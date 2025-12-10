/*
	This file contains code for applying a separable 2D Gaussian filter using CUDA. Both passes of the convolution
	are performed on the GPU.
	Limitations: This algorithm doesn't currently use shared memory, so it could theoretically be made faster
	     for large images.
*/

#pragma once
#include <tira/cuda/error.h>
#include <tira/cuda/callable.h>
#include <tira/functions/filter.h>

namespace tira {
	namespace cuda {

		CUDA_CALLABLE float normal(float x, float mean, float sigma) {
			float c = 1.0f / sqrtf(2.0f * M_PI * sigma * sigma);
			float e = expf(-((x - mean) * (x - mean)) / (2.0f * sigma * sigma));
			return c * e;
		}

		inline std::vector<float> kernel_gaussian_cuda(unsigned int ksize, float sigma) {
			std::vector<float> kernel(ksize);

			const float dx = 1.0f;
			const float halfwidth = static_cast<float>(ksize) * dx / 2.0f;
			const float startx = -halfwidth + (dx / static_cast<float>(2));

			float sum = 0.0f;
			for (unsigned int xi = 0; xi < ksize; ++xi) {
				const float x = startx + static_cast<float>(xi) * dx;
				kernel[xi] = normal(x, 0.0f, sigma);
				sum += kernel[xi];
			}
			if (sum > 0.0f) {
				float inv = 1.0f / sum;
				for (float& v : kernel) v *= inv;
			}
			return kernel;
		}

		/// 1D convolution along the Y (slow) axis
		template<typename T>
		__global__ void kernel_Convolve1DY(const T* source, T* dest,
			const unsigned int width, const unsigned int out_height,
			const float* kernel, const unsigned int kernel_len) {

			unsigned int xi = blockIdx.x * blockDim.x + threadIdx.x;		// calculate the coordinates of the output
			unsigned int yi = blockIdx.y * blockDim.y + threadIdx.y;
			if (yi >= out_height || xi >= width) return;					// return if we are out of bounds of the output image

			T sum = (T)0;														// initialize the running sum
			for (unsigned int ki = 0; ki < kernel_len; ki++) {				// for each sample in the convolution kernel
				sum = sum + kernel[ki] * source[(yi + ki) * width + xi];	// calculate the sample contribution
			}
			dest[yi * width + xi] = sum;									// assign the result to the output image
		}

		/// 1D convolution along the X (fast) axis
		template<typename T>
		__global__ void kernel_Convolve1DX(T* source, T* dest,
			unsigned int width, unsigned int out_width, unsigned int out_height,
			float* kernel, unsigned int kernel_len) {

			unsigned int xi = blockIdx.x * blockDim.x + threadIdx.x;		// calculate the coordinates into the output image
			unsigned int yi = blockIdx.y * blockDim.y + threadIdx.y;
			if (yi >= out_height || xi >= out_width) return;				// return if out of bounds of the output image

			T sum = (T)0;														// initialize the summation
			for (unsigned int ki = 0; ki < kernel_len; ki++) {				// for each sample in the convolution kernel
				sum = sum + kernel[ki] * source[yi * width + xi + ki];		// calculate the sample contribution
			}
			dest[yi * out_width + xi] = sum;								// assign the result to the output image
		}

		template<typename T>
		__global__ void kernel_convolve3_x(const T* source, T* dest, unsigned int width, unsigned int out_width, unsigned int out_height,
			unsigned int out_depth, float* kernel, unsigned int K) {
			unsigned int xi = blockDim.x * blockIdx.x + threadIdx.x;
			unsigned int yi = blockDim.y * blockIdx.y + threadIdx.y;
			unsigned int zi = blockDim.z * blockIdx.z + threadIdx.z;

			if (xi >= out_width || yi >= out_height || zi >= out_depth) return;

			T conv = (T)0;
			for (unsigned int ki = 0; ki < K; ki++)
				conv = conv + source[(zi * out_height + yi) * width + (xi + ki)] * kernel[ki];

			dest[(zi * out_height + yi) * out_width + xi] = conv;
		}

		template<typename T>
		__global__ void kernel_convolve3_y(const T* source, T* dest, unsigned int width, unsigned int height, unsigned int out_height, 
			unsigned int out_depth, float* kernel, unsigned int K) {
			unsigned int xi = blockDim.x * blockIdx.x + threadIdx.x;
			unsigned int yi = blockDim.y * blockIdx.y + threadIdx.y;
			unsigned int zi = blockDim.z * blockIdx.z + threadIdx.z;

			if (xi >= width || yi >= out_height || zi >= out_depth) return;

			T conv = (T)0;
			for (unsigned int ki = 0; ki < K; ki++)
				conv = conv + source[(zi * height + (yi + ki)) * width + xi] * kernel[ki];

			dest[(zi * out_height + yi) * width + xi] = conv;
		}

		template<typename T>
		__global__ void kernel_convolve3_z(const T* source, T* dest, unsigned int width, unsigned int height,
			unsigned int out_depth, float* kernel, unsigned int K) {
			unsigned int xi = blockDim.x * blockIdx.x + threadIdx.x;
			unsigned int yi = blockDim.y * blockIdx.y + threadIdx.y;
			unsigned int zi = blockDim.z * blockIdx.z + threadIdx.z;

			if (xi >= width || yi >= height || zi >= out_depth) return;

			T conv = (T)0;
			for (unsigned int ki = 0; ki < K; ki++)
				conv = conv + source[((zi + ki) * height + yi) * width + xi] * kernel[ki];

			dest[(zi * height + yi) * width + xi] = conv;
		}



		/// <summary>
		/// Convolve a 2D image by a separable Gaussian kernel
		/// </summary>
		/// <typeparam name="T">data type for the input</typeparam>
		/// <param name="source">GPU pointer to the source image array</param>
		/// <param name="width">width of the source image array</param>
		/// <param name="height">height of the source image array</param>
		/// <param name="sigma1">standard deviation of the kernel along the x dimension</param>
		/// <param name="sigma2">standard deviation of the kernel along the y dimension</param>
		/// <param name="out_width">width of the output image after the convolution</param>
		/// <param name="out_height">height of the output image after the convolution</param>
		/// <returns></returns>
		template<typename T>
		T* gaussian_filter2d(const T* source, unsigned int width, unsigned int height,
			float sigma1, float sigma2,
			unsigned int& out_width, unsigned int& out_height) {

			unsigned int window_size_w = (unsigned int)(sigma1 * 6 + 1);		// calculate the window sizes for each kernel
			unsigned int window_size_h = (unsigned int)(sigma2 * 6 + 1);
			out_width = width - (window_size_w - 1);									// calculate the size of the output image
			out_height = height - (window_size_h - 1);

			size_t bytes = sizeof(T) * width * height;							// calculate the number of bytes in the image

			const T* gpu_source;														// create a pointer for the GPU source image
			T* temp_source;

			// determine if the source image is provided on the CPU or GPU
			cudaPointerAttributes attribs;										// create a pointer attribute structure
			HANDLE_ERROR(cudaPointerGetAttributes(&attribs, source));			// get the attributes for the source pointer

			if (attribs.type == cudaMemoryTypeDevice) {							// if the provided pointer is on the device
				gpu_source = source;											// set the gpu_source pointer to source
			}
			else {																// otherwise copy the source image to the GPU
				HANDLE_ERROR(cudaMalloc(&temp_source, bytes));								// allocate space on the GPU for the source image
				HANDLE_ERROR(cudaMemcpy(temp_source, source, bytes, cudaMemcpyHostToDevice));// copy the source image to the GPU
				gpu_source = static_cast<const T*>(temp_source);
			}

			/////////////// Calculate Convolution Kernels
			float* kernel_h = (float*)malloc(sizeof(float) * window_size_h);	// allocate space for the y kernel
			int xh = -(int)(window_size_h / 2);									// calculate the starting coordinate for the kernel
			for (int j = 0; j < window_size_h; j++) {							// for each pixel in the kernel
				kernel_h[j] = normal(xh + j, 0, sigma2);						// calculate the Gaussian value
			}

			// allocate space on the GPU for the y-axis kernel and copy it
			float* gpu_kernel_h;
			HANDLE_ERROR(cudaMalloc(&gpu_kernel_h, sizeof(float) * window_size_h));
			HANDLE_ERROR(cudaMemcpy(gpu_kernel_h, kernel_h, sizeof(float) * window_size_h, cudaMemcpyHostToDevice));


			float* kernel_w = (float*)malloc(sizeof(float) * window_size_w);	// allocate space for the x kernel
			int xw = -(int)(window_size_w / 2);									// calculate the starting coordinate for the kernel
			for (int i = 0; i < window_size_w; i++) {							// for each pixel in the kernel
				kernel_w[i] = normal(xw + i, 0, sigma1);						// calculate the Gaussian value
			}

			// allocate space on the GPU for the x-axis kernel and copy it
			float* gpu_kernel_w;
			HANDLE_ERROR(cudaMalloc(&gpu_kernel_w, sizeof(float) * window_size_w));
			HANDLE_ERROR(cudaMemcpy(gpu_kernel_w, kernel_w, sizeof(float) * window_size_w, cudaMemcpyHostToDevice));
			///////////////End Calculate Convolution Kernels

			// This is a separable convolution and will be done in two passes
			size_t bytes_firstpass = sizeof(T) * width * out_height;			// calculate the number of bytes in the first pass output
			T* gpu_firstpass;													// allocate space on the GPU for the first pass
			HANDLE_ERROR(cudaMalloc(&gpu_firstpass, bytes_firstpass));

			// get the active device properties to calculate the optimal the block size
			int device;
			HANDLE_ERROR(cudaGetDevice(&device));
			cudaDeviceProp props;
			HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
			unsigned int max_threads = props.maxThreadsPerBlock;
			dim3 blockDim = { (unsigned int)std::sqrt(max_threads), (unsigned int)std::sqrt(max_threads) };


			dim3 gridDim = { width / blockDim.x + 1, out_height / blockDim.y + 1 };	// calculate the grid size for the first pass
			// Run the kernel to calculate the first pass of the separable filter. The kernel takes the original source
			// image and outputs an image that is smaller along the y axis.
			kernel_Convolve1DY << <gridDim, blockDim >> > (gpu_source, gpu_firstpass, width, out_height, gpu_kernel_h, window_size_h);

			T* gpu_out;												// calculate the size of the final output image
			size_t out_bytes = sizeof(T) * out_width * out_height;
			HANDLE_ERROR(cudaMalloc(&gpu_out, out_bytes));			// allocate space on the GPU for the final output
			gridDim.x = out_width / blockDim.x + 1;					// re-calculate the grid dimension to match the final pass output
			// Run the kernel to calculate the second pass of the separable Gaussian filter. The kernel takes the first pass
			// as input and produces the final output image on the GPU
			kernel_Convolve1DX << <gridDim, blockDim >> > (gpu_firstpass, gpu_out, width, out_width, out_height, gpu_kernel_w, window_size_w);

			// allocate space for the output image on the CPU and copy the final result
			T* out;
			if (attribs.type == cudaMemoryTypeDevice) {				// if the source was on the device, return an output on the device
				out = gpu_out;
			}
			else {													// otherwise copy the output to the host
				out = (T*)malloc(out_bytes);
				HANDLE_ERROR(cudaMemcpy(out, gpu_out, out_bytes, cudaMemcpyDeviceToHost));
			}

			// free everything
			if (attribs.type == cudaMemoryTypeHost) {				// if the source pointer was on the host, free the interim GPU versions
				HANDLE_ERROR(cudaFree(temp_source));
				HANDLE_ERROR(cudaFree(gpu_out));
			}

			HANDLE_ERROR(cudaFree(gpu_firstpass));
			HANDLE_ERROR(cudaFree(gpu_kernel_w));
			HANDLE_ERROR(cudaFree(gpu_kernel_h));

			return out;

		}

		/// <summary>
		/// Convolve a 3D image by a separable Gaussian kernel
		/// </summary>
		/// <typeparam name="T">data type for the input</typeparam>
		/// <param name="source">GPU pointer to the source image array</param>
		/// <param name="width">width of the source image array</param>
		/// <param name="height">height of the source image array</param>
		/// <param name="sigma_w">standard deviation of the kernel along the x dimension</param>
		/// <param name="sigma_h">standard deviation of the kernel along the y dimension</param>
		/// <param name="sigma_d">standard deviation of the kernel along the z dimension</param>
		/// <param name="out_width">width of the output image after the convolution</param>
		/// <param name="out_height">height of the output image after the convolution</param>
		/// <param name="out_depth">depth of the output image after the convolution</param>
		/// <returns></returns>
		template<typename T>
		T* gaussian_filter3d(const T* source, unsigned int width, unsigned int height, unsigned int depth,
			float sigma_w, float sigma_h, float sigma_d,
			unsigned int window_size_w, unsigned int window_size_h, unsigned int window_size_d,
			unsigned int& out_width, unsigned int& out_height, unsigned int& out_depth) {

			out_width  = width  - (window_size_w - 1);							// calculate the size of the output image (valid convolution)
			out_height = height - (window_size_h - 1);
			out_depth  = depth  - (window_size_d - 1);

			// Calculate the number of bytes in the image
			size_t bytes = sizeof(T) * static_cast<size_t>(width) * static_cast<size_t>(height) * static_cast<size_t>(depth);							

			const T* gpu_source = nullptr;										// create a pointer for the GPU source image
			T* temp_source = nullptr;

			// Determine if the source pointer is on CPU or GPU
			cudaPointerAttributes attribs;										// create a pointer attribute structure
			HANDLE_ERROR(cudaPointerGetAttributes(&attribs, source));			// get the attributes for the source pointer

			if (attribs.type == cudaMemoryTypeDevice) {							// if the provided pointer is on the device
				gpu_source = source;											// set the gpu_source pointer to source
			}
			else {																// otherwise copy the source image to the GPU
				HANDLE_ERROR(cudaMalloc(&temp_source, bytes));								// allocate space on the GPU for the source image
				HANDLE_ERROR(cudaMemcpy(temp_source, source, bytes, cudaMemcpyHostToDevice));// copy the source image to the GPU
				gpu_source = static_cast<const T*>(temp_source);
			}

			// -----  Build 1D Gaussian kernels (host), normalize, then upload to GPU -----
			
			// Z kernel (depth)
			auto kernel_d = kernel_gaussian_cuda(window_size_d, sigma_d);
			float* gpu_kernel_d = nullptr;
			HANDLE_ERROR(cudaMalloc(&gpu_kernel_d, sizeof(float) * window_size_d));
			HANDLE_ERROR(cudaMemcpy(gpu_kernel_d, kernel_d.data(), sizeof(float) * window_size_d, cudaMemcpyHostToDevice));

			// Y kernel (height)
			auto kernel_h = kernel_gaussian_cuda(window_size_h, sigma_h);
			float* gpu_kernel_h = nullptr;
			HANDLE_ERROR(cudaMalloc(&gpu_kernel_h, sizeof(float) * window_size_h));
			HANDLE_ERROR(cudaMemcpy(gpu_kernel_h, kernel_h.data(), sizeof(float) * window_size_h, cudaMemcpyHostToDevice));

			// X kernel (width)
			auto kernel_w = kernel_gaussian_cuda(window_size_w, sigma_w);
			float* gpu_kernel_w = nullptr;
			HANDLE_ERROR(cudaMalloc(&gpu_kernel_w, sizeof(float) * window_size_w));
			HANDLE_ERROR(cudaMemcpy(gpu_kernel_w, kernel_w.data(), sizeof(float) * window_size_w, cudaMemcpyHostToDevice));

			// -------------------- End Build Kernels ----------------------
			 
			// ----- Launch configuration (separable convolution in 3 passes) -----
			int device = 0;
			HANDLE_ERROR(cudaGetDevice(&device));
			cudaDeviceProp props{};
			HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
			unsigned int max_threads = props.maxThreadsPerBlock;
			unsigned int block_xy = static_cast<unsigned int>(std::sqrt(static_cast<float>(max_threads)));
			dim3 blockDim(block_xy, block_xy, 1);


			// --- First pass: Z convolution (reduces depth to out_depth) ---
			T* gpu_firstpass = nullptr;															// allocate space on the GPU for the first pass
			const size_t bytes_firstpass = sizeof(T) * static_cast<size_t>(width) * static_cast<size_t>(height) * static_cast<size_t>(out_depth);
			HANDLE_ERROR(cudaMalloc(&gpu_firstpass, bytes_firstpass));

			dim3 gridDimZ(
				(width + blockDim.x - 1) / blockDim.x,
				(height + blockDim.y - 1) / blockDim.y,
				(out_depth + blockDim.z - 1) / blockDim.z
			);
			kernel_convolve3_z << <gridDimZ, blockDim >> > (gpu_source, gpu_firstpass, width, height, out_depth, gpu_kernel_d, window_size_d);
			HANDLE_ERROR(cudaGetLastError());
			HANDLE_ERROR(cudaDeviceSynchronize());

			if (attribs.type != cudaMemoryTypeDevice) {
				HANDLE_ERROR(cudaFree(temp_source));
				temp_source = nullptr;
			}

			// --- Second pass: Y convolution (reduces height to out_height) ---
			T* gpu_secondpass = nullptr;															// allocate space on the GPU for the second pass
			const size_t bytes_secondpass = sizeof(T) * static_cast<size_t>(width) * static_cast<size_t>(out_height) * static_cast<size_t>(out_depth);
			HANDLE_ERROR(cudaMalloc(&gpu_secondpass, bytes_secondpass));

			// Second pass - Y convolution
			dim3 gridDimY(
				(width + blockDim.x - 1) / blockDim.x,
				(out_height + blockDim.y - 1) / blockDim.y,
				(out_depth + blockDim.z - 1) / blockDim.z
			);
			kernel_convolve3_y << <gridDimY, blockDim >> > (gpu_firstpass, gpu_secondpass, width, height, out_height, out_depth, gpu_kernel_h, window_size_h);
			HANDLE_ERROR(cudaGetLastError());
			HANDLE_ERROR(cudaDeviceSynchronize());
			HANDLE_ERROR(cudaFree(gpu_firstpass));

			// --- Third pass: X convolution (reduces width to out_width) ---
			T* gpu_out = nullptr;																	// calculate the size of the final output image
			const size_t out_bytes = sizeof(T) * static_cast<size_t>(out_width) * static_cast<size_t>(out_height) * static_cast<size_t>(out_depth);
			HANDLE_ERROR(cudaMalloc(&gpu_out, out_bytes));								// allocate space on the GPU for the final output

			// Third/Last pass - X convolution
			dim3 gridDimX(
				(out_width + blockDim.x - 1) / blockDim.x,
				(out_height + blockDim.y - 1) / blockDim.y,
				(out_depth + blockDim.z - 1) / blockDim.z
			);
			kernel_convolve3_x << <gridDimX, blockDim >> > (gpu_secondpass, gpu_out, width, out_width, out_height, out_depth, gpu_kernel_w, window_size_w);
			HANDLE_ERROR(cudaGetLastError());
			HANDLE_ERROR(cudaDeviceSynchronize());
			HANDLE_ERROR(cudaFree(gpu_secondpass));

			// ----- Copy back or return device pointer
			T* out = nullptr;
			if (attribs.type == cudaMemoryTypeDevice) {				// caller gave a device pointer, return a device pointer
				out = gpu_out;
			}
			else {													// otherwise copy the output to the host
				out = static_cast<T*>(malloc(out_bytes));
				HANDLE_ERROR(cudaMemcpy(out, gpu_out, out_bytes, cudaMemcpyDeviceToHost));
				HANDLE_ERROR(cudaFree(gpu_out));
			}

			// ----- Clean up kernels -----
			HANDLE_ERROR(cudaFree(gpu_kernel_w));
			HANDLE_ERROR(cudaFree(gpu_kernel_h));
			HANDLE_ERROR(cudaFree(gpu_kernel_d));

			return out;
		}
	}
}