#pragma once

#include <numbers>
#include <vector>
#include<memory>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <tira/cuda/callable.h>

namespace tira {

	template <typename T>
	CUDA_CALLABLE T normaldist(T mu, T sigma, T x) {
		T sigma_sq = sigma * sigma;
		T pi = static_cast<T>(M_PI);
		T n = static_cast<T>(1) / std::sqrt(static_cast<T>(2) * pi * sigma_sq);

		T exponent = -((x - mu) * (x - mu)) / (static_cast<T>(2) * sigma_sq);
		return n * exp(exponent);
	}
}

namespace tira::cpu {

	/**
	 * @brief Generate a 1D Gaussian kernel sampled on a regular grid.
	 *
	 * @tparam Type      Floating-point type for kernel values
	 * @param ksize		 Number of samples in the kernel
	 * @param mu         Mean of the Gaussian (typically 0)
	 * @param sigma      Standard deviation of the Gaussian
	 * @param dx         Sample spacing of the grid
	 * @return a pointer to an array containing ksize samples of a Gaussian kernel
	 */
	template <typename Type>
	static Type* gaussian1d(unsigned int ksize, Type mu, Type sigma, Type dx) {
		Type* kernel = new Type[ksize];

		// Physical length of the kernel
		const Type halfwidth = static_cast<Type>(ksize) * dx / static_cast<Type>(2);

		// First sample is centered at -halfwidth + dx/2 so the grid is symmetric around mu
		const Type startx = -halfwidth + (dx / static_cast<Type>(2));

		// Sample the Gaussian at regularly spaced positions
		Type sum = static_cast<Type>(0);
		for (unsigned int xi = 0; xi < ksize; ++xi) {
			const Type x = startx + static_cast<Type>(xi) * dx;
			kernel[xi] = tira::normaldist(mu, sigma, x);
			sum += kernel[xi];
		}

		// Normalize the kernel to have sum 1
		if (sum > static_cast<Type>(0)) {
			Type inv_sum = static_cast<Type>(1) / sum;
			for (unsigned i = 0; i < ksize; i++)
				kernel[i] *= inv_sum;
			//for (auto& v : kernel) v *= inv_sum;
		}

		return kernel;
	}

	template <typename ImageType, typename KernelType>
	static ImageType* convolve2(const ImageType* input, const unsigned in_sx, const unsigned in_sy,
		const KernelType* kernel, const unsigned k_sx, const unsigned k_sy,
		unsigned& out_sx, unsigned& out_sy) {

		out_sx = in_sx - (k_sx - 1);
		out_sy = in_sy - (k_sy - 1);
		ImageType* output = (ImageType*)malloc(out_sx * out_sy * sizeof(ImageType));

		for (unsigned yi = 0; yi < out_sy; yi++) {
			for (unsigned xi = 0; xi < out_sx; xi++) {
				ImageType k_sum = (ImageType)0;
				for (unsigned vi = 0; vi < k_sy; vi++) {
					for (unsigned ui = 0; ui < k_sx; ui++) {
						k_sum += kernel[vi * k_sx + ui] * input[(yi + vi) * in_sx + (xi + ui)];
					}
				}
				output[yi * out_sx + xi] = k_sum;
			}
		}
		return output;
	}

	/**
	 * @brief Perform a 3D valid convolution of a volume with a dense kernel.
	 *
	 * @tparam ImageType   Element type of the input/output volume.
	 * @tparam KernelType  Element type of the convolution kernel.
	 * @param input        Pointer to input volume (size in_s0 * in_s1 * in_s2).
	 * @param in_s0        Input size along x.
	 * @param in_s1        Input size along y.
	 * @param in_s2        Input size along z.
	 * @param kernel       Pointer to kernel (size k_s0 * k_s1 * k_s2).
	 * @param k_s0         Kernel size along x.
	 * @param k_s1         Kernel size along y.
	 * @param k_s2         Kernel size along z.
	 * @param[out] out_s0  Output size along x.
	 * @param[out] out_s1  Output size along y.
	 * @param[out] out_s2  Output size along z.
	 * @return ImageType*  Newly allocated output volume; caller must delete[] it.
	 */
	template <typename ImageType, typename KernelType>
	static ImageType* convolve3(const ImageType* input, const unsigned in_s0, const unsigned in_s1,const unsigned in_s2,
		const KernelType* kernel, const unsigned k_s0, const unsigned k_s1, const unsigned k_s2,
		unsigned& out_s0, unsigned& out_s1, unsigned& out_s2) {

		out_s0 = in_s0 - (k_s0 - 1);
		out_s1 = in_s1 - (k_s1 - 1);
		out_s2 = in_s2 - (k_s2 - 1);

		const std::size_t out_size = static_cast<std::size_t>(out_s0) * static_cast<std::size_t>(out_s1) * static_cast<std::size_t>(out_s2);

		auto output = std::make_unique<ImageType[]>(out_size);
		// ImageType* output = (ImageType*)malloc(out_s0 * out_s1 * out_s2 * sizeof(ImageType));

		for (unsigned zi = 0; zi < out_s2; ++zi) {
			for (unsigned yi = 0; yi < out_s1; ++yi) {
				for (unsigned xi = 0; xi < out_s0; ++xi) {

					ImageType k_sum = static_cast<ImageType>(0);

					// Full 3D kernel accumulation
					for (unsigned wi = 0; wi < k_s2; ++wi) {
						for (unsigned vi = 0; vi < k_s1; ++vi) {
							for (unsigned ui = 0; ui < k_s0; ++ui) {
								const unsigned k_idx = wi * k_s1 * k_s0 + vi * k_s0 + ui;
								const unsigned in_idx = (zi + wi) * in_s1 * in_s0 + (yi + vi) * in_s0 + (xi + ui);
								k_sum += kernel[k_idx] * input[in_idx];
							}
						}
					}
					const unsigned out_idx = zi * out_s1 * out_s0 + yi * out_s0 + xi;
					output[out_idx] = k_sum;
				}
			}
		}
		// Hand ownership to caller, must delete[] the returned pointer
		return output.release();
	}

	/**
	 *
	 * @tparam ImageType data type used to represent the image
	 * @tparam KernelType data type used to represent the kernel (default float)
	 * @param input pointer to the input image to be convolved
	 * @param in_sx size of the input image along the fast axis
	 * @param in_sy size of the input image along the slow axis
	 * @param sigma standard deviation of the Gaussian kernel in pixels
	 * @param out_sx size of the output image (consists of the valid portion of the convolution: out_sx <= in_sx)
	 * @param out_sy size of the output image (consists of the valid portion of the convolution: out_sy <= in_sy)
	 * @return a new image consisting of the valid region of the Gaussian convolution, where the size of the image is given by out_sx and out_sy
	 */
	template <typename ImageType, typename KernelType = float>
	static ImageType* gaussian_convolve2(const ImageType* input, const unsigned in_sx, const unsigned in_sy,
		KernelType sigma,
		unsigned& out_sx, unsigned& out_sy) {

		// Calculate a 1D kernel
		unsigned ksize = (unsigned)(6 * sigma);
		KernelType* kernel = gaussian1d(ksize, (KernelType)0, sigma, (KernelType)1);
		for (unsigned i = 0; i < ksize; i++)
			printf("%f", kernel[i]);

		ImageType* dest_x = tira::cpu::convolve2<ImageType, KernelType>(input, in_sx, in_sy, kernel, ksize, 1, out_sx, out_sy);
		ImageType* dest_xy = tira::cpu::convolve2<ImageType, KernelType>(dest_x, out_sx, out_sy, kernel, 1, ksize, out_sx, out_sy);

		delete kernel;
		delete dest_x;
		return dest_xy;
	}
};

#ifdef __CUDACC__
#include <tira/cuda/error.h>

namespace tira::cuda {

	template<typename T>
	__global__ void cuda_convolve1y(const T* source, T* dest,
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
	__global__ void cuda_convolve1x(T* source, T* dest,
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
	template <typename ImageType, typename KernelType = float>
	ImageType* gaussian_convolve2(const ImageType* source, const unsigned in_sx, const unsigned in_sy,
		float sigma1, float sigma2, unsigned int& out_width, unsigned int& out_height) {

		unsigned int window_size_w = (unsigned int)(sigma1 * 6 + 1);		// calculate the window sizes for each kernel
		unsigned int window_size_h = (unsigned int)(sigma2 * 6 + 1);
		out_width = in_sx - window_size_w;									// calculate the size of the output image
		out_height = in_sy - window_size_h;

		size_t bytes = sizeof(ImageType) * in_sx * in_sy;							// calculate the number of bytes in the image

		const ImageType* gpu_source;														// create a pointer for the GPU source image
		ImageType* temp_source;												// temporary image used during the separable convolution

		// determine if the source image is provided on the CPU or GPU
		cudaPointerAttributes attribs;										// create a pointer attribute structure
		HANDLE_ERROR(cudaPointerGetAttributes(&attribs, source));			// get the attributes for the source pointer

		if (attribs.type == cudaMemoryTypeDevice) {							// if the provided pointer is on the device
			gpu_source = source;											// set the gpu_source pointer to source
		}
		else {																// otherwise copy the source image to the GPU
			HANDLE_ERROR(cudaMalloc(&temp_source, bytes));								// allocate space on the GPU for the source image
			HANDLE_ERROR(cudaMemcpy(temp_source, source, bytes, cudaMemcpyHostToDevice));// copy the source image to the GPU
			gpu_source = static_cast<const ImageType*>(temp_source);
		}

		/////////////// Calculate Convolution Kernels

		KernelType* kernel_h = cpu::gaussian1d<KernelType>(window_size_h, 0, sigma2, 1);


		// allocate space on the GPU for the y-axis kernel and copy it
		KernelType* gpu_kernel_h;
		HANDLE_ERROR(cudaMalloc(&gpu_kernel_h, sizeof(KernelType) * window_size_h));
		HANDLE_ERROR(cudaMemcpy(gpu_kernel_h, kernel_h, sizeof(KernelType) * window_size_h, cudaMemcpyHostToDevice));


		KernelType* kernel_w = cpu::gaussian1d<KernelType>(window_size_w, 0, sigma1, 1);

		// allocate space on the GPU for the x-axis kernel and copy it
		KernelType* gpu_kernel_w;
		HANDLE_ERROR(cudaMalloc(&gpu_kernel_w, sizeof(float) * window_size_w));
		HANDLE_ERROR(cudaMemcpy(gpu_kernel_w, kernel_w, sizeof(float) * window_size_w, cudaMemcpyHostToDevice));
		///////////////End Calculate Convolution Kernels

		// This is a separable convolution and will be done in two passes
		size_t bytes_firstpass = sizeof(ImageType) * in_sx * out_height;			// calculate the number of bytes in the first pass output
		ImageType* gpu_firstpass;													// allocate space on the GPU for the first pass
		HANDLE_ERROR(cudaMalloc(&gpu_firstpass, bytes_firstpass));

		// get the active device properties to calculate the optimal the block size
		int device;
		HANDLE_ERROR(cudaGetDevice(&device));
		cudaDeviceProp props;
		HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
		unsigned int max_threads = props.maxThreadsPerBlock;
		dim3 blockDim = { (unsigned int)std::sqrt(max_threads), (unsigned int)std::sqrt(max_threads) };


		dim3 gridDim = { in_sx / blockDim.x + 1, out_height / blockDim.y + 1 };	// calculate the grid size for the first pass
		// Run the kernel to calculate the first pass of the separable filter. The kernel takes the original source
		// image and outputs an image that is smaller along the y axis.
		cuda_convolve1y << <gridDim, blockDim >> > (gpu_source, gpu_firstpass, in_sx, out_height, gpu_kernel_h, window_size_h);

		ImageType* gpu_out;												// calculate the size of the final output image
		size_t out_bytes = sizeof(ImageType) * out_width * out_height;
		HANDLE_ERROR(cudaMalloc(&gpu_out, out_bytes));			// allocate space on the GPU for the final output
		gridDim.x = out_width / blockDim.x + 1;					// re-calculate the grid dimension to match the final pass output
		// Run the kernel to calculate the second pass of the separable Gaussian filter. The kernel takes the first pass
		// as input and produces the final output image on the GPU
		cuda_convolve1x << <gridDim, blockDim >> > (gpu_firstpass, gpu_out, in_sx, out_width, out_height, gpu_kernel_w, window_size_w);

		// allocate space for the output image on the CPU and copy the final result
		ImageType* out;
		if (attribs.type == cudaMemoryTypeDevice) {				// if the source was on the device, return an output on the device
			out = gpu_out;
		}
		else {													// otherwise copy the output to the host
			out = (ImageType*)malloc(out_bytes);
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
}
#endif
