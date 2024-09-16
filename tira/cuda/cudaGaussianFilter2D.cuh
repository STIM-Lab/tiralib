#pragma once

template<typename T>
T* GaussianFilter2D_cpu(T* source, size_t width, size_t height,
	float sigma1, float sigma2,
	size_t& out_width, size_t& out_height) {

	unsigned int window_size_w = (unsigned int)(sigma1 * 6 + 1);
	unsigned int window_size_h = (unsigned int)(sigma2 * 6 + 1);
	out_width = width - window_size_w;
	out_height = height - window_size_h;

	size_t bytes = sizeof(T) * width * height;		// calculate the number of bytes in the image

	T* gpu_source;
	cudaMalloc(&gpu_source, bytes);					// allocate space on the GPU for the source image
	cudaMemcpy(&gpu_source, source, bytes, cudaMemcpyHostToDevice);

	size_t bytes_firstpass = sizeof(T) * width * out_height;
	T* gpu_firstpass;
	cudaMalloc(&gpu_firstpass, bytes_firstpass);
	cudaMemcpy(gpu_firstpass, &gpu_source, bytes_firstpass, cudaMemcpyDeviceToDevice);

	size_t out_bytes = sizeof(T) * out_width * out_height;
	T* gpu_out;
	cudaMalloc(&gpu_out, out_bytes);
	cudaMemcpy(gpu_out, gpu_firstpass, out_bytes, cudaMemcpyDeviceToDevice);

	T* out = (T*)malloc(out_bytes);
	cudaMemcpy(out, gpu_out, out_bytes, cudaMemcpyDeviceToHost);

	return source;

}