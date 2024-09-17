/*
	This file implements an GPU-based eigendecompositions for arrays of 2x2 matrices.
*/

#pragma once
#include <tira/cuda/error.h>
#include <tira/cuda/callable.h>

template<typename T>
CUDA_CALLABLE void eval2D(T* matrix, T& eval0, T& eval1) {
    T d = matrix[0];
    T e = matrix[1];
    T f = matrix[2];
    T g = matrix[3];

    float dpg = d + g;
    float disc = sqrt((4 * e * f) + pow(d - g, 2));
    float a = (dpg + disc) / 2.0f;
    float b = (dpg - disc) / 2.0f;
    eval0 = a < b ? a : b;
    eval1 = a > b ? a : b;
}

template<typename T>
__global__ void kernel_eval2D(T* mats, size_t n, T* evals) {
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    eval2D(&mats[i * 4], evals[i * 2 + 0], evals[i * 2 + 1]);
}

template<typename T>
T* Eigenvalues2D(T* mats, size_t n) {

    T* gpu_mats;
    size_t mats_bytes = sizeof(T) * 4 * n;
    size_t evals_bytes = sizeof(T) * 2 * n;

    // determine if the source image is provided on the CPU or GPU
    cudaPointerAttributes attribs;										// create a pointer attribute structure
    HANDLE_ERROR(cudaPointerGetAttributes(&attribs, mats));			// get the attributes for the source pointer

    if (attribs.type == cudaMemoryTypeDevice) {							// if the provided pointer is on the device
        gpu_mats = mats;									            // set the gpu_source pointer to source
    }
    else {																// otherwise copy the source image to the GPU
        HANDLE_ERROR(cudaMalloc(&gpu_mats, mats_bytes));								// allocate space on the GPU for the source image
        HANDLE_ERROR(cudaMemcpy(gpu_mats, mats, mats_bytes, cudaMemcpyHostToDevice));// copy the source image to the GPU
    }

    // get the active device properties to calculate the optimal the block size
    int device;
    HANDLE_ERROR(cudaGetDevice(&device));
    cudaDeviceProp props;
    HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
    unsigned int max_threads = props.maxThreadsPerBlock;
    dim3 blockDim = max_threads;
    dim3 gridDim = n / blockDim.x + 1;	                                // calculate the grid size for the first pass
    
    T* gpu_evals;
    HANDLE_ERROR(cudaMalloc(&gpu_evals, evals_bytes));
    kernel_eval2D << <gridDim, blockDim >> > (gpu_mats, n, gpu_evals);

    T* evals;
    if (attribs.type == cudaMemoryTypeDevice) {
        evals = gpu_evals;
    }
    else {
        evals = (T*)malloc(sizeof(T) * 2 * n);
        HANDLE_ERROR(cudaMemcpy(evals, gpu_evals, evals_bytes, cudaMemcpyDeviceToHost));
        HANDLE_ERROR(cudaFree(gpu_evals));
        HANDLE_ERROR(cudaFree(gpu_mats));
    }

    return evals;
}