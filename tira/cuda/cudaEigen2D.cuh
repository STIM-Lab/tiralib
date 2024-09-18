/*
	This file implements an GPU-based eigendecompositions for arrays of 2x2 matrices.
*/

#pragma once
#include <tira/cuda/error.h>
#include <tira/cuda/callable.h>

namespace tira {

    namespace cuda {

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

        // small then large
        template<typename T>
        CUDA_CALLABLE void evec2D(T* matrix, T* lambdas, T& theta0, T& theta1) {
            //[a  b]
            //[c  d]

            float a = matrix[0];
            float c = matrix[1];
            float b = matrix[2];
            float d = matrix[3];

            float v[2];
            float len;

            if (c != 0) {
                theta0 = std::atan2(c, lambdas[0] - d);
                theta1 = std::atan2(c, lambdas[1] - d);
            }
            else if (b != 0) {
                theta0 = std::atan2(lambdas[0] - a, b);
                theta1 = std::atan2(lambdas[1] - a, b);
            }
            else {
                theta0 = std::atan2(0, 1);
                theta1 = std::atan2(1, 0);
            }
        }

        template<typename T>
        __global__ void kernel_evec2D(T* mats, T* evals, size_t n, T* evecs) {
            unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
            if (i >= n) return;
            evec2D(&mats[i * 4], &evals[i * 2], evecs[i * 2 + 0], evecs[i * 2 + 1]);
        }

        template<typename T>
        T* Eigenvectors2D(T* mats, T* evals, size_t n) {

            T* gpu_mats;
            T* gpu_evals;
            size_t mats_bytes = sizeof(T) * 4 * n;
            size_t ev_bytes = sizeof(T) * 2 * n;

            // determine if the source image is provided on the CPU or GPU
            cudaPointerAttributes attribs;										// create a pointer attribute structure
            HANDLE_ERROR(cudaPointerGetAttributes(&attribs, mats));			// get the attributes for the source pointer

            if (attribs.type == cudaMemoryTypeDevice) {							// if the provided pointer is on the device
                gpu_mats = mats;									            // set the gpu_source pointer to source
                gpu_evals = evals;
            }
            else {																// otherwise copy the source image to the GPU
                HANDLE_ERROR(cudaMalloc(&gpu_mats, mats_bytes));								// allocate space on the GPU for the source image
                HANDLE_ERROR(cudaMemcpy(gpu_mats, mats, mats_bytes, cudaMemcpyHostToDevice));// copy the source image to the GPU
                HANDLE_ERROR(cudaMalloc(&gpu_evals, ev_bytes));
                HANDLE_ERROR(cudaMemcpy(gpu_evals, evals, ev_bytes, cudaMemcpyHostToDevice));
            }

            // get the active device properties to calculate the optimal the block size
            int device;
            HANDLE_ERROR(cudaGetDevice(&device));
            cudaDeviceProp props;
            HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
            unsigned int max_threads = props.maxThreadsPerBlock;
            dim3 blockDim = max_threads;
            dim3 gridDim = n / blockDim.x + 1;	                                // calculate the grid size for the first pass

            T* gpu_evecs;
            HANDLE_ERROR(cudaMalloc(&gpu_evecs, ev_bytes));
            kernel_evec2D << <gridDim, blockDim >> > (gpu_mats, gpu_evals, n, gpu_evecs);

            T* evecs;
            if (attribs.type == cudaMemoryTypeDevice) {
                evecs = gpu_evecs;
            }
            else {
                evecs = (T*)malloc(ev_bytes);
                HANDLE_ERROR(cudaMemcpy(evecs, gpu_evecs, ev_bytes, cudaMemcpyDeviceToHost));
                HANDLE_ERROR(cudaFree(gpu_evecs));
                HANDLE_ERROR(cudaFree(gpu_evals));
                HANDLE_ERROR(cudaFree(gpu_mats));
            }

            return evecs;

        }
    }
}