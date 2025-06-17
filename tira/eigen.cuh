/*
	This file implements an GPU-based eigendecompositions for arrays of 2x2 matrices.
*/

#pragma once
#include <tira/eigen.h>
#include <tira/cuda/error.h>

namespace tira::cuda {


    template<typename T>
    __global__ void kernel_eval2D(T* mats, size_t n, T* evals) {
        const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) return;
        eval2D(&mats[i * 4], evals[i * 2 + 0], evals[i * 2 + 1]);
    }

    template<typename T>
    T* Eigenvalues2D(T* mats, size_t n, int device) {

        if (device < 0)                                                     // if the device is < zero, run the CPU version
            return cpu::Eigenvalues2D(mats, n);

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
            evals = new T[2 * n];
            HANDLE_ERROR(cudaMemcpy(evals, gpu_evals, evals_bytes, cudaMemcpyDeviceToHost));
            HANDLE_ERROR(cudaFree(gpu_evals));
            HANDLE_ERROR(cudaFree(gpu_mats));
        }

        return evals;
    }

    

    template<typename T>
    __global__ void kernel_evec2Dpolar(T* mats, T* evals, size_t n, T* evecs) {
        const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) return;
        evec2Dpolar(&mats[i * 4], &evals[i * 2], evecs[i * 2 + 0], evecs[i * 2 + 1]);
    }



    template<typename T>
    T* Eigenvectors2DPolar(T* mats, T* evals, size_t n, int device) {

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
        HANDLE_ERROR(cudaGetDevice(&device));
        cudaDeviceProp props;
        HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
        unsigned int max_threads = props.maxThreadsPerBlock;
        dim3 blockDim = max_threads;
        dim3 gridDim = n / blockDim.x + 1;	                                // calculate the grid size for the first pass

        T* gpu_evecs;
        HANDLE_ERROR(cudaMalloc(&gpu_evecs, ev_bytes));
        kernel_evec2Dpolar << <gridDim, blockDim >> > (gpu_mats, gpu_evals, n, gpu_evecs);

        T* evecs;
        if (attribs.type == cudaMemoryTypeDevice) {
            evecs = gpu_evecs;
        }
        else {
            evecs = new T[2 * n];
            HANDLE_ERROR(cudaMemcpy(evecs, gpu_evecs, ev_bytes, cudaMemcpyDeviceToHost));
            HANDLE_ERROR(cudaFree(gpu_evecs));
            HANDLE_ERROR(cudaFree(gpu_evals));
            HANDLE_ERROR(cudaFree(gpu_mats));
        }

        return evecs;
    }

    

    

    template<typename T>
    __global__ void kernel_eval3D(T* mats, size_t n, T* evals) {
        const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) return;
        eval3D(&mats[i * 9], evals[i * 3 + 0], evals[i * 3 + 1], evals[i * 3 + 2]);
    }

    

    

    template<typename T>
    T* Eigenvalues3D(T* mats, size_t n) {

        T* gpu_mats;
        size_t mats_bytes = sizeof(T) * 9 * n;                              // required bytes for storing the tensor
        size_t evals_bytes = sizeof(T) * 3 * n;                             // required bytes for storing eigenvalues

        // determine if the source volume is provided on the CPU or GPU
        cudaPointerAttributes attribs;										// create a pointer attribute structure
        HANDLE_ERROR(cudaPointerGetAttributes(&attribs, mats));			    // get the attributes for the source pointer

        if (attribs.type == cudaMemoryTypeDevice) {							// if the provided pointer is on the device
            gpu_mats = mats;									            // set the gpu_source pointer to source
        }
        else {																// otherwise copy the source image to the GPU
            HANDLE_ERROR(cudaMalloc(&gpu_mats, mats_bytes));								// allocate space on the GPU for the source image
            HANDLE_ERROR(cudaMemcpy(gpu_mats, mats, mats_bytes, cudaMemcpyHostToDevice));   // copy the source image to the GPU
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
        kernel_eval3D << <gridDim, blockDim >> > (gpu_mats, n, gpu_evals);

        if (attribs.type == cudaMemoryTypeDevice) {
            return gpu_evals;
        }
        else {
            T* evals = new T[3 * n];
            HANDLE_ERROR(cudaMemcpy(evals, gpu_evals, evals_bytes, cudaMemcpyDeviceToHost));
            HANDLE_ERROR(cudaFree(gpu_evals));
            HANDLE_ERROR(cudaFree(gpu_mats));
            return evals;
        }
    }

    template<typename T>
    __global__ void kernel_evec3DPolar(T* mats, T* lambda, size_t n, T* evec) {
        const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) return;

        T theta, phi;
        evec3Dpolar(&mats[i * 9], lambda[i * 3 + 1], theta, phi);
        if (isnan(theta) || isnan(phi)) {
            evec[i * 4 + 0] = PI / 2.0;         // acos(0)
            evec[i * 4 + 1] = PI / 2.0;         // std::atan2(1, 0);
        }
        else {
            evec[i * 4 + 0] = theta;
            evec[i * 4 + 1] = phi;
        }

        evec3Dpolar(&mats[i * 9], lambda[i * 3 + 2], theta, phi);

        if (isnan(theta) || isnan(phi)) {
            evec[i * 4 + 2] = PI / 2.0;         // acos(0)
            evec[i * 4 + 3] = 0;                // atan2(0, 0);
        }
        else {
            evec[i * 4 + 2] = theta;
            evec[i * 4 + 3] = phi;
        }
    }

    template<typename T>
    T* Eigenvectors3DPolar(T* mats, T* lambda, const size_t n) {
        T* gpu_mats;
        T* gpu_lambda;
        size_t mats_bytes = sizeof(T) * 9 * n;                              // required bytes for storing the tensor
        size_t evals_bytes = sizeof(T) * 3 * n;                             // required bytes for storing eigenvalues
        size_t evecs_bytes = sizeof(T) * 4 * n;                             // required bytes for storing 2 3D eigenvectors (in polar coordinates)

        // determine if the source volume is provided on the CPU or GPU
        cudaPointerAttributes attribs;										// create a pointer attribute structure
        HANDLE_ERROR(cudaPointerGetAttributes(&attribs, mats));			    // get the attributes for the source pointer

        if (attribs.type == cudaMemoryTypeDevice) {							// if the provided pointer is on the device
            gpu_mats = mats;									            // set the gpu_source pointer to source
            gpu_lambda = lambda;
        }
        else {																// otherwise copy the source image to the GPU
            HANDLE_ERROR(cudaMalloc(&gpu_mats, mats_bytes));								// allocate space on the GPU for the source image
            HANDLE_ERROR(cudaMemcpy(gpu_mats, mats, mats_bytes, cudaMemcpyHostToDevice));   // copy the source image to the GPU
            HANDLE_ERROR(cudaMalloc(&gpu_lambda, evals_bytes));
            HANDLE_ERROR(cudaMemcpy(gpu_lambda, lambda, evals_bytes, cudaMemcpyHostToDevice));
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
        HANDLE_ERROR(cudaMalloc(&gpu_evecs, evecs_bytes));
        kernel_evec3DPolar << <gridDim, blockDim >> > (gpu_mats, gpu_lambda, n, gpu_evecs);

        T* evecs;
        if (attribs.type == cudaMemoryTypeDevice) {
            evecs = gpu_evecs;
        }
        else {
            evecs = new T[4 * n];
            HANDLE_ERROR(cudaMemcpy(evecs, gpu_evecs, evecs_bytes, cudaMemcpyDeviceToHost));
            HANDLE_ERROR(cudaFree(gpu_evecs));
            HANDLE_ERROR(cudaFree(gpu_lambda));
            HANDLE_ERROR(cudaFree(gpu_mats));
        }

        return evecs;
       
    }
}