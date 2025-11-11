/*
	This file implements an GPU-based eigendecompositions for arrays of 2x2 and 3x3 matrices.
*/
#pragma once
#include <tira/eigen.h>
#include <tira/cuda/error.h>

namespace tira::cuda {

    template<typename T>
    __global__ void kernel_eval2_symmetric(const T* mats, const size_t n, T* evals) {
        const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) return;

        T a = mats[i * 4 + 0];
        T b = mats[i * 4 + 1];
        T c = mats[i * 4 + 3];
        eval2_symmetric(a, b, c, evals[i * 2 + 0], evals[i * 2 + 1]);
    }

    template<typename T>
    T* evals2_symmetric(const T* mats, const size_t n, int device) {

        if (device < 0)                                                     // if the device is < zero, run the CPU version
            return cpu::evals2_symmetric(mats, n);

        const T* gpu_mats;
        T* temp_gpu_mats;
        size_t mats_bytes = sizeof(T) * 4 * n;
        size_t evals_bytes = sizeof(T) * 2 * n;

        // determine if the source image is provided on the CPU or GPU
        cudaPointerAttributes attribs;										// create a pointer attribute structure
        HANDLE_ERROR(cudaPointerGetAttributes(&attribs, mats));			// get the attributes for the source pointer

        if (attribs.type == cudaMemoryTypeDevice) {							// if the provided pointer is on the device
            gpu_mats = mats;									            // set the gpu_source pointer to source
        }
        else {																// otherwise copy the source image to the GPU
            HANDLE_ERROR(cudaMalloc(&temp_gpu_mats, mats_bytes));								// allocate space on the GPU for the source image
            HANDLE_ERROR(cudaMemcpy(temp_gpu_mats, mats, mats_bytes, cudaMemcpyHostToDevice));// copy the source image to the GPU
            gpu_mats = static_cast<const T*>(temp_gpu_mats);
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
        kernel_eval2_symmetric << <gridDim, blockDim >> > (gpu_mats, n, gpu_evals);

        T* evals;
        if (attribs.type == cudaMemoryTypeDevice) {
            evals = gpu_evals;
        }
        else {
            evals = new T[2 * n];
            HANDLE_ERROR(cudaMemcpy(evals, gpu_evals, evals_bytes, cudaMemcpyDeviceToHost));
            HANDLE_ERROR(cudaFree(gpu_evals));
            HANDLE_ERROR(cudaFree(temp_gpu_mats));
        }

        return evals;
    }

    template<typename T>
    __global__ void kernel_evec2polar_symmetric(const T* mats, const T* evals, const size_t n, T* evecs) {
        const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) return;
        T a = mats[i * 4 + 0];
        T b = mats[i * 4 + 1];
        T c = mats[i * 4 + 3];
        evec2polar_symmetric(a, b, c, &evals[i * 2], evecs[i * 2 + 0], evecs[i * 2 + 1]);
    }

    template<typename T>
    T* evecs2polar_symmetric(const T* mats, T* evals, size_t n, int device) {

        const T* gpu_mats;
        const T* gpu_evals;
        T* temp_gpu_mats;
        T* temp_gpu_evals;
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
            HANDLE_ERROR(cudaMalloc(&temp_gpu_mats, mats_bytes));								// allocate space on the GPU for the source image
            HANDLE_ERROR(cudaMemcpy(temp_gpu_mats, mats, mats_bytes, cudaMemcpyHostToDevice));// copy the source image to the GPU
            HANDLE_ERROR(cudaMalloc(&temp_gpu_evals, ev_bytes));
            HANDLE_ERROR(cudaMemcpy(temp_gpu_evals, evals, ev_bytes, cudaMemcpyHostToDevice));
            gpu_mats = static_cast<const T*>(temp_gpu_mats);
            gpu_evals = static_cast<const T*>(temp_gpu_evals);
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
        kernel_evec2polar_symmetric << <gridDim, blockDim >> > (gpu_mats, gpu_evals, n, gpu_evecs);

        T* evecs;
        if (attribs.type == cudaMemoryTypeDevice) {
            evecs = gpu_evecs;
        }
        else {
            evecs = new T[2 * n];
            HANDLE_ERROR(cudaMemcpy(evecs, gpu_evecs, ev_bytes, cudaMemcpyDeviceToHost));
            HANDLE_ERROR(cudaFree(gpu_evecs));
            HANDLE_ERROR(cudaFree(temp_gpu_evals));
            HANDLE_ERROR(cudaFree(temp_gpu_mats));
        }

        return evecs;
    }


    // ----- Kernels and wrappers for 3D symmetric matrices -----
    /* Matrices are stored in column-major order per 3x3 block:
         | a  b  d |
         | b  c  e |
         | d  e  f |
       laid out as [a,b,d,  b,c,e,  d,e,f] for each matrix.
    */

    template<typename T>
    __global__ void kernel_eval3_symmetric(const T* mats, const size_t n, T* evals) {
        const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) return;

        const T a = mats[i * 9 + 0];
        const T b = mats[i * 9 + 1];
        const T d = mats[i * 9 + 2];
        const T c = mats[i * 9 + 4];
        const T e = mats[i * 9 + 5];
        const T f = mats[i * 9 + 8];

        T l0, l1, l2;
        tira::eval3_symmetric(a, b, c, d, e, f, l0, l1, l2);
        
        evals[i * 3 + 0] = l0;
        evals[i * 3 + 1] = l1;
        evals[i * 3 + 2] = l2;
    }

    template<typename T>
    T* evals3_symmetric(const T* mats, const size_t n, int device) {
        if (device < 0)                                                     // if the device is < zero, run the CPU version
            return cpu::evals3_symmetric(mats, n);

        // Set up sizes for GPU storage
        const size_t mats_bytes = sizeof(T) * 9 * n;                              // required bytes for storing the tensor
        const size_t evals_bytes = sizeof(T) * 3 * n;                             // required bytes for storing eigenvalues

		// Set up pointers
        const T* gpu_mats;
        T* temp_gpu_mats;

        // Determine if the source volume is provided on the CPU or GPU
        cudaPointerAttributes attribs;										// create a pointer attribute structure
        HANDLE_ERROR(cudaPointerGetAttributes(&attribs, mats));			    // get the attributes for the source pointer

        if (attribs.type == cudaMemoryTypeDevice) {							// if the provided pointer is on the device
            gpu_mats = mats;									            // set the gpu_source pointer to source
        }
        else {																// otherwise copy the source image to the GPU
            HANDLE_ERROR(cudaMalloc(&temp_gpu_mats, mats_bytes));								// allocate space on the GPU for the source image
            HANDLE_ERROR(cudaMemcpy(temp_gpu_mats, mats, mats_bytes, cudaMemcpyHostToDevice));   // copy the source image to the GPU
			gpu_mats = static_cast<const T*>(temp_gpu_mats);
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
        kernel_eval3_symmetric << <gridDim, blockDim >> > (gpu_mats, n, gpu_evals);

        if (attribs.type == cudaMemoryTypeDevice) {
            return gpu_evals;
        }
        else {
            T* evals = new T[3 * n];
            HANDLE_ERROR(cudaMemcpy(evals, gpu_evals, evals_bytes, cudaMemcpyDeviceToHost));
            HANDLE_ERROR(cudaFree(gpu_evals));
            HANDLE_ERROR(cudaFree(temp_gpu_mats));
            return evals;
        }
    }

    template<typename T>
    __global__ void kernel_evec3spherical_symmetric(const T* mats, const T* lambda, const size_t n, T* evecs) {
        const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) return;

        const T a = mats[i * 9 + 0];
        const T b = mats[i * 9 + 1];
        const T d = mats[i * 9 + 2];
        const T c = mats[i * 9 + 4];
        const T e = mats[i * 9 + 5];
        const T f = mats[i * 9 + 8];

        // Now we can safely calculate the eigenvectors
        const T evals[3] = { lambda[i * 3 + 0], lambda[i * 3 + 1], lambda[i * 3 + 2] };
        
        T evec0[3], evec1[3], evec2[3];
        tira::evec3_symmetric(a, b, c, d, e, f, evals, evec0, evec1, evec2);

        // Convert to spherical (theta, phi) coordinates for each eigenvector
        evecs[i * 6 + 0] = atan2(evec0[1], evec0[0]);                      // theta
        evecs[i * 6 + 1] = acos(evec0[2]);                                 // phi

        evecs[i * 6 + 2] = atan2(evec1[1], evec1[0]);                      // theta
        evecs[i * 6 + 3] = acos(evec1[2]);                                 // phi

        evecs[i * 6 + 4] = atan2(evec2[1], evec2[0]);                      // theta
        evecs[i * 6 + 5] = acos(evec2[2]);                                 // phi
    }
    
    template<typename T>
    T* evecs3spherical_symmetric(const T* mats, const T* lambda, const size_t n, int device) {
        if (device < 0)                                                     // if the device is < zero, run the CPU version
            return cpu::evecs3spherical_symmetric(mats, lambda, n);

        const T* gpu_mats;
        const T* gpu_lambda;
        T* temp_gpu_mats;
        T* temp_gpu_lambda;
        size_t mats_bytes = sizeof(T) * 9 * n;                              // required bytes for storing the tensor
        size_t evals_bytes = sizeof(T) * 3 * n;                             // required bytes for storing eigenvalues
        size_t evecs_bytes = sizeof(T) * 6 * n;                             // required bytes for storing 2 3D eigenvectors (in polar coordinates)

        // Determine if the source volume is provided on the CPU or GPU
        cudaPointerAttributes attribs;										// create a pointer attribute structure
        HANDLE_ERROR(cudaPointerGetAttributes(&attribs, mats));			    // get the attributes for the source pointer

        if (attribs.type == cudaMemoryTypeDevice) {							// if the provided pointer is on the device
            gpu_mats = mats;									            // set the gpu_source pointer to source
        }
        else {																// otherwise copy the source image to the GPU
            HANDLE_ERROR(cudaMalloc(&temp_gpu_mats, mats_bytes));								    // allocate space on the GPU for the source image
            HANDLE_ERROR(cudaMemcpy(temp_gpu_mats, mats, mats_bytes, cudaMemcpyHostToDevice));      // copy the source image to the GPU
			gpu_mats = static_cast<const T*>(temp_gpu_mats);
        }

        // Determine if the eigenvalues are provided on the CPU or GPU
        HANDLE_ERROR(cudaPointerGetAttributes(&attribs, lambda));			    // get the attributes for the source pointer
        if (attribs.type == cudaMemoryTypeDevice) {							    // if the provided pointer is on the device
            gpu_lambda = lambda;									            // set the gpu_source pointer to source
        }
        else {																    // otherwise copy the source image to the GPU
            HANDLE_ERROR(cudaMalloc(&temp_gpu_lambda, evals_bytes));
            HANDLE_ERROR(cudaMemcpy(temp_gpu_lambda, lambda, evals_bytes, cudaMemcpyHostToDevice));
            gpu_lambda = static_cast<const T*>(temp_gpu_lambda);
        }

        // Get the active device properties to calculate the optimal the block size
        HANDLE_ERROR(cudaGetDevice(&device));
        cudaDeviceProp props;
        HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
        unsigned int max_threads = props.maxThreadsPerBlock;
        dim3 blockDim = max_threads;
        dim3 gridDim = n / blockDim.x + 1;	                                    // calculate the grid size for the first pass

        T* gpu_evecs;
        HANDLE_ERROR(cudaMalloc(&gpu_evecs, evecs_bytes));
        kernel_evec3spherical_symmetric << <gridDim, blockDim >> > (gpu_mats, gpu_lambda, n, gpu_evecs);

        T* evecs;
        if (attribs.type == cudaMemoryTypeDevice) {
            evecs = gpu_evecs;
        }
        else {
            evecs = new T[6 * n];
            HANDLE_ERROR(cudaMemcpy(evecs, gpu_evecs, evecs_bytes, cudaMemcpyDeviceToHost));
            HANDLE_ERROR(cudaFree(gpu_evecs));
            HANDLE_ERROR(cudaFree(temp_gpu_lambda));
            HANDLE_ERROR(cudaFree(temp_gpu_mats));
        }

        return evecs;
    }
}