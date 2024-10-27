/*
	This file implements an GPU-based eigendecompositions for arrays of 2x2 matrices.
*/

#pragma once
#include <tira/cuda/error.h>
#include <tira/cuda/callable.h>

namespace tira::cuda {
    #define PI 3.14159265358979323846

    template<typename T>
    CUDA_CALLABLE void eval2D(const T* matrix, T& eval0, T& eval1) {
        T d = matrix[0];
        T e = matrix[1];
        T f = matrix[2];
        T g = matrix[3];

        const float dpg = d + g;
        const float disc = sqrt((4 * e * f) + pow(d - g, 2));
        float a = (dpg + disc) / 2.0f;
        float b = (dpg - disc) / 2.0f;

        eval0 = std::abs(a) < std::abs(b) ? a : b;
        eval1 = std::abs(a) > std::abs(b) ? a : b;
    }

    template<typename T>
    __global__ void kernel_eval2D(T* mats, size_t n, T* evals) {
        const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) return;
        eval2D(&mats[i * 4], evals[i * 2 + 0], evals[i * 2 + 1]);
    }


    /// <summary>
    /// CPU code for calculating eigenvalues of a 2D matrix array
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <param name="mats"></param>
    /// <param name="n"></param>
    /// <returns></returns>
    template<typename T>
    T* cpuEigenvalues2D(const T* mats, const size_t n) {

        T* evals = new T[2*n];
        T eval0, eval1;
        for (size_t i = 0; i < n; i++) {
            eval2D(&mats[i * 4], eval0, eval1);
            evals[i * 2 + 0] = eval0;
            evals[i * 2 + 1] = eval1;
        }
        return evals;
    }

    template<typename T>
    T* Eigenvalues2D(T* mats, size_t n, int device) {

        if (device < 0)                                                     // if the device is < zero, run the CPU version
            return cpuEigenvalues2D(mats, n);

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

    // small then large
    template<typename T>
    CUDA_CALLABLE void evec2Dpolar(const T* matrix, const T* lambdas, T& theta0, T& theta1) {
        //[a  b]
        //[c  d]

        const float a = matrix[0];
        const float c = matrix[1];
        const float b = matrix[2];
        const float d = matrix[3];

        if (c != 0) {
            theta0 = std::atan2(c, lambdas[0] - d);
            theta1 = std::atan2(c, lambdas[1] - d);
        }
        else {
            theta0 = std::atan2(lambdas[0] - a, b);
            theta1 = std::atan2(lambdas[1] - a, b);
        }
    }

    template<typename T>
    __global__ void kernel_evec2Dpolar(T* mats, T* evals, size_t n, T* evecs) {
        const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) return;
        evec2Dpolar(&mats[i * 4], &evals[i * 2], evecs[i * 2 + 0], evecs[i * 2 + 1]);
    }

    template<typename T>
    T* cpuEigenvectors2DPolar(const T* mats, const T* evals, const size_t n) {

        T* vecs = new T[2 * n];
        T vec0, vec1;
        for (size_t i = 0; i < n; i++) {
            evec2Dpolar(&mats[i * 4], &evals[i * 2], vec0, vec1);
            vecs[i * 2 + 0] = vec0;
            vecs[i * 2 + 1] = vec1;
        }
        return vecs;
    }

    template<typename T>
    T* Eigenvectors2DPolar(T* mats, T* evals, size_t n, int device) {

        if (device < 0)                                                     // if the device is < zero, run the CPU version
            return cpuEigenvectors2DPolar(mats, evals, n);

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
    CUDA_CALLABLE void eval3D(const T* A, T& eval0, T& eval1, T& eval2) {
        
        T a, b, c, d, e, f, g, h, i;
        a = A[0 * 3 + 0]; b = A[0 * 3 + 1]; c = A[0 * 3 + 2];
        d = A[1 * 3 + 0]; e = A[1 * 3 + 1]; f = A[1 * 3 + 2];
        g = A[2 * 3 + 0]; h = A[2 * 3 + 1]; i = A[2 * 3 + 2];


        // Case: matrix is diagonal
        T p1 = b * b + c * c + f * f;
        if (p1 == 0.0) {
            eval0 = a;
            eval1 = e;
            eval2 = i;
            return;
        }

        T q = (a + e + i) / 3.0;
        T p2 = (a - q) * (a - q) + (e - q) * (e - q) + (i - q) * (i - q) + 2.0 * p1;
        T p = sqrt(p2 / 6.0);
        T pinv = 1.0 / p;
        T B[3][3];
        B[0][0] = pinv * (a - q);
        B[1][1] = pinv * (e - q);
        B[2][2] = pinv * (i - q);
        B[0][1] = pinv * b;
        B[0][2] = pinv * c;
        B[1][0] = pinv * d;
        B[1][2] = pinv * f;
        B[2][0] = pinv * g;
        B[2][1] = pinv * h;
        T detB = B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1]) -
            B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0]) +
            B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0]);
        T r = detB / 2.0;

        // In exact arithmetic for a symmetric matrix - 1 <= r <= 1
        // but computation error can leave it slightly outside this range.

        T phi = 0.0;
        if (r <= -1.0)
            phi = PI / 3.0;
        else if (r > 1.0)
            phi = 0.0;
        else
            phi = acos(r) / 3.0;

        // The eigenvalues satisfy l[0] <= l[1] <= l[2]
        eval2 = q + 2.0 * p * cos(phi);
        eval0 = q + 2.0 * p * cos(phi + (2.0 * PI / 3.0));
        eval1 = 3.0 * q - eval2 - eval0;  // since trace(A) = eig1 + eig2 + eig3
    }

    /// <summary>
    /// CPU code for calculating eigenvalues of a 2D matrix array
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <param name="mats"></param>
    /// <param name="n"></param>
    /// <returns></returns>
    template<typename T>
    T* cpuEigenvalues3D(const T* mats, const size_t n) {

        T* evals = new T[3 * n];
        T eval0, eval1, eval2;
        for (size_t i = 0; i < n; i++) {
            eval3D(&mats[i * 9], eval0, eval1, eval2);
            evals[i * 3 + 0] = eval0;
            evals[i * 3 + 1] = eval1;
            evals[i * 3 + 1] = eval2;
        }
        return evals;
    }

    template<typename T>
    T* Eigenvalues3D(T* mats, size_t n, int device) {
        if (device < 0)                                                     // if the device is < zero, run the CPU version
            return cpuEigenvalues3D(mats, n);
    }
}