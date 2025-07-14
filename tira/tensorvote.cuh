#pragma once
#include "tensorvote.h"
#include "cuda/error.h"
#include "eigen.cuh"

#include <chrono>
#include <glm/glm.hpp>

__global__ static void cuda_stickvote2(glm::mat2* VT, glm::vec2* L, glm::vec2* V, glm::vec2 sigma, unsigned int power, float norm,
    int w, int s0, int s1) {

    int x0 = blockDim.x * blockIdx.x + threadIdx.x;                                       // get the x and y image coordinates for the current thread
    int x1 = blockDim.y * blockIdx.y + threadIdx.y;
    if (x0 >= s0 || x1 >= s1)                                                          // if not within bounds of image, return
        return;

    glm::mat2 Receiver = stickvote2(L, V, sigma, power, norm, w, s0, s1, glm::ivec2(x0, x1));
    VT[x0 * s1 + x1] += Receiver;
}

__global__ static void cuda_platevote2(glm::mat2* VT, glm::vec2* L, glm::vec2 sigma, unsigned int power,
    int w, int s0, int s1, unsigned samples) {

    int x0 = blockDim.x * blockIdx.x + threadIdx.x;                                       // get the x and y image coordinates for the current thread
    int x1 = blockDim.y * blockIdx.y + threadIdx.y;
    if (x0 >= s0 || x1 >= s1)                                                          // if not within bounds of image, return
        return;

    glm::mat2 Receiver = platevote2(L, sigma, w, s0, s1, glm::ivec2(x0, x1), samples);
    VT[x0 * s1 + x1] += Receiver;
}

namespace tira::cuda {

    static void tensorvote2(const float* input_field, float* output_field, unsigned int s0, unsigned int s1, float sigma, float sigma2,
        unsigned int w, unsigned int power, int device, bool STICK, bool PLATE, bool debug, unsigned samples) {

        auto start = std::chrono::high_resolution_clock::now();
        cudaDeviceProp props;
        HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
        auto end = std::chrono::high_resolution_clock::now();
        float t_deviceprops = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        int tensorFieldSize = 4 * s0 * s1;

        start = std::chrono::high_resolution_clock::now();
        //float* L = EigenValues2(input_field, s0 * s1, device);
        float* L = tira::cuda::evals2_symmetric<float>(input_field, s0 * s1, device);
        float* V = tira::cuda::evecs2polar_symmetric(input_field, L, s0 * s1, device);
        end = std::chrono::high_resolution_clock::now();
        float t_eigendecomposition = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        // Declare GPU arrays
        float* gpuOutputField;
        float* gpuV;
        float* gpuL;


        // Allocate GPU arrays
        start = std::chrono::high_resolution_clock::now();
        HANDLE_ERROR(cudaMalloc(&gpuOutputField, tensorFieldSize * sizeof(float)));
        HANDLE_ERROR(cudaMemset(gpuOutputField, 0, tensorFieldSize * sizeof(float)));
        HANDLE_ERROR(cudaMalloc(&gpuV, s0 * s1 * 2 * sizeof(float)));
        HANDLE_ERROR(cudaMalloc(&gpuL, s0 * s1 * 2 * sizeof(float)));
        cudaDeviceSynchronize();
        end = std::chrono::high_resolution_clock::now();
        float t_devicealloc = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        start = std::chrono::high_resolution_clock::now();
        // Copy input arrays
        HANDLE_ERROR(cudaMemcpy(gpuV, V, s0 * s1 * 2 * sizeof(float), cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(gpuL, L, s0 * s1 * 2 * sizeof(float), cudaMemcpyHostToDevice));
        cudaDeviceSynchronize();
        end = std::chrono::high_resolution_clock::now();

        float t_host2device = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        // Specify the CUDA block and grid dimensions
        size_t blockDim = sqrt(props.maxThreadsPerBlock);
        dim3 threads(blockDim, blockDim);
        dim3 blocks(s0 / threads.x + 1, s1 / threads.y + 1);

        float sn = 1.0 / sticknorm2(sigma, sigma2, power);

        if (debug)
            std::cout << "Stick Area: " << sn << std::endl;

        start = std::chrono::high_resolution_clock::now();
        if(STICK)
            cuda_stickvote2 << < blocks, threads >> > ((glm::mat2*)gpuOutputField, (glm::vec2*)gpuL, (glm::vec2*)gpuV, glm::vec2(sigma, sigma2), power, sn, w, s0, s1);
        if (PLATE)
            cuda_platevote2 <<< blocks, threads >>>((glm::mat2*)gpuOutputField, (glm::vec2*)gpuL, glm::vec2(sigma, sigma2), power, w, s0, s1, samples);
        cudaDeviceSynchronize();
        end = std::chrono::high_resolution_clock::now();
        float t_voting = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();


        start = std::chrono::high_resolution_clock::now();
        // Copy the final result back from the GPU
        HANDLE_ERROR(cudaMemcpy(output_field, gpuOutputField, tensorFieldSize * sizeof(float), cudaMemcpyDeviceToHost));
        cudaDeviceSynchronize();
        end = std::chrono::high_resolution_clock::now();
        float t_device2host = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        // Free all of the GPU arrays
        start = std::chrono::high_resolution_clock::now();
        HANDLE_ERROR(cudaFree(gpuOutputField));
        HANDLE_ERROR(cudaFree(gpuV));
        HANDLE_ERROR(cudaFree(gpuL));
        cudaDeviceSynchronize();
        end = std::chrono::high_resolution_clock::now();
        float t_devicefree = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        if (debug) {
            std::cout << "Eigendecomposition:  " << t_eigendecomposition << " ms" << std::endl;
            std::cout << "Voting: " << t_voting << " ms" << std::endl;
            std::cout << "cudaMemcpy (H->D):  " << t_host2device << " ms" << std::endl;
            std::cout << "cudaMemcpy (D->H):  " << t_device2host << " ms" << std::endl;
            std::cout << "cudaMalloc: " << t_devicealloc << " ms" << std::endl;
            std::cout << "cudaFree: " << t_devicefree << " ms" << std::endl;
            std::cout << "cudaDeviceProps: " << t_deviceprops << " ms" << std::endl;
        }
    }



}