/*
	This is Diffusion.
*/

#pragma once
#include <tira/cuda/error.h>
#include <tira/cuda/callable.h>
#include <cuda_runtime.h>

namespace tira::cuda {

    // Kernel for computing the derivative
    template<typename T>
    __global__ void derivative_kernel(T* field, float* new_field, size_t* shape, size_t axis) {
        size_t row = blockIdx.y * blockDim.y + threadIdx.y;
        size_t col = blockIdx.x * blockDim.x + threadIdx.x;

        if (row >= shape[0] || col >= shape[1]) return;
        size_t idx = row * shape[1] + col;

        if (axis == 0) { // Derivative along rows
            if (row == 0) {
                new_field[idx] = 0.5f * field[(row + 1) * shape[1] + col];
            }
            else if (row == shape[0] - 1) {
                new_field[idx] = 0.5f * field[(row - 1) * shape[1] + col];
            }
            else {
                new_field[idx] = 0.5f * (field[(row - 1) * shape[1] + col] + field[(row + 1) * shape[1] + col]);
            }
        }
        else if (axis == 1) { // Derivative along columns
            if (col == 0) {
                new_field[idx] = 0.5f * field[row * shape[1] + (col + 1)];
            }
            else if (col == shape[1] - 1) {
                new_field[idx] = 0.5f * field[row * shape[1] + (col - 1)];
            }
            else {
                new_field[idx] = 0.5f * (field[row * shape[1] + (col - 1)] + field[row * shape[1] + (col + 1)]);
            }
        }
    }

    __global__ void mult(float* in1, float* in2, float* out, size_t* shape) {
        size_t row = blockIdx.y * blockDim.y + threadIdx.y;
        size_t col = blockIdx.x * blockDim.x + threadIdx.x;
        if (row < shape[0] && col < shape[1]) {
            size_t idx = row * shape[1] + col;
            out[idx] = in1[idx] * in2[idx];
        }
    }

    // kernel for adding 4 arrays
    __global__ void add4(float* in1, float* in2, float* in3, float* in4, float* out, size_t* shape) {
        size_t row = blockIdx.y * blockDim.y + threadIdx.y;
        size_t col = blockIdx.x * blockDim.x + threadIdx.x;
        if (row < shape[0] && col < shape[1]) {
            size_t idx = row * shape[1] + col;
            out[idx] = in1[idx] + in2[idx] + in3[idx] + in4[idx];
        }
    }

}
