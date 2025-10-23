#pragma once
#include "tensorvote.h"
#include "cuda/error.h"
#include "eigen.cuh"

#include <chrono>
#include <glm/glm.hpp>

#ifndef TV3_MAX_CONST_NB
#define TV3_MAX_CONST_NB 1536
#endif

namespace tira::tensorvote {
    /*
    *    --------------------------------
                2D Tensor Voting
    *    --------------------------------
    */
    __global__ static void global_stickvote2(glm::mat2* VT, glm::vec2* L, glm::vec2* V, glm::vec2 sigma, unsigned int power, float norm,
        int w, int s0, int s1) {

        int x0 = blockDim.x * blockIdx.x + threadIdx.x;                                       // get the x and y image coordinates for the current thread
        int x1 = blockDim.y * blockIdx.y + threadIdx.y;
        if (x0 >= s0 || x1 >= s1)                                                          // if not within bounds of image, return
            return;

        glm::mat2 Receiver = stickvote2(L, V, sigma, power, norm, w, s0, s1, glm::ivec2(x0, x1));
        VT[x0 * s1 + x1] += Receiver;
    }

    __global__ static void global_platevote2(glm::mat2* VT, glm::vec2* L, glm::vec2 sigma, unsigned int power,
        int w, int s0, int s1, unsigned samples) {

        int x0 = blockDim.x * blockIdx.x + threadIdx.x;                                       // get the x and y image coordinates for the current thread
        int x1 = blockDim.y * blockIdx.y + threadIdx.y;
        if (x0 >= s0 || x1 >= s1)                                                          // if not within bounds of image, return
            return;

        glm::mat2 Receiver = platevote2(L, sigma, w, s0, s1, glm::ivec2(x0, x1), samples);
        VT[x0 * s1 + x1] += Receiver;
    }

    static void tensorvote2_cuda(const float* input_field, float* output_field, unsigned int s0, unsigned int s1, float sigma, float sigma2,
        unsigned int w, unsigned int power, int device, bool STICK, bool PLATE, bool debug, unsigned samples) {

        auto start = std::chrono::high_resolution_clock::now();
        cudaDeviceProp props;
        HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
        auto end = std::chrono::high_resolution_clock::now();
        float t_deviceprops = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        int tensorFieldSize = 4 * s0 * s1;

        start = std::chrono::high_resolution_clock::now();
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
        if (STICK)
            global_stickvote2 << < blocks, threads >> > ((glm::mat2*)gpuOutputField, (glm::vec2*)gpuL, (glm::vec2*)gpuV, glm::vec2(sigma, sigma2), power, sn, w, s0, s1);
        if (PLATE)
            global_platevote2 << < blocks, threads >> > ((glm::mat2*)gpuOutputField, (glm::vec2*)gpuL, glm::vec2(sigma, sigma2), power, w, s0, s1, samples);
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


    /*
    *    --------------------------------
		        3D Tensor Voting
    *    --------------------------------
    */
    struct Neighbor3D_CUDA {
        int du, dv, dw;      // x2 (u), x1 (v), x0 (w) offsets
        float3 d;            // normalized direction
        float c1, c2;        // exp(-l2/s1^2), exp(-l2/s2^2)
    };
    __constant__ Neighbor3D_CUDA d_neighbors_const[TV3_MAX_CONST_NB];

	// Convert Neighbor3D to Neighbor3D_CUDA
    static inline Neighbor3D_CUDA convertNeighbor3D(const Neighbor3D& n) {
        Neighbor3D_CUDA nc;
        nc.du = n.du;   nc.dv = n.dv;   nc.dw = n.dw;
        nc.d = make_float3(n.d.x, n.d.y, n.d.z);
        nc.c1 = n.c1;   nc.c2 = n.c2;
        return nc;
	}

    // Host helper: upload neighbor table to constant or global memory
    struct DeviceNeighbors {
		Neighbor3D_CUDA* d_ptr = nullptr;           // null -> using constant memory
        int count = 0;
        bool used_const = false;
    };

    static DeviceNeighbors upload_neighbors(const std::vector<Neighbor3D>& NB) {
        DeviceNeighbors dn;
        dn.count = (int)NB.size();

        std::vector<Neighbor3D_CUDA> host(NB.size());
        for (size_t i = 0; i < NB.size(); ++i) host[i] = convertNeighbor3D(NB[i]);
        const size_t host_bytes = host.size() * sizeof(Neighbor3D_CUDA);
		// Upload to constant memory if it fits, otherwise to global memory
        if (host_bytes <= sizeof(d_neighbors_const)) {
            cudaMemcpyToSymbol(d_neighbors_const, host.data(), host_bytes, 0, cudaMemcpyHostToDevice);
            dn.used_const = true;
            dn.d_ptr = nullptr;
        }
        else {
            cudaMalloc(&dn.d_ptr, host_bytes);
            cudaMemcpy(dn.d_ptr, host.data(), host_bytes, cudaMemcpyHostToDevice);
            dn.used_const = false;
        }
        return dn;
    }

    static void free_neighbors(DeviceNeighbors& dn) {
        if (!dn.used_const && dn.d_ptr) {
            cudaFree(dn.d_ptr);
            dn.d_ptr = nullptr;
        }
        dn.count = 0;
        dn.used_const = false;
    }

    __global__ static void spherical_to_cart3_kernel(float* Qout,
        const float* V6, size_t N)
    {
        size_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= N) return;
		float theta = V6[i * 6 + 4];
        float phi = V6[i * 6 + 5];
        float st, ct; sincosf(theta, &st, &ct);
		float sp, cp; sincosf(phi, &sp, &cp);

        Qout[i * 3 + 0] = ct * sp;
		Qout[i * 3 + 1] = st * sp;
		Qout[i * 3 + 2] = cp;
    }

    __global__ static void global_stickvote3(glm::mat3* VT, const glm::vec3* L, const glm::vec3* Q, const Neighbor3D_CUDA* d_neighbors_glob,
        int nb_count, int usedConst, unsigned int power, float norm, int s0, int s1, int s2) {

        int x2 = blockDim.x * blockIdx.x + threadIdx.x;                                     // get the x, y, and z volume coordinates for the current thread
        int x1 = blockDim.y * blockIdx.y + threadIdx.y;
        int x0 = blockDim.z * blockIdx.z + threadIdx.z;
        if (x0 >= s0 || x1 >= s1 || x2 >= s2)                                               // if not within bounds of image, return
            return;
		const unsigned base_recv = (unsigned)x0 * s1 * s2 + (unsigned)x1 * s2 + (unsigned)x2;
        
		// Accumulator for the receiver voxel (symmetric 3x3 matrix)
		float m00 = 0.0f, m01 = 0.0f, m02 = 0.0f;
        float m11 = 0.0f, m12 = 0.0f, m22 = 0.0f;

        // Loop over neighbors
        for (int k = 0; k < nb_count; ++k) {
            const Neighbor3D_CUDA& nb = (usedConst ? d_neighbors_const[k] : d_neighbors_glob[k]);

            const int r0 = x0 + nb.dw;
			const int r1 = x1 + nb.dv;
            const int r2 = x2 + nb.du;
			if ((unsigned)r0 >= s0 || (unsigned)r1 >= s1 || (unsigned)r2 >= s2) continue;   // out of bounds

			const unsigned base_voter = (unsigned)r0 * s1 * s2 + (unsigned)r1 * s2 + (unsigned)r2;

            // Read eigenvalues at vote (only l1 and l2)
			const float l1 = L[base_voter].y;
			const float l2 = L[base_voter].z;
            float scale = std::copysignf(fabsf(l2) - fabsf(l1), l2);

            // Eigenvector at voter in spherical angles
            const glm::vec3 q = Q[base_voter];
            const float qx = q.x, qy = q.y, qz = q.z;

            // Dot(q, nb.d)
			const float qTd = qx * nb.d.x + qy * nb.d.y + qz * nb.d.z;
            const float qTd2 = qTd * qTd;

			// Compute the decay function
            float a = 1.0f - qTd2, b = qTd2;
            float ap = a, bp = b;
			for (unsigned i = 1; i < power; ++i) { ap *= a; bp *= b; }      // a^p, b^p
			float eta = (power == 1) 
                ? (nb.c1 * a) + (nb.c2 * b) 
                : (nb.c1 * ap) + (nb.c2 * bp);
			const float term = scale * eta;

            const float rx = qx - 2.0f * qTd * nb.d.x;
            const float ry = qy - 2.0f * qTd * nb.d.y;
            const float rz = qz - 2.0f * qTd * nb.d.z;

            // Accumulate outer product
			m00 += term * rx * rx; m01 += term * rx * ry; m02 += term * rx * rz;
            m11 += term * ry * ry; m12 += term * ry * rz; m22 += term * rz * rz;
        }

        // Write back to VT
		glm::mat3 Receiver;
        Receiver[0][0] = m00; Receiver[0][1] = m01; Receiver[0][2] = m02;
        Receiver[1][0] = m01; Receiver[1][1] = m11; Receiver[1][2] = m12;
        Receiver[2][0] = m02; Receiver[2][1] = m12; Receiver[2][2] = m22;
		VT[base_recv] += Receiver * norm;
    }

    __global__ static void global_platevote3(glm::mat3* VT, glm::vec3* L, glm::vec2 sigma, unsigned int power,
        int w, int s0, int s1, int s2, unsigned samples) {
	    return; // Not implemented yet
    }

    static void tensorvote3_cuda(const float* input_field, float* output_field, unsigned int s0, unsigned int s1, unsigned int s2, float sigma,
        float sigma2, unsigned int w, unsigned int power, int device, bool STICK, bool PLATE, bool debug, unsigned samples) {
        auto start = std::chrono::high_resolution_clock::now();
        cudaDeviceProp props;
        HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
        auto end = std::chrono::high_resolution_clock::now();
        float t_deviceprops = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

		const size_t n_voxels           = (size_t)s0 * s1 * s2;
        const size_t tensorField_bytes  = sizeof(float) * 9 * n_voxels;
        const size_t evals_bytes        = sizeof(float) * 3 * n_voxels;
        const size_t evecs_bytes        = sizeof(float) * 6 * n_voxels;

        // Build neighbor table on host
		start = std::chrono::high_resolution_clock::now();
        const glm::vec2 sig(sigma, sigma2);
		std::vector<Neighbor3D> NB = build_neighbors3d((int)w, sig);
        DeviceNeighbors d_nb = upload_neighbors(NB);
		end = std::chrono::high_resolution_clock::now();
        float t_buildnb = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        if (debug)
			std::cout << "Neighbor table: " << NB.size() << " entries, build+upload time: " << t_buildnb << " ms" << std::endl;

		// Eigendecomposition on GPU
        start = std::chrono::high_resolution_clock::now();
        float* L = tira::cuda::evals3_symmetric(input_field, n_voxels, device);
        float* V = tira::cuda::evecs3spherical_symmetric(input_field, L, n_voxels, device);
        end = std::chrono::high_resolution_clock::now();
        float t_eigendecomposition = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        // Declare GPU arrays
        float* gpuOutputField;
        float* gpuV;
        float* gpuL;

		// Check if input/eigens are already on device
        start = std::chrono::high_resolution_clock::now();
        cudaPointerAttributes attrIn, attrL, attrV;
        HANDLE_ERROR(cudaPointerGetAttributes(&attrL, L));
        HANDLE_ERROR(cudaPointerGetAttributes(&attrV, V));

        if (attrL.type == cudaMemoryTypeDevice) gpuL = L;
        else {
            HANDLE_ERROR(cudaMalloc(&gpuL, evals_bytes));
            HANDLE_ERROR(cudaMemcpy(gpuL, L, evals_bytes, cudaMemcpyHostToDevice));
        }
        if (attrV.type == cudaMemoryTypeDevice) gpuV = V;
        else {
            HANDLE_ERROR(cudaMalloc(&gpuV, evecs_bytes));
            HANDLE_ERROR(cudaMemcpy(gpuV, V, evecs_bytes, cudaMemcpyHostToDevice));
        }
		end = std::chrono::high_resolution_clock::now();
		float t_host2device = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        // Convert spherical angles to cartesian Q
        float* dQ = nullptr;
        HANDLE_ERROR(cudaMalloc(&dQ, 3 * n_voxels * sizeof(float)));
        {
            int t = 256; int g = (int)((n_voxels + t - 1) / t);
            spherical_to_cart3_kernel << <g, t >> > (dQ, (const float*)gpuV, n_voxels);
            HANDLE_ERROR(cudaDeviceSynchronize());
        }

        // Allocate output on device
        start = std::chrono::high_resolution_clock::now();
        HANDLE_ERROR(cudaMalloc(&gpuOutputField, tensorField_bytes));
        HANDLE_ERROR(cudaMemset(gpuOutputField, 0, tensorField_bytes));
        end = std::chrono::high_resolution_clock::now();
        float t_devicealloc = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        // Specify the CUDA block and grid dimensions
        dim3 threads(32, 4, 2);
        dim3 blocks(
            (unsigned)((s2 + threads.x - 1) / threads.x),
            (unsigned)((s1 + threads.y - 1) / threads.y),
            (unsigned)((s0 + threads.z - 1) / threads.z)
        );
		float sn = 1.0f / sticknorm3(sigma, sigma2, power);

        start = std::chrono::high_resolution_clock::now();
        if (STICK)
            global_stickvote3 << <blocks, threads >> > ((glm::mat3*)gpuOutputField, (const glm::vec3*)gpuL, (const glm::vec3*)dQ, 
                d_nb.d_ptr, d_nb.count, d_nb.used_const ? 1 : 0, power, sn, (int)s0, (int)s1, (int)s2);
        if (PLATE)
            global_platevote3 << <blocks, threads >> > ((glm::mat3*)gpuOutputField, (glm::vec3*)gpuL, sig, power, sn, w, s0, s1, s2);
        cudaDeviceSynchronize();
        end = std::chrono::high_resolution_clock::now();
        float t_voting = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

		free_neighbors(d_nb);

        start = std::chrono::high_resolution_clock::now();
        // Copy the final result back from the GPU
		cudaPointerAttributes attrOut;
		HANDLE_ERROR(cudaPointerGetAttributes(&attrOut, output_field));
        if (attrOut.type == cudaMemoryTypeDevice) 
            HANDLE_ERROR(cudaMemcpy(output_field, gpuOutputField, tensorField_bytes, cudaMemcpyDeviceToDevice));
        else
			HANDLE_ERROR(cudaMemcpy(output_field, gpuOutputField, tensorField_bytes, cudaMemcpyDeviceToHost));
        
        cudaDeviceSynchronize();
        end = std::chrono::high_resolution_clock::now();
        float t_device2host = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        // Free all of the GPU arrays
        start = std::chrono::high_resolution_clock::now();
        HANDLE_ERROR(cudaFree(gpuOutputField));
        if (attrL.type == cudaMemoryTypeDevice) HANDLE_ERROR(cudaFree(L));
        else delete[] L;
        if (attrV.type == cudaMemoryTypeDevice) HANDLE_ERROR(cudaFree(V));
		else delete[] V;
		if (gpuL && gpuL != L) HANDLE_ERROR(cudaFree(gpuL));
		if (gpuV && gpuV != V) HANDLE_ERROR(cudaFree(gpuV));
        cudaDeviceSynchronize();
        end = std::chrono::high_resolution_clock::now();
        float t_devicefree = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        if (debug) {
            //std::cout << "Eigendecomposition:  " << t_eigendecomposition << " ms" << std::endl;
            std::cout << "Voting: " << t_voting << " ms" << std::endl;
            std::cout << "cudaMemcpy (H->D):  " << t_host2device << " ms" << std::endl;
            std::cout << "cudaMemcpy (D->H):  " << t_device2host << " ms" << std::endl;
            std::cout << "cudaMalloc: " << t_devicealloc << " ms" << std::endl;
            std::cout << "cudaFree: " << t_devicefree << " ms" << std::endl;
            std::cout << "cudaDeviceProps: " << t_deviceprops << " ms" << std::endl;
        }
    }
}