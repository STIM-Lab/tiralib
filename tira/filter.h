#pragma once

#include <numbers>
#include <vector>
#include<memory>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace tira::cpu {

	template <typename T>
	static T normaldist(T mu, T sigma, T x) {
		T sigma_sq = sigma * sigma;
		T pi = static_cast<T>(M_PI);
		T n = static_cast<T>(1) / std::sqrt(static_cast<T>(2) * pi * sigma_sq);

		T exponent = -((x - mu) * (x - mu)) / (static_cast<T>(2) * sigma_sq);
		return n * exp(exponent);
	}
	
	/**
	 * @brief Generate a 1D Gaussian kernel sampled on a regular grid.
	 *
	 * @tparam T         Floating-point type for kernel values.
	 * @param ksize		 Number of samples in the kernel.
	 * @param mu         Mean of the Gaussian (typically 0).
	 * @param sigma      Standard deviation of the Gaussian.
	 * @param dx         Sample spacing of the grid.
	 * @return std::vector<T> 1D Gaussian kernel of length @p windowSize.
	 */
	template <typename T>
	static std::vector<T> kernel_gaussian(unsigned int ksize, T mu, T sigma, T dx) {
		std::vector<T> kernel(ksize);

		// Physical length of the kernel
		const T halfwidth = static_cast<T>(ksize) * dx / static_cast<T>(2);

		// First sample is centered at -halfwidth + dx/2 so the grid is symmetric around mu
		const T startx = -halfwidth + (dx / static_cast<T>(2));

		// Sample the Gaussian at regularly spaced positions
		T sum = static_cast<T>(0);
		for (unsigned int xi = 0; xi < ksize; ++xi) {
			const T x = startx + static_cast<T>(xi) * dx;
			kernel[xi] = normaldist(mu, sigma, x);
			sum += kernel[xi];
		}

		// Normalize the kernel to have sum 1
		if (sum > static_cast<T>(0)) {
			T inv_sum = static_cast<T>(1) / sum;
			for (auto& v : kernel) v *= inv_sum;
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
};