#pragma once

#include <numbers>

namespace tira::cpu {

	template <typename T>
	static T normaldist(T mu, T sigma, T x) {
		
		T sigma_sq = sigma * sigma;
		T pi = (T)3.14159265358979323846;
		T n = 1 / sqrt(2 * pi * sigma_sq);

		T exponent = -pow(x - mu, 2) / (2 * sigma_sq);
		return n * exp(exponent);
	}


	template <typename T>
	static T* kernel_gaussian(unsigned WindowSize, T mu, T sigma, T dx) {
		T* g = (T*)malloc(WindowSize * sizeof(T));

		T sx = WindowSize * dx;
		T halfwidth = sx / 2;

		T startx = -halfwidth + (dx / 2);

		for (unsigned xi = 0; xi < WindowSize; xi++) {
			T x = startx + xi * dx;
			g[xi] = normaldist(mu, sigma, x);
		}
		return g;
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
};