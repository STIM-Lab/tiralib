#pragma once

#include "tira/solvers.h"

namespace tira::calculus {
	inline std::vector<double> central_difference_coefficients(const unsigned int derivative, unsigned int order) {
		if ((derivative + order) % 2 == 0) order += 1;

		const unsigned int N = order + derivative;		// calculate the number of samples required to achieve the desired order

		std::vector<double> Coefficients;

		std::vector<double> Samples(N);					// allocate a vector that will be used to store sample points

		const int ri = N / 2;
		for (unsigned int ci = 0; ci < N; ci++) {		// calculate the point for each sample
			Samples[ci] = -static_cast<double>(ri) + static_cast<double>(ci);				// store that point in the Samples vector
		}
		/*std::vector<T> A(N * N);
		for (unsigned int ri = 0; ri < N; ri++) {
			for (unsigned int ci = 0; ci < N; ci++) {
				A[ri * N + ci] = pow(Samples[ci], ri);
			}
		}
		std::vector<T> b(N, 0.0);
		b[derivative] = tgamma(derivative + 1);
		std::vector<T> x(N);
		tira::solvers::Ax_b(&A[0], &b[0], &x[0], N);*/

		std::vector<double> x(N);
		tira::solvers::finite_difference_vandermonde(Samples, Coefficients, derivative);

		return Coefficients;
	}

	/// <summary>
	/// Calculate the finite difference coefficients given a derivative and order of accuracy
	/// </summary>
	/// <param name="derivative"></param>
	/// <param name="order"></param>
	/// <returns></returns>
	inline std::vector< std::vector<double> > finite_difference_coefficients(const unsigned int derivative, const unsigned int order) {

		const unsigned int N = order + derivative;		// calculate the number of samples required to achieve the desired order

		std::vector< std::vector<double> > Coefficients;

		std::vector<double> Samples(N);					// allocate a vector that will be used to store sample points

		for (unsigned int ri = 0; ri < N; ri++) {			// for each shifted sample position
			for (unsigned int ci = 0; ci < N; ci++) {		// calculate the point for each sample
				Samples[ci] = -static_cast<double>(ri) + static_cast<double>(ci);				// store that point in the Samples vector
			}

			std::vector<double> x(N);
			tira::solvers::finite_difference_vandermonde(Samples, x, derivative);

			Coefficients.push_back(x);
		}
		return Coefficients;
	}

	inline void printCoefficients(const std::vector<double> Coefficients) {

		for (unsigned int i = 0; i < Coefficients.size(); i++) {
			std::cout << Coefficients[i] << " ";
		}
		std::cout << std::endl;
	}

	inline void printCoefficients(const std::vector< std::vector<double> > C) {
		for (unsigned int i = 0; i < C.size(); i++) {
			printCoefficients(C[i]);
		}
	}
}