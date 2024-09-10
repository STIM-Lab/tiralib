// Eigen is used to solve the linear system required to calculate derivatives at any specified order
#include <iostream>

namespace tira {
	namespace solvers {

		template<typename T>
		void printmat(T* A, int n, int m) {
			for (int mi = 0; mi < m; mi++) {
				for (int ni = 0; ni < n; ni++) {
					std::cout << A[mi * m + ni] << "\t";
				}
				std::cout << std::endl;
			}
		}

		/// <summary>
		/// Simple implementation of the LAPACK gesv function:
		/// DGESV computes the solution to a real system of linear equations
		/// A* X = B,
		/// where A is an N - by - N matrix and X and B are N - by - NRHS matrices.

		/// The LU decomposition with partial pivoting and row interchanges is
		/// used to factor A as
		/// A = P * L * U,
		/// where P is a permutation matrix, L is unit lower triangular, and U is
		/// upper triangular.The factored form of A is then used to solve the
		/// system of equations A* X = B.
		/// </summary>
		/// <typeparam name="T">Data type for the calculation</typeparam>
		/// <param name="n">The number of linear equations, i.e., the order of the matrix A. N >= 0. </param>
		/// <param name="nrhs">The number of right hand sides, i.e., the number of columns of the matrix B. NRHS >= 0. </param>
		/// <param name="a">On entry, the N-by-N coefficient matrix A.</param>
		/// <param name="lda">The leading dimension of the array A.  LDA >= max(1,N)</param>
		/// <param name="ipiv">The pivot indices that define the permutation matrix P; row i of the matrix was interchanged with row IPIV(i) </param>
		/// <param name="b">On entry, the N-by-NRHS matrix of right hand side matrix B. On exit, if INFO = 0, the N - by - NRHS solution matrix X.</param>
		/// <param name="ldb">The leading dimension of the array B.  LDB >= max(1,N)</param>
		/// <param name="info">= 0:  successful exit; < 0:  if INFO = -i, the i - th argument had an illegal value; > 0:  if INFO = i, U(i, i) is exactly zero.The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed.</param>
		template<typename T>
		void gesv(int n, int nrhs, T* a, int lda, int* ipiv, T* b, int ldb, int info) {
			for (int k = 0; k < n - 1; k++) {					// for each pivot
				for (int i = k + 1; i < n; i++) {				// for each equation
					float m = a[i * n + k] / a[k * n + k];		// calculate scale 
					for (int j = k; j < n; j++) {				// for each coefficient
						a[i * n + j] -= m * a[k * n + j];		// subtract
					}
					b[i] -= m * b[k];							// forward elimination
				}
			}
			b[n - 1] = b[n - 1] / a[(n - 1) * n + n - 1];		// calculate
			for (int i = n - 2; i >= 0; i--) {					// for each equation
				float sum = b[i];								// start a sum
				for (int j = i + 1; j < n; j++) {				// for each coefficient
					sum -= a[i * n + j] * b[j];					// add coeffs and knowns
				}
				b[i] = sum / a[i * n + i];						// calculate
			}
		}

		/// <summary>
		/// Solve the linear system Ax=b
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="A">Pointer to an NxN matrix array</param>
		/// <param name="b">Pointer to an N element b vector</param>
		/// <param name="x">Pointer to memory allocated to store the N element solution vector x</param>
		/// <param name="N">Size of the linear system (number of equations/unknowns)</param>
		template<typename T>
		void Ax_b(T* A, T* b, T* x, int N) {

			std::vector<T> AA(A, A + (N * N));
			std::vector<T> bb(b, b + N);
			
			gesv<T>(N, 0, &AA[0], 0, 0, &bb[0], 0, 0);
			memcpy(x, &bb[0], N * sizeof(T));
		}

		/// <summary>
		/// Find the sum of products, recursive
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="nums		">Pointer to an array of N offsets</param>
		/// <param name="k		">Amount of elements in one product</param>
		/// <param name="index		">Index of the excluded element</param>
		/// <param name="sum		">Current sum</param>
		/// <param name="start		">Next index after the last added</param>
		/// <param name="currentProduct  ">Current product in the branch</param>
		/// <param name="selected	">Depth of the current branch</param>
		void combinationProduct(const std::vector<double>& nums, int k, int index, double& sum, int start = 0, double currentProduct = 1, int selected = 0) {
			if (selected == k) {			// leave if the length of product = k
				sum += currentProduct;
				return;
			}
			
			for (int i = start; i < nums.size(); ++i) {
				if (i == index) continue;  // Exclude the element at index j
				combinationProduct(nums, k, index, sum, i + 1, currentProduct * nums[i], selected + 1);
			}
		}

		/// <summary>
		/// Find the sum of products of all possible combinations of k elements, excluding nums[index]
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="nums	">Pointer to an N array  of offsets</param>
		/// <param name="k	">Amount of elements in one product</param>
		/// <param name="index	">Index of the excluded element</param>
		long int sumOfProducts(const std::vector<double>& nums, int k, int index) {
			if (k == 0)			// only 1 way to choose 0 elements
				return 1;
			double sum = 0;
			combinationProduct(nums, k, index, sum);
			return sum;
		}

		/// <summary>
		/// Finds the product of differences between elements in nums and nums[index]
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="nums	">Pointer to an N array  of offsets</param>
		/// <param name="index	">Index of the processed element</param>
		long int prodOfDifferences(const std::vector<double>& nums, unsigned int index) {
			long int prod = 1;
			for (int i = 0; i < nums.size(); i++) {
				if (i != index) {
					prod *= nums[i] - nums[index];
				}
			}
			return prod;
		}

		/// <summary>
		/// Solves Ax=b, where A is a Vandermonde Matrix, and b_derivative = derivative!
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="nums		">Pointer to an N array of geometric progression bases</param>
		/// <param name="coefs		">Pointer to the solution</param>
		/// <param name="derivative	">The calculated derivative</param>
		template<typename T>
		void finite_difference_vandermonde(const std::vector<T>& samples, std::vector<T>& coefs, unsigned int derivative) {

			std::vector<double> dbl_samples(samples.begin(), samples.end());
			std::vector<double> dbl_coef(samples.size());
			float sign = std::pow(-1.0f, derivative);
			long int fact = tgamma(derivative + 1);							// this should only be a problem if the derivative is large
			for (int j = 0; j < samples.size(); j++) {					// for each sample
				double sum = sumOfProducts(dbl_samples, dbl_samples.size() - derivative - 1, j);
				double prod = prodOfDifferences(dbl_samples, j);
				dbl_coef[j] = sign * sum / prod * fact;
			}

			coefs = std::vector<T>(dbl_coef.begin(), dbl_coef.end());
		}
	}
}