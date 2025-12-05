#pragma once

#include <cmath>
#include <tira/cuda/callable.h>

#ifdef __CUDACC__
    #include <tira/cuda/error.h>
#endif

namespace tira::constant {
	constexpr double PI = 3.14159265358979323846;
    template <typename T>
	constexpr T TIRA_EIGEN_EPSILON = static_cast<T>(1e-12);
}

namespace tira::shared {
    template <typename T>
    CUDA_CALLABLE void swap(T& a, T& b) {
        T temp = a;
        a = b;
        b = temp;
    }

    template <typename T>
    CUDA_CALLABLE T sgn(const T& a) {
		return (a > T(0)) - (a < T(0));
    }
    
    template <typename T>
    CUDA_CALLABLE T abs(T a) {
        return (a < (T)0) ? -a : a;
	}

    template <typename T>
    CUDA_CALLABLE T dot3(const T* a, const T* b) {
        // | a0  a1  a2 |
        // | b0  b1  b2 |
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}
    
    template <typename T>
    CUDA_CALLABLE T determinant3(const T a, const T b, const T c, const T d,
        const T e, const T f) {
	    // | a  b  d |
	    // | b  c  e |
	    // | d  e  f |
	    return a * (c * f - e * e) - b * (b * f - e * d) + d * (b * e - c * d);
    }
    
    template <typename T>
    CUDA_CALLABLE void normalize3(T* v) {
        // | v0  v1  v2 |
        T n = sqrt(dot3(v, v));
        if (n > tira::constant::TIRA_EIGEN_EPSILON<T>) {
			T inv_n = T(1) / n;
            v[0] *= inv_n; 
			v[1] *= inv_n;
			v[2] *= inv_n;
        }
	}

    template <typename T>
    CUDA_CALLABLE void cross3(const T* a, const T* b, T* result) {
	    // | a0  a1  a2 |
	    // | b0  b1  b2 |
	    // result = | a1*b2 - a2*b1, a2*b0 - a0*b2, a0*b1 - a1*b0 |
	    result[0] = a[1] * b[2] - a[2] * b[1];
	    result[1] = a[2] * b[0] - a[0] * b[2];
	    result[2] = a[0] * b[1] - a[1] * b[0];
    }

    /// <summary>
    /// Calculate the roots of a quadratic equation where a=1 using stable
    /// assumptions for eigenvalues of 2x2 matrices. In particular, if b = 0
    /// this function will set both roots (eigenvalues) to zero. If the discriminant
    /// is negative, we assume that the negative is round-off error and set the
    /// discriminant to zero.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <param name="b"></param>
    /// <param name="c"></param>
    /// <returns></returns>
    template <typename T>
    CUDA_CALLABLE void quad_root(T b, T c, T& r1, T& r2) {
        T disc = b * b - 4 * c;
        if (b == 0) {
            r1 = 0;
            r2 = 0;
            return;
        }
        if (disc < 0) disc = 0;
        T q = -T(0.5) * (b + sgn(b) * sqrt(disc));
        r1 = q;
        r2 = c / q;
    }

    /// Calculate the eigenvalues of a 2x2 matrix
    template<typename T>
    CUDA_CALLABLE void eval2_symmetric(const T a, const T b, const T c, T& eval0, T& eval1) {
        // | a  b |
        // | b  c |

        // this algorithm uses the trace method to calculate the eigenvalues
        T tr = a + c;                                   // calculate the matrix trace
        T det = a * c - b * b;                          // calculate the determinant

        T l0, l1;
        quad_root(-tr, det, l0, l1);    // find the roots of the quadratic equation - these are the eigenvalues

        eval0 = abs(l0) < abs(l1) ? l0 : l1;  // sort the eigenvalues based on their magnitude
        eval1 = abs(l0) > abs(l1) ? l0 : l1;
    }

    /// Calculate the eigenvectors of a 2x2 matrix associated with the two eigenvalues stored in lambdas.
    /// The eigenvectors are returned in polar coordinates (theta). The input matrices are assumed to be in
    /// column-major format (similar to OpenGL).
    template<typename T>
    CUDA_CALLABLE void evec2polar_symmetric(const T a, const T b, const T c, const T* lambdas, T& theta0, T& theta1) {
        //[a  b]
        //[b  c]

        T l0 = lambdas[0];
        T l1 = lambdas[1];

        T a_l0 = a - l0;
        T a_l1 = a - l1;

        if (abs(a_l0) >= abs(a_l1)) {
            theta1 = atan2(b, a_l0);
            theta0 = theta1 + (tira::constant::PI / 2.0);
            if (theta0 > tira::constant::PI) theta0 -= 2 * tira::constant::PI;
        }
        else {
            theta0 = atan2(b, a_l1);
            theta1 = theta0 + (tira::constant::PI / 2.0);
            if (theta1 > tira::constant::PI) theta1 -= 2 * tira::constant::PI;
        }

    }


    /// Compute a right-handed orthonormal set { U, V, evec }.
    /// The vector evec is guaranteed to be unit-length, in which case
    /// there is no need to worry about a division by zero when computing invLength.
    template<typename T>
    CUDA_CALLABLE void ComputeOrthogonalComplement(const T* evec, T* U, T* V) {
        // Compute a right-handed orthonormal set { U, V, evec }
		const T absX = abs(evec[0]);
		const T absY = abs(evec[1]);
		const T absZ = abs(evec[2]);

        T temp[3];
        if (absY - absX > tira::constant::TIRA_EIGEN_EPSILON<T>) {
            if (absZ - absX > tira::constant::TIRA_EIGEN_EPSILON<T>) {
                temp[0] = T(1); temp[1] = T(0); temp[2] = T(0);         // X is smallest
            }
            else {
            temp[0] = T(0); temp[1] = T(0); temp[2] = T(1);         // Z is smallest
            }
        }
        else {
            if (absZ - absY > tira::constant::TIRA_EIGEN_EPSILON<T>) {
                temp[0] = T(0); temp[1] = T(1); temp[2] = T(0);          // Y is smallest
            }
            else {
                temp[0] = T(0); temp[1] = T(0); temp[2] = T(1);         // Z is smallest
            }
        }
		cross3(evec, temp, U);  // U = evec x temp
        normalize3(U);

		cross3(evec, U, V);     // V = evec x U
    }

    template <typename T>
    CUDA_CALLABLE void evec3_symmetric_0(const T a, const T b, const T c, const T d, const T e, const T f,
        const T eval0, T* evec0) {
        // | a  b  d |
	    // | b  c  e |
	    // | d  e  f |

        // M = A - eval0 * I
        T m00 = a - eval0; T m01 = b;       T m02 = d;
        T m10 = b;       T m11 = c - eval0; T m12 = e;
        T m20 = d;       T m21 = e;       T m22 = f - eval0;

        // Compute cofactors of M
        T c00 = m11 * m22 - m12 * m21;
        T c01 = m12 * m20 - m10 * m22;
        T c02 = m10 * m21 - m11 * m20;
        T c10 = m02 * m21 - m01 * m22;
        T c11 = m00 * m22 - m02 * m20;
        T c12 = m01 * m20 - m00 * m21;
        T c20 = m01 * m12 - m02 * m11;
        T c21 = m02 * m10 - m00 * m12;
        T c22 = m00 * m11 - m01 * m10;

        // Find the row of cofactors with the largest L2 norm
        T r0sqr = c00 * c00 + c01 * c01 + c02 * c02;
        T r1sqr = c10 * c10 + c11 * c11 + c12 * c12;
        T r2sqr = c20 * c20 + c21 * c21 + c22 * c22;

        if (r0sqr >= r1sqr && r0sqr >= r2sqr)
        {
            evec0[0] = c00; evec0[1] = c01; evec0[2] = c02;
        }
        else if (r1sqr >= r0sqr && r1sqr >= r2sqr)
        {
            evec0[0] = c10; evec0[1] = c11; evec0[2] = c12;
        }
        else
        {
            evec0[0] = c20; evec0[1] = c21; evec0[2] = c22;
        }
        normalize3(evec0); // ensure that the eigenvector is unit-length
    }

    template <typename T>
    CUDA_CALLABLE void evec3_symmetric_1(const T a, const T b, const T c, const T d, const T e, const T f,
        const T eval1, const T* evec0, T* evec1) {

        // Compute a right-handed orthonormal set { U, V, evec0 }
        T U[3]; T V[3];
	    ComputeOrthogonalComplement(evec0, U, V);

        // AU, AV
	    const T AU[3] = {   a * U[0] + b * U[1] + d * U[2],
		                    b * U[0] + c * U[1] + e * U[2],
		                    d * U[0] + e * U[1] + f * U[2]  };
	    const T AV[3] = {   a * V[0] + b * V[1] + d * V[2],
		                    b * V[0] + c * V[1] + e * V[2],
		                    d * V[0] + e * V[1] + f * V[2]  };
	    
        // 2x2 projected (A - eval1*I) matrix on {U,V}
        T m00 = dot3(U, AU) - eval1;
	    T m01 = dot3(U, AV);
	    T m11 = dot3(V, AV) - eval1;

        // Choose the largest-length row of M to compute the eigenvector
        T absM00 = abs(m00);
        T absM01 = abs(m01);
        T absM11 = abs(m11);

        if (absM00 >= absM11)
        {
            if (absM00 > tira::constant::TIRA_EIGEN_EPSILON<T> || absM01 > tira::constant::TIRA_EIGEN_EPSILON<T>)
            {
                if (absM00 >= absM01)
                {
                    m01 /= m00;
                    T invLen = (T)1 / sqrt((T)1 + m01 * m01);
                    m00 = invLen;
                    m01 *= invLen;
                }
                else
                {
                    m00 /= m01;
                    T invLen = (T)1 / sqrt((T)1 + m00 * m00);
                    m01 = invLen;
                    m00 *= invLen;
                }
                evec1[0] = m01 * U[0] - m00 * V[0];
                evec1[1] = m01 * U[1] - m00 * V[1];
                evec1[2] = m01 * U[2] - m00 * V[2];
            }
            else
            {
                evec1[0] = U[0]; evec1[1] = U[1]; evec1[2] = U[2];
            }
        }
        else
        {
            if (absM11 > tira::constant::TIRA_EIGEN_EPSILON<T> || absM01 > tira::constant::TIRA_EIGEN_EPSILON<T>)
            {
                if (absM11 >= absM01)
                {
                    m01 /= m11;
                    T invLen = (T)1 / sqrt((T)1 + m01 * m01);
                    m11 = invLen;
                    m01 *= invLen;
                }
                else
                {
                    m11 /= m01;
                    T invLen = (T)1 / sqrt((T)1 + m11 * m11);
                    m01 = invLen;
                    m11 *= invLen;
                }
                evec1[0] = m11 * U[0] - m01 * V[0];
                evec1[1] = m11 * U[1] - m01 * V[1];
                evec1[2] = m11 * U[2] - m01 * V[2];
            }
            else
            {
                evec1[0] = U[0]; evec1[1] = U[1]; evec1[2] = U[2];
            }
        }
    }

    /**
     * This calculation is based on a description in:
     * Smith, "Eigenvalues of a symmetric 3 x 3 matrix," Communications of the ACM 4(4) 168 (April 1961)     *
     * https://doi.org/10.1145/355578.366316
     *
     * @tparam T
     * @param a
     * @param b
     * @param c
     * @param d
     * @param e
     * @param f
     * @param eval0
     * @param eval1
     * @param eval2
     */
    template<typename T>
    CUDA_CALLABLE void eval3_symmetric(T a, T b, T c, T d, T e, T f,
        T& eval0, T& eval1, T& eval2) {
        // | a   b   d |
        // | b   c   e |
        // | d   e   f |

        // determine if the matrix is diagonal
        const T p1 = b * b + d * d + e * e;
        if (p1 < tira::constant::TIRA_EIGEN_EPSILON<T>) {
            eval0 = a;  eval1 = c;  eval2 = f;
            if (eval0 > eval1) swap(eval0, eval1);
            if (eval1 > eval2) swap(eval1, eval2);
            if (eval0 > eval1) swap(eval0, eval1);
            return;
        }
        
		const T tr = a + c + f; // trace of the matrix
        const T q = tr / T(3);

		const T aq = a - q;
		const T cq = c - q;
		const T fq = f - q;

        const T p2 = aq * aq + cq * cq + fq * fq + T(2) * p1;
        T p = sqrt(p2 / T(6));

        // The matrix C = A - q*I is represented by the following, where
        // b00, b11 and b22 are computed after these comments,
        //       +-           -+         +-               -+
        //       | Ba  Bb  Bd  |      1  | a-q   b     d   |
        // B =   | Bb  Bc  Be  |    = _  | b     c-q   e   |
        //       | Bd  Be  Bf  |      p  | d     e     f-q |
        //       +-           -+         +-               -+
        const T p_inv = T(1) / p;

        const T Ba = p_inv * aq;
        const T Bb = p_inv * b;
        const T Bc = p_inv * cq;
        const T Bd = p_inv * d;
        const T Be = p_inv * e;
        const T Bf = p_inv * fq;

        // calculate the determinant of B
        const T det_B = Ba * (Bc * Bf - Be * Be) - Bb * (Bb * Bf - Bd * Be) + Bd * (Bb * Be - Bd * Bc);
        T r = det_B / T(2);

        T phi;
        if (r < T(-1)) phi = T(tira::constant::PI) / T(3);
        else if (r >= T(1)) phi = T(0);
        else phi = acos(r) / T(3);

        eval2 = q + T(2) * p * cos(phi);
        eval0 = q + T(2) * p * cos(phi + (T(2) * T(tira::constant::PI) / T(3)));
        eval1 = tr - eval0 - eval2;
    }

    /// Calculate the eigenvector of a 3x3 matrix associated with the eigenvalue lambda.
    /// The result is returned in polar coordinates (theta, phi).
    /// Avoiding the fragile cross-product method.
    template<typename T>
    CUDA_CALLABLE void evec3_symmetric(T a, T b, T c, T d, T e, T f,
        const T* evals, T* evec0, T* evec1, T* evec2) {
        // Newer version (column-major)
        // | a  b  d |
        // | b  c  e |
        // | d  e  f |

		// Test to see if this matrix is diagonal
        const T upper_diagonal_norm = b * b + d * d + e * e;
        if (upper_diagonal_norm < tira::constant::TIRA_EIGEN_EPSILON<T>) {
			// The evals are pre-sorted, we must find which diagonal element corresponds to which eigenvalue
            // and assign the correct basis vector
            
            // Create pairs of (evals, evec_idx)
            struct EvalIdx { T eval; int index; };
			EvalIdx vi[3] = { {a, 0} , {c, 1}, {f, 2} };
            T basis[3][3] = { 
                {T(1), T(0), T(0)},
                {T(0), T(1), T(0)},
				{T(0), T(0), T(1)}
            };

			// Bubble sort on the pairs based on evals (ascending)
            if (vi[0].eval > vi[1].eval) swap(vi[0], vi[1]);
            if (vi[1].eval > vi[2].eval) swap(vi[1], vi[2]);
            if (vi[0].eval > vi[1].eval) swap(vi[0], vi[1]);

            const int idx0 = vi[0].index;
			const int idx1 = vi[1].index;
			const int idx2 = vi[2].index;

            // Assign the eigenvectors based on the sorted order
			evec0[0] = basis[idx0][0]; evec0[1] = basis[idx0][1]; evec0[2] = basis[idx0][2];
			evec1[0] = basis[idx1][0]; evec1[1] = basis[idx1][1]; evec1[2] = basis[idx1][2];
			evec2[0] = basis[idx2][0]; evec2[1] = basis[idx2][1]; evec2[2] = basis[idx2][2];
            return;
        }
        
        const T q = (a + c + f) / T(3);
        const T p2 = (a - q) * (a - q) + (c - q) * (c - q) + (f - q) * (f - q) + T(2) * upper_diagonal_norm;
        const T p = sqrt(p2 / T(6));
        T halfDet;
        if (p == T(0))
            halfDet = T(0);
        else
            // formula is det(A-qI) / (2 * p^3)
            halfDet = determinant3(a - q, b, c - q, d, e, f - q) / (T(2) * (p*p*p));
        halfDet = (halfDet < T(-1)) ? T(-1) : ((halfDet > T(1)) ? T(1) : halfDet);

        if (halfDet >= T(0)) {
            evec3_symmetric_0(a, b, c, d, e, f, evals[2], evec2);
            evec3_symmetric_1(a, b, c, d, e, f, evals[1], evec2, evec1);
            cross3(evec2, evec1, evec0);        // evec0 = evec2 x evec1
        }
        else {
            evec3_symmetric_0(a, b, c, d, e, f, evals[0], evec0);
            evec3_symmetric_1(a, b, c, d, e, f, evals[1], evec0, evec1);
            cross3(evec0, evec1, evec2);        // evec2 = evec0 x evec1
            normalize3(evec2);
        }
    }

    template <typename T>
    CUDA_CALLABLE void evec3spherical_symmetric(T a, T b, T c, T d, T e, T f,
        const T* evals, T* evec0, T* evec1, T* evec2) {

        T cart_evec0[3];
        T cart_evec1[3];
        T cart_evec2[3];
        evec3_symmetric(a, b, c, d, e, f, evals, cart_evec0, cart_evec1, cart_evec2);

        // Clamp to [-1,1] to avoid NaNs from acos in case of numerical drift
        cart_evec0[2] = (cart_evec0[2] < T(-1)) ? T(-1) : ((cart_evec0[2] > T(1)) ? T(1) : cart_evec0[2]);
        cart_evec1[2] = (cart_evec1[2] < T(-1)) ? T(-1) : ((cart_evec1[2] > T(1)) ? T(1) : cart_evec1[2]);
        cart_evec2[2] = (cart_evec2[2] < T(-1)) ? T(-1) : ((cart_evec2[2] > T(1)) ? T(1) : cart_evec2[2]);

        evec0[0] = atan2(cart_evec0[1], cart_evec0[0]);
        evec0[1] = acos(cart_evec0[2]);

        evec1[0] = atan2(cart_evec1[1], cart_evec1[0]);
        evec1[1] = acos(cart_evec1[2]);

        evec2[0] = atan2(cart_evec2[1], cart_evec2[0]);
        evec2[1] = acos(cart_evec2[2]);
    }
}


/**
 * CPU namespace contains functions that are run completely on the host. All input and output pointers
 * are allocated on the host.
 */
namespace tira::cpu {
    /**
     *
     * @tparam Type is the data type used to represent the input matrices
     * @param mats is a pointer to the start of the matrix array in DRAM
     * @param n is the number of matrices in the array
     * @return a dynamically-allocated array of eigenvalue pairs (bytes = n * 2 * sizeof(Type))
     */
    template<typename Type>
    Type* evals2_symmetric(const Type* mats, const size_t n) {

        Type* evals = new Type[2*n];
        Type eval0, eval1;
        for (size_t i = 0; i < n; i++) {
            Type a = mats[i * 4 + 0];
            Type b = mats[i * 4 + 1];
            Type c = mats[i * 4 + 3];
            shared::eval2_symmetric(a, b, c, eval0, eval1);
            evals[i * 2 + 0] = eval0;
            evals[i * 2 + 1] = eval1;
        }
        return evals;
    }

    /**
	* CPU code for calculating eigenvectors of an array of 2x2 matrices given an
	* array of corresponding eigenvalues.
	* @tparam T data type for the matrix elements (ex. float)
    * @param mats pointer to an array of N 2x2 matrices (2x2xN elements of type T)
	* @param evals pointer to an array of eigenvalues (2xN elements of type T)
	* @param n number of matrices in the array
    */
    template<typename T>
    T* evecs2polar_symmetric(const T* mats, const T* evals, const size_t n) {

        T* vecs = new T[2 * n];
        T vec0, vec1;
        for (size_t i = 0; i < n; i++) {
            T a = mats[i * 4 + 0];
            T b = mats[i * 4 + 1];
            T c = mats[i * 4 + 3];
            shared::evec2polar_symmetric(a, b, c, &evals[i * 2], vec0, vec1);
            vecs[i * 2 + 0] = vec0;
            vecs[i * 2 + 1] = vec1;
        }
        return vecs;
    }

    /**
     * Calculate the eigenvalues for an array of symmetric 3x3 matrices and return the results in a new array.
     * @tparam T data type for the matrix elements
     * @param mats pointer to an array of N 3x3 matrices (3x3xN elements of type T)
     * @param n number of matrices in the array
     * @return a newly allocated array of eigenvalues in ascending order (3xN elements of type T)
     */
    template<typename T>
    T* evals3_symmetric(const T* mats, const size_t n) {

        T* evals = new T[3 * n];
        T eval0, eval1, eval2;
        for (size_t i = 0; i < n; i++) {
			T a = mats[i * 9 + 0];
			T b = mats[i * 9 + 1];
			T d = mats[i * 9 + 2];
			T c = mats[i * 9 + 4];
			T e = mats[i * 9 + 5];
			T f = mats[i * 9 + 8];

            shared::eval3_symmetric(a,b,c,d,e,f, eval0, eval1, eval2);

            evals[i * 3 + 0] = eval0;
            evals[i * 3 + 1] = eval1;
            evals[i * 3 + 2] = eval2;
        }
        return evals;
    }

    /**
     * Calculate the eigenvectors of an array of 3x3 matrices. The eigenvectors are returned in spherical coordinates (theta, phi)
     * where theta (0 < theta < 2pi) is the azimuthal angle and phi (0 < phi < pi) is the polar angle.
     * @tparam T data type for the matrix elements (ex. float)
     * @param mats pointer to an array of 3x3 matrices (3x3xN elements of type T)
     * @param lambda pointer to an array of eigenvalues (3xN elements of type T)
     * @param n number of matrices in the array
     * @return a new array of eigenvectors (3x2xN elements of type T) in spherical coordinates
     */
    template<typename T>
    T* evecs3spherical_symmetric(const T* mats, const T* lambda, const size_t n) {
        // | a  b  d |
        // | b  c  e |
        // | d  e  f |

        T* evecs = new T[3 * 2 * n];
        for (size_t i = 0; i < n; i++) {
            T a = mats[i * 9 + 0];
            T b = mats[i * 9 + 1];
            T d = mats[i * 9 + 2];
            T c = mats[i * 9 + 4];
            T e = mats[i * 9 + 5];
            T f = mats[i * 9 + 8];

			T evec0[3], evec1[3], evec2[3];
            T evals[] = { lambda[i * 3 + 0], lambda[i * 3 + 1], lambda[i * 3 + 2] };
            shared::evec3_symmetric(a, b, c, d, e, f, evals, evec0, evec1, evec2);

            // Clamp to [-1,1] to avoid NaNs from acos in case of numerical drift
            evec0[2] = (evec0[2] < T(-1)) ? T(-1) : ((evec0[2] > T(1)) ? T(1) : evec0[2]);
            evec1[2] = (evec1[2] < T(-1)) ? T(-1) : ((evec1[2] > T(1)) ? T(1) : evec1[2]);
            evec2[2] = (evec2[2] < T(-1)) ? T(-1) : ((evec2[2] > T(1)) ? T(1) : evec2[2]);

            evecs[i * 6 + 0] = std::atan2(evec0[1], evec0[0]);
            evecs[i * 6 + 1] = std::acos(evec0[2]);

            evecs[i * 6 + 2] = std::atan2(evec1[1], evec1[0]);
            evecs[i * 6 + 3] = std::acos(evec1[2]);

            evecs[i * 6 + 4] = std::atan2(evec2[1], evec2[0]);
            evecs[i * 6 + 5] = std::acos(evec2[2]);
        }
        return evecs;
    }
}

/**
 * The CUDA namespace runs everything on the GPU. Input and output pointers can be on either the host or device.
 * If pointers are located on the host, the data will be copied to the currently active CUDA device.
 * This region is only compiled when it's passed to nvcc.
 */
#ifdef __CUDACC__
namespace tira::cuda {

	// ------ Kernels ------

    /**
     *
     * @tparam Type is the data type used to represent the input matrices
     * @param mats is a pointer to the start of the matrix array in DRAM
     * @param n is the number of matrices in the array
     * @return a dynamically-allocated array of eigenvalue pairs (bytes = n * 2 * sizeof(Type))
     */
    template <typename Type>
    __global__ void kernel_eval2_symmetric(const Type* mats, const size_t n, Type* evals) {
        const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) return;

        Type eval0, eval1;
        Type a = mats[i * 4 + 0];
        Type b = mats[i * 4 + 1];
        Type c = mats[i * 4 + 3];
        shared::eval2_symmetric(a, b, c, eval0, eval1);
        evals[i * 2 + 0] = eval0;
        evals[i * 2 + 1] = eval1;
    }

    template <typename Type>
    __global__ void kernel_evec2polar_symmetric(const Type* mats, const Type* evals, const size_t n, Type* evecs) {
        const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) return;
        Type a = mats[i * 4 + 0];
        Type b = mats[i * 4 + 1];
        Type c = mats[i * 4 + 3];
        tira::shared::evec2polar_symmetric(a, b, c, &evals[i * 2], evecs[i * 2 + 0], evecs[i * 2 + 1]);
	}

	template<typename Type>
    __global__ void kernel_eval3_symmetric(const Type* mats, const size_t n, Type* evals) {
        const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) return;

        const Type a = mats[i * 9 + 0];
        const Type b = mats[i * 9 + 1];
        const Type d = mats[i * 9 + 2];
        const Type c = mats[i * 9 + 4];
        const Type e = mats[i * 9 + 5];
        const Type f = mats[i * 9 + 8];

        Type l0, l1, l2;
        tira::shared::eval3_symmetric(a, b, c, d, e, f, l0, l1, l2);

        evals[i * 3 + 0] = l0;
        evals[i * 3 + 1] = l1;
        evals[i * 3 + 2] = l2;
    }

	template<typename Type>
    __global__ void kernel_evec3spherical_symmetric(const Type* mats, const Type* evals, const size_t n, Type* evecs) {
        const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) return;
        const Type a = mats[i * 9 + 0];
        const Type b = mats[i * 9 + 1];
        const Type d = mats[i * 9 + 2];
        const Type c = mats[i * 9 + 4];
        const Type e = mats[i * 9 + 5];
        const Type f = mats[i * 9 + 8];
        tira::shared::evec3spherical_symmetric(a, b, c, d, e, f,
            &evals[i * 3],
            &evecs[i * 6],
            &evecs[i * 6 + 2],
            &evecs[i * 6 + 4]);
	}

	// ------ Host-side launchers ------

    /**
     * @brief Compute eigenvalues of N symmetric 2x2 matrices on the GPU.
     *
     * @tparam Type is the data type used to represent the input matrices
     * @param mats_device Pointer to N 2x2 matrices in device memory.
     * @param n           Number of matrices.
     * @return Type*      Pointer to eigenvalues in device memory (2*N elements).
     *                    Caller is responsible for cudaFree().
     */
    template<typename Type>
    Type* evals2_symmetric(const Type* mats, const size_t n, int device) {
        if (device < 0) return cpu::evals2_symmetric(mats, n);

        const Type* gpu_mats;
        Type* temp_gpu_mats;
        size_t mats_bytes = sizeof(Type) * 4 * n;
        size_t evals_bytes = sizeof(Type) * 2 * n;

        // determine if the source image is provided on the CPU or GPU
        cudaPointerAttributes attribs;										// create a pointer attribute structure
        HANDLE_ERROR(cudaPointerGetAttributes(&attribs, mats));			// get the attributes for the source pointer

        // if the provided pointer is on the device set the gpu_source pointer to source
        if (attribs.type == cudaMemoryTypeDevice)   gpu_mats = mats;

        else {																// otherwise copy the source image to the GPU
            HANDLE_ERROR(cudaMalloc(&temp_gpu_mats, mats_bytes));								// allocate space on the GPU for the source image
            HANDLE_ERROR(cudaMemcpy(temp_gpu_mats, mats, mats_bytes, cudaMemcpyHostToDevice));// copy the source image to the GPU
            gpu_mats = static_cast<const Type*>(temp_gpu_mats);
        }

        // get the active device properties to calculate the optimal the block size
        HANDLE_ERROR(cudaGetDevice(&device));
        cudaDeviceProp props;
        HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
        const unsigned int blockSize = std::min(256u, static_cast<unsigned int>(props.maxThreadsPerBlock));
        dim3 blockDim(blockSize);
        dim3 gridDim(static_cast<unsigned int>((n + blockDim.x - 1) / blockDim.x));

        Type* gpu_evals;
        HANDLE_ERROR(cudaMalloc(&gpu_evals, evals_bytes));
        kernel_eval2_symmetric << <gridDim, blockDim >> > (gpu_mats, n, gpu_evals);
        HANDLE_ERROR(cudaGetLastError());

        Type* evals;
        if (attribs.type == cudaMemoryTypeDevice)   evals = gpu_evals;
        else {
            evals = new Type[2 * n];
            HANDLE_ERROR(cudaMemcpy(evals, gpu_evals, evals_bytes, cudaMemcpyDeviceToHost));
            HANDLE_ERROR(cudaFree(gpu_evals));
            HANDLE_ERROR(cudaFree(temp_gpu_mats));
        }
        return evals;
    }

	// TODO: assert both mats and evals are on the same side (host or device)
    // or independently inspect lambda. Just like the 3D version.
    template<typename Type>
    Type* evecs2polar_symmetric(const Type* mats, const Type* evals, size_t n, int device) {
		if (device < 0)  return cpu::evecs2polar_symmetric(mats, evals, n);

		HANDLE_ERROR(cudaSetDevice(device));

        size_t mats_bytes   = sizeof(Type) * 4 * n;
        size_t ev_bytes     = sizeof(Type) * 2 * n;

        const Type* gpu_mats    = nullptr;
        const Type* gpu_evals   = nullptr;
        Type* temp_gpu_mats     = nullptr;
        Type* temp_gpu_evals    = nullptr;

        // determine if the source image is provided on the CPU or GPU
        cudaPointerAttributes attribs{};										// create a pointer attribute structure
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
            gpu_mats = static_cast<const Type*>(temp_gpu_mats);
            gpu_evals = static_cast<const Type*>(temp_gpu_evals);
        }

        // get the active device properties to calculate the optimal the block size
        HANDLE_ERROR(cudaGetDevice(&device));
        cudaDeviceProp props;
        HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
        const unsigned int blockSize = std::min(256u, static_cast<unsigned int>(props.maxThreadsPerBlock));
        dim3 blockDim(blockSize);
        dim3 gridDim(static_cast<unsigned int>((n + blockDim.x - 1) / blockDim.x));

        Type* gpu_evecs;
        HANDLE_ERROR(cudaMalloc(&gpu_evecs, ev_bytes));
        kernel_evec2polar_symmetric << <gridDim, blockDim >> > (gpu_mats, gpu_evals, n, gpu_evecs);
        HANDLE_ERROR(cudaGetLastError());

        Type* evecs;
        if (attribs.type == cudaMemoryTypeDevice)
            evecs = gpu_evecs;
        else {
            evecs = new Type[2 * n];
            HANDLE_ERROR(cudaMemcpy(evecs, gpu_evecs, ev_bytes, cudaMemcpyDeviceToHost));
            HANDLE_ERROR(cudaFree(gpu_evecs));
            HANDLE_ERROR(cudaFree(temp_gpu_evals));
            HANDLE_ERROR(cudaFree(temp_gpu_mats));
        }
        return evecs;
    }

    template<typename Type>
    Type* evals3_symmetric(const Type* mats, const size_t n, int device) {
        if (device < 0)  return cpu::evals3_symmetric(mats, n);

		HANDLE_ERROR(cudaSetDevice(device));

        // Set up sizes for GPU storage
        const size_t mats_bytes  = sizeof(Type) * 9 * n;                              // required bytes for storing the tensor
        const size_t evals_bytes = sizeof(Type) * 3 * n;                             // required bytes for storing eigenvalues

        // Set up pointers
        const Type* gpu_mats = nullptr;
        Type* temp_gpu_mats  = nullptr;

        // Determine if the source volume is provided on the CPU or GPU
        cudaPointerAttributes attribs{};									// create a pointer attribute structure
        HANDLE_ERROR(cudaPointerGetAttributes(&attribs, mats));			    // get the attributes for the source pointer

		const bool mats_on_device = (attribs.type == cudaMemoryTypeDevice);
        if (mats_on_device)   gpu_mats = mats;
        else {																// otherwise copy the source image to the GPU
            HANDLE_ERROR(cudaMalloc(&temp_gpu_mats, mats_bytes));								// allocate space on the GPU for the source image
            HANDLE_ERROR(cudaMemcpy(temp_gpu_mats, mats, mats_bytes, cudaMemcpyHostToDevice));   // copy the source image to the GPU
            gpu_mats = temp_gpu_mats;
        }

        // get the active device properties to calculate the optimal the block size
        cudaDeviceProp props;
        HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
		const unsigned int blockSize = std::min(256u, static_cast<unsigned int>(props.maxThreadsPerBlock));
        dim3 blockDim(blockSize);
        dim3 gridDim(static_cast<unsigned int>((n + blockDim.x - 1) / blockDim.x));

        Type* gpu_evals = nullptr;
        HANDLE_ERROR(cudaMalloc(&gpu_evals, evals_bytes));
        kernel_eval3_symmetric << <gridDim, blockDim >> > (gpu_mats, n, gpu_evals);
        HANDLE_ERROR(cudaGetLastError());

        if (mats_on_device)
            return gpu_evals;
        else {
            Type* evals = new Type[3 * n];
            HANDLE_ERROR(cudaMemcpy(evals, gpu_evals, evals_bytes, cudaMemcpyDeviceToHost));
            HANDLE_ERROR(cudaFree(gpu_evals));
            if (temp_gpu_mats) HANDLE_ERROR(cudaFree(temp_gpu_mats));
            return evals;
        }
    }

    template<typename Type>
    Type* evecs3spherical_symmetric(const Type* mats, const Type* lambda, const size_t n, int device) {
        if (device < 0)  return cpu::evecs3spherical_symmetric(mats, lambda, n);

		HANDLE_ERROR(cudaSetDevice(device));

        const Type* gpu_mats    = nullptr;
        const Type* gpu_lambda  = nullptr;
        Type* temp_gpu_mats     = nullptr;
        Type* temp_gpu_lambda   = nullptr;
        const size_t mats_bytes  = sizeof(Type) * 9 * n;                              // required bytes for storing the tensor
        const size_t evals_bytes = sizeof(Type) * 3 * n;                             // required bytes for storing eigenvalues
        const size_t evecs_bytes = sizeof(Type) * 6 * n;                             // required bytes for storing 2 3D eigenvectors (in polar coordinates)

        // Determine if the volume is provided on the host or device
        cudaPointerAttributes matsAttr{};										// create a pointer attribute structure
        HANDLE_ERROR(cudaPointerGetAttributes(&matsAttr, mats));			    // get the attributes for the source pointer

		const bool mats_on_device = (matsAttr.type == cudaMemoryTypeDevice);
        if (mats_on_device)  gpu_mats = mats;
        else {																// otherwise copy the source image to the GPU
            HANDLE_ERROR(cudaMalloc(&temp_gpu_mats, mats_bytes));								    // allocate space on the GPU for the source image
            HANDLE_ERROR(cudaMemcpy(temp_gpu_mats, mats, mats_bytes, cudaMemcpyHostToDevice));      // copy the source image to the GPU
            gpu_mats = temp_gpu_mats;
        }

        // Determine if the eigenvalues are provided on the host or device
        cudaPointerAttributes lamAttr{};										// create a pointer attribute structure
        HANDLE_ERROR(cudaPointerGetAttributes(&lamAttr, lambda));			    // get the attributes for the source pointer
        
		const bool lambda_on_device = (lamAttr.type == cudaMemoryTypeDevice);
        if (lambda_on_device)   gpu_lambda = lambda;
        else {																    // otherwise copy the source image to the GPU
            HANDLE_ERROR(cudaMalloc(&temp_gpu_lambda, evals_bytes));
            HANDLE_ERROR(cudaMemcpy(temp_gpu_lambda, lambda, evals_bytes, cudaMemcpyHostToDevice));
            gpu_lambda = temp_gpu_lambda;
        }

        // Set the optimal grid and block size
        cudaDeviceProp props{};
        HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
        const unsigned int blockSize = std::min(256u, static_cast<unsigned int>(props.maxThreadsPerBlock));
        dim3 blockDim(blockSize);
        dim3 gridDim(static_cast<unsigned int>((n + blockDim.x - 1) / blockDim.x));

		// Launche kernel to calculate eigenvectors
        Type* gpu_evecs = nullptr;
        HANDLE_ERROR(cudaMalloc(&gpu_evecs, evecs_bytes));
        kernel_evec3spherical_symmetric << <gridDim, blockDim >> > (gpu_mats, gpu_lambda, n, gpu_evecs);
        HANDLE_ERROR(cudaGetLastError());

        // Decide where to store the results
		const bool keep_on_device = mats_on_device || lambda_on_device;

        // Return / cleanup
        Type* evecs = nullptr;
        if (keep_on_device)
			evecs = gpu_evecs;              // caller cudaFree()
        else {
            evecs = new Type[6 * n];
            HANDLE_ERROR(cudaMemcpy(evecs, gpu_evecs, evecs_bytes, cudaMemcpyDeviceToHost));
            HANDLE_ERROR(cudaFree(gpu_evecs));
        }
		// Free temporary device copies if we allocated them
        if (temp_gpu_mats)    HANDLE_ERROR(cudaFree(temp_gpu_lambda));
        if (temp_gpu_lambda)  HANDLE_ERROR(cudaFree(temp_gpu_mats));

        return evecs;
    }
}
#endif

/**
 * Any functions in this general namespace will look at the associated pointers or device numbers and
 * determine which functions in tira::cpu and tira::cuda are executed. This might be unnecessary.
 */
namespace tira {

}
