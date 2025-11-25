#pragma once

#include <tira/cuda/callable.h>

#define PI 3.14159265358979323846
#define TIRA_EIGEN_EPSILON 1e-12

namespace tira {
    template <typename T>
    CUDA_CALLABLE void swap(T& a, T& b) {
        T temp = a;
        a = b;
        b = temp;
    }

    template <typename T>
    CUDA_CALLABLE T sgn(T& a) {
        if (a > (T)0) return (T)(1);
        if (a < (T)0) return (T)(-1);
        else return (T)0;
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
        if (n > T(TIRA_EIGEN_EPSILON)) {
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
        T q = -0.5 * (b + sgn(b) * sqrt(disc));
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
            theta1 = std::atan2(b, a_l0);
            theta0 = theta1 + (PI / 2.0);
            if (theta0 > PI) theta0 -= 2 * PI;
        }
        else {
            theta0 = std::atan2(b, a_l1);
            theta1 = theta0 + (PI / 2.0);
            if (theta1 > PI) theta1 -= 2 * PI;
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
        if (absY - absX > TIRA_EIGEN_EPSILON) {
            if (absZ - absX > TIRA_EIGEN_EPSILON) {
                temp[0] = 1.0f; temp[1] = 0.0f; temp[2] = 0.0f;         // X is smallest
            }
            else {
            temp[0] = 0.0f; temp[1] = 0.0f; temp[2] = 1.0f;         // Z is smallest
            }
        }
        else {
            if (absZ - absY > TIRA_EIGEN_EPSILON) {
                temp[0] = 0.0f; temp[1] = 1.0f; temp[2] = 0.0f;          // Y is smallest
            }
            else {
                temp[0] = 0.0f; temp[1] = 0.0f; temp[2] = 1.0f;         // Z is smallest
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
            if (absM00 > T(TIRA_EIGEN_EPSILON) || absM01 > T(TIRA_EIGEN_EPSILON))
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
            if (absM11 > T(TIRA_EIGEN_EPSILON) || absM01 > T(TIRA_EIGEN_EPSILON))
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
        if (p1 < T(TIRA_EIGEN_EPSILON)) {
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
        if (r < T(-1)) phi = T(PI) / T(3);
        else if (r >= T(1)) phi = T(0);
        else phi = acos(r) / T(3);

        eval2 = q + T(2) * p * cos(phi);
        eval0 = q + T(2) * p * cos(phi + (T(2) * T(PI) / T(3)));
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
        if (upper_diagonal_norm < T(TIRA_EIGEN_EPSILON)) {
			// The evals are pre-sorted, we must find which diagonal element corresponds to which eigenvalue
            // and assign the correct basis vector
            
            // Create pairs of (evals, evec_idx)
			T val_idx[3][2] = { {a, 0} , {c, 1}, {f, 2} };
			T basis[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };

			// Bubble sort on the pairs based on evals
            if (val_idx[0][0] > val_idx[1][0]) {
                swap(val_idx[0][0], val_idx[1][0]);
                swap(val_idx[0][1], val_idx[1][1]);
			}
            if (val_idx[1][0] > val_idx[2][0]) {
                swap(val_idx[1][0], val_idx[2][0]);
                swap(val_idx[1][1], val_idx[2][1]);
            }
            if (val_idx[0][0] > val_idx[1][0]) {
                swap(val_idx[0][0], val_idx[1][0]);
                swap(val_idx[0][1], val_idx[1][1]);
			}

			// Assign the eigenvectors based on the sorted order
			int idx0 = (int)val_idx[0][1];
			int idx1 = (int)val_idx[1][1];
			int idx2 = (int)val_idx[2][1];

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

        evec0[0] = std::atan2(cart_evec0[1], cart_evec0[0]);
        evec0[1] = std::acos(cart_evec0[2]);

        evec1[0] = std::atan2(cart_evec1[1], cart_evec1[0]);
        evec1[1] = std::acos(cart_evec1[2]);

        evec2[0] = std::atan2(cart_evec2[1], cart_evec2[0]);
        evec2[1] = std::acos(cart_evec2[2]);
    }
}



namespace tira::cpu {

    /// <summary>
    /// CPU code for calculating eigenvalues of a 2D matrix array
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <param name="mats"></param>
    /// <param name="n"></param>
    /// <returns></returns>
    template<typename T>
    T* evals2_symmetric(const T* mats, const size_t n) {

        T* evals = new T[2*n];
        T eval0, eval1;
        for (size_t i = 0; i < n; i++) {
            T a = mats[i * 4 + 0];
            T b = mats[i * 4 + 1];
            T c = mats[i * 4 + 3];
            eval2_symmetric(a, b, c, eval0, eval1);
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
            evec2polar_symmetric(a, b, c, &evals[i * 2], vec0, vec1);
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

            eval3_symmetric(a,b,c,d,e,f, eval0, eval1, eval2);

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
        for (unsigned int i = 0; i < n; i++) {
            T a = mats[i * 9 + 0];
            T b = mats[i * 9 + 1];
            T d = mats[i * 9 + 2];
            T c = mats[i * 9 + 4];
            T e = mats[i * 9 + 5];
            T f = mats[i * 9 + 8];

			// Now we can safely calculate the eigenvectors
			T evec0[3], evec1[3], evec2[3];
            T evals[] = { lambda[i * 3 + 0], lambda[i * 3 + 1], lambda[i * 3 + 2] };
            evec3_symmetric(a, b, c, d, e, f, evals, evec0, evec1, evec2);

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

