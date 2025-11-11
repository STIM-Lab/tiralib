#pragma once

#include <tira/cuda/callable.h>

#define PI 3.14159265358979323846
#define TIRA_EIGEN_EPSILON 1e-12

namespace tira {

    // function used to quickly swap two values
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
        if (n > T(0)) { v[0] /= n; v[1] /= n;   v[2] /= n;  }
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

        eval0 = std::abs(l0) < std::abs(l1) ? l0 : l1;  // sort the eigenvalues based on their magnitude
        eval1 = std::abs(l0) > std::abs(l1) ? l0 : l1;
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

        T abs_a_l0 = std::abs(a_l0);
        T abs_a_l1 = std::abs(a_l1);
        if (abs_a_l0 >= abs_a_l1) {
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
        T inv;
        if (std::fabs(evec[0]) > std::fabs(evec[1])) {
            // The component of maximum absolute value is either evec[0] or evec[2]
		    inv = T(1) / sqrt(evec[0] * evec[0] + evec[2] * evec[2]);
            U[0] = -evec[2] * inv;  U[1] = T(0);    U[2] = evec[0] * inv; // U is orthogonal to evec
        }
        else {
		    // The component of maximum absolute value is either evec[1] or evec[2]
		    inv = T(1) / sqrt(evec[1] * evec[1] + evec[2] * evec[2]);
            U[0] = T(0);    U[1] = evec[2] * inv;   U[2] = -evec[1] * inv;
        }
		cross3(evec, U, V); // V = evec x U
    }

    template <typename T>
    CUDA_CALLABLE void evec3_symmetric_0(const T a, const T b, const T c, const T d, const T e, const T f,
        const T eval0, T* evec0) {
        // | a  b  d |
	    // | b  c  e |
	    // | d  e  f |

        const T upper_diagonal_norm = b * b + d * d + e * e;
        if (upper_diagonal_norm < T(TIRA_EIGEN_EPSILON)) {
            evec0[0] = T(1);    evec0[1] = T(0);    evec0[2] = T(0);
            return;
        }

        // Solve (A - LI) v = 0 via cross products of two rows with largest area
	    T r0[3] = { a - eval0,  b,          d           };
	    T r1[3] = { b,          c - eval0,  e           };
	    T r2[3] = { d,          e,          f - eval0   };

        // calculate the cross-product of each two rows
        T r0xr1[3], r0xr2[3], r1xr2[3];
        cross3(r0, r1, r0xr1);
	    cross3(r0, r2, r0xr2);
	    cross3(r1, r2, r1xr2);

        // dot products - to find out which cross-product has the largest length
		const T d01 = dot3(r0xr1, r0xr1);
		const T d02 = dot3(r0xr2, r0xr2);
		const T d12 = dot3(r1xr2, r1xr2);

        T dmax = d01;
        int imax = 0;
        if (d02 > dmax) { dmax = d02; imax = 1; }
        if (d12 > dmax) { imax = 2; }
        if (imax == 0) {
            // r0xr1 is the largest cross product
			T norm = sqrt(d01);
            evec0[0] = r0xr1[0] / norm;
            evec0[1] = r0xr1[1] / norm;
            evec0[2] = r0xr1[2] / norm;
        }
        else if (imax == 1) {
            T norm = sqrt(d02);
            // r0xr2 is the largest cross product
            evec0[0] = r0xr2[0] / norm;
            evec0[1] = r0xr2[1] / norm;
            evec0[2] = r0xr2[2] / norm;
        }
        else {
            // r1xr2 is the largest cross product
            T norm = sqrt(d12);
            evec0[0] = r1xr2[0] / norm;
            evec0[1] = r1xr2[1] / norm;
            evec0[2] = r1xr2[2] / norm;
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
	    
        // 2x2 projected (A - LI) matrix on {U,V}
        T m00 = dot3(U, AU) - eval1;
	    T m01 = dot3(U, AV);
	    T m11 = dot3(V, AV) - eval1;

        // Choose the largest-length row of M to compute the eigenvector
        T absM00 = std::fabs(m00);
        T absM01 = std::fabs(m01);
        T absM11 = std::fabs(m11);
        T maxAbsComp;

        if (absM00 >= absM11) {
            maxAbsComp = (absM00 >= absM01) ? absM00 : absM01;
            if (maxAbsComp > T(0)) {
                if (absM00 >= absM01) {
                    m01 /= m00;
                    m00 = T(1) / std::sqrt(T(1) + m01 * m01);
                    m01 *= m00;
                }
                else {
                    m00 /= m01;
                    m01 = T(1) / std::sqrt(T(1) + m00 * m00);
                    m00 *= m01;
                }
				//evec1 = subtract3(multiply3(m01, U), multiply3(m00, V));
                evec1[0] = m01 * U[0] - m00 * V[0];
                evec1[1] = m01 * U[1] - m00 * V[1];
                evec1[2] = m01 * U[2] - m00 * V[2];
            }
            else {
                evec1[0] = U[0];    evec1[1] = U[1];    evec1[2] = U[2];
            }
        }
        else {
            maxAbsComp = (absM11 >= absM01) ? absM11 : absM01;
            if (maxAbsComp > T(0)) {
                if (absM11 >= absM01) {
                    m01 /= m11;
                    m11 = T(1) / std::sqrt(T(1) + m01 * m01);
                    m01 *= m11;
                }
                else {
                    m11 /= m01;
                    m01 = T(1) / std::sqrt(T(1) + m11 * m11);
                    m11 *= m01;
                }
                // Subtract(Multiply(m11, U), Multiply(m01, V))
				// evec1 = subtract3(multiply3(m11, U), multiply3(m01, V));
                evec1[0] = m11 * U[0] - m01 * V[0];
                evec1[1] = m11 * U[1] - m01 * V[1];
                evec1[2] = m11 * U[2] - m01 * V[2];

            }
            else {
                evec1[0] = U[0];    evec1[1] = U[1];    evec1[2] = U[2];
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
            if (eval0 > eval1) std::swap(eval0, eval1);
            if (eval1 > eval2) std::swap(eval1, eval2);
            if (eval0 > eval1) std::swap(eval0, eval1);
            return;
        }
        
		const T tr = a + c + f; // trace of the matrix
        const T q = tr / T(3);
        const T p2 = pow(a-q, 2) + pow(c-q, 2) + pow(f-q, 2) + T(2) * p1;
        T p = sqrt(p2 / T(6));
        if (std::isnan(p)) {
            eval0 = -111.0; eval1 = -111.0; eval2 = -111.0;
            return; // Exit and report "error -111"
        }

        // The matrix C = A - q*I is represented by the following, where
        // b00, b11 and b22 are computed after these comments,
        //       +-           -+         +-               -+
        //       | Ba  Bb  Bd  |      1  | a-q   b     d   |
        // B =   | Bb  Bc  Be  |    = _  | b     c-q   e   |
        //       | Bd  Be  Bf  |      p  | d     e     f-q |
        //       +-           -+         +-               -+
        const T p_inv = T(1) / p;
        if (std::isnan(p_inv)) {
            eval0 = -222.0; eval1 = -222.0; eval2 = -222.0;
            return; // Exit and report "error -222"
        }
        const T Ba = p_inv * (a - q);
        const T Bb = p_inv * b;
        const T Bc = p_inv * (c - q);
        const T Bd = p_inv * d;
        const T Be = p_inv * e;
        const T Bf = p_inv * (f - q);

        // calculate the determinant of B
        const T det_B = Ba * (Bc * Bf - Be * Be) - Bb * (Bb * Bf - Bd * Be) + Bd * (Bb * Be - Bd * Bc);
        T r = det_B / T(2);

        T phi;
        if (r < T(-1)) phi = T(PI) / T(3);
        else if (r >= T(1)) phi = T(0);
        else phi = acos(r) / T(3);
        if (std::isnan(phi)) {
            eval0 = -333.0; eval1 = -333.0; eval2 = -333.0;
            return; // Exit and report "error -333"
        }

        eval2 = q + T(2) * p * cos(phi);
        eval0 = q + T(2) * p * cos(phi + (T(2) * T(PI) / T(3)));
        eval1 = tr - eval0 - eval2;
        if (std::isnan(eval0)) {
            eval0 = -444.0; eval1 = -444.0; eval2 = -444.0;
            return; // Exit and report "error -444"
		}
        if (std::isnan(eval2)) {
            eval0 = -555.0; eval1 = -555.0; eval2 = -555.0;
			return; // Exit and report "error -555"
        }
    }

    /**
     * Calculate the eigenvalues of a 3x3 symmetric matrix. This code is adapted from a paper by David Eberly
     * available here:
     * https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf
     *
     * Links to the adopted code are provided and Eberly's comments are appended in the adopted code for the function.
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
    /*template<typename T>
    CUDA_CALLABLE void eval3_symmetric(T a, T b, T c, T d, T e, T f,
        T& eval0, T& eval1, T& eval2) {

        eval3_symmetric_wikipedia(a, b, c, d, e, f, eval0, eval1, eval2);
        return;
	    // | a   b   d |
	    // | b   c   e |
	    // | d   e   f |

        // To guard against floating-point overflow, we precondition the matrix by normalizing by the largest value
        T max0 = (fabs(a) > fabs(b)) ? fabs(a) : fabs(b);              // calculate the largest absolute value in the matrix
        T max1 = (fabs(d) > fabs(c)) ? fabs(d) : fabs(c);
        T max2 = (fabs(e) > fabs(f)) ? fabs(e) : fabs(f);
        T maxElement = (max0 > max1) ? ((max0 > max2) ? max0 : max2) : ((max1 > max2) ? max1 : max2);
        if (maxElement == 0.0) {
            eval0 = 0.0; eval1 = 0.0; eval2 = 0.0; // if the matrix is zero, return zero eigenvalues
            return;
        }

        T invMaxElement = 1.0 / maxElement;                            // normalize the matrix
        a *= invMaxElement; b *= invMaxElement; d *= invMaxElement;
        c *= invMaxElement; e *= invMaxElement; f *= invMaxElement;

        // Case: matrix is diagonal
        T p1 = b * b + d * d + e * e;
        if (p1 == 0.0) {
            eval0 = a * maxElement;
            eval1 = c * maxElement;
            eval2 = f * maxElement;
            if (eval0 > eval2) swap(eval0, eval2);
            if (eval0 > eval1) swap(eval0, eval1);
            if (eval1 > eval2) swap(eval1, eval2);
            return;
        }

        // The PDF defines the matrix B = (A - q*I)/p, where:
        // q = tr(A)/3
        // p = sqrt(tr((A - q*I)^2)/6)

        T q = (a + c + f) / 3.0;            // calculate q = tr(A) / 3

        // The matrix C = A - q*I is represented by the following, where
        // b00, b11 and b22 are computed after these comments,
        //   +-           -+        +-               -+
        //   | b00 a01 a02 |        | a-q   b     d   |
        //   | a01 b11 a12 |    =   | b     c-q   e   |
        //   | a02 a12 b22 |        | d     b     f-q |
        //   +-           -+        +-               -+
        //T tr_C = (a - q) + (c - q) + (f - q);
        //T tr_C2 = pow(a - q, 2) + pow(c - q, 2) + pow(f - q, 2);            // calculate the tr(C^2)
        T p2 = (a - q) * (a - q) + (c - q) * (c - q) + (f - q) * (f - q) + 2.0 * p1;
        T p = sqrt(p2 / 6.0);
        T pinv = 1.0 / p;
        T B[3][3];
        B[0][0] = pinv * (a - q);
        B[1][1] = pinv * (c - q);
        B[2][2] = pinv * (f - q);
        B[0][1] = pinv * b;
        B[0][2] = pinv * d;
        B[1][2] = pinv * e;
        B[1][0] = B[0][1];
        B[2][0] = B[0][2];
        B[2][1] = B[1][2];
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
        eval2 = (q + 2.0 * p * cos(phi))                    * maxElement;
        eval0 = (q + 2.0 * p * cos(phi + (2.0 * PI / 3.0))) * maxElement;
        eval1 = (3.0 * q - eval2 - eval0)                   * maxElement;                // since trace(A) = eig1 + eig2 + eig3
    }
    */

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

		// test to see if this matrix is diagonal, return basis vectors if it is
        const T upper_diagonal_norm = b * b + d * d + e * e;
        if (upper_diagonal_norm < T(TIRA_EIGEN_EPSILON)) {
            evec0[0] = T(1);    evec0[1] = T(0);    evec0[2] = T(0);
            evec1[0] = T(0);    evec1[1] = T(1);    evec1[2] = T(0);
            evec2[0] = T(0);    evec2[1] = T(0);    evec2[2] = T(1);
            return;
        }
        
        const T q = (a + c + f) / T(3);
        const T p2 = pow(a - q, 2) + pow(c - q, 2) + pow(f - q, 2) + T(2) * upper_diagonal_norm;
        const T p = sqrt(p2 / T(6));
        T halfDet;
        if (p == T(0))
            halfDet = T(0);
        else
            // formula is det(A-qI) / (2 * p^3)
            halfDet = determinant3(a - q, b, c - q, d, e, f - q) / (T(2) * pow(p, 3));
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

		// Clamp to [-1,1] to avoid NaNs from acos in case of numerical drift
		evec0[2] = (evec0[2] < T(-1)) ? T(-1) : ((evec0[2] > T(1)) ? T(1) : evec0[2]);
		evec1[2] = (evec1[2] < T(-1)) ? T(-1) : ((evec1[2] > T(1)) ? T(1) : evec1[2]);
		evec2[2] = (evec2[2] < T(-1)) ? T(-1) : ((evec2[2] > T(1)) ? T(1) : evec2[2]);
    }

    template <typename T>
    CUDA_CALLABLE void evec3spherical_symmetric(T a, T b, T c, T d, T e, T f,
        const T* evals, T* evec0, T* evec1, T* evec2) {

        T cart_evec0[3];
        T cart_evec1[3];
        T cart_evec2[3];
        evec3_symmetric(a, b, c, d, e, f, evals, cart_evec0, cart_evec1, cart_evec2);

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
            /*if (eval0 == -111 || eval0 == -222 || eval0 == -333 || eval0 == -444 
                || eval0 == -555 and func_vote) {
                std::cout << "----------------------------------------------" << std::endl;
                std::cout << "Wrong evals at " << i << std::endl;
                std::cout << "eval0: " << eval0 << " eval1: " << eval1 << " eval2: " << eval2 << std::endl;
				std::cout << "Matrix: " << std::endl;
				std::cout << a << " " << b << " " << d << std::endl;
				std::cout << b << " " << c << " " << e << std::endl;
                std::cout << d << " " << e << " " << f << std::endl;
            }*/
			// The new eigenvalues are scaled. Revert the scaling
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

