#pragma once

#include <tira/cuda/callable.h>

#define PI 3.14159265358979323846


// function used to quickly swap two values
template <typename T>
CUDA_CALLABLE void swap(T& a, T& b) {
    T temp = a;
    a = b;
    b = temp;
}

template <typename T>
CUDA_CALLABLE T trace2(const T* matrix) {
    return matrix[0] + matrix[3];
}

template <typename T>
CUDA_CALLABLE T determinant2(const T* matrix) {
    return matrix[0] * matrix[3] - matrix[1] * matrix[2];
}

template <typename T>
CUDA_CALLABLE T determinanat3(const T a, const T b, const T c, const T d, 
    const T e, const T f) {
	// | a  b  d |
	// | b  c  e |
	// | d  e  f |
	return a * (c * f - e * e) - b * (b * f - e * d) + d * (b * e - c * d);
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
template <typename T>
CUDA_CALLABLE int sgn(T x) {
    if (x < 0) return -1;
    else return 1;
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


/// Robustly compute a right-handed orthonormal set { U, V, W }.
/// The vector W is guaranteed to be unit-length, in which case
/// there is no need to worry about a division by zero when computing invLength.
template<typename T>
CUDA_CALLABLE void ComputeOrthogonalComplement(const T* evec, T* U, T* V) {
    T invLength;
    if (std::fabs(evec[0]) > std::fabs(evec[1])) {
        // The component of maximum absolute value is either evec[0] or evec[2]
		invLength = (T)1 / sqrt(evec[0] * evec[0] + evec[2] * evec[2]);
        U[0] = -evec[2] * invLength;
        U[1] = (T)0;
        U[2] = evec[0] * invLength; // U is orthogonal to evec
    }
    else {
		// The component of maximum absolute value is either evec[1] or evec[2]
		invLength = (T)1 / sqrt(evec[1] * evec[1] + evec[2] * evec[2]);
        U[0] = (T)0;
        U[1] = evec[2] * invLength;
        U[2] = - evec[1] * invLength;
    }
    V[0] = evec[1] * U[2] - evec[2] * U[1];
    V[1] = evec[2] * U[0] - evec[0] * U[2];
	V[2] = evec[0] * U[1] - evec[1] * U[0];
}

template <typename T>
CUDA_CALLABLE void ComputeEigenvector0(const T a, const T b, const T c, const T d, const T e, const T f,
    const T eval0, T* evec0) {
    // | a  b  d |
	// | b  c  e |
	// | d  e  f |

	T r0[] = { a - eval0, b, d };
	T r1[] = { b, c - eval0, e };
	T r2[] = { d, e, f - eval0 };

    // calculate the cross-product of each two rows
    T r0xr1[3], r0xr2[3], r1xr2[3];
    cross3(r0, r1, r0xr1);
	cross3(r0, r2, r0xr2);
	cross3(r1, r2, r1xr2);

    // dot products - to find out which cross-product has the largest length
    T d01 = r0xr1[0] * r0xr1[0] + r0xr1[1] * r0xr1[1] + r0xr1[2] * r0xr1[2];
    T d02 = r0xr2[0] * r0xr2[0] + r0xr2[1] * r0xr2[1] + r0xr2[2] * r0xr2[2];
    T d12 = r1xr2[0] * r1xr2[0] + r1xr2[1] * r1xr2[1] + r1xr2[2] * r1xr2[2];

    T* r = r0xr1;
    T norm = sqrt(d01);

    if (d02 > d01 && d02 > d12) {
        r = r0xr2;
        norm = sqrt(d02);
    }
    else if (d12 > d01 && d12 > d02) {
        r = r1xr2;
        norm = sqrt(d12);
    }

    evec0[0] = r[0] / norm;
    evec0[1] = r[1] / norm;
    evec0[2] = r[2] / norm;
}

template <typename T>
CUDA_CALLABLE void ComputeEigenvector1(const T a, const T b, const T c, const T d, const T e, const T f,
    const T eval1, const T* evec0, T* evec1) {
	T* U = new T[3];
	T* V = new T[3];
	ComputeOrthogonalComplement(evec0, U, V);

	T AU[] = {
		a * U[0] + b * U[1] + d * U[2],
		b * U[0] + c * U[1] + e * U[2],
		d * U[0] + e * U[1] + f * U[2]
	};
	T AV[] = {
		a * V[0] + b * V[1] + d * V[2],
		b * V[0] + c * V[1] + e * V[2],
		d * V[0] + e * V[1] + f * V[2]
	};
	T m00 = U[0] * AU[0] + U[1] * AU[1] + U[2] * AU[2] - eval1;
	T m01 = U[0] * AV[0] + U[1] * AV[1] + U[2] * AV[2];
	T m11 = V[0] * AU[0] + V[1] * AU[1] + V[2] * AU[2] - eval1;

    T absM00 = std::fabs(m00);
    T absM01 = std::fabs(m01);
    T absM11 = std::fabs(m11);
    T maxAbsComp;

    // Compare absM00 and absM11 first
    if (absM00 >= absM11) {
        maxAbsComp = (absM00 >= absM01) ? absM00 : absM01;
        if (maxAbsComp > (T)0) {
            if (absM00 >= absM01) {
                m01 /= m00;
                m00 = (T)1 / std::sqrt((T)1 + m01 * m01);
                m01 *= m00;
            }
            else {
                m00 /= m01;
                m01 = (T)1 / std::sqrt((T)1 + m00 * m00);
                m00 *= m01;
            }
            // Subtract(Multiply(m01, U), Multiply(m00, V))
            evec1[0] = m01 * U[0] - m00 * V[0];
            evec1[1] = m01 * U[1] - m00 * V[1];
            evec1[2] = m01 * U[2] - m00 * V[2];

        }
        else evec1 = U;
    }
    else {
        maxAbsComp = (absM11 >= absM01) ? absM11 : absM01;
        if (maxAbsComp > (T)0) {
            if (absM11 >= absM01) {
                m01 /= m11;
                m11 = (T)1 / std::sqrt((T)1 + m01 * m01);
                m01 *= m11;
            }
            else {
                m11 /= m01;
                m01 = (T)1 / std::sqrt((T)1 + m11 * m11);
                m11 *= m01;
            }
            // Subtract(Multiply(m11, U), Multiply(m01, V))
            evec1[0] = m11 * U[0] - m01 * V[0];
            evec1[1] = m11 * U[1] - m01 * V[1];
            evec1[2] = m11 * U[2] - m01 * V[2];

        }
        else evec1 = U;
    }
}
/// Calculate the eigenvalues of a 2x2 matrix
template<typename T>
CUDA_CALLABLE void eval2_symmetric(const T a, const T b, const T c, T& eval0, T& eval1) {
    // | a  b |
    // | b  c |


    // this algorithm uses the trace method to calculate the eigenvalues
    T tr = a + c; //trace2(matrix);              // calculate the matrix trace
    T det = a * c - b * b; //determinant2(matrix);       // calculate the determinant

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

    //const float a = matrix[0];
    //const float b = matrix[1];

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

/// Calculate the eigenvalues of a 3x3 matrix
template<typename T>
CUDA_CALLABLE void eval3_symmetric(const T a, const T b, const T c, const T d, const T e, const T f, 
    T& eval0, T& eval1, T& eval2) {
	// | a   b   d |
	// | b   c   e |
	// | d   e   f |

    // Case: matrix is diagonal
    T p1 = b * b + d * d + e * e;
    if (p1 == 0.0) {
        eval0 = a;
        eval1 = c;
        eval2 = f;
        if (eval0 > eval2) swap(eval0, eval2);
        if (eval0 > eval1) swap(eval0, eval1);
        if (eval1 > eval2) swap(eval1, eval2);
        return;
    }

    T q = (a + c + f) / 3.0;
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
    eval2 = q + 2.0 * p * cos(phi);
    eval0 = q + 2.0 * p * cos(phi + (2.0 * PI / 3.0));
    eval1 = 3.0 * q - eval2 - eval0;                // since trace(A) = eig1 + eig2 + eig3
}

/// Calculate the eigenvector of a 3x3 matrix associated with the eigenvalue lambda.
/// The result is returned in polar coordinates (theta, phi).
/// Avoiding the fragile cross-product method.
template<typename T>
CUDA_CALLABLE void evec3spherical_symmetric(const T a, const T b, const T c, const T d, const T e, const T f,
    const T* evals, T* evec0, T* evec1, T* evec2) {
    // Newer version (column-major)
    // | a  b  d |
    // | b  c  e |
    // | d  e  f |

	T determinant = determinanat3(a, b, c, d, e, f);
    if (determinant >= 0.0) {
		ComputeEigenvector0(a, b, c, d, e, f, evals[0], evec0);
		ComputeEigenvector1(a, b, c, d, e, f, evals[1], evec0, evec1);
		cross3(evec0, evec1, evec2);        // evec2 = evec0 x evec1
    }
    else {
		ComputeEigenvector0(a, b, c, d, e, f, evals[2], evec2);
		ComputeEigenvector1(a, b, c, d, e, f, evals[1], evec2, evec0);
		cross3(evec2, evec1, evec0);        // evec0 = evec2 x evec1
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
    T* eigenvalues2_symmetric(const T* mats, const size_t n) {

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

    /// <summary>
    /// CPU code for calculating eigenvectors of an array of 2x2 matrices given an
    /// array of corresponding eigenvalues.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <param name="mats">is an array of 2x2 matrices</param>
    /// <param name="evals">is an array of 2x1 eigenvalues</param>
    /// <param name="n">is the number of matrices and eigenvalue pairs in the array</param>
    /// <returns></returns>
    template<typename T>
    T* eigenvectors2polar_symmetric(const T* mats, const T* evals, const size_t n) {

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

    /// <summary>
    /// CPU code for calculating eigenvalues of a 3D matrix array
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <param name="mats"></param>
    /// <param name="n"></param>
    /// <returns></returns>
    template<typename T>
    T* eigenvalues3_symmetric(const T* mats, const size_t n) {

        T* evals = new T[3 * n];
        double eval0, eval1, eval2;
        for (size_t i = 0; i < n; i++) {
			double a = mats[i * 9 + 0];
			double b = mats[i * 9 + 1];
			double d = mats[i * 9 + 2];
			double c = mats[i * 9 + 4];
			double e = mats[i * 9 + 5];
			double f = mats[i * 9 + 8];

            eval3_symmetric(a,b,c,d,e,f, eval0, eval1, eval2);
            evals[i * 3 + 0] = eval0;
            evals[i * 3 + 1] = eval1;
            evals[i * 3 + 2] = eval2;
        }
        return evals;
    }

    template<typename T>
    T* eigenvectors3spherical_symmetric(const T* mats, const T* lambda, const size_t n) {
        // | a  b  d |
        // | b  c  e |
        // | d  e  f |
        //T* evecs = new T[4 * n];
        T* evecs = new T[9 * n];
        size_t i = 13;
        double a = mats[i * 9 + 0];
        double b = mats[i * 9 + 1];
        double d = mats[i * 9 + 2];
        double c = mats[i * 9 + 4];
        double e = mats[i * 9 + 5];
        double f = mats[i * 9 + 8];
		double* evec0 = new double[3];
		double* evec1 = new double[3];
		double* evec2 = new double[3];
		double evals[] = { lambda[i * 3 + 0], lambda[i * 3 + 1], lambda[i * 3 + 2] };
		double norm = b * b + d * d + e * e; // norm of the off-diagonal elements
        if (norm > 0.0)
			evec3spherical_symmetric(a, b, c, d, e, f, evals, evec0, evec1, evec2);
        else {
            evec0[0] = 1.0; evec0[1] = 0.0; evec0[2] = 0.0;
			evec1[0] = 0.0; evec1[1] = 1.0; evec1[2] = 0.0;
			evec2[0] = 0.0; evec2[1] = 0.0; evec2[2] = 1.0; // the matrix is diagonal, returns the standard basis vectors
        }

		evecs[i * 9 + 0] = evec0[0];
		evecs[i * 9 + 1] = evec0[1];
		evecs[i * 9 + 2] = evec0[2];        // first eigenvector
		evecs[i * 9 + 3] = evec1[0];
		evecs[i * 9 + 4] = evec1[1];
		evecs[i * 9 + 5] = evec1[2];        // second eigenvector
		evecs[i * 9 + 6] = evec2[0];
		evecs[i * 9 + 7] = evec2[1];
		evecs[i * 9 + 8] = evec2[2];        // third eigenvector
        return evecs;
    }
}

