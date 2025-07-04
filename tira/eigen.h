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
CUDA_CALLABLE void eval3_symmetric(const T* A, T& eval0, T& eval1, T& eval2) {

    T a, b, c, d, e, f, g, h, i;
    a = A[0 * 3 + 0]; b = A[0 * 3 + 1]; c = A[0 * 3 + 2];
    d = A[1 * 3 + 0]; e = A[1 * 3 + 1]; f = A[1 * 3 + 2];
    g = A[2 * 3 + 0]; h = A[2 * 3 + 1]; i = A[2 * 3 + 2];


    // Case: matrix is diagonal
    T p1 = b * b + c * c + f * f;
    if (p1 == 0.0) {
        eval0 = a;
        eval1 = e;
        eval2 = i;
        if (eval0 > eval2) swap(eval0, eval2);
        if (eval0 > eval1) swap(eval0, eval1);
        if (eval1 > eval2) swap(eval1, eval2);
        return;
    }

    T q = (a + e + i) / 3.0;
    T p2 = (a - q) * (a - q) + (e - q) * (e - q) + (i - q) * (i - q) + 2.0 * p1;
    T p = sqrt(p2 / 6.0);
    T pinv = 1.0 / p;
    T B[3][3];
    B[0][0] = pinv * (a - q);
    B[1][1] = pinv * (e - q);
    B[2][2] = pinv * (i - q);
    B[0][1] = pinv * b;
    B[0][2] = pinv * c;
    B[1][0] = pinv * d;
    B[1][2] = pinv * f;
    B[2][0] = pinv * g;
    B[2][1] = pinv * h;
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
    eval1 = 3.0 * q - eval2 - eval0;  // since trace(A) = eig1 + eig2 + eig3
}

/// Calculate the eigenvector of a 3x3 matrix associated with the eigenvalue lambda.
/// The result is returned in polar coordinates (theta, phi).
template<typename T>
CUDA_CALLABLE void evec3spherical_symmetric(const T* matrix, T lambda, T& theta, T& phi) {

    float a, b, c, d, e, f, g, h, i;
    a = matrix[0 * 3 + 0]; b = matrix[0 * 3 + 1]; c = matrix[0 * 3 + 2];
    d = matrix[1 * 3 + 0]; e = matrix[1 * 3 + 1]; f = matrix[1 * 3 + 2];
    g = matrix[2 * 3 + 0]; h = matrix[2 * 3 + 1]; i = matrix[2 * 3 + 2];

    // rows of (A - lambda*I)
    // all the rows multiplied by the eigenvector yield zero vector => eigenvector is prependicular to at least two of the rows
    T row0[] = { a - lambda, b, c };
    T row1[] = { d, e - lambda, f };
    T row2[] = { g, h, i - lambda };

    // calculate the cross-product of each two rows
    // v is parallel to the cross product of two of these rows
    T r0xr1[] = { row0[1] * row1[2] - row0[2] * row1[1],
        row0[2] * row1[0] - row0[0] * row1[2],
        row0[0] * row1[1] - row0[1] * row1[0] };

    T r0xr2[] = { row0[1] * row2[2] - row0[2] * row2[1],
        row0[2] * row2[0] - row0[0] * row2[2],
        row0[0] * row2[1] - row0[1] * row2[0] };

    T r1xr2[] = { row1[1] * row2[2] - row1[2] * row2[1],
        row1[2] * row2[0] - row1[0] * row2[2],
        row1[0] * row2[1] - row1[1] * row2[0] };

    // dot product - to find out which cross-product has the largest length
    float d0 = r0xr1[0] * r0xr1[0] + r0xr1[1] * r0xr1[1] + r0xr1[2] * r0xr1[2];
    float d1 = r0xr2[0] * r0xr2[0] + r0xr2[1] * r0xr2[1] + r0xr2[2] * r0xr2[2];
    float d2 = r1xr2[0] * r1xr2[0] + r1xr2[1] * r1xr2[1] + r1xr2[2] * r1xr2[2];
    int imax = 0;
    float dmax = d0;

    if (d1 > dmax) {
        dmax = d1;
        imax = 1;
    }
    if (d2 > dmax)
        imax = 2;

    T evec[3];
    if (imax == 0.0) {
        evec[0] = r0xr1[0] / sqrt(d0);
        evec[1] = r0xr1[1] / sqrt(d0);
        evec[2] = r0xr1[2] / sqrt(d0);
    }
    else if (imax == 1.0) {
        evec[0] = r0xr2[0] / sqrt(d1);
        evec[1] = r0xr2[1] / sqrt(d1);
        evec[2] = r0xr2[2] / sqrt(d1);
    }
    else {
        evec[0] = r1xr2[0] / sqrt(d2);
        evec[1] = r1xr2[1] / sqrt(d2);
        evec[2] = r1xr2[2] / sqrt(d2);
    }
    theta = acos(evec[2]);
    phi = atan2(evec[1], evec[0]);

};



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
        T eval0, eval1, eval2;
        for (size_t i = 0; i < n; i++) {
            eval3_symmetric(&mats[i * 9], eval0, eval1, eval2);
            evals[i * 3 + 0] = eval0;
            evals[i * 3 + 1] = eval1;
            evals[i * 3 + 2] = eval2;
        }
        return evals;
    }

    template<typename T>
    T* eigenvectors3spherical_symmetric(const T* mats, const T* lambda, const size_t n) {

        T* evec = new T[4 * n];
        for (size_t i = 0; i < n; i++) {
            T theta, phi;
            evec3spherical_symmetric(&mats[i * 9], lambda[i * 3 + 1], theta, phi);
            if (isnan(theta) || isnan(phi)) {
                evec[i * 4 + 0] = acos(0);
                evec[i * 4 + 1] = atan2(1, 0);
            }
            else {
                evec[i * 4 + 0] = theta;
                evec[i * 4 + 1] = phi;
            }

            evec3spherical_symmetric(&mats[i * 9], lambda[i * 3 + 2], theta, phi);

            if (isnan(theta) || isnan(phi)) {
                evec[i * 4 + 2] = acos(1);
                evec[i * 4 + 3] = atan2(0, 0);
            }
            else {
                evec[i * 4 + 2] = theta;
                evec[i * 4 + 3] = phi;
            }
        }
        return evec;
    }
}

