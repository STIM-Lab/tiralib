#pragma once

#include <tira/cuda/callable.h>

#ifdef __CUDACC__
    #include <tira/cuda/error.h>
#endif

#include "tensor.h"
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

#define TV_PI 3.14159265358979323846
#define TV_EPSILON 1e-12
#ifndef TV3_MAX_CONST_NB
#define TV3_MAX_CONST_NB 1536
#endif

namespace tira {
    /**
     * Calculates the decay value for a stick tensor voting field
     *
     * @tparam Type data type used to calculate the decay value
     * @param sin_theta angular term used to determine the angular decay, with maximum values in directions orthogonal to the stick tensor
     * @param length distance between the voter and receiver
     * @param sigma standard deviation indicating the attenuation of the vote signal
     * @return a decay value used to scale a stick tensor vote
     */
    template <typename Type>
    CUDA_CALLABLE static Type stick_decay(Type sin_theta, Type length, Type sigma) {
        Type c = exp(-(length * length) / (sigma * sigma));
        Type tp = sin_theta;
        return c * tp;
    }

    /**
     * Calculate the distance attentuation for the tensor vote, which is just a Gaussian
     * function based on the distance between the voter and receiver.
     *
     * @tparam Type data type for the calculation
     * @param dist distance between the voter and receiver
     * @param sigma standard deviation for the decay
     * @return the magnitude of the decay function
     */
    template <typename Type>
    CUDA_CALLABLE static Type atv_attenuation(Type dist, Type sigma) {
        if (sigma == Type(0)) {
            return (dist > Type(TV_EPSILON)) ? Type(0) : Type(1);
        }
        return exp(-(dist * dist) / (sigma * sigma));
    }

    /**
     * Calculate the rotation matrix R for analytical tensor voting. The rotation matrix
     * is returned as a diagonal matrix in column-major format:
     * [ a  b ]
     * [ b  c ]
     *
     * @tparam Type data type for the calculation
     * @param dx is the x direction from the voter to the receiver
     * @param dy is the y direction from the voter to the receiver
     * @param out_a is the upper-left component of the rotation matrix
     * @param out_b is the upper-right component of the rotation matrix
     * @param out_c is the lower-right component of the rotation matrix
     */
    template <typename Type>
    CUDA_CALLABLE static void atv_rotation(Type dx, Type dy, Type& out_a, Type& out_b, Type& out_c) {
        out_a = 1 - 2 * dx * dx;
        out_b = 0 - 2 * dx * dy;
        out_c = 1 - 2 * dy * dy;
    }

    /**
     * Calculate the tensor produced by a normalized eigenvector. This is the outer product
     * of the vector with itself. The matrix is returned as a diagonal matrix in
     * column-major format:
     * [ a  b ]
     * [ b  c ]
     *
     * @tparam Type data type for the calculation
     * @param evec_x x direction of the eigenvector
     * @param evec_y y direction of the eigenvector
     * @param out_a
     * @param out_b
     * @param out_c
     */
    template <typename Type>
    CUDA_CALLABLE static void atv_normalized_T(Type evec_x, Type evec_y,
        Type& out_a, Type& out_b, Type& out_c) {
        out_a = evec_x * evec_x;
        out_b = evec_x * evec_y;
        out_c = evec_y * evec_y;
    }

    /*
     * The following functions are used to calculate the H matrix in the analytical tensor
     * voting paper.
     */

    /**
     * Calculates the v^T * T_d * v * T_d term of the H matrix given the direction
     * vector v and the normalized matrix T_d.
     * The result is a 2x2 diagonal matrix in column-major form:
     * [ a  b ]
     * [ b  c ]
     *
     * @tparam Type data type used in the calculation
     * @param dx is the x direction from voter to receiver
     * @param dy is the y direction from voter to receiver
     * @param T_a is the upper-left component of T_d
     * @param T_b is the upper-right component of T_d
     * @param T_c is the lower-right component of T_d
     * @param out_a is the upper-left component of the output matrix
     * @param out_b is the upper-right component of the output matrix
     * @param out_c is the lower-right component of the output matrix
     */
    template <typename Type>
    CUDA_CALLABLE void atv_vTvT(Type dx, Type dy, Type T_a, Type T_b, Type T_c,
        Type& out_a, Type& out_b, Type& out_c) {

        Type vTv = dx * (dx * T_a + dy * T_b) + dy * (dx * T_b + dy * T_c);
        out_a = T_a * vTv;
        out_b = T_b * vTv;
        out_c = T_c * vTv;
    }

    /**
     * Calculates the T_d * v * v^T * T_d term of the H matrix given the direction
     * vector v and the normalized matrix T_d.
     * The result is a 2x2 diagonal matrix in column-major form:
     * [ a  b ]
     * [ b  c ]
     * @tparam Type data type used in the calculation
     * @param dx is the x direction from voter to receiver
     * @param dy is the y direction from voter to receiver
     * @param T_a is the upper-left component of T_d
     * @param T_b is the upper-right component of T_d
     * @param T_c is the lower-right component of T_d
     * @param out_a is the upper-left component of the output matrix
     * @param out_b is the upper-right component of the output matrix
     * @param out_c is the lower-right component of the output matrix
     */
    template <typename Type>
    CUDA_CALLABLE void atv_TvvT(Type dx, Type dy, Type T_a, Type T_b, Type T_c,
        Type& out_a, Type& out_b, Type& out_c) {

        /*Type vv = dx * dx + dy * dy;
        Type aa = T_a * T_a;
        Type bb = T_b * T_b;
        Type cc = T_c * T_c;
        Type ab = T_a * T_c;
        Type bc = T_b * T_c;
        out_a = (aa + bb) * vv;
        out_b = (ab + bc) * vv;
        out_c = (bb + cc) * vv;*/
        
        // v = [dx, dy]^T
        // Tv = T * v
        const Type tvx = T_a * dx + T_b * dy;
        const Type tvy = T_b * dx + T_c * dy;

        // T v v^T T = (T v) (T v)^T
        out_a = tvx * tvx;
        out_b = tvx * tvy;
        out_c = tvy * tvy;
    }

    /**
     * Calculates the H matrix for analytical tensor voting.
     * The result is a 2x2 diagonal matrix in column-major form:
     * [ a  b ]
     * [ b  c ]
     *
     * @tparam Type  data type used in the calculation
     * @param lambda0
     * @param lambda1
     * @param ev0x x direction of the smallest eigenvector
     * @param ev0y y direction of the smallest eigenvector
     * @param ev1x x direction of the largest eigenvector
     * @param ev1y y direction of the largest eigenvector
     * @param dx is the x direction from voter to receiver
     * @param dy is the y direction from voter to receiver
     * @param H_a upper-left component of H
     * @param H_b upper-right component of H
     * @param H_c lower-right component of H
     */
    template <typename Type>
    CUDA_CALLABLE static void atv_H(Type lambda0, Type lambda1,
        Type ev0x, Type ev0y, Type ev1x, Type ev1y, Type dx, Type dy,
        Type& H_a, Type& H_b, Type& H_c) {

        // Calculate the normalized tensors T0 and T1
        Type T1_a, T1_b, T1_c;
        atv_normalized_T(ev1x, ev1y, T1_a, T1_b, T1_c);

        Type T0_a, T0_b, T0_c;
        atv_normalized_T(ev0x, ev0y, T0_a, T0_b, T0_c);
        T0_a += T1_a;
        T0_b += T1_b;
        T0_c += T1_c;

        // calculate the internal terms vTvT for each matrix
        Type vT0vT0_a, vT0vT0_b, vT0vT0_c;
        atv_vTvT(dx, dy, T0_a, T0_b, T0_c, vT0vT0_a, vT0vT0_b, vT0vT0_c);
        Type vT1vT1_a, vT1vT1_b, vT1vT1_c;
        atv_vTvT(dx, dy, T1_a, T1_b, T1_c, vT1vT1_a, vT1vT1_b, vT1vT1_c);

        // calculate the internal terms TvvT for each matrix
        Type T0vvT0_a, T0vvT0_b, T0vvT0_c;
        atv_TvvT(dx, dy, T0_a, T0_b, T0_c, T0vvT0_a, T0vvT0_b, T0vvT0_c);
        Type T1vvT1_a, T1vvT1_b, T1vvT1_c;
        atv_TvvT(dx, dy, T1_a, T1_b, T1_c, T1vvT1_a, T1vvT1_b, T1vvT1_c);

        Type s1 = lambda1 - lambda0;
        const Type inv3 = Type(1.0) / Type(3.0);        // K = 1
        H_a = s1 * (T1_a - inv3 * (vT1vT1_a + 2 * T1vvT1_a));
        H_b = s1 * (T1_b - inv3 * (vT1vT1_b + 2 * T1vvT1_b));
        H_c = s1 * (T1_c - inv3 * (vT1vT1_c + 2 * T1vvT1_c));

        Type s0 = lambda0;
        const Type inv4 = Type(1.0) / Type(4.0);        // K = 2
        H_a += s0 * (T0_a - inv4 * (vT0vT0_a + 2 * T0vvT0_a));
        H_b += s0 * (T0_b - inv4 * (vT0vT0_b + 2 * T0vvT0_b));
        H_c += s0 * (T0_c - inv4 * (vT0vT0_c + 2 * T0vvT0_c));
    }

    /**
     * Implements analytical tensor voting given:
     *      1) a pair of eigenvalues such that l0 <= l1
     *      2) a pair of eigenvectors associated with these eigenvalues
     *      3) the position of the receiver relative to the voter (where the voter is at the origin)
     *      4) a sigma value specifying the distance attenuation of the vote
     * This function calculates the resulting vote at the receiver location as a
     * 2x2 diagonal matrix in column-major form:
     * [ a  b ]
     * [ b  c ]
     * @tparam Type data type used in the calculation
     * @param lambda0 smallest eigenvalue
     * @param lambda1 largest eigenvalue
     * @param ev0x x coordinate of the smallest eigenvector
     * @param ev0y y coordinate of the smallest eigenvector
     * @param ev1x x coordinate of the largest eigenvector
     * @param ev1y y coordinate of the largest eigenvector
     * @param rx x coordinate of the receiver relative to the voter location
     * @param ry y coordinate of the receiver relative to the voter location
     * @param sigma distance attenuation parameter (standard deviation of a Gaussian)
     * @param V_a upper-left component of the resulting vote matrix
     * @param V_b upper-right component of the resulting vote matrix
     * @param V_c lower-right component of the resulting vote matrix
     */
    template <typename Type>
    CUDA_CALLABLE static void atv(Type lambda0, Type lambda1,
        Type ev0x, Type ev0y, Type ev1x, Type ev1y, Type rx, Type ry, Type sigma,
        Type& V_a, Type& V_b, Type& V_c) {

        // calculate the distance attenuation function
        Type dist = std::sqrt(rx * rx + ry * ry);
        Type c = atv_attenuation(dist, sigma);

        // calculate the normalized direction vector from voter to receiver
        Type dx, dy;
        if (dist == 0) {
            dx = 0;
            dy = 0;
        }
        else {
            dx = rx / dist;
            dy = ry / dist;
        }


        // calculate the rotation matrix
        Type R_a, R_b, R_c;
        atv_rotation(dx, dy, R_a, R_b, R_c);

        // calculate the H matrix
        Type H_a, H_b, H_c;
        atv_H(lambda0, lambda1, ev0x, ev0y, ev1x, ev1y, dx, dy, H_a, H_b, H_c);

        // calculate the final tensor vote
        Type Ra2 = R_a * R_a;
        Type Rb2 = R_b * R_b;
        Type Rc2 = R_c * R_c;
        Type Rab = R_a * R_b;
        Type Rbc = R_b * R_c;
        Type Rac = R_a * R_c;
        V_a = c * (Ra2 * H_a    + 2 * Rab * H_b     + Rb2 * H_c);
        V_b = c * (Rb2 * H_b    + Rab * H_a         + Rac * H_b     + Rbc * H_c);
        V_c = c * (Rb2 * H_a    + 2 * Rbc * H_b     + Rc2 * H_c);
    }

    /**
     * Calculate the result of analytical tensor voting across a window of voters, where
     * the receiver is at the center. This function sums the contributions of all votes
     * within the window and returns the final vote result at the receiver. The receiver
     * is a 2x2 diagonal matrix in column-major form:
     * [ a  b ]
     * [ b  c ]
     *
     * @tparam Type is the data type for the calculation
     * @param lambdas a pointer to an array of eigenvalues for each sample point in an image
     * @param thetas a pointer to an array of theta values expressing the corresponding
     * eigenvector directions in polar coordinates
     * @param sigma is the decay term describing the fall-off of the vote contribution with distance
     * @param w is the size of the window (usually some factor of sigma)
     * @param s0 image size along the fast axis (usually x)
     * @param s1 image size along the slow axis (usually y)
     * @param x0 first coordinate for the receiver position (center of the window)
     * @param x1 second coordinate for the receiver position (center of the window)
     * @param out_a upper-left component of the summed receiver vote
     * @param out_b upper-right component of the summed receiver vote
     * @param out_c lower-right component of the summed receiver vote
     */
    template <typename Type>
    CUDA_CALLABLE static void atv_window(const Type* lambdas, const Type* thetas, const Type sigma,
        const int w, const int s0, const int s1, const int x0, const int x1,
        Type& out_a, Type& out_b, Type& out_c) {

        // Indices for the symmetric tensor:
        // | a  b |
        // | b  c |
        out_a = out_b = out_c = 0;                           // initialize the receiver value to zero

        const int hw = w / 2;                                // calculate the half window size
        for (int v = -hw; v <= hw; v++) {                    // for each pixel in the window
            const int r0 = static_cast<int>(x0) + v;         // calculate the position relative to the voter
            if (r0 >= 0 && r0 < s0) {
                for (int u = -hw; u <= hw; u++) {

                    const int r1 = static_cast<int>(x1) + u;
                    if (r1 >= 0 && r1 < s1) {
                        const int idx = r0 * s1 * 2 + r1 * 2;
                        const Type theta0 = thetas[idx + 0];               // retrieve the smallest eigenvector
                        const Type evx0 = std::cos(theta0);
                        const Type evy0 = std::sin(theta0);
                        const Type theta1 = thetas[idx + 1];
                        const Type evx1 = std::cos(theta1);
                        const Type evy1 = std::sin(theta1);

                        Type Va, Vb, Vc;

                        //tira::atv<Type>(u, v, sigma, theta, a, b, c); // calculate the stick vote contribution
                        const Type l0 = lambdas[idx + 0];                 // load both eigenvalues associated with the voter
                        const Type l1 = lambdas[idx + 1];
                        tira::atv<Type>(l0, l1, evx0, evy0, evx1, evy1, u, v, sigma, Va, Vb, Vc);

                        out_a += Va;
                        out_b += Vb;
                        out_c += Vc;
                    }
                }
            }
        }
    }

    /**
     * Compute the angular-and-distance modulation (eta) for our stick vote kernel.
     *
     * Kernel:
     *  eta = c1 * (t1^power) + c2 * (t2^power)
     *  where c1 = exp(-dist2 / sigma1^2), c2 = exp(-dist2 / sigma2^2)
     *
     * @tparam Type data type
     * @param dist2 squared distance between voter and receiver
     * @param qTd   dot product between stick direction and voter-to-receiver direction
     * @param sigma1 primary sigma
     * @param sigma2 secondary sigma
     * @param power refinement exponent (use >= 1; if 0, treated as 1)
     * @return eta
     */
    template <typename Type>
    CUDA_CALLABLE static Type tk_stick_eta_2d(Type dist, Type qTd, Type sigma1, Type sigma2, unsigned power) {

        if (power < 1) power = 1;

        const Type c1 = atv_attenuation(dist, sigma1);
        const Type c2 = atv_attenuation(dist, sigma2);

        const Type qTd2 = qTd * qTd;
        const Type t1 = Type(1) - qTd2;
        const Type t2 = qTd2;

        if (power == 1) {
            return c1 * t1 + c2 * t2;
        }

        // Integer-power evaluation (avoids std::pow and matches the style used in your 3D code)
        Type e1 = t1;
        Type e2 = t2;
        for (unsigned i = 1; i < power; ++i) {
            e1 *= t1;
            e2 *= t2;
        }
        return c1 * e1 + c2 * e2;
    }
    
    /**
     * Compute the 2D stick vote tensor contribution at a receiver position relative
     * to the voter, using our tensor-kernel (tk_) stick formulation.
     *
     * Output is a 2x2 symmetric tensor in column-major form:
     * [ Va  Vb ]
     * [ Vb  Vc ]
     *
     * @tparam Type data type used in the calculation
     * @param rx     receiver x relative to voter
     * @param ry     receiver y relative to voter
     * @param sigma1 primary sigma (orthogonal channel)
     * @param sigma2 secondary sigma (colinear channel)
     * @param power  refinement exponent
     * @param theta  stick direction angle (largest eigenvector) in polar coordinates
     * @param Va     upper-left output component
     * @param Vb     upper-right output component
     * @param Vc     lower-right output component
     */
    template <typename Type>
    CUDA_CALLABLE static void tk_stickvote_2d(const Type rx, const Type ry,
        const Type sigma1, const Type sigma2, const unsigned power,
        const Type theta, Type& Va, Type& Vb, Type& Vc)
    {
        // Stick direction q from theta
        const Type qx = cos(theta);
        const Type qy = sin(theta);

        // Distance and normalized direction d from voter to receiver
        const Type dist = std::sqrt(rx * rx + ry * ry);

        Type dx, dy;
        if (dist == Type(0)) {
            dx = Type(0);
            dy = Type(0);
        } else {
            dx = rx / dist;
            dy = ry / dist;
        }

        // Dot product q · d
        const Type qTd = qx * dx + qy * dy;

        // Our-method decay / refinement term
        const Type eta = tk_stick_eta_2d(dist, qTd, sigma1, sigma2, power);

        // Reflection/rotation matrix R = I - 2 d d^T
        Type R_a, R_b, R_c;
        atv_rotation(dx, dy, R_a, R_b, R_c);

        // Rotate/reflect the stick direction to get the vote orientation
        const Type Rq_x = R_a * qx + R_b * qy;
        const Type Rq_y = R_b * qx + R_c * qy;

        // Output vote tensor: (Rq)(Rq)^T scaled by eta
        Va = Rq_x * Rq_x * eta;
        Vb = Rq_x * Rq_y * eta;
        Vc = Rq_y * Rq_y * eta;
    }

    /**
     * Calculate the accumulated value of a tensor voting receiver by integrating votes across a window surrounding
     * the receiver position
     *
     * * The receiver vote is returned as a 2x2 symmetric tensor in column-major form:
     * [ out_a  out_b ]
     * [ out_b  out_c ]
     *
     * @tparam Type data type used for the calculation
     * @param lambdas array of eigenvalues (2 per pixel) corresponding to tensors in the image
     * @param thetas array of eigenvectors in polar coordinates (2 per pixel) corresponding to tensors in the image
     * @param sigma1 primary sigma (orthogonal channel)
     * @param sigma2 secondary sigma (colinear channel)
     * @param w window size used for summing (normally about 6x sigma)
     * @param s0 size of the tensor field along the first (slow) dimension
     * @param s1 size of the tensor field along the second (fast) dimension
     * @param x0 position of the receiver
     * @param x1 position of the receiver
     * @param out_a upper-left  output component
     * @param out_b upper-right output component
     * @param out_c lower-right output component
     */
    template <typename Type>
    CUDA_CALLABLE static void tk_stick_window_2d(const Type* lambdas, const Type* thetas, const Type sigma1, const Type sigma2,
        const unsigned power, const int w, const int s0, const int s1, const int x0, const int x1,
        Type& out_a, Type& out_b, Type& out_c) {

        out_a = out_b = out_c = Type(0);                                    // initialize the receiver value to zero

        const int hw = w / 2;                                               // calculate the half window size
        for (int v = -hw; v <= hw; v++) {                                   // for each pixel in the window
            const int r0 = x0 + v;                                          // calculate the position relative to the voter
            if (r0 >= 0 && r0 < s0) {
                for (int u = -hw; u <= hw; u++) {
                    const int r1 = x1 + u;
                    if (r1 >= 0 && r1 < s1) {
                        const int idx = r0 * s1 * 2 + r1 * 2;
                        const Type theta = thetas[idx + 1];               // retrieve the largest eigenvector

                        Type Va, Vb, Vc;
                        tk_stickvote_2d<Type>((Type)u, (Type)v, sigma1, sigma2, power, theta, Va, Vb, Vc); // calculate the stick vote contribution

                        const Type l0 = lambdas[idx + 0];                 // load both eigenvalues associated with the voter
                        const Type l1 = lambdas[idx + 1];

                        Type scale = fabs(l1) - fabs(l0);                              // calculate the vote scale based on the difference between eigenvalues
                        if (l1 < Type(0)) scale = scale * (-1);                        // TODO: Probably don't need this test since l1 should be larger than l0
                        out_a += scale * Va;                                           // accumulate the receiver vote
                        out_b += scale * Vb;
                        out_c += scale * Vc;
                    }
                }
            }
        }
    }

    /**
     * Calculate the 2D plate tensor contribution using numerical integration.
     *
     * @tparam Type data type used for the calculation
     * @param u is the position of the reciever relative to the tensor casting the vote (voter)
     * @param v is the position of the reciever relative to the tensor casting the vote (voter)
     * @param sigma is the attenuation of the tensor vote field
     * @param out_a upper left value in the symmetric 2x2 output matrix
     * @param out_b upper right value in the symmetric 2x2 output matrix
     * @param out_c lower right value in the symmetric 2x2 output matrix
     * @param n number of samples to use for numerical integration
     */
    template <typename Type>
    CUDA_CALLABLE  static void tk_platevote_numerical(const Type rx, const Type ry, const Type sigma1, const Type sigma2, const unsigned power,
        Type& out_a, Type& out_b, Type& out_c, const unsigned int n = 10) {

        out_a = out_b = out_c = 0;
        if (n == 0) return;

        const Type dtheta = TV_PI / static_cast<Type>(n);

        for (unsigned int i = 0; i < n; i++) {
            const Type theta = dtheta * static_cast<Type>(i);

            Type a, b, c;
            tk_stickvote_2d<Type>(rx, ry, sigma1, sigma2, power, theta, a, b, c);

            out_a += a;
            out_b += b;
            out_c += c;
        }
        const Type inv_n = Type(1) / static_cast<Type>(n);
        out_a *= inv_n;
        out_b *= inv_n;
        out_c *= inv_n;
    }

    /**
     * Accumulate the 2D plate vote tensor kernels at a receiver by integrating TK plate votes over a local window.
     * The contribution from each voter is scaled by lambda0, the smallest eigenvalue of the tensor at that location.
     *
     * @tparam Type data type used for the calculation
     * @param lambdas array of eigenvalues (2 per pixel) corresponding to tensors in the image
     * @param sigma attenuation value for the voting field
     * @param w window size used for summing (normally about 6x sigma)
     * @param s0 size of the tensor field along the first (slow) dimension
     * @param s1 size of the tensor field along the second (fast) dimension
     * @param x0 position of the receiver
     * @param x1 position of the receiver
     * @param out_a upper left value of the 2x2 symmetric tensor vote
     * @param out_b upper right value of the 2x2 symmetric tensor vote
     * @param out_c lower right value of the 2x2 symmetric tensor vote
     * @param samples number of stick tensor samples used for numerical integration
     */
    template <typename Type>
    CUDA_CALLABLE static void tk_platevote_window(const Type* lambdas, Type sigma1, Type sigma2, const unsigned power,
        const int w, const int s0, const int s1, const int x0, const int x1,
        Type& out_a, Type& out_b, Type& out_c, const unsigned samples = 10) {

        out_a = out_b = out_c = Type(0);
        if (samples == 0) return;
        
        const int hw = w / 2;

        for (int v = -hw; v <= hw; v++) {                    
            const int r0 = x0 + v;
            if (r0 >= 0 && r0 < s0) {
                for (int u = -hw; u <= hw; u++) {
                    int r1 = x1 + u;
                    if (r1 >= 0 && r1 < s1) {

                        const int idx = r0 * s1 * 2 + r1 * 2;
                        const Type l0 = lambdas[idx + 0];

                        if (l0 != Type(0)) {
                            Type Va, Vb, Vc;
                            tk_platevote_numerical<Type>(u, v, sigma1, sigma2, power, Va, Vb, Vc, samples);

                            const Type scale = fabsf(l0);
                            out_a += scale * Va;
                            out_b += scale * Vb;
                            out_c += scale * Vc;
                        }
                    }
                }
            }
        }
    }


    /**
    * @brief Precomputed 3D neighbor offset and weights for voting.
    */
    struct Neighbor3D {
            int du, dv, dw;                 // index offsets along x2 (u), x1 (v), and x0 (w)
            glm::vec3 d;                    // normalized direction from sender to receiver (du, dv, dw)
            float l2;                       // squared length of offset: du^2 + dv^2 + dw^2
            float c1;                       // exp(-l2 / sigma1^2)
            float c2;                       // exp(-l2 / sigma2^2)
        };
    
    /**
    * @brief Normalization constant for 3D stick voting
    *
    * @tparam T        Floating point type
    * @param sigma1    Stick-vote sigma along first axis.
    * @param sigma2    Stick-vote sigma along second axis.
    * @param p         Stick-vote power (refinement exponent).
    * @return          Normalization factor for 3D stick votes.
    */
    template <typename Type>
    CUDA_CALLABLE static Type sticknorm3(const Type sigma1, const Type sigma2, const unsigned p) {
        Type pi_term = TV_PI * sqrt(TV_PI) / 2.0;
        Type sig1_3 = sigma1 * sigma1 * sigma1;
        Type sig2_3 = sigma2 * sigma2 * sigma2;
        Type num1 = pow(2, 2 * p + 1) * factorial(p) * factorial(p);
        Type num2 = 2.0;
        Type den1 = factorial(2 * p + 1);
        Type den2 = 2 * p + 1;
        return pi_term * ((sig1_3 * num1 / den1) + (sig2_3 * num2 / den2));
    }

    /**
    * @brief Accumulate one 3D stick vote contribution into a tensor.
    *
    * @param M      3x3 tensor accumulator (modified in place).
    * @param n      Precomputed neighbor direction and decay factors.
    * @param q      Local stick direction (largest eigenvector).
    * @param power  Stick-vote power (refinement exponent).
    * @param scale  Scale from eigenvalues at the voter.
    */
    template<typename Dir, typename Vec>
    CUDA_CALLABLE void stickvote3_accumulate_kernel(float& m00, float& m01, float& m02,
        float& m11, float& m12, float& m22,
        const Dir& d, const float c1, const float c2,
        const Vec& q, const unsigned power, const float scale)
    {
        const float qTd = q.x * d.x + q.y * d.y + q.z * d.z;
        const float qTd2 = qTd * qTd;

        float eta;
        if (power == 1)
            eta = c1 * (1.0f - qTd2) + c2 * qTd2;
        else {
            const float t1 = 1.0f - qTd2;
            const float t2 = qTd2;
            float e1 = t1;
            float e2 = t2;
            for (unsigned i = 1; i < power; ++i) {
                e1 *= t1;
                e2 *= t2;
            }
			eta = c1 * e1 + c2 * e2;
        }

        // Reflected direction
        const float rx = q.x - 2.0f * qTd * d.x;
        const float ry = q.y - 2.0f * qTd * d.y;
        const float rz = q.z - 2.0f * qTd * d.z;
        const float term = scale * eta;

        m00 += term * rx * rx; m01 += term * rx * ry; m02 += term * rx * rz;
        m11 += term * ry * ry; m12 += term * ry * rz; m22 += term * rz * rz;
    }

    /**
     * CPU namespace contains functions that are run completely on the host. All input and output pointers
     * are allocated on the host.
     */
    namespace cpu {

        /**
        * Performs analyticaltensor voting on a 2D image using eigenvalues and eigenvectors.
        * 
        * @tparam Type data type used for the calculation
        * @param t_out is a pointer to the array to be filled with the tensor voting result
        * @param lambdas array of eigenvalues (2 per pixel) corresponding to tensors in the image
        * @param evecs is an array of eigenvectors (2 per pixel) in polar coordinates
        * @param sigma attenuation value for the voting field
        * @param w window size used for summing (normally about 6x sigma)
        * @param shape0 size of the tensor field along the first (slow) dimension
        * @param shape1 size of the tensor field along the second (fast) dimension
        * @param stick is a boolean flag specifying calculation of the stick tensor vote
        * @param plate is a boolean flag specifying calculation of the plate tensor vote
        * @param samples number of stick tensor samples used for numerical integration
        */
        template <typename Type>
        static void tensorvote_atv(Type* t_out, const Type* lambdas, const Type* evecs, Type sigma, const unsigned w,
            const unsigned shape0, const unsigned shape1,
            const bool stick = true, const bool plate = true, const unsigned samples = 10) {

            Type out_a, out_b, out_c;


            //float a, b, c;
            for (int x0 = 0; x0 < shape0; x0++) {
                for (int x1 = 0; x1 < shape1; x1++) {
                    out_a = out_b = out_c = 0;
                    atv_window(lambdas, evecs, sigma, w, shape0, shape1, x0, x1, out_a, out_b, out_c);
                    const unsigned idx = x0 * shape1 * 4 + x1 * 4;
                    /*if (stick) {
                        stickvote_window(lambdas, evecs, sigma, w, shape0, shape1, x0, x1, a, b, c);
                        out_a += a;
                        out_b += b;
                        out_c += c;
                    }

                    if (plate) {
                        platevote_window(lambdas, sigma, w, shape0, shape1, x0, x1, a, b, c, samples);
                        out_a += a;
                        out_b += b;
                        out_c += c;
                    }*/
                    t_out[idx + 0] = out_a;
                    t_out[idx + 1] = out_b;
                    t_out[idx + 2] = out_b;
                    t_out[idx + 3] = out_c;
                }
            }
        }

        /**
        * Tensor-kernel tensor voting on 2D images using eigenvalues and eigenvectors.
        * 
        * @tparam Type data type used for the calculation
        * @param t_out is a pointer to the array to be filled with the tensor voting result
        * @param lambdas array of eigenvalues (2 per pixel) corresponding to tensors in the image
        * @param evecs is an array of eigenvectors (2 per pixel) in polar coordinates
        * @param sigma1 primary sigma (orthogonal channel)
        * @param sigma2 secondary sigma (colinear channel)
        * @param power the refinement term
        * @param w window size used for summing (normally about 6x sigma)
        * @param shape0 size of the tensor field along the first (slow) dimension
        * @param shape1 size of the tensor field along the second (fast) dimension
        * @param stick is a boolean flag specifying calculation of the stick tensor vote
        * @param plate is a boolean flag specifying calculation of the plate tensor vote
        * @param samples number of stick tensor samples used for numerical integration
        */
        template <typename Type>
        static void tensorvote_tk(Type* t_out, const Type* lambdas, const Type* evecs, Type sigma1, Type sigma2, unsigned power,
            const unsigned w, const unsigned shape0, const unsigned shape1,
            const bool stick = true, const bool plate = true, const unsigned samples = 10) {

            Type out_a, out_b, out_c;

            for (int x0 = 0; x0 < shape0; x0++) {
                for (int x1 = 0; x1 < shape1; x1++) {
                    out_a = out_b = out_c = Type(0);

                    Type a, b, c;
                    if (stick) {
                        tk_stick_window_2d(lambdas, evecs, sigma1, sigma2, power, w, shape0, shape1, x0, x1, a, b, c);
                        out_a += a;
                        out_b += b;
                        out_c += c;
                    }

                    if (plate) {
                        tk_platevote_window(lambdas, sigma1, sigma2, power, w, shape0, shape1, x0, x1, a, b, c, samples);
                        out_a += a;
                        out_b += b;
                        out_c += c;
                    }

                    const unsigned idx = x0 * shape1 * 4 + x1 * 4;
                    t_out[idx + 0] = out_a;
                    t_out[idx + 1] = out_b;
                    t_out[idx + 2] = out_b;
                    t_out[idx + 3] = out_c;
                }
            }
        }

        template <typename Type>
        static void tensorvote(Type* t_out, const Type* t_in, const Type sigma1, const Type sigma2, const int power, const unsigned w,
            const unsigned shape0, const unsigned shape1, 
            const bool stick = true, const bool plate = true, const unsigned samples = 10, const bool is_atv = false) {

            // Calculate the eigenvalues for each tensor in the field
            float* lambdas = evals2_symmetric(t_in, shape0 * shape1);

            // Calculate the eigenvectors for each tensor in the field
            float* evecs = evecs2polar_symmetric(t_in, lambdas, shape0 * shape1);

            // Perform tensor voting
            if (is_atv)
                tensorvote_atv(t_out, lambdas, evecs, sigma1, w, shape0, shape1, stick, plate, samples);
            else
                tensorvote_tk(t_out, lambdas, evecs, sigma1, sigma2, power, w, shape0, shape1, stick, plate, samples);
        }





        /**
         * @brief Build all 3D neighbor offsets and Gaussian weights for a window.
         *
         * @param w         Side length of cubic voting window.
         * @param sigma     (sigma1, sigma2) for orthogonal and lateral decay.
         * @return          Vector of precomputed 3D neighbors.
         */
        inline std::vector<Neighbor3D> build_neighbors3d(int w, glm::vec2 sigma) {
            std::vector<Neighbor3D> neighbors;
            neighbors.reserve(static_cast<size_t>(w) * w * w);

            const int hw = w / 2;
            const float invsig1 = (sigma.x > 0) ? 1.0f / (sigma.x * sigma.x) : 0.0f;
            const float invsig2 = (sigma.y > 0) ? 1.0f / (sigma.y * sigma.y) : 0.0f;

            for (int dw = -hw; dw <= hw; dw++) {                         // For each pixel in the window
                for (int dv = -hw; dv <= hw; dv++) {
                    for (int du = -hw; du <= hw; du++) {
                        Neighbor3D n;
                        n.du = du; n.dv = dv; n.dw = dw;
                        n.l2 = float(du * du + dv * dv + dw * dw);
                        if (n.l2 == 0.0f) n.d = glm::vec3(0.0f, 0.0f, 0.0f);
                        else n.d = glm::vec3(float(du), float(dv), float(dw)) / sqrtf(n.l2);

						// Degenerate cases - contribution at voter position when sigma is zero
                        n.c1 = sigma.x == 0.0f ? (n.l2 > TV_EPSILON ? 0.0f : 1.0f) : expf(-n.l2 * invsig1);
                        n.c2 = sigma.y == 0.0f ? (n.l2 > TV_EPSILON ? 0.0f : 1.0f) : expf(-n.l2 * invsig2);
                        neighbors.push_back(n);
                    }
                }
            }
            return neighbors;
        }
        
        /**
         * @brief Compute 3D stick vote tensor at a receiver voxel.
         *
         * Uses precomputed 3D neighbors and largest eigenvectors to accumulate
         * stick votes around a receiver location.
         *
         * @param L      Volume of eigenvalues (L0, L1, L2) per voxel.
         * @param Q      Volume of largest eigenvectors (stick direction) per voxel.
         * @param NB     Precomputed neighbor offsets and decay factors.
         * @param power  Stick-vote power (refinement exponent).
         * @param s0     Volume size along first dimension.
         * @param s1     Volume size along second dimension.
         * @param s2     Volume size along third dimension.
         * @param x      Receiver voxel location index (x0, x1, x2).
         * @return       3x3 stick vote tensor at the receiver.
         */
        inline glm::mat3 stickvote3(const glm::vec3* L, const glm::vec3* Q, const std::vector<Neighbor3D>& NB,
            const unsigned power, const unsigned s0, const unsigned s1, const unsigned s2, const glm::ivec3 x) {

            const int x0 = x[0];
            const int x1 = x[1];
            const int x2 = x[2];

            // Accumulate in symmetric 6-value form
			float m00 = 0.0f, m01 = 0.0f, m02 = 0.0f;
			float m11 = 0.0f, m12 = 0.0f, m22 = 0.0f;

            for (const auto& n : NB) {                         // For each pixel in the window
                const int r0 = x0 + n.dw;
                if ((unsigned)r0 >= s0) continue;
                const int r1 = x1 + n.dv;
                if ((unsigned)r1 >= s1) continue;
                const int r2 = x2 + n.du;
                if ((unsigned)r2 >= s2) continue;

                // Flat index of the voter
                const unsigned base = (unsigned)r0 * s1 * s2 + (unsigned)r1 * s2 + (unsigned)r2;

                // Stick tensor direction (largest eigenvector)
                const glm::vec3 q = Q[base];
                // Eigenvalues (L0 <= L1 <= L2) - L0 not used for stick votes
                const float l1 = L[base].y;
                const float l2 = L[base].z;

                // Calculate the contribution scale factor
                const float scale = std::copysignf(fabsf(l2) - fabsf(l1), l2);

                // Calculate the accumulated contribution of (du,dv,dw) to (x,y,z)
                stickvote3_accumulate_kernel(m00, m01, m02, m11, m12, m22, n.d, n.c1, n.c2, q, power, scale);
            }
			glm::mat3 Votee(0.0f);
			Votee[0][0] = m00; Votee[0][1] = m01; Votee[0][2] = m02;
			Votee[1][0] = m01; Votee[1][1] = m11; Votee[1][2] = m12;
			Votee[2][0] = m02; Votee[2][1] = m12; Votee[2][2] = m22;
            return Votee;
        }

        /**
         * @brief Closed-form 3D plate vote from a single neighbor direction.
         *
         * Uses hypergeometric and beta integrals to analytically integrate
         * stick votes around a plate normal.
         *
         * @param d      Normalized direction from voter to receiver.
         * @param c1     Gaussian weight for sigma1 (exponential term).
         * @param c2     Gaussian weight for sigma2 (exponential term).
         * @param power  Plate-vote power (refinement exponent).
         * @param evec0  Plate normal (smallest eigenvector).
         * @param K0     Precomputed beta integral for J0-like term.
         * @param K1     Precomputed beta integral for J1-like term.
         * @return       3x3 plate vote tensor contribution.
         */
        inline glm::mat3 platevote3_analytic(const glm::vec3& d, float c1, float c2, unsigned power,
            const glm::vec3& evec0, double K0, double K1) {

            glm::vec3 dn = d;

            // Building the local coordinate system (rotate so that plate normal aligns with z-axis)
            glm::vec3 u, v;
            if (std::fabs(evec0.z) < 0.999f)                            // if evec0 is not aligned with z
                u = glm::normalize(glm::cross(glm::vec3(0.0f, 0.0f, 1.0f), evec0));
            else                                                        // if evec0 is already aligned with z, cross with x
                u = glm::normalize(glm::cross(glm::vec3(1.0f, 0.0f, 0.0f), evec0));
            v = glm::cross(evec0, u);

            const glm::mat3 Z(u, v, evec0);                             // rotation matrix to align plate normal with z-axis
            const glm::mat3 Z_trans = glm::transpose(Z);

            const glm::vec3 d_local = Z_trans * dn;                     // now d is aligned with z axis
            const float dx = d_local.x, dy = d_local.y, dz = d_local.z;

            // Calculate the length and angle of the reflected direction matrix on XY plane
            const float alpha = std::sqrt(dx * dx + dy * dy);
            const float a2 = alpha * alpha;
            const float phi = std::atan2(dy, dx);                       // angle of projection in the local xy plane

            // Compute beta and hypergeometric integrals
            const double p_d = static_cast<double>(power);
            const double a2_d = std::min(static_cast<double>(a2), 1.0);
            const double J0 = 0.5 * TV_PI * boost::math::hypergeometric_pFq(std::vector<double>{ -p_d, 1.5 },
                std::vector<double>{ 2.0 }, a2_d);
            const double J1 = TV_PI * boost::math::hypergeometric_pFq(std::vector<double>{ -p_d, 0.5 },
                std::vector<double>{ 1.0 }, a2_d);

            // Terms in the rotated frame
            glm::mat3 A(0.0f), B(0.0f);
            const float tmp_a2 = 1.0f - 2.0f * a2;

            A[0][0] = (tmp_a2 * tmp_a2) * static_cast<float>(J0);
            A[0][2] = -2.0f * alpha * dz * tmp_a2 * static_cast<float>(J0);
            A[2][0] = -2.0f * alpha * dz * tmp_a2 * static_cast<float>(J0);
            A[1][1] = static_cast<float>(J1 - J0);
            A[2][2] = 4.0f * a2 * (dz * dz) * static_cast<float>(J0);

            const float a2p = std::pow(a2, static_cast<float>(power));

            B[0][0] = a2p * (tmp_a2 * tmp_a2) * static_cast<float>(K0);
            B[0][2] = -2.0f * alpha * a2p * dz * tmp_a2 * static_cast<float>(K0);
            B[2][0] = -2.0f * alpha * a2p * dz * tmp_a2 * static_cast<float>(K0);
            B[1][1] = a2p * static_cast<float>(K1 - K0);
            B[2][2] = 4.0f * a2p * a2 * (dz * dz) * static_cast<float>(K0);

            // Initial coordinate system was rotated so the plate normal aligns with the z-axis
            // Build the rotation matrix around Z axis by +/- phi to rotate back to original coordinates
            glm::mat3 Rz(0.0f);
            const float cph = std::cos(phi), sph = std::sin(phi);
            Rz[0][0] = cph;   Rz[0][1] = -sph;
            Rz[1][0] = sph;   Rz[1][1] = cph;
            Rz[2][2] = 1.0f;
            const glm::mat3 Rz_rev = glm::transpose(Rz);

            // Rotate back to original coordinates in the local frame
            const glm::mat3 term_a = Rz_rev * A * Rz;
            const glm::mat3 term_b = Rz_rev * B * Rz;

            // Rotate from local back to global frame
            const glm::mat3 term_a_global = Z * term_a * Z_trans;
            const glm::mat3 term_b_global = Z * term_b * Z_trans;

            // Combine
            glm::mat3 PlateVote(0.0f);
            PlateVote[0][0] = c1 * term_a_global[0][0] + c2 * term_b_global[0][0];
            PlateVote[0][1] = c1 * term_a_global[0][1] + c2 * term_b_global[0][1];
            PlateVote[0][2] = c1 * term_a_global[0][2] + c2 * term_b_global[0][2];
            PlateVote[1][0] = c1 * term_a_global[1][0] + c2 * term_b_global[1][0];
            PlateVote[1][1] = c1 * term_a_global[1][1] + c2 * term_b_global[1][1];
            PlateVote[1][2] = c1 * term_a_global[1][2] + c2 * term_b_global[1][2];
            PlateVote[2][0] = c1 * term_a_global[2][0] + c2 * term_b_global[2][0];
            PlateVote[2][1] = c1 * term_a_global[2][1] + c2 * term_b_global[2][1];
            PlateVote[2][2] = c1 * term_a_global[2][2] + c2 * term_b_global[2][2];
            return PlateVote;
        }

        /**
         * @brief Numerical 3D plate vote via sampled stick integration.
         *
         * Integrates stick votes on a great circle around the plate normal
         * to approximate a plate vote tensor.
         *
         * @param d        Direction from voter to receiver.
         * @param c1       Gaussian weight for sigma1 (exponential term).
         * @param c2       Gaussian weight for sigma2 (exponential term).
         * @param power    Plate-vote power (refinement exponent).
         * @param evec0    Plate normal (smallest eigenvector).
         * @param samples  Number of angular samples for integration.
         * @return         Approximate 3x3 plate vote tensor.
         */
        inline glm::mat3 platevote3_numerical(const glm::vec3& d, float c1, float c2, unsigned power,
            const glm::vec3& evec0, unsigned samples = 20) {

            // A zero-vector has no orientation and cannot vote
            const float evec_len2 = evec0.x * evec0.x + evec0.y * evec0.y + evec0.z * evec0.z;
            glm::mat3 V(0.0f);
            if (samples == 0 or !(evec_len2 > TV_EPSILON)) return V;

            glm::vec3 dn = d;

            // Build an orthonomal basis (u,v) spanning the plane perpendicular to d
            glm::vec3 u, v;
            if (std::fabs(evec0.z) < 0.999f) {
                // Choose any vector not colinear to d
                u = glm::normalize(glm::cross(glm::vec3(0.0f, 0.0f, 1.0f), evec0));
            }
            else {
                u = glm::normalize(glm::cross(glm::vec3(1.0f, 0.0f, 0.0f), evec0));
            }
            v = glm::cross(evec0, u);

            // Integrate over [0, pi] to avoid double counting q and -q
            // (the stick vote is symmetric along the stick axis)
            const float dbeta = float(TV_PI) / float(samples);
			float m00 = 0.0f, m01 = 0.0f, m02 = 0.0f;
			float m11 = 0.0f, m12 = 0.0f, m22 = 0.0f;

            for (unsigned i = 0; i < samples; ++i) {
                const float beta = dbeta * float(i);
				const float cb = cosf(beta);
				const float sb = sinf(beta);
                const glm::vec3 q = cb * u + sb * v;

				stickvote3_accumulate_kernel(m00, m01, m02, m11, m12, m22, dn, c1, c2, q, power, 1.0f);
            }

            // Integral = sum / samples
			const float inv_samples = 1.0f / static_cast<float>(samples);
			m00 *= inv_samples; m01 *= inv_samples; m02 *= inv_samples;
			m11 *= inv_samples; m12 *= inv_samples; m22 *= inv_samples;
			V[0][0] = m00; V[0][1] = m01; V[0][2] = m02;
			V[1][0] = m01; V[1][1] = m11; V[1][2] = m12;
			V[2][0] = m02; V[2][1] = m12; V[2][2] = m22;
            
            return V;
        }

        /**
         * @brief 3D plate vote tensor at a receiver voxel.
         *
         * Combines either analytical or numerical plate votes from all
         * neighbors, scaled by (L1 - L0) at each voter.
         *
         * @param L        Volume of eigenvalues (L0, L1, L2) per voxel.
         * @param Q_small  Volume of smallest eigenvectors (plate normals) per voxel.
         * @param NB       Precomputed neighbor offsets and decay factors.
         * @param power    Plate-vote power (refinement exponent).
         * @param s0       Volume size along first dimension.
         * @param s1       Volume size along second dimension.
         * @param s2       Volume size along third dimension.
         * @param x        Receiver voxel index (x0, x1, x2).
         * @param samples  >0 for numerical integration, 0 for analytic form.
         * @return         3x3 plate vote tensor at the receiver.
         */
        inline glm::mat3 platevote3(const glm::vec3* L, const glm::vec3* Q_small, const std::vector<Neighbor3D>& NB, const unsigned power,
            const unsigned s0, const unsigned s1, const unsigned s2, const glm::ivec3 x, const unsigned samples = 0) {

            const int x0 = x[0];
            const int x1 = x[1];
            const int x2 = x[2];

            glm::mat3 Receiver(0.0f);

            const double p_d = static_cast<double>(power);
            const double K0 = boost::math::beta(0.5, p_d + 1.5);
            const double K1 = boost::math::beta(0.5, p_d + 0.5);

            for (const auto& n : NB) {                         // For each pixel in the window
                const int r0 = x0 + n.dw;
                if ((unsigned)r0 >= s0) continue;
                const int r1 = x1 + n.dv;
                if ((unsigned)r1 >= s1) continue;
                const int r2 = x2 + n.du;
                if ((unsigned)r2 >= s2) continue;

                const unsigned base = (unsigned)r0 * s1 * s2 + (unsigned)r1 * s2 + (unsigned)r2;

                const glm::vec3 evec0 = Q_small[base];
                const float l0 = L[base].x;
                const float l1 = L[base].y;
                float scale = std::copysign(std::abs(l1) - std::abs(l0), l1);

                const glm::vec3 d = n.d;
                const float c1 = n.c1, c2 = n.c2;
                glm::mat3 V;
                if (samples > 0)
                    // Numerical integration of stick votes to form a plate vote
                    V = platevote3_numerical(d, c1, c2, power, evec0, samples);
                else
                    // Analytical closed form solution from direction d
                    V = platevote3_analytic(d, c1, c2, power, evec0, K0, K1);

                Receiver[0][0] += scale * V[0][0]; Receiver[0][1] += scale * V[0][1]; Receiver[0][2] += scale * V[0][2];
                Receiver[1][0] += scale * V[1][0]; Receiver[1][1] += scale * V[1][1]; Receiver[1][2] += scale * V[1][2];
                Receiver[2][0] += scale * V[2][0]; Receiver[2][1] += scale * V[2][1]; Receiver[2][2] += scale * V[2][2];
            }
            return Receiver;
        }

        /**
         * @brief CPU implementation of full 3D tensor voting (stick + plate).
         *
         * For each voxel, accumulates stick and/or plate votes from a cubic
         * neighborhood using precomputed neighbors and eigen-frames.
         *
         * @param VT        Output volume of 3x3 vote tensors.
         * @param L         Volume of eigenvalues (L0, L1, L2) per voxel.
         * @param Q_large   Volume of largest eigenvectors (stick directions).
         * @param Q_small   Volume of smallest eigenvectors (plate normals).
         * @param sigma     (sigma1, sigma2) for stick/plate decay.
         * @param power     Voting power (refinement exponent).
         * @param w         Side length of cubic voting window.
         * @param s0        Volume size along first dimension.
         * @param s1        Volume size along second dimension.
         * @param s2        Volume size along third dimension.
         * @param STICK     If true, include stick votes.
         * @param PLATE     If true, include plate votes.
         * @param samples   Samples for numerical plate (0 = analytic).
         */
        inline void tensorvote3(glm::mat3* VT, const glm::vec3* L, const glm::vec3* Q_large, const glm::vec3* Q_small,
            glm::vec2 sigma, unsigned int power, const unsigned w, const unsigned s0, const unsigned s1, const unsigned s2,
            const bool STICK = true, const bool PLATE = true, const unsigned samples = 20) {
            const float sticknorm = 1.0 / sticknorm3(sigma.x, sigma.y, power);
            //const float platenorm = 1.0; // / TV_PI;           //  Not too sure about this
            // Pre-compute the neighbor offsets and Gaussian factors once for (w, sigmas)
            const auto NB = build_neighbors3d((int)w, sigma);

            // O(N * w^3) complexity
            for (unsigned int x0 = 0; x0 < s0; x0++) {
                for (unsigned int x1 = 0; x1 < s1; x1++) {
                    for (unsigned int x2 = 0; x2 < s2; x2++) {
                        glm::mat3 Vote(0.0f);
                        if (STICK)
                            Vote += sticknorm * stickvote3(L, Q_large, NB, power, s0, s1, s2, glm::ivec3(x0, x1, x2));
                        if (PLATE)
                            Vote += sticknorm * platevote3(L, Q_small, NB, power, s0, s1, s2, glm::ivec3(x0, x1, x2), samples);
                        VT[x0 * s1 * s2 + x1 * s2 + x2] = Vote;
                    }
                }
            }
        }
    }

    /**
     * The CUDA namespace runs everything on the GPU. Input and output pointers can be on either the host or device.
     * If pointers are located on the host, the data will be copied to the currently active CUDA device.
     * This region is only compiled when it's passed to nvcc.
     */

#ifdef __CUDACC__
    namespace cuda {
		
        // ------- 2D Tensor Voting -------

        template <typename Type>
        __global__ static void global_stickvote(Type* t_out, Type* L, Type* V, float sigma, int w, int shape0, int shape1) {

            int x0 = blockDim.x * blockIdx.x + threadIdx.x;                                       // get the x and y image coordinates for the current thread
            int x1 = blockDim.y * blockIdx.y + threadIdx.y;
            if (x0 >= shape0 || x1 >= shape1)                                                          // if not within bounds of image, return
                return;

            float a, b, c;
            stickvote_window(L, V, sigma, w, shape0, shape1, x0, x1, a, b, c);
            t_out[x0 * shape1 * 4 + x1 * 4 + 0] += a;
            t_out[x0 * shape1 * 4 + x1 * 4 + 1] += b;
            t_out[x0 * shape1 * 4 + x1 * 4 + 2] += b;
            t_out[x0 * shape1 * 4 + x1 * 4 + 3] += c;
        }

        template <typename Type>
        __global__ static void global_platevote(Type* t_out, Type* L, float sigma, int w, int s0, int s1, unsigned samples) {

            int x0 = blockDim.x * blockIdx.x + threadIdx.x;                                       // get the x and y image coordinates for the current thread
            int x1 = blockDim.y * blockIdx.y + threadIdx.y;
            if (x0 >= s0 || x1 >= s1)                                                          // if not within bounds of image, return
                return;

            float a, b, c;
            platevote_window(L, sigma, w, s0, s1, x0, x1, a, b, c, samples);

            t_out[x0 * s1 * 4 + x1 * 4 + 0] += a;
            t_out[x0 * s1 * 4 + x1 * 4 + 1] += b;
            t_out[x0 * s1 * 4 + x1 * 4 + 2] += b;
            t_out[x0 * s1 * 4 + x1 * 4 + 3] += c;
        }


        template <typename Type>
        static void tensorvote(Type* t_out, const Type* t_in, const Type sigma, const unsigned w,
            const unsigned shape0, const unsigned shape1,
            const bool stick = true, const bool plate = true, const unsigned samples = 10) {

            /*
            if (device < 0) {
                return cpu::tensorvote2((glm::mat2*)output_field,
                            (glm::vec2*)tira::cpu::evals2_symmetric<float>(input_field, shape0 * shape1),
                            (glm::vec2*)tira::cpu::evecs2polar_symmetric(input_field,
                            tira::cpu::evals2_symmetric<float>(input_field, shape0 * shape1), shape0 * shape1),
                            glm::vec2(sigma, sigma2), power, w, shape0, shape1, STICK, PLATE, samples);
            }
            */


            int device;
            HANDLE_ERROR(cudaGetDevice(&device));

            auto start = std::chrono::high_resolution_clock::now();
            cudaDeviceProp props;
            HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
            auto end = std::chrono::high_resolution_clock::now();
            float t_deviceprops = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

            int tensorFieldSize = 4 * shape0 * shape1;

            start = std::chrono::high_resolution_clock::now();
            float* L = tira::cuda::evals2_symmetric(t_in, shape0 * shape1, device);
            float* V = tira::cuda::evecs2polar_symmetric(t_in, L, shape0 * shape1, device);
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
            HANDLE_ERROR(cudaMalloc(&gpuV, shape0 * shape1 * 2 * sizeof(float)));
            HANDLE_ERROR(cudaMalloc(&gpuL, shape0 * shape1 * 2 * sizeof(float)));
            cudaDeviceSynchronize();
            end = std::chrono::high_resolution_clock::now();
            float t_devicealloc = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

            start = std::chrono::high_resolution_clock::now();
            // Copy input arrays
            HANDLE_ERROR(cudaMemcpy(gpuV, V, shape0 * shape1 * 2 * sizeof(float), cudaMemcpyHostToDevice));
            HANDLE_ERROR(cudaMemcpy(gpuL, L, shape0 * shape1 * 2 * sizeof(float), cudaMemcpyHostToDevice));
            cudaDeviceSynchronize();
            end = std::chrono::high_resolution_clock::now();

            float t_host2device = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

            // Specify the CUDA block and grid dimensions
            size_t blockDim = sqrt(props.maxThreadsPerBlock);
            dim3 threads(blockDim, blockDim);
            dim3 blocks(shape0 / threads.x + 1, shape1 / threads.y + 1);

            // perform tensor voting
            start = std::chrono::high_resolution_clock::now();
            if (stick)
                global_stickvote << < blocks, threads >> > (gpuOutputField, gpuL, gpuV, sigma, w, shape0, shape1);
            if (plate)
                global_platevote << < blocks, threads >> > (gpuOutputField, gpuL, sigma, w, shape0, shape1, samples);
            cudaDeviceSynchronize();
            end = std::chrono::high_resolution_clock::now();
            float t_voting = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();


            start = std::chrono::high_resolution_clock::now();
            // Copy the final result back from the GPU
            HANDLE_ERROR(cudaMemcpy(t_out, gpuOutputField, tensorFieldSize * sizeof(float), cudaMemcpyDeviceToHost));
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

            /*if (debug) {
                std::cout << "Eigendecomposition:  " << t_eigendecomposition << " ms" << std::endl;
                std::cout << "Voting: " << t_voting << " ms" << std::endl;
                std::cout << "cudaMemcpy (H->D):  " << t_host2device << " ms" << std::endl;
                std::cout << "cudaMemcpy (D->H):  " << t_device2host << " ms" << std::endl;
                std::cout << "cudaMalloc: " << t_devicealloc << " ms" << std::endl;
                std::cout << "cudaFree: " << t_devicefree << " ms" << std::endl;
                std::cout << "cudaDeviceProps: " << t_deviceprops << " ms" << std::endl;
            }*/
        }


		// ------- 3D Tensor Voting -------

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
                HANDLE_ERROR(cudaMemcpyToSymbol(d_neighbors_const, host.data(), host_bytes, 0, cudaMemcpyHostToDevice));
                dn.used_const = true;
                dn.d_ptr = nullptr;
            }
            else {
                HANDLE_ERROR(cudaMalloc(&dn.d_ptr, host_bytes));
                HANDLE_ERROR(cudaMemcpy(dn.d_ptr, host.data(), host_bytes, cudaMemcpyHostToDevice));
                dn.used_const = false;
            }
            return dn;
        }

        static void free_neighbors(DeviceNeighbors& dn) {
            if (!dn.used_const && dn.d_ptr) {
                HANDLE_ERROR(cudaFree(dn.d_ptr));
                dn.d_ptr = nullptr;
            }
            dn.count = 0;
            dn.used_const = false;
        }

        __device__ static inline double hyp2f1_neg_p_device(const unsigned int p, const double b, const double c, const double z) {
            double sum = 1.0;
            double term = 1.0;
            const double p_d = static_cast<double>(p);

            for (unsigned int k = 0; k < p; ++k) {
                double k_d = static_cast<double>(k);
                double ratio = ((-p_d + k_d) * (b + k_d)) / ((k_d + 1.0) * (c + k_d));
                term *= ratio * z;
                sum += term;
            }
            return sum;
        }

        __global__ static void spherical_to_cart3_kernel(float* Qout, const float* V6, bool evec_largest, size_t N) {
            size_t i = blockIdx.x * blockDim.x + threadIdx.x;
            if (i >= N) return;

            float theta, phi;
            // Stick vote uses the largest eigenvector (evec2)
            if (evec_largest) {
                theta = V6[i * 6 + 4];
                phi   = V6[i * 6 + 5];
            }
            // Plate vote uses the smallest eigenvector (evec0)
            else {
                theta = V6[i * 6 + 0];
                phi   = V6[i * 6 + 1];
            }

			glm::vec3 q = spherical_to_cartesian(theta, phi);

            Qout[i * 3 + 0] = q.x;
            Qout[i * 3 + 1] = q.y;
            Qout[i * 3 + 2] = q.z;
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
                const float scale = std::copysignf(fabsf(l2) - fabsf(l1), l2);

                // Eigenvector at voter in spherical angles
                const glm::vec3 q = Q[base_voter];
                
                stickvote3_accumulate_kernel(m00, m01, m02, m11, m12, m22,
					nb.d, nb.c1, nb.c2, q, power, scale);
            }

            // Write back to VT
            glm::mat3 Receiver;
            Receiver[0][0] = m00; Receiver[0][1] = m01; Receiver[0][2] = m02;
            Receiver[1][0] = m01; Receiver[1][1] = m11; Receiver[1][2] = m12;
            Receiver[2][0] = m02; Receiver[2][1] = m12; Receiver[2][2] = m22;
            VT[base_recv] += Receiver * norm;
        }

        __device__ static inline void platevote3_numerical_device(float& m00, float& m01, float& m02,
                                                                  float& m11, float& m12, float& m22,
                                                                  const float3 d, const float c1, const float c2,
			                                                      unsigned power, const glm::vec3 evec0, 
                                                                  float scale, const unsigned samples)
        {
            // A zero-vector has no orientation and cannot vote
            const float evec_len2 = evec0.x * evec0.x + evec0.y * evec0.y + evec0.z * evec0.z;
            if (samples == 0 or !(evec_len2 > TV_EPSILON)) return;

            float3 dn = d;

            glm::vec3 u, v;
            if (fabsf(evec0.z) < 0.999f) {
                u = glm::normalize(glm::cross(glm::vec3(0.0f, 0.0f, 1.0f), evec0));
            }
            else {
                u = glm::normalize(glm::cross(glm::vec3(1.0f, 0.0f, 0.0f), evec0));
            }
            v = glm::cross(evec0, u);

            const float dbeta = float(TV_PI) / float(samples);
            float v00 = 0.0f, v01 = 0.0f, v02 = 0.0f;
            float v11 = 0.0f, v12 = 0.0f, v22 = 0.0f;

            for (unsigned i = 0; i < samples; ++i) {
                float beta = dbeta * (float(i));
                float cb = cosf(beta);
                float sb = sinf(beta);
                glm::vec3 q = cb * u + sb * v;
                stickvote3_accumulate_kernel(v00, v01, v02, v11, v12, v22, dn, c1, c2, q, power, 1.0f);
            }

            m00 += scale * v00 / samples; m01 += scale * v01 / samples; m02 += scale * v02 / samples;
            m11 += scale * v11 / samples; m12 += scale * v12 / samples; m22 += scale * v22 / samples;
        }

        __device__ static inline void platevote3_analytical_device(float& m00, float& m01, float& m02,
            float& m11, float& m12, float& m22,
            const float3 d, const float c1, const float c2,
            const glm::vec3 evec0, float scale, const unsigned power,
            const double K0_d, const double K1_d)
        {
            // A zero-vector has no orientation and cannot vote
            const float evec_len2 = evec0.x * evec0.x + evec0.y * evec0.y + evec0.z * evec0.z;
            if (!(evec_len2 > TV_EPSILON)) return;
            float3 dn = d;
            glm::vec3 u, v;
            if (fabsf(evec0.z) < 0.999f) {
                u = glm::normalize(glm::cross(evec0, glm::vec3(0.0f, 0.0f, 1.0f)));
            }
            else {
                u = glm::normalize(glm::cross(evec0, glm::vec3(1.0f, 0.0f, 0.0f)));
            }
            v = glm::cross(evec0, u);

            // Rotation matrix Z = [u v evec0]
            // Zt = transpose(Z)
            const float Zt_00 = u.x, Zt_01 = u.y, Zt_02 = u.z;
            const float Zt_10 = v.x, Zt_11 = v.y, Zt_12 = v.z;
            const float Zt_20 = evec0.x, Zt_21 = evec0.y, Zt_22 = evec0.z;

            // Rotate d into local frame
            const float dx = dn.x * Zt_00 + dn.y * Zt_01 + dn.z * Zt_02;
            const float dy = dn.x * Zt_10 + dn.y * Zt_11 + dn.z * Zt_12;
            const float dz = dn.x * Zt_20 + dn.y * Zt_21 + dn.z * Zt_22;

            // Calculate local angles
            const float alpha = sqrtf(dx * dx + dy * dy);
            const float a2 = alpha * alpha;
            const double a2_d = fminf(static_cast<double>(a2), 1.0); // prevent issues with hyp2f1
            float phi = atan2f(dy, dx);

            // Compute beta and hypergeometric integrals
            const double p_d = static_cast<double>(power);
            const double J0 = 0.5 * TV_PI * hyp2f1_neg_p_device(power, 1.5, 2.0, a2_d);
            const double J1 = TV_PI * hyp2f1_neg_p_device(power, 0.5, 1.0, a2_d);

            // Terms in the rotated frame (A, B)
            const float tmp_a2 = 1.0f - 2.0f * a2;

            float A00 = (tmp_a2 * tmp_a2) * static_cast<float>(J0);
            float A02 = -2.0f * alpha * dz * tmp_a2 * static_cast<float>(J0);
            float A11 = static_cast<float>(J1 - J0);
            float A20 = A02;
            float A22 = 4.0f * a2 * (dz * dz) * static_cast<float>(J0);

            const float a2p = powf(a2, static_cast<float>(power));
            const float K0 = static_cast<float>(K0_d);
            const float K1 = static_cast<float>(K1_d);

            float B00 = a2p * (tmp_a2 * tmp_a2) * K0;
            float B02 = -2.0f * alpha * dz * tmp_a2 * a2p * K0;
            float B11 = a2p * (K1 - K0);
            float B20 = B02;
            float B22 = 4.0f * a2 * (dz * dz) * a2p * K0;

            // Rotate back around Z by phi
            float cph, sph;
            sincosf(phi, &sph, &cph);

            float ta00 = cph * cph * A00 + sph * sph * A11;
            float ta01 = cph * sph * (A00 - A11);
            float ta02 = cph * A02;
            float ta11 = sph * sph * A00 + cph * cph * A11;
            float ta12 = sph * A02;
            float ta22 = A22;

            float tb00 = cph * cph * B00 + sph * sph * B11;
            float tb01 = cph * sph * (B00 - B11);
            float tb02 = cph * B02;
            float tb11 = sph * sph * B00 + cph * cph * B11;
            float tb12 = sph * B02;
            float tb22 = B22;

            // Combine terms in local frame
            float V_loc_00 = c1 * ta00 + c2 * tb00;
            float V_loc_01 = c1 * ta01 + c2 * tb01;
            float V_loc_02 = c1 * ta02 + c2 * tb02;
            float V_loc_11 = c1 * ta11 + c2 * tb11;
            float V_loc_12 = c1 * ta12 + c2 * tb12;
            float V_loc_22 = c1 * ta22 + c2 * tb22;

            // Final rotation: V_global = Z * V_local * Zt
            float T00 = u.x * V_loc_00 + v.x * V_loc_01 + evec0.x * V_loc_02;
            float T01 = u.x * V_loc_01 + v.x * V_loc_11 + evec0.x * V_loc_12;
            float T02 = u.x * V_loc_02 + v.x * V_loc_12 + evec0.x * V_loc_22;

            float T10 = u.y * V_loc_00 + v.y * V_loc_01 + evec0.y * V_loc_02;
            float T11 = u.y * V_loc_01 + v.y * V_loc_11 + evec0.y * V_loc_12;
            float T12 = u.y * V_loc_02 + v.y * V_loc_12 + evec0.y * V_loc_22;

            float T20 = u.z * V_loc_00 + v.z * V_loc_01 + evec0.z * V_loc_02;
            float T21 = u.z * V_loc_01 + v.z * V_loc_11 + evec0.z * V_loc_12;
            float T22 = u.z * V_loc_02 + v.z * V_loc_12 + evec0.z * V_loc_22;

            // Accumulate (T * Zt) * scale
            m00 += scale * (T00 * u.x + T01 * v.x + T02 * evec0.x);
            m01 += scale * (T00 * u.y + T01 * v.y + T02 * evec0.y);
            m02 += scale * (T00 * u.z + T01 * v.z + T02 * evec0.z);
            m11 += scale * (T10 * u.y + T11 * v.y + T12 * evec0.y);
            m12 += scale * (T10 * u.z + T11 * v.z + T12 * evec0.z);
            m22 += scale * (T20 * u.z + T21 * v.z + T22 * evec0.z);

            return;
        }
        __global__ static void global_platevote3(glm::mat3* VT, const glm::vec3* L, const glm::vec3* Q_small, const Neighbor3D_CUDA* d_neigbors_glob,
            int nb_count, int usedConst, unsigned int power, float norm,
            int s0, int s1, int s2, unsigned samples, double K0_d, double K1_d) {
            int x2 = blockDim.x * blockIdx.x + threadIdx.x;
            int x1 = blockDim.y * blockIdx.y + threadIdx.y;
            int x0 = blockDim.z * blockIdx.z + threadIdx.z;
            if (x0 >= s0 || x1 >= s1 || x2 >= s2)
                return;

            const unsigned base_recv = (unsigned)x0 * s1 * s2 + (unsigned)x1 * s2 + (unsigned)x2;

            // Accumulator for the receiver voxel (symmetric 3x3 matrix)
            float m00 = 0.0f, m01 = 0.0f, m02 = 0.0f;
            float m11 = 0.0f, m12 = 0.0f, m22 = 0.0f;

            // Loop over neighbors
            for (int k = 0; k < nb_count; ++k) {
                const Neighbor3D_CUDA& nb = (usedConst ? d_neighbors_const[k] : d_neigbors_glob[k]);
                const int r0 = x0 + nb.dw;
                const int r1 = x1 + nb.dv;
                const int r2 = x2 + nb.du;
                if ((unsigned)r0 >= s0 || (unsigned)r1 >= s1 || (unsigned)r2 >= s2) continue;   // out of bounds
                const unsigned base_voter = (unsigned)r0 * s1 * s2 + (unsigned)r1 * s2 + (unsigned)r2;

                // Read eigenvalues at vote (only l0 and l1 for plate)
                const float l0 = L[base_voter].x;
                const float l1 = L[base_voter].y;
                float scale = std::copysignf(fabsf(l1) - fabsf(l0), l1);

                // Eigenvector at voter in spherical angles
                const glm::vec3 evec0 = Q_small[base_voter];

                // Numerical integration
                if (samples > 0)
                    platevote3_numerical_device(m00, m01, m02, m11, m12, m22,
                        nb.d, nb.c1, nb.c2, power, evec0, scale, samples);
                else
                    platevote3_analytical_device(m00, m01, m02, m11, m12, m22,
                        nb.d, nb.c1, nb.c2, evec0, scale, power, K0_d, K1_d);
            }
            glm::mat3 Receiver;
            Receiver[0][0] = m00; Receiver[0][1] = m01; Receiver[0][2] = m02;
            Receiver[1][0] = m01; Receiver[1][1] = m11; Receiver[1][2] = m12;
            Receiver[2][0] = m02; Receiver[2][1] = m12; Receiver[2][2] = m22;
            VT[base_recv] += Receiver * norm;
        }

        static void tensorvote3(const float* input_field, float* output_field, unsigned int s0, unsigned int s1, unsigned int s2, float* largest_q, float* smallest_q,
            float sigma, float sigma2, unsigned int w, unsigned int power, int device, bool STICK, bool PLATE, bool debug, unsigned samples) {
            
            HANDLE_ERROR(cudaSetDevice(device));

            auto start = std::chrono::high_resolution_clock::now();
            cudaDeviceProp props;
            HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
            auto end = std::chrono::high_resolution_clock::now();
            float t_deviceprops = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

            const size_t n_voxels = (size_t)s0 * s1 * s2;
            const size_t tensorField_bytes = sizeof(float) * 9 * n_voxels;
            const size_t evals_bytes = sizeof(float) * 3 * n_voxels;
            const size_t evecs_bytes = sizeof(float) * 6 * n_voxels;

            // Build neighbor table on host
            start = std::chrono::high_resolution_clock::now();
            const glm::vec2 sig(sigma, sigma2);
            std::vector<Neighbor3D> NB = cpu::build_neighbors3d((int)w, sig);
            DeviceNeighbors d_nb = upload_neighbors(NB);
            end = std::chrono::high_resolution_clock::now();
            float t_buildnb = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (debug)
                std::cout << "Neighbor table: " << NB.size() << " entries, build+upload time: " << t_buildnb << " ms" << std::endl;

            // Eigendecomposition on GPU
            start = std::chrono::high_resolution_clock::now();
            // float* L = tira::cuda::evals3_symmetric(input_field, n_voxels, device);
            // float* V = tira::cuda::evecs3spherical_symmetric(input_field, L, n_voxels, device);
            float* L = tira::cpu::evals3_symmetric(input_field, n_voxels);
            // float* V = tira::cpu::evecs3spherical_symmetric(input_field, L, n_voxels);
            end = std::chrono::high_resolution_clock::now();
            float t_eigendecomposition = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

            // Declare GPU arrays
            float* gpuOutputField;
            // float* gpuV;
            float* gpuL;

            // Check if input/eigens are already on device
            start = std::chrono::high_resolution_clock::now();
            cudaPointerAttributes attrL, attrV;
            HANDLE_ERROR(cudaPointerGetAttributes(&attrL, L));
            // HANDLE_ERROR(cudaPointerGetAttributes(&attrV, V));

            if (attrL.type == cudaMemoryTypeDevice) gpuL = L;
            else {
                HANDLE_ERROR(cudaMalloc(&gpuL, evals_bytes));
                HANDLE_ERROR(cudaMemcpy(gpuL, L, evals_bytes, cudaMemcpyHostToDevice));
            }
            /*if (attrV.type == cudaMemoryTypeDevice) gpuV = V;
            else {
                HANDLE_ERROR(cudaMalloc(&gpuV, evecs_bytes));
                HANDLE_ERROR(cudaMemcpy(gpuV, V, evecs_bytes, cudaMemcpyHostToDevice));
            }*/
            end = std::chrono::high_resolution_clock::now();
            float t_host2device = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

            // Convert spherical angles to cartesian Q
            float* dQ_large = nullptr;
            float* dQ_small = nullptr;
            const size_t cartesian_bytes = sizeof(float) * 3 * n_voxels;
            if (STICK) {
                HANDLE_ERROR(cudaMalloc(&dQ_large, cartesian_bytes));
                HANDLE_ERROR(cudaMemcpy(dQ_large, largest_q, cartesian_bytes, cudaMemcpyHostToDevice));
            }
            if (PLATE) {
                HANDLE_ERROR(cudaMalloc(&dQ_small, cartesian_bytes));
                HANDLE_ERROR(cudaMemcpy(dQ_small, smallest_q, cartesian_bytes, cudaMemcpyHostToDevice));
            }
            /*{
                int t = 256; int g = (int)((n_voxels + t - 1) / t);
                if (STICK)
                    spherical_to_cart3_kernel << <g, t >> > (dQ_large, (const float*)gpuV, true, n_voxels);
                if (PLATE)
                    spherical_to_cart3_kernel << <g, t >> > (dQ_small, (const float*)gpuV, false, n_voxels);
                HANDLE_ERROR(cudaDeviceSynchronize());
            }*/
           
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
            float sticknorm = 1.0f / sticknorm3(sigma, sigma2, power);
            float platenorm = 1.0f; // Not implemented yet

            // Calculate K0 and K1 for analytical plate voting
            const double p_d = static_cast<double>(power);
            const double K0 = boost::math::beta(0.5, p_d + 1.5);
            const double K1 = boost::math::beta(0.5, p_d + 0.5);

            start = std::chrono::high_resolution_clock::now();
            if (STICK)
                global_stickvote3 << <blocks, threads >> > ((glm::mat3*)gpuOutputField, (const glm::vec3*)gpuL, (const glm::vec3*)dQ_large,
                    d_nb.d_ptr, d_nb.count, d_nb.used_const ? 1 : 0, power, sticknorm, (int)s0, (int)s1, (int)s2);
            if (PLATE)
                global_platevote3 << <blocks, threads >> > ((glm::mat3*)gpuOutputField, (const glm::vec3*)gpuL, (const glm::vec3*)dQ_small,
                    d_nb.d_ptr, d_nb.count, d_nb.used_const ? 1 : 0, power, sticknorm, (int)s0, (int)s1, (int)s2, samples, K0, K1);
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
			if (attrL.type == cudaMemoryTypeDevice) HANDLE_ERROR(cudaFree(L));              // This is fragile but works for now
			else delete[] L;                                                                // works for either new[] or cudaMalloc
            // if (attrV.type == cudaMemoryTypeDevice) HANDLE_ERROR(cudaFree(V));
            // else delete[] V;
            if (gpuL && gpuL != L) HANDLE_ERROR(cudaFree(gpuL));
            // if (gpuV && gpuV != V) HANDLE_ERROR(cudaFree(gpuV));
			if (dQ_large) HANDLE_ERROR(cudaFree(dQ_large));
			if (dQ_small) HANDLE_ERROR(cudaFree(dQ_small));
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
    #endif

}
