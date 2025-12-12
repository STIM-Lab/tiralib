#pragma once

#include <tira/cuda/callable.h>

#ifdef __CUDACC__
    #include <tira/cuda/error.h>
#endif

#include "eigen.h"
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <glm/glm.hpp>

#define TV_PI 3.14159265358979323846
#define TIRA_VOTE_EPSILON 1e-12
#ifndef TV3_MAX_CONST_NB
#define TV3_MAX_CONST_NB 1536
#endif

namespace tira {

    template <typename T>
    CUDA_CALLABLE static T decay(T term, T length, T sigma, const unsigned power = 1) {
        T c = exp(-(length * length) / (sigma * sigma));

        T tp = term;
        for (unsigned int pi = 1; pi < power; pi++)
            tp *= term;

        return c * tp;
    }

    template <typename T>
    CUDA_CALLABLE static T FractionalAnisotropy(T l0, T l1, T l2) {
        T l_average = (l0 + l1 + l2) / T(3);
        T numer = (l2 - l_average) * (l2 - l_average) + (l1 - l_average) * (l1 - l_average) + (l0 - l_average) * (l0 - l_average);
        T denom = l0 * l0 + l1 * l1 + l2 * l2;
        return sqrtf(T(3) / T(2) * (numer / denom));
    }

    template <typename T>
    CUDA_CALLABLE static T PlateDecay2D(T length, T sigma) {
        T c = TV_PI * exp(-(length * length) / (sigma * sigma)) / 2.0f;
        return c;
    }

    CUDA_CALLABLE static double factorial(const unsigned n) {
        double fac = 1;
        for (unsigned int i = 1; i <= n; i++)
            fac *= i;
        return fac;
    }

    template <typename T>
    CUDA_CALLABLE static T sticknorm2(const T sigma1, const T sigma2, const unsigned p) {
        T num = TV_PI * factorial(2 * p);
        T ex = pow(2, 2 * p);
        T facp = factorial(p);
        T trig_int = num / (ex * facp * facp);
        return trig_int * (sigma1 * sigma1 + sigma2 * sigma2);
    }

    /**
        * @brief Normalization constant for 3D stick voting
        *
        * @tparam T        Floating point type
        * @param sigma1    Stick-vote sigma along first axis.
        * @param sigma2    Stick-vote sigma along second axis.
        * @param p         Stick-vote power (refinement exponent).
        * @return          Normalization factor for 3D stick votes.
    */
    template <typename T>
    CUDA_CALLABLE static T sticknorm3(const T sigma1, const T sigma2, const unsigned p) {
        T pi_term = TV_PI * sqrt(TV_PI) / 2.0;
        T sig1_3 = sigma1 * sigma1 * sigma1;
        T sig2_3 = sigma2 * sigma2 * sigma2;
        T num1 = pow(2, 2 * p + 1) * factorial(p) * factorial(p);
        T num2 = 2.0;
        T den1 = factorial(2 * p + 1);
        T den2 = 2 * p + 1;
        return pi_term * ((sig1_3 * num1 / den1) + (sig2_3 * num2 / den2));
    }

    /// <summary>
    /// Calculate the stick vote for the relative position (u, v) given the voter eigenvales and eigenvectors
    /// </summary>
    /// <param name="uv">position of the receiver relative to the voter</param>
    /// <param name="sigma">decay value (standard deviation)</param>
    /// <param name="theta">orientation of the voter</param>
    /// <param name="power">refinement term</param>
    /// <returns></returns>
    CUDA_CALLABLE static glm::mat2 stickvote2(const glm::vec2 uv, const glm::vec2 sigma, const float theta, const unsigned power) {

        const float cos_theta = cos(theta);
        const float sin_theta = sin(theta);
        const glm::vec2 q(cos_theta, sin_theta);

        glm::vec2 d = uv;                                       // normalize the direction vector
        const float l = glm::length(d);                               // calculate ell (distance between voter/votee)
        if (l == 0) d = glm::vec2(0, 0);                         // assumes that the voter DOES contribute to itself
        else d = glm::normalize(d);

        const float qTd = glm::dot(q, d);

        // Calculate the decay terms
        float eta1 = sigma[0] == 0.0f ? (l > TIRA_VOTE_EPSILON ? 0.0f : 1.0f) : decay(1 - qTd * qTd, l, sigma[0], power);
        float eta2 = sigma[1] == 0.0f ? (l > TIRA_VOTE_EPSILON ? 0.0f : 1.0f) : decay(qTd * qTd, l, sigma[1], power);

        const glm::mat2 R = glm::mat2(1.0f) - 2.0f * glm::outerProduct(d, d);
        const glm::vec2 Rq = R * q;
        const glm::mat2 RqRq = glm::outerProduct(Rq, Rq);

        return RqRq * (eta1 + eta2);
    }

    /// <summary>
    /// Calculate the receiver vote from an image of eigenvalues and eigenvectors.
    /// </summary>
    /// <param name="L">pointer to an image of eigenvalues</param>
    /// <param name="V">pointer to an image of eigenvectors</param>
    /// <param name="sigma">decay value (standard deviation)</param>
    /// <param name="power">refinement term</param>
    /// <param name="norm"></param>
    /// <param name="w">width of the vote region</param>
    /// <param name="s0">size of the L and V images along the first dimension</param>
    /// <param name="s1">size of the L and V images along the second dimension</param>
    /// <param name="x">position of the receiver</param>
    /// <returns></returns>
    CUDA_CALLABLE static glm::mat2 stickvote2(const glm::vec2* L, const glm::vec2* V, const glm::vec2 sigma, const unsigned power, const float norm,
        const int w, const int s0, const int s1, const glm::ivec2 x) {

        const int x0 = x[0];
        const int x1 = x[1];

        glm::mat2 Votee(0.0f);

        const int hw = w / 2;
        for (int v = -hw; v <= hw; v++) {                    // for each pixel in the window
            const int r0 = static_cast<int>(x0) + v;
            if (r0 >= 0 && r0 < s0) {
                for (int u = -hw; u <= hw; u++) {

                    const int r1 = static_cast<int>(x1) + u;
                    if (r1 >= 0 && r1 < s1) {
                        // calculate the contribution of (u,v) to (x,y)
                        glm::vec2 Vpolar = V[r0 * s1 + r1];
                        const float theta = Vpolar[1];
                        const glm::vec2 uv(u, v);
                        glm::mat2 vote = stickvote2(uv, sigma, theta, power);
                        const float l0 = L[r0 * s1 + r1][0];
                        const float l1 = L[r0 * s1 + r1][1];
                        float scale = fabs(l1) - fabs(l0);
                        if (l1 < 0) scale = scale * (-1);
                        Votee = Votee + scale * vote * norm;
                    }
                }
            }
        }
        return Votee;
    }

    // the plate tensor field is symmetric and not impacted by the refinement term
    // thus -> p = 1
    CUDA_CALLABLE  static glm::mat2 platevote2(glm::vec2 uv, glm::vec2 sigma) {

        // calculate the distance between voter and votee
        const float length = sqrt(uv[0] * uv[0] + uv[1] * uv[1]);
        const float l2 = length * length;
        const float s12 = sigma[0] * sigma[0];
        const float s22 = sigma[1] * sigma[1];
        float e1 = 0;
        if (sigma[0] > 0)
            e1 = expf(-l2 / s12);
        float e2 = 0;
        if (sigma[1] > 0)
            e2 = expf(-l2 / s22);

        const float alpha = atan2f(uv[1], uv[0]);
        const float two_a = 2 * alpha;
        const float cos_2a = cosf(two_a);
        const float sin_2a = sinf(two_a);
        glm::mat2 M;
        M[0][0] = cos_2a + 2;
        M[1][0] = sin_2a;
        M[0][1] = sin_2a;
        M[1][1] = 2 - cos_2a;

        const glm::mat2 I(1.0f);

        const float c = 1.0f / (TV_PI * (s12 + s22));

        return c * (e1 * (I - 0.25f * M) + e2 * (0.25f * M));
    }

    CUDA_CALLABLE  static glm::mat2 platevote2_numerical(const glm::vec2 uv, const glm::vec2 sigma, const unsigned int n = 20) {

        const float dtheta = TV_PI / static_cast<double>(n);
        glm::mat2 V(0.0f);
        for (unsigned int i = 0; i < n; i++) {
            const float theta = dtheta * i;
            V = V + stickvote2(uv, sigma, theta, 1);
        }
        const float norm = (float)1.0f / static_cast<float>(n);
        return V * norm;
    }

    CUDA_CALLABLE static glm::mat2 platevote2(const glm::vec2* L, const glm::vec2 sigma,
        const int w, const int s0, const int s1, const glm::ivec2 x, const unsigned samples = 0) {

        const int x0 = x[0];
        const int x1 = x[1];


        glm::mat2 Receiver(0.0f);

        const int hw = w / 2;

        for (int v = -hw; v <= hw; v++) {                    // for each pixel in the window
            const int r0 = static_cast<int>(x0) + v;
            if (r0 >= 0 && r0 < s0) {
                for (int u = -hw; u <= hw; u++) {
                    int r1 = static_cast<int>(x1) + u;
                    if (r1 >= 0 && r1 < s1) {
                        const float l0 = L[r0 * s1 + r1][0];
                        if (l0 != 0) {
                            const glm::vec2 uv(u, v);
                            if (samples > 0)             // if a sample number is provided, use numerical integration
                                Receiver += fabsf(l0) * platevote2_numerical(uv, sigma, samples);
                            else                         // otherwise use analytical integration (in progress)
                                Receiver += fabsf(l0) * platevote2(uv, sigma);
                        }
                    }
                }
            }
        }
        return Receiver;
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
     * CPU namespace contains functions that are run completely on the host. All input and output pointers
     * are allocated on the host.
     */
    namespace cpu {
        /**
         * @brief Raise a scalar to an integer power (power >= 1).
         *
         * @param t      Base scalar.
         * @param power  Exponent.
         * @return       t^power.
         */
        inline float term_power(float t, unsigned power) {
            float r = t;
            for (unsigned i = 1; i < power; ++i) r *= t;
            return r;
        }

        inline void tensorvote2(glm::mat2* VT, glm::vec2* L, glm::vec2* V, glm::vec2 sigma, unsigned int power, const unsigned w,
            const unsigned s0, const unsigned s1, const bool STICK = true, const bool PLATE = true, const unsigned samples = 0) {

            const float sticknorm = 1.0f / sticknorm2(sigma[0], sigma[1], power);
            for (int x0 = 0; x0 < s0; x0++) {
                for (int x1 = 0; x1 < s1; x1++) {
                    glm::mat2 Vote(0.0f);
                    if (STICK)
                        Vote = Vote + stickvote2(L, V, sigma, power, sticknorm, w, s0, s1, glm::ivec2(x0, x1));
                    if (PLATE)
                        Vote = Vote + platevote2(L, sigma, w, s0, s1, glm::ivec2(x0, x1), samples);
                    VT[x0 * s1 + x1] = Vote;
                }
            }
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
                        n.c1 = sigma.x == 0.0f ? (n.l2 > TIRA_VOTE_EPSILON ? 0.0f : 1.0f) : expf(-n.l2 * invsig1);
                        n.c2 = sigma.y == 0.0f ? (n.l2 > TIRA_VOTE_EPSILON ? 0.0f : 1.0f) : expf(-n.l2 * invsig2);
                        neighbors.push_back(n);
                    }
                }
            }
            return neighbors;
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
        inline void stickvote3_accumulate_kernel(glm::mat3& M, const Neighbor3D& n, const glm::vec3 q, const unsigned power, const float scale) {
            // Calculate the contribution of (du,dv,dw) to (x,y,z)
            const float qTd = q.x * n.d.x + q.y * n.d.y + q.z * n.d.z;
            const float qTd2 = qTd * qTd;
            float eta;
            if (power == 1)
                eta = n.c1 * (1 - qTd2) + n.c2 * qTd2;
            else
                eta = n.c1 * term_power(1 - qTd2, power) + n.c2 * term_power(qTd2, power);
            const float rx = q.x - 2.0f * qTd * n.d.x;
            const float ry = q.y - 2.0f * qTd * n.d.y;
            const float rz = q.z - 2.0f * qTd * n.d.z;
            const float term = scale * eta;
            M[0][0] += term * rx * rx; M[0][1] += term * rx * ry; M[0][2] += term * rx * rz;
            M[1][0] += term * ry * rx; M[1][1] += term * ry * ry; M[1][2] += term * ry * rz;
            M[2][0] += term * rz * rx; M[2][1] += term * rz * ry; M[2][2] += term * rz * rz;
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

            glm::mat3 Votee(0.0f);

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
                stickvote3_accumulate_kernel(Votee, n, q, power, scale);
            }
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
            if (samples == 0 or !(evec_len2 > TIRA_VOTE_EPSILON)) return V;

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
            for (unsigned i = 0; i < samples; ++i) {
                const float beta = dbeta * float(i);
                const glm::vec3 q = std::cos(beta) * u + std::sin(beta) * v;
                stickvote3_accumulate_kernel(V, Neighbor3D{ 0,0,0,dn,0.0f,c1,c2 }, q, power, 1.0f);
            }

            // Integral = sum / samples
            V[0][0] /= samples; V[0][1] /= samples; V[0][2] /= samples;
            V[1][0] /= samples; V[1][1] /= samples; V[1][2] /= samples;
            V[2][0] /= samples; V[2][1] /= samples; V[2][2] /= samples;
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
            const float platenorm = 1.0; // / TV_PI;           //  Not too sure about this
            // Pre-compute the neighbor offsets and Gaussian factors once for (w, sigmas)
            const auto NB = build_neighbors3d((int)w, sigma);

            // O(N * w^3) complexity
            for (int x0 = 0; x0 < s0; x0++) {
                for (int x1 = 0; x1 < s1; x1++) {
                    for (int x2 = 0; x2 < s2; x2++) {
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

        __global__ static void global_stickvote2(glm::mat2* VT, glm::vec2* L, glm::vec2* V, glm::vec2 sigma, unsigned int power, float norm,
            int w, int s0, int s1) {

            int x0 = blockDim.x * blockIdx.x + threadIdx.x;                                       // get the x and y image coordinates for the current thread
            int x1 = blockDim.y * blockIdx.y + threadIdx.y;
            if (x0 >= s0 || x1 >= s1)                                                          // if not within bounds of image, return
                return;

            glm::mat2 Receiver = stickvote2(L, V, sigma, power, norm, w, s0, s1, glm::ivec2(x0, x1));
            VT[x0 * s1 + x1] += Receiver;
        }

        __global__ static void global_platevote2(glm::mat2* VT, glm::vec2* L, glm::vec2 sigma, unsigned int power,
            int w, int s0, int s1, unsigned samples) {

            int x0 = blockDim.x * blockIdx.x + threadIdx.x;                                       // get the x and y image coordinates for the current thread
            int x1 = blockDim.y * blockIdx.y + threadIdx.y;
            if (x0 >= s0 || x1 >= s1)                                                          // if not within bounds of image, return
                return;

            glm::mat2 Receiver = platevote2(L, sigma, w, s0, s1, glm::ivec2(x0, x1), samples);
            VT[x0 * s1 + x1] += Receiver;
        }

        // TODO: assert both mats and evals are on the same side (host or device)
		// Free up L and V after use based on their allocation side (host/device)
        inline void tensorvote2(const float* input_field, float* output_field, unsigned int s0, unsigned int s1, float sigma, float sigma2,
                                unsigned int w, unsigned int power, int device, bool STICK, bool PLATE, bool debug, unsigned samples) {

            if (device < 0) return cpu::tensorvote2((glm::mat2*)output_field,
                                                        (glm::vec2*)tira::cpu::evals2_symmetric<float>(input_field, s0 * s1),
                                                        (glm::vec2*)tira::cpu::evecs2polar_symmetric(input_field,
                                                        tira::cpu::evals2_symmetric<float>(input_field, s0 * s1), s0 * s1),
				                                        glm::vec2(sigma, sigma2), power, w, s0, s1, STICK, PLATE, samples);

            HANDLE_ERROR(cudaSetDevice(device));

            auto start = std::chrono::high_resolution_clock::now();
            cudaDeviceProp props;
            HANDLE_ERROR(cudaGetDeviceProperties(&props, device));
            auto end = std::chrono::high_resolution_clock::now();
            float t_deviceprops = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

            int tensorFieldSize = 4 * s0 * s1;

            start = std::chrono::high_resolution_clock::now();
            float* L = tira::cuda::evals2_symmetric<float>(input_field, s0 * s1, device);
            float* V = tira::cuda::evecs2polar_symmetric(input_field, L, s0 * s1, device);
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
            HANDLE_ERROR(cudaMalloc(&gpuV, s0 * s1 * 2 * sizeof(float)));
            HANDLE_ERROR(cudaMalloc(&gpuL, s0 * s1 * 2 * sizeof(float)));
            cudaDeviceSynchronize();
            end = std::chrono::high_resolution_clock::now();
            float t_devicealloc = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

            start = std::chrono::high_resolution_clock::now();
            // Copy input arrays
            HANDLE_ERROR(cudaMemcpy(gpuV, V, s0 * s1 * 2 * sizeof(float), cudaMemcpyHostToDevice));
            HANDLE_ERROR(cudaMemcpy(gpuL, L, s0 * s1 * 2 * sizeof(float), cudaMemcpyHostToDevice));
            cudaDeviceSynchronize();
            end = std::chrono::high_resolution_clock::now();

            float t_host2device = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

            // Specify the CUDA block and grid dimensions
            size_t blockDim = sqrt(props.maxThreadsPerBlock);
            dim3 threads(blockDim, blockDim);
            dim3 blocks(s0 / threads.x + 1, s1 / threads.y + 1);

            float sn = 1.0 / sticknorm2(sigma, sigma2, power);

            if (debug)
                std::cout << "Stick Area: " << sn << std::endl;

            start = std::chrono::high_resolution_clock::now();
            if (STICK)
                global_stickvote2 << < blocks, threads >> > ((glm::mat2*)gpuOutputField, (glm::vec2*)gpuL, (glm::vec2*)gpuV, glm::vec2(sigma, sigma2), power, sn, w, s0, s1);
            if (PLATE)
                global_platevote2 << < blocks, threads >> > ((glm::mat2*)gpuOutputField, (glm::vec2*)gpuL, glm::vec2(sigma, sigma2), power, w, s0, s1, samples);
            cudaDeviceSynchronize();
            end = std::chrono::high_resolution_clock::now();
            float t_voting = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();


            start = std::chrono::high_resolution_clock::now();
            // Copy the final result back from the GPU
            HANDLE_ERROR(cudaMemcpy(output_field, gpuOutputField, tensorFieldSize * sizeof(float), cudaMemcpyDeviceToHost));
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

            if (debug) {
                std::cout << "Eigendecomposition:  " << t_eigendecomposition << " ms" << std::endl;
                std::cout << "Voting: " << t_voting << " ms" << std::endl;
                std::cout << "cudaMemcpy (H->D):  " << t_host2device << " ms" << std::endl;
                std::cout << "cudaMemcpy (D->H):  " << t_device2host << " ms" << std::endl;
                std::cout << "cudaMalloc: " << t_devicealloc << " ms" << std::endl;
                std::cout << "cudaFree: " << t_devicefree << " ms" << std::endl;
                std::cout << "cudaDeviceProps: " << t_deviceprops << " ms" << std::endl;
            }
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

        __device__ static inline float term_power_device(float t, unsigned power) {
            float r = t;
            for (unsigned i = 1; i < power; ++i) r *= t;
            return r;
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

			const float st = sinf(theta), ct = cosf(theta);
			const float sp = sinf(phi), cp = cosf(phi);

            Qout[i * 3 + 0] = ct * sp;
            Qout[i * 3 + 1] = st * sp;
            Qout[i * 3 + 2] = cp;
        }

        __device__ static inline void stickvote3_accumulate_kernel_device(float& m00, float& m01, float& m02,
                                                                          float& m11, float& m12, float& m22,
                                                                          const float3 d, const float c1, const float c2,
                                                                          const glm::vec3 q, const unsigned power, const float scale)
        {
            const float qTd = q.x * d.x + q.y * d.y + q.z * d.z;
            const float qTd2 = qTd * qTd;
            float eta;
            if (power == 1)
                eta = c1 * (1.0f - qTd2) + c2 * qTd2;
            else
                eta = c1 * term_power_device(1.0f - qTd2, power) + c2 * term_power_device(qTd2, power);
            const float rx = q.x - 2.0f * qTd * d.x;
            const float ry = q.y - 2.0f * qTd * d.y;
            const float rz = q.z - 2.0f * qTd * d.z;
            const float term = scale * eta;

            m00 += term * rx * rx; m01 += term * rx * ry; m02 += term * rx * rz;
            m11 += term * ry * ry; m12 += term * ry * rz; m22 += term * rz * rz;
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
                
                stickvote3_accumulate_kernel_device(m00, m01, m02, m11, m12, m22,
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
            if (samples == 0 or !(evec_len2 > TIRA_VOTE_EPSILON)) return;

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
                stickvote3_accumulate_kernel_device(v00, v01, v02, v11, v12, v22, dn, c1, c2, q, power, 1.0f);
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
            if (!(evec_len2 > TIRA_VOTE_EPSILON)) return;
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

        static void tensorvote3(const float* input_field, float* output_field, unsigned int s0, unsigned int s1, unsigned int s2, float sigma,
            float sigma2, unsigned int w, unsigned int power, int device, bool STICK, bool PLATE, bool debug, unsigned samples) {
            
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
            float* L = tira::cuda::evals3_symmetric(input_field, n_voxels, device);
            float* V = tira::cuda::evecs3spherical_symmetric(input_field, L, n_voxels, device);
            end = std::chrono::high_resolution_clock::now();
            float t_eigendecomposition = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

            // Declare GPU arrays
            float* gpuOutputField;
            float* gpuV;
            float* gpuL;

            // Check if input/eigens are already on device
            start = std::chrono::high_resolution_clock::now();
            cudaPointerAttributes attrL, attrV;
            HANDLE_ERROR(cudaPointerGetAttributes(&attrL, L));
            HANDLE_ERROR(cudaPointerGetAttributes(&attrV, V));

            if (attrL.type == cudaMemoryTypeDevice) gpuL = L;
            else {
                HANDLE_ERROR(cudaMalloc(&gpuL, evals_bytes));
                HANDLE_ERROR(cudaMemcpy(gpuL, L, evals_bytes, cudaMemcpyHostToDevice));
            }
            if (attrV.type == cudaMemoryTypeDevice) gpuV = V;
            else {
                HANDLE_ERROR(cudaMalloc(&gpuV, evecs_bytes));
                HANDLE_ERROR(cudaMemcpy(gpuV, V, evecs_bytes, cudaMemcpyHostToDevice));
            }
            end = std::chrono::high_resolution_clock::now();
            float t_host2device = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

            // Convert spherical angles to cartesian Q
            float* dQ_large = nullptr;
            float* dQ_small = nullptr;
            const size_t cartesian_bytes = sizeof(float) * 3 * n_voxels;
            if (STICK) HANDLE_ERROR(cudaMalloc(&dQ_large, cartesian_bytes));
            if (PLATE) HANDLE_ERROR(cudaMalloc(&dQ_small, cartesian_bytes));
            {
                int t = 256; int g = (int)((n_voxels + t - 1) / t);
                if (STICK)
                    spherical_to_cart3_kernel << <g, t >> > (dQ_large, (const float*)gpuV, true, n_voxels);
                if (PLATE)
                    spherical_to_cart3_kernel << <g, t >> > (dQ_small, (const float*)gpuV, false, n_voxels);
                HANDLE_ERROR(cudaDeviceSynchronize());
            }

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
            if (attrV.type == cudaMemoryTypeDevice) HANDLE_ERROR(cudaFree(V));
            else delete[] V;
            if (gpuL && gpuL != L) HANDLE_ERROR(cudaFree(gpuL));
            if (gpuV && gpuV != V) HANDLE_ERROR(cudaFree(gpuV));
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