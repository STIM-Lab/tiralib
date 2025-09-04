#pragma once

#include <tira/cuda/callable.h>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <glm/glm.hpp>

#define TV_PI 3.14159265358979323846

namespace tira::tensorvote {
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
        T l_average = (l0 + l1 + l2) / 3.0;
        T numer = (l2 - l_average) * (l2 - l_average) + (l1 - l_average) * (l1 - l_average) + (l0 - l_average) * (l0 - l_average);
        T denom = l0 * l0 + l1 * l1 + l2 * l2;
        return std::sqrt(3.0 / 2.0 * (numer / denom));
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
        T ex = std::pow(2, 2 * p);
        T facp = factorial(p);
        T trig_int = num / (ex * facp * facp);
        return trig_int * (sigma1 * sigma1 + sigma2 * sigma2);
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

        float eta1 = 0;
        float eta2 = 0;
        if (sigma[0] > 0)
            eta1 = decay(1 - qTd * qTd, l, sigma[0], power);                       // calculate the decay function
        if (sigma[1] > 0)
            eta2 = decay(qTd * qTd, l, sigma[1], power);

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
        for (int v = -hw; v < hw; v++) {                    // for each pixel in the window
            const int r0 = static_cast<int>(x0) + v;
            if (r0 >= 0 && r0 < s0) {
                for (int u = -hw; u < hw; u++) {

                    const int r1 = static_cast<int>(x1) + u;
                    if (r1 >= 0 && r1 < s1) {
                        // calculate the contribution of (u,v) to (x,y)
                        glm::vec2 Vpolar = V[r0 * s1 + r1];
                        const float theta = Vpolar[1];
                        const glm::vec2 uv(u, v);
                        glm::mat2 vote = stickvote2(uv, sigma, theta, power);
                        const float l0 = L[r0 * s1 + r1][0];
                        const float l1 = L[r0 * s1 + r1][1];
                        float scale = std::abs(l1) - std::abs(l0);
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
            e1 = std::exp(-l2 / s12);
        float e2 = 0;
        if (sigma[1] > 0)
            e2 = std::exp(-l2 / s22);

        const float alpha = std::atan2(uv[1], uv[0]);
        const float two_a = 2 * alpha;
        const float cos_2a = std::cos(two_a);
        const float sin_2a = std::sin(two_a);
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

        for (int v = -hw; v < hw; v++) {                    // for each pixel in the window
            const int r0 = static_cast<int>(x0) + v;
            if (r0 >= 0 && r0 < s0) {
                for (int u = -hw; u < hw; u++) {
                    int r1 = static_cast<int>(x1) + u;
                    if (r1 >= 0 && r1 < s1) {
                        const float l0 = L[r0 * s1 + r1][0];
                        if (l0 != 0) {
                            const glm::vec2 uv(u, v);
                            if (samples > 0)             // if a sample number is provided, use numerical integration
                                Receiver += std::abs(l0) * platevote2_numerical(uv, sigma, samples);
                            else                         // otherwise use analytical integration (in progress)
                                Receiver += std::abs(l0) * platevote2(uv, sigma);
                        }
                    }
                }
            }
        }
        return Receiver;
    }

    /// <summary>
    /// Calculate the 3D stick vote for the relative position (u, v) given the voter eigenvales and eigenvectors
    /// </summary>
    /// <param name="uv">position of the receiver relative to the voter</param>
    /// <param name="sigma">decay value (standard deviation)</param>
    /// <param name="theta">orientation of the voter</param>
    /// <param name="power">refinement term</param>
    /// <returns></returns>
    CUDA_CALLABLE static glm::mat3 stickvote3(const glm::vec3 uvw, const glm::vec2 sigma, const float theta,
        const float phi, const unsigned power) {

        const float cos_theta = cosf(theta);
        const float sin_theta = sinf(theta);
        const float cos_phi = cosf(phi);
        const float sin_phi = sinf(phi);

        const glm::vec3 q(cos_theta * sin_phi, sin_theta * sin_phi, cos_phi);

        glm::vec3 d = uvw;                                                  // normalize the direction vector
        const float l = glm::length(d);                                     // calculate ell (distance between voter/votee)
        if (l == 0) d = glm::vec3(q.x, q.y, q.z);                           // assumes that the voter DOES contribute to itself
        else d = glm::normalize(d);

        const float qTd = glm::dot(q, d);

        float eta1 = 0;
        float eta2 = 0;
        if (sigma[0] > 0)
            eta1 = decay(1 - qTd * qTd, l, sigma[0], power);                       // calculate the decay function
        if (sigma[1] > 0)
            eta2 = decay(qTd * qTd, l, sigma[1], power);

        const glm::mat3 R = glm::mat3(1.0f) - 2.0f * glm::outerProduct(d, d);
        const glm::vec3 Rq = R * q;
        const glm::mat3 RqRq = glm::outerProduct(Rq, Rq);

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
    CUDA_CALLABLE static glm::mat3 stickvote3(const glm::vec3* L, const glm::vec2* V, const glm::vec2 sigma, const unsigned power, const float norm,
        const int w, const unsigned s0, const unsigned s1, const unsigned s2, const glm::ivec3 x) {

        const int x0 = x[0];
        const int x1 = x[1];
        const int x2 = x[2];

        glm::mat3 Votee(0.0f);

        const int hw = w / 2;
        for (int dw = -hw; dw <= hw; dw++) {                         // For each pixel in the window
            const int r0 = x0 + dw;
            if (r0 < 0 || r0 >= (int)s0) continue;

            for (int dv = -hw; dv <= hw; dv++) {
                const int r1 = x1 + dv;
                if (r1 < 0 || r1 >= (int)s1) continue;

                for (int du = -hw; du <= hw; du++) {
                    const int r2 = x2 + du;
                    if (r2 < 0 || r2 >= (int)s2) continue;

					// Flat index of the voter
                    const unsigned base = (unsigned)r0 * s1 * s2 + (unsigned)r1 * s2 + (unsigned)r2;

					// Largest eigenvector in polar coordinates
                    glm::vec2 Vpolar = V[base];
                    const float theta = Vpolar.x;
                    const float phi = Vpolar.y;

                    // Calculate the contribution of (du,dv,dw) to (x,y,z)
                    const glm::vec3 uvw((float)du, (float)dv, (float)dw);
                    glm::mat3 vote = stickvote3(uvw, sigma, theta, phi, power);

                    const float l1 = L[base][1];
                    const float l2 = L[base][2];
                    float scale = std::copysign(std::abs(l2) - std::abs(l1), l2);
                    Votee += scale * vote * norm;
                }
            }
        }
        return Votee;
    }

    CUDA_CALLABLE  static glm::mat3 platevote3(glm::vec3 uvw, glm::vec2 sigma, const unsigned power) {

        // calculate the distance between voter and votee
        const float length = sqrt(uvw[0] * uvw[0] + uvw[1] * uvw[1] + uvw[2] * uvw[2]);
        glm::vec3 d(0.0f);
        if (length != 0) d = glm::normalize(uvw);                 // normalize direction vector

        float dx = d.x;
        float dy = d.y;
        float dz = d.z;

        // calculate the length and angle of the reflected direction matrix on xy-plane
        float alpha = sqrt(dx * dx + dy * dy);
        float a2 = alpha * alpha;
        float phi = std::atan2(dy, dx);

        // build the rotation matrix about the z-axis by phi degrees
        glm::mat3 Rz(0.0f);
        Rz[0][0] = std::cos(phi);   Rz[0][1] = -std::sin(phi);
        Rz[1][0] = std::sin(phi);   Rz[1][1] = std::cos(phi);
        Rz[2][2] = 1.0f;

        // build the rotation matrix about the z-axis by -phi degrees
        glm::mat3 Rz_rev = glm::transpose(Rz);

        // compute beta and hypergeometric integrals
        double p_d = static_cast<double>(power);
        double a2_d = static_cast<double>(a2);
        double J0 = 0.5 * TV_PI * boost::math::hypergeometric_pFq(std::vector<double>{ -p_d, 1.5 },
            std::vector<double>{ 2.0 }, a2_d);
        double J1 = TV_PI * boost::math::hypergeometric_pFq(std::vector<double>{ -p_d, 0.5 },
            std::vector<double>{ 1.0 }, a2_d);
        double K0 = boost::math::beta(0.5, p_d + 1.5);
        double K1 = boost::math::beta(0.5, p_d + 0.5);

        // calcualte each term
        glm::mat3 A(0.0f);
        glm::mat3 B(0.0f);

        float tmp_a2 = 1.0f - 2.0f * a2;

        A[0][0] = (tmp_a2 * tmp_a2) * static_cast<float>(J0);
        A[0][2] = -2.0f * alpha * dz * tmp_a2 * static_cast<float>(J0);
        A[2][0] = -2.0f * alpha * dz * tmp_a2 * static_cast<float>(J0);
        A[1][1] = static_cast<float>(J1 - J0);
        A[2][2] = 4.0f * a2 * (dz * dz) * static_cast<float>(J0);

        float a2p = std::pow(a2, static_cast<float>(power));

        B[0][0] = a2p * (tmp_a2 * tmp_a2) * static_cast<float>(K0);
        B[0][2] = -2.0f * alpha * a2p * dz * tmp_a2 * static_cast<float>(K0);
        B[2][0] = -2.0f * alpha * a2p * dz * tmp_a2 * static_cast<float>(K0);
        B[1][1] = a2p * static_cast<float>(K1 - K0);
        B[2][2] = 4.0f * a2p * a2 * (dz * dz) * static_cast<float>(K0);

        // rotate back to original coordinates
        glm::mat3 term_a = Rz * A * Rz_rev;
        glm::mat3 term_b = Rz * B * Rz_rev;

        // calculate exponential terms
        float e1 = 0.0f, e2 = 0.0f;
        float l2 = length * length;
        float s12 = sigma[0] * sigma[0];
        float s22 = sigma[1] * sigma[1];
        if (sigma[0] > 0)
            e1 = std::exp(-l2 / s12);
        if (sigma[1] > 0)
            e2 = std::exp(-l2 / s22);

        // final integral
        glm::mat3 PlateVote = (e1 * term_a) + (e2 * term_b);

        return PlateVote;
    }

    CUDA_CALLABLE  static glm::mat3 platevote3_numerical(const glm::vec3 uvw, const glm::vec2 sigma, const unsigned power, const unsigned int n = 20) {

        const float dbeta = TV_PI / static_cast<double>(n);
        glm::mat3 V(0.0f);
        for (unsigned int i = 0; i < n; i++) {
            const float beta = dbeta * i;
            V = V + stickvote3(uvw, sigma, asin(1.0f), beta, 1);
        }
        const float norm = (float)1.0f / static_cast<float>(n);
        return V * norm;
    }

    CUDA_CALLABLE static glm::mat3 platevote3(const glm::vec3* L, const glm::vec2 sigma, const unsigned power,
        const int w, const unsigned s0, const unsigned s1, const unsigned s2, const glm::ivec3 x, const unsigned samples = 0) {

        const int x0 = x[0];
        const int x1 = x[1];
        const int x2 = x[2];

        glm::mat3 Receiver(0.0f);

        const int hw = w / 2;

        for (int w = -hw; w < hw; w++) {                         // for each pixel in the window
            const int r0 = static_cast<int>(x0) + w;
            if (r0 >= 0 && r0 < s0) {
                for (int v = -hw; v < hw; v++) {
                    const int r1 = static_cast<int>(x1) + v;
                    if (r1 >= 0 && r1 < s1) {
                        for (int u = -hw; u < hw; u++) {
                            const int r2 = static_cast<int>(x2) + u;
                            if (r2 >= 0 && r2 < s2) {
                                unsigned base = r0 * s1 * s2 + r1 * s2 + r2;
                                const float l1 = L[base][1];
                                const float l2 = L[base][2];
                                float scale = std::copysign(std::abs(l2) - std::abs(l1), l2);
                                if (scale != 0.0) {
                                    const glm::vec3 uvw(u, v, w);
                                    if (samples > 0)             // if a sample number is provided, use numerical integration
                                        Receiver += scale * platevote3_numerical(uvw, sigma, power, samples);
                                    else                         // otherwise use analytical integration (in progress)
                                        Receiver += scale * platevote3(uvw, sigma, power);
                                }
                            }
                        }
                    }
                }
            }
        }
        return Receiver;
    }

    static void tensorvote2_cpu(glm::mat2* VT, glm::vec2* L, glm::vec2* V, glm::vec2 sigma, unsigned int power, const unsigned w,
        const unsigned s0, const unsigned s1, const bool STICK = true, const bool PLATE = true, const unsigned samples = 0) {

        const float sticknorm = 1.0 / sticknorm2(sigma[0], sigma[1], power);
        for (int x0 = 0; x0 < s0; x0++) {
            for (int x1 = 0; x1 < s1; x1++) {
                glm::mat2 Vote(0.0f);
                if(STICK)
                    Vote = Vote + stickvote2(L, V, sigma, power, sticknorm, w, s0, s1, glm::ivec2(x0, x1));
                if(PLATE)
                    Vote = Vote + platevote2(L, sigma, w, s0, s1, glm::ivec2(x0, x1), samples);
                VT[x0 * s1 + x1] = Vote;
            }
        }
    }

    static void tensorvote3_cpu(glm::mat3* VT, const glm::vec3* L, const glm::vec2* V, glm::vec2 sigma, unsigned int power, const unsigned w,
        const unsigned s0, const unsigned s1, const unsigned s2, const bool STICK = true, const bool PLATE = true, const unsigned samples = 0) {
        const float sticknorm = 1.0;    // not yet implemented
        for (int x0 = 0; x0 < s0; x0++) {
            for (int x1 = 0; x1 < s1; x1++) {
                for (int x2 = 0; x2 < s2; x2++) {
                    glm::mat3 Vote(0.0f);
                    if (STICK)
                        Vote = Vote + stickvote3(L, V, sigma, power, sticknorm, w, s0, s1, s2, glm::ivec3(x0, x1, x2));
                    if (PLATE)
                        Vote = Vote + platevote3(L, sigma, power, w, s0, s1, s2, glm::ivec3(x0, x1, x2));
                    VT[x0 * s1 * s2 + x1 * s2 + x2] = Vote;
                }
            }
        }
    }
}