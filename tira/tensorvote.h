#pragma once

#include <tira/cuda/callable.h>
#include <glm/glm.hpp>

#define TV_PI 3.14159265358979323846

template <typename T>
CUDA_CALLABLE static T decay(T term, T length, T sigma, const unsigned power = 1) {
    T c = exp(-(length * length) / (sigma * sigma));

    T tp = term;
    for (unsigned int pi = 1; pi < power; pi++)
        tp *= term;

    return c * tp;
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

CUDA_CALLABLE  static glm::mat2 platevote2(glm::vec2 uv, glm::vec2 sigma) {

    //float length = sqrt(u * u + v * v);                     // calculate the distance between voter and votee
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

void hsa_tensorvote2(const float* input_field, float* output_field, unsigned int s0, unsigned int s1, float sigma, float sigma2,
    unsigned int w, unsigned int power, int device, bool STICK, bool PLATE, bool debug, unsigned samples);

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

namespace tira::cpu {
    static void tensorvote2(glm::mat2* VT, glm::vec2* L, glm::vec2* V, glm::vec2 sigma, unsigned int power,
        const unsigned w, const unsigned s0, const unsigned s1, const bool STICK = true, const bool PLATE = true, const unsigned samples = 0) {

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
}