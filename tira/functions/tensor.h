#include <tira/cuda/callable.h>

#ifdef __CUDACC__
    #include <tira/cuda/error.h>
#endif

#include "eigen.h"
#include <glm/glm.hpp>

namespace tira {
    template <typename T>
    CUDA_CALLABLE static T FractionalAnisotropy(T l0, T l1, T l2) {
        T l_average = (l0 + l1 + l2) / T(3);
        T numer = (l2 - l_average) * (l2 - l_average) + (l1 - l_average) * (l1 - l_average) + (l0 - l_average) * (l0 - l_average);
        T denom = l0 * l0 + l1 * l1 + l2 * l2;
        return sqrtf(T(3) / T(2) * (numer / denom));
    }

    CUDA_CALLABLE static double factorial(const unsigned n) {
        double fac = 1;
        for (unsigned int i = 1; i <= n; i++)
            fac *= i;
        return fac;
    }

    CUDA_CALLABLE glm::vec3 spherical_to_cartesian(const float theta, const float phi) {
        const float st = sinf(theta), ct = cosf(theta);
        const float sp = sinf(phi),   cp = cosf(phi);
        return glm::vec3(ct * sp, st * sp, cp);
    }
}