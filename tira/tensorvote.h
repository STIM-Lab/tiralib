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

	template <typename T>
    CUDA_CALLABLE static T sticknorm3(const T sigma1, const T sigma2, const unsigned p) {
        T pi_term = TV_PI * sqrt(TV_PI) / 2.0;
		T sig1_3 = sigma1 * sigma1 * sigma1;
		T sig2_3 = sigma2 * sigma2 * sigma2;
        T num1 = std::pow(2, 2 * p + 1) * factorial(p) * factorial(p);
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

    static void tensorvote2_cpu(glm::mat2* VT, glm::vec2* L, glm::vec2* V, glm::vec2 sigma, unsigned int power, const unsigned w,
        const unsigned s0, const unsigned s1, const bool STICK = true, const bool PLATE = true, const unsigned samples = 0) {

        const float sticknorm = 1.0 / sticknorm2(sigma[0], sigma[1], power);
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


    struct Neighbor3D {
		int du, dv, dw;                 // index offsets along x2 (u), x1 (v), and x0 (w)
		glm::vec3 d;                    // normalized direction from sender to receiver (du, dv, dw)
		float l2;                       // squared length of offset: du^2 + dv^2 + dw^2
        float c1;                       // exp(-l2 / sigma1^2)
		float c2;                       // exp(-l2 / sigma2^2)
    };

    inline std::vector<Neighbor3D> build_neighbors3d(int w, glm::vec2 sigma) {
        std::vector<Neighbor3D> neighbors;
        neighbors.reserve(static_cast<size_t>(w) * w * w);

        const int hw = w / 2;
		const float invsig1 = (sigma.x > 0) ? 1.0f / (sigma.x * sigma.x) : 0.0f;
		const float invsig2 = (sigma.y > 0) ? 1.0f / (sigma.y * sigma.y) : 0.0f;

        for (int dw = -hw; dw < hw; dw++) {                         // For each pixel in the window
            for (int dv = -hw; dv < hw; dv++) {
                for (int du = -hw; du < hw; du++) {
                    Neighbor3D n;
					n.du = du; n.dv = dv; n.dw = dw;
					n.l2 = float(du * du + dv * dv + dw * dw);
                    if (n.l2 == 0.0f) n.d = glm::vec3(0.0f, 0.0f, 0.0f);
					else n.d = glm::vec3(float(du), float(dv), float(dw)) / sqrtf(n.l2);
                    
					n.c1 = sigma.x > 0 ? expf(-n.l2 * invsig1) : 0.0f;
                    n.c2 = sigma.y > 0 ? expf(-n.l2 * invsig2) : 0.0f;
					neighbors.push_back(n);
                }
            }
        }
        return neighbors;
	}

	inline float term_power(float t, unsigned power) {
        float r = t;
        for (unsigned i = 1; i < power; ++i) r *= t;
		return r;
    }
    CUDA_CALLABLE inline void stickvote3_accumulate_kernel(glm::mat3& M, const Neighbor3D& n, const glm::vec3 q, const unsigned power, const float scale) {
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

    /// <summary>
	/// Accumulate the stick vote for a receiver at position x from the precomputed neighbor offsets and eigenvalues/vectors.
    /// </summary>
    /// <param name="L">pointer to a volume of eigenvalues</param>
    /// <param name="Q">pointer to a volume of largest eigenvectors in cartesian coordinates</param>
	/// <param name="NB">precomputed neighbor offsets and direction vectors</param>
    /// <param name="power">refinement term</param>
    /// <param name="norm"></param>
    /// <param name="w">width of the vote region</param>
    /// <param name="s0">size of the L and V images along the first dimension</param>
    /// <param name="s1">size of the L and V images along the second dimension</param>
    /// <param name="x">position of the receiver</param>
    /// <returns></returns>
    CUDA_CALLABLE static glm::mat3 stickvote3(const glm::vec3* L, const glm::vec3* Q, const std::vector<Neighbor3D>& NB,
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
			const float scale = std::copysign(std::abs(l2) - std::abs(l1), l2);

            // Calculate the accumulated contribution of (du,dv,dw) to (x,y,z)
			stickvote3_accumulate_kernel(Votee, n, q, power, scale);
        }
        return Votee;
    }

    CUDA_CALLABLE  static glm::mat3 platevote3_analytic(const glm::vec3& d, float c1, float c2, unsigned power, double K0, double K1) {

        // Calculate the distance between voter and votee
        glm::vec3 dn = d;
		float len = std::sqrt(dn.x * dn.x + dn.y * dn.y + dn.z * dn.z);
		if (len != 0.0f) dn = dn / len;
		else dn = glm::vec3(0.0f, 0.0f, 0.0f);

        const float dx = dn.x, dy = dn.y, dz = dn.z;

        // Calculate the length and angle of the reflected direction matrix on XY plane
        const float alpha = std::sqrt(dx * dx + dy * dy);
        const float a2 = alpha * alpha;
        const float phi = std::atan2(dy, dx);

        // Build the rotation matrix around Z axis by +/- phi
        glm::mat3 Rz(0.0f);
		const float cph = std::cos(phi), sph = std::sin(phi);
        Rz[0][0] = cph;   Rz[0][1] = -sph;
        Rz[1][0] = sph;   Rz[1][1] = cph;
        Rz[2][2] = 1.0f;
        const glm::mat3 Rz_rev = glm::transpose(Rz);

        // Compute beta and hypergeometric integrals
        const double p_d = static_cast<double>(power);
        const double a2_d = static_cast<double>(a2);
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

        // Rotate back to original coordinates
        const glm::mat3 term_a = Rz * A * Rz_rev;
        const glm::mat3 term_b = Rz * B * Rz_rev;

        // Combine
        glm::mat3 PlateVote(0.0f);
		PlateVote[0][0] = c1 * term_a[0][0] + c2 * term_b[0][0];
		PlateVote[0][1] = c1 * term_a[0][1] + c2 * term_b[0][1];
		PlateVote[0][2] = c1 * term_a[0][2] + c2 * term_b[0][2];
		PlateVote[1][0] = c1 * term_a[1][0] + c2 * term_b[1][0];
		PlateVote[1][1] = c1 * term_a[1][1] + c2 * term_b[1][1];
		PlateVote[1][2] = c1 * term_a[1][2] + c2 * term_b[1][2];
		PlateVote[2][0] = c1 * term_a[2][0] + c2 * term_b[2][0];
		PlateVote[2][1] = c1 * term_a[2][1] + c2 * term_b[2][1];
		PlateVote[2][2] = c1 * term_a[2][2] + c2 * term_b[2][2];
        return PlateVote;
    }

    CUDA_CALLABLE inline glm::mat3 platevote3_numerical(const glm::vec3& d, float c1, float c2, unsigned power, unsigned samples = 20) {
        glm::mat3 V(0.0f);
        if (samples == 0) return V;

		glm::vec3 dn = d;
		float len = std::sqrt(dn.x * dn.x + dn.y * dn.y + dn.z * dn.z);
		if (len != 0.0f) dn /= len; else dn = glm::vec3(0.0f, 0.0f, 0.0f);

		// Build an orthonomal basis (u,v) spanning the plane perpendicular to d
        glm::vec3 u(1, 0, 0), v(0, 1, 0);
        if (d.x != 0.0f || d.y != 0.0f || d.z != 0.0f) {
			// Choose any vector not colinear to d
			glm::vec3 a = (std::fabs(dn.z) < 0.999f) ? glm::vec3(0.0f, 0.0f, 1.0f) : glm::vec3(1.0f, 0.0f, 0.0f);
			u = glm::normalize(glm::cross(a, dn));
			v = glm::cross(dn, u);
		}

		// Integrate over [0, pi] to avoid double counting q and -q
		// (the stick vote is symmetric along the stick axis)
        const float dbeta = float(TV_PI) / float(samples);
        for (unsigned i = 0; i < samples; ++i) {
            const float beta = dbeta * float(i);
			const glm::vec3 q = std::cos(beta) * u + std::sin(beta) * v;
			stickvote3_accumulate_kernel(V, Neighbor3D{ 0,0,0,dn,0.0f,c1,c2 }, q, power, 1.0f);
        }

		// Take the average
        const float invN = 1.0f / float(samples);
		V[0][0] *= invN; V[0][1] *= invN; V[0][2] *= invN;
        V[1][0] *= invN; V[1][1] *= invN; V[1][2] *= invN;
		V[2][0] *= invN; V[2][1] *= invN; V[2][2] *= invN;
        return V;
    }

    CUDA_CALLABLE static glm::mat3 platevote3(const glm::vec3* L, const std::vector<Neighbor3D>& NB, const unsigned power,
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
        
            const float l0 = L[base].x;
            const float l1 = L[base].y;
            float scale = std::copysign(std::abs(l1) - std::abs(l0), l1);
            
            const glm::vec3 d = n.d;
			const float c1 = n.c1, c2 = n.c2;
            glm::mat3 V;
            if (samples > 0)
                // Numerical integration of stick votes to form a plate vote
                V = platevote3_numerical(d, c1, c2, power, samples);
            else
				// Analytical closed form solution from direction d
                V = platevote3_analytic(d, c1, c2, power, K0, K1);
            
			Receiver[0][0] += scale * V[0][0]; Receiver[0][1] += scale * V[0][1]; Receiver[0][2] += scale * V[0][2];
			Receiver[1][0] += scale * V[1][0]; Receiver[1][1] += scale * V[1][1]; Receiver[1][2] += scale * V[1][2];
			Receiver[2][0] += scale * V[2][0]; Receiver[2][1] += scale * V[2][1]; Receiver[2][2] += scale * V[2][2];
        }
        return Receiver;
    }

    static void tensorvote3_cpu(glm::mat3* VT, const glm::vec3* L, const glm::vec3* Q, glm::vec2 sigma, unsigned int power, const unsigned w,
        const unsigned s0, const unsigned s1, const unsigned s2, const bool STICK = true, const bool PLATE = true, const unsigned samples = 20) {
		const float sticknorm = 1.0 / sticknorm3(sigma.x, sigma.y, power);
        const float platenorm = 1.0 / TV_PI;
        // Pre-compute the neighbor offsets and Gaussian factors once for (w, sigmas)
        const auto NB = build_neighbors3d((int)w, sigma);

		// O(N * w^3) complexity
        for (int x0 = 0; x0 < s0; x0++) {
            for (int x1 = 0; x1 < s1; x1++) {
                for (int x2 = 0; x2 < s2; x2++) {
                    glm::mat3 Vote(0.0f);
                    if (STICK)
                        Vote += stickvote3(L, Q, NB, power, s0, s1, s2, glm::ivec3(x0, x1, x2));
                    if (PLATE)
                        Vote += platevote3(L, NB, power, s0, s1, s2, glm::ivec3(x0, x1, x2), samples);
                    VT[x0 * s1 * s2 + x1 * s2 + x2] = sticknorm * Vote;
                }
            }
        }
    }
}