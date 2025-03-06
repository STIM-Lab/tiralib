#pragma once

#include <tira/cuda/callable.h>

#ifdef __CUDACC__
#include <cuda/std/array>
#include <cuda/std/complex>
using namespace cuda::std;
#else
#include <array>
#include <complex>
using namespace std;
#endif

#include <iomanip>
#include <sstream>
#include <chrono>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>


namespace tira {

	template<typename T>
	class vec3 {


	public:
		array<T, 3> _ptr;


		CUDA_CALLABLE vec3() = default;
		CUDA_CALLABLE vec3(vec3 const&) = default;

		CUDA_CALLABLE vec3(T x, T y, T z) {
			_ptr[0] = x;
			_ptr[1] = y;
			_ptr[2] = z;
		}

		

		// copy assignment operator
		CUDA_CALLABLE vec3& operator=(const vec3& b) {
			_ptr = b._ptr;
			return *this;
		}

		// other assignment operators
		CUDA_CALLABLE vec3& operator+=(const vec3& b) {
			_ptr[0] = _ptr[0] + b._ptr[0];
			_ptr[1] = _ptr[1] + b._ptr[1];
			_ptr[2] = _ptr[2] + b._ptr[2];
			return *this;
		}

		// unary operators
		CUDA_CALLABLE vec3 operator-() {
			return vec3(-_ptr[0], -_ptr[1], -_ptr[2]);
		}

		// binary operators
		CUDA_CALLABLE vec3 operator+(vec3 b) {
			return vec3(_ptr[0] + b._ptr[0], _ptr[1] + b._ptr[1], _ptr[2] + b._ptr[2]);
		}		
		CUDA_CALLABLE vec3 operator-(vec3 b) {
			return vec3(_ptr[0] - b._ptr[0], _ptr[1] - b._ptr[1], _ptr[2] - b._ptr[2]);
		}
		CUDA_CALLABLE vec3 operator*(T b) {
			return vec3(_ptr[0] * b, _ptr[1] * b, _ptr[2] * b);
		}
		CUDA_CALLABLE vec3 operator/(T b) {
			return vec3(_ptr[0] / b, _ptr[1] / b, _ptr[2] / b);
		}


		CUDA_CALLABLE T dot(vec3 b) {
			return _ptr[0] * b._ptr[0] + _ptr[1] * b._ptr[1] + _ptr[2] * b._ptr[2];
		}
		CUDA_CALLABLE T norm() {
			return sqrt(_ptr[0] * _ptr[0] + _ptr[1] * _ptr[1] + _ptr[2] * _ptr[2]);
		}
		CUDA_CALLABLE vec3 cross(vec3 b) {
			return vec3(_ptr[1] * b._ptr[2] - _ptr[2] * b._ptr[1],
				_ptr[2] * b._ptr[0] - _ptr[0] * b._ptr[2],
				_ptr[0] * b._ptr[1] - _ptr[1] * b._ptr[0]);
		}
		CUDA_CALLABLE vec3 normalize(T& len) {
			len = norm();
			return vec3(_ptr[0] / len, _ptr[1] / len, _ptr[2] / len);
		}
		CUDA_CALLABLE vec3 normalize() {
			T n;
			return normalize(n);
		}
		CUDA_CALLABLE T& operator[](size_t i) {
			return _ptr[i];
		}
		CUDA_CALLABLE T get(size_t i) {
			T temp = _ptr[i];
			return temp;
		}

	};

	template<typename T>
	class cvec3 : public vec3< complex<typename T> >{

	public:
		CUDA_CALLABLE cvec3() = default;

		CUDA_CALLABLE cvec3(cvec3& b) {
			vec3< complex<T> >::_ptr = b._ptr;
		}

		CUDA_CALLABLE cvec3(complex<T> x, complex<T> y, complex<T> z) {
			vec3< complex<T> >::_ptr[0] = x;
			vec3< complex<T> >::_ptr[1] = y;
			vec3< complex<T> >::_ptr[2] = z;
		}

		// copy assignment operator
		CUDA_CALLABLE cvec3& operator=(const cvec3& b) {
			vec3< complex<T> >::_ptr = b._ptr;
			return *this;
		}

		CUDA_CALLABLE void operator=(const vec3< complex<T> > b) {
			vec3< complex<T> >::_ptr[0] = b._ptr[0];
			vec3< complex<T> >::_ptr[1] = b._ptr[1];
			vec3< complex<T> >::_ptr[2] = b._ptr[2];
		}

		CUDA_CALLABLE complex<T> dot(cvec3 b) {
			return vec3< complex<T> >::_ptr[0] * conj(b._ptr[0]) + vec3< complex<T> >::_ptr[1] * conj(b._ptr[1]) + vec3< complex<T> >::_ptr[2] * conj(b._ptr[2]);
		}
		CUDA_CALLABLE complex<T> dot(vec3<T> b) {
			return vec3< complex<T> >::_ptr[0] * b[0] + vec3< complex<T> >::_ptr[1] * b[1] + vec3< complex<T> >::_ptr[2] * b[2];
		}
		CUDA_CALLABLE T norm() {
			return sqrt((vec3< complex<T> >::_ptr[0] * conj(vec3< complex<T> >::_ptr[0]) + 
				vec3< complex<T> >::_ptr[1] * conj(vec3< complex<T> >::_ptr[1]) + 
				vec3< complex<T> >::_ptr[2] * conj(vec3< complex<T> >::_ptr[2])).real());
		}
	};

	template <typename T>
	class planewave {
		cvec3<T> _E0;			// E vector at the origin (p = 0)
		cvec3<T> _k;				// complex k vector (direction of propagation times wavenumber) in 3D


		/*static glm::tmat3x3<T> _Rot(glm::tvec3<T> from, glm::tvec3<T> to) {
			glm::tvec3<T> from_hat = glm::normalize(from);
			glm::tvec3<T> to_hat = glm::normalize(to);

			glm::tvec3<T> cp = glm::cross(from_hat, to_hat);
			T cp_length_sq = cp[0] + cp[1] + cp[2];
			//if the cross product is zero
			if (cp_length_sq == 0) {
				// calculate the dot product to figure out which degenerate matrix to return
				T dp = glm::dot(from_hat, to_hat);
				if (dp > 0) {		// if the dot product is positive, the to and from vectors are identical
					return glm::tmat3x3<T>(1.0);
				}
				if (dp < 0) {		// if the dot product is negative, then the to and from vectors are pointing in opposite directions
					glm::tmat3x3<T> rot(-1.0);
					rot[0][0] = 1.0;
					return rot;
				}
			}
			glm::tvec3<T> axis = glm::normalize(cp);

			T dot = glm::dot(from_hat, to_hat);							// calculate the dot product (used to differentiate between two angles)
			T arcsin = std::asin(glm::length(cp));
			T angle;
			if (dot >= 0)
				angle = arcsin;									// calculate the angle of rotation (inverse sin of the length of the cross product)
			else
				angle = M_PI - arcsin;
			//T angle = std::asin(glm::length(cp));

			glm::tquat<T> q = glm::angleAxis(angle, axis);		// generate a quaternion from the rotation properties

			return glm::mat3_cast(q);							// generate and return a rotation matrix
		}*/
		/// <summary>
		/// Generate a rotation matrix used to transform coordinates relative to k into global coordinates
		/// </summary>
		/// <returns></returns>
		/*glm::tmat3x3<T> _Rot() {
			glm::tvec3<T> z(0.0, 0.0, 1.0);						// z axis (local coordinates)
			glm::tvec3<T> k(
				_k[0].real(), 
				_k[1].real(), 
				_k[2].real());			// convert k to a vector in global coordinates

			return _Rot(z, k);
		}*/

		/// <summary>
		/// Generate a rotation matrix that transforms global coordinates to coordinates relative to k
		/// </summary>
		/// <returns></returns>
		/*glm::tmat3x3<T> _iRot() {
			glm::tvec3<T> z(0.0, 0.0, 1.0);						// z axis (local coordinates)
			glm::tvec3<T> k(
				_k[0].real(), 
				_k[1].real(), 
				_k[2].real());			// convert k to a vector in global coordinates

			glm::tvec3<T> k_hat = glm::normalize(k);			// convert k to a unit vector

			glm::tvec3<T> cp = glm::cross(z, k_hat);			// calculate the cross product of z and k_hat (giving us the axis of rotation and angle)
			T cp_len_sq = cp[0] * cp[0] + cp[1] * cp[1] + cp[2] * cp[2];
			if (cp_len_sq == 0) {
				T z_dot_khat = glm::dot(z, k_hat);
				if (z_dot_khat > 0) return glm::tmat3x3<T>(1.0);
				if (z_dot_khat < 0) {
					glm::tmat3x3<T> M = glm::tmat3x3<T>(-1.0);
					M[0][0] = 1.0;
					return M;
				}
			}
			
			glm::tvec3<T> axis = glm::normalize(cp);			// calculate the axis of rotation (normalized cross product)
			T dot = glm::dot(z, k_hat);							// calculate the dot product (used to differentiate between two angles)
			T arcsin = std::asin(glm::length(cp));
			T angle;
			if (dot >= 0)
				angle = arcsin;									// calculate the angle of rotation (inverse sin of the length of the cross product)
			else
				angle = M_PI - arcsin;
			
			glm::tquat<T> q = glm::angleAxis(angle, axis);		// generate a quaternion from the rotation properties
			glm::tquat<T> iq = glm::inverse(q);					// calculate the inverse of q

			return glm::mat3_cast(iq);							// generate and return a rotation matrix
		}*/

		/// Multiply a real matrix and a complex vector
		/*glm::tvec3<std::complex<T>> _cMul(glm::tmat3x3<T> M,
			std::complex<T> x, 
			std::complex<T> y, 
			std::complex<T> z) {

			glm::tvec3<T> r(x.real(), y.real(), z.real());				// create a vector from the real coordinates
			glm::tvec3<T> i(x.imag(), y.imag(), z.imag());				// create a vector from the imaginary coordinates

			glm::tvec3<T> R = M * r;											// calculate the real components of the E vector
			glm::tvec3<T> I = M * i;											// calculate the imaginary components of the E vector

			glm::tvec3<std::complex<T>> v(										// create the complex vector from real and imaginary parts
				std::complex<T>(R.x, I.x),
				std::complex<T>(R.y, I.y),
				std::complex<T>(R.z, I.z)
			);
			return v;															// return the complex 3D vector result
		}*/

		/*std::complex<T> _cDot(glm::tvec3<std::complex<T>> a, glm::tvec3<T> b) {
			std::complex<T> c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
			return c;
		}*/

	public:

		/// <summary>
		/// Generate a set of plane waves filling a solid angle (with an optional internal obscuration) using Monte-Carlo sampling (Mersenne Twister).
		/// </summary>
		/// <param name="alpha">External angle to fill (if simulating an objective use asin(NA))</param>
		/// <param name="k">Wavenumber (2pi/lambda)</param>
		/// <param name="Ex">X component of the E vector</param>
		/// <param name="Ey">Y component of the E vector</param>
		/// <param name="N">Number of Monte-Carlo samples</param>
		/// <param name="beta">Optional internal obscuration (reflective optics will generally use beta = 0</param>
		/// <returns>An std::vector of plane waves that fill the desired solid angle</returns>
		/*static std::vector< tira::planewave<T> > SolidAngleMC(T alpha, T kx, T ky, T kz,
		std::complex<T> Ex, std::complex<T> Ey, std::complex<T> Ez, size_t N,
		T beta = 0, glm::vec<3, T> clip = glm::vec<3, T>(0.0, 0.0, 0.0), unsigned long seed = 0) {
			if (beta > alpha) throw "ERROR: beta cannot be larger than alpha";
			T cos_alpha = std::cos(alpha);														// calculate the cosine of the focusing angle
			T cos_beta = std::cos(beta);

			glm::tvec3<T> Z(0.0, 0.0, 1.0);
			glm::tvec3<T> K(kx, ky, kz);
			glm::vec<3, T> D = glm::normalize(K);									// calculate the direction of the k vector
			glm::tmat3x3<T> R = _Rot(Z, D);
			T k = std::sqrt(kx * kx + ky * ky + kz * kz);

			std::mt19937_64 mt{ static_cast<unsigned int>(										// create a mersenne twister random number generator
				seed
				) };
			std::uniform_real_distribution<T> z_dist(cos_alpha, cos_beta);							// create two distributions for z = [0, cos(alpha)] and theta = [0, 2pi]
			std::uniform_real_distribution<T> theta_dist(0.0, 2 * M_PI);
			T theta, phi, x, y, z;
			std::vector< tira::planewave<T> > P;												// allocate space for an array of N plane waves
			tira::planewave<T> p;
			for (unsigned int pi = 0; pi < N; pi++) {


				z = z_dist(mt);																	// generate a random z and theta
				theta = theta_dist(mt);
				phi = std::acos(z);																// calculate phi (project onto the sphere)
				x = std::cos(theta) * std::sin(phi);
				y = std::sin(theta) * std::sin(phi);
				z = std::cos(phi);

				glm::tvec3<T> D_new = R * glm::normalize(glm::tvec3<T>(x, y, z));
				if(glm::length(clip) == 0 || glm::dot(clip, D_new) < 0){
					p._k[0] = D_new[0] * k;
					p._k[1] = D_new[1] * k;
					p._k[2] = D_new[2] * k;
					p._E0[0] = (1.0 / (T)N) * Ex;
					p._E0[1] = (1.0 / (T)N) * Ey;
					p._E0[2] = (1.0 / (T)N) * Ez;
					P.push_back(p);
				}
			}
			return P;
		}*/

		/*static std::vector< tira::planewave<T> > SolidAnglePolar(T alpha, T kx, T ky, T kz,
			std::complex<T> Ex, std::complex<T> Ey, std::complex<T> Ez,
			size_t thetaN = 100, size_t phiN = 50, T beta = 0, glm::vec<3, T> clip = glm::vec<3, T>(0.0, 0.0, 0.0)) {
			size_t N = thetaN * phiN;
			if (beta > alpha) throw "ERROR: beta cannot be larger than alpha";
			T cos_alpha = std::cos(alpha);											// calculate the cosine of the focusing angle
			T cos_beta = std::cos(beta);
			glm::tvec3<T> Z(0.0, 0.0, 1.0);
			glm::tvec3<T> K(kx, ky, kz);
			glm::vec<3, T> D = glm::normalize(K);									// calculate the direction of the k vector
			glm::tmat3x3<T> R = _Rot(Z, D);
			T k = std::sqrt(kx * kx + ky * ky + kz * kz);
			T theta, phi, x, y, z;
			size_t pi;																// 1D index into the array
			T dtheta = 2 * M_PI / (T)thetaN;										// calculate the theta and phi increments
			T dphi = (alpha - beta) / (T)phiN;

			std::vector< tira::planewave<T> > P;									// allocate an array to store the plane waves
			tira::planewave<T> p;
			for (size_t itheta = 0; itheta < thetaN; itheta++) {
				for (size_t iphi = 0; iphi < phiN; iphi++) {
					theta = dtheta * (T)itheta;
					phi = beta + dphi * (T)iphi; 

					x = std::cos(theta) * std::sin(phi);							// convert from spherical to cartesian coordinates
					y = std::sin(theta) * std::sin(phi);
					z = std::cos(phi);
					glm::tvec3<T> D_new = R * glm::normalize(glm::tvec3<T>(x, y, z));
					if (glm::length(clip) == 0 || glm::dot(clip, D_new) < 0) {
						pi = itheta * phiN + iphi;
						p._k[0] = D_new[0] * k;
						p._k[1] = D_new[1] * k;
						p._k[2] = D_new[2] * k;
						p._E0[0] = Ex * std::sin(phi) * dtheta * dphi;
						p._E0[1] = Ey * std::sin(phi) * dtheta * dphi;
						p._E0[2] = Ez * std::sin(phi) * dtheta * dphi;
						P.push_back(p);
					}
				}
			}
			return P;
		}*/

		// Create a plane wave with a given k and E vector. This expects both k and E to be orthogonal.
		CUDA_CALLABLE planewave<T>(complex<T> kx, complex<T> ky, complex<T> kz, complex<T> Ex, complex<T> Ey, complex<T> Ez) {
			_k[0] = kx;
			_k[1] = ky;
			_k[2] = kz;

			

			// save E while correcting for any out-of-plane components
			_E0[0] = Ex;
			_E0[1] = Ey;
			_E0[2] = Ez;
		}

		CUDA_CALLABLE planewave<T>(){
			_k[0] = 0;
			_k[1] = 0;
			_k[2] = 1;
			_E0[0] = 1;
			_E0[1] = 0;
			_E0[2] = 0;
		}

		// copy constructor
		CUDA_CALLABLE planewave<T>(const planewave<T>& t) {
			_E0 = t._E0;
			_k = t._k;
		}

		/// <summary>
		/// Returns the 3-dimensional E vector at the origin
		/// </summary>
		/// <returns></returns>
		CUDA_CALLABLE cvec3<T> E0() {
			return cvec3(_E0[0], _E0[1], _E0[2]);
		}

		CUDA_CALLABLE T I0() {
			return (_E0[0] * conj(_E0[0]) + _E0[1] * conj(_E0[1]) + _E0[2] * conj(_E0[2])).real();
		}

		/// Returns the 3-dimensional E vector at any point in space
		CUDA_CALLABLE cvec3<T> E(T x, T y, T z){
			complex<T> i(0, 1);
			complex<T> k_dot_r = _k[0] * x + _k[1] * y + _k[2] * z;
			complex<T> phase = exp(i * k_dot_r);
			cvec3<T> A = E0();
			cvec3<T> E(A[0] * phase, A[1] * phase, A[2] * phase);
			return E;
		}

		/// <summary>
		/// Return the k vector
		/// </summary>
		/// <returns></returns>
		CUDA_CALLABLE cvec3<T> k() {
			return cvec3(_k[0], _k[1], _k[2]);
		}

		/// <summary>
		/// Return a real k vector (assumes user knows the k vector is real)
		/// </summary>
		/// <returns></returns>
		CUDA_CALLABLE vec3<T> k_real() {
			return vec3(_k[0].real(), _k[1].real(), _k[2].real());
		}

		/// <summary>
		/// Returns the wavenumber
		/// </summary>
		/// <returns></returns>
		CUDA_CALLABLE T k_mag() {
			complex<T> sum = _k[0] * conj(_k[0]) + _k[1] * conj(_k[1]) + _k[2] * conj(_k[2]);
			return std::sqrt(sum.real());
		}

		/// <summary>
		/// Get the real direction of the k vector
		/// </summary>
		/// <returns></returns>
		CUDA_CALLABLE vec3<T> k_dir() {
			T kmag = k_mag();
			return k_real() / kmag;
		}

		/// <summary>
		/// Convert the plane wave to a string.
		/// </summary>
		/// <returns></returns>
		std::string str(int width = 25, int precision = 4) {
			std::stringstream ss;
			vec3<T> s = k_dir();
			ss << "|k|: " << k_mag()<<std::endl;
			ss << "s: "<<std::scientific<<std::setprecision(precision)<<std::setw(width)<< s[0] <<std::setw(width)<< s[1] <<std::setw(width)<< s[2]<<std::endl;
			cvec3<T> E = E0();
			ss << "E: "<<std::scientific<<std::setprecision(precision)<<std::setw(width)<< E[0] <<std::setw(width)<< E[1] <<std::setw(width)<< E[2]<<std::endl;
			return ss.str();
		}

		CUDA_CALLABLE tira::planewave<T> wind(T x, T y, T z) {
			complex<T> k_dot_r = _k[0] * x + _k[1] * y + _k[2] * z;
			complex<T> i(0.0, 1.0);
			complex<T> phase = std::exp(i * k_dot_r);
			return tira::planewave<T>(_k[0], _k[1], _k[2], _E0[0] * phase, _E0[1] * phase, _E0[2] * phase);
		}

		/// Calculate the reflected and transmitted plane waves at a single boundary
		CUDA_CALLABLE planewave reflect(complex<T> nr) {
			vec3<T> ki = k_real();					// k vector
			T k;
			vec3<T> si = ki.normalize(k);			// direction of the k vector
			cvec3<T> Ei = E0();						// incident electric field vector at the origin

			T cos_theta_i = si[2];										// calculate the sine and cosine of the angle of incidence
			T sin_theta_i = std::sqrt(1 - cos_theta_i * cos_theta_i);
			//T theta_i = asin(sin_theta_i);

			// calculate the basis vectors for the plane of incidence
			vec3<T> z(0, 0, 1);
			vec3<T> y;
			if (sin_theta_i == 0) {
				complex<T> rp = (1.0f - nr) / (1.0f + nr);		// compute the Fresnel coefficients (only 2 needed)
				complex<T> tp = 2.0f / (1.0f + nr);

				cvec3<T> Er;
				Er = Ei * rp;				// apply the Fresnel coefficients to set the amplitude of the scattered plane waves
				// generate the reflected and transmitted waves
				return planewave<T>(0, 0, -ki[2], Er[0], Er[1], Er[2]);
			}
			y = vec3<T>(si[0] / sin_theta_i, si[1] / sin_theta_i, 0);

			vec3<T> x = z.cross(y);

			vec3<T> kr(ki[0], ki[1], -ki[2]);		// calculate the reflected k vector

			// use Snell's law to calculate the transmitted angle sine and cosine
			complex<T> sin_theta_t = (1.0f / nr) * sin_theta_i;
			complex<T> cos_theta_t = sqrt(1.0f - sin_theta_t * sin_theta_t);

			// calculate the Fresnel coefficients
			complex<T> sin_tmi = sin_theta_t * cos_theta_i - cos_theta_t * sin_theta_i;
			complex<T> sin_2_tmi = sin_tmi * sin_tmi;
			complex<T> tan_tmi = sin_2_tmi / (complex<T>(1, 0) - sin_2_tmi);
			complex<T> sin_tpi = sin_theta_t * cos_theta_i + cos_theta_t * sin_theta_i;
			complex<T> sin_2_tpi = sin_tpi * sin_tpi;
			complex<T> tan_tpi = sin_2_tpi / (complex<T>(1, 0) - sin_2_tpi);

			complex<T> rs = sin_tmi / sin_tpi;
			complex<T> rp = tan_tmi / tan_tpi;

			// calculate the p component directions for each E vector
			vec3<T> Eip_dir = y * cos_theta_i - z * sin_theta_i;
			vec3<T> Erp_dir = y * cos_theta_i + z * sin_theta_i;


			// calculate the s and t components of each E vector
			complex<T> Ei_s = Ei.dot(x);
			complex<T> Ei_p = Ei.dot(Eip_dir);
			complex<T> Er_s = rs * Ei_s;
			complex<T> Er_p = rp * Ei_p;

			// calculate the E vector for each plane wave
			cvec3<T> Er;
			Er[0] = Erp_dir[0] * Er_p + x[0] * Er_s;
			Er[1] = Erp_dir[1] * Er_p + x[1] * Er_s;
			Er[2] = Erp_dir[2] * Er_p + x[2] * Er_s;

			// generate the reflected and transmitted waves
			return planewave<T>(kr[0], kr[1], kr[2], Er[0], Er[1], Er[2]);
		}

		CUDA_CALLABLE planewave refract(complex<T> nr) {
			vec3<T> ki = k_real();					// k vector
			T k;
			vec3<T> si = ki.normalize(k);			// direction of the k vector
			cvec3<T> Ei = E0();						// incident electric field vector at the origin

			T cos_theta_i = si[2];										// calculate the sine and cosine of the angle of incidence
			T sin_theta_i = std::sqrt(1 - cos_theta_i * cos_theta_i);

			// calculate the basis vectors for the plane of incidence
			vec3<T> z(0, 0, 1);
			vec3<T> y;
			if (sin_theta_i == 0) {
				complex<T> rp = (1.0f - nr) / (1.0f + nr);		// compute the Fresnel coefficients (only 2 needed)
				complex<T> tp = 2.0f / (1.0f + nr);

				cvec3<T> Er;
				Er = Ei * rp;				// apply the Fresnel coefficients to set the amplitude of the scattered plane waves
				cvec3<T> Et;
				Et = Ei * tp;
				// generate the reflected and transmitted waves
				return planewave<T>(0, 0, ki[2] * nr, Et[0], Et[1], Et[2]);
			}
			y = vec3<T>(si[0] / sin_theta_i, si[1] / sin_theta_i, 0);

			vec3<T> x = z.cross(y);

			vec3<T> kr(ki[0], ki[1], -ki[2]);		// calculate the reflected k vector

			// use Snell's law to calculate the transmitted angle sine and cosine
			complex<T> sin_theta_t = (1.0f / nr) * sin_theta_i;
			complex<T> cos_theta_t = sqrt(1.0f - sin_theta_t * sin_theta_t);

			cvec3<T> kt;										// calculate the transmitted k vector
			complex<T> kt_mag = k * nr;
			kt[0] = kt_mag * (y[0] * sin_theta_t + z[0] * cos_theta_t);
			kt[1] = kt_mag * (y[1] * sin_theta_t + z[1] * cos_theta_t);
			kt[2] = kt_mag * (y[2] * sin_theta_t + z[2] * cos_theta_t);

			// calculate the Fresnel coefficients
			complex<T> sin_tmi = sin_theta_t * cos_theta_i - cos_theta_t * sin_theta_i;
			complex<T> cos_tmi = sqrt(complex<T>(1, 0) - sin_tmi * sin_tmi);
			complex<T> sin_2_tmi = sin_tmi * sin_tmi;
			complex<T> tan_tmi = sin_2_tmi / (complex<T>(1, 0) - sin_2_tmi);
			complex<T> sin_tpi = sin_theta_t * cos_theta_i + cos_theta_t * sin_theta_i;
			//complex<T> cos_tpi = std::sqrt(complex<T>(1, 0) - sin_tpi * sin_tpi);
			complex<T> sin_2_tpi = sin_tpi * sin_tpi;
			complex<T> tan_tpi = sin_2_tpi / (complex<T>(1, 0) - sin_2_tpi);

			complex<T> rs = sin_tmi / sin_tpi;
			complex<T> rp = tan_tmi / tan_tpi;
			complex<T> ts = (2.0f * (sin_theta_t * cos_theta_i)) / sin_tpi;
			complex<T> tp = ((2.0f * sin_theta_t * cos_theta_i) / (sin_tpi * cos_tmi));

			// calculate the p component directions for each E vector
			vec3<T> Eip_dir = y * cos_theta_i - z * sin_theta_i;
			vec3<T> Erp_dir = y * cos_theta_i + z * sin_theta_i;
			cvec3<T> Etp_dir(0.0, 0.0, 0.0);
			Etp_dir[0] = y[0] * cos_theta_t - z[0] * sin_theta_t;
			Etp_dir[1] = y[1] * cos_theta_t - z[1] * sin_theta_t;
			Etp_dir[2] = y[2] * cos_theta_t - z[2] * sin_theta_t;

			// calculate the s and t components of each E vector
			complex<T> Ei_s = Ei.dot(x);
			complex<T> Ei_p = Ei.dot(Eip_dir);
			complex<T> Er_s = rs * Ei_s;
			complex<T> Er_p = rp * Ei_p;
			complex<T> Et_s = ts * Ei_s;
			complex<T> Et_p = tp * Ei_p;

			// calculate the E vector for the transmitted wave
			cvec3<T> Et;
			Et[0] = Etp_dir[0] * Et_p + x[0] * Et_s;
			Et[1] = Etp_dir[1] * Et_p + x[1] * Et_s;
			Et[2] = Etp_dir[2] * Et_p + x[2] * Et_s;

			// return the transmitted wave
			return planewave<T>(kt[0], kt[1], kt[2], Et[0], Et[1], Et[2]);
		}

		CUDA_CALLABLE void scatter(complex<T> nr, planewave<T>& r, planewave<T>& t) {
			vec3<T> ki = k_real();					// k vector
			T k;
			vec3<T> si = ki.normalize(k);			// direction of the k vector
			cvec3<T> Ei = E0();						// incident electric field vector at the origin

			T cos_theta_i = si[2];										// calculate the sine and cosine of the angle of incidence
			T sin_theta_i = std::sqrt(1 - cos_theta_i * cos_theta_i);

			// calculate the basis vectors for the plane of incidence
			vec3<T> z(0, 0, 1);						
			vec3<T> y;
			if (sin_theta_i == 0) {
				complex<T> rp = (1.0f - nr) / (1.0f + nr);		// compute the Fresnel coefficients (only 2 needed)
				complex<T> tp = 2.0f / (1.0f + nr);

				cvec3<T> Er;
				Er = Ei * rp;				// apply the Fresnel coefficients to set the amplitude of the scattered plane waves
				cvec3<T> Et;
				Et = Ei * tp;
				// generate the reflected and transmitted waves
				r = planewave<T>(0, 0, -ki[2], Er[0], Er[1], Er[2]);
				t = planewave<T>(0, 0, ki[2] * nr, Et[0], Et[1], Et[2]);
				return;
			}
			y = vec3<T>(si[0] / sin_theta_i, si[1] / sin_theta_i, 0);
			
			vec3<T> x = z.cross(y);

			vec3<T> kr(ki[0], ki[1], -ki[2]);		// calculate the reflected k vector

			// use Snell's law to calculate the transmitted angle sine and cosine
			complex<T> sin_theta_t = (1.0f / nr) * sin_theta_i;	
			complex<T> cos_theta_t = std::sqrt(1.0f - sin_theta_t * sin_theta_t);

			cvec3<T> kt;										// calculate the transmitted k vector
			complex<T> kt_mag = k * nr;
			kt[0] = kt_mag * (y[0] * sin_theta_t + z[0] * cos_theta_t);
			kt[1] = kt_mag * (y[1] * sin_theta_t + z[1] * cos_theta_t);
			kt[2] = kt_mag * (y[2] * sin_theta_t + z[2] * cos_theta_t);

			
			
			// allocate variables to store the Fresnel coefficients
			complex<T> rs, rp, ts, tp;

			// allocate E-vector components for all plane waves
			complex<T> Ei_s, Ei_p, Er_s, Er_p, Et_s, Et_p;

			// calculate the s and t components of the incident E vector
			vec3<T> Eip_dir = y * cos_theta_i - z * sin_theta_i;
			Ei_s = Ei.dot(x);
			Ei_p = Ei.dot(Eip_dir);

			// calculate trigonometric functions required for the s Fresnel coefficients
			complex<T> sin_tmi = sin_theta_t * cos_theta_i - cos_theta_t * sin_theta_i;
			complex<T> sin_tpi = sin_theta_t * cos_theta_i + cos_theta_t * sin_theta_i;

			// calculate the s components for the reflected and transmitted vector
			// (the above if statement handles division by zero)
			rs = sin_tmi / sin_tpi;
			Er_s = rs * Ei_s;
			ts = (2.0f * (sin_theta_t * cos_theta_i)) / sin_tpi;
			Et_s = ts * Ei_s;

			// calculate the trigonometric functions required for the p Fresnel coefficients
			//complex<T> sin_2_tmi = sin_tmi * sin_tmi;
			//complex<T> tan_tmi = sin_tmi / sqrt(complex<T>(1, 0) - sin_2_tmi);
			complex<T> tan_tmi = sin_theta_t * cos_theta_t - sin_theta_i * cos_theta_i;
			complex<T> sin_2_tpi = sin_tpi * sin_tpi;
			if (sin_2_tpi.real() == 1.0f && sin_2_tpi.imag() == 0) {
				Er_p = 0.0f;
			}
			else {
				//complex<T> tan_tpi = sin_2_tpi / (complex<T>(1, 0) - sin_2_tpi);
				complex<T> tan_tpi = sin_theta_t * cos_theta_t + sin_theta_i * cos_theta_i;
				rp = tan_tmi / tan_tpi;
				Er_p = rp * Ei_p;
			}

			complex<T> cos_tmi = std::sqrt(complex<T>(1, 0) - sin_tmi * sin_tmi);
			//tp = ((2.0f * sin_theta_t * cos_theta_i) / (sin_tpi * cos_tmi));
			tp = (2.0f * sin_theta_t * cos_theta_i) / (sin_theta_t * cos_theta_t + sin_theta_i * cos_theta_i);
			Et_p = tp * Ei_p;
			
			// calculate the p component directions for each E vector
			
			vec3<T> Erp_dir = y * cos_theta_i + z * sin_theta_i;
			cvec3<T> Etp_dir(0.0, 0.0, 0.0);
			Etp_dir[0] = y[0] * cos_theta_t - z[0] * sin_theta_t;
			Etp_dir[1] = y[1] * cos_theta_t - z[1] * sin_theta_t;
			Etp_dir[2] = y[2] * cos_theta_t - z[2] * sin_theta_t;

			// calculate the E vector for each plane wave
			cvec3<T> Er;
			Er[0] = Erp_dir[0] * Er_p + x[0] * Er_s;
			Er[1] = Erp_dir[1] * Er_p + x[1] * Er_s;
			Er[2] = Erp_dir[2] * Er_p + x[2] * Er_s;

			cvec3<T> Et;
			Et[0] = Etp_dir[0] * Et_p + x[0] * Et_s;
			Et[1] = Etp_dir[1] * Et_p + x[1] * Et_s;
			Et[2] = Etp_dir[2] * Et_p + x[2] * Et_s;

			// generate the reflected and transmitted waves
			r = planewave<T>(kr[0], kr[1], kr[2], Er[0], Er[1], Er[2]);
			t = planewave<T>(kt[0], kt[1], kt[2], Et[0], Et[1], Et[2]);
		}

		/// @param P is a plane representing the position and orientation of the surface
		/// @param r is the ratio n1/n0 (n0 is real)
		/// @param r is the reflected component of the plane wave
		/// @param t is the transmitted component of the plane wave
		CUDA_CALLABLE void scatter(vec3<T> plane_normal, vec3<T> plane_position, complex<T> nr, planewave<T>& r, planewave<T>& t) {

			// the k component for the input plane wave must be real (n0 is non-absorbing)
#ifndef __CUDACC__
			if(_k[0].imag() != 0.0 || _k[1].imag() != 0.0 || _k[2].imag() != 0.0)
				throw "Cannot scatter a plane wave with an imaginary k-vector";
#endif

			vec3<T> ki(_k[0].real(), _k[1].real(), _k[2].real());	// make the k vector real
			vec3<T> kr;											// the reflected wave can only have a real k vector
			cvec3<T> kt(0.0, 0.0, 0.0);			// the transmitted plane wave can have a complex k vector
			cvec3<T> Ei = E0();;				// All E vectors are complex
			cvec3<T> Er(0.0, 0.0, 0.0);
			cvec3<T> Et(0.0, 0.0, 0.0);

			vec3<T> n = plane_normal.normalize();				// make sure that the plane normal is normalized
			vec3<T> ki_dir = ki.normalize();						// calculate the direction of the k vector

			// swap the plane normal if necessary (it should be pointing in the opposite direction of ki)
			if(ki_dir.dot(n) > 0) n = -n;

			// use Snll's Law to calculate the incident angle
			T cos_theta_i = ki_dir.dot(-n);
			T sin_theta_i = std::sqrt(1 - cos_theta_i * cos_theta_i);
			T theta_i = std::acos(cos_theta_i);

			// define the basis vectors for the calculation (plane of incidence)
			vec3<T> z_hat = -n;
			vec3<T> side = n * ki_dir.dot(n);
			vec3<T> diff = ki_dir - side;
			T diff_len_sq = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];

			// handle the degenerate case where theta_i is 0 (the plane wave hits head-on)
			if(diff_len_sq == 0){
				complex<T> rp = (1.0f - nr) / (1.0f + nr);		// compute the Fresnel coefficients (only 2 needed)
				complex<T> tp = 2.0f / (1.0f + nr);

				kr = -ki;
				kt[0] = ki[0] * nr;			// the transmission vector is just scaled by the refractive index change
				kt[1] = ki[1] * nr;
				kt[2] = ki[2] * nr;

				Er = Ei * rp;				// apply the Fresnel coefficients to set the amplitude of the scattered plane waves
				Et = Ei * tp;

				// calculate the phase offset based on the plane position
				T phase_r = plane_position.dot(ki - kr);
				complex<T> phase_t = 
					plane_position[0] * (ki[0] - kt[0]) +
					plane_position[1] * (ki[1] - kt[1]) +
					plane_position[2] * (ki[2] - kt[2]);
			}
			else{
				T cos_theta_r = cos_theta_i;
				T sin_theta_r = sin_theta_i;
				T theta_r = theta_i;

				complex<T> sin_theta_t = (1.0f/nr) * std::sin(theta_i);	// compute the sin of theta_t using Snell's law
				complex<T> cos_theta_t = std::sqrt(1.0f - sin_theta_t * sin_theta_t);
				complex<T> theta_t = std::asin(sin_theta_t);

				// calculate the remaining basis vectors
				vec3<T> y_hat = (ki_dir - side).normalize();
				vec3<T> x_hat = y_hat.cross(z_hat);

				// calculate the k-vector magnitudes (wave frequency)
				T ki_mag = ki.norm();
				T kr_mag = ki_mag;
				complex<T> kt_mag = ki_mag * nr;

				// calculate the k vector directions
				vec3<T> kr_dir = y_hat * sin_theta_r - z_hat * cos_theta_r;
				cvec3<T> kt_dir(0.0, 0.0, 0.0);
				kt_dir[0] = y_hat[0] * sin_theta_t + z_hat[0] * cos_theta_t;
				kt_dir[1] = y_hat[1] * sin_theta_t + z_hat[1] * cos_theta_t;
				kt_dir[2] = y_hat[2] * sin_theta_t + z_hat[2] * cos_theta_t;

				// calculate the k vectors
				ki = ki_dir * ki_mag;
				kr = kr_dir * kr_mag;
				kt = kt_dir * kt_mag;

				// calculate the Fresnel coefficients
				complex<T> rs = std::sin(theta_t - theta_i) / std::sin(theta_t + theta_i);
				complex<T> rp = std::tan(theta_t - theta_i) / std::tan(theta_t + theta_i);
				complex<T> ts = (2.0f * (sin_theta_t * cos_theta_i)) / std::sin(theta_t + theta_i);
				complex<T> tp = ((2.0f * sin_theta_t * cos_theta_i) / (std::sin(theta_t + theta_i) * std::cos(theta_t - theta_i)));

				// calculate the p component directions for each E vector
				vec3<T> Eip_dir = y_hat * cos_theta_i - z_hat * sin_theta_i;
				vec3<T> Erp_dir = y_hat * cos_theta_r + z_hat * sin_theta_r;
				cvec3<T> Etp_dir(0.0, 0.0, 0.0);
				Etp_dir[0] = y_hat[0] * cos_theta_t - z_hat[0] * sin_theta_t;
				Etp_dir[1] = y_hat[1] * cos_theta_t - z_hat[1] * sin_theta_t;
				Etp_dir[2] = y_hat[2] * cos_theta_t - z_hat[2] * sin_theta_t;
				
				// calculate the s and t components of each E vector
				complex<T> Ei_s = Ei.dot(x_hat);
				complex<T> Ei_p = Ei.dot(Eip_dir);
				complex<T> Er_s = rs * Ei_s;
				complex<T> Er_p = rp * Ei_p;
				complex<T> Et_s = ts * Ei_s;
				complex<T> Et_p = tp * Ei_p;

				// calculate the E vector for each plane wave
				Er[0] = Erp_dir[0] * Er_p + x_hat[0] * Er_s;
				Er[1] = Erp_dir[1] * Er_p + x_hat[1] * Er_s;
				Er[2] = Erp_dir[2] * Er_p + x_hat[2] * Er_s;

				Et[0] = Etp_dir[0] * Et_p + x_hat[0] * Et_s;
				Et[1] = Etp_dir[1] * Et_p + x_hat[1] * Et_s;
				Et[2] = Etp_dir[2] * Et_p + x_hat[2] * Et_s;
			}
			// calculate the phase offset based on the plane positions
			T phase_r = plane_position.dot(ki - kr);
			complex<T> phase_t =
				plane_position[0] * (ki[0] - kt[0]) +
				plane_position[1] * (ki[1] - kt[1]) +
				plane_position[2] * (ki[2] - kt[2]);

			// apply the phase offset
			Er = Er * std::exp(complex<T>(0, 1) * phase_r);
			Et = Et * std::exp(complex<T>(0, 1) * phase_t);

			// generate the reflected and transmitted waves
			r = planewave<T>(kr[0], kr[1], kr[2], Er[0], Er[1], Er[2]);
			t = planewave<T>(kt[0], kt[1], kt[2], Et[0], Et[1], Et[2]);
		}


	};
}
