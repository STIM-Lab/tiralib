#pragma once

#include "glm/glm.hpp"
#include "glm/gtc/quaternion.hpp"

#include <complex>
#include <vector>
#include <iomanip>
#include <chrono>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>


namespace tira {

	template <typename T>
	class planewave {
		std::complex<T> _E0[3];			// E vector at the origin (p = 0)
		std::complex<T> _k[3];			// complex k vector (direction of propagation times wavenumber) in 3D


		static glm::tmat3x3<T> _Rot(glm::tvec3<T> from, glm::tvec3<T> to) {
			glm::tvec3<T> from_hat = glm::normalize(from);
			glm::tvec3<T> to_hat = glm::normalize(to);

			/*if (from_hat[0] == to_hat[0] && from_hat[1] == to_hat[1] && from_hat[2] == to_hat[2])
				return glm::tmat3x3<T>(1.0);
			if (from_hat[0] == -to_hat[0] && from_hat[1] == -to_hat[1] && from_hat[2] == -to_hat[2]) {
				glm::tmat3x3<T> rot(-1.0);
				rot[0][0] = 1.0;
				return rot;
			}*/

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
		}
		/// <summary>
		/// Generate a rotation matrix used to transform coordinates relative to k into global coordinates
		/// </summary>
		/// <returns></returns>
		glm::tmat3x3<T> _Rot() {
			glm::tvec3<T> z(0.0, 0.0, 1.0);						// z axis (local coordinates)
			glm::tvec3<T> k(
				_k[0].real(), 
				_k[1].real(), 
				_k[2].real());			// convert k to a vector in global coordinates

			return _Rot(z, k);
		}

		/// <summary>
		/// Generate a rotation matrix that transforms global coordinates to coordinates relative to k
		/// </summary>
		/// <returns></returns>
		glm::tmat3x3<T> _iRot() {
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
		}

		/// Multiply a real matrix and a complex vector
		glm::tvec3<std::complex<T>> _cMul(glm::tmat3x3<T> M, 
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
		}

		std::complex<T> _cDot(glm::tvec3<std::complex<T>> a, glm::tvec3<T> b){
			std::complex<T> c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
			return c;
		}

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
		static std::vector< tira::planewave<T> > SolidAngleMC(T alpha, T kx, T ky, T kz, 
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
		}

		static std::vector< tira::planewave<T> > SolidAnglePolar(T alpha, T kx, T ky, T kz, 
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
		}

		// Create a plane wave with a given k vector and the corresponding perpendicular components of E.
		/*planewave<T>(std::complex<T> kx, std::complex<T> ky, std::complex<T> kz, std::complex<T> Ex_perp, std::complex<T> Ey_perp) {
			_k[0] = kx;
			_k[1] = ky;
			_k[2] = kz;
			_E0[0] = Ex_perp;
			_E0[1] = Ey_perp;
		}*/

		// Create a plane wave with a given k and E vector. This expects both k and E to be orthogonal.
		planewave<T>(std::complex<T> kx, std::complex<T> ky, std::complex<T> kz, std::complex<T> Ex, std::complex<T> Ey, std::complex<T> Ez) {
			_k[0] = kx;
			_k[1] = ky;
			_k[2] = kz;

			// ORTHOGONALIZE
			// calculate the length of the k vector
			/*std::complex<T> k_mag = sqrt(kx * kx + ky * ky + kz * kz);
			
			// calculate the normalized direction of k
			std::complex<T> dx = kx / k_mag;
			std::complex<T> dy = ky / k_mag;
			std::complex<T> dz = kz / k_mag;

			// calculate the projection of E0 onto k
			std::complex<T> E_dot_k = Ex * dx + Ey * dy + Ez * dz;
			*/

			// save E while correcting for any out-of-plane components
			_E0[0] = Ex;// -dx * E_dot_k;
			_E0[1] = Ey;// -dy * E_dot_k;
			_E0[2] = Ez;// -dz * E_dot_k;
		}

		planewave<T>(){
			_k[0] = 0;
			_k[1] = 0;
			_k[2] = 1;
			_E0[0] = 1;
			_E0[1] = 0;
			_E0[2] = 0;
		}

		/// <summary>
		/// Returns the 3-dimensional E vector at the origin
		/// </summary>
		/// <returns></returns>
		glm::vec<3, std::complex<T>> getE0() {
			//glm::tmat3x3<T> M = _Rot();							// get the matrix mapping local coordinates into global coordinates
			//return _cMul(M, _E0[0], _E0[1], 0.0);				// multiply the local E0 coordinates by the rotation matrix
			return glm::vec<3, std::complex<T> >(_E0[0], _E0[1], _E0[2]);
		}

		/// Returns the 3-dimensional E vector at any point in space
		glm::vec<3, std::complex<T>> getE(T x, T y, T z){
			std::complex<T> i(0, 1);
			std::complex<T> k_dot_r = _k[0] * x + _k[1] * y + _k[2] * z;
			std::complex<T> phase = std::exp(i * k_dot_r);
			glm::vec<3, std::complex<T>> E0 = getE0();
			glm::vec<3, std::complex<T>> E(E0[0] * phase, E0[1] * phase, E0[2] * phase);
			return E;
		}

		/// <summary>
		/// Return the k vector
		/// </summary>
		/// <returns></returns>
		glm::vec<3, std::complex<T>> getK() {
			return glm::vec<3, std::complex<T>>(_k[0], _k[1], _k[2]);
		}

		/// <summary>
		/// Return a real k vector (assumes user knows the k vector is real)
		/// </summary>
		/// <returns></returns>
		glm::vec<3, T> getKreal() {
			return glm::vec<3, T>(_k[0].real(), _k[1].real(), _k[2].real());
		}

		/// <summary>
		/// Returns the wavenumber
		/// </summary>
		/// <returns></returns>
		T getKmag() {
			return std::sqrt(_k[0].real() * _k[0].real() + _k[1].real() * _k[1].real() + _k[2].real() * _k[2].real());
		}

		/// <summary>
		/// Get the real direction of the k vector
		/// </summary>
		/// <returns></returns>
		glm::vec<3, T> getDirection() {
			T kmag = getKmag();
			return getKreal() / kmag;
		}

		/// <summary>
		/// Convert the plane wave to a string.
		/// </summary>
		/// <returns></returns>
		std::string str(int width = 25, int precision = 4) {
			std::stringstream ss;
			glm::vec<3, T> s = getDirection();
			ss << "k: " << getKmag()<<std::endl;
			ss << "s: "<<std::scientific<<std::setprecision(precision)<<std::setw(width)<< s[0] <<std::setw(width)<< s[1] <<std::setw(width)<< s[2]<<std::endl;
			glm::tvec3<std::complex<T>> E = getE0();
			ss << "E: "<<std::scientific<<std::setprecision(precision)<<std::setw(width)<< E[0] <<std::setw(width)<< E[1] <<std::setw(width)<< E[2]<<std::endl;
			return ss.str();
		}

		tira::planewave<T> wind(T x, T y, T z) {
			std::complex<T> k_dot_r = _k[0] * x + _k[1] * y + _k[2] * z;
			std::complex<T> i(0.0, 1.0);
			std::complex<T> phase = std::exp(i * k_dot_r);
			return tira::planewave<T>(_k[0], _k[1], _k[2], _E0[0] * phase, _E0[1] * phase, _E0[2] * phase);
		}

		/// Calculate the reflected and transmitted plane waves at a single boundary

		/// @param P is a plane representing the position and orientation of the surface
		/// @param r is the ratio n1/n0 (n0 is real)
		/// @param r is the reflected component of the plane wave
		/// @param t is the transmitted component of the plane wave
		void scatter(glm::vec<3, T> plane_normal, glm::vec<3, T> plane_position, std::complex<T> nr, planewave<T>& r, planewave<T>& t) {

			// the k component for the input plane wave must be real (n0 is non-absorbing)
			if(_k[0].imag() != 0.0 || _k[1].imag() != 0.0 || _k[2].imag() != 0.0)
				throw "Cannot scatter a plane wave with an imaginary k-vector";

			glm::tvec3<T> ki(_k[0].real(), _k[1].real(), _k[2].real());	// make the k vector real
			glm::tvec3<T> kr;											// the reflected wave can only have a real k vector
			glm::tvec3< std::complex<T> > kt(0.0, 0.0, 0.0);			// the transmitted plane wave can have a complex k vector
			glm::tvec3< std::complex<T> > Ei = getE0();;				// All E vectors are complex
			glm::tvec3< std::complex<T> > Er(0.0, 0.0, 0.0);
			glm::tvec3< std::complex<T> > Et(0.0, 0.0, 0.0);

 			glm::tvec3<T> n = glm::normalize(plane_normal);				// make sure that the plane normal is normalized
			glm::tvec3<T> ki_dir = glm::normalize(ki);						// calculate the direction of the k vector

			// swap the plane normal if necessary (it should be pointing in the opposite direction of ki)
			if(glm::dot(ki_dir, n) > 0) n = -n;

			// use Snll's Law to calculate the transmitted angle
			T cos_theta_i = glm::dot(ki_dir, -n);
			T sin_theta_i = std::sqrt(1 - cos_theta_i * cos_theta_i);
			T theta_i = std::acos(cos_theta_i);

			// define the basis vectors for the calculation (plane of incidence)
			glm::tvec3<T> z_hat = -n;
			glm::tvec3<T> side = n * glm::dot(ki_dir, n);
			glm::tvec3<T> diff = ki_dir - side;
			T diff_len_sq = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];

			// handle the degenerate case where theta_i is 0 (the plane wave hits head-on)
			if(diff_len_sq == 0){
				std::complex<T> rp = (1.0 - nr) / (1.0 + nr);		// compute the Fresnel coefficients (only 2 needed)
				std::complex<T> tp = 2.0 / (1.0 + nr);

				kr = -ki;
				kt[0] = ki[0] * nr;			// the transmission vector is just scaled by the refractive index change
				kt[1] = ki[1] * nr;
				kt[2] = ki[2] * nr;

				Er = Ei * rp;				// apply the Fresnel coefficients to set the amplitude of the scattered plane waves
				Et = Ei * tp;

				// calculate the phase offset based on the plane position
				T phase_r = glm::dot(plane_position, ki - kr);
				std::complex<T> phase_t = 
					plane_position[0] * (ki[0] - kt[0]) +
					plane_position[1] * (ki[1] - kt[1]) +
					plane_position[2] * (ki[2] - kt[2]);
			}
			else{
				T cos_theta_r = cos_theta_i;
				T sin_theta_r = sin_theta_i;
				T theta_r = theta_i;

				std::complex<T> sin_theta_t = (1.0/nr) * std::sin(theta_i);	// compute the sin of theta_t using Snell's law
				std::complex<T> cos_theta_t = std::sqrt(1.0 - sin_theta_t * sin_theta_t);
				std::complex<T> theta_t = std::asin(sin_theta_t);

				// calculate the remaining basis vectors
				glm::tvec3<T> y_hat = glm::normalize(ki_dir - side);
				glm::tvec3<T> x_hat = glm::cross(y_hat, z_hat);

				// calculate the k-vector magnitudes (wave frequency)
				T ki_mag = glm::length(ki);
				T kr_mag = ki_mag;
				std::complex<T> kt_mag = ki_mag * nr;

				// calculate the k vector directions
				glm::tvec3<T> kr_dir = y_hat * sin_theta_r - z_hat * cos_theta_r;
				glm::vec<3, std::complex<T>> kt_dir(0.0, 0.0, 0.0);
				kt_dir[0] = y_hat[0] * sin_theta_t + z_hat[0] * cos_theta_t;
				kt_dir[1] = y_hat[1] * sin_theta_t + z_hat[1] * cos_theta_t;
				kt_dir[2] = y_hat[2] * sin_theta_t + z_hat[2] * cos_theta_t;

				// calculate the k vectors
				ki = ki_dir * ki_mag;
				kr = kr_dir * kr_mag;
				kt = kt_dir * kt_mag;

				// calculate the Fresnel coefficients
				std::complex<T> rs = std::sin(theta_t - theta_i) / std::sin(theta_t + theta_i);
				std::complex<T> rp = std::tan(theta_t - theta_i) / std::tan(theta_t + theta_i);
				std::complex<T> ts = (2.0 * (sin_theta_t * cos_theta_i)) / std::sin(theta_t + theta_i);
				std::complex<T> tp = ((2.0 * sin_theta_t * cos_theta_i) / (std::sin(theta_t + theta_i) * std::cos(theta_t - theta_i)));

				// calculate the p component directions for each E vector
				glm::tvec3<T> Eip_dir = y_hat * cos_theta_i - z_hat * sin_theta_i;
				glm::tvec3<T> Erp_dir = y_hat * cos_theta_r + z_hat * sin_theta_r;
				glm::vec<3, std::complex<T>> Etp_dir(0.0, 0.0, 0.0);
				Etp_dir[0] = y_hat[0] * cos_theta_t - z_hat[0] * sin_theta_t;
				Etp_dir[1] = y_hat[1] * cos_theta_t - z_hat[1] * sin_theta_t;
				Etp_dir[2] = y_hat[2] * cos_theta_t - z_hat[2] * sin_theta_t;
				
				// calculate the s and t components of each E vector
				std::complex<T> Ei_s = _cDot(Ei, x_hat);
				std::complex<T> Ei_p = _cDot(Ei, Eip_dir);
				std::complex<T> Er_s = rs * Ei_s;
				std::complex<T> Er_p = rp * Ei_p;
				std::complex<T> Et_s = ts * Ei_s;
				std::complex<T> Et_p = tp * Ei_p;

				// calculate the E vector for each plane wave
				Er[0] = Erp_dir[0] * Er_p + x_hat[0] * Er_s;
				Er[1] = Erp_dir[1] * Er_p + x_hat[1] * Er_s;
				Er[2] = Erp_dir[2] * Er_p + x_hat[2] * Er_s;

				Et[0] = Etp_dir[0] * Et_p + x_hat[0] * Et_s;
				Et[1] = Etp_dir[1] * Et_p + x_hat[1] * Et_s;
				Et[2] = Etp_dir[2] * Et_p + x_hat[2] * Et_s;
			}
			// calculate the phase offset based on the plane positions
			T phase_r = glm::dot(plane_position, ki - kr);
			std::complex<T> phase_t =
				plane_position[0] * (ki[0] - kt[0]) +
				plane_position[1] * (ki[1] - kt[1]) +
				plane_position[2] * (ki[2] - kt[2]);

			// apply the phase offset
			Er = Er * std::exp(std::complex<T>(0, 1) * phase_r);
			Et = Et * std::exp(std::complex<T>(0, 1) * phase_t);

			// generate the reflected and transmitted waves
			r = planewave<T>(kr[0], kr[1], kr[2], Er[0], Er[1], Er[2]);
			t = planewave<T>(kt[0], kt[1], kt[2], Et[0], Et[1], Et[2]);
		}


	};
}