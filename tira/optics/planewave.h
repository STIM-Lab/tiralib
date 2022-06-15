#pragma once

#include "glm/glm.hpp"
#include "glm/gtc/quaternion.hpp"

#include <complex>


namespace tira {
	template <typename T>
	class planewave {
		std::complex<T> m_E0[2];		// Ex and Ey (complex amplitude of the plane wave)
		std::complex<T> m_k[3];		// complex k vector (direction of propagation) in 3D

		/// <summary>
		/// Generate a rotation matrix used to transform coordinates relative to k into global coordinates
		/// </summary>
		/// <returns></returns>
		glm::tmat3x3<T> Rot() {
			glm::tvec3<T> z(0.0, 0.0, 1.0);						// z axis (local coordinates)
			glm::tvec3<T> k(
				m_k[0].real(), 
				m_k[1].real(), 
				m_k[2].real());			// convert k to a vector in global coordinates

			glm::tvec3<T> k_hat = glm::normalize(k);			// convert k to a unit vector

			glm::tvec3<T> cp = glm::cross(z, k_hat);			// calculate the cross product of z and k_hat (giving us the axis of rotation and angle)
			glm::tvec3<T> axis = glm::normalize(cp);			// calculate the axis of rotation (normalized cross product)
			T angle = std::asin(glm::length(cp));				// calculate the angle of rotation (inverse sin of the length of the cross product)
			
			glm::tquat<T> q = glm::angleAxis(angle, axis);		// generate a quaternion from the rotation properties

			return glm::mat3_cast(q);							// generate and return a rotation matrix
		}

		glm::tvec3<std::complex<T>> cMul(glm::tmat3x3<T> M, 
			std::complex<T> x, 
			std::complex<T> y, 
			std::complex<T> z) {

			glm::tvec3<T> r(m_E0[0].real(), m_E0[1].real(), 0.0);				// create a vector from the real coordinates
			glm::tvec3<T> i(m_E0[1].imag(), m_E0[1].imag(), 0.0);				// create a vector from the imaginary coordinates

			glm::tvec3<T> R = M * r;											// calculate the real components of the E vector
			glm::tvec3<T> I = M * i;											// calculate the imaginary components of the E vector

			glm::tvec3<std::complex<T>> v(										// create the complex vector from real and imaginary parts
				std::complex<T>(R.x, I.x),
				std::complex<T>(R.y, I.y),
				std::complex<T>(R.z, I.z)
			);
			return v;															// return the complex 3D vector result

		}

	public:
		planewave<T>(std::complex<T> kx, std::complex<T> ky, std::complex<T> kz, std::complex<T> Ex, std::complex<T> Ey) {
			m_k[0] = kx;
			m_k[1] = ky;
			m_k[2] = kz;
			m_E0[0] = Ex;
			m_E0[1] = Ey;

		}

		/// <summary>
		/// Returns the 3-dimensional E vector at zero phase
		/// </summary>
		/// <returns></returns>
		glm::vec<3, std::complex<T>> getE0() {

			glm::tmat3x3<T> M = Rot();							// get the matrix mapping local coordinates into global coordinates
			return cMul(M, m_E0[0], m_E0[1], 0.0);				// multiply the local E0 coordinates by the rotation matrix
			//return M * glm::tvec3<T>(1.0, 1.0, 1.0);
		}

		/// <summary>
		/// Convert the plane wave to a string.
		/// </summary>
		/// <returns></returns>
		std::string str() {
			std::stringstream ss;
			ss << "k: [" << m_k[0] << std::endl << m_k[1] << std::endl << m_k[2] << "]" << std::endl;
			ss << std::endl;
			glm::tvec3<std::complex<T>> E = getE0();
			ss << "E: [" << E.x << std::endl << E.y << std::endl << E.z << "]" << std::endl;
			return ss.str();
		}
	};
}