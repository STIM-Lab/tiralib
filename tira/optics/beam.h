#ifndef TIRA_BEAM
#define TIRA_BEAM

//#include "tira/geometry/vec3.h"
#include "tira/optics/planewave.h"
#include <vector>
#include <numbers>
#include <random>

namespace tira{
	namespace optics{
template<typename T>
class beam {

	glm::vec<3, T> m_direction;
	glm::vec<3, T> m_focus;
	glm::vec<3, std::complex<T>> m_E;
	T m_k;



public:
	/// <summary>
	/// Beam constructor
	/// </summary>
	/// <param name="d">Propagation direction of the beam (automatically normalized).</param>
	/// <param name="f">Focal point for the beam.</param>
	/// <param name="pol">Beam polarization (automatically normalized and orthogonalized). Only works for linear polarization, though.</param>
	/// <param name="E">Complex beam amplitude (applied to the polarization to create a linearly polarized beam).</param>
	/// <param name="lambda">Wavelength</param>
	beam(
		vec3<T> d = vec3<T>(0, 0, 1), 
		vec3<T> f = vec3<T>(0, 0, 0),
		vec3<T> pol = vec3<T>(1, 0, 0),
		std::complex<T> E = std::complex<T>(1.0, 0.0),
		T lambda = 1.0
	) {

		m_direction = d.direction();				//normalize and set the plane wave direction
		pol = (pol - d.dot(pol)).direction();			//orthogonalize and normalize the polarization vector
		m_focus = f;								//set the focal point
		m_E[0] = pol[0] * E;						//convert the polarization and complex amplitude to a linearly-polarized E vector
		m_E[1] = pol[1] * E;
		m_E[2] = pol[2] * E;
		m_k = 2.0 * std::numbers::pi / lambda;		//calculate and store the wavenumber
	}

	/// <summary>
	/// Generate a focal spot for the beam using Monte-Carlo sampling.
	/// </summary>
	/// <param name="NA">Numerical aperture of the focusing optics.</param>
	/// <param name="samples">Number of Monte-Carlo samples.</param>
	/// <param name="NA_in">NA of the internal obscuration.</param>
	/// <returns></returns>
	std::vector< planewave<T> > mcfocus(T NA, size_t samples, T NA_in = 0) {

		//create a rotation matrix to orient the beam along the specified direction vector
		tira::matrix_sq<double, 3> R;// = tira::matrix_sq<double, 3>::identity();		//initialize the rotation matrix to the identity matrix
		tira::vec3<double> Z(0, 0, 1);												//create a vector along the Z direction
		tira::quaternion<double> q;
		q.CreateRotation(Z, m_direction);
		R = q.toMatrix3();

		std::vector< planewave<T> > result(samples);								//generate an array of N plane waves
		double alpha = std::asin(NA);												//calculate the angle of the band-pass
		double alpha_in = std::asin(NA_in);											//calculate the angle of the central obscuration
		double cos_alpha = std::cos(alpha);											//calculate the cosine of the band-pass and central obscuration
		double cos_alpha_in = std::cos(alpha_in);

		//random sampling is done using Archimede's method - uniform samples are taken in cylindrical space and projected onto a sphere
		std::random_device rd;														//create a random number generator
		std::default_random_engine eng(rd());
		std::uniform_real_distribution<> theta_dist(0.0, 2.0 * std::numbers::pi);	//create a distribution for theta (in cylindrical coordinates)
		std::uniform_real_distribution<> z_dist(cos_alpha, cos_alpha_in);			//create a uniform distribution for sampling the z-axis

		//generate a guide plane wave
		planewave<T> pw(m_direction[0] * m_k, m_direction[1] * m_k, m_direction[2] * m_k, m_E[0], m_E[1], m_E[2]);
		if (samples == 1) {
			result[0] = pw;
			return result;
		}

		double theta, phi;															//variables to store the coordinates for theta/phi on a sphere
		tira::vec3<T> dir;
		T w = 1.0 / samples;														//integration weight
		for (size_t n = 0; n < samples; n++) {
			theta = theta_dist(eng);												//randomly generate a theta coordinate between 0 and 2pi
			phi = std::acos(z_dist(eng));											//randomly generate a z coordinate based on the NA and project to phi
			dir = R * tira::vec3<double>(1.0, theta, phi).sph2cart();				//convert the value to Cartesian coordinates and store the sample vector
			result[n] = pw.refract(dir);											//refract the plane wave to align with the sample direction
			result[n] = result[n].translate(m_focus[0], m_focus[1], m_focus[2]);	//refocus the plane wave to the focal position
			result[n] = result[n] * w;												//apply the integration weight (1/N)
		}

		return result;																//return the sample array
	}
};
}								//end namespace optics
}								//end namespace tira

#endif
