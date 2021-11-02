#ifndef TIRA_PLANEWAVE_H
#define TIRA_PLANEWAVE_H

#include <string>
#include <sstream>
#include <cmath>

#include "tira/geometry/vec3.h"
#include "tira/geometry/quaternion.h"
//#include "../math/constants.h"
#include "tira/geometry/plane.h"
#include <complex>


namespace tira{
	namespace optics{


		/// evaluate the scalar field produced by a plane wave at a point (x, y, z)

		/// @param x is the x-coordinate of the point
		/// @param y is the y-coordinate of the point
		/// @param z is the z-coordinate of the point
		/// @param A is the amplitude of the plane wave, specifically the field at (0, 0, 0)
		/// @param kx is the k-vector component in the x direction
		/// @param ky is the k-vector component in the y direction
		/// @param kz is the k-vector component in the z direction
		template<typename T>
		std::complex<T> planewave_scalar(T x, T y, T z, std::complex<T> A, T kx, T ky, T kz){
			T d = x * kx + y * ky + z * kz;						//calculate the dot product between k and p = (x, y, z) to find the distance p is along the propagation direction
			std::complex<T> di = std::complex<T>(0, d);		//calculate the phase shift that will have to be applied to propagate the wave distance d
			return A * exp(di);									//multiply the phase term by the amplitude at (0, 0, 0) to propagate the wave to p
		}

		/// evaluate the scalar field produced by a plane wave at several positions

		/// @param field is a pre-allocated block of memory that will store the complex field at all points
		/// @param N is the number of field values to be evaluated
		/// @param x is a set of x coordinates defining positions within the field (NULL implies that all values are zero)
		/// @param y is a set of y coordinates defining positions within the field (NULL implies that all values are zero)
		/// @param z is a set of z coordinates defining positions within the field (NULL implies that all values are zero)
		/// @param A is the amplitude of the plane wave, specifically the field at (0, 0, 0)
		/// @param kx is the k-vector component in the x direction
		/// @param ky is the k-vector component in the y direction
		/// @param kz is the k-vector component in the z direction
		template<typename T>
		void cpu_planewave_scalar(std::complex<T>* field, size_t N, T* x, T* y = NULL, T* z = NULL, std::complex<T> A = 1.0, T kx = 0.0, T ky = 0.0, T kz = 0.0){
			T px, py, pz;
			for(size_t i = 0; i < N; i++){										// for each element in the array
				(x == NULL) ? px = 0 : px = x[i];								// test for NULL values
				(y == NULL) ? py = 0 : py = y[i];
				(z == NULL) ? pz = 0 : pz = z[i];

				field[i] = planewave_scalar(px, py, pz, A, kx, ky, kz);			// call the single-value plane wave function
			}
		}


template<typename T>
class planewave{

protected:

	cvec3<T> m_k;					//k-vector, pointed in propagation direction with magnitude |k| = tau / lambda = 2pi / lambda
	cvec3<T> m_E;					//amplitude (for a scalar plane wave, only E0[0] is used)

	/// Bend a plane wave via refraction, given that the new propagation direction is known
	CUDA_CALLABLE planewave<T> bend(vec3<T> v) const {

		vec3<T> k_real(m_k.get(0).real(), m_k.get(1).real(), m_k.get(2).real());			//force the vector to be real (can only refract real directions)

		vec3<T> kn_hat = v.direction();					//normalize the new k
		vec3<T> k_hat = k_real.direction();				//normalize the current k

		planewave<T> new_p;								//create a new plane wave

		T k_dot_kn = k_hat.dot(kn_hat);					//if kn is equal to k or -k, handle the degenerate case

		//if k . n < 0, then the bend is a reflection
		if(k_dot_kn < 0) k_hat = -k_hat;				//flip k_hat

		if(k_dot_kn == -1){
			new_p.m_k = -m_k;
			new_p.m_E = m_E;
			return new_p;
		}
		else if(k_dot_kn == 1){
			new_p.m_k = m_k;
			new_p.m_E = m_E;
			return new_p;
		}

		vec3<T> r = k_hat.cross(kn_hat);					//compute the rotation vector
		T theta = asin(r.len());							//compute the angle of the rotation about r
		quaternion<T> q;									//create a quaternion to capture the rotation
		q.CreateRotation(theta, r.direction());	
		matrix_sq<T, 3> R = q.toMatrix3();
		vec3< std::complex<T> > E(m_E.get(0), m_E.get(1), m_E.get(2));
		vec3< std::complex<T> > E0n = R * E;					//apply the rotation to E0
		//new_p.m_k = kn_hat * kmag();
		//new_p.m_E = E0n;
		new_p.m_k[0] = kn_hat[0] * kmag();
		new_p.m_k[1] = kn_hat[1] * kmag();
		new_p.m_k[2] = kn_hat[2] * kmag();

		new_p.m_E[0] = E0n[0];
		new_p.m_E[1] = E0n[1];
		new_p.m_E[2] = E0n[2];


		return new_p;
	}

public:

	
	
	///constructor: create a plane wave propagating along k
	CUDA_CALLABLE planewave(std::complex<T> kx, std::complex<T> ky, std::complex<T> kz,
		std::complex<T> Ex, std::complex<T> Ey, std::complex<T> Ez) {

		m_k = cvec3<T>(kx, ky, kz);
		m_E = cvec3<T>(Ex, Ey, Ez);
		force_orthogonal();
	}

	CUDA_CALLABLE planewave() : planewave(0, 0, 1, 1, 0, 0) {}

	//copy constructor
	CUDA_CALLABLE planewave(const planewave& other) {
		m_k = other.m_k;
		m_E = other.m_E;
	}

	/// Assignment operator
	CUDA_CALLABLE planewave& operator=(const planewave& rhs) {
		m_k = rhs.m_k;
		m_E = rhs.m_E;

		return *this;
	}

	/// Forces the k and E vectors to be orthogonal
	CUDA_CALLABLE void force_orthogonal() {

		/*if (m_E.norm2() == 0) return;

		cvec3<T> k_dir = m_k.direction();							//calculate the normalized direction vectors for k and E
		cvec3<T> E_dir = m_E.direction();
		cvec3<T> side = k_dir.cross(E_dir);						//calculate a side vector for projection
		cvec3<T> side_dir = side.direction();					//normalize the side vector
		E_dir = side_dir.cross(k_dir);								//calculate the new E vector direction
		T E_norm = m_E.norm2();
		m_E = E_dir * E_norm;								//apply the new direction to the existing E vector
		*/
	}

	CUDA_CALLABLE cvec3<T> k() {
		return m_k;
	}

	CUDA_CALLABLE cvec3<T> E() {
		return m_E;
	}

	CUDA_CALLABLE cvec3<T> evaluate(T x, T y, T z) {
		
		std::complex<T> k_dot_r = m_k[0] * x + m_k[1] * y + m_k[2] * z;
		std::complex<T> e_k_dot_r = std::exp(std::complex<T>(0, 1) * k_dot_r);

		cvec3<T> result;
		result[0] = m_E[0] * e_k_dot_r;
		result[1] = m_E[1] * e_k_dot_r;
		result[2] = m_E[2] * e_k_dot_r;
		return result;
	}

	CUDA_CALLABLE T kmag() const {
		return std::sqrt(std::real(m_k.get(0) * std::conj(m_k.get(0)) + m_k.get(1) * std::conj(m_k.get(1)) + m_k.get(2) * std::conj(m_k.get(2))));
	}

	/// Return a plane wave with the origin translated by (x, y, z)
	CUDA_CALLABLE planewave<T> translate(T x, T y, T z) const {
		planewave<T> result;
		cvec3<T> k = m_k;
		result.m_k = k;
		std::complex<T> k_dot_r = k[0] * (-x) + k[1] * (-y) + k[2] * (-z);
		std::complex<T> exp_k_dot_r = std::exp(std::complex<T>(0.0, 1.0) * k_dot_r);

		cvec3<T> E = m_E;
		result.m_E[0] = E[0] * exp_k_dot_r;
		result.m_E[1] = E[1] * exp_k_dot_r;
		result.m_E[2] = E[2] * exp_k_dot_r;
		return result;
	}

	///multiplication operator: scale E0
    CUDA_CALLABLE planewave<T>& operator* (const T& rhs) {
		m_E = m_E * rhs;
		return *this;
	}

	///return a plane wave with the applied refractive index (scales the k-vector by the input)
	CUDA_CALLABLE planewave<T> ri(T n) {
		planewave<T> result;
		result.m_E = m_E;
		result.m_k = m_k * n;
		return result;
	}
	CUDA_CALLABLE planewave<T> refract(vec3<T> kn) const {
		return bend(kn);
	}

	/*CUDA_CALLABLE T lambda() const{
		return stim::TAU / k.len();
	}

	CUDA_CALLABLE T kmag() const{
		return k.len();
	}

	CUDA_CALLABLE vec< complex<T> > E(){
		return E0;
	}

	CUDA_CALLABLE vec<T> kvec(){
		return k;
	}

	/// calculate the value of the field produced by the plane wave given a three-dimensional position
	CUDA_CALLABLE vec< complex<T> > pos(T x, T y, T z){
		return pos( stim::vec<T>(x, y, z) );
	}

	CUDA_CALLABLE vec< complex<T> > pos(vec<T> p = vec<T>(0, 0, 0)){
		vec< complex<T> > result;

		T kdp = k.dot(p);
		complex<T> x = complex<T>(0, kdp);
		complex<T> expx = exp(x);

		result[0] = E0[0] * expx;
		result[1] = E0[1] * expx;
		result[2] = E0[2] * expx;

		return result;
	}

	//scales k based on a transition from material ni to material nt
	CUDA_CALLABLE planewave<T> n(T ni, T nt){
		return planewave<T>(k * (nt / ni), E0);
	}

	

	/// Calculate the result of a plane wave hitting an interface between two refractive indices

	/// @param P is a plane representing the position and orientation of the surface
	/// @param n0 is the refractive index outside of the surface (in the direction of the normal)
	/// @param n1 is the refractive index inside the surface (in the direction away from the normal)
	/// @param r is the reflected component of the plane wave
	/// @param t is the transmitted component of the plane wave
	void scatter(stim::plane<T> P, T n0, T n1, planewave<T> &r, planewave<T> &t){
		scatter(P, n1/n0, r, t);
	}*/

	/// Calculate the scattering result when nr = n1/n0

	/// @param P is a plane representing the position and orientation of the surface
	/// @param r is the ratio n1/n0
	/// @param n1 is the refractive index inside the surface (in the direction away from the normal)
	/// @param r is the reflected component of the plane wave
	/// @param t is the transmitted component of the plane wave
	
	int scatter(vec3<T> plane_normal, vec3<T> plane_position, std::complex<T> nr, planewave<T>& r, planewave<T>& t) {
		
		if (m_k[0].imag() != 0.0 || m_k[1].imag() != 0.0 || m_k[2].imag() != 0) {
			std::cout << "ERROR: cannot scatter a plane wave with an imaginary k-vector." << std::endl;
		}

		vec3<T> ki(m_k[0].real(), m_k[1].real(), m_k[2].real());	//force the current k vector to be real
		vec3<T> kr;
		cvec3<T> kt, Ei, Er, Et;

		plane_normal = plane_normal.direction();
		vec3<T> k_dir = ki.direction();								//calculate the direction of the incident plane wave

		//int facing = plane_face(k_dir, plane_normal);				//determine which direction the plane wave is coming in
		if (k_dir.dot(plane_normal) > 0) {							//if the wave hits the back of the plane, invert the plane and nr
			std::cout << "ERROR: k-vector intersects the wrong side of the boundary." << std::endl;
			return -1;												//the plane wave is impacting the wrong side of the surface
		}

		//use Snell's Law to calculate the transmitted angle
		T cos_theta_i = k_dir.dot(-plane_normal);					//compute the cosine of theta_i
		T sin_theta_i = std::sqrt(1 - cos_theta_i * cos_theta_i);
		T theta_i = acos(cos_theta_i);								//compute theta_i

		//handle the degenerate case where theta_i is 0 (the plane wave hits head-on)
		if (theta_i == 0) {
			std::complex<T> rp = (1.0 - nr) / (1.0 + nr);			//compute the Fresnel coefficients
			std::complex<T> tp = 2.0 / (1.0 + nr);

			kr = -ki;												//the reflection vector is the inverse of the incident vector
			kt[0] = ki[0] * nr;
			kt[1] = ki[1] * nr;
			kt[2] = ki[2] * nr;
			
			Er = m_E * rp;											//compute the E vectors based on the Fresnel coefficients
			Et = m_E * tp;

			//calculate the phase offset based on the plane positions
			T phase_r = plane_position.dot(ki - kr);
			std::complex<T> phase_t =
				plane_position[0] * (ki[0] - kt[0]) +
				plane_position[1] * (ki[1] - kt[1]) +
				plane_position[2] * (ki[2] - kt[2]);
		}
		else {
			T cos_theta_r = cos_theta_i;
			T sin_theta_r = sin_theta_i;
			T theta_r = theta_i;

			std::complex<T> sin_theta_t = (1.0/nr) * sin(theta_i);		//compute the sine of theta_t using Snell's law
			std::complex<T> cos_theta_t = std::sqrt(1.0 - sin_theta_t * sin_theta_t);
			std::complex<T> theta_t = asin(sin_theta_t);				//compute the cosine of theta_t

			//Define the basis vectors for the calculation (plane of incidence)
			vec3<T> z_hat = -plane_normal;
			vec3<T> plane_perpendicular = plane_normal * k_dir.dot(plane_normal);
			vec3<T> y_hat = (k_dir - plane_perpendicular).direction();
			vec3<T> x_hat = y_hat.cross(z_hat);

			//calculate the k-vector magnitudes
			T ki_mag = ki.norm2();
			T kr_mag = ki_mag;
			std::complex<T> kt_mag = ki_mag * nr;

			//calculate the k vector directions
			vec3<T> ki_dir = y_hat * sin_theta_i + z_hat * cos_theta_i;
			vec3<T> kr_dir = y_hat * sin_theta_r - z_hat * cos_theta_r;
			cvec3<T> kt_dir;
			kt_dir[0] = y_hat[0] * sin_theta_t + z_hat[0] * cos_theta_t;
			kt_dir[1] = y_hat[1] * sin_theta_t + z_hat[1] * cos_theta_t;
			kt_dir[2] = y_hat[2] * sin_theta_t + z_hat[2] * cos_theta_t;

			//calculate the k vectors
			ki = ki_dir * ki_mag;
			kr = kr_dir * kr_mag;
			kt = kt_dir * kt_mag;

			//calculate the Fresnel coefficients
			std::complex<T> rs = std::sin(theta_t - theta_i) / std::sin(theta_t + theta_i);
			std::complex<T> rp = std::tan(theta_t - theta_i) / std::tan(theta_t + theta_i);
			std::complex<T> ts = (2.0 * (sin_theta_t * cos_theta_i)) / std::sin(theta_t + theta_i);
			std::complex<T> tp = ((2.0 * sin_theta_t * cos_theta_i) / (std::sin(theta_t + theta_i) * std::cos(theta_t - theta_i)));

			//calculate the p component directions for each E vector
			vec3<T> Eip_dir = y_hat * cos_theta_i - z_hat * sin_theta_i;
			vec3<T> Erp_dir = y_hat * cos_theta_r + z_hat * sin_theta_r;
			cvec3<T> Etp_dir;
			Etp_dir[0] = y_hat[0] * cos_theta_t - z_hat[0] * sin_theta_t;
			Etp_dir[1] = y_hat[1] * cos_theta_t - z_hat[1] * sin_theta_t;
			Etp_dir[2] = y_hat[2] * cos_theta_t - z_hat[2] * sin_theta_t;

			//calculate the s and t components of each E vector
			std::complex<T> Ei_s = m_E.dot(x_hat);
			std::complex<T> Ei_p = m_E.dot(Eip_dir);
			std::complex<T> Er_s = rs * Ei_s;
			std::complex<T> Er_p = rp * Ei_p;
			std::complex<T> Et_s = ts * Ei_s;
			std::complex<T> Et_p = tp * Ei_p;

			//calculate the E vector for each plane wave
			Er[0] = Erp_dir[0] * Er_p + x_hat[0] * Er_s;
			Er[1] = Erp_dir[1] * Er_p + x_hat[1] * Er_s;
			Er[2] = Erp_dir[2] * Er_p + x_hat[2] * Er_s;

			Et[0] = Etp_dir[0] * Et_p + x_hat[0] * Et_s;
			Et[1] = Etp_dir[1] * Et_p + x_hat[1] * Et_s;
			Et[2] = Etp_dir[2] * Et_p + x_hat[2] * Et_s;
		}


		//calculate the phase offset based on the plane positions
		T phase_r = plane_position.dot(ki - kr);
		std::complex<T> phase_t =
			plane_position[0] * (ki[0] - kt[0]) +
			plane_position[1] * (ki[1] - kt[1]) +
			plane_position[2] * (ki[2] - kt[2]);

		//apply the phase offset
		Er = Er * std::exp(std::complex<T>(0, 1) * phase_r);
		Et = Et * std::exp(std::complex<T>(0, 1) * phase_t);

		//generate the reflected and transmitted waves
		r = planewave<T>(kr[0], kr[1], kr[2], Er[0], Er[1], Er[2]);
		t = planewave<T>(kt[0], kt[1], kt[2], Et[0], Et[1], Et[2]);

		return 0;
	}

	std::string str()
	{
		std::stringstream ss;
		ss << "k: " << m_k << std::endl;
		ss << "E: " << m_E << std::endl;
		return ss.str();
	}
};					//end planewave class
}					//end namespace optics
}					//end namespace tira

template <typename T>
std::ostream& operator<<(std::ostream& os, tira::optics::planewave<T> p)
{
    os<<p.str();
    return os;
}

#endif