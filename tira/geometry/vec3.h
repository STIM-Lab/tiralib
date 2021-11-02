#ifndef TIRA_VEC3_H
#define TIRA_VEC3_H


#include <stim/cuda/cudatools/callable.h>
#include <complex>
#include <cmath>

#include <sstream>


namespace tira{


/// A class designed to act as a 3D vector with CUDA compatibility
template<typename T>
class vec3{

protected:
	T ptr[3];

public:

	CUDA_CALLABLE vec3(){}

	CUDA_CALLABLE vec3(T v){
		ptr[0] = ptr[1] = ptr[2] = v;
	}

	CUDA_CALLABLE vec3(T x, T y, T z){
		ptr[0] = x;
		ptr[1] = y;
		ptr[2] = z;
	}

	//copy constructor
	CUDA_CALLABLE vec3( const vec3<T>& other){
		ptr[0] = other.ptr[0];
		ptr[1] = other.ptr[1];
		ptr[2] = other.ptr[2];
	}

	//access an element using an index
	CUDA_CALLABLE T& operator[](size_t idx){
		return ptr[idx];
	}
	//read only accessor method
	CUDA_CALLABLE T get(size_t idx) const {
		return ptr[idx];
	}

	CUDA_CALLABLE T* data(){
		return ptr;
	}

/// Casting operator. Creates a new vector with a new type U.
/*	template< typename U >
	CUDA_CALLABLE operator vec3<U>(){
		vec3<U> result((U)ptr[0], (U)ptr[1], (U)ptr[2]);
		//result.ptr[0] = (U)ptr[0];
		//result.ptr[1] = (U)ptr[1];
		//result.ptr[2] = (U)ptr[2];

		return result;
	}*/

	// computes the squared Euclidean length (useful for several operations where only >, =, or < matter)
	CUDA_CALLABLE T len_sq() const{
		return ptr[0] * ptr[0] + ptr[1] * ptr[1] + ptr[2] * ptr[2];
	}

	/// computes the Euclidean length of the vector
	CUDA_CALLABLE T len() const{
		return sqrt(len_sq());
	}

	/// Calculate the L2 norm of the vector, accounting for the possibility that it is complex
	CUDA_CALLABLE T norm2() const{
		T conj_dot = std::real(ptr[0] * std::conj(ptr[0]) + ptr[1] * std::conj(ptr[1]) + ptr[2] * std::conj(ptr[2]));
		return sqrt(conj_dot);
	}

	/// Calculate the normalized direction vector
	CUDA_CALLABLE vec3<T> direction() const {
		vec3<T> result;
		T length = norm2();
		result[0] = ptr[0] / length;
		result[1] = ptr[1] / length;
		result[2] = ptr[2] / length;
		return result;
	}

	/// Convert the vector from cartesian to spherical coordinates (x, y, z -> r, theta, phi where theta = [-PI, PI])
	CUDA_CALLABLE vec3<T> cart2sph() const{
		vec3<T> sph;
		sph.ptr[0] = len();
		sph.ptr[1] = std::atan2(ptr[1], ptr[0]);
		if(sph.ptr[0] == 0)
			sph.ptr[2] = 0;
		else
			sph.ptr[2] = std::acos(ptr[2] / sph.ptr[0]);
		return sph;
	}

	CUDA_CALLABLE vec3<T> cart2cyl() const{
		vec3<T> cyl;
		cyl.ptr[0] = sqrt(pow(ptr[0],2) + pow(ptr[1],2));
		cyl.ptr[1] = atan(ptr[1]/ptr[0]);
		cyl.ptr[2] = ptr[2];

		return cyl;
	}

	/// Convert the vector from cartesian to spherical coordinates (r, theta, phi -> x, y, z where theta = [0, 2*pi])
	CUDA_CALLABLE vec3<T> sph2cart() const{
		vec3<T> cart;
		cart.ptr[0] = ptr[0] * std::cos(ptr[1]) * std::sin(ptr[2]);
		cart.ptr[1] = ptr[0] * std::sin(ptr[1]) * std::sin(ptr[2]);
		cart.ptr[2] = ptr[0] * std::cos(ptr[2]);

		return cart;
	}

	/// Convert the vector from cylindrical to cart coordinates (r, theta, z -> x, y, z where theta = [0, 2*pi])
	CUDA_CALLABLE vec3<T> cyl2cart() const{
		vec3<T> cart;
		cart.ptr[0] = ptr[0] * std::cos(ptr[1]);
		cart.ptr[1] = ptr[0] * std::sin(ptr[1]);
		cart.ptr[2] = ptr[2];

		return cart;
	}

	/// Computes the normalized vector (where each coordinate is divided by the L2 norm)
	CUDA_CALLABLE vec3<T> norm() const{
        vec3<T> result;
        T l = len();						//compute the vector length
        return (*this) / l;
	}

	/// Computes the cross product of a 3-dimensional vector
	CUDA_CALLABLE vec3<T> cross(const vec3<T> rhs) const{

		vec3<T> result;

		result[0] = (ptr[1] * rhs.ptr[2] - ptr[2] * rhs.ptr[1]);
		result[1] = (ptr[2] * rhs.ptr[0] - ptr[0] * rhs.ptr[2]);
		result[2] = (ptr[0] * rhs.ptr[1] - ptr[1] * rhs.ptr[0]);

		return result;
	}

	/// Compute the Euclidean inner (dot) product
    CUDA_CALLABLE T dot(vec3<T> rhs) const{
        return ptr[0] * rhs.ptr[0] + ptr[1] * rhs.ptr[1] + ptr[2] * rhs.ptr[2];
    }

	/// Arithmetic addition operator

    /// @param rhs is the right-hand-side operator for the addition
	CUDA_CALLABLE vec3<T> operator+(vec3<T> rhs) const{
		vec3<T> result;
		result.ptr[0] = ptr[0] + rhs[0];
		result.ptr[1] = ptr[1] + rhs[1];
		result.ptr[2] = ptr[2] + rhs[2];
		return result;
	}

	/// Arithmetic addition to a scalar

	/// @param rhs is the right-hand-side operator for the addition
	CUDA_CALLABLE vec3<T> operator+(T rhs) const{
		vec3<T> result;
		result.ptr[0] = ptr[0] + rhs;
		result.ptr[1] = ptr[1] + rhs;
		result.ptr[2] = ptr[2] + rhs;
		return result;
	}

	/// Arithmetic subtraction operator

	/// @param rhs is the right-hand-side operator for the subtraction
	CUDA_CALLABLE vec3<T> operator-(vec3<T> rhs) const{
		vec3<T> result;
		result.ptr[0] = ptr[0] - rhs[0];
		result.ptr[1] = ptr[1] - rhs[1];
		result.ptr[2] = ptr[2] - rhs[2];
		return result;
	}
	/// Arithmetic subtraction to a scalar

	/// @param rhs is the right-hand-side operator for the addition
	CUDA_CALLABLE vec3<T> operator-(T rhs) const{
		vec3<T> result;
		result.ptr[0] = ptr[0] - rhs;
		result.ptr[1] = ptr[1] - rhs;
		result.ptr[2] = ptr[2] - rhs;
		return result;
	}

	/// Arithmetic scalar multiplication operator

	/// @param rhs is the right-hand-side operator for the subtraction
	CUDA_CALLABLE vec3<T> operator*(T rhs) const{
		vec3<T> result;
		result.ptr[0] = ptr[0] * rhs;
		result.ptr[1] = ptr[1] * rhs;
		result.ptr[2] = ptr[2] * rhs;
		return result;
	}

	/// Arithmetic scalar division operator

	/// @param rhs is the right-hand-side operator for the subtraction
	CUDA_CALLABLE vec3<T> operator/(T rhs) const{
		return (*this) * ((T)1.0/rhs);
	}

	/// Multiplication by a scalar, followed by assignment
	CUDA_CALLABLE vec3<T> operator*=(T rhs){
		ptr[0] = ptr[0] * rhs;
		ptr[1] = ptr[1] * rhs;
		ptr[2] = ptr[2] * rhs;
		return *this;
	}

	/// Addition and assignment
	CUDA_CALLABLE vec3<T> operator+=(vec3<T> rhs){
		ptr[0] = ptr[0] + rhs;
		ptr[1] = ptr[1] + rhs;
		ptr[2] = ptr[2] + rhs;
		return *this;
	}

	/// Assign a scalar to all values
	CUDA_CALLABLE vec3<T> & operator=(T rhs){
		ptr[0] = ptr[0] = rhs;
		ptr[1] = ptr[1] = rhs;
		ptr[2] = ptr[2] = rhs;
		return *this;
	}

	/// Casting and assignment
	template<typename Y>
	CUDA_CALLABLE vec3<T> & operator=(vec3<Y> rhs){
		ptr[0] = (T)rhs.ptr[0];
		ptr[1] = (T)rhs.ptr[1];
		ptr[2] = (T)rhs.ptr[2];
		return *this;
	}

	/// Unary minus (returns the negative of the vector)
	CUDA_CALLABLE vec3<T> operator-() const{
		vec3<T> result;
		result.ptr[0] = -ptr[0];
		result.ptr[1] = -ptr[1];
		result.ptr[2] = -ptr[2];
		return result;
	}

	CUDA_CALLABLE bool operator==(vec3<T> rhs) const{
		if(rhs[0] == ptr[0] && rhs[1] == ptr[1] && rhs[2] == ptr[2])
			return true;
		else
			return false;	
	}

//#ifndef __CUDACC__
	/// Outputs the vector as a string
std::string str() const{
		std::stringstream ss;

		const size_t N = 3;

		ss<<"[";
		for(size_t i=0; i<N; i++)
		{
			ss<<ptr[i];
			if(i != N-1)
				ss<<", ";
		}
		ss<<"]";

		return ss.str();
	}
//#endif

	size_t size(){ return 3; }

	};						//end class vec3

	/// Start Class cvec3 (complex vector class)
	template<typename T>
	class cvec3 : public vec3< std::complex<T> > {

	public:
		CUDA_CALLABLE cvec3() : vec3< std::complex<T> >() {}
		CUDA_CALLABLE cvec3(T x, T y, T z) : vec3< std::complex<T> >(std::complex<T>(x), std::complex<T>(y), std::complex<T>(z)) {}
		CUDA_CALLABLE cvec3(std::complex<T> x, std::complex<T> y, std::complex<T> z) : vec3< std::complex<T> >(x, y, z) {}

		CUDA_CALLABLE cvec3& operator=(const vec3< std::complex<T> > rhs) {
			vec3< std::complex<T> >::ptr[0] = rhs.get(0);
			vec3< std::complex<T> >::ptr[1] = rhs.get(1);
			vec3< std::complex<T> >::ptr[2] = rhs.get(2);
			return *this;
		}
		CUDA_CALLABLE std::complex<T> dot(const vec3<T> rhs) {
			std::complex<T> result =
				vec3< std::complex<T> >::ptr[0] * rhs.get(0) +
				vec3< std::complex<T> >::ptr[1] * rhs.get(1) +
				vec3< std::complex<T> >::ptr[2] * rhs.get(2);

			return result;
		}

		/// @param rhs is the right-hand-side operator for the addition
		CUDA_CALLABLE cvec3<T> operator+(cvec3<T> rhs) const {
			return vec3< std::complex<T> >::operator+(rhs);
		}
		/// @param rhs is the right-hand-side operator for adding a real vector to a complex vector
		CUDA_CALLABLE cvec3<T> operator+(vec3<T> rhs) const {
			cvec3<T> result;
			result.ptr[0] = vec3< std::complex<T> >::ptr[0] + rhs[0];
			result.ptr[1] = vec3< std::complex<T> >::ptr[1] + rhs[1];
			result.ptr[2] = vec3< std::complex<T> >::ptr[2] + rhs[2];
			return result;
		}
	};						//end class cvec3
}							//end namespace tira




/// Multiply a vector by a constant when the vector is on the right hand side
template <typename T>
tira::vec3<T> operator*(T lhs, tira::vec3<T> rhs){
    return rhs * lhs;
}

//stream operator
template<typename T>
std::ostream& operator<<(std::ostream& os, tira::vec3<T> const& rhs){
	os<<rhs.str();
	return os;
}

#endif
