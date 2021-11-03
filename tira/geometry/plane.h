#ifndef TIRA_PLANE_H
#define TIRA_PLANE_H

#include <iostream>
#include <tira/geometry/vec3.h>
#include <tira/cuda/callable.h>
#include <tira/geometry/quaternion.h>


namespace tira
{
template<typename T> class plane;
}

template<typename T>
CUDA_CALLABLE tira::plane<T> operator-(tira::plane<T> v);

namespace tira
{

template <typename T>
class plane
{
	protected:
		vec3<T> P;					//position of a point on the plane
		vec3<T> N;					//plane normal
		vec3<T> U;					//direction vector for U,V coordinates

		///Initializes the plane with standard coordinates.
		///
		CUDA_CALLABLE void init()
		{
			P = vec3<T>(0, 0, 0);
			N = vec3<T>(0, 0, 1);
			U = vec3<T>(1, 0, 0);
		}

	public:
	
		CUDA_CALLABLE plane()
		{
			init();
		}

		CUDA_CALLABLE plane(vec3<T> n, vec3<T> p = vec3<T>(0, 0, 0))
		{
			init();
			P = p;
			rotate(n.direction());
		}

		CUDA_CALLABLE plane(T z_pos)
		{
			init();
			P[2] = z_pos;
		}

		//create a plane from three points (a triangle)
		CUDA_CALLABLE plane(vec3<T> a, vec3<T> b, vec3<T> c)
		{
			init();
			P = c;
			vec3<T> n = (c - a).cross(b - a);
			try
			{
				if(n.len() != 0)
				{
					rotate(n.norm());
				} else {
				 	throw 42;
				}
			}
			catch(int i)
			{
				std::cerr << "No plane can be creates as all points a,b,c lie on a straight line" << std::endl;
			}  
		}
	
		template< typename U >
		CUDA_CALLABLE operator plane<U>()
		{
			plane<U> result(N, P);
			return result;

		}

		CUDA_CALLABLE vec3<T> n()
		{
			return N;
		}

		CUDA_CALLABLE vec3<T> p()
		{
			return P;
		}

		CUDA_CALLABLE vec3<T> u()
		{
			return U;
		}

		//return a 3D coordinate from UV coordinates on the plane
		CUDA_CALLABLE vec3<T> operator()(T u, T v) {
			vec3<T> V = N.cross(U);
			return P + U * u + V * v;
		}

		//return a 3D coordinate from UVW coordinates (assuming that +W is the normal direction)
		CUDA_CALLABLE vec3<T> operator()(T u, T v, T w) {
			return operator()(u, v) + N * w;
		}

		///flip the plane front-to-back
		CUDA_CALLABLE plane<T> flip(){
			plane<T> result = *this;
			result.N = -result.N;
			return result;
		}

		//determines how a vector v intersects the plane (1 = intersects front, 0 = within plane,     -1 = intersects back)
		CUDA_CALLABLE int face(vec3<T> v){
			
			T dprod = v.dot(N);             //get the dot product between v and N

			//conditional returns the appropriate value
			if(dprod < 0)
				return 1;
			else if(dprod > 0)
				return -1;
			else
				return 0;
		}

		//determine on which side of the plane a point lies (1 = front, 0 = on the plane, -1 = bac    k)
		CUDA_CALLABLE int side(vec3<T> p){

			vec3<T> v = p - P;    //get the vector from P to the query point p

			return face(v);
		}

		//compute the component of v that is perpendicular to the plane
		CUDA_CALLABLE vec3<T> perpendicular(vec3<T> v){
			return N * v.dot(N);
		}

		//compute the projection of v in the plane
		CUDA_CALLABLE vec3<T> parallel(vec3<T> v){
			return v - perpendicular(v);
		}

		CUDA_CALLABLE void setU(vec3<T> v)
		{
			U = (parallel(v.norm())).norm();		
		}

		CUDA_CALLABLE void decompose(vec3<T> v, vec3<T>& para, vec3<T>& perp){
			perp = N * v.dot(N);
			para = v - perp;
		}

		//get both the parallel and perpendicular components of a vector v w.r.t. the plane
		CUDA_CALLABLE void project(vec3<T> v, vec3<T> &v_par, vec3<T> &v_perp){

			v_perp = v.dot(N);
			v_par = v - v_perp;
		}

		//compute the reflection of v off of the plane
		CUDA_CALLABLE vec3<T> reflect(vec3<T> v){

			//compute the reflection using N_prime as the plane normal
			vec3<T> par = parallel(v);
			vec3<T> r = (-v) + par * 2;
			return r;

		}

		CUDA_CALLABLE plane<T> operator-()
		{
			plane<T> p = *this;

			//negate the normal vector
			p.N = -p.N;
			return p;
		}

		//output a string
		std::string str(){
			std::stringstream ss;
			ss<<"P: "<<P<<std::endl;
			ss<<"N: "<<N<<std::endl;
			ss<<"U: "<<U;
			return ss.str();
		}


		CUDA_CALLABLE void rotate(vec3<T> n)
		{
			quaternion<T> q;
			q.CreateRotation(N, n);
			matrix_sq<T, 3> M = q.toMatrix3();
			N = M * N;
			U = M * U;

		}

		CUDA_CALLABLE void rotate(vec3<T> n, vec3<T> &Y)
		{
			quaternion<T> q;
			q.CreateRotation(N, n);
			
			N = q.toMatrix3() * N;
			U = q.toMatrix3() * U;
			Y = q.toMatrix3() * Y;

		}

		CUDA_CALLABLE void rotate(vec3<T> n, vec3<T> &X, vec3<T> &Y)
		{
			quaternion<T> q;
			q.CreateRotation(N, n);
			
			N = q.toMatrix3() * N;
			U = q.toMatrix3() * U;
			X = q.toMatrix3() * X;
			Y = q.toMatrix3() * Y;

		}

};
		
		
}
#endif
