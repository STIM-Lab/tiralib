#ifndef TIRA_QUATERNION_H
#define TIRA_QUATERNION_H

#include <tira/geometry/matrix_sq.h>
#include <tira/cuda/callable.h>

namespace tira {

template<typename T>
class quaternion
{
public:
	T w;
	T x;
	T y;
	T z;

	CUDA_CALLABLE void normalize(){

		double length=sqrt(w*w + x*x + y*y + z*z);
		w=w/length;
		x=x/length;
		y=y/length;
		z=z/length;
	}

	//calculate the quaternion length (norm)
	CUDA_CALLABLE T norm() {
		return sqrt(w*w + x * x + y * y + z * z);
	}

	CUDA_CALLABLE void CreateRotation(T theta, T ux, T uy, T uz){

		vec3<T> u(ux, uy, uz);
		CreateRotation(theta, u);		
	}

	CUDA_CALLABLE void CreateRotation(T theta, vec3<T> u){

		vec3<T> u_hat = u.norm();

		//assign the given Euler rotation to this quaternion
		w = (T)cos(theta/2);
		x = u_hat[0]*(T)sin(theta/2);
		y = u_hat[1]*(T)sin(theta/2);
		z = u_hat[2]*(T)sin(theta/2);
	}

	CUDA_CALLABLE void CreateRotation(vec3<T> from, vec3<T> to){

		from = from.norm();
		to = to.norm();
		vec3<T> r = from.cross(to);			//compute the rotation vector
		//T l = r.len();
		//if (l > 1) l = 1;					//we have seen degenerate cases where |r| > 1 (probably due to loss of precision in the cross product)
		//T theta = asin(l);				//compute the angle of the rotation about r
		//deal with a zero vector (both k and kn point in the same direction)
		T cos_theta = from.dot(to);			//calculate the cosine between the two vectors
		if (cos_theta == (T)(-1)) {
			y = -1;
			return;
		}
		T theta = acos(cos_theta);			//calculate the angle between the two vectors
		if(theta == (T)0){
			return;
		}

		//create a quaternion to capture the rotation
		CreateRotation(theta, r.norm());
	}



	CUDA_CALLABLE quaternion<T> operator *(quaternion<T> &param){

		float A, B, C, D, E, F, G, H;


		A = (w + x)*(param.w + param.x);
		B = (z - y)*(param.y - param.z);
		C = (w - x)*(param.y + param.z);
		D = (y + z)*(param.w - param.x);
		E = (x + z)*(param.x + param.y);
		F = (x - z)*(param.x - param.y);
		G = (w + y)*(param.w - param.z);
		H = (w - y)*(param.w + param.z);

		quaternion<T> result;
		result.w = B + (-E - F + G + H) /2;
		result.x = A - (E + F + G + H)/2;
		result.y = C + (E - F + G - H)/2;
		result.z = D + (E - F - G + H)/2;

		return result;
	}
	
	CUDA_CALLABLE matrix_sq<T, 3> toMatrix3(){

		matrix_sq<T, 3> R;

		T s, wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2;
		s = sqrt(norm());
		xx = x * x;		xy = x * y;		xz = x * z;
		yy = y * y;		yz = y * z;
		zz = z * z;
		wx = w * x;		wy = w * y;		wz = w * z;

		R(0, 0) = 1 - 2 * s * (yy + zz);
		R(0, 1) = 2 * s * (xy - wz);
		R(0, 2) = 2 * s * (xz + wy);
		R(1, 0) = 2 * s * (xy + wz);
		R(1, 1) = 1 - 2 * s * (xx + zz);
		R(1, 2) = 2 * s * (yz - wx);
		R(2, 0) = 2 * s * (xz - wy);
		R(2, 1) = 2 * s * (yz + wx);
		R(2, 2) = 1 - 2 * s * (xx + yy);

		return R;

	    /*T wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2;


	    // calculate coefficients
	    x2 = x + x; y2 = y + y;
	    z2 = z + z;
	    xx = x * x2; xy = x * y2; xz = x * z2;
	    yy = y * y2; yz = y * z2; zz = z * z2;
	    wx = w * x2; wy = w * y2; wz = w * z2;

		result(0, 0) = 1 - (yy + zz);
		result(0, 1) = xy - wz;

		result(0, 2) = xz + wy;

		result(1, 0) = xy + wz;
		result(1, 1) = 1 - (xx + zz);

		result(1, 2) = yz - wx;

		result(2, 0) = xz - wy;
		result(2, 1) = yz + wx;

		result(2, 2) = 1 - (xx + yy);
		
		return result;
		*/
	}

	CUDA_CALLABLE matrix_sq<T, 4> toMatrix4(){
		matrix_sq<T, 4> R;
	    T s, wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2;
		s = sqrt(norm());
		xx = x * x;		xy = x * y;		xz = x * z;
		yy = y * y;		yz = y * z;
		zz = z * z;
		wx = w * x;		wy = w * y;		wz = w * z;

		R(0, 0) = 1 - 2 * s * (yy + zz);
		R(0, 1) = 2 * s * (xy - wz);
		R(0, 2) = 2 * s * (xz + wy);
		R(1, 0) = 2 * s * (xy + wz);
		R(1, 1) = 1 - 2 * s * (xx + zz);
		R(1, 2) = 2 * s * (yz - wx);
		R(2, 0) = 2 * s * (xz - wy);
		R(2, 1) = 2 * s * (yz + wx);
		R(2, 2) = 1 - 2 * s * (xx + yy);

		R(0, 3) = 0;
		R(1, 3) = 0;
		R(2, 3) = 0;
		R(3, 0) = 0;
		R(3, 1) = 0;
		R(3, 2) = 0;
		R(3, 3) = 1;

	    // calculate coefficients
	    /*x2 = x + x; y2 = y + y;
	    z2 = z + z;
	    xx = x * x2; xy = x * y2; xz = x * z2;
	    yy = y * y2; yz = y * z2; zz = z * z2;
	    wx = w * x2; wy = w * y2; wz = w * z2;

		result(0, 0) = 1 - (yy + zz);
		result(0, 1) = xy - wz;

		result(0, 2) = xz + wy;

		result(1, 0) = xy + wz;
		result(1, 1) = 1 - (xx + zz);

		result(1, 2) = yz - wx;

		result(2, 0) = xz - wy;
		result(2, 1) = yz + wx;

		result(2, 2) = 1 - (xx + yy);

		result(3, 3) = 1;*/

		return R;
	}


	CUDA_CALLABLE quaternion(){
		w=0; x=0; y=0; z=0;
	}

	CUDA_CALLABLE quaternion(T c, T i, T j, T k){
		w=c;  x=i;  y=j;  z=k;
	}

	// create a pure quaternion from a vector
	CUDA_CALLABLE quaternion(vec3<T> v){
		w = 0; x = v[0]; y = v[1]; z = v[2];
	}

	CUDA_CALLABLE quaternion<T> conj(){
		quaternion<T> c;
		c.w = w;
		c.x = -x;
		c.y = -y;
		c.z = -z;
		return c;
	}

};

}	//end rts namespace


#endif
