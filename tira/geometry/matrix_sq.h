#ifndef TIRA_MATRIX_SQ_H
#define TIRA_MATRIX_SQ_H

//#include "rts/vector.h"
#include <string.h>
#include <iostream>
#include <tira/geometry/vec3.h>
#include <tira/cuda/callable.h>

namespace tira {

template <class T, int N>
struct matrix_sq
{
	//the matrix will be stored in column-major order (compatible with OpenGL)
	T M[N*N];

	CUDA_CALLABLE matrix_sq()
	{
		for(int r=0; r<N; r++)
			for(int c=0; c<N; c++)
				if(r == c)
					(*this)(r, c) = 1;
				else
					(*this)(r, c) = 0;
	}

	CUDA_CALLABLE matrix_sq(T rhs[N*N])
	{
		memcpy(M,rhs, sizeof(T)*N*N);
	}

	CUDA_CALLABLE matrix_sq<T,N> set(T rhs[N*N])
	{
		memcpy(M, rhs, sizeof(T)*N*N);
		return *this;
	}

	//create a symmetric matrix given the rhs values, given in column-major order
	CUDA_CALLABLE void setsym(T rhs[(N*N+N)/2]){
		const size_t L = (N*N+N)/2;		//store the number of values

		size_t r, c;
		r = c = 0;
		for(size_t i = 0; i < L; i++){ 				//for each value
			if(r == c) M[c * N + r] = rhs[i];
			else M[c*N + r] = M[r * N + c] = rhs[i];
			r++;
			if(r == N) r = ++c;
		}
	}

	CUDA_CALLABLE T& operator()(int row, int col)
	{
		return M[col * N + row];
	}

	CUDA_CALLABLE matrix_sq<T, N> operator=(T rhs)
	{
		int Nsq = N*N;
		for(int i=0; i<Nsq; i++)
			M[i] = rhs;

		return *this;
	}
	
	// M - rhs*I
	CUDA_CALLABLE matrix_sq<T, N> operator-(T rhs)
	{
		for(int i=0; i<N; i++)
			for(int j=0 ; j<N; j++)
				if(i == j)
					M[i*N+j] -= rhs;
		return *this;
	}
	
	// element wise divide
	CUDA_CALLABLE matrix_sq<T, N> operator/(T rhs)
	{
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++)
				M[i*N+j] = M[i*N+j]/rhs;
		return *this;
	}

	template<typename Y>
	CUDA_CALLABLE vec3<Y> operator*(vec3<Y> rhs){
		vec3<Y> result(0, 0, 0);
		for(int r=0; r<3; r++)
			for(int c=0; c<3; c++)
				result[r] += (*this)(r, c) * rhs[c];

		return result;
	}

	std::string toStr()
	{
		std::stringstream ss;

		for(int r = 0; r < N; r++)
		{
			ss << "| ";
			for(int c=0; c<N; c++)
			{
				ss << (*this)(r, c) << " ";
			}
			ss << "|" << std::endl;
		}

		return ss.str();
	}

        std::string toTensor()
        {
                std::stringstream ss; 

                ss << (*this)(0, 0) << " ";
                ss << (*this)(0, 1) << " ";
                ss << (*this)(1, 1) << " ";
                ss << (*this)(2, 0) << " ";
                ss << (*this)(2, 1) << " ";
                ss << (*this)(2, 2); 

                return ss.str();
        }


	static matrix_sq<T, N> identity() {
		matrix_sq<T, N> I;
		I = 0;
		for (size_t i = 0; i < N; i++)
			I.M[i * N + i] = 1;
		return I;
	}
};

}	//end namespace rts

template <typename T, int N>
std::ostream& operator<<(std::ostream& os, tira::matrix_sq<T, N> M)
{
    os<<M.toStr();
    return os;
}

//#if __GNUC__ > 3 && __GNUC_MINOR__ > 7
//template<class T, int N> using rtsMatrix = rts::matrix<T, N>;
//#endif

#endif
