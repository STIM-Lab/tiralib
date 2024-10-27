#pragma once

#include <cmath>

#include "geometry.h"

namespace tira{

	template <typename T>
	class cube : public geometry<T>{
	protected:
		void genGeometry(){
			T v[] = {
				-0.5, -0.5, -0.5,
				-0.5, 0.5, -0.5,
				0.5, 0.5, -0.5,
				0.5, -0.5, -0.5,
				-0.5, -0.5, 0.5,
				-0.5, 0.5, 0.5,
				0.5, 0.5, 0.5,
				0.5, -0.5, 0.5,
				0.5, -0.5, -0.5,
				0.5, -0.5, 0.5,
				0.5, 0.5, 0.5,
				0.5, 0.5, -0.5,
				-0.5, -0.5, -0.5,
				-0.5, -0.5, 0.5,
				-0.5, 0.5, 0.5,
				-0.5, 0.5, -0.5,
				-0.5, -0.5, -0.5,
				-0.5, -0.5, 0.5,
				0.5, -0.5, 0.5,
				0.5, -0.5, -0.5,
				-0.5, 0.5, -0.5,
				-0.5, 0.5, 0.5,
				0.5, 0.5, 0.5,
				0.5, 0.5, -0.5
			};
			geometry<T>::_vertices.insert(geometry<T>::_vertices.begin(), std::begin(v), std::end(v));

			T n[] = {
				0.0, 0.0, -1.0,
				0.0, 0.0, -1.0,
				0.0, 0.0, -1.0,
				0.0, 0.0, -1.0,
				0.0, 0.0, 1.0,
				0.0, 0.0, 1.0,
				0.0, 0.0, 1.0,
				0.0, 0.0, 1.0,
				1.0, 0.0, 0.0,
				1.0, 0.0, 0.0,
				1.0, 0.0, 0.0,
				1.0, 0.0, 0.0,
				-1.0, 0.0, 0.0,
				-1.0, 0.0, 0.0,
				-1.0, 0.0, 0.0,
				-1.0, 0.0, 0.0,
				0.0, -1.0, 0.0,
				0.0, -1.0, 0.0,
				0.0, -1.0, 0.0,
				0.0, -1.0, 0.0,
				0.0, 1.0, 0.0,
				0.0, 1.0, 0.0,
				0.0, 1.0, 0.0,
				0.0, 1.0, 0.0
			};
			geometry<T>::_normals.insert(geometry<T>::_normals.begin(), std::begin(n), std::end(n));

			T t[] = {
				// negative z plane
				0.0, 0.0,
				0.0, 1.0,
				1.0, 1.0,
				1.0, 0.0,
				//positive z plane
				0.0, 0.0,
				0.0, 1.0,
				1.0, 1.0,
				1.0, 0.0,
				//positive x plane
				0.0, 0.0,
				0.0, 1.0,
				1.0, 1.0,
				1.0, 0.0,
				//negative x plane
				0.0, 0.0,
				0.0, 1.0,
				1.0, 1.0,
				1.0, 0.0,
				//negative y plane
				0.0, 0.0,
				0.0, 1.0,
				1.0, 1.0,
				1.0, 0.0,
				//positive y plane
				0.0, 0.0,
				0.0, 1.0,
				1.0, 1.0,
				1.0, 0.0
			};
			geometry<T>::_texcoords.insert(geometry<T>::_texcoords.begin(), std::begin(t), std::end(t));

			int i[] = {
				0, 1, 2,
				0, 2, 3,

				4, 6, 5,
				4, 7, 6,

				8, 10, 9,
				8, 11, 10,

				12, 13, 14,
				12, 14, 15,

				16, 18, 17,
				16, 19, 18,

				20, 21, 22,
				20, 22, 23
			};
			geometry<T>::_indices.insert(geometry<T>::_indices.begin(), std::begin(i), std::end(i));
		}

		public:
			cube() {
				genGeometry();
			}
	};
}