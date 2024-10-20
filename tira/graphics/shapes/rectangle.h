#pragma once

#include "geometry.h"

namespace tira{

	template <typename T>
	class rectangle : public geometry<T>{
	protected:
		void genGeometry(){
			T v[] = {
				-0.5, -0.5, 0.0,
				-0.5, 0.5, 0.0,
				0.5, 0.5, 0.0,
				0.5, -0.5, 0.0
			};
			geometry<T>::_vertices.insert(geometry<T>::_vertices.begin(), std::begin(v), std::end(v));

			T n[] = {
				0.0, 0.0, -1.0,
				0.0, 0.0, -1.0,
				0.0, 0.0, -1.0,
				0.0, 0.0, -1.0
			};
			geometry<T>::_normals.insert(geometry<T>::_normals.begin(), std::begin(n), std::end(n));

			T t[] = {
				0.0, 0.0,
				0.0, 1.0,
				1.0, 1.0,
				1.0, 0.0
			};
			geometry<T>::_texcoords.insert(geometry<T>::_texcoords.begin(), std::begin(t), std::end(t));

			int i[] = {
				0, 1, 2,
				0, 2, 3
			};
			geometry<T>::_indices.insert(geometry<T>::_indices.begin(), std::begin(i), std::end(i));
		}

		public:
			rectangle() {
				genGeometry();
			}
	};
}