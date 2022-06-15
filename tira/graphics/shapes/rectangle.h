#pragma once

#include <vector>
#include "shape.h"

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
			geometry<T>::m_vertices.insert(geometry<T>::m_vertices.begin(), std::begin(v), std::end(v));

			T n[] = {
				0.0, 0.0, -1.0,
				0.0, 0.0, -1.0,
				0.0, 0.0, -1.0,
				0.0, 0.0, -1.0
			};
			geometry<T>::m_normals.insert(geometry<T>::m_normals.begin(), std::begin(n), std::end(n));

			T t[] = {
				0.0, 0.0,
				0.0, 1.0,
				1.0, 1.0,
				1.0, 0.0
			};
			geometry<T>::m_texcoords.insert(geometry<T>::m_texcoords.begin(), std::begin(t), std::end(t));

			int i[] = {
				0, 1, 2,
				0, 2, 3
			};
			geometry<T>::m_indices.insert(geometry<T>::m_indices.begin(), std::begin(i), std::end(i));
		}

		public:
			rectangle() {
				genGeometry();
			}
	};
}