#pragma once

#include <cmath>

#include "geometry.h"

namespace tira{

	template <typename T>
	class circle : public geometry<T>{
	protected:
		void genGeometry(size_t points){
			std::vector<T> v(points * 3 + 3);							// generate a vector to store 3D points representing the circle

			// calculate vertices
			float d_rad = 2 * M_PI / points;
			T theta, x, y;

			for(size_t p = 0; p < points; p++){				// for each point
				theta = d_rad * p;							// calculate the theta coordinate (angle around the circle) for the point
				x = 0.5 * cos(theta);
				y = 0.5 * sin(theta);
				v[p * 3 + 0] = x;
				v[p * 3 + 1] = y;
				v[p * 3 + 2] = 0;
			}
			v[points * 3 + 0] = 0.0;						// add the center point
			v[points * 3 + 1] = 0.0;
			v[points * 3 + 2] = 0.0;

			geometry<T>::m_vertices.insert(geometry<T>::m_vertices.begin(), std::begin(v), std::end(v));

			std::vector<T> n(points * 3 + 3);
			for(size_t p = 0; p < points; p++){
				n[p * 3 + 0] = 0;
				n[p * 3 + 1] = 0;
				n[p * 3 + 2] = 1;
			}
			n[points * 3 + 0] = 0.0;						// add the center normal
			n[points * 3 + 1] = 0.0;
			n[points * 3 + 2] = 1.0;
			geometry<T>::m_normals.insert(geometry<T>::m_normals.begin(), std::begin(n), std::end(n));
			std::vector<T> t(points * 2 + 2);
			for(size_t p = 0; p < points; p++){
				t[p * 2 + 0] = v[p * 3 + 0] + 0.5;
				t[p * 2 + 1] = v[p * 3 + 1] + 0.5;
			}
			t[points * 2 + 0] = 0.5;						// add the center texture coordinate
			t[points * 2 + 1] = 0.5;

			geometry<T>::m_texcoords.insert(geometry<T>::m_texcoords.begin(), std::begin(t), std::end(t));

			// calculate the indices
			std::vector<int> i(points * 3);
			for(size_t p = 0; p < points - 1; p++){
				i[p * 3 + 0] = (int)p;
				i[p * 3 + 1] = (int)points;
				i[p * 3 + 2] = (int)p + 1;
			}
			i[(points - 1) * 3 + 0] = (int)points - 1;
			i[(points - 1) * 3 + 1] = (int)points;
			i[(points - 1) * 3 + 2] = (int)0;
			geometry<T>::m_indices.insert(geometry<T>::m_indices.begin(), std::begin(i), std::end(i));
		}

		public:
			circle(size_t points = 10) {
				genGeometry(points);
			}
	};
}