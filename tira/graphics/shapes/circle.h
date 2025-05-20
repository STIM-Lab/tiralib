#pragma once

#include <cmath>

#include "mesh_triangle.h"

namespace tira{

	template <typename T>
	trimesh<T> circle(size_t points){
        trimesh<T> S;
		std::vector<T> v(points * 3 + 3);							// generate a vector to store 3D points representing the circle

		// calculate vertices
		float d_rad = 2 * 3.14159265358979323846 / points;
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

		S.vertices(v);

		std::vector<T> n(points * 3 + 3);
		for(size_t p = 0; p < points; p++){
			n[p * 3 + 0] = 0;
			n[p * 3 + 1] = 0;
			n[p * 3 + 2] = 1;
		}
		n[points * 3 + 0] = 0.0;						// add the center normal
		n[points * 3 + 1] = 0.0;
		n[points * 3 + 2] = 1.0;
        S.normals(n);

		// texture coordinates will be polar coordinates
		std::vector<T> t(points * 2 + 2);
		for(size_t p = 0; p < points; p++){
			t[p * 2 + 0] = 0.5;
			t[p * 2 + 1] = d_rad * p;
		}
		t[points * 2 + 0] = 0.0;						// add the center texture coordinate
		t[points * 2 + 1] = 0.0;
        S.texcoords(t);


		// calculate the indices
		std::vector<unsigned int> i(points * 3);
		for(size_t p = 0; p < points - 1; p++){
			i[p * 3 + 0] = (unsigned int)p;
			i[p * 3 + 1] = (unsigned int)points;
			i[p * 3 + 2] = (unsigned int)p + 1;
		}
		i[(points - 1) * 3 + 0] = (unsigned int)points - 1;
		i[(points - 1) * 3 + 1] = (unsigned int)points;
		i[(points - 1) * 3 + 2] = (unsigned int)0;
        S.indices(i);
        return S;
	}
}