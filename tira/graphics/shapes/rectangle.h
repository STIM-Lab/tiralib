#pragma once

#include "mesh_triangle.h"

namespace tira{

	template <typename T>
	trimesh<T> rectangle(){
          	trimesh<T> R;
			std::vector<T> v = {
				-0.5, -0.5, 0.0,
				-0.5, 0.5, 0.0,
				0.5, 0.5, 0.0,
				0.5, -0.5, 0.0
			};
            R.vertices(v);

			std::vector<T> n = {
				0.0, 0.0, -1.0,
				0.0, 0.0, -1.0,
				0.0, 0.0, -1.0,
				0.0, 0.0, -1.0
			};
            R.normals(n);

			std::vector<T> t = {
				0.0, 0.0,
				0.0, 1.0,
				1.0, 1.0,
				1.0, 0.0
			};
            R.texcoords(t);

			std::vector<unsigned int> i = {
				0, 1, 2,
				0, 2, 3
			};
            R.indices(i);

            return R;
		}
}