#pragma once

#include "tmesh.h"

namespace tira{

	static tmesh rectangle(){
		tmesh R;
			std::vector<float> v = {
				-0.5, -0.5, 0.0,
				-0.5, 0.5, 0.0,
				0.5, 0.5, 0.0,
				0.5, -0.5, 0.0
			};
            R.Vertices(&v[0], v.size() / 3);

			std::vector<float> n = {
				0.0, 0.0, -1.0,
				0.0, 0.0, -1.0,
				0.0, 0.0, -1.0,
				0.0, 0.0, -1.0
			};
            R.Normals(&n[0], n.size() / 3);

			std::vector<float> t = {
				0.0, 0.0,
				0.0, 1.0,
				1.0, 1.0,
				1.0, 0.0
			};
            R.TexCoords(&t[0], n.size() / 2, 2);

			std::vector<unsigned int> i = {
				0, 1, 2,
				0, 2, 3
			};
            R.Indices(i);

            return R;
		}
}