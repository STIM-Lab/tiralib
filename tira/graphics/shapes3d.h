#include <vector>
#include "tira/graphics/glGeometry.h"

namespace tira {

	/// Generate the vertices and index arrays for a cube, including normals
	template <typename T>
	void GenerateCube(std::vector<T> &v, std::vector<unsigned int> &i) {
		v.clear();			//resize the array to match the number of required vertex components
		T av[] = {
			// vertex (x, y, z), normal(x, y, z)

			// negative z plane
			-0.5, -0.5, -0.5, 0.0, 0.0, -1.0, 0.0, 0.0,
			-0.5, 0.5, -0.5, 0.0, 0.0, -1.0, 0.0, 1.0,
			0.5, 0.5, -0.5, 0.0, 0.0, -1.0, 1.0, 1.0,
			0.5, -0.5, -0.5, 0.0, 0.0, -1.0, 1.0, 0.0,
			//positive z plane
			-0.5, -0.5, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0,
			-0.5, 0.5, 0.5, 0.0, 0.0, 1.0, 0.0, 1.0,
			0.5, 0.5, 0.5, 0.0, 0.0, 1.0, 1.0, 1.0,
			0.5, -0.5, 0.5, 0.0, 0.0, 1.0, 1.0, 0.0,
			//positive x plane
			0.5, -0.5, -0.5, 1.0, 0.0, 0.0, 0.0, 0.0,
			0.5, -0.5, 0.5, 1.0, 0.0, 0.0, 0.0, 1.0,
			0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0,
			0.5, 0.5, -0.5, 1.0, 0.0, 0.0, 1.0, 0.0,
			//negative x plane
			-0.5, -0.5, -0.5, -1.0, 0.0, 0.0, 0.0, 0.0,
			-0.5, -0.5, 0.5, -1.0, 0.0, 0.0, 0.0, 1.0,
			-0.5, 0.5, 0.5, -1.0, 0.0, 0.0, 1.0, 1.0,
			-0.5, 0.5, -0.5, -1.0, 0.0, 0.0, 1.0, 0.0,
			//negative y plane
			-0.5, -0.5, -0.5, 0.0, -1.0, 0.0, 0.0, 0.0,
			-0.5, -0.5, 0.5, 0.0, -1.0, 0.0, 0.0, 1.0,
			0.5, -0.5, 0.5, 0.0, -1.0, 0.0, 1.0, 1.0,
			0.5, -0.5, -0.5, 0.0, -1.0, 0.0, 1.0, 0.0,
			//positive y plane
			-0.5, 0.5, -0.5, 0.0, 1.0, 0.0, 0.0, 0.0,
			-0.5, 0.5, 0.5, 0.0, 1.0, 0.0, 0.0, 1.0,
			0.5, 0.5, 0.5, 0.0, 1.0, 0.0, 1.0, 1.0,
			0.5, 0.5, -0.5, 0.0, 1.0, 0.0, 1.0, 0.0
		};
		v.insert(v.begin(), std::begin(av), std::end(av));

		i.clear();
		int ai[] = {
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
		i.insert(i.begin(), std::begin(ai), std::end(ai));
	}

	/// Generate the vertices and index arrays for a cube, including normals
	template <typename T>
	void GenerateSphere(std::vector<T>& v, std::vector<unsigned int>& idx, unsigned int stacks, unsigned int sectors) {
		v.clear();
		const float PI = 3.14159;
		//std::vector<float>().swap(vertices);
		//std::vector<float>().swap(normals);
		//std::vector<float>().swap(texCoords);

		float radius = 0.5f;
		float x, y, z, xy;                              // vertex position
		float nx, ny, nz, lengthInv = 1.0f / radius;    // vertex normal
		float s, t;                                     // vertex texCoord

		float sectorStep = 2 * PI / sectors;
		float stackStep = PI / stacks;
		float sectorAngle, stackAngle;

		for (int i = 0; i <= stacks; ++i)
		{
			stackAngle = PI / 2 - i * stackStep;        // starting from pi/2 to -pi/2
			xy = radius * cosf(stackAngle);             // r * cos(u)
			z = radius * sinf(stackAngle);              // r * sin(u)

			// add (sectorCount+1) vertices per stack
			// the first and last vertices have same position and normal, but different tex coords
			for (int j = 0; j <= sectors; ++j)
			{
				sectorAngle = j * sectorStep;           // starting from 0 to 2pi

				// vertex position (x, y, z)
				x = xy * cosf(sectorAngle);             // r * cos(u) * cos(v)
				y = xy * sinf(sectorAngle);             // r * cos(u) * sin(v)
				v.push_back(x);
				v.push_back(y);
				v.push_back(z);

				// normalized vertex normal (nx, ny, nz)
				nx = x * lengthInv;
				ny = y * lengthInv;
				nz = z * lengthInv;
				v.push_back(nx);
				v.push_back(ny);
				v.push_back(nz);

				// vertex tex coord (s, t) range between [0, 1]
				s = (float)j / sectors;
				t = (float)i / stacks;
				v.push_back(s);
				v.push_back(t);
			}
		}
		// generate CCW index list of sphere triangles
		// k1--k1+1
		// |  / |
		// | /  |
		// k2--k2+1
		idx.clear();
		//std::vector<int> lineIndices;
		int k1, k2;
		for (int i = 0; i < stacks; ++i)
		{
			k1 = i * (sectors + 1);     // beginning of current stack
			k2 = k1 + sectors + 1;      // beginning of next stack

			for (int j = 0; j < sectors; ++j, ++k1, ++k2)
			{
				// 2 triangles per sector excluding first and last stacks
				// k1 => k2 => k1+1
				if (i != 0)
				{
					idx.push_back(k1);
					idx.push_back(k2);
					idx.push_back(k1 + 1);
				}

				// k1+1 => k2 => k2+1
				if (i != (stacks - 1))
				{
					idx.push_back(k1 + 1);
					idx.push_back(k2);
					idx.push_back(k2 + 1);
				}

				/*// store indices for lines
				// vertical lines for all stacks, k1 => k2
				lineIndices.push_back(k1);
				lineIndices.push_back(k2);
				if (i != 0)  // horizontal lines except 1st stack, k1 => k+1
				{
					lineIndices.push_back(k1);
					lineIndices.push_back(k1 + 1);
				}*/
			}
		}
	}
	
}