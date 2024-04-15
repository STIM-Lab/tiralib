#pragma once

#include <cmath>

#include "geometry.h"

namespace tira {
	template <typename T>
	class sphere : public geometry<T>{
		void genGeometry(unsigned int stacks, unsigned int sectors) {

			float radius = 0.5f;
			float x, y, z, xy;                              // vertex position
			float nx, ny, nz, lengthInv = 1.0f / radius;    // vertex normal
			float s, t;                                     // vertex texCoord

			float sectorStep = 2 * 3.14159265358979323846 / sectors;
			float stackStep = 3.14159265358979323846 / stacks;
			float sectorAngle, stackAngle;

			// add the first vertex (top of the sphere)
			geometry<T>::m_vertices.push_back(0);
			geometry<T>::m_vertices.push_back(0);
			geometry<T>::m_vertices.push_back(radius);

			geometry<T>::m_normals.push_back(0);
			geometry<T>::m_normals.push_back(0);
			geometry<T>::m_normals.push_back(1);

			geometry<T>::m_texcoords.push_back(0);
			geometry<T>::m_texcoords.push_back(0);

			for (int i = 1; i < stacks; ++i)
			{
				stackAngle = 3.14159265358979323846 / 2 - i * stackStep;      // starting from pi/2 to -pi/2
				xy = radius * cosf(stackAngle);             // r * cos(u)
				z = radius * sinf(stackAngle);              // r * sin(u)

				// add (sectorCount+1) vertices per stack
				// the first and last vertices have same position and normal, but different tex coords
				for (int j = 0; j < sectors; ++j)
				{
					sectorAngle = j * sectorStep;           // starting from 0 to 2pi

					// vertex position (x, y, z)
					x = xy * cosf(sectorAngle);             // r * cos(u) * cos(v)
					y = xy * sinf(sectorAngle);             // r * cos(u) * sin(v)
					geometry<T>::m_vertices.push_back(x);
					geometry<T>::m_vertices.push_back(y);
					geometry<T>::m_vertices.push_back(z);

					// normalized vertex normal (nx, ny, nz)
					nx = x * lengthInv;
					ny = y * lengthInv;
					nz = z * lengthInv;
					geometry<T>::m_normals.push_back(nx);
					geometry<T>::m_normals.push_back(ny);
					geometry<T>::m_normals.push_back(nz);

					// vertex tex coord (s, t) range between [0, 1]
					s = (float)j / sectors;
					t = (float)i / stacks;
					geometry<T>::m_texcoords.push_back(s);
					geometry<T>::m_texcoords.push_back(t);
				}
			}

			// add the first vertex (top of the sphere)
			geometry<T>::m_vertices.push_back(0);
			geometry<T>::m_vertices.push_back(0);
			geometry<T>::m_vertices.push_back(-radius);

			geometry<T>::m_normals.push_back(0);
			geometry<T>::m_normals.push_back(0);
			geometry<T>::m_normals.push_back(-1);

			geometry<T>::m_texcoords.push_back(1);
			geometry<T>::m_texcoords.push_back(1);

			// generate CCW index list of sphere triangles
			// k1--k1+1
			// |  / |
			// | /  |
			// k2--k2+1
			//idx.clear();
			//std::vector<int> lineIndices;
			int k1, k2;
			for (int i = 0; i < stacks; ++i)
			{
				if (i == 0)
					k1 = 0;
				else
					k1 = (i - 1) * sectors + 1;     // beginning of current stack

				k2 = k1 + sectors;      // beginning of next stack

				for (int j = 0; j < sectors; ++j, ++k1, ++k2)
				{
					// 2 triangles per sector excluding first and last stacks
					// k1 => k2 => k1+1
					if (i == 0) {
						geometry<T>::m_indices.push_back(0);
						geometry<T>::m_indices.push_back(j + 1);

						if (j == sectors - 1)
							geometry<T>::m_indices.push_back(1);
						else
							geometry<T>::m_indices.push_back(j + 2);
					}
					else if (i == (stacks - 1)) {
						int last = geometry<T>::m_vertices.size() / 3 - 1;
						geometry<T>::m_indices.push_back(k1);
						geometry<T>::m_indices.push_back(last);
						if (j == sectors - 1)
							geometry<T>::m_indices.push_back(last - sectors);
						else
							geometry<T>::m_indices.push_back(k1 + 1);
					}
					else
					{
						if (j == sectors - 1) {
							geometry<T>::m_indices.push_back(k1);
							geometry<T>::m_indices.push_back(k2);
							geometry<T>::m_indices.push_back(k1 - sectors + 1);
							geometry<T>::m_indices.push_back(k1 - sectors + 1);
							geometry<T>::m_indices.push_back(k2);
							geometry<T>::m_indices.push_back(k2 - sectors + 1);
						}
						else {
							geometry<T>::m_indices.push_back(k1);
							geometry<T>::m_indices.push_back(k2);
							geometry<T>::m_indices.push_back(k1 + 1);
							geometry<T>::m_indices.push_back(k1 + 1);
							geometry<T>::m_indices.push_back(k2);
							geometry<T>::m_indices.push_back(k2 + 1);
						}
					}

					// k1+1 => k2 => k2+1
					//if (i != (stacks - 1))
					//{
					//	geometry<T>::m_indices.push_back(k1 + 1);
					//	geometry<T>::m_indices.push_back(k2);
					//	geometry<T>::m_indices.push_back(k2 + 1);
					//}

				}
			}
		} // end GenerateSphere()

	public:
		sphere(unsigned int stacks = 10, unsigned int slices = 10) {
			genGeometry(stacks, slices);
		}
	};
}