#pragma once

#include <cmath>

#include "tmesh.h"

namespace tira {
	static tmesh sphere(unsigned int stacks, unsigned int sectors){

			float radius = 0.5f;
			float x, y, z, xy;                              // vertex position
			float nx, ny, nz, lengthInv = 1.0f / radius;    // vertex normal
			float s, t;                                     // vertex texCoord

			float sectorStep = 2 * 3.14159265358979323846 / sectors;
			float stackStep = 3.14159265358979323846 / stacks;
			float sectorAngle, stackAngle;

            std::vector<float> vertices;
			// add the first vertex (top of the sphere)
			vertices.push_back(0);
			vertices.push_back(0);
			vertices.push_back(radius);

            std::vector<float> normals;
			normals.push_back(0);
			normals.push_back(0);
			normals.push_back(1);

            std::vector<float> texcoords;
			texcoords.push_back(0);
			texcoords.push_back(0);

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
					vertices.push_back(x);
					vertices.push_back(y);
					vertices.push_back(z);

					// normalized vertex normal (nx, ny, nz)
					nx = x * lengthInv;
					ny = y * lengthInv;
					nz = z * lengthInv;
					normals.push_back(nx);
					normals.push_back(ny);
					normals.push_back(nz);

					// vertex tex coord (s, t) range between [0, 1]
					s = (float)j / sectors;
					t = (float)i / stacks;
					texcoords.push_back(s);
					texcoords.push_back(t);
				}
			}

			// add the first vertex (top of the sphere)
			vertices.push_back(0);
			vertices.push_back(0);
			vertices.push_back(-radius);

			normals.push_back(0);
			normals.push_back(0);
			normals.push_back(-1);

			texcoords.push_back(1);
			texcoords.push_back(1);

			// generate CCW index list of sphere triangles
			// k1--k1+1
			// |  / |
			// | /  |
			// k2--k2+1
			//idx.clear();
			//std::vector<int> lineIndices;
            std::vector<unsigned int> indices;
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
						indices.push_back(0);
						indices.push_back(j + 1);

						if (j == sectors - 1)
							indices.push_back(1);
						else
							indices.push_back(j + 2);
					}
					else if (i == (stacks - 1)) {
						int last = vertices.size() / 3 - 1;
						indices.push_back(k1);
						indices.push_back(last);
						if (j == sectors - 1)
							indices.push_back(last - sectors);
						else
							indices.push_back(k1 + 1);
					}
					else
					{
						if (j == sectors - 1) {
							indices.push_back(k1);
							indices.push_back(k2);
							indices.push_back(k1 - sectors + 1);
							indices.push_back(k1 - sectors + 1);
							indices.push_back(k2);
							indices.push_back(k2 - sectors + 1);
						}
						else {
							indices.push_back(k1);
							indices.push_back(k2);
							indices.push_back(k1 + 1);
							indices.push_back(k1 + 1);
							indices.push_back(k2);
							indices.push_back(k2 + 1);
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
			tmesh S;
            S.Vertices(&vertices[0], vertices.size() / 3);
            S.Normals(&normals[0], normals.size() / 3);
            S.TexCoords(&texcoords[0], texcoords.size()/2, 2);
            S.Indices(indices);
            return S;
		} // end GenerateSphere()
}