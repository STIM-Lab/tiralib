#pragma once

#include <cmath>

#include "geometry.h"

namespace tira {
    template <typename T>
    class cylinder : public geometry<T> {
		// Function for generating cylinder geometry
        void genGeometry(unsigned int stacks, unsigned int sectors, T radius, T height) {

			T dz = height / stacks;						// change along the z axis
			T dtheta = 2 * M_PI / sectors;				// change in theta

			T zstart = - height / 2.0f;					// start position (bottom end cap)

			// CREATE THE BOTTOM END CAP
			geometry<T>::m_vertices.push_back(0);		// add the center vertex of the bottom end cap
			geometry<T>::m_vertices.push_back(0);
			geometry<T>::m_vertices.push_back(zstart);
			geometry<T>::m_normals.push_back(0);
			geometry<T>::m_normals.push_back(0);
			geometry<T>::m_normals.push_back(-1);
			geometry<T>::m_texcoords.push_back(0);
			geometry<T>::m_texcoords.push_back(0);

			for (unsigned int ti = 0; ti < sectors; ti++) {		// generate vertices for the triangle fan
				T theta = dtheta * ti;
				T x = radius * cos(theta);
				T y = radius * sin(theta);
				geometry<T>::m_vertices.push_back(x);
				geometry<T>::m_vertices.push_back(y);
				geometry<T>::m_vertices.push_back(zstart);

				geometry<T>::m_normals.push_back(0);
				geometry<T>::m_normals.push_back(0);
				geometry<T>::m_normals.push_back(-1);

				geometry<T>::m_texcoords.push_back((double)ti / (double)sectors);
				geometry<T>::m_texcoords.push_back(0);
			}
			for (unsigned int ti = 0; ti < sectors; ti++) {
				geometry<T>::m_indices.push_back(ti + 1);
				geometry<T>::m_indices.push_back(0);				
				if (ti == sectors - 1) geometry<T>::m_indices.push_back(1);
				else geometry<T>::m_indices.push_back(ti + 2);
			}

			// CREATE THE TUBE BODY

			for (unsigned int si = 1; si <= stacks; si++) {
				T z = zstart + si * dz;								// calculate the z coordinate of the current stack
				for (unsigned int ti = 0; ti < sectors; ti++) {
					T theta = dtheta * ti;
					T cos_theta = cos(theta);
					T sin_theta = sin(theta);
					T x = radius * cos_theta;
					T y = radius * sin_theta;

					geometry<T>::m_vertices.push_back(x);
					geometry<T>::m_vertices.push_back(y);
					geometry<T>::m_vertices.push_back(z);

					geometry<T>::m_normals.push_back(cos_theta);
					geometry<T>::m_normals.push_back(sin_theta);
					geometry<T>::m_normals.push_back(0);

					geometry<T>::m_texcoords.push_back((double)ti / (double)sectors);
					geometry<T>::m_texcoords.push_back((si * dz) / height);
				}
			}
			for (unsigned int si = 0; si < stacks; si++) {
				unsigned int s1_start = si * sectors + 1;
				unsigned int s2_start = s1_start + sectors;
				for (unsigned int ti = 0; ti < sectors; ti++) {
					geometry<T>::m_indices.push_back(s1_start + ti);
					if (ti == sectors - 1) geometry<T>::m_indices.push_back(s1_start);
					else geometry<T>::m_indices.push_back(s1_start + ti + 1);					
					geometry<T>::m_indices.push_back(s2_start + ti);

					geometry<T>::m_indices.push_back(s2_start + ti);
					if (ti == sectors - 1) geometry<T>::m_indices.push_back(s1_start);
					else geometry<T>::m_indices.push_back(s1_start + ti + 1);
					if (ti == sectors - 1) geometry<T>::m_indices.push_back(s2_start);
					else geometry<T>::m_indices.push_back(s2_start + ti + 1);
				}
			}

			// CREATE TOP END CAP
			geometry<T>::m_vertices.push_back(0);		// add the center vertex of the bottom end cap
			geometry<T>::m_vertices.push_back(0);
			geometry<T>::m_vertices.push_back(-zstart);
			geometry<T>::m_normals.push_back(0);
			geometry<T>::m_normals.push_back(0);
			geometry<T>::m_normals.push_back(1);
			geometry<T>::m_texcoords.push_back(1);
			geometry<T>::m_texcoords.push_back(0);

			unsigned int si_start = stacks * sectors + 1;
			unsigned int endcap = geometry<T>::m_vertices.size() / 3 - 1;
			for (unsigned int ti = 0; ti < sectors; ti++) {
				geometry<T>::m_indices.push_back(endcap);
				geometry<T>::m_indices.push_back(si_start + ti);
				if (ti == sectors - 1) geometry<T>::m_indices.push_back(si_start);
				else geometry<T>::m_indices.push_back(si_start + ti + 1);

			}

            /*//const float PI = 3.14159;

			// Declare variables for vertex position, normal, and texture coordinate
            float x, y, z, xy;											// vertex position
            float nx, ny, nz, lengthInv = 1.0f / radius;				// vertex normal
            float s, t;													// vertex texCoord

			// Calculate sector and stack steps based on the number of sectors and stacks
			float sectorStep = 2 * M_PI / sectors;
			float stackStep = height / stacks;
			float sectorAngle, stackHeight;


			// Loop through the number of stacks and sectors to generate vertices
			for (int i = 0; i <= stacks; ++i)
			{
				stackHeight = i * stackStep;        // starting from 0 to height

				// Add (sectorCount+1) vertices per stack
				// The first and last vertices have same position and normal, but different tex coords
				for (int j = 0; j <= sectors; ++j)
				{
					// starting from 0 to 2pi
					sectorAngle = j * sectorStep;				

					// vertex position (x, y, z)
					x = radius * cosf(sectorAngle);
					y = radius * sinf(sectorAngle);
					z = stackHeight;

					// Add vertex position to the vertex list
					geometry<T>::m_vertices.push_back(x);
					geometry<T>::m_vertices.push_back(y);
					geometry<T>::m_vertices.push_back(z);

					// Calculate normalized vertex normal (nx, ny, nz)
					nx = x / radius;
					ny = y / radius;
					nz = 0.0f;

					// Add vertex normal to the normal list
					geometry<T>::m_normals.push_back(nx);
					geometry<T>::m_normals.push_back(ny);
					geometry<T>::m_normals.push_back(nz);

					// Calculate vertex tex coord (s, t) range between [0, 1]
					s = (float)j / sectors;
					t = (float)i / stacks;

					// Add vertex texture coordinate to the texture coordinate list
					geometry<T>::m_texcoords.push_back(s);
					geometry<T>::m_texcoords.push_back(t);
				}
			}
			// generate CCW index list of cylinder triangles
			// k1--k1+1
			// |  / |
			// | /  |
			// k2--k2+1

			int k1, k2;
			for (int i = 0; i < stacks; ++i)
			{
				
				for (int j = 0; j < sectors; ++j, ++k1, ++k2)
				{
					k1 = i * (sectors + 1) + j;			// beginning of current stack
					k2 = k1 + sectors + 1;				// beginning of next stack


					// Add indices for two triangles that make up a side face of the cylinder
					geometry<T>::m_indices.push_back(k1);			// first vertex of first triangle
					geometry<T>::m_indices.push_back(k2);			// first vertex of second triangle
					geometry<T>::m_indices.push_back(k1 + 1);		// second vertex of first triangle
					
					
					geometry<T>::m_indices.push_back(k1 + 1);		// second vertex of first triangle
					geometry<T>::m_indices.push_back(k2);			// first vertex of second triangle
					geometry<T>::m_indices.push_back(k2 + 1);		// second vertex of second triangle
					
					// Add triangles for bottom faces
					if (i == 0 && j < sectors - 1)
					{
						geometry<T>::m_indices.push_back(k1);
						geometry<T>::m_indices.push_back(k1 + 1);
						geometry<T>::m_indices.push_back(0);
					}
					else if (i == stacks - 1 && j < sectors - 1)
					{
						geometry<T>::m_indices.push_back(k2 + 1);
						geometry<T>::m_indices.push_back(k2);
						geometry<T>::m_indices.push_back((stacks + 1) * (sectors + 1) - 1);
					}
				}
			}*/
		} // end GenerateCylinder()

	public:
		cylinder(unsigned int stacks = 10, unsigned int slices = 10, float radius = 0.5f, float height = 1.0f) {
			genGeometry(stacks, slices, radius, height);
		}
	};



}