#pragma once

#include <cmath>

#include "geometry.h"

namespace tira {

	static tmesh cylinder(unsigned int stacks = 1, unsigned int sectors = 10, float radius = 0.5, float height = 1) {

		std::vector<float> vertices;
		std::vector<float> normals;
		std::vector<float> texcoords;
		std::vector<unsigned int> indices;

		float dz = height / stacks;						// change along the z axis
		float dtheta = 2 * 3.14159265358979323846 / sectors;				// change in theta

		float zstart = - height / 2.0f;					// start position (bottom end cap)

		// CREATE THE BOTTOM END CAP
		vertices.push_back(0);		// add the center vertex of the bottom end cap
		vertices.push_back(0);
		vertices.push_back(zstart);
		normals.push_back(0);
		normals.push_back(0);
		normals.push_back(-1);
		texcoords.push_back(0);
		texcoords.push_back(0);

		for (unsigned int ti = 0; ti < sectors; ti++) {		// generate vertices for the triangle fan
			float theta = dtheta * ti;
			float x = radius * cos(theta);
			float y = radius * sin(theta);
			vertices.push_back(x);
			vertices.push_back(y);
			vertices.push_back(zstart);

			normals.push_back(0);
			normals.push_back(0);
			normals.push_back(-1);

			texcoords.push_back((double)ti / (double)sectors);
			texcoords.push_back(0);
		}
		for (unsigned int ti = 0; ti < sectors; ti++) {
			indices.push_back(ti + 1);
			indices.push_back(0);				
			if (ti == sectors - 1) indices.push_back(1);
			else indices.push_back(ti + 2);
		}

		// CREATE THE TUBE BODY

		for (unsigned int si = 1; si <= stacks; si++) {
			float z = zstart + si * dz;								// calculate the z coordinate of the current stack
			for (unsigned int ti = 0; ti < sectors; ti++) {
				float theta = dtheta * ti;
				float cos_theta = cos(theta);
				float sin_theta = sin(theta);
				float x = radius * cos_theta;
				float y = radius * sin_theta;

				vertices.push_back(x);
				vertices.push_back(y);
				vertices.push_back(z);

				normals.push_back(cos_theta);
				normals.push_back(sin_theta);
				normals.push_back(0);

				texcoords.push_back((double)ti / (double)sectors);
				texcoords.push_back((si * dz) / height);
			}
		}
		for (unsigned int si = 0; si < stacks; si++) {
			unsigned int s1_start = si * sectors + 1;
			unsigned int s2_start = s1_start + sectors;
			for (unsigned int ti = 0; ti < sectors; ti++) {
				indices.push_back(s1_start + ti);
				if (ti == sectors - 1) indices.push_back(s1_start);
				else indices.push_back(s1_start + ti + 1);					
				indices.push_back(s2_start + ti);

				indices.push_back(s2_start + ti);
				if (ti == sectors - 1) indices.push_back(s1_start);
				else indices.push_back(s1_start + ti + 1);
				if (ti == sectors - 1) indices.push_back(s2_start);
				else indices.push_back(s2_start + ti + 1);
			}
		}

		// CREATE TOP END CAP
		vertices.push_back(0);		// add the center vertex of the bottom end cap
		vertices.push_back(0);
		vertices.push_back(-zstart);
		normals.push_back(0);
		normals.push_back(0);
		normals.push_back(1);
		texcoords.push_back(1);
		texcoords.push_back(0);

		unsigned int si_start = stacks * sectors + 1;
		unsigned int endcap = vertices.size() / 3 - 1;
		for (unsigned int ti = 0; ti < sectors; ti++) {
			indices.push_back(endcap);
			indices.push_back(si_start + ti);
			if (ti == sectors - 1) indices.push_back(si_start);
			else indices.push_back(si_start + ti + 1);

		}

		tmesh C;
		C.Vertices(&vertices[0], vertices.size() / 3);
		C.Normals(&normals[0], normals.size() / 3);
		C.TexCoords(&texcoords[0], texcoords.size() / 2, 2);
		C.Indices(indices);
		return C;
	} // end GenerateCylinder()

}