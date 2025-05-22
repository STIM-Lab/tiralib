#pragma once
#include "icosphere.h"

namespace tira {

	template<typename Type>
	class Superquadric : public Icosphere<Type> {
	protected:


		Type signpow(Type x, Type exponent) {
			if (x < 0) return -pow(abs(x), exponent);
			else return pow(abs(x), exponent);
		}

		void sq_vertex(Type alpha, Type beta, Type theta, Type phi, Type& x, Type& y, Type& z) {

			Type cos_phi = cos(phi);
			Type sin_theta = sin(theta);
			Type sin_phi = sin(phi);
			Type cos_theta = cos(theta);

			x = signpow(cos_phi, beta);
			y = -signpow(sin_theta, alpha) * signpow(sin_phi, beta);
			z = signpow(cos_theta, alpha) * signpow(sin_phi, beta);
		}

	public:
		Superquadric(Type l0, Type l1, Type l2, Type gamma, Type radius = 1.0f, unsigned int subdiv = 4, bool smooth = false) : Icosphere<Type>(radius, subdiv, smooth) {
			
			// calculate the linear and planar anisotropy
			Type suml = l0 + l1 + l2;
			Type Cl = (l0 - l1) / suml;
			Type Cp = 2 * (l1 - l2) / suml;

			size_t Nv = Icosphere<Type>::_vertices.size() / 3;
			for (size_t vi = 0; vi < Nv; vi++) {
				Type x = Icosphere<Type>::_vertices[vi * 3 + 0];
				Type y = Icosphere<Type>::_vertices[vi * 3 + 1];
				Type z = Icosphere<Type>::_vertices[vi * 3 + 2];

				Type theta = atan2(y, x);
				Type phi = atan2(sqrt(x * x + y * y), z);

				Type nx, ny, nz;

				if (Cl >= Cp) {
					Type alpha = pow(1 - Cp, gamma);
					Type beta = pow(1 - Cl, gamma);
					sq_vertex(alpha, beta, theta, phi, x, y, z);
					sq_vertex(2.0f - alpha, 2.0f - beta, theta, phi, nx, ny, nz);
				}
				else {
					Type alpha = pow(1 - Cl, gamma);
					Type beta = pow(1 - Cp, gamma);
					sq_vertex(alpha, beta, theta, phi, z, y, x);
					y = -y;
					sq_vertex(2.0f - alpha, 2.0f - beta, theta, phi, nz, ny, nx);
					ny = -ny;					
				}

				Icosphere<Type>::_vertices[vi * 3 + 0] = l0 * x * radius;
				Icosphere<Type>::_vertices[vi * 3 + 1] = l1 * y * radius;
				Icosphere<Type>::_vertices[vi * 3 + 2] = l2 * z * radius;

				Type sx = 1.0f / l0;
				Type sy = 1.0f / l1;
				Type sz = 1.0f / l2;
				Icosphere<Type>::_normals[vi * 3 + 0] = sx * nx;
				Icosphere<Type>::_normals[vi * 3 + 1] = sy * ny;
				Icosphere<Type>::_normals[vi * 3 + 2] = sz * nz;

			}
			Icosphere<Type>::buildInterleavedVertices();
		}
	};

	template<typename Type>
	trimesh<Type> superquadric(Type l0, Type l1, Type l2, Type gamma, Type radius = 1.0f, unsigned int subdiv = 4, bool smooth = false) {

		trimesh<Type> S;

		Superquadric s(l0, l1, l2, gamma, radius, subdiv, smooth);
		S.vertices(s._vertices);
		S.normals(s._normals);
		S.texcoords(s._texcoords);
		S.indices(s._indices);
		return S;
	}
}