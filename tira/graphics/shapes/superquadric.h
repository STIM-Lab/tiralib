#include "icosphere.h"

namespace tira {
	class superquadric : public icosphere {
	protected:


		float signpow(float x, float exponent) {
			if (x < 0) return -pow(abs(x), exponent);
			else return pow(abs(x), exponent);
		}

		void sq_vertex(float alpha, float beta, float theta, float phi, float& x, float& y, float& z) {

			float cos_phi = cos(phi);
			float sin_theta = sin(theta);
			float sin_phi = sin(phi);
			float cos_theta = cos(theta);

			x = signpow(cos_phi, beta);
			y = -signpow(sin_theta, alpha) * signpow(sin_phi, beta);
			z = signpow(cos_theta, alpha) * signpow(sin_phi, beta);
		}

	public:
		superquadric(float l0, float l1, float l2, float gamma, float radius = 1.0f, unsigned int subdiv = 4, bool smooth = false) : icosphere(radius, subdiv, smooth) {
			
			// calculate the linear and planar anisotropy
			float suml = l0 + l1 + l2;
			float Cl = (l0 - l1) / suml;
			float Cp = 2 * (l1 - l2) / suml;

			size_t Nv = vertices.size() / 3;
			for (size_t vi = 0; vi < Nv; vi++) {
				float x = vertices[vi * 3 + 0];
				float y = vertices[vi * 3 + 1];
				float z = vertices[vi * 3 + 2];

				float theta = atan2(y, x);
				float phi = atan2(sqrt(x * x + y * y), z);

				float nx, ny, nz;

				if (Cl >= Cp) {
					float alpha = pow(1 - Cp, gamma);
					float beta = pow(1 - Cl, gamma);
					sq_vertex(alpha, beta, theta, phi, x, y, z);
					sq_vertex(2.0f - alpha, 2.0f - beta, theta, phi, nx, ny, nz);
				}
				else {
					float alpha = pow(1 - Cl, gamma);
					float beta = pow(1 - Cp, gamma);
					sq_vertex(alpha, beta, theta, phi, z, y, x);
					y = -y;
					sq_vertex(2.0f - alpha, 2.0f - beta, theta, phi, nz, ny, nx);
					ny = -ny;					
				}

				vertices[vi * 3 + 0] = l0 * x * radius;
				vertices[vi * 3 + 1] = l1 * y * radius;
				vertices[vi * 3 + 2] = l2 * z * radius;

				float sx = 1.0f / l0;
				float sy = 1.0f / l1;
				float sz = 1.0f / l2;
				normals[vi * 3 + 0] = sx * nx;
				normals[vi * 3 + 1] = sy * ny;
				normals[vi * 3 + 2] = sz * nz;

			}
			buildInterleavedVertices();
		}
	};
}