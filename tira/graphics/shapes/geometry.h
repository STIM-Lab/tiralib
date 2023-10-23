#pragma once
#include <vector>
#include <string>
#include <fstream>

namespace tira{
	template <typename T>
	class geometry {
	protected:
		std::vector<T> m_vertices;
		std::vector<T> m_normals;
		std::vector<T> m_texcoords;
		std::vector<unsigned int> m_indices;
		size_t vertex_dim = 3;
		size_t normal_dim = 3;
		size_t texture_dim = 2;

	public:
		std::vector<T> getVertices() { return m_vertices; }
		std::vector<T> getNormals() { return m_normals; }
		std::vector<T> getTexCoords() { return m_texcoords; }
		std::vector<unsigned int> getIndices() { return m_indices; }
		size_t getNumVertices() { return m_vertices.size() / vertex_dim; }
		size_t getNumTriangles() { return m_indices.size() / 3; }
		std::vector<T> getInterleavedVertices() {
			// create a vector to hold the interleaved vertices
			std::vector<T> interleaved(m_vertices.size() +
				m_normals.size() +
				m_texcoords.size());


			size_t stride = vertex_dim + normal_dim + texture_dim;
			for (size_t i = 0; i < getNumVertices(); i++) {
				for (size_t v = 0; v < vertex_dim; v++)
					interleaved[i * stride + v] = m_vertices[i * vertex_dim + v];
				for (size_t n = 0; n < normal_dim; n++)
					interleaved[i * stride + vertex_dim + n] = m_normals[i * normal_dim + n];
				for (size_t t = 0; t < texture_dim; t++)
					interleaved[i * stride + vertex_dim + normal_dim + t] = m_texcoords[i * texture_dim + t];
			}

			return interleaved;
		}

		// save an OBJ file for this geometry
		void obj(std::string filename) {
			std::ofstream f(filename);										// open a file for writing
			size_t nV = getNumVertices();									// get the number of vertices
			for (size_t vi = 0; vi < nV; vi++) {							// for each vertex
				f << "v ";													// output the initial "v" to indicate the vertex
				for (size_t ci = 0; ci < vertex_dim; ci++) {				// for each coordinate
					f << m_vertices[3 * vi + ci];							// output the coordinate value
					if (ci < vertex_dim - 1) f << " ";						// output a space (if this isn't the last coordinate)
				}
				if (vi < nV - 1)f << std::endl;								// move to the next line
			}

			if (m_texcoords.size()) {										// if there are any texture coordinates
				f << std::endl;
				for (size_t ti = 0; ti < nV; ti++) {						// for each vertex
					f << "vt ";												// output the initial "v" to indicate the vertex
					for (size_t ci = 0; ci < texture_dim; ci++) {			// for each coordinate
						f << m_texcoords[texture_dim * ti + ci];						// output the coordinate value
						if (ci < texture_dim - 1) f << " ";					// output a space (if this isn't the last coordinate)
					}
					if (ti < nV - 1)f << std::endl;							// move to the next line
				}
			}

			if (m_normals.size()) {										// if there are any texture coordinates
				f << std::endl;
				for (size_t ni = 0; ni < nV; ni++) {						// for each vertex
					f << "vn ";												// output the initial "v" to indicate the vertex
					for (size_t ci = 0; ci < normal_dim; ci++) {			// for each coordinate
						f << m_normals[normal_dim * ni + ci];						// output the coordinate value
						if (ci < normal_dim - 1) f << " ";					// output a space (if this isn't the last coordinate)
					}
					if (ni < nV - 1)f << std::endl;							// move to the next line
				}
			}

			f << std::endl;
			size_t nT = getNumTriangles();								//get the number of triangles
			for (size_t ti = 0; ti <nT; ti++) {					// for each vertex
				f << "f ";												// output the initial "v" to indicate the vertex
				for (size_t i = 0; i < 3; i++) {						// for each point of the triangle
					f << m_indices[3 * ti + i] + 1;						// output the coordinate value
					if (i < 2) f << " ";					// output a space (if this isn't the last coordinate)
				}
				if (ti < nT - 1)f << std::endl;							// move to the next line
			}
			f.close();														// close the OBJ file
		}
	};
}