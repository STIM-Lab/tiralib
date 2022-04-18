#pragma once
#include <vector>

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
		size_t getNumVertices() { return m_vertices.size() / 3; }
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
	};
}