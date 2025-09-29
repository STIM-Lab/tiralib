 #pragma once
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <functional>
#include <exception>
#include <iostream>

#include <glm/glm.hpp>

/* The TIRA trimesh class is designed to store data for efficiently render triangular meshes
	using graphics APIs. It provides details for storing vertices, normals, and texture coordinates.
 */

namespace tira{

	/**
	 * @brief The TIRA tmesh class is designed to store data for efficiently rendering and transferring triangular
	 * 				meshes through a variety of APIs. It provides a framework for storing vertex positions, normals, and
	 * 				texture coordinates.
	 * 
	 */
	class tmesh {
	protected:

		/**
		 * @brief Stores the coordinates defining vertex positions, and assumed to be 3 coordinates for a 3D object.
		 */
		std::vector<glm::vec3> m_vertices;

		/**
		 * @brief Stores the coordinates defining vertex normals, and assumed to be 3 coordinates for a 3D object.
		 */
		std::vector<glm::vec3> m_normals;

		/**
		 * @brief Stores the coordinates definining vertex texture coordinates. Texture coordinates are generally stored as
		 * 			3D coordinates, however several helper functions allow this to be reduced when outputting the data.
		 * 
		 */
		std::vector<glm::vec3> m_texcoords;

		/**
		 * @brief Stores the indices defining the geometric connectivity through the 3D mesh. Every three vertices specify a single triangle.
		 * 
		 */
		std::vector<unsigned int> m_indices;

	public:
		// deprecated
		//std::vector<glm::vec3> getVertices() const { return m_vertices; }
		//std::vector<glm::vec3> getNormals() const { return m_normals; }
		//std::vector<glm::vec3> getTexCoords() const { return m_texcoords; }
		//std::vector<unsigned int> getIndices() const { return m_indices; }
		//unsigned int getVertexDim() const { return _vertex_dim; }
		//unsigned int getNormalDim() const { return _normal_dim; }
		//unsigned int getTextureDim() const { return _texture_dim; }

		/**
		 * @brief Returns an std::vector of 3D coordinates associated with vertex positions.
		 * 
		 * @return std::vector<glm::vec3> is an array of vec3 positions
		 */
		std::vector<glm::vec3> Vertices() const { return m_vertices; }

		/**
		 * @brief Returns the vertex positions as a flattened 1D array of floating point values. This makes the structure
		 * 			compatible with most APIs that do not work specifically with glm::vec3 structures.
		 * 
		 * @return std::vector<float> is the flattened array of positions as floating point values
		 */
		std::vector<float> RavelVertices() {
			std::vector<float> result(m_vertices.size() * sizeof(glm::vec3));
			memcpy(&result[0], &m_vertices[0], m_vertices.size() * sizeof(glm::vec3));
			return result;
		}

		/**
		 * @brief Sets the vertex positions for this mesh.
		 * 
		 * @param vertices is an array of vertex positions that will replace any existing vertex data.
		 */
		void Vertices(std::vector<glm::vec3>& vertices) { m_vertices = vertices; }

		/**
		 * @brief Sets the vertex positions for this mesh using a flattened floating point array.
		 * 
		 * @param v is a pointer to the vertex position data
		 * @param N is the number of vertices specified in this array
		 */
		void Vertices(float* v, size_t N) {
			m_vertices.resize(N);
			memcpy(&m_vertices[0], v, N * sizeof(glm::vec3));
		}

		/**
		 * @brief Returns an std::vector of 3D coordinates associated with vertex normals.
		 * 
		 * @return std::vector<glm::vec3> is an array of vec3 normals for each vertex
		 */
		std::vector<glm::vec3> Normals() const { return m_normals; }

		/**
		 * @brief Returns the vertex normals as a flattened 1D array of floating point values. This makes the structure
		 * 			compatible with most APIs that do not work specifically with glm::vec3 structures.
		 * 
		 * @return std::vector<float> is the flattened array of normals as floating point values
		 */
		std::vector<float> RavelNormals() {
			std::vector<float> result(m_normals.size() * sizeof(glm::vec3));
			memcpy(&result[0], &m_normals[0], m_normals.size() * sizeof(glm::vec3));
			return result;
		}

		/**
		 * @brief Sets the vertex normals for this mesh.
		 * 
		 * @param normals is an array of vertex normals that will replace any existing vertex data.
		 */
		void Normals(std::vector<glm::vec3>& normals) { m_normals = normals; }

		/**
		 * @brief Sets the vertex normals for this mesh using a flattened floating point array.
		 * 
		 * @param n is a pointer to the vertex normal data
		 * @param N is the number of vertices specified in this array
		 */
		void Normals(float* n, size_t N) {
			m_normals.resize(N);
			memcpy(&m_normals[0], n, N * sizeof(glm::vec3));
		}

		/**
		 * @brief Sets the vertex texture coordinates for this mesh.
		 * 
		 * @return std::vector<glm::vec3> is an array of vec3 texture coordinates for each vertex
		 */
		std::vector<glm::vec3> TexCoords() const { return m_texcoords; }

		/**
		 * @brief Returns the vertex texture coordinates as a flattened 1D array of floating point values. This makes the structure
		 * 			compatible with most APIs that do not work specifically with glm::vec3 structures.
		 * 
		 * @return std::vector<float> is the flattened array of texture coordinates as floating point values
		 */
		std::vector<float> RavelTexCoords() {
			std::vector<float> result(m_texcoords.size() * sizeof(glm::vec3));
			memcpy(&result[0], &m_texcoords[0], m_texcoords.size() * sizeof(glm::vec3));
			return result;
		}

		/**
		 * @brief Sets the vertex texture coordinates for this mesh.
		 * 
		 * @param texcoords is an array of vertex texture coordinates that will replace any existing vertex data.
		 */
		void TexCoords(std::vector<glm::vec3>& texcoords) { m_texcoords = texcoords; }

		/**
		 * @brief Sets the vertex texture coordinates for this mesh using a flattened floating point array.
		 * 
		 * @param tc is a pointer to the vertex texture coordinate data
		 * @param N is the number of vertices
		 * @param coords is the number of texture coordinates for each vertex 
		 * 				(this will be stored internally as 3 with any additional vertices set to zero)
		 */
		void TexCoords(float* tc, size_t N, int coords = 3) {
			m_texcoords.clear();
			m_texcoords.resize(N);

			// handle cases where there are fewer than 3 texture coordinates
			if(coords == 3)
				memcpy(&m_texcoords[0], tc, N * sizeof(float) * coords);
			else if (coords == 2) {
				for (size_t ti = 0; ti < N; ti++) {
					m_texcoords[ti].x = tc[ti * 2 + 0];
					m_texcoords[ti].y = tc[ti * 2 + 1];
				}
			}
			else if (coords == 1) {
				for (size_t ti = 0; ti < N; ti++) {
					m_texcoords[ti].x = tc[ti];
				}
			}
			else
				throw std::runtime_error("ERROR: invalid number of texture coordinates");
		}

		/**
		 * @brief Returns the indices describing the mesh connectivity, where every three vertices represent a triangle.
		 * 
		 * @return std::vector<unsigned int> is the index array
		 */
		std::vector<unsigned int> Indices() const { return m_indices; }

		/**
		 * @brief Sets the index array for the current mesh
		 * 
		 * @param indices is the list of vertex indices defining the connectivity for the mesh, which will replace
		 * 					any current indices
		 */
		void Indices(std::vector<unsigned int>& indices) { m_indices = indices; }

	
		/**
		 * @brief Returns the number of triangles in the mesh.
		 * 
		 * @return unsigned int is the total number of triangles in the mesh
		 */
		unsigned int NumTriangles() const { return m_indices.size() / 3; }

		/**
		 * @brief Returns the number of vertices in the mesh
		 * 
		 * @return unsigned int is the total number of vertices in the triangle mesh
		 */
		unsigned int NumVertices() const { return m_vertices.size(); }

		size_t NumIndices() const { return m_indices.size(); }

		/**
		 * @brief Returns the total number of bytes needed to store the entire mesh
		 * 
		 * @return size_t is the total size of the mesh (in bytes)
		 */
		size_t Bytes() {
			size_t nbytes = 0;
			nbytes += m_vertices.size() * sizeof(glm::mat3);
			nbytes += m_normals.size() * sizeof(glm::mat3);
			nbytes += m_texcoords.size() * sizeof(glm::mat3);
			nbytes += m_indices.size() * sizeof(unsigned int);
			return nbytes;
		}

		/**
		 * @brief Returns a vector of interleaved vertex positions, normals, and texture coordinates. This data is commonly
		 * 			used for graphics APIs such as OpenGL for specifying geometry and vertex attributes.
		 * 
		 * @param texcoords is the number of values used to specify the texture coordinates at each vertex. This value is
		 * 					stored internally as 3 texture coordinates per vertex, but many applications rely on fewer. Coordinates
		 * 					will be returned in lexicographic order (x first, then xy, then xyz).
		 * @return std::vector<float> is the array of interleaved vertex data
		 */
		std::vector<float> Interleave(int texcoords = 3) {

			size_t n_vertices = NumVertices();
			size_t stride = 3 + 3 + texcoords;

			// create a vector to hold the interleaved vertices
			std::vector<float> interleaved(n_vertices * stride);


			//size_t stride = _vertex_dim + _normal_dim + _texture_dim;
			for (size_t i = 0; i < NumVertices(); i++) {

				interleaved[i * stride + 0] = m_vertices[i][0];
				interleaved[i * stride + 1] = m_vertices[i][1];
				interleaved[i * stride + 2] = m_vertices[i][2];

				interleaved[i * stride + 3] = m_normals[i][0];
				interleaved[i * stride + 4] = m_normals[i][1];
				interleaved[i * stride + 5] = m_normals[i][2];

				for(size_t tci = 0; tci < texcoords; tci++)
					interleaved[i * stride + 6 + tci] = m_texcoords[i][tci];
			}

			return interleaved;
		}

		/**
		 * @brief Merge two triangle meshes into a single mesh.
		 * 
		 * @param rhs is the second triangle mesh to be merged with this one.
		 * @return tmesh is the merged triangle mesh.
		 */
		tmesh Merge(tmesh rhs) {

			tmesh merged;


			// merge vertices
			merged.m_vertices.reserve(m_vertices.size() + rhs.m_vertices.size()); // preallocate memory
			merged.m_vertices.insert(merged.m_vertices.end(), m_vertices.begin(), m_vertices.end());
			merged.m_vertices.insert(merged.m_vertices.end(), rhs.m_vertices.begin(), rhs.m_vertices.end());

			// merge normals
			merged.m_normals.reserve(m_normals.size() + rhs.m_normals.size()); // preallocate memory
			merged.m_normals.insert(merged.m_normals.end(), m_normals.begin(), m_normals.end());
			merged.m_normals.insert(merged.m_normals.end(), rhs.m_normals.begin(), rhs.m_normals.end());

			// merge texcoords
			merged.m_texcoords.reserve(m_texcoords.size() + rhs.m_texcoords.size()); // preallocate memory
			merged.m_texcoords.insert(merged.m_texcoords.end(), m_texcoords.begin(), m_texcoords.end());
			merged.m_texcoords.insert(merged.m_texcoords.end(), rhs.m_texcoords.begin(), rhs.m_texcoords.end());

			// merge indices
			merged.m_indices.reserve(m_indices.size() + rhs.m_indices.size()); // preallocate memory
			merged.m_indices.insert(merged.m_indices.end(), m_indices.begin(), m_indices.end());
			merged.m_indices.insert(merged.m_indices.end(), rhs.m_indices.begin(), rhs.m_indices.end());

			// increment the indices from the second mesh
			size_t V = NumVertices();
			for (size_t i = m_indices.size(); i < merged.m_indices.size(); i++)
				merged.m_indices[i] += V;

			return merged;
		}

		/**
		 * @brief Assignment operator
		 * 
		 * @param rhs is the triangle mesh to be assigned
		 * @return tmesh& is a reference to a new triangle mesh that is a copy of the right hand side value.
		 */
		tmesh& operator=(const tmesh rhs) {
			m_vertices = rhs.m_vertices;
			m_normals = rhs.m_normals;
			m_texcoords = rhs.m_texcoords;
			m_indices = rhs.m_indices;
			return *this;
		}

		/**
		 * @brief Construct a new triangle mesh by making a copy of an existing mesh
		 * 
		 * @param c the existing mesh to copy.
		 */
		tmesh(const tmesh& c) {
			m_vertices = c.m_vertices;
			m_normals = c.m_normals;
			m_texcoords = c.m_texcoords;
			m_indices = c.m_indices;
		}

		tmesh() {}

		/**
		 * @brief Transforms all vertices in the model by a specified matrix and returns the transformed model.
		 * 
		 * @param matrix is the transformation matrix that will be applied to all vertices.
		 * @param component is the component that is transformed (0 transforms positions and normals, 1 transforms texture coordinates)
		 * @return tmesh is the transformed mesh
		 */
		tmesh Transform(glm::mat3 matrix, int component = 0) {
			tmesh result(*this);

			size_t V = NumVertices();
			for (size_t vi = 0; vi < V; vi++) {
				if (component == 0) {
					result.m_vertices[vi] = matrix * result.m_vertices[vi];
					result.m_normals[vi] = matrix * result.m_normals[vi];
				}
				else if (component == 1) {
					result.m_texcoords[vi] = matrix * result.m_texcoords[vi];
				}
				else
					throw std::runtime_error("ERROR: unknown component type for transformation");
			}
			return result;
			
		}

		/**
		 * @brief Translates all vertices by the specified coordinates.
		 * 
		 * @param coords is a vector specifying the direction and magnitude of the transformation
		 * @param component is the component that is transformed (0 transforms positions, 1 transforms texture coordinates)
		 * @return tmesh is the transformed mesh
		 */
		tmesh Translate(glm::vec3 coords, int component = 0) {
			tmesh result(*this);
			size_t V = NumVertices();
			for (size_t vi = 0; vi < V; vi++) {
				if(component == 0)
					result.m_vertices[vi] += coords;
				else if(component == 1)
					result.m_texcoords[vi] += coords;
			}

			return result;
		}


		/**
		 * @brief Scales all vertices by the specified coordinates.
		 * 
		 * @param coords is a vector specifying the magnitude of the scaling along each dimension
		 * @param component is the component that is transformed (0 transforms positions, 1 transforms texture coordinates)
		 * @return tmesh is the transformed mesh
		 */
		tmesh Scale(glm::vec3 coords, int component = 0) {
			
			tmesh result(*this);
			size_t V = NumVertices();
			
			for (size_t vi = 0; vi < V; vi++) {
				if(component == 0)
					result.m_vertices[vi] *= coords;
				else if(component == 1)
					result.m_texcoords[vi] *= coords;
			}
			return result;
		}

		/**
		 * @brief Creates a mesh composed of tiled replicas of the current mesh.
		 * 
		 * @param d is the direction along which the current mesh is replicated (how each new mesh is translated)
		 * @param N is the number of copies to make
		 * @return tmesh is the new tiled mesh
		 */
		tmesh Tile(glm::vec3 d, size_t N) const {
			const unsigned int nV = NumVertices();
			tmesh result;
			result.m_vertices.resize(N * m_vertices.size());			// reserve space in the output mesh for all sub-elements
			result.m_normals.resize(N * m_normals.size());
			result.m_texcoords.resize(N * m_texcoords.size());
			result.m_indices.resize(N * m_indices.size());

			tmesh temp = (*this);
			for(unsigned int ni = 0; ni < N; ni++) {
				size_t v_offset = ni * m_vertices.size();
				size_t n_offset = ni * m_normals.size();
				size_t t_offset = ni * m_texcoords.size();
				size_t i_offset = ni * m_indices.size();
				std::copy(temp.m_vertices.begin(), temp.m_vertices.end(), result.m_vertices.begin() + v_offset);
				std::copy(temp.m_normals.begin(), temp.m_normals.end(), result.m_normals.begin() + n_offset);
				std::copy(temp.m_texcoords.begin(), temp.m_texcoords.end(), result.m_texcoords.begin() + t_offset);
				std::transform(temp.m_indices.begin(), temp.m_indices.end(), result.m_indices.begin() + i_offset, std::bind(std::plus<unsigned int>(), std::placeholders::_1, ni * nV));

				temp = temp.Translate(d);
			}

			return result;
		}
	};
}	// end namespace tira
