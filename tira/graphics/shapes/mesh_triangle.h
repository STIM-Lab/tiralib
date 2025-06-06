 #pragma once
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <functional>
#include <exception>
#include <iostream>

/* The TIRA trimesh class is designed to store data for efficiently render triangular meshes
	using graphics APIs. It provides details for storing vertices, normals, and texture coordinates.
 */

namespace tira{
	template <typename T>
	class trimesh {
	protected:
		std::vector<T> _vertices;
		std::vector<T> _normals;
		std::vector<T> _texcoords;
		std::vector<unsigned int> _indices;
		unsigned int _vertex_dim = 3;
		unsigned int _normal_dim = 3;
		unsigned int _texture_dim = 2;

		/// <summary>
		/// Tests the compatibility of a merge between this mesh and the requested trimesh.
		/// The number of vertex, texture, and normal coordinates must be the same for it to work.
		/// </summary>
		/// <param name="merging"></param>
		/// <returns></returns>
		bool _check_merge(trimesh& merging) {
			if (_vertex_dim != merging._vertex_dim) {
				std::cout << "ERROR: incompatible trimesh merge, A has " << _vertex_dim << " vertex coordinates and B has " << merging._vertex_dim << std::endl;
				return false;
			}
			if (_normal_dim != merging._normal_dim) {
				std::cout << "ERROR: incompatible trimesh merge, A has " << _normal_dim << " normal coordinates and B has " << merging._normal_dim << std::endl;
				return false;
			}
			if (_texture_dim != merging._texture_dim) {
				std::cout << "ERROR: incompatible trimesh merge, A has " << _texture_dim << " texture coordinates and B has " << merging._texture_dim << std::endl;
				return false;
			}
			return true;
		}

	public:
		// deprecated
		std::vector<T> getVertices() const { return _vertices; }
		std::vector<T> getNormals() const { return _normals; }
		std::vector<T> getTexCoords() const { return _texcoords; }
		std::vector<unsigned int> getIndices() const { return _indices; }
		unsigned int getVertexDim() const { return _vertex_dim; }
		unsigned int getNormalDim() const { return _normal_dim; }
		unsigned int getTextureDim() const { return _texture_dim; }

		// get and set member functions
		std::vector<T> vertices() const { return _vertices; }
		void vertices(std::vector<T>& vertices) { _vertices = vertices; }

		std::vector<T> normals() const { return _normals; }
		void normals(std::vector<T>& normals) { _normals = normals; }

		std::vector<T> texcoords() const { return _texcoords; }
		void texcoords(std::vector<T>& texcoords) { _texcoords = texcoords; }

		std::vector<unsigned int> indices() const { return _indices; }
		void indices(std::vector<unsigned int>& indices) { _indices = indices; }

		unsigned int vdim() const { return _vertex_dim; }
		unsigned int ndim() const { return _normal_dim; }
		unsigned int tdim() const { return _texture_dim; }
		unsigned int num() const { return _indices.size() / 3; }
		unsigned int num_v() const { return _vertices.size() / _vertex_dim; }
		unsigned int num_n() const { return _normals.size() /  _normal_dim; }
		unsigned int num_t() const { return _texcoords.size() /  _texture_dim; }

		size_t getNumVertices() const { return _vertices.size() / _vertex_dim; }
		size_t getNumTriangles() const { return _indices.size() / 3; }

		size_t bytes() {
			size_t nbytes = 0;
			nbytes += _vertices.size() * sizeof(T);
			nbytes += _normals.size() * sizeof(T);
			nbytes += _texcoords.size() * sizeof(T);
			nbytes += _indices.size() * sizeof(unsigned int);
			return nbytes;
		}

		// Returns a vector of vertices with the vertex positions, normals, and texture coordinates interleaved
		std::vector<T> getInterleavedVertices() {
			// create a vector to hold the interleaved vertices
			std::vector<T> interleaved(_vertices.size() +
				_normals.size() +
				_texcoords.size());


			size_t stride = _vertex_dim + _normal_dim + _texture_dim;
			for (size_t i = 0; i < getNumVertices(); i++) {
				for (size_t v = 0; v < _vertex_dim; v++)
					interleaved[i * stride + v] = _vertices[i * _vertex_dim + v];
				for (size_t n = 0; n < _normal_dim; n++)
					interleaved[i * stride + _vertex_dim + n] = _normals[i * _normal_dim + n];
				for (size_t t = 0; t < _texture_dim; t++)
					interleaved[i * stride + _vertex_dim + _normal_dim + t] = _texcoords[i * _texture_dim + t];
			}

			return interleaved;
		}

		/// <summary>
		/// Merge two geometry objects and return the result
		/// </summary>
		/// <param name="rhs"></param>
		/// <returns></returns>
		trimesh<T> merge(trimesh<T> rhs) {
			if (!_check_merge(rhs))								// if the meshes have different numbers of coordinates, throw an exception
				throw std::invalid_argument("Incompatible Mesh Parameters");

			trimesh<T> merged;
			merged._vertex_dim = _vertex_dim;
			merged._normal_dim = _normal_dim;
			merged._texture_dim = _texture_dim;

			// merge vertices
			merged._vertices.reserve(_vertices.size() + rhs._vertices.size()); // preallocate memory
			merged._vertices.insert(merged._vertices.end(), _vertices.begin(), _vertices.end());
			merged._vertices.insert(merged._vertices.end(), rhs._vertices.begin(), rhs._vertices.end());

			// merge normals
			merged._normals.reserve(_normals.size() + rhs._normals.size()); // preallocate memory
			merged._normals.insert(merged._normals.end(), _normals.begin(), _normals.end());
			merged._normals.insert(merged._normals.end(), rhs._normals.begin(), rhs._normals.end());

			// merge texcoords
			merged._texcoords.reserve(_texcoords.size() + rhs._texcoords.size()); // preallocate memory
			merged._texcoords.insert(merged._texcoords.end(), _texcoords.begin(), _texcoords.end());
			merged._texcoords.insert(merged._texcoords.end(), rhs._texcoords.begin(), rhs._texcoords.end());

			// merge indices
			merged._indices.reserve(_indices.size() + rhs._indices.size()); // preallocate memory
			merged._indices.insert(merged._indices.end(), _indices.begin(), _indices.end());
			merged._indices.insert(merged._indices.end(), rhs._indices.begin(), rhs._indices.end());

			// increment the indices from the second mesh
			size_t V = getNumVertices();
			for (size_t i = _indices.size(); i < merged._indices.size(); i++)
				merged._indices[i] += V;

			return merged;
		}

		trimesh<T>& operator=(const trimesh<T> rhs) {
			_vertex_dim = rhs._vertex_dim;
			_normal_dim = rhs._normal_dim;
			_texture_dim = rhs._texture_dim;
			_vertices = rhs._vertices;
			_normals = rhs._normals;
			_texcoords = rhs._texcoords;
			_indices = rhs._indices;
			return *this;
		}

		trimesh(const trimesh& c) {
			_vertex_dim = c._vertex_dim;
			_normal_dim = c._normal_dim;
			_texture_dim = c._texture_dim;
			_vertices = c._vertices;
			_normals = c._normals;
			_texcoords = c._texcoords;
			_indices = c._indices;
		}

		trimesh() {}

		trimesh<T> translate(std::vector<T> coords, int component = 0) {
			if ( (component == 0 && coords.size() > _vertex_dim) || (component == 1 && coords.size() > _texture_dim))
				throw std::invalid_argument("Too many translation coordinates provided");

			trimesh<T> result(*this);
			size_t V = getNumVertices();
			for (size_t d = 0; d < coords.size(); d++) {
				if (coords[d] != 0) {
					for (size_t vi = 0; vi < V; vi++) {
						if(component == 0)
							result._vertices[vi * _vertex_dim + d] += coords[d];
						else if(component == 1)
							result._texcoords[vi * _texture_dim + d] += coords[d];
					}
				}
			}
			return result;
		}


		trimesh<T> scale(std::vector<T> coords, int component = 0) {
			if ( (component == 0 && coords.size() > _vertex_dim) || (component == 1 && coords.size() > _texture_dim))
				throw std::invalid_argument("Too many scale coordinates provided");

			trimesh<T> result(*this);
			size_t V = getNumVertices();
			for (size_t d = 0; d < coords.size(); d++) {
				if (coords[d] != 1) {
					for (size_t vi = 0; vi < V; vi++) {
						if(component == 0)
							result._vertices[vi * _vertex_dim + d] *= coords[d];
						else if(component == 1)
							result._texcoords[vi * _texture_dim + d] *= coords[d];
					}
				}
			}
			return result;
		}

		trimesh<T> tile(std::vector<T> d, size_t N) const {
			const unsigned int nV = getNumVertices();
			trimesh<T> result;
			result._vertices.resize(N * _vertices.size());			// reserve space in the output mesh for all sub-elements
			result._normals.resize(N * _normals.size());
			result._texcoords.resize(N * _texcoords.size());
			result._indices.resize(N * _indices.size());

			trimesh<T> temp = (*this);
			for(unsigned int ni = 0; ni < N; ni++) {
				size_t v_offset = ni * _vertices.size();
				size_t n_offset = ni * _normals.size();
				size_t t_offset = ni * _texcoords.size();
				size_t i_offset = ni * _indices.size();
				std::copy(temp._vertices.begin(), temp._vertices.end(), result._vertices.begin() + v_offset);
				std::copy(temp._normals.begin(), temp._normals.end(), result._normals.begin() + n_offset);
				std::copy(temp._texcoords.begin(), temp._texcoords.end(), result._texcoords.begin() + t_offset);
				std::transform(temp._indices.begin(), temp._indices.end(), result._indices.begin() + i_offset, std::bind(std::plus<unsigned int>(), std::placeholders::_1, ni * nV));

				temp = temp.translate(d);
			}

			result._normal_dim = _normal_dim;
			result._texture_dim = _texture_dim;
			result._vertex_dim = _vertex_dim;

			return result;
		}

		/*// save an OBJ file for this geometry
		void obj(std::string filename) {
			std::ofstream f(filename);										// open a file for writing
			if (!f) throw "ERROR: unable to open OBJ file for writing";
			size_t nV = getNumVertices();									// get the number of vertices
			for (size_t vi = 0; vi < nV; vi++) {							// for each vertex
				f << "v ";													// output the initial "v" to indicate the vertex
				for (size_t ci = 0; ci < _vertex_dim; ci++) {				// for each coordinate
					f << _vertices[3 * vi + ci];							// output the coordinate value
					if (ci < _vertex_dim - 1) f << " ";						// output a space (if this isn't the last coordinate)
				}
				if (vi < nV - 1)f << std::endl;								// move to the next line
			}

			if (_texcoords.size()) {										// if there are any texture coordinates
				f << std::endl;
				for (size_t ti = 0; ti < nV; ti++) {						// for each vertex
					f << "vt ";												// output the initial "v" to indicate the vertex
					for (size_t ci = 0; ci < _texture_dim; ci++) {			// for each coordinate
						f << _texcoords[_texture_dim * ti + ci];						// output the coordinate value
						if (ci < _texture_dim - 1) f << " ";					// output a space (if this isn't the last coordinate)
					}
					if (ti < nV - 1)f << std::endl;							// move to the next line
				}
			}

			if (_normals.size()) {										// if there are any texture coordinates
				f << std::endl;
				for (size_t ni = 0; ni < nV; ni++) {						// for each vertex
					f << "vn ";												// output the initial "v" to indicate the vertex
					for (size_t ci = 0; ci < _normal_dim; ci++) {			// for each coordinate
						f << _normals[_normal_dim * ni + ci];						// output the coordinate value
						if (ci < _normal_dim - 1) f << " ";					// output a space (if this isn't the last coordinate)
					}
					if (ni < nV - 1)f << std::endl;							// move to the next line
				}
			}

			f << std::endl;
			size_t nT = getNumTriangles();								//get the number of triangles
			for (size_t ti = 0; ti <nT; ti++) {					// for each vertex
				f << "f ";												// output the initial "v" to indicate the vertex
				for (size_t i = 0; i < 3; i++) {						// for each point of the triangle
					f << _indices[3 * ti + i] + 1;						// output the coordinate value
					if (i < 2) f << " ";					// output a space (if this isn't the last coordinate)
				}
				if (ti < nT - 1)f << std::endl;							// move to the next line
			}
			f.close();														// close the OBJ file
		}*/
	};
}	// end namespace tira
