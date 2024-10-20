#pragma once
#include <vector>
#include <string>
#include <fstream>

namespace tira{
	template <typename T>
	class geometry {
	protected:
		std::vector<T> _vertices;
		std::vector<T> _normals;
		std::vector<T> _texcoords;
		std::vector<unsigned int> _indices;
		unsigned int _vertex_dim = 3;
		unsigned int _normal_dim = 3;
		unsigned int _texture_dim = 2;

		/// <summary>
		/// Tests the compatibility of a merge between geometry mesh and the merging geometry.
		/// The number of vertex, texture, and normal coordinates must be the same for it to work.
		/// </summary>
		/// <param name="merging"></param>
		/// <returns></returns>
		bool _check_merge(geometry& merging) {
			if (_vertex_dim != merging._vertex_dim) {
				std::cout << "ERROR: incompatible geometry merge, A has " << _vertex_dim << " vertex coordinates and B has " << merging._vertex_dim << std::endl;
				return false;
			}
			if (_normal_dim != merging._normal_dim) {
				std::cout << "ERROR: incompatible geometry merge, A has " << _normal_dim << " normal coordinates and B has " << merging._normal_dim << std::endl;
				return false;
			}
			if (_texture_dim != merging._texture_dim) {
				std::cout << "ERROR: incompatible geometry merge, A has " << _texture_dim << " texture coordinates and B has " << merging._texture_dim << std::endl;
				return false;
			}
			return true;
		}

	public:
		std::vector<T> getVertices() { return _vertices; }
		std::vector<T> getNormals() { return _normals; }
		std::vector<T> getTexCoords() { return _texcoords; }
		unsigned int getVertexDim() { return _vertex_dim; }
		unsigned int getNormalDim() { return _normal_dim; }
		unsigned int getTextureDim() { return _texture_dim; }
		std::vector<unsigned int> getIndices() { return _indices; }
		size_t getNumVertices() { return _vertices.size() / _vertex_dim; }
		size_t getNumTriangles() { return _indices.size() / 3; }

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
		geometry<T> merge(geometry<T> rhs) {
			if (!_check_merge(rhs))								// if the meshes have different numbers of coordinates, throw an exception
				throw std::exception("Incompatible Mesh Parameters");

			geometry<T> merged;
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

		/// <summary>
		/// Combine two geometry objects and return the result
		/// </summary>
		/// <param name="rhs"></param>
		/// <returns></returns>
		geometry<T> operator+(geometry<T> rhs) {

		}

		geometry(const geometry& c) {
			_vertex_dim = c._vertex_dim;
			_normal_dim = c._normal_dim;
			_texture_dim = c._texture_dim;
			_vertices = c._vertices;
			_normals = c._normals;
			_texcoords = c._texcoords;
			_indices = c._indices;
		}

		geometry() {}

		geometry<T> translate(std::vector<T> coords) {
			if (coords.size() > _vertex_dim)
				throw std::exception("Too many translation coordinates provided");

			geometry<T> result(*this);
			size_t V = getNumVertices();
			for (size_t d = 0; d < coords.size(); d++) {
				if (coords[d] != 0) {
					for (size_t vi = 0; vi < V; vi++) {
						result._vertices[vi * _vertex_dim + d] += coords[d];
					}
				}
			}
			return result;
		}

		geometry<T> scale(std::vector<T> coords) {
			if (coords.size() > _vertex_dim)
				throw std::exception("Too many translation coordinates provided");

			geometry<T> result(*this);
			size_t V = getNumVertices();
			for (size_t d = 0; d < coords.size(); d++) {
				if (coords[d] != 1) {
					for (size_t vi = 0; vi < V; vi++) {
						result._vertices[vi * _vertex_dim + d] *= coords[d];
					}
				}
			}
			return result;
		}

		// save an OBJ file for this geometry
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
		}
	};
}