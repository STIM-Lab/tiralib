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
    class linemesh {
    protected:
        std::vector<T> _vertices;
        std::vector<T> _normals;
        std::vector<T> _texcoords;
        std::vector< std::vector<unsigned int> > _strips;
        unsigned int _vertex_dim = 3;
        unsigned int _normal_dim = 3;
        unsigned int _texture_dim = 1;

        std::vector<T> _line_entry(T* array, size_t line_idx, unsigned ndims) {
            size_t n_vertices = _strips[line_idx].size();                    // get the number of vertices in the line

            std::vector<T> result(n_vertices * _vertex_dim);                 // allocate space for the output vertex vector

            for (size_t vi = 0; vi < n_vertices; vi++) {
                size_t vertex_idx = _strips[line_idx][vi];                    // get the ID for the current vertex
                for (size_t di = 0; di < _vertex_dim; di++) {
                    result[vi * _vertex_dim + di] = array[vertex_idx * _vertex_dim + di];
                }
            }
            return result;
        }

    public:

        // get and set member functions
        std::vector<T> vertices() const { return _vertices; }
        void vertices(std::vector<T>& vertices) { _vertices = vertices; }
        /// Returns the vertices associated with a given line index
        std::vector<T> vertices(unsigned line_idx) {
            return _line_entry(_vertices.data(), line_idx, _vertex_dim);
        }

        std::vector<T> normals() const { return _normals; }
        void normals(std::vector<T>& normals) { _normals = normals; }
        std::vector<T> normals(unsigned line_idx) {
            return _line_entry(_normals.data(), line_idx, _normal_dim);
        }

        std::vector<T> texcoords() const { return _texcoords; }
        void texcoords(std::vector<T>& texcoords) { _texcoords = texcoords; }
        std::vector<T> texcoords(unsigned line_idx) {
            return _line_entry(_texcoords.data(), line_idx, _texture_dim);
        }

        std::vector<unsigned int> indices(unsigned line_idx) const { return _strips[line_idx]; }
        void indices(std::vector<unsigned int>& indices, unsigned line_idx) { _strips.push_back(indices); }

        unsigned int vdim() const { return _vertex_dim; }
        unsigned int ndim() const { return _normal_dim; }
        unsigned int tdim() const { return _texture_dim; }
        unsigned int num_lines() const { return _strips.size(); }
        unsigned int num_v() const { return _vertices.size() / _vertex_dim; }
        unsigned int num_n() const { return _normals.size() /  _normal_dim; }
        unsigned int num_t() const { return _texcoords.size() /  _texture_dim; }

        size_t getNumVertices() const { return _vertices.size() / _vertex_dim; }
    };    // end linemesh class
}         // end namespace tira
