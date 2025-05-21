#include <vector>
#include <string>
#include <sstream>

#include "mesh_triangle.h"
#include "mesh_line.h"

namespace tira{
    // Enumeration for primitive types that will be rendered to the OBJ object
    enum OBJenum {
        OBJ_POINTS,
        OBJ_LINES,
        OBJ_LINE_STRIP,
        OBJ_LINE_LOOP,
        OBJ_TRIANGLES,
        OBJ_TRIANGLE_STRIP,
        OBJ_TRIANGLE_FAN,
        OBJ_QUADS,
        OBJ_QUAD_STRIP,
        OBJ_POLYGON
    };

    template <typename T>
    class obj {
    public:
        

    protected:

        typedef std::vector<T> _vec;

        std::vector<_vec> _vertices;
        std::vector<_vec> _normals;
        std::vector<_vec> _texcoords;

        struct _index {
            size_t v { 0 };
            size_t t { 0 };
            size_t n { 0 };
        };

        typedef std::vector<_index> _entry;

        std::vector<_entry> _faces;
        std::vector<_entry> _lines;
        std::vector<_entry> _points;

        // immediate mode members
        bool _drawing = false;                  // flags whether or not a primitive is actively being drawn
        bool _active_normal = false;            // flags whether or not the current vertex will have a normal
        bool _active_texcoord = false;          // flags whether or not the current vertex will have a texture coordinate

        _entry _active_primitive;



        OBJenum _drawtype = OBJ_TRIANGLES;
        void _Vertex() {
            // add the appropriate indices to the active primitive
            _index new_idx;
            new_idx.v = _vertices.size();
            if (_active_normal == true) new_idx.n = _normals.size();
            if (_active_texcoord == true) new_idx.t = _texcoords.size();
            _active_primitive.push_back(new_idx);
        }


        /// Write a vector entry to the OBJ file stream. This is includes the coordinates for a vertex,
        /// normal, or texture coordinate.
        static void _write_vector_entry(std::stringstream& ss, std::vector<_vec>& _array, std::string prefix) {
            for (size_t i = 0; i < _array.size(); i++) {                        // for each vector
                ss << prefix;                                                   // output the prefix (normally "v", "vt", or "vn")
                for (size_t di = 0; di < _array[i].size(); di++) {              // for each coordinate
                    ss << " " << _array[i][di];                                 // output the coordinate value
                }
                ss << std::endl;
            }
        }

        static void _write_index(std::stringstream& ss, _index e) {
            ss << e.v;                                                          // the vertex component is required
            if (e.t != 0 || e.n != 0) ss << "/";                                // if either a texture or normal is provided, output a slash
            if (e.t != 0) ss << e.t;                                            // if a texture coordinate is provided, output the index
            if (e.n != 0) ss << "/" << e.n;                                     // if a normal is provided, output a slash and the index
        }

        void _write_faces(std::stringstream& ss) {
            for (size_t fi = 0; fi < _faces.size(); fi++) {                     // for each face
                ss << "f";                                                      // write an "f"
                for (size_t vi = 0; vi < _faces[fi].size(); vi++) {             // for each vertex
                    ss << " ";                                                  // write a space
                    _write_index(ss, _faces[fi][vi]);                        // write the index
                }
                ss << std::endl;
            }
        }
        void _write_lines(std::stringstream& ss) {
            for (size_t li = 0; li < _lines.size(); li++) {                     // for each face
                ss << "l";                                                      // write an "f"
                for (size_t vi = 0; vi < _lines[li].size(); vi++) {             // for each vertex
                    ss << " ";                                                  // write a space
                    _write_index(ss, _lines[li][vi]);                        // write the index
                }
                ss << std::endl;
            }
        }

    public:
        template<typename D>
        void add(trimesh<D> mesh) {
            // Calculate the offsets required to insert this geometry into the current structure
            const size_t v_offset = _vertices.size() + 1;
            const size_t n_offset = _normals.size() + 1;
            const size_t t_offset = _texcoords.size() + 1;

            // Get the vertex, normal, and texture coordinates from the mesh structure
            std::vector<D> Vertices = mesh.vertices();
            std::vector<D> Normals = mesh.normals();
            std::vector<D> TexCoords = mesh.texcoords();

            // The trimesh structure stores a vertex, normal, and texture coordinate for each
            // vertex, so we only have to loop through each of these arrays once.
            for (size_t vi = 0; vi < mesh.num_v(); vi++) {                  // for each vertex

                // Allocate space for each coordinate of the vertex and fill that structure with the
                // provided coordinate values from the mesh.
                _vec v(mesh.vdim());
                for (size_t di = 0; di < mesh.vdim(); di++) {
                    v[di] = Vertices[vi * mesh.vdim() + di];
                }

                // Do this same pre-allocation and assignment for the normals and texture coordinates
                _vec n(mesh.ndim());
                for (size_t di = 0; di < mesh.ndim(); di++) {
                    n[di] = Normals[vi * mesh.ndim() + di];
                }
                _vec t(mesh.tdim());
                for (size_t di = 0; di < mesh.tdim(); di++) {
                    t[di] = TexCoords[vi * mesh.tdim() + di];
                }

                // Push the coordinates for each vertex component to the end of the associated vector
                // in the OBJ structure.
                _vertices.push_back(v);
                _normals.push_back(n);
                _texcoords.push_back(t);
            }

            const std::vector<unsigned int> Faces = mesh.indices();

            const size_t num_triangles = Faces.size() / 3;                  // calculate the number of triangles
            for (size_t fi = 0; fi < num_triangles; fi++) {             // for each triangle in the trimesh

                // Create an OBJ entry - in this case a triangle that has three vertices
                _entry obj_triangle(3);


                for (size_t vi = 0; vi < 3; vi++) {                     // for each vertex in the triangle
                    unsigned i = Faces[fi * 3 + vi];                        // store the vertex index from the trimesh class

                    // The trimesh structure doesn't keep separate indices for each vertex component,
                    // so making that compatible with an OBJ file requires duplicating the indices for
                    // each texture coordinate and normal. They will be adjusted by current offsets from
                    // geometry currently placed in this OBJ structure.
                    _index new_idx = {
                        i + v_offset,
                        i + t_offset,
                        i + n_offset
                    };

                    obj_triangle[vi] = new_idx;                         // add the new vertex to the OBJ entry for the triangle
                }
                _faces.push_back(obj_triangle);                         // push the triangle to the faces array
            }
        }

        void save(const std::string& filename) {
            std::stringstream ss;

            // Write all of the vertex entries
            _write_vector_entry(ss, _vertices, "v");
            _write_vector_entry(ss, _normals, "vn");
            _write_vector_entry(ss, _texcoords, "vt");

            _write_faces(ss);                                             // write the faces
            _write_lines(ss);


            std::ofstream ofs(filename);
            ofs.write(ss.str().c_str(), ss.str().size());
            ofs.close();
        }

        /// OpenGL Overloads - This section contains functions that are designed to mimic
        /// OpenGL's "immediate mode".



        void Begin(OBJenum mode) {
            if (_drawing == true)
                throw std::runtime_error("Currently drawing, call End() first");
            _drawing = true;
            _drawtype = mode;
        }

        void End() {
            switch (_drawtype) {
            case OBJ_TRIANGLES: {
                if (_active_primitive.size() % 3 != 0)
                    throw std::runtime_error("OBJ_TRIANGLES must contain a number of vertices divisible by 3");
                const size_t num_triangles = _active_primitive.size() / 3;
                for (size_t ti = 0; ti < num_triangles; ti++) {
                    _entry new_triangle(3);
                    for (size_t vi = 0; vi < 3; vi++) {
                        new_triangle[vi] = _active_primitive[ti * 3 + vi];
                    }
                    _faces.push_back(new_triangle);
                }
                break;
            }
            case OBJ_LINE_STRIP: {
                if (_active_primitive.size() < 2)
                    throw std::runtime_error("OBJ_LINE_STRIP must have at least 2 vertices");
                _lines.push_back(_active_primitive);
                break;
            }
            default:
                throw std::runtime_error("Primitive type not supported");
                break;
            }

            _drawing = false;                       // turn off the drawing flag
            _active_normal = false;                 // deactivate the normal
            _active_texcoord = false;               // deactivate the texture coordinate
            _active_primitive.clear();              // clear the active primitive
        }

        /// Immediate Mode methods for assigning normals to the current vertex
        void Normal(T x) {
            _vec new_normal = { x };                // create a new normal and push it into the normals array
            _normals.push_back(new_normal);
            _active_normal = true;
        }
        void Normal(T x, T y) {
            _vec new_normal = { x, y };
            _normals.push_back(new_normal);
            _active_normal = true;
        }
        void Normal(T x, T y, T z) {
            _vec new_normal = { x, y, z };
            _normals.push_back(new_normal);
            _active_normal = true;
        }

        /// Immediate mode methods for assigning texture coordinates to the current vertex
        void TexCoord(T u) {
            _vec new_texcoord = { u };
            _texcoords.push_back(new_texcoord);
            _active_texcoord = true;
        }
        void TexCoord(T u, T v) {
            _vec new_texcoord = { u, v };
            _texcoords.push_back(new_texcoord);
            _active_texcoord = true;
        }
        void TexCoord(T u, T v, T w) {
            _vec new_texcoord = { u, v, w };
            _texcoords.push_back(new_texcoord);
            _active_texcoord = true;
        }

        void Vertex(T x) {
            _vec new_vertex = { x };
            _vertices.push_back(new_vertex);
            _Vertex();
        }
        void Vertex(T x, T y) {
            _vec new_vertex = { x, y };
            _vertices.push_back(new_vertex);
            _Vertex();
        }
        void Vertex(T x, T y, T z) {
            _vec new_vertex = { x, y, z };
            _vertices.push_back(new_vertex);
            _Vertex();
        }
        void Vertex(T x, T y, T z, T w) {
            _vec new_vertex = { x, y, z, w };
            _vertices.push_back(new_vertex);
            _Vertex();
        }
    };
}    // end namespace tira
