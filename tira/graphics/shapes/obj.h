#include <vector>
#include <string>
#include <sstream>

#include "mesh_triangle.h"

namespace tira{
    template <typename T>
    class obj {
        protected:
            std::vector< trimesh<T> > _tmeshes;

        public:
            void write_entries(std::stringstream& s, std::vector<T> a, std::string prefix, unsigned int dims) {
                for (size_t i = 0; i < a.size(); i+=dims) {
                    s << prefix << " ";
                    for(size_t di = 0; di < dims; di++)
                        s << " " << a[i + di];
                    s << std::endl;
                }
            }

            void write_vertices(std::stringstream& ss) {
                for (size_t tmi = 0; tmi < _tmeshes.size(); tmi++) {
                    std::vector<T> v = _tmeshes[tmi].vertices();
                    write_entries(ss, v, "v", _tmeshes[tmi].vdim());
                }
            }

            void write_normals(std::stringstream& ss) {
                for (size_t tmi = 0; tmi < _tmeshes.size(); tmi++) {
                    std::vector<T> n = _tmeshes[tmi].normals();
                    write_entries(ss, n, "vn", _tmeshes[tmi].ndim());
                }
            }

            void write_texcoords(std::stringstream& ss) {
                for (size_t tmi = 0; tmi < _tmeshes.size(); tmi++) {
                    std::vector<T> t = _tmeshes[tmi].texcoords();
                    write_entries(ss, t, "vt", _tmeshes[tmi].tdim());
                }
            }

            void write_indices(std::stringstream& ss) {
                size_t offset = 1;

                for (size_t tmi = 0; tmi < _tmeshes.size(); tmi++) {
                    std::vector<unsigned int> i = _tmeshes[tmi].indices();
                    unsigned int N = i.size() / 3;                            // calculate the number of triangles
                    for(unsigned int ti = 0; ti < N; ti++) {
                        ss << "f ";
                        for(unsigned int vi = 0; vi < 3; vi++) {
                            ss << i[ti * 3 + vi] + offset << "/";                        // output the vertex index
                            ss << i[ti * 3 + vi] + offset<< "/";
                            ss << i[ti * 3 + vi] + offset;
                            if(vi < 2) ss << " ";
                        }
                        ss << std::endl;
                    }
                    offset += _tmeshes[tmi].num_v();
                }
            }

            void add(trimesh<float> t) {
                _tmeshes.push_back(t);
            }

            void save(const std::string& filename) {
                std::stringstream ss;

                // write all of the vertex entries
                write_vertices(ss);
                write_normals(ss);
                write_texcoords(ss);
                write_indices(ss);

                std::ofstream ofs(filename);
                ofs.write(ss.str().c_str(), ss.str().size());
                ofs.close();
            }



    };
}    // end namespace tira
