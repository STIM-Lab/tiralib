#include <vector>
#include <string>
#include <sstream>

#include <tira/shapes/tmesh.h>
#include <tira/shapes/lmesh.h>

namespace tira{
    
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

    struct mtl {
        std::string name;
        float Ka[3] = { 1.0f, 0.0f, 0.0f };         // ambient color
        float Kd[3] = { 1.0f, 0.0f, 0.0f };         // diffuse color
        float Ks[3] = { 0.0f, 0.0f, 0.0f };         // specular color
        float Ns = 0.0f;                            // specular exponent
        float d = 1.0f;                             // dissolve (1.0 is opaque)
    };

    /**
     * @brief This class defines an interface for loading and saving Wavefront OBJ files. The interface
     *          replicates OpenGL immediate mode as a convenient way to specify geomtry.
     */
    class obj {
        

    protected:

        /**
         * @brief Vector data type used to store point positions for vertices, normals, and texture coordinates
         */
        typedef std::vector<float> _vec;

        /**
         * @brief Vector storing a list of vertex positions
         */
        std::vector<_vec> m_vertices;

        /**
         * @brief Vector storing a list of normal directions
         */
        std::vector<_vec> m_normals;

        /**
         * @brief Vector storing a list of texture coordinates
         */
        std::vector<_vec> m_texcoords;


        std::vector<mtl> m_materials;

        /**
         * @brief The _index structure stores indices for a single vertex including the position, texture coordinates, and normals. An
         *          index of 0 indicates that the specified component isn't used.
         */
        struct _index {
            size_t v { 0 };
            size_t t { 0 };
            size_t n { 0 };
        };

        /**
         * @brief A _primitive is an std::vector of indices specifying a single primitive, such as a face, line, or point set
         */
        typedef std::vector<_index> _primitive;

        struct _object {
            std::string name;
            size_t material;                    // specifies the material index used for this object (-1 is no material)
            std::vector<_primitive> faces;
            std::vector<_primitive> lines;
            std::vector<_primitive> points;
        };


        std::vector<_object> m_objects;

        /**
         * @brief Flag specifies whether or not a primitive is actively being drawn.
         */
        bool m_drawing = false;

        /**
         * @brief Flag specifies whether or not the current vertex will have a normal
         */
        bool m_active_normal = false;

        /**
         * @brief Flag specifies whether or not the current vertex will have texture coordinates
         */
        bool m_active_texcoord = false;

        /**
         * @brief Specifies the index of the current active object
         */
        size_t m_active_object;

        /**
         * @brief Stores an active primitive that is currently being drawn
         */
        _primitive m_active_primitive;


        /**
         * @brief Stores the type for the primitive that is currently being drawn
         */
        OBJenum m_drawtype = OBJ_TRIANGLES;

        /**
         * @brief Helper function adds the most recently specified vertex components (position, normal, texture coord) to the current
         *          active primitive
         */
        void m_Vertex() {
            _index new_idx;
            new_idx.v = m_vertices.size();
            if (m_active_normal == true) new_idx.n = m_normals.size();
            if (m_active_texcoord == true) new_idx.t = m_texcoords.size();
            m_active_primitive.push_back(new_idx);
        }


        /**
         * @brief Writes all coordinates in an array to an std::stringstream using the specified prefix. This is used to write
         *          vertex components such as position, normal, and texture coordinates using the Wavefront OBJ file format.
         * 
         * @param ss is the std::stringstream that the coordinates will be written to
         * @param _array is an array containing the coordinates to be written
         * @param prefix is the prefix placed before each coordinate that signifies the vertex component it belongs to (ex. "vt" specifies that the
         *                  coordinate is associated with a texture)
         */
        static void m_WriteCoordinateArray(std::stringstream& ss, std::vector<_vec>& _array, std::string prefix) {
            for (size_t i = 0; i < _array.size(); i++) {                        // for each vector
                ss << prefix;                                                   // output the prefix (normally "v", "vt", or "vn")
                for (size_t di = 0; di < _array[i].size(); di++) {              // for each coordinate
                    ss << " " << _array[i][di];                                 // output the coordinate value
                }
                ss << std::endl;
            }
        }

        /**
         * @brief Write an index item out as a string following the Wavefront OBJ file format (ex. 1/2/3)
         * 
         * @param ss an std::stringstream to write the index to
         * @param e is the index to write
         */
        static void m_WriteIndex(std::stringstream& ss, _index e) {
            ss << e.v;                                                          // the vertex component is required
            if (e.t != 0 || e.n != 0) ss << "/";                                // if either a texture or normal is provided, output a slash
            if (e.t != 0) ss << e.t;                                            // if a texture coordinate is provided, output the index
            if (e.n != 0) ss << "/" << e.n;                                     // if a normal is provided, output a slash and the index
        }

        /**
         * @brief Writes all of the faces associated with an object using the Wavefront OBJ file format
         * 
         * @param ss is an std::stringstream to write the face to
         * @param idx is the index of the object to write
         */
        void m_WriteFaces(std::stringstream& ss, size_t idx = 0) {
            for (size_t fi = 0; fi < m_objects[idx].faces.size(); fi++) {                     // for each face
                ss << "f";                                                      // write an "f"
                for (size_t vi = 0; vi < m_objects[idx].faces[fi].size(); vi++) {             // for each vertex
                    ss << " ";                                                  // write a space
                    m_WriteIndex(ss, m_objects[idx].faces[fi][vi]);                        // write the index
                }
                ss << std::endl;
            }
        }

        /**
         * @brief Write all of the lines associated with an object using the Wavefront OBJ file format
         * 
         * @param ss is an std::stringstream to write the face to
         * @param idx is the index of the object to write
         */
        void m_WriteLines(std::stringstream& ss, size_t idx = 0) {
            for (size_t li = 0; li < m_objects[idx].lines.size(); li++) {                     // for each face
                ss << "l";                                                      // write an "f"
                for (size_t vi = 0; vi < m_objects[idx].lines[li].size(); vi++) {             // for each vertex
                    ss << " ";                                                  // write a space
                    m_WriteIndex(ss, m_objects[idx].lines[li][vi]);                        // write the index
                }
                ss << std::endl;
            }
        }

        /**
         * @brief Writes a Wavefront MTL file containing all materials in the current object.
         * 
         * @param filename is the name of the MTL file to be saved
         */
        void m_SaveMtl(std::string filename) {
            std::stringstream ss;

            for (size_t i = 0; i < m_materials.size(); i++) {
                ss << "newmtl " << m_materials[i].name << std::endl;
                ss << "Ka " << m_materials[i].Ka[0] << " " << m_materials[i].Ka[1] << " " << m_materials[i].Ka[2] << std::endl;
                ss << "Kd " << m_materials[i].Kd[0] << " " << m_materials[i].Kd[1] << " " << m_materials[i].Kd[2] << std::endl;
                ss << "Ks " << m_materials[i].Ks[0] << " " << m_materials[i].Ks[1] << " " << m_materials[i].Ks[2] << std::endl;
                ss << "Ns " << m_materials[i].Ns << std::endl;
                ss << "d " << m_materials[i].d << std::endl;
                ss << std::endl;
            }

            std::ofstream ofs(filename);
            ofs.write(ss.str().c_str(), ss.str().size());
            ofs.close();
        }

    public:

        /**
         * @brief Construct a new OBJ object. This new OBJ structure initializes a single object named "default"
         *          and sets the current active object to this index.
         */
        obj() {
            m_objects.resize(1);
            m_objects[0].name = "default";
            m_objects[0].material = (size_t)-1;
            m_active_object = 0;
        }

        /**
         * @brief Adds a triangle mesh to the current tira::obj structure. The mesh is added as a separate object
         *          that can be assigned a specific material.
         * 
         * @param mesh is the triangular mesh that will be added
         * @param name is the name assigned to the mesh object in the OBJ file
         * @param material_idx is the index of the material in the internal material array ( (size_t)-1 assigns no material)
         */
        void AddMesh(tmesh mesh, std::string name, size_t material_idx = (size_t)-1) {

            // Calculate the offsets required to insert this geometry into the current structure
            const size_t v_offset = m_vertices.size() + 1;
            const size_t n_offset = m_normals.size() + 1;
            const size_t t_offset = m_texcoords.size() + 1;

            // Get the vertex, normal, and texture coordinates from the mesh structure
            std::vector<float> Vertices = mesh.RavelVertices();
            std::vector<float> Normals = mesh.RavelNormals();
            std::vector<float> TexCoords = mesh.RavelTexCoords();

            // The trimesh structure stores a vertex, normal, and texture coordinate for each
            // vertex, so we only have to loop through each of these arrays once.
            for (size_t vi = 0; vi < mesh.NumVertices(); vi++) {                  // for each vertex

                // Allocate space for each coordinate of the vertex and fill that structure with the
                // provided coordinate values from the mesh.
                _vec v(3);
                for (size_t di = 0; di < 3; di++) {
                    v[di] = Vertices[vi * 3 + di];
                }

                // Do this same pre-allocation and assignment for the normals and texture coordinates
                _vec n(3);
                for (size_t di = 0; di < 3; di++) {
                    n[di] = Normals[vi * 3 + di];
                }
                _vec t(3);
                for (size_t di = 0; di < 3; di++) {
                    t[di] = TexCoords[vi * 3 + di];
                }

                // Push the coordinates for each vertex component to the end of the associated vector
                // in the OBJ structure.
                m_vertices.push_back(v);
                m_normals.push_back(n);
                m_texcoords.push_back(t);
            }

            const std::vector<unsigned int> Faces = mesh.Indices();

            const size_t num_triangles = Faces.size() / 3;                  // calculate the number of triangles
            _object new_object;
            new_object.name = name;
            new_object.material = material_idx;
            for (size_t fi = 0; fi < num_triangles; fi++) {             // for each triangle in the trimesh

                // Create an OBJ entry - in this case a triangle that has three vertices
                _primitive obj_triangle(3);


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
                new_object.faces.push_back(obj_triangle);               // push the triangle to the faces array
            }
            m_objects.push_back(new_object);
        }

        /**
         * @brief Adds a triangle mesh to the current tira::obj structure with its own new material. The mesh is added as a separate object
         *          that can be assigned a specific material.
         * 
         * @param mesh is the triangular mesh that will be added
         * @param name is the name assigned to the mesh object in the OBJ file
         * @param material is the material structure defining material properties associated with this mesh
         */
        void AddMesh(tmesh mesh, std::string name, mtl material) {

            m_materials.push_back(material);                            // add the specified material to the m_materials array
            AddMesh(mesh, name, m_materials.size() - 1);                // call the main AddMesh function
        }

        /**
         * @brief Saves the current state of the OBJ structure as a Wavefront OBJ file.
         * 
         * @param filename is the name of the file that this OBJ structure will be saved to
         */
        void Save(const std::string& filename) {
            std::stringstream ss;

            if (m_materials.size() > 0) {
                ss << "mtllib " << filename + ".mtl" << std::endl << std::endl;
                m_SaveMtl(filename + ".mtl");
            }

            // Write all of the vertex entries
            m_WriteCoordinateArray(ss, m_vertices, "v");
            m_WriteCoordinateArray(ss, m_normals, "vn");
            m_WriteCoordinateArray(ss, m_texcoords, "vt");

            // for each object in the OBJ structure, write the object and its associated primitives
            for(size_t i = 0; i < m_objects.size(); i++) {
                ss << std::endl;
                if (m_objects[i].material < (size_t)-1) {       // if the object has an associated material
                    ss << "usemtl " << m_materials[m_objects[i].material].name << std::endl;
                }
                ss << "o " << m_objects[i].name << std::endl;
                m_WriteFaces(ss, i);
                m_WriteLines(ss, i);
            }

            std::ofstream ofs(filename);
            ofs.write(ss.str().c_str(), ss.str().size());
            ofs.close();
        }

        void Load(const std::string& filename) {

            std::ifstream infile(filename);

        }

        
        /**
         * @brief OpenGL immediate mode overload function that begins a new primitive
         * 
         * @param mode primitive to begin (ex. OBJ_LINES, OBJ_TRIANGLES, etc.)
         */
        void Begin(OBJenum mode) {
            if (m_drawing == true)
                throw std::runtime_error("Currently drawing, call End() first");
            m_drawing = true;
            m_drawtype = mode;
        }

        /**
         * @brief OpenGL immediate mode overload function that ends drawing a primitive that was previously started with Begin().
         * 
         */
        void End() {
            switch (m_drawtype) {
            case OBJ_TRIANGLES: {
                if (m_active_primitive.size() % 3 != 0)
                    throw std::runtime_error("OBJ_TRIANGLES must contain a number of vertices divisible by 3");
                const size_t num_triangles = m_active_primitive.size() / 3;
                for (size_t ti = 0; ti < num_triangles; ti++) {
                    _primitive new_triangle(3);
                    for (size_t vi = 0; vi < 3; vi++) {
                        new_triangle[vi] = m_active_primitive[ti * 3 + vi];
                    }
                    m_objects[m_active_object].faces.push_back(new_triangle);
                }
                break;
            }
            case OBJ_LINE_STRIP: {
                if (m_active_primitive.size() < 2)
                    throw std::runtime_error("OBJ_LINE_STRIP must have at least 2 vertices");
                m_objects[m_active_object].lines.push_back(m_active_primitive);
                break;
            }
            default:
                throw std::runtime_error("Primitive type not supported");
                break;
            }

            m_drawing = false;                       // turn off the drawing flag
            m_active_normal = false;                 // deactivate the normal
            m_active_texcoord = false;               // deactivate the texture coordinate
            m_active_primitive.clear();              // clear the active primitive
        }

        /**
         * @brief OpenGL immediate mode overload that sets the vertex normal for any vertices that are specified after the call.
         * 
         * @param x coordinate for the vertex normal
         * @param y coordinate for the vertex normal
         * @param z coordinate for the vertex normal
         */
        void Normal(float x, float y, float z) {
            _vec new_normal = { x, y, z };
            m_normals.push_back(new_normal);
            m_active_normal = true;
        }

        /**
         * @brief OpenGL immediate mode overload that sets the vertex texture coordinate for any vertices that are specified after the call.
         * 
         * @param u coordinate for the vertex texture
         * @param v coordinate for the vertex texture
         * @param w coordinate for the vertex texture
         */
        void TexCoord(float u, float v, float w) {
            _vec new_texcoord = { u, v, w };
            m_texcoords.push_back(new_texcoord);
            m_active_texcoord = true;
        }

        /**
         * @brief OpenGL immediate mode overload that sets the vertex position for any vertices that are specified after the call.
         * 
         * @param x coordinate for the vertex position
         * @param y coordinate for the vertex position
         * @param z coordinate for the vertex position
         */
        void Vertex(float x, float y, float z) {
            _vec new_vertex = { x, y, z };
            m_vertices.push_back(new_vertex);
            m_Vertex();
        }
    };
}
