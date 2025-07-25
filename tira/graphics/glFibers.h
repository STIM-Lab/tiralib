#include <vector>
#include <tira/graphics/camera.h>

#include <tira/graphics/glVertexBuffer.h>
#include <tira/graphics/glVertexBufferLayout.h>
#include <tira/graphics/glVertexArray.h>
#include <tira/graphics/glIndexBuffer.h>
#include <tira/graphics/glShader.h>



namespace tira {

    /**
     * @brief      This class provides the framework for rendering fiber-like structures using OpenGL.
     * 
     * This aims to provide all of the infrastructure for rendering large amounts of fiber-like data such
     * as streamlines, centerlines, medial axes, and networks of interconnected curve-like structures. It provides
     * options for rendering style (ex. lines, tubes, and billboards), and allows the user to specify per-point
     * attributes that can be displayed via the fiber radius or color.
     *
     * @tparam     PointAttributeType  Structure or data type used to specify attributes associated with every point of the fiber.
     * This can be a single attribute (ex. float for radii) or a structure containing multiple data types and variables. The
     * formatting for this structure is specified using a glVertexBufferLayout class.
     */
    template <typename PointAttributeType = float, typename FiberAttributeType = float>
    class glFibers {

        // A single Point consists of (1) a 3D vertex and (2) user-specified per-point attributes.
        // The point forms the basis for a fiber, and additional functions will allow the user to
        // specify which attributes will be used for features (ex. radius, color) when rendering.
        class _Point : public glm::vec3 {
        public:
            PointAttributeType r;                           // stores the attribute associated with this point

            _Point& operator=(const glm::vec3& rhs) {       // use operator= to set coordinates from a vec3
                x = rhs.x;
                y = rhs.y;
                z = rhs.z;
                return *this;
            }
        };


        // A Fiber consists of a vector of Points that are all connected from first to last
        class _Fiber : public std::vector<_Point> {
        public:
            FiberAttributeType attribute;
        };

    public:

        // The RenderMode enumeration allows the user to specify how the fibers will be rendered.
        enum RenderMode {LINES, TUBES, BILLBOARDS};

    protected:
        
        /**
         * Describes how the fibers will be rendered using OpenGL. Initial render types will include
         * GL_LINES and tubes. I would like to implement billboards soon.
         */
        RenderMode _mode;                       
        
        /**
         * Vector that stores all fibers in the structure. This array will be used to build vertex arrays
         * and buffers to facilitate rendering.
         */
        std::vector< _Fiber > _fibers;

        /**
         * Pair that stores an aligned bounding box that surrounds the entire set of fibers.
         * The first element is the lower-left (smallest) corner and the second element is the
         * upper-right (largest) corner.
         */
        std::pair<glm::vec3, glm::vec3> _aabb;

        std::vector<glVertexBuffer> _vbuffers;

        //glVertexArray _varray;
        glVertexBufferLayout _layout;
        bool _buffers_valid;                // flag specifies that the buffers are up to date and accurately represent the fiber data

        glShader _shader;



        void _update_aabb(glm::vec3 p) {
            if (_aabb.first.x > p.x) _aabb.first.x = p.x;
            if (_aabb.first.y > p.y) _aabb.first.y = p.y;
            if (_aabb.first.z > p.z) _aabb.first.z = p.z;
            if (_aabb.second.x < p.x) _aabb.second.x = p.x;
            if (_aabb.second.y < p.y) _aabb.second.y = p.y;
            if (_aabb.second.z < p.z) _aabb.second.z = p.z;
        }

        void _clear_buffers() {
            for (size_t bi = 0; bi < _vbuffers.size(); bi++) {
                _vbuffers[bi].Destroy();
            }
        }

        void _generate_line_buffers() {
            _clear_buffers();                   // destroy existing OpenGL vertex buffers


            for (size_t fi = 0; fi < _fibers.size(); fi++) {    // iterate through each fiber

                std::vector<glm::vec3> varray;


                glVertexBuffer new_buffer(&_fibers[fi][0], _fibers[fi].size() * sizeof(_Point));

                _vbuffers.push_back(new_buffer);
                //_count.push_back(_fibers[fi].size());

            }
            _buffers_valid = true;              // the buffers are now valid and can be rendered
        }

        void _generate_buffers() {
            switch (_mode) {
            case RenderMode::LINES:
                _generate_line_buffers();
                break;
            default:
                throw std::runtime_error("glLines ERROR: only LINES render mode is implemented");
            }
        }

        bool _validate() {
            if (!_buffers_valid)
                _generate_buffers();

            size_t point_bytes = sizeof(_Point);
            size_t layout_bytes = _layout.bytes();

            if (layout_bytes != point_bytes)
                throw std::runtime_error("glFibers ERROR: vertex attribute layout doesn't match the size of the specified template type");
            return true;
        }

    public:

        void clear() {
            _clear_buffers();
            _buffers_valid = false;                                                     // OpenGL buffers have not been allocated or set
            _aabb.first = glm::vec3(std::numeric_limits<float>::infinity());      // initialize the bounding box to infinity
            _aabb.second = glm::vec3(-std::numeric_limits<float>::infinity());
            _fibers.clear();
            _vbuffers.clear();
        }

        glFibers() {
            clear();
            _mode = RenderMode::LINES;                                          // default to rendering lines
            _shader.CreateShader(glShaderStrings::vf_brewer);
        }

        glm::vec3 center() {
            return (_aabb.first + _aabb.second)/2.0f;
        }

        float length() {
            return glm::length(_aabb.second - _aabb.first);
        }

        /**
         * Add a fiber to the glFibers structure based on separate vectors representing the vertices and point attributes.
         * @param vertices
         * @param attributes
         * @return
         */
        size_t add_fiber(std::vector<glm::vec3> vertices, std::vector<PointAttributeType> attributes = {0}) {

            _Fiber new_fiber;                                       // generate a new internal fiber

            for (size_t vi = 0; vi < vertices.size(); vi++) {       // for each vertex in the vertex list
                _Point new_p;                                       // create a new point for the vertex
                new_p = vertices[vi];                               // set the coordinates of the point to the vertex coordinates
                if (vi >= attributes.size())                        // if there are more vertices than attributes
                    new_p.r = attributes.back();                    // use the last radius in the attributes vector
                else
                    new_p.r = attributes[vi];                       // otherwise assign the supplied attribute
                new_fiber.push_back(new_p);                         // push the new point into the new fiber
                _update_aabb(new_p);                                // update the axis aligned bounding box to include the new point
            }
            _fibers.push_back(new_fiber);                           // add the new fiber to the fiber list

            _buffers_valid = false;                                 // falsify the buffers (since we've added new data)
            return _fibers.size() - 1;
        }

        void render(glm::mat4 View, glm::mat4 Proj) {

            _validate();            // validate that the data structure is ready for rendering

            for (size_t bi = 0; bi < _vbuffers.size(); bi++) {
                
                _vbuffers[bi].Bind();
                _shader.Bind();
                _shader.SetUniformMat4f("V", View);
                _shader.SetUniformMat4f("P", Proj);
                GLERROR(glEnableVertexAttribArray(0));
                GLERROR(glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(_Point), (const void*)0));
                _layout.Bind();
                glDrawArrays(GL_LINE_STRIP, 0, _fibers[bi].size());
            }
        }

        void vertex_attributes(std::vector<GLenum> types) {
            _layout = glVertexBufferLayout();
            _layout.Push<float>(3);

            for (size_t i = 0; i < types.size(); i++)
                _layout.Push(1, types[i]);            
        }


    };
}