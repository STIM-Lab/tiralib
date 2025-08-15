#pragma once

#include <vector>
#include <tira/graphics/camera.h>
#include <tira/fiber.h>

#include <tira/graphics/glVertexBuffer.h>
#include <tira/graphics/glVertexBufferLayout.h>
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
    //template <typename VertexAttributeType = float, typename FiberAttributeType = float>
    template<int VertexAttributeNumber = 1, int FiberAttributeNumber = 1>
    class glFibers {

    public:

        // The RenderMode enumeration allows the user to specify how the fibers will be rendered.
        enum RenderMode {LINES, TUBES, BILLBOARDS};
        typedef std::array<float, VertexAttributeNumber> VertexAttributes;
        typedef std::array<float, FiberAttributeNumber> FiberAttributes;

    protected:
        
        /**
         * Describes how the fibers will be rendered using OpenGL. Initial render types will include
         * GL_LINES and tubes. I would like to implement billboards soon.
         */
        RenderMode m_mode;                       
        
        /**
         * Vector that stores all fibers in the structure. This array will be used to build vertex arrays
         * and buffers to facilitate rendering.
         */
        std::vector< fiber< VertexAttributes > > m_fibers;

        /**
         * Pair that stores an aligned bounding box that surrounds the entire set of fibers.
         * The first element is the lower-left (smallest) corner and the second element is the
         * upper-right (largest) corner.
         */
        std::pair<glm::vec3, glm::vec3> m_aabb;

        std::array< std::pair<float, float>, VertexAttributeNumber> m_extrema;

        std::vector<glVertexBuffer> m_vbuffers;

        //glVertexArray _varray;
        glVertexBufferLayout m_layout;
        bool m_buffers_valid;                // flag specifies that the buffers are up to date and accurately represent the fiber data

        inline static const std::string brewer_pts =
        #include "src/vascuvis/brewer_pts.shader"
        ;

        inline static const std::string cmap_vertex_string =
        #include "src/vascuvis/cmap_vertex.shader"
        ;

        inline static const std::string cmap_fragment_string =
        #include "src/vascuvis/cmap_fragment.shader"
        ;

        inline static const std::string fiber_id_string =
        #include "src/vascuvis/fiber_id.shader"
        ;

        glShader m_shader;
        glShader m_fiberid_shader;
        int m_cmap_attribute;                // attribute to be colormapped

        void _assemble_shader(std::string& vertex, std::string& fragment) {

            vertex = "#version 330 core\n";    // shader header and required version
            vertex += brewer_pts;                                           // set of Brewer colormap points
            vertex += "layout(location = 0) in vec4 position;\n";           // layout location for the vertex positions
            vertex += "layout(location = " + std::to_string(m_cmap_attribute + 1) + ") in float value;\n";
            vertex += cmap_vertex_string;

            fragment = "#version 330 core\n";
            fragment += cmap_fragment_string;
        }

        void _generate_cmap_shader() {
            std::string vertex_shader_string, fragment_shader_string;
            _assemble_shader(vertex_shader_string, fragment_shader_string);
            std::cout << vertex_shader_string << std::endl << std::endl << fragment_shader_string << std::endl;
            m_shader.CreateShader(vertex_shader_string, fragment_shader_string);
        }

        void _update_aabb(glm::vec3 p) {
            if (m_aabb.first.x > p.x) m_aabb.first.x = p.x;
            if (m_aabb.first.y > p.y) m_aabb.first.y = p.y;
            if (m_aabb.first.z > p.z) m_aabb.first.z = p.z;
            if (m_aabb.second.x < p.x) m_aabb.second.x = p.x;
            if (m_aabb.second.y < p.y) m_aabb.second.y = p.y;
            if (m_aabb.second.z < p.z) m_aabb.second.z = p.z;
        }

        void _update_extrema(VertexAttributes a) {
            for (size_t ai = 0; ai < VertexAttributeNumber; ai++) {
                if (m_extrema[ai].first > a[ai]) m_extrema[ai].first = a[ai];
                if (m_extrema[ai].second < a[ai]) m_extrema[ai].second = a[ai];
            }
        }

        void _clear_buffers() {
            for (size_t bi = 0; bi < m_vbuffers.size(); bi++) {
                m_vbuffers[bi].Destroy();
            }
        }

        void _clear_extrema() {
            m_aabb.first = glm::vec3(std::numeric_limits<float>::infinity());      // initialize the bounding box to infinity
            m_aabb.second = glm::vec3(-std::numeric_limits<float>::infinity());
            for (size_t ai = 0; ai < VertexAttributeNumber; ai++) {
                m_extrema[ai].first = std::numeric_limits<float>::infinity();
                m_extrema[ai].second = -std::numeric_limits<float>::infinity();
            }
        }

        void _generate_line_buffers() {
            _clear_buffers();                   // destroy existing OpenGL vertex buffers


            for (size_t fi = 0; fi < m_fibers.size(); fi++) {    // iterate through each fiber

                std::vector<glm::vec3> varray;


                glVertexBuffer new_buffer(&m_fibers[fi][0], m_fibers[fi].size() * sizeof(vertex<VertexAttributes>));

                m_vbuffers.push_back(new_buffer);
                //_count.push_back(_fibers[fi].size());

            }
            m_buffers_valid = true;              // the buffers are now valid and can be rendered
        }

        void _generate_buffers() {
            switch (m_mode) {
            case RenderMode::LINES:
                _generate_line_buffers();
                break;
            default:
                throw std::runtime_error("glLines ERROR: only LINES render mode is implemented");
            }
        }

        bool _validate() {
            if (!m_buffers_valid)
                _generate_buffers();

            size_t point_bytes = sizeof(vertex<VertexAttributes>);
            size_t layout_bytes = m_layout.bytes();

            if (layout_bytes != point_bytes)
                throw std::runtime_error("glFibers ERROR: vertex attribute layout doesn't match the size of the specified template type");
            return true;
        }

    public:

        void clear() {
            _clear_buffers();
            m_buffers_valid = false;                                                     // OpenGL buffers have not been allocated or set
            _clear_extrema();
            m_fibers.clear();
            m_vbuffers.clear();
        }

        glFibers() {
            clear();
            m_cmap_attribute = 0;
            m_mode = RenderMode::LINES;                                          // default to rendering lines
            m_layout.Push<float>(3);                                         // push the first attribute (vertex position) as 3xfloat32
            for (size_t ai = 0; ai < VertexAttributeNumber; ai++)               // each additional attribute will be a single float32
                m_layout.Push<float>(1);
            _generate_cmap_shader();
            m_fiberid_shader.CreateShader(fiber_id_string);
        }

        glm::vec3 Center() {
            return (m_aabb.first + m_aabb.second)/2.0f;
        }

        float Length() {
            return glm::length(m_aabb.second - m_aabb.first);
        }

        /**
         * Add a fiber to the glFibers structure based on separate vectors representing the vertices and point attributes.
         * @param vertices
         * @param attributes
         * @return
         */
        size_t AddFiber(std::vector<glm::vec3> positions, std::vector<VertexAttributes> attributes) {

            if (positions.size() != attributes.size())
                throw std::runtime_error("glFibers ERROR: the vectors of positions and attributes have different sizes");

            fiber<VertexAttributes> new_fiber;                                            // generate a new internal fiber

            for (size_t vi = 0; vi < positions.size(); vi++) {           // for each vertex in the vertex list
                VertexAttributes va;
                va = attributes[vi];                                    // assign the supplied attribute
                vertex<VertexAttributes> new_v(positions[vi], va);    // create a new point for the vertex
                new_fiber.push_back(new_v);                             // push the new point into the new fiber
                _update_aabb(new_v);                                    // update the axis aligned bounding box to include the new point
                _update_extrema(va);
            }
            m_fibers.push_back(new_fiber);                               // add the new fiber to the fiber list

            m_buffers_valid = false;                                     // falsify the buffers (since we've added new data)
            return m_fibers.size() - 1;
        }

        void CmapAttribute(int ai) {
            if (ai < 0 || ai >= VertexAttributeNumber)
                throw std::runtime_error("glFibers ERROR: vertex attribute index out of range");
            m_cmap_attribute = ai;
            _generate_cmap_shader();                             // changing attributes requires regenerating the shader
        }

        int CmapAttribute() { return m_cmap_attribute; }

        float Min(size_t attribute) {
            return m_extrema[attribute].first;
        }
        float Max(size_t attribute) {
            return m_extrema[attribute].second;
        }

        void Render(glm::mat4 View, glm::mat4 Proj, std::vector<float> cmap_override = {}) {

            _validate();            // validate that the data structure is ready for rendering

            for (size_t bi = 0; bi < m_vbuffers.size(); bi++) {

                m_vbuffers[bi].Bind();
                m_shader.Bind();
                m_shader.SetUniformMat4f("V", View);
                m_shader.SetUniformMat4f("P", Proj);
                if (cmap_override.size() >= 1)
                    m_shader.SetUniform1f("vmin", cmap_override[0]);
                else
                    m_shader.SetUniform1f("vmin", m_extrema[m_cmap_attribute].first);

                if (cmap_override.size() >= 2)
                    m_shader.SetUniform1f("vmax", cmap_override[1]);
                else
                    m_shader.SetUniform1f("vmax", m_extrema[m_cmap_attribute].second);
                GLERROR(glEnableVertexAttribArray(0));
                GLERROR(glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vertex<VertexAttributes>), (const void*)0));
                m_layout.Bind();
                glDrawArrays(GL_LINE_STRIP, 0, m_fibers[bi].size());
            }
        }

        void RenderFiberID(glm::mat4 View, glm::mat4 Proj) {
            _validate();            // validate that the data structure is ready for rendering

            for (size_t bi = 0; bi < m_vbuffers.size(); bi++) {

                m_vbuffers[bi].Bind();
                m_fiberid_shader.Bind();
                m_fiberid_shader.SetUniformMat4f("V", View);
                m_fiberid_shader.SetUniformMat4f("P", Proj);
                m_fiberid_shader.SetUniform1i("id", bi);


                GLERROR(glEnableVertexAttribArray(0));
                GLERROR(glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vertex<VertexAttributes>), (const void*)0));
                m_layout.Bind();
                glDrawArrays(GL_LINE_STRIP, 0, m_fibers[bi].size());
            }

            // handle ID here:
            size_t hit = 0;
        }




    };
}
