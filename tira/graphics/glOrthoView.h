#include <string>

#include "glMaterial.h"
#include "glVolume.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

namespace tira {


    template<typename Type>
    class glOrthoView : public glVolume<Type> {

        std::string _vertexsource =                                  // Source code for the default vertex shader
            "# version 330 core\n"

            "layout(location = 0) in vec3 vertices;\n"
            "layout(location = 2) in vec2 texcoords;\n"

            "uniform mat4 MVP;\n"
            "uniform mat4 T;"

            "out vec4 vertex_color;\n"
            "out vec4 vertex_texcoord;\n"

            "void main() {\n"
            "    gl_Position = MVP * vec4(vertices.x, vertices.y, vertices.z, 1.0f);\n"
            "    vertex_texcoord = T * vec4(texcoords.x, texcoords.y, 0.0, 1.0f);\n"
            "};\n";

        std::string _fragmentsource =
            "# version 330 core\n"

            "layout(location = 0) out vec4 color;\n"

            "in vec4 vertex_color;\n"
            "in vec4 vertex_texcoord;\n"
            "uniform sampler3D texmap;\n"

            "void main() {\n"
            "    color = texture(texmap, vertex_texcoord.xyz);\n"
            "};\n";

    protected:
        glMaterial _material;
        glGeometry _planes[3];
        glm::vec3 _normalized_slices = glm::vec3(0.0f);             // slice positions given in normalized coordinates (0.0, 1.0)
        float _viewport_aspect = 1.0f;
        float _viewport_width = 1.0f;
        GLenum _filter_type = GL_NEAREST;
        GLenum _internal_format = GL_RGB;

        void _update_uniform_viewport_width() {
            glm::vec3 dims = dimensions();
            float max_width = std::max(dims.x, dims.z);
            float max_height = std::max(dims.y, dims.z);
            _viewport_width = max_width;

            if (max_height > _viewport_width / _viewport_aspect)
                _viewport_width = max_height * _viewport_aspect;
        }

    public:

        void init() {

            // create a shader from the hard-coded shader strings (above)
            _material.CreateShader(_vertexsource, _fragmentsource);
            glGetError();

            // upload the current volume as a texture map
            if (this->Size() > 0)
                _material.SetTexture("texmap", this->getTexture());

            glGetError();

            // generate the geometry representing three planes
            _planes[0] = glGeometry::GenerateRectangle<float>();
            _planes[1] = glGeometry::GenerateRectangle<float>();
            _planes[2] = glGeometry::GenerateRectangle<float>();
        }

        glm::vec3 dimensions() {
            glm::vec3 s(this->dx() * this->X(),
                        this->dy() * this->Y(),
                        this->dz() * this->Z());
            return s;
        }

         /**
         * Overloads the glVolume function for generating a dummy RGB cube for testing.
         * @param X is the number of samples along X for the color cube
         * @param Y is the number of samples along Y for the color cube
         * @param Z is the number of samples along Z for the color cube
         * @param boxes number of boxes in a grid (default is 1, producing a solid color cube)
         */
        void generate_rgb(unsigned int X = 32, unsigned int Y = 32, unsigned int Z = 32, unsigned int boxes = 1) {
            glVolume<Type>::generate_rgb(X, Y, Z, boxes);
            this->_texture.SetFilter(GL_NEAREST);
            _material.SetTexture("texmap", this->_texture);
            _update_uniform_viewport_width();
            glGetError();
        }

        /**
         * Overloads the glVolume function for generating a single-channel grid for testing.
         * @param X is the number of samples along X
         * @param Y is the number of samples along Y
         * @param Z is the number of samples along Z
         * @param boxes number of boxes in a grid (default is 1, producing a solid color cube)
         */
        void generate_grid(unsigned int X = 32, unsigned int Y = 32, unsigned int Z = 32, unsigned int boxes = 1) {
            glVolume<Type>::generate_grid(X, Y, Z, boxes);
            _update_uniform_viewport_width();
        }

        /**
         * Returns the dimensions of the corresponding axis slice.
         * @param i is the axis identifier (0 = X, 1 = Y, 2 = Z)
         * @return a glm::vec2 storing the dimensions of a slice along the specified axis
         */
        glm::vec2 axis_shape(int i) {
            if (i == 0) return glm::vec2(this->dz() * this->Z(), this->dy() * this->Y());
            else if (i == 1) return glm::vec2(this->dx() * this->X(), this->dz() * this->Z());
            else if (i == 2) return glm::vec2(this->dx() * this->X(), this->dy() * this->Y());
            else throw std::runtime_error("Axis out of range");
        }

        glm::vec3 slice_positions() {
            glm::vec3 sp = _normalized_slices * dimensions();
            return sp;
        }

        void slice_positions(glm::vec3 s) {
            _normalized_slices = s / dimensions();
        }

        glm::vec3 slice_positions_normalized() {
            return _normalized_slices;
        }
        glm::vec3 slice_positions_normalized(glm::vec3 n) {
            _normalized_slices = n / dimensions();
            return _normalized_slices;
        }

        /**
         * Renders an orthographic slice to a viewport with the specified dimensions. This allows the user
         * to override the currently specified aspect ratio, potentially distorting the slice. The viewport dimensions
         * are specified in the same units as the volume.
         * @param axis is the axis identifier (0 = X, 1 = Y, 2 = Z)
         * @param viewport_dimensions size of the viewport (in volume units) that the slice will be rendered to
         */
        void render_slice(const int axis, glm::vec2 viewport_dimensions) {

            glm::vec2 slice_dimensions = axis_shape(axis);

            const glm::mat4 I(1.0f);
            const glm::mat4 Projection = glm::ortho(-viewport_dimensions[0]/2.0f,
                                                viewport_dimensions[0]/2.0f,
                                                -viewport_dimensions[1]/2.0f,
                                                viewport_dimensions[1]/2.0f);

            glm::mat4 Model = I;
            glm::mat4 Scale = glm::scale(Model, glm::vec3(slice_dimensions[0], slice_dimensions[1], 1));
            glm::mat4 MVP = Projection * Scale;


            glm::mat4 T = I;

            if (axis == 0) {
                T = glm::rotate(T, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                T = glm::translate(T, glm::vec3(-1.0f, 0.0f, _normalized_slices.x));

            }
            if (axis == 1) {
                T = glm::rotate(T, glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
                T = glm::translate(T, glm::vec3(0.0f, 0.0f, -_normalized_slices.y));
            }
            if (axis == 2) {
                T = glm::translate(T, glm::vec3(0.0f, 0.0f, _normalized_slices.z));
            }
            _update_uniform_viewport_width();
            this->_texture.SetFilter(GL_NEAREST);
            _material.Begin();
            _material.SetTexture("texmap", this->getTexture());

            _material.SetUniformMat4f("MVP", MVP);
            _material.SetUniformMat4f("T", T);
            _planes[axis].Draw();

            _material.End();
        }

        /**
         * Set the aspect ratio for the viewport used to render slices. This ensures that a constant viewport size is
         * used for sequentially rendered slices, making the slice sizes comparable in traditional visualizations.
         * These parameters can be overridden by specifying the overriding values in in render_slice(...).
         * @param aspect_ratio is the aspect ratio of the destination viewport
         */
        void aspect(float aspect_ratio) {
            _viewport_aspect = aspect_ratio;
            _update_uniform_viewport_width();
        }

        void render_slice(int slice, float viewport_width, float viewport_aspect) {
            glm::vec2 viewport_dimensions(viewport_width, viewport_width / viewport_aspect);
            render_slice(slice, viewport_dimensions);
        }

        void render_slice(int slice, int viewport_pixel_width, int viewport_pixel_height) {

            float viewport_aspect = (float)viewport_pixel_width / (float)viewport_pixel_height;

            glm::vec2 slice_dimensions = axis_shape(slice);

            float slice_aspect = slice_dimensions[0] / slice_dimensions[1];

            float viewport_width;
            if (viewport_aspect > slice_aspect) viewport_width = slice_dimensions[1] * viewport_aspect;
            else viewport_width = slice_dimensions[0];

            render_slice(slice, viewport_width, viewport_aspect);

        }

        void render_slice(int slice) {
            render_slice(slice, _viewport_width, _viewport_aspect);
        }

    };
};
