 # pragma once

#include <tira/graphics/glTexture.h>

namespace tira {
    /**
     * Class extends the glTexture class into an OpenGL frame buffer to facilitate off-screen rendering to textures. The
     * constructor takes the size of the frame buffer as well as the desired format, initializing an empty texture that
     * can be used as a render target for OpenGL calls. Using the frame buffer as a render target should just require binding
     * it using the Bind() method. Calling Unbind() will return to default behavior (ex. the default framebuffer using GLFW).
     */
    class glFramebuffer : public glTexture{
        protected:

        GLuint m_framebuffer_id;
        std::vector<GLenum> m_draw_buffers;

        public:
        /**
         * The default constructor doesn't create a frame buffer - it just sets the ID to zero so that we know no buffer is created
         */
        glFramebuffer() {
            m_framebuffer_id = 0;
        }

        /**
         * Constructor takes the frame buffer size and the internal storage format (see the options for internalFormat in
         * the glTexImage2D function.
         * @param width width of the frame buffer and associated texture
         * @param height height of the frame buffer and associated texture
         * @param format internal format used to represent the frame buffer and associated texture on the GPU (see the internalFormat options for glTexImage2D)
         */
        glFramebuffer(int width, int height, int format) : glTexture(width, height, 0, format) {
            GLBREAK(glGenFramebuffers(1, &m_framebuffer_id));
            GLBREAK(glBindFramebuffer(GL_FRAMEBUFFER, m_framebuffer_id));
            GLBREAK(glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_TextureID, 0));
            m_draw_buffers.push_back(GL_COLOR_ATTACHMENT0);
            GLBREAK(glDrawBuffers(m_draw_buffers.size(), &m_draw_buffers[0]));
        }

        /**
         * Make this frame buffer the current render target for all OpenGL render calls
         */
        void Bind() {
            glBindFramebuffer(GL_FRAMEBUFFER, m_framebuffer_id);
        }

        /**
         * Return to the default rendering target (ex. the default GLFW framebuffer)
         */
        void Unbind() {
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
        }

    };


}