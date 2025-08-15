 # pragma once

#include <tira/graphics/glTexture.h>

namespace tira {

    class glFramebuffer : public glTexture{
        protected:

        GLuint m_framebuffer_id;
        std::vector<GLenum> m_draw_buffers;

        public:


        glFramebuffer() {
            m_framebuffer_id = 0;
        }

        glFramebuffer(int width, int height, int format) : glTexture(width, height, 0, format) {
            GLBREAK(glGenFramebuffers(1, &m_framebuffer_id));
            GLBREAK(glBindFramebuffer(GL_FRAMEBUFFER, m_framebuffer_id));
            GLBREAK(glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, m_TextureID, 0));
            m_draw_buffers.push_back(GL_COLOR_ATTACHMENT0);
            GLBREAK(glDrawBuffers(m_draw_buffers.size(), &m_draw_buffers[0]));
        }

        void Bind() {
            glBindFramebuffer(GL_FRAMEBUFFER, m_framebuffer_id);
        }

        void Unbind() {
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
        }




    };


}