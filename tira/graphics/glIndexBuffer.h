#pragma once
#include "glErrorHandler.h"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"

namespace tira {

	/// <summary>
	/// Defines the variables required to create an OpenGL index buffer.
	/// </summary>
	class glIndexBuffer {
	private:
		unsigned int m_BufferID;
		unsigned int m_Count;			// Number of indices 
	public:
		glIndexBuffer() : m_Count(0), m_BufferID(0) {
			//GLCALL(glGenBuffers(1, &m_BufferID));
		}
		// Number of element count
		glIndexBuffer(const unsigned int* data, unsigned int count) : m_Count(count) {
			//GLCALL(glGenBuffers(1, &m_BufferID));
			SetBuffer(data, count);
		}
		~glIndexBuffer() {
			//GLCALL(glDeleteBuffers(1, &m_BufferID));
		}

		void Destroy() {
			GLERROR(glDeleteBuffers(1, &m_BufferID));
		}
		void SetBuffer(const unsigned int* data, unsigned int count) {
			m_Count = count;
			GLERROR(glGenBuffers(1, &m_BufferID));
			GLERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_BufferID));                  // What will be drawn
			GLERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_Count * sizeof(unsigned int), data, GL_STATIC_DRAW));
		}
		void Bind() const {
			GLERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_BufferID));
		}

		void Unbind() const {
			GLERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));
		}
		inline unsigned int GetCount() const { return m_Count; }
	};
}