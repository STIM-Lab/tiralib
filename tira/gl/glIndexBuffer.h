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
		glIndexBuffer() : m_Count(0) {
			GLCALL(glGenBuffers(1, &m_BufferID));
		}
		// Number of element count
		glIndexBuffer(const unsigned int* data, unsigned int count) : m_Count(count) {
			GLCALL(glGenBuffers(1, &m_BufferID));
			SetBuffer(data, count);
		}
		~glIndexBuffer() {
			GLCALL(glDeleteBuffers(1, &m_BufferID));
		}
		void SetBuffer(const unsigned int* data, unsigned int count) {
			m_Count = count;
			GLCALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_BufferID));                  // What will be drawn
			GLCALL(glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_Count * sizeof(unsigned int), data, GL_STATIC_DRAW));
		}
		void Bind() const {
			GLCALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_BufferID));
		}

		void Unbind() const {
			GLCALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));
		}
		inline unsigned int GetCount() const { return m_Count; }
	};
}
//namespace tira {
//	glIndexBuffer::glIndexBuffer() : m_Count(0) {
//		GLCALL(glGenBuffers(1, &m_BufferID));
//	}
//	glIndexBuffer::glIndexBuffer(const unsigned int* data, unsigned int count) : m_Count(count) {
//		GLCALL(glGenBuffers(1, &m_BufferID));
//		SetBuffer(data, count);
//	}
//
//	// Delete the IndexBuffer
//	glIndexBuffer::~glIndexBuffer() {
//		GLCALL(glDeleteBuffers(1, &m_BufferID));
//	}
//
//	void glIndexBuffer::SetBuffer(const unsigned int* data, unsigned int count) {
//		m_Count = count;
//		GLCALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_BufferID));                  // What will be drawn
//		GLCALL(glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_Count * sizeof(unsigned int), data, GL_STATIC_DRAW));
//	}
//
//	void glIndexBuffer::Bind() const {
//		GLCALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_BufferID));
//	}
//
//	void glIndexBuffer::Unbind() const {
//		GLCALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));
//	}
//}