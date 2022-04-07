#pragma once

#include "glErrorHandler.h"

namespace tira {

	/// <summary>
	/// Variables required to build and enable an OpenGL vertex buffer
	/// </summary>
	class glVertexBuffer {
	private:
		unsigned int m_BufferID;
	public:
		/// <summary>
		/// Constructor creates an object storing individual OpenGL vertex coordinates
		/// </summary>
		/// <param name="data">Pointer to a CPU-based array containing vertex information</param>
		/// <param name="size">Size of the vertex buffer (in bytes)</param>
		glVertexBuffer() : m_BufferID(0) {}
		// size usually bytes
		glVertexBuffer(const void* data, unsigned int size) {
			SetBuffer(data, size);
		}
		~glVertexBuffer() {
			GLCALL(glDeleteBuffers(1, &m_BufferID));
		}

		void SetBuffer(const void* data, unsigned int size) {
			GLCALL(glGenBuffers(1, &m_BufferID));
			GLCALL(glBindBuffer(GL_ARRAY_BUFFER, m_BufferID));                  // What will be drawn
			GLCALL(glBufferData(GL_ARRAY_BUFFER, size, data, GL_STATIC_DRAW));
		}
		void Bind() const {
			GLCALL(glBindBuffer(GL_ARRAY_BUFFER, m_BufferID));
		}
		void Unbind() const {
			GLCALL(glBindBuffer(GL_ARRAY_BUFFER, 0));
		}
	};
}

//namespace tira {
//	glVertexBuffer::glVertexBuffer() : m_BufferID(0) {}
//	glVertexBuffer::glVertexBuffer(const void* data, unsigned int size) {
//		SetBuffer(data, size);
//	}
//
//	void glVertexBuffer::SetBuffer(const void* data, unsigned int size) {
//		GLCALL(glGenBuffers(1, &m_BufferID));
//		GLCALL(glBindBuffer(GL_ARRAY_BUFFER, m_BufferID));                  // What will be drawn
//		GLCALL(glBufferData(GL_ARRAY_BUFFER, size, data, GL_STATIC_DRAW));
//	}
//
//	// Delete the VertexBuffer
//	glVertexBuffer::~glVertexBuffer() {
//		GLCALL(glDeleteBuffers(1, &m_BufferID));
//	}
//
//	void glVertexBuffer::Bind() const {
//		GLCALL(glBindBuffer(GL_ARRAY_BUFFER, m_BufferID));
//	}
//
//	void glVertexBuffer::Unbind() const {
//		GLCALL(glBindBuffer(GL_ARRAY_BUFFER, 0));
//	}
//}