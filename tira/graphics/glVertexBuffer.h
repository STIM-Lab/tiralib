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
			//GLCALL(glDeleteBuffers(1, &m_BufferID));
		}

		void Destroy() {
			GLERROR(glDeleteBuffers(1, &m_BufferID));
		}

		void SetBuffer(const void* data, unsigned int size) {
			GLERROR(glGenBuffers(1, &m_BufferID));
			GLERROR(glBindBuffer(GL_ARRAY_BUFFER, m_BufferID));                  // What will be drawn
			GLERROR(glBufferData(GL_ARRAY_BUFFER, size, data, GL_STATIC_DRAW));
		}
		void Bind() const {
			GLERROR(glBindBuffer(GL_ARRAY_BUFFER, m_BufferID));
		}
		void Unbind() const {
			GLERROR(glBindBuffer(GL_ARRAY_BUFFER, 0));
		}
	};
}