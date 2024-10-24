#pragma once
#include "glErrorHandler.h"

namespace tira {

	/// <summary>
	/// Defines the variables required to create an OpenGL index buffer.
	/// </summary>
	class glIndexBuffer {

		unsigned int _BufferID;
		unsigned int _Count;			// Number of indices
	public:
		glIndexBuffer() : _BufferID(0), _Count(0) { }
		// Number of element count
		glIndexBuffer(const unsigned int* data, const unsigned int count) : _Count(count) {
			SetBuffer(data, count);
		}
		~glIndexBuffer() { }

		void Destroy() const {
			GLERROR(glDeleteBuffers(1, &_BufferID));
		}
		void SetBuffer(const unsigned int* data, const unsigned int count) {
			_Count = count;
			GLERROR(glGenBuffers(1, &_BufferID));
			GLERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _BufferID));                  // What will be drawn
			GLERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER, _Count * sizeof(unsigned int), data, GL_STATIC_DRAW));
		}
		void Bind() const {
			GLERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _BufferID));
		}

		static void Unbind() {
			GLERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));
		}
		inline unsigned int GetCount() const { return _Count; }
	};
}