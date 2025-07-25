#pragma once

#include "glErrorHandler.h"

namespace tira {

	/** An OpenGL vertex buffer object (VBO) stores all of the per-vertex data
	 *  in an interleaved format. Generally the coordinates for each point are
	 *  stored first, followed by per-vertex data such as colors and texture
	 *  coordinates. This is a wrapper class that handles the creation,
	 *  destruction, and binding of vertex buffers. It is used in conjunction
	 *  with the glVertexBufferLayout class, which describes the format for the
	 *  data stored in the buffer.
	 * @brief      A wrapper class for an OpenGL Vertex Buffer Object
	 */
	class glVertexBuffer {
	private:
		unsigned int _id;
	public:
		
		/**
		 * @brief      Constructs an object that can hold (but does not contain) an OpenGL vertex buffer.
		 */
		glVertexBuffer() : _id(0) {}

		/**
		 * @brief      Constructs and initializes an OpenGL vertex buffer with the provided data.
		 *
		 * @param[in]  data  A pointer to the data that will be copied to the initialized vertex buffer.
		 * 					 This data is handed over directly to the OpenGL driver (no copy is kept local
		 * 					 to the class).
		 * @param[in]  size  The size of the data in bytes, and also the size used to allocate the buffer
		 * 					 on the GPU.
		 */
		glVertexBuffer(const void* data, unsigned int size) :_id(0) {
			SetBuffer(data, size);
		}

		/** The constructor for this class doesn't really do anything. We can't actually delete the buffer
		 *  because copies of it can be made made regularly when passing the object as a non-reference
		 *  parameter. The buffer can actually be deleted from OpenGL control using the Destroy() function.
		 * @brief      Destructor for the class - doesn't really do anything.
		 */		
		~glVertexBuffer() {}

		/**
		 * @brief      Destroys the buffer and deletes its contents from the GPU.
		 */
		void Destroy() {
			GLERROR(glDeleteBuffers(1, &_id));			// send the command to OpenGL to remove the buffer
			_id = 0;									// set the ID to 0, which indicates an unallocated buffer
		}

		/**
		 * @brief      Ask OpenGL for a new buffer ID (if we don't already have one), and allocate and set the associated data
		 *
		 * @param[in]  data  A pointer to the data that will be copied to the initialized vertex buffer.
		 * 					 This data is handed over directly to the OpenGL driver (no copy is kept local
		 * 					 to the class).
		 * @param[in]  size  The size of the data in bytes, and also the size used to allocate the buffer
		 * 					 on the GPU.
		 */
		void SetBuffer(const void* data, unsigned int size) {
			if (_id == 0) {										// if the buffer is unallocated
				GLERROR(glGenBuffers(1, &_id));					// ask OpenGL for a new buffer ID
				GLERROR(glBindBuffer(GL_ARRAY_BUFFER, _id));    // bind the buffer (one can be bound at a time)
				GLERROR(glBufferData(GL_ARRAY_BUFFER, size, data, GL_STATIC_DRAW));	// allocate and copy the data to the GPU
			}
			else
				throw std::runtime_error("glVertexBuffer has already been allocated - implement code to deal with this!");
		}

		/**
		 * @brief      Bind the buffer to use it for rendering.
		 */
		void Bind() const {
			GLERROR(glBindBuffer(GL_ARRAY_BUFFER, _id));
		}

		/**
		 * @brief      Unbind the buffer (can't think of too many reasons to do this).
		 */
		void Unbind() const {
			GLERROR(glBindBuffer(GL_ARRAY_BUFFER, 0));
		}
	};
}