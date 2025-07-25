# pragma once
#include <vector>
#include "glErrorHandler.h"

namespace tira {

	/**
	 * @brief      Temporary parameter used to identify a data type and convert it to an OpenGL enumeration
	 *
	 * @tparam     T     Data type to be specified
	 */
	template<typename T>
	struct identity{ typedef T type; };

	/**
	 * @brief      Specifies a single vertex buffer element (position, color, etc.)
	 */
	struct glVertexBufferElement {

		/**
		 * OpenGL enumeration representing the data type for this element.
		 */
		unsigned int type;

		/**
		 *  Number of components of the associated data type in this element. For example
		 *  a 3D vertex will have a count of 3.
		 */
		unsigned int count;

		/**
		 *  Flag indicates whether or not the element is normalized. This is a parameter
		 *  required by OpenGL, but will almost always be false;
		 */
		unsigned char normalized;

		/**
		 * @brief      Returns the size of the data type associated with an OpenGL enumeration.
		 *
		 * @param[in]  type  OpenGL enumeration representing the data type
		 *
		 * @return     The size of the type.
		 */
		static unsigned int GetSizeOfType(unsigned int type) {
			switch (type) {
			case GL_FLOAT: return 4;
			case GL_UNSIGNED_INT: return 4;
			case GL_UNSIGNED_BYTE: return 1;
			}
			ASSERT(false);
			return 0;
		}
	};

	/**
	 * @brief      This class stores the data describing the layout of a vertex buffer.
	 */
	class glVertexBufferLayout {
	private:

		/**
		 *  Vector that stores the individual elements of a vertex buffer.
		 */
		std::vector<glVertexBufferElement> _elements;

		/**
		 *  Stores the size of a vertex (in bytes) based on the current layout.
		 */
		unsigned int _stride;

		/**
		 * @brief      Adds a 32-bit floating point element to the layout.
		 *
		 * @param[in]  count      Number of floating point values in this element
		 * @param[in]  <unnamed>  identity<float> parameter flags this as a floating point element
		 */
		void Push(unsigned int count, identity<float>){
			_elements.push_back({ GL_FLOAT, count, GL_FALSE });
			_stride += count * glVertexBufferElement::GetSizeOfType(GL_FLOAT);
		}

		/**
		 * @brief      Adds a 32-bit unsigned integer element to the layout.
		 *
		 * @param[in]  count      Number of unsigned int values in this element
		 * @param[in]  <unnamed>  identity<unsigned int> parameter flags this as an unsigned int element
		 */
		void Push(unsigned int count, identity<unsigned int>) {
			_elements.push_back({ GL_UNSIGNED_INT, count, GL_FALSE });
			_stride += count * glVertexBufferElement::GetSizeOfType(GL_UNSIGNED_INT);
		}

		/**
		 * @brief      Adds an 8-bit unsigned integer element to the layout.
		 *
		 * @param[in]  count      Number of unsigned char values in this element
		 * @param[in]  <unnamed>  identity<unsigned char> parameter flags this as an unsigned char element
		 */
		void Push(unsigned int count, identity<unsigned char>) {
			_elements.push_back({ GL_UNSIGNED_BYTE, count, GL_TRUE });
			_stride += count * glVertexBufferElement::GetSizeOfType(GL_UNSIGNED_BYTE);
		}
	public:

		/**
		 * @brief      Constructor initializes the stride to zero (a vertex with this layout starts with no data)
		 */
		glVertexBufferLayout() : _stride(0) {};


		/**
		 * @brief      Push an element into the current layout by specifying the data type and number of values
		 *
		 * @param[in]  count  Number of values in the element (ex. a 3D vertex will have a count of 3)
		 *
		 * @tparam     T      Dummy variable specifies the data type as a parameter in order to call the relevant Push sub-function
		 */
		template<typename T>
		void Push(unsigned int count){
			Push(count, identity<T>());
		}

		/**
		 * @brief      Push an element into the current layout by specifying its OpenGL enumeration type (ex. GL_FLOAT)
		 *
		 * @param[in]  count  Number of values in the element (ex. a 3D vertex will have a count of 3)
		 *
		 * @param[in]  type	  OpenGL enumeration type
		 */
		void Push(unsigned int count, GLenum type) {
			_elements.push_back({ type, count, GL_FALSE });
			_stride += count * glVertexBufferElement::GetSizeOfType(type);
		}

		/**
		 * @brief      Returns a vector containing all elements in the layout
		 *
		 * @return     Vector of elements (generally passed to a vertex buffer)
		 */
		inline const std::vector<glVertexBufferElement> GetElements() const { return _elements; }

		/**
		 * @brief      Returns the size (in bytes) of a single vertex using this layout
		 *
		 * @return     The number of bytes in a relevant vertex
		 */
		inline unsigned int GetStride() const { return _stride; }

		/**
		 * @brief      Bind the layout before rendering the target vertex buffer object
		 */
		void Bind(){

			// Set Layout
			size_t offset = 0;
			for (unsigned int i = 0; i < _elements.size(); i++) {
				const auto& element = _elements[i];
				GLERROR(glEnableVertexAttribArray(i));
				GLERROR(glVertexAttribPointer(i, element.count, element.type, element.normalized, _stride, (const void*)offset));
				offset += element.count * glVertexBufferElement::GetSizeOfType(element.type);
			}
		}

		size_t bytes() { return _stride; }

	};
}