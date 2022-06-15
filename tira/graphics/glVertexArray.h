# pragma once
#include "glVertexBuffer.h"
#include "glErrorHandler.h"

namespace tira {
	class glVertexBufferLayout;

	class glVertexArray {
	private:
		unsigned int m_ArrayID;
	public:
		glVertexArray() : m_ArrayID(0) {
			//GLCALL(glGenVertexArrays(1, &m_ArrayID));
		}
		~glVertexArray() {
			//GLCALL(glDeleteVertexArrays(1, &m_ArrayID))
		}

		void Destroy() {
			GLERROR(glDeleteVertexArrays(1, &m_ArrayID));
		}

		void AddBuffer(const glVertexBuffer& vb, const glVertexBufferLayout& layout, size_t offset) {
			GLERROR(glGenVertexArrays(1, &m_ArrayID));
			Bind();
			vb.Bind();
			const auto& elements = layout.GetElements();

			// Set Layout
			for (unsigned int i = 0; i < elements.size(); i++) {
				const auto& element = elements[i];
				GLERROR(glEnableVertexAttribArray(i));
				GLERROR(glVertexAttribPointer(i, element.count, element.type, element.normalized, layout.GetStride(), (const void*)offset));
				offset += element.count * glVertexBufferElement::GetSizeOfType(element.type);
			}
		}

		void Bind() const {
			GLERROR(glBindVertexArray(m_ArrayID));
		}

		void Unbind() const {
			GLERROR(glBindVertexArray(0));
		}
	};
}