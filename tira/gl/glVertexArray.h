# pragma once
#include "glVertexBuffer.h"
#include "glErrorHandler.h"

namespace tira {
	class glVertexBufferLayout;

	class glVertexArray {
	private:
		unsigned int m_ArrayID;
	public:
		glVertexArray() {
			GLCALL(glGenVertexArrays(1, &m_ArrayID));
		}
		~glVertexArray() {
			GLCALL(glDeleteVertexArrays(1, &m_ArrayID))
		}

		void AddBuffer(const glVertexBuffer& vb, const glVertexBufferLayout& layout, unsigned long offset) {
			Bind();
			vb.Bind();
			const auto& elements = layout.GetElements();

			// Set Layout
			for (unsigned int i = 0; i < elements.size(); i++) {
				const auto& element = elements[i];
				GLCALL(glEnableVertexAttribArray(i));
				GLCALL(glVertexAttribPointer(i, element.count, element.type, element.normalized, layout.GetStride(), (const void*)offset));
				offset += element.count * glVertexBufferElement::GetSizeOfType(element.type);
			}
		}

		void Bind() const {
			GLCALL(glBindVertexArray(m_ArrayID));
		}

		void Unbind() const {
			GLCALL(glBindVertexArray(0));
		}
	};
}
//
//namespace tira {
//	glVertexArray::glVertexArray() {
//		GLCALL(glGenVertexArrays(1, &m_ArrayID));
//	}
//
//	glVertexArray::~glVertexArray() {
//		GLCALL(glDeleteVertexArrays(1, &m_ArrayID))
//	}
//
//	void glVertexArray::AddBuffer(const glVertexBuffer& vb, const glVertexBufferLayout& layout, unsigned long offset) {
//		Bind();
//		vb.Bind();
//		const auto& elements = layout.GetElements();
//
//		// Set Layout
//		for (unsigned int i = 0; i < elements.size(); i++) {
//			const auto& element = elements[i];
//			GLCALL(glEnableVertexAttribArray(i));
//			GLCALL(glVertexAttribPointer(i, element.count, element.type, element.normalized, layout.GetStride(), (const void*)offset));
//			offset += element.count * glVertexBufferElement::GetSizeOfType(element.type);
//		}
//	}
//
//	void glVertexArray::Bind() const {
//		GLCALL(glBindVertexArray(m_ArrayID));
//	}
//
//	void glVertexArray::Unbind() const {
//		GLCALL(glBindVertexArray(0));
//	}
//}