# pragma once
#include <vector>
#include "glErrorHandler.h"

namespace tira {
	template<typename T>
	struct identity{ typedef T type; };

	struct glVertexBufferElement {
		unsigned int type;
		unsigned int count;
		unsigned char normalized;

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

	class glVertexBufferLayout {
	private:
		std::vector<glVertexBufferElement> m_Elements;
		unsigned int m_Stride;

		void Push(unsigned int count, identity<float>){
			m_Elements.push_back({ GL_FLOAT, count, GL_FALSE });
			m_Stride += count * glVertexBufferElement::GetSizeOfType(GL_FLOAT);
		}

		void Push(unsigned int count, identity<unsigned int>) {
			m_Elements.push_back({ GL_UNSIGNED_INT, count, GL_FALSE });
			m_Stride += count * glVertexBufferElement::GetSizeOfType(GL_UNSIGNED_INT);
		}

		void Push(unsigned int count, identity<unsigned char>) {
			m_Elements.push_back({ GL_UNSIGNED_BYTE, count, GL_TRUE });
			m_Stride += count * glVertexBufferElement::GetSizeOfType(GL_UNSIGNED_BYTE);
		}
	public:
		glVertexBufferLayout() : m_Stride(0) {};


		template<typename T>
		void Push(unsigned int count){
			Push(count, identity<T>());
		}

		inline const std::vector<glVertexBufferElement> GetElements() const { return m_Elements; }
		inline unsigned int GetStride() const { return m_Stride; }

	};
}