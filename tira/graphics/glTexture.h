# pragma once

#include "glGeometry.h"

namespace tira {

	//enum glTextureFilterType { TextureFilterLinear, TextureFilterNearest };

	class glTexture {
		GLuint m_TextureID;		// texture OpenGL name
		//const unsigned char* m_LocalBuffer;	// pointer to pixels in CPU memory
		int m_Width, m_Height;			// width and height (2D and 3D textures)
		int m_Depth;					// depth = 0 for a 2D texture
		GLenum m_TextureType;			// GL_TEXTURE_1D, GL_TEXTURE_2D, GL_TEXTURE_3D, etc.

	public:

		/// Helper function returns an OpenGL data type constant given a specified data type in the template
		template <typename D>
		static GLenum type2gl() {
			if (typeid(D) == typeid(char)) return GL_BYTE;
			if (typeid(D) == typeid(unsigned char)) return GL_UNSIGNED_BYTE;
			if (typeid(D) == typeid(short)) return GL_SHORT;
			if (typeid(D) == typeid(unsigned short)) return GL_UNSIGNED_SHORT;
			if (typeid(D) == typeid(int)) return GL_INT;
			if (typeid(D) == typeid(unsigned int)) return GL_UNSIGNED_INT;
			if (typeid(D) == typeid(float)) return GL_FLOAT;
			else {
				std::cout << "ERROR: data type is not currently supported (may require OpenGL extensions)" << std::endl;
				exit(1);
			}
		}

		/// <summary>
		/// Helper function that returns an OpenGL format constant given a specified number of channels
		/// </summary>
		/// <param name="c">Number of color channels</param>
		/// <returns></returns>
		static GLenum format2gl(unsigned int c) {
			if (c == 1) return GL_RED;
			if (c == 2) return GL_RG;
			if (c == 3) return GL_RGB;
			if (c == 4) return GL_RGBA;
			else {
				std::cout << "ERROR: texture format does not currently support this number of channels" << std::endl;
				exit(1);
			}
		}

		glTexture() {
			m_Depth = m_Width = m_Height = 0;
			m_TextureType = GL_TEXTURE_2D;
			m_TextureID = 0;
		}

		void AssignImage(const unsigned char* bytes,
			int width, int height, int depth,
			GLint internalFormat,
			GLenum externalFormat,
			GLenum externalDataType) {

			m_Width = width;
			m_Height = height;
			m_Depth = depth;

			if (m_TextureID == 0) {
				std::cout << "Texture ID: " << m_TextureID << std::endl;
				GLERROR(glGenTextures(1, &m_TextureID));	// generate a new OpenGL texture name

			}

			if (depth == 0) m_TextureType = GL_TEXTURE_2D;	// set the texture type based on input parameters
			else m_TextureType = GL_TEXTURE_3D;


			GLERROR(glBindTexture(m_TextureType, m_TextureID));	// bind the texture
			glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
			GLBREAK(glTexParameteri(m_TextureType, GL_TEXTURE_MIN_FILTER, GL_LINEAR));
			GLBREAK(glTexParameteri(m_TextureType, GL_TEXTURE_MAG_FILTER, GL_LINEAR));
			GLBREAK(glTexParameteri(m_TextureType, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE));
			GLBREAK(glTexParameteri(m_TextureType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE));
			if (m_TextureType == GL_TEXTURE_3D)
				GLBREAK(glTexParameteri(m_TextureType, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE));

			if (m_TextureType == GL_TEXTURE_2D) {
				GLBREAK(glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, m_Width, m_Height, 0,
					externalFormat,
					externalDataType, bytes));
			}
			else if (m_TextureType == GL_TEXTURE_3D) {
				GLBREAK(glTexImage3D(GL_TEXTURE_3D, 0, internalFormat, 
					m_Width, m_Height, m_Depth, 0,
					externalFormat,			// Format of the pixel data
					externalDataType, bytes));
			}

			GLBREAK(glBindTexture(m_TextureType, 0));	// unbind the texture
		}
		glTexture(const unsigned char* bytes,
			int width, int height, int depth,
			GLenum internalFormat,
			GLenum externalFormat,
			GLenum externalDataType) : glTexture() {
			AssignImage(bytes, width, height, depth, internalFormat, externalFormat,externalDataType);
		}										   
			
		

		void Destroy() {
			GLBREAK(glDeleteTextures(1, &m_TextureID));
		}

		~glTexture() {
			
		}
		
		void Bind() const {
			GLBREAK(glBindTexture(m_TextureType, m_TextureID));	// bind the texture
		}
		void Unbind() const {
			GLBREAK(glBindTexture(m_TextureType, 0));
		}

		void SetFilter(GLenum filter_type) {
			GLBREAK(glBindTexture(m_TextureType, m_TextureID));	// bind the texture
			GLBREAK(glTexParameteri(m_TextureType, GL_TEXTURE_MIN_FILTER, filter_type));
			GLBREAK(glTexParameteri(m_TextureType, GL_TEXTURE_MAG_FILTER, filter_type));
		}

		inline int Width() { return m_Width; }
		inline int Height() { return m_Height; }
		inline int Depth() { return m_Depth; }
		inline GLuint ID() { return m_TextureID; }
	};

}