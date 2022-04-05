# pragma once

#include "glGeometry.h"
#include "octls_image.h"

namespace tira {

	class glTexture {
		unsigned int m_TextureID;		// texture OpenGL name
		//const unsigned char* m_LocalBuffer;	// pointer to pixels in CPU memory
		int m_Width, m_Height;			// width and height (2D and 3D textures)
		int m_Depth;					// depth = 1 for a 2D texture
		GLenum m_TextureType;			// GL_TEXTURE_1D, GL_TEXTURE_2D, GL_TEXTURE_3D, etc.

	public:
		glTexture() {
			m_Depth = m_Width = m_Height = 0;
			m_TextureType = GL_TEXTURE_2D;

			GLCALL(glGenTextures(1, &m_TextureID));	// generate a new OpenGL texture name
		}

		void AssignImage(const unsigned char* bytes,
			int width, int height, int depth,
			GLenum internalFormat,
			GLenum externalFormat,
			GLenum externalDataType) {

			m_Width = width;
			m_Height = height;
			m_Depth = depth;

			if (depth == 0) m_TextureType = GL_TEXTURE_2D;	// set the texture type based on input parameters
			else m_TextureType = GL_TEXTURE_3D;

			GLCALL(glBindTexture(m_TextureType, m_TextureID));	// bind the texture

			GLCALL(glTexParameteri(m_TextureType, GL_TEXTURE_MIN_FILTER, GL_LINEAR));
			GLCALL(glTexParameteri(m_TextureType, GL_TEXTURE_MAG_FILTER, GL_LINEAR));
			GLCALL(glTexParameteri(m_TextureType, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE));
			GLCALL(glTexParameteri(m_TextureType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE));
			if (m_TextureType == GL_TEXTURE_3D)
				GLCALL(glTexParameteri(m_TextureType, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE));

			if (m_TextureType == GL_TEXTURE_2D) {
				GLCALL(glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, m_Width, m_Height, 0,
					externalFormat,
					externalDataType, bytes));
			}
			else if (m_TextureType == GL_TEXTURE_3D) {
				GLCALL(glTexImage3D(GL_TEXTURE_3D, 0, internalFormat, m_Width, m_Height, m_Depth, 0,
					externalFormat,			// Format of the pixel data
					externalDataType, bytes));
			}

			GLCALL(glBindTexture(m_TextureType, 0));	// unbind the texture
		}
		glTexture(const unsigned char* bytes,
			int width, int height, int depth,
			GLenum internalFormat,
			GLenum externalFormat,
			GLenum externalDataType) {
			GLCALL(glGenTextures(1, &m_TextureID));	// generate a new OpenGL texture name
			AssignImage(bytes, width, height, depth, internalFormat, externalFormat,externalDataType);
		}										   
			
		~glTexture() {
			GLCALL(glDeleteTextures(1, &m_TextureID));
		}
		
		void Bind() const {
			GLCALL(glActiveTexture(GL_TEXTURE0));						// 0-8 or 0-32
			GLCALL(glBindTexture(m_TextureType, m_TextureID));
		}
		void Unbind() const {
			GLCALL(glBindTexture(m_TextureType, 0));
		}

		inline int Width() { return m_Width; }
		inline int Height() { return m_Height; }
		inline int Depth() { return m_Depth; }
	};

}

//namespace tira {
//
//	void glTexture::AssignImage(const unsigned char* bytes,
//		int width, int height, int depth,
//		GLenum internalFormat,
//		GLenum externalFormat,
//		GLenum externalDataType) {
//
//		m_Width = width;
//		m_Height = height;
//		m_Depth = depth;
//
//		if (depth == 0) m_TextureType = GL_TEXTURE_2D;	// set the texture type based on input parameters
//		else m_TextureType = GL_TEXTURE_3D;
//
//		GLCALL(glBindTexture(m_TextureType, m_TextureID));	// bind the texture
//
//		GLCALL(glTexParameteri(m_TextureType, GL_TEXTURE_MIN_FILTER, GL_LINEAR));
//		GLCALL(glTexParameteri(m_TextureType, GL_TEXTURE_MAG_FILTER, GL_LINEAR));
//		GLCALL(glTexParameteri(m_TextureType, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE));
//		GLCALL(glTexParameteri(m_TextureType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE));
//		if (m_TextureType == GL_TEXTURE_3D)
//			GLCALL(glTexParameteri(m_TextureType, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE));
//
//		if (m_TextureType == GL_TEXTURE_2D) {
//			GLCALL(glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, m_Width, m_Height, 0,
//				externalFormat,
//				externalDataType, bytes));
//		}
//		else if (m_TextureType == GL_TEXTURE_3D) {
//			GLCALL(glTexImage3D(GL_TEXTURE_3D, 0, internalFormat, m_Width, m_Height, m_Depth, 0,
//				externalFormat,			// Format of the pixel data
//				externalDataType, bytes));
//		}
//
//		GLCALL(glBindTexture(m_TextureType, 0));	// unbind the texture
//	}
//
//	glTexture::glTexture() {
//		m_Depth = m_Width = m_Height = 0;
//		m_TextureType = GL_TEXTURE_2D;
//
//		GLCALL(glGenTextures(1, &m_TextureID));	// generate a new OpenGL texture name
//	}
//
//	glTexture::glTexture(const unsigned char* bytes,
//		int width, int height, int depth,
//		GLenum internalFormat,
//		GLenum externalFormat,
//		GLenum externalDataType) {
//		GLCALL(glGenTextures(1, &m_TextureID));	// generate a new OpenGL texture name
//		AssignImage(bytes, width, height, depth);
//	}
//
//	glTexture::~glTexture() {
//		GLCALL(glDeleteTextures(1, &m_TextureID));
//	}
//
//	void glTexture::Bind() const {
//		GLCALL(glActiveTexture(GL_TEXTURE0));						// 0-8 or 0-32
//		GLCALL(glBindTexture(m_TextureType, m_TextureID));
//	}
//
//	void glTexture::Unbind() const {
//		GLCALL(glBindTexture(m_TextureType, 0));
//	}
//
//}