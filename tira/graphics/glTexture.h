# pragma once

#include <tira/graphics/glGeometry.h>
#include <tira/image.h>


namespace tira {

	//enum glTextureFilterType { TextureFilterLinear, TextureFilterNearest };

	class glTexture {

	protected:
		GLuint m_TextureID;				// texture OpenGL name
		int m_Width, m_Height;			// width and height (2D and 3D textures)
		int m_Depth;					// depth = 0 for a 2D texture
		GLenum m_TextureType;			// GL_TEXTURE_1D, GL_TEXTURE_2D, GL_TEXTURE_3D, etc.
		GLint m_internal_format;		// internal format of the texture on the GPU


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


		/**
		 * Helper function that converts a specified number of color channels into an OpenGL constant
		 * @param c number of color channels
		 * @return OpenGL constant representing that number of channels
		 */
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

		unsigned int _glformat2channels(GLenum transferFormat) {
			switch (transferFormat) {
			case GL_RED:
			case GL_RED_INTEGER:
				return 1;
			case GL_RG:
			case GL_RG_INTEGER:
				return 2;
			case GL_RGB:
			case GL_RGB_INTEGER:
				return 3;
			case GL_RGBA:
			case GL_RGBA_INTEGER:
				return 4;
			default:
				throw std::runtime_error("glTexture ERROR: unable to find the number of channels from the provided transfer format");
			}
		}

		/**
		 * Apparently even empty textures have to be specified with a transfer format that is compatible with the specified internal
		 * format. For example, if you try to create an empty texture with an internal format of GL_R32I, the transfer format has to be
		 * an integral type (even though no data is transferred). This function provides a compatible transfer format given common
		 * internal formats.
		 * @param internalFormat internal format that requires a compatible transfer format type
		 * @return a valid transfer format for the provided internal format
		 */
		GLenum _get_compatible_transfer_format(GLint internalFormat) {
			switch (internalFormat) {
			case GL_R8I:
			case GL_R8UI:
			case GL_R16I:
			case GL_R16UI:
			case GL_R32I:
			case GL_R32UI:
				return GL_RED_INTEGER;
			case GL_RG8I:
			case GL_RG8UI:
			case GL_RG16I:
			case GL_RG16UI:
			case GL_RG32I:
			case GL_RG32UI:
				return GL_RG_INTEGER;
			case GL_RGB8I:
			case GL_RGB8UI:
			case GL_RGB16I:
			case GL_RGB16UI:
			case GL_RGB32I:
			case GL_RGB32UI:
				return GL_RGB_INTEGER;
			case GL_RGBA8I:
			case GL_RGBA8UI:
			case GL_RGBA16I:
			case GL_RGBA16UI:
			case GL_RGBA32I:
			case GL_RGBA32UI:
				return GL_RGBA_INTEGER;
			case GL_RED:
			case GL_R8:
			case GL_R16:
			case GL_R16F:
			case GL_R32F:
				return GL_RED;
			case GL_RG:
			case GL_RG8:
			case GL_RG16:
			case GL_RG16F:
			case GL_RG32F:
				return GL_RG;
			case GL_RGB:
			case GL_RGB8:
			case GL_RGB16:
			case GL_RGB16F:
			case GL_RGB32F:
				return GL_RGB;
			case GL_RGBA:
			case GL_RGBA8:
			case GL_RGBA16:
			case GL_RGBA16F:
			case GL_RGBA32F:
				return GL_RGBA;
			default:
				throw std::runtime_error("glTexture ERROR: unable to find a transfer format compatible with the specified internal format");
			}
		}


	public:


		glTexture() {
			m_Depth = m_Width = m_Height = 0;
			m_TextureType = GL_TEXTURE_2D;
			m_TextureID = 0;
			m_internal_format = 3;
		}

		/**
		 * Construct a glTexture object storing a texture with the specified parameters and data
		 * @param bytes	pointer to the data used to fill the texture
		 * @param width width of the texture in pixels
		 * @param height height of the texture in pixels
		 * @param depth depth of the texture in pixels (0 for a 2D texture)
		 * @param internalFormat format used to store the texture on the GPU (usually the number of color channels but can also be a sized internal format)
		 * @param externalFormat format of the data provided in the bytes array
		 * @param externalDataType data type used to store the data in the bytes array
		 */
		glTexture(const unsigned char* bytes,
		          const int width, const int height, const int depth,
		          const GLint internalFormat,
		          const GLenum externalFormat,
		          const GLenum externalDataType) : glTexture() {
			AssignImage(bytes, width, height, depth, internalFormat, externalFormat,externalDataType);
		}

		/**
		 * Constructor generates an empty texture on the GPU with the specified size and internal format
		 * @param width width of the texture to be created
		 * @param height height of the texture to be created
		 * @param depth depth of a 3D texture (0 for 2D textures)
		 * @param internalFormat OpenGL enumeration describing how the texture is represented on the GPU (see the internalFormat options for glTexImage2D)
		 */
		glTexture(const int width, const int height, const int depth, const GLenum internalFormat) {
			GLenum transferDataType = GL_UNSIGNED_BYTE;					// I think GL_UNSIGNED_BYTE is compatible with any format
			GLenum transferFormat = _get_compatible_transfer_format(internalFormat);

			AssignImage(NULL, width, height, depth, internalFormat, transferFormat, transferDataType);
		}

		/**
		 * Resizes a texture that has already been created.
		 * @param width new width of the texture
		 * @param height new height of the texture
		 * @param depth new depth of the texture (0 will create a 2D texture, any positive integer will create a 3D texture)
		 */
		void Resize(const int width, const int height, const int depth = 0) {
			Bind();
			GLenum transferFormat = _get_compatible_transfer_format(m_internal_format);
			GLenum transferDataType = GL_UNSIGNED_BYTE;			// I think this works in all cases

			AssignImage(NULL, width, height, depth, m_internal_format, transferFormat, transferDataType);
		}

		/**
		 * Main helper function that creates or re-allocates texture memory on the GPU
		 * @param bytes pointer to the data in main memory that will be used to fill the texture (NULL will create a texture without copying data - contents are undefined)
		 * @param width with of the texture to allocate
		 * @param height height of the texture to allocate
		 * @param depth depth of the texture to allocate - 0 creates a 2D texture and any positive integer specifies a 3D texture
		 * @param internalFormat internal format used to store the texture on the GPU
		 * @param externalFormat format for the data stored in bytes used to transfer texture information from main memory
		 * @param externalDataType data type for the contents of bytes
		 */
		void AssignImage(const void* bytes,
		                 int width, int height, int depth,
		                 GLint internalFormat,
		                 GLenum externalFormat,
		                 GLenum externalDataType) {

			m_Width = width;
			m_Height = height;
			m_Depth = depth;

			if (m_TextureID == 0) {
				GLERROR(glGenTextures(1, &m_TextureID));	// generate a new OpenGL texture name
			}

			if (depth == 0) m_TextureType = GL_TEXTURE_2D;	// set the texture type based on input parameters
			else m_TextureType = GL_TEXTURE_3D;


			GLERROR(glBindTexture(m_TextureType, m_TextureID));	// bind the texture
			glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
			GLBREAK(glTexParameteri(m_TextureType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE));
			if (m_TextureType == GL_TEXTURE_3D) {
				GLBREAK(glTexParameteri(m_TextureType, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE));
			}

			if (m_TextureType == GL_TEXTURE_2D) {
				m_internal_format = internalFormat;
				GLBREAK(glTexImage2D(GL_TEXTURE_2D, 0, m_internal_format, m_Width, m_Height, 0,
					externalFormat,
					externalDataType, bytes));
			}
			else {											//if (m_TextureType == GL_TEXTURE_3D)
				m_internal_format = internalFormat;
				GLBREAK(glTexImage3D(GL_TEXTURE_3D, 0, m_internal_format,
					m_Width, m_Height, m_Depth, 0,
					externalFormat,			// Format of the pixel data
					externalDataType, bytes));
			}

			GLBREAK(glBindTexture(m_TextureType, 0));	// unbind the texture
		}

		/**
		 * Deletes a texture and invalidates the texture name
		 */
		void Destroy() const {
			GLBREAK(glDeleteTextures(1, &m_TextureID));
		}

		/**
		 * Destroys the current texture object (doesn't actually do anything to the texture if it is already allocated - use Destroy() for that)
		 */
		~glTexture() {}

		/**
		 * Binds a texture so that it can be used as a shader element
		 */
		void Bind() const {
			GLBREAK(glBindTexture(m_TextureType, m_TextureID));	// bind the texture
		}

		/**
		 * Unbinds a texture - this usually isn't necessary since binding a new texture unbinds previous textures
		 */
		void Unbind() const {
			GLBREAK(glBindTexture(m_TextureType, 0));
		}

		/**
		 * Sets the type of texture filtering used to access the texture
		 * @param filter_type
		 */
		void SetFilter(const GLenum filter_type) const {
			GLBREAK(glBindTexture(m_TextureType, m_TextureID));	// bind the texture
			GLBREAK(glTexParameteri(m_TextureType, GL_TEXTURE_MIN_FILTER, filter_type));
			GLBREAK(glTexParameteri(m_TextureType, GL_TEXTURE_MAG_FILTER, filter_type));
		}

		GLenum internalFormat() {
			GLint format;
			glGetTextureLevelParameteriv(
				m_TextureType,
				0, //e.g. 0
				GL_TEXTURE_INTERNAL_FORMAT,
				&format
			);
			return (GLenum)format;
		}

		/**
		 * Returns the width of the current texture
		 * @return texture width as an int
		 */
		int Width() const { return m_Width; }

		/**
		 * Returns the height of the current texture
		 * @return texture height as an int
		 */
		int Height() const { return m_Height; }

		/**
		 * Returns the depth of the current texture
		 * @return texture depth as an int (returns 0 if the texture is 2D)
		 */
		int Depth() const { return m_Depth; }

		/**
		 * Returns the OpenGL name of the texture
		 * @return texture name (identifier) as an integer (names provided by OpenGL are > 0)
		 */
		GLuint ID() const { return m_TextureID; }

		/**
		 * Tests to see if the current texture exists (if it has a name >0)
		 * @return boolean indicating that the texture is valid (true) or invalid (false)
		 */
		bool valid() {
			if (m_TextureID > 0) return true;
			return false;
		}

		/**
		 * Returns the texture data as a tira::image object
		 * @tparam T data type used to represent the texture
		 * @return the texture represented as a tira::image, useful for saving the texture data as a normal image for debugging
		 */
		template<typename T>
		image<T> Image() {
			Bind();
			GLenum transfer_type = type2gl<T>();											// get the appropriate data type used in the transfer
			GLenum transfer_format = _get_compatible_transfer_format(m_internal_format);	// get the appropriate transfer format
			unsigned int channels = _glformat2channels(transfer_format);

			if (m_Depth != 0) {
				throw std::runtime_error("glTexture ERROR: a 2D image is requested, however the specified texture is 3D");
			}

			image<T> texture_image(m_Width, m_Height, channels);
			glGetTexImage(GL_TEXTURE_2D, 0, transfer_format, transfer_type, texture_image.data());
			return texture_image;
		}
	};

}