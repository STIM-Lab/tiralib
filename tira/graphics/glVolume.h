#pragma once

#include "../volume.h"
#include "glTexture.h"

namespace tira{

	template<class T>
	class glVolume : public volume<T> {

	protected:
		glTexture _texture;

	public:
		void generate_rgb(unsigned int X = 32, unsigned int Y = 32, unsigned int Z = 32, unsigned int boxes = 1){
			volume<T>::generate_rgb(X, Y, Z, boxes);
			GLint external_type = glTexture::type2gl<T>();
			_texture.AssignImage((unsigned char*)&volume<T>::m_data[0], X, Y, Z, GL_RGB, GL_RGB, external_type);
		}

		void generate_grid(unsigned int X = 32, unsigned int Y = 32, unsigned int Z = 32, unsigned int boxes = 1){
			volume<T>::generate_grid(X, Y, Z, boxes);
			GLint external_type = glTexture::type2gl<T>();
			_texture.AssignImage((unsigned char*)&volume<T>::m_data[0], X, Y, Z, GL_RED, GL_RED, external_type);
		}

		void setData(void* volData, unsigned int X, unsigned int Y, unsigned int Z, GLint internalFormat, GLenum externalFormat, GLenum externalDataType){
			_texture.AssignImage((unsigned char*)volData, X, Y, Z, internalFormat, externalFormat, externalDataType);
		}

		void update_texture() {
			GLenum internalFormat = glTexture::format2gl(this->C());
			GLenum externalFormat = glTexture::format2gl(this->C());
			GLenum externalDataType = glTexture::type2gl<T>();
			_texture.AssignImage((unsigned char*)this->Data(), this->X(), this->Y(), this->Z(), internalFormat, externalFormat, externalDataType);
		}

		template<typename D = T>
		void load_npy(std::string filename)
		{
			volume<T>::template load_npy<D>(filename);
			size_t xn = volume<T>::X();
			size_t yn = volume<T>::Y();
			size_t zn = volume<T>::Z();
			size_t cn = volume<T>::C();

			GLint internalFormat = glTexture::format2gl(cn);
			GLenum externalFormat = glTexture::format2gl(cn);
			GLenum externalDataType = glTexture::type2gl<T>();
			setData((unsigned char*)&volume<T>::_data[0], xn, yn, zn, internalFormat, externalFormat, externalDataType);
		}

		void setFilter(GLenum filter_type){
			_texture.SetFilter(filter_type);
		}

		GLuint getTextureID(){
			return _texture.ID();
		}

		void Bind(){
			_texture.Bind();
		}

		void Unbind(){
			_texture.Unbind();
		}

		glTexture getTexture() {
			return _texture;
		}

	};
}