#pragma once

#include "volume.h"
#include "glTexture.h"

namespace tira{

	class glVolume : public volume {

	protected:
		glTexture _texture;

	public:
		void genRGB(unsigned int X = 32, unsigned int Y = 32, unsigned int Z = 32, unsigned int boxes = 1){
			volume::genRGB(X, Y, Z, boxes);
			_texture.AssignImage((unsigned char*)_data, _dims[0], _dims[1], _dims[2], GL_RGB, GL_RGB, GL_UNSIGNED_BYTE);
		}

		void genGrid(unsigned int X = 32, unsigned int Y = 32, unsigned int Z = 32, unsigned int boxes = 1){
			volume::genGrid(X, Y, Z, boxes);
			_texture.AssignImage((unsigned char*)_data, _dims[0], _dims[1], _dims[2], GL_RED, GL_RED, GL_UNSIGNED_BYTE);
		}

		void setData(void* volData, unsigned int X, unsigned int Y, unsigned int Z, GLenum internalFormat, GLenum externalFormat, GLenum externalDataType){
			_texture.AssignImage((unsigned char*)volData, X, Y, Z, internalFormat, externalFormat, externalDataType);
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

	};
}