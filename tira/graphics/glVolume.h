#pragma once

#include "volume.h"
#include "glTexture.h"

namespace tira{

	class glVolume : public volume {

	protected:
		glTexture _texture;

	public:
		void genRGB(unsigned int X = 32, unsigned int Y = 32, unsigned int Z = 32){
			volume::genRGB(X, Y, Z);
			_texture.AssignImage((unsigned char*)_data, _dims[0], _dims[1], _dims[2], GL_RGB, GL_RGB, GL_UNSIGNED_BYTE);
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