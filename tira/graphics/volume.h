#pragma once

#include "glm/glm.hpp"
#include "glm/gtx/quaternion.hpp"

class volume{

protected:
	void* _data;
	size_t _dims[3];
	glm::vec3 _position;
	glm::vec3 _scale;
	glm::quat _orientation;

	/// Generate an RGB cube with the specified resolution
	void genRGB(unsigned int X = 32, unsigned int Y = 32, unsigned int Z = 32, unsigned int boxes = 1){
		_dims[0] = X;
		_dims[1] = Y;
		_dims[2] = Z;

		float box_size = 1.0f / (float)boxes;
		float pixel_size[] = {1.0f/X, 1.0f/Y, 1.0f/Z};

		unsigned char* char_data = (unsigned char*)malloc(X * Y * Z * sizeof(unsigned char)*3);							// allocate space for the cube
		//memset(char_data, 255, X * Y * Z*3);
		size_t idx_z, idx_y, idx;
		float x, y, z;
		unsigned char r;
		unsigned char g;
		unsigned char b;
		bool x_box;
		for(size_t zi = 0; zi < Z; zi++){
			idx_z = zi * X * Y * 3;
			
			for(size_t yi = 0; yi < Y; yi++){
				idx_y = idx_z + yi * X * 3;
				for(size_t xi = 0; xi < X; xi++){
					idx = idx_y + xi * 3;
					x = xi * pixel_size[0];
					if((unsigned int)(x/box_size) % 2)
						x_box = false;
					else
						x_box = true;
					if(x_box){
						r = (double)xi / (double)X * 255;
						g = (double)yi / (double)Y * 255;
						b = (double)zi / (double)Z * 255;
						char_data[idx + 0] = r;
						char_data[idx + 1] = g;
						char_data[idx + 2] = b;
					}
					else{
						char_data[idx + 0] = 0;
						char_data[idx + 1] = 0;
						char_data[idx + 2] = 0;
					}
				}
			}
		}
		_data = (void*)char_data;
	}

public:
	/// <summary>
	/// Default constructor
	/// </summary>
	volume() {
		_data = NULL;
		_dims[0] = 0;
		_dims[1] = 0;
		_dims[2] = 0;
		_position = glm::vec3(0, 0, 0);
		_scale = glm::vec3(1, 1, 1);
		_orientation = glm::quat(1, 0, 0, 0);
	}

	void setScale(float x, float y, float z){
		_scale = glm::vec3(x, y, z);
	}
	void setPosition(float x, float y, float z){
		_position = glm::vec3(x, y, z);
	}
	void rotate(float angle, float x, float y, float z){
		float sa2 = sin(angle / 2.0f);
		float ca2 = cos(angle/2.0f);
		glm::quat rot(ca2, x * sa2, y * sa2, z * sa2);
		_orientation = _orientation * rot;
	}
	void setEulerAngles(float x, float y, float z){
		_orientation = glm::quat(glm::vec3(x, y, z));
	}

	glm::mat4 getAlignmentMatrix(){
		
		// generate an alignment matrix that maps a [-0.5, 0.5] volume space to world space
		glm::mat4 result(1.0f);
		result = glm::scale(result, glm::vec3(1.0, 1.0, 1.0)/_scale);
		result = glm::translate(result, _position);
		result = result * glm::toMat4(_orientation);

		return result;
	}


};