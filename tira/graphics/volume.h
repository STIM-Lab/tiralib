#pragma once

#include "glm/glm.hpp"
#include "glm/gtx/quaternion.hpp"

#define XOR(a,b) (((a)&&(!(b)))||(((!a))&&(b)))

class volume{

protected:
	void* _data;
	size_t _dims[3];
	glm::vec3 _position;
	glm::vec3 _scale;
	glm::quat _orientation;

	bool grid(float x, float y, float z, unsigned int boxes){
		bool x_box, y_box, z_box;
		float box_size = 1.0f / (float)boxes;
		if((unsigned int)(z/box_size) % 2)
			z_box = false;
		else
			z_box = true;
			
		if((unsigned int)(y/box_size) % 2)
			y_box = false;
		else
			y_box = true;

		if((unsigned int)(x/box_size) % 2)
			x_box = false;
		else
			x_box = true;

		return XOR(XOR(x_box, y_box), z_box);
	}

	void genGrid(unsigned int X = 32, unsigned int Y = 32, unsigned int Z = 32, unsigned int boxes = 1){
		_dims[0] = X;
		_dims[1] = Y;
		_dims[2] = Z;

		float pixel_size[] = {1.0f/X, 1.0f/Y, 1.0f/Z};

		unsigned char* char_data = (unsigned char*)malloc(X * Y * Z * sizeof(unsigned char));							// allocate space for the cube
		size_t idx_z, idx_y, idx;
		float x, y, z;
		bool x_box, y_box, z_box;
		for(size_t zi = 0; zi < Z; zi++){
			idx_z = zi * X * Y;
			z = zi * pixel_size[0];
			
			for(size_t yi = 0; yi < Y; yi++){
				idx_y = idx_z + yi * X;
				y = yi * pixel_size[0];
				
				for(size_t xi = 0; xi < X; xi++){
					idx = idx_y + xi;
					x = xi * pixel_size[0];
					
					if(grid(x, y, z, boxes)){
						char_data[idx] = sqrt(x*x + y*y + z*z) / sqrt(3) * 255;
					}
					else{
						char_data[idx] = 0;
					}
				}
			}
		}
		_data = (void*)char_data;
	}

	/// Generate an RGB cube with the specified resolution
	void genRGB(unsigned int X = 32, unsigned int Y = 32, unsigned int Z = 32, unsigned int boxes = 1){
		_dims[0] = X;
		_dims[1] = Y;
		_dims[2] = Z;

		
		float pixel_size[] = {1.0f/X, 1.0f/Y, 1.0f/Z};

		unsigned char* char_data = (unsigned char*)malloc(X * Y * Z * sizeof(unsigned char)*3);							// allocate space for the cube
		//memset(char_data, 255, X * Y * Z*3);
		size_t idx_z, idx_y, idx;
		float x, y, z;
		unsigned char r;
		unsigned char g;
		unsigned char b;
		bool x_box, y_box, z_box;
		for(size_t zi = 0; zi < Z; zi++){
			idx_z = zi * X * Y * 3;
			z = zi * pixel_size[0];
			
			for(size_t yi = 0; yi < Y; yi++){
				idx_y = idx_z + yi * X * 3;
				y = yi * pixel_size[0];
				
				for(size_t xi = 0; xi < X; xi++){
					idx = idx_y + xi * 3;
					x = xi * pixel_size[0];
					
					if(grid(x, y, z, boxes)){
						char_data[idx + 0] = x * 255;
						char_data[idx + 1] = y * 255;
						char_data[idx + 2] = z * 255;
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
	glm::vec3 getScale(){
		return _scale;
	}
	void setPosition(float x, float y, float z){
		_position = glm::vec3(x, y, z);
	}
	glm::vec3 getPosition(){
		return _position;
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
	glm::vec3 getEulerAngles(){
		return glm::eulerAngles(_orientation);
	}


	// Create a matrix that transforms world-space coordinates to the local [-0.5, 0.5] volume coordinates
	glm::mat4 getMatrix(){
		
		glm::mat4 result(1.0f);
		result = glm::scale(result, glm::vec3(1.0, 1.0, 1.0)/_scale);
		result = glm::translate(result, _position);
		result = result * glm::toMat4(_orientation);

		return result;
	}


};