#pragma once

//#include <stim/math/vector.h>
//#include <stim/math/quaternion.h>
//#include <stim/math/matrix_sq.h>
#include "glm/glm.hpp"
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/quaternion.hpp>

#include <ostream>
#include <iostream>

namespace tira{

class camera
{
	glm::vec3 d;
	glm::vec3 p;
	glm::vec3 up;
	float focus;		//focal length of the camera
	float fov;

	//private function makes sure that the up vector is orthogonal to the direction vector and both are normalized
	void stabalize()
	{
		glm::vec3 side = glm::cross(up, d);
		up = glm::cross(d, side);
		up = glm::normalize(up);
		d = glm::normalize(d);
	}

public:
	void setPosition(glm::vec3 pos)
	{
		p = pos;
	}
	void setPosition(float x, float y, float z){setPosition(glm::vec3(x, y, z));}

	void setFocalDistance(float distance){focus = distance;}
	void setFOV(float field_of_view){fov = field_of_view;}

	void LookAt(glm::vec3 pos)
	{
		//find the new direction
		d = pos - p;

		//find the distance from the look-at point to the current position
		focus = glm::length(d);

		//stabalize the camera
		stabalize();
	}
	void LookAt(float px, float py, float pz){LookAt(glm::vec3(px, py, pz));}
	void LookAt(glm::vec3 pos, glm::vec3 new_up){up = new_up; LookAt(pos);}
	void LookAt(float px, float py, float pz, float ux, float uy, float uz){LookAt(glm::vec3(px, py, pz), glm::vec3(ux, uy, uz));}
	void LookAtDolly(float lx, float ly, float lz)
	{
		//find the current focus point
		glm::vec3 f = p + focus*d;
		glm::vec3 T = glm::vec3(lx, ly, lz) - f;
		p = p + T;
	}

	void Dolly(glm::vec3 direction)
	{
		p = p+direction;
	}
	void Dolly(float x, float y, float z){Dolly(glm::vec3(x, y, z));}

	/// Push the camera along the view direction
	void Push(float delta)
	{
		if(delta > focus)
			delta = focus;

		focus -= delta;

		Dolly(d*delta);
	}

	void Pan(float theta_x, float theta_y, float theta_z)
	{
		

		
		glm::quat qx = glm::angleAxis(theta_x, up);			//x rotation is around the up axis

		
		glm::vec3 side = glm::cross(up, d);					//y rotation is around the side axis
		glm::quat qy = glm::angleAxis(theta_y, side);

		glm::quat qz = glm::angleAxis(theta_z, d);			//z rotation is around the direction vector

		
		glm::quat final = qz * qy * qx;						//combine the rotations in x, y, z order

		
		glm::mat3 rotation = glm::mat3_cast(final);			//get the rotation matrix

		
		d = rotation*d;										//apply the rotation
		up = rotation*up;

		
		stabalize();										//stabalize the camera

	}
	void Pan(float theta_x){Pan(theta_x, 0, 0);}
	void Tilt(float theta_y){Pan(0, theta_y, 0);}
	void Twist(float theta_z){Pan(0, 0, theta_z);}

	void Zoom(float delta)
	{
		fov -= delta;
		if(fov < 0.5)
			fov = 0.5;
		if(fov > 180)
			fov = 180;
	}

	void OrbitFocus(float theta_x, float theta_y)
	{
		//find the focal point
		glm::vec3 focal_point = p + focus*d;

		//center the coordinate system on the focal point
		glm::vec3 centered = p - (focal_point - glm::vec3(0, 0, 0));

		//create the x rotation (around the up vector)
		glm::quat qx = glm::angleAxis(theta_x, up);
		centered = glm::vec3(0, 0, 0) + glm::mat3_cast(qx)*(centered - glm::vec3(0, 0, 0));

		//get a side vector for theta_y rotation
		glm::vec3 side = glm::normalize(glm::cross(up, glm::vec3(0, 0, 0) - centered));
		glm::quat qy = glm::angleAxis(theta_y, side);
		centered = glm::vec3(0, 0, 0) + glm::mat3_cast(qy)*(centered - glm::vec3(0, 0, 0));

		//re-position the camera
		p = centered + (focal_point - glm::vec3(0, 0, 0));

		//make sure we are looking at the focal point
		LookAt(focal_point);

		//stabalize the camera
		stabalize();

	}

	void Slide(float u, float v)
	{
		glm::vec3 V = glm::normalize(up);
		glm::vec3 U = glm::normalize(glm::cross(up, d));

		p = p + (V * v) + (U * u);
	}

	//accessor methods
	glm::vec3 getPosition(){return p;}
	glm::vec3 getUp(){return up;}
	glm::vec3 getDirection(){return d;}
	glm::vec3 getLookAt(){return p + focus*d;}
	float getFOV(){return fov;}

	//output the camera settings
	void print(std::ostream& output)
	{
		output<<"Position: "<<glm::to_string(p)<<std::endl;

	}
	friend std::ostream& operator<<(std::ostream& out, const camera& c)
	{
		out<<"Position: "<<glm::to_string(c.p)<<std::endl;
		out<<"Direction: "<<glm::to_string(c.d)<<std::endl;
		out<<"Up: "<<glm::to_string(c.up)<<std::endl;
		out<<"Focal Distance: "<<c.focus<<std::endl;
		return out;
	}

	//constructor
	camera()
	{
		p = glm::vec3(0, 0, 0);
		d = glm::vec3(0, 0, 1);
		up = glm::vec3(0, 1, 0);
		focus = 1;
		fov = 60;

	}

	/// Outputs the camera information as a string
	std::string str(){
		std::stringstream ss;
		ss<<glm::to_string(p)<<"----->"<<glm::to_string(p + d * focus);
		return ss.str();
	}
};

}
