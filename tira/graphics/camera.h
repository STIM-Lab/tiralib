#pragma once

#define GLM_FORCE_LEFT_HANDED
#define GLM_ENABLE_EXPERIMENTAL

#include <tira/cuda/callable.h>

#include "glm/glm.hpp"
#include <glm/gtx/string_cast.hpp>

#include <glm/gtc/quaternion.hpp>

#include <ostream>
#include <sstream>
#include <numbers>

namespace tira{

class camera
{
	glm::vec3 _view;					// view vector (direction the camera is looking)
	glm::vec3 _position;				// camera location (eye point)
	glm::vec3 _up;						// vector pointing upwards (from the camera location)
	float _focus;						// focal length of the camera
	float _fov;							// camera field of view
	float _imagedist;					// distance to a unit-sized image plane

	//private function that updates the up vector to ensure that it is orthogonal to the view vector
	CUDA_CALLABLE void _stabalize() {
		const glm::vec3 side = glm::cross(_up, _view);
		_up = glm::cross(_view, side);
		_up = glm::normalize(_up);
		_view = glm::normalize(_view);
	}

	//private function that updates the VIEW vector to ensure that it is orthogonal to the UP vector
	CUDA_CALLABLE void _stabalize_view() {
		const glm::vec3 side = glm::cross(_up, _view);
		_up = glm::cross(_view, side);
		_up = glm::normalize(_up);
		_view = glm::normalize(_view);
	}

	// private function that calculates the distance to an image plane of size 1 given an FOV (in degrees)
	CUDA_CALLABLE void _update_imagedist() {
		const float fov_radians = _fov * static_cast<float>(std::numbers::pi) / 180.0f;     // convert the field of view to radians
		const float tan_theta = std::tan(fov_radians / 2.0f);								// calculate the tangent of the fov half-angle
		_imagedist = 0.5f * (1.0f / tan_theta);                                             // return the distance to the plane
	}

public:
	CUDA_CALLABLE void position(const glm::vec3 pos)	{
		_position = pos;
	}
	CUDA_CALLABLE void position(const float x, const float y, const float z){ position(glm::vec3(x, y, z));}
	CUDA_CALLABLE void up(const glm::vec3 up) {
		_up = up;
		_stabalize_view();
	}
	CUDA_CALLABLE void focal_length(const float distance){_focus = distance;}
	CUDA_CALLABLE void fov(const float field_of_view) {
		_fov = field_of_view;
		_update_imagedist();
	}

	CUDA_CALLABLE void lookat(const glm::vec3 pos) {
		_view = pos - _position;								// calculate the new direction		
		_focus = glm::length(_view);							// calculate the distance from the look-at point to the current position		
		_stabalize();											// stabalize the camera
	}
	CUDA_CALLABLE void lookat(const float px, const float py, const float pz) { lookat(glm::vec3(px, py, pz)); }
	CUDA_CALLABLE void lookat(const glm::vec3 pos, const glm::vec3 new_up) {_up = new_up; lookat(pos);}
	CUDA_CALLABLE void lookat(const float px, const float py, const float pz,
		const float ux, const float uy, const float uz) {
		lookat(glm::vec3(px, py, pz), glm::vec3(ux, uy, uz));
	}
	CUDA_CALLABLE glm::vec3 ray(const float x, const float y) const {
		const glm::vec3 side = glm::cross(_view, _up);
		const glm::vec3 result = x * side + y * _up + _imagedist * _view;
		return glm::normalize(result);
	}

	// translates the camera position based on the specified vector
	CUDA_CALLABLE void dolly(const glm::vec3 direction)	{
		_position = _position + direction;
	}
	CUDA_CALLABLE void dolly(const float x, const float y, const float z) { dolly(glm::vec3(x, y, z)); }

	/// <summary>
	/// Move the camera a specified distance from the view position
	/// </summary>
	/// <param name="d"></param>
	CUDA_CALLABLE void distance(const float d) {
		const glm::vec3 focal_point = _position + _focus * _view;

		_position = focal_point - _view * d;
		_focus = d;
	}

	// Push the camera along the view direction
	CUDA_CALLABLE void push(float delta) {
		if(delta > _focus)
			delta = _focus;

		_focus -= delta;

		dolly(_view * delta);
	}

	CUDA_CALLABLE void pan(const float theta_x, const float theta_y, const float theta_z) {
		
		const glm::quat qx = glm::angleAxis(theta_x, _up);				//x rotation is around the up axis
		const glm::vec3 side = glm::cross(_up, _view);				//y rotation is around the side axis
		const glm::quat qy = glm::angleAxis(theta_y, side);
		const glm::quat qz = glm::angleAxis(theta_z, _view);				//z rotation is around the direction vector
		const glm::quat final = qz * qy * qx;								//combine the rotations in x, y, z order
		const glm::mat3 rotation = glm::mat3_cast(final);					//get the rotation matrix
		_view = rotation * _view;									//apply the rotation
		_up = rotation * _up;

		_stabalize();											//stabalize the camera

	}
	CUDA_CALLABLE void pan(const float theta_x){ pan(theta_x, 0, 0); }
	CUDA_CALLABLE void tilt(const float theta_y){ pan(0, theta_y, 0); }
	CUDA_CALLABLE void twist(const float theta_z){ pan(0, 0, theta_z); }

	// adjusts the camera fov
	CUDA_CALLABLE void zoom(const float delta) {
		float new_fov = _fov - delta;
		if(_fov < 0.5)
			new_fov = 0.5;
		if(_fov > 180)
			new_fov = 180;
		fov(new_fov);
	}

	// orbits the camera around the focal point by the specified angles (in radians)
	CUDA_CALLABLE void orbit(const float theta_x, const float theta_y) {
		//find the focal point
		const glm::vec3 focal_point = _position + _focus * _view;

		//center the coordinate system on the focal point
		glm::vec3 centered = _position - (focal_point - glm::vec3(0, 0, 0));

		//create the x rotation (around the up vector)
		const glm::quat qx = glm::angleAxis(theta_x, _up);
		centered = glm::vec3(0, 0, 0) + glm::mat3_cast(qx)*(centered - glm::vec3(0, 0, 0));

		//get a side vector for theta_y rotation
		const glm::vec3 side = glm::normalize(glm::cross(_up, glm::vec3(0, 0, 0) - centered));
		const glm::quat qy = glm::angleAxis(theta_y, side);
		centered = glm::vec3(0, 0, 0) + glm::mat3_cast(qy)*(centered - glm::vec3(0, 0, 0));

		//re-position the camera
		_position = centered + (focal_point - glm::vec3(0, 0, 0));

		//make sure we are looking at the focal point
		lookat(focal_point);

		//stabalize the camera
		_stabalize();

	}

	CUDA_CALLABLE void slide(const float u, const float v) {
		const glm::vec3 V = glm::normalize(_up);
		const glm::vec3 U = glm::normalize(glm::cross(_up, _view));

		_position = _position + (V * v) + (U * u);
	}

	//accessor methods
	CUDA_CALLABLE glm::vec3 position() const { return _position; }
	CUDA_CALLABLE glm::vec3 up() const { return _up; }
	CUDA_CALLABLE glm::vec3 view() const { return _view; }
	CUDA_CALLABLE glm::vec3 lookat() const { return _position + _focus * _view; }
	CUDA_CALLABLE float fov() const { return _fov; }

	// returns a vector pointing along the positive x-axis from the viewpoint of the camera
	CUDA_CALLABLE glm::vec3 side() const {
		return glm::cross(_up, _view);
	}

	CUDA_CALLABLE glm::mat4 viewmatrix() const {
		return glm::lookAt(_position, _position + _view, _up);
	}

	/// <summary>
	/// Generates an orthographic projection matrix such that the image plane (based on the FOV and focal distance) is entirely
	/// withing the view frustum. The depth of field is twice the focal distance
	/// </summary>
	/// <returns></returns>
	CUDA_CALLABLE glm::mat4 orthomatrix(float aspect = 1.0f, float zoom_factor = 1.0f, float UP = 0.0f, float RIGHT = 0.0f) const {
		const float fov_rad = glm::radians(_fov);
		const float view_width = _focus * std::tan(fov_rad / 2.0f);
		const float dof = 2 * _focus;
		return glm::ortho(-view_width * aspect / zoom_factor + RIGHT, view_width * aspect / zoom_factor + RIGHT, 
			-view_width / zoom_factor + UP, view_width / zoom_factor + UP, 0.0f, dof);
	}

	CUDA_CALLABLE glm::mat4 perspectivematrix(const float aspect = 1.0f, const float depthfactor = 2.0f) const {
		const float fov_rad = glm::radians(_fov);
		return glm::perspective(fov_rad, aspect, _focus/ depthfactor, depthfactor * _focus);
	}

	//output the camera settings
	void print(std::ostream& output) const {
		output<<"Position: "<<glm::to_string(_position)<<std::endl;
	}
	friend std::ostream& operator<<(std::ostream& out, const camera& c)	{
		out<<"Position: "<<glm::to_string(c._position)<<std::endl;
		out<<"Direction: "<<glm::to_string(c._view)<<std::endl;
		out<<"Up: "<<glm::to_string(c._up)<<std::endl;
		out<<"Focal Distance: "<<c._focus<<std::endl;
		return out;
	}

	//constructor
	CUDA_CALLABLE camera() {
		position(glm::vec3(0, 0, 0));
		lookat(glm::vec3(0, 0, 1));
		up(glm::vec3(0, 1, 0));
		fov(60);

	}

	/// Outputs the camera information as a string
	std::string str() const {
		std::stringstream ss;
		ss<<glm::to_string(_position)<<"----->"<<glm::to_string(_position + _view * _focus);
		return ss.str();
	}
};

}
