#pragma once

#include "glm/glm.hpp"

class volume{

protected:
	void* m_data;
	size_t m_num_pixels[3];
	glm::vec3 m_position;
	glm::vec3 m_scale;
//	glm::quaternion m_orientation;

public:
	glm::mat4 getAlignmentMatrix(){
		glm::mat4 result;
		// generate an alignment matrix and return it
		return result;
	}
};