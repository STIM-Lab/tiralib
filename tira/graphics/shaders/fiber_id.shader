R"(

# shader vertex

#version 330 core

layout(location = 0) in vec4 position;
uniform mat4 V;							// view matrix
uniform mat4 P;							// projection matrix


void main() {
	gl_Position = P * V * position;
}

# shader fragment

#version 330 core

layout(location = 0) out int id_val;
uniform int id;

void main() {
	id_val = id;
}
)" 
