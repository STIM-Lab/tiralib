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

uniform vec4 c = vec4(1.0, 0.0, 1.0, 1.0);
layout(location = 0) out vec4 color;

void main() {
	color = c;
}
)"