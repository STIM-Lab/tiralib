 R"(
# shader vertex

#version 330 core

#define NUM_PTS 11

vec4 brewer[NUM_PTS] = vec4[](
    vec4(0.192157f, 0.211765f, 0.584314f, 1.0f),
	vec4(0.270588f, 0.458824f, 0.705882f, 1.0f),
	vec4(0.454902f, 0.678431f, 0.819608f, 1.0f),
	vec4(0.670588f, 0.85098f, 0.913725f, 1.0f),
	vec4(0.878431f, 0.952941f, 0.972549f, 1.0f),
	vec4(1.0f, 1.0f, 0.74902f, 1.0f),
	vec4(0.996078f, 0.878431f, 0.564706f, 1.0f),
	vec4(0.992157f, 0.682353f, 0.380392f, 1.0f),
	vec4(0.956863f, 0.427451f, 0.262745f, 1.0f),
	vec4(0.843137f, 0.188235f, 0.152941f, 1.0f),
	vec4(0.647059f, 0.0f, 0.14902f, 1.0f)
);

layout(location = 0) in vec4 position;
layout(location = 1) in float value;

uniform mat4 V;							// view matrix
uniform mat4 P;							// projection matrix

uniform float vmin = 1.0f;				// minimum value for the color map
uniform float vmax = 7.0f;				// maximum value for the color map

out vec4 v_color;

void main() {
	gl_Position = P * V * position;			// apply the view and projection transforms to the vertex position

	float v = (value - vmin) / (vmax - vmin);		// normalize the value within the specified range
	v = clamp(v, 0.0f, 1.0f);				// clamp the result to [0, 1] to avoid overruns

	float c = v * (NUM_PTS - 1);
	int c_floor = int(floor(c));
	float m = c - c_floor;
	v_color = mix(brewer[c_floor], brewer[c_floor + 1], m);
	v_color.a = 1.0f;
}

# shader fragment

#version 330 core

layout(location = 0) out vec4 color;

in vec4 v_color;

void main() {
	color = v_color;
}
)"