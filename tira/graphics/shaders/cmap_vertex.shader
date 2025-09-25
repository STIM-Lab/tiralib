 R"(
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
)"