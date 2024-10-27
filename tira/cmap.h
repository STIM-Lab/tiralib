#pragma once

enum class ColorMap { Grey, Brewer, BrewerBlk, Magma, RainbowCycle };

static float  GREY_PTS[11 * 4] = { 0 * 1.0f / (11 - 1),	0 * 1.0f / (11 - 1),	0 * 1.0f / (11 - 1), 1.0f,
													1 * 1.0f / (11 - 1),	1 * 1.0f / (11 - 1),	1 * 1.0f / (11 - 1), 1.0f,
													1 * 1.0f / (11 - 1),	2 * 1.0f / (11 - 1),	2 * 1.0f / (11 - 1), 1.0f,
													3 * 1.0f / (11 - 1),	3 * 1.0f / (11 - 1),	3 * 1.0f / (11 - 1), 1.0f,
													4 * 1.0f / (11 - 1),	4 * 1.0f / (11 - 1),	4 * 1.0f / (11 - 1), 1.0f,
													5 * 1.0f / (11 - 1),	5 * 1.0f / (11 - 1),	5 * 1.0f / (11 - 1), 1.0f,
													6 * 1.0f / (11 - 1),	6 * 1.0f / (11 - 1),	6 * 1.0f / (11 - 1), 1.0f,
													7 * 1.0f / (11 - 1),	7 * 1.0f / (11 - 1),	7 * 1.0f / (11 - 1), 1.0f,
													8 * 1.0f / (11 - 1),	8 * 1.0f / (11 - 1),	8 * 1.0f / (11 - 1), 1.0f,
													9 * 1.0f / (11 - 1),	9 * 1.0f / (11 - 1),	9 * 1.0f / (11 - 1), 1.0f,
													10 * 1.0f / (11 - 1),	10 * 1.0f / (11 - 1),	10 * 1.0f / (11 - 1), 1.0f };

static float  BREWER_PTS[11 * 4] = { 0.192157f, 0.211765f, 0.584314f, 1.0f,
													0.270588f, 0.458824f, 0.705882f, 1.0f,
													0.454902f, 0.678431f, 0.819608f, 1.0f,
													0.670588f, 0.85098f, 0.913725f, 1.0f,
													0.878431f, 0.952941f, 0.972549f, 1.0f,
													1.0f, 1.0f, 0.74902f, 1.0f,
													0.996078f, 0.878431f, 0.564706f, 1.0f,
													0.992157f, 0.682353f, 0.380392f, 1.0f,
													0.956863f, 0.427451f, 0.262745f, 1.0f,
													0.843137f, 0.188235f, 0.152941f, 1.0f,
													0.647059f, 0.0f, 0.14902f, 1.0f };

static float  BREWERBLK_PTS[11 * 4] = { 0.298, 0.443, 1.000, 1.0f,
													0.176, 0.322, 0.878, 1.0f,
													0.059, 0.204, 0.761, 1.0f,
													0.000, 0.118, 0.561, 1.0f,
													0.000, 0.059, 0.278, 1.0f,
													0.000, 0.000, 0.000, 1.0f,
													0.314, 0.004, 0.020, 1.0f,
													0.624, 0.008, 0.039, 1.0f,
													0.824, 0.067, 0.106, 1.0f,
													0.906, 0.180, 0.216, 1.0f,
													0.988, 0.290, 0.325, 1.0f };


static float  MAGMA_PTS[11 * 4] = { 0.001462f, 0.000466f, 0.013866f, 1.0f,
													0.078815f, 0.054184f, 0.211667f, 1.0f,
													0.232077f, 0.059889f, 0.437695f, 1.0f,
													0.390384f, 0.100379f, 0.501864f, 1.0f,
													0.550287f, 0.161158f, 0.505719f, 1.0f,
													0.716387f, 0.214982f, 0.47529f, 1.0f,
													0.868793f, 0.287728f, 0.409303f, 1.0f,
													0.967671f, 0.439703f, 0.35981f, 1.0f,
													0.994738f, 0.62435f, 0.427397f, 1.0f,
													0.99568f, 0.812706f, 0.572645f, 1.0f,
													0.987053f, 0.991438f, 0.749504f, 1.0f };

static float  RAINBOWCYCLE_PTS[7 * 4] = { 1.0f, 0.0f, 0.0f, 1.0f,
													0.0f, 1.0f, 0.0f, 1.0f,
													0.0f, 0.0f, 1.0f, 1.0f,
													1.0f, 0.0f, 0.0f, 1.0f,
													0.0f, 1.0f, 0.0f, 1.0f,
													0.0f, 0.0f, 1.0f, 1.0f,
													1.0f, 0.0f, 0.0f, 1.0f };

namespace tira::cmap {
	template<typename T>
	void cmap(T val, T vmin, T vmax, unsigned char& r, unsigned char& g, unsigned char& b, ColorMap colormap) {
		T a = 0.5;
		if (vmin != vmax)
			a = (val - vmin) / (vmax - vmin);

		a = (std::max)(0.0f, a);																// clamp the normalized value to between zero and 1
		a = (std::min)(1.0f, a);



		float* ctrlPts;
		unsigned int numPts;
		switch (colormap) {
		case ColorMap::Grey:
			ctrlPts = GREY_PTS;
			numPts = 11;
			break;
		case ColorMap::Brewer:
			ctrlPts = BREWER_PTS;
			numPts = 11;
			break;
		case ColorMap::BrewerBlk:
			ctrlPts = BREWERBLK_PTS;
			numPts = 11;
			break;
		case ColorMap::Magma:
			ctrlPts = MAGMA_PTS;
			numPts = 11;
			break;
		case ColorMap::RainbowCycle:
			ctrlPts = RAINBOWCYCLE_PTS;
			numPts = 7;
			break;
		default:
			throw std::runtime_error("Invalid Colormap");
		}
		const float c = a * static_cast<float>(numPts - 1);								// get the real value to interpolate control points
		const int c_floor = static_cast<int>(c);										// calculate the control point index
		const float m = c - static_cast<float>(c_floor);								// calculate the fractional part of the control point index

		float fr = ctrlPts[c_floor * 4 + 0];											// use a LUT to find the "low" color value
		float fg = ctrlPts[c_floor * 4 + 1];
		float fb = ctrlPts[c_floor * 4 + 2];
		if (c_floor != numPts - 1) {										// if there is a fractional component, interpolate
			fr = fr * (1.0f - m) + ctrlPts[(c_floor + 1) * 4 + 0] * m;
			fg = fg * (1.0f - m) + ctrlPts[(c_floor + 1) * 4 + 1] * m;
			fb = fb * (1.0f - m) + ctrlPts[(c_floor + 1) * 4 + 2] * m;
			if (fb > 1)
				std::cout << "ERROR" << std::endl;
		}
		r = 255 * fr;											// save the resulting color in the output image
		g = 255 * fg;
		b = 255 * fb;
	}

}