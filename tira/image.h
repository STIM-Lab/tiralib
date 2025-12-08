#pragma once

#define cimg_use_png
#define cimg_use_jpeg
#include "extern/CImg.h"
#undef min

#include <stack>
#include <algorithm>

#include <tira/field.h>
#include <tira/cmap.h>



namespace tira {
	/// This static class provides the STIM interface for loading, saving, and storing 2D images.
	/// Data is stored in an interleaved (BIP) format (default for saving and loading is RGB).

	template <class T>
	class image : public field<T>{		// the image class extends field

	protected:

		void init(std::vector<size_t> image_shape) {
			if (field<T>::m_shape.size() != 0) {													// if this image has already been allocated
				field<T>::m_shape.clear();														// clear the shape and data vectors to start from scratch
				field<T>::m_data.clear();
			}
			field<T>::m_shape.push_back(image_shape[1]);
			field<T>::m_shape.push_back(image_shape[0]);
			if (image_shape.size() == 3)
				field<T>::m_shape.push_back(image_shape[2]);

			field<T>::m_Allocate();
		}

		/// <summary>
		/// Allocate space for an empty image
		/// </summary>
		/// <param name="x">Width of the image (fast axis)</param>
		/// <param name="y">Height of the image (slow axis)</param>
		/// <param name="c">Number of color channels</param>
		void init(size_t x, size_t y, size_t c = 0) {
			std::vector<size_t> image_shape;
			image_shape.push_back(x);
			image_shape.push_back(y);
			if (c > 1) image_shape.push_back(c);
			init(image_shape);
		}

		/// <summary>
		/// Calculates the 1D array offset given spatial and color coordinates
		/// </summary>
		/// <param name="x">X coordinate of the image (fast axis)</param>
		/// <param name="y">Y coordinate of the image (slow axis)</param>
		/// <param name="c">Color channel</param>
		/// <returns></returns>
		inline size_t idx_offset(size_t x, size_t y, size_t c = 0) const {
			size_t idx =  y * C() * X() + x * C() + c;		// y * C * X + x * C + c
			return idx;
		}

		T& at(size_t x, size_t y, size_t c = 0) {
			size_t idx = idx_offset(x, y, c);
			return field<T>::m_data[idx];
		}

		/// <summary>
		/// Calculate a distance field from an input contour phi = 0
		/// </summary>
		/// <returns></returns>
		tira::image<T> dist(tira::image<int>& binary_boundary) {

			const int w = width();
			const int h = height();



			std::vector<std::tuple<int, int>> neighbors;				// vector stores a template for 4-connected indices
			neighbors.emplace_back(0, 1);
			neighbors.emplace_back(1, 0);
			neighbors.emplace_back(-1, 0);
			neighbors.emplace_back(0, -1);


			// indentifying boundary cells
			for (int y = 1; y < h - 1; y++) {						// for every row in the image
				for (int x = 1; x < w - 1; x++) {					// for every column in the image
					for (int k = 0; k < neighbors.size(); k++) {		// for every neighbor

						int nx = x + get<0>(neighbors[k]);				// calculate the x coordinate of the neighbor cell
						int ny = y + get<1>(neighbors[k]);				// calculate the y coordinate of the neighbor cell

						if (at(x, y) * at(nx, ny) <= 0) {				// if the product of the current cell and neighboring cell is negative (it is a boundary cell)
							binary_boundary(x, y) = 1;					// this cell is a boundary cell
							binary_boundary(nx, ny) = 1;				// the neighboring cell is ALSO a boundary cell
						}
					}
				}
			}



			tira::image<float> dist(w, h);										// create an image to store the distance field
			dist = 9999.0f;																// initialize the distance field to a very large value


			// calculate the distance for all boundary cells to the contour
			for (int y = 1; y < h - 1; y++) {											// for every row in the image
				for (int x = 1; x < w - 1; x++) {										// for every column in the image
					if (binary_boundary(x, y)) {											// if the pixel (x, y) is in the boundary
						for (int k = 0; k < neighbors.size(); k++) {						// for every neighbor

							int nx = x + get<0>(neighbors[k]);								// calculate the x coordinate of the neighbor cell
							int ny = y + get<1>(neighbors[k]);								// calculate the y coordinate of the neighbor cell
							if (binary_boundary(nx, ny)) {												// if the neighboring cell (nx, ny) is ALSO in the boundary
								float da = (abs(at(x, y))) / (abs(at(nx, ny) - at(x, y)));				// calculate distance from pixel(x,y) to contour da
								float db = (abs(at(nx, ny))) / (abs(at(nx, ny) - at(x, y)));			// calculate distance from neighbor to contour db
								float dist_xy = dist(x, y);
								float dist_nxny = dist(nx, ny);
								dist(x, y) = std::min(dist_xy, da);									// minimum between distance and large boundary value of pixel (x,y)
								dist(nx, ny) = std::min(dist_nxny, db);								// minimum between distance and large boundary value of neighbor (nx, ny)
							}
						}
					}
				}
			}
			// at this point, the distance image has the correct distances in cells surround the boundary

			std::vector<float> dist1d;											// create a 1d representation of the distance field
			dist1d.resize(h * w);

			int row = w;

			// copy the distance field into the 1d distance field
			for (int y = 0; y < h; y++) {
				for (int x = 0; x < w; x++) {
					dist1d[y * w + x] = dist(x, y);
				}
			}

			// initializing fast sweeiping algorithm 
			constexpr int NSweeps = 4;

			//// sweep directions { start, end, step }
			const int dirX[NSweeps][3] = { {0, w - 1, 1} , {w - 1, 0, -1}, {w - 1, 0, -1}, {0, w - 1, 1} };
			const int dirY[NSweeps][3] = { {0, h - 1, 1}, {0, h - 1, 1}, {h - 1, 0, -1}, {h - 1, 0, -1} };
			double aa[2];

			for (int s = 0; s < NSweeps; s++) {

				for (int iy = dirY[s][0]; dirY[s][2] * iy <= dirY[s][1]; iy += dirY[s][2]) {
					for (int ix = dirX[s][0]; dirX[s][2] * ix <= dirX[s][1]; ix += dirX[s][2]) {
						constexpr double f = 1.0;
						constexpr double dx = 1.0;

						const int gridPos = iy * row + ix;

						if (iy == 0 || iy == (h - 1)) {                    // calculation for ymin
							if (iy == 0) {
								aa[1] = dist1d[(iy + 1) * row + ix];
							}
							if (iy == (h - 1)) {
								aa[1] = dist1d[(iy - 1) * row + ix];
							}
						}
						else {
							aa[1] = dist1d[(iy - 1) * row + ix] < dist1d[(iy + 1) * row + ix] ? dist1d[(iy - 1) * row + ix] : dist1d[(iy + 1) * row + ix];
						}

						if (ix == 0 || ix == (w - 1)) {                    // calculation for xmin
							if (ix == 0) {
								aa[0] = dist1d[iy * row + (ix + 1)];
							}
							if (ix == (w - 1)) {
								aa[0] = dist1d[iy * row + (ix - 1)];
							}
						}
						else {
							aa[0] = dist1d[iy * row + (ix - 1)] < dist1d[iy * row + (ix + 1)] ? dist1d[iy * row + (ix - 1)] : dist1d[iy * row + (ix + 1)];
						}

						const double a = aa[0];
						const double b = aa[1];
						const double d_new = (fabs(a - b) < f * dx ? (a + b + sqrt(2.0 * f * f * dx * dx - (a - b) * (a - b))) * 0.5 : std::fminf(a, b) + f * dx);

						dist1d[gridPos] = d_new;
					}
				}
			}

			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					dist(x, y) = dist1d[y * w + x];
				}

			}

			return dist;
		}
		
		
	public:
		inline size_t Y() const { return field<T>::m_shape[0]; }
		inline size_t X() const { return field<T>::m_shape[1]; }
		inline size_t C() const {
			if(field<T>::m_shape.size() <=2) return 1;
			return field<T>::m_shape[2];
		}

		/// <summary>
		/// Default constructor initializes an empty image (0x0 with 0 channels)
		/// </summary>
		image() : field<T>() {}			//initialize all variables, don't allocate any memory

		image(std::vector<size_t> shape) {
			field<T>::m_shape = shape;
			field<T>::m_Allocate();
		}

		/// <summary>
		/// Constructor creates an image from a field
		/// </summary>
		/// <param name="F"></param>
		image(field<T> F) : image(F.Shape()) {
			//field<T>::m_shape = F.shape();
			//field<T>::m_Allocate();
			memcpy(&field<T>::m_data[0], F.Data(), F.Bytes());
		}

		

		/// <summary>
		/// Create a new image from scratch given a number of samples and channels
		/// </summary>
		/// <param name="x">size of the image along the X (fast) axis</param>
		/// <param name="y">size of the image along the Y (slow) axis</param>
		/// <param name="c">number of channels</param>
		image(size_t x, size_t y = 1, size_t c = 0) {
			init(x, y, c);
		}

		/// <summary>
		/// Create a new image with the data given in 'data'
		/// </summary>
		/// <param name="data">pointer to the raw image data</param>
		/// <param name="x">size of the image along the X (fast) axis</param>
		/// <param name="y">size of the image along the Y (slow) axis</param>
		/// <param name="c">number of channels (channels are interleaved)</param>
		image(T* data, size_t x, size_t y, size_t c = 0) : image(x, y, c) {
			memcpy(&field<T>::m_data[0], data, field<T>::Bytes());									// use memcpy to copy the raw data into the field array
		}

		/// <summary>
		/// Copy constructor - duplicates an image object
		/// </summary>
		/// <param name="I">image object to be copied</param>
		image(const image<T>& I) : image(I.X(), I.Y(), I.C()) {
			memcpy(&field<T>::m_data[0], &I.m_data[0], field<T>::Bytes());
		}

		/// <summary>
		/// Creates a new image from an image file
		/// </summary>
		/// <param name="filename">Name of the file to load (uses CImg)</param>
		image(std::string filename) {
			load(filename);
		}

		/// Destructor
		~image() {
			
		}

		/// <summary>
		/// Return the image width
		/// </summary>
		/// <returns>Width of the image (fast dimension)</returns>
		size_t width() const {
			return X();
		}

		/// <summary>
		/// Return the height of the image
		/// </summary>
		/// <returns>Height of the image (slow dimension)</returns>
		size_t height() const {
			return Y();
		}

		/// <summary>
		/// Returns the number of channels in the iamge
		/// </summary>
		/// <returns>Number of channels</returns>
		size_t channels() const {
			return C();
		}

		/// <summary>
		/// Assignment operator, assigns one image to another
		/// </summary>
		/// <param name="I"></param>
		/// <returns>pointer to the current object</returns>
		image<T>& operator=(const image<T>& I) {
			if (&I == this)													// handle self-assignment
				return *this;
			init(I.X(), I.Y(), I.C());
			memcpy(&field<T>::m_data[0], &I.m_data[0], field<T>::Bytes());		// copy the data
			return *this;													// return a pointer to the current object
		}

		
		image<T> operator+(T rhs) {
			size_t N = field<T>::Size();							//calculate the total number of values in the image
			image<T> r(this->Shape());								//allocate space for the resulting image
			for (size_t n = 0; n < N; n++)
				r.m_data[n] = field<T>::m_data[n] + rhs;				//add the individual pixels
			return r;												//return the summed result
		}

		image<T> operator-(const image<T>& rhs) const {
			size_t N = field<T>::Size();                           //allocate space for the resulting image
			image<T> r(this->Shape());

			for (size_t n = 0; n < N; n++)                        
				r.m_data[n] = field<T>::m_data[n] - rhs.m_data[n];

			return r;                                              // return the result
		}

		image<T> operator*(T rhs) {
			size_t N = field<T>::Size();							//calculate the total number of values in the image
			image<T> r(this->Shape());								//allocate space for the resulting image
			for (size_t n = 0; n < N; n++)
				r.m_data[n] = field<T>::m_data[n] * rhs;				//add the individual pixels
			return r;												//return the summed result
		}

		image<T> operator/(T rhs) {
			size_t N = field<T>::Size();							//calculate the total number of values in the image
			image<T> r(this->Shape());								//allocate space for the resulting image
			for (size_t n = 0; n < N; n++)
				r.m_data[n] = field<T>::m_data[n] / rhs;				//add the individual pixels
			return r;												//return the summed result
		}

		image<T> operator*(image<T> rhs) {							// point-wise multiplication
			if(X()!=rhs.X() || Y()!=rhs.Y())
				throw std::runtime_error("Images dimensions are incompatible");

			if(C() == rhs.C()) {					// if both images have the same number of color channels
				tira::image<T> result(this->Shape());
				for(size_t i = 0; i < field<T>::Size(); i++) {
					result.m_data[i] = field<T>::m_data[i] * rhs.m_data[i];
				}
				return result;
			}
			else if(C() == 1) {
				tira::image<T> result(X(), Y(), rhs.C());
				for(size_t yi = 0; yi < Y(); yi++) {
					for(size_t xi = 0; xi < X(); xi++) {
						for(size_t ci = 0; ci < rhs.C(); ci++) {
							result(xi, yi, ci) = at(xi, yi) * rhs(xi, yi, ci);
						}
					}
				}
				return result;
			}
			else if(rhs.C() == 1) {
				tira::image<T> result(this->Shape());
				for(size_t yi = 0; yi < Y(); yi++) {
					for(size_t xi = 0; xi < X(); xi++) {
						for(size_t ci = 0; ci < C(); ci++) {
							result(xi, yi, ci) = at(xi, yi, ci) * rhs(xi, yi);
						}
					}
				}
				return result;
			}
			else {
				throw std::runtime_error("Number of color channels are incompatible");
			}
		}

		image<T> operator/(image<T> rhs) {							// point-wise multiplication
			if (X() != rhs.X() || Y() != rhs.Y())
				throw std::runtime_error("Images dimensions are incompatible");

			if (C() == rhs.C()) {					// if both images have the same number of color channels
				tira::image<T> result(this->Shape());
				for (size_t i = 0; i < field<T>::Size(); i++) {
					result.m_data[i] = field<T>::m_data[i] / rhs.m_data[i];
				}
				return result;
			}
			else if (C() == 1) {
				tira::image<T> result(X(), Y(), rhs.C());
				for (size_t yi = 0; yi < Y(); yi++) {
					for (size_t xi = 0; xi < X(); xi++) {
						for (size_t ci = 0; ci < rhs.C(); ci++) {
							result(xi, yi, ci) = at(xi, yi) / rhs(xi, yi, ci);
						}
					}
				}
				return result;
			}
			else if (rhs.C() == 1) {
				tira::image<T> result(X(), Y(), C());
				for (size_t yi = 0; yi < Y(); yi++) {
					for (size_t xi = 0; xi < X(); xi++) {
						for (size_t ci = 0; ci < C(); ci++) {
							result(xi, yi, ci) = at(xi, yi, ci) / rhs(xi, yi);
						}
					}
				}
				return result;
			}
			else {
				throw std::runtime_error("Number of color channels are incompatible");
			}
		}

		image<T> operator-(T rhs) {
			size_t N = field<T>::size();							//calculate the total number of values in the image
			image<T> r(this->shape());								//allocate space for the resulting image
			for (size_t n = 0; n < N; n++)
				r.m_data[n] = field<T>::m_data[n] - rhs;				//add the individual pixels
			return r;												//return the summed result
		}

		image<T> operator-() {
			return static_cast< image<T> >(field<T>::operator-());
		}

		/// <summary>
		/// Set all values in the image to a single constant
		/// </summary>
		/// <param name="v">Constant that all elements will be set to</param>
		image<T>& operator=(T v) {														//set all elements of the image to a given value v
			const size_t N = field<T>::Size();
			for(size_t n = 0; n < N; n++)
				field<T>::m_data[n] = v;
			return *this;
		}

		/// <summary>
		/// Adds two images
		/// </summary>
		/// <param name="rhs"></param>
		/// <returns></returns>
		image<T> operator+(image<T> rhs) {
			size_t N = field<T>::Size();							//calculate the total number of values in the image
			image<T> r(this->Shape());								//allocate space for the resulting image
			for (size_t n = 0; n < N; n++)
				r.m_data[n] = field<T>::m_data[n] + rhs.m_data[n];		//add the individual pixels
			return r;												//return the inverted image
		}

		image<T> abs() {
			image<T> result = field<T>::Abs();
			return result;
		}

		T& operator()(size_t x, size_t y, size_t c = 0) {
			return at(x, y, c);
		}

		image<T> clamp(T min, T max) {
			size_t N = field<T>::Size();
			image<T> r(this->Shape());
			for (size_t n = 0; n < N; n++) {
				if (field<T>::m_data[n] < min) r.m_data[n] = min;
				else if (field<T>::m_data[n] > max) r.m_data[n] = max;
				else r.m_data[n] = field<T>::m_data[n];
			}

			return r;
		}

		/// <summary>
		/// Cast data types
		/// </summary>
		/// <typeparam name="V"></typeparam>
		template<typename V>
		operator image<V>() {
			image<V> r(this->Shape());					//create a new image
			std::copy(&field<T>::m_data[0], &field<T>::m_data[0] + field<T>::Size(), r.Data());		//copy and cast the data
			return r;									//return the new image
		}

		/// <summary>
		/// Save the image as a file (uses the CImg library)
		/// </summary>
		/// <param name="fname">name of the file</param>
		void save(const std::string fname) {
			cimg_library::CImg<T> cimg(static_cast<unsigned int>(X()), static_cast<unsigned int>(Y()), 1, static_cast<unsigned int>(C()));
			get_noninterleaved(cimg.data());
			cimg.save(fname.c_str());
		}

		/// <summary>
		/// Load an image file (uses the CImg library)
		/// </summary>
		/// <param name="filename">Name of the file to load</param>
		void load(const std::string filename) {
			cimg_library::CImg<T> cimg(filename.c_str());									//create a CImg object for the image file
			set_noninterleaved(cimg.data(), cimg.width(), cimg.height(), cimg.spectrum());
		}

		/// <summary>
		/// Returns a channel of the current image as an independent one-channel image
		/// </summary>
		/// <param name="c">channel to return</param>
		/// <returns></returns>
		image<T> channel(const size_t c) const {
			image<T> r(X(), Y());											//create a new single-channel image
			for (size_t x = 0; x < X(); x++) {
				for (size_t y = 0; y < Y(); y++) {
					size_t source_i = idx_offset(x, y, c);
					size_t dest_i = r.idx_offset(x, y, 0);
					r.m_data[dest_i] = field<T>::m_data[source_i];
				}
			}
			return r;
		}

		/// <summary>
		/// Set the specified image channel to the data provided in src
		/// </summary>
		/// <param name="src">raw data used as the channel source</param>
		/// <param name="c">channel that the raw data is assigned to</param>
		void channel(T* src, const size_t c) {
			for (size_t y = 0; y < Y(); y++) {
				for (size_t x = 0; x < X(); x++) {
					field<T>::m_data[idx_offset(x, y, c)] = src[y * X() + x];
				}
			}
		}

		/// <summary>
		/// Set the all pixels in the specified channel to a single value
		/// </summary>
		/// <param name="val">Value that all pixels in the channel will be set to</param>
		/// <param name="c">Channel to be modified</param>
		void channel(T val, size_t c) {
			for (size_t y = 0; y < Y(); y++) {
				for (size_t x = 0; x < X(); x++) {
					field<T>::m_data[idx_offset(x, y, c)] = val;
				}
			}
		}

		/// <summary>
		/// Splits an image into sub-channels and returns each one independently in an std::vector of images
		/// </summary>
		/// <returns></returns>
		std::vector< image<T> > split() const {
			std::vector< image<T> > r;			//create an image array
			r.resize(C());						//create images for each channel

			for (size_t c = 0; c < C(); c++) {	//for each channel
				r[c] = channel(c);				//copy the channel image to the array
			}
			return r;
		}

		/// <summary>
		/// Merge a set of single-channel images into a multi-channel image
		/// </summary>
		/// <param name="list"></param>
		void merge(std::vector< image<T> >& list) {
			size_t x = list[0].width();				//calculate the size of the image
			size_t y = list[0].height();
			init(x, y, list.size());			//re-allocate the image
			for (size_t c = 0; c < list.size(); c++)		//for each channel
				channel(&list[c].channel(0).m_data[0], c);	//insert the channel into the output image
		}

		/// <summary>
		/// Returns the number of foreground pixels (given a specific background value)
		/// </summary>
		/// <param name="background">Value used as background</param>
		/// <returns>Number of foreground pixels</returns>
		size_t foreground(T background = 0) {

			size_t N = X() * Y() * C();

			size_t num_foreground = 0;
			for (size_t n = 0; n < N; n++)
				if (field<T>::m_data[n] != background) num_foreground++;

			return num_foreground;	//return the number of foreground values
		}

		/// <summary>
		/// Returns the number of non-zero values in the image
		/// </summary>
		/// <returns>Number of non-zero values</returns>
		size_t nnz() {
			return foreground(0);
		}

		/// <summary>
		/// Returns the indices of all non-zero values
		/// </summary>
		/// <returns></returns>
		std::vector<size_t> sparse_idx(T background = 0) {

			std::vector<size_t> s;								// allocate an array
			s.resize(foreground(background));					// allocate space in the array

			size_t N = field<T>::size();

			size_t i = 0;
			for (size_t n = 0; n < N; n++) {
				if (field<T>::m_data[n] != background) {
					s[i] = n;
					i++;
				}
			}

			return s;			//return the index list
		}

		/// <summary>
		/// Returns the maximum value in the image
		/// </summary>
		/// <returns></returns>
		T maxv() {
			T max_val = field<T>::m_data[0];				//initialize the maximum value to the first one
			size_t N = field<T>::Size();	//get the number of pixels

			for (size_t n = 0; n < N; n++) {		//for every value

				if (field<T>::m_data[n] > max_val) {			//if the value is higher than the current max
					max_val = field<T>::m_data[n];
				}
			}
			return max_val;
		}


		/// <summary>
		/// Returns the minimum value in the image
		/// </summary>
		/// <returns></returns>
		T minv() {
			T min_val = field<T>::m_data[0];				//initialize the maximum value to the first one
			size_t N = field<T>::Size();	//get the number of pixels

			for (size_t n = 0; n < N; n++) {		//for every value
				if (field<T>::m_data[n] < min_val) {			//if the value is higher than the current max
					min_val = field<T>::m_data[n];
				}
			}

			return min_val;
		}

		/// <summary>
		/// Stretches the contrast of an image such that the highest and lowest values fall within those specified
		/// </summary>
		/// <param name="low">Lowest value allowed in the image</param>
		/// <param name="high">Highest value allowed in the image</param>
		/// <returns></returns>
		image<T> stretch_contrast(T low, T high) {
			T maxval = maxv();
			T minval = minv();
			image<T> result = *this;				//create a new image for output
			if (maxval == minval) {					//if the minimum and maximum values are the same, return an image composed of low
				result = low;
				return result;
			}

			size_t N = field<T>::Size();						//get the number of values in the image
			T range = maxval - minval;			//calculate the current range of the image
			T desired_range = high - low;		//calculate the desired range of the image
			for (size_t n = 0; n < N; n++) {		//for each element in the image
				result.m_data[n] = desired_range * (field<T>::m_data[n] - minval) / range + low;
			}
			return result;
		}

		/// <summary>
		/// Add a border to an image
		/// </summary>
		/// <param name="w">Width of the border in pixels</param>
		/// <param name="value">Value used to generate the border</param>
		/// <returns></returns>
		image<T> border(const size_t w, const T value) {
			image<T> result(width() + w * 2, height() + w * 2, channels());						//create an output image
			result = value;														//assign the border value to all pixels in the new image
			for (size_t y = 0; y < height(); y++) {								//for each pixel in the original image
				for (size_t x = 0; x < width(); x++) {
					size_t n = (y + w) * (width() + w * 2) + x + w;				//calculate the index of the corresponding pixel in the result image
					size_t n0 = idx_offset(x, y);										//calculate the index for this pixel in the original image
					result.m_data[n] = field<T>::m_data[n0];									// copy the original image to the result image afer the border area
				}
			}
			return result;
		}

		image<T> border(std::vector<size_t> width, T val=0) {
			return static_cast< image<T> >(field<T>::border(width, val));
		}

		/// <summary>
		/// Generates a border by replicating edge pixels
		/// </summary>
		/// <param name="w"> width of the border (padding) to create</param>
		/// <returns></returns>
		image<T> border_replicate(const size_t w) {

			image<T> result(width() + w * 2, height() + w * 2, channels());						//create an output image
													//assign the border value to all pixels in the new image
			for (size_t y = 0; y < height(); y++) {								//for each pixel in the original image
				for (size_t x = 0; x < width(); x++) {
					size_t n = (y + w) * (width() + w * 2) + x + w;				//calculate the index of the corresponding pixel in the result image
					size_t n0 = idx_offset(x, y);										//calculate the index for this pixel in the original image
					result.data()[n] = field<T>::m_data[n0];									// copy the original image to the result image afer the border area
				}
			}
			const size_t l = w;
			const size_t r = w + width() - 1;
			const size_t t = w;
			const size_t b = w + height() - 1;
			for (size_t y = 0; y < w; y++) {
				for (size_t x = l; x <= r; x++) {
					result(x, y) = result(x, t);						//pad the top
				}
			}
			for (size_t y = b + 1; y < result.height(); ++y) {
				for (size_t x = l; x <= r; x++) {
					result(x, y) = result(x, b);	//pad the bottom
				}
			}
			for (size_t y = t; y <= b; y++) {
				for (size_t x = 0; x < l; x++) {
					result(x, y) = result(l, y);						//pad the left
				}
			}
			for (size_t y = t; y <= b; y++) {
				for (size_t x = r + 1; x < result.width(); ++x) {
					result(x, y) = result(r, y);		//pad the right
				}
			}
			for (size_t y = 0; y < t; y++) {
				for (size_t x = 0; x < l; x++) {
					result(x, y) = result(l, t);						//pad the top left
				}
			}
			for (size_t y = 0; y < t; y++) {
				for (size_t x = r + 1; x < result.width(); ++x) {
					result(x, y) = result(r, t);		//pad the top right
				}
			}
			for (size_t y = b + 1; y < result.height(); ++y) {
				for (size_t x = 0; x < l; x++) {
					result(x, y) = result(l, b);		//pad the bottom left
				}
			}
			for (size_t y = b + 1; y < result.height(); ++y) {
				for (size_t x = r + 1; x < result.width(); ++x) {
					result(x, y) = result(r, b);		//pad the bottom right
				}
			}
			return result;
		}

		image<T> border_remove(const size_t w) {
			return crop(w, w, X() - 2 * w, Y() - 2 * w);
		}

		/// <summary>
		/// Crops an image to a desired size
		/// </summary>
		/// <param name="x0">X boundary of the cropped image</param>
		/// <param name="y0">Y boundary of the cropped image</param>
		/// <param name="w">Width of the cropped image</param>
		/// <param name="h">Height of the cropped image</param>
		/// <returns></returns>
		image<T> crop(size_t x0, size_t y0, size_t w, size_t h) {
			if (x0 + w > width() || y0 + h > height()) {
				std::cout << "ERROR: cropped image contains an invalid region." << std::endl;
				exit(1);
			}
			image<T> result(w, h, C());								//create the output cropped image

			const size_t line_bytes = w * C() * sizeof(T);				//calculate the number of bytes in a line
			for (size_t yi = 0; yi < h; yi++) {						//for each row in the cropped image
				size_t srci = (y0 + yi) * X() * C() + x0 * C();			//calculate the source index
				size_t dsti = yi * w * C();								//calculate the destination index
				memcpy(&result.m_data[dsti], &field<T>::m_data[srci], line_bytes);	//copy the data
			}
			return result;
		}

		/// <summary>
		/// Convolves the image by a 2D mask and returns the result
		/// </summary>
		/// <param name="mask"></param>
		/// <returns></returns>
		template <typename D>
		image<T> convolve2(image<D> mask) const {
			image<T> result(X() - (mask.X() - 1), Y() - (mask.Y() - 1), C());		// output image will be smaller than the input (only valid region returned)

			const size_t width = result.width();
			const size_t height = result.height();
			const size_t channels = result.channels();
			const size_t kwidth = mask.width();
			const size_t kheight = mask.height();
			T ival, kval;
			for (size_t yi = 0; yi < height; yi++) {
				for (size_t xi = 0; xi < width; xi++) {
					for (size_t ci = 0; ci < channels; ci++) {
						T sum = static_cast<T>(0);
						for (size_t vi = 0; vi < kheight; vi++) {
							for (size_t ui = 0; ui < kwidth; ui++) {
								ival = field<T>::m_data[idx_offset(xi + ui, yi + vi, ci)];
								kval = mask(ui, vi, 0);
								sum += ival * kval;
							}
						}
						result(xi, yi, ci) = sum;
					}
				}
			}
			return result;
		}


		T normal(T x, T sigma) {
			return std::exp(-(x*x)/(2*sigma*sigma)) / std::sqrt(2.0 * std::numbers::pi * sigma * sigma);
		}

		image<T> gaussian_filter(T sigma, T sigma2=-1, bool border=true) {
			if(sigma2 < 0) sigma2 = sigma;

			int w_sigma = ceil(6 * std::max<T>(sigma, 1) + 1);					// calculate the window sizes for each blur kernel
			int w_sigma2 = ceil(6 * std::max<T>(sigma2, 1) + 1);

			image<T> blur(w_sigma, 1);							// create the images representing the blur kernels
			image<T> blur2(1, w_sigma2);

			for(int wi=0; wi < w_sigma; wi++) {				// fill each kernel with the appropriate normalize Gaussian values
				T x = -w_sigma / 2 + wi;
				blur(wi, 0) = normal(x, sigma);
			}

			for(int wi=0; wi < w_sigma2; wi++) {
				T x = -w_sigma2 / 2 + wi;
				blur2(wi, 0) = normal(x, sigma2);
			}

			std::vector<size_t> border_size = { static_cast<size_t>(w_sigma / 2), static_cast<size_t>(w_sigma2 / 2) };
			//image<T> result1 = field<T>::border(border_size);
			image<T> result1 = field<T>::Border(border_size);
			image<T> result2 = result1.convolve2(blur);
			image<T> result3 = result2.convolve2(blur2);

			return result3;		
		}


		/// <summary>
		/// Copies the non-interleaved image data to the specified pointer
		/// </summary>
		/// <param name="dest">destination of the copy (assumes the memory has been allocated)</param>
		void get_noninterleaved(T* dest) {
			//for each channel
			for (size_t y = 0; y < Y(); y++) {
				for (size_t x = 0; x < X(); x++) {
					for (size_t c = 0; c < C(); c++) {
						dest[c * Y() * X() + y * X() + x] = field<T>::m_data[idx_offset(x, y, c)];
					}
				}
			}
		}

		/* I don't know what this function is supposed to do and I didn't write it (DAVID). Currently I believe it literally does nothing.
		image<T> bin(size_t bx, size_t by) {
			size_t nbx = X() / bx;												// calculate the number of bins along x and y
			size_t nby = Y() / by;

			image<T> result(nbx, nby, C());												// generate an image to store the binned result
			result = 0;

			for(size_t yb = 0; yb < nby; yb++){
				for (size_t xb = 0; xb < nbx; xb++) {
					for (size_t yi = 0; yi < by; yi++) {
						for (size_t xi = 0; xi < bx; xi++) {
							for (size_t ci = 0; ci < C(); ci++) {
								result(xb, yb, ci) += at(xb * bx + xi, yb * by + yi, ci);
							}
						}
					}
				}
			}
			return result;
		}
		*/

		void set_noninterleaved(T* data, const size_t width, const size_t height, const size_t chan) {
			init(width, height, chan);

			//for each channel
			for (size_t y = 0; y < Y(); y++) {
				for (size_t x = 0; x < X(); x++) {
					for (size_t c = 0; c < C(); c++) {
						field<T>::m_data[idx_offset(x, y, c)] = data[c * Y() * X() + y * X() + x];
					}
				}
			}
		}

		template<typename D = T>
		void load_npy(std::string filename) {
			field<T>::template LoadNpy<D>(filename);										// load the numpy file using the tira::field class
			if (field<T>::m_shape.size() == 2)										// if the numpy array is only 2D, add a color channel of size 1
				field<T>::m_shape.push_back(1);
		}

		image<T> dist() {
			tira::image<int> boundary(width(), height());
			return dist(boundary);
		}

		image<T> sdf() {

			tira::image<int> boundary(width(), height());
			image<T> distance = dist(boundary);
			int width = distance.width();
			int height = distance.height();

			std::vector<float>SDF;
			SDF.resize(height * width);

			for (int y = 0; y < height; y++)
			{
				for (int x = 0; x < width; x++)
				{
					SDF[y * width + x] = distance(x, y);
				}

			}

			std::vector <float> frozenCells;
			frozenCells.resize(height * width);


			// initializing frozencells
			for (int i = 0; i < height; i++)
			{
				for (int j = 0; j < width; j++)
				{
					frozenCells[i * width + j] = boundary(j, i);
				}
			}

			// turn the whole input distance field to negative
			for (int i = 0; i < SDF.size(); i++) {
				SDF[i] = -1 * SDF[i];
			}

			const int row = width;
			const int nx = width - 1;
			const int ny = height - 1;
			int ix = 0, iy = 0;

			// initialize a std::tuple to store the values of frozen cell
			std::stack<std::tuple<int, int>> stack = {};

			// find the first unfrozen cell
			int gridPos = 0;

			while (frozenCells[gridPos]) {
				ix += (ix < nx ? 1 : 0);
				iy += (iy < ny ? 1 : 0);
				gridPos = row * iy + ix;
			}
			stack.push({ ix, iy });
			// a simple pixel flood
			while (stack.size()) {
				std::tuple<int, int> idsPair = stack.top();
				stack.pop();
				ix = std::get<0>(idsPair);
				iy = std::get<1>(idsPair);
				gridPos = row * iy + ix;
				if (!frozenCells[gridPos]) {
					const double val = -1.0 * SDF[gridPos];
					SDF[gridPos] = val;
					frozenCells[gridPos] = true; // freeze cell when done
					if (ix > 0) {
						stack.push({ ix - 1, iy });
					}
					if (ix < nx) {
						stack.push({ ix + 1, iy });
					}
					if (iy > 0) {
						stack.push({ ix, iy - 1 });
					}
					if (iy < ny) {
						stack.push({ ix, iy + 1 });
					}
				}
			}


			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					distance(x, y) = SDF[y * width + x];
				}
			}
			return distance;
		}

		/// <summary>
		/// Generate a color map of the image.
		/// </summary>
		/// <param name="minval"></param>
		/// <param name="maxval"></param>
		/// <param name="colormap">type of colormap to apply</param>
		/// <returns></returns>
		image<unsigned char> cmap(T minval, T maxval, const ColorMap colormap) {

			if (C() > 1) {																	// throw an exception if the image is more than one channel
				throw "Cannot create a color map from images with more than one channel!";
			}
			image<unsigned char> color_result(X(), Y(), 3);									// create the new color image
			for (size_t i = 0; i < field<T>::m_data.size(); i++) {
				cmap::colormap(field<T>::m_data[i], minval, maxval, color_result.Data()[i * 3 + 0], color_result.Data()[i * 3 + 1], color_result.Data()[i * 3 + 2], colormap);
			}
			return color_result;
		}

		/// <summary>
		/// Generate a color map of the image.
		/// </summary>
		/// <returns></returns>
		image<unsigned char> cmap(ColorMap colormap = ColorMap::Brewer) {
			return cmap(minv(), maxv(), colormap);
		}


	};	// end class image

}	// end namespace tira