#pragma once

#include "extern/CImg.h"

#include <tira/field.h>


namespace tira {
	/// This static class provides the STIM interface for loading, saving, and storing 2D images.
	/// Data is stored in an interleaved (BIP) format (default for saving and loading is RGB).

	template <class T>
	class image : public field<T>{		// the image class extends field

	protected:
		inline size_t C() const { return field<T>::_shape[0]; }
		inline size_t X() const { return field<T>::_shape[1]; }
		inline size_t Y() const { return field<T>::_shape[2]; }

		/// <summary>
		/// Allocate space for an empty image
		/// </summary>
		/// <param name="x">Width of the image (fast axis)</param>
		/// <param name="y">Height of the image (slow axis)</param>
		/// <param name="c">Number of color channels</param>
		void init(size_t x, size_t y, size_t c = 1) {
			field<T>::_shape.push_back(c);
			field<T>::_shape.push_back(x);
			field<T>::_shape.push_back(y);

			field<T>::allocate();
		}

		/// <summary>
		/// Calculates the 1D array offset given spatial and color coordinates
		/// </summary>
		/// <param name="x">X coordinate of the image (fast axis)</param>
		/// <param name="y">Y coordinate of the image (slow axis)</param>
		/// <param name="c">Color channel</param>
		/// <returns></returns>
		inline size_t idx_offset(size_t x, size_t y, size_t c = 0) const {
			return y * C() * X() + x * C() + c;		// y * C * X + x * C + c
		}
		
		
	public:

		/// <summary>
		/// Default constructor initializes an empty image (0x0 with 0 channels)
		/// </summary>
		image() : field<T>() {}			//initialize all variables, don't allocate any memory

		/// <summary>
		/// Create a new image from scratch given a number of samples and channels
		/// </summary>
		/// <param name="x">size of the image along the X (fast) axis</param>
		/// <param name="y">size of the image along the Y (slow) axis</param>
		/// <param name="c">number of channels</param>
		image(size_t x, size_t y = 1, size_t c = 1) {
			init(x, y, c);
		}

		/// <summary>
		/// Create a new image with the data given in 'data'
		/// </summary>
		/// <param name="data">pointer to the raw image data</param>
		/// <param name="x">size of the image along the X (fast) axis</param>
		/// <param name="y">size of the image along the Y (slow) axis</param>
		/// <param name="c">number of channels (channels are interleaved)</param>
		image(T* data, size_t x, size_t y, size_t c = 1) : image(x, y, c) {
			memcpy(&field<T>::_data[0], data, field<T>::bytes());									// use memcpy to copy the raw data into the field array
		}

		/// <summary>
		/// Copy constructor - duplicates an image object
		/// </summary>
		/// <param name="I">image object to be copied</param>
		image(const image<T>& I) : image(I.X(), I.Y(), I.C()) {
			memcpy(&field<T>::_data[0], &I._data[0], field<T>::bytes());
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
		size_t width() {
			return X();
		}

		/// <summary>
		/// Return the height of the image
		/// </summary>
		/// <returns>Height of the image (slow dimension)</returns>
		size_t height() {
			return Y();
		}

		/// <summary>
		/// Returns the number of channels in the iamge
		/// </summary>
		/// <returns>Number of channels</returns>
		size_t channels() {
			return C();
		}

		/// <summary>
		/// Returns the number of raw values in the image
		/// </summary>
		/// <returns></returns>
		size_t size() { return C() * X() * Y(); }

		/// <summary>
		/// Assignment operator, assigns one image to another
		/// </summary>
		/// <param name="I"></param>
		/// <returns>pointer to the current object</returns>
		image<T>& operator=(const image<T>& I) {
			if (&I == this)													// handle self-assignment
				return *this;
			init(I.X(), I.Y(), I.C());										// allocate space and set shape variables
			memcpy(&field<T>::_data[0], &I._data[0], field<T>::bytes());		// copy the data
			return *this;													// return a pointer to the current object
		}

		/// <summary>
		/// Set all values in the image to a single constant
		/// </summary>
		/// <param name="v">Constant that all elements will be set to</param>
		image<T> operator=(T v) {														//set all elements of the image to a given value v
			size_t N = size();
			std::fill(field<T>::_data.begin(), field<T>::_data.end(), v);
			return *this;
		}

		/// <summary>
		/// Adds two images
		/// </summary>
		/// <param name="rhs"></param>
		/// <returns></returns>
		image<T> operator+(image<T> rhs) {
			size_t N = size();										//calculate the total number of values in the image
			image<T> r(X(), Y(), C());								//allocate space for the resulting image
			for (size_t n = 0; n < N; n++)
				r._data[n] = field<T>::_data[n] + rhs._data[n];		//perform the inversion
			return r;												//return the inverted image
		}

		T& operator()(size_t x, size_t y, size_t c = 0) {
			return field<T>::_data[idx_offset(x, y, c)];
		}

		/// <summary>
		/// Cast data types
		/// </summary>
		/// <typeparam name="V"></typeparam>
		template<typename V>
		operator image<V>() {
			image<V> r(X(), Y(), C());					//create a new image
			std::copy(&field<T>::_data[0], &field<T>::_data[0] + size(), r.data());		//copy and cast the data
			return r;									//return the new image
		}

		/// <summary>
		/// Save the image as a file (uses the CImg library)
		/// </summary>
		/// <param name="fname">name of the file</param>
		void save(std::string fname) {
			cimg_library::CImg<T> cimg((unsigned int)X(), (unsigned int)Y(), 1, (unsigned int)C());
			get_noninterleaved(cimg.data());
			cimg.save(fname.c_str());
		}

		/// <summary>
		/// Load an image file (uses the CImg library)
		/// </summary>
		/// <param name="filename">Name of the file to load</param>
		void load(std::string filename) {
			cimg_library::CImg<T> cimg(filename.c_str());									//create a CImg object for the image file
			set_noninterleaved(cimg.data(), cimg.width(), cimg.height(), cimg.spectrum());
		}

		/// <summary>
		/// Returns a channel of the current image as an independent one-channel image
		/// </summary>
		/// <param name="c">channel to return</param>
		/// <returns></returns>
		image<T> channel(size_t c) const {
			image<T> r(X(), Y());											//create a new single-channel image
			for (size_t x = 0; x < X(); x++) {
				for (size_t y = 0; y < Y(); y++) {
					r._data[r.idx_offset(x, y, 0)] = field<T>::_data[idx_offset(x, y, c)];
				}
			}
			return r;
		}

		/// <summary>
		/// Set the specified image channel to the data provided in src
		/// </summary>
		/// <param name="src">raw data used as the channel source</param>
		/// <param name="c">channel that the raw data is assigned to</param>
		void channel(T* src, size_t c) {
			size_t x, y;
			for (y = 0; y < Y(); y++) {
				for (x = 0; x < X(); x++) {
					field<T>::_data[idx_offset(x, y, c)] = src[y * X() + x];
				}
			}
		}

		/// <summary>
		/// Set the all pixels in the specified channel to a single value
		/// </summary>
		/// <param name="val">Value that all pixels in the channel will be set to</param>
		/// <param name="c">Channel to be modified</param>
		void channel(T val, size_t c) {
			size_t x, y;
			for (y = 0; y < Y(); y++) {
				for (x = 0; x < X(); x++) {
					field<T>::_data[idx_offset(x, y, c)] = val;
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
				channel(&list[c].channel(0)._data[0], c);	//insert the channel into the output image
		}

		/// <summary>
		/// Returns a pointer to the raw image data
		/// </summary>
		/// <returns>Pointer to the raw image data</returns>
		T* data() {
			return &field<T>::_data[0];
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
				if (field<T>::_data[n] != background) num_foreground++;

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

			size_t N = size();

			size_t i = 0;
			for (size_t n = 0; n < N; n++) {
				if (field<T>::_data[n] != background) {
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
			T max_val = field<T>::_data[0];				//initialize the maximum value to the first one
			size_t N = size();	//get the number of pixels

			for (size_t n = 0; n < N; n++) {		//for every value

				if (field<T>::_data[n] > max_val) {			//if the value is higher than the current max
					max_val = field<T>::_data[n];
				}
			}
			return max_val;
		}

		/// <summary>
		/// Returns the minimum value in the image
		/// </summary>
		/// <returns></returns>
		T minv() {
			T min_val = field<T>::_data[0];				//initialize the maximum value to the first one
			size_t N = size();	//get the number of pixels

			for (size_t n = 0; n < N; n++) {		//for every value
				if (field<T>::_data[n] < min_val) {			//if the value is higher than the current max
					min_val = field<T>::_data[n];
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

			size_t N = size();						//get the number of values in the image
			T range = maxval - minval;			//calculate the current range of the image
			T desired_range = high - low;		//calculate the desired range of the image
			for (size_t n = 0; n < N; n++) {		//for each element in the image
				result._data[n] = desired_range * (field<T>::_data[n] - minval) / range + low;
			}
			return result;
		}

		/// <summary>
		/// Add a border to an image
		/// </summary>
		/// <param name="w">Width of the border in pixels</param>
		/// <param name="value">Value used to generate the border</param>
		/// <returns></returns>
		image<T> border(size_t w, T value = 0) {
			image<T> result(width() + w * 2, height() + w * 2, channels());						//create an output image
			result = value;														//assign the border value to all pixels in the new image
			for (size_t y = 0; y < height(); y++) {								//for each pixel in the original image
				for (size_t x = 0; x < width(); x++) {
					size_t n = (y + w) * (width() + w * 2) + x + w;				//calculate the index of the corresponding pixel in the result image
					size_t n0 = idx_offset(x, y);										//calculate the index for this pixel in the original image
					result._data[n] = field<T>::_data[n0];									// copy the original image to the result image afer the border area
				}
			}
			return result;
		}

		/// <summary>
		/// Generates a border by replicating edge pixels
		/// </summary>
		/// <param name="p">Width of the border (padding) to create</param>
		/// <returns></returns>
		image<T> border_replicate(size_t w) {
			image<T> result(width() + w * 2, height() + w * 2, channels());						//create an output image
			result = 0;
			//result = value;														//assign the border value to all pixels in the new image
			for (size_t y = 0; y < height(); y++) {								//for each pixel in the original image
				for (size_t x = 0; x < width(); x++) {
					size_t n = (y + w) * (width() + w * 2) + x + w;				//calculate the index of the corresponding pixel in the result image
					size_t n0 = idx_offset(x, y);										//calculate the index for this pixel in the original image
					result.data()[n] = field<T>::_data[n0];									// copy the original image to the result image afer the border area
				}
			}
			size_t l = w;
			size_t r = w + width() - 1;
			size_t t = w;
			size_t b = w + height() - 1;
			for (size_t y = 0; y < w; y++) for (size_t x = l; x <= r; x++) result(x, y) = result(x, t);						//pad the top
			for (size_t y = b + 1; y < result.height(); y++) for (size_t x = l; x <= r; x++) result(x, y) = result(x, b);	//pad the bottom
			for (size_t y = t; y <= b; y++) for (size_t x = 0; x < l; x++) result(x, y) = result(l, y);						//pad the left
			for (size_t y = t; y <= b; y++) for (size_t x = r + 1; x < result.width(); x++) result(x, y) = result(r, y);		//pad the right
			for (size_t y = 0; y < t; y++) for (size_t x = 0; x < l; x++) result(x, y) = result(l, t);						//pad the top left
			for (size_t y = 0; y < t; y++) for (size_t x = r + 1; x < result.width(); x++) result(x, y) = result(r, t);		//pad the top right
			for (size_t y = b + 1; y < result.height(); y++) for (size_t x = 0; x < l; x++) result(x, y) = result(l, b);		//pad the bottom left
			for (size_t y = b + 1; y < result.height(); y++) for (size_t x = r + 1; x < result.width(); x++) result(x, y) = result(r, b);		//pad the bottom right
			return result;
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

			size_t srci;
			size_t dsti;
			size_t line_bytes = w * C() * sizeof(T);				//calculate the number of bytes in a line
			for (size_t yi = 0; yi < h; yi++) {						//for each row in the cropped image
				srci = (y0 + yi) * X() * C() + x0 * C();			//calculate the source index
				dsti = yi * w * C();								//calculate the destination index
				memcpy(&result._data[dsti], &field<T>::_data[srci], line_bytes);	//copy the data
			}
			return result;
		}

		/// <summary>
		/// Convolves the image by a 2D mask and returns the result
		/// </summary>
		/// <param name="mask"></param>
		/// <returns></returns>
		image<T> convolve2(image<T> mask) {
			image<T> result(X() - (mask.X() - 1), Y() - (mask.Y() - 1), C());		// output image will be smaller than the input (only valid region returned)

			T sum;
			for (size_t yi = 0; yi < result.height(); yi++) {
				for (size_t xi = 0; xi < result.width(); xi++) {
					for (size_t ci = 0; ci < result.channels(); ci++) {
						sum = (T)0;
						for (size_t vi = 0; vi < mask.height(); vi++) {
							for (size_t ui = 0; ui < mask.width(); ui++) {
								sum += field<T>::_data[idx_offset(xi + ui, yi + vi, ci)] * mask(ui, vi, 0);
							}
						}
						result(xi, yi, ci) = sum;
					}
				}
			}
			return result;
			//std::cout<<"ERROR tira::image::convolve2 - function has been broken, and shouldn't really be in here."<<std::endl;
			//exit(1);
		}

		/// <summary>
		/// Copies the non-interleaved image data to the specified pointer
		/// </summary>
		/// <param name="data">destination of the copy (assumes the memory has been allocated)</param>
		void get_noninterleaved(T* dest) {
			//for each channel
			for (size_t y = 0; y < Y(); y++) {
				for (size_t x = 0; x < X(); x++) {
					for (size_t c = 0; c < C(); c++) {
						dest[c * Y() * X() + y * X() + x] = field<T>::_data[idx_offset(x, y, c)];
					}
				}
			}
		}

		void set_noninterleaved(T* data, size_t width, size_t height, size_t chan) {
			init(width, height, chan);

			//for each channel
			for (size_t y = 0; y < Y(); y++) {
				for (size_t x = 0; x < X(); x++) {
					for (size_t c = 0; c < C(); c++) {
						field<T>::_data[idx_offset(x, y, c)] = data[c * Y() * X() + y * X() + x];
					}
				}
			}
		}


	};	// end class image

}	// end namespace tira