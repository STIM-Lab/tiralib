#ifndef TIRA_IMAGE_H
#define TIRA_IMAGE_H

#ifdef _WIN32
#undef max
#endif

#include <vector>
#include <iostream>
#include <limits>							//use limits and remove the MIN and MAX macros
#include <typeinfo>
#include <fstream>
#include <cstring>
#include <tira/image/CImg.h>


#include <tira/filename.h>

namespace tira{
/// This static class provides the STIM interface for loading, saving, and storing 2D images.
/// Data is stored in an interleaved (BIP) format (default for saving and loading is RGB).

template <class T>
class image{

protected:

	T* img;										//pointer to the image data (interleaved RGB for color)
	size_t R[3];
	bool interleaved = true;					//by default the data is interleaved

	inline size_t X() const { return R[1]; }
	inline size_t Y() const { return R[2]; }
	inline size_t C() const { return R[0]; }

	void init(){								//initializes all variables, assumes no memory is allocated
		memset(R, 0, sizeof(size_t) * 3);		//set the resolution and number of channels to zero
		img = NULL;
	}

	void unalloc(){								//frees any resources associated with the image
		if(img)	free(img);						//if memory has been allocated, free it
	}


	void clear(){								//clears all image data
		unalloc();								//unallocate previous memory
		init();									//re-initialize the variables
	}

	void allocate(){
		unalloc();
		img = (T*) malloc( sizeof(T) * R[0] * R[1] * R[2] );	//allocate memory
		if (img == NULL) {
			std::cout << "tira::image ERROR - failed to allocate memory for image" << std::endl;
			exit(1);
		}
	}

	void allocate(size_t x, size_t y, size_t c){	//allocate memory based on the resolution
		R[0] = c; R[1] = x; R[2] = y;				//set the resolution
		allocate();									//allocate memory
	}

	inline size_t idx(size_t x, size_t y, size_t c = 0) const {
		return y * R[0] * R[1] + x * R[0] + c;
	}

	/// Returns the value for "white" based on the dynamic range (assumes white is 1.0 for floating point images)
	T white(){
		if (typeid(T) == typeid(double) || typeid(T) == typeid(float))
			return (T)1.0;
		else
			return std::numeric_limits<T>::max();
	}

public:

	size_t bytes() { return size() * sizeof(T); }

	/// Default constructor - creates an empty image object
	image(){ init(); }							//initialize all variables to zero, don't allocate any memory

	/// Constructor with a filename - loads the specified file
	image(std::string filename){				//constructor initialize the image with an image file
		init();
		load(filename);
	}

	/// Create a new image from scratch given a number of samples and channels
	image(size_t x, size_t y = 1, size_t c = 1){
		init();
		allocate(x, y, c);
	}

	/// Create a new image with the data given in 'data'
	image(T* data, size_t x, size_t y, size_t c = 1){
		init();
		allocate(x, y, c);
		memcpy(img, data, bytes());
	}

	/// Copy constructor - duplicates an image object
	image(const image<T>& I){
		init();
		allocate(I.X(), I.Y(), I.C());
		memcpy(img, I.img, bytes());
	}

	/// Destructor - clear memory
	~image(){
		free(img);
	}

	image<T>& operator=(const image<T>& I){
		if(&I == this)									//handle self-assignment
			return *this;
		init();
		allocate(I.X(), I.Y(), I.C());
		memcpy(img, I.img, bytes());
		return *this;
	}

	//determines if a filename represents a valid file format that can be loaded/saved
	static bool test_filename(std::string f) {
		filename fname = f;
		std::string ext = fname.extension();
		if (ext == "bmp" ||
			ext == "jpg" ||
			ext == "png" ||
			ext == "pbm" ||
			ext == "tif" )
			return true;
		else
			return false;
	}

	//save a Netpbm file
	void load_netpbm(std::string filename) {
		std::ifstream infile(filename.c_str(), std::ios::in | std::ios::binary);		//open an output file
		if (!infile) {
			std::cout << "Error opening input file in image::load_netpbm()" << std::endl;
			exit(1);
		}

		size_t nc;													//allocate space for the number of channels
		char format[2];												//allocate space to hold the image format tag
		infile.read(format, 2);										//read the image format tag
		infile.seekg(1, std::ios::cur);								//skip the newline character

		if (format[0] != 'P') {
			std::cout << "Error in image::load_netpbm() - file format tag is invalid: " << format[0] << format[1] << std::endl;
			exit(1);
		}
		if (format[1] == '5') nc = 1;								//get the number of channels from the format flag
		else if (format[1] == '6') nc = 3;
		else {
			std::cout << "Error in image::load_netpbm() - file format tag is invalid: " << format[0] << format[1] << std::endl;
			exit(1);
		}

		unsigned char c;								//stores a character
		while (infile.peek() == '#') {					//if the next character indicates the start of a comment
			while (true) {
				c = infile.get();
				if (c == 0x0A) break;
			}
		}

		std::string sw;									//create a string to store the width of the image
		while(true){
			c = infile.get();							//get a single character
			if (c == ' ') break;						//exit if we've encountered a space
			sw.push_back(c);							//push the character on to the string
		}
		size_t w = atoi(sw.c_str());					//convert the string into an integer

		std::string sh;
		while (true) {
			c = infile.get();
			if (c == 0x0A) break;
			sh.push_back(c);
		}

		while (true) {									//skip the maximum value
			c = infile.get();
			if (c == 0x0A) break;
		}
		size_t h = atoi(sh.c_str());					//convert the string into an integer

		allocate(w, h, nc);													//allocate space for the image
		unsigned char* buffer = (unsigned char*)malloc(w * h * nc);			//create a buffer to store the read data
		infile.read((char*)buffer, size());									//copy the binary data from the file to the image
		infile.close();														//close the file
		for (size_t n = 0; n < size(); n++) img[n] = (T)buffer[n];			//copy the buffer data into the image
		free(buffer);														//free the buffer array
	}
	

	//Copy N data points from source to dest, casting while doing so
	template<typename S, typename D>
	void type_copy(S* source, D* dest, size_t N) {
		if (typeid(S) == typeid(D))						//if both types are the same
			memcpy(dest, source, N * sizeof(S));		//just use a memcpy
		for (size_t n = 0; n < N; n++)					//otherwise, iterate through each element
			dest[n] = (D)source[n];							//copy and cast
	}
	/// Load an image from a file
	void load(std::string filename){
		//Use CImg to load the file
		cimg_library::CImg<T> cimg(filename.c_str());	//create a CImg object for the image file

		set_noninterleaved(cimg.data(), cimg.width(), cimg.height(), cimg.spectrum());
	}



	//save a Netpbm file
	void save_netpbm(std::string filename) {
		std::ofstream outfile(filename.c_str(), std::ios::out | std::ios::binary);		//open an output file
		if(!outfile) {
			std::cout << "Error generating output file in image::save_netpbm()" << std::endl;
			exit(1);
		}
		if (sizeof(T) != 1) {
			std::cout << "Error in image::save_netpbm() - data type must be 8-bit integer." << std::endl;
			exit(1);
		}
		std::string format;
		if (channels() == 1) outfile << "P5" << (char)0x0A;			//output P5 if the file is grayscale
		else if (channels() == 3) outfile << "P6" << (char)0x0A;		//output P6 if the file is color
		else {
			std::cout << "Error in image::save_netpbm() - data must be grayscale or RGB." << std::endl;
			exit(1);
		}
		size_t w = width();
		size_t h = height();
		outfile << w << " " << h << (char)0x0A;			//save the width and height
		outfile << "255" << (char)0x0A;								//output the maximum value
		outfile.write((const char*)img, size());			//write the binary data
		outfile.close();
	}

	//save a file
	void save(std::string fname){
		filename file(fname);
		if (file.extension() == "raw" || file.extension() == "") {
			std::ofstream outfile(fname.c_str(), std::ios::binary);
			outfile.write((char*)img, sizeof(T) * R[0] * R[1] * R[2]);
			outfile.close();
		}
		else {
			cimg_library::CImg<T> cimg((unsigned int)X(), (unsigned int)Y(), 1, (unsigned int)C());
			get_noninterleaved(cimg.data());
			//get_interleaved_rgb(cimg.data());
			cimg.save(fname.c_str());
		}
	}

	/// Returns an image cast to the specified format
	template<typename U>
	image<U> convert() {
		
		image<U> new_image(R[1], R[2], R[0]);					//create a new image with the destination data type

		size_t ni = R[0] * R[1] * R[2];							//calculate the number of data points in the image

		double inmax = (std::numeric_limits<T>::max)();				//get the maximum value for the input image
		double outmax = (std::numeric_limits<U>::max)();				//get the maximum value for the output image
		for (size_t i = 0; i < ni; i++) {							//for each pixel in the image
			if (img[i] > outmax) new_image(i) = outmax;			//if the source pixel is greater than the maximum destination pixel, set the output to maximum
			else new_image(i) = img[i];							//otherwise, copy the source value and cast it to the destination value		
		}
		return new_image;
	}

	void set_interleaved(T* buffer, size_t width, size_t height, size_t channels){
		allocate(width, height, channels);
		memcpy(img, buffer, bytes());
	}

	//create an image from an interleaved buffer
	void set_interleaved_rgb(T* buffer, size_t width, size_t height){
		set_interleaved(buffer, width, height, 3);
	}

	void set_interleaved_bgr(T* buffer, size_t width, size_t height){
		allocate(width, height, 3);
		T value;
		size_t i;
		for(size_t c = 0; c < C(); c++){								//copy directly
			for(size_t y = 0; y < Y(); y++){
				for(size_t x = 0; x < X(); x++){
					i = y * X() * C() + x * C() + (2-c);
					value = buffer[i];
					img[idx(x, y, c)] = value;
				}
			}
		}
	}

	void set_interleaved(T* buffer, size_t width, size_t height){
		set_interleaved_rgb(buffer, width, height);
	}

	//copies data in the given channel order as a non-interleaved image
	void set_noninterleaved(T* data, size_t width, size_t height, size_t chan) {
		allocate(width, height, chan);

		//for each channel
		for (size_t y = 0; y < Y(); y++) {
			for (size_t x = 0; x < X(); x++) {
				for (size_t c = 0; c < C(); c++) {
					img[idx(x, y, c)] = data[c * Y() * X() + y * X() + x];
				}
			}
		}
	}

	void get_interleaved_bgr(T* data){
		//for each channel
		T* source;
		if (C() == 3) {
			source = img;														//if the image has 3 channels, interleave all three
		}
		else if (C() == 4) {
			source = img + X() * Y();											//if the image has 4 channels, skip the alpha channel
		}
		else {
			throw std::runtime_error("ERROR: a BGR image must be 3 or 4 channels");	//throw an error if any other number of channels is provided
		}
		for(size_t y = 0; y < Y(); y++){
			for(size_t x = 0; x < X(); x++){
				for(size_t c = 0; c < 3; c++){
					data[y * X() * 3 + x * 3 + (2-c)] = source[y*3*R[1] + x*3 + c];
				}
			}
		}
	}

	void get_interleaved_rgb(T* data){
		memcpy(data, img, bytes());
	}

	//copies data in the given channel order as a non-interleaved image
	void get_noninterleaved(T* data){
		//for each channel
		for(size_t y = 0; y < Y(); y++){
			for(size_t x = 0; x < X(); x++){
				for(size_t c = 0; c < C(); c++){
					data[c * Y() * X() + y * X() + x] = img[idx(x, y, c)];
				}
			}
		}
	}

	/// Return an image representing a specified channel
	/// @param c is the channel to be returned
	image<T> channel(size_t c) const {		
		image<T> r(X(), Y(), 1);				//create a new image
		for(size_t x = 0; x < X(); x++){
			for(size_t y = 0; y < Y(); y++){
				r.img[r.idx(x, y, 0)] = img[idx(x, y, c)];
			}
		}
		return r;
	}

	/// Returns an std::vector containing each channel as a separate image
	std::vector< image<T> > split() const {
		std::vector< image<T> > r;			//create an image array
		r.resize(C());						//create images for each channel

		for (size_t c = 0; c < C(); c++) {	//for each channel
			r[c] = channel(c);				//copy the channel image to the array
		}
		return r;
	}

	/// Merge a series of single-channel images into a multi-channel image
	void merge(std::vector< image<T> >& list) {
		size_t x = list[0].width();				//calculate the size of the image
		size_t y = list[0].height();
		allocate(x, y, list.size());			//re-allocate the image
		for (size_t c = 0; c < list.size(); c++)		//for each channel
			set_channel(list[c].channel(0).data(), c);	//insert the channel into the output image
	}

	T& operator()(size_t x, size_t y, size_t c = 0){
		return img[idx(x, y, c)];
	}

	/// This function returns a pixel reference based on a 1D index into the image
	T& operator()(size_t i) {
		return img[i];
	}

	/// Set all elements in the image to a given scalar value

	/// @param v is the value used to set all values in the image
	void set_all(T v) {														//set all elements of the image to a given value v
		size_t N = size();
		for (size_t n = 0; n < N; n++) img[n] = v;
	}
	image<T> operator=(T v){
		set_all(v);
		return *this;
	}

	/// invert the image, given a specified maximum value (ex. maxval = 255, I' = 255 - I)
	/*image<T> invert(T maxval) {
		image<T> result(width(), height(), channels());		//create a new image
		size_t N = size();									//get the number of elements in the image
		for (size_t n = 0; n < N; n++)
			result.data()[n] = maxval - img[n];				//perform the inversion and save the result to the new image
		return result;
	}*/

	/// Stretch the contrast of the image such that the minimum and maximum intensity match the given values
	image<T> stretch(T low, T high) {
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
			result.data()[n] = desired_range * (img[n] - minval) / range + low;
		}
		return result;
	}

	/// Add a border of width w with the given value around the image
	/// @param w specifies the total size of the border
	/// @param T is the pixel value (all channels will be the same)
	image<T> border(size_t w, T value = 0) {
		image<T> result(width() + w * 2, height() + w * 2, channels());						//create an output image
		result = value;														//assign the border value to all pixels in the new image
		for (size_t y = 0; y < height(); y++) {								//for each pixel in the original image
			for (size_t x = 0; x < width(); x++) {
				size_t n = (y + w) * (width() + w * 2) + x + w;				//calculate the index of the corresponding pixel in the result image
				size_t n0 = idx(x,y);										//calculate the index for this pixel in the original image
				result.data()[n] = img[n0];									// copy the original image to the result image afer the border area
			}
		}
		return result;
	}

	/// Adds curcular padding for the specified number of pixels - in this case replicating the boundary pixels
	image<T> pad_replicate(size_t p) {
		image<T> result(width() + p * 2, height() + p * 2, channels());						//create an output image
		result = 0;
		//result = value;														//assign the border value to all pixels in the new image
		for (size_t y = 0; y < height(); y++) {								//for each pixel in the original image
			for (size_t x = 0; x < width(); x++) {
				size_t n = (y + p) * (width() + p * 2) + x + p;				//calculate the index of the corresponding pixel in the result image
				size_t n0 = idx(x, y);										//calculate the index for this pixel in the original image
				result.data()[n] = img[n0];									// copy the original image to the result image afer the border area
			}
		}
		size_t l = p;
		size_t r = p + width() - 1;
		size_t t = p;
		size_t b = p + height() - 1;
		for (size_t y = 0; y < p; y++) for (size_t x = l; x <= r; x++) result(x, y) = result(x, t);						//pad the top
		for (size_t y = b + 1; y < result.height(); y++) for (size_t x = l; x <= r; x++) result(x, y) = result(x, b);	//pad the bottom
		for (size_t y = t; y <= b; y++) for (size_t x = 0; x < l; x++) result(x, y) = result(l, y);						//pad the left
		for (size_t y = t; y <= b; y++) for (size_t x = r+1; x < result.width(); x++) result(x, y) = result(r, y);		//pad the right
		for (size_t y = 0; y < t; y++) for (size_t x = 0; x < l; x++) result(x, y) = result(l, t);						//pad the top left
		for (size_t y = 0; y < t; y++) for (size_t x = r+1; x < result.width(); x++) result(x, y) = result(r, t);		//pad the top right
		for (size_t y = b+1; y < result.height(); y++) for (size_t x = 0; x < l; x++) result(x, y) = result(l, b);		//pad the bottom left
		for (size_t y = b+1; y < result.height(); y++) for (size_t x = r + 1; x < result.width(); x++) result(x, y) = result(r, b);		//pad the bottom right
		return result;
	}

	/// Copy the given data to the specified channel

	/// @param c is the channel number that the data will be copied to
	/// @param buffer is a pointer to the image to be copied to channel c

	void set_channel(T* buffer, size_t c){
		size_t x, y;
		for(y = 0; y < Y(); y++){
			for(x = 0; x < X(); x++){
				img[idx(x, y, c)] = buffer[y * X() + x];
			}
		}
	}

	/// Set the specified channel to a constant value

	/// @param c is the channel number that the data will be copied to
	/// @param buffer is a pointer to the image to be copied to channel c

	void set_channel(T val, size_t c) {
		size_t x, y;
		for (y = 0; y < Y(); y++) {
			for (x = 0; x < X(); x++) {
				img[idx(x, y, c)] = val;
			}
		}
	}

	size_t channels() const{
		return C();
	}

	size_t width() const{
		return X();
	}

	size_t height() const{
		return Y();
	}

	T* data(){
		return img;
	}

	//returns the size (number of values) of the image
	size_t size(){ return C() * X() * Y(); }

	/// Returns the number of nonzero values
	size_t nnz(){

		size_t N = X() * Y() * C();

		size_t nz = 0;
		for(size_t n = 0; n < N; n++)
			if(img[n] != 0) nz++;

		return nz;	//return the number of nonzero pixels

	}

	//this function returns indices of pixels that have nonzero values
	std::vector<size_t> sparse_idx(){

		std::vector<size_t> s;				//allocate an array
		s.resize(nnz());					//allocate space in the array

		size_t N = size();
		//size_t C = channels();

		//T* ptr = img.data();				//get a pointer to the image data

		size_t i = 0;
		for(size_t n = 0; n < N; n++){
			if(img[n] != 0){
				s[i] = n;
				i++;
			}
		}

		return s;			//return the index list
	}


	/// Returns the maximum pixel value in the image
	T maxv(){
		T max_val = img[0];				//initialize the maximum value to the first one
		size_t N = size();	//get the number of pixels

		for (size_t n=0; n<N; n++){		//for every value

			if (img[n] > max_val){			//if the value is higher than the current max
				max_val = img[n];
			}
		}
		return max_val;
	}

	/// Returns the maximum pixel value in the image
	T minv(){
		T min_val = img[0];				//initialize the maximum value to the first one
		size_t N = size();	//get the number of pixels

		for (size_t n=0; n<N; n++){		//for every value
			if (img[n] < min_val){			//if the value is higher than the current max
				min_val = img[n];
			}
		}

		return min_val;
	}

	/// Invert an image by calculating I1 = alpha - I0, where alpha is the maximum image value
	image<T> invert(T white_val){
		size_t N = size();						//calculate the total number of values in the image
		image<T> r(X(), Y(), C());				//allocate space for the resulting image
		for(size_t n = 0; n < N; n++)
			r.img[n] = white_val - img[n];		//perform the inversion

		return r;								//return the inverted image
	}

	image<T> crop(size_t x0, size_t y0, size_t w, size_t h){
		if(x0 + w > width() || y0 + h > height()){
			std::cout<<"ERROR: cropped image contains an invalid region."<<std::endl;
			exit(1);
		}
		image<T> result(w, h, C());								//create the output cropped image

		size_t srci;
		size_t dsti;
		size_t line_bytes = w * C();							//calculate the number of bytes in a line
		for (size_t yi = 0; yi < h; yi++) {						//for each row in the cropped image
			srci = (y0 + yi) * X() * C() + x0 * C();			//calculate the source index
			dsti = yi * w * C();								//calculate the destination index
			memcpy(&result.img[dsti], &img[srci], line_bytes);	//copy the data
		}
		return result;
	}

	//crop regions given by an array of 1D index values
	std::vector< image<T> > crop_idx(size_t w, size_t h, std::vector<size_t> idx) {
		std::vector< image<T> > result(idx.size());										//create an array of image files to return
		for (size_t i = 0; i < idx.size(); i++) {										//for each specified index point
			size_t y = idx[i] / X();													//calculate the y coordinate from the 1D index (center of ROI)
			size_t x = idx[i] - y * X();												//calculate the x coordinate (center of ROI)
			y -= w / 2;																	//update x and y values to reflect the lower corner of the ROI
			x -= h / 2;
			result[i] = crop(x, y, w, h);												//get the cropped image and store it in the result array
		}
		return result;
	}

	//operator functions
	image<T> operator+(image<T> rhs) {
		size_t N = size();						//calculate the total number of values in the image
		image<T> r(X(), Y(), C());				//allocate space for the resulting image
		for (size_t n = 0; n < N; n++)
			r.img[n] = img[n] + rhs.img[n];		//perform the inversion
		return r;								//return the inverted image
	}

	image<T> srgb2lab(){
		std::cout<<"ERROR tira::image::srgb2lab - function has been broken, re-implement."<<std::endl;
		exit(1);
	}

	image<T> convolve2(image<T> mask){
		image<T> result(X() - (mask.X() - 1), Y() - (mask.Y() - 1), C());		// output image will be smaller than the input (only valid region returned)

		T sum;
		for (size_t yi = 0; yi < result.height(); yi++) {
			for (size_t xi = 0; xi < result.width(); xi++) {
				for (size_t ci = 0; ci < result.channels(); ci++) {
					sum = (T)0;
					for (size_t vi = 0; vi < mask.height(); vi++) {
						for (size_t ui = 0; ui < mask.width(); ui++) {
							sum += img[idx(xi + ui, yi + vi, ci)] * mask(ui, vi, 0);
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


	image<T> rotate(float angle, float cx, float cy){
		std::cout<<"ERROR tira::image::rotate - function has been broken, and shouldn't really be in here."<<std::endl;
		exit(1);
	}

	// leila's code for non_interleaving data in 3D
	//create an data set from an interleaved buffer
	void set_interleaved3(T* buffer, size_t width, size_t height, size_t depth, size_t channels = 3){
		std::cout<<"ERROR tira::image::set_interleaved3 - tira::image no longer supports 3D images."<<std::endl;
		exit(1);
	}

	/// Casting operator, casts every value in an image to a different data type V
	template<typename V>
	operator image<V>() {
		image<V> r(X(), Y(), C());					//create a new image
		std::copy(img, img + size(), r.data());		//copy and cast the data
		return r;									//return the new image
	}

};

};		//end namespace tira


#endif
