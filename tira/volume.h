#pragma once

#include <tira/image.h>

#include <iomanip>

// Added by HELIA (test)
#include <filesystem>
// end library test

namespace tira {
	/// This static class provides the STIM interface for loading, saving, and storing 2D images.
	/// Data is stored in an interleaved (BIP) format (default for saving and loading is RGB).

	template <class T>
	class volume : public field<T> {		// the image class extends field

	protected:

		

	
		/// <summary>
		/// Allocate space for an empty volume
		/// </summary>
		/// <param name="x">X size (fastest axis)</param>
		/// <param name="y">Y size (slower axis)</param>
		/// <param name="z">Z size (slowest axis)</param>
		/// <param name="c">Number of color channels</param>
		void init(size_t x, size_t y, size_t z, size_t c = 1) {
			field<T>::_shape.push_back(z);
			field<T>::_shape.push_back(y);
			field<T>::_shape.push_back(x);
			field<T>::_shape.push_back(c);

			field<T>::allocate();
		}

		/// <summary>
		/// Calculate the 1D array coordinate given spatial and color coordinates
		/// </summary>
		/// <param name="x">X axis position (fastest axis)</param>
		/// <param name="y">Y axis position (slower axis)</param>
		/// <param name="z">Z axis position (slowest axis)</param>
		/// <param name="c">Color channel</param>
		/// <returns></returns>
		inline size_t idx_offset(size_t x, size_t y, size_t z, size_t c = 0) const {
			return z * C() * X() * Y() + y * C() * X() + x * C() + c;		// z * C * X * Y + y * C * X + x * C + c
		}

		/// <summary>
		/// Binary XOR function
		/// </summary>
		/// <typeparam name="T"></typeparam>
		//bool XOR(bool a, bool b) {
		//	return (((a) && (!(b))) || (((!a)) && (b)));
		//}

		/// <summary>
		/// Evaluates the grid status of a given coordinate (used to generate grids)
		/// </summary>
		/// <param name="x"></param>
		/// <param name="y"></param>
		/// <param name="z"></param>
		/// <param name="boxes"></param>
		/// <returns></returns>
		bool grid(float x, float y, float z, unsigned int boxes) {
			bool x_box, y_box, z_box;
			float box_size = 1.0f / (float)boxes;
			if ((unsigned int)(z / box_size) % 2)
				z_box = false;
			else
				z_box = true;

			if ((unsigned int)(y / box_size) % 2)
				y_box = false;
			else
				y_box = true;

			if ((unsigned int)(x / box_size) % 2)
				x_box = false;
			else
				x_box = true;

			return (x_box ^ y_box) ^ z_box;
		}

		

	public:

		/// <summary>
		/// Default constructor initializes an empty volume (0x0x0 with 0 channels)
		/// </summary>
		volume() : field<T>() {}			//initialize all variables, don't allocate any memory

		/// <summary>
		/// Create a new volume from scratch given a number of samples and channels
		/// </summary>
		/// <param name="x">size of the image along the X (fast) axis</param>
		/// <param name="y">size of the image along the Y (slow) axis</param>
		/// <param name="z">size of the image along the Z (slowest) axis</param>
		/// <param name="c">number of channels</param>
		volume(size_t x, size_t y = 1, size_t z = 1, size_t c = 1) {
			init(x, y, z, c);
		}

		/// <summary>
		/// Create a new volume with the data given in a pointer 'data'
		/// </summary>
		/// <param name="data">pointer to the raw image data</param>
		/// <param name="x">size of the image along the X (fast) axis</param>
		/// <param name="y">size of the image along the Y (slow) axis</param>
		/// <param name="z">size of the image along the Z (slowest) axis</param>
		/// <param name="c">number of channels (channels are interleaved)</param>
		volume(T* data, size_t x, size_t y, size_t z, size_t c = 1) : volume(x, y, z, c) {
			memcpy(&field<T>::_data[0], data, field<T>::bytes());									// use memcpy to copy the raw data into the field array
		}

		/// <summary>
		/// Copy constructor - duplicates a volume object
		/// </summary>
		/// <param name="I">image object to be copied</param>
		volume(const volume<T>& V) : volume(V.X(), V.Y(), V.Z(), V.C()) {
			memcpy(&field<T>::_data[0], &V._data[0], field<T>::bytes());
		}

		/// <summary>
		/// Creates a new image from a stack of image files
		/// </summary>
		/// <param name="filename">Name of the file to load (uses CImg)</param>
		volume(std::string filename) {
			load(filename);
		}

		/// Destructor
		~volume() {

		}

		// access methods for the volume size and number of channels
		inline size_t Z() const { return field<T>::_shape[0]; }
		inline size_t Y() const { return field<T>::_shape[1]; }
		inline size_t X() const { return field<T>::_shape[2]; }
		inline size_t C() const { return field<T>::_shape[3]; }

		void generate_grid(unsigned int X = 32, unsigned int Y = 32, unsigned int Z = 32, unsigned int boxes = 1) {
			init(X, Y, Z, 1);

			float pixel_size[] = { 1.0f / X, 1.0f / Y, 1.0f / Z };

			//unsigned char* char_data = (unsigned char*)malloc(X * Y * Z * sizeof(unsigned char));							// allocate space for the cube
			size_t idx_z, idx_y, idx;
			float x, y, z;
			bool x_box, y_box, z_box;
			for (size_t zi = 0; zi < Z; zi++) {
				idx_z = zi * X * Y;
				z = zi * pixel_size[0];

				for (size_t yi = 0; yi < Y; yi++) {
					idx_y = idx_z + yi * X;
					y = yi * pixel_size[0];

					for (size_t xi = 0; xi < X; xi++) {
						idx = idx_y + xi;
						x = xi * pixel_size[0];

						if (grid(x, y, z, boxes)) {
							field<T>::_data[idx] = sqrt(x * x + y * y + z * z) / sqrt(3) * 255;
						}
						else {
							field<T>::_data[idx] = 0;
						}
					}
				}
			}
		}

		void generate_rgb(unsigned int X = 32, unsigned int Y = 32, unsigned int Z = 32, unsigned int boxes = 1) {
			init(X, Y, Z, 3);


			float pixel_size[] = { 1.0f / X, 1.0f / Y, 1.0f / Z };

			size_t idx_z, idx_y, idx;
			float x, y, z;
			unsigned char r;
			unsigned char g;
			unsigned char b;
			bool x_box, y_box, z_box;
			for (size_t zi = 0; zi < Z; zi++) {
				idx_z = zi * X * Y * 3;
				z = zi * pixel_size[2];

				for (size_t yi = 0; yi < Y; yi++) {
					idx_y = idx_z + yi * X * 3;
					y = yi * pixel_size[1];

					for (size_t xi = 0; xi < X; xi++) {
						idx = idx_y + xi * 3;
						x = xi * pixel_size[0];

						if (grid(x, y, z, boxes)) {
							field<T>::_data[idx + 0] = x * 255;
							field<T>::_data[idx + 1] = y * 255;
							field<T>::_data[idx + 2] = z * 255;
						}
						else {
							field<T>::_data[idx + 0] = 0;
							field<T>::_data[idx + 1] = 0;
							field<T>::_data[idx + 2] = 0;
						}
					}
				}
			}
		}

		/// <summary>
		/// Retrieve a slice of the volume as a 2D image
		/// </summary>
		/// <param name="i">Slice to return</param>
		/// <param name="z">Axis specifying the slice orientation</param>
		/// <returns></returns>
		tira::image<T> slice(size_t i, unsigned int axis = 2) {
			size_t width, height;
			size_t min_x = 0;	size_t max_x = X() - 1;		
			size_t min_y = 0;	size_t max_y = Y() - 1;		
			size_t min_z = 0;	size_t max_z = Z() - 1;		
			if (axis == 0) {
				width = Y();
				height = Z();
				min_x = i; max_x = i;
			}
			else if (axis == 1) {
				width = X();
				height = Z();
				min_y = i; max_y = i;
			}
			else {
				width = X();
				height = Y();
				min_z = i; max_z = i;
			}

			size_t sx = max_x - min_x;
			size_t sy = max_y - min_y;
			size_t sz = max_z - min_z;
			size_t sc = C();

			tira::image<T> result(width, height, C());				// create an image at the appropriate size with the correct
																	// number of color channels
			size_t idx_x, idx_y, idx_z, idx3;
			size_t idx2 = 0;
			T value;
			for (size_t zi = 0; zi <= sz; zi++) {
				idx_z = (zi + min_z) * Y() * X() * C();				// z index into the 3D volume
				for (size_t yi = 0; yi <= sy; yi++) {
					idx_y = (yi + min_y) * X() * C();				// y index into the 3D volume
					for (size_t xi = 0; xi <= sx; xi++) {
						idx_x = (xi + min_x) * C();					// x index into the 3D volume
						for (size_t ci = 0; ci < C(); ci++) {
							idx3 = idx_x + idx_y + idx_z + ci;
							value = field<T>::_data[idx3];
							result.data()[idx2] = value;
							idx2++;
						}
					}
				}
			}
			return result;
		}

		void save(std::string prefix, std::string format = "bmp") {

			for (size_t zi = 0; zi < Z(); zi++) {							// for each z slice
				
				std::stringstream ss;										// get the string representation of the current number
				ss << std::setfill('0') << std::setw(size_t(log10(Z()))) << zi;
				std::string num = ss.str();

				std::string filename = prefix + num + "." + format;			// get the file name for this slice

				tira::image<T> slice_image = slice(zi);						// get an image of the slice
				slice_image.save(filename);
			}
		}

		/// <summary>
		/// Loads a volume from a stack of images specified by a file mask
		/// </summary>
		/// <param name="file_mask"></param>
		void load(std::vector<std::string> file_names) {

			tira::image<T> test(file_names[0]);

			init(test.width(), test.height(), file_names.size(), test.channels());			// allocate space for the new shape					
	
			// copies data from each image file to the destination data
			for (size_t zi = 0; zi < Z(); zi++) {											// iterate through each file

				tira::image<T> img(file_names[zi]);											// load the image file

				for (size_t yi = 0; yi < Y(); yi++) {
					for (size_t xi = 0; xi < X(); xi++) {
						for (size_t c = 0; c < C(); c++) {
							size_t volume_index = idx_offset(xi, yi, zi, c);
							field<T>::_data[volume_index] = img(xi, yi, c);
						}
					}
				}
			}

		}

		void load(std::string file_mask) {
			if (field<T>::_shape.size() != 0) {							// if a volume is already allocated, clear it
				field<T>::_shape.clear();
				field<T>::_data.clear();
			}
		// generate a list of file names from the mask	

			std::vector<std::string> file_names;

			// ------------- using filesystem ---------------------- //
			
			// obtaining the path string from input
			std::size_t pos = file_mask.rfind("*");												// position of * in file_mask
			std::string path(file_mask.substr(0, pos));											// gets the path of the folder containing the images
			std::string ext(file_mask.substr(pos + 1, file_mask.length() - pos));				// stores the file format (jpg, bmp, ppm, ...)

			// check if the file format is valid
			if (ext != ".bmp" &&
				ext != ".png" &&
				ext != ".pbm" &&
				ext != ".tif") {

				std::cout << "File format not valid" << std::endl;

			}

			for (auto& p : std::filesystem::directory_iterator(path))							// iterate through each file and stores the ones with 
			{																					// desired extension in file_name vector
				if (p.path().extension() == ext) {
					file_names.push_back(path + p.path().filename().string());
				}
			}
			if(file_names.size() != 0)
				load(file_names);
			else {
				std::cout << "ERROR: No valid files found." << std::endl;
			}
		}

		/// <summary>
		/// Return a string for describing the volume on the command line
		/// </summary>
		/// <returns>A string summarizing the volume</returns>
		std::string str() {
			std::stringstream ss;
			ss << "X: " << X() << "     Y: " << Y() << "     Z: " << Z() << std::endl;
			ss << "Colors: " << C() << std::endl;
			ss << "data type (bytes): " << sizeof(T) << std::endl;
			return ss.str();
		}


	};	// end class volume

	

}	// end namespace tira