#pragma once

#include <tira/image.h>

#include <iomanip>
#include <filesystem>


namespace tira {
	/// This static class provides the STIM interface for loading, saving, and storing 2D images.
	/// Data is stored in an interleaved (BIP) format (default for saving and loading is RGB).

	template <class T>
	class volume : public field<T> {		// the image class extends field

	protected:

		std::vector<double> _spacing;

	
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

			field<T>::_allocate();

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

		T& at(size_t x, size_t y, size_t z, size_t c = 0) {
			return field<T>::_data[idx_offset(x, y, z, c)];
		}

		tira::volume<float> _dist(tira::volume<int>& binary_boundary) {

			// create registers to hold the size
			int w = X();
			int h = Y();
			int l = Z();

			// resize the binary boundary grid
			binary_boundary.resize(field<T>::_shape);

			std::vector<std::tuple<int, int, int>> neighbors;				// vector stores a template for 4-connected indices
			neighbors.emplace_back(0, 0, 1);
			neighbors.emplace_back(0, 1, 0);
			neighbors.emplace_back(1, 0, 0);
			neighbors.emplace_back(-1, 0, 0);
			neighbors.emplace_back(0, -1, 0);
			neighbors.emplace_back(0, 0, -1);



			// indentifying boundary cells
			for (int y = 1; y < h - 1; y++) {						// for every row in the image
				for (int x = 1; x < w - 1; x++) {					// for every column in the image
					for (int z = 1; z < l - 1; z++) {					// for every length in the image
						for (int k = 0; k < neighbors.size(); k++) {		// for every neighbor

							int nx = x + get<0>(neighbors[k]);				// calculate the x coordinate of the neighbor cell
							int ny = y + get<1>(neighbors[k]);				// calculate the y coordinate of the neighbor cell
							int nz = z + get<2>(neighbors[k]);				// calculate the z coordinate of the neighbor cell

							if (at(x, y, z) * at(nx, ny, nz) <= 0) {				// if the product of the current cell and neighboring cell is negative (it is a boundary cell)
								binary_boundary(x, y, z) = 1;					// this cell is a boundary cell
								binary_boundary(nx, ny, nz) = 1;				// the neighboring cell is ALSO a boundary cell
							}
						}
					}
				}
			}


			tira::volume<float> dist(w, h, l);										// create an image to store the distance field
			dist = 9999.0f;																// initialize the distance field to a very large value



			// calculate the distance for all boundary cells to the contour
			for (int y = 1; y < h - 1; y++) {											// for every row in the image
				for (int x = 1; x < w - 1; x++) {										// for every column in the image
					for (int z = 1; z < l - 1; z++) {										// for every length in the image
						if (binary_boundary(x, y, z)) {											// if the pixel (x, y,z) is in the boundary
							for (int k = 0; k < neighbors.size(); k++) {						// for every neighbor

								int nx = x + get<0>(neighbors[k]);								// calculate the x coordinate of the neighbor cell
								int ny = y + get<1>(neighbors[k]);								// calculate the y coordinate of the neighbor cell
								int nz = z + get<2>(neighbors[k]);								// calculate the z coordinate of the neighbor cell
								if (binary_boundary(nx, ny, nz)) {												// if the neighboring cell (nx, ny) is ALSO in the boundary
									float da = (abs(at(x, y, z))) / (abs(at(nx, ny, nz) - at(x, y, z)));			// calculate distance from pixel(x,y) to contour da
									float db = (abs(at(nx, ny, nz))) / (abs(at(nx, ny, nz) - at(x, y, z)));			// calculate distance from neighbor to contour db
									dist(x, y, z) = std::min(dist(x, y, z), da);									// minimum between distance and large boundary value of pixel (x,y)
									dist(nx, ny, nz) = std::min(dist(nx, ny, nz), db);								// minimum between distance and large boundary value of neighbor (nx, ny)
								}
							}
						}
					}
				}
			}


			// at this point, the distance image has the correct distances in cells surround the boundary

			std::vector<float> dist1d;											// create a 1d representation of the distance field
			dist1d.resize(h * w * l);

			int row = w;

			// copy the distance field into the 1d distance field
			for (int y = 0; y < h; y++) {
				for (int x = 0; x < w; x++) {
					for (int z = 0; z < l; z++) {
						dist1d[z * (w * h) + y * w + x] = dist(x, y, z);
					}
				}
			}

			// initializing fast sweeiping algorithm 
			const int NSweeps = 8;

			//// sweep directions { start, end, step }
			const int dirX[NSweeps][3] = { {0, w - 1, 1} , {w - 1, 0, -1}, {w - 1, 0, -1}, {0, w - 1, 1}, {0, w - 1, 1} , {w - 1, 0, -1}, {w - 1, 0, -1}, {0, w - 1, 1} };
			const int dirY[NSweeps][3] = { {0, h - 1, 1}, {0, h - 1, 1}, {h - 1, 0, -1}, {h - 1, 0, -1}, {0, h - 1, 1}, {0, h - 1, 1}, {h - 1, 0, -1}, {h - 1, 0, -1} };
			const int dirZ[NSweeps][3] = { {0, l - 1, 1}, {0, l - 1, 1}, {0, l - 1, 1}, {0, l - 1, 1}, {l - 1, 0, -1}, {l - 1, 0, -1}, {l - 1, 0, -1}, {l - 1, 0, -1} };
			double aa[3], tmp, eps = 1e-6;
			double d_new, a, b;
			int s, ix, iy, iz, gridPos;
			const double dx = 1.0, f = 1.0;

			for (s = 0; s < NSweeps; s++) {

				for (iy = dirY[s][0]; dirY[s][2] * iy <= dirY[s][1]; iy += dirY[s][2]) {
					for (ix = dirX[s][0]; dirX[s][2] * ix <= dirX[s][1]; ix += dirX[s][2]) {
						for (iz = dirZ[s][0]; dirZ[s][2] * iz <= dirZ[s][1]; iz += dirZ[s][2]) {

							gridPos = ((iz * h + iy) * row + ix);

							if (iy == 0 || iy == (h - 1)) {                    // calculation for ymin
								if (iy == 0) {
									aa[1] = dist1d[(iz * h + (iy + 1)) * row + ix];
								}
								if (iy == (h - 1)) {
									aa[1] = dist1d[(iz * h + (iy - 1)) * row + ix];
								}
							}
							else {
								aa[1] = dist1d[(iz * h + (iy - 1)) * row + ix] < dist1d[(iz * h + (iy + 1)) * row + ix] ? dist1d[(iz * h + (iy - 1)) * row + ix] : dist1d[(iz * h + (iy + 1)) * row + ix];
							}

							if (ix == 0 || ix == (w - 1)) {                    // calculation for xmin
								if (ix == 0) {
									aa[0] = dist1d[(iz * h + iy) * row + (ix + 1)];
								}
								if (ix == (w - 1)) {
									aa[0] = dist1d[(iz * h + iy) * row + (ix - 1)];
								}
							}
							else {
								aa[0] = dist1d[(iz * h + iy) * row + (ix - 1)] < dist1d[(iz * h + iy) * row + (ix + 1)] ? dist1d[(iz * h + iy) * row + (ix - 1)] : dist1d[(iz * h + iy) * row + (ix + 1)];
							}

							if (iz == 0 || iz == (l - 1)) {                    // calculation for Zmin
								if (iz == 0) {
									aa[2] = dist1d[((iz + 1) * h + iy) * row + ix];
								}
								if (iz == (l - 1)) {
									aa[2] = dist1d[((iz - 1) * h + iy) * row + ix];
								}
							}
							else {
								aa[2] = dist1d[((iz - 1) * h + iy) * row + ix] < dist1d[((iz + 1) * h + iy) * row + ix] ? dist1d[((iz - 1) * h + iy) * row + ix] : dist1d[((iz + 1) * h + iy) * row + ix];
							}


							// simple bubble sort
							if (aa[0] > aa[1]) { tmp = aa[0]; aa[0] = aa[1]; aa[1] = tmp; }
							if (aa[1] > aa[2]) { tmp = aa[1]; aa[1] = aa[2]; aa[2] = tmp; }
							if (aa[0] > aa[1]) { tmp = aa[0]; aa[0] = aa[1]; aa[1] = tmp; }

							double d_curr = aa[0] + dx * f; // just a linear equation with the first neighbor value
							double d_new;
							if (d_curr <= (aa[1] + eps)) {
								d_new = d_curr; // accept the solution
							}
							else {
								// quadratic equation with coefficients involving 2 neighbor values aa
								double a = 2.0;
								double b = -2.0 * (aa[0] + aa[1]);
								double c = aa[0] * aa[0] + aa[1] * aa[1] - dx * dx * f * f;
								double D = sqrt(b * b - 4.0 * a * c);
								// choose the minimal root
								d_curr = ((-b + D) > (-b - D) ? (-b + D) : (-b - D)) / (2.0 * a);

								if (d_curr <= (aa[2] + eps))
									d_new = d_curr; // accept the solution
								else {
									// quadratic equation with coefficients involving all 3 neighbor values aa
									a = 3.0;
									b = -2.0 * (aa[0] + aa[1] + aa[2]);
									c = aa[0] * aa[0] + aa[1] * aa[1] + aa[2] * aa[2] - dx * dx * f * f;
									D = sqrt(b * b - 4.0 * a * c);
									// choose the minimal root
									d_new = ((-b + D) > (-b - D) ? (-b + D) : (-b - D)) / (2.0 * a);
								}
							}
							// update if d_new is smaller
							dist1d[gridPos] = dist1d[gridPos] < d_new ? dist1d[gridPos] : d_new;

						}
					}
				}
			}

			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					for (int z = 0; z < l; z++)
					{
						dist(x, y, z) = dist1d[z * (w * h) + y * w + x];
					}
				}

			}


			return dist;
		}

	public:

		/// <summary>
		/// Default constructor initializes an empty volume (0x0x0 with 0 channels)
		/// </summary>
		volume() : field<T>() {
			_spacing = { 1.0f, 1.0f, 1.0f };
		}			//initialize all variables, don't allocate any memory

		/// <summary>
		/// Create a new volume from scratch given a number of samples and channels
		/// </summary>
		/// <param name="x">size of the image along the X (fast) axis</param>
		/// <param name="y">size of the image along the Y (slow) axis</param>
		/// <param name="z">size of the image along the Z (slowest) axis</param>
		/// <param name="c">number of channels</param>
		volume(size_t x, size_t y, size_t z, size_t c = 1) : volume() {
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

		volume(std::vector<size_t> shape) : field<T>(shape) {
			_spacing = { 1.0f, 1.0f, 1.0f };
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
			_spacing = { 1.0f, 1.0f, 1.0f };
		}

		/// Destructor
		~volume() {

		}

		// access methods for the volume size and number of channels
		inline size_t Z() const { return field<T>::_shape[0]; }
		inline size_t Y() const { return field<T>::_shape[1]; }
		inline size_t X() const { return field<T>::_shape[2]; }
		inline size_t C() const { return field<T>::_shape[3]; }

		inline double dx() const { return _spacing[0]; }
		inline double dy() const { return _spacing[1]; }
		inline double dz() const { return _spacing[2]; }

		inline double sx() const { return X() * _spacing[0]; }
		inline double sy() const { return Y() * _spacing[1]; }
		inline double sz() const { return Z() * _spacing[2]; }
		inline double smax() const { return std::max(sx(), std::max(sy(), sz())); }

		/// <summary>
		/// Returns a channel of the current image as an independent one-channel image
		/// </summary>
		/// <param name="c">channel to return</param>
		/// <returns></returns>
		volume<T> channel(const size_t c) const {
			volume<T> r(X(), Y(), Z());							// create a new single-channel image
			T* dp = r.data();									// create a pointer to the raw result
			size_t S = X() * Y() * X();							// number of voxels
			for (size_t i = 0; i < S; i++) {					// for each voxel
				dp[i] = field<T>::_data[i * C() + c];			// copy the value corresponding to the desired channel
			}
			return r;
		}

		T minv() {
			T m = *std::min_element(field<T>::_data.begin(), field<T>::_data.end());
			return m;
		}

		T maxv() {
			T m = *std::max_element(field<T>::_data.begin(), field<T>::_data.end());
			return m;
		}

		void spacing(double x, double y, double z) {
			std::vector<double> spacing = { x, y, z };
			_spacing = spacing;
		}

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


		//calculate gradient along dx (3D)
		tira::volume<float> gradient_dx()
		{
			tira::volume<float>output(X(), Y(), Z());

			for (int j = 0; j < X(); j++)
			{
				for (int i = 0; i < Y(); i++)
				{
					for (int k = 0; k < Z(); k++)
					{
						int j_left = j - 1;
						int j_right = j + 1;
						if (j_left < 0)
						{
							j_left = 0;
							j_right = 1;
							double dist_grad = at(j_right, i, k) - at(j_left, i, k);

						}

						else if (j_right >= X())
						{
							j_right = X() - 1;
							j_left = j_right - 1;
							double dist_grad = at(j_right,i, k) - at(j_left, i, k);

						}


						double dist_grad = (at(j_right, i, k) - at(j_left, i, k)) / 2.0f;

						output(j, i, k) = dist_grad;
					}
				}
			}

			return output;
		}

		// calculate gradient along dy(3D)
		tira::volume<float > gradient_dy()
		{
			tira::volume<float> output(X(), Y(), Z());

			for (int j = 0; j < X(); j++)
			{
				for (int i = 0; i < Y(); i++)
				{
					for (int k = 0; k < Z(); k++)
					{
						int i_left = i - 1;
						int i_right = i + 1;
						if (i_left < 0)
						{
							i_left = 0;
							i_right = 1;
							double dist_grad = at(j, i_right, k) - at(j, i_left, k);

						}

						else if (i_right >= Y())
						{
							i_right = Y() - 1;
							i_left = i_right -1 ;
							double dist_grad = at(j, i_right, k) - at(j, i_left, k);

						}


						double dist_grad = (at(j, i_right, k) - at(j, i_left, k)) / 2.0f;

						output(j, i, k) = dist_grad;
					}
				}
			}

			return output;
		}

		tira::volume<float>gradient_dz()
		{
			//volume<T>dxdy(X(), Y(), Z());

			tira::volume<float>output(X(), Y(), Z());

			for (int j = 0; j < X(); j++)
			{
				for (int i = 0; i < Y(); i++)
				{
					for (int k = 0; k < Z(); k++)
					{
						int k_left = k - 1;
						int k_right = k + 1;
						if (k_left < 0)
						{
							k_left = 0;
							k_right = 1;
							double dist_grad = at(j, i, k_right) - at(j, i, k_left);
						}

						else if (k_right >= Z())
						{
							k_right = Z() - 1;
							k_left = k_right - 1;
							double dist_grad = at(j, i, k_right) - at(j, i, k_left);
						}


						double dist_grad = (at(j, i, k_right) - at(j, i, k_left)) / 2.0f;
						output(j, i, k) = dist_grad;
					}
				}
			}

			return output;
		}


		/// <summary>
		/// Convolves the image by a 3D mask and returns the result
		/// </summary>
		/// <param name="mask"></param>
		/// <returns></returns>

		tira::volume<float>convolve3D(tira::volume<float> mask) {

			tira::volume<float>result(X() - (mask.X() - 1), Y() - (mask.Y() - 1), Z() - (mask.Z() - 1));		// output image will be smaller than the input (only valid region returned)

			float sum;
			for (size_t yi = 0; yi < result.Y(); yi++) {
				for (size_t xi = 0; xi < result.X(); xi++) {
					for (size_t zi = 0; zi < result.Z(); zi++) {
						sum = 0;
						for (size_t vi = 0; vi < mask.Y(); vi++) {
							for (size_t ui = 0; ui < mask.X(); ui++) {
								for (size_t wi = 0; wi < mask.Z(); wi++) {

									sum += at((xi + ui), (yi + vi), (zi + wi)) * mask(ui, vi, wi);
									
								}
							}
						}
						
						result(xi, yi, zi) = sum;
					}

				}
			}

			return result;
		}


		/// <summary>
		/// Convolves the image by a 1D mask and returns the result
		/// </summary>
		/// <param name="mask"></param>
		/// <returns></returns>
		
		tira::volume<float> convolve3D_separate(float* K, int k_size) {
	
			tira::volume<float>result_x(X(), Y() - (k_size - 1), Z());		// output image will be smaller than the input (only valid region returned)

			float sum;
			for (size_t yi = 0; yi < result_x.Y(); yi++) {
				for (size_t xi = 0; xi < result_x.X(); xi++) {
					for (size_t zi = 0; zi < result_x.Z(); zi++) {
						sum = 0;
						for (size_t ui = 0; ui < k_size; ui++) {
							sum += at(xi, (yi + ui), zi) * K[ui];
						}
						result_x(xi, yi, zi) = sum;
					}
				}
			}

			tira::volume<float>result_y(result_x.X() - (k_size - 1), result_x.Y(), result_x.Z());		// output image will be smaller than the input (only valid region returned)

			float sum1;
			for (size_t yi = 0; yi < result_y.Y(); yi++) {
				for (size_t xi = 0; xi < result_y.X(); xi++) {
					for (size_t zi = 0; zi < result_y.Z(); zi++) {
						sum1 = 0;
						for (size_t ui = 0; ui < k_size; ui++) {
							sum1 += result_x((xi + ui), yi, zi) * K[ui];
						}
						result_y(xi, yi, zi) = sum1;
					}
				}
			}

			tira::volume<float>result_z(result_y.X(), result_y.Y(), result_y.Z() - (k_size - 1));		// output image will be smaller than the input (only valid region returned)

			float sum2;
			for (size_t yi = 0; yi < result_z.Y(); yi++) {
				for (size_t xi = 0; xi < result_z.X(); xi++) {
					for (size_t zi = 0; zi < result_z.Z(); zi++) {
						sum2 = 0;
						for (size_t ui = 0; ui < k_size; ui++) {
							sum2 += result_y(xi, yi, (zi + ui)) * K[ui];
						}
						result_z(xi, yi, zi) = sum2;
					}
				}
			}

			return result_z;

		}


		/// <summary>
		/// Add a border to an image
		/// </summary>
		/// <param name="w">Width of the border in pixels</param>
		/// <param name="value">Value used to generate the border</param>
		/// <returns></returns>
		
		tira::volume<float>border_3D(size_t w) {

			tira::volume<float> result(X() + w * 2, Y() + w * 2, Z() + w * 2);

			result = 30;
			//result = value;														//assign the border value to all pixels in the new image
			for (size_t y = 0; y < Y(); y++) {								//for each pixel in the original image
				for (size_t x = 0; x < X(); x++) {
					for (size_t z = 0; z < Z(); z++) {
						size_t n = ((z + w) * (X() + w * 2) * (Y() + w * 2)) + ((y + w) * (X() + w * 2)) + (x + w);				//calculate the index of the corresponding pixel in the result image
						size_t n0 = idx_offset(x, y, z);										//calculate the index for this pixel in the original image
						result.data()[n] = field<T>::_data[n0];									// copy the original image to the result image afer the border area
					}
				}
			}
			return result;

		}
		
		/// <summary>
		/// Generates a border by replicating edge pixels
		/// </summary>
		/// <param name="p">Width of the border (padding) to create</param>
		/// <returns></returns>


		tira::volume<float>border_replicate_3D(size_t w) {

			tira::volume<float> result(X() + w * 2, Y() + w * 2, Z() + w * 2);

			result = 0;
			//result = value;														//assign the border value to all pixels in the new image
			for (size_t y = 0; y < Y(); y++) {								//for each pixel in the original image
				for (size_t x = 0; x < X(); x++) {
					for (size_t z = 0; z < Z(); z++) {
						size_t n = ((z + w) * (X() + w * 2) * (Z() + w * 2)) + ((y + w) * (Z() + w * 2)) + (x + w);				//calculate the index of the corresponding pixel in the result image
						size_t n0 = idx_offset(x, y, z);										//calculate the index for this pixel in the original image
						result.data()[n] = field<T>::_data[n0];									// copy the original image to the result image afer the border area
					}
				}
			}
			size_t l = w;
			size_t r = w + X() - 1;
			size_t t = w;
			size_t b = w + Y() - 1;
			size_t len = w;
			size_t len2 = w + Z() - 1;

			for (size_t y = l; y <= b; y++) for (size_t x = 0; x <= r; x++) for (size_t z = 0; z < w; z++) result(x, y, z) = result(x, t, len2);			    //pad the top
			for (size_t y = l; y <= b; y++) for (size_t x = l; x < r; x++) for (size_t z = len2 + 1; z < result.Z(); z++) result(x, y ,z) = result(x, b, len);	//pad the bottom					
			for (size_t y = 0; y < w; y++) for (size_t x = l; x < r; x++) for (size_t z = l; z <= len2; z++) result(x, y, z) = result(x, t, t);               //pad the left	
			for (size_t y = b + 1; y < result.Y(); y++) for (size_t x = l; x <= r; x++) for (size_t z = l; z <= len2; z++) result(x, y, z) = result(x, t, t);  //pad the right
			for (size_t y = 0; y < w; y++) for (size_t x = 0; x < r; x++) for (size_t z = 0; z < w; z++) result(x, y, z) = result(x, b, len);              //pad the top left	
			for (size_t y = b+1; y < result.Y(); y++) for (size_t x = 0; x < r; x++) for (size_t z = 0; z < w; z++) result(x, y, z) = result(x, b, len);    //pad the top right		
			for (size_t y = 0; y < w; y++) for (size_t x = l; x < r; x++) for (size_t z = len2 + 1; z < result.Z(); z++) result(x, y ,z) = result(x, b, len);  //pad the bottom left		
			for (size_t y = b + 1; y < result.Y(); y++) for (size_t x = l; x < result.X(); x++) for (size_t z = len2 + 1; z < result.Z(); z++) result(x, y, z) = result(x, b, len);//pad the bottom right



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
			if (pos == std::string::npos) {
				std::cout << "Cannot compose a volume from a single image." << std::endl;
				exit(1);
			}
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
		/// Load a volume from an NPY file. This function makes sure that the volume has four channels (even if there is only one color channel)
		/// </summary>
		/// <typeparam name="D"></typeparam>
		/// <param name="filename"></param>
		template<typename D = T>
		void load_npy(std::string filename) {
			field<T>::template load_npy<D>(filename);										// load the numpy file using the tira::field class
			if (field<T>::_shape.size() == 3)										// if the numpy array is only 2D, add a color channel of size 1
				field<T>::_shape.push_back(1);
		}

		T& operator()(size_t x, size_t y, size_t z, size_t c = 0) {
			return at(x, y, z, c);
		}

		/// <summary>
		/// Set all values in the image to a single constant
		/// </summary>
		/// <param name="v">Constant that all elements will be set to</param>
		field<T> operator=(T v) {														//set all elements of the image to a given value v
			return field<T>::operator=(v);
		}

		

		void resize(size_t x, size_t y, size_t z, size_t c = 0) {
			std::vector<size_t> shape = { z, y, x, c };
			tira::field<T>::resize(shape);
		}

		void resize(std::vector<size_t> shape) {
			field<T>::resize(shape);
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
		

		volume<T> dist() {
			volume<int> boundary;
			return _dist(boundary);
		}

		volume<T> sdf() {
			volume<int> boundary;
			volume<T> D = _dist(boundary);

			int width = D.X();
			int height = D.Y();
			int length = D.Z();

			std::vector<float>SDF;
			SDF.resize(height * width * length);

			for (int y = 0; y < height; y++)
			{
				for (int x = 0; x < width; x++)
				{
					for (int z = 0; z < length; z++)
					{
						SDF[((z * height) + y) * width + x] = D(x, y, z);
					}
				}

			}

			std::vector <float> frozenCells;
			frozenCells.resize(height * width * length);


			// initializing frozencells
			for (int y = 0; y < height; y++)
			{
				for (int x = 0; x < width; x++)
				{
					for (int z = 0; z < length; z++)
					{
						frozenCells[((z * height) + y) * width + x] = boundary(x, y, z);
					}
				}
			}

			// turn the whole input distance field to negative
			for (int i = 0; i < SDF.size(); i++) {
				SDF[i] = -1 * SDF[i];
			}


			double val; int gridPos;
			const int row = width;
			const int nx = width - 1;
			const int ny = height - 1;
			const int nz = length - 1;
			int ix = 0, iy = 0, iz = 0;

			// initialize a std::tuple to store the values of frozen cell
			std::stack<std::tuple<int, int, int>> stack = {};

			std::tuple<int, int, int> idsPair;

			// find the first unfrozen cell
			gridPos = 0;


			while (frozenCells[gridPos]) {
				ix += (ix < nx ? 1 : 0);
				iy += (iy < ny ? 1 : 0);
				iz += (iz < nz ? 1 : 0);
				gridPos = (iz * height + iy) * row + ix;
			}
			stack.push({ ix, iy, iz });
			// a simple pixel flood
			while (stack.size()) {
				idsPair = stack.top();
				stack.pop();
				ix = std::get<0>(idsPair);
				iy = std::get<1>(idsPair);
				iz = std::get<2>(idsPair);
				gridPos = (iz * height + iy) * row + ix;
				if (!frozenCells[gridPos]) {
					val = -1.0 * SDF[gridPos];
					SDF[gridPos] = val;
					frozenCells[gridPos] = true; // freeze cell when done
					if (ix > 0) {
						stack.push({ ix - 1, iy , iz });
					}
					if (ix < nx) {
						stack.push({ ix + 1, iy, iz });
					}
					if (iy > 0) {
						stack.push({ ix, iy - 1, iz });
					}
					if (iy < ny) {
						stack.push({ ix, iy + 1, iz });
					}
					if (iz > 0) {
						stack.push({ ix, iy, iz - 1 });
					}
					if (iz < nz) {
						stack.push({ ix, iy, iz + 1 });
					}
				}
			}


			for (int y = 0; y < height; y++)
			{
				for (int x = 0; x < width; x++)
				{
					for (int z = 0; z < length; z++)
					{
						D(x, y, z) = SDF[((z * height) + y) * width + x];
					}
				}

			}
			return D;
		}

		/// <summary>
		/// Generate a color map of the image.
		/// </summary>
		/// <param name="minval"></param>
		/// <param name="maxval"></param>
		/// <param name="colormap">type of colormap to apply</param>
		/// <returns></returns>
		volume<unsigned char> cmap(T minval, T maxval, const ColorMap colormap) {

			if (C() > 1) {																	// throw an exception if the image is more than one channel
				throw "Cannot create a color map from images with more than one channel!";
			}
			volume<unsigned char> color_result(X(), Y(), Z(), 3);									// create the new color image
			float a;																		// create a normalized value as a reference for the color map
			for (size_t i = 0; i < field<T>::_data.size(); i++) {
				cmap::cmap(field<T>::_data[i], minval, maxval, color_result.data()[i * 3 + 0], color_result.data()[i * 3 + 1], color_result.data()[i * 3 + 2], colormap);
			}
			return color_result;
		}

		/// <summary>
		/// Generate a color map of the image.
		/// </summary>
		/// <returns></returns>
		volume<unsigned char> cmap(ColorMap colormap = ColorMap::Brewer) {
			return cmap(minv(), maxv(), colormap);
		}

	};	// end class volume

	

}	// end namespace tira
