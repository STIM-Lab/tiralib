#pragma once

#include <tira/image.h>
#include <tira/progressbar.h>

#include <iomanip>
#include <filesystem>
#include <set>


namespace tira {
	/**
	 * @brief This class extends tira::field in the case where the field is 3D. This provides several convenient simplifications
	 * and addresses the fact that 3D volumes are usually indexed in (x, y, z) order.
	 * @tparam T
	 */
	template <class T>
	class volume : public field<T> {		// the image class extends field

	protected:

		std::vector<double> _spacing;

		/**
		 * @brief Private initialization function sets the shape of the volume and allocates memory to store the data points
		 * @param x is the size of the volume along the fastest dimension
		 * @param y is the size of the volume along the middle dimension
		 * @param z is the size of the volume along the slowest dimension
		 * @param c is the number of color channels, which is stored interleaved (colors are a really fast dimension)
		 */
		void m_Init(size_t x, size_t y, size_t z, size_t c = 1) {
			field<T>::m_shape.push_back(z);
			field<T>::m_shape.push_back(y);
			field<T>::m_shape.push_back(x);
			field<T>::m_shape.push_back(c);

			field<T>::m_Allocate();
		}

		/**
		 * @brief Evaluates the value at a specified z position for a grid of boxes. This function is used
		 *			to generate a basic "dummy" volume for testing.
		 * @param xi is the x index of the volume to test
		 * @param yi is the y index of the volume to test
		 * @param zi is the z index of the volume to test
		 * @param boxes is the number of boxes along a dimension
		 * @return is a binary value specifying whether or not a box exists at the current location
		 */
		static bool m_Grid(const float xi, const float yi, const float zi, const unsigned int boxes) {
			bool x_box, y_box, z_box;
			const float box_size = 1.0f / (float)boxes;
			if ((unsigned int)(zi / box_size) % 2)
				z_box = false;
			else
				z_box = true;

			if ((unsigned int)(yi / box_size) % 2)
				y_box = false;
			else
				y_box = true;

			if ((unsigned int)(xi / box_size) % 2)
				x_box = false;
			else
				x_box = true;

			return (x_box ^ y_box) ^ z_box;
		}

		/**
		 * @brief Retrieve a copy of the value at (x, y, z, c)
		 * @param x is the index of the value in the x dimension
		 * @param y is the index of the value in the y dimension
		 * @param z is the index of the value in the z dimension
		 * @param c is the index of the color value
		 * @return a copy of the value at the specified index coordinates
		 */
		T m_Get(size_t x, size_t y, size_t z, size_t c = 0) {
			size_t i = Idx(x, y, z, c);
			return field<T>::m_data[i];
		}

		volume<float> _dist(volume<int>& binary_boundary) {
			// create registers to hold the size
			int w = X();
			int h = Y();
			int l = Z();

			// resize the binary boundary grid
			binary_boundary.resize(field<T>::m_shape);

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

							if (m_Get(x, y, z) * m_Get(nx, ny, nz) <= 0) {				// if the product of the current cell and neighboring cell is negative (it is a boundary cell)
								binary_boundary(x, y, z) = 1;					// this cell is a boundary cell
								binary_boundary(nx, ny, nz) = 1;				// the neighboring cell is ALSO a boundary cell
							}
						}
					}
				}
			}


			volume<float> dist(w, h, l);										// create an image to store the distance field
			const float bignum = 9999.0f;
			dist = bignum;																// initialize the distance field to a very large value
			/*for (int y = 0; y < h; y++) {
				for (int x = 0; x < w; x++) {
					for (int z = 0; z < l; z++) {
						if (binary_boundary(x, y, z) == 1) {
							dist(x, y, z) = 0.0f;  // Ensure contour is strictly zero
						}
						else {
							dist(x, y, z) = bignum;  // Initialize all others to large value
						}
					}
				}
			}*/



			// calculate the distance for all boundary cells to the contour
			for (int y = 1; y < h - 1; y++) {											// for every row in the image
				for (int x = 1; x < w - 1; x++) {										// for every column in the image
					for (int z = 1; z < l - 1; z++) {										// for every length in the image
						if (binary_boundary(x, y, z)) {											// if the pixel (x, y,z) is in the boundary
							for (int k = 0; k < neighbors.size(); k++) {						// for every neighbor

								int nx = x + get<0>(neighbors[k]);								// calculate the x coordinate of the neighbor cell
								int ny = y + get<1>(neighbors[k]);								// calculate the y coordinate of the neighbor cell
								int nz = z + get<2>(neighbors[k]);								// calculate the z coordinate of the neighbor cell
								if (binary_boundary(nx, ny, nz)) {							// if the neighboring cell (nx, ny) is ALSO in the boundary
									float p_dist = abs(m_Get(x, y, z));
									float n_dist = abs(m_Get(nx, ny, nz));

									float da, db;
									if (p_dist == n_dist) {											// handle division by zero
										da = bignum;												// if the denominator is zero, we initialize the distance to its maximum
										db = bignum;
									}
									else {
										da = (abs(m_Get(x, y, z))) / (abs(m_Get(nx, ny, nz) - m_Get(x, y, z)));
										db = (abs(m_Get(nx, ny, nz))) / (abs(m_Get(nx, ny, nz) - m_Get(x, y, z)));
									}
									dist(x, y, z) = std::min(dist(x, y, z), da);									// minimum between distance and large boundary value of pixel (x,y)
									dist(nx, ny, nz) = std::min(dist(nx, ny, nz), db);					// minimum between distance and large boundary value of neighbor (nx, ny)
								}
							}
						}
					}
				}
			}


			// at this point, the distance image has the correct distances in cells surround the boundary
			

			FastSweep(dist);


			return dist;
		}

	public:
		/**
		 * @brief Default constructor creates a volume with isotropic spacing between voxels
		 */
		volume() : field<T>() {
			_spacing = { 1.0f, 1.0f, 1.0f };
		}			//initialize all variables, don't allocate any memory

		/**
		 * @brief Constructor builds a volume with the specified size, number of color channels, and spacing between voxels.
		 * @param nx is the number of pixels along the x dimension (fastest)
		 * @param ny is the number of pixels along the y dimension
		 * @param nz is the number of pixels along the z dimension
		 * @param nc is the number of color channels
		 * @param spacing is the spacing between voxels (default is 1.0 along all dimensions)
		 */
		volume(size_t nx, size_t ny, size_t nz, size_t nc = 1, std::vector<double> spacing = {1.0, 1.0, 1.0}) : volume() {
			m_Init(nx, ny, nz, nc);
			_spacing = spacing;
		}

		/**
		 * @brief Constructor creates a volume from a provided array
		 * @param data is a pointer to the data that will be copied into the volume
		 * @param nx is the number of voxels along the x-dimension (fastest)
		 * @param ny is the number of voxels along the y-dimension
		 * @param nz is the number of voxels along the z-dimension (slowest)
		 * @param nc is the number of interleaved color channels
		 * @param spacing is the spatial distance between samples (default is 1.0 along all dimensions)
		 */
		volume(T* data, size_t nx, size_t ny, size_t nz, size_t nc = 1, std::vector<double> spacing = {1.0, 1.0, 1.0}) : volume(nx, ny, nz, nc) {
			memcpy(&field<T>::m_data[0], data, field<T>::Bytes());									// use memcpy to copy the raw data into the field array
			_spacing = spacing;
		}

		/**
		 * @brief Generate a new volume from a tira::field shape. I don't really like this function because it doesn't obfuscate the main
		 * difference between a tira::volume and a tira::field, which is the order of dimensions. However it's a good way to create a volume
		 * that is the same size as a field.
		 * @param shape is a vector specifying the shape of the volume from slowest-to-fastest order: {z, y, x, c}
		 * @param spacing is the spacing between each voxel in (x, y, z) order (which is why I don't like this constructor)
		 */
		volume(std::vector<size_t> shape, std::vector<double> spacing = {1.0, 1.0, 1.0}) : field<T>(shape) {
			_spacing = spacing;
		}

		/**
		 * @brief Copy constructor that creates a volume from a given volume
		 * @param V is the volume to copy
		 */
		volume(const volume& V) : volume(V.X(), V.Y(), V.Z(), V.C()) {
			memcpy(&field<T>::m_data[0], &V.m_data[0], field<T>::Bytes());
			_spacing = V._spacing;
		}


		/**
		 * @brief Creates a new volume from a stack of images or a 3D image file (NPY)
		 * @param filemask is a mask used to glob a list of files to construct the volume
		 */
		volume(std::string filemask) {
			if (filemask.ends_with(".npy")) {
				load_npy(filemask);
			}
			else {
				std::vector<std::string> file_list(1);
				file_list[0] = filemask;
				Load(file_list);
			}
			_spacing = { 1.0f, 1.0f, 1.0f };
		}

		volume(std::vector<std::string> file_list) {
			Load(file_list);
			_spacing = { 1.0f, 1.0f, 1.0f };
		}

		/**
		 * @brief Is the default constructor (all of the member variables are smart so this can be simple)
		 */
		~volume() { }

		/**
		 * @brief Assignment operator fills the current field with the data from another field
		 * @param other is the right-hand-side operator in the assignment and the field containing the data to be copied
		 * @return a reference to the current field after the data has been copied from other
		 */
		volume& operator=(const volume& other) {
			field<T>::operator=(other);
			_spacing = other._spacing;
			return *this;
		}

		/**
		 * @brief Get the number of samples along the z-axis (slowest)
		 * @return the number of samples along the z-axis
		 */
		size_t Z() const { return field<T>::m_shape[0]; }

		/**
		 * @brief Get the number of samples along the y-axis
		 * @return the number of samples along the y-axis
		 */
		size_t Y() const { return field<T>::m_shape[1]; }

		/**
		 * @brief Get the number of samples along the x-axis (fastest)
		 * @return the number of samples along the x-axis
		 */
		size_t X() const { return field<T>::m_shape[2]; }

		/**
		 * @brief Get the number of color channels in the volume
		 * @return the number of color channels
		 */
		size_t C() const { return field<T>::m_shape[3]; }

		/**
		 * @brief Get the spacing between samples along the x-axis
		 * @return the spacing between samples along the x-axis
		 */
		double dx() const { return _spacing[0]; }

		/**
		 * @brief Get the spacing between samples along the y-axis
		 * @return the spacing between samples along the y-axis
		 */
		double dy() const { return _spacing[1]; }

		/**
		 * @brief Get the spacing between samples along the z-axis
		 * @return the spacing between samples along the z-axis
		 */
		double dz() const { return _spacing[2]; }

		/**
		 * @brief Calculates the size of the entire volume along the x-axis based on the number of samples and spacing between them
		 * @return the size of the volume along the x axis
		 */
		double sx() const { return X() * _spacing[0]; }

		/**
		 * @brief Calculates the size of the entire volume along the y-axis based on the number of samples and spacing between them
		 * @return the size of the volume along the y axis
		 */
		double sy() const { return Y() * _spacing[1]; }

		/**
		 * @brief Calculates the size of the entire volume along the z-axis based on the number of samples and spacing between them
		 * @return the size of the volume along the z axis
		 */
		double sz() const { return Z() * _spacing[2]; }

		/**
		 * @brief Calculate a 1D index into the data array given a set of spatial coordinates and a color
		 * @param x is the x index into the volume
		 * @param y is the y index into the volume
		 * @param z is the z index into the volume
		 * @param c is the color channel
		 * @return a one-dimensional index providing the location of the specified value in m_data
		 */
		size_t Idx(size_t x, size_t y, size_t z, size_t c = 0) const {
			size_t i = z * C() * X() * Y() + y * C() * X() + x * C() + c;
			return i;		// z * C * X * Y + y * C * X + x * C + c
		}

		/**
		 * @brief Run the fast sweeping algorithm to expand a distance field through an entire volume
		 * @param contour_distance is the initial distance field (usually specified near the contour surface)
		 */
		void FastSweep(volume<float>& contour_distance) {

			// Get volume dimensions
			int x = X();  // Number of voxels in x-direction
			int y = Y();  // Number of voxels in y-direction
			int z = Z();  // Number of voxels in z-direction

			int row = x;  // Row stride in flattened array is width (x)

			// Allocate 1D array to store the 3D volume
			std::vector<float> dist1d(x * y * z);

			// Copy 3D data into 1D array for easier access during sweeping
			for (int j = 0; j < y; j++) {
				for (int i = 0; i < x; i++) {
					for (int k = 0; k < z; k++) {
						dist1d[k * (x * y) + j * x + i] = contour_distance(i, j, k);
					}
				}
			}

			const int NSweeps = 8;  // Number of sweeping passes

			// Sweep directions along each axis for 8 directional passes
			const int dirX[NSweeps][3] = {
				{0, x - 1, 1}, {x - 1, 0, -1}, {x - 1, 0, -1}, {0, x - 1, 1},
				{0, x - 1, 1}, {x - 1, 0, -1}, {x - 1, 0, -1}, {0, x - 1, 1}
			};
			const int dirY[NSweeps][3] = {
				{0, y - 1, 1}, {0, y - 1, 1}, {y - 1, 0, -1}, {y - 1, 0, -1},
				{0, y - 1, 1}, {0, y - 1, 1}, {y - 1, 0, -1}, {y - 1, 0, -1}
			};
			const int dirZ[NSweeps][3] = {
				{0, z - 1, 1}, {0, z - 1, 1}, {0, z - 1, 1}, {0, z - 1, 1},
				{z - 1, 0, -1}, {z - 1, 0, -1}, {z - 1, 0, -1}, {z - 1, 0, -1}
			};

			// Temporary variables
			double aa[3];         // Stores nearest neighbor distances in x, y, z
			double tmp;           // For swapping
			double eps = 1e-6;    // Small epsilon to avoid floating point precision issues
			double d_new;         // Newly computed distance

			// Spacing in each direction
			const float dx = _spacing[0];
			const float dy = _spacing[1];
			const float dz = _spacing[2];
			const double f = 1.0; // RHS of Eikonal equation ∇d = f

			// Perform 8 directional sweeps
			for (int s = 0; s < NSweeps; s++) {
				for (int iy = dirY[s][0]; dirY[s][2] * iy <= dirY[s][1]; iy += dirY[s][2]) {
					for (int ix = dirX[s][0]; dirX[s][2] * ix <= dirX[s][1]; ix += dirX[s][2]) {
						for (int iz = dirZ[s][0]; dirZ[s][2] * iz <= dirZ[s][1]; iz += dirZ[s][2]) {

							// Compute 1D index of current voxel
							int gridPos = ((iz * y + iy) * row + ix);

							// Neighbor in y-direction
							if (iy == 0 || iy == (y - 1)) { // If on the y-boundary (top or bottom edge)
								aa[1] = dist1d[(iz * y + (iy == 0 ? iy + 1 : iy - 1)) * row + ix]; // Use the inner neighbor along y
							}
							else {
								aa[1] = std::min( // Otherwise, use the min of two neighbors in the y-direction
									dist1d[(iz * y + (iy - 1)) * row + ix],  // Neighbor in -y
									dist1d[(iz * y + (iy + 1)) * row + ix]   // Neighbor in +y
								);
							}

							// Neighbor in x-direction
							if (ix == 0 || ix == (x - 1)) { // If on the x-boundary (left or right edge)
								aa[0] = dist1d[(iz * y + iy) * row + (ix == 0 ? ix + 1 : ix - 1)]; // Use the inner neighbor along x
							}
							else {
								aa[0] = std::min( // Otherwise, use the min of two neighbors in the x-direction
									dist1d[(iz * y + iy) * row + (ix - 1)],  // Neighbor in -x
									dist1d[(iz * y + iy) * row + (ix + 1)]   // Neighbor in +x
								);
							}

							// Neighbor in z-direction
							if (iz == 0 || iz == (z - 1)) { // If on the z-boundary (front or back face)
								aa[2] = dist1d[((iz == 0 ? iz + 1 : iz - 1) * y + iy) * row + ix]; // Use the inner neighbor along z
							}
							else {
								aa[2] = std::min( // Otherwise, use the min of two neighbors in the z-direction
									dist1d[((iz - 1) * y + iy) * row + ix],  // Neighbor in -z
									dist1d[((iz + 1) * y + iy) * row + ix]   // Neighbor in +z
								);
							}

							// Sort neighbors in ascending order (so aa[0] <= aa[1] <= aa[2])
							if (aa[0] > aa[1]) { tmp = aa[0]; aa[0] = aa[1]; aa[1] = tmp; } // Swap if needed
							if (aa[1] > aa[2]) { tmp = aa[1]; aa[1] = aa[2]; aa[2] = tmp; } // Swap if needed
							if (aa[0] > aa[1]) { tmp = aa[0]; aa[0] = aa[1]; aa[1] = tmp; } // Swap if needed

							// Conservative estimate using smallest spacing and smallest neighbor
							double min_spacing = std::min({ dx, dy, dz }); // Choose the smallest spacing
							double d_curr = aa[0] + f * min_spacing;       // Initial guess: smallest neighbor + step

							//double d_curr = aa[0] + f / std::sqrt(1.0 / (dx*dx) + 1.0 / (dy*dy) + 1.0 / (dz*dz));


							// If the guess is less than next neighbor (valid), accept it
							if (d_curr <= (aa[1] + eps)) {
								d_new = d_curr; // Use 1 neighbor solution
							}
							else {
								// Solve quadratic using first two neighbors (x and y)
								double a = 1.0 / (dx * dx) + 1.0 / (dy * dy); // Coefficient a
								double b = -2.0 * (aa[0] / (dx * dx) + aa[1] / (dy * dy)); // Coefficient b
								double c = (aa[0] * aa[0]) / (dx * dx) + (aa[1] * aa[1]) / (dy * dy) - f * f; // Coefficient c

								double D = sqrt(b * b - 4.0 * a * c); // Discriminant
								d_curr = ((-b + D) > (-b - D) ? (-b + D) : (-b - D)) / (2.0 * a); // Choose smaller positive root

								// If still valid (less than next neighbor)
								if (d_curr <= (aa[2] + eps)) {
									d_new = d_curr; // Use 2 neighbor solution
								}
								else {
									// Solve quadratic using all three neighbors (x, y, z)
									a = 1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz); // Update a
									b = -2.0 * (aa[0] / (dx * dx) + aa[1] / (dy * dy) + aa[2] / (dz * dz)); // Update b
									c = (aa[0] * aa[0]) / (dx * dx) + (aa[1] * aa[1]) / (dy * dy) + (aa[2] * aa[2]) / (dz * dz) - f * f; // Update c

									D = sqrt(b * b - 4.0 * a * c); // Recompute discriminant
									d_new = ((-b + D) > (-b - D) ? (-b + D) : (-b - D)) / (2.0 * a); // Choose smaller positive root
								}
							}


							// Update voxel if new distance is smaller
							//dist1d[gridPos] = std::min(dist1d[gridPos], static_cast<float>(d_new));
							float abs_old = std::abs(dist1d[gridPos]);
							float sign = (dist1d[gridPos] >= 0.0f) ? 1.0f : -1.0f;
							dist1d[gridPos] = sign * std::min(abs_old, static_cast<float>(d_new));



						}
					}
				}
			}

			// Copy result from 1D array back to 3D volume
			for (int j = 0; j < y; j++) {
				for (int i = 0; i < x; i++) {
					for (int k = 0; k < z; k++) {
						dist(i, j, k) = dist1d[k * (x * y) + j * x + i];
					}
				}
			}

		}

		inline double smax() const { return std::max(sx(), std::max(sy(), sz())); }

		/// <summary>
		/// Returns a single channel of the current volume as an independent one-channel volume
		/// </summary>
		/// <param name="c">channel to return</param>
		/// <returns></returns>
		volume<T> channel(const size_t c) const {
			volume<T> r(X(), Y(), Z());							// create a new single-channel image
			T* dp = r.Data();									// create a pointer to the raw result
			size_t S = X() * Y() * Z();							// number of voxels
			for (size_t i = 0; i < S; i++) {					// for each voxel
				dp[i] = field<T>::m_data[i * C() + c];			// copy the value corresponding to the desired channel
			}
			return r;
		}

		/// <summary>
		/// Cast data types
		/// </summary>
		/// <typeparam name="V"></typeparam>
		template<typename V>
		operator volume<V>() {
			volume<V> r(this->Shape());					//create a new image
			std::copy(this->m_data.begin(), this->m_data.begin() + this->m_data.size(), r.Data());		//copy and cast the data
			return r;									//return the new image
		}

		T minv() {
			T m = *std::min_element(field<T>::m_data.begin(), field<T>::m_data.end());
			return m;
		}

		T maxv() {
			T m = *std::max_element(field<T>::m_data.begin(), field<T>::m_data.end());
			return m;
		}

		void spacing(double x, double y, double z) {
			std::vector<double> spacing = { x, y, z };
			_spacing = spacing;
		}

		std::vector<double> spacing() {
			return _spacing;
		}

		void generate_grid(unsigned int X = 32, unsigned int Y = 32, unsigned int Z = 32, unsigned int boxes = 1) {
			m_Init(X, Y, Z, 1);

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

						if (m_Grid(x, y, z, boxes)) {
							field<T>::m_data[idx] = sqrt(x * x + y * y + z * z) / sqrt(3) * 255;
						}
						else {
							field<T>::m_data[idx] = 0;
						}
					}
				}
			}
		}

		void generate_rgb(unsigned int X = 32, unsigned int Y = 32, unsigned int Z = 32, unsigned int boxes = 1) {
			m_Init(X, Y, Z, 3);


			float pixel_size[] = { 1.0f / X, 1.0f / Y, 1.0f / Z };

			size_t idx_z, idx_y, idx;
			float x, y, z;
			for (size_t zi = 0; zi < Z; zi++) {
				idx_z = zi * X * Y * 3;
				z = zi * pixel_size[2];

				for (size_t yi = 0; yi < Y; yi++) {
					idx_y = idx_z + yi * X * 3;
					y = yi * pixel_size[1];

					for (size_t xi = 0; xi < X; xi++) {
						idx = idx_y + xi * 3;
						x = xi * pixel_size[0];

						if (m_Grid(x, y, z, boxes)) {
							field<T>::m_data[idx + 0] = x * 255;
							field<T>::m_data[idx + 1] = y * 255;
							field<T>::m_data[idx + 2] = z * 255;
						}
						else {
							field<T>::m_data[idx + 0] = 0;
							field<T>::m_data[idx + 1] = 0;
							field<T>::m_data[idx + 2] = 0;
						}
					}
				}
			}
		}

		// This function writes an OBJ file where for every nonzero voxel 
        // in the input volume a cube is drawn. Each cube is defined by 8 vertices 
        // and 6 faces.
		void convertToCubeObj(const std::string& obj_file_name) {
			// Open the OBJ file for writing.
			std::ofstream obj_file(obj_file_name);
			if (!obj_file.is_open()) {
				std::cerr << "ERROR: Unable to open OBJ file for writing.\n";
				return;
			}

			int vertex_counter = 1; // OBJ vertices are 1-indexed.

			// Loop over every voxel in the volume.
			for (int y = 0; y < this->Y(); ++y) {
				for (int x = 0; x < this->X(); ++x) {
					for (int z = 0; z < this->Z(); ++z) {
						// Check if the voxel is non-zero
						if (this->at(x, y, z) > 0) {
							// Voxel center coordinates and half-size
							double cx = static_cast<double>(x);
							double cy = static_cast<double>(y);
							double cz = static_cast<double>(z);
							double hs = 0.5;

							// Define the 8 vertices of the cube
							double v[8][3] = {
								{cx - hs, cy - hs, cz - hs}, {cx + hs, cy - hs, cz - hs},
								{cx + hs, cy + hs, cz - hs}, {cx - hs, cy + hs, cz - hs},
								{cx - hs, cy - hs, cz + hs}, {cx + hs, cy - hs, cz + hs},
								{cx + hs, cy + hs, cz + hs}, {cx - hs, cy + hs, cz + hs}
							};

							// Write the vertices to the OBJ file
							for (int i = 0; i < 8; ++i) {
								obj_file << "v " << v[i][0] << " " << v[i][1] << " " << v[i][2] << "\n";
							}

							int base = vertex_counter;

							// Define the cube faces (6 faces, 4 vertices each)
							int faces[6][4] = {
								{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
								{3, 2, 6, 7}, {0, 3, 7, 4}, {1, 2, 6, 5}
							};

							// Write the faces to the OBJ file
							for (const auto& face : faces) {
								obj_file << "f " << base + face[0] << " " << base + face[1]
									<< " " << base + face[2] << " " << base + face[3] << "\n";
							}

							// Increment the vertex counter for the next cube
							vertex_counter += 8;
						}
					}
				}
			}

			obj_file.close();
			std::cout << "Cube OBJ file written successfully: " << obj_file_name << std::endl;
		}


		/// <summary>
		/// Retrieve a slice of the volume as a 2D image
		/// </summary>
		/// <param name="i">Slice to return</param>
		/// <param name="axis">Axis specifying the slice orientation</param>
		/// <returns></returns>
		image<T> slice(size_t i, unsigned int axis = 2) {
			size_t rx, ry;
			size_t sx, sy, sz;

			// calculate the size of the output slice image
			if (axis == 0) {
				rx = Y();
				ry = Z();
			}
			else if (axis == 1) {
				rx = X();
				ry = Z();

				image<T> result(rx, ry, C());					// allocate the output slice image

				for (size_t zi = 0; zi < Z(); zi++) {
					for (size_t xi = 0; xi < X(); xi++) {
						for (size_t ci = 0; ci < C(); ci++) {
							result(xi, zi, ci) = m_Get(xi, i, zi, ci);
						}
					}
				}
				return result;
			}
			else if (axis == 2) {
				rx = X();
				ry = Y();
				image<T> result(rx, ry, C());					// allocate the output slice image

				for (size_t yi = 0; yi < Y(); yi++) {
					for (size_t xi = 0; xi < X(); xi++) {
						for (size_t ci = 0; ci < C(); ci++) {
							result(xi, yi, ci) = m_Get(xi, yi, i, ci);
						}
					}
				}
				return result;
			}
		}

		tira::volume<float> derivative(unsigned int axis, unsigned int d, unsigned int order, bool print_coefs = false){
			tira::field<float> F = this->Derivative(axis, d, order, print_coefs);			// derivative of the field
			tira::volume<float> D(F.Shape(), _spacing);										// create a new volume
			memcpy(D.Data(), F.Data(), D.Bytes());											// copy the field to the volume
			// field axes: {Z, Y, X, C} => {0, 1, 2, 3}
			// spacing: {dx, dy, dz} => {0, 1, 2}
			int spacing_idx = 2 - static_cast<int>(axis);									// axis 2->dx, 1->dy, 0->dz
			D = D * static_cast<float>(1 / D._spacing[spacing_idx]);						// scale the gradient by the axis spacing
			return D;
		}

		tira::volume<float> gradmag(unsigned int order) {
			tira::field<float> F = field<float>::gradmag(order);
			tira::volume<float> R(F.Data(), X(), Y(), Z());
			return R;
		}

		/// <summary>
		/// Calculates the sign() function at each point in the volume and returns the result.
		/// </summary>
		/// <returns></returns>
		tira::volume<float> sign() {
			tira::field<float> F = field<float>::Sign();
			tira::volume<float> R(F.Data(), X(), Y(), Z());
			return R;
		}

		
		tira::volume<T> Phi_to_sdf() {
			tira::volume<T> phi = *this;								// original field
			tira::volume<T> sdf(phi.X(), phi.Y(), phi.Z(), phi.C(), _spacing);				// output SDF

			float BIG = 99999;
			sdf = BIG;

			int X = phi.X();
			int Y = phi.Y();
			int Z = phi.Z();

			
			// generate the template to detect voxels that are adjacent to the boundary at phi = 0
			std::vector<std::tuple<int, int, int>> neighbors;
			neighbors.emplace_back(-1, -1, -1);
			neighbors.emplace_back(-1, -1, 0);
			neighbors.emplace_back(-1, -1, 1);
			neighbors.emplace_back(-1, 0, -1);
			neighbors.emplace_back(-1, 0, 0);
			neighbors.emplace_back(-1, 0, 1);
			neighbors.emplace_back(-1, 1, -1);
			neighbors.emplace_back(-1, 1, 0);
			neighbors.emplace_back(-1, 1, 1);

			neighbors.emplace_back(0, -1, -1);
			neighbors.emplace_back(0, -1, 0);
			neighbors.emplace_back(0, -1, 1);
			neighbors.emplace_back(0, 0, -1);
			neighbors.emplace_back(0, 0, 0);
			neighbors.emplace_back(0, 0, 1);
			neighbors.emplace_back(0, 1, -1);
			neighbors.emplace_back(0, 1, 0);
			neighbors.emplace_back(0, 1, 1);

			neighbors.emplace_back(1, -1, -1);
			neighbors.emplace_back(1, -1, 0);
			neighbors.emplace_back(1, -1, 1);
			neighbors.emplace_back(1, 0, -1);
			neighbors.emplace_back(1, 0, 0);
			neighbors.emplace_back(1, 0, 1);
			neighbors.emplace_back(1, 1, -1);
			neighbors.emplace_back(1, 1, 0);
			neighbors.emplace_back(1, 1, 1);

			unsigned int band_cells = 0;
			// identifying boundary cells where phi changes sign
			for (int yi = 0; yi < Y; yi++) {							// for every row in the volume
				for (int xi = 0; xi < X; xi++) {						// for every column in the volume
					for (int zi = 0; zi < Z; zi++) {					// for every depth slice in the volume
						float d_min = sdf(xi, yi, zi);			// get the current approximated distance from the SDF
						for (int k = 0; k < neighbors.size(); k++) {	// for every neighbor

							int nxi = xi + get<0>(neighbors[k]);		// compute x index of neighbor
							int nyi = yi + get<1>(neighbors[k]);		// compute y index of neighbor
							int nzi = zi + get<2>(neighbors[k]);		// compute z index of neighbor

							// test if (nxi, nyi, nzi) is inside the image
							if (nxi >= 0 && nyi >= 0 && nzi >= 0 &&
								nxi < X && nyi < Y && nzi < Z) {

								float phi_x = phi(xi, yi, zi);				// get the phi values corresponding to x and n
								float phi_n = phi(nxi, nyi, nzi);

								if (phi_x * phi_n < 0.0f) {							// if x and n are on opposite sides of the boundary
									float phi_space = std::abs(phi_x - phi_n);		// calculate the spacing between the current cell and its neighbor in phi
									float phi_frac = std::abs(phi_x) / phi_space;	// calculate the fraction of the distance between the current and neighbor points

									float x = (float)xi * _spacing[0];
									float y = (float)yi * _spacing[1];
									float z = (float)zi * _spacing[2];
									float nx = (float)nxi * _spacing[0];
									float ny = (float)nyi * _spacing[1];
									float nz = (float)nzi * _spacing[2];

									float voxel_space = std::sqrt((x - nx) * (x - nx) + (y - ny) * (y - ny) + (z - nz) * (z - nz));
									float d = phi_frac * voxel_space;			// calculate the distance to the surface
									if (d < std::abs(d_min)) {					// if the current point is closer to the surface than the previous guess
										d_min = d;
									}
									band_cells++;
								}
							}
						}
						sdf(xi, yi, zi) = d_min;
					}
				}
			}
			FastSweep(sdf);

			tira::volume<float> sign_phi = phi.sign();
			sdf = sdf * sign_phi;


			return sdf;

			//return dist;
		}

		/**
		 * Calculate the gradient of the volume using second-order finite differences for internal voxels and first-order
		 * finite differences for external voxels. This is the standard for most mathematical libraries, and may have something
		 * to do with preventing unusually large values at the edges.
		 * @param axis axis along which the gradient is calculated (X = 0, Y = 1, Z = 2)
		 * @return a volume containing the specified gradient
		 */
		volume gradient(const unsigned int axis) {
			volume output(X(), Y(), Z());

			double d = _spacing[axis];

			size_t i[3];
			size_t S[3] = {X(), Y(), Z()};
			for (i[2] = 0; i[2] < Z(); i[2]++) {
				for (i[1] = 0; i[1] < Y(); i[1]++) {
					for (i[0] = 0; i[0] < X(); i[0]++) {

						// if we are at the left edge of the axis dimension, use the forward difference
						if (i[axis] == 0) {
							std::vector p_this = {i[2], i[1], i[0]};
							std::vector p_right = {i[2], i[1], i[0]};
							p_right[2 - axis] += 1;

							T a = field<T>::Read(p_this);
							T b = field<T>::Read(p_right);

							output(i[0], i[1], i[2]) = (b - a) / d;
						}
						// if we are at the right edge of the axis dimension, use the backward difference
						else if (i[axis] == S[axis] - 1) {
							std::vector p_this = {i[2], i[1], i[0]};
							std::vector p_left = {i[2], i[1], i[0]};
							p_left[2 - axis] -= 1;

							T b = field<T>::Read(p_this);
							T a = field<T>::Read(p_left);
							output(i[0], i[1], i[2]) = (b - a) / d;
						}
						else {
							std::vector p_left = {i[2], i[1], i[0]};
							p_left[2 - axis] -= 1;

							std::vector p_right = {i[2], i[1], i[0]};
							p_right[2 - axis] += 1;

							T a = field<T>::Read(p_left);
							T b = field<T>::Read(p_right);

							output(i[0], i[1], i[2]) = (b - a) / (2.0 * d);
						}
					}
				}
			}
			return output;
		}
	
		// Calculate gradient along dx (3D) with spacing consideration
		tira::volume<float> gradient_dx()
		{
			tira::volume<float> output(X(), Y(), Z());

			double dx = _spacing[0]; // Spacing along X-direction

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
							output(j, i, k) = (m_Get(j_right, i, k) - m_Get(j_left, i, k)) / dx;
						}
						else if (j_right >= X())
						{
							j_right = X() - 1;
							j_left = j_right - 1;
							output(j, i, k) = (m_Get(j_right, i, k) - m_Get(j_left, i, k)) / dx;
						}
						else
						{
							output(j, i, k) = (m_Get(j_right, i, k) - m_Get(j_left, i, k)) / (2.0 * dx);
						}
					}
				}
			}

			return output;
		}


		// Calculate gradient along dy (3D) with spacing consideration
		tira::volume<float> gradient_dy()
		{
			tira::volume<float> output(X(), Y(), Z());

			double dy = _spacing[1]; // Spacing along Y-direction

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
							output(j, i, k) = (m_Get(j, i_right, k) - m_Get(j, i_left, k)) / dy;
						}
						else if (i_right >= Y())
						{
							i_right = Y() - 1;
							i_left = i_right - 1;
							output(j, i, k) = (m_Get(j, i_right, k) - m_Get(j, i_left, k)) / dy;
						}
						else
						{
							output(j, i, k) = (m_Get(j, i_right, k) - m_Get(j, i_left, k)) / (2.0 * dy);
						}
					}
				}
			}

			return output;
		}


		// Calculate gradient along dz (3D) with spacing consideration
		tira::volume<float> gradient_dz()
		{
			tira::volume<float> output(X(), Y(), Z());

			double dz = _spacing[2]; // Spacing along Z-direction

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
							output(j, i, k) = (m_Get(j, i, k_right) - m_Get(j, i, k_left)) / dz;
						}
						else if (k_right >= Z())
						{
							k_right = Z() - 1;
							k_left = k_right - 1;
							output(j, i, k) = (m_Get(j, i, k_right) - m_Get(j, i, k_left)) / dz;
						}
						else
						{
							output(j, i, k) = (m_Get(j, i, k_right) - m_Get(j, i, k_left)) / (2.0 * dz);
						}
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
		template <typename D>
		tira::volume<T>convolve3D(tira::volume<D> mask) {

			tira::volume<T>result(X() - (mask.X() - 1), Y() - (mask.Y() - 1), Z() - (mask.Z() - 1));		// output image will be smaller than the input (only valid region returned)

			for (size_t yi = 0; yi < result.Y(); yi++) {
				for (size_t xi = 0; xi < result.X(); xi++) {
					for (size_t zi = 0; zi < result.Z(); zi++) {
						T sum = static_cast<T>(0);
						for (size_t vi = 0; vi < mask.Y(); vi++) {
							for (size_t ui = 0; ui < mask.X(); ui++) {
								for (size_t wi = 0; wi < mask.Z(); wi++) {

									sum += m_Get((xi + ui), (yi + vi), (zi + wi)) * mask(ui, vi, wi);
									
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
							sum += m_Get(xi, (yi + ui), zi) * K[ui];
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

		// Convolves the volume with a Gaussian kernel and returns the result. This function
		//    adjusts the sigma value to account for the pixel size (sigma is given in terms of the pixel spacing)
		tira::volume<float> gaussian_filter(float sigma, int window_factor, bool pad = false) {

			// calculate the kernel sizes along each dimension (in voxel units)
			float sigma_x = sigma / _spacing[0];
			int kernel_size_x = window_factor * sigma_x;
			if (kernel_size_x % 2 == 0) kernel_size_x++;

			float sigma_y = sigma / _spacing[1];
			int kernel_size_y = window_factor * sigma_y;
			if (kernel_size_y % 2 == 0) kernel_size_y++;

			float sigma_z = sigma / _spacing[2];
			int kernel_size_z = window_factor * sigma_z;
			if (kernel_size_z % 2 == 0) kernel_size_z++;

			// If padding is requested, apply it and re-call the filter with pad=false
			if (pad) {
				int rx = kernel_size_x / 2;
				int ry = kernel_size_y / 2;
				int rz = kernel_size_z / 2;
				tira::volume<float> padded = this->border_separable(rx, ry, rz, 0.0f);
				return padded.gaussian_filter(sigma, window_factor, false);
			}

			// ----- X Axis -----
			float miu_x = static_cast<float>(kernel_size_x) / 2.0f;
			float* kernel_x = (float*)malloc(kernel_size_x * sizeof(float));
			for (int xi = 0; xi < kernel_size_x; xi++) {
				int u = 2 * sigma_x * sigma_x;
				kernel_x[xi] = 1.0f / sqrt(u * 3.14159265358979323846f) * exp(-(xi - miu_x) * (xi - miu_x) / u);
			}
			// normalize kernel_x
			float sum_x = 0.0f;
			for (int i = 0; i < kernel_size_x; i++) sum_x += kernel_x[i];
			for (int i = 0; i < kernel_size_x; i++) kernel_x[i] /= sum_x;

			tira::volume<float> result_x(X() - (kernel_size_x - 1), Y(), Z());
			for (size_t yi = 0; yi < result_x.Y(); yi++) {
				for (size_t xi = 0; xi < result_x.X(); xi++) {
					for (size_t zi = 0; zi < result_x.Z(); zi++) {
						float sum = 0;
						for (size_t ui = 0; ui < kernel_size_x; ui++) {
							sum += m_Get(xi + ui, yi, zi) * kernel_x[ui];
						}
						result_x(xi, yi, zi) = sum;
					}
				}
			}
			free(kernel_x);

			// ----- Y Axis -----
			float miu_y = static_cast<float>(kernel_size_y) / 2.0f;
			float* kernel_y = (float*)malloc(kernel_size_y * sizeof(float));
			for (int yi = 0; yi < kernel_size_y; yi++) {
				int u = 2 * sigma_y * sigma_y;
				kernel_y[yi] = 1.0f / sqrt(u * 3.14159265358979323846f) * exp(-(yi - miu_y) * (yi - miu_y) / u);
			}
			// normalize kernel_y
			float sum_y = 0.0f;
			for (int i = 0; i < kernel_size_y; i++) sum_y += kernel_y[i];
			for (int i = 0; i < kernel_size_y; i++) kernel_y[i] /= sum_y;

			tira::volume<float> result_y(result_x.X(), result_x.Y() - (kernel_size_y - 1), result_x.Z());
			for (size_t yi = 0; yi < result_y.Y(); yi++) {
				for (size_t xi = 0; xi < result_y.X(); xi++) {
					for (size_t zi = 0; zi < result_y.Z(); zi++) {
						float sum = 0;
						for (size_t ui = 0; ui < kernel_size_y; ui++) {
							sum += result_x(xi, yi + ui, zi) * kernel_y[ui];
						}
						result_y(xi, yi, zi) = sum;
					}
				}
			}
			free(kernel_y);

			// ----- Z Axis -----
			float miu_z = static_cast<float>(kernel_size_z) / 2.0f;
			float* kernel_z = (float*)malloc(kernel_size_z * sizeof(float));
			for (int zi = 0; zi < kernel_size_z; zi++) {
				int u = 2 * sigma_z * sigma_z;
				kernel_z[zi] = 1.0f / sqrt(u * 3.14159265358979323846f) * exp(-(zi - miu_z) * (zi - miu_z) / u);
			}
			// normalize kernel_z
			float sum_z = 0.0f;
			for (int i = 0; i < kernel_size_z; i++) sum_z += kernel_z[i];
			for (int i = 0; i < kernel_size_z; i++) kernel_z[i] /= sum_z;

			tira::volume<float> result_z(result_y.X(), result_y.Y(), result_y.Z() - (kernel_size_z - 1));
			for (size_t yi = 0; yi < result_z.Y(); yi++) {
				for (size_t xi = 0; xi < result_z.X(); xi++) {
					for (size_t zi = 0; zi < result_z.Z(); zi++) {
						float sum = 0;
						for (size_t ui = 0; ui < kernel_size_z; ui++) {
							sum += result_y(xi, yi, zi + ui) * kernel_z[ui];
						}
						result_z(xi, yi, zi) = sum;
					}
				}
			}
			free(kernel_z);

			return result_z;
		}


		tira::volume<float> reinitialize_osher(float d_tau, int num_iters) {

			tira::volume<float> phi = *this;

			const int nx = phi.X();
			const int ny = phi.Y();
			const int nz = phi.Z();

			const float dx = _spacing[0];
			const float dy = _spacing[1];
			const float dz = _spacing[2];

			// Initialize psi as a copy of φ to evolve it
			tira::volume<float> psi = phi;

			// Iterate over  time t
			for (int iter = 0; iter < num_iters; ++iter) {
				for (int z = 1; z < nz - 1; ++z) {
					for (int y = 1; y < ny - 1; ++y) {
						for (int x = 1; x < nx - 1; ++x) {

							// Compute sign function of φ 
							float val = phi(x, y, z);
							float sgn = val / std::sqrt(val * val + 1e-8f);  // sgn_ε(φ)

							// Compute forward and backward differences for psi
							float Dx_plus = (psi(x + 1, y, z) - psi(x, y, z)) / dx;
							float Dx_minus = (psi(x, y, z) - psi(x - 1, y, z)) / dx;
							float Dy_plus = (psi(x, y + 1, z) - psi(x, y, z)) / dy;
							float Dy_minus = (psi(x, y, z) - psi(x, y - 1, z)) / dy;
							float Dz_plus = (psi(x, y, z + 1) - psi(x, y, z)) / dz;
							float Dz_minus = (psi(x, y, z) - psi(x, y, z - 1)) / dz;

							// Godunov's scheme: compute |∇psi| using one-sided upwind differences
							float Dx_pos, Dx_neg;
							if (Dx_minus > 0.0f) Dx_pos = Dx_minus * Dx_minus;
							else Dx_pos = 0.0f;
							if (Dx_plus < 0.0f) Dx_pos += Dx_plus * Dx_plus;

							if (Dx_minus < 0.0f) Dx_neg = Dx_minus * Dx_minus;
							else Dx_neg = 0.0f;
							if (Dx_plus > 0.0f) Dx_neg += Dx_plus * Dx_plus;

							float Dy_pos, Dy_neg;
							if (Dy_minus > 0.0f) Dy_pos = Dy_minus * Dy_minus;
							else Dy_pos = 0.0f;
							if (Dy_plus < 0.0f) Dy_pos += Dy_plus * Dy_plus;

							if (Dy_minus < 0.0f) Dy_neg = Dy_minus * Dy_minus;
							else Dy_neg = 0.0f;
							if (Dy_plus > 0.0f) Dy_neg += Dy_plus * Dy_plus;

							float Dz_pos, Dz_neg;
							if (Dz_minus > 0.0f) Dz_pos = Dz_minus * Dz_minus;
							else Dz_pos = 0.0f;
							if (Dz_plus < 0.0f) Dz_pos += Dz_plus * Dz_plus;

							if (Dz_minus < 0.0f) Dz_neg = Dz_minus * Dz_minus;
							else Dz_neg = 0.0f;
							if (Dz_plus > 0.0f) Dz_neg += Dz_plus * Dz_plus;

							float grad_psi = 0.0f;

							// Select correct upwind norm based on sign of φ
							if (sgn > 0.0f) {
								grad_psi = std::sqrt(Dx_pos + Dy_pos + Dz_pos);
							}
							else {
								grad_psi = std::sqrt(Dx_neg + Dy_neg + Dz_neg);
							}

							// Compute update term: sgn(φ) * (|∇psi| - 1)
							float H = grad_psi - 1.0f;

							// Euler update step
							psi(x, y, z) = psi(x, y, z) - d_tau * sgn * H;
						}
					}
				}
			}

			// Return the evolved (signed distance approximation)
			return psi;
		}

		/// <summary>
        /// Upsamples the current volume by the given factor using trilinear interpolation.
        /// Each dimension is scaled by 'factor', and new voxel values are interpolated
        /// <param name="factor">Upsampling factor (e.g., 2 for doubling resolution)</param>
        /// <returns>A new volume with higher resolution and interpolated values</returns>


		tira::volume<float> Trilinear_interpolation(int factor) {
			int new_x = X() * factor;
			int new_y = Y() * factor;
			int new_z = Z() * factor;

			tira::volume<float> upsampled(new_x, new_y, new_z);

			for (int z = 0; z < new_z; z++) {
				float zf = (float)z / factor;
				int z0 = (int)zf;
				int z1 = z0 + 1;
				if (z1 >= Z()) z1 = z0;

				float wz = zf - z0;

				for (int y = 0; y < new_y; y++) {
					float yf = (float)y / factor;
					int y0 = (int)yf;
					int y1 = y0 + 1;
					if (y1 >= Y()) y1 = y0;

					float wy = yf - y0;

					for (int x = 0; x < new_x; x++) {
						float xf = (float)x / factor;
						int x0 = (int)xf;
						int x1 = x0 + 1;
						if (x1 >= X()) x1 = x0;

						float wx = xf - x0;

						// Trilinear interpolation
						float c000 = m_Get(x0, y0, z0);
						float c100 = m_Get(x1, y0, z0);
						float c010 = m_Get(x0, y1, z0);
						float c110 = m_Get(x1, y1, z0);
						float c001 = m_Get(x0, y0, z1);
						float c101 = m_Get(x1, y0, z1);
						float c011 = m_Get(x0, y1, z1);
						float c111 = m_Get(x1, y1, z1);

						float c00 = c000 * (1 - wx) + c100 * wx;
						float c10 = c010 * (1 - wx) + c110 * wx;
						float c01 = c001 * (1 - wx) + c101 * wx;
						float c11 = c011 * (1 - wx) + c111 * wx;

						float c0 = c00 * (1 - wy) + c10 * wy;
						float c1 = c01 * (1 - wy) + c11 * wy;

						float value = c0 * (1 - wz) + c1 * wz;

						upsampled(x, y, z) = value;
					}
				}
			}

			return upsampled;
		}


		//Cubic interpolation using Catmull–Rom spline
		tira::volume<float> Tricubic_interpolation(int factor) {
			int new_x = X() * factor;                      // Scale original X dimension by factor
			int new_y = Y() * factor;                      // Scale original Y dimension by factor
			int new_z = Z() * factor;                      // Scale original Z dimension by factor

			tira::volume<float> upsampled(new_x, new_y, new_z); // Create the upsampled volume with new dimensions

			// Cubic interpolation using Catmull–Rom spline
			auto cubicInterpolate = [&](const float p[4], float t) {
				float a = -0.5f * p[0] + 1.5f * p[1] - 1.5f * p[2] + 0.5f * p[3]; // Compute coefficient a
				float b = p[0] - 2.5f * p[1] + 2.0f * p[2] - 0.5f * p[3];           // Compute coefficient b
				float c = -0.5f * p[0] + 0.5f * p[2];                                // Compute coefficient c
				float d = p[1];                                                      // Coefficient d is the central value
				return ((a * t + b) * t + c) * t + d;                                // Evaluate cubic polynomial at t
				};

			for (int z = 0; z < new_z; z++) {                   // Loop through each z coordinate in the upsampled volume
				float zf = (float)z / factor;                 // Convert integer z to float 
				int z_int = (int)zf;                          // Integer part of zf (floor)
				float dz = zf - z_int;                        // Fractional part of z for interpolation

				for (int y = 0; y < new_y; y++) {               // Loop through each y coordinate in the upsampled volume
					float yf = (float)y / factor;             // Convert integer y to float 
					int y_int = (int)yf;                      // Integer part of yf (floor)
					float dy = yf - y_int;                    // Fractional part of y for interpolation

					for (int x = 0; x < new_x; x++) {           // Loop through each x coordinate in the upsampled volume
						float xf = (float)x / factor;         // Convert integer x to float
						int x_int = (int)xf;                  // Integer part of xf (floor)
						float dx = xf - x_int;                // Fractional part of x for interpolation

						float interpz[4]; //  hold intermediate results from y-interpolation along z

						// Interpolate over the 4 z-neighbors (from z_int - 1 to z_int + 2)
						for (int iz = 0; iz < 4; iz++) {
							int z_idx = z_int - 1 + iz;       // Compute original z index for neighbor
							// ensure it lies within [0, Z()-1]
							if (z_idx < 0) {
								z_idx = 0;
							}
							else if (z_idx >= Z()) {
								z_idx = Z() - 1;
							}

							float interpy[4]; // hold intermediate results from x-interpolation along y

							// Interpolate over the 4 y-neighbors (from y_int - 1 to y_int + 2)
							for (int iy = 0; iy < 4; iy++) {
								int y_idx = y_int - 1 + iy;   // Compute original y index for neighbor
								// ensure it lies within [0, Y()-1]
								if (y_idx < 0) {
									y_idx = 0;
								}
								else if (y_idx >= Y()) {
									y_idx = Y() - 1;
								}
								float p[4]; // store values from 4 x-neighbors

								// Interpolate over the 4 x-neighbors (from x_int - 1 to x_int + 2)
								for (int ix = 0; ix < 4; ix++) {
									int x_idx = x_int - 1 + ix; // Compute original x index for neighbor
									// ensure it lies within [0, X()-1]
									if (x_idx < 0) {
										x_idx = 0;
									}
									else if (x_idx >= X()) {
										x_idx = X() - 1;
									}
									p[ix] = m_Get(x_idx, y_idx, z_idx); // Retrieve voxel value from original volume
								}
								interpy[iy] = cubicInterpolate(p, dx); // Interpolate along x-axis for fixed y and z
							}
							interpz[iz] = cubicInterpolate(interpy, dy); // Interpolate along y-axis for fixed z
						}
						float value = cubicInterpolate(interpz, dz);  // Final interpolation along z-axis
						upsampled(x, y, z) = value;                    // Assign the computed value to the upsampled volume
					}
				}
			}
			upsampled._spacing[0] /= factor;
			upsampled._spacing[1] /= factor;
			upsampled._spacing[2] /= factor;
			return upsampled; // Return the upsampled volume
		}



		/// <summary>
		/// Add a border to an image
		/// </summary>
		/// <param name="w">Width of the border in pixels</param>
		/// <param name="value">Value used to generate the border</param>
		/// <returns></returns>
		tira::volume<T>border(size_t w, T value) {
			tira::volume<float> result(X() + w * 2, Y() + w * 2, Z() + w * 2);
			result = value;													// assign the border value to all pixels in the new volume
 														
			for (size_t y = 0; y < Y(); y++) {								//for each pixel in the original image
				for (size_t x = 0; x < X(); x++) {
					for (size_t z = 0; z < Z(); z++) {
						size_t n = ((z + w) * (X() + w * 2) * (Y() + w * 2)) + ((y + w) * (X() + w * 2)) + (x + w);				//calculate the index of the corresponding pixel in the result image
						size_t n0 = Idx(x, y, z);										//calculate the index for this pixel in the original image
						result.Data()[n] = field<T>::m_data[n0];									// copy the original image to the result image afer the border area
					}
				}
			}
			return result;
		}

		/// Add a border around the volume, specifying the width for each dimension separately
		//tira::volume<T> border(size_t wx, size_t wy, size_t wz, T value) {

		//}

		
		tira::volume<T> border_separable(size_t wx, size_t wy, size_t wz, T value) {
			tira::volume<T> result(X() + 2 * wx, Y() + 2 * wy, Z() + 2 * wz);
			result = value; // Fill entire volume with padding value

			for (size_t y = 0; y < Y(); ++y) {
				for (size_t x = 0; x < X(); ++x) {
					for (size_t z = 0; z < Z(); ++z) {
						result(x + wx, y + wy, z + wz) = (*this)(x, y, z);
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
		tira::volume<T>border_replicate_3D(size_t w) {

			tira::volume<float> result(X() + w * 2, Y() + w * 2, Z() + w * 2);

			result = 0;
			//result = value;														//assign the border value to all pixels in the new image
			for (size_t y = 0; y < Y(); y++) {								//for each pixel in the original image
				for (size_t x = 0; x < X(); x++) {
					for (size_t z = 0; z < Z(); z++) {
						size_t n = ((z + w) * (X() + w * 2) * (Z() + w * 2)) + ((y + w) * (Z() + w * 2)) + (x + w);				//calculate the index of the corresponding pixel in the result image
						size_t n0 = Idx(x, y, z);										//calculate the index for this pixel in the original image
						result.Data()[n] = field<T>::m_data[n0];									// copy the original image to the result image afer the border area
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

		tira::volume<T>border_remove(size_t w) {
			return crop(w, w, w, X() - 2 * w, Y() - 2 * w, Z() - 2 * w);
		}

		/// <summary>
		/// Crops a volume to a desired size
		/// </summary>
		/// <param name="x0">X boundary of the cropped image</param>
		/// <param name="y0">Y boundary of the cropped image</param>
		/// <param name="z0">Z boundary of the cropped image</param>
		/// <param name="w">Width of the cropped image</param>
		/// <param name="h">Height of the cropped image</param>
		/// <param name="d">Depth of the cropped image</param>
		/// <returns></returns>
		tira::volume<T> crop(size_t x0, size_t y0, size_t z0, size_t w, size_t h, size_t d) {
			if (x0 + w > X() || y0 + h > Y() || z0 + d > Z()) {
				std::cout << "ERROR: cropped volume contains an invalid region." << std::endl;
				exit(1);
			}
			tira::volume<T> result(w, h, d, C());												// create the output cropped image

			const size_t line_bytes = w * C() * sizeof(T);										// calculate the number of bytes in a slice
			for (size_t zi = 0; zi < d; zi++) {													// loop through each depth slice
				for (size_t yi = 0; yi < h; yi++) {												// loop through each row
					size_t srci = ((z0 + zi) * Y() * X() + (y0 + yi) * X() + x0) * C();			// calculate the source index
					size_t dsti = (zi * h * w  + yi * w) * C();									// calculate the destination index
					memcpy(&result.m_data[dsti], &field<T>::m_data[srci], line_bytes);			// copy the data
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

		/**
		 * @brief Load a list of files and assemble them as a volume
		 * @param file_names is a vector containing all of the file names to load
		 */
		void Load(std::vector<std::string> file_names, const bool progress = true) {

			tira::image<T> test(file_names[0]);

			m_Init(test.width(), test.height(), file_names.size(), test.channels());			// allocate space for the new shape					
	
			// copies data from each image file to the destination data
			for (size_t zi = 0; zi < Z(); zi++) {											// iterate through each file

				tira::image<T> img(file_names[zi]);											// load the image file

				for (size_t yi = 0; yi < Y(); yi++) {
					for (size_t xi = 0; xi < X(); xi++) {
						for (size_t c = 0; c < C(); c++) {
							size_t volume_index = Idx(xi, yi, zi, c);
							field<T>::m_data[volume_index] = img(xi, yi, c);
						}
					}
				}
				if (progress)															// write a progress bar to the console
					progressbar((float)zi / (float)Z());
			}
			std::cout << std::endl;
		}

		void LoadStack(std::string directory) {
			if (field<T>::m_shape.size() != 0) {							// if a volume is already allocated, clear it
				field<T>::m_shape.clear();
				field<T>::m_data.clear();
			}
			// generate a list of file names from the mask

			// ------------- using filesystem ---------------------- //
			
			// obtaining the path string from input
			/*std::size_t pos = file_mask.rfind("*");												// position of * in file_mask
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
			*/
			std::set<std::filesystem::path> sorted_files;
			for (auto& p : std::filesystem::directory_iterator(directory))							// iterate through each file and stores the ones with
			{																					// desired extension in file_name vector
				//if (p.path().extension() == ext) {
				sorted_files.insert(directory + "/" + p.path().filename().string());
				//}
			}
			std::vector<std::string> file_names;
			for(auto &filename : sorted_files) {
				file_names.push_back(filename.string());
			}
			if(file_names.size() != 0)
				Load(file_names);
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
			field<T>::template LoadNpy<D>(filename);										// load the numpy file using the tira::field class
			if (field<T>::m_shape.size() == 3)										// if the numpy array is only 2D, add a color channel of size 1
				field<T>::m_shape.push_back(1);
		}

		/// Save the volume as a Numpy file. If there is only one channel, this function will squeeze the dimensions.
		template<typename D = T>
		void save_npy(std::string filename) {
			std::vector<size_t> squeezed_shape = field<T>::m_shape;
			if (squeezed_shape.size() == 4 && field<T>::m_shape[3] == 1) {
				squeezed_shape.pop_back();
			}
			field<T>::template SaveNpy<D>(filename, squeezed_shape);
		}

		T& operator()(size_t x, size_t y, size_t z, size_t c = 0) {
			return field<T>::m_data[Idx(x, y, z, c)];
		}

		/// <summary>
		/// Set all values in the image to a single constant
		/// </summary>
		/// <param name="v">Constant that all elements will be set to</param>
		field<T> operator=(T v) {														//set all elements of the image to a given value v
			return field<T>::operator=(v);
		}

		tira::volume<T> operator*(T rhs) {
			size_t N = field<T>::Size();									// calculate the total number of values in the volume
			tira::volume<T> r(this->Shape(), this->_spacing);								// allocate space for the resulting image
			for (size_t n = 0; n < N; n++)
				r.m_data[n] = field<T>::m_data[n] * rhs;						// multiply the individual samples
			return r;														// return the final result
		}

		tira::volume<T> operator+(T rhs) {
			size_t N = field<T>::Size();									// calculate the total number of values in the volume
			tira::volume<T> r(this->Shape(), this->_spacing);								// allocate space for the resulting image
			for (size_t n = 0; n < N; n++)
				r.m_data[n] = field<T>::m_data[n] + rhs;						// add the individual pixels
			return r;														// return the summed result
		}

		tira::volume<T> operator-(T rhs) {
			size_t N = field<T>::Size();									// calculate the total number of values in the volume
			tira::volume<T> r(this->Shape(), this->_spacing);								// allocate space for the resulting image
			for (size_t n = 0; n < N; n++)
				r.m_data[n] = field<T>::m_data[n] - rhs;						// add the individual pixels
			return r;														// return the summed result
		}

		volume<T> operator*(volume<T> rhs) {							// point-wise multiplication
			if (X() != rhs.X() || Y() != rhs.Y())
				throw std::runtime_error("Images dimensions are incompatible");

			if (C() == rhs.C()) {					// if both images have the same number of color channels
				tira::volume<T> result(this->Shape(), this->_spacing);
				for (size_t i = 0; i < field<T>::Size(); i++) {
					result.m_data[i] = field<T>::m_data[i] * rhs.m_data[i];
				}
				return result;
			}
			else if (C() == 1) {
				tira::volume<T> result(X(), Y(), Z(), rhs.C(), this->_spacing);
				for (size_t zi = 0; zi < Z(); zi++) {
					for (size_t yi = 0; yi < Y(); yi++) {
						for (size_t xi = 0; xi < X(); xi++) {
							for (size_t ci = 0; ci < rhs.C(); ci++) {
								result(xi, yi, zi, ci) = m_Get(xi, yi, zi) * rhs(xi, yi, zi, ci);
							}
						}
					}
				}
				return result;
			}
			else if (rhs.C() == 1) {
				tira::volume<T> result(this->Shape(), this->_spacing);
				for (size_t zi = 0; zi < Z(); zi++) {
					for (size_t yi = 0; yi < Y(); yi++) {
						for (size_t xi = 0; xi < X(); xi++) {
							for (size_t ci = 0; ci < C(); ci++) {
								result(xi, yi, zi, ci) = m_Get(xi, yi, zi, ci) * rhs(xi, yi, zi);
							}
						}
					}
				}
				return result;
			}
			else {
				throw std::runtime_error("Number of color channels are incompatible");
			}
		}

		volume<T> operator+(volume<T> rhs) {							// point-wise multiplication
			if (X() != rhs.X() || Y() != rhs.Y())
				throw std::runtime_error("Images dimensions are incompatible");

			tira::volume<T> result(this->Shape(), this->_spacing);						// create the output
			for (size_t zi = 0; zi < Z(); zi++) {
				for (size_t yi = 0; yi < Y(); yi++) {
					for (size_t xi = 0; xi < X(); xi++) {
						for (size_t ci = 0; ci < C(); ci++) {
							result(xi, yi, zi, ci) = m_Get(xi, yi, zi, ci) + rhs(xi, yi, zi);
						}
					}
				}
			}
			return result;
		}

		volume<T> operator-(volume<T> rhs) {							// point-wise multiplication
			if (X() != rhs.X() || Y() != rhs.Y())
				throw std::runtime_error("Images dimensions are incompatible");

			tira::volume<T> result(this->Shape(), this->_spacing);						// create the output
			for (size_t zi = 0; zi < Z(); zi++) {
				for (size_t yi = 0; yi < Y(); yi++) {
					for (size_t xi = 0; xi < X(); xi++) {
						for (size_t ci = 0; ci < C(); ci++) {
							result(xi, yi, zi, ci) = m_Get(xi, yi, zi, ci) - rhs(xi, yi, zi);
						}
					}
				}
			}
			return result;
		}

		volume<T> operator/(volume<T> rhs) {							// point-wise multiplication
			if (X() != rhs.X() || Y() != rhs.Y())
				throw std::runtime_error("Images dimensions are incompatible");

			tira::volume<T> result(this->Shape(), this->_spacing);						// create the output
			for (size_t zi = 0; zi < Z(); zi++) {
				for (size_t yi = 0; yi < Y(); yi++) {
					for (size_t xi = 0; xi < X(); xi++) {
						for (size_t ci = 0; ci < C(); ci++) {
							result(xi, yi, zi, ci) = m_Get(xi, yi, zi, ci) / rhs(xi, yi, zi);
						}
					}
				}
			}
			return result;
		}


		void resize(size_t x, size_t y, size_t z, size_t c = 0) {
			std::vector<size_t> shape = { z, y, x, c };
			tira::field<T>::Resize(shape);
		}

		void resize(std::vector<size_t> shape) {
			for (size_t d = shape.size(); d < 4; d++)
				shape.push_back(1);
			field<T>::Resize(shape);
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

		/// Calculates the signed distance field from a binary image
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
			for (size_t i = 0; i < field<T>::m_data.size(); i++) {
				cmap::colormap(field<T>::m_data[i], minval, maxval, color_result.Data()[i * 3 + 0], color_result.Data()[i * 3 + 1], color_result.Data()[i * 3 + 2], colormap);
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
