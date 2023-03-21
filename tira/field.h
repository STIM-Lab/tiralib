#pragma once

#include <typeinfo>

#include "extern/libnpy/npy.hpp"

namespace tira {

	template <class T>
	class field {
	protected:
		std::vector<T> _data;							// pointer to raw field data
		std::vector<size_t> _shape;						// dimension sizes, specified as the slowest dimension first and the fastest dimension last (C ordering)

		/// <summary>
		/// Allocate data in the _data vector based on the values in _shape
		/// </summary>
		void allocate() {
			if (_shape.size() == 0) _data.resize(0);
			size_t s = _shape[0];
			for (size_t di = 1; di < _shape.size(); di++)
				s *= _shape[di];
			_data.resize(s);
		}

		void setShape(std::vector<size_t> S) {
			_shape = S;
			allocate();
		}

		size_t idx_offset(size_t d, size_t i) const {				// base function calculates the index offset for the last dimension (called by recursive version)
			size_t off = i;
			for (size_t di = d + 1; di < _shape.size(); di++) {
				off *= _shape[di];
			}
			return off;
		}
		/*template<typename... D>
		size_t idx_offset(size_t d, size_t i, D... more) const {	// recursive function for finding the index from an arbitrary set of coordinates
			size_t off = i;
			for (size_t di = d + 1; di < _shape.size(); di++) {
				off *= _shape[di];
			}
			return off + idx_offset(d + 1, more...);
		}*/

		/// <summary>
		/// Calculate the 2D index into the _data array given a set of coordinates using image-based indexing (specified fast-to-slow)
		/// </summary>
		/// <param name="coords"></param>
		/// <returns></returns>
		size_t idx_f2s(std::vector<size_t> coords) {

		}

		

	public:

		field() {}											// default constructor sets up an empty field

		/*template<typename... D>
		field(const D... d) {
			setShape(d...);
			allocate();										// allocate space to store the field
		}*/
		field(std::vector<size_t> shape) {
			_shape = shape;
			allocate();
		}

		/// <summary>
		/// Copy constructor
		/// </summary>
		/// <param name="other"></param>
		field(const field<T>& other) {
			_shape = other._shape;
			_data = other._data;
		}

		/// <summary>
		/// Convert a 1D index to an ND field coordinate
		/// </summary>
		/// <param name="idx"></param>
		/// <returns></returns>
		inline void coord(size_t idx, std::vector<size_t> &c) const {

			size_t nC = c.size();
			size_t step;
			size_t i, j;
			for (i = 0; i < nC; i++) {
				step = 1;
				for (j = i + 1; j < nC; j++) {
					step *= _shape[j];
				}
				c[i] = idx / step;
				idx -= c[i] * step;
			}
		}

		size_t idx(const std::vector<size_t> &coord) const {
			size_t i = 0;

			size_t step;
			for (size_t c = 0; c < coord.size(); c++) {
				step = 1;
				for (size_t j = c + 1; j < coord.size(); j++) {
					step *= _shape[j];
				}
				i += coord[c] * step;
			}
			return i;
		}

		template<typename D = T>
		void load_npy(std::string filename) {					// fill the field with data from an npy file
			std::vector<unsigned long> shape;
			std::vector<D> data;
			bool fortran_order;
			try {
				npy::LoadArrayFromNumpy<D>(filename, shape, fortran_order, data);	// load NPY array and metadata
			}
			catch(const std::runtime_error &e){
				std::cout << e.what() << std::endl;
				exit(1);
			}

			// calculate the number of NPY array elements D can fit in a single field object T
			size_t D_per_T = sizeof(T) / sizeof(D);
			if (D_per_T * sizeof(D) != sizeof(T)) {
				std::cout << "ERROR: Cannot divide the field type <" << typeid(T).name() << "> (" << sizeof(T) << " bytes) into NumPy elements <" << typeid(D).name() << "> (" << sizeof(D) << " bytes)" << std::endl;
				exit(1);
			}

			size_t Dims_per_T;
			if (D_per_T == 1) Dims_per_T = 0;									// if the size of the field type matches individual array elements, just duplicate the field dimensions
			else {																// the field type is composed of multiple array elements, so assume these are broken up into dimensions at the end of the array
				size_t total_elements = 1;										// start with one element
				for (size_t di = shape.size() - 1; di >= 0; di--) {				// go backwards through the NumPy array dimensions to figure out how many dimensions compose a field type
					total_elements *= shape[di];
					if (total_elements == D_per_T) {							// if the total number of elements matches the number in the field type
						Dims_per_T = shape.size() - di;							// save the number of dimensions used to represent the type
						break;
					}
				}
				if (total_elements != D_per_T) {
					std::cout << "ERROR: Unable to match field type <" << typeid(T).name() << "> with the trailing dimensions of the NumPy array" << std::endl;
					exit(1);
				}
			}

			size_t FieldDims = shape.size() - Dims_per_T;						// calculate the total number of dimensions in the field
			std::vector<size_t> new_shape(shape.begin(), shape.begin() + FieldDims);	// calculate the actual shape of the field (composed of elements T, each composed of sub-elements D)

			if (fortran_order) {													// if the numpy array uses fortran ordering, than the indices are flipped
				std::reverse(new_shape.begin(), new_shape.end());
				setShape(new_shape);
			}
			else {
				setShape(new_shape);													// set the new shape of the field
			}
			allocate();																	// allocate space for the field

			memcpy(&_data[0], &data[0], bytes());									// copy the data from the

		}
		
		void save_npy(const std::string& filename) {
			bool fortran_order = false;
			std::vector<unsigned long> shape(_shape.size());
			for (int i = 0; i < shape.size(); i++)
				shape[i] = _shape[i];

			npy::SaveArrayAsNumpy(filename, fortran_order, shape.size(), (const unsigned long *) & shape[0], &_data[0]);
		}


		size_t bytes() const {
			return size() * sizeof(T);
		}

		std::vector<size_t> shape() const { return _shape; }

		/// <summary>
		/// Re-shape the current field based on the specified parameters
		/// </summary>
		/// <param name="new_shape"></param>
		void reshape(std::vector<size_t> new_shape) {
			size_t n = new_shape[0];
			for (size_t i = 1; i < new_shape.size(); i++)
				n *= new_shape[i];
			if (n != size()) {
				std::cout << "ERROR: new shape (";
				for (size_t i = 0; i < new_shape.size(); i++)
					std::cout << new_shape[i];
				std::cout << ") and old shape (";
				for (size_t i = 0; i < _shape.size(); i++)
					std::cout << _shape[i];
				std::cout << ") must have the same number of elements." << std::endl;
			}
			_shape = new_shape;
		}

		// return the total number of elements in the field
		size_t size() const {							
			
			return _data.size();
		}

		/// <summary>
		/// Set all values in the image to a single constant
		/// </summary>
		/// <param name="v">Constant that all elements will be set to</param>
		field<T> operator=(T v) {														//set all elements of the image to a given value v
			size_t N = field<T>::size();
			std::fill(field<T>::_data.begin(), field<T>::_data.end(), v);
			return *this;
		}

		/*template<typename... D>
		T& operator()(size_t x, D... more) {							// returns a reference to the indexed value
			return _data[idx_offset(0, x, more...)];
		}*/

		T& operator()(std::vector<size_t> x) {
			return _data[idx(x)];
		}

		T& operator()(size_t i) {
			return _data[i];
		}

		

		/// <summary>
		/// Read-only version of the accessor method
		/// </summary>
		/// <param name="x"></param>
		/// <returns></returns>
		T read(const std::vector<size_t> &x) const {
			size_t i = idx(x);
			if (i >= _data.size()) {
				std::cout << "ERROR: index out of range: i = " << i << " (" << _data.size() << " max)" << std::endl;
				exit(1);
			}
			return _data[i];
		}

		T read(const size_t i) const {
			if (i >= _data.size()) {
				std::cout << "ERROR: index out of range: i = " << i << " (" << _data.size() << " max)" << std::endl;
				exit(1);
			}
			return _data[i];
		}


		/// <summary>
		/// Returns the convolution of the current field with a given kernel
		/// </summary>
		/// <param name="mask"></param>
		/// <returns></returns>
		template<class C>
		field<T> convolve(const field<C> K) const {

			// make sure the number of dimensions in k matches the number of dimensions in the current field
			if (K.shape().size() > _shape.size()) {
				std::cout << "error: convolution kernel (d = " << K.shape().size() << ") has more dimensions than the field (d = " << _shape.size() << ")" << std::endl;
				exit(1);
			}

			// calculate the size of the output image and allocate space
			std::vector<size_t> result_shape = _shape;												// start with the size of the source image
			for (size_t i = 0; i < result_shape.size(); i++) result_shape[i] -= K.shape()[i] - 1;	// offset the size by the size of the kernel (K-1)
			field<T> result(result_shape);															// generate a new field to store the convolved image

			std::vector<size_t> k_coord(_shape.size());
			std::vector<size_t> s_coord(_shape.size());
			std::vector<size_t> r_coord(_shape.size());
			T sum;																	// create a T value to accumulate the convolution result for a single pixel
			// FOR each pixel in the output (result) image
			for (size_t ri = 0; ri < result.size(); ri++) {							// for each element in the result field
				sum = (T)0;
				result.coord(ri, r_coord);										// get the coordinate in the result array
				// FOR each pixel in the kernel
				for (size_t ki = 0; ki < K.size(); ki++) {							// for each element in the kernel
					K.coord(ki, k_coord);										// get the coordinate in the kernel
					
					// CALCULATE the coordinate into the source image
					for (size_t ci = 0; ci < s_coord.size(); ci++) {				// for each coordinate
						s_coord[ci] = r_coord[ci] + k_coord[ci];					// offset the coordinate by the kernel coordinate
					}
					size_t s_idx = idx(s_coord);									// calculate the index into the source image
					T s_val = _data[s_idx];											// get the value from the source field
					C kval = K.read(k_coord);										// get the value of the kernel
					sum = sum + kval * s_val;										// accumulate the product of the kernel and the source
				}
				result(r_coord) = sum;												// store the accumulated value in the result field
			}

			return result;															// return the result field

		}

		

		template<class C>
		field<T> convolve2(const field<C> K) const {

			// make sure all of the dimensions are correct
			if (_shape.size() != 2) {
				std::cout << "ERROR: cannot perform a 2D convolution on a field that is not 2D (D = " << _shape.size() << ")" << std::endl;
				exit(1);
			}
			if (K.shape().size() != 2) {
				std::cout << "ERROR: cannot perform a 2D convolution with a kernel that is not 2D (D = " << K.shape().size() << ")" << std::endl;
				exit(1);
			}

			// calculate the size of the output image (only valid region is returned) and allocate the result field
			std::vector<size_t> result_shape = _shape;												
			size_t RY = _shape[0] - (K.shape()[0] - 1);
			size_t RX = _shape[1] - (K.shape()[1] - 1);
			field<T> result({ RY, RX });

			// store the dimensions of the kernel and source image
			size_t KY = K.shape()[0];
			size_t KX = K.shape()[1];
			size_t SY = _shape[0];
			size_t SX = _shape[1];

			size_t yi, xi, kxi, kyi, result_i, source_i, k_i;			// allocate indices
			T input_v, result_v;
			C k_v;

			// FOR each pixel along the Y direction in the output (result) field
			for (yi = 0; yi < RY; yi++) {
				// FOR each pixel along the X direction in the output (result) field
				for (xi = 0; xi < RX; xi++) {
					result_i = yi * RX + xi;							// calculate the pixel index at the current coordinate
					result_v = (T)0.0;									// initialize the convolved value at this pixel to zero
					// FOR each pixel along the X direction in the kernel
					for (kyi = 0; kyi < KY; kyi++) {
						// FOR each pixel along the Y direction in the kernel
						for (kxi = 0; kxi < KX; kxi++) {
							k_i = kyi * KX + kxi;						// calculate the kernel index
							source_i = (yi + kyi) * SX + (xi + kxi);	// calculate the corresponding index into the source image
							k_v = K.read(k_i);							// read the value of the kernel
							input_v = read(source_i);					// read the value of the source pixel
							result_v += input_v *k_v;					// accumulate the result
						}
					}
					result._data[result_i] = result_v;					// store the result
				}
			}
			return result;

		}

		void resize(std::vector<size_t> new_size) {
			setShape(new_size);
		}


		T* data() {
			return &_data[0];
		}
	};
}