#pragma once

//external libraries

#include <cmath>
#include <numbers>
#include <typeinfo>

#include "extern/libnpy/npy.hpp"
#include "tira/calculus.h"

namespace tira {

	template <class T>
	class field {
	protected:
		std::vector<T> _data;							// pointer to raw field data
		std::vector<size_t> _shape;						// dimension sizes, specified as the slowest dimension first and the fastest dimension last (C ordering)

		/// <summary>
		/// Allocate data in the _data vector based on the values in _shape
		/// </summary>
		void _allocate() {
			if (_shape.size() == 0) _data.resize(0);
			size_t s = _shape[0];
			for (size_t di = 1; di < _shape.size(); di++)
				s *= _shape[di];
			_data.resize(s);
		}

		void _setShape(std::vector<size_t> S) {
			_shape = S;
			_allocate();
		}

		// test to see if two shape vectors are compatible (same values and same dimensions)
		bool _sameshape(std::vector<size_t> S) {
			if (_shape.size() != S.size()) return false;
			for (size_t i = 0; i < _shape.size(); i++) {
				if (_shape[i] != S[i]) return false;
			}
			return true;
		}

		size_t _idx_offset(size_t d, size_t i) const {				// base function calculates the index offset for the last dimension (called by recursive version)
			size_t off = i;
			for (size_t di = d + 1; di < _shape.size(); di++) {
				off *= _shape[di];
			}
			return off;
		}

		

		

		/// <summary>
		/// Calculate the partial derivative of the field along the specified axis and return a pointer to the resulting data.
		/// </summary>
		/// <param name="axis">Axis along which the partial derivative is calculated</param>
		/// <param name="d">Derivative (ex. 2 for second derivative)</param>
		/// <param name="order">Order of accuracy (requires order+d sample points)</param>
		/// <returns>Pointer to the derivative data in an array that is the same format as the current field</returns>
		T* _derivative_ptr(unsigned int axis, unsigned int d, unsigned int order, bool print_coefs = false) {

			if (axis >= _shape.size()) throw "ERROR: axis out of range of field";
			std::vector< std::vector<double> > C = tira::calculus::finite_difference_coefficients(d, order);		// calculate the list of finite difference coefficients
			if(print_coefs) tira::calculus::printCoefficients(C);

			T* derivative = new T[_data.size()];												// allocate a dynamic array for the derivative data

			int S = C.size();											// store the number of sample points
			int S2 = S / 2;												// store the number of sample points on one side of the template
			int ci;														// current template to use
			T accum;
			int offset;
			std::vector<size_t> c(_shape.size());
			std::vector<size_t> cs(_shape.size());
			iterator i = begin();
			while (i != end()) {
				
				i.coord(&c[0]);								// get the coordinate c for the current iterator location
				if (c[axis] < S2) {							// if the current point is within half of the window size from the starting edge, select a shifted template
					ci = c[axis];
					offset = -ci;
				}
				else if (c[axis] >= _shape[axis] - S2) {	// if the current point is within half of the window size from the trailing edge, select a shifted template
					ci = S - (_shape[axis] - c[axis]);
					offset = -ci;
				}
				else {										// we have enough data on both sides to use the center template
					ci = S2;
					offset = -ci;
				}
				accum = 0;									// initialize the sum for the template to zero

				for (int si = 0; si < S; si++) {			// for each sample in the finite difference template
					cs = c;									// initialize the first template location to the current coordinate
					cs[axis] += si + offset;				// adjust based on which template is being used
					T val = C[ci][si] * read(cs);			// calculate the product of the image and template
					accum += val;							// save the result in the sum
				}
				derivative[idx(c)] = accum;
				++i;														// for each point in the field
			}
			return derivative;

		}


		

	public:

		struct iterator {

			// iterator tags, used to make this more efficient for standard library algorithms
			using iterator_category = std::bidirectional_iterator_tag;
			using difference_type = std::ptrdiff_t;
			using value_type = T;
			using pointer = T*;
			using reference = T&;

			iterator() {
				_ptr = NULL;
				//_coord.assign(0);
				_field = NULL;
			}

			iterator(iterator &c) {
				_ptr = c._ptr;
				_coord = c._coord;
				_field = c._field;
			}

			iterator(field<T>* f) {
				_field = f;
				_ptr = _field->data();
				_coord = std::vector<size_t>(_field->shape().size(), 0);
			}
			
			iterator(field<T>* f, pointer ptr) {
				_field = f;
				_ptr = ptr;
			}

			iterator& operator=(iterator& c) {												// copy assignment operator
				_ptr = c._ptr;
				_coord = c._coord;
				_field = c._field;
			}

			reference operator*() const { return *_ptr; }								// dereference operator
			pointer operator->() { return _ptr; }

			iterator& operator++() {													// increment operator

				_ptr++;

				int D = _field->ndims();
				_coord[D - 1]++;

				//if the position overflows, update
				for (int di = D - 1; di > 0; di--) {
					if (_coord[di] >= _field->_shape[di]) {
						_coord[di] = 0;
						_coord[di - 1]++;
					}
					else break;
				}
				return *this;
			}							
			iterator& operator++(int) { iterator tmp = *this; ++(*this); return tmp; }	// increment operator (prefix)

			friend bool operator== (const iterator& a, const iterator& b) { 
				return a._ptr == b._ptr;
			};
			friend bool operator!= (const iterator& a, const iterator& b) { 
				return a._ptr != b._ptr;
			};
			std::vector<size_t> coord() { return _coord; }
			void coord(size_t* c) {										// copy the coordinate to a pre-allocated array
				memcpy(c, &_coord[0], _coord.size() * sizeof(size_t)); 
			}			
			size_t coord(int axis) { return _coord[axis]; }

			size_t idx() { return _ptr - _field->data(); }


		private:
			pointer _ptr;												// pointer directly to the current element
			std::vector<size_t> _coord;									// N-dimensional coordinate into the field
			field<T>* _field;											// pointer to the field object being iterated over
		};

		field() {}											// default constructor sets up an empty field

		field(std::vector<size_t> shape) {
			_shape = shape;
			_allocate();
		}

		field(std::vector<size_t> shape, T val) : field(shape) {
			size_t N = field<T>::size();
			std::fill(field<T>::_data.begin(), field<T>::_data.end(), val);
		}

		field(std::vector<size_t> shape, T* ptr) {
			_setShape(shape);													// set the new shape of the field
			_allocate();																// allocate space for the field
			memcpy(&_data[0], ptr, bytes());										// copy the data from the							
		}

		field(std::string npy_filename) {
			load_npy(npy_filename);
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

		unsigned int ndims() { return _shape.size(); }

		template<typename D = T>
		void load_npy(const std::string filename) {					// fill the field with data from an npy file
			std::vector<unsigned long> shape;
			std::vector<D> data;
			bool fortran_order;
			try {
				npy::LoadArrayFromNumpy<D>(filename, shape, fortran_order, data);	// load NPY array and metadata
			}
			catch(const std::runtime_error &e){
				std::cout << "ERROR loading NumPy file: " << e.what() << std::endl;
				exit(1);
			}

			// calculate the number of NPY array elements D can fit in a single field object T
			const size_t D_per_T = sizeof(T) / sizeof(D);
			if (D_per_T * sizeof(D) != sizeof(T)) {
				std::cout << "ERROR: Cannot divide the field type <" << typeid(T).name() << "> (" << sizeof(T) << " bytes) into NumPy elements <" << typeid(D).name() << "> (" << sizeof(D) << " bytes)" << std::endl;
				exit(1);
			}

			size_t Dims_per_T = 0;
			if (D_per_T != 1) {																// the field type is composed of multiple array elements, so assume these are broken up into dimensions at the end of the array
				size_t total_elements = 1;										// start with one element
				//for (size_t di = shape.size() - 1; di >= 0; di--) {				// go backwards through the NumPy array dimensions to figure out how many dimensions compose a field type
				for(size_t di = 0; di < shape.size(); di++) {
					const size_t din = shape.size() - di - 1;
					total_elements *= shape[din];
					if (total_elements == D_per_T) {							// if the total number of elements matches the number in the field type
						Dims_per_T = shape.size() - din;							// save the number of dimensions used to represent the type
						break;
					}
				}
				if (total_elements != D_per_T) {
					std::cout << "ERROR: Unable to match field type <" << typeid(T).name() << "> with the trailing dimensions of the NumPy array" << std::endl;
					exit(1);
				}
			}

			const size_t FieldDims = shape.size() - Dims_per_T;						// calculate the total number of dimensions in the field
			std::vector<size_t> new_shape(shape.begin(), shape.begin() + FieldDims);	// calculate the actual shape of the field (composed of elements T, each composed of sub-elements D)

			if (fortran_order) {													// if the numpy array uses fortran ordering, than the indices are flipped
				std::ranges::reverse(new_shape.begin(), new_shape.end());
				_setShape(new_shape);
			}
			else {
				_setShape(new_shape);													// set the new shape of the field
			}
			_allocate();																	// allocate space for the field

			memcpy(reinterpret_cast<void*>(&_data[0]), reinterpret_cast<void*>(&data[0]), bytes());									// copy the data from the
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

		field<T> operator+(field<T> rhs) {
			field<T> result(_shape);
			size_t N = field<T>::size();
			for(size_t ni = 0; ni < N; ni++) {
				result._data[ni] = _data[ni] + rhs._data[ni];
			}
			return result;
		}

		// Binary subtraction (subtract one field from another)
		field<T> operator-(field<T> rhs) {
			if(!_sameshape(rhs._shape))
				throw std::runtime_error("Cannot subtract fields of different shapes");

			field<T> result(_shape);
			size_t N = field<T>::size();
			for(size_t ni = 0; ni < N; ni++) {
				result._data[ni] = _data[ni] - rhs._data[ni];
			}
			return result;
		}

		// Unary negation operator (multiply the field by -1)
		field<T> operator-() {
			field<T> result(_shape);
			for(size_t i = 0; i < _data.size(); i++)
				result._data[i] = -_data[i];
			return result;
		}

		field<T> operator*(field<T> rhs) {
			field<T> result(_shape);
			size_t N = field<T>::size();
			for(size_t ni = 0; ni < N; ni++) {
				result._data[ni] = _data[ni] * rhs._data[ni];
			}
			return result;
		}

		field<T> operator/(field<T> rhs) {
			field<T> result(_shape);
			size_t N = field<T>::size();
			for(size_t ni = 0; ni < N; ni++) {
				result._data[ni] = _data[ni] / rhs._data[ni];
			}
			return result;
		}

		field<T> operator/(T rhs) {
			field<T> result(_shape);
			size_t N = field<T>::size();
			for (size_t ni = 0; ni < N; ni++) {
				result._data[ni] = _data[ni] / rhs;
			}
			return result;
		}

		field<T> operator*(T rhs) {
			field<T> result(_shape);
			size_t N = field<T>::size();
			for(size_t ni = 0; ni < N; ni++) {
				result._data[ni] = _data[ni] * rhs;
			}
			return result;
		}

		field<T> abs() {
			field<T> result(_shape);
			size_t N = field<T>::size();
			for(size_t ni = 0; ni < N; ni++) {
				result._data[ni] = std::abs(_data[ni]);
			}
			return result;
		}

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
				throw "ERROR: index out of range:";
			}
			return _data[i];
		}

		T read(const size_t i) const {
			if (i >= _data.size()) {
				throw "ERROR: index out of range";
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
				result.coord(ri, r_coord);											// get the coordinate in the result array
				// FOR each pixel in the kernel
				for (size_t ki = 0; ki < K.size(); ki++) {							// for each element in the kernel
					K.coord(ki, k_coord);											// get the coordinate in the kernel
					
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
			std::vector<size_t> R = { RY, RX };
			field<T> result(R);

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

		/// <summary>
		/// Pads the field with a constant value
		/// </summary>
		/// <param name="dist	">Vector of paddings in each direction</param>
		/// <param name="val	">The value for the padded values</param>
		field<T> border(std::vector<size_t> dist, T val=0) {
			std::vector<size_t> new_shape = _shape;
			size_t D = _shape.size();

			for (size_t si = 0; si < dist.size(); si++) {
				new_shape[si] += 2 * dist[si];				// calculate new shape 
			}

			field<T> result(new_shape, val);				// initialize result field with val as outer values

			std::vector<size_t> new_coord(D);
			for (size_t ui = 0; ui < size(); ui++) {
				coord(ui, new_coord);						// find new coordinates
				for (int si = 0; si < dist.size(); si++) {
					new_coord[si] += dist[si];
				}
				size_t new_index = result.idx(new_coord);	//find new index
				result._data[new_index] = _data[ui];		// replace the values from the initial field
			}
			return result;
		}

		/// <summary>
		/// Pads the field with a constant value
		/// </summary>
		/// <param name="dist	">Padding length</param>
		/// <param name="val	">The value for the padded values</param>
		field<T> border(size_t dist, T val = 0) {
			std::vector<size_t> dist_vect = { dist, dist };
			return border(dist_vect, val);
		}

		/// <summary>
		/// Generates a border around the field by replicating edge values
		/// </summary>
		/// <param name="dist">A vector of padding distances in each axis direction</param>
		/// <returns></returns>
		field<T> border_replicate(std::vector<size_t> dist) {
			std::vector<size_t> new_shape = _shape;
			size_t D = _shape.size();
			std::vector<size_t> coord_i(D);					// initialize memory for the coordinates in the cycle

			for (size_t si = 0; si < D; si++) {
				new_shape[si] += 2 * dist[si];				// calculate new shape 
			}

			field<T> result(new_shape);						// initialize result field 

			size_t total_size = result.size();
			for (size_t index = 0; index < total_size; index++) {				// for each element in resulting field
				
				result.coord(index, coord_i);									// convert index to coordinates
				for (size_t di = 0; di < D; di++) {
					coord_i[di] = std::max(dist[di], coord_i[di]);				// bound in range [dist; n + 2 * dist)
					coord_i[di] -= dist[di];									// move the range to [0; n + dist)
					coord_i[di] = std::min(_shape[di] - 1, coord_i[di]);		// bound in range [0; n)
				}	
				result._data[index] = _data[idx(coord_i)];						
			}
			return result;
		}

		/// <summary>
		/// Generates a border around the field with constant padding distance by replicating edge values
		/// </summary>
		/// <param name="dist">Length of padding</param>
		/// <returns></returns>
		field<T> border_replicate(size_t dist) {
			std::vector<size_t> dist_vect(_shape.size(), dist);
			return border_replicate(dist_vect);
		}


		field<T> crop(std::vector<size_t> min_coord, std::vector<size_t> max_coord) const {

			size_t D = _shape.size();												// get the number of dimensions
			std::vector<size_t> new_dims(D);
			for (size_t di = 0; di < D; di++) {										// calculate the size of the output field
				new_dims[di] = max_coord[di] - min_coord[di];
			}

			field<T> result(new_dims);
			std::vector<size_t> s_coord(_shape.size());								// stores the source field coordinate
			std::vector<size_t> r_coord(_shape.size());								// stores the result field coordinate
			// FOR each pixel in the output (result) image
			for (size_t ri = 0; ri < result.size(); ri++) {							// for each element in the result field
				result.coord(ri, r_coord);											// get the coordinate in the result array

				// CALCULATE the coordinate into the source image
				for (size_t ci = 0; ci < s_coord.size(); ci++) {				// for each coordinate
					s_coord[ci] = r_coord[ci] + min_coord[ci];					// offset the coordinate by the kernel coordinate
				}
				size_t s_idx = idx(s_coord);									// calculate the index into the source field
				T s_val = _data[s_idx];											// get the value from the source field

				result(r_coord) = s_val;											// store the accumulated value in the result field
			}

			return result;															// return the result field
		}



		void resize(std::vector<size_t> new_size) {
			_setShape(new_size);
		}


		T* data() {
			return &_data[0];
		}

		iterator begin() { 
			return iterator(this);
		}
		iterator end() { 
			return iterator(this, data() + size());
			//return iterator(this);
		}

		field<T> derivative(unsigned int axis, unsigned int d, unsigned int order, bool print_coefs = false) {

			T* ptr = _derivative_ptr(axis, d, order, print_coefs);

			field<T> result(_shape, ptr);

			delete ptr;
			return result;
		}

		static std::vector<double> fd_coefficients(unsigned int d, unsigned int order) {
			return tira::calculus::central_difference_coefficients(d, order);
		}

		field<T> central_derivative(unsigned int axis, unsigned int d, unsigned int order, bool pad = true, bool print_coefs = false) {

			std::vector<double> C = tira::calculus::central_difference_coefficients(d, order);
			if(print_coefs) tira::calculus::printCoefficients(C);

			std::vector<size_t> kernel_shape(_shape.size(), 1);
			kernel_shape[axis] = C.size();
			tira::field<T> kernel(kernel_shape);

			for(size_t ci = 0; ci < C.size(); ci++) {
				kernel(ci) = C[ci];
			}

			tira::field<T> D;
			if (pad) {
				std::vector<size_t> dist = kernel_shape;
				for (int si = 0; si < dist.size(); si++) {
					dist[si] /= 2;							// divide the shape of kernel by 2 to get the needed padding
				}
				D = border_replicate(dist).convolve(kernel);
			}
			else {
				D = convolve(kernel);
			}
			return D;
		}

		
	};
}