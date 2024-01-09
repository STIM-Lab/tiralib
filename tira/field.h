#pragma once

//external libraries

// Eigen is used to solve the linear system required to calculate derivatives at any specified order
#include <Eigen/Dense>


#include <cmath>
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

		/// <summary>
		/// Calculate the finite difference coefficients given a derivative and set of sample points
		/// </summary>
		/// <param name="derivative"></param>
		/// <param name="samples"></param>
		/// <returns></returns>
		Eigen::VectorX<T> finite_difference_coefficients(unsigned int derivative, Eigen::VectorX<T> samples) {

			unsigned int N = samples.size();

			Eigen::MatrixX<T> S(N, N);
			for (unsigned int ri = 0; ri < N; ri++) {
				for (unsigned int ci = 0; ci < N; ci++) {
					S(ri, ci) = pow(samples[ci], ri);
				}
			}

			Eigen::VectorX<T> b = Eigen::VectorX<T>::Zero(N);
			b(derivative) = tgamma(derivative + 1);

			return S.colPivHouseholderQr().solve(b);
		}

		/// <summary>
		/// Calculate the finite difference coefficients given a derivative and order of accuracy
		/// </summary>
		/// <param name="derivative"></param>
		/// <param name="order"></param>
		/// <returns></returns>
		std::vector< std::vector<T> > finite_difference_coefficients(unsigned int derivative, unsigned int order) {

			unsigned int N = order + derivative;		// calculate the number of samples required to achieve the desired order

			std::vector< std::vector<T> > Coefficients;

			Eigen::VectorX<T> Samples(N);				// allocate a vector that will be used to store sample points

			for (int ri = 0; ri < N; ri++) {			// for each shifted sample position
				for (int ci = 0; ci < N; ci++) {		// calculate the point for each sample
					Samples(ci) = -ri + ci;				// store that point in the Samples vector
				}
				std::vector<T> c(N);
				Eigen::Map< Eigen::VectorX<T> >(&c[0], N) = finite_difference_coefficients(derivative, Samples);
				Coefficients.push_back(c);
			}
			return Coefficients;
		}

		/// <summary>
		/// Calculate the partial derivative of the field along the specified axis and return a pointer to the resulting data.
		/// </summary>
		/// <param name="axis">Axis along which the partial derivative is calculated</param>
		/// <param name="d">Derivative (ex. 2 for second derivative)</param>
		/// <param name="order">Order of accuracy (requires order+d sample points)</param>
		/// <returns>Pointer to the derivative data in an array that is the same format as the current field</returns>
		T* derivative_ptr(unsigned int axis, unsigned int d, unsigned int order) {

			if (axis >= _shape.size()) throw "ERROR: axis out of range of field";
			std::vector< std::vector<T> > C = finite_difference_coefficients(d, order);		// calculate the list of finite difference coefficients

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
				
				//c = i.coord();											// calculate the coordinate of the current axis
				i.coord(&c[0]);
				if (c[axis] < S2) {										// if the current point is within half of the window size from the edge, select a shifted template
					ci = c[axis];
					offset = -ci;
				}
				else if (c[axis] >= _shape[axis] - S2) {				// if the current point is within half of the window size from the edge, select a shifted template
					ci = S - (_shape[axis] - c[axis]);
					offset = -ci;
				}
				else {
					ci = S2;
					offset = -ci;
				}
				accum = 0;

				for (int si = 0; si < S; si++) {						// for each sample in the finite difference template
					cs = c;
					cs[axis] += si + offset;
					accum += C[ci][si] * read(cs);
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

			iterator(iterator &c) {
				//_ptr = c._ptr;
				//_idx = c._idx;
				_coord = c._coord;
				//_shape = c._shape;
				//_N = c._N;
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
				//_coord = std::vector<size_t>(_field->shape().size(), 0);
			}

			iterator& operator=(iterator& c) {												// copy assignment operator
				_ptr = c._ptr;
				//_idx = c._idx;
				_coord = c._coord;
				_field = c._field;
				//_shape = c._shape;
				//_N = c._N;
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
			allocate();
		}

		field(std::vector<size_t> shape, T* ptr) {
			setShape(shape);													// set the new shape of the field		
			allocate();																// allocate space for the field
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
		void load_npy(std::string filename) {					// fill the field with data from an npy file
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

		iterator begin() { 
			return iterator(this);
		}
		iterator end() { 
			return iterator(this, data() + size());
			//return iterator(this);
		}

		field<T> derivative(unsigned int axis, unsigned int d, unsigned int order) {

			T* ptr = derivative_ptr(axis, d, order);

			field<T> result(_shape, ptr);

			delete ptr;
			return result;


		}

		
	};
}