#pragma once

#include <cmath>
#include <typeinfo>
#include <fstream>
#include <numbers>
#include "extern/libnpy/npy.hpp"
#include "tira/calculus.h"

namespace tira {
	/**
	 * The field class provides a container for an N-dimensional uniform grid. This is used as a base class for more specialized
	 * grids such as images (tira::image) and volumes (tira::volume). We try to put general functions here that are common to
	 * uniform grids, such as gradients, convolutions, and get/set methods.
	 *
	 * The field class uses an std::vector as its foundation.
	 *
	 * @brief This class represents a generic n-dimensional uniform grid and provides common methods for processing that type of
	 * data. This includes calculating gradients, convolutions, and saving/loading data using general file formats like
	 * RAW and NumPy.
	 * @tparam T is the data type used to represent the data at each sample point
	 */
	template <class T>
	class field {
	protected:
		/**
		 * @brief The m_data attribute stores the raw data as an std::vector. Data is stored "slowest dimension first" such that
		 * iterating through this array would cause the corresponding array indices for the first dimensions to move the slowest
		 */
		std::vector<T> m_data;

		/**
		 * @brief The m_shape attribute stores the number of samples in each grid dimension. The dimensions are stored "slowest first"
		 * such that, when iterating across the m_data array the indices for the first dimensions move the slowest. This is the traditional
		 * order used for C/C++ and Python
		 */
		std::vector<size_t> m_shape;


		/**
		 * @brief Allocates the memory required to store a field based on the sample sizes stored in m_shape
		 */
		void m_Allocate() {
			if (m_shape.size() == 0) m_data.resize(0);
			size_t s = m_shape[0];
			for (size_t di = 1; di < m_shape.size(); di++)
				s *= m_shape[di];
			m_data.resize(s);
		}

		/**
		 * @brief Sets the shape of the field based on an array of inputs and allocates the necessary space to store the data
		 * @param S is an std::vector specifying the number of samples stored in the field for each dimension
		 */
		void m_SetShape(const std::vector<size_t> S) {
			m_shape = S;
			m_Allocate();
		}

		/**
		 * @brief Determines whether or not the specified shape vector is identical to the current shape. This is used to
		 * test whether or not element-to-element operations can be performed between two fields.
		 * @param S is the shape that will be compared to the current object
		 * @return
		 */
		bool m_IsSameShape(const std::vector<size_t> S) const {
			if (m_shape.size() != S.size()) return false;			// if the number of dimensions are different, the fields are not the same shape
			for (size_t i = 0; i < m_shape.size(); i++) {			// for each dimension
				if (m_shape[i] != S[i]) return false;				// if the size of the dimensions are different, the fields are not the same shape
			}
			return true;											// if all shape values are identical, the fields are the same shape
		}


		/**
		 * @brief Calculate the partial derivative of the field along the specified axis and return a pointer to the resulting data.
		 * @param axis is the axis along which the partial derivative is calculated
		 * @param d is the derivative to calculate (ex. 2 for a second derivative)
		 * @param order is the order of accuracy (requires order + d sample points)
		 * @param print_coefs prints the coefficients to the console
		 * @return a pointer to the derivative data in an array that is the same format as the current field
		 */
		T* m_Derivative(unsigned int axis, unsigned int d, unsigned int order, bool print_coefs = false) {

			if (axis >= m_shape.size()) throw "ERROR: axis out of range of field";
			std::vector< std::vector<double> > C = tira::calculus::finite_difference_coefficients(d, order);		// calculate the list of finite difference coefficients
			if(print_coefs) tira::calculus::printCoefficients(C);

			T* derivative = new T[m_data.size()];												// allocate a dynamic array for the derivative data

			int S = C.size();											// store the number of sample points
			int S2 = S / 2;												// store the number of sample points on one side of the template
			int ci;														// current template to use
			T accum;
			int offset;
			std::vector<size_t> c(m_shape.size());
			std::vector<size_t> cs(m_shape.size());
			iterator i = Begin();
			while (i != End()) {
				
				i.coord(&c[0]);								// get the coordinate c for the current iterator location
				if (c[axis] < S2) {							// if the current point is within half of the window size from the starting edge, select a shifted template
					ci = c[axis];
					offset = -ci;
				}
				else if (c[axis] >= m_shape[axis] - S2) {	// if the current point is within half of the window size from the trailing edge, select a shifted template
					ci = S - (m_shape[axis] - c[axis]);
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
					T val = C[ci][si] * Read(cs);			// calculate the product of the image and template
					accum += val;							// save the result in the sum
				}
				derivative[Idx(c)] = accum;
				++i;														// for each point in the field
			}
			return derivative;

		}


		

	public:
		/**
		 * @brief This structure defines an iterator that facilitates traversal of points in the field
		 */
		struct iterator {

			// iterator tags, used to make this more efficient for standard library algorithms
			using iterator_category = std::bidirectional_iterator_tag;
			using difference_type = std::ptrdiff_t;
			using value_type = T;
			using pointer = T*;
			using reference = T&;

			/**
			 * Constructor initializes an empty iterator that is unassociated with a specific field or sample point
			 */
			iterator() {
				m_ptr = NULL;
				m_field = NULL;
			}

			/**
			 * @brief Copy constructor generates a new iterator that is a copy of an existing iterator c
			 * @param c is an existing iterator to be copied
			 */
			iterator(iterator &c) {
				m_ptr = c.m_ptr;
				m_coord = c.m_coord;
				m_field = c.m_field;
			}

			/**
			 * @brief Create an iterator associated with an existing field f. The iterator will initially point to the first sample in the field (0, 0, ...)
			 * @param f is the field that this iterator will be associated with
			 */
			explicit iterator(field* f) {
				m_field = f;
				m_ptr = m_field->Data();
				m_coord = std::vector<size_t>(m_field->Shape().size(), 0);
			}

			/**
			 * @brief Create an iterator that points to a specified location within an existing field
			 * @param f is the field that this iterator will be associated with
			 * @param ptr is a pointer to the value in the field
			 */
			iterator(field* f, pointer ptr) {
				m_field = f;
				m_ptr = ptr;
			}

			/**
			 * @brief Copy assignment operator creates a new copy of an existing iterator during an assignment operation.
			 * @param c is the iterator to copy during assignment
			 * @return a reference to the new iterator that will be stored in the left-hand-side operand of an assignment operation
			 */
			iterator& operator=(const iterator& c) {
				if (this != &c) {
					m_ptr = c.m_ptr;
					m_coord = c.m_coord;
					m_field = c.m_field;
				}
				return (*this);
			}

			/**
			 * @brief Defines the dereferencing operator for a field::iterator
			 * @return a reference to the value that the current iterator is pointing to
			 */
			reference operator*() const { return *m_ptr; }

			/**
			 * @brief Defines the pointer operator
			 * @return a reference to the value that the current iterator is pointing to
			 */
			pointer operator->() { return m_ptr; }

			/**
			 * @brief Prefix increment operator, increments the iterator lexicographically through the field by one sample point.
			 * @return a reference to the current iterator after it is incremented
			 */
			iterator& operator++() {

				++m_ptr;

				const int D = m_field->Ndims();
				++m_coord[D - 1];

				//if the position overflows, update
				for (int di = D - 1; di > 0; --di) {
					if (m_coord[di] >= m_field->m_shape[di]) {
						m_coord[di] = 0;
						++m_coord[di - 1];
					}
					else break;
				}
				return *this;
			}

			/**
			 * @brief Postfix increment operator, increments the iterator lexicographically through the field by one sample point.
			 * @return Returns a reference to the iterator after it has been incremented
			 */
			iterator operator++(int) {
				iterator tmp = *this;
				++(*this);
				return tmp;
			}

			/**
			 * @brief Boolean comparison operator, compares if two iterators are pointing to the same sample point within the same field
			 * @param a is the left-hand-side iterator to compare
			 * @param b is the right-hand-side iterator to compare
			 * @return a boolean value that is true if both iterators are pointing to the same location and false otherwise
			 */
			friend bool operator== (const iterator& a, const iterator& b) {
				return a.m_ptr == b.m_ptr;
			};

			/**
			 * @brief Boolean comparison operator, compares if two iterators are NOT pointing to the same sample point within the same field
			 * @param a is the left-hand-side iterator to compare
			 * @param b is the right-hand-side iterator to compare
			 * @return a boolean value that is true if both iterators are pointing to the same location and false otherwise
			 */
			friend bool operator!= (const iterator& a, const iterator& b) {
				return a.m_ptr != b.m_ptr;
			};

			/**
			 * @brief Converts the current iterator to a valid N-dimensional coordinate within the field
			 * @return a valid N-D coordinate as an std::vector
			 */
			std::vector<size_t> coord() { return m_coord; }

			/**
			 * @brief Copies the coordinate for the current iterator position into a pre-allocated array
			 * @param c is the location of the array to which the current iterator position will be copied
			 */
			void coord(size_t* c) const {
				memcpy(c, &m_coord[0], m_coord.size() * sizeof(size_t));
			}

			/**
			 * @brief Returns the coordinate of the current iterator position for a single axis
			 * @param axis is the axis coordinate that will be returned
			 * @return the coordinate for the specified axis
			 */
			size_t coord(int axis) { return m_coord[axis]; }

			/**
			 * @brief Returns a one-dimensional index to the data point referenced by the iterator. This 1D index is the
			 * index value that can be used to directly access the element in the private field std::vector.
			 * @return a one-dimensional index into the field vector
			 */
			size_t idx() { return m_ptr - m_field->Data(); }


		private:
			/**
			 * @brief Pointer that directly points to the current element in the field
			 */
			pointer m_ptr;

			/**
			 * @brief N-dimensional coordinate into the field that is actively maintained as the iterator moves
			 */
			std::vector<size_t> m_coord;

			/**
			 * @brief Pointer to the field object that is being iterated across
			 */
			field<T>* m_field;
		};

		/**
		 * @brief Default constructor that creates an empty field with no data points or dimensions
		 */
		field() {}

		/**
		 * @brief Constructor creates a field with the specified shape and allocates the necessary space. The contents of the
		 * field are not defined at creation.
		 * @param shape is the shape of the field to be created
		 */
		explicit field(const std::vector<size_t> shape) {
			m_shape = shape;
			m_Allocate();
		}

		/**
		 * @brief Constructor that creates a field with the specified shape and assignes the value val to all points
		 * @param shape is the shape of the field to be created
		 * @param val is the value at all sample points
		 */
		field(std::vector<size_t> shape, T val) : field(shape) {
			std::fill(m_data.begin(), m_data.end(), val);
		}

		/**
		 * @brief Constructor creates a field with the specified shape and fills it with a provided set of data points.
		 * @param shape is the shape of the field to be created
		 * @param ptr is a pointer to the values that will be copied into the field (ptr must point to an appropriate number of data points to fill the field)
		 */
		field(const std::vector<size_t> shape, const T* ptr) {
			m_SetShape(shape);													// set the new shape of the field
			m_Allocate();															// allocate space for the field
			memcpy(&m_data[0], ptr, Bytes());										// copy the data from the
		}

		/**
		 * @brief Constructor creates a new field using the contents of a NumPy file
		 * @param npy_filename is the name of the NumPy file to be loaded
		 */
		explicit field(const std::string npy_filename) {
			LoadNpy(npy_filename);
		}

		/**
		 * @brief Copy constructor copies an existing field into the current field (including the shape and all data)
		 * @param other is the constructor to be copied
		 */
		field(const field<T>& other) {
			m_shape = other.m_shape;
			m_data = other.m_data;
		}

		/**
		 * @brief Assignment operator fills the current field with the data from another field
		 * @param other is the right-hand-side operator in the assignment and the field containing the data to be copied
		 * @return a reference to the current field after the data has been copied from other
		 */
		field& operator=(const field& other) {
			m_data = other.m_data;
			m_shape = other.m_shape;
			return *this;
		}

		/**
		 * @brief Assignment operator sets all values in the field to the specified value
		 * @param v the value that will be assigned to all values in the field
		 * @return a field with the specified value at all sample points
		 */
		field& operator=(const T& v) {														//set all elements of the image to a given value v
			std::fill(m_data.begin(), m_data.end(), v);
			return *this;
		}

		/**
		 * @brief Addition operator performs an element-by-element addition of two fields that have the same shape.
		 * @param rhs is the field that will be added to this field
		 * @return a new field that is the sum of both fields
		 */
		field operator+(field<T> rhs) {
			field result(m_shape);
			size_t N = Size();
			for(size_t ni = 0; ni < N; ni++) {
				result.m_data[ni] = m_data[ni] + rhs.m_data[ni];
			}
			return result;
		}

		/**
		 * @brief Subtraction operator performs an element-by-element addition of two fields that have the same shape.
		 * @param rhs is the field that will be subtracted from this field
		 * @return a new field that is the difference of both fields
		 */
		field<T> operator-(field<T> rhs) {
			if(!m_IsSameShape(rhs.m_shape))
				throw std::runtime_error("Cannot subtract fields of different shapes");

			field result(m_shape);
			size_t N = Size();
			for(size_t ni = 0; ni < N; ni++) {
				result.m_data[ni] = m_data[ni] - rhs.m_data[ni];
			}
			return result;
		}

		/**
		 * @brief Unitary negative operator that calculates the negative of each value in the field
		 * @return a new field where each element is the negative of the corresponding element in the current field
		 */
		field operator-() {
			field result(m_shape);
			for(size_t i = 0; i < m_data.size(); i++)
				result.m_data[i] = -m_data[i];
			return result;
		}

		/**
		 * @brief Performs a piecewise multiplication between two fields
		 * @param rhs is the right-hand-side operator in the element-wise multiplication
		 * @return a new field that is the element-by-element product of the current field and rhs
		 */
		field operator*(field rhs) {
			field result(m_shape);
			size_t N = Size();
			for(size_t ni = 0; ni < N; ni++) {
				result.m_data[ni] = m_data[ni] * rhs.m_data[ni];
			}
			return result;
		}

		/**
		 * @brief Performs a piecewise division between two fields
		 * @param rhs is the right-hand-side operator in the element-wise division
		 * @return a new field that is the current field divided by rhs
		 */
		field operator/(field rhs) {
			field result(m_shape);
			size_t N = Size();
			for(size_t ni = 0; ni < N; ni++) {
				result.m_data[ni] = m_data[ni] / rhs.m_data[ni];
			}
			return result;
		}

		/**
		 * @brief Divides the entire field by a scalar value
		 * @param rhs is the scalar value that will be divided into every element of the field
		 * @return a new field where each element of the current field is divided by rhs
		 */
		field operator/(T rhs) {
			field result(m_shape);
			size_t N = Size();
			for (size_t ni = 0; ni < N; ni++) {
				result.m_data[ni] = m_data[ni] / rhs;
			}
			return result;
		}

		/**
		 * @brief Multiplies the entire field by a scalar value
		 * @param rhs is the scalar value that will be multiplied by every element of the field
		 * @return a new field where each element of the current field is multiplied by rhs
		 */
		field operator*(T rhs) {
			field result(m_shape);
			size_t N = Size();
			for(size_t ni = 0; ni < N; ni++) {
				result.m_data[ni] = m_data[ni] * rhs;
			}
			return result;
		}

		/**
		 * @brief Access the specified element of the field
		 * @param x N-dimensional coordinate of the desired element
		 * @return a reference to the value specified by the N-dimensional coordinate x
		 */
		T& operator()(std::vector<size_t> x) {
			return m_data[Idx(x)];
		}

		/**
		 * @brief Access the specified element of the field using a one-dimensional index
		 * @param i is the one-dimensional index used to access an element in m_data
		 * @return a reference to the corresponding value in the scalar field
		 */
		T& operator()(size_t i) {
			return m_data[i];
		}

		/**
		 * @brief Calculates the absolute value of each element in the field
		 * @return a new field where each point is the absolute value of the corresponding point in this field
		 */
		field Abs() {
			field result(m_shape);
			size_t N = Size();
			for(size_t ni = 0; ni < N; ni++) {
				result.m_data[ni] = std::abs(m_data[ni]);
			}
			return result;
		}

		/**
		 * @brief Retrieve an iterator to the first lexicographic element of the field
		 * @return an iterator to this field at coordinate (0, 0, ...)
		 */
		iterator Begin() {
			return iterator(this);
		}


		/**
		 * @brief Adds a border around the field with the specified width and value. The width of the border is specified
		 * for each dimension, so that different sized borders can be implemented.
		 * @param width is the width of the border that will be added in each dimension
		 * @param val is the value used to fill the border
		 * @return a new field with the specified border added
		 */
		field Border(const std::vector<size_t> width, T val=0) {
			std::vector<size_t> new_shape = m_shape;
			size_t D = m_shape.size();

			for (size_t si = 0; si < width.size(); si++) {
				new_shape[si] += 2 * width[si];				// calculate new shape
			}

			field<T> result(new_shape, val);				// initialize result field with val as outer values

			std::vector<size_t> new_coord(D);
			for (size_t ui = 0; ui < Size(); ui++) {
				Coord(ui, new_coord);						// find new coordinates
				for (int si = 0; si < width.size(); si++) {
					new_coord[si] += width[si];
				}
				size_t new_index = result.Idx(new_coord);	//find new index
				result.m_data[new_index] = m_data[ui];		// replace the values from the initial field
			}
			return result;
		}

		/**
		 * @brief Adds a border around the entire field with the specified width and value.
		 * @param width is the width of the border that will be added in each dimension
		 * @param val is the value used to fill the border
		 * @return a new field with the specified border added
		 */
		field Border(const size_t width, T val = 0) {
			std::vector dist_vect = { width, width };
			return Border(dist_vect, val);
		}

		/**
		 * @brief Adds a border that replicates the value of the nearest point on the edge
		 * @param width is the size of the border along each dimension
		 * @return a new field with the replicated border added
		 */
		field BorderReplicate(std::vector<size_t> width) {
			std::vector<size_t> new_shape = m_shape;
			size_t D = m_shape.size();
			std::vector<size_t> coord_i(D);					// initialize memory for the coordinates in the cycle

			for (size_t si = 0; si < D; si++) {
				new_shape[si] += 2 * width[si];				// calculate new shape
			}

			field result(new_shape);						// initialize result field

			size_t total_size = result.Size();
			for (size_t index = 0; index < total_size; index++) {				// for each element in resulting field

				result.Coord(index, coord_i);									// convert index to coordinates
				for (size_t di = 0; di < D; di++) {
					coord_i[di] = std::max(width[di], coord_i[di]);				// bound in range [dist; n + 2 * dist)
					coord_i[di] -= width[di];									// move the range to [0; n + dist)
					coord_i[di] = std::min(m_shape[di] - 1, coord_i[di]);		// bound in range [0; n)
				}
				result.m_data[index] = m_data[Idx(coord_i)];
			}
			return result;
		}

		/**
		 * @brief Adds a border that replicates the value of the nearest point on the edge
		 * @param width is the size of the border that will be added around the entire field
		 * @return
		 */
		field<T> BorderReplicate(size_t width) {
			std::vector dist_vect(m_shape.size(), width);
			return BorderReplicate(dist_vect);
		}

		/**
		 * @brief Returns the number of bytes required to store the field data.
		 * @return the number of bytes storing the field data
		 */
		size_t Bytes() const {
			return Size() * sizeof(T);
		}

		/**
		 * @brief Returns a const pointer to the raw sample data so that it can be used in const functions for read-only access
		 * @return a const pointer to the raw sample data
		 */
		const T* ConstData() const {
			return &m_data[0];
		}




		/**
		 * @brief Convolves the current field with a kernel field, returning only the valid "center" portion of the result.
		 * @tparam C is the data type for the kernel field. This type should be compatible with T for any operations involved in
		 * the convolution (specifically addition and multiplication).
		 * @param K is a field representing the kernel that will be convolved with the current field
		 * @return a new field representing the convolution of the current field with K. Note that only the valid central part
		 * of the field will be returned (the result will be smaller than the current field).
		 */
		template<class C>
		field Convolve(const field<C> K) const {

			// make sure the number of dimensions in k matches the number of dimensions in the current field
			if (K.shape().size() > m_shape.size()) {
				std::cout << "error: convolution kernel (d = " << K.shape().size() << ") has more dimensions than the field (d = " << m_shape.size() << ")" << std::endl;
				exit(1);
			}

			// calculate the size of the output image and allocate space
			std::vector<size_t> result_shape = m_shape;												// start with the size of the source image
			for (size_t i = 0; i < result_shape.size(); i++) result_shape[i] -= K.shape()[i] - 1;	// offset the size by the size of the kernel (K-1)
			field result(result_shape);															// generate a new field to store the convolved image

			std::vector<size_t> k_coord(m_shape.size());
			std::vector<size_t> s_coord(m_shape.size());
			std::vector<size_t> r_coord(m_shape.size());
			T sum;																	// create a T value to accumulate the convolution result for a single pixel
			// FOR each pixel in the output (result) image
			for (size_t ri = 0; ri < result.Size(); ri++) {							// for each element in the result field
				sum = (T)0;
				result.Coord(ri, r_coord);											// get the coordinate in the result array
				// FOR each pixel in the kernel
				for (size_t ki = 0; ki < K.size(); ki++) {							// for each element in the kernel
					K.coord(ki, k_coord);											// get the coordinate in the kernel

					// CALCULATE the coordinate into the source image
					for (size_t ci = 0; ci < s_coord.size(); ci++) {				// for each coordinate
						s_coord[ci] = r_coord[ci] + k_coord[ci];					// offset the coordinate by the kernel coordinate
					}
					size_t s_idx = Idx(s_coord);									// calculate the index into the source image
					T s_val = m_data[s_idx];											// get the value from the source field
					C kval = K.read(k_coord);										// get the value of the kernel
					sum = sum + kval * s_val;										// accumulate the product of the kernel and the source
				}
				result(r_coord) = sum;												// store the accumulated value in the result field
			}

			return result;															// return the result field
		}


		/**
		 * @brief Converts a one-dimensional index into the data vector to an N-dimensional coordinate to the sample point
		 * @param idx is the one-dimensional index into linear memory for the field
		 * @param c is a coordinate that will be written to
		 */
		void Coord(size_t idx, std::vector<size_t> &c) const {

			size_t nC = c.size();
			size_t step;
			size_t i, j;
			for (i = 0; i < nC; i++) {
				step = 1;
				for (j = i + 1; j < nC; j++) {
					step *= m_shape[j];
				}
				c[i] = idx / step;
				idx -= c[i] * step;
			}
		}


		/**
		 * @brief Crop a contiguous region from the field and returns it as a new field
		 * @param min_coord is the coordinate of the minimum corner of the cropped region
		 * @param max_coord is the coordinate of the maximum corner of the cropped region
		 * @return
		 */
		field Crop(std::vector<size_t> min_coord, std::vector<size_t> max_coord) const {

			size_t D = m_shape.size();												// get the number of dimensions
			std::vector<size_t> new_dims(D);
			for (size_t di = 0; di < D; di++) {										// calculate the size of the output field
				new_dims[di] = max_coord[di] - min_coord[di];
			}

			field<T> result(new_dims);
			std::vector<size_t> s_coord(m_shape.size());								// stores the source field coordinate
			std::vector<size_t> r_coord(m_shape.size());								// stores the result field coordinate
			// FOR each pixel in the output (result) image
			for (size_t ri = 0; ri < result.Size(); ri++) {							// for each element in the result field
				result.Coord(ri, r_coord);											// get the coordinate in the result array

				// CALCULATE the coordinate into the source image
				for (size_t ci = 0; ci < s_coord.size(); ci++) {				// for each coordinate
					s_coord[ci] = r_coord[ci] + min_coord[ci];					// offset the coordinate by the kernel coordinate
				}
				size_t s_idx = Idx(s_coord);									// calculate the index into the source field
				T s_val = m_data[s_idx];											// get the value from the source field

				result(r_coord) = s_val;											// store the accumulated value in the result field
			}

			return result;															// return the result field
		}

		/**
		 * @brief Returns a pointer to the raw sample data
		 * @return a pointer to the raw sample data
		 */
		T* Data() {
			return &m_data[0];
		}

		/**
		 * @brief Calculates the partial derivative of the field with respect to the specified dimension
		 * @param axis is the dimension along which the partial derivative is computed
		 * @param d is the derivative to calculate (ex. 2 is the second derivative)
		 * @param order is the order of accuracy to use in the calculation. This order will be used for all sample points,
		 * including those at the edges
		 * @param print_coefs is a flag that prints the finite difference coefficients to the console for debugging
		 * @return a new field representing the partial derivative of the current field
		 */
		field Derivative(const unsigned int axis, const unsigned int d, const unsigned int order, const bool print_coefs = false) {

			T* ptr = m_Derivative(axis, d, order, print_coefs);

			field result(m_shape, ptr);

			delete ptr;
			return result;
		}

		/**
		 * @brief Retrieve an iterator just past the last element of the field
		 * @return an iterator to this field pointing just past the last element
		 */
		iterator End() {
			return iterator(this, Data() + Size());
		}

		/**
		 * @brief Converts the specified coordinate into a one-dimensional index that can be used to directly access elements in the m_data vector
		 * @param coord is the N-dimensional coordinate that will be converted into an index
		 * @return an index to the element specified by the coord vector
		 */
		size_t Idx(const std::vector<size_t> &coord) const {
			size_t i = 0;

			size_t step;
			for (size_t c = 0; c < coord.size(); c++) {
				step = 1;
				for (size_t j = c + 1; j < coord.size(); j++) {
					step *= m_shape[j];
				}
				i += coord[c] * step;
			}
			return i;
		}


		/**
		 * @brief Loads a NumPy file and stores the resulting array in the field.
		 *
		 * Note that NumPy files only support basic data types, so this
		 * doesn't work if the sample points represent complex data structures. However, the D typename can be used to load data structures that are
		 * composed of a single several basic data types. For example, a field where each element T is a 3x3 matrix consisting of 32-bit floating
		 * point numbers can be loaded by specifying D as float.
		 *
		 * This function uses an externaly library cnpy (https://github.com/rogersce/cnpy) to load the file.
		 *
		 * @tparam D is the data type used to store values in the NumPy file
		 * @param filename is the name of the NumPy file to load
		 */
		template<typename D = T>
		void LoadNpy(const std::string filename) {					// fill the field with data from an npy file
			std::vector<unsigned long> shape;
			std::vector<D> data;
			bool fortran_order;
			try {
				npy::LoadArrayFromNumpy<D>(filename, shape, fortran_order, data);	// load NPY array and metadata
			}
			catch(const std::runtime_error &e){
				std::cout << "ERROR loading NumPy file "<<filename<<": " << e.what() << std::endl;
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
				m_SetShape(new_shape);
			}
			else {
				m_SetShape(new_shape);													// set the new shape of the field
			}
			m_Allocate();																	// allocate space for the field

			memcpy(reinterpret_cast<void*>(&m_data[0]), reinterpret_cast<void*>(&data[0]), Bytes());									// copy the data from the
		}

		/**
		 * @brief Calculates the average field value
		 * @return a scalar value representing the mean field value
		 */
		T Mean() {
			T s = Sum();
			T m = s / Size();
			return m;
		}

		/**
		 * @brief Returns the number of dimensions in the field
		 * @return the number of dimensions
		 */
		unsigned int Ndims() const { return m_shape.size(); }

		/**
		 * @brief Returns a copy of the element referenced by the given coordinate x. Use this function if you do not want
		 * to change or update the specified value
		 * @param x is the coordinate of the element to be copied
		 * @return a non-referenced copy of the element
		 */
		T Read(const std::vector<size_t> &x) const {
			size_t i = Idx(x);
			if (i >= m_data.size()) {
				throw "ERROR: index out of range:";
			}
			return m_data[i];
		}

		/**
		 * @brief Returns a copy of the element referenced by the given coordinate x. Use this function if you do not want
		 * to change or update the specified value
		 * @param i is the one-dimensional index of the element to be copied
		 * @return a non-referenced copy of the element
		 */
		T Read(const size_t i) const {
			if (i >= m_data.size()) {
				throw "ERROR: index out of range";
			}
			return m_data[i];
		}

		/**
		 * @brief Reshapes the current field by changing the m_shape attribute. This doesn't change the data in the field, and the new shape must be able
		 * to accomodate the data of the current field.
		 * @param new_shape is the desired shape for the field
		 */
		void Reshape(std::vector<size_t> new_shape) {
			size_t n = new_shape[0];
			for (size_t i = 1; i < new_shape.size(); i++)
				n *= new_shape[i];
			if (n != Size()) {
				std::cout << "ERROR: new shape (";
				for (size_t i = 0; i < new_shape.size(); i++)
					std::cout << new_shape[i];
				std::cout << ") and old shape (";
				for (size_t i = 0; i < m_shape.size(); i++)
					std::cout << m_shape[i];
				std::cout << ") must have the same number of elements." << std::endl;
			}
			m_shape = new_shape;
		}

		/**
		 * @brief Change the size and shape of the field and reallocate the data
		 * @param new_size is a vector specifying the new shape of the field
		 */
		void Resize(const std::vector<size_t> new_size) {
			m_SetShape(new_size);
		}

		/**
		 * @brief Saves a NumPy file and stores the resulting array in the field.
		 *
		 * Note that NumPy files only support basic data types, so this doesn't work if the sample points represent complex data structures.
		 * However, the D typename can be used to load data structures that are composed of a single several basic data types. For example,
		 * a field where each element T is a 3x3 matrix consisting of 32-bit floating point numbers can be saved by specifying D as float.
		 *
		 * This function uses an externaly library cnpy (https://github.com/rogersce/cnpy) to load the file.
		 *
		 * @tparam D is the data type used to store values in the NumPy file
		 * @param filename is the name of the NumPy file to save
		 */
		template<typename D = T>
		void SaveNpy(const std::string& filename) {
			std::vector<size_t> destshape = m_shape;
			SaveNpy<D>(filename, destshape);
		}

		/**
		 * @brief Saves the current field data to a NumPy file, allowing the user to specify the dimensions of the saved grid.
		 * @tparam D is the basic data type used for NumPy (see previous save/load functions)
		 * @param filename is the file name for the NumPy file that will be produced
		 * @param dest_shape is the shape of the output NumPy grid
		 */
		template<typename D = T>
		void SaveNpy(const std::string& filename, std::vector<size_t> dest_shape) {
			bool fortran_order = false;										// default to standard (C) order
			std::vector cast_dest_shape(dest_shape.begin(), dest_shape.end());

			if (sizeof(D) < sizeof(T))										// if the data type stored is smaller than the grid data type
				cast_dest_shape.push_back(sizeof(T) / sizeof(D));			// add another dimension to account for this
			try {
				npy::SaveArrayAsNumpy(filename, fortran_order, cast_dest_shape.size(), (const unsigned long*)&cast_dest_shape[0], (D*)(&m_data[0]));
			}
			catch (std::out_of_range&) {
				throw std::runtime_error("field ERROR: data type not recognized by external code npy.hpp");
			}
		}

		/**
		 * @brief Saves the field as a direct copy to disk
		 * @param filename is the name of the RAW file that will be produced
		 */
		void SaveRaw(const std::string& filename) {
			std::ofstream outfile;
			outfile.open(filename, std::ios::app | std::ios::binary);
			outfile.write((char*)&m_data[0], Bytes());
		}




		/**
		 * @brief Returns the current shape of the field
		 * @return the shape of the field as an N-dimensional std::vector
		 */
		std::vector<size_t> Shape() const { return m_shape; }


		/**
		 * @brief Returns the result of the sgn() function for each sample point as a field:
		 *          {	-1		x < 0
		 * sgn(x) = {	0		x = 0
		 *			{	1		x > 0
		 * @return a new field with the result of sgn(x) at each sample point
		 */
		field Sign() {
			field R(m_shape);
			for (size_t xi = 0; xi < Size(); xi++) {
				if (m_data[xi] < 0) R.m_data[xi] = -1;
				else if (m_data[xi] > 0) R.m_data[xi] = 1;
				else R.m_data[xi] = 0;
			}
			return R;
		}


		/**
		 * @brief Returns the total number of sample points in the field across all dimensions
		 * @return the total number of samples (the size of the m_data vector)
		 */
		size_t Size() const {
			
			return m_data.size();
		}







		

		/*template<class C>
		field<T> convolve2(const field<C> K) const {

			// make sure all of the dimensions are correct
			if (m_shape.size() != 2) {
				std::cout << "ERROR: cannot perform a 2D convolution on a field that is not 2D (D = " << m_shape.size() << ")" << std::endl;
				exit(1);
			}
			if (K.shape().size() != 2) {
				std::cout << "ERROR: cannot perform a 2D convolution with a kernel that is not 2D (D = " << K.shape().size() << ")" << std::endl;
				exit(1);
			}

			// calculate the size of the output image (only valid region is returned) and allocate the result field
			std::vector<size_t> result_shape = m_shape;
			size_t RY = m_shape[0] - (K.shape()[0] - 1);
			size_t RX = m_shape[1] - (K.shape()[1] - 1);
			std::vector<size_t> R = { RY, RX };
			field<T> result(R);

			// store the dimensions of the kernel and source image
			size_t KY = K.shape()[0];
			size_t KX = K.shape()[1];
			size_t SY = m_shape[0];
			size_t SX = m_shape[1];

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
							input_v = Read(source_i);					// read the value of the source pixel
							result_v += input_v *k_v;					// accumulate the result
						}
					}
					result.m_data[result_i] = result_v;					// store the result
				}
			}
			return result;

		}
		*/















		field<T> gradmag(unsigned int order) {

			T* sum_sq = new T[Size()]{};							// allocate space to store the sum of squares of each derivative (initialize to zero)

			for (size_t di = 0; di < m_shape.size(); di++) {			// for each dimension
				if (m_shape[di] == 1) continue;						// if the dimension is singular, ignore it
				T* part_deriv = m_Derivative(di, 1, order);		// calculate the partial derivative along dimension di
				for (size_t xi = 0; xi < Size(); xi++) {			// for each point in the field
					sum_sq[xi] += part_deriv[xi] * part_deriv[xi];	// calculate the square and sum
				}
			}
			
			// calculate the square root for each element
			for (size_t xi = 0; xi < Size(); xi++) {			// for each point in the field
				sum_sq[xi] = std::sqrt(sum_sq[xi]);				// calculate the square root
			}
			field<T> result(m_shape, sum_sq);
			delete sum_sq;
			return result;
		}



		/**
		 * @brief Sum of all of the elements in the field
		 * @return a scalar value representing the sum of all sample points
		 */
		T Sum() {
			T s = 0;
			for (size_t xi = 0; xi < Size(); xi++) {
				s += m_data[xi];
			}
			return s;
		}



		//static std::vector<double> fd_coefficients(unsigned int d, unsigned int order) {
		//	return tira::calculus::central_difference_coefficients(d, order);
		//}

		/*field<T> central_derivative(unsigned int axis, unsigned int d, unsigned int order, bool pad = true, bool print_coefs = false) {

			std::vector<double> C = tira::calculus::central_difference_coefficients(d, order);
			if(print_coefs) tira::calculus::printCoefficients(C);

			std::vector<size_t> kernel_shape(m_shape.size(), 1);
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
				D = BorderReplicate(dist).convolve(kernel);
			}
			else {
				D = convolve(kernel);
			}
			return D;
		}*/

		
	};
}