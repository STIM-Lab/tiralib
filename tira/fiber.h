#pragma once

#include <vector>

#include <glm/glm.hpp>

namespace tira {

	/**
	 * @brief      Class describes a vertex consisting of a 3D coordinate and other user-defined attributes.
	 *
	 * @tparam     VertexAttribute  Data type for user-defined attributes tstored with each vertex position
	 */
	template <typename VertexAttribute = float>
	class vertex : public glm::vec3 {

		/**
		 * User-defined attribute stored at each vertex position (ex. radius, color, etc)
		 */
		VertexAttribute _va;
	public:

		/**
		 * @brief      Constructs a new vertex from a 3D coordinate and attribute
		 *
		 * @param[in]  p     { parameter_description }
		 * @param[in]  r     { parameter_description }
		 */
		vertex(glm::vec3 p, VertexAttribute r) : glm::vec3(p) { _va = r; }

		void va(VertexAttribute r) { _va = r; }
		VertexAttribute va() { return _va; }
	};

	/**
	 * @brief      Class defines a fiber as an array of sequential vertices.
	 *
	 * @tparam     VertexAttribute  Data type for additional attributes stored with each vertex position
	 */
	template <typename VertexAttribute = float>
	class fiber : public std::vector< vertex<VertexAttribute> > {


	protected:
		std::vector< float > _lv;		// array storing the parameterized length values at each vertex

		float _gaussian(float d, float sigma) {
			float n = 1.0f / std::sqrt(2.0f * std::numbers::pi * sigma * sigma);
			float y = -(d * d) / (2 * sigma * sigma);
			return n * std::exp(y);
		}

	public:

		fiber() : std::vector< vertex<VertexAttribute> >() {}

		/**
		 * @brief      Creates a vertex from a position and attribute, and inserts it at the end of the fiber
		 *
		 * @param[in]  p     3D coordinate providing the spatial position of the vertex
		 * @param[in]  r     user-defined attribute associated with this vertex (ex. radius)
		 */
		void push_back(glm::vec3 p, VertexAttribute r) {
			vertex<VertexAttribute> v(p, r);						// generate a new vertex from the provided position and attribute

			// update the LV vector to store the length at the new vertex
			if (this->size() == 0) _lv.push_back(0);				// if the current fiber is empty, the first vertex is l(v) = 0
			else {
				float d = glm::length(p - this->back()) ;			// otherwise l(v) = l(v_{n-1}) + |v_{n} - v_{n-1}|
				_lv.push_back(d + _lv.back());
			}

			std::vector< vertex<VertexAttribute> >::push_back(v);	// push the vertex into the fiber
		}

		void push_back(vertex<VertexAttribute> v) {
			push_back(glm::vec3(v), v.va());
		}

		fiber smooth_gaussian(float sigma, bool anchor_endpoints = true) {

			fiber smoothed;
			for (size_t vi=0; vi<this->size(); vi++) {
				if (anchor_endpoints && (vi == 0 || vi == this->size() - 1)) {
					vertex v = this->at(vi);
					smoothed.push_back(v);
				}
				else {
					glm::vec3 p(0.0f);
					float g_energy = 0.0f;
					for (size_t pi=0; pi<this->size(); pi++) {
						float d = _lv[vi] - _lv[pi];
						p += _gaussian(d, sigma) * this->at(pi);
						g_energy += _gaussian(d, sigma);
					}
					smoothed.push_back(p / g_energy, this->at(vi).va());
				}
			}
			return smoothed;
		}

	};
}
