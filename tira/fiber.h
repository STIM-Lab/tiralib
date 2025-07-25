#pragma once

#include <vector>

#include <glm/glm.hpp>

namespace tira {

	/**
	 * @brief      Class describes a vertex consisting of a 3D coordinate and other user-defined attributes.
	 *
	 * @tparam     VertexAttribute  Data type for user-defined attributes tstored with each vertex position
	 */
	template <typename VertexAttribute>
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

	public:

		fiber() : std::vector< vertex<VertexAttribute> >() {}

		/**
		 * @brief      Calculates the length of the fiber using linear interpolation between points.
		 *
		 * @return     Current fiber length
		 */
		float length() {

			float l = 0.0f;
			for (size_t vi = 1; vi < this->size(); vi++) {
				l += glm::length(this->at(vi) - this->at(vi + 1));
			}

			return l;
		}

		/**
		 * @brief      Creates a vertex from a position and attribute, and inserts it at the end of the fiber
		 *
		 * @param[in]  p     3D coordinate providing the spatial position of the vertex
		 * @param[in]  r     user-defined attribute associated with this vertex (ex. radius)
		 */
		void push_back(glm::vec3 p, VertexAttribute r) {
			vertex<VertexAttribute> v(p, r);
			std::vector< vertex<VertexAttribute> >::push_back(v);
		}

	};
}
