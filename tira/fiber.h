#pragma once

#include <algorithm>
#include <vector>
#include <array>
#include <iostream>
#include <numbers>

#include <glm/glm.hpp>

namespace tira {

	/**
	 * @brief      Class describes a vertex consisting of a 3D coordinate and other user-defined attributes.
	 *
	 * @tparam     VertexAttribute  Data type for user-defined attributes tstored with each vertex position
	 */
	template <typename VertexAttributeType = float>
	class vertex : public glm::vec3 {

		/**
		 * User-defined attribute stored at each vertex position (ex. radius, color, etc)
		 */
		VertexAttributeType m_va;
	public:

		vertex() {}

		/**
		 * @brief      Constructs a new vertex from a 3D coordinate and attribute
		 *
		 * @param[in]  p     { parameter_description }
		 */
		vertex(const glm::vec3& p) : glm::vec3(p) { m_va = 0; }

		/**
		 * @brief      Constructs a new vertex from a 3D coordinate and attribute
		 *
		 * @param[in]  p     { parameter_description }
		 * @param[in]  r     { parameter_description }
		 */
		vertex(glm::vec3 p, VertexAttributeType r) : glm::vec3(p) { m_va = r; }

		vertex(const vertex& v) : glm::vec3(v) { m_va = v.m_va; }

		vertex& operator=(const vertex& v) {
			glm::vec3::operator=(v);
			m_va = v.m_va;
			return *this;
		}


		void Attribute(VertexAttributeType r) { m_va = r; }
		VertexAttributeType Attribute() const { return m_va; }
	};

	/**
	 * @brief      Class defines a fiber as an array of sequential vertices.
	 *
	 * @tparam     VertexAttribute  Data type for additional attributes stored with each vertex position
	 */
	template <typename VertexAttributeType = float>
	class fiber : public std::vector< vertex<VertexAttributeType> > {


	protected:
		std::vector< float > m_lparam;		// array storing the parameterized length values at each vertex

		float m_Gaussian(float d, float sigma) {
			float n = 1.0f / std::sqrt(2.0f * std::numbers::pi * sigma * sigma);
			float y = -(d * d) / (2 * sigma * sigma);
			return n * std::exp(y);
		}

		fiber<float> m_Derivative() {
			fiber<float> d_dt;

			glm::vec3 p0 = this->at(0);
			glm::vec3 p1 = this->at(1);
			float l0_1 = glm::length(p1 - p0);

			if (l0_1 != 0)
				d_dt.AddLastVertex(vertex<float>((p1 - p0) / l0_1, l0_1));
			else
				d_dt.AddLastVertex(vertex<float>(glm::vec3(0.0f), 0));

			float l1_2, l0_2;
			glm::vec3 p2;
			for (size_t pi = 1; pi < this->size() - 1; pi++) {
				p0 = this->at(pi - 1);
				p1 = this->at(pi);
				p2 = this->at(pi + 1);

				l0_1 = glm::length(p1 - p0);
				l1_2 = glm::length(p2 - p1);
				l0_2 = l0_1 + l1_2;
				if (l0_2 != 0)
					d_dt.AddLastVertex((p2 - p0) / l0_2, l0_2);
				else
					d_dt.AddLastVertex(vertex<float>(glm::vec3(0.0f), 0));

			}

			p1 = this->at(this->size() - 2);
			p2 = this->at(this->size() - 1);
			l1_2 = glm::length(p2 - p1);
			if (l1_2 != 0)
				d_dt.AddLastVertex(vertex<float>((p2 - p1) / l1_2, l1_2));
			else
				d_dt.AddLastVertex(vertex<float>(glm::vec3(0.0f), 0));

			return d_dt;
		}

		/*bool m_TestDuplicates(vertex<VertexAttribute> a, vertex<VertexAttribute> b) {
			if (a.x == b.x && a.y == b.y && a.z == b.z)
				return true;
			return false;
		}*/

	public:

		fiber() : std::vector< vertex<VertexAttributeType> >() {}

		/**
		 * @brief      Creates a vertex from a position and attribute, and inserts it at the end of the fiber
		 *
		 * @param[in]  p     3D coordinate providing the spatial position of the vertex
		 * @param[in]  r     user-defined attribute associated with this vertex (ex. radius)
		 */
		void AddLastVertex(glm::vec3 p, VertexAttributeType r) {
			vertex<VertexAttributeType> v(p, r);						// generate a new vertex from the provided position and attribute

			// update the LV vector to store the length at the new vertex
			if (this->size() == 0) m_lparam.push_back(0);				// if the current fiber is empty, the first vertex is l(v) = 0
			else {
				float d = glm::length(p - this->back());			// otherwise l(v) = l(v_{n-1}) + |v_{n} - v_{n-1}|
				m_lparam.push_back(d + m_lparam.back());
			}

			std::vector< vertex<VertexAttributeType> >::push_back(v);	// push the vertex into the fiber
		}

		void AddLastVertex(vertex<VertexAttributeType> v) {
			AddLastVertex(glm::vec3(v), v.Attribute());
		}

		bool BoundingBox(glm::vec3& aabb_min, glm::vec3& aabb_max) {
			if (this->size() == 0) return false;						// if the fiber is empty, return false (AABB is invalid)

			aabb_min = this->at(0);									// initialize the corners with the first point
			aabb_max = aabb_min;
			for (size_t pi = 1; pi < this->size(); pi++) {				// for each of the remaining points in the fiber
				glm::vec3 p = this->at(pi);
				if (p.x < aabb_min.x) aabb_min.x = p.x;
				else if (p.x > aabb_max.x) aabb_max.x = p.x;

				if (p.y < aabb_min.y) aabb_min.y = p.y;
				else if (p.y > aabb_max.y) aabb_max.y = p.y;

				if (p.z < aabb_min.z) aabb_min.z = p.z;
				else if (p.z > aabb_max.z) aabb_max.z = p.z;
			}
			return true;
		}

		float Length() {
			float l = 0;
			for (size_t vi = 1; vi < this->size(); vi++) {
				l += glm::length(this->at(vi) - this->at(vi - 1));
			}
			return l;
		}

		float ChordLength() {
			if (this->size() < 2) return 0.0f;  // not enough points to compute a chord
			return glm::length(this->back() - this->front());
		}

		/**
		 * @brief Smooth a fiber by applying a Gaussian kernel to its vertex coordinates
		 * @param sigma is the standard deviation for the smoothing kernel
		 * @return a smoothed version of the current fiber
		 */
		fiber Smooth(float sigma) {

			fiber smoothed_fiber;

			// add the first vertex to the smoothed fiber unchanged - the end points will remain the same
			smoothed_fiber.AddLastVertex(this->at(0));

			// pre-calculate the Gaussian window size and number of vertices
			int window = (int)(sigma * 3);
			int num_vertices = (int)this->size();

			// for each internal vertex of the fiber
			for (int vi = 1; vi < num_vertices - 1; vi++) {
				float gaussian_integral = 0.0f;
				glm::vec3 weighted_point(0.0f);

				// for every point within range of the current point to be smoothed
				for (int wi = -window; wi < window; wi++) {
					int fid = vi + wi;							// calculate the ID of the windowed vertex

					// if the windowed vertex is within the fiber, apply the appropriate weights
					if (fid >= 0 && fid < num_vertices) {
						float d = std::abs(m_lparam[vi] - m_lparam[fid]);
						float gaussian_value = m_Gaussian(d, sigma);
						gaussian_integral += gaussian_value;
						weighted_point += gaussian_value * this->at(fid);
					}
				}
				glm::vec3 smoothed_point = weighted_point / gaussian_integral;
				smoothed_fiber.AddLastVertex(smoothed_point, this->at(vi).Attribute());
			}

			// add the last point to the smoothed fiber unchanged
			smoothed_fiber.AddLastVertex((this->at(0)));

			// return the smoothed version of the fiber
			return smoothed_fiber;

		}

		size_t RemoveDuplicates() {
			size_t size_before = this->size();

			auto last = std::unique(this->begin(), this->end());
			this->erase(last, this->end());

			size_t size_after = this->size();

			return size_before - size_after;
		}

		std::vector<glm::vec3> Derivative() {
			fiber<float> d_dt = m_Derivative();

			std::vector<glm::vec3> result;
			for (size_t vi = 0; vi < this->size(); vi++) {
				result.push_back(d_dt[vi]);
			}
			return result;
		}

		/**
		 * The curvature is calculated using partial derivatives:
		 * \f[
		 * \kappa = \frac{\sqrt{(z''y'-y''z')^2 + (x''z'-z''x')^2 + (y''x' - x''y')^2}}{\sqrt{(x'^2+y'^2+z'^2)^3}}
		 * \f]
		 * this is calculated by evaluating the following terms:
		 * \f[
		 * \kappa = \frac{\sqrt{ a^2 + b^2 + c^2 }}{\sqrt{d^3}}
		 * \f]
		 * @return
		 */
		std::vector<float> Curvature() {
			std::vector<float> kappa;
			fiber<float> d_dt = m_Derivative();
			fiber<float> dd_dt = d_dt.m_Derivative();
			for (size_t vi = 0; vi < this->size(); vi++) {
				float xp = d_dt[vi].x;
				float yp = d_dt[vi].y;
				float zp = d_dt[vi].z;

				float xpp = dd_dt[vi].x;
				float ypp = dd_dt[vi].y;
				float zpp = dd_dt[vi].z;

				float a = zpp * yp - ypp * zp;
				float b = xpp * zp - zpp * xp;
				float c = ypp * xp - xpp * yp;
				float d = xp * xp + yp * yp + zp * zp;

				float num = std::sqrt(a * a + b * b + c * c);
				float denom = std::sqrt(d * d * d);

				kappa.push_back(num / denom);
			}
			return kappa;
		}



	};
}