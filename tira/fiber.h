#pragma once

#include <algorithm>
#include <vector>
#include <array>
#include <iostream>

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

		fiber<float> _d1() {
			fiber<float> d_dt;

			glm::vec3 p0 = this->at(0);
			glm::vec3 p1 = this->at(1);
			float l0_1 = glm::length(p1 - p0);

			if (l0_1 != 0)
				d_dt.push_back(vertex<float>((p1 - p0) / l0_1, l0_1));
			else
				d_dt.push_back(vertex<float>(glm::vec3(0.0f), 0));

			float l1_2, l0_2;
			glm::vec3 p2;
			for (size_t pi=1; pi<this->size()-1; pi++) {
				p0 = this->at(pi-1);
				p1 = this->at(pi);
				p2 = this->at(pi+1);

				l0_1 = glm::length(p1 - p0);
				l1_2 = glm::length(p2 - p1);
				l0_2 = l0_1 + l1_2;
				if (l0_2 != 0)
					d_dt.push_back((p2 - p0) / l0_2, l0_2);
				else
					d_dt.push_back(vertex<float>(glm::vec3(0.0f), 0));

			}

			p1 = this->at(this->size()-2);
			p2 = this->at(this->size()-1);
			l1_2 = glm::length(p2 - p1);
			if (l1_2 != 0)
				d_dt.push_back(vertex<float>((p2 - p1) / l1_2, l1_2));
			else
				d_dt.push_back(vertex<float>(glm::vec3(0.0f), 0));

			return d_dt;
		}

		bool _test_duplicates(vertex<VertexAttribute> a, vertex<VertexAttribute> b) {
			if (a.x == b.x && a.y == b.y && a.z == b.z)
				return true;
			return false;
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

		float length() {
			float l = 0;
			for (size_t vi=1; vi<this->size(); vi++) {
				l+= glm::length(this->at(vi) - this->at(vi-1));
			}
			return l;
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

		size_t remove_duplicates() {
			size_t size_before = this->size();

			auto last = std::unique(this->begin(), this->end());
			this->erase(last, this->end());

			size_t size_after = this->size();

			return size_before - size_after;
		}

		std::vector<glm::vec3> derivative() {
			fiber<float> d_dt = _d1();

			std::vector<glm::vec3> result;
			for (size_t vi=0; vi<this->size(); vi++) {
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
		std::vector<float> curvature() {
			std::vector<float> kappa;
			fiber<float> d_dt = _d1();
			fiber<float> dd_dt = d_dt._d1();
			for (size_t vi=0; vi<this->size(); vi++) {
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

				float num = std::sqrt( a*a + b*b + c*c );
				float denom = std::sqrt(d * d * d);

				kappa.push_back( num / denom );
			}
			return kappa;
		}



	};
}
