#pragma once

#include <vector>

#include <glm/glm.hpp>

namespace tira {

	template <typename VertexAttribute>
	class vertex : public glm::vec3 {
	public:
		VertexAttribute a;

		vertex(glm::vec3 p, VertexAttribute r) : glm::vec3(p) {
			a = r;
		}
	};

	template <typename VertexAttribute = float>
	class fiber : public std::vector< vertex<VertexAttribute> > {


	protected:

		//std::vector< vertex<VertexAttribute> > _vertices;
		float _length;
		bool _length_valid;

	public:

		fiber() : std::vector< vertex<VertexAttribute> >() {
			_length = 0;
			_length_valid = false;
		}

		float length() {
			if (_length_valid) return _length;

			_length = 0.0f;
			for (size_t vi = 1; vi < this->size(); vi++) {
				_length += glm::length(this->at(vi) - this->at(vi + 1));
			}

			_length_valid = true;
			return _length;
		}

		void push_back(glm::vec3 p, VertexAttribute r) {
			push_back(p, r);
		}

	};
}
