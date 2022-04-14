#include "shapes3d.h"

namespace tira{
    template <typename T>
	glGeometry glGenerateCube(){
		std::vector<T> v;
		std::vector<unsigned int> i;
		GenerateCube(v, i);
		glGeometry cube;

		glVertexBufferLayout layout;
		layout.Push<float>(3);
		layout.Push<float>(3);
		layout.Push<float>(2);
		cube.AddMesh(&v[0], v.size() * sizeof(T), layout, &i[0], i.size(), 0);

		return cube;
	}

	template <typename T>
	glGeometry glGenerateSphere(unsigned int stacks, unsigned int slices) {
		std::vector<T> v;
		std::vector<unsigned int> i;
		GenerateSphere(v, i, stacks, slices);
		glGeometry sphere;

		glVertexBufferLayout layout;
		layout.Push<float>(3);
		layout.Push<float>(3);
		layout.Push<float>(2);
		sphere.AddMesh(&v[0], v.size() * sizeof(T), layout, &i[0], i.size(), 0);

		return sphere;
	}
}