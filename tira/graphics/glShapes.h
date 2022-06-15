#pragma once

#include "shapes.h"

namespace tira{

	template <typename T>
	glGeometry glGenerateRectangle() {
		rectangle<T> r;
		std::vector<T> v = r.getInterleavedVertices();									// generate a vector of interleaved vertices
		std::vector<unsigned int> i = r.getIndices();									// generate a vector of indices
		glGeometry rectangle;

		glVertexBufferLayout layout;
		layout.Push<T>(3);																// add three vertex components to the layout
		layout.Push<T>(3);																// add three normal components to the layout
		layout.Push<T>(2);																// add two texture components to the layout
		rectangle.AddMesh(&v[0], v.size() * sizeof(T), layout, &i[0], i.size(), 0);

		return rectangle;
	}

	template <typename T>
	glGeometry glGenerateCircle(unsigned int slices = 10) {
		circle<T> c(slices);
		std::vector<T> v = c.getInterleavedVertices();									// generate a vector of interleaved vertices
		std::vector<unsigned int> i = c.getIndices();									// generate a vector of indices
		glGeometry circle;

		glVertexBufferLayout layout;
		layout.Push<T>(3);																// add three vertex components to the layout
		layout.Push<T>(3);																// add three normal components to the layout
		layout.Push<T>(2);																// add two texture components to the layout
		circle.AddMesh(&v[0], v.size() * sizeof(T), layout, &i[0], i.size(), 0);

		return circle;
	}

    template <typename T>
	glGeometry glGenerateCube(){
		cube<T> c;
		std::vector<T> v = c.getInterleavedVertices();
		std::vector<unsigned int> i = c.getIndices();
		glGeometry cube;

		glVertexBufferLayout layout;
		layout.Push<T>(3);
		layout.Push<T>(3);
		layout.Push<T>(2);
		cube.AddMesh(&v[0], v.size() * sizeof(T), layout, &i[0], i.size(), 0);

		return cube;
	}

	template <typename T>
	glGeometry glGenerateSphere(unsigned int stacks, unsigned int slices) {
		sphere<T> s(stacks, slices);
		std::vector<T> v = s.getInterleavedVertices();
		std::vector<unsigned int> i = s.getIndices();

		glGeometry sphere;

		glVertexBufferLayout layout;
		layout.Push<T>(3);
		layout.Push<T>(3);
		layout.Push<T>(2);
		sphere.AddMesh(&v[0], v.size() * sizeof(T), layout, &i[0], i.size(), 0);

		return sphere;
	}

	template <typename T>
	glGeometry glGenerateIcosahedron() {
		tira::icosahedron<T> ico(0.5f);
		glGeometry icosahedron;

		glVertexBufferLayout layout;
		layout.Push<T>(3);
		layout.Push<T>(3);
		layout.Push<T>(2);
		icosahedron.AddMesh(ico.getInterleavedVertices(), ico.getInterleavedVertexCount() * 8 * sizeof(T), layout, ico.getIndices(), ico.getIndexCount(), 0);

		return icosahedron;
	}

	template <typename T>
	glGeometry glGenerateIcosphere(unsigned int subdiv, bool smooth = false) {
		tira::icosphere ico(0.5f, subdiv, smooth);
		glGeometry icosphere;

		glVertexBufferLayout layout;
		layout.Push<T>(3);
		layout.Push<T>(3);
		layout.Push<T>(2);
		icosphere.AddMesh(ico.getInterleavedVertices(), ico.getInterleavedVertexCount() * 8 * sizeof(T), layout, ico.getIndices(), ico.getIndexCount(), 0);

		return icosphere;
	}
}