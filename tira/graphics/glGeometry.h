# pragma once

#include "glErrorHandler.h"
#include "glVertexBuffer.h"
#include "glVertexBufferLayout.h"
#include "glIndexBuffer.h"
#include "glVertexArray.h"
#include "glShader.h"
#include "glTexture.h"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "shapes.h"

namespace tira {

	/// <summary>
	/// glGeometry class encodes the variables required to send vertex information to the OpenGL renderer. This
	/// includes a vertex array, vertex buffer, and index array.
	/// </summary>
	class glGeometry {
		glVertexBuffer vb;
		glVertexArray va;
		glIndexBuffer ib;

	public:
		glGeometry() {}

		void Unbind() {
			va.Unbind();
			vb.Unbind();
			ib.Unbind();
		}

		/// <summary>
		/// Add a vertex list and corresponding indices as a mesh to the geometry class.
		/// </summary>
		/// <param name="vertices">Vertex data stored in a CPU-based array</param>
		/// <param name="bytes">Size of the vertex data in bytes</param>
		/// <param name="layout"></param>
		/// <param name="indices">Array of indices specifying the order of vertices in the mesh</param>
		/// <param name="count">Number of indices</param>
		/// <param name="offset"></param>
		void AddMesh(const void* vertices, unsigned int bytes,
			glVertexBufferLayout layout,
			const unsigned int* indices, unsigned int count, unsigned int offset) {
			vb.SetBuffer(vertices, bytes);
			ib.SetBuffer(indices, count);
			va.AddBuffer(vb, layout, offset);
			//va.Bind();
			//ib.Bind();
		}

		void Draw() {
			va.Bind();              // Bind index array
			ib.Bind();              // Bind index buffer
			GLERROR(glDrawElements(GL_TRIANGLES, ib.GetCount(), GL_UNSIGNED_INT, nullptr));
		}

		void Destroy() {
			va.Destroy();
			vb.Destroy();
			ib.Destroy();
		}

		/// Static functions for generating common shapes

		/// <summary>
		/// Generates a rectangle centered at (x, y, z) = 0 with a width and height of 1
		/// </summary>
		/// <typeparam name="T">Data type used to represent vertices, normals, and tetxure coordinates.</typeparam>
		/// <returns></returns>
		template <typename T>
		static glGeometry GenerateRectangle() {
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
		static glGeometry GenerateCircle(unsigned int slices = 10) {
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
		static glGeometry GenerateCube(){
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
		static glGeometry GenerateSphere(unsigned int stacks, unsigned int slices) {
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
		static glGeometry GenerateIcosahedron() {
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
		static glGeometry GenerateIcosphere(unsigned int subdiv = 4, bool smooth = false) {
			tira::icosphere<float> ico(0.5f, subdiv, smooth);
			glGeometry icosphere;

			glVertexBufferLayout layout;
			layout.Push<T>(3);
			layout.Push<T>(3);
			layout.Push<T>(2);
			icosphere.AddMesh(ico.getInterleavedVertices(), ico.getInterleavedVertexCount() * 8 * sizeof(T), layout, ico.getIndices(), ico.getIndexCount(), 0);

			return icosphere;
		}

		


		template <typename T>
		static glGeometry GenerateSuperquadric(T l0, T l1, T l2, T gamma, unsigned int subdiv = 4, bool smooth = false) {
			tira::superquadric sq(l0, l1, l2, gamma, 0.5f, subdiv, smooth);

			glGeometry superquadric;

			glVertexBufferLayout layout;
			layout.Push<T>(3);
			layout.Push<T>(3);
			layout.Push<T>(2);
			superquadric.AddMesh(sq.getInterleavedVertices(), sq.getInterleavedVertexCount() * 8 * sizeof(T), layout, sq.getIndices(), sq.getIndexCount(), 0);

			return superquadric;
		}

		template <typename T>
		static glGeometry GenerateCylinder(unsigned int stacks, unsigned int slices) {
			cylinder<T> s(stacks, slices);
			std::vector<T> v = s.getInterleavedVertices();
			std::vector<unsigned int> i = s.getIndices();

			glGeometry cylinder;

			glVertexBufferLayout layout;
			layout.Push<T>(3);
			layout.Push<T>(3);
			layout.Push<T>(2);
			cylinder.AddMesh(&v[0], v.size() * sizeof(T), layout, &i[0], i.size(), 0);

			return cylinder;
		}
	};
}