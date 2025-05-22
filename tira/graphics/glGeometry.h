# pragma once

#include "glErrorHandler.h"
#include "glVertexBuffer.h"
#include "glVertexBufferLayout.h"
#include "glIndexBuffer.h"
#include "glVertexArray.h"
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


		glGeometry() {}

		template<typename T>
		glGeometry(trimesh<T> mesh) {
			std::vector<T> v = mesh.getInterleavedVertices();									// generate a vector of interleaved vertices
			std::vector<unsigned int> i = mesh.getIndices();									// generate a vector of indices
			glVertexBufferLayout layout;
			layout.Push<T>(mesh.getVertexDim());													// add the number of vertex components to the layout
			layout.Push<T>(mesh.getNormalDim());													// add the number of normal components to the layout
			layout.Push<T>(mesh.getTextureDim());													// add the number of texture components to the layout
			AddMesh(&v[0], v.size() * sizeof(T), layout, &i[0], i.size(), 0);
		}


		/// Static functions for generating common shapes

		/// <summary>
		/// Generates a rectangle centered at (x, y, z) = 0 with a width and height of 1
		/// </summary>
		/// <typeparam name="T">Data type used to represent vertices, normals, and tetxure coordinates.</typeparam>
		/// <returns></returns>
		template <typename T>
		static glGeometry GenerateRectangle() {
			trimesh<T> r = rectangle<T>();
			std::vector<T> v = r.getInterleavedVertices();									// generate a vector of interleaved vertices
			const std::vector<unsigned int> i = r.getIndices();									// generate a vector of indices
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

			trimesh<T> c = circle<T>(slices);

			std::vector<T> v = c.getInterleavedVertices();									// generate a vector of interleaved vertices
			std::vector<unsigned int> i = c.getIndices();									// generate a vector of indices


			glVertexBufferLayout layout;
			layout.Push<T>(3);																// add three vertex components to the layout
			layout.Push<T>(3);																// add three normal components to the layout
			layout.Push<T>(2);																// add two texture components to the layout

			glGeometry circle;
			circle.AddMesh(&v[0], v.size() * sizeof(T), layout, &i[0], i.size(), 0);
			return circle;
		}

	    template <typename T>
		static glGeometry GenerateCube(){
			trimesh<T> c = circle<T>();
			//cube<T> c;
			std::vector<T> v = c.getInterleavedVertices();
			std::vector<unsigned int> i = c.getIndices();


			glVertexBufferLayout layout;
			layout.Push<T>(3);
			layout.Push<T>(3);
			layout.Push<T>(2);

			glGeometry cube;
			cube.AddMesh(&v[0], v.size() * sizeof(T), layout, &i[0], i.size(), 0);
			return cube;
		}

		template <typename T>
		static glGeometry GenerateSphere(unsigned int stacks, unsigned int slices) {
			//sphere<T> s(stacks, slices);
			trimesh<T> s = sphere<T>(slices, stacks);
			std::vector<T> v = s.getInterleavedVertices();
			std::vector<unsigned int> i = s.getIndices();


			glVertexBufferLayout layout;
			layout.Push<T>(3);
			layout.Push<T>(3);
			layout.Push<T>(2);

			glGeometry sphere;
			sphere.AddMesh(&v[0], v.size() * sizeof(T), layout, &i[0], i.size(), 0);
			return sphere;
		}

		template <typename T>
		static glGeometry GenerateIcosahedron() {
			//tira::icosahedron<T> ico(0.5f);
			//
			trimesh<T> ico = icosahedron<T>();

			glVertexBufferLayout layout;
			layout.Push<T>(3);
			layout.Push<T>(3);
			layout.Push<T>(2);

			glGeometry icosahedron;
			icosahedron.AddMesh(ico.getInterleavedVertices(), ico.getInterleavedVertexCount() * 8 * sizeof(T), layout, ico.getIndices(), ico.getIndexCount(), 0);
			return icosahedron;
		}

		template <typename T>
		static glGeometry GenerateIcosphere(unsigned int subdiv = 4, bool smooth = false) {
			trimesh<T> ico = icosphere(0.5f, subdiv, smooth);


			glVertexBufferLayout layout;
			layout.Push<T>(3);
			layout.Push<T>(3);
			layout.Push<T>(2);

			glGeometry icosphere;
			icosphere.AddMesh(ico.getInterleavedVertices(), ico.getInterleavedVertexCount() * 8 * sizeof(T), layout, ico.getIndices(), ico.getIndexCount(), 0);
			return icosphere;
		}

		


		template <typename T>
		static glGeometry GenerateSuperquadric(T l0, T l1, T l2, T gamma, unsigned int subdiv = 4, bool smooth = false) {

			trimesh<T> sq = superquadric<T>(l0, l1, l2, gamma, subdiv, smooth);

			glVertexBufferLayout layout;
			layout.Push<T>(3);
			layout.Push<T>(3);
			layout.Push<T>(2);

			glGeometry superquadric;
			superquadric.AddMesh(sq.getInterleavedVertices(), sq.getInterleavedVertexCount() * 8 * sizeof(T), layout, sq.getIndices(), sq.getIndexCount(), 0);
			return superquadric;
		}

		template <typename T>
		static glGeometry GenerateCylinder(unsigned int stacks, unsigned int slices) {

			trimesh<T> c = cylinder<T>(slices, stacks);

			std::vector<T> v = c.getInterleavedVertices();
			std::vector<unsigned int> i = c.getIndices();



			glVertexBufferLayout layout;
			layout.Push<T>(3);
			layout.Push<T>(3);
			layout.Push<T>(2);

			glGeometry cylinder;
			cylinder.AddMesh(&v[0], v.size() * sizeof(T), layout, &i[0], i.size(), 0);
			return cylinder;
		}
	};
}