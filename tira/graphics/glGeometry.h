# pragma once

#include "glErrorHandler.h"
#include "glVertexBuffer.h"
#include "glVertexBufferLayout.h"
#include "glIndexBuffer.h"
#include "glVertexArray.h"
#include <tira/shapes.h>

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

		glGeometry(tmesh mesh) {
			std::vector<float> v = mesh.Interleave(2);									// generate a vector of interleaved vertices
			std::vector<unsigned int> i = mesh.Indices();									// generate a vector of indices
			glVertexBufferLayout layout;
			layout.Push<float>(3);													// add the number of vertex components to the layout
			layout.Push<float>(3);													// add the number of normal components to the layout
			layout.Push<float>(2);													// add the number of texture components to the layout
			size_t stride = (3 + 3 + 2) * sizeof(float);
			AddMesh(&v[0], mesh.NumVertices() * stride, layout, &i[0], i.size(), 0);
		}


		/// Static functions for generating common shapes

		/// <summary>
		/// Generates a rectangle centered at (x, y, z) = 0 with a width and height of 1
		/// </summary>
		/// <typeparam name="T">Data type used to represent vertices, normals, and tetxure coordinates.</typeparam>
		/// <returns></returns>
		static glGeometry GenerateRectangle() {
			tmesh r = rectangle();
			std::vector<float> v = r.Interleave(2);									// generate a vector of interleaved vertices
			const std::vector<unsigned int> i = r.Indices();									// generate a vector of indices
			glGeometry rectangle;

			glVertexBufferLayout layout;
			layout.Push<float>(3);																// add three vertex components to the layout
			layout.Push<float>(3);																// add three normal components to the layout
			layout.Push<float>(2);																// add two texture components to the layout
			size_t stride = (3 + 3 + 2) * sizeof(float);
			rectangle.AddMesh(&v[0], r.NumVertices() * stride, layout, &i[0], i.size(), 0);

			return rectangle;
		}

		static glGeometry GenerateCircle(unsigned int slices = 10) {

			tmesh c = circle(slices);

			std::vector<float> v = c.Interleave(2);									// generate a vector of interleaved vertices
			std::vector<unsigned int> i = c.Indices();									// generate a vector of indices


			glVertexBufferLayout layout;
			layout.Push<float>(3);																// add three vertex components to the layout
			layout.Push<float>(3);																// add three normal components to the layout
			layout.Push<float>(2);																// add two texture components to the layout
			size_t stride = (3 + 3 + 2) * sizeof(float);

			glGeometry circle;
			circle.AddMesh(&v[0], c.NumVertices() * stride, layout, &i[0], i.size(), 0);
			return circle;
		}

		static glGeometry GenerateCube(){

			tmesh c = cube();
			//cube<T> c;
			std::vector<float> v = c.Interleave(2);
			std::vector<unsigned int> i = c.Indices();


			glVertexBufferLayout layout;
			layout.Push<float>(3);
			layout.Push<float>(3);
			layout.Push<float>(2);
			size_t stride = (3 + 3 + 2) * sizeof(float);

			glGeometry cube;
			cube.AddMesh(&v[0], c.NumVertices() * stride, layout, &i[0], i.size(), 0);
			return cube;
		}

		template <typename T>
		static glGeometry GenerateSphere(unsigned int stacks, unsigned int slices) {
			//sphere<T> s(stacks, slices);
			tmesh s = sphere(stacks, slices);
			s = sphere(slices, stacks);
			std::vector<float> v = s.Interleave(2);
			std::vector<unsigned int> i = s.Indices();


			glVertexBufferLayout layout;
			layout.Push<float>(3);
			layout.Push<float>(3);
			layout.Push<float>(2);
			size_t stride = (3 + 3 + 2) * sizeof(float);

			glGeometry sphere;
			sphere.AddMesh(&v[0], s.NumVertices() * stride, layout, &i[0], i.size(), 0);
			return sphere;
		}

		static glGeometry GenerateIcosahedron() {
			//tira::icosahedron<T> ico(0.5f);
			//
			tira::tmesh ico = icosohedron();
			std::vector<float> v = ico.Interleave(2);

			glVertexBufferLayout layout;
			layout.Push<float>(3);
			layout.Push<float>(3);
			layout.Push<float>(2);
			size_t stride = (3 + 3 + 2) * sizeof(float);

			glGeometry icosahedron;
			icosahedron.AddMesh(&v[0], ico.NumVertices() * stride, layout, &ico.Indices()[0], ico.NumIndices(), 0);
			return icosahedron;
		}

		static glGeometry GenerateIcosphere(unsigned int subdiv = 4, bool smooth = false) {
			tmesh ico = icosphere(0.5f, subdiv, smooth);
			std::vector<float> v = ico.Interleave(2);

			glVertexBufferLayout layout;
			layout.Push<float>(3);
			layout.Push<float>(3);
			layout.Push<float>(2);
			size_t stride = (3 + 3 + 2) * sizeof(float);

			glGeometry icosphere;
			icosphere.AddMesh(&v[0], ico.NumVertices() * stride, layout, &ico.Indices()[0], ico.NumIndices(), 0);
			return icosphere;
		}


		static glGeometry GenerateSuperquadric(float l0, float l1, float l2, float gamma, unsigned int subdiv = 4, bool smooth = false) {

			tmesh sq = superquadric(l0, l1, l2, gamma, subdiv, smooth);
			std::vector<float> v = sq.Interleave(2);

			glVertexBufferLayout layout;
			layout.Push<float>(3);
			layout.Push<float>(3);
			layout.Push<float>(2);
			size_t stride = (3 + 3 + 2) * sizeof(float);


			glGeometry superquadric;
			superquadric.AddMesh(&v[0], sq.NumVertices() * stride, layout, &sq.Indices()[0], sq.NumIndices(), 0);
			return superquadric;
		}

		static glGeometry GenerateCylinder(unsigned int stacks, unsigned int slices) {

			tmesh c = cylinder(slices, stacks);

			std::vector<float> v = c.Interleave(2);
			std::vector<unsigned int> i = c.Indices();



			glVertexBufferLayout layout;
			layout.Push<float>(3);
			layout.Push<float>(3);
			layout.Push<float>(2);
			size_t stride = (3 + 3 + 2) * sizeof(float);

			glGeometry cylinder;
			cylinder.AddMesh(&v[0], c.NumVertices() * stride, layout, &i[0], i.size(), 0);
			return cylinder;
		}
	};
}