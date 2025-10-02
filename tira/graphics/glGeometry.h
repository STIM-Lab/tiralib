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

		void AddMesh(tmesh mesh, int texcoords) {
			std::vector<float> v = mesh.Interleave(texcoords);									// generate a vector of interleaved vertices
			std::vector<unsigned int> i = mesh.Indices();									// generate a vector of indices
			glVertexBufferLayout layout;
			layout.Push<float>(3);													// add the number of vertex components to the layout
			layout.Push<float>(3);													// add the number of normal components to the layout
			layout.Push<float>(texcoords);													// add the number of texture components to the layout
			AddMesh(&v[0], v.size() * sizeof(float), layout, &i[0], i.size(), 0);
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

		glGeometry(tmesh mesh, int texcoords = 2) {
			AddMesh(mesh, texcoords);
		}


		/// Static functions for generating common shapes

		/// <summary>
		/// Generates a rectangle centered at (x, y, z) = 0 with a width and height of 1
		/// </summary>
		/// <typeparam name="T">Data type used to represent vertices, normals, and tetxure coordinates.</typeparam>
		/// <returns></returns>
		static glGeometry GenerateRectangle() {
			tmesh r = rectangle();
			glGeometry rectangle;
			rectangle.AddMesh(r, 2);

			return rectangle;
		}

		static glGeometry GenerateCircle(unsigned int slices = 10) {

			tmesh c = circle(slices);
			glGeometry circle(c, 2);
			return circle;
		}

		static glGeometry GenerateCube(){
			tmesh c = cube();
			glGeometry cube(c, 2);
			return cube;
		}

		static glGeometry GenerateSphere(unsigned int stacks, unsigned int slices) {
			tmesh s = sphere(slices, stacks);
			glGeometry sphere(s, 2);
			return sphere;
		}

		static glGeometry GenerateIcosahedron() {

			tmesh ico = icosahedron();
			glGeometry icosahedron(ico);
			return icosahedron;
		}

		static glGeometry GenerateIcosphere(unsigned int subdiv = 4, bool smooth = false) {
			tmesh ico = icosphere(0.5f, subdiv, smooth);
			glGeometry icosphere(ico, 2);
			return icosphere;
		}


		static glGeometry GenerateSuperquadric(float l0, float l1, float l2, float gamma, unsigned int subdiv = 4, bool smooth = false) {

			tmesh sq = superquadric(l0, l1, l2, gamma, subdiv, smooth);
			glGeometry superquadric(sq);
			return superquadric;
		}

		static glGeometry GenerateCylinder(unsigned int stacks, unsigned int slices) {

			tmesh c = cylinder(slices, stacks);
			glGeometry cylinder(c, 2);
			return cylinder;
		}
	};
}