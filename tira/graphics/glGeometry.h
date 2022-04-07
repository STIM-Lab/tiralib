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
			ib.Bind();
		}

		void Draw() {
			va.Bind();              // Bind index array
			ib.Bind();              // Bind index buffer
			GLCALL(glDrawElements(GL_TRIANGLES, ib.GetCount(), GL_UNSIGNED_INT, nullptr));
		}
	};
}
//namespace tira {
//	glGeometry::glGeometry() {}
//
//	void glGeometry::AddMesh(const void* vertices, unsigned int bytes,
//		glVertexBufferLayout layout,
//		const unsigned int* indices, unsigned int count, unsigned int offset) {
//		vb.SetBuffer(vertices, bytes);
//		ib.SetBuffer(indices, count);
//		va.AddBuffer(vb, layout, offset);
//		ib.Bind();
//	}
//
//	void glGeometry::Unbind() {
//		va.Unbind();
//		vb.Unbind();
//		ib.Unbind();
//	}
//
//	void glGeometry::Draw() {
//		va.Bind();              // Bind index array
//		ib.Bind();              // Bind index buffer
//		GLCALL(glDrawElements(GL_TRIANGLES, ib.GetCount(), GL_UNSIGNED_INT, nullptr));
//	}
//}