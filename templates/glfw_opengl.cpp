#include <GL/glew.h>
#include <GLFW/glfw3.h> // Will drag system OpenGL headers

#include "tira/graphics_gl.h"

GLFWwindow* window;                                     // pointer to the GLFW window that will be created (used in GLFW calls to request properties)
const char* glsl_version = "#version 130";              // specify the version of GLSL
tira::camera camera;


void glfw_error_callback(int error, const char* description)
{
	fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

std::string vertex_shader_source =
"#version 330 core\n"
"layout(location = 0) in vec3 aPos;\n"
"uniform mat4 PVM;\n "
"out vec4 vertexColor;\n"
"void main() {\n"
"	gl_Position = PVM * vec4(aPos, 1.0);\n"
"	vertexColor = vec4(1.0, 0.0, 0.0, 1.0);\n"
"}\n";

std::string fragment_shader_source =
"#version 330 core\n"
"out vec4 FragColor;\n"
"in vec4 vertexColor;\n"
"void main() {\n"
"	FragColor = vertexColor;\n"
"}\n";

GLFWwindow* InitGLFW() {
	GLFWwindow* window;

	// Setup window
	glfwSetErrorCallback(glfw_error_callback);
	if (!glfwInit())
		return NULL;

	// GL 3.0 + GLSL 130
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

	// Create window with graphics context
	window = glfwCreateWindow(1600, 1200, "GLFW+OpenGL3 Hello World Program", NULL, NULL);
	if (window == NULL)
		return NULL;
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1); // Enable vsync
	return window;
}

void InitGLEW() {
	GLenum err = glewInit();
	if (GLEW_OK != err) {
		/* Problem: glewInit failed, something is seriously wrong. */
		fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
		exit(1);
	}
}

int main(int argc, char** argv) {
	// Initialize OpenGL
	window = InitGLFW();                                // create a GLFW window
	InitGLEW();

	// Enable OpenGL environment parameters
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	
	camera.setPosition(0, 0, 1);
	camera.LookAt(0, 0, 0);

	glm::mat4 Mview = camera.getMatrix();
	

	tira::glGeometry square = tira::glGeometry::GenerateRectangle<float>();	// create a square
	tira::glShader shader(vertex_shader_source, fragment_shader_source);

	

	while (!glfwWindowShouldClose(window)){

		
		glfwPollEvents();													// Poll and handle events (inputs, window resize, etc.)

		int display_w, display_h;                                           // size of the frame buffer (openGL display)
		glfwGetFramebufferSize(window, &display_w, &display_h);

		float aspect = (float)display_w / (float)display_h;
		glm::mat4 Mprojection;
		if (aspect > 1) {
			Mprojection = glm::ortho(-aspect, aspect, -1.0f, 1.0f, -10.0f, 10.0f);
		}
		else {
			Mprojection = glm::ortho(-1.0f, 1.0f, -1.0f/aspect, 1.0f/aspect, -10.0f, 10.0f);
		}

		glViewport(0, 0, display_w, display_h);								// specifies the area of the window where OpenGL can render
		
		glClearColor(0, 0, 0, 0);
		
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		
		/// Render Something Here
		shader.Bind();
		shader.SetUniformMat4f("PVM", Mprojection * Mview);
		
		square.Draw();


		glfwSwapBuffers(window);

	}

	glfwDestroyWindow(window);                                      // Destroy the GLFW rendering window
	glfwTerminate();                                                // Terminate GLFW

}