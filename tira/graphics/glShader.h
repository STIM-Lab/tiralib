# pragma once
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "glErrorHandler.h"
#include <glm/glm.hpp>

namespace tira {
	/// This static class provides the STIM interface for loading, saving, and storing 2D images.
	/// Data is stored in an interleaved (BIP) format (default for saving and loading is RGB).
	struct ShaderProgramSource {
		std::string VertexSource;       // Public field so use Capital. 
		std::string FragmentSource;
	};

	struct glShaderUniform{

		std::string name;
		GLenum type;
		GLint location;

		GLint size() {
			switch (type) {
			case GL_FLOAT:
				return sizeof(GLfloat);
			case GL_FLOAT_VEC2:
				return 2 * sizeof(GLfloat);
			case GL_FLOAT_VEC3:
				return 3 * sizeof(GLfloat);
			case GL_FLOAT_VEC4:
				return 4 * sizeof(GLfloat);
			case GL_DOUBLE:
				return sizeof(GLdouble);
			case GL_DOUBLE_VEC2:
				return 2 * sizeof(GLdouble);
			case GL_DOUBLE_VEC3:
				return 3 * sizeof(GLdouble);
			case GL_DOUBLE_VEC4:
				return 4 * sizeof(GLdouble);
			case GL_INT:
				return sizeof(GLint);
			case GL_INT_VEC2:
				return 2 * sizeof(GLint);
			case GL_INT_VEC3:
				return 3 * sizeof(GLint);
			case GL_INT_VEC4:
				return 4 * sizeof(GLint);
			case GL_UNSIGNED_INT:
				return sizeof(GLuint);
			case GL_UNSIGNED_INT_VEC2:
				return 2 * sizeof(GLuint);
			case GL_UNSIGNED_INT_VEC3:
				return 3 * sizeof(GLuint);
			case GL_UNSIGNED_INT_VEC4:
				return 4 * sizeof(GLuint);
			case GL_BOOL:
				return sizeof(GLboolean);
			case GL_BOOL_VEC2:
				return 2 * sizeof(GLboolean);
			case GL_BOOL_VEC3:
				return 3 * sizeof(GLboolean);
			case GL_BOOL_VEC4:
				return 4 * sizeof(GLboolean);
			case GL_FLOAT_MAT2:
				return 2 * 2 * sizeof(GLfloat);
			case GL_FLOAT_MAT3:
				return 3 * 3 * sizeof(GLfloat);
			case GL_FLOAT_MAT4:
				return 4 * 4 * sizeof(GLfloat);
			case GL_FLOAT_MAT2x3:
				return 2 * 3 * sizeof(GLfloat);
			case GL_FLOAT_MAT2x4:
				return 2 * 4 * sizeof(GLfloat);
			case GL_FLOAT_MAT3x2:
				return 3 * 2 * sizeof(GLfloat);
			case GL_FLOAT_MAT3x4:
				return 3 * 4 * sizeof(GLfloat);
			case GL_FLOAT_MAT4x2:
				return 4 * 2 * sizeof(GLfloat);
			case GL_FLOAT_MAT4x3:
				return 4 * 3 * sizeof(GLfloat);
			case GL_DOUBLE_MAT2:
				return 2 * 2 * sizeof(GLdouble);
			case GL_DOUBLE_MAT3:
				return 3 * 3 * sizeof(GLdouble);
			case GL_DOUBLE_MAT4:
				return 4 * 4 * sizeof(GLdouble);
			case GL_DOUBLE_MAT2x3:
				return 2 * 3 * sizeof(GLdouble);
			case GL_DOUBLE_MAT2x4:
				return 2 * 4 * sizeof(GLdouble);
			case GL_DOUBLE_MAT3x2:
				return 3 * 2 * sizeof(GLdouble);
			case GL_DOUBLE_MAT3x4:
				return 3 * 4 * sizeof(GLdouble);
			case GL_DOUBLE_MAT4x2:
				return 4 * 2 * sizeof(GLdouble);
			case GL_DOUBLE_MAT4x3:
				return 4 * 3 * sizeof(GLdouble);
			case GL_SAMPLER_1D:
			case GL_SAMPLER_2D:
			case GL_SAMPLER_3D:
			case GL_SAMPLER_CUBE:
			case GL_SAMPLER_1D_SHADOW:
			case GL_SAMPLER_2D_SHADOW:
			case GL_SAMPLER_1D_ARRAY:
			case GL_SAMPLER_2D_ARRAY:
			case GL_SAMPLER_1D_ARRAY_SHADOW:
			case GL_SAMPLER_2D_ARRAY_SHADOW:
			case GL_SAMPLER_2D_MULTISAMPLE:
			case GL_SAMPLER_2D_MULTISAMPLE_ARRAY:
			case GL_SAMPLER_CUBE_SHADOW:
			case GL_SAMPLER_BUFFER:
			case GL_SAMPLER_2D_RECT:
			case GL_SAMPLER_2D_RECT_SHADOW:
			case GL_INT_SAMPLER_1D:
			case GL_INT_SAMPLER_2D:
			case GL_INT_SAMPLER_3D:
			case GL_INT_SAMPLER_CUBE:
			case GL_INT_SAMPLER_1D_ARRAY:
			case GL_INT_SAMPLER_2D_ARRAY:
			case GL_INT_SAMPLER_2D_MULTISAMPLE:
			case GL_INT_SAMPLER_2D_MULTISAMPLE_ARRAY:
			case GL_INT_SAMPLER_BUFFER:
			case GL_INT_SAMPLER_2D_RECT:
			case GL_UNSIGNED_INT_SAMPLER_1D:
			case GL_UNSIGNED_INT_SAMPLER_2D:
			case GL_UNSIGNED_INT_SAMPLER_3D:
			case GL_UNSIGNED_INT_SAMPLER_CUBE:
			case GL_UNSIGNED_INT_SAMPLER_1D_ARRAY:
			case GL_UNSIGNED_INT_SAMPLER_2D_ARRAY:
			case GL_UNSIGNED_INT_SAMPLER_2D_MULTISAMPLE:
			case GL_UNSIGNED_INT_SAMPLER_2D_MULTISAMPLE_ARRAY:
			case GL_UNSIGNED_INT_SAMPLER_BUFFER:
			case GL_UNSIGNED_INT_SAMPLER_2D_RECT:
			case GL_IMAGE_1D:
			case GL_IMAGE_2D:
			case GL_IMAGE_3D:
			case GL_IMAGE_2D_RECT:
			case GL_IMAGE_CUBE:
			case GL_IMAGE_BUFFER:
			case GL_IMAGE_1D_ARRAY:
			case GL_IMAGE_2D_ARRAY:
			case GL_IMAGE_2D_MULTISAMPLE:
			case GL_IMAGE_2D_MULTISAMPLE_ARRAY:
			case GL_INT_IMAGE_1D:
			case GL_INT_IMAGE_2D:
			case GL_INT_IMAGE_3D:
			case GL_INT_IMAGE_2D_RECT:
			case GL_INT_IMAGE_CUBE:
			case GL_INT_IMAGE_BUFFER:
			case GL_INT_IMAGE_1D_ARRAY:
			case GL_INT_IMAGE_2D_ARRAY:
			case GL_INT_IMAGE_2D_MULTISAMPLE:
			case GL_INT_IMAGE_2D_MULTISAMPLE_ARRAY:
			case GL_UNSIGNED_INT_IMAGE_1D:
			case GL_UNSIGNED_INT_IMAGE_2D:
			case GL_UNSIGNED_INT_IMAGE_3D:
			case GL_UNSIGNED_INT_IMAGE_2D_RECT:
			case GL_UNSIGNED_INT_IMAGE_CUBE:
			case GL_UNSIGNED_INT_IMAGE_BUFFER:
			case GL_UNSIGNED_INT_IMAGE_1D_ARRAY:
			case GL_UNSIGNED_INT_IMAGE_2D_ARRAY:
			case GL_UNSIGNED_INT_IMAGE_2D_MULTISAMPLE:
			case GL_UNSIGNED_INT_IMAGE_2D_MULTISAMPLE_ARRAY:
			case GL_UNSIGNED_INT_ATOMIC_COUNTER:
			default:
				return 0;
			}
		}
	};


	class glShader {
	protected:
		std::string m_FilePath;
		unsigned int m_ShaderID;
		mutable std::unordered_map<std::string, glShaderUniform> m_UniformCache;
		// Caching for uniforms
	public:
		/// <summary>
		/// Load a file containing both a fragment and vertex shader
		/// The input file should have the shaders labeled with the following:
		/// # shader vertex
		/// ....
		/// # shader fragment
		/// ....
		/// </summary>
		/// <param name="filepath">File path and name</param>
		glShader() : m_ShaderID(0) {}
		glShader(const std::string& filepath) {
			m_ShaderID = 0;
			LoadShader(filepath);
		}
		glShader(const std::string& vertexSource, const std::string& fragmentSource) {
			m_ShaderID = 0;
			CreateShader(vertexSource, fragmentSource);
		}

		void CreateShader(const std::string& vertexShader, const std::string& fragmentShader) {
			if (m_ShaderID != 0)
				glDeleteProgram(m_ShaderID);
			m_ShaderID = glCreateProgram();
			unsigned int vs = CompileShader(GL_VERTEX_SHADER, vertexShader);
			unsigned int fs = CompileShader(GL_FRAGMENT_SHADER, fragmentShader);
			GLERROR(glAttachShader(m_ShaderID, vs));
			GLERROR(glAttachShader(m_ShaderID, fs));
			GLERROR(glLinkProgram(m_ShaderID));
			GLERROR(glValidateProgram(m_ShaderID));

			GLERROR(glDeleteShader(vs));
			GLERROR(glDeleteShader(fs));
			CacheUniforms();							// cache the new uniform variables
		}

		void LoadShader(const std::string& filepath) {
			ShaderProgramSource source = ParseShader(filepath);
			CreateShader(source.VertexSource, source.FragmentSource);
		}

		void Bind() const {
			GLERROR(glUseProgram(m_ShaderID));
		}

		void Unbind() const {
			GLERROR(glUseProgram(0));
		}

		// Set uniforms
		void SetUniform1f(const std::string& name, float value) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniform1f(location, value));
		}
		void SetUniform2f(const std::string& name, float v0, float v1) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniform2f(location, v0, v1));
		}
		void SetUniform3f(const std::string& name, float v0, float v1, float v2) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniform3f(location, v0, v1, v2));
		}
		void SetUniform4f(const std::string& name, float v0, float v1, float v2, float v3) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniform4f(location, v0, v1, v2, v3));
		}

		void SetUniform1i(const std::string& name, int value) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniform1i(location, value));
		}
		void SetUniform2i(const std::string& name, int v0, int v1) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniform2i(location, v0, v1));
		}
		void SetUniform3i(const std::string& name, int v0, int v1, int v2) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniform3i(location, v0, v1, v2));
		}
		void SetUniform4i(const std::string& name, int v0, int v1, int v2, int v3) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniform4i(location, v0, v1, v2, v3));
		}

		void SetUniform1ui(const std::string& name, unsigned int value) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniform1ui(location, value));
		}
		void SetUniform2ui(const std::string& name, unsigned int v0, unsigned int v1) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniform2ui(location, v0, v1));
		}
		void SetUniform3ui(const std::string& name, unsigned int v0, unsigned int v1, unsigned int v2) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniform3ui(location, v0, v1, v2));
		}
		void SetUniform4ui(const std::string& name, unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniform4ui(location, v0, v1, v2, v3));
		}

		inline void SetUniformMat4f(const std::string& name, const glm::mat4& matrix) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniformMatrix4fv(location, 1, GL_FALSE, &matrix[0][0]));
		}
		inline void SetUniformMat3f(const std::string& name, const glm::mat3& matrix) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniformMatrix3fv(location, 1, GL_FALSE, &matrix[0][0]));
		}
		inline void SetUniformMat2f(const std::string& name, const glm::mat2& matrix) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniformMatrix2fv(location, 1, GL_FALSE, &matrix[0][0]));
		}
		inline void SetUniformMat2x3f(const std::string& name, const glm::mat2x3& matrix) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniformMatrix2x3fv(location, 1, GL_FALSE, &matrix[0][0]));
		}
		inline void SetUniformMat3x2f(const std::string& name, const glm::mat3x2& matrix) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniformMatrix3x2fv(location, 1, GL_FALSE, &matrix[0][0]));
		}
		inline void SetUniformMat2x4f(const std::string& name, const glm::mat2x4& matrix) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniformMatrix2x4fv(location, 1, GL_FALSE, &matrix[0][0]));
		}
		inline void SetUniformMat4x2f(const std::string& name, const glm::mat4x2& matrix) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniformMatrix4x2fv(location, 1, GL_FALSE, &matrix[0][0]));
		}
		inline void SetUniformMat3x4f(const std::string& name, const glm::mat3x4& matrix) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniformMatrix3x4fv(location, 1, GL_FALSE, &matrix[0][0]));
		}
		inline void SetUniformMat4x3f(const std::string& name, const glm::mat4x3& matrix) {
			GLint location = GetUniformLocation(name);
			GLERROR(glUniformMatrix4x3fv(location, 1, GL_FALSE, &matrix[0][0]));
		}

		std::vector< glShaderUniform > GetUniformList(){
			if (m_ShaderID == 0) return std::vector<glShaderUniform>();	// if there isn't a shader program, return an empty vector
			GLint count;												// stores the number of uniforms
			glGetProgramiv(m_ShaderID, GL_ACTIVE_UNIFORMS, &count);		// get the number of uniforms

			GLint s; 							// size of the variable
			GLenum type; 							// type of the variable (float, vec3 or mat4, etc)

			const GLsizei bufSize = 16; 			// maximum name length
			GLchar name[bufSize]; 					// variable name in GLSL
			GLsizei length; 						// name length

			std::vector< glShaderUniform > uniforms(count);
			for(int i = 0; i < count; i++){
				glGetActiveUniform(m_ShaderID, (GLuint)i, bufSize, &length, &s, &type, name);
				uniforms[i].name = std::string(name);
				uniforms[i].type = type;
				uniforms[i].location = glGetUniformLocation(m_ShaderID, name);
			}
			return uniforms;
		}

	protected:

		ShaderProgramSource ParseShader(const std::string& filepath) {             // Slow
			std::ifstream stream(filepath);
			enum class ShaderType {
				NONE = -1, VERTEX = 0, FRAGMENT = 1
			};
			std::string line;
			std::stringstream ss[2];
			ShaderType type = ShaderType::NONE;
			while (getline(stream, line)) {
				if (line.find("# shader") != std::string::npos) {
					if (line.find("vertex") != std::string::npos) {
						// set mode to vertex
						type = ShaderType::VERTEX;
					}
					else if (line.find("fragment") != std::string::npos) {
						// set mode to fragment
						type = ShaderType::FRAGMENT;
					}
				}
				else {
					ss[(int)type] << line << "\n";
				}
			}
			return { ss[0].str(), ss[1].str() };
		}
		unsigned int CompileShader(unsigned int type, const std::string& source) {
			unsigned int id = glCreateShader(type);
			const char* src = source.c_str();
			GLERROR(glShaderSource(id, 1, &src, nullptr));
			GLERROR(glCompileShader(id));

			int result;
			glGetShaderiv(id, GL_COMPILE_STATUS, &result);
			if (result == GL_FALSE) {
				int length;
				GLERROR(glGetShaderiv(id, GL_INFO_LOG_LENGTH, &length));
				#ifdef _MSC_VER
				char* message = (char*)_malloca(length * sizeof(char));
				#else
				char* message = (char*)alloca(length * sizeof(char));
				#endif
				GLERROR(glGetShaderInfoLog(id, length, &length, message));
				std::cout << "Failed to compile " << (type == GL_VERTEX_SHADER ? "vertex" : "fragment") << "shader" << std::endl;
				std::cout << message << std::endl;
				GLERROR(glDeleteShader(id));
				return 0;
			}

			return id;
			//TODO: Error handling
		}
		
		GLint GetUniformLocation(const std::string& name) const {
			if (m_UniformCache.find(name) != m_UniformCache.end())
				return m_UniformCache[name].location;
			return -1;
		}

		void CacheUniforms(){
			std::vector< glShaderUniform > uniform_list = GetUniformList();
			m_UniformCache.clear();
			for(size_t i = 0; i < uniform_list.size(); i++){
				m_UniformCache[uniform_list[i].name] = uniform_list[i];
			}
		}

	};
}