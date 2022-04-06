# pragma once
#include <unordered_map>
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


	class glShader {
	private:
		std::string m_FilePath;
		unsigned int m_ShaderID;
		mutable std::unordered_map<std::string, GLint> m_UniformLocationCache;
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
		//~glShader() {
		//	GLCALL(glDeleteProgram(m_ShaderID));
		//}

		void LoadShader(const std::string& filepath) {
			ShaderProgramSource source = ParseShader(filepath);
			m_ShaderID = CreateShader(source.VertexSource, source.FragmentSource);
		}

		void Bind() const {
			GLCALL(glUseProgram(m_ShaderID));
		}

		void Unbind() const {
			GLCALL(glUseProgram(0));
		}

		// Set uniforms
		// Mimic larfer project
		void SetUniform1i(const std::string& name, float value) {
			GLint location = GetUniformLocation(name);
			GLCALL(glUniform1i(location, value));
		}
		void SetUniform1f(const std::string& name, float value) {
			GLint location = GetUniformLocation(name);
			GLCALL(glUniform1f(location, value));
		}
		void SetUniform3f(const std::string& name, float v0, float v1, float v2) {
			GLint location = GetUniformLocation(name);
			GLCALL(glUniform3f(location, v0, v1, v2));
		}
		void SetUniform4f(const std::string& name, float v0, float v1, float v2, float v3) {
			GLint location = GetUniformLocation(name);
			GLCALL(glUniform4f(location, v0, v1, v2, v3));
		}
		inline void SetUniformMat4f(const std::string& name, const glm::mat4& matrix) {
			GLint location = GetUniformLocation(name);
			GLCALL(glUniformMatrix4fv(location, 1, GL_FALSE, &matrix[0][0]));
		}

	private:

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
			GLCALL(glShaderSource(id, 1, &src, nullptr));
			GLCALL(glCompileShader(id));

			int result;
			glGetShaderiv(id, GL_COMPILE_STATUS, &result);
			if (result == GL_FALSE) {
				int length;
				GLCALL(glGetShaderiv(id, GL_INFO_LOG_LENGTH, &length));
				#ifdef _MSC_VER
				char* message = (char*)_malloca(length * sizeof(char));
				#else
				char* message = (char*)alloca(length * sizeof(char));
				#endif
				GLCALL(glGetShaderInfoLog(id, length, &length, message));
				std::cout << "Failed to compile " << (type == GL_VERTEX_SHADER ? "vertex" : "fragment") << "shader" << std::endl;
				std::cout << message << std::endl;
				GLCALL(glDeleteShader(id));
				return 0;
			}

			return id;
			//TODO: Error handling
		}
		unsigned int CreateShader(const std::string& vertexShader, const std::string& fragmentShader) {
			unsigned int program = glCreateProgram();
			unsigned int vs = CompileShader(GL_VERTEX_SHADER, vertexShader);
			unsigned int fs = CompileShader(GL_FRAGMENT_SHADER, fragmentShader);
			GLCALL(glAttachShader(program, vs));
			GLCALL(glAttachShader(program, fs));
			GLCALL(glLinkProgram(program));
			GLCALL(glValidateProgram(program));

			GLCALL(glDeleteShader(vs));
			GLCALL(glDeleteShader(fs));

			return program;
		}
		GLint GetUniformLocation(const std::string& name) const {
			if (m_UniformLocationCache.find(name) != m_UniformLocationCache.end())
				return m_UniformLocationCache[name];
			GLint location = glGetUniformLocation(m_ShaderID, name.c_str());
			m_UniformLocationCache[name] = location;
			return location;
		}

	};
}
//
//#include "glShader.h"
//#include <iostream>
//#include <fstream>
//#include <string>
//#include <sstream>
//
//#include "glGeometry.h"
//
//namespace tira {
//    glShader::glShader() : m_ShaderID(0) {}
//
//    glShader::glShader(const std::string& filepath) {
//        m_ShaderID = 0;
//        LoadShader(filepath);
//    }
//
//    glShader::~glShader() {
//        GLCALL(glDeleteProgram(m_ShaderID));
//    }
//
//    void glShader::LoadShader(const std::string& filepath) {
//        ShaderProgramSource source = ParseShader(filepath);
//        m_ShaderID = CreateShader(source.VertexSource, source.FragmentSource);
//    }
//
//    // Deal with Basic.shader file 
//    ShaderProgramSource glShader::ParseShader(const std::string& filepath) {             // Slow
//        std::ifstream stream(filepath);
//        enum class ShaderType {
//            NONE = -1, VERTEX = 0, FRAGMENT = 1
//        };
//        std::string line;
//        std::stringstream ss[2];
//        ShaderType type = ShaderType::NONE;
//        while (getline(stream, line)) {
//            if (line.find("# shader") != std::string::npos) {
//                if (line.find("vertex") != std::string::npos) {
//                    // set mode to vertex
//                    type = ShaderType::VERTEX;
//                }
//                else if (line.find("fragment") != std::string::npos) {
//                    // set mode to fragment
//                    type = ShaderType::FRAGMENT;
//                }
//            }
//            else {
//                ss[(int)type] << line << "\n";
//            }
//        }
//        return { ss[0].str(), ss[1].str() };
//    }
//
//    unsigned int glShader::CompileShader(unsigned int type, const std::string& source) {
//        unsigned int id = glCreateShader(type);
//        const char* src = source.c_str();
//        GLCALL(glShaderSource(id, 1, &src, nullptr));
//        GLCALL(glCompileShader(id));
//
//        int result;
//        glGetShaderiv(id, GL_COMPILE_STATUS, &result);
//        if (result == GL_FALSE) {
//            int length;
//            GLCALL(glGetShaderiv(id, GL_INFO_LOG_LENGTH, &length));
//            char* message = (char*)_malloca(length * sizeof(char));
//            GLCALL(glGetShaderInfoLog(id, length, &length, message));
//            std::cout << "Failed to compile " << (type == GL_VERTEX_SHADER ? "vertex" : "fragment") << "shader" << std::endl;
//            std::cout << message << std::endl;
//            GLCALL(glDeleteShader(id));
//            return 0;
//        }
//
//        return id;
//        //TODO: Error handling
//    }
//
//    unsigned int glShader::CreateShader(const std::string& vertexShader, const std::string& fragmentShader) {
//        unsigned int program = glCreateProgram();
//        unsigned int vs = CompileShader(GL_VERTEX_SHADER, vertexShader);
//        unsigned int fs = CompileShader(GL_FRAGMENT_SHADER, fragmentShader);
//        GLCALL(glAttachShader(program, vs));
//        GLCALL(glAttachShader(program, fs));
//        GLCALL(glLinkProgram(program));
//        GLCALL(glValidateProgram(program));
//
//        GLCALL(glDeleteShader(vs));
//        GLCALL(glDeleteShader(fs));
//
//        return program;
//    }
//
//    GLint glShader::GetUniformLocation(const std::string& name) const {
//        if (m_UniformLocationCache.find(name) != m_UniformLocationCache.end())
//            return m_UniformLocationCache[name];
//        GLint location = glGetUniformLocation(m_ShaderID, name.c_str());
//        m_UniformLocationCache[name] = location;
//        return location;
//    }
//
//
//    void glShader::Bind() const {
//        GLCALL(glUseProgram(m_ShaderID));
//    }
//
//    void glShader::Unbind() const {
//        GLCALL(glUseProgram(0));
//    }
//
//    void glShader::SetUniform1i(const std::string& name, float value) {
//        GLint location = GetUniformLocation(name);
//        GLCALL(glUniform1i(location, value));
//    }
//
//    void glShader::SetUniform1f(const std::string& name, float value) {
//        GLint location = GetUniformLocation(name);
//        GLCALL(glUniform1f(location, value));
//    }
//
//    void glShader::SetUniform3f(const std::string& name, float v0, float v1, float v2) {
//        GLint location = GetUniformLocation(name);
//        GLCALL(glUniform3f(location, v0, v1, v2));
//    }
//    void glShader::SetUniform4f(const std::string& name, float v0, float v1, float v2, float v3) {
//        GLint location = GetUniformLocation(name);
//        GLCALL(glUniform4f(location, v0, v1, v2, v3));
//    }
//    void glShader::SetUniformMat4f(const std::string& name, const glm::mat4& matrix) {
//        GLint location = GetUniformLocation(name);
//        GLCALL(glUniformMatrix4fv(location, 1, GL_FALSE, &matrix[0][0]));
//    }
//}
//
