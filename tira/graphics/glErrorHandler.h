# pragma once

#ifdef _MSC_VER
#define DEBUG_BREAK __debugbreak()
#else
#include <signal.h>
#define DEBUG_BREAK raise(SIGTRAP)
#endif

#include <GL/glew.h>
#include <iostream>

#define ASSERT(x) if (!(x)) DEBUG_BREAK;
#define GLBREAK(x) tira::glClearError();\
    x;\
    ASSERT(tira::glLogCall(#x, __FILE__, __LINE__))

#define GLERROR(x) tira::glClearError();\
    x;\
    tira::glLogCall(#x, __FILE__, __LINE__)

namespace tira {
    // For multiple errors. Loop until all errors found. 
    inline void glClearError() {
        while (glGetError() != GL_NO_ERROR);            
    }

    inline std::string glErrorString(GLenum error) {
        if(error == GL_INVALID_VALUE) return "GL_INVALID_VALUE";
        else if(error == GL_INVALID_ENUM) return "GL_INVALID_ENUM";
        else if(error == GL_INVALID_OPERATION) return "GL_INVALID_OPERATION";
        else if(error == GL_OUT_OF_MEMORY) return "GL_OUT_OF_MEMORY";
        else if(error == GL_STACK_OVERFLOW) return "GL_STACK_OVERFLOW";
        else if(error == GL_STACK_UNDERFLOW) return "GL_STACK_UNDERFLOW";
        else if(error == GL_INVALID_FRAMEBUFFER_OPERATION) return "GL_INVALID_FRAMEBUFFER_OPERATION";
        else if(error == GL_CONTEXT_LOST) return "GL_CONTEXT_LOST";
        else return "UNKNOWN ERROR??";
    }
        // For multiple errors. Loop until all errors found. 
    inline bool glLogCall(const char* function, const char* file, int line) {
        while (GLenum error = glGetError()) {
            std::cout << "[OpenGL Error] (" << error << ", "<<glErrorString(error)<<"):" << function << " " << file << ":" << line << std::endl;
            return false;
        }
        return true;
    }
}