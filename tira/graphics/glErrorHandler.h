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
#define GLCALL(x) tira::glClearError();\
    x;\
    ASSERT(tira::glLogCall(#x, __FILE__, __LINE__))

namespace tira {
    // For multiple errors. Loop until all errors found. 
    inline void glClearError() {
        while (glGetError() != GL_NO_ERROR);            
    }
        // For multiple errors. Loop until all errors found. 
    inline bool glLogCall(const char* function, const char* file, int line) {
        while (GLenum error = glGetError()) {
            std::cout << "[OpenGL Error] (" << error << "):" << function << " " << file << ":" << line << std::endl;
            return false;
        }
        return true;
    }
}