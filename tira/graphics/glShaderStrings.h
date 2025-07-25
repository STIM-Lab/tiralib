namespace tira {
    /**
     * Data structure containing a series of string constants that can be used to quickly assemble basic
     * shaders. The string data is loaded from individual .shader source files written in GLSL, and any
     * loaded strings will be directly compiled into the source code and do not need to be loaded independently
     * during runtime.
     */
    struct glShaderStrings {
        /**
         * @brief Basic vertex and fragment shader program that requires a view and projection matrix and renders all pixels as the specified color.
         * @param V (uniform mat4) is a view matrix used to transform vertices
         * @param P (uniform mat4) is a projection matrix used to transform the vertices
         * @param c (uniform vec4) is a color used to render all fragments (default = [1.0, 0.0, 1.0, 1.0])
         */
        inline static const std::string vf_basic =
        #include "shaders/vf_basic.shader"
        ;

        inline static const std::string vf_brewer =
        #include "shaders/vf_brewer.shader"
        ;
    };
}



