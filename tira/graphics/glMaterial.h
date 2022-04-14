#include "tira/graphics/glShader.h"
#include "tira/graphics/glTexture.h"
#include "tira/image/image.h"

namespace tira{

    class glMaterial : public glShader {

        

        struct glTextureUnit {
            glTexture texture;
            glShaderUniform uniform;
        };
        // Unordered map contains information for the texture uniform (in the shader) and texture map (on the GPU)
        // keyed to the uniform name. This will be updated anytime a new shader is loaded, and the textures will
        // only be destroyed if they no longer exist as variables in the shader.
        std::unordered_map<std::string, glTextureUnit> m_TextureUnits;

        bool isSampler(GLenum uniform_type){
            switch(uniform_type){
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
                    return true;
                default:
                    return false;
            };
        }

        /// Create and return a vector containing the uniform details about every texture sampler in the shader
        std::vector<glShaderUniform> GetSamplerList(){
            std::vector<glShaderUniform> sampler_list;

            // for each uniform in the shader
            for(std::unordered_map<std::string, glShaderUniform>::iterator i = m_UniformCache.begin();
                i != m_UniformCache.end(); i++){
                    if(isSampler((*i).second.type))             // if the uniform is a texture sampler
                        sampler_list.push_back((*i).second);    // store it in the texture list
            }
            return sampler_list;                                // return a list of all texture sampler uniforms
        }

        /// Create a list of appropriate textures for all samplers in the shader. If textures already exist,
        /// remove any that are unused.
        void RefreshTextures(){
            std::unordered_map<std::string, glTextureUnit> new_texture_units;

            std::vector<glShaderUniform> sampler_list = GetSamplerList();   // get a list of all samplers in the shader
            glUseProgram(m_ShaderID);                                       // bund the current shader program to add textures
            for(size_t i = 0; i < sampler_list.size(); i++){                // for each sampler in the shader
                // look for an existing texture unit with the sampler name
                std::unordered_map<std::string, glTextureUnit>::iterator u = m_TextureUnits.find(sampler_list[i].name);
                if(u == m_TextureUnits.end()){                              // if the name doesn't exist
                    glTextureUnit tu;
                    tu.uniform = sampler_list[i];
                    new_texture_units[sampler_list[i].name] = tu;
                }
                else{                                                       // if the texture unit already exists
                    new_texture_units[sampler_list[i].name] = (*u).second;  // copy it to the new texture unit map
                    m_TextureUnits.erase(u);                                // erase the old copy of the texture unit
                }
                glUniform1i(sampler_list[i].location, i);                   // assign the texture unit to the uniform
            }

            //erase any remaining texture units that are unused in the current shader
            for(std::unordered_map<std::string, glTextureUnit>::iterator i = m_TextureUnits.begin();
                i != m_TextureUnits.end(); i++){
                    (*i).second.texture.Destroy();
            }
            m_TextureUnits.clear();
            m_TextureUnits = new_texture_units;
        }

        public:
        glMaterial() : glShader(){
            
        }

        /// <summary>
        /// Constructor to create a new material based on a shader source file.
        /// </summary>
        /// <param name="filename">Source code containing the vertex and fragment shaders</param>
        glMaterial(std::string filename) : glShader(filename){
            RefreshTextures();                              // identify all texture samplers so that textures can be loaded separately
        }

        /// <summary>
        /// Load a shader from a single source file separated by "# shader vertex" and "# shader fragment" blocks.
        /// </summary>
        /// <param name="filepath">Name of the source file.</param>
        void LoadShader(const std::string& filepath) {
            glShader::LoadShader(filepath);                 // call the glShader function to load a new shader
            RefreshTextures();                              // identify all texture samplers so that textures can be loaded separately
        }

        /// <summary>
        /// Create a shader program by compiling source code for a vertex and fragment shader.
        /// </summary>
        /// <param name="vertex_source">String containing the vertex shader source code.</param>
        /// <param name="fragment_source">String containing the fragment shader source code.</param>
        void CreateShader(std::string vertex_source, std::string fragment_source) {
            glShader::CreateShader(vertex_source, fragment_source);
            RefreshTextures();
        }

        /// <summary>
        /// Begin using the material (all geometry rendered between Begin() and End() will use this material).
        /// </summary>
        void Begin(){
            int i = 0;
            for(std::unordered_map<std::string, glTextureUnit>::iterator tu = m_TextureUnits.begin();
                tu != m_TextureUnits.end(); tu++){
                    glActiveTexture(GL_TEXTURE0 + i);
                    (*tu).second.texture.Bind();
                    i++;
            }
            Bind();
        }

        /// <summary>
        /// End using the material (all geometry rendered between Begin() and End() will use this material).
        /// </summary>
        void End(){
            Unbind();
        }

        /// <summary>
        /// Assign a texture map to an existing uniform variable.
        /// </summary>
        /// <typeparam name="T">Data type used to represent the input image.</typeparam>
        /// <param name="name">Uniform variable name for the associated texture sampler.</param>
        /// <param name="teximage">Image object containing the texture image.</param>
        /// <param name="internalFormat">Format used to represent the texture internally on the GPU.</param>
        /// <param name="filter_type">Filter type used for texture sampling (GL_LINEAR, GL_NEAREST).</param>
        template <typename T>
        void SetTexture(std::string name, tira::image<T> teximage, GLenum internalFormat = GL_RGB, GLenum filter_type = GL_LINEAR) {
            std::unordered_map<std::string, glTextureUnit>::iterator i = m_TextureUnits.find(name); // create an iterator
            if(i == m_TextureUnits.end()) return;               // if the texture name doesn't exist, return

            GLenum externalFormat;
            size_t numChannels = teximage.channels();
            if (numChannels == 1) externalFormat = GL_RED;
            else if (numChannels == 2) externalFormat = GL_RG;
            else if (numChannels == 3) externalFormat = GL_RGB;
            else if (numChannels == 4) externalFormat = GL_RGBA;
            else externalFormat = GL_RGBA;
            (*i).second.texture.AssignImage(teximage.data(),
                                            teximage.width(),
                                            teximage.height(),
                                            0,
                                            internalFormat,
                                            externalFormat,
                                            GL_UNSIGNED_BYTE);
            (*i).second.texture.SetFilter(filter_type);
        }

        /// <summary>
        /// Return a glTexture object representing a texture with a given sampler name
        /// </summary>
        /// <param name="sampler_name"></param>
        /// <returns></returns>
        glTexture GetTexture(std::string sampler_name) {
            return m_TextureUnits[sampler_name].texture;
        }
    };
}