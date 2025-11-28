// GPU - Low-level GPU API implementation

#include "gpu.h"
#include "gl_utils.h"
#include <core/md_log.h>

namespace gpu {

// Texture management
Texture create_texture_2d(const TextureDesc& desc) {
    Texture tex;
    tex.width = desc.width;
    tex.height = desc.height;
    tex.depth = 1;
    
    glGenTextures(1, &tex.id);
    glBindTexture(GL_TEXTURE_2D, tex.id);
    
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, desc.min_filter);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, desc.mag_filter);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, desc.wrap_s);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, desc.wrap_t);
    
    // Use glTexStorage2D for immutable storage if no data provided, otherwise use glTexImage2D
    if (desc.data == nullptr) {
        glTexStorage2D(GL_TEXTURE_2D, 1, desc.internal_format, desc.width, desc.height);
    } else {
        glTexImage2D(GL_TEXTURE_2D, 0, desc.internal_format, desc.width, desc.height, 0, 
                     desc.format, desc.type, desc.data);
    }
    
    glBindTexture(GL_TEXTURE_2D, 0);
    
    return tex;
}

Texture create_texture_3d(const TextureDesc& desc) {
    Texture tex;
    tex.width = desc.width;
    tex.height = desc.height;
    tex.depth = desc.depth;
    
    glGenTextures(1, &tex.id);
    glBindTexture(GL_TEXTURE_3D, tex.id);
    
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, desc.min_filter);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, desc.mag_filter);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, desc.wrap_s);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, desc.wrap_t);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, desc.wrap_r);
    
    glTexImage3D(GL_TEXTURE_3D, 0, desc.internal_format, desc.width, desc.height, desc.depth, 0,
                 desc.format, desc.type, desc.data);
    
    glBindTexture(GL_TEXTURE_3D, 0);
    
    return tex;
}

void update_texture_2d(Texture& tex, const void* data, int width, int height) {
    if (tex.id == 0) return;
    
    glBindTexture(GL_TEXTURE_2D, tex.id);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, data);
    glBindTexture(GL_TEXTURE_2D, 0);
    
    tex.width = width;
    tex.height = height;
}

void update_texture_3d(Texture& tex, const void* data, int width, int height, int depth) {
    if (tex.id == 0) return;
    
    glBindTexture(GL_TEXTURE_3D, tex.id);
    glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, width, height, depth, GL_RGBA, GL_UNSIGNED_BYTE, data);
    glBindTexture(GL_TEXTURE_3D, 0);
    
    tex.width = width;
    tex.height = height;
    tex.depth = depth;
}

void destroy_texture(Texture& tex) {
    if (tex.id != 0) {
        glDeleteTextures(1, &tex.id);
        tex.id = 0;
        tex.width = 0;
        tex.height = 0;
        tex.depth = 0;
    }
}

void bind_texture(const Texture& tex, uint32_t slot) {
    glActiveTexture(GL_TEXTURE0 + slot);
    if (tex.depth > 1) {
        glBindTexture(GL_TEXTURE_3D, tex.id);
    } else {
        glBindTexture(GL_TEXTURE_2D, tex.id);
    }
}

void unbind_texture(uint32_t slot) {
    glActiveTexture(GL_TEXTURE0 + slot);
    glBindTexture(GL_TEXTURE_2D, 0);
    glBindTexture(GL_TEXTURE_3D, 0);
}

// Shader management
Shader create_shader(const ShaderSource& source) {
    Shader shader;
    shader.program = glCreateProgram();
    
    // Compile vertex shader
    if (source.vertex) {
        uint32_t vs = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vs, 1, &source.vertex, nullptr);
        glCompileShader(vs);
        
        // Check compilation
        int success;
        glGetShaderiv(vs, GL_COMPILE_STATUS, &success);
        if (!success) {
            char log[512];
            glGetShaderInfoLog(vs, 512, nullptr, log);
            MD_LOG_ERROR("Vertex shader compilation failed: %s", log);
        }
        glAttachShader(shader.program, vs);
        glDeleteShader(vs);
    }
    
    // Compile fragment shader
    if (source.fragment) {
        uint32_t fs = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fs, 1, &source.fragment, nullptr);
        glCompileShader(fs);
        
        // Check compilation
        int success;
        glGetShaderiv(fs, GL_COMPILE_STATUS, &success);
        if (!success) {
            char log[512];
            glGetShaderInfoLog(fs, 512, nullptr, log);
            MD_LOG_ERROR("Fragment shader compilation failed: %s", log);
        }
        glAttachShader(shader.program, fs);
        glDeleteShader(fs);
    }
    
    // Compile geometry shader
    if (source.geometry) {
        uint32_t gs = glCreateShader(GL_GEOMETRY_SHADER);
        glShaderSource(gs, 1, &source.geometry, nullptr);
        glCompileShader(gs);
        
        // Check compilation
        int success;
        glGetShaderiv(gs, GL_COMPILE_STATUS, &success);
        if (!success) {
            char log[512];
            glGetShaderInfoLog(gs, 512, nullptr, log);
            MD_LOG_ERROR("Geometry shader compilation failed: %s", log);
        }
        glAttachShader(shader.program, gs);
        glDeleteShader(gs);
    }
    
    // Link program
    glLinkProgram(shader.program);
    
    // Check linking
    int success;
    glGetProgramiv(shader.program, GL_LINK_STATUS, &success);
    if (!success) {
        char log[512];
        glGetProgramInfoLog(shader.program, 512, nullptr, log);
        MD_LOG_ERROR("Shader program linking failed: %s", log);
    }
    
    return shader;
}

void destroy_shader(Shader& shader) {
    if (shader.program != 0) {
        glDeleteProgram(shader.program);
        shader.program = 0;
    }
}

void use_shader(const Shader& shader) {
    glUseProgram(shader.program);
}

void set_uniform(const Shader& shader, const char* name, int value) {
    int loc = glGetUniformLocation(shader.program, name);
    if (loc != -1) {
        glUniform1i(loc, value);
    }
}

void set_uniform(const Shader& shader, const char* name, float value) {
    int loc = glGetUniformLocation(shader.program, name);
    if (loc != -1) {
        glUniform1f(loc, value);
    }
}

void set_uniform(const Shader& shader, const char* name, const vec2_t& value) {
    int loc = glGetUniformLocation(shader.program, name);
    if (loc != -1) {
        glUniform2f(loc, value.x, value.y);
    }
}

void set_uniform(const Shader& shader, const char* name, const vec3_t& value) {
    int loc = glGetUniformLocation(shader.program, name);
    if (loc != -1) {
        glUniform3f(loc, value.x, value.y, value.z);
    }
}

void set_uniform(const Shader& shader, const char* name, const vec4_t& value) {
    int loc = glGetUniformLocation(shader.program, name);
    if (loc != -1) {
        glUniform4f(loc, value.x, value.y, value.z, value.w);
    }
}

void set_uniform(const Shader& shader, const char* name, const mat4_t& value) {
    int loc = glGetUniformLocation(shader.program, name);
    if (loc != -1) {
        glUniformMatrix4fv(loc, 1, GL_FALSE, &value.elem[0][0]);
    }
}

// Framebuffer management
Framebuffer create_framebuffer(int width, int height) {
    Framebuffer fb;
    fb.width = width;
    fb.height = height;
    
    glGenFramebuffers(1, &fb.fbo);
    
    return fb;
}

void destroy_framebuffer(Framebuffer& fb) {
    if (fb.fbo != 0) {
        glDeleteFramebuffers(1, &fb.fbo);
        fb.fbo = 0;
        fb.width = 0;
        fb.height = 0;
    }
}

void bind_framebuffer(const Framebuffer& fb) {
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fb.fbo);
}

void unbind_framebuffer() {
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
}

void attach_texture_to_framebuffer(const Framebuffer& fb, const Texture& tex, GLenum attachment) {
    glBindFramebuffer(GL_FRAMEBUFFER, fb.fbo);
    if (tex.depth > 1) {
        glFramebufferTexture(GL_FRAMEBUFFER, attachment, tex.id, 0);
    } else {
        glFramebufferTexture2D(GL_FRAMEBUFFER, attachment, GL_TEXTURE_2D, tex.id, 0);
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void set_draw_buffers(const GLenum* attachments, int count) {
    glDrawBuffers(count, attachments);
}

// Vertex Array Object management
VertexArray create_vertex_array() {
    VertexArray vao;
    glGenVertexArrays(1, &vao.vao);
    return vao;
}

void destroy_vertex_array(VertexArray& vao) {
    if (vao.vao != 0) {
        glDeleteVertexArrays(1, &vao.vao);
        vao.vao = 0;
    }
}

void bind_vertex_array(const VertexArray& vao) {
    glBindVertexArray(vao.vao);
}

void unbind_vertex_array() {
    glBindVertexArray(0);
}

// Rendering state
void set_viewport(int x, int y, int width, int height) {
    glViewport(x, y, width, height);
}

void set_scissor(int x, int y, int width, int height) {
    glScissor(x, y, width, height);
}

void enable_blend(GLenum src, GLenum dst) {
    glEnable(GL_BLEND);
    glBlendFunc(src, dst);
}

void disable_blend() {
    glDisable(GL_BLEND);
}

void enable_depth_test(GLenum func) {
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(func);
}

void disable_depth_test() {
    glDisable(GL_DEPTH_TEST);
}

void enable_cull_face(GLenum mode) {
    glEnable(GL_CULL_FACE);
    glCullFace(mode);
}

void disable_cull_face() {
    glDisable(GL_CULL_FACE);
}

void enable_scissor_test() {
    glEnable(GL_SCISSOR_TEST);
}

void disable_scissor_test() {
    glDisable(GL_SCISSOR_TEST);
}

void set_depth_mask(bool enable) {
    glDepthMask(enable ? GL_TRUE : GL_FALSE);
}

void set_color_mask(bool r, bool g, bool b, bool a) {
    glColorMask(r ? GL_TRUE : GL_FALSE, g ? GL_TRUE : GL_FALSE, 
                b ? GL_TRUE : GL_FALSE, a ? GL_TRUE : GL_FALSE);
}

void clear_color(float r, float g, float b, float a) {
    glClearColor(r, g, b, a);
}

void clear(GLbitfield mask) {
    glClear(mask);
}

// Drawing
void draw_arrays(GLenum mode, int first, int count) {
    glDrawArrays(mode, first, count);
}

void draw_elements(GLenum mode, int count, GLenum type, const void* indices) {
    glDrawElements(mode, count, type, indices);
}

} // namespace gpu
