// GPU - Low-level GPU API abstraction
// Handles textures, shaders, framebuffers, and other GPU resources

#pragma once

#include "gl.h"
#include <core/md_vec_math.h>
#include <stdint.h>

namespace gpu {

// Texture management
struct Texture {
    uint32_t id = 0;
    int width = 0;
    int height = 0;
    int depth = 0;  // For 3D textures
};

struct TextureDesc {
    int width = 0;
    int height = 0;
    int depth = 1;  // 1 for 2D, >1 for 3D
    GLenum internal_format = GL_RGBA8;
    GLenum format = GL_RGBA;
    GLenum type = GL_UNSIGNED_BYTE;
    GLenum min_filter = GL_LINEAR;
    GLenum mag_filter = GL_LINEAR;
    GLenum wrap_s = GL_CLAMP_TO_EDGE;
    GLenum wrap_t = GL_CLAMP_TO_EDGE;
    GLenum wrap_r = GL_CLAMP_TO_EDGE;
    const void* data = nullptr;
};

Texture create_texture_2d(const TextureDesc& desc);
Texture create_texture_3d(const TextureDesc& desc);
void update_texture_2d(Texture& tex, const void* data, int width, int height);
void update_texture_3d(Texture& tex, const void* data, int width, int height, int depth);
void destroy_texture(Texture& tex);
void bind_texture(const Texture& tex, uint32_t slot = 0);
void unbind_texture(uint32_t slot = 0);

// Shader management
struct Shader {
    uint32_t program = 0;
};

struct ShaderSource {
    const char* vertex = nullptr;
    const char* fragment = nullptr;
    const char* geometry = nullptr;
};

Shader create_shader(const ShaderSource& source);
void destroy_shader(Shader& shader);
void use_shader(const Shader& shader);
void set_uniform(const Shader& shader, const char* name, int value);
void set_uniform(const Shader& shader, const char* name, float value);
void set_uniform(const Shader& shader, const char* name, const vec2_t& value);
void set_uniform(const Shader& shader, const char* name, const vec3_t& value);
void set_uniform(const Shader& shader, const char* name, const vec4_t& value);
void set_uniform(const Shader& shader, const char* name, const mat4_t& value);

// Framebuffer management
struct Framebuffer {
    uint32_t fbo = 0;
    uint32_t width = 0;
    uint32_t height = 0;
};

Framebuffer create_framebuffer(int width, int height);
void destroy_framebuffer(Framebuffer& fb);
void bind_framebuffer(const Framebuffer& fb);
void unbind_framebuffer();
void attach_texture_to_framebuffer(const Framebuffer& fb, const Texture& tex, GLenum attachment);
void set_draw_buffers(const GLenum* attachments, int count);

// Vertex Array Object management
struct VertexArray {
    uint32_t vao = 0;
};

VertexArray create_vertex_array();
void destroy_vertex_array(VertexArray& vao);
void bind_vertex_array(const VertexArray& vao);
void unbind_vertex_array();

// Rendering state
void set_viewport(int x, int y, int width, int height);
void set_scissor(int x, int y, int width, int height);
void enable_blend(GLenum src = GL_SRC_ALPHA, GLenum dst = GL_ONE_MINUS_SRC_ALPHA);
void disable_blend();
void enable_depth_test(GLenum func = GL_LESS);
void disable_depth_test();
void enable_cull_face(GLenum mode = GL_BACK);
void disable_cull_face();
void enable_scissor_test();
void disable_scissor_test();
void set_depth_mask(bool enable);
void set_color_mask(bool r, bool g, bool b, bool a);
void clear_color(float r, float g, float b, float a);
void clear(GLbitfield mask);

// Drawing
void draw_arrays(GLenum mode, int first, int count);
void draw_elements(GLenum mode, int count, GLenum type, const void* indices);

} // namespace gpu
