#pragma once

#include "gl.h"
#include <core/md_str.h>

namespace gl {

bool get_shader_compile_error(char* buffer, int max_length, GLuint shader);
bool get_program_link_error(char* buffer, int max_length, GLuint program);
GLuint compile_shader_from_source(md_str_t source, GLenum shader_type, md_str_t defines = {}, md_str_t base_include_dir = {});
bool attach_link_detach(GLuint program, const GLuint shaders[], int num_shaders);

// A 2D texture of the size and format, made again when it is of another, linearly filtered and clamped
bool init_texture_2D(GLuint* texture, int width, int height, GLenum format);
bool free_texture(GLuint* texture);
// The data of all of a 2D texture (RGBA8, RGBA32F, R32F or RGB10_A2)
bool set_texture_2D_data(GLuint texture, int level, const void* data, GLenum format);

}  // namespace gl
