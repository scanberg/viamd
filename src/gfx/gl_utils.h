#pragma once

#include <core/gl.h>

namespace gl {

bool get_shader_compile_error(char* buffer, int max_length, GLuint shader);
bool get_program_link_error(char* buffer, int max_length, GLuint program);

}  // namespace gl
