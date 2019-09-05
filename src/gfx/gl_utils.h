#pragma once

#include <core/string_types.h>
#include "gl.h"

namespace gl {

bool get_shader_compile_error(char* buffer, int max_length, GLuint shader);
bool get_program_link_error(char* buffer, int max_length, GLuint program);
GLuint compile_shader_from_file(CString filename, GLenum shader_type);
bool attach_link_detach(GLuint program, Array<const GLuint> shaders);
bool attach_link_detach_with_transform_feedback(GLuint program, Array<const GLuint> shaders, Array<const char*> varyings, GLenum buffer_capture_mode);

}  // namespace gl
