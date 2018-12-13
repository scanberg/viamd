#include "gl_utils.h"

bool gl::get_shader_compile_error(char * buffer, int max_length, GLuint shader) {
	GLint success = 0;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
	if (success == GL_FALSE) {
		int length;
		glGetShaderInfoLog(shader, max_length, &length, buffer);
		return true;
	}
	else {
		return false;
	}
}

bool gl::get_program_link_error(char * buffer, int max_length, GLuint program) {
	GLint success = 0;
	glGetProgramiv(program, GL_LINK_STATUS, &success);
	if (success == GL_FALSE) {
		int length;
		glGetProgramInfoLog(program, max_length, &length, buffer);
		return true;
	}
	else {
		return false;
	}
}
