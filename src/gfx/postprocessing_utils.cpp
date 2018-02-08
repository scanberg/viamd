#include "postprocessing_utils.h"

namespace postprocessing {

// This is for the fullscreen quad
GLuint vao;
GLuint vbo;
GLuint v_shader_fs_quad;

static const char* v_shader_src_fs_quad = R"(
#version 150 core

out vec2 tc;

void main() {
	uint idx = uint(gl_VertexID) % 3U;
	gl_Position = vec4(
		(float( idx     &1U)) * 4.0 - 1.0,
		(float((idx>>1U)&1U)) * 4.0 - 1.0,
		0, 1.0);
	tc = gl_Position.xy * 0.5 + 0.5;
}
)";

namespace ssao {
GLuint fbo_depth_linear;
GLuint fbo_hbao;

GLuint tex_random;
GLuint tex_depth_linear;
GLuint tex_blur;
GLuint tex_result;

GLuint prog_depth_linearize;
GLuint prog_hbao;
GLuint prog_blur;
GLuint prog_blur_vert;
GLuint prog_blur_horiz;

static const char* f_shader_src = R"(
#version 150 core

in vec2 tc;
in vec4 col;

out vec4 out_frag;

void main() {
	out_frag = col;
}
)";

void initialize() {
	/*
    v_shader = glCreateShader(GL_VERTEX_SHADER);
    f_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(v_shader, 1, &v_shader_src, 0);
    glShaderSource(f_shader, 1, &f_shader_src, 0);
    glCompileShader(v_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, v_shader)) {
        printf("Error while compiling immediate vertex shader:\n%s\n", buffer);
    }
    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        printf("Error while compiling immediate fragment shader:\n%s\n", buffer);
    }

    program = glCreateProgram();
    glAttachShader(program, v_shader);
    glAttachShader(program, f_shader);
    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, program)) {
        printf("Error while linking immediate program:\n%s\n", buffer);
    }

    glDetachShader(program, v_shader);
    glDetachShader(program, f_shader);

    glDeleteShader(v_shader);
    glDeleteShader(f_shader);

    attrib_loc_pos = glGetAttribLocation(program, "in_pos");
    attrib_loc_tc = glGetAttribLocation(program, "in_tc");
    attrib_loc_col = glGetAttribLocation(program, "in_col");
    uniform_loc_mvp = glGetUniformLocation(program, "u_mvp");

    glGenBuffers(1, &vbo);
    glGenBuffers(1, &ibo);

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)0);

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_UNSIGNED_SHORT, GL_TRUE, sizeof(Vertex), (const GLvoid*)12);

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (const GLvoid*)16);

    glBindVertexArray(0);

    vertices.reserve(100000);
	*/
}

void shutdown() {
	/*
    if (vbo) glDeleteBuffers(1, &vbo);
    if (ibo) glDeleteBuffers(1, &ibo);
    if (vao) glDeleteVertexArrays(1, &vao);
    if (program) glDeleteProgram(program);
	*/
}

}  // namespace ssao

void ssao(GLuint depth_tex, float strength);

}  // namespace postprocessing