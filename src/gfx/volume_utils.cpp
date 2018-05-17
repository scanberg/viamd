#include "volume_utils.h"
#include <core/common.h>
#include <core/log.h>
#include <core/math_utils.h>
#include <gfx/gl_utils.h>
#include <stdio.h>

namespace volume {

static GLuint vao = 0;
static GLuint vbo = 0;

static GLuint prog_volume_renderer = 0;

static GLint attrib_loc_pos = -1;

static GLint uniform_loc_view_matrix = -1;
static GLint uniform_loc_proj_matrix = -1;
static GLint uniform_loc_tex_depth = -1;
static GLint uniform_loc_tex_volume = -1;
static GLint uniform_loc_color = -1;
static GLint uniform_loc_scale = -1;
static GLint uniform_loc_view_to_vol_matrix = -1;

static const char* v_shader_src_volume_renderer = R"(
#version 150 core

in vec3 in_pos;

uniform mat4 u_view_mat;
uniform mat4 u_proj_mat;

out vec3 view_pos;

void main() {
	view_pos = (u_view_mat * vec4(in_pos)).xyz;
	gl_Position = u_proj_mat * vec4(view_pos, 1);
}
)";

static const char* f_shader_src_volume_renderer = R"(
#version 150 core

uniform sampler2D u_tex_depth;
uniform sampler3D u_tex_volume;
uniform vec4	  u_color;
uniform float	  u_scale;
uniform mat4	  u_view_to_vol_mat;

in  vec3 view_pos;
out vec4 out_frag;

vec4 fetch_voxel(vec3 tc) {
	float a = texture(u_tex_volume, tc).x * u_scale;
	return vec4(u_color.xyz, a);
}

const float REF = 150.0;

void main() {

	out_frag = vec4(1,0,0,1);
	return;

	vec3 p = view_pos;
	vec3 d = normalize(view_pos);

	vec4 result = vec4(0);
	
	for (float d = 0; d < 100.0; d += step) {
		p += d * step;
		vec3 tc	= view_to_tex_mat * vec4(p, 1);
		if (any(lessThan(tc, vec3(0))) || any(greaterThan(tc, vec3(1)))) break;
		vec4 rgba = fetch_voxel(tc);
		
		rgba.a     = 1.0 - pow(1.0 - rgba.a, step * REF);
		result.rgb = result.rgb + (1.0 - result.a) * rgba.a * rgba.rgb;
		result.a   = result.a + (1.0 - result.a) * rgba.a;
		if (result.a > 0.99) break;
	}

	out_frag = result;
}
)";

void initialize() {
    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    GLuint v_shader = glCreateShader(GL_VERTEX_SHADER);
    GLuint f_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(v_shader, 1, &v_shader_src_volume_renderer, 0);
    glShaderSource(f_shader, 1, &f_shader_src_volume_renderer, 0);
    glCompileShader(v_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, v_shader)) {
        LOG_ERROR("Error while compiling volume renderer vertex shader:\n%s\n", buffer);
    }
    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        LOG_ERROR("Error while compiling volume renderer fragment shader:\n%s\n", buffer);
    }

    prog_volume_renderer = glCreateProgram();
    glAttachShader(prog_volume_renderer, v_shader);
    glAttachShader(prog_volume_renderer, f_shader);
    glLinkProgram(prog_volume_renderer);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, prog_volume_renderer)) {
        LOG_ERROR("Error while linking volume renderer program:\n%s\n", buffer);
    }

    glDetachShader(prog_volume_renderer, v_shader);
    glDetachShader(prog_volume_renderer, f_shader);

    glDeleteShader(v_shader);
    glDeleteShader(f_shader);

    attrib_loc_pos = glGetAttribLocation(prog_volume_renderer, "in_pos");

    uniform_loc_view_matrix = glGetUniformLocation(prog_volume_renderer, "u_view_mat");
    uniform_loc_proj_matrix = glGetUniformLocation(prog_volume_renderer, "u_proj_mat");
    uniform_loc_tex_depth = glGetUniformLocation(prog_volume_renderer, "u_tex_depth");
    uniform_loc_tex_volume = glGetUniformLocation(prog_volume_renderer, "u_tex_volume");
    uniform_loc_color = glGetUniformLocation(prog_volume_renderer, "u_color");
    uniform_loc_scale = glGetUniformLocation(prog_volume_renderer, "u_scale");
    uniform_loc_view_to_vol_matrix = glGetUniformLocation(prog_volume_renderer, "u_view_to_vol_mat");

    // From here:
    // https://stackoverflow.com/questions/28375338/cube-using-single-gl-triangle-strip
    static const GLfloat cube_strip[] = {
        0.f, 1.f, 1.f,  // Front-top-left
        1.f, 1.f, 1.f,  // Front-top-right
        0.f, 0.f, 1.f,  // Front-bottom-left
        1.f, 0.f, 1.f,  // Front-bottom-right
        1.f, 0.f, 0.f,  // Back-bottom-right
        1.f, 1.f, 1.f,  // Front-top-right
        1.f, 1.f, 0.f,  // Back-top-right
        0.f, 1.f, 1.f,  // Front-top-left
        0.f, 1.f, 0.f,  // Back-top-left
        0.f, 0.f, 1.f,  // Front-bottom-left
        0.f, 0.f, 0.f,  // Back-bottom-left
        1.f, 0.f, 0.f,  // Back-bottom-right
        0.f, 1.f, 0.f,  // Back-top-left
        1.f, 1.f, 0.f   // Back-top-right
    };

    if (!vbo) {
        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(cube_strip), cube_strip, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    if (!vao) {
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        if (attrib_loc_pos != -1) {
            glEnableVertexAttribArray(attrib_loc_pos);
            glVertexAttribPointer(attrib_loc_pos, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid*)0);
        }
        glBindVertexArray(0);
    }
}

void shutdown() {}

GLuint create_volume_texture() {
    GLuint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_3D, tex);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_3D, 0);
    return tex;
}

void free_volume_texture(GLuint texture) {
    if (glIsTexture(texture)) glDeleteTextures(1, &texture);
}

void set_volume_texture_data(GLuint texture, const Volume& volume) {
    if (glIsTexture(texture)) {
        glTexImage3D(GL_TEXTURE_3D, 0, GL_R8, volume.dim.x, volume.dim.y, volume.dim.z, 0, GL_RED, GL_FLOAT, volume.voxel_data.data);
    }
}

void render_volume_texture(GLuint volume_texture, GLuint depth_texture, const mat4& basis, const mat4& view_matrix, const mat4& proj_matrix,
                           vec3 color, float opacity_scale) {

    const mat4 view_to_vol_matrix = basis * math::inverse(view_matrix);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, depth_texture);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_3D, volume_texture);

    glUseProgram(prog_volume_renderer);

    glUniform1i(uniform_loc_tex_depth, 0);
    glUniform1i(uniform_loc_tex_volume, 1);
    glUniform3fv(uniform_loc_color, 1, &color[0]);
    glUniform1f(uniform_loc_scale, opacity_scale);
    glUniformMatrix4fv(uniform_loc_view_matrix, 1, GL_FALSE, &view_matrix[0][0]);
    glUniformMatrix4fv(uniform_loc_proj_matrix, 1, GL_FALSE, &proj_matrix[0][0]);
    glUniformMatrix4fv(uniform_loc_view_to_vol_matrix, 1, GL_FALSE, &view_to_vol_matrix[0][0]);

    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 42);
    glBindVertexArray(0);
}

}  // namespace volume
