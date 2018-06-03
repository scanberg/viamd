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

static GLint uniform_loc_model_view_proj_matrix = -1;
static GLint uniform_loc_tex_depth = -1;
static GLint uniform_loc_tex_volume = -1;
static GLint uniform_loc_color = -1;
static GLint uniform_loc_scale = -1;
static GLint uniform_loc_inv_res = -1;
static GLint uniform_loc_view_to_model_matrix = -1;
static GLint uniform_loc_model_to_tex_matrix = -1;
static GLint uniform_loc_inv_proj_matrix = -1;

static const char* v_shader_src_volume_renderer = R"(
#version 150 core

in vec3 in_pos;

uniform mat4 u_view_to_model_mat;
uniform mat4 u_model_view_proj_mat;

out vec3 model_pos;
out vec3 model_eye;
out vec3 color;

void main() {
	color = mix(vec3(1,0,0), vec3(0,1,0), float(gl_VertexID) / 42.0);
	model_pos = in_pos;
	model_eye = u_view_to_model_mat[3].xyz;
	gl_Position = u_model_view_proj_mat * vec4(in_pos, 1);
}
)";

static const char* f_shader_src_volume_renderer = R"(
#version 150 core

uniform sampler2D u_tex_depth;
uniform sampler3D u_tex_volume;
uniform vec3	  u_color;
uniform float	  u_scale;
uniform vec2	  u_inv_res;
uniform mat4      u_view_to_model_mat;
uniform mat4	  u_model_to_tex_mat;
uniform mat4	  u_inv_proj_mat;

in  vec3 model_pos;
in  vec3 model_eye;
in  vec3 color;

out vec4 out_frag;

// Modified version from source found here:
// https://gamedev.stackexchange.com/questions/18436/most-efficient-aabb-vs-ray-collision-algorithms#18459

bool ray_vs_aabb(out float t_entry, out float t_exit, in vec3 ori, in vec3 dir, in vec3 min_box, in vec3 max_box) {
	vec3 dir_frac = 1.0 / dir;

	vec3 tv_1 = (min_box - ori) * dir_frac;
	vec3 tv_2 = (max_box - ori) * dir_frac;

	vec3 tv_min = min(tv_1, tv_2);
	vec3 tv_max = max(tv_1, tv_2); 

	float t_min = max(max(tv_min.x, tv_min.y), tv_min.z);
	float t_max = min(min(tv_max.x, tv_max.y), tv_max.z);

	if (t_max < 0 || t_min > t_max) {
		return false;
	}

	t_entry = t_min;
	t_exit  = t_max;

	return true;
} 

vec4 depth_to_view_coord(vec2 tc, float depth) {
    vec4 clip_coord = vec4(vec3(tc, depth) * 2.0 - 1.0, 1.0);
    vec4 view_coord = u_inv_proj_mat * clip_coord;
    return view_coord / view_coord.w;
}

vec4 fetch_voxel(vec3 tc) {
	float a = min(texture(u_tex_volume, tc).x * u_scale, 1.0);
	return vec4(mix(vec3(0), u_color, min(a * 10.0, 1.0)), a);
}

const float REF = 150.0;
const float step = 0.005;

void main() {
	// Do everything in model space
	vec3 ori = model_eye;
    vec3 dir = normalize(model_pos - model_eye);

	float t_entry, t_exit;
	if (!ray_vs_aabb(t_entry, t_exit, ori, dir, vec3(0), vec3(1))) discard;
	t_entry = max(0, t_entry);

	float depth = texelFetch(u_tex_depth, ivec2(gl_FragCoord.xy), 0).x;
	if (depth < 1.0) {
		vec3 stop_pos = (u_view_to_model_mat * depth_to_view_coord(gl_FragCoord.xy * u_inv_res, depth)).xyz;
		t_exit = min(t_exit, dot(stop_pos - ori, dir));
	}

	vec4 result = vec4(0);
	for (float t = t_entry + step; t < t_exit; t += step) {
		vec3 pos  = ori + dir * t;
		vec4 rgba = fetch_voxel((u_model_to_tex_mat * vec4(pos, 1)).xyz);
		rgba.a    = 1.0 - pow(1.0 - rgba.a, step * REF);
		rgba.rgb *= rgba.a;

		result += (1.0 - result.a) * rgba;
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

    uniform_loc_model_view_proj_matrix = glGetUniformLocation(prog_volume_renderer, "u_model_view_proj_mat");
    uniform_loc_tex_depth = glGetUniformLocation(prog_volume_renderer, "u_tex_depth");
    uniform_loc_tex_volume = glGetUniformLocation(prog_volume_renderer, "u_tex_volume");
    uniform_loc_color = glGetUniformLocation(prog_volume_renderer, "u_color");
    uniform_loc_scale = glGetUniformLocation(prog_volume_renderer, "u_scale");
    uniform_loc_inv_res = glGetUniformLocation(prog_volume_renderer, "u_inv_res");
    uniform_loc_view_to_model_matrix = glGetUniformLocation(prog_volume_renderer, "u_view_to_model_mat");
    uniform_loc_model_to_tex_matrix = glGetUniformLocation(prog_volume_renderer, "u_model_to_tex_mat");
    uniform_loc_inv_proj_matrix = glGetUniformLocation(prog_volume_renderer, "u_inv_proj_mat");

    // From here:
    // https://stackoverflow.com/questions/28375338/cube-using-single-gl-triangle-strip
    static const float cube_strip[] = {0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1,
                                       0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0};

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

void create_volume_texture(GLuint* texture, const ivec3& dim) {
    ASSERT(texture);
    if (*texture != 0 && glIsTexture(*texture)) {
        glDeleteTextures(1, texture);
    }
    glGenTextures(1, texture);
    glBindTexture(GL_TEXTURE_3D, *texture);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, dim.x, dim.y, dim.z, 0, GL_RED, GL_FLOAT, nullptr);
    glBindTexture(GL_TEXTURE_3D, 0);
}

void free_volume_texture(GLuint texture) {
    if (glIsTexture(texture)) glDeleteTextures(1, &texture);
}

void set_volume_texture_data(GLuint texture, ivec3 dim, void* data) {
    if (glIsTexture(texture)) {
        glBindTexture(GL_TEXTURE_3D, texture);
        glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, dim.x, dim.y, dim.z, GL_RED, GL_FLOAT, data);
        // glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, volume.dim.x, volume.dim.y, volume.dim.z, 0, GL_RED, GL_FLOAT, volume.voxel_data.data);
        glBindTexture(GL_TEXTURE_3D, 0);
    }
}

mat4 compute_model_to_world_matrix(const vec3& min_world_aabb, const vec3& max_world_aabb) {
    vec3 ext = max_world_aabb - min_world_aabb;
    return mat4(vec4(ext.x, 0, 0, 0), vec4(0, ext.y, 0, 0), vec4(0, 0, ext.z, 0), vec4(min_world_aabb, 1));
}

mat4 compute_texture_to_model_matrix(const ivec3& dim) {
    return mat4(1);
    const vec3 cell_ext = 1.f / vec3(dim);
    const vec3 scl = 1.f + cell_ext;
    return math::inverse(mat4(vec4(scl.x, 0, 0, 0), vec4(0, scl.y, 0, 0), vec4(0, 0, scl.z, 0), vec4(-0.5f * cell_ext, 1.f)));
}

void save_volume_to_file(const Volume& volume, const char* file) {
    FILE* fs = fopen(file, "wb");
    if (!fs) {
        LOG_ERROR("Could not open file %s", fs);
    }
    fwrite(volume.voxel_data.data, sizeof(Volume::VoxelDataType), volume.voxel_data.count, fs);
    fclose(fs);
}

void render_volume_texture(GLuint volume_texture, GLuint depth_texture, const mat4& texture_matrix, const mat4& model_matrix, const mat4& view_matrix,
                           const mat4& proj_matrix, vec3 color, float opacity_scale) {
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    const mat4 model_to_view_matrix = view_matrix * model_matrix;
    const mat4 view_to_model_matrix = math::inverse(model_to_view_matrix);
    const mat4 model_view_proj_matrix = proj_matrix * model_to_view_matrix;
    const mat4 model_to_tex_matrix = math::inverse(texture_matrix);
    const mat4 inv_proj_matrix = math::inverse(proj_matrix);
    const vec2 inv_res = vec2(1.f / (float)(viewport[2]), 1.f / (float)(viewport[3]));

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glDepthMask(GL_TRUE);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, depth_texture);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_3D, volume_texture);

    glUseProgram(prog_volume_renderer);

    glUniform1i(uniform_loc_tex_depth, 0);
    glUniform1i(uniform_loc_tex_volume, 1);
    glUniform3fv(uniform_loc_color, 1, &color[0]);
    glUniform1f(uniform_loc_scale, opacity_scale);
    glUniform2fv(uniform_loc_inv_res, 1, &inv_res[0]);
    glUniformMatrix4fv(uniform_loc_model_view_proj_matrix, 1, GL_FALSE, &model_view_proj_matrix[0][0]);
    glUniformMatrix4fv(uniform_loc_view_to_model_matrix, 1, GL_FALSE, &view_to_model_matrix[0][0]);
    glUniformMatrix4fv(uniform_loc_model_to_tex_matrix, 1, GL_FALSE, &model_to_tex_matrix[0][0]);
    glUniformMatrix4fv(uniform_loc_inv_proj_matrix, 1, GL_FALSE, &inv_proj_matrix[0][0]);

    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 42);
    glBindVertexArray(0);

    glUseProgram(0);

    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
}

}  // namespace volume
