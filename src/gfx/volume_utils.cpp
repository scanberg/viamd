#include "volume_utils.h"
#include <core/common.h>
#include <core/log.h>
#include <core/math_utils.h>
#include <gfx/gl_utils.h>
#include <stdio.h>

#include <fstream>

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
static GLint uniform_loc_inv_res = -1;
static GLint uniform_loc_view_to_vol_matrix = -1;
static GLint uniform_loc_inv_proj_matrix = -1;

static const char* v_shader_src_volume_renderer = R"(
#version 150 core

in vec3 in_pos;

uniform mat4 u_view_mat;
uniform mat4 u_proj_mat;

out vec3 view_pos;

void main() {
	view_pos = (u_view_mat * vec4(in_pos, 1)).xyz;
	gl_Position = u_proj_mat * vec4(view_pos, 1);
}
)";

static const char* f_shader_src_volume_renderer = R"(
#version 150 core

uniform sampler2D u_tex_depth;
uniform sampler3D u_tex_volume;
uniform vec3	  u_color;
uniform float	  u_scale;
uniform vec2	  u_inv_res;
uniform mat4	  u_view_to_vol_mat;
uniform mat4	  u_inv_proj_mat;

in  vec3 view_pos;
out vec4 out_frag;

/*
float ray_vs_aabb(vec3 o, vec3 d, vec3 min_box, vec3 max_box) {
	vec3 dirfrac = 1.f / 
	dirfrac.x = 1.0f / r.dir.x;
	dirfrac.y = 1.0f / r.dir.y;
	dirfrac.z = 1.0f / r.dir.z;
	// lb is the corner of AABB with minimal coordinates - left bottom, rt is maximal corner
	// r.org is origin of ray
	float t1 = (lb.x - r.org.x)*dirfrac.x;
	float t2 = (rt.x - r.org.x)*dirfrac.x;
	float t3 = (lb.y - r.org.y)*dirfrac.y;
	float t4 = (rt.y - r.org.y)*dirfrac.y;
	float t5 = (lb.z - r.org.z)*dirfrac.z;
	float t6 = (rt.z - r.org.z)*dirfrac.z;

	float tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
	float tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));

	// if tmax < 0, ray (line) is intersecting AABB, but the whole AABB is behind us
	if (tmax < 0)
	{
		t = tmax;
		return false;
	}

	// if tmin > tmax, ray doesn't intersect AABB
	if (tmin > tmax)
	{
		t = tmax;
		return false;
	}

	t = tmin;
	return true;
}
*/

vec4 depth_to_view_coord(vec2 tc, float depth) {
    vec4 clip_coord = vec4(vec3(tc, depth) * 2.0 - 1.0, 1.0);
    vec4 view_coord = u_inv_proj_mat * clip_coord;
    return view_coord / view_coord.w;
}

vec4 fetch_voxel(vec3 tc) {
	float a = texture(u_tex_volume, tc).x * u_scale;
	return vec4(u_color, a);
}

const float REF = 150.0;
const float step = 0.01;

void main() {
    vec3 dir = normalize(view_pos);
	vec3 pos = view_pos + dir * step;

	float depth = texelFetch(u_tex_depth, ivec2(gl_FragCoord.xy), 0).x;
	float view_z = depth_to_view_coord(gl_FragCoord.xy * u_inv_res, depth).z;

//	vec3 texCoord = (u_view_to_vol_mat * vec4(pos, 1)).xyz;
//	out_frag = vec4(mix(vec3(1, 0, 0), vec3(0, 1, 0), lessThan(texCoord, vec3(1.0))), 1.0);
//	return;

	vec4 result = vec4(0);
	for (float t = 0; t < 100.0; t += step) {
		if (pos.z < view_z) break;
		vec3 tc	= (u_view_to_vol_mat * vec4(pos, 1)).xyz;
		if (any(lessThan(tc, vec3(0))) || any(greaterThan(tc, vec3(1)))) break;
		vec4 rgba = fetch_voxel(tc);
		
		rgba.a     = 1.0 - pow(1.0 - rgba.a, step * REF);

		rgba.rgb *= rgba.a;
		result += (1.0 - result.a) * rgba;

		//result.rgb = result.rgb + (1.0 - result.a) * rgba.a * rgba.rgb;
		//result.a   = result.a + (1.0 - result.a) * rgba.a;
		if (result.a > 0.99) break;
		pos += dir * t;
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

GLuint create_volume_texture() {
    GLuint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_3D, tex);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
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

mat4 compute_model_to_world_matrix(const vec3& min_world_aabb, const vec3& max_world_aabb) {
    vec3 ext = max_world_aabb - min_world_aabb;
    return mat4(vec4(ext.x, 0, 0, 0), vec4(0, ext.y, 0, 0), vec4(0, 0, ext.z, 0), vec4(min_world_aabb, 1));
}

mat4 compute_data_to_model_matrix(const ivec3& dim) { return mat4(1.f); }

void dump_volume(const Volume& volume, const char* file) {
    std::ofstream of(file, std::ios::binary);
    if (!of.is_open()) return;

    of.write(reinterpret_cast<char*>(volume.voxel_data.data), volume.voxel_data.size_in_bytes());
}

void render_volume_texture(GLuint volume_texture, GLuint depth_texture, const mat4& basis, const mat4& view_matrix, const mat4& proj_matrix,
                           vec3 color, float opacity_scale) {
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    const mat4 vol_to_view_matrix = view_matrix * basis;
    const mat4 view_to_vol_matrix = math::inverse(vol_to_view_matrix);
    const mat4 inv_proj_matrix = math::inverse(proj_matrix);
    const vec2 inv_res = vec2(1.f / (float)(viewport[2]), 1.f / (float)(viewport[3]));

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
    glUniformMatrix4fv(uniform_loc_view_matrix, 1, GL_FALSE, &vol_to_view_matrix[0][0]);
    glUniformMatrix4fv(uniform_loc_proj_matrix, 1, GL_FALSE, &proj_matrix[0][0]);
    glUniformMatrix4fv(uniform_loc_view_to_vol_matrix, 1, GL_FALSE, &view_to_vol_matrix[0][0]);
    glUniformMatrix4fv(uniform_loc_inv_proj_matrix, 1, GL_FALSE, &inv_proj_matrix[0][0]);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_DEPTH_CLAMP);
    glDepthFunc(GL_LESS);
    glDepthMask(GL_TRUE);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 42);
    glBindVertexArray(0);

    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_CLAMP);
}

}  // namespace volume
