//#define NOMINMAX

#include "raytracing_utils.h"
#include <core/log.h>
#include <gfx/gl_utils.h>
#include <gfx/immediate_draw_utils.h>

struct Volume {
    DynamicArray<uint32> voxel_data{};
    ivec3 dim = {256, 256, 256};
    vec3 min_box = {0, 0, 0};
    vec3 max_box = {0, 0, 0};
    vec3 voxel_ext = {0, 0, 0};
};

namespace render {

static GLuint vao = 0;
static GLuint vbo = 0;

static const char* v_shader_fs_quad_src = R"(
#version 150 core

out vec2 uv;

void main() {
	uint idx = uint(gl_VertexID) % 3U;
	gl_Position = vec4(
		(float( idx     &1U)) * 4.0 - 1.0,
		(float((idx>>1U)&1U)) * 4.0 - 1.0,
		0, 1.0);
	uv = gl_Position.xy * 0.5 + 0.5;
}
)";

namespace cone_trace {

static Volume volume;

static struct {
    GLuint program = 0;
    GLuint voxel_texture;

    struct {
        GLint depth_tex = -1;
        GLint normal_tex = -1;
        GLint color_alpha_tex = -1;
        GLint f0_smoothness_tex = -1;
        GLint voxel_tex = -1;
        GLint voxel_grid_min = -1;
        GLint voxel_grid_size = -1;
        GLint voxel_dimensions = -1;
        GLint voxel_extent = -1;
        GLint indirect_diffuse_scale = -1;
        GLint indirect_specular_scale = -1;
        GLint ambient_occlusion_scale = -1;
        GLint cone_angle = -1;
        GLint inv_view_mat = -1;
        GLint inv_view_proj_mat = -1;
        GLint world_space_camera = -1;
    } uniform_location;

} gl;

// Modified version of
// https://github.com/Cigg/Voxel-Cone-Tracing

static const char* f_shader_src = R"(
#version 330 core

in vec2 uv;
out vec4 frag_color;

// Textures
uniform sampler2D u_depth_texture;
uniform sampler2D u_normal_texture;
uniform sampler2D u_color_alpha_texture;
uniform sampler2D u_f0_smoothness_texture;

// Voxel stuff
uniform sampler3D u_voxel_texture;
uniform vec3 u_voxel_grid_world_size;
uniform vec3 u_voxel_grid_world_min;
uniform ivec3 u_voxel_dimensions;
uniform float u_voxel_extent;

// Scaling factors
uniform float u_indirect_diffuse_scale = 1.f;
uniform float u_indirect_specular_scale = 1.f;
uniform float u_ambient_occlusion_scale = 1.f;
uniform float u_cone_angle = 0.08;

// View parameters
uniform mat4 u_inv_view_mat;
uniform mat4 u_inv_view_proj_mat;
uniform vec3 u_world_space_camera;

const float MAX_DIST = 1000.0;
const float ALPHA_THRESH = 0.95;

// 6 60 degree cone
const int NUM_CONES = 6;
const vec3 cone_directions[6] = vec3[]
(                            vec3(0, 1, 0),
                            vec3(0, 0.5, 0.866025),
                            vec3(0.823639, 0.5, 0.267617),
                            vec3(0.509037, 0.5, -0.700629),
                            vec3(-0.509037, 0.5, -0.700629),
                            vec3(-0.823639, 0.5, 0.267617)
                            );
const float cone_weights[6] = float[](0.25, 0.15, 0.15, 0.15, 0.15, 0.15);

// // 5 90 degree cones
// const int NUM_CONES = 5;
// vec3 coneDirections[5] = vec3[]
// (                            vec3(0, 1, 0),
//                             vec3(0, 0.707, 0.707),
//                             vec3(0, 0.707, -0.707),
//                             vec3(0.707, 0.707, 0),
//                             vec3(-0.707, 0.707, 0)
//                             );
// float coneWeights[5] = float[](0.28, 0.18, 0.18, 0.18, 0.18);

vec4 sample_voxels(vec3 world_position, float lod) {
    vec3 tc = (world_position - u_voxel_grid_world_min) / (u_voxel_grid_world_size);
	//if (any(lessThan(tc, vec3(0))) || any(greaterThan(tc, vec3(1)))) return vec4(-1);
	//tc = tc + vec3(0.5) / vec3(u_voxel_dimensions);
    return textureLod(u_voxel_texture, tc, lod);
}

// Third argument to say how long between steps?
vec4 cone_trace(vec3 world_position, vec3 world_normal, vec3 direction, float tan_half_angle, out float occlusion) {
    
    // lod level 0 mipmap is full size, level 1 is half that size and so on
    float lod = 0.0;
    vec4 rgba = vec4(0);
    occlusion = 0.0;

    float dist = u_voxel_extent; // Start one voxel away to avoid self occlusion
    vec3 start_pos = world_position + world_normal * u_voxel_extent; // Plus move away slightly in the normal direction to avoid
                                                                     // self occlusion in flat surfaces

    while(dist < MAX_DIST && rgba.a < ALPHA_THRESH) {
        // smallest sample diameter possible is the voxel size
        float diameter = max(u_voxel_extent, 2.0 * tan_half_angle * dist);
        float lod_level = log2(diameter / u_voxel_extent);
        vec4 voxel_color = sample_voxels(start_pos + dist * direction, lod_level);
		//if (voxel_color.a < 0) break;

        // front-to-back compositing
        float a = (1.0 - rgba.a);
        rgba += a * voxel_color;
        //color += a * voxel_color.rgb;
        //alpha += a * voxel_color.a;
        //occlusion += a * voxelColor.a;
        occlusion += (a * voxel_color.a) / (1.0 + 0.03 * diameter);
        //dist += diameter * 0.5; // smoother
        dist += diameter; // faster but misses more voxels
    }

    return rgba;
}

vec4 indirect_light(in vec3 world_position, in vec3 world_normal, in mat3 tangent_to_world, out float occlusion_out) {
    vec4 color = vec4(0);
    occlusion_out = 0.0;

    for(int i = 0; i < NUM_CONES; i++) {
        float occlusion = 0.0;
        // 60 degree cones -> tan(30) = 0.577
        // 90 degree cones -> tan(45) = 1.0
        color += cone_weights[i] * cone_trace(world_position, world_normal, tangent_to_world * cone_directions[i], 0.577, occlusion);
        occlusion_out += cone_weights[i] * occlusion;
    }

    occlusion_out = 1.0 - occlusion_out;

    return color;
}

mat3 compute_ON_basis(in vec3 v1) {
    vec3 v0;
    vec3 v2;
    float d0 = dot(vec3(1,0,0), v1);
    float d1 = dot(vec3(0,0,1), v1);
    if (d0 < d1) {
        v0 = normalize(vec3(1,0,0) - v1 * d0);
    } else {
        v0 = normalize(vec3(0,0,1) - v1 * d1);
    }
    v2 = cross(v0, v1);
    return mat3(v0, v1, v2);
}

vec4 depth_to_world_coord(vec2 tex_coord, float depth) {
    vec4 clip_coord = vec4(vec3(tex_coord, depth) * 2.0 - 1.0, 1.0);
    vec4 world_coord = u_inv_view_proj_mat * clip_coord;
    return world_coord / world_coord.w;
}

// https://aras-p.info/texts/CompactNormalStorage.html
vec3 decode_normal(vec2 enc) {
    vec2 fenc = enc*4-2;
    float f = dot(fenc,fenc);
    float g = sqrt(1-f/4.0);
    vec3 n;
    n.xy = fenc*g;
    n.z = 1-f/2.0;
    return n;
}

vec3 fresnel(vec3 f0, float H_dot_V) {
    //const float n1 = 1.0;
    //const float n2 = 1.5;
    //const float R0 = pow((n1-n2)/(n1+n2), 2);
    return f0 + (1.0 - f0)*pow(1.0 - H_dot_V, 5);
}

vec3 shade(vec3 albedo, float alpha, vec3 f0, float smoothness, vec3 P, vec3 V, vec3 N) {
    const vec3 env_radiance = vec3(0.5);
    const vec3 dir_radiance = vec3(0.5);
    const vec3 L = normalize(vec3(1));
    const float spec_exp = 10.0;

    mat3 tangent_to_world = compute_ON_basis(N);

    float N_dot_V = max(0.0, dot(N, V));
    vec3 R = -V + 2.0 * dot(N, V) * N;
    vec3 H = normalize(L + V);
    float H_dot_V = max(0.0, dot(H, V));
    float N_dot_H = max(0.0, dot(N, H));
    float N_dot_L = max(0.0, dot(N, L));
    vec3 fresnel_direct = fresnel(f0, H_dot_V);
    vec3 fresnel_indirect = fresnel(f0, N_dot_V);
	float tan_half_angle = tan(mix(45.0, 0.0, smoothness));

    float diffuse_occlusion;
    vec3 direct_diffuse = env_radiance + N_dot_L * dir_radiance;
    vec3 indirect_diffuse = indirect_light(P, N, tangent_to_world, diffuse_occlusion).rgb;

    float specular_occlusion;
    vec3 direct_specular = dir_radiance * pow(N_dot_H, spec_exp);
    vec3 indirect_specular = cone_trace(P, N, R, tan_half_angle, specular_occlusion).rgb; 

    vec3 result = vec3(0);
    result += albedo * (direct_diffuse + u_indirect_diffuse_scale * indirect_diffuse) * pow(diffuse_occlusion, u_ambient_occlusion_scale * 0.5);
    result += direct_specular * fresnel_direct + u_indirect_specular_scale * indirect_specular * fresnel_indirect;

    return vec3(smoothness);
}

void main() {
    float depth = texture(u_depth_texture, uv).r;
	if (depth == 1.0) discard;

    vec2 encoded_normal = texture(u_normal_texture, uv).rg;
	vec4 albedo_alpha = texture(u_color_alpha_texture, uv);
	vec4 f0_smoothness = texture(u_f0_smoothness_texture, uv);
	vec3 albedo = albedo_alpha.rgb;
	float alpha = albedo_alpha.a;
	vec3 f0 = f0_smoothness.rgb;
	float smoothness = f0_smoothness.a;

    vec3 world_position = depth_to_world_coord(uv, depth).xyz;
    vec3 world_normal = mat3(u_inv_view_mat) * decode_normal(encoded_normal);
    vec3 world_eye = u_world_space_camera;

	vec3 P = world_position;
	vec3 V = normalize(world_eye - world_position);
	vec3 N = world_normal;

    frag_color = vec4(shade(albedo, alpha, f0, smoothness, P, V, N), 1);

    /*

	vec3 L = normalize(vec3(1,1,1));
	vec3 N = world_normal;
	vec3 E = normalize(world_eye - world_position);
	float cos_term = max(0, dot(L, N));

    mat3 tangent_to_world = compute_ON_basis(world_normal);
  
    // Calculate diffuse light
    vec3 diffuse_contribution = vec3(0);
    {
		vec3 direct_diffuse = cos_term * vec3(1,1,1);

        // Indirect diffuse light
        float occlusion = 0.0;
        vec3 indirect_diffuse = indirect_light(world_position, world_normal, tangent_to_world, occlusion).rgb;

        // Sum direct and indirect diffuse light and tweak a little bit
        // occlusion = min(1.0, 1.5 * occlusion); // Make occlusion brighter
		// TODO: apply occlusion
        diffuse_contribution = (direct_diffuse + u_indirect_diffuse_scale * indirect_diffuse) * albedo;

        diffuse_contribution = diffuse_contribution * u_ambient_occlusion_scale * occlusion;
    }
    
    // Calculate specular light
    vec3 specular_contribution = vec3(0);
    {
        // 0.2 = 22.6 degrees, 0.1 = 11.4 degrees, 0.07 = 8 degrees angle
        const float cone_angle = 0.07;

        vec3 E_refl = -E + 2.0 * dot(N, E) * N;

        // Maybe fix so that the cone doesnt trace below the plane defined by the surface normal.
        // For example so that the floor doesnt reflect itself when looking at it with a small angle
        float occlusion;
        vec4 traced_specular = cone_trace(world_position, world_normal, E_refl, cone_angle, occlusion); 
        specular_contribution = u_indirect_specular_scale * traced_specular.rgb;
    }
    */

    //frag_color = vec4(diffuse_contribution + specular_contribution, 1.0);
}
)";

static void initialize() {
    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    GLuint v_shader = glCreateShader(GL_VERTEX_SHADER);
    GLuint f_shader = glCreateShader(GL_FRAGMENT_SHADER);

    glShaderSource(v_shader, 1, &v_shader_fs_quad_src, 0);
    glShaderSource(f_shader, 1, &f_shader_src, 0);

    glCompileShader(v_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, v_shader)) {
        LOG_ERROR("Compiling cone_tracing vertex shader:\n%s\n", buffer);
    }

    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        LOG_ERROR("Compiling cone_tracing fragment shader:\n%s\n", buffer);
    }

    gl.program = glCreateProgram();
    glAttachShader(gl.program, v_shader);
    glAttachShader(gl.program, f_shader);
    glLinkProgram(gl.program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, gl.program)) {
        LOG_ERROR("Linking cone_tracing program:\n%s\n", buffer);
    }

    glDetachShader(gl.program, v_shader);
    glDetachShader(gl.program, f_shader);

    glDeleteShader(v_shader);
    glDeleteShader(f_shader);

    gl.uniform_location.depth_tex = glGetUniformLocation(gl.program, "u_depth_texture");
    gl.uniform_location.normal_tex = glGetUniformLocation(gl.program, "u_normal_texture");
    gl.uniform_location.color_alpha_tex = glGetUniformLocation(gl.program, "u_color_alpha_texture");
    gl.uniform_location.f0_smoothness_tex = glGetUniformLocation(gl.program, "u_f0_smoothness_texture");
    gl.uniform_location.voxel_tex = glGetUniformLocation(gl.program, "u_voxel_texture");
    gl.uniform_location.voxel_grid_min = glGetUniformLocation(gl.program, "u_voxel_grid_world_min");
    gl.uniform_location.voxel_grid_size = glGetUniformLocation(gl.program, "u_voxel_grid_world_size");
    gl.uniform_location.voxel_dimensions = glGetUniformLocation(gl.program, "u_voxel_dimensions");
    gl.uniform_location.voxel_extent = glGetUniformLocation(gl.program, "u_voxel_extent");
    gl.uniform_location.indirect_diffuse_scale = glGetUniformLocation(gl.program, "u_indirect_diffuse_scale");
    gl.uniform_location.indirect_specular_scale = glGetUniformLocation(gl.program, "u_indirect_specular_scale");
    gl.uniform_location.ambient_occlusion_scale = glGetUniformLocation(gl.program, "u_ambient_occlusion_scale");
    gl.uniform_location.cone_angle = glGetUniformLocation(gl.program, "u_cone_angle");
    gl.uniform_location.inv_view_mat = glGetUniformLocation(gl.program, "u_inv_view_mat");
    gl.uniform_location.inv_view_proj_mat = glGetUniformLocation(gl.program, "u_inv_view_proj_mat");
    gl.uniform_location.world_space_camera = glGetUniformLocation(gl.program, "u_world_space_camera");

    glGenTextures(1, &gl.voxel_texture);
    glBindTexture(GL_TEXTURE_3D, gl.voxel_texture);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);
    glBindTexture(GL_TEXTURE_3D, 0);
}

static void shutdown() {
    if (gl.program) {
        glDeleteProgram(gl.program);
        gl.program = 0;
    }
    if (gl.voxel_texture) {
        glDeleteTextures(1, &gl.voxel_texture);
        gl.voxel_texture = 0;
    }
}

}  // namespace cone_trace

void initialize() {
    if (!vao) glGenVertexArrays(1, &vao);
    if (!vbo) glGenBuffers(1, &vbo);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, 12, nullptr, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid*)0);
    glBindVertexArray(0);

    cone_trace::initialize();
}

void shutdown() { cone_trace::shutdown(); }

inline ivec3 compute_voxel_coord(const Volume& data, const vec3& coord) {
    return math::clamp(ivec3((coord - data.min_box) / data.voxel_ext), ivec3(0), data.dim - 1);
}

inline int compute_voxel_idx(const ivec3& res, const ivec3& coord) { return coord.z * res.x * res.y + coord.y * res.x + coord.x; }

inline int compute_voxel_idx(const Volume& data, const vec3& coord) { return compute_voxel_idx(data.dim, compute_voxel_coord(data, coord)); }

inline uint32 accumulate_voxel_color(uint32 current_color, uint32 new_color, float counter) {
    // @TODO: Implement proper color blending
    vec4 c = math::convert_color(current_color);
    vec4 n = math::convert_color(new_color);

    c = (counter * c + n) / (counter + 1.f);
    return math::convert_color(c);
}

void voxelize_scene(Array<const vec3> atom_pos, Array<const float> atom_radii, Array<const uint32> atom_colors, ivec3 resolution, vec3 min_box,
                    vec3 max_box) {
    ASSERT(atom_pos.count == atom_radii.count);
    ASSERT(atom_pos.count == atom_colors.count);
    const int32 N = atom_pos.count;

    if (min_box == vec3(0, 0, 0) && max_box == vec3(0, 0, 0)) {
        min_box = vec3(FLT_MAX);
        max_box = vec3(-FLT_MAX);
        for (int32 i = 0; i < N; i++) {
            min_box = math::min(min_box, atom_pos[i] - atom_radii[i]);
            max_box = math::max(max_box, atom_pos[i] + atom_radii[i]);
        }
    }

    cone_trace::volume.dim = resolution;
    cone_trace::volume.voxel_data.resize(resolution.x * resolution.y * resolution.z);
    cone_trace::volume.voxel_data.set_mem_to_zero();
    cone_trace::volume.min_box = min_box;
    cone_trace::volume.max_box = max_box;
    cone_trace::volume.voxel_ext = (max_box - min_box) / vec3(resolution);

    // For running mean
    DynamicArray<float> counter(cone_trace::volume.voxel_data.count, 1);

    for (int32 i = 0; i < N; i++) {
        const auto& pos = atom_pos[i];
        const auto& rad = atom_radii[i];
        const auto& col = atom_colors[i];
        const float r2 = rad * rad;
        ivec3 min_cc = compute_voxel_coord(cone_trace::volume, pos - rad);
        ivec3 max_cc = compute_voxel_coord(cone_trace::volume, pos + rad);
        ivec3 cc;
        for (cc.z = min_cc.z; cc.z <= max_cc.z; cc.z++) {
            for (cc.y = min_cc.y; cc.y <= max_cc.y; cc.y++) {
                for (cc.x = min_cc.x; cc.x <= max_cc.x; cc.x++) {
                    vec3 min_voxel = cone_trace::volume.min_box + vec3(cc) * cone_trace::volume.voxel_ext;
                    vec3 max_voxel = min_voxel + cone_trace::volume.voxel_ext;
                    vec3 clamped_pos = math::clamp(pos, min_voxel, max_voxel);
                    vec3 d = clamped_pos - pos;

                    if (dot(d, d) < r2) {
                        int voxel_idx = compute_voxel_idx(cone_trace::volume.dim, cc);
                        cone_trace::volume.voxel_data[voxel_idx] =
                            accumulate_voxel_color(cone_trace::volume.voxel_data[voxel_idx], col, counter[voxel_idx]);
                        counter[voxel_idx]++;
                    }
                }
            }
        }
    }

    glBindTexture(GL_TEXTURE_3D, cone_trace::gl.voxel_texture);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA8, cone_trace::volume.dim.x, cone_trace::volume.dim.y, cone_trace::volume.dim.z, 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, cone_trace::volume.voxel_data.data);
    glGenerateMipmap(GL_TEXTURE_3D);
    glBindTexture(GL_TEXTURE_3D, 0);
}  // namespace render

void draw_voxelized_scene(const mat4& view_mat, const mat4& proj_mat) {
    immediate::set_view_matrix(view_mat);
    immediate::set_proj_matrix(proj_mat);

    for (int32 z = 0; z < cone_trace::volume.dim.z; z++) {
        for (int32 y = 0; y < cone_trace::volume.dim.y; y++) {
            for (int32 x = 0; x < cone_trace::volume.dim.x; x++) {
                int32 i = compute_voxel_idx(cone_trace::volume.dim, ivec3(x, y, z));
                if (cone_trace::volume.voxel_data[i] > 0) {
                    vec3 min_box = cone_trace::volume.min_box + vec3(x, y, z) * cone_trace::volume.voxel_ext;
                    vec3 max_box = min_box + cone_trace::volume.voxel_ext;
                    immediate::draw_aabb_lines(min_box, max_box, cone_trace::volume.voxel_data[i]);
                }
            }
        }
    }

    immediate::flush();
}

void cone_trace_scene(GLuint depth_tex, GLuint normal_tex, GLuint color_alpha_tex, GLuint f0_smoothness_tex, const mat4& view_mat,
                      const mat4& proj_mat, float indirect_diffuse_scale, float indirect_specular_scale, float ambient_occlusion_scale) {

    mat4 inv_view_proj_mat = math::inverse(proj_mat * view_mat);
    mat4 inv_view_mat = math::inverse(view_mat);
    vec3 world_space_camera = inv_view_mat * vec4(0, 0, 0, 1);
    // printf("cam: %.2f %.2f %.2f\n", world_space_camera.x, world_space_camera.y, world_space_camera.z);
    vec3 voxel_grid_min = cone_trace::volume.min_box;
    vec3 voxel_grid_ext = cone_trace::volume.max_box - cone_trace::volume.min_box;
    float voxel_ext = math::max(math::max(cone_trace::volume.voxel_ext.x, cone_trace::volume.voxel_ext.y), cone_trace::volume.voxel_ext.z);

    const float cone_angle = 0.07;  // 0.2 = 22.6 degrees, 0.1 = 11.4 degrees, 0.07 = 8 degrees angle

    glUseProgram(cone_trace::gl.program);

    glUniform1i(cone_trace::gl.uniform_location.depth_tex, 0);
    glUniform1i(cone_trace::gl.uniform_location.normal_tex, 1);
    glUniform1i(cone_trace::gl.uniform_location.color_alpha_tex, 2);
    glUniform1i(cone_trace::gl.uniform_location.f0_smoothness_tex, 3);
    glUniform1i(cone_trace::gl.uniform_location.voxel_tex, 4);
    glUniform3fv(cone_trace::gl.uniform_location.voxel_grid_min, 1, &voxel_grid_min[0]);
    glUniform3fv(cone_trace::gl.uniform_location.voxel_grid_size, 1, &voxel_grid_ext[0]);
    glUniform3iv(cone_trace::gl.uniform_location.voxel_dimensions, 1, &cone_trace::volume.dim[0]);
    glUniform1f(cone_trace::gl.uniform_location.voxel_extent, voxel_ext);
    glUniform1f(cone_trace::gl.uniform_location.indirect_diffuse_scale, indirect_diffuse_scale);
    glUniform1f(cone_trace::gl.uniform_location.indirect_specular_scale, indirect_specular_scale);
    glUniform1f(cone_trace::gl.uniform_location.ambient_occlusion_scale, ambient_occlusion_scale);
    glUniformMatrix4fv(cone_trace::gl.uniform_location.inv_view_mat, 1, GL_FALSE, &inv_view_mat[0][0]);
    glUniformMatrix4fv(cone_trace::gl.uniform_location.inv_view_proj_mat, 1, GL_FALSE, &inv_view_proj_mat[0][0]);
    glUniform3fv(cone_trace::gl.uniform_location.world_space_camera, 1, &world_space_camera[0]);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, depth_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, normal_tex);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, color_alpha_tex);

    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, f0_smoothness_tex);

    glActiveTexture(GL_TEXTURE4);
    glBindTexture(GL_TEXTURE_3D, cone_trace::gl.voxel_texture);

    glDisable(GL_DEPTH_TEST);
    // glEnable(GL_BLEND);
    // glBlendFunc(GL_ONE, GL_ONE);
    // glColorMask(1, 1, 1, 0);
    glDepthMask(GL_FALSE);

    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);

    // glDisable(GL_BLEND);
    // glBlendFunc(GL_ONE, GL_ZERO);
    // glColorMask(1, 1, 1, 1);
    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);
}

}  // namespace render
