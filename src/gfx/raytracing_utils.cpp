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
        GLint voxel_tex = -1;
        GLint voxel_grid_min = -1;
        GLint voxel_grid_size = -1;
        GLint voxel_dimensions = -1;
        GLint voxel_extent = -1;
        GLint use_indirect_diffuse = -1;
        GLint use_indirect_specular = -1;
        GLint use_ambient_occlusion = -1;
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

// Voxel stuff
uniform sampler3D u_voxel_texture;
uniform vec3 u_voxel_grid_world_size;
uniform vec3 u_voxel_grid_world_min;
uniform ivec3 u_voxel_dimensions;
uniform float u_voxel_extent;

// Toggle "booleans"
uniform float u_use_indirect_diffuse;
uniform float u_use_indirect_specular;
uniform float u_use_ambient_occluision;

// View parameters
uniform mat4 u_inv_view_mat;
uniform mat4 u_inv_view_proj_mat;
uniform vec3 u_world_space_camera;

const float MAX_DIST = 100.0;
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
    // Why not z here???
    vec3 tc = (world_position - u_voxel_grid_world_min) / u_voxel_grid_world_size;
	vec3 offset = vec3(0.5) / vec3(u_voxel_dimensions);
	tc = tc + offset;
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

        // front-to-back compositing
        float a = (1.0 - rgba.a);
        rgba += a * voxel_color;
        //color += a * voxel_color.rgb;
        //alpha += a * voxel_color.a;
        //occlusion += a * voxelColor.a;
        occlusion += (a * voxel_color.a) / (1.0 + 0.03 * diameter);
        dist += diameter * 0.5; // smoother
        //dist += diameter; // faster but misses more voxels
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

void main() {
    float depth = texture(u_depth_texture, uv).r;
	if (depth == 1.0) discard;

    vec2 encoded_normal = texture(u_normal_texture, uv).rg;
    vec3 world_position = depth_to_world_coord(uv, depth).xyz;
    vec3 world_normal = mat3(u_inv_view_mat) * decode_normal(encoded_normal);
    vec3 world_eye = u_world_space_camera;

    mat3 tangent_to_world = compute_ON_basis(world_normal);
  
    // Calculate diffuse light
    vec3 diffuse_contribution;
    {
        // Indirect diffuse light
        float occlusion = 0.0;
        vec3 indirect_diffuse = indirect_light(world_position, world_normal, tangent_to_world, occlusion).rgb;

        // Sum direct and indirect diffuse light and tweak a little bit
        occlusion = min(1.0, 1.5 * occlusion); // Make occlusion brighter
        diffuse_contribution = u_use_indirect_diffuse * 2.0 * occlusion * 4.0 * indirect_diffuse;
    }
    
    // Calculate specular light
    vec3 specular_contribution;
    {
        // 0.2 = 22.6 degrees, 0.1 = 11.4 degrees, 0.07 = 8 degrees angle
        const float cone_angle = 0.07;
        vec3 reflect_dir = normalize(world_eye - 2.0 * dot(world_eye, world_normal) * world_normal);

        // Maybe fix so that the cone doesnt trace below the plane defined by the surface normal.
        // For example so that the floor doesnt reflect itself when looking at it with a small angle
        float occlusion;
        vec4 traced_specular = cone_trace(world_position, world_normal, reflect_dir, cone_angle, occlusion); 
        specular_contribution = u_use_indirect_specular * 2.0 * traced_specular.rgb;
    }

    frag_color = vec4(diffuse_contribution * 0.4, 1.0);
	//frag_color = vec4(sample_voxels(world_position, 0).rgb, 1.0);
	//frag_color = vec4(tc, 1.0);
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
    gl.uniform_location.voxel_tex = glGetUniformLocation(gl.program, "u_voxel_texture");
    gl.uniform_location.voxel_grid_min = glGetUniformLocation(gl.program, "u_voxel_grid_world_min");
    gl.uniform_location.voxel_grid_size = glGetUniformLocation(gl.program, "u_voxel_grid_world_size");
    gl.uniform_location.voxel_dimensions = glGetUniformLocation(gl.program, "u_voxel_dimensions");
    gl.uniform_location.voxel_extent = glGetUniformLocation(gl.program, "u_voxel_extent");
    gl.uniform_location.use_indirect_diffuse = glGetUniformLocation(gl.program, "u_use_indirect_diffuse");
    gl.uniform_location.use_indirect_specular = glGetUniformLocation(gl.program, "u_use_indirect_specular");
    gl.uniform_location.use_ambient_occlusion = glGetUniformLocation(gl.program, "u_use_ambient_occlusion");
    gl.uniform_location.inv_view_mat = glGetUniformLocation(gl.program, "u_inv_view_mat");
    gl.uniform_location.inv_view_proj_mat = glGetUniformLocation(gl.program, "u_inv_view_proj_mat");
    gl.uniform_location.world_space_camera = glGetUniformLocation(gl.program, "u_world_space_camera");

    glGenTextures(1, &gl.voxel_texture);
    glBindTexture(GL_TEXTURE_3D, gl.voxel_texture);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
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
    DynamicArray<float> counter(cone_trace::volume.voxel_data.count, 0);

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
                    immediate::draw_aabb(min_box, max_box, cone_trace::volume.voxel_data[i], false);
                }
            }
        }
    }

    immediate::flush();
}

void cone_trace_scene(GLuint depth_tex, GLuint normal_tex, const mat4& view_mat, const mat4& proj_mat) {

    float use_indirect_diffuse = 1.f;
    float use_indirect_specular = 1.f;
    float use_ambient_occlusion = 1.f;

    mat4 inv_view_proj_mat = math::inverse(proj_mat * view_mat);
    mat4 inv_view_mat = math::inverse(view_mat);
    vec3 world_space_camera = inv_view_mat * vec4(0, 0, 0, 1);
    vec3 voxel_grid_min = cone_trace::volume.min_box;
    vec3 voxel_grid_ext = cone_trace::volume.max_box - cone_trace::volume.min_box;
    float voxel_ext = math::max(math::max(cone_trace::volume.voxel_ext.x, cone_trace::volume.voxel_ext.y), cone_trace::volume.voxel_ext.z);

    glUseProgram(cone_trace::gl.program);

    glUniform1i(cone_trace::gl.uniform_location.depth_tex, 0);
    glUniform1i(cone_trace::gl.uniform_location.normal_tex, 1);
    glUniform1i(cone_trace::gl.uniform_location.voxel_tex, 2);
    glUniform3fv(cone_trace::gl.uniform_location.voxel_grid_min, 1, &voxel_grid_min[0]);
    glUniform3fv(cone_trace::gl.uniform_location.voxel_grid_size, 1, &voxel_grid_ext[0]);
    glUniform3iv(cone_trace::gl.uniform_location.voxel_dimensions, 1, &cone_trace::volume.dim[0]);
    glUniform1f(cone_trace::gl.uniform_location.voxel_extent, voxel_ext);
    glUniform1f(cone_trace::gl.uniform_location.use_indirect_diffuse, use_indirect_diffuse);
    glUniform1f(cone_trace::gl.uniform_location.use_indirect_specular, use_indirect_specular);
    glUniform1f(cone_trace::gl.uniform_location.use_ambient_occlusion, use_ambient_occlusion);
    glUniformMatrix4fv(cone_trace::gl.uniform_location.inv_view_mat, 1, GL_FALSE, &inv_view_mat[0][0]);
    glUniformMatrix4fv(cone_trace::gl.uniform_location.inv_view_proj_mat, 1, GL_FALSE, &inv_view_proj_mat[0][0]);
    glUniform3fv(cone_trace::gl.uniform_location.world_space_camera, 1, &world_space_camera[0]);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, depth_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, normal_tex);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_3D, cone_trace::gl.voxel_texture);

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);
    glColorMask(1, 1, 1, 0);
    glDepthMask(GL_FALSE);

    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);

    glDisable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);
    glColorMask(1, 1, 1, 1);
    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);
}

}  // namespace render
