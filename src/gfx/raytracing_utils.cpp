//#define NOMINMAX

#include "raytracing_utils.h"
#include <core/gl.h>
#include <core/log.h>
#include <gfx/gl_utils.h>
#include <gfx/immediate_draw_utils.h>

struct VoxelData {
    DynamicArray<uint32> data{};
    ivec3 dim = {256, 256, 256};
    vec3 min_box = {0, 0, 0};
    vec3 max_box = {0, 0, 0};
    vec3 voxel_ext = {0, 0, 0};
};

static VoxelData voxel_data;

namespace render {

static GLuint vbo = 0;
static GLsizeiptr vbo_size = MEGABYTES(4);

inline void set_vbo_data(const void* data, GLsizeiptr size_in_bytes) {
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    if (vbo_size > size_in_bytes) {
        glBufferSubData(GL_ARRAY_BUFFER, 0, size_in_bytes, data);
    } else {
        glBufferData(GL_ARRAY_BUFFER, size_in_bytes, data, GL_STREAM_DRAW);
        vbo_size = size_in_bytes;
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

namespace shit {
static GLuint v_shader = 0;
static GLuint f_shader = 0;
static GLuint program = 0;

static const char* v_shader_src = R"(
#version 150 core
void main() {

}
)";

static const char* f_shader_src = R"(
#version 150 core
void main() {

}
)";

static void initialize() {
    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    v_shader = glCreateShader(GL_VERTEX_SHADER);
    f_shader = glCreateShader(GL_FRAGMENT_SHADER);

    glShaderSource(v_shader, 1, &v_shader_src, 0);
    glShaderSource(f_shader, 1, &f_shader_src, 0);

    glCompileShader(v_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, v_shader)) {
        LOG_ERROR("Compiling vdw vertex shader:\n%s\n", buffer);
    }

    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        LOG_ERROR("Compiling vdw fragment shader:\n%s\n", buffer);
    }

    program = glCreateProgram();
    glAttachShader(program, v_shader);
    glAttachShader(program, f_shader);
    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, program)) {
        LOG_ERROR("Linking vdw program:\n%s\n", buffer);
    }

    glDetachShader(program, v_shader);
    glDetachShader(program, f_shader);

    glDeleteShader(v_shader);
    glDeleteShader(f_shader);
}

static void shutdown() {
    if (program) glDeleteProgram(program);
}

}  // namespace shit

void initialize() {

    if (!vbo) {
        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, vbo_size, nullptr, GL_STREAM_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    shit::initialize();
}

void shutdown() {
    if (vbo) glDeleteBuffers(1, &vbo);

    shit::shutdown();
}

inline ivec3 compute_voxel_coord(const VoxelData& data, const vec3& coord) {
    return math::clamp(ivec3((coord - data.min_box) / data.voxel_ext), ivec3(0), data.dim - 1);
}

inline int compute_voxel_idx(const ivec3& res, const ivec3& coord) { return coord.z * res.x * res.y + coord.y * res.x + coord.x; }

inline int compute_voxel_idx(const VoxelData& data, const vec3& coord) { return compute_voxel_idx(data.dim, compute_voxel_coord(data, coord)); }

inline uint32 accumulate_voxel_color(uint32 current_color, uint32 new_color, uint32 counter) {
    // @TODO: Implement proper color blending
    return new_color;
}

void voxelize_scene(Array<const vec3> atom_pos, Array<const uint32> atom_color, ivec3 resolution, vec3 min_box, vec3 max_box) {
    ASSERT(atom_pos.count == atom_color.count);

    if (min_box == vec3(0, 0, 0) && max_box == vec3(0, 0, 0)) {
        min_box = vec3(FLT_MAX);
        max_box = vec3(-FLT_MAX);
        for (const auto& p : atom_pos) {
            min_box = math::min(min_box, p);
            max_box = math::max(max_box, p);
        }
    }

    voxel_data.dim = resolution;
    voxel_data.data.resize(resolution.x * resolution.y * resolution.z);
    voxel_data.data.set_mem_to_zero();
    voxel_data.min_box = min_box;
    voxel_data.max_box = max_box;
    voxel_data.voxel_ext = (max_box - min_box) / vec3(resolution);

    // For running mean
    DynamicArray<uint8> counter(voxel_data.data.count, 0);

    for (int32 i = 0; i < atom_pos.count; i++) {
        const auto& pos = atom_pos[i];
        const auto& col = atom_color[i];
        ivec3 coord = compute_voxel_coord(voxel_data, pos);
        int voxel_idx = compute_voxel_idx(voxel_data.dim, coord);
        voxel_data.data[voxel_idx] = accumulate_voxel_color(voxel_data.data[voxel_idx], col, counter[voxel_idx]);
    }
}

void draw_voxelized_scene(const mat4& view_mat, const mat4& proj_mat) {
    immediate::set_view_matrix(view_mat);
    immediate::set_proj_matrix(proj_mat);

    for (int32 z = 0; z < voxel_data.dim.z; z++) {
        for (int32 y = 0; y < voxel_data.dim.y; y++) {
            for (int32 x = 0; x < voxel_data.dim.x; x++) {
                int32 i = compute_voxel_idx(voxel_data.dim, ivec3(x, y, z));
                if (voxel_data.data[i] > 0) {
                    vec3 min_box = voxel_data.min_box + vec3(x, y, z) * voxel_data.voxel_ext;
                    vec3 max_box = min_box + voxel_data.voxel_ext;
                    immediate::draw_aabb(min_box, max_box, voxel_data.data[i], true);
                }
            }
        }
    }

    immediate::flush();
}

void cone_trace_scene(const mat4& view_mat, const mat4& proj_mat) {}

void draw_spatial_hash_cells(const spatialhash::Frame& frame, const mat4& view_mat, const mat4& proj_mat) {
    immediate::set_view_matrix(view_mat);
    immediate::set_proj_matrix(proj_mat);

    for (int32 z = 0; z < frame.cell_count.z; z++) {
        for (int32 y = 0; y < frame.cell_count.y; y++) {
            for (int32 x = 0; x < frame.cell_count.x; x++) {
                int32 i = frame.cell_count.x * frame.cell_count.y * z + frame.cell_count.x * y + x;
                if (frame.cells[i].count > 0) {
                    vec3 min_box = frame.min_box + vec3(x, y, z) * frame.cell_ext;
                    vec3 max_box = min_box + frame.cell_ext;
                    immediate::draw_aabb(min_box, max_box, immediate::COLOR_GREEN);
                }
            }
        }
    }

    immediate::flush();
}

}  // namespace render
