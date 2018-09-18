//#define NOMINMAX

#include "raytracing_utils.h"
#include <core/gl.h>
#include <core/log.h>
#include <gfx/gl_utils.h>
#include <gfx/immediate_draw_utils.h>

struct VoxelData {
    ivec3 vol_dim = {256, 256, 256};
    DynamicArray<uint32> voxel_data;
};

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

void voxelize_scene(Array<const vec3> atom_pos, Array<const uint32> atom_color, vec3 min_box, vec3 max_box) {}

void draw_voxelized_scene(const mat4& view_mat, const mat4& proj_mat) {}

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
