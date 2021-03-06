﻿#include "image.h"
#include "molecule_draw.h"

#include <core/common.h>
#include <core/log.h>
#include <core/string_utils.h>
#include <core/spatial_hash.h>
#include <gfx/gl_utils.h>
#include <mol/molecule_utils.h>
#include <gfx/immediate_draw_utils.h>

#define PUSH_GPU_SECTION(lbl)                                                                       \
    {                                                                                               \
        if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); \
    }
#define POP_GPU_SECTION()                       \
    {                                           \
        if (glPopDebugGroup) glPopDebugGroup(); \
    }

static bool is_orthographic_proj_matrix(const mat4& M) { return M[2][3] == 0.0f; }

namespace draw {
static GLuint vao = 0;
static GLuint tex[4] = {};
static GLint major_version = 0;
static GLint minor_version = 0;

namespace vdw {
static GLuint program_persp = 0;
static GLuint program_ortho = 0;
static GLuint program_depth_prepass = 0;

static GLuint vbo_icosa = 0;
static GLuint ibo_icosa = 0;

static void initialize() {
    {
        GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/vdw.vert", GL_VERTEX_SHADER);
        GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/vdw.geom", GL_GEOMETRY_SHADER);
        GLuint f_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/vdw.frag", GL_FRAGMENT_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(g_shader);
            glDeleteShader(f_shader);
        };

        if (!program_persp) program_persp = glCreateProgram();
        const GLuint shaders[] = {v_shader, g_shader, f_shader};
        gl::attach_link_detach(program_persp, shaders);
    }
    {
        GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/vdw_ortho.vert", GL_VERTEX_SHADER);
        GLuint f_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/vdw_ortho.frag", GL_FRAGMENT_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(f_shader);
        };

        if (!program_ortho) program_ortho = glCreateProgram();
        const GLuint shaders[] = {v_shader, f_shader};
        gl::attach_link_detach(program_ortho, shaders);
    }
    {
        GLuint shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/vdw_depth_prepass.vert", GL_VERTEX_SHADER);
        defer { glDeleteShader(shader); };

        if (!program_depth_prepass) program_depth_prepass = glCreateProgram();
        gl::attach_link_detach(program_depth_prepass, {&shader, 1});
    }
    {
        constexpr float X = .525731112119133606;
        constexpr float Z = .850650808352039932;
        constexpr GLfloat vertices[12][3] = 
        {    
            {-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},    
            {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},    
            {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0} 
        };

        constexpr GLubyte indices[20][3] = 
        { 
            {0,4,1},  {0,9,4},  {9,5,4},  {4,5,8},  {4,8,1},    
            {8,10,1}, {8,3,10}, {5,3,8},  {5,2,3},  {2,7,3},    
            {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6}, 
            {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5},  {7,2,11} 
        };

        if (!vbo_icosa) {
            glGenBuffers(1, &vbo_icosa);
            glBindBuffer(GL_ARRAY_BUFFER, vbo_icosa);
            glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
        }
        if (!ibo_icosa) {
            glGenBuffers(1, &ibo_icosa);
            glBindBuffer(GL_ARRAY_BUFFER, ibo_icosa);
            glBufferData(GL_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
        }
    }
}

static void shutdown() {
    if (program_persp) glDeleteProgram(program_persp);
    if (program_ortho) glDeleteProgram(program_ortho);
}

}  // namespace vdw

namespace licorice {
static GLuint program_persp = 0;
static GLuint program_ortho = 0;

static void initialize() {
    GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/licorice.vert", GL_VERTEX_SHADER);
    GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/licorice.geom", GL_GEOMETRY_SHADER);
    GLuint f_shader_persp = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/licorice.frag", GL_FRAGMENT_SHADER);
    GLuint f_shader_ortho = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/licorice.frag", GL_FRAGMENT_SHADER, "#define ORTHO");
    defer {
        glDeleteShader(v_shader);
        glDeleteShader(g_shader);
        glDeleteShader(f_shader_persp);
        glDeleteShader(f_shader_ortho);
    };

    if (!program_persp) program_persp = glCreateProgram();
    {
        const GLuint shaders[] = {v_shader, g_shader, f_shader_persp};
        gl::attach_link_detach(program_persp, shaders);
    }

    if (!program_ortho) program_ortho = glCreateProgram();
    {
        const GLuint shaders[] = {v_shader, g_shader, f_shader_ortho};
        gl::attach_link_detach(program_ortho, shaders);
    }
}

static void shutdown() {
    if (program_persp) glDeleteProgram(program_persp);
    if (program_ortho) glDeleteProgram(program_ortho);
}

}  // namespace licorice

namespace ribbon {
static GLuint program = 0;

void intitialize() {
    GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/ribbons.vert", GL_VERTEX_SHADER);
    GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/ribbons.geom", GL_GEOMETRY_SHADER);
    GLuint f_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/ribbons.frag", GL_FRAGMENT_SHADER);
    defer {
        glDeleteShader(v_shader);
        glDeleteShader(g_shader);
        glDeleteShader(f_shader);
    };

    if (!program) program = glCreateProgram();
    const GLuint shaders[] = {v_shader, g_shader, f_shader};
    gl::attach_link_detach(program, shaders);
}

void shutdown() {
    if (program) glDeleteProgram(program);
}

}  // namespace ribbon

namespace cartoon {
static GLuint program = 0;

void intitialize() {
    GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/cartoon.vert", GL_VERTEX_SHADER);
    GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/cartoon.geom", GL_GEOMETRY_SHADER);
    GLuint f_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/cartoon.frag", GL_FRAGMENT_SHADER);
    defer {
        glDeleteShader(v_shader);
        glDeleteShader(g_shader);
        glDeleteShader(f_shader);
    };

    if (!program) program = glCreateProgram();
    const GLuint shaders[] = {v_shader, g_shader, f_shader};
    gl::attach_link_detach(program, shaders);
}

void shutdown() {
    if (program) glDeleteProgram(program);
}

}  // namespace cartoon

namespace backbone_spline {
static GLuint extract_control_points_program = 0;
static GLuint compute_spline_program = 0;
static GLuint draw_spline_program = 0;

void initialize() {
    {
        GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/backbone_control_points.vert", GL_VERTEX_SHADER);
        GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/backbone_control_points.geom", GL_GEOMETRY_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(g_shader);
        };

        const GLuint shaders[] = {v_shader, g_shader};
        const GLchar* feedback_varyings[] = {"out_control_point",     "out_support_vector_xy", "out_support_vector_z_tangent_vector_x",
                                             "out_tangent_vector_yz", "out_classification",    "out_atom_index"};
        if (!extract_control_points_program) extract_control_points_program = glCreateProgram();
        gl::attach_link_detach_with_transform_feedback(extract_control_points_program, shaders, feedback_varyings, GL_INTERLEAVED_ATTRIBS);
    }

    {
        GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/backbone_compute_spline.vert", GL_VERTEX_SHADER);
        GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/backbone_compute_spline.geom", GL_GEOMETRY_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(g_shader);
        };

        const GLuint shaders[] = {v_shader, g_shader};
        const GLchar* feedback_varyings[] = {"out_control_point",     "out_support_vector_xy", "out_support_vector_z_tangent_vector_x",
                                             "out_tangent_vector_yz", "out_classification",    "out_atom_index"};
        if (!compute_spline_program) compute_spline_program = glCreateProgram();
        gl::attach_link_detach_with_transform_feedback(compute_spline_program, shaders, feedback_varyings, GL_INTERLEAVED_ATTRIBS);
    }

    {
        GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/backbone_draw_spline.vert", GL_VERTEX_SHADER);
        GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/backbone_draw_spline.geom", GL_GEOMETRY_SHADER);
        GLuint f_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/backbone_draw_spline.frag", GL_FRAGMENT_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(g_shader);
            glDeleteShader(f_shader);
        };

        const GLuint shaders[] = {v_shader, g_shader, f_shader};
        if (!draw_spline_program) draw_spline_program = glCreateProgram();
        gl::attach_link_detach(draw_spline_program, shaders);
    }
}

void shutdown() {
    if (extract_control_points_program) glDeleteProgram(extract_control_points_program);
    if (compute_spline_program) glDeleteProgram(compute_spline_program);
    if (draw_spline_program) glDeleteProgram(draw_spline_program);
}
}  // namespace backbone_spline

namespace pbc_view_velocity {
static GLuint program = 0;
void initialize() {
    GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/compute_pbc_view_velocity.vert", GL_VERTEX_SHADER);
    defer { glDeleteShader(v_shader); };

    const GLuint shaders[] = {v_shader};
    const GLchar* feedback_varyings[] = {"out_view_velocity"};
    if (!program) program = glCreateProgram();
    gl::attach_link_detach_with_transform_feedback(program, shaders, feedback_varyings, GL_INTERLEAVED_ATTRIBS);
}

void shutdown() {
    if (program) glDeleteProgram(program);
}
}  // namespace pbc_view_velocity

namespace lean_and_mean {
namespace vdw {
static GLuint program_persp = 0;
static GLuint program_ortho = 0;

static void initialize() {
    {
        GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/lean_and_mean/vdw.vert", GL_VERTEX_SHADER);
        GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/lean_and_mean/vdw.geom", GL_GEOMETRY_SHADER);
        GLuint f_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/lean_and_mean/vdw.frag", GL_FRAGMENT_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(g_shader);
            glDeleteShader(f_shader);
        };

        if (!program_persp) program_persp = glCreateProgram();
        const GLuint shaders[] = {v_shader, g_shader, f_shader};
        gl::attach_link_detach(program_persp, shaders);
    }
    {
        GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/lean_and_mean/vdw_ortho.vert", GL_VERTEX_SHADER);
        GLuint f_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/lean_and_mean/vdw_ortho.frag", GL_FRAGMENT_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(f_shader);
        };

        if (!program_ortho) program_ortho = glCreateProgram();
        const GLuint shaders[] = {v_shader, f_shader};
        gl::attach_link_detach(program_ortho, shaders);
    }
}

static void shutdown() {
    if (program_persp) glDeleteProgram(program_persp);
    if (program_ortho) glDeleteProgram(program_ortho);
}
}  // namespace vdw

namespace licorice {
static GLuint program = 0;

static void initialize() {
    GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/lean_and_mean/licorice.vert", GL_VERTEX_SHADER);
    GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/lean_and_mean/licorice.geom", GL_GEOMETRY_SHADER);
    GLuint f_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/lean_and_mean/licorice.frag", GL_FRAGMENT_SHADER);
    defer {
        glDeleteShader(v_shader);
        glDeleteShader(g_shader);
        glDeleteShader(f_shader);
    };

    if (!program) program = glCreateProgram();
    const GLuint shaders[] = {v_shader, g_shader, f_shader};
    gl::attach_link_detach(program, shaders);
}

static void shutdown() {
    if (program) glDeleteProgram(program);
}

}  // namespace licorice

namespace ribbon {
static GLuint program = 0;

void intitialize() {
    GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/lean_and_mean/ribbons.vert", GL_VERTEX_SHADER);
    GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/lean_and_mean/ribbons.geom", GL_GEOMETRY_SHADER);
    GLuint f_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/lean_and_mean/ribbons.frag", GL_FRAGMENT_SHADER);
    defer {
        glDeleteShader(v_shader);
        glDeleteShader(g_shader);
        glDeleteShader(f_shader);
    };

    if (!program) program = glCreateProgram();
    const GLuint shaders[] = {v_shader, g_shader, f_shader};
    gl::attach_link_detach(program, shaders);
}

void shutdown() {
    if (program) glDeleteProgram(program);
}

}  // namespace ribbon
}  // namespace lean_and_mean

namespace culling {
namespace aabb {
static GLuint program_compute = 0;
static GLuint program_draw = 0;
static GLuint program_cull = 0;
static GLuint program_fill_draw_multi_arrays_cmd = 0;
static GLuint program_fill_draw_cmd = 0;
static GLuint indirect_cmd_buffer = 0;
static GLuint indirect_cmd_count = 0;

void intitialize() {
    {
        GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/culling/compute_aabb.vert", GL_VERTEX_SHADER);
        defer { glDeleteShader(v_shader); };

        if (!program_compute) program_compute = glCreateProgram();
        const GLchar* feedback_varyings[] = {"out_aabb_center", "out_aabb_extent"};
        const GLuint shaders[] = {v_shader};
        gl::attach_link_detach_with_transform_feedback(program_compute, shaders, feedback_varyings, GL_INTERLEAVED_ATTRIBS);
    }
    {
        GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/culling/draw_aabb.vert", GL_VERTEX_SHADER);
        GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/culling/draw_aabb.geom", GL_GEOMETRY_SHADER);
        GLuint f_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/culling/draw_aabb.frag", GL_FRAGMENT_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(g_shader);
            glDeleteShader(f_shader);
        };

        if (!program_draw) program_draw = glCreateProgram();
        const GLuint shaders[] = {v_shader, g_shader, f_shader};
        gl::attach_link_detach(program_draw, shaders);
    }
    {
        GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/culling/draw_aabb.vert", GL_VERTEX_SHADER);
        GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/culling/draw_aabb.geom", GL_GEOMETRY_SHADER);
        GLuint f_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/culling/cull_aabb.frag", GL_FRAGMENT_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(g_shader);
            glDeleteShader(f_shader);
        };

        if (!program_cull) program_cull = glCreateProgram();
        const GLuint shaders[] = {v_shader, g_shader, f_shader};
        gl::attach_link_detach(program_cull, shaders);
    }
    {
        GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/culling/fill_indirect_draw_elements.vert", GL_VERTEX_SHADER);
        defer { glDeleteShader(v_shader); };

        if (!program_fill_draw_cmd) program_fill_draw_cmd = glCreateProgram();
        const GLuint shaders[] = {v_shader};
        gl::attach_link_detach(program_fill_draw_cmd, shaders);
    }
    {
        GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/culling/fill_indirect_draw_multi_arrays.vert", GL_VERTEX_SHADER);
        defer { glDeleteShader(v_shader); };

        if (!program_fill_draw_multi_arrays_cmd) program_fill_draw_multi_arrays_cmd = glCreateProgram();
        const GLuint shaders[] = {v_shader};
        gl::attach_link_detach(program_fill_draw_multi_arrays_cmd, shaders);
    }
    if (!indirect_cmd_buffer) {
        glGenBuffers(1, &indirect_cmd_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, indirect_cmd_buffer);
        glBufferData(GL_ARRAY_BUFFER, MEGABYTES(1), NULL, GL_DYNAMIC_COPY);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    if (!indirect_cmd_count) {
        glGenBuffers(1, &indirect_cmd_count);
        glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, indirect_cmd_count);
        glBufferData(GL_ATOMIC_COUNTER_BUFFER, sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
        glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);
    }
}

void shutdown() {
    if (program_compute) glDeleteProgram(program_compute);
    if (program_draw) glDeleteProgram(program_draw);
    if (program_cull) glDeleteProgram(program_cull);
    if (indirect_cmd_buffer) glDeleteBuffers(1, &indirect_cmd_buffer);
    if (indirect_cmd_count) glDeleteBuffers(1, &indirect_cmd_count);
}

}  // namespace aabb
}  // namespace culling

namespace scan {
static struct {
    GLuint prefixsum = 0;
    GLuint offsets = 0;
    GLuint combine = 0;
} program;

static GLuint offset_buffer = 0;
static uvec3 max_work_group_count = {0, 0, 0};

void initialize() {
    {
        GLuint shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/sdf/scan.comp", GL_COMPUTE_SHADER, "#define TASK TASK_SUM");
        defer { glDeleteShader(shader); };

        if (!program.prefixsum) program.prefixsum = glCreateProgram();
        const GLuint shaders[] = {shader};
        gl::attach_link_detach(program.prefixsum, shaders);
    }
    {
        GLuint shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/sdf/scan.comp", GL_COMPUTE_SHADER, "#define TASK TASK_OFFSETS");
        defer { glDeleteShader(shader); };

        if (!program.offsets) program.offsets = glCreateProgram();
        const GLuint shaders[] = {shader};
        gl::attach_link_detach(program.offsets, shaders);
    }
    {
        GLuint shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/sdf/scan.comp", GL_COMPUTE_SHADER, "#define TASK TASK_COMBINE");
        defer { glDeleteShader(shader); };

        if (!program.combine) program.combine = glCreateProgram();
        const GLuint shaders[] = {shader};
        gl::attach_link_detach(program.combine, shaders);
    }

    if (!offset_buffer) {
        glGenBuffers(1, &offset_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, offset_buffer);
        glBufferData(GL_ARRAY_BUFFER, MEGABYTES(2), 0, GL_DYNAMIC_COPY);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 0, (GLint*)&max_work_group_count[0]);
}

void shutdown() {
    glDeleteProgram(program.prefixsum);
    glDeleteProgram(program.offsets);
    glDeleteProgram(program.combine);
}

static void scan(const GLuint input, const GLuint output, GLuint elements);
}  // namespace scan

namespace sdf {

static GLuint program_draw = 0;
static GLuint program_create_work_indirect = 0;
static GLuint program_distance_spheres = 0;
static GLuint cube_vbo = 0;
static GLuint ubo = 0;

static struct {
    GLuint program_bin = 0;
    GLuint program_sort = 0;

    GLuint buffer_cell_offset = 0;
    GLuint buffer_cell_count = 0;
    GLuint buffer_global_cell_idx = 0;
    GLuint buffer_local_cell_idx = 0;
    GLuint buffer_sphere = 0;

    AABB aabb = {};
    uvec3 dim = {0, 0, 0};
    uint32_t sphere_count = 0;
    uint32_t num_cells = 0;
} spatial_hash;

static struct {
    GLuint texture = 0;
    AABB aabb = {};
    uvec3 dim = {0, 0, 0};
} sdf_volume;

#define SPHERE_BINNING_GROUP_COUNT 512
#define SPHERE_WRITE_COMPRESSED_GROUP_COUNT 512

void initialize() {
    if (!sdf_volume.texture) {
        glGenTextures(1, &sdf_volume.texture);
    }

    if (!cube_vbo) {
        // https://stackoverflow.com/questions/28375338/cube-using-single-gl-triangle-strip
        const uint8_t cube_strip[42] = {0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1,
                                        0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0};
        glGenBuffers(1, &cube_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, cube_vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(cube_strip), cube_strip, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    if (!ubo) {
        glGenBuffers(1, &ubo);
        glBindBuffer(GL_ARRAY_BUFFER, ubo);
        glBufferData(GL_ARRAY_BUFFER, KILOBYTES(1), 0, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    if (!spatial_hash.buffer_cell_offset) {
        glGenBuffers(1, &spatial_hash.buffer_cell_offset);
    }

    if (!spatial_hash.buffer_cell_count) {
        glGenBuffers(1, &spatial_hash.buffer_cell_count);
    }

    if (!spatial_hash.buffer_sphere) {
        glGenBuffers(1, &spatial_hash.buffer_sphere);
    }

    if (!spatial_hash.buffer_global_cell_idx) {
        glGenBuffers(1, &spatial_hash.buffer_global_cell_idx);
    }

    if (!spatial_hash.buffer_local_cell_idx) {
        glGenBuffers(1, &spatial_hash.buffer_local_cell_idx);
    }

    {
        GLuint shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/sdf/bin.comp", GL_COMPUTE_SHADER,
                                                     "#define GROUP_SIZE " TOSTRING(SPHERE_BINNING_GROUP_COUNT));
        defer { glDeleteShader(shader); };

        if (!spatial_hash.program_bin) spatial_hash.program_bin = glCreateProgram();
        const GLuint shaders[] = {shader};
        gl::attach_link_detach(spatial_hash.program_bin, shaders);
    }
    {
        GLuint shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/sdf/sort.comp", GL_COMPUTE_SHADER,
                                                     "#define GROUP_SIZE " TOSTRING(SPHERE_WRITE_COMPRESSED_GROUP_COUNT));
        defer { glDeleteShader(shader); };

        if (!spatial_hash.program_sort) spatial_hash.program_sort = glCreateProgram();
        const GLuint shaders[] = {shader};
        gl::attach_link_detach(spatial_hash.program_sort, shaders);
    }
    {
        GLuint shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/sdf/fill_populated_voxels_indirect.comp", GL_COMPUTE_SHADER);
        defer { glDeleteShader(shader); };

        if (!program_create_work_indirect) program_create_work_indirect = glCreateProgram();
        const GLuint shaders[] = {shader};
        gl::attach_link_detach(program_create_work_indirect, shaders);
    }
    {
        GLuint shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/sdf/write_distance_spheres.comp", GL_COMPUTE_SHADER);
        defer { glDeleteShader(shader); };

        if (!program_distance_spheres) program_distance_spheres = glCreateProgram();
        const GLuint shaders[] = {shader};
        gl::attach_link_detach(program_distance_spheres, shaders);
    }
    /*
    {
        GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/sdf/raycaster.vert", GL_VERTEX_SHADER);
        GLuint f_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/sdf/raycaster.frag", GL_FRAGMENT_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(f_shader);
        };

        if (!program_draw) program_draw = glCreateProgram();
        const GLuint shaders[] = {v_shader, f_shader};
        gl::attach_link_detach(program_draw, shaders);
    }
    */
}

void shutdown() {}

}  // namespace sdf

void initialize() {
    if (!major_version) glGetIntegerv(GL_MAJOR_VERSION, &major_version);
    if (!minor_version) glGetIntegerv(GL_MINOR_VERSION, &minor_version);

    if (major_version < 3 && minor_version < 3) {
        LOG_ERROR("Molecule draw requires version OpenGL 3.3 or higher.");
        return;
    }

    if (!vao) glGenVertexArrays(1, &vao);
    if (!tex[0]) glGenTextures(4, tex);

    vdw::initialize();
    licorice::initialize();
    ribbon::intitialize();
    cartoon::intitialize();
    backbone_spline::initialize();
    lean_and_mean::vdw::initialize();
    lean_and_mean::licorice::initialize();
    lean_and_mean::ribbon::intitialize();
    pbc_view_velocity::initialize();

    if (major_version >= 4 && minor_version >= 3) {
        culling::aabb::intitialize();
        sdf::initialize();
        scan::initialize();
    }
}

void shutdown() {
    if (vao) glDeleteVertexArrays(1, &vao);
    if (tex[0]) glDeleteTextures(4, tex);
    vdw::shutdown();
    licorice::shutdown();
    ribbon::shutdown();
    cartoon::shutdown();
    backbone_spline::shutdown();
    lean_and_mean::vdw::shutdown();
    lean_and_mean::licorice::shutdown();
    lean_and_mean::ribbon::shutdown();
    pbc_view_velocity::shutdown();
    culling::aabb::shutdown();
    sdf::shutdown();
    scan::shutdown();
}

void draw_vdw(GLuint atom_position_buffer, GLuint atom_radius_buffer, GLuint atom_color_buffer, GLuint atom_view_velocity_buffer, i32 atom_count,
              const ViewParam& view_param, float radius_scale) {
    ASSERT(glIsBuffer(atom_position_buffer));
    ASSERT(glIsBuffer(atom_radius_buffer));
    ASSERT(glIsBuffer(atom_color_buffer));
    ASSERT(glIsBuffer(atom_view_velocity_buffer));

    const mat4 curr_view_to_prev_clip_mat = view_param.matrix.previous.view_proj_jittered * view_param.matrix.inverse.view;
    const vec2 res = view_param.resolution;
    const vec4 jitter_uv = vec4(view_param.jitter.current / res, view_param.jitter.previous / res);
    static uint32_t frame = 0;
    frame++;

    // Depth pre-pass
#if 0
    {
        PUSH_GPU_SECTION("VDW PREPASS")
        const GLuint prog = vdw::program_depth_prepass;
        const mat4 mvp_mat = view_param.matrix.current.view_proj_jittered;

        glUseProgram(prog);

        glUniform1i(glGetUniformLocation(prog, "u_buf_pos"), 0);
        glUniform1i(glGetUniformLocation(prog, "u_buf_rad"), 1);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_mvp_mat"), 1, GL_FALSE, &mvp_mat[0][0]);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_BUFFER, tex[0]);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, atom_position_buffer);
    
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_BUFFER, tex[1]);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, atom_radius_buffer);

        glBindVertexArray(vao);
        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, vdw::vbo_icosa);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (const GLvoid*)0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vdw::ibo_icosa);

        glDepthFunc(GL_LESS);
        glColorMask(0, 0, 0, 0); // Output nothing but depth
        glDepthMask(1);

        glDrawElementsInstanced(GL_TRIANGLES, 60, GL_UNSIGNED_BYTE, 0, atom_count);

        glColorMask(1, 1, 1, 1);
        POP_GPU_SECTION()
    }
#endif

    const bool ortho = is_orthographic_proj_matrix(view_param.matrix.current.proj);
    GLuint prog = ortho ? vdw::program_ortho : vdw::program_persp;
    glUseProgram(prog);

    if (ortho) {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_BUFFER, tex[0]);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, atom_position_buffer);

        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_BUFFER, tex[1]);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, atom_radius_buffer);

        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_BUFFER, tex[2]);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color_buffer);

        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_BUFFER, tex[3]);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, atom_view_velocity_buffer);

        glBindBuffer(GL_TEXTURE_BUFFER, 0);
        glActiveTexture(GL_TEXTURE0);

        // Uniforms
        glUniform1i(glGetUniformLocation(prog, "u_buf_position"), 0);
        glUniform1i(glGetUniformLocation(prog, "u_buf_radius"), 1);
        glUniform1i(glGetUniformLocation(prog, "u_buf_color"), 2);
        glUniform1i(glGetUniformLocation(prog, "u_buf_view_velocity"), 3);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_mat"), 1, GL_FALSE, &view_param.matrix.current.view[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.proj_jittered[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_curr_view_to_prev_clip_mat"), 1, GL_FALSE, &curr_view_to_prev_clip_mat[0][0]);
        glUniform1f(glGetUniformLocation(prog, "u_radius_scale"), radius_scale);
        glUniform4fv(glGetUniformLocation(prog, "u_jitter_uv"), 1, &jitter_uv[0]);

        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLES, 0, atom_count * 3);
        glBindVertexArray(0);
    } else {
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, atom_position_buffer);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(AtomPosition), (const GLvoid*)0);
        glEnableVertexAttribArray(0);

        glBindBuffer(GL_ARRAY_BUFFER, atom_radius_buffer);
        glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(AtomRadius), (const GLvoid*)0);
        glEnableVertexAttribArray(1);

        glBindBuffer(GL_ARRAY_BUFFER, atom_color_buffer);
        glVertexAttribPointer(2, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(AtomColor), (const GLvoid*)0);
        glEnableVertexAttribArray(2);

        glBindBuffer(GL_ARRAY_BUFFER, atom_view_velocity_buffer);
        glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(AtomVelocity), (const GLvoid*)0);
        glEnableVertexAttribArray(3);

        glBindBuffer(GL_ARRAY_BUFFER, 0);

        // Uniforms
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_mat"), 1, GL_FALSE, &view_param.matrix.current.view[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.proj_jittered[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_inv_proj_mat"), 1, GL_FALSE, &view_param.matrix.inverse.proj_jittered[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_curr_view_to_prev_clip_mat"), 1, GL_FALSE, &curr_view_to_prev_clip_mat[0][0]);
        glUniform1f(glGetUniformLocation(prog, "u_radius_scale"), radius_scale);
        glUniform4fv(glGetUniformLocation(prog, "u_jitter_uv"), 1, &jitter_uv[0]);
        glUniform1ui(glGetUniformLocation(prog, "u_frame"), frame);

        glDrawArrays(GL_POINTS, 0, atom_count);
        glBindVertexArray(0);
    }
    glUseProgram(0);
}

void draw_licorice(GLuint atom_position_buffer, GLuint atom_color_buffer, GLuint atom_velocity_buffer, GLuint bond_buffer, i32 bond_count,
                   const ViewParam& view_param, float radius_scale) {
    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, atom_position_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(AtomPosition), (const GLvoid*)0);

    glBindBuffer(GL_ARRAY_BUFFER, atom_color_buffer);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(AtomColor), (const GLvoid*)0);

    glBindBuffer(GL_ARRAY_BUFFER, atom_velocity_buffer);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(AtomVelocity), (const GLvoid*)0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bond_buffer);

    const mat4 curr_view_to_prev_clip_mat = view_param.matrix.previous.view_proj_jittered * view_param.matrix.inverse.view;
    const vec2 res = view_param.resolution;
    const vec4 jitter_uv = vec4(view_param.jitter.current / res, view_param.jitter.previous / res);

    const bool ortho = is_orthographic_proj_matrix(view_param.matrix.current.proj_jittered);
    const GLuint prog = ortho ? licorice::program_ortho : licorice::program_persp;
    glUseProgram(prog);

    glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_mat"), 1, GL_FALSE, &view_param.matrix.current.view[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.proj_jittered[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_curr_view_to_prev_clip_mat"), 1, GL_FALSE, &curr_view_to_prev_clip_mat[0][0]);
    glUniform1f(glGetUniformLocation(prog, "u_radius"), 0.25f * radius_scale);
    glUniform4fv(glGetUniformLocation(prog, "u_jitter_uv"), 1, &jitter_uv[0]);

    glDrawElements(GL_LINES, bond_count * 2, GL_UNSIGNED_INT, 0);
    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void compute_backbone_control_points(GLuint dst_buffer, GLuint atom_position_buffer, GLuint backbone_index_buffer, int num_backbone_indices,
                                     GLuint ramachandran_tex) {
    glEnable(GL_RASTERIZER_DISCARD);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, atom_position_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(AtomPosition), (const GLvoid*)0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, backbone_index_buffer);

    glUseProgram(backbone_spline::extract_control_points_program);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, ramachandran_tex);
    //@TODO: Replace ramachandran texture to some vertex attribute which holds the type of secondary structure

    glUniform1i(glGetUniformLocation(backbone_spline::extract_control_points_program, "u_ramachandran_tex"), 0);

    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, dst_buffer);
    glBeginTransformFeedback(GL_POINTS);
    glDrawElements(GL_TRIANGLES_ADJACENCY, num_backbone_indices, GL_UNSIGNED_INT, 0);
    glEndTransformFeedback();

    glUseProgram(0);

    glDisable(GL_RASTERIZER_DISCARD);
}

void compute_backbone_spline(GLuint dst_buffer, GLuint control_point_buffer, GLuint control_point_index_buffer, int num_control_point_indices,
                             float tension) {
    glEnable(GL_RASTERIZER_DISCARD);
    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFFU);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, control_point_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, position));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, support_vector));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, classification));
    glEnableVertexAttribArray(3);
    glVertexAttribIPointer(3, 1, GL_UNSIGNED_INT, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, atom_index));

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, control_point_index_buffer);

    glUseProgram(backbone_spline::compute_spline_program);

    glUniform1f(glGetUniformLocation(backbone_spline::compute_spline_program, "u_tension"), tension);

    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, dst_buffer);
    glBeginTransformFeedback(GL_POINTS);
    glDrawElements(GL_LINE_STRIP_ADJACENCY, num_control_point_indices, GL_UNSIGNED_INT, 0);
    glEndTransformFeedback();

    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);
    glDisable(GL_RASTERIZER_DISCARD);
}

void compute_pbc_view_velocity(GLuint dst_buffer, GLuint position_buffer, GLuint old_position_buffer, int count, const ViewParam& view_param,
                               const vec3& box_ext) {
    glEnable(GL_RASTERIZER_DISCARD);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, position_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(AtomPosition), (const GLvoid*)0);

    glBindBuffer(GL_ARRAY_BUFFER, old_position_buffer);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(AtomPosition), (const GLvoid*)0);

    glUseProgram(pbc_view_velocity::program);

    glUniformMatrix4fv(glGetUniformLocation(pbc_view_velocity::program, "u_view_mat"), 1, GL_FALSE, &view_param.matrix.current.view[0][0]);
    glUniform3fv(glGetUniformLocation(pbc_view_velocity::program, "u_box_ext"), 1, &box_ext[0]);

    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, dst_buffer);
    glBeginTransformFeedback(GL_POINTS);
    glDrawArrays(GL_POINTS, 0, count);
    glEndTransformFeedback();

    glUseProgram(0);

    glDisable(GL_RASTERIZER_DISCARD);
}

void draw_spline(GLuint spline_buffer, GLuint spline_index_buffer, i32 num_spline_indices, const ViewParam& view_param, uint32_t s_color,
                 uint32_t v_color, uint32_t t_color) {
    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFFU);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, spline_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, position));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, support_vector));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, tangent_vector));

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, spline_index_buffer);

    GLuint prog = backbone_spline::draw_spline_program;
    glUseProgram(prog);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_proj_mat"), 1, GL_FALSE, &(view_param.matrix.current.view_proj_jittered)[0][0]);
    glUniform4fv(glGetUniformLocation(prog, "u_s_color"), 1, &math::convert_color(s_color)[0]);
    glUniform4fv(glGetUniformLocation(prog, "u_v_color"), 1, &math::convert_color(v_color)[0]);
    glUniform4fv(glGetUniformLocation(prog, "u_t_color"), 1, &math::convert_color(t_color)[0]);

    glDrawElements(GL_LINE_STRIP, num_spline_indices, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);
}

void draw_ribbons(GLuint spline_buffer, GLuint spline_index_buffer, GLuint atom_color_buffer, GLuint atom_velocity_buffer, i32 num_spline_indices,
                  const ViewParam& view_param) {
    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFFU);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, spline_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, position));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, support_vector));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, tangent_vector));
    glEnableVertexAttribArray(3);
    glVertexAttribIPointer(3, 1, GL_UNSIGNED_INT, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, atom_index));

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, spline_index_buffer);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, tex[0]);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color_buffer);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, tex[1]);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, atom_velocity_buffer);

    const vec2 res = view_param.resolution;
    const vec4 jitter_uv = vec4(view_param.jitter.current / res, view_param.jitter.previous / res);

    glUseProgram(ribbon::program);

    // Uniforms
    glUniform1i(glGetUniformLocation(ribbon::program, "u_atom_color_buffer"), 0);
    glUniform1i(glGetUniformLocation(ribbon::program, "u_atom_velocity_buffer"), 1);
    glUniformMatrix4fv(glGetUniformLocation(ribbon::program, "u_normal_mat"), 1, GL_FALSE, &view_param.matrix.current.norm[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(ribbon::program, "u_view_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.view_proj_jittered[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(ribbon::program, "u_prev_view_proj_mat"), 1, GL_FALSE,
                       &view_param.matrix.previous.view_proj_jittered[0][0]);
    glUniform4fv(glGetUniformLocation(ribbon::program, "u_jitter_uv"), 1, &jitter_uv[0]);

    glDrawElements(GL_LINE_STRIP, num_spline_indices, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);
}

void draw_cartoon(GLuint spline_buffer, GLuint spline_index_buffer, GLuint atom_color_buffer, i32 num_spline_indices, const ViewParam& view_param) {
    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFFU);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, spline_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, position));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, support_vector));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, tangent_vector));
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, classification));
    glEnableVertexAttribArray(4);
    glVertexAttribIPointer(4, 1, GL_UNSIGNED_INT, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, atom_index));

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, spline_index_buffer);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, tex[0]);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color_buffer);

    GLuint prog = cartoon::program;
    glUseProgram(prog);
    glUniform1i(glGetUniformLocation(prog, "u_atom_color_buffer"), 0);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_normal_mat"), 1, GL_FALSE, &view_param.matrix.current.norm[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.view_proj_jittered[0][0]);
    glDrawElements(GL_LINE_STRIP, num_spline_indices, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);
}

namespace lean_and_mean {

void draw_vdw(GLuint atom_position_buffer, GLuint atom_radius_buffer, GLuint atom_color_buffer, GLuint atom_mask_buffer, i32 atom_count,
              const ViewParam& view_param, float radius_scale, vec4 color, uint32_t mask) {
    ASSERT(glIsBuffer(atom_position_buffer));
    ASSERT(glIsBuffer(atom_radius_buffer));
    ASSERT(glIsBuffer(atom_color_buffer));
    ASSERT(glIsBuffer(atom_mask_buffer));

    const vec2 res = view_param.resolution;
    const vec4 jitter_uv = vec4(view_param.jitter.current / res, view_param.jitter.previous / res);

    const bool ortho = is_orthographic_proj_matrix(view_param.matrix.current.proj_jittered);
    const GLuint prog = ortho ? vdw::program_ortho : vdw::program_persp;
    glUseProgram(prog);

    if (ortho) {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_BUFFER, tex[0]);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, atom_position_buffer);

        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_BUFFER, tex[1]);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, atom_radius_buffer);

        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_BUFFER, tex[2]);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color_buffer);

        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_BUFFER, tex[3]);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, atom_mask_buffer);

        glBindBuffer(GL_TEXTURE_BUFFER, 0);
        glActiveTexture(GL_TEXTURE0);

        // Uniforms
        glUniform1i(glGetUniformLocation(prog, "u_buf_position"), 0);
        glUniform1i(glGetUniformLocation(prog, "u_buf_radius"), 1);
        glUniform1i(glGetUniformLocation(prog, "u_buf_color"), 2);
        glUniform1i(glGetUniformLocation(prog, "u_buf_mask"), 3);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_mat"), 1, GL_FALSE, &view_param.matrix.current.view[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.proj_jittered[0][0]);
        glUniform1f(glGetUniformLocation(prog, "u_radius_scale"), radius_scale);
        glUniform1ui(glGetUniformLocation(prog, "u_mask"), mask);
        glUniform4fv(glGetUniformLocation(prog, "u_color"), 1, &color[0]);

        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLES, 0, atom_count * 3);
        glBindVertexArray(0);
    } else {
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, atom_position_buffer);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(AtomPosition), (const GLvoid*)0);
        glEnableVertexAttribArray(0);

        glBindBuffer(GL_ARRAY_BUFFER, atom_radius_buffer);
        glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(AtomRadius), (const GLvoid*)0);
        glEnableVertexAttribArray(1);

        glBindBuffer(GL_ARRAY_BUFFER, atom_color_buffer);
        glVertexAttribPointer(2, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(AtomColor), (const GLvoid*)0);
        glEnableVertexAttribArray(2);

        glBindBuffer(GL_ARRAY_BUFFER, atom_mask_buffer);
        glVertexAttribIPointer(3, 1, GL_UNSIGNED_BYTE, sizeof(AtomMask), (const GLvoid*)0);
        glEnableVertexAttribArray(3);

        // Uniforms
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_mat"), 1, GL_FALSE, &view_param.matrix.current.view[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.proj_jittered[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_inv_proj_mat"), 1, GL_FALSE, &view_param.matrix.inverse.proj_jittered[0][0]);
        glUniform4fv(glGetUniformLocation(prog, "u_jitter_uv"), 1, &jitter_uv[0]);
        glUniform4fv(glGetUniformLocation(prog, "u_color"), 1, &color[0]);
        glUniform1f(glGetUniformLocation(prog, "u_radius_scale"), radius_scale);
        glUniform1ui(glGetUniformLocation(prog, "u_mask"), mask);

        glDrawArrays(GL_POINTS, 0, atom_count);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }
    glUseProgram(0);
}

void draw_licorice(GLuint atom_position_buffer, GLuint atom_color_buffer, GLuint atom_mask_buffer, GLuint bond_buffer, i32 bond_count,
                   const ViewParam& view_param, float radius_scale, vec4 color, uint32_t mask) {
    ASSERT(glIsBuffer(atom_position_buffer));
    ASSERT(glIsBuffer(atom_color_buffer));
    ASSERT(glIsBuffer(atom_mask_buffer));
    ASSERT(glIsBuffer(bond_buffer));

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, atom_position_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(AtomPosition), (const GLvoid*)0);

    glBindBuffer(GL_ARRAY_BUFFER, atom_color_buffer);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(AtomColor), (const GLvoid*)0);

    glBindBuffer(GL_ARRAY_BUFFER, atom_mask_buffer);
    glEnableVertexAttribArray(2);
    glVertexAttribIPointer(2, 1, GL_UNSIGNED_BYTE, sizeof(AtomMask), (const GLvoid*)0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bond_buffer);

    GLuint prog = licorice::program;
    glUseProgram(prog);

    glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_mat"), 1, GL_FALSE, &view_param.matrix.current.view[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.proj_jittered[0][0]);
    glUniform4fv(glGetUniformLocation(prog, "u_color"), 1, &color[0]);
    glUniform1ui(glGetUniformLocation(prog, "u_mask"), mask);
    glUniform1f(glGetUniformLocation(prog, "u_radius"), 0.25f * radius_scale);

    glDrawElements(GL_LINES, bond_count * 2, GL_UNSIGNED_INT, 0);
    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void draw_ribbons(GLuint spline_buffer, GLuint spline_index_buffer, GLuint atom_color_buffer, GLuint atom_mask_buffer, i32 num_spline_indices,
                  const ViewParam& view_param, float scale, vec4 color, uint32_t mask) {
    ASSERT(glIsBuffer(spline_buffer));
    ASSERT(glIsBuffer(spline_index_buffer));
    ASSERT(glIsBuffer(atom_color_buffer));
    ASSERT(glIsBuffer(atom_mask_buffer));

    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFFU);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, spline_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, position));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, support_vector));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_SHORT, GL_TRUE, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, tangent_vector));
    glEnableVertexAttribArray(3);
    glVertexAttribIPointer(3, 1, GL_UNSIGNED_INT, sizeof(ControlPoint), (const GLvoid*)offsetof(ControlPoint, atom_index));

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, spline_index_buffer);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, tex[0]);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color_buffer);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, tex[1]);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, atom_mask_buffer);

    const vec2 ribbon_scale = vec2(1.0f, 0.1f) * scale;

    GLuint prog = ribbon::program;
    glUseProgram(prog);

    // Uniforms
    glUniform1i(glGetUniformLocation(prog, "u_atom_color_buffer"), 0);
    glUniform1i(glGetUniformLocation(prog, "u_atom_mask_buffer"), 1);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.view_proj_jittered[0][0]);
    glUniform4fv(glGetUniformLocation(prog, "u_color"), 1, &color[0]);
    glUniform2fv(glGetUniformLocation(prog, "u_scale"), 1, &ribbon_scale[0]);
    glUniform1ui(glGetUniformLocation(prog, "u_mask"), mask);

    glDrawElements(GL_LINE_STRIP, num_spline_indices, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);
}

}  // namespace lean_and_mean

namespace culling {
void compute_residue_aabbs(GLuint aabb_buffer, GLuint atom_position_buffer, GLuint atom_radii_buffer, GLuint atom_range_buffer, int32_t count) {
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, tex[0]);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, atom_position_buffer);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, tex[1]);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, atom_radii_buffer);

    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, atom_range_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribIPointer(0, 1, GL_INT, 8, 0);
    glEnableVertexAttribArray(1);
    glVertexAttribIPointer(1, 1, GL_INT, 8, (const void*)4);

    glUseProgram(aabb::program_compute);

    glUniform1i(glGetUniformLocation(aabb::program_compute, "buf_atom_pos"), 0);
    glUniform1i(glGetUniformLocation(aabb::program_compute, "buf_atom_rad"), 1);

    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, aabb_buffer);
    glBeginTransformFeedback(GL_POINTS);
    glDrawArrays(GL_POINTS, 0, count);
    glEndTransformFeedback();

    glUseProgram(0);
    glBindVertexArray(0);
}

static void initialize_cmd_buffer(int count) {
    // Do we need to resize cmd buffer?
    const int CMD_SIZE = 4;
    GLint current_size = 0;
    glBindBuffer(GL_ARRAY_BUFFER, aabb::indirect_cmd_buffer);
    glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &current_size);
    const GLsizeiptr required_size = count * CMD_SIZE * sizeof(uint32_t);
    if (required_size > current_size) {
        glBufferData(GL_ARRAY_BUFFER, required_size, nullptr, GL_DYNAMIC_COPY);
    }
}

static void clear_cmd_count() {
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, aabb::indirect_cmd_count);
    glClearBufferData(GL_ATOMIC_COUNTER_BUFFER, GL_R32UI, GL_RED_INTEGER, GL_UNSIGNED_INT, 0);
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);
}

static uint32_t fill_cmd_buffer(GLuint visibility_buffer, int count) {
    initialize_cmd_buffer(count);
    clear_cmd_count();

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, visibility_buffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, aabb::indirect_cmd_buffer);
    glBindBufferBase(GL_ATOMIC_COUNTER_BUFFER, 2, aabb::indirect_cmd_count);

    /*
    typedef  struct {
        uint  count;
        uint  instanceCount;
        uint  first;
        uint  baseInstance;
    } DrawArraysIndirectCommand;
    */
    glBindVertexArray(vao);

    // Fill indirect_cmd_buffer
    glUseProgram(aabb::program_fill_draw_cmd);
    glDrawArrays(GL_POINTS, 0, count);
    glUseProgram(0);

    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

    glBindVertexArray(0);

    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, aabb::indirect_cmd_count);
    uint32_t cmd_count;
    glGetBufferSubData(GL_ATOMIC_COUNTER_BUFFER, 0, 4, &cmd_count);
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);

    return cmd_count;
}

static uint32_t fill_cmd_buffer(GLuint visibility_buffer, GLuint chunk_buffer, int chunk_count) {
    initialize_cmd_buffer(chunk_count);
    clear_cmd_count();

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, visibility_buffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, aabb::indirect_cmd_buffer);
    glBindBufferBase(GL_ATOMIC_COUNTER_BUFFER, 2, aabb::indirect_cmd_count);

    /*
    typedef  struct {
        uint  count;
        uint  instanceCount;
        uint  first;
        uint  baseInstance;
    } DrawArraysIndirectCommand;
    */
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, chunk_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribIPointer(0, 1, GL_UNSIGNED_INT, 8, 0);
    glEnableVertexAttribArray(1);
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_INT, 8, (const void*)4);

    // Fill indirect_cmd_buffer
    glUseProgram(aabb::program_fill_draw_multi_arrays_cmd);
    glDrawArrays(GL_POINTS, 0, chunk_count);
    glUseProgram(0);

    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glBindVertexArray(0);

    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, aabb::indirect_cmd_count);
    uint32_t cmd_count;
    glGetBufferSubData(GL_ATOMIC_COUNTER_BUFFER, 0, 4, &cmd_count);
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);

    return cmd_count;
}

void draw_aabbs(GLuint aabb_buffer, const ViewParam& view_param, int32_t count) {
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, aabb_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 24, 0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 24, (const void*)12);

    glUseProgram(aabb::program_draw);

    glUniformMatrix4fv(glGetUniformLocation(aabb::program_draw, "u_view_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.view_proj_jittered[0][0]);
    glDrawArrays(GL_POINTS, 0, count);

    glUseProgram(0);
    glBindVertexArray(0);
}

void draw_culled_aabbs(GLuint visibility_buffer, GLuint aabb_buffer, const ViewParam& view_param, GLsizei count) {
    fill_cmd_buffer(visibility_buffer, count);
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, aabb::indirect_cmd_count);
    uint32_t draw_size;
    glGetBufferSubData(GL_ATOMIC_COUNTER_BUFFER, 0, 4, &draw_size);
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);

    // printf("visible count %u\n", draw_size);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, aabb_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 24, 0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 24, (const void*)12);

    glUseProgram(aabb::program_draw);
    glUniformMatrix4fv(glGetUniformLocation(aabb::program_draw, "u_view_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.view_proj_jittered[0][0]);

    glBindBuffer(GL_DRAW_INDIRECT_BUFFER, aabb::indirect_cmd_buffer);
    glMultiDrawArraysIndirect(GL_POINTS, 0, draw_size, 16);
    glBindBuffer(GL_DRAW_INDIRECT_BUFFER, 0);
    glUseProgram(0);
    glBindVertexArray(0);
}

void cull_aabbs(GLuint visibility_buffer, GLuint aabb_buffer, const ViewParam& view_param, GLsizei count) {
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, aabb_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 24, 0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 24, (const void*)12);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, visibility_buffer);
    glClearBufferData(GL_SHADER_STORAGE_BUFFER, GL_R32I, GL_RED, GL_INT, 0);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, visibility_buffer);
    glUseProgram(aabb::program_cull);

    glUniformMatrix4fv(glGetUniformLocation(aabb::program_cull, "u_view_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.view_proj_jittered[0][0]);
    glDrawArrays(GL_POINTS, 0, count);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    glUseProgram(0);
    glBindVertexArray(0);
}

void draw_vdw_indirect(GLuint atom_position_buffer, GLuint atom_radius_buffer, GLuint atom_color_buffer, GLuint atom_view_velocity_buffer,
                       GLuint indirect_cmd_buffer, i32 cmd_count, const ViewParam& view_param, float radius_scale) {
    ASSERT(glIsBuffer(atom_position_buffer));
    ASSERT(glIsBuffer(atom_radius_buffer));
    ASSERT(glIsBuffer(atom_color_buffer));
    ASSERT(glIsBuffer(atom_view_velocity_buffer));
    ASSERT(glIsBuffer(indirect_cmd_buffer));

    const mat4 curr_view_to_prev_clip_mat = view_param.matrix.previous.view_proj_jittered * view_param.matrix.inverse.view;
    const vec2 res = view_param.resolution;
    const vec4 jitter_uv = vec4(view_param.jitter.current / res, view_param.jitter.previous / res);
    static uint32_t frame = 0;
    frame++;

    const bool ortho = is_orthographic_proj_matrix(view_param.matrix.current.proj);
    GLuint prog = ortho ? vdw::program_ortho : vdw::program_persp;
    glUseProgram(prog);

    if (ortho) {
        /*
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_BUFFER, tex[0]);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, atom_position_buffer);

        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_BUFFER, tex[1]);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, atom_radius_buffer);

        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_BUFFER, tex[2]);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color_buffer);

        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_BUFFER, tex[3]);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, atom_view_velocity_buffer);

        glBindBuffer(GL_TEXTURE_BUFFER, 0);
        glActiveTexture(GL_TEXTURE0);

        // Uniforms
        glUniform1i(glGetUniformLocation(prog, "u_buf_position"), 0);
        glUniform1i(glGetUniformLocation(prog, "u_buf_radius"), 1);
        glUniform1i(glGetUniformLocation(prog, "u_buf_color"), 2);
        glUniform1i(glGetUniformLocation(prog, "u_buf_view_velocity"), 3);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_mat"), 1, GL_FALSE, &view_param.matrix.current.view[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.proj_jittered[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_curr_view_to_prev_clip_mat"), 1, GL_FALSE, &curr_view_to_prev_clip_mat[0][0]);
        glUniform1f(glGetUniformLocation(prog, "u_radius_scale"), radius_scale);
        glUniform4fv(glGetUniformLocation(prog, "u_jitter_uv"), 1, &jitter_uv[0]);

        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLES, 0, atom_count * 3);
        glBindVertexArray(0);
        */
    } else {
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, atom_position_buffer);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(AtomPosition), (const GLvoid*)0);
        glEnableVertexAttribArray(0);

        glBindBuffer(GL_ARRAY_BUFFER, atom_radius_buffer);
        glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(AtomRadius), (const GLvoid*)0);
        glEnableVertexAttribArray(1);

        glBindBuffer(GL_ARRAY_BUFFER, atom_color_buffer);
        glVertexAttribPointer(2, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(AtomColor), (const GLvoid*)0);
        glEnableVertexAttribArray(2);

        glBindBuffer(GL_ARRAY_BUFFER, atom_view_velocity_buffer);
        glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(AtomVelocity), (const GLvoid*)0);
        glEnableVertexAttribArray(3);

        glBindBuffer(GL_ARRAY_BUFFER, 0);

        // Uniforms
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_mat"), 1, GL_FALSE, &view_param.matrix.current.view[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.proj_jittered[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_inv_proj_mat"), 1, GL_FALSE, &view_param.matrix.inverse.proj_jittered[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_curr_view_to_prev_clip_mat"), 1, GL_FALSE, &curr_view_to_prev_clip_mat[0][0]);
        glUniform1f(glGetUniformLocation(prog, "u_radius_scale"), radius_scale);
        glUniform4fv(glGetUniformLocation(prog, "u_jitter_uv"), 1, &jitter_uv[0]);
        glUniform1ui(glGetUniformLocation(prog, "u_frame"), frame);

        glBindBuffer(GL_DRAW_INDIRECT_BUFFER, indirect_cmd_buffer);
        glMultiDrawArraysIndirect(GL_POINTS, 0, cmd_count, 16);
        glBindBuffer(GL_DRAW_INDIRECT_BUFFER, 0);

        glBindVertexArray(0);
    }
    glUseProgram(0);
}

void draw_culled_vdw(GLuint atom_position_buffer, GLuint atom_radius_buffer, GLuint atom_color_buffer, GLuint atom_view_velocity_buffer,
                     GLuint chunk_visibility_buffer, GLuint chunk_buffer, i32 chunk_count, const ViewParam& view_param, float radius_scale) {
    const uint32_t cmd_count = fill_cmd_buffer(chunk_visibility_buffer, chunk_buffer, chunk_count);

    draw_vdw_indirect(atom_position_buffer, atom_radius_buffer, atom_color_buffer, atom_view_velocity_buffer, aabb::indirect_cmd_buffer, cmd_count,
                      view_param, radius_scale);
}

}  // namespace culling

namespace sdf {

#include <core/spatial_hash.h>

template <typename T>
T div_up(T x, T div) {
    return (x + div - T(1)) / div;
}

struct OptimalSDFVolume {
    AABB aabb {};
    uvec3 dim {};
};

OptimalSDFVolume compute_optimal_sdf_volume(const AABB& aabb, uint32_t max_dim) {
    OptimalSDFVolume vol;
    const vec3 ext = aabb.ext();
    const float max_ext = math::max(ext.x, math::max(ext.y, ext.z));

    vec3 out_ext = vec3(max_ext);
    uvec3 out_dim = uvec3(max_dim);
    for (int i = 0; i < 3; i++) {
        if (ext[i] == max_ext) continue;
        float half_ext = max_ext * 0.5f;
        while (ext[i] < half_ext) {
            out_ext[i] = half_ext;
            out_dim[i] = out_dim[i] >> 2;
            half_ext = half_ext * 0.5f;
        }
    }

    const vec3 center = (aabb.min + aabb.max) * 0.5f;
    const vec3 half_ext = out_ext * 0.5f;
    vol.aabb.min = center - half_ext;
    vol.aabb.max = center + half_ext;
    vol.dim = out_dim;

    return vol;
}

static vec3 calc_voxel_extent(const AABB& aabb, const uvec3& dim) { return aabb.ext() / vec3(dim); }

static u32 calc_num_voxels(const uvec3& dim) { return dim.x * dim.y * dim.z; }

static int8_t encode_distance(float d, float voxel_ext) {
    const float n = d / (voxel_ext * 4.0f);
    return (int8_t)math::round(math::clamp(n, -1.0f, 1.0f) * 127.0f);
}

#define SDF_TEXTURE_FORMAT GL_R8_SNORM

void init_sdf_volume(const AABB& aabb, uvec3 dim, float max_val) {
    sdf_volume.aabb = aabb;
    if (sdf_volume.dim != dim) {
        sdf_volume.dim = dim;

        // @NOTE: Compute log2 (integer version) to find the amount of mipmaps required
        uint32_t mips = 1;
        {
            uint32_t max_dim = math::max(dim.x, math::max(dim.y, dim.z));
            while (max_dim >>= 1) ++mips;
        }

        if (sdf_volume.texture) glDeleteTextures(1, &sdf_volume.texture);
        glGenTextures(1, &sdf_volume.texture);
        glBindTexture(GL_TEXTURE_3D, sdf_volume.texture);
        glTexStorage3D(GL_TEXTURE_3D, mips, SDF_TEXTURE_FORMAT, dim.x, dim.y, dim.z);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_3D, 0);
    }

    PUSH_GPU_SECTION("Clear SDF Volume")
#if SDF_TEXTURE_FORMAT == GL_R32F
    const float value = 1.0f;
    glClearTexImage(sdf_volume.texture, 0, GL_RED, GL_FLOAT, &value);
#else
    const int8_t value = 127;
    glClearTexImage(sdf_volume.texture, 0, GL_RED, GL_BYTE, &value);
#endif
    POP_GPU_SECTION()
}

void init_spatial_hash(const AABB& aabb, uvec3 dim, uint32_t element_count) {
    spatial_hash.aabb = aabb;
    if (spatial_hash.dim != dim) {
        spatial_hash.dim = dim;
        spatial_hash.num_cells = div_up(dim.x * dim.y * dim.z, 4U) * 4U; // Round upwards to nearest 4, scan requires this alignment since it operates on uvec4
        glBindBuffer(GL_ARRAY_BUFFER, spatial_hash.buffer_cell_offset);
        glBufferData(GL_ARRAY_BUFFER, spatial_hash.num_cells * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
        glBindBuffer(GL_ARRAY_BUFFER, spatial_hash.buffer_cell_count);
        glBufferData(GL_ARRAY_BUFFER, spatial_hash.num_cells * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    if (spatial_hash.sphere_count != element_count) {
        spatial_hash.sphere_count = element_count;
        glBindBuffer(GL_ARRAY_BUFFER, spatial_hash.buffer_sphere);
        glBufferData(GL_ARRAY_BUFFER, element_count * sizeof(vec4), NULL, GL_DYNAMIC_COPY);
        glBindBuffer(GL_ARRAY_BUFFER, spatial_hash.buffer_global_cell_idx);
        glBufferData(GL_ARRAY_BUFFER, element_count * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
        glBindBuffer(GL_ARRAY_BUFFER, spatial_hash.buffer_local_cell_idx);
        glBufferData(GL_ARRAY_BUFFER, element_count * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    PUSH_GPU_SECTION("Clear Spatial Hash Buffers")
    glBindBuffer(GL_ARRAY_BUFFER, spatial_hash.buffer_cell_offset);
    glClearBufferData(GL_ARRAY_BUFFER, GL_R32UI, GL_RED_INTEGER, GL_UNSIGNED_INT, 0);

    glBindBuffer(GL_ARRAY_BUFFER, spatial_hash.buffer_cell_count);
    glClearBufferData(GL_ARRAY_BUFFER, GL_R32UI, GL_RED_INTEGER, GL_UNSIGNED_INT, 0);

    glBindBuffer(GL_ARRAY_BUFFER, spatial_hash.buffer_global_cell_idx);
    glClearBufferData(GL_ARRAY_BUFFER, GL_R32UI, GL_RED_INTEGER, GL_UNSIGNED_INT, 0);

    glBindBuffer(GL_ARRAY_BUFFER, spatial_hash.buffer_local_cell_idx);
    glClearBufferData(GL_ARRAY_BUFFER, GL_R32UI, GL_RED_INTEGER, GL_UNSIGNED_INT, 0);

    glBindBuffer(GL_ARRAY_BUFFER, spatial_hash.buffer_sphere);
    glClearBufferData(GL_ARRAY_BUFFER, GL_RGBA32F, GL_RED, GL_FLOAT, 0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    POP_GPU_SECTION()
}

void bin_spheres(const GLuint position_buffer, uint32_t count, const AABB& sphere_aabb) {
    PUSH_GPU_SECTION("Sphere Binning")
    glUseProgram(spatial_hash.program_bin);

    struct Block {
        vec3 aabb_min;
        float _p0;
        vec3 cell_ext;
        float _p1;
        uvec3 cell_dim;
        uint32_t num_spheres;
    } block;

    block.aabb_min = spatial_hash.aabb.min;
    block.cell_ext = calc_voxel_extent(spatial_hash.aabb, spatial_hash.dim);
    block.cell_dim = spatial_hash.dim;
    block.num_spheres = count;

    glBindBuffer(GL_UNIFORM_BUFFER, ubo);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(Block), &block);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, position_buffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, spatial_hash.buffer_local_cell_idx);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, spatial_hash.buffer_cell_count);
    glUniformBlockBinding(spatial_hash.program_bin, 0, 0);

    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    int num_groups = div_up(count, (u32)SPHERE_BINNING_GROUP_COUNT);
    glDispatchCompute(num_groups, 1, 1);

    glUseProgram(0);
    POP_GPU_SECTION()
}

void compute_sphere_offsets() {
    PUSH_GPU_SECTION("Compute Spatial Hash Offsets")
    scan::scan(spatial_hash.buffer_cell_count, spatial_hash.buffer_cell_offset, spatial_hash.num_cells);
    POP_GPU_SECTION()
}

void write_sorted_spheres(const GLuint position_buffer, const GLuint radius_buffer, i32 count, const AABB& sphere_aabb, float max_radius) {
    PUSH_GPU_SECTION("Sort Spheres")
    glUseProgram(spatial_hash.program_sort);

    struct Block {
        vec3 aabb_min;
        float p0;
        vec3 cell_ext;
        float max_radius;
        uvec3 cell_dim;
        uint32_t num_elements;
    } block;

    block.aabb_min = spatial_hash.aabb.min;
    block.cell_ext = calc_voxel_extent(spatial_hash.aabb, spatial_hash.dim);
    block.max_radius = max_radius;
    block.cell_dim = spatial_hash.dim;
    block.num_elements = count;

    glBindBuffer(GL_UNIFORM_BUFFER, ubo);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(Block), &block);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, position_buffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, radius_buffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, spatial_hash.buffer_cell_offset);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, spatial_hash.buffer_cell_count);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, spatial_hash.buffer_local_cell_idx);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, spatial_hash.buffer_sphere);
    glUniformBlockBinding(spatial_hash.program_sort, 0, 0);

    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    int num_groups = div_up(count, SPHERE_WRITE_COMPRESSED_GROUP_COUNT);
    glDispatchCompute(num_groups, 1, 1);

    glUseProgram(0);
    POP_GPU_SECTION()
}

void probe_work() {
    PUSH_GPU_SECTION("Create Work Indirect");
   // glUseProgram(program_create_work_indirect);

    struct Block {
        vec3 voxel_ext;
        vec3 cell_ext;
        uvec3 cell_dim;
        float search_radius;
    } block;

    POP_GPU_SECTION();
}

void write_sphere_distances(GLuint level, float max_sphere_radius, float max_distance) {
    PUSH_GPU_SECTION("Write Distances Spheres");
    glUseProgram(program_distance_spheres);

    struct Block {
        vec3 aabb_min;
        float p0;
        vec3 voxel_ext;
        float max_radius;
        vec3 cell_ext;
        float max_distance;
        uvec3 cell_dim;
        float search_radius;
    } block;

    block.aabb_min = spatial_hash.aabb.min;
    block.voxel_ext = calc_voxel_extent(sdf_volume.aabb, sdf_volume.dim);
    block.max_radius = max_sphere_radius;
    block.cell_ext = calc_voxel_extent(spatial_hash.aabb, spatial_hash.dim);
    block.max_distance = max_distance;
    block.cell_dim = spatial_hash.dim;
    block.search_radius = max_distance + max_sphere_radius;

    glBindBuffer(GL_UNIFORM_BUFFER, ubo);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(Block), &block);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, spatial_hash.buffer_cell_offset);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, spatial_hash.buffer_cell_count);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, spatial_hash.buffer_sphere);
    glUniformBlockBinding(program_distance_spheres, 0, 0);

#if SDF_TEXTURE_FORMAT == GL_R32F
    glBindImageTexture(0, sdf_volume.texture, level, GL_FALSE, 0, GL_WRITE_ONLY, GL_R32F);
#else
    glBindImageTexture(0, sdf_volume.texture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_R8_SNORM);
#endif

    uvec3 num_groups = sdf_volume.dim / 4U;
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glDispatchCompute(num_groups.x, num_groups.y, num_groups.z);
    glUseProgram(0);
    POP_GPU_SECTION();
}

void compute_vdw_sdf(const GLuint atom_pos_buffer, const GLuint atom_rad_buffer, i32 atom_count, const AABB& sphere_aabb, float max_sphere_radius) {
    const uint32_t max_dim = 128U;
    const uint32_t downsample_factor = 32U;
    //const float32 max_sphere_radius = 2.0f;

    auto optimal_vol = compute_optimal_sdf_volume(sphere_aabb, max_dim);
    const vec3 voxel_ext = optimal_vol.aabb.ext() / vec3(optimal_vol.dim);
    const f32 max_distance = voxel_ext.x * 1.0f;
    const f32 spatial_hash_radius = (max_distance + max_sphere_radius) / 1.5f;  // @NOTE: ~search diameter / 3.0 according to nvidias presentation
    const uvec3 spatial_hash_dim = optimal_vol.aabb.ext() / spatial_hash_radius;

    init_sdf_volume(optimal_vol.aabb, optimal_vol.dim, 1.0f);
    init_spatial_hash(optimal_vol.aabb, spatial_hash_dim, atom_count);

    glBindBuffer(GL_UNIFORM_BUFFER, ubo);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo);

    bin_spheres(atom_pos_buffer, atom_count, sphere_aabb);
    compute_sphere_offsets();
    probe_work();
    write_sorted_spheres(atom_pos_buffer, atom_rad_buffer, atom_count, sphere_aabb, max_sphere_radius);

    #if 0
    {
        const uint32_t num_cells = spatial_hash.cell_count;
        uint32_t* cell_offset_data = (uint32_t*)MALLOC(num_cells * sizeof(uint32_t));
        uint32_t* cell_count_data = (uint32_t*)MALLOC(num_cells * sizeof(uint32_t));
        vec3* pos_data = (vec3*)MALLOC(atom_count * sizeof(vec3));
        vec4* sphere_data = (vec4*)MALLOC(atom_count * sizeof(vec4));
        defer {
            FREE(cell_offset_data);
            FREE(cell_count_data);
            FREE(pos_data);
            FREE(sphere_data);
        };

        glBindBuffer(GL_ARRAY_BUFFER, spatial_hash.buffer_cell_offset);
        glGetBufferSubData(GL_ARRAY_BUFFER, 0, num_cells * sizeof(uint32_t), cell_offset_data);

        glBindBuffer(GL_ARRAY_BUFFER, spatial_hash.buffer_cell_count);
        glGetBufferSubData(GL_ARRAY_BUFFER, 0, num_cells * sizeof(uint32_t), cell_count_data);

        glBindBuffer(GL_ARRAY_BUFFER, atom_pos_buffer);
        glGetBufferSubData(GL_ARRAY_BUFFER, 0, atom_count * sizeof(vec3), pos_data);

        glBindBuffer(GL_ARRAY_BUFFER, spatial_hash.buffer_sphere);
        glGetBufferSubData(GL_ARRAY_BUFFER, 0, atom_count * sizeof(vec4), sphere_data);
        printf("DONE!");
    }
    #endif

    write_sphere_distances(0, max_sphere_radius, max_distance);
}

void draw_sdf_debug(const ViewParam& view_param) {
    immediate::set_model_view_matrix(view_param.matrix.current.view);
    immediate::set_proj_matrix(view_param.matrix.current.proj_jittered);
    immediate::draw_box_wireframe(spatial_hash.aabb.min, spatial_hash.aabb.max, immediate::COLOR_CYAN);
    ivec3 dim = spatial_hash.dim;
    vec3 ext = spatial_hash.aabb.ext() / vec3(dim);
    for (int z = 0; z < dim.z; z++) {
        for (int y = 0; y < dim.y; y++) {
            for (int x = 0; x < dim.x; x++) {
                const vec3 min_aabb = spatial_hash.aabb.min + ext * vec3(x, y, z);
                const vec3 max_aabb = min_aabb + ext;
                immediate::draw_box_wireframe(min_aabb, max_aabb, immediate::COLOR_BLACK);
            }
        }
    }
    immediate::flush();
}

void draw_sdf(const ViewParam& view_param) {
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, cube_vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_UNSIGNED_BYTE, GL_FALSE, 0, (const void*)0);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, sdf_volume.texture);

    glUseProgram(program_draw);

    glUniformMatrix4fv(glGetUniformLocation(program_draw, "u_model_view_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.view_proj_jittered[0][0]);
    glUniform4fv(glGetUniformLocation(program_draw, "u_model_eye_pos"), 1, &view_param.matrix.inverse.view[3][0]);
    glUniform1i(glGetUniformLocation(program_draw, "u_tex_sdf"), 0);

    glDrawArrays(GL_TRIANGLE_STRIP, 0, 42);
    glUseProgram(0);

    glBindTexture(GL_TEXTURE_3D, 0);
    glBindVertexArray(0);
}

}  // namespace sdf
}  // namespace draw

namespace draw::scan {

static const size_t GROUPSIZE = 512;
static const size_t BATCH_ELEMENTS = GROUPSIZE * 4;

inline static GLuint snapdiv(GLuint input, GLuint align) { return (input + align - 1) / align; }

static void scan(const GLuint input, const GLuint output, GLuint elements) {
    PUSH_GPU_SECTION("Scan")
    ASSERT((elements % 4) == 0);
    ASSERT(elements < (GLuint64)BATCH_ELEMENTS * BATCH_ELEMENTS);

    glUseProgram(program.prefixsum);
    glUniform1ui(0, elements);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, input);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, output);

    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

    GLuint groups = snapdiv(elements, BATCH_ELEMENTS);
    ASSERT(groups <= max_work_group_count[0]);
    glDispatchCompute(groups, 1, 1);

    if (groups > 1) {
        GLuint groupcombines = snapdiv(groups, BATCH_ELEMENTS);
        ASSERT(groupcombines <= BATCH_ELEMENTS);

        glUseProgram(program.offsets);
        glUniform1ui(0, elements);

        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, input);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, output);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, offset_buffer);

        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        ASSERT(groupcombines <= max_work_group_count[0]);
        glDispatchCompute(groupcombines, 1, 1);

        if (groupcombines > 1) {
            ASSERT(false);
            // @TODO: Compute offsets for groups tier 2,
            // This is only required when elements exceed 2048 * 2048 = ‭4 194 304‬
        }

        glUseProgram(program.combine);
        glUniform1ui(0, elements);

        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, offset_buffer);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, output);

        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        GLuint groups = snapdiv(elements, GROUPSIZE);
        ASSERT(groups < max_work_group_count[0]);
        glDispatchCompute(groups, 1, 1);
    }

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 1);
    POP_GPU_SECTION()
}

void test_scan() {
    const size_t count = 10000;

    GLuint buffers[2];
    glGenBuffers(2, buffers);

    glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
    glBufferData(GL_ARRAY_BUFFER, count * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);

    void* ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    if (ptr) {
        uint32_t* data = (uint32_t*)ptr;
        for (size_t i = 0; i < count; i++) {
            data[i] = i % 10 == 0 ? 1 : 0;
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }

    glBindBuffer(GL_ARRAY_BUFFER, buffers[1]);
    glBufferData(GL_ARRAY_BUFFER, count * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    GLuint query_id[2];
    glGenQueries(2, query_id);

    glQueryCounter(query_id[0], GL_TIMESTAMP);
    scan(buffers[0], buffers[1], count);
    glQueryCounter(query_id[1], GL_TIMESTAMP);

    unsigned int stopTimerAvailable = 0;
    while (!stopTimerAvailable) {
        glGetQueryObjectuiv(query_id[1], GL_QUERY_RESULT_AVAILABLE, &stopTimerAvailable);
    }

    // get query results
    uint64_t t0, t1;
    glGetQueryObjectui64v(query_id[0], GL_QUERY_RESULT, &t0);
    glGetQueryObjectui64v(query_id[1], GL_QUERY_RESULT, &t1);

    LOG_NOTE("Time taken to do prefixsum: %f ms\n", (t1 - t0) / 1000000.0);

    glBindBuffer(GL_ARRAY_BUFFER, buffers[1]);
    ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);
    if (ptr) {
        uint32_t* data = (uint32_t*)ptr;
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }

    glDeleteBuffers(2, buffers);
    glDeleteQueries(2, query_id);
}
}  // namespace draw::scan