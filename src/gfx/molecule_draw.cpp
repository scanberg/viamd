#include "image.h"
#include "molecule_draw.h"

#include <core/common.h>
#include <core/log.h>
#include <core/string_utils.h>
#include <gfx/gl_utils.h>
#include <mol/molecule_utils.h>

static bool is_orthographic_proj_matrix(const mat4& M) { return math::length2(vec2(M[3])) > 0.0f; }

namespace draw {
static GLuint vao = 0;

namespace vdw {
static GLuint program_persp = 0;
static GLuint program_ortho = 0;
static GLuint tex_position = 0;
static GLuint tex_radius = 0;
static GLuint tex_color = 0;
static GLuint tex_view_vel = 0;

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

    if (!tex_position) glGenTextures(4, &tex_position);
}

static void shutdown() {
    if (program_persp) glDeleteProgram(program_persp);
    if (program_ortho) glDeleteProgram(program_ortho);
    if (tex_position) glDeleteTextures(4, &tex_position);
}

}  // namespace vdw

namespace licorice {
static GLuint program = 0;

static GLint uniform_loc_view_mat = -1;
static GLint uniform_loc_proj_mat = -1;
static GLint uniform_loc_curr_view_to_prev_clip_mat = -1;
static GLint uniform_loc_radius = -1;
static GLint uniform_loc_jitter_uv = -1;

static void initialize() {
    GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/licorice.vert", GL_VERTEX_SHADER);
    GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/licorice.geom", GL_GEOMETRY_SHADER);
    GLuint f_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/licorice.frag", GL_FRAGMENT_SHADER);
    defer {
        glDeleteShader(v_shader);
        glDeleteShader(g_shader);
        glDeleteShader(f_shader);
    };

    if (!program) program = glCreateProgram();
    const GLuint shaders[] = {v_shader, g_shader, f_shader};
    gl::attach_link_detach(program, shaders);

    uniform_loc_view_mat = glGetUniformLocation(program, "u_view_mat");
    uniform_loc_proj_mat = glGetUniformLocation(program, "u_proj_mat");
    uniform_loc_curr_view_to_prev_clip_mat = glGetUniformLocation(program, "u_curr_view_to_prev_clip_mat");
    uniform_loc_radius = glGetUniformLocation(program, "u_radius");
    uniform_loc_jitter_uv = glGetUniformLocation(program, "u_jitter_uv");
}

static void shutdown() {
    if (program) glDeleteProgram(program);
}

}  // namespace licorice

namespace ribbon {
static GLuint program = 0;
static GLuint atom_color_tex = 0;
static GLuint atom_velocity_tex = 0;

static GLint uniform_loc_atom_color_buffer = -1;
static GLint uniform_loc_atom_velocity_buffer = -1;
static GLint uniform_loc_normal_mat = -1;
static GLint uniform_loc_view_proj_mat = -1;
static GLint uniform_loc_prev_view_proj_mat = -1;
static GLint uniform_loc_scale = -1;
static GLint uniform_loc_jitter_uv = -1;

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

    uniform_loc_atom_color_buffer = glGetUniformLocation(program, "u_atom_color_buffer");
    uniform_loc_atom_velocity_buffer = glGetUniformLocation(program, "u_atom_velocity_buffer");
    uniform_loc_normal_mat = glGetUniformLocation(program, "u_normal_mat");
    uniform_loc_view_proj_mat = glGetUniformLocation(program, "u_view_proj_mat");
    uniform_loc_prev_view_proj_mat = glGetUniformLocation(program, "u_prev_view_proj_mat");
    uniform_loc_scale = glGetUniformLocation(program, "u_scale");
    uniform_loc_jitter_uv = glGetUniformLocation(program, "u_jitter_uv");

    if (!atom_color_tex) {
        glGenTextures(1, &atom_color_tex);
    }

    if (!atom_velocity_tex) {
        glGenTextures(1, &atom_velocity_tex);
    }
}

void shutdown() {
    if (program) glDeleteProgram(program);
    if (atom_color_tex) glDeleteTextures(1, &atom_color_tex);
    if (atom_velocity_tex) glDeleteTextures(1, &atom_velocity_tex);
}

}  // namespace ribbon

namespace cartoon {
static GLuint program = 0;
static GLuint atom_color_tex = 0;
static GLuint atom_velocity_tex = 0;

static GLint uniform_loc_atom_color_buffer = -1;
static GLint uniform_loc_atom_velocity_buffer = -1;
static GLint uniform_loc_normal_mat = -1;
static GLint uniform_loc_view_proj_mat = -1;
static GLint uniform_loc_prev_view_proj_mat = -1;
static GLint uniform_loc_scale = -1;
static GLint uniform_loc_atom_color_tex = -1;

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

    uniform_loc_atom_color_buffer = glGetUniformLocation(program, "u_atom_color_buffer");
    uniform_loc_atom_velocity_buffer = glGetUniformLocation(program, "u_atom_velocity_buffer");
    uniform_loc_normal_mat = glGetUniformLocation(program, "u_normal_mat");
    uniform_loc_view_proj_mat = glGetUniformLocation(program, "u_view_proj_mat");
    uniform_loc_prev_view_proj_mat = glGetUniformLocation(program, "u_prev_view_proj_mat");
    uniform_loc_scale = glGetUniformLocation(program, "u_scale");
    uniform_loc_atom_color_tex = glGetUniformLocation(program, "u_atom_color_tex");

    if (!atom_color_tex) {
        glGenTextures(1, &atom_color_tex);
    }

    if (!atom_velocity_tex) {
        glGenTextures(1, &atom_velocity_tex);
    }
}

void shutdown() {
    if (program) glDeleteProgram(program);
    if (atom_color_tex) glDeleteTextures(1, &atom_color_tex);
    if (atom_velocity_tex) glDeleteTextures(1, &atom_velocity_tex);
}

}  // namespace cartoon

namespace backbone_spline {
static GLuint extract_control_points_program = 0;
static GLuint compute_spline_program = 0;
static GLuint draw_spline_program = 0;

static GLint uniform_loc_ramachandran_tex = 0;

static GLint uniform_loc_tension = 0;

static GLint uniform_loc_view_proj_mat = 0;
static GLint uniform_loc_s_color = 0;
static GLint uniform_loc_v_color = 0;
static GLint uniform_loc_t_color = 0;

void initialize() {
    {
        GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/backbone_control_points.vert", GL_VERTEX_SHADER);
        GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/backbone_control_points.geom", GL_GEOMETRY_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(g_shader);
        };

        const GLuint shaders[] = {v_shader, g_shader};
        const GLchar* feedback_varyings[] = {"out_control_point", "out_support_vector_xy", "out_support_vector_z_tangent_vector_x", "out_tangent_vector_yz", "out_classification", "out_atom_index"};
        if (!extract_control_points_program) extract_control_points_program = glCreateProgram();
        gl::attach_link_detach_with_transform_feedback(extract_control_points_program, shaders, feedback_varyings, GL_INTERLEAVED_ATTRIBS);

        uniform_loc_ramachandran_tex = glGetUniformLocation(extract_control_points_program, "u_ramachandran_tex");
    }

    {
        GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/backbone_compute_spline.vert", GL_VERTEX_SHADER);
        GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/backbone_compute_spline.geom", GL_GEOMETRY_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(g_shader);
        };

        const GLuint shaders[] = {v_shader, g_shader};
        const GLchar* feedback_varyings[] = {"out_control_point", "out_support_vector_xy", "out_support_vector_z_tangent_vector_x", "out_tangent_vector_yz", "out_classification", "out_atom_index"};
        if (!compute_spline_program) compute_spline_program = glCreateProgram();
        gl::attach_link_detach_with_transform_feedback(compute_spline_program, shaders, feedback_varyings, GL_INTERLEAVED_ATTRIBS);

        // uniform_loc_num_subdivisions = glGetUniformLocation(compute_spline_program, "u_num_subdivisions");
        uniform_loc_tension = glGetUniformLocation(compute_spline_program, "u_tension");
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

        uniform_loc_view_proj_mat = glGetUniformLocation(draw_spline_program, "u_view_proj_mat");
        uniform_loc_s_color = glGetUniformLocation(draw_spline_program, "u_s_color");
        uniform_loc_v_color = glGetUniformLocation(draw_spline_program, "u_v_color");
        uniform_loc_t_color = glGetUniformLocation(draw_spline_program, "u_t_color");
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
static GLuint program = 0;

static GLint uniform_loc_view_mat = -1;
static GLint uniform_loc_proj_mat = -1;
static GLint uniform_loc_inv_proj_mat = -1;
static GLint uniform_loc_jitter_uv = -1;
static GLint uniform_loc_radius_scale = -1;
static GLint uniform_loc_color = -1;
static GLint uniform_loc_mask = -1;

static void initialize() {
    GLuint v_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/lean_and_mean/vdw.vert", GL_VERTEX_SHADER);
    GLuint g_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/lean_and_mean/vdw.geom", GL_GEOMETRY_SHADER);
    GLuint f_shader = gl::compile_shader_from_file(VIAMD_SHADER_DIR "/lean_and_mean/vdw.frag", GL_FRAGMENT_SHADER);
    defer {
        glDeleteShader(v_shader);
        glDeleteShader(g_shader);
        glDeleteShader(f_shader);
    };

    if (!program) program = glCreateProgram();
    const GLuint shaders[] = {v_shader, g_shader, f_shader};
    gl::attach_link_detach(program, shaders);

    uniform_loc_view_mat = glGetUniformLocation(program, "u_view_mat");
    uniform_loc_proj_mat = glGetUniformLocation(program, "u_proj_mat");
    uniform_loc_inv_proj_mat = glGetUniformLocation(program, "u_inv_proj_mat");
    uniform_loc_jitter_uv = glGetUniformLocation(program, "u_jitter_uv");
    uniform_loc_radius_scale = glGetUniformLocation(program, "u_radius_scale");
    uniform_loc_color = glGetUniformLocation(program, "u_color");
    uniform_loc_mask = glGetUniformLocation(program, "u_mask");
}

static void shutdown() {
    if (program) glDeleteProgram(program);
}
}  // namespace vdw

namespace licorice {
static GLuint program = 0;

static GLint uniform_loc_view_mat = -1;
static GLint uniform_loc_proj_mat = -1;
static GLint uniform_loc_radius = -1;
static GLint uniform_loc_mask = -1;
static GLint uniform_loc_color = -1;

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

    uniform_loc_view_mat = glGetUniformLocation(program, "u_view_mat");
    uniform_loc_proj_mat = glGetUniformLocation(program, "u_proj_mat");
    uniform_loc_radius = glGetUniformLocation(program, "u_radius");
    uniform_loc_mask = glGetUniformLocation(program, "u_mask");
    uniform_loc_color = glGetUniformLocation(program, "u_color");
}

static void shutdown() {
    if (program) glDeleteProgram(program);
}

}  // namespace licorice

namespace ribbon {
static GLuint program = 0;
static GLuint atom_color_tex = 0;
static GLuint atom_mask_tex = 0;

static GLint uniform_loc_atom_color_buffer = -1;
static GLint uniform_loc_atom_mask_buffer = -1;
static GLint uniform_loc_view_proj_mat = -1;
static GLint uniform_loc_scale = -1;
static GLint uniform_loc_color = -1;
static GLint uniform_loc_mask = -1;

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

    uniform_loc_atom_color_buffer = glGetUniformLocation(program, "u_atom_color_buffer");
    uniform_loc_atom_mask_buffer = glGetUniformLocation(program, "u_atom_mask_buffer");
    uniform_loc_view_proj_mat = glGetUniformLocation(program, "u_view_proj_mat");
    uniform_loc_scale = glGetUniformLocation(program, "u_scale");
    uniform_loc_color = glGetUniformLocation(program, "u_color");
    uniform_loc_mask = glGetUniformLocation(program, "u_mask");

    if (!atom_color_tex) {
        glGenTextures(1, &atom_color_tex);
    }

    if (!atom_mask_tex) {
        glGenTextures(1, &atom_mask_tex);
    }
}

void shutdown() {
    if (program) glDeleteProgram(program);
    if (atom_color_tex) glDeleteTextures(1, &atom_color_tex);
    if (atom_mask_tex) glDeleteTextures(1, &atom_mask_tex);
}

}  // namespace ribbon

}  // namespace lean_and_mean

void initialize() {
    if (!vao) glGenVertexArrays(1, &vao);
    vdw::initialize();
    licorice::initialize();
    ribbon::intitialize();
    cartoon::intitialize();
    backbone_spline::initialize();
    lean_and_mean::vdw::initialize();
    lean_and_mean::licorice::initialize();
    lean_and_mean::ribbon::intitialize();
    pbc_view_velocity::initialize();
}

void shutdown() {
    if (vao) glDeleteVertexArrays(1, &vao);
    vdw::shutdown();
    licorice::shutdown();
    ribbon::shutdown();
    cartoon::shutdown();
    backbone_spline::shutdown();
    lean_and_mean::vdw::shutdown();
    lean_and_mean::licorice::shutdown();
    lean_and_mean::ribbon::shutdown();
    pbc_view_velocity::shutdown();
}

void draw_vdw(GLuint atom_position_buffer, GLuint atom_radius_buffer, GLuint atom_color_buffer, GLuint atom_view_velocity_buffer, int32 atom_count, const ViewParam& view_param, float radius_scale) {
    ASSERT(glIsBuffer(atom_position_buffer));
    ASSERT(glIsBuffer(atom_radius_buffer));
    ASSERT(glIsBuffer(atom_color_buffer));
    ASSERT(glIsBuffer(atom_view_velocity_buffer));

    const mat4 curr_view_to_prev_clip_mat = view_param.matrix.previous.view_proj_jittered * view_param.matrix.inverse.view;
    const vec2 res = view_param.resolution;
    const vec4 jitter_uv = vec4(view_param.jitter.current / res, view_param.jitter.previous / res);
    static uint32 frame = 0;
    frame++;

    const bool ortho = is_orthographic_proj_matrix(view_param.matrix.current.proj);

    if (ortho) {
        glBindTexture(GL_TEXTURE_BUFFER, vdw::tex_position);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, atom_position_buffer);
        glUniform1i(glGetUniformLocation(vdw::program_ortho, "u_buf_position"), 0);

        glBindTexture(GL_TEXTURE_BUFFER, vdw::tex_radius);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, atom_radius_buffer);
        glUniform1i(glGetUniformLocation(vdw::program_ortho, "u_buf_radius"), 1);

        glBindTexture(GL_TEXTURE_BUFFER, vdw::tex_color);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color_buffer);
        glUniform1i(glGetUniformLocation(vdw::program_ortho, "u_buf_color"), 2);

        glBindTexture(GL_TEXTURE_BUFFER, vdw::tex_view_vel);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, atom_view_velocity_buffer);
        glUniform1i(glGetUniformLocation(vdw::program_ortho, "u_buf_view_velocity"), 3);

        glBindBuffer(GL_TEXTURE_BUFFER, 0);
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

        glUseProgram(vdw::program_persp);

        // Uniforms
        glUniformMatrix4fv(glGetUniformLocation(vdw::program_persp, "u_view_mat"), 1, GL_FALSE, &view_param.matrix.current.view[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(vdw::program_persp, "u_proj_mat"), 1, GL_FALSE, &view_param.matrix.current.proj_jittered[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(vdw::program_persp, "u_inv_proj_mat"), 1, GL_FALSE, &view_param.matrix.inverse.proj_jittered[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(vdw::program_persp, "u_curr_view_to_prev_clip_mat"), 1, GL_FALSE, &curr_view_to_prev_clip_mat[0][0]);
        glUniform1f(glGetUniformLocation(vdw::program_persp, "u_radius_scale"), radius_scale);
        glUniform4fv(glGetUniformLocation(vdw::program_persp, "u_jitter_uv"), 1, &jitter_uv[0]);
        glUniform1ui(glGetUniformLocation(vdw::program_persp, "u_frame"), frame);

        glDrawArrays(GL_POINTS, 0, atom_count);
        glBindVertexArray(0);
    }
    glUseProgram(0);
}

void draw_licorice(GLuint atom_position_buffer, GLuint atom_color_buffer, GLuint atom_velocity_buffer, GLuint bond_buffer, int32 bond_count, const ViewParam& view_param, float radius_scale) {
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

    glUseProgram(licorice::program);
    glUniformMatrix4fv(licorice::uniform_loc_view_mat, 1, GL_FALSE, &view_param.matrix.current.view[0][0]);
    glUniformMatrix4fv(licorice::uniform_loc_proj_mat, 1, GL_FALSE, &view_param.matrix.current.proj_jittered[0][0]);
    glUniformMatrix4fv(licorice::uniform_loc_curr_view_to_prev_clip_mat, 1, GL_FALSE, &curr_view_to_prev_clip_mat[0][0]);
    glUniform1f(licorice::uniform_loc_radius, 0.25f * radius_scale);
    glUniform4fv(licorice::uniform_loc_jitter_uv, 1, &jitter_uv[0]);

    glDrawElements(GL_LINES, bond_count * 2, GL_UNSIGNED_INT, 0);
    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void compute_backbone_control_points(GLuint dst_buffer, GLuint atom_position_buffer, GLuint backbone_index_buffer, int num_backbone_indices, GLuint ramachandran_tex) {
    glEnable(GL_RASTERIZER_DISCARD);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, atom_position_buffer);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(AtomPosition), (const GLvoid*)0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, backbone_index_buffer);

    glUseProgram(backbone_spline::extract_control_points_program);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, ramachandran_tex);
    glUniform1i(backbone_spline::uniform_loc_ramachandran_tex, 0);

    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, dst_buffer);
    glBeginTransformFeedback(GL_POINTS);
    glDrawElements(GL_TRIANGLES_ADJACENCY, num_backbone_indices, GL_UNSIGNED_INT, 0);
    glEndTransformFeedback();

    glUseProgram(0);

    glDisable(GL_RASTERIZER_DISCARD);
}

void compute_backbone_spline(GLuint dst_buffer, GLuint control_point_buffer, GLuint control_point_index_buffer, int num_control_point_indices, float tension) {
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

    // glUniform1i(backbone_spline::uniform_loc_num_subdivisions, num_subdivisions);
    glUniform1f(backbone_spline::uniform_loc_tension, tension);

    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, dst_buffer);
    glBeginTransformFeedback(GL_POINTS);
    glDrawElements(GL_LINE_STRIP_ADJACENCY, num_control_point_indices, GL_UNSIGNED_INT, 0);
    glEndTransformFeedback();

    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);
    glDisable(GL_RASTERIZER_DISCARD);
}

void compute_pbc_view_velocity(GLuint dst_buffer, GLuint position_buffer, GLuint old_position_buffer, int count, const ViewParam& view_param, const vec3& box_ext) {
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

void draw_spline(GLuint spline_buffer, GLuint spline_index_buffer, int32 num_spline_indices, const ViewParam& view_param, uint32 s_color, uint32 v_color, uint32 t_color) {
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

    glUseProgram(backbone_spline::draw_spline_program);
    glUniformMatrix4fv(backbone_spline::uniform_loc_view_proj_mat, 1, GL_FALSE, &(view_param.matrix.current.view_proj_jittered)[0][0]);
    glUniform4fv(backbone_spline::uniform_loc_s_color, 1, &math::convert_color(s_color)[0]);
    glUniform4fv(backbone_spline::uniform_loc_v_color, 1, &math::convert_color(v_color)[0]);
    glUniform4fv(backbone_spline::uniform_loc_t_color, 1, &math::convert_color(t_color)[0]);
    glDrawElements(GL_LINE_STRIP, num_spline_indices, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);
}

void draw_ribbons(GLuint spline_buffer, GLuint spline_index_buffer, GLuint atom_color_buffer, GLuint atom_velocity_buffer, int32 num_spline_indices, const ViewParam& view_param) {
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
    glBindTexture(GL_TEXTURE_BUFFER, ribbon::atom_color_tex);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color_buffer);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ribbon::atom_velocity_tex);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, atom_velocity_buffer);

    const vec2 res = view_param.resolution;
    const vec4 jitter_uv = vec4(view_param.jitter.current / res, view_param.jitter.previous / res);

    glUseProgram(ribbon::program);

    // Uniforms
    glUniform1i(ribbon::uniform_loc_atom_color_buffer, 0);
    glUniform1i(ribbon::uniform_loc_atom_velocity_buffer, 1);
    glUniformMatrix4fv(ribbon::uniform_loc_normal_mat, 1, GL_FALSE, &view_param.matrix.current.norm[0][0]);
    glUniformMatrix4fv(ribbon::uniform_loc_view_proj_mat, 1, GL_FALSE, &view_param.matrix.current.view_proj_jittered[0][0]);
    glUniformMatrix4fv(ribbon::uniform_loc_prev_view_proj_mat, 1, GL_FALSE, &view_param.matrix.previous.view_proj_jittered[0][0]);
    glUniform4fv(ribbon::uniform_loc_jitter_uv, 1, &jitter_uv[0]);

    glDrawElements(GL_LINE_STRIP, num_spline_indices, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);
}

void draw_cartoon(GLuint spline_buffer, GLuint spline_index_buffer, GLuint atom_color_buffer, int32 num_spline_indices, const ViewParam& view_param) {
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
    glBindTexture(GL_TEXTURE_BUFFER, cartoon::atom_color_tex);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color_buffer);

    glUseProgram(cartoon::program);
    glUniform1i(cartoon::uniform_loc_atom_color_tex, 0);
    glUniformMatrix4fv(cartoon::uniform_loc_normal_mat, 1, GL_FALSE, &view_param.matrix.current.norm[0][0]);
    glUniformMatrix4fv(cartoon::uniform_loc_view_proj_mat, 1, GL_FALSE, &view_param.matrix.current.view_proj_jittered[0][0]);
    glDrawElements(GL_LINE_STRIP, num_spline_indices, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);
}

namespace lean_and_mean {

void draw_vdw(GLuint atom_position_buffer, GLuint atom_radius_buffer, GLuint atom_color_buffer, GLuint atom_mask_buffer, int32 atom_count, const ViewParam& view_param, float radius_scale, vec4 color,
              uint32 mask) {
    ASSERT(glIsBuffer(atom_position_buffer));
    ASSERT(glIsBuffer(atom_radius_buffer));
    ASSERT(glIsBuffer(atom_color_buffer));
    ASSERT(glIsBuffer(atom_mask_buffer));

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

    const vec2 res = view_param.resolution;
    const vec4 jitter_uv = vec4(view_param.jitter.current / res, view_param.jitter.previous / res);

    glUseProgram(vdw::program);

    // Uniforms
    glUniformMatrix4fv(vdw::uniform_loc_view_mat, 1, GL_FALSE, &view_param.matrix.current.view[0][0]);
    glUniformMatrix4fv(vdw::uniform_loc_proj_mat, 1, GL_FALSE, &view_param.matrix.current.proj_jittered[0][0]);
    glUniformMatrix4fv(vdw::uniform_loc_inv_proj_mat, 1, GL_FALSE, &view_param.matrix.inverse.proj_jittered[0][0]);
    glUniform4fv(vdw::uniform_loc_jitter_uv, 1, &jitter_uv[0]);
    glUniform4fv(vdw::uniform_loc_color, 1, &color[0]);
    glUniform1f(vdw::uniform_loc_radius_scale, radius_scale);
    glUniform1ui(vdw::uniform_loc_mask, mask);

    glDrawArrays(GL_POINTS, 0, atom_count);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    glUseProgram(0);
    glBindVertexArray(0);
}

void draw_licorice(GLuint atom_position_buffer, GLuint atom_color_buffer, GLuint atom_mask_buffer, GLuint bond_buffer, int32 bond_count, const ViewParam& view_param, float radius_scale, vec4 color,
                   uint32 mask) {
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

    glUseProgram(licorice::program);
    glUniformMatrix4fv(licorice::uniform_loc_view_mat, 1, GL_FALSE, &view_param.matrix.current.view[0][0]);
    glUniformMatrix4fv(licorice::uniform_loc_proj_mat, 1, GL_FALSE, &view_param.matrix.current.proj_jittered[0][0]);
    glUniform4fv(licorice::uniform_loc_color, 1, &color[0]);
    glUniform1ui(licorice::uniform_loc_mask, mask);
    glUniform1f(licorice::uniform_loc_radius, 0.25f * radius_scale);

    glDrawElements(GL_LINES, bond_count * 2, GL_UNSIGNED_INT, 0);
    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void draw_ribbons(GLuint spline_buffer, GLuint spline_index_buffer, GLuint atom_color_buffer, GLuint atom_mask_buffer, int32 num_spline_indices, const ViewParam& view_param, float scale, vec4 color,
                  uint32 mask) {
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
    glBindTexture(GL_TEXTURE_BUFFER, ribbon::atom_color_tex);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color_buffer);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ribbon::atom_mask_tex);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, atom_mask_buffer);

    glUseProgram(ribbon::program);

    const vec2 ribbon_scale = vec2(1.0f, 0.1f) * scale;

    // Uniforms
    glUniformMatrix4fv(ribbon::uniform_loc_view_proj_mat, 1, GL_FALSE, &view_param.matrix.current.view_proj_jittered[0][0]);
    glUniform4fv(ribbon::uniform_loc_color, 1, &color[0]);
    glUniform2fv(ribbon::uniform_loc_scale, 1, &ribbon_scale[0]);
    glUniform1ui(ribbon::uniform_loc_mask, mask);
    glUniform1i(ribbon::uniform_loc_atom_color_buffer, 0);
    glUniform1i(ribbon::uniform_loc_atom_mask_buffer, 1);

    glDrawElements(GL_LINE_STRIP, num_spline_indices, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);
}

}  // namespace lean_and_mean

}  // namespace draw
