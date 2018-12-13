#include "immediate_draw_utils.h"
#include <core/types.h>
#include <core/array_types.h>
#include <core/gl.h>
#include <core/log.h>
#include <core/math_utils.h>
#include <gfx/gl_utils.h>

namespace immediate {

struct Vertex {
    vec3 position = {0, 0, 0};
    vec3 normal = {0, 0, 1};
    vec2 uv = {0, 0};
    uint32 color = DEFAULT_COLOR;
};

using Index = uint16;

struct DrawCommand {
    Index offset;
    Index count;
    GLenum primitive_type;
    GLuint program;
    int view_matrix_idx = -1;
    int proj_matrix_idx = -1;
    int material_idx = -1;
};

static DynamicArray<mat4> matrix_stack;
static DynamicArray<Material> material_stack;

static DynamicArray<DrawCommand> commands;
static DynamicArray<Vertex> vertices;
static DynamicArray<Index> indices;

static GLuint vbo = 0;
static GLuint ibo = 0;
static GLuint vao = 0;
static GLuint ubo_material = 0;
static GLuint default_tex = 0;

static GLuint program = 0;

static GLint uniform_loc_mvp_matrix = -1;
static GLint uniform_loc_normal_matrix = -1;
static GLint uniform_loc_point_size = -1;
static GLint uniform_block_index_material = -1;

static int curr_view_matrix_idx = -1;
static int curr_proj_matrix_idx = -1;
static int curr_material_idx = -1;

static const char* v_shader_src = R"(
#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

struct Material {
	vec3 f0;
	float smoothness;
	vec2 uv_scale;
	float _pad_0;
	float _pad_1;
};

layout(std140) uniform u_material_buffer {
	Material material;
};

uniform mat4 u_mvp_matrix;
uniform mat3 u_normal_matrix;
uniform float u_point_size = 1.f;

layout(location = 0) in vec3 in_position;
layout(location = 1) in vec3 in_normal;
layout(location = 2) in vec2 in_uv;
layout(location = 3) in vec4 in_color;

out vec3 normal;
out vec2 uv;
out vec4 color;

void main() {
	gl_Position = u_mvp_matrix * vec4(in_position, 1);
    gl_PointSize = max(u_point_size, 400.f / gl_Position.w);
	normal = u_normal_matrix * in_normal;
	uv = in_uv;
	color = in_color;
}
)";

static const char* f_shader_src = R"(
#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

struct Material {
	vec3 f0;
	float smoothness;
	vec2 uv_scale;
	float _pad_0;
	float _pad_1;
};

layout(std140) uniform u_material_buffer {
	Material material;
};

uniform sampler2D u_base_color_texture;

in vec3 normal;
in vec2 uv;
in vec4 color;

layout(location = 0) out vec4 out_color_alpha;
layout(location = 1) out vec4 out_f0_smoothness;
layout(location = 2) out vec4 out_normal;

vec4 encode_normal (vec3 n) {
    float p = sqrt(n.z*8+8);
    return vec4(n.xy/p + 0.5,0,0);
}

void main() {
	out_color_alpha = texture(u_base_color_texture, uv) * color;
	out_f0_smoothness = vec4(material.f0, material.smoothness);
	out_normal = encode_normal(normalize(normal));
}
)";

static inline void append_draw_command(Index offset, Index count, GLenum primitive_type) {
    if (commands.count > 0 && commands.back().primitive_type == primitive_type) {
        commands.back().count += count;
    } else {
        ASSERT(curr_view_matrix_idx > -1, "Immediate Mode View Matrix not set!");
        ASSERT(curr_proj_matrix_idx > -1, "Immediate Mode Proj Matrix not set!");
        ASSERT(curr_material_idx > -1, "Material not set!");

        DrawCommand cmd{offset, count, primitive_type, program, curr_view_matrix_idx, curr_proj_matrix_idx, curr_material_idx};
        commands.push_back(cmd);
    }
}

void initialize() {
    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    GLuint v_shader = glCreateShader(GL_VERTEX_SHADER);
    GLuint f_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(v_shader, 1, &v_shader_src, 0);
    glShaderSource(f_shader, 1, &f_shader_src, 0);
    glCompileShader(v_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, v_shader)) {
        LOG_ERROR("Error while compiling immediate vertex shader:\n%s\n", buffer);
    }
    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        LOG_ERROR("Error while compiling immediate fragment shader:\n%s\n", buffer);
    }

    program = glCreateProgram();
    glAttachShader(program, v_shader);
    glAttachShader(program, f_shader);
    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, program)) {
        LOG_ERROR("Error while linking immediate program:\n%s\n", buffer);
    }

    glDetachShader(program, v_shader);
    glDetachShader(program, f_shader);

    glDeleteShader(v_shader);
    glDeleteShader(f_shader);

    uniform_loc_mvp_matrix = glGetUniformLocation(program, "u_mvp_matrix");
    uniform_loc_normal_matrix = glGetUniformLocation(program, "u_normal_matrix");
    uniform_loc_point_size = glGetUniformLocation(program, "u_point_size");
    uniform_block_index_material = glGetUniformBlockIndex(program, "u_material_buffer");

    glGenBuffers(1, &vbo);
    glGenBuffers(1, &ibo);

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)offsetof(Vertex, position));

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)offsetof(Vertex, normal));

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)offsetof(Vertex, uv));

    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (const GLvoid*)offsetof(Vertex, color));

    glBindVertexArray(0);

    if (!ubo_material) glGenBuffers(1, &ubo_material);
    glBindBuffer(GL_UNIFORM_BUFFER, ubo_material);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(Material), nullptr, GL_DYNAMIC_DRAW);

    constexpr uint32 pixel_data = 0xffffffff;
    if (!default_tex) glGenTextures(1, &default_tex);
    glBindTexture(GL_TEXTURE_2D, default_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, &pixel_data);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    vertices.reserve(100000);
    indices.reserve(100000);
}

void shutdown() {
    if (vbo) glDeleteBuffers(1, &vbo);
    if (ibo) glDeleteBuffers(1, &ibo);
    if (vao) glDeleteVertexArrays(1, &vao);
    if (program) glDeleteProgram(program);
}

void set_view_matrix(const mat4& model_view_matrix) {
    curr_view_matrix_idx = (int)matrix_stack.count;
    matrix_stack.push_back(model_view_matrix);
}

void set_proj_matrix(const mat4& proj_matrix) {
    curr_proj_matrix_idx = (int)matrix_stack.count;
    matrix_stack.push_back(proj_matrix);
}

void set_material(const Material& material) {
    curr_material_idx = (int)material_stack.count;
    material_stack.push_back(material);
}

void flush() {
    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, vertices.count * sizeof(Vertex), vertices.data, GL_STREAM_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.count * sizeof(Index), indices.data, GL_STREAM_DRAW);

    glEnable(GL_PROGRAM_POINT_SIZE);
    // glEnable(GL_BLEND);
    // glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glUseProgram(program);
    // glPointSize(20.f);
    glLineWidth(1.f);

    glUniform1f(uniform_loc_point_size, 4.f);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_material);
    glUniformBlockBinding(program, uniform_block_index_material, 0);

    int current_view_matrix_idx = -999;
    int current_proj_matrix_idx = -999;
    int current_material_idx = -999;

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, default_tex);

    for (const auto& cmd : commands) {
        bool update_view = false;
        bool update_mvp = false;
        if (cmd.view_matrix_idx != current_view_matrix_idx) {
            current_view_matrix_idx = cmd.view_matrix_idx;
            update_view = true;
            update_mvp = true;
        }
        if (cmd.proj_matrix_idx != current_proj_matrix_idx) {
            current_proj_matrix_idx = cmd.proj_matrix_idx;
            update_mvp = true;
        }

        if (update_view) {
            mat3 normal_matrix = mat3(math::transpose(math::inverse(matrix_stack[cmd.view_matrix_idx])));
            glUniformMatrix3fv(uniform_loc_normal_matrix, 1, GL_FALSE, &normal_matrix[0][0]);
        }
        if (update_mvp) {
            mat4 mvp_matrix = matrix_stack[cmd.proj_matrix_idx] * matrix_stack[cmd.view_matrix_idx];
            glUniformMatrix4fv(uniform_loc_mvp_matrix, 1, GL_FALSE, &mvp_matrix[0][0]);
        }
        if (cmd.material_idx != current_material_idx) {
            current_material_idx = cmd.material_idx;
            const Material& material = cmd.material_idx == -1 ? DEFAULT_MATERIAL : material_stack[cmd.material_idx];
            glBindBuffer(GL_UNIFORM_BUFFER, ubo_material);
            glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(Material), &material);
            glBindBuffer(GL_UNIFORM_BUFFER, 0);
            if (material.texture_id != 0) {
                glBindTexture(GL_TEXTURE_2D, material.texture_id);
            } else {
                glBindTexture(GL_TEXTURE_2D, default_tex);
            }
        }

        glDrawElements(cmd.primitive_type, cmd.count, sizeof(Index) == 2 ? GL_UNSIGNED_SHORT : GL_UNSIGNED_INT,
                       reinterpret_cast<const void*>(cmd.offset * sizeof(Index)));
    }

    glBindVertexArray(0);
    glUseProgram(0);

    glDisable(GL_PROGRAM_POINT_SIZE);
    // glDisable(GL_BLEND);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    vertices.clear();
    indices.clear();
    commands.clear();
    matrix_stack.clear();
    material_stack.clear();
    curr_view_matrix_idx = -1;
    curr_proj_matrix_idx = -1;
    curr_material_idx = -1;
}

// PRIMITIVES
void draw_point(const vec3& pos, uint32 color) {
    const Index idx = (Index)vertices.count;

    vertices.push_back({pos, vec3(0, 0, 1), {0, 0}, color});
    indices.push_back(idx);

    append_draw_command(idx, 1, GL_POINTS);
}

void draw_line(const vec3& from, const vec3& to, uint32 color) {
    const Index idx = (Index)vertices.count;

    vertices.push_back({from, {0, 0, 1}, {0, 0}, color});
    vertices.push_back({to, {0, 0, 1}, {0, 1}, color});

    indices.push_back(idx);
    indices.push_back(idx + 1);

    append_draw_command(idx, 2, GL_LINES);
}

void draw_triangle(const vec3& p0, const vec3& p1, const vec3& p2, uint32 color) {
    const Index idx = (Index)vertices.count;
    const vec3 normal = math::normalize(math::cross(p1 - p0, p2 - p0));

    vertices.push_back({p0, normal, {0, 0}, color});
    vertices.push_back({p1, normal, {1, 0}, color});
    vertices.push_back({p2, normal, {0, 1}, color});

    indices.push_back(idx);
    indices.push_back(idx + 1);
    indices.push_back(idx + 2);

    append_draw_command(idx, 3, GL_TRIANGLES);
}

void draw_plane(const vec3& center, const vec3& vec_u, const vec3& vec_v, uint32 color) {
    const Index idx = (Index)vertices.count;
    const vec3 normal = math::normalize(math::cross(vec_u, vec_v));

    vertices.push_back({{center - vec_u + vec_v}, normal, {0, 1}});
    vertices.push_back({{center - vec_u - vec_v}, normal, {0, 0}});
    vertices.push_back({{center + vec_u + vec_v}, normal, {1, 1}});
    vertices.push_back({{center + vec_u - vec_v}, normal, {1, 0}});

    indices.push_back(idx);
    indices.push_back(idx + 1);
    indices.push_back(idx + 2);
    indices.push_back(idx + 2);
    indices.push_back(idx + 1);
    indices.push_back(idx + 3);

    append_draw_command(idx, 6, GL_TRIANGLES);
}

/*
void draw_aabb(const vec3& min_box, const vec3& max_box) {
    const Index idx = (Index)vertices.count;
    const vec3 normal = {0, 0, 1};

    // @ TODO: This is incorrect and needs to be fixed

    vertices.push_back({{max_box[0], max_box[1], max_box[2]}, normal});
    vertices.push_back({{min_box[0], max_box[1], max_box[2]}, normal});
    vertices.push_back({{max_box[0], max_box[1], min_box[2]}, normal});
    vertices.push_back({{min_box[0], max_box[1], min_box[2]}, normal});

    vertices.push_back({{max_box[0], min_box[1], max_box[2]}, normal});
    vertices.push_back({{min_box[0], min_box[1], max_box[2]}, normal});
    vertices.push_back({{max_box[0], min_box[1], min_box[2]}, normal});
    vertices.push_back({{min_box[0], min_box[1], min_box[2]}, normal});

    indices.push_back(idx + 4 - 1);
    indices.push_back(idx + 3 - 1);
    indices.push_back(idx + 7 - 1);

    indices.push_back(idx + 8 - 1);
    indices.push_back(idx + 5 - 1);
    indices.push_back(idx + 3 - 1);

    indices.push_back(idx + 1 - 1);
    indices.push_back(idx + 4 - 1);
    indices.push_back(idx + 2 - 1);

    indices.push_back(idx + 7 - 1);
    indices.push_back(idx + 6 - 1);
    indices.push_back(idx + 5 - 1);

    indices.push_back(idx + 2 - 1);
    indices.push_back(idx + 1 - 1);

    append_draw_command(idx, 14, GL_TRIANGLE_STRIP);
}
*/

void draw_aabb_lines(const vec3& min_box, const vec3& max_box, uint32 color) {
    // Z = min
    draw_line(vec3(min_box[0], min_box[1], min_box[2]), vec3(max_box[0], min_box[1], min_box[2]), color);
    draw_line(vec3(min_box[0], min_box[1], min_box[2]), vec3(min_box[0], max_box[1], min_box[2]), color);
    draw_line(vec3(max_box[0], min_box[1], min_box[2]), vec3(max_box[0], max_box[1], min_box[2]), color);
    draw_line(vec3(min_box[0], max_box[1], min_box[2]), vec3(max_box[0], max_box[1], min_box[2]), color);

    // Z = max
    draw_line(vec3(min_box[0], min_box[1], max_box[2]), vec3(max_box[0], min_box[1], max_box[2]), color);
    draw_line(vec3(min_box[0], min_box[1], max_box[2]), vec3(min_box[0], max_box[1], max_box[2]), color);
    draw_line(vec3(max_box[0], min_box[1], max_box[2]), vec3(max_box[0], max_box[1], max_box[2]), color);
    draw_line(vec3(min_box[0], max_box[1], max_box[2]), vec3(max_box[0], max_box[1], max_box[2]), color);

    // Z min max
    draw_line(vec3(min_box[0], min_box[1], min_box[2]), vec3(min_box[0], min_box[1], max_box[2]), color);
    draw_line(vec3(min_box[0], max_box[1], min_box[2]), vec3(min_box[0], max_box[1], max_box[2]), color);
    draw_line(vec3(max_box[0], min_box[1], min_box[2]), vec3(max_box[0], min_box[1], max_box[2]), color);
    draw_line(vec3(max_box[0], max_box[1], min_box[2]), vec3(max_box[0], max_box[1], max_box[2]), color);
}

void draw_basis(const mat4& basis, const float scale, uint32 x_color, uint32 y_color, uint32 z_color) {
    const vec3 o = vec3(basis[3]);
    const vec3 x = o + vec3(basis[0]) * scale;
    const vec3 y = o + vec3(basis[1]) * scale;
    const vec3 z = o + vec3(basis[2]) * scale;

    draw_line(o, x, x_color);
    draw_line(o, y, y_color);
    draw_line(o, z, z_color);
}

}  // namespace immediate
