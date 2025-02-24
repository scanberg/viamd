#include <gfx/immediate_draw_utils.h>

#include <gfx/gl_utils.h>
#include <color_utils.h>

#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>

#include <stddef.h>

namespace immediate {

using Index = uint32_t;

struct DrawCommand {
    uint32_t offset = 0;
    uint32_t count = 0;
    int32_t view_matrix_idx = -1;
    int32_t proj_matrix_idx = -1;
    GLenum primitive_type;
};

static md_array(mat4_t) matrix_stack;

static md_array(DrawCommand) commands;
static md_array(Vertex) vertices;
static md_array(Index) indices;

static md_allocator_i* arena;

static GLuint vbo = 0;
static GLuint ibo = 0;
static GLuint vao = 0;
static GLuint default_tex = 0;

static GLuint program = 0;

static GLint uniform_loc_mvp_matrix = -1;
static GLint uniform_loc_normal_matrix = -1;
static GLint uniform_loc_uv_scale = -1;
static GLint uniform_loc_point_size = -1;

static int32_t curr_view_matrix_idx = -1;
static int32_t curr_proj_matrix_idx = -1;

static const char* v_shader_src = R"(
#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_mvp_matrix;
//uniform mat3 u_normal_matrix;
//uniform vec2 u_uv_scale = vec2(1,1);
uniform float u_point_size = 1.f;

layout(location = 0) in vec3 in_position;
layout(location = 1) in vec4 in_color;

out vec4 color;

void main() {
	gl_Position = u_mvp_matrix * vec4(in_position, 1);
    gl_PointSize = max(u_point_size, 200.f / gl_Position.w);
	color = in_color;
}
)";

static const char* f_shader_src = R"(
#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

in vec4 color;

layout(location = 0) out vec4 out_color_alpha;
layout(location = 1) out vec4 out_normal;

vec4 encode_normal (vec3 n) {
    float p = sqrt(n.z*8+8);
    return vec4(n.xy/p + 0.5,0,0);
}

void main() {
	out_color_alpha = color;
	out_normal = encode_normal(vec3(0,0,1));
}
)";

static inline void append_draw_command(size_t count, GLenum primitive_type) {
    const uint32_t max_size = (sizeof(Index) == 2 ? 0xFFFFU : 0xFFFFFFFFU);
    ASSERT(md_array_size(indices) + count < max_size);
    // Can we append data to previous draw command?
    DrawCommand* last_cmd = md_array_last(commands);
    if (last_cmd &&
        last_cmd->primitive_type == primitive_type &&
        last_cmd->view_matrix_idx == curr_view_matrix_idx &&
        last_cmd->proj_matrix_idx == curr_proj_matrix_idx) {
        size_t capacity = max_size - last_cmd->count;

        if (count < capacity) {
            last_cmd->count += (uint32_t)count;
            return;
        }
        else {
            last_cmd->count += (uint32_t)capacity;
            count -= capacity;
        }
    }
    ASSERT(curr_view_matrix_idx > -1 && "Immediate Mode View Matrix not set!");
    ASSERT(curr_proj_matrix_idx > -1 && "Immediate Mode Proj Matrix not set!");

    const size_t offset = md_array_size(indices) - count;
    DrawCommand cmd {(uint32_t)offset, (uint32_t)count, curr_view_matrix_idx, curr_proj_matrix_idx, primitive_type};
    
    md_array_push(commands, cmd, arena);
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
        MD_LOG_ERROR("Error while compiling immediate vertex shader:\n%s\n", buffer);
    }
    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        MD_LOG_ERROR("Error while compiling immediate fragment shader:\n%s\n", buffer);
    }

    program = glCreateProgram();
    glAttachShader(program, v_shader);
    glAttachShader(program, f_shader);
    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, program)) {
        MD_LOG_ERROR("Error while linking immediate program:\n%s\n", buffer);
    }

    glDetachShader(program, v_shader);
    glDetachShader(program, f_shader);

    glDeleteShader(v_shader);
    glDeleteShader(f_shader);

    uniform_loc_mvp_matrix = glGetUniformLocation(program, "u_mvp_matrix");
    uniform_loc_normal_matrix = glGetUniformLocation(program, "u_normal_matrix");
    uniform_loc_uv_scale = glGetUniformLocation(program, "u_uv_scale");
    uniform_loc_point_size = glGetUniformLocation(program, "u_point_size");

    glGenBuffers(1, &vbo);
    glGenBuffers(1, &ibo);

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)offsetof(Vertex, coord));

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (const GLvoid*)offsetof(Vertex, color));

    glBindVertexArray(0);

    constexpr uint32_t pixel_data = 0xffffffff;
    if (!default_tex) glGenTextures(1, &default_tex);
    glBindTexture(GL_TEXTURE_2D, default_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, &pixel_data);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (!arena) {
        arena = md_vm_arena_create(GIGABYTES(1));
    }

    md_array_ensure(vertices, 100000, arena);
    md_array_ensure(indices,  100000, arena);
    md_array_ensure(matrix_stack, 10, arena);
}

void shutdown() {
    if (vbo) glDeleteBuffers(1, &vbo);
    if (ibo) glDeleteBuffers(1, &ibo);
    if (vao) glDeleteVertexArrays(1, &vao);
    if (program) glDeleteProgram(program);
    md_vm_arena_destroy(arena);
}

void set_model_view_matrix(const mat4_t& model_view_matrix) {
    curr_view_matrix_idx = (int)md_array_size(matrix_stack);
    md_array_push(matrix_stack, model_view_matrix, arena);
}

void set_proj_matrix(const mat4_t& proj_matrix) {
    curr_proj_matrix_idx = (int)md_array_size(matrix_stack);
    md_array_push(matrix_stack, proj_matrix, arena);
}

void render() {
    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, md_array_bytes(vertices), vertices, GL_STREAM_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, md_array_bytes(indices), indices, GL_STREAM_DRAW);

    glEnable(GL_PROGRAM_POINT_SIZE);
    //glEnable(GL_BLEND);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glUseProgram(program);
    glLineWidth(1.f);

    glUniform1f(uniform_loc_point_size, 1.f);

    int current_view_matrix_idx = -1;
    int current_proj_matrix_idx = -1;

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, default_tex);

    const GLenum index_type = sizeof(Index) == 2 ? GL_UNSIGNED_SHORT : GL_UNSIGNED_INT;

    for (size_t i = 0; i < md_array_size(commands); ++i) {
        const auto& cmd = commands[i];
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
            mat3_t normal_matrix = mat3_from_mat4(mat4_transpose(mat4_inverse(matrix_stack[cmd.view_matrix_idx])));
            glUniformMatrix3fv(uniform_loc_normal_matrix, 1, GL_FALSE, &normal_matrix.elem[0][0]);
        }
        if (update_mvp) {
            mat4_t mvp_matrix = mat4_mul(matrix_stack[cmd.proj_matrix_idx], matrix_stack[cmd.view_matrix_idx]);
            glUniformMatrix4fv(uniform_loc_mvp_matrix, 1, GL_FALSE, &mvp_matrix.elem[0][0]);
        }
        // @TODO: Enable textures to be bound...

        glDrawElements(cmd.primitive_type, cmd.count, index_type, (const void*)(cmd.offset * sizeof(Index)));
    }

    glBindVertexArray(0);
    glUseProgram(0);

    glDisable(GL_PROGRAM_POINT_SIZE);
    //glDisable(GL_BLEND);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    md_array_shrink(vertices, 0);
    md_array_shrink(indices,  0);
    md_array_shrink(commands, 0);
    md_array_shrink(matrix_stack, 0);
    curr_view_matrix_idx = -1;
    curr_proj_matrix_idx = -1;
}

// PRIMITIVES
void draw_point(vec3_t pos, uint32_t color) {
    const Index idx = (Index)md_array_size(vertices);

    Vertex v = {pos, color};
    md_array_push(vertices, v,  arena);
    md_array_push(indices, idx, arena);

    append_draw_command(1, GL_POINTS);
}

void draw_points_v(const Vertex verts[], size_t count, vec4_t color_mult) {
    const Index idx = (Index)md_array_size(vertices);
    md_array_resize(vertices, md_array_size(vertices) + count, arena);
    Vertex* v = md_array_end(vertices) - count;

    ASSERT(v);

    for (size_t i = 0; i < count; ++i) {
        vec4_t c = convert_color(verts[i].color) * color_mult;
        v[i].coord = verts[i].coord;
        v[i].color = convert_color(c);
    }

    for (Index i = idx; i < idx + count; i += 1) {
        md_array_push(indices, i, arena);
    }

    append_draw_command(count, GL_POINTS);
}

void draw_line(vec3_t from, vec3_t to, uint32_t color) {
    const Index idx = (Index)md_array_size(vertices);

    Vertex v[2] = {
        {from, color},
        {to,   color}
    };
    md_array_push_array(vertices, v, 2, arena);

    md_array_push(indices, idx + 0, arena);
    md_array_push(indices, idx + 1, arena);

    append_draw_command(2, GL_LINES);
}

void draw_lines_v(const Vertex verts[], size_t count, vec4_t color_mult) {
    if ((count & 1) == 1) {
        MD_LOG_DEBUG("Expected even number of vertices as in function %s", __FUNCTION__);
    }

    // Make sure count is an even number
    count = count - (count & 1);

    const Index idx = (Index)md_array_size(vertices);
    md_array_resize(vertices, md_array_size(vertices) + count, arena);
    Vertex* v = md_array_end(vertices) - count;

    ASSERT(v);

    for (size_t i = 0; i < count; ++i) {
        vec4_t c = convert_color(verts[i].color) * color_mult;
        v[i].coord = verts[i].coord;
        v[i].color = convert_color(c);
    }

    for (Index i = idx; i < idx + count; i += 2) {
        md_array_push(indices, i + 0, arena);
        md_array_push(indices, i + 1, arena);
    }

    append_draw_command(count, GL_LINES);
}

void draw_triangle(vec3_t p0, vec3_t p1, vec3_t p2, uint32_t color) {
    const Index idx = (Index)md_array_size(vertices);
    //const vec3_t normal = vec3_normalize(vec3_cross(vec3_sub(p1, p0), vec3_sub(p2, p0)));

    Vertex v[3] = {
        {p0, color},
        {p1, color},
        {p2, color},
    };
    md_array_push_array(vertices, v, 3, arena);

    md_array_push(indices, idx + 0, arena);
    md_array_push(indices, idx + 1, arena);
    md_array_push(indices, idx + 2, arena);

    append_draw_command(3, GL_TRIANGLES);
}

void draw_triangles_v(const Vertex verts[], size_t count, vec4_t color_mult) {
    if ((count % 3) != 0) {
        MD_LOG_DEBUG("Expected number of vertices to be divisible by 3 in function %s", __FUNCTION__);
    }

    // Make sure count is a correct number here
    count = count - (count % 3);

    const Index idx = (Index)md_array_size(vertices);
    md_array_resize(vertices, md_array_size(vertices) + count, arena);
    Vertex* v = md_array_end(vertices) - count;

    ASSERT(v);

    for (size_t i = 0; i < count; ++i) {
        vec4_t c = convert_color(verts[i].color) * color_mult;
        v[i].coord = verts[i].coord;
        v[i].color = convert_color(c);
    }

    for (Index i = idx; i < idx + count; i += 3) {
        md_array_push(indices, i + 0, arena);
        md_array_push(indices, i + 1, arena);
        md_array_push(indices, i + 2, arena);
    }

    append_draw_command(count, GL_TRIANGLES);
}

void draw_plane(vec3_t center, vec3_t u, vec3_t v, uint32_t color) {
    const Index idx = (Index)md_array_size(vertices);
    //const vec3_t normal = vec3_normalize(vec3_cross(u, v));

    Vertex vert[4] = {
        {vec3_sub(center, vec3_add(u, v)), color},
        {vec3_sub(center, vec3_sub(u, v)), color},
        {vec3_add(center, vec3_add(u, v)), color},
        {vec3_add(center, vec3_sub(u, v)), color},
    };
    md_array_push_array(vertices, vert, 4, arena);

    md_array_push(indices, idx + 0, arena);
    md_array_push(indices, idx + 1, arena);
    md_array_push(indices, idx + 2, arena);
    md_array_push(indices, idx + 2, arena);
    md_array_push(indices, idx + 1, arena);
    md_array_push(indices, idx + 3, arena);

    append_draw_command(6, GL_TRIANGLES);
}

void draw_plane_wireframe(vec3_t center, vec3_t vec_u, vec3_t vec_v, uint32_t color, int segments_u, int segments_v) {
    ASSERT(segments_u > 0);
    ASSERT(segments_v > 0);

    //const vec3_t normal = vec3_normalize(vec3_cross(vec_u, vec_v));

    for (int i = 0; i <= segments_u; i++) {
        const float t = -1.0f + 2.0f * ((float)i / (float)segments_u);
        const vec3_t u = vec3_mul_f(vec_u, t);
        draw_line(vec3_sub(center, vec3_add(vec_v, u)), vec3_add(center, vec3_add(vec_v, u)), color);
    }

    for (int i = 0; i <= segments_v; i++) {
        const float t = -1.0f + 2.0f * ((float)i / (float)segments_v);
        const vec3_t v = vec3_mul_f(vec_v, t);
        draw_line(vec3_sub(center, vec3_add(vec_u, v)), vec3_add(center, vec3_add(vec_u, v)), color);
    }
}

/*
void draw_aabb(vec3_t min_box, vec3_t max_box) {
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

void draw_box_wireframe(vec3_t min_box, vec3_t max_box, uint32_t color) {
    // Z = min
    draw_line(vec3_t{min_box.elem[0], min_box.elem[1], min_box.elem[2]}, vec3_t{max_box.elem[0], min_box.elem[1], min_box.elem[2]}, color);
    draw_line(vec3_t{min_box.elem[0], min_box.elem[1], min_box.elem[2]}, vec3_t{min_box.elem[0], max_box.elem[1], min_box.elem[2]}, color);
    draw_line(vec3_t{max_box.elem[0], min_box.elem[1], min_box.elem[2]}, vec3_t{max_box.elem[0], max_box.elem[1], min_box.elem[2]}, color);
    draw_line(vec3_t{min_box.elem[0], max_box.elem[1], min_box.elem[2]}, vec3_t{max_box.elem[0], max_box.elem[1], min_box.elem[2]}, color);

    // Z = max
    draw_line(vec3_t{min_box.elem[0], min_box.elem[1], max_box.elem[2]}, vec3_t{max_box.elem[0], min_box.elem[1], max_box.elem[2]}, color);
    draw_line(vec3_t{min_box.elem[0], min_box.elem[1], max_box.elem[2]}, vec3_t{min_box.elem[0], max_box.elem[1], max_box.elem[2]}, color);
    draw_line(vec3_t{max_box.elem[0], min_box.elem[1], max_box.elem[2]}, vec3_t{max_box.elem[0], max_box.elem[1], max_box.elem[2]}, color);
    draw_line(vec3_t{min_box.elem[0], max_box.elem[1], max_box.elem[2]}, vec3_t{max_box.elem[0], max_box.elem[1], max_box.elem[2]}, color);

    // Z min max
    draw_line(vec3_t{min_box.elem[0], min_box.elem[1], min_box.elem[2]}, vec3_t{min_box.elem[0], min_box.elem[1], max_box.elem[2]}, color);
    draw_line(vec3_t{min_box.elem[0], max_box.elem[1], min_box.elem[2]}, vec3_t{min_box.elem[0], max_box.elem[1], max_box.elem[2]}, color);
    draw_line(vec3_t{max_box.elem[0], min_box.elem[1], min_box.elem[2]}, vec3_t{max_box.elem[0], min_box.elem[1], max_box.elem[2]}, color);
    draw_line(vec3_t{max_box.elem[0], max_box.elem[1], min_box.elem[2]}, vec3_t{max_box.elem[0], max_box.elem[1], max_box.elem[2]}, color);
}

void draw_box_wireframe(vec3_t min_box, vec3_t max_box, vec4_t color) {
    draw_box_wireframe(min_box, max_box, convert_color(color));
}

void draw_box_wireframe(vec3_t min_box, vec3_t max_box, mat4_t model_matrix, uint32_t color) {
    const mat3_t R = mat3_from_mat4(model_matrix);
    const vec3_t trans = vec3_from_vec4(model_matrix.col[3]);
    // Z = min
    draw_line(vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], min_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], min_box.elem[1], min_box.elem[2]})), color);
    draw_line(vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], min_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], max_box.elem[1], min_box.elem[2]})), color);
    draw_line(vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], min_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], max_box.elem[1], min_box.elem[2]})), color);
    draw_line(vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], max_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], max_box.elem[1], min_box.elem[2]})), color);
                                                                                                                                                                    
    // Z = max                                                                                                                                                      
    draw_line(vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], min_box.elem[1], max_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], min_box.elem[1], max_box.elem[2]})), color);
    draw_line(vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], min_box.elem[1], max_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], max_box.elem[1], max_box.elem[2]})), color);
    draw_line(vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], min_box.elem[1], max_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], max_box.elem[1], max_box.elem[2]})), color);
    draw_line(vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], max_box.elem[1], max_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], max_box.elem[1], max_box.elem[2]})), color);
                                                                                                                                                                     
    // Z min to max                                                                                                                                                 
    draw_line(vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], min_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], min_box.elem[1], max_box.elem[2]})), color);
    draw_line(vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], max_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], max_box.elem[1], max_box.elem[2]})), color);
    draw_line(vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], min_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], min_box.elem[1], max_box.elem[2]})), color);
    draw_line(vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], max_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], max_box.elem[1], max_box.elem[2]})), color);
}

void draw_box_wireframe(vec3_t min_box, vec3_t max_box, mat4_t model_matrix, vec4_t color) {
    draw_box_wireframe(min_box, max_box, model_matrix, convert_color(color));
}

void draw_basis(mat4_t basis, const float scale, uint32_t x_color, uint32_t y_color, uint32_t z_color) {
    const vec3_t o = vec3_from_vec4(basis.col[3]);
    const vec3_t x = vec3_add(o, vec3_mul_f(vec3_from_vec4(basis.col[0]), scale));
    const vec3_t y = vec3_add(o, vec3_mul_f(vec3_from_vec4(basis.col[1]), scale));
    const vec3_t z = vec3_add(o, vec3_mul_f(vec3_from_vec4(basis.col[2]), scale));

    draw_line(o, x, x_color);
    draw_line(o, y, y_color);
    draw_line(o, z, z_color);
}

void draw_basis(mat4_t basis, const float scale, vec4_t x_color, vec4_t y_color, vec4_t z_color) {
    draw_basis(basis, scale, convert_color(x_color), convert_color(y_color), convert_color(z_color));
}

}  // namespace immediate
