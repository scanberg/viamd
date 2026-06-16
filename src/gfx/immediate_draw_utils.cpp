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
    int32_t model_matrix_idx = -1;
    uint32_t picking_base_idx = 0;
    GLenum primitive_type;
};

struct PassParams {
    mat4_t view = mat4_ident();
    mat4_t proj = mat4_ident();
    uint32_t picking_base_idx = 0;
};

struct Queue {
    const char* label = nullptr;

    md_allocator_i* arena = nullptr;
    md_array(mat4_t) model_matrices = nullptr;
    md_array(DrawCommand) commands = nullptr;
    md_array(Vertex) vertices = nullptr;
    md_array(Index) indices = nullptr;

    int32_t curr_model_matrix_idx = -1;
    uint32_t curr_picking_base_idx = 0;
};

static GLuint vbo = 0;
static GLuint ibo = 0;
static GLuint vao = 0;
static GLuint default_tex = 0;

static GLuint program = 0;
static GLuint program_shaded = 0;

static GLint uniform_loc_mvp_matrix = -1;
static GLint uniform_loc_normal_matrix = -1;
static GLint uniform_loc_uv_scale = -1;
static GLint uniform_loc_point_size = -1;
static GLint uniform_loc_picking_base_idx = -1;

static const char* v_shader_src = R"(
#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_mvp_matrix;
uniform mat3 u_normal_matrix;
uniform float u_point_size = 1.f;
uniform uint u_picking_base_idx = 0u;

layout(location = 0) in vec3 in_position;
layout(location = 1) in vec4 in_color;
layout(location = 2) in vec3 in_normal;
layout(location = 3) in uint in_picking_idx;

out vec4 color;
out vec3 normal;
flat out uint picking_idx;

void main() {
    gl_Position = u_mvp_matrix * vec4(in_position, 1.0);
    gl_PointSize = max(u_point_size, 200.0 / gl_Position.w);
    color = in_color;
    normal = u_normal_matrix * in_normal;
    picking_idx = (in_picking_idx == 0xFFFFFFFFu) ? 0xFFFFFFFFu : in_picking_idx + u_picking_base_idx;
}
)";

static const char* f_shader_src = R"(
#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

in vec4 color;
in vec3 normal;
flat in uint picking_idx;

layout(location = 0) out vec4 out_color_alpha;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec2 out_velocity;
layout(location = 3) out vec4 out_picking_idx;

vec4 encode_normal(vec3 n) {
    float p = sqrt(n.z * 8.0 + 8.0);
    return vec4(n.xy / p + 0.5, 0.0, 0.0);
}

vec4 encode_u32(uint index) {
    return vec4(
        (index & 0x000000FFU) >> 0U,
        (index & 0x0000FF00U) >> 8U,
        (index & 0x00FF0000U) >> 16U,
        (index & 0xFF000000U) >> 24U) / 255.0;
}

void main() {
    out_color_alpha = color;
    out_normal = encode_normal(normal);
    out_velocity = vec2(0,0);
    out_picking_idx = encode_u32(picking_idx);
}
)";

static inline void append_draw_command(Queue& ctx, size_t count, GLenum primitive_type) {
    const uint32_t max_size = (sizeof(Index) == 2 ? 0xFFFFU : 0xFFFFFFFFU);
    ASSERT(md_array_size(ctx.indices) + count < max_size);
    DrawCommand* last_cmd = md_array_last(ctx.commands);
    if (last_cmd &&
        last_cmd->primitive_type == primitive_type &&
        last_cmd->model_matrix_idx == ctx.curr_model_matrix_idx &&
        last_cmd->picking_base_idx == ctx.curr_picking_base_idx) {
        size_t capacity = max_size - last_cmd->count;
        if (count < capacity) {
            last_cmd->count += (uint32_t)count;
            return;
        }
        last_cmd->count += (uint32_t)capacity;
        count -= capacity;
    }

    ASSERT(ctx.curr_model_matrix_idx > -1 && "Immediate draw model matrix not set");
    const size_t offset = md_array_size(ctx.indices) - count;
    DrawCommand cmd {(uint32_t)offset, (uint32_t)count, ctx.curr_model_matrix_idx, ctx.curr_picking_base_idx, primitive_type};
    md_array_push(ctx.commands, cmd, ctx.arena);
}

static void init_queue_storage(Queue& q, md_allocator_i* arena) {
    MEMSET(&q, 0, sizeof(Queue));
    q.arena = arena;
    md_array_ensure(q.vertices, 100000, q.arena);
    md_array_ensure(q.indices, 100000, q.arena);
    md_array_ensure(q.model_matrices, 10, q.arena);
    md_array_ensure(q.commands, 32, q.arena);
    q.curr_model_matrix_idx = -1;
    q.curr_picking_base_idx = 0;
}

void encoder_set_model(Queue& ctx, const mat4_t& model_matrix);
void encoder_set_picking_base_idx(Queue& ctx, uint32_t base_idx);
void encoder_submit(Queue& ctx, const mat4_t& view, const mat4_t& proj);
void draw_point(Queue& ctx, vec3_t pos, uint32_t color, uint32_t picking_idx);
void draw_points_v(Queue& ctx, const Vertex verts[], size_t count, vec4_t color_mult);
void draw_line(Queue& ctx, vec3_t from, vec3_t to, uint32_t color);
void draw_lines_v(Queue& ctx, const Vertex verts[], size_t count, vec4_t color_mult);
void draw_triangle(Queue& ctx, vec3_t p0, vec3_t p1, vec3_t p2, uint32_t color, uint32_t picking_idx);
void draw_triangles_v(Queue& ctx, const Vertex verts[], size_t count, vec4_t color_mult);
void draw_plane(Queue& ctx, vec3_t center, vec3_t u, vec3_t v, uint32_t color, uint32_t picking_idx);
void draw_plane_wireframe(Queue& ctx, vec3_t center, vec3_t vec_u, vec3_t vec_v, uint32_t color, int segments_u, int segments_v);
void draw_box_wireframe(Queue& ctx, vec3_t min_box, vec3_t max_box, uint32_t color);
void draw_box_wireframe(Queue& ctx, vec3_t min_box, vec3_t max_box, vec4_t color);
void draw_box_wireframe(Queue& ctx, vec3_t min_box, vec3_t max_box, mat4_t model_matrix, uint32_t color);
void draw_box_wireframe(Queue& ctx, vec3_t min_box, vec3_t max_box, mat4_t model_matrix, vec4_t color);
void draw_basis(Queue& ctx, mat4_t basis, const float scale, uint32_t x_color, uint32_t y_color, uint32_t z_color);
void draw_basis(Queue& ctx, mat4_t basis, const float scale, vec4_t x_color, vec4_t y_color, vec4_t z_color);
void draw_sphere(Queue& ctx, vec3_t center, float radius, uint32_t color, uint32_t picking_idx, int stacks, int slices);
void draw_sphere_wireframe(Queue& ctx, vec3_t center, float radius, uint32_t color, int stacks, int slices);
void draw_cylinder(Queue& ctx, vec3_t from, vec3_t to, float radius, uint32_t color, uint32_t picking_idx, int segments);
void draw_cylinder_wireframe(Queue& ctx, vec3_t from, vec3_t to, float radius, uint32_t color, int segments);
void draw_cone(Queue& ctx, vec3_t base, vec3_t tip, float radius, uint32_t color, uint32_t picking_idx, int segments);
void draw_cone_wireframe(Queue& ctx, vec3_t base, vec3_t tip, float radius, uint32_t color, int segments);
void draw_box(Queue& ctx, vec3_t mn, vec3_t mx, uint32_t color, uint32_t picking_idx);

// All queue APIs require an explicit queue; this helper enforces that contract
// and lazily injects an identity model matrix for first use.
static Queue& require_queue(Queue* queue) {
    ASSERT(queue);
    if (queue->curr_model_matrix_idx < 0) {
        queue->curr_model_matrix_idx = (int32_t)md_array_size(queue->model_matrices);
        md_array_push(queue->model_matrices, mat4_ident(), queue->arena);
    }
    return *queue;
}

void set_model(Queue* queue, const mat4_t& model_mat) {
    Queue& q = require_queue(queue);
    encoder_set_model(q, model_mat);
}

void set_picking_base_idx(Queue* queue, uint32_t base_idx) {
    Queue& q = require_queue(queue);
    encoder_set_picking_base_idx(q, base_idx);
}

void point(Queue* queue, vec3_t pos, uint32_t color, uint32_t picking_idx) {
    Queue& q = require_queue(queue);
    draw_point(q, pos, color, picking_idx);
}

void line(Queue* queue, vec3_t from, vec3_t to, uint32_t color) {
    Queue& q = require_queue(queue);
    draw_line(q, from, to, color);
}

void triangle(Queue* queue, vec3_t v0, vec3_t v1, vec3_t v2, uint32_t color, uint32_t picking_idx) {
    Queue& q = require_queue(queue);
    draw_triangle(q, v0, v1, v2, color, picking_idx);
}

void triangle_wireframe(Queue* queue, vec3_t v0, vec3_t v1, vec3_t v2, uint32_t color) {
    Queue& q = require_queue(queue);
    draw_line(q, v0, v1, color);
    draw_line(q, v1, v2, color);
    draw_line(q, v2, v0, color);
}

void sphere(Queue* queue, vec3_t center, float radius, uint32_t color, uint32_t picking_idx, int stacks, int slices) {
    Queue& q = require_queue(queue);
    draw_sphere(q, center, radius, color, picking_idx, stacks, slices);
}

void sphere_wireframe(Queue* queue, vec3_t center, float radius, uint32_t color, int stacks, int slices) {
    Queue& q = require_queue(queue);
    draw_sphere_wireframe(q, center, radius, color, stacks, slices);
}

void cylinder(Queue* queue, vec3_t from, vec3_t to, float radius, uint32_t color, uint32_t picking_idx, int segments) {
    Queue& q = require_queue(queue);
    draw_cylinder(q, from, to, radius, color, picking_idx, segments);
}

void cylinder_wireframe(Queue* queue, vec3_t from, vec3_t to, float radius, uint32_t color, int segments) {
    Queue& q = require_queue(queue);
    draw_cylinder_wireframe(q, from, to, radius, color, segments);
}

void cone(Queue* queue, vec3_t base, vec3_t tip, float radius, uint32_t color, uint32_t picking_idx, int segments) {
    Queue& q = require_queue(queue);
    draw_cone(q, base, tip, radius, color, picking_idx, segments);
}

void cone_wireframe(Queue* queue, vec3_t base, vec3_t tip, float radius, uint32_t color, int segments) {
    Queue& q = require_queue(queue);
    draw_cone_wireframe(q, base, tip, radius, color, segments);
}

void capsule(Queue* queue, vec3_t from, vec3_t to, float radius, uint32_t color, uint32_t picking_idx, int segments) {
    Queue& q = require_queue(queue);
    draw_cylinder(q, from, to, radius, color, picking_idx, segments);
    draw_sphere(q, from, radius, color, picking_idx, segments, segments * 2);
    draw_sphere(q, to, radius, color, picking_idx, segments, segments * 2);
}

void capsule_wireframe(Queue* queue, vec3_t from, vec3_t to, float radius, uint32_t color, int segments) {
    Queue& q = require_queue(queue);
    draw_cylinder_wireframe(q, from, to, radius, color, segments);
    draw_sphere_wireframe(q, from, radius, color, segments, segments * 2);
    draw_sphere_wireframe(q, to, radius, color, segments, segments * 2);
}

void box(Queue* queue, vec3_t min_box, vec3_t max_box, uint32_t color, uint32_t picking_idx) {
    Queue& q = require_queue(queue);
    draw_box(q, min_box, max_box, color, picking_idx);
}

void points(Queue* queue, const Vertex verts[], size_t count, vec4_t color_mult) {
    Queue& q = require_queue(queue);
    draw_points_v(q, verts, count, color_mult);
}

void lines(Queue* queue, const Vertex verts[], size_t count, vec4_t color_mult) {
    Queue& q = require_queue(queue);
    draw_lines_v(q, verts, count, color_mult);
}

void triangles(Queue* queue, const Vertex verts[], size_t count, vec4_t color_mult) {
    Queue& q = require_queue(queue);
    draw_triangles_v(q, verts, count, color_mult);
}

void triangles_wireframe(Queue* queue, const Vertex verts[], size_t count, vec4_t color_mult) {
    Queue& q = require_queue(queue);
    if (!verts) return;
    for (size_t i = 0; i + 2 < count; i += 3) {
        const uint32_t c0 = u32_from_vec4(vec4_mul(vec4_from_u32(verts[i + 0].color), color_mult));
        const uint32_t c1 = u32_from_vec4(vec4_mul(vec4_from_u32(verts[i + 1].color), color_mult));
        const uint32_t c2 = u32_from_vec4(vec4_mul(vec4_from_u32(verts[i + 2].color), color_mult));
        draw_line(q, verts[i + 0].coord, verts[i + 1].coord, c0);
        draw_line(q, verts[i + 1].coord, verts[i + 2].coord, c1);
        draw_line(q, verts[i + 2].coord, verts[i + 0].coord, c2);
    }
}

void plane(Queue* queue, vec3_t center, vec3_t plane_u, vec3_t plane_v, uint32_t color, uint32_t picking_idx) {
    Queue& q = require_queue(queue);
    draw_plane(q, center, plane_u, plane_v, color, picking_idx);
}

void plane_wireframe(Queue* queue, vec3_t center, vec3_t plane_u, vec3_t plane_v, uint32_t color, int segments_u, int segments_v) {
    Queue& q = require_queue(queue);
    draw_plane_wireframe(q, center, plane_u, plane_v, color, segments_u, segments_v);
}

void box_wireframe(Queue* queue, vec3_t min_box, vec3_t max_box, uint32_t color) {
    Queue& q = require_queue(queue);
    draw_box_wireframe(q, min_box, max_box, color);
}

void box_wireframe(Queue* queue, vec3_t min_box, vec3_t max_box, vec4_t color) {
    Queue& q = require_queue(queue);
    draw_box_wireframe(q, min_box, max_box, color);
}

void box_wireframe(Queue* queue, vec3_t min_box, vec3_t max_box, mat4_t model_matrix, uint32_t color) {
    Queue& q = require_queue(queue);
    draw_box_wireframe(q, min_box, max_box, model_matrix, color);
}

void box_wireframe(Queue* queue, vec3_t min_box, vec3_t max_box, mat4_t model_matrix, vec4_t color) {
    Queue& q = require_queue(queue);
    draw_box_wireframe(q, min_box, max_box, model_matrix, color);
}

void basis(Queue* queue, mat4_t basis_mat, float scale, uint32_t x_color, uint32_t y_color, uint32_t z_color) {
    Queue& q = require_queue(queue);
    draw_basis(q, basis_mat, scale, x_color, y_color, z_color);
}

void basis(Queue* queue, mat4_t basis_mat, float scale, vec4_t x_color, vec4_t y_color, vec4_t z_color) {
    Queue& q = require_queue(queue);
    draw_basis(q, basis_mat, scale, x_color, y_color, z_color);
}

Scope::Scope(Queue* supplied_queue, const char* label) {
    (void)label;
    ASSERT(supplied_queue);
    queue = supplied_queue;
    // Begin scope: push a fresh identity model transform in the queue stream.
    set_model(queue, mat4_ident());
    set_picking_base_idx(queue, 0);
}

Scope::~Scope() {
    if (queue) {
        // End scope: reset back to identity for the following commands.
        set_model(queue, mat4_ident());
        set_picking_base_idx(queue, 0);
    }
    queue = nullptr;
}

Scope::Scope(Scope&& other) noexcept {
    queue = other.queue;
    other.queue = nullptr;
}

Scope& Scope::operator=(Scope&& other) noexcept {
    if (this != &other) {
        queue = other.queue;
        other.queue = nullptr;
    }
    return *this;
}

void set_model(Scope& s, const mat4_t& model_mat) { set_model(s.queue, model_mat); }
void set_picking_base_idx(Scope& s, uint32_t base_idx) { set_picking_base_idx(s.queue, base_idx); }

Queue* queue_create(const char* label) {
    md_allocator_i* heap = md_get_heap_allocator();
    Queue* queue = (Queue*)md_alloc(heap, sizeof(Queue));
    md_allocator_i* arena = md_arena_allocator_create(heap, MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE);
    init_queue_storage(*queue, arena);
    queue->label = label;
    return queue;
}

void queue_destroy(Queue* queue) {
    if (!queue) return;
    if (queue->arena) {
        md_arena_allocator_destroy(queue->arena);
    }
    md_free(md_get_heap_allocator(), queue, sizeof(Queue));
}

void queue_reset(Queue* queue) {
    if (!queue) return;
    md_array_shrink(queue->vertices, 0);
    md_array_shrink(queue->indices, 0);
    md_array_shrink(queue->commands, 0);
    md_array_shrink(queue->model_matrices, 0);
    queue->curr_model_matrix_idx = -1;
    queue->curr_picking_base_idx = 0;
}

void render(Queue* queue, const RenderParams& params) {
    if (!queue) return;

    if (md_array_size(queue->commands) == 0) {
        return;
    }

    encoder_submit(*queue, params.view, params.proj);
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
    uniform_loc_picking_base_idx = glGetUniformLocation(program, "u_picking_base_idx");

    glGenBuffers(1, &vbo);
    glGenBuffers(1, &ibo);

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)offsetof(Vertex, coord));

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (const GLvoid*)offsetof(Vertex, color));

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)offsetof(Vertex, normal));

    glEnableVertexAttribArray(3);
    glVertexAttribIPointer(3, 1, GL_UNSIGNED_INT, sizeof(Vertex), (const GLvoid*)offsetof(Vertex, picking_idx));

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
}

void shutdown() {
    if (vbo) glDeleteBuffers(1, &vbo);
    if (ibo) glDeleteBuffers(1, &ibo);
    if (vao) glDeleteVertexArrays(1, &vao);
    if (program) glDeleteProgram(program);
    if (program_shaded) glDeleteProgram(program_shaded);
}

void encoder_set_model(Queue& ctx, const mat4_t& model_matrix) {
    ctx.curr_model_matrix_idx = (int)md_array_size(ctx.model_matrices);
    md_array_push(ctx.model_matrices, model_matrix, ctx.arena);
}

void encoder_set_picking_base_idx(Queue& ctx, uint32_t base_idx) {
    ctx.curr_picking_base_idx = base_idx;
}

void encoder_submit(Queue& ctx, const mat4_t& view, const mat4_t& proj) {
	size_t num_commands = md_array_size(ctx.commands);
	if (num_commands == 0) {
		return;
	}

	glBindVertexArray(vao);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, md_array_bytes(ctx.vertices), ctx.vertices, GL_STREAM_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, md_array_bytes(ctx.indices), ctx.indices, GL_STREAM_DRAW);

	glEnable(GL_PROGRAM_POINT_SIZE);
	glLineWidth(1.f);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, default_tex);

	const GLenum index_type = sizeof(Index) == 2 ? GL_UNSIGNED_SHORT : GL_UNSIGNED_INT;
	const mat4_t vp_matrix = mat4_mul(proj, view);

	int current_model_matrix_idx = -1;

    glUseProgram(program);

	for (size_t i = 0; i < num_commands; ++i) {
		const auto& cmd = ctx.commands[i];

		if (cmd.model_matrix_idx != current_model_matrix_idx) {
			current_model_matrix_idx = cmd.model_matrix_idx;
			const mat4_t mv_matrix = mat4_mul(view, ctx.model_matrices[cmd.model_matrix_idx]);
			mat3_t normal_matrix = mat3_from_mat4(mat4_transpose(mat4_inverse(mv_matrix)));
			mat4_t mvp_matrix = mat4_mul(vp_matrix, ctx.model_matrices[cmd.model_matrix_idx]);

		    glUniformMatrix3fv(uniform_loc_normal_matrix, 1, GL_FALSE, &normal_matrix.elem[0][0]);
			glUniformMatrix4fv(uniform_loc_mvp_matrix, 1, GL_FALSE, &mvp_matrix.elem[0][0]);
		}

		glUniform1ui(uniform_loc_picking_base_idx, cmd.picking_base_idx);

		glDrawElements(cmd.primitive_type, cmd.count, index_type, (const void*)(cmd.offset * sizeof(Index)));
	}

	glBindVertexArray(0);
	glUseProgram(0);

	glDisable(GL_PROGRAM_POINT_SIZE);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	md_array_shrink(ctx.vertices, 0);
	md_array_shrink(ctx.indices,  0);
	md_array_shrink(ctx.commands, 0);
	md_array_shrink(ctx.model_matrices, 0);
	ctx.curr_model_matrix_idx = -1;
}

// PRIMITIVES
void draw_point(Queue& ctx, vec3_t pos, uint32_t color, uint32_t picking_idx) {
    const Index idx = (Index)md_array_size(ctx.vertices);

    Vertex v = {pos, color, {0,0,1}, picking_idx};
    md_array_push(ctx.vertices, v,  ctx.arena);
    md_array_push(ctx.indices, idx, ctx.arena);

    append_draw_command(ctx, 1, GL_POINTS);
}

void draw_points_v(Queue& ctx, const Vertex verts[], size_t count, vec4_t color_mult) {
    const Index idx = (Index)md_array_size(ctx.vertices);
    md_array_resize(ctx.vertices, md_array_size(ctx.vertices) + count, ctx.arena);
    Vertex* v = md_array_end(ctx.vertices) - count;

    ASSERT(v);

    for (size_t i = 0; i < count; ++i) {
        vec4_t c = convert_color(verts[i].color) * color_mult;
        v[i].coord = verts[i].coord;
        v[i].color = convert_color(c);
        v[i].normal = verts[i].normal;
        v[i].picking_idx = verts[i].picking_idx;
    }

    for (Index i = idx; i < idx + count; i += 1) {
        md_array_push(ctx.indices, i, ctx.arena);
    }

    append_draw_command(ctx, count, GL_POINTS);
}

void draw_line(Queue& ctx, vec3_t from, vec3_t to, uint32_t color) {
    const Index idx = (Index)md_array_size(ctx.vertices);

    Vertex v[2] = {
        {from, color, {0,0,1}, 0xFFFFFFFF},
        {to,   color, {0,0,1}, 0xFFFFFFFF}
    };
    md_array_push_array(ctx.vertices, v, 2, ctx.arena);

    md_array_push(ctx.indices, idx + 0, ctx.arena);
    md_array_push(ctx.indices, idx + 1, ctx.arena);

    append_draw_command(ctx, 2, GL_LINES);
}

void draw_lines_v(Queue& ctx, const Vertex verts[], size_t count, vec4_t color_mult) {
    if ((count & 1) == 1) {
        MD_LOG_DEBUG("Expected even number of vertices as in function %s", __FUNCTION__);
    }

    // Make sure count is an even number
    count = count - (count & 1);

    const Index idx = (Index)md_array_size(ctx.vertices);
    md_array_resize(ctx.vertices, md_array_size(ctx.vertices) + count, ctx.arena);
    Vertex* v = md_array_end(ctx.vertices) - count;

    ASSERT(v);

    for (size_t i = 0; i < count; ++i) {
        vec4_t c = convert_color(verts[i].color) * color_mult;
        v[i].coord = verts[i].coord;
        v[i].color = convert_color(c);
        v[i].normal = verts[i].normal;
        v[i].picking_idx = 0xFFFFFFFFu;
    }

    for (Index i = idx; i < idx + count; i += 2) {
        md_array_push(ctx.indices, i + 0, ctx.arena);
        md_array_push(ctx.indices, i + 1, ctx.arena);
    }

    append_draw_command(ctx, count, GL_LINES);
}

void draw_triangle(Queue& ctx, vec3_t p0, vec3_t p1, vec3_t p2, uint32_t color, uint32_t picking_idx) {
    const Index idx = (Index)md_array_size(ctx.vertices);
    const vec3_t normal = vec3_normalize(vec3_cross(vec3_sub(p1, p0), vec3_sub(p2, p0)));

    Vertex v[3] = {
        {p0, color, normal, picking_idx},
        {p1, color, normal, picking_idx},
        {p2, color, normal, picking_idx},
    };
    md_array_push_array(ctx.vertices, v, 3, ctx.arena);

    md_array_push(ctx.indices, idx + 0, ctx.arena);
    md_array_push(ctx.indices, idx + 1, ctx.arena);
    md_array_push(ctx.indices, idx + 2, ctx.arena);

    append_draw_command(ctx, 3, GL_TRIANGLES);
}

void draw_triangles_v(Queue& ctx, const Vertex verts[], size_t count, vec4_t color_mult) {
    if ((count % 3) != 0) {
        MD_LOG_DEBUG("Expected number of vertices to be divisible by 3 in function %s", __FUNCTION__);
    }

    // Make sure count is a correct number here
    count = count - (count % 3);

    const Index idx = (Index)md_array_size(ctx.vertices);
    md_array_resize(ctx.vertices, md_array_size(ctx.vertices) + count, ctx.arena);
    Vertex* v = md_array_end(ctx.vertices) - count;

    ASSERT(v);

    for (size_t i = 0; i < count; ++i) {
        vec4_t c = convert_color(verts[i].color) * color_mult;
        v[i].coord = verts[i].coord;
        v[i].color = convert_color(c);
        v[i].normal = verts[i].normal;
        v[i].picking_idx = verts[i].picking_idx;
    }

    for (Index i = idx; i < idx + count; i += 3) {
        md_array_push(ctx.indices, i + 0, ctx.arena);
        md_array_push(ctx.indices, i + 1, ctx.arena);
        md_array_push(ctx.indices, i + 2, ctx.arena);
    }

    append_draw_command(ctx, count, GL_TRIANGLES);
}

void draw_plane(Queue& ctx, vec3_t center, vec3_t u, vec3_t v, uint32_t color, uint32_t picking_idx) {
    const Index idx = (Index)md_array_size(ctx.vertices);
    const vec3_t normal = vec3_normalize(vec3_cross(u, v));

    Vertex vert[4] = {
        {vec3_sub(center, vec3_add(u, v)), color, normal, picking_idx},
        {vec3_sub(center, vec3_sub(u, v)), color, normal, picking_idx},
        {vec3_add(center, vec3_add(u, v)), color, normal, picking_idx},
        {vec3_add(center, vec3_sub(u, v)), color, normal, picking_idx},
    };
    md_array_push_array(ctx.vertices, vert, 4, ctx.arena);

    md_array_push(ctx.indices, idx + 0, ctx.arena);
    md_array_push(ctx.indices, idx + 1, ctx.arena);
    md_array_push(ctx.indices, idx + 2, ctx.arena);
    md_array_push(ctx.indices, idx + 2, ctx.arena);
    md_array_push(ctx.indices, idx + 1, ctx.arena);
    md_array_push(ctx.indices, idx + 3, ctx.arena);

    append_draw_command(ctx, 6, GL_TRIANGLES);
}

void draw_plane_wireframe(Queue& ctx, vec3_t center, vec3_t vec_u, vec3_t vec_v, uint32_t color, int segments_u, int segments_v) {
    ASSERT(segments_u > 0);
    ASSERT(segments_v > 0);

    //const vec3_t normal = vec3_normalize(vec3_cross(vec_u, vec_v));

    for (int i = 0; i <= segments_u; i++) {
        const float t = -1.0f + 2.0f * ((float)i / (float)segments_u);
        const vec3_t u = vec3_mul1(vec_u, t);
        draw_line(ctx, vec3_sub(center, vec3_add(vec_v, u)), vec3_add(center, vec3_add(vec_v, u)), color);
    }

    for (int i = 0; i <= segments_v; i++) {
        const float t = -1.0f + 2.0f * ((float)i / (float)segments_v);
        const vec3_t v = vec3_mul1(vec_v, t);
        draw_line(ctx, vec3_sub(center, vec3_add(vec_u, v)), vec3_add(center, vec3_add(vec_u, v)), color);
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

    append_draw_command(ctx, idx, 14, GL_TRIANGLE_STRIP);
}
*/

void draw_box_wireframe(Queue& ctx, vec3_t min_box, vec3_t max_box, uint32_t color) {
    // Z = min
    draw_line(ctx, vec3_t{min_box.elem[0], min_box.elem[1], min_box.elem[2]}, vec3_t{max_box.elem[0], min_box.elem[1], min_box.elem[2]}, color);
    draw_line(ctx, vec3_t{min_box.elem[0], min_box.elem[1], min_box.elem[2]}, vec3_t{min_box.elem[0], max_box.elem[1], min_box.elem[2]}, color);
    draw_line(ctx, vec3_t{max_box.elem[0], min_box.elem[1], min_box.elem[2]}, vec3_t{max_box.elem[0], max_box.elem[1], min_box.elem[2]}, color);
    draw_line(ctx, vec3_t{min_box.elem[0], max_box.elem[1], min_box.elem[2]}, vec3_t{max_box.elem[0], max_box.elem[1], min_box.elem[2]}, color);

    // Z = max
    draw_line(ctx, vec3_t{min_box.elem[0], min_box.elem[1], max_box.elem[2]}, vec3_t{max_box.elem[0], min_box.elem[1], max_box.elem[2]}, color);
    draw_line(ctx, vec3_t{min_box.elem[0], min_box.elem[1], max_box.elem[2]}, vec3_t{min_box.elem[0], max_box.elem[1], max_box.elem[2]}, color);
    draw_line(ctx, vec3_t{max_box.elem[0], min_box.elem[1], max_box.elem[2]}, vec3_t{max_box.elem[0], max_box.elem[1], max_box.elem[2]}, color);
    draw_line(ctx, vec3_t{min_box.elem[0], max_box.elem[1], max_box.elem[2]}, vec3_t{max_box.elem[0], max_box.elem[1], max_box.elem[2]}, color);

    // Z min max
    draw_line(ctx, vec3_t{min_box.elem[0], min_box.elem[1], min_box.elem[2]}, vec3_t{min_box.elem[0], min_box.elem[1], max_box.elem[2]}, color);
    draw_line(ctx, vec3_t{min_box.elem[0], max_box.elem[1], min_box.elem[2]}, vec3_t{min_box.elem[0], max_box.elem[1], max_box.elem[2]}, color);
    draw_line(ctx, vec3_t{max_box.elem[0], min_box.elem[1], min_box.elem[2]}, vec3_t{max_box.elem[0], min_box.elem[1], max_box.elem[2]}, color);
    draw_line(ctx, vec3_t{max_box.elem[0], max_box.elem[1], min_box.elem[2]}, vec3_t{max_box.elem[0], max_box.elem[1], max_box.elem[2]}, color);
}

void draw_box_wireframe(Queue& ctx, vec3_t min_box, vec3_t max_box, vec4_t color) {
    draw_box_wireframe(ctx, min_box, max_box, convert_color(color));
}

void draw_box_wireframe(Queue& ctx, vec3_t min_box, vec3_t max_box, mat4_t model_matrix, uint32_t color) {
    const mat3_t R = mat3_from_mat4(model_matrix);
    const vec3_t trans = vec3_from_vec4(model_matrix.col[3]);
    // Z = min
    draw_line(ctx, vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], min_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], min_box.elem[1], min_box.elem[2]})), color);
    draw_line(ctx, vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], min_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], max_box.elem[1], min_box.elem[2]})), color);
    draw_line(ctx, vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], min_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], max_box.elem[1], min_box.elem[2]})), color);
    draw_line(ctx, vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], max_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], max_box.elem[1], min_box.elem[2]})), color);
                                                                                                                                                                    
    // Z = max                                                                                                                                                      
    draw_line(ctx, vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], min_box.elem[1], max_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], min_box.elem[1], max_box.elem[2]})), color);
    draw_line(ctx, vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], min_box.elem[1], max_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], max_box.elem[1], max_box.elem[2]})), color);
    draw_line(ctx, vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], min_box.elem[1], max_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], max_box.elem[1], max_box.elem[2]})), color);
    draw_line(ctx, vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], max_box.elem[1], max_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], max_box.elem[1], max_box.elem[2]})), color);
                                                                                                                                                                     
    // Z min to max                                                                                                                                                 
    draw_line(ctx, vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], min_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], min_box.elem[1], max_box.elem[2]})), color);
    draw_line(ctx, vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], max_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{min_box.elem[0], max_box.elem[1], max_box.elem[2]})), color);
    draw_line(ctx, vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], min_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], min_box.elem[1], max_box.elem[2]})), color);
    draw_line(ctx, vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], max_box.elem[1], min_box.elem[2]})), vec3_add(trans, mat3_mul_vec3(R, vec3_t{max_box.elem[0], max_box.elem[1], max_box.elem[2]})), color);
}

void draw_box_wireframe(Queue& ctx, vec3_t min_box, vec3_t max_box, mat4_t model_matrix, vec4_t color) {
    draw_box_wireframe(ctx, min_box, max_box, model_matrix, convert_color(color));
}

void draw_basis(Queue& ctx, mat4_t basis, const float scale, uint32_t x_color, uint32_t y_color, uint32_t z_color) {
    const vec3_t o = vec3_from_vec4(basis.col[3]);
    const vec3_t x = vec3_add(o, vec3_mul1(vec3_from_vec4(basis.col[0]), scale));
    const vec3_t y = vec3_add(o, vec3_mul1(vec3_from_vec4(basis.col[1]), scale));
    const vec3_t z = vec3_add(o, vec3_mul1(vec3_from_vec4(basis.col[2]), scale));

    draw_line(ctx, o, x, x_color);
    draw_line(ctx, o, y, y_color);
    draw_line(ctx, o, z, z_color);
}

void draw_basis(Queue& ctx, mat4_t basis, const float scale, vec4_t x_color, vec4_t y_color, vec4_t z_color) {
    draw_basis(ctx, basis, scale, convert_color(x_color), convert_color(y_color), convert_color(z_color));
}

// Build an orthonormal tangent frame given a normalized axis direction.
// Returns two vectors (out_u, out_v) perpendicular to 'axis' and to each other.
static void build_tangent_frame(vec3_t axis, vec3_t& out_u, vec3_t& out_v) {
    vec3_t up = (fabsf(axis.x) < 0.9f) ? vec3_t{1, 0, 0} : vec3_t{0, 1, 0};
    out_u = vec3_normalize(vec3_cross(up, axis));
    out_v = vec3_cross(axis, out_u);
}

void draw_sphere(Queue& ctx, vec3_t center, float radius, uint32_t color, uint32_t picking_idx, int stacks, int slices) {
    ASSERT(stacks >= 2 && slices >= 3);
    const float pi  = 3.14159265358979323846f;
    const float tau = 2.0f * pi;

    // Emit a quad (two triangles) directly so we can assign smooth per-vertex normals.
    // On a unit sphere the outward normal at a point equals the point's position relative to center.
    auto emit_quad = [&](vec3_t n00, vec3_t n10, vec3_t n01, vec3_t n11) {
        vec3_t p00 = vec3_add(center, vec3_mul1(n00, radius));
        vec3_t p10 = vec3_add(center, vec3_mul1(n10, radius));
        vec3_t p01 = vec3_add(center, vec3_mul1(n01, radius));
        vec3_t p11 = vec3_add(center, vec3_mul1(n11, radius));
        const Index base = (Index)md_array_size(ctx.vertices);
        Vertex v[4] = {
            {p00, color, n00, picking_idx},
            {p10, color, n10, picking_idx},
            {p01, color, n01, picking_idx},
            {p11, color, n11, picking_idx},
        };
        md_array_push_array(ctx.vertices, v, 4, ctx.arena);
        // tri 0: p00, p11, p10
        md_array_push(ctx.indices, base + 0, ctx.arena);
        md_array_push(ctx.indices, base + 3, ctx.arena);
        md_array_push(ctx.indices, base + 1, ctx.arena);
        // tri 1: p00, p01, p11
        md_array_push(ctx.indices, base + 0, ctx.arena);
        md_array_push(ctx.indices, base + 2, ctx.arena);
        md_array_push(ctx.indices, base + 3, ctx.arena);
        append_draw_command(ctx, 6, GL_TRIANGLES);
    };

    for (int i = 0; i < stacks; ++i) {
        const float phi0 = pi * ((float) i      / (float)stacks) - pi * 0.5f;
        const float phi1 = pi * ((float)(i + 1) / (float)stacks) - pi * 0.5f;

        const float cp0 = cosf(phi0), sp0 = sinf(phi0);
        const float cp1 = cosf(phi1), sp1 = sinf(phi1);

        for (int j = 0; j < slices; ++j) {
            const float th0 = tau * ((float) j      / (float)slices);
            const float th1 = tau * ((float)(j + 1) / (float)slices);

            const float ct0 = cosf(th0), st0 = sinf(th0);
            const float ct1 = cosf(th1), st1 = sinf(th1);

            // Unit-sphere positions == smooth normals
            vec3_t n00 = {cp0*ct0, cp0*st0, sp0};
            vec3_t n10 = {cp1*ct0, cp1*st0, sp1};
            vec3_t n01 = {cp0*ct1, cp0*st1, sp0};
            vec3_t n11 = {cp1*ct1, cp1*st1, sp1};

            emit_quad(n00, n10, n01, n11);
        }
    }
}

void draw_sphere_wireframe(Queue& ctx, vec3_t center, float radius, uint32_t color, int stacks, int slices) {
    ASSERT(stacks >= 2 && slices >= 3);
    const float pi  = 3.14159265358979323846f;
    const float tau = 2.0f * pi;

    // Latitude rings
    for (int i = 1; i < stacks; ++i) {
        const float phi = pi * ((float)i / (float)stacks) - pi * 0.5f;
        const float cp  = cosf(phi), sp = sinf(phi);
        for (int j = 0; j < slices; ++j) {
            const float th0 = tau * ((float) j      / (float)slices);
            const float th1 = tau * ((float)(j + 1) / (float)slices);
            vec3_t a = vec3_add(center, vec3_mul1({cp*cosf(th0), cp*sinf(th0), sp}, radius));
            vec3_t b = vec3_add(center, vec3_mul1({cp*cosf(th1), cp*sinf(th1), sp}, radius));
            draw_line(ctx, a, b, color);
        }
    }

    // Longitude arcs
    for (int j = 0; j < slices; ++j) {
        const float th = tau * ((float)j / (float)slices);
        const float ct = cosf(th), st = sinf(th);
        for (int i = 0; i < stacks; ++i) {
            const float phi0 = pi * ((float) i      / (float)stacks) - pi * 0.5f;
            const float phi1 = pi * ((float)(i + 1) / (float)stacks) - pi * 0.5f;
            vec3_t a = vec3_add(center, vec3_mul1({cosf(phi0)*ct, cosf(phi0)*st, sinf(phi0)}, radius));
            vec3_t b = vec3_add(center, vec3_mul1({cosf(phi1)*ct, cosf(phi1)*st, sinf(phi1)}, radius));
            draw_line(ctx, a, b, color);
        }
    }
}

void draw_cylinder(Queue& ctx, vec3_t from, vec3_t to, float radius, uint32_t color, uint32_t picking_idx, int segments) {
    ASSERT(segments >= 3);
    const float tau = 6.28318530717958647692f;

    vec3_t axis = vec3_normalize(vec3_sub(to, from));
    vec3_t u, v;
    build_tangent_frame(axis, u, v);

    for (int i = 0; i < segments; ++i) {
        const float th0 = tau * ((float) i      / (float)segments);
        const float th1 = tau * ((float)(i + 1) / (float)segments);

        const float c0 = cosf(th0), s0 = sinf(th0);
        const float c1 = cosf(th1), s1 = sinf(th1);

        // Smooth outward normals for the side: radial direction only (no axis component)
        vec3_t sn0 = vec3_add(vec3_mul1(u, c0), vec3_mul1(v, s0));
        vec3_t sn1 = vec3_add(vec3_mul1(u, c1), vec3_mul1(v, s1));

        vec3_t a0 = vec3_add(from, vec3_mul1(sn0, radius));
        vec3_t a1 = vec3_add(from, vec3_mul1(sn1, radius));
        vec3_t b0 = vec3_add(to,   vec3_mul1(sn0, radius));
        vec3_t b1 = vec3_add(to,   vec3_mul1(sn1, radius));

        // Side quad – smooth normals interpolated across the quad
        {
            const Index base = (Index)md_array_size(ctx.vertices);
            Vertex sv[4] = {
                {a0, color, sn0, picking_idx},
                {a1, color, sn1, picking_idx},
                {b0, color, sn0, picking_idx},
                {b1, color, sn1, picking_idx},
            };
            md_array_push_array(ctx.vertices, sv, 4, ctx.arena);
            md_array_push(ctx.indices, base + 0, ctx.arena);
            md_array_push(ctx.indices, base + 3, ctx.arena);
            md_array_push(ctx.indices, base + 2, ctx.arena);
            md_array_push(ctx.indices, base + 0, ctx.arena);
            md_array_push(ctx.indices, base + 1, ctx.arena);
            md_array_push(ctx.indices, base + 3, ctx.arena);
            append_draw_command(ctx, 6, GL_TRIANGLES);
        }

        // Bottom cap – flat normal pointing away from 'to'
        {
            vec3_t cap_n = vec3_mul1(axis, -1.f);
            const Index base = (Index)md_array_size(ctx.vertices);
            Vertex cv[3] = {
                {from, color, cap_n, picking_idx},
                {a1,   color, cap_n, picking_idx},
                {a0,   color, cap_n, picking_idx},
            };
            md_array_push_array(ctx.vertices, cv, 3, ctx.arena);
            md_array_push(ctx.indices, base + 0, ctx.arena);
            md_array_push(ctx.indices, base + 1, ctx.arena);
            md_array_push(ctx.indices, base + 2, ctx.arena);
            append_draw_command(ctx, 3, GL_TRIANGLES);
        }

        // Top cap – flat normal pointing away from 'from'
        {
            const Index base = (Index)md_array_size(ctx.vertices);
            Vertex cv[3] = {
                {to, color, axis, picking_idx},
                {b0, color, axis, picking_idx},
                {b1, color, axis, picking_idx},
            };
            md_array_push_array(ctx.vertices, cv, 3, ctx.arena);
            md_array_push(ctx.indices, base + 0, ctx.arena);
            md_array_push(ctx.indices, base + 1, ctx.arena);
            md_array_push(ctx.indices, base + 2, ctx.arena);
            append_draw_command(ctx, 3, GL_TRIANGLES);
        }
    }
}

void draw_cylinder_wireframe(Queue& ctx, vec3_t from, vec3_t to, float radius, uint32_t color, int segments) {
    ASSERT(segments >= 3);
    const float tau = 6.28318530717958647692f;

    vec3_t axis = vec3_normalize(vec3_sub(to, from));
    vec3_t u, v;
    build_tangent_frame(axis, u, v);

    for (int i = 0; i < segments; ++i) {
        const float th0 = tau * ((float) i      / (float)segments);
        const float th1 = tau * ((float)(i + 1) / (float)segments);

        vec3_t r0 = vec3_add(vec3_mul1(u, radius * cosf(th0)), vec3_mul1(v, radius * sinf(th0)));
        vec3_t r1 = vec3_add(vec3_mul1(u, radius * cosf(th1)), vec3_mul1(v, radius * sinf(th1)));

        // Bottom ring
        draw_line(ctx, vec3_add(from, r0), vec3_add(from, r1), color);
        // Top ring
        draw_line(ctx, vec3_add(to, r0), vec3_add(to, r1), color);
        // Side edge
        draw_line(ctx, vec3_add(from, r0), vec3_add(to, r0), color);
    }
}

void draw_cone(Queue& ctx, vec3_t base, vec3_t tip, float radius, uint32_t color, uint32_t picking_idx, int segments) {
    ASSERT(segments >= 3);
    const float tau = 6.28318530717958647692f;

    vec3_t axis_vec = vec3_sub(tip, base);
    float height    = vec3_length(axis_vec);
    vec3_t axis     = vec3_mul1(axis_vec, 1.f / height);
    vec3_t u, v;
    build_tangent_frame(axis, u, v);

    // Slant length and trig for cone normal
    const float slant = sqrtf(height * height + radius * radius);
    const float cos_a = height / slant;   // cos of half-angle
    const float sin_a = radius / slant;   // sin of half-angle

    vec3_t base_cap_n = vec3_mul1(axis, -1.f);

    for (int i = 0; i < segments; ++i) {
        const float th0 = tau * ((float) i      / (float)segments);
        const float th1 = tau * ((float)(i + 1) / (float)segments);

        const float c0 = cosf(th0), s0 = sinf(th0);
        const float c1 = cosf(th1), s1 = sinf(th1);

        // Radial direction at each angle
        vec3_t rad0 = vec3_add(vec3_mul1(u, c0), vec3_mul1(v, s0));
        vec3_t rad1 = vec3_add(vec3_mul1(u, c1), vec3_mul1(v, s1));

        // Smooth outward lateral normals: tilt radial inward by the cone half-angle
        vec3_t ln0 = vec3_sub(vec3_mul1(rad0, cos_a), vec3_mul1(axis, sin_a));
        vec3_t ln1 = vec3_sub(vec3_mul1(rad1, cos_a), vec3_mul1(axis, sin_a));

        vec3_t a0 = vec3_add(base, vec3_mul1(rad0, radius));
        vec3_t a1 = vec3_add(base, vec3_mul1(rad1, radius));

        // Lateral face – smooth normals at base ring; use average at apex for continuity
        vec3_t apex_n = vec3_normalize(vec3_add(ln0, ln1));
        {
            const Index base_idx = (Index)md_array_size(ctx.vertices);
            Vertex sv[3] = {
                {a0,  color, ln0,    picking_idx},
                {tip, color, apex_n, picking_idx},
                {a1,  color, ln1,    picking_idx},
            };
            md_array_push_array(ctx.vertices, sv, 3, ctx.arena);
            md_array_push(ctx.indices, base_idx + 0, ctx.arena);
            md_array_push(ctx.indices, base_idx + 2, ctx.arena);
            md_array_push(ctx.indices, base_idx + 1, ctx.arena);
            append_draw_command(ctx, 3, GL_TRIANGLES);
        }

        // Base cap – flat normal pointing away from tip
        {
            const Index base_idx = (Index)md_array_size(ctx.vertices);
            Vertex cv[3] = {
                {base, color, base_cap_n, picking_idx},
                {a0,   color, base_cap_n, picking_idx},
                {a1,   color, base_cap_n, picking_idx},
            };
            md_array_push_array(ctx.vertices, cv, 3, ctx.arena);
            md_array_push(ctx.indices, base_idx + 0, ctx.arena);
            md_array_push(ctx.indices, base_idx + 2, ctx.arena);
            md_array_push(ctx.indices, base_idx + 1, ctx.arena);
            append_draw_command(ctx, 3, GL_TRIANGLES);
        }
    }
}

void draw_cone_wireframe(Queue& ctx, vec3_t base, vec3_t tip, float radius, uint32_t color, int segments) {
    ASSERT(segments >= 3);
    const float tau = 6.28318530717958647692f;

    vec3_t axis = vec3_normalize(vec3_sub(tip, base));
    vec3_t u, v;
    build_tangent_frame(axis, u, v);

    for (int i = 0; i < segments; ++i) {
        const float th0 = tau * ((float) i      / (float)segments);
        const float th1 = tau * ((float)(i + 1) / (float)segments);

        vec3_t r0 = vec3_add(vec3_mul1(u, radius * cosf(th0)), vec3_mul1(v, radius * sinf(th0)));
        vec3_t r1 = vec3_add(vec3_mul1(u, radius * cosf(th1)), vec3_mul1(v, radius * sinf(th1)));

        // Base ring
        draw_line(ctx, vec3_add(base, r0), vec3_add(base, r1), color);
        // Lateral edge
        draw_line(ctx, vec3_add(base, r0), tip, color);
    }
}

void draw_box(Queue& ctx, vec3_t mn, vec3_t mx, uint32_t color, uint32_t picking_idx) {
    const vec3_t p000 = {mn.x, mn.y, mn.z};
    const vec3_t p100 = {mx.x, mn.y, mn.z};
    const vec3_t p010 = {mn.x, mx.y, mn.z};
    const vec3_t p110 = {mx.x, mx.y, mn.z};
    const vec3_t p001 = {mn.x, mn.y, mx.z};
    const vec3_t p101 = {mx.x, mn.y, mx.z};
    const vec3_t p011 = {mn.x, mx.y, mx.z};
    const vec3_t p111 = {mx.x, mx.y, mx.z};

    // -Z face
    draw_triangle(ctx, p000, p110, p100, color, picking_idx);
    draw_triangle(ctx, p000, p010, p110, color, picking_idx);
    // +Z face
    draw_triangle(ctx, p001, p101, p111, color, picking_idx);
    draw_triangle(ctx, p001, p111, p011, color, picking_idx);
    // -Y face
    draw_triangle(ctx, p000, p100, p101, color, picking_idx);
    draw_triangle(ctx, p000, p101, p001, color, picking_idx);
    // +Y face
    draw_triangle(ctx, p010, p111, p110, color, picking_idx);
    draw_triangle(ctx, p010, p011, p111, color, picking_idx);
    // -X face
    draw_triangle(ctx, p000, p001, p011, color, picking_idx);
    draw_triangle(ctx, p000, p011, p010, color, picking_idx);
    // +X face
    draw_triangle(ctx, p100, p110, p111, color, picking_idx);
    draw_triangle(ctx, p100, p111, p101, color, picking_idx);
}

}  // namespace immediate
