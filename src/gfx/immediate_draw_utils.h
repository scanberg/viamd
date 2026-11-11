#pragma once

#include <core/md_vec.h>
#include <stdint.h>

namespace immediate {

constexpr uint32_t COLOR_WHITE = 0xffffffff;
constexpr uint32_t COLOR_BLACK = 0xff000000;
constexpr uint32_t COLOR_RED = 0xff0000ff;
constexpr uint32_t COLOR_GREEN = 0xff00ff00;
constexpr uint32_t COLOR_BLUE = 0xffff0000;
constexpr uint32_t COLOR_YELLOW = 0xff00ffff;
constexpr uint32_t COLOR_MAGENTA = 0xffff00ff;
constexpr uint32_t COLOR_CYAN = 0xffffff00;
constexpr uint32_t COLOR_GRAY = 0xff888888;

constexpr uint32_t DEFAULT_COLOR = COLOR_BLACK;

struct Vertex {
    md_f32x3_t   coord;
    uint32_t color;
	md_f32x3_t   normal = {0,0,1};
    uint32_t picking_idx = 0xFFFFFFFF;
};

struct RenderParams {
    md_f32x4x4_t view = MD_MAT4_IDENT_INIT;
    md_f32x4x4_t proj = MD_MAT4_IDENT_INIT;
};

struct Queue;
struct Command;

// Scope is optional semantic sugar around a queue-backed recording span.
// It does not own command storage; it only marks begin/end and resets model state.
struct Scope {
    Scope(Queue* queue, const char* label = nullptr);
    ~Scope();

    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;
    Scope(Scope&& other) noexcept;
    Scope& operator=(Scope&& other) noexcept;

	// Support implicit cast of Scope to queue to allow for submit to scope syntax: e.g. `immediate::line(scope, ...)`
	operator Queue*() { return queue; }

    Queue* queue = nullptr;
};

void initialize();
void shutdown();

Queue* queue_create(const char* label = nullptr);
void queue_destroy(Queue* queue);
void queue_reset(Queue* queue);

void render(Queue* queue, const RenderParams& params);

// Direct queue API: emit commands straight into queue storage without creating a Scope.
void set_model(Queue* queue, const md_f32x4x4_t& model_mat);
void set_picking_base_idx(Queue* queue, uint32_t base_idx);

void point(Queue* queue, md_f32x3_t pos, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);
void line(Queue* queue, md_f32x3_t from, md_f32x3_t to, uint32_t color = DEFAULT_COLOR);
void triangle(Queue* queue, md_f32x3_t v0, md_f32x3_t v1, md_f32x3_t v2, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);
void triangle_wireframe(Queue* queue, md_f32x3_t v0, md_f32x3_t v1, md_f32x3_t v2, uint32_t color = DEFAULT_COLOR);

void sphere(Queue* queue, md_f32x3_t center, float radius, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF, int stacks = 12, int slices = 16);
void sphere_wireframe(Queue* queue, md_f32x3_t center, float radius, uint32_t color = DEFAULT_COLOR, int stacks = 12, int slices = 16);

void cylinder(Queue* queue, md_f32x3_t from, md_f32x3_t to, float radius, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF, int segments = 16);
void cylinder_wireframe(Queue* queue, md_f32x3_t from, md_f32x3_t to, float radius, uint32_t color = DEFAULT_COLOR, int segments = 16);

void cone(Queue* queue, md_f32x3_t base, md_f32x3_t tip, float radius, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF, int segments = 16);
void cone_wireframe(Queue* queue, md_f32x3_t base, md_f32x3_t tip, float radius, uint32_t color = DEFAULT_COLOR, int segments = 16);

// A solid arrow from 'from' to 'to': a shaft and a conical head, all in world units, as one closed surface.
// The head is shortened to the arrow's length if it is longer; its radius is at least the shaft's.
void arrow(Queue* queue, md_f32x3_t from, md_f32x3_t to, float shaft_radius, float head_radius, float head_length, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF, int segments = 32);

// The side of a truncated cone from 'from' to 'to', its radius going linearly from radius_from to radius_to, open at
// both ends (closed by whatever joins it there, e.g. spheres of the same radii)
void frustum(Queue* queue, md_f32x3_t from, md_f32x3_t to, float radius_from, float radius_to, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF, int segments = 16);

void capsule(Queue* queue, md_f32x3_t from, md_f32x3_t to, float radius, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF, int segments = 16);
void capsule_wireframe(Queue* queue, md_f32x3_t from, md_f32x3_t to, float radius, uint32_t color = DEFAULT_COLOR, int segments = 16);

void box(Queue* queue, md_f32x3_t min_box, md_f32x3_t max_box, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);

void points(Queue* queue, const Vertex verts[], size_t count, md_f32x4_t color_mult = {1,1,1,1});
void lines(Queue* queue, const Vertex verts[], size_t count, md_f32x4_t color_mult = {1,1,1,1});
void triangles(Queue* queue, const Vertex verts[], size_t count, md_f32x4_t color_mult = {1,1,1,1});
void triangles_wireframe(Queue* queue, const Vertex verts[], size_t count, md_f32x4_t color_mult = { 1,1,1,1 });

/*
__________________
/        ^(v)     /
/        /        /
/     (c). -----> /
/              (u)/
/_________________/

Draws a plane given a center point and two support vectors.
*/

void plane(Queue* queue, md_f32x3_t center, md_f32x3_t plane_u, md_f32x3_t plane_v, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);
void plane_wireframe(Queue* queue, md_f32x3_t center, md_f32x3_t plane_u, md_f32x3_t plane_v, uint32_t color = DEFAULT_COLOR, int segments_u = 4, int segments_v = 4);

void box_wireframe(Queue* queue, md_f32x3_t min_box, md_f32x3_t max_box, uint32_t color = DEFAULT_COLOR);
void box_wireframe(Queue* queue, md_f32x3_t min_box, md_f32x3_t max_box, md_f32x4_t color);
void box_wireframe(Queue* queue, md_f32x3_t min_box, md_f32x3_t max_box, md_f32x4x4_t model_matrix, uint32_t color = DEFAULT_COLOR);
void box_wireframe(Queue* queue, md_f32x3_t min_box, md_f32x3_t max_box, md_f32x4x4_t model_matrix, md_f32x4_t color);

void basis(Queue* queue, md_f32x4x4_t basis, float scale = 1.f, uint32_t x_color = COLOR_RED, uint32_t y_color = COLOR_GREEN, uint32_t z_color = COLOR_BLUE);
void basis(Queue* queue, md_f32x4x4_t basis, float scale = 1.f, md_f32x4_t x_color = {1,0,0,1}, md_f32x4_t y_color = {0,1,0,1}, md_f32x4_t z_color = {0,0,1,1});

}  // namespace immediate
