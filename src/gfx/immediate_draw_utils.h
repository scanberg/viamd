#pragma once

#include <core/types.h>
#include <core/vector_types.h>

namespace immediate {

void initialize();
void shutdown();

void set_view_matrix(const mat4& model_view_mat);
void set_proj_matrix(const mat4& proj_mat);
void flush();

constexpr uint32 COLOR_WHITE = 0xffffffff;
constexpr uint32 COLOR_BLACK = 0xff000000;
constexpr uint32 COLOR_RED = 0xff0000ff;
constexpr uint32 COLOR_GREEN = 0xff00ff00;
constexpr uint32 COLOR_BLUE = 0xffff0000;
constexpr uint32 COLOR_YELLOW = 0xff00ffff;
constexpr uint32 COLOR_MAGENTA = 0xffff00ff;
constexpr uint32 COLOR_CYAN = 0xffffff00;

// Primitives
void draw_point(const float pos[3], const uint32 color = COLOR_WHITE);
inline void draw_point(vec3 pos, const uint32 color = COLOR_WHITE) { draw_point(&pos[0], color); }

void draw_line(const float from[3], const float to[3], const uint32 color = COLOR_WHITE);
inline void draw_line(vec3 from, vec3 to, const uint32 color = COLOR_WHITE) { draw_line(&from[0], &to[0], color); }

void draw_triangle(const float v0[3], const float v1[3], const float v2[3], const uint32 color = COLOR_WHITE);
inline void draw_triangle(vec3 v0, vec3 v1, vec3 v2, const uint32 color = COLOR_WHITE) { draw_triangle(&v0[0], &v1[0], &v2[0], color); }

// void draw_quad(const float v0[3], const float v1[3], const float v2[3], const float v3[3], const uint32 color = COLOR_WHITE);
// inline void draw_quad(vec3 v0, vec3 v1, vec3 v2, vec3 v3, const uint32 color = COLOR_WHITE) {
//	draw_quad(&v0[0], &v1[0], &v2[0], &v3[0], color);
//}

// Composits
// void draw_wired_sphere(const float pos[3], const float radius, const uint32 color = COLOR_WHITE);
void draw_aabb(const float min_box[3], const float max_box[3], const uint32 color = COLOR_WHITE, bool filled = false);
inline void draw_aabb(vec3 min_box, vec3 max_box, const uint32 color = COLOR_WHITE, bool filled = false) {
    draw_aabb(&min_box[0], &max_box[0], color, filled);
}

void draw_basis(const mat4& basis, const float scale = 1.f, const uint32 x_axis_color = COLOR_RED, const uint32 y_axis_color = COLOR_GREEN,
                const uint32 z_axis_color = COLOR_BLUE);

// Advanced
// void draw_sphere_glyph(const float pos[3], const float radius, const uint32 color);
// void draw_capsule(const float v0[3], const float v1[3], const uint32 color);

}  // namespace immediate
