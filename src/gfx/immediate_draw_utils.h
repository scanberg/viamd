#pragma once

#include <core/types.h>

namespace immediate {

void initialize();
void shutdown();

void set_view_matrix(const mat4& model_view_mat);
void set_proj_matrix(const mat4& proj_mat);
void flush();

constexpr unsigned char COLOR_WHITE[4] = { 255,255,255,255 };
constexpr unsigned char COLOR_BLACK[4] = { 0,0,0,255 };
constexpr unsigned char COLOR_RED[4] = { 255,0,0,255 };
constexpr unsigned char COLOR_GREEN[4] = { 0,255,0,255 };
constexpr unsigned char COLOR_BLUE[4] = { 0,0,255,255 };
constexpr unsigned char COLOR_YELLOW[4] = { 255,255,0,255 };
constexpr unsigned char COLOR_MAGENTA[4] = { 255,0,255,255 };
constexpr unsigned char COLOR_CYAN[4] = { 0,255,255,255 };

// Primitives
void draw_point(const float pos[3], const unsigned char color[4] = COLOR_WHITE);
inline void draw_point(vec3 pos, const unsigned char color[4] = COLOR_WHITE) {
	draw_point(&pos[0], color);
}

void draw_line(const float from[3], const float to[3], const unsigned char color[4] = COLOR_WHITE);
inline void draw_line(vec3 from, vec3 to, const unsigned char color[4] = COLOR_WHITE) {
	draw_line(&from[0], &to[0], color);
}

void draw_triangle(const float v0[3], const float v1[3], const float v2[3], const unsigned char color[4] = COLOR_WHITE);
inline void draw_triangle(vec3 v0, vec3 v1, vec3 v2, const unsigned char color[4] = COLOR_WHITE) {
	draw_triangle(&v0[0], &v1[0], &v2[0], color);
}

void draw_quad(const float v0[3], const float v1[3], const float v2[3], const float v3[3], const unsigned char color[4] = COLOR_WHITE);
inline void draw_quad(vec3 v0, vec3 v1, vec3 v2, vec3 v3, const unsigned char color[4] = COLOR_WHITE) {
	draw_quad(&v0[0], &v1[0], &v2[0], &v3[0], color);
}

// Composits
void draw_wired_sphere(const float pos[3], const float radius, const unsigned char color[4] = COLOR_WHITE);
void draw_axis_aligned_box(const float min_box[3], const float max_box[3], const unsigned char color[4] = COLOR_WHITE);

// Advanced
void draw_sphere_glyph(const float pos[3], const float radius, const unsigned char color[4]);
void draw_capsule(const float v0[3], const float v1[3], const unsigned char color[4]);

}