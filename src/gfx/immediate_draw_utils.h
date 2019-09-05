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

constexpr uint32 DEFAULT_COLOR = COLOR_WHITE;

// 3D Primitives
void draw_point(const vec3& pos, uint32 color = DEFAULT_COLOR);
void draw_line(const vec3& from, const vec3& to, uint32 color = DEFAULT_COLOR);
void draw_triangle(const vec3& v0, const vec3& v1, const vec3& v2, uint32 color = DEFAULT_COLOR);

/*
         __________________
        /        ^(v)     /
       /        /        /
      /     (c). -----> /
     /              (u)/
    /_________________/

        Draws a plane given a center point and two support vectors.
*/
void draw_plane(const vec3& center, const vec3& plane_u, const vec3& plane_v, uint32 color = DEFAULT_COLOR);
void draw_plane_wireframe(const vec3& center, const vec3& plane_u, const vec3& plane_v, uint32 color = DEFAULT_COLOR, int segments_u = 4, int segments_v = 4);

// Composits
void draw_box_wireframe(const vec3& min_box, const vec3& max_box, uint32 color = DEFAULT_COLOR);
void draw_box_wireframe(const vec3& min_box, const vec3& max_box, const mat4& model_matrix, uint32 color = DEFAULT_COLOR);
void draw_basis(const mat4& basis, const float scale = 1.f, uint32 x_color = COLOR_RED, uint32 y_color = COLOR_GREEN, uint32 z_color = COLOR_BLUE);

// Advanced
// void draw_sphere_glyph(const float pos[3], const float radius, const uint32 color);
// void draw_capsule(const float v0[3], const float v1[3], const uint32 color);

}  // namespace immediate
