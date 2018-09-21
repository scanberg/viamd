#pragma once

#include <core/types.h>
#include <core/vector_types.h>

namespace immediate {

void initialize();
void shutdown();

struct Material {
    vec3 base_color = {1, 1, 1};
    float alpha = 1;
    vec3 f0 = {0.04f, 0.04f, 0.04f};
    float smoothness = 0;
    vec2 uv_scale = {1, 1};
    uint32 texture_id = 0;
    float _pad;
};

void set_view_matrix(const mat4& model_view_mat);
void set_proj_matrix(const mat4& proj_mat);
void set_material(const Material& mat);

void flush();

constexpr vec3 COLOR_WHITE = {1, 1, 1};
constexpr vec3 COLOR_BLACK = {0, 0, 0};
constexpr vec3 COLOR_RED = {1, 0, 0};
constexpr vec3 COLOR_GREEN = {0, 1, 0};
constexpr vec3 COLOR_BLUE = {0, 0, 1};
constexpr vec3 COLOR_YELLOW = {1, 1, 0};
constexpr vec3 COLOR_MAGENTA = {1, 0, 1};
constexpr vec3 COLOR_CYAN = {0, 1, 1};

constexpr Material MATERIAL_ROUGH_WHITE = {COLOR_WHITE};
constexpr Material MATERIAL_ROUGH_BLACK = {COLOR_BLACK};
constexpr Material MATERIAL_ROUGH_RED = {COLOR_RED};
constexpr Material MATERIAL_ROUGH_GREEN = {COLOR_GREEN};
constexpr Material MATERIAL_ROUGH_BLUE = {COLOR_BLUE};
constexpr Material MATERIAL_ROUGH_YELLOW = {COLOR_YELLOW};
constexpr Material MATERIAL_ROUGH_MAGENTA = {COLOR_MAGENTA};
constexpr Material MATERIAL_ROUGH_CYAN = {COLOR_CYAN};

constexpr Material MATERIAL_GLOSSY_WHITE = {COLOR_WHITE, 1, {0.04f, 0.04f, 0.04f}, 0.5};
constexpr Material MATERIAL_GLOSSY_BLACK = {COLOR_BLACK, 1, {0.04f, 0.04f, 0.04f}, 0.5};
constexpr Material MATERIAL_GLOSSY_RED = {COLOR_RED, 1, {0.04f, 0.04f, 0.04f}, 0.5};
constexpr Material MATERIAL_GLOSSY_GREEN = {COLOR_GREEN, 1, {0.04f, 0.04f, 0.04f}, 0.5};
constexpr Material MATERIAL_GLOSSY_BLUE = {COLOR_BLUE, 1, {0.04f, 0.04f, 0.04f}, 0.5};
constexpr Material MATERIAL_GLOSSY_YELLOW = {COLOR_YELLOW, 1, {0.04f, 0.04f, 0.04f}, 0.5};
constexpr Material MATERIAL_GLOSSY_MAGENTA = {COLOR_MAGENTA, 1, {0.04f, 0.04f, 0.04f}, 0.5};
constexpr Material MATERIAL_GLOSSY_CYAN = {COLOR_CYAN, 1, {0.04f, 0.04f, 0.04f}, 0.5};

// Primitives

// 2D
void draw_point(const vec3& pos, const vec3& color = COLOR_WHITE);
void draw_line(const vec3& from, const vec3& to, const vec3& color = COLOR_WHITE);

// 3D
void draw_triangle(const vec3& v0, const vec3& v1, const vec3& v2);

/*
         __________________
        /        ^(v)     /
       /        /        /
      /     (c). -----> /
     /              (u)/
    /_________________/

        Draws a plane given a center point and two support vectors.
*/
void draw_plane(const vec3& center, const vec3& plane_u, const vec3& plane_v);

// Composits
void draw_aabb(const vec3& min_box, const vec3& max_box);
void draw_aabb_lines(const vec3& min_box, const vec3& max_box);
void draw_basis(const mat4& basis, const float scale = 1.f, const vec3& x_axis_color = COLOR_RED, const vec3& y_axis_color = COLOR_GREEN,
                const vec3& z_axis_color = COLOR_BLUE);

// Advanced
// void draw_sphere_glyph(const float pos[3], const float radius, const uint32 color);
// void draw_capsule(const float v0[3], const float v1[3], const uint32 color);

}  // namespace immediate
