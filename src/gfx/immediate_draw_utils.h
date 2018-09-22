#pragma once

#include <core/types.h>
#include <core/vector_types.h>

namespace immediate {

void initialize();
void shutdown();

struct Material {
    vec4 color_alpha = {1, 1, 1, 1};
    vec3 f0 = {0.04f, 0.04f, 0.04f};
    float smoothness = 0;
    vec2 uv_scale = {1, 1};
    uint32 texture_id = 0;
    float _pad = 0;
};

void set_view_matrix(const mat4& model_view_mat);
void set_proj_matrix(const mat4& proj_mat);
void set_material(const Material& mat);

void flush();

const vec4 COLOR_WHITE = {1, 1, 1, 1};
const vec4 COLOR_BLACK = {0, 0, 0, 1};
const vec4 COLOR_RED = {1, 0, 0, 1};
const vec4 COLOR_GREEN = {0, 1, 0, 1};
const vec4 COLOR_BLUE = {0, 0, 1, 1};
const vec4 COLOR_YELLOW = {1, 1, 0, 1};
const vec4 COLOR_MAGENTA = {1, 0, 1, 1};
const vec4 COLOR_CYAN = {0, 1, 1, 1};

const Material MATERIAL_ROUGH_WHITE = {COLOR_WHITE};
const Material MATERIAL_ROUGH_BLACK = {COLOR_BLACK};
const Material MATERIAL_ROUGH_RED = {COLOR_RED};
const Material MATERIAL_ROUGH_GREEN = {COLOR_GREEN};
const Material MATERIAL_ROUGH_BLUE = {COLOR_BLUE};
const Material MATERIAL_ROUGH_YELLOW = {COLOR_YELLOW};
const Material MATERIAL_ROUGH_MAGENTA = {COLOR_MAGENTA};
const Material MATERIAL_ROUGH_CYAN = {COLOR_CYAN};

const Material MATERIAL_GLOSSY_WHITE = {COLOR_WHITE, {0.04f, 0.04f, 0.04f}, 0.5};
const Material MATERIAL_GLOSSY_BLACK = {COLOR_BLACK, {0.04f, 0.04f, 0.04f}, 0.5};
const Material MATERIAL_GLOSSY_RED = {COLOR_RED, {0.04f, 0.04f, 0.04f}, 0.5};
const Material MATERIAL_GLOSSY_GREEN = {COLOR_GREEN, {0.04f, 0.04f, 0.04f}, 0.5};
const Material MATERIAL_GLOSSY_BLUE = {COLOR_BLUE, {0.04f, 0.04f, 0.04f}, 0.5};
const Material MATERIAL_GLOSSY_YELLOW = {COLOR_YELLOW, {0.04f, 0.04f, 0.04f}, 0.5};
const Material MATERIAL_GLOSSY_MAGENTA = {COLOR_MAGENTA, {0.04f, 0.04f, 0.04f}, 0.5};
const Material MATERIAL_GLOSSY_CYAN = {COLOR_CYAN, {0.04f, 0.04f, 0.04f}, 0.5};

// Primitives

void draw_point(const vec3& pos);
void draw_line(const vec3& from, const vec3& to);
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
void draw_basis(const mat4& basis, const float scale = 1.f);

// Advanced
// void draw_sphere_glyph(const float pos[3], const float radius, const uint32 color);
// void draw_capsule(const float v0[3], const float v1[3], const uint32 color);

}  // namespace immediate
