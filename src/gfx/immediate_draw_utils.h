#pragma once

#include <core/md_vec_math.h>

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
    vec3_t   coord;
    uint32_t color;
	vec3_t   normal = {0,0,1};
    uint32_t picking_idx = 0xFFFFFFFF;
};

// All radiance values are in linear light units.
// light_dir is a world-space direction pointing *towards* the light source.
struct LightingDesc {
    vec3_t light_dir    = {0.57735026f, 0.57735026f, 0.57735026f};  // normalized world-space direction towards light
    vec3_t dir_radiance = {10.f, 10.f, 10.f};                       // directional light colour * intensity
    vec3_t env_radiance = {1.f,  1.f,  1.f};                        // ambient/environment radiance
    float  roughness    = 0.4f;
    float  F0           = 0.04f;                                    // Fresnel reflectance at normal incidence
};

struct RenderParams {
    mat4_t view = mat4_ident();
    mat4_t proj = mat4_ident();
    uint32_t picking_base_idx = 0;
    bool shaded = false;
    const LightingDesc* lighting = nullptr;
};

struct Queue;
struct Command;

struct Scope {
    Scope(const char* label = nullptr);
    ~Scope();

    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;
    Scope(Scope&& other) noexcept;
    Scope& operator=(Scope&& other) noexcept;

    struct Encoder;
    Encoder* encoder = nullptr;
};

void initialize();
void shutdown();

Queue* queue_create(const char* label = nullptr);
void queue_destroy(Queue* queue);
void queue_reset(Queue* queue);
void queue_submit(Queue* queue, Scope& scope);

void render(Queue* queue, const RenderParams& params);

void set_model(Scope& s, const mat4_t& model_mat);
void set_picking_base_idx(Scope& s, uint32_t base_idx);

void submit(Queue* queue, Scope& scope);

void point(Scope& s, vec3_t pos, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);
void line(Scope& s, vec3_t from, vec3_t to, uint32_t color = DEFAULT_COLOR);
void triangle(Scope& s, vec3_t v0, vec3_t v1, vec3_t v2, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);
void triangle_wireframe(Scope& s, vec3_t v0, vec3_t v1, vec3_t v2, uint32_t color = DEFAULT_COLOR);

// Sphere centered at 'center' with given radius.
void sphere(Scope& s, vec3_t center, float radius, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF, int stacks = 12, int slices = 16);
void sphere_wireframe(Scope& s, vec3_t center, float radius, uint32_t color = DEFAULT_COLOR, int stacks = 12, int slices = 16);

// Cylinder from 'from' to 'to' with given radius. Cap geometry is included for the filled variant.
void cylinder(Scope& s, vec3_t from, vec3_t to, float radius, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF, int segments = 16);
void cylinder_wireframe(Scope& s, vec3_t from, vec3_t to, float radius, uint32_t color = DEFAULT_COLOR, int segments = 16);

// Cone with apex at 'tip' and base centered at 'base'.
void cone(Scope& s, vec3_t base, vec3_t tip, float radius, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF, int segments = 16);
void cone_wireframe(Scope& s, vec3_t base, vec3_t tip, float radius, uint32_t color = DEFAULT_COLOR, int segments = 16);

void capsule(Scope& s, vec3_t from, vec3_t to, float radius, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF, int segments = 16);
void capsule_wireframe(Scope& s, vec3_t from, vec3_t to, float radius, uint32_t color = DEFAULT_COLOR, int segments = 16);

// Solid box given axis-aligned min/max corners.
void box(Scope& s, vec3_t min_box, vec3_t max_box, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);

// Batch
void points(Scope& s, const Vertex verts[], size_t count, vec4_t color_mult = {1,1,1,1});
void lines(Scope& s, const Vertex verts[], size_t count, vec4_t color_mult = {1,1,1,1});
void triangles(Scope& s, const Vertex verts[], size_t count, vec4_t color_mult = {1,1,1,1});
void triangles_wireframe(Scope& s, const Vertex verts[], size_t count, vec4_t color_mult = { 1,1,1,1 });

/*
         __________________
        /        ^(v)     /
       /        /        /
      /     (c). -----> /
     /              (u)/
    /_________________/

        Draws a plane given a center point and two support vectors.
*/
void plane(Scope& s, vec3_t center, vec3_t plane_u, vec3_t plane_v, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);
void plane_wireframe(Scope& s, vec3_t center, vec3_t plane_u, vec3_t plane_v, uint32_t color = DEFAULT_COLOR, int segments_u = 4, int segments_v = 4);

// Composits
void box_wireframe(Scope& s, vec3_t min_box, vec3_t max_box, uint32_t color = DEFAULT_COLOR);
void box_wireframe(Scope& s, vec3_t min_box, vec3_t max_box, vec4_t color);
void box_wireframe(Scope& s, vec3_t min_box, vec3_t max_box, mat4_t model_matrix, uint32_t color = DEFAULT_COLOR);
void box_wireframe(Scope& s, vec3_t min_box, vec3_t max_box, mat4_t model_matrix, vec4_t color);

void basis(Scope& s, mat4_t basis, float scale = 1.f, uint32_t x_color = COLOR_RED, uint32_t y_color = COLOR_GREEN, uint32_t z_color = COLOR_BLUE);
void basis(Scope& s, mat4_t basis, float scale = 1.f, vec4_t x_color = {1,0,0,1}, vec4_t y_color = {0,1,0,1}, vec4_t z_color = {0,0,1,1});

}  // namespace immediate
