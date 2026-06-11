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

struct Encoder {
    void* impl = nullptr;
};

struct LightingDesc;

struct RenderParams {
    mat4_t view = mat4_ident();
    mat4_t proj = mat4_ident();
    uint32_t picking_base_idx = 0;
    bool shaded = false;
    const LightingDesc* lighting = nullptr;
};

struct GfxContext;
struct Queue;

struct Scope {
    Scope(const char* label = nullptr);
    ~Scope();

    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;
    Scope(Scope&& other) noexcept;
    Scope& operator=(Scope&& other) noexcept;

    Encoder& encoder();

private:
    Encoder enc = {};
    bool active = false;

    friend void queue_submit(Queue* q, Scope& s);
};

void initialize();
void shutdown();

GfxContext* context_create();
void context_destroy(GfxContext* ctx);
Queue* queue_create(GfxContext* ctx, const char* label = nullptr);
void queue_destroy(Queue* queue);
void queue_reset(Queue* queue);
void queue_submit(Queue* queue, Scope& scope);
void render(Queue* queue, const RenderParams& params);

// Encoder-oriented API (preferred)
Encoder encoder_begin(const char* label = nullptr);
void encoder_end(Encoder& e);
void encoder_discard(Encoder& e);
void encoder_set_model(Encoder& e, const mat4_t& model_mat);
void encoder_set_picking_base_idx(Encoder& e, uint32_t base_idx);
void encoder_submit(Encoder& e, const mat4_t& view, const mat4_t& proj, uint32_t picking_base_idx = 0);
void encoder_submit_shaded(Encoder& e, const mat4_t& view, const mat4_t& proj, const LightingDesc& lighting, uint32_t picking_base_idx = 0);

void point(Encoder& e, vec3_t pos, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);
void line(Encoder& e, vec3_t from, vec3_t to, uint32_t color = DEFAULT_COLOR);
void triangle(Encoder& e, vec3_t v0, vec3_t v1, vec3_t v2, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);
void point(Scope& s, vec3_t pos, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);
void line(Scope& s, vec3_t from, vec3_t to, uint32_t color = DEFAULT_COLOR);
void triangle(Scope& s, vec3_t v0, vec3_t v1, vec3_t v2, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);

// Lighting descriptor for render_shaded().
// All radiance values are in linear light units.
// light_dir is a world-space direction pointing *towards* the light source.
struct LightingDesc {
    vec3_t light_dir    = {0.57735026f, 0.57735026f, 0.57735026f};  // normalized world-space direction towards light
    vec3_t dir_radiance = {10.f, 10.f, 10.f};                       // directional light colour * intensity
    vec3_t env_radiance = {1.f,  1.f,  1.f};                        // ambient/environment radiance
    float  roughness    = 0.4f;
    float  F0           = 0.04f;                                     // Fresnel reflectance at normal incidence
};

// Like render(), but triangle draw commands are shaded using the provided lighting descriptor.
// Lines and points are drawn unlit. The same GBuffer MRT layout (locations 0-3) is written.
// 3D Primitives
void draw_point(Encoder& e, vec3_t pos, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);
void draw_line(Encoder& e, vec3_t from, vec3_t to, uint32_t color = DEFAULT_COLOR);
void draw_triangle(Encoder& e, vec3_t v0, vec3_t v1, vec3_t v2, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);

// Sphere centered at 'center' with given radius.
void draw_sphere(Encoder& e, vec3_t center, float radius, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF, int stacks = 12, int slices = 16);
void draw_sphere_wireframe(Encoder& e, vec3_t center, float radius, uint32_t color = DEFAULT_COLOR, int stacks = 12, int slices = 16);

// Cylinder from 'from' to 'to' with given radius. Cap geometry is included for the filled variant.
void draw_cylinder(Encoder& e, vec3_t from, vec3_t to, float radius, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF, int segments = 16);
void draw_cylinder_wireframe(Encoder& e, vec3_t from, vec3_t to, float radius, uint32_t color = DEFAULT_COLOR, int segments = 16);

// Cone with apex at 'tip' and base centered at 'base'.
void draw_cone(Encoder& e, vec3_t base, vec3_t tip, float radius, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF, int segments = 16);
void draw_cone_wireframe(Encoder& e, vec3_t base, vec3_t tip, float radius, uint32_t color = DEFAULT_COLOR, int segments = 16);

// Solid box given axis-aligned min/max corners.
void draw_box(Encoder& e, vec3_t min_box, vec3_t max_box, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);

// Batch
void draw_points_v(Encoder& e, const Vertex verts[], size_t count, vec4_t color_mult = {1,1,1,1});
void draw_lines_v(Encoder& e, const Vertex verts[], size_t count, vec4_t color_mult = {1,1,1,1});
void draw_triangles_v(Encoder& e, const Vertex verts[], size_t count, vec4_t color_mult = {1,1,1,1});

/*
         __________________
        /        ^(v)     /
       /        /        /
      /     (c). -----> /
     /              (u)/
    /_________________/

        Draws a plane given a center point and two support vectors.
*/
void draw_plane(Encoder& e, vec3_t center, vec3_t plane_u, vec3_t plane_v, uint32_t color = DEFAULT_COLOR, uint32_t picking_idx = 0xFFFFFFFF);
void draw_plane_wireframe(Encoder& e, vec3_t center, vec3_t plane_u, vec3_t plane_v, uint32_t color = DEFAULT_COLOR, int segments_u = 4, int segments_v = 4);

// Composits
void draw_box_wireframe(Encoder& e, vec3_t min_box, vec3_t max_box, uint32_t color = DEFAULT_COLOR);
void draw_box_wireframe(Encoder& e, vec3_t min_box, vec3_t max_box, vec4_t color);
void draw_box_wireframe(Encoder& e, vec3_t min_box, vec3_t max_box, mat4_t model_matrix, uint32_t color = DEFAULT_COLOR);
void draw_box_wireframe(Encoder& e, vec3_t min_box, vec3_t max_box, mat4_t model_matrix, vec4_t color);

void draw_basis(Encoder& e, mat4_t basis, float scale = 1.f, uint32_t x_color = COLOR_RED, uint32_t y_color = COLOR_GREEN, uint32_t z_color = COLOR_BLUE);
void draw_basis(Encoder& e, mat4_t basis, float scale = 1.f, vec4_t x_color = {1,0,0,1}, vec4_t y_color = {0,1,0,1}, vec4_t z_color = {0,0,1,1});

// Advanced
// void draw_sphere_glyph(const float pos[3], const float radius, const uint32 color);
// void draw_capsule(const float v0[3], const float v1[3], const uint32 color);

}  // namespace immediate
