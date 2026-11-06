#pragma once

#include <core/md_vec_math.h>
#include "camera.h"

/*
struct FpsControllerState {
    struct {
        bool forward_button = false;
        bool backward_button = false;
        bool left_button = false;
        bool right_button = false;
        vec2 mouse_vel_ndc = {0, 0};
        float delta_time = 0.f;
    } input;

    struct {
        float move_speed = 1.f;
        float rotation_speed = 1.f;
    } params;
};
*/
struct TrackballControllerInput {
    bool rotate_button = false;
    bool pan_button = false;
    bool dolly_button = false;
    float dolly_delta = 0;
    vec2_t mouse_coord_prev = {0, 0};
    vec2_t mouse_coord_curr = {0, 0};
    vec2_t screen_size = {0, 0};
    float fov_y = 0.5f;
    vec2_t ndc_center = {0, 0};     // Where in the view the look-at point is seen (see ViewRegion), which rotation is about
};

struct TrackballControllerParam {
    float pan_scale = 1.0f;
    float pan_exponent = 1.f;
    float dolly_drag_scale = 0.01f;
    float dolly_drag_exponent = 1.1f;
    float dolly_delta_scale = 0.1f;
    float dolly_delta_exponent = 1.1f;
    float min_distance = 1.f;
    float max_distance = 1000.f;
};

/*
struct TrackballControllerState {
    struct {
        bool rotate_button = false;
        bool pan_button = false;
        bool dolly_button = false;
        float dolly_delta = 0;
        vec2 mouse_coord_prev = {0, 0};
        vec2 mouse_coord_curr = {0, 0};
        vec2 screen_size = {0, 0};
        float fov_y = 0.5f;
    } input;

    struct {
        float pan_scale = 1.0f;
        float pan_exponent = 1.f;
        float dolly_drag_scale = 0.01f;
        float dolly_drag_exponent = 1.1f;
        float dolly_delta_scale = 0.1f;
        float dolly_delta_exponent = 1.1f;
        float min_distance = 1.f;
        float max_distance = 1000.f;
    } params;

    // Distance to the rotational focus point
    float distance = 14.f;
};
*/

enum TrackballFlags_ {
	TrackballFlags_None                 = 0x0,
    TrackballFlags_PanEnabled           = 0x1,
    TrackballFlags_RotateEnabled        = 0x2,
    TrackballFlags_DollyEnabled         = 0x4,
    TrackballFlags_PanReturnsTrue       = 0x10,
    TrackballFlags_RotateReturnsTrue    = 0x20,
    TrackballFlags_DollyReturnsTrue     = 0x40,
    
    TrackballFlags_EnableAllInteractions     = TrackballFlags_PanEnabled | TrackballFlags_RotateEnabled | TrackballFlags_DollyEnabled,
    TrackballFlags_AnyInteractionReturnsTrue = TrackballFlags_PanReturnsTrue | TrackballFlags_RotateReturnsTrue | TrackballFlags_DollyReturnsTrue,
};

typedef uint32_t TrackballFlags;

void camera_trackball(Camera* camera, vec2_t prev_ndc, vec2_t curr_ndc);
void camera_move(Camera* camera, vec3_t vec);
vec3_t camera_get_look_at(const ViewTransform& transform);
vec3_t camera_position_from_look_at(const vec3_t& look_at, const quat_t& orientation, float distance);

void camera_interpolate_look_at(vec3_t* out_pos, quat_t* out_ori, float* out_dist, vec3_t in_pos[2], quat_t in_ori[2], float in_dist[2], double t);

mat4_t camera_world_to_view_matrix(const ViewTransform& transform);
mat4_t camera_view_to_world_matrix(const ViewTransform& transform);

mat4_t camera_view_to_clip_matrix_persp(const Camera& camera, float aspect_ratio);
mat4_t camera_clip_to_view_matrix_persp(const Camera& camera, float aspect_ratio);

mat4_t camera_view_to_clip_matrix_persp(const Camera& camera, int width, int height, float texel_offset_x, float texel_offset_y);
mat4_t camera_clip_to_view_matrix_persp(const Camera& camera, int width, int height, float texel_offset_x, float texel_offset_y);

mat4_t camera_view_to_clip_matrix_ortho(float left, float right, float bottom, float top);
mat4_t camera_clip_to_view_matrix_ortho(float left, float right, float bottom, float top);

mat4_t camera_view_to_clip_matrix_ortho(float left, float right, float bottom, float top, float near, float far);
mat4_t camera_clip_to_view_matrix_ortho(float left, float right, float bottom, float top, float near, float far);

// @TODO: Fix the name to something more descriptive. This modifies the position, orientation and distance using a trackball modality
bool camera_controller_trackball(ViewTransform* transform, const TrackballControllerInput& input, const TrackballControllerParam& param = {}, TrackballFlags flags = -1);

//void camera_controller_fps(Camera* camera, const FpsControllerState& state);

ViewTransform compute_optimal_view(const vec3_t& center, const vec3_t& half_ext, const mat3_t& basis = mat3_ident(), float distance_scale = 3.0f);

// A rectangle of a view, such as the part of it that windows leave uncovered: min and max in the view's
// normalized device coordinates (x to the right, y up, -1 to 1 across the view), and the view's width over
// its height. The default is all of a square view.
//
// What is framed in a region is seen at its center by shifting the projection there: the projection matrix
// is premultiplied by mat4_translate(center.x, center.y, 0), which moves the picture over in clip space. The
// camera still looks at what it frames, so it rotates about it and zooms towards it.
struct ViewRegion {
    vec2_t min = {-1, -1};
    vec2_t max = { 1,  1};
    float  aspect_ratio = 1.0f;
};

static inline vec2_t view_region_center(const ViewRegion& region) { return (region.min + region.max) * 0.5f; }

// A good default view of a set of atoms: the whole system, or a structure within it. x, y, z hold all
// num_atoms atoms, as drawn (no periodic treatment); the view frames those listed in indices, or all of them
// when indices is null. cell_A is the unit cell basis (columns a, b, c), or null for none - it is only used
// to judge whether the world axes mean anything. The atoms are framed to fit in region, when the projection
// is shifted to its center (see ViewRegion). See the implementation for the choices this makes.
ViewTransform camera_compute_default_view(const vec3_t* xyz, size_t num_atoms, const int32_t* indices, size_t count, const mat3_t* cell_A, float fov_y, const ViewRegion& region = {});

// The camera distance at which all atoms (those in indices, or all count of them) fit in a view from
// orientation, looking at look_at: in the largest square in the middle of region, with the projection
// shifted to its center (see ViewRegion).
float camera_fit_distance(const vec3_t* xyz, const int32_t* indices, size_t count, vec3_t look_at, quat_t orientation, float fov_y, const ViewRegion& region = {});

// Lazy stupid procedure on top of interpolate_look_at
void camera_animate(ViewTransform* current, const ViewTransform& target, double dt, double target_factor = 0.12f);

// A shift of the projection (see ViewRegion) animated towards target, at the pace of camera_animate
void camera_animate_shift(vec2_t* current, vec2_t target, double dt, double target_factor = 0.12f);