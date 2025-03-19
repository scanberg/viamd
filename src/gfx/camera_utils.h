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
    TrackballFlags_PanEnabled           = 0x1,
    TrackballFlags_RotateEnabled        = 0x2,
    TrackballFlags_DollyEnabled         = 0x4,
    TrackballFlags_PanReturnsTrue       = 0x10,
    TrackballFlags_RotateReturnsTrue    = 0x20,
    TrackballFlags_DollyReturnsTrue     = 0x40,
    
    TrackballFlags_EnableAllInteractions     = TrackballFlags_PanEnabled | TrackballFlags_RotateEnabled | TrackballFlags_DollyEnabled,
    TrackballFlags_AnyInteractionReturnsTrue = TrackballFlags_PanReturnsTrue | TrackballFlags_RotateReturnsTrue | TrackballFlags_DollyReturnsTrue
};

typedef uint32_t TrackballFlags;

void camera_trackball(Camera* camera, vec2_t prev_ndc, vec2_t curr_ndc);
void camera_move(Camera* camera, vec3_t vec);

void camera_interpolate_look_at(vec3_t* out_pos, quat_t* out_ori, float* out_dist, vec3_t in_pos[2], quat_t in_ori[2], float in_dist[2], double t);

mat4_t camera_world_to_view_matrix(const Camera& camera);
mat4_t camera_view_to_world_matrix(const Camera& camera);

mat4_t camera_perspective_projection_matrix(const Camera& camera, float aspect_ratio);
mat4_t camera_inverse_perspective_projection_matrix(const Camera& camera, float aspect_ratio);

mat4_t camera_perspective_projection_matrix(const Camera& camera, int width, int height, float texel_offset_x, float texel_offset_y);
mat4_t camera_inverse_perspective_projection_matrix(const Camera& camera, int width, int height, float texel_offset_x, float texel_offset_y);

mat4_t camera_orthographic_projection_matrix(float left, float right, float bottom, float top);
mat4_t camera_inverse_orthographic_projection_matrix(float left, float right, float bottom, float top);

mat4_t camera_orthographic_projection_matrix(float left, float right, float bottom, float top, float near, float far);
mat4_t camera_inverse_orthographic_projection_matrix(float left, float right, float bottom, float top, float near, float far);

// @TODO: Fix the name to something more descriptive. This modifies the position, orientation and distance using a trackball modality
bool camera_controller_trackball(vec3_t* position, quat_t* orientation, float* distance, TrackballControllerInput input, TrackballControllerParam param = {}, TrackballFlags flags = -1);

//void camera_controller_fps(Camera* camera, const FpsControllerState& state);

void camera_compute_optimal_view(vec3_t* out_position, quat_t* out_orientation, float* out_distance, mat3_t in_basis, vec3_t in_min_ext, vec3_t in_max_ext, float in_distance_scale = 3.0f);

// Lazy stupid procedure on top of interpolate_look_at
void camera_animate(Camera* camera, quat_t target_ori, vec3_t target_pos, float target_dist, double dt, double target_factor = 0.12f);