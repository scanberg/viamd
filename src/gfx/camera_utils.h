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
// How far in front of the camera a point is, along the viewing direction (negative behind it).
// This is the depth that depth of field focuses at, and for the look-at point it equals the distance.
float camera_depth_of_point(const ViewTransform& transform, vec3_t point);
// Turns the camera to look at 'point' with the eye fixed, keeping its roll. The distance becomes the eye-to-point distance.
// Returns false (leaving the transform unchanged) if the point is on top of the eye.
bool camera_aim_at(ViewTransform* transform, vec3_t point);
vec3_t camera_position_from_look_at(const vec3_t& look_at, const quat_t& orientation, float distance);

void camera_interpolate_look_at(vec3_t* out_pos, quat_t* out_ori, float* out_dist, vec3_t in_pos[2], quat_t in_ori[2], float in_dist[2], double t);

// Pose of the camera at 'time' along a path through keyframes, which must be sorted by time (count > 0).
// The camera holds the first/last keyframe outside their range and moves smoothly (continuous velocity)
// through the ones in between, so it does not stop at each of them. The look-at point follows a
// Catmull-Rom spline, distance and fov a monotone cubic (no overshoot) and the orientation a spherical
// cubic Bezier along the shortest rotation. The start and end of the path ease in and out.
// A keyframe with spin_turns adds whole turns around the look-at point over the segment that ends at it.
// The easing of a key shapes the segment that ends at it. With loop, the path is cyclic: the last key must
// be in the same pose as the first, and the velocity is then continuous across the seam.
// follow_now is where the follow target is now. The look-at point of a key with follow is then kept relative to
// the target, blended with the plain look-at across segments between keys that do and do not follow. Without it,
// follow is ignored.
// With upright, the camera is kept level: its up is turned as close to *upright as it can be about the direction it
// looks in, then by the keys' roll (interpolated). A spin around the camera's up is then a spin around *upright.
void camera_keyframes_evaluate(ViewTransform* out_transform, float* out_fov_y, const CameraKeyframe* keys, size_t count, double time, bool loop = false, const vec3_t* follow_now = nullptr, const vec3_t* key_follow_now = nullptr, const vec3_t* upright = nullptr, quat_t* out_spin = nullptr);

// Turns the camera about the direction it looks in, so that its up is as close to 'up' as it can be, then by 'roll'
// radians (positive leans the camera's up to the left). Where it looks from and at stay. Looking along 'up', the
// camera's right is kept instead.
void camera_level(ViewTransform* transform, vec3_t up, float roll = 0.0f);
// The roll that camera_level would need to give this orientation, in radians in [-pi, pi]
float camera_roll(const ViewTransform& transform, vec3_t up);

// A value that is keyed over time, with the easing of the key that ends each segment (eases[i] shapes the
// stretch from i - 1 to i; eases[0] is not used). Times must be strictly increasing. Holds the first and
// last value outside of them.
double keyed_curve_evaluate(const double* times, const double* values, const KeyEase* eases, size_t n, double time);

// The trajectory frame at 'time', given by the keyframes that have one (use_frame). Returns false if
// none do. The keys need not be sorted. The frame can go backward from one key to the next. It passes through
// each key, holds outside the first and last, and the speed is continuous through the ones in between (a
// monotone cubic, so it never overshoots a key, and it stops where it turns around).
bool camera_keyframes_evaluate_frame(double* out_frame, const CameraKeyframe* keys, size_t count, double time);

// Trajectory frame at 'time' through two anchors, (t_begin, f_begin) and (t_end, f_end), and the keys that
// have a frame. The anchors always take part, so a key in between bends the speed of the trajectory; a key
// at or beyond an anchor takes its place. Before the first point and after the last the frame is held.
double camera_keyframes_frame_with_anchors(const CameraKeyframe* keys, size_t count, double t_begin, double f_begin, double t_end, double f_end, double time);

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

// A good default view of a set of atoms: the whole system, or a structure within it. x, y, z hold all
// num_atoms atoms, as drawn (no periodic treatment); the view frames those listed in indices, or all of them
// when indices is null. cell_A is the unit cell basis (columns a, b, c), or null for none - it is only used
// to judge whether the world axes mean anything. See the implementation for the choices this makes.
ViewTransform camera_compute_default_view(const vec3_t* xyz, size_t num_atoms, const int32_t* indices, size_t count, const mat3_t* cell_A, float fov_y);

// The camera distance at which all atoms (those in indices, or all count of them) fit in a view from
// orientation, looking at look_at.
float camera_fit_distance(const vec3_t* xyz, const int32_t* indices, size_t count, vec3_t look_at, quat_t orientation, float fov_y);

// Lazy stupid procedure on top of interpolate_look_at
void camera_animate(ViewTransform* current, const ViewTransform& target, double dt, double target_factor = 0.12f);