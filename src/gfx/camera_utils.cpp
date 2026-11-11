#include "camera_utils.h"
#include <core/md_common.h>
#include <core/md_vec_math.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>

static md_f32x4_t projection_extents(float fov_y, int width, int height, float texel_offset_x, float texel_offset_y) {
    const float aspect_ratio = (float)width / (float)height;
    const float half_h = tanf(fov_y * 0.5f);
    const float half_w = aspect_ratio * half_h;
    const float texel_size_x = half_w / (float)(0.5f * width);
    const float texel_size_y = half_h / (float)(0.5f * height);
    const float jitter_x = texel_size_x * texel_offset_x;
    const float jitter_y = texel_size_y * texel_offset_y;

    // xy = frustum extents at distance 1, zw = jitter at distance 1
    return {half_w, half_h, jitter_x, jitter_y};
}

md_f32x4x4_t camera_view_to_world_matrix(const ViewTransform& transform) {
    md_f32x4x4_t M = md_f32x4x4_from_quat(transform.orientation);
    M.col[3] = md_f32x4_from_f32x3(transform.position, 1);
    return M;
}

md_f32x4x4_t camera_world_to_view_matrix(const ViewTransform& transform) {
    const md_f32x4x4_t R = md_f32x4x4_from_quat(md_quat_conj(transform.orientation));
    const md_f32x4x4_t T = md_f32x4x4_translate(-transform.position.x, -transform.position.y, -transform.position.z);
    return R * T;
}

md_f32x4x4_t camera_view_to_clip_matrix_persp(const Camera& camera, float aspect_ratio) {
    return md_f32x4x4_persp(camera.fov_y, aspect_ratio, camera.near_plane, camera.far_plane);
}

md_f32x4x4_t camera_clip_to_view_matrix_persp(const Camera& camera, float aspect_ratio) {
    return md_f32x4x4_persp_inv(camera.fov_y, aspect_ratio, camera.near_plane, camera.far_plane);
}

md_f32x4x4_t camera_view_to_clip_matrix_persp(const Camera& camera, int width, int height, float texel_offset_x, float texel_offset_y) {
    const md_f32x4_t ext = projection_extents(camera.fov_y, width, height, texel_offset_x, texel_offset_y);

    const float cn = camera.near_plane;
    const float cf = camera.far_plane;
    const float xm = ext.z - ext.x;
    const float xp = ext.z + ext.x;
    const float ym = ext.w - ext.y;
    const float yp = ext.w + ext.y;

    return md_f32x4x4_frustum(xm * cn, xp * cn, ym * cn, yp * cn, cn, cf);
}

md_f32x4x4_t camera_clip_to_view_matrix_persp(const Camera& camera, int width, int height, float texel_offset_x, float texel_offset_y) {
    const md_f32x4_t ext = projection_extents(camera.fov_y, width, height, texel_offset_x, texel_offset_y);

    const float cn = camera.near_plane;
    const float cf = camera.far_plane;
    const float xm = ext.z - ext.x;
    const float xp = ext.z + ext.x;
    const float ym = ext.w - ext.y;
    const float yp = ext.w + ext.y;

    return md_f32x4x4_frustum_inv(xm * cn, xp * cn, ym * cn, yp * cn, cn, cf);
}

md_f32x4x4_t camera_view_to_clip_matrix_ortho(float l, float r, float b, float t) {
    return md_f32x4x4_ortho_2d(l, r, b, t);
}

md_f32x4x4_t camera_clip_to_view_matrix_ortho(float l, float r, float b, float t) {
    return md_f32x4x4_ortho_2d_inv(l, r, b, t);
}

md_f32x4x4_t camera_view_to_clip_matrix_ortho(float l, float r, float b, float t, float n, float f) {
    return md_f32x4x4_ortho(l, r, b, t, n, f);
}

md_f32x4x4_t camera_clip_to_view_matrix_ortho(float l, float r, float b, float t, float n, float f) {
    return md_f32x4x4_ortho_inv(l, r, b, t, n, f);
}

md_f32x4x4_t look_at(md_f32x3_t look_from, md_f32x3_t look_at, md_f32x3_t look_up) {
    const md_f32x3_t f = md_f32x3_normalize(look_at - look_from);
    const md_f32x3_t s = md_f32x3_normalize(md_f32x3_cross(f, look_up));
    const md_f32x3_t u = md_f32x3_cross(s, f);
	const md_f32x4x4_t M = {
            s.x, u.x, -f.x, 0.0f,
            s.y, u.y, -f.y, 0.0f,
            s.z, u.z, -f.z, 0.0f,
            -md_f32x3_dot(s, look_from), -md_f32x3_dot(u, look_from), md_f32x3_dot(f, look_from), 1.0f,
    };
    return M;
}

static inline float project_to_sphere(float r, md_f32x2_t v) {
    const float d = md_f32x2_length(v);
    if (d < r * 0.70710678118654752440f) {
        // On sphere
        return sqrtf(r * r - d * d);
    } else {
        // On hyperbola
        float t = r / 1.41421356237309504880f;
        return t * t / d;
    }
}

static inline md_quat_t trackball(md_f32x2_t prev_ndc, md_f32x2_t curr_ndc) {
    static const float TRACKBALLSIZE = 0.8f;
    if (md_f32x2_dot(curr_ndc - prev_ndc, curr_ndc - prev_ndc) == 0.f) return md_quat_t{1, 0, 0, 0};

    const md_f32x3_t p1 = md_f32x3_from_f32x2(prev_ndc, project_to_sphere(TRACKBALLSIZE, prev_ndc));
    const md_f32x3_t p2 = md_f32x3_from_f32x2(curr_ndc, project_to_sphere(TRACKBALLSIZE, curr_ndc));

    const md_f32x3_t axis = md_f32x3_normalize(md_f32x3_cross(p2, p1));
    float t = MD_CLAMP(md_f32x3_length(p1 - p2) / (2.0f * TRACKBALLSIZE), -1.f, 1.f);
    float angle = 2.f * asinf(t);

    return md_quat_axis_angle(axis, angle);
}

void camera_trackball(Camera* camera, md_f32x2_t prev_ndc, md_f32x2_t curr_ndc) {
    MD_ASSERT(camera);
    const md_quat_t q = trackball(prev_ndc, curr_ndc);
    camera->orientation = camera->orientation * q;
}

void camera_move(Camera* camera, md_f32x3_t t) {
    MD_ASSERT(camera);
    camera->position = camera->position + camera->orientation * t;
}

md_f32x3_t camera_get_look_at(const ViewTransform& transform) {
    return transform.position - transform.orientation * md_f32x3_t{0, 0, transform.distance};
}

md_f32x3_t camera_position_from_look_at(const md_f32x3_t& look_at, const md_quat_t& orientation, float distance) {
    return look_at + orientation * md_f32x3_t{0, 0, distance};
}

// high precision version
static inline md_f32x3_t highp_quat_vec3_mul(const md_quat_t& q, const md_f32x3_t& v) {
    double u[4] = {q.x, q.y, q.z, q.w};
    double t[3] = {
        2.0 * (u[1] * (double)v.z - (double)v.y * u[2]),
        2.0 * (u[2] * (double)v.x - (double)v.z * u[0]),
        2.0 * (u[0] * (double)v.y - (double)v.x * u[1]),
    };
    double res[3] = {
        v.x + t[0] * u[3] + (u[1] * t[2] - t[1] * u[2]),
        v.y + t[1] * u[3] + (u[2] * t[0] - t[2] * u[0]),
        v.z + t[2] * u[3] + (u[0] * t[1] - t[0] * u[1]),
    };
    return md_f32x3_set((float)res[0], (float)res[1], (float)res[2]);
}

// We want to interpolate along an arc which is formed by maintaining a distance to the look_at position and smoothly interpolating the orientation,
// We linearly interpolate a look_at position which is implicitly defined by position, orientation and distance
// There is some precision errors creeping into the posision because we transform back and forth to look at using the orientation
void camera_interpolate_look_at(md_f32x3_t* out_pos, md_quat_t* out_ori, float* out_dist, md_f32x3_t in_pos[2], md_quat_t in_ori[2], float in_dist[2], double t) {

    // This is to combat floating point inaccuracies when interpolating in this arc like fashion
    // If we are close enough, we just set it to that
    const double quat_epsilon = 1.0e-7;
    const double pos_epsilon = 1.0e-7;
    const double dist_epsilon = 1.0e-7;

    md_quat_t ori;
    md_f32x3_t pos;
    float dist;

    const float qd = md_quat_dot(in_ori[0], in_ori[1]);
    if (fabsf(qd) > 1.0 - quat_epsilon) {
        ori = in_ori[1];
    } else {
        // Due to the inherent rotational duality of quaternions we want to make sure we rotate along the shortest 'path'
        md_quat_t qa = qd < 0.0f ? md_quat_conj(in_ori[0]) : in_ori[0];
        md_quat_t qb = in_ori[1];
        ori = md_quat_normalize(md_quat_slerp(qa, qb, (float)t));
    }

    if (fabs((double)in_dist[0] - (double)in_dist[1]) < dist_epsilon) {
        dist = in_dist[1];
    } else {
        dist = (float)lerp(in_dist[0], in_dist[1], t);
    }

    double dx = in_pos[0].x - in_pos[1].x;
    double dy = in_pos[0].y - in_pos[1].y;
    double dz = in_pos[0].z - in_pos[1].z;
    double d2 = sqrt(dx*dx + dy*dy + dz*dz);
    if (d2 < pos_epsilon) {
        pos = in_pos[1];
    } else {
        md_f32x3_t l0 = in_pos[0] - highp_quat_vec3_mul(in_ori[0], md_f32x3_set(0, 0, in_dist[0]));
        md_f32x3_t l1 = in_pos[1] - highp_quat_vec3_mul(in_ori[1], md_f32x3_set(0, 0, in_dist[1]));
        md_f32x3_t look_at = {
            lerp(l0.x, l1.x, t),
            lerp(l0.y, l1.y, t),
            lerp(l0.z, l1.z, t),
        };
        pos = camera_position_from_look_at(look_at, ori, dist);
    }

    *out_pos = pos;
    *out_ori = ori;
    *out_dist = dist;
}

bool camera_controller_trackball(ViewTransform* transform, const TrackballControllerInput& input, const TrackballControllerParam& param, TrackballFlags flags) {
    MD_ASSERT(transform);

    const md_f32x2_t half_res = input.screen_size * 0.5f;
    const md_f32x2_t ndc_prev = (md_f32x2_t{input.mouse_coord_prev.x, input.screen_size.y - input.mouse_coord_prev.y} - half_res) / half_res;
    const md_f32x2_t ndc_curr = (md_f32x2_t{input.mouse_coord_curr.x, input.screen_size.y - input.mouse_coord_curr.y} - half_res) / half_res;
    const md_f32x2_t mouse_coord_delta = input.mouse_coord_curr - input.mouse_coord_prev;
    const bool mouse_move = mouse_coord_delta != md_f32x2_t{0, 0};

    if ((flags & TrackballFlags_RotateEnabled) && input.rotate_button && mouse_move) {
        const md_quat_t q = trackball(ndc_prev - input.ndc_center, ndc_curr - input.ndc_center);
        const md_f32x3_t look_at = camera_get_look_at(*transform);
        transform->orientation = md_quat_normalize(transform->orientation * q);
        transform->position = camera_position_from_look_at(look_at, transform->orientation, transform->distance);
        if (flags & TrackballFlags_RotateReturnsTrue) return true;
    } else if ((flags & TrackballFlags_PanEnabled) && input.pan_button && mouse_move) {
        const float aspect_ratio = input.screen_size.x / input.screen_size.y;
        const float scl = tanf(input.fov_y * 0.5f);
        const md_f32x2_t delta = (ndc_curr - ndc_prev) * md_f32x2_t{1, -1} * md_f32x2_t{aspect_ratio * scl, scl};
        const md_f32x3_t move = transform->orientation * md_f32x3_t{-delta.x, delta.y, 0} * powf(transform->distance * param.pan_scale, param.pan_exponent);
        transform->position = transform->position + move;
        if (flags & TrackballFlags_PanReturnsTrue) return true;
    } else if ((flags & TrackballFlags_DollyEnabled) && ((input.dolly_button && mouse_move) || input.dolly_delta != 0.f)) {
        float delta = -(input.mouse_coord_curr.y - input.mouse_coord_prev.y) * powf(transform->distance * param.dolly_drag_scale, param.dolly_drag_exponent);
        delta -= input.dolly_delta * powf(transform->distance * param.dolly_delta_scale, param.dolly_delta_exponent);
        const md_f32x3_t look_at = camera_get_look_at(*transform);
        transform->distance = MD_CLAMP(transform->distance + delta, param.min_distance, param.max_distance);
        transform->position = camera_position_from_look_at(look_at, transform->orientation, transform->distance);
        if (flags & TrackballFlags_DollyReturnsTrue) return true;
    }
    return false;
}

ViewTransform compute_optimal_view(const md_f32x3_t& center, const md_f32x3_t& half_ext, const md_f32x3x3_t& basis, float distance_scale) {
    const float len = MD_MAX(md_f32x3_length(half_ext * 0.5f), 5.0f);
    const float max_ext = MD_MAX(MD_MAX(half_ext.x, half_ext.y), half_ext.z);
    const float min_ext = MD_MIN(MD_MIN(half_ext.x, half_ext.y), half_ext.z);

	md_f32x3_t up  = md_f32x3_set(0, 0, 1);
    md_f32x3_t dir = md_f32x3_normalize(md_f32x3_set(0.6f, 0.5f, 1.0f));
    if (max_ext > 0.0f) {
        const float aniso_ext = max_ext / min_ext;

        // We want to align the view such that we the longest axis of the aabb align with the X-axis, the mid axis with the Y-axis
        int l[3] = { 0, 1, 2 };
        auto swap = [](int& a, int& b) { int t = a; a = b; b = t; };
        if (aniso_ext > 1.1f) {
            // The aabb is not uniform, so we sort the axes by length
            if (half_ext[l[0]] < half_ext[l[1]]) swap(l[0], l[1]);
            if (half_ext[l[1]] < half_ext[l[2]]) swap(l[1], l[2]);
            if (half_ext[l[0]] < half_ext[l[1]]) swap(l[0], l[1]);
            // Now the axes are sorted with respect to the length l[0] > l[1] > l[2]
        }

               up    = basis[l[1]];
        md_f32x3_t right = basis[l[0]];
        md_f32x3_t out   = basis[l[2]];
        dir = md_f32x3_normalize(right * dir.x + up * dir.y + out * dir.z);
    }

    const md_f32x3_t pos = center + dir * len * distance_scale;

    ViewTransform result = {
        .orientation = md_quat_from_f32x4x4(md_f32x4x4_look_at(pos, center, up)),
        .position = pos,
        .distance = md_f32x3_length(pos - center),
	};
    
    return result;
}

// We use an exponential interpolation of the deltas with a common factor
static double animate_interpolation_factor(double dt, double target_factor) {
    dt = MD_CLAMP(dt, 1.0 / 1000.0, 1.0 / 20.0);
    const double INV_TARGET_DT = 100.0;
    return target_factor * dt * INV_TARGET_DT;
}

// Within these of the target, an animation is at the target: the exponential approach would go on in the last bits
// for ever, and a view that never stops moving is never still enough to keep (lazy rendering). Far below a pixel.
static const double SNAP_SHIFT    = 1.0e-5;     // of the shift, in clip space (a pixel is about 1e-3)
static const double SNAP_RELATIVE = 1.0e-5;     // of the position and the distance, relative to the distance
static const double SNAP_ANGLE    = 1.0e-5;     // of the orientation, in radians (a pixel is about 1e-4 at most)

void camera_animate_shift(md_f32x2_t* current, md_f32x2_t target, double dt, double target_factor) {
    MD_ASSERT(current);
    const float t = (float)MD_CLAMP(animate_interpolation_factor(dt, target_factor), 0.0, 1.0);
    *current = *current + (target - *current) * t;
    if (fabs((double)current->x - (double)target.x) < SNAP_SHIFT && fabs((double)current->y - (double)target.y) < SNAP_SHIFT) {
        *current = target;
    }
}

void camera_animate(ViewTransform* current, const ViewTransform& target, double dt, double target_factor) {
    MD_ASSERT(current);

    double interpolation_factor = animate_interpolation_factor(dt, target_factor);

    md_f32x3_t pos[2] = {current->position, target.position};
    md_quat_t ori[2] = {current->orientation, target.orientation};
    float dist[2] = {current->distance, target.distance};
    camera_interpolate_look_at(&current->position, &current->orientation, &current->distance, pos, ori, dist, interpolation_factor);

    // At the target when within a fraction of a pixel of it. The angle between the orientations by the length of their
    // difference (about half the angle), which floats resolve where the dot product does not.
    const double scale = MD_MAX((double)target.distance, 1.0e-3);
    double d_pos = 0.0, dot = 0.0, d_ori_minus = 0.0, d_ori_plus = 0.0;
    for (int i = 0; i < 3; ++i) {
        const double d = (double)current->position.elem[i] - (double)target.position.elem[i];
        d_pos += d * d;
    }
    for (int i = 0; i < 4; ++i) {
        const double a = current->orientation.elem[i], b = target.orientation.elem[i];
        dot += a * b;
        d_ori_minus += (a - b) * (a - b);
        d_ori_plus  += (a + b) * (a + b);
    }
    const double angle = 2.0 * sqrt(dot < 0.0 ? d_ori_plus : d_ori_minus);
    if (sqrt(d_pos) < SNAP_RELATIVE * scale && fabs((double)current->distance - (double)target.distance) < SNAP_RELATIVE * scale &&
        angle < SNAP_ANGLE) {
        *current = target;
    }
}

// -----------------------------------------------------------------------------------------------------------
// Default view
//
// What makes a good default view depends on whether the world axes mean anything.
//
// * A system that fills its periodic cell (spans it along at least two of the cell axes) was built in a frame
//   that means something: membranes and slabs are set up with their normal along Z and the box edges lie
//   along the axes. Z is up, and the whole system is shown in the world view: X towards the viewer and
//   slightly to the left, seen a little from above. A structure inside such a system keeps Z up and only
//   turns about it to show its broad side. A cell the atoms do not fill - a crystallographic cell around a
//   protein, a placeholder CRYST1 - says nothing about the frame and is ignored.
// * A large, flat system without a cell saying so - a membrane patch or a surface - is a slab: its normal is
//   up (snapped to a world axis when close to one) and it is seen from the side like the world view.
// * Anything else is a molecule whose coordinates carry no preferred direction. A small one gets the view in
//   which the most of its atoms can be seen: sampled directions are scored by the visible share of each
//   atom, in the spirit of viewpoint entropy. That finds face-on for planar molecules and avoids looking
//   down bonds or along rows where atoms eclipse each other. A large one is too dense for visibility to
//   discriminate - only its surface shows from any side - so it is shown by its shape: broad side towards
//   the viewer, long axis across the screen.
//
// Signs and in-plane rotations that the shape leaves open are settled by the world view, so the same input
// always gives the same view. In every case the view is centered on the projected extent and the camera is
// backed off until every atom fits, in the region of the view it is made for.

static const float  DV_AZIMUTH          = 25.0f * (3.14159265f / 180.0f); // Swing from the face towards screen right
static const float  DV_ELEVATION        = 20.0f * (3.14159265f / 180.0f); // Swing up from the horizontal
static const float  DV_ATOM_RADIUS      = 1.0f;
static const float  DV_FILL             = 0.9f;   // Fraction of the view the fitted atoms may fill
static const float  DV_MIN_DISTANCE     = 15.0f;
static const float  DV_CELL_SPAN_MIN    = 0.75f;  // Span of the atoms along a cell axis, in cell lengths, that counts as filling it
static const float  DV_CELL_SPAN_MAX    = 1.5f;   // Beyond this the cell is too small to be the box of the system (placeholder)
static const float  DV_SLAB_FLATNESS    = 0.6f;   // A slab's thickness is at most this fraction of its middle extent...
static const float  DV_SLAB_MIN_SIZE    = 50.0f;  // ...and its middle extent is at least this long (Ångström)
static const float  DV_SNAP_COS         = 0.966f; // cos(15 deg): a slab normal this close to a world axis is snapped to it
static const float  DV_ANISOTROPY       = 1.2f;   // Extent ratio above which an axis of a shape is trusted
static const size_t DV_SEARCH_MAX_ATOMS = 1500;   // Above this, a molecule is shown by its shape rather than by visibility
static const int    DV_SEARCH_DIRS      = 192;

static const md_f32x3_t DV_X = {1, 0, 0};
static const md_f32x3_t DV_Y = {0, 1, 0};
static const md_f32x3_t DV_Z = {0, 0, 1};

struct DvFrame {
    md_f32x3_t face;  // Towards the viewer, before the swing
    md_f32x3_t right;
    md_f32x3_t up;
    bool   swing;
};

// The part of v orthogonal to the unit vector n, normalized. False when v is (nearly) parallel to n.
static bool dv_reject(md_f32x3_t* out, md_f32x3_t v, md_f32x3_t n) {
    const md_f32x3_t r = md_f32x3_sub(v, md_f32x3_mul1(n, md_f32x3_dot(v, n)));
    const float  l = md_f32x3_length(r);
    if (l < 1.0e-3f) return false;
    *out = md_f32x3_div1(r, l);
    return true;
}

static md_f32x3_t dv_orient(md_f32x3_t v, md_f32x3_t pref) {
    return md_f32x3_dot(v, pref) < 0.0f ? md_f32x3_mul1(v, -1.0f) : v;
}

static md_f32x3_t dv_swing(const DvFrame& f) {
    if (!f.swing) return f.face;
    const float ca = cosf(DV_AZIMUTH),   sa = sinf(DV_AZIMUTH);
    const float ce = cosf(DV_ELEVATION), se = sinf(DV_ELEVATION);
    return md_f32x3_normalize(md_f32x3_add(md_f32x3_add(md_f32x3_mul1(f.face, ce * ca), md_f32x3_mul1(f.right, ce * sa)), md_f32x3_mul1(f.up, se)));
}

// The world view's camera frame (b towards the camera, s screen right, u screen up). It has a component along
// every world axis, which is what makes it a tie breaker for vectors lying exactly along one.
static void dv_world_camera(md_f32x3_t* b, md_f32x3_t* s, md_f32x3_t* u) {
    const DvFrame world = {DV_X, DV_Y, DV_Z, true};
    *b = dv_swing(world);
    *s = md_f32x3_normalize(md_f32x3_cross(DV_Z, *b));
    *u = md_f32x3_cross(*b, *s);
}

// A unit vector orthogonal to the unit vector n, as close to the world X (or Y) as possible
static md_f32x3_t dv_perp(md_f32x3_t n) {
    md_f32x3_t v;
    if (!dv_reject(&v, DV_X, n)) dv_reject(&v, DV_Y, n);
    return v;
}

static md_f32x3_t dv_load(const md_f32x3_t* xyz, const int32_t* indices, size_t i) {
    const size_t idx = indices ? (size_t)indices[i] : i;
    return xyz[idx];
}

// Does the system fill its cell, i.e. span it along at least two cell axes?
static bool dv_fills_cell(const md_f32x3_t* xyz, size_t num_atoms, const md_f32x3x3_t& A) {
    if (!num_atoms || fabsf(md_f32x3x3_determinant(A)) < 1.0e-6f) return false;
    const md_f32x3x3_t I = md_f32x3x3_inverse(A);
    md_f32x3_t fmin = md_f32x3_set1( FLT_MAX);
    md_f32x3_t fmax = md_f32x3_set1(-FLT_MAX);
    for (size_t i = 0; i < num_atoms; ++i) {
        const md_f32x3_t f = md_f32x3x3_mul_f32x3(I, xyz[i]);
        fmin = md_f32x3_min(fmin, f);
        fmax = md_f32x3_max(fmax, f);
    }
    int filled = 0;
    for (int k = 0; k < 3; ++k) {
        const float span = fmax.elem[k] - fmin.elem[k];
        filled += (DV_CELL_SPAN_MIN <= span && span <= DV_CELL_SPAN_MAX) ? 1 : 0;
    }
    return filled >= 2;
}

// With up given: the face that shows the broad side of a shape with covariance C, when it has one about up.
// Otherwise the world X, as seen from the world view.
static md_f32x3_t dv_face_about_up(md_f32x3_t up, const md_f32x3x3_t* C) {
    md_f32x3_t pref_b, pref_s, pref_u;
    dv_world_camera(&pref_b, &pref_s, &pref_u);

    const md_f32x3_t e1 = dv_perp(up);
    const md_f32x3_t e2 = md_f32x3_cross(up, e1);
    md_f32x3_t face = e1;
    if (C) {
        const float a = md_f32x3_dot(e1, md_f32x3x3_mul_f32x3(*C, e1));
        const float b = md_f32x3_dot(e1, md_f32x3x3_mul_f32x3(*C, e2));
        const float c = md_f32x3_dot(e2, md_f32x3x3_mul_f32x3(*C, e2));
        const float h = 0.5f * (a + c);
        const float r = sqrtf(MD_MAX(0.0f, 0.25f * (a - c) * (a - c) + b * b));
        const float l_major = h + r;
        const float l_minor = h - r;
        if (l_major > l_minor * DV_ANISOTROPY * DV_ANISOTROPY) {
            // Eigenvector of the minor eigenvalue; of the two equivalent forms take the better conditioned
            const md_f32x2_t v0 = {b, l_minor - a};
            const md_f32x2_t v1 = {l_minor - c, b};
            const md_f32x2_t v  = (v0.x * v0.x + v0.y * v0.y > v1.x * v1.x + v1.y * v1.y) ? v0 : v1;
            face = md_f32x3_normalize(md_f32x3_add(md_f32x3_mul1(e1, v.x), md_f32x3_mul1(e2, v.y)));
        }
    }
    return dv_orient(face, pref_b);
}

// Score of a view along d (towards the camera): the sum over atoms of the square root of their visible
// share of the projected disc. Concave, so that seeing every atom partly beats seeing some fully and others
// not at all, and never saturating, so that a view with less overlap always scores higher. Orthographic,
// rasterized with sphere depth.
struct DvRaster {
    int      dim;
    float    px;
    float    R;
    float*   zbuf;
    int32_t* ibuf;
    int32_t* disc;
    int32_t* vis;
};

static float dv_visibility(DvRaster& r, const md_f32x3_t* p, size_t n, md_f32x3_t d) {
    const md_f32x3_t s = dv_perp(d);
    const md_f32x3_t u = md_f32x3_cross(d, s);
    const int    dim = r.dim;
    const float  rad = DV_ATOM_RADIUS;

    for (int k = 0; k < dim * dim; ++k) { r.zbuf[k] = -FLT_MAX; r.ibuf[k] = -1; }
    for (size_t i = 0; i < n; ++i) { r.disc[i] = 0; r.vis[i] = 0; }

    for (size_t i = 0; i < n; ++i) {
        const float sx = md_f32x3_dot(p[i], s);
        const float sy = md_f32x3_dot(p[i], u);
        const float sz = md_f32x3_dot(p[i], d);
        const int x0 = MD_MAX(0,       (int)floorf((sx - rad + r.R) / r.px));
        const int x1 = MD_MIN(dim - 1, (int)floorf((sx + rad + r.R) / r.px));
        const int y0 = MD_MAX(0,       (int)floorf((sy - rad + r.R) / r.px));
        const int y1 = MD_MIN(dim - 1, (int)floorf((sy + rad + r.R) / r.px));
        for (int py = y0; py <= y1; ++py) {
            const float dy = (py + 0.5f) * r.px - r.R - sy;
            for (int px = x0; px <= x1; ++px) {
                const float dx = (px + 0.5f) * r.px - r.R - sx;
                const float rho2 = dx * dx + dy * dy;
                if (rho2 > rad * rad) continue;
                r.disc[i] += 1;
                const float depth = sz + sqrtf(rad * rad - rho2);
                const int   k = py * dim + px;
                if (depth > r.zbuf[k]) {
                    r.zbuf[k] = depth;
                    r.ibuf[k] = (int32_t)i;
                }
            }
        }
    }
    for (int k = 0; k < dim * dim; ++k) {
        if (r.ibuf[k] >= 0) r.vis[r.ibuf[k]] += 1;
    }
    float score = 0.0f;
    for (size_t i = 0; i < n; ++i) {
        if (r.disc[i] > 0) score += sqrtf((float)r.vis[i] / (float)r.disc[i]);
    }
    return score;
}

// The direction (towards the camera) from which the most of the atoms p (relative to their center) are seen.
static md_f32x3_t dv_best_direction(const md_f32x3_t* p, size_t n, const md_f32x3_t axis[3]) {
    md_f32x3_t pref_b, pref_s, pref_u;
    dv_world_camera(&pref_b, &pref_s, &pref_u);

    float R = 0.0f;
    for (size_t i = 0; i < n; ++i) R = MD_MAX(R, md_f32x3_length(p[i]));
    R += DV_ATOM_RADIUS;

    DvRaster r = {};
    r.dim  = MD_CLAMP((int)ceilf(2.0f * R / (DV_ATOM_RADIUS / 3.0f)), 32, 256);
    r.px   = 2.0f * R / (float)r.dim;
    r.R    = R;
    r.zbuf = (float*)  malloc(sizeof(float)   * r.dim * r.dim);
    r.ibuf = (int32_t*)malloc(sizeof(int32_t) * r.dim * r.dim);
    r.disc = (int32_t*)malloc(sizeof(int32_t) * n);
    r.vis  = (int32_t*)malloc(sizeof(int32_t) * n);

    // The principal axes first, with a small bonus: an exact face-on or edge-on view beats a sampled
    // direction a few degrees off it that scores the same
    md_f32x3_t best = pref_b;
    float  best_score = -1.0f;
    for (int k = 0; k < 6; ++k) {
        const md_f32x3_t d = md_f32x3_mul1(axis[k / 2], (k & 1) ? -1.0f : 1.0f);
        const float  score = dv_visibility(r, p, n, d) * 1.01f;
        if (score > best_score) { best_score = score; best = d; }
    }
    // Then the sphere, evenly (Fibonacci)
    const float golden_angle = 2.39996323f;
    for (int k = 0; k < DV_SEARCH_DIRS; ++k) {
        const float cz = 1.0f - (2.0f * k + 1.0f) / (float)DV_SEARCH_DIRS;
        const float sr = sqrtf(MD_MAX(0.0f, 1.0f - cz * cz));
        const md_f32x3_t d = md_f32x3_set(sr * cosf(golden_angle * k), sr * sinf(golden_angle * k), cz);
        const float  score = dv_visibility(r, p, n, d);
        if (score > best_score) { best_score = score; best = d; }
    }
    // Symmetric molecules score the same from both sides: then take the side facing the world view
    if (md_f32x3_dot(best, pref_b) < 0.0f) {
        const md_f32x3_t flip = md_f32x3_mul1(best, -1.0f);
        if (dv_visibility(r, p, n, flip) >= best_score * 0.995f) best = flip;
    }

    free(r.zbuf);
    free(r.ibuf);
    free(r.disc);
    free(r.vis);
    return best;
}

// Given the view direction b: right along the long axis of the projection, up following from it, upright
// with respect to the world view where the sign is free.
static DvFrame dv_frame_about_face(md_f32x3_t b, const md_f32x3_t* p, size_t n) {
    md_f32x3_t pref_b, pref_s, pref_u;
    dv_world_camera(&pref_b, &pref_s, &pref_u);

    const md_f32x3_t e1 = dv_perp(b);
    const md_f32x3_t e2 = md_f32x3_cross(b, e1);
    float a = 0, c = 0, bb = 0;
    for (size_t i = 0; i < n; ++i) {
        const float x = md_f32x3_dot(p[i], e1);
        const float y = md_f32x3_dot(p[i], e2);
        a += x * x; bb += x * y; c += y * y;
    }
    const float h = 0.5f * (a + c);
    const float r = sqrtf(MD_MAX(0.0f, 0.25f * (a - c) * (a - c) + bb * bb));
    const float l_major = h + r;
    const float l_minor = h - r;

    DvFrame f = {b, e1, e2, false};
    if (l_major > l_minor * DV_ANISOTROPY * DV_ANISOTROPY) {
        const md_f32x2_t v0 = {bb, l_major - a};
        const md_f32x2_t v1 = {l_major - c, bb};
        const md_f32x2_t v  = (v0.x * v0.x + v0.y * v0.y > v1.x * v1.x + v1.y * v1.y) ? v0 : v1;
        f.right = dv_orient(md_f32x3_normalize(md_f32x3_add(md_f32x3_mul1(e1, v.x), md_f32x3_mul1(e2, v.y))), pref_s);
        f.up    = md_f32x3_cross(b, f.right);
    } else if (dv_reject(&f.up, pref_u, b)) {
        f.right = md_f32x3_cross(f.up, b);
    } else {
        dv_reject(&f.right, pref_s, b);
        f.up = md_f32x3_cross(b, f.right);
    }
    if (md_f32x3_dot(f.up, pref_u) < 0.0f) {
        // Upright takes priority over the sign of right: half a turn about the view direction
        f.up    = md_f32x3_mul1(f.up,    -1.0f);
        f.right = md_f32x3_mul1(f.right, -1.0f);
    }
    return f;
}

float camera_fit_distance(const md_f32x3_t* xyz, const int32_t* indices, size_t count, md_f32x3_t look_at, md_quat_t orientation, float fov_y, const ViewRegion& region) {
    const md_f32x3_t s = md_quat_mul_f32x3(orientation, DV_X);
    const md_f32x3_t u = md_quat_mul_f32x3(orientation, DV_Y);
    const md_f32x3_t b = md_quat_mul_f32x3(orientation, DV_Z);
    // The shift moves the picture over without changing it, so only the size of the region matters: half the
    // side of the largest square in it, as the tangent of an angle off the view axis
    const md_f32x2_t half_ext = (region.max - region.min) * 0.5f;
    const float tan_half_side = MD_MIN(half_ext.x * region.aspect_ratio, half_ext.y) * tanf(fov_y * 0.5f) * DV_FILL;
    float dist = DV_MIN_DISTANCE;
    for (size_t i = 0; i < count; ++i) {
        const md_f32x3_t p = md_f32x3_sub(dv_load(xyz, indices, i), look_at);
        const float lateral = MD_MAX(fabsf(md_f32x3_dot(p, s)), fabsf(md_f32x3_dot(p, u))) + DV_ATOM_RADIUS;
        dist = MD_MAX(dist, md_f32x3_dot(p, b) + lateral / tan_half_side);
    }
    return dist;
}

ViewTransform camera_compute_default_view(const md_f32x3_t* xyz, size_t num_atoms, const int32_t* indices, size_t count, const md_f32x3x3_t* cell_A, float fov_y, const ViewRegion& region) {
    ViewTransform result = {};
    if (!indices) count = num_atoms;
    if (!count) return result;

    // Shape: covariance and principal axes about the mean, extents along them
    md_f32x3_t mean = {0, 0, 0};
    for (size_t i = 0; i < count; ++i) mean = md_f32x3_add(mean, dv_load(xyz, indices, i));
    mean = md_f32x3_div1(mean, (float)count);
    const md_f32x3x3_t C   = md_f32x3x3_covariance_matrix(xyz, nullptr, indices, count, mean);
    const md_f32x3x3_t PCA = md_f32x3x3_orthonormalize(md_f32x3x3_extract_rotation(md_f32x3x3_eigen(C).vectors));
    const md_f32x3x3_t basis = md_f32x3x3_transpose(PCA); // Axis i is column i of basis (row i of PCA)

    md_f32x3_t pmin = md_f32x3_set1( FLT_MAX);
    md_f32x3_t pmax = md_f32x3_set1(-FLT_MAX);
    for (size_t i = 0; i < count; ++i) {
        const md_f32x3_t q = md_f32x3x3_mul_f32x3(PCA, dv_load(xyz, indices, i));
        pmin = md_f32x3_min(pmin, q);
        pmax = md_f32x3_max(pmax, q);
    }
    // Sorted longest first; the atom radius keeps a flat or linear shape from reading as infinitely anisotropic
    int l[3] = {0, 1, 2};
    float e[3];
    for (int k = 0; k < 3; ++k) e[k] = (pmax.elem[k] - pmin.elem[k]) + 2.0f * DV_ATOM_RADIUS;
    for (int pass = 0; pass < 2; ++pass)
        for (int k = 0; k < 2; ++k)
            if (e[l[k]] < e[l[k + 1]]) { int t = l[k]; l[k] = l[k + 1]; l[k + 1] = t; }
    const md_f32x3_t axis[3] = {basis.col[l[0]], basis.col[l[1]], basis.col[l[2]]};
    const float  ext[3]  = {e[l[0]], e[l[1]], e[l[2]]};

    md_f32x3_t pref_b, pref_s, pref_u;
    dv_world_camera(&pref_b, &pref_s, &pref_u);

    DvFrame f = {DV_X, DV_Y, DV_Z, true};
    const bool whole = (indices == nullptr) || count == num_atoms;
    const bool slab  = ext[2] <= ext[1] * DV_SLAB_FLATNESS && ext[1] >= DV_SLAB_MIN_SIZE;

    if (cell_A && dv_fills_cell(xyz, num_atoms, *cell_A)) {
        // The world frame means something: Z up
        if (!whole) f.face = dv_face_about_up(DV_Z, &C);
    } else if (slab) {
        // Membrane or surface: its normal is up, seen from the side
        md_f32x3_t up = axis[2];
        const md_f32x3_t world[3] = {DV_X, DV_Y, DV_Z};
        for (int k = 0; k < 3; ++k) {
            if (fabsf(md_f32x3_dot(up, world[k])) >= DV_SNAP_COS) up = world[k];
        }
        f.up   = dv_orient(up, DV_Z);
        f.face = dv_face_about_up(f.up, &C);
    } else if (count <= DV_SEARCH_MAX_ATOMS) {
        // Small molecule: the view that shows the most of it
        md_f32x3_t* p = (md_f32x3_t*)malloc(sizeof(md_f32x3_t) * count);
        for (size_t i = 0; i < count; ++i) p[i] = md_f32x3_sub(dv_load(xyz, indices, i), mean);
        const md_f32x3_t b = dv_best_direction(p, count, axis);
        f = dv_frame_about_face(b, p, count);
        free(p);
    } else if (ext[0] > ext[2] * DV_ANISOTROPY) {
        // Large molecule: by its shape
        f.face  = dv_orient(axis[2], pref_b);
        f.right = dv_orient(axis[0], pref_s);
        f.up    = md_f32x3_cross(f.face, f.right);
        if (md_f32x3_dot(f.up, pref_u) < 0.0f) {
            f.up    = md_f32x3_mul1(f.up,    -1.0f);
            f.right = md_f32x3_mul1(f.right, -1.0f);
        }
    }
    // else: a large, round molecule shows the same from any side - the world view

    if (f.swing) {
        // Complete the frame from face and up (right handed: face = right x up)
        f.right = md_f32x3_normalize(md_f32x3_cross(f.up, f.face));
        f.up    = md_f32x3_cross(f.face, f.right);
    }

    // Camera frame: b from the look-at towards the camera, s screen right, u screen up
    const md_f32x3_t b = dv_swing(f);
    const md_f32x3_t s = md_f32x3_normalize(md_f32x3_cross(f.up, b));
    const md_f32x3_t u = md_f32x3_cross(b, s);

    // Center on the middle of the extent as seen in that frame
    md_f32x3_t vmin = md_f32x3_set1( FLT_MAX);
    md_f32x3_t vmax = md_f32x3_set1(-FLT_MAX);
    for (size_t i = 0; i < count; ++i) {
        const md_f32x3_t q = md_f32x3_sub(dv_load(xyz, indices, i), mean);
        const md_f32x3_t v = md_f32x3_set(md_f32x3_dot(q, s), md_f32x3_dot(q, u), md_f32x3_dot(q, b));
        vmin = md_f32x3_min(vmin, v);
        vmax = md_f32x3_max(vmax, v);
    }
    const md_f32x3_t mid    = md_f32x3_mul1(md_f32x3_add(vmin, vmax), 0.5f);
    const md_f32x3_t center = md_f32x3_add(mean, md_f32x3_add(md_f32x3_add(md_f32x3_mul1(s, mid.x), md_f32x3_mul1(u, mid.y)), md_f32x3_mul1(b, mid.z)));

    const md_quat_t orientation = md_quat_from_f32x4x4(md_f32x4x4_look_at(md_f32x3_add(center, b), center, u));
    const float  dist = camera_fit_distance(xyz, indices, count, center, orientation, fov_y, region);

    result.orientation = orientation;
    result.position    = md_f32x3_add(center, md_f32x3_mul1(b, dist));
    result.distance    = dist;
    return result;
}
