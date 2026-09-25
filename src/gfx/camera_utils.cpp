#include "camera_utils.h"
#include <core/md_common.h>
#include <core/md_vec_math.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>

static vec4_t projection_extents(float fov_y, int width, int height, float texel_offset_x, float texel_offset_y) {
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

mat4_t camera_view_to_world_matrix(const ViewTransform& transform) {
    mat4_t M = mat4_from_quat(transform.orientation);
    M.col[3] = vec4_from_vec3(transform.position, 1);
    return M;
}

mat4_t camera_world_to_view_matrix(const ViewTransform& transform) {
    const mat4_t R = mat4_from_quat(quat_conj(transform.orientation));
    const mat4_t T = mat4_translate(-transform.position.x, -transform.position.y, -transform.position.z);
    return R * T;
}

mat4_t camera_view_to_clip_matrix_persp(const Camera& camera, float aspect_ratio) {
    return mat4_persp(camera.fov_y, aspect_ratio, camera.near_plane, camera.far_plane);
}

mat4_t camera_clip_to_view_matrix_persp(const Camera& camera, float aspect_ratio) {
    return mat4_persp_inv(camera.fov_y, aspect_ratio, camera.near_plane, camera.far_plane);
}

mat4_t camera_view_to_clip_matrix_persp(const Camera& camera, int width, int height, float texel_offset_x, float texel_offset_y) {
    const vec4_t ext = projection_extents(camera.fov_y, width, height, texel_offset_x, texel_offset_y);

    const float cn = camera.near_plane;
    const float cf = camera.far_plane;
    const float xm = ext.z - ext.x;
    const float xp = ext.z + ext.x;
    const float ym = ext.w - ext.y;
    const float yp = ext.w + ext.y;

    return mat4_frustum(xm * cn, xp * cn, ym * cn, yp * cn, cn, cf);
}

mat4_t camera_clip_to_view_matrix_persp(const Camera& camera, int width, int height, float texel_offset_x, float texel_offset_y) {
    const vec4_t ext = projection_extents(camera.fov_y, width, height, texel_offset_x, texel_offset_y);

    const float cn = camera.near_plane;
    const float cf = camera.far_plane;
    const float xm = ext.z - ext.x;
    const float xp = ext.z + ext.x;
    const float ym = ext.w - ext.y;
    const float yp = ext.w + ext.y;

    return mat4_frustum_inv(xm * cn, xp * cn, ym * cn, yp * cn, cn, cf);
}

mat4_t camera_view_to_clip_matrix_ortho(float l, float r, float b, float t) {
    return mat4_ortho_2d(l, r, b, t);
}

mat4_t camera_clip_to_view_matrix_ortho(float l, float r, float b, float t) {
    return mat4_ortho_2d_inv(l, r, b, t);
}

mat4_t camera_view_to_clip_matrix_ortho(float l, float r, float b, float t, float n, float f) {
    return mat4_ortho(l, r, b, t, n, f);
}

mat4_t camera_clip_to_view_matrix_ortho(float l, float r, float b, float t, float n, float f) {
    return mat4_ortho_inv(l, r, b, t, n, f);
}

mat4_t look_at(vec3_t look_from, vec3_t look_at, vec3_t look_up) {
    const vec3_t f = vec3_normalize(look_at - look_from);
    const vec3_t s = vec3_normalize(vec3_cross(f, look_up));
    const vec3_t u = vec3_cross(s, f);
	const mat4_t M = {
            s.x, u.x, -f.x, 0.0f,
            s.y, u.y, -f.y, 0.0f,
            s.z, u.z, -f.z, 0.0f,
            -vec3_dot(s, look_from), -vec3_dot(u, look_from), vec3_dot(f, look_from), 1.0f,
    };
    return M;
}

static inline float project_to_sphere(float r, vec2_t v) {
    const float d = vec2_length(v);
    if (d < r * 0.70710678118654752440f) {
        // On sphere
        return sqrtf(r * r - d * d);
    } else {
        // On hyperbola
        float t = r / 1.41421356237309504880f;
        return t * t / d;
    }
}

static inline quat_t trackball(vec2_t prev_ndc, vec2_t curr_ndc) {
    static const float TRACKBALLSIZE = 0.8f;
    if (vec2_dot(curr_ndc - prev_ndc, curr_ndc - prev_ndc) == 0.f) return quat_t{1, 0, 0, 0};

    const vec3_t p1 = vec3_from_vec2(prev_ndc, project_to_sphere(TRACKBALLSIZE, prev_ndc));
    const vec3_t p2 = vec3_from_vec2(curr_ndc, project_to_sphere(TRACKBALLSIZE, curr_ndc));

    const vec3_t axis = vec3_normalize(vec3_cross(p2, p1));
    float t = CLAMP(vec3_length(p1 - p2) / (2.0f * TRACKBALLSIZE), -1.f, 1.f);
    float angle = 2.f * asinf(t);

    return quat_axis_angle(axis, angle);
}

void camera_trackball(Camera* camera, vec2_t prev_ndc, vec2_t curr_ndc) {
    ASSERT(camera);
    const quat_t q = trackball(prev_ndc, curr_ndc);
    camera->orientation = camera->orientation * q;
}

void camera_move(Camera* camera, vec3_t t) {
    ASSERT(camera);
    camera->position = camera->position + camera->orientation * t;
}

vec3_t camera_get_look_at(const ViewTransform& transform) {
    return transform.position - transform.orientation * vec3_t{0, 0, transform.distance};
}

vec3_t camera_position_from_look_at(const vec3_t& look_at, const quat_t& orientation, float distance) {
    return look_at + orientation * vec3_t{0, 0, distance};
}

// high precision version
static inline vec3_t highp_quat_vec3_mul(const quat_t& q, const vec3_t& v) {
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
    return vec3_set((float)res[0], (float)res[1], (float)res[2]);
}

// We want to interpolate along an arc which is formed by maintaining a distance to the look_at position and smoothly interpolating the orientation,
// We linearly interpolate a look_at position which is implicitly defined by position, orientation and distance
// There is some precision errors creeping into the posision because we transform back and forth to look at using the orientation
void camera_interpolate_look_at(vec3_t* out_pos, quat_t* out_ori, float* out_dist, vec3_t in_pos[2], quat_t in_ori[2], float in_dist[2], double t) {

    // This is to combat floating point inaccuracies when interpolating in this arc like fashion
    // If we are close enough, we just set it to that
    const double quat_epsilon = 1.0e-7;
    const double pos_epsilon = 1.0e-7;
    const double dist_epsilon = 1.0e-7;

    quat_t ori;
    vec3_t pos;
    float dist;

    const float qd = quat_dot(in_ori[0], in_ori[1]);
    if (fabsf(qd) > 1.0 - quat_epsilon) {
        ori = in_ori[1];
    } else {
        // Due to the inherent rotational duality of quaternions we want to make sure we rotate along the shortest 'path'
        quat_t qa = qd < 0.0f ? quat_conj(in_ori[0]) : in_ori[0];
        quat_t qb = in_ori[1];
        ori = quat_normalize(quat_slerp(qa, qb, (float)t));
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
        vec3_t l0 = in_pos[0] - highp_quat_vec3_mul(in_ori[0], vec3_set(0, 0, in_dist[0]));
        vec3_t l1 = in_pos[1] - highp_quat_vec3_mul(in_ori[1], vec3_set(0, 0, in_dist[1]));
        vec3_t look_at = {
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
    ASSERT(transform);

    const vec2_t half_res = input.screen_size * 0.5f;
    const vec2_t ndc_prev = (vec2_t{input.mouse_coord_prev.x, input.screen_size.y - input.mouse_coord_prev.y} - half_res) / half_res;
    const vec2_t ndc_curr = (vec2_t{input.mouse_coord_curr.x, input.screen_size.y - input.mouse_coord_curr.y} - half_res) / half_res;
    const vec2_t mouse_coord_delta = input.mouse_coord_curr - input.mouse_coord_prev;
    const bool mouse_move = mouse_coord_delta != vec2_t{0, 0};

    if ((flags & TrackballFlags_RotateEnabled) && input.rotate_button && mouse_move) {
        const quat_t q = trackball(ndc_prev, ndc_curr);
        const vec3_t look_at = camera_get_look_at(*transform);
        transform->orientation = quat_normalize(transform->orientation * q);
        transform->position = camera_position_from_look_at(look_at, transform->orientation, transform->distance);
        if (flags & TrackballFlags_RotateReturnsTrue) return true;
    } else if ((flags & TrackballFlags_PanEnabled) && input.pan_button && mouse_move) {
        const float aspect_ratio = input.screen_size.x / input.screen_size.y;
        const float scl = tanf(input.fov_y * 0.5f);
        const vec2_t delta = (ndc_curr - ndc_prev) * vec2_t{1, -1} * vec2_t{aspect_ratio * scl, scl};
        const vec3_t move = transform->orientation * vec3_t{-delta.x, delta.y, 0} * powf(transform->distance * param.pan_scale, param.pan_exponent);
        transform->position = transform->position + move;
        if (flags & TrackballFlags_PanReturnsTrue) return true;
    } else if ((flags & TrackballFlags_DollyEnabled) && ((input.dolly_button && mouse_move) || input.dolly_delta != 0.f)) {
        float delta = -(input.mouse_coord_curr.y - input.mouse_coord_prev.y) * powf(transform->distance * param.dolly_drag_scale, param.dolly_drag_exponent);
        delta -= input.dolly_delta * powf(transform->distance * param.dolly_delta_scale, param.dolly_delta_exponent);
        const vec3_t look_at = camera_get_look_at(*transform);
        transform->distance = CLAMP(transform->distance + delta, param.min_distance, param.max_distance);
        transform->position = camera_position_from_look_at(look_at, transform->orientation, transform->distance);
        if (flags & TrackballFlags_DollyReturnsTrue) return true;
    }
    return false;
}

ViewTransform compute_optimal_view(const vec3_t& center, const vec3_t& half_ext, const mat3_t& basis, float distance_scale) {
    const float len = MAX(vec3_length(half_ext * 0.5f), 5.0f);
    const float max_ext = MAX(MAX(half_ext.x, half_ext.y), half_ext.z);
    const float min_ext = MIN(MIN(half_ext.x, half_ext.y), half_ext.z);

	vec3_t up  = vec3_set(0, 0, 1);
    vec3_t dir = vec3_normalize(vec3_set(0.6f, 0.5f, 1.0f));
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
        vec3_t right = basis[l[0]];
        vec3_t out   = basis[l[2]];
        dir = vec3_normalize(right * dir.x + up * dir.y + out * dir.z);
    }

    const vec3_t pos = center + dir * len * distance_scale;

    ViewTransform result = {
        .orientation = quat_from_mat4(mat4_look_at(pos, center, up)),
        .position = pos,
        .distance = vec3_length(pos - center),
	};
    
    return result;
}

void camera_animate(ViewTransform* current, const ViewTransform& target, double dt, double target_factor) {
    ASSERT(current);

    dt = CLAMP(dt, 1.0 / 1000.0, 1.0 / 20.0);

    // We use an exponential interpolation of the deltas with a common factor
    const double INV_TARGET_DT = 100.0;
    double interpolation_factor = target_factor * dt * INV_TARGET_DT;

    vec3_t pos[2] = {current->position, target.position};
    quat_t ori[2] = {current->orientation, target.orientation};
    float dist[2] = {current->distance, target.distance};
    camera_interpolate_look_at(&current->position, &current->orientation, &current->distance, pos, ori, dist, interpolation_factor);
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
// backed off until every atom fits.

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

static const vec3_t DV_X = {1, 0, 0};
static const vec3_t DV_Y = {0, 1, 0};
static const vec3_t DV_Z = {0, 0, 1};

struct DvFrame {
    vec3_t face;  // Towards the viewer, before the swing
    vec3_t right;
    vec3_t up;
    bool   swing;
};

// The part of v orthogonal to the unit vector n, normalized. False when v is (nearly) parallel to n.
static bool dv_reject(vec3_t* out, vec3_t v, vec3_t n) {
    const vec3_t r = vec3_sub(v, vec3_mul1(n, vec3_dot(v, n)));
    const float  l = vec3_length(r);
    if (l < 1.0e-3f) return false;
    *out = vec3_div1(r, l);
    return true;
}

static vec3_t dv_orient(vec3_t v, vec3_t pref) {
    return vec3_dot(v, pref) < 0.0f ? vec3_mul1(v, -1.0f) : v;
}

static vec3_t dv_swing(const DvFrame& f) {
    if (!f.swing) return f.face;
    const float ca = cosf(DV_AZIMUTH),   sa = sinf(DV_AZIMUTH);
    const float ce = cosf(DV_ELEVATION), se = sinf(DV_ELEVATION);
    return vec3_normalize(vec3_add(vec3_add(vec3_mul1(f.face, ce * ca), vec3_mul1(f.right, ce * sa)), vec3_mul1(f.up, se)));
}

// The world view's camera frame (b towards the camera, s screen right, u screen up). It has a component along
// every world axis, which is what makes it a tie breaker for vectors lying exactly along one.
static void dv_world_camera(vec3_t* b, vec3_t* s, vec3_t* u) {
    const DvFrame world = {DV_X, DV_Y, DV_Z, true};
    *b = dv_swing(world);
    *s = vec3_normalize(vec3_cross(DV_Z, *b));
    *u = vec3_cross(*b, *s);
}

// A unit vector orthogonal to the unit vector n, as close to the world X (or Y) as possible
static vec3_t dv_perp(vec3_t n) {
    vec3_t v;
    if (!dv_reject(&v, DV_X, n)) dv_reject(&v, DV_Y, n);
    return v;
}

static vec3_t dv_load(const vec3_t* xyz, const int32_t* indices, size_t i) {
    const size_t idx = indices ? (size_t)indices[i] : i;
    return xyz[idx];
}

// Does the system fill its cell, i.e. span it along at least two cell axes?
static bool dv_fills_cell(const vec3_t* xyz, size_t num_atoms, const mat3_t& A) {
    if (!num_atoms || fabsf(mat3_determinant(A)) < 1.0e-6f) return false;
    const mat3_t I = mat3_inverse(A);
    vec3_t fmin = vec3_set1( FLT_MAX);
    vec3_t fmax = vec3_set1(-FLT_MAX);
    for (size_t i = 0; i < num_atoms; ++i) {
        const vec3_t f = mat3_mul_vec3(I, xyz[i]);
        fmin = vec3_min(fmin, f);
        fmax = vec3_max(fmax, f);
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
static vec3_t dv_face_about_up(vec3_t up, const mat3_t* C) {
    vec3_t pref_b, pref_s, pref_u;
    dv_world_camera(&pref_b, &pref_s, &pref_u);

    const vec3_t e1 = dv_perp(up);
    const vec3_t e2 = vec3_cross(up, e1);
    vec3_t face = e1;
    if (C) {
        const float a = vec3_dot(e1, mat3_mul_vec3(*C, e1));
        const float b = vec3_dot(e1, mat3_mul_vec3(*C, e2));
        const float c = vec3_dot(e2, mat3_mul_vec3(*C, e2));
        const float h = 0.5f * (a + c);
        const float r = sqrtf(MAX(0.0f, 0.25f * (a - c) * (a - c) + b * b));
        const float l_major = h + r;
        const float l_minor = h - r;
        if (l_major > l_minor * DV_ANISOTROPY * DV_ANISOTROPY) {
            // Eigenvector of the minor eigenvalue; of the two equivalent forms take the better conditioned
            const vec2_t v0 = {b, l_minor - a};
            const vec2_t v1 = {l_minor - c, b};
            const vec2_t v  = (v0.x * v0.x + v0.y * v0.y > v1.x * v1.x + v1.y * v1.y) ? v0 : v1;
            face = vec3_normalize(vec3_add(vec3_mul1(e1, v.x), vec3_mul1(e2, v.y)));
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

static float dv_visibility(DvRaster& r, const vec3_t* p, size_t n, vec3_t d) {
    const vec3_t s = dv_perp(d);
    const vec3_t u = vec3_cross(d, s);
    const int    dim = r.dim;
    const float  rad = DV_ATOM_RADIUS;

    for (int k = 0; k < dim * dim; ++k) { r.zbuf[k] = -FLT_MAX; r.ibuf[k] = -1; }
    for (size_t i = 0; i < n; ++i) { r.disc[i] = 0; r.vis[i] = 0; }

    for (size_t i = 0; i < n; ++i) {
        const float sx = vec3_dot(p[i], s);
        const float sy = vec3_dot(p[i], u);
        const float sz = vec3_dot(p[i], d);
        const int x0 = MAX(0,       (int)floorf((sx - rad + r.R) / r.px));
        const int x1 = MIN(dim - 1, (int)floorf((sx + rad + r.R) / r.px));
        const int y0 = MAX(0,       (int)floorf((sy - rad + r.R) / r.px));
        const int y1 = MIN(dim - 1, (int)floorf((sy + rad + r.R) / r.px));
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
static vec3_t dv_best_direction(const vec3_t* p, size_t n, const vec3_t axis[3]) {
    vec3_t pref_b, pref_s, pref_u;
    dv_world_camera(&pref_b, &pref_s, &pref_u);

    float R = 0.0f;
    for (size_t i = 0; i < n; ++i) R = MAX(R, vec3_length(p[i]));
    R += DV_ATOM_RADIUS;

    DvRaster r = {};
    r.dim  = CLAMP((int)ceilf(2.0f * R / (DV_ATOM_RADIUS / 3.0f)), 32, 256);
    r.px   = 2.0f * R / (float)r.dim;
    r.R    = R;
    r.zbuf = (float*)  malloc(sizeof(float)   * r.dim * r.dim);
    r.ibuf = (int32_t*)malloc(sizeof(int32_t) * r.dim * r.dim);
    r.disc = (int32_t*)malloc(sizeof(int32_t) * n);
    r.vis  = (int32_t*)malloc(sizeof(int32_t) * n);

    // The principal axes first, with a small bonus: an exact face-on or edge-on view beats a sampled
    // direction a few degrees off it that scores the same
    vec3_t best = pref_b;
    float  best_score = -1.0f;
    for (int k = 0; k < 6; ++k) {
        const vec3_t d = vec3_mul1(axis[k / 2], (k & 1) ? -1.0f : 1.0f);
        const float  score = dv_visibility(r, p, n, d) * 1.01f;
        if (score > best_score) { best_score = score; best = d; }
    }
    // Then the sphere, evenly (Fibonacci)
    const float golden_angle = 2.39996323f;
    for (int k = 0; k < DV_SEARCH_DIRS; ++k) {
        const float cz = 1.0f - (2.0f * k + 1.0f) / (float)DV_SEARCH_DIRS;
        const float sr = sqrtf(MAX(0.0f, 1.0f - cz * cz));
        const vec3_t d = vec3_set(sr * cosf(golden_angle * k), sr * sinf(golden_angle * k), cz);
        const float  score = dv_visibility(r, p, n, d);
        if (score > best_score) { best_score = score; best = d; }
    }
    // Symmetric molecules score the same from both sides: then take the side facing the world view
    if (vec3_dot(best, pref_b) < 0.0f) {
        const vec3_t flip = vec3_mul1(best, -1.0f);
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
static DvFrame dv_frame_about_face(vec3_t b, const vec3_t* p, size_t n) {
    vec3_t pref_b, pref_s, pref_u;
    dv_world_camera(&pref_b, &pref_s, &pref_u);

    const vec3_t e1 = dv_perp(b);
    const vec3_t e2 = vec3_cross(b, e1);
    float a = 0, c = 0, bb = 0;
    for (size_t i = 0; i < n; ++i) {
        const float x = vec3_dot(p[i], e1);
        const float y = vec3_dot(p[i], e2);
        a += x * x; bb += x * y; c += y * y;
    }
    const float h = 0.5f * (a + c);
    const float r = sqrtf(MAX(0.0f, 0.25f * (a - c) * (a - c) + bb * bb));
    const float l_major = h + r;
    const float l_minor = h - r;

    DvFrame f = {b, e1, e2, false};
    if (l_major > l_minor * DV_ANISOTROPY * DV_ANISOTROPY) {
        const vec2_t v0 = {bb, l_major - a};
        const vec2_t v1 = {l_major - c, bb};
        const vec2_t v  = (v0.x * v0.x + v0.y * v0.y > v1.x * v1.x + v1.y * v1.y) ? v0 : v1;
        f.right = dv_orient(vec3_normalize(vec3_add(vec3_mul1(e1, v.x), vec3_mul1(e2, v.y))), pref_s);
        f.up    = vec3_cross(b, f.right);
    } else if (dv_reject(&f.up, pref_u, b)) {
        f.right = vec3_cross(f.up, b);
    } else {
        dv_reject(&f.right, pref_s, b);
        f.up = vec3_cross(b, f.right);
    }
    if (vec3_dot(f.up, pref_u) < 0.0f) {
        // Upright takes priority over the sign of right: half a turn about the view direction
        f.up    = vec3_mul1(f.up,    -1.0f);
        f.right = vec3_mul1(f.right, -1.0f);
    }
    return f;
}

float camera_fit_distance(const vec3_t* xyz, const int32_t* indices, size_t count, vec3_t look_at, quat_t orientation, float fov_y) {
    const vec3_t s = quat_mul_vec3(orientation, DV_X);
    const vec3_t u = quat_mul_vec3(orientation, DV_Y);
    const vec3_t b = quat_mul_vec3(orientation, DV_Z);
    const float tan_half_fov = tanf(fov_y * 0.5f) * DV_FILL;
    float dist = DV_MIN_DISTANCE;
    for (size_t i = 0; i < count; ++i) {
        const vec3_t p = vec3_sub(dv_load(xyz, indices, i), look_at);
        const float lateral = MAX(fabsf(vec3_dot(p, s)), fabsf(vec3_dot(p, u))) + DV_ATOM_RADIUS;
        dist = MAX(dist, vec3_dot(p, b) + lateral / tan_half_fov);
    }
    return dist;
}

ViewTransform camera_compute_default_view(const vec3_t* xyz, size_t num_atoms, const int32_t* indices, size_t count, const mat3_t* cell_A, float fov_y) {
    ViewTransform result = {};
    if (!indices) count = num_atoms;
    if (!count) return result;

    // Shape: covariance and principal axes about the mean, extents along them
    vec3_t mean = {0, 0, 0};
    for (size_t i = 0; i < count; ++i) mean = vec3_add(mean, dv_load(xyz, indices, i));
    mean = vec3_div1(mean, (float)count);
    const mat3_t C   = mat3_covariance_matrix(xyz, nullptr, indices, count, mean);
    const mat3_t PCA = mat3_orthonormalize(mat3_extract_rotation(mat3_eigen(C).vectors));
    const mat3_t basis = mat3_transpose(PCA); // Axis i is column i of basis (row i of PCA)

    vec3_t pmin = vec3_set1( FLT_MAX);
    vec3_t pmax = vec3_set1(-FLT_MAX);
    for (size_t i = 0; i < count; ++i) {
        const vec3_t q = mat3_mul_vec3(PCA, dv_load(xyz, indices, i));
        pmin = vec3_min(pmin, q);
        pmax = vec3_max(pmax, q);
    }
    // Sorted longest first; the atom radius keeps a flat or linear shape from reading as infinitely anisotropic
    int l[3] = {0, 1, 2};
    float e[3];
    for (int k = 0; k < 3; ++k) e[k] = (pmax.elem[k] - pmin.elem[k]) + 2.0f * DV_ATOM_RADIUS;
    for (int pass = 0; pass < 2; ++pass)
        for (int k = 0; k < 2; ++k)
            if (e[l[k]] < e[l[k + 1]]) { int t = l[k]; l[k] = l[k + 1]; l[k + 1] = t; }
    const vec3_t axis[3] = {basis.col[l[0]], basis.col[l[1]], basis.col[l[2]]};
    const float  ext[3]  = {e[l[0]], e[l[1]], e[l[2]]};

    vec3_t pref_b, pref_s, pref_u;
    dv_world_camera(&pref_b, &pref_s, &pref_u);

    DvFrame f = {DV_X, DV_Y, DV_Z, true};
    const bool whole = (indices == nullptr) || count == num_atoms;
    const bool slab  = ext[2] <= ext[1] * DV_SLAB_FLATNESS && ext[1] >= DV_SLAB_MIN_SIZE;

    if (cell_A && dv_fills_cell(xyz, num_atoms, *cell_A)) {
        // The world frame means something: Z up
        if (!whole) f.face = dv_face_about_up(DV_Z, &C);
    } else if (slab) {
        // Membrane or surface: its normal is up, seen from the side
        vec3_t up = axis[2];
        const vec3_t world[3] = {DV_X, DV_Y, DV_Z};
        for (int k = 0; k < 3; ++k) {
            if (fabsf(vec3_dot(up, world[k])) >= DV_SNAP_COS) up = world[k];
        }
        f.up   = dv_orient(up, DV_Z);
        f.face = dv_face_about_up(f.up, &C);
    } else if (count <= DV_SEARCH_MAX_ATOMS) {
        // Small molecule: the view that shows the most of it
        vec3_t* p = (vec3_t*)malloc(sizeof(vec3_t) * count);
        for (size_t i = 0; i < count; ++i) p[i] = vec3_sub(dv_load(xyz, indices, i), mean);
        const vec3_t b = dv_best_direction(p, count, axis);
        f = dv_frame_about_face(b, p, count);
        free(p);
    } else if (ext[0] > ext[2] * DV_ANISOTROPY) {
        // Large molecule: by its shape
        f.face  = dv_orient(axis[2], pref_b);
        f.right = dv_orient(axis[0], pref_s);
        f.up    = vec3_cross(f.face, f.right);
        if (vec3_dot(f.up, pref_u) < 0.0f) {
            f.up    = vec3_mul1(f.up,    -1.0f);
            f.right = vec3_mul1(f.right, -1.0f);
        }
    }
    // else: a large, round molecule shows the same from any side - the world view

    if (f.swing) {
        // Complete the frame from face and up (right handed: face = right x up)
        f.right = vec3_normalize(vec3_cross(f.up, f.face));
        f.up    = vec3_cross(f.face, f.right);
    }

    // Camera frame: b from the look-at towards the camera, s screen right, u screen up
    const vec3_t b = dv_swing(f);
    const vec3_t s = vec3_normalize(vec3_cross(f.up, b));
    const vec3_t u = vec3_cross(b, s);

    // Center on the middle of the extent as seen in that frame
    vec3_t vmin = vec3_set1( FLT_MAX);
    vec3_t vmax = vec3_set1(-FLT_MAX);
    for (size_t i = 0; i < count; ++i) {
        const vec3_t q = vec3_sub(dv_load(xyz, indices, i), mean);
        const vec3_t v = vec3_set(vec3_dot(q, s), vec3_dot(q, u), vec3_dot(q, b));
        vmin = vec3_min(vmin, v);
        vmax = vec3_max(vmax, v);
    }
    const vec3_t mid    = vec3_mul1(vec3_add(vmin, vmax), 0.5f);
    const vec3_t center = vec3_add(mean, vec3_add(vec3_add(vec3_mul1(s, mid.x), vec3_mul1(u, mid.y)), vec3_mul1(b, mid.z)));

    const quat_t orientation = quat_from_mat4(mat4_look_at(vec3_add(center, b), center, u));
    const float  dist = camera_fit_distance(xyz, indices, count, center, orientation, fov_y);

    result.orientation = orientation;
    result.position    = vec3_add(center, vec3_mul1(b, dist));
    result.distance    = dist;
    return result;
}
