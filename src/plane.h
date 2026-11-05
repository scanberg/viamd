#pragma once

// A plane (Plane: a vec4) and what is computed with one, apart from the application: the axes laid out on it, the
// polygon it cuts out of a box, and the plane that fits a set of points. Used for the planes of the electronic
// structure representations (es_planes.h).

#include <core/md_vec_math.h>
#include <stddef.h>

// A plane in space, n . x + d = 0 with n of unit length (world, Angstrom), as a vec4. As a cut it keeps the side n
// points to, n . x + d >= 0; flipping it is negating all four. It has no position or turn within itself: what is laid
// out on it (points, a gizmo) takes plane_frame's axes, and a gizmo sits where the plane cuts what it is drawn over.
struct Plane {
    vec3_t normal = {0, 0, 1};
    float  d = 0.0f;
};
static_assert(sizeof(Plane) == sizeof(vec4_t), "a plane is a vec4");

static inline vec4_t plane_vec4(const Plane& p) { return vec4_set(p.normal.x, p.normal.y, p.normal.z, p.d); }
static inline float  plane_distance(const Plane& p, vec3_t x) { return vec3_dot(p.normal, x) + p.d; }
static inline vec3_t plane_project(const Plane& p, vec3_t x) { return vec3_sub(x, vec3_mul1(p.normal, plane_distance(p, x))); }
static inline void   plane_flip(Plane* p) { p->normal = vec3_mul1(p->normal, -1.0f); p->d = -p->d; }
// The plane through 'point' facing 'normal' (any length but 0)
static inline Plane  plane_from_point_normal(vec3_t point, vec3_t normal) {
    Plane p;
    p.normal = vec3_normalize(normal);
    p.d = -vec3_dot(p.normal, point);
    return p;
}
// Axes in the plane, u and v, with u x v = n: u is 'right' made perpendicular to n (world x, or else y, when 'right'
// is too close to n for that), v = n x u. The turn of the plane about n is not the plane's: it is taken from what the
// axes are for, such as the camera's right for a gizmo, so that they line up with the screen.
void plane_frame(vec3_t* u, vec3_t* v, const Plane& p, vec3_t right);

// The polygon a plane cuts out of a box, the unit cube of 'box_to_world' (affine): up to 6 corners in world space,
// in order around it. Returns how many (0 when it misses the box).
int plane_box_polygon(vec3_t out[6], const Plane& plane, const mat4_t& box_to_world);

// How plane_fit placed the plane: what the points determined of it
enum class PlaneFitKind {
    None,       // no points: the plane as it was
    Point,      // one point, or all at one place: moved through it, facing as it did
    Line,       // along a line: the plane through the line nearest the previous one
    Loose,      // the two least spread directions alike: the normal between them, nearest the previous one
    Kept,       // spread alike every way (a tetrahedron, a blob): moved through the centroid, facing as it did
    Fit,        // a plane: the least-squares one (exact for three points)
};

struct PlaneFitInfo {
    PlaneFitKind kind = PlaneFitKind::None;
    float rms = 0.0f;   // of the points' distances to the plane given
};

// The two smallest spreads of the points (eigenvalues of their scatter) closer than this, as a ratio, leave the normal
// undetermined between their directions; all three closer than this, everywhere. At 0.5 the plane along the second
// direction fits only twice as badly as along the first (as rms distance, 1.4 times).
#define PLANE_FIT_LOOSE_RATIO 0.5

// The plane through a set of points, turned as little from 'prev' as they allow, and on the same side (the normal
// within 90 degrees of prev's, so what it cuts away stays on the side it was). Always through their centroid; the
// normal (PlaneFitKind):
// - points spread in two directions much more than in the third: along the third, the least-squares plane, which
//   minimises the sum of squared perpendicular distances (exact for three points not in a line);
// - points along a line: the plane through the line nearest to prev (its normal: prev's without the part along the
//   line);
// - the two least spread directions alike (PLANE_FIT_LOOSE_RATIO): any normal between them fits about as well, and
//   the one nearest prev's is taken, rather than one that turns with the smallest change in the points;
// - spread alike in all three: prev's.
// No points: prev.
Plane plane_fit(const Plane& prev, const vec3_t* points, size_t count, PlaneFitInfo* info = nullptr);
