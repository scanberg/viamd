#include <plane.h>

#include <math.h>

static inline vec3_t transform_point(const mat4_t& m, vec3_t p) {
    const vec4_t r = mat4_mul_vec4(m, vec4_set(p.x, p.y, p.z, 1.0f));
    return vec3_set(r.x, r.y, r.z);
}

void plane_frame(vec3_t* u, vec3_t* v, const Plane& p, vec3_t right) {
    const vec3_t n = p.normal;
    // 'right' without its part along n; the next of x, y when it is (nearly) along n: of x and y one is always at least
    // 45 degrees off n (n.x^2 + n.y^2 <= 1), so the last is never too short
    const vec3_t candidates[3] = {right, vec3_set(1, 0, 0), vec3_set(0, 1, 0)};
    vec3_t r = vec3_set(1, 0, 0);
    for (int i = 0; i < 3; ++i) {
        const float len = vec3_length(candidates[i]);
        if (!(len > 0.0f)) continue;
        const vec3_t c = vec3_mul1(candidates[i], 1.0f / len);
        r = vec3_sub(c, vec3_mul1(n, vec3_dot(c, n)));
        if (vec3_length(r) > 0.1f) break;
    }
    *u = vec3_normalize(r);
    *v = vec3_cross(n, *u);
}

int plane_box_polygon(vec3_t out[6], const Plane& plane, const mat4_t& box_to_world) {
    vec3_t c[8];
    float dist[8];
    for (int i = 0; i < 8; ++i) {
        c[i] = transform_point(box_to_world, vec3_set((float)(i & 1), (float)((i >> 1) & 1), (float)((i >> 2) & 1)));
        dist[i] = plane_distance(plane, c[i]);
    }
    // the 12 edges as pairs of corners that differ in one bit
    static const int edge[12][2] = {{0,1},{2,3},{4,5},{6,7},{0,2},{1,3},{4,6},{5,7},{0,4},{1,5},{2,6},{3,7}};
    vec3_t pts[12];
    int n = 0;
    for (int e = 0; e < 12; ++e) {
        const float da = dist[edge[e][0]], db = dist[edge[e][1]];
        if ((da > 0.0f) == (db > 0.0f)) continue;
        const float t = da / (da - db);
        pts[n++] = vec3_lerp(c[edge[e][0]], c[edge[e][1]], t);
    }
    if (n < 3) return 0;
    // in order around their centre, by angle in the plane
    vec3_t centre = vec3_zero();
    for (int i = 0; i < n; ++i) centre = vec3_add(centre, pts[i]);
    centre = vec3_mul1(centre, 1.0f / (float)n);
    vec3_t u, v;
    plane_frame(&u, &v, plane, vec3_set(1, 0, 0));
    float ang[12];
    for (int i = 0; i < n; ++i) {
        const vec3_t d = vec3_sub(pts[i], centre);
        ang[i] = atan2f(vec3_dot(d, v), vec3_dot(d, u));
    }
    for (int i = 1; i < n; ++i) {
        for (int j = i; j > 0 && ang[j] < ang[j - 1]; --j) {
            const float ta = ang[j]; ang[j] = ang[j - 1]; ang[j - 1] = ta;
            const vec3_t tp = pts[j]; pts[j] = pts[j - 1]; pts[j - 1] = tp;
        }
    }
    // an edge through a corner gives the corner twice
    int m = 0;
    for (int i = 0; i < n && m < 6; ++i) {
        if (m > 0 && vec3_distance(pts[i], out[m - 1]) < 1.0e-5f) continue;
        out[m++] = pts[i];
    }
    if (m > 1 && vec3_distance(out[0], out[m - 1]) < 1.0e-5f) --m;
    return m >= 3 ? m : 0;
}

// The eigenvectors of a symmetric 3x3 matrix by Jacobi rotations, in double: the columns of V, with the eigenvalues
// on the diagonal of A after. A covariance of atoms is as well conditioned as their geometry, and three atoms make it
// singular, exactly where the normal comes from: in double that is the exact plane to well below what is drawn.
static void symmetric_eigen3(double A[3][3], double V[3][3]) {
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) V[i][j] = i == j ? 1.0 : 0.0;
    for (int sweep = 0; sweep < 32; ++sweep) {
        const double off = A[0][1] * A[0][1] + A[0][2] * A[0][2] + A[1][2] * A[1][2];
        const double diag = A[0][0] * A[0][0] + A[1][1] * A[1][1] + A[2][2] * A[2][2];
        if (off <= 1.0e-30 * diag || off == 0.0) break;
        for (int p = 0; p < 2; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                if (A[p][q] == 0.0) continue;
                const double theta = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1.0));
                const double c = 1.0 / sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < 3; ++k) {   // A <- A J
                    const double akp = A[k][p], akq = A[k][q];
                    A[k][p] = c * akp - s * akq;
                    A[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k) {   // A <- J^T A
                    const double apk = A[p][k], aqk = A[q][k];
                    A[p][k] = c * apk - s * aqk;
                    A[q][k] = s * apk + c * aqk;
                }
                for (int k = 0; k < 3; ++k) {   // V <- V J
                    const double vkp = V[k][p], vkq = V[k][q];
                    V[k][p] = c * vkp - s * vkq;
                    V[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
}

// The unit vector in the span of a and b (orthonormal) nearest to v, or 0 when v is across both
static vec3_t nearest_in_span(vec3_t v, vec3_t a, vec3_t b) {
    const vec3_t p = vec3_add(vec3_mul1(a, vec3_dot(v, a)), vec3_mul1(b, vec3_dot(v, b)));
    const float len = vec3_length(p);
    return len > 1.0e-4f ? vec3_mul1(p, 1.0f / len) : vec3_zero();
}

Plane plane_fit(const Plane& prev, const vec3_t* points, size_t count, PlaneFitInfo* info) {
    PlaneFitInfo fit = {};
    if (!points || count == 0) {
        if (info) *info = fit;
        return prev;
    }
    double c[3] = {0, 0, 0};
    for (size_t i = 0; i < count; ++i) {
        c[0] += points[i].x; c[1] += points[i].y; c[2] += points[i].z;
    }
    for (int k = 0; k < 3; ++k) c[k] /= (double)count;
    double A[3][3] = {};
    for (size_t i = 0; i < count; ++i) {
        const double d[3] = {points[i].x - c[0], points[i].y - c[1], points[i].z - c[2]};
        for (int r = 0; r < 3; ++r) for (int s = 0; s < 3; ++s) A[r][s] += d[r] * d[s];
    }
    const vec3_t centroid = vec3_set((float)c[0], (float)c[1], (float)c[2]);
    const double spread = A[0][0] + A[1][1] + A[2][2];   // sum of squared distances from the centroid
    vec3_t n = prev.normal;
    if (spread <= 1.0e-8 * (double)count) {
        // all at one place (to 1e-4 A rms)
        fit.kind = PlaneFitKind::Point;
    } else {
        double V[3][3];
        symmetric_eigen3(A, V);
        int order[3] = {0, 1, 2};   // by eigenvalue, largest first
        for (int i = 0; i < 3; ++i) for (int j = i + 1; j < 3; ++j) if (A[order[j]][order[j]] > A[order[i]][order[i]]) { const int t = order[i]; order[i] = order[j]; order[j] = t; }
        const double l0 = A[order[0]][order[0]], l1 = A[order[1]][order[1]], l2 = fmax(A[order[2]][order[2]], 0.0);
        const auto column = [&V](int k) { return vec3_normalize(vec3_set((float)V[0][k], (float)V[1][k], (float)V[2][k])); };
        const vec3_t e1 = column(order[1]), e2 = column(order[2]);
        if (l1 <= 1.0e-6 * l0) {
            // along a line (off it by under a thousandth of its spread along it): the plane through it nearest to prev
            fit.kind = PlaneFitKind::Line;
            n = nearest_in_span(prev.normal, e1, e2);
            // prev across the line: any plane through it is as near, one is taken
            if (vec3_length(n) == 0.0f) n = e2;
        } else if (l2 > PLANE_FIT_LOOSE_RATIO * l0) {
            fit.kind = PlaneFitKind::Kept;
        } else if (l2 > PLANE_FIT_LOOSE_RATIO * l1) {
            fit.kind = PlaneFitKind::Loose;
            n = nearest_in_span(prev.normal, e1, e2);
            // prev along the most spread direction: no nearest, the least-squares one
            if (vec3_length(n) == 0.0f) n = e2;
        } else {
            fit.kind = PlaneFitKind::Fit;
            n = e2;
        }
    }
    if (vec3_dot(n, prev.normal) < 0.0f) n = vec3_mul1(n, -1.0f);
    const Plane plane = plane_from_point_normal(centroid, n);
    // the rms distance to the plane given
    double sum = 0.0;
    for (size_t i = 0; i < count; ++i) {
        const double dist = (double)plane_distance(plane, points[i]);
        sum += dist * dist;
    }
    fit.rms = (float)sqrt(sum / (double)count);
    if (info) *info = fit;
    return plane;
}
