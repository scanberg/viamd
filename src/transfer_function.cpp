#include <transfer_function.h>

#include <core/md_hash.h>

#include <math.h>
#include <string.h>

static inline float clamp01(float x) { return x < 0.0f ? 0.0f : (x > 1.0f ? 1.0f : x); }

float transfer_function_param(const TransferFunction& tf, float value) {
    const float ext = tf.range_end - tf.range_beg;
    if (!(fabsf(ext) > 0.0f)) return 0.5f;
    return clamp01((value - tf.range_beg) / ext);
}

float transfer_function_value(const TransferFunction& tf, float t) {
    return tf.range_beg + t * (tf.range_end - tf.range_beg);
}

float transfer_function_curve(const TransferFunction& tf, float t) {
    const int n = tf.num_points;
    if (n <= 0) return 0.0f;
    const TransferFunctionPoint* p = tf.points;
    if (t <= p[0].t) return p[0].alpha;
    if (t >= p[n - 1].t) return p[n - 1].alpha;
    for (int i = 1; i < n; ++i) {
        if (t <= p[i].t) {
            const float dt = p[i].t - p[i - 1].t;
            const float s = dt > 0.0f ? (t - p[i - 1].t) / dt : 1.0f;
            return p[i - 1].alpha + s * (p[i].alpha - p[i - 1].alpha);
        }
    }
    return p[n - 1].alpha;
}

void transfer_function_bake_alpha(const TransferFunction& tf, float* out, int n, int stride) {
    if (!out || n <= 0) return;
    const float opacity = clamp01(tf.opacity);
    for (int i = 0; i < n; ++i) {
        const float t = n > 1 ? (float)i / (float)(n - 1) : 0.0f;
        out[(size_t)i * (size_t)stride] = opacity * clamp01(transfer_function_curve(tf, t));
    }
}

int transfer_function_insert(TransferFunction* tf, float t, float alpha) {
    if (!tf || tf->num_points >= TRANSFER_FUNCTION_MAX_POINTS) return -1;
    t = clamp01(t);
    alpha = clamp01(alpha);
    int i = 0;
    while (i < tf->num_points && tf->points[i].t <= t) ++i;
    memmove(tf->points + i + 1, tf->points + i, sizeof(TransferFunctionPoint) * (size_t)(tf->num_points - i));
    tf->points[i] = {t, alpha};
    tf->num_points += 1;
    return i;
}

bool transfer_function_remove(TransferFunction* tf, int i) {
    if (!tf || tf->num_points <= 1 || i < 0 || i >= tf->num_points) return false;
    memmove(tf->points + i, tf->points + i + 1, sizeof(TransferFunctionPoint) * (size_t)(tf->num_points - i - 1));
    tf->num_points -= 1;
    return true;
}

void transfer_function_move(TransferFunction* tf, int i, float t, float alpha) {
    if (!tf || i < 0 || i >= tf->num_points) return;
    const float lo = i > 0 ? tf->points[i - 1].t : 0.0f;
    const float hi = i + 1 < tf->num_points ? tf->points[i + 1].t : 1.0f;
    t = clamp01(t);
    tf->points[i].t = t < lo ? lo : (t > hi ? hi : t);
    tf->points[i].alpha = clamp01(alpha);
}

void transfer_function_set_range(TransferFunction* tf, float beg, float end) {
    if (!tf) return;
    if (tf->symmetric) {
        const float r = fmaxf(fabsf(beg), fabsf(end));
        tf->range_beg = -r;
        tf->range_end = r;
    } else {
        tf->range_beg = beg;
        tf->range_end = end;
    }
}

uint64_t transfer_function_hash(const TransferFunction& tf) {
    uint64_t h = md_hash64(&tf.colormap, sizeof(tf.colormap), 0);
    const float f[3] = {tf.range_beg, tf.range_end, tf.opacity};
    h = md_hash64(f, sizeof(f), h);
    h = md_hash64(&tf.num_points, sizeof(tf.num_points), h);
    const int n = tf.num_points < 0 ? 0 : (tf.num_points > TRANSFER_FUNCTION_MAX_POINTS ? TRANSFER_FUNCTION_MAX_POINTS : tf.num_points);
    return md_hash64(tf.points, sizeof(TransferFunctionPoint) * (size_t)n, h);
}

void transfer_function_ramp(TransferFunction* tf, float from) {
    if (!tf) return;
    from = clamp01(from);
    tf->num_points = 0;
    if (from > 0.0f) tf->points[tf->num_points++] = {0.0f, 0.0f};
    tf->points[tf->num_points++] = {from, 0.0f};
    tf->points[tf->num_points++] = {1.0f, 1.0f};
}

void transfer_function_ramp_magnitude(TransferFunction* tf, float from) {
    if (!tf) return;
    const float h = 0.5f * clamp01(from);
    tf->num_points = 0;
    tf->points[tf->num_points++] = {0.0f, 1.0f};
    tf->points[tf->num_points++] = {0.5f - h, 0.0f};
    if (h > 0.0f) tf->points[tf->num_points++] = {0.5f + h, 0.0f};
    tf->points[tf->num_points++] = {1.0f, 1.0f};
}
