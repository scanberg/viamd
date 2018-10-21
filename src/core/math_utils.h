#pragma once

#define NOMINMAX
//#ifndef GLM_FORCE_SSE2
//#define GLM_FORCE_SSE2
//#endif

//#ifndef GLM_FORCE_ALIGNED
//#define GLM_FORCE_ALIGNED
//#endif

#include <core/types.h>
#include <core/vector_types.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/random.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/spline.hpp>
#include <stdlib.h>

namespace math {
constexpr float PI = 3.14159265358979323846264338327950288f;
constexpr float SQRT_TWO = 1.41421356237309504880f;
constexpr float EPSILON = 1.192092896e-07f;
constexpr float FLOAT_MAX = 3.402823466e+38f;
constexpr float RAD_TO_DEG = 180.f / PI;
constexpr float DEG_TO_RAD = PI / 180.f;

// Core
using glm::abs;
using glm::ceil;
using glm::clamp;
using glm::exp;
using glm::floor;
using glm::fract;
using glm::inversesqrt;
using glm::max;
using glm::min;
using glm::mod;
using glm::modf;
using glm::pow;
using glm::round;
using glm::sign;
using glm::sqrt;
using glm::step;

// Trigonometry
using glm::acos;
using glm::acosh;
using glm::asin;
using glm::asinh;
using glm::atan;
using glm::atanh;
using glm::cos;
using glm::radians;
using glm::sin;
using glm::tan;
using glm::tanh;

// Vector
using glm::cross;
using glm::distance;
using glm::dot;
using glm::length;
template <typename T>
inline float length2(const T& v) {
    return dot(v, v);
}
using glm::distance;
template <typename T>
inline float distance2(const T& a, const T& b) {
    T c = a - b;
    return dot(c, c);
}
using glm::normalize;

template <int N, typename T, glm::qualifier Q>
T angle(glm::vec<N, T, Q> const& a, glm::vec<N, T, Q> const& b) {
    return acos(dot(normalize(a), normalize(b)));
}

template <int N, typename T, glm::qualifier Q>
T angle(glm::vec<N, T, Q> const& a, glm::vec<N, T, Q> const& b, glm::vec<N, T, Q> const& c) {
    return angle(a - b, c - b);
}

inline float dihedral_angle(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3) {
    vec3 b1 = p1 - p0;
    vec3 b2 = p2 - p1;
    vec3 b3 = p3 - p2;
    vec3 c1 = math::cross(b1, b2);
    vec3 c2 = math::cross(b2, b3);
    return glm::atan(glm::dot(glm::cross(c1, c2), glm::normalize(b2)), glm::dot(c1, c2));
}

inline float dihedral_angle(const vec3 p[4]) { return dihedral_angle(p[0], p[1], p[2], p[3]); }

// Matrix
using glm::determinant;
using glm::inverse;
using glm::transpose;

// Quat
using glm::conjugate;

// Cast
using glm::mat3_cast;
using glm::mat4_cast;
using glm::quat_cast;

// Interpolation
template <typename T, typename V>
T lerp(T const& a, T const& b, V t) {
    return ((V)1.0 - t) * a + t * b;
}

template <typename T, typename V>
T mix(T const& a, T const& b, V t) {
    return lerp(a, b, t);
}

template <typename T, typename V>
T catmull_rom(T const& v1, T const& v2, T const& v3, T const& v4, V s) {
    return glm::catmullRom(v1, v2, v3, v4, s);
}

template <typename T, typename V>
T cubic(T const& v1, T const& v2, T const& v3, T const& v4, V s) {
    return glm::cubic(v1, v2, v3, v4, s);
}

template <typename T, typename V>
T hermite(T const& v1, T const& t1, T const& v2, T const& t2, V s) {
    return glm::cubic(v1, t1, v2, t2, s);
}

template <typename T, typename V>
static T spline(const T& p0, const T& p1, const T& p2, const T& p3, V s, V tension = (V)0.5) {
    T v0 = (p2 - p0) * tension;
    T v1 = (p3 - p1) * tension;
    V s2 = s * s;
    V s3 = s * s2;

    return ((V)2.0 * p1 - (V)2.0 * p2 + v0 + v1) * s3 + (-(V)3.0 * p1 + (V)3.0 * p2 - (V)2.0 * v0 - v1) * s2 + v0 * s + p1;
}

template <typename T, typename V>
static T spline_tangent(const T& p0, const T& p1, const T& p2, const T& p3, V s, V tension = (V)0.5) {
	T v0 = (p2 - p0) * tension;
	T v1 = (p3 - p1) * tension;

	// f(t) = (2p1 - 2p2 + v0 + v1)s^3 + (-3p1 + 3p2 - 2v0 - v1)s^2 + v0s + p1;
	// df(t)/dt = (2p1 - 2p2 + v0 + v1)*3s^2 + (-3p1 + 3p2 - 2v0 - v1)*2s + v0;
	return ((V)2.0 * p1 - (V)2.0 * p2 + v0 + v1) * (V)3.0 * s * s + (-(V)3.0 * p1 + (V)3.0 * p2 - (V)2.0 * v0 - v1) * (V)2.0 * s + v0;
}

static inline glm_vec4 spline(const glm_vec4 p0, const glm_vec4 p1, const glm_vec4 p2, const glm_vec4 p3, float s, float tension = 0.5f) {
    const glm_vec4 v_t = _mm_set_ps1(tension);

    const glm_vec4 v0 = glm_vec4_mul(glm_vec4_sub(p2, p0), v_t);
    const glm_vec4 v1 = glm_vec4_mul(glm_vec4_sub(p3, p1), v_t);

    // (2p1  - 2p2 + v0 + v1) s3 + (-3p1 + 3p2 - 2v0 - v1) s2 + v0 s + p1
    // a = 2(p1 - p2) + (v0 + v1)
    // b = 3(p2 - p1) - (2v0 + v1) 

    const glm_vec4 a = glm_vec4_add(glm_vec4_mul(_mm_set_ps1(2), glm_vec4_sub(p1, p2)), glm_vec4_add(v0, v1));
    const glm_vec4 b = glm_vec4_sub(glm_vec4_mul(_mm_set_ps1(3), glm_vec4_sub(p2, p1)), glm_vec4_add(glm_vec4_mul(_mm_set_ps1(2), v0), v1));
    
    const glm_vec4 res_0 = glm_vec4_add(glm_vec4_mul(a, _mm_set_ps1(s*s*s)), glm_vec4_mul(b, _mm_set_ps1(s*s)));
    const glm_vec4 res_1 = glm_vec4_add(glm_vec4_mul(v0, _mm_set_ps1(s)), p1);
    return glm_vec4_add(res_0, res_1);
}

// Quaternion
template <typename T, glm::qualifier Q>
glm::tquat<T, Q> angle_axis(T const& angle, glm::vec<3, T, Q> const& axis) {
    return glm::angleAxis(angle, axis);
}

// Random
inline float rnd() { return (float)rand() / (float)RAND_MAX; }
inline void set_rnd_seed(unsigned int seed) { srand(seed); }

// Color

// http://lolengine.net/blog/2013/07/27/rgb-to-hsv-in-glsl

inline vec3 rgb_to_hsv(vec3 c) {
    vec4 K = vec4(0.0f, -1.0f / 3.0f, 2.0f / 3.0f, -1.0f);
    vec4 p = mix(vec4(c.b, c.g, K.w, K.z), vec4(c.g, c.b, K.x, K.y), step(c.b, c.g));
    vec4 q = mix(vec4(p.x, p.y, p.w, c.r), vec4(c.r, p.y, p.z, p.x), step(p.x, c.r));

    float d = q.x - min(q.w, q.y);
    float e = 1.0e-10f;
    return vec3(abs(q.z + (q.w - q.y) / (6.0f * d + e)), d / (q.x + e), q.x);
}

inline vec3 hsv_to_rgb(vec3 c) {
    vec4 K = vec4(1.0f, 2.0f / 3.0f, 1.0f / 3.0f, 3.0f);
    vec3 p = abs(fract(vec3(c.x) + vec3(K)) * 6.0f - vec3(K.w));
    return c.z * mix(vec3(K.x), clamp(p - vec3(K.x), 0.0f, 1.0f), c.y);
}

constexpr float HCLgamma = 3;
constexpr float HCLy0 = 100;
constexpr float HCLmaxL = 0.530454533953517f;  // == exp(HCLgamma / HCLy0) - 0.5

inline vec3 hcl_to_rgb(vec3 HCL) {
    vec3 RGB = vec3(0);
    if (HCL.z != 0) {
        float H = HCL.x;
        float C = HCL.y;
        float L = HCL.z * HCLmaxL;
        float Q = exp((1 - C / (2 * L)) * (HCLgamma / HCLy0));
        float U = (2 * L - C) / (2 * Q - 1);
        float V = C / Q;
        float T = tan((H + min(fract(2 * H) / 4.f, fract(-2 * H) / 8.f)) * PI * 2);
        H *= 6;
        if (H <= 1) {
            RGB.r = 1;
            RGB.g = T / (1 + T);
        } else if (H <= 2) {
            RGB.r = (1 + T) / T;
            RGB.g = 1;
        } else if (H <= 3) {
            RGB.g = 1;
            RGB.b = 1 + T;
        } else if (H <= 4) {
            RGB.g = 1 / (1 + T);
            RGB.b = 1;
        } else if (H <= 5) {
            RGB.r = -1 / T;
            RGB.b = 1;
        } else {
            RGB.r = 1;
            RGB.b = -T;
        }
        RGB = RGB * V + U;
    }
    return RGB;
}

inline vec3 rgb_to_hcl(vec3 rgb) {
    vec3 HCL;
    float H = 0;
    float U = min(rgb.r, min(rgb.g, rgb.b));
    float V = max(rgb.r, max(rgb.g, rgb.b));
    float Q = HCLgamma / HCLy0;
    HCL.y = V - U;
    if (HCL.y != 0) {
        H = atan(rgb.g - rgb.b, rgb.r - rgb.g) / PI;
        Q *= U / V;
    }
    Q = exp(Q);
    HCL.x = fract(H / 2.f - min(fract(H), fract(-H)) / 6.f);
    HCL.y *= Q;
    HCL.z = lerp(-U, V, Q) / (HCLmaxL * 2);
    return HCL;
}

inline vec4 convert_color(uint32 color) { return glm::unpackUnorm4x8(color); }
inline uint32 convert_color(vec4 color) { return glm::packUnorm4x8(color); }

}  // namespace math
