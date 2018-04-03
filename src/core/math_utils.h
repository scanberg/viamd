#pragma once

#include <core/types.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/random.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/spline.hpp>
#include <float.h>
#include <stdlib.h>

namespace math {
constexpr float PI = glm::pi<float>();
constexpr float EPSILON = glm::epsilon<float>();
constexpr float FLOAT_MAX = FLT_MAX;
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

// Matrix
using glm::determinant;
using glm::inverse;
using glm::transpose;

// Cast
using glm::mat3_cast;
using glm::mat4_cast;
using glm::quat_cast;

// Interpolation
using glm::mix;
template <typename T, typename V>
T lerp(T const& a, T const& b, V t) {
    return glm::mix(a, b, t);
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

// Quaternion
template <typename T, glm::qualifier Q>
glm::tquat<T, Q> angle_axis(T const& angle, glm::vec<3, T, Q> const& axis) {
    return glm::angleAxis(angle, axis);
}

// Random
inline float rnd() { return rand() / (float)RAND_MAX; }
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

}  // namespace math