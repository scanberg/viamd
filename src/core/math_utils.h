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
using glm::min;
using glm::max;
using glm::clamp;
using glm::abs;
using glm::floor;
using glm::ceil;
using glm::fract;
using glm::sign;
using glm::mod;
using glm::modf;
using glm::step;
using glm::exp;
using glm::pow;
using glm::sqrt;
using glm::inversesqrt;

// Trigonometry
using glm::cos;
using glm::acos;
using glm::acosh;
using glm::sin;
using glm::asin;
using glm::asinh;
using glm::tan;
using glm::tanh;
using glm::atan;
using glm::atanh;
using glm::radians;

// Vector
using glm::dot;
using glm::cross;
using glm::distance;
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

// Matrix
using glm::transpose;
using glm::inverse;
using glm::determinant;

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

// Quaternion
template <typename T, glm::qualifier Q>
glm::tquat<T, Q> angle_axis(T const& angle, glm::vec<3, T, Q> const& axis) {
    return glm::angleAxis(angle, axis);
}

// Random
inline float rnd() { return rand() / (float)RAND_MAX; }
inline void set_rnd_seed(unsigned int seed) { srand(seed); }

// Source from here
// https://github.com/gegaryfa/Color-conversions

inline vec3 rgb_to_xyz(vec3 rgb) {
    float r = rgb[0];
    float g = rgb[1];
    float b = rgb[2];

    if (r > 0.04045f)
        r = pow(((r + 0.055f) / 1.055f), 2.4f);
    else
        r = r / 12.92f;

    if (g > 0.04045f)
        g = pow(((g + 0.055f) / 1.055f), 2.4f);
    else
        g = g / 12.92f;

    if (b > 0.04045f)
        b = pow(((b + 0.055f) / 1.055f), 2.4f);
    else
        b = b / 12.92f;

    // Scale or don't scale with 100?
    r = r * 100;
    g = g * 100;
    b = b * 100;

    // Observer. = 2°, Illuminant = D65
    float x = r * 0.4124f + g * 0.3576f + b * 0.1805f;
    float y = r * 0.2126f + g * 0.7152f + b * 0.0722f;
    float z = r * 0.0193f + g * 0.1192f + b * 0.9505f;

    return {x, y, z};
}

inline vec3 xyz_to_rgb(vec3 xyz) {
    float x = xyz[0] / 100.f;
    float y = xyz[1] / 100.f;
    float z = xyz[2] / 100.f;

    float r = x * 3.2406f + (y * -1.5372f) + z * (-0.4986f);
    float g = x * (-0.9689f) + y * 1.8758f + z * 0.0415f;
    float b = x * 0.0557f + y * (-0.2040f) + z * 1.0570f;

    if (r > 0.0031308f)
        r = 1.055f * powf(r, (1.0f / 2.4f)) - 0.055f;
    else
        r = 12.92f * r;

    if (g > 0.0031308f)
        g = 1.055f * powf(g, (1.0f / 2.4f)) - 0.055f;
    else
        g = 12.92f * g;

    if (b > 0.0031308f)
        b = 1.055f * powf(b, (1.0f / 2.4f)) - 0.055f;
    else
        b = 12.92f * b;

    return {r, g, b};
}

// D65 Reference white point in CIEXYZ
constexpr float ref_X = 95.047f;
constexpr float ref_Y = 100.0f;
constexpr float ref_Z = 108.883f;

inline vec3 xyz_to_lab(vec3 xyz) {
    float x = (xyz[0] / ref_X);  // ref_X = 95.047
    float y = (xyz[1] / ref_Y);  // ref_Y = 100.0
    float z = (xyz[2] / ref_Z);  // ref_Z = 108.883

    if (x > 0.008856f)
        x = pow(x, (1.f / 3.f));
    else
        x = (7.787f * x) + (16.f / 116.f);

    if (y > 0.008856f)
        y = pow(y, (1.f / 3.f));
    else
        y = (7.787f * y) + (16.f / 116.f);

    if (z > 0.008856f)
        z = pow(z, (1.f / 3.f));
    else
        z = (7.787f * z) + (16.f / 116.f);

    float L = (116.f * y) - 16.f;
    float a = 500.f * (x - y);
    float b = 200.f * (y - z);

    return {L, a, b};
}

inline vec3 lab_to_xyz(vec3 Lab) {
    float y = (Lab[0] + 16.f) / 116.f;
    float x = (Lab[1] / 500.f) + y;
    float z = y - (Lab[2] / 200.f);

    if (powf(y, 3.f) > 0.008856f)
        y = powf(y, 3.f);
    else
        y = (y - (16.f / 116.f)) / 7.787f;

    if (powf(x, 3.f) > 0.008856f)
        x = powf(x, 3.f);
    else
        x = (x - (16.f / 116.f)) / 7.787f;

    if (powf(z, 3.f) > 0.008856f)
        z = powf(z, 3.f);
    else
        z = (z - (16.f / 116.f)) / 7.787f;

    x = ref_X * x;  // ref_X =  95.047     Observer= 2°, Illuminant= D65
    y = ref_Y * x;  // ref_Y = 100.000
    z = ref_Z * x;  // ref_Z = 108.883

    return {x, y, z};
}

inline vec3 rgb_to_lab(vec3 rgb) {
	return xyz_to_lab(rgb_to_xyz(rgb));
}

inline vec3 lab_to_rgb(vec3 lab) {
	return xyz_to_rgb(lab_to_xyz(lab));
}

inline vec3 lab_to_lch(vec3 lab) {
	float l = lab[0];
	float a = lab[1];
	float b = lab[2];
    float c = sqrt(a * a + b * b);
    float h = mod(atan(b, a) + 2.f * PI, 2.f * PI);
    return {l, c, h};
};

inline vec3 lch_to_lab(vec3 lch) {
	float l = lch[0];
	float c = lch[1];
	float h = lch[2];
	float a = c * cos(h);
	float b = c * sin(h);
	return {l, a, b};
}

inline float lab_color_difference(vec3 lab1, vec3 lab2) {
    return distance(lab1, lab2);
}

inline float rgb_color_lab_difference(vec3 rgb1, vec3 rgb2) {
	return lab_color_difference(rgb_to_lab(rgb1), rgb_to_lab(rgb2));
}

}