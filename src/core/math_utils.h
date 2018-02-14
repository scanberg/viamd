#pragma once

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

	// Core
	using glm::min;
	using glm::max;
	using glm::clamp;
	using glm::abs;
	using glm::floor;
	using glm::ceil;
	using glm::fract;
	using glm::sign;
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

	// Casts
	using glm::mat3_cast;
	using glm::mat4_cast;
	using glm::quat_cast;

	// Interpolation
	using glm::mix;
	template<typename T, typename V>
	T lerp(T const& a, T const& b, V t) { return glm::mix(a,b,t); }

	template<typename T, typename V>
	T catmull_rom(T const& v1, T const& v2, T const& v3, T const& v4, V s) { return glm::catmullRom(v1, v2, v3, v4, s); }

	template<typename T, typename V>
	T cubic(T const& v1, T const& v2, T const& v3, T const& v4, V s) { return glm::cubic(v1, v2, v3, v4, s); }

	template<typename T, typename V>
	T hermite(T const& v1, T const& t1, T const& v2, T const& t2, V s) { return glm::cubic(v1, v2, v3, v4, s); }

	// Quaternion
	template<typename T, glm::qualifier Q>
	glm::tquat<T, Q> angle_axis(T const& angle, glm::vec<3, T, Q> const& axis) { return glm::angleAxis(angle, axis); }

	// Random
	inline float rnd() { return rand() / (float)RAND_MAX; }
	inline void set_rnd_seed(unsigned int seed) { srand(seed); }

}