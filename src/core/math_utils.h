#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/random.hpp>
#include <float.h>

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

	// Vector algebra
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

	using glm::mat3_cast;
	using glm::mat4_cast;
	using glm::quat_cast;

	template<typename T, glm::qualifier Q>
	glm::tquat<T, Q> angle_axis(T const& angle, glm::vec<3, T, Q> const& v) { return glm::angleAxis(angle, v); }

	// Random
	inline float rnd() { return glm::linearRand(0.f, 1.f); }

}