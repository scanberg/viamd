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

	using glm::dot;
	using glm::cross;
	using glm::abs;
	using glm::floor;
	using glm::ceil;
	using glm::clamp;
	using glm::exp;
	using glm::pow;
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
	using glm::distance;
	using glm::sqrt;
	using glm::inversesqrt;

	using glm::transpose;
	using glm::inverse;

	using glm::linearRand;
	using glm::sphericalRand;
	using glm::ballRand;
	inline float rnd() { return glm::linearRand(0.f, 1.f); }

}