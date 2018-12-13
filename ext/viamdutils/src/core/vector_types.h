#pragma once

//#ifndef GLM_FORCE_SSE2
//#define GLM_FORCE_SSE2
//#endif

/*
#ifndef GLM_FORCE_ALIGNED
#define GLM_FORCE_ALIGNED
#endif
*/

#define GLM_FORCE_DEFAULT_ALIGNED_GENTYPES

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/mat2x2.hpp>
#include <glm/mat2x3.hpp>
#include <glm/mat2x4.hpp>
#include <glm/mat3x2.hpp>
#include <glm/mat3x3.hpp>
#include <glm/mat3x4.hpp>
#include <glm/mat4x2.hpp>
#include <glm/mat4x3.hpp>
#include <glm/mat4x4.hpp>

using quat = glm::tquat<float>;
using vec2 = glm::tvec2<float>;
using vec3 = glm::tvec3<float>;
using vec4 = glm::tvec4<float>;

using f32vec2 = glm::tvec2<float>;
using f32vec3 = glm::tvec3<float>;
using f32vec4 = glm::tvec4<float>;

using dvec2 = glm::tvec2<double>;
using dvec3 = glm::tvec3<double>;
using dvec4 = glm::tvec4<double>;

using f64vec2 = glm::tvec2<double>;
using f64vec3 = glm::tvec3<double>;
using f64vec4 = glm::tvec4<double>;

using ivec2 = glm::tvec2<int32_t>;
using ivec3 = glm::tvec3<int32_t>;
using ivec4 = glm::tvec4<int32_t>;

using i8vec2 = glm::tvec2<int8_t>;
using i8vec3 = glm::tvec3<int8_t>;
using i8vec4 = glm::tvec4<int8_t>;

using i16vec2 = glm::tvec2<int16_t>;
using i16vec3 = glm::tvec3<int16_t>;
using i16vec4 = glm::tvec4<int16_t>;

using i32vec2 = glm::tvec2<int32_t>;
using i32vec3 = glm::tvec3<int32_t>;
using i32vec4 = glm::tvec4<int32_t>;

using i64vec2 = glm::tvec2<int64_t>;
using i64vec3 = glm::tvec3<int64_t>;
using i64vec4 = glm::tvec4<int64_t>;

using uvec2 = glm::tvec2<uint32_t>;
using uvec3 = glm::tvec3<uint32_t>;
using uvec4 = glm::tvec4<uint32_t>;

using u8vec2 = glm::tvec2<uint8_t>;
using u8vec3 = glm::tvec3<uint8_t>;
using u8vec4 = glm::tvec4<uint8_t>;

using u16vec2 = glm::tvec2<uint16_t>;
using u16vec3 = glm::tvec3<uint16_t>;
using u16vec4 = glm::tvec4<uint16_t>;

using u32vec2 = glm::tvec2<uint32_t>;
using u32vec3 = glm::tvec3<uint32_t>;
using u32vec4 = glm::tvec4<uint32_t>;

using u64vec2 = glm::tvec2<uint64_t>;
using u64vec3 = glm::tvec3<uint64_t>;
using u64vec4 = glm::tvec4<uint64_t>;

using mat2 = glm::tmat2x2<float>;
using mat3 = glm::tmat3x3<float>;
using mat4 = glm::tmat4x4<float>;
using mat2x2 = glm::tmat2x2<float>;
using mat2x3 = glm::tmat2x3<float>;
using mat2x4 = glm::tmat2x4<float>;
using mat3x2 = glm::tmat3x2<float>;
using mat3x3 = glm::tmat3x3<float>;
using mat3x4 = glm::tmat3x4<float>;
using mat4x2 = glm::tmat4x2<float>;
using mat4x3 = glm::tmat4x3<float>;
using mat4x4 = glm::tmat4x4<float>;

using Range = vec2;
