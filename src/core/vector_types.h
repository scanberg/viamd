#pragma once

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <glm/gtc/quaternion.hpp>
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
using glm::vec2;
using glm::vec3;
using glm::vec4;

using f32vec2 = glm::tvec2<float>;
using f32vec3 = glm::tvec3<float>;
using f32vec4 = glm::tvec4<float>;

using glm::dvec2;
using glm::dvec3;
using glm::dvec4;

using f64vec2 = glm::tvec2<double>;
using f64vec3 = glm::tvec3<double>;
using f64vec4 = glm::tvec4<double>;

using glm::ivec2;
using glm::ivec3;
using glm::ivec4;

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

using glm::uvec2;
using glm::uvec3;
using glm::uvec4;

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

using glm::mat2;
using glm::mat3;
using glm::mat4;

using glm::mat2x2;
using glm::mat2x3;
using glm::mat2x4;

using glm::mat3x2;
using glm::mat3x3;
using glm::mat3x4;

using glm::mat4x2;
using glm::mat4x3;
using glm::mat4x4;

using Range = vec2;
