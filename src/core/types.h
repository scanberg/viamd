#pragma once

#include <stdint.h>
//#include <glm/fwd.hpp>
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

typedef int8_t int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;

typedef uint8_t uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

typedef float float32;
typedef double float64;

using glm::vec2;
using glm::vec3;
using glm::vec4;
using glm::quat;

using glm::f32vec2;
using glm::f32vec3;
using glm::f32vec4;

using glm::dvec2;
using glm::dvec3;
using glm::dvec4;

using glm::f64vec2;
using glm::f64vec3;
using glm::f64vec4;

using glm::ivec2;
using glm::ivec3;
using glm::ivec4;

using glm::i8vec2;
using glm::i8vec3;
using glm::i8vec4;

using glm::i16vec2;
using glm::i16vec3;
using glm::i16vec4;

using glm::i32vec2;
using glm::i32vec3;
using glm::i32vec4;

using glm::i64vec2;
using glm::i64vec3;
using glm::i64vec4;

using glm::uvec2;
using glm::uvec3;
using glm::uvec4;

using glm::u8vec2;
using glm::u8vec3;
using glm::u8vec4;

using glm::u16vec2;
using glm::u16vec3;
using glm::u16vec4;

using glm::u32vec2;
using glm::u32vec3;
using glm::u32vec4;

using glm::u64vec2;
using glm::u64vec3;
using glm::u64vec4;

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

//  @TODO: Someday remove dependency on GLM and implement proper SIMD compatible vector types (No templates) using the principles bellow
/*
#define VS

#if defined(CLANG_OR_GCC)
typedef float float2 __attribute__((ext_vector_type(2)));
typedef float float3 __attribute__((ext_vector_type(3)));
typedef float float4 __attribute__((ext_vector_type(4)));
#elif defined(VS)
struct float2 {
	float x, y;
	// xx
	__declspec(property(get = get_xx, put = put_xx)) float2 xx;
	inline float2 get_xx() const { return { x, x }; }
	inline void put_xx(float2 v) { *this = v.xx; }
	// xy
	__declspec(property(get = get_xy, put = put_xy)) float2 xy;
	inline float2 get_xy() const { return { x, y }; }
	inline void put_xy(float2 v) { *this = v.xy; }
	// yx
	__declspec(property(get = get_yx, put = put_yx)) float2 yx;
	inline float2 get_yx() const { return { y, x }; }
	inline void put_yx(float2 v) { *this = v.yx; }
	// yy
	__declspec(property(get = get_yy, put = put_yy)) float2 yy;
	inline float2 get_yy() const { return { y, y }; }
	inline void put_yy(float2 v) { *this = v.yy; }
};
struct float3 {
	float x, y, z;
	// xx
	__declspec(property(get = get_xx, put = put_xx)) float2 xx;
	inline float2 get_xx() const { return { x, x }; }
	inline void put_xx(float2 v) { this->x = v.x; this->y = v.x; }
	// xy
	__declspec(property(get = get_xy, put = put_xy)) float2 xy;
	inline float2 get_xy() const { return { x, y }; }
	inline void put_xy(float2 v) { this->x = v.x; this->y = v.y; }
	// yx
	__declspec(property(get = get_yx, put = put_yx)) float2 yx;
	inline float2 get_yx() const { return { y, x }; }
	inline void put_yx(float2 v) { this->x = v.y; this->y = v.x; }
	// yy
	__declspec(property(get = get_yy, put = put_yy)) float2 yy;
	inline float2 get_yy() const { return { y, y }; }
	inline void put_yy(float2 v) { this->x = v.y; this->y = v.y; }

	// xxx
	__declspec(property(get = get_xxx, put = put_xxx)) float3 xxx;
	inline float3 get_xxx() const { return { x, x, x }; }
	inline void put_xxx(float3 v) { *this = {v.x, v.x, v.x}; }
	// xxy
	__declspec(property(get = get_xxy, put = put_xxy)) float3 xxy;
	inline float3 get_xxy() const { return { x, x, y }; }
	inline void put_xxy(float3 v) { *this = { v.x, v.x, v.y }; }
	// xxz
	__declspec(property(get = get_xxz, put = put_xxz)) float3 xxz;
	inline float3 get_xxz() const { return { x, x, z }; }
	inline void put_xxz(float3 v) { *this = { v.x, v.x, v.z }; }

	// xyx
	__declspec(property(get = get_xyx, put = put_xyx)) float3 xyx;
	inline float3 get_xyx() const { return { x, y, x }; }
	inline void put_xyx(float3 v) { *this = { v.x, v.y, v.x }; }
	// xyy
	__declspec(property(get = get_xyy, put = put_xyy)) float3 xyy;
	inline float3 get_xyy() const { return { x, y, y }; }
	inline void put_xyy(float3 v) { *this = { v.x, v.y, v.y }; }
	// xyz
	__declspec(property(get = get_xyz, put = put_xyz)) float3 xyz;
	inline float3 get_xyz() const { return { x, y, z }; }
	inline void put_xyz(float3 v) { *this = { v.x, v.y, v.z }; }

	// xzx
	__declspec(property(get = get_xzx, put = put_xzx)) float3 xzx;
	inline float3 get_xzx() const { return { x, z, x }; }
	inline void put_xzx(float3 v) { *this = { v.x, v.z, v.x }; }
	// xzy
	__declspec(property(get = get_xzy, put = put_xzy)) float3 xzy;
	inline float3 get_xzy() const { return { x, z, y }; }
	inline void put_xzy(float3 v) { *this = { v.x, v.z, v.y }; }
	// xzz
	__declspec(property(get = get_xzz, put = put_xzz)) float3 xzz;
	inline float3 get_xzz() const { return { x, z, z }; }
	inline void put_xzz(float3 v) { *this = { v.x, v.z, v.z }; }

	// yxx
	__declspec(property(get = get_yxx, put = put_yxx)) float3 yxx;
	inline float3 get_yxx() const { return { y, x, x }; }
	inline void put_yxx(float3 v) { *this = { v.y, v.x, v.x }; }
	// yxy
	__declspec(property(get = get_yxy, put = put_yxy)) float3 yxy;
	inline float3 get_yxy() const { return { y, x, y }; }
	inline void put_yxy(float3 v) { *this = { v.y, v.x, v.y }; }
	// yxz
	__declspec(property(get = get_yxz, put = put_yxz)) float3 yxz;
	inline float3 get_yxz() const { return { y, x, z }; }
	inline void put_yxz(float3 v) { *this = { v.y, v.x, v.z }; }

	// yyx
	__declspec(property(get = get_yyx, put = put_yyx)) float3 yyx;
	inline float3 get_yyx() const { return { y, y, x }; }
	inline void put_yyx(float3 v) { *this = { v.y, v.y, v.x }; }
	// yyy
	__declspec(property(get = get_yyy, put = put_yyy)) float3 yyy;
	inline float3 get_yyy() const { return { y, y, y }; }
	inline void put_yyy(float3 v) { *this = { v.y, v.y, v.y }; }
	// yyz
	__declspec(property(get = get_yyz, put = put_yyz)) float3 yyz;
	inline float3 get_yyz() const { return { y, y, z }; }
	inline void put_yyz(float3 v) { *this = { v.y, v.y, v.z }; }

	// yzx
	__declspec(property(get = get_yzx, put = put_yzx)) float3 yzx;
	inline float3 get_yzx() const { return { y, z, x }; }
	inline void put_yzx(float3 v) { *this = { v.y, v.z, v.x }; }
	// yzy
	__declspec(property(get = get_yzy, put = put_yzy)) float3 yzy;
	inline float3 get_yzy() const { return { y, z, y }; }
	inline void put_yzy(float3 v) { *this = { v.y, v.z, v.y }; }
	// yzz
	__declspec(property(get = get_yzz, put = put_yzz)) float3 yzz;
	inline float3 get_yzz() const { return { y, z, z }; }
	inline void put_yzz(float3 v) { *this = { v.y, v.z, v.z }; }

	// zxx
	__declspec(property(get = get_zxx, put = put_zxx)) float3 zxx;
	inline float3 get_zxx() const { return { z, x, x }; }
	inline void put_zxx(float3 v) { *this = { v.z, v.x, v.x }; }
	// zxy
	__declspec(property(get = get_zxy, put = put_zxy)) float3 zxy;
	inline float3 get_zxy() const { return { z, x, y }; }
	inline void put_zxy(float3 v) { *this = { v.z, v.x, v.y }; }
	// zxz
	__declspec(property(get = get_zxz, put = put_zxz)) float3 zxz;
	inline float3 get_zxz() const { return { z, x, z }; }
	inline void put_zxz(float3 v) { *this = { v.z, v.x, v.z }; }

	// zyx
	__declspec(property(get = get_zyx, put = put_zyx)) float3 zyyx;
	inline float3 get_zyx() const { return { z, y, x }; }
	inline void put_zyx(float3 v) { *this = { v.z, v.y, v.x }; }
	// zyy
	__declspec(property(get = get_zyy, put = put_zyy)) float3 zyy;
	inline float3 get_zyy() const { return { z, y, y }; }
	inline void put_zyy(float3 v) { *this = { v.z, v.y, v.y }; }
	// zyz
	__declspec(property(get = get_zyz, put = put_zyz)) float3 zyz;
	inline float3 get_zyz() const { return { z, y, z }; }
	inline void put_zyz(float3 v) { *this = { v.z, v.y, v.z }; }

	// zzx
	__declspec(property(get = get_zzx, put = put_zzx)) float3 zzx;
	inline float3 get_zzx() const { return { z, z, x }; }
	inline void put_zzx(float3 v) { *this = { v.z, v.z, v.x }; }
	// zzy
	__declspec(property(get = get_zzy, put = put_yzy)) float3 zzy;
	inline float3 get_zzy() const { return { z, z, y }; }
	inline void put_zzy(float3 v) { *this = { v.z, v.z, v.y }; }
	// zzz
	__declspec(property(get = get_zzz, put = put_zzz)) float3 zzz;
	inline float3 get_zzz() const { return { z, z, z }; }
	inline void put_zzz(float3 v) { *this = { v.z, v.z, v.z }; }
};

struct float4 {
	float x, y, z, w;

	// TODO: Swizzle permutations 
};
#else
	#error
#endif

struct float3x3 {
	float3 col[3];
};

struct float4x4 {
	float4 col[4];
};
*/