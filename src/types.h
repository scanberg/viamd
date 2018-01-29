#pragma once

#include <stdint.h>

#define Kilobytes(x) (x << 10)
#define Megabytes(x) (Kilobytes(x) << 10)
#define Gigabytes(x) (Megabytes(x) << 10)
#define UNUSED(x) (void)(x)

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

struct f32vec2;
struct f32vec3;
struct f32vec4;

struct i32vec2;
struct i32vec3;
struct i32vec4;
