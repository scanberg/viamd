#pragma once

#define VIAMD_USE_SIMD 1

#include <stdint.h>
#ifdef WIN32
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

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

struct F32Vec2 {
    union {
        struct { float32 x,y; };
        struct { float32 r,g; };
        struct { float32 s,t; };
    };

    F32Vec2& operator += (float32 v) {
        x += v; y += v;
    }

    F32Vec2& operator -= (float32 v) {
        x -= v; y -= v;
    }

    F32Vec2& operator *= (float32 v) {
        x *= v; y *= v;
        return *this;
    }

    F32Vec2& operator /= (float32 v) {
        x /= v; y /= v;
        return *this;
    }

    float32& operator[](int i) { return (&x) + i; }
};

struct F32Vec3 {
    union {
        struct { float32 x,y,z; };
        struct { float32 r,g,b; };
        struct { float32 s,t,u; };
#ifdef VIAMD_USE_SIMD
        __m128 simd;
#endif
    };

    F32Vec3& operator += (float32 v) {
        x += v; y += v; z += v;
    }

    F32Vec3& operator -= (float32 v) {
        x -= v; y -= v; z -= v;
    }

    F32Vec3& operator *= (float32 v) {
        x *= v; y *= v; z *= v;
        return *this;
    }

    F32Vec3& operator /= (float32 v) {
        x /= v; y /= v; z /= v;
        return *this;
    }
};

struct F32Vec4 {
    union {
        struct { float32 x,y,z,w; };
        struct { float32 r,g,b,a; };
        struct { float32 s,t,u,v; };
#ifdef VIAMD_USE_SIMD
        __m128 simd;
#endif
    };

    F32Vec4& operator += (float32 v) {
        x += v; y += v; z += v; w += v;
    }

    F32Vec4& operator -= (float32 v) {
        x -= v; y -= v; z -= v; w -= v;
    }

    F32Vec4& operator *= (float32 v) {
        x *= v; y *= v; z *= v; w *= v;
        return *this;
    }
    F32Vec4& operator /= (float32 v) {
        x /= v; y /= v; z /= v; w *= v;
        return *this;
    }
};

struct I32Vec2 {
    union {
        struct { int32 x,y; };
        struct { int32 r,g; };
        struct { int32 s,t; };
    };

    I32Vec2& operator += (int32 v) {
        x += v; y += v;
        return *this;
    }

    I32Vec2& operator -= (int32 v) {
        x -= v; y -= v;
        return *this;
    }

    I32Vec2& operator *= (int32 v) {
        x *= v; y *= v;
        return *this;
    }

    I32Vec2& operator /= (int32 v) {
        x /= v; y /= v;
        return *this;
    }
};

struct I32Vec3 {
    union {
        struct { int32 x,y,z; };
        struct { int32 r,g,b; };
        struct { int32 s,t,u; };
#ifdef VIAMD_USE_SIMD
        __m128i simd;
#endif
    };

    I32Vec3& operator += (int32 v) {
        x += v; y += v; z += v;
        return *this;
    }

    I32Vec3& operator -= (int32 v) {
        x -= v; y -= v; z -= v;
        return *this;
    }

    I32Vec3& operator *= (int32 v) {
        x *= v; y *= v; z *= v;
        return *this;
    }

    I32Vec3& operator /= (int32 v) {
        x /= v; y /= v; z /= v;
        return *this;
    }
};

struct I32Vec4 {
    union {
        struct { int32 x,y,z,w; };
        struct { int32 r,g,b,a; };
        struct { int32 s,t,u,v; };
#ifdef VIAMD_USE_SIMD
        __m128i simd;
#endif
    };

    I32Vec4& operator += (int32 v) {
        x += v; y += v; z += v; w += v;
        return *this;
    }

    I32Vec4& operator -= (int32 v) {
        x -= v; y -= v; z -= v; w += v;
        return *this;
    }

    I32Vec4& operator *= (int32 v) {
        x *= v; y *= v; z *= v; w *= v;
        return *this;
    }
    I32Vec4& operator /= (int32 v) {
        x /= v; y /= v; z /= v; w *= v;
        return *this;
    }
};

struct F32Mat2 {
    F32Vec2 col[2];
};


// using glm::vec2;
// using glm::vec3;
// using glm::vec4;
// using glm::quat;

// using glm::ivec2;
// using glm::ivec3;
// using glm::ivec4;

// using glm::uvec2;
// using glm::uvec3;
// using glm::uvec4;

// using glm::mat2;
// using glm::mat3;
// using glm::mat4;


