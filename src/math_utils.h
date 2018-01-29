#pragma once

#define VIAMD_USE_SIMD 1

#include "types.h"

#ifdef WIN32
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

struct f32vec2 {
    f32vec2() = default;
    f32vec2(float32 val) : x(val), y(val) {}
    f32vec2(float32 _x, float32 _y) : x(_x), y(_y) {}

    union {
        struct {
            float32 x, y;
        };
        struct {
            float32 r, g;
        };
        struct {
            float32 s, t;
        };
    };

    f32vec2 &operator+=(float32 val) {
        x += val;
        y += val;
        return *this;
    }

    f32vec2 &operator-=(float32 val) {
        x -= val;
        y -= val;
        return *this;
    }

    f32vec2 &operator*=(float32 val) {
        x *= val;
        y *= val;
        return *this;
    }

    f32vec2 &operator/=(float32 val) {
        x /= val;
        y /= val;
        return *this;
    }

    float32 &operator[](int i) { return (&x)[i]; }
};

struct f32vec3 {
    f32vec3() = default;
    f32vec3(float32 val) : x(val), y(val), z(val) {}
    f32vec3(float32 _x, float32 _y, float32 _z) : x(_x), y(_y), z(_z) {}
    f32vec3(__m128 vec) : m128(vec) {}

    union {
        struct {
            float32 x, y, z;
        };
        struct {
            float32 r, g, b;
        };
        struct {
            float32 s, t, u;
        };
#ifdef VIAMD_USE_SIMD
        __m128 m128;
#endif
    };

    f32vec3 &operator+=(float32 val) {
        x += val;
        y += val;
        z += val;
        return *this;
    }

    f32vec3 &operator-=(float32 val) {
        x -= val;
        y -= val;
        z -= val;
        return *this;
    }

    f32vec3 &operator*=(float32 val) {
        x *= val;
        y *= val;
        z *= val;
        return *this;
    }

    f32vec3 &operator/=(float32 val) {
        x /= val;
        y /= val;
        z /= val;
        return *this;
    }

    float32 &operator[](int i) { return (&x)[i]; }
#ifdef VIAMD_USE_SIMD
    operator __m128() const { return m128; }
#endif
};

struct f32vec4 {
    f32vec4() = default;
    f32vec4(float32 val) : x(val), y(val), z(val), w(val) {}
    f32vec4(float32 _x, float32 _y, float32 _z, float32 _w) : x(_x), y(_y), z(_z), w(_w) {}
    f32vec4(__m128 vec) : m128(vec) {}

    union {
        struct {
            float32 x, y, z, w;
        };
        struct {
            float32 r, g, b, a;
        };
        struct {
            float32 s, t, u, v;
        };
#ifdef VIAMD_USE_SIMD
        __m128 m128;
#endif
    };

    f32vec4 &operator+=(float32 val) {
        x += val;
        y += val;
        z += val;
        w += val;
        return *this;
    }

    f32vec4 &operator-=(float32 val) {
        x -= val;
        y -= val;
        z -= val;
        w -= val;
        return *this;
    }

    f32vec4 &operator*=(float32 val) {
        x *= val;
        y *= val;
        z *= val;
        w *= val;
        return *this;
    }
    f32vec4 &operator/=(float32 val) {
        x /= val;
        y /= val;
        z /= val;
        w *= val;
        return *this;
    }

    float32 &operator[](int i) { return (&x)[i]; }
#ifdef VIAMD_USE_SIMD
    operator __m128() const { return m128; }
#endif
};

struct i32vec2 {
    union {
        struct {
            int32 x, y;
        };
        struct {
            int32 r, g;
        };
        struct {
            int32 s, t;
        };
    };

    i32vec2 &operator+=(int32 val) {
        x += val;
        y += val;
        return *this;
    }

    i32vec2 &operator-=(int32 val) {
        x -= val;
        y -= val;
        return *this;
    }

    i32vec2 &operator*=(int32 val) {
        x *= val;
        y *= val;
        return *this;
    }

    i32vec2 &operator/=(int32 val) {
        x /= val;
        y /= val;
        return *this;
    }
};

struct i32vec3 {
    union {
        struct {
            int32 x, y, z;
        };
        struct {
            int32 r, g, b;
        };
        struct {
            int32 s, t, u;
        };
#ifdef VIAMD_USE_SIMD
        __m128i m128i;
#endif
    };

    i32vec3 &operator+=(int32 val) {
        x += val;
        y += val;
        z += val;
        return *this;
    }

    i32vec3 &operator-=(int32 val) {
        x -= val;
        y -= val;
        z -= val;
        return *this;
    }

    i32vec3 &operator*=(int32 val) {
        x *= val;
        y *= val;
        z *= val;
        return *this;
    }

    i32vec3 &operator/=(int32 val) {
        x /= val;
        y /= val;
        z /= val;
        return *this;
    }
};

struct i32vec4 {
    union {
        struct {
            int32 x, y, z, w;
        };
        struct {
            int32 r, g, b, a;
        };
        struct {
            int32 s, t, u, v;
        };
#ifdef VIAMD_USE_SIMD
        __m128i m128i;
#endif
    };

    i32vec4 &operator+=(int32 val) {
        x += val;
        y += val;
        z += val;
        w += val;
        return *this;
    }

    i32vec4 &operator-=(int32 val) {
        x -= val;
        y -= val;
        z -= val;
        w += val;
        return *this;
    }

    i32vec4 &operator*=(int32 val) {
        x *= val;
        y *= val;
        z *= val;
        w *= val;
        return *this;
    }
    i32vec4 &operator/=(int32 val) {
        x /= val;
        y /= val;
        z /= val;
        w *= val;
        return *this;
    }
};

struct f32mat2 {
    f32mat2() = default;
	f32mat2(float32 s) {
		col[0] = { s,0 };
		col[1] = { 0,s };
	}
    f32mat2(f32vec2 c0, f32vec2 c1) {
        col[0] = c0;
        col[1] = c1;
    }
	f32mat2(
		float32 c00, float32 c01,
		float32 c10, float32 c11) {
		col[0] = { c00, c01 };
		col[1] = { c10, c11 };
	}
    f32vec2 col[2];

	f32vec2& operator[](int i) { return col[i]; }
};

struct f32mat3 {
	f32mat3() = default;
	f32mat3(float32 s) {
		col[0] = { s,0,0 };
		col[1] = { 0,s,0 };
		col[2] = { 0,s,0 };
	}
	f32mat3(f32vec3 c0, f32vec3 c1, f32vec3 c2) {
		col[0] = c0;
		col[1] = c1;
		col[2] = c2;
	}
	f32mat3(
		float32 c00, float32 c01, float32 c02,
		float32 c10, float32 c11, float32 c12,
		float32 c20, float32 c21, float32 c22) {
		col[0] = { c00, c01, c02 };
		col[1] = { c10, c11, c12 };
		col[2] = { c20, c21, c22 };
	}
	f32vec3 col[3];

	f32vec3& operator[](int i) { return col[i]; }
};

struct f32mat4 {
	f32mat4() = default;
	f32mat4(float32 s) {
		col[0] = { s,0,0,0 };
		col[1] = { 0,s,0,0 };
		col[2] = { 0,0,s,0 };
		col[3] = { 0,0,0,s };
	}
	f32mat4(f32vec4 c0, f32vec4 c1, f32vec4 c2, f32vec4 c3) {
		col[0] = c0;
		col[1] = c1;
		col[2] = c2;
		col[3] = c3;
	}
	f32mat4(
		float32 c00, float32 c01, float32 c02, float32 c03,
		float32 c10, float32 c11, float32 c12, float32 c13,
		float32 c20, float32 c21, float32 c22, float32 c23,
		float32 c30, float32 c31, float32 c32, float32 c33) {
		col[0] = { c00, c01, c02, c03 };
		col[1] = { c10, c11, c12, c13 };
		col[2] = { c20, c21, c22, c23 };
		col[3] = { c30, c31, c32, c33 };
	}
	f32vec4 col[3];

	f32vec4& operator[](int i) { return col[i]; }
};

namespace math {

// *** DOT PRODUCT ***
inline float32 dot(f32vec2 a, f32vec2 b) { return a.x * b.x + a.y * b.y; }

inline float32 dot(f32vec3 a, f32vec3 b) { return _mm_dp_ps(a.m128, b.m128, 0xff).m128_f32[0]; }

inline float32 dot(f32vec4 a, f32vec4 b) { return _mm_dp_ps(a.m128, b.m128, 0xff).m128_f32[0]; }

inline int32 dot(i32vec2 a, i32vec2 b) { return a.x * b.x + a.y * b.y; }

inline int32 dot(i32vec3 a, i32vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

inline int32 dot(i32vec4 a, i32vec4 b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }

// *** CROSS PRODUCT ***
// TODO: Possibly reverse order? Use right-hand
inline f32vec3 cross(f32vec3 a, f32vec3 b) {
    __m128 result = _mm_sub_ps(_mm_mul_ps(b.m128, _mm_shuffle_ps(a.m128, a.m128, _MM_SHUFFLE(3, 0, 2, 1))),
                               _mm_mul_ps(a.m128, _mm_shuffle_ps(b.m128, b.m128, _MM_SHUFFLE(3, 0, 2, 1))));
    return _mm_shuffle_ps(result, result, _MM_SHUFFLE(3, 0, 2, 1));
}

// OPERATORS FLOAT
inline f32vec2 operator+(f32vec2 a, f32vec2 b) { return {a.x + b.x, a.y + b.y}; }
inline f32vec3 operator+(f32vec3 a, f32vec3 b) { return _mm_add_ps(a, b); }
inline f32vec4 operator+(f32vec4 a, f32vec4 b) { return _mm_add_ps(a, b); }

inline f32vec2 operator-(f32vec2 a, f32vec2 b) { return {a.x - b.x, a.y - b.y}; }
inline f32vec3 operator-(f32vec3 a, f32vec3 b) { return _mm_sub_ps(a, b); }
inline f32vec4 operator-(f32vec4 a, f32vec4 b) { return _mm_sub_ps(a, b); }

inline f32vec2 operator*(f32vec2 a, f32vec2 b) { return {a.x * b.x, a.y * b.y}; }
inline f32vec3 operator*(f32vec3 a, f32vec3 b) { return _mm_mul_ps(a, b); }
inline f32vec4 operator*(f32vec4 a, f32vec4 b) { return _mm_mul_ps(a, b); }

inline f32vec2 operator*(f32vec2 v, float32 s) { return {v.x * s, v.y * s}; }
inline f32vec3 operator*(f32vec3 v, float32 s) { return _mm_mul_ps(v, _mm_set_ps1(s)); }
inline f32vec4 operator*(f32vec4 v, float32 s) { return _mm_mul_ps(v, _mm_set_ps1(s)); }

inline f32vec2 operator*(float32 s, f32vec2 v) { return v * s; }
inline f32vec3 operator*(float32 s, f32vec3 v) { return v * s; }
inline f32vec4 operator*(float32 s, f32vec4 v) { return v * s; }

inline f32vec2 operator/(f32vec2 a, f32vec2 b) { return {a.x / b.x, a.y / b.y}; }
inline f32vec3 operator/(f32vec3 a, f32vec3 b) { return _mm_div_ps(a, b); }
inline f32vec4 operator/(f32vec4 a, f32vec4 b) { return _mm_div_ps(a, b); }

inline f32vec2 operator/(f32vec2 v, float32 s) { return {v.x / s, v.y / s}; }
inline f32vec3 operator/(f32vec3 v, float32 s) { return _mm_div_ps(v, _mm_set_ps1(s)); }
inline f32vec4 operator/(f32vec4 v, float32 s) { return _mm_div_ps(v, _mm_set_ps1(s)); }

inline f32vec2 operator/(float32 s, f32vec2 v) { return {s / v.x, s / v.y}; }
inline f32vec3 operator/(float32 s, f32vec3 v) { return _mm_div_ps(_mm_set_ps1(s), v); }
inline f32vec4 operator/(float32 s, f32vec4 v) { return _mm_div_ps(_mm_set_ps1(s), v); }

// MATRIX OPERATIONS

template<typename T>
void swap(T& x, T& y) {
	T tmp = x;
	x = y;
	y = tmp;
}

inline f32mat2 transpose(f32mat2 m) {
	swap(m[0][1], m[1][0]);
	return m;
}

inline f32mat3 transpose(f32mat3 m) {
	// Swap non diag-elements
	swap(m[0][1], m[1][0]);
	swap(m[0][2], m[2][0]);
	swap(m[1][2], m[2][1]);

	return m;
}

inline f32mat4 transpose(f32mat4 m)
{
	__m128& col0 = m.col[0].m128;
	__m128& col1 = m.col[1].m128;
	__m128& col2 = m.col[2].m128;
	__m128& col3 = m.col[3].m128;

	__m128 tmp3, tmp2, tmp1, tmp0;

	tmp0 = _mm_shuffle_ps(col0, col1, 0x44);      
	tmp2 = _mm_shuffle_ps(col0, col1, 0xEE);          
	tmp1 = _mm_shuffle_ps(col2, col3, 0x44);          
	tmp3 = _mm_shuffle_ps(col2, col3, 0xEE);          

	col0 = _mm_shuffle_ps(tmp0, tmp1, 0x88);
	col1 = _mm_shuffle_ps(tmp0, tmp1, 0xDD);
	col2 = _mm_shuffle_ps(tmp2, tmp3, 0x88);
	col3 = _mm_shuffle_ps(tmp2, tmp3, 0xDD);

	return m;
}

inline f32mat2 operator*(f32mat2 a, f32mat2 b) {
	b = transpose(b);
	return { a[0] * b[0], a[1] * b[1] };
}

inline f32mat3 operator*(f32mat3 a, f32mat3 b) {
	b = transpose(b);
	return { a[0] * b[0], a[1] * b[1], a[2] * b[2] };
}

inline f32mat4 operator*(f32mat4 a, f32mat4 b) {
	b = transpose(b);
	return { a[0] * b[0], a[1] * b[1], a[2] * b[2], a[3] * b[3] };
}

inline f32vec2 operator*(f32mat2 m, f32vec2 v) {
	return { m[0][0] * v[0] + m[1][0] * v[1], m[0][1] * v[0] + m[1][1] * v[1] };
}

// OPERATORS INT
i32vec2 operator+(const i32vec2 &, const i32vec2 &);
i32vec3 operator+(const i32vec3 &, const i32vec3 &);
i32vec4 operator+(const i32vec4 &, const i32vec4 &);

i32vec2 operator-(const i32vec2 &, const i32vec2 &);
i32vec3 operator-(const i32vec3 &, const i32vec3 &);
i32vec4 operator-(const i32vec4 &, const i32vec4 &);

i32vec2 operator*(const i32vec2 &, const i32vec2 &);
i32vec3 operator*(const i32vec3 &, const i32vec3 &);
i32vec4 operator*(const i32vec4 &, const i32vec4 &);

i32vec2 operator*(const i32vec2 &, const int32 &);
i32vec3 operator*(const i32vec3 &, const int32 &);
i32vec4 operator*(const i32vec4 &, const int32 &);

i32vec2 operator*(const int32 &, const i32vec2 &);
i32vec3 operator*(const int32 &, const i32vec3 &);
i32vec4 operator*(const int32 &, const i32vec4 &);

i32vec2 operator/(const i32vec2 &, const i32vec2 &);
i32vec3 operator/(const i32vec3 &, const i32vec3 &);
i32vec4 operator/(const i32vec4 &, const i32vec4 &);

i32vec2 operator/(const int32 &, const i32vec2 &);
i32vec3 operator/(const int32 &, const i32vec3 &);
i32vec4 operator/(const int32 &, const i32vec4 &);

i32vec2 operator/(const i32vec2 &, const int32 &);
i32vec3 operator/(const i32vec3 &, const int32 &);
i32vec4 operator/(const i32vec4 &, const int32 &);

}  // namespace math