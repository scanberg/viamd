#include "math_utils.h"

constexpr int mask_all = 0xFF;

namespace math {

float32 dot(const F32Vec2& a, const F32Vec2& b) {
    return a.x * b.x + a.y * b.y;
}
float32 dot(const F32Vec3& a, const F32Vec3& b) {
    F32Vec3 res;
    res.m128 = _mm_dp_ps(a.m128, b.m128, mask_all);
    return res[0];
}
/*
float32 dot(const F32Vec4& a, const F32Vec4& b) {
    F32Vec4 res = _mm_dp_ps(a, b, mask_all);
    return res[0];
}

int32 dot(const I32Vec2& a, const I32Vec2& b) {

}
int32 dot(const I32Vec3& a, const I32Vec3& b) {

}
int32 dot(const I32Vec4& a, const I32Vec4& b) {

}

F32Vec3 cross(const F32Vec3& a, const F32Vec3& b) {

}

// FLOAT
F32Vec2 operator + (const F32Vec2& a, const F32Vec2& b) {

}
F32Vec3 operator + (const F32Vec3& a, const F32Vec3& b) {

}
F32Vec4 operator + (const F32Vec4& a, const F32Vec4& b) {

}

F32Vec2 operator - (const F32Vec2& a, const F32Vec2& b) {

}
F32Vec3 operator - (const F32Vec3& a, const F32Vec3& b) {

}
F32Vec4 operator - (const F32Vec4& a, const F32Vec4& b) {

}

F32Vec2 operator * (const F32Vec2& a, const F32Vec2& b) {

}
F32Vec3 operator * (const F32Vec3& a, const F32Vec3& b) {

}
F32Vec4 operator * (const F32Vec4& a, const F32Vec4& b) {

}

F32Vec2 operator * (const F32Vec2& a, const int32& b) {

}
F32Vec3 operator * (const F32Vec3& a, const int32& b) {

}
F32Vec4 operator * (const F32Vec4& a, const int32& b) {

}

F32Vec2 operator * (const int32& a, const F32Vec2& b) {

}
F32Vec3 operator * (const int32& a, const F32Vec3& b) {

}
F32Vec4 operator * (const int32& a, const F32Vec4& b) {

}

F32Vec2 operator / (const F32Vec2& a, const F32Vec2& b) {

}
F32Vec3 operator / (const F32Vec3& a, const F32Vec3& b) {

}
F32Vec4 operator / (const F32Vec4& a, const F32Vec4& b) {

}

F32Vec2 operator / (const int32& a, const F32Vec2& b) {

}
F32Vec3 operator / (const int32& a, const F32Vec3& b) {

}
F32Vec4 operator / (const int32& a, const F32Vec4& b) {

}

F32Vec2 operator / (const F32Vec2& a, const int32& b) {

}
F32Vec3 operator / (const F32Vec3& a, const int32& b) {

}
F32Vec4 operator / (const F32Vec4& a, const int32& b) {

}

// INTEGER
I32Vec2 operator + (const I32Vec2& a, const I32Vec2& b) {

}
I32Vec3 operator + (const I32Vec3& a, const I32Vec3& b) {

}
I32Vec4 operator + (const I32Vec4& a, const I32Vec4& b) {

}

I32Vec2 operator - (const I32Vec2& a, const I32Vec2& b) {

}
I32Vec3 operator - (const I32Vec3& a, const I32Vec3& b) {

}
I32Vec4 operator - (const I32Vec4& a, const I32Vec4& b) {

}

I32Vec2 operator * (const I32Vec2& a, const I32Vec2& b) {

}
I32Vec3 operator * (const I32Vec3& a, const I32Vec3& b) {

}
I32Vec4 operator * (const I32Vec4& a, const I32Vec4& b) {

}

I32Vec2 operator * (const I32Vec2& a, const int32& b) {

}
I32Vec3 operator * (const I32Vec3& a, const int32& b) {

}
I32Vec4 operator * (const I32Vec4& a, const int32& b) {

}

I32Vec2 operator * (const int32& a, const I32Vec2& b) {

}
I32Vec3 operator * (const int32& a, const I32Vec3& b) {

}
I32Vec4 operator * (const int32& a, const I32Vec4& b) {

}

I32Vec2 operator / (const I32Vec2& a, const I32Vec2& b) {

}
I32Vec3 operator / (const I32Vec3& a, const I32Vec3& b) {

}
I32Vec4 operator / (const I32Vec4& a, const I32Vec4& b) {

}

I32Vec2 operator / (const int32& a, const I32Vec2& b) {

}
I32Vec3 operator / (const int32& a, const I32Vec3& b) {

}
I32Vec4 operator / (const int32& a, const I32Vec4& b) {

}

I32Vec2 operator / (const I32Vec2& a, const int32& b) {

}
I32Vec3 operator / (const I32Vec3& a, const int32& b) {

}
I32Vec4 operator / (const I32Vec4& a, const int32& b) {

}
*/
}