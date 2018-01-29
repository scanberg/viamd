#pragma once

#include "types.h"

namespace math {

float32 dot(const F32Vec2&, const F32Vec2&);
float32 dot(const F32Vec3&, const F32Vec3&);
float32 dot(const F32Vec4&, const F32Vec4&);

int32 dot(const I32Vec2&, const I32Vec2&);
int32 dot(const I32Vec3&, const I32Vec3&);
int32 dot(const I32Vec4&, const I32Vec4&);

F32Vec3 cross(const F32Vec3&, const F32Vec3&);

// FLOAT
F32Vec2 operator + (const F32Vec2&, const F32Vec2&);
F32Vec3 operator + (const F32Vec3&, const F32Vec3&);
F32Vec4 operator + (const F32Vec4&, const F32Vec4&);

F32Vec2 operator - (const F32Vec2&, const F32Vec2&);
F32Vec3 operator - (const F32Vec3&, const F32Vec3&);
F32Vec4 operator - (const F32Vec4&, const F32Vec4&);

F32Vec2 operator * (const F32Vec2&, const F32Vec2&);
F32Vec3 operator * (const F32Vec3&, const F32Vec3&);
F32Vec4 operator * (const F32Vec4&, const F32Vec4&);

F32Vec2 operator * (const F32Vec2&, const int32&);
F32Vec3 operator * (const F32Vec3&, const int32&);
F32Vec4 operator * (const F32Vec4&, const int32&);

F32Vec2 operator * (const int32&, const F32Vec2&);
F32Vec3 operator * (const int32&, const F32Vec3&);
F32Vec4 operator * (const int32&, const F32Vec4&);

F32Vec2 operator / (const F32Vec2&, const F32Vec2&);
F32Vec3 operator / (const F32Vec3&, const F32Vec3&);
F32Vec4 operator / (const F32Vec4&, const F32Vec4&);

F32Vec2 operator / (const int32&, const F32Vec2&);
F32Vec3 operator / (const int32&, const F32Vec3&);
F32Vec4 operator / (const int32&, const F32Vec4&);

F32Vec2 operator / (const F32Vec2&, const int32&);
F32Vec3 operator / (const F32Vec3&, const int32&);
F32Vec4 operator / (const F32Vec4&, const int32&);

// INTEGER
I32Vec2 operator + (const I32Vec2&, const I32Vec2&);
I32Vec3 operator + (const I32Vec3&, const I32Vec3&);
I32Vec4 operator + (const I32Vec4&, const I32Vec4&);

I32Vec2 operator - (const I32Vec2&, const I32Vec2&);
I32Vec3 operator - (const I32Vec3&, const I32Vec3&);
I32Vec4 operator - (const I32Vec4&, const I32Vec4&);

I32Vec2 operator * (const I32Vec2&, const I32Vec2&);
I32Vec3 operator * (const I32Vec3&, const I32Vec3&);
I32Vec4 operator * (const I32Vec4&, const I32Vec4&);

I32Vec2 operator * (const I32Vec2&, const int32&);
I32Vec3 operator * (const I32Vec3&, const int32&);
I32Vec4 operator * (const I32Vec4&, const int32&);

I32Vec2 operator * (const int32&, const I32Vec2&);
I32Vec3 operator * (const int32&, const I32Vec3&);
I32Vec4 operator * (const int32&, const I32Vec4&);

I32Vec2 operator / (const I32Vec2&, const I32Vec2&);
I32Vec3 operator / (const I32Vec3&, const I32Vec3&);
I32Vec4 operator / (const I32Vec4&, const I32Vec4&);

I32Vec2 operator / (const int32&, const I32Vec2&);
I32Vec3 operator / (const int32&, const I32Vec3&);
I32Vec4 operator / (const int32&, const I32Vec4&);

I32Vec2 operator / (const I32Vec2&, const int32&);
I32Vec3 operator / (const I32Vec3&, const int32&);
I32Vec4 operator / (const I32Vec4&, const int32&);

}