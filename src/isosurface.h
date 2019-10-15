#pragma once

#include <core/vector_types.h>

struct IsoSurfaces {
    static constexpr int MaxCount = 8;
    float values[MaxCount] = {};
    vec4 colors[MaxCount] = {};
    int count = 0;
};

bool insert(IsoSurfaces& surface, float value, const vec4& color);
void clear(IsoSurfaces& surface);
void sort(IsoSurfaces& surface);