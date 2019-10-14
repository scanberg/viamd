#pragma once

#include <core/vector_types.h>

struct IsoSurface {
    static constexpr int MaxCount = 8;
    float values[MaxCount] = {};
    vec4 colors[MaxCount] = {};
    int count = 0;
};

bool insert(IsoSurface& surface, float value, const vec4& color);
void clear(IsoSurface& surface);
void sort(IsoSurface& surface);