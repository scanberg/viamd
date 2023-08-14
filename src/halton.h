#pragma once

// https://en.wikipedia.org/wiki/Halton_sequence

static inline float halton(int index, int base) {
    float f = 1;
    float r = 0;
    const float ifb = 1.f / base;
    while (index > 0) {
        f = f * ifb;
        r = r + f * (float)(index % base);
        index = (int)(index * ifb);
    }
    return r;
}
