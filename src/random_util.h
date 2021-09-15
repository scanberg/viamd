inline float halton(int index, int base) {
    float f = 1;
    float r = 0;
    const float ifb = 1.f / base;
    while (index > 0) {
        f = f * ifb;
        r = r + f * fmodf((float)index, (float)base);
        index = (int)(index * ifb);
    }
    return r;
}

inline void generate_halton_sequence(float* dst, int count, int base) {
    for (int i = 0; i < count; i++) {
        dst[i] = halton(i + 1, base);
    }
}

inline void generate_halton_sequence(vec2_t* dst, int count, int base_x, int base_y) {
    for (int i = 0; i < count; i++) {
        dst[i].x = halton(i + 1, base_x);
        dst[i].y = halton(i + 1, base_y);
    }
}