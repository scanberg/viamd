#include "image.h"
#include <core/common.h>
#include <core/vector_types.h>
#include <core/math_utils.h>
#include <stb_image.h>

bool init_image(Image* img, int32 width, int32 height, uint32 color) {
    ASSERT(img);
    ASSERT(width > 0);
    ASSERT(height > 0);

    uint32* data = (uint32*)MALLOC(width * height * sizeof(uint32));
    if (!data) return false;

    free_image(img);
    img->width = width;
    img->height = height;
    img->data = data;
    return true;
}

void free_image(Image* img) {
    ASSERT(img);
    if (img->data) {
        FREE(img->data);
    }
    img->width = 0;
    img->height = 0;
    img->data = nullptr;
}

bool read_image(Image* img, CString filename) {
    ASSERT(img);

    StringBuffer<512> zstr = filename;
    int x, y, channels;
    uint8* data = stbi_load(zstr.cstr(), &x, &y, &channels, 4);
    if (!data) return false;

    free_image(img);
    img->width = x;
    img->height = y;
    img->data = (uint32*)data;

    return true;
}

// All this is ported and stolen from here and needs to be verified
// http://blog.ivank.net/fastest-gaussian-blur.html

void boxes_for_gauss_3(int* box_w, int n, float sigma) {  // standard deviation, number of boxes
    ASSERT(box_w);
    float wIdeal = math::sqrt((12 * sigma * sigma / n) + 1);  // Ideal averaging filter width
    int wl = (int)wIdeal;
    if (wl % 2 == 0) wl--;
    int wu = wl + 2;

    float mIdeal = (12 * sigma * sigma - n * wl * wl - 4 * n * wl - 3 * n) / (-4 * wl - 4);
    int m = (int)(mIdeal + 0.5f);
    // var sigmaActual = Math.sqrt( (m*wl*wl + (n-m)*wu*wu - n)/12 );

    for (int i = 0; i < n; i++) box_w[i] = (i < m ? wl : wu);
}

void box_blur_h(uint32* src, uint32* dst, int w, int h, int r) {
    float iarr = 1.f / (r + r + 1);
    for (int i = 0; i < h; i++) {
        int ti = i * w;
        int li = ti;
        int ri = ti + r;
        const vec4 fv = math::convert_color(src[ti]);
        const vec4 lv = math::convert_color(src[ti + w - 1]);
        vec4 val = (float)(r + 1) * fv;

        for (int j = 0; j < r; j++) val += math::convert_color(src[ti + j]);
        for (int j = 0; j <= r; j++) {
            val += math::convert_color(src[ri++]) - fv;
            dst[ti++] = math::convert_color(val * iarr);
        }
        for (int j = r + 1; j < w - r; j++) {
            val += math::convert_color(src[ri++]) - math::convert_color(src[li++]);
            dst[ti++] = math::convert_color(val * iarr);
        }
        for (int j = w - r; j < w; j++) {
            val += lv - math::convert_color(src[li++]);
            dst[ti++] = math::convert_color(val * iarr);
        }
    }
}

void box_blur_v(uint32* src, uint32* dst, int w, int h, int r) {
    float iarr = 1.f / (r + r + 1);
    for (int i = 0; i < w; i++) {
        int ti = i;
        int li = ti;
        int ri = ti + r * w;
        const vec4 fv = math::convert_color(src[ti]);
        const vec4 lv = math::convert_color(src[ti + w * (h - 1)]);
        vec4 val = (float)(r + 1) * fv;
        for (int j = 0; j < r; j++) val += math::convert_color(src[ti + j * w]);
        for (int j = 0; j <= r; j++) {
            val += math::convert_color(src[ri]) - fv;
            dst[ti] = math::convert_color(val * iarr);
            ri += w;
            ti += w;
        }
        for (int j = r + 1; j < h - r; j++) {
            val += math::convert_color(src[ri]) - math::convert_color(src[li]);
            dst[ti] = math::convert_color(val * iarr);
            li += w;
            ri += w;
            ti += w;
        }
        for (int j = h - r; j < h; j++) {
            val += lv - math::convert_color(src[li]);
            dst[ti] = math::convert_color(val * iarr);
            li += w;
            ti += w;
        }
    }
}

void box_blur(Image* src, Image* dst, int32 radius) {
    ASSERT(src);
    ASSERT(dst);
    ASSERT(src->width == dst->width);
    ASSERT(src->height == dst->height);

    const int w = src->width;
    const int h = src->height;
    const int r = radius;

    memcpy(dst->data, src->data, w * h * sizeof(uint32));
    box_blur_h(dst->data, src->data, w, h, r);
    box_blur_v(src->data, dst->data, w, h, r);

    /*
for (int y = 0; y < h; y++) {
    const int32 y_m = (y == 0) ? 0 : y - 1;
    const int32 y_p = (y == h - 1) ? h - 1 : y + 1;
    for (int x = 0; x < w; x++) {
        const int32 x_m = (x == 0) ? 0 : x - 1;
        const int32 x_p = (x == w - 1) ? w - 1 : x + 1;

        const vec4 val = 1.f / 9.f *
                         (math::convert_color(src->data[y_m * w + x_m]) + math::convert_color(src->data[y_m * w + x]) + math::convert_color(src->data[y_m * w + x_p]) +
                          math::convert_color(src->data[y * w + x_m]) + math::convert_color(src->data[y * w + x]) + math::convert_color(src->data[y * w + x_p]) +
                          math::convert_color(src->data[y_p * w + x_m]) + math::convert_color(src->data[y_p * w + x]) + math::convert_color(src->data[y_p * w + x_p]));

        dst->data[y * w + x] = math::convert_color(val);
    }
}
    */
}

void gaussian_blur(Image* src, Image* dst, int32 radius) {
    int box_w[3];
    boxes_for_gauss_3(box_w, 3, radius);

    box_blur(src, dst, box_w[0] / 2);
    box_blur(dst, src, box_w[1] / 2);
    box_blur(src, dst, box_w[2] / 2);
}
