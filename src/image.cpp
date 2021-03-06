#include "image.h"
#include <core/common.h>
#include <core/log.h>
#include <core/file.h>
#include <core/vector_types.h>
#include <core/math_utils.h>
#include <stb_image.h>
#include <stb_image_write.h>

bool init_image(Image* img, i32 width, i32 height) {
    ASSERT(img);
    ASSERT(width > 0);
    ASSERT(height > 0);

    u32* data = (u32*)MALLOC(width * height * sizeof(u32));
    if (!data) return false;

    free_image(img);
    img->width = width;
    img->height = height;
    img->data = data;
    return true;
}

bool init_image(Image* img, const Image& other) {
    ASSERT(img);
    ASSERT(other);

    if (!init_image(img, other.width, other.height)) {
        return false;
    }
    memcpy(img->data, other.data, other.width * other.height * sizeof(u32));
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

bool read_image(Image* img, CStringView filename) {
    ASSERT(img);

    StringBuffer<512> zstr = filename;
    int x, y, channels;
    u8* data = stbi_load(zstr.cstr(), &x, &y, &channels, 4);
    if (!data) return false;

    free_image(img);
    img->width = x;
    img->height = y;
    img->data = (u32*)data;

    return true;
}

static void write_func(void* context, void* data, int size) {
    ASSERT(context);
    FILE* file = (FILE*)context;
    fwrite(data, 1, size, file);
}

static FILE* open_file(CStringView filename) {
    FILE* file = fopen(filename, "wb");
    if (!file) {
        LOG_ERROR("Failed to open file")
    }
    return file;
}

bool write_image_jpg(const Image& img, CStringView filename, int quality) {
    FILE* file = open_file(filename);
    defer { fclose(file); };
    return file && stbi_write_jpg_to_func(write_func, file, img.width, img.height, 4, img.data, quality) != 0;
}

bool write_image_png(const Image& img, CStringView filename) {
    FILE* file = open_file(filename);
    defer { fclose(file); };
    return file && stbi_write_png_to_func(write_func, file, img.width, img.height, 4, img.data, img.width * sizeof(u32)) != 0;
}

bool write_image_bmp(const Image& img, CStringView filename) {
    FILE* file = open_file(filename);
    defer { fclose(file); };
    return stbi_write_bmp_to_func(write_func, file, img.width, img.height, 4, img.data) != 0;
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

void box_blur_h(u32* src, u32* dst, int w, int h, int r) {
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

void box_blur_v(u32* src, u32* dst, int w, int h, int r) {
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

void box_blur(u32* data_in, u32* data_out, i32 w, i32 h, i32 r) {
    ASSERT(data_in);
    ASSERT(data_out);

    box_blur_h(data_out, data_in, w, h, r);
    box_blur_v(data_in, data_out, w, h, r);
}

void gaussian_blur(Image* img, i32 radius) {
    ASSERT(img);
    const int w = img->width;
    const int h = img->height;

    int box_w[3];
    boxes_for_gauss_3(box_w, 3, (float)radius);
    u32* tmp_data = (u32*)TMP_MALLOC(img->width * img->height * sizeof(u32));
    ASSERT(tmp_data);
    memcpy(tmp_data, img->data, w * h * sizeof(u32));

    box_blur(tmp_data, img->data, w, h, box_w[0] / 2);
    box_blur(img->data, tmp_data, w, h, box_w[1] / 2);
    box_blur(tmp_data, img->data, w, h, box_w[2] / 2);
}
