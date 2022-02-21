#ifndef  _CRT_SECURE_NO_WARNINGS
#define  _CRT_SECURE_NO_WARNINGS
#endif // ! _CRT_SECURE_NO_WARNINGS

#include "image.h"

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_file.h>

#define STBI_MALLOC(sz)                     md_alloc(default_temp_allocator, sz)
#define STBI_REALLOC_SIZED(p,oldsz,newsz)   md_realloc(default_temp_allocator, p, oldsz, newsz)
#define STBI_FREE(p)                        {}

#define STBIW_MALLOC(sz)                     md_alloc(default_temp_allocator, sz)
#define STBIW_REALLOC_SIZED(p,oldsz,newsz)   md_realloc(default_temp_allocator, p, oldsz, newsz)
#define STBIW_FREE(p)                        {}

#include <stb_image.h>
#include <stb_image_write.h>

#include <string.h>
#include <math.h>

#include "color_utils.h"

bool init_image(image_t* img, int32_t width, int32_t height, md_allocator_i* alloc) {
    ASSERT(img);
    ASSERT(width > 0);
    ASSERT(height > 0);

    if (img->data) {
        md_print(MD_LOG_TYPE_DEBUG, "Image struct is not empty, possibly leaking memory here");
    }

    uint32_t* data = (uint32_t*)md_alloc(alloc, width * height * sizeof(uint32_t));
    if (!data) return false;

    img->width = width;
    img->height = height;
    img->data = data;
    return true;
}

bool init_image(image_t* img, const image_t other, md_allocator_i* alloc) {
    ASSERT(img);

    if (!init_image(img, other.width, other.height, alloc)) {
        return false;
    }

    memcpy(img->data, other.data, other.width * other.height * sizeof(uint32_t));
    return true;
}

void free_image(image_t* img, md_allocator_i* alloc) {
    ASSERT(img);
    if (img->data) {
        md_free(alloc, img->data, img->width * img->height * sizeof(uint32_t));
    }
    img->width = 0;
    img->height = 0;
    img->data = nullptr;
}

bool read_image(image_t* img, const char* filename, md_allocator_i* alloc) {
    ASSERT(img);

    if (img->data) {
        md_print(MD_LOG_TYPE_DEBUG, "POSSIBLY LEAKING MEMORY HERE");
        ASSERT(false);
    }

    int x, y, channels;
    uint8_t* tmp_data = stbi_load(filename, &x, &y, &channels, 4);
    if (!tmp_data) return false;

    void *img_data = md_alloc(alloc, x * y * 4);
    memcpy(img_data, tmp_data, x * y * 4);

    img->width = x;
    img->height = y;
    img->data = (uint32_t*)img_data;

    return true;
}

static void write_func(void* context, void* data, int size) {
    ASSERT(context);
    FILE* file = (FILE*)context;
    fwrite(data, 1, size, file);
}

static FILE* open_file(const char* filename) {
    FILE* file = fopen(filename, "wb");
    if (!file) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to open file");
    }
    return file;
}

bool write_image_jpg(const image_t img, const char* filename, int quality) {
    FILE* file = open_file(filename);
    bool result = false;
    if (file) {
        result = stbi_write_jpg_to_func(write_func, file, img.width, img.height, 4, img.data, quality) != 0;
        fclose(file);
    }
    return result;
}

bool write_image_png(const image_t img, const char* filename) {
    FILE* file = open_file(filename);
    bool result = false;
    if (file) {
        result = stbi_write_png_to_func(write_func, file, img.width, img.height, 4, img.data, img.width * sizeof(uint32_t)) != 0;
        fclose(file);
    }
    return result;
}

bool write_image_bmp(const image_t img, const char* filename) {
    FILE* file = open_file(filename);
    bool result = false;
    if (file) {
        result = stbi_write_bmp_to_func(write_func, file, img.width, img.height, 4, img.data) != 0;
        fclose(file);
    }
    return result;
}

// All this is ported and stolen from here and needs to be verified
// http://blog.ivank.net/fastest-gaussian-blur.html

void boxes_for_gauss_3(int* box_w, int n, float sigma) {  // standard deviation, number of boxes
    ASSERT(box_w);
    float wIdeal = sqrtf((12 * sigma * sigma / n) + 1);  // Ideal averaging filter width
    int wl = (int)wIdeal;
    if (wl % 2 == 0) wl--;
    int wu = wl + 2;

    float mIdeal = (12 * sigma * sigma - n * wl * wl - 4 * n * wl - 3 * n) / (-4 * wl - 4);
    int m = (int)(mIdeal + 0.5f);
    // var sigmaActual = Math.sqrt( (m*wl*wl + (n-m)*wu*wu - n)/12 );

    for (int i = 0; i < n; i++) box_w[i] = (i < m ? wl : wu);
}

void box_blur_h(uint32_t* src, uint32_t* dst, int w, int h, int r) {
    float iarr = 1.f / (r + r + 1);
    for (int i = 0; i < h; i++) {
        int ti = i * w;
        int li = ti;
        int ri = ti + r;
        const vec4_t fv = convert_color(src[ti]);
        const vec4_t lv = convert_color(src[ti + w - 1]);
        vec4_t val = (float)(r + 1) * fv;

        for (int j = 0; j < r; j++) val = val + convert_color(src[ti + j]);
        for (int j = 0; j <= r; j++) {
            val = val + convert_color(src[ri++]) - fv;
            dst[ti++] = convert_color(val * iarr);
        }
        for (int j = r + 1; j < w - r; j++) {
            val = val + convert_color(src[ri++]) - convert_color(src[li++]);
            dst[ti++] = convert_color(val * iarr);
        }
        for (int j = w - r; j < w; j++) {
            val = val + lv - convert_color(src[li++]);
            dst[ti++] = convert_color(val * iarr);
        }
    }
}

void box_blur_v(uint32_t* src, uint32_t* dst, int w, int h, int r) {
    float iarr = 1.f / (r + r + 1);
    for (int i = 0; i < w; i++) {
        int ti = i;
        int li = ti;
        int ri = ti + r * w;
        const vec4_t fv = convert_color(src[ti]);
        const vec4_t lv = convert_color(src[ti + w * (h - 1)]);
        vec4_t val = (float)(r + 1) * fv;
        for (int j = 0; j < r; j++) val = val + convert_color(src[ti + j * w]);
        for (int j = 0; j <= r; j++) {
            val = val + convert_color(src[ri]) - fv;
            dst[ti] = convert_color(val * iarr);
            ri += w;
            ti += w;
        }
        for (int j = r + 1; j < h - r; j++) {
            val = val + convert_color(src[ri]) - convert_color(src[li]);
            dst[ti] = convert_color(val * iarr);
            li += w;
            ri += w;
            ti += w;
        }
        for (int j = h - r; j < h; j++) {
            val = val + lv - convert_color(src[li]);
            dst[ti] = convert_color(val * iarr);
            li += w;
            ti += w;
        }
    }
}

void box_blur(uint32_t* data_in, uint32_t* data_out, int32_t w, int32_t h, int32_t r) {
    ASSERT(data_in);
    ASSERT(data_out);

    box_blur_h(data_out, data_in, w, h, r);
    box_blur_v(data_in, data_out, w, h, r);
}

void gaussian_blur(image_t* img, int32_t radius) {
    ASSERT(img);
    const int w = img->width;
    const int h = img->height;

    int box_w[3];
    boxes_for_gauss_3(box_w, 3, (float)radius);
    int64_t tmp_size = img->width * img->height * sizeof(uint32_t);
    uint32_t* tmp_data = (uint32_t*)md_alloc(default_temp_allocator, tmp_size);
    ASSERT(tmp_data);
    memcpy(tmp_data, img->data, w * h * sizeof(uint32_t));

    box_blur(tmp_data, img->data, w, h, box_w[0] / 2);
    box_blur(img->data, tmp_data, w, h, box_w[1] / 2);
    box_blur(tmp_data, img->data, w, h, box_w[2] / 2);

    md_free(default_temp_allocator, tmp_data, tmp_size);
}
