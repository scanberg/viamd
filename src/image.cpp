#include <image.h>

#include <core/md_log.h>
#include <core/md_os.h>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_PNG
#define STBI_ONLY_JPEG
#include <stb_image.h>

uint8_t* image_decode_rgba(const void* data, size_t size, int* width, int* height) {
    int channels = 0;
    return stbi_load_from_memory((const stbi_uc*)data, (int)size, width, height, &channels, 4);
}

void image_free(void* pixels) {
    stbi_image_free(pixels);
}

static void write_func(void* context, void* data, int size) {
    ASSERT(context);
    md_file_t* file = (md_file_t*)context;
    md_file_write(*file, data, size);
}

static inline md_file_t open_file(str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        MD_LOG_ERROR("Failed to open file");
    }
    return file;
}

bool image_write_jpg(str_t filename, const void* rgba, int width, int height, int quality) {
    md_file_t file = open_file(filename);
    bool result = false;
    if (md_file_valid(file)) {
        result = stbi_write_jpg_to_func(write_func, &file, width, height, 4, rgba, quality) != 0;
        md_file_close(&file);
    }
    return result;
}

bool image_write_png(str_t filename, const void* rgba, int width, int height) {
    md_file_t file = open_file(filename);
    bool result = false;
    if (md_file_valid(file)) {
        result = stbi_write_png_to_func(write_func, &file, width, height, 4, rgba, width * sizeof(uint32_t)) != 0;
        md_file_close(&file);
    }
    return result;
}

bool image_write_bmp(str_t filename, const void* rgba, int width, int height) {
    md_file_t file = open_file(filename);
    bool result = false;
    if (md_file_valid(file)) {
        result = stbi_write_bmp_to_func(write_func, &file, width, height, 4, rgba) != 0;
        md_file_close(&file);
    }
    return result;
}

bool image_alpha_bounds(const uint8_t* rgba, int width, int height, uint8_t threshold, int* x0, int* y0, int* x1, int* y1) {
    int minx = width, miny = height, maxx = -1, maxy = -1;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (rgba[((size_t)y * (size_t)width + (size_t)x) * 4 + 3] > threshold) {
                if (x < minx) minx = x;
                if (x > maxx) maxx = x;
                if (y < miny) miny = y;
                if (y > maxy) maxy = y;
            }
        }
    }
    if (maxx < 0) return false;
    *x0 = minx;
    *y0 = miny;
    *x1 = maxx + 1;
    *y1 = maxy + 1;
    return true;
}

void image_bleed_transparent(uint8_t* rgba, int width, int height) {
    double sum[3] = {0.0, 0.0, 0.0}, weight = 0.0;
    const size_t n = (size_t)width * (size_t)height;
    for (size_t i = 0; i < n; ++i) {
        const double a = (double)rgba[i * 4 + 3];
        for (int c = 0; c < 3; ++c) sum[c] += a * (double)rgba[i * 4 + c];
        weight += a;
    }
    if (weight <= 0.0) return;
    uint8_t avg[3];
    for (int c = 0; c < 3; ++c) avg[c] = (uint8_t)(sum[c] / weight + 0.5);
    for (size_t i = 0; i < n; ++i) {
        if (rgba[i * 4 + 3] == 0) {
            for (int c = 0; c < 3; ++c) rgba[i * 4 + c] = avg[c];
        }
    }
}

void image_replace_light_grey(uint8_t* rgba, int width, int height, uint8_t min_value, uint8_t tolerance, uint8_t r, uint8_t g, uint8_t b) {
    const size_t n = (size_t)width * (size_t)height;
    for (size_t i = 0; i < n; ++i) {
        uint8_t* p = &rgba[i * 4];
        if (p[3] == 0) continue;
        const uint8_t hi = p[0] > p[1] ? (p[0] > p[2] ? p[0] : p[2]) : (p[1] > p[2] ? p[1] : p[2]);
        const uint8_t lo = p[0] < p[1] ? (p[0] < p[2] ? p[0] : p[2]) : (p[1] < p[2] ? p[1] : p[2]);
        if (hi - lo <= tolerance && lo >= min_value) {
            p[0] = r;
            p[1] = g;
            p[2] = b;
        }
    }
}
