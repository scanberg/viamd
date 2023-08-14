#ifndef  _CRT_SECURE_NO_WARNINGS
#define  _CRT_SECURE_NO_WARNINGS
#endif // ! _CRT_SECURE_NO_WARNINGS

#include <image.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>

#include <color_utils.h>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include <stb_image.h>
#include <stb_image_write.h>

#include <math.h>

bool image_init(image_t* img, int32_t width, int32_t height, md_allocator_i* alloc) {
    ASSERT(img);
    ASSERT(width > 0);
    ASSERT(height > 0);

    if (img->data) {
        MD_LOG_DEBUG("Image struct is not empty, possibly leaking memory here");
    }

    uint32_t* data = (uint32_t*)md_alloc(alloc, width * height * sizeof(uint32_t));
    if (!data) {
        MD_LOG_ERROR("Failed to allocate memory");
        return false;
    }

    img->width = width;
    img->height = height;
    img->data = data;
    
    return true;
}

bool image_copy(image_t* img, const image_t* other, md_allocator_i* alloc) {
    ASSERT(img);
    ASSERT(other);

    if (!image_init(img, other->width, other->height, alloc)) {
        return false;
    }

    MEMCPY(img->data, other->data, other->width * other->height * sizeof(uint32_t));
    return true;
}

void image_free(image_t* img, md_allocator_i* alloc) {
    ASSERT(img);
    if (img->data) {
        md_free(alloc, img->data, img->width * img->height * sizeof(uint32_t));
    }
    img->width = 0;
    img->height = 0;
    img->data = nullptr;
}

bool image_read(image_t* img, str_t filename, md_allocator_i* alloc) {
    ASSERT(img);

    if (img->data) {
        MD_LOG_DEBUG("Image struct is not empty, possibly leaking memory here");
    }

    // Ensure zero terminated cstr
    filename = str_copy(filename, default_temp_allocator);

    int x, y, channels;
    uint8_t* tmp_data = stbi_load(filename.ptr, &x, &y, &channels, 4);
    if (!tmp_data) return false;

    void *img_data = md_alloc(alloc, x * y * 4);
    memcpy(img_data, tmp_data, x * y * 4);

    img->width  = x;
    img->height = y;
    img->data = (uint32_t*)img_data;

    return true;
}

static void write_func(void* context, void* data, int size) {
    ASSERT(context);
    FILE* file = (FILE*)context;
    fwrite(data, 1, size, file);
}

static inline FILE* open_file(str_t filename) {
    FILE* file = (FILE*)md_file_open(filename, MD_FILE_WRITE | MD_FILE_BINARY);
    if (!file) {
        MD_LOG_ERROR("Failed to open file");
    }
    return file;
}

bool image_write_jpg(const image_t* img, str_t filename, int quality) {
    FILE* file = open_file(filename);
    bool result = false;
    if (file) {
        result = stbi_write_jpg_to_func(write_func, file, img->width, img->height, 4, img->data, quality) != 0;
        fclose(file);
    }
    return result;
}

bool image_write_png(const image_t* img, str_t filename) {
    FILE* file = open_file(filename);
    bool result = false;
    if (file) {
        result = stbi_write_png_to_func(write_func, file, img->width, img->height, 4, img->data, img->width * sizeof(uint32_t)) != 0;
        fclose(file);
    }
    return result;
}

bool image_write_bmp(const image_t* img, str_t filename) {
    FILE* file = open_file(filename);
    bool result = false;
    if (file) {
        result = stbi_write_bmp_to_func(write_func, file, img->width, img->height, 4, img->data) != 0;
        fclose(file);
    }
    return result;
}

// All this is ported and stolen from here and needs to be verified
// http://blog.ivank.net/fastest-gaussian-blur.html

static void boxes_for_gauss_3(int* box_w, int n, float sigma) {  // standard deviation, number of boxes
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

static void box_blur_h(uint32_t* src, uint32_t* dst, int w, int h, int r) {
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

static void box_blur_v(uint32_t* src, uint32_t* dst, int w, int h, int r) {
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

static void box_blur(uint32_t* data_in, uint32_t* data_out, int32_t w, int32_t h, int32_t r) {
    ASSERT(data_in);
    ASSERT(data_out);

    box_blur_h(data_out, data_in, w, h, r);
    box_blur_v(data_in, data_out, w, h, r);
}

void image_gaussian_blur(image_t* img, int32_t radius) {
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
