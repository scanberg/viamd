#pragma once

#include <core/md_str.h>

struct md_allocator_i;

typedef struct image_t {
    int width;
    int height;
    unsigned int* data;
} image_t;

bool image_init(image_t* img, int width, int height, struct md_allocator_i* alloc);
void image_free(image_t* img, struct md_allocator_i* alloc);

bool image_copy(image_t* img, const image_t* other, struct md_allocator_i* alloc);

bool image_read(image_t* img, str_t filename, struct md_allocator_i* alloc);

// @NOTE: quality is between 0-100
bool image_write_jpg(const image_t* img, str_t filename, int quality);
bool image_write_png(const image_t* img, str_t filename);
bool image_write_bmp(const image_t* img, str_t filename);

void image_gaussian_blur(image_t* img, int kernel_width_in_pixels = 4);
