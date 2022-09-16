#pragma once

#include <stdint.h>
#include <core/md_str.h>

struct md_allocator_i;

typedef struct image_t {
    int32_t width;
    int32_t height;
    uint32_t* data;
} image_t;

bool init_image(image_t* img, int32_t width, int32_t height, struct md_allocator_i* alloc);
bool init_image(image_t* img, image_t other, struct md_allocator_i* alloc);
void free_image(image_t* img, struct md_allocator_i* alloc);

bool read_image(image_t* img, str_t filename, struct md_allocator_i* alloc);

// @NOTE: quality is between 0-100
bool write_image_jpg(const image_t img, str_t filename, int quality);
bool write_image_png(const image_t img, str_t filename);
bool write_image_bmp(const image_t img, str_t filename);

void gaussian_blur(image_t* img, int32_t kernel_width = 4);
