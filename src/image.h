#pragma once

#include <core/types.h>
#include <core/string_types.h>

struct Image {
    i32 width = 0;
    i32 height = 0;
    u32* data = nullptr;
    operator bool() const { return data != nullptr; }
};

bool init_image(Image* img, i32 width, i32 height);
bool init_image(Image* img, const Image& other);
void free_image(Image* img);

bool read_image(Image* img, CStringView filename);
// @NOTE: quality is between 0-100
bool write_image_jpg(const Image& img, CStringView filename, int quality);
bool write_image_png(const Image& img, CStringView filename);
bool write_image_bmp(const Image& img, CStringView filename);

void gaussian_blur(Image* img, i32 kernel_width = 4);
