#pragma once

#include <core/types.h>
#include <core/string_types.h>

struct Image {
    int32 width = 0;
    int32 height = 0;
    uint32* data = nullptr;
};

bool init_image(Image* img, int32 width, int32 height, uint32 color = 0xffffffff);
void free_image(Image* img);

bool read_image(Image* img, CString filename);

void box_blur(Image* src, Image* dst, int32 kernel_width = 4);
void gaussian_blur(Image* src, Image* dst, int32 kernel_width = 4);
