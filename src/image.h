#pragma once

#include <core/types.h>
#include <core/string_types.h>

struct Image {
    int32 width = 0;
    int32 height = 0;
    uint32* data = nullptr;
    operator bool() const { return data != nullptr; }
};

bool init_image(Image* img, int32 width, int32 height);
bool init_image(Image* img, const Image& other);
void free_image(Image* img);
bool read_image(Image* img, CString filename);

void gaussian_blur(Image* img, int32 kernel_width = 4);
