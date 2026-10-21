#pragma once

#include <core/md_str.h>

enum ImageFormat {
    ImageFormat_Unknown = 0,
    ImageFormat_JPG,
    ImageFormat_PNG,
    ImageFormat_BMP,
};

// The image format a file extension (without the dot) names, case-insensitive.
// "jpg" and "jpeg" are the same format: the macOS save panel proposes ".jpeg" for "jpg"
ImageFormat image_format_from_ext(str_t ext);

// @NOTE: JPG quality is between 0-100
// rgba: 32-bit per pixel, 8-bit per channel, rgba
bool image_write_jpg(str_t filename, const void* rgba, int width, int height, int quality);
bool image_write_png(str_t filename, const void* rgba, int width, int height);
bool image_write_bmp(str_t filename, const void* rgba, int width, int height);
