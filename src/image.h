#pragma once

#include <core/md_str.h>

// @NOTE: JPG quality is between 0-100
// rgba: 32-bit per pixel, 8-bit per channel, rgba
bool image_write_jpg(str_t filename, const void* rgba, int width, int height, int quality);
bool image_write_png(str_t filename, const void* rgba, int width, int height);
bool image_write_bmp(str_t filename, const void* rgba, int width, int height);
