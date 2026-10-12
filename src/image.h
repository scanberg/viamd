#pragma once

#include <core/md_str.h>

#include <stdint.h>
#include <stddef.h>

// @NOTE: JPG quality is between 0-100
// rgba: 32-bit per pixel, 8-bit per channel, rgba
bool image_write_jpg(str_t filename, const void* rgba, int width, int height, int quality);
bool image_write_png(str_t filename, const void* rgba, int width, int height);
bool image_write_bmp(str_t filename, const void* rgba, int width, int height);

// Decodes an image held in memory (png, jpg ...) to rgba, 8 bits per channel. The pixels are freed with image_free.
// Returns null if it is not an image.
uint8_t* image_decode_rgba(const void* data, size_t size, int* width, int* height);
void image_free(void* pixels);

// The smallest rectangle [x0, x1) x [y0, y1) that holds every pixel with an alpha above 'threshold'. False when there is none.
bool image_alpha_bounds(const uint8_t* rgba, int width, int height, uint8_t threshold, int* x0, int* y0, int* x1, int* y1);

// Gives the visible pixels that are grey (their channels differ by at most 'tolerance') and at least 'min_value' bright the
// colour (r, g, b), keeping their alpha. Coloured and dark pixels are left as they are.
void image_replace_light_grey(uint8_t* rgba, int width, int height, uint8_t min_value, uint8_t tolerance, uint8_t r, uint8_t g, uint8_t b);

// Gives the pixels that are fully transparent the colour that the rest has on average, so that smoothing between them and
// the visible pixels (a texture drawn smaller, or with mipmaps) does not leave a dark fringe.
void image_bleed_transparent(uint8_t* rgba, int width, int height);
