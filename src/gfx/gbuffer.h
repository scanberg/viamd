#pragma once

#include "gl.h"

#include <cstdint>

#define GL_COLOR_ATTACHMENT_COLOR        GL_COLOR_ATTACHMENT0
#define GL_COLOR_ATTACHMENT_NORMAL       GL_COLOR_ATTACHMENT1
#define GL_COLOR_ATTACHMENT_VELOCITY     GL_COLOR_ATTACHMENT2
#define GL_COLOR_ATTACHMENT_PICKING      GL_COLOR_ATTACHMENT3
#define GL_COLOR_ATTACHMENT_TRANSPARENCY GL_COLOR_ATTACHMENT4

enum class GBufferFlags : uint32_t {
    None = 0,
    Depth = 1 << 0,
    Color = 1 << 1,
    Normal = 1 << 2,
    Velocity = 1 << 3,
    Picking = 1 << 4,
    Transparency = 1 << 5,
    History = 1 << 6,
};

struct GBuffer {
    struct {
        uint32_t depth = 0;
        uint32_t color = 0;
        uint32_t normal = 0;
        uint32_t velocity = 0;
        uint32_t picking = 0;
        uint32_t transparency = 0;
        uint32_t transparency_hdr = 0;  // RGBA16F, premultiplied linear radiance blended over the scene before tone mapping (isosurfaces). Not attached to fbo.
        uint32_t history = 0;       // TAA history written this frame, at the output size
        uint32_t history_prev = 0;  // TAA history of the previous frame; swap the two after post-processing
    } tex;

    uint32_t fbo = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t out_width = 0;         // what post-processing makes of it: the size of the image shown, and of the history
    uint32_t out_height = 0;
};

// out_width / out_height: the size of the image post-processing makes, when it is upscaled from the G-buffer's;
// 0: the G-buffer's own
void gbuffer_init(GBuffer* gbuf, int width, int height, int out_width = 0, int out_height = 0);
void gbuffer_clear(GBuffer* gbuf);
void gbuffer_free(GBuffer* gbuf);
