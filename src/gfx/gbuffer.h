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
        uint32_t history = 0;
    } tex;

    uint32_t fbo = 0;
    uint32_t width = 0;
    uint32_t height = 0;
};

void gbuffer_init(GBuffer* gbuf, int width, int height);
void gbuffer_clear(GBuffer* gbuf);
void gbuffer_free(GBuffer* gbuf);
