#include <gfx/gbuffer.h>

#include <gfx/gl_utils.h>

#include <core/md_vec_math.h>

void gbuffer_init(GBuffer* gbuf, int width, int height) {
    ASSERT(gbuf);

    bool attach_textures_deferred = false;
    if (!gbuf->fbo) {
        glGenFramebuffers(1, &gbuf->fbo);
        attach_textures_deferred = true;
    }

    if (!gbuf->tex.depth) glGenTextures(1, &gbuf->tex.depth);
    if (!gbuf->tex.color) glGenTextures(1, &gbuf->tex.color);
    if (!gbuf->tex.normal) glGenTextures(1, &gbuf->tex.normal);
    if (!gbuf->tex.velocity) glGenTextures(1, &gbuf->tex.velocity);
    if (!gbuf->tex.transparency) glGenTextures(1, &gbuf->tex.transparency);
    if (!gbuf->tex.picking) glGenTextures(1, &gbuf->tex.picking);
    if (!gbuf->tex.history) glGenTextures(1, &gbuf->tex.history);

    glBindTexture(GL_TEXTURE_2D, gbuf->tex.depth);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH24_STENCIL8, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->tex.color);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->tex.normal);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16, width, height, 0, GL_RG, GL_UNSIGNED_SHORT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->tex.velocity);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16F, width, height, 0, GL_RG, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->tex.picking);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->tex.transparency);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    // Temporal history is not attached to the G-buffer FBO, but it follows the same size/lifetime.
    glBindTexture(GL_TEXTURE_2D, gbuf->tex.history);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R11F_G11F_B10F, width, height, 0, GL_RGB, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, 0);

    gbuf->width = width;
    gbuf->height = height;

    const GLenum draw_buffers[] = {
        GL_COLOR_ATTACHMENT_COLOR,
        GL_COLOR_ATTACHMENT_NORMAL,
        GL_COLOR_ATTACHMENT_VELOCITY,
        GL_COLOR_ATTACHMENT_TRANSPARENCY,
        GL_COLOR_ATTACHMENT_PICKING,
    };

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf->fbo);
    if (attach_textures_deferred) {
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, gbuf->tex.depth, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_STENCIL_ATTACHMENT, GL_TEXTURE_2D, gbuf->tex.depth, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_COLOR, GL_TEXTURE_2D, gbuf->tex.color, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_NORMAL, GL_TEXTURE_2D, gbuf->tex.normal, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_VELOCITY, GL_TEXTURE_2D, gbuf->tex.velocity, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_PICKING, GL_TEXTURE_2D, gbuf->tex.picking, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_TRANSPARENCY, GL_TEXTURE_2D, gbuf->tex.transparency, 0);
    }

    ASSERT(glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
    glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
    glClearColor(0, 0, 0, 0);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

    // Clear history deterministically so first TAA frame has defined previous data.
    GLuint clear_fbo = 0;
    glGenFramebuffers(1, &clear_fbo);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, clear_fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gbuf->tex.history, 0);
    ASSERT(glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);
    glDeleteFramebuffers(1, &clear_fbo);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
}

void gbuffer_clear(GBuffer* gbuffer) {
    ASSERT(gbuffer);

    const GLenum draw_buffers[] = {
        GL_COLOR_ATTACHMENT_COLOR,
        GL_COLOR_ATTACHMENT_NORMAL,
        GL_COLOR_ATTACHMENT_VELOCITY,
        GL_COLOR_ATTACHMENT_PICKING,
        GL_COLOR_ATTACHMENT_TRANSPARENCY,
    };

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuffer->fbo);
    glViewport(0, 0, gbuffer->width, gbuffer->height);

    glDepthMask(1);
    glColorMask(1, 1, 1, 1);
    glStencilMask(0xFF);
    const vec4_t zero = {0, 0, 0, 0};
    const vec4_t picking = {1, 1, 1, 1};

    glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
    glClearBufferfv(GL_COLOR, 0, zero.elem);
    glClearBufferfv(GL_COLOR, 1, zero.elem);
    glClearBufferfv(GL_COLOR, 2, zero.elem);
    glClearBufferfv(GL_COLOR, 3, picking.elem);
    glClearBufferfv(GL_COLOR, 4, zero.elem);
    glClearBufferfi(GL_DEPTH_STENCIL, 0, 1.0f, 0x01);
}

void gbuffer_free(GBuffer* gbuf) {
    ASSERT(gbuf);
    if (gbuf->fbo) glDeleteFramebuffers(1, &gbuf->fbo);
    if (gbuf->tex.depth) glDeleteTextures(1, &gbuf->tex.depth);
    if (gbuf->tex.color) glDeleteTextures(1, &gbuf->tex.color);
    if (gbuf->tex.normal) glDeleteTextures(1, &gbuf->tex.normal);
    if (gbuf->tex.velocity) glDeleteTextures(1, &gbuf->tex.velocity);
    if (gbuf->tex.transparency) glDeleteTextures(1, &gbuf->tex.transparency);
    if (gbuf->tex.picking) glDeleteTextures(1, &gbuf->tex.picking);
    if (gbuf->tex.history) glDeleteTextures(1, &gbuf->tex.history);
}
