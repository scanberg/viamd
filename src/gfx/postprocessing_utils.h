#pragma once

#include "gl.h"
#include "view_param.h"

#include <core/md_vec_math.h>

enum GbufferTexture : uint32_t {
    GBUFFER_TEX_COLOR = 0,
    GBUFFER_TEX_NORMAL = 1,
    GBUFFER_TEX_VELOCITY = 2,
    GBUFFER_TEX_PICKING = 3,
    GBUFFER_TEX_TRANSPARENCY = 4,
};

enum GBufferMask : uint32_t {
    GBUFFER_MASK_NONE         = 0,

    // COLOR OUTPUT BUFFERS
    GBUFFER_MASK_COLOR        = 1 << GBUFFER_TEX_COLOR,
    GBUFFER_MASK_NORMAL       = 1 << GBUFFER_TEX_NORMAL,
    GBUFFER_MASK_VELOCITY     = 1 << GBUFFER_TEX_VELOCITY,
    GBUFFER_MASK_PICKING      = 1 << GBUFFER_TEX_PICKING,
    GBUFFER_MASK_TRANSPARENCY = 1 << GBUFFER_TEX_TRANSPARENCY,

    GBUFFER_MASK_DEPTH        = 1 << 16,
    GBUFFER_MASK_STENCIL      = 1 << 17,

    GBUFFER_MASK_ALL          = 0xFFFFFFFFu
};

ENUM_FLAGS(GBufferMask)

#define GL_COLOR_ATTACHMENT_COLOR        (GL_COLOR_ATTACHMENT0 + GBUFFER_TEX_COLOR)
#define GL_COLOR_ATTACHMENT_NORMAL       (GL_COLOR_ATTACHMENT0 + GBUFFER_TEX_NORMAL)
#define GL_COLOR_ATTACHMENT_VELOCITY     (GL_COLOR_ATTACHMENT0 + GBUFFER_TEX_VELOCITY)
#define GL_COLOR_ATTACHMENT_PICKING      (GL_COLOR_ATTACHMENT0 + GBUFFER_TEX_PICKING)
#define GL_COLOR_ATTACHMENT_TRANSPARENCY (GL_COLOR_ATTACHMENT0 + GBUFFER_TEX_TRANSPARENCY)

// Poor fit perhaps
struct GBuffer {
    struct {
        uint32_t depth = 0;
        uint32_t color = 0;
        uint32_t normal = 0;
        uint32_t velocity = 0;
        uint32_t picking = 0;
        uint32_t transparency = 0;
    } tex;

    uint32_t fbo = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    GBufferMask mask = GBUFFER_MASK_NONE;
};

namespace postprocessing {

void initialize(int width, int height);
void shutdown();

typedef int Tonemapping;
enum Tonemapping_ {
    Tonemapping_Passthrough,
    Tonemapping_ExposureGamma,
    Tonemapping_Filmic,
    Tonemapping_ACES,
};

struct Descriptor {
    struct {
        vec3_t color = {20.f, 20.f, 20.f};
    } background;

#if 0
    struct {
        bool enabled = true;
        float clip_point = 1.0f;
    } bloom;
#endif

    struct {
        bool enabled = true;
        Tonemapping mode = Tonemapping_ACES;
        float exposure = 1.0f;
        float gamma = 2.4f;
    } tonemapping;

    struct {
        bool enabled = true;
        float radius = 6.0f;
        float intensity = 3.0f;
        float bias = 0.1f;
    } ambient_occlusion;

    struct {
        bool enabled = true;
        float focus_depth = 0.5f;
        float focus_scale = 10.f;
    } depth_of_field;

    struct {
        bool enabled = true;
    } fxaa;

    struct {
        bool enabled = true;
        float feedback_min = 0.80f;
        float feedback_max = 0.95f;
        struct {
            bool enabled = true;
            float motion_scale = 0.5f;
        } motion_blur;
    } temporal_aa;

    struct {
        bool enabled = true;
        float weight = 1.0f;
    } sharpen;

    struct {
        GLuint depth = 0;
        GLuint color = 0;
        GLuint normal = 0;
        GLuint velocity = 0;
        GLuint transparency = 0;
    } input_textures;
};

void apply_tonemapping(GLuint color_tex, Tonemapping tonemapping, float exposure = 1.0f, float gamma = 2.4f);

void shade_and_postprocess(const Descriptor& desc, const ViewParam& view_param);

void blit_static_velocity(GLuint tex_depth, const ViewParam& view_param);

void scale_hsv(GLuint color_tex, vec3_t hsv_scale);

void blit_texture(GLuint tex);
void blit_color(vec4_t color);

void blur_texture_gaussian(GLuint tex, int num_passes = 1);
void blur_texture_box(GLuint tex, int num_passes = 1);

}  // namespace postprocessing

void gbuffer_clear(const GBuffer& gbuf, GBufferMask mask = GBUFFER_MASK_ALL);
void gbuffer_bind(const GBuffer& gbuf, GBufferMask mask = GBUFFER_MASK_ALL);

void gbuffer_init(GBuffer* gbuf, int width, int height, GBufferMask mask = GBUFFER_MASK_ALL);
void gbuffer_free(GBuffer* gbuf);
