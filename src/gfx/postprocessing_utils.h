#pragma once

#include "gl.h"
#include "view_param.h"

#include <core/md_vec_math.h>

#define INVALID_PICKING_IDX (~0U)

#define GL_COLOR_ATTACHMENT_COLOR        GL_COLOR_ATTACHMENT0
#define GL_COLOR_ATTACHMENT_NORMAL       GL_COLOR_ATTACHMENT1
#define GL_COLOR_ATTACHMENT_VELOCITY     GL_COLOR_ATTACHMENT2
#define GL_COLOR_ATTACHMENT_PICKING      GL_COLOR_ATTACHMENT3
#define GL_COLOR_ATTACHMENT_TRANSPARENCY GL_COLOR_ATTACHMENT4

// Poor fit perhaps
struct GBuffer {
    struct {
        uint32_t depth = 0;
        uint32_t color = 0;
        uint32_t normal = 0;
        uint32_t velocity = 0;
        uint32_t transparency = 0;
        uint32_t picking = 0;
        uint32_t fbo = 0;
    } deferred;

    struct {
        // @NOTE: Many of each, we submit the read and use it some frame(s) later
        // This means that we read with N-1 frames latency
        uint32_t color[2] = {};
        uint32_t depth[2] = {};
        uint32_t frame = 0;
    } pbo_picking;

    uint32_t width = 0;
    uint32_t height = 0;
};

struct PickingData {
    uint32_t idx = INVALID_PICKING_IDX;
    float depth = 1.0f;
    vec3_t world_coord = {0, 0, 0};
    vec2_t screen_coord = {0, 0};
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
        vec4_t color = {20.f, 20.f, 20.f, 1.0f};
    } background;

    struct {
        bool enabled = true;
        float clip_point = 1.0f;
    } bloom;

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
        float feedback_min = 0.88f;
        float feedback_max = 0.97f;
        struct {
            bool enabled = true;
            float motion_scale = 0.5f;
        } motion_blur;
    } temporal_reprojection;

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

void clear_gbuffer(GBuffer* gbuf);
void init_gbuffer(GBuffer* gbuf, int width, int height);
void destroy_gbuffer(GBuffer* gbuf);
PickingData read_picking_data(GBuffer* fbo, int x, int y);
