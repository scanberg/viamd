/*-----------------------------------------------------------------------
  Copyright (c) 2014, NVIDIA. All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
   * Neither the name of its contributors may be used to endorse
     or promote products derived from this software without specific
     prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
  OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
-----------------------------------------------------------------------*/

// Depth linearization and the proj_info setup are based on NVIDIAs gl_ssao sample and are covered by the notice above

#include <gfx/postprocessing_utils.h>

#include <core/md_str.h>
#include <core/md_log.h>

#include <gfx/gl_utils.h>

#include <float.h>
#include <string.h>

#include <shaders.inl>

#include <gfx/profiler.h>

struct GLResetState {
	GLint fbo = 0;
	GLint draw_buffers[8] = { 0 };
	GLint viewport[4] = { 0 };
	GLint scissor_rect[4] = { 0 };
};

static void record_gl_reset_state(GLResetState* state) {
	ASSERT(state);
	glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &state->fbo);
	glGetIntegerv(GL_VIEWPORT, state->viewport);
	glGetIntegerv(GL_SCISSOR_BOX, state->scissor_rect);
	for (size_t i = 0; i < ARRAY_SIZE(state->draw_buffers); ++i) {
		glGetIntegerv(GL_DRAW_BUFFER0 + i, &state->draw_buffers[i]);
	}
}

static void reset_gl_state(const GLResetState& state) {
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, state.fbo);
	if (state.fbo == 0) {
		// For the default framebuffer glDrawBuffers only accepts NONE, FRONT_LEFT,
		// FRONT_RIGHT, BACK_LEFT and BACK_RIGHT -- FRONT, BACK, LEFT, RIGHT and
		// FRONT_AND_BACK are not accepted, even though glDrawBuffer takes them and
		// the GL_DRAW_BUFFER0 query reports them straight back (glDrawBuffer(GL_BACK)
		// reads back as GL_BACK). Feeding that value to glDrawBuffers is INVALID_ENUM
		// on a strict driver, so restore it through glDrawBuffer instead.
		glDrawBuffer((GLenum)state.draw_buffers[0]);
	} else {
		glDrawBuffers((int)ARRAY_SIZE(state.draw_buffers), (const GLenum*)state.draw_buffers);
	}
	glViewport(state.viewport[0], state.viewport[1], state.viewport[2], state.viewport[3]);
	glScissor(state.scissor_rect[0], state.scissor_rect[1], state.scissor_rect[2], state.scissor_rect[3]);
}

namespace postprocessing {

// @TODO: Use half-res render targets for SSAO
// @TODO: Use shared textures for all postprocessing operations
// @TODO: Use some kind of unified pipeline for all post processing operations

// The display transform (shaders/tonemap/display.frag): x = exposure * L, display = curve(x) / curve(W), then the sRGB
// encoding with the gamma setting on top. L is in units of reference white (postprocess_pipeline::lighting). The scene
// is composited over the background there, after the curve, by its coverage: the background is not part of the scene.
namespace display {

enum Curve {
    Curve_None    = 0,  // clip to [0, 1]
    Curve_Neutral = 1,  // Khronos PBR Neutral
    Curve_ACES    = 2,  // Stephen Hill's fit of the ACES RRT + ODT
    Curve_Hable   = 3,  // John Hable's Uncharted 2 curve
};

struct Transform {
    Curve curve = Curve_None;
    float exposure = 1.0f;
    float gamma = 2.2f;
};

static Transform transform(const postprocess_pipeline::Settings& settings) {
    Transform t = {};
    // Basic shading shows the colours as they are, and no tone mapping means no curve: either way only the encoding
    if (settings.shading.mode == postprocess_pipeline::ShadingMode_Basic || !settings.tonemap.enabled) return t;
    switch (settings.tonemap.mode) {
        case postprocess_pipeline::Tonemapper_Neutral: t.curve = Curve_Neutral; break;
        case postprocess_pipeline::Tonemapper_ACES:    t.curve = Curve_ACES;    break;
        case postprocess_pipeline::Tonemapper_Filmic:  t.curve = Curve_Hable;   break;
        default:                                       t.curve = Curve_None;    break;
    }
    t.exposure = MAX(settings.tonemap.exposure, 1.0e-6f);
    t.gamma = CLAMP(settings.tonemap.gamma, 0.1f, 10.0f);
    return t;
}

static inline float srgb_decode(float x) {
    x = CLAMP(x, 0.0f, 1.0f);
    return (x <= 0.04045f) ? (x / 12.92f) : powf((x + 0.055f) / 1.055f, 2.4f);
}

// The display referred linear value of the background, the sRGB colour c, which the encoding (with its gamma
// adjustment) shows as c
static vec3_t background(vec3_t c, const Transform& t) {
    const float k = t.gamma / 2.2f;
    return {srgb_decode(powf(CLAMP(c.x, 0.0f, 1.0f), k)), srgb_decode(powf(CLAMP(c.y, 0.0f, 1.0f), k)), srgb_decode(powf(CLAMP(c.z, 0.0f, 1.0f), k))};
}

}  // namespace display

static inline bool is_orthographic_proj_matrix(const float P[4][4]) { return P[2][3] == 0.0f; }

// A visibility estimated in screen space at half resolution, and its history over frames (see filter_visibility)
struct Gl_Visibility {
    GLuint tex[2] = {};             // half res RG32F (visibility, linear depth), ping-pong for the separable blur
    GLuint tex_history[2] = {};     // half res RG32F (frame count + visibility, linear depth), ping-pong over frames
    GLuint tex_visibility = 0;      // full res R16F: the upsampled visibility, which the compose pass applies
    int history_curr = 0;           // index of the history written last
    bool history_valid = false;
    int history_width = 0;
    int history_height = 0;
    int history_variant = 0;
    int history_mode = 0;           // estimator the history was accumulated with
    mat4_t history_view = {};       // view and projection the history was written with
    mat4_t history_proj = {};
    uint32_t frame = 0;             // drives the per-frame permutation of the sample pattern
};

static struct {
    GLuint vao = 0;
    GLuint v_shader_fs_quad = 0;
    bool programs_ready = false;    // programs are compiled once; initialize() on resize only reallocates targets
    uint32_t tex_width = 0;         // the size of the G-buffer, and of the targets up to the anti-aliasing
    uint32_t tex_height = 0;
    uint32_t out_width = 0;         // the size of the image shown, when it is upscaled to it: of the TAA's output and on
    uint32_t out_height = 0;

    struct {
        GLuint fbo = 0;
        GLuint scratch_fbo = 0;
        GLuint tex_rgba8 = 0;
        GLuint tex_half_dof = 0;        // half res RGBA16F: colour + signed CoC (DOF prepass)
        GLuint tex_color[2] = {0, 0};   // HDR ping-pong (R11F_G11F_B10F): compose, SSAO, DOF
        GLuint tex_ldr[2] = {0, 0};     // LDR ping-pong (RGB10_A2): everything after tone mapping
        GLuint tex_coverage[2] = {0, 0};// R8: coverage of the geometry, sharp and through the depth of field
        GLuint tex_history_prev = 0;    // at the output size
        GLuint out_fbo = 0;             // the targets at the output size, when upscaling: the TAA's and what follows it
        GLuint tex_out[2] = {0, 0};     // LDR ping-pong (RGB10_A2) at the output size
    } rt;

    struct {
        GLuint fbo = 0;
        GLuint tex_tilemax = 0;
        GLuint tex_neighbormax = 0;
        uint32_t tex_width = 0;
        uint32_t tex_height = 0;
    } velocity;

    struct {
        GLuint fbo = 0;
        GLuint texture = 0;
        int levels = 0;
        struct {
            GLuint program_persp = 0;
            GLuint program_ortho = 0;
        } linearize;
        struct {
            GLuint program = 0;
        } downsample;
    } linear_depth;

    // The screen-space visibilities: of the surroundings (the ambient occlusion) and of the key light (the shadows).
    // Each is estimated at half resolution, then filtered over time and space and upsampled by the same passes.
    struct {
        GLuint fbo = 0;
        Gl_Visibility ao;
        Gl_Visibility shadow;
        // Index 0 = orthographic, 1 = perspective
        GLuint program_ao[2] = {};          // ssao/ssao.frag:          point samples ("Performance")
        GLuint program_ao_quality[2] = {};  // ssao/ssao_bitmask.frag:  slices with a visibility bitmask ("Quality")
        GLuint program_shadow[2] = {};      // shadows/screen_space.frag: a slice along the light, a bitmask over its chord
        GLuint program_temporal[2] = {};    // the filters, for both
        GLuint program_blur[2] = {};
        GLuint program_upsample[2] = {};
    } ssao;

    struct {
        GLuint fbo = 0;
        GLuint tex_half[2] = {};        // half res RGBA16F: far (filled) / near (premultiplied) layers for the composite
        GLuint tex_tile[2] = {};        // tile res RGBA16F, see apply_dof
        GLuint tex_src[3] = {};         // half res RGBA16F: near / far (mipmapped) / CoC (mipmapped) sources (dof_prepass)
        GLuint tex_far_pull = 0;        // half res RGBA16F, mipmapped: far layer as gathered, holes at foreground (dof_fill)
        GLuint tex_far_push = 0;        // quarter res RGBA16F, mipmapped: level i holds the filled far layer of pull level i + 1
        GLuint tex_band_src = 0;        // RGBA16F, level i: source of near band i + 1 (dof_near_down)
        GLuint tex_band_coc = 0;        // R16F, level i: its CoC
        GLuint tex_band_acc = 0;        // RGBA16F, level i: near field of bands i + 1 and up (dof_near)
        int src_width = 0;
        int src_height = 0;
        int pull_levels = 0;
        int band_top = 0;               // coarsest near band
        int band_width = 0;             // size of level 0 of the band textures (pyramid level 1, padded)
        int band_height = 0;
        GLuint program_prepass = 0;
        GLuint program_tile = 0;
        GLuint program_gather = 0;
        GLuint program_near_down = 0;
        GLuint program_near = 0;
        GLuint program_fill = 0;
        GLuint program_composite = 0;
    } bokeh_dof;

    struct {
        GLuint program = 0;
        struct {
            GLint mode = -1;
            GLint tex_color = -1;
        } uniform_loc;
    } tonemapping;

    struct {
        GLuint program = 0;
        struct {
            GLint tex_rgba = -1;
        } uniform_loc;
    } luma;

    struct {
        GLuint program = 0;
        struct {
            GLint tex_rgbl = -1;
            GLint rcp_res  = -1;
            GLint tc_scl   = -1;
        } uniform_loc;
    } fxaa;

    // The temporal resolve: [with motion blur][upscaling from a G-buffer smaller than the image shown]
    struct Gl_TemporalProgram {
        GLuint program = 0;
        struct {
            GLint tex_linear_depth = -1;
            GLint tex_main = -1;
            GLint tex_prev = -1;
            GLint tex_vel = -1;
            GLint tex_vel_neighbormax = -1;
            GLint texel_size = -1;
            GLint in_size = -1;
            GLint out_texel_size = -1;
            GLint time = -1;
            GLint feedback_min = -1;
            GLint feedback_max = -1;
            GLint motion_scale = -1;
            GLint jitter_uv = -1;
            GLint pixel_size = -1;
        } uniform_loc;
    } temporal[2][2];

    // An image to the output size, a frame on its own (Catmull-Rom): where the image is upscaled without TAA
    struct {
        GLuint program = 0;
        GLint  loc_tex = -1;
        GLint  loc_src_texel = -1;
        GLint  loc_src_size = -1;
    } upscale;
} gl;

static constexpr str_t v_shader_src_fs_quad = STR_LIT(
R"(
#version 150 core

out vec2 tc;

uniform vec2 u_tc_scl = vec2(1,1);

void main() {
	uint idx = uint(gl_VertexID) % 3U;
	gl_Position = vec4(
		(float( idx     &1U)) * 4.0 - 1.0,
		(float((idx>>1U)&1U)) * 4.0 - 1.0,
		0, 1.0);
	tc = (gl_Position.xy * 0.5 + 0.5) * u_tc_scl;
}
)");

static constexpr str_t f_shader_src_linearize_depth = STR_LIT(
R"(
#ifndef PERSPECTIVE
#define PERSPECTIVE 1
#endif

// z_n * z_f,  z_n - z_f,  z_f, *not used*
uniform vec4 u_clip_info;
uniform sampler2D u_tex_depth;

float ReconstructCSZ(float d, vec4 clip_info) {
#if PERSPECTIVE
    return (clip_info[0] / (d*clip_info[1] + clip_info[2]));
#else
    return (clip_info[1] + clip_info[2] - d*clip_info[1]);
#endif
}

out vec4 out_frag;

void main() {
  float d = texelFetch(u_tex_depth, ivec2(gl_FragCoord.xy), 0).x;
  out_frag = vec4(ReconstructCSZ(d, u_clip_info));
}
)");


static GLuint setup_program_from_source(str_t name, str_t f_shader_src, str_t defines = {}) {
    GLuint f_shader = gl::compile_shader_from_source(f_shader_src, GL_FRAGMENT_SHADER, defines);
    GLuint program = 0;

    if (f_shader) {
        char buffer[1024];
        program = glCreateProgram();

        glAttachShader(program, gl.v_shader_fs_quad);
        glAttachShader(program, f_shader);
        glLinkProgram(program);
        if (gl::get_program_link_error(buffer, sizeof(buffer), program)) {
            MD_LOG_ERROR("Error while linking %.*s program:\n%s", (int)name.len, name.ptr, buffer);
            glDeleteProgram(program);
            return 0;
        }

        glDetachShader(program, gl.v_shader_fs_quad);
        glDetachShader(program, f_shader);
        glDeleteShader(f_shader);
    }

    return program;
}

static void ensure_texture_2d(GLuint* tex, GLint internal_format, int width, int height, GLenum format, GLenum type, GLenum min_filter = GL_LINEAR, GLenum mag_filter = GL_LINEAR) {
    ASSERT(tex);
    if (!*tex) {
        glGenTextures(1, tex);
    }

    glBindTexture(GL_TEXTURE_2D, *tex);
    glTexImage2D(GL_TEXTURE_2D, 0, internal_format, width, height, 0, format, type, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min_filter);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag_filter);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);
}

// Texture with a full mip chain (filled by glGenerateMipmap)
static void ensure_texture_2d_mipmapped(GLuint* tex, GLint internal_format, int width, int height, GLenum format, GLenum type) {
    ASSERT(tex);
    if (!*tex) {
        glGenTextures(1, tex);
    }

    int levels = 1;
    while ((width >> levels) > 0 || (height >> levels) > 0) ++levels;

    glBindTexture(GL_TEXTURE_2D, *tex);
    for (int level = 0; level < levels; ++level) {
        glTexImage2D(GL_TEXTURE_2D, level, internal_format, MAX(width >> level, 1), MAX(height >> level, 1), 0, format, type, nullptr);
    }
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, levels - 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);
}

static void ensure_framebuffer(GLuint* fbo) {
    ASSERT(fbo);
    if (!*fbo) {
        glGenFramebuffers(1, fbo);
    }
}

// Best-effort GPU copy path. Drivers may schedule this on a copy engine when available.
static void copy_texture_2d(GLuint dst_tex, GLuint src_tex, int width, int height) {
    ASSERT(glIsTexture(dst_tex));
    ASSERT(glIsTexture(src_tex));

    if (glCopyImageSubData) {
        glCopyImageSubData(src_tex, GL_TEXTURE_2D, 0, 0, 0, 0,
                           dst_tex, GL_TEXTURE_2D, 0, 0, 0, 0,
                           width, height, 1);
        return;
    }

    ensure_framebuffer(&gl.rt.scratch_fbo);

    GLint prev_read_fbo = 0;
    GLint prev_tex = 0;
    glGetIntegerv(GL_READ_FRAMEBUFFER_BINDING, &prev_read_fbo);
    glGetIntegerv(GL_TEXTURE_BINDING_2D, &prev_tex);

    glBindFramebuffer(GL_READ_FRAMEBUFFER, gl.rt.scratch_fbo);
    glFramebufferTexture2D(GL_READ_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, src_tex, 0);

    glBindTexture(GL_TEXTURE_2D, dst_tex);
    glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, width, height);

    glBindFramebuffer(GL_READ_FRAMEBUFFER, prev_read_fbo);
    glBindTexture(GL_TEXTURE_2D, (GLuint)prev_tex);
}

namespace ssao {
// Scale-free SSAO, see shaders/ssao/ssao.frag. None of these constants depend on the scale of the scene: the sampling
// range is defined in screen space and the occlusion test is angular, so there is no world-space radius to tune.
static constexpr int   MAX_MIP        = 5;      // must match AO_MAX_MIP in ssao.frag (and SH_MAX_MIP in screen_space.frag)
static constexpr float R_MIN_PX       = 2.0f;   // smallest sampling distance, full-res pixels
static constexpr float R_MAX_FRACTION = 0.1f;   // largest sampling distance, fraction of the viewport height

// The screen-space shadows (shaders/shadows/screen_space.frag) reach further: a cavity is in the shadow of its rim
// as far below it as the light's slant carries it
static constexpr float SHADOW_R_MIN_PX       = 1.0f;
static constexpr float SHADOW_R_MAX_FRACTION = 0.3f;

static GLuint setup_variant(str_t name, str_t src, bool perspective) {
    return setup_program_from_source(name, src, perspective ? STR_LIT("#define AO_PERSPECTIVE 1") : STR_LIT("#define AO_PERSPECTIVE 0"));
}

static GLuint* all_programs[] = {gl.ssao.program_ao, gl.ssao.program_ao_quality, gl.ssao.program_shadow,
                                 gl.ssao.program_temporal, gl.ssao.program_blur, gl.ssao.program_upsample};

static void delete_programs() {
    for (GLuint* programs : all_programs) {
        for (int i = 0; i < 2; ++i) {
            if (programs[i]) glDeleteProgram(programs[i]);
            programs[i] = 0;
        }
    }
}

void initialize_programs() {
    const str_t ao_src       = {(const char*)ssao_frag, ssao_frag_size};
    const str_t quality_src  = {(const char*)ssao_bitmask_frag, ssao_bitmask_frag_size};
    const str_t shadow_src   = {(const char*)screen_space_frag, screen_space_frag_size};
    const str_t temporal_src = {(const char*)ssao_temporal_frag, ssao_temporal_frag_size};
    const str_t blur_src     = {(const char*)blur_frag, blur_frag_size};
    const str_t upsample_src = {(const char*)upsample_frag, upsample_frag_size};
    delete_programs();
    for (int i = 0; i < 2; ++i) {
        gl.ssao.program_ao[i]         = setup_variant(STR_LIT("ssao"),          ao_src,       i == 1);
        gl.ssao.program_ao_quality[i] = setup_variant(STR_LIT("ssao quality"),  quality_src,  i == 1);
        gl.ssao.program_shadow[i]     = setup_variant(STR_LIT("shadows"),       shadow_src,   i == 1);
        gl.ssao.program_temporal[i]   = setup_variant(STR_LIT("ssao temporal"), temporal_src, i == 1);
        gl.ssao.program_blur[i]       = setup_variant(STR_LIT("ssao blur"),     blur_src,     i == 1);
        gl.ssao.program_upsample[i]   = setup_variant(STR_LIT("ssao upsample"), upsample_src, i == 1);
    }
}

static void initialize_targets(Gl_Visibility* v, int width, int height) {
    // Must have the dimensions of mip level 1 of the linear depth texture.
    // RG32F: the depth channel feeds the bilateral weights and the history test, which work at the scale of a pixel
    // footprint. Half floats quantise depth to about 1e-3 of its value, more than a footprint at 4K, and the camera
    // jitter then moves depths across quantisation steps from frame to frame.
    const int half_w = MAX(width / 2, 1);
    const int half_h = MAX(height / 2, 1);
    for (int i = 0; i < 2; ++i) {
        ensure_texture_2d(&v->tex[i],         GL_RG32F, half_w, half_h, GL_RG, GL_FLOAT, GL_NEAREST, GL_NEAREST);
        ensure_texture_2d(&v->tex_history[i], GL_RG32F, half_w, half_h, GL_RG, GL_FLOAT, GL_NEAREST, GL_NEAREST);
    }
    ensure_texture_2d(&v->tex_visibility, GL_R16F, width, height, GL_RED, GL_FLOAT, GL_NEAREST, GL_NEAREST);
    v->history_valid = false;
}

static void shutdown(Gl_Visibility* v) {
    if (v->tex[0]) glDeleteTextures(2, v->tex);
    if (v->tex_history[0]) glDeleteTextures(2, v->tex_history);
    if (v->tex_visibility) glDeleteTextures(1, &v->tex_visibility);
    v->tex_visibility = 0;
    v->tex[0] = v->tex[1] = 0;
    v->tex_history[0] = v->tex_history[1] = 0;
    v->history_valid = false;
}

void initialize_targets(int width, int height) {
    initialize_targets(&gl.ssao.ao, width, height);
    initialize_targets(&gl.ssao.shadow, width, height);
}

void shutdown() {
    shutdown(&gl.ssao.ao);
    shutdown(&gl.ssao.shadow);
    delete_programs();
}

}  // namespace ssao

namespace fxaa {
void initialize() {
    gl.luma.program = setup_program_from_source(STR_LIT("luma"), {(const char*)luma_frag, luma_frag_size});
    gl.luma.uniform_loc.tex_rgba = glGetUniformLocation(gl.luma.program, "u_tex_rgba");

    str_t defines = STR_INIT("#define FXAA_PC 1\n#define FXAA_GLSL_130 1\n#define FXAA_QUALITY__PRESET 12");
    gl.fxaa.program = setup_program_from_source(STR_LIT("fxaa"), {(const char*)fxaa_frag, fxaa_frag_size}, defines);
    gl.fxaa.uniform_loc.tex_rgbl = glGetUniformLocation(gl.fxaa.program, "u_tex_rgbl");
    gl.fxaa.uniform_loc.rcp_res  = glGetUniformLocation(gl.fxaa.program, "u_rcp_res");
    gl.fxaa.uniform_loc.tc_scl   = glGetUniformLocation(gl.fxaa.program, "u_tc_scl");

}

void shutdown() {
    if (gl.luma.program) glDeleteProgram(gl.luma.program);
    if (gl.fxaa.program) glDeleteProgram(gl.fxaa.program);
}
}

namespace compose {

struct ubo_data_t {
    vec4_t proj_info;
    float  linear_colors;   // 1: ColorSpace_Linear
    float  light_strength;  // Settings::shading.light
    float  use_shadow;      // 1: u_texture_shadow holds the visibility of the key light
    float  time;
    vec3_t sky_radiance;
    float  roughness;
    vec3_t dir_radiance;
    float  F0;
    vec3_t light_dir;
    float  basic;           // 1: ShadingMode_Basic
    vec3_t ground_radiance;
    float  basic_ambient;
    float  basic_specular;
    float  use_ao;
    vec2_t inv_viewport;
};
static_assert(sizeof(ubo_data_t) == 112, "ubo_data_t must match the std140 layout of UniformData in compose_deferred.frag");
static_assert(offsetof(ubo_data_t, inv_viewport) == 104, "std140 aligns a vec2 to 8 bytes");

static struct {
    struct {
        GLuint program = 0;
        struct {
            GLint uniform_data = -1;
            GLint texture_depth = -1;
            GLint texture_color = -1;
            GLint texture_normal = -1;
            GLint texture_ao = -1;
            GLint texture_shadow = -1;
        } uniform_loc;
    } persp, ortho;
    GLuint ubo = 0;
} compose;

static void init_variant(decltype(compose.persp)& v, str_t defines) {
    v.program = setup_program_from_source(STR_LIT("compose deferred"), {(const char*)compose_deferred_frag, compose_deferred_frag_size}, defines);
    v.uniform_loc.texture_depth = glGetUniformLocation(v.program, "u_texture_depth");
    v.uniform_loc.texture_color = glGetUniformLocation(v.program, "u_texture_color");
    v.uniform_loc.texture_normal = glGetUniformLocation(v.program, "u_texture_normal");
    v.uniform_loc.texture_ao = glGetUniformLocation(v.program, "u_texture_ao");
    v.uniform_loc.texture_shadow = glGetUniformLocation(v.program, "u_texture_shadow");
    v.uniform_loc.uniform_data = glGetUniformBlockIndex(v.program, "UniformData");
}

void initialize() {
    init_variant(compose.persp, STR_LIT("#define COMPOSE_PROJECTION_PERSPECTIVE 1\n"));
    init_variant(compose.ortho, STR_LIT("#define COMPOSE_PROJECTION_PERSPECTIVE 0\n"));
    
    if (compose.ubo == 0) {
        glGenBuffers(1, &compose.ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, compose.ubo);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(ubo_data_t), 0, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
    }
}

void shutdown() {
    if (compose.persp.program) glDeleteProgram(compose.persp.program);
    if (compose.ortho.program) glDeleteProgram(compose.ortho.program);
    if (compose.ubo)        glDeleteBuffers(1, &compose.ubo);
    compose.ubo = 0;
}
}  // namespace compose

namespace tonemapping {

static struct {
    GLuint program = 0;
    struct {
        GLint texture = -1;
        GLint texture_coverage = -1;
        GLint texture_iso = -1;
        GLint use_iso = -1;
        GLint background = -1;
        GLint curve = -1;
        GLint exposure = -1;
        GLint gamma = -1;
    } uniform_loc;
} display_transform;

static struct {
    GLuint program_forward = 0;
    GLuint program_inverse = 0;
    struct {
        GLint texture = -1;
    } uniform_loc;
} fast_reversible;

void initialize() {
    {
        display_transform.program = setup_program_from_source(STR_LIT("Display transform"), {(const char*)display_frag, display_frag_size});
        display_transform.uniform_loc.texture          = glGetUniformLocation(display_transform.program, "u_texture");
        display_transform.uniform_loc.texture_coverage = glGetUniformLocation(display_transform.program, "u_texture_coverage");
        display_transform.uniform_loc.texture_iso      = glGetUniformLocation(display_transform.program, "u_texture_iso");
        display_transform.uniform_loc.use_iso          = glGetUniformLocation(display_transform.program, "u_use_iso");
        display_transform.uniform_loc.background       = glGetUniformLocation(display_transform.program, "u_background");
        display_transform.uniform_loc.curve            = glGetUniformLocation(display_transform.program, "u_curve");
        display_transform.uniform_loc.exposure         = glGetUniformLocation(display_transform.program, "u_exposure");
        display_transform.uniform_loc.gamma            = glGetUniformLocation(display_transform.program, "u_gamma");
    }
    {
        // Fast Reversible (For AA) (Credits to Brian Karis: http://graphicrants.blogspot.com/2013/12/tone-mapping.html)
        fast_reversible.program_forward = setup_program_from_source(STR_LIT("Fast Reversible"), {(const char*)fast_reversible_frag, fast_reversible_frag_size}, STR_LIT("#define USE_INVERSE 0"));
        fast_reversible.program_inverse = setup_program_from_source(STR_LIT("Fast Reversible"), {(const char*)fast_reversible_frag, fast_reversible_frag_size}, STR_LIT("#define USE_INVERSE 1"));
        fast_reversible.uniform_loc.texture = glGetUniformLocation(fast_reversible.program_forward, "u_texture");
    }
}

void shutdown() {
    if (display_transform.program) glDeleteProgram(display_transform.program);
    if (fast_reversible.program_forward) glDeleteProgram(fast_reversible.program_forward);
    if (fast_reversible.program_inverse) glDeleteProgram(fast_reversible.program_inverse);
}

}  // namespace tonemapping

namespace dof {
// Half-res depth of field:
//   prepass (half-res colour + CoC, gather sources) -> tiles -> far gather -> near bands -> far fill -> full-res composite
// The CoC is aperture * (1 - focus / z), with the aperture given as a fraction of the viewport height, which makes the
// blur independent of both the scene scale and the output resolution.
// Everything is deterministic (no noise for TAA to average), see dof_gather.frag and dof_near.frag.
static constexpr int   TILE             = 8;      // half-res pixels per tile side, must match TILE in the shaders
static constexpr float MAX_COC_FRACTION = 0.04f;  // CoC clamp, fraction of the viewport height
static constexpr int   MAX_BANDS        = 8;

static GLuint setup(str_t name, const unsigned char* src, size_t size, str_t defines = {}) {
    return setup_program_from_source(name, {(const char*)src, size}, defines);
}

// Coarsest near band for a viewport height: band L holds CoCs up to 2^(L + 1.5) half-res px (see dof_near.frag)
static int band_top(int height) {
    const float max_coc_half = MAX_COC_FRACTION * (float)height * 0.5f;
    const int top = (int)ceilf(log2f(MAX(max_coc_half, 1.0f)) - 1.5f);
    return CLAMP(top, 1, MAX_BANDS - 1);
}

// A texture with exactly the given number of levels (w, h halved per level)
static void ensure_texture_2d_levels(GLuint* tex, GLint internal_format, int width, int height, int levels, GLenum format, GLenum type, GLenum min_filter) {
    if (!*tex) glGenTextures(1, tex);
    glBindTexture(GL_TEXTURE_2D, *tex);
    for (int level = 0; level < levels; ++level) {
        glTexImage2D(GL_TEXTURE_2D, level, internal_format, MAX(width >> level, 1), MAX(height >> level, 1), 0, format, type, nullptr);
    }
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, levels - 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min_filter);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);
}

static GLuint* all_programs[] = {&gl.bokeh_dof.program_prepass, &gl.bokeh_dof.program_tile, &gl.bokeh_dof.program_gather, &gl.bokeh_dof.program_near_down,
                                 &gl.bokeh_dof.program_near, &gl.bokeh_dof.program_fill, &gl.bokeh_dof.program_composite};

void initialize_programs() {
    for (GLuint* prog : all_programs) {
        if (*prog) glDeleteProgram(*prog);
        *prog = 0;
    }
    gl.bokeh_dof.program_prepass   = setup(STR_LIT("DOF prepass"),   dof_prepass_frag,   dof_prepass_frag_size);
    gl.bokeh_dof.program_tile      = setup(STR_LIT("DOF tile"),      dof_tile_frag,      dof_tile_frag_size);
    gl.bokeh_dof.program_gather    = setup(STR_LIT("DOF gather"),    dof_gather_frag,    dof_gather_frag_size);
    gl.bokeh_dof.program_near_down = setup(STR_LIT("DOF near down"), dof_near_down_frag, dof_near_down_frag_size);
    gl.bokeh_dof.program_near      = setup(STR_LIT("DOF near"),      dof_near_frag,      dof_near_frag_size);
    gl.bokeh_dof.program_fill      = setup(STR_LIT("DOF fill"),      dof_fill_frag,      dof_fill_frag_size);
    gl.bokeh_dof.program_composite = setup(STR_LIT("DOF composite"), dof_composite_frag, dof_composite_frag_size);
}

void initialize_targets(int32_t width, int32_t height) {
    const int half_w = MAX(width / 2, 1);
    const int half_h = MAX(height / 2, 1);
    const int tile_w = DIV_UP(half_w, TILE);
    const int tile_h = DIV_UP(half_h, TILE);
    for (int i = 0; i < 2; ++i) {
        ensure_texture_2d(&gl.bokeh_dof.tex_half[i], GL_RGBA16F, half_w, half_h, GL_RGBA, GL_FLOAT, GL_LINEAR, GL_LINEAR);
    }
    for (int i = 0; i < 2; ++i) {
        ensure_texture_2d(&gl.bokeh_dof.tex_tile[i], GL_RGBA16F, tile_w, tile_h, GL_RGBA, GL_FLOAT, GL_NEAREST, GL_NEAREST);
    }
    ensure_texture_2d(&gl.bokeh_dof.tex_src[0], GL_RGBA16F, half_w, half_h, GL_RGBA, GL_FLOAT, GL_NEAREST, GL_NEAREST);
    for (int i = 1; i < 3; ++i) {
        ensure_texture_2d_mipmapped(&gl.bokeh_dof.tex_src[i], GL_RGBA16F, half_w, half_h, GL_RGBA, GL_FLOAT);
    }

    // Far fill
    ensure_texture_2d_mipmapped(&gl.bokeh_dof.tex_far_pull, GL_RGBA16F, half_w, half_h, GL_RGBA, GL_FLOAT);
    ensure_texture_2d_mipmapped(&gl.bokeh_dof.tex_far_push, GL_RGBA16F, MAX(half_w / 2, 1), MAX(half_h / 2, 1), GL_RGBA, GL_FLOAT);
    // The push reads one level of it at a time (base = max level), bilinear within that level
    glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.tex_far_push);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);
    int levels = 1;
    while ((half_w >> levels) > 0 || (half_h >> levels) > 0) ++levels;
    gl.bokeh_dof.pull_levels = levels;

    // Near bands 1..top at pyramid levels 1..top, padded to a multiple of the coarsest texel so every level covers
    // the whole viewport
    const int top = band_top(height);
    const int band_w = DIV_UP(half_w, 1 << top) << (top - 1);
    const int band_h = DIV_UP(half_h, 1 << top) << (top - 1);
    ensure_texture_2d_levels(&gl.bokeh_dof.tex_band_src, GL_RGBA16F, band_w, band_h, top, GL_RGBA, GL_FLOAT, GL_NEAREST_MIPMAP_NEAREST);
    ensure_texture_2d_levels(&gl.bokeh_dof.tex_band_coc, GL_R16F,    band_w, band_h, top, GL_RED,  GL_FLOAT, GL_NEAREST_MIPMAP_NEAREST);
    ensure_texture_2d_levels(&gl.bokeh_dof.tex_band_acc, GL_RGBA16F, band_w, band_h, top, GL_RGBA, GL_FLOAT, GL_LINEAR_MIPMAP_NEAREST);
    gl.bokeh_dof.band_top    = top;
    gl.bokeh_dof.band_width  = band_w;
    gl.bokeh_dof.band_height = band_h;

    gl.bokeh_dof.src_width  = half_w;
    gl.bokeh_dof.src_height = half_h;
}

void shutdown() {
    for (GLuint* prog : all_programs) {
        if (*prog) glDeleteProgram(*prog);
        *prog = 0;
    }
    GLuint* textures[] = {&gl.bokeh_dof.tex_half[0], &gl.bokeh_dof.tex_half[1], &gl.bokeh_dof.tex_tile[0], &gl.bokeh_dof.tex_tile[1],
                          &gl.bokeh_dof.tex_src[0], &gl.bokeh_dof.tex_src[1], &gl.bokeh_dof.tex_src[2], &gl.bokeh_dof.tex_far_pull,
                          &gl.bokeh_dof.tex_far_push, &gl.bokeh_dof.tex_band_src, &gl.bokeh_dof.tex_band_coc, &gl.bokeh_dof.tex_band_acc};
    for (GLuint* tex : textures) {
        if (*tex) glDeleteTextures(1, tex);
        *tex = 0;
    }
}
}  // namespace dof

namespace blit {
static GLuint program_tex = 0;
static GLuint program_tex_dither = 0;
static GLuint program_marks = 0;
static GLuint program_coverage = 0;
static GLint uniform_loc_texture = -1;

// The coverage of the geometry: 1 where the G-buffer holds a surface (its colour has alpha), 0 for the background
constexpr str_t f_shader_src_coverage = STR_LIT(R"(
#version 150 core

uniform sampler2D u_texture;

out vec4 out_frag;

void main() {
    out_frag = vec4(texelFetch(u_texture, ivec2(gl_FragCoord.xy), 0).a > 0.0 ? 1.0 : 0.0);
}
)");

constexpr str_t f_shader_src_tex = STR_LIT(R"(
#version 150 core

uniform sampler2D u_texture;

out vec4 out_frag;

void main() {
    out_frag = texelFetch(u_texture, ivec2(gl_FragCoord.xy), 0);
}
)");

// Final copy to the (8-bit) output with triangular dither of +-1 LSB. Quantisation happens only here, so this replaces
// the large per-pixel albedo noise compose used to add against banding.
constexpr str_t f_shader_src_tex_dither = STR_LIT(R"(
#version 150 core

uniform sampler2D u_texture;
uniform float u_time;

out vec4 out_frag;

float hash(vec2 p) {
    vec3 p3 = fract(vec3(p.xyx) * 0.1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}

void main() {
    vec4 c = texelFetch(u_texture, ivec2(gl_FragCoord.xy), 0);
    vec2 p = gl_FragCoord.xy + u_time * 61.0;
    float d = hash(p) + hash(p + vec2(17.31, 41.17)) - 1.0;   // triangular in [-1, 1]
    out_frag = vec4(c.rgb + d / 255.0, c.a);
}
)");

// The final copy with the marks of the selection and the highlight put over it (Inputs::marks, Settings::marks).
// The image is anti-aliased, the marks are not: they are a pixel soft, by a tent over the 3x3 pixels about each, which
// four bilinear taps between the pixels make. The outline is the band in which the tent of where the marks are at all
// is neither 0 nor 1: two pixels wide, and with taps further out (u_outline_step) as many times wider.
constexpr str_t f_shader_src_marks = STR_LIT(R"(
#version 410 core

uniform sampler2D u_texture;
uniform sampler2D u_marks;          // the output's size, filtered: selection in view, at all; highlight in view, at all
uniform sampler2D u_marks_group;    // the highlight's group: 1 + it, 8 more in view, 0 none (in 1/255)
uniform float u_time;

uniform vec4  u_visible[2];         // selection, highlight
uniform vec4  u_hidden[2];
uniform vec4  u_outline[2];
uniform vec3  u_groups[7];          // the highlight's groups after the first
uniform float u_highlight;          // how much of the highlight shows
uniform int   u_outline_step;

out vec4 out_frag;

float hash(vec2 p) {
    vec3 p3 = fract(vec3(p.xyx) * 0.1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}

// The mean of the four 2x2 blocks of pixels about this one, r pixels out diagonally
vec4 marks_about(vec2 px, float r, vec2 inv_size) {
    return 0.25 * (texture(u_marks, (px + vec2(-r, -r)) * inv_size) +
                   texture(u_marks, (px + vec2( r, -r)) * inv_size) +
                   texture(u_marks, (px + vec2(-r,  r)) * inv_size) +
                   texture(u_marks, (px + vec2( r,  r)) * inv_size));
}

// The highlight's group at the pixel, of the 3x3 about it where it has none (by its tent its marks reach that far):
// the one in view first, of those the last
int group_about(vec2 px, vec2 inv_size) {
    float g = texelFetch(u_marks_group, ivec2(px), 0).r;
    if (g == 0.0) {
        vec4 m = max(max(textureGather(u_marks_group, (px + vec2(-0.49, -0.49)) * inv_size),
                         textureGather(u_marks_group, (px + vec2( 0.49, -0.49)) * inv_size)),
                     max(textureGather(u_marks_group, (px + vec2(-0.49,  0.49)) * inv_size),
                         textureGather(u_marks_group, (px + vec2( 0.49,  0.49)) * inv_size)));
        g = max(max(m.x, m.y), max(m.z, m.w));
    }
    int code = int(g * 255.0 + 0.5);
    return code > 0 ? (code - 1) % 8 : 0;
}

void main() {
    vec4 c = texelFetch(u_texture, ivec2(gl_FragCoord.xy), 0);
    vec2 inv_size = 1.0 / vec2(textureSize(u_marks, 0));

    // The tent of the marks about the pixel, and for the outline, wider with the taps further out
    vec4 m = marks_about(gl_FragCoord.xy, 0.5, inv_size);
    vec4 ring = m;
    if (u_outline_step > 1) {
        ring = 0.5 * (m + marks_about(gl_FragCoord.xy, float(u_outline_step) - 0.5, inv_size));
    }

    if (max(max(ring.y, ring.w), max(m.y, m.w)) > 0.0) {
        vec2 vis    = m.xz;
        vec2 hidden = max(m.yw - m.xz, 0.0);
        vec2 edge   = clamp(2.0 - 2.0 * abs(2.0 * ring.yw - 1.0), 0.0, 1.0);

        // The highlight over the selection, in the colours of its group
        vec3 visible[2] = vec3[2](u_visible[0].rgb, u_visible[1].rgb);
        vec3 hidden_c[2] = vec3[2](u_hidden[0].rgb, u_hidden[1].rgb);
        vec3 outline[2] = vec3[2](u_outline[0].rgb, u_outline[1].rgb);
        if (max(ring.w, m.w) > 0.0) {
            int g = group_about(gl_FragCoord.xy, inv_size);
            if (g > 0) {
                visible[1] = hidden_c[1] = outline[1] = u_groups[g - 1];
            }
        }
        for (int i = 0; i < 2; ++i) {
            float k = i == 0 ? 1.0 : u_highlight;
            c.rgb = mix(c.rgb, visible[i],  u_visible[i].a * vis[i] * k);
            c.rgb = mix(c.rgb, hidden_c[i], u_hidden[i].a * hidden[i] * k);
            c.rgb = mix(c.rgb, outline[i],  u_outline[i].a * edge[i] * k);
        }
    }

    vec2 p = gl_FragCoord.xy + u_time * 61.0;
    float d = hash(p) + hash(p + vec2(17.31, 41.17)) - 1.0;   // triangular in [-1, 1]
    out_frag = vec4(c.rgb + d / 255.0, c.a);
}
)");

void initialize() {
    program_tex = setup_program_from_source(STR_LIT("blit texture"), f_shader_src_tex);
    uniform_loc_texture = glGetUniformLocation(program_tex, "u_texture");

    program_tex_dither = setup_program_from_source(STR_LIT("blit texture dither"), f_shader_src_tex_dither);

    program_marks = setup_program_from_source(STR_LIT("present marks"), f_shader_src_marks);

    program_coverage = setup_program_from_source(STR_LIT("blit coverage"), f_shader_src_coverage);

}

void shutdown() {
    if (program_tex) glDeleteProgram(program_tex);
    if (program_tex_dither) glDeleteProgram(program_tex_dither);
    if (program_marks) glDeleteProgram(program_marks);
    if (program_coverage) glDeleteProgram(program_coverage);
    program_tex = program_tex_dither = program_marks = program_coverage = 0;
}
}  // namespace blit

namespace blur {
static GLuint program_gaussian = 0;
static GLuint program_box = 0;
static GLint uniform_loc_texture = -1;
static GLint uniform_loc_inv_res_dir = -1;

constexpr str_t f_shader_src_gaussian = STR_LIT(R"(
#version 150 core

#define KERNEL_RADIUS 5

uniform sampler2D u_texture;
uniform vec2      u_inv_res_dir;

in vec2 tc;
out vec4 out_frag;

float blur_weight(float r) {
    const float sigma = KERNEL_RADIUS * 0.5;
    const float falloff = 1.0 / (2.0*sigma*sigma);
    float w = exp2(-r*r*falloff);
    return w;
}

void main() {
    vec2 uv = tc;
    vec4  c_tot = texture(u_texture, uv);
    float w_tot = 1.0;

    for (float r = 1; r <= KERNEL_RADIUS; ++r) {
        float w = blur_weight(r);
        vec4  c = texture(u_texture, uv + u_inv_res_dir * r);
        c_tot += c * w;
        w_tot += w;
    }
    for (float r = 1; r <= KERNEL_RADIUS; ++r) {
        float w = blur_weight(r);
        vec4  c = texture(u_texture, uv - u_inv_res_dir * r);
        c_tot += c * w;
        w_tot += w;
    }

    out_frag = c_tot / w_tot;
}
)");

constexpr str_t f_shader_src_box = STR_LIT(R"(
#version 150 core

uniform sampler2D u_texture;
out vec4 out_frag;

void main() {
    vec4 c = vec4(0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2(-1, -1), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2( 0, -1), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2(+1, -1), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2(-1,  0), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2( 0,  0), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2(+1,  0), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2(-1, +1), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2( 0, +1), 0);
    c += texelFetch(u_texture, ivec2(gl_FragCoord.xy) + ivec2(+1, +1), 0);

    out_frag = c / 9.0;
}
)");

void initialize() {
    program_gaussian = setup_program_from_source(STR_LIT("gaussian blur"), f_shader_src_gaussian);
    uniform_loc_texture = glGetUniformLocation(program_gaussian, "u_texture");
    uniform_loc_inv_res_dir = glGetUniformLocation(program_gaussian, "u_inv_res_dir");

    program_box = setup_program_from_source(STR_LIT("box blur"), f_shader_src_box);
}

void shutdown() {
    if (program_gaussian) glDeleteProgram(program_gaussian);
    if (program_box) glDeleteProgram(program_box);
}
}  // namespace blit

namespace velocity {
#define VEL_TILE_SIZE 8

struct {
    GLuint program = 0;
    struct {
		GLint tex_depth = -1;
        GLint curr_clip_to_prev_clip_mat = -1;
        GLint jitter_uv = -1;
    } uniform_loc;
} blit_velocity;

struct {
    GLuint program = 0;
    struct {
        GLint tex_vel = -1;
        GLint tex_linear_depth = -1;
    } uniform_loc;
} blit_tilemax;

struct {
    GLuint program = 0;
    struct {
        GLint tex_vel = -1;
        GLint tex_linear_depth = -1;
        GLint tex_vel_texel_size = -1;
    } uniform_loc;
} blit_neighbormax;

struct {
    GLuint program = 0;
    struct {
        GLint tex_vel = -1;
        GLint texel_size = -1;
    } uniform_loc;
} blit_dilate;

void initialize_programs() {
    {
        blit_velocity.program = setup_program_from_source(STR_LIT("screen-space velocity"), {(const char*)blit_velocity_frag, blit_velocity_frag_size});
		blit_velocity.uniform_loc.tex_depth = glGetUniformLocation(blit_velocity.program, "u_tex_depth");
        blit_velocity.uniform_loc.curr_clip_to_prev_clip_mat = glGetUniformLocation(blit_velocity.program, "u_curr_clip_to_prev_clip_mat");
        blit_velocity.uniform_loc.jitter_uv = glGetUniformLocation(blit_velocity.program, "u_jitter_uv");

    }
    {
        str_t defines = STR_INIT("#define TILE_SIZE " STRINGIFY_VAL(VEL_TILE_SIZE));
        blit_tilemax.program = setup_program_from_source(STR_LIT("tilemax"), {(const char*)blit_tilemax_frag, blit_tilemax_frag_size}, defines);
        blit_tilemax.uniform_loc.tex_vel = glGetUniformLocation(blit_tilemax.program, "u_tex_vel");
        blit_tilemax.uniform_loc.tex_linear_depth = glGetUniformLocation(blit_tilemax.program, "u_tex_linear_depth");
    }
    {
        blit_neighbormax.program = setup_program_from_source(STR_LIT("neighbormax"), {(const char*)blit_neighbormax_frag, blit_neighbormax_frag_size});
        blit_neighbormax.uniform_loc.tex_vel = glGetUniformLocation(blit_neighbormax.program, "u_tex_vel");
        blit_neighbormax.uniform_loc.tex_linear_depth = glGetUniformLocation(blit_neighbormax.program, "u_tex_linear_depth");
        blit_neighbormax.uniform_loc.tex_vel_texel_size = glGetUniformLocation(blit_neighbormax.program, "u_tex_vel_texel_size");
    }
    {
        blit_dilate.program = setup_program_from_source(STR_LIT("velocity dilate"), { (const char*)blit_velocity_dilate_frag, blit_velocity_dilate_frag_size });
        blit_dilate.uniform_loc.tex_vel = glGetUniformLocation(blit_dilate.program, "u_tex_vel");
        blit_dilate.uniform_loc.texel_size = glGetUniformLocation(blit_dilate.program, "u_texel_size");
    }
}

void initialize_targets(int32_t width, int32_t height) {
    if (!gl.velocity.tex_tilemax) {
        glGenTextures(1, &gl.velocity.tex_tilemax);
    }

    if (!gl.velocity.tex_neighbormax) {
        glGenTextures(1, &gl.velocity.tex_neighbormax);
    }

    // Clamp to at least one texel: a zero-sized attachment makes the FBO incomplete.
    gl.velocity.tex_width  = MAX(width  / VEL_TILE_SIZE, 1);
    gl.velocity.tex_height = MAX(height / VEL_TILE_SIZE, 1);

    glBindTexture(GL_TEXTURE_2D, gl.velocity.tex_tilemax);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16F, gl.velocity.tex_width, gl.velocity.tex_height, 0, GL_RG, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    glBindTexture(GL_TEXTURE_2D, gl.velocity.tex_neighbormax);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16F, gl.velocity.tex_width, gl.velocity.tex_height, 0, GL_RG, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    ensure_framebuffer(&gl.velocity.fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, gl.velocity.fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.velocity.tex_tilemax, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gl.velocity.tex_neighbormax, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void shutdown() {
    if (blit_velocity.program) glDeleteProgram(blit_velocity.program);
    if (blit_tilemax.program) glDeleteProgram(blit_tilemax.program);
    if (blit_neighbormax.program) glDeleteProgram(blit_neighbormax.program);
    if (blit_dilate.program) glDeleteProgram(blit_dilate.program);
    if (gl.velocity.tex_tilemax) glDeleteTextures(1, &gl.velocity.tex_tilemax);
    if (gl.velocity.tex_neighbormax) glDeleteTextures(1, &gl.velocity.tex_neighbormax);
    blit_velocity.program = blit_tilemax.program = blit_neighbormax.program = blit_dilate.program = 0;
    gl.velocity.tex_tilemax = gl.velocity.tex_neighbormax = 0;
}
}  // namespace velocity

namespace temporal {
static void get_locations(decltype(gl.temporal[0][0]) p) {  // a reference to one of them
    const GLuint prog = p.program;
    p.uniform_loc.tex_linear_depth    = glGetUniformLocation(prog, "u_tex_linear_depth");
    p.uniform_loc.tex_main            = glGetUniformLocation(prog, "u_tex_main");
    p.uniform_loc.tex_prev            = glGetUniformLocation(prog, "u_tex_prev");
    p.uniform_loc.tex_vel             = glGetUniformLocation(prog, "u_tex_vel");
    p.uniform_loc.tex_vel_neighbormax = glGetUniformLocation(prog, "u_tex_vel_neighbormax");
    p.uniform_loc.texel_size          = glGetUniformLocation(prog, "u_texel_size");
    p.uniform_loc.in_size             = glGetUniformLocation(prog, "u_in_size");
    p.uniform_loc.out_texel_size      = glGetUniformLocation(prog, "u_out_texel_size");
    p.uniform_loc.jitter_uv           = glGetUniformLocation(prog, "u_jitter_uv");
    p.uniform_loc.time                = glGetUniformLocation(prog, "u_time");
    p.uniform_loc.feedback_min        = glGetUniformLocation(prog, "u_feedback_min");
    p.uniform_loc.feedback_max        = glGetUniformLocation(prog, "u_feedback_max");
    p.uniform_loc.motion_scale        = glGetUniformLocation(prog, "u_motion_scale");
    p.uniform_loc.pixel_size          = glGetUniformLocation(prog, "u_pixel_size");
}

// Catmull-Rom from nine bilinear taps (the four of a row of weights folded into two where they are the same sign):
// sharp, but each frame on its own, for upscaling without TAA. Clamped at 0 for its overshoot.
static constexpr str_t f_shader_src_upscale = STR_LIT(R"(
#version 150 core
uniform sampler2D u_tex;
uniform vec4 u_src_texel;   // 1 / size, size of u_tex
uniform vec2 u_src_size;    // of the image in it: its corner (the texture can be larger, shared with larger views)
in vec2 tc;
out vec4 out_frag;

void main() {
    vec2 p   = tc * u_src_size;
    vec2 t1  = floor(p - 0.5) + 0.5;
    vec2 f   = p - t1;
    vec2 w0  = f * (-0.5 + f * (1.0 - 0.5 * f));
    vec2 w1  = 1.0 + f * f * (-2.5 + 1.5 * f);
    vec2 w2  = f * (0.5 + f * (2.0 - 1.5 * f));
    vec2 w3  = f * f * (-0.5 + 0.5 * f);
    vec2 w12 = w1 + w2;
    vec2 t0  = (t1 - 1.0) * u_src_texel.xy;
    vec2 t3  = (t1 + 2.0) * u_src_texel.xy;
    vec2 t12 = (t1 + w2 / w12) * u_src_texel.xy;
    vec4 c =
        (texture(u_tex, vec2(t0.x,  t0.y)) * w0.x + texture(u_tex, vec2(t12.x, t0.y)) * w12.x + texture(u_tex, vec2(t3.x,  t0.y)) * w3.x) * w0.y +
        (texture(u_tex, vec2(t0.x, t12.y)) * w0.x + texture(u_tex, vec2(t12.x, t12.y)) * w12.x + texture(u_tex, vec2(t3.x, t12.y)) * w3.x) * w12.y +
        (texture(u_tex, vec2(t0.x,  t3.y)) * w0.x + texture(u_tex, vec2(t12.x, t3.y)) * w12.x + texture(u_tex, vec2(t3.x,  t3.y)) * w3.x) * w3.y;
    out_frag = max(c, vec4(0.0));
}
)");

void initialize() {
    const str_t src = {(const char*)temporal_frag, temporal_frag_size};
    gl.temporal[1][0].program = setup_program_from_source(STR_LIT("temporal aa + motion-blur"), src);
    gl.temporal[0][0].program = setup_program_from_source(STR_LIT("temporal aa"), src, STR_LIT("#define USE_MOTION_BLUR 0\n"));
    gl.temporal[1][1].program = setup_program_from_source(STR_LIT("temporal aa upscale + motion-blur"), src, STR_LIT("#define UPSCALE 1\n"));
    gl.temporal[0][1].program = setup_program_from_source(STR_LIT("temporal aa upscale"), src, STR_LIT("#define USE_MOTION_BLUR 0\n#define UPSCALE 1\n"));
    for (auto& row : gl.temporal) {
        for (auto& p : row) get_locations(p);
    }

    gl.upscale.program       = setup_program_from_source(STR_LIT("upscale"), f_shader_src_upscale);
    gl.upscale.loc_tex       = glGetUniformLocation(gl.upscale.program, "u_tex");
    gl.upscale.loc_src_texel = glGetUniformLocation(gl.upscale.program, "u_src_texel");
    gl.upscale.loc_src_size  = glGetUniformLocation(gl.upscale.program, "u_src_size");
}

void shutdown() {
    for (auto& row : gl.temporal) {
        for (auto& p : row) {
            if (p.program) glDeleteProgram(p.program);
            p.program = 0;
        }
    }
    if (gl.upscale.program) glDeleteProgram(gl.upscale.program);
    gl.upscale.program = 0;
}
}  // namespace temporal

namespace sharpen {
static GLuint program = 0;
void initialize() {
    constexpr str_t f_shader_src_sharpen = STR_LIT(
 R"(#version 150 core

    uniform sampler2D u_tex;
    uniform float u_weight;
    out vec4 out_frag;

    void main() {
        vec3 cc = texelFetch(u_tex, ivec2(gl_FragCoord.xy), 0).rgb;
        vec3 cl = texelFetch(u_tex, ivec2(gl_FragCoord.xy) + ivec2(-1, 0), 0).rgb;
        vec3 ct = texelFetch(u_tex, ivec2(gl_FragCoord.xy) + ivec2( 0, 1), 0).rgb;
        vec3 cr = texelFetch(u_tex, ivec2(gl_FragCoord.xy) + ivec2( 1, 0), 0).rgb;
        vec3 cb = texelFetch(u_tex, ivec2(gl_FragCoord.xy) + ivec2( 0,-1), 0).rgb;

        float pos_weight = 1.0 + u_weight;
        float neg_weight = -u_weight * 0.25;
        out_frag = vec4(vec3(pos_weight * cc + neg_weight * (cl + ct + cr + cb)), 1.0);
    })");
    program = setup_program_from_source(STR_LIT("sharpen"), f_shader_src_sharpen);
}

void sharpen(GLuint in_texture, float weight = 1.0f) {
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, in_texture);

    glUseProgram(program);
    glUniform1i(glGetUniformLocation(sharpen::program, "u_tex"), 0);
    glUniform1f(glGetUniformLocation(sharpen::program, "u_weight"), weight);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);

    glUseProgram(0);
}

void shutdown() {
    if (program) glDeleteProgram(program);
}
}

static void initialize_programs() {
    if (!gl.vao) glGenVertexArrays(1, &gl.vao);

    gl.v_shader_fs_quad = gl::compile_shader_from_source(v_shader_src_fs_quad, GL_VERTEX_SHADER);

    gl.linear_depth.linearize.program_persp = setup_program_from_source(STR_LIT("linearize depth persp"), f_shader_src_linearize_depth, STR_LIT("#version 150 core\n#define PERSPECTIVE 1"));
    gl.linear_depth.linearize.program_ortho = setup_program_from_source(STR_LIT("linearize depth ortho"), f_shader_src_linearize_depth, STR_LIT("#version 150 core\n#define PERSPECTIVE 0"));
    gl.linear_depth.downsample.program = setup_program_from_source(STR_LIT("linear depth downsample"), {(const char*)depth_downsample_frag, depth_downsample_frag_size});

    ssao::initialize_programs();
    dof::initialize_programs();
    velocity::initialize_programs();
    tonemapping::initialize();
    temporal::initialize();
    blit::initialize();
    blur::initialize();
    sharpen::initialize();
    compose::initialize();
    fxaa::initialize();

    gl.programs_ready = true;
}

// Called at startup and on every resize. Programs used to be recompiled (and the old ones leaked) on every call;
// now they are created once and only the size dependent targets are reallocated here.
// width / height: the G-buffer's size; out_width / out_height: the image shown, when it is upscaled to it (0: the same)
void initialize(int width, int height, int out_width = 0, int out_height = 0) {
    out_width  = out_width  > 0 ? out_width  : width;
    out_height = out_height > 0 ? out_height : height;
    if (!gl.programs_ready) {
        initialize_programs();
    }

    if (gl.linear_depth.texture)
        glDeleteTextures(1, &gl.linear_depth.texture);
    
    {
        // R32F: R16F cannot hold the far plane (1e5 > 65504) and its ~11 bit mantissa quantises depth to about a pixel
        // footprint at Retina resolutions. Mip levels 1..MAX_MIP are rotated-grid subsamples used by the SSAO.
        int max_dim = MAX(width, height);
        int levels = 1;
        while ((max_dim >> levels) > 0 && levels < ssao::MAX_MIP + 1) ++levels;
        gl.linear_depth.levels = levels;

        glGenTextures(1, &gl.linear_depth.texture);
        glBindTexture(GL_TEXTURE_2D, gl.linear_depth.texture);
        glTexStorage2D(GL_TEXTURE_2D, levels, GL_R32F, width, height);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, levels - 1);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    ensure_framebuffer(&gl.rt.fbo);
    ensure_framebuffer(&gl.rt.scratch_fbo);
    ensure_framebuffer(&gl.rt.out_fbo);
    ensure_framebuffer(&gl.linear_depth.fbo);
    ensure_framebuffer(&gl.ssao.fbo);
    ensure_framebuffer(&gl.bokeh_dof.fbo);

    {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.linear_depth.fbo);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.linear_depth.texture, 0);
        GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            MD_LOG_ERROR("Something went wrong in creating framebuffer for depth linearization");
        }
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    }

    // Internal transient targets.
    // HDR stages use R11F_G11F_B10F, fine for linear HDR. Its 6/6/5 bit mantissas are too coarse for the display
    // referred values after tone mapping (steps of 1/128 and 1/64 near white: visible banding, and the TAA
    // exponential history can get stuck), so the LDR stages and the TAA history use RGB10_A2 at the same 32 bpp.
    ensure_texture_2d(&gl.rt.tex_color[0], GL_R11F_G11F_B10F, width, height, GL_RGB, GL_FLOAT);
    ensure_texture_2d(&gl.rt.tex_color[1], GL_R11F_G11F_B10F, width, height, GL_RGB, GL_FLOAT);
    ensure_texture_2d(&gl.rt.tex_ldr[0], GL_RGB10_A2, width, height, GL_RGBA, GL_UNSIGNED_INT_2_10_10_10_REV);
    ensure_texture_2d(&gl.rt.tex_ldr[1], GL_RGB10_A2, width, height, GL_RGBA, GL_UNSIGNED_INT_2_10_10_10_REV);
    ensure_texture_2d(&gl.rt.tex_coverage[0], GL_R8, width, height, GL_RED, GL_UNSIGNED_BYTE);
    ensure_texture_2d(&gl.rt.tex_coverage[1], GL_R8, width, height, GL_RED, GL_UNSIGNED_BYTE);
    ensure_texture_2d(&gl.rt.tex_history_prev, GL_RGB10_A2, out_width, out_height, GL_RGBA, GL_UNSIGNED_INT_2_10_10_10_REV);
    ensure_texture_2d(&gl.rt.tex_out[0], GL_RGB10_A2, out_width, out_height, GL_RGBA, GL_UNSIGNED_INT_2_10_10_10_REV);
    ensure_texture_2d(&gl.rt.tex_out[1], GL_RGB10_A2, out_width, out_height, GL_RGBA, GL_UNSIGNED_INT_2_10_10_10_REV);

    ensure_texture_2d(&gl.rt.tex_rgba8, GL_RGBA8, width, height, GL_RGBA, GL_UNSIGNED_BYTE);
    ensure_texture_2d(&gl.rt.tex_half_dof, GL_RGBA16F, MAX(width / 2, 1), MAX(height / 2, 1), GL_RGBA, GL_FLOAT);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.rt.fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.rt.tex_color[0], 0);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gl.rt.tex_color[1], 0);

    GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
    if (status != GL_FRAMEBUFFER_COMPLETE) {
        MD_LOG_ERROR("Something went wrong in creating framebuffer for targets");
    }

    GLenum buffers[] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
    glDrawBuffers(2, buffers);
    glClear(GL_COLOR_BUFFER_BIT);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

    gl.tex_width = width;
    gl.tex_height = height;
    gl.out_width = out_width;
    gl.out_height = out_height;

    ssao::initialize_targets(width, height);
    dof::initialize_targets(width, height);
    velocity::initialize_targets(width, height);
}

void shutdown() {
    ssao::shutdown();
    dof::shutdown();
    velocity::shutdown();
    tonemapping::shutdown();
    temporal::shutdown();
    blit::shutdown();
    blur::shutdown();
    sharpen::shutdown();
    compose::shutdown();
    fxaa::shutdown();

    if (gl.linear_depth.linearize.program_persp) glDeleteProgram(gl.linear_depth.linearize.program_persp);
    if (gl.linear_depth.linearize.program_ortho) glDeleteProgram(gl.linear_depth.linearize.program_ortho);
    if (gl.linear_depth.downsample.program)      glDeleteProgram(gl.linear_depth.downsample.program);
    if (gl.linear_depth.texture) glDeleteTextures(1, &gl.linear_depth.texture);
    gl.linear_depth.linearize.program_persp = gl.linear_depth.linearize.program_ortho = gl.linear_depth.downsample.program = 0;
    gl.linear_depth.texture = 0;

    if (gl.vao) glDeleteVertexArrays(1, &gl.vao);
    if (gl.v_shader_fs_quad) glDeleteShader(gl.v_shader_fs_quad);
    gl.vao = 0;
    gl.v_shader_fs_quad = 0;
    gl.programs_ready = false;
    if (gl.rt.fbo) glDeleteFramebuffers(1, &gl.rt.fbo);
    if (gl.rt.scratch_fbo) glDeleteFramebuffers(1, &gl.rt.scratch_fbo);
    if (gl.linear_depth.fbo) glDeleteFramebuffers(1, &gl.linear_depth.fbo);
    if (gl.velocity.fbo) glDeleteFramebuffers(1, &gl.velocity.fbo);
    if (gl.ssao.fbo) glDeleteFramebuffers(1, &gl.ssao.fbo);
    if (gl.bokeh_dof.fbo) glDeleteFramebuffers(1, &gl.bokeh_dof.fbo);
    if (gl.rt.tex_color[0]) glDeleteTextures(2, gl.rt.tex_color);
    if (gl.rt.tex_ldr[0]) glDeleteTextures(2, gl.rt.tex_ldr);
    if (gl.rt.tex_coverage[0]) glDeleteTextures(2, gl.rt.tex_coverage);
    gl.rt.tex_coverage[0] = gl.rt.tex_coverage[1] = 0;
    if (gl.rt.tex_history_prev) glDeleteTextures(1, &gl.rt.tex_history_prev);
    if (gl.rt.tex_rgba8) glDeleteTextures(1, &gl.rt.tex_rgba8);
    if (gl.rt.tex_half_dof) glDeleteTextures(1, &gl.rt.tex_half_dof);
}

void compute_linear_depth(GLuint depth_tex, float near_plane, float far_plane, bool orthographic = false) {
    const vec4_t clip_info {near_plane * far_plane, near_plane - far_plane, far_plane, 0};

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, depth_tex);

    GLuint program = orthographic ? gl.linear_depth.linearize.program_ortho : gl.linear_depth.linearize.program_persp;
    glUseProgram(program);
    glUniform1i(glGetUniformLocation (program, "u_tex_depth"), 0);
    glUniform4fv(glGetUniformLocation(program, "u_clip_info"), 1, &clip_info.x);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
}

void downsample_depth(GLuint linear_depth_tex, int src_lod) {
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);

    glUseProgram(gl.linear_depth.downsample.program);
    glUniform1i(glGetUniformLocation(gl.linear_depth.downsample.program, "u_tex_linear_depth"), 0);
    glUniform1i(glGetUniformLocation(gl.linear_depth.downsample.program, "u_src_lod"), src_lod);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
}

static vec4_t compute_proj_info(const float proj_mat[4][4]);

// What the passes of a screen-space visibility share for a frame: the view, the sizes, the units
struct VisibilityFrame {
    int width = 0;              // full res (the viewport)
    int height = 0;
    int half_w = 0;
    int half_h = 0;
    int variant = 0;            // 0 = orthographic, 1 = perspective
    vec4_t proj_info = {};
    float px_scale = 0;         // world size of a full-res pixel at depth 1 (persp) or absolute (ortho)
    float z_max = 0;            // background pixels hold the far clip distance: at or beyond this they are skipped
};

static VisibilityFrame visibility_frame(const GLResetState& reset_state, const ViewParam& view_param, float z_far) {
    VisibilityFrame f;
    f.width  = reset_state.viewport[2];
    f.height = reset_state.viewport[3];
    f.half_w = MAX(f.width  / 2, 1);
    f.half_h = MAX(f.height / 2, 1);
    const float (*proj_mat)[4] = view_param.matrix.curr.proj.elem;
    f.variant   = is_orthographic_proj_matrix(proj_mat) ? 0 : 1;
    f.proj_info = compute_proj_info(proj_mat);
    f.px_scale  = f.proj_info.y / (float)f.height;
    f.z_max     = z_far * 0.99f;
    return f;
}

// The history of v is only usable if it was written in the directly preceding frame of this view, at this size and
// with this estimator (switching modes shows the other estimator at once instead of blending into it)
static bool visibility_history_valid(const Gl_Visibility& v, const VisibilityFrame& f, const ViewParam& view_param, bool temporal, int mode) {
    return temporal && v.history_valid &&
        v.history_width == f.half_w && v.history_height == f.half_h && v.history_variant == f.variant &&
        v.history_mode == mode &&
        memcmp(&v.history_view, &view_param.matrix.prev.view, sizeof(mat4_t)) == 0 &&
        memcmp(&v.history_proj, &view_param.matrix.prev.proj, sizeof(mat4_t)) == 0;
}

// Uniforms every estimator takes, with its program bound and the linear depth on unit 0, the normals on unit 1
static void set_estimator_uniforms(GLuint program, const VisibilityFrame& f, float r_min, float r_max, int frame) {
    glUniform1i(glGetUniformLocation(program, "u_tex_linear_depth"), 0);
    glUniform1i(glGetUniformLocation(program, "u_tex_normal"), 1);
    glUniform4fv(glGetUniformLocation(program, "u_proj_info"), 1, &f.proj_info.x);
    glUniform2f(glGetUniformLocation(program, "u_full_res"), (float)f.width, (float)f.height);
    glUniform1f(glGetUniformLocation(program, "u_px_scale"), f.px_scale);
    glUniform1f(glGetUniformLocation(program, "u_r_min"), r_min);
    glUniform1f(glGetUniformLocation(program, "u_r_max"), r_max);
    glUniform1f(glGetUniformLocation(program, "u_z_max"), f.z_max);
    glUniform1i(glGetUniformLocation(program, "u_frame"), frame);
}

// Binds the half-res target of the estimator of v, and its inputs: the linear depth on unit 0, the normals on unit 1
static void begin_estimator(Gl_Visibility* v, const VisibilityFrame& f, GLuint linear_depth_tex, GLuint normal_tex) {
    glBindVertexArray(gl.vao);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.ssao.fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, v->tex[0], 0);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glViewport(0, 0, f.half_w, f.half_h);
    glScissor(0, 0, f.half_w, f.half_h);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, normal_tex);
}

// What follows the estimator of a visibility, which has written (visibility, linear depth) into v->tex[0] at half res:
// with temporal set (the estimator's pattern permuted by the frame), the accumulation in the history
// (shaders/ssao/ssao_temporal.frag), reprojected with velocity_tex; then a separable bilateral blur (blur.frag) and
// a joint bilateral upsample (upsample.frag) to full res, into v->tex_visibility, which the compose pass applies.
// Without temporal the pattern is fixed and the result stable.
static void filter_visibility(Gl_Visibility* v, const VisibilityFrame& f, GLuint linear_depth_tex, GLuint normal_tex, GLuint velocity_tex,
                              const ViewParam& view_param, bool temporal, bool history_valid, int mode) {
    const int variant = f.variant;

    glBindVertexArray(gl.vao);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.ssao.fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glViewport(0, 0, f.half_w, f.half_h);
    glScissor(0, 0, f.half_w, f.half_h);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, normal_tex);

    // Index into v->tex of the latest result
    int src = 0;

    if (temporal) {
        PUSH_GPU_SECTION("Temporal")
        const int hist_prev = v->history_curr;
        const int hist_curr = hist_prev ^ 1;

        // This frame's view space -> the previous frame's. The AO passes use view space with z into the screen (positive
        // linear depth), i.e. the GL view space mirrored in z: S * M * S with S = diag(1, 1, -1, 1).
        mat4_t curr_to_prev = view_param.matrix.prev.view * view_param.matrix.inv.view;
        for (int c = 0; c < 4; ++c) {
            for (int r = 0; r < 4; ++r) {
                if ((c == 2) != (r == 2)) curr_to_prev.elem[c][r] = -curr_to_prev.elem[c][r];
            }
        }
        const vec4_t proj_info_prev = compute_proj_info(view_param.matrix.prev.proj.elem);
        const vec2_t jitter_delta   = view_param.jitter.curr - view_param.jitter.prev;

        const GLuint program = gl.ssao.program_temporal[variant];
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, v->tex_history[hist_curr], 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, v->tex[1], 0);
        const GLenum buffers[2] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
        glDrawBuffers(2, buffers);

        glUseProgram(program);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, v->tex[0]);
        // normal_tex is bound to unit 1
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, v->tex_history[hist_prev]);
        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_2D, velocity_tex);
        glUniform1i(glGetUniformLocation(program, "u_tex_ao"), 0);
        glUniform1i(glGetUniformLocation(program, "u_tex_normal"), 1);
        glUniform1i(glGetUniformLocation(program, "u_tex_history"), 2);
        glUniform1i(glGetUniformLocation(program, "u_tex_velocity"), 3);
        glUniform4fv(glGetUniformLocation(program, "u_proj_info"), 1, &f.proj_info.x);
        glUniform4fv(glGetUniformLocation(program, "u_proj_info_prev"), 1, &proj_info_prev.x);
        glUniformMatrix4fv(glGetUniformLocation(program, "u_curr_to_prev"), 1, GL_FALSE, &curr_to_prev.elem[0][0]);
        glUniform2f(glGetUniformLocation(program, "u_jitter_delta"), jitter_delta.x, jitter_delta.y);
        glUniform2f(glGetUniformLocation(program, "u_full_res"), (float)f.width, (float)f.height);
        glUniform1f(glGetUniformLocation(program, "u_px_scale"), f.px_scale);
        glUniform1f(glGetUniformLocation(program, "u_z_max"), f.z_max);
        glUniform1i(glGetUniformLocation(program, "u_history_valid"), history_valid ? 1 : 0);
        glDrawArrays(GL_TRIANGLES, 0, 3);

        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, 0, 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        glActiveTexture(GL_TEXTURE0);

        v->history_curr    = hist_curr;
        v->history_valid   = true;
        v->history_width   = f.half_w;
        v->history_height  = f.half_h;
        v->history_variant = variant;
        v->history_mode    = mode;
        v->history_view    = view_param.matrix.curr.view;
        v->history_proj    = view_param.matrix.curr.proj;
        src = 1;
        POP_GPU_SECTION()
    } else {
        v->history_valid = false;
    }

    PUSH_GPU_SECTION("Blur")
    {
        const GLuint program = gl.ssao.program_blur[variant];
        glUseProgram(program);
        glUniform1i(glGetUniformLocation(program, "u_tex"), 0);
        glUniform1i(glGetUniformLocation(program, "u_tex_normal"), 1);
        glUniform4fv(glGetUniformLocation(program, "u_proj_info"), 1, &f.proj_info.x);
        glUniform2f(glGetUniformLocation(program, "u_full_res"), (float)f.width, (float)f.height);
        glUniform1f(glGetUniformLocation(program, "u_px_scale"), f.px_scale);
        glUniform1f(glGetUniformLocation(program, "u_z_max"), f.z_max);
        // The temporal filter already averages every pixel over the whole pattern; the blur then only has to cover
        // one period for the frames right after a disocclusion
        glUniform1i(glGetUniformLocation(program, "u_radius"), temporal ? 2 : 4);
        const GLint loc_dir = glGetUniformLocation(program, "u_dir");

        // normal_tex is still bound to unit 1
        glActiveTexture(GL_TEXTURE0);
        for (int pass = 0; pass < 2; ++pass) {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, v->tex[src ^ 1], 0);
            glBindTexture(GL_TEXTURE_2D, v->tex[src]);
            glUniform2i(loc_dir, pass == 0 ? 1 : 0, pass == 0 ? 0 : 1);
            glDrawArrays(GL_TRIANGLES, 0, 3);
            src ^= 1;
        }
    }
    POP_GPU_SECTION()

    // Into v->tex_visibility at full res, for the compose pass to apply to the light it occludes
    PUSH_GPU_SECTION("Upsample")
    {
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, v->tex_visibility, 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        glViewport(0, 0, f.width, f.height);
        glScissor(0, 0, f.width, f.height);
        const GLuint program = gl.ssao.program_upsample[variant];
        glUseProgram(program);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, v->tex[src]);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, linear_depth_tex);
        glUniform1i(glGetUniformLocation(program, "u_tex_ao"), 0);
        glUniform1i(glGetUniformLocation(program, "u_tex_linear_depth"), 1);
        glUniform1f(glGetUniformLocation(program, "u_px_scale"), f.px_scale);
        glUniform1f(glGetUniformLocation(program, "u_z_max"), f.z_max);
        glDrawArrays(GL_TRIANGLES, 0, 3);
    }
    POP_GPU_SECTION()

    glActiveTexture(GL_TEXTURE0);
}

// Half-res scale-free SSAO, upsampled to full res into gl.ssao.ao.tex_visibility, which the compose pass applies to the
// light of the surroundings. z_far is the far clip distance: background pixels hold exactly that linear depth and are
// skipped. With temporal set, the sample pattern is permuted every frame and accumulated in a history, reprojected
// with velocity_tex (filter_visibility). mode selects the estimator; everything after it is shared.
void compute_ssao(GLuint linear_depth_tex, GLuint normal_tex, GLuint velocity_tex, const ViewParam& view_param, float z_far, float intensity, bool temporal,
                  postprocess_pipeline::SsaoMode mode) {
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(normal_tex));
    ASSERT(!temporal || glIsTexture(velocity_tex));

    GLResetState reset_state = {};
    record_gl_reset_state(&reset_state);

    Gl_Visibility* v = &gl.ssao.ao;
    const VisibilityFrame f = visibility_frame(reset_state, view_param, z_far);
    const float r_max = MAX(ssao::R_MAX_FRACTION * (float)f.height, 2.0f * ssao::R_MIN_PX);
    const bool history_valid = visibility_history_valid(*v, f, view_param, temporal, (int)mode);
    const int frame = temporal ? (int)(v->frame++ & 15) : 0;

    PUSH_GPU_SECTION("AO")
    {
        const GLuint program = mode == postprocess_pipeline::SsaoMode_Quality ? gl.ssao.program_ao_quality[f.variant] : gl.ssao.program_ao[f.variant];
        begin_estimator(v, f, linear_depth_tex, normal_tex);
        glUseProgram(program);
        set_estimator_uniforms(program, f, ssao::R_MIN_PX, r_max, frame);
        glUniform1f(glGetUniformLocation(program, "u_intensity"), MAX(intensity, 0.0f));
        glDrawArrays(GL_TRIANGLES, 0, 3);
    }
    POP_GPU_SECTION()

    filter_visibility(v, f, linear_depth_tex, normal_tex, velocity_tex, view_param, temporal, history_valid, (int)mode);

    reset_gl_state(reset_state);
    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(0);
}

// Half-res screen-space shadows of the key light (shaders/shadows/screen_space.frag), filtered as the ambient occlusion
// is (filter_visibility), into gl.ssao.shadow.tex_visibility, which the compose pass applies to the key light.
// light_dir: toward the light, in the view space (GL: z toward the viewer); light_radius: its angular radius (radians).
void compute_shadows(GLuint linear_depth_tex, GLuint normal_tex, GLuint velocity_tex, const ViewParam& view_param, float z_far, vec3_t light_dir, float light_radius, bool temporal) {
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(normal_tex));
    ASSERT(!temporal || glIsTexture(velocity_tex));

    GLResetState reset_state = {};
    record_gl_reset_state(&reset_state);

    Gl_Visibility* v = &gl.ssao.shadow;
    const VisibilityFrame f = visibility_frame(reset_state, view_param, z_far);
    const float r_max = MAX(ssao::SHADOW_R_MAX_FRACTION * (float)f.height, 2.0f * ssao::SHADOW_R_MIN_PX);
    const bool history_valid = visibility_history_valid(*v, f, view_param, temporal, 0);
    const int frame = temporal ? (int)(v->frame++ & 15) : 0;

    // The passes' view space has z into the screen
    const vec3_t L = vec3_normalize(vec3_t{light_dir.x, light_dir.y, -light_dir.z});

    PUSH_GPU_SECTION("Shadow Rays")
    {
        const GLuint program = gl.ssao.program_shadow[f.variant];
        begin_estimator(v, f, linear_depth_tex, normal_tex);
        glUseProgram(program);
        set_estimator_uniforms(program, f, ssao::SHADOW_R_MIN_PX, r_max, frame);
        glUniform3f(glGetUniformLocation(program, "u_light_dir"), L.x, L.y, L.z);
        glUniform1f(glGetUniformLocation(program, "u_light_radius"), MAX(light_radius, 1.0e-3f));
        glDrawArrays(GL_TRIANGLES, 0, 3);
    }
    POP_GPU_SECTION()

    filter_visibility(v, f, linear_depth_tex, normal_tex, velocity_tex, view_param, temporal, history_valid, 0);

    reset_gl_state(reset_state);
    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(0);
}

// View position = ((uv * xy + zw) * z, z) (perspective) or (uv * xy + zw, z) (orthographic), z = linear depth.
// The y offset terms carry the same sign as the x ones. The NVIDIA original flips it, which only shows with an
// off-centre frustum: with the TAA jitter that misplaced every reconstructed position by twice the y jitter.
static vec4_t compute_proj_info(const float proj_mat[4][4]) {
    if (!is_orthographic_proj_matrix(proj_mat)) {
        return vec4_t{
            2.0f / (proj_mat[0][0]),
            2.0f / (proj_mat[1][1]),
            -(1.0f - proj_mat[2][0]) / proj_mat[0][0],
            -(1.0f - proj_mat[2][1]) / proj_mat[1][1]
        };
    }

    return vec4_t{
        2.0f / (proj_mat[0][0]),
        2.0f / (proj_mat[1][1]),
        -(1.0f + proj_mat[3][0]) / proj_mat[0][0],
        -(1.0f + proj_mat[3][1]) / proj_mat[1][1]
    };
}

// ao_tex: the visibility of the surroundings, the ambient occlusion (0: none); shadow_tex: the visibility of the key
// light, its shadows (0: none)
static void compose_deferred(GLuint linear_depth_tex, GLuint color_tex, GLuint normal_tex, GLuint ao_tex, GLuint shadow_tex, const float proj_mat[4][4], bool orthographic, int width, int height, const postprocess_pipeline::Settings& settings, float time) {
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(color_tex));
    ASSERT(glIsTexture(normal_tex));

    const bool basic = settings.shading.mode == postprocess_pipeline::ShadingMode_Basic;


    const postprocess_pipeline::Lighting light = postprocess_pipeline::lighting(settings.shading.mode);
    const float roughness = 0.4f;
    const float F0 = 0.04f;

    compose::ubo_data_t data = {
        .proj_info = compute_proj_info(proj_mat),
        .linear_colors = settings.shading.colors == postprocess_pipeline::ColorSpace_Linear ? 1.0f : 0.0f,
        .light_strength = MAX(settings.shading.light, 0.0f),
        .use_shadow = shadow_tex ? 1.0f : 0.0f,
        .time = time,
        .sky_radiance = light.sky_radiance * MAX(settings.shading.environment, 0.0f),
        .roughness = roughness,
        .dir_radiance = light.light_irradiance,
        .F0 = F0,
        .light_dir = light.light_dir,
        .basic = basic ? 1.0f : 0.0f,
        .ground_radiance = light.ground_radiance * MAX(settings.shading.environment, 0.0f),
        .basic_ambient = CLAMP(settings.shading.ambient, 0.0f, 1.0f),
        .basic_specular = MAX(settings.shading.specular, 0.0f),
        .use_ao = ao_tex ? 1.0f : 0.0f,
        .inv_viewport = {1.0f / (float)MAX(width, 1), 1.0f / (float)MAX(height, 1)},
    };

    glBindBuffer(GL_UNIFORM_BUFFER, compose::compose.ubo);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(data), &data);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, color_tex);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, normal_tex);
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, ao_tex ? ao_tex : gl.ssao.ao.tex_visibility);  // a valid texture either way: use_ao gates it
    glActiveTexture(GL_TEXTURE4);
    glBindTexture(GL_TEXTURE_2D, shadow_tex ? shadow_tex : gl.ssao.shadow.tex_visibility);  // use_shadow gates it
    glActiveTexture(GL_TEXTURE0);

    auto& variant = orthographic ? compose::compose.ortho : compose::compose.persp;
    GLuint program = variant.program;

    glUseProgram(program);

    glBindBufferBase(GL_UNIFORM_BUFFER, 0, compose::compose.ubo);
    glUniformBlockBinding(program, variant.uniform_loc.uniform_data, 0);
    glUniform1i (variant.uniform_loc.texture_depth, 0);
    glUniform1i (variant.uniform_loc.texture_color, 1);
    glUniform1i (variant.uniform_loc.texture_normal, 2);
    glUniform1i (variant.uniform_loc.texture_ao, 3);
    glUniform1i (variant.uniform_loc.texture_shadow, 4);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);

    glUseProgram(0);
}

// Depth of field from color_tex into the currently bound color target.
// aperture is the CoC of an object at infinity as a fraction of the viewport height (0.01 = 1% of the view height).
void apply_dof(GLuint linear_depth_tex, GLuint color_tex, float focus_depth, float aperture) {
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(color_tex));

    GLResetState reset_state = {};
    record_gl_reset_state(&reset_state);

    const int width  = reset_state.viewport[2];
    const int height = reset_state.viewport[3];
    const int half_w = MAX(width  / 2, 1);
    const int half_h = MAX(height / 2, 1);
    const int tile_w = DIV_UP(half_w, dof::TILE);
    const int tile_h = DIV_UP(half_h, dof::TILE);

    const float aperture_px = MAX(aperture, 0.0f) * (float)height;   // full-res px
    const float max_coc_px  = dof::MAX_COC_FRACTION * (float)height;
    const float focus       = MAX(focus_depth, 1.0e-3f);

    glBindVertexArray(gl.vao);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.bokeh_dof.fbo);

    PUSH_GPU_SECTION("DOF prepass")
    {
        // Half-res colour + CoC, and the mipmapped sources the gather integrates over
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.rt.tex_half_dof, 0);
        for (int i = 0; i < 3; ++i) {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1 + i, GL_TEXTURE_2D, gl.bokeh_dof.tex_src[i], 0);
        }
        const GLenum buffers[4] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3};
        glDrawBuffers(4, buffers);

        // The mip chains also average texels outside the viewport (when the targets are larger than it): those hold
        // nothing that spreads
        const GLfloat zero[4] = {0, 0, 0, 0};
        glViewport(0, 0, gl.bokeh_dof.src_width, gl.bokeh_dof.src_height);
        glScissor(0, 0, gl.bokeh_dof.src_width, gl.bokeh_dof.src_height);
        for (int i = 1; i < 4; ++i) {
            glClearBufferfv(GL_COLOR, i, zero);
        }

        const GLuint program = gl.bokeh_dof.program_prepass;
        glViewport(0, 0, half_w, half_h);
        glScissor(0, 0, half_w, half_h);
        glUseProgram(program);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, color_tex);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, linear_depth_tex);
        glUniform1i(glGetUniformLocation(program, "u_tex_color"), 0);
        glUniform1i(glGetUniformLocation(program, "u_tex_linear_depth"), 1);
        glUniform1f(glGetUniformLocation(program, "u_focus"), focus);
        glUniform1f(glGetUniformLocation(program, "u_aperture_px"), aperture_px * 0.5f);   // half-res px
        glUniform1f(glGetUniformLocation(program, "u_max_coc_px"), max_coc_px * 0.5f);
        glDrawArrays(GL_TRIANGLES, 0, 3);

        for (int i = 0; i < 3; ++i) {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1 + i, GL_TEXTURE_2D, 0, 0);
        }
        glDrawBuffer(GL_COLOR_ATTACHMENT0);

        // The far gather integrates over the mip chains (the near bands make their own levels)
        glActiveTexture(GL_TEXTURE0);
        for (int i = 1; i < 3; ++i) {
            glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.tex_src[i]);
            glGenerateMipmap(GL_TEXTURE_2D);
        }
    }
    POP_GPU_SECTION()

    PUSH_GPU_SECTION("DOF tiles")
    {
        // tex_tile[0]: (max near |CoC|, max far CoC, min near CoC) per tile
        // tex_tile[1]: (far gather radius, near-field reach, min near CoC of the 3x3 tiles around)
        const GLuint program = gl.bokeh_dof.program_tile;
        const GLint loc_dilate = glGetUniformLocation(program, "u_dilate");
        // Foreground blur reaches up to the CoC clamp; the dilation has to search that far (in tiles)
        const int reach_tiles = (int)ceilf(max_coc_px * 0.5f / (float)dof::TILE);
        glViewport(0, 0, tile_w, tile_h);
        glScissor(0, 0, tile_w, tile_h);
        glUseProgram(program);
        glUniform1i(glGetUniformLocation(program, "u_tex"), 0);
        glUniform1i(glGetUniformLocation(program, "u_tex_near_src"), 1);
        glUniform1i(glGetUniformLocation(program, "u_tex_coc_src"), 2);
        glUniform1i(glGetUniformLocation(program, "u_reach_tiles"), reach_tiles);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.tex_src[0]);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.tex_src[2]);
        glActiveTexture(GL_TEXTURE0);

        const GLuint dst[2] = {gl.bokeh_dof.tex_tile[0], gl.bokeh_dof.tex_tile[1]};
        const GLuint src[2] = {gl.rt.tex_half_dof,       gl.bokeh_dof.tex_tile[0]};
        for (int pass = 0; pass < 2; ++pass) {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, dst[pass], 0);
            glBindTexture(GL_TEXTURE_2D, src[pass]);
            glUniform1i(loc_dilate, pass);
            glDrawArrays(GL_TRIANGLES, 0, 3);
        }
    }
    POP_GPU_SECTION()

    PUSH_GPU_SECTION("DOF gather")
    {
        // Far field, holes where the centre belongs to the foreground (filled by the fill pass)
        glViewport(0, 0, half_w, half_h);
        glScissor(0, 0, half_w, half_h);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.bokeh_dof.tex_far_pull, 0);

        const GLuint program = gl.bokeh_dof.program_gather;
        glUseProgram(program);
        const GLuint textures[4] = {gl.rt.tex_half_dof, gl.bokeh_dof.tex_tile[1], gl.bokeh_dof.tex_src[1], gl.bokeh_dof.tex_src[2]};
        const char* names[4] = {"u_tex", "u_tex_tile", "u_tex_far_src", "u_tex_coc_src"};
        for (int i = 0; i < 4; ++i) {
            glActiveTexture(GL_TEXTURE0 + i);
            glBindTexture(GL_TEXTURE_2D, textures[i]);
            glUniform1i(glGetUniformLocation(program, names[i]), i);
        }
        glDrawArrays(GL_TRIANGLES, 0, 3);
    }
    POP_GPU_SECTION()

    PUSH_GPU_SECTION("DOF near")
    {
        // Near field in bands of CoC, band L gathered at pyramid level L, see dof_near.frag
        const int top = gl.bokeh_dof.band_top;
        {
            // Band sources, level L of the pyramid in level L - 1 of the band textures
            const GLuint program = gl.bokeh_dof.program_near_down;
            glUseProgram(program);
            glUniform1i(glGetUniformLocation(program, "u_tex_near_src"), 0);
            glUniform1i(glGetUniformLocation(program, "u_tex_coc_src"), 1);
            glUniform1i(glGetUniformLocation(program, "u_top"), top);
            const GLint loc_level = glGetUniformLocation(program, "u_level");
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.tex_src[0]);
            glActiveTexture(GL_TEXTURE1);
            glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.tex_src[2]);
            const GLenum two_buffers[2] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
            glDrawBuffers(2, two_buffers);
            for (int level = 1; level <= top; ++level) {
                glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.bokeh_dof.tex_band_src, level - 1);
                glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gl.bokeh_dof.tex_band_coc, level - 1);
                const int w = gl.bokeh_dof.band_width  >> (level - 1);
                const int h = gl.bokeh_dof.band_height >> (level - 1);
                glViewport(0, 0, w, h);
                glScissor(0, 0, w, h);
                glUniform1i(loc_level, level);
                glDrawArrays(GL_TRIANGLES, 0, 3);
            }
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, 0, 0);
            glDrawBuffer(GL_COLOR_ATTACHMENT0);
        }
        {
            // Gather each band and add the levels above, coarsest first; the last pass writes the half-res near layer.
            // Level L is written to level L - 1 of tex_band_acc while level L of it (the levels above) is read: the
            // sampled levels are restricted to that one.
            const GLuint program = gl.bokeh_dof.program_near;
            glUseProgram(program);
            glUniform1i(glGetUniformLocation(program, "u_tex_src"), 0);
            glUniform1i(glGetUniformLocation(program, "u_tex_coc"), 1);
            glUniform1i(glGetUniformLocation(program, "u_tex_up"), 2);
            glUniform1i(glGetUniformLocation(program, "u_tex_tile"), 3);
            glUniform1i(glGetUniformLocation(program, "u_top"), top);
            const GLint loc_level   = glGetUniformLocation(program, "u_level");
            const GLint loc_src_lod = glGetUniformLocation(program, "u_src_lod");
            const GLint loc_has_up  = glGetUniformLocation(program, "u_has_up");
            glActiveTexture(GL_TEXTURE3);
            glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.tex_tile[1]);
            for (int level = top; level >= 0; --level) {
                if (level > 0) {
                    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.bokeh_dof.tex_band_acc, level - 1);
                    const int w = gl.bokeh_dof.band_width  >> (level - 1);
                    const int h = gl.bokeh_dof.band_height >> (level - 1);
                    glViewport(0, 0, w, h);
                    glScissor(0, 0, w, h);
                } else {
                    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.bokeh_dof.tex_half[1], 0);
                    glViewport(0, 0, half_w, half_h);
                    glScissor(0, 0, half_w, half_h);
                }
                glActiveTexture(GL_TEXTURE0);
                glBindTexture(GL_TEXTURE_2D, level > 0 ? gl.bokeh_dof.tex_band_src : gl.bokeh_dof.tex_src[0]);
                glActiveTexture(GL_TEXTURE1);
                glBindTexture(GL_TEXTURE_2D, level > 0 ? gl.bokeh_dof.tex_band_coc : gl.bokeh_dof.tex_src[2]);
                glActiveTexture(GL_TEXTURE2);
                if (level < top) {
                    glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.tex_band_acc);
                    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, level);
                    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL,  level);
                } else {
                    glBindTexture(GL_TEXTURE_2D, 0);
                }
                glUniform1i(loc_level, level);
                glUniform1i(loc_src_lod, level > 0 ? level - 1 : 0);
                glUniform1i(loc_has_up, level < top ? 1 : 0);
                glDrawArrays(GL_TRIANGLES, 0, 3);
            }
            glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.tex_band_acc);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL,  top - 1);
            glBindTexture(GL_TEXTURE_2D, 0);
            glActiveTexture(GL_TEXTURE0);
        }
    }
    POP_GPU_SECTION()

    PUSH_GPU_SECTION("DOF fill")
    {
        // Pull-push fill of the far layer behind the foreground, see dof_fill.frag
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.tex_far_pull);
        glGenerateMipmap(GL_TEXTURE_2D);

        const GLuint program = gl.bokeh_dof.program_fill;
        glUseProgram(program);
        glUniform1i(glGetUniformLocation(program, "u_tex_pull"), 0);
        glUniform1i(glGetUniformLocation(program, "u_tex_coarse"), 1);
        const GLint loc_level = glGetUniformLocation(program, "u_level");
        const GLint loc_top   = glGetUniformLocation(program, "u_top");

        // Push level l is written to tex_far_push level l - 1 (to tex_half[0] for l = 0) and reads level l + 1, the one
        // above it. Restricting the sampled levels of tex_far_push to that one keeps it from being read where it is
        // written.
        glActiveTexture(GL_TEXTURE1);
        const int top = gl.bokeh_dof.pull_levels - 1;
        for (int level = top; level >= 0; --level) {
            const int w = MAX(gl.bokeh_dof.src_width  >> level, 1);
            const int h = MAX(gl.bokeh_dof.src_height >> level, 1);
            if (level > 0) {
                glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.bokeh_dof.tex_far_push, level - 1);
            } else {
                glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.bokeh_dof.tex_half[0], 0);
            }
            if (level < top) {
                glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.tex_far_push);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, level);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL,  level);
            } else {
                glBindTexture(GL_TEXTURE_2D, 0);   // nothing above the top level
            }
            glViewport(0, 0, level > 0 ? w : half_w, level > 0 ? h : half_h);
            glScissor (0, 0, level > 0 ? w : half_w, level > 0 ? h : half_h);
            glUniform1i(loc_level, level);
            glUniform1i(loc_top, level == top ? 1 : 0);
            glDrawArrays(GL_TRIANGLES, 0, 3);
        }
        glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.tex_far_push);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL,  MAX(top - 1, 0));
        glBindTexture(GL_TEXTURE_2D, 0);
        glActiveTexture(GL_TEXTURE0);
    }
    POP_GPU_SECTION()

    PUSH_GPU_SECTION("DOF composite")
    {
        reset_gl_state(reset_state);
        const GLuint program = gl.bokeh_dof.program_composite;
        glUseProgram(program);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, color_tex);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, linear_depth_tex);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.tex_half[0]);
        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_2D, gl.bokeh_dof.tex_half[1]);
        glUniform1i(glGetUniformLocation(program, "u_tex_color"), 0);
        glUniform1i(glGetUniformLocation(program, "u_tex_linear_depth"), 1);
        glUniform1i(glGetUniformLocation(program, "u_tex_far"), 2);
        glUniform1i(glGetUniformLocation(program, "u_tex_near"), 3);
        glUniform1f(glGetUniformLocation(program, "u_focus"), focus);
        glUniform1f(glGetUniformLocation(program, "u_aperture_px"), aperture_px);
        glUniform1f(glGetUniformLocation(program, "u_max_coc_px"), max_coc_px);
        glDrawArrays(GL_TRIANGLES, 0, 3);
    }
    POP_GPU_SECTION()

    glUseProgram(0);
    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(0);
}

// color_tex: premultiplied radiance of the scene. coverage_tex (r): coverage of the geometry. iso_tex (a, optional): the
// isosurfaces' coverage. background: display referred linear, put behind the scene after the curve.
void apply_display_transform(GLuint color_tex, GLuint coverage_tex, GLuint iso_tex, const display::Transform& t, vec3_t background) {
    ASSERT(glIsTexture(color_tex));
    ASSERT(glIsTexture(coverage_tex));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, color_tex);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, coverage_tex);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, iso_tex ? iso_tex : coverage_tex);     // a valid texture either way
    glActiveTexture(GL_TEXTURE0);

    const auto& p = tonemapping::display_transform;
    glUseProgram(p.program);
    glUniform1i(p.uniform_loc.texture, 0);
    glUniform1i(p.uniform_loc.texture_coverage, 1);
    glUniform1i(p.uniform_loc.texture_iso, 2);
    glUniform1i(p.uniform_loc.use_iso, iso_tex ? 1 : 0);
    glUniform3fv(p.uniform_loc.background, 1, background.elem);
    glUniform1i(p.uniform_loc.curve, (int)t.curve);
    glUniform1f(p.uniform_loc.exposure, t.exposure);
    glUniform1f(p.uniform_loc.gamma, t.gamma);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void apply_aa_tonemapping(GLuint color_tex) {
    ASSERT(glIsTexture(color_tex));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glUseProgram(tonemapping::fast_reversible.program_forward);
    glUniform1i(tonemapping::fast_reversible.uniform_loc.texture, 0);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void apply_inverse_aa_tonemapping(GLuint color_tex) {
    ASSERT(glIsTexture(color_tex));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glUseProgram(tonemapping::fast_reversible.program_inverse);
    glUniform1i(tonemapping::fast_reversible.uniform_loc.texture, 0);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void blit_static_velocity(GLuint depth_tex, const ViewParam& view_param) {
    //mat4_t curr_clip_to_prev_clip_mat = view_param.matrix.previous.view_proj * view_param.matrix.inverse.view_proj;
    mat4_t curr_clip_to_prev_clip_mat = view_param.matrix.prev.proj * view_param.matrix.prev.view * view_param.matrix.inv.view * view_param.matrix.inv.proj;

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, depth_tex);

    const vec2_t res = view_param.resolution;
    const vec2_t jitter_uv_cur = view_param.jitter.curr / res;
    const vec2_t jitter_uv_prev = view_param.jitter.prev / res;
    vec4_t jitter_uv = {jitter_uv_cur.x, jitter_uv_cur.y, jitter_uv_prev.x, jitter_uv_prev.y};

    glUseProgram(velocity::blit_velocity.program);
	glUniform1i(velocity::blit_velocity.uniform_loc.tex_depth, 0);
    glUniformMatrix4fv(velocity::blit_velocity.uniform_loc.curr_clip_to_prev_clip_mat, 1, GL_FALSE, &curr_clip_to_prev_clip_mat.elem[0][0]);
    glUniform4fv(velocity::blit_velocity.uniform_loc.jitter_uv, 1, jitter_uv.elem);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void blit_tilemax(GLuint velocity_tex, GLuint linear_depth_tex) {
    ASSERT(glIsTexture(velocity_tex));
    ASSERT(glIsTexture(linear_depth_tex));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, velocity_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);

    glUseProgram(velocity::blit_tilemax.program);
    glUniform1i(velocity::blit_tilemax.uniform_loc.tex_vel, 0);
    glUniform1i(velocity::blit_tilemax.uniform_loc.tex_linear_depth, 1);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void blit_neighbormax(GLuint velocity_tex, GLuint linear_depth_tex, int tex_width, int tex_height) {
    ASSERT(glIsTexture(velocity_tex));
    ASSERT(glIsTexture(linear_depth_tex));
    const vec2_t texel_size = {1.f / tex_width, 1.f / tex_height};

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, velocity_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);

    glUseProgram(velocity::blit_neighbormax.program);
    glUniform1i(velocity::blit_neighbormax.uniform_loc.tex_vel, 0);
    glUniform1i(velocity::blit_neighbormax.uniform_loc.tex_linear_depth, 1);
    glUniform2fv(velocity::blit_neighbormax.uniform_loc.tex_vel_texel_size, 1, &texel_size.x);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void blit_velocity_dilate(GLuint velocity_tex, int tex_width, int tex_height) {
    ASSERT(glIsTexture(velocity_tex));
    const vec2_t texel_size = {1.f / tex_width, 1.f / tex_height};

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, velocity_tex);

    glUseProgram(velocity::blit_dilate.program);
    glUniform1i(velocity::blit_dilate.uniform_loc.tex_vel, 0);
    glUniform2fv(velocity::blit_dilate.uniform_loc.texel_size, 1, &texel_size.x);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

// The temporal resolve into the bound targets: the history (location 0) and the image shown (location 1), of the
// viewport's size, out_width x out_height. The current frame is width x height, in the corner of the pipeline's textures
// (which can be larger, shared with larger views); upscaling when it is smaller than the output.
void apply_temporal_aa(GLuint linear_depth_tex, GLuint color_tex, GLuint prev_history_tex, GLuint velocity_tex, GLuint velocity_neighbormax_tex, const vec2_t& curr_jitter,
                       const vec2_t& prev_jitter, float feedback_min, float feedback_max, float motion_scale, float time,
                       vec2_t pixel_size, int width, int height, int out_width, int out_height) {
    const bool upscale = width != out_width || height != out_height;
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(color_tex));
    ASSERT(glIsTexture(velocity_tex));
    ASSERT(glIsTexture(velocity_neighbormax_tex));

    const vec2_t res = {(float)gl.tex_width, (float)gl.tex_height};
    const vec2_t inv_res = 1.0f / res;
    const vec4_t texel_size = vec4_t{inv_res.x, inv_res.y, res.x, res.y};
    const vec2_t in_size = {(float)width, (float)height};
    const vec4_t out_texel_size = vec4_t{1.0f / (float)out_width, 1.0f / (float)out_height, (float)out_width, (float)out_height};
    const vec2_t jitter_uv_curr = curr_jitter / res;
    const vec2_t jitter_uv_prev = prev_jitter / res;
    const vec4_t jitter_uv = vec4_t{jitter_uv_curr.x, jitter_uv_curr.y, jitter_uv_prev.x, jitter_uv_prev.y};

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, prev_history_tex ? prev_history_tex : color_tex);

    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, velocity_tex);

    glActiveTexture(GL_TEXTURE4);
    glBindTexture(GL_TEXTURE_2D, velocity_neighbormax_tex);

    const auto& p = gl.temporal[motion_scale != 0.f ? 1 : 0][upscale ? 1 : 0];
    glUseProgram(p.program);
    glUniform1i(p.uniform_loc.tex_linear_depth, 0);
    glUniform1i(p.uniform_loc.tex_main, 1);
    glUniform1i(p.uniform_loc.tex_prev, 2);
    glUniform1i(p.uniform_loc.tex_vel, 3);
    glUniform1i(p.uniform_loc.tex_vel_neighbormax, 4);
    glUniform4fv(p.uniform_loc.texel_size, 1, &texel_size.x);
    glUniform2fv(p.uniform_loc.in_size, 1, in_size.elem);
    glUniform4fv(p.uniform_loc.out_texel_size, 1, &out_texel_size.x);
    glUniform4fv(p.uniform_loc.jitter_uv, 1, &jitter_uv.x);
    glUniform1f(p.uniform_loc.time, time);
    glUniform1f(p.uniform_loc.feedback_min, feedback_min);
    glUniform1f(p.uniform_loc.feedback_max, feedback_max);
    glUniform1f(p.uniform_loc.motion_scale, motion_scale);
    glUniform2fv(p.uniform_loc.pixel_size, 1, pixel_size.elem);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

// The source (src_width x src_height in the corner of a texture of the pipeline's size) to the viewport, Catmull-Rom
static void upscale_texture(GLuint src, int src_width, int src_height) {
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, src);
    glUseProgram(gl.upscale.program);
    glUniform1i(gl.upscale.loc_tex, 0);
    glUniform4f(gl.upscale.loc_src_texel, 1.0f / (float)gl.tex_width, 1.0f / (float)gl.tex_height, (float)gl.tex_width, (float)gl.tex_height);
    glUniform2f(gl.upscale.loc_src_size, (float)src_width, (float)src_height);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void blit_texture(GLuint tex) {
    ASSERT(glIsTexture(tex));
    glUseProgram(blit::program_tex);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);
    glUniform1i(blit::uniform_loc_texture, 0);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

static void blit_texture_dither(GLuint tex, float time) {
    ASSERT(glIsTexture(tex));
    glUseProgram(blit::program_tex_dither);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);
    glUniform1i(glGetUniformLocation(blit::program_tex_dither, "u_texture"), 0);
    glUniform1f(glGetUniformLocation(blit::program_tex_dither, "u_time"), time);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

// The final copy, with the marks over it (Settings::marks)
static void present_marks(GLuint tex, GLuint marks, GLuint marks_group, const postprocess_pipeline::Settings& settings, float time) {
    ASSERT(glIsTexture(tex));
    ASSERT(glIsTexture(marks));
    ASSERT(glIsTexture(marks_group));
    const GLuint p = blit::program_marks;
    glUseProgram(p);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, marks);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, marks_group);
    glUniform1i(glGetUniformLocation(p, "u_texture"), 0);
    glUniform1i(glGetUniformLocation(p, "u_marks"), 1);
    glUniform1i(glGetUniformLocation(p, "u_marks_group"), 2);
    glUniform3fv(glGetUniformLocation(p, "u_groups"), 7, settings.marks.groups[0].elem);
    glUniform1f(glGetUniformLocation(p, "u_time"), time);
    glUniform4fv(glGetUniformLocation(p, "u_visible"), 2, settings.marks.visible[0].elem);
    glUniform4fv(glGetUniformLocation(p, "u_hidden"),  2, settings.marks.hidden[0].elem);
    glUniform4fv(glGetUniformLocation(p, "u_outline"), 2, settings.marks.outline[0].elem);
    glUniform1f(glGetUniformLocation(p, "u_highlight"), settings.marks.highlight);
    glUniform1i(glGetUniformLocation(p, "u_outline_step"), MAX(settings.marks.outline_step, 1));
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glActiveTexture(GL_TEXTURE0);
    glUseProgram(0);
}

// The coverage of the geometry in the G-buffer colour gbuffer_color, into r
static void blit_coverage(GLuint gbuffer_color) {
    ASSERT(glIsTexture(gbuffer_color));
    glUseProgram(blit::program_coverage);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, gbuffer_color);
    glUniform1i(glGetUniformLocation(blit::program_coverage, "u_texture"), 0);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void blur_texture_gaussian(GLuint tex, int num_passes) {
    ASSERT(glIsTexture(tex));
    ASSERT(num_passes > 0);

    GLResetState reset_state = {};
    record_gl_reset_state(&reset_state);

    GLint w, h;

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);

    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

    glBindVertexArray(gl.vao);

    glUseProgram(blur::program_gaussian);
    glUniform1i(blur::uniform_loc_texture, 0);

    glViewport(0, 0, w, h);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.rt.scratch_fbo);

    for (int i = 0; i < num_passes; ++i) {
        glBindTexture(GL_TEXTURE_2D, tex);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.rt.tex_rgba8, 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        glUniform2f(blur::uniform_loc_inv_res_dir, 1.0f / w, 0.0f);
        glDrawArrays(GL_TRIANGLES, 0, 3);

        glBindTexture(GL_TEXTURE_2D, gl.rt.tex_rgba8);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex, 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        glUniform2f(blur::uniform_loc_inv_res_dir, 0.0f, 1.0f / h);
        glDrawArrays(GL_TRIANGLES, 0, 3);
    }

    glUseProgram(0);
    glBindVertexArray(0);

    reset_gl_state(reset_state);
}

void blur_texture_box(GLuint tex, int num_passes) {
    ASSERT(glIsTexture(tex));
    ASSERT(num_passes > 0);

    GLResetState reset_state = {};
    record_gl_reset_state(&reset_state);

    GLint w, h;

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);

    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

    glBindVertexArray(gl.vao);

    glUseProgram(blur::program_box);

    glViewport(0, 0, w, h);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.rt.scratch_fbo);

    for (int i = 0; i < num_passes; ++i) {
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.rt.tex_rgba8, 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        glDrawArrays(GL_TRIANGLES, 0, 3);

        glBindTexture(GL_TEXTURE_2D, gl.rt.tex_rgba8);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex, 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        glDrawArrays(GL_TRIANGLES, 0, 3);
    }

    glUseProgram(0);
    glBindVertexArray(0);

    reset_gl_state(reset_state);
}

static void compute_luma(GLuint tex) {
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);

    glBindVertexArray(gl.vao);
    glUseProgram(gl.luma.program);
    glUniform1i(gl.luma.uniform_loc.tex_rgba, 0);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glUseProgram(0);
    glBindVertexArray(0);
}

static void compute_fxaa(GLuint tex_rgbl, int width, int height) {
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex_rgbl);

    int w, h;
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

    float scl_x = (float)width  / (float)w;
    float scl_y = (float)height / (float)h;

    glBindVertexArray(gl.vao);
    glUseProgram(gl.fxaa.program);
    glUniform2f(gl.fxaa.uniform_loc.tc_scl, scl_x, scl_y);
    glUniform1i(gl.fxaa.uniform_loc.tex_rgbl, 0);
    glUniform2f(gl.fxaa.uniform_loc.rcp_res, 1.0f / (float)w, 1.0f / (float)h);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glUseProgram(0);
    glBindVertexArray(0);
}

void execute(const postprocess_pipeline::Inputs& in, const postprocess_pipeline::Settings& settings, const ViewParam& view_param) {
    ASSERT(glIsTexture(in.depth));
    ASSERT(glIsTexture(in.color));
    ASSERT(glIsTexture(in.normal));
    if (settings.taa.enabled) {
        ASSERT(glIsTexture(in.velocity));
    }

    const bool do_compose = true;
    const bool do_velocity = settings.taa.enabled;
    const bool do_ssao = settings.ssao.enabled;
    const bool do_shadows = settings.shadows.enabled && settings.shading.light > 0.0f;
    const bool do_dof = settings.dof.enabled;
    const bool do_transparency = in.transparency != 0;
    const bool do_transparency_hdr = in.transparency_hdr != 0;
    const bool do_fxaa = settings.fxaa.enabled;
    const bool do_taa = settings.taa.enabled && do_velocity && in.history;
    const bool do_sharpen = settings.sharpen.enabled;
    const bool do_present = true;

    // For seeding noise
    static float time = 0.f;
    time = time + 0.01f;
    if (time > 100.f) time -= 100.f;

    GLResetState reset_state = {};
    record_gl_reset_state(&reset_state);

    // The image shown is the viewport's size; the G-buffer is that or, upscaling, smaller (in.width / in.height)
    const int out_width  = reset_state.viewport[2];
    const int out_height = reset_state.viewport[3];
    int width  = in.width  > 0 ? in.width  : out_width;
    int height = in.height > 0 ? in.height : out_height;
    const bool upscale = width != out_width || height != out_height;

    // The targets are shared by the views (the main one and the component windows'), of different sizes: here they only
    // grow, and each view uses their corner. The main view sets them to its size when it changes (initialize()).
    if (width > (int)gl.tex_width || height > (int)gl.tex_height || out_width > (int)gl.out_width || out_height > (int)gl.out_height) {
        initialize(MAX(width, (int)gl.tex_width), MAX(height, (int)gl.tex_height), MAX(out_width, (int)gl.out_width), MAX(out_height, (int)gl.out_height));
    }

    // Postprocessing targets have no depth attachment – disable depth test so no
    // driver discards fragments silently when the G-buffer pass leaves it enabled.
    GLboolean last_depth_test = glIsEnabled(GL_DEPTH_TEST);
    glDisable(GL_DEPTH_TEST);

    glEnable(GL_SCISSOR_TEST);
    glScissor (0, 0, width, height);
    glViewport(0, 0, width, height);
    glBindVertexArray(gl.vao);

    const auto near_dist = view_param.clip_planes.near;
    const auto far_dist = view_param.clip_planes.far;
    const auto ortho = is_orthographic_proj_matrix(view_param.matrix.curr.proj.elem);

    PUSH_GPU_SECTION("Linearize Depth")
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.linear_depth.fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.linear_depth.texture, 0);
    compute_linear_depth(in.depth, near_dist, far_dist, ortho);
    POP_GPU_SECTION()

    const GLenum draw_buffers[2] = {
        GL_COLOR_ATTACHMENT0,
        GL_COLOR_ATTACHMENT1,
    };

    if (do_ssao || do_shadows) {
        PUSH_GPU_SECTION("Generate Linear Depth Mipmaps") {
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.linear_depth.fbo);
            glDrawBuffer(GL_COLOR_ATTACHMENT0);
            for (int level = 1; level < gl.linear_depth.levels; ++level) {
                const int w = MAX((int)gl.tex_width  >> level, 1);
                const int h = MAX((int)gl.tex_height >> level, 1);
                glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.linear_depth.texture, level);
                glViewport(0, 0, w, h);
                glScissor(0, 0, w, h);
                downsample_depth(gl.linear_depth.texture, level - 1);
            }
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.linear_depth.texture, 0);
        }
        POP_GPU_SECTION()
    }

    if (do_velocity) {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.velocity.fbo);
        glViewport(0, 0, gl.velocity.tex_width, gl.velocity.tex_height);
        glScissor(0, 0, gl.velocity.tex_width, gl.velocity.tex_height);

        PUSH_GPU_SECTION("Velocity: Tilemax") {
            glDrawBuffer(GL_COLOR_ATTACHMENT0);
            blit_tilemax(in.velocity, gl.linear_depth.texture);
        }
        POP_GPU_SECTION()

        PUSH_GPU_SECTION("Velocity: Neighbormax") {
            glDrawBuffer(GL_COLOR_ATTACHMENT1);
            blit_neighbormax(gl.velocity.tex_tilemax, gl.linear_depth.texture, gl.velocity.tex_width, gl.velocity.tex_height);
        }
        POP_GPU_SECTION()
    }
     
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.rt.fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.rt.tex_color[0], 0);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gl.rt.tex_color[1], 0);
    // The history is the output's size: upscaling, it goes to the targets of that size (to_output_size below)
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, upscale ? 0 : in.history, 0);
    glViewport(0, 0, width, height);
    glScissor(0, 0, width, height);
    glDrawBuffers(2, draw_buffers);
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);

    // Ping-pong between attachment 0 and 1 of gl.rt.fbo. Which textures sit there changes at the tone mapping pass
    // (HDR pair -> LDR pair), so the source texture is tracked explicitly.
    const GLuint* pair = gl.rt.tex_color;
    GLenum dst_buffer = GL_COLOR_ATTACHMENT0;
    GLuint src_texture = 0;

    auto swap_target = [&]() {
        src_texture = pair[dst_buffer - GL_COLOR_ATTACHMENT0];
        dst_buffer = dst_buffer == GL_COLOR_ATTACHMENT0 ? GL_COLOR_ATTACHMENT1 : GL_COLOR_ATTACHMENT0;
    };
    glDrawBuffer(dst_buffer);

    // The visibilities first, the ambient occlusion's and the shadows': the compose pass applies them to the light they
    // occlude, the surroundings' and the key light's
    if (do_ssao) {
        PUSH_GPU_SECTION("SSAO")
        const bool temporal = settings.ssao.temporal && in.velocity != 0;
        compute_ssao(gl.linear_depth.texture, in.normal, in.velocity, view_param, far_dist, settings.ssao.intensity, temporal, settings.ssao.mode);
        POP_GPU_SECTION()
    }

    if (do_shadows) {
        PUSH_GPU_SECTION("Shadows")
        const postprocess_pipeline::Lighting light = postprocess_pipeline::lighting(settings.shading.mode);
        const bool temporal = settings.shadows.temporal && in.velocity != 0;
        compute_shadows(gl.linear_depth.texture, in.normal, in.velocity, view_param, far_dist, light.light_dir, (float)DEG_TO_RAD(settings.shadows.softness), temporal);
        POP_GPU_SECTION()
    }

    if (do_compose) {
        PUSH_GPU_SECTION("Compose")
        compose_deferred(gl.linear_depth.texture, in.color, in.normal, do_ssao ? gl.ssao.ao.tex_visibility : 0, do_shadows ? gl.ssao.shadow.tex_visibility : 0,
                         view_param.matrix.curr.proj.elem, ortho, width, height, settings, time);
        POP_GPU_SECTION()
    } else {
        // Seed transient pipeline from external color input when compose stage is disabled.
        blit_texture(in.color);
    }

    // HDR effects -----------------------------------------------------------
    if (do_dof) {
        swap_target();
        glDrawBuffer(dst_buffer);
        PUSH_GPU_SECTION("DOF")
        apply_dof(gl.linear_depth.texture, src_texture, settings.dof.focus_depth, settings.dof.aperture);
        POP_GPU_SECTION()
    }

    // After DOF: the isosurfaces have no depth of their own for it to work with, so they stay sharp
    if (do_transparency_hdr) {
        PUSH_GPU_SECTION("Add HDR Transparency")
        glDrawBuffer(dst_buffer);
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
        blit_texture(in.transparency_hdr);
        glDisable(GL_BLEND);
        POP_GPU_SECTION()
    }

    // HDR -> LDR boundary ---------------------------------------------------

    // The coverage of the geometry, which puts the scene over the background: sharp from the G-buffer, and blurred as
    // the colour was when the depth of field is on
    GLuint coverage_tex = gl.rt.tex_coverage[0];
    {
        PUSH_GPU_SECTION("Coverage")
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.rt.scratch_fbo);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.rt.tex_coverage[0], 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        blit_coverage(in.color);
        if (do_dof) {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.rt.tex_coverage[1], 0);
            apply_dof(gl.linear_depth.texture, gl.rt.tex_coverage[0], settings.dof.focus_depth, settings.dof.aperture);
            coverage_tex = gl.rt.tex_coverage[1];
        }
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.rt.fbo);
        POP_GPU_SECTION()
    }

    // Always: tone mapping off (and basic shading) is the display transform without a curve, not no transform
    {
        PUSH_GPU_SECTION("Display transform")
        swap_target();
        // From here on the values are display referred: switch the ping-pong pair to the RGB10_A2 targets
        pair = gl.rt.tex_ldr;
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.rt.tex_ldr[0], 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gl.rt.tex_ldr[1], 0);
        glDrawBuffer(dst_buffer);
        const display::Transform t = display::transform(settings);
        apply_display_transform(src_texture, coverage_tex, do_transparency_hdr ? in.transparency_hdr : 0, t, display::background(settings.background_color, t));
        POP_GPU_SECTION()
    }

    // LDR effects -----------------------------------------------------------

    if (do_transparency) {
        PUSH_GPU_SECTION("Add Transparency")
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);    // premultiplied
        blit_texture(in.transparency);
        glDisable(GL_BLEND);
        POP_GPU_SECTION()
    }

    if (do_fxaa) {
        swap_target();
        PUSH_GPU_SECTION("Luma")
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.rt.scratch_fbo);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.rt.tex_rgba8, 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        compute_luma(src_texture);
        POP_GPU_SECTION()
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.rt.fbo);
        glDrawBuffer(dst_buffer);
        PUSH_GPU_SECTION("FXAA")
        compute_fxaa(gl.rt.tex_rgba8, width, height);
        POP_GPU_SECTION()
    }

    // From here on, upscaling: the targets at the output size. The ping-pong goes on in them.
    auto to_output_size = [&]() {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.rt.out_fbo);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.rt.tex_out[0], 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gl.rt.tex_out[1], 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, in.history, 0);
        glViewport(0, 0, out_width, out_height);
        glScissor(0, 0, out_width, out_height);
        pair = gl.rt.tex_out;
        dst_buffer = GL_COLOR_ATTACHMENT0;
    };

    if (do_taa) {
        GLuint history_prev = in.history_prev;
        if (!history_prev) {
            copy_texture_2d(gl.rt.tex_history_prev, in.history, out_width, out_height);
            history_prev = gl.rt.tex_history_prev;
        }

        swap_target();
        if (upscale) to_output_size();
        const float feedback_min = settings.taa.feedback_min;
        const float feedback_max = settings.taa.feedback_max;
        const float motion_scale = settings.taa.motion_blur.enabled ? settings.taa.motion_blur.motion_scale : 0.f;

        // The temporal shader writes two outputs:
        //   location 0 (out_buff) = temporal resolve, used as history next frame
        //   location 1 (out_frag) = motion-blur composited result, used for display
        // When motion blur is active we must enable both draw buffers so location 1
        // reaches tex_taa_frag; otherwise a single buffer suffices (out_buff == out_frag).
        PUSH_GPU_SECTION("Temporal AA")
        const GLenum taa_buffers[2] = { GL_COLOR_ATTACHMENT2, dst_buffer };
        glDrawBuffers(2, taa_buffers);

        // A pixel's height in view space: x at linear depth 1 (perspective, y = 1) or anywhere (orthographic, y = 0).
        // The prefilter measures steps in depth by it.
        const float (*P)[4] = view_param.matrix.curr.proj.elem;
        const vec2_t pixel_size = {2.0f / (P[1][1] * (float)height), is_orthographic_proj_matrix(P) ? 0.0f : 1.0f};

        apply_temporal_aa(gl.linear_depth.texture, src_texture, history_prev, in.velocity, gl.velocity.tex_neighbormax, view_param.jitter.curr, view_param.jitter.prev,
                          feedback_min, feedback_max, motion_scale, time, pixel_size, width, height, out_width, out_height);

        POP_GPU_SECTION()
    } else if (upscale) {
        // Each frame on its own (FXAA or none)
        PUSH_GPU_SECTION("Upscale")
        swap_target();
        to_output_size();
        glDrawBuffer(dst_buffer);
        upscale_texture(src_texture, width, height);
        POP_GPU_SECTION()
    }

    if (do_sharpen) {
        PUSH_GPU_SECTION("Sharpen")
        swap_target();
        glDrawBuffer(dst_buffer);
        sharpen::sharpen(src_texture, settings.sharpen.weight);
        POP_GPU_SECTION()
    }

    // Activate backbuffer or whatever was bound before
    reset_gl_state(reset_state);
    glDisable(GL_SCISSOR_TEST);

    if (do_present) {
        swap_target();
        glDepthMask(0);
        if (!settings.present) {
            // Kept to be presented later (postprocess_pipeline::present): as it is, without the dither
            blit_texture(src_texture);
        } else if (in.marks && in.marks_group) {
            PUSH_GPU_SECTION("Present marks")
            present_marks(src_texture, in.marks, in.marks_group, settings, time);
            POP_GPU_SECTION()
        } else {
            blit_texture_dither(src_texture, time);
        }
    }

    if (last_depth_test) glEnable(GL_DEPTH_TEST);
    glDepthMask(1);
    glColorMask(1, 1, 1, 1);
}

}  // namespace postprocessing

namespace postprocess_pipeline {

void initialize(int width, int height, int out_width, int out_height) {
    // glTexStorage2D rejects zero dimensions outright (GL_INVALID_VALUE), so guard the
    // 0x0 startup/iconified case here as well. See gbuffer_init().
    postprocessing::initialize(MAX(width, 1), MAX(height, 1), MAX(out_width, 0), MAX(out_height, 0));
}

void shutdown() {
    postprocessing::shutdown();
}


Lighting lighting(ShadingMode mode) {
    // Physically based: the key light from the upper right, 55 degrees off the view axis, which models the shapes,
    // and light from all around, nearly uniform: a sky a little brighter than the ground, which keeps a hint of up
    // and down. The balance is that of the earlier shading (a uniform environment of radiance 6 against a key light
    // of 10, in its units): the side away from the key light keeps about 0.37 of the lit side. A white surface
    // facing the viewer:
    //   key: 0.87 * cos(55 deg) = 0.50, environment: (0.59 + 0.49) / 2 = 0.54
    // which, less what the specular reflection takes off the diffuse (Fresnel, about 3.5%), is a diffuse reflection
    // of 1: reference white. The specular reflection adds about 0.02 (Neutral's offset takes off 0.04). Facing the
    // key light it is about 1.4. Settings::shading.light and .environment scale the two.
    // Basic: from the upper right too, but 27 degrees off the view axis, so the surfaces facing the viewer keep
    // close to their full colour. The radiances are not used.
    Lighting l = {
        .light_dir        = vec3_normalize(vec3_set(1.0f, 1.0f, 1.0f)),
        .light_irradiance = vec3_set1(0.87f * (float)PI),
        .sky_radiance     = vec3_set1(0.59f),
        .ground_radiance  = vec3_set1(0.49f),
    };
    if (mode == ShadingMode_Basic) {
        l.light_dir = vec3_normalize(vec3_set(0.3f, 0.4f, 1.0f));
    }
    return l;
}

void execute(const Inputs& in, const Settings& settings, const ViewParam& view) {
    ASSERT(glIsTexture(in.depth));
    ASSERT(glIsTexture(in.color));
    ASSERT(glIsTexture(in.normal));
    if (settings.taa.enabled) {
        ASSERT(glIsTexture(in.velocity));
    }

    postprocessing::execute(in, settings, view);
}

void present(GLuint image, GLuint marks, GLuint marks_group, const Settings& settings) {
    ASSERT(glIsTexture(image));
    // For the dither, as execute() seeds it
    static float time = 0.f;
    time = time + 0.01f;
    if (time > 100.f) time -= 100.f;

    const GLboolean last_depth_test = glIsEnabled(GL_DEPTH_TEST);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glDepthMask(0);
    if (marks && marks_group) {
        PUSH_GPU_SECTION("Present marks")
        postprocessing::present_marks(image, marks, marks_group, settings, time);
        POP_GPU_SECTION()
    } else {
        postprocessing::blit_texture_dither(image, time);
    }
    if (last_depth_test) glEnable(GL_DEPTH_TEST);
    glDepthMask(1);
}

}  // namespace postprocess_pipeline
