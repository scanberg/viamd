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

#include <shaders.inl>

#define PUSH_GPU_SECTION(lbl)                                                                       \
    {                                                                                               \
        if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); \
    }
#define POP_GPU_SECTION()                       \
    {                                           \
        if (glPopDebugGroup) glPopDebugGroup(); \
    }

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

typedef int Tonemapping;
enum Tonemapping_ {
    Tonemapping_Passthrough,
    Tonemapping_ExposureGamma,
    Tonemapping_Filmic,
    Tonemapping_ACES,
};

static Tonemapping to_legacy_tonemapper(postprocess_pipeline::Tonemapper tm) {
    switch (tm) {
        case postprocess_pipeline::Tonemapper_ExposureGamma: return Tonemapping_ExposureGamma;
        case postprocess_pipeline::Tonemapper_Filmic: return Tonemapping_Filmic;
        case postprocess_pipeline::Tonemapper_ACES: return Tonemapping_ACES;
        case postprocess_pipeline::Tonemapper_Passthrough:
        default: return Tonemapping_Passthrough;
    }
}

static inline bool is_orthographic_proj_matrix(const float P[4][4]) { return P[2][3] == 0.0f; }

static struct {
    GLuint vao = 0;
    GLuint v_shader_fs_quad = 0;
    bool programs_ready = false;    // programs are compiled once; initialize() on resize only reallocates targets
    uint32_t tex_width = 0;
    uint32_t tex_height = 0;

    struct {
        GLuint fbo = 0;
        GLuint scratch_fbo = 0;
        GLuint tex_rgba8 = 0;
        GLuint tex_half_dof = 0;        // half res RGBA16F: colour + signed CoC (DOF prepass)
        GLuint tex_color[2] = {0, 0};   // HDR ping-pong (R11F_G11F_B10F): compose, SSAO, DOF
        GLuint tex_ldr[2] = {0, 0};     // LDR ping-pong (RGB10_A2): everything after tone mapping
        GLuint tex_history_prev = 0;
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

    struct {
        GLuint fbo = 0;
        GLuint tex[2] = {};             // half res RG16F (visibility, linear depth), ping-pong for the separable blur
        // Index 0 = orthographic, 1 = perspective
        GLuint program_ao[2] = {};
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

    struct {
        struct {
            GLuint program = 0;
            struct {
                GLint tex_linear_depth = -1;
                GLint tex_main = -1;
                GLint tex_prev = -1;
                GLint tex_vel = -1;
                GLint tex_vel_neighbormax = -1;
                GLint texel_size = -1;
                GLint time = -1;
                GLint feedback_min = -1;
                GLint feedback_max = -1;
                GLint motion_scale = -1;
                GLint jitter_uv = -1;
            } uniform_loc;
        } with_motion_blur;
        struct {
            GLuint program = 0;
            struct {
                GLint tex_linear_depth = -1;
                GLint tex_main = -1;
                GLint tex_prev = -1;
                GLint tex_vel = -1;
                GLint texel_size = -1;
                GLint time = -1;
                GLint feedback_min = -1;
                GLint feedback_max = -1;
                GLint motion_scale = -1;
                GLint jitter_uv = -1;
            } uniform_loc;
        } no_motion_blur;
    } temporal;
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
static constexpr int   MAX_MIP        = 5;      // must match AO_MAX_MIP in ssao.frag
static constexpr float R_MIN_PX       = 2.0f;   // smallest sampling distance, full-res pixels
static constexpr float R_MAX_FRACTION = 0.1f;   // largest sampling distance, fraction of the viewport height

static GLuint setup_variant(str_t name, str_t src, bool perspective) {
    return setup_program_from_source(name, src, perspective ? STR_LIT("#define AO_PERSPECTIVE 1") : STR_LIT("#define AO_PERSPECTIVE 0"));
}

void initialize_programs() {
    const str_t ao_src       = {(const char*)ssao_frag, ssao_frag_size};
    const str_t blur_src     = {(const char*)blur_frag, blur_frag_size};
    const str_t upsample_src = {(const char*)upsample_frag, upsample_frag_size};
    for (int i = 0; i < 2; ++i) {
        if (gl.ssao.program_ao[i])       glDeleteProgram(gl.ssao.program_ao[i]);
        if (gl.ssao.program_blur[i])     glDeleteProgram(gl.ssao.program_blur[i]);
        if (gl.ssao.program_upsample[i]) glDeleteProgram(gl.ssao.program_upsample[i]);
        gl.ssao.program_ao[i]       = setup_variant(STR_LIT("ssao"),          ao_src,       i == 1);
        gl.ssao.program_blur[i]     = setup_variant(STR_LIT("ssao blur"),     blur_src,     i == 1);
        gl.ssao.program_upsample[i] = setup_variant(STR_LIT("ssao upsample"), upsample_src, i == 1);
    }
}

void initialize_targets(int width, int height) {
    // Must have the dimensions of mip level 1 of the linear depth texture
    const int half_w = MAX(width / 2, 1);
    const int half_h = MAX(height / 2, 1);
    ensure_texture_2d(&gl.ssao.tex[0], GL_RG16F, half_w, half_h, GL_RG, GL_FLOAT, GL_NEAREST, GL_NEAREST);
    ensure_texture_2d(&gl.ssao.tex[1], GL_RG16F, half_w, half_h, GL_RG, GL_FLOAT, GL_NEAREST, GL_NEAREST);
}

void shutdown() {
    if (gl.ssao.tex[0]) glDeleteTextures(2, gl.ssao.tex);
    gl.ssao.tex[0] = gl.ssao.tex[1] = 0;
    for (int i = 0; i < 2; ++i) {
        if (gl.ssao.program_ao[i])       glDeleteProgram(gl.ssao.program_ao[i]);
        if (gl.ssao.program_blur[i])     glDeleteProgram(gl.ssao.program_blur[i]);
        if (gl.ssao.program_upsample[i]) glDeleteProgram(gl.ssao.program_upsample[i]);
        gl.ssao.program_ao[i] = gl.ssao.program_blur[i] = gl.ssao.program_upsample[i] = 0;
    }
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

namespace highlight {

static struct {
    GLuint program = 0;
    GLuint selection_texture = 0;
    struct {
        GLint texture_atom_idx = -1;
        GLint buffer_selection = -1;
        GLint highlight = -1;
        GLint selection = -1;
        GLint outline = -1;
    } uniform_loc;
} highlight;

void initialize() {
    highlight.program = setup_program_from_source(STR_LIT("highlight"), {(const char*)highlight_frag, highlight_frag_size});
    if (!highlight.selection_texture) glGenTextures(1, &highlight.selection_texture);
    highlight.uniform_loc.texture_atom_idx = glGetUniformLocation(highlight.program, "u_texture_atom_idx");
    highlight.uniform_loc.buffer_selection = glGetUniformLocation(highlight.program, "u_buffer_selection");
    highlight.uniform_loc.highlight = glGetUniformLocation(highlight.program, "u_highlight");
    highlight.uniform_loc.selection = glGetUniformLocation(highlight.program, "u_selection");
    highlight.uniform_loc.outline = glGetUniformLocation(highlight.program, "u_outline");
}

void shutdown() {
    if (highlight.program) glDeleteProgram(highlight.program);
}
}  // namespace highlight

namespace hsv {

static struct {
    GLuint program = 0;
    struct {
        GLint texture_color = -1;
        GLint hsv_scale = -1;
    } uniform_loc;
} gl;

void initialize() {
    gl.program = setup_program_from_source(STR_LIT("scale hsv"), {(const char*)scale_hsv_frag, scale_hsv_frag_size});
    gl.uniform_loc.texture_color = glGetUniformLocation(gl.program, "u_texture_atom_color");
    gl.uniform_loc.hsv_scale = glGetUniformLocation(gl.program, "u_hsv_scale");
}

void shutdown() {
    if (gl.program) glDeleteProgram(gl.program);
}
}  // namespace hsv

namespace compose {

struct ubo_data_t {
    vec4_t proj_info;
    vec3_t bg_color;
    float  time;
    vec3_t env_radiance;
    float  roughness;
    vec3_t dir_radiance;
    float  F0;
    vec3_t light_dir;
};

static struct {
    struct {
        GLuint program = 0;
        struct {
            GLint uniform_data = -1;
            GLint texture_depth = -1;
            GLint texture_color = -1;
            GLint texture_normal = -1;
        } uniform_loc;
    } persp, ortho;
    GLuint ubo = 0;
} compose;

static void init_variant(decltype(compose.persp)& v, str_t defines) {
    v.program = setup_program_from_source(STR_LIT("compose deferred"), {(const char*)compose_deferred_frag, compose_deferred_frag_size}, defines);
    v.uniform_loc.texture_depth = glGetUniformLocation(v.program, "u_texture_depth");
    v.uniform_loc.texture_color = glGetUniformLocation(v.program, "u_texture_color");
    v.uniform_loc.texture_normal = glGetUniformLocation(v.program, "u_texture_normal");
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
    } uniform_loc;
} passthrough;

static struct {
    GLuint program = 0;
    struct {
        GLint texture = -1;
        GLint exposure = -1;
        GLint gamma = -1;
    } uniform_loc;
} exposure_gamma;

static struct {
    GLuint program = 0;
    struct {
        GLint texture = -1;
        GLint exposure = -1;
        GLint gamma = -1;
    } uniform_loc;
} filmic;

static struct {
    GLuint program = 0;
    struct {
        GLint texture = -1;
        GLint exposure = -1;
        GLint gamma = -1;
    } uniform_loc;
} aces;

static struct {
    GLuint program_forward = 0;
    GLuint program_inverse = 0;
    struct {
        GLint texture = -1;
    } uniform_loc;
} fast_reversible;

void initialize() {
    {
        // PASSTHROUGH
        passthrough.program = setup_program_from_source(STR_LIT("Passthrough"), {(const char*)passthrough_frag, passthrough_frag_size});
        passthrough.uniform_loc.texture = glGetUniformLocation(passthrough.program, "u_texture");
    }
    {
        // EXPOSURE GAMMA
        exposure_gamma.program = setup_program_from_source(STR_LIT("Exposure Gamma"), {(const char*)exposure_gamma_frag, exposure_gamma_frag_size});
        exposure_gamma.uniform_loc.texture = glGetUniformLocation(exposure_gamma.program, "u_texture");
        exposure_gamma.uniform_loc.exposure = glGetUniformLocation(exposure_gamma.program, "u_exposure");
        exposure_gamma.uniform_loc.gamma = glGetUniformLocation(exposure_gamma.program, "u_gamma");
    }
    {
        // UNCHARTED
        filmic.program = setup_program_from_source(STR_LIT("Filmic"), {(const char*)uncharted_frag, uncharted_frag_size});
        filmic.uniform_loc.texture = glGetUniformLocation(filmic.program, "u_texture");
        filmic.uniform_loc.exposure = glGetUniformLocation(filmic.program, "u_exposure");
        filmic.uniform_loc.gamma = glGetUniformLocation(filmic.program, "u_gamma");
    }
    {
        // ACES
        aces.program = setup_program_from_source(STR_LIT("ACES"), {(const char*)aces_frag, aces_frag_size});
        aces.uniform_loc.texture = glGetUniformLocation(aces.program, "u_texture");
        aces.uniform_loc.exposure = glGetUniformLocation(aces.program, "u_exposure");
        aces.uniform_loc.gamma = glGetUniformLocation(aces.program, "u_gamma");
    }
    {
        // Fast Reversible (For AA) (Credits to Brian Karis: http://graphicrants.blogspot.com/2013/12/tone-mapping.html)
        fast_reversible.program_forward = setup_program_from_source(STR_LIT("Fast Reversible"), {(const char*)fast_reversible_frag, fast_reversible_frag_size}, STR_LIT("#define USE_INVERSE 0"));
        fast_reversible.program_inverse = setup_program_from_source(STR_LIT("Fast Reversible"), {(const char*)fast_reversible_frag, fast_reversible_frag_size}, STR_LIT("#define USE_INVERSE 1"));
        fast_reversible.uniform_loc.texture = glGetUniformLocation(fast_reversible.program_forward, "u_texture");
    }
}

void shutdown() {
    if (passthrough.program) glDeleteProgram(passthrough.program);
    if (exposure_gamma.program) glDeleteProgram(exposure_gamma.program);
    if (filmic.program) glDeleteProgram(filmic.program);
    if (aces.program) glDeleteProgram(aces.program);
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
static GLuint program_tex_covered = 0;
static GLuint program_col = 0;
static GLint uniform_loc_texture = -1;
static GLint uniform_loc_color = -1;

// A premultiplied texture, covered by the coverage (alpha) of another: what lies under it
constexpr str_t f_shader_src_tex_covered = STR_LIT(R"(
#version 150 core

uniform sampler2D u_texture;
uniform sampler2D u_cover;

out vec4 out_frag;

void main() {
    ivec2 px = ivec2(gl_FragCoord.xy);
    out_frag = texelFetch(u_texture, px, 0) * (1.0 - clamp(texelFetch(u_cover, px, 0).a, 0.0, 1.0));
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

// Premultiplied: what it writes goes into the transparency buffer, which holds premultiplied colour
constexpr str_t f_shader_src_col = STR_LIT(R"(
#version 150 core

uniform vec4 u_color;
out vec4 out_frag;

void main() {
	float a = clamp(u_color.a, 0.0, 1.0);
	out_frag = vec4(u_color.rgb * a, a);
}
)");

void initialize() {
    program_tex = setup_program_from_source(STR_LIT("blit texture"), f_shader_src_tex);
    uniform_loc_texture = glGetUniformLocation(program_tex, "u_texture");

    program_tex_dither = setup_program_from_source(STR_LIT("blit texture dither"), f_shader_src_tex_dither);

    program_tex_covered = setup_program_from_source(STR_LIT("blit texture covered"), f_shader_src_tex_covered);

    program_col = setup_program_from_source(STR_LIT("blit color"), f_shader_src_col);
    uniform_loc_color = glGetUniformLocation(program_col, "u_color");
}

void shutdown() {
    if (program_tex) glDeleteProgram(program_tex);
    if (program_tex_dither) glDeleteProgram(program_tex_dither);
    if (program_tex_covered) glDeleteProgram(program_tex_covered);
    if (program_col) glDeleteProgram(program_col);
    program_tex = program_tex_dither = program_tex_covered = program_col = 0;
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
void initialize() {
    {
        gl.temporal.with_motion_blur.program = setup_program_from_source(STR_LIT("temporal aa + motion-blur"), {(const char*)temporal_frag, temporal_frag_size});
        gl.temporal.no_motion_blur.program   = setup_program_from_source(STR_LIT("temporal aa"), {(const char*)temporal_frag, temporal_frag_size}, STR_LIT("#define USE_MOTION_BLUR 0\n"));

        gl.temporal.with_motion_blur.uniform_loc.tex_linear_depth = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_tex_linear_depth");
        gl.temporal.with_motion_blur.uniform_loc.tex_main = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_tex_main");
        gl.temporal.with_motion_blur.uniform_loc.tex_prev = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_tex_prev");
        gl.temporal.with_motion_blur.uniform_loc.tex_vel = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_tex_vel");
        gl.temporal.with_motion_blur.uniform_loc.tex_vel_neighbormax = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_tex_vel_neighbormax");
        gl.temporal.with_motion_blur.uniform_loc.texel_size = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_texel_size");
        gl.temporal.with_motion_blur.uniform_loc.jitter_uv = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_jitter_uv");
        gl.temporal.with_motion_blur.uniform_loc.time = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_time");
        gl.temporal.with_motion_blur.uniform_loc.feedback_min = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_feedback_min");
        gl.temporal.with_motion_blur.uniform_loc.feedback_max = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_feedback_max");
        gl.temporal.with_motion_blur.uniform_loc.motion_scale = glGetUniformLocation(gl.temporal.with_motion_blur.program, "u_motion_scale");

        gl.temporal.no_motion_blur.uniform_loc.tex_linear_depth = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_tex_linear_depth");
        gl.temporal.no_motion_blur.uniform_loc.tex_main = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_tex_main");
        gl.temporal.no_motion_blur.uniform_loc.tex_prev = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_tex_prev");
        gl.temporal.no_motion_blur.uniform_loc.tex_vel = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_tex_vel");
        gl.temporal.no_motion_blur.uniform_loc.texel_size = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_texel_size");
        gl.temporal.no_motion_blur.uniform_loc.jitter_uv = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_jitter_uv");
        gl.temporal.no_motion_blur.uniform_loc.time = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_time");
        gl.temporal.no_motion_blur.uniform_loc.feedback_min = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_feedback_min");
        gl.temporal.no_motion_blur.uniform_loc.feedback_max = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_feedback_max");
        gl.temporal.no_motion_blur.uniform_loc.motion_scale = glGetUniformLocation(gl.temporal.no_motion_blur.program, "u_motion_scale");
    }
}

void shutdown() {
    if (gl.temporal.with_motion_blur.program) glDeleteProgram(gl.temporal.with_motion_blur.program);
    if (gl.temporal.no_motion_blur.program)   glDeleteProgram(gl.temporal.no_motion_blur.program);
    gl.temporal.with_motion_blur.program = gl.temporal.no_motion_blur.program = 0;
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
    highlight::initialize();
    hsv::initialize();
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
void initialize(int width, int height) {
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
    ensure_texture_2d(&gl.rt.tex_history_prev, GL_RGB10_A2, width, height, GL_RGBA, GL_UNSIGNED_INT_2_10_10_10_REV);

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

    ssao::initialize_targets(width, height);
    dof::initialize_targets(width, height);
    velocity::initialize_targets(width, height);
}

void shutdown() {
    ssao::shutdown();
    dof::shutdown();
    velocity::shutdown();
    highlight::shutdown();
    hsv::shutdown();
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

// Half-res scale-free SSAO, multiplied into the currently bound color target.
// z_far is the far clip distance: background pixels hold exactly that linear depth and are skipped.
// frame only rotates the sample pattern; pass 0 unless something (TAA) integrates over frames.
void compute_ssao(GLuint linear_depth_tex, GLuint normal_tex, const float proj_mat[4][4], float z_far, float intensity, int frame) {
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(normal_tex));

    GLResetState reset_state = {};
    record_gl_reset_state(&reset_state);

    const int width  = reset_state.viewport[2];
    const int height = reset_state.viewport[3];
    const int half_w = MAX(width  / 2, 1);
    const int half_h = MAX(height / 2, 1);

    const int    variant   = is_orthographic_proj_matrix(proj_mat) ? 0 : 1;
    const vec4_t proj_info = compute_proj_info(proj_mat);
    const float  px_scale  = proj_info.y / (float)height;  // world size of a pixel at depth 1 (persp) or absolute (ortho)
    const float  z_max     = z_far * 0.99f;
    const float  r_max     = MAX(ssao::R_MAX_FRACTION * (float)height, 2.0f * ssao::R_MIN_PX);

    glBindVertexArray(gl.vao);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.ssao.fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glViewport(0, 0, half_w, half_h);
    glScissor(0, 0, half_w, half_h);

    PUSH_GPU_SECTION("AO")
    {
        const GLuint program = gl.ssao.program_ao[variant];
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.ssao.tex[0], 0);
        glUseProgram(program);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, linear_depth_tex);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, normal_tex);
        glUniform1i(glGetUniformLocation(program, "u_tex_linear_depth"), 0);
        glUniform1i(glGetUniformLocation(program, "u_tex_normal"), 1);
        glUniform4fv(glGetUniformLocation(program, "u_proj_info"), 1, &proj_info.x);
        glUniform2f(glGetUniformLocation(program, "u_full_res"), (float)width, (float)height);
        glUniform1f(glGetUniformLocation(program, "u_px_scale"), px_scale);
        glUniform1f(glGetUniformLocation(program, "u_r_min"), ssao::R_MIN_PX);
        glUniform1f(glGetUniformLocation(program, "u_r_max"), r_max);
        glUniform1f(glGetUniformLocation(program, "u_intensity"), MAX(intensity, 0.0f));
        glUniform1f(glGetUniformLocation(program, "u_z_max"), z_max);
        glUniform1i(glGetUniformLocation(program, "u_frame"), frame);
        glDrawArrays(GL_TRIANGLES, 0, 3);
    }
    POP_GPU_SECTION()

    PUSH_GPU_SECTION("Blur")
    {
        const GLuint program = gl.ssao.program_blur[variant];
        glUseProgram(program);
        glUniform1i(glGetUniformLocation(program, "u_tex"), 0);
        glUniform1i(glGetUniformLocation(program, "u_tex_normal"), 1);
        glUniform4fv(glGetUniformLocation(program, "u_proj_info"), 1, &proj_info.x);
        glUniform2f(glGetUniformLocation(program, "u_full_res"), (float)width, (float)height);
        glUniform1f(glGetUniformLocation(program, "u_px_scale"), px_scale);
        glUniform1f(glGetUniformLocation(program, "u_z_max"), z_max);
        const GLint loc_dir = glGetUniformLocation(program, "u_dir");

        // normal_tex is still bound to unit 1
        glActiveTexture(GL_TEXTURE0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.ssao.tex[1], 0);
        glBindTexture(GL_TEXTURE_2D, gl.ssao.tex[0]);
        glUniform2i(loc_dir, 1, 0);
        glDrawArrays(GL_TRIANGLES, 0, 3);

        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.ssao.tex[0], 0);
        glBindTexture(GL_TEXTURE_2D, gl.ssao.tex[1]);
        glUniform2i(loc_dir, 0, 1);
        glDrawArrays(GL_TRIANGLES, 0, 3);
    }
    POP_GPU_SECTION()

    PUSH_GPU_SECTION("Upsample + apply")
    {
        reset_gl_state(reset_state);
        const GLuint program = gl.ssao.program_upsample[variant];
        glUseProgram(program);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, gl.ssao.tex[0]);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, linear_depth_tex);
        glUniform1i(glGetUniformLocation(program, "u_tex_ao"), 0);
        glUniform1i(glGetUniformLocation(program, "u_tex_linear_depth"), 1);
        glUniform1f(glGetUniformLocation(program, "u_px_scale"), px_scale);
        glUniform1f(glGetUniformLocation(program, "u_z_max"), z_max);

        glEnable(GL_BLEND);
        glBlendFunc(GL_ZERO, GL_SRC_COLOR);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glDisable(GL_BLEND);
    }
    POP_GPU_SECTION()

    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(0);
}

static vec4_t compute_proj_info(const float proj_mat[4][4]) {
    if (!is_orthographic_proj_matrix(proj_mat)) {
        return vec4_t{
            2.0f / (proj_mat[0][0]),
            2.0f / (proj_mat[1][1]),
            -(1.0f - proj_mat[2][0]) / proj_mat[0][0],
            -(1.0f + proj_mat[2][1]) / proj_mat[1][1]
        };
    }

    return vec4_t{
        2.0f / (proj_mat[0][0]),
        2.0f / (proj_mat[1][1]),
        -(1.0f + proj_mat[3][0]) / proj_mat[0][0],
        -(1.0f - proj_mat[3][1]) / proj_mat[1][1]
    };
}

static void compose_deferred(GLuint linear_depth_tex, GLuint color_tex, GLuint normal_tex, const float proj_mat[4][4], bool orthographic, const vec3_t bg_color, float time) {
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(color_tex));
    ASSERT(glIsTexture(normal_tex));

    const vec3_t env_radiance = bg_color * 0.25f;
    const vec3_t dir_radiance = {10, 10, 10};
    const float roughness = 0.4f;
    const float F0 = 0.04f;
    const vec3_t L = {0.57735026918962576451f, 0.57735026918962576451f, 0.57735026918962576451f}; // 1.0 / sqrt(3)

    compose::ubo_data_t data = {
        .proj_info = compute_proj_info(proj_mat),
        .bg_color = bg_color,
        .time = time,
        .env_radiance = env_radiance,
        .roughness = roughness,
        .dir_radiance = dir_radiance,
        .F0 = F0,
        .light_dir = L,
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

    auto& variant = orthographic ? compose::compose.ortho : compose::compose.persp;
    GLuint program = variant.program;

    glUseProgram(program);

    glBindBufferBase(GL_UNIFORM_BUFFER, 0, compose::compose.ubo);
    glUniformBlockBinding(program, variant.uniform_loc.uniform_data, 0);
    glUniform1i (variant.uniform_loc.texture_depth, 0);
    glUniform1i (variant.uniform_loc.texture_color, 1);
    glUniform1i (variant.uniform_loc.texture_normal, 2);

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);

    glUseProgram(0);
}

void highlight_selection(GLuint atom_idx_tex, GLuint selection_buffer, const vec3_t& highlight, const vec3_t& selection, const vec3_t& outline) {
    ASSERT(glIsTexture(atom_idx_tex));
    ASSERT(glIsBuffer(selection_buffer));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, atom_idx_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, highlight::highlight.selection_texture);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, selection_buffer);

    glUseProgram(highlight::highlight.program);
    glUniform1i(highlight::highlight.uniform_loc.texture_atom_idx, 0);
    glUniform1i(highlight::highlight.uniform_loc.buffer_selection, 1);
    glUniform3fv(highlight::highlight.uniform_loc.highlight, 1, &highlight.x);
    glUniform3fv(highlight::highlight.uniform_loc.selection, 1, &selection.x);
    glUniform3fv(highlight::highlight.uniform_loc.outline, 1, &outline.x);
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

void apply_tonemapping(GLuint color_tex, Tonemapping tonemapping, float exposure, float gamma) {
    ASSERT(glIsTexture(color_tex));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    switch (tonemapping) {
        case Tonemapping_ExposureGamma:
            glUseProgram(tonemapping::exposure_gamma.program);
            glUniform1i(tonemapping::exposure_gamma.uniform_loc.texture, 0);
            glUniform1f(tonemapping::exposure_gamma.uniform_loc.exposure, exposure);
            glUniform1f(tonemapping::exposure_gamma.uniform_loc.gamma, gamma);
            break;
        case Tonemapping_Filmic:
            glUseProgram(tonemapping::filmic.program);
            glUniform1i(tonemapping::filmic.uniform_loc.texture, 0);
            glUniform1f(tonemapping::filmic.uniform_loc.exposure, exposure);
            glUniform1f(tonemapping::filmic.uniform_loc.gamma, gamma);
            break;
        case Tonemapping_ACES:
            glUseProgram(tonemapping::aces.program);
            glUniform1i(tonemapping::aces.uniform_loc.texture, 0);
            glUniform1f(tonemapping::aces.uniform_loc.exposure, exposure);
            glUniform1f(tonemapping::aces.uniform_loc.gamma, gamma);
            break;
        case Tonemapping_Passthrough:
        default:
            glUseProgram(tonemapping::passthrough.program);
            glUniform1i(tonemapping::passthrough.uniform_loc.texture, 0);
            break;
    }

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

void apply_temporal_aa(GLuint linear_depth_tex, GLuint color_tex, GLuint prev_history_tex, GLuint velocity_tex, GLuint velocity_neighbormax_tex, const vec2_t& curr_jitter,
                       const vec2_t& prev_jitter, float feedback_min, float feedback_max, float motion_scale, float time) {
    ASSERT(glIsTexture(linear_depth_tex));
    ASSERT(glIsTexture(color_tex));
    ASSERT(glIsTexture(velocity_tex));
    ASSERT(glIsTexture(velocity_neighbormax_tex));

    const vec2_t res = {(float)gl.tex_width, (float)gl.tex_height};
    const vec2_t inv_res = 1.0f / res;
    const vec4_t texel_size = vec4_t{inv_res.x, inv_res.y, res.x, res.y};
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

    if (motion_scale != 0.f) {
        glUseProgram(gl.temporal.with_motion_blur.program);

        glUniform1i(gl.temporal.with_motion_blur.uniform_loc.tex_linear_depth, 0);
        glUniform1i(gl.temporal.with_motion_blur.uniform_loc.tex_main, 1);
        glUniform1i(gl.temporal.with_motion_blur.uniform_loc.tex_prev, 2);
        glUniform1i(gl.temporal.with_motion_blur.uniform_loc.tex_vel, 3);
        glUniform1i(gl.temporal.with_motion_blur.uniform_loc.tex_vel_neighbormax, 4);

        glUniform4fv(gl.temporal.with_motion_blur.uniform_loc.texel_size, 1, &texel_size.x);
        glUniform4fv(gl.temporal.with_motion_blur.uniform_loc.jitter_uv, 1, &jitter_uv.x);
        glUniform1f(gl.temporal.with_motion_blur.uniform_loc.time, time);
        glUniform1f(gl.temporal.with_motion_blur.uniform_loc.feedback_min, feedback_min);
        glUniform1f(gl.temporal.with_motion_blur.uniform_loc.feedback_max, feedback_max);
        glUniform1f(gl.temporal.with_motion_blur.uniform_loc.motion_scale, motion_scale);
    } else {
        glUseProgram(gl.temporal.no_motion_blur.program);

        glUniform1i(gl.temporal.no_motion_blur.uniform_loc.tex_linear_depth, 0);
        glUniform1i(gl.temporal.no_motion_blur.uniform_loc.tex_main, 1);
        glUniform1i(gl.temporal.no_motion_blur.uniform_loc.tex_prev, 2);
        glUniform1i(gl.temporal.no_motion_blur.uniform_loc.tex_vel, 3);

        glUniform4fv(gl.temporal.no_motion_blur.uniform_loc.texel_size, 1, &texel_size.x);
        glUniform4fv(gl.temporal.no_motion_blur.uniform_loc.jitter_uv, 1, &jitter_uv.x);
        glUniform1f(gl.temporal.no_motion_blur.uniform_loc.time, time);
        glUniform1f(gl.temporal.no_motion_blur.uniform_loc.feedback_min, feedback_min);
        glUniform1f(gl.temporal.no_motion_blur.uniform_loc.feedback_max, feedback_max);
        glUniform1f(gl.temporal.no_motion_blur.uniform_loc.motion_scale, motion_scale);
    }

    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void scale_hsv(GLuint color_tex, vec3_t hsv_scale) {
    GLResetState reset_state = {};
    record_gl_reset_state(&reset_state);

    GLint w, h;

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH,  &w);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

    glBindVertexArray(gl.vao);

    glViewport(0, 0, w, h);
    glScissor(0, 0, w, h);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.rt.scratch_fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.rt.tex_rgba8, 0);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glUseProgram(hsv::gl.program);
    glUniform1i(hsv::gl.uniform_loc.texture_color, 0);
    glUniform3fv(hsv::gl.uniform_loc.hsv_scale, 1, &hsv_scale.x);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, color_tex, 0);
    glDrawBuffer(GL_COLOR_ATTACHMENT1);
    blit_texture(gl.rt.tex_rgba8);

    reset_gl_state(reset_state);
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

// tex scaled by (1 - alpha of cover), both premultiplied
static void blit_texture_covered(GLuint tex, GLuint cover) {
    ASSERT(glIsTexture(tex));
    ASSERT(glIsTexture(cover));
    glUseProgram(blit::program_tex_covered);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, cover);
    glActiveTexture(GL_TEXTURE0);
    glUniform1i(glGetUniformLocation(blit::program_tex_covered, "u_texture"), 0);
    glUniform1i(glGetUniformLocation(blit::program_tex_covered, "u_cover"), 1);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void blit_color(vec4_t color) {
    glUseProgram(blit::program_col);
    glUniform4fv(blit::uniform_loc_color, 1, &color.x);
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
    const bool do_tonemap = settings.tonemap.enabled;
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

    // Rotates the SSAO sample pattern only when TAA is there to integrate it; otherwise the image is stable.
    static int frame_index = 0;
    frame_index = (frame_index + 1) & 1023;
    const int noise_frame = do_taa ? frame_index : 0;

    GLResetState reset_state = {};
    record_gl_reset_state(&reset_state);

    int width = reset_state.viewport[2];
    int height = reset_state.viewport[3];

    if (width > (int)gl.tex_width || height > (int)gl.tex_height) {
        initialize(width, height);
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

    if (do_ssao) {
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
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, in.history, 0);
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

    if (do_compose) {
        PUSH_GPU_SECTION("Compose")
        compose_deferred(gl.linear_depth.texture, in.color, in.normal, view_param.matrix.curr.proj.elem, ortho, settings.background_color, time);
        POP_GPU_SECTION()
    } else {
        // Seed transient pipeline from external color input when compose stage is disabled.
        blit_texture(in.color);
    }

    if (do_ssao) {
        PUSH_GPU_SECTION("SSAO")
        compute_ssao(gl.linear_depth.texture, in.normal, view_param.matrix.curr.proj.elem, far_dist, settings.ssao.intensity, noise_frame);
        POP_GPU_SECTION()
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

    if (do_tonemap) {
        PUSH_GPU_SECTION("Tonemapping")
        swap_target();
        // From here on the values are display referred: switch the ping-pong pair to the RGB10_A2 targets
        pair = gl.rt.tex_ldr;
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.rt.tex_ldr[0], 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gl.rt.tex_ldr[1], 0);
        glDrawBuffer(dst_buffer);
        Tonemapping tonemapper = settings.tonemap.enabled ? to_legacy_tonemapper(settings.tonemap.mode) : Tonemapping_Passthrough;
        apply_tonemapping(src_texture, tonemapper, settings.tonemap.exposure, settings.tonemap.gamma);
        POP_GPU_SECTION()
    }

    // LDR effects -----------------------------------------------------------

    if (do_transparency) {
        PUSH_GPU_SECTION("Add Transparency")
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);    // premultiplied
        if (in.transparency_under_hdr && do_transparency_hdr) {
            blit_texture_covered(in.transparency, in.transparency_hdr);
        } else {
            blit_texture(in.transparency);
        }
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

    if (do_taa) {
        GLuint history_prev = in.history_prev;
        if (!history_prev) {
            copy_texture_2d(gl.rt.tex_history_prev, in.history, width, height);
            history_prev = gl.rt.tex_history_prev;
        }

        swap_target();
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

        apply_temporal_aa(gl.linear_depth.texture, src_texture, history_prev, in.velocity, gl.velocity.tex_neighbormax, view_param.jitter.curr, view_param.jitter.prev,
                          feedback_min, feedback_max, motion_scale, time);

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
        blit_texture_dither(src_texture, time);
    }

    if (last_depth_test) glEnable(GL_DEPTH_TEST);
    glDepthMask(1);
    glColorMask(1, 1, 1, 1);
}

}  // namespace postprocessing

namespace postprocess_pipeline {

void initialize(int width, int height) {
    // glTexStorage2D rejects zero dimensions outright (GL_INVALID_VALUE), so guard the
    // 0x0 startup/iconified case here as well. See gbuffer_init().
    postprocessing::initialize(MAX(width, 1), MAX(height, 1));
}

void shutdown() {
    postprocessing::shutdown();
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

}  // namespace postprocess_pipeline
