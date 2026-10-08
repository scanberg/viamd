#include "volumerender_utils.h"

#include <gfx/gl.h>
#include <gfx/gl_utils.h>
#include <gfx/postprocessing_utils.h>
#include <color_utils.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_vec_math.h>

#include <implot.h>

#include <shaders.inl>

#define PUSH_GPU_SECTION(lbl)                                                                       \
    {                                                                                               \
        if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); \
    }
#define POP_GPU_SECTION()                       \
    {                                           \
        if (glPopDebugGroup) glPopDebugGroup(); \
    }

static constexpr str_t v_shader_src_fs_quad = STR_LIT(
    R"(
#version 410 core

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

namespace volume {

static struct {
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint ubo = 0;
    GLuint fbo = 0;
    GLuint ssbo = 0;

    GLuint tex_entry  = 0;
    GLuint tex_exit   = 0;
    GLuint tex_result = 0;

    uint32_t width  = 0;
    uint32_t height = 0;

    struct {
        GLuint entry_exit = 0;
        GLuint entry_exit_depth = 0;
        GLuint dvr_only = 0;
        GLuint splat_color = 0;
    } program;

    struct {
        int major = 0;
        int minor = 0;
    } version;
} gl;

struct UniformData {
    mat4_t view_to_model_mat;
    mat4_t model_to_view_mat;
    mat4_t inv_proj_mat;
    mat4_t model_view_proj_mat;

    vec2_t inv_res;
    float  time;
    float  gamma;

    vec3_t clip_volume_min;
    float  tf_min;
    vec3_t clip_volume_max;
    float  tf_inv_ext;

    vec3_t gradient_spacing_world_space;
    float  exposure;
    mat4_t gradient_spacing_tex_space;

    vec3_t env_radiance;
    float  roughness;
    vec3_t dir_radiance;
    float  F0;
};

// -----------------------------------------------------------------------------
// GPU timings
// -----------------------------------------------------------------------------

static constexpr int    TIMER_QUERY_COUNT     = 64;
static constexpr int    TIMER_PUBLISH_FRAMES  = 30;

static struct {
    struct {
        GLuint id = 0;
        int    stage = 0;
        bool   pending = false;
    } query[TIMER_QUERY_COUNT];

    int    active = -1;     // the query recording right now: GL_TIME_ELAPSED queries cannot nest
    double accum_ms[TimingStage_Count] = {};
    int    accum_frames = 0;
    GpuTimings published = {};
} timing;

static void timer_begin(TimingStage stage) {
    if (timing.active != -1) return;
    for (int i = 0; i < TIMER_QUERY_COUNT; ++i) {
        if (timing.query[i].id && !timing.query[i].pending) {
            timing.query[i].stage = stage;
            timing.active = i;
            glBeginQuery(GL_TIME_ELAPSED, timing.query[i].id);
            return;
        }
    }
}

static void timer_end() {
    if (timing.active == -1) return;
    glEndQuery(GL_TIME_ELAPSED);
    timing.query[timing.active].pending = true;
    timing.active = -1;
}

void timings_new_frame() {
    for (int i = 0; i < TIMER_QUERY_COUNT; ++i) {
        if (!timing.query[i].pending) continue;
        GLint available = 0;
        glGetQueryObjectiv(timing.query[i].id, GL_QUERY_RESULT_AVAILABLE, &available);
        if (!available) continue;
        GLuint64 ns = 0;
        glGetQueryObjectui64v(timing.query[i].id, GL_QUERY_RESULT, &ns);
        timing.accum_ms[timing.query[i].stage] += (double)ns * 1.0e-6;
        timing.query[i].pending = false;
    }

    timing.accum_frames += 1;
    if (timing.accum_frames >= TIMER_PUBLISH_FRAMES) {
        GpuTimings t = {};
        for (int i = 0; i < TimingStage_Count; ++i) {
            t.ms[i] = (float)(timing.accum_ms[i] / timing.accum_frames);
            t.total_ms += t.ms[i];
            timing.accum_ms[i] = 0.0;
        }
        timing.published = t;
        timing.accum_frames = 0;
    }
}

GpuTimings timings_get() {
    return timing.published;
}

const char* timing_stage_name(TimingStage stage) {
    switch (stage) {
    case TimingStage_BlockMinMax: return "Block min/max";
    case TimingStage_EntryExit:   return "Entry / exit";
    case TimingStage_Raycast:     return "Raycast";
    default: return "";
    }
}

// -----------------------------------------------------------------------------
// Isosurfaces
// -----------------------------------------------------------------------------

static constexpr int   ISO_MAX_COUNT = 8;
// Optical densities are given per 1/150 world unit: the extinction per world unit inside a surface is
// tau * ISO_OPTICAL_SCALE. Kept from the previous raycaster so existing settings look the same.
static constexpr float ISO_OPTICAL_SCALE = 150.0f;
static constexpr float ISO_SAMPLES_PER_VOXEL = 1.0f;

enum IsoVariant {
    IsoVariant_Uniform,
    IsoVariant_ColorVolume,
    IsoVariant_Field,
    IsoVariant_Count
};

struct IsoProgram {
    GLuint program = 0;
    GLint  loc_iso_values = -1;
    GLint  loc_iso_colors = -1;
    GLint  loc_iso_tau = -1;
    GLint  loc_iso_count = -1;
    GLint  loc_tex_volume = -1;
    GLint  loc_tex_depth = -1;
    GLint  loc_tex_color_volume = -1;
    GLint  loc_tex_field = -1;
    GLint  loc_tex_field_colormap = -1;
    GLint  loc_tex_entry = -1;
    GLint  loc_tex_exit = -1;
    GLint  block_index = -1;
};

// std140, mirrors IsoUniforms in isosurface.frag
struct IsoUniformData {
    mat4_t clip_to_model;
    mat4_t model_to_view;
    mat4_t grad_offsets;

    vec3_t clip_min;
    float  samples_per_voxel;
    vec3_t clip_max;
    float  use_depth;

    vec3_t env_radiance;
    float  roughness;
    vec3_t dir_radiance;
    float  F0;
    vec3_t light_dir;
    float  field_beg;

    vec2_t inv_res;
    float  field_inv_ext;
    float  optical_scale;

    float  use_proxy;
    float  entry_from_near;
    float  pad0;
    float  pad1;
};
static_assert(sizeof(IsoUniformData) == 3 * 64 + 7 * 16, "IsoUniformData must match the std140 layout of IsoUniforms");

// Block min/max grids: what the isosurface rays skip empty space with. One per density volume, kept
// here keyed by the volume's texture and rebuilt when notify_data_changed() says its texels changed.
static constexpr int BLOCK_MIN_SIZE     = 8;    // voxels per block side, at least
static constexpr int BLOCK_MAX_PER_AXIS = 32;   // larger volumes get larger blocks: the proxy draws one box per block
static constexpr int BLOCK_CACHE_SIZE   = 64;

struct BlockGrid {
    GLuint   volume = 0;
    GLuint   minmax = 0;            // RG32F, one texel per block
    uint64_t version = 0;           // of the volume's texels
    uint64_t built_version = 0;     // the version minmax was built from
    int      dim[3] = {};           // of the volume when built
    int      block_dim[3] = {};
    int      block_size = 0;
    uint64_t last_used = 0;
};

static struct {
    IsoProgram prog[IsoVariant_Count];
    GLuint ubo = 0;
    GLuint fbo = 0;
    GLuint vao = 0;     // empty: full screen triangle and proxy boxes both come from gl_VertexID

    BlockGrid grid[BLOCK_CACHE_SIZE];
    uint64_t  use_counter = 0;

    struct {
        GLuint program = 0;
        GLint  loc_volume = -1;
        GLint  loc_layer = -1;
        GLint  loc_block_size = -1;
        GLuint fbo = 0;
    } minmax;

    struct {
        GLuint program = 0;
        GLint  loc_minmax = -1;
        GLint  loc_model_to_clip = -1;
        GLint  loc_block_dim = -1;
        GLint  loc_block_ext = -1;
        GLint  loc_clip_min = -1;
        GLint  loc_clip_max = -1;
        GLint  loc_iso_values = -1;
        GLint  loc_iso_tau = -1;
        GLint  loc_iso_count = -1;
        GLuint fbo = 0;
        GLuint tex_entry = 0;       // DEPTH_COMPONENT32F, nearest proxy depth
        GLuint tex_exit = 0;        // DEPTH_COMPONENT32F, farthest proxy depth
        int    width = 0;
        int    height = 0;
    } proxy;
} iso;

static void iso_program_setup(IsoProgram* p, GLuint v_shader, str_t defines) {
    GLuint f_shader = gl::compile_shader_from_source({(const char*)isosurface_frag, isosurface_frag_size}, GL_FRAGMENT_SHADER, defines);
    if (!f_shader) {
        MD_LOG_ERROR("Isosurface shader compilation failed (" STR_FMT "), keeping the previous program", STR_ARG(defines));
        return;
    }
    if (!p->program) p->program = glCreateProgram();
    const GLuint shaders[] = {v_shader, f_shader};
    gl::attach_link_detach(p->program, shaders, (int)ARRAY_SIZE(shaders));
    glDeleteShader(f_shader);

    const GLuint prog = p->program;
    p->loc_iso_values         = glGetUniformLocation(prog, "u_iso_values");
    p->loc_iso_colors         = glGetUniformLocation(prog, "u_iso_colors");
    p->loc_iso_tau            = glGetUniformLocation(prog, "u_iso_tau");
    p->loc_iso_count          = glGetUniformLocation(prog, "u_iso_count");
    p->loc_tex_volume         = glGetUniformLocation(prog, "u_tex_volume");
    p->loc_tex_depth          = glGetUniformLocation(prog, "u_tex_depth");
    p->loc_tex_color_volume   = glGetUniformLocation(prog, "u_tex_color_volume");
    p->loc_tex_field          = glGetUniformLocation(prog, "u_tex_field");
    p->loc_tex_field_colormap = glGetUniformLocation(prog, "u_tex_field_colormap");
    p->loc_tex_entry          = glGetUniformLocation(prog, "u_tex_entry");
    p->loc_tex_exit           = glGetUniformLocation(prog, "u_tex_exit");
    p->block_index            = glGetUniformBlockIndex(prog, "IsoUniforms");
}

static void iso_initialize() {
    GLuint v_shader = gl::compile_shader_from_source(v_shader_src_fs_quad, GL_VERTEX_SHADER);
    if (!v_shader) {
        MD_LOG_ERROR("Isosurface vertex shader compilation failed");
        return;
    }
    iso_program_setup(&iso.prog[IsoVariant_Uniform],     v_shader, STR_LIT(""));
    iso_program_setup(&iso.prog[IsoVariant_ColorVolume], v_shader, STR_LIT("#define USE_COLOR_VOLUME"));
    iso_program_setup(&iso.prog[IsoVariant_Field],       v_shader, STR_LIT("#define USE_FIELD"));

    {
        GLuint f_shader = gl::compile_shader_from_source({(const char*)block_minmax_frag, block_minmax_frag_size}, GL_FRAGMENT_SHADER);
        if (f_shader) {
            if (!iso.minmax.program) iso.minmax.program = glCreateProgram();
            const GLuint shaders[] = {v_shader, f_shader};
            gl::attach_link_detach(iso.minmax.program, shaders, (int)ARRAY_SIZE(shaders));
            glDeleteShader(f_shader);
            iso.minmax.loc_volume     = glGetUniformLocation(iso.minmax.program, "u_tex_volume");
            iso.minmax.loc_layer      = glGetUniformLocation(iso.minmax.program, "u_layer");
            iso.minmax.loc_block_size = glGetUniformLocation(iso.minmax.program, "u_block_size");
        } else {
            MD_LOG_ERROR("Block min/max shader compilation failed, isosurfaces will not skip empty space");
        }
    }
    glDeleteShader(v_shader);

    {
        GLuint v_proxy = gl::compile_shader_from_source({(const char*)block_proxy_vert, block_proxy_vert_size}, GL_VERTEX_SHADER);
        GLuint f_proxy = gl::compile_shader_from_source({(const char*)block_proxy_frag, block_proxy_frag_size}, GL_FRAGMENT_SHADER);
        if (v_proxy && f_proxy) {
            if (!iso.proxy.program) iso.proxy.program = glCreateProgram();
            const GLuint shaders[] = {v_proxy, f_proxy};
            gl::attach_link_detach(iso.proxy.program, shaders, (int)ARRAY_SIZE(shaders));
            const GLuint prog = iso.proxy.program;
            iso.proxy.loc_minmax        = glGetUniformLocation(prog, "u_tex_minmax");
            iso.proxy.loc_model_to_clip = glGetUniformLocation(prog, "u_model_to_clip");
            iso.proxy.loc_block_dim     = glGetUniformLocation(prog, "u_block_dim");
            iso.proxy.loc_block_ext     = glGetUniformLocation(prog, "u_block_ext");
            iso.proxy.loc_clip_min      = glGetUniformLocation(prog, "u_clip_min");
            iso.proxy.loc_clip_max      = glGetUniformLocation(prog, "u_clip_max");
            iso.proxy.loc_iso_values    = glGetUniformLocation(prog, "u_iso_values");
            iso.proxy.loc_iso_tau       = glGetUniformLocation(prog, "u_iso_tau");
            iso.proxy.loc_iso_count     = glGetUniformLocation(prog, "u_iso_count");
        } else {
            MD_LOG_ERROR("Block proxy shader compilation failed, isosurfaces will not skip empty space");
        }
        if (v_proxy) glDeleteShader(v_proxy);
        if (f_proxy) glDeleteShader(f_proxy);
    }

    if (!iso.ubo) {
        glGenBuffers(1, &iso.ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, iso.ubo);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(IsoUniformData), 0, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
    }
    if (!iso.fbo) {
        glGenFramebuffers(1, &iso.fbo);
    }
    if (!iso.vao) {
        glGenVertexArrays(1, &iso.vao);
    }
    if (!iso.minmax.fbo) {
        glGenFramebuffers(1, &iso.minmax.fbo);
    }
    if (!iso.proxy.fbo) {
        glGenFramebuffers(1, &iso.proxy.fbo);
    }
}

static BlockGrid* block_grid_find(GLuint volume) {
    for (int i = 0; i < BLOCK_CACHE_SIZE; ++i) {
        if (iso.grid[i].volume == volume) return &iso.grid[i];
    }
    return NULL;
}

static BlockGrid* block_grid_acquire(GLuint volume) {
    if (BlockGrid* g = block_grid_find(volume)) {
        return g;
    }
    // A free slot, else the one used longest ago
    BlockGrid* slot = NULL;
    for (int i = 0; i < BLOCK_CACHE_SIZE; ++i) {
        if (!iso.grid[i].volume) {
            slot = &iso.grid[i];
            break;
        }
        if (!slot || iso.grid[i].last_used < slot->last_used) {
            slot = &iso.grid[i];
        }
    }
    const GLuint minmax = slot->minmax;     // the texture object is reused, its storage respecified on build
    *slot = BlockGrid{};
    slot->minmax    = minmax;
    slot->volume    = volume;
    slot->version   = 1;
    slot->last_used = ++iso.use_counter;
    return slot;
}

void notify_data_changed(uint32_t volume_texture) {
    if (!volume_texture) return;
    block_grid_acquire(volume_texture)->version += 1;
}

uint64_t data_version(uint32_t volume_texture) {
    const BlockGrid* g = block_grid_find(volume_texture);
    return g ? g->version : 0;
}

// The block grid of the volume, (re)built if its texels changed since. NULL when it cannot be made.
static const BlockGrid* block_grid_update(GLuint volume, const int dim[3]) {
    if (!iso.minmax.program || !iso.minmax.fbo) return NULL;

    BlockGrid* g = block_grid_acquire(volume);
    g->last_used = ++iso.use_counter;

    const bool same_dim = g->dim[0] == dim[0] && g->dim[1] == dim[1] && g->dim[2] == dim[2];
    if (g->minmax && same_dim && g->built_version == g->version) {
        return g;
    }

    int bs = BLOCK_MIN_SIZE;
    while (DIV_UP(MAX(dim[0], MAX(dim[1], dim[2])), bs) > BLOCK_MAX_PER_AXIS) {
        bs *= 2;
    }
    const int bd[3] = { DIV_UP(dim[0], bs), DIV_UP(dim[1], bs), DIV_UP(dim[2], bs) };

    if (!g->minmax) {
        glGenTextures(1, &g->minmax);
    }
    if (!same_dim || g->block_size != bs || g->built_version == 0) {
        glBindTexture(GL_TEXTURE_3D, g->minmax);
        glTexImage3D(GL_TEXTURE_3D, 0, GL_RG32F, bd[0], bd[1], bd[2], 0, GL_RG, GL_FLOAT, NULL);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_3D, 0);
    }

    PUSH_GPU_SECTION("ISO BLOCK MIN/MAX")
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, iso.minmax.fbo);
    glViewport(0, 0, bd[0], bd[1]);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glDisable(GL_SCISSOR_TEST);

    glUseProgram(iso.minmax.program);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, volume);
    glUniform1i(iso.minmax.loc_volume, 0);
    glUniform1i(iso.minmax.loc_block_size, bs);
    glBindVertexArray(iso.vao);

    bool complete = true;
    timer_begin(TimingStage_BlockMinMax);
    for (int z = 0; z < bd[2]; ++z) {
        glFramebufferTextureLayer(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, g->minmax, 0, z);
        if (z == 0) {
            glDrawBuffer(GL_COLOR_ATTACHMENT0);
            const GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
            if (status != GL_FRAMEBUFFER_COMPLETE) {
                MD_LOG_ERROR("Block min/max framebuffer is incomplete (0x%04X)", (unsigned int)status);
                complete = false;
                break;
            }
        }
        glUniform1i(iso.minmax.loc_layer, z);
        glDrawArrays(GL_TRIANGLES, 0, 3);
    }
    timer_end();

    glFramebufferTextureLayer(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, 0, 0, 0);
    glBindVertexArray(0);
    glUseProgram(0);
    POP_GPU_SECTION()

    if (!complete) {
        return NULL;
    }

    MEMCPY(g->dim, dim, sizeof(g->dim));
    MEMCPY(g->block_dim, bd, sizeof(g->block_dim));
    g->block_size = bs;
    g->built_version = g->version;
    return g;
}

// Whether any corner of the clip box lies at or behind the near plane. Then rays start at the near plane
// rather than at the proxy's nearest depth, which cannot see the faces of a box around the camera.
static bool near_plane_cuts(const mat4_t& model_to_clip, vec3_t lo, vec3_t hi) {
    for (int i = 0; i < 8; ++i) {
        const vec4_t p = { (i & 1) ? hi.x : lo.x, (i & 2) ? hi.y : lo.y, (i & 4) ? hi.z : lo.z, 1.0f };
        const vec4_t c = mat4_mul_vec4(model_to_clip, p);
        if (c.z <= -c.w) return true;
    }
    return false;
}

// Rasterizes the boxes of the blocks that can hold a surface into the entry (nearest) and exit (farthest)
// depth textures
static bool proxy_render(const BlockGrid& g, const mat4_t& model_to_clip, vec3_t clip_min, vec3_t clip_max,
                         const float* values, const float* tau, int count, int width, int height) {
    if (!iso.proxy.program || !iso.proxy.fbo) return false;

    if (iso.proxy.width < width || iso.proxy.height < height) {
        iso.proxy.width  = MAX(iso.proxy.width, width);
        iso.proxy.height = MAX(iso.proxy.height, height);
        GLuint* texs[2] = { &iso.proxy.tex_entry, &iso.proxy.tex_exit };
        for (GLuint* t : texs) {
            if (!*t) glGenTextures(1, t);
            glBindTexture(GL_TEXTURE_2D, *t);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32F, iso.proxy.width, iso.proxy.height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
        }
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    const vec3_t block_ext = {
        (float)g.block_size / (float)g.dim[0],
        (float)g.block_size / (float)g.dim[1],
        (float)g.block_size / (float)g.dim[2],
    };
    const int instances = g.block_dim[0] * g.block_dim[1] * g.block_dim[2];

    PUSH_GPU_SECTION("ISO PROXY ENTRY / EXIT")
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, iso.proxy.fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, iso.proxy.tex_entry, 0);
    glDrawBuffer(GL_NONE);
    const GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
    if (status != GL_FRAMEBUFFER_COMPLETE) {
        MD_LOG_ERROR("Isosurface proxy framebuffer is incomplete (0x%04X)", (unsigned int)status);
        POP_GPU_SECTION()
        return false;
    }

    glViewport(0, 0, width, height);
    glDisable(GL_SCISSOR_TEST);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);    // both faces: the nearest and the farthest of every box count
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);

    glUseProgram(iso.proxy.program);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, g.minmax);
    glUniform1i(iso.proxy.loc_minmax, 0);
    glUniformMatrix4fv(iso.proxy.loc_model_to_clip, 1, GL_FALSE, &model_to_clip.elem[0][0]);
    glUniform3i(iso.proxy.loc_block_dim, g.block_dim[0], g.block_dim[1], g.block_dim[2]);
    glUniform3fv(iso.proxy.loc_block_ext, 1, block_ext.elem);
    glUniform3fv(iso.proxy.loc_clip_min, 1, clip_min.elem);
    glUniform3fv(iso.proxy.loc_clip_max, 1, clip_max.elem);
    glUniform1fv(iso.proxy.loc_iso_values, count, values);
    glUniform1fv(iso.proxy.loc_iso_tau, count, tau);
    glUniform1i(iso.proxy.loc_iso_count, count);
    glBindVertexArray(iso.vao);

    timer_begin(TimingStage_EntryExit);
    glClearDepth(1.0);
    glClear(GL_DEPTH_BUFFER_BIT);
    glDepthFunc(GL_LESS);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 14, instances);

    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, iso.proxy.tex_exit, 0);
    glClearDepth(0.0);
    glClear(GL_DEPTH_BUFFER_BIT);
    glDepthFunc(GL_GREATER);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 14, instances);
    timer_end();

    glClearDepth(1.0);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
    glBindVertexArray(0);
    glUseProgram(0);
    POP_GPU_SECTION()
    return true;
}

// The draw target and the state a pass changes, captured so the caller's are put back afterwards
struct SavedState {
    GLint fbo = 0;
    GLint viewport[4] = {};
    GLint draw_buffer[8] = {};
    GLint draw_buffer_count = 0;
    GLboolean depth_test = GL_FALSE;
    GLboolean blend = GL_FALSE;
    GLboolean cull_face = GL_FALSE;
    GLboolean scissor_test = GL_FALSE;
    GLboolean depth_mask = GL_TRUE;
    GLint depth_func = GL_LESS;
    GLint blend_src_rgb = GL_ONE;
    GLint blend_dst_rgb = GL_ZERO;
    GLint blend_src_alpha = GL_ONE;
    GLint blend_dst_alpha = GL_ZERO;
};

static SavedState save_state() {
    SavedState s = {};
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &s.fbo);
    glGetIntegerv(GL_VIEWPORT, s.viewport);
    for (int i = 0; i < 8; ++i) {
        glGetIntegerv(GL_DRAW_BUFFER0 + i, &s.draw_buffer[i]);
        if (s.draw_buffer[i] != GL_NONE) {
            s.draw_buffer_count = i + 1;
        }
    }
    s.depth_test = glIsEnabled(GL_DEPTH_TEST);
    s.blend      = glIsEnabled(GL_BLEND);
    s.cull_face  = glIsEnabled(GL_CULL_FACE);
    s.scissor_test = glIsEnabled(GL_SCISSOR_TEST);
    glGetBooleanv(GL_DEPTH_WRITEMASK, &s.depth_mask);
    glGetIntegerv(GL_DEPTH_FUNC, &s.depth_func);
    glGetIntegerv(GL_BLEND_SRC_RGB,   &s.blend_src_rgb);
    glGetIntegerv(GL_BLEND_DST_RGB,   &s.blend_dst_rgb);
    glGetIntegerv(GL_BLEND_SRC_ALPHA, &s.blend_src_alpha);
    glGetIntegerv(GL_BLEND_DST_ALPHA, &s.blend_dst_alpha);
    return s;
}

static void restore_state(const SavedState& s) {
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, s.fbo);
    glViewport(s.viewport[0], s.viewport[1], s.viewport[2], s.viewport[3]);
    // For the default framebuffer glDrawBuffers rejects the FRONT/BACK/LEFT/RIGHT tokens that
    // glDrawBuffer accepts and that GL_DRAW_BUFFER0 reports back, so it goes through glDrawBuffer there.
    // See reset_gl_state() in postprocessing_utils.cpp for the same fix.
    if (s.fbo == 0) {
        glDrawBuffer(s.draw_buffer_count > 0 ? (GLenum)s.draw_buffer[0] : GL_NONE);
    } else if (s.draw_buffer_count > 0) {
        glDrawBuffers(s.draw_buffer_count, (const GLenum*)s.draw_buffer);
    } else {
        glDrawBuffer(GL_NONE);
    }
    if (s.depth_test) glEnable(GL_DEPTH_TEST); else glDisable(GL_DEPTH_TEST);
    if (s.blend)      glEnable(GL_BLEND);      else glDisable(GL_BLEND);
    if (s.cull_face)  glEnable(GL_CULL_FACE);  else glDisable(GL_CULL_FACE);
    if (s.scissor_test) glEnable(GL_SCISSOR_TEST); else glDisable(GL_SCISSOR_TEST);
    glDepthMask(s.depth_mask);
    glDepthFunc(s.depth_func);
    glBlendFuncSeparate(s.blend_src_rgb, s.blend_dst_rgb, s.blend_src_alpha, s.blend_dst_alpha);
}

// Binds desc's colour target, or keeps the caller's framebuffer when there is none
static bool bind_color_target(GLuint fbo, uint32_t color, uint32_t width, uint32_t height, bool clear) {
    if (!color) return true;
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, color, 0);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
    if (status != GL_FRAMEBUFFER_COMPLETE) {
        MD_LOG_ERROR("Volume render target framebuffer is incomplete (0x%04X)", (unsigned int)status);
        return false;
    }
    glViewport(0, 0, width, height);
    if (clear) {
        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT);
    }
    return true;
}

bool render_isosurfaces(const IsoRenderDesc& desc) {
    const int count = CLAMP((int)desc.iso.count, 0, ISO_MAX_COUNT);
    if (!desc.texture.density_volume || (count > 0 && (!desc.iso.values || !desc.iso.colors)) || desc.render_target.width == 0 || desc.render_target.height == 0) {
        return false;
    }

    // A field takes the place of the colour volume: the two are alternatives, never both
    const bool use_field        = desc.iso.use_field && desc.texture.field_volume && desc.texture.field_colormap;
    const bool use_color_volume = !use_field && desc.iso.use_color_volume && desc.texture.color_volume;
    const IsoProgram& p = iso.prog[use_field ? IsoVariant_Field : use_color_volume ? IsoVariant_ColorVolume : IsoVariant_Uniform];

    int dim[3] = {};
    if (!p.program || !gl::get_texture_dim(dim, desc.texture.density_volume) || dim[0] <= 0 || dim[1] <= 0 || dim[2] <= 0) {
        return false;
    }

    float  values[ISO_MAX_COUNT] = {};
    vec4_t colors[ISO_MAX_COUNT] = {};
    float  tau[ISO_MAX_COUNT] = {};
    MEMCPY(values, desc.iso.values, count * sizeof(float));
    MEMCPY(colors, desc.iso.colors, count * sizeof(vec4_t));
    if (desc.iso.optical_densities) {
        MEMCPY(tau, desc.iso.optical_densities, count * sizeof(float));
    }

    const mat4_t model_to_view = desc.matrix.view * desc.matrix.model;
    const mat4_t view_to_model = mat4_inverse(model_to_view);
    const mat4_t model_to_clip = desc.matrix.proj * model_to_view;
    const mat4_t clip_to_model = mat4_inverse(model_to_clip);

    // One voxel, the smallest of the three axes, along each view axis: the gradient taken over those
    // offsets is the view space normal directly
    const vec3_t voxel_ext = {
        vec3_length(vec3_from_vec4(desc.matrix.model.col[0])) / (float)dim[0],
        vec3_length(vec3_from_vec4(desc.matrix.model.col[1])) / (float)dim[1],
        vec3_length(vec3_from_vec4(desc.matrix.model.col[2])) / (float)dim[2],
    };
    const float h = MAX(1.0e-6f, MIN(voxel_ext.x, MIN(voxel_ext.y, voxel_ext.z)));

    const float n1 = 1.0f;
    const float n2 = desc.shading.ior;
    const float field_ext = desc.field.range_end - desc.field.range_beg;

    IsoUniformData data = {};
    data.clip_to_model     = clip_to_model;
    data.model_to_view     = model_to_view;
    data.grad_offsets      = view_to_model * mat4_scale(h, h, h);
    data.clip_min          = desc.clip_volume.min;
    data.samples_per_voxel = ISO_SAMPLES_PER_VOXEL;
    data.clip_max          = desc.clip_volume.max;
    data.use_depth         = desc.render_target.depth ? 1.0f : 0.0f;
    data.env_radiance      = desc.shading.env_radiance;
    data.roughness         = desc.shading.roughness;
    data.dir_radiance      = desc.shading.dir_radiance;
    data.F0                = powf((n1 - n2) / (n1 + n2), 2.0f);
    data.light_dir         = vec3_normalize(vec3_set(1.0f, 1.0f, 1.0f));   // as in compose_deferred
    data.field_beg         = desc.field.range_beg;
    data.inv_res           = {1.0f / (float)desc.render_target.width, 1.0f / (float)desc.render_target.height};
    data.field_inv_ext     = field_ext != 0.0f ? 1.0f / field_ext : 0.0f;
    data.optical_scale     = ISO_OPTICAL_SCALE;
    data.entry_from_near   = near_plane_cuts(model_to_clip, desc.clip_volume.min, desc.clip_volume.max) ? 1.0f : 0.0f;

    const SavedState saved = save_state();

    PUSH_GPU_SECTION("ISOSURFACES")
    // Empty space: the blocks that can hold one of these surfaces, drawn as boxes for the span of every ray
    bool use_proxy = false;
    if (count > 0) {
        if (const BlockGrid* g = block_grid_update(desc.texture.density_volume, dim)) {
            use_proxy = proxy_render(*g, model_to_clip, desc.clip_volume.min, desc.clip_volume.max, values, tau, count,
                                     (int)desc.render_target.width, (int)desc.render_target.height);
        }
    }
    data.use_proxy = use_proxy ? 1.0f : 0.0f;

    // The caller's framebuffer is the target when there is no colour texture
    if (!desc.render_target.color) {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, saved.fbo);
        glViewport(saved.viewport[0], saved.viewport[1], saved.viewport[2], saved.viewport[3]);
        if (saved.fbo == 0) {
            glDrawBuffer(saved.draw_buffer_count > 0 ? (GLenum)saved.draw_buffer[0] : GL_NONE);
        } else if (saved.draw_buffer_count > 0) {
            glDrawBuffers(saved.draw_buffer_count, (const GLenum*)saved.draw_buffer);
        }
    }

    const bool bound = bind_color_target(iso.fbo, desc.render_target.color, desc.render_target.width, desc.render_target.height, desc.render_target.clear_color);
    if (bound && count > 0) {
        glBindBuffer(GL_UNIFORM_BUFFER, iso.ubo);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(IsoUniformData), &data);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, iso.ubo);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_3D, desc.texture.density_volume);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, desc.render_target.depth);
        if (use_color_volume) {
            glActiveTexture(GL_TEXTURE2);
            glBindTexture(GL_TEXTURE_3D, desc.texture.color_volume);
        }
        if (use_field) {
            glActiveTexture(GL_TEXTURE3);
            glBindTexture(GL_TEXTURE_3D, desc.texture.field_volume);
            glActiveTexture(GL_TEXTURE4);
            glBindTexture(GL_TEXTURE_2D, desc.texture.field_colormap);
        }
        glActiveTexture(GL_TEXTURE5);
        glBindTexture(GL_TEXTURE_2D, use_proxy ? iso.proxy.tex_entry : 0);
        glActiveTexture(GL_TEXTURE6);
        glBindTexture(GL_TEXTURE_2D, use_proxy ? iso.proxy.tex_exit : 0);
        glActiveTexture(GL_TEXTURE0);

        glUseProgram(p.program);
        glUniformBlockBinding(p.program, p.block_index, 0);
        glUniform1i(p.loc_tex_volume, 0);
        glUniform1i(p.loc_tex_depth, 1);
        if (p.loc_tex_color_volume != -1)   glUniform1i(p.loc_tex_color_volume, 2);
        if (p.loc_tex_field != -1)          glUniform1i(p.loc_tex_field, 3);
        if (p.loc_tex_field_colormap != -1) glUniform1i(p.loc_tex_field_colormap, 4);
        glUniform1i(p.loc_tex_entry, 5);
        glUniform1i(p.loc_tex_exit,  6);
        glUniform1fv(p.loc_iso_values, count, values);
        glUniform4fv(p.loc_iso_colors, count, (const float*)colors);
        glUniform1fv(p.loc_iso_tau,    count, tau);
        glUniform1i (p.loc_iso_count,  count);

        // Premultiplied radiance over whatever is in the target
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
        glDisable(GL_SCISSOR_TEST);
        glDepthMask(GL_FALSE);
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

        glBindVertexArray(iso.vao);
        timer_begin(TimingStage_Raycast);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        timer_end();
        glBindVertexArray(0);
        glUseProgram(0);
    }
    POP_GPU_SECTION()

    restore_state(saved);
    return bound;
}

void initialize() {
    if (!timing.query[0].id) {
        for (int i = 0; i < TIMER_QUERY_COUNT; ++i) {
            glGenQueries(1, &timing.query[i].id);
        }
    }

    iso_initialize();

    glGetIntegerv(GL_MAJOR_VERSION, (GLint*)&gl.version.major);
    glGetIntegerv(GL_MINOR_VERSION, (GLint*)&gl.version.minor);

    GLuint v_shader_entry_exit          = gl::compile_shader_from_source({(const char*)entryexit_vert, entryexit_vert_size}, GL_VERTEX_SHADER);
    GLuint f_shader_entry_exit          = gl::compile_shader_from_source({(const char*)entryexit_frag, entryexit_frag_size}, GL_FRAGMENT_SHADER);
    GLuint f_shader_entry_exit_depth    = gl::compile_shader_from_source({(const char*)entryexit_frag, entryexit_frag_size}, GL_FRAGMENT_SHADER, STR_LIT("#define SAMPLE_DEPTH"));

    GLuint v_shader_vol                 = gl::compile_shader_from_source(v_shader_src_fs_quad, GL_VERTEX_SHADER);
    GLuint f_shader_dvr_only            = gl::compile_shader_from_source({(const char*)raycaster_frag, raycaster_frag_size}, GL_FRAGMENT_SHADER, STR_LIT("#define INCLUDE_DVR"));

    defer {
        glDeleteShader(v_shader_vol);
        glDeleteShader(v_shader_entry_exit);
        glDeleteShader(f_shader_entry_exit);
        glDeleteShader(f_shader_entry_exit_depth);
        glDeleteShader(f_shader_dvr_only);
    };

    if (v_shader_entry_exit == 0 || v_shader_vol == 0 || f_shader_entry_exit == 0 || f_shader_entry_exit_depth == 0 || f_shader_dvr_only == 0) {
        MD_LOG_ERROR("shader compilation failed, shader program for raycasting will not be updated");
        return;
    }
    
    if (gl.version.major >= 4 && gl.version.minor >= 3) {
        GLuint c_shader_splat_color = gl::compile_shader_from_source({ (const char*)splat_color_comp, splat_color_comp_size }, GL_COMPUTE_SHADER);
        if (c_shader_splat_color == 0) {
            MD_LOG_ERROR("shader compilation failed, shader program for splat color computation will not be updated");
            return;
        }
        if (!gl.program.splat_color) gl.program.splat_color = glCreateProgram();
        gl::attach_link_detach(gl.program.splat_color, &c_shader_splat_color, 1);
        glDeleteShader(c_shader_splat_color);
    }

    if (!gl.program.entry_exit) gl.program.entry_exit = glCreateProgram();
    if (!gl.program.entry_exit_depth) gl.program.entry_exit_depth = glCreateProgram();
    if (!gl.program.dvr_only) gl.program.dvr_only = glCreateProgram();

    {
        const GLuint shaders[] = {v_shader_entry_exit, f_shader_entry_exit};
        gl::attach_link_detach(gl.program.entry_exit, shaders, (int)ARRAY_SIZE(shaders));
    }
    {
        const GLuint shaders[] = {v_shader_entry_exit, f_shader_entry_exit_depth};
        gl::attach_link_detach(gl.program.entry_exit_depth, shaders, (int)ARRAY_SIZE(shaders));
    }
    {
        const GLuint shaders[] = {v_shader_vol, f_shader_dvr_only};
        gl::attach_link_detach(gl.program.dvr_only, shaders, (int)ARRAY_SIZE(shaders));
    }


    if (!gl.vbo) {
        // https://stackoverflow.com/questions/28375338/cube-using-single-gl-triangle-strip
        constexpr uint8_t cube_strip[42] = {0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1,
                                            0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0};
        glGenBuffers(1, &gl.vbo);
        glBindBuffer(GL_ARRAY_BUFFER, gl.vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(cube_strip), cube_strip, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    if (!gl.vao) {
        glGenVertexArrays(1, &gl.vao);
        glBindVertexArray(gl.vao);
        glBindBuffer(GL_ARRAY_BUFFER, gl.vbo);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_UNSIGNED_BYTE, GL_FALSE, 0, (const GLvoid*)0);
        glBindVertexArray(0);
    }

    if (!gl.ubo) {
        glGenBuffers(1, &gl.ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, gl.ubo);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(UniformData), 0, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
    }

    if (!gl.tex_entry) {
        glGenTextures(1, &gl.tex_entry);
    }

    if (!gl.tex_exit) {
        glGenTextures(1, &gl.tex_exit);
    }

    if (!gl.tex_result) {
        glGenTextures(1, &gl.tex_result);
    }

    if (!gl.fbo) {
        glGenFramebuffers(1, &gl.fbo);
    }

    if (!gl.ssbo) {
        glGenBuffers(1, &gl.ssbo);
    }
}

void shutdown() {}

mat4_t compute_model_to_world_matrix(vec3_t min_world_aabb, vec3_t max_world_aabb) {
    vec3_t ext = max_world_aabb - min_world_aabb;
    vec3_t off = min_world_aabb;
    return {
            ext.x, 0, 0, 0,
            0, ext.y, 0, 0,
            0, 0, ext.z, 0,
            off.x, off.y, off.z, 1,
        };
}

mat4_t compute_world_to_model_matrix(vec3_t min_world_aabb, vec3_t max_world_aabb) {
    vec3_t ext = max_world_aabb - min_world_aabb;
    vec3_t off = min_world_aabb;
    return {
            1.0f / ext.x, 0, 0, 0,
            0, 1.0f / ext.y, 0, 0,
            0, 0, 1.0f / ext.z, 0,
            -off.x, -off.y, -off.z, 1,
        };
}

mat4_t compute_model_to_texture_matrix(int dim_x, int dim_y, int dim_z) {
    (void)dim_x;
    (void)dim_y;
    (void)dim_z;
    return mat4_ident();
}

mat4_t compute_texture_to_model_matrix(int dim_x, int dim_y, int dim_z) {
    (void)dim_x;
    (void)dim_y;
    (void)dim_z;
    return mat4_ident();
}

void compute_transfer_function_texture_simple(uint32_t* tex, int colormap, float alpha_scale, int res) {
    ASSERT(tex);
    if (res <= 0) {
        MD_LOG_ERROR("Bad input resolution");
        return;
    }

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    size_t bytes = sizeof(uint32_t) * res;
    uint32_t* pixel_data = (uint32_t*)md_temp_alloc(temp, bytes);

    // Update colormap texture
    for (int i = 0; i < res; ++i) {
        float t = CLAMP((float)i / (float)(res - 1), 0.0f, 1.0f);
        ImVec4 col = ImPlot::SampleColormap(t, colormap);

        col.w = MIN(160 * t*t, 0.341176f) * alpha_scale;
        col.w = CLAMP(col.w, 0.0f, 1.0f);

        pixel_data[i] = ImGui::ColorConvertFloat4ToU32(col);
    }

    gl::init_texture_2D(tex, res, 1, GL_RGBA8);
    gl::set_texture_2D_data(*tex, 0, pixel_data, GL_RGBA8);

}

void compute_transfer_function_texture(uint32_t* tex, int colormap, ramp_type_t type, float ramp_scale, float ramp_period, int res) {
    ASSERT(tex);
    if (res <= 0) {
        MD_LOG_ERROR("Bad input resolution");
        return;
    }

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    size_t bytes = sizeof(uint32_t) * res;
    uint32_t* pixel_data = (uint32_t*)md_temp_alloc(temp, bytes);

    const float s = ramp_scale;
    const float p = ramp_period;

    // Update colormap texture
    for (int i = 0; i < res; ++i) {
        float t = (float)i / (float)(res - 1);
        ImVec4 col = ImPlot::SampleColormap(t, colormap);

        // Remap linear t which has a profile of '/' to '\/' to put the alpha ramp in the center (0.5)
        // y(x) = abs((1.0 - c) * (c - x)) Paste this into a grapher and you will see!alpha_origin;
        switch (type) {
        case RAMP_TYPE_SAWTOOTH:
            t = s * fmodf(t, p);
            break;
        case RAMP_TYPE_TRIANGLE:
            t = s * (2.0f/p) * fabsf(fmodf(t-p*0.5f,p)-p*0.5f);
            break;
        default:
            break;
        }

        t = CLAMP(t, 0.0f, 1.0f);
        col.w = t;

        pixel_data[i] = ImGui::ColorConvertFloat4ToU32(col);
    }

    gl::init_texture_2D(tex, res, 1, GL_RGBA8);
    gl::set_texture_2D_data(*tex, 0, pixel_data, GL_RGBA8);

}

static void splat_point_color_volume_GPU(uint32_t vol_texture, const int volume_dim[3], const float voxel_spacing[3], const float world_to_model[4][4], const float index_to_world[4][4], const vec4_t* point_xyzw, const uint32_t* point_color, size_t point_count, float power) {
    PUSH_GPU_SECTION("SPLAT COLOR VOLUME")
    glUseProgram(gl.program.splat_color);

    glUniformMatrix4fv(glGetUniformLocation(gl.program.splat_color, "u_world_to_model"), 1, GL_FALSE, (const float*)world_to_model);
    glUniformMatrix4fv(glGetUniformLocation(gl.program.splat_color, "u_voxel_to_world"), 1, GL_FALSE, (const float*)index_to_world);
    glUniform3iv(glGetUniformLocation(gl.program.splat_color, "u_volume_dim"),    1, volume_dim);
    glUniform3fv(glGetUniformLocation(gl.program.splat_color, "u_voxel_spacing"), 1, voxel_spacing);
    glUniform1i (glGetUniformLocation(gl.program.splat_color, "u_num_points"), (int)point_count);
    glUniform1f(glGetUniformLocation(gl.program.splat_color, "u_power"), power);

    size_t point_xyzw_offset    = 0;
	size_t point_xyzw_size      = ALIGN_TO(sizeof(vec4_t) * point_count, 256);
	size_t point_color_offset   = point_xyzw_offset + point_xyzw_size;
    size_t point_color_size     = ALIGN_TO(sizeof(uint32_t) * point_count, 256);
	size_t total_buffer_size    = ALIGN_TO(point_xyzw_size + point_color_size, 256);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, gl.ssbo);
    glBufferData(GL_SHADER_STORAGE_BUFFER, total_buffer_size , NULL, GL_DYNAMIC_DRAW);

    glBufferSubData(GL_SHADER_STORAGE_BUFFER, point_xyzw_offset,  sizeof(vec4_t)   * point_count, point_xyzw);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, point_color_offset, sizeof(uint32_t) * point_count, point_color);

    // Bind two contiguous vec4 arrays as separate shader storage buffer bindings for position/radius and color
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 0, gl.ssbo, point_xyzw_offset,  point_xyzw_size);  // point_xyzw
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 1, gl.ssbo, point_color_offset, point_color_size); // point_color
    glBindImageTexture(0, vol_texture, 0, GL_TRUE, 0, GL_WRITE_ONLY, GL_RGBA8);

    const int num_wg[3] = {
        DIV_UP((int)volume_dim[0], 8),
        DIV_UP((int)volume_dim[1], 8),
        DIV_UP((int)volume_dim[2], 8),
    };

    glDispatchCompute(num_wg[0], num_wg[1], num_wg[2]);
    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    POP_GPU_SECTION()
}

static void splat_point_color_volume_CPU(uint32_t vol_texture, const int volume_dim[3], const float voxel_spacing[3], const float world_to_model[4][4], const float index_to_world[4][4], const vec4_t* point_xyzw, const uint32_t* point_color, size_t point_count, float power) {
    ASSERT(volume_dim[0] % 4 == 0);
    ASSERT(volume_dim[1] % 4 == 0);
    ASSERT(volume_dim[2] % 4 == 0);

    (void)voxel_spacing;
    (void)world_to_model;

    md_temp_scope_t temp_scope = md_temp_begin();
    defer { md_temp_end(temp_scope); };

    size_t bytes = sizeof(vec4_t) * volume_dim[0] * volume_dim[1] * volume_dim[2];
    vec4_t* result = (vec4_t*)md_temp_alloc(temp_scope, bytes);
    ASSERT(result);

    mat4_t index_to_world_mat = mat4_load((const float*)index_to_world);

    for (int z = 0; z < volume_dim[2]; ++z) {
        for (int y = 0; y < volume_dim[1]; ++y) {
            for (int x = 0; x < volume_dim[0]; ++x) {
                vec4_t voxel_pos_index = {(float)x, (float)y, (float)z, 1.0f};
                vec4_t voxel_pos_world = mat4_mul_vec4(index_to_world_mat, voxel_pos_index);
                vec4_t acc = {0};
                float  sum = 0.0f;

                for (size_t i = 0; i < point_count; ++i) {
                    vec4_t point_world = point_xyzw[i];
                    float sigma = point_world.w;
                    point_world.w = 1.0f; // Ignore radius for distance calculation, we will use it as part of the influence factor instead
                    float dist2 = vec4_distance_squared(voxel_pos_world, point_world);

                    float inv2sig2 = 1.0f / (2.0f * sigma * sigma);
                    float w = expf(-dist2 * inv2sig2 * power);
                    acc += w * convert_color(point_color[i]);
                    sum += w;
                }

                vec4_t color = (sum > 0.0f) ? acc / sum : vec4_set(1.0f, 1.0f, 1.0f, 0.0f);
                acc /= sum;
                int linear_idx = x + y * volume_dim[0] + z * volume_dim[0] * volume_dim[1];
                result[linear_idx] = color;
            }
        }
    }

    // Write result into volume texture
    glBindTexture(GL_TEXTURE_3D, vol_texture);
    glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, volume_dim[0], volume_dim[1], volume_dim[2], GL_RGBA, GL_FLOAT, result);
    glBindTexture(GL_TEXTURE_3D, 0);

}

// vol_origin is the origin of the volume in world space
// voxel_spacing is the spacing of voxels in world space
void compute_point_color_volume(uint32_t vol_texture, const int volume_dim[3], const float voxel_spacing[3], const float world_to_model[4][4], const float index_to_world[4][4], const vec4_t* point_xyzw, const uint32_t* point_color, size_t point_count, double power) {
    if (!glIsTexture(vol_texture)) {
        MD_LOG_ERROR("Invalid volume texture");
        return;
    }

    int gl_major, gl_minor;
    glGetIntegerv(GL_MAJOR_VERSION, &gl_major);
    glGetIntegerv(GL_MINOR_VERSION, &gl_minor);

    if (gl_major > 4 || (gl_major == 4 && gl_minor >= 3)) {
        splat_point_color_volume_GPU(vol_texture, volume_dim, voxel_spacing, world_to_model, index_to_world, point_xyzw, point_color, point_count, (float)power);
    } else {
        splat_point_color_volume_CPU(vol_texture, volume_dim, voxel_spacing, world_to_model, index_to_world, point_xyzw, point_color, point_count, (float)power);
    }
}

// What the shared entry/exit + raycasting path takes: the union of the two public descriptions,
// with exactly one of dvr / iso enabled.
struct RaycastDesc {
    struct {
        uint32_t depth = 0;
        uint32_t color = 0;
        uint32_t width = 0;
        uint32_t height = 0;
        bool clear_color = false;
    } render_target;

    struct {
        uint32_t density_volume = 0;
        uint32_t color_volume   = 0;
        uint32_t transfer_function = 0;
        uint32_t field_volume   = 0;
        uint32_t field_colormap = 0;
    } texture;

    struct {
        mat4_t model = {};
        mat4_t view = {};
        mat4_t proj = {};
        mat4_t inv_proj = {};
    } matrix;

    struct {
        vec3_t min = {0, 0, 0};
        vec3_t max = {1, 1, 1};
    } clip_volume;

    struct {
        bool enabled = false;
    } temporal;

    struct {
        bool enabled = false;
        size_t count = 0;
        const float* values = NULL;
        const vec4_t* colors = NULL;
        const float* optical_densities = NULL;
        bool use_color_volume = false;
        bool use_field = false;
    } iso;

    struct {
        bool enabled = false;
        float min_tf_value = 0.0f;
        float max_tf_value = 1.0f;
    } dvr;

    struct {
        float range_beg = 0.0f;
        float range_end = 1.0f;
    } field;

    struct {
        vec3_t env_radiance = {0,0,0};
        float roughness = 0.4f;
        vec3_t dir_radiance = {1,1,1};
        float ior = 1.5f;
        float exposure = 1.0f;
        float gamma = 2.2f;
    } shading;

    vec3_t voxel_spacing = {};
};

static void render_raycast(const RaycastDesc& desc) {
    if (!desc.dvr.enabled && !desc.iso.enabled) return;

    int    iso_count = CLAMP((int)desc.iso.count, 0, 8);
    float  iso_values[8];
    vec4_t iso_colors[8];
    float  iso_optical_densities[8] = { 0 };

    MEMCPY(iso_values, desc.iso.values, iso_count * sizeof(float));
    MEMCPY(iso_colors, desc.iso.colors, iso_count * sizeof(vec4_t));
    if (desc.iso.optical_densities) {
        MEMCPY(iso_optical_densities, desc.iso.optical_densities, iso_count * sizeof(float));
    }

    // Sort on iso value
    for (int i = 0; i < iso_count - 1; ++i) {
        for (int j = i + 1; j < iso_count; ++j) {
            if (iso_values[j] < iso_values[i]) {
                float  val_tmp = iso_values[i];
                vec4_t col_tmp = iso_colors[i];
                float  od_tmp = iso_optical_densities[i];
                iso_values[i] = iso_values[j];
                iso_colors[i] = iso_colors[j];
                iso_optical_densities[i] = iso_optical_densities[j];
                iso_values[j] = val_tmp;
                iso_colors[j] = col_tmp;
                iso_optical_densities[j] = od_tmp;
            }
        }
    }

    // For the default framebuffer glDrawBuffers rejects the FRONT/BACK/LEFT/RIGHT
    // tokens that glDrawBuffer accepts and that GL_DRAW_BUFFER0 reports back, so the
    // restore has to go through glDrawBuffer there. See reset_gl_state() in
    // postprocessing_utils.cpp for the same fix.
    auto restore_draw_buffers = [](GLint fbo, const GLint* buffers, GLint count) {
        if (fbo == 0) {
            glDrawBuffer(count > 0 ? (GLenum)buffers[0] : GL_NONE);
        } else if (count > 0) {
            glDrawBuffers(count, (const GLenum*)buffers);
        } else {
            glDrawBuffer(GL_NONE);
        }
    };

    GLint bound_fbo;
    GLint bound_viewport[4];
    GLint bound_draw_buffer[8] = {0};
    GLint bound_draw_buffer_count = 0;
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &bound_fbo);
    glGetIntegerv(GL_VIEWPORT, bound_viewport);
    for (int i = 0; i < 8; ++i) {
        glGetIntegerv(GL_DRAW_BUFFER0 + i, &bound_draw_buffer[i]);
        // @NOTE: Assume that its tightly packed and if we stumple upon a zero draw buffer index, we enterpret that as the 'end'
        if (bound_draw_buffer[i] != GL_NONE) {
            bound_draw_buffer_count = i + 1;
        }
    }

    if (gl.width < desc.render_target.width ||
        gl.height < desc.render_target.height)
    {
        gl.width = desc.render_target.width;
        gl.height = desc.render_target.height;
        gl::init_texture_2D(&gl.tex_entry,  gl.width, gl.height, GL_RGB16);
        gl::init_texture_2D(&gl.tex_exit,   gl.width, gl.height, GL_RGB16);
        gl::init_texture_2D(&gl.tex_result, gl.width, gl.height, GL_RGBA8);
    }

    const mat4_t model_to_view_matrix = mat4_mul(desc.matrix.view, desc.matrix.model);

    static float time = 0.0f;
    time += 1.0f / 100.0f;
    if (time > 100.0) time -= 100.0f;
    if (!desc.temporal.enabled) {
        time = 0.0f;
    }

    float tf_min = desc.dvr.min_tf_value;
    float tf_max = desc.dvr.max_tf_value;
    float tf_ext = tf_max - tf_min;
    float inv_tf_ext = tf_ext == 0 ? 1.0f : 1.0f / tf_ext;

    const float n1 = 1.0f;
    const float n2 = desc.shading.ior;
    const float F0 = powf((n1-n2)/(n1+n2), 2.0f);

    UniformData data;
    data.view_to_model_mat = mat4_inverse(model_to_view_matrix);
    data.model_to_view_mat = model_to_view_matrix;
    data.inv_proj_mat      = desc.matrix.inv_proj;
    data.model_view_proj_mat = desc.matrix.proj * model_to_view_matrix;
    data.inv_res = {1.f / (float)(desc.render_target.width), 1.f / (float)(desc.render_target.height)};
    data.time = time;
    data.gamma = desc.shading.gamma;
    data.clip_volume_min = desc.clip_volume.min;
    data.tf_min = tf_min;
    data.clip_volume_max = desc.clip_volume.max;
    data.tf_inv_ext = inv_tf_ext;
    data.gradient_spacing_world_space = desc.voxel_spacing;
    data.exposure = desc.shading.exposure;
    data.gradient_spacing_tex_space = data.view_to_model_mat * mat4_scale(desc.voxel_spacing.x, desc.voxel_spacing.y, desc.voxel_spacing.z);
    data.env_radiance = desc.shading.env_radiance;
    data.roughness = desc.shading.roughness;
    data.dir_radiance = desc.shading.dir_radiance;
    data.F0 = F0;

    glBindBuffer(GL_UNIFORM_BUFFER, gl.ubo);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(UniformData), &data);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    bool use_depth = desc.render_target.depth;

    if (use_depth) {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, desc.render_target.depth);
    }

    glBindBufferBase(GL_UNIFORM_BUFFER, 0, gl.ubo);
    glBindVertexArray(gl.vao);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.tex_entry, 0);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gl.tex_exit,  0);
    const GLuint draw_bufs[] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
    glDrawBuffers(2, draw_bufs);

    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);

    glViewport(0, 0, desc.render_target.width, desc.render_target.height);

    glEnable(GL_CULL_FACE);
    {
        PUSH_GPU_SECTION("VOLUME ENTRY / EXIT");
        
        const GLuint prog = use_depth ? gl.program.entry_exit_depth : gl.program.entry_exit;
        const GLint uniform_block_index = glGetUniformBlockIndex(prog, "UniformData");
        const GLint uniform_loc_tex_depth = glGetUniformLocation(prog, "u_tex_depth");

        glUseProgram(prog);
        glUniform1i(uniform_loc_tex_depth, 0);
        glUniformBlockBinding(prog, uniform_block_index, 0);

        glCullFace(GL_FRONT);
        timer_begin(TimingStage_EntryExit);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 42);
        timer_end();
        POP_GPU_SECTION()
    }

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, gl.tex_entry);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, gl.tex_exit);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_3D, desc.texture.density_volume);

    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, desc.texture.transfer_function);

    // A field takes the place of the colour volume: the two are alternatives, never both
    const bool use_field        = desc.iso.enabled && desc.iso.use_field && desc.texture.field_volume && desc.texture.field_colormap;
    const bool use_color_volume = desc.iso.enabled && !use_field && desc.iso.use_color_volume;

    if (use_color_volume) {
        glActiveTexture(GL_TEXTURE4);
        glBindTexture(GL_TEXTURE_3D, desc.texture.color_volume);
    }
    if (use_field) {
        glActiveTexture(GL_TEXTURE5);
        glBindTexture(GL_TEXTURE_3D, desc.texture.field_volume);
        glActiveTexture(GL_TEXTURE6);
        glBindTexture(GL_TEXTURE_2D, desc.texture.field_colormap);
    }

    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.tex_result, 0);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    if (desc.render_target.color) {
        ASSERT(glIsTexture(desc.render_target.color));
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, desc.render_target.color, 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        GLenum vol_status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
        if (vol_status != GL_FRAMEBUFFER_COMPLETE) {
            MD_LOG_ERROR("Volume render target framebuffer is incomplete (0x%04X)", (unsigned int)vol_status);
        }
        if (desc.render_target.clear_color) {
            glClearColor(0, 0, 0, 0);
            glClear(GL_COLOR_BUFFER_BIT);
        }
    } else {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, bound_fbo);
        glViewport(bound_viewport[0], bound_viewport[1], bound_viewport[2], bound_viewport[3]);
        restore_draw_buffers(bound_fbo, bound_draw_buffer, bound_draw_buffer_count);
    }

    PUSH_GPU_SECTION("VOLUME RAYCASTING")
    {
        GLuint vol_prog = desc.dvr.enabled ? gl.program.dvr_only : 0;

        if (vol_prog == 0) {
            MD_LOG_DEBUG("No raycasting shader program available for the current render description, skipping raycasting");
            glBindVertexArray(0);
            glDisable(GL_BLEND);
            glEnable(GL_DEPTH_TEST);
            glCullFace(GL_BACK);
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, bound_fbo);
            glViewport(bound_viewport[0], bound_viewport[1], bound_viewport[2], bound_viewport[3]);
            restore_draw_buffers(bound_fbo, bound_draw_buffer, bound_draw_buffer_count);
            POP_GPU_SECTION()
            return;
        }

        const GLint uniform_block_index             = glGetUniformBlockIndex(vol_prog, "UniformData");
        const GLint uniform_loc_tex_entry           = glGetUniformLocation(vol_prog, "u_tex_entry");
        const GLint uniform_loc_tex_exit            = glGetUniformLocation(vol_prog, "u_tex_exit");
        const GLint uniform_loc_tex_tf              = glGetUniformLocation(vol_prog, "u_tex_tf");
        const GLint uniform_loc_tex_density_volume  = glGetUniformLocation(vol_prog, "u_tex_density_volume");
        const GLint uniform_loc_tex_color_volume    = glGetUniformLocation(vol_prog, "u_tex_color_volume");
        const GLint uniform_loc_iso_values          = glGetUniformLocation(vol_prog, "u_iso.values");
        const GLint uniform_loc_iso_colors          = glGetUniformLocation(vol_prog, "u_iso.colors");
        const GLint uniform_loc_iso_optical_densities = glGetUniformLocation(vol_prog, "u_iso.optical_densities");
        const GLint uniform_loc_iso_count           = glGetUniformLocation(vol_prog, "u_iso.count");

        glUseProgram(vol_prog);

        glUniform1i(uniform_loc_tex_entry, 0);
        glUniform1i(uniform_loc_tex_exit,  1);
        glUniform1i(uniform_loc_tex_density_volume, 2);
        glUniform1i(uniform_loc_tex_tf, 3);
        glUniform1i(uniform_loc_tex_color_volume, 4);
        if (use_field) {
            const float ext = desc.field.range_end - desc.field.range_beg;
            glUniform1i(glGetUniformLocation(vol_prog, "u_tex_field"), 5);
            glUniform1i(glGetUniformLocation(vol_prog, "u_tex_field_colormap"), 6);
            glUniform2f(glGetUniformLocation(vol_prog, "u_field_range"), desc.field.range_beg, ext != 0.0f ? 1.0f / ext : 0.0f);
        }
        glUniform1fv(uniform_loc_iso_values, (GLsizei)iso_count, (const float*)iso_values);
        glUniform4fv(uniform_loc_iso_colors, (GLsizei)iso_count, (const float*)iso_colors);
        glUniform1fv(uniform_loc_iso_optical_densities, (GLsizei)iso_count, (const float*)iso_optical_densities);
        glUniform1i(uniform_loc_iso_count, (int)iso_count);
        glUniformBlockBinding(vol_prog, uniform_block_index, 0);

        timer_begin(TimingStage_Raycast);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        timer_end();

        glBindVertexArray(0);
        glUseProgram(0);
    }
    POP_GPU_SECTION()

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glCullFace(GL_BACK);

    if (desc.render_target.color) {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, bound_fbo);
        glViewport(bound_viewport[0], bound_viewport[1], bound_viewport[2], bound_viewport[3]);
        restore_draw_buffers(bound_fbo, bound_draw_buffer, bound_draw_buffer_count);
    }

}

void render_dvr(const DvrRenderDesc& desc) {
    if (!desc.texture.density_volume || !desc.texture.transfer_function) return;

    RaycastDesc rc = {};
    rc.render_target.depth       = desc.render_target.depth;
    rc.render_target.color       = desc.render_target.color;
    rc.render_target.width       = desc.render_target.width;
    rc.render_target.height      = desc.render_target.height;
    rc.render_target.clear_color = desc.render_target.clear_color;
    rc.texture.density_volume    = desc.texture.density_volume;
    rc.texture.transfer_function = desc.texture.transfer_function;
    rc.matrix.model    = desc.matrix.model;
    rc.matrix.view     = desc.matrix.view;
    rc.matrix.proj     = desc.matrix.proj;
    rc.matrix.inv_proj = desc.matrix.inv_proj;
    rc.clip_volume.min = desc.clip_volume.min;
    rc.clip_volume.max = desc.clip_volume.max;
    rc.temporal.enabled = desc.temporal.enabled;
    rc.dvr.enabled = true;
    rc.dvr.min_tf_value = desc.tf.min_value;
    rc.dvr.max_tf_value = desc.tf.max_value;
    rc.voxel_spacing = desc.voxel_spacing;
    render_raycast(rc);
}

}  // namespace volume
