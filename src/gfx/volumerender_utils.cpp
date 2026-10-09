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
    GLuint vao = 0;     // empty: the full screen triangle comes from gl_VertexID
    GLuint ubo = 0;
    GLuint fbo = 0;
    GLuint ssbo = 0;

    struct {
        GLuint dvr = 0;
        GLuint splat_color = 0;
    } program;

    struct {
        GLint tex_volume = -1;
        GLint tex_depth = -1;
        GLint tex_tf = -1;
        GLint block_index = -1;
    } dvr_loc;

    struct {
        int major = 0;
        int minor = 0;
    } version;
} gl;

// std140, mirrors DvrUniforms in dvr.frag
struct DvrUniformData {
    mat4_t clip_to_model;
    vec3_t clip_min;
    float  tf_min;
    vec3_t clip_max;
    float  tf_inv_ext;
    vec2_t inv_res;
    float  time;
    float  use_depth;
};
static_assert(sizeof(DvrUniformData) == 64 + 3 * 16, "DvrUniformData must match the std140 layout of DvrUniforms");

// -----------------------------------------------------------------------------
// GPU timings
// -----------------------------------------------------------------------------

// A redraw of the orbital grid alone issues three queries per panel, and they come back a few frames late
static constexpr int    TIMER_QUERY_COUNT     = 256;
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

enum IsoVariant {
    IsoVariant_Uniform,
    IsoVariant_ColorVolume,
    IsoVariant_Field,
    IsoVariant_Count
};

static const char* iso_variant_defines[IsoVariant_Count] = {
    "",
    "#define USE_COLOR_VOLUME",
    "#define USE_FIELD",
};

// How the rays find the surfaces: isosurface_fast.frag samples once per voxel, isosurface.frag intersects
// every cell exactly. Both read the same IsoUniforms.
enum IsoMode {
    IsoMode_Fast,
    IsoMode_Exact,
    IsoMode_Count
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
    GLint  loc_tex_minmax = -1;
    GLint  block_index = -1;
};

// std140, mirrors IsoUniforms in isosurface.frag
struct IsoUniformData {
    mat4_t clip_to_model;
    mat4_t model_to_view;
    mat4_t grad_offsets;

    vec3_t clip_min;
    float  block_size;      // exact: voxels per block side of the block grid; fast: samples per voxel
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
    float  use_blocks;
    float  pad;
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
    IsoProgram prog[IsoMode_Count][IsoVariant_Count];          // compiled on first use
    bool       prog_failed[IsoMode_Count][IsoVariant_Count];
    GLuint ubo = 0;
    GLuint fbo = 0;
    GLuint vao = 0;     // empty: full screen triangle and proxy boxes both come from gl_VertexID

    BlockGrid grid[BLOCK_CACHE_SIZE];
    uint64_t  use_counter = 0;
    uint64_t  version_counter = 0;  // versions are unique across volumes and evictions, so a (texture, version) pair never repeats

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

static void iso_program_setup(IsoProgram* p, GLuint v_shader, str_t source, str_t defines) {
    GLuint f_shader = gl::compile_shader_from_source(source, GL_FRAGMENT_SHADER, defines);
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
    p->loc_tex_minmax         = glGetUniformLocation(prog, "u_tex_minmax");
    p->block_index            = glGetUniformBlockIndex(prog, "IsoUniforms");
}

static void iso_program_compile(int mode, int variant, GLuint v_shader) {
    const str_t source = (mode == IsoMode_Exact) ? str_t{(const char*)isosurface_frag,      isosurface_frag_size}
                                                 : str_t{(const char*)isosurface_fast_frag, isosurface_fast_frag_size};
    iso_program_setup(&iso.prog[mode][variant], v_shader, source, str_from_cstr(iso_variant_defines[variant]));
    iso.prog_failed[mode][variant] = iso.prog[mode][variant].program == 0;
}

// The program of a mode and variant, compiled the first time it is asked for; NULL if it does not compile
static const IsoProgram* iso_program_get(int mode, int variant) {
    IsoProgram& p = iso.prog[mode][variant];
    if (!p.program && !iso.prog_failed[mode][variant]) {
        GLuint v_shader = gl::compile_shader_from_source(v_shader_src_fs_quad, GL_VERTEX_SHADER);
        if (v_shader) {
            iso_program_compile(mode, variant, v_shader);
            glDeleteShader(v_shader);
        } else {
            iso.prog_failed[mode][variant] = true;
        }
    }
    return p.program ? &p : NULL;
}

static void iso_initialize() {
    GLuint v_shader = gl::compile_shader_from_source(v_shader_src_fs_quad, GL_VERTEX_SHADER);
    if (!v_shader) {
        MD_LOG_ERROR("Isosurface vertex shader compilation failed");
        return;
    }
    // Programs already made are rebuilt (a reinitialize picks up shader changes); the others wait until used
    for (int m = 0; m < IsoMode_Count; ++m) {
        for (int v = 0; v < IsoVariant_Count; ++v) {
            iso.prog_failed[m][v] = false;
            if (iso.prog[m][v].program) {
                iso_program_compile(m, v, v_shader);
            }
        }
    }

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
    slot->version   = ++iso.version_counter;
    slot->last_used = ++iso.use_counter;
    return slot;
}

void notify_data_changed(uint32_t volume_texture) {
    if (!volume_texture) return;
    block_grid_acquire(volume_texture)->version = ++iso.version_counter;
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
    // Boxes that reach past the far plane (or between the eye and the near plane) are not clipped there
    // but clamped to it, so the farthest depth is the far plane rather than a front face inside the box
    const GLboolean depth_clamp = glIsEnabled(GL_DEPTH_CLAMP);
    glEnable(GL_DEPTH_CLAMP);

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
    if (!depth_clamp) glDisable(GL_DEPTH_CLAMP);
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
    const int variant = use_field ? IsoVariant_Field : use_color_volume ? IsoVariant_ColorVolume : IsoVariant_Uniform;
    int mode = desc.iso.exact ? IsoMode_Exact : IsoMode_Fast;
    const IsoProgram* pp = iso_program_get(mode, variant);
    if (!pp) {
        // The other mode rather than nothing
        mode = (mode == IsoMode_Exact) ? IsoMode_Fast : IsoMode_Exact;
        pp = iso_program_get(mode, variant);
    }

    int dim[3] = {};
    if (!pp || !gl::get_texture_dim(dim, desc.texture.density_volume) || dim[0] <= 0 || dim[1] <= 0 || dim[2] <= 0) {
        return false;
    }
    const IsoProgram& p = *pp;

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
    data.block_size        = 1.0f;
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
    // Empty space: the blocks that can hold one of these surfaces, drawn as boxes for the span of every ray,
    // and walked by the raycaster to cross the others in one step
    bool use_proxy = false;
    const BlockGrid* grid = count > 0 ? block_grid_update(desc.texture.density_volume, dim) : NULL;
    if (grid) {
        use_proxy = proxy_render(*grid, model_to_clip, desc.clip_volume.min, desc.clip_volume.max, values, tau, count,
                                 (int)desc.render_target.width, (int)desc.render_target.height);
        data.block_size = (float)grid->block_size;
    }
    data.use_proxy  = use_proxy ? 1.0f : 0.0f;
    data.use_blocks = grid ? 1.0f : 0.0f;
    if (mode == IsoMode_Fast) {
        data.block_size = 1.0f;     // samples per voxel
    }

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
        glActiveTexture(GL_TEXTURE7);
        glBindTexture(GL_TEXTURE_3D, grid ? grid->minmax : 0);
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
        glUniform1i(p.loc_tex_minmax, 7);
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

        // The proxy depths are rendered into again by the next call: not left bound for sampling
        glActiveTexture(GL_TEXTURE5);
        glBindTexture(GL_TEXTURE_2D, 0);
        glActiveTexture(GL_TEXTURE6);
        glBindTexture(GL_TEXTURE_2D, 0);
        glActiveTexture(GL_TEXTURE7);
        glBindTexture(GL_TEXTURE_3D, 0);
        glActiveTexture(GL_TEXTURE0);
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

    {
        GLuint v_shader = gl::compile_shader_from_source(v_shader_src_fs_quad, GL_VERTEX_SHADER);
        GLuint f_shader = gl::compile_shader_from_source({(const char*)dvr_frag, dvr_frag_size}, GL_FRAGMENT_SHADER);
        if (v_shader && f_shader) {
            if (!gl.program.dvr) gl.program.dvr = glCreateProgram();
            const GLuint shaders[] = {v_shader, f_shader};
            gl::attach_link_detach(gl.program.dvr, shaders, (int)ARRAY_SIZE(shaders));
            gl.dvr_loc.tex_volume  = glGetUniformLocation(gl.program.dvr, "u_tex_volume");
            gl.dvr_loc.tex_depth   = glGetUniformLocation(gl.program.dvr, "u_tex_depth");
            gl.dvr_loc.tex_tf      = glGetUniformLocation(gl.program.dvr, "u_tex_tf");
            gl.dvr_loc.block_index = glGetUniformBlockIndex(gl.program.dvr, "DvrUniforms");
        } else {
            MD_LOG_ERROR("DVR shader compilation failed, the DVR program will not be updated");
        }
        if (v_shader) glDeleteShader(v_shader);
        if (f_shader) glDeleteShader(f_shader);
    }

    if (gl.version.major > 4 || (gl.version.major == 4 && gl.version.minor >= 3)) {
        GLuint c_shader_splat_color = gl::compile_shader_from_source({ (const char*)splat_color_comp, splat_color_comp_size }, GL_COMPUTE_SHADER);
        if (c_shader_splat_color) {
            if (!gl.program.splat_color) gl.program.splat_color = glCreateProgram();
            gl::attach_link_detach(gl.program.splat_color, &c_shader_splat_color, 1);
            glDeleteShader(c_shader_splat_color);
        } else {
            MD_LOG_ERROR("shader compilation failed, shader program for splat color computation will not be updated");
        }
    }

    if (!gl.vao) {
        glGenVertexArrays(1, &gl.vao);
    }

    if (!gl.ubo) {
        glGenBuffers(1, &gl.ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, gl.ubo);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(DvrUniformData), 0, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
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

static void splat_point_color_volume_GPU(uint32_t vol_texture, const int volume_dim[3], const float index_to_world[4][4], const vec4_t* point_xyzw, const uint32_t* point_color, size_t point_count, float power) {
    PUSH_GPU_SECTION("SPLAT COLOR VOLUME")
    glUseProgram(gl.program.splat_color);

    glUniformMatrix4fv(glGetUniformLocation(gl.program.splat_color, "u_voxel_to_world"), 1, GL_FALSE, (const float*)index_to_world);
    glUniform3iv(glGetUniformLocation(gl.program.splat_color, "u_volume_dim"), 1, volume_dim);
    glUniform1i (glGetUniformLocation(gl.program.splat_color, "u_num_points"), (int)point_count);
    glUniform1f (glGetUniformLocation(gl.program.splat_color, "u_power"), power);

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

// The CPU twin of splat_color.comp (used where there are no compute shaders, GL < 4.3): normalised
// Gaussian weights evaluated relative to the largest one, over the points that can matter in each 8^3
// block of voxels.
static void splat_point_color_volume_CPU(uint32_t vol_texture, const int volume_dim[3], const float index_to_world[4][4], const vec4_t* point_xyzw, const uint32_t* point_color, size_t point_count, float power) {
    const int dim_x = volume_dim[0];
    const int dim_y = volume_dim[1];
    const int dim_z = volume_dim[2];
    if (dim_x <= 0 || dim_y <= 0 || dim_z <= 0 || point_count == 0) return;

    static constexpr int   BLOCK  = 8;
    static constexpr float CUTOFF = 20.0f;  // a point weighing less than exp(-CUTOFF) of the best one is left out

    md_temp_scope_t temp_scope = md_temp_begin();
    defer { md_temp_end(temp_scope); };

    uint32_t* result = (uint32_t*)md_temp_alloc(temp_scope, sizeof(uint32_t) * (size_t)dim_x * dim_y * dim_z);
    vec4_t*   points = (vec4_t*)  md_temp_alloc(temp_scope, sizeof(vec4_t) * point_count);    // xyz, power / (2 sigma^2)
    vec4_t*   colors = (vec4_t*)  md_temp_alloc(temp_scope, sizeof(vec4_t) * point_count);
    uint32_t* cand   = (uint32_t*)md_temp_alloc(temp_scope, sizeof(uint32_t) * point_count);
    if (!result || !points || !colors || !cand) return;

    for (size_t i = 0; i < point_count; ++i) {
        const float sigma = MAX(point_xyzw[i].w, 0.1f);
        points[i] = vec4_set(point_xyzw[i].x, point_xyzw[i].y, point_xyzw[i].z, power / (2.0f * sigma * sigma));
        colors[i] = convert_color(point_color[i]);
    }

    const mat4_t M = mat4_load((const float*)index_to_world);

    for (int bz = 0; bz < dim_z; bz += BLOCK) {
        for (int by = 0; by < dim_y; by += BLOCK) {
            for (int bx = 0; bx < dim_x; bx += BLOCK) {
                const int lo[3] = { bx, by, bz };
                const int hi[3] = { MIN(bx + BLOCK, dim_x) - 1, MIN(by + BLOCK, dim_y) - 1, MIN(bz + BLOCK, dim_z) - 1 };

                // World space bounds of the voxel centres of the block (the volume may be rotated)
                vec3_t bmin = vec3_set( FLT_MAX,  FLT_MAX,  FLT_MAX);
                vec3_t bmax = vec3_set(-FLT_MAX, -FLT_MAX, -FLT_MAX);
                for (int c = 0; c < 8; ++c) {
                    const vec4_t idx = { (float)((c & 1) ? hi[0] : lo[0]), (float)((c & 2) ? hi[1] : lo[1]), (float)((c & 4) ? hi[2] : lo[2]), 1.0f };
                    const vec3_t p = vec3_from_vec4(mat4_mul_vec4(M, idx));
                    bmin = vec3_min(bmin, p);
                    bmax = vec3_max(bmax, p);
                }

                // The best exponent some point is guaranteed everywhere in the block, then the points that can matter
                float e_best = -FLT_MAX;
                for (size_t i = 0; i < point_count; ++i) {
                    const vec3_t p   = vec3_from_vec4(points[i]);
                    const vec3_t far = vec3_max(vec3_abs(vec3_sub(p, bmin)), vec3_abs(vec3_sub(p, bmax)));
                    e_best = MAX(e_best, -vec3_dot(far, far) * points[i].w);
                }
                uint32_t num_cand = 0;
                for (size_t i = 0; i < point_count; ++i) {
                    const vec3_t p    = vec3_from_vec4(points[i]);
                    const vec3_t near = vec3_sub(p, vec3_clamp(p, bmin, bmax));
                    if (-vec3_dot(near, near) * points[i].w >= e_best - CUTOFF) {
                        cand[num_cand++] = (uint32_t)i;
                    }
                }

                for (int z = lo[2]; z <= hi[2]; ++z) {
                    for (int y = lo[1]; y <= hi[1]; ++y) {
                        for (int x = lo[0]; x <= hi[0]; ++x) {
                            const vec3_t xw = vec3_from_vec4(mat4_mul_vec4(M, vec4_set((float)x, (float)y, (float)z, 1.0f)));
                            // Streamed log-sum-exp: weights relative to the largest exponent so far
                            float  m   = -FLT_MAX;
                            float  sum = 0.0f;
                            vec4_t acc = {0, 0, 0, 0};
                            for (uint32_t j = 0; j < num_cand; ++j) {
                                const uint32_t i = cand[j];
                                const vec3_t d = vec3_sub(xw, vec3_from_vec4(points[i]));
                                const float  e = -vec3_dot(d, d) * points[i].w;
                                if (e > m) {
                                    const float scale = expf(m - e);
                                    sum = sum * scale + 1.0f;
                                    acc = vec4_add(vec4_mul1(acc, scale), colors[i]);
                                    m = e;
                                } else {
                                    const float w = expf(e - m);
                                    sum += w;
                                    acc = vec4_add(acc, vec4_mul1(colors[i], w));
                                }
                            }
                            const vec4_t color = (sum > 0.0f) ? vec4_mul1(acc, 1.0f / sum) : vec4_set(1.0f, 1.0f, 1.0f, 0.0f);
                            result[(size_t)x + (size_t)y * dim_x + (size_t)z * dim_x * dim_y] = convert_color(color);
                        }
                    }
                }
            }
        }
    }

    glBindTexture(GL_TEXTURE_3D, vol_texture);
    glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, dim_x, dim_y, dim_z, GL_RGBA, GL_UNSIGNED_BYTE, result);
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

    (void)voxel_spacing;
    (void)world_to_model;
    if ((gl_major > 4 || (gl_major == 4 && gl_minor >= 3)) && gl.program.splat_color) {
        splat_point_color_volume_GPU(vol_texture, volume_dim, index_to_world, point_xyzw, point_color, point_count, (float)power);
    } else {
        splat_point_color_volume_CPU(vol_texture, volume_dim, index_to_world, point_xyzw, point_color, point_count, (float)power);
    }
}

// -----------------------------------------------------------------------------
// Direct volume rendering
// -----------------------------------------------------------------------------

void render_dvr(const DvrRenderDesc& desc) {
    if (!desc.texture.density_volume || !desc.texture.transfer_function || !gl.program.dvr ||
        desc.render_target.width == 0 || desc.render_target.height == 0) {
        return;
    }

    // Moves the jitter pattern every frame when something (TAA) integrates it, keeps it still otherwise
    static float time = 0.0f;
    time = desc.temporal.enabled ? fmodf(time + 0.01f, 100.0f) : 0.0f;

    const float tf_ext = desc.tf.max_value - desc.tf.min_value;
    const mat4_t model_to_clip = desc.matrix.proj * desc.matrix.view * desc.matrix.model;

    DvrUniformData data = {};
    data.clip_to_model = mat4_inverse(model_to_clip);
    data.clip_min   = desc.clip_volume.min;
    data.tf_min     = desc.tf.min_value;
    data.clip_max   = desc.clip_volume.max;
    data.tf_inv_ext = tf_ext != 0.0f ? 1.0f / tf_ext : 1.0f;
    data.inv_res    = {1.0f / (float)desc.render_target.width, 1.0f / (float)desc.render_target.height};
    data.time       = time;
    data.use_depth  = desc.render_target.depth ? 1.0f : 0.0f;

    const SavedState saved = save_state();

    PUSH_GPU_SECTION("DVR")
    if (bind_color_target(gl.fbo, desc.render_target.color, desc.render_target.width, desc.render_target.height, desc.render_target.clear_color)) {
        glBindBuffer(GL_UNIFORM_BUFFER, gl.ubo);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(DvrUniformData), &data);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, gl.ubo);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_3D, desc.texture.density_volume);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, desc.render_target.depth);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, desc.texture.transfer_function);
        glActiveTexture(GL_TEXTURE0);

        glUseProgram(gl.program.dvr);
        glUniformBlockBinding(gl.program.dvr, gl.dvr_loc.block_index, 0);
        glUniform1i(gl.dvr_loc.tex_volume, 0);
        glUniform1i(gl.dvr_loc.tex_depth, 1);
        glUniform1i(gl.dvr_loc.tex_tf, 2);

        // Premultiplied colour over whatever is in the target
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
        glDisable(GL_SCISSOR_TEST);
        glDepthMask(GL_FALSE);
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

        glBindVertexArray(gl.vao);
        timer_begin(TimingStage_Raycast);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        timer_end();
        glBindVertexArray(0);
        glUseProgram(0);
    }
    POP_GPU_SECTION()

    restore_state(saved);
}

}  // namespace volume
