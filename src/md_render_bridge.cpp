#include "md_render_bridge.h"

#include <viamd.h>
#include <representation.h>
#include <color_utils.h>

#include <core/md_log.h>
#include <core/md_hash.h>
#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_common.h>
#include <core/md_arena_allocator.h>

#include <gfx/gl.h>
#include <gfx/gl_utils.h>
#include <imgui.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if MD_ENABLE_GPU

#include <md_render.h>
#include <md_render_mol.h>
#include <md_compose.h>
#include <md_post.h>
#include <core/md_gpu_gl.h>
#include <gfx/postprocessing_utils.h>

// Pairs of atoms drawn as sticks (md_render_mol_rep_t.pairs): by kind, each with its id (what a pick of it reports)
struct PairSet {
    md_gpu_addr_t pairs = 0, ids = 0;
    uint32_t      cap = 0;
    uint32_t      offset[MD_RENDER_MOL_BOND_KIND_COUNT + 1] = {};
};

static uint32_t pair_count(const PairSet& p) { return p.offset[MD_RENDER_MOL_BOND_KIND_COUNT]; }

static void free_pairs(md_gpu_stream_t s, PairSet* p) {
    if (p->pairs) md_gpu_free(s, p->pairs);
    if (p->ids)   md_gpu_free(s, p->ids);
    *p = {};
}

// A representation's arrays on the GPU: its atoms' colours and which it shows, as md_gl has them (rep_colors), and the
// profile of its cartoon
struct BridgeRep {
    uint32_t      id = 0;               // Representation::node.id
    md_gpu_addr_t color = 0;            // uint32 RGBA8 per atom
    md_gpu_addr_t visible = 0;          // bitset over the atoms, uint32 words
    uint64_t      visible_hash = 0;     // of what it holds
    uint32_t      atom_count = 0;
    bool          seen = false;         // still a representation (render)
    md_gpu_addr_t profile = 0;          // Cartoon: md_render_mol_cartoon_profile's, of profile_key
    uint64_t      profile_key = 0;
    PairSet       pairs;                // Bonds: its pairs, each one's id its index in the representation's bond_pairs
};

static void free_rep(md_gpu_stream_t s, BridgeRep* r) {
    if (r->color)   md_gpu_free(s, r->color);
    if (r->visible) md_gpu_free(s, r->visible);
    if (r->profile) md_gpu_free(s, r->profile);
    free_pairs(s, &r->pairs);
}

// A texture md_gpu writes and GL reads, and the GL texture of it. The two take turns with it: md_gpu writes it between
// acquire and release, GL reads it between begin_read and end_read, and md_gpu writes it again only after GL said so
// (gl_done)
struct SharedTex {
    md_gpu_texture_t    tex = nullptr;
    md_gpu_gl_texture_t gl  = {};
    bool                imported = false;
    bool                released = false;   // released by md_gpu, not yet read by GL
};

struct MdRenderBridge {
    bool        ok = false;
    char        status[256] = "";

    md_gpu_device_t      device = nullptr;
    md_gpu_stream_t      stream = nullptr;     // graphics
    md_render_ctx_t      render = nullptr;
    md_render_mol_ctx_t  mol    = nullptr;
    md_render_view_t     view   = nullptr;

    // The system: its topology, the per-frame block drawn, the positions, the marks (selection, highlight, its groups)
    bool                 has_topo = false;
    md_render_mol_topo_t topo = {};
    md_gpu_addr_t        dynamic = 0;
    md_gpu_addr_t        positions = 0;
    md_gpu_addr_t        marks = 0;            // uint8 per atom (md_render_marks_write)
    uint32_t             atom_count = 0;
    uint32_t             pending = 0;          // MolBit_ flags not yet taken up
    bool                 updated = false;      // the per-frame block holds the positions
    int64_t              regroup_frame = -1;   // the trajectory frame of the last regroup
    double               time = 0.0;           // md_render's: one per update of the positions, the velocities' unit
    // The bonds of the frame shown while they are tracked (md_bond_delta_t, md_chem_delta_t): md_render_mol's block
    // of them (md_render_mol_bonds_write), the system's that are not the frame's hidden, the frame's own drawn
    md_gpu_addr_t        bonds_block = 0;
    uint64_t             bonds_cap = 0;
    md_render_mol_bonds_t bonds = {};
    bool                 frame_bonds = false;
    md_gpu_addr_t        secondary = 0;        // int32 per protein backbone segment
    uint32_t             secondary_count = 0;
    uint64_t             secondary_hash = 0;
    uint64_t             content[MD_RENDER_CONTENT_COUNT] = {};

    md_array(BridgeRep)         reps = nullptr;
    md_array(md_render_batch_t) batches = nullptr;

    // The shapes (emit_shapes): the points of their capsules, a colour each and a radius scale a capsule, on the GPU
    // and as they were uploaded; the dipole of each arrow, by the arrows' order
    md_gpu_addr_t         shape_pos = 0, shape_col = 0, shape_scale = 0;
    uint32_t              shape_cap = 0;          // points the buffers hold (and scales)
    md_array(md_float3_t) shape_pos_cpu = nullptr;
    md_array(uint32_t)    shape_col_cpu = nullptr;
    md_array(float)       shape_scale_cpu = nullptr;
    struct Dipole { md_attribute_id_t key; uint32_t index; };
    md_array(Dipole)      dipoles = nullptr;

    // The parts of the assembly shown (place_ops): the atoms of each (0: all) and its copies' transforms
    struct Part { md_gpu_addr_t atoms; uint32_t c0, c1; };
    md_array(Part) parts = nullptr;
    int            parts_of = -1;           // the assembly, -1: none
    uint64_t       parts_topology = 0;

    // The frame: md_compose's scene and md_post's frame at the scene's size, the histories at the image's, the image
    // kept (the frame, or the history written last); shared with GL, the image shown and the scene's depth
    md_compose_ctx_t compose = nullptr;
    md_post_ctx_t    post = nullptr;
    uint32_t  width = 0, height = 0;            // the scene's (the G-buffer's)
    uint32_t  out_width = 0, out_height = 0;    // the image's
    md_gpu_texture_t scene = nullptr, frame = nullptr;
    md_gpu_texture_t history[2] = {};
    int       hist = 0;                         // the history written last
    bool      hist_valid = false;               // it holds the image of the frame before, for the TAA to go on from
    md_gpu_texture_t image = nullptr;
    SharedTex depth, out;
    bool      exported = false;     // depth holds the scene drawn last
    bool      still = true;         // md_render did not rasterise this frame: nothing moved
    md_render_frame_t last = {};    // the frame drawn last (its batches are `batches`), for the masks alone
    bool      has_last = false;
    int       last_ui_frame = -1;   // the ImGui frame it was drawn in
    uint32_t  presents = 0;         // seeds the dither

    // GL: the image shown into the framebuffer; the depth as GL reads it (depth_gl)
    GLuint present_program = 0, depth_program = 0, vao = 0;
    GLint  ploc_image = -1, ploc_size = -1, dloc_depth = -1, dloc_size = -1;
    GLuint depth_fbo = 0, depth_tex = 0;
    uint32_t depth_tex_w = 0, depth_tex_h = 0;
    bool   depth_tex_valid = false;

    // md_render's answers to the picks, the newest last
    static constexpr int PICKS = 8;
    md_render_pick_t picks[PICKS] = {};
    uint32_t         num_picks = 0;
};

static void* gl_loader(const char* name) {
    GL3WglProc p = gl3wGetProcAddress(name);
    void* v = nullptr;
    static_assert(sizeof(p) == sizeof(v), "function and object pointers of one size");
    memcpy(&v, &p, sizeof(v));
    return v;
}

static void set_status(MdRenderBridge* b, const char* fmt, const char* detail = nullptr) {
    snprintf(b->status, sizeof(b->status), fmt, detail ? detail : "");
}

// ---- GL ----------------------------------------------------------------------------------------------------------

static const char* fullscreen_vs = R"(
#version 330 core
void main() {
    vec2 p = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
    gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}
)";

// The image shown, as md_post presented it (sRGB-encoded, dithered), into the framebuffer: its rows from the top, GL's
// from the bottom
static const char* present_fs = R"(
#ifdef RECT
uniform sampler2DRect u_image;
#define FETCH(t, p) texelFetch(t, p)
#else
uniform sampler2D u_image;
#define FETCH(t, p) texelFetch(t, p, 0)
#endif
uniform ivec2 u_size;
out vec4 out_color;
void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    p.y = u_size.y - 1 - p.y;
    out_color = vec4(FETCH(u_image, p).rgb, 1.0);
}
)";

// md_render's depth (reverse Z, 0 where nothing; rows from the top) as GL's window depth (rows from the bottom):
// md_render's projection is md_gl's with its depth reversed, d = (1 - ndc z) / 2
static const char* depth_fs = R"(
#ifdef RECT
uniform sampler2DRect u_depth;
#define FETCH(t, p) texelFetch(t, p)
#else
uniform sampler2D u_depth;
#define FETCH(t, p) texelFetch(t, p, 0)
#endif
uniform ivec2 u_size;
out vec4 out_depth;
void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    p.y = u_size.y - 1 - p.y;
    float d = FETCH(u_depth, p).r;
    out_depth = vec4(d > 0.0 ? 1.0 - d : 1.0);
}
)";

static GLuint gl_program(const char* fs_body, bool rect) {
    char fs_src[4096];
    snprintf(fs_src, sizeof(fs_src), "#version 330 core\n%s%s", rect ? "#define RECT 1\n" : "", fs_body);
    GLuint vs = gl::compile_shader_from_source(md_str_from_cstr(fullscreen_vs), GL_VERTEX_SHADER);
    GLuint fs = gl::compile_shader_from_source(md_str_from_cstr(fs_src), GL_FRAGMENT_SHADER);
    GLuint program = 0;
    if (vs && fs) {
        program = glCreateProgram();
        const GLuint shaders[] = {vs, fs};
        if (!gl::attach_link_detach(program, shaders, 2)) {
            glDeleteProgram(program);
            program = 0;
        }
    }
    if (vs) glDeleteShader(vs);
    if (fs) glDeleteShader(fs);
    return program;
}

static bool gl_init(MdRenderBridge* b, bool rect) {
    b->present_program = gl_program(present_fs, rect);
    b->depth_program   = gl_program(depth_fs, rect);
    if (!b->present_program || !b->depth_program) return false;
    b->ploc_image = glGetUniformLocation(b->present_program, "u_image");
    b->ploc_size  = glGetUniformLocation(b->present_program, "u_size");
    b->dloc_depth = glGetUniformLocation(b->depth_program, "u_depth");
    b->dloc_size  = glGetUniformLocation(b->depth_program, "u_size");
    glGenVertexArrays(1, &b->vao);
    glGenFramebuffers(1, &b->depth_fbo);
    return true;
}

// ---- Shared textures ---------------------------------------------------------------------------------------------

// Before md_gpu writes s again: GL is done with it (read or not)
static void gl_done(SharedTex* s) {
    if (!s->released) return;
    md_gpu_gl_begin_read(&s->gl);
    md_gpu_gl_end_read(&s->gl);
    s->released = false;
}

static void share_write_begin(MdRenderBridge* b, SharedTex* s) {
    gl_done(s);
    md_gpu_texture_share_acquire(b->stream, s->tex);
}

static void share_write_end(MdRenderBridge* b, SharedTex* s) {
    md_gpu_texture_share_release(b->stream, s->tex);
    s->released = true;
}

// GL reads what md_gpu released: between gl_read_begin and gl_read_end
static void gl_read_begin(SharedTex* s) {
    if (s->released) md_gpu_gl_begin_read(&s->gl);
}

static void gl_read_end(SharedTex* s) {
    if (!s->released) return;
    md_gpu_gl_end_read(&s->gl);
    s->released = false;
}

static void shared_free(SharedTex* s) {
    gl_done(s);
    if (s->imported) md_gpu_gl_release(&s->gl);
    if (s->tex) md_gpu_texture_destroy(s->tex);
    *s = SharedTex{};
}

static bool shared_make(MdRenderBridge* b, SharedTex* s, md_gpu_format_t format, uint32_t width, uint32_t height, const char* label) {
    md_gpu_texture_desc_t d = {};
    d.type   = MD_GPU_TEX_2D;
    d.format = format;
    d.usage  = MD_GPU_TEX_STORAGE | MD_GPU_TEX_SHARED;
    d.width  = width;
    d.height = height;
    d.label  = label;
    s->tex = md_gpu_texture_create(b->stream, &d);
    if (!s->tex) {
        set_status(b, "Could not make a texture shared with OpenGL: %s", md_gpu_last_error());
        return false;
    }
    if (!md_gpu_gl_import(&s->gl, s->tex, gl_loader)) {
        set_status(b, "OpenGL could not import md_gpu's texture (see the log): %s", md_gpu_last_error());
        return false;
    }
    s->imported = true;
    return true;
}

static md_gpu_texture_t gpu_texture(MdRenderBridge* b, md_gpu_format_t format, md_gpu_tex_usage_t usage, uint32_t w, uint32_t h, const char* label) {
    md_gpu_texture_desc_t d = {};
    d.type = MD_GPU_TEX_2D; d.format = format; d.usage = usage; d.width = w; d.height = h; d.label = label;
    return md_gpu_texture_create(b->stream, &d);
}

static void free_targets(MdRenderBridge* b) {
    md_gpu_texture_t* t[] = {&b->scene, &b->frame, &b->history[0], &b->history[1]};
    for (md_gpu_texture_t* x : t) {
        if (*x) md_gpu_texture_destroy(*x);
        *x = nullptr;
    }
    shared_free(&b->depth);
    shared_free(&b->out);
    b->image = nullptr;
    b->hist_valid = false;
    b->exported = false;
    b->width = b->height = b->out_width = b->out_height = 0;
}

// The frame's textures at the scene's size (width, height) and the image's (out_width, out_height)
static bool ensure_targets(MdRenderBridge* b, uint32_t width, uint32_t height, uint32_t out_width, uint32_t out_height) {
    if (b->scene && b->width == width && b->height == height && b->out_width == out_width && b->out_height == out_height) return true;
    free_targets(b);
    const md_gpu_tex_usage_t SS = MD_GPU_TEX_STORAGE | MD_GPU_TEX_SAMPLED;
    b->scene      = gpu_texture(b, MD_GPU_FORMAT_RGBA16_FLOAT, SS, width, height, "main view scene");
    b->frame      = gpu_texture(b, MD_GPU_FORMAT_RGBA16_FLOAT, SS, width, height, "main view frame");
    b->history[0] = gpu_texture(b, MD_GPU_FORMAT_RGBA16_FLOAT, SS, out_width, out_height, "main view history 0");
    b->history[1] = gpu_texture(b, MD_GPU_FORMAT_RGBA16_FLOAT, SS, out_width, out_height, "main view history 1");
    if (!b->scene || !b->frame || !b->history[0] || !b->history[1]) {
        set_status(b, "Could not make the main view's targets: %s", md_gpu_last_error());
        free_targets(b);
        return false;
    }
    if (!shared_make(b, &b->depth, MD_GPU_FORMAT_R32_FLOAT, width, height, "main view depth") ||
        !shared_make(b, &b->out, MD_GPU_FORMAT_RGBA8_UNORM, out_width, out_height, "main view image")) {
        free_targets(b);
        return false;
    }
    if (!b->present_program && !gl_init(b, b->out.gl.target != GL_TEXTURE_2D)) {
        set_status(b, "The main view's GL shaders did not compile%s");
        free_targets(b);
        return false;
    }
    b->width = width; b->height = height; b->out_width = out_width; b->out_height = out_height;
    return true;
}

// ---- The system --------------------------------------------------------------------------------------------------

static md_gpu_addr_t realloc_gpu(MdRenderBridge* b, md_gpu_addr_t old, size_t size) {
    if (old) md_gpu_free(b->stream, old);
    return size ? md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, size).gpu : 0;
}

static void free_parts(MdRenderBridge* b);

// Pairs of atoms (atom_a[i], atom_b[i]) of the kinds of their flags (md_bond_flags_t; null: covalent) into p, by kind,
// each with its index i as its id; those with an atom the system has not left out
static void upload_pairs(MdRenderBridge* b, PairSet* p, const uint32_t* atom_a, const uint32_t* atom_b, const uint32_t* flags, size_t count,
                         md_allocator_i* alloc) {
    const uint32_t n_atoms = b->atom_count;
    auto kind_of = [flags](size_t i) -> uint32_t { return flags ? MD_RENDER_MOL_BOND_KIND(md_bond_kind((md_bond_flags_t)flags[i])) : 0u; };
    auto valid = [&](size_t i) { return atom_a[i] < n_atoms && atom_b[i] < n_atoms; };
    uint32_t off[MD_RENDER_MOL_BOND_KIND_COUNT + 1] = {};
    for (size_t i = 0; i < count; ++i) if (valid(i)) off[kind_of(i) + 1] += 1;
    for (uint32_t k = 0; k < MD_RENDER_MOL_BOND_KIND_COUNT; ++k) off[k + 1] += off[k];
    const uint32_t num = off[MD_RENDER_MOL_BOND_KIND_COUNT];
    memset(p->offset, 0, sizeof(p->offset));
    if (num > p->cap) {
        p->pairs = realloc_gpu(b, p->pairs, (size_t)num * 8);
        p->ids   = realloc_gpu(b, p->ids, (size_t)num * 4);
        p->cap   = p->pairs && p->ids ? num : 0;
    }
    if (num == 0 || !p->cap) return;
    md_temp_scope_t temp = md_temp_begin_in(alloc);
    uint32_t* pairs = (uint32_t*)md_vm_arena_push(alloc, (size_t)num * 8);
    uint32_t* ids   = (uint32_t*)md_vm_arena_push(alloc, (size_t)num * 4);
    uint32_t at[MD_RENDER_MOL_BOND_KIND_COUNT];
    memcpy(at, off, sizeof(at));
    for (size_t i = 0; i < count; ++i) {
        if (!valid(i)) continue;
        const uint32_t j = at[kind_of(i)]++;
        pairs[2 * j] = atom_a[i];
        pairs[2 * j + 1] = atom_b[i];
        ids[j] = (uint32_t)i;
    }
    md_gpu_upload(b->stream, p->pairs, pairs, (size_t)num * 8);
    md_gpu_upload(b->stream, p->ids, ids, (size_t)num * 4);
    md_temp_end(temp);
    memcpy(p->offset, off, sizeof(off));
}

static void free_system(MdRenderBridge* b) {
    if (b->has_topo) md_gpu_free(b->stream, b->topo.base);
    if (b->dynamic)   md_gpu_free(b->stream, b->dynamic);
    if (b->positions) md_gpu_free(b->stream, b->positions);
    if (b->marks)     md_gpu_free(b->stream, b->marks);
    if (b->secondary) md_gpu_free(b->stream, b->secondary);
    if (b->bonds_block) md_gpu_free(b->stream, b->bonds_block);
    b->bonds_block = 0;
    b->bonds_cap = 0;
    b->frame_bonds = false;
    for (size_t i = 0; i < md_array_size(b->reps); ++i) free_rep(b->stream, &b->reps[i]);
    md_array_shrink(b->reps, 0);
    free_parts(b);
    b->has_topo = false;
    b->topo = {};
    b->dynamic = b->positions = b->marks = b->secondary = 0;
    b->secondary_count = 0;
    b->secondary_hash = 0;
    b->atom_count = 0;
    b->updated = false;
    b->regroup_frame = -1;
}

static size_t bitset_words32(uint32_t count) { return ((size_t)count + 63) / 64 * 2; }

// A set over [0, count) as words of 32 bits: those of 64 in the order they are stored
static void dense_words(uint64_t* out, uint32_t count, md_bitset_t set) {
    const size_t n = ((size_t)count + 63) / 64;
    for (size_t k = 0; k < n; ++k) out[k] = md_bitset_word(set, (uint32_t)(k * 64));
}

static void upload_set(MdRenderBridge* b, md_gpu_addr_t dst, uint32_t count, md_bitset_t set, md_allocator_i* alloc) {
    const size_t n = ((size_t)count + 63) / 64;
    md_temp_scope_t temp = md_temp_begin_in(alloc);
    uint64_t* words = (uint64_t*)md_vm_arena_push(alloc, n * sizeof(uint64_t));
    dense_words(words, count, set);
    md_gpu_upload(b->stream, dst, words, n * sizeof(uint64_t));
    md_temp_end(temp);
}

// The topology of the system shown, made again when it changed
static bool ensure_system(MdRenderBridge* b, ApplicationState* state) {
    const md_system_t& sys = state->mold.sys;
    const uint32_t n = (uint32_t)sys.atom.count;
    if (b->has_topo && b->atom_count == n && b->topo.topology_version == sys.topology_version) return true;
    free_system(b);
    if (n == 0) return false;

    md_allocator_i* alloc = state->allocator.frame;
    md_temp_scope_t temp = md_temp_begin_in(alloc);
    md_render_mol_plan_t plan = {};
    bool ok = md_render_mol_plan(&plan, &sys, nullptr, temp);
    if (ok) {
        const md_gpu_addr_t base = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, plan.size).gpu;
        ok = base && md_render_mol_write(b->stream, &b->topo, &plan, base, b->content);
        if (!ok && base) md_gpu_free(b->stream, base);
    }
    md_temp_end(temp);
    if (!ok) {
        set_status(b, "md_render_mol could not plan the system (see the log)%s");
        return false;
    }
    b->has_topo  = true;
    b->atom_count = n;
    b->dynamic   = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, b->topo.dynamic_size_velocity).gpu;    // velocities derived
    b->positions = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, (size_t)n * 12).gpu;
    b->marks     = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, n).gpu;
    if (!b->dynamic || !b->positions || !b->marks) {
        set_status(b, "Out of GPU memory for the system%s");
        free_system(b);
        return false;
    }
    md_gpu_memset(b->stream, b->marks, 0, n);
    // All of it, then
    b->pending |= MolBit_DirtyPosition | MolBit_DirtyFlags | MolBit_DirtyBonds | MolBit_ResetBackboneHistory;
    // The representations' arrays are of the atoms: given again for the new system
    flag_all_representations_as_dirty(state);
    return true;
}

static md_render_mol_atoms_t atoms_of(const MdRenderBridge* b, const ApplicationState* state) {
    md_render_mol_atoms_t a = {};
    a.position  = b->positions;
    a.marks     = b->marks;
    md_render_mol_frame_cell(a.cell, &state->mold.state.unitcell);
    // The velocities: the displacement since the last update at the minimum image, per update (time counts them),
    // so that what moved across the cell moves by what is left of it rather than streaking over the system
    a.velocity_dt = 1.0f;
    a.bonds = b->frame_bonds ? &b->bonds : nullptr;
    return a;
}

// The bonds of the frame shown while they are tracked (bond_tracking_delta, bond_tracking_chem), as md_gl gets them
// (update_md_buffers, md_gl_mol_set_bond_orders): the system's broken and those of the region of the frame's chemistry
// hidden, the region's and those formed drawn as the frame has them, with their orders. Picks report md_gl's: a bond
// of the system by its index, one formed the system's bond count past its index in the delta; a ring of the region
// the system's ring count past its index in the chemistry's.
static void update_frame_bonds(MdRenderBridge* b, ApplicationState* state) {
    const md_bond_delta_t* delta = bond_tracking_delta(state);
    const md_chem_delta_t* chem = bond_tracking_chem(state);
    b->frame_bonds = false;
    b->content[MD_RENDER_CONTENT_GEOMETRY] += 1;
    if (!delta || (delta->num_formed == 0 && md_bitset_empty(delta->broken))) return;
    md_allocator_i* alloc = state->allocator.frame;
    md_temp_scope_t temp = md_temp_begin_in(alloc);
    md_render_mol_bonds_plan_t plan;
    if (md_render_mol_bonds_plan(&plan, &state->mold.sys, nullptr, &state->mold.state, delta, chem, temp)) {
        if (plan.size > b->bonds_cap) {
            b->bonds_block = realloc_gpu(b, b->bonds_block, (size_t)plan.size);
            b->bonds_cap = b->bonds_block ? plan.size : 0;
        }
        b->frame_bonds = b->bonds_block && md_render_mol_bonds_write(b->stream, &b->bonds, &plan, b->bonds_block, b->content);
    }
    md_temp_end(temp);
}

// What changed since the scene was last made: positions, the marks, the bonds of the frame
static bool take_up(MdRenderBridge* b, ApplicationState* state) {
    const uint32_t dirty = b->pending;
    b->pending = 0;
    // The frame's bonds' sides are the update's: one follows
    if (dirty & MolBit_DirtyBonds) update_frame_bonds(b, state);
    if (dirty & MolBit_DirtyFlags) {
        // The selection and the highlight, and the group of the highlight each atom is in (of the last group it is in,
        // as update_md_buffers has it): the colour it is marked in
        md_allocator_i* alloc = state->allocator.frame;
        md_temp_scope_t temp = md_temp_begin_in(alloc);
        uint8_t* marks = (uint8_t*)md_vm_arena_push(alloc, b->atom_count);
        md_render_marks_write(marks, b->atom_count, md_bitmask_view(&state->selection.selection_mask), md_bitmask_view(&state->selection.highlight_mask));
        const auto& groups = state->selection.highlight_groups;
        for (size_t g = 0; g < groups.count; ++g) {
            md_render_marks_tag(marks, b->atom_count, groups.sets[g], (uint32_t)(g % HIGHLIGHT_GROUP_COLORS));
        }
        md_gpu_upload(b->stream, b->marks, marks, b->atom_count);
        md_temp_end(temp);
        b->content[MD_RENDER_CONTENT_STATE] += 1;
    }
    if ((dirty & (MolBit_DirtyPosition | MolBit_DirtyRadius | MolBit_ResetBackboneHistory | MolBit_DirtyBonds)) || !b->updated) {
        if (!md_render_mol_upload_state(b->stream, b->positions, &state->mold.state, b->content)) {
            set_status(b, "md_render_mol could not upload the positions (see the log)%s");
            return false;
        }
        md_render_mol_update_flags_t flags = 0;
        // The small molecules regrouped on a new frame of the trajectory, not between frames
        const int64_t frame = md_state_has_frame(&state->mold.state) ? (int64_t)md_state_frame_nearest(&state->mold.state) : 0;
        if (frame != b->regroup_frame || !b->updated) {
            flags |= MD_RENDER_MOL_UPDATE_REGROUP;
            b->regroup_frame = frame;
        }
        if (dirty & MolBit_ResetBackboneHistory) flags |= MD_RENDER_MOL_UPDATE_RESET_BACKBONE;
        // A jump (a seek, another trajectory) is not a motion: at rest
        if (dirty & (MolBit_ClearVelocity | MolBit_ResetBackboneHistory)) flags |= MD_RENDER_MOL_UPDATE_RESET_VELOCITY;
        const md_render_mol_atoms_t atoms = atoms_of(b, state);
        if (!md_render_mol_update(b->mol, b->stream, &b->topo, b->dynamic, &atoms, flags, b->content)) {
            set_status(b, "md_render_mol could not update the positions (see the log)%s");
            return false;
        }
        b->updated = true;
        b->time += 1.0;
    }
    return true;
}

// The secondary structure shown, for the cartoons' profiles: its hash, 0 for none (all coil)
static uint64_t update_secondary(MdRenderBridge* b, ApplicationState* state) {
    const md_system_t& sys = state->mold.sys;
    const uint32_t count = (uint32_t)sys.protein_backbone.segment.count;
    const md_secondary_structure_t* ss = count ? displayed_secondary_structure(state) : nullptr;
    if (!ss) return 0;
    const uint64_t hash = md_hash64(ss, count * sizeof(md_secondary_structure_t), 1) | 1;
    if (count != b->secondary_count) {
        b->secondary = realloc_gpu(b, b->secondary, count * sizeof(int32_t));
        b->secondary_count = count;
        b->secondary_hash = 0;
    }
    if (hash != b->secondary_hash && b->secondary) {
        static_assert(sizeof(md_secondary_structure_t) == sizeof(int32_t), "an int32 per segment");
        md_gpu_upload(b->stream, b->secondary, ss, count * sizeof(int32_t));
        b->secondary_hash = hash;
    }
    return b->secondary ? hash : 0;
}

// ---- The representations -----------------------------------------------------------------------------------------

static BridgeRep* find_rep(MdRenderBridge* b, uint32_t id) {
    for (size_t i = 0; i < md_array_size(b->reps); ++i) {
        if (b->reps[i].id == id) return &b->reps[i];
    }
    return nullptr;
}

static bool type_drawn(RepresentationType t) {
    return t == RepresentationType::SpaceFill || t == RepresentationType::Licorice || t == RepresentationType::BallAndStick ||
           t == RepresentationType::Ribbons || t == RepresentationType::Cartoon || t == RepresentationType::Bonds;
}

static bool active(const ApplicationState& state) {
    const MdRenderBridge* b = state.md_render_bridge;
    return b && b->ok;
}

// The shapes md_render draws (emit_shapes)
static bool shape_shown(const ApplicationState& state, const Representation& rep) {
    if (rep.type == RepresentationType::DipoleMoment) return representation_shown(state, rep);
    if (rep.type == RepresentationType::SimulationBox) return representation_shown(state, rep) && rep.type_is_valid;
    return false;
}

// Whether md_render draws the representation (as it is now: shown, of the system's atoms)
static bool draws(const ApplicationState& state, const Representation& rep) {
    if (!active(state) || !state.md_render_bridge->has_topo) return false;
    if (rep.type == RepresentationType::DipoleMoment || rep.type == RepresentationType::SimulationBox) return shape_shown(state, rep);
    if (!type_drawn(rep.type) || !rep.type_is_valid || !representation_shown(state, rep)) return false;
    const BridgeRep* r = find_rep(state.md_render_bridge, rep.node.id);
    return r && r->color && r->visible && r->atom_count == state.md_render_bridge->atom_count;
}

static BridgeRep* find_or_add_rep(MdRenderBridge* b, uint32_t id, md_allocator_i* alloc) {
    BridgeRep* r = find_rep(b, id);
    if (r) return r;
    BridgeRep nr = {};
    nr.id = id;
    md_array_push(b->reps, nr, alloc);
    return &b->reps[md_array_size(b->reps) - 1];
}

void md_render_bridge_rep_pairs(ApplicationState* state, const Representation& rep, const uint32_t* atom_a, const uint32_t* atom_b,
                                const uint32_t* flags, size_t count) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!b || !b->ok || !b->has_topo || rep.type != RepresentationType::Bonds) return;
    BridgeRep* r = find_or_add_rep(b, rep.node.id, state->allocator.persistent);
    // Each pair's index the one it picks as (representation_weak_bond_draw_ops)
    upload_pairs(b, &r->pairs, atom_a, atom_b, flags, count, state->allocator.frame);
    b->content[MD_RENDER_CONTENT_GEOMETRY] += 1;
}

void md_render_bridge_rep_colors(ApplicationState* state, const Representation& rep, const uint32_t* colors, size_t count) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!b || !b->ok || !type_drawn(rep.type)) return;
    if (!b->has_topo || count != b->atom_count || !colors) return;
    BridgeRep* r = find_or_add_rep(b, rep.node.id, state->allocator.persistent);
    if (r->atom_count != count) {
        r->color   = realloc_gpu(b, r->color, count * sizeof(uint32_t));
        r->visible = realloc_gpu(b, r->visible, bitset_words32((uint32_t)count) * 4);
        r->visible_hash = 0;
        r->atom_count = (uint32_t)count;
    }
    if (!r->color || !r->visible) return;
    md_gpu_upload(b->stream, r->color, colors, count * sizeof(uint32_t));
    b->content[MD_RENDER_CONTENT_APPEARANCE] += 1;
    // What it shows is geometry: rasterised again only when it changed
    const uint64_t hash = md_bitset_hash64(rep.atom_mask, count) | 1;
    if (hash != r->visible_hash) {
        upload_set(b, r->visible, (uint32_t)count, rep.atom_mask, state->allocator.frame);
        r->visible_hash = hash;
        b->content[MD_RENDER_CONTENT_GEOMETRY] += 1;
    }
}

// The style of each kind of bond of a representation: md_render_mol's, its coordinate bonds solid where it has them so
static void bond_styles(md_render_mol_bond_style_t out[MD_RENDER_MOL_BOND_KIND_COUNT], const Representation& rep) {
    for (int k = 0; k < MD_RENDER_MOL_BOND_KIND_COUNT; ++k) out[k] = md_render_mol_bond_styles[k];
    if (!rep.coordinate_bonds_dashed) {
        md_render_mol_bond_style_t& c = out[MD_RENDER_MOL_BOND_KIND(MD_BOND_KIND_COORDINATE)];
        c.radius = 1.0f;
        c.dash = 0.0f;
        c.gap = 0.0f;
    }
    if (rep.bond_color == BondColorMode::Uniform) {
        uint32_t color = convert_color(rep.bond_base_color);
        if (rep.tint_scale > 0.0f || rep.saturation < 1.0f) {
            tint_colors(&color, 1, convert_color(rep.tint_color), rep.tint_scale, rep.saturation);
        }
        for (int k = 0; k < MD_RENDER_MOL_BOND_KIND_COUNT; ++k) out[k].color = color;
    }
}

// The style of each kind of the pairs of a Bonds representation, as md_gl's dashes (md_gl.h): by kind md_render_mol's
// (those that are not covalent dashed), none, or all (the covalent, elastic and metallic ones in a style of their own)
static void pair_styles(md_render_mol_bond_style_t out[MD_RENDER_MOL_BOND_KIND_COUNT], const Representation& rep) {
    for (int k = 0; k < MD_RENDER_MOL_BOND_KIND_COUNT; ++k) {
        out[k] = md_render_mol_bond_styles[k];
        if (rep.bond_dashes == BondDashes::None) {
            out[k].radius = 1.0f;
            out[k].dash = out[k].gap = 0.0f;
        } else if (rep.bond_dashes == BondDashes::All && out[k].dash == 0.0f) {
            out[k].radius = 1.0f;
            out[k].dash = 2.4f;
            out[k].gap = 0.6f;
        }
    }
    if (rep.bond_color == BondColorMode::Uniform) {
        uint32_t color = convert_color(rep.bond_base_color);
        if (rep.tint_scale > 0.0f || rep.saturation < 1.0f) {
            tint_colors(&color, 1, convert_color(rep.tint_color), rep.tint_scale, rep.saturation);
        }
        for (int k = 0; k < MD_RENDER_MOL_BOND_KIND_COUNT; ++k) out[k].color = color;
    }
}

static uint32_t bond_kinds(const Representation& rep) {
    uint32_t kinds = MD_RENDER_MOL_BONDS(MD_BOND_KIND_COVALENT) | MD_RENDER_MOL_BONDS(MD_BOND_KIND_COORDINATE) | MD_RENDER_MOL_BONDS(MD_BOND_KIND_METALLIC);
    if (rep.elastic_bonds) kinds |= MD_RENDER_MOL_BONDS(MD_BOND_KIND_ELASTIC);
    if (rep.interaction_bonds) {
        kinds |= MD_RENDER_MOL_BONDS(MD_BOND_KIND_CONTACT) | MD_RENDER_MOL_BONDS(MD_BOND_KIND_HYDROGEN) |
                 MD_RENDER_MOL_BONDS(MD_BOND_KIND_IONIC) | MD_RENDER_MOL_BONDS(MD_BOND_KIND_HALOGEN);
    }
    return kinds;
}

// The profile of a cartoon: of the secondary structure shown and its scales
static md_gpu_addr_t cartoon_profile(MdRenderBridge* b, BridgeRep* r, const Representation& rep, uint64_t secondary) {
    const float scale[3] = {rep.scale.x, rep.scale.z, rep.scale.y};     // coil, helix, sheet
    uint64_t key = md_hash64(scale, sizeof(scale), secondary) | 1;
    if (!r->profile) {
        r->profile = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, md_render_mol_cartoon_size(&b->topo)).gpu;
        r->profile_key = 0;
        if (!r->profile) return 0;
    }
    if (key != r->profile_key) {
        if (!md_render_mol_cartoon_profile(b->mol, b->stream, &b->topo, r->profile, secondary ? b->secondary : 0, scale, b->content)) return 0;
        r->profile_key = key;
    }
    return r->profile;
}

static bool emit(MdRenderBridge* b, const md_render_mol_atoms_t& atoms, const md_render_mol_rep_t& mr, md_allocator_i* alloc) {
    const size_t at = md_array_size(b->batches);
    uint32_t cap = (uint32_t)(md_array_capacity(b->batches) - at);
    md_array_resize(b->batches, at + cap, alloc);
    uint32_t n = md_render_mol_emit(&b->topo, b->dynamic, &atoms, &mr, b->batches + at, cap);
    if (n > cap) {
        md_array_resize(b->batches, at + n, alloc);
        n = md_render_mol_emit(&b->topo, b->dynamic, &atoms, &mr, b->batches + at, n);
    }
    md_array_resize(b->batches, at + n, alloc);
    return n > 0;
}

// md_gl's view as md_render's camera: the same view, the projection's depth reversed (d = (1 - ndc z) / 2, so that
// md_gl's depth is 1 - d) and its jitter given apart, in pixels (x right, y down)
static void camera_of(md_render_camera_t* cam, float jitter[2], const ApplicationState& state, uint32_t width, uint32_t height) {
    const ViewParam& p = state.view.param;
    memcpy(cam->view, &p.matrix.curr.view.elem[0][0], sizeof(cam->view));
    const md_float4x4_t& P = p.matrix.curr.proj_no_jitter;
    const md_float4x4_t& J = p.matrix.curr.proj;
    for (int c = 0; c < 4; ++c) {
        for (int r = 0; r < 4; ++r) cam->proj[c * 4 + r] = P.elem[c][r];
        cam->proj[c * 4 + 2] = 0.5f * (P.elem[c][3] - P.elem[c][2]);
    }
    // The shift of the image by the jitter, in NDC: of a perspective by column 2 (times z over w = -z), of an
    // orthographic one by column 3; either way that of the point at z = -1, whose w is 1
    const float dx = (P.elem[2][0] - J.elem[2][0]) + (J.elem[3][0] - P.elem[3][0]);
    const float dy = (P.elem[2][1] - J.elem[2][1]) + (J.elem[3][1] - P.elem[3][1]);
    jitter[0] =  dx * 0.5f * (float)width;
    jitter[1] = -dy * 0.5f * (float)height;
}

// ---- Assemblies --------------------------------------------------------------------------------------------------

static void free_parts(MdRenderBridge* b) {
    for (size_t i = 0; i < md_array_size(b->parts); ++i) {
        if (b->parts[i].atoms) md_gpu_free(b->stream, b->parts[i].atoms);
    }
    md_array_shrink(b->parts, 0);
    b->parts_of = -1;
    b->parts_topology = 0;
}

// The parts of the assembly shown, as place_ops places md_gl's: each the atoms of its instances (of all of them: every
// atom, 0) in the copies its operators make. None while the view's copies are not of it (yet): then as without one.
static void update_parts(MdRenderBridge* b, ApplicationState* state) {
    const md_system_t& sys = state->mold.sys;
    const md_assembly_data_t& a = sys.assembly;
    const int assembly = state->view.assembly;
    if (assembly < 0 || (size_t)assembly >= a.count) {
        if (md_array_size(b->parts) || b->parts_of >= 0) free_parts(b);
        return;
    }
    const size_t part_beg = a.group_offset[assembly], part_end = a.group_offset[assembly + 1];
    const size_t num_offsets = md_array_size(state->view.assembly_copy_offset);
    const size_t num_parts = num_offsets > 0 ? num_offsets - 1 : 0;
    if (num_parts == 0 || num_parts != part_end - part_beg || num_parts >= 0xFFFF) {
        free_parts(b);
        return;
    }
    if (b->parts_of == assembly && b->parts_topology == sys.topology_version && md_array_size(b->parts) == num_parts) return;
    free_parts(b);
    const uint32_t n = b->atom_count;
    const size_t num_words = ((size_t)n + 63) / 64;
    md_allocator_i* alloc = state->allocator.frame;
    md_temp_scope_t temp = md_temp_begin_in(alloc);
    uint64_t* words = (uint64_t*)md_vm_arena_push(alloc, num_words * sizeof(uint64_t));
    for (size_t p = part_beg; p < part_end; ++p) {
        const AssemblyGroup g = assembly_group(a, p);
        MdRenderBridge::Part part = {};
        part.c0 = state->view.assembly_copy_offset[p - part_beg];
        part.c1 = state->view.assembly_copy_offset[p - part_beg + 1];
        if ((size_t)(g.inst_end - g.inst) < sys.instance.count) {
            memset(words, 0, num_words * sizeof(uint64_t));
            for (const md_instance_idx_t* k = g.inst; k < g.inst_end; ++k) {
                const md_urange_t r = md_system_instance_atom_range(&sys, (size_t)*k);
                for (uint32_t i = r.beg; i < r.end && i < n; ++i) words[i >> 6] |= 1ull << (i & 63);
            }
            part.atoms = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, num_words * sizeof(uint64_t)).gpu;
            if (part.atoms) md_gpu_upload(b->stream, part.atoms, words, num_words * sizeof(uint64_t));
        }
        md_array_push(b->parts, part, state->allocator.persistent);
    }
    md_temp_end(temp);
    b->parts_of = assembly;
    b->parts_topology = sys.topology_version;
    b->content[MD_RENDER_CONTENT_GEOMETRY] += 1;
}

// A representation in the copies of the scene: with an assembly shown once per part, of its atoms in its copies (the
// part in the low bits of the key); else in its periodic images
static void emit_placed(MdRenderBridge* b, const ApplicationState* state, const md_render_mol_atoms_t& atoms, md_render_mol_rep_t mr,
                        const Representation& rep, md_allocator_i* alloc) {
    const size_t num_parts = md_array_size(b->parts);
    if (num_parts > 0) {
        for (size_t p = 0; p < num_parts; ++p) {
            const MdRenderBridge::Part& part = b->parts[p];
            if (part.c1 <= part.c0 || part.c1 > md_array_size(state->view.assembly_transforms)) continue;
            md_render_mol_atoms_t pa = atoms;
            pa.visible = part.atoms;
            md_render_mol_rep_t pr = mr;
            pr.key = mr.key | (uint64_t)(p + 1);
            pr.transform = &state->view.assembly_transforms[part.c0].elem[0][0];
            pr.transform_count = part.c1 - part.c0;
            emit(b, pa, pr, alloc);
        }
        return;
    }
    const size_t copies = md_array_size(state->view.periodic_transforms);
    if (rep.periodic_images && copies > 1) {
        mr.transform = &state->view.periodic_transforms[0].elem[0][0];
        mr.transform_count = (uint32_t)copies;
    }
    emit(b, atoms, mr, alloc);
}

// The transform of the copy a hit is in, of the representation of key (emit_placed), if it is in one
static bool copy_transform(const MdRenderBridge* b, const ApplicationState& state, uint64_t key, uint32_t instance, md_float4x4_t* out) {
    const uint32_t part = (uint32_t)(key & 0xFFFF);
    if (part > 0) {
        if (part > md_array_size(b->parts)) return false;
        const size_t c = (size_t)b->parts[part - 1].c0 + instance;
        if (c >= md_array_size(state.view.assembly_transforms)) return false;
        *out = state.view.assembly_transforms[c];
        return true;
    }
    if (instance < md_array_size(state.view.periodic_transforms) && md_array_size(state.view.periodic_transforms) > 1) {
        *out = state.view.periodic_transforms[instance];
        return true;
    }
    return false;
}

// ---- Shapes ------------------------------------------------------------------------------------------------------

// The dipoles' arrows and the cell's edges: capsules of md_render's over points of their own. Their keys are past the
// representations' (below 2^56), with a low byte that is none of md_render_mol's parts; a box's carry its
// representation. Shafts and heads hold the arrows of every dipole representation, in the order of `dipoles`.
static constexpr uint64_t SHAPE_KEY = 0xFFull << 56;
static uint64_t shape_key(uint32_t rep_id, uint32_t what) { return SHAPE_KEY | ((uint64_t)(rep_id + 1) << 16) | ((uint64_t)what << 8) | 0xFFu; }
enum { SHAPE_SHAFTS, SHAPE_HEADS, SHAPE_CELL, SHAPE_IMAGES, SHAPE_AXES, SHAPE_ORIGIN };
static bool is_shape_key(uint64_t key) { return (key >> 56) == 0xFF && (key & 0xFF) == 0xFF; }

// A batch of capsules over the shapes' points [first, first + 2 count), with their radius scales from `scale` (or
// none: -1)
struct ShapeBatch {
    uint64_t key;
    uint32_t first, count;
    int64_t  scale;
    md_render_capsule_style_t style;
    const float* transform;
    uint32_t     transform_count;
};

// The batches of the shapes this frame, into b->batches. The points are uploaded when they changed (in place, the
// content's counters bumped, or where there is room for them), so that a frame of the same ones draws nothing.
static void emit_shapes(MdRenderBridge* b, ApplicationState* state, md_allocator_i* alloc) {
    md_allocator_i* temp = state->allocator.frame;
    md_array(md_float3_t) pos = nullptr;
    md_array(uint32_t)    col = nullptr;
    md_array(float)       scale = nullptr;
    md_array(ShapeBatch)  sb = nullptr;
    md_array_shrink(b->dipoles, 0);
    const size_t num_reps = md_array_size(state->representation.reps);

    // The arrows: the shafts, then the heads, of the radii their scales hold
    md_array(DipoleArrow) arrows = nullptr;
    md_array(uint32_t)    arrow_col = nullptr;
    for (size_t i = 0; i < num_reps; ++i) {
        const Representation& rep = state->representation.reps[i];
        DipoleArrow a;
        if (rep.type != RepresentationType::DipoleMoment || !shape_shown(*state, rep) || !representation_dipole_arrow(&a, *state, rep)) continue;
        md_array_push(arrows, a, temp);
        md_array_push(arrow_col, convert_color(rep.dipole.color), temp);
        md_array_push(b->dipoles, (MdRenderBridge::Dipole{rep.dipole.dipole_key, rep.dipole.dipole_index}), state->allocator.persistent);
    }
    const uint32_t num_arrows = (uint32_t)md_array_size(arrows);
    if (num_arrows) {
        for (int head = 0; head < 2; ++head) {
            ShapeBatch x = {};
            x.key   = shape_key(0, head ? SHAPE_HEADS : SHAPE_SHAFTS);
            x.first = (uint32_t)md_array_size(pos);
            x.count = num_arrows;
            x.scale = (int64_t)md_array_size(scale);
            x.style.radius = 1.0f;
            x.style.flat   = true;
            x.style.tapered = head != 0;                    // to a point
            for (uint32_t k = 0; k < num_arrows; ++k) {
                const DipoleArrow& a = arrows[k];
                md_array_push(pos, head ? a.neck : a.beg, temp);
                md_array_push(pos, head ? a.tip : a.neck, temp);
                md_array_push(col, arrow_col[k], temp);
                md_array_push(col, arrow_col[k], temp);
                md_array_push(scale, head ? a.head_radius : a.shaft_radius, temp);
            }
            md_array_push(sb, x, temp);
        }
    }

    // The cells: of each box its own and those of the periodic images, their edges as wide on screen all along, joined
    // by their round ends (main.cpp's draw_simulation_box, which md_gl draws). With axes, those of its own cell along
    // a, b and c from the origin red, green and blue, flat-ended: the ends of the edges they meet, and the origin's
    // ball of the box's colour, close them. The axes and the ball are a hundredth wider than the edges, so that the
    // edges of the images they lie on never show through them (md_gl draws its own cell first, and the depth test
    // keeps the first of equal depths; md_render's raster has no order).
    const float AXIS_WIDER = 1.01f;
    const md_unitcell_t& cell = state->mold.state.unitcell;
    if (md_unitcell_flags(&cell) != 0) {
        md_float3x3_t A = {0};
        md_unitcell_A_extract_float(A.elem, &cell);
        const md_float4x4_t M = state->operations.state_rotation * md_float4x4_from_float3x3(A);
        md_float3_t corner[8];       // bit k of i: the coordinate along cell vector k
        for (int i = 0; i < 8; ++i) {
            corner[i] = md_float3_from_float4(md_float4x4_mul_float4(M, md_float4_set((float)(i & 1), (float)((i >> 1) & 1), (float)((i >> 2) & 1), 1.0f)));
        }
        const uint32_t axis_color[3] = { convert_color(md_float4_set(0.85f, 0.15f, 0.15f, 1.0f)), convert_color(md_float4_set(0.15f, 0.70f, 0.15f, 1.0f)), convert_color(md_float4_set(0.15f, 0.30f, 0.90f, 1.0f)) };
        // The thickness in points, as md_gl's: the view is the whole window, DisplaySize in points; never thinner than
        // 1.2 pixels shown (0.6 a radius), whatever the render scale; 0: a pixel
        const float h_render = (float)MD_MAX(state->gbuffer.height, 1U);
        const float h_px     = (float)MD_MAX(state->gbuffer.out_height, 1U);
        const float h_pt     = ImGui::GetIO().DisplaySize.y > 0.0f ? ImGui::GetIO().DisplaySize.y : h_px;
        const size_t num_xf  = md_array_size(state->view.periodic_transforms);
        for (size_t i = 0; i < num_reps; ++i) {
            const Representation& rep = state->representation.reps[i];
            if (rep.type != RepresentationType::SimulationBox || !shape_shown(*state, rep)) continue;
            const uint32_t color = convert_color(md_float4_set(rep.box.color.x, rep.box.color.y, rep.box.color.z, 1.0f));
            const uint32_t base = (uint32_t)md_array_size(pos);
            auto edge = [&](int i0, int i1, uint32_t c) {
                md_array_push(pos, corner[i0], temp);
                md_array_push(pos, corner[i1], temp);
                md_array_push(col, c, temp);
                md_array_push(col, c, temp);
            };
            for (int k = 0; k < 3; ++k) edge(0, 1 << k, color);                         // base + 0: from the origin
            for (int c = 1; c < 8; ++c) for (int k = 0; k < 3; ++k) if (!(c & (1 << k))) edge(c, c | (1 << k), color);
            edge(0, 0, color);                                                          // base + 24: the origin
            for (int k = 0; k < 3; ++k) edge(0, 1 << k, axis_color[k]);                 // base + 26: the axes

            md_render_capsule_style_t st = {};
            st.radius_px = rep.box.thickness > 0.0f ? MD_MAX(0.5f * rep.box.thickness * h_render / h_pt, 0.6f * h_render / h_px) : 0.5f * h_render / h_px;
            ShapeBatch x = {};
            x.scale = -1;
            x.style = st;
            x.key   = shape_key(rep.node.id, SHAPE_CELL);
            x.first = rep.box.axes ? base + 6 : base;
            x.count = rep.box.axes ? 9 : 12;
            md_array_push(sb, x, temp);
            if (rep.periodic_images && num_xf > 1) {
                x.key   = shape_key(rep.node.id, SHAPE_IMAGES);
                x.first = base;
                x.count = 12;
                x.transform = &state->view.periodic_transforms[1].elem[0][0];
                x.transform_count = (uint32_t)(num_xf - 1);
                md_array_push(sb, x, temp);
            }
            if (rep.box.axes) {
                x.key   = shape_key(rep.node.id, SHAPE_ORIGIN);
                x.first = base + 24;
                x.count = 1;
                x.style.radius_px = st.radius_px * AXIS_WIDER;
                x.transform = nullptr;
                x.transform_count = 0;
                md_array_push(sb, x, temp);
                x.key   = shape_key(rep.node.id, SHAPE_AXES);
                x.first = base + 26;
                x.count = 3;
                x.style.flat = true;
                md_array_push(sb, x, temp);
            }
        }
    }

    const uint32_t np = (uint32_t)md_array_size(pos);
    if (np == 0) return;
    // Uploaded when changed: in place with the counter of the content bumped, or into larger buffers (new addresses)
    if (np > b->shape_cap) {
        b->shape_pos   = realloc_gpu(b, b->shape_pos, (size_t)np * sizeof(md_float3_t));
        b->shape_col   = realloc_gpu(b, b->shape_col, (size_t)np * sizeof(uint32_t));
        b->shape_scale = realloc_gpu(b, b->shape_scale, (size_t)np * sizeof(float));
        b->shape_cap   = b->shape_pos && b->shape_col && b->shape_scale ? np : 0;
        md_array_shrink(b->shape_pos_cpu, 0);
        md_array_shrink(b->shape_col_cpu, 0);
        md_array_shrink(b->shape_scale_cpu, 0);
        if (!b->shape_cap) return;
    }
    md_allocator_i* keep = state->allocator.persistent;
    const size_t ns = md_array_size(scale);
    const bool geo = md_array_size(b->shape_pos_cpu) != np || memcmp(b->shape_pos_cpu, pos, np * sizeof(md_float3_t)) ||
                     md_array_size(b->shape_scale_cpu) != ns || (ns && memcmp(b->shape_scale_cpu, scale, ns * sizeof(float)));
    const bool app = md_array_size(b->shape_col_cpu) != np || memcmp(b->shape_col_cpu, col, np * sizeof(uint32_t));
    if (geo) {
        md_gpu_upload(b->stream, b->shape_pos, pos, np * sizeof(md_float3_t));
        if (ns) md_gpu_upload(b->stream, b->shape_scale, scale, ns * sizeof(float));
        md_array_resize(b->shape_pos_cpu, np, keep);
        memcpy(b->shape_pos_cpu, pos, np * sizeof(md_float3_t));
        md_array_resize(b->shape_scale_cpu, ns, keep);
        if (ns) memcpy(b->shape_scale_cpu, scale, ns * sizeof(float));
        b->content[MD_RENDER_CONTENT_GEOMETRY] += 1;
    }
    if (app) {
        md_gpu_upload(b->stream, b->shape_col, col, np * sizeof(uint32_t));
        md_array_resize(b->shape_col_cpu, np, keep);
        memcpy(b->shape_col_cpu, col, np * sizeof(uint32_t));
        b->content[MD_RENDER_CONTENT_APPEARANCE] += 1;
    }

    for (size_t i = 0; i < md_array_size(sb); ++i) {
        const ShapeBatch& x = sb[i];
        md_render_batch_t bt = {};
        bt.key   = x.key;
        bt.kind  = MD_RENDER_CAPSULES;
        bt.count = x.count;
        bt.points.count    = 2 * x.count;
        bt.points.position = b->shape_pos + (uint64_t)x.first * sizeof(md_float3_t);
        bt.points.color    = b->shape_col + (uint64_t)x.first * sizeof(uint32_t);
        bt.color.mode      = MD_RENDER_COLOR_POINTS;
        bt.color.sharpness = 1.0f;
        bt.style.capsules  = x.style;
        if (x.scale >= 0) bt.style.capsules.radius_scale = b->shape_scale + (uint64_t)x.scale * sizeof(float);
        bt.transform       = x.transform;
        bt.transform_count = x.transform_count;
        md_array_push(b->batches, bt, alloc);
    }
}

// ---- The frame ---------------------------------------------------------------------------------------------------

// md_compose's shading of the settings', the compose pass's that was (compose_deferred: its roughness and F0)
static void shading_of(md_compose_shading_t* sh, const postprocess_pipeline::Settings& s) {
    const postprocess_pipeline::Lighting l = postprocess_pipeline::lighting(s.shading.mode);
    const float env = MD_MAX(s.shading.environment, 0.0f);
    *sh = {};
    sh->mode = s.shading.mode == postprocess_pipeline::ShadingMode_Basic ? MD_COMPOSE_SHADING_BASIC : MD_COMPOSE_SHADING_PBR;
    sh->linear_colors = s.shading.colors == postprocess_pipeline::ColorSpace_Linear;
    for (int k = 0; k < 3; ++k) {
        sh->light_dir[k]       = l.light_dir.elem[k];
        sh->dir_radiance[k]    = l.light_irradiance.elem[k];
        sh->sky_radiance[k]    = l.sky_radiance.elem[k] * env;
        sh->ground_radiance[k] = l.ground_radiance.elem[k] * env;
    }
    sh->light_strength = MD_MAX(s.shading.light, 0.0f);
    sh->roughness = 0.4f;
    sh->f0        = 0.04f;
    sh->ambient   = MD_CLAMP(s.shading.ambient, 0.0f, 1.0f);
    sh->specular  = MD_MAX(s.shading.specular, 0.0f);
}

// md_post's display transform of the settings' (display::transform of the post-processing that was): basic shading
// and no tone mapping are the encoding alone
static void display_of(md_post_display_t* d, const postprocess_pipeline::Settings& s) {
    d->curve = MD_POST_CURVE_NONE;
    if (s.shading.mode != postprocess_pipeline::ShadingMode_Basic && s.tonemap.enabled) {
        switch (s.tonemap.mode) {
        case postprocess_pipeline::Tonemapper_Neutral: d->curve = MD_POST_CURVE_NEUTRAL; break;
        case postprocess_pipeline::Tonemapper_ACES:    d->curve = MD_POST_CURVE_ACES;    break;
        case postprocess_pipeline::Tonemapper_Filmic:  d->curve = MD_POST_CURVE_HABLE;   break;
        default: break;
        }
    }
    d->exposure = s.tonemap.exposure;
    d->gamma    = s.tonemap.gamma;
    for (int k = 0; k < 3; ++k) d->background[k] = s.background_color.elem[k];
}

void md_render_bridge_render(ApplicationState* state, const postprocess_pipeline::Settings& settings) {
    MdRenderBridge* b = state->md_render_bridge;
    if (b) b->still = true;
    if (!b || !b->ok) return;
    b->has_last = false;
    md_array_shrink(b->batches, 0);
    const uint32_t width = state->gbuffer.width, height = state->gbuffer.height;
    const uint32_t out_width = state->gbuffer.out_width ? state->gbuffer.out_width : width;
    const uint32_t out_height = state->gbuffer.out_height ? state->gbuffer.out_height : height;
    if (width == 0 || height == 0) return;
    if (!ensure_targets(b, width, height, out_width, out_height)) {
        b->ok = false;
        return;
    }

    // The system's representations and shapes, when there is one (else the background alone)
    if (ensure_system(b, state) && take_up(b, state)) {
    // Forget the representations that are no more
    for (size_t i = 0; i < md_array_size(b->reps); ++i) b->reps[i].seen = false;
    const md_render_mol_atoms_t atoms = atoms_of(b, state);
    const uint64_t secondary = update_secondary(b, state);
    update_parts(b, state);
    md_allocator_i* alloc = state->allocator.persistent;     // the batches' array is kept

    md_render_mol_bond_style_t styles[64][MD_RENDER_MOL_BOND_KIND_COUNT];
    size_t num_styles = 0;
    const float dim = md_bitmask_empty(&state->selection.selection_mask) ? 0.0f : 1.0f - MD_CLAMP(state->selection.color.saturation, 0.0f, 1.0f);
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        const Representation& rep = state->representation.reps[i];
        BridgeRep* r = find_rep(b, rep.node.id);
        if (r) r->seen = true;
        if (!draws(*state, rep) || !r) continue;

        md_render_mol_rep_t mr = {};
        mr.key   = ((uint64_t)rep.node.id + 1) << 16;   // and the part of an assembly in the low bits (emit_placed)
        mr.color.mode = MD_RENDER_COLOR_POINTS;
        // What is not selected, dimmed while something is (to the selection colour's saturation)
        mr.color.desaturate_unmarked = dim;
        mr.color.sharpness = rep.bond_color == BondColorMode::SmoothAtom ? 1.0f - MD_CLAMP(rep.bond_sharpness, 0.0f, 1.0f) : 1.0f;
        mr.atom_color = r->color;
        mr.visible = r->visible;
        mr.chain_gaps = rep.chain_gaps;
        // While the system's bonds are being inferred, every type as spheres (representation_draw_op)
        const RepresentationType type = rep.sphere_fallback ? RepresentationType::SpaceFill : rep.type;
        switch (type) {
        case RepresentationType::SpaceFill:
            mr.type = MD_RENDER_MOL_REP_SPACEFILL;
            mr.params.spacefill.radius_scale = rep.sphere_fallback ? 1.0f : rep.scale.x;
            break;
        case RepresentationType::Licorice:
            mr.type = MD_RENDER_MOL_REP_LICORICE;
            mr.params.licorice.radius = 0.2f * rep.scale.x;                  // md_gl's scales
            mr.bond_orders = rep.bond_orders;
            break;
        case RepresentationType::BallAndStick:
            mr.type = MD_RENDER_MOL_REP_BALL_AND_STICK;
            mr.params.ball_and_stick.ball_scale = 0.2f * rep.scale.x;
            mr.params.ball_and_stick.stick_radius = 0.2f * rep.scale.y;
            mr.bond_orders = rep.bond_orders;
            break;
        case RepresentationType::Ribbons:
            mr.type = MD_RENDER_MOL_REP_RIBBONS;
            mr.params.ribbons.half_width = rep.scale.x;
            mr.params.ribbons.half_thickness = 0.1f * rep.scale.y;
            break;
        case RepresentationType::Cartoon: {
            mr.type = MD_RENDER_MOL_REP_CARTOON;
            const float w = rep.scale.w > 0.0f ? rep.scale.w : 1.0f;
            mr.params.cartoon.profile = cartoon_profile(b, r, rep, secondary);
            mr.params.cartoon.band_radius = MD_MAX(rep.nucleic_band_radius, 0.01f) * w;
            mr.params.cartoon.ring_radius = MD_MAX(rep.nucleic_ring_radius, 0.01f) * w;
            mr.params.cartoon.ring_plate  = 0.6f;
            mr.params.cartoon.lipid_scale = w;
            mr.params.cartoon.split_heads = rep.two_head_spheres;
            break;
        }
        case RepresentationType::Bonds:
            // The pairs its filter picks (md_render_bridge_rep_pairs), as thick as a licorice's sticks of the same scale
            if (pair_count(r->pairs) == 0) continue;
            mr.type = MD_RENDER_MOL_REP_PAIRS;
            mr.params.pairs.radius = 0.2f * rep.scale.x;
            mr.pairs = r->pairs.pairs;
            mr.pair_id = r->pairs.ids;
            mr.pair_offset = r->pairs.offset;
            break;
        default:
            continue;
        }
        if (type == RepresentationType::Licorice || type == RepresentationType::BallAndStick || type == RepresentationType::Bonds) {
            if (num_styles == MD_COUNTOF(styles)) continue;
            if (type == RepresentationType::Bonds) {
                pair_styles(styles[num_styles], rep);
            } else {
                bond_styles(styles[num_styles], rep);
                mr.bond_kinds = bond_kinds(rep);
            }
            mr.bond_style = styles[num_styles++];
        }
        emit_placed(b, state, atoms, mr, rep, alloc);
    }
    for (size_t i = 0; i < md_array_size(b->reps);) {
        if (!b->reps[i].seen) {
            free_rep(b->stream, &b->reps[i]);
            md_array_swap_back_and_pop(b->reps, i);
        } else {
            ++i;
        }
    }
    emit_shapes(b, state, alloc);
    }

    md_render_frame_t f = {};
    f.width  = width;
    f.height = height;
    camera_of(&f.camera, f.jitter_offset, *state, width, height);
    f.jitter = MD_RENDER_JITTER_GIVEN;
    f.time   = b->time;
    f.batches = b->batches;
    f.batch_count = (uint32_t)md_array_size(b->batches);
    for (int i = 0; i < MD_RENDER_CONTENT_COUNT; ++i) f.content[i] = b->content[i];
    md_render_status_t st = {};
    if (!md_render_view_draw(b->stream, b->view, &f, &st)) {
        set_status(b, "md_render could not draw the frame (see the log)%s");
        return;
    }
    b->still = (st.work & MD_RENDER_WORK_RASTER) == 0;
    b->last = f;
    b->has_last = true;
    b->last_ui_frame = ImGui::GetFrameCount();

    // Composed, the display transform, the TAA into the history at the image's size (the frame alone when there is
    // none to go on from)
    md_compose_desc_t cd = {};
    cd.camera = f.camera;
    shading_of(&cd.shading, settings);
    cd.ao = settings.ssao.enabled ? MD_MAX(settings.ssao.intensity, 0.0f) : 0.0f;
    cd.ao_mode = settings.ssao.mode == postprocess_pipeline::SsaoMode_Quality ? MD_COMPOSE_AO_QUALITY : MD_COMPOSE_AO_PERFORMANCE;
    cd.frame  = (uint32_t)ImGui::GetFrameCount();
    cd.target = b->scene;
    md_post_display_t dd = {};
    display_of(&dd, settings);
    dd.scene  = b->scene;
    dd.target = b->frame;
    bool ok = md_compose(b->stream, b->compose, b->view, &cd) && md_post_display(b->stream, b->post, &dd);
    if (ok && settings.taa.enabled) {
        md_post_resolve_t rd = {};
        rd.view = b->view;
        rd.camera = f.camera;
        rd.jitter[0] = st.jitter[0];
        rd.jitter[1] = st.jitter[1];
        rd.frame = b->frame;
        rd.history_prev = b->hist_valid ? b->history[b->hist] : nullptr;
        rd.history = b->history[b->hist ^ 1];
        rd.feedback_min = settings.taa.feedback_min;
        rd.feedback_max = settings.taa.feedback_max;
        ok = md_post_resolve(b->stream, b->post, &rd);
        b->hist ^= 1;
        b->hist_valid = ok;
        b->image = ok ? b->history[b->hist] : nullptr;
    } else {
        b->hist_valid = false;
        b->image = ok ? b->frame : nullptr;
    }
    if (!ok) {
        set_status(b, "The main view's frame could not be made (see the log)%s");
        return;
    }

    // The depth, for GL to draw behind the scene (depth_gl)
    share_write_begin(b, &b->depth);
    md_render_export_t ex = {b->depth.tex, nullptr, nullptr};
    b->exported = md_render_view_export(b->stream, b->view, &ex);
    share_write_end(b, &b->depth);
    b->depth_tex_valid = false;
    snprintf(b->status, sizeof(b->status), "Drawing %u batches, %u x %u, shown at %u x %u", f.batch_count, width, height, out_width, out_height);
}

bool md_render_bridge_present(ApplicationState* state, const postprocess_pipeline::Settings* marks) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!b || !b->ok || !b->image || !b->out.tex) return false;
    // The marks changed without the scene being made: md_render's masks again (only they are made)
    if (marks && b->has_last && b->last_ui_frame != ImGui::GetFrameCount() && b->last.width == b->width && b->last.height == b->height) {
        if (take_up(b, state)) {
            md_render_frame_t f = b->last;
            for (int i = 0; i < MD_RENDER_CONTENT_COUNT; ++i) f.content[i] = b->content[i];
            if (md_render_view_draw(b->stream, b->view, &f, nullptr)) {
                b->last = f;
                b->last_ui_frame = ImGui::GetFrameCount();
            }
        }
    }
    md_post_marks_t m = {};
    if (marks) {
        const auto& mk = marks->marks;
        for (int i = 0; i < 2; ++i) {
            for (int k = 0; k < 4; ++k) {
                m.visible[i][k] = mk.visible[i].elem[k];
                m.hidden[i][k]  = mk.hidden[i].elem[k];
                m.outline[i][k] = mk.outline[i].elem[k];
            }
        }
        for (int g = 0; g < 7; ++g) for (int k = 0; k < 3; ++k) m.groups[g][k] = mk.groups[g].elem[k];
        m.highlight  = mk.highlight;
        m.outline_px = (uint32_t)MD_MAX(mk.outline_step, 1);
    }
    md_post_present_t pd = {};
    pd.image = b->image;
    pd.marks_view = marks && b->has_last ? b->view : nullptr;
    pd.marks  = marks ? &m : nullptr;
    pd.target = b->out.tex;
    pd.frame  = b->presents++;
    share_write_begin(b, &b->out);
    const bool ok = md_post_present(b->stream, b->post, &pd);
    share_write_end(b, &b->out);
    if (!ok) return false;

    // Into the framebuffer bound, over the viewport the caller set
    gl_read_begin(&b->out);
    const GLboolean depth_test = glIsEnabled(GL_DEPTH_TEST), blend = glIsEnabled(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glUseProgram(b->present_program);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(b->out.gl.target, b->out.gl.name);
    glUniform1i(b->ploc_image, 0);
    glUniform2i(b->ploc_size, (GLint)b->out_width, (GLint)b->out_height);
    glBindVertexArray(b->vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glBindTexture(b->out.gl.target, 0);
    glUseProgram(0);
    if (depth_test) glEnable(GL_DEPTH_TEST);
    if (blend) glEnable(GL_BLEND);
    gl_read_end(&b->out);
    return true;
}

uint32_t md_render_bridge_depth_gl(ApplicationState* state) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!b || !b->ok || !b->exported || !b->depth_program) return 0;
    if (b->depth_tex_valid) return b->depth_tex;
    if (!b->depth_tex || b->depth_tex_w != b->width || b->depth_tex_h != b->height) {
        if (b->depth_tex) glDeleteTextures(1, &b->depth_tex);
        glGenTextures(1, &b->depth_tex);
        glBindTexture(GL_TEXTURE_2D, b->depth_tex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, (GLsizei)b->width, (GLsizei)b->height, 0, GL_RED, GL_FLOAT, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glBindTexture(GL_TEXTURE_2D, 0);
        b->depth_tex_w = b->width;
        b->depth_tex_h = b->height;
    }
    GLint fbo = 0, viewport[4];
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &fbo);
    glGetIntegerv(GL_VIEWPORT, viewport);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, b->depth_fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, b->depth_tex, 0);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glViewport(0, 0, (GLsizei)b->width, (GLsizei)b->height);
    const GLboolean depth_test = glIsEnabled(GL_DEPTH_TEST), blend = glIsEnabled(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    gl_read_begin(&b->depth);
    glUseProgram(b->depth_program);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(b->depth.gl.target, b->depth.gl.name);
    glUniform1i(b->dloc_depth, 0);
    glUniform2i(b->dloc_size, (GLint)b->width, (GLint)b->height);
    glBindVertexArray(b->vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glBindTexture(b->depth.gl.target, 0);
    glUseProgram(0);
    gl_read_end(&b->depth);
    if (depth_test) glEnable(GL_DEPTH_TEST);
    if (blend) glEnable(GL_BLEND);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, (GLuint)fbo);
    glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
    b->depth_tex_valid = true;
    return b->depth_tex;
}

// ---- Picking ------------------------------------------------------------------------------------------------------

static void collect_picks(MdRenderBridge* b) {
    md_render_pick_t p;
    while (md_render_view_pick_result(b->view, &p)) {
        if (b->num_picks == MdRenderBridge::PICKS) {
            memmove(b->picks, b->picks + 1, sizeof(b->picks[0]) * (MdRenderBridge::PICKS - 1));
            b->num_picks -= 1;
        }
        b->picks[b->num_picks++] = p;
    }
}

void md_render_bridge_pick(ApplicationState* state, float x, float y) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!b || !b->ok || !b->has_last) return;
    collect_picks(b);
    const int px = (int)x, py = (int)y;
    if (px < 0 || py < 0 || px >= (int)b->width || py >= (int)b->height) return;
    md_render_view_pick(b->stream, b->view, (uint32_t)px, (uint32_t)py);
}

// The hit as the domain's: its index in the domain's range of the frame
static void set_hit(const ApplicationState& state, PickingHit* hit, PickingDomainID domain, uint64_t key, uint32_t local) {
    const PickingSpace* space = picking_handler_find_space(state.picking_handler, hit->frame_idx);
    const PickingRange* range = space ? picking_space_find_range(*space, domain, key) : nullptr;
    hit->domain    = domain;
    hit->key       = key;
    hit->local_idx = local;
    hit->raw_idx   = range && local < range->end - range->beg ? range->beg + local : INVALID_PICKING_IDX;
}

// Whether a hit on the stick between a pair of atoms is in its middle half (md_gl's bond.frag: the bond), not a quarter at
// an end (its atom), with the pair's ends where the hit's copy has them
static bool pair_middle(const MdRenderBridge* b, const ApplicationState& state, const md_render_mol_hit_t& mh, md_atom_pair_t pair) {
    const md_float3_t* xyz = state.mold.state.xyz;
    const md_float3_t p = {mh.position[0], mh.position[1], mh.position[2]};
    md_float3_t a = xyz[pair.idx[0]], e = xyz[pair.idx[1]];
    md_float4x4_t T;
    if (copy_transform(b, state, mh.rep_key, mh.instance, &T)) {
        a = md_float3_from_float4(md_float4x4_mul_float4(T, md_float4_from_float3(a, 1.0f)));
        e = md_float3_from_float4(md_float4x4_mul_float4(T, md_float4_from_float3(e, 1.0f)));
    }
    const md_float3_t ab = e - a;
    const float len2 = md_float3_dot(ab, ab);
    const float t = len2 > 0.0f ? md_float3_dot(p - a, ab) / len2 : 0.5f;
    return t >= 0.25f && t <= 0.75f;
}

bool md_render_bridge_hit(ApplicationState* state, PickingHit* hit) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!hit || !b || !b->ok) return false;
    collect_picks(b);
    if (b->num_picks == 0) return false;
    // The newest answer, of the pixel asked about a frame or two ago (as the GL picking's readback was)
    const md_render_pick_t* found = &b->picks[b->num_picks - 1];
    *hit = {};
    hit->frame_idx     = state->picking_handler.frame_idx;
    hit->raw_idx       = INVALID_PICKING_IDX;
    // Where, as the GL picking had it: the pixel's centre, rows from the bottom, and GL's window depth (md_render's
    // is reversed)
    hit->surface_coord = {(float)found->x + 0.5f, (float)b->height - (float)found->y - 0.5f};
    hit->depth         = found->key ? 1.0f - found->depth : 1.0f;
    if (!found->key) return true;
    // A dipole's arrow picks as md_gl's (its group's range, keyed by its vector attribute); the cells' edges as none
    if (found && is_shape_key(found->key)) {
        const uint32_t what = (uint32_t)((found->key >> 8) & 0xFF);
        if ((what == SHAPE_SHAFTS || what == SHAPE_HEADS) && found->element < md_array_size(b->dipoles)) {
            const MdRenderBridge::Dipole& d = b->dipoles[found->element];
            set_hit(*state, hit, PickingDomain_Dipole, d.key, d.index);
        }
        return true;     // the cells' edges as nothing to pick
    }
    md_render_mol_hit_t mh = {};
    // Of a system no more (or not known): nothing, rather than what md_render drew there taken for something else
    if (!b->has_topo || !md_render_mol_decode_pick(found, &mh) || mh.kind == MD_RENDER_MOL_HIT_NONE) return true;
    const md_system_t& sys = state->mold.sys;
    switch (mh.kind) {
    case MD_RENDER_MOL_HIT_BOND: {
        // As md_gl's: the middle half of a bond picks as the bond, the quarters at its ends as their atoms
        if (!pair_middle(b, *state, mh, sys.bond.pairs[mh.bond])) {
            set_hit(*state, hit, PickingDomain_Atom, 0, mh.atom);
        } else if (mh.bond < sys.bond.count) {
            set_hit(*state, hit, PickingDomain_Bond, 0, mh.bond);
        } else {
            // An interaction over the sticks of a licorice: of the representation's range of them
            set_hit(*state, hit, PickingDomain_Bond, (mh.rep_key >> 16) - 1, (uint32_t)(mh.bond - sys.bond.count));
        }
        break;
    }
    case MD_RENDER_MOL_HIT_PAIR: {
        // As md_gl's, the middle half the bond, the quarters at its ends their atoms. A pair of a Bonds representation
        // is a weak bond of its range (its index in bond_pairs); otherwise one of the frame's own bonds (its sticks,
        // multiple sticks or partner): md_gl's index, of the system's or past it of one formed (update_frame_bonds)
        const uint32_t rep_id = (uint32_t)((mh.rep_key >> 16) - 1);
        const Representation* rep = nullptr;
        for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
            if (state->representation.reps[i].node.id == rep_id) rep = &state->representation.reps[i];
        }
        const md_bond_delta_t* delta = bond_tracking_delta(state);
        md_atom_pair_t pair = {};
        uint64_t key = 0;
        uint32_t local = mh.bond;
        if (rep && rep->type == RepresentationType::Bonds && mh.bond < md_array_size(rep->bond_pairs)) {
            pair = rep->bond_pairs[mh.bond];
            key = rep_id;
        } else if (rep && rep->type != RepresentationType::Bonds && mh.bond < sys.bond.count) {
            pair = sys.bond.pairs[mh.bond];
        } else if (rep && rep->type != RepresentationType::Bonds && delta && mh.bond - sys.bond.count < delta->num_formed) {
            pair = delta->formed_pair[mh.bond - sys.bond.count];
        } else {
            set_hit(*state, hit, PickingDomain_Atom, 0, mh.atom);
            break;
        }
        if (pair_middle(b, *state, mh, pair)) {
            set_hit(*state, hit, PickingDomain_Bond, key, local);
        } else {
            set_hit(*state, hit, PickingDomain_Atom, 0, mh.atom);
        }
        break;
    }
    case MD_RENDER_MOL_HIT_BACKBONE:
        set_hit(*state, hit, PickingDomain_BackboneSegment, 0, mh.segment);
        break;
    case MD_RENDER_MOL_HIT_RING:
        set_hit(*state, hit, PickingDomain_Ring, 0, mh.ring);
        break;
    case MD_RENDER_MOL_HIT_CHAIN_GAP:
        set_hit(*state, hit, PickingDomain_ChainGap, 0, mh.residue);
        break;
    default:
        // An atom, or what stands for one: a lipid's sphere or tube, a nucleotide's part, a residue's sphere
        set_hit(*state, hit, PickingDomain_Atom, 0, mh.atom);
        break;
    }
    return true;
}

void md_render_bridge_buffers(ApplicationState* state, uint32_t dirty) {
    MdRenderBridge* b = state->md_render_bridge;
    if (b) b->pending |= dirty;
}

// ---- Life --------------------------------------------------------------------------------------------------------

static bool create(ApplicationState* state) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!state->gpu_device) {
        set_status(b, "No GPU device (md_gpu)%s");
        return false;
    }
    md_gpu_device_info_t info = {};
    md_gpu_device_info(state->gpu_device, &info);
    if (!info.supports_graphics) {
        set_status(b, "The GPU device does not draw (md_gpu graphics)%s");
        return false;
    }
    if (!info.supports_sharing) {
        set_status(b, "The GPU device cannot share textures with OpenGL%s");
        return false;
    }
    b->device = state->gpu_device;
    b->stream = md_gpu_stream_default(b->device, MD_GPU_STREAM_GRAPHICS);
    b->render = b->stream ? md_render_create(b->device, nullptr) : nullptr;
    b->mol    = b->render ? md_render_mol_create(b->device, nullptr) : nullptr;
    b->view   = b->mol ? md_render_view_create(b->render, "main view") : nullptr;
    b->compose = b->view ? md_compose_create(b->device) : nullptr;
    b->post    = b->compose ? md_post_create(b->device) : nullptr;
    if (!b->post) {
        set_status(b, "md_render, md_compose or md_post could not be made: %s", md_gpu_last_error());
        return false;
    }
    // Whether GL can read what md_gpu writes: the targets made once (and GL's shaders compiled)
    if (!ensure_targets(b, 16, 16, 16, 16)) return false;
    set_status(b, "Ready%s");
    return true;
}

bool md_render_bridge_create(ApplicationState* state) {
    if (!state->md_render_bridge) {
        state->md_render_bridge = new MdRenderBridge();
        state->md_render_bridge->ok = create(state);
    }
    MdRenderBridge* b = state->md_render_bridge;
    if (!b->ok) return false;
    b->pending |= MolBit_DirtyPosition | MolBit_DirtyFlags;
    b->updated = false;
    flag_all_representations_as_dirty(state);
    return true;
}

void md_render_bridge_destroy(ApplicationState* state) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!b) return;
    if (b->stream) {
        free_system(b);
        if (b->shape_pos)   md_gpu_free(b->stream, b->shape_pos);
        if (b->shape_col)   md_gpu_free(b->stream, b->shape_col);
        if (b->shape_scale) md_gpu_free(b->stream, b->shape_scale);
        free_targets(b);
        if (b->view) md_render_view_destroy(b->stream, b->view);
        md_gpu_stream_sync(b->stream);
    }
    if (b->post) md_post_destroy(b->post);
    if (b->compose) md_compose_destroy(b->compose);
    if (b->mol) md_render_mol_destroy(b->mol);
    if (b->render) md_render_destroy(b->render);
    if (b->present_program) glDeleteProgram(b->present_program);
    if (b->depth_program) glDeleteProgram(b->depth_program);
    if (b->vao) glDeleteVertexArrays(1, &b->vao);
    if (b->depth_fbo) glDeleteFramebuffers(1, &b->depth_fbo);
    if (b->depth_tex) glDeleteTextures(1, &b->depth_tex);
    md_array_free(b->reps, state->allocator.persistent);
    md_array_free(b->batches, state->allocator.persistent);
    md_array_free(b->parts, state->allocator.persistent);
    md_array_free(b->shape_pos_cpu, state->allocator.persistent);
    md_array_free(b->shape_col_cpu, state->allocator.persistent);
    md_array_free(b->shape_scale_cpu, state->allocator.persistent);
    md_array_free(b->dipoles, state->allocator.persistent);
    delete b;
    state->md_render_bridge = nullptr;
}

const char* md_render_bridge_status(const ApplicationState& state) {
    const MdRenderBridge* b = state.md_render_bridge;
    if (!b) return "Not made";
    return b->status[0] ? b->status : (b->ok ? "Ready" : "Failed");
}

#else

bool md_render_bridge_create(ApplicationState*) { return false; }
void md_render_bridge_destroy(ApplicationState*) {}
const char* md_render_bridge_status(const ApplicationState&) { return "Built without md_gpu (MD_ENABLE_GPU): the main view is not drawn"; }
void md_render_bridge_pick(ApplicationState*, float, float) {}
bool md_render_bridge_hit(ApplicationState*, PickingHit*) { return false; }
void md_render_bridge_rep_colors(ApplicationState*, const Representation&, const uint32_t*, size_t) {}
void md_render_bridge_rep_pairs(ApplicationState*, const Representation&, const uint32_t*, const uint32_t*, const uint32_t*, size_t) {}
void md_render_bridge_buffers(ApplicationState*, uint32_t) {}
void md_render_bridge_render(ApplicationState*, const postprocess_pipeline::Settings&) {}
bool md_render_bridge_present(ApplicationState*, const postprocess_pipeline::Settings*) { return false; }
uint32_t md_render_bridge_depth_gl(ApplicationState*) { return 0; }

#endif
