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
#include <core/md_gpu_gl.h>

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

// A texture md_render's targets are exported into, and the GL texture of it
struct SharedTex {
    md_gpu_texture_t    tex = nullptr;
    md_gpu_gl_texture_t gl  = {};
    bool                imported = false;
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

    // The parts of the assembly shown (place_ops): the atoms of each (0: all) and its copies' transforms
    struct Part { md_gpu_addr_t atoms; uint32_t c0, c1; };
    md_array(Part) parts = nullptr;
    int            parts_of = -1;           // the assembly, -1: none
    uint64_t       parts_topology = 0;

    // The export, shared with GL
    uint32_t  width = 0, height = 0;
    SharedTex depth, color, surface;
    bool      exported = false;     // the textures hold a frame
    // GL reads the textures between its begin_read and end_read, md_gpu writes them between acquire and release, in
    // turn: 0 nothing to read, 1 released by md_gpu (GL to begin), 2 GL reading (to end, before md_gpu acquires)
    int       gl_read = 0;
    bool      still = true;         // md_render did not rasterise this frame: nothing moved
    md_render_frame_t last = {};    // the frame drawn last (its batches are `batches`), for the marks alone
    bool      has_last = false;
    int       last_ui_frame = -1;   // the ImGui frame it was drawn in

    // The GL merges: into the G-buffer, into the marks
    GLuint program = 0;
    GLuint vao = 0;
    GLint  loc_depth = -1, loc_color = -1, loc_surface = -1, loc_size = -1, loc_dim = -1, loc_still = -1, loc_index = -1;
    GLuint marks_program = 0;
    GLint  mloc_depth = -1, mloc_color = -1, mloc_scene = -1, mloc_size = -1, mloc_out_size = -1;

    // md_render's answers to the picks, the newest last, for the GL picking's hits of its index to find
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

// The G-buffer's encodings, as the snippet md_gl's shaders write with (main.cpp, shader_output_snippet)
static const char* merge_vs = R"(
#version 330 core
void main() {
    vec2 p = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
    gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}
)";

static const char* merge_fs = R"(
#ifdef RECT
uniform sampler2DRect u_depth;
uniform sampler2DRect u_color;
uniform sampler2DRect u_surface;
#define FETCH(t, p) texelFetch(t, p)
#else
uniform sampler2D u_depth;
uniform sampler2D u_color;
uniform sampler2D u_surface;
#define FETCH(t, p) texelFetch(t, p, 0)
#endif
uniform ivec2 u_size;
uniform float u_dim;      // the saturation of what is not selected; 1 while nothing is
uniform float u_still;    // 1: md_render did not rasterise this frame, nothing of it moved
uniform vec4  u_index;    // the picking index of what md_render drew (PickingDomain_MdRender), as encode_index has it

layout(location = 0) out vec4 out_color;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec4 out_velocity;
layout(location = 3) out vec4 out_atom_index;

vec3 oct_decode(vec2 e) {
    vec3 n = vec3(e, 1.0 - abs(e.x) - abs(e.y));
    if (n.z < 0.0) n.xy = (1.0 - abs(n.yx)) * vec2(n.x >= 0.0 ? 1.0 : -1.0, n.y >= 0.0 ? 1.0 : -1.0);
    return normalize(n);
}

vec2 encode_normal(vec3 n) {
   float p = sqrt(n.z * 8.0 + 8.0);
   return n.xy / p + 0.5;
}

void main() {
    // md_render's rows from the top, GL's from the bottom
    ivec2 p = ivec2(gl_FragCoord.xy);
    p.y = u_size.y - 1 - p.y;
    float d = FETCH(u_depth, p).r;
    if (d <= 0.0) discard;
    vec4 c = FETCH(u_color, p);
    vec4 s = FETCH(u_surface, p);
    uint bits = uint(c.a * 255.0 + 0.5);
    vec3 rgb = c.rgb;
    if ((bits & 1u) == 0u) {
        rgb = mix(vec3(dot(rgb, vec3(0.2126, 0.7152, 0.0722))), rgb, u_dim);
    }
    // md_render's projection is md_gl's with its depth reversed: d = (1 - ndc z) / 2
    gl_FragDepth   = 1.0 - d;
    out_color      = vec4(rgb, 1.0);
    out_normal     = vec4(encode_normal(oct_decode(s.xy)), 0.0, 0.0);
    // Pixels, y down, to the G-buffer's: of the uv, y up
    out_velocity   = u_still > 0.5 ? vec4(0.0) : vec4(s.z / float(u_size.x), -s.w / float(u_size.y), 0.0, 0.0);
    out_atom_index = u_index;
}
)";

// The marks as md_gl's marks pass has them (main.cpp, shader_output_snippet_marks), blended by max into its targets:
// of the selected 1 in r in view and in g there at all, of the highlighted the same in b and a, and the highlight's
// group + 1, 8 more in view. md_render's visible surface is in view where it is the scene's (the G-buffer's depth is
// what the merge wrote there, or nearer). Of the size of the image presented: md_render's texel under each pixel.
static const char* marks_fs = R"(
#ifdef RECT
uniform sampler2DRect u_depth;
uniform sampler2DRect u_color;
#define FETCH(t, p) texelFetch(t, p)
#else
uniform sampler2D u_depth;
uniform sampler2D u_color;
#define FETCH(t, p) texelFetch(t, p, 0)
#endif
uniform sampler2D u_scene;    // the G-buffer's depth
uniform ivec2 u_size;         // md_render's, the G-buffer's
uniform vec2  u_out_size;     // the marks'

layout(location = 0) out vec4 out_marks;
layout(location = 1) out vec4 out_group;

void main() {
    ivec2 p = clamp(ivec2(gl_FragCoord.xy / u_out_size * vec2(u_size)), ivec2(0), u_size - 1);
    float scene = texelFetch(u_scene, p, 0).r;
    p.y = u_size.y - 1 - p.y;
    uint a = uint(FETCH(u_color, p).a * 255.0 + 0.5);
    uint bits = a & 15u;
    if (bits == 0u) discard;
    float d = FETCH(u_depth, p).r;
    bool in_view = d > 0.0 && 1.0 - d <= scene + 1e-5;
    float sel_view = (bits & 1u) != 0u && in_view ? 1.0 : 0.0;
    float hl_view  = (bits & 4u) != 0u && in_view ? 1.0 : 0.0;
    float sel_any  = (bits & 2u) != 0u ? 1.0 : 0.0;
    float hl_any   = (bits & 8u) != 0u ? 1.0 : 0.0;
    out_marks = vec4(sel_view, sel_any, hl_view, hl_any);
    uint group = 1u + (a >> 4u) + (hl_view > 0.0 ? 8u : 0u);
    out_group = vec4(hl_any > 0.0 ? float(group) / 255.0 : 0.0, 0.0, 0.0, 0.0);
}
)";

static GLuint gl_program(const char* fs_body, bool rect) {
    char fs_src[4096];
    snprintf(fs_src, sizeof(fs_src), "#version 330 core\n%s%s", rect ? "#define RECT 1\n" : "", fs_body);
    GLuint vs = gl::compile_shader_from_source(md_str_from_cstr(merge_vs), GL_VERTEX_SHADER);
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
    b->program = gl_program(merge_fs, rect);
    b->marks_program = gl_program(marks_fs, rect);
    if (!b->program || !b->marks_program) return false;
    b->loc_depth   = glGetUniformLocation(b->program, "u_depth");
    b->loc_color   = glGetUniformLocation(b->program, "u_color");
    b->loc_surface = glGetUniformLocation(b->program, "u_surface");
    b->loc_size    = glGetUniformLocation(b->program, "u_size");
    b->loc_dim     = glGetUniformLocation(b->program, "u_dim");
    b->loc_still   = glGetUniformLocation(b->program, "u_still");
    b->loc_index   = glGetUniformLocation(b->program, "u_index");
    b->mloc_depth    = glGetUniformLocation(b->marks_program, "u_depth");
    b->mloc_color    = glGetUniformLocation(b->marks_program, "u_color");
    b->mloc_scene    = glGetUniformLocation(b->marks_program, "u_scene");
    b->mloc_size     = glGetUniformLocation(b->marks_program, "u_size");
    b->mloc_out_size = glGetUniformLocation(b->marks_program, "u_out_size");
    glGenVertexArrays(1, &b->vao);
    return true;
}

// ---- Shared textures ---------------------------------------------------------------------------------------------

// GL begins to read what md_gpu released
static void begin_reads(MdRenderBridge* b) {
    if (b->gl_read != 1) return;
    SharedTex* t[3] = {&b->depth, &b->color, &b->surface};
    for (SharedTex* s : t) md_gpu_gl_begin_read(&s->gl);
    b->gl_read = 2;
}

// GL is done with what was exported (read or not): the two take turns with each texture, and md_gpu writes it again
// only after GL said so. Before md_gpu acquires them again.
static void settle_reads(MdRenderBridge* b) {
    begin_reads(b);
    if (b->gl_read != 2) return;
    SharedTex* t[3] = {&b->depth, &b->color, &b->surface};
    for (SharedTex* s : t) md_gpu_gl_end_read(&s->gl);
    b->gl_read = 0;
}

static void shared_free(MdRenderBridge* b, SharedTex* s) {
    if (s->imported) md_gpu_gl_release(&s->gl);
    if (s->tex) md_gpu_texture_destroy(s->tex);
    *s = SharedTex{};
    (void)b;
}

static bool shared_make(MdRenderBridge* b, SharedTex* s, md_gpu_format_t format, const char* label) {
    md_gpu_texture_desc_t d = {};
    d.type   = MD_GPU_TEX_2D;
    d.format = format;
    d.usage  = MD_GPU_TEX_STORAGE | MD_GPU_TEX_SHARED;
    d.width  = b->width;
    d.height = b->height;
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

static bool shared_resize(MdRenderBridge* b, uint32_t width, uint32_t height) {
    if (b->depth.tex && b->width == width && b->height == height) return true;
    settle_reads(b);
    shared_free(b, &b->depth);
    shared_free(b, &b->color);
    shared_free(b, &b->surface);
    b->exported = false;
    b->width = width;
    b->height = height;
    if (!shared_make(b, &b->depth, MD_GPU_FORMAT_R32_FLOAT, "md_render export depth") ||
        !shared_make(b, &b->color, MD_GPU_FORMAT_RGBA8_UNORM, "md_render export color") ||
        !shared_make(b, &b->surface, MD_GPU_FORMAT_RGBA16_FLOAT, "md_render export surface")) {
        return false;
    }
    if (!b->program && !gl_init(b, b->depth.gl.target != GL_TEXTURE_2D)) {
        set_status(b, "The merge shader did not compile%s");
        return false;
    }
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
    return state.md_render.enabled && b && b->ok && b->has_topo;
}

bool md_render_bridge_active(const ApplicationState& state) { return active(state); }

bool md_render_bridge_draws(const ApplicationState& state, const Representation& rep) {
    if (!active(state)) return false;
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
    if (!state->md_render.enabled || !b || !b->ok || !b->has_topo || rep.type != RepresentationType::Bonds) return;
    BridgeRep* r = find_or_add_rep(b, rep.node.id, state->allocator.persistent);
    // Each pair's index the one it picks as (representation_weak_bond_draw_ops)
    upload_pairs(b, &r->pairs, atom_a, atom_b, flags, count, state->allocator.frame);
    b->content[MD_RENDER_CONTENT_GEOMETRY] += 1;
}

void md_render_bridge_rep_colors(ApplicationState* state, const Representation& rep, const uint32_t* colors, size_t count) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!state->md_render.enabled || !b || !b->ok || !type_drawn(rep.type)) return;
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

// ---- The frame ---------------------------------------------------------------------------------------------------

// The view's targets into the shared textures: md_gpu acquires them (after GL's last read, settle_reads), writes and
// releases them to GL
static void export_frame(MdRenderBridge* b) {
    SharedTex* t[3] = {&b->depth, &b->color, &b->surface};
    for (SharedTex* s : t) md_gpu_texture_share_acquire(b->stream, s->tex);
    md_render_export_t ex = {b->depth.tex, b->color.tex, b->surface.tex};
    b->exported = md_render_view_export(b->stream, b->view, &ex);
    for (SharedTex* s : t) md_gpu_texture_share_release(b->stream, s->tex);
    b->gl_read = 1;
}

void md_render_bridge_render(ApplicationState* state) {
    MdRenderBridge* b = state->md_render_bridge;
    if (b) b->still = true;
    if (!state->md_render.enabled || !b || !b->ok) return;
    settle_reads(b);
    b->has_last = false;
    md_array_shrink(b->batches, 0);

    if (!ensure_system(b, state)) return;
    if (!take_up(b, state)) return;

    // Forget the representations that are no more
    for (size_t i = 0; i < md_array_size(b->reps); ++i) b->reps[i].seen = false;
    const md_render_mol_atoms_t atoms = atoms_of(b, state);
    const uint64_t secondary = update_secondary(b, state);
    update_parts(b, state);
    md_allocator_i* alloc = state->allocator.persistent;     // the batches' array is kept

    md_render_mol_bond_style_t styles[64][MD_RENDER_MOL_BOND_KIND_COUNT];
    size_t num_styles = 0;
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        const Representation& rep = state->representation.reps[i];
        BridgeRep* r = find_rep(b, rep.node.id);
        if (r) r->seen = true;
        if (!md_render_bridge_draws(*state, rep) || !r) continue;

        md_render_mol_rep_t mr = {};
        mr.key   = ((uint64_t)rep.node.id + 1) << 16;   // and the part of an assembly in the low bits (emit_placed)
        mr.color.mode = MD_RENDER_COLOR_POINTS;
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

    const uint32_t width = state->gbuffer.width, height = state->gbuffer.height;
    if (width == 0 || height == 0 || md_array_size(b->batches) == 0) {
        b->exported = false;
        return;
    }
    if (!shared_resize(b, width, height)) {
        b->ok = false;
        return;
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
        b->exported = false;
        return;
    }
    b->still = (st.work & MD_RENDER_WORK_RASTER) == 0;
    b->last = f;
    b->has_last = true;
    b->last_ui_frame = ImGui::GetFrameCount();
    export_frame(b);
    snprintf(b->status, sizeof(b->status), "Drawing %u batches, %u x %u", f.batch_count, width, height);
}

// The marks of the frame drawn last into md_gl's marks pass (render, after draw_representations_marks): drawn by md_render
// again first when the scene was not made this frame (only its masks: the marks changed)
void md_render_bridge_marks(ApplicationState* state) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!active(*state) || !b->has_last || !b->marks_program || b->last.width != b->width || b->last.height != b->height) return;
    if (b->last_ui_frame != ImGui::GetFrameCount()) {
        settle_reads(b);
        if (!take_up(b, state)) return;
        md_render_frame_t f = b->last;
        for (int i = 0; i < MD_RENDER_CONTENT_COUNT; ++i) f.content[i] = b->content[i];
        if (!md_render_view_draw(b->stream, b->view, &f, nullptr)) return;
        b->last = f;
        b->last_ui_frame = ImGui::GetFrameCount();
        export_frame(b);
    }
    if (!b->exported) return;
    begin_reads(b);

    // Into the marks pass's targets, by max as it set (gbuffer_marks_begin)
    glUseProgram(b->marks_program);
    const SharedTex* t[2] = {&b->depth, &b->color};
    const GLint locs[2] = {b->mloc_depth, b->mloc_color};
    for (int i = 0; i < 2; ++i) {
        glActiveTexture(GL_TEXTURE0 + i);
        glBindTexture(t[i]->gl.target, t[i]->gl.name);
        glUniform1i(locs[i], i);
    }
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, state->gbuffer.tex.depth);
    glUniform1i(b->mloc_scene, 2);
    glUniform2i(b->mloc_size, (GLint)b->width, (GLint)b->height);
    glUniform2f(b->mloc_out_size, (float)state->gbuffer.out_width, (float)state->gbuffer.out_height);
    GLboolean cull = glIsEnabled(GL_CULL_FACE);
    glDisable(GL_CULL_FACE);
    glBindVertexArray(b->vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    if (cull) glEnable(GL_CULL_FACE);
    for (int i = 2; i >= 0; --i) {
        glActiveTexture(GL_TEXTURE0 + i);
        glBindTexture(i < 2 ? t[i]->gl.target : GL_TEXTURE_2D, 0);
    }
    glUseProgram(0);
}

void md_render_bridge_merge(ApplicationState* state) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!state->md_render.enabled || !b || !b->ok || !b->program) return;
    if (!b->exported) return;
    // Read until md_gpu next writes them (settle_reads): the marks read them too
    begin_reads(b);
    SharedTex* t[3] = {&b->depth, &b->color, &b->surface};

    const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY, GL_COLOR_ATTACHMENT_PICKING};
    glDrawBuffers((GLsizei)MD_COUNTOF(draw_buffers), draw_buffers);
    glViewport(0, 0, (GLsizei)b->width, (GLsizei)b->height);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glDepthMask(GL_TRUE);
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glDisable(GL_BLEND);

    glUseProgram(b->program);
    const GLint units[3] = {0, 1, 2};
    const GLint locs[3] = {b->loc_depth, b->loc_color, b->loc_surface};
    for (int i = 0; i < 3; ++i) {
        glActiveTexture(GL_TEXTURE0 + units[i]);
        glBindTexture(t[i]->gl.target, t[i]->gl.name);
        glUniform1i(locs[i], units[i]);
    }
    const bool selection = !md_bitmask_empty(&state->selection.selection_mask);
    glUniform2i(b->loc_size, (GLint)b->width, (GLint)b->height);
    glUniform1f(b->loc_dim, selection ? MD_CLAMP(state->selection.color.saturation, 0.0f, 1.0f) : 1.0f);
    glUniform1f(b->loc_still, b->still ? 1.0f : 0.0f);
    {
        // Picked as what md_render found there (md_render_bridge_resolve_hit)
        const PickingSpace* space = picking_handler_current_space(&state->picking_handler);
        const PickingRange* range = space ? picking_space_find_range(*space, PickingDomain_MdRender, 0) : nullptr;
        const uint32_t idx = range ? range->beg : INVALID_PICKING_IDX;
        glUniform4f(b->loc_index, (float)(idx & 0xFF) / 255.0f, (float)((idx >> 8) & 0xFF) / 255.0f, (float)((idx >> 16) & 0xFF) / 255.0f, (float)(idx >> 24) / 255.0f);
    }
    glBindVertexArray(b->vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    for (int i = 2; i >= 0; --i) {
        glActiveTexture(GL_TEXTURE0 + units[i]);
        glBindTexture(t[i]->gl.target, 0);
    }
    glUseProgram(0);
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
    if (!active(*state) || !b->exported) return;
    collect_picks(b);
    const int px = (int)x, py = (int)y;
    if (px < 0 || py < 0 || px >= (int)b->width || py >= (int)b->height) return;
    md_render_view_pick(b->stream, b->view, (uint32_t)px, b->height - 1 - (uint32_t)py);     // rows from the top
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

void md_render_bridge_resolve_hit(ApplicationState* state, PickingHit* hit) {
    if (!hit || hit->domain != PickingDomain_MdRender) return;
    MdRenderBridge* b = state->md_render_bridge;
    const md_render_pick_t* found = nullptr;
    if (active(*state)) {
        collect_picks(b);
        const uint32_t px = (uint32_t)hit->surface_coord.x, py = b->height - 1 - (uint32_t)hit->surface_coord.y;
        for (uint32_t i = b->num_picks; i-- > 0;) {
            if (b->picks[i].x == px && b->picks[i].y == py) { found = &b->picks[i]; break; }
        }
    }
    md_render_mol_hit_t mh = {};
    if (!found || !md_render_mol_decode_pick(found, &mh) || mh.kind == MD_RENDER_MOL_HIT_NONE) {
        // Not known yet (or no more): nothing, rather than what md_render drew there taken for something else
        hit->domain = 0;
        hit->key = 0;
        hit->local_idx = 0;
        hit->raw_idx = INVALID_PICKING_IDX;
        return;
    }
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
    if (!b->view) {
        set_status(b, "md_render could not be made: %s", md_gpu_last_error());
        return false;
    }
    // Whether GL can read what md_gpu writes: a texture made and imported once (and the merge's shader compiled)
    if (!shared_resize(b, 16, 16)) return false;
    set_status(b, "Ready%s");
    return true;
}

bool md_render_bridge_enable(ApplicationState* state, bool enable) {
    state->md_render.enabled = enable;
    if (!enable) {
        // Read what was exported, so that the textures can be written again when it is turned on
        if (state->md_render_bridge) settle_reads(state->md_render_bridge);
        flag_all_representations_as_dirty(state);
        return true;
    }
    if (!state->md_render_bridge) {
        state->md_render_bridge = new MdRenderBridge();
        state->md_render_bridge->ok = create(state);
        if (!state->md_render_bridge->ok) {
            MD_LOG_INFO("md_render bridge: %s", state->md_render_bridge->status);
        }
    }
    MdRenderBridge* b = state->md_render_bridge;
    if (!b->ok) {
        state->md_render.enabled = false;
        return false;
    }
    b->pending |= MolBit_DirtyPosition | MolBit_DirtyFlags;
    b->updated = false;
    flag_all_representations_as_dirty(state);
    return true;
}

void md_render_bridge_destroy(ApplicationState* state) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!b) return;
    settle_reads(b);
    if (b->stream) {
        free_system(b);
        shared_free(b, &b->depth);
        shared_free(b, &b->color);
        shared_free(b, &b->surface);
        if (b->view) md_render_view_destroy(b->stream, b->view);
        md_gpu_stream_sync(b->stream);
    }
    if (b->mol) md_render_mol_destroy(b->mol);
    if (b->render) md_render_destroy(b->render);
    if (b->program) glDeleteProgram(b->program);
    if (b->vao) glDeleteVertexArrays(1, &b->vao);
    md_array_free(b->reps, state->allocator.persistent);
    md_array_free(b->batches, state->allocator.persistent);
    md_array_free(b->parts, state->allocator.persistent);
    delete b;
    state->md_render_bridge = nullptr;
}

const char* md_render_bridge_status(const ApplicationState& state) {
    const MdRenderBridge* b = state.md_render_bridge;
    if (!b) return "Off";
    return b->status[0] ? b->status : (state.md_render.enabled ? "On" : "Off");
}

#else

bool md_render_bridge_enable(ApplicationState* state, bool enable) { (void)state; return !enable; }
void md_render_bridge_destroy(ApplicationState*) {}
const char* md_render_bridge_status(const ApplicationState&) { return "Built without md_gpu (MD_ENABLE_GPU)"; }
bool md_render_bridge_active(const ApplicationState&) { return false; }
void md_render_bridge_marks(ApplicationState*) {}
bool md_render_bridge_draws(const ApplicationState&, const Representation&) { return false; }
void md_render_bridge_pick(ApplicationState*, float, float) {}
void md_render_bridge_resolve_hit(ApplicationState*, PickingHit*) {}
void md_render_bridge_rep_colors(ApplicationState*, const Representation&, const uint32_t*, size_t) {}
void md_render_bridge_rep_pairs(ApplicationState*, const Representation&, const uint32_t*, const uint32_t*, const uint32_t*, size_t) {}
void md_render_bridge_buffers(ApplicationState*, uint32_t) {}
void md_render_bridge_render(ApplicationState*) {}
void md_render_bridge_merge(ApplicationState*) {}

#endif
