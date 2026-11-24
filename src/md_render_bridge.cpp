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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if MD_ENABLE_GPU

#include <md_render.h>
#include <md_render_mol.h>
#include <core/md_gpu_gl.h>

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
};

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

    // The system: its topology, the per-frame block drawn, the positions, the selection and highlight
    bool                 has_topo = false;
    md_render_mol_topo_t topo = {};
    md_gpu_addr_t        dynamic = 0;
    md_gpu_addr_t        positions = 0;
    md_gpu_addr_t        selection = 0;
    md_gpu_addr_t        highlight = 0;
    uint32_t             atom_count = 0;
    uint32_t             pending = 0;          // MolBit_ flags not yet taken up
    bool                 updated = false;      // the per-frame block holds the positions
    int64_t              regroup_frame = -1;   // the trajectory frame of the last regroup
    md_gpu_addr_t        secondary = 0;        // int32 per protein backbone segment
    uint32_t             secondary_count = 0;
    uint64_t             secondary_hash = 0;
    uint64_t             content[MD_RENDER_CONTENT_COUNT] = {};

    md_array(BridgeRep)         reps = nullptr;
    md_array(md_render_batch_t) batches = nullptr;

    // The export, shared with GL
    uint32_t  width = 0, height = 0;
    SharedTex depth, color, surface;
    bool      exported = false;     // the textures hold a frame
    bool      fresh = false;        // exported this frame: GL reads after md_gpu's release (begin_read, end_read)
    bool      still = true;         // md_render did not rasterise this frame: nothing moved

    // The GL merge
    GLuint program = 0;
    GLuint vao = 0;
    GLint  loc_depth = -1, loc_color = -1, loc_surface = -1, loc_size = -1, loc_dim = -1, loc_still = -1;
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
    out_atom_index = vec4(1.0);     // INVALID_PICKING_IDX
}
)";

static bool gl_init(MdRenderBridge* b, bool rect) {
    char defines[64];
    snprintf(defines, sizeof(defines), "%s", rect ? "#define RECT 1\n" : "");
    char fs_src[4096];
    snprintf(fs_src, sizeof(fs_src), "#version 330 core\n%s%s", defines, merge_fs);
    GLuint vs = gl::compile_shader_from_source(md_str_from_cstr(merge_vs), GL_VERTEX_SHADER);
    GLuint fs = gl::compile_shader_from_source(md_str_from_cstr(fs_src), GL_FRAGMENT_SHADER);
    if (!vs || !fs) {
        if (vs) glDeleteShader(vs);
        if (fs) glDeleteShader(fs);
        return false;
    }
    b->program = glCreateProgram();
    const GLuint shaders[] = {vs, fs};
    const bool linked = gl::attach_link_detach(b->program, shaders, 2);
    glDeleteShader(vs);
    glDeleteShader(fs);
    if (!linked) {
        glDeleteProgram(b->program);
        b->program = 0;
        return false;
    }
    b->loc_depth   = glGetUniformLocation(b->program, "u_depth");
    b->loc_color   = glGetUniformLocation(b->program, "u_color");
    b->loc_surface = glGetUniformLocation(b->program, "u_surface");
    b->loc_size    = glGetUniformLocation(b->program, "u_size");
    b->loc_dim     = glGetUniformLocation(b->program, "u_dim");
    b->loc_still   = glGetUniformLocation(b->program, "u_still");
    glGenVertexArrays(1, &b->vao);
    return true;
}

// ---- Shared textures ---------------------------------------------------------------------------------------------

// A frame exported and not yet read by GL is read (without drawing): the two take turns with each texture, and
// md_gpu writes it again only after GL said it is done
static void settle_reads(MdRenderBridge* b) {
    if (!b->fresh) return;
    SharedTex* t[3] = {&b->depth, &b->color, &b->surface};
    for (SharedTex* s : t) md_gpu_gl_begin_read(&s->gl);
    for (SharedTex* s : t) md_gpu_gl_end_read(&s->gl);
    b->fresh = false;
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

static void free_system(MdRenderBridge* b) {
    if (b->has_topo) md_gpu_free(b->stream, b->topo.base);
    if (b->dynamic)   md_gpu_free(b->stream, b->dynamic);
    if (b->positions) md_gpu_free(b->stream, b->positions);
    if (b->selection) md_gpu_free(b->stream, b->selection);
    if (b->highlight) md_gpu_free(b->stream, b->highlight);
    if (b->secondary) md_gpu_free(b->stream, b->secondary);
    for (size_t i = 0; i < md_array_size(b->reps); ++i) {
        if (b->reps[i].color)   md_gpu_free(b->stream, b->reps[i].color);
        if (b->reps[i].visible) md_gpu_free(b->stream, b->reps[i].visible);
        if (b->reps[i].profile) md_gpu_free(b->stream, b->reps[i].profile);
    }
    md_array_shrink(b->reps, 0);
    b->has_topo = false;
    b->topo = {};
    b->dynamic = b->positions = b->selection = b->highlight = b->secondary = 0;
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
    b->dynamic   = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, b->topo.dynamic_size).gpu;
    b->positions = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, (size_t)n * 12).gpu;
    b->selection = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, bitset_words32(n) * 4).gpu;
    b->highlight = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, bitset_words32(n) * 4).gpu;
    if (!b->dynamic || !b->positions || !b->selection || !b->highlight) {
        set_status(b, "Out of GPU memory for the system%s");
        free_system(b);
        return false;
    }
    md_gpu_memset(b->stream, b->selection, 0, bitset_words32(n) * 4);
    md_gpu_memset(b->stream, b->highlight, 0, bitset_words32(n) * 4);
    // All of it, then
    b->pending |= MolBit_DirtyPosition | MolBit_DirtyFlags | MolBit_ResetBackboneHistory;
    // The representations' arrays are of the atoms: given again for the new system
    flag_all_representations_as_dirty(state);
    return true;
}

static md_render_mol_atoms_t atoms_of(const MdRenderBridge* b, const ApplicationState* state) {
    md_render_mol_atoms_t a = {};
    a.position  = b->positions;
    a.selection = b->selection;
    a.highlight = b->highlight;
    md_render_mol_frame_cell(a.cell, &state->mold.state.unitcell);
    return a;
}

// What changed since the scene was last made: positions, the selection and highlight
static bool take_up(MdRenderBridge* b, ApplicationState* state) {
    const uint32_t dirty = b->pending;
    b->pending = 0;
    if (dirty & MolBit_DirtyFlags) {
        upload_set(b, b->selection, b->atom_count, md_bitmask_view(&state->selection.selection_mask), state->allocator.frame);
        upload_set(b, b->highlight, b->atom_count, md_bitmask_view(&state->selection.highlight_mask), state->allocator.frame);
        b->content[MD_RENDER_CONTENT_STATE] += 1;
    }
    if ((dirty & (MolBit_DirtyPosition | MolBit_DirtyRadius | MolBit_ResetBackboneHistory)) || !b->updated) {
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
        const md_render_mol_atoms_t atoms = atoms_of(b, state);
        if (!md_render_mol_update(b->mol, b->stream, &b->topo, b->dynamic, &atoms, flags, b->content)) {
            set_status(b, "md_render_mol could not update the positions (see the log)%s");
            return false;
        }
        b->updated = true;
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
           t == RepresentationType::Ribbons || t == RepresentationType::Cartoon;
}

static bool active(const ApplicationState& state) {
    const MdRenderBridge* b = state.md_render_bridge;
    return state.md_render.enabled && b && b->ok && b->has_topo;
}

bool md_render_bridge_draws(const ApplicationState& state, const Representation& rep) {
    if (!active(state)) return false;
    if (!type_drawn(rep.type) || rep.sphere_fallback || !rep.type_is_valid || !representation_shown(state, rep)) return false;
    // md_gl's for now: the parts of an assembly shown, and bonds formed and broken over the run
    if (state.view.assembly >= 0) return false;
    if ((rep.type == RepresentationType::Licorice || rep.type == RepresentationType::BallAndStick) && bond_tracking_delta(&state)) return false;
    const BridgeRep* r = find_rep(state.md_render_bridge, rep.node.id);
    return r && r->color && r->visible && r->atom_count == state.md_render_bridge->atom_count;
}

void md_render_bridge_rep_colors(ApplicationState* state, const Representation& rep, const uint32_t* colors, size_t count) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!state->md_render.enabled || !b || !b->ok || !type_drawn(rep.type)) return;
    if (!b->has_topo || count != b->atom_count || !colors) return;
    BridgeRep* r = find_rep(b, rep.node.id);
    if (!r) {
        BridgeRep nr = {};
        nr.id = rep.node.id;
        md_array_push(b->reps, nr, state->allocator.persistent);
        r = &b->reps[md_array_size(b->reps) - 1];
    }
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

// ---- The frame ---------------------------------------------------------------------------------------------------

void md_render_bridge_render(ApplicationState* state) {
    MdRenderBridge* b = state->md_render_bridge;
    if (b) b->still = true;
    if (!state->md_render.enabled || !b || !b->ok) return;
    settle_reads(b);
    md_array_shrink(b->batches, 0);

    if (!ensure_system(b, state)) return;
    if (!take_up(b, state)) return;

    // Forget the representations that are no more
    for (size_t i = 0; i < md_array_size(b->reps); ++i) b->reps[i].seen = false;
    const md_render_mol_atoms_t atoms = atoms_of(b, state);
    const uint64_t secondary = update_secondary(b, state);
    md_allocator_i* alloc = state->allocator.persistent;     // the batches' array is kept

    md_render_mol_bond_style_t styles[64][MD_RENDER_MOL_BOND_KIND_COUNT];
    size_t num_styles = 0;
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        const Representation& rep = state->representation.reps[i];
        BridgeRep* r = find_rep(b, rep.node.id);
        if (r) r->seen = true;
        if (!md_render_bridge_draws(*state, rep) || !r) continue;

        md_render_mol_rep_t mr = {};
        mr.key   = (uint64_t)rep.node.id + 1;
        mr.color.mode = MD_RENDER_COLOR_POINTS;
        mr.color.sharpness = rep.bond_color == BondColorMode::SmoothAtom ? 1.0f - MD_CLAMP(rep.bond_sharpness, 0.0f, 1.0f) : 1.0f;
        mr.atom_color = r->color;
        mr.visible = r->visible;
        mr.chain_gaps = rep.chain_gaps;
        switch (rep.type) {
        case RepresentationType::SpaceFill:
            mr.type = MD_RENDER_MOL_REP_SPACEFILL;
            mr.params.spacefill.radius_scale = rep.scale.x;
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
        default:
            continue;
        }
        if (rep.type == RepresentationType::Licorice || rep.type == RepresentationType::BallAndStick) {
            if (num_styles == MD_COUNTOF(styles)) continue;
            bond_styles(styles[num_styles], rep);
            mr.bond_style = styles[num_styles++];
            mr.bond_kinds = bond_kinds(rep);
        }
        const size_t copies = md_array_size(state->view.periodic_transforms);
        if (rep.periodic_images && copies > 1) {
            mr.transform = &state->view.periodic_transforms[0].elem[0][0];
            mr.transform_count = (uint32_t)copies;
        }
        emit(b, atoms, mr, alloc);
    }
    for (size_t i = 0; i < md_array_size(b->reps);) {
        if (!b->reps[i].seen) {
            BridgeRep& r = b->reps[i];
            if (r.color)   md_gpu_free(b->stream, r.color);
            if (r.visible) md_gpu_free(b->stream, r.visible);
            if (r.profile) md_gpu_free(b->stream, r.profile);
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

    // Every frame the scene is made, even when md_render did nothing: GL reads the textures only between its
    // begin_read and end_read, which order its reading before md_gpu's next writing
    SharedTex* t[3] = {&b->depth, &b->color, &b->surface};
    for (SharedTex* s : t) md_gpu_texture_share_acquire(b->stream, s->tex);
    md_render_export_t ex = {b->depth.tex, b->color.tex, b->surface.tex};
    const bool ok = md_render_view_export(b->stream, b->view, &ex);
    for (SharedTex* s : t) md_gpu_texture_share_release(b->stream, s->tex);
    b->fresh = true;
    b->exported = ok;
    snprintf(b->status, sizeof(b->status), "Drawing %u batches, %u x %u", f.batch_count, width, height);
}

void md_render_bridge_merge(ApplicationState* state) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!state->md_render.enabled || !b || !b->ok || !b->program) return;
    if (!b->exported) {
        settle_reads(b);
        return;
    }
    SharedTex* t[3] = {&b->depth, &b->color, &b->surface};
    if (b->fresh) {
        for (SharedTex* s : t) md_gpu_gl_begin_read(&s->gl);
    }

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
    glBindVertexArray(b->vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    for (int i = 2; i >= 0; --i) {
        glActiveTexture(GL_TEXTURE0 + units[i]);
        glBindTexture(t[i]->gl.target, 0);
    }
    glUseProgram(0);

    if (b->fresh) {
        for (SharedTex* s : t) md_gpu_gl_end_read(&s->gl);
        b->fresh = false;
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
bool md_render_bridge_draws(const ApplicationState&, const Representation&) { return false; }
void md_render_bridge_rep_colors(ApplicationState*, const Representation&, const uint32_t*, size_t) {}
void md_render_bridge_buffers(ApplicationState*, uint32_t) {}
void md_render_bridge_render(ApplicationState*) {}
void md_render_bridge_merge(ApplicationState*) {}

#endif
