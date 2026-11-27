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

#include <imgui.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if MD_ENABLE_GPU

#include <md_render.h>
#include <md_render_mol.h>
#include <md_compose.h>
#include <md_post.h>
#include <gfx/postprocessing_utils.h>
#include <gfx/isolines.h>
#include <gpu_volume.h>
#include <iso_meshes.h>
#include <es_planes.h>
#include <color_scale.h>
#include <transfer_function.h>
#include <transfer_function_widget.h>
#include <implot.h>

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
    // kept (the frame, or the history written last); the image shown
    md_compose_ctx_t compose = nullptr;
    md_post_ctx_t    post = nullptr;
    uint32_t  width = 0, height = 0;            // the scene's (the G-buffer's)
    uint32_t  out_width = 0, out_height = 0;    // the image's
    md_gpu_texture_t scene = nullptr, frame = nullptr;
    md_gpu_texture_t history[2] = {};
    md_gpu_texture_t shown = nullptr;           // the TAA's image with the motion blur, while it is on (else the history)
    md_compose_history_t vis_history = nullptr; // the shadows' and the occlusion's, while they are filtered over time
    int       hist = 0;                         // the history written last
    bool      hist_valid = false;               // it holds the image of the frame before, for the TAA to go on from
    md_gpu_texture_t image = nullptr;
    md_gpu_texture_t out = nullptr;             // the image shown: the image kept presented, the overlay over it
    bool      still = true;         // md_render did not rasterise this frame: nothing moved
    md_render_frame_t last = {};    // the frame drawn last (its batches are `batches`), for the masks alone
    bool      has_last = false;
    int       last_ui_frame = -1;   // the ImGui frame it was drawn in
    uint32_t  presents = 0;         // seeds the dither

    // The overlay drawn last (draw_overlay), for picking it: its items, the camera and the size it was drawn at; and
    // the pixel asked about last (md_render_bridge_pick: the scene's, rows from the top)
    md_array(OverlayItem) overlay = nullptr;
    md_render_camera_t    overlay_camera = {};
    float                 pick_x = -1.0f, pick_y = -1.0f;

    // The lookup tables of the electronic structure representations (lut_texture): their transfer functions and colour
    // maps, each kept while it is drawn
    struct Lut { uint64_t key; md_gpu_texture_t tex; bool used; };
    md_array(Lut) luts = nullptr;
    bool          es_dropped_logged = false;    // that some were left out (collect_electronic_structure), once

    // The colours of the vertices of the isosurfaces drawn as meshes (emit_meshes), a representation's each, by what
    // they are of, kept while drawn; what the meshes were as drawn last (a new one may be where an old one was)
    struct MeshColors { uint32_t rep_id; uint64_t key; md_gpu_addr_t color; uint32_t cap; bool used; };
    md_array(MeshColors) mesh_colors = nullptr;
    uint64_t             meshes_key = 0;
    md_float4x4_t        mesh_model = {};     // their bohr to the world

    // The glyph representations' (emit_glyphs), by representation
    struct Glyphs {
        uint32_t      rep_id = 0;
        uint32_t      cap = 0;              // glyphs the buffers hold
        md_gpu_addr_t points = 0, radius = 0, colors = 0, attr = 0, marks = 0;
        uint32_t      count = 0, version = 0, color_version = 0;
        bool          cylinders = false, valid = false, used = false;
        uint32_t      marked = UINT32_MAX;  // the glyph marked as the highlight
    };
    md_array(Glyphs) glyphs = nullptr;


    // md_render's answers to the picks, the newest last
    static constexpr int PICKS = 8;
    md_render_pick_t picks[PICKS] = {};
    uint32_t         num_picks = 0;
};

static void set_status(MdRenderBridge* b, const char* fmt, const char* detail = nullptr) {
    snprintf(b->status, sizeof(b->status), fmt, detail ? detail : "");
}

static md_gpu_texture_t gpu_texture(MdRenderBridge* b, md_gpu_format_t format, md_gpu_tex_usage_t usage, uint32_t w, uint32_t h, const char* label) {
    md_gpu_texture_desc_t d = {};
    d.type = MD_GPU_TEX_2D; d.format = format; d.usage = usage; d.width = w; d.height = h; d.label = label;
    return md_gpu_texture_create(b->stream, &d);
}

static void free_targets(MdRenderBridge* b) {
    md_gpu_texture_t* t[] = {&b->scene, &b->frame, &b->history[0], &b->history[1], &b->shown};
    for (md_gpu_texture_t* x : t) {
        if (*x) md_gpu_texture_destroy(*x);
        *x = nullptr;
    }
    if (b->out) md_gpu_texture_destroy(b->out);
    b->out = nullptr;
    b->image = nullptr;
    b->hist_valid = false;
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
    b->shown      = gpu_texture(b, MD_GPU_FORMAT_RGBA16_FLOAT, SS, out_width, out_height, "main view shown");
    if (!b->scene || !b->frame || !b->history[0] || !b->history[1] || !b->shown) {
        set_status(b, "Could not make the main view's targets: %s", md_gpu_last_error());
        free_targets(b);
        return false;
    }
    // Written by md_post_present, drawn over by md_post_overlay, drawn by ImGui
    b->out = gpu_texture(b, MD_GPU_FORMAT_RGBA8_UNORM, MD_GPU_TEX_STORAGE | MD_GPU_TEX_RENDER_TARGET | MD_GPU_TEX_SAMPLED, out_width, out_height,
                         "main view image");
    if (!b->out) {
        set_status(b, "Could not make the main view's image: %s", md_gpu_last_error());
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
static md_gpu_addr_t cartoon_profile(MdRenderBridge* b, md_gpu_addr_t* profile, uint64_t* profile_key, const Representation& rep, uint64_t secondary) {
    const float scale[3] = {rep.scale.x, rep.scale.z, rep.scale.y};     // coil, helix, sheet
    uint64_t key = md_hash64(scale, sizeof(scale), secondary) | 1;
    if (!*profile) {
        *profile = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, md_render_mol_cartoon_size(&b->topo)).gpu;
        *profile_key = 0;
        if (!*profile) return 0;
    }
    if (key != *profile_key) {
        if (!md_render_mol_cartoon_profile(b->mol, b->stream, &b->topo, *profile, secondary ? b->secondary : 0, scale, b->content)) return 0;
        *profile_key = key;
    }
    return *profile;
}

// A representation of atoms as md_render_mol's: its type, its parameters and its bonds' styles (into styles, which emit
// reads), its cartoon's profile (in *profile, of *profile_key). False of a type of something else, and of Bonds (the
// pairs of the main view's: md_render_bridge_rep_pairs). Its key, colours and visibility are the caller's.
static bool mol_rep_params(md_render_mol_rep_t* mr, md_render_mol_bond_style_t styles[MD_RENDER_MOL_BOND_KIND_COUNT], MdRenderBridge* b,
                           const Representation& rep, md_gpu_addr_t* profile, uint64_t* profile_key, uint64_t secondary) {
    mr->color.sharpness = rep.bond_color == BondColorMode::SmoothAtom ? 1.0f - MD_CLAMP(rep.bond_sharpness, 0.0f, 1.0f) : 1.0f;
    mr->chain_gaps = rep.chain_gaps;
    // While the system's bonds are being inferred, every type as spheres (representation_draw_op)
    const RepresentationType type = rep.sphere_fallback ? RepresentationType::SpaceFill : rep.type;
    switch (type) {
    case RepresentationType::SpaceFill:
        mr->type = MD_RENDER_MOL_REP_SPACEFILL;
        mr->params.spacefill.radius_scale = rep.sphere_fallback ? 1.0f : rep.scale.x;
        return true;
    case RepresentationType::Licorice:
        mr->type = MD_RENDER_MOL_REP_LICORICE;
        mr->params.licorice.radius = 0.2f * rep.scale.x;                  // md_gl's scales
        mr->bond_orders = rep.bond_orders;
        break;
    case RepresentationType::BallAndStick:
        mr->type = MD_RENDER_MOL_REP_BALL_AND_STICK;
        mr->params.ball_and_stick.ball_scale = 0.2f * rep.scale.x;
        mr->params.ball_and_stick.stick_radius = 0.2f * rep.scale.y;
        mr->bond_orders = rep.bond_orders;
        break;
    case RepresentationType::Ribbons:
        mr->type = MD_RENDER_MOL_REP_RIBBONS;
        mr->params.ribbons.half_width = rep.scale.x;
        mr->params.ribbons.half_thickness = 0.1f * rep.scale.y;
        return true;
    case RepresentationType::Cartoon: {
        mr->type = MD_RENDER_MOL_REP_CARTOON;
        const float w = rep.scale.w > 0.0f ? rep.scale.w : 1.0f;
        mr->params.cartoon.profile = cartoon_profile(b, profile, profile_key, rep, secondary);
        mr->params.cartoon.band_radius = MD_MAX(rep.nucleic_band_radius, 0.01f) * w;
        mr->params.cartoon.ring_radius = MD_MAX(rep.nucleic_ring_radius, 0.01f) * w;
        mr->params.cartoon.ring_plate  = 0.6f;
        mr->params.cartoon.lipid_scale = w;
        mr->params.cartoon.split_heads = rep.two_head_spheres;
        return true;
    }
    default:
        return false;
    }
    // the sticks
    bond_styles(styles, rep);
    mr->bond_kinds = bond_kinds(rep);
    mr->bond_style = styles;
    return true;
}

// A representation's batches, after those in *out
static bool emit(MdRenderBridge* b, md_render_batch_t** out, const md_render_mol_atoms_t& atoms, const md_render_mol_rep_t& mr, md_allocator_i* alloc) {
    md_render_batch_t*& o = *out;
    const size_t at = md_array_size(o);
    uint32_t cap = (uint32_t)(md_array_capacity(o) - at);
    md_array_resize(o, at + cap, alloc);
    uint32_t n = md_render_mol_emit(&b->topo, b->dynamic, &atoms, &mr, o + at, cap);
    if (n > cap) {
        md_array_resize(o, at + n, alloc);
        n = md_render_mol_emit(&b->topo, b->dynamic, &atoms, &mr, o + at, n);
    }
    md_array_resize(o, at + n, alloc);
    return n > 0;
}

// md_gl's view as md_render's camera: the same view, the projection's depth reversed (d = (1 - ndc z) / 2, so that
// md_gl's depth is 1 - d) and its jitter given apart, in pixels (x right, y down)
static void camera_from(md_render_camera_t* cam, const md_float4x4_t& V, const md_float4x4_t& P) {
    md_memcpy(cam->view, &V.elem[0][0], sizeof(cam->view));
    for (int c = 0; c < 4; ++c) {
        for (int r = 0; r < 4; ++r) cam->proj[c * 4 + r] = P.elem[c][r];
        cam->proj[c * 4 + 2] = 0.5f * (P.elem[c][3] - P.elem[c][2]);
    }
}

static void camera_of(md_render_camera_t* cam, float jitter[2], const ApplicationState& state, uint32_t width, uint32_t height) {
    const ViewParam& p = state.view.param;
    const md_float4x4_t& P = p.matrix.curr.proj_no_jitter;
    const md_float4x4_t& J = p.matrix.curr.proj;
    camera_from(cam, p.matrix.curr.view, P);
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
            emit(b, &b->batches, pa, pr, alloc);
        }
        return;
    }
    const size_t copies = md_array_size(state->view.periodic_transforms);
    if (rep.periodic_images && copies > 1) {
        mr.transform = &state->view.periodic_transforms[0].elem[0][0];
        mr.transform_count = (uint32_t)copies;
    }
    emit(b, &b->batches, atoms, mr, alloc);
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
enum { SHAPE_SHAFTS, SHAPE_HEADS, SHAPE_CELL, SHAPE_IMAGES, SHAPE_AXES, SHAPE_ORIGIN, SHAPE_GLYPHS };
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
        const float h_render = (float)MD_MAX(state->target.height, 1U);
        const float h_px     = (float)MD_MAX(state->target.out_height, 1U);
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

// ---- Glyphs -------------------------------------------------------------------------------------------------------

// The glyph under the cursor, in the view or in a table (glyph_hover_ui first, as it is marked in place of the other)
static void glyph_hovered(const ApplicationState& state, uint32_t* rep_id, uint32_t* point) {
    const bool ui = state.glyph_hover_ui.rep_id != 0;
    *rep_id = ui ? state.glyph_hover_ui.rep_id : state.glyph_hover.rep_id;
    *point  = ui ? state.glyph_hover_ui.point : state.glyph_hover.point;
}

static void free_glyphs(MdRenderBridge* b, MdRenderBridge::Glyphs* x) {
    const md_gpu_addr_t a[5] = {x->points, x->radius, x->colors, x->attr, x->marks};
    for (md_gpu_addr_t p : a) if (p) md_gpu_free(b->stream, p);
    *x = {};
}

// The marks of the glyphs: the hovered one's highlight (state 1), the others none. Uploaded when it moved, the content's
// state counted.
static void glyph_marks(MdRenderBridge* b, const ApplicationState& state) {
    uint32_t rep_id, point;
    glyph_hovered(state, &rep_id, &point);
    for (size_t i = 0; i < md_array_size(b->glyphs); ++i) {
        MdRenderBridge::Glyphs& x = b->glyphs[i];
        if (!x.valid) continue;
        const uint32_t want = rep_id && x.rep_id == rep_id && point < x.count ? point : UINT32_MAX;
        if (want == x.marked) continue;
        md_temp_scope_t temp = md_temp_begin_in(state.allocator.frame);
        const size_t n = MD_ALIGN_TO(x.count, 4);
        uint8_t* m = (uint8_t*)md_vm_arena_push(state.allocator.frame, n);
        md_memset(m, 0, n);
        if (want != UINT32_MAX) m[want] = 2;
        md_gpu_upload(b->stream, x.marks, m, n);
        md_temp_end(temp);
        x.marked = want;
        b->content[MD_RENDER_CONTENT_STATE] += 1;
    }
}

// The batches of the glyph representations shown (glyphs.h): spheres, or flat-ended cylinders with their glyph as the
// attribute of both their points, over buffers of each representation's own, uploaded as their versions change.
// Picked as their points (md_render_bridge_hit), the hovered one marked as the highlight (glyph_marks).
static void emit_glyphs(MdRenderBridge* b, ApplicationState* state, md_allocator_i* alloc) {
    for (size_t i = 0; i < md_array_size(b->glyphs); ++i) b->glyphs[i].used = false;
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        const Representation& rep = state->representation.reps[i];
        GlyphDraw d;
        if (rep.type != RepresentationType::Glyphs || !representation_shown(*state, rep) || !rep.type_is_valid || !glyphs_draw_data(&d, rep.glyphs)) continue;
        MdRenderBridge::Glyphs* x = nullptr;
        for (size_t k = 0; k < md_array_size(b->glyphs) && !x; ++k) if (b->glyphs[k].rep_id == rep.node.id) x = &b->glyphs[k];
        if (!x) {
            md_array_push(b->glyphs, MdRenderBridge::Glyphs{}, state->allocator.persistent);
            x = &b->glyphs[md_array_size(b->glyphs) - 1];
            x->rep_id = rep.node.id;
        }
        x->used = true;
        const uint32_t n = d.count, per = d.cylinders ? 2 : 1;
        if (n > x->cap) {
            const uint32_t rep_id = x->rep_id;
            free_glyphs(b, x);
            x->rep_id = rep_id;
            x->used   = true;
            x->points = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, (size_t)n * 2 * 4 * sizeof(float)).gpu;
            x->radius = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, (size_t)n * sizeof(float)).gpu;
            x->colors = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, (size_t)n * sizeof(uint32_t)).gpu;
            x->attr   = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, (size_t)n * 2 * sizeof(uint32_t)).gpu;
            x->marks  = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, MD_ALIGN_TO(n, 4)).gpu;
            if (!x->points || !x->radius || !x->colors || !x->attr || !x->marks) {
                free_glyphs(b, x);
                x->rep_id = rep.node.id;
                x->used   = true;
                continue;
            }
            x->cap = n;
            // A cylinder's two points are of its glyph
            md_temp_scope_t temp = md_temp_begin_in(state->allocator.frame);
            uint32_t* attr = (uint32_t*)md_vm_arena_push(state->allocator.frame, (size_t)n * 2 * sizeof(uint32_t));
            for (uint32_t j = 0; j < 2 * n; ++j) attr[j] = j / 2;
            md_gpu_upload(b->stream, x->attr, attr, (size_t)n * 2 * sizeof(uint32_t));
            md_temp_end(temp);
            md_gpu_memset(b->stream, x->marks, 0, MD_ALIGN_TO(n, 4));
        }
        if (!x->valid || d.version != x->version || d.cylinders != x->cylinders || n != x->count) {
            md_gpu_upload(b->stream, x->points, d.points, (size_t)n * per * 4 * sizeof(float));
            md_gpu_upload(b->stream, x->radius, d.radius, (size_t)n * sizeof(float));
            b->content[MD_RENDER_CONTENT_GEOMETRY] += 1;
        }
        if (!x->valid || d.color_version != x->color_version || n != x->count) {
            md_gpu_upload(b->stream, x->colors, d.colors, (size_t)n * sizeof(uint32_t));
            b->content[MD_RENDER_CONTENT_APPEARANCE] += 1;
        }
        if (n != x->count) {
            md_gpu_memset(b->stream, x->marks, 0, MD_ALIGN_TO(n, 4));
            x->marked = UINT32_MAX;
        }
        x->version = d.version;
        x->color_version = d.color_version;
        x->cylinders = d.cylinders;
        x->count = n;
        x->valid = true;

        md_render_batch_t bt = {};
        bt.key   = shape_key(rep.node.id, SHAPE_GLYPHS);
        bt.count = n;
        bt.points.count    = n * per;
        bt.points.position = x->points;
        bt.points.stride   = 16;
        bt.points.color    = x->colors;
        bt.points.marks    = x->marks;
        bt.color.mode      = MD_RENDER_COLOR_POINTS;
        bt.color.sharpness = 1.0f;
        bt.color.unlit     = d.unlit;
        if (d.cylinders) {
            bt.kind = MD_RENDER_CAPSULES;
            bt.points.attr_index = x->attr;
            bt.points.attr_count = n;
            bt.style.capsules.radius = 1.0f;
            bt.style.capsules.radius_scale = x->radius;
            bt.style.capsules.flat = true;
        } else {
            bt.kind = MD_RENDER_SPHERES;
            bt.points.radius = x->points + 12;
            bt.points.radius_stride = 16;
            bt.style.spheres.radius_scale = 1.0f;
        }
        md_array_push(b->batches, bt, alloc);
    }
    // Of representations no more
    for (size_t i = 0; i < md_array_size(b->glyphs);) {
        if (b->glyphs[i].used) { ++i; continue; }
        free_glyphs(b, &b->glyphs[i]);
        md_array_swap_back_and_pop(b->glyphs, i);
    }
    glyph_marks(b, *state);
}

// ---- The isosurfaces drawn as meshes ------------------------------------------------------------------------------

// Their batches' keys: shapes', a surface of a representation's meshes each (not picked: the cells' edges' are not)
enum { SHAPE_MESH = 0x10 };

// The entries of a lookup table (a transfer function, a colour map: lut_texture)
static constexpr int LUT_RES = 256;

static uint32_t rgba8_of(md_float4_t c) {
    uint32_t v = 0;
    for (int k = 0; k < 4; ++k) v |= (uint32_t)(MD_CLAMP(c.elem[k], 0.0f, 1.0f) * 255.0f + 0.5f) << (8 * k);
    return v;
}

static float srgb_to_linear1(float c) {
    return c <= 0.04045f ? c / 12.92f : powf((c + 0.055f) / 1.055f, 2.4f);
}

// The colour volume's (RGBA8, its dims, texture space t) at t, trilinear between its texels' centres
static md_float4_t color_volume_at(const uint32_t* v, const int dim[3], md_float3_t t) {
    float f[3];
    int   i0[3], i1[3];
    for (int k = 0; k < 3; ++k) {
        const float x = MD_CLAMP(t.elem[k] * (float)dim[k] - 0.5f, 0.0f, (float)(dim[k] - 1));
        i0[k] = (int)x;
        i1[k] = MD_MIN(i0[k] + 1, dim[k] - 1);
        f[k] = x - (float)i0[k];
    }
    md_float4_t c = {0, 0, 0, 0};
    for (int n = 0; n < 8; ++n) {
        const int x = (n & 1) ? i1[0] : i0[0], y = (n & 2) ? i1[1] : i0[1], z = (n & 4) ? i1[2] : i0[2];
        const float w = ((n & 1) ? f[0] : 1.0f - f[0]) * ((n & 2) ? f[1] : 1.0f - f[1]) * ((n & 4) ? f[2] : 1.0f - f[2]);
        const uint32_t p = v[((size_t)z * dim[1] + y) * dim[0] + x];
        for (int k = 0; k < 4; ++k) c.elem[k] += w * (float)((p >> (8 * k)) & 0xFF) / 255.0f;
    }
    return c;
}

// A representation's vertex colours, of `key` (made from what it is of, on the host, when they are not), on the GPU
static md_gpu_addr_t mesh_colors(MdRenderBridge* b, uint32_t rep_id, uint64_t key, uint32_t count, md_allocator_i* alloc,
                                 void (*fill)(uint32_t* out, uint32_t count, const void* user), const void* user) {
    MdRenderBridge::MeshColors* mc = nullptr;
    for (size_t i = 0; i < md_array_size(b->mesh_colors); ++i) {
        if (b->mesh_colors[i].rep_id == rep_id) mc = &b->mesh_colors[i];
    }
    if (!mc) {
        md_array_push(b->mesh_colors, (MdRenderBridge::MeshColors{rep_id, 0, 0, 0, false}), alloc);
        mc = &b->mesh_colors[md_array_size(b->mesh_colors) - 1];
    }
    mc->used = true;
    if (mc->key == key && mc->color) return mc->color;
    if (mc->cap < count) {
        if (mc->color) md_gpu_free(b->stream, mc->color);
        mc->color = md_gpu_malloc(b->stream, MD_GPU_MEM_DEVICE, sizeof(uint32_t) * count).gpu;
        mc->cap = mc->color ? count : 0;
    }
    uint32_t* px = (uint32_t*)malloc(sizeof(uint32_t) * MD_MAX(count, 1u));
    if (!mc->color || !px) {
        free(px);
        mc->key = 0;
        return 0;
    }
    fill(px, count, user);
    const bool ok = md_gpu_upload(b->stream, mc->color, px, sizeof(uint32_t) * count);
    free(px);
    mc->key = ok ? key : 0;
    return ok ? mc->color : 0;
}

struct FieldColors { const float* values; float beg, inv_ext; uint32_t lut[LUT_RES]; };
static void fill_field_colors(uint32_t* out, uint32_t count, const void* user) {
    const FieldColors& f = *(const FieldColors*)user;
    for (uint32_t i = 0; i < count; ++i) {
        const float t = MD_CLAMP((f.values[i] - f.beg) * f.inv_ext, 0.0f, 1.0f);
        out[i] = f.lut[(int)(t * (float)(LUT_RES - 1) + 0.5f)];
    }
}

struct AtomColors { const float* xyz; const uint32_t* vol; int dim[3]; md_float4x4_t world_to_texture; };
static void fill_atom_colors(uint32_t* out, uint32_t count, const void* user) {
    const AtomColors& a = *(const AtomColors*)user;
    const float scl = (float)BOHR_TO_ANGSTROM;
    for (uint32_t i = 0; i < count; ++i) {
        const md_float4_t w = md_float4_set(a.xyz[3 * i] * scl, a.xyz[3 * i + 1] * scl, a.xyz[3 * i + 2] * scl, 1.0f);
        const md_float3_t t = md_float3_from_float4(md_float4x4_mul_float4(a.world_to_texture, w));
        md_float4_t c = color_volume_at(a.vol, a.dim, t);
        c.w = 1.0f;
        out[i] = rgba8_of(c);
    }
}

// The isosurfaces of the electronic structure representations drawn as meshes (es_mesh_mode): transparent triangles
// of md_render's, an isovalue's sheet a batch, of its colour and coverage; by a field, its colour map's colour at each
// vertex (the coverage the surface's); by atom, the colour volume's there times the surface's colour. Without their
// optical density inside, and the field's isolines.
static void emit_meshes(MdRenderBridge* b, ApplicationState* state, md_allocator_i* alloc) {
    // read when the frame is drawn
    b->mesh_model = iso_meshes_model_matrix();
    uint64_t drawn = 0;
    for (size_t i = 0; i < md_array_size(b->mesh_colors); ++i) b->mesh_colors[i].used = false;
    for (size_t r = 0; r < md_array_size(state->representation.reps); ++r) {
        const Representation& rep = state->representation.reps[r];
        if (rep.type != RepresentationType::ElectronicStructure || !representation_shown(*state, rep) || !rep.type_is_valid) continue;
        const ElectronicStructureRepresentation& es = rep.electronic_structure;
        if (es.display != ElectronicStructureDisplay::Isosurface || !es_mesh_mode(state, es)) continue;
        EsMeshDraw d;
        if (!electronic_structure_mesh_draw(&d, es)) continue;
        md_gpu_stream_wait(b->stream, d.ready);

        // The vertices' colours: by the field, by atom, or none (the surfaces' own)
        md_gpu_addr_t colors = 0;
        bool tinted = false;
        if (d.values) {
            FieldColors fc;
            fc.values = d.values;
            fc.beg = es.field_map.range_beg;
            const float ext = es.field_map.range_end - es.field_map.range_beg;
            fc.inv_ext = ext != 0.0f ? 1.0f / ext : 0.0f;
            for (int i = 0; i < LUT_RES; ++i) {
                ImVec4 c = ImPlot::SampleColormap((float)i / (float)(LUT_RES - 1), es.field_map.colormap);
                c.w = 1.0f;
                fc.lut[i] = ImGui::ColorConvertFloat4ToU32(c);
            }
            uint64_t key = md_hash64(&d.values_version, sizeof(d.values_version), 0xf1e1d);
            const float m[3] = {(float)es.field_map.colormap, es.field_map.range_beg, es.field_map.range_end};
            key = md_hash64(m, sizeof(m), key);
            colors = mesh_colors(b, rep.node.id, key, d.num_vertices, alloc, fill_field_colors, &fc);
        } else if (es.coloring == SurfaceColoring::AtomColors && es.color_rgba && d.xyz) {
            AtomColors ac;
            ac.xyz = d.xyz;
            ac.vol = es.color_rgba;
            for (int k = 0; k < 3; ++k) ac.dim[k] = es.color_vol.dim[k];
            ac.world_to_texture = md_float4x4_inverse(es.color_vol.texture_to_world);
            uint64_t key = md_hash64(&es.col_hash, sizeof(es.col_hash), d.version);
            key = md_hash64(&ac.vol, sizeof(ac.vol), key);
            colors = mesh_colors(b, rep.node.id, key, d.num_vertices, alloc, fill_atom_colors, &ac);
            tinted = colors != 0;
        }
        drawn = md_hash64(&d.version, sizeof(d.version), drawn);
        drawn = md_hash64(&colors, sizeof(colors), drawn);

        for (uint32_t i = 0; i < d.num_parts; ++i) {
            const EsMeshDraw::Part& p = d.parts[i];
            const float alpha = MD_CLAMP(p.color.w, 0.0f, 1.0f);
            if (alpha <= 0.0f || p.index_count < 3) continue;
            md_render_batch_t bt = {};
            bt.key = shape_key(rep.node.id, SHAPE_MESH + i);
            bt.kind = MD_RENDER_TRIANGLES;
            bt.count = p.index_count / 3;
            bt.points.count = d.num_vertices;
            bt.points.position = d.position;
            bt.points.stride = 16;
            bt.index = d.index + (md_gpu_addr_t)p.index_offset * sizeof(uint32_t);
            bt.style.triangles.normal = d.normal;
            bt.style.triangles.normal_stride = 16;
            if (colors) {
                bt.points.color = colors;
                bt.color.mode = MD_RENDER_COLOR_POINTS;
                if (tinted) {
                    for (int k = 0; k < 3; ++k) bt.color.tint[k] = srgb_to_linear1(MD_CLAMP(p.color.elem[k], 0.0f, 1.0f));
                    bt.color.tint[3] = 1.0f;
                }
            } else {
                bt.color.mode = MD_RENDER_COLOR_UNIFORM;
                bt.color.uniform = rgba8_of(md_float4_set(p.color.x, p.color.y, p.color.z, 1.0f));
            }
            bt.color.transparency = alpha < 1.0f ? 1.0f - alpha : 0.0f;
            bt.transform = &b->mesh_model.elem[0][0];
            bt.transform_count = 1;
            md_array_push(b->batches, bt, alloc);
        }
    }
    // the colours of what is no more drawn
    for (size_t i = 0; i < md_array_size(b->mesh_colors);) {
        if (!b->mesh_colors[i].used) {
            if (b->mesh_colors[i].color) md_gpu_free(b->stream, b->mesh_colors[i].color);
            md_array_swap_back_and_pop(b->mesh_colors, i);
        } else {
            ++i;
        }
    }
    // A mesh or its colours made where another was: its contents changed under the same addresses
    if (drawn != b->meshes_key) {
        b->meshes_key = drawn;
        b->content[MD_RENDER_CONTENT_GEOMETRY] += 1;
        b->content[MD_RENDER_CONTENT_APPEARANCE] += 1;
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

// ---- The electronic structure representations ---------------------------------------------------------------------

// A lookup table of LUT_RES x 1 RGBA8 texels by what it is of (key), made from px (filled by the caller on a miss: px
// null asks whether it is there), and kept while it is drawn (luts_sweep)
static md_gpu_texture_t lut_texture(MdRenderBridge* b, uint64_t key, const uint32_t* px, md_allocator_i* alloc) {
    for (size_t i = 0; i < md_array_size(b->luts); ++i) {
        if (b->luts[i].key == key) {
            b->luts[i].used = true;
            return b->luts[i].tex;
        }
    }
    if (!px) return nullptr;
    md_gpu_texture_desc_t d = {};
    d.type = MD_GPU_TEX_2D;
    d.format = MD_GPU_FORMAT_RGBA8_UNORM;
    d.usage = MD_GPU_TEX_SAMPLED;
    d.width = LUT_RES;
    d.height = 1;
    d.depth_or_layers = 1;
    d.label = "Lookup table";
    md_gpu_texture_t tex = md_gpu_texture_create(b->stream, &d);
    if (!tex || !md_gpu_upload_texture(b->stream, tex, nullptr, px, sizeof(uint32_t) * LUT_RES)) {
        if (tex) md_gpu_texture_destroy(tex);
        return nullptr;
    }
    md_array_push(b->luts, (MdRenderBridge::Lut{key, tex, true}), alloc);
    return tex;
}

// The tables not drawn since the last sweep released
static void luts_sweep(MdRenderBridge* b) {
    for (size_t i = 0; i < md_array_size(b->luts);) {
        if (!b->luts[i].used) {
            md_gpu_texture_destroy(b->luts[i].tex);
            md_array_swap_back_and_pop(b->luts, i);
        } else {
            b->luts[i].used = false;
            ++i;
        }
    }
}

// A transfer function as md_compose reads it: the colour map's sRGB colour, the opacity of a slab 1/150 of the box thick
static md_gpu_texture_t tf_texture(MdRenderBridge* b, const TransferFunction& tf, md_allocator_i* alloc) {
    const uint64_t key = md_hash64_combine(transfer_function_hash(tf), 0x7f0a11ce);
    if (md_gpu_texture_t t = lut_texture(b, key, nullptr, alloc)) return t;
    float f[LUT_RES * 4];
    transfer_function_bake_rgba(tf, f, LUT_RES);
    uint32_t px[LUT_RES];
    for (int i = 0; i < LUT_RES; ++i) {
        uint32_t c = 0;
        for (int k = 0; k < 4; ++k) c |= (uint32_t)(MD_CLAMP(f[4 * i + k], 0.0f, 1.0f) * 255.0f + 0.5f) << (8 * k);
        px[i] = c;
    }
    return lut_texture(b, key, px, alloc);
}

// A colour map (ImPlot's), opaque
static md_gpu_texture_t colormap_texture(MdRenderBridge* b, int colormap, md_allocator_i* alloc) {
    const uint64_t key = md_hash64(&colormap, sizeof(colormap), 0xc0105ca1e);
    if (md_gpu_texture_t t = lut_texture(b, key, nullptr, alloc)) return t;
    uint32_t px[LUT_RES];
    for (int i = 0; i < LUT_RES; ++i) {
        ImVec4 c = ImPlot::SampleColormap((float)i / (float)(LUT_RES - 1), colormap);
        c.w = 1.0f;
        px[i] = ImGui::ColorConvertFloat4ToU32(c);
    }
    return lut_texture(b, key, px, alloc);
}

static void mat_store(float out[16], const md_float4x4_t& m) {
    md_memcpy(out, &m.elem[0][0], sizeof(float) * 16);
}

// What the frame reads of a volume: written by the evaluations' stream, waited for on the bridge's
static md_gpu_texture_t volume_take(MdRenderBridge* b, const GpuVolume& v) {
    if (v.tex) md_gpu_stream_wait(b->stream, v.ready);
    return v.tex;
}

// The electronic structure representations shown, as md_compose draws them: the raycast isosurfaces and the volume
// renderings as its volumes (at most MD_COMPOSE_MAX_VOLUMES), the planes as its slices (MD_COMPOSE_MAX_SLICES), what is
// past those left out. Of a volume not there yet, nothing; of a slice, its edge. The isosurfaces drawn as meshes are
// not drawn yet.
static void collect_electronic_structure(MdRenderBridge* b, ApplicationState* state, md_compose_volume_t* vols, uint32_t* num_vols,
                                         md_compose_slice_t* slices, uint32_t* num_slices) {
    md_allocator_i* alloc = state->allocator.persistent;
    *num_vols = *num_slices = 0;
    // a point of the GUI in the target's pixels (an isoline's width); a slice's edge 1.5 points wide, as the outlines
    const float h_pt = ImGui::GetIO().DisplaySize.y;
    const float px_per_point = h_pt > 0.0f ? (float)b->height / h_pt : 1.0f;
    const float edge_px = 0.75f * ImGui::GetIO().DisplayFramebufferScale.x * (float)b->width / (float)MD_MAX(state->app.framebuffer.width, 1);
    bool dropped = false;

    for (size_t r = 0; r < md_array_size(state->representation.reps); ++r) {
        Representation& rep = state->representation.reps[r];
        if (rep.type != RepresentationType::ElectronicStructure || !representation_shown(*state, rep) || !rep.type_is_valid) continue;
        ElectronicStructureRepresentation& es = rep.electronic_structure;
        const GpuVolume dv = gpu_volume_get(es.density_vol.gpu);

        if (es.display == ElectronicStructureDisplay::Isosurface || es.display == ElectronicStructureDisplay::Volume) {
            const bool iso = es.display == ElectronicStructureDisplay::Isosurface;
            if (!dv.tex || (iso && es_mesh_mode(state, es))) continue;
            if (*num_vols == MD_COMPOSE_MAX_VOLUMES) {
                dropped = true;
                continue;
            }
            md_compose_volume_t& v = vols[(*num_vols)++];
            v = {};
            v.field = volume_take(b, dv);
            mat_store(v.model, dv.texture_to_world);
            for (int k = 0; k < 3; ++k) {
                v.clip_min[k] = 0.0f;
                v.clip_max[k] = 1.0f;
            }
            if (iso) {
                IsoDesc d;
                electronic_structure_iso_desc_init(&d, es);
                v.iso_count = (uint32_t)MD_MIN(d.count, (size_t)MD_COMPOSE_MAX_ISO);
                for (uint32_t i = 0; i < v.iso_count; ++i) {
                    v.iso_values[i] = d.values[i];
                    for (int k = 0; k < 4; ++k) v.iso_colors[i][k] = d.colors[i].elem[k];
                    v.iso_densities[i] = d.optical_densities[i];
                }
                if (es.coloring == SurfaceColoring::AtomColors) {
                    v.color_volume = volume_take(b, gpu_volume_get(es.color_vol.gpu));
                } else if (es.coloring == SurfaceColoring::Field) {
                    v.color_field = volume_take(b, gpu_volume_get(es.field_vol.gpu));
                    if (v.color_field) {
                        v.color_map = colormap_texture(b, es.field_map.colormap, alloc);
                        if (!v.color_map) v.color_field = nullptr;
                        v.color_range[0] = es.field_map.range_beg;
                        v.color_range[1] = es.field_map.range_end;
                        const IsolineDesc il = color_scale_isoline_desc(es.field_isolines, px_per_point);
                        v.isoline_count = MD_MIN(il.count, (uint32_t)MD_COMPOSE_MAX_ISOLINES);
                        for (uint32_t i = 0; i < v.isoline_count; ++i) v.isoline_values[i] = il.values[i];
                        v.isoline_hover = il.hover;
                        v.isoline_hover_value = il.hover_value;
                        v.isoline_width = il.width;
                    }
                }
                v.normals = dv.normals;
            } else {
                md_float4_t planes[ES_MAX_PLANES];
                v.plane_count = (uint32_t)MD_MIN(es_volume_clip_planes(planes, es), MD_COMPOSE_MAX_PLANES);
                for (uint32_t i = 0; i < v.plane_count; ++i) {
                    for (int k = 0; k < 4; ++k) v.planes[i][k] = planes[i].elem[k];
                }
                v.transfer = tf_texture(b, es.volume.tf, alloc);
                v.transfer_range[0] = es.volume.tf.range_beg;
                v.transfer_range[1] = es.volume.tf.range_end;
                if (!v.transfer) *num_vols -= 1;
            }
            continue;
        }

        if (es.display != ElectronicStructureDisplay::Planes || !es.has_box) continue;
        // Each slice from its own image where it is of the plane as it is, else from the volume; its edge alone while
        // there is neither
        es_slices_conform(&es);
        const bool on_planes = es_plane_images_supported(state, es);
        md_gpu_texture_t transfer = nullptr;
        for (int i = 0; i < es.slices.count && i < ES_MAX_PLANES; ++i) {
            const ElectronicStructurePlane& p = es.slices.items[i];
            EsSliceFrame f;
            if (!es_slice_frame(&f, p)) continue;
            const md_float4_t edge = es_slice_edge_color(state, rep.node.id, i, p.enabled);
            float opacity = p.enabled ? MD_CLAMP(p.opacity, 0.0f, 1.0f) : 0.0f;
            GpuVolume field = {};
            if (opacity > 0.0f) {
                const ElectronicStructurePlaneImage& im = es.plane_images[i];
                if (on_planes && im.current) field = gpu_volume_get(im.img[im.front].gpu);
                if (!field.tex) field = dv;
                if (!transfer) transfer = tf_texture(b, es.volume.tf, alloc);
                if (!field.tex || !transfer) opacity = 0.0f;
            }
            if (opacity <= 0.0f && edge.w <= 0.0f) continue;
            if (*num_slices == MD_COMPOSE_MAX_SLICES) {
                dropped = true;
                continue;
            }
            md_compose_slice_t& s = slices[(*num_slices)++];
            s = {};
            mat_store(s.rect, md_float4x4_inverse(es_slice_rect_from_world(f)));
            if (opacity > 0.0f) {
                s.field = volume_take(b, field);
                mat_store(s.field_from_world, md_float4x4_inverse(field.texture_to_world));
                s.transfer = transfer;
                s.transfer_range[0] = es.volume.tf.range_beg;
                s.transfer_range[1] = es.volume.tf.range_end;
            }
            s.opacity = opacity;
            for (int k = 0; k < 4; ++k) s.edge[k] = edge.elem[k];
            s.edge_px = edge_px;
        }
    }
    if (dropped && !b->es_dropped_logged) {
        MD_LOG_INFO("The main view draws at most %d volumes and %d planes of the electronic structure representations: the rest are left out",
                    MD_COMPOSE_MAX_VOLUMES, MD_COMPOSE_MAX_SLICES);
    }
    b->es_dropped_logged = dropped;
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
    const uint32_t width = state->target.width, height = state->target.height;
    const uint32_t out_width = state->target.out_width ? state->target.out_width : width;
    const uint32_t out_height = state->target.out_height ? state->target.out_height : height;
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

    md_render_mol_bond_style_t styles[MD_RENDER_MOL_BOND_KIND_COUNT];     // a representation's, read by its emit
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
        mr.atom_color = r->color;
        mr.visible = r->visible;
        if (rep.type == RepresentationType::Bonds && !rep.sphere_fallback) {
            // The pairs its filter picks (md_render_bridge_rep_pairs), as thick as a licorice's sticks of the same scale
            if (pair_count(r->pairs) == 0) continue;
            mr.type = MD_RENDER_MOL_REP_PAIRS;
            mr.color.sharpness = rep.bond_color == BondColorMode::SmoothAtom ? 1.0f - MD_CLAMP(rep.bond_sharpness, 0.0f, 1.0f) : 1.0f;
            mr.chain_gaps = rep.chain_gaps;
            mr.params.pairs.radius = 0.2f * rep.scale.x;
            mr.pairs = r->pairs.pairs;
            mr.pair_id = r->pairs.ids;
            mr.pair_offset = r->pairs.offset;
            pair_styles(styles, rep);
            mr.bond_style = styles;
        } else if (!mol_rep_params(&mr, styles, b, rep, &r->profile, &r->profile_key, secondary)) {
            continue;
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
    emit_glyphs(b, state, alloc);
    emit_meshes(b, state, alloc);
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
    md_compose_volume_t vols[MD_COMPOSE_MAX_VOLUMES];
    md_compose_slice_t  slices[MD_COMPOSE_MAX_SLICES];
    uint32_t num_vols = 0, num_slices = 0;
    collect_electronic_structure(b, state, vols, &num_vols, slices, &num_slices);
    md_compose_desc_t cd = {};
    cd.camera = f.camera;
    cd.volumes = vols;
    cd.volume_count = num_vols;
    cd.slices = slices;
    cd.slice_count = num_slices;
    shading_of(&cd.shading, settings);
    cd.ao = settings.ssao.enabled ? MD_MAX(settings.ssao.intensity, 0.0f) : 0.0f;
    cd.ao_mode = settings.ssao.mode == postprocess_pipeline::SsaoMode_Quality ? MD_COMPOSE_AO_QUALITY : MD_COMPOSE_AO_PERFORMANCE;
    // The samples turned by the frame where the TAA averages them, else still (no flicker)
    cd.frame  = settings.taa.enabled ? (uint32_t)ImGui::GetFrameCount() : 0u;
    // the volumes' rays through the samples the raster's surfaces are of
    cd.jitter[0] = st.jitter[0];
    cd.jitter[1] = st.jitter[1];
    // The key light's shadows, found in the depth (md_compose): the light's size, their softness
    cd.shadows = settings.shadows.enabled;
    cd.light_radius = (float)MD_DEG_TO_RAD(MD_MAX(settings.shadows.softness, 0.0f));
    // Both accumulated over frames (md_compose's history), else their sample patterns held still. The history goes
    // with the filter: one made again starts over rather than from frames long past.
    if (settings.occlusion_temporal) {
        if (!b->vis_history) b->vis_history = md_compose_history_create();
    } else if (b->vis_history) {
        md_compose_history_destroy(b->stream, b->vis_history);
        b->vis_history = nullptr;
    }
    cd.history = b->vis_history;
    cd.target = b->scene;
    md_post_display_t dd = {};
    display_of(&dd, settings);
    dd.scene  = b->scene;
    dd.target = b->frame;
    bool ok = md_compose(b->stream, b->compose, b->view, &cd) && md_post_display(b->stream, b->post, &dd);
    luts_sweep(b);
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
        // The motion blur into the image shown; the history goes on without it
        const bool blur = settings.taa.motion_blur.enabled && settings.taa.motion_blur.motion_scale > 0.0f;
        rd.motion_blur = blur ? settings.taa.motion_blur.motion_scale : 0.0f;
        rd.shown = blur ? b->shown : nullptr;
        rd.seed = (uint32_t)ImGui::GetFrameCount();
        ok = md_post_resolve(b->stream, b->post, &rd);
        b->hist ^= 1;
        b->hist_valid = ok;
        b->image = ok ? (blur ? b->shown : b->history[b->hist]) : nullptr;
    } else {
        b->hist_valid = false;
        b->image = ok ? b->frame : nullptr;
    }
    if (!ok) {
        set_status(b, "The main view's frame could not be made (see the log)%s");
        return;
    }

    snprintf(b->status, sizeof(b->status), "Drawing %u batches, %u x %u, shown at %u x %u", f.batch_count, width, height, out_width, out_height);
}

// ---- The overlay -------------------------------------------------------------------------------------------------

// Over the image (md_post_overlay), of the camera the scene was drawn with and hidden behind its depth where asked;
// kept for picking it
// The overlay's items as md_post_overlay's, in alloc
static md_post_overlay_item_t* overlay_items(const OverlayItem* src, size_t n, md_allocator_i* alloc) {
    md_post_overlay_item_t* items = (md_post_overlay_item_t*)md_vm_arena_push(alloc, MD_MAX(n, (size_t)1) * sizeof(md_post_overlay_item_t));
    for (size_t i = 0; i < n; ++i) {
        const OverlayItem& it = src[i];
        md_post_overlay_item_t& x = items[i];
        x = {};
        x.kind   = it.shape == OverlayShape::Disc ? MD_POST_OVERLAY_DISC : it.shape == OverlayShape::Segment ? MD_POST_OVERLAY_SEGMENT : MD_POST_OVERLAY_TRIANGLE;
        x.flags  = (it.pixels ? (uint32_t)MD_POST_OVERLAY_PIXELS : 0u) | (it.depth ? (uint32_t)MD_POST_OVERLAY_DEPTH : 0u);
        md_memcpy(x.p0, &it.p[0], sizeof(x.p0));
        md_memcpy(x.p1, &it.p[1], sizeof(x.p1));
        md_memcpy(x.p2, &it.p[2], sizeof(x.p2));
        x.size   = it.size;
        x.width  = it.width;
        x.color  = it.color;
        x.color2 = it.color2;
    }
    return items;
}

static bool draw_overlay(MdRenderBridge* b, ApplicationState* state, const Overlay* o) {
    md_array_shrink(b->overlay, 0);
    const size_t n = o ? md_array_size(o->items) : 0;
    if (n == 0 || !b->has_last) return true;
    md_array_push_array(b->overlay, o->items, n, state->allocator.persistent);
    b->overlay_camera = b->last.camera;
    md_temp_scope_t temp = md_temp_begin_in(state->allocator.frame);
    md_post_overlay_item_t* items = overlay_items(o->items, n, state->allocator.frame);
    md_post_overlay_t ov = {};
    ov.target = b->out;
    ov.camera = b->overlay_camera;
    ov.items  = items;
    ov.count  = (uint32_t)n;
    ov.depth  = md_render_view_targets(b->view).depth;
    const bool ok = md_post_overlay(b->stream, b->post, &ov);
    md_temp_end(temp);
    return ok;
}

// The item of the overlay drawn last under the point (the image's pixels, rows from the top) of those that pick: the
// last drawn of them, as it is on top. depth: its first point's (md_render's, reverse Z).
static const OverlayItem* overlay_pick(const MdRenderBridge* b, float x, float y, float* depth) {
    const float* V = b->overlay_camera.view;
    const float* P = b->overlay_camera.proj;
    const float W = (float)b->out_width, H = (float)b->out_height;
    const md_float2_t q = {x, y};
    // a point's pixel, its clip w (> 0 in front) and depth
    auto project = [&](md_float3_t p, float* w, float* d) {
        float v[4], c[4];
        for (int r = 0; r < 4; ++r) v[r] = V[r] * p.x + V[4 + r] * p.y + V[8 + r] * p.z + V[12 + r];
        for (int r = 0; r < 4; ++r) c[r] = P[r] * v[0] + P[4 + r] * v[1] + P[8 + r] * v[2] + P[12 + r] * v[3];
        *w = c[3];
        *d = c[3] != 0.0f ? c[2] / c[3] : 0.0f;
        return md_float2_t{(c[0] / c[3] * 0.5f + 0.5f) * W, (0.5f - c[1] / c[3] * 0.5f) * H};
    };
    auto px_size = [&](const OverlayItem& it, float w) { return it.pixels ? it.size : it.size * P[5] * 0.5f * H / w; };
    for (size_t i = md_array_size(b->overlay); i-- > 0;) {
        const OverlayItem& it = b->overlay[i];
        if (!it.pick_domain) continue;
        float w0, d0, w1, d1, w2, d2;
        const md_float2_t a = project(it.p[0], &w0, &d0);
        if (!(w0 > 0.0f)) continue;
        bool in = false;
        if (it.shape == OverlayShape::Disc) {
            in = md_float2_length(md_float2_sub(q, a)) <= px_size(it, w0);
        } else if (it.shape == OverlayShape::Segment) {
            const md_float2_t e = project(it.p[1], &w1, &d1);
            if (!(w1 > 0.0f)) continue;
            const md_float2_t ab = md_float2_sub(e, a);
            const float l2 = md_float2_dot(ab, ab);
            const float t = l2 > 0.0f ? MD_CLAMP(md_float2_dot(md_float2_sub(q, a), ab) / l2, 0.0f, 1.0f) : 0.0f;
            const md_float2_t c = md_float2_add(a, md_float2_mul1(ab, t));
            in = md_float2_length(md_float2_sub(q, c)) <= px_size(it, w0) * (1.0f - t) + px_size(it, w1) * t;
        } else {
            const md_float2_t e = project(it.p[1], &w1, &d1), f = project(it.p[2], &w2, &d2);
            if (!(w1 > 0.0f) || !(w2 > 0.0f)) continue;
            auto side = [](md_float2_t p, md_float2_t u, md_float2_t v) { return (v.x - u.x) * (p.y - u.y) - (v.y - u.y) * (p.x - u.x); };
            const float s0 = side(q, a, e), s1 = side(q, e, f), s2 = side(q, f, a);
            in = (s0 >= 0.0f && s1 >= 0.0f && s2 >= 0.0f) || (s0 <= 0.0f && s1 <= 0.0f && s2 <= 0.0f);
        }
        if (in) {
            *depth = d0;
            return &it;
        }
    }
    return nullptr;
}

// md_post's marks of the settings' (marks_settings)
static void marks_of(md_post_marks_t* m, const postprocess_pipeline::Settings& s) {
    const auto& mk = s.marks;
    *m = {};
    for (int i = 0; i < 2; ++i) {
        for (int k = 0; k < 4; ++k) {
            m->visible[i][k] = mk.visible[i].elem[k];
            m->hidden[i][k]  = mk.hidden[i].elem[k];
            m->outline[i][k] = mk.outline[i].elem[k];
        }
    }
    for (int g = 0; g < 7; ++g) for (int k = 0; k < 3; ++k) m->groups[g][k] = mk.groups[g].elem[k];
    m->highlight  = mk.highlight;
    m->outline_px = (uint32_t)MD_MAX(mk.outline_step, 1);
}

bool md_render_bridge_present(ApplicationState* state, const postprocess_pipeline::Settings* marks, const Overlay* overlay) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!b || !b->ok || !b->image || !b->out) return false;
    // The marks changed without the scene being made: md_render's masks again (only they are made)
    if (marks && b->has_last && b->last_ui_frame != ImGui::GetFrameCount() && b->last.width == b->width && b->last.height == b->height) {
        if (take_up(b, state)) {
            glyph_marks(b, *state);
            md_render_frame_t f = b->last;
            for (int i = 0; i < MD_RENDER_CONTENT_COUNT; ++i) f.content[i] = b->content[i];
            if (md_render_view_draw(b->stream, b->view, &f, nullptr)) {
                b->last = f;
                b->last_ui_frame = ImGui::GetFrameCount();
            }
        }
    }
    md_post_marks_t m = {};
    if (marks) marks_of(&m, *marks);
    md_post_present_t pd = {};
    pd.image = b->image;
    pd.marks_view = marks && b->has_last ? b->view : nullptr;
    pd.marks  = marks ? &m : nullptr;
    pd.target = b->out;
    pd.frame  = b->presents++;
    return md_post_present(b->stream, b->post, &pd) && draw_overlay(b, state, overlay);
}

md_gpu_texture_t md_render_bridge_image(const ApplicationState* state) {
    const MdRenderBridge* b = state->md_render_bridge;
    return b && b->ok && b->image ? b->out : nullptr;
}

bool md_render_bridge_read_image(ApplicationState* state, uint32_t* rgba, uint32_t width, uint32_t height) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!b || !b->ok || !b->image || !b->out || b->out_width != width || b->out_height != height) return false;
    const size_t bytes = (size_t)width * height * 4;
    md_gpu_mem_t rb = md_gpu_malloc(b->stream, MD_GPU_MEM_HOST_READ, bytes);
    const bool ok = rb.cpu && md_gpu_copy_from_texture(b->stream, rb.gpu, b->out, nullptr);
    md_gpu_stream_sync(b->stream);
    if (ok) memcpy(rgba, rb.cpu, bytes);
    if (rb.gpu) md_gpu_free(b->stream, rb.gpu);
    return ok;
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
    if (!b || !b->ok) return;
    b->pick_x = x;
    b->pick_y = y;
    if (!b->has_last) return;
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
struct MdRenderPanels;
static bool panel_rep_transform(const MdRenderPanels* p, uint64_t rep_key, md_float4x4_t* out);

// The transform of the copy a hit is in: of the main view's copies (panels null, emit_placed), or of a panel's
// representation
static bool hit_transform(const MdRenderBridge* b, const MdRenderPanels* panels, const ApplicationState& state, const md_render_mol_hit_t& mh,
                          md_float4x4_t* out) {
    return panels ? panel_rep_transform(panels, mh.rep_key, out) : copy_transform(b, state, mh.rep_key, mh.instance, out);
}

static bool pair_middle(const MdRenderBridge* b, const MdRenderPanels* panels, const ApplicationState& state, const md_render_mol_hit_t& mh, md_atom_pair_t pair) {
    const md_float3_t* xyz = state.mold.state.xyz;
    const md_float3_t p = {mh.position[0], mh.position[1], mh.position[2]};
    md_float3_t a = xyz[pair.idx[0]], e = xyz[pair.idx[1]];
    md_float4x4_t T;
    if (hit_transform(b, panels, state, mh, &T)) {
        a = md_float3_from_float4(md_float4x4_mul_float4(T, md_float4_from_float3(a, 1.0f)));
        e = md_float3_from_float4(md_float4x4_mul_float4(T, md_float4_from_float3(e, 1.0f)));
    }
    const md_float3_t ab = e - a;
    const float len2 = md_float3_dot(ab, ab);
    const float t = len2 > 0.0f ? md_float3_dot(p - a, ab) / len2 : 0.5f;
    return t >= 0.25f && t <= 0.75f;
}

// The hit of an answer of md_render's on what md_render_mol drew (an atom, a bond, a segment of the cartoon, a ring, a
// gap), in their domains as md_gl's were; nothing of a system no more. Of the main view's representations (panels null)
// or a panel set's.
static void mol_hit(const MdRenderBridge* b, const MdRenderPanels* panels, ApplicationState* state, const md_render_pick_t& found, PickingHit* hit) {
    md_render_mol_hit_t mh = {};
    // Of a system no more (or not known): nothing, rather than what md_render drew there taken for something else
    if (!b->has_topo || !md_render_mol_decode_pick(&found, &mh) || mh.kind == MD_RENDER_MOL_HIT_NONE) return;
    const md_system_t& sys = state->mold.sys;
    switch (mh.kind) {
    case MD_RENDER_MOL_HIT_BOND: {
        // As md_gl's: the middle half of a bond picks as the bond, the quarters at its ends as their atoms. An
        // interaction over the sticks of a licorice (past the system's bonds), of the representation's range of them.
        if (mh.bond >= sys.bond.count) {
            set_hit(*state, hit, PickingDomain_Bond, (mh.rep_key >> 16) - 1, (uint32_t)(mh.bond - sys.bond.count));
        } else if (pair_middle(b, panels, *state, mh, sys.bond.pairs[mh.bond])) {
            set_hit(*state, hit, PickingDomain_Bond, 0, mh.bond);
        } else {
            set_hit(*state, hit, PickingDomain_Atom, 0, mh.atom);
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
        if (pair_middle(b, panels, *state, mh, pair)) {
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

bool md_render_bridge_hit(ApplicationState* state, PickingHit* hit) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!hit || !b || !b->ok) return false;
    // The overlay first, as it is over the scene: as it is drawn now, of the pixel asked about now
    if (b->pick_x >= 0.0f && md_array_size(b->overlay) && b->width && b->height) {
        float d = 0.0f;
        const float sx = (float)b->out_width / (float)b->width, sy = (float)b->out_height / (float)b->height;
        if (const OverlayItem* it = overlay_pick(b, b->pick_x * sx, b->pick_y * sy, &d)) {
            *hit = {};
            hit->frame_idx     = state->picking_handler.frame_idx;
            hit->surface_coord = {b->pick_x, (float)b->height - b->pick_y};
            hit->depth         = 1.0f - d;
            set_hit(*state, hit, it->pick_domain, it->pick_key, it->pick_index);
            return true;
        }
    }
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
        } else if (what == SHAPE_GLYPHS) {
            // A glyph as its point, of its representation's range
            set_hit(*state, hit, PickingDomain_Glyph, ((found->key >> 16) & 0xFFFFFFFFFFull) - 1, found->element);
        }
        return true;     // the cells' edges and the meshes as nothing to pick
    }
    mol_hit(b, nullptr, state, *found, hit);
    return true;
}

void md_render_bridge_buffers(ApplicationState* state, uint32_t dirty) {
    MdRenderBridge* b = state->md_render_bridge;
    if (b) b->pending |= dirty;
}

// ---- Panels ------------------------------------------------------------------------------------------------------

static constexpr uint32_t PANEL_MAX      = 16;
static constexpr uint32_t PANEL_SAMPLES  = 16;               // a still panel's: it settles over as many frames
static constexpr uint64_t PANEL_REP_BASE = 0xFF0001;         // a representation's key: (base + its index) << 16, of no
                                                             // representation of the main view's, of no part
static constexpr uint64_t PANEL_MESH     = 0xFEull << 56;    // the meshes' (their low byte none of md_render_mol's)
static_assert(PANEL_TRANSFER_RES == LUT_RES, "a transfer function's texels");

struct MdRenderPanels {
    uint32_t width = 0, height = 0, count = 0, columns = 0;     // as made
    bool     shared = true;                     // the panels in one view (view[0]), else each in its own
    md_render_view_t view[PANEL_MAX] = {};
    uint64_t content[PANEL_MAX][MD_RENDER_CONTENT_COUNT] = {};  // a view's own counters, past the system's
    // Each panel's image as it settles: its histories (the one written last, and whether it holds the panel), what its
    // surfaces were when it was made
    md_gpu_texture_t history[PANEL_MAX][2] = {};
    int      hist[PANEL_MAX] = {};
    bool     settled[PANEL_MAX] = {};
    uint64_t surfaces[PANEL_MAX] = {};
    // A panel's on the way: composed, displayed, presented (in turn, one panel after another)
    md_gpu_texture_t scene = nullptr, frame = nullptr, image = nullptr;
    // Their images side by side, which ImGui draws
    md_gpu_texture_t atlas = nullptr;
    uint32_t  atlas_w = 0, atlas_h = 0;
    bool      has_image = false;
    // The molecule: its sets of colours, and of each representation which atoms it shows (of its set's, and of its own),
    // its transform; what they were made of. The cartoons' profiles, by what they are of, kept while drawn.
    struct Colors { md_gpu_addr_t color; uint32_t cap; uint64_t key; md_array(uint64_t) shown; };
    struct Rep    { md_gpu_addr_t visible; uint32_t cap; uint64_t key; md_float4x4_t xf; bool has_xf; };
    struct Profile { uint64_t key; md_gpu_addr_t addr; uint64_t made; bool used; };
    md_array(Colors)  colors = nullptr;
    md_array(Rep)     reps = nullptr;
    md_array(Profile) profiles = nullptr;
    uint32_t atom_count = 0;
    uint64_t topology = 0;                      // the system's they are of (its topology_version)
    uint64_t look_key = 0, present_key = 0;     // the shading and display; the marks and the overlay
    bool     too_many_logged = false;
    md_array(md_render_batch_t) batches = nullptr;
    md_float4x4_t mesh_model = {};              // the meshes' bohr to the world
    // The view asked about last (PANEL_MAX: none) and its newest answer
    uint32_t pick_view = PANEL_MAX;
    bool     has_pick = false;
    md_render_pick_t pick = {};
};

// The transform of a panel representation's hit (its key's), if it has one
static bool panel_rep_transform(const MdRenderPanels* p, uint64_t rep_key, md_float4x4_t* out) {
    const uint64_t i = (rep_key >> 16) - PANEL_REP_BASE;
    if ((rep_key >> 16) < PANEL_REP_BASE || i >= md_array_size(p->reps) || !p->reps[i].has_xf) return false;
    *out = p->reps[i].xf;
    return true;
}

// The views but the first (each panel's own), gone
static void panels_free_views(MdRenderBridge* b, MdRenderPanels* p, uint32_t first) {
    for (uint32_t v = first; v < PANEL_MAX; ++v) {
        if (p->view[v]) md_render_view_destroy(b->stream, p->view[v]);
        p->view[v] = nullptr;
        if (p->pick_view == v) p->pick_view = PANEL_MAX;
    }
}

static void panels_free_targets(MdRenderBridge* b, MdRenderPanels* p) {
    for (uint32_t i = 0; i < PANEL_MAX; ++i) {
        for (int k = 0; k < 2; ++k) {
            if (p->history[i][k]) md_gpu_texture_destroy(p->history[i][k]);
            p->history[i][k] = nullptr;
        }
        p->settled[i] = false;
        p->surfaces[i] = 0;
    }
    md_gpu_texture_t* t[] = {&p->scene, &p->frame, &p->image};
    for (md_gpu_texture_t* x : t) {
        if (*x) md_gpu_texture_destroy(*x);
        *x = nullptr;
    }
    if (p->atlas) md_gpu_texture_destroy(p->atlas);
    p->atlas = nullptr;
    p->has_image = false;
    p->width = p->height = p->count = p->columns = p->atlas_w = p->atlas_h = 0;
}

// The textures of panels of the size, so many, in so many columns
static bool panels_ensure(MdRenderBridge* b, MdRenderPanels* p, uint32_t w, uint32_t h, uint32_t count, uint32_t columns) {
    if (p->scene && p->width == w && p->height == h && p->count == count && p->columns == columns) return true;
    panels_free_targets(b, p);
    const md_gpu_tex_usage_t SS = MD_GPU_TEX_STORAGE | MD_GPU_TEX_SAMPLED;
    p->scene = gpu_texture(b, MD_GPU_FORMAT_RGBA16_FLOAT, SS, w, h, "panel scene");
    p->frame = gpu_texture(b, MD_GPU_FORMAT_RGBA16_FLOAT, SS, w, h, "panel frame");
    p->image = gpu_texture(b, MD_GPU_FORMAT_RGBA8_UNORM, MD_GPU_TEX_STORAGE | MD_GPU_TEX_RENDER_TARGET, w, h, "panel image");
    bool ok = p->scene && p->frame && p->image;
    for (uint32_t i = 0; i < count && ok; ++i) {
        for (int k = 0; k < 2; ++k) {
            p->history[i][k] = gpu_texture(b, MD_GPU_FORMAT_RGBA16_FLOAT, SS, w, h, "panel history");
            ok = ok && p->history[i][k];
        }
    }
    const uint32_t rows = (count + columns - 1) / columns;
    p->atlas = ok ? gpu_texture(b, MD_GPU_FORMAT_RGBA8_UNORM, MD_GPU_TEX_SAMPLED, w * columns, h * rows, "panels") : nullptr;
    ok = ok && p->atlas;
    if (!ok) {
        MD_LOG_ERROR("The panels' textures (%u of %u x %u) could not be made: %s", count, w, h, md_gpu_last_error());
        panels_free_targets(b, p);
        return false;
    }
    p->atlas_w = w * columns;
    p->atlas_h = h * rows;
    p->width = w;
    p->height = h;
    p->count = count;
    p->columns = columns;
    // The images begin again from their views' first samples, which they would not if their rasters had settled
    for (uint32_t v = 0; v < PANEL_MAX; ++v) p->content[v][MD_RENDER_CONTENT_APPEARANCE] += 1;
    return true;
}

// The molecule's sets of colours (taken when they changed: then also which atoms they show), and each representation's
// atoms shown (of its set's, and of its own), uploaded when they changed. False when they cannot be.
static bool panels_molecule(MdRenderBridge* b, MdRenderPanels* p, const PanelsDesc& d, md_allocator_i* alloc, md_allocator_i* temp_alloc) {
    const uint32_t n = b->atom_count;
    const size_t words = ((size_t)n + 63) / 64;
    // Of another system: all of it again, and the profiles (of its topology's size) gone
    if (p->atom_count != n || p->topology != b->topo.topology_version) {
        for (size_t i = 0; i < md_array_size(p->colors); ++i) p->colors[i].key = 0;
        for (size_t i = 0; i < md_array_size(p->reps); ++i) p->reps[i].key = 0;
        for (size_t i = 0; i < md_array_size(p->profiles); ++i) if (p->profiles[i].addr) md_gpu_free(b->stream, p->profiles[i].addr);
        md_array_shrink(p->profiles, 0);
        p->atom_count = n;
        p->topology = b->topo.topology_version;
    }
    // A set's colours, and its atoms not transparent
    while (md_array_size(p->colors) < d.color_count) md_array_push(p->colors, (MdRenderPanels::Colors{}), alloc);
    for (uint32_t c = 0; c < d.color_count; ++c) {
        MdRenderPanels::Colors& cs = p->colors[c];
        const uint32_t* src = d.colors[c];
        if (!src) return false;
        const uint64_t key = md_hash64(src, (size_t)n * 4, 0xc010125) | 1;
        if (key == cs.key) continue;
        if (cs.cap < n) {
            cs.color = realloc_gpu(b, cs.color, (size_t)n * 4);
            cs.cap = cs.color ? n : 0;
            if (!cs.color) return false;
        }
        md_gpu_upload(b->stream, cs.color, src, (size_t)n * 4);
        md_array_resize(cs.shown, words, alloc);
        md_memset(cs.shown, 0, words * 8);
        for (uint32_t i = 0; i < n; ++i) {
            if (src[i] >> 24) cs.shown[i >> 6] |= 1ull << (i & 63);
        }
        cs.key = key;
        for (uint32_t v = 0; v < PANEL_MAX; ++v) p->content[v][MD_RENDER_CONTENT_APPEARANCE] += 1;
    }
    // A representation's atoms: of its set's, those of its own
    while (md_array_size(p->reps) > d.rep_count) {
        if (md_array_last(p->reps)->visible) md_gpu_free(b->stream, md_array_last(p->reps)->visible);
        md_array_pop(p->reps);
    }
    while (md_array_size(p->reps) < d.rep_count) md_array_push(p->reps, (MdRenderPanels::Rep{}), alloc);
    for (uint32_t r = 0; r < d.rep_count; ++r) {
        const PanelRep& pr = d.reps[r];
        MdRenderPanels::Rep& rp = p->reps[r];
        rp.has_xf = pr.transform != nullptr;
        rp.xf = pr.transform ? *pr.transform : md_float4x4_ident();
        if (pr.colors >= d.color_count) return false;
        const MdRenderPanels::Colors& cs = p->colors[pr.colors];
        uint64_t key = md_hash64(&cs.key, sizeof(cs.key), 0x815181e);
        if (pr.atoms) key = md_bitset_hash64(*pr.atoms, key);
        key |= 1;
        if (key == rp.key && rp.visible) continue;
        if (rp.cap < n) {
            rp.visible = realloc_gpu(b, rp.visible, bitset_words32(n) * 4);
            rp.cap = rp.visible ? n : 0;
            if (!rp.visible) return false;
        }
        md_temp_scope_t temp = md_temp_begin_in(temp_alloc);
        uint64_t* bits = (uint64_t*)md_vm_arena_push(temp_alloc, MD_MAX(words, (size_t)1) * 8);
        for (size_t k = 0; k < words; ++k) {
            bits[k] = cs.shown[k] & (pr.atoms ? md_bitset_word(*pr.atoms, (uint32_t)(k * 64)) : ~0ull);
        }
        md_gpu_upload(b->stream, rp.visible, bits, words * 8);
        md_temp_end(temp);
        rp.key = key;
        for (uint32_t v = 0; v < PANEL_MAX; ++v) p->content[v][MD_RENDER_CONTENT_GEOMETRY] += 1;
    }
    return true;
}

// A cartoon's profile for the panels: one of those kept, of what it is of, made when none is
static md_gpu_addr_t panel_profile(MdRenderBridge* b, MdRenderPanels* p, const Representation& rep, uint64_t secondary, md_allocator_i* alloc,
                                   md_render_mol_rep_t* mr, md_render_mol_bond_style_t styles[MD_RENDER_MOL_BOND_KIND_COUNT]) {
    const float scale[3] = {rep.scale.x, rep.scale.z, rep.scale.y};
    const uint64_t key = md_hash64(scale, sizeof(scale), secondary) | 1;     // cartoon_profile's
    size_t at = md_array_size(p->profiles);
    for (size_t i = 0; i < md_array_size(p->profiles); ++i) {
        if (p->profiles[i].key == key) at = i;
    }
    if (at == md_array_size(p->profiles)) md_array_push(p->profiles, (MdRenderPanels::Profile{key, 0, 0, false}), alloc);
    MdRenderPanels::Profile& pf = p->profiles[at];
    pf.used = true;
    mol_rep_params(mr, styles, b, rep, &pf.addr, &pf.made, secondary);
    return pf.addr;
}

// A panel volume's transfer function's lookup table's key (lut_texture)
static uint64_t panel_transfer_key(const uint32_t* transfer) {
    return md_hash64(transfer, sizeof(uint32_t) * PANEL_TRANSFER_RES, 0x7a2f);
}

// Whether a panel's meshes are there to draw
static bool panel_has_meshes(const PanelSurfaces& s) {
    const IsoMeshes* m = s.meshes;
    return m && m->gpu_xyz && m->gpu_normal && m->gpu_index && m->num_vertices && s.sheet_count && s.sheets;
}

// What a panel's surfaces are: its image made again when it changed
static uint64_t panel_surfaces_key(const PanelSurfaces& s) {
    uint64_t k = 0x5afe5;
    const uint32_t nv = s.volumes ? MD_MIN(s.volume_count, (uint32_t)MD_COMPOSE_MAX_VOLUMES) : 0;
    for (uint32_t j = 0; j < nv; ++j) {
        const PanelVolume& v = s.volumes[j];
        const uint64_t ver = v.volume ? gpu_volume_version(v.volume) : 0;
        k = md_hash64(&v.volume, sizeof(v.volume), k);
        k = md_hash64(&ver, sizeof(ver), k);
        const uint32_t ni = v.iso_values && v.iso_colors ? MD_MIN(v.iso_count, (uint32_t)MD_COMPOSE_MAX_ISO) : 0;
        k = md_hash64(&ni, sizeof(ni), k);
        if (ni) {
            k = md_hash64(v.iso_values, sizeof(float) * ni, k);
            k = md_hash64(v.iso_colors, sizeof(md_float4_t) * ni, k);
        }
        if (v.transfer) k = md_hash64(v.transfer, sizeof(uint32_t) * PANEL_TRANSFER_RES, k);
        const float box[8] = {v.transfer_range[0], v.transfer_range[1], v.clip_min.x, v.clip_min.y, v.clip_min.z, v.clip_max.x, v.clip_max.y, v.clip_max.z};
        k = md_hash64(box, sizeof(box), k);
    }
    if (panel_has_meshes(s)) {
        k = md_hash64(&s.meshes->version, sizeof(s.meshes->version), k);
        for (uint32_t j = 0; j < s.sheet_count; ++j) {
            const PanelSheet& sh = s.sheets[j];
            const float f[6] = {(float)sh.index_offset, (float)sh.index_count, sh.color.x, sh.color.y, sh.color.z, sh.color.w};
            k = md_hash64(f, sizeof(f), k);
        }
    }
    return k | 1;
}

// A panel's meshes, after the batches in *out: a sheet a batch of triangles, in its colour and coverage
static void panel_meshes(MdRenderBridge* b, MdRenderPanels* p, const PanelSurfaces& s, md_render_batch_t** out, md_allocator_i* alloc) {
    if (!panel_has_meshes(s)) return;
    const IsoMeshes& m = *s.meshes;
    md_gpu_stream_wait(b->stream, m.gpu_ready);
    for (uint32_t j = 0; j < s.sheet_count; ++j) {
        const PanelSheet& sh = s.sheets[j];
        const float alpha = MD_CLAMP(sh.color.w, 0.0f, 1.0f);
        if (alpha <= 0.0f || sh.index_count < 3 || sh.index_offset + sh.index_count > m.num_indices) continue;
        md_render_batch_t bt = {};
        bt.key = PANEL_MESH | ((uint64_t)j << 8) | 0xFFu;
        bt.kind = MD_RENDER_TRIANGLES;
        bt.count = sh.index_count / 3;
        bt.points.count = m.num_vertices;
        bt.points.position = m.gpu_xyz;
        bt.points.stride = 16;
        bt.index = m.gpu_index + (md_gpu_addr_t)sh.index_offset * sizeof(uint32_t);
        bt.style.triangles.normal = m.gpu_normal;
        bt.style.triangles.normal_stride = 16;
        bt.color.mode = MD_RENDER_COLOR_UNIFORM;
        bt.color.uniform = rgba8_of(md_float4_set(sh.color.x, sh.color.y, sh.color.z, 1.0f));
        bt.color.transparency = alpha < 1.0f ? 1.0f - alpha : 0.0f;
        bt.transform = &p->mesh_model.elem[0][0];
        bt.transform_count = 1;
        md_array_push(*out, bt, alloc);
    }
}

// Panel i composed over the raster of its view as it is (st), displayed and settled into its history: a sample more,
// or the first
static bool panel_make(MdRenderBridge* b, MdRenderPanels* p, uint32_t i, md_render_view_t view, const md_render_camera_t& cam,
                       const md_render_status_t& st, const PanelSurfaces& s, const md_compose_shading_t& sh, const md_post_display_t& disp,
                       md_allocator_i* alloc) {
    md_compose_volume_t vols[MD_COMPOSE_MAX_VOLUMES] = {};
    uint32_t nv = 0;
    for (uint32_t j = 0; s.volumes && j < s.volume_count && nv < MD_COMPOSE_MAX_VOLUMES; ++j) {
        const PanelVolume& pv = s.volumes[j];
        const uint32_t ni = pv.iso_values && pv.iso_colors ? MD_MIN(pv.iso_count, (uint32_t)MD_COMPOSE_MAX_ISO) : 0;
        const GpuVolume gv = pv.volume ? gpu_volume_get(pv.volume) : GpuVolume{};
        if (!gv.tex || (!ni && !pv.transfer)) continue;
        md_compose_volume_t& vol = vols[nv];
        vol.field = volume_take(b, gv);
        mat_store(vol.model, gv.texture_to_world);
        for (int k = 0; k < 3; ++k) {
            vol.clip_min[k] = MD_CLAMP(pv.clip_min.elem[k], 0.0f, 1.0f);
            vol.clip_max[k] = MD_CLAMP(pv.clip_max.elem[k], 0.0f, 1.0f);
        }
        vol.iso_count = ni;
        for (uint32_t q = 0; q < ni; ++q) {
            vol.iso_values[q] = pv.iso_values[q];
            for (int k = 0; k < 4; ++k) vol.iso_colors[q][k] = pv.iso_colors[q].elem[k];
        }
        vol.normals = gv.normals;
        if (pv.transfer) {
            vol.transfer = lut_texture(b, panel_transfer_key(pv.transfer), pv.transfer, alloc);
            vol.transfer_range[0] = pv.transfer_range[0];
            vol.transfer_range[1] = pv.transfer_range[1];
            if (!vol.transfer && !ni) continue;
        }
        nv += 1;
    }
    md_compose_desc_t cd = {};
    cd.camera = cam;
    cd.shading = sh;
    cd.volumes = vols;
    cd.volume_count = nv;
    cd.frame = st.sample;
    cd.jitter[0] = st.jitter[0];
    cd.jitter[1] = st.jitter[1];
    cd.target = p->scene;
    md_post_display_t dd = disp;
    dd.scene = p->scene;
    dd.target = p->frame;
    // The mean of the samples: the history of those before weighed by their number
    const bool go_on = st.sample > 0 && p->settled[i];
    md_post_resolve_t rd = {};
    rd.view = view;
    rd.camera = cam;
    rd.jitter[0] = st.jitter[0];
    rd.jitter[1] = st.jitter[1];
    rd.frame = p->frame;
    rd.history_prev = go_on ? p->history[i][p->hist[i]] : nullptr;
    rd.history = p->history[i][p->hist[i] ^ 1];
    rd.feedback_min = rd.feedback_max = go_on ? (float)st.sample / (float)(st.sample + 1) : 0.0f;
    const bool ok = md_compose(b->stream, b->compose, view, &cd) && md_post_display(b->stream, b->post, &dd) && md_post_resolve(b->stream, b->post, &rd);
    p->hist[i] ^= 1;
    p->settled[i] = ok;
    return ok;
}

// Panel i's image, with the marks of its view over it (null: none) and the overlay (ov's: its items and camera; null:
// none), into its place in the atlas
static bool panel_present(MdRenderBridge* b, MdRenderPanels* p, uint32_t i, md_render_view_t view, const md_post_marks_t* marks,
                          md_post_overlay_t* ov) {
    md_post_present_t pd = {};
    pd.image = p->history[i][p->hist[i]];
    pd.marks_view = marks ? view : nullptr;
    pd.marks = marks;
    pd.target = p->image;
    pd.frame = b->presents++;
    if (!md_post_present(b->stream, b->post, &pd)) return false;
    if (ov) {
        ov->target = p->image;
        ov->depth  = md_render_view_targets(view).depth;
        if (!md_post_overlay(b->stream, b->post, ov)) return false;
    }
    md_gpu_tex_region_t src = {};
    src.extent[0] = p->width;
    src.extent[1] = p->height;
    src.extent[2] = 1;
    md_gpu_tex_region_t dst = {};
    dst.offset[0] = (i % p->columns) * p->width;
    dst.offset[1] = (i / p->columns) * p->height;
    dst.extent[0] = p->width;
    dst.extent[1] = p->height;
    dst.extent[2] = 1;
    return md_gpu_copy_texture(b->stream, p->atlas, &dst, p->image, &src);
}

bool md_render_bridge_panels_draw(ApplicationState* state, MdRenderPanels** pp, const PanelsDesc& d) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!b || !b->ok || !pp) return false;
    // The system as the main view took it up: its topology, positions and marks
    const md_system_t& sys = state->mold.sys;
    if (!b->has_topo || !b->updated || b->atom_count != sys.atom.count || b->topo.topology_version != sys.topology_version) return false;
    const uint32_t count = MD_MIN(d.count, PANEL_MAX);
    if (count == 0 || d.width == 0 || d.height == 0 || !d.settings || (d.rep_count && !d.reps)) return false;
    if (!*pp) *pp = new MdRenderPanels();
    MdRenderPanels* p = *pp;
    if (!panels_ensure(b, p, d.width, d.height, count, MD_CLAMP(d.columns, 1u, count))) return false;
    md_allocator_i* alloc = state->allocator.persistent;
    md_allocator_i* temp_alloc = state->allocator.frame;
    if (!panels_molecule(b, p, d, alloc, temp_alloc)) return false;

    // The look: when it changed every panel is made again; the marks and the overlay: presented again
    md_compose_shading_t sh;
    shading_of(&sh, *d.settings);
    md_post_display_t disp = {};
    display_of(&disp, *d.settings);
    {
        const float look[] = {(float)sh.mode, (float)sh.linear_colors, sh.light_dir[0], sh.light_dir[1], sh.light_dir[2],
                              sh.dir_radiance[0], sh.dir_radiance[1], sh.dir_radiance[2], sh.sky_radiance[0], sh.sky_radiance[1],
                              sh.sky_radiance[2], sh.ground_radiance[0], sh.ground_radiance[1], sh.ground_radiance[2],
                              sh.light_strength, sh.roughness, sh.f0, sh.ambient, sh.specular, (float)disp.curve, disp.exposure,
                              disp.gamma, disp.background[0], disp.background[1], disp.background[2]};
        const uint64_t key = md_hash64(look, sizeof(look), 0x100c) | 1;
        if (key != p->look_key) {
            p->look_key = key;
            for (uint32_t v = 0; v < PANEL_MAX; ++v) p->content[v][MD_RENDER_CONTENT_APPEARANCE] += 1;
        }
    }
    md_post_marks_t mk = {};
    if (d.marks) marks_of(&mk, *d.settings);
    md_render_camera_t cam;
    camera_from(&cam, d.view, d.proj);
    uint64_t present_key = md_hash64(&mk, sizeof(mk), d.marks ? 2 : 1);
    if (d.overlay_count) {
        present_key = md_hash64(d.overlay, sizeof(OverlayItem) * d.overlay_count, present_key);
        present_key = md_hash64(&cam, sizeof(cam), present_key);
    }
    const bool present_changed = present_key != p->present_key;
    p->present_key = present_key;

    // The molecule: each representation, of its atoms in its set's colours (the same marks as the main view's)
    md_array_shrink(p->batches, 0);
    for (size_t i = 0; i < md_array_size(p->profiles); ++i) p->profiles[i].used = false;
    const md_render_mol_atoms_t atoms = atoms_of(b, state);
    const uint64_t secondary = update_secondary(b, state);
    for (uint32_t r = 0; r < d.rep_count; ++r) {
        const PanelRep& pr = d.reps[r];
        if (!pr.rep) continue;
        // md_render takes at most 4095 batches a frame; a representation emits up to about a dozen
        if (md_array_size(p->batches) > 4000) {
            if (!p->too_many_logged) MD_LOG_INFO("A view draws %u representations of the molecule: those past %u are left out", d.rep_count, r);
            p->too_many_logged = true;
            break;
        }
        md_render_mol_rep_t mr = {};
        md_render_mol_bond_style_t styles[MD_RENDER_MOL_BOND_KIND_COUNT];
        if (pr.rep->type == RepresentationType::Cartoon && !pr.rep->sphere_fallback) {
            panel_profile(b, p, *pr.rep, secondary, alloc, &mr, styles);
        } else {
            md_gpu_addr_t none = 0;
            uint64_t none_key = 0;
            if (!mol_rep_params(&mr, styles, b, *pr.rep, &none, &none_key, secondary)) continue;
        }
        mr.key = (PANEL_REP_BASE + r) << 16;
        mr.color.mode = MD_RENDER_COLOR_POINTS;
        mr.atom_color = p->colors[pr.colors].color;
        mr.visible = p->reps[r].visible;
        if (pr.transform) {
            mr.transform = &p->reps[r].xf.elem[0][0];
            mr.transform_count = 1;
        }
        emit(b, &p->batches, atoms, mr, alloc);
    }
    for (size_t i = 0; i < md_array_size(p->profiles);) {
        if (!p->profiles[i].used) {
            if (p->profiles[i].addr) md_gpu_free(b->stream, p->profiles[i].addr);
            md_array_swap_back_and_pop(p->profiles, i);
        } else {
            ++i;
        }
    }
    const size_t mol_count = md_array_size(p->batches);

    // One view for all of them, unless they have meshes of their own
    bool shared = true;
    for (uint32_t i = 0; i < count; ++i) shared = shared && !panel_has_meshes(d.panels[i]);
    if (shared != p->shared) {
        for (uint32_t i = 0; i < PANEL_MAX; ++i) p->settled[i] = false;
        p->shared = shared;
    }
    panels_free_views(b, p, shared ? 1 : count);    // those not drawn (the first goes on as panel 0's)
    p->mesh_model = iso_meshes_model_matrix();
    // A panel whose surfaces changed made again from its first sample, with the others of its view
    for (uint32_t i = 0; i < count; ++i) {
        const uint64_t key = panel_surfaces_key(d.panels[i]);
        if (key == p->surfaces[i]) continue;
        p->surfaces[i] = key;
        uint64_t* c = p->content[shared ? 0 : i];
        c[MD_RENDER_CONTENT_APPEARANCE] += 1;
        if (!shared) c[MD_RENDER_CONTENT_GEOMETRY] += 1;
    }
    // Their transfer functions kept while they are shown (luts_sweep), whether made again this frame or not
    for (uint32_t i = 0; i < count; ++i) {
        const PanelSurfaces& s = d.panels[i];
        for (uint32_t j = 0; s.volumes && j < s.volume_count; ++j) {
            if (s.volumes[j].transfer) lut_texture(b, panel_transfer_key(s.volumes[j].transfer), nullptr, alloc);
        }
    }

    // The overlay, over each panel presented
    md_temp_scope_t ov_temp = md_temp_begin_in(temp_alloc);
    md_post_overlay_t ov = {};
    if (d.overlay_count && d.overlay) {
        ov.camera = cam;
        ov.items = overlay_items(d.overlay, d.overlay_count, temp_alloc);
        ov.count = d.overlay_count;
    }

    const uint32_t num_views = shared ? 1 : count;
    bool presented = false, ok = true;
    for (uint32_t v = 0; v < num_views && ok; ++v) {
        md_render_batch_t* batches = p->batches;
        size_t num_batches = mol_count;
        md_temp_scope_t temp = md_temp_begin_in(temp_alloc);
        if (!shared) {
            md_array(md_render_batch_t) own = nullptr;
            md_array_push_array(own, p->batches, mol_count, temp_alloc);
            panel_meshes(b, p, d.panels[v], &own, temp_alloc);
            batches = own;
            num_batches = md_array_size(own);
        }
        if (!p->view[v]) p->view[v] = md_render_view_create(b->render, shared ? "panels" : "panel");
        md_render_frame_t f = {};
        f.width = p->width;
        f.height = p->height;
        f.camera = cam;
        f.batches = batches;
        f.batch_count = (uint32_t)num_batches;
        f.jitter = MD_RENDER_JITTER_STATIC;
        f.samples = PANEL_SAMPLES;
        f.time = b->time;
        for (int k = 0; k < MD_RENDER_CONTENT_COUNT; ++k) f.content[k] = b->content[k] + p->content[v][k];
        md_render_status_t st = {};
        ok = p->view[v] && md_render_view_draw(b->stream, p->view[v], &f, &st);
        md_temp_end(temp);
        if (!ok) break;
        const bool made = (st.work & (MD_RENDER_WORK_RASTER | MD_RENDER_WORK_RESOLVE)) != 0;
        const bool marked = (st.work & MD_RENDER_WORK_MASK) != 0 || present_changed || !p->has_image;
        const uint32_t i0 = shared ? 0 : v, i1 = shared ? count : v + 1;
        for (uint32_t i = i0; i < i1 && ok; ++i) {
            const bool make = made || !p->settled[i];
            if (!make && !marked) continue;
            if (make) ok = panel_make(b, p, i, p->view[v], cam, st, d.panels[i], sh, disp, alloc);
            if (!ok) break;
            presented = true;
            ok = panel_present(b, p, i, p->view[v], d.marks ? &mk : nullptr, ov.count ? &ov : nullptr);
        }
    }
    md_temp_end(ov_temp);
    if (!ok) {
        MD_LOG_ERROR("The panels could not be drawn (see the log)");
        return false;
    }
    if (presented) p->has_image = true;
    return p->has_image;
}

bool md_render_bridge_panels_image(const MdRenderPanels* p, uint32_t i, md_gpu_texture_t* texture, md_float2_t* uv0, md_float2_t* uv1) {
    if (!p || !p->has_image || i >= p->count) return false;
    const float W = (float)p->atlas_w, H = (float)p->atlas_h;
    const float x0 = (float)((i % p->columns) * p->width), y0 = (float)((i / p->columns) * p->height);
    *texture = p->atlas;
    *uv0 = {x0 / W, y0 / H};
    *uv1 = {(x0 + (float)p->width) / W, (y0 + (float)p->height) / H};
    return true;
}

static void panels_collect(MdRenderPanels* p) {
    if (p->pick_view >= PANEL_MAX || !p->view[p->pick_view]) return;
    md_render_pick_t q;
    while (md_render_view_pick_result(p->view[p->pick_view], &q)) {
        p->pick = q;
        p->has_pick = true;
    }
}

void md_render_bridge_panels_pick(ApplicationState* state, MdRenderPanels* p, uint32_t i, float x, float y) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!b || !b->ok || !p || i >= p->count) return;
    const uint32_t v = p->shared ? 0 : i;
    if (!p->view[v] || x < 0.0f || y < 0.0f || x >= (float)p->width || y >= (float)p->height) return;
    if (v != p->pick_view) {
        // the answers of the view asked about before are of another panel
        panels_collect(p);
        p->has_pick = false;
        p->pick_view = v;
    }
    md_render_view_pick(b->stream, p->view[v], (uint32_t)x, (uint32_t)y);
}

bool md_render_bridge_panels_hit(ApplicationState* state, MdRenderPanels* p, PickingHit* hit) {
    MdRenderBridge* b = state->md_render_bridge;
    if (!b || !b->ok || !p || !hit) return false;
    panels_collect(p);
    if (!p->has_pick) return false;
    const md_render_pick_t& found = p->pick;
    *hit = {};
    hit->frame_idx     = state->picking_handler.frame_idx;
    hit->raw_idx       = INVALID_PICKING_IDX;
    hit->surface_coord = {(float)found.x + 0.5f, (float)p->height - (float)found.y - 0.5f};
    hit->depth         = found.key ? 1.0f - found.depth : 1.0f;
    // The meshes as nothing to pick
    if (found.key && (found.key >> 56) != (PANEL_MESH >> 56)) mol_hit(b, p, state, found, hit);
    return true;
}

void md_render_bridge_panels_surface_hit(PickingHit* hit, ApplicationState* state, MdRenderPanels* panels, uint32_t i,
                                         const InteractionSurfaceState& surface, uint64_t source, const md_float4x4_t& clip_to_world,
                                         int width, int height) {
    const ImVec2 scale = ImGui::GetIO().DisplayFramebufferScale;
    md_render_bridge_panels_pick(state, panels, i, surface.mouse_local.x * scale.x, surface.mouse_local.y * scale.y);
    *hit = {};
    if (md_render_bridge_panels_hit(state, panels, hit)) {
        hit->source = source;
        hit->screen_coord = surface.mouse_local;
        const md_float4_t viewport = {0, 0, (float)width, (float)height};
        hit->world_pos = md_float4x4_unproject({hit->surface_coord.x, hit->surface_coord.y, hit->depth}, clip_to_world, viewport);
    }
}

void md_render_bridge_panels_destroy(ApplicationState* state, MdRenderPanels** pp) {
    if (!pp || !*pp) return;
    MdRenderPanels* p = *pp;
    MdRenderBridge* b = state->md_render_bridge;
    md_allocator_i* alloc = state->allocator.persistent;
    if (b && b->stream) {
        panels_free_views(b, p, 0);
        panels_free_targets(b, p);
        for (size_t i = 0; i < md_array_size(p->colors); ++i) if (p->colors[i].color) md_gpu_free(b->stream, p->colors[i].color);
        for (size_t i = 0; i < md_array_size(p->reps); ++i) if (p->reps[i].visible) md_gpu_free(b->stream, p->reps[i].visible);
        for (size_t i = 0; i < md_array_size(p->profiles); ++i) if (p->profiles[i].addr) md_gpu_free(b->stream, p->profiles[i].addr);
    }
    for (size_t i = 0; i < md_array_size(p->colors); ++i) md_array_free(p->colors[i].shown, alloc);
    md_array_free(p->colors, alloc);
    md_array_free(p->reps, alloc);
    md_array_free(p->profiles, alloc);
    md_array_free(p->batches, alloc);
    delete p;
    *pp = nullptr;
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
    // Whether the targets can be made: once
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
        for (size_t i = 0; i < md_array_size(b->glyphs); ++i) free_glyphs(b, &b->glyphs[i]);
        free_targets(b);
        md_compose_history_destroy(b->stream, b->vis_history);
        b->vis_history = nullptr;
        if (b->view) md_render_view_destroy(b->stream, b->view);
        md_gpu_stream_sync(b->stream);
    }
    if (b->post) md_post_destroy(b->post);
    if (b->compose) md_compose_destroy(b->compose);
    if (b->mol) md_render_mol_destroy(b->mol);
    if (b->render) md_render_destroy(b->render);
    md_array_free(b->reps, state->allocator.persistent);
    md_array_free(b->batches, state->allocator.persistent);
    md_array_free(b->glyphs, state->allocator.persistent);
    md_array_free(b->overlay, state->allocator.persistent);
    md_array_free(b->parts, state->allocator.persistent);
    md_array_free(b->shape_pos_cpu, state->allocator.persistent);
    md_array_free(b->shape_col_cpu, state->allocator.persistent);
    md_array_free(b->shape_scale_cpu, state->allocator.persistent);
    md_array_free(b->dipoles, state->allocator.persistent);
    for (size_t i = 0; i < md_array_size(b->luts); ++i) md_gpu_texture_destroy(b->luts[i].tex);
    md_array_free(b->luts, state->allocator.persistent);
    for (size_t i = 0; i < md_array_size(b->mesh_colors); ++i) {
        if (b->stream && b->mesh_colors[i].color) md_gpu_free(b->stream, b->mesh_colors[i].color);
    }
    md_array_free(b->mesh_colors, state->allocator.persistent);
    delete b;
    state->md_render_bridge = nullptr;
}

const char* md_render_bridge_status(const ApplicationState& state) {
    const MdRenderBridge* b = state.md_render_bridge;
    if (!b) return "Not made";
    return b->status[0] ? b->status : (b->ok ? "Ready" : "Failed");
}

uint32_t md_render_bridge_split_heads(const ApplicationState* state) {
    const MdRenderBridge* b = state->md_render_bridge;
    return b && b->has_topo ? b->topo.lipid_sphere_count[1] : 0;
}

#endif
