#define IMGUI_DEFINE_MATH_OPERATORS

#include <representation.h>

#include <md_util.h>

#include <core/md_log.h>
#include <core/md_bitset.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>

#include <color_utils.h>
#include <color_scale.h>
#include <serialization_utils.h>
#include <imgui_widgets.h>
#include <property_grid.h>

#include <gfx/gl_utils.h>
#include <iso_meshes.h>
#include <script_suggestions.h>
#include <system_scene.h>

#include <imgui.h>
#include <imgui_internal.h>

#include <stdio.h>
#include <string.h>

// ---------------------------------------------------------------------------------------------------------
// The tree, as representation_tree.h sees it
// ---------------------------------------------------------------------------------------------------------

// The nodes of the representations, then those of the groups. Valid until either array changes.
struct TreeNodes {
    RepresentationNode** nodes = nullptr;
    size_t count = 0;
    size_t num_reps = 0;
};

// Reading the tree does not change it, which is all the const_cast is for: the tree functions take the same
// array of nodes whether they move them or not
static TreeNodes tree_nodes(const SystemContext& sys, md_allocator_i* alloc) {
    const size_t num_reps   = md_array_size(sys.representation.reps);
    const size_t num_groups = md_array_size(sys.representation.groups);
    TreeNodes t;
    t.count = num_reps + num_groups;
    t.num_reps = num_reps;
    t.nodes = (RepresentationNode**)md_alloc(alloc, sizeof(RepresentationNode*) * MAX(t.count, (size_t)1));
    for (size_t i = 0; i < num_reps; ++i) {
        t.nodes[i] = const_cast<RepresentationNode*>(&sys.representation.reps[i].node);
    }
    for (size_t i = 0; i < num_groups; ++i) {
        t.nodes[num_reps + i] = const_cast<RepresentationNode*>(&sys.representation.groups[i].node);
    }
    return t;
}

// Kept in the order of the window, so that whatever walks the arrays - drawing, the legends, saving - meets
// them in the order the user put them in
template <typename T>
static void sort_by_order(T* items, size_t count) {
    for (size_t i = 1; i < count; ++i) {
        if (items[i - 1].node.order <= items[i].node.order) continue;
        T item = items[i];
        size_t j = i;
        while (j > 0 && items[j - 1].node.order > item.node.order) {
            items[j] = items[j - 1];
            j -= 1;
        }
        items[j] = item;
    }
}

// Pointers into the arrays do not survive this: the items are moved to their places
static void tree_normalize(SystemContext& sys) {
    md_temp_scope_t temp = md_temp_begin();
    TreeNodes t = tree_nodes(sys, md_temp_allocator(temp));
    representation_tree_normalize(t.nodes, t.count);
    md_temp_end(temp);

    sort_by_order(sys.representation.reps,   md_array_size(sys.representation.reps));
    sort_by_order(sys.representation.groups, md_array_size(sys.representation.groups));
}

static int find_rep_idx(const SystemContext& sys, uint32_t id) {
    for (size_t i = 0; i < md_array_size(sys.representation.reps); ++i) {
        if (sys.representation.reps[i].node.id == id) return (int)i;
    }
    return -1;
}

static int find_group_idx(const SystemContext& sys, uint32_t id) {
    if (id == 0) return -1;
    for (size_t i = 0; i < md_array_size(sys.representation.groups); ++i) {
        if (sys.representation.groups[i].node.id == id) return (int)i;
    }
    return -1;
}

Representation* find_representation(SystemContext& sys, uint32_t id) {
    const int i = find_rep_idx(sys, id);
    return i >= 0 ? &sys.representation.reps[i] : nullptr;
}

RepresentationGroup* find_representation_group(SystemContext& sys, uint32_t id) {
    const int i = find_group_idx(sys, id);
    return i >= 0 ? &sys.representation.groups[i] : nullptr;
}

const RepresentationGroup* find_representation_group(const SystemContext& sys, uint32_t id) {
    const int i = find_group_idx(sys, id);
    return i >= 0 ? &sys.representation.groups[i] : nullptr;
}

static bool group_has_filter(const RepresentationGroup& g) {
    return !str_empty(str_trim(str_from_cstr(g.filter.text)));
}

// Whether every group from 'parent' outwards is enabled. A parent that names no group is the top level,
// as it is to the tree.
static bool groups_enabled(const SystemContext& sys, uint32_t parent) {
    const size_t num_groups = md_array_size(sys.representation.groups);
    for (size_t steps = 0; parent != 0 && steps <= num_groups; ++steps) {
        const RepresentationGroup* g = find_representation_group(sys, parent);
        if (!g) break;
        if (!g->enabled) return false;
        parent = g->node.parent;
    }
    return true;
}

bool representation_shown(const SystemContext& sys, const Representation& rep) {
    return sys.scene.visible && rep.enabled && groups_enabled(sys, rep.node.parent);
}

bool representation_group_shown(const SystemContext& sys, const RepresentationGroup& group) {
    return sys.scene.visible && group.enabled && groups_enabled(sys, group.node.parent);
}

RepresentationContext representation_context(const SystemContext& sys, uint32_t parent) {
    RepresentationContext ctx;
    const size_t num_groups = md_array_size(sys.representation.groups);
    for (size_t steps = 0; parent != 0 && steps <= num_groups; ++steps) {
        const RepresentationGroup* g = find_representation_group(sys, parent);
        if (!g) break;
        if (group_has_filter(*g)) {
            // The innermost one is all it takes: its own filter was evaluated within the context around it,
            // and it carries that context's validity and dependence on the frame along
            ctx.group   = g;
            ctx.atoms   = &g->atom_mask;
            ctx.valid   = g->filter.result.success;
            ctx.dynamic = g->filter.result.is_dynamic;
            return ctx;
        }
        parent = g->node.parent;
    }
    return ctx;
}

bool representation_is_dynamic(const SystemContext& sys, const Representation& rep) {
    // The result has the context's in it as of the last evaluation; the context is asked as well for one that has
    // come to depend on the frame since, which flags what is within it for evaluation anyway
    return rep.filter.result.is_dynamic || representation_context(sys, rep.node.parent).dynamic;
}

script::FilterCtx representation_filter_context(const SystemContext& sys, const RepresentationContext& ctx) {
    script::FilterCtx fctx = filter_context(sys);
    fctx.atoms = ctx.atoms;
    return fctx;
}

// Everything within a group is evaluated in a context that just changed, or may have: the groups' contexts
// and the representations' atoms. For a representation it is just itself.
static void flag_within_as_dirty(SystemContext& sys, bool group, uint32_t id) {
    if (!group) {
        if (Representation* rep = find_representation(sys, id)) flag_representation_as_dirty(rep);
        sys.representation.atom_visibility_mask_dirty = true;
        return;
    }

    md_temp_scope_t temp = md_temp_begin();
    TreeNodes t = tree_nodes(sys, md_temp_allocator(temp));
    for (size_t k = 0; k < t.count; ++k) {
        const RepresentationNode* node = t.nodes[k];
        const bool self = node->group && node->id == id;
        if (!self && !representation_tree_is_within(t.nodes, t.count, k, id)) continue;
        if (k < t.num_reps) {
            flag_representation_as_dirty(&sys.representation.reps[k]);
        } else {
            sys.representation.groups[k - t.num_reps].filt_is_dirty = true;
        }
    }
    md_temp_end(temp);
    sys.representation.atom_visibility_mask_dirty = true;
}

// The ids of what is within the group, representations and groups apart, in window order
static void collect_within(md_array(uint32_t)* rep_ids, md_array(uint32_t)* group_ids, const SystemContext& sys, uint32_t id, md_allocator_i* alloc) {
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    TreeNodes t = tree_nodes(sys, md_temp_allocator(temp));
    for (size_t k = 0; k < t.count; ++k) {
        if (!representation_tree_is_within(t.nodes, t.count, k, id)) continue;
        if (k < t.num_reps) {
            if (rep_ids) md_array_push(*rep_ids, t.nodes[k]->id, alloc);
        } else {
            if (group_ids) md_array_push(*group_ids, t.nodes[k]->id, alloc);
        }
    }
    md_temp_end(temp);
}

bool representation_can_move(const SystemContext& sys, bool group, uint32_t id, bool target_group, uint32_t target_id, RepresentationDrop drop) {
    md_temp_scope_t temp = md_temp_begin();
    TreeNodes t = tree_nodes(sys, md_temp_allocator(temp));
    const int idx    = representation_tree_find(t.nodes, t.count, group, id);
    const int target = target_id ? representation_tree_find(t.nodes, t.count, target_group, target_id) : -1;
    const bool ok = idx >= 0 && (target_id == 0 || target >= 0) && representation_tree_can_move(t.nodes, t.count, (size_t)idx, target, drop);
    md_temp_end(temp);
    return ok;
}

bool representation_move(SystemContext& sys, bool group, uint32_t id, bool target_group, uint32_t target_id, RepresentationDrop drop) {
    md_temp_scope_t temp = md_temp_begin();
    TreeNodes t = tree_nodes(sys, md_temp_allocator(temp));
    const int idx    = representation_tree_find(t.nodes, t.count, group, id);
    const int target = target_id ? representation_tree_find(t.nodes, t.count, target_group, target_id) : -1;
    bool moved = false;
    uint32_t old_parent = 0, new_parent = 0;
    if (idx >= 0 && (target_id == 0 || target >= 0)) {
        old_parent = t.nodes[idx]->parent;
        moved = representation_tree_move(t.nodes, t.count, (size_t)idx, target, drop);
        new_parent = t.nodes[idx]->parent;
    }
    md_temp_end(temp);

    if (!moved) return false;

    sort_by_order(sys.representation.reps,   md_array_size(sys.representation.reps));
    sort_by_order(sys.representation.groups, md_array_size(sys.representation.groups));

    // Into another group is into another context, and maybe out of sight. Among its siblings it is neither.
    if (old_parent != new_parent) {
        flag_within_as_dirty(sys, group, id);
    }
    return true;
}

// ---------------------------------------------------------------------------------------------------------
// Representations
// ---------------------------------------------------------------------------------------------------------

static void init_representation(SystemContext& sys, Representation* rep) {
    ApplicationState* state = sys.app;
#if EXPERIMENTAL_GFX_API
    rep->gfx_rep = md_gfx_rep_create(sys.data.sys.atom.count);
#endif
    rep->md_rep = md_gl_rep_create(sys.data.gl_mol);
    // Of the glyphs of the system loaded now: one of a system before it goes
    md_gl_rep_destroy(rep->glyph_rep);
    rep->glyph_rep = cartoon_glyphs_rep_create(sys.data.cartoon_glyphs);
    rep->atom_mask = md_bitset_empty_set();

    // Default to the first per atom field the system offers, if it offers any - unless the
    // representation already names one the system has. This runs again for every representation
    // when a system is loaded, which is AFTER a workspace's representations were read, and for a
    // clone: neither may lose the attribute it was coloured by, or how. There is no list to consult:
    // the attribute table is the list.
    if (!md_attributes_get(&sys.data.sys.attributes, rep->atom_attribute.key)) {
        md_attribute_id_t first_attribute = MD_ATTRIBUTE_INVALID;
        if (atom_attribute_query(&first_attribute, 1, sys.data.sys) > 0) {
            atom_attribute_select(&rep->atom_attribute, first_attribute, sys.data.sys);
        }
    }

    // A system loaded under a representation that colours by a field: the field's values were of the
    // old one, whatever the frame says
    surface_field_invalidate(&rep->electronic_structure.field_vol);
    // and its meshes of the old one's basis
    electronic_structure_mesh_free(state, &rep->electronic_structure);

    flag_representation_as_dirty(rep);
}

// What a copy must not share with its original: the GPU side of it and the evaluated atoms. The clone
// evaluates its own, and must neither write into nor free the original's (each frees its own on removal).
static void clear_owned(Representation* rep) {
    rep->md_rep = {};
    rep->glyph_rep = {};
    rep->atom_mask = {};
    rep->bond_pairs = nullptr;
    rep->electronic_structure.density_vol.tex_id = 0;
    rep->electronic_structure.color_vol.tex_id   = 0;
    rep->electronic_structure.field_vol          = SurfaceFieldVolume{};
    rep->electronic_structure.iso_mesh           = nullptr;
}

// Appends a copy, at the place in the tree it is given, and does not normalize: that is the caller's to do
// once it has placed all it is copying. Returns its id.
static uint32_t push_copy(SystemContext& sys, const Representation& src, uint32_t parent, int32_t order) {
    ApplicationState* state = sys.app;
    // A copy first: src may well be an element of the array the push is about to grow
    Representation copy = src;
    clear_owned(&copy);
    copy.node.id     = sys.representation.next_rep_id++;
    copy.node.parent = parent;
    copy.node.order  = order;
    copy.node.group  = false;
    md_array_push(sys.representation.reps, copy, state->allocator.persistent);
    init_representation(sys, md_array_last(sys.representation.reps));
    return copy.node.id;
}

bool representation_can_add(const SystemContext& sys, size_t count) {
    return md_array_size(sys.representation.reps) + count <= REPRESENTATION_MAX;
}

// Said when one is not made for the limit
static void log_representation_limit() {
    VIAMD_LOG_ERROR("A system has at most %zu representations; remove one to make another", REPRESENTATION_MAX);
}

Representation* create_representation(SystemContext& sys, RepresentationType type, ColorMapping color_mapping, str_t filter, uint32_t group) {
    ApplicationState* state = sys.app;
    ASSERT(state);
    if (!representation_can_add(sys)) {
        log_representation_limit();
        return nullptr;
    }
    md_array_push(sys.representation.reps, Representation(), state->allocator.persistent);
    Representation* rep = md_array_last(sys.representation.reps);
    rep->node.id     = sys.representation.next_rep_id++;
    rep->node.parent = group;
    rep->node.order  = INT32_MAX;   // the last in its group
    rep->type = type;
    rep->color_mapping = color_mapping;
    if (!str_empty(filter)) {
        str_copy_to_char_buf(rep->filter.text, sizeof(rep->filter.text), filter);
    }
    // Opens on the HOMO, derived from the occupations rather than handed over by whoever loaded
    // the file. -1 (nothing occupied) clamps to 0, which is the only orbital there is to show.
    OrbitalFrontier frontier = {};
    es_orbital_frontier(&frontier, sys.data.sys, es_path::alpha_occupation);
    rep->electronic_structure.orbital_idx = MAX(frontier.homo_idx, 0);
    init_representation(sys, rep);

    const uint32_t id = rep->node.id;
    tree_normalize(sys);
    return find_representation(sys, id);
}

Representation* clone_representation(SystemContext& sys, const Representation& rep) {
    if (!representation_can_add(sys)) {
        log_representation_limit();
        return nullptr;
    }
    // Taken before anything moves: rep may be an element of the array
    const Representation src = rep;
    tree_normalize(sys);
    // Just after the original, which is the original's order + 1 once normalized. One that is not in the
    // array goes last.
    const Representation* orig = find_representation(sys, src.node.id);
    const uint32_t id = orig ? push_copy(sys, src, orig->node.parent, orig->node.order + 1)
                             : push_copy(sys, src, src.node.parent, INT32_MAX);
    tree_normalize(sys);
    sys.representation.atom_visibility_mask_dirty = true;
    return find_representation(sys, id);
}

void remove_representation(SystemContext& sys, size_t idx) {
    ApplicationState* state = sys.app;
    ASSERT(state);
    ASSERT(idx < md_array_size(sys.representation.reps));
    auto& rep = sys.representation.reps[idx];
    md_bitset_free(&rep.atom_mask, state->allocator.persistent);
    md_array_free(rep.bond_pairs, state->allocator.persistent);
    md_gl_rep_destroy(rep.md_rep);
    md_gl_rep_destroy(rep.glyph_rep);
    // A readback queued for this representation's volume would land in a texture that no longer
    // exists, so let the queue run out first.
    gpu_volume_jobs_drain(state);
    if (rep.electronic_structure.density_vol.tex_id) gl::free_texture(&rep.electronic_structure.density_vol.tex_id);
    if (rep.electronic_structure.color_vol.tex_id)   gl::free_texture(&rep.electronic_structure.color_vol.tex_id);
    surface_field_free(&rep.electronic_structure.field_vol);
    electronic_structure_mesh_free(state, &rep.electronic_structure);

    // In order: the array is the order of the window
    const size_t count = md_array_size(sys.representation.reps);
    for (size_t i = idx; i + 1 < count; ++i) {
        sys.representation.reps[i] = sys.representation.reps[i + 1];
    }
    md_array_pop(sys.representation.reps);
    recompute_atom_visibility_mask(sys);
}

static void free_group(SystemContext& sys, size_t idx) {
    ApplicationState* state = sys.app;
    ASSERT(idx < md_array_size(sys.representation.groups));
    md_bitset_free(&sys.representation.groups[idx].atom_mask, state->allocator.persistent);
    const size_t count = md_array_size(sys.representation.groups);
    for (size_t i = idx; i + 1 < count; ++i) {
        sys.representation.groups[i] = sys.representation.groups[i + 1];
    }
    md_array_pop(sys.representation.groups);
}

void remove_all_representations(SystemContext& sys) {
    while (md_array_size(sys.representation.reps) > 0) {
        remove_representation(sys, md_array_size(sys.representation.reps) - 1);
    }
    while (md_array_size(sys.representation.groups) > 0) {
        free_group(sys, md_array_size(sys.representation.groups) - 1);
    }
}

void recompute_atom_visibility_mask(SystemContext& sys) {
    auto& mask = sys.representation.visibility_mask;

    md_bitmask_clear(&mask);
    for (size_t i = 0; i < md_array_size(sys.representation.reps); ++i) {
        const Representation& rep = sys.representation.reps[i];
        if (!representation_shown(sys, rep)) continue;
        md_bitmask_or(&mask, rep.atom_mask);
    }
    sys.representation.visibility_mask_hash = md_bitset_hash64(md_bitmask_view(&mask), 0);
}

// ---------------------------------------------------------------------------------------------------------
// Groups
// ---------------------------------------------------------------------------------------------------------

RepresentationGroup* create_representation_group(SystemContext& sys, str_t name, uint32_t parent) {
    ApplicationState* state = sys.app;
    ASSERT(state);
    RepresentationGroup group = {};
    group.node.id     = sys.representation.next_group_id++;
    group.node.parent = parent;
    group.node.order  = INT32_MAX;
    group.node.group  = true;
    str_copy_to_char_buf(group.name, sizeof(group.name), str_empty(name) ? STR_LIT("group") : name);
    md_array_push(sys.representation.groups, group, state->allocator.persistent);
    tree_normalize(sys);
    return find_representation_group(sys, group.node.id);
}

size_t representation_group_rep_count(const SystemContext& sys, uint32_t id) {
    md_temp_scope_t temp = md_temp_begin();
    md_array(uint32_t) rep_ids = 0;
    md_array(uint32_t) group_ids = 0;
    collect_within(&rep_ids, &group_ids, sys, id, md_temp_allocator(temp));
    const size_t n = md_array_size(rep_ids);
    md_temp_end(temp);
    return n;
}

void clone_representation_group(SystemContext& sys, uint32_t id) {
    ApplicationState* state = sys.app;
    ASSERT(state);
    tree_normalize(sys);
    const RepresentationGroup* orig = find_representation_group(sys, id);
    if (!orig) return;

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    md_array(uint32_t) rep_ids = 0;
    md_array(uint32_t) group_ids = 0;
    collect_within(&rep_ids, &group_ids, sys, id, alloc);

    // All of it or none: a copy without some of what it held would not be a copy
    if (!representation_can_add(sys, md_array_size(rep_ids))) {
        log_representation_limit();
        md_temp_end(temp);
        return;
    }

    // The groups first, so that every parent within the copy has its new id before anything refers to it.
    // The copy of the group itself goes just after the original and what it holds; what is within it keeps
    // the orders it had, which are only ever compared with each other.
    const size_t num_groups = md_array_size(group_ids) + 1;
    uint32_t* old_ids = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * num_groups);
    uint32_t* new_ids = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * num_groups);
    old_ids[0] = id;
    for (size_t i = 1; i < num_groups; ++i) old_ids[i] = group_ids[i - 1];

    auto map_parent = [&](uint32_t parent) -> uint32_t {
        for (size_t i = 0; i < num_groups; ++i) {
            if (old_ids[i] == parent) return new_ids[i];
        }
        return parent;
    };

    for (size_t i = 0; i < num_groups; ++i) {
        new_ids[i] = sys.representation.next_group_id++;
    }
    for (size_t i = 0; i < num_groups; ++i) {
        const RepresentationGroup* src = find_representation_group(sys, old_ids[i]);
        if (!src) continue;
        RepresentationGroup copy = *src;
        copy.atom_mask = {};            // its own, evaluated again
        copy.filt_is_dirty = true;
        copy.node.id = new_ids[i];
        if (i == 0) {
            copy.node.order = src->node.order + 1;
        } else {
            copy.node.parent = map_parent(src->node.parent);
        }
        md_array_push(sys.representation.groups, copy, state->allocator.persistent);
    }
    for (size_t i = 0; i < md_array_size(rep_ids); ++i) {
        const Representation* src = find_representation(sys, rep_ids[i]);
        if (!src) continue;
        push_copy(sys, *src, map_parent(src->node.parent), src->node.order);
    }

    md_temp_end(temp);
    tree_normalize(sys);
    sys.representation.atom_visibility_mask_dirty = true;
}

void remove_representation_group(SystemContext& sys, uint32_t id) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    md_array(uint32_t) rep_ids = 0;
    md_array(uint32_t) group_ids = 0;
    collect_within(&rep_ids, &group_ids, sys, id, alloc);
    md_array_push(group_ids, id, alloc);

    for (size_t i = 0; i < md_array_size(rep_ids); ++i) {
        const int idx = find_rep_idx(sys, rep_ids[i]);
        if (idx >= 0) remove_representation(sys, (size_t)idx);
    }
    for (size_t i = 0; i < md_array_size(group_ids); ++i) {
        const int idx = find_group_idx(sys, group_ids[i]);
        if (idx >= 0) free_group(sys, (size_t)idx);
    }

    md_temp_end(temp);
    tree_normalize(sys);
    recompute_atom_visibility_mask(sys);
}

void ungroup_representation_group(SystemContext& sys, uint32_t id) {
    // Normalized, everything within the group is numbered between the group and its next sibling: given the
    // group's parent, what it holds sorts exactly into its place
    tree_normalize(sys);
    const RepresentationGroup* g = find_representation_group(sys, id);
    if (!g) return;
    const uint32_t parent = g->node.parent;

    for (size_t i = 0; i < md_array_size(sys.representation.reps); ++i) {
        Representation& rep = sys.representation.reps[i];
        if (rep.node.parent != id) continue;
        rep.node.parent = parent;
        flag_within_as_dirty(sys, false, rep.node.id);
    }
    for (size_t i = 0; i < md_array_size(sys.representation.groups); ++i) {
        RepresentationGroup& child = sys.representation.groups[i];
        if (child.node.parent != id) continue;
        child.node.parent = parent;
        flag_within_as_dirty(sys, true, child.node.id);
    }

    const int idx = find_group_idx(sys, id);
    if (idx >= 0) free_group(sys, (size_t)idx);
    tree_normalize(sys);
    sys.representation.atom_visibility_mask_dirty = true;
}

// A new group in the place of a representation or a group, which goes into it
static void wrap_in_group(SystemContext& sys, bool group, uint32_t id) {
    ApplicationState* state = sys.app;
    tree_normalize(sys);
    const RepresentationNode* node = nullptr;
    if (group) {
        if (const RepresentationGroup* g = find_representation_group(sys, id)) node = &g->node;
    } else {
        if (const Representation* r = find_representation(sys, id)) node = &r->node;
    }
    if (!node) return;

    RepresentationGroup wrapper = {};
    wrapper.node.id     = sys.representation.next_group_id++;
    wrapper.node.parent = node->parent;
    wrapper.node.order  = node->order - 1;
    wrapper.node.group  = true;
    snprintf(wrapper.name, sizeof(wrapper.name), "group");
    md_array_push(sys.representation.groups, wrapper, state->allocator.persistent);

    // The push may have moved the groups: the node is looked up again
    if (group) {
        if (RepresentationGroup* g = find_representation_group(sys, id)) g->node.parent = wrapper.node.id;
    } else {
        if (Representation* r = find_representation(sys, id)) r->node.parent = wrapper.node.id;
    }
    // A group without a filter changes no context
    tree_normalize(sys);
}

// ---------------------------------------------------------------------------------------------------------
// Evaluation
// ---------------------------------------------------------------------------------------------------------

// The context of a group: its own filter evaluated within the context it is in, which comes first
static void update_group(SystemContext& sys, size_t idx, size_t depth) {
    ApplicationState* state = sys.app;
    RepresentationGroup* g = &sys.representation.groups[idx];
    if (!g->filt_is_dirty) return;
    if (depth > md_array_size(sys.representation.groups)) return;    // A cycle; the tree breaks those

    const int parent_idx = find_group_idx(sys, g->node.parent);
    if (parent_idx >= 0) {
        update_group(sys, (size_t)parent_idx, depth + 1);
    }

    if (!group_has_filter(*g)) {
        md_bitset_free(&g->atom_mask, state->allocator.persistent);
        g->filter.result = {};
        g->filter.result.success = true;
    } else {
        // Without a system there is nothing to evaluate against; it is done once there is one
        if (sys.data.sys.atom.count == 0) return;

        const RepresentationContext ctx = representation_context(sys, g->node.parent);
        if (!ctx.valid) {
            md_bitset_free(&g->atom_mask, state->allocator.persistent);
            g->filter.result = {};
            // A context that follows the frame can come back valid in another frame, and this with it
            g->filter.result.is_dynamic = ctx.dynamic;
            snprintf(g->filter.result.error, sizeof(g->filter.result.error), "The filter of the group '%s' it is in is not valid", ctx.group ? ctx.group->name : "");
        } else {
            filter_evaluate(&g->filter, &g->atom_mask, state->script.language, representation_filter_context(sys, ctx), state->allocator.persistent);
            g->filter.result.is_dynamic |= ctx.dynamic;
        }
    }
    g->filt_is_dirty = false;
}

bool representation_uses_atom_colors(const Representation& rep) {
    switch (rep.type) {
        case RepresentationType::ElectronicStructure:
            return rep.electronic_structure.coloring == SurfaceColoring::AtomColors;
        case RepresentationType::DipoleMoment:
        case RepresentationType::SimulationBox:
            return false;
        default:
            return true;
    }
}

bool representation_particle_type_valid(const SystemContext& sc, RepresentationType type) {
    const md_system_t& sys = sc.data.sys;
    switch (type) {
    case RepresentationType::SpaceFill:     return sys.atom.count > 0;
    case RepresentationType::Licorice:      return sys.bond.count > 0;
    case RepresentationType::BallAndStick:  return sys.atom.count > 0;
    case RepresentationType::Ribbons:       return sys.protein_backbone.range.count > 0;
    // Proteins, and the lipids and nucleic acids of the glyphs
    case RepresentationType::Cartoon:       return sys.protein_backbone.range.count > 0 || cartoon_glyphs_valid(sc.data.cartoon_glyphs);
    default:                                return false;
    }
}

md_gl_draw_op_t representation_draw_op(const Representation& rep) {
    md_gl_draw_op_t op = {
        .type = (md_gl_rep_type_t)rep.type,
        .args = {},
        .rep = rep.md_rep,
    };
    switch (rep.type) {
    case RepresentationType::SpaceFill:
        op.args.space_fill.radius_scale = rep.scale.x;
        break;
    case RepresentationType::Licorice:
        op.args.licorice.radius = rep.scale.x;
        op.args.licorice.color_mode = (md_gl_bond_mode_t)rep.bond_color;
        op.args.licorice.sharpness = rep.bond_sharpness;
        op.args.licorice.uniform_color = convert_color(rep.bond_base_color);
        op.args.licorice.bond_orders = rep.bond_orders;
        op.args.licorice.chain_gaps = rep.chain_gaps;
        op.args.licorice.elastic_bonds = rep.elastic_bonds;
        if (rep.tint_scale > 0.0f || rep.saturation < 1.0f) {
            tint_colors(&op.args.licorice.uniform_color, 1, convert_color(rep.tint_color), rep.tint_scale, rep.saturation);
        }
        break;
    case RepresentationType::BallAndStick:
        op.args.ball_and_stick.ball_scale = rep.scale.x;
        op.args.ball_and_stick.stick_radius = rep.scale.y;
        op.args.ball_and_stick.color_mode = (md_gl_bond_mode_t)rep.bond_color;
        op.args.ball_and_stick.sharpness = rep.bond_sharpness;
        op.args.ball_and_stick.uniform_color = convert_color(rep.bond_base_color);
        op.args.ball_and_stick.bond_orders = rep.bond_orders;
        op.args.ball_and_stick.chain_gaps = rep.chain_gaps;
        op.args.ball_and_stick.elastic_bonds = rep.elastic_bonds;
        if (rep.tint_scale > 0.0f || rep.saturation < 1.0f) {
            tint_colors(&op.args.ball_and_stick.uniform_color, 1, convert_color(rep.tint_color), rep.tint_scale, rep.saturation);
        }
        break;
    case RepresentationType::Ribbons:
        op.args.ribbons.width_scale = rep.scale.x;
        op.args.ribbons.thickness_scale = rep.scale.y;
        op.args.ribbons.chain_gaps = rep.chain_gaps;
        break;
    case RepresentationType::Cartoon:
        op.args.cartoon.coil_scale = rep.scale.x;
        op.args.cartoon.sheet_scale = rep.scale.y;
        op.args.cartoon.helix_scale = rep.scale.z;
        op.args.cartoon.chain_gaps = rep.chain_gaps;
        break;
    default:
        ASSERT(false);
        break;
    }
    return op;
}

// The draws of the glyphs of a Cartoon (see representation.h)
void representation_glyph_draw_ops(md_array(md_gl_draw_op_t)* ops, const SystemContext& sys, const Representation& rep, md_allocator_i* alloc) {
    if (!cartoon_glyphs_valid(sys.data.cartoon_glyphs) || !rep.glyph_rep.id || !rep.glyphs_shown || !sys.data.cartoon_glyphs.positions_valid) return;
    const float scale = rep.scale.w > 0.0f ? rep.scale.w : 1.0f;

    // The heads; the other points (control points, corners of slabs) have no sphere, and nucleic acids none at all
    if (sys.data.cartoon_glyphs.num_spheres > 0) {
        md_gl_draw_op_t spheres = {
            .type = MD_GL_REP_SPACE_FILL,
            .args = {},
            .rep = rep.glyph_rep,
            .copies = {},
        };
        spheres.args.space_fill.radius_scale = scale;
        md_array_push(*ops, spheres, alloc);
    }

    if (sys.data.cartoon_glyphs.num_tubes > 0) {
        md_gl_draw_op_t tubes = {
            .type = MD_GL_REP_TUBES,
            .args = {},
            .rep = rep.glyph_rep,
            .copies = {},
        };
        tubes.args.tubes.radius_scale = scale;
        // The bands of the nucleic acids are made with the default radius, a representation has its own
        tubes.args.tubes.band_scale = MAX(rep.nucleic_band_radius, 0.01f) / CARTOON_NUCLEIC_BAND_RADIUS;
        md_array_push(*ops, tubes, alloc);
    }

    if (sys.data.cartoon_glyphs.num_slabs > 0) {
        md_gl_draw_op_t slabs = {
            .type = MD_GL_REP_SLABS,
            .args = {},
            .rep = rep.glyph_rep,
            .copies = {},
        };
        // The rims (and the sticks, slabs of two corners) of the radius of the representation, the plate within a rim
        // recessed to 60% of its thickness
        slabs.args.slabs.radius = MAX(rep.nucleic_ring_radius, 0.01f) * scale;
        slabs.args.slabs.plate  = 0.6f;
        md_array_push(*ops, slabs, alloc);
    }
}

// The colours of all atoms of the system for a representation which colours by atom, by its colour mapping, then towards
// its tint and saturation (see representation.h)
void representation_atom_colors(uint32_t* colors, size_t num_atoms, SystemContext& sc, Representation* rep, const md_bitset_t* shown) {
    ApplicationState* state = sc.app;
    const md_system_t& sys = sc.data.sys;
    md_allocator_i* frame_alloc = state->allocator.frame;
    const size_t bytes = num_atoms * sizeof(uint32_t);

    switch (rep->color_mapping) {
    case ColorMapping::Uniform:
        color_atoms_uniform(colors, num_atoms, convert_color(rep->base_color));
        break;
    case ColorMapping::Type:
        color_atoms_type(colors, num_atoms, sys);
        break;
    case ColorMapping::Serial:
        color_atoms_idx(colors, num_atoms, sys);
        break;
    case ColorMapping::ResName:
        color_atoms_res_name(colors, num_atoms, sys);
        break;
    case ColorMapping::ResType:
        color_atoms_res_type(colors, num_atoms, sys);
        break;
    case ColorMapping::ResSeqId:
        color_atoms_res_seq_id(colors, num_atoms, sys);
        break;
    case ColorMapping::ResIndex:
        color_atoms_res_idx(colors, num_atoms, sys);
        break;
    case ColorMapping::ChainId:
        color_atoms_chain_id(colors, num_atoms, sys);
        break;
    case ColorMapping::ChainIndex:
        color_atoms_chain_idx(colors, num_atoms, sys);
        break;
    case ColorMapping::PartIndex:
        color_atoms_part_idx(colors, num_atoms, sys);
        break;
    case ColorMapping::SecondaryStructure: {
        SecondaryStructurePalette palette = {
            .coil  = convert_color(rep->secondary_structure.color_coil),
            .helix = convert_color(rep->secondary_structure.color_helix),
            .sheet = convert_color(rep->secondary_structure.color_sheet),
        };
        color_atoms_secondary_structure(colors, num_atoms, sys, displayed_secondary_structure(sc), palette);
        break;
    }
    case ColorMapping::Attribute: {
        AtomAttributeColoring& coloring = rep->atom_attribute;
        const md_attribute_t* attr = md_attributes_get(&sys.attributes, coloring.key);
        size_t num_extracted = 0;
        float* values = nullptr;

        if (attr) {
            values = (float*)md_vm_arena_push(frame_alloc, sizeof(float) * num_atoms);

            // Fixing the variant axis hands back exactly the atom axis, so there is no
            // offset arithmetic here to get wrong. A field with no variant axis is rank 1
            // and takes no indices at all.
            const uint32_t variant = (uint32_t)CLAMP(coloring.variant_idx, 0, MAX(atom_attribute_variant_count(attr) - 1, 0));
            const md_attribute_slice_t slice = attr->format.rank > 1 ? md_attribute_slice_1(variant) : md_attribute_slice_all();
            num_extracted = md_attribute_extract_f32(values, num_atoms, attr, slice, md_unit_none());
        }

        if (num_extracted == num_atoms) {
            // The span of the atoms shown, measured again only when the values or the atoms
            // changed: it scans every variant of the attribute
            uint64_t span_hash = md_hash64(&coloring.key, sizeof(coloring.key), 1);
            const uint64_t version = md_attributes_version(&sys.attributes, coloring.key);
            span_hash = md_hash64(&version, sizeof(version), span_hash);
            span_hash = shown ? md_bitset_hash64(*shown, span_hash) : md_hash64(&num_atoms, sizeof(num_atoms), span_hash);
            if (span_hash != coloring.span_hash) {
                coloring.span = atom_attribute_span(attr, shown);
                coloring.span_hash = span_hash;
            }
            color_scale_update_range(&coloring.scale, coloring.span);

            for (size_t i = 0; i < num_atoms; ++i) {
                // An atom without a value is not a point on the ramp: NAN would reach the
                // colormap lookup as an index, so it gets a neutral grey of its own instead
                colors[i] = atom_attribute_value_absent(values[i]) ? IM_COL32(128, 128, 128, 255) : color_scale_color_u32(coloring.scale, values[i]);
            }
        } else {
            if (attr) {
                MD_LOG_DEBUG("Failed to extract values for the selected atom attribute");
            }
            MEMSET(colors, 0xFFFFFFFFu, bytes);
        }
        break;
    }
    default:
        ASSERT(false);
        break;
    }

    if (rep->tint_scale > 0.0f || rep->saturation < 1.0f) {
        uint32_t tint_color = convert_color(rep->tint_color);
        tint_colors(colors, num_atoms, tint_color, rep->tint_scale, rep->saturation);
    }
}

static void update_representation(SystemContext& sc, Representation* rep) {
    ApplicationState* state = sc.app;
    ASSERT(state);
    ASSERT(rep);

    // A hidden one keeps what needs doing for when it is shown again
    if (!rep->needs_update) return;
    if (!representation_shown(sc, *rep)) return;

    const auto& sys = sc.data.sys;
    size_t num_atoms = md_system_atom_count(&sys);

    md_allocator_i* frame_alloc = state->allocator.frame;
    md_temp_scope_t temp = md_temp_begin_in(frame_alloc);
    defer { md_temp_end(temp); };

    bool again = false;     // not all done: updated again next frame
    uint32_t* colors = 0;
    if (representation_uses_atom_colors(*rep)) {
        colors = (uint32_t*)md_vm_arena_push(frame_alloc, sizeof(uint32_t) * num_atoms);

        // The atoms first: the colours are applied through them, and a colour scale's range can
        // follow the values of exactly the atoms shown
        if (rep->dynamic_evaluation && representation_is_dynamic(sc, *rep)) {
            rep->filt_is_dirty = true;
        }
        if (rep->filt_is_dirty) {
            // Within the context of the groups it is in. A context that is not valid leaves nothing to
            // evaluate within, which is what an invalid filter of its own leaves too.
            const RepresentationContext ctx = representation_context(sc, rep->node.parent);
            if (ctx.valid && rep->type == RepresentationType::Bonds) {
                // The bonds, and the atoms at their ends: those are what is coloured, and shown
                filter_evaluate_bonds(&rep->filter, &rep->bond_pairs, state->script.language, representation_filter_context(sc, ctx), state->allocator.persistent);
                md_bitmask_t ends = md_bitmask_create(0, (uint32_t)num_atoms, frame_alloc);
                for (size_t i = 0; i < md_array_size(rep->bond_pairs); ++i) {
                    md_bitmask_set_bit(&ends, (uint32_t)rep->bond_pairs[i].idx[0]);
                    md_bitmask_set_bit(&ends, (uint32_t)rep->bond_pairs[i].idx[1]);
                }
                md_bitset_free(&rep->atom_mask, state->allocator.persistent);
                rep->atom_mask = md_bitmask_extract(&ends, state->allocator.persistent);
                rep->filter.result.is_dynamic |= ctx.dynamic;
            } else if (ctx.valid) {
                filter_evaluate(&rep->filter, &rep->atom_mask, state->script.language, representation_filter_context(sc, ctx), state->allocator.persistent);
                // What it picks follows the frame if the context it picks within does
                rep->filter.result.is_dynamic |= ctx.dynamic;
            } else {
                md_bitset_free(&rep->atom_mask, state->allocator.persistent);
                md_array_shrink(rep->bond_pairs, 0);
                rep->filter.result = {};
                rep->filter.result.is_dynamic = ctx.dynamic;     // it may come back valid in another frame
                snprintf(rep->filter.result.error, sizeof(rep->filter.result.error), "The filter of the group '%s', the context of this one, is not valid", ctx.group ? ctx.group->name : "");
            }
            rep->filt_is_dirty = false;
            sc.representation.atom_visibility_mask_dirty = true;
        }

        // Of the atoms shown, for a colour scale which spans their values
        representation_atom_colors(colors, num_atoms, sc, rep, rep->filter.result.success ? &rep->atom_mask : nullptr);
    }

    switch (rep->type) {
    case RepresentationType::SpaceFill:
    case RepresentationType::Licorice:
    case RepresentationType::BallAndStick:
    case RepresentationType::Ribbons:
    case RepresentationType::Cartoon:
        rep->type_is_valid = representation_particle_type_valid(sc, rep->type);
        break;
    case RepresentationType::ElectronicStructure:
        // Evaluated from the quantum chemistry of the first system alone, until it is per system (the GPU data and the
        // evaluations are of the first): of another, nothing is drawn rather than the first's orbitals on its atoms
        rep->type_is_valid = &sc == &state->system_context && electronic_structure_source_supported(es_source_mask(sys), rep->electronic_structure.source);
        if (rep->type_is_valid) {
            // what it waits for without blocking, asked again next frame
            again = !electronic_structure_representation_update(state, rep, colors);
        }
        break;
    case RepresentationType::DipoleMoment:
        rep->type_is_valid = md_attributes_get(&sc.data.sys.attributes, rep->dipole.dipole_key) != NULL;
        break;
    case RepresentationType::SimulationBox:
        rep->type_is_valid = md_unitcell_flags(&sc.data.state.unitcell) != 0;
        break;
    case RepresentationType::Bonds:
        rep->type_is_valid = sys.atom.count > 0;
        break;
    default:
        ASSERT(false);
        break;
    }

    if (colors) {
        if (rep->filter.result.success) {
            filter_colors(colors, num_atoms, rep->atom_mask);
            sc.representation.atom_visibility_mask_dirty = true;
            md_gl_rep_set_atom_colors(rep->md_rep, 0, (uint32_t)num_atoms, colors, 0);
            representation_set_weak_bonds(sc, rep);
            rep->glyphs_shown = false;
            if (rep->type == RepresentationType::Cartoon) {
                rep->glyphs_shown = cartoon_glyphs_rep_set_colors(sc.data.cartoon_glyphs, rep->glyph_rep, colors, num_atoms, rep->two_head_spheres);
            }
            sc.representation.gl_version += 1;

#if EXPERIMENTAL_GFX_API
            md_gfx_rep_attr_t attributes = {};
            attributes.spacefill.radius_scale = 1.0f;
            md_gfx_rep_set_type_and_attr(rep->gfx_rep, MD_GFX_REP_TYPE_SPACEFILL, &attributes);
            md_gfx_rep_set_color(rep->gfx_rep, 0, (uint32_t)mol.atom.count, (md_gfx_color_t*)colors, 0);
#endif
        }
    }

    rep->needs_update = again;
}

void representation_set_weak_bonds(SystemContext& sc, Representation* rep) {
    ApplicationState* state = sc.app;
    ASSERT(state && rep);
    if (!rep->md_rep.id) return;
    const md_system_t& sys = sc.data.sys;

    const md_atom_pair_t* pairs = nullptr;
    size_t count = 0;
    if (rep->type == RepresentationType::Licorice || rep->type == RepresentationType::BallAndStick) {
        pairs = sys.bond.pairs ? sys.bond.pairs + sys.bond.count : nullptr;
        count = sys.bond.noncovalent_count;
    } else if (rep->type == RepresentationType::Bonds) {
        pairs = rep->bond_pairs;
        count = md_array_size(rep->bond_pairs);
    }

    md_temp_scope_t temp = md_temp_begin_in(state->allocator.frame);
    defer { md_temp_end(temp); };
    uint32_t* atom_a = count ? (uint32_t*)md_vm_arena_push(state->allocator.frame, sizeof(uint32_t) * count) : nullptr;
    uint32_t* atom_b = count ? (uint32_t*)md_vm_arena_push(state->allocator.frame, sizeof(uint32_t) * count) : nullptr;
    for (size_t i = 0; i < count; ++i) {
        atom_a[i] = (uint32_t)pairs[i].idx[0];
        atom_b[i] = (uint32_t)pairs[i].idx[1];
    }
    md_gl_rep_set_weak_bonds(rep->md_rep, (uint32_t)count, atom_a, atom_b, nullptr);
}

// The dashes of a weak bond drawn dashed, and the gaps between them, as multiples of its radius
#define WEAK_BOND_DASH_LENGTH 2.0f
#define WEAK_BOND_GAP_LENGTH  1.5f

void representation_weak_bond_draw_ops(md_array(md_gl_draw_op_t)* ops, const SystemContext& sys, const Representation& rep, md_allocator_i* alloc) {
    ASSERT(ops);
    md_gl_draw_op_t op = {
        .type = MD_GL_REP_WEAK_BONDS,
        .args = {},
        .rep = rep.md_rep,
    };
    bool dashed = true;
    if (rep.type == RepresentationType::Licorice || rep.type == RepresentationType::BallAndStick) {
        // The non-covalent bonds of the system among its atoms, thinner than its sticks, in their colours
        if (!rep.noncovalent_bonds || sys.data.sys.bond.noncovalent_count == 0) return;
        const md_gl_draw_op_t sticks = representation_draw_op(rep);
        const bool licorice = rep.type == RepresentationType::Licorice;
        op.args.weak_bonds.radius = 0.5f * 0.2f * (licorice ? rep.scale.x : rep.scale.y);
        op.args.weak_bonds.color_mode    = licorice ? sticks.args.licorice.color_mode    : sticks.args.ball_and_stick.color_mode;
        op.args.weak_bonds.sharpness     = licorice ? sticks.args.licorice.sharpness     : sticks.args.ball_and_stick.sharpness;
        op.args.weak_bonds.uniform_color = licorice ? sticks.args.licorice.uniform_color : sticks.args.ball_and_stick.uniform_color;
    } else if (rep.type == RepresentationType::Bonds) {
        if (md_array_size(rep.bond_pairs) == 0) return;
        // As thick as the sticks of a licorice of the same scale
        op.args.weak_bonds.radius = 0.2f * rep.scale.x;
        op.args.weak_bonds.color_mode = (md_gl_bond_mode_t)rep.bond_color;
        op.args.weak_bonds.sharpness = rep.bond_sharpness;
        op.args.weak_bonds.uniform_color = convert_color(rep.bond_base_color);
        if (rep.tint_scale > 0.0f || rep.saturation < 1.0f) {
            tint_colors(&op.args.weak_bonds.uniform_color, 1, convert_color(rep.tint_color), rep.tint_scale, rep.saturation);
        }
        dashed = rep.bonds_dashed;
    } else {
        return;
    }
    if (dashed) {
        op.args.weak_bonds.dash_length = WEAK_BOND_DASH_LENGTH * op.args.weak_bonds.radius;
        op.args.weak_bonds.gap_length  = WEAK_BOND_GAP_LENGTH  * op.args.weak_bonds.radius;
    } else {
        // Solid: md_gl draws two dashes at least, so one too long for any bond with a gap too small to see
        op.args.weak_bonds.dash_length = 1.0e6f;
        op.args.weak_bonds.gap_length  = 1.0e-6f;
    }
    md_array_push(*ops, op, alloc);
}

void update_all_representations(SystemContext& sys) {
    // The contexts first, so that every representation is evaluated within an up to date one. Hidden groups
    // wait: what is within them is hidden too, and nothing evaluated there would be seen.
    for (size_t i = 0; i < md_array_size(sys.representation.groups); ++i) {
        if (representation_group_shown(sys, sys.representation.groups[i])) {
            update_group(sys, i, 0);
        }
    }
    for (size_t i = 0; i < md_array_size(sys.representation.reps); ++i) {
        Representation* rep = &sys.representation.reps[i];
        // The field it is coloured by is evaluated over frames, shown or not
        if (rep->type == RepresentationType::ElectronicStructure) electronic_structure_field_poll(rep);
        update_representation(sys, rep);
    }
}

void init_all_representations(SystemContext& sys) {
    ApplicationState* state = sys.app;
    for (size_t i = 0; i < md_array_size(sys.representation.reps); ++i) {
        init_representation(sys, &sys.representation.reps[i]);
    }
    // The contexts were selections of the system that was there before
    for (size_t i = 0; i < md_array_size(sys.representation.groups); ++i) {
        RepresentationGroup& g = sys.representation.groups[i];
        md_bitset_free(&g.atom_mask, state->allocator.persistent);
        g.filt_is_dirty = true;
    }
    tree_normalize(sys);
}

void flag_representation_as_dirty(Representation* rep) {
    ASSERT(rep);
    rep->filt_is_dirty = true;
    rep->needs_update  = true;
    rep->electronic_structure.vol_hash = 0;
}

bool representation_can_show_bond(const SystemContext& sys, const Representation& rep) {
    return rep.type == RepresentationType::Bonds && representation_context(sys, rep.node.parent).atoms == nullptr;
}

bool representation_shows_bond(const Representation& rep, md_atom_pair_t pair) {
    for (size_t i = 0; i < md_array_size(rep.bond_pairs); ++i) {
        const md_atom_pair_t p = rep.bond_pairs[i];
        if ((p.idx[0] == pair.idx[0] && p.idx[1] == pair.idx[1]) || (p.idx[0] == pair.idx[1] && p.idx[1] == pair.idx[0])) return true;
    }
    return false;
}

bool representation_show_bond(SystemContext& sys, md_atom_pair_t pair, md_bond_flags_t flags, Representation* into) {
    ApplicationState* app = sys.app;
    ASSERT(app);
    char expr[64];
    snprintf(expr, sizeof(expr), "bond(%d, %d)", pair.idx[0] + 1, pair.idx[1] + 1);

    if (into) {
        ASSERT(representation_can_show_bond(*app, *into));
        char text[sizeof(into->filter.text)];
        const bool empty = str_empty(str_trim(str_from_cstr(into->filter.text)));
        const int n = empty ? snprintf(text, sizeof(text), "%s", expr) : snprintf(text, sizeof(text), "(%s) or %s", into->filter.text, expr);
        if (n < 0 || (size_t)n >= sizeof(text)) return false;
        MEMCPY(into->filter.text, text, (size_t)n + 1);
        flag_representation_as_dirty(into);
        return true;
    }

    // Distinct from what draws the bond already: one colour, dashed as a weak bond is, and a stick over that of a
    // licorice (whose non-covalent bonds are half as thick as its sticks)
    Representation* rep = create_representation(sys, RepresentationType::Bonds, ColorMapping::Uniform, str_from_cstr(expr));
    if (!rep) return false;
    const bool noncovalent = md_bond_is_noncovalent(flags);
    snprintf(rep->name, sizeof(rep->name), "bond %d-%d", pair.idx[0] + 1, pair.idx[1] + 1);
    rep->bonds_dashed    = noncovalent;
    rep->scale.x         = noncovalent ? 1.0f : 1.25f;
    rep->bond_color      = BondColorMode::Uniform;
    rep->bond_base_color = {1.0f, 0.55f, 0.1f, 1.0f};
    rep->open_properties = true;
    app->show_representation_window = true;
    return true;
}

// --- New representations from the context menu of the viewport ---

namespace {

// What a target is a whole of, at the atom it is found from
enum class TargetUnit { Residue, Part, Chain, Molecule, Lipid };

// "ALA 23": the name and number the file gives the residue
void residue_label(char* buf, size_t cap, const md_system_t& sys, size_t res_idx) {
    const str_t name = md_system_residue_name(&sys, res_idx);
    snprintf(buf, cap, STR_FMT " %d", STR_ARG(name), (int)md_system_residue_seq_id(&sys, res_idx));
}

int64_t lipid_of_residue(const md_system_t& sys, size_t res_idx) {
    const md_lipid_data_t& lipids = sys.lipid;
    for (size_t i = 0; i < lipids.count; ++i) {
        const int32_t* res = md_lipid_residues(&lipids, i);
        const size_t n = md_lipid_residue_count(&lipids, i);
        for (size_t k = 0; k < n; ++k) {
            if ((size_t)res[k] == res_idx) return (int64_t)i;
        }
    }
    return -1;
}

// The unit which holds the atom, by an id which tells it apart from the others of its kind, -1 for none: then the atom
// falls back to a finer unit (target_unit_set does the same)
int64_t target_unit_id(TargetUnit unit, size_t atom, const md_system_t& sys) {
    const md_residue_idx_t res = md_system_residue_find_by_atom_idx(&sys, atom);
    switch (unit) {
    case TargetUnit::Residue:  return res;
    case TargetUnit::Part:     return md_system_part_find_by_atom_idx(&sys, atom);
    case TargetUnit::Chain: {
        const md_chain_idx_t chain = res >= 0 ? md_system_chain_find_by_residue_idx(&sys, (size_t)res) : -1;
        return chain >= 0 ? chain : -2 - (int64_t)md_system_part_find_by_atom_idx(&sys, atom);
    }
    case TargetUnit::Molecule: return md_system_molecule_find_by_atom_idx(&sys, atom);
    case TargetUnit::Lipid:    return res >= 0 ? lipid_of_residue(sys, (size_t)res) : -1;
    }
    return -1;
}

void target_unit_set(md_bitmask_t* mask, TargetUnit unit, size_t atom, const md_system_t& sys) {
    const md_residue_idx_t res = md_system_residue_find_by_atom_idx(&sys, atom);
    switch (unit) {
    case TargetUnit::Lipid: {
        const int64_t lipid = res >= 0 ? lipid_of_residue(sys, (size_t)res) : -1;
        if (lipid < 0) break;
        const int32_t* lres = md_lipid_residues(&sys.lipid, (size_t)lipid);
        for (size_t k = 0; k < md_lipid_residue_count(&sys.lipid, (size_t)lipid); ++k) {
            const md_urange_t r = md_system_residue_atom_range(&sys, (size_t)lres[k]);
            if (r.beg < r.end) md_bitmask_set_range(mask, r.beg, r.end);
        }
        return;
    }
    case TargetUnit::Molecule: {
        const md_molecule_idx_t mol = md_system_molecule_find_by_atom_idx(&sys, atom);
        if (mol < 0) break;
        mask_set_molecule(mask, (size_t)mol, sys);
        return;
    }
    case TargetUnit::Chain:
        if (res < 0) break;
        mask_set_chain_of_res(mask, (size_t)res, sys);
        return;
    case TargetUnit::Part: {
        const md_part_idx_t part = md_system_part_find_by_atom_idx(&sys, atom);
        if (part < 0) break;
        const md_urange_t r = md_system_part_atom_range(&sys, (size_t)part);
        if (r.beg < r.end) md_bitmask_set_range(mask, r.beg, r.end);
        return;
    }
    case TargetUnit::Residue:
        break;
    }
    // The residue, or the atom where it is in none
    if (res >= 0) {
        const md_urange_t r = md_system_residue_atom_range(&sys, (size_t)res);
        if (r.beg < r.end) {
            md_bitmask_set_range(mask, r.beg, r.end);
            return;
        }
    }
    md_bitmask_set_bit(mask, (uint32_t)atom);
}

// What the menu calls the unit which holds the atom: "ALA 23", "chain A", "POPC 112", "SOL" (a molecule by its type)
void target_unit_label(char* buf, size_t cap, TargetUnit unit, size_t atom, const md_system_t& sys) {
    const md_residue_idx_t res = md_system_residue_find_by_atom_idx(&sys, atom);
    switch (unit) {
    case TargetUnit::Chain:
    case TargetUnit::Part: {
        const md_chain_idx_t chain = res >= 0 ? md_system_chain_find_by_residue_idx(&sys, (size_t)res) : -1;
        const md_part_idx_t  part  = md_system_part_find_by_atom_idx(&sys, atom);
        const md_part_kind_t kind  = part >= 0 ? md_system_part_kind(&sys, (size_t)part) : MD_PART_KIND_UNKNOWN;
        // A strand by its chain, as the files call it; a ligand by its residue, as a single one is
        if (chain >= 0 && (unit == TargetUnit::Chain || md_part_kind_is_polymer(kind))) {
            const str_t id = md_system_chain_id(&sys, (size_t)chain);
            snprintf(buf, cap, "chain " STR_FMT, STR_ARG(id));
            return;
        }
        if (part >= 0 && md_part_kind_is_polymer(kind)) {
            // As polymer(i) counts them
            int num = 0;
            for (size_t i = 0; i <= (size_t)part; ++i) num += md_part_kind_is_polymer(md_system_part_kind(&sys, i));
            snprintf(buf, cap, "polymer %d", num);
            return;
        }
        if (part >= 0) {
            // A ligand or a glycan of several residues: by the first of them
            const md_urange_t r = md_system_part_residue_range(&sys, (size_t)part);
            if (r.end - r.beg > 1) {
                char first[48];
                residue_label(first, sizeof(first), sys, r.beg);
                snprintf(buf, cap, "%s and %d more", first, (int)(r.end - r.beg - 1));
                return;
            }
        }
        break;
    }
    case TargetUnit::Molecule: {
        const md_molecule_idx_t mol = md_system_molecule_find_by_atom_idx(&sys, atom);
        const md_molecule_type_idx_t type = mol >= 0 ? md_system_molecule_type_idx(&sys, (size_t)mol) : -1;
        const str_t name = type >= 0 ? md_system_molecule_type_name(&sys, (size_t)type) : str_t{};
        if (!str_empty(name)) {
            snprintf(buf, cap, STR_FMT, STR_ARG(name));
            return;
        }
        if (res < 0 && mol >= 0) {
            snprintf(buf, cap, "molecule %d", (int)mol + 1);
            return;
        }
        break;
    }
    case TargetUnit::Lipid:
    case TargetUnit::Residue:
        break;
    }
    if (res >= 0) {
        residue_label(buf, cap, sys, (size_t)res);
    } else {
        snprintf(buf, cap, "atom %d", (int)atom + 1);
    }
}

const char* target_unit_plural(TargetUnit unit) {
    switch (unit) {
    case TargetUnit::Residue:  return "residues";
    case TargetUnit::Part:     return "parts";
    case TargetUnit::Chain:    return "chains";
    case TargetUnit::Molecule: return "molecules";
    case TargetUnit::Lipid:    return "lipids";
    }
    return "";
}

// The unit a hit on the atom is taken as, at Atom granularity: what one most likely means by pointing at it. A cartoon
// is a strand, the atoms of a protein or nucleic acid a residue of it, a lipid the molecule (all of its residues), a
// ligand or a glycan its part, unless it is a residue within a polymer (a modified one).
TargetUnit target_unit_guess(size_t atom, bool cartoon, const md_system_t& sys) {
    const md_residue_idx_t res = md_system_residue_find_by_atom_idx(&sys, atom);
    const md_part_idx_t part = md_system_part_find_by_atom_idx(&sys, atom);
    const md_part_kind_t part_kind = part >= 0 ? md_system_part_kind(&sys, (size_t)part) : MD_PART_KIND_UNKNOWN;
    if (res < 0) {
        // Nothing but the bonds to go by (an XYZ file): the molecule
        return TargetUnit::Molecule;
    }
    switch (md_system_residue_kind(&sys, (size_t)res)) {
    case MD_RESIDUE_KIND_LIPID:
        return TargetUnit::Lipid;
    case MD_RESIDUE_KIND_AMINO_ACID:
    case MD_RESIDUE_KIND_NUCLEOTIDE:
        return (cartoon && part >= 0 && md_part_kind_is_polymer(part_kind)) ? TargetUnit::Part : TargetUnit::Residue;
    case MD_RESIDUE_KIND_WATER:
    case MD_RESIDUE_KIND_ION:
        return TargetUnit::Residue;
    default:
        if (cartoon && part >= 0 && md_part_kind_is_polymer(part_kind)) return TargetUnit::Part;
        // A ligand or a glycan, a part of its own; a residue within a strand by itself
        return (part >= 0 && !md_part_kind_is_polymer(part_kind) && part_kind != MD_PART_KIND_WATER) ? TargetUnit::Part : TargetUnit::Residue;
    }
}

// All the residues which share the name of that of the atom, with the number of them. Of a residue no polymer holds.
size_t target_same_name(md_bitmask_t* mask, size_t atom, const md_system_t& sys) {
    const md_residue_idx_t res = md_system_residue_find_by_atom_idx(&sys, atom);
    if (res < 0) return 0;
    const md_residue_kind_t kind = md_system_residue_kind(&sys, (size_t)res);
    if (kind == MD_RESIDUE_KIND_AMINO_ACID || kind == MD_RESIDUE_KIND_NUCLEOTIDE) return 0;
    const md_part_idx_t part = md_system_part_find_by_atom_idx(&sys, atom);
    if (part >= 0 && md_part_kind_is_polymer(md_system_part_kind(&sys, (size_t)part))) return 0;

    const str_t name = md_system_residue_name(&sys, (size_t)res);
    if (str_empty(name)) return 0;
    size_t num = 0;
    for (size_t i = 0; i < md_system_residue_count(&sys); ++i) {
        if (!str_eq(md_system_residue_name(&sys, i), name)) continue;
        const md_urange_t r = md_system_residue_atom_range(&sys, i);
        if (r.beg < r.end) {
            md_bitmask_set_range(mask, r.beg, r.end);
            num += 1;
        }
    }
    return num;
}

// What a target can be shown as, and what it is suggested to be shown as: by what it is, and not as something which
// draws most of it already, which a new one of the same type would only coincide with
void target_suggest(RepresentationTarget* t, const SystemContext& sc) {
    const md_system_t& sys = sc.data.sys;
    const md_bitset_t atoms = t->atoms;
    const size_t num_atoms = md_bitset_count(atoms);

    // What the cartoon would draw of it: the residues of a protein backbone, and the points of the glyphs (lipids and
    // nucleic acids)
    size_t num_backbone = 0;
    const md_protein_backbone_data_t& bb = sys.protein_backbone;
    for (size_t i = 0; i < bb.segment.count; ++i) {
        const md_atom_idx_t ca = bb.segment.atoms[i].ca;
        num_backbone += ca >= 0 && md_bitset_test(atoms, (uint32_t)ca);
    }
    size_t num_glyphs = 0;
    const CartoonGlyphs& glyphs = sc.data.cartoon_glyphs;
    if (cartoon_glyphs_valid(glyphs)) {
        for (uint32_t i = 0; i < glyphs.num_points; ++i) {
            const int32_t a = glyphs.point_color_atom[i];
            num_glyphs += a >= 0 && md_bitset_test(atoms, (uint32_t)a);
        }
    }
    t->ribbons = num_backbone > 0 && representation_particle_type_valid(sc, RepresentationType::Ribbons);
    t->cartoon = (num_backbone > 0 || num_glyphs > 0) && representation_particle_type_valid(sc, RepresentationType::Cartoon);

    // Mostly of a polymer: the residues of a protein or a nucleic acid
    size_t num_polymer_atoms = 0;
    size_t num_amino_atoms = 0;
    for (size_t i = 0; i < md_system_residue_count(&sys); ++i) {
        const md_urange_t r = md_system_residue_atom_range(&sys, i);
        if (r.beg >= r.end || !md_bitset_test_any_range(atoms, r.beg, r.end)) continue;
        const md_residue_kind_t kind = md_system_residue_kind(&sys, i);
        const size_t n = md_bitset_count_range(atoms, r.beg, r.end);
        num_amino_atoms   += kind == MD_RESIDUE_KIND_AMINO_ACID ? n : 0;
        num_polymer_atoms += (kind == MD_RESIDUE_KIND_AMINO_ACID || kind == MD_RESIDUE_KIND_NUCLEOTIDE) ? n : 0;
    }

    // By preference: a strand (or more) as a cartoon, a lipid as its glyph; the residues of a polymer, which a cartoon
    // likely shows already, as sticks; a ligand as balls and sticks; something large (solvent, a membrane) as spheres
    RepresentationType prefs[4];
    size_t num_prefs = 0;
    const bool lipids = num_glyphs > 0 && num_backbone == 0 && num_polymer_atoms == 0;
    if (t->cartoon && (num_backbone >= 8 || lipids || (num_glyphs >= 8))) {
        prefs[num_prefs++] = RepresentationType::Cartoon;
        prefs[num_prefs++] = RepresentationType::Licorice;
        prefs[num_prefs++] = RepresentationType::SpaceFill;
        prefs[num_prefs++] = RepresentationType::BallAndStick;
    } else if (num_atoms > 20000) {
        prefs[num_prefs++] = RepresentationType::SpaceFill;
        prefs[num_prefs++] = RepresentationType::Licorice;
        prefs[num_prefs++] = RepresentationType::BallAndStick;
    } else if (2 * num_polymer_atoms > num_atoms) {
        prefs[num_prefs++] = RepresentationType::Licorice;
        prefs[num_prefs++] = RepresentationType::BallAndStick;
        prefs[num_prefs++] = RepresentationType::SpaceFill;
    } else {
        prefs[num_prefs++] = RepresentationType::BallAndStick;
        prefs[num_prefs++] = RepresentationType::SpaceFill;
        prefs[num_prefs++] = RepresentationType::Licorice;
    }

    // The types which draw most of it already
    bool drawn[(int)RepresentationType::Count] = {};
    md_temp_scope_t temp = md_temp_begin();
    for (size_t i = 0; i < md_array_size(sc.representation.reps); ++i) {
        const Representation& rep = sc.representation.reps[i];
        if ((int)rep.type > (int)RepresentationType::Cartoon || !representation_shown(sc, rep)) continue;
        if (drawn[(int)rep.type]) continue;
        const md_bitset_t both = md_bitset_and(rep.atom_mask, atoms, md_temp_allocator(temp));
        drawn[(int)rep.type] = 2 * md_bitset_count(both) > num_atoms;
    }
    md_temp_end(temp);

    // The first not drawn already, else the first that can be shown at all
    t->suggested_type = RepresentationType::SpaceFill;
    bool found = false;
    for (int pass = 0; pass < 2 && !found; ++pass) {
        for (size_t i = 0; i < num_prefs; ++i) {
            if ((pass == 1 || !drawn[(int)prefs[i]]) && representation_target_can_show(sc, *t, prefs[i])) {
                t->suggested_type = prefs[i];
                found = true;
                break;
            }
        }
    }
    t->suggested_color = (t->suggested_type == RepresentationType::Cartoon)
        ? (2 * num_amino_atoms > num_atoms ? ColorMapping::SecondaryStructure : ColorMapping::ResType)
        : ColorMapping::Type;
}

// Done with the mask: its atoms the target's, a single atom none
bool target_finish(RepresentationTarget* t, md_bitmask_t* mask, const SystemContext& sys) {
    const bool ok = md_bitmask_count(mask) > 1;
    if (ok) {
        t->atoms = md_bitmask_extract(mask, md_get_heap_allocator());
        target_suggest(t, sys);
    }
    md_bitmask_free(mask);
    return ok;
}

}  // namespace

bool representation_targets_from_hit(RepresentationTarget* target, RepresentationTarget* all, const SystemContext& sc, const PickingHit& hit) {
    const ApplicationState& app = *sc.app;
    ASSERT(target);
    representation_target_free(target);
    if (all) representation_target_free(all);

    const md_system_t& sys = sc.data.sys;
    if (sys.atom.count == 0) return false;

    // The atoms the hit stands for, and whether it is of a cartoon (a strand rather than a residue of it)
    int32_t seeds[64];
    size_t num_seeds = 0;
    bool cartoon = false;
    auto push = [&](int64_t atom) {
        if (atom >= 0 && (size_t)atom < sys.atom.count && num_seeds < ARRAY_SIZE(seeds)) seeds[num_seeds++] = (int32_t)atom;
    };
    auto push_residue = [&](int64_t res) {
        if (res >= 0 && (size_t)res < sys.residue.count) push(md_system_residue_atom_range(&sys, (size_t)res).beg);
    };

    if (hit.domain == PickingDomain_Atom) {
        push(hit.local_idx);
    } else if (hit.domain == PickingDomain_Bond) {
        md_atom_pair_t pair;
        md_bond_flags_t flags;
        bool formed;
        if (picked_bond(sc, hit.local_idx, &pair, &flags, &formed)) {
            push(pair.idx[0]);
            push(pair.idx[1]);
        }
    } else if (hit.domain == PickingDomain_Ring) {
        const int32_t* ring_atoms = nullptr;
        size_t ring_count = 0;
        md_ring_flags_t ring_flags = MD_RING_FLAG_NONE;
        if (picked_ring(sc, hit.local_idx, &ring_atoms, &ring_count, &ring_flags)) {
            for (size_t i = 0; i < ring_count; ++i) push(ring_atoms[i]);
        }
    } else if (hit.domain == PickingDomain_BackboneSegment) {
        const md_protein_backbone_data_t& bb = sys.protein_backbone;
        if (hit.local_idx < bb.segment.count) {
            cartoon = true;
            const md_residue_idx_t res = bb.segment.res_idx ? bb.segment.res_idx[hit.local_idx] : -1;
            if (res >= 0) push_residue(res);
            else push(bb.segment.atoms[hit.local_idx].ca);
        }
    } else if (hit.domain == PickingDomain_CartoonGlyph) {
        const CartoonGlyphs& glyphs = sc.data.cartoon_glyphs;
        const uint32_t point = cartoon_glyphs_picked_point(glyphs, hit.local_idx);
        if (point < glyphs.num_points) {
            cartoon = true;
            const int32_t res = cartoon_glyphs_point_residue(glyphs, point);
            if (res >= 0) push_residue(res);
            else push(glyphs.point_color_atom[point]);
        }
    } else if (hit.domain == PickingDomain_ChainGap) {
        cartoon = true;
        push_residue(hit.local_idx);
        push_residue((int64_t)hit.local_idx + 1);
    }
    if (num_seeds == 0) return false;

    md_bitmask_t mask = md_bitmask_create(0, (uint32_t)sys.atom.count, md_get_heap_allocator());

    const SelectionGranularity granularity = app.selection.granularity;
    TargetUnit unit = TargetUnit::Residue;
    if (granularity != SelectionGranularity::Atom) {
        // What a click there selects, which the highlight showed
        switch (granularity) {
        case SelectionGranularity::Part:     unit = md_system_part_count(&sys) ? TargetUnit::Part : TargetUnit::Residue; break;
        case SelectionGranularity::Chain:    unit = TargetUnit::Chain; break;
        case SelectionGranularity::Molecule: unit = sys.molecule.count ? TargetUnit::Molecule : TargetUnit::Part; break;
        default:                             unit = TargetUnit::Residue; break;
        }
        for (size_t i = 0; i < num_seeds; ++i) {
            mask_set_atom_by_selection_granularity(&mask, (size_t)seeds[i], granularity, sys);
        }
    } else {
        unit = target_unit_guess((size_t)seeds[0], cartoon, sys);
        for (size_t i = 0; i < num_seeds; ++i) {
            target_unit_set(&mask, unit, (size_t)seeds[i], sys);
        }
    }

    // Its name: that of the unit, or how many of them there are (a bond between two residues)
    int64_t ids[ARRAY_SIZE(seeds)];
    size_t num_ids = 0;
    for (size_t i = 0; i < num_seeds; ++i) {
        const int64_t id = target_unit_id(unit, (size_t)seeds[i], sys);
        bool found = false;
        for (size_t k = 0; k < num_ids; ++k) found |= ids[k] == id;
        if (!found) ids[num_ids++] = id;
    }
    if (num_ids > 1) {
        snprintf(target->label, sizeof(target->label), "%d %s", (int)num_ids, target_unit_plural(unit));
    } else {
        target_unit_label(target->label, sizeof(target->label), unit, (size_t)seeds[0], sys);
    }

    // All of its name, where it is one of several (and the one alone, where it is a single atom: an ion)
    RepresentationTarget same = {};
    bool same_ok = false;
    {
        md_bitmask_t same_mask = md_bitmask_create(0, (uint32_t)sys.atom.count, md_get_heap_allocator());
        const size_t num = target_same_name(&same_mask, (size_t)seeds[0], sys);
        const md_residue_idx_t res = md_system_residue_find_by_atom_idx(&sys, (size_t)seeds[0]);
        if (num > 1 && res >= 0) {
            const str_t name = md_system_residue_name(&sys, (size_t)res);
            snprintf(same.label, sizeof(same.label), "all " STR_FMT, STR_ARG(name));
            same_ok = target_finish(&same, &same_mask, sc);
        } else {
            md_bitmask_free(&same_mask);
        }
    }

    bool ok = target_finish(target, &mask, sc);
    if (ok && same_ok && md_bitset_equal(target->atoms, same.atoms)) {
        // The only one of its name: it is all of them
        representation_target_free(&same);
        same_ok = false;
    }
    if (!ok && same_ok) {
        // A single atom, such as an ion: all of its name instead
        *target = same;
        same = {};
        same_ok = false;
        ok = true;
    }
    if (all && same_ok) {
        *all = same;
    } else {
        representation_target_free(&same);
    }
    return ok;
}

bool representation_target_from_selection(RepresentationTarget* target, const SystemContext& sys) {
    ASSERT(target);
    representation_target_free(target);
    if (md_bitmask_count(&sys.selection.selection_mask) < 2) return false;
    target->atoms = md_bitmask_extract(&sys.selection.selection_mask, md_get_heap_allocator());
    snprintf(target->label, sizeof(target->label), "selection");
    target->selection = true;
    target_suggest(target, sys);
    return true;
}

void representation_target_free(RepresentationTarget* target) {
    ASSERT(target);
    md_bitset_free(&target->atoms, md_get_heap_allocator());
    *target = RepresentationTarget{};
}

bool representation_target_can_show(const SystemContext& sys, const RepresentationTarget& target, RepresentationType type) {
    switch (type) {
    case RepresentationType::SpaceFill:
    case RepresentationType::Licorice:
    case RepresentationType::BallAndStick:
        return representation_particle_type_valid(sys, type);
    case RepresentationType::Ribbons: return target.ribbons;
    case RepresentationType::Cartoon: return target.cartoon;
    default: return false;
    }
}

Representation* create_representation_of_target(SystemContext& sys, const RepresentationTarget& target, RepresentationType type) {
    ApplicationState* app = sys.app;
    ASSERT(app);
    if (md_bitset_empty(target.atoms)) return nullptr;

    // The most telling way of writing it which fits
    md_temp_scope_t temp = md_temp_begin();
    char filter[sizeof(AtomFilter::text)] = "";
    const md_array(ScriptSuggestion) suggestions = script_suggest_selection(app->script.language, STR_LIT("x"), target.atoms, sys.data.sys, md_temp_allocator(temp));
    for (size_t i = 0; i < md_array_size(suggestions); ++i) {
        str_t expr = suggestions[i].text;
        if (!str_begins_with(expr, STR_LIT("x = "))) continue;
        expr = str_trim(str_substr(expr, 4));
        if (str_ends_with(expr, STR_LIT(";"))) expr = str_trim(str_substr(expr, 0, expr.len - 1));
        if (!str_empty(expr) && expr.len < sizeof(filter)) {
            str_copy_to_char_buf(filter, sizeof(filter), expr);
            break;
        }
    }
    md_temp_end(temp);
    if (filter[0] == '\0') {
        VIAMD_LOG_INFO("'%s' is too scattered to write as the filter of a representation", target.label);
        return nullptr;
    }

    ColorMapping color = ColorMapping::Type;
    if (type == RepresentationType::Cartoon || type == RepresentationType::Ribbons) {
        // As suggested where it was suggested a cartoon, else by the secondary structure of a protein
        color = target.suggested_type == RepresentationType::Cartoon ? target.suggested_color
              : (target.ribbons ? ColorMapping::SecondaryStructure : ColorMapping::ResType);
    }
    Representation* rep = create_representation(sys, type, color, str_from_cstr(filter));
    if (!rep) return nullptr;
    // A selection by the expression, where that is short enough to read as a name
    if (target.selection && strlen(filter) < 32) {
        snprintf(rep->name, sizeof(rep->name), "%s", filter);
    } else {
        snprintf(rep->name, sizeof(rep->name), "%s", target.label);
    }
    // Opened, to be styled, when the window lists it
    rep->open_properties = true;
    return rep;
}

void flag_all_representations_as_dirty(SystemContext& sys) {
    for (size_t i = 0; i < md_array_size(sys.representation.groups); ++i) {
        sys.representation.groups[i].filt_is_dirty = true;
    }
    for (size_t i = 0; i < md_array_size(sys.representation.reps); ++i) {
        flag_representation_as_dirty(&sys.representation.reps[i]);
    }
}

void flag_dynamic_representations_as_dirty(SystemContext& sys) {
    // A context that follows the frame is evaluated again; what is within it decides for itself whether to
    // follow, by its own auto-update
    for (size_t i = 0; i < md_array_size(sys.representation.groups); ++i) {
        RepresentationGroup& g = sys.representation.groups[i];
        if (g.filter.result.is_dynamic && representation_group_shown(sys, g)) {
            g.filt_is_dirty = true;
        }
    }
    for (size_t i = 0; i < md_array_size(sys.representation.reps); ++i) {
        Representation& rep = sys.representation.reps[i];
        if (!representation_shown(sys, rep)) continue;

        const bool atoms_follow  = rep.dynamic_evaluation && representation_is_dynamic(sys, rep);
        // The volume of an electronic structure is evaluated around the atoms where they are; it is
        // evaluated again only if the frame actually moved them (vol_hash)
        const bool colors_follow = rep.color_mapping == ColorMapping::SecondaryStructure ||
                                   rep.type == RepresentationType::ElectronicStructure ||
                                   (rep.dynamic_evaluation && rep.color_mapping == ColorMapping::Attribute);
        if (atoms_follow) {
            flag_representation_as_dirty(&rep);
        } else if (colors_follow) {
            rep.needs_update = true;
        }
    }
}

// The nucleotides on a nucleic backbone, atomistic or coarse grained: those the cartoon shows
static size_t num_nucleotides_on_backbone(const md_system_t& sys) {
    size_t num = 0;
    for (size_t i = 0; i < sys.nucleic_backbone.segment.count; ++i) {
        num += sys.nucleic_backbone.segment.atoms[i].c4 >= 0;
    }
    return num;
}

// What a system starts out shown with: by what is in it
void create_default_representations(SystemContext& sc) {
    ApplicationState* state = sc.app;
    bool amino_acid_present = false;
    bool nucleic_present = false;
    bool lipid_present = false;
    bool lipid_cartoon = false;     // The lipids as their cartoon (a sphere for each head, a tube for each tail)
    size_t num_nucleotides = 0;
    bool ion_present = false;
    bool water_present = false;
    bool ligand_present = false;
    size_t num_coarse_grained = 0;
    // Any electronic structure, not orbitals in particular: a file can carry density properties and
    // no SCF block, and its density properties are as much something to show as orbitals would be.
    // Only the first system's quantum chemistry is evaluated so far (update_representation)
    const bool electronic_structure_present = &sc == &state->system_context && es_source_mask(sc.data.sys) != 0;


    if (sc.data.sys.atom.count > 3'000'000) {
        VIAMD_LOG_INFO("Large system detected, creating default representation for all atoms");
        Representation* rep = create_representation(sc, RepresentationType::SpaceFill, ColorMapping::Type);
        snprintf(rep->name, sizeof(rep->name), "default");
        goto done;
    }

    if (sc.data.sys.residue.count == 0) {
        // No residues present
        Representation* rep = create_representation(sc, RepresentationType::BallAndStick, ColorMapping::Type);
        snprintf(rep->name, sizeof(rep->name), "default");
        goto done;
    }

    // What is present, by what the residues are. The selections of the representations below are by the same.
    for (size_t i = 0; i < sc.data.sys.residue.count; ++i) {
        switch (md_system_residue_kind(&sc.data.sys, i)) {
        case MD_RESIDUE_KIND_AMINO_ACID: amino_acid_present = true; break;
        case MD_RESIDUE_KIND_NUCLEOTIDE: nucleic_present = true; num_nucleotides += 1; break;
        case MD_RESIDUE_KIND_ION:          ion_present = true;        break;
        case MD_RESIDUE_KIND_WATER:        water_present = true;      break;
        case MD_RESIDUE_KIND_LIPID:        lipid_present = true;      break;
        default:                           ligand_present = true;     break;
        }
    }
    for (size_t i = 0; i < sc.data.sys.atom.count; ++i) {
        num_coarse_grained += md_atom_particle_kind(&sc.data.sys.atom, i) == MD_PARTICLE_BEAD;
    }

    // Lipids are shown as their cartoon where their heads and tails are known (the atoms of 'head or tail'), and the
    // rest of them (lipids by name whose anatomy is not known, as without bonds) as ligands
    if (lipid_present) {
        const md_lipid_data_t& lipids = sc.data.sys.lipid;
        size_t num_lipid_residues = 0;
        size_t num_known = 0;
        for (size_t i = 0; i < sc.data.sys.residue.count; ++i) {
            num_lipid_residues += md_system_residue_kind(&sc.data.sys, i) == MD_RESIDUE_KIND_LIPID;
        }
        for (size_t i = 0; i < lipids.count; ++i) {
            if (lipids.head_offset[i + 1] > lipids.head_offset[i]) num_known += md_lipid_residue_count(&lipids, i);
        }
        lipid_cartoon = num_known > 0 && cartoon_glyphs_valid(sc.data.cartoon_glyphs);
        if (!lipid_cartoon || num_known < num_lipid_residues) ligand_present = true;
    }
    if (lipid_cartoon) {
        Representation* lip = create_representation(sc, RepresentationType::Cartoon, ColorMapping::ResType, STR_LIT("lipid"));
        snprintf(lip->name, sizeof(lip->name), "lipid");
    }

    // Coarse grained when most of it is. A handful of beads, or atoms a loader could not assign an
    // element to, should not take the protein, ligand and water representations from the rest.
    if (2 * num_coarse_grained > sc.data.sys.atom.count) {
        // Proteins of beads have a trace backbone and its secondary structure (CG-DSSP): a cartoon. So do nucleic acids
        // of beads where every nucleotide is on a backbone, their bases with the rings of ideal ones. The water as a
        // representation of its own, hidden unless there is nothing else, as for atoms: the water beads of a solvated
        // membrane are most of it. The rest as beads.
        const bool protein_trace = amino_acid_present && sc.data.sys.protein_backbone.segment.count > 0;
        const size_t num_on_backbone = num_nucleotides_on_backbone(sc.data.sys);
        const bool nucleic_cartoon = nucleic_present && num_on_backbone >= 10 && num_on_backbone == num_nucleotides && cartoon_glyphs_valid(sc.data.cartoon_glyphs);
        if (protein_trace) {
            Representation* prot = create_representation(sc, RepresentationType::Cartoon, ColorMapping::SecondaryStructure, STR_LIT("protein"));
            snprintf(prot->name, sizeof(prot->name), "protein");
        }
        if (nucleic_cartoon) {
            Representation* nucl = create_representation(sc, RepresentationType::Cartoon, ColorMapping::ResType, STR_LIT("nucleic"));
            snprintf(nucl->name, sizeof(nucl->name), "nucleic");
        }
        const char* shown[4];
        int num_shown = 0;
        if (protein_trace)   shown[num_shown++] = "protein";
        if (nucleic_cartoon) shown[num_shown++] = "nucleic";
        if (lipid_cartoon)   shown[num_shown++] = state->script.language == script::Language::V2 ? "particle(head) or particle(tail)" : "head or tail";
        if (water_present)   shown[num_shown++] = "water";
        char rest[160] = "";
        if (num_shown == 1 && !strchr(shown[0], ' ')) {
            snprintf(rest, sizeof(rest), "not %s", shown[0]);
        } else if (num_shown > 0) {
            int len = snprintf(rest, sizeof(rest), "not (");
            for (int i = 0; i < num_shown; ++i) {
                len += snprintf(rest + len, sizeof(rest) - len, "%s%s", i > 0 ? " or " : "", shown[i]);
            }
            snprintf(rest + len, sizeof(rest) - len, ")");
        }
        Representation* rep = create_representation(sc, RepresentationType::SpaceFill, ColorMapping::Type, str_from_cstr(rest));
        snprintf(rep->name, sizeof(rep->name), "default");
        if (water_present) {
            Representation* water = create_representation(sc, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("water"));
            snprintf(water->name, sizeof(water->name), "water");
            water->enabled = !amino_acid_present && !nucleic_present && !ligand_present && !lipid_present;
        }
        goto done;
    }

    if (amino_acid_present) {
        RepresentationType type = RepresentationType::Cartoon;
        ColorMapping color = ColorMapping::SecondaryStructure;

        // Several strands are told apart by color; a short single strand (or none: free amino acids) is shown atom by atom
        const md_system_t& sys = sc.data.sys;
        size_t num_chains = 0;
        size_t chain_res_count = 0;
        for (size_t i = 0; i < md_system_part_count(&sys); ++i) {
            if (md_system_part_kind(&sys, i) == MD_PART_KIND_PEPTIDE) {
                if (num_chains++ == 0) {
                    const md_urange_t r = md_system_part_residue_range(&sys, i);
                    chain_res_count = r.end - r.beg;
                }
            }
        }
        if (num_chains > 1) {
            // By the chain ids the data gives, else by the strands themselves
            color = md_system_chain_count(&sys) > 0 ? ColorMapping::ChainId : ColorMapping::PartIndex;
        } else if (chain_res_count < 20) {
            type = RepresentationType::BallAndStick;
            color = ColorMapping::Type;
        }

        Representation* prot = create_representation(sc, type, color, STR_LIT("protein"));
        snprintf(prot->name, sizeof(prot->name), "protein");
    }
    if (nucleic_present) {
        // A cartoon where every nucleotide is on a backbone (and there are some of them), else atom by atom: the cartoon
        // shows nothing of a nucleotide off a backbone
        const size_t num_on_backbone = num_nucleotides_on_backbone(sc.data.sys);
        const bool cartoon = num_on_backbone >= 10 && num_on_backbone == num_nucleotides && cartoon_glyphs_valid(sc.data.cartoon_glyphs);
        Representation* nucl = create_representation(sc, cartoon ? RepresentationType::Cartoon : RepresentationType::BallAndStick, cartoon ? ColorMapping::ResType : ColorMapping::Type, STR_LIT("nucleic"));
        snprintf(nucl->name, sizeof(nucl->name), "nucleic");
    }
    if (ion_present) {
        Representation* ion = create_representation(sc, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("ion"));
        snprintf(ion->name, sizeof(ion->name), "ion");
    }
    if (ligand_present) {
        // The lipids whose anatomy is known are their cartoon: VIAL's heads and tails are units, all of them particles(...)
        const bool v2 = state->script.language == script::Language::V2;
        const str_t filter = !lipid_cartoon ? STR_LIT("not (protein or nucleic or water or ion)")
                           : v2 ? STR_LIT("not (protein or nucleic or water or ion or particle(head) or particle(tail))")
                                : STR_LIT("not (protein or nucleic or water or ion or head or tail)");
        Representation* ligand = create_representation(sc, RepresentationType::BallAndStick, ColorMapping::Type, filter);
        snprintf(ligand->name, sizeof(ligand->name), "ligand");
    }
    if (water_present) {
        Representation* water = create_representation(sc, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("water"));
        water->scale.x = 0.5f;
        snprintf(water->name, sizeof(water->name), "water");
        water->enabled = false;
        if (!amino_acid_present && !nucleic_present && !ligand_present && !lipid_present) {
            water->enabled = true;
        }
    }

done:
    // The cell, where there is one, hidden: there to be shown, for the coarse grained and the large too
    if (md_unitcell_flags(&sc.data.state.unitcell) != 0) {
        Representation* box = create_representation(sc, RepresentationType::SimulationBox);
        snprintf(box->name, sizeof(box->name), "box");
        box->enabled = false;
    }

    if (electronic_structure_present) {
        Representation* rep = create_representation(sc, RepresentationType::ElectronicStructure);
        snprintf(rep->name, sizeof(rep->name), "electronic structure");
        rep->enabled = true;
        electronic_structure_select_available_source(&rep->electronic_structure, sc.data.sys);


		// ONE representation, for the ground state moment. A group each meant a VeloxChem file with
		// transition dipoles opened with a representation per group, and with the electric,
		// magnetic and velocity sets that is three nobody asked for on top of the one they wanted.
		// The ground state is what a dipole means to someone who has not said otherwise; the rest
		// are a representation away, and the index slider covers the excited states within each.
		DipoleGroup groups[16];
		size_t num_groups = MIN(dipole_groups_gather(groups, ARRAY_SIZE(groups), sc.data.sys), ARRAY_SIZE(groups));
        for (size_t i = 0; i < num_groups; ++i) {
            // The group name is the identity here, the same string the path spells.
            if (!str_eq(groups[i].label, STR_LIT("ground_state"))) continue;

            vec3_t vec = {0, 0, 0};
            if (!dipole_moment_read(&vec, nullptr, sc.data.sys, groups[i].key, 0)) continue;
            if (vec3_length(vec) <= 1e-3f) continue;

            Representation* dipole_rep = create_representation(sc, RepresentationType::DipoleMoment);
            dipole_rep->dipole.dipole_key   = groups[i].key;
            dipole_rep->dipole.dipole_index = 0;

            // Lower case, like every other auto created representation - "protein", "water",
            // "electronic structure". dipole_label_pretty title cases for menus and tooltips,
            // which is a different job and stays as it is.
            snprintf(dipole_rep->name, sizeof(dipole_rep->name), "dipole moment");
            dipole_rep->enabled = true;
            break;
        }
    }

    recompute_atom_visibility_mask(sc);
}

// ---------------------------------------------------------------------------------------------------------
// Workspace
// ---------------------------------------------------------------------------------------------------------

bool representation_read_entry(Representation* rep, str_t ident, str_t arg) {
    if (str_eq(ident, STR_LIT("Name"))) {
        viamd::extract_to_char_buf(rep->name, sizeof(rep->name), arg);
    } else if (str_eq(ident, STR_LIT("Filter"))) {
        viamd::extract_to_char_buf(rep->filter.text, sizeof(rep->filter.text), arg);
    } else if (str_eq(ident, STR_LIT("Enabled"))) {
        viamd::extract_bool(rep->enabled, arg);
    } else if (str_eq(ident, STR_LIT("Type"))) {
        viamd::extract_enum(rep->type, arg, (int)RepresentationType::Count);
    } else if (str_eq(ident, STR_LIT("ColorMapping"))) {
        viamd::extract_enum(rep->color_mapping, arg, (int)ColorMapping::Count);
    } else if (str_eq(ident, STR_LIT("StaticColor")) || str_eq(ident, STR_LIT("BaseColor"))) {
        viamd::extract_vec4(rep->base_color, arg);
    } else if (str_eq(ident, STR_LIT("Saturation"))) {
        viamd::extract_flt(rep->saturation, arg);
    } else if (str_eq(ident, STR_LIT("TintColor"))) {
        viamd::extract_vec4(rep->tint_color, arg);
    } else if (str_eq(ident, STR_LIT("TintScale"))) {
        viamd::extract_flt(rep->tint_scale, arg);
    } else if (str_eq(ident, STR_LIT("SecondaryStructureColorUnknown"))) {
        viamd::extract_vec4(rep->secondary_structure.color_unknown, arg);
    } else if (str_eq(ident, STR_LIT("SecondaryStructureColorCoil"))) {
        viamd::extract_vec4(rep->secondary_structure.color_coil, arg);
    } else if (str_eq(ident, STR_LIT("SecondaryStructureColorHelix"))) {
        viamd::extract_vec4(rep->secondary_structure.color_helix, arg);
    } else if (str_eq(ident, STR_LIT("SecondaryStructureColorSheet"))) {
        viamd::extract_vec4(rep->secondary_structure.color_sheet, arg);
    } else if (str_eq(ident, STR_LIT("BondColor"))) {
        viamd::extract_enum(rep->bond_color, arg, (int)BondColorMode::Count);
    } else if (str_eq(ident, STR_LIT("BondSharpness"))) {
        viamd::extract_flt(rep->bond_sharpness, arg);
    } else if (str_eq(ident, STR_LIT("BondBaseColor"))) {
        viamd::extract_vec4(rep->bond_base_color, arg);
    } else if (str_eq(ident, STR_LIT("BondOrders"))) {
        viamd::extract_bool(rep->bond_orders, arg);
    } else if (str_eq(ident, STR_LIT("BondsDashed"))) {
        viamd::extract_bool(rep->bonds_dashed, arg);
    } else if (str_eq(ident, STR_LIT("ChainGaps"))) {
        viamd::extract_bool(rep->chain_gaps, arg);
    } else if (str_eq(ident, STR_LIT("ElasticBonds"))) {
        viamd::extract_bool(rep->elastic_bonds, arg);
    } else if (str_eq(ident, STR_LIT("NoncovalentBonds"))) {
        viamd::extract_bool(rep->noncovalent_bonds, arg);
    } else if (str_eq(ident, STR_LIT("TwoHeadSpheres"))) {
        viamd::extract_bool(rep->two_head_spheres, arg);
    } else if (str_eq(ident, STR_LIT("NucleicBandRadius"))) {
        viamd::extract_flt(rep->nucleic_band_radius, arg);
    } else if (str_eq(ident, STR_LIT("NucleicRingRadius"))) {
        viamd::extract_flt(rep->nucleic_ring_radius, arg);
    } else if (str_eq(ident, STR_LIT("Radius"))) {
        // DEPRECATED
        viamd::extract_flt(rep->scale.x, arg);
    } else if (str_eq(ident, STR_LIT("Tension"))) {
        // DEPRECATED
    } else if (str_eq(ident, STR_LIT("Width"))) {
        viamd::extract_flt(rep->scale.x, arg);
    } else if (str_eq(ident, STR_LIT("Thickness"))) {
        viamd::extract_flt(rep->scale.y, arg);
    } else if (str_eq(ident, STR_LIT("Param"))) {
        viamd::extract_vec4(rep->scale, arg);
    } else if (str_eq(ident, STR_LIT("AutoUpdate"))) {
        viamd::extract_bool(rep->dynamic_evaluation, arg);
    } else if (str_eq(ident, STR_LIT("DynamicEval"))) {
        // DEPRECATED, and not read: the window used to clear it whenever it drew a representation whose filter
        // did not depend on the frame, so a stored 0 mostly says nothing about what the user chose. Read, it
        // would keep a representation from following a context that does (a group within(...)) without a word.
        // AutoUpdate is the choice alone.
    } else if (str_eq(ident, STR_LIT("AtomicPropertyPath"))) {
        // Ids are a function of the path: this resolves before the dataset is loaded
        rep->atom_attribute.key = md_attributes_id_from_path(arg);
    } else if (str_eq(ident, STR_LIT("AtomicPropertyVariant"))) {
        viamd::extract_int(rep->atom_attribute.variant_idx, arg);
    } else if (str_eq(ident, STR_LIT("AtomicPropertyColormap"))) {
        viamd::extract_int(rep->atom_attribute.scale.colormap, arg);
    } else if (str_eq(ident, STR_LIT("AtomicPropertyRange"))) {
        // A range that was written is one somebody set: a workspace from before the range could
        // follow the values has no AtomicPropertyAutoRange, and keeps it. One that has, has it
        // after this line and says for itself.
        float r[2];
        if (viamd::extract_flt_vec(r, 2, arg)) {
            rep->atom_attribute.scale.range_beg = r[0];
            rep->atom_attribute.scale.range_end = r[1];
            rep->atom_attribute.scale.auto_range = false;
        }
    } else if (str_eq(ident, STR_LIT("AtomicPropertyAutoRange"))) {
        viamd::extract_bool(rep->atom_attribute.scale.auto_range, arg);
    } else if (str_eq(ident, STR_LIT("AtomicPropertySymmetric"))) {
        viamd::extract_bool(rep->atom_attribute.scale.symmetric, arg);
    } else if (str_eq(ident, STR_LIT("AtomicPropertyLegend"))) {
        viamd::extract_bool(rep->atom_attribute.scale.show_legend, arg);
    } else if (str_eq(ident, STR_LIT("AtomicPropertyDataRange")) || str_eq(ident, STR_LIT("AtomicPropertySymmetricZero"))) {
        // Written by earlier versions: the span is measured from the data now, and the old
        // symmetric flag never reached the colours
    } else if (str_eq(ident, STR_LIT("DipolePath"))) {
        rep->dipole.dipole_key = md_attributes_id_from_path(arg);
    } else if (str_eq(ident, STR_LIT("DipoleIndex"))) {
        int i;
        if (viamd::extract_int(i, arg) && i >= 0) rep->dipole.dipole_index = (uint32_t)i;
    } else if (str_eq(ident, STR_LIT("DipoleColor"))) {
        viamd::extract_vec4(rep->dipole.color, arg);
    } else if (str_eq(ident, STR_LIT("DipoleOffset"))) {
        viamd::extract_vec3(rep->dipole.offset, arg);
    } else if (str_eq(ident, STR_LIT("DipoleScale"))) {
        viamd::extract_dbl(rep->dipole.scale, arg);
    } else if (str_eq(ident, STR_LIT("DipoleRadius"))) {
        viamd::extract_flt(rep->dipole.radius, arg);
    } else if (str_eq(ident, STR_LIT("PeriodicImages"))) {
        viamd::extract_bool(rep->periodic_images, arg);
    } else if (str_eq(ident, STR_LIT("BoxColor"))) {
        viamd::extract_vec4(rep->box.color, arg);
    } else if (str_eq(ident, STR_LIT("BoxThickness"))) {
        viamd::extract_flt(rep->box.thickness, arg);
    } else if (str_eq(ident, STR_LIT("BoxAxes"))) {
        viamd::extract_bool(rep->box.axes, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureMoIdx")) || str_eq(ident, STR_LIT("ElectronicStructureOrbitalIdx"))) {
        viamd::extract_int(rep->electronic_structure.orbital_idx, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureNtoIdx")) || str_eq(ident, STR_LIT("ElectronicStructureExcitedStateIdx"))) {
        viamd::extract_int(rep->electronic_structure.excited_state_idx, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureNtoLambdaIdx"))) {
        viamd::extract_int(rep->electronic_structure.nto_lambda_idx, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureDensityPropertyPath"))) {
        // The id is a function of the path, so this resolves without the dataset being
        // loaded yet - and names the same property after a reload, which the index it
        // replaces did not.
        rep->electronic_structure.density_property_key = md_attributes_id_from_path(arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureDensityPropertyIdx"))) {
        // DEPRECATED. A position in a list that only existed while a particular file was
        // open; there is nothing to resolve it against. Workspaces written before the
        // path was stored fall back to the first property available.
    } else if (str_eq(ident, STR_LIT("ElectronicStructureRes"))) {
        int res;
        if (viamd::extract_int(res, arg)) rep->electronic_structure.resolution = (VolumeResolution)res;
    } else if (str_eq(ident, STR_LIT("ElectronicStructureType"))) {
        int type;
        if (viamd::extract_int(type, arg)) electronic_structure_set_legacy_type(&rep->electronic_structure, type);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureSource"))) {
        int source;
        if (viamd::extract_int(source, arg)) {
            rep->electronic_structure.source = (ElectronicStructureSource)source;
            electronic_structure_set_source_defaults(&rep->electronic_structure);
        }
    } else if (str_eq(ident, STR_LIT("ElectronicStructureField"))) {
        int field;
        if (viamd::extract_int(field, arg)) electronic_structure_set_legacy_field(&rep->electronic_structure, field);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureUseMagnitude"))) {
        viamd::extract_bool(rep->electronic_structure.use_magnitude, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureSpin"))) {
        int spin;
        if (viamd::extract_int(spin, arg)) rep->electronic_structure.spin = (ElectronicStructureSpin)spin;
    } else if (str_eq(ident, STR_LIT("ElectronicStructureNtoComponent"))) {
        int component;
        if (viamd::extract_int(component, arg)) rep->electronic_structure.nto_component = (ElectronicStructureNtoComponent)component;
    } else if (str_eq(ident, STR_LIT("ElectronicStructureTransitionDensityComponent"))) {
        int component;
        if (viamd::extract_int(component, arg)) rep->electronic_structure.transition_density_component = (ElectronicStructureTransitionDensityComponent)component;
    } else if (str_eq(ident, STR_LIT("ElectronicStructureIso"))) {
        viamd::extract_dbl(rep->electronic_structure.iso_value, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureColoring"))) {
        int c;
        if (viamd::extract_int(c, arg) && c >= 0 && c < (int)SurfaceColoring::Count) rep->electronic_structure.coloring = (SurfaceColoring)c;
    } else if (str_eq(ident, STR_LIT("ElectronicStructureFieldKind"))) {
        int k;
        if (viamd::extract_int(k, arg) && k >= 0 && k < (int)SurfaceFieldKind::Count) rep->electronic_structure.field_kind = (SurfaceFieldKind)k;
    } else if (str_eq(ident, STR_LIT("ElectronicStructureFieldColormap"))) {
        viamd::extract_int(rep->electronic_structure.field_map.colormap, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureFieldRange"))) {
        float r[2];
        if (viamd::extract_flt_vec(r, 2, arg)) { rep->electronic_structure.field_map.range_beg = r[0]; rep->electronic_structure.field_map.range_end = r[1]; }
    } else if (str_eq(ident, STR_LIT("ElectronicStructureFieldSymmetric"))) {
        viamd::extract_bool(rep->electronic_structure.field_map.symmetric, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureFieldAutoRange"))) {
        viamd::extract_bool(rep->electronic_structure.field_map.auto_range, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureFieldLegend"))) {
        viamd::extract_bool(rep->electronic_structure.field_map.show_legend, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureFieldIsolinesShow"))) {
        viamd::extract_bool(rep->electronic_structure.field_isolines.show, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureFieldIsolineWidth"))) {
        viamd::extract_flt(rep->electronic_structure.field_isolines.width, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureFieldIsolineCount"))) {
        int n;
        if (viamd::extract_int(n, arg)) {
            rep->electronic_structure.field_isolines.count = CLAMP(n, 0, ISOLINES_MAX);
        }
    } else if (str_eq(ident, STR_LIT("ElectronicStructureTintPos"))) {
        viamd::extract_vec4(rep->electronic_structure.tint_psi_pos, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureTintNeg"))) {
        viamd::extract_vec4(rep->electronic_structure.tint_psi_neg, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureTintDen"))) {
        viamd::extract_vec4(rep->electronic_structure.tint_den, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureColPos"))) {
        viamd::extract_vec4(rep->electronic_structure.col_psi_pos, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureColNeg"))) {
        viamd::extract_vec4(rep->electronic_structure.col_psi_neg, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureColDen"))) {
        viamd::extract_vec4(rep->electronic_structure.col_den, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureColAtt"))) {
        viamd::extract_vec4(rep->electronic_structure.col_att, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureColDet"))) {
        viamd::extract_vec4(rep->electronic_structure.col_det, arg);
    } else if (str_eq(ident, STR_LIT("ElectronicStructureDensityPropertyIsoCount"))) {
        int n;
        if (viamd::extract_int(n, arg)) {
            rep->electronic_structure.density_property.num_isos = CLAMP(n, 0, (int)ARRAY_SIZE(rep->electronic_structure.density_property.values));
        }
    } else if (str_begins_with(ident, STR_LIT("ElectronicStructureFieldIsoline"))) {
        // ElectronicStructureFieldIsoline<i>: the value of line i, in the field's own unit
        const str_t num = str_substr(ident, STR_LIT("ElectronicStructureFieldIsoline").len, SIZE_MAX);
        int i = 0;
        bool ok = num.len > 0;
        for (size_t k = 0; k < num.len; ++k) {
            if (num.ptr[k] < '0' || num.ptr[k] > '9') { ok = false; break; }
            i = i * 10 + (num.ptr[k] - '0');
        }
        if (ok && i < ISOLINES_MAX) {
            viamd::extract_flt(rep->electronic_structure.field_isolines.values[i], arg);
        }
    } else {
        for (int i = 0; i < (int)ARRAY_SIZE(rep->electronic_structure.density_property.values); ++i) {
            char key[64];
            snprintf(key, sizeof(key), "ElectronicStructureDensityPropertyIsoValue%d", i);
            if (str_eq(ident, str_from_cstr(key))) {
                viamd::extract_dbl(rep->electronic_structure.density_property.values[i], arg);
                return true;
            }
            snprintf(key, sizeof(key), "ElectronicStructureDensityPropertyIsoColor%d", i);
            if (str_eq(ident, str_from_cstr(key))) {
                viamd::extract_vec4(rep->electronic_structure.density_property.colors[i], arg);
                return true;
            }
        }
        return false;
    }
    return true;
}

// The entries of one section, read before what they are for is made: the system it is of is one of them (System=)
struct SectionEntries {
    md_array(str_t) ident = 0;
    md_array(str_t) arg = 0;
    int system = -1;    // The id it names, -1 for none
};

static void read_section_entries(SectionEntries* out, viamd::deserialization_state_t& state, md_allocator_i* alloc) {
    str_t ident, arg;
    while (viamd::next_entry(ident, arg, state)) {
        if (str_eq(ident, STR_LIT("System"))) {
            viamd::extract_int(out->system, arg);
            continue;
        }
        md_array_push(out->ident, ident, alloc);
        md_array_push(out->arg, arg, alloc);
    }
}

static void deserialize_representation(SystemContext& sys, const SectionEntries& entries) {
    if (!representation_can_add(sys)) {
        // Beyond as many as a system has: its entries are passed over
        VIAMD_LOG_ERROR("Workspace: more representations than a system has (at most %zu); the rest are left out", REPRESENTATION_MAX);
        return;
    }
    Representation* rep = create_representation(sys);
    for (size_t i = 0; i < md_array_size(entries.ident); ++i) {
        const str_t ident = entries.ident[i];
        const str_t arg = entries.arg[i];
        if (str_eq(ident, STR_LIT("Group"))) {
            // By the id of a [RepresentationGroup] before it. One that names nothing leaves it at the top level.
            int group = 0;
            if (viamd::extract_int(group, arg) && group > 0 && find_representation_group(sys, (uint32_t)group)) {
                rep->node.parent = (uint32_t)group;
            }
        } else {
            representation_read_entry(rep, ident, arg);
        }
    }
    // Last in the group it says it is in: the sections are in the order of the window
    tree_normalize(sys);
}

void representation_write_entries(viamd::serialization_state_t& state, const SystemContext& sys, const Representation& rep) {
    const md_attributes_t* attributes = &sys.data.sys.attributes;

    viamd::write_str(state,  STR_LIT("Name"), str_from_cstr(rep.name));
    viamd::write_str(state,  STR_LIT("Filter"), str_from_cstr(rep.filter.text));
    viamd::write_bool(state, STR_LIT("Enabled"), rep.enabled);
    viamd::write_int(state,  STR_LIT("Type"), (int)rep.type);
    viamd::write_int(state,  STR_LIT("ColorMapping"), (int)rep.color_mapping);
    viamd::write_vec4(state, STR_LIT("BaseColor"), rep.base_color);
    viamd::write_flt(state,  STR_LIT("Saturation"), rep.saturation);
    viamd::write_vec4(state, STR_LIT("TintColor"), rep.tint_color);
    viamd::write_flt(state,  STR_LIT("TintScale"), rep.tint_scale);
    viamd::write_vec4(state, STR_LIT("SecondaryStructureColorUnknown"), rep.secondary_structure.color_unknown);
    viamd::write_vec4(state, STR_LIT("SecondaryStructureColorCoil"),    rep.secondary_structure.color_coil);
    viamd::write_vec4(state, STR_LIT("SecondaryStructureColorHelix"),   rep.secondary_structure.color_helix);
    viamd::write_vec4(state, STR_LIT("SecondaryStructureColorSheet"),   rep.secondary_structure.color_sheet);
    viamd::write_int(state,  STR_LIT("BondColor"), (int)rep.bond_color);
    viamd::write_flt(state,  STR_LIT("BondSharpness"), rep.bond_sharpness);
    viamd::write_vec4(state, STR_LIT("BondBaseColor"), rep.bond_base_color);
    viamd::write_bool(state, STR_LIT("BondOrders"), rep.bond_orders);
    if (rep.type == RepresentationType::Bonds) {
        viamd::write_bool(state, STR_LIT("BondsDashed"), rep.bonds_dashed);
    }
    viamd::write_bool(state, STR_LIT("ChainGaps"), rep.chain_gaps);
    viamd::write_bool(state, STR_LIT("ElasticBonds"), rep.elastic_bonds);
    viamd::write_bool(state, STR_LIT("NoncovalentBonds"), rep.noncovalent_bonds);
    viamd::write_bool(state, STR_LIT("TwoHeadSpheres"), rep.two_head_spheres);
    viamd::write_flt(state,  STR_LIT("NucleicBandRadius"), rep.nucleic_band_radius);
    viamd::write_flt(state,  STR_LIT("NucleicRingRadius"), rep.nucleic_ring_radius);

    viamd::write_vec4(state, STR_LIT("Param"), rep.scale);
    viamd::write_bool(state, STR_LIT("PeriodicImages"), rep.periodic_images);
    viamd::write_bool(state, STR_LIT("AutoUpdate"), rep.dynamic_evaluation);

    // Attributes by path: an id is a hash of it, and the path is what a reader can see
    if (const md_attribute_t* prop = md_attributes_get(attributes, rep.atom_attribute.key)) {
        viamd::write_str(state, STR_LIT("AtomicPropertyPath"), prop->path);
        viamd::write_int(state, STR_LIT("AtomicPropertyVariant"), rep.atom_attribute.variant_idx);
        const ColorScale& scale = rep.atom_attribute.scale;
        const float range[2] = { scale.range_beg, scale.range_end };
        viamd::write_int(state, STR_LIT("AtomicPropertyColormap"), scale.colormap);
        viamd::write_flt_vec(state, STR_LIT("AtomicPropertyRange"), range, 2);
        // After the range, which on its own reads as a range somebody set
        viamd::write_bool(state, STR_LIT("AtomicPropertyAutoRange"), scale.auto_range);
        viamd::write_bool(state, STR_LIT("AtomicPropertySymmetric"), scale.symmetric);
        viamd::write_bool(state, STR_LIT("AtomicPropertyLegend"), scale.show_legend);
    }

    if (rep.type == RepresentationType::DipoleMoment) {
        if (const md_attribute_t* dip = md_attributes_get(attributes, rep.dipole.dipole_key)) {
            viamd::write_str(state, STR_LIT("DipolePath"), dip->path);
        }
        viamd::write_int(state,  STR_LIT("DipoleIndex"), (int64_t)rep.dipole.dipole_index);
        viamd::write_vec4(state, STR_LIT("DipoleColor"), rep.dipole.color);
        viamd::write_vec3(state, STR_LIT("DipoleOffset"), rep.dipole.offset);
        viamd::write_dbl(state,  STR_LIT("DipoleScale"), rep.dipole.scale);
        viamd::write_flt(state,  STR_LIT("DipoleRadius"), rep.dipole.radius);
    }

    if (rep.type == RepresentationType::SimulationBox) {
        viamd::write_vec4(state, STR_LIT("BoxColor"), rep.box.color);
        viamd::write_flt(state,  STR_LIT("BoxThickness"), rep.box.thickness);
        viamd::write_bool(state, STR_LIT("BoxAxes"), rep.box.axes);
    }

    if (rep.type == RepresentationType::ElectronicStructure) {
        viamd::write_int(state,  STR_LIT("ElectronicStructureMoIdx"),    rep.electronic_structure.orbital_idx);
        viamd::write_int(state,  STR_LIT("ElectronicStructureNtoIdx"),   rep.electronic_structure.excited_state_idx);
        viamd::write_int(state,  STR_LIT("ElectronicStructureNtoLambdaIdx"), rep.electronic_structure.nto_lambda_idx);
        if (const md_attribute_t* prop = md_attributes_get(attributes, rep.electronic_structure.density_property_key)) {
            viamd::write_str(state, STR_LIT("ElectronicStructureDensityPropertyPath"), prop->path);
        }
        viamd::write_int(state,  STR_LIT("ElectronicStructureType"),     electronic_structure_legacy_type(rep.electronic_structure));
        viamd::write_int(state,  STR_LIT("ElectronicStructureSource"),   (int)rep.electronic_structure.source);
        viamd::write_int(state,  STR_LIT("ElectronicStructureField"),    (int)electronic_structure_legacy_field(rep.electronic_structure));
        viamd::write_bool(state, STR_LIT("ElectronicStructureUseMagnitude"), rep.electronic_structure.use_magnitude);
        viamd::write_int(state,  STR_LIT("ElectronicStructureSpin"),     (int)rep.electronic_structure.spin);
        viamd::write_int(state,  STR_LIT("ElectronicStructureNtoComponent"), (int)rep.electronic_structure.nto_component);
        viamd::write_int(state,  STR_LIT("ElectronicStructureTransitionDensityComponent"), (int)rep.electronic_structure.transition_density_component);
        viamd::write_int(state,  STR_LIT("ElectronicStructureRes"), (int)rep.electronic_structure.resolution);
        viamd::write_dbl(state,  STR_LIT("ElectronicStructureIso"),      rep.electronic_structure.iso_value);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureColPos"),   rep.electronic_structure.col_psi_pos);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureColNeg"),   rep.electronic_structure.col_psi_neg);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureColDen"),   rep.electronic_structure.col_den);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureColAtt"),   rep.electronic_structure.col_att);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureColDet"),   rep.electronic_structure.col_det);
        viamd::write_int(state,  STR_LIT("ElectronicStructureColoring"), (int)rep.electronic_structure.coloring);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureTintPos"),  rep.electronic_structure.tint_psi_pos);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureTintNeg"),  rep.electronic_structure.tint_psi_neg);
        viamd::write_vec4(state, STR_LIT("ElectronicStructureTintDen"),  rep.electronic_structure.tint_den);
        if (rep.electronic_structure.coloring == SurfaceColoring::Field) {
            const float range[2] = { rep.electronic_structure.field_map.range_beg, rep.electronic_structure.field_map.range_end };
            viamd::write_int(state,      STR_LIT("ElectronicStructureFieldKind"),      (int)rep.electronic_structure.field_kind);
            viamd::write_int(state,      STR_LIT("ElectronicStructureFieldColormap"),  rep.electronic_structure.field_map.colormap);
            viamd::write_flt_vec(state,  STR_LIT("ElectronicStructureFieldRange"),     range, 2);
            viamd::write_bool(state,     STR_LIT("ElectronicStructureFieldSymmetric"), rep.electronic_structure.field_map.symmetric);
            viamd::write_bool(state,     STR_LIT("ElectronicStructureFieldAutoRange"), rep.electronic_structure.field_map.auto_range);
            viamd::write_bool(state,     STR_LIT("ElectronicStructureFieldLegend"),    rep.electronic_structure.field_map.show_legend);
            const ColorScaleIsolines& lines = rep.electronic_structure.field_isolines;
            viamd::write_bool(state,     STR_LIT("ElectronicStructureFieldIsolinesShow"),  lines.show);
            viamd::write_flt(state,      STR_LIT("ElectronicStructureFieldIsolineWidth"),  lines.width);
            viamd::write_int(state,      STR_LIT("ElectronicStructureFieldIsolineCount"),  lines.count);
            for (int j = 0; j < lines.count; ++j) {
                char key[64];
                snprintf(key, sizeof(key), "ElectronicStructureFieldIsoline%d", j);
                viamd::write_flt(state, str_from_cstr(key), lines.values[j]);
            }
        }
        viamd::write_int(state,  STR_LIT("ElectronicStructureDensityPropertyIsoCount"), rep.electronic_structure.density_property.num_isos);
        for (int j = 0; j < rep.electronic_structure.density_property.num_isos; ++j) {
            char key[64];
            snprintf(key, sizeof(key), "ElectronicStructureDensityPropertyIsoValue%d", j);
            viamd::write_dbl(state, str_from_cstr(key), rep.electronic_structure.density_property.values[j]);
            snprintf(key, sizeof(key), "ElectronicStructureDensityPropertyIsoColor%d", j);
            viamd::write_vec4(state, str_from_cstr(key), rep.electronic_structure.density_property.colors[j]);
        }
    }
}

static void serialize_representation(viamd::serialization_state_t& state, const SystemContext& sys, const Representation& rep) {
    viamd::write_section_header(state, STR_LIT("Representation"));
    viamd::write_int(state, STR_LIT("System"), (int)sys.id);
    if (find_representation_group(sys, rep.node.parent)) {
        viamd::write_int(state, STR_LIT("Group"), rep.node.parent);
    }
    representation_write_entries(state, sys, rep);
}

static void deserialize_group(SystemContext& sys, const SectionEntries& entries) {
    ApplicationState* data = sys.app;
    RepresentationGroup group = {};
    int id = 0;
    int parent = 0;
    for (size_t i = 0; i < md_array_size(entries.ident); ++i) {
        const str_t ident = entries.ident[i];
        const str_t arg = entries.arg[i];
        if (str_eq(ident, STR_LIT("Id"))) {
            viamd::extract_int(id, arg);
        } else if (str_eq(ident, STR_LIT("Group"))) {
            viamd::extract_int(parent, arg);
        } else if (str_eq(ident, STR_LIT("Name"))) {
            viamd::extract_to_char_buf(group.name, sizeof(group.name), arg);
        } else if (str_eq(ident, STR_LIT("Filter"))) {
            viamd::extract_to_char_buf(group.filter.text, sizeof(group.filter.text), arg);
        } else if (str_eq(ident, STR_LIT("Enabled"))) {
            viamd::extract_bool(group.enabled, arg);
        }
    }

    // The id it was saved with is what the sections after it refer to it by, so it is kept - unless another
    // group already has it, which only a workspace edited by hand can make happen
    if (id > 0 && !find_representation_group(sys, (uint32_t)id)) {
        group.node.id = (uint32_t)id;
        sys.representation.next_group_id = MAX(sys.representation.next_group_id, (uint32_t)id + 1);
    } else {
        group.node.id = sys.representation.next_group_id++;
    }
    // A group comes before what is within it, so the group it is in has been read already
    if (parent > 0 && find_representation_group(sys, (uint32_t)parent)) {
        group.node.parent = (uint32_t)parent;
    }
    group.node.order = INT32_MAX;
    group.node.group = true;
    group.filt_is_dirty = true;
    md_array_push(sys.representation.groups, group, data->allocator.persistent);
    tree_normalize(sys);
}

static void serialize_group(viamd::serialization_state_t& state, const SystemContext& sys, const RepresentationGroup& group) {
    viamd::write_section_header(state, STR_LIT("RepresentationGroup"));
    viamd::write_int(state, STR_LIT("System"), (int)sys.id);
    viamd::write_int(state, STR_LIT("Id"), group.node.id);
    if (find_representation_group(sys, group.node.parent)) {
        viamd::write_int(state, STR_LIT("Group"), group.node.parent);
    }
    viamd::write_str(state,  STR_LIT("Name"), str_from_cstr(group.name));
    viamd::write_str(state,  STR_LIT("Filter"), str_from_cstr(group.filter.text));
    viamd::write_bool(state, STR_LIT("Enabled"), group.enabled);
}

// Depth first, the order of the window: a group is written before what is within it, which is what reading
// it back relies on
static void serialize_level(viamd::serialization_state_t& state, const SystemContext& sys, const TreeNodes& t, uint32_t parent, size_t depth) {
    if (depth > t.count) return;
    md_temp_scope_t temp = md_temp_begin();
    int* kids = (int*)md_temp_alloc(temp, sizeof(int) * MAX(t.count, (size_t)1));
    const size_t n = representation_tree_children(kids, t.count, t.nodes, t.count, parent);
    for (size_t i = 0; i < n; ++i) {
        const size_t k = (size_t)kids[i];
        if (k < t.num_reps) {
            serialize_representation(state, sys, sys.representation.reps[k]);
        } else {
            const RepresentationGroup& group = sys.representation.groups[k - t.num_reps];
            serialize_group(state, sys, group);
            serialize_level(state, sys, t, group.node.id, depth + 1);
        }
    }
    md_temp_end(temp);
}

void representations_serialize(viamd::serialization_state_t& state, const SystemContext& sys) {
    md_temp_scope_t temp = md_temp_begin();
    const TreeNodes t = tree_nodes(sys, md_temp_allocator(temp));
    serialize_level(state, sys, t, 0, 0);
    md_temp_end(temp);
}

bool representations_deserialize(str_t section, viamd::deserialization_state_t& state, WorkspaceSystemFn system_of, void* user) {
    const bool rep = str_eq(section, STR_LIT("Representation"));
    const bool group = str_eq(section, STR_LIT("RepresentationGroup"));
    if (!rep && !group) return false;

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    SectionEntries entries = {};
    read_section_entries(&entries, state, md_temp_allocator(temp));
    // Of a system the workspace does not load: left out
    SystemContext* sys = system_of(user, entries.system);
    if (!sys) return true;
    if (rep) deserialize_representation(*sys, entries);
    else     deserialize_group(*sys, entries);
    return true;
}

// ---------------------------------------------------------------------------------------------------------
// The properties of an electronic structure representation
// ---------------------------------------------------------------------------------------------------------

// How the isosurfaces of an electronic structure representation are coloured, as rows of its Color section: the
// colouring, and for a field the field, its colour scale and its isolines. The atoms' colours are not offered for a
// density property; a field is offered whenever the dataset has what one is made from.
static bool draw_surface_coloring(ElectronicStructureRepresentation& es, const md_system_t& sys, bool allow_atom_colors, bool mesh_mode) {
    bool update_rep = false;

    bool any_field = false;
    for (int k = 0; k < (int)SurfaceFieldKind::Count; ++k) {
        any_field |= surface_field_available((SurfaceFieldKind)k, sys);
    }

    prop::row("coloring");
    if (ImGui::BeginCombo("##coloring", surface_coloring_str[(int)es.coloring])) {
        for (int i = 0; i < (int)SurfaceColoring::Count; ++i) {
            const SurfaceColoring c = (SurfaceColoring)i;
            const bool enabled = (c == SurfaceColoring::AtomColors) ? allow_atom_colors :
                                 (c == SurfaceColoring::Field)      ? any_field : true;
            ImGui::BeginDisabled(!enabled);
            if (ImGui::Selectable(surface_coloring_str[i], es.coloring == c)) {
                es.coloring = c;
                // A field this dataset can make, if the one the representation names is not
                if (c == SurfaceColoring::Field && !surface_field_available(es.field_kind, sys)) {
                    const SurfaceFieldKind first = surface_field_first_available(sys);
                    surface_field_scale_for_kind(&es.field_map, es.field_kind, first);
                    es.field_kind = first;
                }
                update_rep = true;
            }
            ImGui::EndDisabled();
            if (!enabled && c == SurfaceColoring::Field) {
                ImGui::SetItemTooltip("Nothing in this dataset to make a field from:\n"
                                      "the electrostatic potential needs a basis, a ground state density and the QM atoms' nuclear charges,\n"
                                      "the local ionization energy a basis and canonical orbitals with their occupations and energies,\n"
                                      "the local electron attachment energy the same and at least one virtual orbital below zero,\n"
                                      "the embedding potential classical charges (atom/charge) on atoms outside the QM region");
            }
        }
        ImGui::EndCombo();
    }

    if (es.coloring == SurfaceColoring::AtomColors) {
        const double min_power = 2.0;
        const double max_power = 20.0;
        prop::row("splat power", "How sharply the colours of the atoms fall off around them: the power of the gaussians they are splatted with");
        update_rep |= ImGui::SliderScalar("##splat_power", ImGuiDataType_Double, &es.gaussian_splatting_power, &min_power, &max_power, "%.2f");
    }

    if (es.coloring != SurfaceColoring::Field) {
        return update_rep;
    }

    // Where the field stands, on its label: how much of it has been evaluated, and where
    const SurfaceFieldVolume& vol = es.field_vol;
    char status[160] = "";
    if (!mesh_mode && surface_field_busy(vol)) {
        snprintf(status, sizeof(status), "evaluating, %zu of %zu voxels (%dx%dx%d) so far", vol.num_evaluated, vol.num_voxels,
                 vol.grid.dim[0], vol.grid.dim[1], vol.grid.dim[2]);
    } else if (!mesh_mode && vol.num_surface_samples > 0) {
        snprintf(status, sizeof(status), "%zu of %zu voxels (%dx%dx%d) evaluated on the %s", vol.num_evaluated, vol.num_voxels,
                 vol.grid.dim[0], vol.grid.dim[1], vol.grid.dim[2], vol.on_gpu ? "GPU" : "CPU");
    }
    prop::row("field", status);
    if (ImGui::BeginCombo("##field", surface_field_kind_info(es.field_kind).label)) {
        for (int k = 0; k < (int)SurfaceFieldKind::Count; ++k) {
            const SurfaceFieldKind kind = (SurfaceFieldKind)k;
            ImGui::BeginDisabled(!surface_field_available(kind, sys));
            if (ImGui::Selectable(surface_field_kind_info(kind).label, es.field_kind == kind)) {
                if (es.field_kind != kind) surface_field_scale_for_kind(&es.field_map, es.field_kind, kind);
                es.field_kind = kind;
                update_rep = true;
            }
            ImGui::EndDisabled();
        }
        ImGui::EndCombo();
    }
    if (es.field_kind == SurfaceFieldKind::EmbeddingPotential) {
        // What this dataset's potential is made of, since a polarizable run may or may not have kept its
        // induced dipoles. Counted only while the tooltip shows: it reads a column over every atom.
        if (ImGui::BeginItemTooltip()) {
            ImGui::TextUnformatted("The electrostatic potential of the environment on the surface: its charges, the dipoles and\n"
                                   "quadrupoles where the potential has them, and for a polarizable embedding the dipoles the\n"
                                   "calculation induced, where it kept them - the potential the QM region was embedded in.\n"
                                   "The QM region is left out, whatever charges a topology gives its atoms.");
            // What this dataset's potential is made of: whatever it carries, as it carries it
            const SurfaceFieldEmbeddingContent c = surface_field_embedding_content(sys);
            ImGui::Separator();
            ImGui::TextUnformatted("This dataset:");
            if (c.charges)     ImGui::BulletText("charges on %zu sites", c.charges);
            if (c.dipoles)     ImGui::BulletText("permanent dipoles on %zu sites", c.dipoles);
            if (c.quadrupoles) ImGui::BulletText("quadrupoles on %zu sites", c.quadrupoles);
            if (c.induced)     ImGui::BulletText("induced dipoles on %zu polarizable sites", c.induced);
            if (c.induced) {
                ImGui::TextUnformatted("The total potential the QM region was embedded in.");
            } else if (c.polarizable) {
                ImGui::Text("The PERMANENT part only: %zu sites are polarizable, but the file keeps no induced\n"
                            "dipoles, so their part of what the QM region was embedded in is missing.", c.polarizable);
            }
            ImGui::EndTooltip();
        }
    } else if (es.field_kind == SurfaceFieldKind::ElectrostaticPotential) {
        ImGui::SetItemTooltip("The electrostatic potential of the QM region on the surface: its nuclei and its ground state\n"
                              "electron density, exact over the basis. Negative (red) where the electrons dominate - lone pairs,\n"
                              "pi systems, where an electrophile is drawn - and positive (blue) where the nuclei are less screened.\n"
                              "Mapped onto an isodensity surface (0.001-0.002 au) it is the molecular electrostatic potential.\n"
                              "An embedding's environment is not part of it.");
    } else if (es.field_kind == SurfaceFieldKind::AverageLocalIonizationEnergy) {
        ImGui::SetItemTooltip("The average local ionization energy of the QM region on the surface: the orbital energies of its\n"
                              "occupied orbitals, weighted by how much each contributes to the electron density there,\n"
                              "sum n |e| |phi|^2 / sum n |phi|^2. Lowest (red) where the least tightly held electrons are - the\n"
                              "sites an electrophile or a radical attacks - and never below the HOMO's |e|. Meant for an\n"
                              "isodensity surface (0.001-0.002 au). Kohn-Sham orbital energies give lower values than Hartree-Fock.");
    } else if (es.field_kind == SurfaceFieldKind::LocalElectronAttachmentEnergy) {
        ImGui::SetItemTooltip("The local electron attachment energy of the QM region on the surface: the energies of its virtual\n"
                              "orbitals below zero, the ones that bind an added electron, weighted by their densities against the\n"
                              "electron density. Never positive; lowest where a nucleophile attacks or a halogen bond forms (the\n"
                              "sigma hole). Defined for hybrid DFT orbitals and read on the 0.004 au isodensity surface; it depends\n"
                              "strongly on the method and basis, and is 0 everywhere with Hartree-Fock virtuals.");
    }

    // The colours are applied per pixel as the surface is shaded, so a change to the scale needs
    // nothing re-evaluated: it is not an update of the representation
    // With mesh rendering the field is evaluated at the vertices of the meshes, not on a volume
    color_scale_draw_controls(&es.field_map, mesh_mode ? electronic_structure_mesh_field_span(es) : surface_field_span(vol),
                              surface_field_kind_info(es.field_kind).unit, "on surface", "its 1st to 99th percentile on the surface");
    // Drawn as the surface is shaded, like the colours: nothing to re-evaluate
    color_scale_draw_isoline_controls(&es.field_isolines, es.field_map, surface_field_kind_info(es.field_kind).unit);

    return update_rep;
}


// The properties of an electronic structure representation: what the volume is of (Source), the surfaces drawn from
// it (Surface) and how they are coloured (Color)
static bool draw_representations_window_electronic_structure(SystemContext& sys, Representation& rep) {
    ApplicationState* state = sys.app;
    ElectronicStructureRepresentation& es = rep.electronic_structure;

    // The key is an opaque hash now, so it is no longer worth showing beside the label: a property
    // with no label is named by its path's last segment already, which is the readable thing.
    auto density_property_combo_label = [](char* buf, size_t buf_cap, const DensityProperty* prop) {
        ASSERT(buf);
        ASSERT(buf_cap > 0);
        if (!prop || str_empty(prop->label)) {
            snprintf(buf, buf_cap, "None");
            return;
        }
        snprintf(buf, buf_cap, "%.*s", (int)prop->label.len, prop->label.ptr);
    };

    bool update_rep = false;
    const bool mesh_mode = state->settings.isosurface_rendering == IsosurfaceRendering_Mesh;
    const ElectronicStructureSourceFlags source_mask = es_source_mask(sys.data.sys);

    if (prop::section("Source")) {
        prop::row("source");
        if (ImGui::BeginCombo("##source", electronic_structure_source_str[(int)es.source])) {
            for (int n = 0; n < (int)ElectronicStructureSource::Count; n++) {
                ElectronicStructureSource source = (ElectronicStructureSource)n;
                bool is_selected = (es.source == source);
                bool disabled = !electronic_structure_source_supported(source_mask, source);

                if (disabled) ImGui::PushDisabled();
                if (ImGui::Selectable(electronic_structure_source_str[n], is_selected)) {
                    es.source = source;
                    electronic_structure_set_source_defaults(&es);
                    update_rep = true;
                }
                if (disabled) ImGui::PopDisabled();

                if (is_selected) {
                    ImGui::SetItemDefaultFocus();
                }
            }
            ImGui::EndCombo();
        }

        if (electronic_structure_is_density_property(es)) {
            // Gathered here rather than kept in a list somewhere: it is a handful of entries borrowed
            // from the attribute table, and a stored copy would be one more thing to invalidate on load.
            DensityProperty density_props[64];
            const int num_density_props = (int)MIN(density_properties_gather(density_props, ARRAY_SIZE(density_props), sys.data.sys), ARRAY_SIZE(density_props));
            if (num_density_props > 0) {
                // A key that names nothing in the current list - a fresh representation, or a workspace
                // saved against a file whose properties differ - falls back to the first one rather
                // than drawing nothing.
                const DensityProperty* selected = nullptr;
                for (int n = 0; n < num_density_props; ++n) {
                    if (density_props[n].key == es.density_property_key) {
                        selected = &density_props[n];
                        break;
                    }
                }
                if (!selected) {
                    selected = &density_props[0];
                    es.density_property_key = selected->key;
                    update_rep = true;
                }

                char preview[256];
                density_property_combo_label(preview, sizeof(preview), selected);
                prop::row("property");
                if (ImGui::BeginCombo("##property", preview)) {
                    for (int n = 0; n < num_density_props; ++n) {
                        char label[256];
                        density_property_combo_label(label, sizeof(label), &density_props[n]);
                        const bool is_selected = (es.density_property_key == density_props[n].key);
                        if (ImGui::Selectable(label, is_selected)) {
                            if (es.density_property_key != density_props[n].key) {
                                es.density_property_key = density_props[n].key;
                                update_rep = true;
                            }
                        }

                        if (is_selected) {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }
            } else {
                prop::hint("no density properties in this dataset");
            }
        }

        if (electronic_structure_uses_spin(es)) {
            // Asked of the table: a second, DISTINCT set of coefficients. This used to be true for a
            // restricted file too, because beta is published there as an alias of alpha and the fan-in
            // only checked that the column existed - so the selector appeared and toggled between two
            // names for one datum.
            const bool has_beta = es_has_distinct_beta_orbitals(sys.data.sys);
            if (es.source == ElectronicStructureSource::MolecularOrbital) {
                if (es.spin != ElectronicStructureSpin::Alpha && es.spin != ElectronicStructureSpin::Beta) {
                    es.spin = ElectronicStructureSpin::Alpha;
                    update_rep = true;
                }
                if (has_beta) {
                    const char* spin_options[] = { "Alpha", "Beta" };
                    int spin_idx = es.spin == ElectronicStructureSpin::Beta ? 1 : 0;
                    prop::row("spin");
                    if (ImGui::Combo("##spin", &spin_idx, spin_options, IM_ARRAYSIZE(spin_options))) {
                        es.spin = spin_idx == 1 ? ElectronicStructureSpin::Beta : ElectronicStructureSpin::Alpha;
                        update_rep = true;
                    }
                }
            }
            else if (es.source == ElectronicStructureSource::ElectronDensity) {
                if (has_beta) {
                    const ElectronicStructureSpin spin_options[] = {
                        ElectronicStructureSpin::Total,
                        ElectronicStructureSpin::Alpha,
                        ElectronicStructureSpin::Beta,
                        ElectronicStructureSpin::Difference,
                    };
                    const char* spin_labels[] = { "Total", "Alpha", "Beta", "Difference" };
                    int spin_idx = 0;
                    for (int i = 0; i < (int)IM_ARRAYSIZE(spin_options); ++i) {
                        if (es.spin == spin_options[i]) {
                            spin_idx = i;
                            break;
                        }
                    }
                    if (es.spin != spin_options[spin_idx]) {
                        es.spin = spin_options[spin_idx];
                        update_rep = true;
                    }
                    prop::row("spin");
                    if (ImGui::Combo("##spin", &spin_idx, spin_labels, IM_ARRAYSIZE(spin_labels))) {
                        es.spin = spin_options[spin_idx];
                        if (es.spin != ElectronicStructureSpin::Difference) {
                            es.use_magnitude = false;
                        }
                        update_rep = true;
                    }
                } else if (es.spin != ElectronicStructureSpin::Total) {
                    es.spin = ElectronicStructureSpin::Total;
                    es.use_magnitude = false;
                    update_rep = true;
                }
            }
        }

        if (es.source == ElectronicStructureSource::NaturalTransitionOrbital) {
            prop::row("component");
            if (ImGui::Combo("##component", (int*)&es.nto_component, electronic_structure_nto_component_str, IM_ARRAYSIZE(electronic_structure_nto_component_str))) {
                update_rep = true;
            }
        }

        if (es.source == ElectronicStructureSource::TransitionDensity) {
            prop::row("component");
            if (ImGui::Combo("##component", (int*)&es.transition_density_component, electronic_structure_transition_density_component_str, IM_ARRAYSIZE(electronic_structure_transition_density_component_str))) {
                update_rep = true;
            }
        }

        if (electronic_structure_uses_orbital_idx(es)) {
            // The (homo)/(lumo) markers sit at different indices for the two spins in an unrestricted
            // calculation, and the rendered volume already follows es.spin - so the frontier has to
            // follow it too, or the entry marked (homo) under spin=Beta is alpha's HOMO index.
            //
            // Formatted here rather than read out of a precomputed list: it is an index, an energy and
            // two markers, and a copy of it would be one more thing to keep in step with the table.
            const bool use_beta = (es.spin == ElectronicStructureSpin::Beta) && es_has_distinct_beta_orbitals(sys.data.sys);

            size_t num_orbitals = 0;
            OrbitalFrontier frontier = {};
            if (es_orbital_extent(sys.data.sys, &num_orbitals, nullptr) && num_orbitals > 0) {
                es_orbital_frontier(&frontier, sys.data.sys, use_beta ? es_path::beta_occupation : es_path::alpha_occupation);

                // Energies alongside the index, so the list can be read as an orbital diagram. A column
                // that does not span the orbital set is not one to index with, so it is dropped whole
                // rather than bounds checked per row - the list then falls back to bare indices.
                size_t num_energies = 0;
                const double* energy = es_orbital_energies(&num_energies, sys.data.sys, use_beta ? es_path::beta_energy : es_path::alpha_energy);
                if (num_energies != num_orbitals) energy = nullptr;

                // The UI font is monospaced (Dejavu Sans Mono), so a fixed field width per column is all
                // the alignment this needs: right aligned index, then a fixed number of decimals, which
                // puts every decimal point of the list on one column. Widths are passed in rather than
                // baked into the format so the preview - a single line, nothing to align against - can
                // ask for the same text unpadded.
                int idx_width = 1;
                for (size_t m = num_orbitals; m >= 10; m /= 10) idx_width++;

                auto orbital_label = [&](char* buf, size_t cap, int n, int idx_w, int ene_w) {
                    const char* mark = (n == frontier.homo_idx) ? " (homo)" : (n == frontier.lumo_idx) ? " (lumo)" : "";
                    if (energy) {
                        snprintf(buf, cap, "%*i  %*.4f%s", idx_w, n + 1, ene_w, energy[n], mark);
                    } else {
                        snprintf(buf, cap, "%*i%s", idx_w, n + 1, mark);
                    }
                };

                es.orbital_idx = CLAMP(es.orbital_idx, 0, (int)num_orbitals - 1);

                char preview[48];
                orbital_label(preview, sizeof(preview), es.orbital_idx, 0, 0);
                prop::row("orbital", energy ? "The orbitals by index and energy (Ha), highest first" : nullptr);
                if (ImGui::BeginCombo("##orbital", preview)) {
                    // Highest index first, so the lowest lying orbital sits at the bottom - the
                    // conventional orbital diagram, and the same order the Molecular Orbitals table
                    // defaults to. Canonical MOs come out of the SCF in ascending energy, so descending
                    // index IS descending energy; sorting on the energies instead would only differ for
                    // a non-aufbau set, and would then disagree with that table.
                    for (int n = (int)num_orbitals - 1; n >= 0; n--) {
                        bool is_selected = (es.orbital_idx == n);
                        char label[48];
                        orbital_label(label, sizeof(label), n, idx_width, 10);
                        if (ImGui::Selectable(label, is_selected)) {
                            if (es.orbital_idx != n) {
                                update_rep = true;
                            }
                            es.orbital_idx = n;
                        }

                        if (is_selected) {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }
            }
        }

        if (electronic_structure_uses_excited_state_idx(es)) {
            const int num_excited_states = (int)es_excited_state_count(sys.data.sys);
            if (num_excited_states > 0) {
                es.excited_state_idx = CLAMP(es.excited_state_idx, 0, num_excited_states - 1);
                char preview[16];
                snprintf(preview, sizeof(preview), "%i", es.excited_state_idx + 1);
                prop::row("state", "The excited state");
                if (ImGui::BeginCombo("##state", preview)) {
                    for (int n = 0; n < num_excited_states; n++) {
                        const bool is_selected = (es.excited_state_idx == n);
                        char label[16];
                        snprintf(label, sizeof(label), "%i", n + 1);
                        if (ImGui::Selectable(label, is_selected)) {
                            if (es.excited_state_idx != n) {
                                update_rep = true;
                            }
                            es.excited_state_idx = n;
                        }

                        if (is_selected) {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }
                if (electronic_structure_uses_nto_lambda_idx(es)) {
                    // The same 1e-3 the rest of the tree stops at, and the same place the zero padded
                    // tail of the ragged lambda axis stops it.
                    const double LAMBDA_CUTOFF = 1.0e-3;
                    double lambdas[32];
                    const bool nto_ready = es_nto_ready(state, (size_t)es.excited_state_idx);
                    const int num_lambdas = nto_ready ? (int)es_nto_lambdas(lambdas, ARRAY_SIZE(lambdas), sys.data.sys,
                                                                            (size_t)es.excited_state_idx, LAMBDA_CUTOFF) : 0;
                    if (!nto_ready) {
                        prop::row("\xce\xbb", "The natural transition orbital pair, by its weight");
                        ImGui::TextDisabled("computing...");
                    }
                    if (num_lambdas > 0) {
                        es.nto_lambda_idx = CLAMP(es.nto_lambda_idx, 0, num_lambdas - 1);

                        auto lambda_label = [&lambdas](char* buf, size_t cap, int n) {
                            // The bytes rather than a u8 literal: u8"" is const char8_t* in C++20 and
                            // a format string has to be a plain char array. \xce\xbb is U+03BB, lambda.
                            snprintf(buf, cap, "\xce\xbb[%i] (%.3f)", n + 1, lambdas[n]);
                        };

                        char preview_lambda[32];
                        lambda_label(preview_lambda, sizeof(preview_lambda), es.nto_lambda_idx);
                        prop::row("\xce\xbb", "The natural transition orbital pair, by its weight");
                        if (ImGui::BeginCombo("##lambda", preview_lambda)) {
                            for (int n = 0; n < num_lambdas; n++) {
                                const bool is_selected = (es.nto_lambda_idx == n);
                                char label[32];
                                lambda_label(label, sizeof(label), n);
                                if (ImGui::Selectable(label, is_selected)) {
                                    if (es.nto_lambda_idx != n) {
                                        update_rep = true;
                                    }
                                    es.nto_lambda_idx = n;
                                }

                                if (is_selected) {
                                    ImGui::SetItemDefaultFocus();
                                }
                            }
                            ImGui::EndCombo();
                        }
                    }
                }
            }
        }

        if (electronic_structure_uses_magnitude_toggle(es)) {
            const bool density = es.source == ElectronicStructureSource::ElectronDensity;
            prop::row("magnitude", density ? (const char*)u8"The surfaces of |ρ| instead of ρ" : (const char*)u8"The surfaces of |Ψ| instead of Ψ");
            update_rep |= ImGui::Checkbox("##magnitude", &es.use_magnitude);
        }
        prop::section_end();
    }

    const double min_tau = 0.0;
    const double max_tau = 1.0;
    const double iso_min = 1.0e-2;
    const double iso_max = 5.0;

    if (prop::section("Surface")) {
        if (electronic_structure_is_density_property(es)) {
            // Its isosurfaces as a list: the value and the colour of each, at least one of them
            const int cap = (int)ARRAY_SIZE(es.density_property.values);
            const int num_isos = CLAMP(es.density_property.num_isos, 1, cap);
            if (num_isos != es.density_property.num_isos) {
                es.density_property.num_isos = num_isos;
                update_rep = true;
            }
            const ImGuiStyle& style = ImGui::GetStyle();
            const float btn_w = ImGui::GetFrameHeight();
            const float gap   = style.ItemInnerSpacing.x;
            int remove = -1;
            for (int i = 0; i < es.density_property.num_isos; ++i) {
                ImGui::PushID(i);
                prop::row(i == 0 ? "isovalues" : "");
                ImGui::SetNextItemWidth(MAX(prop::control_width() - 2.0f * (btn_w + gap), 1.0f));
                update_rep |= ImGui::DragScalar("##value", ImGuiDataType_Double, &es.density_property.values[i], 0.001f, NULL, NULL, "%.6g");
                ImGui::SameLine(0.0f, gap);
                update_rep |= ImGui::ColorEdit4("##color", es.density_property.colors[i].elem, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel | ImGuiColorEditFlags_AlphaPreviewHalf);
                ImGui::SameLine(0.0f, gap);
                ImGui::BeginDisabled(es.density_property.num_isos <= 1);
                if (ImGui::DeleteButton(ICON_FA_XMARK, ImVec2(btn_w, 0))) remove = i;
                ImGui::EndDisabled();
                ImGui::SetItemTooltip("Remove the isosurface");
                ImGui::PopID();
            }
            if (remove >= 0) {
                for (int i = remove; i + 1 < es.density_property.num_isos; ++i) {
                    es.density_property.values[i] = es.density_property.values[i + 1];
                    es.density_property.colors[i] = es.density_property.colors[i + 1];
                }
                es.density_property.num_isos -= 1;
                update_rep = true;
            }
            prop::row("");
            ImGui::BeginDisabled(es.density_property.num_isos >= cap);
            if (ImGui::Button(ICON_FA_PLUS, ImVec2(btn_w, 0))) {
                const int n = es.density_property.num_isos;
                es.density_property.values[n] = es.density_property.values[n - 1] + 0.05;
                es.density_property.colors[n] = electronic_structure_density_property_default_color(es, n);
                es.density_property.num_isos += 1;
                update_rep = true;
            }
            ImGui::EndDisabled();
            ImGui::SetItemTooltip("Add an isosurface");
            prop::row("optical density", "Optical density shared by all the isosurfaces: how much they absorb of what they enclose");
            ImGui::SliderScalar("##optical_density", ImGuiDataType_Double, &es.iso_optical_density, &min_tau, &max_tau, "%.4f", ImGuiSliderFlags_Logarithmic);
        } else {
            // The band a field is evaluated in follows the isovalue
            prop::row("isovalue", electronic_structure_iso_value_tooltip(es));
            update_rep |= ImGui::SliderScalar("##isovalue", ImGuiDataType_Double, &es.iso_value, &iso_min, &iso_max, "%.4g", ImGuiSliderFlags_Logarithmic) && es.coloring == SurfaceColoring::Field;
            prop::row("optical density", "Optical density of the isosurfaces: how much they absorb of what they enclose");
            ImGui::SliderScalar("##optical_density", ImGuiDataType_Double, &es.iso_optical_density, &min_tau, &max_tau, "%.4f", ImGuiSliderFlags_Logarithmic);
        }
        prop::row("resolution", mesh_mode ? "The size of the cells the surfaces are extracted at" : "The resolution of the volume the surfaces are drawn from");
        if (ImGui::Combo("##resolution", (int*)&es.resolution, volume_resolution_str, IM_ARRAYSIZE(volume_resolution_str))) {
            update_rep = true;
        }
        prop::section_end();
    }

    if (prop::section("Color")) {
        const bool density_property = electronic_structure_is_density_property(es);
        update_rep |= draw_surface_coloring(es, sys.data.sys, !density_property, mesh_mode);
        if (es.coloring == SurfaceColoring::AtomColors) {
            update_rep |= draw_representation_coloring(sys, rep);
        }

        // The colours of the surfaces, or the tint of what colours them; a property's are per isovalue (Surface)
        if (!density_property) {
            const bool tint = es.coloring != SurfaceColoring::Uniform;
            const bool attachment = es.source == ElectronicStructureSource::TransitionDensity && es.transition_density_component == ElectronicStructureTransitionDensityComponent::Attachment;
            const bool detachment = es.source == ElectronicStructureSource::TransitionDensity && es.transition_density_component == ElectronicStructureTransitionDensityComponent::Detachment;
            const char* tint_tooltip = "Multiplies the colours of the surface, and its alpha is the surface's opacity";
            if (electronic_structure_is_signed(es)) {
                prop::row(tint ? (const char*)u8"surface tint +" : "positive", tint ? tint_tooltip : nullptr);
                prop::color("##pos", tint ? es.tint_psi_pos.elem : es.col_psi_pos.elem, true);
                prop::row(tint ? (const char*)u8"surface tint −" : "negative", tint ? tint_tooltip : nullptr);
                prop::color("##neg", tint ? es.tint_psi_neg.elem : es.col_psi_neg.elem, true);
            } else {
                vec4_t& col = attachment ? (tint ? es.tint_att : es.col_att) :
                              detachment ? (tint ? es.tint_det : es.col_det) :
                                           (tint ? es.tint_den : es.col_den);
                prop::row(tint ? "surface tint" : "color", tint ? tint_tooltip : nullptr);
                prop::color("##color", col.elem, true);
            }
        }
        prop::section_end();
    }

    return update_rep;
}

// ---------------------------------------------------------------------------------------------------------
// The Scene window (it was the Representations window): the systems loaded, and how each is shown
// ---------------------------------------------------------------------------------------------------------

#define REPRESENTATION_NODE_DND "VIAMD_REPRESENTATION_NODE"

// A representation or a group, by what names it. Also the drag and drop payload.
struct NodeRef {
    uint32_t id = 0;        // 0 for none
    bool group = false;
    uint32_t system = SYSTEM_ID_NONE;   // The id of the system it is of: dragged, it is dropped within that one only
};

// What the window does to the tree. Applied once all of it has been drawn: while the rows are drawn the
// arrays are being walked, and nothing may add to them, remove from them or reorder them.
struct WindowActions {
    struct {
        NodeRef node;
        NodeRef target;     // id 0: the top level
        RepresentationDrop drop = RepresentationDrop_IntoLast;
    } move;
    NodeRef clone;
    NodeRef remove;
    NodeRef wrap;           // into a new group, in its place
    uint32_t ungroup = 0;
    bool add_representation = false;
    bool add_group = false;
    uint32_t add_into = 0;  // the group to add into, 0 for the top level
    bool remove_all = false;
};

// The tree as it is drawn this frame
struct WindowTree {
    TreeNodes t;
    float button_size = 0.0f;
    float button_spacing = 0.0f;
    float handle_width = 0.0f;  // The same on every row, whatever its icon, so that the names line up
};

static void apply_window_actions(SystemContext& sys, const WindowActions& act) {
    if (act.remove_all) {
        remove_all_representations(sys);
        return;
    }
    if (act.move.node.id) {
        representation_move(sys, act.move.node.group, act.move.node.id, act.move.target.group, act.move.target.id, act.move.drop);
    }
    if (act.clone.id) {
        if (act.clone.group) {
            clone_representation_group(sys, act.clone.id);
        } else if (const Representation* rep = find_representation(sys, act.clone.id)) {
            clone_representation(sys, *rep);
        }
    }
    if (act.remove.id) {
        if (act.remove.group) {
            remove_representation_group(sys, act.remove.id);
        } else {
            const int idx = find_rep_idx(sys, act.remove.id);
            if (idx >= 0) remove_representation(sys, (size_t)idx);
        }
    }
    if (act.wrap.id) {
        wrap_in_group(sys, act.wrap.group, act.wrap.id);
    }
    if (act.ungroup) {
        ungroup_representation_group(sys, act.ungroup);
    }
    if (act.add_representation) {
        create_representation(sys, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT(""), act.add_into);
    }
    if (act.add_group) {
        create_representation_group(sys, STR_LIT("group"), act.add_into);
    }
}

// The menu of a row, or of the empty space below them (node.id 0)
// Whether node can be duplicated: the representations it adds fit (REPRESENTATION_MAX)
static bool can_duplicate(const SystemContext& sys, NodeRef node) {
    return representation_can_add(sys, node.group ? representation_group_rep_count(sys, node.id) : 1);
}

static const char* limit_tip() {
    static char buf[96];
    snprintf(buf, sizeof(buf), "A system has at most %zu representations", REPRESENTATION_MAX);
    return buf;
}

static void draw_tree_menu(const SystemContext& sys, NodeRef node, WindowActions* act) {
    const bool group = node.id && node.group;
    if (node.id == 0 || group) {
        if (ImGui::MenuItem(ICON_FA_PLUS " Add representation", nullptr, false, representation_can_add(sys))) {
            act->add_representation = true;
            act->add_into = group ? node.id : 0;
        }
        if (ImGui::MenuItem(ICON_FA_FOLDER_PLUS " Add group")) {
            act->add_group = true;
            act->add_into = group ? node.id : 0;
        }
    }
    if (node.id == 0) return;

    if (group) ImGui::Separator();
    if (ImGui::MenuItem(ICON_FA_OBJECT_GROUP " Group")) {
        act->wrap = node;
    }
    ImGui::SetItemTooltip("Put it in a new group, in its place");
    if (group) {
        if (ImGui::MenuItem(ICON_FA_OBJECT_UNGROUP " Ungroup")) {
            act->ungroup = node.id;
        }
        ImGui::SetItemTooltip("Remove the group and keep what is in it, in its place");
    }
    ImGui::Separator();
    if (ImGui::MenuItem(ICON_FA_COPY " Duplicate", nullptr, false, can_duplicate(sys, node))) {
        act->clone = node;
    }
    if (!can_duplicate(sys, node)) ImGui::SetItemTooltip("%s", limit_tip());
    if (ImGui::MenuItem(ICON_FA_XMARK " Remove")) {
        act->remove = node;
    }
}

#define DRAG_HANDLE_HINT "Drag to move it: between others, or onto a group to put it in"

// The handle a row is dragged by, the first thing on it. Its icon is the row's to draw (draw_handle_icon): a
// group's says whether it is open, which is known only once its tree node is drawn.
struct DragHandle {
    ImVec2 min;
    ImVec2 size;
    bool hot = false;       // Hovered or held
    bool tooltip = false;   // Hovered and not being dragged
};

static DragHandle draw_drag_handle(const WindowTree& tree, NodeRef node, const char* name) {
    DragHandle h;
    h.size = {tree.handle_width, tree.button_size};
    ImGui::InvisibleButton("##handle", h.size);
    h.min = ImGui::GetItemRectMin();
    const bool hovered = ImGui::IsItemHovered();
    h.hot = hovered || ImGui::IsItemActive();
    if (h.hot) {
        ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeAll);
    }

    if (ImGui::BeginDragDropSource(ImGuiDragDropFlags_None)) {
        ImGui::SetDragDropPayload(REPRESENTATION_NODE_DND, &node, sizeof(node));
        if (node.group) ImGui::Text("%s %s", ICON_FA_FOLDER, name);
        else            ImGui::TextUnformatted(name);
        ImGui::EndDragDropSource();
    } else {
        h.tooltip = hovered;
    }
    return h;
}

static void draw_handle_icon(const DragHandle& h, const char* icon, ImU32 col) {
    const ImVec2 icon_size = ImGui::CalcTextSize(icon);
    ImGui::GetWindowDrawList()->AddText(ImFloor(h.min + (h.size - icon_size) * 0.5f), col, icon);
}

// Where on a row something dropped goes: above it, below it or into it, by the height the mouse is at
static RepresentationDrop drop_zone(float t, bool group, bool open) {
    if (group) {
        if (t < 0.25f) return RepresentationDrop_Before;
        if (t > 0.75f) return open ? RepresentationDrop_IntoFirst : RepresentationDrop_After;
        return RepresentationDrop_IntoLast;
    }
    return t < 0.5f ? RepresentationDrop_Before : RepresentationDrop_After;
}

// The whole row as a place to drop on. x0 is where the row's content starts, which is where a line
// placing something beside it starts too.
static void draw_drop_target(const SystemContext& sys, NodeRef node, bool open, const ImRect& row, float x0, WindowActions* act) {
    if (!ImGui::BeginDragDropTargetCustom(row, ImGui::GetID("##drop"))) return;

    const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(REPRESENTATION_NODE_DND, ImGuiDragDropFlags_AcceptBeforeDelivery | ImGuiDragDropFlags_AcceptNoDrawDefaultRect);
    if (payload && payload->DataSize == (int)sizeof(NodeRef)) {
        NodeRef src;
        MEMCPY(&src, payload->Data, sizeof(src));

        const float t = (ImGui::GetMousePos().y - row.Min.y) / MAX(row.GetHeight(), 1.0f);
        const RepresentationDrop drop = drop_zone(t, node.group, open);

        if (src.system == sys.id && representation_can_move(sys, src.group, src.id, node.group, node.id, drop)) {
            ImDrawList* dl = ImGui::GetWindowDrawList();
            const ImU32 col = ImGui::GetColorU32(ImGuiCol_DragDropTarget);
            switch (drop) {
            case RepresentationDrop_Before:
                dl->AddLine({x0, row.Min.y}, {row.Max.x, row.Min.y}, col, 2.0f);
                break;
            case RepresentationDrop_After:
                dl->AddLine({x0, row.Max.y}, {row.Max.x, row.Max.y}, col, 2.0f);
                break;
            case RepresentationDrop_IntoFirst:
                dl->AddLine({x0 + ImGui::GetStyle().IndentSpacing, row.Max.y}, {row.Max.x, row.Max.y}, col, 2.0f);
                break;
            case RepresentationDrop_IntoLast:
                dl->AddRect({x0, row.Min.y}, row.Max, col, ImGui::GetStyle().FrameRounding, 0, 2.0f);
                break;
            }
            if (payload->IsDelivery()) {
                act->move.node   = src;
                act->move.target = node;
                act->move.drop   = drop;
            }
        }
    }
    ImGui::EndDragDropTarget();
}

// What was hidden missed what happened meanwhile - the frame changing, above all - as hidden representations
// are not updated. Brought up to date, but not evaluated again for nothing: an electronic structure's volume
// is evaluated again only if the frame did change (vol_hash).
static void refresh_shown_again(SystemContext& sys, NodeRef node) {
    md_temp_scope_t temp = md_temp_begin();
    TreeNodes t = tree_nodes(sys, md_temp_allocator(temp));
    for (size_t k = 0; k < t.count; ++k) {
        const bool self = t.nodes[k]->group == node.group && t.nodes[k]->id == node.id;
        if (!self && !(node.group && representation_tree_is_within(t.nodes, t.count, k, node.id))) continue;
        if (k < t.num_reps) {
            Representation& rep = sys.representation.reps[k];
            rep.needs_update = true;
            // Not one whose auto-update is off: what it picks was left where it is on purpose
            if (rep.dynamic_evaluation && representation_is_dynamic(sys, rep)) rep.filt_is_dirty = true;
        } else {
            RepresentationGroup& g = sys.representation.groups[k - t.num_reps];
            if (g.filter.result.is_dynamic) g.filt_is_dirty = true;
        }
    }
    md_temp_end(temp);
}

// The buttons at the end of a row: show/hide, duplicate, remove. Placed from the right edge.
static void draw_row_buttons(const WindowTree& tree, bool* enabled, float right_x, NodeRef node, SystemContext& sys, WindowActions* act) {
    const ImVec2 btn_size = {tree.button_size, tree.button_size};
    const float total = (tree.button_size + tree.button_spacing) * 3;
    ImGui::SameLine(right_x - total, tree.button_spacing);
    if (ImGui::Button(*enabled ? ICON_FA_EYE : ICON_FA_EYE_SLASH, btn_size)) {
        *enabled = !*enabled;
        sys.representation.atom_visibility_mask_dirty = true;
        if (*enabled) {
            refresh_shown_again(sys, node);
        }
    }
    ImGui::SetItemTooltip(node.group ? "Show/Hide the group and everything in it" : "Show/Hide");
    ImGui::SameLine(0, tree.button_spacing);
    const bool duplicable = can_duplicate(sys, node);
    ImGui::BeginDisabled(!duplicable);
    if (ImGui::Button(ICON_FA_COPY, btn_size)) {
        act->clone = node;
    }
    ImGui::EndDisabled();
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
        if (!duplicable)     ImGui::SetTooltip("%s", limit_tip());
        else if (node.group) ImGui::SetTooltip("Duplicate the group and everything in it");
        else                 ImGui::SetTooltip("Duplicate");
    }
    ImGui::SameLine(0, tree.button_spacing);
    if (ImGui::DeleteButton(ICON_FA_XMARK, btn_size)) {
        act->remove = node;
    }
    ImGui::SetItemTooltip(node.group ? "Remove the group and everything in it" : "Remove");
}

static void draw_tree_level(SystemContext& sys, const WindowTree& tree, uint32_t parent, size_t depth, WindowActions* act);

// How the atoms of a representation are coloured, as rows of its Color section: the mapping and what it takes - a
// colour, an attribute through a colour scale, the colours of the secondary structures - and for sticks the bonds.
// Also what an electronic structure coloured by its atoms takes its colours from.
bool draw_representation_coloring(SystemContext& sys, Representation& rep) {
    bool update_rep = false;

    prop::row("coloring");
    if (ImGui::BeginCombo("##atom_coloring", color_mapping_str[(int)rep.color_mapping])) {
        for (ColorMapping mapping : color_mapping_order) {
            if (ImGui::Selectable(color_mapping_str[(int)mapping], rep.color_mapping == mapping)) {
                rep.color_mapping = mapping;
                update_rep = true;
            }
        }
        ImGui::EndCombo();
    }
    if (rep.color_mapping == ColorMapping::ResType) {
        ImGui::SetItemTooltip("The colors of the residue types, set in the Residue Types of the System window");
    } else if ((rep.color_mapping == ColorMapping::ChainId || rep.color_mapping == ColorMapping::ChainIndex) && md_system_chain_count(&sys.data.sys) == 0) {
        ImGui::SetItemTooltip("The data gives no chain ids, so the system has no chains: Part Index tells its strands apart");
    } else if (rep.color_mapping == ColorMapping::PartIndex) {
        ImGui::SetItemTooltip("One color per part: each strand of a polymer, each ligand, each water");
    }

    if (rep.color_mapping == ColorMapping::Uniform) {
        prop::row("color");
        update_rep |= prop::color("##base_color", rep.base_color.elem);
    } else if (rep.color_mapping == ColorMapping::Attribute) {
        // The list of per atom fields IS the system's attribute table under atom/.
        // Queried here rather than cached anywhere, so it cannot disagree with the data.
        md_attribute_id_t prop_ids[64];
        size_t num_props = MIN(atom_attribute_query(prop_ids, ARRAY_SIZE(prop_ids), sys.data.sys), ARRAY_SIZE(prop_ids));

        const md_attributes_t& attributes = sys.data.sys.attributes;
        const md_attribute_t* selected_prop = md_attributes_get(&attributes, rep.atom_attribute.key);

        // A key the table no longer holds - a reload which dropped that field - falls
        // back to the first available rather than leaving the representation blank.
        if (!selected_prop && num_props > 0) {
            atom_attribute_select(&rep.atom_attribute, prop_ids[0], sys.data.sys);
            selected_prop = md_attributes_get(&attributes, rep.atom_attribute.key);
            update_rep = true;
        }

        if (num_props > 0 && selected_prop) {
            prop::row("attribute");
            if (ImGui::BeginCombo("##attribute", atom_attribute_label(selected_prop).ptr)) {
                for (size_t i = 0; i < num_props; ++i) {
                    const md_attribute_t* attr = md_attributes_get(&attributes, prop_ids[i]);
                    if (!attr) continue;
                    bool selected = prop_ids[i] == rep.atom_attribute.key;
                    if (ImGui::Selectable(atom_attribute_label(attr).ptr, selected)) {
                        atom_attribute_select(&rep.atom_attribute, prop_ids[i], sys.data.sys);
                        update_rep = true;
                    }
                }
                ImGui::EndCombo();
            }

            const int num_variants = atom_attribute_variant_count(selected_prop);
            if (num_variants > 1) {
                int idx = rep.atom_attribute.variant_idx + 1;
                prop::row("index");
                if (ImGui::SliderInt("##index", &idx, 1, num_variants)) {
                    update_rep = true;
                }
                rep.atom_attribute.variant_idx = CLAMP(idx - 1, 0, num_variants - 1);
            }

            // The same scale an isosurface coloured by a field has. Its span is what
            // the colours were last made from, the atoms of this representation.
            const bool multiple = num_variants > 1;
            update_rep |= color_scale_draw_controls(&rep.atom_attribute.scale, rep.atom_attribute.span, selected_prop->unit, "shown atoms",
                                                    multiple ? "the smallest to the largest value of the atoms shown, at every index"
                                                             : "the smallest to the largest value of the atoms shown");
        } else {
            prop::hint("no per atom attributes in this dataset");
        }
    } else if (rep.color_mapping == ColorMapping::SecondaryStructure) {
        prop::row("unknown");
        update_rep |= prop::color("##ss_unknown", rep.secondary_structure.color_unknown.elem);
        prop::row("coil");
        update_rep |= prop::color("##ss_coil", rep.secondary_structure.color_coil.elem);
        prop::row("helix");
        update_rep |= prop::color("##ss_helix", rep.secondary_structure.color_helix.elem);
        prop::row("sheet");
        update_rep |= prop::color("##ss_sheet", rep.secondary_structure.color_sheet.elem);
    }

    if (rep.type == RepresentationType::Licorice || rep.type == RepresentationType::BallAndStick || rep.type == RepresentationType::Bonds) {
        // How the bonds are coloured, and what that takes
        prop::row("bonds");
        update_rep |= ImGui::Combo("##bond_color", (int*)(&rep.bond_color), bond_color_mode_str, IM_ARRAYSIZE(bond_color_mode_str));
        if (rep.bond_color == BondColorMode::SmoothAtom) {
            prop::row("sharpness", "How sharply a bond turns from the colour of one atom to the other's");
            ImGui::SliderFloat("##bond_sharpness", &rep.bond_sharpness, 0.0f, 1.0f);
        } else if (rep.bond_color == BondColorMode::Uniform) {
            prop::row("bond color");
            prop::color("##bond_base_color", rep.bond_base_color.elem);
        }
    }

    // The last step of the colours, whatever made them: towards the tint, by its amount, then their saturation
    // (tint_colors, in OKLab). The colour is faint while the amount leaves it nothing to do.
    const float gap = ImGui::GetStyle().ItemInnerSpacing.x;
    const float swatch_w = 3.0f * ImGui::GetFrameHeight();
    prop::row("tint", "How far the colours go towards the tint colour: 0 not at all, 1 its hue and chroma");
    ImGui::SetNextItemWidth(MAX(prop::control_width() - swatch_w - gap, 1.0f));
    update_rep |= ImGui::SliderFloat("##tint_scale", &rep.tint_scale, 0.0f, 1.0f, "%.2f", ImGuiSliderFlags_AlwaysClamp);
    ImGui::SameLine(0.0f, gap);
    ImGui::SetNextItemWidth(swatch_w);
    if (rep.tint_scale <= 0.0f) ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.35f);
    update_rep |= prop::color("##tint_color", rep.tint_color.elem);
    if (rep.tint_scale <= 0.0f) ImGui::PopStyleVar();
    prop::row("saturation", "How much of their chroma the colours keep: 1 all of it, 0 grey");
    update_rep |= ImGui::SliderFloat("##saturation", &rep.saturation, 0.0f, 1.0f, "%.2f", ImGuiSliderFlags_AlwaysClamp);

    return update_rep;
}

// Whether a chain of the system is broken somewhere (MD_RESIDUE_FLAG_BREAK_AFTER): residues missing in a crystal
// structure, which the representations with chain gaps draw a dashed line over
static bool system_has_chain_breaks(const md_system_t& sys) {
    for (size_t i = 0; i < sys.residue.count; ++i) {
        if (md_system_residue_flags(&sys, i) & MD_RESIDUE_FLAG_BREAK_AFTER) return true;
    }
    return false;
}

// Whether the system has elastic bonds (an elastic network, restraints), kept per topology: a system may have millions
// of bonds, and the properties of a representation are drawn every frame
static bool system_has_elastic_bonds(const md_system_t& sys) {
    static uint64_t version = 0;
    static bool     result  = false;
    if (version != sys.topology_version || version == 0) {
        version = sys.topology_version;
        result = false;
        for (size_t i = 0; i < sys.bond.count && sys.bond.flags && !result; ++i) {
            result = md_bond_kind(sys.bond.flags[i]) == MD_BOND_KIND_ELASTIC;
        }
    }
    return result;
}

// The sizes of what a representation draws, as rows of its Geometry section: scales of the defaults, so 1 is as is
bool draw_representation_geometry(SystemContext& sys, Representation& rep) {
    bool update_rep = false;
    const char* scale_fmt = (const char*)u8"%.2f×";
    switch (rep.type) {
    case RepresentationType::SpaceFill:
        prop::row("radius");
        update_rep |= ImGui::SliderFloat("##radius", &rep.scale[0], 0.1f, 4.f, scale_fmt);
        break;
    case RepresentationType::Licorice:
        prop::row("radius");
        update_rep |= ImGui::SliderFloat("##radius", &rep.scale[0], 0.1f, 4.0f, scale_fmt);
        break;
    case RepresentationType::Bonds:
        prop::row("radius", "As the sticks of a licorice of the same scale");
        update_rep |= ImGui::SliderFloat("##radius", &rep.scale[0], 0.1f, 4.0f, scale_fmt);
        prop::row("dashed", "Drawn as dashes rather than sticks, as weak bonds often are");
        update_rep |= ImGui::Checkbox("##bonds_dashed", &rep.bonds_dashed);
        break;
    case RepresentationType::BallAndStick:
        prop::row("balls");
        update_rep |= ImGui::SliderFloat("##balls", &rep.scale[0], 0.1f, 4.f, scale_fmt);
        prop::row("sticks");
        update_rep |= ImGui::SliderFloat("##sticks", &rep.scale[1], 0.1f, 4.f, scale_fmt);
        break;
    case RepresentationType::Ribbons:
        prop::row("width");
        update_rep |= ImGui::SliderFloat("##width", &rep.scale[0], 0.1f, 3.f, scale_fmt);
        prop::row("thickness");
        update_rep |= ImGui::SliderFloat("##thickness", &rep.scale[1], 0.1f, 3.f, scale_fmt);
        break;
    case RepresentationType::Cartoon:
        prop::row("coil");
        update_rep |= ImGui::SliderFloat("##coil", &rep.scale[0], 0.1f, 3.f, scale_fmt);
        prop::row("sheet");
        update_rep |= ImGui::SliderFloat("##sheet", &rep.scale[1], 0.1f, 3.f, scale_fmt);
        prop::row("helix");
        update_rep |= ImGui::SliderFloat("##helix", &rep.scale[2], 0.1f, 3.f, scale_fmt);
        if (cartoon_glyphs_valid(sys.data.cartoon_glyphs)) {
            prop::row("lipid, nucleic", "Thickness of the heads and tails of the lipids and of the backbones and bases of the nucleic acids");
            update_rep |= ImGui::SliderFloat("##lipid_nucleic", &rep.scale[3], 0.1f, 3.f, scale_fmt);
            if (sys.data.sys.nucleic_backbone.range.count > 0) {
                prop::row("nucleic backbone", "Radius of the band along the backbone of the nucleic acids: half its thickness, its width stays (Å, times the lipid, nucleic scale)");
                update_rep |= ImGui::SliderFloat("##nucleic_backbone", &rep.nucleic_band_radius, 0.1f, 1.5f, "%.2f Å", ImGuiSliderFlags_AlwaysClamp);
            }
            if (sys.data.cartoon_glyphs.num_slabs > 0) {
                prop::row("nucleic rings", "Radius of the rims of the sugar and base rings, and of the stick from each sugar to its base (Å, times the lipid, nucleic scale)");
                update_rep |= ImGui::SliderFloat("##nucleic_rings", &rep.nucleic_ring_radius, 0.1f, 1.0f, "%.2f Å", ImGuiSliderFlags_AlwaysClamp);
            }
            if (sys.data.cartoon_glyphs.num_split_heads > 0) {
                prop::row("two head spheres", "The head of a phospholipid as two spheres, split at its phosphate: the head group with the phosphate, and the backbone the tails hang off");
                update_rep |= ImGui::Checkbox("##two_head_spheres", &rep.two_head_spheres);
            }
        }
        break;
    default:
        break;
    }
    // Shown where the system has a break to draw it over
    const bool gaps = rep.type == RepresentationType::Cartoon || rep.type == RepresentationType::Ribbons ||
                      rep.type == RepresentationType::Licorice || rep.type == RepresentationType::BallAndStick;
    if (gaps && system_has_chain_breaks(sys.data.sys)) {
        prop::row("chain gaps", "A dashed line over each break of a chain, where residues are missing (or its pieces are not bonded)");
        update_rep |= ImGui::Checkbox("##chain_gaps", &rep.chain_gaps);
    }
    if (rep.type == RepresentationType::Licorice || rep.type == RepresentationType::BallAndStick) {
        prop::row("bond orders", "Double and triple bonds as two and three sticks, aromatic rings with a circle (dashed in fused aromatic systems, whose sextets are not assigned to particular rings), delocalized bonds with a dashed line");
        update_rep |= ImGui::Checkbox("##bond_orders", &rep.bond_orders);
        if (sys.data.sys.bond.noncovalent_count > 0) {
            prop::row("non-covalent bonds", "The non-covalent bonds of the system among its atoms (hydrogen bonds and interactions a file gives or a user marks), dashed and thinner than the sticks. Off to leave them to a Bonds representation");
            update_rep |= ImGui::Checkbox("##noncovalent_bonds", &rep.noncovalent_bonds);
        }
        if (system_has_elastic_bonds(sys.data.sys)) {
            prop::row("elastic bonds", "The springs of the force field which are not bonds of the molecules: the elastic network of a coarse grained protein or nucleic acid, restraints. Drawn as plain bonds when shown");
            update_rep |= ImGui::Checkbox("##elastic_bonds", &rep.elastic_bonds);
        }
    }
    return update_rep;
}

bool draw_representation_type(SystemContext& sys, Representation& rep, bool atoms_only) {
    bool update_rep = false;
    // @TODO: Only display the representations which can be used for the current dataset
    prop::row("type");
    if (!rep.type_is_valid) ImGui::PushInvalid();
    if (ImGui::BeginCombo("##type", representation_type_str[(int)rep.type])) {
        const int count = atoms_only ? (int)RepresentationType::Cartoon + 1 : (int)RepresentationType::Count;
        for (int i = 0; i < count; ++i) {
            if (i == (int)RepresentationType::ElectronicStructure) {
                // Enlisted when the system carries ANY electronic structure source - orbitals are
                // one, and a file without an SCF block can still carry density properties
                if (es_source_mask(sys.data.sys) == 0) continue;
            }
            if (ImGui::Selectable(representation_type_str[(int)i], i == (int)rep.type)) {
                rep.type = (RepresentationType)i;
                if (rep.type == RepresentationType::ElectronicStructure) {
                    // The default source is orbitals, which this system may not have
                    electronic_structure_select_available_source(&rep.electronic_structure, sys.data.sys);
                }
                update_rep = true;
            }
        }
        ImGui::EndCombo();
    }
    if (!rep.type_is_valid) ImGui::PopInvalid();
    return update_rep;
}

// The properties of a representation, below its row when it is open: what it is and what it shows on rows of their
// own, then its sections. Returns whether it needs updating.
static bool draw_representation_properties(SystemContext& sys, Representation& rep) {
    ApplicationState* state = sys.app;
    bool update_rep = false;
    prop::begin("representation");

    update_rep |= draw_representation_type(sys, rep, false);

    // Of what md_gl draws and the box, where there is a cell to repeat in. Drawn so every frame: not an update.
    const bool bonds = rep.type == RepresentationType::Bonds;
    const bool periodic = ((RepresentationType::SpaceFill <= rep.type && rep.type <= RepresentationType::Cartoon) || rep.type == RepresentationType::SimulationBox || bonds) &&
                          md_unitcell_flags(&sys.data.state.unitcell) != 0;
    const bool atoms = RepresentationType::SpaceFill <= rep.type && rep.type <= RepresentationType::Cartoon;

    if (atoms || bonds) {
        // Completed within the context it is evaluated in: what is not there would not match. Either takes atoms or
        // bonds, the one as the other: the atoms of bonds, the bonds among atoms.
        const RepresentationContext ctx = representation_context(sys, rep.node.parent);
        const script::FilterCtx fctx = representation_filter_context(sys, ctx);
        prop::row("filter");
        if (filter_input("##filter", &rep.filter, state->script.language, fctx, 0, nullptr, bonds ? "all bonds (or bonds, or atoms: the bonds among them)" : "all (or atoms, or bonds: their atoms)")) {
            rep.filt_is_dirty = true;
            update_rep = true;
        }
        if (rep.filter.result.success && ImGui::IsItemHovered() && !filter_completions_open()) {
            if (ctx.group) {
                ImGui::SetTooltip("Within the context of the group '%s': %s", ctx.group->name, ctx.group->filter.text);
            } else if (bonds) {
                ImGui::SetTooltip("%zu bonds.\nA set of bonds: bond(kind='hydrogen'), bond(origin='user'), bond(bond_delta()).\n"
                                  "A set of atoms: the bonds among them, of the system.", md_array_size(rep.bond_pairs));
            }
        }

        if (periodic) {
            prop::row("periodic images", "Repeated in the periodic images of the system, as many as it is shown in (Shown as, under the system)");
            ImGui::Checkbox("##periodic_images", &rep.periodic_images);
        }
    }

    // Shown only where there is something to follow; the choice is kept either way, so that it holds
    // again for a filter edited back into depending on the frame
    if (representation_uses_atom_colors(rep) && (representation_is_dynamic(sys, rep) || rep.color_mapping == ColorMapping::Attribute)) {
        prop::row("auto-update", "Evaluate the filter and the colours again as the frame changes");
        update_rep |= ImGui::Checkbox("##auto_update", &rep.dynamic_evaluation);
        if (!rep.dynamic_evaluation) {
            ImGui::SameLine();
            if (ImGui::Button("update")) {
                rep.filt_is_dirty = true;
                update_rep = true;
            }
        }
    }

    switch (rep.type) {
    case RepresentationType::ElectronicStructure:
        update_rep |= draw_representations_window_electronic_structure(sys, rep);
        break;
    case RepresentationType::DipoleMoment:
        if (prop::section("Dipole")) {
            // The combo lists GROUPS, which is what the file publishes. A group holding one
            // moment per excited state is then picked apart by the index widget below, rather
            // than flattening every (group, element) pair into one list nobody can scan.
            DipoleGroup groups[16];
            const size_t num_groups = MIN(dipole_groups_gather(groups, ARRAY_SIZE(groups), sys.data.sys), ARRAY_SIZE(groups));

            // Asked of the key the representation holds, not searched for in the list: a key
            // that names nothing - a workspace pointed at a file without that group - simply
            // leaves the preview empty and the index widget away.
            DipoleGroup selected = {};
            const bool has_selection = dipole_group_from_key(&selected, sys.data.sys, rep.dipole.dipole_key);

            // A group can come back shorter than it was (fewer excited states), which would
            // otherwise leave the representation pointing past its end.
            if (has_selection && rep.dipole.dipole_index >= selected.count) {
                rep.dipole.dipole_index = 0;
                update_rep = true;
            }

            char preview[64] = "";
            if (has_selection) dipole_label_pretty(preview, sizeof(preview), selected.label);

            prop::row("dipole");
            if (ImGui::BeginCombo("##dipole", preview)) {
                for (size_t i = 0; i < num_groups; ++i) {
                    char lbl[64];
                    dipole_label_pretty(lbl, sizeof(lbl), groups[i].label);
                    const bool is_selected = groups[i].key == rep.dipole.dipole_key;
                    if (ImGui::Selectable(lbl, is_selected)) {
                        if (!is_selected) {
                            rep.dipole.dipole_key   = groups[i].key;
                            rep.dipole.dipole_index = 0;   // a different group; the old element means nothing in it
                            update_rep = true;
                        }
                    }
                    if (is_selected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }

            // Only an array of moments has an index to choose. One moment is addressed by the
            // group alone, and a slider from 1 to 1 is a control that cannot do anything.
            if (has_selection && selected.count > 1) {
                int index = (int)rep.dipole.dipole_index + 1;
                prop::row("index");
                if (ImGui::SliderInt("##index", &index, 1, (int)selected.count)) {
                    rep.dipole.dipole_index = (uint32_t)CLAMP(index - 1, 0, (int)selected.count - 1);
                    update_rep = true;
                }
            }
            if (num_groups == 0) {
                prop::hint("no dipole moments in this dataset");
            }

            // Drawn from the attribute every frame, so none of these is an update of the representation
            const double dipole_min = 0.01;
            const double dipole_max = 100.0;
            prop::row("color");
            prop::color("##color", rep.dipole.color.elem, true);
            prop::row("scale");
            ImGui::SliderScalar("##scale", ImGuiDataType_Double, &rep.dipole.scale, &dipole_min, &dipole_max, "%.3g", ImGuiSliderFlags_Logarithmic);
            prop::row("radius");
            ImGui::SliderFloat("##radius", &rep.dipole.radius, 0.01f, 0.5f, "%.2f");
            prop::row("offset");
            ImGui::DragFloat3("##offset", rep.dipole.offset.elem, 0.01f, -100.0f, 100.0f, "%.2f");
            prop::section_end();
        }
        break;
    case RepresentationType::SimulationBox:
        if (prop::section("Box")) {
            // Drawn from the cell every frame, so none of these is an update of the representation
            prop::row("color");
            prop::color("##color", rep.box.color.elem);
            prop::row("thickness", "Diameter of the edges, shaded as the rest of the scene and never thinner than a pixel.\n"
                                   "0: lines of a pixel, unshaded.");
            ImGui::SliderFloat("##thickness", &rep.box.thickness, 0.0f, 2.0f, rep.box.thickness > 0.0f ? "%.2f Å" : "line", ImGuiSliderFlags_AlwaysClamp);
            prop::row("axes", "The edges along a, b and c from the origin of the cell in red, green and blue");
            ImGui::Checkbox("##axes", &rep.box.axes);
            if (periodic) {
                prop::row("periodic images", "Repeated in the periodic images of the system, as many as it is shown in (Shown as, under the system)");
                ImGui::Checkbox("##periodic_images", &rep.periodic_images);
            }
            if (!rep.type_is_valid) {
                prop::hint("the system has no unit cell");
            }
            prop::section_end();
        }
        break;
    default:
        if (prop::section("Color")) {
            update_rep |= draw_representation_coloring(sys, rep);
            prop::section_end();
        }
        if (prop::section("Geometry", false)) {
            update_rep |= draw_representation_geometry(sys, rep);
            prop::section_end();
        }
        break;
    }

    prop::end();
    return update_rep;
}

// Where a row is: both kinds lay out from the left edge (the handle, the arrow, the name) and from the
// right (the buttons)
struct RowLayout {
    ImVec2 min;         // screen position the row starts at
    float right_x;      // window local x of the right edge, for SameLine
    float right_screen; // screen x of the right edge
};

static RowLayout row_begin() {
    RowLayout r;
    r.min = ImGui::GetCursorScreenPos();
    const float avail = ImGui::GetContentRegionAvail().x;
    r.right_x = ImGui::GetCursorPosX() + avail;
    r.right_screen = r.min.x + avail;
    return r;
}

static void draw_representation_row(SystemContext& sys, const WindowTree& tree, Representation& rep, WindowActions* act) {
    const NodeRef node = {rep.node.id, false, sys.id};
    ImGui::PushID("rep");
    ImGui::PushID((int)rep.node.id);

    const RowLayout row = row_begin();
    const float pad = 3.0f;

    const DragHandle handle = draw_drag_handle(tree, node, rep.name);
    draw_handle_icon(handle, ICON_FA_GRIP_VERTICAL, ImGui::GetColorU32(handle.hot ? ImGuiCol_Text : ImGuiCol_TextDisabled));
    if (handle.tooltip) ImGui::SetTooltip(DRAG_HANDLE_HINT);
    ImGui::SameLine(0, 0);

    if (rep.open_properties) {
        ImGui::SetNextItemOpen(true);
        rep.open_properties = false;
    }
    ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0, pad));
    const bool open = ImGui::TreeNodeEx("##node", ImGuiTreeNodeFlags_FramePadding);
    ImGui::PopStyleVar();

    // An invalid filter, or one in a context that is not valid, is said on the row too: the properties
    // are often closed
    ImGui::SameLine();
    const float buttons_width = (tree.button_size + tree.button_spacing) * 3;
    ImGui::SetNextItemWidth(MAX(row.right_screen - ImGui::GetCursorScreenPos().x - buttons_width - tree.button_spacing, 1.0f));
    const bool flag_invalid = representation_uses_atom_colors(rep) && !rep.filter.result.success && !rep.filt_is_dirty && rep.filter.text[0] != '\0';
    if (flag_invalid) ImGui::PushInvalid();
    ImGui::InputText("##name", rep.name, sizeof(rep.name));
    if (flag_invalid) {
        ImGui::PopInvalid();
        ImGui::SetItemTooltip("%s", rep.filter.result.error);
    }

    draw_row_buttons(tree, &rep.enabled, row.right_x, node, sys, act);

    const ImRect row_rect = {row.min, {row.right_screen, ImGui::GetItemRectMax().y}};
    draw_drop_target(sys, node, open, row_rect, row.min.x, act);

    if (ImGui::IsMouseReleased(ImGuiMouseButton_Right) && ImGui::IsWindowHovered() && row_rect.Contains(ImGui::GetMousePos())) {
        ImGui::OpenPopup("##menu");
    }
    if (ImGui::BeginPopup("##menu")) {
        draw_tree_menu(sys, node, act);
        ImGui::EndPopup();
    }

    if (open) {
        if (draw_representation_properties(sys, rep)) {
            flag_representation_as_dirty(&rep);
        }
        ImGui::Spacing();
        ImGui::TreePop();
    }

    ImGui::PopID();
    ImGui::PopID();
}

// The filter of a group, the first line within it. A grid of its own, not the one of the properties of the
// representations: that column follows the widest label drawn, so opening a representation within the group
// would move the group's context along with it. Its own, per group, keeps it where it is whatever is open.
// Completed within the context the group is in, as it is evaluated.
static void draw_group_filter(SystemContext& sys, RepresentationGroup& group) {
    ApplicationState* state = sys.app;
    char grid[32];
    snprintf(grid, sizeof(grid), "group_context_%u", group.node.id);
    prop::begin(grid);
    prop::row("context");

    const RepresentationContext outer = representation_context(sys, group.node.parent);
    const bool edited = filter_input("##context", &group.filter, state->script.language, representation_filter_context(sys, outer), 0, nullptr, "none: what is within sees the whole system");

    // An error is filter_input's to say
    const bool has_filter = group_has_filter(group);
    if (ImGui::IsItemHovered() && !filter_completions_open() && (!has_filter || group.filter.result.success)) {
        char within[96] = "";
        if (outer.group) snprintf(within, sizeof(within), ", itself within the group '%s'", outer.group->name);
        if (has_filter && !group.filt_is_dirty) {
            ImGui::SetTooltip("The filters within this group are evaluated within the %llu atoms this picks%s:\n"
                              "residue(1) is its first residue, all is all of it - as the left of 'in' in a script.",
                              (unsigned long long)md_bitset_count(group.atom_mask), within);
        } else if (!has_filter) {
            ImGui::SetTooltip("A filter here is the context of the filters within the group: they are evaluated within the atoms it\n"
                              "picks, as the left of 'in' in a script. Without one, the group only groups%s.", within);
        }
    }
    prop::end();

    if (edited) {
        flag_within_as_dirty(sys, true, group.node.id);
    }
}

static void draw_group_row(SystemContext& sys, const WindowTree& tree, RepresentationGroup& group, size_t depth, WindowActions* act) {
    const NodeRef node = {group.node.id, true, sys.id};
    ImGui::PushID("group");
    ImGui::PushID((int)group.node.id);

    const RowLayout row = row_begin();
    const float pad = 3.0f;

    // The handle is the folder, so that the name starts where the names of what is in it do, one indent in
    const DragHandle handle = draw_drag_handle(tree, node, group.name);
    ImGui::SameLine(0, 0);

    // A group is open when it is made, which is when it is to be filled
    ImGui::SetNextItemOpen(true, ImGuiCond_Once);
    ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0, pad));
    const bool open = ImGui::TreeNodeEx("##node", ImGuiTreeNodeFlags_FramePadding);
    ImGui::PopStyleVar();

    // The folder says that it is a group. In the accent colour it has a context of its own, in red one that does
    // not evaluate; the context itself is in the tooltip, for when the group is closed.
    const bool has_filter = group_has_filter(group);
    const bool invalid = has_filter && !group.filter.result.success && !group.filt_is_dirty;
    const ImU32 icon_col = invalid     ? ImGui::GetColorU32(ImVec4(0.9f, 0.3f, 0.3f, 1.0f))
                         : has_filter  ? ImGui::GetColorU32(ImGuiCol_CheckMark)
                         : ImGui::GetColorU32(handle.hot ? ImGuiCol_Text : ImGuiCol_TextDisabled);
    draw_handle_icon(handle, open ? ICON_FA_FOLDER_OPEN : ICON_FA_FOLDER, icon_col);
    if (handle.tooltip) {
        if (invalid)         ImGui::SetTooltip("Context: %s\n%s\n\n" DRAG_HANDLE_HINT, group.filter.text, group.filter.result.error);
        else if (has_filter) ImGui::SetTooltip("Context: %s\n\n" DRAG_HANDLE_HINT, group.filter.text);
        else                 ImGui::SetTooltip("A group without a context of its own\n\n" DRAG_HANDLE_HINT);
    }

    ImGui::SameLine();
    const float buttons_width = (tree.button_size + tree.button_spacing) * 3;
    ImGui::SetNextItemWidth(MAX(row.right_screen - ImGui::GetCursorScreenPos().x - buttons_width - tree.button_spacing, 1.0f));
    ImGui::InputText("##name", group.name, sizeof(group.name));

    draw_row_buttons(tree, &group.enabled, row.right_x, node, sys, act);

    const ImRect row_rect = {row.min, {row.right_screen, ImGui::GetItemRectMax().y}};
    draw_drop_target(sys, node, open, row_rect, row.min.x, act);

    if (ImGui::IsMouseReleased(ImGuiMouseButton_Right) && ImGui::IsWindowHovered() && row_rect.Contains(ImGui::GetMousePos())) {
        ImGui::OpenPopup("##menu");
    }
    if (ImGui::BeginPopup("##menu")) {
        draw_tree_menu(sys, node, act);
        ImGui::EndPopup();
    }

    if (open) {
        draw_group_filter(sys, group);

        // What is within a hidden group is drawn dimmed: it is hidden too, whatever its own eye says
        const bool dim = !group.enabled;
        if (dim) ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
        draw_tree_level(sys, tree, group.node.id, depth + 1, act);
        if (dim) ImGui::PopStyleVar();

        ImGui::TreePop();
    }

    ImGui::PopID();
    ImGui::PopID();
}

static void draw_tree_level(SystemContext& sys, const WindowTree& tree, uint32_t parent, size_t depth, WindowActions* act) {
    ApplicationState* state = sys.app;
    const TreeNodes& t = tree.t;
    if (depth > t.count) return;    // A cycle; the tree breaks those when it is normalized

    int* kids = (int*)md_alloc(state->allocator.frame, sizeof(int) * MAX(t.count, (size_t)1));
    const size_t n = representation_tree_children(kids, t.count, t.nodes, t.count, parent);
    for (size_t i = 0; i < n; ++i) {
        const size_t k = (size_t)kids[i];
        if (k < t.num_reps) {
            draw_representation_row(sys, tree, sys.representation.reps[k], act);
        } else {
            draw_group_row(sys, tree, sys.representation.groups[k - t.num_reps], depth, act);
        }
    }
}

// The representations of a system: a heading with what adds to them and removes them all, then their tree
static void draw_representations_section(SystemContext& sys, const WindowTree& tree, WindowActions* act) {
    ImGui::PushID("representations");
    const float btn = tree.button_size;
    const float spacing = tree.button_spacing;
    const float right_x = ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x;

    ImGui::AlignTextToFramePadding();
    ImGui::TextDisabled("Representations");
    // How many of as many as it can have, once that comes near
    const size_t num_reps = md_array_size(sys.representation.reps);
    if (num_reps + 4 >= REPRESENTATION_MAX) {
        ImGui::SameLine();
        if (num_reps >= REPRESENTATION_MAX) ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.2f, 1.0f), "%zu / %zu", num_reps, REPRESENTATION_MAX);
        else                                ImGui::TextDisabled("%zu / %zu", num_reps, REPRESENTATION_MAX);
        ImGui::SetItemTooltip("%s", limit_tip());
    }

    ImGui::SameLine(right_x - (btn + spacing) * 3, spacing);
    if (ImGui::Button(ICON_FA_FOLDER_PLUS, {btn, btn})) {
        act->add_group = true;
    }
    ImGui::SetItemTooltip("Add a group. What is in it is shown and hidden with it, and its filter is the\n"
                          "context the filters within it are evaluated in.");
    ImGui::SameLine(0, spacing);
    ImGui::BeginDisabled(!representation_can_add(sys));
    if (ImGui::Button(ICON_FA_PLUS, {btn, btn})) {
        act->add_representation = true;
    }
    ImGui::EndDisabled();
    ImGui::SetItemTooltip("%s", representation_can_add(sys) ? "Add a representation" : limit_tip());

    // Removing them all, behind a question: it cannot be taken back
    ImGui::SameLine(0, spacing);
    ImGui::BeginDisabled(md_array_size(sys.representation.reps) == 0 && md_array_size(sys.representation.groups) == 0);
    if (ImGui::Button(ICON_FA_TRASH_CAN, {btn, btn})) {
        ImGui::OpenPopup("##remove_all");
    }
    ImGui::EndDisabled();
    ImGui::SetItemTooltip("Remove every representation and group");
    if (ImGui::BeginPopup("##remove_all")) {
        ImGui::TextUnformatted("Remove every representation and group?");
        if (ImGui::DeleteButton("remove all")) {
            act->remove_all = true;
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("cancel")) {
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
    ImGui::PopID();

    if (tree.t.count == 0) {
        ImGui::TextDisabled("None yet. Add one with " ICON_FA_PLUS ".");
    }
    // What is within a hidden system is drawn dimmed: it is hidden too, whatever its own eye says
    const bool dim = !sys.scene.visible;
    if (dim) ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
    draw_tree_level(sys, tree, 0, 0, act);
    if (dim) ImGui::PopStyleVar();
}

// One system in the window: its header and, under it, one level in: its files, what is done to its coordinates and
// what it is shown as, then its representations. Closed, it lists none of them. The changes to it, in act, are applied
// once the window is drawn.
static void draw_system_in_window(SystemContext& sys, const WindowTree& tree, WindowActions* act, bool* remove) {
    if (draw_system_header(sys, remove)) {
        ImGui::Indent();
        draw_system_files(sys);
        draw_system_operations(sys);
        draw_representations_section(sys, tree, act);
        ImGui::Unindent();
    }
}

// The space below the rows of the last system listed: dropping there puts it last at its top level, and a right click
// adds to it
static void draw_below_rows(SystemContext& sys, WindowActions* act) {
    const ImVec2 rest = {MAX(ImGui::GetContentRegionAvail().x, 1.0f), MAX(ImGui::GetContentRegionAvail().y, ImGui::GetFrameHeight())};
    ImGui::InvisibleButton("##below", rest);
    if (ImGui::BeginDragDropTarget()) {
        const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(REPRESENTATION_NODE_DND, ImGuiDragDropFlags_AcceptBeforeDelivery | ImGuiDragDropFlags_AcceptNoDrawDefaultRect);
        if (payload && payload->DataSize == (int)sizeof(NodeRef)) {
            NodeRef src;
            MEMCPY(&src, payload->Data, sizeof(src));
            if (src.system == sys.id && representation_can_move(sys, src.group, src.id, false, 0, RepresentationDrop_IntoLast)) {
                const ImVec2 min = ImGui::GetItemRectMin();
                ImGui::GetWindowDrawList()->AddLine(min, {ImGui::GetItemRectMax().x, min.y}, ImGui::GetColorU32(ImGuiCol_DragDropTarget), 2.0f);
                if (payload->IsDelivery()) {
                    act->move.node   = src;
                    act->move.target = NodeRef{};
                    act->move.drop   = RepresentationDrop_IntoLast;
                }
            }
        }
        ImGui::EndDragDropTarget();
    }
    if (ImGui::BeginPopupContextItem("##below_menu")) {
        draw_tree_menu(sys, NodeRef{}, act);
        ImGui::EndPopup();
    }
}

void draw_representations_window(ApplicationState* state) {
    if (!state->show_representation_window) return;

    ImGui::SetNextWindowSize({300,200}, ImGuiCond_FirstUseEver);
    // The systems loaded and how each is shown. Its open state is kept in a workspace under its former name
    // (Representations), so that one written before it was renamed opens it as it was.
    if (!ImGui::Begin("Scene", &state->show_representation_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
        ImGui::End();
        return;
    }

    // Every system, in its order, each with what is changed of it. The arrays are only read while the rows are drawn
    // - every change to them waits in its actions - so the nodes gathered here hold until the end.
    struct Listed {
        SystemContext* sys;
        WindowTree tree;
        WindowActions act;
        bool remove;
    };
    Listed listed[SYSTEM_MAX];
    size_t num_listed = 0;
    const float button_size = ImGui::GetFontSize() + 3.0f * 2;
    const float handle_width = MAX(ImGui::CalcTextSize(ICON_FA_GRIP_VERTICAL).x,
                               MAX(ImGui::CalcTextSize(ICON_FA_FOLDER).x, ImGui::CalcTextSize(ICON_FA_FOLDER_OPEN).x)) + 6.0f;
    for (size_t i = 0; i < state->num_systems; ++i) {
        SystemContext& sys = *state->systems[i];
        // An empty system is not listed, but for representations which are there (a workspace whose files did not
        // load): still listed, to be removed or kept for the next
        const bool loaded = system_loaded(sys);
        if (!loaded && md_array_size(sys.representation.reps) == 0 && md_array_size(sys.representation.groups) == 0) continue;
        Listed& l = listed[num_listed++];
        l.sys = &sys;
        l.tree.t = tree_nodes(sys, state->allocator.frame);
        l.tree.button_size = button_size;
        l.tree.button_spacing = 2.0f;
        l.tree.handle_width = handle_width;
        l.act = WindowActions{};
        l.remove = false;
    }

    if (system_count(*state) == 0) {
        ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
        ImGui::TextWrapped("Nothing loaded. Open files with File > Open File..., or drop them on the window.");
        ImGui::PopStyleColor();
    }
    for (size_t i = 0; i < num_listed; ++i) {
        Listed& l = listed[i];
        ImGui::PushID((int)l.sys->id);
        if (system_loaded(*l.sys)) {
            draw_system_in_window(*l.sys, l.tree, &l.act, &l.remove);
        } else {
            ImGui::Spacing();
            draw_representations_section(*l.sys, l.tree, &l.act);
        }
        if (i + 1 == num_listed) draw_below_rows(*l.sys, &l.act);
        ImGui::PopID();
    }

    ImGui::End();

    for (size_t i = 0; i < num_listed; ++i) {
        Listed& l = listed[i];
        if (l.remove) {
            // The rows were drawn from what it held; nothing else is done to them. Taken out of the list, but for the
            // first, which stays empty for the next to be loaded into.
            system_remove(*l.sys);
            system_release(state, l.sys);
            continue;
        }
        apply_window_actions(*l.sys, l.act);
    }
}

// The legends of the representations that colour by a value: the field of an isosurface, or the
// attribute of its atoms. One per representation at most, since one representation colours by one
// thing, and only while what it shows exists - a legend of nothing would mislead.
// A legend is named by its representation's id, so it stays with the representation when they are reordered
void draw_representation_legends(ApplicationState& state) {
    const size_t num_reps = md_array_size(state.representation.reps);
    int slot = 0;
    for (size_t i = 0; i < num_reps; ++i) {
        Representation& rep = state.representation.reps[i];
        // The value under the cursor of a legend is of this frame only: none unless its legend is drawn and hovered
        rep.electronic_structure.field_isolines.hover = false;
        if (!representation_shown(state.system_context, rep)) continue;

        if (rep.type == RepresentationType::ElectronicStructure) {
            ElectronicStructureRepresentation& es = rep.electronic_structure;
            const bool has_field = state.settings.isosurface_rendering == IsosurfaceRendering_Mesh ? electronic_structure_mesh_field_span(es).valid : es.field_vol.tex_id != 0;
            if (es.coloring == SurfaceColoring::Field && es.field_map.show_legend && has_field) {
                color_scale_draw_legend(es.field_map, surface_field_kind_info(es.field_kind).unit, surface_field_kind_info(es.field_kind).label, rep.name, (int)rep.node.id, slot++, &es.field_isolines);
                continue;
            }
        }

        if (representation_uses_atom_colors(rep) && rep.color_mapping == ColorMapping::Attribute && rep.atom_attribute.scale.show_legend) {
            const md_attribute_t* attr = md_attributes_get(&state.mold.sys.attributes, rep.atom_attribute.key);
            if (!attr) continue;
            char label[128];
            atom_attribute_legend_label(label, sizeof(label), attr, rep.atom_attribute.variant_idx);
            color_scale_draw_legend(rep.atom_attribute.scale, attr->unit, label, rep.name, (int)rep.node.id, slot++);
        }
    }
}

// ## Of every system

void update_all_representations(ApplicationState* app) {
    for (size_t i = 0; i < app->num_systems; ++i) update_all_representations(*app->systems[i]);
}

void flag_all_representations_as_dirty(ApplicationState* app) {
    for (size_t i = 0; i < app->num_systems; ++i) flag_all_representations_as_dirty(*app->systems[i]);
}

void flag_dynamic_representations_as_dirty(ApplicationState* app) {
    for (size_t i = 0; i < app->num_systems; ++i) flag_dynamic_representations_as_dirty(*app->systems[i]);
}
