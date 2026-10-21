#define IMGUI_DEFINE_MATH_OPERATORS

#include <representation.h>

#include <md_filter.h>
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
static TreeNodes tree_nodes(const ApplicationState& app, md_allocator_i* alloc) {
    const size_t num_reps   = md_array_size(app.representation.reps);
    const size_t num_groups = md_array_size(app.representation.groups);
    TreeNodes t;
    t.count = num_reps + num_groups;
    t.num_reps = num_reps;
    t.nodes = (RepresentationNode**)md_alloc(alloc, sizeof(RepresentationNode*) * MAX(t.count, (size_t)1));
    for (size_t i = 0; i < num_reps; ++i) {
        t.nodes[i] = const_cast<RepresentationNode*>(&app.representation.reps[i].node);
    }
    for (size_t i = 0; i < num_groups; ++i) {
        t.nodes[num_reps + i] = const_cast<RepresentationNode*>(&app.representation.groups[i].node);
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
static void tree_normalize(ApplicationState* app) {
    md_temp_scope_t temp = md_temp_begin();
    TreeNodes t = tree_nodes(*app, md_temp_allocator(temp));
    representation_tree_normalize(t.nodes, t.count);
    md_temp_end(temp);

    sort_by_order(app->representation.reps,   md_array_size(app->representation.reps));
    sort_by_order(app->representation.groups, md_array_size(app->representation.groups));
}

static int find_rep_idx(const ApplicationState& app, uint32_t id) {
    for (size_t i = 0; i < md_array_size(app.representation.reps); ++i) {
        if (app.representation.reps[i].node.id == id) return (int)i;
    }
    return -1;
}

static int find_group_idx(const ApplicationState& app, uint32_t id) {
    if (id == 0) return -1;
    for (size_t i = 0; i < md_array_size(app.representation.groups); ++i) {
        if (app.representation.groups[i].node.id == id) return (int)i;
    }
    return -1;
}

Representation* find_representation(ApplicationState* app, uint32_t id) {
    const int i = find_rep_idx(*app, id);
    return i >= 0 ? &app->representation.reps[i] : nullptr;
}

RepresentationGroup* find_representation_group(ApplicationState* app, uint32_t id) {
    const int i = find_group_idx(*app, id);
    return i >= 0 ? &app->representation.groups[i] : nullptr;
}

const RepresentationGroup* find_representation_group(const ApplicationState& app, uint32_t id) {
    const int i = find_group_idx(app, id);
    return i >= 0 ? &app.representation.groups[i] : nullptr;
}

static bool group_has_filter(const RepresentationGroup& g) {
    return !str_empty(str_trim(str_from_cstr(g.filter.text)));
}

// Whether every group from 'parent' outwards is enabled. A parent that names no group is the top level,
// as it is to the tree.
static bool groups_enabled(const ApplicationState& app, uint32_t parent) {
    const size_t num_groups = md_array_size(app.representation.groups);
    for (size_t steps = 0; parent != 0 && steps <= num_groups; ++steps) {
        const RepresentationGroup* g = find_representation_group(app, parent);
        if (!g) break;
        if (!g->enabled) return false;
        parent = g->node.parent;
    }
    return true;
}

bool representation_shown(const ApplicationState& app, const Representation& rep) {
    return rep.enabled && groups_enabled(app, rep.node.parent);
}

bool representation_group_shown(const ApplicationState& app, const RepresentationGroup& group) {
    return group.enabled && groups_enabled(app, group.node.parent);
}

RepresentationContext representation_context(const ApplicationState& app, uint32_t parent) {
    RepresentationContext ctx;
    const size_t num_groups = md_array_size(app.representation.groups);
    for (size_t steps = 0; parent != 0 && steps <= num_groups; ++steps) {
        const RepresentationGroup* g = find_representation_group(app, parent);
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

bool representation_is_dynamic(const ApplicationState& app, const Representation& rep) {
    // The result has the context's in it as of the last evaluation; the context is asked as well for one that has
    // come to depend on the frame since, which flags what is within it for evaluation anyway
    return rep.filter.result.is_dynamic || representation_context(app, rep.node.parent).dynamic;
}

md_filter_context_t representation_filter_context(const ApplicationState& app, const RepresentationContext& ctx) {
    md_filter_context_t fctx = filter_context(&app);
    fctx.atoms = ctx.atoms;
    return fctx;
}

// Everything within a group is evaluated in a context that just changed, or may have: the groups' contexts
// and the representations' atoms. For a representation it is just itself.
static void flag_within_as_dirty(ApplicationState* app, bool group, uint32_t id) {
    if (!group) {
        if (Representation* rep = find_representation(app, id)) flag_representation_as_dirty(rep);
        app->representation.atom_visibility_mask_dirty = true;
        return;
    }

    md_temp_scope_t temp = md_temp_begin();
    TreeNodes t = tree_nodes(*app, md_temp_allocator(temp));
    for (size_t k = 0; k < t.count; ++k) {
        const RepresentationNode* node = t.nodes[k];
        const bool self = node->group && node->id == id;
        if (!self && !representation_tree_is_within(t.nodes, t.count, k, id)) continue;
        if (k < t.num_reps) {
            flag_representation_as_dirty(&app->representation.reps[k]);
        } else {
            app->representation.groups[k - t.num_reps].filt_is_dirty = true;
        }
    }
    md_temp_end(temp);
    app->representation.atom_visibility_mask_dirty = true;
}

// The ids of what is within the group, representations and groups apart, in window order
static void collect_within(md_array(uint32_t)* rep_ids, md_array(uint32_t)* group_ids, const ApplicationState& app, uint32_t id, md_allocator_i* alloc) {
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    TreeNodes t = tree_nodes(app, md_temp_allocator(temp));
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

bool representation_can_move(const ApplicationState& app, bool group, uint32_t id, bool target_group, uint32_t target_id, RepresentationDrop drop) {
    md_temp_scope_t temp = md_temp_begin();
    TreeNodes t = tree_nodes(app, md_temp_allocator(temp));
    const int idx    = representation_tree_find(t.nodes, t.count, group, id);
    const int target = target_id ? representation_tree_find(t.nodes, t.count, target_group, target_id) : -1;
    const bool ok = idx >= 0 && (target_id == 0 || target >= 0) && representation_tree_can_move(t.nodes, t.count, (size_t)idx, target, drop);
    md_temp_end(temp);
    return ok;
}

bool representation_move(ApplicationState* app, bool group, uint32_t id, bool target_group, uint32_t target_id, RepresentationDrop drop) {
    md_temp_scope_t temp = md_temp_begin();
    TreeNodes t = tree_nodes(*app, md_temp_allocator(temp));
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

    sort_by_order(app->representation.reps,   md_array_size(app->representation.reps));
    sort_by_order(app->representation.groups, md_array_size(app->representation.groups));

    // Into another group is into another context, and maybe out of sight. Among its siblings it is neither.
    if (old_parent != new_parent) {
        flag_within_as_dirty(app, group, id);
    }
    return true;
}

// ---------------------------------------------------------------------------------------------------------
// Representations
// ---------------------------------------------------------------------------------------------------------

static void init_representation(ApplicationState* state, Representation* rep) {
#if EXPERIMENTAL_GFX_API
    rep->gfx_rep = md_gfx_rep_create(state->mold.sys.atom.count);
#endif
    rep->md_rep = md_gl_rep_create(state->mold.gl_mol);
    // Of the glyphs of the system loaded now: one of a system before it goes
    md_gl_rep_destroy(rep->glyph_rep);
    rep->glyph_rep = cartoon_glyphs_rep_create(state->mold.cartoon_glyphs);
    rep->atom_mask = md_bitset_empty_set();

    // Default to the first per atom field the system offers, if it offers any - unless the
    // representation already names one the system has. This runs again for every representation
    // when a system is loaded, which is AFTER a workspace's representations were read, and for a
    // clone: neither may lose the attribute it was coloured by, or how. There is no list to consult:
    // the attribute table is the list.
    if (!md_attributes_get(&state->mold.sys.attributes, rep->atom_attribute.key)) {
        md_attribute_id_t first_attribute = MD_ATTRIBUTE_INVALID;
        if (atom_attribute_query(&first_attribute, 1, state->mold.sys) > 0) {
            atom_attribute_select(&rep->atom_attribute, first_attribute, state->mold.sys);
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
    rep->electronic_structure.density_vol.tex_id = 0;
    rep->electronic_structure.color_vol.tex_id   = 0;
    rep->electronic_structure.field_vol          = SurfaceFieldVolume{};
    rep->electronic_structure.iso_mesh           = nullptr;
}

// Appends a copy, at the place in the tree it is given, and does not normalize: that is the caller's to do
// once it has placed all it is copying. Returns its id.
static uint32_t push_copy(ApplicationState* state, const Representation& src, uint32_t parent, int32_t order) {
    // A copy first: src may well be an element of the array the push is about to grow
    Representation copy = src;
    clear_owned(&copy);
    copy.node.id     = state->representation.next_rep_id++;
    copy.node.parent = parent;
    copy.node.order  = order;
    copy.node.group  = false;
    md_array_push(state->representation.reps, copy, state->allocator.persistent);
    init_representation(state, md_array_last(state->representation.reps));
    return copy.node.id;
}

Representation* create_representation(ApplicationState* state, RepresentationType type, ColorMapping color_mapping, str_t filter, uint32_t group) {
    ASSERT(state);
    md_array_push(state->representation.reps, Representation(), state->allocator.persistent);
    Representation* rep = md_array_last(state->representation.reps);
    rep->node.id     = state->representation.next_rep_id++;
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
    es_orbital_frontier(&frontier, state->mold.sys, es_path::alpha_occupation);
    rep->electronic_structure.orbital_idx = MAX(frontier.homo_idx, 0);
    init_representation(state, rep);

    const uint32_t id = rep->node.id;
    tree_normalize(state);
    return find_representation(state, id);
}

Representation* clone_representation(ApplicationState* state, const Representation& rep) {
    ASSERT(state);
    // Taken before anything moves: rep may be an element of the array
    const Representation src = rep;
    tree_normalize(state);
    // Just after the original, which is the original's order + 1 once normalized. One that is not in the
    // array goes last.
    const Representation* orig = find_representation(state, src.node.id);
    const uint32_t id = orig ? push_copy(state, src, orig->node.parent, orig->node.order + 1)
                             : push_copy(state, src, src.node.parent, INT32_MAX);
    tree_normalize(state);
    state->representation.atom_visibility_mask_dirty = true;
    return find_representation(state, id);
}

void remove_representation(ApplicationState* state, size_t idx) {
    ASSERT(state);
    ASSERT(idx < md_array_size(state->representation.reps));
    auto& rep = state->representation.reps[idx];
    md_bitset_free(&rep.atom_mask, state->allocator.persistent);
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
    const size_t count = md_array_size(state->representation.reps);
    for (size_t i = idx; i + 1 < count; ++i) {
        state->representation.reps[i] = state->representation.reps[i + 1];
    }
    md_array_pop(state->representation.reps);
    recompute_atom_visibility_mask(state);
}

static void free_group(ApplicationState* state, size_t idx) {
    ASSERT(idx < md_array_size(state->representation.groups));
    md_bitset_free(&state->representation.groups[idx].atom_mask, state->allocator.persistent);
    const size_t count = md_array_size(state->representation.groups);
    for (size_t i = idx; i + 1 < count; ++i) {
        state->representation.groups[i] = state->representation.groups[i + 1];
    }
    md_array_pop(state->representation.groups);
}

void remove_all_representations(ApplicationState* state) {
    while (md_array_size(state->representation.reps) > 0) {
        remove_representation(state, md_array_size(state->representation.reps) - 1);
    }
    while (md_array_size(state->representation.groups) > 0) {
        free_group(state, md_array_size(state->representation.groups) - 1);
    }
}

void recompute_atom_visibility_mask(ApplicationState* state) {
    ASSERT(state);
    auto& mask = state->representation.visibility_mask;

    md_bitmask_clear(&mask);
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        const Representation& rep = state->representation.reps[i];
        if (!representation_shown(*state, rep)) continue;
        md_bitmask_or(&mask, rep.atom_mask);
    }
    state->representation.visibility_mask_hash = md_bitset_hash64(md_bitmask_view(&mask), 0);
}

// ---------------------------------------------------------------------------------------------------------
// Groups
// ---------------------------------------------------------------------------------------------------------

RepresentationGroup* create_representation_group(ApplicationState* state, str_t name, uint32_t parent) {
    ASSERT(state);
    RepresentationGroup group = {};
    group.node.id     = state->representation.next_group_id++;
    group.node.parent = parent;
    group.node.order  = INT32_MAX;
    group.node.group  = true;
    str_copy_to_char_buf(group.name, sizeof(group.name), str_empty(name) ? STR_LIT("group") : name);
    md_array_push(state->representation.groups, group, state->allocator.persistent);
    tree_normalize(state);
    return find_representation_group(state, group.node.id);
}

void clone_representation_group(ApplicationState* state, uint32_t id) {
    ASSERT(state);
    tree_normalize(state);
    const RepresentationGroup* orig = find_representation_group(state, id);
    if (!orig) return;

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    md_array(uint32_t) rep_ids = 0;
    md_array(uint32_t) group_ids = 0;
    collect_within(&rep_ids, &group_ids, *state, id, alloc);

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
        new_ids[i] = state->representation.next_group_id++;
    }
    for (size_t i = 0; i < num_groups; ++i) {
        const RepresentationGroup* src = find_representation_group(state, old_ids[i]);
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
        md_array_push(state->representation.groups, copy, state->allocator.persistent);
    }
    for (size_t i = 0; i < md_array_size(rep_ids); ++i) {
        const Representation* src = find_representation(state, rep_ids[i]);
        if (!src) continue;
        push_copy(state, *src, map_parent(src->node.parent), src->node.order);
    }

    md_temp_end(temp);
    tree_normalize(state);
    state->representation.atom_visibility_mask_dirty = true;
}

void remove_representation_group(ApplicationState* state, uint32_t id) {
    ASSERT(state);
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    md_array(uint32_t) rep_ids = 0;
    md_array(uint32_t) group_ids = 0;
    collect_within(&rep_ids, &group_ids, *state, id, alloc);
    md_array_push(group_ids, id, alloc);

    for (size_t i = 0; i < md_array_size(rep_ids); ++i) {
        const int idx = find_rep_idx(*state, rep_ids[i]);
        if (idx >= 0) remove_representation(state, (size_t)idx);
    }
    for (size_t i = 0; i < md_array_size(group_ids); ++i) {
        const int idx = find_group_idx(*state, group_ids[i]);
        if (idx >= 0) free_group(state, (size_t)idx);
    }

    md_temp_end(temp);
    tree_normalize(state);
    recompute_atom_visibility_mask(state);
}

void ungroup_representation_group(ApplicationState* state, uint32_t id) {
    ASSERT(state);
    // Normalized, everything within the group is numbered between the group and its next sibling: given the
    // group's parent, what it holds sorts exactly into its place
    tree_normalize(state);
    const RepresentationGroup* g = find_representation_group(state, id);
    if (!g) return;
    const uint32_t parent = g->node.parent;

    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        Representation& rep = state->representation.reps[i];
        if (rep.node.parent != id) continue;
        rep.node.parent = parent;
        flag_within_as_dirty(state, false, rep.node.id);
    }
    for (size_t i = 0; i < md_array_size(state->representation.groups); ++i) {
        RepresentationGroup& child = state->representation.groups[i];
        if (child.node.parent != id) continue;
        child.node.parent = parent;
        flag_within_as_dirty(state, true, child.node.id);
    }

    const int idx = find_group_idx(*state, id);
    if (idx >= 0) free_group(state, (size_t)idx);
    tree_normalize(state);
    state->representation.atom_visibility_mask_dirty = true;
}

// A new group in the place of a representation or a group, which goes into it
static void wrap_in_group(ApplicationState* state, bool group, uint32_t id) {
    tree_normalize(state);
    const RepresentationNode* node = nullptr;
    if (group) {
        if (const RepresentationGroup* g = find_representation_group(state, id)) node = &g->node;
    } else {
        if (const Representation* r = find_representation(state, id)) node = &r->node;
    }
    if (!node) return;

    RepresentationGroup wrapper = {};
    wrapper.node.id     = state->representation.next_group_id++;
    wrapper.node.parent = node->parent;
    wrapper.node.order  = node->order - 1;
    wrapper.node.group  = true;
    snprintf(wrapper.name, sizeof(wrapper.name), "group");
    md_array_push(state->representation.groups, wrapper, state->allocator.persistent);

    // The push may have moved the groups: the node is looked up again
    if (group) {
        if (RepresentationGroup* g = find_representation_group(state, id)) g->node.parent = wrapper.node.id;
    } else {
        if (Representation* r = find_representation(state, id)) r->node.parent = wrapper.node.id;
    }
    // A group without a filter changes no context
    tree_normalize(state);
}

// ---------------------------------------------------------------------------------------------------------
// Evaluation
// ---------------------------------------------------------------------------------------------------------

// The context of a group: its own filter evaluated within the context it is in, which comes first
static void update_group(ApplicationState* state, size_t idx, size_t depth) {
    RepresentationGroup* g = &state->representation.groups[idx];
    if (!g->filt_is_dirty) return;
    if (depth > md_array_size(state->representation.groups)) return;    // A cycle; the tree breaks those

    const int parent_idx = find_group_idx(*state, g->node.parent);
    if (parent_idx >= 0) {
        update_group(state, (size_t)parent_idx, depth + 1);
    }

    if (!group_has_filter(*g)) {
        md_bitset_free(&g->atom_mask, state->allocator.persistent);
        g->filter.result = {};
        g->filter.result.success = true;
    } else {
        // Without a system there is nothing to evaluate against; it is done once there is one
        if (state->mold.sys.atom.count == 0) return;

        const RepresentationContext ctx = representation_context(*state, g->node.parent);
        if (!ctx.valid) {
            md_bitset_free(&g->atom_mask, state->allocator.persistent);
            g->filter.result = {};
            // A context that follows the frame can come back valid in another frame, and this with it
            g->filter.result.is_dynamic = ctx.dynamic;
            snprintf(g->filter.result.error, sizeof(g->filter.result.error), "The filter of the group '%s' it is in is not valid", ctx.group ? ctx.group->name : "");
        } else {
            filter_evaluate(&g->filter, &g->atom_mask, representation_filter_context(*state, ctx), state->allocator.persistent);
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

bool representation_atom_type_valid(const ApplicationState& state, RepresentationType type) {
    const md_system_t& sys = state.mold.sys;
    switch (type) {
    case RepresentationType::SpaceFill:     return sys.atom.count > 0;
    case RepresentationType::Licorice:      return sys.bond.count > 0;
    case RepresentationType::BallAndStick:  return sys.atom.count > 0;
    case RepresentationType::Ribbons:       return sys.protein_backbone.range.count > 0;
    // Proteins, and the lipids and nucleic acids of the glyphs
    case RepresentationType::Cartoon:       return sys.protein_backbone.range.count > 0 || cartoon_glyphs_valid(state.mold.cartoon_glyphs);
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
void representation_glyph_draw_ops(md_array(md_gl_draw_op_t)* ops, const ApplicationState& state, const Representation& rep, md_allocator_i* alloc) {
    if (!cartoon_glyphs_valid(state.mold.cartoon_glyphs) || !rep.glyph_rep.id || !rep.glyphs_shown || !state.mold.cartoon_glyphs.positions_valid) return;
    const float scale = rep.scale.w > 0.0f ? rep.scale.w : 1.0f;

    // The heads; the other points (control points, corners of slabs) have no sphere, and nucleic acids none at all
    if (state.mold.cartoon_glyphs.num_spheres > 0) {
        md_gl_draw_op_t spheres = {
            .type = MD_GL_REP_SPACE_FILL,
            .args = {},
            .rep = rep.glyph_rep,
            .copies = {},
        };
        spheres.args.space_fill.radius_scale = scale;
        md_array_push(*ops, spheres, alloc);
    }

    if (state.mold.cartoon_glyphs.num_tubes > 0) {
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

    if (state.mold.cartoon_glyphs.num_slabs > 0) {
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
void representation_atom_colors(uint32_t* colors, size_t num_atoms, ApplicationState* state, Representation* rep, const md_bitset_t* shown) {
    const md_system_t& sys = state->mold.sys;
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
    case ColorMapping::CompName:
        color_atoms_comp_name(colors, num_atoms, sys);
        break;
    case ColorMapping::CompType:
        color_atoms_comp_type(colors, num_atoms, sys);
        break;
    case ColorMapping::CompSeqId:
        color_atoms_comp_seq_id(colors, num_atoms, sys);
        break;
    case ColorMapping::CompIndex:
        color_atoms_comp_idx(colors, num_atoms, sys);
        break;
    case ColorMapping::InstId:
        color_atoms_inst_id(colors, num_atoms, sys);
        break;
    case ColorMapping::InstIndex:
        color_atoms_inst_idx(colors, num_atoms, sys);
        break;
    case ColorMapping::SecondaryStructure: {
        SecondaryStructurePalette palette = {
            .coil  = convert_color(rep->secondary_structure.color_coil),
            .helix = convert_color(rep->secondary_structure.color_helix),
            .sheet = convert_color(rep->secondary_structure.color_sheet),
        };
        color_atoms_secondary_structure(colors, num_atoms, sys, displayed_secondary_structure(state), palette);
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

static void update_representation(ApplicationState* state, Representation* rep) {
    ASSERT(state);
    ASSERT(rep);

    // A hidden one keeps what needs doing for when it is shown again
    if (!rep->needs_update) return;
    if (!representation_shown(*state, *rep)) return;

    const auto& sys = state->mold.sys;
    size_t num_atoms = md_system_atom_count(&sys);

    md_allocator_i* frame_alloc = state->allocator.frame;
    md_temp_scope_t temp = md_temp_begin_in(frame_alloc);
    defer { md_temp_end(temp); };

    uint32_t* colors = 0;
    if (representation_uses_atom_colors(*rep)) {
        colors = (uint32_t*)md_vm_arena_push(frame_alloc, sizeof(uint32_t) * num_atoms);

        // The atoms first: the colours are applied through them, and a colour scale's range can
        // follow the values of exactly the atoms shown
        if (rep->dynamic_evaluation && representation_is_dynamic(*state, *rep)) {
            rep->filt_is_dirty = true;
        }
        if (rep->filt_is_dirty) {
            // Within the context of the groups it is in. A context that is not valid leaves nothing to
            // evaluate within, which is what an invalid filter of its own leaves too.
            const RepresentationContext ctx = representation_context(*state, rep->node.parent);
            if (ctx.valid) {
                filter_evaluate(&rep->filter, &rep->atom_mask, representation_filter_context(*state, ctx), state->allocator.persistent);
                // What it picks follows the frame if the context it picks within does
                rep->filter.result.is_dynamic |= ctx.dynamic;
            } else {
                md_bitset_free(&rep->atom_mask, state->allocator.persistent);
                rep->filter.result = {};
                rep->filter.result.is_dynamic = ctx.dynamic;     // it may come back valid in another frame
                snprintf(rep->filter.result.error, sizeof(rep->filter.result.error), "The filter of the group '%s', the context of this one, is not valid", ctx.group ? ctx.group->name : "");
            }
            rep->filt_is_dirty = false;
            state->representation.atom_visibility_mask_dirty = true;
        }

        // Of the atoms shown, for a colour scale which spans their values
        representation_atom_colors(colors, num_atoms, state, rep, rep->filter.result.success ? &rep->atom_mask : nullptr);
    }

    switch (rep->type) {
    case RepresentationType::SpaceFill:
    case RepresentationType::Licorice:
    case RepresentationType::BallAndStick:
    case RepresentationType::Ribbons:
    case RepresentationType::Cartoon:
        rep->type_is_valid = representation_atom_type_valid(*state, rep->type);
        break;
    case RepresentationType::ElectronicStructure:
        rep->type_is_valid = electronic_structure_source_supported(es_source_mask(sys), rep->electronic_structure.source);
        if (rep->type_is_valid) {
            electronic_structure_representation_update(state, rep, colors);
        }
        break;
    case RepresentationType::DipoleMoment:
        rep->type_is_valid = md_attributes_get(&state->mold.sys.attributes, rep->dipole.dipole_key) != NULL;
        break;
    case RepresentationType::SimulationBox:
        rep->type_is_valid = md_unitcell_flags(&state->mold.state.unitcell) != 0;
        break;
    default:
        ASSERT(false);
        break;
    }

    if (colors) {
        if (rep->filter.result.success) {
            filter_colors(colors, num_atoms, rep->atom_mask);
            state->representation.atom_visibility_mask_dirty = true;
            md_gl_rep_set_atom_colors(rep->md_rep, 0, (uint32_t)num_atoms, colors, 0);
            rep->glyphs_shown = false;
            if (rep->type == RepresentationType::Cartoon) {
                rep->glyphs_shown = cartoon_glyphs_rep_set_colors(state->mold.cartoon_glyphs, rep->glyph_rep, colors, num_atoms, rep->two_head_spheres);
            }
            state->representation.gl_version += 1;

#if EXPERIMENTAL_GFX_API
            md_gfx_rep_attr_t attributes = {};
            attributes.spacefill.radius_scale = 1.0f;
            md_gfx_rep_set_type_and_attr(rep->gfx_rep, MD_GFX_REP_TYPE_SPACEFILL, &attributes);
            md_gfx_rep_set_color(rep->gfx_rep, 0, (uint32_t)mol.atom.count, (md_gfx_color_t*)colors, 0);
#endif
        }
    }

    rep->needs_update = false;
}

void update_all_representations(ApplicationState* state) {
    // The contexts first, so that every representation is evaluated within an up to date one. Hidden groups
    // wait: what is within them is hidden too, and nothing evaluated there would be seen.
    for (size_t i = 0; i < md_array_size(state->representation.groups); ++i) {
        if (representation_group_shown(*state, state->representation.groups[i])) {
            update_group(state, i, 0);
        }
    }
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        update_representation(state, &state->representation.reps[i]);
    }
}

void init_all_representations(ApplicationState* state) {
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        init_representation(state, &state->representation.reps[i]);
    }
    // The contexts were selections of the system that was there before
    for (size_t i = 0; i < md_array_size(state->representation.groups); ++i) {
        RepresentationGroup& g = state->representation.groups[i];
        md_bitset_free(&g.atom_mask, state->allocator.persistent);
        g.filt_is_dirty = true;
    }
    tree_normalize(state);
}

void flag_representation_as_dirty(Representation* rep) {
    ASSERT(rep);
    rep->filt_is_dirty = true;
    rep->needs_update  = true;
    rep->electronic_structure.vol_hash = 0;
}

void flag_all_representations_as_dirty(ApplicationState* state) {
    ASSERT(state);
    for (size_t i = 0; i < md_array_size(state->representation.groups); ++i) {
        state->representation.groups[i].filt_is_dirty = true;
    }
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        flag_representation_as_dirty(&state->representation.reps[i]);
    }
}

void flag_dynamic_representations_as_dirty(ApplicationState* state) {
    ASSERT(state);
    // A context that follows the frame is evaluated again; what is within it decides for itself whether to
    // follow, by its own auto-update
    for (size_t i = 0; i < md_array_size(state->representation.groups); ++i) {
        RepresentationGroup& g = state->representation.groups[i];
        if (g.filter.result.is_dynamic && representation_group_shown(*state, g)) {
            g.filt_is_dirty = true;
        }
    }
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        Representation& rep = state->representation.reps[i];
        if (!representation_shown(*state, rep)) continue;

        const bool atoms_follow  = rep.dynamic_evaluation && representation_is_dynamic(*state, rep);
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
void create_default_representations(ApplicationState* state) {
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
    const bool electronic_structure_present = es_source_mask(state->mold.sys) != 0;


    if (state->mold.sys.atom.count > 3'000'000) {
        VIAMD_LOG_INFO("Large system detected, creating default representation for all atoms");
        Representation* rep = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("all"));
        snprintf(rep->name, sizeof(rep->name), "default");
        goto done;
    }

    if (state->mold.sys.component.count == 0) {
        // No residues present
        Representation* rep = create_representation(state, RepresentationType::BallAndStick, ColorMapping::Type, STR_LIT("all"));
        snprintf(rep->name, sizeof(rep->name), "default");
        goto done;
    }

    // What is present, by what the components are. The selections of the representations below are by the same.
    for (size_t i = 0; i < state->mold.sys.component.count; ++i) {
        switch (md_component_kind(&state->mold.sys.component, i)) {
        case MD_COMPONENT_KIND_AMINO_ACID: amino_acid_present = true; break;
        case MD_COMPONENT_KIND_NUCLEOTIDE: nucleic_present = true; num_nucleotides += 1; break;
        case MD_COMPONENT_KIND_ION:        ion_present = true;        break;
        case MD_COMPONENT_KIND_WATER:      water_present = true;      break;
        case MD_COMPONENT_KIND_LIPID:      lipid_present = true;      break;
        default:                           ligand_present = true;     break;
        }
    }
    for (size_t i = 0; i < state->mold.sys.atom.count; ++i) {
        num_coarse_grained += md_atom_particle_kind(&state->mold.sys.atom, i) == MD_PARTICLE_BEAD;
    }

    // Lipids are shown as their cartoon where their heads and tails are known (the atoms of 'head or tail'), and the
    // rest of them (lipids by name whose anatomy is not known, as without bonds) as ligands
    if (lipid_present) {
        const md_lipid_data_t& lipids = state->mold.sys.lipid;
        size_t num_lipid_comps = 0;
        size_t num_known = 0;
        for (size_t i = 0; i < state->mold.sys.component.count; ++i) {
            num_lipid_comps += md_component_kind(&state->mold.sys.component, i) == MD_COMPONENT_KIND_LIPID;
        }
        for (size_t i = 0; i < lipids.count; ++i) {
            if (lipids.head_offset[i + 1] > lipids.head_offset[i]) num_known += md_index_range_size(&lipids.comps, i);
        }
        lipid_cartoon = num_known > 0 && cartoon_glyphs_valid(state->mold.cartoon_glyphs);
        if (!lipid_cartoon || num_known < num_lipid_comps) ligand_present = true;
    }
    if (lipid_cartoon) {
        Representation* lip = create_representation(state, RepresentationType::Cartoon, ColorMapping::CompType, STR_LIT("lipid"));
        snprintf(lip->name, sizeof(lip->name), "lipid");
    }

    // Coarse grained when most of it is. A handful of beads, or atoms a loader could not assign an
    // element to, should not take the protein, ligand and water representations from the rest.
    if (2 * num_coarse_grained > state->mold.sys.atom.count) {
        // Proteins of beads have a trace backbone and its secondary structure (CG-DSSP): a cartoon. So do nucleic acids
        // of beads where every nucleotide is on a backbone, their bases with the rings of ideal ones. The rest as
        // beads.
        const bool protein_trace = amino_acid_present && state->mold.sys.protein_backbone.segment.count > 0;
        const size_t num_on_backbone = num_nucleotides_on_backbone(state->mold.sys);
        const bool nucleic_cartoon = nucleic_present && num_on_backbone >= 10 && num_on_backbone == num_nucleotides && cartoon_glyphs_valid(state->mold.cartoon_glyphs);
        if (protein_trace) {
            Representation* prot = create_representation(state, RepresentationType::Cartoon, ColorMapping::SecondaryStructure, STR_LIT("protein"));
            snprintf(prot->name, sizeof(prot->name), "protein");
        }
        if (nucleic_cartoon) {
            Representation* nucl = create_representation(state, RepresentationType::Cartoon, ColorMapping::CompType, STR_LIT("nucleic"));
            snprintf(nucl->name, sizeof(nucl->name), "nucleic");
        }
        const char* shown[3];
        int num_shown = 0;
        if (protein_trace)   shown[num_shown++] = "protein";
        if (nucleic_cartoon) shown[num_shown++] = "nucleic";
        if (lipid_cartoon)   shown[num_shown++] = "head or tail";
        char rest[64] = "all";
        if (num_shown == 1 && !strchr(shown[0], ' ')) {
            snprintf(rest, sizeof(rest), "not %s", shown[0]);
        } else if (num_shown > 0) {
            int len = snprintf(rest, sizeof(rest), "not (");
            for (int i = 0; i < num_shown; ++i) {
                len += snprintf(rest + len, sizeof(rest) - len, "%s%s", i > 0 ? " or " : "", shown[i]);
            }
            snprintf(rest + len, sizeof(rest) - len, ")");
        }
        Representation* rep = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Type, str_from_cstr(rest));
        snprintf(rep->name, sizeof(rep->name), "default");
        goto done;
    }

    if (amino_acid_present) {
        RepresentationType type = RepresentationType::Cartoon;
        ColorMapping color = ColorMapping::SecondaryStructure;

        // Several chains are told apart by color; a short single chain (or none: free amino acids) is shown atom by atom
        size_t num_chains = 0;
        size_t chain_res_count = 0;
        for (size_t i = 0; i < state->mold.sys.instance.count; ++i) {
            if (md_system_instance_entity_kind(&state->mold.sys, i) == MD_ENTITY_KIND_PEPTIDE) {
                if (num_chains++ == 0) chain_res_count = md_instance_comp_count(&state->mold.sys.instance, i);
            }
        }
        if (num_chains > 1) {
            color = ColorMapping::InstId;
        } else if (chain_res_count < 20) {
            type = RepresentationType::BallAndStick;
            color = ColorMapping::Type;
        }

        Representation* prot = create_representation(state, type, color, STR_LIT("protein"));
        snprintf(prot->name, sizeof(prot->name), "protein");
    }
    if (nucleic_present) {
        // A cartoon where every nucleotide is on a backbone (and there are some of them), else atom by atom: the cartoon
        // shows nothing of a nucleotide off a backbone
        const size_t num_on_backbone = num_nucleotides_on_backbone(state->mold.sys);
        const bool cartoon = num_on_backbone >= 10 && num_on_backbone == num_nucleotides && cartoon_glyphs_valid(state->mold.cartoon_glyphs);
        Representation* nucl = create_representation(state, cartoon ? RepresentationType::Cartoon : RepresentationType::BallAndStick, cartoon ? ColorMapping::CompType : ColorMapping::Type, STR_LIT("nucleic"));
        snprintf(nucl->name, sizeof(nucl->name), "nucleic");
    }
    if (ion_present) {
        Representation* ion = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("ion"));
        snprintf(ion->name, sizeof(ion->name), "ion");
    }
    if (ligand_present) {
        Representation* ligand = create_representation(state, RepresentationType::BallAndStick, ColorMapping::Type, lipid_cartoon ? STR_LIT("not (protein or nucleic or water or ion or head or tail)") : STR_LIT("not (protein or nucleic or water or ion)"));
        snprintf(ligand->name, sizeof(ligand->name), "ligand");
    }
    if (water_present) {
        Representation* water = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("water"));
        water->scale.x = 0.5f;
        snprintf(water->name, sizeof(water->name), "water");
        water->enabled = false;
        if (!amino_acid_present && !nucleic_present && !ligand_present && !lipid_present) {
            water->enabled = true;
        }
    }

done:
    // The cell, where there is one, hidden: there to be shown, for the coarse grained and the large too
    if (md_unitcell_flags(&state->mold.state.unitcell) != 0) {
        Representation* box = create_representation(state, RepresentationType::SimulationBox);
        snprintf(box->name, sizeof(box->name), "box");
        box->enabled = false;
    }

    if (electronic_structure_present) {
        Representation* rep = create_representation(state, RepresentationType::ElectronicStructure);
        snprintf(rep->name, sizeof(rep->name), "electronic structure");
        rep->enabled = true;
        electronic_structure_select_available_source(&rep->electronic_structure, state->mold.sys);


		// ONE representation, for the ground state moment. A group each meant a VeloxChem file with
		// transition dipoles opened with a representation per group, and with the electric,
		// magnetic and velocity sets that is three nobody asked for on top of the one they wanted.
		// The ground state is what a dipole means to someone who has not said otherwise; the rest
		// are a representation away, and the index slider covers the excited states within each.
		DipoleGroup groups[16];
		size_t num_groups = MIN(dipole_groups_gather(groups, ARRAY_SIZE(groups), state->mold.sys), ARRAY_SIZE(groups));
        for (size_t i = 0; i < num_groups; ++i) {
            // The group name is the identity here, the same string the path spells.
            if (!str_eq(groups[i].label, STR_LIT("ground_state"))) continue;

            vec3_t vec = {0, 0, 0};
            if (!dipole_moment_read(&vec, nullptr, state->mold.sys, groups[i].key, 0)) continue;
            if (vec3_length(vec) <= 1e-3f) continue;

            Representation* dipole_rep = create_representation(state, RepresentationType::DipoleMoment);
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

    recompute_atom_visibility_mask(state);
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
    } else if (str_eq(ident, STR_LIT("Opacity"))) {
        viamd::extract_flt(rep->opacity, arg);
    } else if (str_eq(ident, STR_LIT("IOR"))) {
        viamd::extract_flt(rep->ior, arg);
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
    } else if (str_eq(ident, STR_LIT("ChainGaps"))) {
        viamd::extract_bool(rep->chain_gaps, arg);
    } else if (str_eq(ident, STR_LIT("ElasticBonds"))) {
        viamd::extract_bool(rep->elastic_bonds, arg);
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

static void deserialize_representation(ApplicationState* data, viamd::deserialization_state_t& state) {
    Representation* rep = create_representation(data);
    str_t ident, arg;
    while (viamd::next_entry(ident, arg, state)) {
        if (str_eq(ident, STR_LIT("Group"))) {
            // By the id of a [RepresentationGroup] before it. One that names nothing leaves it at the top level.
            int group = 0;
            if (viamd::extract_int(group, arg) && group > 0 && find_representation_group(data, (uint32_t)group)) {
                rep->node.parent = (uint32_t)group;
            }
        } else {
            representation_read_entry(rep, ident, arg);
        }
    }
    // Last in the group it says it is in: the sections are in the order of the window
    tree_normalize(data);
}

void representation_write_entries(viamd::serialization_state_t& state, const ApplicationState* app_state, const Representation& rep) {
    const md_attributes_t* attributes = &app_state->mold.sys.attributes;

    viamd::write_str(state,  STR_LIT("Name"), str_from_cstr(rep.name));
    viamd::write_str(state,  STR_LIT("Filter"), str_from_cstr(rep.filter.text));
    viamd::write_bool(state, STR_LIT("Enabled"), rep.enabled);
    viamd::write_int(state,  STR_LIT("Type"), (int)rep.type);
    viamd::write_int(state,  STR_LIT("ColorMapping"), (int)rep.color_mapping);
    viamd::write_vec4(state, STR_LIT("BaseColor"), rep.base_color);
    viamd::write_flt(state,  STR_LIT("Saturation"), rep.saturation);
    viamd::write_flt(state,  STR_LIT("Opacity"), rep.opacity);
    viamd::write_flt(state,  STR_LIT("IOR"), rep.ior);
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
    viamd::write_bool(state, STR_LIT("ChainGaps"), rep.chain_gaps);
    viamd::write_bool(state, STR_LIT("ElasticBonds"), rep.elastic_bonds);
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

static void serialize_representation(viamd::serialization_state_t& state, const ApplicationState* app_state, const Representation& rep) {
    viamd::write_section_header(state, STR_LIT("Representation"));
    if (find_representation_group(*app_state, rep.node.parent)) {
        viamd::write_int(state, STR_LIT("Group"), rep.node.parent);
    }
    representation_write_entries(state, app_state, rep);
}

static void deserialize_group(ApplicationState* data, viamd::deserialization_state_t& state) {
    RepresentationGroup group = {};
    int id = 0;
    int parent = 0;
    str_t ident, arg;
    while (viamd::next_entry(ident, arg, state)) {
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
    if (id > 0 && !find_representation_group(data, (uint32_t)id)) {
        group.node.id = (uint32_t)id;
        data->representation.next_group_id = MAX(data->representation.next_group_id, (uint32_t)id + 1);
    } else {
        group.node.id = data->representation.next_group_id++;
    }
    // A group comes before what is within it, so the group it is in has been read already
    if (parent > 0 && find_representation_group(data, (uint32_t)parent)) {
        group.node.parent = (uint32_t)parent;
    }
    group.node.order = INT32_MAX;
    group.node.group = true;
    group.filt_is_dirty = true;
    md_array_push(data->representation.groups, group, data->allocator.persistent);
    tree_normalize(data);
}

static void serialize_group(viamd::serialization_state_t& state, const ApplicationState* app, const RepresentationGroup& group) {
    viamd::write_section_header(state, STR_LIT("RepresentationGroup"));
    viamd::write_int(state, STR_LIT("Id"), group.node.id);
    if (find_representation_group(*app, group.node.parent)) {
        viamd::write_int(state, STR_LIT("Group"), group.node.parent);
    }
    viamd::write_str(state,  STR_LIT("Name"), str_from_cstr(group.name));
    viamd::write_str(state,  STR_LIT("Filter"), str_from_cstr(group.filter.text));
    viamd::write_bool(state, STR_LIT("Enabled"), group.enabled);
}

// Depth first, the order of the window: a group is written before what is within it, which is what reading
// it back relies on
static void serialize_level(viamd::serialization_state_t& state, const ApplicationState* app, const TreeNodes& t, uint32_t parent, size_t depth) {
    if (depth > t.count) return;
    md_temp_scope_t temp = md_temp_begin();
    int* kids = (int*)md_temp_alloc(temp, sizeof(int) * MAX(t.count, (size_t)1));
    const size_t n = representation_tree_children(kids, t.count, t.nodes, t.count, parent);
    for (size_t i = 0; i < n; ++i) {
        const size_t k = (size_t)kids[i];
        if (k < t.num_reps) {
            serialize_representation(state, app, app->representation.reps[k]);
        } else {
            const RepresentationGroup& group = app->representation.groups[k - t.num_reps];
            serialize_group(state, app, group);
            serialize_level(state, app, t, group.node.id, depth + 1);
        }
    }
    md_temp_end(temp);
}

void representations_serialize(viamd::serialization_state_t& state, const ApplicationState* app) {
    md_temp_scope_t temp = md_temp_begin();
    const TreeNodes t = tree_nodes(*app, md_temp_allocator(temp));
    serialize_level(state, app, t, 0, 0);
    md_temp_end(temp);
}

bool representations_deserialize(ApplicationState* app, str_t section, viamd::deserialization_state_t& state) {
    if (str_eq(section, STR_LIT("Representation"))) {
        deserialize_representation(app, state);
        return true;
    }
    if (str_eq(section, STR_LIT("RepresentationGroup"))) {
        deserialize_group(app, state);
        return true;
    }
    return false;
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
                    es.field_kind = surface_field_first_available(sys);
                }
                update_rep = true;
            }
            ImGui::EndDisabled();
            if (!enabled && c == SurfaceColoring::Field) {
                ImGui::SetItemTooltip("Nothing in this dataset to make a field from:\n"
                                      "the electrostatic potential needs a basis, a ground state density and the QM atoms' nuclear charges,\n"
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
    if (!mesh_mode && vol.num_surface_samples > 0) {
        snprintf(status, sizeof(status), "%zu of %zu voxels (%dx%dx%d) evaluated on the %s", vol.num_evaluated, vol.num_voxels,
                 vol.grid.dim[0], vol.grid.dim[1], vol.grid.dim[2], vol.on_gpu ? "GPU" : "CPU");
    }
    prop::row("field", status);
    if (ImGui::BeginCombo("##field", surface_field_kind_str[(int)es.field_kind])) {
        for (int k = 0; k < (int)SurfaceFieldKind::Count; ++k) {
            const SurfaceFieldKind kind = (SurfaceFieldKind)k;
            ImGui::BeginDisabled(!surface_field_available(kind, sys));
            if (ImGui::Selectable(surface_field_kind_str[k], es.field_kind == kind)) {
                es.field_kind = kind;
                update_rep = true;
            }
            ImGui::EndDisabled();
        }
        ImGui::EndCombo();
    }
    if (es.field_kind == SurfaceFieldKind::EmbeddingPotential) {
        ImGui::SetItemTooltip("The electrostatic potential of the environment's classical multipoles on the surface: the charges,\n"
                              "and the dipoles and quadrupoles where the potential has them.\n"
                              "The QM region is left out, whatever charges a topology gives its atoms.\n"
                              "For a polarizable embedding this is the potential of its PERMANENT multipoles:\n"
                              "the induced dipoles are solved for during the calculation and not stored.");
    } else if (es.field_kind == SurfaceFieldKind::ElectrostaticPotential) {
        ImGui::SetItemTooltip("The electrostatic potential of the QM region on the surface: its nuclei and its ground state\n"
                              "electron density, exact over the basis. Negative (red) where the electrons dominate - lone pairs,\n"
                              "pi systems, where an electrophile is drawn - and positive (blue) where the nuclei are less screened.\n"
                              "Mapped onto an isodensity surface (0.001-0.002 au) it is the molecular electrostatic potential.\n"
                              "An embedding's environment is not part of it.");
    }

    // The colours are applied per pixel as the surface is shaded, so a change to the scale needs
    // nothing re-evaluated: it is not an update of the representation
    // With mesh rendering the field is evaluated at the vertices of the meshes, not on a volume
    color_scale_draw_controls(&es.field_map, mesh_mode ? electronic_structure_mesh_field_span(es) : surface_field_span(vol),
                              surface_field_unit(es.field_kind), "on surface", "its 1st to 99th percentile on the surface");
    // Drawn as the surface is shaded, like the colours: nothing to re-evaluate
    color_scale_draw_isoline_controls(&es.field_isolines, es.field_map, surface_field_unit(es.field_kind));

    return update_rep;
}


// The properties of an electronic structure representation: what the volume is of (Source), the surfaces drawn from
// it (Surface) and how they are coloured (Color)
static bool draw_representations_window_electronic_structure(ApplicationState* state, Representation& rep) {
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
    const ElectronicStructureSourceFlags source_mask = es_source_mask(state->mold.sys);

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
            const int num_density_props = (int)MIN(density_properties_gather(density_props, ARRAY_SIZE(density_props), state->mold.sys), ARRAY_SIZE(density_props));
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
            const bool has_beta = es_has_distinct_beta_orbitals(state->mold.sys);
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
            const bool use_beta = (es.spin == ElectronicStructureSpin::Beta) && es_has_distinct_beta_orbitals(state->mold.sys);

            size_t num_orbitals = 0;
            OrbitalFrontier frontier = {};
            if (es_orbital_extent(state->mold.sys, &num_orbitals, nullptr) && num_orbitals > 0) {
                es_orbital_frontier(&frontier, state->mold.sys, use_beta ? es_path::beta_occupation : es_path::alpha_occupation);

                // Energies alongside the index, so the list can be read as an orbital diagram. A column
                // that does not span the orbital set is not one to index with, so it is dropped whole
                // rather than bounds checked per row - the list then falls back to bare indices.
                size_t num_energies = 0;
                const double* energy = es_orbital_energies(&num_energies, state->mold.sys, use_beta ? es_path::beta_energy : es_path::alpha_energy);
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
            const int num_excited_states = (int)es_excited_state_count(state->mold.sys);
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
                    const int num_lambdas = (int)es_nto_lambdas(lambdas, ARRAY_SIZE(lambdas), state->mold.sys,
                                                                (size_t)es.excited_state_idx, LAMBDA_CUTOFF);
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
        update_rep |= draw_surface_coloring(es, state->mold.sys, !density_property, mesh_mode);
        if (es.coloring == SurfaceColoring::AtomColors) {
            update_rep |= draw_representation_coloring(state, rep);
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
// The Representations window
// ---------------------------------------------------------------------------------------------------------

#define REPRESENTATION_NODE_DND "VIAMD_REPRESENTATION_NODE"

// A representation or a group, by what names it. Also the drag and drop payload.
struct NodeRef {
    uint32_t id = 0;        // 0 for none
    bool group = false;
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

static void apply_window_actions(ApplicationState* state, const WindowActions& act) {
    if (act.remove_all) {
        remove_all_representations(state);
        return;
    }
    if (act.move.node.id) {
        representation_move(state, act.move.node.group, act.move.node.id, act.move.target.group, act.move.target.id, act.move.drop);
    }
    if (act.clone.id) {
        if (act.clone.group) {
            clone_representation_group(state, act.clone.id);
        } else if (const Representation* rep = find_representation(state, act.clone.id)) {
            clone_representation(state, *rep);
        }
    }
    if (act.remove.id) {
        if (act.remove.group) {
            remove_representation_group(state, act.remove.id);
        } else {
            const int idx = find_rep_idx(*state, act.remove.id);
            if (idx >= 0) remove_representation(state, (size_t)idx);
        }
    }
    if (act.wrap.id) {
        wrap_in_group(state, act.wrap.group, act.wrap.id);
    }
    if (act.ungroup) {
        ungroup_representation_group(state, act.ungroup);
    }
    if (act.add_representation) {
        create_representation(state, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("all"), act.add_into);
    }
    if (act.add_group) {
        create_representation_group(state, STR_LIT("group"), act.add_into);
    }
}

// The menu of a row, or of the empty space below them (node.id 0)
static void draw_tree_menu(NodeRef node, WindowActions* act) {
    const bool group = node.id && node.group;
    if (node.id == 0 || group) {
        if (ImGui::MenuItem(ICON_FA_PLUS " Add representation")) {
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
    if (ImGui::MenuItem(ICON_FA_COPY " Duplicate")) {
        act->clone = node;
    }
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
static void draw_drop_target(const ApplicationState& state, NodeRef node, bool open, const ImRect& row, float x0, WindowActions* act) {
    if (!ImGui::BeginDragDropTargetCustom(row, ImGui::GetID("##drop"))) return;

    const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(REPRESENTATION_NODE_DND, ImGuiDragDropFlags_AcceptBeforeDelivery | ImGuiDragDropFlags_AcceptNoDrawDefaultRect);
    if (payload && payload->DataSize == (int)sizeof(NodeRef)) {
        NodeRef src;
        MEMCPY(&src, payload->Data, sizeof(src));

        const float t = (ImGui::GetMousePos().y - row.Min.y) / MAX(row.GetHeight(), 1.0f);
        const RepresentationDrop drop = drop_zone(t, node.group, open);

        if (representation_can_move(state, src.group, src.id, node.group, node.id, drop)) {
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
static void refresh_shown_again(ApplicationState* state, NodeRef node) {
    md_temp_scope_t temp = md_temp_begin();
    TreeNodes t = tree_nodes(*state, md_temp_allocator(temp));
    for (size_t k = 0; k < t.count; ++k) {
        const bool self = t.nodes[k]->group == node.group && t.nodes[k]->id == node.id;
        if (!self && !(node.group && representation_tree_is_within(t.nodes, t.count, k, node.id))) continue;
        if (k < t.num_reps) {
            Representation& rep = state->representation.reps[k];
            rep.needs_update = true;
            // Not one whose auto-update is off: what it picks was left where it is on purpose
            if (rep.dynamic_evaluation && representation_is_dynamic(*state, rep)) rep.filt_is_dirty = true;
        } else {
            RepresentationGroup& g = state->representation.groups[k - t.num_reps];
            if (g.filter.result.is_dynamic) g.filt_is_dirty = true;
        }
    }
    md_temp_end(temp);
}

// The buttons at the end of a row: show/hide, duplicate, remove. Placed from the right edge.
static void draw_row_buttons(const WindowTree& tree, bool* enabled, float right_x, NodeRef node, ApplicationState* state, WindowActions* act) {
    const ImVec2 btn_size = {tree.button_size, tree.button_size};
    const float total = (tree.button_size + tree.button_spacing) * 3;
    ImGui::SameLine(right_x - total, tree.button_spacing);
    if (ImGui::Button(*enabled ? ICON_FA_EYE : ICON_FA_EYE_SLASH, btn_size)) {
        *enabled = !*enabled;
        state->representation.atom_visibility_mask_dirty = true;
        if (*enabled) {
            refresh_shown_again(state, node);
        }
    }
    ImGui::SetItemTooltip(node.group ? "Show/Hide the group and everything in it" : "Show/Hide");
    ImGui::SameLine(0, tree.button_spacing);
    if (ImGui::Button(ICON_FA_COPY, btn_size)) {
        act->clone = node;
    }
    ImGui::SetItemTooltip(node.group ? "Duplicate the group and everything in it" : "Duplicate");
    ImGui::SameLine(0, tree.button_spacing);
    if (ImGui::DeleteButton(ICON_FA_XMARK, btn_size)) {
        act->remove = node;
    }
    ImGui::SetItemTooltip(node.group ? "Remove the group and everything in it" : "Remove");
}

static void draw_tree_level(ApplicationState* state, const WindowTree& tree, uint32_t parent, size_t depth, WindowActions* act);

// How the atoms of a representation are coloured, as rows of its Color section: the mapping and what it takes - a
// colour, an attribute through a colour scale, the colours of the secondary structures - and for sticks the bonds.
// Also what an electronic structure coloured by its atoms takes its colours from.
bool draw_representation_coloring(ApplicationState* state, Representation& rep) {
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
    if (rep.color_mapping == ColorMapping::CompType) {
        ImGui::SetItemTooltip("The colors of the residue types, set in the Residue Types of the System window");
    }

    if (rep.color_mapping == ColorMapping::Uniform) {
        prop::row("color");
        update_rep |= prop::color("##base_color", rep.base_color.elem);
    } else if (rep.color_mapping == ColorMapping::Attribute) {
        // The list of per atom fields IS the system's attribute table under atom/.
        // Queried here rather than cached anywhere, so it cannot disagree with the data.
        md_attribute_id_t prop_ids[64];
        size_t num_props = MIN(atom_attribute_query(prop_ids, ARRAY_SIZE(prop_ids), state->mold.sys), ARRAY_SIZE(prop_ids));

        const md_attributes_t& attributes = state->mold.sys.attributes;
        const md_attribute_t* selected_prop = md_attributes_get(&attributes, rep.atom_attribute.key);

        // A key the table no longer holds - a reload which dropped that field - falls
        // back to the first available rather than leaving the representation blank.
        if (!selected_prop && num_props > 0) {
            atom_attribute_select(&rep.atom_attribute, prop_ids[0], state->mold.sys);
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
                        atom_attribute_select(&rep.atom_attribute, prop_ids[i], state->mold.sys);
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

    if (rep.type == RepresentationType::Licorice || rep.type == RepresentationType::BallAndStick) {
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

    // Of md_gl's representations: below 1 they are drawn with the isosurface meshes' transparency
    if (RepresentationType::SpaceFill <= rep.type && rep.type <= RepresentationType::Cartoon) {
        prop::row("opacity", "How much of what lies behind it it hides: 1 opaque. Below that it is drawn transparent, with "
                             "the isosurfaces, and is then not picked and casts no ambient occlusion or shadows");
        update_rep |= ImGui::SliderFloat("##opacity", &rep.opacity, 0.0f, 1.0f, "%.2f", ImGuiSliderFlags_AlwaysClamp);
        // Only while transparent
        ImGui::BeginDisabled(rep.opacity >= 1.0f);
        prop::row("refraction", "Index of refraction: how much of the light the transparent surface reflects, the more toward its\n"
                                "silhouette, which it outlines (Fresnel). 1: none, an even veil. 1.33 water, 1.5 glass (the isosurfaces'),\n"
                                "higher a brighter rim and more of the surroundings mirrored.");
        update_rep |= ImGui::SliderFloat("##ior", &rep.ior, 1.0f, 3.0f, "%.2f", ImGuiSliderFlags_AlwaysClamp);
        ImGui::EndDisabled();
    }

    return update_rep;
}

// Whether a chain of the system is broken somewhere (MD_COMPONENT_FLAG_BREAK_AFTER): residues missing in a crystal
// structure, which the representations with chain gaps draw a dashed line over
static bool system_has_chain_breaks(const md_system_t& sys) {
    for (size_t i = 0; i < sys.component.count; ++i) {
        if (md_component_flags(&sys.component, i) & MD_COMPONENT_FLAG_BREAK_AFTER) return true;
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
            result = (sys.bond.flags[i] & MD_BOND_FLAG_ELASTIC) != 0;
        }
    }
    return result;
}

// The sizes of what a representation draws, as rows of its Geometry section: scales of the defaults, so 1 is as is
bool draw_representation_geometry(ApplicationState* state, Representation& rep) {
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
        if (cartoon_glyphs_valid(state->mold.cartoon_glyphs)) {
            prop::row("lipid, nucleic", "Thickness of the heads and tails of the lipids and of the backbones and bases of the nucleic acids");
            update_rep |= ImGui::SliderFloat("##lipid_nucleic", &rep.scale[3], 0.1f, 3.f, scale_fmt);
            if (state->mold.sys.nucleic_backbone.range.count > 0) {
                prop::row("nucleic backbone", "Radius of the band along the backbone of the nucleic acids: half its thickness, its width stays (Å, times the lipid, nucleic scale)");
                update_rep |= ImGui::SliderFloat("##nucleic_backbone", &rep.nucleic_band_radius, 0.1f, 1.5f, "%.2f Å", ImGuiSliderFlags_AlwaysClamp);
            }
            if (state->mold.cartoon_glyphs.num_slabs > 0) {
                prop::row("nucleic rings", "Radius of the rims of the sugar and base rings, and of the stick from each sugar to its base (Å, times the lipid, nucleic scale)");
                update_rep |= ImGui::SliderFloat("##nucleic_rings", &rep.nucleic_ring_radius, 0.1f, 1.0f, "%.2f Å", ImGuiSliderFlags_AlwaysClamp);
            }
            if (state->mold.cartoon_glyphs.num_split_heads > 0) {
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
    if (gaps && system_has_chain_breaks(state->mold.sys)) {
        prop::row("chain gaps", "A dashed line over each break of a chain, where residues are missing (or its pieces are not bonded)");
        update_rep |= ImGui::Checkbox("##chain_gaps", &rep.chain_gaps);
    }
    if (rep.type == RepresentationType::Licorice || rep.type == RepresentationType::BallAndStick) {
        prop::row("bond orders", "Double and triple bonds as two and three sticks, aromatic rings with a circle (dashed in fused aromatic systems, whose sextets are not assigned to particular rings), delocalized bonds with a dashed line");
        update_rep |= ImGui::Checkbox("##bond_orders", &rep.bond_orders);
        if (system_has_elastic_bonds(state->mold.sys)) {
            prop::row("elastic bonds", "The springs of the force field which are not bonds of the molecules: the elastic network of a coarse grained protein or nucleic acid, restraints. Drawn as plain bonds when shown");
            update_rep |= ImGui::Checkbox("##elastic_bonds", &rep.elastic_bonds);
        }
    }
    return update_rep;
}

bool draw_representation_type(ApplicationState* state, Representation& rep, bool atoms_only) {
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
                if (es_source_mask(state->mold.sys) == 0) continue;
            }
            if (ImGui::Selectable(representation_type_str[(int)i], i == (int)rep.type)) {
                rep.type = (RepresentationType)i;
                if (rep.type == RepresentationType::ElectronicStructure) {
                    // The default source is orbitals, which this system may not have
                    electronic_structure_select_available_source(&rep.electronic_structure, state->mold.sys);
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
static bool draw_representation_properties(ApplicationState* state, Representation& rep) {
    bool update_rep = false;
    prop::begin("representation");

    update_rep |= draw_representation_type(state, rep, false);

    // Of what md_gl draws and the box, where there is a cell to repeat in. Drawn so every frame: not an update.
    const bool periodic = ((RepresentationType::SpaceFill <= rep.type && rep.type <= RepresentationType::Cartoon) || rep.type == RepresentationType::SimulationBox) &&
                          md_unitcell_flags(&state->mold.state.unitcell) != 0;
    const bool atoms = RepresentationType::SpaceFill <= rep.type && rep.type <= RepresentationType::Cartoon;

    if (atoms) {
        // Completed within the context it is evaluated in: what is not there would not match
        const RepresentationContext ctx = representation_context(*state, rep.node.parent);
        prop::row("filter");
        if (filter_input("##filter", &rep.filter, representation_filter_context(*state, ctx))) {
            rep.filt_is_dirty = true;
            update_rep = true;
        }
        if (ctx.group && rep.filter.result.success && ImGui::IsItemHovered() && !filter_completions_open()) {
            ImGui::SetTooltip("Within the context of the group '%s': %s", ctx.group->name, ctx.group->filter.text);
        }

        if (periodic) {
            prop::row("periodic images", "Repeated in the periodic images of the view, as many as Visuals > View > Periodic Images shows");
            ImGui::Checkbox("##periodic_images", &rep.periodic_images);
        }
    }

    // Shown only where there is something to follow; the choice is kept either way, so that it holds
    // again for a filter edited back into depending on the frame
    if (representation_uses_atom_colors(rep) && (representation_is_dynamic(*state, rep) || rep.color_mapping == ColorMapping::Attribute)) {
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
        update_rep |= draw_representations_window_electronic_structure(state, rep);
        break;
    case RepresentationType::DipoleMoment:
        if (prop::section("Dipole")) {
            // The combo lists GROUPS, which is what the file publishes. A group holding one
            // moment per excited state is then picked apart by the index widget below, rather
            // than flattening every (group, element) pair into one list nobody can scan.
            DipoleGroup groups[16];
            const size_t num_groups = MIN(dipole_groups_gather(groups, ARRAY_SIZE(groups), state->mold.sys), ARRAY_SIZE(groups));

            // Asked of the key the representation holds, not searched for in the list: a key
            // that names nothing - a workspace pointed at a file without that group - simply
            // leaves the preview empty and the index widget away.
            DipoleGroup selected = {};
            const bool has_selection = dipole_group_from_key(&selected, state->mold.sys, rep.dipole.dipole_key);

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
                prop::row("periodic images", "Repeated in the periodic images of the view, as many as Visuals > View > Periodic Images shows");
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
            update_rep |= draw_representation_coloring(state, rep);
            prop::section_end();
        }
        if (prop::section("Geometry", false)) {
            update_rep |= draw_representation_geometry(state, rep);
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

static void draw_representation_row(ApplicationState* state, const WindowTree& tree, Representation& rep, WindowActions* act) {
    const NodeRef node = {rep.node.id, false};
    ImGui::PushID("rep");
    ImGui::PushID((int)rep.node.id);

    const RowLayout row = row_begin();
    const float pad = 3.0f;

    const DragHandle handle = draw_drag_handle(tree, node, rep.name);
    draw_handle_icon(handle, ICON_FA_GRIP_VERTICAL, ImGui::GetColorU32(handle.hot ? ImGuiCol_Text : ImGuiCol_TextDisabled));
    if (handle.tooltip) ImGui::SetTooltip(DRAG_HANDLE_HINT);
    ImGui::SameLine(0, 0);

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

    draw_row_buttons(tree, &rep.enabled, row.right_x, node, state, act);

    const ImRect row_rect = {row.min, {row.right_screen, ImGui::GetItemRectMax().y}};
    draw_drop_target(*state, node, open, row_rect, row.min.x, act);

    if (ImGui::IsMouseReleased(ImGuiMouseButton_Right) && ImGui::IsWindowHovered() && row_rect.Contains(ImGui::GetMousePos())) {
        ImGui::OpenPopup("##menu");
    }
    if (ImGui::BeginPopup("##menu")) {
        draw_tree_menu(node, act);
        ImGui::EndPopup();
    }

    if (open) {
        if (draw_representation_properties(state, rep)) {
            flag_representation_as_dirty(&rep);
        }
        ImGui::Spacing();
        ImGui::TreePop();
    }

    ImGui::PopID();
    ImGui::PopID();
}

// The filter of a group, the first line within it, a row of the same grid as the properties of the representations.
// Completed within the context the group is in, as it is evaluated.
static void draw_group_filter(ApplicationState* state, RepresentationGroup& group) {
    prop::begin("representation");
    prop::row("context");

    const RepresentationContext outer = representation_context(*state, group.node.parent);
    const bool edited = filter_input("##context", &group.filter, representation_filter_context(*state, outer), 0, nullptr, "none: what is within sees the whole system");

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
        flag_within_as_dirty(state, true, group.node.id);
    }
}

static void draw_group_row(ApplicationState* state, const WindowTree& tree, RepresentationGroup& group, size_t depth, WindowActions* act) {
    const NodeRef node = {group.node.id, true};
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

    draw_row_buttons(tree, &group.enabled, row.right_x, node, state, act);

    const ImRect row_rect = {row.min, {row.right_screen, ImGui::GetItemRectMax().y}};
    draw_drop_target(*state, node, open, row_rect, row.min.x, act);

    if (ImGui::IsMouseReleased(ImGuiMouseButton_Right) && ImGui::IsWindowHovered() && row_rect.Contains(ImGui::GetMousePos())) {
        ImGui::OpenPopup("##menu");
    }
    if (ImGui::BeginPopup("##menu")) {
        draw_tree_menu(node, act);
        ImGui::EndPopup();
    }

    if (open) {
        draw_group_filter(state, group);

        // What is within a hidden group is drawn dimmed: it is hidden too, whatever its own eye says
        const bool dim = !group.enabled;
        if (dim) ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
        draw_tree_level(state, tree, group.node.id, depth + 1, act);
        if (dim) ImGui::PopStyleVar();

        ImGui::TreePop();
    }

    ImGui::PopID();
    ImGui::PopID();
}

static void draw_tree_level(ApplicationState* state, const WindowTree& tree, uint32_t parent, size_t depth, WindowActions* act) {
    const TreeNodes& t = tree.t;
    if (depth > t.count) return;    // A cycle; the tree breaks those when it is normalized

    int* kids = (int*)md_alloc(state->allocator.frame, sizeof(int) * MAX(t.count, (size_t)1));
    const size_t n = representation_tree_children(kids, t.count, t.nodes, t.count, parent);
    for (size_t i = 0; i < n; ++i) {
        const size_t k = (size_t)kids[i];
        if (k < t.num_reps) {
            draw_representation_row(state, tree, state->representation.reps[k], act);
        } else {
            draw_group_row(state, tree, state->representation.groups[k - t.num_reps], depth, act);
        }
    }
}

void draw_representations_window(ApplicationState* state) {
    if (!state->representation.show_window) return;

    ImGui::SetNextWindowSize({300,200}, ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Representations", &state->representation.show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
        ImGui::End();
        return;
    }

    WindowActions act;

    if (ImGui::Button(ICON_FA_FOLDER_PLUS)) {
        act.add_group = true;
    }
    ImGui::SetItemTooltip("Add a group. What is in it is shown and hidden with it, and its filter is the\n"
                          "context the filters within it are evaluated in.");
    ImGui::SameLine();
    if (ImGui::Button(ICON_FA_PLUS)) {
        act.add_representation = true;
    }
    ImGui::SetItemTooltip("Add a representation");

    // Removing everything at the far end, behind a question: it cannot be taken back
    ImGui::SameLine();
    const float trash_w = ImGui::CalcTextSize(ICON_FA_TRASH_CAN).x + 2.0f * ImGui::GetStyle().FramePadding.x;
    ImGui::SetCursorPosX(ImGui::GetCursorPosX() + MAX(ImGui::GetContentRegionAvail().x - trash_w, 0.0f));
    ImGui::BeginDisabled(md_array_size(state->representation.reps) == 0 && md_array_size(state->representation.groups) == 0);
    if (ImGui::Button(ICON_FA_TRASH_CAN)) {
        ImGui::OpenPopup("##remove_all");
    }
    ImGui::EndDisabled();
    ImGui::SetItemTooltip("Remove every representation and group");
    if (ImGui::BeginPopup("##remove_all")) {
        ImGui::TextUnformatted("Remove every representation and group?");
        if (ImGui::DeleteButton("remove all")) {
            act.remove_all = true;
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("cancel")) {
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }

    ImGui::Spacing();
    ImGui::Separator();

    // The arrays are only read while the rows are drawn - every change to them waits in 'act' - so the
    // nodes gathered here hold until the end
    WindowTree tree;
    tree.t = tree_nodes(*state, state->allocator.frame);
    tree.button_size = ImGui::GetFontSize() + 3.0f * 2;
    tree.button_spacing = 2.0f;
    tree.handle_width = MAX(ImGui::CalcTextSize(ICON_FA_GRIP_VERTICAL).x,
                        MAX(ImGui::CalcTextSize(ICON_FA_FOLDER).x, ImGui::CalcTextSize(ICON_FA_FOLDER_OPEN).x)) + 6.0f;

    if (tree.t.count == 0) {
        ImGui::TextDisabled("Nothing to show. Add a representation with " ICON_FA_PLUS ".");
    }
    draw_tree_level(state, tree, 0, 0, &act);

    // The space below the rows: dropping there puts it last at the top level, and a right click adds
    const ImVec2 rest = {MAX(ImGui::GetContentRegionAvail().x, 1.0f), MAX(ImGui::GetContentRegionAvail().y, ImGui::GetFrameHeight())};
    ImGui::InvisibleButton("##below", rest);
    if (ImGui::BeginDragDropTarget()) {
        const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(REPRESENTATION_NODE_DND, ImGuiDragDropFlags_AcceptBeforeDelivery | ImGuiDragDropFlags_AcceptNoDrawDefaultRect);
        if (payload && payload->DataSize == (int)sizeof(NodeRef)) {
            NodeRef src;
            MEMCPY(&src, payload->Data, sizeof(src));
            if (representation_can_move(*state, src.group, src.id, false, 0, RepresentationDrop_IntoLast)) {
                const ImVec2 min = ImGui::GetItemRectMin();
                ImGui::GetWindowDrawList()->AddLine(min, {ImGui::GetItemRectMax().x, min.y}, ImGui::GetColorU32(ImGuiCol_DragDropTarget), 2.0f);
                if (payload->IsDelivery()) {
                    act.move.node   = src;
                    act.move.target = NodeRef{};
                    act.move.drop   = RepresentationDrop_IntoLast;
                }
            }
        }
        ImGui::EndDragDropTarget();
    }
    if (ImGui::BeginPopupContextItem("##below_menu")) {
        draw_tree_menu(NodeRef{}, &act);
        ImGui::EndPopup();
    }

    ImGui::End();

    apply_window_actions(state, act);
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
        if (!representation_shown(state, rep)) continue;

        if (rep.type == RepresentationType::ElectronicStructure) {
            ElectronicStructureRepresentation& es = rep.electronic_structure;
            const bool has_field = state.settings.isosurface_rendering == IsosurfaceRendering_Mesh ? electronic_structure_mesh_field_span(es).valid : es.field_vol.tex_id != 0;
            if (es.coloring == SurfaceColoring::Field && es.field_map.show_legend && has_field) {
                color_scale_draw_legend(es.field_map, surface_field_unit(es.field_kind), surface_field_kind_str[(int)es.field_kind], rep.name, (int)rep.node.id, slot++, &es.field_isolines);
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
