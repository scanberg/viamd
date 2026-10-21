#pragma once

// The representations - how the atoms, the electronic structure and the dipoles are drawn - and the groups
// they can be put in: creating and removing them, where they are in the tree of groups, evaluating them,
// keeping them in a workspace, and the Representations window that edits all of it.
//
// The tree itself is plain data on each representation and group (representation_tree.h). Two rules hang
// off it:
//
//   shown    A representation is shown when it is enabled and so is every group around it.
//   context  A group's filter, when it has one, is the context of everything within it: their filters are
//            evaluated within the atoms the group picks (md_filter_context_t::atoms), and so is the filter of a
//            group within it. The context of a representation is that of the innermost group around it that
//            has a filter; without one it is the whole system. A filter is completed within its context too.

#include <viamd.h>

namespace viamd {
struct serialization_state_t;
struct deserialization_state_t;
}

// --- Representations ---

// Appended at the end of 'group' (0 for the top level)
Representation* create_representation(ApplicationState* app, RepresentationType type = RepresentationType::SpaceFill, ColorMapping color_mapping = ColorMapping::Type, str_t filter = STR_LIT("all"), uint32_t group = 0);
// Placed just after the original, in the same group
Representation* clone_representation(ApplicationState* app, const Representation& rep);
void remove_representation(ApplicationState* app, size_t idx);
// The representations and the groups
void remove_all_representations(ApplicationState* app);
void create_default_representations(ApplicationState* app);

Representation* find_representation(ApplicationState* app, uint32_t id);

// --- Groups ---

// Appended at the end of 'parent' (0 for the top level)
RepresentationGroup* create_representation_group(ApplicationState* app, str_t name = STR_LIT("group"), uint32_t parent = 0);
// The group and everything within it, placed just after the original
void clone_representation_group(ApplicationState* app, uint32_t id);
// The group and everything within it
void remove_representation_group(ApplicationState* app, uint32_t id);
// The group goes, and what it held takes its place, in the group around it
void ungroup_representation_group(ApplicationState* app, uint32_t id);

RepresentationGroup* find_representation_group(ApplicationState* app, uint32_t id);
const RepresentationGroup* find_representation_group(const ApplicationState& app, uint32_t id);

// --- Tree ---

bool representation_shown(const ApplicationState& app, const Representation& rep);
bool representation_group_shown(const ApplicationState& app, const RepresentationGroup& group);

// The context something within the group 'parent' is evaluated in: the innermost group out from it that has
// a filter, or none (atoms NULL) for the whole system
struct RepresentationContext {
    const RepresentationGroup* group = nullptr;     // whose filter it is, NULL for the whole system
    const md_bitset_t* atoms = nullptr;             // NULL for the whole system
    bool valid = true;                              // false while that group's filter, or one around it, is not
    bool dynamic = false;                           // it depends on the frame
};
RepresentationContext representation_context(const ApplicationState& app, uint32_t parent);
// What a filter in that context is evaluated and completed against: the application's (filter_context), within its atoms
md_filter_context_t representation_filter_context(const ApplicationState& app, const RepresentationContext& ctx);

// Moves a representation (group false) or a group, with everything within it, to where 'drop' puts it
// relative to the target - a representation or a group by id, or the top level (target_id 0, which takes only
// RepresentationDrop_IntoLast). False when it cannot go there: a group cannot go anywhere within itself.
bool representation_move(ApplicationState* app, bool group, uint32_t id, bool target_group, uint32_t target_id, RepresentationDrop drop);
bool representation_can_move(const ApplicationState& app, bool group, uint32_t id, bool target_group, uint32_t target_id, RepresentationDrop drop);

// --- Evaluation ---

// For a system that was just loaded: GPU resources made again, every filter evaluated again
void init_all_representations(ApplicationState* app);
// Once per frame: the contexts of the groups, then the representations, that need it and are shown
void update_all_representations(ApplicationState* app);

void flag_representation_as_dirty(Representation* rep);
// The representations and the contexts of the groups: what they are evaluated against changed (the script, the system)
void flag_all_representations_as_dirty(ApplicationState* app);
// The frame changed: what follows it is evaluated again. That is a filter or a context which depends on the
// frame (unless its auto-update is off), colours by secondary structure or by an attribute, and an
// electronic structure, whose volume moves with the atoms.
void flag_dynamic_representations_as_dirty(ApplicationState* app);

// Whether the atoms of a representation depend on the frame: its own filter, or the context it is in
bool representation_is_dynamic(const ApplicationState& app, const Representation& rep);
bool representation_uses_atom_colors(const Representation& rep);

// The union of the atoms of the shown representations, visibility_mask
void recompute_atom_visibility_mask(ApplicationState* app);

// --- Representations of atoms, wherever they are drawn ---
// What the Representations window and others (the reference structures of the Density Volume window) share of the
// representations md_gl draws, SpaceFill to Cartoon: whether one can be drawn, its draw op, its colours, its
// properties in the window and in a workspace.

// Whether the system has what a type draws: atoms, bonds, a backbone
bool representation_atom_type_valid(const ApplicationState& app, RepresentationType type);

// The draw op of a representation of atoms: its type, arguments and colours (md_rep). The copies and the subset are
// the caller's.
md_gl_draw_op_t representation_draw_op(const Representation& rep);
// The draw ops of the glyphs of a Cartoon, the lipids and nucleic acids (cartoon_glyphs.h), pushed to ops: their
// spheres, tubes and slabs, thickened by the w of its scale, in its colours of the glyphs (glyph_rep). None while it
// shows none or their positions are not derived (cartoon_glyphs_update_positions). The copies and the subset (of the
// glyphs) are the caller's.
void representation_glyph_draw_ops(md_array(md_gl_draw_op_t)* ops, const ApplicationState& app, const Representation& rep, md_allocator_i* alloc);

// The colours of all num_atoms atoms of the system for a representation which colours by atom
// (representation_uses_atom_colors): by its colour mapping, then towards its tint and saturation. 'shown' are the
// atoms a colour scale of an attribute spans the values of (NULL for all), whose span it keeps (atom_attribute).
void representation_atom_colors(uint32_t* colors, size_t num_atoms, ApplicationState* app, Representation* rep, const md_bitset_t* shown);

// --- Workspace ---

// [RepresentationGroup] and [Representation] sections, depth first: a group before what is within it
void representations_serialize(viamd::serialization_state_t& state, const ApplicationState* app);
// The entries of a representation in a section, all but where it is in the tree (Group), and one of them read: false
// when it is not one of them
void representation_write_entries(viamd::serialization_state_t& state, const ApplicationState* app, const Representation& rep);
bool representation_read_entry(Representation* rep, str_t ident, str_t arg);
// True when the section was one of the above, which it has then read
bool representations_deserialize(ApplicationState* app, str_t section, viamd::deserialization_state_t& state);

// --- UI ---

void draw_representations_window(ApplicationState* app);
// Rows of the properties of a representation, in a property grid (property_grid.h), as the Representations window
// draws them: its type (atoms_only: SpaceFill to Cartoon), the rows of its Color section, those of its Geometry
// section. Each returns whether it needs updating.
bool draw_representation_type(ApplicationState* app, Representation& rep, bool atoms_only);
bool draw_representation_coloring(ApplicationState* app, Representation& rep);
bool draw_representation_geometry(ApplicationState* app, Representation& rep);
// The legends of the shown representations that colour by a value. A field's legend also edits the isolines of
// its surfaces and sets the transient one (ColorScaleIsolines), so this is drawn before the frame is rendered.
void draw_representation_legends(ApplicationState& app);
