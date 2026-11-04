#pragma once

// LAYER GROUPS
//
// The layers of the Ramachandran and Shapespace windows can be put in groups, as representations are, so that they show
// one subset of the system without each repeating it. A group's subset - a filter, or the selection - is the context of
// the layers in it: a layer's own filter is evaluated within the atoms it picks, as the left hand side of 'in' is within
// its right in a script (residue(1) is the group's first residue), and a layer without a filter shows the whole of it.
// What a layer is over in time, and how it is drawn, stay the layer's: a group is typically the same residues over
// several periods. A group without a subset only groups. One level: a group is not put in a group.
//
// A group also shows or hides what it holds, and keeps it together in the list. Where the groups and layers are in the
// list is a tree (representation_tree.h), its nodes on the layers and groups themselves. What is drawn first is what is
// listed first.

#include <representation_tree.h>
#include <filter_input.h>

#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_str.h>

#include <imgui.h>
#include <imgui_internal.h>

#include <stddef.h>
#include <stdint.h>

struct ApplicationState;

namespace viamd {
struct serialization_state_t;
}

constexpr int LAYER_MAX_GROUPS = 16;

struct LayerGroup {
    RepresentationNode node = {0, 0, 0, true};
    char  name[64] = "";            // Empty: named after its subset (layer_group_label)
    bool  visible = true;           // What it holds is shown only while it is
    bool  open = true;              // Listed with what it holds

    AtomFilter filter = {};         // The subset. Empty, and not following the selection: no context of its own.
    bool  follow_selection = false; // The selection instead of the filter

    // ### Derived, not saved
    md_bitset_t atoms = {};         // The context: what the filter picks, or the selection
    md_array(md_bitset_t) sets = 0; // The structures the filter gives, when they are asked for (layer_group_update)
    uint64_t key = 0;               // Of what they were found from
    uint64_t hash = 0;              // Of what was found: what is within it is found again when this changes
    bool  valid = true;             // The filter evaluated (or there is none)
};

static inline bool layer_group_has_context(const LayerGroup& g) {
    return g.follow_selection || g.filter.text[0] != '\0';
}

// The context found again when what it is of changes: the filter, the script's names it may use, the selection it
// follows, the system (generation). follow_frame: a filter that depends on the coordinates follows the frame shown;
// otherwise it is taken in the frame it is evaluated in, as a density over frames could not follow it. want_sets: the
// structures of the filter too, one set each (a group following the selection has one: all of it).
void layer_group_update(LayerGroup* g, const ApplicationState* app, uint64_t generation, bool follow_frame, bool want_sets);

void layer_group_free(LayerGroup* g);

// What is within g is evaluated in: the system as shown, and g's atoms when it has a subset. g may be NULL.
script::FilterCtx layer_group_filter_context(const ApplicationState* app, const LayerGroup* g);

// Its name, or what it is of
void layer_group_label(char* buf, size_t cap, const LayerGroup& g);

// The filter of a layer taken out of a group, so that it shows what it did: its own within the group's ("(a) in union(b)"),
// or the group's for one without. Unchanged for a group without a filter. Returns false when it does not fit.
bool layer_group_compose_filter(char* layer_filter, size_t cap, const LayerGroup& g);

// Whether a layer taken out of g can be given a filter that shows what it does within g (layer_group_compose_filter):
// not one that follows the selection within a filter, nor one of its own within a group following the selection or
// within a filter that depends on the coordinates ('in' takes its right side when it compiles)
bool layer_group_keeps_filter(const char* layer_filter, bool layer_follows_selection, const LayerGroup& g);

// ## The tree
//
// The nodes of a window's layers and groups, as representation_tree takes them. One level: what puts a group within a
// group is refused.

bool layer_tree_can_move(RepresentationNode* const nodes[], size_t count, size_t idx, int target, RepresentationDrop drop);
bool layer_tree_move(RepresentationNode* nodes[], size_t count, size_t idx, int target, RepresentationDrop drop);

// The order of a new node: first within the group 'parent' (or at the top), or last
int32_t layer_tree_first_order(RepresentationNode* const nodes[], size_t count, uint32_t parent);
constexpr int32_t LAYER_TREE_LAST = INT32_MAX;

// ## UI

// Dragged between the rows of one window's list: 'window' tells the lists apart
struct LayerNodeRef {
    uint32_t window = 0;
    uint32_t id = 0;
    bool     group = false;
};
#define LAYER_NODE_DND "LAYER_NODE"

// The last item, a row of the list, as something to drag
void layer_tree_drag_source(LayerNodeRef ref, const char* label);

// The row as a place to drop a layer or group of the same window on: above it, below it, or into it for a group (open:
// whether it is listed with what it holds). x0 is where the row's content starts. Returns true on the drop, with what
// was dropped and where.
bool layer_tree_drop_target(const ImRect& row, float x0, LayerNodeRef target, bool open, RepresentationNode* const nodes[], size_t count,
                            LayerNodeRef* dropped, RepresentationDrop* drop);

// A place to drop at the end of the group 'parent' (0: of the list): its empty row, or the space below the list.
// Returns true on the drop, with what was dropped; it goes last in 'parent'.
bool layer_tree_drop_target_end(const ImRect& rect, float x0, uint32_t window, uint32_t parent, RepresentationNode* const nodes[], size_t count,
                                LayerNodeRef* dropped);

enum LayerGroupAction : int {
    LayerGroupAction_None = 0,
    LayerGroupAction_Duplicate,
    LayerGroupAction_Ungroup,
    LayerGroupAction_Remove,
};

struct LayerGroupRow {
    bool clicked = false;       // Selected
    bool hovered = false;
    LayerGroupAction action = LayerGroupAction_None;
    ImRect rect;                // Of the row, to drop on
    float  x0 = 0.0f;           // Where its content starts
};

// A group's row: open or not, shown or not, its name and what it is of (summary, one line). Its menu (duplicate,
// ungroup, remove) is the window's to act on. ungroup_warning, if given: what ungrouping will not keep.
LayerGroupRow layer_group_row(LayerGroup& g, uint32_t window, bool selected, const char* summary, bool invalid, const char* ungroup_warning = nullptr);

// The inspector's rows of a group: its name and its subset. 'what' names what the window's layers are of ("residues",
// "structures"). Returns true when the subset changed.
bool layer_group_inspector(LayerGroup& g, const ApplicationState* app, const char* what);

// ## Workspace

// Its section's entries; id: what the layers in it name it by (Group=id)
void layer_group_serialize(viamd::serialization_state_t& state, str_t section, const LayerGroup& g);
// One entry of a group's section. Its file id is written to *file_id.
void layer_group_deserialize_entry(LayerGroup* g, str_t ident, str_t arg, int* file_id);
