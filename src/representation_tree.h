#pragma once

// Where the representations sit in the Representations window: a tree whose inner nodes are groups.
//
// Each representation and each group carries its own place in it - the group it is in and an order
// among what else is in there - and nothing else describes the tree. There is no second structure to
// keep in step with the two arrays, and nothing to rebuild when one of them changes.
//
// 'order' is a sort key between siblings. representation_tree_normalize renumbers it to the position
// in the window, depth first and two apart, which keeps the siblings in order and turns every placement
// into picking a number between two others:
//
//     before X                  X.order - 1
//     after X                   X.order + 1    - and after what X holds: that sorts among X's own
//                                                children, which are not X's siblings
//     first within group G      G.order + 1    - G's first child is G.order + 2
//     last within group G       INT32_MAX
//
// It also means that a node lifted out of a group into the group around it lands where the group was,
// in the order it had: everything within G is numbered between G and G's next sibling.
//
// This unit is free of ImGui and of ApplicationState so that it can be tested; representation.h keeps
// the tree on the representations and groups and draws it.

#include <stdint.h>
#include <stddef.h>

struct RepresentationNode {
    uint32_t id     = 0;      // Never 0. Groups and representations are numbered apart, so a node is
                              // named by its id together with what it is.
    uint32_t parent = 0;      // The id of the group it is in, 0 at the top level
    int32_t  order  = 0;      // Where it is among its siblings, see above
    bool     group  = false;
};

enum RepresentationDrop : int {
    RepresentationDrop_Before,      // the target's sibling, just before it
    RepresentationDrop_After,       // the target's sibling, just after it and what it holds
    RepresentationDrop_IntoFirst,   // the first within the target group
    RepresentationDrop_IntoLast,    // the last within the target group, or at the top level
};

// The index of the node with this id, -1 for none
int representation_tree_find(RepresentationNode* const nodes[], size_t count, bool group, uint32_t id);

// Whether node idx is somewhere within the group 'group_id', at any depth. Nothing is within itself.
bool representation_tree_is_within(RepresentationNode* const nodes[], size_t count, size_t idx, uint32_t group_id);

// The groups around node idx, the innermost first. Writes at most cap and returns how many there are.
// A parent that names no group ends the walk there: the node is at the top level as far as it goes.
size_t representation_tree_ancestors(int out[], size_t cap, RepresentationNode* const nodes[], size_t count, size_t idx);

// What is directly in the group 'parent' (0 for the top level), in order. A node whose parent names no
// group is at the top level. Writes at most cap and returns how many there are.
size_t representation_tree_children(int out[], size_t cap, RepresentationNode* const nodes[], size_t count, uint32_t parent);

// Brings the tree to its canonical form: a parent that names no group, or that would put a group within
// itself, becomes the top level, and 'order' is renumbered to the depth first position, two apart. The
// relative order of siblings is kept, ties broken by their position in 'nodes'.
void representation_tree_normalize(RepresentationNode* nodes[], size_t count);

// Whether node idx can go where 'drop' puts it relative to 'target'. target -1 is the top level, which
// takes only RepresentationDrop_IntoLast. A group cannot go anywhere within itself, and nothing goes into
// something that is not a group.
bool representation_tree_can_move(RepresentationNode* const nodes[], size_t count, size_t idx, int target, RepresentationDrop drop);

// Moves node idx, with everything in it for a group, to where 'drop' puts it relative to 'target', and
// leaves the tree normalized. Returns false and changes nothing when it cannot go there.
bool representation_tree_move(RepresentationNode* nodes[], size_t count, size_t idx, int target, RepresentationDrop drop);
