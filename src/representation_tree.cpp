#include <representation_tree.h>

#include <core/md_common.h>
#include <core/md_allocator.h>

int representation_tree_find(RepresentationNode* const nodes[], size_t count, bool group, uint32_t id) {
    if (id == 0) return -1;
    for (size_t i = 0; i < count; ++i) {
        if (nodes[i]->group == group && nodes[i]->id == id) return (int)i;
    }
    return -1;
}

// The group a node is in as far as the tree goes: one that names no group is at the top level
static uint32_t effective_parent(RepresentationNode* const nodes[], size_t count, const RepresentationNode* node) {
    return representation_tree_find(nodes, count, true, node->parent) >= 0 ? node->parent : 0;
}

bool representation_tree_is_within(RepresentationNode* const nodes[], size_t count, size_t idx, uint32_t group_id) {
    if (idx >= count || group_id == 0) return false;
    uint32_t p = nodes[idx]->parent;
    // Bounded by the number of nodes, so that a cycle - which normalize breaks, but which this does not
    // assume has been - cannot keep it going
    for (size_t steps = 0; p != 0 && steps < count; ++steps) {
        if (p == group_id) return true;
        const int g = representation_tree_find(nodes, count, true, p);
        if (g < 0) break;
        p = nodes[g]->parent;
    }
    return false;
}

size_t representation_tree_ancestors(int out[], size_t cap, RepresentationNode* const nodes[], size_t count, size_t idx) {
    if (idx >= count) return 0;
    size_t n = 0;
    uint32_t p = nodes[idx]->parent;
    for (size_t steps = 0; p != 0 && steps < count; ++steps) {
        const int g = representation_tree_find(nodes, count, true, p);
        if (g < 0) break;
        if (out && n < cap) out[n] = g;
        n += 1;
        p = nodes[g]->parent;
    }
    return n;
}

static bool sorts_before(RepresentationNode* const nodes[], int a, int b) {
    return nodes[a]->order < nodes[b]->order || (nodes[a]->order == nodes[b]->order && a < b);
}

size_t representation_tree_children(int out[], size_t cap, RepresentationNode* const nodes[], size_t count, uint32_t parent) {
    if (count == 0) return 0;

    md_temp_scope_t temp = md_temp_begin();
    int* kids = (int*)md_temp_alloc(temp, sizeof(int) * count);

    size_t n = 0;
    for (size_t i = 0; i < count; ++i) {
        if (effective_parent(nodes, count, nodes[i]) == parent) {
            // Insertion sort as they come: there are a handful
            size_t j = n++;
            while (j > 0 && sorts_before(nodes, (int)i, kids[j - 1])) {
                kids[j] = kids[j - 1];
                j -= 1;
            }
            kids[j] = (int)i;
        }
    }
    for (size_t i = 0; i < n && i < cap; ++i) {
        out[i] = kids[i];
    }

    md_temp_end(temp);
    return n;
}

// Every level reads the orders of its own siblings, all of them, before it writes any: a sibling is
// numbered only once the whole level has been sorted, and what is below it is only read after that.
static void renumber(RepresentationNode* nodes[], size_t count, uint32_t parent, int32_t* next, size_t depth) {
    if (depth > count) return;  // A cycle. normalize breaks them before it gets here.

    md_temp_scope_t temp = md_temp_begin();
    int* kids = (int*)md_temp_alloc(temp, sizeof(int) * count);
    const size_t n = representation_tree_children(kids, count, nodes, count, parent);
    for (size_t i = 0; i < n; ++i) {
        RepresentationNode* node = nodes[kids[i]];
        node->order = *next;
        *next += 2;
        if (node->group) {
            renumber(nodes, count, node->id, next, depth + 1);
        }
    }
    md_temp_end(temp);
}

void representation_tree_normalize(RepresentationNode* nodes[], size_t count) {
    for (size_t i = 0; i < count; ++i) {
        nodes[i]->parent = effective_parent(nodes, count, nodes[i]);
    }

    // A group that is its own ancestor: the walk up from it comes back to it. The one found first is
    // taken out to the top level, which puts the rest of its cycle below it.
    for (size_t i = 0; i < count; ++i) {
        if (!nodes[i]->group) continue;
        uint32_t p = nodes[i]->parent;
        for (size_t steps = 0; p != 0 && steps <= count; ++steps) {
            if (p == nodes[i]->id) {
                nodes[i]->parent = 0;
                break;
            }
            const int g = representation_tree_find(nodes, count, true, p);
            if (g < 0) break;
            p = nodes[g]->parent;
        }
    }

    int32_t next = 0;
    renumber(nodes, count, 0, &next, 0);
}

// Where a node ends up when dropped on 'target': the group it goes into
static uint32_t drop_parent(RepresentationNode* const nodes[], size_t count, int target, RepresentationDrop drop) {
    if (target < 0) return 0;
    const RepresentationNode* t = nodes[target];
    return (drop == RepresentationDrop_IntoFirst || drop == RepresentationDrop_IntoLast) ? t->id : effective_parent(nodes, count, t);
}

bool representation_tree_can_move(RepresentationNode* const nodes[], size_t count, size_t idx, int target, RepresentationDrop drop) {
    if (idx >= count || target >= (int)count) return false;
    if (target < 0) {
        return drop == RepresentationDrop_IntoLast;
    }
    if ((size_t)target == idx) return false;

    const bool into = drop == RepresentationDrop_IntoFirst || drop == RepresentationDrop_IntoLast;
    if (into && !nodes[target]->group) return false;

    const RepresentationNode* node = nodes[idx];
    if (node->group) {
        const uint32_t parent = drop_parent(nodes, count, target, drop);
        if (parent == node->id) return false;
        const int g = representation_tree_find(nodes, count, true, parent);
        if (g >= 0 && representation_tree_is_within(nodes, count, (size_t)g, node->id)) return false;
    }
    return true;
}

bool representation_tree_move(RepresentationNode* nodes[], size_t count, size_t idx, int target, RepresentationDrop drop) {
    if (!representation_tree_can_move(nodes, count, idx, target, drop)) return false;

    // The numbers picked below rely on the orders being two apart
    representation_tree_normalize(nodes, count);

    RepresentationNode* node = nodes[idx];
    const uint32_t parent = drop_parent(nodes, count, target, drop);
    int32_t order = INT32_MAX;
    if (target >= 0) {
        const int32_t t = nodes[target]->order;
        switch (drop) {
        case RepresentationDrop_Before:    order = t - 1; break;
        case RepresentationDrop_After:     order = t + 1; break;
        case RepresentationDrop_IntoFirst: order = t + 1; break;
        case RepresentationDrop_IntoLast:  order = INT32_MAX; break;
        }
    }
    node->parent = parent;
    node->order  = order;

    representation_tree_normalize(nodes, count);
    return true;
}
