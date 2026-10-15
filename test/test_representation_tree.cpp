#include "utest.h"

#include <representation_tree.h>

#include <core/md_common.h>

#include <stdio.h>
#include <string.h>

/* The Representations window's tree is nothing but (parent, order) on each node. These pin down what the
 * window relies on: siblings in order, every placement a drag can ask for, and that a group can never be
 * dragged into itself - which would detach it and everything in it from the top level. */

namespace {

struct Tree {
    RepresentationNode  storage[32];
    RepresentationNode* nodes[32];
    size_t count = 0;

    Tree() {
        for (size_t i = 0; i < 32; ++i) nodes[i] = &storage[i];
    }

    // Returns the index of the node, which is what the tree functions take
    size_t add(bool group, uint32_t id, uint32_t parent, int32_t order) {
        storage[count] = RepresentationNode{id, parent, order, group};
        return count++;
    }
    size_t rep(uint32_t id, uint32_t parent = 0, int32_t order = 0)   { return add(false, id, parent, order); }
    size_t group(uint32_t id, uint32_t parent = 0, int32_t order = 0) { return add(true,  id, parent, order); }

    int find_rep(uint32_t id)   const { return representation_tree_find(nodes, count, false, id); }
    int find_group(uint32_t id) const { return representation_tree_find(nodes, count, true,  id); }
};

// The tree as the window lists it: "r1 g1[r2 g2[r3]] r4"
static void describe_level(char* buf, size_t cap, size_t* len, const Tree& t, uint32_t parent, int depth) {
    int kids[32];
    const size_t n = representation_tree_children(kids, 32, t.nodes, t.count, parent);
    for (size_t i = 0; i < n && depth < 32; ++i) {
        const RepresentationNode* node = t.nodes[kids[i]];
        *len += (size_t)snprintf(buf + *len, cap - *len, "%s%c%u", i ? " " : "", node->group ? 'g' : 'r', node->id);
        if (node->group) {
            *len += (size_t)snprintf(buf + *len, cap - *len, "[");
            describe_level(buf, cap, len, t, node->id, depth + 1);
            *len += (size_t)snprintf(buf + *len, cap - *len, "]");
        }
    }
}

static const char* describe(const Tree& t) {
    static char buf[512];
    size_t len = 0;
    buf[0] = '\0';
    describe_level(buf, sizeof(buf), &len, t, 0, 0);
    return buf;
}

// The orders of the nodes in window order, which normalize makes 0, 2, 4, ...
static bool orders_are_depth_first(const Tree& t) {
    int32_t expected = 0;
    // Walk the window order by repeatedly taking the smallest order: after normalize that IS the window order
    bool used[32] = {};
    for (size_t k = 0; k < t.count; ++k) {
        int best = -1;
        for (size_t i = 0; i < t.count; ++i) {
            if (!used[i] && (best < 0 || t.nodes[i]->order < t.nodes[best]->order)) best = (int)i;
        }
        used[best] = true;
        if (t.nodes[best]->order != expected) return false;
        expected += 2;
    }
    return true;
}

}  // namespace

UTEST(representation_tree, siblings_by_order_ties_by_position) {
    Tree t;
    t.rep(1, 0, 4);
    t.rep(2, 0, 2);
    t.rep(3, 0, 4);   // the same order as r1, and after it in the array
    t.rep(4, 0, -7);
    EXPECT_STREQ("r4 r2 r1 r3", describe(t));
}

UTEST(representation_tree, ids_of_groups_and_representations_are_apart) {
    Tree t;
    t.group(1);
    t.rep(1, 1);
    EXPECT_EQ(0, t.find_group(1));
    EXPECT_EQ(1, t.find_rep(1));
    EXPECT_EQ(-1, t.find_group(2));
    EXPECT_EQ(-1, t.find_rep(0));     // 0 is the top level, not a node
    EXPECT_STREQ("g1[r1]", describe(t));
}

UTEST(representation_tree, normalize_numbers_depth_first_two_apart) {
    Tree t;
    t.rep(1, 0, 5);
    t.group(1, 0, 1);
    t.rep(2, 1, 9);
    t.rep(3, 1, 2);
    t.group(2, 1, 100);
    t.rep(4, 2, 0);

    representation_tree_normalize(t.nodes, t.count);
    EXPECT_STREQ("g1[r3 r2 g2[r4]] r1", describe(t));
    EXPECT_TRUE(orders_are_depth_first(t));
    EXPECT_EQ(0,  t.nodes[t.find_group(1)]->order);
    EXPECT_EQ(2,  t.nodes[t.find_rep(3)]->order);
    EXPECT_EQ(6,  t.nodes[t.find_group(2)]->order);
    EXPECT_EQ(8,  t.nodes[t.find_rep(4)]->order);
    EXPECT_EQ(10, t.nodes[t.find_rep(1)]->order);

    // Twice is once
    representation_tree_normalize(t.nodes, t.count);
    EXPECT_STREQ("g1[r3 r2 g2[r4]] r1", describe(t));
    EXPECT_TRUE(orders_are_depth_first(t));
}

UTEST(representation_tree, normalize_takes_missing_parents_and_cycles_to_the_top) {
    Tree t;
    t.rep(1, 7, 0);     // group 7 does not exist: a workspace written by hand, or a group removed under it
    t.group(1, 2, 1);   // g1 in g2 and g2 in g1: neither is reachable from the top level
    t.group(2, 1, 2);
    t.rep(2, 2, 3);

    representation_tree_normalize(t.nodes, t.count);
    EXPECT_EQ(0u, t.nodes[t.find_rep(1)]->parent);
    // The first group of the cycle comes out, the rest of it hangs below as it did
    EXPECT_STREQ("r1 g1[g2[r2]]", describe(t));
    EXPECT_TRUE(orders_are_depth_first(t));
}

UTEST(representation_tree, ancestors_innermost_first) {
    Tree t;
    t.group(1);
    t.group(2, 1);
    t.group(3, 2);
    const size_t r = t.rep(1, 3);
    t.rep(2);

    int out[8];
    const size_t n = representation_tree_ancestors(out, 8, t.nodes, t.count, r);
    ASSERT_EQ((size_t)3, n);
    EXPECT_EQ(t.find_group(3), out[0]);
    EXPECT_EQ(t.find_group(2), out[1]);
    EXPECT_EQ(t.find_group(1), out[2]);

    // A short buffer still says how many there are
    EXPECT_EQ((size_t)3, representation_tree_ancestors(out, 1, t.nodes, t.count, r));
    EXPECT_EQ((size_t)0, representation_tree_ancestors(out, 8, t.nodes, t.count, (size_t)t.find_rep(2)));

    EXPECT_TRUE (representation_tree_is_within(t.nodes, t.count, r, 1));
    EXPECT_TRUE (representation_tree_is_within(t.nodes, t.count, r, 3));
    EXPECT_FALSE(representation_tree_is_within(t.nodes, t.count, (size_t)t.find_group(3), 3));   // not within itself
    EXPECT_FALSE(representation_tree_is_within(t.nodes, t.count, (size_t)t.find_group(1), 2));
    EXPECT_FALSE(representation_tree_is_within(t.nodes, t.count, r, 0));
}

UTEST(representation_tree, every_placement_a_drop_can_ask_for) {
    Tree t;
    t.rep(1);
    t.group(1);
    t.rep(2, 1);
    t.rep(3, 1);
    t.rep(4);
    representation_tree_normalize(t.nodes, t.count);
    ASSERT_STREQ("r1 g1[r2 r3] r4", describe(t));

    const size_t r1 = (size_t)t.find_rep(1);
    const size_t r4 = (size_t)t.find_rep(4);
    const int g1 = t.find_group(1);
    const int r2 = t.find_rep(2);
    const int r3 = t.find_rep(3);

    EXPECT_TRUE(representation_tree_move(t.nodes, t.count, r4, (int)r1, RepresentationDrop_Before));
    EXPECT_STREQ("r4 r1 g1[r2 r3]", describe(t));

    EXPECT_TRUE(representation_tree_move(t.nodes, t.count, r4, r2, RepresentationDrop_After));
    EXPECT_STREQ("r1 g1[r2 r4 r3]", describe(t));

    EXPECT_TRUE(representation_tree_move(t.nodes, t.count, r1, g1, RepresentationDrop_IntoFirst));
    EXPECT_STREQ("g1[r1 r2 r4 r3]", describe(t));

    EXPECT_TRUE(representation_tree_move(t.nodes, t.count, r1, g1, RepresentationDrop_IntoLast));
    EXPECT_STREQ("g1[r2 r4 r3 r1]", describe(t));

    // Out to the top level, at the end
    EXPECT_TRUE(representation_tree_move(t.nodes, t.count, (size_t)r3, -1, RepresentationDrop_IntoLast));
    EXPECT_STREQ("g1[r2 r4 r1] r3", describe(t));

    // Before the group, from inside it
    EXPECT_TRUE(representation_tree_move(t.nodes, t.count, r4, g1, RepresentationDrop_Before));
    EXPECT_STREQ("r4 g1[r2 r1] r3", describe(t));

    EXPECT_TRUE(orders_are_depth_first(t));
}

UTEST(representation_tree, after_a_group_passes_over_what_it_holds) {
    Tree t;
    t.rep(1);
    t.group(1);
    t.rep(2, 1);
    t.group(2, 1);
    t.rep(3, 2);
    t.rep(4);
    representation_tree_normalize(t.nodes, t.count);
    ASSERT_STREQ("r1 g1[r2 g2[r3]] r4", describe(t));

    EXPECT_TRUE(representation_tree_move(t.nodes, t.count, (size_t)t.find_rep(1), t.find_group(1), RepresentationDrop_After));
    EXPECT_STREQ("g1[r2 g2[r3]] r1 r4", describe(t));

    // And a group moves with what it holds
    EXPECT_TRUE(representation_tree_move(t.nodes, t.count, (size_t)t.find_group(2), t.find_rep(4), RepresentationDrop_After));
    EXPECT_STREQ("g1[r2] r1 r4 g2[r3]", describe(t));
    EXPECT_TRUE(orders_are_depth_first(t));
}

UTEST(representation_tree, a_group_never_goes_within_itself) {
    Tree t;
    t.group(1);
    t.rep(1, 1);
    t.group(2, 1);
    t.rep(2, 2);
    t.rep(3);
    representation_tree_normalize(t.nodes, t.count);
    ASSERT_STREQ("g1[r1 g2[r2]] r3", describe(t));

    const size_t g1 = (size_t)t.find_group(1);
    const int g2 = t.find_group(2);
    const int r1 = t.find_rep(1);
    const int r2 = t.find_rep(2);

    struct { int target; RepresentationDrop drop; } refused[] = {
        { (int)g1, RepresentationDrop_IntoLast  },  // into itself
        { (int)g1, RepresentationDrop_Before    },  // onto itself
        { g2,      RepresentationDrop_IntoFirst },  // into a group within it
        { r2,      RepresentationDrop_Before    },  // beside something within it, so within it too
        { r1,      RepresentationDrop_After     },
    };
    for (size_t i = 0; i < ARRAY_SIZE(refused); ++i) {
        EXPECT_FALSE(representation_tree_can_move(t.nodes, t.count, g1, refused[i].target, refused[i].drop));
        EXPECT_FALSE(representation_tree_move(t.nodes, t.count, g1, refused[i].target, refused[i].drop));
        EXPECT_STREQ("g1[r1 g2[r2]] r3", describe(t));
    }

    // Nothing goes into a representation, and the top level takes only the end
    EXPECT_FALSE(representation_tree_can_move(t.nodes, t.count, (size_t)t.find_rep(3), r1, RepresentationDrop_IntoFirst));
    EXPECT_FALSE(representation_tree_can_move(t.nodes, t.count, (size_t)t.find_rep(3), -1, RepresentationDrop_Before));

    // The other way round is fine: the inner group out past the outer one
    EXPECT_TRUE(representation_tree_move(t.nodes, t.count, (size_t)g2, t.find_rep(3), RepresentationDrop_After));
    EXPECT_STREQ("g1[r1] r3 g2[r2]", describe(t));
}

UTEST(representation_tree, lifted_out_of_a_group_is_where_the_group_was) {
    // What ungrouping does: what the group holds goes to the group around it, the group goes away
    Tree t;
    t.rep(1);
    t.group(1);
    t.rep(2, 1);
    t.group(2, 1);
    t.rep(3, 2);
    t.rep(4, 1);
    t.rep(5);
    representation_tree_normalize(t.nodes, t.count);
    ASSERT_STREQ("r1 g1[r2 g2[r3] r4] r5", describe(t));

    for (size_t i = 0; i < t.count; ++i) {
        if (t.nodes[i]->parent == 1) t.nodes[i]->parent = 0;
    }
    // The group itself would be removed; here it is left empty, in its place
    representation_tree_normalize(t.nodes, t.count);
    EXPECT_STREQ("r1 g1[] r2 g2[r3] r4 r5", describe(t));
}
