#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <core/md_str.h>
#include <md_fragment.h>

struct md_system_t;
struct md_allocator_i;

// The fragments of a system (sys->fragment, md_fragment.h) as a tree to show, over the levels the system
// already has:
//
//   set        molecules alike (one template), two or more: "water x24"
//   molecule   one molecule of several fragments: named by its residue, "protein", "nucleic acid", or its formula
//   residue    a residue of a molecule of several residues, with the fragments it is the first of
//   fragment   a leaf: a fragment, whole
//
// A level of one child is not made: a molecule of one fragment is that fragment, a residue which is one fragment is
// a leaf. A fragment over residues (a peptide group, a phosphate, a disulfide) is shown whole under the first of them,
// its name saying the next ("peptide, with ALA 4"): every atom is in one leaf, so the tree's quantities, which add up
// over its atoms, are right, but a residue's node is its fragments, not its atoms (of a protein: its alpha carbon,
// its side chain and the peptide group to the next residue). Names repeated under one parent are numbered
// ("phenyl 1", "phenyl 2").
//
// The nodes are in pre-order, a parent before its children, and their atoms are ordered so that every node is one
// contiguous range of atom_idx, its children's ranges inside its own.

#define FRAG_TREE_NAME_CAP 48

enum FragTreeKind : uint8_t {
    FRAG_TREE_SET = 0,
    FRAG_TREE_MOLECULE,
    FRAG_TREE_RESIDUE,
    FRAG_TREE_FRAGMENT,
};

struct FragmentTree {
    size_t    count = 0;
    uint32_t* parent = nullptr;         // [count] MD_FRAGMENT_NONE for a root
    uint8_t*  kind = nullptr;           // [count] FragTreeKind
    uint8_t*  frag_kind = nullptr;      // [count] md_fragment_kind_t of a leaf (of its fragment), MD_FRAGMENT_KIND_COUNT else
    char*     name = nullptr;           // [count * FRAG_TREE_NAME_CAP]
    uint32_t* atom_beg = nullptr;       // [count] the node's atoms are atom_idx[atom_beg .. atom_end)
    uint32_t* atom_end = nullptr;
    int32_t*  atom_idx = nullptr;
    size_t    num_atom_idx = 0;
    uint32_t* atom_fragment = nullptr;  // [num_atoms] the leaf of each atom of the system, MD_FRAGMENT_NONE for none
    size_t    num_atoms = 0;
};

// Builds the tree of the molecules which have an atom in 'include' (an atom count of flags, NULL for all) into
// 'alloc', which owns it. Needs the system's fragments (md_util_system_infer_fragments).
bool frag_tree_build(FragmentTree* out, const struct md_system_t* sys, const bool* include, struct md_allocator_i* alloc);

static inline md_str_t frag_tree_name(const FragmentTree* t, size_t idx) {
    return idx < t->count ? md_str_from_cstr(t->name + idx * FRAG_TREE_NAME_CAP) : md_str_t{};
}

static inline bool frag_tree_is_leaf(const FragmentTree* t, size_t idx) {
    // Children follow their parent directly, so a node is a leaf when the next one is not its child
    return idx < t->count && (idx + 1 >= t->count || t->parent[idx + 1] != (uint32_t)idx);
}
