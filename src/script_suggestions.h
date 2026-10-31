#pragma once

#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_str.h>

#include <stddef.h>
#include <stdint.h>

struct md_system_t;
struct md_allocator_i;

// The lines the Script menu of the viewport's context menu offers for what is selected: the selection, to keep in a
// variable, and the distance, angle or dihedral of atoms picked one by one.
//
// A line is only offered when it does what it says. The expression of a selection is evaluated against the system
// (script::filter) and has to pick exactly the selected atoms. A measurement has to take the picked atoms and no others:
// each of its arguments is evaluated in the context the line puts it in and has to pick its atom alone, and the line
// has to compile. Whatever does not hold is dropped, so a form which does not fit the system - chain() of a molecule
// which is not a chain, resid() where the sequence numbers repeat, an atom name that is not unique - never reaches
// the menu, rather than failing once it is in the script.
//
// Nothing here knows of the application: the system and the atoms in, the lines out, which is what lets
// test/test_script_suggestions.cpp check them against real files.

struct ScriptSuggestion {
    str_t text;     // The line, ready to append to the script: 'sel1 = name("CA") and chain("A");'
    str_t hint;     // What the line gives, which it does not show itself: "11 residues", "one value for each of 12 ALA"
};

// Ways of writing the selected atoms, the most telling first: by what they are (water(), name("CA") and chain("A")),
// by the chains, residues and atom names which hold them, by their positions in the system, and as atom indices.
// ident is the name the lines assign to. Empty for an empty selection.
md_array(ScriptSuggestion) script_suggest_selection(str_t ident, md_bitset_t sel, const md_system_t& sys, md_allocator_i* alloc);

// The distance (2 atoms), angle (3) or dihedral (4) of the atoms, in the order they are given (the order they were
// picked: the vertex of an angle is the second). By index; then by the residues which hold them, as the system counts
// them and as the file numbers them; and, where they share a residue, the same measurement in every residue of its
// name. Empty for any other number of atoms.
md_array(ScriptSuggestion) script_suggest_measurement(str_t ident, const int32_t* atoms, size_t count, const md_system_t& sys, md_allocator_i* alloc);
