#pragma once

#include <core/types.h>
#include <core/array.h>
#include <core/string_utils.h>
#include <mol/element.h>
#include <string.h>

#ifdef WIN32
#pragma warning(disable:4996) // strncpy instead of strncpy_s (windows)
#endif

struct Bond {
	int32 idx_a;
	int32 idx_b;
};

using Label = StringBuffer<8>;

struct Residue {
	Label id;
	int32 beg_atom_idx;
	int32 end_atom_idx;
	int32 chain_idx = -1;
};

struct Chain {
	Label id;
	int32 beg_res_idx;
	int32 end_res_idx;
};

// Interface to access molecular data
struct MoleculeStructure {
	Array<vec3>		atom_positions;
	Array<Element>	atom_elements;
	Array<Label>	atom_labels;
	Array<int32>	atom_residue_indices;

	Array<Bond>		bonds;
	Array<Residue>	residues;
	Array<Chain>	chains;
};

enum MoleculeStructureAllocationFlags {
	MOL_POSITIONS = BIT(0),
	MOL_ELEMENTS = BIT(1),
	MOL_LABELS = BIT(2),
	MOL_RESIDUE_INDICES = BIT(3),
	MOL_ALL = 0xffffffff
};

MoleculeStructure allocate_molecule_structure(int atom_count, int bond_count, int residue_count, int chain_count, MoleculeStructureAllocationFlags alloc_flags = MOL_ALL, Allocator* allocator = &default_alloc);

// This is a bit risky
void free_molecule_structure(MoleculeStructure& mol, Allocator* allocator = &default_alloc);