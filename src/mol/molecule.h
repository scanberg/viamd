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
struct MoleculeInterface {
	Array<vec3>		atom_positions;
	Array<Element>	atom_elements;
	Array<Label>	atom_labels;
	Array<int32>	atom_residue_indices;

	Array<Bond>		bonds;
	Array<Residue>	residues;
	Array<Chain>	chains;
};

struct MoleculeData {
	DynamicArray<vec3>		atom_position_data;
	DynamicArray<Element>	atom_element_data;
	DynamicArray<Label>		atom_label_data;
	DynamicArray<int32>		atom_residue_index_data;

	DynamicArray<Bond>		bond_data;
	DynamicArray<Residue>	residue_data;
	DynamicArray<Chain>		chain_data;

	operator MoleculeInterface() {
		return { atom_position_data, atom_element_data, atom_label_data, atom_residue_index_data, bond_data, residue_data, chain_data };
	}
};

enum MoleculeStructureAllocationFlags {
	MOL_POSITIONS = BIT(0),
	MOL_ELEMENTS = BIT(1),
	MOL_LABELS = BIT(2),
	MOL_RESIDUE_INDICES = BIT(3),
	MOL_ALL = 0xffffffff
};

MoleculeInterface allocate_molecule_structure(int atom_count, int bond_count, int residue_count, int chain_count, MoleculeStructureAllocationFlags alloc_flags = MOL_ALL, Allocator* allocator = &default_alloc);

// This is a bit risky
void free_molecule_structure(MoleculeInterface& mol, Allocator* allocator = &default_alloc);