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

using Label    = StringBuffer<8>;
using AtomIdx  = int32;
using ResIdx   = int16;
using ChainIdx = int16;

struct Residue {
	Label id;
	AtomIdx beg_atom_idx = 0;
	AtomIdx end_atom_idx = 0;
	ChainIdx chain_idx = -1;
};

struct Chain {
	Label id;
	ResIdx beg_res_idx;
	ResIdx end_res_idx;
};

struct BackboneSegment {
	int32 ca_idx;
	//int32 ha_idx;
	//int32 cb_idx;
	int32 n_idx;
	int32 c_idx;
	int32 o_idx;
};


// Interface to access molecular data
struct MoleculeStructure {
	Array<vec3>		atom_positions;
	Array<Element>	atom_elements;
	Array<Label>	atom_labels;
	Array<ResIdx>	atom_residue_indices;

	Array<Bond>		bonds;
	Array<Residue>	residues;
	Array<Chain>	chains;

	// If this is not zero in length it shall have the same length as residues
	Array<BackboneSegment> backbone_segments;
};

MoleculeStructure* allocate_molecule_structure(int num_atoms, int num_bonds, int num_residues, int num_chains, int num_backbone_segments);
void free_molecule_structure(MoleculeStructure* mol);