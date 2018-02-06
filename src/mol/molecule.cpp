#include "molecule.h"

MoleculeStructure allocate_molecule_structure(int atom_count, int bond_count, int residue_count, int chain_count, MoleculeStructureAllocationFlags alloc_flags, Allocator& allocator)
{
	MoleculeStructure mol;
	int64 alloc_size = 0;

	if (atom_count > 0) {
		if (alloc_flags & MOL_POSITIONS) alloc_size += atom_count * sizeof(vec3);
		if (alloc_flags & MOL_ELEMENTS) alloc_size += atom_count * sizeof(Element);
		if (alloc_flags & MOL_LABELS) alloc_size += atom_count * sizeof(Label);
		if (alloc_flags & MOL_RESIDUE_INDICES) alloc_size += atom_count * sizeof(int32);
	}
	if (bond_count > 0) {
		alloc_size += bond_count * sizeof(Bond);
	}
	if (residue_count > 0) {
		alloc_size += residue_count * sizeof(Residue);
	}
	if (chain_count > 0) {
		alloc_size += chain_count * sizeof(Chain);
	}

	void* data = allocator.alloc(alloc_size);

	mol.atom_positions =		{ (vec3*)data, alloc_flags & MOL_POSITIONS ? atom_count : 0 };
	mol.atom_elements =			{ (Element*)(mol.atom_positions.end()), alloc_flags & MOL_ELEMENTS ? atom_count : 0 };
	mol.atom_labels =			{ (Label*)(mol.atom_elements.end()), alloc_flags & MOL_LABELS ? atom_count : 0 };
	mol.atom_residue_indices =	{ (int32*)(mol.atom_labels.end()), alloc_flags & MOL_RESIDUE_INDICES ? atom_count : 0 };
	
	mol.bonds =		{ (Bond*)(mol.atom_residue_indices.end()), bond_count };
	mol.residues =	{ (Residue*)(mol.bonds.end()), residue_count };
	mol.chains =	{ (Chain*)(mol.residues.end()), chain_count };

	return mol;
}

void free_molecule_structure(MoleculeStructure& mol, Allocator& alloc) {
	alloc.free(mol.atom_positions.data);
}