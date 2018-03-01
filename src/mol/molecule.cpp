#include "molecule.h"

MoleculeStructure* allocate_molecule_structure(int num_atoms, int num_bonds, int num_residues, int num_chains, int num_backbone_segments) {
	MoleculeStructure* mol = (MoleculeStructure*)MALLOC(sizeof(MoleculeStructure));
	int64 alloc_size = 0;

	alloc_size += num_atoms * (sizeof(vec3) + sizeof(Element) + sizeof(Label) + sizeof(ResIdx));
	alloc_size += num_bonds * sizeof(Bond);
	alloc_size += num_residues * sizeof(Residue);
	alloc_size += num_chains * sizeof(Chain);
	alloc_size += num_backbone_segments * sizeof(BackboneSegment);

	void* data = MALLOC(alloc_size);

	mol->atom_positions =		{ (vec3*)data, num_atoms };
	mol->atom_elements =		{ (Element*)(mol->atom_positions.end()), num_atoms };
	mol->atom_labels =			{ (Label*)(mol->atom_elements.end()), num_atoms };
	mol->atom_residue_indices =	{ (ResIdx*)(mol->atom_labels.end()), num_atoms };
	
	mol->bonds =			 { (Bond*)(mol->atom_residue_indices.end()), num_bonds };
	mol->residues =			 { (Residue*)(mol->bonds.end()), num_residues };
	mol->chains =			 { (Chain*)(mol->residues.end()), num_chains };
	mol->backbone_segments = { (BackboneSegment*)(mol->chains.end()), num_backbone_segments };

	return mol;
}

void free_molecule_structure(MoleculeStructure* mol) {
	ASSERT(mol);
	FREE(mol->atom_positions.beg());
	FREE(mol);
}