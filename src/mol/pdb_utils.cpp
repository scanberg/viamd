#include "pdb_utils.h"
#include <mol/element.h>
#include <core/string_utils.h>

static inline bool valid_line(CString line, uint32 options) {
	if ((options & PDB_READ_ATOM) && compare_n(line, "ATOM", 4))
		return true;
	else if ((options & PDB_READ_HETATM) && compare_n(line, "HETATM", 6))
		return true;

	return false;
}

PdbResult load_pdb_from_file(const char* filename, PdbLoadParams params, Allocator& alloc) {
	return parse_pdb_from_string(read_textfile(filename), params, alloc);
}

PdbResult parse_pdb_from_string(CString pdb_string, PdbLoadParams params, Allocator& alloc) {
	CString line;

	int current_res_id = -1;
	char current_chain_id = -1;

	struct PdbAtom {
		vec3 position;
		Label label {};
		Element element;
		int32 residue_idx;
		float occupancy = 0;
		float temp_factor = 0;
	};

	DynamicArray<PdbAtom> atoms;
	DynamicArray<Residue> residues;
	DynamicArray<Chain> chains;
	DynamicArray<Bond> bonds;

	while (extract_line(line, pdb_string)) {
		// printf("%s\n", line.data);
		if (valid_line(line, params)) {
			PdbAtom atom;
			atom.position.x = to_float(line.substr(30, 8));
			atom.position.y = to_float(line.substr(38, 8));
			atom.position.z = to_float(line.substr(46, 8));
			if (line.count > 60) atom.occupancy = to_float(line.substr(54, 6));
			if (line.count > 66) atom.temp_factor = to_float(line.substr(60, 6));
			atom.label = trim(line.substr(12, 4));

			// Try to determine element from optional element column first, then name
			atom.element = element::get_from_string(line.substr(76, 2));
			if (atom.element == Element::Unknown) {
				atom.element = element::get_from_string(atom.label);
			}

			auto res_id = to_int(line.substr(22, 4));
			auto chain_id = line[21];

			// New Chain
			if (current_chain_id != chain_id) {
				current_chain_id = chain_id;
				Chain chain;
				chain.beg_res_idx = static_cast<int>(residues.count);
				chain.end_res_idx = chain.beg_res_idx;
				chain.id = chain_id;
				chains.push_back(chain);
			}

			// New Residue
			if (res_id != current_res_id) {
				current_res_id = res_id;

				Residue residue;
				residue.beg_atom_idx = static_cast<int>(atoms.count);
				residue.end_atom_idx = residue.beg_atom_idx;
				copy(String(residue.id.data, Label::MAX_LENGTH-1), trim(line.substr(17, 3)));
				residues.push_back(residue);

				// TODO: Match against Amino Acid?
			}
			residues.back().end_atom_idx++;

			// Add Atom
			atoms.push_back(atom);
		}
	}

	// TODO: Compute bonds?


	PdbStructure pdb;
	MoleculeStructure& mol = pdb;
	mol = allocate_molecule_structure(atoms.count, bonds.count, residues.count, chains.count, MOL_ALL, alloc);

	// Copy data into molecule
	memcpy(mol.residues.data, residues.data, residues.count * sizeof(Residue));
	memcpy(mol.chains.data, chains.data, chains.size() * sizeof(Chain));
	memcpy(mol.bonds.data, bonds.data, bonds.size() * sizeof(Bond));

	// Atom data
	for (int i = 0; i < atoms.size(); i++) {
		mol.atom_positions[i] = atoms[i].position;
		mol.atom_elements[i] = atoms[i].element;
		mol.atom_labels[i] = atoms[i].label;
	}

	return { true, pdb };
}