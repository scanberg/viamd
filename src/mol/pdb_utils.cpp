#include "pdb_utils.h"
#include <mol/element.h>
#include <mol/molecule_utils.h>
#include <core/string_utils.h>

static inline bool valid_line(CString line, uint32 options) {
	if ((options & PDB_READ_ATOM) && compare_n(line, "ATOM", 4))
		return true;
	else if ((options & PDB_READ_HETATM) && compare_n(line, "HETATM", 6))
		return true;

	return false;
}

PdbResult load_pdb_from_file(const char* filename, PdbLoadParams params, Allocator* alloc) {
	return parse_pdb_from_string(read_textfile(filename), params, alloc);
}

PdbResult parse_pdb_from_string(CString pdb_string, PdbLoadParams params, Allocator* alloc) {
	DynamicArray<vec3> positions;
	DynamicArray<Label> labels;
	DynamicArray<Element> elements;
	DynamicArray<int32> residue_indices;
	DynamicArray<float> occupancies;
	DynamicArray<float> temp_factors;
	DynamicArray<Residue> residues;
	DynamicArray<Chain> chains;
	DynamicArray<Bond> bonds;

	int current_res_id = -1;
	char current_chain_id = -1;
	int num_atoms = 0;
	CString line;
	while (extract_line(line, pdb_string)) {
		if (valid_line(line, params)) {
			labels.push_back(trim(line.substr(12, 4)));

			positions.push_back(vec3(to_float(line.substr(30, 8)), to_float(line.substr(38, 8)), to_float(line.substr(46, 8))));

			if (line.count > 60) {
				occupancies.push_back(to_float(line.substr(54, 6)));
			}

			if (line.count > 66) {
				temp_factors.push_back(to_float(line.substr(60, 6)));
			}

			// Try to determine element from optional element column first, then name
			Element elem = element::get_from_string(line.substr(76, 2));
			if (elem == Element::Unknown) elem = element::get_from_string(labels.back());
			elements.push_back(elem);

			residue_indices.push_back((int)residues.size());

			auto res_id = to_int(line.substr(22, 4));
			char chain_id = line[21];

			// New Chain
			if (current_chain_id != chain_id) {
				current_chain_id = chain_id;
				Chain chain;
				chain.beg_res_idx = (int)residues.count;
				chain.end_res_idx = chain.beg_res_idx;
				chain.id = chain_id;
				chains.push_back(chain);
			}

			// New Residue
			if (res_id != current_res_id) {
				current_res_id = res_id;
				Residue residue;
				residue.beg_atom_idx = num_atoms;
				residue.end_atom_idx = residue.beg_atom_idx;
				copy(String(residue.id.beg(), Label::MAX_LENGTH-1), trim(line.substr(17, 3)));
				residues.push_back(residue);
			}
			residues.back().end_atom_idx++;

			// Add Atom
			num_atoms++;
		}
	}

	bonds = compute_covalent_bonds(positions, elements, residues);

	PdbStructure pdb;
	MoleculeStructure& mol = pdb;
	mol = allocate_molecule_structure(num_atoms, bonds.size(), residues.size(), chains.size(), MOL_ALL, alloc);

	// Copy data into molecule
	memcpy(mol.atom_positions.data, positions.data, positions.size() * sizeof(vec3));
	memcpy(mol.atom_elements.data, elements.data, elements.size() * sizeof(Element));
	memcpy(mol.atom_labels.data, labels.data, labels.size() * sizeof(Label));
	memcpy(mol.atom_residue_indices.data, residue_indices.data, residue_indices.size() * sizeof(int32));

	memcpy(mol.residues.data, residues.data, residues.size() * sizeof(Residue));
	memcpy(mol.chains.data, chains.data, chains.size() * sizeof(Chain));
	memcpy(mol.bonds.data, bonds.data, bonds.size() * sizeof(Bond));

	return { true, pdb };
}