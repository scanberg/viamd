#include "gro_utils.h"
#include <mol/element.h>
#include <core/string_utils.h>

GroResult load_gro_from_file(const char* filename, Allocator& alloc) {
	return parse_gro_from_string(read_textfile(filename), alloc);
}

GroResult parse_gro_from_string(CString gro_string, Allocator& alloc) {

	CString header;
	CString length;

	extract_line(header, gro_string);
	extract_line(length, gro_string);

	int count = to_int(length);

    if (count == 0) {
        return {false};
    }

    vec3 pos, vel, box;
    char atom_name[6], res_name[6];
	atom_name[5] = '\0';
	res_name[5] = '\0';
    int atom_idx, res_idx;
    int cur_res = -1;

    struct GroAtom {
		vec3 position;
		vec3 velocity;
		Label label {};
		Element element;
		int32 residue_idx;
	};

	DynamicArray<GroAtom> atoms;
	DynamicArray<Residue> residues;
	DynamicArray<Chain> chains;
	DynamicArray<Bond> bonds;

	char buffer[256] = {};
	String line(buffer, 256);
    for (int i = 0; i < count; ++i) {
        // Get line first and then scanf the line to avoid bug when velocities are not present in data
		copy_line(line, gro_string);
        auto result =
            sscanf(line, "%5d%5c%5c%5d%8f%8f%8f%8f%8f%8f", &res_idx, res_name, atom_name,
                    &atom_idx, &pos.x, &pos.y, &pos.z, &vel.x, &vel.y, &vel.z);

        if (result > 0) {
            if (cur_res != res_idx) {
                cur_res = res_idx;
				//auto amino = aminoacid::getFromString(res_name_trim);
				//if (amino != AminoAcid::Unknown) {
				//	mol.pushStructure<structure::AminoAcid>(static_cast<int>(residues.size()), amino);
				//}
				Residue res {res_name, i, 0};
                residues.push_back(res);
            }
            residues.back().end_atom_idx++;

			auto elem = element::get_from_string(trim(CString(atom_name)));

			atoms.push_back({pos, vel, atom_name, elem, res_idx});

			//gro.atom_labels.push_back(trim(CString(atom_name)));
			//gro.atom_elements.push_back(elem);
			//gro.atom_positions.push_back(pos);
            //mol.pushAtom<field::Label, field::Element, field::Position, field::Velocity>(atom_name_trim, elem, pos, vel);
        }
    }

	//std::getline(ss, buffer); // Get simulation box
	copy_line(line, gro_string);
    sscanf(line, "%8f %8f %8f", &box.x, &box.y, &box.z);

    // Convert from nm to ångström
    for (auto& atom : atoms) {
        atom.position *= 10.f;
    }
    box *= 10.f;

	GroStructure gro;
	MoleculeStructure& mol = gro;
	mol = allocate_molecule_structure(atoms.size(), bonds.size(), residues.size(), chains.size(), MOL_ALL, alloc);

	// Copy data into molecule
	memcpy(mol.residues.data, residues.data, residues.size() * sizeof(Residue));
	memcpy(mol.chains.data, chains.data, chains.size() * sizeof(Chain));
	memcpy(mol.bonds.data, bonds.data, bonds.size() * sizeof(Bond));

	// Atom data
	for (int i = 0; i < atoms.size(); i++) {
		mol.atom_positions[i] = atoms[i].position;
		mol.atom_elements[i] = atoms[i].element;
		mol.atom_labels[i] = atoms[i].label;
	}

    //computeBonds(mol);
	//computeBackbones(mol);
    gro.box = mat3(vec3(box.x, 0, 0), vec3(0, box.y, 0), vec3(0, 0, box.z));

	return { true, gro };
}