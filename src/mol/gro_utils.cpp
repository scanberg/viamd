#include "gro_utils.h"
#include <core/string_utils.h>
#include <mol/element.h>
#include <mol/molecule_utils.h>

MoleculeStructure* allocate_and_load_gro_from_file(const char* filename) {
	String txt = allocate_and_read_textfile(filename);
	auto mol = allocate_and_parse_gro_from_string(txt);
	FREE(txt);
	return mol;
}

MoleculeStructure* allocate_and_parse_gro_from_string(CString gro_string) {

	CString header;
	CString length;

	extract_line(header, gro_string);
	extract_line(length, gro_string);

	int num_atoms = to_int(length);

    if (num_atoms == 0) {
        return nullptr;
    }

	int res_count = 0;
    int cur_res = -1;

	/*
    struct GroAtom {
		vec3 position;
		vec3 velocity;
		Label label {};
		Element element;
		int32 residue_idx;
	};
	*/

	DynamicArray<vec3> positions;
	DynamicArray<vec3> velocities;
	DynamicArray<Label> labels;
	DynamicArray<Element> elements;
	DynamicArray<ResIdx> residue_indices;
	//DynamicArray<GroAtom> atoms;
	DynamicArray<Residue> residues;

	char buffer[256] = {};
	String line(buffer, 256);
    for (int i = 0; i < num_atoms; ++i) {
		vec3 pos, vel;
		int atom_idx, res_idx;
		char atom_name[8] = {};
		char res_name[8] = {};

        // Get line first and then scanf the line to avoid bug when velocities are not present in data
		copy_line(line, gro_string);
        auto result =
            sscanf(line, "%5d%5c%5c%5d%8f%8f%8f%8f%8f%8f", &res_idx, res_name, atom_name,
                    &atom_idx, &pos.x, &pos.y, &pos.z, &vel.x, &vel.y, &vel.z);
        if (result > 0) {
			if (cur_res != res_idx) {
				cur_res = res_idx;
				res_count = (int)residues.count;
				//auto amino = aminoacid::getFromString(res_name_trim);
				//if (amino != AminoAcid::Unknown) {
				//	mol.pushStructure<structure::AminoAcid>(static_cast<int>(residues.size()), amino);
				//}
				CString res_name_trim = trim(CString(res_name));
				Residue res{ res_name_trim, i, i };
				residues.push_back(res);
			}
			residues.back().end_atom_idx++;

			CString atom_name_trim = trim(CString(atom_name));

            auto elem = element::get_from_string(atom_name_trim);
			positions.push_back(pos);
			velocities.push_back(vel);
			labels.push_back(atom_name_trim);
			elements.push_back(elem);
			residue_indices.push_back((ResIdx)res_count);



			//atoms.push_back({pos, vel, atom_name, elem, res_idx});

			//gro.atom_labels.push_back(trim(CString(atom_name)));
			//gro.atom_elements.push_back(elem);
			//gro.atom_positions.push_back(pos);
            //mol.pushAtom<field::Label, field::Element, field::Position, field::Velocity>(atom_name_trim, elem, pos, vel);
        }
    }

	//std::getline(ss, buffer); // Get simulation box
	vec3 box;
	copy_line(line, gro_string);
    sscanf(line, "%8f %8f %8f", &box.x, &box.y, &box.z);

    // Convert from nm to ångström
    for (auto& p : positions) {
        p *= 10.f;
    }
	for (auto& v : velocities) {
		v *= 10.f;
	}
    box *= 10.f;

	DynamicArray<BackboneSegment> backbone_segments = compute_backbone_segments(residues, labels);
	DynamicArray<Bond> bonds = compute_covalent_bonds(positions, elements, residues);
	DynamicArray<Chain> chains = compute_chains(residues, bonds, residue_indices);

   	for (int c = 0; c < chains.count; c++) {
	   	for (int i = chains[c].beg_res_idx; i < chains[c].end_res_idx; i++) {
	   		residues[i].chain_idx = (ChainIdx)c;
	   	}
	}

	MoleculeStructure* mol = allocate_molecule_structure(num_atoms, (int)bonds.size(), (int)residues.size(), (int)chains.size(), (int)backbone_segments.size());

	// Copy data into molecule
	memcpy(mol->atom_positions.data, positions.data, positions.size_in_bytes());
	memcpy(mol->atom_elements.data, elements.data, elements.size_in_bytes());
	memcpy(mol->atom_labels.data, labels.data, labels.size_in_bytes());
	memcpy(mol->atom_residue_indices.data, residue_indices.data, residue_indices.size_in_bytes());

	memcpy(mol->residues.data, residues.data, residues.size_in_bytes());
	memcpy(mol->chains.data, chains.data, chains.size_in_bytes());
	memcpy(mol->bonds.data, bonds.data, bonds.size_in_bytes());
	memcpy(mol->backbone_segments.data, backbone_segments.data, backbone_segments.size_in_bytes());

    //gro.box = mat3(vec3(box.x, 0, 0), vec3(0, box.y, 0), vec3(0, 0, box.z));

	return mol;
}