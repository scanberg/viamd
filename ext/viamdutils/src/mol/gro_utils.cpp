#include "gro_utils.h"
#include <core/string_utils.h>
#include <mol/element.h>
#include <mol/molecule_utils.h>
#include <mol/hydrogen_bond.h>

bool allocate_and_load_gro_from_file(MoleculeStructure* mol, const char* filename) {
    String txt = allocate_and_read_textfile(filename);
    auto res = allocate_and_parse_gro_from_string(mol, txt);
    FREE(txt);
    return res;
}

bool allocate_and_parse_gro_from_string(MoleculeStructure* mol, CString gro_string) {
    CString header;
    CString length;

    extract_line(header, gro_string);
    extract_line(length, gro_string);

    int num_atoms = to_int(length);

    if (num_atoms == 0) {
        return false;
    }

    int res_count = 0;
    int cur_res = -1;

    DynamicArray<vec3> positions;
    DynamicArray<vec3> velocities;
    DynamicArray<Label> labels;
    DynamicArray<Element> elements;
    DynamicArray<ResIdx> residue_indices;
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
            sscanf(line, "%5d%5c%5c%5d%8f%8f%8f%8f%8f%8f", &res_idx, res_name, atom_name, &atom_idx, &pos.x, &pos.y, &pos.z, &vel.x, &vel.y, &vel.z);
        if (result > 0) {
            if (cur_res != res_idx) {
                cur_res = res_idx;
                res_count = (int)residues.count;
                CString res_name_trim = trim(CString(res_name));
                Residue res{};
                res.name = res_name_trim;
                res.id = res_idx;
                res.chain_idx = 0;
                res.atom_idx = {i, i};
                residues.push_back(res);
            }
            residues.back().atom_idx.end++;

            CString atom_name_trim = trim(CString(atom_name));
            CString element_str = atom_name_trim;

            if (is_amino_acid(residues.back())) {
                // If we have an amino acid, we can assume its an organic element with just one letter. C/N/H/O?
                element_str = element_str.substr(0, 1);
            }
            Element elem = element::get_from_string(element_str);

            positions.push_back(pos);
            velocities.push_back(vel);
            labels.push_back(atom_name_trim);
            elements.push_back(elem);
            residue_indices.push_back((ResIdx)res_count);
        }
    }

    vec3 box{};
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
    DynamicArray<Bond> covalent_bonds = compute_covalent_bonds(residues, residue_indices, positions, elements);
    DynamicArray<Chain> chains = compute_chains(residues, covalent_bonds, residue_indices);
    DynamicArray<HydrogenBondDonor> donors = hydrogen_bond::compute_donors(labels);
    DynamicArray<HydrogenBondAcceptor> acceptors = hydrogen_bond::compute_acceptors(elements);

    for (ChainIdx c = 0; c < chains.count; c++) {
        for (auto i = chains[c].res_idx.beg; i < chains[c].res_idx.end; i++) {
            residues[i].chain_idx = c;
        }
    }

    init_molecule_structure(mol, num_atoms, (int32)covalent_bonds.size(), (int32)residues.size(), (int32)chains.size(),
                            (int32)backbone_segments.size(), (int32)donors.size(), (int32)acceptors.size());

    // Copy data into molecule
    memcpy(mol->atom.positions, positions.data, positions.size_in_bytes());
    memcpy(mol->atom.elements, elements.data, elements.size_in_bytes());
    memcpy(mol->atom.labels, labels.data, labels.size_in_bytes());
    memcpy(mol->atom.residue_indices, residue_indices.data, residue_indices.size_in_bytes());

    memcpy(mol->residues.data, residues.data, residues.size_in_bytes());
    memcpy(mol->chains.data, chains.data, chains.size_in_bytes());
    memcpy(mol->covalent_bonds.data, covalent_bonds.data, covalent_bonds.size_in_bytes());
    memcpy(mol->backbone_segments.data, backbone_segments.data, backbone_segments.size_in_bytes());
    memcpy(mol->hydrogen_bond.donors.data, donors.data, donors.size_in_bytes());
    memcpy(mol->hydrogen_bond.acceptors.data, acceptors.data, acceptors.size_in_bytes());

    return true;
}
