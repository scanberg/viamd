#include "pdb_utils.h"
#include <mol/element.h>
#include <mol/molecule_utils.h>
#include <mol/trajectory_utils.h>
#include <core/string_utils.h>

static inline bool valid_line(CString line, uint32 options) {
    if ((options & PDB_READ_ATOM) && compare_n(line, "ATOM", 4))
        return true;
    else if ((options & PDB_READ_HETATM) && compare_n(line, "HETATM", 6))
        return true;

    return false;
}

bool allocate_and_load_pdb_from_file(MoleculeDynamic* md, const char* filename, PdbLoadParams params) {
    String txt = allocate_and_read_textfile(filename);
    auto res = allocate_and_parse_pdb_from_string(md, txt, params);
    FREE(txt);
    return res;
}

bool allocate_and_parse_pdb_from_string(MoleculeDynamic* md, CString pdb_string, PdbLoadParams params) {
    free_molecule_structure(&md->molecule);
    free_trajectory(&md->trajectory);

    DynamicArray<vec3> positions;
    DynamicArray<Label> labels;
    DynamicArray<Element> elements;
    DynamicArray<ResIdx> residue_indices;
    DynamicArray<float> occupancies;
    DynamicArray<float> temp_factors;
    DynamicArray<Residue> residues;
    DynamicArray<Chain> chains;

    int current_res_id = -1;
    char current_chain_id = -1;
    int num_atoms = 0;
    int num_models = 0;
    int num_frames = 0;
    mat3 box(0);
    CString line;
    while (extract_line(line, pdb_string)) {
        if (compare_n(line, "CRYST1", 6)) {
            vec3 dim(to_float(line.substr(6, 9)), to_float(line.substr(15, 9)), to_float(line.substr(24, 9)));
            vec3 angles(to_float(line.substr(33, 7)), to_float(line.substr(40, 7)), to_float(line.substr(47, 7)));
            box[0].x = dim.x;
            box[1].y = dim.y;
            box[2].z = dim.z;
        } else if (compare_n(line, "MODEL", 5)) {

        } else if (compare_n(line, "ENDMDL", 6)) {
            if (params & PDB_TREAT_MODELS_AS_FRAMES) {
                /*
                if (!md->molecule) {
                        bonds = compute_covalent_bonds(positions, elements, residues);
                        backbone_segments = compute_backbone_segments(residues, labels);

                        init_molecule_structure(&md->molecule, num_atoms, (int32)bonds.count, (int32)residues.count, (int32)chains.count,
                (int32)backbone_segments.count);

                        // Copy data into molecule
                        memcpy(md->molecule.atom_positions.data, positions.data, positions.size() * sizeof(vec3));
                        memcpy(md->molecule.atom_elements.data, elements.data, elements.size() * sizeof(Element));
                        memcpy(md->molecule.atom_labels.data, labels.data, labels.size() * sizeof(Label));
                        memcpy(md->molecule.atom_residue_indices.data, residue_indices.data, residue_indices.size() * sizeof(ResIdx));

                        memcpy(md->molecule.residues.data, residues.data, residues.size() * sizeof(Residue));
                        memcpy(md->molecule.chains.data, chains.data, chains.size() * sizeof(Chain));
                        memcpy(md->molecule.bonds.data, bonds.data, bonds.size() * sizeof(Bond));
                        memcpy(md->molecule.backbone_segments.data, backbone_segments.data, backbone_segments.size() * sizeof(BackboneSegment));
                }
                */

                num_frames++;
            } else {
                ASSERT(false);
                num_models++;
                num_atoms = 0;
                current_res_id = -1;
                current_chain_id = -1;

                positions.clear();
                labels.clear();
                elements.clear();
                residue_indices.clear();
                occupancies.clear();
                temp_factors.clear();
                residues.clear();
                chains.clear();
            }
        } else if (compare_n(line, "TER", 3)) {
            if (params & PDB_TREAT_MODELS_AS_FRAMES && num_frames > 0) continue;
            current_chain_id = line[21];
            Chain chain;
            chain.beg_res_idx = (ResIdx)residues.size();
            chain.end_res_idx = chain.beg_res_idx;
            chain.id = current_chain_id;
            chains.push_back(chain);
        } else if (valid_line(line, params)) {
            positions.push_back(vec3(to_float(line.substr(30, 8)), to_float(line.substr(38, 8)), to_float(line.substr(46, 8))));
            if (params & PDB_TREAT_MODELS_AS_FRAMES && num_frames > 0) continue;
            labels.push_back(trim(line.substr(12, 4)));
            if (line.count > 60) {
                occupancies.push_back(to_float(line.substr(54, 6)));
            }
            if (line.count > 66) {
                temp_factors.push_back(to_float(line.substr(60, 6)));
            }

            // Try to determine element from optional element column first, then name

            Element elem = Element::Unknown;
            if (line.count >= 78) {
                element::get_from_string(line.substr(76, 2));
            }
            if (elem == Element::Unknown) elem = element::get_from_string(labels.back());
            elements.push_back(elem);

            auto res_id = to_int(line.substr(22, 4));
            char chain_id = line[21];

            // New Chain
            if (current_chain_id != chain_id) {
                current_chain_id = chain_id;
                Chain chain;
                chain.beg_res_idx = (ResIdx)residues.size();
                chain.end_res_idx = chain.beg_res_idx;
                chain.id = chain_id;
                chains.push_back(chain);
            }

            // New Residue
            if (res_id != current_res_id) {
                current_res_id = res_id;
                Residue residue;
                residue.id = res_id;
                residue.beg_atom_idx = num_atoms;
                residue.end_atom_idx = residue.beg_atom_idx;
                residue.chain_idx = (ChainIdx)(chains.size() - 1);
                copy(String(residue.name.beg(), Label::MAX_LENGTH - 1), trim(line.substr(17, 3)));
                residues.push_back(residue);
                chains.back().end_res_idx++;
            }
            residues.back().end_atom_idx++;

            residue_indices.push_back((ResIdx)(residues.size() - 1));

            // Add Atom
            num_atoms++;
        }
    }

    if (!md->molecule) {
        auto mol_pos = positions.sub_array(0, num_atoms);
        auto covalent_bonds = compute_covalent_bonds(mol_pos, elements, residues);
        auto backbone_segments = compute_backbone_segments(residues, labels);
        auto donors = hydrogen_bond::compute_donors(labels);
        auto acceptors = hydrogen_bond::compute_acceptors(elements);

        init_molecule_structure(&md->molecule, num_atoms, (int32)covalent_bonds.count, (int32)residues.count, (int32)chains.count,
                                (int32)backbone_segments.count, (int32)donors.count, (int32)acceptors.count);

        // Copy data into molecule
        memcpy(md->molecule.atom_positions.data, mol_pos.data, mol_pos.size_in_bytes());
        memcpy(md->molecule.atom_elements.data, elements.data, elements.size_in_bytes());
        memcpy(md->molecule.atom_labels.data, labels.data, labels.size_in_bytes());
        memcpy(md->molecule.atom_residue_indices.data, residue_indices.data, residue_indices.size_in_bytes());

        memcpy(md->molecule.residues.data, residues.data, residues.size_in_bytes());
        memcpy(md->molecule.chains.data, chains.data, chains.size_in_bytes());
        memcpy(md->molecule.covalent_bonds.data, covalent_bonds.data, covalent_bonds.size_in_bytes());
        memcpy(md->molecule.backbone_segments.data, backbone_segments.data, backbone_segments.size_in_bytes());
        memcpy(md->molecule.hydrogen_bond.donors.data, donors.data, donors.size_in_bytes());
        memcpy(md->molecule.hydrogen_bond.acceptors.data, acceptors.data, acceptors.size_in_bytes());
    }

    if (num_frames > 0) {
        init_trajectory(&md->trajectory, num_atoms, num_frames);
        // COPY POSITION DATA

        ASSERT(positions.count > 0);
        memcpy(md->trajectory.position_data.data, positions.data, sizeof(vec3) * positions.count);

        for (int i = 0; i < md->trajectory.num_frames; i++) {
            int index = i;
            float time = 0;
            Array<vec3> atom_positions{md->trajectory.position_data.beg() + i * num_atoms, num_atoms};
            md->trajectory.frame_buffer[i] = {index, time, box, atom_positions};
        }
    }

    return true;
}
