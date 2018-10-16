#pragma once

#include <core/types.h>
#include <core/array_types.h>
#include <core/vector_types.h>
#include <core/string_utils.h>
#include <mol/element.h>

using Label = StringBuffer<8>;
using AtomIdx = int32;
using ResIdx = int32;
using ChainIdx = int32;

struct Bond {
    int32 idx_a;
    int32 idx_b;
};

struct BackboneSegment {
    AtomIdx ca_idx;
    AtomIdx n_idx;
    AtomIdx c_idx;
    AtomIdx o_idx;
};

struct HydrogenBondDonor {
    AtomIdx donor_idx = 0;
    AtomIdx hydro_idx = 0;
};

typedef AtomIdx HydrogenBondAcceptor;

struct Residue {
    Label name = "";
    ResIdx id = -1;
    AtomIdx beg_atom_idx = 0;
    AtomIdx end_atom_idx = 0;
    ChainIdx chain_idx = -1;
};

struct Chain {
    Label id = "";
    ResIdx beg_res_idx = 0;
    ResIdx end_res_idx = 0;
};

// Interface to access molecular data
struct MoleculeStructure {
    Array<vec3> atom_positions{};
    Array<Element> atom_elements{};
    Array<Label> atom_labels{};
    Array<ResIdx> atom_residue_indices{};

    Array<Bond> covalent_bonds{};
    Array<Residue> residues{};
    Array<Chain> chains{};

    // If this is not zero in length it should have the same length as residues
    Array<BackboneSegment> backbone_segments{};

    struct {
        Array<HydrogenBondDonor> donors{};
        Array<HydrogenBondAcceptor> acceptors{};
    } hydrogen_bond;

    operator bool() const { return atom_positions.count > 0; }
};

// Chain func

inline AtomIdx get_atom_beg_idx(MoleculeStructure& mol, Chain chain) {
    ASSERT(0 <= chain.beg_res_idx && chain.beg_res_idx < mol.residues.count);
    return mol.residues[chain.beg_res_idx].beg_atom_idx;
}

inline AtomIdx get_atom_beg_idx(const MoleculeStructure& mol, Chain chain) {
    ASSERT(0 <= chain.beg_res_idx && chain.beg_res_idx < mol.residues.count);
    return mol.residues[chain.beg_res_idx].beg_atom_idx;
}

inline AtomIdx get_atom_end_idx(MoleculeStructure& mol, Chain chain) {
    ASSERT(0 < chain.end_res_idx && chain.end_res_idx <= mol.residues.count);
    return mol.residues[chain.end_res_idx - 1].end_atom_idx;
}

inline AtomIdx get_atom_end_idx(const MoleculeStructure& mol, Chain chain) {
    ASSERT(0 < chain.end_res_idx && chain.end_res_idx <= mol.residues.count);
    return mol.residues[chain.end_res_idx - 1].end_atom_idx;
}

inline Chain get_chain(MoleculeStructure& mol, ChainIdx idx) {
    ASSERT(0 <= idx && idx < mol.chains.count);
    return mol.chains[idx];
}

inline Chain get_chain(const MoleculeStructure& mol, ChainIdx idx) {
    ASSERT(0 <= idx && idx < mol.chains.count);
    return mol.chains[idx];
}

inline Array<BackboneSegment> get_backbone(MoleculeStructure& mol, Chain chain) {
    ASSERT(0 < chain.end_res_idx && chain.end_res_idx <= mol.residues.count);
    if (mol.backbone_segments.count == 0) return {};
    return {mol.backbone_segments.beg() + chain.beg_res_idx, mol.backbone_segments.beg() + chain.end_res_idx};
}

inline Array<const BackboneSegment> get_backbone(const MoleculeStructure& mol, Chain chain) {
    ASSERT(0 < chain.end_res_idx && chain.end_res_idx <= mol.residues.count);
    if (mol.backbone_segments.count == 0) return {};
    return {mol.backbone_segments.beg() + chain.beg_res_idx, mol.backbone_segments.beg() + chain.end_res_idx};
}

inline Array<Residue> get_residues(MoleculeStructure& mol, Chain chain) {
    return mol.residues.sub_array(chain.beg_res_idx, chain.end_res_idx - chain.beg_res_idx);
}

inline Array<const Residue> get_residues(const MoleculeStructure& mol, Chain chain) {
    return mol.residues.sub_array(chain.beg_res_idx, chain.end_res_idx - chain.beg_res_idx);
}

inline Array<vec3> get_positions(MoleculeStructure& mol, Chain chain) {
    auto beg_atom_idx = get_atom_beg_idx(mol, chain);
    auto end_atom_idx = get_atom_end_idx(mol, chain);
    return mol.atom_positions.sub_array(beg_atom_idx, end_atom_idx - beg_atom_idx);
}

inline Array<const vec3> get_positions(const MoleculeStructure& mol, Chain chain) {
    auto beg_atom_idx = get_atom_beg_idx(mol, chain);
    auto end_atom_idx = get_atom_end_idx(mol, chain);
    return mol.atom_positions.sub_array(beg_atom_idx, end_atom_idx - beg_atom_idx);
}

inline Array<Element> get_elements(MoleculeStructure& mol, Chain chain) {
    auto beg_atom_idx = mol.residues[chain.beg_res_idx].beg_atom_idx;
    auto end_atom_idx = mol.residues[chain.end_res_idx - 1].end_atom_idx;
    return mol.atom_elements.sub_array(beg_atom_idx, end_atom_idx - beg_atom_idx);
}

inline Array<const Element> get_elements(const MoleculeStructure& mol, Chain chain) {
    auto beg_atom_idx = mol.residues[chain.beg_res_idx].beg_atom_idx;
    auto end_atom_idx = mol.residues[chain.end_res_idx - 1].end_atom_idx;
    return mol.atom_elements.sub_array(beg_atom_idx, end_atom_idx - beg_atom_idx);
}

inline Array<Label> get_labels(MoleculeStructure& mol, Chain chain) {
    auto beg_atom_idx = mol.residues[chain.beg_res_idx].beg_atom_idx;
    auto end_atom_idx = mol.residues[chain.end_res_idx - 1].end_atom_idx;
    return mol.atom_labels.sub_array(beg_atom_idx, end_atom_idx - beg_atom_idx);
}

inline Array<const Label> get_labels(const MoleculeStructure& mol, Chain chain) {
    auto beg_atom_idx = mol.residues[chain.beg_res_idx].beg_atom_idx;
    auto end_atom_idx = mol.residues[chain.end_res_idx - 1].end_atom_idx;
    return mol.atom_labels.sub_array(beg_atom_idx, end_atom_idx - beg_atom_idx);
}

// Res func
inline Array<vec3> get_positions(MoleculeStructure& mol, Residue res) {
    return mol.atom_positions.sub_array(res.beg_atom_idx, res.end_atom_idx - res.beg_atom_idx);
}

inline Array<const vec3> get_positions(const MoleculeStructure& mol, Residue res) {
    return mol.atom_positions.sub_array(res.beg_atom_idx, res.end_atom_idx - res.beg_atom_idx);
}

inline Array<Element> get_elements(MoleculeStructure& mol, Residue res) {
    return mol.atom_elements.sub_array(res.beg_atom_idx, res.end_atom_idx - res.beg_atom_idx);
}

inline Array<const Element> get_elements(const MoleculeStructure& mol, Residue res) {
    return mol.atom_elements.sub_array(res.beg_atom_idx, res.end_atom_idx - res.beg_atom_idx);
}

inline Array<Label> get_labels(MoleculeStructure& mol, Residue res) {
    return mol.atom_labels.sub_array(res.beg_atom_idx, res.end_atom_idx - res.beg_atom_idx);
}

inline Array<const Label> get_labels(const MoleculeStructure& mol, Residue res) {
    return mol.atom_labels.sub_array(res.beg_atom_idx, res.end_atom_idx - res.beg_atom_idx);
}

bool init_molecule_structure(MoleculeStructure* mol, int32 num_atoms, int32 num_bonds, int32 num_residues, int32 num_chains,
                             int32 num_backbone_segments = 0, int32 num_hydrogen_bond_donors = 0, int32 num_hydrogen_bond_acceptors = 0);
void free_molecule_structure(MoleculeStructure* mol);
