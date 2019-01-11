#pragma once

#include <core/math_utils.h>
#include "image.h"

enum class ColorMapping { STATIC_COLOR, CPK, RES_ID, RES_INDEX, CHAIN_ID, CHAIN_INDEX, SECONDARY_STRUCTURE };

inline vec4 color_from_hash(uint32 hash) {
    constexpr float CHROMA = 0.45f;
    constexpr float LUMINANCE = 0.90f;
    constexpr uint32 MOD = 21;
    constexpr float SCL = 1.f / (float)MOD;

    return vec4(math::hcl_to_rgb(vec3((hash % MOD) * SCL, CHROMA, LUMINANCE)), 1);
}

void color_atoms_cpk(Array<uint32> dst_atom_colors, Array<const Element> elements);
void color_atoms_residue_id(Array<uint32> dst_atom_colors, Array<const Residue> residues);
void color_atoms_residue_index(Array<uint32> dst_atom_colors, Array<const Residue> residues);
void color_atoms_chain_id(Array<uint32> dst_atom_colors, Array<const Chain> chains, Array<const Residue> residues);
void color_atoms_chain_index(Array<uint32> dst_atom_colors, Array<const Chain> chains, Array<const Residue> residues);

void color_atoms_backbone_angles(Array<uint32> dst_atom_colors, Array<const Residue> residues, Array<const BackboneSequence> bb_seq, Array<const vec2> bb_angles, const Image& ramachandran_color_map);
