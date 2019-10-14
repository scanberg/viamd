#pragma once

#include <core/types.h>
#include <core/array_types.h>
#include <core/bitfield.h>
#include <core/math_utils.h>
#include "image.h"

enum class ColorMapping { Uniform, Cpk, ResId, ResIndex, ChainId, ChainIndex, SecondaryStructure };

inline vec4 color_from_hash(uint32 hash) {
    constexpr float chroma = 0.8f;
    constexpr float luminance = 0.90f;
    constexpr uint32 mod = 21;
    const float hue = (hash % mod) / (float)mod;
    const vec3 rgb = math::hcl_to_rgb(vec3(hue, chroma, luminance));

    return vec4(rgb, 1);
}

void color_atoms_uniform(ArrayView<uint32> dst_atom_colors, const vec4& color);
void color_atoms_cpk(ArrayView<uint32> dst_atom_colors, ArrayView<const Element> elements);
void color_atoms_residue_id(ArrayView<uint32> dst_atom_colors, ArrayView<const Residue> residues);
void color_atoms_residue_index(ArrayView<uint32> dst_atom_colors, ArrayView<const Residue> residues);
void color_atoms_chain_id(ArrayView<uint32> dst_atom_colors, ArrayView<const Chain> chains);
void color_atoms_chain_index(ArrayView<uint32> dst_atom_colors, ArrayView<const Chain> chains);

void color_atoms_backbone_angles(ArrayView<uint32> dst_atom_colors, ArrayView<const Residue> residues, ArrayView<const BackboneSequence> bb_seq, ArrayView<const BackboneAngle> bb_angles,
                                 const Image& ramachandran_color_map);

void filter_colors(ArrayView<uint32> colors, Bitfield mask);
void desaturate_colors(ArrayView<uint32> colors, Bitfield mask, float scale);

vec3 magma_color_scale(float t);
vec3 inferno_color_scale(float t);
vec3 plasma_color_scale(float t);
vec3 viridis_color_scale(float t);

vec3 orange_color_scale(float t);
vec3 green_color_scale(float t);