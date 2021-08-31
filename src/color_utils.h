#pragma once

#include <core/types.h>
#include <core/array_types.h>
#include <core/bitfield.h>
#include <core/math_utils.h>

#include <md_molecule.h>
#include "image.h"

enum class ColorMapping { Uniform, Cpk, ResId, ResIndex, ChainId, ChainIndex, SecondaryStructure };

inline vec4 color_from_hash(u32 hash) {
    constexpr float chroma = 0.8f;
    constexpr float luminance = 1.0f;
    constexpr u32 mod = 21;
    const float hue = (hash % mod) / (float)mod;
    const vec3 rgb = math::hcl_to_rgb(vec3(hue, chroma, luminance));

    return vec4(rgb, 1);
}

void color_atoms_uniform(Array<u32> dst_atom_colors, const vec4& color);
void color_atoms_cpk(Array<u32> dst_atom_colors, const md_molecule_t& mol);
void color_atoms_residue_id(Array<u32> dst_atom_colors, const md_molecule_t& mol);
void color_atoms_residue_index(Array<u32> dst_atom_colors, const md_molecule_t& mol);
void color_atoms_chain_id(Array<u32> dst_atom_colors, const md_molecule_t& mol);
void color_atoms_chain_index(Array<u32> dst_atom_colors, const md_molecule_t& mol);
void color_atoms_secondary_structure(Array<u32> dst_atom_colors, const md_molecule_t& mol);

void color_atoms_backbone_angles(Array<u32> dst_atom_colors, const md_molecule_t& mol, const Image& ramachandran_color_map);

void filter_colors(Array<u32> colors, Bitfield mask);
//void desaturate_colors(Array<u32> colors, Bitfield mask, float scale);

vec3 magma_color_scale(float t);
vec3 inferno_color_scale(float t);
vec3 plasma_color_scale(float t);
vec3 viridis_color_scale(float t);

vec3 orange_color_scale(float t);
vec3 green_color_scale(float t);

vec4 qualitative_color_scale(int idx);