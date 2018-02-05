#pragma once

#include <core/allocator.h>

enum class ColorMapping {
    STATIC_COLOR,
    CPK,
    RES_ID,
    RES_INDEX,
    CHAIN_ID,
    CHAIN_INDEX
};

void transform_positions(Array<vec3> positions, const mat4 transformation);

Array<vec3>     compute_bonds(const Array<vec3> atom_pos, const Array<Element> atom_elem, Allocator& alloc = default_alloc);
Array<Backbone> compute_backbones(const Array<Residue> residues, const Array<Bond> bonds, Allocator& alloc = default_alloc);
Array<float>    compute_atom_radii(const Array<Element> elements, Allocator& alloc = default_alloc);
Array<uint32>   compute_atom_colors(MoleculeStructure& mol, ColorMapping mapping);