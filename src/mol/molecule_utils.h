#pragma once

#include <core/allocator.h>
#include <core/array.h>
#include <mol/molecule.h>

enum class ColorMapping { STATIC_COLOR, CPK, RES_ID, RES_INDEX, CHAIN_ID, CHAIN_INDEX };

// Tangent and binormal is perhaps redundant
struct SplineSegment {
	vec3 position;
	vec3 tangent;
	vec3 normal;
	vec3 binormal;
};

void transform_positions(Array<vec3> positions, const mat4& transformation);
void compute_bounding_box(vec3* min_box, vec3* max_box, const Array<vec3> positions);

void linear_interpolation_periodic(Array<vec3> positions, const Array<vec3> prev_pos, const Array<vec3> next_pos, float t, mat3 sim_box);
void linear_interpolation(Array<vec3> positions, const Array<vec3> prev_pos, const Array<vec3> next_pos, float t);

DynamicArray<Bond> compute_covalent_bonds(const Array<vec3> atom_pos, const Array<Element> atom_elem, const Array<Residue> residues = {});
DynamicArray<Chain> compute_chains(const Array<Residue> residue, const Array<Bond> bonds, const Array<int32> atom_residue_indices = {});
DynamicArray<BackboneSegment> compute_backbone(const Chain& chain, const Array<Residue> residues, const Array<Label> atom_labels);
DynamicArray<SplineSegment> compute_spline(const Array<vec3> atom_pos, const Array<BackboneSegment>& backbone, int num_subdivisions = 1);

DynamicArray<float> compute_atom_radii(const Array<Element> elements);
void compute_atom_radii(Array<float> radii_dst, const Array<Element> elements);

DynamicArray<uint32> compute_atom_colors(const MoleculeStructure& mol, ColorMapping mapping);
void compute_atom_colors(Array<uint32> color_dst, const MoleculeStructure& mol, ColorMapping mapping);

bool is_amino_acid(Residue res);

namespace draw {
void initialize();
void shutdown();
void draw_vdw(const Array<vec3> atom_positions, const Array<float> atom_radii, const Array<uint32> atom_colors, const mat4& view_mat,
              const mat4& proj_mat, float radii_scale = 1.f);
void draw_licorice(const Array<vec3> atom_positions, const Array<Bond> atom_bonds, const Array<uint32> atom_colors, const mat4& view_mat,
                   const mat4& proj_mat, float radii_scale = 1.f);
void draw_backbone(const Array<BackboneSegment> backbone, const Array<vec3> atom_positions, const mat4& view_mat, const mat4& proj_mat);
void draw_spline(const Array<SplineSegment> spline, const mat4& view_mat, const mat4& proj_mat);
}  // namespace draw
