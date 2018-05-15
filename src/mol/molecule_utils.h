#pragma once

#include <core/array.h>
#include <mol/molecule.h>
#include <mol/trajectory.h>
#include <mol/molecule_dynamic.h>
#include <mol/aminoacid.h>

enum class ColorMapping { STATIC_COLOR, CPK, RES_ID, RES_INDEX, CHAIN_ID, CHAIN_INDEX };

typedef bool (*FilterCommandFunc)(Array<bool> mask, const MoleculeDynamic& dynamic, Array<const CString> args);

struct FilterCommand {
    StringBuffer<16> keyword{};
    FilterCommandFunc func = nullptr;
};

enum class RamachandranConformationClassification {
    None,
    BetaHigh,
    BetaMid,
    BetaLow,
    AlphaHigh,
    AlphaMid,
    AlphaLow,
    LeftAlphaHigh,
    LeftAlphaMid,
    LeftAlphaLow,
    PMid,
    PLow
};

struct Volume {
    DynamicArray<uint8> data;
    ivec3 dim;
};

// tangent AND binormal is perhaps redundant
struct SplineSegment {
    vec3 position;
    vec3 tangent;
    vec3 normal;
    vec3 binormal;

    uint32 index;
    uint32 color;
};

struct BackboneAngles {
    float omega;
    float phi;
    float psi;
};

struct BackboneAnglesTrajectory {
    int num_segments = 0;
    int num_frames = 0;
    Array<BackboneAngles> angle_data{};
};

inline Array<BackboneAngles> get_backbone_angles(BackboneAnglesTrajectory& backbone_angle_traj, int frame_index) {
    if (backbone_angle_traj.angle_data.count == 0 || backbone_angle_traj.num_segments == 0) return {};
    ASSERT(frame_index < backbone_angle_traj.angle_data.count / backbone_angle_traj.num_segments);
    return Array<BackboneAngles>(&backbone_angle_traj.angle_data[frame_index * backbone_angle_traj.num_segments], backbone_angle_traj.num_segments);
}

inline Array<BackboneAngles> get_backbone_angles(BackboneAnglesTrajectory& backbone_angle_traj, int frame_offset, int frame_count) {
    if (backbone_angle_traj.angle_data.count == 0 || backbone_angle_traj.num_segments == 0) return {};
    int32 num_frames = (int32)backbone_angle_traj.angle_data.count / backbone_angle_traj.num_segments;
    ASSERT(frame_offset < num_frames);
    ASSERT(frame_offset + frame_count <= num_frames);
    return backbone_angle_traj.angle_data.sub_array(frame_offset * backbone_angle_traj.num_segments, frame_count * backbone_angle_traj.num_segments);
}

inline int32 get_backbone_angles_trajectory_current_frame_count(const BackboneAnglesTrajectory& backbone_angle_traj) {
    if (backbone_angle_traj.angle_data.count == 0 || backbone_angle_traj.num_segments == 0) return 0;
    return (int32)backbone_angle_traj.angle_data.count / backbone_angle_traj.num_segments;
}

inline Array<BackboneAngles> get_backbone_angles(BackboneAnglesTrajectory& backbone_angle_traj, int frame_index, Chain chain) {
    return get_backbone_angles(backbone_angle_traj, frame_index).sub_array(chain.beg_res_idx, chain.end_res_idx - chain.beg_res_idx);
}

void transform_positions(Array<vec3> positions, const mat4& transformation);
void compute_bounding_box(vec3* min_box, vec3* max_box, Array<const vec3> positions);

void linear_interpolation_periodic(Array<vec3> positions, Array<const vec3> prev_pos, Array<const vec3> next_pos, float t, mat3 sim_box);
void linear_interpolation(Array<vec3> positions, Array<const vec3> prev_pos, Array<const vec3> next_pos, float t);
void spline_interpolation_periodic(Array<vec3> positions, Array<const vec3> pos0, Array<const vec3> pos1, Array<const vec3> pos2,
                                   Array<const vec3> pos3, float t, mat3 sim_box);
void spline_interpolation(Array<vec3> positions, Array<const vec3> pos0, Array<const vec3> pos1, Array<const vec3> pos2, Array<const vec3> pos3,
                          float t);

DynamicArray<Bond> compute_covalent_bonds(Array<const vec3> atom_pos, Array<const Element> atom_elem, Array<const Residue> residues = {});
DynamicArray<Chain> compute_chains(Array<const Residue> residue, Array<const Bond> bonds, Array<const ResIdx> atom_residue_indices = {});
DynamicArray<BackboneSegment> compute_backbone_segments(Array<const Residue> residues, Array<const Label> atom_labels);
DynamicArray<SplineSegment> compute_spline(Array<const vec3> atom_pos, Array<const uint32> colors, Array<const BackboneSegment> backbone,
                                           int32 num_subdivisions = 1, float tension = 0.5f);

// Computes the dihedral angles within the backbone:
// omega = dihedral(CA[i-1], C[i-1], N[i], CA[i])
// phi   = dihedral( C[i-1], N[i],  CA[i],  C[i])
// psi   = dihedral( N[i],  CA[i],   C[i],  N[i+1])
// As seen here https://en.wikipedia.org/wiki/Ramachandran_plot.
DynamicArray<BackboneAngles> compute_backbone_angles(Array<const vec3> atom_pos, Array<const BackboneSegment> backbone_segments);
void compute_backbone_angles(Array<BackboneAngles> dst, Array<const vec3> atom_pos, Array<const BackboneSegment> backbone_segments);

void init_backbone_angles_trajectory(BackboneAnglesTrajectory* data, const MoleculeDynamic& dynamic);
void free_backbone_angles_trajectory(BackboneAnglesTrajectory* data);
void compute_backbone_angles_trajectory(BackboneAnglesTrajectory* bb_angle_traj, const MoleculeDynamic& dynamic);

DynamicArray<float> compute_atom_radii(Array<const Element> elements);
void compute_atom_radii(Array<float> radii_dst, Array<const Element> elements);

DynamicArray<uint32> compute_atom_colors(const MoleculeStructure& mol, ColorMapping mapping, uint32 static_color = 0xffffffff);
void compute_atom_colors(Array<uint32> color_dst, const MoleculeStructure& mol, ColorMapping mapping, uint32 static_color = 0xffffffff);

Volume compute_occupancy_volume(Array<vec3> atom_pos);
void compute_occupancy_volume(Volume* volume, Array<vec3> atom_pos);

Volume compute_occupancy_volume(Array<vec3> atom_pos, vec3 min_box, vec3 max_box);
void compute_occupancy_volume(Volume* volume, Array<vec3> atom_pos, vec3 min_box, vec3 max_box);

struct HydrogenBondDonor {
	AtomIdx donor_idx = 0;
	AtomIdx hydro_idx = 0;
};

typedef AtomIdx HydrogenBondAcceptor;

struct HydrogenBondCandidates {
	DynamicArray<HydrogenBondDonor> donors;
	DynamicArray<HydrogenBondAcceptor> acceptors;
};

struct HydrogenBond {
	AtomIdx acc_idx = 0;
	AtomIdx don_idx = 0;
};

struct HydrogenBondTrajectory {
	DynamicArray<HydrogenBond> bond_data{};
	DynamicArray<Array<HydrogenBond>> frame_bonds{};
};

void compute_hydrogen_bond_candidates(HydrogenBondCandidates* candidates, const MoleculeStructure& mol);
inline HydrogenBondCandidates compute_hydrogen_bond_candidates(const MoleculeStructure& mol) {
	HydrogenBondCandidates candidates;
	compute_hydrogen_bond_candidates(&candidates, mol);
	return candidates;
}

void compute_hydrogen_bonds(DynamicArray<HydrogenBond>* bonds, const MoleculeDynamic& dyn, int32 frame_idx, float dist_cutoff = 2.f, float angle_cutoff = 3.14159265f * 0.25f, const HydrogenBondCandidates* candidates = nullptr);
DynamicArray<HydrogenBond> compute_hydrogen_bonds(const MoleculeDynamic& dyn, int32 frame_idx, float dist_cutoff = 2.f, float angle_cutoff = 3.14159265f * 0.25f, const HydrogenBondCandidates* candidates = nullptr);

HydrogenBondTrajectory compute_hydrogen_bonds_trajectory(const MoleculeDynamic& dyn, float dist_cutoff, float angle_cutoff);

// bool filter_valid(CString filter);
// bool filter_colors(Array<uint32> color_dst, const MoleculeStructure& mol, CString filter);

namespace filter {
void initialize();
void shutdown();
bool compute_filter_mask(Array<bool> mask, const MoleculeDynamic& mol, CString filter);
void filter_colors(Array<uint32> colors, Array<bool> mask);
}  // namespace filter

inline bool is_amino_acid(Residue res) { return aminoacid::get_from_string(res.name) != AminoAcid::Unknown; }

namespace draw {
void initialize();
void shutdown();
void draw_vdw(Array<const vec3> atom_positions, Array<const float> atom_radii, Array<const uint32> atom_colors, const mat4& view_mat,
              const mat4& proj_mat, float radii_scale = 1.f);
void draw_licorice(Array<const vec3> atom_positions, Array<const Bond> atom_bonds, Array<const uint32> atom_colors, const mat4& view_mat,
                   const mat4& proj_mat, float radii_scale = 1.f);
void draw_ribbons(Array<const BackboneSegment> backbone_segments, Array<const Chain> chains, Array<const vec3> atom_positions,
                  Array<const uint32> atom_colors, const mat4& view_mat, const mat4& proj_mat, int num_subdivisions = 8, float tension = 0.5f,
                  float width_scale = 1.f, float thickness_scale = 1.f);

// DEBUG
void draw_backbone(Array<const BackboneSegment> backbone, Array<const vec3> atom_positions, const mat4& view_mat, const mat4& proj_mat);
void draw_spline(Array<const SplineSegment> spline, const mat4& view_mat, const mat4& proj_mat);

}  // namespace draw

namespace ramachandran {
void initialize();
void shutdown();

void clear_accumulation_texture();
// Radius is given as percentage of normalized texture space coordinates (1.0 = 1% of texture width and height)
void compute_accumulation_texture(Array<const BackboneAngles> angles, vec4 color, float radius = 1.f, float outline = 0.f);
uint32 get_accumulation_texture();
uint32 get_segmentation_texture();
}  // namespace ramachandran
