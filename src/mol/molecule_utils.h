#pragma once

#include <core/array.h>
#include <mol/molecule.h>
#include <mol/trajectory.h>
#include <mol/molecule_dynamic.h>
#include <mol/aminoacid.h>

enum class ColorMapping { STATIC_COLOR, CPK, RES_ID, RES_INDEX, CHAIN_ID, CHAIN_INDEX };

typedef bool (*FilterCommandFunc)(Array<bool> mask, const MoleculeDynamic& dynamic, const Array<CString> args);

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
    return backbone_angle_traj.angle_data.sub_array(frame_offset * backbone_angle_traj.num_segments,
                                                    frame_count * backbone_angle_traj.num_segments);
}

inline int32 get_backbone_angles_trajectory_current_frame_count(const BackboneAnglesTrajectory& backbone_angle_traj) {
	if (backbone_angle_traj.angle_data.count == 0 || backbone_angle_traj.num_segments == 0) return 0;
	return (int32)backbone_angle_traj.angle_data.count / backbone_angle_traj.num_segments;
}

inline Array<BackboneAngles> get_backbone_angles(BackboneAnglesTrajectory& backbone_angle_traj, int frame_index, Chain chain) {
    return get_backbone_angles(backbone_angle_traj, frame_index).sub_array(chain.beg_res_idx, chain.end_res_idx - chain.beg_res_idx);
}

void transform_positions(Array<vec3> positions, const mat4& transformation);
void compute_bounding_box(vec3* min_box, vec3* max_box, const Array<vec3> positions);

void linear_interpolation_periodic(Array<vec3> positions, const Array<vec3> prev_pos, const Array<vec3> next_pos, float t, mat3 sim_box);
void linear_interpolation(Array<vec3> positions, const Array<vec3> prev_pos, const Array<vec3> next_pos, float t);

inline float dihedral_angle(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3) {
    vec3 b1 = p1 - p0;
    vec3 b2 = p2 - p1;
    vec3 b3 = p3 - p2;
    vec3 c1 = glm::cross(b1, b2);
    vec3 c2 = glm::cross(b2, b3);
    return glm::atan(glm::dot(glm::cross(c1, c2), glm::normalize(b2)), glm::dot(c1, c2));
}

inline float dihedral_angle(const vec3 p[4]) { return dihedral_angle(p[0], p[1], p[2], p[3]); }

DynamicArray<Bond> compute_covalent_bonds(const Array<vec3> atom_pos, const Array<Element> atom_elem, const Array<Residue> residues = {});
DynamicArray<Chain> compute_chains(const Array<Residue> residue, const Array<Bond> bonds, const Array<ResIdx> atom_residue_indices = {});
DynamicArray<BackboneSegment> compute_backbone_segments(const Array<Residue> residues, const Array<Label> atom_labels);
DynamicArray<SplineSegment> compute_spline(const Array<vec3> atom_pos, const Array<uint32> colors, const Array<BackboneSegment>& backbone,
                                           int32 num_subdivisions = 1);

// Computes the dihedral angles within the backbone:
// omega = dihedral(CA[i-1], C[i-1], N[i], CA[i])
// phi   = dihedral( C[i-1], N[i],  CA[i],  C[i])
// psi   = dihedral( N[i],  CA[i],   C[i],  N[i+1])
// As seen here https://en.wikipedia.org/wiki/Ramachandran_plot.
DynamicArray<BackboneAngles> compute_backbone_angles(const Array<vec3> atom_pos, const Array<BackboneSegment> backbone_segments);
void compute_backbone_angles(Array<BackboneAngles> dst, const Array<vec3> atom_pos, const Array<BackboneSegment> backbone_segments);

void init_backbone_angles_trajectory(BackboneAnglesTrajectory* data, const MoleculeDynamic& dynamic);
void free_backbone_angles_trajectory(BackboneAnglesTrajectory* data);
void compute_backbone_angles_trajectory(BackboneAnglesTrajectory* bb_angle_traj, const MoleculeDynamic& dynamic);

DynamicArray<float> compute_atom_radii(const Array<Element> elements);
void compute_atom_radii(Array<float> radii_dst, const Array<Element> elements);

DynamicArray<uint32> compute_atom_colors(const MoleculeStructure& mol, ColorMapping mapping, uint32 static_color = 0xffffffff);
void compute_atom_colors(Array<uint32> color_dst, const MoleculeStructure& mol, ColorMapping mapping, uint32 static_color = 0xffffffff);

Volume compute_occupancy_volume(Array<vec3> atom_pos);
void compute_occupancy_volume(Volume* volume, Array<vec3> atom_pos);

Volume compute_occupancy_volume(Array<vec3> atom_pos, vec3 min_box, vec3 max_box);
void compute_occupancy_volume(Volume* volume, Array<vec3> atom_pos, vec3 min_box, vec3 max_box);

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
void draw_vdw(const Array<vec3> atom_positions, const Array<float> atom_radii, const Array<uint32> atom_colors, const mat4& view_mat,
              const mat4& proj_mat, float radii_scale = 1.f);
void draw_licorice(const Array<vec3> atom_positions, const Array<Bond> atom_bonds, const Array<uint32> atom_colors, const mat4& view_mat,
                   const mat4& proj_mat, float radii_scale = 1.f);
void draw_ribbons(const Array<BackboneSegment> backbone_segments, const Array<Chain> chains, const Array<vec3> atom_positions,
                  const Array<uint32> atom_colors, const mat4& view_mat, const mat4& proj_mat, int num_subdivisions = 8);

// DEBUG
void draw_backbone(const Array<BackboneSegment> backbone, const Array<vec3> atom_positions, const mat4& view_mat, const mat4& proj_mat);
void draw_spline(const Array<SplineSegment> spline, const mat4& view_mat, const mat4& proj_mat);

}  // namespace draw

namespace ramachandran {
void initialize();
void shutdown();
// Radius is given as percentage of normalized texture space coordinates (1.0 = 1% of texture width and height)
void compute_accumulation_texture(const Array<BackboneAngles> angles, const Array<BackboneAngles> highlighted_angles, float radius = 1.f,
                                  float opacity = 1.f);
uint32 get_accumulation_texture();
uint32 get_segmentation_texture();
}  // namespace ramachandran
