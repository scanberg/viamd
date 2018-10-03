#pragma once

#include <core/array_types.h>
#include <core/math_utils.h>
#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>
#include <mol/molecule_dynamic.h>
#include <mol/aminoacid.h>

enum class ColorMapping { STATIC_COLOR, CPK, RES_ID, RES_INDEX, CHAIN_ID, CHAIN_INDEX };

struct DynamicBasis {
    AtomIdx origin_idx = -1;
    AtomIdx x_idx = -1;
    AtomIdx y_idx = -1;
    AtomIdx z_idx = -1;

    vec3 extent;
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
#ifdef DEBUG
    int32 num_frames = (int32)backbone_angle_traj.angle_data.count / backbone_angle_traj.num_segments;
    ASSERT(frame_offset < num_frames);
    ASSERT(frame_offset + frame_count <= num_frames);
#endif
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
void compute_bounding_box(vec3* min_box, vec3* max_box, Array<const vec3> positions, Array<const float> radii = {});
vec3 compute_com(Array<const vec3> positions, Array<const float> masses = {});

void linear_interpolation_periodic(Array<vec3> positions, Array<const vec3> prev_pos, Array<const vec3> next_pos, float t, mat3 sim_box);
void linear_interpolation(Array<vec3> positions, Array<const vec3> prev_pos, Array<const vec3> next_pos, float t);
void spline_interpolation_periodic(Array<vec3> positions, Array<const vec3> pos0, Array<const vec3> pos1, Array<const vec3> pos2,
                                   Array<const vec3> pos3, float t, mat3 sim_box);
void spline_interpolation(Array<vec3> positions, Array<const vec3> pos0, Array<const vec3> pos1, Array<const vec3> pos2, Array<const vec3> pos3,
                          float t);

DynamicArray<Bond> compute_covalent_bonds(Array<const vec3> atom_pos, Array<const Element> atom_elem, Array<const ResIdx> atom_res_idx = {});
DynamicArray<Chain> compute_chains(Array<const Residue> residue, Array<const Bond> bonds, Array<const ResIdx> atom_res_idx = {});
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

inline bool is_amino_acid(Residue res) { return aminoacid::get_from_string(res.name) != AminoAcid::Unknown; }

// Compute a linear transform which fits a set of points_a onto a set of points_b in least squares fashion.
mat4 compute_linear_transform(Array<const vec3> pos_frame_a, Array<const vec3> pos_frame_b);
mat4 compute_linear_transform(Array<const vec3> pos_frame_a, Array<const vec3> pos_frame_b, Array<const float> mass);
void compute_RS(mat3* R, mat3* S, Array<const vec3> x0, Array<const vec3> x, Array<const float> mass = {});

// 3x3 matrix operations
void diagonalize(const mat3& M, mat3* Q, mat3* D);
void decompose(const mat3& M, mat3* R, mat3* S);
