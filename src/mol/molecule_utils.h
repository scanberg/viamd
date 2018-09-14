#pragma once

#include <core/array_types.h>
#include <core/math_utils.h>
#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>
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

struct DynamicBasis {
    AtomIdx origin_idx = -1;
    AtomIdx x_idx = -1;
    AtomIdx y_idx = -1;
    AtomIdx z_idx = -1;

    vec3 extent;
};

struct HydrogenBond {
    AtomIdx acc_idx = 0;
    AtomIdx don_idx = 0;
    AtomIdx hyd_idx = 0;
};

struct HydrogenBondTrajectory {
    DynamicArray<HydrogenBond> bond_data{};
    DynamicArray<Array<HydrogenBond>> frame_bonds{};
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
void compute_bounding_box(vec3* min_box, vec3* max_box, Array<const vec3> positions);
vec3 compute_com(Array<const vec3> positions, Array<const float> masses = {});

/*
inline mat(const DynamicBasis& basis, Array<const vec3> atom_positions) {
    mat4 mat;
    vec3 x = atom_positions[basis.x_idx] - atom_positions[basis.origin_idx];
    vec3 y = atom_positions[basis.y_idx] - atom_positions[basis.origin_idx];
    vec3 z = atom_positions[basis.z_idx] - atom_positions[basis.origin_idx];

    return mat;
}
*/

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

namespace hydrogen_bond {
int32 compute_acceptors(DynamicArray<HydrogenBondAcceptor>* acceptors, Array<const Element> elements);
DynamicArray<HydrogenBondAcceptor> compute_acceptors(Array<const Element> elements);

int32 compute_donors(DynamicArray<HydrogenBondAcceptor>* acceptors, Array<const Label> labels);
DynamicArray<HydrogenBondDonor> compute_donors(Array<const Label> labels);

int32 compute_bonds(DynamicArray<HydrogenBond>* bonds, Array<const HydrogenBondDonor> donors, Array<const HydrogenBondAcceptor> acceptors,
                    Array<const vec3> atom_positions, float dist_cutoff = 3.f, float angle_cutoff = 20.f * math::DEG_TO_RAD);
DynamicArray<HydrogenBond> compute_bonds(Array<const HydrogenBondDonor> donors, Array<const HydrogenBondAcceptor> acceptors,
                                         Array<const vec3> atom_positions, float dist_cutoff = 3.f, float angle_cutoff = 20.f * math::DEG_TO_RAD);

void compute_bonds_trajectory(HydrogenBondTrajectory* hbt, const MoleculeDynamic& dyn, float dist_cutoff, float angle_cutoff);
HydrogenBondTrajectory compute_bonds_trajectory(const MoleculeDynamic& dyn, float dist_cutoff, float angle_cutoff);
}  // namespace hydrogen_bond

// bool filter_valid(CString filter);
// bool filter_colors(Array<uint32> color_dst, const MoleculeStructure& mol, CString filter);

namespace filter {
void initialize();
void shutdown();
bool compute_filter_mask(Array<bool> mask, const MoleculeDynamic& dynamic, CString filter);
void filter_colors(Array<uint32> colors, Array<bool> mask);

template <typename T>
void extract_filtered_data(DynamicArray<T>* dst_data, Array<const T> src_data, Array<const bool> mask) {
    ASSERT(dst_data);
    dst_data->clear();
    for (int32 i = 0; i < mask.count; i++) {
        if (mask[i]) dst_data->push_back(src_data[i]);
    }
}

template <typename T>
DynamicArray<T> extract_filtered_data(Array<const T> data, Array<const bool> mask) {
    DynamicArray<T> result{};
    extract_filtered_data(&result, data, mask);
    return result;
}

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

// Compute a linear transform which fits a set of points_a onto a set of points_b in least squares fashion.
mat4 compute_linear_transform(Array<const vec3> pos_frame_a, Array<const vec3> pos_frame_b);
mat4 compute_linear_transform(Array<const vec3> pos_frame_a, Array<const vec3> pos_frame_b, Array<const float> mass);
void compute_RS(mat3* R, mat3* S, Array<const vec3> x0, Array<const vec3> x, Array<const float> mass = {});

// 3x3 matrix operations
void diagonalize(const mat3& M, mat3* Q, mat3* D);
void decompose(const mat3& M, mat3* R, mat3* S);
