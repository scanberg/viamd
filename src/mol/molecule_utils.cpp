//#define NOMINMAX

#include "molecule_utils.h"
#include <core/common.h>
#include <core/hash.h>
#include <core/log.h>
#include <mol/trajectory_utils.h>
#include <mol/spatial_hash.h>
#include <gfx/gl_utils.h>
#include <gfx/immediate_draw_utils.h>

//#include <imgui.h>
#include <ctype.h>
#include <fstream>

void transform_positions(Array<vec3> positions, const mat4& transformation) {
    for (auto& p : positions) {
        p = vec3(transformation * vec4(p, 1));
    }
}

void compute_bounding_box(vec3* min_box, vec3* max_box, Array<const vec3> positions) {
    ASSERT(min_box);
    ASSERT(max_box);

    if (positions.count == 0) {
        *min_box = *max_box = vec3(0);
    }

    *min_box = *max_box = positions.data[0];
    for (int64 i = 0; i < positions.count; i++) {
        *min_box = math::min(*min_box, positions.data[i]);
        *max_box = math::max(*max_box, positions.data[i]);
    }
}

vec3 compute_com(Array<const vec3> positions, Array<const float> masses) {
    if (positions.count == 0) return {0, 0, 0};
    if (positions.count == 1) return positions[0];

    vec3 com{0};
    if (masses.count == 0) {
        for (const auto& p : positions) {
            com += p;
        }
        com = com / (float)positions.count;
    } else {
        ASSERT(masses.count == positions.count);
        vec3 pos_mass_sum{0};
        float mass_sum{0};
        for (int32 i = 0; i < positions.count; i++) {
            pos_mass_sum += positions[i] * masses[i];
            mass_sum += masses[i];
        }
        com = pos_mass_sum / mass_sum;
    }

    return com;
}

inline bool periodic_jump(const vec3& p_prev, const vec3& p_next, const vec3& half_box) {
    const vec3 abs_delta = math::abs(p_next - p_prev);
    if (abs_delta.x > half_box.x) return true;
    if (abs_delta.y > half_box.y) return true;
    if (abs_delta.z > half_box.z) return true;
    return false;
}

void linear_interpolation(Array<vec3> positions, Array<const vec3> prev_pos, Array<const vec3> next_pos, float t) {
    ASSERT(prev_pos.count == positions.count);
    ASSERT(next_pos.count == positions.count);

    for (int i = 0; i < positions.count; i++) {
        positions[i] = math::mix(prev_pos[i], next_pos[i], t);
    }
}

// @TODO: Fix this, is it possible in theory to get a good interpolation between frames with periodicity without modifying source data?
// @PERFORMANCE: VECTORIZE THE LIVING SHIET OUT OF THIS
void linear_interpolation_periodic(Array<vec3> positions, Array<const vec3> prev_pos, Array<const vec3> next_pos, float t, mat3 sim_box) {
    ASSERT(prev_pos.count == positions.count);
    ASSERT(next_pos.count == positions.count);

    const vec3 full_box_ext = sim_box * vec3(1);
    const vec3 half_box_ext = full_box_ext * 0.5f;

    for (int i = 0; i < positions.count; i++) {
        vec3 next = next_pos[i];
        vec3 prev = prev_pos[i];

        if (periodic_jump(prev, next, half_box_ext)) {
            // Atom moved across periodic boundry we cannot apply linearly interpolate directly
            positions[i] = t < 0.5f ? prev : next;
        } else {
            positions[i] = math::mix(prev, next, t);
        }
    }
}

void spline_interpolation_periodic(Array<vec3> positions, Array<const vec3> pos0, Array<const vec3> pos1, Array<const vec3> pos2,
                                   Array<const vec3> pos3, float t, mat3 sim_box) {
    ASSERT(pos0.count == positions.count);
    ASSERT(pos1.count == positions.count);
    ASSERT(pos2.count == positions.count);
    ASSERT(pos3.count == positions.count);

    const vec3 full_box_ext = sim_box * vec3(1);
    const vec3 half_box_ext = full_box_ext * 0.5f;

    for (int i = 0; i < positions.count; i++) {
        vec3 p0 = pos0[i];
        vec3 p1 = pos1[i];
        vec3 p2 = pos2[i];
        vec3 p3 = pos3[i];

        if (periodic_jump(p0, p1, half_box_ext)) p0 = p1;
        if (periodic_jump(p2, p3, half_box_ext)) p3 = p2;

        if (periodic_jump(p1, p2, half_box_ext)) {
            positions[i] = t < 0.5f ? p1 : p2;
        } else {
            positions[i] = math::spline(p0, p1, p2, p3, t);
        }
    }
}

void spline_interpolation(Array<vec3> positions, Array<const vec3> pos0, Array<const vec3> pos1, Array<const vec3> pos2, Array<const vec3> pos3,
                          float t) {
    ASSERT(pos0.count == positions.count);
    ASSERT(pos1.count == positions.count);
    ASSERT(pos2.count == positions.count);
    ASSERT(pos3.count == positions.count);

    for (int i = 0; i < positions.count; i++) {
        vec3 p0 = pos0[i];
        vec3 p1 = pos1[i];
        vec3 p2 = pos2[i];
        vec3 p3 = pos3[i];

        positions[i] = math::spline(p0, p1, p2, p3, t);
    }
}

// Computes covalent bonds between a set of atoms with given positions and elements.
// The approach is inspired by the technique used in NGL (https://github.com/arose/ngl)
DynamicArray<Bond> compute_covalent_bonds(Array<const vec3> atom_pos, Array<const Element> atom_elem, Array<const Residue> residues) {
    ASSERT(atom_pos.count == atom_elem.count);

    DynamicArray<Bond> bonds;

    auto try_and_create_atom_bond = [& pos = atom_pos, &elem = atom_elem, &bonds = bonds](int i, int j) -> bool {
        auto d = element::covalent_radius(elem[i]) + element::covalent_radius(elem[j]);
        auto d1 = d + 0.3f;
        auto d2 = d - 0.5f;
        auto v = pos[i] - pos[j];
        auto dist2 = glm::dot(v, v);

        if (dist2 < (d1 * d1) && dist2 > (d2 * d2)) {
            bonds.push_back({i, j});
            return true;
        } else
            return false;
    };

    if (residues) {
        // Create connections within residues
        for (int64 i = 0; i < residues.count; i++) {
            const Residue& res = residues[i];
            for (int atom_i = res.beg_atom_idx; atom_i < res.end_atom_idx - 1; atom_i++) {
                for (int atom_j = atom_i + 1; atom_j < res.end_atom_idx; atom_j++) {
                    try_and_create_atom_bond(atom_i, atom_j);
                }
            }
        }

        // @TODO:
        // Create connections between consecutive residues (Is it enough and correct to only check
        // consecutive residues?)

        for (int res_i = 0; res_i < residues.count - 1; res_i++) {
            auto res_j = res_i + 1;
            auto atom_i_beg = residues[res_i].beg_atom_idx;
            auto atom_i_end = residues[res_i].end_atom_idx;
            auto atom_j_beg = residues[res_j].beg_atom_idx;
            auto atom_j_end = residues[res_j].end_atom_idx;

            for (int atom_i = atom_i_beg; atom_i < atom_i_end; atom_i++) {
                for (int atom_j = atom_j_beg; atom_j < atom_j_end; atom_j++) {
                    try_and_create_atom_bond(atom_i, atom_j);
                }
            }
        }
    } else {
        // Brute force N^2 check
        // @TODO: Use spatial hash
        auto atom_count = atom_pos.count;
        for (int atom_i = 0; atom_i < atom_count - 1; ++atom_i) {
            for (int atom_j = 1; atom_j < atom_count; ++atom_j) {
                try_and_create_atom_bond(atom_i, atom_j);
            }
        }
    }

    return bonds;
}

// @NOTE this method is sub-optimal and can surely be improved...
// Residues should have no more than 2 potential connections to other residues.
DynamicArray<Chain> compute_chains(Array<const Residue> residues, Array<const Bond> bonds, Array<const ResIdx> atom_residue_indices) {
    DynamicArray<Bond> residue_bonds;

    if (atom_residue_indices) {
        for (const auto& bond : bonds) {
            if (atom_residue_indices[bond.idx_a] != atom_residue_indices[bond.idx_b]) {
                residue_bonds.push_back({atom_residue_indices[bond.idx_a], atom_residue_indices[bond.idx_b]});
            }
        }
    } else {
        ASSERT(false, "Not implemented Yeti");
    }

    if (residue_bonds.count == 0) {
        // No residue bonds, No chains.
        return {};
    }

    DynamicArray<int> residue_chains(residues.count, -1);

    if (residue_bonds.count > 0) {
        int curr_chain_idx = 0;
        int res_bond_idx = 0;
        for (int i = 0; i < residues.count; i++) {
            if (residue_chains[i] == -1) residue_chains[i] = curr_chain_idx++;
            for (; res_bond_idx < residue_bonds.count; res_bond_idx++) {
                const auto& res_bond = residue_bonds[res_bond_idx];
                if (i == res_bond.idx_a) {
                    residue_chains[res_bond.idx_b] = residue_chains[res_bond.idx_a];
                } else if (res_bond.idx_a > i)
                    break;
            }
        }
    }

    DynamicArray<Chain> chains;
    int curr_chain_idx = -1;
    for (int i = 0; i < residue_chains.count; i++) {
        if (residue_chains[i] != curr_chain_idx) {
            curr_chain_idx = residue_chains[i];
            Label lbl;
            snprintf(lbl.beg(), Label::MAX_LENGTH, "C%i", curr_chain_idx);
            chains.push_back({lbl, (ResIdx)i, (ResIdx)i});
        }
        chains.back().end_res_idx++;
    }

    return chains;
}

template <int64 N>
bool match(const Label& lbl, const char (&cstr)[N]) {
    for (int64 i = 0; i < N; i++) {
        if (tolower(lbl[i]) != tolower(cstr[i])) return false;
    }
    return true;
}

DynamicArray<BackboneSegment> compute_backbone_segments(Array<const Residue> residues, Array<const Label> atom_labels) {
    DynamicArray<BackboneSegment> segments;
    for (auto& res : residues) {
        auto ca_idx = -1;
        auto n_idx = -1;
        auto c_idx = -1;
        auto o_idx = -1;
        if (is_amino_acid(res)) {
            // find atoms
            for (int32 i = res.beg_atom_idx; i < res.end_atom_idx; i++) {
                const auto& lbl = atom_labels[i];
                if (ca_idx == -1 && match(lbl, "CA")) ca_idx = i;
                if (n_idx == -1 && match(lbl, "N")) n_idx = i;
                if (c_idx == -1 && match(lbl, "C")) c_idx = i;
                if (o_idx == -1 && match(lbl, "O")) o_idx = i;
            }

            // Could not match "O"
            if (o_idx == -1) {
                // Pick first atom containing O after C atom
                for (int32 i = c_idx; i < res.end_atom_idx; i++) {
                    const auto& lbl = atom_labels[i];
                    if (lbl[0] == 'o' || lbl[0] == 'O') o_idx = i;
                }
            }

            if (ca_idx == -1 || n_idx == -1 || c_idx == -1 || o_idx == -1) {
                LOG_ERROR("Could not identify all backbone indices for residue %s.\n", res.name.beg());
            }
            segments.push_back({ca_idx, n_idx, c_idx, o_idx});
        } else {
            segments.push_back({-1, -1, -1, -1});
        }
    }

    return segments;
}

DynamicArray<SplineSegment> compute_spline(Array<const vec3> atom_pos, Array<const uint32> colors, Array<const BackboneSegment> backbone,
                                           int num_subdivisions, float tension) {
    if (backbone.count < 4) return {};

    DynamicArray<vec3> p_tmp;
    DynamicArray<vec3> o_tmp;
    DynamicArray<vec3> c_tmp;
    DynamicArray<int> ca_idx;

    auto d_p0 = atom_pos[backbone[1].ca_idx] - atom_pos[backbone[0].ca_idx];
    p_tmp.push_back(atom_pos[backbone[0].ca_idx] - d_p0);

    auto d_o0 = atom_pos[backbone[1].o_idx] - atom_pos[backbone[0].o_idx];
    o_tmp.push_back(atom_pos[backbone[0].o_idx] - d_o0);

    auto d_c0 = atom_pos[backbone[1].c_idx] - atom_pos[backbone[0].c_idx];
    c_tmp.push_back(atom_pos[backbone[0].c_idx] - d_c0);

    ca_idx.push_back(backbone[0].ca_idx);

    const int size = (int)(backbone.count);
    for (auto i = 0; i < size; i++) {
        p_tmp.push_back(atom_pos[backbone[i].ca_idx]);
        o_tmp.push_back(atom_pos[backbone[i].o_idx]);
        c_tmp.push_back(atom_pos[backbone[i].c_idx]);
        ca_idx.push_back(backbone[i].ca_idx);
    }

    auto d_pn = atom_pos[backbone[size - 1].ca_idx] - atom_pos[backbone[size - 2].ca_idx];
    p_tmp.push_back(atom_pos[backbone[size - 1].ca_idx] + d_pn);
    p_tmp.push_back(p_tmp.back() + d_pn);

    auto d_on = atom_pos[backbone[size - 1].o_idx] - atom_pos[backbone[size - 2].o_idx];
    o_tmp.push_back(atom_pos[backbone[size - 1].o_idx] + d_on);
    o_tmp.push_back(o_tmp.back() + d_on);

    auto d_cn = atom_pos[backbone[size - 1].c_idx] - atom_pos[backbone[size - 2].c_idx];
    c_tmp.push_back(atom_pos[backbone[size - 1].c_idx] + d_cn);
    c_tmp.push_back(c_tmp.back() + d_cn);

    ca_idx.push_back(backbone[size - 1].ca_idx);
    ca_idx.push_back(backbone[size - 1].ca_idx);

    // NEEDED?
    for (int64 i = 1; i < o_tmp.size(); i++) {
        vec3 v0 = o_tmp[i - 1] - c_tmp[i - 1];
        vec3 v1 = o_tmp[i] - c_tmp[i];

        if (glm::dot(v0, v1) < 0) {
            o_tmp[i] = c_tmp[i] - v1;
        }
    }

    DynamicArray<SplineSegment> segments;

    for (int64 i = 1; i < p_tmp.size() - 2; i++) {
        auto p0 = p_tmp[i - 1];
        auto p1 = p_tmp[i];
        auto p2 = p_tmp[i + 1];
        auto p3 = p_tmp[i + 2];

        auto o0 = o_tmp[i - 1];
        auto o1 = o_tmp[i];
        auto o2 = o_tmp[i + 1];
        auto o3 = o_tmp[i + 2];

        auto c0 = c_tmp[i - 1];
        auto c1 = c_tmp[i];
        auto c2 = c_tmp[i + 1];
        auto c3 = c_tmp[i + 2];

        uint32 idx = ca_idx[i];
        uint32 color = colors[idx];

        auto count = (i < (p_tmp.size() - 3)) ? num_subdivisions : num_subdivisions + 1;
        for (int n = 0; n < count; n++) {
            auto t = n / (float)(num_subdivisions);

            vec3 p = math::spline(p0, p1, p2, p3, t, tension);
            vec3 o = math::spline(o0, o1, o2, o3, t, tension);
            vec3 c = math::spline(c0, c1, c2, c3, t, tension);

            vec3 v_dir = math::normalize(o - c);

            const float eps = 0.0001f;
            float d0 = math::max(0.f, t - eps);
            float d1 = math::min(t + eps, 1.f);

            vec3 tangent = math::normalize(math::spline(p0, p1, p2, p3, d1, tension) - math::spline(p0, p1, p2, p3, d0, tension));
            vec3 normal = math::normalize(math::cross(v_dir, tangent));
            vec3 binormal = math::normalize(math::cross(tangent, normal));
            // vec3 normal = v_dir;
            // vec3 binormal = math::normalize(math::cross(tangent, normal));

            segments.push_back({p, tangent, normal, binormal, idx, color});
        }
    }

    return segments;
}

DynamicArray<BackboneAngles> compute_backbone_angles(Array<const vec3> pos, Array<const BackboneSegment> backbone) {
    if (backbone.count == 0) return {};
    DynamicArray<BackboneAngles> angles(backbone.count);
    compute_backbone_angles(angles, pos, backbone);
    return angles;
}

void compute_backbone_angles(Array<BackboneAngles> dst, Array<const vec3> pos, Array<const BackboneSegment> backbone_segments) {
    ASSERT(dst.count >= backbone_segments.count);
    float omega, phi, psi;

    omega = 0;
    phi = 0;
    psi = math::dihedral_angle(pos[backbone_segments[0].n_idx], pos[backbone_segments[0].ca_idx], pos[backbone_segments[0].c_idx],
                               pos[backbone_segments[1].n_idx]);
    dst[0] = {omega, phi, psi};

    for (int64 i = 1; i < backbone_segments.count - 1; i++) {
        omega = math::dihedral_angle(pos[backbone_segments[i - 1].ca_idx], pos[backbone_segments[i - 1].c_idx], pos[backbone_segments[i].n_idx],
                                     pos[backbone_segments[i].ca_idx]);
        phi = math::dihedral_angle(pos[backbone_segments[i - 1].c_idx], pos[backbone_segments[i].n_idx], pos[backbone_segments[i].ca_idx],
                                   pos[backbone_segments[i].c_idx]);
        psi = math::dihedral_angle(pos[backbone_segments[i].n_idx], pos[backbone_segments[i].ca_idx], pos[backbone_segments[i].c_idx],
                                   pos[backbone_segments[i + 1].n_idx]);
        dst[i] = {omega, phi, psi};
    }

    auto N = backbone_segments.count - 1;
    omega = math::dihedral_angle(pos[backbone_segments[N - 1].ca_idx], pos[backbone_segments[N - 1].c_idx], pos[backbone_segments[N].n_idx],
                                 pos[backbone_segments[N].ca_idx]);
    phi = math::dihedral_angle(pos[backbone_segments[N - 1].c_idx], pos[backbone_segments[N].n_idx], pos[backbone_segments[N].ca_idx],
                               pos[backbone_segments[N].c_idx]);
    psi = 0;
    dst[N] = {omega, phi, psi};
}

void init_backbone_angles_trajectory(BackboneAnglesTrajectory* data, const MoleculeDynamic& dynamic) {
    ASSERT(data);
    if (!dynamic.molecule || !dynamic.trajectory) return;

    if (data->angle_data) {
        FREE(data->angle_data.data);
    }

    int32 alloc_count = (int32)dynamic.molecule.backbone_segments.count * (int32)dynamic.trajectory.frame_buffer.count;
    data->num_segments = (int32)dynamic.molecule.backbone_segments.count;
    data->num_frames = 0;
    data->angle_data = {(BackboneAngles*)CALLOC(alloc_count, sizeof(BackboneAngles)), alloc_count};
}

void free_backbone_angles_trajectory(BackboneAnglesTrajectory* data) {
    ASSERT(data);
    if (data->angle_data) {
        FREE(data->angle_data.data);
        *data = {};
    }
}

void compute_backbone_angles_trajectory(BackboneAnglesTrajectory* data, const MoleculeDynamic& dynamic) {
    ASSERT(dynamic.trajectory && dynamic.molecule);
    if (dynamic.trajectory.num_frames == 0 || dynamic.molecule.backbone_segments.count == 0) return;

    //@NOTE: Trajectory may be loading while this is taking place, therefore read num_frames once and stick to that
    int32 traj_num_frames = dynamic.trajectory.num_frames;

    // @NOTE: If we are up to date, no need to compute anything
    if (traj_num_frames == data->num_frames) {
        return;
    }

    // @TODO: parallelize?
    // @NOTE: Only compute data for indices which are new
    for (int32 f_idx = data->num_frames; f_idx < traj_num_frames; f_idx++) {
        Array<const vec3> frame_pos = get_trajectory_positions(dynamic.trajectory, f_idx);
        Array<BackboneAngles> frame_angles = get_backbone_angles(*data, f_idx);
        for (const Chain& c : dynamic.molecule.chains) {
            auto bb_segments = get_backbone(dynamic.molecule, c);
            auto bb_angles = frame_angles.sub_array(c.beg_res_idx, c.end_res_idx - c.beg_res_idx);
            compute_backbone_angles(bb_angles, frame_pos, bb_segments);
        }
    }
}

DynamicArray<float> compute_atom_radii(Array<const Element> elements) {
    DynamicArray<float> radii(elements.count, 0);
    compute_atom_radii(radii, elements);
    return radii;
}

void compute_atom_radii(Array<float> radii_dst, Array<const Element> elements) {
    ASSERT(radii_dst.count <= elements.count);
    for (int64 i = 0; i < radii_dst.count; i++) {
        radii_dst[i] = element::vdw_radius(elements[i]);
    }
}

DynamicArray<uint32> compute_atom_colors(const MoleculeStructure& mol, ColorMapping mapping, uint32 static_color) {
    DynamicArray<uint32> colors(mol.atom_elements.count, 0xFFFFFFFF);
    compute_atom_colors(colors, mol, mapping, static_color);
    return colors;
}

static inline vec3 compute_color(uint32 hash) {
    constexpr float CHROMA = 0.45f;
    constexpr float LUMINANCE = 0.90f;
    constexpr int32 MOD = 21;
    constexpr float SCL = 1.f / (float)MOD;

    return math::hcl_to_rgb(vec3((hash % MOD) * SCL, CHROMA, LUMINANCE));
}

void compute_atom_colors(Array<uint32> color_dst, const MoleculeStructure& mol, ColorMapping mapping, uint32 static_color) {
    // @TODO: Implement more mappings

    // CPK
    switch (mapping) {
        case ColorMapping::STATIC_COLOR:
            for (int64 i = 0; i < color_dst.count; i++) {
                color_dst[i] = static_color;
            }
            break;
        case ColorMapping::CPK:
            for (int64 i = 0; i < color_dst.count; i++) {
                color_dst[i] = element::color(mol.atom_elements[i]);
            }
            break;
        case ColorMapping::RES_ID:
            // Color based on residues, not unique by any means.
            // Perhaps use predefined colors if residues are Amino acids
            for (int64 i = 0; i < color_dst.count; i++) {
                if (i < mol.atom_residue_indices.count) {
                    const auto& res = mol.residues[mol.atom_residue_indices[i]];
                    vec3 c = compute_color(hash::crc32(res.name.operator CString()));
                    unsigned char color[4];
                    color[0] = (unsigned char)(c.x * 255);
                    color[1] = (unsigned char)(c.y * 255);
                    color[2] = (unsigned char)(c.z * 255);
                    color[3] = (unsigned char)255;
                    color_dst[i] = *(uint32*)(color);
                }
            }
            break;
        case ColorMapping::RES_INDEX:
            for (int64 i = 0; i < color_dst.count; i++) {
                if (i < mol.atom_residue_indices.count) {
                    vec3 c = compute_color(mol.atom_residue_indices[i]);
                    unsigned char color[4];
                    color[0] = (unsigned char)(c.x * 255);
                    color[1] = (unsigned char)(c.y * 255);
                    color[2] = (unsigned char)(c.z * 255);
                    color[3] = (unsigned char)(255);
                    color_dst[i] = *(uint32*)(color);
                }
            }
            break;
        case ColorMapping::CHAIN_ID:
            for (int64 i = 0; i < color_dst.count; i++) {
                if (i < mol.atom_residue_indices.count) {
                    const auto& res = mol.residues[mol.atom_residue_indices[i]];
                    if (res.chain_idx < mol.chains.count) {
                        vec3 c = compute_color(hash::crc32(res.name.operator CString()));
                        unsigned char color[4];
                        color[0] = (unsigned char)(c.x * 255);
                        color[1] = (unsigned char)(c.y * 255);
                        color[2] = (unsigned char)(c.z * 255);
                        color[3] = (unsigned char)(255);
                        color_dst[i] = *(uint32*)(color);
                    }
                }
            }
        case ColorMapping::CHAIN_INDEX:
            for (int64 i = 0; i < color_dst.count; i++) {
                if (i < mol.atom_residue_indices.count) {
                    const auto& res = mol.residues[mol.atom_residue_indices[i]];
                    if (res.chain_idx < mol.chains.count) {
                        vec3 c = compute_color(res.chain_idx);
                        unsigned char color[4];
                        color[0] = (unsigned char)(c.x * 255);
                        color[1] = (unsigned char)(c.y * 255);
                        color[2] = (unsigned char)(c.z * 255);
                        color[3] = (unsigned char)(255);
                        color_dst[i] = *(uint32*)(color);
                    }
                }
            }
        default:
            break;
    }
}

namespace hydrogen_bond {

// Computes the potential donors given a set of atom labels.
// OH and NH atoms are assumed to be donors if the concecutive atom is marked with 'H' for Hydrogen.
int32 compute_donors(DynamicArray<HydrogenBondDonor>* donors, Array<const Label> labels) {
    ASSERT(donors);
    int32 pre_count = (int32)donors->count;
    const int32 num_labels = (int32)labels.count;
    for (int32 i = 0; i < num_labels; i++) {
        if (compare_n(labels[i], "OH", 2) || compare_n(labels[i], "NH", 2)) {
            if (i + 1 < num_labels && compare_n(labels[i + 1], "H", 1)) {
                donors->push_back({i, i + 1});
            }
            if (i + 2 < num_labels && compare_n(labels[i + 2], "H", 1)) {
                donors->push_back({i, i + 2});
            }
        }
    }
    return (int32)donors->count - pre_count;
}

DynamicArray<HydrogenBondDonor> compute_donors(Array<const Label> labels) {
    DynamicArray<HydrogenBondDonor> donors;
    compute_donors(&donors, labels);
    return donors;
}

// Computes the potential acceptors given a set of atom elements.
// This essentially just a filter on atom element which extracts Oxygen and Nitrogen
int32 compute_acceptors(DynamicArray<HydrogenBondAcceptor>* acceptors, Array<const Element> elements) {
    ASSERT(acceptors);
    const int32 pre_count = (int32)acceptors->count;
    for (int32 i = 0; i < (int32)elements.count; i++) {
        if (elements[i] == Element::O || elements[i] == Element::N) {
            acceptors->push_back(i);
        }
    }
    return (int32)acceptors->count - pre_count;
}

DynamicArray<HydrogenBondAcceptor> compute_acceptors(Array<const Element> elements) {
    DynamicArray<HydrogenBondAcceptor> acceptors;
    compute_acceptors(&acceptors, elements);
    return acceptors;
}

// Computes hydrogen bonds given a certain set of potential donors, acceptors and atomic positions from a frame.
// The distance cutoff sets the distance from bonds to potential acceptors.
//

int32 compute_bonds(DynamicArray<HydrogenBond>* bonds, Array<const HydrogenBondDonor> donors, Array<const HydrogenBondAcceptor> acceptors,
                    Array<const vec3> atom_positions, float dist_cutoff, float angle_cutoff) {
    const int32 num_acceptors = (int32)acceptors.count;
    if (!num_acceptors) return 0;

    DynamicArray<vec3> acceptor_pos(num_acceptors);
    DynamicArray<AtomIdx> acceptor_idx(num_acceptors);
    for (int32 i = 0; i < num_acceptors; i++) {
        acceptor_pos[i] = atom_positions[acceptors[i]];
        acceptor_idx[i] = acceptors[i];
    }

    int32 pre_count = (int32)bonds->count;
    spatialhash::Frame frame = spatialhash::compute_frame(acceptor_pos, vec3(dist_cutoff));
    for (const auto& don : donors) {
        vec3 donor_pos = atom_positions[don.donor_idx];
        vec3 hydro_pos = atom_positions[don.hydro_idx];
        spatialhash::for_each_within(frame, hydro_pos, dist_cutoff,
                                     [bonds, &donor_pos, &hydro_pos, &acceptor_idx, &don, angle_cutoff](int32 idx, const vec3& pos) {
                                         AtomIdx g_idx = acceptor_idx[idx];
                                         if (g_idx == don.donor_idx) return;
                                         const vec3 a = hydro_pos - donor_pos;
                                         const vec3 b = pos - hydro_pos;
                                         if (math::angle(a, b) < angle_cutoff) {
                                             bonds->push_back({g_idx, don.donor_idx, don.hydro_idx});
                                         }
                                     });
    }
    return (int32)bonds->count - pre_count;
}

DynamicArray<HydrogenBond> compute_bonds(Array<const HydrogenBondDonor> donors, Array<const HydrogenBondAcceptor> acceptors,
                                         Array<const vec3> atom_positions, float dist_cutoff, float angle_cutoff) {
    DynamicArray<HydrogenBond> bonds;
    compute_bonds(&bonds, donors, acceptors, atom_positions, dist_cutoff, angle_cutoff);
    return bonds;
}

void compute_bonds_trajectory(HydrogenBondTrajectory* hbt, const MoleculeDynamic& dyn, float max_dist, float max_angle) {
    ASSERT(hbt);
    for (int32 i = 0; i < dyn.trajectory.num_frames; i++) {
        Array<HydrogenBond> frame_bonds{(HydrogenBond*)(hbt->bond_data.end()), int64(0)};
        Array<const vec3> atom_positions = get_trajectory_positions(dyn.trajectory, i);
        frame_bonds.count = compute_bonds(&hbt->bond_data, dyn.molecule.hydrogen_bond.donors, dyn.molecule.hydrogen_bond.acceptors, atom_positions,
                                          max_dist, max_angle);
    }
}

HydrogenBondTrajectory compute_bonds_trajectory(const MoleculeDynamic& dyn, float max_dist, float max_angle) {
    HydrogenBondTrajectory hbt;
    compute_bonds_trajectory(&hbt, dyn, max_dist, max_angle);
    return hbt;
}
}  // namespace hydrogen_bond

namespace filter {
static DynamicArray<FilterCommand> filter_commands;

static bool is_modifier(CString str) {
    if (compare(str, "and", true)) return true;
    if (compare(str, "or", true)) return true;
    return false;
}

static bool is_keyword(CString str) {
    if (compare(str, "and", true)) return true;
    if (compare(str, "or", true)) return true;
    if (compare(str, "not", true)) return true;
    return false;
}

int32 count_parentheses(CString str) {
    int beg_parentheses_count = 0;
    int end_parentheses_count = 0;

    for (int64 i = 0; i < str.count; i++) {
        if (str[i] == '(') beg_parentheses_count++;
        if (str[i] == ')') end_parentheses_count++;
    }
    return beg_parentheses_count - end_parentheses_count;
}

CString extract_parenthesis(CString str) {
    const char* beg = str.beg();
    while (beg != str.end() && *beg != '(') beg++;
    if (beg == str.end()) return {};

    const char* end = beg + 1;
    int count = 1;
    while (end++ != str.end() && count > 0) {
        if (*end == '(') count++;
        if (*end == ')') count--;
    }
    if (end == str.end()) return {};
    return CString(beg, end);
}

DynamicArray<CString> extract_chunks(CString str) {
    DynamicArray<CString> chunks;

    const char* beg = str.beg();
    while (beg != str.end()) {
        if (*beg == '(') {
            CString par = extract_parenthesis(CString(beg, str.end()));
            chunks.push_back({par.beg(), par.end()});  // Exclude actual parentheses
            beg = par.end();
        } else if (*beg != ' ') {
            const char* end = beg;
            while (end != str.end() && *end != ' ') end++;
            chunks.push_back(CString(beg, end));
            beg = end;
        } else
            beg++;
    }

    DynamicArray<CString> big_chunks;

    CString* chunk = chunks.beg();
    while (chunk != chunks.end()) {
        if (chunk->front() == '(') {
            big_chunks.push_back(*chunk);
            chunk++;
        } else if (is_keyword(*chunk)) {
            big_chunks.push_back(*chunk);
            chunk++;
        } else {
            CString* beg_chunk = chunk;
            CString* end_chunk = chunk + 1;
            while (end_chunk != chunks.end() && !is_modifier(*end_chunk) && end_chunk->front() != '(') end_chunk++;
            big_chunks.push_back({beg_chunk->beg(), (end_chunk - 1)->end()});
            chunk = end_chunk;
        }
    }

    return big_chunks;
}

FilterCommand* find_filter_command(CString command) {
    for (auto& f : filter_commands) {
        if (compare(command, f.keyword)) return &f;
    }
    return nullptr;
}

void combine_mask_and(Array<bool> dst, Array<bool> src_a, Array<bool> src_b, bool state_not) {
    if (state_not) {
        for (int i = 0; i < dst.count; i++) dst[i] = src_a[i] & !src_b[i];
    } else {
        for (int i = 0; i < dst.count; i++) dst[i] = src_a[i] & src_b[i];
    }
}

void combine_mask_or(Array<bool> dst, Array<bool> src_a, Array<bool> src_b, bool state_not) {
    if (state_not) {
        for (int i = 0; i < dst.count; i++) dst[i] = src_a[i] | !src_b[i];
    } else {
        for (int i = 0; i < dst.count; i++) dst[i] = src_a[i] | src_b[i];
    }
}

bool internal_filter_mask(Array<bool> mask, const MoleculeDynamic& dyn, CString filter) {
    DynamicArray<CString> chunks = extract_chunks(filter);
    DynamicArray<bool> chunk_mask(mask.count);

    bool state_and = true;
    bool state_or = false;
    bool state_not = false;

    for (const auto& chunk : chunks) {
        if (compare(chunk, "and", true)) {
            state_and = true;
        } else if (compare(chunk, "or", true)) {
            state_or = true;
        } else if (compare(chunk, "not", true)) {
            state_not = true;
        } else {
            if (chunk.front() == '(') {
                ASSERT(chunk.back() == ')');
                if (!internal_filter_mask(chunk_mask, dyn, CString(chunk.beg() + 1, chunk.end() - 1))) return false;
            } else {
                auto tokens = ctokenize(chunk);
                auto cmd = find_filter_command(tokens[0]);
                if (!cmd) {
                    StringBuffer<32> buf = tokens[0];
                    LOG_ERROR("Could not match command: '%s'\n", buf.beg());
                    return false;
                }

                auto args = tokens.sub_array(1);

                while (args.count > 0 && compare(args[0], "not", true)) {
                    state_not = !state_not;
                    args = args.sub_array(1);
                }

                if (!cmd->func(chunk_mask, dyn, args)) {
                    StringBuffer<32> buf = tokens[0];
                    LOG_ERROR("Could not parse command: '%s' with arguments: ", buf.beg());
                    for (const auto& arg : args) {
                        buf = arg;
                        printf("'%s'", buf.beg());
                    }
                    return false;
                }
            }

            if (state_and)
                combine_mask_and(mask, mask, chunk_mask, state_not);
            else if (state_or)
                combine_mask_or(mask, mask, chunk_mask, state_not);

            state_and = false;
            state_or = false;
            state_not = false;
        }

        if (state_and && state_or) {
            LOG_ERROR("Cannot use both 'and' and 'or' to combine filter options\n");
            return false;
        }
    }

    return true;
}

bool compute_filter_mask(Array<bool> mask, const MoleculeDynamic& dyn, CString filter) {
    auto count = dyn.molecule.atom_elements.count;
    ASSERT(count == dyn.molecule.atom_labels.count);
    ASSERT(count == dyn.molecule.atom_positions.count);
    ASSERT(count == dyn.molecule.atom_residue_indices.count);
    ASSERT(count == mask.count);

    if (count_parentheses(filter) != 0) {
        LOG_ERROR("Unmatched parentheses\n");
        return false;
    }

    memset(mask.data, 1, mask.count);
    return internal_filter_mask(mask, dyn, filter);
}

void filter_colors(Array<uint32> colors, Array<bool> mask) {
    ASSERT(colors.count == mask.count);
    for (int i = 0; i < colors.count; i++) {
        if (mask[i])
            colors[i] |= 0xff000000;
        else
            colors[i] &= ~0xff000000;
    }
}

void initialize() {

    /*
            all
            water
            aminoacid
            backbone?
                        protein

            name
            element
            atomicnumber
            atom
            residue
            resname
            resid
            chain
            chainid
    */

    auto filter_amino_acid = [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString>) {
        memset(mask.data, 0, mask.count);
        for (const auto& res : dyn.molecule.residues) {
            if (is_amino_acid(res)) {
                memset(mask.data + res.beg_atom_idx, 1, (res.end_atom_idx - res.beg_atom_idx));
            }
        }
        return true;
    };

    filter_commands.push_back({"all", [](Array<bool> mask, const MoleculeDynamic&, Array<const CString>) {
                                   memset(mask.data, 1, mask.count);
                                   return true;
                               }});
    filter_commands.push_back({"water", [](Array<bool>, const MoleculeDynamic&, Array<const CString>) { return true; }});  // NOT DONE
    filter_commands.push_back({"aminoacid", filter_amino_acid});
    filter_commands.push_back({"backbone", [](Array<bool>, const MoleculeDynamic&, Array<const CString>) { return true; }});  // NOT DONE
    filter_commands.push_back({"protein", filter_amino_acid});

    filter_commands.push_back({"name", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   if (args.count == 0) return false;

                                   for (int i = 0; i < dyn.molecule.atom_labels.count; i++) {
                                       mask[i] = false;
                                       for (const auto& arg : args) {
                                           if (compare(dyn.molecule.atom_labels[i], arg)) {
                                               mask[i] = true;
                                               break;
                                           }
                                       }
                                   }
                                   return true;
                               }});

    filter_commands.push_back({"label", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   return find_filter_command("label")->func(mask, dyn, args);
                               }});

    filter_commands.push_back({"element", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   Array<Element> elements = {(Element*)(TMP_MALLOC(args.count * sizeof(Element))), args.count};
                                   for (int i = 0; i < elements.count; i++) {
                                       elements[i] = element::get_from_string(args[i]);
                                       if (elements[i] == Element::Unknown) return false;
                                   }

                                   for (int i = 0; i < dyn.molecule.atom_elements.count; i++) {
                                       mask[i] = false;
                                       for (const auto& ele : elements) {
                                           if (dyn.molecule.atom_elements[i] == ele) {
                                               mask[i] = true;
                                               break;
                                           }
                                       }
                                   }
                                   TMP_FREE(elements.data);
                                   return true;
                               }});

    filter_commands.push_back({"atomicnumber", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   DynamicArray<IntRange> ranges;
                                   if (!extract_ranges(&ranges, args)) return false;
                                   for (int i = 0; i < dyn.molecule.atom_elements.count; i++) {
                                       int atomnr = (int)dyn.molecule.atom_elements[i];
                                       mask[i] = false;
                                       for (auto range : ranges) {
                                           if (range.x <= atomnr && atomnr <= range.y) {
                                               mask[i] = true;
                                               break;
                                           }
                                       }
                                   }
                                   return true;
                               }});

    filter_commands.push_back({"atom", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   DynamicArray<IntRange> ranges;
                                   if (!extract_ranges(&ranges, args)) return false;
                                   memset(mask.data, 0, mask.size_in_bytes());
                                   for (auto range : ranges) {
                                       range.x = math::clamp(range.x - 1, 0, (int32)dyn.molecule.atom_positions.count - 1);
                                       range.y = math::clamp(range.y - 1, 0, (int32)dyn.molecule.atom_positions.count - 1);
                                       if (range.x == range.y)
                                           mask[range.x] = true;
                                       else
                                           memset(mask.data + range.x, 1, range.y - range.x + 1);
                                   }
                                   return true;
                               }});

    filter_commands.push_back({"residue", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   DynamicArray<IntRange> ranges;
                                   if (!extract_ranges(&ranges, args)) return false;
                                   memset(mask.data, 0, mask.size_in_bytes());
                                   for (auto range : ranges) {
                                       range.x = math::clamp(range.x, 0, (int32)dyn.molecule.residues.count - 1);
                                       range.y = math::clamp(range.y, 0, (int32)dyn.molecule.residues.count - 1);
                                       for (int i = range.x; i <= range.y; i++) {
                                           int beg = dyn.molecule.residues[i].beg_atom_idx;
                                           int end = dyn.molecule.residues[i].end_atom_idx;
                                           memset(mask.data + beg, 1, end - beg);
                                       }
                                   }
                                   /*
   for (int i = 0; i < args.count; i++) {
       auto res = to_int(args[i]);
       if (!res.success) return false;
       int res_idx = res.value;
       if (res_idx < 0 || dyn.molecule.residues.count <= res_idx) return false;
       int beg = dyn.molecule.residues[res_idx].beg_atom_idx;
       int end = dyn.molecule.residues[res_idx].end_atom_idx;
       memset(mask.data + beg, 1, end - beg);
   }
                                   */
                                   return true;
                               }});

    filter_commands.push_back({"resname", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   memset(mask.data, 0, mask.count);
                                   for (int i = 0; i < args.count; i++) {
                                       for (const auto& res : dyn.molecule.residues) {
                                           if (compare(args[i], res.name)) {
                                               int beg = res.beg_atom_idx;
                                               int end = res.end_atom_idx;
                                               memset(mask.data + beg, 1, end - beg);
                                           }
                                       }
                                   }
                                   return true;
                               }});

    filter_commands.push_back({"resid", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   DynamicArray<IntRange> ranges;
                                   if (!extract_ranges(&ranges, args)) return false;
                                   memset(mask.data, 0, mask.size_in_bytes());
                                   for (auto range : ranges) {
                                       range.x = math::clamp(range.x - 1, 0, (int32)dyn.molecule.residues.count - 1);
                                       range.y = math::clamp(range.y - 1, 0, (int32)dyn.molecule.residues.count - 1);
                                       for (int i = range.x; i <= range.y; i++) {
                                           int beg = dyn.molecule.residues[i].beg_atom_idx;
                                           int end = dyn.molecule.residues[i].end_atom_idx;
                                           memset(mask.data + beg, 1, end - beg);
                                       }
                                   }
                                   return true;
                               }});

    filter_commands.push_back({"chain", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   DynamicArray<IntRange> ranges;
                                   if (!extract_ranges(&ranges, args)) return false;
                                   memset(mask.data, 0, mask.size_in_bytes());
                                   for (auto range : ranges) {
                                       range.x = math::clamp(range.x - 1, 0, (int32)dyn.molecule.residues.count - 1);
                                       range.y = math::clamp(range.y - 1, 0, (int32)dyn.molecule.residues.count - 1);
                                       for (int i = range.x; i <= range.y; i++) {
                                           Chain chain = get_chain(dyn.molecule, (ChainIdx)i);
                                           int beg = get_atom_beg_idx(dyn.molecule, chain);
                                           int end = get_atom_end_idx(dyn.molecule, chain);
                                           memset(mask.data + beg, 1, end - beg);
                                       }
                                   }
                                   return true;
                               }});

    filter_commands.push_back({"chainid", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   memset(mask.data, 0, mask.count);
                                   for (int i = 0; i < args.count; i++) {
                                       for (const auto& chain : dyn.molecule.chains) {
                                           if (compare(args[i], chain.id)) {
                                               int beg = get_atom_beg_idx(dyn.molecule, chain);
                                               int end = get_atom_end_idx(dyn.molecule, chain);
                                               memset(mask.data + beg, 1, end - beg);
                                           }
                                       }
                                   }
                                   return true;
                               }});
}

void shutdown() {}

}  // namespace filter

namespace draw {

static GLuint instanced_quad_vao = 0;
static GLuint instanced_quad_ibo = 0;

static GLuint vbo = 0;
static GLsizeiptr vbo_size = MEGABYTES(4);

inline void draw_instanced_quads(int num_instances) {
    glBindVertexArray(instanced_quad_vao);
    glDrawElementsInstanced(GL_TRIANGLE_STRIP, 4, GL_UNSIGNED_BYTE, 0, num_instances);
    glBindVertexArray(0);
}

inline void set_vbo_data(const void* data, GLsizeiptr size_in_bytes) {
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    if (vbo_size > size_in_bytes) {
        glBufferSubData(GL_ARRAY_BUFFER, 0, size_in_bytes, data);
    } else {
        glBufferData(GL_ARRAY_BUFFER, size_in_bytes, data, GL_STREAM_DRAW);
        vbo_size = size_in_bytes;
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

namespace vdw {
static GLuint v_shader = 0;
static GLuint f_shader = 0;
static GLuint program = 0;

static GLuint buf_position_radius = 0;
static GLuint buf_color = 0;
static GLuint buf_picking = 0;
static GLuint tex_position_radius = 0;
static GLuint tex_color = 0;
static GLuint tex_picking = 0;

static GLint uniform_loc_view_mat = -1;
static GLint uniform_loc_proj_mat = -1;
static GLint uniform_loc_inv_proj_mat = -1;
static GLint uniform_loc_fov = -1;
static GLint uniform_loc_tex_pos_rad = -1;
static GLint uniform_loc_tex_color = -1;
static GLint uniform_loc_tex_picking = -1;

static const char* v_shader_src = R"(
#version 150 core

uniform mat4 u_view_mat;
uniform mat4 u_proj_mat;
uniform mat4 u_inv_proj_mat;

uniform samplerBuffer u_tex_pos_rad;
uniform samplerBuffer u_tex_color;
uniform samplerBuffer u_tex_picking;

out Fragment {
    flat vec4 color;
    flat vec4 view_sphere;
	flat vec4 picking_color;
    smooth vec4 view_coord;
} out_frag;

// From Inigo Quilez!
void proj_sphere(in vec4 sphere, 
				 in float fle,
				 out vec2 axis_a,
				 out vec2 axis_b,
				 out vec2 center) {
	vec3  o = sphere.xyz;
    float r2 = sphere.w*sphere.w;
	float z2 = o.z*o.z;	
	float l2 = dot(o,o);
	
	// axis
	axis_a = fle*sqrt(-r2*(r2-l2)/((l2-z2)*(r2-z2)*(r2-z2)))*vec2( o.x,o.y);
	axis_b = fle*sqrt(-r2*(r2-l2)/((l2-z2)*(r2-z2)*(r2-l2)))*vec2(-o.y,o.x);
	center = -fle*o.z*o.xy/(z2-r2);
}

void main() {
	int VID = gl_VertexID;
	int IID = gl_InstanceID;
	vec2 uv = vec2(VID / 2, VID % 2) * 2.0 - 1.0; 

	vec4 pos_rad = texelFetch(u_tex_pos_rad, IID);
	vec4 color	 = texelFetch(u_tex_color, IID);
	vec4 picking = texelFetch(u_tex_picking, IID);

	vec3 pos = pos_rad.xyz;
	float rad = pos_rad.w;

    vec4 view_coord = u_view_mat * vec4(pos, 1.0);
    float len = length(view_coord.xyz);
    vec3 view_dir = view_coord.xyz / len;

    out_frag.color = color;
    out_frag.view_sphere = vec4(view_coord.xyz, rad);
	out_frag.picking_color = picking;

	// Focal length
	float fle = u_proj_mat[1][1];
	// 1.0 / aspect_ratio
	float inv_ar = u_proj_mat[0][0] / u_proj_mat[1][1];

	vec2 axis_a;
	vec2 axis_b;
	vec2 center;
	proj_sphere(out_frag.view_sphere, fle, axis_a, axis_b, center);

	// Bias depth with radius of sphere to make sure we only write depth values greater than this
	view_coord.z += rad;

	vec2 xy = (center + axis_a * uv.x + axis_b * uv.y) * vec2(inv_ar, 1.0);
	float z = -u_proj_mat[2][2] - u_proj_mat[3][2] / view_coord.z;
	vec4 pc = vec4(xy, z, 1);

	// Compute view_coordinate which is used for ray-sphere intersection in fragment shader
	out_frag.view_coord = u_inv_proj_mat * pc;
	out_frag.view_coord = out_frag.view_coord / out_frag.view_coord.w;

    gl_Position = pc;
}
)";

static const char* f_shader_src = R"(
#version 150 core
#extension GL_ARB_conservative_depth : enable
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_proj_mat;
uniform float u_exposure = 1.0;

in Fragment {
    flat vec4 color;
    flat vec4 view_sphere;
	flat vec4 picking_color;
    smooth vec4 view_coord;
} in_frag;

#ifdef GL_EXT_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif
layout(location = 0) out vec4 out_color;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec4 out_picking_id;

// https://aras-p.info/texts/CompactNormalStorage.html
vec4 encode_normal (vec3 n) {
    float p = sqrt(n.z*8+8);
    return vec4(n.xy/p + 0.5,0,0);
}

void main() {
    vec3 center = in_frag.view_sphere.xyz;
    float radius = in_frag.view_sphere.w;
    vec3 view_dir = -normalize(in_frag.view_coord.xyz);

    vec3 m = -center;
    vec3 d = -view_dir;
    float r = radius;
    float b = dot(m, d);
    float c = dot(m, m) - r*r;
    float discr = b*b - c;
    if (discr < 0.0) discard;
    float t = -b -sqrt(discr);

    vec3 view_hit = d * t;
    vec3 view_normal = (view_hit - center) / radius;

    gl_FragDepth = (-u_proj_mat[2][2] - u_proj_mat[3][2] / view_hit.z) * 0.5 + 0.5;
    out_color = in_frag.color;
	out_normal = encode_normal(view_normal);
	out_picking_id = in_frag.picking_color;
}
)";

static void initialize() {
    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    v_shader = glCreateShader(GL_VERTEX_SHADER);
    f_shader = glCreateShader(GL_FRAGMENT_SHADER);

    glShaderSource(v_shader, 1, &v_shader_src, 0);
    glShaderSource(f_shader, 1, &f_shader_src, 0);

    glCompileShader(v_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, v_shader)) {
        LOG_ERROR("Compiling vdw vertex shader:\n%s\n", buffer);
    }

    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        LOG_ERROR("Compiling vdw fragment shader:\n%s\n", buffer);
    }

    program = glCreateProgram();
    glAttachShader(program, v_shader);
    glAttachShader(program, f_shader);
    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, program)) {
        LOG_ERROR("Linking vdw program:\n%s\n", buffer);
    }

    glDetachShader(program, v_shader);
    glDetachShader(program, f_shader);

    glDeleteShader(v_shader);
    glDeleteShader(f_shader);

    if (!buf_position_radius) {
        glGenBuffers(1, &buf_position_radius);
        glBindBuffer(GL_ARRAY_BUFFER, buf_position_radius);
        glBufferData(GL_ARRAY_BUFFER, MEGABYTES(20), 0, GL_STREAM_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    if (!buf_color) {
        glGenBuffers(1, &buf_color);
        glBindBuffer(GL_ARRAY_BUFFER, buf_color);
        glBufferData(GL_ARRAY_BUFFER, MEGABYTES(5), 0, GL_STREAM_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    if (!buf_picking) {
        glGenBuffers(1, &buf_picking);
        glBindBuffer(GL_ARRAY_BUFFER, buf_picking);
        glBufferData(GL_ARRAY_BUFFER, MEGABYTES(5), 0, GL_STREAM_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    if (!tex_position_radius) glGenTextures(1, &tex_position_radius);
    if (!tex_color) glGenTextures(1, &tex_color);
    if (!tex_picking) glGenTextures(1, &tex_picking);

    uniform_loc_view_mat = glGetUniformLocation(program, "u_view_mat");
    uniform_loc_proj_mat = glGetUniformLocation(program, "u_proj_mat");
    uniform_loc_inv_proj_mat = glGetUniformLocation(program, "u_inv_proj_mat");
    uniform_loc_fov = glGetUniformLocation(program, "u_fov");
    uniform_loc_tex_pos_rad = glGetUniformLocation(vdw::program, "u_tex_pos_rad");
    uniform_loc_tex_color = glGetUniformLocation(vdw::program, "u_tex_color");
    uniform_loc_tex_picking = glGetUniformLocation(vdw::program, "u_tex_picking");
}

static void shutdown() {
    if (program) glDeleteProgram(program);

    if (buf_position_radius) glDeleteBuffers(1, &buf_position_radius);
    if (buf_color) glDeleteBuffers(1, &buf_color);
    if (!buf_picking) glDeleteBuffers(1, &buf_picking);

    if (tex_position_radius) glDeleteTextures(1, &tex_position_radius);
    if (tex_color) glDeleteTextures(1, &tex_color);
    if (buf_picking) glDeleteTextures(1, &buf_picking);
}

}  // namespace vdw

namespace licorice {
struct Vertex {
    vec3 position;
    uint32 color;
};

static GLuint vao = 0;
static GLuint ibo = 0;

static GLuint v_shader = 0;
static GLuint g_shader = 0;
static GLuint f_shader = 0;
static GLuint program = 0;

static GLint attrib_loc_pos = -1;
static GLint attrib_loc_col = -1;
static GLint uniform_loc_view_mat = -1;
static GLint uniform_loc_proj_mat = -1;
static GLint uniform_loc_radius = -1;

static const char* v_shader_src = R"(
#version 150 core

uniform mat4 u_view_mat;

in vec3	 v_position;
in vec4  v_color;

out Vertex {
    flat vec4 color;
	flat uint picking_id;
} out_vert;

void main() {
    gl_Position = u_view_mat * vec4(v_position, 1.0);
    out_vert.color = v_color;
	out_vert.picking_id = uint(gl_VertexID);
}
)";

static const char* g_shader_src = R"(
#version 150 core

uniform mat4 u_proj_mat;
uniform float u_radius = 1.0;

layout (lines) in;
layout (triangle_strip, max_vertices = 24) out;

in Vertex {
    flat vec4 color;
	flat uint picking_id;
} in_vert[];

out Fragment {
    flat vec4 color[2];
	flat vec4 picking_color[2];
    smooth vec3 view_pos;

	flat vec4  capsule_center_radius;
	flat vec4  capsule_axis_length;
} out_frag;

vec4 pack_u32(uint data) {
	return vec4(
        (data & uint(0x000000FF)) >> 0,
        (data & uint(0x0000FF00)) >> 8,
        (data & uint(0x00FF0000)) >> 16,
        (data & uint(0xFF000000)) >> 24) / 255.0;
}

vec4 prismoid[8];

void emit_vertex(int a){
    out_frag.view_pos = prismoid[a].xyz;
	gl_Position = u_proj_mat * prismoid[a];
    EmitVertex();
}

void emit(int a, int b, int c, int d)
{
    emit_vertex(a);
    emit_vertex(b);
    emit_vertex(c);
    emit_vertex(d);
    EndPrimitive(); 
}

vec3 get_ortho_vec(vec3 v, vec3 A, vec3 B){
    if(abs(1-dot(v,A))>0.001){
        return normalize(cross(v,A));
    }else{
        return normalize(cross(v,B));
    }
}

void main()
{
    if (in_vert[0].color.a == 0 || in_vert[1].color.a == 0) {
        EndPrimitive();
        return;
    }

    // Compute orientation vectors for the two connecting faces:
    vec3 p0 = gl_in[0].gl_Position.xyz;
    vec3 p1 = gl_in[1].gl_Position.xyz;
	float r = u_radius;
	float l = distance(p0, p1);
	vec3 a = (p1 - p0) / l;
	vec3 c = (p0 + p1) * 0.5;

	out_frag.color[0] = in_vert[0].color;
	out_frag.color[1] = in_vert[1].color;

    out_frag.picking_color[0] = pack_u32(in_vert[0].picking_id);
    out_frag.picking_color[1] = pack_u32(in_vert[1].picking_id);

	out_frag.capsule_center_radius = vec4(c, r);
	out_frag.capsule_axis_length = vec4(a, l);

    // Extend end points to properly fit the sphere caps
    p0 -= a * r;
    p1 += a * r;

	vec3 B = vec3(0,0,1);
	vec3 A = vec3(1,0,0);
    vec3 o = get_ortho_vec(a,A,B);

    // Declare scratch variables for basis vectors:
    vec3 i,j,k;

    // Compute vertices of prismoid:
    j = a; i = o; k = normalize(cross(i, j)); i = normalize(cross(k, j)); ; i *= r; k *= r;
    prismoid[0] = vec4(p0 + i + k, 1);
    prismoid[1] = vec4(p0 + i - k, 1);
    prismoid[2] = vec4(p0 - i - k, 1);
    prismoid[3] = vec4(p0 - i + k, 1);
    prismoid[4] = vec4(p1 + i + k, 1);
    prismoid[5] = vec4(p1 + i - k, 1);
    prismoid[6] = vec4(p1 - i - k, 1);
    prismoid[7] = vec4(p1 - i + k, 1);

    // Emit the six faces of the prismoid:
    emit(0,1,3,2); emit(5,4,6,7);
    emit(4,5,0,1); emit(3,2,7,6);
    emit(0,3,4,7); emit(2,1,6,5);
}
)";

static const char* f_shader_src = R"(
#version 150 core
#extension GL_ARB_conservative_depth : enable
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_proj_mat;
uniform float u_exposure = 1.0;
uniform float u_radius = 1.0;

in Fragment {
    flat vec4 color[2];
	flat vec4 picking_color[2];
    smooth vec3 view_pos;

	flat vec4  capsule_center_radius;
	flat vec4  capsule_axis_length;
} in_frag;

#ifdef GL_EXT_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif
layout(location = 0) out vec4 out_color;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec4 out_picking_id;

// Source from Ingo Quilez (https://www.shadertoy.com/view/Xt3SzX)
float intersect_capsule(in vec3 ro, in vec3 rd, in vec3 cc, in vec3 ca, float cr,
                      float ch, out vec3 normal, out int side)  // cc center, ca orientation axis, cr radius, ch height
{
    vec3 oc = ro - cc;
    ch *= 0.5;

    float card = dot(ca, rd);
    float caoc = dot(ca, oc);

    float a = 1.0 - card * card;
    float b = dot(oc, rd) - caoc * card;
    float c = dot(oc, oc) - caoc * caoc - cr * cr;
    float h = b * b - a * c;
    if (h < 0.0) return -1.0;
    float t = (-b - sqrt(h)) / a;

    float y = caoc + t * card;
    side = int(y > 0);

    // body
    if (abs(y) < ch) {
        normal = normalize(oc + t * rd - ca * y);
        return t;
    }

    // caps
    float sy = sign(y);
    oc = ro - (cc + sy * ca * ch);
    b = dot(rd, oc);
    c = dot(oc, oc) - cr * cr;
    h = b * b - c;
    if (h > 0.0) {
        t = -b - sqrt(h);
        normal = normalize(oc + rd * t);
        return t;
    }

    return -1.0;
}

vec4 encode_normal (vec3 n) {
    float p = sqrt(n.z*8+8);
    return vec4(n.xy/p + 0.5,0,0);
}

void main() {
    vec3 ro = vec3(0);
    vec3 rd = normalize(in_frag.view_pos);
	vec3 cc = in_frag.capsule_center_radius.xyz;
	float cr = in_frag.capsule_center_radius.w;
	vec3 ca = in_frag.capsule_axis_length.xyz;
    float ch = in_frag.capsule_axis_length.w;

    vec3 normal;
    int side;
    float t = intersect_capsule(ro, rd, cc, ca, cr, ch, normal, side);
    if (t < 0.0) {
        discard;
        return;
    }

    vec3 pos = ro + t*rd;
    vec4 color = in_frag.color[side];
	vec4 picking_color = in_frag.picking_color[side];

    vec4 coord = vec4(0, 0, pos.z, 1);
    coord = u_proj_mat * coord;
    coord = coord / coord.w;

    gl_FragDepth = coord.z * 0.5 + 0.5;
    out_color = color;
	out_normal = encode_normal(normal);
	out_picking_id = picking_color;
}
)";

static void initialize() {
    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    v_shader = glCreateShader(GL_VERTEX_SHADER);
    g_shader = glCreateShader(GL_GEOMETRY_SHADER);
    f_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(v_shader, 1, &v_shader_src, 0);
    glShaderSource(g_shader, 1, &g_shader_src, 0);
    glShaderSource(f_shader, 1, &f_shader_src, 0);

    glCompileShader(v_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, v_shader)) {
        LOG_ERROR("Compiling licorice vertex shader:\n%s\n", buffer);
    }
    glCompileShader(g_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, g_shader)) {
        LOG_ERROR("Compiling licorice geometry shader:\n%s\n", buffer);
    }
    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        LOG_ERROR("Compiling licorice fragment shader:\n%s\n", buffer);
    }

    program = glCreateProgram();
    glAttachShader(program, v_shader);
    glAttachShader(program, g_shader);
    glAttachShader(program, f_shader);
    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, program)) {
        LOG_ERROR("Linking licorice program:\n%s\n", buffer);
    }

    glDetachShader(program, v_shader);
    glDetachShader(program, g_shader);
    glDetachShader(program, f_shader);

    glDeleteShader(v_shader);
    glDeleteShader(g_shader);
    glDeleteShader(f_shader);

    attrib_loc_pos = glGetAttribLocation(program, "v_position");
    attrib_loc_col = glGetAttribLocation(program, "v_color");
    uniform_loc_view_mat = glGetUniformLocation(program, "u_view_mat");
    uniform_loc_proj_mat = glGetUniformLocation(program, "u_proj_mat");
    uniform_loc_radius = glGetUniformLocation(program, "u_radius");

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    glEnableVertexAttribArray(attrib_loc_pos);
    glVertexAttribPointer(attrib_loc_pos, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)0);

    glEnableVertexAttribArray(attrib_loc_col);
    glVertexAttribPointer(attrib_loc_col, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (const GLvoid*)12);

    glBindVertexArray(0);

    glGenBuffers(1, &ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, MEGABYTES(4), nullptr, GL_STREAM_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

static void shutdown() {
    if (vao) glDeleteVertexArrays(1, &vao);
    if (program) glDeleteProgram(program);
    if (ibo) glDeleteBuffers(1, &ibo);
}

}  // namespace licorice

namespace ribbons {

struct Vertex {
    vec3 position;
    vec3 tangent;
    vec3 normal;
    unsigned int color = 0xffffffff;
    unsigned int picking_id = 0xffffffff;
};

static GLuint vao = 0;
static GLuint program = 0;

static GLint attrib_loc_position = -1;
static GLint attrib_loc_tangent = -1;
static GLint attrib_loc_normal = -1;
static GLint attrib_loc_color = -1;
static GLint attrib_loc_picking = -1;

static GLint uniform_loc_view_mat = -1;
static GLint uniform_loc_proj_mat = -1;
static GLint uniform_loc_scale_x = -1;
static GLint uniform_loc_scale_y = -1;

static const char* v_shader_src = R"(
#version 150 core
in vec3 v_position;
in vec3 v_tangent;
in vec3 v_normal;
in vec4 v_color;
in uint v_picking_id;

out Vertex {
    vec3 tangent;
    vec3 normal;
    vec4 color;
    vec4 picking_color;
} out_vert;

vec4 pack_u32(uint data) {
    return vec4(
        (data & uint(0x000000FF)) >> 0,
        (data & uint(0x0000FF00)) >> 8,
        (data & uint(0x00FF0000)) >> 16,
        (data & uint(0xFF000000)) >> 24) / 255.0;
}
 
void main() {
    out_vert.tangent = v_tangent;
    out_vert.normal = v_normal;
    out_vert.color = v_color;
    out_vert.picking_color = pack_u32(v_picking_id);
    gl_Position = vec4(v_position, 1);
} 
)";

static const char* g_shader_src = R"(
#version 150 core

#define CIRCLE_RES 8

layout(lines) in;
layout(triangle_strip, max_vertices = 56) out;

uniform mat4 u_view_mat;
uniform mat4 u_proj_mat;
uniform float u_scale_x = 1.0;
uniform float u_scale_y = 0.1;

in Vertex {
    vec3 tangent;
    vec3 normal;
    vec4 color;
    vec4 picking_color;
} in_vert[];

out Fragment {
    smooth vec4 color;
    smooth vec3 view_normal;
    flat vec4 picking_color;
} out_frag;

void emit(mat4 mat, int input_idx, vec4 v, vec3 n) {
    vec4 view_coord = u_view_mat * mat * v;
    mat3 norm_mat = inverse(transpose(mat3(u_view_mat) * mat3(mat)));
    out_frag.view_normal = normalize(norm_mat * n);
    out_frag.color = in_vert[input_idx].color;
    out_frag.picking_color = in_vert[input_idx].picking_color;
    gl_Position = u_proj_mat * view_coord;
    EmitVertex();
}

mat4 compute_mat(vec3 tangent, vec3 normal, vec3 translation) {
    vec3 binormal = normalize(cross(tangent, normal));
    return mat4(vec4(binormal, 0), vec4(normal, 0), vec4(tangent, 0), vec4(translation, 1));
}

void main() {
    if (in_vert[0].color.a == 0) {
        EndPrimitive();
        return;
    }

    vec3 pos[2];
    pos[0] = gl_in[0].gl_Position.xyz;
    pos[1] = gl_in[1].gl_Position.xyz;

    vec3 t[2];
    t[0] = normalize(in_vert[0].tangent);
    t[1] = normalize(in_vert[1].tangent);

    vec3 n[2];
    n[0] = normalize(in_vert[0].normal);
    n[1] = normalize(in_vert[1].normal);

    mat4 mat[2];
    mat[0] = compute_mat(t[0], n[0], pos[0]);
    mat[1] = compute_mat(t[1], n[1], pos[1]);

	// BOTTOM
	emit(mat[0], 0, vec4( 1 * u_scale_x, -1 * u_scale_y, 0, 1), vec3(0,-1,0));
	emit(mat[0], 0, vec4(-1 * u_scale_x, -1 * u_scale_y, 0, 1), vec3(0,-1,0));
	emit(mat[1], 1, vec4( 1 * u_scale_x, -1 * u_scale_y, 0, 1), vec3(0,-1,0));
	emit(mat[1], 1, vec4(-1 * u_scale_x, -1 * u_scale_y, 0, 1), vec3(0,-1,0));
	EndPrimitive();

	// TOP
	emit(mat[0], 0, vec4(-1 * u_scale_x, 1 * u_scale_y, 0, 1), vec3(0,1,0));
	emit(mat[0], 0, vec4( 1 * u_scale_x, 1 * u_scale_y, 0, 1), vec3(0,1,0));
	emit(mat[1], 1, vec4(-1 * u_scale_x, 1 * u_scale_y, 0, 1), vec3(0,1,0));
	emit(mat[1], 1, vec4( 1 * u_scale_x, 1 * u_scale_y, 0, 1), vec3(0,1,0));
	EndPrimitive();

	// LEFT
	emit(mat[0], 0, vec4(-1 * u_scale_x, -1 * u_scale_y, 0, 1), vec3(-1,0,0));
	emit(mat[0], 0, vec4(-1 * u_scale_x,  1 * u_scale_y, 0, 1), vec3(-1,0,0));
	emit(mat[1], 1, vec4(-1 * u_scale_x, -1 * u_scale_y, 0, 1), vec3(-1,0,0));
	emit(mat[1], 1, vec4(-1 * u_scale_x,  1 * u_scale_y, 0, 1), vec3(-1,0,0));
	EndPrimitive();

	// RIGHT
	emit(mat[0], 0, vec4(1 * u_scale_x,  1 * u_scale_y, 0, 1), vec3(1,0,0));
	emit(mat[0], 0, vec4(1 * u_scale_x, -1 * u_scale_y, 0, 1), vec3(1,0,0));
	emit(mat[1], 1, vec4(1 * u_scale_x,  1 * u_scale_y, 0, 1), vec3(1,0,0));
	emit(mat[1], 1, vec4(1 * u_scale_x, -1 * u_scale_y, 0, 1), vec3(1,0,0));
    EndPrimitive();
}
)";

static const char* f_shader_src = R"(
#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

in Fragment {
    smooth vec4 color;
    smooth vec3 view_normal;
	flat vec4 picking_color;
} in_frag;

layout(location = 0) out vec4 out_color;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec4 out_picking_id;

vec4 encode_normal (vec3 n) {
    float p = sqrt(n.z*8+8);
    return vec4(n.xy/p + 0.5,0,0);
}

void main() {
    out_color = in_frag.color;
	out_normal = encode_normal(in_frag.view_normal);
    out_picking_id = in_frag.picking_color;
}
)";

void intitialize() {
    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    GLuint v_shader = glCreateShader(GL_VERTEX_SHADER);
    GLuint g_shader = glCreateShader(GL_GEOMETRY_SHADER);
    GLuint f_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(v_shader, 1, &v_shader_src, 0);
    glShaderSource(g_shader, 1, &g_shader_src, 0);
    glShaderSource(f_shader, 1, &f_shader_src, 0);

    glCompileShader(v_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, v_shader)) {
        LOG_ERROR("Compiling ribbons vertex shader:\n%s\n", buffer);
    }
    glCompileShader(g_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, g_shader)) {
        LOG_ERROR("Compiling ribbons geometry shader:\n%s\n", buffer);
    }
    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        LOG_ERROR("Compiling ribbons fragment shader:\n%s\n", buffer);
    }

    program = glCreateProgram();
    glAttachShader(program, v_shader);
    glAttachShader(program, g_shader);
    glAttachShader(program, f_shader);
    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, program)) {
        LOG_ERROR("Linking ribbons program:\n%s\n", buffer);
    }

    glDetachShader(program, v_shader);
    glDetachShader(program, g_shader);
    glDetachShader(program, f_shader);

    glDeleteShader(v_shader);
    glDeleteShader(g_shader);
    glDeleteShader(f_shader);

    attrib_loc_position = glGetAttribLocation(program, "v_position");
    attrib_loc_tangent = glGetAttribLocation(program, "v_tangent");
    attrib_loc_normal = glGetAttribLocation(program, "v_normal");
    attrib_loc_color = glGetAttribLocation(program, "v_color");
    attrib_loc_picking = glGetAttribLocation(program, "v_picking_id");

    uniform_loc_view_mat = glGetUniformLocation(program, "u_view_mat");
    uniform_loc_proj_mat = glGetUniformLocation(program, "u_proj_mat");
    uniform_loc_scale_x = glGetUniformLocation(program, "u_scale_x");
    uniform_loc_scale_y = glGetUniformLocation(program, "u_scale_y");

    if (!vao) {
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);

        glEnableVertexAttribArray(attrib_loc_position);
        glVertexAttribPointer(attrib_loc_position, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)0);

        glEnableVertexAttribArray(attrib_loc_tangent);
        glVertexAttribPointer(attrib_loc_tangent, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)12);

        glEnableVertexAttribArray(attrib_loc_normal);
        glVertexAttribPointer(attrib_loc_normal, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)24);

        glEnableVertexAttribArray(attrib_loc_color);
        glVertexAttribPointer(attrib_loc_color, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (const GLvoid*)36);

        glEnableVertexAttribArray(attrib_loc_picking);
        glVertexAttribIPointer(attrib_loc_picking, 1, GL_UNSIGNED_INT, sizeof(Vertex), (const GLvoid*)40);

        glBindVertexArray(0);
    }
}

void shutdown() {
    if (vao) glDeleteVertexArrays(1, &vao);
}

}  // namespace ribbons

void initialize() {
    if (!instanced_quad_ibo) {
        const unsigned char data[4] = {0, 1, 2, 3};
        glGenBuffers(1, &instanced_quad_ibo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, instanced_quad_ibo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 4, data, GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }

    if (!instanced_quad_vao) {
        glGenVertexArrays(1, &instanced_quad_vao);
        glBindVertexArray(instanced_quad_vao);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, instanced_quad_ibo);
        glBindVertexArray(0);
    }

    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, vbo_size, nullptr, GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    vdw::initialize();
    licorice::initialize();
    ribbons::intitialize();
    ramachandran::initialize();
}

void shutdown() {
    if (instanced_quad_vao) glDeleteVertexArrays(1, &instanced_quad_vao);
    if (instanced_quad_ibo) glDeleteBuffers(1, &instanced_quad_ibo);
    if (vbo) glDeleteBuffers(1, &vbo);

    vdw::shutdown();
    licorice::shutdown();
    ribbons::shutdown();
    ramachandran::shutdown();
}

void draw_vdw(Array<const vec3> atom_positions, Array<const float> atom_radii, Array<const uint32> atom_colors, const mat4& view_mat,
              const mat4& proj_mat, float radii_scale) {
    uint32 count = (uint32)atom_positions.count;
    ASSERT(count == atom_radii.count && count == atom_colors.count);

    mat4 inv_proj_mat = math::inverse(proj_mat);

    glBindBuffer(GL_ARRAY_BUFFER, vdw::buf_position_radius);
    vec4* gpu_pos_rad = (vec4*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    glBindBuffer(GL_ARRAY_BUFFER, vdw::buf_color);
    uint32* gpu_color = (uint32*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    glBindBuffer(GL_ARRAY_BUFFER, vdw::buf_picking);
    uint32* gpu_picking = (uint32*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

    unsigned int draw_count = 0;
    // DISCARD ANY ZERO RADII OR ZERO COLOR ALPHA ATOMS HERE
    for (uint32 i = 0; i < count; i++) {
        if (atom_radii[i] <= 0.f) continue;
        if ((atom_colors[i] & 0xff000000) == 0) continue;
        gpu_pos_rad[draw_count] = vec4(atom_positions[i], atom_radii[i] * radii_scale);
        gpu_color[draw_count] = atom_colors[i];
        gpu_picking[draw_count] = i;
        draw_count++;
    }

    glUnmapBuffer(GL_ARRAY_BUFFER);
    glBindBuffer(GL_ARRAY_BUFFER, vdw::buf_color);
    glUnmapBuffer(GL_ARRAY_BUFFER);
    glBindBuffer(GL_ARRAY_BUFFER, vdw::buf_position_radius);
    glUnmapBuffer(GL_ARRAY_BUFFER);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // glEnable(GL_DEPTH_TEST);

    glUseProgram(vdw::program);

    // Texture 0
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, vdw::tex_position_radius);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32F, vdw::buf_position_radius);

    // Texture 1
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, vdw::tex_color);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, vdw::buf_color);

    // Texture 2
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_BUFFER, vdw::tex_picking);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, vdw::buf_picking);

    // Uniforms
    glUniform1i(vdw::uniform_loc_tex_pos_rad, 0);
    glUniform1i(vdw::uniform_loc_tex_color, 1);
    glUniform1i(vdw::uniform_loc_tex_picking, 2);
    glUniformMatrix4fv(vdw::uniform_loc_view_mat, 1, GL_FALSE, &view_mat[0][0]);
    glUniformMatrix4fv(vdw::uniform_loc_proj_mat, 1, GL_FALSE, &proj_mat[0][0]);
    glUniformMatrix4fv(vdw::uniform_loc_inv_proj_mat, 1, GL_FALSE, &inv_proj_mat[0][0]);

    // Draw
    draw_instanced_quads(draw_count);

    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // glDisable(GL_DEPTH_TEST);
}

void draw_licorice(Array<const vec3> atom_positions, Array<const Bond> atom_bonds, Array<const uint32> atom_colors, const mat4& view_mat,
                   const mat4& proj_mat, float radii_scale) {
    ASSERT(atom_positions.count == atom_colors.count);

    const auto num_bytes = atom_positions.count * sizeof(licorice::Vertex);
    licorice::Vertex* data = (licorice::Vertex*)TMP_MALLOC(num_bytes);

    for (int64_t i = 0; i < atom_positions.count; i++) {
        data[i].position = atom_positions[i];
        data[i].color = atom_colors[i];
    }

    draw::set_vbo_data(data, atom_positions.count * sizeof(licorice::Vertex));
    TMP_FREE(data);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, licorice::ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, atom_bonds.count * sizeof(Bond), atom_bonds.data, GL_STREAM_DRAW);

    // glEnable(GL_DEPTH_TEST);

    glBindVertexArray(licorice::vao);
    glUseProgram(licorice::program);
    glUniformMatrix4fv(licorice::uniform_loc_view_mat, 1, GL_FALSE, &view_mat[0][0]);
    glUniformMatrix4fv(licorice::uniform_loc_proj_mat, 1, GL_FALSE, &proj_mat[0][0]);
    glUniform1f(licorice::uniform_loc_radius, 0.25f * radii_scale);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, licorice::ibo);
    glDrawElements(GL_LINES, (GLsizei)atom_bonds.count * 2, GL_UNSIGNED_INT, (const void*)0);
    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    // glDisable(GL_DEPTH_TEST);
}

void draw_ribbons(Array<const BackboneSegment> backbone_segments, Array<const Chain> chains, Array<const vec3> atom_positions,
                  Array<const uint32> atom_colors, const mat4& view_mat, const mat4& proj_mat, int num_subdivisions, float tension, float width_scale,
                  float thickness_scale) {
    if (backbone_segments.count == 0) return;
    if (chains.count == 0) return;
    if (atom_positions.count == 0) return;
    if (atom_colors.count == 0) return;

    struct DrawInfo {
        int32 offset;
        int32 count;
    };

    DynamicArray<DrawInfo> draw_data;
    DynamicArray<BackboneSegment> visible_segments;
    DynamicArray<SplineSegment> spline_segments;
    DynamicArray<ribbons::Vertex> vertices;

    int32 offset = 0;
    for (const auto& c : chains) {
        visible_segments.clear();
        for (int i = c.beg_res_idx; i < c.end_res_idx; i++) {
            int ca_idx = backbone_segments[i].ca_idx;
            if (ca_idx == -1 || (atom_colors[ca_idx] & 0xff000000) == 0)
                break;
            else {
                visible_segments.push_back(backbone_segments[i]);
            }
        }
        // Only do this if all segments within a chain was visible
        if (visible_segments.size() == (c.end_res_idx - c.beg_res_idx)) {
            auto splines = compute_spline(atom_positions, atom_colors, visible_segments, num_subdivisions, tension);
            // spline_segments.append(splines);
            for (const auto& s : splines) {
                vertices.push_back({s.position, s.tangent, s.normal, s.color, s.index});
            }
            draw_data.push_back({offset, (int32)splines.count});
            offset += (int32)splines.count;
        }
    }

    set_vbo_data(vertices.data, vertices.size_in_bytes());

    // glEnable(GL_DEPTH_TEST);

    glBindVertexArray(ribbons::vao);
    glUseProgram(ribbons::program);
    glUniformMatrix4fv(ribbons::uniform_loc_view_mat, 1, GL_FALSE, &view_mat[0][0]);
    glUniformMatrix4fv(ribbons::uniform_loc_proj_mat, 1, GL_FALSE, &proj_mat[0][0]);
    glUniform1f(ribbons::uniform_loc_scale_x, 0.5f * width_scale);
    glUniform1f(ribbons::uniform_loc_scale_y, 0.1f * thickness_scale);
    for (const auto& di : draw_data) {
        glDrawArrays(GL_LINE_STRIP, di.offset, di.count);
    }
    glUseProgram(0);
    glBindVertexArray(0);

    // glDisable(GL_DEPTH_TEST);
}

void draw_backbone(Array<const BackboneSegment> backbone_segments, Array<const vec3> atom_positions, const mat4& view_mat, const mat4& proj_mat) {
    immediate::set_view_matrix(view_mat);
    immediate::set_proj_matrix(proj_mat);

    for (const auto& seg : backbone_segments) {
        if (seg.ca_idx > -1) immediate::draw_point(atom_positions[seg.ca_idx], immediate::COLOR_BLACK);
        if (seg.c_idx > -1) immediate::draw_point(atom_positions[seg.c_idx], immediate::COLOR_GREEN);
        if (seg.o_idx > -1) immediate::draw_point(atom_positions[seg.o_idx], immediate::COLOR_RED);
    }

    for (int i = 1; i < backbone_segments.count; i++) {
        if (backbone_segments[i].ca_idx > -1 && backbone_segments[i - 1].ca_idx > -1)
            immediate::draw_line(atom_positions[backbone_segments[i - 1].ca_idx], atom_positions[backbone_segments[i].ca_idx],
                                 immediate::COLOR_WHITE);
    }

    immediate::flush();
}

void draw_spline(Array<const SplineSegment> spline, const mat4& view_mat, const mat4& proj_mat) {
    immediate::set_view_matrix(view_mat);
    immediate::set_proj_matrix(proj_mat);

    for (const auto& seg : spline) {
        immediate::draw_line(seg.position, seg.position + seg.tangent, immediate::COLOR_BLUE);
        immediate::draw_line(seg.position, seg.position + seg.normal, immediate::COLOR_GREEN);
        immediate::draw_line(seg.position, seg.position + seg.binormal, immediate::COLOR_RED);
    }

    for (int64 i = 1; i < spline.count; i++) {
        immediate::draw_line(spline[i - 1].position, spline[i].position, immediate::COLOR_WHITE);
    }

    immediate::flush();
}

}  // namespace draw

namespace ramachandran {

// Segmentation texture data
constexpr int seg_width = 36;
constexpr int seg_height = 36;
constexpr GLenum seg_data_format = GL_BGR;
/* UGLY
constexpr unsigned char seg_data[] =
    R"(  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P                PP PP PP PP PP                 P  P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?
?  ?  ?  ?  P  P                                            P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P  P  P P  P  ?  ?  ?  ?  ?
?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P                                      P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?
P  P  P                                      P  P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P  P P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?
?  ?  ?  ?  ?  ?  ?  ?  P  P  P         P  P  P  P  P  P                   P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  PP  P  P  P
P  P  P  P  P                   P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  PP  P  P  P  P  ?  ?  P  P                   P  P  ?  ?
?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P  P  P  PP  P  ?  ?  ?  ?  P  P  P                   P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P  P
P  P  P  PP  P  ?  ?  ?  ?  P  P  P                   P  P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P            P  P  P  ?  ?  ?  P  P  P  P P
P  P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P            P  P  ?  ?  ?  ?  ?  ?  P  P                P  P  P P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P
P  P             P  P  ?  ?  ?  ?  ?  P  P  P               P  P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P  P             P  P  P  ?  ?  ?  ?  P  P
P               P  P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P  P          P  P  P  ?  ?  ?  ?  ?  P  P               P  P  P  P  ?  ?  ?  ?  ?
?  ?  ?  ?  ?  ?  ?  P  P  P  P          P  P  P  ?  ?  ?  P  P  P               P  P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P  P          P
P  P  P  ?  ?  ?  P  P               P  P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P  P       P  P  P  ?  P  ?  ?  P  P               P  P
P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P  P  P    P  P  P  P  P  P  P  P  P               P  P  P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?
?  ?  ?  P  P  P    P  P  P  P  P  P  P  P  P               P  P  P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P  P P  P  P  P  P  P  ?  ?
?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P  P                                         P  P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P
P  P                                      P  P  P  ?  ?  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P  P P  P  P  P  P  P  ?  ?  ?  ?  ?  ?  ?
?  ?  ?  ?  ?  ?  ?  ?  ?  P  P                                      P  P  P  P  P  P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P  P P  P  P P
P  P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P  P  P  P  P                                                P  P  P P  P  P   P  P  P  P P  ?  P  P  P
P  P  P  P  P                                                   P  P  P  P  P  P  P  ?  ?  P P  P  P  P  P  P  P P  P  P  P  P  ?  ?  P  P  P  P  P  P
P                      PP PP PP PP PP                                P  P  P  P  P  P  ?  ?  ?  ?  P  P  P  P                      PP PP PP PP PP P  ?
?  ?  P  ?  ?  ?  ?  ?  P  P  P  P  P                   PP PP ?? PP PP PP                             P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P  P  P PP
PP ?? PP PP PP                    P  P  P  P  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P                   PP PP ?? ?? PP PP                    P  P  P
?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  ?  P  P                   PP PP PP PP PP PP                    P  P  ?)";
*/
constexpr unsigned char seg_data[] =
    R"(z™Ì UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌz™Ìz™Ìz™Ì­¿Ø­¿Øÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ·ÌÈ·ÌÈ·ÌÈ·ÌÈ·ÌÈÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Ø­¿Ø­¿Øz™Ìz™Ì UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌz™Ìz™Ìz™Ì­¿Ø­¿Øÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Ø­¿Øz™Ìz™Ìz™Ì UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌz™Ìz™Ìz™Ì­¿Ø­¿Ø­¿Ø­¿Øÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Øz™Ìz™Ìz™Ìz™Ìz™Ì UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌz™Ìz™Ìz™Ìz™Ì­¿Ø­¿Øÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Ø­¿Øz™Ìz™Ìz™Ìz™Ì UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌz™Ìz™Ì­¿Ø­¿Ø­¿Øÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Ø­¿Ø­¿Øz™Ìz™Ìz™Ì UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌz™Ìz™Ìz™Ì­¿Ø­¿Ø­¿Øÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Ø­¿Øz™Ìz™Ìz™Ìz™Ì UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌ UÌz™Ìz™Ìz™Ìz™Ì­¿Ø­¿Ø­¿ØÿÿÿÿÿÿÿÿÿÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Øz™Ìz™Ìz™Ìz™Ìz™Ì UÌ UÌ UÌ UÌ UÌ UÌ UÌz™Ìz™Ìz™Ìz™Ìz™Ìz™Ì­¿Ø­¿ØÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Øz™Ìz™Ìz™Ìz™Ìz™Ìz™Ì UÌz™Ìz™Ì UÌz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ì­¿Ø­¿ØÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌžfÌžfÌÂ·ÌÂ·ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Øz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ì­¿Ø­¿Ø­¿Ø­¿Ø­¿ØÌÂ·ÌÂ·ÌžfÌžfÌžfÌžfÌÂ·ÌÂ·ÌÂ·ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Ø­¿Øz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ì­¿Ø­¿Ø­¿Ø­¿Ø­¿Ø­¿Ø­¿ØÌÂ·ÌÂ·ÌžfÌžfÌžfÌžfÌÂ·ÌÂ·ÌÂ·ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Ø­¿Ø­¿Øz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ì­¿Ø­¿ØÿÿÿÿÿÿÿÿÿÿÿÿÌÂ·ÌÂ·ÌÂ·ÌžfÌžfÌžfÌÂ·ÌÂ·ÌÂ·ÌÂ·ÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Ø­¿Ø­¿Ø­¿Øz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ì­¿Ø­¿ØÿÿÿÿÿÿÿÿÿÿÿÿÌÂ·ÌÂ·ÌžfÌžfÌžfÌžfÌžfÌžfÌÂ·ÌÂ·ÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Ø­¿Ø¾Ì·¾Ì·•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz¾Ì·¾Ì·¾Ì·ÿÿÿÿÿÿÿÿÿÿÿÿÌÂ·ÌÂ·ÌžfÌžfÌp ÌžfÌžfÌÂ·ÌÂ·ÌÂ·ÿÿÿÿÿÿÿÿÿÿÿÿ¾Ì·¾Ì·¾Ì·¾Ì·•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz¾Ì·¾Ì·¾Ì·ÿÿÿÿÿÿÿÿÿÿÿÿÌÂ·ÌÂ·ÌÂ·ÌžfÌp Ìp ÌžfÌÂ·ÌÂ·ÌÂ·ÿÿÿÿÿÿÿÿÿÿÿÿ¾Ì·¾Ì·¾Ì·¾Ì·•Ìz•Ìz•Ìz•Ìz•Ìz•ÌzDÌ •Ìz•Ìz•Ìz•Ìz•Ìz¾Ì·¾Ì·¾Ì·ÿÿÿÿÿÿÿÿÿÌÂ·ÌÂ·ÌÂ·ÌžfÌžfÌžfÌžfÌžfÌÂ·ÌÂ·ÿÿÿÿÿÿÿÿÿÿÿÿ¾Ì·¾Ì·¾Ì·¾Ì·•Ìz•Ìz•Ìz•Ìz•ÌzDÌ DÌ DÌ DÌ •Ìz•Ìz•Ìz¾Ì·¾Ì·¾Ì·¾Ì·ÿÿÿÿÿÿÿÿÿÌÂ·ÌÂ·ÌÂ·ÌžfÌžfÌžfÌÂ·ÌÂ·ÌÂ·ÿÿÿÿÿÿÿÿÿÿÿÿ¾Ì·¾Ì·¾Ì·¾Ì·•Ìz•Ìz•Ìz•Ìz•ÌzDÌ DÌ DÌ DÌ DÌ •Ìz•Ìz•Ìz¾Ì·¾Ì·¾Ì·ÿÿÿÿÿÿÿÿÿÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌžfÌžfÌžfÌÂ·ÌÂ·ÿÿÿÿÿÿÿÿÿÿÿÿ¾Ì·¾Ì·¾Ì·¾Ì·•Ìz•Ìz•Ìz•ÌzDÌ DÌ DÌ DÌ DÌ DÌ DÌ •Ìz•Ìz•Ìz¾Ì·¾Ì·¾Ì·ÿÿÿÿÿÿÌÂ·ÌÂ·ÌÂ·ÌžfÌÂ·ÌžfÌžfÌÂ·ÌÂ·ÿÿÿÿÿÿÿÿÿÿÿÿ¾Ì·¾Ì·¾Ì·¾Ì·•Ìz•Ìz•Ìz•Ìz•ÌzDÌ DÌ DÌ DÌ DÌ DÌ DÌ •Ìz•Ìz¾Ì·¾Ì·¾Ì·¾Ì·ÿÿÿÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÿÿÿÿÿÿÿÿÿÿÿÿ¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·•Ìz•Ìz•Ìz•ÌzDÌ DÌ DÌ DÌ DÌ DÌ DÌ DÌ •Ìz•Ìz¾Ì·¾Ì·¾Ì·ÿÿÿÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÌÂ·ÿÿÿÿÿÿÿÿÿÿÿÿ¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·•Ìz•Ìz•Ìz•Ìz•ÌzDÌ DÌ DÌ DÌ DÌ DÌ DÌ •Ìz•Ìz•Ìz¾Ì·¾Ì·¾Ì·ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·•Ìz•Ìz•Ìz•Ìz•ÌzDÌ DÌ DÌ DÌ DÌ DÌ DÌ •Ìz•Ìz¾Ì·¾Ì·¾Ì·ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ¾Ì·¾Ì·¾Ì·¾Ì·•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•ÌzDÌ DÌ DÌ DÌ DÌ DÌ •Ìz•Ìz•Ìz¾Ì·¾Ì·¾Ì·ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ¾Ì·¾Ì·¾Ì·•Ìz•Ìz¾Ì·•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•ÌzDÌ DÌ DÌ DÌ •Ìz•Ìz•Ìz¾Ì·¾Ì·¾Ì·ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz¾Ì·¾Ì·ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz¾Ì·¾Ì·¾Ì·ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ¾Ì·¾Ì·¾Ì·­¿Ø¾Ì·¾Ì·¾Ì·¾Ì·•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz•Ìz¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Ø­¿Ø¾Ì·¾Ì·¾Ì·­¿Ø­¿Ø­¿Ø­¿Ø¾Ì·•Ìz¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Ø­¿Ø­¿Ø­¿Ø­¿Ø­¿Øz™Ìz™Ì­¿Ø¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·¾Ì·ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Ø­¿Ø­¿Ø­¿Øz™Ìz™Ì­¿Ø­¿Ø­¿Ø­¿Ø­¿Ø­¿Ø­¿Øÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ·ÌÈ·ÌÈ·ÌÈ·ÌÈ·ÌÈÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Ø­¿Ø­¿Ø­¿Ø­¿Øz™Ìz™Ìz™Ìz™Ì­¿Ø­¿Ø­¿Ø­¿Øÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ·ÌÈ·ÌÈ·ÌÈ·ÌÈ·ÌÈÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Øz™Ìz™Ìz™Ì­¿Øz™Ìz™Ìz™Ìz™Ìz™Ì­¿Ø­¿Ø­¿Ø­¿Ø­¿Øÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ·ÌÈ·ÌÈQÌ··ÌÈ·ÌÈ·ÌÈÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Øz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ì­¿Ø­¿Ø­¿Ø­¿Øÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ·ÌÈ·ÌÈQÌ··ÌÈ·ÌÈ·ÌÈÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Ø­¿Ø­¿Øz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ì­¿Ø­¿Øÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ·ÌÈ·ÌÈQÌ·QÌ··ÌÈ·ÌÈÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Ø­¿Øz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ìz™Ì­¿Ø­¿Øÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ·ÌÈ·ÌÈ·ÌÈ·ÌÈ·ÌÈ·ÌÈÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ­¿Ø­¿Øz™Ì)";
// Accumulation texture data
constexpr int acc_width = 2048;
constexpr int acc_height = 2048;

static GLuint segmentation_tex = 0;
static GLuint accumulation_tex = 0;
static GLuint coord_tex = 0;
static GLuint coord_buf = 0;
static GLuint fbo = 0;
static GLuint program = 0;

static GLint uniform_loc_coord_tex = -1;
static GLint uniform_loc_instance_offset = -1;
static GLint uniform_loc_radius = -1;
static GLint uniform_loc_color = -1;
static GLint uniform_loc_outline = -1;

GLuint get_accumulation_texture() { return accumulation_tex; }
GLuint get_segmentation_texture() { return segmentation_tex; }

// @NOTE: This should generate a quad with a certain size in texture coordinates
constexpr const char* v_shader_src = R"(
#version 150 core

uniform int u_instance_offset = 0;
uniform samplerBuffer u_tex_coord;
uniform float u_radius;
out vec2 uv;

void main() {
	int VID = gl_VertexID;
	int IID = gl_InstanceID + u_instance_offset;

	vec2 coord = texelFetch(u_tex_coord, IID).xy;
	uv = vec2(VID / 2, VID % 2) * 2.0 - 1.0; 

	gl_Position = vec4(coord * 2.0 - 1.0 + uv * u_radius, 0, 1);
}
)";

// @NOTE: Do some radial falloff based on uv coordinate
constexpr const char* f_shader_src = R"(
#version 150 core

uniform vec4 u_color;
uniform float u_outline;
in vec2 uv;
out vec4 out_frag;

float step_dist(float edge, float dist) {
	float factor = 1.0; // <-- value can be played around with a bit
	float mask = step(edge, dist);
	float step_w = factor * length(vec2(dFdx(mask), dFdy(mask)));
	return smoothstep(-step_w/2.0, step_w/2.0, mask);
}

void main() {
	float dist = sqrt(dot(uv, uv));
	if (u_outline > 0) {
		vec4 inner_color = u_color.rgba;
		vec4 rim_color = vec4(0,0,0,1);
		vec4 outer_color = vec4(0,0,0,0);
		out_frag = mix(inner_color, mix(rim_color, outer_color, step_dist(1.0, dist)), step_dist(1.0 - u_outline, dist));
	} else {
		float falloff = max(0, 1.0 - dist);
		out_frag = vec4(u_color.rgb, u_color.a * falloff);	
	}
}
)";

void initialize() {
    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    GLuint v_shader = glCreateShader(GL_VERTEX_SHADER);
    GLuint f_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(v_shader, 1, &v_shader_src, 0);
    glShaderSource(f_shader, 1, &f_shader_src, 0);

    glCompileShader(v_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, v_shader)) {
        LOG_ERROR("Compiling ramachandran vertex shader:\n%s\n", buffer);
    }
    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        LOG_ERROR("Compiling ramachandran fragment shader:\n%s\n", buffer);
    }

    program = glCreateProgram();
    glAttachShader(program, v_shader);
    glAttachShader(program, f_shader);
    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, program)) {
        LOG_ERROR("Linking ramachandran program:\n%s\n", buffer);
    }

    glDetachShader(program, v_shader);
    glDetachShader(program, f_shader);

    glDeleteShader(v_shader);
    glDeleteShader(f_shader);

    uniform_loc_coord_tex = glGetUniformLocation(program, "u_coord_tex");
    uniform_loc_instance_offset = glGetUniformLocation(program, "u_instance_offset");
    uniform_loc_radius = glGetUniformLocation(program, "u_radius");
    uniform_loc_color = glGetUniformLocation(program, "u_color");
    uniform_loc_outline = glGetUniformLocation(program, "u_outline");

    if (!segmentation_tex) {
        glGenTextures(1, &segmentation_tex);
        glBindTexture(GL_TEXTURE_2D, segmentation_tex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, seg_width, seg_height, 0, seg_data_format, GL_UNSIGNED_BYTE, seg_data);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    if (!accumulation_tex) {
        glGenTextures(1, &accumulation_tex);
        glBindTexture(GL_TEXTURE_2D, accumulation_tex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, acc_width, acc_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    if (!coord_buf) {
        glGenBuffers(1, &coord_buf);
    }

    if (!coord_tex) {
        glGenTextures(1, &coord_tex);
    }

    if (!fbo) {
        glGenFramebuffers(1, &fbo);
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, accumulation_tex, 0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }
}

void shutdown() {
    if (segmentation_tex) glDeleteTextures(1, &segmentation_tex);
    if (accumulation_tex) glDeleteTextures(1, &accumulation_tex);
    if (coord_buf) glDeleteBuffers(1, &coord_buf);
    if (coord_tex) glDeleteTextures(1, &coord_tex);
    if (fbo) glDeleteFramebuffers(1, &fbo);
}

void clear_accumulation_texture() {
    GLint last_viewport[4];
    glGetIntegerv(GL_VIEWPORT, last_viewport);

    glViewport(0, 0, acc_width, acc_height);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);

    glViewport(last_viewport[0], last_viewport[1], (GLsizei)last_viewport[2], (GLsizei)last_viewport[3]);
}

void compute_accumulation_texture(Array<const BackboneAngles> angles, vec4 color, float radius, float outline) {
    struct Coord {
        unsigned short x, y;
    };

    // Use fast scratch memory here
    Coord* coords = (Coord*)TMP_MALLOC((angles.count) * sizeof(Coord));

    constexpr float ONE_OVER_TWO_PI = 1.f / (2.f * math::PI);

    int32 count = 0;
    for (const auto& angle : angles) {
        if (angle.phi == 0 || angle.psi == 0) continue;
        vec2 coord = vec2(angle.phi, angle.psi) * ONE_OVER_TWO_PI + 0.5f;  // [-PI, PI] -> [0, 1]
        coord.y = 1.f - coord.y;
        coords[count].x = (unsigned short)(coord.x * 0xffff);
        coords[count].y = (unsigned short)(coord.y * 0xffff);
        count++;
    }

    draw::set_vbo_data(coords, count * 2 * sizeof(unsigned short));

    TMP_FREE(coords);

    // Backup GL state
    GLint last_polygon_mode[2];
    glGetIntegerv(GL_POLYGON_MODE, last_polygon_mode);
    GLint last_viewport[4];
    glGetIntegerv(GL_VIEWPORT, last_viewport);
    GLenum last_blend_src_rgb;
    glGetIntegerv(GL_BLEND_SRC_RGB, (GLint*)&last_blend_src_rgb);
    GLenum last_blend_dst_rgb;
    glGetIntegerv(GL_BLEND_DST_RGB, (GLint*)&last_blend_dst_rgb);
    GLenum last_blend_src_alpha;
    glGetIntegerv(GL_BLEND_SRC_ALPHA, (GLint*)&last_blend_src_alpha);
    GLenum last_blend_dst_alpha;
    glGetIntegerv(GL_BLEND_DST_ALPHA, (GLint*)&last_blend_dst_alpha);
    GLenum last_blend_equation_rgb;
    glGetIntegerv(GL_BLEND_EQUATION_RGB, (GLint*)&last_blend_equation_rgb);
    GLenum last_blend_equation_alpha;
    glGetIntegerv(GL_BLEND_EQUATION_ALPHA, (GLint*)&last_blend_equation_alpha);
    GLboolean last_enable_blend = glIsEnabled(GL_BLEND);
    GLboolean last_enable_cull_face = glIsEnabled(GL_CULL_FACE);
    GLboolean last_enable_depth_test = glIsEnabled(GL_DEPTH_TEST);
    GLboolean last_enable_scissor_test = glIsEnabled(GL_SCISSOR_TEST);

    // RENDER TO ACCUMULATION TEXTURE

    glViewport(0, 0, acc_width, acc_height);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // Texture 0
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, coord_tex);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RG16, draw::vbo);

    glUseProgram(program);
    glUniform1i(uniform_loc_coord_tex, 0);

    // Draw
    glUniform1f(uniform_loc_radius, radius * 0.01f);
    glUniform1i(uniform_loc_instance_offset, 0);
    glUniform4fv(uniform_loc_color, 1, &color[0]);
    glUniform1f(uniform_loc_outline, outline);
    draw::draw_instanced_quads(count);

    glUseProgram(0);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

    // Restore modified GL state
    glBlendEquationSeparate(last_blend_equation_rgb, last_blend_equation_alpha);
    glBlendFuncSeparate(last_blend_src_rgb, last_blend_dst_rgb, last_blend_src_alpha, last_blend_dst_alpha);
    if (last_enable_blend)
        glEnable(GL_BLEND);
    else
        glDisable(GL_BLEND);
    if (last_enable_cull_face)
        glEnable(GL_CULL_FACE);
    else
        glDisable(GL_CULL_FACE);
    if (last_enable_depth_test)
        glEnable(GL_DEPTH_TEST);
    else
        glDisable(GL_DEPTH_TEST);
    if (last_enable_scissor_test)
        glEnable(GL_SCISSOR_TEST);
    else
        glDisable(GL_SCISSOR_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, last_polygon_mode[0]);
    glViewport(last_viewport[0], last_viewport[1], (GLsizei)last_viewport[2], (GLsizei)last_viewport[3]);
}

}  // namespace ramachandran

mat4 compute_linear_transform(Array<const vec3> pos_frame_a, Array<const vec3> pos_frame_b) {
    ASSERT(pos_frame_a.count == pos_frame_b.count);

    const vec3 com_a = compute_com(pos_frame_a);
    DynamicArray<vec3> q(pos_frame_a.count);
    for (int32 i = 0; i < q.count; i++) {
        q[i] = pos_frame_a[i] - com_a;
    }

    const vec3 com_b = compute_com(pos_frame_b);
    DynamicArray<vec3> p(pos_frame_b.count);
    for (int32 i = 0; i < p.count; i++) {
        p[i] = pos_frame_b[i] - com_b;
    }

    mat3 Apq{0};
    for (int32 i = 0; i < p.count; i++) {
        Apq[0][0] += p[i].x * q[i].x;
        Apq[0][1] += p[i].y * q[i].x;
        Apq[0][2] += p[i].z * q[i].x;
        Apq[1][0] += p[i].x * q[i].y;
        Apq[1][1] += p[i].y * q[i].y;
        Apq[1][2] += p[i].z * q[i].y;
        Apq[2][0] += p[i].x * q[i].z;
        Apq[2][1] += p[i].y * q[i].z;
        Apq[2][2] += p[i].z * q[i].z;
    }

    mat3 Aqq{0};
    for (int32 i = 0; i < q.count; i++) {
        Aqq[0][0] += q[i].x * q[i].x;
        Aqq[0][1] += q[i].y * q[i].x;
        Aqq[0][2] += q[i].z * q[i].x;
        Aqq[1][0] += q[i].x * q[i].y;
        Aqq[1][1] += q[i].y * q[i].y;
        Aqq[1][2] += q[i].z * q[i].y;
        Aqq[2][0] += q[i].x * q[i].z;
        Aqq[2][1] += q[i].y * q[i].z;
        Aqq[2][2] += q[i].z * q[i].z;
    }

    mat4 result = Apq / Aqq;
    result[3] = vec4(com_b, 1);
    return result;
}

mat4 compute_linear_transform(Array<const vec3> pos_frame_a, Array<const vec3> pos_frame_b, Array<const float> mass) {
    ASSERT(pos_frame_a.count == pos_frame_b.count);
    ASSERT(mass.count == pos_frame_a.count);

    const vec3 com_a = compute_com(pos_frame_a, mass);
    DynamicArray<vec3> q(pos_frame_a.count);
    for (int32 i = 0; i < q.count; i++) {
        q[i] = pos_frame_a[i] - com_a;
    }

    const vec3 com_b = compute_com(pos_frame_b, mass);
    DynamicArray<vec3> p(pos_frame_b.count);
    for (int32 i = 0; i < p.count; i++) {
        p[i] = pos_frame_b[i] - com_b;
    }

    mat3 Apq{0};
    for (int32 i = 0; i < p.count; i++) {
        Apq[0][0] += mass[i] * p[i].x * q[i].x;
        Apq[0][1] += mass[i] * p[i].y * q[i].x;
        Apq[0][2] += mass[i] * p[i].z * q[i].x;
        Apq[1][0] += mass[i] * p[i].x * q[i].y;
        Apq[1][1] += mass[i] * p[i].y * q[i].y;
        Apq[1][2] += mass[i] * p[i].z * q[i].y;
        Apq[2][0] += mass[i] * p[i].x * q[i].z;
        Apq[2][1] += mass[i] * p[i].y * q[i].z;
        Apq[2][2] += mass[i] * p[i].z * q[i].z;
    }

    mat3 Aqq{0};
    for (int32 i = 0; i < q.count; i++) {
        Aqq[0][0] += mass[i] * q[i].x * q[i].x;
        Aqq[0][1] += mass[i] * q[i].y * q[i].x;
        Aqq[0][2] += mass[i] * q[i].z * q[i].x;
        Aqq[1][0] += mass[i] * q[i].x * q[i].y;
        Aqq[1][1] += mass[i] * q[i].y * q[i].y;
        Aqq[1][2] += mass[i] * q[i].z * q[i].y;
        Aqq[2][0] += mass[i] * q[i].x * q[i].z;
        Aqq[2][1] += mass[i] * q[i].y * q[i].z;
        Aqq[2][2] += mass[i] * q[i].z * q[i].z;
    }

    mat4 result = Apq / Aqq;
    result[3] = vec4(com_b, 1);
    return result;
}

void compute_RS(mat3* R, mat3* S, Array<const vec3> x0, Array<const vec3> x, Array<const float> m) {
    ASSERT(x0.count == x.count);
    ASSERT(m.count == x0.count);

    const vec3 com_x0 = compute_com(x0, m);
    DynamicArray<vec3> q(x0.count);
    for (int32 i = 0; i < q.count; i++) {
        q[i] = x0[i] - com_x0;
    }

    const vec3 com_x = compute_com(x, m);
    DynamicArray<vec3> p(x.count);
    for (int32 i = 0; i < p.count; i++) {
        p[i] = x[i] - com_x;
    }

    mat3 Apq{0};
    for (int32 i = 0; i < p.count; i++) {
        Apq[0][0] += m[i] * p[i].x * q[i].x;
        Apq[0][1] += m[i] * p[i].y * q[i].x;
        Apq[0][2] += m[i] * p[i].z * q[i].x;
        Apq[1][0] += m[i] * p[i].x * q[i].y;
        Apq[1][1] += m[i] * p[i].y * q[i].y;
        Apq[1][2] += m[i] * p[i].z * q[i].y;
        Apq[2][0] += m[i] * p[i].x * q[i].z;
        Apq[2][1] += m[i] * p[i].y * q[i].z;
        Apq[2][2] += m[i] * p[i].z * q[i].z;
    }

    mat3 Q, D;
    diagonalize(math::transpose(Apq) * Apq, &Q, &D);
    D[0][0] = sqrtf(D[0][0]);
    D[1][1] = sqrtf(D[1][1]);
    D[2][2] = sqrtf(D[2][2]);

    *S = Q * D * math::inverse(Q);
    *R = Apq * math::inverse(*S);
}

// from here https://stackoverflow.com/questions/4372224/fast-method-for-computing-3x3-symmetric-matrix-spectral-decomposition
// Slightly modified version of  Stan Melax's code for 3x3 matrix diagonalization (Thanks Stan!)
// source: http://www.melax.com/diag.html?attredirects=0
void Diagonalize(const float (&A)[3][3], float (&Q)[3][3], float (&D)[3][3]) {
    // A must be a symmetric matrix.
    // returns Q and D such that
    // Diagonal matrix D = QT * A * Q;  and  A = Q*D*QT
    const int maxsteps = 24;  // certainly wont need that many.
    int k0, k1, k2;
    float o[3], m[3];
    float q[4] = {0.0f, 0.0f, 0.0f, 1.0f};
    float jr[4];
    float sqw, sqx, sqy, sqz;
    float tmp1, tmp2, mq;
    float AQ[3][3];
    float thet, sgn, t, c;
    for (int i = 0; i < maxsteps; ++i) {
        // quat to matrix
        sqx = q[0] * q[0];
        sqy = q[1] * q[1];
        sqz = q[2] * q[2];
        sqw = q[3] * q[3];
        Q[0][0] = (sqx - sqy - sqz + sqw);
        Q[1][1] = (-sqx + sqy - sqz + sqw);
        Q[2][2] = (-sqx - sqy + sqz + sqw);
        tmp1 = q[0] * q[1];
        tmp2 = q[2] * q[3];
        Q[1][0] = 2.0f * (tmp1 + tmp2);
        Q[0][1] = 2.0f * (tmp1 - tmp2);
        tmp1 = q[0] * q[2];
        tmp2 = q[1] * q[3];
        Q[2][0] = 2.0f * (tmp1 - tmp2);
        Q[0][2] = 2.0f * (tmp1 + tmp2);
        tmp1 = q[1] * q[2];
        tmp2 = q[0] * q[3];
        Q[2][1] = 2.0f * (tmp1 + tmp2);
        Q[1][2] = 2.0f * (tmp1 - tmp2);

        // AQ = A * Q
        AQ[0][0] = Q[0][0] * A[0][0] + Q[1][0] * A[0][1] + Q[2][0] * A[0][2];
        AQ[0][1] = Q[0][1] * A[0][0] + Q[1][1] * A[0][1] + Q[2][1] * A[0][2];
        AQ[0][2] = Q[0][2] * A[0][0] + Q[1][2] * A[0][1] + Q[2][2] * A[0][2];
        AQ[1][0] = Q[0][0] * A[0][1] + Q[1][0] * A[1][1] + Q[2][0] * A[1][2];
        AQ[1][1] = Q[0][1] * A[0][1] + Q[1][1] * A[1][1] + Q[2][1] * A[1][2];
        AQ[1][2] = Q[0][2] * A[0][1] + Q[1][2] * A[1][1] + Q[2][2] * A[1][2];
        AQ[2][0] = Q[0][0] * A[0][2] + Q[1][0] * A[1][2] + Q[2][0] * A[2][2];
        AQ[2][1] = Q[0][1] * A[0][2] + Q[1][1] * A[1][2] + Q[2][1] * A[2][2];
        AQ[2][2] = Q[0][2] * A[0][2] + Q[1][2] * A[1][2] + Q[2][2] * A[2][2];
        // D = Qt * AQ
        D[0][0] = AQ[0][0] * Q[0][0] + AQ[1][0] * Q[1][0] + AQ[2][0] * Q[2][0];
        D[0][1] = AQ[0][0] * Q[0][1] + AQ[1][0] * Q[1][1] + AQ[2][0] * Q[2][1];
        D[0][2] = AQ[0][0] * Q[0][2] + AQ[1][0] * Q[1][2] + AQ[2][0] * Q[2][2];
        D[1][0] = AQ[0][1] * Q[0][0] + AQ[1][1] * Q[1][0] + AQ[2][1] * Q[2][0];
        D[1][1] = AQ[0][1] * Q[0][1] + AQ[1][1] * Q[1][1] + AQ[2][1] * Q[2][1];
        D[1][2] = AQ[0][1] * Q[0][2] + AQ[1][1] * Q[1][2] + AQ[2][1] * Q[2][2];
        D[2][0] = AQ[0][2] * Q[0][0] + AQ[1][2] * Q[1][0] + AQ[2][2] * Q[2][0];
        D[2][1] = AQ[0][2] * Q[0][1] + AQ[1][2] * Q[1][1] + AQ[2][2] * Q[2][1];
        D[2][2] = AQ[0][2] * Q[0][2] + AQ[1][2] * Q[1][2] + AQ[2][2] * Q[2][2];
        o[0] = D[1][2];
        o[1] = D[0][2];
        o[2] = D[0][1];
        m[0] = fabs(o[0]);
        m[1] = fabs(o[1]);
        m[2] = fabs(o[2]);

        k0 = (m[0] > m[1] && m[0] > m[2]) ? 0 : (m[1] > m[2]) ? 1 : 2;  // index of largest element of offdiag
        k1 = (k0 + 1) % 3;
        k2 = (k0 + 2) % 3;
        if (o[k0] == 0.0f) {
            break;  // diagonal already
        }
        thet = (D[k2][k2] - D[k1][k1]) / (2.0f * o[k0]);
        sgn = (thet > 0.0f) ? 1.0f : -1.0f;
        thet *= sgn;                                                             // make it positive
        t = sgn / (thet + ((thet < 1.E6f) ? sqrtf(thet * thet + 1.0f) : thet));  // sign(T)/(|T|+sqrt(T^2+1))
        c = 1.0f / sqrtf(t * t + 1.0f);                                          //  c= 1/(t^2+1) , t=s/c
        if (c == 1.0f) {
            break;  // no room for improvement - reached machine precision.
        }
        jr[0] = jr[1] = jr[2] = jr[3] = 0.0f;
        jr[k0] = sgn * sqrtf((1.0f - c) / 2.0f);  // using 1/2 angle identity sin(a/2) = sqrt((1-cos(a))/2)
        jr[k0] *= -1.0f;                          // since our quat-to-matrix convention was for v*M instead of M*v
        jr[3] = sqrtf(1.0f - jr[k0] * jr[k0]);
        if (jr[3] == 1.0f) {
            break;  // reached limits of floating point precision
        }
        q[0] = (q[3] * jr[0] + q[0] * jr[3] + q[1] * jr[2] - q[2] * jr[1]);
        q[1] = (q[3] * jr[1] - q[0] * jr[2] + q[1] * jr[3] + q[2] * jr[0]);
        q[2] = (q[3] * jr[2] + q[0] * jr[1] - q[1] * jr[0] + q[2] * jr[3]);
        q[3] = (q[3] * jr[3] - q[0] * jr[0] - q[1] * jr[1] - q[2] * jr[2]);
        mq = sqrtf(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
        q[0] /= mq;
        q[1] /= mq;
        q[2] /= mq;
        q[3] /= mq;
    }
}

void diagonalize(const mat3& M, mat3* Q, mat3* D) {
    ASSERT(Q);
    ASSERT(D);
    Diagonalize((const float(&)[3][3])M, (float(&)[3][3]) * Q, (float(&)[3][3]) * D);
}

void decompose(const mat3& M, mat3* R, mat3* S) {
    ASSERT(R);
    ASSERT(S);
    mat3 Q, D;
    mat3 AtA = math::transpose(M) * M;
    diagonalize(AtA, &Q, &D);
    float det = math::determinant(AtA);
    D[0][0] = sqrtf(D[0][0]);
    D[1][1] = sqrtf(D[1][1]);
    D[2][2] = sqrtf(D[2][2]);
    *S = Q * D * math::inverse(Q);
    *R = M * math::inverse(*S);
}
