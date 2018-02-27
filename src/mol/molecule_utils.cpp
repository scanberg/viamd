#include "molecule_utils.h"
#include <gfx/gl_utils.h>
#include <gfx/immediate_draw_utils.h>
#include <core/hash.h>
#include <core/common.h>
#include <core/math_utils.h>
#include <mol/aminoacid.h>
#include <mol/trajectory_utils.h>
#include <imgui.h>
#include <ctype.h>

void transform_positions(Array<vec3> positions, const mat4& transformation) {
    for (auto& p : positions) {
        p = vec3(transformation * vec4(p, 1));
    }
}

void compute_bounding_box(vec3* min_box, vec3* max_box, const Array<vec3> positions) {
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

// @TODO: Fix this, is it possible in theory to get a good interpolation between frames with periodicity without modifying source data?
void linear_interpolation_periodic(Array<vec3> positions, const Array<vec3> prev_pos, const Array<vec3> next_pos, float t, mat3 sim_box) {
    ASSERT(prev_pos.count == positions.count);
    ASSERT(next_pos.count == positions.count);

    vec3 full_box_ext = sim_box * vec3(1);
    vec3 half_box_ext = full_box_ext * 0.5f;

    for (int i = 0; i < positions.count; i++) {
        vec3 next = next_pos[i];
        vec3 prev = prev_pos[i];

        vec3 delta = next - prev;
        vec3 sign = math::sign(delta);

        vec3 abs_delta = math::abs(delta);
        if (abs_delta.x > half_box_ext.x) next.x -= full_box_ext.x * sign.x;
        if (abs_delta.y > half_box_ext.y) next.y -= full_box_ext.y * sign.y;
        if (abs_delta.z > half_box_ext.z) next.z -= full_box_ext.z * sign.z;

        // next -= math::step(half_box_ext, math::abs(delta)) * delta;

        positions[i] = math::mix(prev, next, t);
    }
}

void linear_interpolation(Array<vec3> positions, const Array<vec3> prev_pos, const Array<vec3> next_pos, float t) {
    ASSERT(prev_pos.count == positions.count);
    ASSERT(next_pos.count == positions.count);

    for (int i = 0; i < positions.count; i++) {
        positions[i] = math::mix(prev_pos[i], next_pos[i], t);
    }
}

DynamicArray<Bond> compute_covalent_bonds(const Array<vec3> atom_pos, const Array<Element> atom_elem, const Array<Residue> residues) {
    ASSERT(atom_pos.count == atom_elem.count);

    DynamicArray<Bond> bonds;

    auto try_and_create_atom_bond = [& pos = atom_pos, &elem = atom_elem, &bonds = bonds ](int i, int j)->bool {
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

/*
static inline bool idx_in_residue(int idx, const Residue& res) { return res.beg_atom_idx <= idx && idx < res.end_atom_idx; }

static bool residues_are_connected(Residue res_a, Residue res_b, const Array<Bond> bonds) {
        for (const auto& bond : bonds) {
                if (idx_in_residue(bond.idx_a, res_a) && idx_in_residue(bond.idx_b, res_b)) return true;
                if (idx_in_residue(bond.idx_b, res_a) && idx_in_residue(bond.idx_a, res_b)) return true;
        }
        return false;
}
*/

// @NOTE this method is sub-optimal and can surely be improved...
// Residues should have no more than 2 potential connections to other residues.
DynamicArray<Chain> compute_chains(const Array<Residue> residues, const Array<Bond> bonds, const Array<int32> atom_residue_indices) {
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

    DynamicArray<int> residue_chains(residues.count, -1);

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

    DynamicArray<Chain> chains;
    curr_chain_idx = -1;
    for (int i = 0; i < residue_chains.count; i++) {
        if (residue_chains[i] != curr_chain_idx) {
            curr_chain_idx = residue_chains[i];
            Label lbl;
            snprintf(lbl.beg(), Label::MAX_LENGTH, "C%i", curr_chain_idx);
            chains.push_back({lbl, i, i});
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

DynamicArray<BackboneSegment> compute_backbone(const Chain& chain, const Array<Residue> residues, const Array<Label> atom_labels) {
    DynamicArray<BackboneSegment> backbones;
    for (int32 res_idx = chain.beg_res_idx; res_idx < chain.end_res_idx; res_idx++) {
        const auto& residue = residues[res_idx];
        if (is_amino_acid(residue) == false) continue;
        // find atoms
        auto ca_idx = -1;
        auto ha_idx = -1;
        auto cb_idx = -1;
        auto n_idx = -1;
        auto c_idx = -1;
        auto o_idx = -1;
        for (int32 i = residue.beg_atom_idx; i < residue.end_atom_idx; i++) {
            const auto& lbl = atom_labels[i];
            if (ca_idx == -1 && match(lbl, "CA")) ca_idx = i;
            if (ha_idx == -1 && match(lbl, "HA")) ha_idx = i;
            if (cb_idx == -1 && match(lbl, "CB")) cb_idx = i;
            if (n_idx == -1 && match(lbl, "N")) n_idx = i;
            if (c_idx == -1 && match(lbl, "C")) c_idx = i;
            if (o_idx == -1 && match(lbl, "O")) o_idx = i;
        }
        if (ca_idx == -1) {
            printf("No CA label found for residue[%i]: %s.\n", res_idx, residues[res_idx].id.beg());
        }
        if (ha_idx == -1) {
            printf("No HA label found for residue[%i]: %s.\n", res_idx, residues[res_idx].id.beg());
        }
        if (cb_idx == -1) {
            printf("No CB label found for residue[%i]: %s.\n", res_idx, residues[res_idx].id.beg());
        }
        if (n_idx == -1) {
            printf("No N label found for residue[%i]: %s.\n", res_idx, residues[res_idx].id.beg());
        }
        if (c_idx == -1) {
            printf("No C label found for residue[%i]: %s.\n", res_idx, residues[res_idx].id.beg());
        }
        if (o_idx == -1) {
            printf("No O label found for residue[%i]: %s.\n", res_idx, residues[res_idx].id.beg());
        }

        backbones.push_back({ca_idx, ha_idx, cb_idx, n_idx, c_idx, o_idx});
    }

    return backbones;
}

DynamicArray<SplineSegment> compute_spline(const Array<vec3> atom_pos, const Array<BackboneSegment>& backbone, int num_subdivisions) {
    if (backbone.count < 4) return {};

    // @TODO: Use C -> O vector for orientation

    DynamicArray<vec3> p_tmp;
    DynamicArray<vec3> o_tmp;
    DynamicArray<vec3> c_tmp;
    DynamicArray<int> ca_idx;

    auto d_p0 = atom_pos[backbone[1].ca_idx] - atom_pos[backbone[0].ca_idx];
    p_tmp.push_back(atom_pos[backbone[0].ca_idx] - d_p0);

    auto d_o0 = atom_pos[backbone[1].o_idx] - atom_pos[backbone[0].o_idx];
    o_tmp.push_back(atom_pos[backbone[0].o_idx] - d_o0);

    if (backbone[0].c_idx > -1 && backbone[1].c_idx > -1) {
        auto d_c0 = atom_pos[backbone[1].c_idx] - atom_pos[backbone[0].c_idx];
        c_tmp.push_back(atom_pos[backbone[0].c_idx] - d_c0);
    } else {
        // This should never happen
        c_tmp.push_back(atom_pos[backbone[0].ca_idx] - d_p0);
    }

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

    const float tension = 0.5f;

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

        for (int n = 0; n < num_subdivisions; n++) {
            auto t = n / (float)(num_subdivisions);

            vec3 p = math::spline(p0, p1, p2, p3, t, tension);
            vec3 o = math::spline(o0, o1, o2, o3, t, tension);
            vec3 c = math::spline(c0, c1, c2, c3, t, tension);

            vec3 v_dir = math::normalize(o - c);

            const float eps = 0.0001f;
            float d0 = math::max(0.f, t - eps);
            float d1 = math::min(t + eps, 1.f);

            vec3 tangent = math::normalize(math::spline(p0, p1, p2, p3, d1, tension) - math::spline(p0, p1, p2, p3, d0, tension));
            // vec3 binormal = math::normalize(math::cross(v_dir, tangent));
            // vec3 normal = math::normalize(math::cross(tangent, binormal));
            vec3 normal = math::normalize(v_dir);
            vec3 binormal = math::normalize(math::cross(tangent, normal));

            segments.push_back({p, tangent, normal, binormal});
        }
    }

    return segments;
}

DynamicArray<BackboneAngles> compute_backbone_angles(const Array<vec3> pos, const Array<BackboneSegment> backbone) {
	if (backbone.count == 0) return {};
	DynamicArray<BackboneAngles> angles(backbone.count);
	compute_backbone_angles(angles, pos, backbone);
	return angles;
}

void compute_backbone_angles(Array<BackboneAngles> dst, const Array<vec3> pos, const Array<BackboneSegment> backbone) {
	ASSERT(dst.count >= backbone.count);
	float omega, phi, psi;

	omega = 0;
	phi = 0;
	psi = compute_dihedral_angle(pos[backbone[0].n_idx], pos[backbone[0].ca_idx], pos[backbone[0].c_idx], pos[backbone[1].n_idx]);
	dst[0] = { omega, phi, psi };

	for (int64 i = 1; i < backbone.count - 1; i++) {
		omega = compute_dihedral_angle(pos[backbone[i - 1].ca_idx], pos[backbone[i - 1].c_idx], pos[backbone[i].n_idx], pos[backbone[i].ca_idx]);
		phi = compute_dihedral_angle(pos[backbone[i - 1].c_idx], pos[backbone[i].n_idx], pos[backbone[i].ca_idx], pos[backbone[i].c_idx]);
		psi = compute_dihedral_angle(pos[backbone[i].n_idx], pos[backbone[i].ca_idx], pos[backbone[i].c_idx], pos[backbone[i + 1].n_idx]);
		dst[i] = { omega, phi, psi };
	}

	auto N = backbone.count - 1;
	omega = compute_dihedral_angle(pos[backbone[N - 1].ca_idx], pos[backbone[N - 1].c_idx], pos[backbone[N].n_idx], pos[backbone[N].ca_idx]);
	phi = compute_dihedral_angle(pos[backbone[N - 1].c_idx], pos[backbone[N].n_idx], pos[backbone[N].ca_idx], pos[backbone[N].c_idx]);
	psi = 0;
	dst[N] = { omega, phi, psi };
}

BackboneAnglesTrajectory compute_backbone_angles_trajectory(const Trajectory& trajectory, const Array<BackboneSegment> backbone) {
	if (trajectory.num_frames == 0 || backbone.count == 0) return {};

	BackboneAnglesTrajectory bat;
	bat.num_frames = trajectory.num_frames;
	bat.num_segments = backbone.count;
	bat.angle_data = DynamicArray<BackboneAngles>(bat.num_frames * bat.num_segments);

	// @TODO: parallelize?
	for (int f_idx = 0; f_idx < trajectory.num_frames; f_idx++) {
		auto pos = get_trajectory_positions(trajectory, f_idx);
		auto b_angles = get_backbone_angles(bat, f_idx);
		compute_backbone_angles(b_angles, pos, backbone);
	}
}

DynamicArray<float> compute_atom_radii(const Array<Element> elements) {
    DynamicArray<float> radii(elements.count, 0);
    compute_atom_radii(radii, elements);
    return radii;
}

void compute_atom_radii(Array<float> radii_dst, const Array<Element> elements) {
    ASSERT(radii_dst.count <= elements.count);
    for (int64 i = 0; i < radii_dst.count; i++) {
        radii_dst[i] = element::vdw_radius(elements[i]);
    }
}

DynamicArray<uint32> compute_atom_colors(const MoleculeStructure& mol, ColorMapping mapping) {
    DynamicArray<uint32> colors(mol.atom_elements.count, 0xFFFFFFFF);
    compute_atom_colors(colors, mol, mapping);
    return colors;
}

void compute_atom_colors(Array<uint32> color_dst, const MoleculeStructure& mol, ColorMapping mapping) {
    // @TODO: Implement more mappings

    // CPK
    switch (mapping) {
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
                    unsigned int h = hash::crc32(res.id.beg(), Label::MAX_LENGTH);
                    float hue = (h % 32) / 32.f;
                    vec3 c = math::hcl_to_rgb(vec3(hue, 0.8f, 0.8f));
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
                    unsigned int h = hash::crc32(mol.atom_residue_indices[i]);
                    float hue = (h % 15) / 15.f;
                    vec3 c = math::hcl_to_rgb(vec3(hue, 0.8f, 0.8f));
                    unsigned char color[4];
                    color[0] = (unsigned char)(c.x * 255);
                    color[1] = (unsigned char)(c.y * 255);
                    color[2] = (unsigned char)(c.z * 255);
                    color[3] = (unsigned char)(255);
                    color_dst[i] = *(uint32*)(color);
                }
            }
            break;
        case ColorMapping::CHAIN_INDEX:
            for (int64 i = 0; i < color_dst.count; i++) {
                if (i < mol.atom_residue_indices.count) {
                    const auto& res = mol.residues[mol.atom_residue_indices[i]];
                    if (res.chain_idx < mol.chains.count) {
                        unsigned int h = hash::crc32(res.chain_idx);
                        float hue = (h % 32) / 32.f;
                        vec3 c = math::hcl_to_rgb(vec3(hue, 0.8f, 0.8f));
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

bool is_amino_acid(Residue res) { return aminoacid::get_from_string(res.id) != AminoAcid::Unknown; }

namespace draw {

static constexpr int VERTEX_BUFFER_SIZE = MEGABYTES(4);
static GLuint empty_vao = 0;
static GLuint vbo = 0;

namespace vdw {
static GLuint v_shader = 0;
static GLuint f_shader = 0;
static GLuint program = 0;

static GLuint vao = 0;
static GLuint ibo = 0;
static GLuint buf_position_radius = 0;
static GLuint buf_color = 0;
static GLuint tex_position_radius = 0;
static GLuint tex_color = 0;

static GLint uniform_loc_view_mat = -1;
static GLint uniform_loc_proj_mat = -1;
static GLint uniform_loc_inv_proj_mat = -1;
static GLint uniform_loc_fov = -1;
static GLint uniform_loc_tex_pos_rad = -1;
static GLint uniform_loc_tex_color = -1;

static const char* v_shader_src = R"(
#version 150 core

uniform mat4 u_view_mat;
uniform mat4 u_proj_mat;
uniform mat4 u_inv_proj_mat;

uniform samplerBuffer u_tex_pos_rad;
uniform samplerBuffer u_tex_color;

out Fragment {
    flat vec4 color;
    flat vec4 view_sphere;
	flat vec4 picking_color;
    smooth vec4 view_coord;
} out_frag;

vec4 pack_u32(uint data) {
	return vec4(
        (data & uint(0x000000FF)) >> 0,
        (data & uint(0x0000FF00)) >> 8,
        (data & uint(0x00FF0000)) >> 16,
        (data & uint(0xFF000000)) >> 24) / 255.0;
}

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
	vec4 color = texelFetch(u_tex_color, IID);

	vec3 pos = pos_rad.xyz;
	float rad = pos_rad.w;

    vec4 view_coord = u_view_mat * vec4(pos, 1.0);
    float len = length(view_coord.xyz);
    vec3 view_dir = view_coord.xyz / len;

    out_frag.color = color;
    out_frag.view_sphere = vec4(view_coord.xyz, rad);
	out_frag.picking_color = pack_u32(uint(IID));

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
layout(location = 1) out vec4 out_picking_id;

float fresnel(float H_dot_V) {   
    const float n1 = 1.0;
    const float n2 = 1.5;
    const float R0 = pow((n1-n2)/(n1+n2), 2);

    return R0 + (1.0 - R0)*pow(1.0 - H_dot_V, 5);
}

vec3 srgb_to_rgb_approx(vec3 srgb) {
    return pow(srgb, vec3(2.2));
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
    vec4 color = in_frag.color;

    // Compute Color
    const vec3 env_radiance = vec3(1.0);
    const vec3 dir_radiance = vec3(10.0);
    const vec3 L = normalize(vec3(1));
    const float spec_exp = 50.0;

    vec3 N = view_normal;
    vec3 V = view_dir;
    vec3 H = normalize(L + V);
    float H_dot_V = max(0.0, dot(H, V));
    float N_dot_H = max(0.0, dot(N, H));
    float N_dot_L = max(0.0, dot(N, L));
    float fr = fresnel(H_dot_V);

    vec3 diffuse = color.rgb * (env_radiance + N_dot_L * dir_radiance);
    vec3 specular = dir_radiance * pow(N_dot_H, spec_exp);

    color.rgb = mix(diffuse, specular, fr);

    gl_FragDepth = (-u_proj_mat[2][2] - u_proj_mat[3][2] / view_hit.z) * 0.5 + 0.5;
    out_color = vec4(color.rgb, color.a);
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
        printf("Error while compiling vdw vertex shader:\n%s\n", buffer);
    }

    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        printf("Error while compiling vdw fragment shader:\n%s\n", buffer);
    }

    program = glCreateProgram();
    glAttachShader(program, v_shader);
    glAttachShader(program, f_shader);
    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, program)) {
        printf("Error while linking vdw program:\n%s\n", buffer);
    }

    glDetachShader(program, v_shader);
    glDetachShader(program, f_shader);

    glDeleteShader(v_shader);
    glDeleteShader(f_shader);

    if (!vao) glGenVertexArrays(1, &vao);

    if (!ibo) {
        const unsigned char data[4] = {0, 1, 2, 3};
        glGenBuffers(1, &ibo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 4, data, GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }

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

    if (!tex_position_radius) glGenTextures(1, &tex_position_radius);

    if (!tex_color) glGenTextures(1, &tex_color);

    glBindVertexArray(vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBindVertexArray(0);

    uniform_loc_view_mat = glGetUniformLocation(program, "u_view_mat");
    uniform_loc_proj_mat = glGetUniformLocation(program, "u_proj_mat");
    uniform_loc_inv_proj_mat = glGetUniformLocation(program, "u_inv_proj_mat");
    uniform_loc_fov = glGetUniformLocation(program, "u_fov");
    uniform_loc_tex_pos_rad = glGetUniformLocation(vdw::program, "u_tex_pos_rad");
    uniform_loc_tex_color = glGetUniformLocation(vdw::program, "u_tex_color");
}

static void shutdown() {
    if (program) glDeleteProgram(program);
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
static GLint uniform_loc_radius_scl = -1;

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
uniform float u_radius_scl = 1.0;

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
	float r = 1.0 * u_radius_scl;
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
uniform float u_radius_scl = 1.0;

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
layout(location = 1) out vec4 out_picking_id;

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

vec4 pack_u32(uint data) {
	return vec4(
        (data & uint(0x000000FF)) >> 0,
        (data & uint(0x0000FF00)) >> 8,
        (data & uint(0x00FF0000)) >> 16,
        (data & uint(0xFF000000)) >> 24) / 255.0;
}

float fresnel(float H_dot_V) {   
    const float n1 = 1.0;
    const float n2 = 1.5;
    const float R0 = pow((n1-n2)/(n1+n2), 2);

    return R0 + (1.0 - R0)*pow(1.0 - H_dot_V, 5);
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

    // Compute Color
    const vec3 env_radiance = vec3(1.0);
    const vec3 dir_radiance = vec3(10.0);
    const vec3 L = normalize(vec3(1));
    const float spec_exp = 50.0;

    vec3 N = normal;
    vec3 V = -rd;
    vec3 H = normalize(L + V);
    float H_dot_V = max(0.0, dot(H, V));
    float N_dot_H = max(0.0, dot(N, H));
    float N_dot_L = max(0.0, dot(N, L));
    float fr = fresnel(H_dot_V);

    vec3 diffuse = color.rgb * (env_radiance + N_dot_L * dir_radiance);
    vec3 specular = dir_radiance * pow(N_dot_H, spec_exp);

    color.rgb = mix(diffuse, specular, fr);

    vec4 coord = vec4(0, 0, pos.z, 1);
    coord = u_proj_mat * coord;
    coord = coord / coord.w;

    gl_FragDepth = coord.z * 0.5 + 0.5;
    out_color = vec4(color.rgb, color.a);
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
        printf("Error while compiling licorice vertex shader:\n%s\n", buffer);
    }
    glCompileShader(g_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, g_shader)) {
        printf("Error while compiling licorice geometry shader:\n%s\n", buffer);
    }
    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        printf("Error while compiling licorice fragment shader:\n%s\n", buffer);
    }

    program = glCreateProgram();
    glAttachShader(program, v_shader);
    glAttachShader(program, g_shader);
    glAttachShader(program, f_shader);
    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, program)) {
        printf("Error while linking licorice program:\n%s\n", buffer);
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
    uniform_loc_radius_scl = glGetUniformLocation(program, "u_radius_scl");

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
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, VERTEX_BUFFER_SIZE, nullptr, GL_STREAM_DRAW);
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

#define CIRCLE_RES 32

layout(lines) in;
layout(triangle_strip, max_vertices = 256) out;

uniform mat4 u_view_mat;
uniform mat4 u_proj_mat;
uniform float scale_x = 1.0;
uniform float scale_y = 0.1;

in Vertex {
    vec3 tangent;
    vec3 normal;
    vec4 color;
    vec4 picking_color;
} in_vert[];

out Fragment {
    smooth vec3 view_position;
    smooth vec3 view_normal;
    smooth vec4 color;
    flat vec4 picking_color;
} out_frag;

void emit(mat4 mat, int input_idx, vec4 v) {
    vec4 view_coord = u_view_mat * mat * v;
    mat3 norm_mat = inverse(transpose(mat3(u_view_mat) * mat3(mat)));
    out_frag.view_position = view_coord.xyz;
    out_frag.view_normal = normalize(norm_mat * v.xyz);
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

    const float delta_angle = 6.28318530718 / float(CIRCLE_RES);
    for (int u = 0; u <= CIRCLE_RES; u++) {
        float angle = delta_angle * u;
        vec4 v = vec4(scale_x * cos(angle), scale_y * sin(angle), 0, 1);
        // disc 0
        // FIX THIS COLOR INTERPOLATION
        emit(mat[0], 0, v);
        // disc 1
        emit(mat[1], 1, v);
    }
    EndPrimitive();
}
)";

static const char* f_shader_src = R"(
#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

in Fragment {
	smooth vec3 view_position;
    smooth vec3 view_normal;
    smooth vec4 color;
	flat vec4 picking_color;
} in_frag;

layout(location = 0) out vec4 out_color;
layout(location = 1) out vec4 out_picking_id;

void main() {
    vec4 color = in_frag.color;
    vec3 V = -normalize(in_frag.view_position);
	vec3 N = normalize(in_frag.view_normal);
	vec3 L = normalize(vec3(1,1,1));

    color.rgb = color.rgb * max(0.0, dot(N, L));

    out_color = vec4(color.rgb, color.a);
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
        printf("Error while compiling ribbons vertex shader:\n%s\n", buffer);
    }
    glCompileShader(g_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, g_shader)) {
        printf("Error while compiling ribbons geometry shader:\n%s\n", buffer);
    }
    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        printf("Error while compiling ribbons fragment shader:\n%s\n", buffer);
    }

    program = glCreateProgram();
    glAttachShader(program, v_shader);
    glAttachShader(program, g_shader);
    glAttachShader(program, f_shader);
    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, program)) {
        printf("Error while linking ribbons program:\n%s\n", buffer);
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

}

void initialize() {
    glGenVertexArrays(1, &empty_vao);
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, VERTEX_BUFFER_SIZE, nullptr, GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    vdw::initialize();
    licorice::initialize();
	ribbons::intitialize();
}

void shutdown() {
    if (empty_vao) glDeleteVertexArrays(1, &empty_vao);
    if (vbo) glDeleteBuffers(1, &vbo);

    vdw::shutdown();
    licorice::shutdown();
	ribbons::shutdown();
}

void draw_vdw(const Array<vec3> atom_positions, const Array<float> atom_radii, const Array<uint32> atom_colors, const mat4& view_mat,
              const mat4& proj_mat, float radii_scale) {
    int64_t count = atom_positions.count;
    ASSERT(count == atom_radii.count && count == atom_colors.count);

    mat4 inv_proj_mat = math::inverse(proj_mat);

    unsigned int draw_count = 0;
    static bool initialized = false;

    // if (!initialized) {
    glBindBuffer(GL_ARRAY_BUFFER, vdw::buf_position_radius);
    vec4* gpu_pos_rad = (vec4*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    glBindBuffer(GL_ARRAY_BUFFER, vdw::buf_color);
    uint32* gpu_color = (uint32*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

    // @ TODO: DISCARD ANY ZERO RADII OR ZERO COLOR ALPHA ATOMS HERE
    for (int64_t i = 0; i < count; i++) {
        if (atom_radii[i] <= 0.f) continue;
        if ((atom_colors[i] & 0xff000000) == 0) continue;
        gpu_pos_rad[i] = vec4(atom_positions[i], atom_radii[i] * radii_scale);
        gpu_color[i] = atom_colors[i];
        draw_count++;
    }

    glUnmapBuffer(GL_ARRAY_BUFFER);
    glBindBuffer(GL_ARRAY_BUFFER, vdw::buf_position_radius);
    glUnmapBuffer(GL_ARRAY_BUFFER);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    initialized = true;
    //}

    glEnable(GL_DEPTH_TEST);

    glBindVertexArray(vdw::vao);
    glUseProgram(vdw::program);

    // Texture 0
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, vdw::tex_position_radius);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32F, vdw::buf_position_radius);

    // Texture 1
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, vdw::tex_color);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, vdw::buf_color);

    // Uniforms
    glUniform1i(vdw::uniform_loc_tex_pos_rad, 0);
    glUniform1i(vdw::uniform_loc_tex_color, 1);
    glUniformMatrix4fv(vdw::uniform_loc_view_mat, 1, GL_FALSE, &view_mat[0][0]);
    glUniformMatrix4fv(vdw::uniform_loc_proj_mat, 1, GL_FALSE, &proj_mat[0][0]);
    glUniformMatrix4fv(vdw::uniform_loc_inv_proj_mat, 1, GL_FALSE, &inv_proj_mat[0][0]);

    // Draw call
    glDrawElementsInstanced(GL_TRIANGLE_STRIP, 4, GL_UNSIGNED_BYTE, 0, draw_count);

    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glDisable(GL_DEPTH_TEST);
}

void draw_licorice(const Array<vec3> atom_positions, const Array<Bond> atom_bonds, const Array<uint32> atom_colors, const mat4& view_mat,
                   const mat4& proj_mat, float radii_scale) {
    int64_t count = atom_positions.count;
    ASSERT(count == atom_colors.count);
    ASSERT(count * sizeof(licorice::Vertex) < VERTEX_BUFFER_SIZE);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    licorice::Vertex* data = (licorice::Vertex*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (int64_t i = 0; i < count; i++) {
        data[i].position = atom_positions[i];
        data[i].color = atom_colors[i];
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, licorice::ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, atom_bonds.count * sizeof(Bond), atom_bonds.data, GL_STREAM_DRAW);

    glEnable(GL_DEPTH_TEST);

    glBindVertexArray(licorice::vao);
    glUseProgram(licorice::program);
    glUniformMatrix4fv(licorice::uniform_loc_view_mat, 1, GL_FALSE, &view_mat[0][0]);
    glUniformMatrix4fv(licorice::uniform_loc_proj_mat, 1, GL_FALSE, &proj_mat[0][0]);
    glUniform1f(licorice::uniform_loc_radius_scl, radii_scale);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, licorice::ibo);
    glDrawElements(GL_LINES, (GLsizei)atom_bonds.count * 2, GL_UNSIGNED_INT, (const void*)0);
    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glDisable(GL_DEPTH_TEST);
}

void draw_ribbons(const Array<SplineSegment> spline, const mat4& view_mat, const mat4& proj_mat) {
	ASSERT(spline.count * sizeof(ribbons::Vertex) < VERTEX_BUFFER_SIZE);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	ribbons::Vertex* data = (ribbons::Vertex*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
	for (int64_t i = 0; i < spline.count; i++) {
		data[i].position = spline[i].position;
		data[i].tangent = spline[i].tangent;
		data[i].normal = spline[i].normal;
		data[i].color = *((uint32*)immediate::COLOR_WHITE);
		data[i].picking_id = 0xffffffff;
	}
	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glEnable(GL_DEPTH_TEST);

	glBindVertexArray(ribbons::vao);
	glUseProgram(ribbons::program);
	glUniformMatrix4fv(ribbons::uniform_loc_view_mat, 1, GL_FALSE, &view_mat[0][0]);
	glUniformMatrix4fv(ribbons::uniform_loc_proj_mat, 1, GL_FALSE, &proj_mat[0][0]);
	glUniform1f(ribbons::uniform_loc_scale_x, 2.0f);
	glUniform1f(ribbons::uniform_loc_scale_y, 0.5f);
	glDrawArrays(GL_LINE_STRIP, 0, spline.count);
	glUseProgram(0);
	glBindVertexArray(0);

	glDisable(GL_DEPTH_TEST);
}

void draw_backbone(const Array<BackboneSegment> backbone, const Array<vec3> atom_positions, const mat4& view_mat, const mat4& proj_mat) {
    immediate::set_view_matrix(view_mat);
    immediate::set_proj_matrix(proj_mat);

    for (const auto& seg : backbone) {
        if (seg.ca_idx > -1) immediate::draw_point(atom_positions[seg.ca_idx], immediate::COLOR_BLACK);
        if (seg.c_idx > -1) immediate::draw_point(atom_positions[seg.c_idx], immediate::COLOR_GREEN);
        if (seg.o_idx > -1) immediate::draw_point(atom_positions[seg.o_idx], immediate::COLOR_RED);
    }

    for (int i = 1; i < backbone.count; i++) {
        if (backbone[i].ca_idx > -1 && backbone[i - 1].ca_idx > -1)
            immediate::draw_line(atom_positions[backbone[i - 1].ca_idx], atom_positions[backbone[i].ca_idx], immediate::COLOR_WHITE);
    }

    immediate::flush();
}

void draw_spline(const Array<SplineSegment> spline, const mat4& view_mat, const mat4& proj_mat) {
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