#include <core/array_types.h>
#include <core/hash.h>
#include <mol/molecule_structure.h>
#include "color_utils.h"

void color_atoms_cpk(Array<uint32> dst_atom_colors, Array<const Element> elements) {
    for (int64 i = 0; i < dst_atom_colors.count; i++) {
        dst_atom_colors[i] = element::color(elements[i]);
    }
}

void color_atoms_residue_id(Array<uint32> dst_atom_colors, Array<const Residue> residues) {
    memset_array(dst_atom_colors, 0xffffffff);
    for (const auto& res : residues) {
        const uint32 color = math::convert_color(color_from_hash(hash::crc32(CString(res.name))));
        memset_array(dst_atom_colors, color, res.atom_range);
    }
}
void color_atoms_residue_index(Array<uint32> dst_atom_colors, Array<const Residue> residues) {
    memset_array(dst_atom_colors, 0xffffffff);
    for (int64 i = 0; i < residues.count; i++) {
        const uint32 color = math::convert_color(color_from_hash(hash::crc32(i)));
        memset_array(dst_atom_colors, color, residues[i].atom_range);
    }
}
void color_atoms_chain_id(Array<uint32> dst_atom_colors, Array<const Chain> chains) {
    memset_array(dst_atom_colors, 0xffffffff);
    for (const auto& chain : chains) {
        const uint32 color = math::convert_color(color_from_hash(hash::crc32(CString(chain.id))));
        memset_array(dst_atom_colors, color, chain.atom_range);
    }
}
void color_atoms_chain_index(Array<uint32> dst_atom_colors, Array<const Chain> chains) {
    memset_array(dst_atom_colors, 0xffffffff);
    for (int64 i = 0; i < chains.count; i++) {
        const uint32 color = math::convert_color(color_from_hash(hash::crc32(i)));
        memset_array(dst_atom_colors, color, chains[i].atom_range);
    }
}

inline uint32 lerp_pixel(const Image& color_map, const vec2& coords) {
    const int x0 = math::clamp((int32)(coords.x * color_map.width), 0, color_map.width - 1);
    const int y0 = math::clamp((int32)(coords.y * color_map.height), 0, color_map.height - 1);
    const int x1 = math::min(x0 + 1, color_map.width - 1);
    const int y1 = math::min(y0 + 1, color_map.height - 1);
    const vec2 f = math::fract(coords);
    const vec4 cx0 = math::mix(math::convert_color(color_map.data[y0 * color_map.width + x0]), math::convert_color(color_map.data[y0 * color_map.width + x1]), f.x);
    const vec4 cx1 = math::mix(math::convert_color(color_map.data[y1 * color_map.width + x0]), math::convert_color(color_map.data[y1 * color_map.width + x1]), f.x);
    return math::convert_color(math::mix(cx0, cx1, f.y));
}

void color_atoms_backbone_angles(Array<uint32> dst_atom_colors, Array<const Residue> residues, Array<const BackboneSequence> bb_seq, Array<const vec2> bb_angles, const Image& color_map) {
    memset_array(dst_atom_colors, 0xffffffff);

    if (color_map.width == 0 || color_map.height == 0) return;
    const float one_over_two_pi = 1.f / (2.f * math::PI);

    for (const auto& seq : bb_seq) {
        if (seq.end - seq.beg < 2) continue;

        for (int64 i = seq.beg + 1; i < seq.end - 1; i++) {
            const vec2 coord = vec2(0, 1) + vec2(1, -1) * (bb_angles[i] * one_over_two_pi + 0.5f);
            const uint32 color = lerp_pixel(color_map, coord);
            memset_array(dst_atom_colors, color, residues[i].atom_range);
        }

        // Do first and last segment explicitly since it lacks adjacent [next] amino acid to properly compute phi and psi.
        {
            const auto dst_i = seq.beg;
            const auto src_i = math::min(seq.end - 1, dst_i + 1);
            const uint32 color = dst_atom_colors[residues[src_i].atom_range.beg];  // copy color from previous residue
            memset_array(dst_atom_colors, color, residues[dst_i].atom_range);
        }
        {
            const auto dst_i = math::max(0, seq.end - 1);
            const auto src_i = math::max(0, dst_i - 1);
            const uint32 color = dst_atom_colors[residues[src_i].atom_range.beg];  // copy color from previous residue
            memset_array(dst_atom_colors, color, residues[dst_i].atom_range);
        }
    }
}
