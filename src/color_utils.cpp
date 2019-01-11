#include <core/array_types.h>
#include <core/hash.h>
#include <mol/molecule_structure.h>
#include "color_utils.h"

inline void memset32(void* buf, uint32_t val, int32_t count) {
    for (int32_t i = 0; i < count; i++) {
        ((uint32_t*)buf)[i] = val;
    }
}

void color_atoms_cpk(Array<uint32> dst_atom_colors, Array<const Element> elements) {
    for (int64 i = 0; i < dst_atom_colors.count; i++) {
        dst_atom_colors[i] = element::color(elements[i]);
    }
}

void color_atoms_residue_id(Array<uint32> dst_atom_colors, Array<const Residue> residues) {
    memset32(dst_atom_colors.data, 0xffffffff, dst_atom_colors.size());
    for (const auto& res : residues) {
        const uint32 color = math::convert_color(color_from_hash(hash::crc32(res.name.operator CString())));
        memset32(dst_atom_colors.data + res.atom_idx.beg, color, res.atom_idx.end - res.atom_idx.beg);
    }
}
void color_atoms_residue_index(Array<uint32> dst_atom_colors, Array<const Residue> residues) {
    memset32(dst_atom_colors.data, 0xffffffff, dst_atom_colors.size());
    for (int64 i = 0; i < residues.count; i++) {
        const uint32 color = math::convert_color(color_from_hash(hash::crc32(i)));
        memset32(dst_atom_colors.data + residues[i].atom_idx.beg, color, residues[i].atom_idx.end - residues[i].atom_idx.beg);
    }
}
void color_atoms_chain_id(Array<uint32> dst_atom_colors, Array<const Chain> chains, Array<const Residue> residues) {
    memset32(dst_atom_colors.data, 0xffffffff, dst_atom_colors.size());
    for (const auto& chain : chains) {
        const uint32 color = math::convert_color(color_from_hash(hash::crc32(chain.id.operator CString())));
        const auto beg_idx = residues[chain.res_idx.beg].atom_idx.beg;
        const auto end_idx = residues[chain.res_idx.end - 1].atom_idx.end;
        memset32(dst_atom_colors.data + beg_idx, color, end_idx - beg_idx);
    }
}
void color_atoms_chain_index(Array<uint32> dst_atom_colors, Array<const Chain> chains, Array<const Residue> residues) {
    memset32(dst_atom_colors.data, 0xffffffff, dst_atom_colors.size());
    for (int64 i = 0; i < chains.count; i++) {
        const uint32 color = math::convert_color(color_from_hash(hash::crc32(i)));
        const auto beg_idx = residues[chains[i].res_idx.beg].atom_idx.beg;
        const auto end_idx = residues[chains[i].res_idx.end - 1].atom_idx.end;
        memset32(dst_atom_colors.data + beg_idx, color, end_idx - beg_idx);
    }
}

inline uint32 fetch_pixel(const Image& color_map, const vec2& coords) {
    int x = math::clamp((int32)(coords.x * color_map.width), 0, math::max(0, color_map.width - 1));
    int y = math::clamp((int32)(coords.y * color_map.height), 0, math::max(0, color_map.height - 1));
    return color_map.data[y * color_map.width + x];
}

void color_atoms_backbone_angles(Array<uint32> dst_atom_colors, Array<const Residue> residues, Array<const BackboneSequence> bb_seq, Array<const vec2> bb_angles, const Image& color_map) {
    memset32(dst_atom_colors.data, 0xffffffff, dst_atom_colors.size());
    for (const auto& seq : bb_seq) {
        for (int64 i = seq.beg; i < seq.end; i++) {
            const float one_over_two_pi = 1.f / (2.f * math::PI);
            const vec2 coord = vec2(0, 1) + vec2(1, -1) * (bb_angles[i] * one_over_two_pi + 0.5f);
            const uint32 color = fetch_pixel(color_map, bb_angles[i] * one_over_two_pi + 0.5f);
            memset32(dst_atom_colors.data + residues[i].atom_idx.beg, color, residues[i].atom_idx.end - residues[i].atom_idx.beg);
        }
    }
}
