#include "color_utils.h"

#include <core/md_hash.h>
#include <md_molecule.h>
#include <md_util.h>

static void set_colors(uint32_t* colors, size_t count, uint32_t color) {
    if (color == 0xFFFFFFFFU) {
        MEMSET(colors, 0xFF, sizeof(uint32_t) * count);
    } else {
        for (size_t i = 0; i < count; ++i) {
            colors[i] = color;
        }
    }
}

static float u32_to_hue(uint32_t u32) {
    const float scl = 1.618033988749f / 21.0f;
    const float hue = (float)u32 * scl;
    return fractf(hue);
}

static uint32_t u32_to_color(uint32_t u32) {
	const float chroma = 0.9f;
	const float luminance = 1.0f;
	const vec3_t hcl = {u32_to_hue(u32), chroma, luminance};
	const vec3_t rgb  = hcl_to_rgb(hcl);
	return convert_color({rgb.x, rgb.y, rgb.z, 1.0f});
}

void color_atoms_uniform(uint32_t* colors, size_t count, vec4_t color, const md_bitfield_t* mask) {
    if (mask) {
        const uint32_t u32_color = convert_color(color);
        md_bitfield_iter_t it = md_bitfield_iter_create(mask);
        while (md_bitfield_iter_next(&it)) {
            colors[md_bitfield_iter_idx(&it)] = u32_color;
        }
    } else {
        set_colors(colors, count, convert_color(color));
    }
}

const float chroma = 0.9f;
const float luminance = 1.0f;

void color_atoms_cpk(uint32_t* colors, size_t count, const md_molecule_t& mol) {
    for (size_t i = 0; i < count; i++) {
        colors[i] = mol.atom.element ? md_util_element_cpk_color(mol.atom.element[i]) : 0xFFFFFFFFU;
    }
}

void color_atoms_type(uint32_t* colors, size_t count, const md_molecule_t& mol) {
    for (size_t i = 0; i < count; ++i) {
        const uint32_t u32 = md_hash32(mol.atom.type[i].buf, mol.atom.type[i].len, 0);
        colors[i] = u32_to_color(u32);
    }
}

void color_atoms_idx(uint32_t* colors, size_t count, const md_molecule_t&) {
    for (size_t i = 0; i < count; ++i) {
        colors[i] = u32_to_color((uint32_t)i);
    }
}

void color_atoms_res_name(uint32_t* colors, size_t count, const md_molecule_t& mol) {
    set_colors(colors, count, 0xFFFFFFFFU);
    for (size_t i = 0; i < mol.residue.count; i++) {
        str_t str = mol.residue.name[i];
        const uint32_t u32 = md_hash32(str.ptr, str.len, 0);
        const uint32_t color = u32_to_color(u32);
        md_range_t range = md_residue_atom_range(mol.residue, i);
        set_colors(colors + range.beg, range.end - range.beg, color);
    }
}

void color_atoms_res_id(uint32_t* colors, size_t count, const md_molecule_t& mol) {
    set_colors(colors, count, 0xFFFFFFFFU);
    for (size_t i = 0; i < mol.residue.count; i++) {
        const uint32_t color = u32_to_color(mol.residue.id[i]);
        md_range_t range = md_residue_atom_range(mol.residue, i);
        set_colors(colors + range.beg, range.end - range.beg, color);
    }
}
void color_atoms_res_idx(uint32_t* colors, size_t count, const md_molecule_t& mol) {
    set_colors(colors, count, 0xFFFFFFFFU);
    for (size_t i = 0; i < mol.residue.count; i++) {
        const uint32_t color = u32_to_color((uint32_t)i);
        md_range_t range = md_residue_atom_range(mol.residue, i);
        set_colors(colors + range.beg, range.end - range.beg, color);
    }
}

void color_atoms_chain_id(uint32_t* colors, size_t count, const md_molecule_t& mol) {
    set_colors(colors, count, 0xFFFFFFFFU);
    for (size_t i = 0; i < mol.chain.count; i++) {
        str_t str = mol.chain.id[i];
        uint32_t u32 = 0;
        for (size_t j = 0; j < str.len; ++j) {
            u32 += str.ptr[j];
        }
        const uint32_t color = u32_to_color(u32);
        md_range_t range = md_chain_atom_range(mol.chain, i);
        set_colors(colors + range.beg, range.end - range.beg, color);
    }
}
void color_atoms_chain_idx(uint32_t* colors, size_t count, const md_molecule_t& mol) {
    set_colors(colors, count, 0xFFFFFFFFU);
    for (size_t i = 0; i < mol.chain.count; i++) {
        const uint32_t color = u32_to_color((uint32_t)i);
        md_range_t range = md_chain_atom_range(mol.chain, i);
        set_colors(colors + range.beg, range.end - range.beg, color);
    }
}

void color_atoms_sec_str(uint32_t* colors, size_t count, const md_molecule_t& mol) {
    const uint32_t color_unknown = 0x22222222;
    const uint32_t color_coil    = 0xDDDDDDDD;
    const uint32_t color_helix   = 0xFF22DD22;
    const uint32_t color_sheet   = 0xFFDD2222;

    set_colors(colors, count, color_unknown);
    if (mol.backbone.secondary_structure) {
        for (size_t i = 0; i < mol.backbone.count; i++) {
            const vec4_t w = convert_color((uint32_t)mol.backbone.secondary_structure[i]);
            const vec4_t rgba = w.x * convert_color(color_coil) + w.y * convert_color(color_helix) + w.z * convert_color(color_sheet);
            const uint32_t color = convert_color(rgba);
            md_residue_idx_t res_idx = mol.backbone.residue_idx[i];
            md_range_t range = md_residue_atom_range(mol.residue, res_idx);
            set_colors(colors + range.beg, range.end - range.beg, color);
        }
    }
}

void filter_colors(uint32_t* colors, size_t num_colors, const md_bitfield_t* mask) {
    for (size_t i = 0; i < num_colors; ++i) {
        colors[i] &= 0x00FFFFFFU;
    }

    md_bitfield_iter_t it = md_bitfield_iter_create(mask);
    while (md_bitfield_iter_next(&it)) {
        uint64_t idx = md_bitfield_iter_idx(&it);
        colors[idx] = 0xFF000000U | (colors[idx] & 0x00FFFFFFU);
    }
}

void desaturate_colors(uint32_t* colors, const md_bitfield_t* mask, float scale) {
    int64_t beg_bit = mask->beg_bit;
    int64_t end_bit = mask->end_bit;
    while ((beg_bit = md_bitfield_scan(mask, beg_bit, end_bit)) != 0) {
        int64_t i = beg_bit - 1;
        vec4_t rgba = convert_color(colors[i]);
        vec3_t hsv = rgb_to_hsv(vec3_from_vec4(rgba));
        hsv.y *= scale;
        rgba = vec4_from_vec3(hsv_to_rgb(hsv), rgba.w);
        colors[i] = convert_color(rgba);
    }
}
