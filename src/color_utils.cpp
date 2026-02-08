#include "color_utils.h"

#include <core/md_hash.h>
#include <md_system.h>
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

#if 0
static float u32_to_hue(uint32_t u32) {
    const uint32_t mod = 24;
    // The u32 is assumed to be random between [0, UINT_MAX[
    return (float)(u32 % mod) / (float)mod;
}

static uint32_t u32_to_color(uint32_t u32) {
	const float chroma = 0.9f;
	const float luminance = 1.0f;
	const float hue = u32_to_hue(u32);
	const vec3_t lch = {hue, chroma, luminance};
	const vec3_t rgb  = oklch_to_srgb(lch);
	return convert_color({rgb.x, rgb.y, rgb.z, 1.0f});
}
#endif

// https://github.com/uekstrom/colorwheel
static const vec3_t wheel_colors[] = {
       {9.65246549e-01, 8.13827873e-01, 3.54995481e-02},
       {1.00000000e+00, 7.54400671e-01, 5.96369855e-08},
       {1.00000000e+00, 6.89975121e-01, 3.70155717e-07},
       {1.00000000e+00, 6.24055810e-01, 2.89885028e-07},
       {1.00000000e+00, 5.56388222e-01, 2.37011569e-02},
       {1.00000000e+00, 4.88310276e-01, 1.07342654e-01},
       {1.00000000e+00, 4.14139317e-01, 1.45399871e-01},
       {9.96733205e-01, 3.32898008e-01, 1.83204421e-01},
       {9.55860966e-01, 2.67691368e-01, 2.29196088e-01},
       {9.16812299e-01, 2.29007900e-01, 3.16986609e-01},
       {8.74115199e-01, 2.12576460e-01, 4.00747900e-01},
       {8.26373747e-01, 2.15992904e-01, 4.79921197e-01},
       {7.73442962e-01, 2.35149069e-01, 5.55273021e-01},
       {7.13797072e-01, 2.68304699e-01, 6.22764084e-01},
       {6.46284947e-01, 3.14257459e-01, 6.70179364e-01},
       {5.75826839e-01, 3.55350643e-01, 7.13045820e-01},
       {4.97736406e-01, 3.92601178e-01, 7.47614002e-01},
       {4.07784381e-01, 4.29442504e-01, 7.68234303e-01},
       {2.99426780e-01, 4.68016769e-01, 7.71232221e-01},
       {1.61515156e-01, 5.11441924e-01, 7.57959971e-01},
       {0.00000000e+00, 5.61134848e-01, 7.38655511e-01},
       {0.00000000e+00, 6.06702393e-01, 7.07289110e-01},
       {0.00000000e+00, 6.48146935e-01, 6.65693542e-01},
       {1.29865957e-02, 6.80804928e-01, 6.03421641e-01},
       {1.89447215e-01, 7.08723329e-01, 5.34681691e-01},
       {3.15812031e-01, 7.34579462e-01, 4.62348754e-01},
       {4.21580101e-01, 7.60995139e-01, 3.86226481e-01},
       {5.21205993e-01, 7.90186326e-01, 3.12525366e-01},
       {6.34918907e-01, 8.12483129e-01, 2.36750831e-01},
       {7.50249701e-01, 8.27893673e-01, 1.66742811e-01},
       {8.62790551e-01, 8.34388326e-01, 1.05388198e-01},
       {9.65246549e-01, 8.13827873e-01, 3.54995481e-02}
};

static inline vec3_t sample_colormap(const vec3_t* colormap, size_t colormap_size, float t) {
    const float scaled_t = t * (colormap_size - 1);
    const size_t idx0 = (size_t)floorf(scaled_t);
    const size_t idx1 = (size_t)ceilf(scaled_t);
    const float local_t = scaled_t - idx0;
    return vec3_lerp(colormap[idx0], colormap[idx1], local_t);
}

static inline uint32_t u32_to_color_sequential(uint32_t u32, uint32_t num_colors) {
    float t = (float)(u32 % num_colors) / (float)num_colors;
    vec3_t rgb = sample_colormap(wheel_colors, ARRAY_SIZE(wheel_colors), t);
    return convert_color({ rgb.x, rgb.y, rgb.z, 1.0f });
}

static inline uint32_t u32_to_color_qualitative(uint32_t u32) {
	// The golden ratio conjugate is used to ensure good distribution of colors even for non-random u32 values (e.g. sequential indices)
    const float t = fractf((float)u32 * 1.618033988749f);
    vec3_t rgb = sample_colormap(wheel_colors, ARRAY_SIZE(wheel_colors), t);
    return convert_color({ rgb.x, rgb.y, rgb.z, 1.0f });
}

static inline uint32_t u32_to_color(uint32_t u32) {
	return u32_to_color_qualitative(u32);
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

void color_atoms_cpk(uint32_t* colors, size_t count, const md_system_t& sys) {
    for (size_t i = 0; i < count; i++) {
        md_atomic_number_t z = md_atom_atomic_number(&sys.atom, i);
        colors[i] = z ? md_util_element_cpk_color(z) : 0xFFFFFFFFU;
    }
}

void color_atoms_type(uint32_t* colors, size_t count, const md_system_t& sys) {
    for (size_t i = 0; i < count; ++i) {
		md_atom_type_idx_t idx = md_atom_type_idx(&sys.atom, i);
        colors[i] = md_atom_type_color(&sys.atom.type, idx);
    }
}

void color_atoms_idx(uint32_t* colors, size_t count, const md_system_t&) {
    for (size_t i = 0; i < count; ++i) {
        colors[i] = u32_to_color((uint32_t)i);
    }
}

void color_atoms_comp_name(uint32_t* colors, size_t count, const md_system_t& sys) {
    set_colors(colors, count, 0xFFFFFFFFU);
    for (size_t i = 0; i < sys.comp.count; i++) {
        str_t str = sys.comp.name[i];
        const uint32_t u32 = md_hash32(str.ptr, str.len, 0);
        const uint32_t color = u32_to_color(u32);
        md_urange_t range = md_comp_atom_range(&sys.comp, i);
        set_colors(colors + range.beg, range.end - range.beg, color);
    }
}

void color_atoms_comp_seq_id(uint32_t* colors, size_t count, const md_system_t& sys) {
    set_colors(colors, count, 0xFFFFFFFFU);
    for (size_t i = 0; i < sys.comp.count; i++) {
        const uint32_t color = u32_to_color(sys.comp.seq_id[i]);
        md_urange_t range = md_comp_atom_range(&sys.comp, i);
        set_colors(colors + range.beg, range.end - range.beg, color);
    }
}
void color_atoms_comp_idx(uint32_t* colors, size_t count, const md_system_t& sys) {
    set_colors(colors, count, 0xFFFFFFFFU);
    for (size_t i = 0; i < sys.comp.count; i++) {
        const uint32_t color = u32_to_color((uint32_t)i);
        md_urange_t range = md_comp_atom_range(&sys.comp, i);
        set_colors(colors + range.beg, range.end - range.beg, color);
    }
}

void color_atoms_inst_id(uint32_t* colors, size_t count, const md_system_t& sys) {
    set_colors(colors, count, 0xFFFFFFFFU);
    for (size_t i = 0; i < sys.inst.count; i++) {
        str_t str = sys.inst.id[i];
        uint32_t u32 = 0;
        for (size_t j = 0; j < str.len; ++j) {
            u32 += str.ptr[j];
        }
        const uint32_t color = u32_to_color(u32);
        md_urange_t range = md_system_inst_atom_range(&sys, i);
        set_colors(colors + range.beg, range.end - range.beg, color);
    }
}
void color_atoms_inst_idx(uint32_t* colors, size_t count, const md_system_t& sys) {
    set_colors(colors, count, 0xFFFFFFFFU);
    for (size_t i = 0; i < sys.inst.count; i++) {
        const uint32_t color = u32_to_color((uint32_t)i);
        md_urange_t range = md_system_inst_atom_range(&sys, i);
        set_colors(colors + range.beg, range.end - range.beg, color);
    }
}

void color_atoms_sec_str(uint32_t* colors, size_t count, const md_system_t& sys) {
    const uint32_t color_unknown = 0x22222222;
    const uint32_t color_coil    = 0xDDDDDDDD;
    const uint32_t color_helix   = 0xFF22DD22;
    const uint32_t color_sheet   = 0xFFDD2222;

    if (sys.protein_backbone.segment.secondary_structure) {
        for (size_t i = 0; i < sys.protein_backbone.segment.count; i++) {
            md_secondary_structure_t ss = sys.protein_backbone.segment.secondary_structure[i];
            uint32_t color = color_unknown;
            switch (ss) {
                case MD_SECONDARY_STRUCTURE_HELIX_310:
                case MD_SECONDARY_STRUCTURE_HELIX_ALPHA:
                case MD_SECONDARY_STRUCTURE_HELIX_PI:
                    color = color_helix;
                    break;
                case MD_SECONDARY_STRUCTURE_BETA_SHEET:
                case MD_SECONDARY_STRUCTURE_BETA_BRIDGE:
                    color = color_sheet;
                    break;
                case MD_SECONDARY_STRUCTURE_COIL:
                case MD_SECONDARY_STRUCTURE_TURN:
                case MD_SECONDARY_STRUCTURE_BEND:
                case MD_SECONDARY_STRUCTURE_UNKNOWN:
                default:
                    color = color_coil;
                    break;
            }
            md_comp_idx_t res_idx = sys.protein_backbone.segment.comp_idx[i];
            md_urange_t range = md_comp_atom_range(&sys.comp, res_idx);
            set_colors(colors + range.beg, range.end - range.beg, color);
        }
    } else {
        set_colors(colors, count, color_unknown);
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

void scale_saturation(uint32_t* colors, const md_bitfield_t* mask, float scale) {
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

void scale_saturation(uint32_t* colors, size_t count, float scale) {
    for (size_t i = 0; i < count; ++i) {
        vec4_t rgba = convert_color(colors[i]);
        vec3_t hsv = rgb_to_hsv(vec3_from_vec4(rgba));
        hsv.y *= scale;
        rgba = vec4_from_vec3(hsv_to_rgb(hsv), rgba.w);
        colors[i] = convert_color(rgba);
    }
}
