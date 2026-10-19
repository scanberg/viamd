#include "color_utils.h"

#include <core/md_allocator.h>
#include <core/md_hash.h>
#include <md_system.h>
#include <md_util.h>

#include <math.h>
#include <vector>

static void set_colors(uint32_t* colors, size_t count, uint32_t color) {
    if (color == 0xFFFFFFFFU) {
        MEMSET(colors, 0xFF, sizeof(uint32_t) * count);
    } else {
        for (size_t i = 0; i < count; ++i) {
            colors[i] = color;
        }
    }
}


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
    // The multiplication is done in integer arithmetic (Fibonacci hashing): 0x9E3779B9 = 2^32 / phi, and the wrap-around modulo 2^32
    // is exactly the fractional part of u32 * phi. Doing this in float breaks for u32 >= 2^24 (e.g. hashes), where the product has no
    // fractional bits left and every value maps to t = 0 (the same color).
    const uint32_t h = u32 * 0x9E3779B9u;
    const float t = (float)(h >> 8) * (1.0f / 16777216.0f);
    vec3_t rgb = sample_colormap(wheel_colors, ARRAY_SIZE(wheel_colors), t);
    return convert_color({ rgb.x, rgb.y, rgb.z, 1.0f });
}

static inline uint32_t u32_to_color(uint32_t u32) {
	return u32_to_color_qualitative(u32);
}

void color_atoms_uniform(uint32_t* colors, size_t count, uint32_t color) {
    set_colors(colors, count, color);
}

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
    for (size_t i = 0; i < sys.component.count; i++) {
        str_t str = sys.component.name[i];
        const uint32_t u32 = md_hash32(str.ptr, str.len, 0);
        const uint32_t color = u32_to_color(u32);
        md_urange_t range = md_component_atom_range(&sys.component, i);
        set_colors(colors + range.beg, range.end - range.beg, color);
    }
}

void color_atoms_comp_type(uint32_t* colors, size_t count, const md_system_t& sys) {
    set_colors(colors, count, 0xFFFFFFFFU);
    for (size_t i = 0; i < sys.component.count; i++) {
        const md_component_type_idx_t type_idx = md_component_type_idx(&sys.component, i);
        if (type_idx < 0) continue;
        const uint32_t color = md_component_type_color(&sys.component.type, (size_t)type_idx) | 0xFF000000U;
        md_urange_t range = md_component_atom_range(&sys.component, i);
        range.end = MIN(range.end, (uint32_t)count);
        if (range.beg < range.end) set_colors(colors + range.beg, range.end - range.beg, color);
    }
}

void color_atoms_comp_seq_id(uint32_t* colors, size_t count, const md_system_t& sys) {
    set_colors(colors, count, 0xFFFFFFFFU);
    for (size_t i = 0; i < sys.component.count; i++) {
        const uint32_t color = u32_to_color(sys.component.seq_id[i]);
        md_urange_t range = md_component_atom_range(&sys.component, i);
        set_colors(colors + range.beg, range.end - range.beg, color);
    }
}
void color_atoms_comp_idx(uint32_t* colors, size_t count, const md_system_t& sys) {
    set_colors(colors, count, 0xFFFFFFFFU);
    for (size_t i = 0; i < sys.component.count; i++) {
        const uint32_t color = u32_to_color((uint32_t)i);
        md_urange_t range = md_component_atom_range(&sys.component, i);
        set_colors(colors + range.beg, range.end - range.beg, color);
    }
}

void color_atoms_inst_id(uint32_t* colors, size_t count, const md_system_t& sys) {
    set_colors(colors, count, 0xFFFFFFFFU);
    for (size_t i = 0; i < sys.instance.count; i++) {
        str_t str = sys.instance.id[i];
        uint32_t u32 = 0;
        for (size_t j = 0; j < str.len; ++j) {
            u32 += str.ptr[j];
        }
        const uint32_t color = u32_to_color(u32);
        md_urange_t range = md_system_instance_atom_range(&sys, i);
        set_colors(colors + range.beg, range.end - range.beg, color);
    }
}
void color_atoms_inst_idx(uint32_t* colors, size_t count, const md_system_t& sys) {
    set_colors(colors, count, 0xFFFFFFFFU);
    for (size_t i = 0; i < sys.instance.count; i++) {
        const uint32_t color = u32_to_color((uint32_t)i);
        md_urange_t range = md_system_instance_atom_range(&sys, i);
        set_colors(colors + range.beg, range.end - range.beg, color);
    }
}

void color_atoms_secondary_structure(uint32_t* colors, size_t count, const md_system_t& sys, const md_secondary_structure_t* secondary_structure, const SecondaryStructurePalette& palette) {
    const uint32_t color_unknown = palette.unknown;
    const uint32_t color_coil    = palette.coil;
    const uint32_t color_helix   = palette.helix;
    const uint32_t color_sheet   = palette.sheet;

    if (secondary_structure) {
        for (size_t i = 0; i < sys.protein_backbone.range.count; ++i) {
            uint32_t range_beg = sys.protein_backbone.range.offset[i];
            uint32_t range_end = sys.protein_backbone.range.offset[i+1];
            uint32_t range_ext = range_end - range_beg;

            md_temp_scope_t temp_scope = md_temp_begin();
            defer { md_temp_end(temp_scope); };

            uint32_t* seg_colors = md_temp_alloc_array(temp_scope, uint32_t, range_ext);
            for (size_t j = 0; j < range_ext; ++j) {
                md_secondary_structure_t ss = secondary_structure[range_beg + j];
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
                seg_colors[j] = color;
            }

            // Fill isolated coils within sheets to sheet color for better visualization
            for (size_t j = 0; j < range_ext; ++j) {
                uint32_t color = seg_colors[j];
#if 1
                if (j > 0 && j + 1 < range_ext) {
                    if (color == color_coil && seg_colors[j-1] == color_sheet && seg_colors[j+1] == color_sheet) {
                        color = color_sheet;
                    }
                }
#endif
                md_component_idx_t comp_idx = sys.protein_backbone.segment.comp_idx[range_beg + j];
                md_urange_t range = md_component_atom_range(&sys.component, comp_idx);
                set_colors(colors + range.beg, range.end - range.beg, color);
            }
        }
    } else {
        set_colors(colors, count, color_unknown);
    }
}

void filter_colors(uint32_t* colors, size_t num_colors, md_bitset_t mask) {
    for (size_t i = 0; i < num_colors; ++i) {
        colors[i] &= 0x00FFFFFFU;
    }

    md_bitset_iter_t it = md_bitset_iter(mask);
    while (md_bitset_iter_next(&it)) {
        uint64_t idx = it.idx;
        colors[idx] = 0xFF000000U | (colors[idx] & 0x00FFFFFFU);
    }
}

void scale_saturation(uint32_t* colors, md_bitset_t mask, float scale) {
    md_bitset_iter_t it = md_bitset_iter(mask);
    while (md_bitset_iter_next(&it)) {
        const uint32_t i = it.idx;
        colors[i] = scale_saturation(colors[i], scale);
    }
}

void scale_saturation(uint32_t* colors, size_t count, float scale) {
    for (size_t i = 0; i < count; ++i) {
        colors[i] = scale_saturation(colors[i], scale);
    }
}

void tint_colors(uint32_t* colors, size_t count, uint32_t tint_color, float tint_strength, float saturation) {
    const float t = CLAMP(tint_strength, 0.0f, 1.0f);
    const float s = CLAMP(saturation, 0.0f, 1.0f);

    vec4_t tint_rgba = convert_color(tint_color);
    vec3_t tint_lab = linear_srgb_to_oklab(vec3_from_vec4(tint_rgba));

    const float L_coeff  = 0.1f * t;
    const float ab_coeff = 1.0f * t;

    for (size_t i = 0; i < count; ++i) {
        vec4_t rgba = convert_color(colors[i]);
        vec3_t lab = linear_srgb_to_oklab(vec3_from_vec4(rgba));

        // First apply tint.
        lab.x = lab.x * (1.0f - L_coeff)  + tint_lab.x * L_coeff;
        lab.y = lab.y * (1.0f - ab_coeff) + tint_lab.y * ab_coeff;
        lab.z = lab.z * (1.0f - ab_coeff) + tint_lab.z * ab_coeff;

        // Then apply saturation independently by scaling chroma in OKLab.
        lab.y *= s;
        lab.z *= s;

        vec3_t out_rgb = oklab_to_linear_srgb(lab);
        colors[i] = convert_color(vec4_from_vec3(out_rgb, rgba.w));
    }
}

void tint_colors(vec4_t* colors, size_t count, vec4_t tint_color, float tint_strength, float saturation) {
    const float t = CLAMP(tint_strength, 0.0f, 1.0f);
    const float s = CLAMP(saturation, 0.0f, 1.0f);
    vec3_t tint_lab = linear_srgb_to_oklab(vec3_from_vec4(tint_color));
    const float L_coeff = 0.1f * t;
    const float ab_coeff = 1.0f * t;
    for (size_t i = 0; i < count; ++i) {
        vec4_t rgba = colors[i];
        vec3_t lab = linear_srgb_to_oklab(vec3_from_vec4(rgba));
        // First apply tint.
        lab.x = lab.x * (1.0f - L_coeff) + tint_lab.x * L_coeff;
        lab.y = lab.y * (1.0f - ab_coeff) + tint_lab.y * ab_coeff;
        lab.z = lab.z * (1.0f - ab_coeff) + tint_lab.z * ab_coeff;
        // Then apply saturation independently by scaling chroma in OKLab.
        lab.y *= s;
        lab.z *= s;
        vec3_t out_rgb = oklab_to_linear_srgb(lab);
        colors[i] = vec4_from_vec3(out_rgb, rgba.w);
    }
}
// ## Optimized categorical colors

vec3_t oklch_to_srgb_mapped(float L, float C, float h) {
    L = CLAMP(L, 0.0f, 1.0f);
    C = MAX(C, 0.0f);
    const float ang = h * 2.0f * (float)PI;
    const float ca = cosf(ang);
    const float sa = sinf(ang);
    const auto linear = [&](float c) { return oklab_to_linear_srgb(vec3_t{L, c * ca, c * sa}); };
    const auto in_gamut = [](vec3_t v) {
        const float eps = 1.0e-4f;
        return v.x >= -eps && v.y >= -eps && v.z >= -eps && v.x <= 1.0f + eps && v.y <= 1.0f + eps && v.z <= 1.0f + eps;
    };

    vec3_t lin = linear(C);
    if (!in_gamut(lin)) {
        // The grey of the lightness is in gamut, the chroma in between is found by bisection
        float lo = 0.0f, hi = C;
        lin = linear(0.0f);
        for (int i = 0; i < 24; ++i) {
            const float mid = 0.5f * (lo + hi);
            const vec3_t m = linear(mid);
            if (in_gamut(m)) { lo = mid; lin = m; } else { hi = mid; }
        }
    }
    return linear_srgb_to_srgb(lin);
}

namespace {

// Machado, Oliveira and Fernandes 2009: deuteranopia of severity 1, in linear sRGB
constexpr float deutan_matrix[3][3] = {
    { 0.367322f, 0.860646f, -0.227968f},
    { 0.280085f, 0.672501f,  0.047413f},
    {-0.011820f, 0.042940f,  0.968881f},
};

struct OptColor {
    float L = 0, C = 0, h = 0;  // The color shown, OKLCH (h in turns)
    vec3_t srgb = {};           // The color shown
    vec3_t lab = {};            // Its OKLab
    vec3_t cvd = {};            // The OKLab of what a deuteranope sees
};

inline vec3_t shown_from_stored(uint32_t stored, bool linear_input) {
    const vec4_t c = convert_color(stored);
    const vec3_t rgb = {c.x, c.y, c.z};
    // Taken as a linear reflectance, the value is shown as its encoding
    return linear_input ? linear_srgb_to_srgb(rgb) : rgb;
}

inline uint32_t stored_from_shown(vec3_t srgb, bool linear_input) {
    const vec3_t rgb = linear_input ? srgb_to_linear_srgb(srgb) : srgb;
    return convert_color(vec4_t{rgb.x, rgb.y, rgb.z, 1.0f});
}

void set_shown(OptColor& c, vec3_t srgb) {
    c.srgb = srgb;
    c.lab = srgb_to_oklab(srgb);
    const vec3_t lin = srgb_to_linear_srgb(srgb);
    vec3_t d;
    d.x = deutan_matrix[0][0] * lin.x + deutan_matrix[0][1] * lin.y + deutan_matrix[0][2] * lin.z;
    d.y = deutan_matrix[1][0] * lin.x + deutan_matrix[1][1] * lin.y + deutan_matrix[1][2] * lin.z;
    d.z = deutan_matrix[2][0] * lin.x + deutan_matrix[2][1] * lin.y + deutan_matrix[2][2] * lin.z;
    c.cvd = linear_srgb_to_oklab(vec3_clamp_f(d, 0.0f, 1.0f));
}

void set_lch(OptColor& c, float L, float C, float h) {
    c.L = L; c.C = C; c.h = h - floorf(h);
    set_shown(c, oklch_to_srgb_mapped(c.L, c.C, c.h));
}

inline float distance(vec3_t a, vec3_t b) {
    const vec3_t d = a - b;
    return sqrtf(d.x * d.x + d.y * d.y + d.z * d.z);
}

// A small deterministic generator, the same colors for the same input on every platform
struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed * 0x9E3779B97F4A7C15ull + 0x632BE59BD9B4E019ull) {}
    uint64_t next() { s ^= s >> 12; s ^= s << 25; s ^= s >> 27; return s * 0x2545F4914F6CDD1Dull; }
    float uniform() { return (float)(next() >> 40) * (1.0f / 16777216.0f); }
    float normal() {
        const float u1 = MAX(uniform(), 1.0e-7f);
        const float u2 = uniform();
        return sqrtf(-2.0f * logf(u1)) * cosf(2.0f * (float)PI * u2);
    }
};

}  // namespace

ColorOptimizeResult optimize_colors(ColorOptimizeItem* items, size_t count, const ColorOptimizeDesc& desc) {
    ColorOptimizeResult result = {};
    if (!items || count == 0) return result;

    const float Lmin = CLAMP(MIN(desc.lightness_min, desc.lightness_max), 0.0f, 1.0f);
    const float Lmax = CLAMP(MAX(desc.lightness_min, desc.lightness_max), 0.0f, 1.0f);
    const float Cmin = MAX(0.0f, MIN(desc.chroma_min, desc.chroma_max));
    const float Cmax = MAX(desc.chroma_min, desc.chroma_max);
    const float grey_chroma = 0.03f;    // Below this a color counts as a grey, and one which keeps its hue stays grey
    const float hue_slack   = 0.02f;    // How far (in turns) one which keeps its hue may turn
    const float tau         = 0.01f;    // The softness of the smallest distance which is pushed up
    const float wc          = MAX(desc.cvd_weight, 0.0f);

    std::vector<OptColor> c(count);
    std::vector<float> h0(count, 0.0f);
    std::vector<uint8_t> grey(count, 0);
    float wmax = 0.0f;
    for (size_t i = 0; i < count; ++i) wmax = MAX(wmax, items[i].weight);
    if (wmax <= 0.0f) wmax = 1.0f;
    std::vector<float> wq(count);   // The weights, normalized and flattened: a rare item may come closer, not much closer
    for (size_t i = 0; i < count; ++i) wq[i] = powf(CLAMP(items[i].weight / wmax, 1.0e-6f, 1.0f), 0.25f);

    // The band of each item: one of a large share of the weight is held lighter and calmer
    std::vector<float> lo_L(count, Lmin), hi_C(count, Cmax);
    if (desc.calm_dominant) {
        double wsum = 0.0;
        for (size_t i = 0; i < count; ++i) wsum += MAX(items[i].weight, 0.0f);
        for (size_t i = 0; i < count && wsum > 0.0; ++i) {
            const float share = (float)(MAX(items[i].weight, 0.0f) / wsum);
            const float dominance = MIN(1.0f, 2.0f * share);   // From none to full as the share goes to a half
            lo_L[i] = Lmin + (Lmax - Lmin) * 0.8f * dominance;
            hi_C[i] = MAX(Cmin, Cmax * (1.0f - 0.6f * dominance));
        }
    }

    std::vector<size_t> moving;
    for (size_t i = 0; i < count; ++i) {
        const vec3_t shown = shown_from_stored(items[i].color, desc.linear_input);
        set_shown(c[i], shown);
        const vec3_t lab = c[i].lab;
        c[i].L = lab.x;
        c[i].C = sqrtf(lab.y * lab.y + lab.z * lab.z);
        c[i].h = atan2f(lab.z, lab.y) / (2.0f * (float)PI);
        c[i].h -= floorf(c[i].h);
        if (items[i].fixed) continue;

        moving.push_back(i);
        h0[i]   = c[i].h;
        grey[i] = items[i].keep_hue && c[i].C < grey_chroma;
        float L = CLAMP(c[i].L, lo_L[i], Lmax);
        float C = items[i].keep_hue ? (grey[i] ? c[i].C : CLAMP(c[i].C, grey_chroma, MAX(hi_C[i], c[i].C))) : CLAMP(c[i].C, Cmin, hi_C[i]);
        float h = c[i].h;
        // Items which start out alike (types of one hashed color, or of one element) are spread first
        for (size_t j = 0; j < i; ++j) {
            if (distance(lab, c[j].lab) < 0.02f) {
                if (items[i].keep_hue) L = CLAMP(L + ((moving.size() & 1) ? 0.06f : -0.06f), lo_L[i], Lmax);
                else                   h = h + 0.618034f * (float)moving.size();
                break;
            }
        }
        set_lch(c[i], L, C, h);
    }
    result.num_moved = moving.size();

    // The pairs which can change are those with an item which moves. Each adds exp(-e / tau) to what is minimized, e
    // being its distance (normal vision and deuteranopia) over how much the pair matters: a soft minimum.
    const auto pair_term = [&](const OptColor& a, const OptColor& b, size_t i, size_t j) -> double {
        const float d = (distance(a.lab, b.lab) + wc * distance(a.cvd, b.cvd)) / (1.0f + wc);
        const float importance = 0.5f + 0.5f * wq[i] * wq[j];
        return exp(-(double)(d / importance) / (double)tau);
    };

    if (!moving.empty()) {
        std::vector<double> term(count * count, 0.0);
        for (size_t i = 0; i < count; ++i) {
            for (size_t j = i + 1; j < count; ++j) {
                if (items[i].fixed && items[j].fixed) continue;
                term[i * count + j] = term[j * count + i] = pair_term(c[i], c[j], i, j);
            }
        }

        Rng rng(desc.seed);
        std::vector<double> trial(count, 0.0);
        const uint64_t total = MIN((uint64_t)desc.moves_per_item * moving.size(), (uint64_t)120000);
        for (uint64_t m = 0; m < total; ++m) {
            const size_t k = moving[rng.next() % moving.size()];
            const float a = 1.0f - 0.85f * (float)m / (float)total;    // Steps shrink as it settles

            // Drawn before they are clamped: CLAMP evaluates its arguments more than once
            const float nL = rng.normal(), nC = rng.normal(), nh = rng.normal();
            const float L0 = c[k].L + nL * 0.04f * a;
            float L = CLAMP(L0, lo_L[k], Lmax);
            float C, h;
            if (items[k].keep_hue) {
                const float C0 = c[k].C + nC * (grey[k] ? 0.01f : 0.025f) * a;
                const float Chi = grey[k] ? grey_chroma : MAX(hi_C[k], c[k].C);
                C = CLAMP(C0, grey[k] ? 0.0f : grey_chroma, Chi);
                float dh = c[k].h + nh * 0.01f * a - h0[k];
                dh -= roundf(dh);
                h = h0[k] + CLAMP(dh, -hue_slack, hue_slack);
            } else {
                const float C0 = c[k].C + nC * 0.025f * a;
                C = CLAMP(C0, Cmin, hi_C[k]);
                h = c[k].h + nh * 0.06f * a;
            }
            OptColor cand = c[k];
            set_lch(cand, L, C, h);

            double before = 0.0, after = 0.0;
            for (size_t j = 0; j < count; ++j) {
                if (j == k) continue;
                before += term[k * count + j];
                trial[j] = pair_term(cand, c[j], k, j);
                after += trial[j];
            }
            if (after <= before) {
                c[k] = cand;
                for (size_t j = 0; j < count; ++j) {
                    if (j == k) continue;
                    term[k * count + j] = term[j * count + k] = trial[j];
                }
            }
        }

        for (size_t i : moving) {
            items[i].color = stored_from_shown(c[i].srgb, desc.linear_input);
            set_shown(c[i], shown_from_stored(items[i].color, desc.linear_input));  // What is reported is what is stored
        }
    }

    float min_d = 1.0e9f;
    for (size_t i = 0; i < count; ++i) {
        for (size_t j = i + 1; j < count; ++j) {
            if (items[i].fixed && items[j].fixed) continue;
            min_d = MIN(min_d, distance(c[i].lab, c[j].lab));
        }
    }
    result.min_distance = min_d < 1.0e9f ? min_d : 0.0f;
    return result;
}
