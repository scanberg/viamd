#include "utest.h"

#include <color_utils.h>

#include <md_system.h>
#include <md_util.h>
#include <md_pdb.h>
#include <core/md_arena_allocator.h>

#include <math.h>
#include <initializer_list>     // The loops over braced lists: libc++ does not bring it in through the others

// The colour space conversions are header only and are what the representation colouring and the
// colour optimizer (optimize_colors) are built on. Round trips are the useful property: a conversion pair that
// does not return where it started will drift every time a colour passes through it.

namespace {

const vec3_t samples[] = {
    {0.00f, 0.00f, 0.00f},
    {1.00f, 1.00f, 1.00f},
    {0.50f, 0.50f, 0.50f},
    {1.00f, 0.00f, 0.00f},
    {0.00f, 1.00f, 0.00f},
    {0.00f, 0.00f, 1.00f},
    {0.20f, 0.60f, 0.35f},
    {0.85f, 0.42f, 0.10f},
    {0.05f, 0.10f, 0.90f},
    {0.73f, 0.73f, 0.20f},
};

}  // namespace

#define EXPECT_VEC3_NEAR(got, want, tol)                    \
    do {                                                    \
        const vec3_t _g = (got);                            \
        const vec3_t _w = (want);                           \
        EXPECT_NEAR((double)_w.x, (double)_g.x, (tol));     \
        EXPECT_NEAR((double)_w.y, (double)_g.y, (tol));     \
        EXPECT_NEAR((double)_w.z, (double)_g.z, (tol));     \
    } while (0)

UTEST(viamd_color, srgb_oklab_round_trip) {
    for (size_t i = 0; i < ARRAY_SIZE(samples); ++i) {
        const vec3_t back = oklab_to_srgb(srgb_to_oklab(samples[i]));
        EXPECT_VEC3_NEAR(back, samples[i], 1.0e-4);
    }
}

UTEST(viamd_color, oklab_anchors) {
    // Oklab is built so that white sits at L = 1 with no chroma and black at the origin. If these
    // move, every lightness band and contrast target expressed in Oklab moves with them.
    const vec3_t white = srgb_to_oklab(vec3_t{1, 1, 1});
    EXPECT_NEAR(1.0, (double)white.x, 1.0e-3);
    EXPECT_NEAR(0.0, (double)white.y, 1.0e-3);
    EXPECT_NEAR(0.0, (double)white.z, 1.0e-3);

    const vec3_t black = srgb_to_oklab(vec3_t{0, 0, 0});
    EXPECT_NEAR(0.0, (double)black.x, 1.0e-3);
    EXPECT_NEAR(0.0, (double)black.y, 1.0e-3);
    EXPECT_NEAR(0.0, (double)black.z, 1.0e-3);

    // Grey has no chroma either, and its lightness rises with the grey level
    const vec3_t mid = srgb_to_oklab(vec3_t{0.5f, 0.5f, 0.5f});
    EXPECT_NEAR(0.0, (double)mid.y, 1.0e-3);
    EXPECT_NEAR(0.0, (double)mid.z, 1.0e-3);
    EXPECT_GT(mid.x, black.x);
    EXPECT_LT(mid.x, white.x);
}

UTEST(viamd_color, linear_srgb_round_trip) {
    for (size_t i = 0; i < ARRAY_SIZE(samples); ++i) {
        const vec3_t back = linear_srgb_to_srgb(srgb_to_linear_srgb(samples[i]));
        EXPECT_VEC3_NEAR(back, samples[i], 1.0e-5);
    }
}

UTEST(viamd_color, xyz_lab_round_trip) {
    for (size_t i = 0; i < ARRAY_SIZE(samples); ++i) {
        const vec3_t xyz  = rgb_to_XYZ(samples[i]);
        const vec3_t back = XYZ_to_rgb(Lab_to_XYZ(XYZ_to_Lab(xyz)));
        EXPECT_VEC3_NEAR(back, samples[i], 1.0e-3);
    }
}

UTEST(viamd_color, rgb_hsv_round_trip) {
    for (size_t i = 0; i < ARRAY_SIZE(samples); ++i) {
        const vec3_t back = hsv_to_rgb(rgb_to_hsv(samples[i]));
        EXPECT_VEC3_NEAR(back, samples[i], 1.0e-5);
    }
}

UTEST(viamd_color, rgb_hcl_round_trip_is_approximate) {
    // Unlike the Oklab and XYZ paths above, which come back to within 1e-4, this HCL is a cheap
    // perceptual approximation and does not round trip exactly: the worst sample here is off by
    // about 0.015, and black does not return to black. That is the conversion's real behaviour
    // rather than a defect, so it is written down - if a colour path needs to survive repeated
    // conversion, Oklab is the one to use.
    for (size_t i = 0; i < ARRAY_SIZE(samples); ++i) {
        const vec3_t back = hcl_to_rgb(rgb_to_hcl(samples[i]));
        EXPECT_VEC3_NEAR(back, samples[i], 0.02);
    }
}

UTEST(viamd_color, oklch_hue_sweep_stays_in_gamut) {
    // A hue walked at fixed lightness and chroma is handed straight to the renderer, so the
    // conversion has to come back inside the display gamut.
    for (int i = 0; i < 64; ++i) {
        const float hue = (float)i / 64.0f;
        const vec3_t rgb = oklch_to_srgb(vec3_t{0.6f, 0.08f, hue});
        for (int c = 0; c < 3; ++c) {
            EXPECT_GE(rgb.elem[c], 0.0f);
            EXPECT_LE(rgb.elem[c], 1.0f);
        }
    }
}

// Res Type: every atom takes the color of the type of its residue, which is that of the standard residue the type is a
// form of (HIE is one of HIS) until it is given one of its own
UTEST(viamd_color, res_type) {
    static const char* pdb =
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O\n"
        "ATOM      5  N   HIE A   2       3.332   1.536   0.000  1.00  0.00           N\n"
        "ATOM      6  CA  HIE A   2       3.970   2.846   0.000  1.00  0.00           C\n"
        "ATOM      7  C   HIE A   2       5.486   2.705   0.000  1.00  0.00           C\n"
        "ATOM      8  O   HIE A   2       6.009   1.593   0.000  1.00  0.00           O\n"
        "HETATM    9  O   HOH A 101      20.000  20.000  20.000  1.00  0.00           O\n";

    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, str_from_cstr(pdb), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    ASSERT_EQ(3, (int)sys.residue.count);

    uint32_t colors[9] = {};
    color_atoms_res_type(colors, 9, sys);
    for (int a = 0; a < 4; ++a) EXPECT_EQ(md_monomer_color(MD_MONOMER_ALA), colors[a]);
    for (int a = 4; a < 8; ++a) EXPECT_EQ(md_monomer_color(MD_MONOMER_HIS), colors[a]);
    EXPECT_EQ(md_residue_type_color(&sys.residue.type, (size_t)md_residue_type_idx(&sys.residue, 2)) | 0xFF000000u, colors[8]);

    // A color of its own, opaque whatever alpha it was given
    const md_residue_type_idx_t his = md_residue_type_idx(&sys.residue, 1);
    sys.residue.type.color[his] = 0x00102030u;
    color_atoms_res_type(colors, 9, sys);
    for (int a = 4; a < 8; ++a) EXPECT_EQ(0xFF102030u, colors[a]);

    md_arena_allocator_destroy(alloc);
}

// The gamut mapping lowers the chroma of a color out of gamut, and keeps its lightness and hue
UTEST(viamd_color, oklch_mapped_keeps_lightness_and_hue) {
    for (int i = 0; i < 36; ++i) {
        const float h = (float)i / 36.0f;
        for (float L : {0.3f, 0.6f, 0.9f}) {
            const vec3_t rgb = oklch_to_srgb_mapped(L, 0.4f, h);    // Far out of gamut at every hue
            for (int c = 0; c < 3; ++c) {
                EXPECT_GE(rgb.elem[c], 0.0f);
                EXPECT_LE(rgb.elem[c], 1.0f);
            }
            const vec3_t lab = srgb_to_oklab(rgb);
            EXPECT_NEAR((double)L, (double)lab.x, 2.0e-3);
            const float C = sqrtf(lab.y * lab.y + lab.z * lab.z);
            if (C > 0.02f) {
                float dh = atan2f(lab.z, lab.y) / (2.0f * (float)PI) - h;
                dh -= roundf(dh);
                EXPECT_LT(fabsf(dh), 0.01f);
            }
        }
    }
}

namespace {

float oklab_distance(uint32_t a, uint32_t b, bool linear_input) {
    vec3_t ca = vec3_from_vec4(convert_color(a)), cb = vec3_from_vec4(convert_color(b));
    if (linear_input) { ca = linear_srgb_to_srgb(ca); cb = linear_srgb_to_srgb(cb); }
    const vec3_t d = srgb_to_oklab(ca) - srgb_to_oklab(cb);
    return sqrtf(d.x * d.x + d.y * d.y + d.z * d.z);
}

}  // namespace

// Twelve categories of one color come apart, and the same input gives the same colors
UTEST(viamd_color, optimize_colors_spreads_alike) {
    ColorOptimizeItem items[12];
    for (auto& it : items) it.color = 0xFF8080C0u;
    const ColorOptimizeResult r = optimize_colors(items, 12);
    EXPECT_EQ((size_t)12, r.num_moved);
    EXPECT_GT(r.min_distance, 0.08f);
    float min_d = 1.0f;
    for (int i = 0; i < 12; ++i) {
        EXPECT_EQ(0xFF000000u, items[i].color & 0xFF000000u);
        for (int j = i + 1; j < 12; ++j) min_d = MIN(min_d, oklab_distance(items[i].color, items[j].color, false));
    }
    EXPECT_NEAR((double)r.min_distance, (double)min_d, 1.0e-4);

    ColorOptimizeItem again[12];
    for (auto& it : again) it.color = 0xFF8080C0u;
    optimize_colors(again, 12);
    for (int i = 0; i < 12; ++i) EXPECT_EQ(items[i].color, again[i].color);
}

// Fixed ones keep their colors, the others keep away from them, and stay within the band as they are shown
UTEST(viamd_color, optimize_colors_anchors_and_band) {
    const uint32_t anchors[3] = {0xFF2020F0u, 0xFFF04646u, 0xFF838383u};   // CPK oxygen, nitrogen, carbon (ABGR)
    for (bool linear : {false, true}) {
        ColorOptimizeItem items[9];
        for (int i = 0; i < 3; ++i) { items[i].color = anchors[i]; items[i].fixed = true; }
        for (int i = 3; i < 9; ++i) items[i].color = anchors[i % 3];        // Start on top of the anchors
        ColorOptimizeDesc desc;
        desc.linear_input = linear;
        const ColorOptimizeResult r = optimize_colors(items, 9, desc);
        EXPECT_EQ((size_t)6, r.num_moved);
        for (int i = 0; i < 3; ++i) EXPECT_EQ(anchors[i], items[i].color);
        for (int i = 3; i < 9; ++i) {
            for (int j = 0; j < 9; ++j) {
                if (j != i) EXPECT_GT(oklab_distance(items[i].color, items[j].color, linear), 0.06f);
            }
            vec3_t shown = vec3_from_vec4(convert_color(items[i].color));
            if (linear) shown = linear_srgb_to_srgb(shown);
            const vec3_t lab = srgb_to_oklab(shown);
            EXPECT_GE(lab.x, desc.lightness_min - 0.01f);
            EXPECT_LE(lab.x, desc.lightness_max + 0.01f);
        }
    }
}

// Pulled colors which stand apart already keep their colors exactly, a grey among them too, and pulled items of one
// color stay of one color. The free ones keep away from them.
UTEST(viamd_color, optimize_colors_pulled_stay) {
    const uint32_t cpk[4] = {0xFF2020F0u, 0xFFF04646u, 0xFF838383u, 0xFF30C8FFu};   // Oxygen, nitrogen, carbon, sulfur (ABGR)
    ColorOptimizeItem items[12];
    for (int i = 0; i < 8; ++i) { items[i].color = cpk[i % 4]; items[i].pull = 0.5f; items[i].weight = 10.0f; }
    for (int i = 8; i < 12; ++i) items[i].color = cpk[2];                           // Free, start on top of carbon
    const ColorOptimizeResult r = optimize_colors(items, 12);
    for (int i = 0; i < 8; ++i) EXPECT_EQ(cpk[i % 4], items[i].color);
    EXPECT_EQ((size_t)4, r.num_moved);
    EXPECT_GT(r.min_distance, 0.06f);
}

// Pulled colors too close to tell apart move, but not far, and the harder the pull the less. Pulled items of one color
// move as one.
UTEST(viamd_color, optimize_colors_pulled_drift) {
    const uint32_t reds[2] = {0xFF2020F0u, 0xFF2828E0u};
    float moved[2] = {};
    for (int p = 0; p < 2; ++p) {
        ColorOptimizeItem items[3];
        for (int i = 0; i < 3; ++i) { items[i].color = reds[i & 1]; items[i].pull = p == 0 ? 0.3f : 1.5f; }
        ColorOptimizeDesc desc;
        optimize_colors(items, 3, desc);
        EXPECT_EQ(items[0].color, items[2].color);
        for (int i = 0; i < 2; ++i) {
            const float d = oklab_distance(items[i].color, reds[i], false);
            EXPECT_LE(d, desc.max_drift + 0.005f);
            moved[p] += d;
        }
        if (p == 0) EXPECT_GT(oklab_distance(items[0].color, items[1].color, false), 0.06f);
    }
    EXPECT_GT(moved[0], 0.04f);
    EXPECT_LT(moved[1], moved[0]);
}

// Rare categories of one look are told apart as much as common ones: gold, sulfur and phosphorus among far more carbon
// and hydrogen, as viamd shows them (linear input), come to the distance which is enough. The colors are those of
// mdlib's CPK table before it was tuned for separation (mdlib f7ee83e), where gold and sulfur were too close: the
// table's own colors no longer make the case.
UTEST(viamd_color, optimize_colors_rare_pulled_apart) {
    // Au, S, P, C, H, O, N (ABGR)
    const uint32_t colors[7] = {0xFF23D1FFu, 0xFF0AF0F0u, 0xFF00A5FFu, 0xFF838383u, 0xFFFFFFFFu, 0xFF0000F0u, 0xFFF04646u};
    const float count[7] = {144, 60, 12, 600, 900, 40, 10};
    ColorOptimizeItem items[7];
    uint32_t cpk[7];
    for (int i = 0; i < 7; ++i) {
        cpk[i] = items[i].color = colors[i];
        items[i].weight = count[i];
        items[i].pull = 0.5f;
    }
    ColorOptimizeDesc desc;
    desc.linear_input = true;
    EXPECT_LT(oklab_distance(cpk[0], cpk[1], true), 0.06f);     // Gold and sulfur, the closest
    const ColorOptimizeResult r = optimize_colors(items, 7, desc);
    EXPECT_GT(r.min_distance, desc.enough_distance - 0.005f);
    for (int i = 0; i < 7; ++i) EXPECT_LE(oklab_distance(items[i].color, cpk[i], true), desc.max_drift + 0.005f);
    for (int i = 3; i < 7; ++i) EXPECT_EQ(cpk[i], items[i].color);    // Carbon, hydrogen, oxygen and nitrogen stand apart already
}

// A category of most of the weight (the solvent) is lighter and calmer than the rare ones around it
UTEST(viamd_color, optimize_colors_calm_dominant) {
    ColorOptimizeItem items[8];
    for (int i = 0; i < 8; ++i) { items[i].color = 0xFF8080C0u; items[i].weight = 1.0f; }
    items[0].weight = 20.0f;
    optimize_colors(items, 8);
    const vec3_t dom = srgb_to_oklab(vec3_from_vec4(convert_color(items[0].color)));
    const float dom_C = sqrtf(dom.y * dom.y + dom.z * dom.z);
    ColorOptimizeDesc d;
    EXPECT_GE(dom.x, d.lightness_min + 0.8f * (d.lightness_max - d.lightness_min) * 0.99f);
    EXPECT_LE(dom_C, d.chroma_max * 0.4f + 0.005f);
}
