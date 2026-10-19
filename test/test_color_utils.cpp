#include "utest.h"

#include <color_utils.h>

#include <md_system.h>
#include <md_util.h>
#include <md_pdb.h>
#include <core/md_arena_allocator.h>

#include <math.h>

// The colour space conversions are header only and are what the representation colouring and the
// Oklab colour optimizer are built on. Round trips are the useful property: a conversion pair that
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
    // The colour optimizer walks hue at fixed lightness and chroma, and hands the result straight
    // to the renderer, so the conversion has to come back inside the display gamut.
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
UTEST(viamd_color, comp_type) {
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
    ASSERT_EQ(3, (int)sys.component.count);

    uint32_t colors[9] = {};
    color_atoms_comp_type(colors, 9, sys);
    for (int a = 0; a < 4; ++a) EXPECT_EQ(md_monomer_color(MD_MONOMER_ALA), colors[a]);
    for (int a = 4; a < 8; ++a) EXPECT_EQ(md_monomer_color(MD_MONOMER_HIS), colors[a]);
    EXPECT_EQ(md_component_type_color(&sys.component.type, (size_t)md_component_type_idx(&sys.component, 2)) | 0xFF000000u, colors[8]);

    // A color of its own, opaque whatever alpha it was given
    const md_component_type_idx_t his = md_component_type_idx(&sys.component, 1);
    sys.component.type.color[his] = 0x00102030u;
    color_atoms_comp_type(colors, 9, sys);
    for (int a = 4; a < 8; ++a) EXPECT_EQ(0xFF102030u, colors[a]);

    md_arena_allocator_destroy(alloc);
}
