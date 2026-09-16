#include "utest.h"

#include <color_utils.h>

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
