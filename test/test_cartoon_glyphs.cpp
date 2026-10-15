#include "utest.h"

#include <cartoon_glyphs.h>

#include <md_system.h>
#include <md_util.h>
#include <md_pdb.h>
#include <md_tpr.h>
#include <core/md_arena_allocator.h>

/* The cartoon of the lipids is built on the CPU (the points, which atoms each stands for, the tubes through them) and
 * placed on the GPU. Without a GL molecule to derive the points from, cartoon_glyphs_init builds the first part and
 * stops, which is what these check: the anatomy of the lipids turned into spheres and tubes. */

static uint32_t count_points(const CartoonGlyphs& g, CartoonGlyphKind kind, CartoonGlyphVariant variant, bool sphere) {
    uint32_t n = 0;
    for (uint32_t i = 0; i < g.num_points; ++i) {
        if (g.point_kind[i] == kind && g.point_variant[i] == variant && (g.point_radius[i] > 0.0f) == sphere) ++n;
    }
    return n;
}

static bool point_has_atom_named(const CartoonGlyphs& g, const md_system_t& sys, uint32_t p, const char* name) {
    for (uint32_t j = g.point_offset[p]; j < g.point_offset[p + 1]; ++j) {
        if (str_eq_cstr(md_atom_name(&sys.atom, g.point_atom[j]), name)) return true;
    }
    return false;
}

static bool point_has_element(const CartoonGlyphs& g, const md_system_t& sys, uint32_t p, md_atomic_number_t z) {
    for (uint32_t j = g.point_offset[p]; j < g.point_offset[p + 1]; ++j) {
        if (md_atom_atomic_number(&sys.atom, g.point_atom[j]) == z) return true;
    }
    return false;
}

UTEST(cartoon_glyphs, dppc_heads_split_at_the_phosphate) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, STR_LIT(VIAMD_TEST_DATA_DIR "/dppc64.pdb"), MD_PDB_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    ASSERT_EQ(64, (int)sys.lipid.count);

    CartoonGlyphs g = {};
    cartoon_glyphs_init(&g, sys, md_gl_mol_t{}, alloc);
    EXPECT_EQ(0u, g.gl_mol.id);

    // Every head splits: one sphere for it, or the phosphocholine and the glycerol with its esters
    EXPECT_EQ(64u, g.num_split_heads);
    EXPECT_EQ(64u, count_points(g, CartoonGlyph_LipidHead, CartoonGlyphVariant_OneHeadSphere, true));
    EXPECT_EQ(64u, count_points(g, CartoonGlyph_LipidHead, CartoonGlyphVariant_TwoHeadSpheres, true));
    EXPECT_EQ(64u, count_points(g, CartoonGlyph_LipidBackbone, CartoonGlyphVariant_TwoHeadSpheres, true));
    // The tubes start at the backbone in either case, at a point without a sphere
    EXPECT_EQ(64u, count_points(g, CartoonGlyph_LipidBackbone, CartoonGlyphVariant_Always, false));

    for (uint32_t i = 0; i < g.num_points; ++i) {
        const uint32_t n = g.point_offset[i + 1] - g.point_offset[i];
        if (g.point_variant[i] == CartoonGlyphVariant_OneHeadSphere) {
            EXPECT_EQ(20u, n);
        } else if (g.point_variant[i] == CartoonGlyphVariant_TwoHeadSpheres) {
            // Choline, its O, the P and the two free O of the phosphate; the O to the glycerol, the glycerol, the esters
            EXPECT_EQ(10u, n);
            const bool outer = g.point_kind[i] == CartoonGlyph_LipidHead;
            EXPECT_EQ(outer, point_has_element(g, sys, i, 15));  // P
            EXPECT_EQ(outer, point_has_element(g, sys, i, 7));   // N
        }
    }

    // Two tubes per lipid, each from a backbone root
    ASSERT_EQ(128u, g.num_tubes);
    for (uint32_t t = 0; t < g.num_tubes; ++t) {
        const uint32_t first = g.tube_point[g.tube_offset[t]];
        EXPECT_EQ((uint8_t)CartoonGlyph_LipidBackbone, g.point_kind[first]);
        EXPECT_EQ((uint8_t)CartoonGlyphVariant_Always, g.point_variant[first]);
        EXPECT_GE(g.tube_offset[t + 1] - g.tube_offset[t], 8u);
        for (uint32_t k = g.tube_offset[t] + 1; k < g.tube_offset[t + 1]; ++k) {
            EXPECT_EQ((uint8_t)CartoonGlyph_LipidTail, g.point_kind[g.tube_point[k]]);
            EXPECT_EQ(0.0f, g.point_radius[g.tube_point[k]]);
        }
    }

    // A tube picks as its control point, after the points
    EXPECT_EQ(g.num_points + g.num_tube_points, cartoon_glyphs_picking_count(g));
    EXPECT_EQ(g.tube_point[0], cartoon_glyphs_picked_point(g, g.num_points));
    EXPECT_EQ(UINT32_MAX, cartoon_glyphs_picked_point(g, cartoon_glyphs_picking_count(g)));

    cartoon_glyphs_free(&g);
    md_arena_allocator_destroy(alloc);
}

UTEST(cartoon_glyphs, martini_phospholipids_split_and_sterols_do_not) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &st, STR_LIT(VIAMD_TEST_DATA_DIR "/martini/martini3_membrane.tpr")));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL & ~MD_UTIL_INFER_BOND_BIT));

    CartoonGlyphs g = {};
    cartoon_glyphs_init(&g, sys, md_gl_mol_t{}, alloc);

    // POPC, POPE and POPS at their PO4: the head group bead and PO4, GL1 and GL2. CHOL has no phosphate.
    EXPECT_EQ(32u + 16u + 8u, g.num_split_heads);
    for (uint32_t i = 0; i < g.num_points; ++i) {
        if (g.point_variant[i] != CartoonGlyphVariant_TwoHeadSpheres) continue;
        const uint32_t n = g.point_offset[i + 1] - g.point_offset[i];
        EXPECT_EQ(2u, n);
        const bool outer = g.point_kind[i] == CartoonGlyph_LipidHead;
        EXPECT_EQ(outer, point_has_atom_named(g, sys, i, "PO4"));
        EXPECT_EQ(!outer, point_has_atom_named(g, sys, i, "GL1"));
    }
    // The heads of the sterols are one sphere, whatever a representation asks
    EXPECT_EQ(16u, count_points(g, CartoonGlyph_LipidHead, CartoonGlyphVariant_Always, true));

    cartoon_glyphs_free(&g);
    md_arena_allocator_destroy(alloc);
}
