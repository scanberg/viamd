#include "utest.h"

#include <cartoon_glyphs.h>

#include <md_system.h>
#include <md_util.h>
#include <md_pdb.h>
#include <md_tpr.h>
#include <md_gro.h>
#include <md_xyz.h>
#include <core/md_arena_allocator.h>

#include <float.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>

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
    // A head sphere per lipid, and the two of each split head
    EXPECT_EQ(64u * 3u, g.num_spheres);

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

// The optimal LINCS cholesterol of Martini 2 carries the mass of its rings on dummies which do not interact (VOH V1 V2
// V3), its ring beads virtual sites of them, so the spine of its skeleton runs over the dummies to the side chain. It is
// drawn as the cholesterol of 8 beads is: the slab of its ring beads, then the side chain C1 C2, and no control point on
// a dummy. Here one of each, of the same beads where they are (the dummies on the beads they stand in for, as a
// structure built from the cholesterol of 8 beads has them).
UTEST(cartoon_glyphs, martini_cholesterol_of_dummies_drawn_as_of_beads) {
    static const char gro[] =
        "CHOL of 12 and of 8\n"
        "   20\n"
        "    1CHOL   VOH    1   9.702  18.783  12.773\n"
        "    1CHOL    V2    2   9.932  18.814  13.207\n"
        "    1CHOL    V3    3   9.919  18.555  13.286\n"
        "    1CHOL    V1    4  10.081  18.381  13.768\n"
        "    1CHOL   ROH    5   9.702  18.783  12.773\n"
        "    1CHOL    R1    6   9.905  18.703  12.954\n"
        "    1CHOL    R2    7   9.932  18.814  13.207\n"
        "    1CHOL    R3    8   9.919  18.555  13.286\n"
        "    1CHOL    R4    9  10.034  18.695  13.515\n"
        "    1CHOL    R5   10  10.117  18.521  13.415\n"
        "    1CHOL    C1   11  10.081  18.381  13.768\n"
        "    1CHOL    C2   12   9.921  18.432  14.041\n"
        "    2CHOL   ROH   13   6.702  18.783  12.773\n"
        "    2CHOL    R1   14   6.905  18.703  12.954\n"
        "    2CHOL    R2   15   6.932  18.814  13.207\n"
        "    2CHOL    R3   16   6.919  18.555  13.286\n"
        "    2CHOL    R4   17   7.034  18.695  13.515\n"
        "    2CHOL    R5   18   7.117  18.521  13.415\n"
        "    2CHOL    C1   19   7.081  18.381  13.768\n"
        "    2CHOL    C2   20   6.921  18.432  14.041\n"
        "  20.00000  30.00000  30.00000\n";
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &st, str_from_cstr(gro)));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    ASSERT_EQ(2u, (uint32_t)sys.lipid.count);

    CartoonGlyphs g = {};
    cartoon_glyphs_init(&g, sys, md_gl_mol_t{}, alloc);
    ASSERT_EQ(2u, g.num_tubes);
    // Alike: the head, the sections of the slab, C1 and C2
    EXPECT_EQ(g.tube_offset[2] - g.tube_offset[1], g.tube_offset[1] - g.tube_offset[0]);
    for (uint32_t t = 0; t < 2; ++t) {
        const uint32_t beg = g.tube_offset[t], end = g.tube_offset[t + 1];
        ASSERT_GE(end - beg, 3u);
        EXPECT_TRUE(point_has_atom_named(g, sys, g.tube_point[end - 1], "C2"));
        EXPECT_TRUE(point_has_atom_named(g, sys, g.tube_point[end - 2], "C1"));
        for (uint32_t k = beg; k < end; ++k) {
            const uint32_t p = g.tube_point[k];
            const bool dummy = g.point_offset[p + 1] - g.point_offset[p] == 1 && md_atom_name(&sys.atom, g.point_atom[g.point_offset[p]]).ptr[0] == 'V';
            EXPECT_FALSE(dummy);
        }
    }

    cartoon_glyphs_free(&g);
    md_arena_allocator_destroy(alloc);
}

UTEST(cartoon_glyphs, nucleic_acids_are_a_band_with_sugar_and_base_slabs) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, STR_LIT(VIAMD_TEST_DATA_DIR "/nucl-dna.gro")));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    const md_nucleic_backbone_data_t& bb = sys.nucleic_backbone;
    ASSERT_GT(bb.segment.count, (size_t)100);

    // A sugar and a ring for a pyrimidine, two for a purine, and the stick between the sugar and the base
    uint32_t expected_slabs = 0;
    uint32_t nucleotides = 0;
    for (size_t s = 0; s < bb.segment.count; ++s) {
        if (bb.segment.atoms[s].c4 < 0) continue;
        const str_t name = md_component_name(&sys.component, bb.segment.comp_idx[s]);
        const bool purine = str_eq_cstr(name, "DA") || str_eq_cstr(name, "DG");
        expected_slabs += purine ? 4 : 3;
        nucleotides += 1;
    }

    CartoonGlyphs g = {};
    cartoon_glyphs_init(&g, sys, md_gl_mol_t{}, alloc);
    EXPECT_EQ(expected_slabs, g.num_slabs);
    // A control point per nucleotide, which faces its base and which the band passes through
    EXPECT_EQ(nucleotides, g.num_tube_points);
    // No spheres: the points are control points and corners
    EXPECT_EQ(0u, g.num_spheres);
    for (uint32_t t = 0; t < g.num_tubes; ++t) {
        EXPECT_EQ(1, (int)g.tube_through[t]);
    }
    for (uint32_t i = 0; i < g.num_tube_points; ++i) {
        EXPECT_GT(g.tube_flat[i], 0.0f);
        ASSERT_LT(g.tube_face[i], g.num_points);
        EXPECT_EQ((uint8_t)CartoonGlyph_NucleicBase, g.point_kind[g.tube_face[i]]);
    }
    uint32_t sticks = 0;
    for (uint32_t s = 0; s < g.num_slabs; ++s) {
        const uint32_t n = g.slab_offset[s + 1] - g.slab_offset[s];
        EXPECT_TRUE(n == 2 || n == 5 || n == 6);
        // Corners without spheres
        for (uint32_t k = g.slab_offset[s]; k < g.slab_offset[s + 1]; ++k) {
            EXPECT_EQ(0.0f, g.point_radius[g.slab_point[k]]);
        }
        if (n == 2) {
            // A stick: a capsule from the C1' of the sugar to the glycosidic N of the base, each a corner of its ring
            const uint32_t c1 = g.slab_point[g.slab_offset[s] + 0];
            const uint32_t n9 = g.slab_point[g.slab_offset[s] + 1];
            EXPECT_EQ((uint8_t)CartoonGlyph_NucleicSugar, g.point_kind[c1]);
            EXPECT_EQ((uint8_t)CartoonGlyph_NucleicBase,  g.point_kind[n9]);
            EXPECT_TRUE(point_has_atom_named(g, sys, c1, "C1'"));
            EXPECT_TRUE(point_has_element(g, sys, n9, 7));
            sticks += 1;
        }
    }
    EXPECT_EQ(nucleotides, sticks);

    // A slab picks as its first corner, after the points and the control points
    const uint32_t first_slab = g.num_points + g.num_tube_points;
    EXPECT_EQ(g.slab_point[0], cartoon_glyphs_picked_point(g, first_slab));
    EXPECT_EQ(UINT32_MAX, cartoon_glyphs_picked_point(g, first_slab + g.num_slabs));

    cartoon_glyphs_free(&g);
    md_arena_allocator_destroy(alloc);
}

// The atoms of the corners of each slab, each slab a key of its atoms in order of index, sorted: the slabs as a set
static md_array(uint64_t) slab_keys(const CartoonGlyphs& g, md_allocator_i* alloc) {
    md_array(uint64_t) keys = nullptr;
    for (uint32_t s = 0; s < g.num_slabs; ++s) {
        int32_t atoms[8];
        const uint32_t m = g.slab_offset[s + 1] - g.slab_offset[s];
        for (uint32_t k = 0; k < m && k < 8; ++k) atoms[k] = g.point_atom[g.point_offset[g.slab_point[g.slab_offset[s] + k]]];
        for (uint32_t i = 1; i < m && i < 8; ++i) {
            for (uint32_t j = i; j > 0 && atoms[j - 1] > atoms[j]; --j) {
                const int32_t t = atoms[j]; atoms[j] = atoms[j - 1]; atoms[j - 1] = t;
            }
        }
        uint64_t key = 1469598103934665603ull;
        for (uint32_t k = 0; k < m && k < 8; ++k) key = (key ^ (uint64_t)(uint32_t)atoms[k]) * 1099511628211ull;
        md_array_push(keys, key, alloc);
    }
    std::sort(keys, keys + md_array_size(keys));
    return keys;
}

// A nucleic acid from an XYZ file, the element and the position of each atom: without components, its backbones are found
// from the bonds alone (md_util_system_infer) and the rings of each nucleotide from its C1'. The same slabs as of the
// same atoms read with their residues, and the same band, though no point stands for a component.
UTEST(cartoon_glyphs, nucleic_acids_without_components) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));
    md_system_t ref = { .alloc = alloc };
    md_system_state_t ref_st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&ref, &ref_st, STR_LIT(VIAMD_TEST_DATA_DIR "/nucl-dna.gro")));
    ASSERT_TRUE(md_util_system_infer(&ref, &ref_st, MD_UTIL_INFER_ALL));

    char* text = (char*)md_alloc(alloc, 64 * (ref.atom.count + 1));
    size_t len = (size_t)snprintf(text, 64, "%zu\nxyz\n", ref.atom.count);
    for (size_t i = 0; i < ref.atom.count; ++i) {
        const str_t sym = md_atomic_number_symbol(md_atom_atomic_number(&ref.atom, i));
        const vec3_t p = md_state_coord(&ref_st, i);
        len += (size_t)snprintf(text + len, 64, "%.*s %.4f %.4f %.4f\n", (int)sym.len, sym.ptr, p.x, p.y, p.z);
    }
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_xyz_system_init_from_str(&sys, &st, str_t{ text, len }, MD_XYZ_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    ASSERT_EQ((size_t)0, sys.component.count);
    ASSERT_EQ(ref.nucleic_backbone.segment.count, sys.nucleic_backbone.segment.count);

    CartoonGlyphs g_ref = {};
    CartoonGlyphs g = {};
    cartoon_glyphs_init(&g_ref, ref, md_gl_mol_t{}, alloc);
    cartoon_glyphs_init(&g, sys, md_gl_mol_t{}, alloc);
    ASSERT_GT(g_ref.num_slabs, 0u);
    EXPECT_EQ(g_ref.num_points, g.num_points);
    EXPECT_EQ(g_ref.num_tubes, g.num_tubes);
    EXPECT_EQ(g_ref.num_tube_points, g.num_tube_points);
    EXPECT_EQ(g_ref.num_slabs, g.num_slabs);

    md_array(uint64_t) keys_ref = slab_keys(g_ref, alloc);
    md_array(uint64_t) keys = slab_keys(g, alloc);
    ASSERT_EQ(md_array_size(keys_ref), md_array_size(keys));
    EXPECT_EQ(0, memcmp(keys_ref, keys, md_array_bytes(keys)));

    for (uint32_t i = 0; i < g.num_points; ++i) EXPECT_EQ(-1, g.point_comp[i]);
    for (uint32_t i = 0; i < g.num_tube_points; ++i) {
        ASSERT_LT(g.tube_face[i], g.num_points);
        EXPECT_EQ((uint8_t)CartoonGlyph_NucleicBase, g.point_kind[g.tube_face[i]]);
    }

    cartoon_glyphs_free(&g);
    cartoon_glyphs_free(&g_ref);
    md_arena_allocator_destroy(alloc);
}

// Where a point is in a frame: the weighted sum of its atoms, each in its periodic image nearest to its reference
static vec3_t point_position(const CartoonGlyphs& g, const md_system_state_t& st, uint32_t p) {
    const vec3_t ref = st.xyz[g.point_ref[p]];
    vec3_t pos = vec3_zero();
    for (uint32_t j = g.point_offset[p]; j < g.point_offset[p + 1]; ++j) {
        vec3_t d = vec3_sub(st.xyz[g.point_atom[j]], ref);
        if (md_unitcell_flags(&st.unitcell)) md_util_min_image_vec3(&d, 1, &st.unitcell);
        pos = vec3_add(pos, vec3_mul1(vec3_add(ref, d), g.point_weight[j]));
    }
    return pos;
}

// The nucleic acids of Martini: a band through the BB2 of each nucleotide, and its base with the rings of an ideal base
// (two for a purine) and a stick from the band to it; no sugar. The rings come out where the beads of the frame put
// them: their edges as long as the bonds of a ring, flat, the glycosidic N a corner of them.
UTEST(cartoon_glyphs, martini_nucleic_acids_have_the_rings_of_ideal_bases) {
    static const char* names[] = { "martini2_dna", "martini2_rna", "martini3_rna" };
    for (size_t n = 0; n < ARRAY_SIZE(names); ++n) {
        md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
        md_system_t sys = { .alloc = alloc };
        md_system_state_t st = { .alloc = alloc };
        char path[256];
        snprintf(path, sizeof(path), VIAMD_TEST_DATA_DIR "/martini/%s.gro", names[n]);
        ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, str_from_cstr(path)));
        ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
        const md_nucleic_backbone_data_t& bb = sys.nucleic_backbone;
        ASSERT_EQ((size_t)2, bb.range.count);

        uint32_t expected_slabs = 0;
        for (size_t s = 0; s < bb.segment.count; ++s) {
            const str_t name = md_component_name(&sys.component, bb.segment.comp_idx[s]);
            const bool purine = str_find_char(NULL, name, 'A') || str_find_char(NULL, name, 'G');
            expected_slabs += purine ? 3 : 2;
        }

        CartoonGlyphs g = {};
        cartoon_glyphs_init(&g, sys, md_gl_mol_t{}, alloc);
        EXPECT_EQ_MSG(expected_slabs, g.num_slabs, names[n]);
        EXPECT_EQ(2u, g.num_tubes);
        EXPECT_EQ((uint32_t)bb.segment.count, g.num_tube_points);
        EXPECT_EQ(0u, g.num_spheres);

        // Every point a weighted sum of its atoms: the weights sum to 1
        for (uint32_t i = 0; i < g.num_points; ++i) {
            float sum = 0.0f;
            for (uint32_t j = g.point_offset[i]; j < g.point_offset[i + 1]; ++j) sum += g.point_weight[j];
            EXPECT_NEAR(1.0f, sum, 1e-3f);
        }
        // The band through the BB2, facing the glycosidic N
        for (uint32_t i = 0; i < g.num_tube_points; ++i) {
            const uint32_t p = g.tube_point[i];
            EXPECT_EQ((uint8_t)CartoonGlyph_NucleicBackbone, g.point_kind[p]);
            EXPECT_TRUE(point_has_atom_named(g, sys, p, "BB2"));
            ASSERT_LT(g.tube_face[i], g.num_points);
            EXPECT_EQ((uint8_t)CartoonGlyph_NucleicBase, g.point_kind[g.tube_face[i]]);
        }

        uint32_t sticks = 0;
        float shortest_edge = FLT_MAX, longest_edge = 0.0f, worst_plane = 0.0f, longest_stick = 0.0f;
        for (uint32_t s = 0; s < g.num_slabs; ++s) {
            const uint32_t beg = g.slab_offset[s], m = g.slab_offset[s + 1] - beg;
            ASSERT_TRUE(m == 2 || m == 5 || m == 6);
            if (m == 2) {
                // From the band to the glycosidic N, the corner of a ring
                EXPECT_EQ((uint8_t)CartoonGlyph_NucleicBackbone, g.point_kind[g.slab_point[beg]]);
                EXPECT_EQ((uint8_t)CartoonGlyph_NucleicBase, g.point_kind[g.slab_point[beg + 1]]);
                longest_stick = MAX(longest_stick, vec3_distance(point_position(g, st, g.slab_point[beg]), point_position(g, st, g.slab_point[beg + 1])));
                sticks += 1;
                continue;
            }
            vec3_t c[6];
            vec3_t mid = vec3_zero();
            for (uint32_t k = 0; k < m; ++k) {
                EXPECT_EQ((uint8_t)CartoonGlyph_NucleicBase, g.point_kind[g.slab_point[beg + k]]);
                c[k] = point_position(g, st, g.slab_point[beg + k]);
                mid = vec3_add(mid, vec3_mul1(c[k], 1.0f / (float)m));
            }
            for (uint32_t k = 0; k < m; ++k) {
                const float d = vec3_distance(c[k], c[(k + 1) % m]);
                shortest_edge = MIN(shortest_edge, d);
                longest_edge = MAX(longest_edge, d);
            }
            const vec3_t normal = vec3_normalize(vec3_cross(vec3_sub(c[0], mid), vec3_sub(c[1], mid)));
            for (uint32_t k = 0; k < m; ++k) worst_plane = MAX(worst_plane, fabsf(vec3_dot(vec3_sub(c[k], mid), normal)));
        }
        EXPECT_EQ((uint32_t)bb.segment.count, sticks);
        // The bonds of the rings of the ideal bases are 1.31 to 1.45 long; the purines of Martini 2 flex a little
        EXPECT_GT(shortest_edge, 1.26f);
        EXPECT_LT(longest_edge, 1.50f);
        EXPECT_LT(worst_plane, 0.05f);
        // Within the band and the base: shorter than the longest bond of a Martini nucleotide is drawn (md_gl culls a
        // capsule longer than that as across the periodic boundary)
        EXPECT_LT(longest_stick, 5.0f);

        cartoon_glyphs_free(&g);
        md_arena_allocator_destroy(alloc);
    }
}
