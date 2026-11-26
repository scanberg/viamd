#include "utest.h"
#include <density_cache.h>
#include <loader.h>
#include <task_system.h>

#include <md_system.h>
#include <md_qm.h>
#include <md_gto.h>
#include <md_gto_dmat.h>
#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_unit.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

// The pool, once for the process (it has no re-entrant initialize): whichever test brings it up first
static void density_test_task_system() {
    if (!task_system::executor()) {
        task_system::initialize(4);
        atexit([]() { task_system::shutdown(); });
    }
}

// A density read once is kept as its packed triangle (the symmetric part), handed out again while what it is made
// of keeps its version, read again once it changes (the old one staying with who holds it), and gone on a clear.
UTEST(density_cache, keeps_the_symmetric_triangle_until_its_version_changes) {
    md_allocator_i* heap = md_get_heap_allocator();
    md_system_t sys = {.alloc = heap};
    sys.attributes.alloc = heap;
    const md_str_t path = MD_STR_LIT("orbital/total/density");
    const double m[9] = { 1.0, 2.0, 3.0,
                          4.0, 5.0, 6.0,
                          7.0, 8.0, 9.0 };
    const md_attribute_id_t id = md_qm_publish_matrix(&sys, path, MD_STR_LIT("Density"), md_unit_none(), m, 3, 3);
    ASSERT_NE(MD_ATTRIBUTE_INVALID, id);

    DensityCache* cache = density_cache_create();
    std::shared_ptr<const PackedDensity> a = density_cache_packed(cache, sys, path, nullptr);
    ASSERT_TRUE(a != nullptr);
    ASSERT_EQ(3u, a->dim);
    const double sym[6] = { 1.0, 3.0, 5.0, 5.0, 7.0, 9.0 };   // (0,0) (0,1) (0,2) (1,1) (1,2) (2,2)
    for (int k = 0; k < 6; ++k) EXPECT_EQ(sym[k], a->tri[k]);
    EXPECT_EQ(sizeof(double) * 6, density_cache_bytes(cache));

    std::shared_ptr<const PackedDensity> b = density_cache_packed(cache, sys, path, nullptr);
    EXPECT_EQ(a.get(), b.get());

    // New contents: read again, and the one held before is left as it was
    double* data = (double*)md_attributes_data(&sys.attributes, id, MD_ATTRIBUTE_TYPE_F64);
    ASSERT_TRUE(data != nullptr);
    data[1] = 4.0;      // (0,1) and (1,0) now both 4
    md_attributes_touch(&sys.attributes, id);
    std::shared_ptr<const PackedDensity> c = density_cache_packed(cache, sys, path, nullptr);
    ASSERT_TRUE(c != nullptr);
    EXPECT_NE(a.get(), c.get());
    EXPECT_EQ(4.0, c->tri[1]);
    EXPECT_EQ(3.0, a->tri[1]);
    EXPECT_EQ(sizeof(double) * 6, density_cache_bytes(cache));    // one entry per path

    // Without a cache: read, not kept
    std::shared_ptr<const PackedDensity> d = density_cache_packed(nullptr, sys, path, nullptr);
    ASSERT_TRUE(d != nullptr);
    EXPECT_NE(c.get(), d.get());
    EXPECT_EQ(0, memcmp(c->tri, d->tri, sizeof(double) * 6));

    EXPECT_TRUE(density_cache_packed(cache, sys, MD_STR_LIT("orbital/alpha/density"), nullptr) == nullptr);

    density_cache_clear(cache);
    EXPECT_EQ(0u, density_cache_bytes(cache));
    EXPECT_EQ(4.0, c->tri[1]);      // still held
    density_cache_destroy(cache);
    md_system_free(&sys);
}

// A file's density (h2o in cc-pVDZ from Molden: computed from its orbitals on demand), kept as the triangle of the
// square md_qm hands out
UTEST(density_cache, molden_density_is_the_square_ones_triangle) {
    const md_str_t file = MD_STR_LIT(VIAMD_TEST_DATA_DIR "/molden/h2o_ccpvdz.molden");
    loader::LoaderState reader = {};
    loader::init(&reader, file);
    md_allocator_i* heap = md_get_heap_allocator();
    md_system_t sys = {.alloc = heap};
    md_system_state_t state = {.alloc = heap};
    ASSERT_TRUE(loader::load(&sys, &state, file, reader));

    md_str_t path = MD_STR_LIT("orbital/total/density");
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, path);
    if (!attr) {
        path = MD_STR_LIT("orbital/alpha/density");
        attr = md_attributes_find(&sys.attributes, path);
    }
    ASSERT_TRUE(attr != nullptr);
    const size_t n = md_qm_extract_symmetric_f64(nullptr, 0, &sys.attributes, attr, md_attribute_slice_all(), nullptr);
    ASSERT_GT(n, 0u);
    double* full = (double*)md_alloc(heap, sizeof(double) * n * n);
    ASSERT_EQ(n, md_qm_extract_symmetric_f64(full, n * n, &sys.attributes, attr, md_attribute_slice_all(), nullptr));

    DensityCache* cache = density_cache_create();
    std::shared_ptr<const PackedDensity> p = density_cache_packed(cache, sys, path, nullptr);
    ASSERT_TRUE(p != nullptr);
    ASSERT_EQ(n, p->dim);
    size_t k = 0, differ = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j, ++k) differ += p->tri[k] != 0.5 * (full[i * n + j] + full[j * n + i]);
    }
    EXPECT_EQ(0u, differ);

    density_cache_destroy(cache);
    md_free(heap, full, sizeof(double) * n * n);
    md_system_state_free(&state);
    md_system_free(&sys);
}

// A density's tiles: built from its orbitals once (the same as md_gto_dmat_t's own build), handed out again while the
// orbitals and the atoms' places hold, built again once an atom has moved more than half the margin
UTEST(density_cache, tiles_are_built_once_and_follow_the_atoms) {
    density_test_task_system();
    const md_str_t file = MD_STR_LIT(VIAMD_TEST_DATA_DIR "/molden/h2o_ccpvdz.molden");
    loader::LoaderState reader = {};
    loader::init(&reader, file);
    md_allocator_i* heap = md_get_heap_allocator();
    md_system_t sys = {.alloc = heap};
    md_system_state_t state = {.alloc = heap};
    ASSERT_TRUE(loader::load(&sys, &state, file, reader));
    md_gto_basis_t basis = {};
    ASSERT_TRUE(md_gto_basis_extract_attributes(&basis, &sys.attributes, heap));
    const size_t na = md_gto_basis_num_atoms(&basis);
    float* xyz = (float*)md_alloc(heap, sizeof(float) * 3 * na);
    double* c = (double*)md_alloc(heap, sizeof(double) * 3 * na);
    const md_attribute_t* ca = md_attributes_find(&sys.attributes, MD_STR_LIT("qm/atom/coordinate"));     // Angstrom
    ASSERT_TRUE(ca != nullptr);
    ASSERT_EQ(3 * na, md_attribute_extract_f64(c, 3 * na, &sys.attributes, ca, md_attribute_slice_all(), md_unit_none()));
    for (size_t i = 0; i < 3 * na; ++i) xyz[i] = (float)(c[i] * 1.8897261246257702);

    md_str_t path = MD_STR_LIT("orbital/total/density");
    if (!md_attributes_find(&sys.attributes, path)) path = MD_STR_LIT("orbital/alpha/density");
    DensityCache* cache = density_cache_create();
    std::shared_ptr<const DensityTiles> a = density_cache_tiles(cache, sys, path, &basis, xyz, 0);
    ASSERT_TRUE(a != nullptr);
    EXPECT_EQ(DENSITY_TILES_CUTOFF, a->dm.cutoff);
    EXPECT_EQ(DENSITY_TILES_MARGIN, a->dm.margin);
    EXPECT_EQ(sizeof(double) * 256 * a->dm.num_tiles, density_cache_bytes(cache));

    // the same as a build of its own
    double* rows = nullptr; double* w = nullptr; size_t nao = 0;
    const size_t nmo = md_qm_density_orbitals(&rows, &w, &nao, &sys, path, heap);
    ASSERT_GT(nmo, 0u);
    md_gto_dmat_t ref = {};
    const md_gto_dmat_desc_t desc = { .basis = &basis, .atom_xyz = xyz, .cutoff = DENSITY_TILES_CUTOFF, .margin = DENSITY_TILES_MARGIN };
    ASSERT_TRUE(md_gto_dmat_init(&ref, &desc, heap));
    ASSERT_TRUE(md_gto_dmat_build(&ref, rows, nao, 1, w, nmo, nullptr));
    ASSERT_EQ(ref.num_tiles, a->dm.num_tiles);
    EXPECT_EQ(0, memcmp(ref.tiles, a->dm.tiles, sizeof(double) * 256 * ref.num_tiles));

    // within half the margin: the same tiles; past it: built again
    xyz[0] += 0.4f;
    EXPECT_EQ(a.get(), density_cache_tiles(cache, sys, path, &basis, xyz, 0).get());
    xyz[0] += 0.2f;
    std::shared_ptr<const DensityTiles> b = density_cache_tiles(cache, sys, path, &basis, xyz, 0);
    ASSERT_TRUE(b != nullptr);
    EXPECT_NE(a.get(), b.get());
    EXPECT_EQ(sizeof(double) * 256 * b->dm.num_tiles, density_cache_bytes(cache));     // one entry per path

    // no orbitals for that density: none
    EXPECT_TRUE(density_cache_tiles(cache, sys, MD_STR_LIT("orbital/no/density"), &basis, xyz, 0) == nullptr);

    density_cache_clear(cache);
    EXPECT_EQ(0u, density_cache_bytes(cache));
    density_cache_destroy(cache);
    md_gto_dmat_free(&ref);
    md_free(heap, rows, sizeof(double) * nmo * nao);
    md_free(heap, w, sizeof(double) * nmo);
    md_free(heap, xyz, sizeof(float) * 3 * na);
    md_free(heap, c, sizeof(double) * 3 * na);
    md_gto_basis_free(&basis, heap);
    md_system_state_free(&state);
    md_system_free(&sys);
}

// The read at load (density_cache_preload): on the pool, into the cache, where a consumer then finds it
UTEST(density_cache, preload_reads_it_on_the_pool) {
    density_test_task_system();
    md_allocator_i* heap = md_get_heap_allocator();
    md_system_t sys = {.alloc = heap};
    sys.attributes.alloc = heap;
    const md_str_t path = MD_STR_LIT("orbital/total/density");
    const double m[4] = { 2.0, 0.5, 0.5, 1.0 };
    ASSERT_NE(MD_ATTRIBUTE_INVALID, md_qm_publish_matrix(&sys, path, MD_STR_LIT("Density"), md_unit_none(), m, 2, 2));

    DensityCache* cache = density_cache_create();
    density_cache_preload(cache, sys, path);
    task_system::wait_all();
    EXPECT_EQ(sizeof(double) * 3, density_cache_bytes(cache));
    std::shared_ptr<const PackedDensity> a = density_cache_packed(cache, sys, path, nullptr);
    ASSERT_TRUE(a != nullptr);
    EXPECT_EQ(0.5, a->tri[1]);
    EXPECT_EQ(a.get(), density_cache_packed(cache, sys, path, nullptr).get());

    density_cache_preload(cache, sys, MD_STR_LIT("orbital/alpha/density"));     // none: nothing started
    density_cache_destroy(cache);
    md_system_free(&sys);
}
