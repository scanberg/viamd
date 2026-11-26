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
#include <core/md_vec_math.h>

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

// The build at load (density_cache_preload): on the pool, into the cache, where a consumer at the same places then
// finds it
UTEST(density_cache, preload_builds_the_tiles_on_the_pool) {
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
    md_float3_t* xyz = (md_float3_t*)md_alloc(heap, sizeof(md_float3_t) * na);
    double* c = (double*)md_alloc(heap, sizeof(double) * 3 * na);
    const md_attribute_t* ca = md_attributes_find(&sys.attributes, MD_STR_LIT("qm/atom/coordinate"));     // Angstrom
    ASSERT_TRUE(ca != nullptr);
    ASSERT_EQ(3 * na, md_attribute_extract_f64(c, 3 * na, &sys.attributes, ca, md_attribute_slice_all(), md_unit_none()));
    for (size_t i = 0; i < na; ++i) {
        xyz[i] = {(float)(c[3 * i + 0] * 1.8897261246257702), (float)(c[3 * i + 1] * 1.8897261246257702), (float)(c[3 * i + 2] * 1.8897261246257702)};
    }

    md_str_t path = MD_STR_LIT("orbital/total/density");
    if (!md_attributes_find(&sys.attributes, path)) path = MD_STR_LIT("orbital/alpha/density");
    DensityCache* cache = density_cache_create();
    density_cache_preload(cache, sys, path, xyz, na);
    task_system::wait_all();
    const size_t bytes = density_cache_bytes(cache);
    EXPECT_GT(bytes, 0u);
    std::shared_ptr<const DensityTiles> a = density_cache_tiles(cache, sys, path, &basis, xyz[0].elem, sizeof(md_float3_t));
    ASSERT_TRUE(a != nullptr);
    EXPECT_EQ(sizeof(double) * 256 * a->dm.num_tiles, bytes);     // the preloaded one, not a second
    EXPECT_EQ(a.get(), density_cache_tiles(cache, sys, path, &basis, xyz[0].elem, sizeof(md_float3_t)).get());

    density_cache_clear(cache);
    density_cache_preload(cache, sys, MD_STR_LIT("orbital/no/density"), xyz, na);      // no orbitals for it: nothing kept
    task_system::wait_all();
    EXPECT_EQ(0u, density_cache_bytes(cache));
    density_cache_preload(cache, sys, path, nullptr, 0);                               // no places: nothing started
    task_system::wait_all();
    EXPECT_EQ(0u, density_cache_bytes(cache));

    density_cache_destroy(cache);
    md_free(heap, xyz, sizeof(md_float3_t) * na);
    md_free(heap, c, sizeof(double) * 3 * na);
    md_gto_basis_free(&basis, heap);
    md_system_state_free(&state);
    md_system_free(&sys);
}
