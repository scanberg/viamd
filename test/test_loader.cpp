#include "utest.h"

#include <loader.h>

#include <core/md_str.h>
#include <core/md_arena_allocator.h>
#include <md_system.h>

/* loader is the table that decides which reader gets a file. It is three parallel arrays indexed by
 * LoaderType plus a lookup, and the arrays are maintained by hand behind #if MD_HDF5 guards - so
 * the failure mode is not a crash but a silent shift, where adding a format one place and not
 * another makes every entry after it describe the wrong reader. Most of what follows checks
 * that the three arrays still agree with each other and with the enum. */

UTEST(viamd_loader, every_type_has_a_name_and_an_extension) {
    /* Undefined is the one entry allowed an empty extension - it is what a lookup returns when it
     * recognises nothing. */
    EXPECT_TRUE(str_empty(loader::type_ext(LoaderType_Undefined)));

    for (LoaderType t = LoaderType_Undefined + 1; t < LoaderType_COUNT; ++t) {
        EXPECT_FALSE(str_empty(loader::type_name(t)));
        EXPECT_FALSE(str_empty(loader::type_ext(t)));
        /* A type that can load nothing would never be reachable. */
        EXPECT_NE(0u, loader::type_flags(t) & (LoaderFlag_System | LoaderFlag_Trajectory | LoaderFlag_Supplemental));
    }
}

UTEST(viamd_loader, extensions_are_unique) {
    /* Two types claiming one extension makes type_from_ext return whichever comes first, and the
     * other format becomes unreachable without anyone noticing. */
    for (LoaderType a = LoaderType_Undefined + 1; a < LoaderType_COUNT; ++a) {
        for (LoaderType b = a + 1; b < LoaderType_COUNT; ++b) {
            EXPECT_FALSE(str_eq_ignore_case(loader::type_ext(a), loader::type_ext(b)));
        }
    }
}

UTEST(viamd_loader, an_extension_round_trips_to_its_own_type) {
    /* The property that keeps the tables aligned: whatever extension a type advertises must look
     * that same type back up. A shifted array fails this immediately. */
    for (LoaderType t = LoaderType_Undefined + 1; t < LoaderType_COUNT; ++t) {
        EXPECT_EQ(t, loader::type_from_ext(loader::type_ext(t)));
    }
}

UTEST(viamd_loader, extension_lookup_ignores_case_and_rejects_the_unknown) {
    EXPECT_EQ((LoaderType)LoaderType_PDB, loader::type_from_ext(STR_LIT("pdb")));
    EXPECT_EQ((LoaderType)LoaderType_PDB, loader::type_from_ext(STR_LIT("PDB")));
    EXPECT_EQ((LoaderType)LoaderType_PDB, loader::type_from_ext(STR_LIT("Pdb")));
    EXPECT_EQ((LoaderType)LoaderType_GRO, loader::type_from_ext(STR_LIT("GRO")));

    EXPECT_EQ((LoaderType)LoaderType_Undefined, loader::type_from_ext(STR_LIT("")));
    EXPECT_EQ((LoaderType)LoaderType_Undefined, loader::type_from_ext(STR_LIT("txt")));
    /* A prefix of a real extension is not a match. */
    EXPECT_EQ((LoaderType)LoaderType_Undefined, loader::type_from_ext(STR_LIT("pd")));
    EXPECT_EQ((LoaderType)LoaderType_Undefined, loader::type_from_ext(STR_LIT("pdbx")));
}

UTEST(viamd_loader, out_of_range_types_fall_back_rather_than_read_past_the_table) {
    /* type_name/type_ext/type_flags are handed values that came from a file or a saved workspace,
     * so a stale or corrupt one has to land on Undefined instead of indexing off the end. */
    const LoaderType bad = LoaderType_COUNT;
    EXPECT_TRUE(str_eq(loader::type_name(bad),  loader::type_name(LoaderType_Undefined)));
    EXPECT_TRUE(str_eq(loader::type_ext(bad),   loader::type_ext(LoaderType_Undefined)));
    EXPECT_EQ(loader::type_flags(bad),          loader::type_flags(LoaderType_Undefined));

    const LoaderType way_off = (LoaderType)0x7FFFFFFFu;
    EXPECT_TRUE(str_eq(loader::type_name(way_off), loader::type_name(LoaderType_Undefined)));
    EXPECT_EQ(loader::type_flags(way_off),         loader::type_flags(LoaderType_Undefined));
}

UTEST(viamd_loader, molecular_and_quantum_formats_are_labelled_as_such) {
    /* The MM/QM flags drive which parts of the UI offer a file, so a format landing in the wrong
     * camp is a visible bug with an invisible cause. */
    EXPECT_NE(0u, loader::type_flags(LoaderType_PDB)    & LoaderFlag_MM);
    EXPECT_EQ(0u, loader::type_flags(LoaderType_PDB)    & LoaderFlag_QM);
    EXPECT_NE(0u, loader::type_flags(LoaderType_MOLDEN) & LoaderFlag_QM);
    EXPECT_EQ(0u, loader::type_flags(LoaderType_MOLDEN) & LoaderFlag_MM);
#if MD_HDF5
    EXPECT_NE(0u, loader::type_flags(LoaderType_TREXIO) & LoaderFlag_QM);
    EXPECT_EQ(0u, loader::type_flags(LoaderType_TREXIO) & LoaderFlag_MM);
    /* VeloxChem h5 is the one format that is both: a system with coordinates and a trajectory, and
     * quantum data alongside. */
    EXPECT_NE(0u, loader::type_flags(LoaderType_VLX_H5) & LoaderFlag_QM);
    EXPECT_NE(0u, loader::type_flags(LoaderType_VLX_H5) & LoaderFlag_MM);
#endif
}

UTEST(viamd_loader, ase_trajectory_contains_a_system_and_a_run) {
    const str_t path = STR_LIT(VIAMD_ASE_TRAJ_TEST_DATA_DIR "/ase_fixed.traj");
    loader::LoaderState reader = {};
    loader::init(&reader, path);
    ASSERT_EQ((LoaderType)LoaderType_ASE_TRAJ, reader.type);
    EXPECT_NE(0u, reader.flags & LoaderFlag_System);
    EXPECT_NE(0u, reader.flags & LoaderFlag_Trajectory);

    md_allocator_i* alloc = md_get_heap_allocator();
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(loader::load(&sys, &state, path, reader));
    EXPECT_EQ(3u, sys.atom.count);
    ASSERT_TRUE(loader::publish_run(&sys, path, STR_LIT("run/ase"), 0));
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/ase/atom/position")) != nullptr);
    md_system_state_free(&state);
    md_system_free(&sys);
}

UTEST(viamd_loader, trajectory_only_formats_do_not_claim_to_be_systems) {
    /* A trajectory has no topology of its own, so offering one as a system would load atoms with no
     * elements or names. */
    const LoaderType traj_only[] = { LoaderType_XTC, LoaderType_TRR, LoaderType_DCD, LoaderType_LAMMPSTRJ };
    for (size_t i = 0; i < ARRAY_SIZE(traj_only); ++i) {
        EXPECT_NE(0u, loader::type_flags(traj_only[i]) & LoaderFlag_Trajectory);
        EXPECT_EQ(0u, loader::type_flags(traj_only[i]) & LoaderFlag_System);
    }
}

/* ---- init(), which is where the extension is only the first guess ---------------------------- */

UTEST(viamd_loader, init_picks_a_reader_from_the_path) {
    loader::LoaderState state = {};

    loader::init(&state, STR_LIT(VIAMD_TEST_DATA_DIR "/molden/h2o_ccpvdz.molden"));
    EXPECT_EQ((LoaderType)LoaderType_MOLDEN, state.type);
    EXPECT_NE(0u, state.flags & LoaderFlag_QM);

    loader::init(&state, STR_LIT(VIAMD_TEST_DATA_DIR "/1ALA-560ns.pdb"));
    EXPECT_EQ((LoaderType)LoaderType_PDB, state.type);

    /* An unrecognised extension on a file that is not a molden file leaves the type undefined and
     * asks for the dialogue rather than guessing a reader - init() is also called on whatever a
     * user drags in, and picking a reader by accident is worse than asking. */
    loader::init(&state, STR_LIT(VIAMD_TEST_DATA_DIR "/script.txt"));
    EXPECT_EQ((LoaderType)LoaderType_Undefined, state.type);
    EXPECT_NE(0u, state.flags & LoaderFlag_RequiresDialogue);

    /* A path with no extension at all, and one that does not exist, both have to be survivable and
     * reach the same answer rather than asserting somewhere inside a reader. */
    loader::init(&state, STR_LIT(VIAMD_TEST_DATA_DIR "/there_is_no_such_file"));
    EXPECT_EQ((LoaderType)LoaderType_Undefined, state.type);
    EXPECT_NE(0u, state.flags & LoaderFlag_RequiresDialogue);
    loader::init(&state, STR_LIT(""));
    EXPECT_EQ((LoaderType)LoaderType_Undefined, state.type);
}

#if MD_HDF5
UTEST(viamd_loader, trexio_is_recognised_by_content_not_only_by_extension) {
    /* .h5 is not one format: VeloxChem and TREXIO both use it, and they are told apart by looking
     * inside. Here the extension already says trexio, so this checks the other half - that the
     * content check agrees rather than overriding a correct guess. */
    loader::LoaderState state = {};
    loader::init(&state, STR_LIT(VIAMD_TEST_DATA_DIR "/trexio/h2o.trexio"));
    EXPECT_EQ((LoaderType)LoaderType_TREXIO, state.type);
    EXPECT_NE(0u, state.flags & LoaderFlag_QM);
}
#endif

#if MD_HDF5
UTEST(viamd_loader, h5md_is_a_system_and_its_own_trajectory) {
    /* One file, both halves: load() gives the system from it and publish_run() the positions, which
     * is the path a structure file with frames of its own takes (init_trajectory_data). */
    const str_t path = STR_LIT(VIAMD_TEST_DATA_DIR "/h5md/peptide_tip3p.h5md");
    loader::LoaderState state = {};
    loader::init(&state, path);
    EXPECT_EQ((LoaderType)LoaderType_H5MD, state.type);
    const LoaderFlags expected = LoaderFlag_System | LoaderFlag_Trajectory | LoaderFlag_MM | LoaderFlag_Topology;
    EXPECT_EQ(expected, state.flags & expected);

    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {};
    sys.alloc = arena;
    md_system_state_t sys_state = {};
    sys_state.alloc = arena;
    ASSERT_TRUE(loader::load(&sys, &sys_state, path, state));
    EXPECT_EQ(1577u, sys.atom.count);
    EXPECT_LT(0u, sys.bond.count);
    EXPECT_TRUE(loader::publish_run(&sys, path, STR_LIT("run/h5md"), 0));
    const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/h5md/time"));
    ASSERT_TRUE(time != NULL);
    EXPECT_EQ(5u, time->format.shape[0]);
    md_vm_arena_destroy(arena);
}
#endif

#if MD_HDF5
UTEST(viamd_loader, a_veloxchem_h5_stays_veloxchem) {
    /* The mirror of the case above: a .h5 that is not TREXIO must not be rerouted by the content
     * check. Getting this backwards would send every VeloxChem file to the wrong reader. */
    loader::LoaderState state = {};
    loader::init(&state, STR_LIT(VIAMD_TEST_DATA_DIR "/vlx/h2o.h5"));
    EXPECT_EQ((LoaderType)LoaderType_VLX_H5, state.type);
    EXPECT_NE(0u, state.flags & LoaderFlag_QM);
}
#endif
