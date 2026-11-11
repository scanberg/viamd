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
    EXPECT_TRUE(md_str_empty(loader::type_ext(LoaderType_Undefined)));

    for (LoaderType t = LoaderType_Undefined + 1; t < LoaderType_COUNT; ++t) {
        EXPECT_FALSE(md_str_empty(loader::type_name(t)));
        EXPECT_FALSE(md_str_empty(loader::type_ext(t)));
        /* A type that can load nothing would never be reachable. Coordinates without atoms are loaded for a
         * topology. */
        EXPECT_NE(0u, loader::type_flags(t) & (LoaderFlag_System | LoaderFlag_Trajectory | LoaderFlag_Supplemental | LoaderFlag_Coordinates));
    }
}

UTEST(viamd_loader, extensions_are_unique) {
    /* Two types claiming one extension makes type_from_ext return whichever comes first, and the
     * other format becomes unreachable without anyone noticing. */
    for (LoaderType a = LoaderType_Undefined + 1; a < LoaderType_COUNT; ++a) {
        for (LoaderType b = a + 1; b < LoaderType_COUNT; ++b) {
            EXPECT_FALSE(md_str_eq_ignore_case(loader::type_ext(a), loader::type_ext(b)));
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
    EXPECT_EQ((LoaderType)LoaderType_PDB, loader::type_from_ext(MD_STR_LIT("pdb")));
    EXPECT_EQ((LoaderType)LoaderType_PDB, loader::type_from_ext(MD_STR_LIT("PDB")));
    EXPECT_EQ((LoaderType)LoaderType_PDB, loader::type_from_ext(MD_STR_LIT("Pdb")));
    EXPECT_EQ((LoaderType)LoaderType_GRO, loader::type_from_ext(MD_STR_LIT("GRO")));

    EXPECT_EQ((LoaderType)LoaderType_Undefined, loader::type_from_ext(MD_STR_LIT("")));
    EXPECT_EQ((LoaderType)LoaderType_Undefined, loader::type_from_ext(MD_STR_LIT("txt")));
    /* A prefix of a real extension is not a match. */
    EXPECT_EQ((LoaderType)LoaderType_Undefined, loader::type_from_ext(MD_STR_LIT("pd")));
    EXPECT_EQ((LoaderType)LoaderType_Undefined, loader::type_from_ext(MD_STR_LIT("pdbx")));
}

UTEST(viamd_loader, out_of_range_types_fall_back_rather_than_read_past_the_table) {
    /* type_name/type_ext/type_flags are handed values that came from a file or a saved workspace,
     * so a stale or corrupt one has to land on Undefined instead of indexing off the end. */
    const LoaderType bad = LoaderType_COUNT;
    EXPECT_TRUE(md_str_eq(loader::type_name(bad),  loader::type_name(LoaderType_Undefined)));
    EXPECT_TRUE(md_str_eq(loader::type_ext(bad),   loader::type_ext(LoaderType_Undefined)));
    EXPECT_EQ(loader::type_flags(bad),          loader::type_flags(LoaderType_Undefined));

    const LoaderType way_off = (LoaderType)0x7FFFFFFFu;
    EXPECT_TRUE(md_str_eq(loader::type_name(way_off), loader::type_name(LoaderType_Undefined)));
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
    const md_str_t path = MD_STR_LIT(VIAMD_ASE_TRAJ_TEST_DATA_DIR "/ase_fixed.traj");
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
    ASSERT_TRUE(loader::publish_run(&sys, path, MD_STR_LIT("run/ase"), 0));
    EXPECT_TRUE(md_attributes_find(&sys.attributes, MD_STR_LIT("run/ase/atom/position")) != nullptr);
    md_system_state_free(&state);
    md_system_free(&sys);
}

UTEST(viamd_loader, trajectory_only_formats_do_not_claim_to_be_systems) {
    /* A trajectory has no topology of its own, so offering one as a system would load atoms with no
     * elements or names. */
    const LoaderType traj_only[] = { LoaderType_XTC, LoaderType_TRR, LoaderType_DCD, LoaderType_LAMMPSTRJ };
    for (size_t i = 0; i < MD_COUNTOF(traj_only); ++i) {
        EXPECT_NE(0u, loader::type_flags(traj_only[i]) & LoaderFlag_Trajectory);
        EXPECT_EQ(0u, loader::type_flags(traj_only[i]) & LoaderFlag_System);
    }
}

/* ---- init(), which is where the extension is only the first guess ---------------------------- */

UTEST(viamd_loader, init_picks_a_reader_from_the_path) {
    loader::LoaderState state = {};

    loader::init(&state, MD_STR_LIT(VIAMD_TEST_DATA_DIR "/molden/h2o_ccpvdz.molden"));
    EXPECT_EQ((LoaderType)LoaderType_MOLDEN, state.type);
    EXPECT_NE(0u, state.flags & LoaderFlag_QM);

    loader::init(&state, MD_STR_LIT(VIAMD_TEST_DATA_DIR "/1ALA-560ns.pdb"));
    EXPECT_EQ((LoaderType)LoaderType_PDB, state.type);

    /* An unrecognised extension on a file that is not a molden file leaves the type undefined and
     * asks for the dialogue rather than guessing a reader - init() is also called on whatever a
     * user drags in, and picking a reader by accident is worse than asking. */
    loader::init(&state, MD_STR_LIT(VIAMD_TEST_DATA_DIR "/script.txt"));
    EXPECT_EQ((LoaderType)LoaderType_Undefined, state.type);
    EXPECT_NE(0u, state.flags & LoaderFlag_RequiresDialogue);

    /* A path with no extension at all, and one that does not exist, both have to be survivable and
     * reach the same answer rather than asserting somewhere inside a reader. */
    loader::init(&state, MD_STR_LIT(VIAMD_TEST_DATA_DIR "/there_is_no_such_file"));
    EXPECT_EQ((LoaderType)LoaderType_Undefined, state.type);
    EXPECT_NE(0u, state.flags & LoaderFlag_RequiresDialogue);
    loader::init(&state, MD_STR_LIT(""));
    EXPECT_EQ((LoaderType)LoaderType_Undefined, state.type);
}

#if MD_HDF5
UTEST(viamd_loader, trexio_is_recognised_by_content_not_only_by_extension) {
    /* .h5 is not one format: VeloxChem and TREXIO both use it, and they are told apart by looking
     * inside. Here the extension already says trexio, so this checks the other half - that the
     * content check agrees rather than overriding a correct guess. */
    loader::LoaderState state = {};
    loader::init(&state, MD_STR_LIT(VIAMD_TEST_DATA_DIR "/trexio/h2o.trexio"));
    EXPECT_EQ((LoaderType)LoaderType_TREXIO, state.type);
    EXPECT_NE(0u, state.flags & LoaderFlag_QM);
}
#endif

#if MD_HDF5
UTEST(viamd_loader, h5md_is_a_system_and_its_own_trajectory) {
    /* One file, both halves: load() gives the system from it and publish_run() the positions, which
     * is the path a structure file with frames of its own takes (init_trajectory_data). */
    const md_str_t path = MD_STR_LIT(VIAMD_TEST_DATA_DIR "/h5md/peptide_tip3p.h5md");
    loader::LoaderState state = {};
    loader::init(&state, path);
    EXPECT_EQ((LoaderType)LoaderType_H5MD, state.type);
    const LoaderFlags expected = LoaderFlag_System | LoaderFlag_Trajectory | LoaderFlag_MM | LoaderFlag_Topology;
    EXPECT_EQ(expected, state.flags & expected);

    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {};
    sys.alloc = arena;
    md_system_state_t sys_state = {};
    sys_state.alloc = arena;
    ASSERT_TRUE(loader::load(&sys, &sys_state, path, state));
    EXPECT_EQ(1577u, sys.atom.count);
    EXPECT_LT(0u, sys.bond.count);
    EXPECT_TRUE(loader::publish_run(&sys, path, MD_STR_LIT("run/h5md"), 0));
    const md_attribute_t* time = md_attributes_find(&sys.attributes, MD_STR_LIT("run/h5md/time"));
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
    loader::init(&state, MD_STR_LIT(VIAMD_TEST_DATA_DIR "/vlx/h2o.h5"));
    EXPECT_EQ((LoaderType)LoaderType_VLX_H5, state.type);
    EXPECT_NE(0u, state.flags & LoaderFlag_QM);
}
#endif

/* ---- A topology without coordinates, and coordinates without atoms ------------------------------ */

#include <core/md_os.h>
#include <md_util.h>
#include <string.h>
#include <stdio.h>

static bool write_file(const char* path, const char* text) {
    md_file_t f = {};
    if (!md_file_open(&f, md_str_from_cstr(path), MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) return false;
    md_file_write(f, text, strlen(text));
    md_file_close(&f);
    return true;
}

/* Two waters, as Amber writes them */
static const char* water_prmtop =
    "%VERSION  VERSION_STAMP = V0001.000  DATE = 11/02/26  12:00:00\n"
    "%FLAG TITLE\n%FORMAT(20a4)\nwater\n"
    "%FLAG POINTERS\n%FORMAT(10I8)\n"
    "       6       2       6       0       0       0       0       0       0       0\n"
    "       0       2       0       0       0       1       0       0       2       0\n"
    "       0       0       0       0       0       0       0       0       3       0\n"
    "       0\n"
    "%FLAG ATOM_NAME\n%FORMAT(20a4)\nO   H1  H2  O   H1  H2  \n"
    "%FLAG CHARGE\n%FORMAT(5E16.8)\n"
    " -1.51973982E+01  7.59869910E+00  7.59869910E+00 -1.51973982E+01  7.59869910E+00\n"
    "  7.59869910E+00\n"
    "%FLAG ATOMIC_NUMBER\n%FORMAT(10I8)\n       8       1       1       8       1       1\n"
    "%FLAG MASS\n%FORMAT(5E16.8)\n"
    "  1.60000000E+01  1.00800000E+00  1.00800000E+00  1.60000000E+01  1.00800000E+00\n"
    "  1.00800000E+00\n"
    "%FLAG RESIDUE_LABEL\n%FORMAT(20a4)\nWAT WAT \n"
    "%FLAG RESIDUE_POINTER\n%FORMAT(10I8)\n       1       4\n"
    "%FLAG BONDS_INC_HYDROGEN\n%FORMAT(10I8)\n"
    "       0       3       1       0       6       1       3       6       2       9\n"
    "      12       1       9      15       1      12      15       2\n"
    "%FLAG BONDS_WITHOUT_HYDROGEN\n%FORMAT(10I8)\n\n";

static const char* water_rst7 =
    "water\n"
    "     6\n"
    "   1.0000000   1.0000000   1.0000000   1.9572000   1.0000000   1.0000000\n"
    "   0.7600000   1.9260000   1.0000000   5.0000000   5.0000000   5.0000000\n"
    "   5.9572000   5.0000000   5.0000000   4.7600000   5.9260000   5.0000000\n";

static const char* water_pdb =
    "ATOM      1  O   WAT     1       2.000   2.000   2.000  1.00  0.00           O\n"
    "ATOM      2  H1  WAT     1       2.957   2.000   2.000  1.00  0.00           H\n"
    "ATOM      3  H2  WAT     1       1.760   2.926   2.000  1.00  0.00           H\n"
    "ATOM      4  O   WAT     2       6.000   6.000   6.000  1.00  0.00           O\n"
    "ATOM      5  H1  WAT     2       6.957   6.000   6.000  1.00  0.00           H\n"
    "ATOM      6  H2  WAT     2       5.760   6.926   6.000  1.00  0.00           H\n"
    "END\n";

static const char* water_lammpstrj =
    "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n6\nITEM: BOX BOUNDS pp pp pp\n0 20\n0 20\n0 20\nITEM: ATOMS id type x y z\n"
    "1 1 3.0 3.0 3.0\n2 2 3.957 3.0 3.0\n3 2 2.76 3.926 3.0\n4 1 7.0 7.0 7.0\n5 2 7.957 7.0 7.0\n6 2 6.76 7.926 7.0\n"
    "ITEM: TIMESTEP\n100\nITEM: NUMBER OF ATOMS\n6\nITEM: BOX BOUNDS pp pp pp\n0 20\n0 20\n0 20\nITEM: ATOMS id type x y z\n"
    "1 1 3.5 3.0 3.0\n2 2 4.457 3.0 3.0\n3 2 3.26 3.926 3.0\n4 1 7.5 7.0 7.0\n5 2 8.457 7.0 7.0\n6 2 7.26 7.926 7.0\n";

UTEST(viamd_loader, amber_files_are_told_by_their_contents) {
    ASSERT_TRUE(write_file("viamd_test_water.top", water_prmtop));
    ASSERT_TRUE(write_file("prmtop", water_prmtop));
    ASSERT_TRUE(write_file("inpcrd", water_rst7));
    ASSERT_TRUE(write_file("viamd_test_water.coor", water_pdb));

    loader::LoaderState s = {};
    /* .top is GROMACS' too: the contents say it is Amber's */
    loader::init(&s, MD_STR_LIT("viamd_test_water.top"));
    EXPECT_EQ((LoaderType)LoaderType_AMBER_PRMTOP, s.type);
    EXPECT_NE(0u, s.flags & LoaderFlag_NeedsCoordinates);
    EXPECT_EQ(6u, loader::atom_count(MD_STR_LIT("viamd_test_water.top"), s));
    /* Amber's own names have no extension */
    loader::init(&s, MD_STR_LIT("prmtop"));
    EXPECT_EQ((LoaderType)LoaderType_AMBER_PRMTOP, s.type);
    loader::init(&s, MD_STR_LIT("inpcrd"));
    EXPECT_EQ((LoaderType)LoaderType_AMBER_RESTART, s.type);
    EXPECT_NE(0u, s.flags & LoaderFlag_Coordinates);
    EXPECT_EQ(6u, loader::atom_count(MD_STR_LIT("inpcrd"), s));
    /* NAMD writes a PDB file as .coor when binaryoutput is off */
    loader::init(&s, MD_STR_LIT("viamd_test_water.coor"));
    EXPECT_EQ((LoaderType)LoaderType_PDB, s.type);

    remove("viamd_test_water.top"); remove("prmtop"); remove("inpcrd"); remove("viamd_test_water.coor");
}

UTEST(viamd_loader, a_topology_loads_with_coordinates_from_any_of_their_files) {
    ASSERT_TRUE(write_file("viamd_test_water.prmtop", water_prmtop));
    ASSERT_TRUE(write_file("viamd_test_water.rst7", water_rst7));
    ASSERT_TRUE(write_file("viamd_test_water.pdb", water_pdb));
    ASSERT_TRUE(write_file("viamd_test_water.lammpstrj", water_lammpstrj));
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));

    loader::LoaderState ls = {};
    loader::init(&ls, MD_STR_LIT("viamd_test_water.prmtop"));
    ASSERT_EQ((LoaderType)LoaderType_AMBER_PRMTOP, ls.type);

    struct { const char* file; float x0; } coords[] = {
        { "viamd_test_water.rst7", 1.0f },          // a restart
        { "viamd_test_water.pdb", 2.0f },           // a structure file of the same atoms
        { "viamd_test_water.lammpstrj", 3.0f },     // the first frame of a trajectory
    };
    for (size_t i = 0; i < MD_COUNTOF(coords); ++i) {
        md_system_t sys = {};
        sys.alloc = arena;
        md_system_state_t state = {};
        state.alloc = arena;
        ls.coordinates = coords[i].file;
        ASSERT_TRUE(loader::load(&sys, &state, MD_STR_LIT("viamd_test_water.prmtop"), ls));
        EXPECT_EQ(6u, sys.atom.count);
        EXPECT_EQ(6u, sys.bond.count);
        ASSERT_TRUE(md_system_state_has_coords(&state));
        EXPECT_NEAR(coords[i].x0, state.xyz[0].x, 1e-4f);
        /* The run the trajectory was read through is gone again */
        EXPECT_TRUE(md_attributes_find(&sys.attributes, MD_STR_LIT("run/__coordinates/time")) == nullptr);
        md_system_free(&sys);
    }

    /* Without coordinates, or with those of other atoms, there is no system */
    md_system_t sys = {};
    sys.alloc = arena;
    md_system_state_t state = {};
    state.alloc = arena;
    ls.coordinates = nullptr;
    EXPECT_FALSE(loader::load(&sys, &state, MD_STR_LIT("viamd_test_water.prmtop"), ls));
    ASSERT_TRUE(write_file("viamd_test_water.rst7", "water\n     3\n   1.0000000   1.0000000   1.0000000   1.9572000   1.0000000   1.0000000\n   0.7600000   1.9260000   1.0000000\n"));
    ls.coordinates = "viamd_test_water.rst7";
    EXPECT_FALSE(loader::load(&sys, &state, MD_STR_LIT("viamd_test_water.prmtop"), ls));
    md_system_free(&sys);

    /* Applied to the structure file's system instead: the topology fills it in */
    md_system_t pdb = {};
    pdb.alloc = arena;
    md_system_state_t pdb_state = {};
    pdb_state.alloc = arena;
    loader::LoaderState ps = {};
    loader::init(&ps, MD_STR_LIT("viamd_test_water.pdb"));
    ASSERT_TRUE(loader::load(&pdb, &pdb_state, MD_STR_LIT("viamd_test_water.pdb"), ps));
    md_util_system_infer(&pdb, &pdb_state, MD_UTIL_INFER_ALL, NULL);
    ASSERT_TRUE(loader::load_supplemental(&pdb, MD_STR_LIT("viamd_test_water.prmtop"), ls, MD_STR_LIT("")));
    EXPECT_EQ(6u, pdb.bond.count);
    EXPECT_TRUE(md_attributes_find(&pdb.attributes, MD_STR_LIT("atom/charge")) != nullptr);
    md_system_free(&pdb);

    md_arena_allocator_destroy(arena);
    remove("viamd_test_water.prmtop"); remove("viamd_test_water.rst7"); remove("viamd_test_water.pdb"); remove("viamd_test_water.lammpstrj");
}

UTEST(viamd_loader, trajectories_say_how_many_atoms_they_hold) {
    /* What lets a trajectory loaded before its topology wait for it rather than fail against a system of other atoms */
    loader::LoaderState s = {};
    loader::init(&s, MD_STR_LIT(VIAMD_TEST_DATA_DIR "/catalyst.xtc"));
    EXPECT_EQ(1336u, loader::atom_count(MD_STR_LIT(VIAMD_TEST_DATA_DIR "/catalyst.xtc"), s));
    loader::init(&s, MD_STR_LIT(VIAMD_TEST_DATA_DIR "/tryptophan-md.trr"));
    EXPECT_EQ(6495u, loader::atom_count(MD_STR_LIT(VIAMD_TEST_DATA_DIR "/tryptophan-md.trr"), s));
}
