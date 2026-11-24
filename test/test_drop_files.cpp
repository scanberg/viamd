#include "utest.h"

#include <drop_files.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>
#include <core/md_file.h>
#include <string.h>
#include <md_system.h>

#include <stdio.h>
#include <string>

/* What matches what among the files of a drop and the systems loaded, on real files (drop_files.h): a topology and the
 * files its coordinates come from, dropped together in any order or onto what is loaded, make one system; coordinates
 * of other atoms make another. The systems are joined here as the application joins them (md_system_merge). */

namespace {

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

// A plan as text: per system its action and its files in load order, by name
static std::string plan_str(const PlacementPlan& p, const DropFiles& d) {
    std::string s;
    static const char* actions[] = { "join", "new", "copy", "ask" };
    for (size_t k = 0; k < p.num_systems; ++k) {
        if (k) s += " | ";
        s += actions[p.systems[k].action];
        for (size_t i = 0; i < p.systems[k].num_entries; ++i) {
            const char* path = d.files[p.systems[k].entries[i].file].entry.path;
            const char* slash = strrchr(path, '/');
            s += " ";
            s += slash ? slash + 1 : path;
        }
    }
    for (size_t i = 0; i < p.num_refused; ++i) s += " refused";
    return s;
}

struct Drop {
    DropFiles* d = new DropFiles();
    PlacementPlan* plan = new PlacementPlan();
    Drop(const DropSystem* loaded = nullptr, size_t num_loaded = 0) { drop_files_init(d, loaded, num_loaded); }
    ~Drop() { drop_files_free(d); delete d; delete plan; }
    void add(const char* path) {
        loader::LoaderState ls = {};
        loader::init(&ls, md_str_from_cstr(path));
        drop_files_add(d, md_str_from_cstr(path), ls);
    }
    std::string place() {
        drop_files_plan(plan, d, 8);
        return plan_str(*plan, *d);
    }
};

// The system a placement of new files makes, joined as the application joins them
static bool build(md_system_t* sys, md_system_state_t* state, Drop& drop, const PlacementSystem& ps) {
    for (size_t i = 0; i < ps.num_entries; ++i) {
        const size_t f = ps.entries[i].file;
        const DropFile* df = drop_files_system(drop.d, f);
        if (!df) return false;
        if (drop.d->files[f].kind == PlacementKind_Trajectory && sys->atom.count) {
            // Its frames the run, which a system without coordinates shows the first of
            if (!md_system_state_has_coords(state) && !md_system_merge(sys, state, &df->sys, &df->state)) return false;
            if (!loader::publish_run(sys, md_str_from_cstr(df->entry.path), MD_STR_LIT("run/a"), 0)) return false;
        } else if (!md_system_merge(sys, state, &df->sys, &df->state)) {
            return false;
        }
    }
    return true;
}

struct Files {
    Files() {
        write_file("viamd_test_drop.prmtop", water_prmtop);
        write_file("viamd_test_drop.rst7", water_rst7);
        write_file("viamd_test_drop.pdb", water_pdb);
        write_file("viamd_test_drop.lammpstrj", water_lammpstrj);
        write_file("viamd_test_drop_other.rst7", "water\n     3\n   1.0000000   1.0000000   1.0000000   1.9572000   1.0000000   1.0000000\n   0.7600000   1.9260000   1.0000000\n");
    }
    ~Files() {
        remove("viamd_test_drop.prmtop"); remove("viamd_test_drop.rst7"); remove("viamd_test_drop.pdb");
        remove("viamd_test_drop.lammpstrj"); remove("viamd_test_drop_other.rst7");
    }
};

}  // namespace

UTEST(viamd_drop_files, a_topology_and_its_coordinates_in_any_order_are_one_system) {
    Files files;
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    struct { const char* coords; float x0; const char* plan; } cases[] = {
        { "viamd_test_drop.rst7", 1.0f, "new viamd_test_drop.prmtop viamd_test_drop.rst7" },               // a restart
        { "viamd_test_drop.pdb", 2.0f, "new viamd_test_drop.prmtop viamd_test_drop.pdb" },                 // a structure of its atoms
        { "viamd_test_drop.lammpstrj", 3.0f, "new viamd_test_drop.prmtop viamd_test_drop.lammpstrj" },     // a trajectory: its run
    };
    for (size_t i = 0; i < MD_COUNTOF(cases); ++i) {
        for (int order = 0; order < 2; ++order) {
            Drop drop;
            if (order == 0) { drop.add("viamd_test_drop.prmtop"); drop.add(cases[i].coords); }
            else            { drop.add(cases[i].coords); drop.add("viamd_test_drop.prmtop"); }
            const std::string got = drop.place();
            EXPECT_STREQ(cases[i].plan, got.c_str());
            ASSERT_EQ(1u, drop.plan->num_systems);

            md_system_t sys = { .alloc = arena };
            md_system_state_t state = { .alloc = arena };
            ASSERT_TRUE(build(&sys, &state, drop, drop.plan->systems[0]));
            EXPECT_EQ(6u, sys.atom.count);
            EXPECT_EQ(6u, md_bond_list_count(&sys.given_bond));
            // The topology's names
            EXPECT_TRUE(md_str_eq_cstr(md_system_residue_name(&sys, 0), "WAT"));
            ASSERT_TRUE(md_system_state_has_coords(&state));
            EXPECT_NEAR(cases[i].x0, state.xyz[0].x, 1e-4f);
            md_system_free(&sys);
        }
    }
    md_arena_allocator_destroy(arena);
}

UTEST(viamd_drop_files, coordinates_of_other_atoms_are_a_system_of_their_own) {
    Files files;
    Drop drop;
    drop.add("viamd_test_drop.prmtop");
    drop.add("viamd_test_drop_other.rst7");
    const std::string got = drop.place();
    EXPECT_STREQ("new viamd_test_drop.prmtop | new viamd_test_drop_other.rst7", got.c_str());
}

UTEST(viamd_drop_files, dropped_onto_what_is_loaded) {
    Files files;
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));

    // A topology alone, loaded: its coordinates join it, as does a trajectory of its atoms
    md_system_t top = { .alloc = arena };
    md_system_state_t top_state = { .alloc = arena };
    loader::LoaderState ls = {};
    loader::init(&ls, MD_STR_LIT("viamd_test_drop.prmtop"));
    ASSERT_EQ(MD_LOAD_OK, loader::load(&top, &top_state, MD_STR_LIT("viamd_test_drop.prmtop"), ls));
    EXPECT_FALSE(md_system_state_has_coords(&top_state));
    const DropSystem loaded[] = { { 7, &top, &top_state, {} } };
    {
        Drop drop(loaded, 1);
        drop.add("viamd_test_drop.rst7");
        const std::string got = drop.place();
        EXPECT_STREQ("join viamd_test_drop.rst7", got.c_str());
        EXPECT_EQ(7u, drop.plan->systems[0].system);
    }
    {
        Drop drop(loaded, 1);
        drop.add("viamd_test_drop.lammpstrj");
        const std::string got = drop.place();
        EXPECT_STREQ("join viamd_test_drop.lammpstrj", got.c_str());
    }

    // A structure loaded: the topology joins it; coordinates would replace its own, so they go into a copy of it
    md_system_t pdb = { .alloc = arena };
    md_system_state_t pdb_state = { .alloc = arena };
    loader::init(&ls, MD_STR_LIT("viamd_test_drop.pdb"));
    ASSERT_EQ(MD_LOAD_OK, loader::load(&pdb, &pdb_state, MD_STR_LIT("viamd_test_drop.pdb"), ls));
    const DropSystem with_pdb[] = { { 3, &pdb, &pdb_state, {} } };
    {
        Drop drop(with_pdb, 1);
        drop.add("viamd_test_drop.prmtop");
        const std::string got = drop.place();
        EXPECT_STREQ("join viamd_test_drop.prmtop", got.c_str());
    }
    {
        Drop drop(with_pdb, 1);
        drop.add("viamd_test_drop.rst7");
        const std::string got = drop.place();
        EXPECT_STREQ("copy viamd_test_drop.rst7", got.c_str());
    }
    {
        // Of other atoms: a system of its own
        Drop drop(with_pdb, 1);
        drop.add("viamd_test_drop_other.rst7");
        const std::string got = drop.place();
        EXPECT_STREQ("new viamd_test_drop_other.rst7", got.c_str());
    }

    md_system_free(&top);
    md_system_free(&pdb);
    md_arena_allocator_destroy(arena);
}

UTEST(viamd_drop_files, a_trajectory_alone_is_a_system_of_points) {
    Drop drop;
    drop.add(VIAMD_TEST_DATA_DIR "/catalyst.xtc");
    const std::string got = drop.place();
    EXPECT_STREQ("new catalyst.xtc", got.c_str());

    // With its structure: one system, the trajectory its run
    Drop both;
    both.add(VIAMD_TEST_DATA_DIR "/catalyst.xtc");
    both.add(VIAMD_TEST_DATA_DIR "/catalyst.gro");
    const std::string got_both = both.place();
    EXPECT_STREQ("new catalyst.gro catalyst.xtc", got_both.c_str());
}
