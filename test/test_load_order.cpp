#include "utest.h"

#include <load_order.h>

#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_os.h>

#include <stdio.h>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

/* A topology without coordinates and the files its coordinates come from are two halves of a system, and the
 * order they are loaded in is the user's. These tests load every order and hold that they all end in the same
 * system. The application is a stand in: it keeps what was loaded in which role, and knows the atom counts of the
 * files by name - the files need not be there, the loader tells their kinds from their extensions. */

namespace {

struct App {
    std::map<std::string, size_t> counts;   // atoms per file name
    size_t atoms = 0;
    std::string built_from;                 // the file that built the system
    std::string coordinates;                // its coordinates, for a topology
    std::vector<std::string> supplements;
    std::string run;
    FileList list = {};
    PendingFiles pending = {};
    md_allocator_i* alloc = md_get_heap_allocator();

    ~App() { file_list_free(&list, alloc); }

    static std::string name(md_str_t path) { return std::string(path.ptr, path.len); }
    size_t count(const std::string& n) const { auto it = counts.find(n); return it == counts.end() ? 0 : it->second; }

    LoadResult load(md_str_t path, const loader::LoaderState& ls, FileRole role) {
        const std::string n = name(path);
        switch (role) {
        case FileRole_System:
            if (ls.flags & LoaderFlag_NeedsCoordinates) {
                if (!ls.coordinates || !ls.coordinates[0]) return LoadResult_Pending;
                if (count(n) != count(ls.coordinates)) return LoadResult_Failed;
                coordinates = ls.coordinates;
            } else {
                coordinates.clear();
            }
            atoms = count(n);
            built_from = n;
            supplements.clear();
            run.clear();
            file_list_record(&list, file_entry_from(path, ls), FileRole_System, alloc);
            return LoadResult_System;
        case FileRole_Trajectory:
            if (atoms == 0 || count(n) != atoms) return LoadResult_Failed;
            run = n;
            file_list_record(&list, file_entry_from(path, ls), FileRole_Trajectory, alloc);
            return LoadResult_Trajectory;
        case FileRole_Supplement:
            if (atoms == 0 || count(n) != atoms) return LoadResult_Failed;
            supplements.push_back(n);
            file_list_record(&list, file_entry_from(path, ls), FileRole_Supplement, alloc);
            return LoadResult_Supplement;
        default:
            return LoadResult_Failed;
        }
    }

    LoadOps ops() {
        LoadOps o = {};
        o.user = this;
        o.load = [](void* u, md_str_t p, const loader::LoaderState& ls, FileRole r, bool) { return ((App*)u)->load(p, ls, r); };
        o.system_atoms = [](void* u) { return ((App*)u)->atoms; };
        o.run_frames = [](void* u) { return ((App*)u)->run.empty() ? (size_t)0 : (size_t)10; };
        o.files = [](void* u) { return (const FileList*)&((App*)u)->list; };
        o.atom_count = [](void* u, md_str_t p, const loader::LoaderState&) { return ((App*)u)->count(name(p)); };
        return o;
    }

    LoadResult open(const std::string& path) {
        loader::LoaderState ls = {};
        loader::init(&ls, md_str_from_cstr(path.c_str()), nullptr);
        return load_in_order(&pending, ops(), md_str_from_cstr(path.c_str()), ls);
    }

    // What there is to see once the files are loaded
    std::string str() const {
        std::string s = built_from + (coordinates.empty() ? "" : "(" + coordinates + ")");
        for (const auto& sup : supplements) s += " +" + sup;
        if (!run.empty()) s += " run:" + run;
        if (pending.topology.path[0] || pending.coordinates.path[0] || pending.trajectory.path[0]) s += " (waiting)";
        return s;
    }
};

// Every order of the files, from nothing loaded and from a system of other atoms (viamd's default dataset)
void every_order(int* utest_result, std::vector<std::string> files, const std::map<std::string, size_t>& counts, const char* expected) {
    std::sort(files.begin(), files.end());
    do {
        for (int start = 0; start < 2; ++start) {
            App app;
            app.counts = counts;
            app.counts["default.pdb"] = 153;
            std::string order;
            if (start) {
                app.open("default.pdb");
                order = "default.pdb ";
            }
            for (const auto& f : files) {
                app.open(f);
                order += f + " ";
            }
            const std::string got = app.str();
            if (got != expected) printf("order: %s\n", order.c_str());
            EXPECT_STREQ(expected, got.c_str());
        }
    } while (std::next_permutation(files.begin(), files.end()));
}

}  // namespace

UTEST(viamd_load_order, amber_topology_restart_and_trajectory_in_any_order) {
    every_order(utest_result, { "sys.prmtop", "eq.rst7", "md.nc" }, { {"sys.prmtop", 12}, {"eq.rst7", 12}, {"md.nc", 12} },
        "sys.prmtop(eq.rst7) run:md.nc");
}

UTEST(viamd_load_order, a_topology_and_a_trajectory_in_any_order) {
    every_order(utest_result, { "sys.prmtop", "md.nc" }, { {"sys.prmtop", 12}, {"md.nc", 12} },
        "sys.prmtop(md.nc) run:md.nc");
    every_order(utest_result, { "sys.psf", "md.dcd" }, { {"sys.psf", 30}, {"md.dcd", 30} },
        "sys.psf(md.dcd) run:md.dcd");
}

UTEST(viamd_load_order, psf_binary_coordinates_and_dcd_in_any_order) {
    // A .coor is told by its contents (it may be a PDB file): this one is there, binary, of 30 atoms
    {
        md_file_t f = {};
        ASSERT_TRUE(md_file_open(&f, MD_STR_LIT("out.coor"), MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
        const uint32_t n = 30;
        md_file_write(f, &n, 4);
        const double zero[90] = {};
        md_file_write(f, zero, sizeof(zero));
        md_file_close(&f);
    }
    every_order(utest_result, { "sys.psf", "out.coor", "md.dcd" }, { {"sys.psf", 30}, {"out.coor", 30}, {"md.dcd", 30} },
        "sys.psf(out.coor) run:md.dcd");
    remove("out.coor");
}

UTEST(viamd_load_order, a_structure_file_builds_and_the_topology_fills_it_in_in_any_order) {
    every_order(utest_result, { "sys.pdb", "sys.prmtop", "md.nc" }, { {"sys.pdb", 12}, {"sys.prmtop", 12}, {"md.nc", 12} },
        "sys.pdb +sys.prmtop run:md.nc");
    every_order(utest_result, { "sys.pdb", "sys.psf", "md.dcd" }, { {"sys.pdb", 30}, {"sys.psf", 30}, {"md.dcd", 30} },
        "sys.pdb +sys.psf run:md.dcd");
}

UTEST(viamd_load_order, halves_of_other_systems_wait_and_change_nothing) {
    App app;
    app.counts = { {"default.pdb", 153}, {"sys.prmtop", 12}, {"other.rst7", 99} };
    app.open("default.pdb");
    EXPECT_EQ(LoadResult_Pending, app.open("sys.prmtop"));
    EXPECT_EQ(LoadResult_Pending, app.open("other.rst7"));
    const std::string got = app.str();
    EXPECT_STREQ("default.pdb (waiting)", got.c_str());
}
