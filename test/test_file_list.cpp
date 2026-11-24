#include "utest.h"

#include <file_list.h>
#include <serialization_utils.h>

#include <core/md_allocator.h>
#include <core/md_writer.h>
#include <core/md_file.h>
#include <core/md_path.h>

#include <string.h>
#include <string>

/* A system is an ordered list of files - the first builds it, the rest fill it in - and so is each of its runs -
 * the first opens it, the rest join it. These tests hold the rules that keep the list what is loaded, how a file
 * is read from its place, and the workspace format: [System], [Run] and [Input] tied together by ids. */

namespace {

FileEntry file(const char* path, LoaderType type = LoaderType_Undefined) {
    FileEntry e = {};
    snprintf(e.path, sizeof(e.path), "%s", path);
    e.type = type;
    return e;
}

struct List {
    FileList list = {};
    md_allocator_i* alloc = md_get_heap_allocator();
    ~List() { file_list_free(&list, alloc); }
    bool record(const FileEntry& e, FileRole role) { return file_list_record(&list, e, role, alloc); }

    // "path:role:run" per entry, the run as -1 for the system's own
    std::string str() const {
        std::string s;
        static const char* names[] = { "None", "System", "Supplement", "Trajectory", "Series" };
        for (size_t i = 0; i < file_list_count(&list); ++i) {
            const FileEntry* e = file_list_get(&list, i);
            if (!s.empty()) s += " ";
            s += std::string(e->path) + ":" + names[file_list_role(&list, i)] + ":" + std::to_string(e->run) + (e->missing ? "!" : "");
        }
        return s;
    }
};

}  // namespace

UTEST(viamd_file_list, a_file_on_its_own_takes_the_role_of_its_loader) {
    // What can supplement the loaded system does, else what can build one does, and a trajectory only otherwise
    EXPECT_EQ(FileRole_Series,     file_role(LoaderFlag_Supplemental | LoaderFlag_Temporal));
    EXPECT_EQ(FileRole_Supplement, file_role(LoaderFlag_Supplemental | LoaderFlag_System | LoaderFlag_Trajectory));
    EXPECT_EQ(FileRole_System,     file_role(LoaderFlag_System | LoaderFlag_Trajectory));
    EXPECT_EQ(FileRole_Trajectory, file_role(LoaderFlag_Trajectory));
    EXPECT_EQ(FileRole_None,       file_role(LoaderFlag_None));

    EXPECT_EQ(FileRole_System,     file_role(loader::type_flags(LoaderType_PDB)));
    EXPECT_EQ(FileRole_System,     file_role(loader::type_flags(LoaderType_TPR)));
    EXPECT_EQ(FileRole_Trajectory, file_role(loader::type_flags(LoaderType_XTC)));
    EXPECT_EQ(FileRole_Supplement, file_role(loader::type_flags(LoaderType_ITP)));
    EXPECT_EQ(FileRole_Series,     file_role(loader::type_flags(LoaderType_EDR)));
    EXPECT_EQ(FileRole_Series,     file_role(loader::type_flags(LoaderType_XVG)));
}

UTEST(viamd_file_list, the_role_is_the_place) {
    List l;
    EXPECT_TRUE(l.record(file("a.gro"), FileRole_System));
    EXPECT_TRUE(l.record(file("a.xtc"), FileRole_Trajectory));
    EXPECT_TRUE(l.record(file("a.top"), FileRole_Supplement));    // after the run was opened, but the system's
    EXPECT_TRUE(l.record(file("a.edr"), FileRole_Series));
    EXPECT_TRUE(l.record(file("b.xvg"), FileRole_Series));
    { const std::string got_ = l.str(); EXPECT_STREQ("a.gro:System:-1 a.xtc:Trajectory:0 a.top:Supplement:-1 a.edr:Series:0 b.xvg:Series:0", got_.c_str()); }
    EXPECT_STREQ("a.gro", file_list_system_file(&l.list)->path);
    EXPECT_STREQ("a.xtc", file_list_run_file(&l.list)->path);
}

UTEST(viamd_file_list, building_starts_over) {
    List l;
    l.record(file("a.gro"), FileRole_System);
    l.record(file("a.xtc"), FileRole_Trajectory);
    l.list.atom_count = 10;
    l.record(file("b.pdb"), FileRole_System);
    { const std::string got_ = l.str(); EXPECT_STREQ("b.pdb:System:-1", got_.c_str()); }
    EXPECT_EQ(0u, md_array_size(l.list.runs));
    EXPECT_EQ(0u, l.list.atom_count);
}

UTEST(viamd_file_list, a_trajectory_opens_a_run_and_closing_runs_drops_their_files) {
    List l;
    l.record(file("a.pdb"), FileRole_System);
    l.record(file("a.pdb"), FileRole_Trajectory);   // a multi model PDB: the system's file, and its run's
    l.record(file("a.xvg"), FileRole_Series);
    l.record(file("a.itp"), FileRole_Supplement);
    { const std::string got_ = l.str(); EXPECT_STREQ("a.pdb:System:-1 a.pdb:Trajectory:0 a.xvg:Series:0 a.itp:Supplement:-1", got_.c_str()); }
    EXPECT_STREQ("a.pdb", file_list_run_file(&l.list)->path);

    file_list_drop_runs(&l.list);
    { const std::string got_ = l.str(); EXPECT_STREQ("a.pdb:System:-1 a.itp:Supplement:-1", got_.c_str()); }
    EXPECT_EQ(0u, md_array_size(l.list.runs));
}

UTEST(viamd_file_list, a_series_needs_a_run) {
    List l;
    l.record(file("a.gro"), FileRole_System);
    EXPECT_FALSE(l.record(file("a.edr"), FileRole_Series));
    { const std::string got_ = l.str(); EXPECT_STREQ("a.gro:System:-1", got_.c_str()); }
}

UTEST(viamd_file_list, a_file_loaded_again_stays_in_its_place) {
    List l;
    l.record(file("a.gro"), FileRole_System);
    l.record(file("a.itp", LoaderType_ITP), FileRole_Supplement);
    l.record(file("a.xtc"), FileRole_Trajectory);
    FileEntry gone = file("a.csv", LoaderType_CSV);
    gone.missing = true;
    l.record(gone, FileRole_Series);
    l.record(file("a.itp", LoaderType_ITP), FileRole_Supplement);
    l.record(file("a.csv", LoaderType_CSV), FileRole_Series);   // put back: no longer missing
    { const std::string got_ = l.str(); EXPECT_STREQ("a.gro:System:-1 a.itp:Supplement:-1 a.xtc:Trajectory:0 a.csv:Series:0", got_.c_str()); }
}

// ## Workspace

namespace {

struct Written {
    viamd::serialization_state_t s = {};
    Written(const char* workspace, md_allocator_i* alloc) {
        s.filename = md_str_from_cstr(workspace);
        s.sw = md_str_writer_create(alloc);
    }
    ~Written() { md_str_writer_free(&s.sw); }
    md_str_t text() const { return md_str_writer_to_str(&s.sw); }
};

struct Read {
    WorkspaceFiles files = {};
    SystemPlan plan = {};
    md_allocator_i* alloc = md_get_heap_allocator();
    // Of the system with the id system_id, the first the workspace names by default
    Read(md_str_t text, md_str_t folder, size_t max_runs = 1, int system_id = -1) {
        viamd::deserialization_state_t d = {};
        d.filename = MD_STR_LIT("test.via");
        d.text = text;
        md_str_t section;
        while (viamd::next_section_header(section, d)) {
            workspace_files_read(&files, d, folder, alloc);
        }
        if (system_id < 0) system_id = md_array_size(files.systems) ? files.systems[0].id : 0;
        workspace_files_plan(&plan, &files, system_id, max_runs, alloc);
    }
    ~Read() {
        md_array_free(files.systems, alloc); md_array_free(files.runs, alloc); md_array_free(files.inputs, alloc);
        md_array_free(plan.inputs, alloc);   md_array_free(plan.frame_counts, alloc);
    }
    std::string str() const {
        std::string s;
        static const char* names[] = { "None", "System", "Supplement", "Trajectory", "Series" };
        for (size_t i = 0; i < md_array_size(plan.inputs); ++i) {
            md_str_t name = md_str_from_cstr(plan.inputs[i].entry.path);
            name = md_path_filename(name);
            if (!s.empty()) s += " ";
            s += std::string(name.ptr, name.len) + ":" + names[plan.inputs[i].role] + ":" + std::to_string(plan.inputs[i].run);
        }
        return s;
    }
};

// A path as the workspace reader resolves it: canonical (on Windows a path from the root gets the current drive),
// or as given where it cannot be made so
std::string resolved(const char* path) {
    char buf[1024];
    const size_t len = md_path_write_canonical(buf, sizeof(buf), md_str_from_cstr(path));
    return len ? std::string(buf, len) : std::string(path);
}

bool contains(md_str_t text, const char* what) {
    size_t loc;
    return ((loc = md_str_find_str(text, md_str_from_cstr(what))) != MD_STR_NPOS);
}

}  // namespace

UTEST(viamd_file_list, workspace_round_trip) {
    md_allocator_i* alloc = md_get_heap_allocator();

    // A canonical folder ends with its separator
    char folder_buf[1024];
    const size_t folder_len = md_path_write_canonical(folder_buf, sizeof(folder_buf), MD_STR_LIT(VIAMD_TEST_DATA_DIR));
    ASSERT_GT(folder_len, 0u);
    const md_str_t folder = {folder_buf, folder_len};
    ASSERT_TRUE(md_str_ends_with(folder, MD_STR_LIT("/")));
    auto at = [&](const char* name) { return std::string(folder.ptr, folder.len) + name; };

    List l;
    l.record(file(at("catalyst.gro").c_str(), LoaderType_GRO), FileRole_System);
    l.record(file(at("catalyst.xtc").c_str(), LoaderType_XTC), FileRole_Trajectory);
    FileEntry lammps = file(at("sub dir/system.data").c_str(), LoaderType_LAMMPSDATA);    // not there, and with a space
    lammps.explicit_type = true;
    snprintf(lammps.arg, sizeof(lammps.arg), "atom-ID atom-type x y z");
    l.record(lammps, FileRole_Supplement);
    FileEntry gone = file(at("gone.edr").c_str(), LoaderType_EDR);
    gone.missing = true;
    l.record(gone, FileRole_Series);
    l.list.atom_count = 1234;
    l.list.runs[0].frame_count = 56;

    const std::string workspace = at("viamd_test.via"), data = at("sub dir/system.data");
    Written w(workspace.c_str(), alloc);
    EXPECT_EQ(8, file_list_serialize(w.s, &l.list, 3, 7));    // the next free run id
    const md_str_t text = w.text();

    // Ids as given, owners by id, paths relative to the workspace, the loader only where it was chosen
    EXPECT_TRUE(contains(text, "[System]\nId=3\nAtomCount=1234\n"));
    EXPECT_TRUE(contains(text, "[Run]\nId=7\nSystem=3\nFrameCount=56\n"));
    EXPECT_TRUE(contains(text, "[Input]\nSystem=3\nPath=./catalyst.gro\n"));
    EXPECT_TRUE(contains(text, "[Input]\nRun=7\nPath=./catalyst.xtc\n"));
    EXPECT_TRUE(contains(text, "[Input]\nSystem=3\nPath=./sub dir/system.data\nLoader=data\nArgument=atom-ID atom-type x y z\n"));
    EXPECT_TRUE(contains(text, "[Input]\nRun=7\nPath=./gone.edr\n"));
    EXPECT_FALSE(contains(text, "Loader=gro"));

    // Planned per owner: the system's own first, then its run - the topology no longer comes after the trajectory
    Read r(text, folder);
    { const std::string got_ = r.str(); EXPECT_STREQ("catalyst.gro:System:-1 system.data:Supplement:-1 catalyst.xtc:Trajectory:0 gone.edr:Series:0", got_.c_str()); }
    EXPECT_EQ(1234u, r.plan.atom_count);
    ASSERT_EQ(1u, md_array_size(r.plan.frame_counts));
    EXPECT_EQ(56u, r.plan.frame_counts[0]);
    EXPECT_STREQ(data.c_str(), r.plan.inputs[1].entry.path);   // a file that is not there keeps its path
    EXPECT_TRUE(r.plan.inputs[1].entry.explicit_type);
    EXPECT_EQ(LoaderType_LAMMPSDATA, r.plan.inputs[1].entry.type);
    EXPECT_STREQ("atom-ID atom-type x y z", r.plan.inputs[1].entry.arg);
    EXPECT_EQ(0u, r.plan.skipped_runs + r.plan.skipped_inputs);
}

UTEST(viamd_file_list, the_order_of_the_sections_does_not_matter) {
    const md_str_t text = MD_STR_LIT(
        "[Input]\nRun=5\nPath=/d/b.edr\n"
        "[Input]\nSystem=2\nPath=/d/a.top\n"
        "[Run]\nId=5\nSystem=2\n"
        "[Input]\nRun=5\nPath=/d/a.xtc\n"
        "[Input]\nSystem=2\nPath=/d/a.gro\n"
        "[System]\nId=2\n");
    Read r(text, MD_STR_LIT("/"));
    // Within an owner the order is the load order; the first of each is what builds or opens it
    { const std::string got_ = r.str(); EXPECT_STREQ("a.top:System:-1 a.gro:Supplement:-1 b.edr:Trajectory:0 a.xtc:Series:0", got_.c_str()); }
}

UTEST(viamd_file_list, each_system_is_planned_on_its_own) {
    const md_str_t text = MD_STR_LIT(
        "[System]\nId=0\n[System]\nId=1\n"
        "[Run]\nId=0\nSystem=0\n[Run]\nId=1\nSystem=0\n[Run]\nId=2\nSystem=1\n"
        "[Input]\nSystem=0\nPath=/d/a.gro\n"
        "[Input]\nRun=0\nPath=/d/rep1.xtc\n"
        "[Input]\nRun=1\nPath=/d/rep2.xtc\n"
        "[Input]\nSystem=1\nPath=/d/b.gro\n"
        "[Input]\nRun=2\nPath=/d/b.xtc\n"
        "[Input]\nRun=9\nPath=/d/dangling.xtc\n"
        "[Input]\nSystem=0\nRun=0\nPath=/d/both.xtc\n"
        "[Input]\nPath=/d/none.xtc\n");
    Read one(text, MD_STR_LIT("/"));
    { const std::string got_ = one.str(); EXPECT_STREQ("a.gro:System:-1 rep1.xtc:Trajectory:0", got_.c_str()); }
    EXPECT_EQ(1u, one.plan.skipped_runs);      // of the first system
    EXPECT_EQ(3u, one.plan.skipped_inputs);    // dangling, both owners, none

    // The second system, with its own files and run; what belongs nowhere was counted with the first
    Read second(text, MD_STR_LIT("/"), 1, 1);
    { const std::string got_ = second.str(); EXPECT_STREQ("b.gro:System:-1 b.xtc:Trajectory:0", got_.c_str()); }
    EXPECT_EQ(0u, second.plan.skipped_runs);
    EXPECT_EQ(0u, second.plan.skipped_inputs);

    // A system the workspace does not name has no files
    Read none(text, MD_STR_LIT("/"), 1, 7);
    EXPECT_EQ(0u, md_array_size(none.plan.inputs));

    // A build that opens several runs takes them in order
    Read two(text, MD_STR_LIT("/"), 4);
    { const std::string got_ = two.str(); EXPECT_STREQ("a.gro:System:-1 rep1.xtc:Trajectory:0 rep2.xtc:Trajectory:1", got_.c_str()); }
    EXPECT_EQ(0u, two.plan.skipped_runs);
}

UTEST(viamd_file_list, version_1_files_load_in_their_slots_order) {
    // Keys in any order: structure, trajectory, energies, series - and an empty slot is no file
    const md_str_t text = MD_STR_LIT(
        "[Files]\n"
        "SeriesFile=b.xvg\n"
        "TrajectoryFile=traj/a.xtc\n"
        "EnergyFile=a.edr\n"
        "MoleculeFile=a.gro\n"
        "SeriesFile=/abs/c.csv\n"
        "CoarseGrained=1\n");
    Read r(text, MD_STR_LIT("/ws/"));
    { const std::string got_ = r.str(); EXPECT_STREQ("a.gro:System:-1 a.xtc:Trajectory:0 a.edr:Series:0 b.xvg:Series:0 c.csv:Series:0", got_.c_str()); }
    const std::string traj = resolved("/ws/traj/a.xtc"), csv = resolved("/abs/c.csv");   // kept: the checks take pointers
    EXPECT_STREQ(traj.c_str(), r.plan.inputs[1].entry.path);
    EXPECT_STREQ(csv.c_str(), r.plan.inputs[4].entry.path);
    EXPECT_FALSE(r.plan.inputs[0].opens_own_run);   // it has a trajectory of its own

    // Without a trajectory the structure file's own frames were the run, and the series were along them
    Read own(MD_STR_LIT("[Files]\nMoleculeFile=a.pdb\nTrajectoryFile=\nSeriesFile=a.xvg\n"), MD_STR_LIT("/ws/"));
    { const std::string got_ = own.str(); EXPECT_STREQ("a.pdb:System:-1 a.xvg:Series:0", got_.c_str()); }
    EXPECT_TRUE(own.plan.inputs[0].opens_own_run);

    // Saved with nothing loaded
    Read empty(MD_STR_LIT("[Files]\nMoleculeFile=\nTrajectoryFile=\n"), MD_STR_LIT("/ws/"));
    EXPECT_EQ(0u, md_array_size(empty.plan.inputs));
}

UTEST(viamd_file_list, an_unknown_loader_falls_back_to_the_extension) {
    Read r(MD_STR_LIT("[System]\nId=0\n[Input]\nSystem=0\nPath=/x/file.h5\nLoader=nosuchformat\n"), MD_STR_LIT("/"));
    ASSERT_EQ(1u, md_array_size(r.plan.inputs));
    EXPECT_FALSE(r.plan.inputs[0].entry.explicit_type);
    EXPECT_EQ(LoaderType_Undefined, r.plan.inputs[0].entry.type);
}

// A workspace from before a topology was a system of its own names the coordinates with it: they are a file of the
// system after it, but for a trajectory, which gave its first frame and is the run already
UTEST(viamd_file_list, the_coordinates_of_a_topology_are_a_file_of_its_system) {
    // Files that are there, so that their paths resolve: a workspace folder with the topology and its coordinates
    const char* names[] = { "viamd_test_sys.prmtop", "viamd_test_md.nc", "viamd_test_eq.rst7" };
    for (const char* n : names) {
        md_file_t f = {};
        ASSERT_TRUE(md_file_open(&f, md_str_from_cstr(n), MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
        md_file_close(&f);
    }
    char folder_buf[1024];
    const size_t folder_len = md_path_write_canonical(folder_buf, sizeof(folder_buf), MD_STR_LIT("./"));
    ASSERT_GT(folder_len, 0u);
    std::string folder_s(folder_buf, folder_len);
    if (folder_s.back() != '/') folder_s += '/';
    const md_str_t folder = {folder_s.c_str(), folder_s.size()};

    Read restart(MD_STR_LIT("[System]\nId=0\n[Input]\nSystem=0\nPath=./viamd_test_sys.prmtop\nCoordinates=./viamd_test_eq.rst7\n"), folder);
    { const std::string got_ = restart.str(); EXPECT_STREQ("viamd_test_sys.prmtop:System:-1 viamd_test_eq.rst7:Supplement:-1", got_.c_str()); }

    Read trajectory(MD_STR_LIT("[System]\nId=0\n[Run]\nId=0\nSystem=0\n[Input]\nSystem=0\nPath=./viamd_test_sys.prmtop\nCoordinates=./viamd_test_md.nc\n"
                               "[Input]\nRun=0\nPath=./viamd_test_md.nc\n"), folder);
    { const std::string got_ = trajectory.str(); EXPECT_STREQ("viamd_test_sys.prmtop:System:-1 viamd_test_md.nc:Trajectory:0", got_.c_str()); }

    for (const char* n : names) remove(n);
}

UTEST(viamd_file_list, a_file_moved_out_takes_the_files_of_the_run_it_opens) {
    List l;
    l.record(file("a.prmtop"), FileRole_System);
    l.record(file("a.rst7"), FileRole_Supplement);
    l.record(file("a.nc"), FileRole_Trajectory);
    l.record(file("a.mdout"), FileRole_Series);
    List rest, moved;

    // Coordinates joined to a topology go alone, and the topology is what is left
    EXPECT_EQ(1u, file_list_split(&rest.list, &moved.list, &l.list, 1, l.alloc));
    { const std::string got_ = rest.str();  EXPECT_STREQ("a.prmtop:System:-1 a.nc:Trajectory:0 a.mdout:Series:0", got_.c_str()); }
    { const std::string got_ = moved.str(); EXPECT_STREQ("a.rst7:System:-1", got_.c_str()); }

    // A trajectory takes the series along it, and the run goes with them
    EXPECT_EQ(2u, file_list_split(&rest.list, &moved.list, &l.list, 2, l.alloc));
    { const std::string got_ = rest.str();  EXPECT_STREQ("a.prmtop:System:-1 a.rst7:Supplement:-1", got_.c_str()); }
    { const std::string got_ = moved.str(); EXPECT_STREQ("a.nc:Trajectory:0 a.mdout:Series:0", got_.c_str()); }
    EXPECT_EQ(0u, md_array_size(rest.list.runs));
    EXPECT_EQ(1u, md_array_size(moved.list.runs));

    // The file that built the system: what is left builds it
    EXPECT_EQ(1u, file_list_split(&rest.list, &moved.list, &l.list, 0, l.alloc));
    { const std::string got_ = rest.str();  EXPECT_STREQ("a.rst7:System:-1 a.nc:Trajectory:0 a.mdout:Series:0", got_.c_str()); }

    // A series alone is of no run where it goes
    EXPECT_EQ(1u, file_list_split(&rest.list, &moved.list, &l.list, 3, l.alloc));
    { const std::string got_ = rest.str();  EXPECT_STREQ("a.prmtop:System:-1 a.rst7:Supplement:-1 a.nc:Trajectory:0", got_.c_str()); }
    { const std::string got_ = moved.str(); EXPECT_STREQ("a.mdout:System:-1", got_.c_str()); }
}

UTEST(viamd_file_list, a_structure_of_several_frames_moved_out_takes_its_run) {
    List l;
    l.record(file("a.pdb"), FileRole_System);
    l.record(file("a.pdb"), FileRole_Trajectory);
    l.record(file("a.xvg"), FileRole_Series);
    l.record(file("a.ndx"), FileRole_Supplement);
    List rest, moved;
    EXPECT_EQ(3u, file_list_split(&rest.list, &moved.list, &l.list, 0, l.alloc));
    { const std::string got_ = rest.str();  EXPECT_STREQ("a.ndx:System:-1", got_.c_str()); }
    { const std::string got_ = moved.str(); EXPECT_STREQ("a.pdb:System:-1 a.pdb:Trajectory:0 a.xvg:Series:0", got_.c_str()); }
}

UTEST(viamd_file_list, runs_close_one_at_a_time) {
    List l;
    l.record(file("a.gro"), FileRole_System);
    l.record(file("r1.xtc"), FileRole_Trajectory);
    l.record(file("r1.edr"), FileRole_Series);
    l.record(file("r2.xtc"), FileRole_Trajectory);
    l.record(file("r2.edr"), FileRole_Series);
    // A series along a run named, the first, after the second was opened
    FileEntry along_first = file("r1.xvg");
    along_first.run = 0;
    l.record(along_first, FileRole_Series);
    { const std::string got_ = l.str(); EXPECT_STREQ("a.gro:System:-1 r1.xtc:Trajectory:0 r1.edr:Series:0 r2.xtc:Trajectory:1 r2.edr:Series:1 r1.xvg:Series:0", got_.c_str()); }
    EXPECT_EQ(0, file_list_run_of(&l.list, MD_STR_LIT("r1.xtc")));
    EXPECT_EQ(1, file_list_run_of(&l.list, MD_STR_LIT("r2.xtc")));
    EXPECT_EQ(-1, file_list_run_of(&l.list, MD_STR_LIT("r1.edr")));
    EXPECT_STREQ("r2.xtc", file_list_run_file_of(&l.list, 1)->path);

    // The first closed: its files go, and the second takes its place
    file_list_remove_run(&l.list, 0);
    { const std::string got_ = l.str(); EXPECT_STREQ("a.gro:System:-1 r2.xtc:Trajectory:0 r2.edr:Series:0", got_.c_str()); }
    EXPECT_EQ(1u, md_array_size(l.list.runs));
}
