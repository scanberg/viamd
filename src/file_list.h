#pragma once

#include <loader.h>

#include <core/md_str.h>
#include <core/md_array.h>

struct md_allocator_i;

namespace viamd {
struct serialization_state_t;
struct deserialization_state_t;
}

// #file_list
//
// The files a system is built from. A system is an ordered list of files: the first builds it, and each
// one after fills it in (a topology, QM data for its atoms). A run belongs to a system and is an ordered list
// of files too: the first opens it with its frames, and each one after joins it (energies, xvg and csv
// series). Nothing is typed by slot - there is no "the molecule file" or "the energy file" - only which
// system or run a file belongs to, and its place there. That is also where its data goes in mdlib: into the
// system's attribute table, or under the run's prefix.
//
// How a file is read follows from its place (file_list_role), not from what is loaded when it is read, so
// loading the same files into the same places builds the same system.
//
// In a workspace (version 2):
//
//     [System]            Id, AtomCount
//     [Run]               Id, System, FrameCount
//     [Input]             System or Run, Path, Loader (when chosen), Argument (when any)
//
// Ids are local to the file: written densely at each save, and only resolved within it. Every section which
// refers to a system or a run does so by these ids. AtomCount and FrameCount tell whether the files are still
// the ones the workspace was saved with: what is stored by atom index or by frame is not applied to others.
// The files load per system - its own in order, then each of its runs, each in order - so a run is opened
// on the complete system.
//
// Free of the application on purpose, so that the rules are tested on their own (test/test_file_list.cpp).

// What a file is to its system
enum FileRole {
    FileRole_None = 0,
    FileRole_System,      // Builds the system: the first file of a system
    FileRole_Supplement,  // Fills in the system: every later file of the system
    FileRole_Trajectory,  // Opens a run with its frames: the first file of a run
    FileRole_Series,      // Joins a run, sampled along it: every later file of a run
};

// The role a file loaded on its own takes (dropped, or opened from the menu), from its loader flags: what can
// supplement the loaded system does (a VeloxChem file of the same atoms), else what can build one does, and a
// trajectory only otherwise.
FileRole file_role(LoaderFlags flags);

struct FileEntry {
    char       path[1024] = "";                 // Canonical
    LoaderType type = LoaderType_Undefined;     // The loader the file was read with
    bool       explicit_type = false;           // The loader was chosen by the user rather than told from the file
    char       arg[256] = "";                   // The loader's argument (the LAMMPS atom format), empty when none
    int        run = -1;                        // The run it belongs to (index into FileList::runs), -1 for the system's own
    bool       missing = false;                 // Named by a workspace but not loaded: kept, so a save does not drop it
};

struct RunFiles {
    size_t frame_count = 0;                     // Of the run, as loaded or as the workspace stored it
};

struct FileList {
    md_array(FileEntry) entries = 0;            // In load order: the system's own and those of its runs
    md_array(RunFiles)  runs = 0;               // In the order they were opened
    size_t atom_count = 0;                      // Of the system, as loaded or as the workspace stored it
};

static inline size_t           file_list_count(const FileList* list)            { return md_array_size(list->entries); }
static inline const FileEntry* file_list_get(const FileList* list, size_t i)    { return i < md_array_size(list->entries) ? &list->entries[i] : NULL; }

void file_list_free (FileList* list, struct md_allocator_i* alloc);
void file_list_clear(FileList* list);

// The role of entry i, from its place: first of the system, later in the system, first of its run, later in it
FileRole file_list_role(const FileList* list, size_t i);

// A file was loaded in this role, or was named by a workspace and could not be (entry.missing). Building
// starts the list over. A trajectory opens a new run, which the file then belongs to. A series joins the last
// run, and is not recorded without one. A supplement or a series identical to one already in its place (path,
// loader and argument) is not added twice. Returns false when nothing was recorded.
bool file_list_record(FileList* list, FileEntry entry, FileRole role, struct md_allocator_i* alloc);

// The runs were closed: their files go with them
void file_list_drop_runs(FileList* list);

// The file that built the system, NULL when there is none
const FileEntry* file_list_system_file(const FileList* list);

// The file that opened the last run, or else the one that built the system, NULL when there is neither
const FileEntry* file_list_run_file(const FileList* list);

// ## Workspace

// Writes the [System] section of the list with the given id, a [Run] for each of its runs with ids from
// first_run_id, and an [Input] per file in load order. Paths are relative to the workspace (state.filename) so
// that a dataset can be moved as a whole, and absolute where no relative path exists (another volume).
// Returns the next free run id.
int file_list_serialize(viamd::serialization_state_t& state, const FileList* list, int system_id, int first_run_id);

// What a workspace says about its systems and their files, read section by section and resolved once all is read
struct WorkspaceFiles {
    struct System { int id; size_t atom_count; };
    struct Run    { int id; int system; size_t frame_count; };
    struct Input  {
        FileEntry entry;
        int       system;           // Ids as written, -1 where not given
        int       run;
        FileRole  role;             // None: from its place. Given for a [Files] section, which had slots instead
        bool      opens_own_run;    // See SystemPlan::Input
    };
    md_array(System) systems = 0;
    md_array(Run)    runs = 0;
    md_array(Input)  inputs = 0;
    bool legacy = false;    // From a [Files] section: written before systems and runs
};

// Reads the section the state is at if it is one of [System], [Run], [Input], or [Files] (a workspace written
// before version 2), and returns whether it was. A relative path is resolved against folder, the workspace's.
bool workspace_files_read(WorkspaceFiles* files, viamd::deserialization_state_t& state, str_t folder, struct md_allocator_i* alloc);

// One system as it loads: its files in order with the role each is read in
struct SystemPlan {
    struct Input {
        FileEntry entry;
        FileRole  role;
        int       run;              // The run's index within the system, -1 for the system's own
        bool      opens_own_run;    // The system's file opens its run too when it holds several frames (legacy)
    };
    md_array(Input)  inputs = 0;
    size_t           atom_count = 0;            // As stored, 0 when not
    md_array(size_t) frame_counts = 0;          // Per run, as stored, 0 when not
    size_t skipped_systems = 0;                 // Named, but not loaded: more than this build loads
    size_t skipped_runs = 0;
    size_t skipped_inputs = 0;                  // Belonging to no system or run the workspace names
};

// The plan for the first system the workspace names, with at most max_runs of its runs. What is left out is
// counted in the plan, so that it can be reported.
void workspace_files_plan(SystemPlan* plan, const WorkspaceFiles* files, size_t max_runs, struct md_allocator_i* alloc);
