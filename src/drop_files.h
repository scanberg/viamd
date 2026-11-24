#pragma once

#include <file_placement.h>
#include <file_list.h>
#include <loader.h>

#include <md_system.h>

// #drop_files
//
// The files of one drop (or one Open), matched against each other and against the systems loaded, for file_placement:
// what each file is (drop_file_kind), and whether it is of a system's atoms, which mdlib answers - md_system_match for a
// system of its own, the format's check for data of a system's atoms and for a series along a run.
//
// A file which is a system of its own is loaded when it is first matched, on an arena of its own, and kept: it is the
// system the file is, which loading it then takes (load_data's loaded) rather than reading the file again. A trajectory
// is the system of its atoms by their count, at its first frame; matched for a series of the drop, its frames are
// published as a run of that system too.
//
// Free of the application, so that what matches what is tested on real files (test/test_drop_files.cpp).

// A system loaded, as a drop sees it
struct DropSystem {
    uint32_t                 id = PLACEMENT_NONE;
    const md_system_t*       sys = nullptr;
    const md_system_state_t* state = nullptr;
    md_str_t                 run = {};          // "run/<name>" of its run, empty for none
};

struct DropFile {
    FileEntry           entry;
    loader::LoaderState ls;
    PlacementKind       kind = PlacementKind_System;
    uint32_t            target = PLACEMENT_AUTO;

    // The system it is, loaded when first matched: 1 loaded, -1 it would not load, 0 not tried
    int                 loaded = 0;
    md_allocator_i*     arena = nullptr;
    md_system_t         sys = {};
    md_system_state_t   state = {};
    bool                run_published = false;  // A trajectory's frames as the run DROP_RUN of sys
};

struct DropFiles {
    DropFile          files[PLACEMENT_MAX_FILES];
    size_t            num_files = 0;
    size_t            num_skipped = 0;          // Past PLACEMENT_MAX_FILES
    const DropSystem* loaded = nullptr;
    size_t            num_loaded = 0;
};

// What a file of these loader flags is to placement
PlacementKind drop_file_kind(LoaderFlags flags);

// The systems loaded stay as they are while the drop is planned (pointers into them are kept)
void drop_files_init(DropFiles* drop, const DropSystem* loaded, size_t num_loaded);

// A file of the drop, its loader told; target a loaded system's id, PLACEMENT_AUTO or PLACEMENT_NEW
void drop_files_add(DropFiles* drop, md_str_t path, const loader::LoaderState& ls, uint32_t target = PLACEMENT_AUTO);

// Whether file is of onto's atoms (PlacementQuery)
PlacementMatch drop_files_match(DropFiles* drop, size_t file, PlacementRef onto);

// The plan of the drop, at most max_systems systems in all (file_placement_plan)
void drop_files_plan(PlacementPlan* plan, DropFiles* drop, size_t max_systems);

// The system the file is, loaded now if it was not: NULL when it does not load, or is no system of its own
const DropFile* drop_files_system(DropFiles* drop, size_t file);

void drop_files_free(DropFiles* drop);
