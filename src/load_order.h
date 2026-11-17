#pragma once

#include <file_list.h>
#include <loader.h>

#include <core/md_str.h>

// #load_order
//
// A file loaded on its own (dropped, or opened from the menu) takes the role its loader gives it - except that a
// topology without coordinates (an Amber prmtop, a PSF) and the coordinates for it (a restart, NAMD binary
// coordinates, a structure file of the same atoms, a trajectory) are two halves of one system, and come in either
// order. Whichever comes first waits (PendingFiles); the second builds the system:
//
//   topology, then restart or trajectory    the system from the topology with those coordinates; a trajectory then
//   restart or trajectory, then topology    opens as the run too
//   structure file and topology, any order  the structure builds the system, the topology fills it in
//   new coordinates                         a system built from a topology is built again with them, its run opened
//                                           again
//
// A topology of as many atoms as the loaded system fills it in, and a trajectory of other atoms than the loaded
// system's (as its header tells) waits for its topology rather than failing against the system.
//
// Free of the application, so that the orders are tested on their own (test/test_load_order.cpp): what loading does
// to the application is handed in (LoadOps).

// What load_data_from_file made of a file
enum LoadResult {
    LoadResult_Failed = 0,
    LoadResult_System,      // A new system, with its run when the file holds several frames
    LoadResult_Trajectory,  // A new run for the loaded system
    LoadResult_Supplement,  // Added to the loaded system, or along its run: the system and the run are as they were
    LoadResult_Pending,     // Half of a system (a topology without coordinates, or coordinates without one): kept until the other half is loaded
};

// Halves of a system waiting for the other half. Empty paths when none waits.
struct PendingFiles {
    FileEntry topology;
    FileEntry coordinates;
    FileEntry trajectory;
};

// What loading does, in the application or in a test
struct LoadOps {
    void* user = nullptr;
    // Loads the file in the role (load_data). A topology with load_state.coordinates set builds the system with them.
    LoadResult (*load)(void* user, md_str_t path, const loader::LoaderState& load_state, FileRole role, bool opens_own_run) = nullptr;
    size_t (*system_atoms)(void* user) = nullptr;   // 0 without a system
    size_t (*run_frames)(void* user) = nullptr;     // 0 without a run
    const FileList* (*files)(void* user) = nullptr; // the files the loaded system was built from
    // The atom count of a file as its header tells (loader::atom_count when not given)
    size_t (*atom_count)(void* user, md_str_t path, const loader::LoaderState& load_state) = nullptr;
};

// The loader state of an entry, as it was loaded
loader::LoaderState load_order_loader_state(const FileEntry& e);

// Loads a file on its own, in the order rules above
LoadResult load_in_order(PendingFiles* pending, const LoadOps& ops, md_str_t path, const loader::LoaderState& load_state);
