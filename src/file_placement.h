#pragma once

#include <stddef.h>
#include <stdint.h>

// #file_placement
//
// Where the files of one drop (or one Open) go: into which loaded system, into a new one, or nowhere, and in what order
// each system loads them. A pure function of what the files are, the systems loaded, and whether a file is of a
// system's atoms - which the caller answers (mdlib's md_system_match for a system of its own, a header's atom count for
// a trajectory, time ranges for a series) - so that the rules are tested on their own (test/test_file_placement.cpp).
//
// THE RULES
//   Matching      A file matches a system when it is of its atoms. It is a candidate of the system when it matches and
//                 replaces nothing there (its coordinates over those the system has, complete bonds over complete
//                 bonds, an attribute of the same path): PlacementMatch_Join. PlacementMatch_Override matches but would
//                 replace.
//   Groups        The files of one drop are matched against each other before the loaded systems, in three passes by
//                 what depends on what: the files which are systems of their own (topologies first, so that their
//                 names stand, then the other systems, then trajectories) form the systems the drop makes; then the
//                 data of a system's atoms (an index file); then the series along a run.
//   Placement     A system the drop makes joins the one loaded system it is a candidate of. Of none, it is a new
//                 system; of several, the user picks one of them or a new system (PlacementAction_Ask). One that
//                 matches a loaded system but would replace something there is never joined to it: it is a new system
//                 made as a copy of that one, the drop's files applied after (PlacementAction_Copy). Placement never
//                 replaces anything of a loaded system.
//   Alone         What cannot stand alone - data of a system's atoms, a series - is refused when it has no candidate.
//                 A trajectory alone is a system of points.
//   Explicit      A file dropped onto a system (target) goes into it if it matches it, replacing what it does, and is
//                 refused otherwise; one dropped onto the space below the systems (PLACEMENT_NEW) makes a new system,
//                 whatever it matches.
//   The cap       A new system past max_systems is refused; a choice offers a new system only while there is room.

constexpr uint32_t PLACEMENT_NONE = UINT32_MAX;         // No system
constexpr uint32_t PLACEMENT_AUTO = UINT32_MAX - 1;     // A file's target: placed by the rules
constexpr uint32_t PLACEMENT_NEW  = UINT32_MAX - 2;     // A file's target: a new system

constexpr size_t PLACEMENT_MAX_FILES      = 32;         // Of one drop
constexpr size_t PLACEMENT_MAX_CANDIDATES = 16;

// What a file is to placement
enum PlacementKind {
    PlacementKind_System,       // A system of its own: atoms (a structure, a restart, QM data, an itp's molecule types)
    PlacementKind_Topology,     // A system of its own whose names stand over a structure file's (prmtop, psf, top, tpr)
    PlacementKind_Trajectory,   // Frames: a run of a system, or a system of points alone
    PlacementKind_Data,         // Data of a system's atoms, which cannot stand alone (an index file)
    PlacementKind_Series,       // Data along a run, which cannot stand alone (edr, xvg, csv, Amber energies)
};

enum PlacementMatch {
    PlacementMatch_No = 0,      // Not of its atoms (or, for a series, not along that run)
    PlacementMatch_Join,        // Of its atoms, and replaces nothing there
    PlacementMatch_Override,    // Of its atoms, but replaces something there
};

// What a file is matched onto: a loaded system (and a run of it, for a series), or the system or run a file of the
// drop makes
struct PlacementRef {
    uint32_t system = PLACEMENT_NONE;   // A loaded system's id
    int32_t  run    = -1;               // Its run, for a series
    int32_t  file   = -1;               // Else a file of the drop: the system it builds, or the run a trajectory opens
};

struct PlacementFile {
    PlacementKind kind   = PlacementKind_System;
    uint32_t      target = PLACEMENT_AUTO;      // A loaded system's id, PLACEMENT_AUTO or PLACEMENT_NEW
};

struct PlacementLoaded {
    uint32_t id = PLACEMENT_NONE;
    size_t   num_runs = 0;
};

// Whether file (an index into the drop) is of onto's atoms. For a series onto is a run.
struct PlacementQuery {
    void* user = nullptr;
    PlacementMatch (*match)(void* user, size_t file, PlacementRef onto) = nullptr;
};

enum PlacementAction {
    PlacementAction_Join,       // Into the loaded system 'system'
    PlacementAction_New,        // A new system
    PlacementAction_Copy,       // A new system, a copy of the loaded system 'system' (its files but its runs and their
                                // series), the files applied after
    PlacementAction_Ask,        // The user picks among candidates (and a new system where new_allowed)
};

enum PlacementReason {
    PlacementReason_NeedsSystem,    // Data of a system's atoms, of no system
    PlacementReason_NeedsRun,       // A series, along no run
    PlacementReason_NotOfTarget,    // Dropped onto a system whose atoms it is not of
    PlacementReason_Cap,            // A new system, past max_systems
    PlacementReason_TooMany,        // Past PLACEMENT_MAX_FILES
};

// A file in the order its system loads it
struct PlacementEntry {
    size_t       file = 0;
    PlacementRef run;           // Its run: a trajectory opens one (file is itself); a series joins a loaded system's run
                                // or the run a trajectory of the drop opens. Neither: of the system.
};

struct PlacementSystem {
    PlacementAction action = PlacementAction_New;
    uint32_t        system = PLACEMENT_NONE;    // Join, Copy
    PlacementRef    candidates[PLACEMENT_MAX_CANDIDATES];  // Ask: loaded systems, the systems the drop makes, or runs
    size_t          num_candidates = 0;
    bool            new_allowed = false;        // Ask: a new system is offered too
    PlacementEntry  entries[PLACEMENT_MAX_FILES];
    size_t          num_entries = 0;
};

struct PlacementRefusal {
    size_t          file = 0;
    PlacementReason reason = PlacementReason_NeedsSystem;
};

struct PlacementPlan {
    PlacementSystem  systems[PLACEMENT_MAX_FILES];
    size_t           num_systems = 0;
    PlacementRefusal refused[PLACEMENT_MAX_FILES];
    size_t           num_refused = 0;
};

// The plan for the files of one drop onto the systems loaded, at most max_systems of them in all
void file_placement_plan(PlacementPlan* plan, const PlacementFile* files, size_t num_files, const PlacementLoaded* loaded, size_t num_loaded, size_t max_systems, PlacementQuery query);
