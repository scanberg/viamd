#pragma once

// THE SYSTEM IN THE SCENE
//
// The Scene window lists what is loaded by system: the system first, then its files, the operations done to its
// coordinates before they are drawn, the copies of it that are drawn (periodic images or an assembly), then its
// representations. This unit is the part of that which is the system's own (ApplicationState::system): its rows in the
// window, what the menus and the viewport do to it, and its section of a workspace.
//
// The operations are the coordinate operations of apply_state_operations and the tracking of bond changes. Every
// system lists all of them, opted into one by one, in one fixed order, the order they apply in - none can be moved
// before another, as no other order is valid. Turned off, an operation keeps its settings.
//
// What is done once, to the frame shown (centering, wrapping, making whole now, recomputing the bonds), is an action
// of the system, not an operation it has.

#include <stdint.h>
#include <stddef.h>

#include <core/md_str.h>
#include <file_list.h>

struct ApplicationState;
struct SystemContext;

namespace viamd {
struct serialization_state_t;
struct deserialization_state_t;
}

enum SystemOperation : uint32_t {
    SystemOperation_Center     = 1u << 0,   // operations.recenter: its target in the middle of the cell, its orientation kept or not
    SystemOperation_Wrap       = 1u << 1,   // operations.apply_pbc: every atom into the cell
    SystemOperation_MakeWhole  = 1u << 2,   // operations.unwrap_structures: the molecules split across the cell joined
    SystemOperation_TrackBonds = 1u << 3,   // operations.bond_tracking: the bonds formed and broken in each frame
};
constexpr int SYSTEM_OPERATION_COUNT = 4;

// The operations in the order they are listed and applied, by index
SystemOperation system_operation_at(int idx);
const char* system_operation_name(SystemOperation op);

// Whether the operation is on
bool system_operation_enabled(const SystemContext& sys, SystemOperation op);

// Turns it on or off. Turned on, it is applied to the frame shown at once.
void system_set_operation_enabled(SystemContext& sys, SystemOperation op, bool enabled);

// Every operation off, its settings back to the defaults: as a new system has them
void system_reset_operations(SystemContext& sys);

// Applied once, to the frame shown, in the order that keeps them valid together (apply_state_operations). Nothing
// while an assembly is shown.
void system_apply_now(SystemContext& sys, bool center, bool wrap, bool make_whole);

// The assembly shown, -1 for none: its copies at once, and the view frames it
void system_set_assembly(SystemContext& sys, int assembly);

// The display name of the system: the file that built it
void system_label(char* buf, size_t cap, const SystemContext& sys);

// ## UI

// The header of the system in the Scene window: the system, shown or not, its menu of the operations to add, the
// actions, and removing it, which is asked first and then left to the caller (*remove), once the window is drawn.
// Returns whether it is open: what is the system's is listed under it then, indented, by the caller.
bool draw_system_header(SystemContext& sys, bool* remove);

// The files of the system, one line each: what it is, its name and what it holds in short. The System window has them
// in full (its Files), and this links to it.
void draw_system_files(SystemContext& sys);

// The system's operations, in their order, folded under one row that says which are on, then what it is shown as (its
// copies): under its files
void draw_system_operations(SystemContext& sys);

// Unloads the system: its files and data, its representations, its operations and what it is shown as
void system_remove(SystemContext& sys);

// ## Files

// What a file is to its system ("Structure", "Trajectory", "Energies", ...), and what it holds in short ("500 frames,
// 0 - 4990 ps"): the Scene window's files and the System window's
const char* system_file_role_label(const FileEntry& e, FileRole role);
void system_file_details(char* buf, size_t cap, const SystemContext& sys, size_t i);

// ## Workspace

// [SystemScene] of the system, by its id, which its [System] has too (file_list_serialize)
void system_scene_serialize(viamd::serialization_state_t& state, const SystemContext& sys);

// One entry of [SystemScene] or of the [Operations] of a workspace written before it, read into app. Returns false
// for an entry it does not read: System, which the caller matches against the system loaded, and what refers to the
// loaded system (Assembly, RecenterTarget), which the caller keeps until it is loaded.
bool system_scene_read_entry(SystemContext& sys, str_t ident, str_t arg);
