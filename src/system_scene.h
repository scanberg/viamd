#pragma once

// THE SYSTEM IN THE SCENE
//
// The Representations window lists what is loaded by system: the system first, then the operations done to its
// coordinates before they are drawn, then the copies of it that are drawn (periodic images or an assembly), then its
// representations. This unit is the part of that which is the system's own (ApplicationState::system): its rows in the
// window, what the menus and the viewport do to it, and its section of a workspace.
//
// The operations are the coordinate operations of apply_state_operations and the tracking of bond changes. A system
// has each of them or not, and lists those it has in one fixed order, the order they apply in - none can be moved
// before another, as no other order is valid. One it has is turned on and off by its own flag, which keeps its
// settings while it is off. Removing it turns it off and puts its settings back to their defaults.
//
// What is done once, to the frame shown (centering, wrapping, making whole now, recomputing the bonds), is an action
// of the system, not an operation it has.

#include <stdint.h>
#include <stddef.h>

#include <core/md_str.h>

struct ApplicationState;

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

// Whether the system has the operation, and whether it is on
bool system_has_operation(const ApplicationState& app, SystemOperation op);
bool system_operation_enabled(const ApplicationState& app, SystemOperation op);

// Turns it on or off (it is had either way once turned on). Turned on, it is applied to the frame shown at once.
void system_set_operation_enabled(ApplicationState* app, SystemOperation op, bool enabled);

// Added on, or removed: turned off, its settings back to the defaults
void system_add_operation(ApplicationState* app, SystemOperation op);
void system_remove_operation(ApplicationState* app, SystemOperation op);

// What was turned on by other means than the window (the Operations menu, the viewport, a workspace) is had: once a
// frame
void system_operations_sync(ApplicationState* app);

// Applied once, to the frame shown, in the order that keeps them valid together (apply_state_operations). Nothing
// while an assembly is shown.
void system_apply_now(ApplicationState* app, bool center, bool wrap, bool make_whole);

// The assembly shown, -1 for none: its copies at once, and the view frames it
void system_set_assembly(ApplicationState* app, int assembly);

// The display name of the system: the file that built it
void system_label(char* buf, size_t cap, const ApplicationState& app);

// ## UI

// The rows of the system in the Representations window: its header (the system, shown or not, the operations to add
// and the actions), and when it is open, its operations and its copies. Returns whether it is open: its
// representations are listed under it then.
bool draw_system_rows(ApplicationState* app);

// The state of the tracking of bond changes in the frame shown, and its advanced settings: under the operation in the
// window, and in the Operations menu
void draw_bond_tracking_details(ApplicationState* app);

// ## Workspace

// [SystemScene], of the system with the id the files of the workspace give it (file_list_serialize)
void system_scene_serialize(viamd::serialization_state_t& state, const ApplicationState* app, int system_id);

// One entry of [SystemScene] or of the [Operations] of a workspace written before it, read into app. Returns false
// for an entry it does not read: System, which the caller matches against the system loaded, and what refers to the
// loaded system (Assembly, RecenterTarget), which the caller keeps until it is loaded.
bool system_scene_read_entry(ApplicationState* app, str_t ident, str_t arg);

// After the entries of a workspace written before [SystemScene] were read: it has the operations that are on, and
// those with settings of their own (a target, a scope)
void system_scene_from_legacy(ApplicationState* app);
