#pragma once

#include <event.h>

namespace viamd {

// These are defined as compile time hashes as we want other components in the system to
// Define their own EventTypes that may be sent to the event system in a similar fashion
// There may be a collision at some point, but its not very likely

enum : EventType {
	EventType_ViamdInitialize 	   			= HASH_STR_LIT("VIAMD Initialize"),
	EventType_ViamdShutdown 	   			= HASH_STR_LIT("VIAMD Shutdown"),
	EventType_ViamdFrameTick 	   			= HASH_STR_LIT("VIAMD Frame Tick"),				// This is called once per frame

	EventType_ViamdRenderOpaque				= HASH_STR_LIT("VIAMD Render Opaque"),			// Render opaque geometry into GBuffer
	EventType_ViamdRenderTransparent		= HASH_STR_LIT("VIAMD Render Transparent"),		// Render transparent geometry (Occurs after opaque render pass)
	// Before the main view is rendered, every frame: a handler that draws into it (RenderOpaque, RenderTransparent) folds
	// what it draws from into the key (payload RenderKey), as it changes without input (a hover, a result arriving). The
	// view is made again when the key changes, and else shown as it was made (lazy rendering).
	EventType_ViamdRenderKey				= HASH_STR_LIT("VIAMD Render Key"),

	// While the Windows menu is open, every frame: a component adds its windows to the request (payload
	// WindowMenuRequest) with window_menu_add. Where each is listed is its path, see window_menu.h
	EventType_ViamdWindowMenuRequest		= HASH_STR_LIT("VIAMD Window Menu Request"),

	// Workspace. Serialize: write your sections. Loading a workspace is three steps (see #workspace in
	// viamd.cpp): DeserializeBegin - reset what you keep in a workspace to its defaults, since a
	// section missing from the file means "the defaults"; Deserialize, once per section nobody else
	// reads - read it if it is yours; DeserializeEnd - the files are loaded, apply what needed them.
	// All three carry the DeserializationState (its filename identifies the workspace).
	EventType_ViamdSerialize				= HASH_STR_LIT("VIAMD Serialize Workspace"),
	EventType_ViamdDeserializeBegin			= HASH_STR_LIT("VIAMD Deserialize Workspace Begin"),
	EventType_ViamdDeserialize				= HASH_STR_LIT("VIAMD Deserialize Workspace"),
	EventType_ViamdDeserializeEnd			= HASH_STR_LIT("VIAMD Deserialize Workspace End"),

	EventType_ViamdSystemInit				= HASH_STR_LIT("VIAMD System Initialize"),		// Called when a system is initialized
	EventType_ViamdSystemFree				= HASH_STR_LIT("VIAMD System Free"),			// Called when a system is freed
	
    EventType_ViamdLoadData					= HASH_STR_LIT("VIAMD Load Data"),				// Called when a data is loaded

	EventType_ViamdTrajectoryInit			= HASH_STR_LIT("VIAMD Trajectory Initialize"),	// Called when a trajectory is initialized
	EventType_ViamdTrajectoryFree			= HASH_STR_LIT("VIAMD Trajectory Free"),		// Called when a trajectory is freed

	EventType_ViamdHighlightMaskChanged		= HASH_STR_LIT("VIAMD Highlight Mask Changed"),
	EventType_ViamdSelectionMaskChanged		= HASH_STR_LIT("VIAMD Selection Mask Changed"),
	EventType_ViamdRepresentationChanged	= HASH_STR_LIT("VIAMD Representation Changed"),		// Called when representations have been modified

    // System state corresponds to the transient portion of the system i.e. atom coordinates + unitcell, which can be modified by a script or by a trajectory frame change.
    EventType_ViamdSystemStateChanged		= HASH_STR_LIT("VIAMD System State Changed"),		// Called when system state has been modified (e.g. by a script or by a trajectory frame change)


	// Picking events
	// This reserve range is to give each event handler a chance to reserve a picking range before a frame is rendered.
	EventType_ViamdPickingRangeReserve				= HASH_STR_LIT("VIAMD Picking Range Reserve"),
	EventType_ViamdPickingHit						= HASH_STR_LIT("VIAMD Picking Hit"),
	EventType_ViamdPickingTooltipTextRequest		= HASH_STR_LIT("VIAMD Picking Tooltip Text Request"),

	// InteractionSurface
    EventType_ViamdInteractionSurface				= HASH_STR_LIT("VIAMD Interaction Surface"),
    EventType_ViamdViewFit							= HASH_STR_LIT("VIAMD View Fit"),

	// The script's evaluation, state.script.eval (payload ApplicationState).
	//   Init   one has been made, taking over what an earlier one held
	//   Free   it is about to be freed, and is still state.script.eval: stop and wait for whatever reads it on other
	//          threads, and drop what points into it
	//   Begin  it is about to be run over frames. Until End its attributes change (script::attributes) and are not
	//          read; what it writes per frame for any thread to read goes on being read (script::eval_pairs,
	//          script::eval_output_frames)
	//   End    the run has ended, finished or interrupted: its attributes are what it published
	EventType_ViamdScriptEvalInit					= HASH_STR_LIT("VIAMD Script Eval Initialize"),
	EventType_ViamdScriptEvalFree					= HASH_STR_LIT("VIAMD Script Eval Free"),
	EventType_ViamdScriptEvalBegin					= HASH_STR_LIT("VIAMD Script Eval Begin"),
	EventType_ViamdScriptEvalEnd					= HASH_STR_LIT("VIAMD Script Eval End"),

	// Every frame, for the task overlay (draw_async_task_window in main.cpp): a handler with long work under way that the
	// task system does not show adds a line to the request (payload TaskProgressRequest)
	EventType_ViamdTaskProgressRequest				= HASH_STR_LIT("VIAMD Task Progress Request"),
};

enum : EventPayloadType {
	EventPayloadType_Undefined					= 0,
	EventPayloadType_ApplicationState			= HASH_STR_LIT("Payload Application State"),
	EventPayloadType_Representation				= HASH_STR_LIT("Payload Representation"),
	EventPayloadType_SerializationState			= HASH_STR_LIT("Payload Serialization State"),
	EventPayloadType_DeserializationState		= HASH_STR_LIT("Payload Deserialization State"),
	EventPayloadType_PickingSpace				= HASH_STR_LIT("Payload Picking Space"),
	EventPayloadType_PickingHit					= HASH_STR_LIT("Payload Picking Hit"),
	EventPayloadType_PickingTooltipTextRequest	= HASH_STR_LIT("Payload Picking Tooltip Text Request"),
    EventPayloadType_LoadData					= HASH_STR_LIT("Payload Load Data"),
    EventPayloadType_InteractionSurfaceEvent    = HASH_STR_LIT("Payload Interaction Surface Event"),
	EventPayloadType_ViewFitRequest				= HASH_STR_LIT("Payload View Fit Request"),
	EventPayloadType_RenderKey					= HASH_STR_LIT("Payload Render Key"),
	EventPayloadType_TaskProgressRequest		= HASH_STR_LIT("Payload Task Progress Request"),
	EventPayloadType_WindowMenuRequest			= HASH_STR_LIT("Payload Window Menu Request"),
};

}
