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

	EventType_ViamdWindowDrawMenu			= HASH_STR_LIT("VIAMD Draw Window Menu"),		// Draw the menu

	EventType_ViamdSerialize				= HASH_STR_LIT("VIAMD Serialize Workspace"),
	EventType_ViamdDeserialize				= HASH_STR_LIT("VIAMD Deserialize Workspace"),

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

	EventType_ViamdRepresentationInfoFill					= HASH_STR_LIT("VIAMD Representation Info Fill"),
	EventType_ViamdRepresentationEvalElectronicStructure	= HASH_STR_LIT("VIAMD Representation Eval ElectronicStructure"),
	EventType_ViamdRepresentationEvalAtomProperty			= HASH_STR_LIT("VIAMD Representation Eval AtomProperty"),

	// Picking events
	// This reserve range is to give each event handler a chance to reserve a picking range before a frame is rendered.
	EventType_ViamdPickingRangeReserve				= HASH_STR_LIT("VIAMD Picking Range Reserve"),
	EventType_ViamdPickingHit						= HASH_STR_LIT("VIAMD Picking Hit"),
	EventType_ViamdPickingTooltipTextRequest		= HASH_STR_LIT("VIAMD Picking Tooltip Text Request"),

	// InteractionSurface
    EventType_ViamdInteractionSurface				= HASH_STR_LIT("VIAMD Interaction Surface"),
};

enum : EventPayloadType {
	EventPayloadType_Undefined					= 0,
	EventPayloadType_ApplicationState			= HASH_STR_LIT("Payload Application State"),
	EventPayloadType_RepresentationInfo			= HASH_STR_LIT("Payload Representation Info"),
	EventPayloadType_Representation				= HASH_STR_LIT("Payload Representation"),
	EventPayloadType_SerializationState			= HASH_STR_LIT("Payload Serialization State"),
	EventPayloadType_DeserializationState		= HASH_STR_LIT("Payload Deserialization State"),
	EventPayloadType_EvalElectronicStructure	= HASH_STR_LIT("Payload Eval ElectronicStructure"),
	EventPayloadType_EvalAtomProperty			= HASH_STR_LIT("Payload Eval AtomProperty"),
	EventPayloadType_PickingSpace				= HASH_STR_LIT("Payload Picking Space"),
	EventPayloadType_PickingHit					= HASH_STR_LIT("Payload Picking Hit"),
	EventPayloadType_PickingTooltipTextRequest	= HASH_STR_LIT("Payload Picking Tooltip Text Request"),
    EventPayloadType_LoadData					= HASH_STR_LIT("Payload Load Data"),
    EventPayloadType_InteractionSurfaceEvent    = HASH_STR_LIT("Payload Interaction Surface Event"),
};

}
