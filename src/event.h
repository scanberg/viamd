#pragma once

#include <stdint.h>
#include <stddef.h>
#include <core/md_hash.h>

namespace viamd {
	
typedef uint32_t EventType;
typedef uint32_t EventPayloadType;

// These are defined as compile time hashes as we want other components in the system to
// Define their own EventTypes that may be sent to the event system in a similar fashion
// There may be a collision at some point, but its not very likely

enum : EventType {
	EventType_ViamdInitialize 	   		= HASH_STR_LIT("VIAMD Initialize"),
	EventType_ViamdShutdown 	   		= HASH_STR_LIT("VIAMD Shutdown"),
	EventType_ViamdFrameTick 	   		= HASH_STR_LIT("VIAMD Frame Tick"),				// This is called once per frame

	EventType_ViamdRenderOpaque			= HASH_STR_LIT("VIAMD Render Opaque"),			// Render opaque geometry into GBuffer
	EventType_ViamdRenderTransparent	= HASH_STR_LIT("VIAMD Render Transparent"),		// Render transparent geometry (Occurs after opaque render pass)

	EventType_ViamdWindowDrawMenu		= HASH_STR_LIT("VIAMD Draw Window Menu"),		// Draw the menu

	EventType_ViamdSerialize			= HASH_STR_LIT("VIAMD Serialize Workspace"),
	EventType_ViamdDeserialize			= HASH_STR_LIT("VIAMD Deserialize Workspace"),

	EventType_ViamdTopologyInit			= HASH_STR_LIT("VIAMD Topology Initialize"),	// Called when topology is initialized
	EventType_ViamdTopologyFree			= HASH_STR_LIT("VIAMD Topology Free"),			// Called when topology is freed

	EventType_ViamdTrajectoryInit		= HASH_STR_LIT("VIAMD Trajectory Initialize"),	// Called when a trajectory is initialized
	EventType_ViamdTrajectoryFree		= HASH_STR_LIT("VIAMD Trajectory Free"),		// Called when a trajectory is freed

	//EventType_ViamdFileOpen 	   		= HASH_STR_LIT("VIAMD Open File"),

	//EventType_AnimationFrameChange		= HASH_STR_LIT("Animation Frame Changed"),

	EventType_ViamdHoverMaskChanged			= HASH_STR_LIT("VIAMD Hover Mask Changed"),
	EventType_ViamdSelectionMaskChanged		= HASH_STR_LIT("VIAMD Selection Mask Changed"),
	EventType_ViamdRepresentationChanged	= HASH_STR_LIT("VIAMD Representation Changed"),		// Called when representations have been modified

    // System state corresponds to the transient portion of the system i.e. atom coordinates + unitcell, which can be modified by a script or by a trajectory frame change.
    EventType_ViamdSystemStateChanged		= HASH_STR_LIT("VIAMD System State Changed"),		// Called when system state has been modified (e.g. by a script or by a trajectory frame change)

	EventType_ViamdRepresentationInfoFill					= HASH_STR_LIT("VIAMD Representation Info Fill"),
	EventType_ViamdRepresentationEvalElectronicStructure	= HASH_STR_LIT("VIAMD Representation Eval ElectronicStructure"),
	EventType_ViamdRepresentationEvalAtomProperty			= HASH_STR_LIT("VIAMD Representation Eval AtomProperty"),

	//EventType_PropertyHovered 			= HASH_STR_LIT("Property Hovered"),
	//EventType_ScriptChange				= HASH_STR_LIT("Script Changed"),
	//EventType_ScriptEvalStart			= HASH_STR_LIT("Script Evaluation Started"),
	//EventType_ScriptEvalComplete		= HASH_STR_LIT("Script Evaluation Completed"),
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
};

struct Event {
	EventType	type;
	EventPayloadType payload_type;
	uint64_t    timestamp;
	const void* payload;
};

// This represents the interface to implement if you want to define a component that can process events
struct EventHandler {
    virtual void process_events(const Event* events, size_t num_events) = 0;
};

void event_system_register_handler(EventHandler& event_handler);

// Queues up an event to be processed (Prefer this, unless the event has to be processed now)
// Events are processed in batches by each registered event
void event_system_enqueue_event(EventType type, EventPayloadType payload_type = EventPayloadType_Undefined, const void* payload = 0, uint64_t delay_in_ms = 0);

// This immediately broadcasts an event in the system
void event_system_broadcast_event(EventType type, EventPayloadType payload_type = EventPayloadType_Undefined, const void* payload = 0);

// Call once per frame to process the queued up events
void event_system_process_event_queue();

}
