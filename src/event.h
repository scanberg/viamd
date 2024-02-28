#pragma once

#include <stdint.h>
#include <stddef.h>
#include <core/md_hash.h>

namespace viamd {
	
typedef uint32_t EventType;

// These are defined as compile time hashes as we want other components in the system to
// Define their own EventTypes that may be sent to the event system in a similar fashion
// There may be a collision at some point, but its not very likely

enum : EventType {
	EventType_ViamdInitialize 	   		= HASH_STR_LIT("VIAMD Initialize"),
	EventType_ViamdShutdown 	   		= HASH_STR_LIT("VIAMD Shutdown"),
	EventType_ViamdFrameTick 	   		= HASH_STR_LIT("VIAMD Frame Tick"),				// This is called once per frame
	EventType_ViamdDrawMenu				= HASH_STR_LIT("VIAMD Draw Menu"),				// Draw the menu

	EventType_ViamdSerialize			= HASH_STR_LIT("VIAMD Serialize Workspace"),
	EventType_ViamdDeserialize			= HASH_STR_LIT("VIAMD Deserialize Workspace"),

	EventType_ViamdTopologyInit			= HASH_STR_LIT("VIAMD Topology Initialize"),	// Called when topology is initialized
	EventType_ViamdTopologyFree			= HASH_STR_LIT("VIAMD Topology Free"),			// Called when topology is freed

	EventType_ViamdTrajectoryInit		= HASH_STR_LIT("VIAMD Trajectory Initialize"),	// Called when a trajectory is initialized
	EventType_ViamdTrajectoryFree		= HASH_STR_LIT("VIAMD Trajectory Free"),		// Called when a trajectory is freed

	//EventType_ViamdFileOpen 	   		= HASH_STR_LIT("VIAMD Open File"),

	//EventType_AnimationFrameChange		= HASH_STR_LIT("Animation Frame Changed"),

	EventType_HoverMaskChanged			= HASH_STR_LIT("Hover Mask Changed"),
	EventType_SelectionMaskChanged		= HASH_STR_LIT("Selection Mask Changed"),

	//EventType_RepresentationChange		= HASH_STR_LIT("Representation Changed"),

	//EventType_PropertyHovered 			= HASH_STR_LIT("Property Hovered"),
	//EventType_ScriptChange				= HASH_STR_LIT("Script Changed"),
	//EventType_ScriptEvalStart			= HASH_STR_LIT("Script Evaluation Started"),
	//EventType_ScriptEvalComplete		= HASH_STR_LIT("Script Evaluation Completed"),
};

struct Event {
	EventType	type;
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
void event_system_enqueue_event(EventType type, const void* payload = 0, uint64_t delay_in_ms = 0);

// This immediately broadcasts an event in the system
void event_system_broadcast_event(EventType type, const void* payload = 0);

// Call once per frame to process the queued up events
void event_system_process_event_queue();

}
