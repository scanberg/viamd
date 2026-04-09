#pragma once

#include <stdint.h>
#include <stddef.h>
#include <core/md_hash.h>

namespace viamd {
	
typedef uint32_t EventType;
typedef uint32_t EventPayloadType;

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

// Register an event handler to receive events submitted to the system
void event_system_register_handler(EventHandler& event_handler);

// Queues up an event to be processed (Prefer this, unless the event has to be processed now)
// Events are processed in batches by each registered event handler
void event_system_enqueue_event(EventType type, EventPayloadType payload_type = 0, const void* payload = 0, uint64_t delay_in_ms = 0);

// This immediately broadcasts an event in the system
void event_system_broadcast_event(EventType type, EventPayloadType payload_type = 0, const void* payload = 0);

// Call once per frame to process the queued up events
void event_system_process_event_queue();

}
