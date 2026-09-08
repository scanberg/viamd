#include "event.h"
#include <core/md_array.h>
#include <core/md_os.h>
#include <core/md_allocator.h>

#include <algorithm>

namespace viamd {

struct event_system_t {
	md_array(Event) event_queue = 0;
	md_array(EventHandler*) event_handlers = 0;
};

static event_system_t event_system = {};

void event_system_register_handler(EventHandler& handler) {
	md_array_push(event_system.event_handlers, &handler, md_get_heap_allocator());
}

void event_system_enqueue_event(EventType type, EventPayloadType payload_type, const void* payload, uint64_t delay_in_ms) {
	md_timestamp_t time_now = md_time_now();

	// md_time_now() counts platform ticks - nanoseconds on unix, QPC ticks on windows - so a delay
	// in milliseconds has to be converted rather than added straight on. Adding it raw made a
	// requested second come out as a microsecond on one platform and a tenth of that on the other.
	const double ticks_per_ms = 1.0 / md_time_as_milliseconds(1);
	const md_timestamp_t delay_in_ticks = (md_timestamp_t)((double)delay_in_ms * ticks_per_ms);

	Event e = {
		.type = type,
		.payload_type = payload_type,
		.timestamp = (uint64_t)(time_now + delay_in_ticks),
		.payload = payload,
	};
	md_array_push(event_system.event_queue, e, md_get_heap_allocator());
}

void event_system_broadcast_event(EventType type, EventPayloadType payload_type, const void* payload) {
	Event e = {
		.type = type,
		.payload_type = payload_type,
		.timestamp = (uint64_t)md_time_now(),
		.payload = payload
	};

	for (size_t i = 0; i < md_array_size(event_system.event_handlers); ++i) {
		event_system.event_handlers[i]->process_events(&e, 1);
	}
}

void event_system_process_event_queue() {
	size_t num_events = md_array_size(event_system.event_queue);
	if (num_events == 0) {
		return;
	}

	// Stable, so events enqueued within the same tick keep the order they were enqueued in. Several
	// of the load sequences enqueue a run of events back to back and depend on that order.
	Event* beg = event_system.event_queue;
	std::stable_sort(beg, beg + num_events, [](const Event& a, const Event& b) {
		return a.timestamp < b.timestamp;
	});
	
	uint64_t time_now = md_time_now();
	size_t num_events_to_process = num_events;
	for (size_t i = 0; i < num_events; ++i) {
		if (event_system.event_queue[i].timestamp > time_now) {
			num_events_to_process = i;
			break;
		}
	}

	if (num_events_to_process) {
		for (size_t i = 0; i < md_array_size(event_system.event_handlers); ++i) {
			EventHandler* handler = event_system.event_handlers[i];
			handler->process_events(event_system.event_queue, num_events_to_process);
		}
		size_t num_events_left = num_events - num_events_to_process;
		if (num_events_left) {
			MEMCPY(beg, beg + num_events_to_process, sizeof(Event) * num_events_left);
		}
		md_array_shrink(event_system.event_queue, num_events_left);
	}
}

}
