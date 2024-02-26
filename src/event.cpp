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
	md_array_push(event_system.event_handlers, &handler, md_heap_allocator);
}

void event_system_enqueue_event(EventType type, const void* payload, uint64_t delay_in_ms) {
	md_timestamp_t time_now = md_time_current();
	Event e = {
		.type = type,
		.timestamp = time_now + delay_in_ms,
		.payload = payload,
	};
	md_array_push(event_system.event_queue, e, md_heap_allocator);
}

void event_system_broadcast_event(EventType type, const void* payload) {
	Event e = {
		.type = type,
		.timestamp = (uint64_t)md_time_current(),
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

	Event* beg = event_system.event_queue;
	std::sort(beg, beg + num_events, [](const Event& a, const Event& b) {
		return a.timestamp > b.timestamp;
	});
	
	uint64_t time_now = md_time_current();
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
