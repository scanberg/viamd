#include "event.h"
#include <core/md_array.h>
#include <core/md_os.h>
#include <core/md_allocator.h>

#include <algorithm>

namespace viamd {

struct event_system_t {
	md_array(Event) queue = 0;
	md_array(EventHandler*) handlers = 0;
    md_allocator_i* alloc = md_get_heap_allocator();
};

static event_system_t event_system = {};

void event_system_register_handler(EventHandler& handler) {
	md_array_push(event_system.handlers, &handler, event_system.alloc);
}

void event_system_enqueue_event(EventType type, EventPayloadType payload_type, const void* payload, uint64_t delay_in_ms) {
	md_tick_t time_now = md_tick_now();

	// md_tick_now() counts platform ticks - nanoseconds on unix, QueryPerformanceCounter ticks
	// on windows - so a delay given in milliseconds has to be converted rather than added raw.
	// Added raw, a requested second came out as a microsecond on unix, and as a different wrong
	// answer on windows, where the tick rate varies by machine. mdlib only offers the
	// ticks-to-ms direction, so the rate is derived from it once.
	static const double ms_per_tick = md_tick_to_milliseconds(1);
	const md_tick_t delay_in_ticks = (md_tick_t)((double)delay_in_ms / ms_per_tick);

	Event e = {
		.type = type,
		.payload_type = payload_type,
		.process_time = (uint64_t)(time_now + delay_in_ticks),
		.payload = payload,
	};
	md_array_push(event_system.queue, e, event_system.alloc);
}

void event_system_broadcast_event(EventType type, EventPayloadType payload_type, const void* payload) {
	Event e = {
		.type = type,
		.payload_type = payload_type,
		.process_time = (uint64_t)md_tick_now(),
		.payload = payload
	};

	for (size_t i = 0; i < md_array_size(event_system.handlers); ++i) {
		event_system.handlers[i]->process_events(&e, 1);
	}
}

void event_system_process_event_queue() {
	size_t num_events = md_array_size(event_system.queue);
	if (num_events == 0) {
		return;
	}

	Event* beg = event_system.queue;
	std::sort(beg, beg + num_events, [](const Event& a, const Event& b) {
		return a.process_time < b.process_time;
	});
	
	uint64_t time_now = md_tick_now();
	size_t num_events_to_process = num_events;
	for (size_t i = 0; i < num_events; ++i) {
		if (event_system.queue[i].process_time > time_now) {
			num_events_to_process = i;
			break;
		}
	}

	if (num_events_to_process) {
		for (size_t i = 0; i < md_array_size(event_system.handlers); ++i) {
			EventHandler* handler = event_system.handlers[i];
			handler->process_events(event_system.queue, num_events_to_process);
		}
		size_t num_events_left = num_events - num_events_to_process;
		if (num_events_left) {
			MEMCPY(beg, beg + num_events_to_process, sizeof(Event) * num_events_left);
		}
		md_array_shrink(event_system.queue, num_events_left);
	}
}

}
