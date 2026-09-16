#include "utest.h"

#include <event.h>

#include <core/md_os.h>

#include <vector>

// The event system is global and has no way to unregister a handler, so a handler registered by one
// test outlives it. Every test here therefore uses its own event type and its handler ignores
// anything else, which is also how the real components behave.

namespace {

struct recorder_t : viamd::EventHandler {
    viamd::EventType watch = 0;
    std::vector<viamd::EventType> seen;
    std::vector<const void*> payloads;
    int batches = 0;

    void process_events(const viamd::Event* events, size_t num_events) final {
        bool any = false;
        for (size_t i = 0; i < num_events; ++i) {
            if ((events[i].type & 0xFFFF0000u) != (watch & 0xFFFF0000u)) continue;
            seen.push_back(events[i].type);
            payloads.push_back(events[i].payload);
            any = true;
        }
        if (any) batches += 1;
    }

    void reset(viamd::EventType type) {
        watch = type;
        seen.clear();
        payloads.clear();
        batches = 0;
    }
};

// The recorders have to outlive every test. Registration is permanent - there is no way to take a
// handler back out - so a handler on a test's stack would be dispatched to long after that stack
// frame is gone, and the crash lands in whichever later test happens to broadcast first.
recorder_t& primary() {
    static recorder_t r;
    return r;
}

recorder_t& secondary() {
    static recorder_t r;
    return r;
}

void begin(viamd::EventType type) {
    static bool registered = false;
    if (!registered) {
        viamd::event_system_register_handler(primary());
        viamd::event_system_register_handler(secondary());
        registered = true;
    }
    primary().reset(type);
    secondary().reset(type);
    // Leave nothing queued behind from an earlier test
    viamd::event_system_process_event_queue();
    primary().reset(type);
    secondary().reset(type);
}

}  // namespace

UTEST(viamd_event, broadcast_reaches_every_handler) {
    const viamd::EventType type = 0x00010001u;

    begin(type);
    recorder_t& a = primary();
    recorder_t& b = secondary();

    int payload = 7;
    viamd::event_system_broadcast_event(type, 0, &payload);

    ASSERT_EQ(1u, (uint32_t)a.seen.size());
    ASSERT_EQ(1u, (uint32_t)b.seen.size());
    EXPECT_EQ(type, a.seen[0]);
    EXPECT_EQ((const void*)&payload, a.payloads[0]);
    EXPECT_EQ((const void*)&payload, b.payloads[0]);
}

UTEST(viamd_event, queued_events_are_delivered_once) {
    const viamd::EventType base = 0x00020000u;

    begin(base);
    recorder_t& r = primary();

    viamd::event_system_enqueue_event(base | 1);
    viamd::event_system_enqueue_event(base | 2);
    viamd::event_system_enqueue_event(base | 3);

    viamd::event_system_process_event_queue();
    ASSERT_EQ(3u, (uint32_t)r.seen.size());
    EXPECT_EQ(1, r.batches);   // one batch, not three calls

    // Draining is complete: processing again delivers nothing
    viamd::event_system_process_event_queue();
    EXPECT_EQ(3u, (uint32_t)r.seen.size());
}

UTEST(viamd_event, queued_events_keep_their_order) {
    // Load sequences enqueue a run of events back to back and depend on the order, which means the
    // queue has to be sorted stably - events enqueued within the same tick share a timestamp.
    const viamd::EventType base = 0x00030000u;

    begin(base);
    recorder_t& r = primary();

    for (uint32_t i = 0; i < 32; ++i) {
        viamd::event_system_enqueue_event(base | i);
    }
    viamd::event_system_process_event_queue();

    ASSERT_EQ(32u, (uint32_t)r.seen.size());
    for (uint32_t i = 0; i < 32; ++i) {
        EXPECT_EQ(base | i, r.seen[i]);
    }
}

UTEST(viamd_event, a_delayed_event_actually_waits) {
    // The delay is given in milliseconds while md_time_now() counts platform ticks, so it has to be
    // converted. Added raw, a requested second came out as a microsecond and the event fired at
    // once - on windows it was a different wrong answer again, since the tick rate differs.
    const viamd::EventType base = 0x00040000u;

    begin(base);
    recorder_t& r = primary();

    viamd::event_system_enqueue_event(base | 1, 0, 0, 10000);   // ten seconds out
    viamd::event_system_process_event_queue();
    EXPECT_EQ(0u, (uint32_t)r.seen.size());

    md_thread_sleep(20);
    viamd::event_system_process_event_queue();
    EXPECT_EQ(0u, (uint32_t)r.seen.size());

    // A short delay does come through once it is due, which shows the units are not merely large
    viamd::event_system_enqueue_event(base | 2, 0, 0, 5);
    md_thread_sleep(30);
    viamd::event_system_process_event_queue();
    ASSERT_EQ(1u, (uint32_t)r.seen.size());
    EXPECT_EQ(base | 2, r.seen[0]);
}

UTEST(viamd_event, a_due_event_does_not_hold_back_a_pending_one) {
    // The queue is processed up to the first event that is not due yet, so a long delayed event
    // must not stall the ones behind it.
    const viamd::EventType base = 0x00050000u;

    begin(base);
    recorder_t& r = primary();

    viamd::event_system_enqueue_event(base | 1, 0, 0, 60000);
    viamd::event_system_enqueue_event(base | 2);
    viamd::event_system_process_event_queue();

    ASSERT_EQ(1u, (uint32_t)r.seen.size());
    EXPECT_EQ(base | 2, r.seen[0]);
}
