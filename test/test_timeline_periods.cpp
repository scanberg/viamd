#include "utest.h"

#include <timeline_periods.h>
#include <serialization_utils.h>

#include <core/md_allocator.h>
#include <core/md_writer.h>

#include <string.h>

UTEST(timeline_periods, frames_fixed_and_following) {
    int beg = -1, end = -1;
    TimelinePeriod r = {};

    r.beg = 10; r.end = 20;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, 1001));
    EXPECT_EQ(10, beg); EXPECT_EQ(20, end);

    // An end dragged past the other: the lower one is the beginning
    r.beg = 20; r.end = 10;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, 1001));
    EXPECT_EQ(10, beg); EXPECT_EQ(20, end);

    // Beyond the end of the trajectory: clamped, or nothing when it is all beyond
    r.beg = 990; r.end = 1200;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, 1001));
    EXPECT_EQ(990, beg); EXPECT_EQ(1000, end);
    r.beg = 2000; r.end = 2100;
    EXPECT_FALSE(timeline_period_frames(&beg, &end, r, 0.0, 1001));

    // Following: the last 60 frames, from the current frame rounded
    r.follow = true;
    r.beg = -60; r.end = 0;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 100.4, 1001));
    EXPECT_EQ(40, beg); EXPECT_EQ(100, end);
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 30.0, 1001));
    EXPECT_EQ(0, beg); EXPECT_EQ(30, end);

    // Looking ahead from the last frame covers nothing
    r.beg = 10; r.end = 60;
    EXPECT_FALSE(timeline_period_frames(&beg, &end, r, 1000.0, 1001));

    EXPECT_FALSE(timeline_period_frames(&beg, &end, r, 0.0, 0));
}

UTEST(timeline_periods, follow_keeps_the_frames_it_covers) {
    TimelinePeriod r = {};
    r.beg = 200; r.end = 100;

    timeline_period_set_follow(&r, true, 150.0, 1001);
    EXPECT_TRUE(r.follow);
    EXPECT_EQ(-50, r.beg);
    EXPECT_EQ(50, r.end);

    timeline_period_set_follow(&r, false, 980.0, 1001);
    EXPECT_FALSE(r.follow);
    EXPECT_EQ(930, r.beg);
    EXPECT_EQ(1000, r.end);     // 1030 is past the last frame
}

UTEST(timeline_periods, add_names_colours_and_remove) {
    TimelinePeriods rs = {};

    TimelinePeriod* a = timeline_periods_add(&rs, nullptr, 0, 100, false);
    TimelinePeriod* b = timeline_periods_add(&rs, nullptr, 300, 200, false);
    TimelinePeriod* m = timeline_periods_add(&rs, nullptr, 50, 50, false);
    ASSERT_TRUE(a && b && m);
    EXPECT_STREQ("Period 1", rs.items[0].name);
    EXPECT_STREQ("Period 2", rs.items[1].name);
    EXPECT_STREQ("Marker 1", rs.items[2].name);
    EXPECT_EQ(200, rs.items[1].beg);
    EXPECT_EQ(300, rs.items[1].end);
    EXPECT_TRUE(timeline_period_is_marker(rs.items[2]));

    // Each its own colour while the palette lasts
    EXPECT_FALSE(rs.items[0].color.x == rs.items[1].color.x && rs.items[0].color.y == rs.items[1].color.y && rs.items[0].color.z == rs.items[1].color.z);

    const uint32_t id_b = rs.items[1].id;
    rs.legacy_filter_id = id_b;
    timeline_periods_remove(&rs, id_b);
    EXPECT_EQ(2, rs.count);
    EXPECT_EQ(0u, rs.legacy_filter_id);
    EXPECT_EQ(-1, timeline_periods_index(&rs, id_b));
    EXPECT_STREQ("Marker 1", rs.items[1].name);

    // The name freed is the next one given
    TimelinePeriod* c = timeline_periods_add(&rs, nullptr, 10, 20, false);
    ASSERT_TRUE(c != nullptr);
    EXPECT_STREQ("Period 2", c->name);

    for (int i = rs.count; i < TIMELINE_MAX_PERIODS; ++i) {
        ASSERT_TRUE(timeline_periods_add(&rs, "x", 0, 1, false) != nullptr);
    }
    EXPECT_TRUE(timeline_periods_add(&rs, "one too many", 0, 1, false) == nullptr);
}

UTEST(timeline_periods, draw_order_puts_short_periods_on_top) {
    TimelinePeriods rs = {};
    timeline_periods_add(&rs, "Production",  200, 1000, false);  // 0
    timeline_periods_add(&rs, "Event",       430, 470,  false);  // 1
    timeline_periods_add(&rs, "Equilibrate", 0,   200,  false);  // 2
    timeline_periods_add(&rs, "Same length", 600, 800,  false);  // 3, as long as Equilibrate but younger
    timeline_periods_add(&rs, "Marker",      450, 450,  false);  // 4
    timeline_periods_add(&rs, "Beyond",      5000, 6000, false); // 5, covers nothing

    int order[TIMELINE_MAX_PERIODS];
    const int n = timeline_periods_draw_order(order, &rs, 0.0, 1001);
    ASSERT_EQ(5, n);
    EXPECT_EQ(0, order[0]);
    EXPECT_EQ(2, order[1]);
    EXPECT_EQ(3, order[2]);
    EXPECT_EQ(1, order[3]);
    EXPECT_EQ(4, order[4]);

    int depth[TIMELINE_MAX_PERIODS];
    timeline_periods_nesting(depth, order, n, &rs, 0.0, 1001);
    EXPECT_EQ(0, depth[0]);     // Production
    EXPECT_EQ(0, depth[1]);     // Equilibrate only touches Production
    EXPECT_EQ(1, depth[2]);     // Same length lies on Production
    EXPECT_EQ(1, depth[3]);     // Event lies on Production
    EXPECT_EQ(2, depth[4]);     // The marker on Production and Event
}

UTEST(timeline_periods, hit_test) {
    TimelinePeriodBar bars[4];
    for (auto& b : bars) { b.y0 = 0; b.y1 = 20; b.visible = true; }
    bars[0].x0 = 0;   bars[0].x1 = 100;     // A
    bars[1].x0 = 100; bars[1].x1 = 300;     // B, beginning where A ends
    bars[2].x0 = 150; bars[2].x1 = 190;     // C, on B
    bars[3].x0 = bars[3].x1 = 260; bars[3].marker = true;    // A marker, on B

    const int top_first[4] = { 3, 2, 0, 1 };
    const float reach = 7.0f, marker_reach = 6.0f;

    // At the boundary between A and B, the side decides
    TimelinePeriodHit h = timeline_periods_hit_test(bars, top_first, 4, 98.0f, 10.0f, reach, marker_reach);
    EXPECT_EQ(0, h.index); EXPECT_EQ(TimelinePeriodPart_End, h.part);
    h = timeline_periods_hit_test(bars, top_first, 4, 102.0f, 10.0f, reach, marker_reach);
    EXPECT_EQ(1, h.index); EXPECT_EQ(TimelinePeriodPart_Beg, h.part);

    // The period on top
    h = timeline_periods_hit_test(bars, top_first, 4, 170.0f, 10.0f, reach, marker_reach);
    EXPECT_EQ(2, h.index); EXPECT_EQ(TimelinePeriodPart_Body, h.part);
    h = timeline_periods_hit_test(bars, top_first, 4, 152.0f, 10.0f, reach, marker_reach);
    EXPECT_EQ(2, h.index); EXPECT_EQ(TimelinePeriodPart_Beg, h.part);
    h = timeline_periods_hit_test(bars, top_first, 4, 220.0f, 10.0f, reach, marker_reach);
    EXPECT_EQ(1, h.index); EXPECT_EQ(TimelinePeriodPart_Body, h.part);

    // Markers before anything else in reach
    h = timeline_periods_hit_test(bars, top_first, 4, 263.0f, 10.0f, reach, marker_reach);
    EXPECT_EQ(3, h.index); EXPECT_EQ(TimelinePeriodPart_Body, h.part);

    // Outside the lane
    h = timeline_periods_hit_test(bars, top_first, 4, 170.0f, 40.0f, reach, marker_reach);
    EXPECT_EQ(-1, h.index);

    // The middle of a narrow bar moves it, its ends resize it
    TimelinePeriodBar narrow[1];
    narrow[0].x0 = 100; narrow[0].x1 = 109; narrow[0].y1 = 20; narrow[0].visible = true;
    const int one[1] = { 0 };
    h = timeline_periods_hit_test(narrow, one, 1, 104.5f, 10.0f, reach, marker_reach);
    EXPECT_EQ(TimelinePeriodPart_Body, h.part);
    h = timeline_periods_hit_test(narrow, one, 1, 101.0f, 10.0f, reach, marker_reach);
    EXPECT_EQ(TimelinePeriodPart_Beg, h.part);
}

UTEST(timeline_periods, workspace_round_trip) {
    TimelinePeriods rs = {};
    TimelinePeriod* a = timeline_periods_add(&rs, "Unfolding", 430, 470, false);
    strcpy(a->note, "The helix unwinds here.\nCompare with the second run.");
    TimelinePeriod* b = timeline_periods_add(&rs, "Recent", -60, 0, true);
    b->color = vec4_set(0.25f, 0.5f, 0.75f, 1.0f);
    timeline_periods_add(&rs, nullptr, 12, 12, false);

    viamd::serialization_state_t w = {};
	md_str_writer_init(&w.sw, md_get_heap_allocator());
    timeline_periods_serialize(w, &rs);

    viamd::deserialization_state_t r = {};
    r.text = md_str_writer_to_str(&w.sw);

    TimelinePeriods out = {};
    str_t section;
    int sections = 0;
    while (viamd::next_section_header(section, r)) {
        ASSERT_TRUE(str_eq(section, STR_LIT("TimelinePeriod")));
        ASSERT_TRUE(timeline_periods_deserialize_period(r, &out));
        sections += 1;
    }
    md_str_writer_free(&w.sw);

    ASSERT_EQ(3, sections);
    ASSERT_EQ(3, out.count);
    EXPECT_STREQ("Unfolding", out.items[0].name);
    EXPECT_STREQ(a->note, out.items[0].note);
    EXPECT_EQ(430, out.items[0].beg);
    EXPECT_EQ(470, out.items[0].end);
    EXPECT_FALSE(out.items[0].follow);

    EXPECT_STREQ("Recent", out.items[1].name);
    EXPECT_TRUE(out.items[1].follow);
    EXPECT_EQ(-60, out.items[1].beg);
    EXPECT_EQ(0, out.items[1].end);
    EXPECT_NEAR(0.25f, out.items[1].color.x, 1e-6f);
    EXPECT_NEAR(0.75f, out.items[1].color.z, 1e-6f);
    EXPECT_EQ(0u, out.legacy_filter_id);

    EXPECT_STREQ("Marker 1", out.items[2].name);
    EXPECT_TRUE(timeline_period_is_marker(out.items[2]));

    // Each keeps its id, which what refers to it was saved with
    for (int i = 0; i < 3; ++i) {
        EXPECT_EQ(rs.items[i].id, out.items[i].id);
    }
}

UTEST(timeline_periods, ids_survive_a_save_and_a_load) {
    // Ids with a gap, as after a period was removed
    TimelinePeriods rs = {};
    timeline_periods_add(&rs, "A", 0, 10, false);
    const uint32_t gone = timeline_periods_add(&rs, "B", 10, 20, false)->id;
    const uint32_t c = timeline_periods_add(&rs, "C", 20, 30, false)->id;
    timeline_periods_remove(&rs, gone);

    viamd::serialization_state_t w = {};
    md_str_writer_init(&w.sw, md_get_heap_allocator());
    timeline_periods_serialize(w, &rs);

    // Read into periods which already have one, of C's id: C gets another, and new ones never take a saved one
    TimelinePeriods out = {};
    TimelinePeriod* holder = timeline_periods_add(&out, "Already there", 0, 1, false);
    holder->id = c;
    out.next_id = c + 1;
    const int before = out.count;

    viamd::deserialization_state_t r = {};
    r.text = md_str_writer_to_str(&w.sw);
    str_t section;
    while (viamd::next_section_header(section, r)) {
        ASSERT_TRUE(timeline_periods_deserialize_period(r, &out));
    }
    md_str_writer_free(&w.sw);

    ASSERT_EQ(before + 2, out.count);
    const TimelinePeriod& a = out.items[before];
    const TimelinePeriod& cc = out.items[before + 1];
    EXPECT_EQ(rs.items[0].id, a.id);            // free: kept
    EXPECT_NE(c, cc.id);                        // taken: another
    EXPECT_STREQ("C", cc.name);
    for (int i = 0; i < out.count; ++i) {
        for (int j = i + 1; j < out.count; ++j) {
            EXPECT_NE(out.items[i].id, out.items[j].id);
        }
    }
    TimelinePeriod* next = timeline_periods_add(&out, "Next", 0, 1, false);
    ASSERT_TRUE(next != nullptr);
    for (int i = 0; i < out.count - 1; ++i) {
        EXPECT_NE(out.items[i].id, next->id);
    }
}
