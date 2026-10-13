#include "utest.h"

#include <timeline_ranges.h>
#include <serialization_utils.h>

#include <core/md_allocator.h>
#include <core/md_str_builder.h>

#include <string.h>

UTEST(timeline_ranges, frames_fixed_and_following) {
    int beg = -1, end = -1;
    TimelineRange r = {};

    r.beg = 10; r.end = 20;
    ASSERT_TRUE(timeline_range_frames(&beg, &end, r, 0.0, 1001));
    EXPECT_EQ(10, beg); EXPECT_EQ(20, end);

    // An end dragged past the other: the lower one is the beginning
    r.beg = 20; r.end = 10;
    ASSERT_TRUE(timeline_range_frames(&beg, &end, r, 0.0, 1001));
    EXPECT_EQ(10, beg); EXPECT_EQ(20, end);

    // Beyond the end of the trajectory: clamped, or nothing when it is all beyond
    r.beg = 990; r.end = 1200;
    ASSERT_TRUE(timeline_range_frames(&beg, &end, r, 0.0, 1001));
    EXPECT_EQ(990, beg); EXPECT_EQ(1000, end);
    r.beg = 2000; r.end = 2100;
    EXPECT_FALSE(timeline_range_frames(&beg, &end, r, 0.0, 1001));

    // Following: the last 60 frames, from the current frame rounded
    r.follow = true;
    r.beg = -60; r.end = 0;
    ASSERT_TRUE(timeline_range_frames(&beg, &end, r, 100.4, 1001));
    EXPECT_EQ(40, beg); EXPECT_EQ(100, end);
    ASSERT_TRUE(timeline_range_frames(&beg, &end, r, 30.0, 1001));
    EXPECT_EQ(0, beg); EXPECT_EQ(30, end);

    // Looking ahead from the last frame covers nothing
    r.beg = 10; r.end = 60;
    EXPECT_FALSE(timeline_range_frames(&beg, &end, r, 1000.0, 1001));

    EXPECT_FALSE(timeline_range_frames(&beg, &end, r, 0.0, 0));
}

UTEST(timeline_ranges, follow_keeps_the_frames_it_covers) {
    TimelineRange r = {};
    r.beg = 200; r.end = 100;

    timeline_range_set_follow(&r, true, 150.0, 1001);
    EXPECT_TRUE(r.follow);
    EXPECT_EQ(-50, r.beg);
    EXPECT_EQ(50, r.end);

    timeline_range_set_follow(&r, false, 980.0, 1001);
    EXPECT_FALSE(r.follow);
    EXPECT_EQ(930, r.beg);
    EXPECT_EQ(1000, r.end);     // 1030 is past the last frame
}

UTEST(timeline_ranges, add_names_colours_and_remove) {
    TimelineRanges rs = {};

    TimelineRange* a = timeline_ranges_add(&rs, nullptr, 0, 100, false);
    TimelineRange* b = timeline_ranges_add(&rs, nullptr, 300, 200, false);
    TimelineRange* m = timeline_ranges_add(&rs, nullptr, 50, 50, false);
    ASSERT_TRUE(a && b && m);
    EXPECT_STREQ("Range 1", rs.items[0].name);
    EXPECT_STREQ("Range 2", rs.items[1].name);
    EXPECT_STREQ("Marker 1", rs.items[2].name);
    EXPECT_EQ(200, rs.items[1].beg);
    EXPECT_EQ(300, rs.items[1].end);
    EXPECT_TRUE(timeline_range_is_marker(rs.items[2]));

    // Each its own colour while the palette lasts
    EXPECT_FALSE(rs.items[0].color.x == rs.items[1].color.x && rs.items[0].color.y == rs.items[1].color.y && rs.items[0].color.z == rs.items[1].color.z);

    const uint32_t id_b = rs.items[1].id;
    rs.filter_id = id_b;
    timeline_ranges_remove(&rs, id_b);
    EXPECT_EQ(2, rs.count);
    EXPECT_EQ(0u, rs.filter_id);
    EXPECT_EQ(-1, timeline_ranges_index(&rs, id_b));
    EXPECT_STREQ("Marker 1", rs.items[1].name);

    // The name freed is the next one given
    TimelineRange* c = timeline_ranges_add(&rs, nullptr, 10, 20, false);
    ASSERT_TRUE(c != nullptr);
    EXPECT_STREQ("Range 2", c->name);

    for (int i = rs.count; i < TIMELINE_MAX_RANGES; ++i) {
        ASSERT_TRUE(timeline_ranges_add(&rs, "x", 0, 1, false) != nullptr);
    }
    EXPECT_TRUE(timeline_ranges_add(&rs, "one too many", 0, 1, false) == nullptr);
}

UTEST(timeline_ranges, filter_follows_the_marked_range) {
    TimelineRanges rs = {};
    int beg, end;
    EXPECT_FALSE(timeline_ranges_filter_frames(&beg, &end, &rs, 0.0, 100));

    TimelineRange* r = timeline_ranges_add(&rs, "Window", -5, 5, true);
    rs.filter_id = r->id;
    ASSERT_TRUE(timeline_ranges_filter_frames(&beg, &end, &rs, 50.0, 100));
    EXPECT_EQ(45, beg); EXPECT_EQ(55, end);
}

UTEST(timeline_ranges, draw_order_puts_short_ranges_on_top) {
    TimelineRanges rs = {};
    timeline_ranges_add(&rs, "Production",  200, 1000, false);  // 0
    timeline_ranges_add(&rs, "Event",       430, 470,  false);  // 1
    timeline_ranges_add(&rs, "Equilibrate", 0,   200,  false);  // 2
    timeline_ranges_add(&rs, "Same length", 600, 800,  false);  // 3, as long as Equilibrate but younger
    timeline_ranges_add(&rs, "Marker",      450, 450,  false);  // 4
    timeline_ranges_add(&rs, "Beyond",      5000, 6000, false); // 5, covers nothing

    int order[TIMELINE_MAX_RANGES];
    const int n = timeline_ranges_draw_order(order, &rs, 0.0, 1001);
    ASSERT_EQ(5, n);
    EXPECT_EQ(0, order[0]);
    EXPECT_EQ(2, order[1]);
    EXPECT_EQ(3, order[2]);
    EXPECT_EQ(1, order[3]);
    EXPECT_EQ(4, order[4]);

    int depth[TIMELINE_MAX_RANGES];
    timeline_ranges_nesting(depth, order, n, &rs, 0.0, 1001);
    EXPECT_EQ(0, depth[0]);     // Production
    EXPECT_EQ(0, depth[1]);     // Equilibrate only touches Production
    EXPECT_EQ(1, depth[2]);     // Same length lies on Production
    EXPECT_EQ(1, depth[3]);     // Event lies on Production
    EXPECT_EQ(2, depth[4]);     // The marker on Production and Event
}

UTEST(timeline_ranges, hit_test) {
    TimelineRangeBar bars[4];
    for (auto& b : bars) { b.y0 = 0; b.y1 = 20; b.visible = true; }
    bars[0].x0 = 0;   bars[0].x1 = 100;     // A
    bars[1].x0 = 100; bars[1].x1 = 300;     // B, beginning where A ends
    bars[2].x0 = 150; bars[2].x1 = 190;     // C, on B
    bars[3].x0 = bars[3].x1 = 260; bars[3].marker = true;    // A marker, on B

    const int top_first[4] = { 3, 2, 0, 1 };
    const float reach = 7.0f, marker_reach = 6.0f;

    // At the boundary between A and B, the side decides
    TimelineRangeHit h = timeline_ranges_hit_test(bars, top_first, 4, 98.0f, 10.0f, reach, marker_reach);
    EXPECT_EQ(0, h.index); EXPECT_EQ(TimelineRangePart_End, h.part);
    h = timeline_ranges_hit_test(bars, top_first, 4, 102.0f, 10.0f, reach, marker_reach);
    EXPECT_EQ(1, h.index); EXPECT_EQ(TimelineRangePart_Beg, h.part);

    // The range on top
    h = timeline_ranges_hit_test(bars, top_first, 4, 170.0f, 10.0f, reach, marker_reach);
    EXPECT_EQ(2, h.index); EXPECT_EQ(TimelineRangePart_Body, h.part);
    h = timeline_ranges_hit_test(bars, top_first, 4, 152.0f, 10.0f, reach, marker_reach);
    EXPECT_EQ(2, h.index); EXPECT_EQ(TimelineRangePart_Beg, h.part);
    h = timeline_ranges_hit_test(bars, top_first, 4, 220.0f, 10.0f, reach, marker_reach);
    EXPECT_EQ(1, h.index); EXPECT_EQ(TimelineRangePart_Body, h.part);

    // Markers before anything else in reach
    h = timeline_ranges_hit_test(bars, top_first, 4, 263.0f, 10.0f, reach, marker_reach);
    EXPECT_EQ(3, h.index); EXPECT_EQ(TimelineRangePart_Body, h.part);

    // Outside the lane
    h = timeline_ranges_hit_test(bars, top_first, 4, 170.0f, 40.0f, reach, marker_reach);
    EXPECT_EQ(-1, h.index);

    // The middle of a narrow bar moves it, its ends resize it
    TimelineRangeBar narrow[1];
    narrow[0].x0 = 100; narrow[0].x1 = 109; narrow[0].y1 = 20; narrow[0].visible = true;
    const int one[1] = { 0 };
    h = timeline_ranges_hit_test(narrow, one, 1, 104.5f, 10.0f, reach, marker_reach);
    EXPECT_EQ(TimelineRangePart_Body, h.part);
    h = timeline_ranges_hit_test(narrow, one, 1, 101.0f, 10.0f, reach, marker_reach);
    EXPECT_EQ(TimelineRangePart_Beg, h.part);
}

UTEST(timeline_ranges, workspace_round_trip) {
    TimelineRanges rs = {};
    TimelineRange* a = timeline_ranges_add(&rs, "Unfolding", 430, 470, false);
    strcpy(a->note, "The helix unwinds here.\nCompare with the second run.");
    TimelineRange* b = timeline_ranges_add(&rs, "Recent", -60, 0, true);
    b->color = vec4_set(0.25f, 0.5f, 0.75f, 1.0f);
    rs.filter_id = b->id;
    timeline_ranges_add(&rs, nullptr, 12, 12, false);

    viamd::serialization_state_t w = {};
    md_strb_init(&w.sb, md_get_heap_allocator());
    timeline_ranges_serialize(w, &rs);

    viamd::deserialization_state_t r = {};
    r.text = md_strb_to_str(w.sb);

    TimelineRanges out = {};
    str_t section;
    int sections = 0;
    while (viamd::next_section_header(section, r)) {
        ASSERT_TRUE(str_eq(section, STR_LIT("TimelineRange")));
        ASSERT_TRUE(timeline_ranges_deserialize_range(r, &out));
        sections += 1;
    }
    md_strb_free(&w.sb);

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
    EXPECT_EQ(out.items[1].id, out.filter_id);

    EXPECT_STREQ("Marker 1", out.items[2].name);
    EXPECT_TRUE(timeline_range_is_marker(out.items[2]));
}
