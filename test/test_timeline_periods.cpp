#include "utest.h"

#include <timeline_periods.h>
#include <serialization_utils.h>

#include <core/md_allocator.h>
#include <core/md_writer.h>
#include <core/md_vec_math.h>

#include <string.h>

// Frames every 10 ps from 0 to 10000 ps
static double g_times[1001];
static ClockAxis axis_10ps(double offset = 0.0, double scale = 1.0) {
    for (int i = 0; i < 1001; ++i) g_times[i] = 10.0 * i;
    ClockAxis a;
    a.times = g_times;
    a.count = 1001;
    a.offset = offset;
    a.scale = scale;
    return a;
}

UTEST(timeline_periods, frames_fixed_and_following) {
    const ClockAxis a = axis_10ps();
    int beg = -1, end = -1;
    TimelinePeriod r = {};

    r.beg = 100; r.end = 200;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, a));
    EXPECT_EQ(10, beg); EXPECT_EQ(20, end);

    // An end dragged past the other: the lower one is the beginning
    r.beg = 200; r.end = 100;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, a));
    EXPECT_EQ(10, beg); EXPECT_EQ(20, end);

    // Ends between frames: the frames within
    r.beg = 95; r.end = 205;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, a));
    EXPECT_EQ(10, beg); EXPECT_EQ(20, end);
    r.beg = 101; r.end = 199;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, a));
    EXPECT_EQ(11, beg); EXPECT_EQ(19, end);

    // Beyond the last frame: clamped, or nothing when it is all beyond
    r.beg = 9900; r.end = 12000;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, a));
    EXPECT_EQ(990, beg); EXPECT_EQ(1000, end);
    r.beg = 20000; r.end = 21000;
    EXPECT_FALSE(timeline_period_frames(&beg, &end, r, 0.0, a));
    r.beg = -500; r.end = -100;
    EXPECT_FALSE(timeline_period_frames(&beg, &end, r, 0.0, a));

    // Following: the last 600 ps, from the current time
    r.follow = true;
    r.beg = -600; r.end = 0;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 1000.0, a));
    EXPECT_EQ(40, beg); EXPECT_EQ(100, end);
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 300.0, a));
    EXPECT_EQ(0, beg); EXPECT_EQ(30, end);

    // Looking ahead from the last frame covers nothing
    r.beg = 100; r.end = 600;
    EXPECT_FALSE(timeline_period_frames(&beg, &end, r, 10000.0, a));

    EXPECT_FALSE(timeline_period_frames(&beg, &end, r, 0.0, ClockAxis{}));
}

UTEST(timeline_periods, frames_of_runs_sampled_differently) {
    // The same time of a run every 10 ps and one every 40 ps, the latter 1000 ps later on the clock and in ns
    const ClockAxis fine = axis_10ps();
    double coarse_times[101];
    for (int i = 0; i < 101; ++i) coarse_times[i] = 0.04 * i;     // ns
    ClockAxis coarse;
    coarse.times = coarse_times;
    coarse.count = 101;
    coarse.offset = 1.0;        // ns
    coarse.scale = 1000.0;      // into ps

    TimelinePeriod r = {};
    r.beg = 2000; r.end = 3000;
    int beg, end;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, fine));
    EXPECT_EQ(200, beg); EXPECT_EQ(300, end);
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, coarse));
    EXPECT_EQ(25, beg); EXPECT_EQ(50, end);

    // A marker: the frame at its time, and of a run without a frame there the nearest
    r.beg = r.end = 2000;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, fine));
    EXPECT_EQ(200, beg); EXPECT_EQ(200, end);
    r.beg = r.end = 2015;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, coarse));
    EXPECT_EQ(25, beg); EXPECT_EQ(25, end);
    r.beg = r.end = 2025;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, coarse));
    EXPECT_EQ(26, beg); EXPECT_EQ(26, end);

    // Narrower than a frame of the coarse run, between two of its frames: the one nearest its middle
    r.beg = 2042; r.end = 2058;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, coarse));
    EXPECT_EQ(26, beg); EXPECT_EQ(26, end);

    // Before the coarse run begins
    r.beg = 0; r.end = 500;
    EXPECT_FALSE(timeline_period_frames(&beg, &end, r, 0.0, coarse));
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, fine));
    EXPECT_EQ(0, beg); EXPECT_EQ(50, end);
}

UTEST(timeline_periods, frames_of_ends_at_frames_computed_alike) {
    // Ends made from frames, as a workspace from before the clock has them, cover those frames exactly, also where the
    // times are not exact in binary
    double times[50];
    for (int i = 0; i < 50; ++i) times[i] = 0.1 * i;
    ClockAxis a;
    a.times = times;
    a.count = 50;
    a.offset = 0.3;
    a.scale = 3.0;
    TimelinePeriod r = {};
    r.beg = 7;
    r.end = 21;
    timeline_period_frames_to_time(&r, a);
    int beg, end;
    ASSERT_TRUE(timeline_period_frames(&beg, &end, r, 0.0, a));
    EXPECT_EQ(7, beg);
    EXPECT_EQ(21, end);

    // Following: offsets in frames become offsets of as many mean frame times
    TimelinePeriod f = {};
    f.follow = true;
    f.beg = -10;
    f.end = 0;
    timeline_period_frames_to_time(&f, a);
    EXPECT_NEAR(-3.0, f.beg, 1e-9);
    EXPECT_NEAR(0.0, f.end, 1e-9);
}

UTEST(timeline_periods, follow_centers_and_keeps_the_extent) {
    TimelinePeriod r = {};
    r.beg = 2000; r.end = 1000;

    // Centered on the current time, wherever that is
    timeline_period_set_follow(&r, true, 4000.0, 10.0, 0.0, 10000.0);
    EXPECT_TRUE(r.follow);
    EXPECT_EQ(-500.0, r.beg);
    EXPECT_EQ(500.0, r.end);

    // Fixed again, it keeps the time it covers, within the clock's span
    timeline_period_set_follow(&r, false, 9800.0, 10.0, 0.0, 10000.0);
    EXPECT_FALSE(r.follow);
    EXPECT_EQ(9300.0, r.beg);
    EXPECT_EQ(10000.0, r.end);     // 10300 is past the end

    // An extent of an odd number of steps leans into the past
    r.beg = 100; r.end = 190;
    timeline_period_set_follow(&r, true, 0.0, 10.0, 0.0, 10000.0);
    EXPECT_EQ(-50.0, r.beg);
    EXPECT_EQ(40.0, r.end);

    // A marker follows on the current time
    r.follow = false;
    r.beg = 700; r.end = 700;
    timeline_period_set_follow(&r, true, 3000.0, 10.0, 0.0, 10000.0);
    EXPECT_EQ(0.0, r.beg);
    EXPECT_EQ(0.0, r.end);
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
    const int n = timeline_periods_draw_order(order, &rs, 0.0, 0.0, 1000.0);
    ASSERT_EQ(5, n);
    EXPECT_EQ(0, order[0]);
    EXPECT_EQ(2, order[1]);
    EXPECT_EQ(3, order[2]);
    EXPECT_EQ(1, order[3]);
    EXPECT_EQ(4, order[4]);

    int depth[TIMELINE_MAX_PERIODS];
    timeline_periods_nesting(depth, order, n, &rs, 0.0, 0.0, 1000.0);
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
    TimelinePeriod* a = timeline_periods_add(&rs, "Unfolding", 430.25, 470, false);
    strcpy(a->note, "The helix unwinds here.\nCompare with the second run.");
    TimelinePeriod* b = timeline_periods_add(&rs, "Recent", -60, 0, true);
    b->color = md_float4_set(0.25f, 0.5f, 0.75f, 1.0f);
    timeline_periods_add(&rs, nullptr, 12, 12, false);

    viamd::serialization_state_t w = {};
	md_str_writer_init(&w.sw, md_get_heap_allocator());
    timeline_periods_serialize(w, &rs, "ps");

    viamd::deserialization_state_t r = {};
    r.text = md_str_writer_to_str(&w.sw);

    TimelinePeriods out = {};
    md_str_t section;
    int sections = 0;
    while (viamd::next_section_header(section, r)) {
        ASSERT_TRUE(md_str_eq(section, MD_STR_LIT("TimelinePeriod")));
        char unit[16];
        bool in_frames = true;
        ASSERT_TRUE(timeline_periods_deserialize_period(r, &out, unit, sizeof(unit), &in_frames) != nullptr);
        EXPECT_STREQ("ps", unit);
        EXPECT_FALSE(in_frames);
        sections += 1;
    }
    md_str_writer_free(&w.sw);

    ASSERT_EQ(3, sections);
    ASSERT_EQ(3, out.count);
    EXPECT_STREQ("Unfolding", out.items[0].name);
    EXPECT_STREQ(a->note, out.items[0].note);
    EXPECT_EQ(430.25, out.items[0].beg);
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
    timeline_periods_serialize(w, &rs, "ps");

    // Read into periods which already have one, of C's id: C gets another, and new ones never take a saved one
    TimelinePeriods out = {};
    TimelinePeriod* holder = timeline_periods_add(&out, "Already there", 0, 1, false);
    holder->id = c;
    out.next_id = c + 1;
    const int before = out.count;

    viamd::deserialization_state_t r = {};
    r.text = md_str_writer_to_str(&w.sw);
    md_str_t section;
    while (viamd::next_section_header(section, r)) {
        char unit[16];
        bool in_frames;
        ASSERT_TRUE(timeline_periods_deserialize_period(r, &out, unit, sizeof(unit), &in_frames) != nullptr);
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

UTEST(timeline_periods, reads_frames_from_before_the_clock) {
    const char* text =
        "[TimelinePeriod]\n"
        "Id=3\n"
        "Name=Production\n"
        "Follow=0\n"
        "Frames=200,1000\n";
    viamd::deserialization_state_t r = {};
    r.text = md_str_from_cstr(text);
    TimelinePeriods out = {};
    md_str_t section;
    ASSERT_TRUE(viamd::next_section_header(section, r));
    char unit[16];
    bool in_frames = false;
    const TimelinePeriod* p = timeline_periods_deserialize_period(r, &out, unit, sizeof(unit), &in_frames);
    ASSERT_TRUE(p != nullptr);
    EXPECT_TRUE(in_frames);
    EXPECT_STREQ("", unit);
    EXPECT_EQ(200.0, p->beg);
    EXPECT_EQ(1000.0, p->end);
    EXPECT_EQ(3u, p->id);
}
