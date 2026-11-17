#include "utest.h"

#include "../svg_core.h"

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#include <stdio.h>
#include <string.h>

// Every kind of command the flow diagram draws, written to a file and read back
UTEST(svg, export) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));

    vg_scene_t scene = {};
    vg_scene_init(&scene, md_float2_t{ 320.0f, 180.0f }, arena);

    vg_add_rect_filled(&scene, md_float2_t{ 12.0f, 32.0f }, md_float2_t{ 88.0f, 72.0f }, VG_COLOR_RGBA(0, 128, 255, 200), 6.0f);
    vg_add_rect(&scene, md_float2_t{ 100.0f, 32.0f }, md_float2_t{ 140.0f, 72.0f }, VG_COLOR_RGB(0, 0, 0), 0.0f, 1.5f);

    vg_style_t gradient_style = vg_style_stroke(VG_COLOR_RGB(0, 0, 0), 6.0f);
    gradient_style.stroke = vg_paint_linear_gradient(md_float2_t{ 40.0f, 120.0f }, md_float2_t{ 164.0f, 84.0f }, VG_COLOR_RGB(0, 96, 180), VG_COLOR_RGB(255, 180, 0));
    vg_add_bezier_cubic_styled(&scene, md_float2_t{ 40.0f, 120.0f }, md_float2_t{ 80.0f, 84.0f }, md_float2_t{ 124.0f, 156.0f }, md_float2_t{ 164.0f, 112.0f }, &gradient_style);

    vg_path_move_to(&scene, md_float2_t{ 196.0f, 112.0f });
    vg_path_line_to(&scene, md_float2_t{ 236.0f, 112.0f });
    vg_path_bezier_cubic_curve_to(&scene, md_float2_t{ 240.0f, 120.0f }, md_float2_t{ 248.0f, 140.0f }, md_float2_t{ 252.0f, 148.0f });
    vg_path_close(&scene);
    const vg_style_t fill_style = { vg_paint_color(VG_COLOR_RGBA(40, 40, 40, 180)), vg_paint_none(), 1.0f, 4.0f, VG_LINE_CAP_BUTT, VG_LINE_JOIN_MITER, VG_FILL_RULE_NON_ZERO };
    vg_path_fill_styled(&scene, &fill_style);

    vg_add_text(&scene, md_float2_t{ 12.0f, 144.0f }, 14.0f, VG_COLOR_RGB(20, 20, 20), MD_STR_LIT("Vector <export>"));
    vg_add_text_rotated(&scene, md_float2_t{ 300.0f, 90.0f }, 12.0f, VG_COLOR_RGB(20, 20, 20), VG_TEXT_ANCHOR_MIDDLE, -90.0f, MD_STR_LIT("up"));

    const char* path = "viamd_test_svg_export.svg";
    ASSERT_TRUE(vg_scene_write_svg_file(&scene, md_str_from_cstr(path)));

    char buf[8192] = {0};
    FILE* f = fopen(path, "rb");
    ASSERT_TRUE(f != NULL);
    const size_t len = fread(buf, 1, sizeof(buf) - 1, f);
    fclose(f);
    remove(path);
    ASSERT_TRUE(len > 0);

    EXPECT_TRUE(strstr(buf, "<svg") != NULL);
    EXPECT_TRUE(strstr(buf, "<rect") != NULL);
    EXPECT_TRUE(strstr(buf, "rx=\"6\"") != NULL);
    EXPECT_TRUE(strstr(buf, "<path d=\"M 40 120 C") != NULL);
    EXPECT_TRUE(strstr(buf, "<path d=\"M 196 112 L 236 112 C") != NULL);
    EXPECT_TRUE(strstr(buf, "<linearGradient") != NULL);
    EXPECT_TRUE(strstr(buf, "url(#mdvg-grad-") != NULL);
    EXPECT_TRUE(strstr(buf, "Vector &lt;export&gt;") != NULL);
    EXPECT_TRUE(strstr(buf, "text-anchor=\"middle\"") != NULL);
    EXPECT_TRUE(strstr(buf, "rotate(-90") != NULL);

    md_arena_allocator_destroy(arena);
}
