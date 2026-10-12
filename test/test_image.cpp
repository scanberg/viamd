#include "utest.h"

#include <image.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

UTEST(viamd_image, a_written_png_decodes_to_the_same_pixels) {
    const int W = 3, H = 2;
    const uint8_t src[W * H * 4] = {
        255, 0, 0, 255,   0, 255, 0, 128,   0, 0, 255, 0,
        10, 20, 30, 255,  40, 50, 60, 200,  70, 80, 90, 1,
    };

    char path[] = "/tmp/viamd_image_testXXXXXX";
    const int fd = mkstemp(path);
    ASSERT_TRUE(fd >= 0);
    fclose(fdopen(fd, "wb"));
    ASSERT_TRUE(image_write_png(str_from_cstr(path), src, W, H));

    FILE* f = fopen(path, "rb");
    ASSERT_TRUE(f != nullptr);
    std::vector<uint8_t> bytes;
    uint8_t buf[256];
    size_t n;
    while ((n = fread(buf, 1, sizeof(buf), f)) > 0) bytes.insert(bytes.end(), buf, buf + n);
    fclose(f);
    remove(path);

    int w = 0, h = 0;
    uint8_t* pixels = image_decode_rgba(bytes.data(), bytes.size(), &w, &h);
    ASSERT_TRUE(pixels != nullptr);
    EXPECT_EQ(W, w);
    EXPECT_EQ(H, h);
    EXPECT_EQ(0, memcmp(src, pixels, sizeof(src)));
    image_free(pixels);
}

UTEST(viamd_image, something_that_is_not_an_image_is_refused) {
    const char text[] = "this is not a picture";
    int w = 0, h = 0;
    EXPECT_TRUE(image_decode_rgba(text, sizeof(text), &w, &h) == nullptr);
}

UTEST(viamd_image, the_bounds_of_the_visible_part_leave_out_a_transparent_margin) {
    const int W = 5, H = 4;
    std::vector<uint8_t> img((size_t)W * H * 4, 0);
    auto set = [&](int x, int y, uint8_t a) { img[((size_t)y * W + x) * 4 + 3] = a; };
    set(1, 1, 255);
    set(3, 2, 100);
    int x0, y0, x1, y1;
    ASSERT_TRUE(image_alpha_bounds(img.data(), W, H, 8, &x0, &y0, &x1, &y1));
    EXPECT_EQ(1, x0);
    EXPECT_EQ(1, y0);
    EXPECT_EQ(4, x1);
    EXPECT_EQ(3, y1);
}

UTEST(viamd_image, an_alpha_at_or_below_the_threshold_is_not_visible) {
    std::vector<uint8_t> img(4 * 4, 0);
    img[3] = 8;
    int x0, y0, x1, y1;
    EXPECT_FALSE(image_alpha_bounds(img.data(), 2, 2, 8, &x0, &y0, &x1, &y1));
    img[3] = 9;
    EXPECT_TRUE(image_alpha_bounds(img.data(), 2, 2, 8, &x0, &y0, &x1, &y1));
}

UTEST(viamd_image, transparent_pixels_take_the_average_colour_of_the_visible_ones) {
    uint8_t img[3 * 4] = {
        200, 100, 0, 255,
        0, 0, 0, 0,
        100, 200, 40, 255,
    };
    image_bleed_transparent(img, 3, 1);
    EXPECT_EQ(150, (int)img[4]);
    EXPECT_EQ(150, (int)img[5]);
    EXPECT_EQ(20, (int)img[6]);
    EXPECT_EQ(0, (int)img[7]);       /* it stays transparent */
    EXPECT_EQ(200, (int)img[0]);     /* the visible ones are not touched */
}

UTEST(viamd_image, an_image_with_nothing_visible_is_left_as_it_is) {
    uint8_t img[8] = {5, 6, 7, 0, 8, 9, 10, 0};
    image_bleed_transparent(img, 2, 1);
    EXPECT_EQ(5, (int)img[0]);
    EXPECT_EQ(9, (int)img[5]);
}

UTEST(viamd_image, light_grey_becomes_the_colour_asked_for_and_keeps_its_alpha) {
    uint8_t img[4 * 4] = {
        224, 224, 224, 255,     /* light grey: replaced */
        226, 222, 224, 90,      /* nearly grey, partly transparent: replaced, alpha kept */
        200, 40, 90, 255,       /* coloured: kept */
        60, 60, 60, 255,        /* dark grey: kept */
    };
    image_replace_light_grey(img, 4, 1, 128, 24, 0, 0, 0);
    EXPECT_EQ(0, (int)img[0]); EXPECT_EQ(0, (int)img[1]); EXPECT_EQ(0, (int)img[2]); EXPECT_EQ(255, (int)img[3]);
    EXPECT_EQ(0, (int)img[4]); EXPECT_EQ(0, (int)img[6]); EXPECT_EQ(90, (int)img[7]);
    EXPECT_EQ(200, (int)img[8]); EXPECT_EQ(40, (int)img[9]); EXPECT_EQ(90, (int)img[10]);
    EXPECT_EQ(60, (int)img[12]); EXPECT_EQ(60, (int)img[13]);
}

UTEST(viamd_image, transparent_pixels_are_not_recoloured) {
    uint8_t img[4] = {230, 230, 230, 0};
    image_replace_light_grey(img, 1, 1, 128, 24, 0, 0, 0);
    EXPECT_EQ(230, (int)img[0]);
}
