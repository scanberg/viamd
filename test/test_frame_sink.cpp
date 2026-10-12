#include "utest.h"

#include <frame_sink.h>

#include <core/md_os.h>

// The sink spawns processes and writes to a folder, which is done with POSIX calls here. Windows has its
// own implementation in frame_sink.cpp that these tests do not cover.
#if !defined(_WIN32)

// The decoder is implemented in image.cpp
#include <stb_image.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include <string>
#include <vector>

namespace {

std::string make_temp_dir() {
    char tmpl[] = "/tmp/viamd_sink_XXXXXX";
    return mkdtemp(tmpl) ? std::string(tmpl) : std::string();
}

bool write_text_file(const std::string& path, const std::string& text) {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) return false;
    const bool ok = fwrite(text.data(), 1, text.size(), f) == text.size();
    fclose(f);
    return ok;
}

bool read_binary_file(const std::string& path, std::vector<uint8_t>& out) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    uint8_t buf[4096];
    size_t n;
    out.clear();
    while ((n = fread(buf, 1, sizeof(buf), f)) > 0) out.insert(out.end(), buf, buf + n);
    fclose(f);
    return true;
}

// Stands in for ffmpeg: whatever arrives on stdin is stored where ffmpeg would have written the video
std::string make_fake_ffmpeg(const std::string& dir, const char* tail = "") {
    const std::string path = dir + "/fake_ffmpeg.sh";
    write_text_file(path, std::string("#!/bin/sh\nfor last; do :; done\ncat > \"$last\"\n") + tail + "\n");
    chmod(path.c_str(), 0755);
    return path;
}

str_t to_str(const std::string& s) { return str_t{s.c_str(), s.size()}; }

frame_sink::Status wait_until_done(frame_sink::Sink* sink) {
    frame_sink::Status st = frame_sink::status(sink);
    for (int i = 0; i < 2000 && !st.done; ++i) {
        md_thread_sleep(5);
        st = frame_sink::status(sink);
    }
    return st;
}

void remove_dir(const std::string& dir) {
    const std::string cmd = "rm -rf '" + dir + "'";
    if (system(cmd.c_str()) != 0) { /* nothing else to do */ }
}

}  // namespace

UTEST(viamd_frame_sink, png_sequence_is_written_upright_and_in_full) {
    const std::string dir = make_temp_dir();
    ASSERT_FALSE(dir.empty());

    const int W = 8, H = 4, N = 10;   // more frames than there are buffers, so acquire has to wait
    frame_sink::Desc desc;
    desc.kind = frame_sink::Kind::PngSequence;
    desc.dir = to_str(dir);
    desc.prefix = STR_LIT("shot");
    desc.width = W;
    desc.height = H;
    desc.max_in_flight = 3;

    char err[256] = {};
    frame_sink::Sink* sink = frame_sink::create(desc, err, sizeof(err));
    ASSERT_TRUE(sink != nullptr);

    for (int i = 0; i < N; ++i) {
        uint8_t* px = frame_sink::acquire(sink, true);
        ASSERT_TRUE(px != nullptr);
        memset(px, 0, (size_t)W * H * 4);
        for (int a = 0; a < W * H; ++a) px[a * 4 + 3] = 255;
        // OpenGL hands over the bottom row first: that row carries the frame number in red
        for (int x = 0; x < W; ++x) px[x * 4 + 0] = (uint8_t)(10 + i);
        frame_sink::submit(sink, px, i);
    }
    frame_sink::close(sink);

    const frame_sink::Status st = wait_until_done(sink);
    EXPECT_TRUE(st.done);
    EXPECT_TRUE(st.ok);
    EXPECT_EQ(N, st.written);
    EXPECT_EQ(0, st.failed);

    for (int i = 0; i < N; ++i) {
        char path[256];
        snprintf(path, sizeof(path), "%s/shot_%05d.png", dir.c_str(), i);
        int w = 0, h = 0, c = 0;
        uint8_t* img = stbi_load(path, &w, &h, &c, 4);
        ASSERT_TRUE(img != nullptr);
        EXPECT_EQ(W, w);
        EXPECT_EQ(H, h);
        EXPECT_EQ(0, (int)img[0]);                              // top row of the picture
        EXPECT_EQ(10 + i, (int)img[((H - 1) * W) * 4 + 0]);     // the bottom row ends up at the bottom
        stbi_image_free(img);
    }

    frame_sink::destroy(sink);
    remove_dir(dir);
}

UTEST(viamd_frame_sink, ffmpeg_arguments_describe_the_raw_input_and_the_output) {
    frame_sink::Desc desc;
    desc.kind = frame_sink::Kind::Ffmpeg;
    desc.dir = STR_LIT("/out dir");
    desc.prefix = STR_LIT("clip");
    desc.width = 1001;
    desc.height = 563;
    desc.fps = 30.0f;
    desc.crf = 20;

    std::vector<std::string> a;
    frame_sink::ffmpeg_arguments(desc, a);
    ASSERT_TRUE(a.size() > 10);
    EXPECT_STREQ("ffmpeg", a[0].c_str());
    EXPECT_STREQ("/out dir/clip.mp4", a.back().c_str());

    auto value_after = [&](const char* flag) -> std::string {
        for (size_t i = 0; i + 1 < a.size(); ++i) if (a[i] == flag) return a[i + 1];
        return "<missing>";
    };
    EXPECT_STREQ("rawvideo", value_after("-f").c_str());
    EXPECT_STREQ("rgba", value_after("-pixel_format").c_str());
    EXPECT_STREQ("1001x563", value_after("-video_size").c_str());
    EXPECT_STREQ("30", value_after("-framerate").c_str());
    EXPECT_STREQ("-", value_after("-i").c_str());
    EXPECT_STREQ("20", value_after("-crf").c_str());

    // Odd sizes have to be rounded to even for yuv420p, and the picture has to be turned right way up
    const std::string vf = value_after("-vf");
    EXPECT_TRUE(vf.find("vflip") != std::string::npos);
    EXPECT_TRUE(vf.find("trunc(iw/2)*2") != std::string::npos);

    desc.ffmpeg = STR_LIT("/opt/ff/ffmpeg");
    frame_sink::ffmpeg_arguments(desc, a);
    EXPECT_STREQ("/opt/ff/ffmpeg", a[0].c_str());
}

UTEST(viamd_frame_sink, ffmpeg_receives_every_frame_in_order) {
    const std::string dir = make_temp_dir();
    ASSERT_FALSE(dir.empty());
    const std::string fake = make_fake_ffmpeg(dir);

    const int W = 8, H = 4, N = 12;
    const size_t frame_bytes = (size_t)W * H * 4;

    frame_sink::Desc desc;
    desc.kind = frame_sink::Kind::Ffmpeg;
    desc.dir = to_str(dir);
    desc.prefix = STR_LIT("movie");
    desc.width = W;
    desc.height = H;
    desc.max_in_flight = 3;
    desc.ffmpeg = to_str(fake);

    char err[256] = {};
    frame_sink::Sink* sink = frame_sink::create(desc, err, sizeof(err));
    ASSERT_TRUE(sink != nullptr);

    for (int i = 0; i < N; ++i) {
        uint8_t* px = frame_sink::acquire(sink, true);
        ASSERT_TRUE(px != nullptr);
        memset(px, i + 1, frame_bytes);
        frame_sink::submit(sink, px, i);
    }
    frame_sink::close(sink);

    const frame_sink::Status st = wait_until_done(sink);
    EXPECT_TRUE(st.done);
    EXPECT_TRUE(st.ok);
    EXPECT_EQ(N, st.written);

    std::vector<uint8_t> out;
    ASSERT_TRUE(read_binary_file(dir + "/movie.mp4", out));
    ASSERT_EQ((size_t)N * frame_bytes, out.size());
    for (int i = 0; i < N; ++i) {
        EXPECT_EQ(i + 1, (int)out[(size_t)i * frame_bytes]);
        EXPECT_EQ(i + 1, (int)out[(size_t)(i + 1) * frame_bytes - 1]);
    }

    frame_sink::destroy(sink);
    remove_dir(dir);
}

UTEST(viamd_frame_sink, ffmpeg_failing_is_reported) {
    const std::string dir = make_temp_dir();
    ASSERT_FALSE(dir.empty());
    const std::string fake = make_fake_ffmpeg(dir, "exit 3");

    frame_sink::Desc desc;
    desc.kind = frame_sink::Kind::Ffmpeg;
    desc.dir = to_str(dir);
    desc.prefix = STR_LIT("movie");
    desc.width = 4;
    desc.height = 4;
    desc.ffmpeg = to_str(fake);

    char err[256] = {};
    frame_sink::Sink* sink = frame_sink::create(desc, err, sizeof(err));
    ASSERT_TRUE(sink != nullptr);
    uint8_t* px = frame_sink::acquire(sink, true);
    ASSERT_TRUE(px != nullptr);
    memset(px, 0, 4 * 4 * 4);
    frame_sink::submit(sink, px, 0);
    frame_sink::close(sink);

    const frame_sink::Status st = wait_until_done(sink);
    EXPECT_TRUE(st.done);
    EXPECT_FALSE(st.ok);
    EXPECT_TRUE(strstr(st.message, "code 3") != nullptr);

    frame_sink::destroy(sink);
    remove_dir(dir);
}

UTEST(viamd_frame_sink, a_missing_ffmpeg_is_an_error_not_a_crash) {
    const std::string dir = make_temp_dir();
    ASSERT_FALSE(dir.empty());

    frame_sink::Desc desc;
    desc.kind = frame_sink::Kind::Ffmpeg;
    desc.dir = to_str(dir);
    desc.prefix = STR_LIT("movie");
    desc.width = 4;
    desc.height = 4;
    desc.ffmpeg = STR_LIT("/nonexistent/path/to/ffmpeg");

    char err[256] = {};
    EXPECT_TRUE(frame_sink::create(desc, err, sizeof(err)) == nullptr);
    EXPECT_TRUE(err[0] != '\0');

    remove_dir(dir);
}

UTEST(viamd_frame_sink, cancel_returns_every_buffer_and_finishes) {
    const std::string dir = make_temp_dir();
    ASSERT_FALSE(dir.empty());

    frame_sink::Desc desc;
    desc.kind = frame_sink::Kind::PngSequence;
    desc.dir = to_str(dir);
    desc.prefix = STR_LIT("c");
    desc.width = 64;
    desc.height = 64;
    desc.max_in_flight = 4;

    char err[256] = {};
    frame_sink::Sink* sink = frame_sink::create(desc, err, sizeof(err));
    ASSERT_TRUE(sink != nullptr);

    for (int i = 0; i < 4; ++i) {
        uint8_t* px = frame_sink::acquire(sink, true);
        ASSERT_TRUE(px != nullptr);
        memset(px, 7, 64 * 64 * 4);
        frame_sink::submit(sink, px, i);
    }
    frame_sink::cancel(sink);

    const frame_sink::Status st = wait_until_done(sink);
    EXPECT_TRUE(st.done);
    EXPECT_LE(st.written + st.failed, 4);

    frame_sink::destroy(sink);
    remove_dir(dir);
}

UTEST(viamd_frame_sink, no_buffer_is_handed_out_twice) {
    const std::string dir = make_temp_dir();
    ASSERT_FALSE(dir.empty());

    frame_sink::Desc desc;
    desc.kind = frame_sink::Kind::PngSequence;
    desc.dir = to_str(dir);
    desc.prefix = STR_LIT("p");
    desc.width = 4;
    desc.height = 4;
    desc.max_in_flight = 2;

    char err[256] = {};
    frame_sink::Sink* sink = frame_sink::create(desc, err, sizeof(err));
    ASSERT_TRUE(sink != nullptr);

    uint8_t* a = frame_sink::acquire(sink, false);
    uint8_t* b = frame_sink::acquire(sink, false);
    EXPECT_TRUE(a != nullptr);
    EXPECT_TRUE(b != nullptr);
    EXPECT_TRUE(a != b);
    EXPECT_TRUE(frame_sink::acquire(sink, false) == nullptr);   // the pool is exhausted, the caller must wait

    frame_sink::cancel(sink);
    frame_sink::destroy(sink);
    remove_dir(dir);
}

#endif  // !_WIN32

UTEST(viamd_frame_sink, the_codec_sets_the_encoder_and_the_file) {
    frame_sink::Desc desc;
    desc.kind = frame_sink::Kind::Ffmpeg;
    desc.dir = STR_LIT("/out");
    desc.prefix = STR_LIT("clip");
    desc.width = 640;
    desc.height = 480;
    desc.crf = 24;

    auto value_after = [](const std::vector<std::string>& a, const char* flag) -> std::string {
        for (size_t i = 0; i + 1 < a.size(); ++i) if (a[i] == flag) return a[i + 1];
        return "<missing>";
    };
    auto has = [](const std::vector<std::string>& a, const char* flag) {
        for (const std::string& s : a) if (s == flag) return true;
        return false;
    };

    std::vector<std::string> a;
    desc.codec = frame_sink::Codec::H264;
    frame_sink::ffmpeg_arguments(desc, a);
    EXPECT_STREQ("libx264", value_after(a, "-c:v").c_str());
    EXPECT_STREQ("/out/clip.mp4", a.back().c_str());
    EXPECT_TRUE(has(a, "-movflags"));

    desc.codec = frame_sink::Codec::H265;
    frame_sink::ffmpeg_arguments(desc, a);
    EXPECT_STREQ("libx265", value_after(a, "-c:v").c_str());
    EXPECT_STREQ("hvc1", value_after(a, "-tag:v").c_str());   /* so that QuickTime plays it */
    EXPECT_STREQ("24", value_after(a, "-crf").c_str());
    EXPECT_STREQ("/out/clip.mp4", a.back().c_str());

    desc.codec = frame_sink::Codec::Vp9;
    frame_sink::ffmpeg_arguments(desc, a);
    EXPECT_STREQ("libvpx-vp9", value_after(a, "-c:v").c_str());
    EXPECT_STREQ("0", value_after(a, "-b:v").c_str());        /* the crf is the only quality control */
    EXPECT_FALSE(has(a, "-movflags"));                         /* an mp4 option */
    EXPECT_STREQ("/out/clip.webm", a.back().c_str());
    EXPECT_STREQ("webm", frame_sink::file_extension(frame_sink::Codec::Vp9));
    EXPECT_STREQ("mp4", frame_sink::file_extension(frame_sink::Codec::H265));
}
