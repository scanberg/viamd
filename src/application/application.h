#pragma once

#include <core/md_platform.h>
#include <core/md_str.h>
#include <core/md_vec_math.h>

#include "IconsFontAwesome5.h"

namespace application {

typedef vec2_t Coordinate;

struct Context {
    struct {
        const char* title;
        int width, height;
        bool vsync;
        bool should_close;
        void* ptr;
    } window;

    struct {
        int width;
        int height;
    } framebuffer;

    struct {
        struct {
            int major;
            int minor;
            int revision;
        } version;
    } gl_info;

    struct {
        uint64_t delta_ns;
        uint64_t total_ns;

        double delta_s;
        double total_s;
    } timing;
};

// Context
bool initialize(Context* ctx, int64_t width, int64_t height, const char* title);
void shutdown(Context* ctx);
void update(Context* ctx);
void render_imgui(Context* ctx);
void swap_buffers(Context* ctx);

// File Dialog
typedef uint32_t FileDialogFlags;

enum FileDialogFlags_ {
    FileDialogFlags_Open = 0x1,
    FileDialogFlags_Save = 0x2,
    FileDialogFlags_Directory = 0x4
};

struct FileDialogResult {
    enum Result { Ok, Cancel };
    Result result;
    int64_t path_len;
    char path[512];
};

FileDialogResult file_dialog(FileDialogFlags flags, str_t default_path = {}, str_t filter = {});

}  // namespace application
