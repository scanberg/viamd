#pragma once

#include <core/md_platform.h>
#include <core/md_str.h>
#include <core/md_vec_math.h>

#include "IconsFontAwesome6.h"

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
    FileDialog_Open = 0x1,
    FileDialog_Save = 0x2,
    FileDialog_Dir  = 0x4
};


// Opens a file system file dialogue which prompts the user to either open or save files/directories.
// returns true if successful and false if it fails or is canceled by the user.
// path_buf is a pointer to a string buffer which the null-terminated path is written to
// path_cap is the capacity of the string buffer
// flags represents the type of dialogue to be opened, e.g. FileDialog_Save to save a file, FileDialog_Open | FileDialog_Dir to open a directory.
// filter is a null-terminated string containing comma separated extensions to be applied as a filter for the files: e.g. "jpg,png,bmp" to limit the scope of files to files with endings .jpg, .png or .bmp
bool file_dialog(char* path_buf, int64_t path_cap, FileDialogFlags flags, const char* filter = 0);

}  // namespace application
