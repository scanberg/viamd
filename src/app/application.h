#pragma once

#include <stdint.h>
#include <stddef.h>

#include <core/md_str.h>
#include <core/md_gpu.h>

namespace application {

typedef void (*FileDropCallback)(size_t file_count, const md_str_t file_paths[], void* user_data);

struct Context {
    struct {
        md_str_t title;
        unsigned int width, height;
        float scale_factor;
        bool vsync;
        bool should_close;
        void* ptr;
    } window;

    struct {
        uint32_t width;
        uint32_t height;
    } framebuffer;

    // What everything is drawn with, and how the window is shown: by md_gpu (a surface), or through OpenGL (gl)
    struct {
        md_gpu_device_t device;
        md_gpu_stream_t stream;     // GRAPHICS: draws the window and presents it
        bool gl;
    } gpu;

    struct {
        double delta_s;
        double total_s;
    } timing;

    struct {
        FileDropCallback callback;
        void* user_data;
    } file_drop;
};

// Context
bool initialize(Context* ctx, size_t width, size_t height, md_str_t title);
void shutdown(Context* ctx);
void update(Context* ctx);
// ImGui's frame into the window (and the other windows it has)
void render_imgui(Context* ctx);
void swap_buffers(Context* ctx);
// The window as drawn this frame, the GUI and all, after render_imgui: RGBA8, rows from the top, of the framebuffer's
// size (width, height). False when it cannot be read.
bool read_frame(Context* ctx, uint32_t* rgba, uint32_t width, uint32_t height);

// File Dialog
typedef unsigned int FileDialogFlag;

enum {
    FileDialogFlag_Open = 0x1,
    FileDialogFlag_Save = 0x2,
    FileDialogFlag_Dir  = 0x4
};

// Opens a file system file dialogue which prompts the user to either open or save files/directories.
// returns true if successful and false if it fails or is canceled by the user.
// path_buf is a pointer to a string buffer which the null-terminated path is written to
// path_cap is the capacity of the string buffer
// flags represents the type of dialogue to be opened, e.g. FileDialog_Save to save a file, FileDialog_Open | FileDialog_Dir to open a directory.
// filter is a null-terminated string containing comma separated extensions to be applied as a filter for the files: e.g. "jpg,png,bmp" to limit the scope of files to files with endings .jpg, .png or .bmp
bool file_dialog(char* path_buf, size_t path_cap, FileDialogFlag flags, md_str_t filter = {});

}  // namespace application
