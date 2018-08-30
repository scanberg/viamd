#pragma once

#ifdef __APPLE__
#include "TargetConditionals.h"
#ifdef TARGET_OS_MAC
#define OS_MAC_OSX
#endif
#elif defined _WIN32 || defined _WIN64
#define OS_WINDOWS
#elif defined __linux__
#define OS_LINUX
#else
#error
#endif

#include <core/keys.h>
#include <core/types.h>
#include <core/array.h>
#include <core/string_utils.h>

namespace platform {

constexpr int MAX_KEYS = 512;
constexpr int MAX_MOUSE_BUTTONS = 8;

struct Context {
    struct {
        void* ptr;
        bool should_close;
        const char* title;
        bool vsync;
        int32 width, height;
    } window;

    struct {
        int32 width;
        int32 height;
    } framebuffer;

    struct {
        struct {
            bool down[MAX_KEYS];
            bool hit[MAX_KEYS];
            bool release[MAX_KEYS];
        } key;

        struct {
            bool down[MAX_MOUSE_BUTTONS];
            bool hit[MAX_MOUSE_BUTTONS];
            bool release[MAX_MOUSE_BUTTONS];

            vec2 coord_curr;
            vec2 coord_prev;

            vec2 ndc_curr;
            vec2 ndc_prev;

            vec2 velocity;
            vec2 scroll;
        } mouse;
    } input;

    struct {
        uint64 delta_ns;
        uint64 total_ns;

        float32 dt;
        float64 total_s;
    } timing;
};

void initialize(Context* ctx, int32 width, int32 height, const char* title);
void shutdown(Context* ctx);
void update(Context* ctx);
void swap_buffers(Context* ctx);

typedef StringBuffer<512> Path;

struct DirEntry {
    enum Type { File, Dir, Link, Unknown };
    Type type;
    Path name;
};

typedef uint64 Timestamp;

// Platform specific
void sleep(int32 milliseconds);
DynamicArray<DirEntry> list_directory(CString dir_path);
CString get_cwd();

Timestamp get_time();
float compute_delta_ms(Timestamp t0, Timestamp t1);

struct FileDialogResult {
    enum Result { FILE_OK, FILE_CANCEL };
    Path path;
    Result result;
};

typedef int FileDialogFlags;
enum FileDialogFlags_ { FileDialogFlags_Open = BIT(0), FileDialogFlags_Save = BIT(1), FileDialogFlags_Directory = BIT(2) };

FileDialogResult file_dialog(FileDialogFlags flags, CString default_path = {}, CString filter = {});

}  // namespace platform
