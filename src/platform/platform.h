#pragma once

/*
@TODO:
        - Fix delta time computation for unix systems, this is not correct at al
        - Implement directory watch (+file watch) for different systems (win32, osx and linux)
        - Make sure that filters work with all systems for file_dialog
*/

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

#include "keys.h"
#include <core/common.h>
#include <core/types.h>
#include <core/array_types.h>
#include <core/string_types.h>

namespace platform {

const int MAX_KEYS = 512;
const int MAX_MOUSE_BUTTONS = 8;

struct Coordinate {
    float32 x;
    float32 y;
};

inline bool operator==(const Coordinate& a, const Coordinate& b) { return a.x == b.x && a.y == b.y; }
inline bool operator!=(const Coordinate& a, const Coordinate& b) { return a.x != b.x || a.y != b.y; }
inline Coordinate operator+(const Coordinate& a, const Coordinate& b) { return {a.x + b.x, a.y + b.y}; }
inline Coordinate operator-(const Coordinate& a, const Coordinate& b) { return {a.x - b.x, a.y - b.y}; }

struct Context {
    struct {
        const char* title;
        int32 width, height;
        bool vsync;
        bool should_close;
        void* ptr;
    } window;

    struct {
        int32 width;
        int32 height;
    } framebuffer;

    struct {
        struct {
            int major;
            int minor;
            int revision;
        } version;
    } gl_info;

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
			bool clicked[MAX_MOUSE_BUTTONS]; // This implies that the user did not move the mouse during the hit and the release of the button.

            Coordinate win_coord;  // Window coordinates
            Coordinate win_delta;

            Coordinate ndc_coord;  // Normalized device coordinates
			Coordinate ndc_delta;

            float32 scroll_delta;

            bool moving;
        } mouse;
    } input;

    struct {
        uint64 delta_ns;
        uint64 total_ns;

        float32 delta_s;
        float64 total_s;
    } timing;
};

// Context
bool initialize(Context* ctx, int32 width, int32 height, const char* title);
void shutdown(Context* ctx);
void update(Context* ctx);
void render_imgui(Context* ctx);
void swap_buffers(Context* ctx);

// Timing
typedef uint64 Timestamp;

void sleep(int32 milliseconds);
Timestamp get_time();
float compute_delta_ms(Timestamp t0, Timestamp t1);

// Filesystem
typedef StringBuffer<512> Path;
typedef int32 FileDialogFlags;
// typedef void (*DirectoryWatchCallback)(const DirectoryEntry entry, void* usr_data);

enum FileDialogFlags_ { FileDialogFlags_Open = BIT(0), FileDialogFlags_Save = BIT(1), FileDialogFlags_Directory = BIT(2) };

struct FileDialogResult {
    enum Result { Ok, Cancel };
    Path path;
    Result result;
};

struct DirectoryEntry {
    enum Type { File, Dir, Link, Unknown };
    Type type;
    Path name;
};

CString get_cwd();
DynamicArray<DirectoryEntry> list_directory(CString directory);
FileDialogResult file_dialog(FileDialogFlags flags, CString default_path = {}, CString filter = {});
/*
bool add_directory_watch(CString directory, DirectoryWatchCallback, void* usr_data = NULL);
bool remove_directory_watch(CString directory);
*/
}  // namespace platform
