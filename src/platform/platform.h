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
//#include <core/common.h>
//#include <core/types.h>
//#include <core/array_types.h>
//#include <core/string_types.h>
#include <core/md_str.h>

#include "IconsFontAwesome5.h"

namespace platform {

const int MAX_KEYS = 512;
const int MAX_MOUSE_BUTTONS = 8;

struct Coordinate {
    float x;
    float y;
};

inline bool operator==(const Coordinate& a, const Coordinate& b) { return a.x == b.x && a.y == b.y; }
inline bool operator!=(const Coordinate& a, const Coordinate& b) { return a.x != b.x || a.y != b.y; }
inline Coordinate operator+(const Coordinate& a, const Coordinate& b) { return {a.x + b.x, a.y + b.y}; }
inline Coordinate operator-(const Coordinate& a, const Coordinate& b) { return {a.x - b.x, a.y - b.y}; }

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

            float scroll_delta;

            bool moving;
        } mouse;
    } input;

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

// Timing
typedef uint64_t Timestamp;

void sleep(int64_t milliseconds);
Timestamp get_time();
double compute_delta_ms(Timestamp t0, Timestamp t1);
double compute_delta_s(Timestamp t0, Timestamp t1);

// Filesystem
typedef uint32_t FileDialogFlags;

enum FileDialogFlags_ { FileDialogFlags_Open = 0x1, FileDialogFlags_Save = 0x2, FileDialogFlags_Directory = 0x4 };

struct FileDialogResult {
    enum Result { Ok, Cancel };
    Result result;
    int64_t path_len;
    char path[512];
};

/*
struct DirectoryEntry {
    enum Type { File, Dir, Link, Unknown };
    Type type = Unknown;
    Path name = {};
};

str_t get_cwd();
DynamicArray<DirectoryEntry> list_directory(CStringView directory);
*/
FileDialogResult file_dialog(FileDialogFlags flags, str_t default_path = {}, str_t filter = {});

/*
bool add_directory_watch(CString directory, DirectoryWatchCallback, void* usr_data = NULL);
bool remove_directory_watch(CString directory);
*/
}  // namespace platform
