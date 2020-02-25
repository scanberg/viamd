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

        float delta_s;
        float total_s;
    } timing;
};

// Context
bool initialize(Context* ctx, i32 width, i32 height, const char* title);
void shutdown(Context* ctx);
void update(Context* ctx);
void render_imgui(Context* ctx);
void swap_buffers(Context* ctx);

// Timing
typedef uint64_t Timestamp;

void sleep(int32_t milliseconds);
Timestamp get_time();
float compute_delta_ms(Timestamp t0, Timestamp t1);

// Filesystem
typedef StringBuffer<512> Path;
typedef int32_t FileDialogFlags;

enum FileDialogFlags_ { FileDialogFlags_Open = BIT(0), FileDialogFlags_Save = BIT(1), FileDialogFlags_Directory = BIT(2) };

struct FileDialogResult {
    enum Result { Ok, Cancel };
    Result result;
    Path path;
};

struct DirectoryEntry {
    enum Type { File, Dir, Link, Unknown };
    Type type = Unknown;
    Path name = {};
};

CStringView get_cwd();
DynamicArray<DirectoryEntry> list_directory(CStringView directory);
FileDialogResult file_dialog(FileDialogFlags flags, CStringView default_path = {}, CStringView filter = {});

// Atomic operations
int32_t atomic_fetch_and_add(volatile int32_t* ptr, int32_t add);
uint32_t atomic_fetch_and_add(volatile uint32_t* ptr, uint32_t add);

int64_t atomic_fetch_and_add(volatile int64_t* ptr, int64_t add);
uint64_t atomic_fetch_and_add(volatile uint64_t* ptr, uint64_t add);

/*
bool add_directory_watch(CString directory, DirectoryWatchCallback, void* usr_data = NULL);
bool remove_directory_watch(CString directory);
*/
}  // namespace platform
