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

#include "common.h"
#include "keys.h"
#include "types.h"
#include "array_types.h"
#include "string_types.h"

namespace platform {

const int MAX_KEYS = 512;
const int MAX_MOUSE_BUTTONS = 8;

typedef StringBuffer<512> Path;
typedef uint64 Timestamp;

typedef int32 FileDialogFlags;
enum FileDialogFlags_ { FileDialogFlags_Open = BIT(0), FileDialogFlags_Save = BIT(1), FileDialogFlags_Directory = BIT(2) };

struct FileDialogResult {
    enum Result { FILE_OK, FILE_CANCEL };
    Path path;
    Result result;
};

struct DirectoryEntry {
    enum Type { File, Dir, Link, Unknown };
    Type type;
    Path name;
};

struct Coordinate {
    float32 x;
    float32 y;
};

inline bool operator==(const Coordinate& a, const Coordinate& b) { return a.x == b.x && a.y == b.y; }
inline bool operator!=(const Coordinate& a, const Coordinate& b) { return a.x != b.x || a.y != b.y; }

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

            Coordinate coord;      // Window coordinates
            Coordinate ndc_coord;  // Normalized device coordinates
            /*
struct {
    Coordinate window;
    struct {
        Coordinate current;
        Coordinate previous;
    } window;

    struct {
        Coordinate current;
        Coordinate previous;
    } ndc;
} coordinate;
            */
            float32 scroll_delta;
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
void swap_buffers(Context* ctx);

// Timing
void sleep(int32 milliseconds);
Timestamp get_time();
float compute_delta_ms(Timestamp t0, Timestamp t1);

// Filesystem
DynamicArray<DirectoryEntry> list_directory(CString directory);
CString get_cwd();
FileDialogResult file_dialog(FileDialogFlags flags, CString default_path = {}, CString filter = {});

}  // namespace platform
