#pragma once

#ifdef __APPLE__
    #include "TargetConditionals.h"
    #ifdef TARGET_OS_MAC
        #define OS_MAC
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
        int width, height;
    } window;

	struct {
		int width;
		int height;
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
        float dt;
        uint64 delta_ns;
        uint64 total_ns;

		double total_s;
    } timing;
};

void initialize(Context* ctx, int width, int height, const char* title);
void shutdown(Context* ctx);
void update(Context* ctx);
void swap_buffers(Context* ctx);

typedef StringBuffer<512> Path;

struct DirEntry {
    enum Type {
        File,
        Dir,
        Link,
        Unknown
    };
    Type type;
    Path name;
};

// Filter example "jpg|png|tga"
DynamicArray<DirEntry> list_directory(CString dir_path, CString filter = {});
CString get_cwd();

Path open_file_dialog(CString filter = {});

}  // namespace platform

/*
Window* create_window(int width, int height, const char* window_title);
void	destroy_window(Window* window);
void	set_window_should_close(Window* window, bool value);
bool	window_in_focus(Window* window);
bool	window_should_close(Window* window);
void	get_framebuffer_size(Window* window, int* width, int* height);
*/
/*
InputState* get_input_state();
double get_delta_time();
void swap_buffers(Window* window);
void set_vsync(bool value);

void* perm_alloc(size_t);
void perm_free(void*);

void* scratch_malloc(size_t);
void scratch_free(void*);

}  // namespace platform

/*
struct KeyboardEvent {
    InputState* state;
    Key::Key_t key;

};

struct MouseEvent {

};
*/

// struct FileEvent {
//    CString path;
//};

/*
typedef void (*FileCallback)(const FileEvent& event);

void register_file_event_callback(FileCallback callback);
    


}
*/
