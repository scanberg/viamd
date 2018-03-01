#pragma once

#include <core/keys.h>
#include <core/types.h>
#include <core/allocator.h>
//#include <core/string_utils.h>

struct DefaultAllocator;
struct FrameAllocator;

namespace platform {

constexpr int MAX_KEYS = 512;
constexpr int MAX_MOUSE_BUTTONS = 8;

struct Window;

struct InputState {
    bool key_down[MAX_KEYS];
    bool key_hit[MAX_KEYS];
    bool key_release[MAX_KEYS];

    bool mouse_down[MAX_MOUSE_BUTTONS];
    bool mouse_hit[MAX_MOUSE_BUTTONS];
    bool mouse_release[MAX_MOUSE_BUTTONS];

    vec2 mouse_screen_coords;
    vec2 prev_mouse_screen_coords;
    vec2 mouse_ndc_coords;
    vec2 prev_mouse_ndc_coords;
    vec2 mouse_velocity;
    vec2 mouse_scroll;
};

void initialize();
Window* create_window(int width, int height, const char* window_title);
void destroy_window(Window* window);
void set_window_should_close(Window* window, bool value);
void shutdown();
void update();
bool window_in_focus(Window* window);
bool window_should_close(Window* window);
void get_framebuffer_size(Window* window, int* width, int* height);
InputState* get_input_state();
double get_delta_time();
void swap_buffers(Window* window);
void set_vsync(bool value);

DefaultAllocator* default_allocator();
FrameAllocator* frame_allocator();

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
