#include <core/keys.h>
#include <core/types.h>
#include <GL/gl3w.h>

namespace platform {

constexpr int MaxKeys = 512;
constexpr int MaxMouseButtons = 8;

struct Window;

struct InputState {
    bool key_down[MaxKeys];
    bool key_hit[MaxKeys];
    bool key_release[MaxKeys];

    bool mouse_down[MaxMouseButtons];
    bool mouse_hit[MaxMouseButtons];
    bool mouse_release[MaxMouseButtons];

    ivec2 mouse_coord;
	ivec2 mouse_velocity;
    ivec2 mouse_scroll;
};

void        initialize();
Window*     create_window(int width, int height, const char* window_title);
void        destroy_window(Window* window);
void		set_window_should_close(Window* window, bool value);
void        shutdown();
void        update();
bool        window_in_focus(Window* window);
bool        window_should_close(Window* window);
void        get_framebuffer_size(Window* window, int* width, int* height);
InputState* get_input_state();
void        swap_buffers(Window* window);

}


