#include "keys.h"
#include "types.h"
#include <GL/gl3w.h>

namespace platform {

struct Window;

struct InputState {
    bool key_down[512];
    bool key_hit[512];
    bool key_release[512];

    bool mouse_down[8];
    bool mouse_hit[8];
    bool mouse_release[8];

    ivec2 mouse_coord;
	ivec2 mouse_velocity;
    ivec2 mouse_scroll;
};

void        initialize();
Window*     create_window(int width, int height, const char* window_title);
void        destroy_window(Window* window);
void        shutdown();
void        update();
bool        window_in_focus(Window* window);
bool        window_should_close(Window* window);
void        get_framebuffer_size(Window* window, int* width, int* height);
InputState* get_input_state();
void        swap_buffers(Window* window);

}


