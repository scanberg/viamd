#include <core/platform.h>
#include <core/gl.h>
#include <GLFW/glfw3.h>
#ifdef _WIN32
#undef APIENTRY
#define GLFW_EXPOSE_NATIVE_WIN32
#define GLFW_EXPOSE_NATIVE_WGL
#include <GLFW/glfw3native.h>
#endif
#include <imgui.h>
#include <stdio.h>
#include <stdlib.h>

namespace platform {

struct Window {
    GLFWwindow* glfw_window = NULL;
};

// TODO: Perhaps extend this to array of windows
static Window     g_window;
static InputState g_input_state;

// Data
static double       g_time = 0.0f;
static GLuint       g_font_texture = 0;
static int          g_shader_handle = 0, g_vert_handle = 0, g_frag_handle = 0;
static int          g_attrib_location_tex = 0, g_attrib_location_proj_mat = 0;
static int          g_attrib_location_position = 0, g_attrib_location_uv = 0, g_attrib_location_color = 0;
static unsigned int g_vbo_handle = 0, g_vao_handle = 0, g_elements_handle = 0;

static void error_callback(int error, const char* description) {
    fprintf(stderr, "Error %d: %s\n", error, description);
}

static void gl_callback(
	GLenum source,
	GLenum type,
	GLuint id,
	GLenum severity,
	GLsizei length,
	const GLchar* message,
	const void* userParam ) {
	(void)source; (void)type; (void)id;
	(void)severity; (void)length; (void)userParam;

	if (severity >= GL_DEBUG_SEVERITY_LOW)
		fprintf(stderr, "%s\n\n", message);
	if (severity == GL_DEBUG_SEVERITY_HIGH) {
		fprintf(stderr, "Aborting...\n");
		abort();
	}
}

static void mouse_button_callback(GLFWwindow*, int button, int action, int /*mods*/) {
    if (action == GLFW_PRESS) {
        g_input_state.mouse_down[button] = true;
        g_input_state.mouse_hit[button] = true;
    } else if (action == GLFW_RELEASE) {
        g_input_state.mouse_down[button] = false;
        g_input_state.mouse_release[button] = true;
    }
}

static void mouse_scroll_callback(GLFWwindow*, double xoffset, double yoffset) {
    ImGuiIO& io = ImGui::GetIO();

    g_input_state.mouse_scroll.y = (float)yoffset;
    io.MouseWheelH += (float)xoffset;
    io.MouseWheel += (float)yoffset;
}

static void key_callback(GLFWwindow*, int key, int, int action, int mods) {
    ImGuiIO& io = ImGui::GetIO();

    if (action == GLFW_PRESS) {
        io.KeysDown[key] = true;
        g_input_state.key_down[key] = true;
        g_input_state.key_hit[key] = true;
    }
    if (action == GLFW_RELEASE) {
        io.KeysDown[key] = false;
        g_input_state.key_down[key] = false;
    }

    // IMGUI
    (void)mods; // Modifiers are not reliable across systems
    io.KeyCtrl = io.KeysDown[GLFW_KEY_LEFT_CONTROL] || io.KeysDown[GLFW_KEY_RIGHT_CONTROL];
    io.KeyShift = io.KeysDown[GLFW_KEY_LEFT_SHIFT] || io.KeysDown[GLFW_KEY_RIGHT_SHIFT];
    io.KeyAlt = io.KeysDown[GLFW_KEY_LEFT_ALT] || io.KeysDown[GLFW_KEY_RIGHT_ALT];
    io.KeySuper = io.KeysDown[GLFW_KEY_LEFT_SUPER] || io.KeysDown[GLFW_KEY_RIGHT_SUPER];
}

static void char_callback(GLFWwindow*, unsigned int c) {
    ImGuiIO& io = ImGui::GetIO();
    if (c > 0 && c < 0x10000)
        io.AddInputCharacter((unsigned short)c);
}


// This is the main rendering function that you have to implement and provide to ImGui (via setting up 'RenderDrawListsFn' in the ImGuiIO structure)
// Note that this implementation is little overcomplicated because we are saving/setting up/restoring every OpenGL state explicitly, in order to be able to run within any OpenGL engine that doesn't do so. 
// If text or lines are blurry when integrating ImGui in your engine: in your Render function, try translating your projection matrix by (0.5f,0.5f) or (0.375f,0.375f)
static void imgui_render_draw_lists(ImDrawData* draw_data)
{
    // Avoid rendering when minimized, scale coordinates for retina displays (screen coordinates != framebuffer coordinates)
    ImGuiIO& io = ImGui::GetIO();
    int fb_width = (int)(io.DisplaySize.x * io.DisplayFramebufferScale.x);
    int fb_height = (int)(io.DisplaySize.y * io.DisplayFramebufferScale.y);
    if (fb_width == 0 || fb_height == 0)
        return;
    draw_data->ScaleClipRects(io.DisplayFramebufferScale);

    // Backup GL state
    GLenum last_active_texture; glGetIntegerv(GL_ACTIVE_TEXTURE, (GLint*)&last_active_texture);
    glActiveTexture(GL_TEXTURE0);
    GLint last_program; glGetIntegerv(GL_CURRENT_PROGRAM, &last_program);
    GLint last_texture; glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture);
    GLint last_sampler; glGetIntegerv(GL_SAMPLER_BINDING, &last_sampler);
    GLint last_array_buffer; glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &last_array_buffer);
    GLint last_element_array_buffer; glGetIntegerv(GL_ELEMENT_ARRAY_BUFFER_BINDING, &last_element_array_buffer);
    GLint last_vertex_array; glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &last_vertex_array);
    GLint last_polygon_mode[2]; glGetIntegerv(GL_POLYGON_MODE, last_polygon_mode);
    GLint last_viewport[4]; glGetIntegerv(GL_VIEWPORT, last_viewport);
    GLint last_scissor_box[4]; glGetIntegerv(GL_SCISSOR_BOX, last_scissor_box);
    GLenum last_blend_src_rgb; glGetIntegerv(GL_BLEND_SRC_RGB, (GLint*)&last_blend_src_rgb);
    GLenum last_blend_dst_rgb; glGetIntegerv(GL_BLEND_DST_RGB, (GLint*)&last_blend_dst_rgb);
    GLenum last_blend_src_alpha; glGetIntegerv(GL_BLEND_SRC_ALPHA, (GLint*)&last_blend_src_alpha);
    GLenum last_blend_dst_alpha; glGetIntegerv(GL_BLEND_DST_ALPHA, (GLint*)&last_blend_dst_alpha);
    GLenum last_blend_equation_rgb; glGetIntegerv(GL_BLEND_EQUATION_RGB, (GLint*)&last_blend_equation_rgb);
    GLenum last_blend_equation_alpha; glGetIntegerv(GL_BLEND_EQUATION_ALPHA, (GLint*)&last_blend_equation_alpha);
    GLboolean last_enable_blend = glIsEnabled(GL_BLEND);
    GLboolean last_enable_cull_face = glIsEnabled(GL_CULL_FACE);
    GLboolean last_enable_depth_test = glIsEnabled(GL_DEPTH_TEST);
    GLboolean last_enable_scissor_test = glIsEnabled(GL_SCISSOR_TEST);

    // Setup render state: alpha-blending enabled, no face culling, no depth testing, scissor enabled, polygon fill
    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_SCISSOR_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // Setup viewport, orthographic projection matrix
    glViewport(0, 0, (GLsizei)fb_width, (GLsizei)fb_height);
    const float ortho_projection[4][4] =
    {
        { 2.0f/io.DisplaySize.x, 0.0f,                   0.0f, 0.0f },
        { 0.0f,                  2.0f/-io.DisplaySize.y, 0.0f, 0.0f },
        { 0.0f,                  0.0f,                  -1.0f, 0.0f },
        {-1.0f,                  1.0f,                   0.0f, 1.0f },
    };
    glUseProgram(g_shader_handle);
    glUniform1i(g_attrib_location_tex, 0);
    glUniformMatrix4fv(g_attrib_location_proj_mat, 1, GL_FALSE, &ortho_projection[0][0]);
    glBindVertexArray(g_vao_handle);
    glBindSampler(0, 0); // Rely on combined texture/sampler state.

    for (int n = 0; n < draw_data->CmdListsCount; n++)
    {
        const ImDrawList* cmd_list = draw_data->CmdLists[n];
        const ImDrawIdx* idx_buffer_offset = 0;

        glBindBuffer(GL_ARRAY_BUFFER, g_vbo_handle);
        glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)cmd_list->VtxBuffer.Size * sizeof(ImDrawVert), (const GLvoid*)cmd_list->VtxBuffer.Data, GL_STREAM_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_elements_handle);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, (GLsizeiptr)cmd_list->IdxBuffer.Size * sizeof(ImDrawIdx), (const GLvoid*)cmd_list->IdxBuffer.Data, GL_STREAM_DRAW);

        for (int cmd_i = 0; cmd_i < cmd_list->CmdBuffer.Size; cmd_i++)
        {
            const ImDrawCmd* pcmd = &cmd_list->CmdBuffer[cmd_i];
            if (pcmd->UserCallback)
            {
                pcmd->UserCallback(cmd_list, pcmd);
            }
            else
            {
                glBindTexture(GL_TEXTURE_2D, (GLuint)(intptr_t)pcmd->TextureId);
                glScissor((int)pcmd->ClipRect.x, (int)(fb_height - pcmd->ClipRect.w), (int)(pcmd->ClipRect.z - pcmd->ClipRect.x), (int)(pcmd->ClipRect.w - pcmd->ClipRect.y));
                glDrawElements(GL_TRIANGLES, (GLsizei)pcmd->ElemCount, sizeof(ImDrawIdx) == 2 ? GL_UNSIGNED_SHORT : GL_UNSIGNED_INT, idx_buffer_offset);
            }
            idx_buffer_offset += pcmd->ElemCount;
        }
    }

    // Restore modified GL state
    glUseProgram(last_program);
    glBindTexture(GL_TEXTURE_2D, last_texture);
    glBindSampler(0, last_sampler);
    glActiveTexture(last_active_texture);
    glBindVertexArray(last_vertex_array);
    glBindBuffer(GL_ARRAY_BUFFER, last_array_buffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, last_element_array_buffer);
    glBlendEquationSeparate(last_blend_equation_rgb, last_blend_equation_alpha);
    glBlendFuncSeparate(last_blend_src_rgb, last_blend_dst_rgb, last_blend_src_alpha, last_blend_dst_alpha);
    if (last_enable_blend) glEnable(GL_BLEND); else glDisable(GL_BLEND);
    if (last_enable_cull_face) glEnable(GL_CULL_FACE); else glDisable(GL_CULL_FACE);
    if (last_enable_depth_test) glEnable(GL_DEPTH_TEST); else glDisable(GL_DEPTH_TEST);
    if (last_enable_scissor_test) glEnable(GL_SCISSOR_TEST); else glDisable(GL_SCISSOR_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, last_polygon_mode[0]);
    glViewport(last_viewport[0], last_viewport[1], (GLsizei)last_viewport[2], (GLsizei)last_viewport[3]);
    glScissor(last_scissor_box[0], last_scissor_box[1], (GLsizei)last_scissor_box[2], (GLsizei)last_scissor_box[3]);
}

static const char* imgui_get_clipboard_text(void* user_data) {
    return glfwGetClipboardString((GLFWwindow*)user_data);
}

static void imgui_set_clipboard_text(void* user_data, const char* text) {
    glfwSetClipboardString((GLFWwindow*)user_data, text);
}

static bool imgui_create_fonts_texture() {
    // Build texture atlas
    ImGuiIO& io = ImGui::GetIO();
    unsigned char* pixels;
    int width, height;
    io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);   // Load as RGBA 32-bits (75% of the memory is wasted, but default font is so small) because it is more likely to be compatible with user's existing shaders. If your ImTextureId represent a higher-level concept than just a GL texture id, consider calling GetTexDataAsAlpha8() instead to save on GPU memory.

    // Upload texture to graphics system
    GLint last_texture;
    glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture);
    glGenTextures(1, &g_font_texture);
    glBindTexture(GL_TEXTURE_2D, g_font_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

    // Store our identifier
    io.Fonts->TexID = (void *)(intptr_t)g_font_texture;

    // Restore state
    glBindTexture(GL_TEXTURE_2D, last_texture);

    return true;
}

bool imgui_create_device_objects() {
    // Backup GL state
    GLint last_texture, last_array_buffer, last_vertex_array;
    glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture);
    glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &last_array_buffer);
    glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &last_vertex_array);

    const GLchar *vertex_shader =
        "#version 150\n"
        "uniform mat4 ProjMtx;\n"
        "in vec2 Position;\n"
        "in vec2 UV;\n"
        "in vec4 Color;\n"
        "out vec2 Frag_UV;\n"
        "out vec4 Frag_Color;\n"
        "void main()\n"
        "{\n"
        "   Frag_UV = UV;\n"
        "   Frag_Color = Color;\n"
        "   gl_Position = ProjMtx * vec4(Position.xy,0,1);\n"
        "}\n";

    const GLchar* fragment_shader =
        "#version 150\n"
        "uniform sampler2D Texture;\n"
        "in vec2 Frag_UV;\n"
        "in vec4 Frag_Color;\n"
        "out vec4 Out_Color;\n"
        "void main()\n"
        "{\n"
        "   Out_Color = Frag_Color * texture( Texture, Frag_UV.st);\n"
        "}\n";

    g_shader_handle = glCreateProgram();
    g_vert_handle = glCreateShader(GL_VERTEX_SHADER);
    g_frag_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(g_vert_handle, 1, &vertex_shader, 0);
    glShaderSource(g_frag_handle, 1, &fragment_shader, 0);
    glCompileShader(g_vert_handle);
    glCompileShader(g_frag_handle);
    glAttachShader(g_shader_handle, g_vert_handle);
    glAttachShader(g_shader_handle, g_frag_handle);
    glLinkProgram(g_shader_handle);

    g_attrib_location_tex = glGetUniformLocation(g_shader_handle, "Texture");
    g_attrib_location_proj_mat = glGetUniformLocation(g_shader_handle, "ProjMtx");
    g_attrib_location_position = glGetAttribLocation(g_shader_handle, "Position");
    g_attrib_location_uv = glGetAttribLocation(g_shader_handle, "UV");
    g_attrib_location_color = glGetAttribLocation(g_shader_handle, "Color");

    glGenBuffers(1, &g_vbo_handle);
    glGenBuffers(1, &g_elements_handle);

    glGenVertexArrays(1, &g_vao_handle);
    glBindVertexArray(g_vao_handle);
    glBindBuffer(GL_ARRAY_BUFFER, g_vbo_handle);
    glEnableVertexAttribArray(g_attrib_location_position);
    glEnableVertexAttribArray(g_attrib_location_uv);
    glEnableVertexAttribArray(g_attrib_location_color);

    glVertexAttribPointer(g_attrib_location_position, 2, GL_FLOAT, GL_FALSE, sizeof(ImDrawVert), (GLvoid*)IM_OFFSETOF(ImDrawVert, pos));
    glVertexAttribPointer(g_attrib_location_uv, 2, GL_FLOAT, GL_FALSE, sizeof(ImDrawVert), (GLvoid*)IM_OFFSETOF(ImDrawVert, uv));
    glVertexAttribPointer(g_attrib_location_color, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(ImDrawVert), (GLvoid*)IM_OFFSETOF(ImDrawVert, col));

    imgui_create_fonts_texture();

    // Restore modified GL state
    glBindTexture(GL_TEXTURE_2D, last_texture);
    glBindBuffer(GL_ARRAY_BUFFER, last_array_buffer);
    glBindVertexArray(last_vertex_array);

    return true;
}

static void imgui_invalidate_device_objects() {
    if (g_vao_handle) glDeleteVertexArrays(1, &g_vao_handle);
    if (g_vbo_handle) glDeleteBuffers(1, &g_vbo_handle);
    if (g_elements_handle) glDeleteBuffers(1, &g_elements_handle);
    g_vao_handle = g_vbo_handle = g_elements_handle = 0;

    if (g_shader_handle && g_vert_handle) glDetachShader(g_shader_handle, g_vert_handle);
    if (g_vert_handle) glDeleteShader(g_vert_handle);
    g_vert_handle = 0;

    if (g_shader_handle && g_frag_handle) glDetachShader(g_shader_handle, g_frag_handle);
    if (g_frag_handle) glDeleteShader(g_frag_handle);
    g_frag_handle = 0;

    if (g_shader_handle) glDeleteProgram(g_shader_handle);
    g_shader_handle = 0;

    if (g_font_texture)
    {
        glDeleteTextures(1, &g_font_texture);
        ImGui::GetIO().Fonts->TexID = 0;
        g_font_texture = 0;
    }
}

bool imgui_init(GLFWwindow* window, bool install_callbacks) {
    g_window.glfw_window = window;

    ImGuiIO& io = ImGui::GetIO();
    io.KeyMap[ImGuiKey_Tab] = Key::KEY_TAB;                         // Keyboard mapping. ImGui will use those indices to peek into the io.KeyDown[] array.
    io.KeyMap[ImGuiKey_LeftArrow] = Key::KEY_LEFT;
    io.KeyMap[ImGuiKey_RightArrow] = Key::KEY_RIGHT;
    io.KeyMap[ImGuiKey_UpArrow] = Key::KEY_UP;
    io.KeyMap[ImGuiKey_DownArrow] = Key::KEY_DOWN;
    io.KeyMap[ImGuiKey_PageUp] = Key::KEY_PAGE_UP;
    io.KeyMap[ImGuiKey_PageDown] = Key::KEY_PAGE_DOWN;
    io.KeyMap[ImGuiKey_Home] = Key::KEY_HOME;
    io.KeyMap[ImGuiKey_End] = Key::KEY_END;
    io.KeyMap[ImGuiKey_Insert] = Key::KEY_INSERT;
    io.KeyMap[ImGuiKey_Delete] = Key::KEY_DELETE;
    io.KeyMap[ImGuiKey_Backspace] = Key::KEY_BACKSPACE;
    io.KeyMap[ImGuiKey_Enter] = Key::KEY_ENTER;
    io.KeyMap[ImGuiKey_Escape] = Key::KEY_ESCAPE;
    io.KeyMap[ImGuiKey_A] = Key::KEY_A;
    io.KeyMap[ImGuiKey_C] = Key::KEY_C;
    io.KeyMap[ImGuiKey_V] = Key::KEY_V;
    io.KeyMap[ImGuiKey_X] = Key::KEY_X;
    io.KeyMap[ImGuiKey_Y] = Key::KEY_Y;
    io.KeyMap[ImGuiKey_Z] = Key::KEY_Z;

    io.RenderDrawListsFn = imgui_render_draw_lists;       // Alternatively you can set this to NULL and call ImGui::GetDrawData() after ImGui::Render() to get the same ImDrawData pointer.
    io.SetClipboardTextFn = imgui_set_clipboard_text;
    io.GetClipboardTextFn = imgui_get_clipboard_text;
    io.ClipboardUserData = g_window.glfw_window;
#ifdef _WIN32
    io.ImeWindowHandle = glfwGetWin32Window(g_window.glfw_window);
#endif

    if (install_callbacks)
    {
        glfwSetMouseButtonCallback(window, mouse_button_callback);
        glfwSetScrollCallback(window, mouse_scroll_callback);
        glfwSetKeyCallback(window, key_callback);
        glfwSetCharCallback(window, char_callback);
    }

    return true;
}

static void imgui_shutdown()
{
    imgui_invalidate_device_objects();
}

static void imgui_new_frame()
{
    if (!g_font_texture)
        imgui_create_device_objects();

    ImGuiIO& io = ImGui::GetIO();

    // Setup display size (every frame to accommodate for window resizing)
    int w, h;
    int display_w, display_h;
    glfwGetWindowSize(g_window.glfw_window, &w, &h);
    glfwGetFramebufferSize(g_window.glfw_window, &display_w, &display_h);
    io.DisplaySize = ImVec2((float)w, (float)h);
    io.DisplayFramebufferScale = ImVec2(w > 0 ? ((float)display_w / w) : 0, h > 0 ? ((float)display_h / h) : 0);

    // Setup time step
    double current_time =  glfwGetTime();
    io.DeltaTime = g_time > 0.0 ? (float)(current_time - g_time) : (float)(1.0f/60.0f);
    g_time = current_time;

    // Setup inputs
    // (we already got mouse wheel, keyboard keys & characters from glfw callbacks polled in glfwPollEvents())
    if (glfwGetWindowAttrib(g_window.glfw_window, GLFW_FOCUSED)) {
        if (io.WantMoveMouse) {
            glfwSetCursorPos(g_window.glfw_window, (double)io.MousePos.x, (double)io.MousePos.y);   // Set mouse position if requested by io.WantMoveMouse flag (used when io.NavMovesTrue is enabled by user and using directional navigation)
        }
        else {
            double mouse_x, mouse_y;
            glfwGetCursorPos(g_window.glfw_window, &mouse_x, &mouse_y);
            io.MousePos = ImVec2((float)mouse_x, (float)mouse_y);
        }
    }
    else {
        io.MousePos = ImVec2(-FLT_MAX,-FLT_MAX);
    }

    for (int i = 0; i < 3; i++) {
        // If a mouse press event came, always pass it as "mouse held this frame", so we don't miss click-release events that are shorter than 1 frame.
        io.MouseDown[i] = g_input_state.mouse_hit[i] || glfwGetMouseButton(g_window.glfw_window, i) != 0;
        //g_input_state.mouse_hit[i] = false;
    }

    // Hide OS mouse cursor if ImGui is drawing it
    glfwSetInputMode(g_window.glfw_window, GLFW_CURSOR, io.MouseDrawCursor ? GLFW_CURSOR_HIDDEN : GLFW_CURSOR_NORMAL);

    // Start the frame. This call will update the io.WantCaptureMouse, io.WantCaptureKeyboard flag that you can use to dispatch inputs (or not) to your application.
    ImGui::NewFrame();
}


void initialize() {
    if (!glfwInit()) {
        // TODO Throw critical error
        error_callback(1, "Error while initializing. Terminating.");
        exit(1);
    }
    glfwSetErrorCallback(error_callback);
}

Window* create_window(int width, int height, const char* window_title) {
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    GLFWwindow* window = glfwCreateWindow(width, height, window_title, NULL, NULL);
    if (!window) {
        error_callback(2, "Could not create window. Terminating.");
        return nullptr;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    gl3wInit();

	if (glDebugMessageCallback) {
		glEnable(GL_DEBUG_OUTPUT);
		glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
		glDebugMessageCallback(gl_callback, nullptr);
		glDebugMessageControl(
			GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, NULL, true
		);
	}

	ImGui::CreateContext();
    imgui_init(window, true);

	double x, y;
	glfwGetCursorPos(window, &x, &y);
	g_input_state.prev_mouse_screen_coords = { x,y };
	g_input_state.mouse_screen_coords = { x,y };

	int w, h;
	glfwGetFramebufferSize(window, &w, &h);
	vec2 half_res = vec2(w, h) * 0.5f;
	g_input_state.mouse_ndc_coords = (vec2(x, h-y) - half_res) / half_res;
	g_input_state.prev_mouse_ndc_coords = g_input_state.mouse_ndc_coords;

    g_window.glfw_window = window;
    return &g_window; 
}

void destroy_window(Window* window) {
    glfwDestroyWindow(window->glfw_window);
}

void set_window_should_close(Window* window, bool value) {
	glfwSetWindowShouldClose(window->glfw_window, (int)value);
}

void shutdown() {
    imgui_shutdown();
	ImGui::DestroyContext();
    glfwTerminate();
}

void update() {
    // Reset hit states
    memset(g_input_state.key_hit, 0, MAX_KEYS);
    memset(g_input_state.key_release, 0, MAX_KEYS);
    memset(g_input_state.mouse_hit, 0, MAX_MOUSE_BUTTONS);
    memset(g_input_state.mouse_release, 0, MAX_MOUSE_BUTTONS);
	g_input_state.mouse_scroll = { 0,0 };

    glfwPollEvents();
    double x, y;
    glfwGetCursorPos(g_window.glfw_window, &x, &y);
	int w, h;
	glfwGetFramebufferSize(g_window.glfw_window, &w, &h);

	vec2 new_coord{ x,y };
	g_input_state.mouse_velocity = new_coord - g_input_state.mouse_screen_coords;
	g_input_state.prev_mouse_screen_coords = g_input_state.mouse_screen_coords;
	g_input_state.mouse_screen_coords = new_coord;

	vec2 half_res = vec2(w, h) * 0.5f;
	g_input_state.prev_mouse_ndc_coords = g_input_state.mouse_ndc_coords;
	g_input_state.mouse_ndc_coords = (vec2(x, h-y) - half_res) / half_res;

    imgui_new_frame();
}

bool window_in_focus(Window* window) {
    return glfwGetWindowAttrib(window->glfw_window, GLFW_FOCUSED);
}

bool window_should_close(Window* window) {
    return glfwWindowShouldClose(window->glfw_window);
}

void get_framebuffer_size(Window* window, int* width, int* height) {
    glfwGetFramebufferSize(window->glfw_window, width, height);
}

InputState* get_input_state() {
    return &g_input_state;
}

double get_delta_time() {
	return ImGui::GetIO().DeltaTime;
}

void swap_buffers(Window* window) {
    glfwSwapBuffers(window->glfw_window);
}

void set_vsync(bool value) {
	glfwSwapInterval((int)value);
}

}