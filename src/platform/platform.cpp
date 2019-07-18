#include "platform.h"
#include <core/gl.h>
#include <core/log.h>
#include <core/string_utils.h>
#include <GLFW/glfw3.h>
#ifdef OS_WINDOWS
//#undef APIENTRY
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN 1
#endif
#define GLFW_EXPOSE_NATIVE_WIN32
#define GLFW_EXPOSE_NATIVE_WGL
#include <GLFW/glfw3native.h>
#endif
#include <imgui.h>
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
//#include "cousine_font.inl"
#include "droid_sans.inl"
#include <stdio.h>
#include <stdlib.h>

#include <nfd.h>

#define ARRAY_SIZE(x) (sizeof(x) / sizeof(x[0]))

namespace platform {

// Data
static struct {
    Context internal_ctx{};

    struct {
        Path cwd{};
    } file_system;
} data;


static void error_callback(int error, const char* description) { LOG_ERROR("%d: %s\n", error, description); }

static void APIENTRY gl_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam) {
    (void)source;
    (void)type;
    (void)id;
    (void)severity;
    (void)length;
    (void)userParam;

    if (severity > GL_DEBUG_SEVERITY_LOW) {
        LOG_WARNING(message);
    }
    if (severity == GL_DEBUG_SEVERITY_HIGH) {
        LOG_ERROR("A SEVERE GL ERROR HAS OCCURED: %s", message);
        abort();
    }
}

static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (action == GLFW_PRESS) {
        data.internal_ctx.input.mouse.down[button] = true;
        data.internal_ctx.input.mouse.hit[button] = true;
    } else if (action == GLFW_RELEASE) {
        data.internal_ctx.input.mouse.down[button] = false;
        data.internal_ctx.input.mouse.release[button] = true;
    }

    ImGui_ImplGlfw_MouseButtonCallback(window, button, action, mods);
}

static void mouse_scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    data.internal_ctx.input.mouse.scroll_delta = (float)yoffset;
    ImGui_ImplGlfw_ScrollCallback(window, xoffset, yoffset);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS) {
        data.internal_ctx.input.key.down[key] = true;
        data.internal_ctx.input.key.hit[key] = true;
    }
    if (action == GLFW_RELEASE) {
        data.internal_ctx.input.key.down[key] = false;
    }

    ImGui_ImplGlfw_KeyCallback(window, key, scancode, action, mods);
}

static void char_callback(GLFWwindow* window, unsigned int c) { ImGui_ImplGlfw_CharCallback(window, c); }

bool initialize(Context* ctx, int32 width, int32 height, const char* title) {
    if (!glfwInit()) {
        // TODO Throw critical error
        error_callback(1, "Error while initializing glfw.");
        return false;
    }
    glfwSetErrorCallback(error_callback);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef OS_MAC_OSX
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    GLFWwindow* window = glfwCreateWindow(width, height, title, NULL, NULL);
    if (!window) {
        error_callback(2, "Could not create glfw window.");
        return false;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    if (gl3wInit() != GL3W_OK) {
        error_callback(3, "Could not load gl functions.");
        return false;
    }

    if (glDebugMessageCallback) {
        glEnable(GL_DEBUG_OUTPUT);
        glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
        glDebugMessageCallback(gl_callback, NULL);
        glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, NULL, true);
    }

    glfwGetVersion(&data.internal_ctx.gl_info.version.major, &data.internal_ctx.gl_info.version.minor, &data.internal_ctx.gl_info.version.revision);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard Controls
    // io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls
    // io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;  // Enable Docking
    // io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;  // Enable Multi-Viewport / Platform Windows
    // io.ConfigFlags |= ImGuiConfigFlags_ViewportsNoTaskBarIcons;
    // io.ConfigFlags |= ImGuiConfigFlags_ViewportsNoMerge;
    io.ConfigFlags |= ImGuiConfigFlags_DpiEnableScaleFonts;
    //io.ConfigDockingWithShift = true;
	io.ConfigWindowsMoveFromTitleBarOnly = true;

    // default range is 0x0020 - 0x00FF.
    // In our case we BAKE ALL!
    // @TODO: reduce this to only the used characters... for some release...
    static const ImWchar ranges[] = {0x0020, 0xFFFF, 0};
    ImFontConfig config;
    config.OversampleV = 4;
    config.OversampleH = 4;
    ImGui::GetIO().Fonts->AddFontFromMemoryCompressedTTF((void*)DroidSans_compressed_data, DroidSans_compressed_size, 14.f, &config, ranges);

    ImGui_ImplGlfw_InitForOpenGL(window, false);
    ImGui_ImplOpenGL3_Init("#version 150");

    data.internal_ctx.window.ptr = window;
    data.internal_ctx.window.title = title;
    data.internal_ctx.window.width = width;
    data.internal_ctx.window.height = height;
    data.internal_ctx.window.vsync = true;

    int w, h;
    glfwGetFramebufferSize(window, &w, &h);
    data.internal_ctx.framebuffer.width = w;
    data.internal_ctx.framebuffer.height = h;

    double x, y;
    glfwGetCursorPos(window, &x, &y);
    Coordinate win_coord = {(float)x, (float)y};
    data.internal_ctx.input.mouse.win_coord = win_coord;

    const float half_res_x = width * 0.5f;
    const float half_res_y = height * 0.5f;
    Coordinate ndc_coord = {(win_coord.x - half_res_x) / half_res_x, ((height - win_coord.y) - half_res_y) / half_res_y};
    data.internal_ctx.input.mouse.ndc_delta = ndc_coord - data.internal_ctx.input.mouse.ndc_coord;
    data.internal_ctx.input.mouse.ndc_coord = ndc_coord;

    data.internal_ctx.input.mouse.moving = false;

    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, mouse_scroll_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetCharCallback(window, char_callback);

    memcpy(ctx, &data.internal_ctx, sizeof(Context));

    return true;
}

void shutdown(Context* ctx) {
    glfwDestroyWindow((GLFWwindow*)data.internal_ctx.window.ptr);
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwTerminate();

    ctx->window.ptr = nullptr;
    ctx->window.title = "";
    ctx->window.width = 0;
    ctx->window.height = 0;
}

void update(Context* ctx) {
    // Reset hit states
    memset(data.internal_ctx.input.key.hit, 0, MAX_KEYS);
    memset(data.internal_ctx.input.key.release, 0, MAX_KEYS);
    memset(data.internal_ctx.input.mouse.hit, 0, MAX_MOUSE_BUTTONS);
    memset(data.internal_ctx.input.mouse.release, 0, MAX_MOUSE_BUTTONS);
    data.internal_ctx.input.mouse.scroll_delta = 0;

    glfwPollEvents();

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    if (ctx->window.width != data.internal_ctx.window.width || ctx->window.height != data.internal_ctx.window.height) {
        glfwSetWindowSize((GLFWwindow*)ctx->window.ptr, ctx->window.width, ctx->window.height);
        data.internal_ctx.window.width = ctx->window.width;
        data.internal_ctx.window.height = ctx->window.height;
    }
    int w, h;
    glfwGetFramebufferSize((GLFWwindow*)data.internal_ctx.window.ptr, &w, &h);
    data.internal_ctx.framebuffer.width = w;
    data.internal_ctx.framebuffer.height = h;

    glfwGetWindowSize((GLFWwindow*)data.internal_ctx.window.ptr, &w, &h);
    data.internal_ctx.window.width = w;
    data.internal_ctx.window.height = h;

    if (ctx->window.vsync != data.internal_ctx.window.vsync) {
        data.internal_ctx.window.vsync = ctx->window.vsync;
        glfwSwapInterval((int)ctx->window.vsync);
    }

    data.internal_ctx.window.should_close = (bool)glfwWindowShouldClose((GLFWwindow*)data.internal_ctx.window.ptr);

    double x, y;
    glfwGetCursorPos((GLFWwindow*)data.internal_ctx.window.ptr, &x, &y);
    Coordinate win_coord{(float)x, (float)y};

    // If user has modified value, set the mouse pointer to that value
    if (ctx->input.mouse.win_coord != data.internal_ctx.input.mouse.win_coord) {
        win_coord = ctx->input.mouse.win_coord;
    }

    data.internal_ctx.input.mouse.win_delta = win_coord - data.internal_ctx.input.mouse.win_coord;
    data.internal_ctx.input.mouse.win_coord = win_coord;

    constexpr Coordinate zero_coord{0, 0};
    data.internal_ctx.input.mouse.moving = data.internal_ctx.input.mouse.win_delta != zero_coord;

    const float half_res_x = w * 0.5f;
    const float half_res_y = h * 0.5f;
    data.internal_ctx.input.mouse.ndc_coord = {(win_coord.x - half_res_x) / half_res_x, ((h - win_coord.y) - half_res_y) / half_res_y};

    // Any mouse key hit
    bool any_mouse_key_hit = false;
    for (int i = 0; i < MAX_MOUSE_BUTTONS; i++) {
        if (data.internal_ctx.input.mouse.hit[i]) {
            any_mouse_key_hit = true;
            break;
        }
    }

    static Coordinate mouse_click_coord = win_coord;
    if (any_mouse_key_hit) {
        mouse_click_coord = win_coord;
    }

    for (int i = 0; i < MAX_MOUSE_BUTTONS; i++) {
        data.internal_ctx.input.mouse.clicked[i] = (mouse_click_coord == win_coord) && data.internal_ctx.input.mouse.release[i];
    }

    double t = glfwGetTime();
    data.internal_ctx.timing.delta_s = (float)(t - data.internal_ctx.timing.total_s);
    data.internal_ctx.timing.total_s = t;

    memcpy(ctx, &data.internal_ctx, sizeof(Context));
}

void render_imgui(Context* ctx) {
    (void)ctx;
    GLFWwindow* window = (GLFWwindow*)data.internal_ctx.window.ptr;

    ImGui::Render();
    glfwMakeContextCurrent(window);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    // Update and Render additional Platform Windows
    if (ImGui::GetIO().ConfigFlags & ImGuiConfigFlags_ViewportsEnable) {
        ImGui::UpdatePlatformWindows();
        ImGui::RenderPlatformWindowsDefault();
    }

    glfwMakeContextCurrent(window);
}

void swap_buffers(Context* ctx) { glfwSwapBuffers((GLFWwindow*)ctx->window.ptr); }

FileDialogResult file_dialog(FileDialogFlags flags, CString default_path, CString filter) {
    // Create zero terminated strings
    Path path_buf = default_path;
    convert_backslashes(path_buf);
    StringBuffer<256> filter_buf = filter;

    nfdchar_t* out_path = NULL;
    nfdresult_t result = NFD_ERROR;

    if (flags & FileDialogFlags_Open) {
        result = NFD_OpenDialog(filter_buf.cstr(), path_buf.cstr(), &out_path);
    } else if (flags & FileDialogFlags_Save) {
        result = NFD_SaveDialog(filter_buf.cstr(), path_buf.cstr(), &out_path);
    }
    defer {
        if (out_path) free(out_path);
    };

    if (result == NFD_OKAY) {
        Path res_path = out_path;
        convert_backslashes(res_path);
        return {res_path, FileDialogResult::FILE_OK};
    } else if (result == NFD_CANCEL) {
        return {{}, FileDialogResult::FILE_CANCEL};
    }
    LOG_ERROR("%s\n", NFD_GetError());
    return {{}, FileDialogResult::FILE_CANCEL};
}

#ifdef OS_WINDOWS
#include "platform_win32.inl"
#elif defined OS_MAC_OSX
#include "platform_osx.inl"
#elif defined OS_LINUX
#include "platform_linux.inl"
#endif

}  // namespace platform
