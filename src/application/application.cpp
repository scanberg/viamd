#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "application.h"

#include <core/md_log.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_platform.h>
#include <core/md_str.h>

#include "gfx/gl.h"
#include <GLFW/glfw3.h>
#include <nfd.h>

#include <imgui.h>
#include <implot.h>
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

// Compressed fonts
#include "dejavu_sans_mono.inl"
#include "fa_solid.inl"
#include "IconsFontAwesome6.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

namespace application {

// Data
static struct {
    Context internal_ctx{};
} data;

static void error_callback(int error, const char* description) { MD_LOG_ERROR("%d: %s\n", error, description); }

static void APIENTRY gl_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message,
                                 const void* userParam) {
    (void)source;
    (void)type;
    (void)id;
    (void)severity;
    (void)length;
    (void)userParam;

    if (severity == GL_DEBUG_SEVERITY_HIGH) {
        MD_LOG_ERROR("A SEVERE GL ERROR HAS OCCURED: %s", message);
        ASSERT(false);
    } else {
        //MD_LOG_INFO("%s", message);
    }
}

bool initialize(Context* ctx, int width, int height, const char* title) {
    if (!glfwInit()) {
        // TODO Throw critical error
        MD_LOG_ERROR("Error while initializing glfw.");
        return false;
    }
    glfwSetErrorCallback(error_callback);

    if (width == 0 && height == 0) {
        int count;
        GLFWmonitor** monitors = glfwGetMonitors(&count);
        if (count > 0) {
            int pos_x, pos_y, dim_x, dim_y;
            glfwGetMonitorWorkarea(monitors[0], &pos_x, &pos_y, &dim_x, &dim_y);
            width  = (int)(dim_x * 0.9);
            height = (int)(dim_y * 0.9);
        }
    }

#if MD_PLATFORM_OSX
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    GLFWwindow* window = glfwCreateWindow((int)width, (int)height, title, NULL, NULL);
    if (!window) {
        MD_LOG_ERROR("Could not create glfw window.");
        return false;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    if (gl3wInit() != GL3W_OK) {
        MD_LOG_ERROR("Could not load gl functions.");
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
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard Controls
    // io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;  // Enable Docking
#if !MD_PLATFORM_LINUX
    //io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;  // Enable Multi-Viewport / Platform Windows
#endif
    // io.ConfigFlags |= ImGuiConfigFlags_ViewportsNoTaskBarIcons;
    // io.ConfigFlags |= ImGuiConfigFlags_ViewportsNoMerge;
    // io.ConfigFlags |= ImGuiConfigFlags_DpiEnableScaleFonts;
    // io.ConfigDockingWithShift = true;
    io.ConfigWindowsMoveFromTitleBarOnly = true;

    float xscale, yscale;
    glfwGetWindowContentScale(window, &xscale, &yscale);
    const float dpi_scale = (xscale + yscale) * 0.5f;

    // default range is 0x0020 - 0x00FF.
    // Added some greek letters
    const ImWchar ranges_characters[] = {0x0020, 0x00FF, 0x03C6, 0x03C8, 0x2074, 0x207b, 0};
    const ImWchar ranges_icons[] = {ICON_MIN_FA, ICON_MAX_FA, 0};
    const float font_size[] = {16, 18, 20, 24, 30};
    const char* font_names[] = {"DejaVu Sans Mono 16", "DejaVu Sans Mono 18", "DejaVu Sans Mono 20", "DejaVu Sans Mono 24", "DejaVu Sans Mono 30"};

    STATIC_ASSERT(ARRAY_SIZE(font_size) == ARRAY_SIZE(font_names), "font_size and font_names must have the same size");
    
    for (size_t i = 0; i < ARRAY_SIZE(font_size); ++i) {
        ImFontConfig config;
        snprintf(config.Name, sizeof(config.Name), "%s", font_names[i]);
        config.OversampleV = 2;
        config.OversampleH = 3;
        config.PixelSnapH = true;

        const float size = font_size[i];

        // CHARACTERS
        config.RasterizerMultiply = 1.f;
        ImGui::GetIO().Fonts->AddFontFromMemoryCompressedTTF((void*)dejavu_sans_mono_compressed_data, dejavu_sans_mono_compressed_size, size, &config, ranges_characters);

        // ICONS
        const float scl = 0.75f; // 0.875f;   // We scale this a bit to better fit within buttons and such.
        config.RasterizerMultiply = 0.9f;
        config.MergeMode = true;
        ImGui::GetIO().Fonts->AddFontFromMemoryCompressedTTF((void*)fa_solid_compressed_data, fa_solid_compressed_size, size * scl, &config, ranges_icons);
    }

    io.Fonts->Build();
    io.FontDefault = io.Fonts->Fonts[1]; // Set default to 18px

    ImGui_ImplGlfw_InitForOpenGL(window, false);
    ImGui_ImplOpenGL3_Init("#version 150");

    data.internal_ctx.window.ptr = window;
    data.internal_ctx.window.title = title;
    data.internal_ctx.window.width = (int)width;
    data.internal_ctx.window.height = (int)height;
    data.internal_ctx.window.vsync = true;

    int w, h;
    glfwGetFramebufferSize(window, &w, &h);
    data.internal_ctx.framebuffer.width = w;
    data.internal_ctx.framebuffer.height = h;

    glfwSetMouseButtonCallback(window, ImGui_ImplGlfw_MouseButtonCallback);
    glfwSetScrollCallback(window, ImGui_ImplGlfw_ScrollCallback);
    glfwSetKeyCallback(window, ImGui_ImplGlfw_KeyCallback);
    glfwSetCharCallback(window, ImGui_ImplGlfw_CharCallback);

    glfwSetWindowUserPointer(window, &data.internal_ctx);

    GLFWdropfun drop_cb = [](GLFWwindow* window, int num_files, const char** paths) {
        Context* ctx = (Context*)glfwGetWindowUserPointer(window);
        ASSERT(ctx);
        
        for (int i = 0; i < num_files; ++i) {
            MD_LOG_DEBUG("User dropped file: '%s'", paths[i]);
        }
        
        if (ctx->file_drop.callback) {
            ctx->file_drop.callback(num_files, paths, ctx->file_drop.user_data);
        }
    };
    glfwSetDropCallback(window, drop_cb);

    memcpy(ctx, &data.internal_ctx, sizeof(Context));

    return true;
}

void shutdown(Context* ctx) {
    glfwDestroyWindow((GLFWwindow*)data.internal_ctx.window.ptr);
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    glfwTerminate();

    ctx->window.ptr = nullptr;
    ctx->window.title = "";
    ctx->window.width = 0;
    ctx->window.height = 0;
}

void update(Context* ctx) {
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

    if (ctx->file_drop.callback != data.internal_ctx.file_drop.callback) {
        data.internal_ctx.file_drop.callback = ctx->file_drop.callback;
    }

    if (ctx->file_drop.user_data != data.internal_ctx.file_drop.user_data) {
        data.internal_ctx.file_drop.user_data = ctx->file_drop.user_data;
    }

    data.internal_ctx.window.should_close = (bool)glfwWindowShouldClose((GLFWwindow*)data.internal_ctx.window.ptr);

    double t = glfwGetTime();
    data.internal_ctx.timing.delta_s = (t - data.internal_ctx.timing.total_s);
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

bool file_dialog(char* str_buf, int str_cap, FileDialogFlags flags, const char* filter) {    
    nfdchar_t* out_path = NULL;
    defer { if (out_path) free(out_path); };

    nfdresult_t result = NFD_ERROR;

    const char* default_path = 0;

    if (flags & FileDialog_Open) {
        result = NFD_OpenDialog(filter, default_path, &out_path);
    } else if (flags & FileDialog_Save) {
        result = NFD_SaveDialog(filter, default_path, &out_path);
    }

    if (result == NFD_OKAY) {
        char ext[16] = {0};
        if (flags & FileDialog_Save) {
            str_t path = str_from_cstr(out_path);
            str_t pext  = extract_ext(path);
            if (str_empty(pext) && filter) {
                // get ext from supplied filter (first match)
                pext.ptr = filter;
                const char* delim = strchr(filter, ',');
                if (delim) {
                    pext.len = delim - filter;
                } else {
                    pext.len = strlen(filter);
                }
                snprintf(ext, sizeof(ext), ".%.*s", (int)pext.len, pext.ptr);
            }
        }
        snprintf(str_buf, str_cap, "%s%s", out_path, ext);
        convert_backslashes(str_buf, str_cap);
        return true;
    } else if (result == NFD_ERROR) {
        MD_LOG_ERROR("%s\n", NFD_GetError());
    }
    /* fallthrough for NFD_CANCEL */
    return false;
}

}  // namespace application
