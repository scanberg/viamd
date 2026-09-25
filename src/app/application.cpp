#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <app/application.h>

#include <core/md_log.h>
#include <core/md_common.h>
#include <core/md_platform.h>
#include <core/md_str.h>
#include <core/md_allocator.h>

#include <gfx/gl.h>
#include <GLFW/glfw3.h>

#if MD_PLATFORM_WINDOWS
#define GLFW_EXPOSE_NATIVE_WIN32
#include <GLFW/glfw3native.h>
#elif MD_PLATFORM_LINUX
#define GLFW_EXPOSE_NATIVE_X11
#include <GLFW/glfw3native.h>
#endif
#include <nfd.h>
#include <nfd_glfw3.h>

#include <imgui.h>
#include <implot.h>

#include <app/imgui_impl_glfw.h>
#include <app/imgui_impl_opengl3.h>

// Compressed fonts
#include <app/dejavu_sans_mono.inl>
#include <app/fa_solid.inl>
#include <app/IconsFontAwesome6.h>

#include <stdio.h> // snprintf
#include <stdint.h> // uint64_t
#include <stdlib.h> // free

namespace application {

// Data
static struct {
    Context internal_ctx{};
} data;

static void error_callback(int error, const char* description) { MD_LOG_ERROR("%d: %s\n", error, description); }

static const char* gl_debug_source_str(GLenum v) {
    switch (v) {
    case GL_DEBUG_SOURCE_API:             return "API";
    case GL_DEBUG_SOURCE_WINDOW_SYSTEM:   return "WINDOW_SYSTEM";
    case GL_DEBUG_SOURCE_SHADER_COMPILER: return "SHADER_COMPILER";
    case GL_DEBUG_SOURCE_THIRD_PARTY:     return "THIRD_PARTY";
    case GL_DEBUG_SOURCE_APPLICATION:     return "APPLICATION";
    default:                              return "OTHER";
    }
}

static const char* gl_debug_type_str(GLenum v) {
    switch (v) {
    case GL_DEBUG_TYPE_ERROR:               return "ERROR";
    case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR: return "DEPRECATED";
    case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR:  return "UNDEFINED";
    case GL_DEBUG_TYPE_PORTABILITY:         return "PORTABILITY";
    case GL_DEBUG_TYPE_PERFORMANCE:         return "PERFORMANCE";
    case GL_DEBUG_TYPE_MARKER:              return "MARKER";
    case GL_DEBUG_TYPE_PUSH_GROUP:          return "PUSH_GROUP";
    case GL_DEBUG_TYPE_POP_GROUP:           return "POP_GROUP";
    default:                                return "OTHER";
    }
}

static const char* gl_debug_severity_str(GLenum v) {
    switch (v) {
    case GL_DEBUG_SEVERITY_HIGH:         return "HIGH";
    case GL_DEBUG_SEVERITY_MEDIUM:       return "MEDIUM";
    case GL_DEBUG_SEVERITY_LOW:          return "LOW";
    case GL_DEBUG_SEVERITY_NOTIFICATION: return "NOTIFICATION";
    default:                             return "UNKNOWN";
    }
}

// KHR_debug leaves the behaviour of calling any GL or window system function from
// inside the debug callback undefined, so everything below works from the message
// itself plus state we track ourselves. Context identification is logged once at
// startup instead (see initialize()).
#define GL_DEBUG_GROUP_STACK_SIZE 8
#define GL_DEBUG_GROUP_LABEL_LEN  48

static struct {
    char     group[GL_DEBUG_GROUP_STACK_SIZE][GL_DEBUG_GROUP_LABEL_LEN];
    int      depth;
    GLuint   last_id;
    uint64_t repeat;
} gl_debug_state = {};

// Renders the active PUSH_GPU_SECTION labels as "G-buffer > Postprocessing", which is
// what tells us where in the frame the message came from.
static void gl_debug_group_path(char* buf, size_t cap) {
    size_t off = 0;
    buf[0] = '\0';
    int depth = MIN(gl_debug_state.depth, GL_DEBUG_GROUP_STACK_SIZE);
    for (int i = 0; i < depth; ++i) {
        int n = snprintf(buf + off, cap - off, "%s%s", (i > 0) ? " > " : "", gl_debug_state.group[i]);
        if (n < 0 || (size_t)n >= cap - off) break;
        off += (size_t)n;
    }
}

static bool is_power_of_ten(uint64_t v) {
    uint64_t p = 1;
    while (p < v) {
        uint64_t next = p * 10;
        if (next < p) return false;  // overflow
        p = next;
    }
    return p == v;
}

static void APIENTRY gl_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message,
                                 const void* userParam) {
    (void)userParam;

    // glPushDebugGroup / glPopDebugGroup are reported through this same callback, so
    // mirror the stack here. Depth is tracked past the array bound so pushes and pops
    // stay balanced even if we stop recording labels.
    if (type == GL_DEBUG_TYPE_PUSH_GROUP) {
        if (gl_debug_state.depth >= 0 && gl_debug_state.depth < GL_DEBUG_GROUP_STACK_SIZE) {
            char* dst = gl_debug_state.group[gl_debug_state.depth];
            if (length > 0) {
                snprintf(dst, GL_DEBUG_GROUP_LABEL_LEN, "%.*s", (int)length, message);
            } else {
                snprintf(dst, GL_DEBUG_GROUP_LABEL_LEN, "%s", message);
            }
        }
        gl_debug_state.depth += 1;
        return;
    }
    if (type == GL_DEBUG_TYPE_POP_GROUP) {
        if (gl_debug_state.depth > 0) gl_debug_state.depth -= 1;
        return;
    }

    // Drivers disagree on severity -- the same GL_INVALID_ENUM is reported HIGH by one
    // vendor and LOW by another -- so classify on type and log everything that is not
    // notification spam. A bug report from a machine we do not have is close to
    // useless without the source/type/id triplet and the frame location.
    if (severity == GL_DEBUG_SEVERITY_NOTIFICATION) {
        return;
    }

    // A bad call inside the render loop repeats every frame. Report the 1st, 10th,
    // 100th ... occurrence of a given id and drop the rest, so a persistent error
    // stays visible without burying everything else.
    if (id == gl_debug_state.last_id) {
        gl_debug_state.repeat += 1;
        if (!is_power_of_ten(gl_debug_state.repeat)) return;
    } else {
        gl_debug_state.last_id = id;
        gl_debug_state.repeat  = 1;
    }

    char where[256];
    gl_debug_group_path(where, sizeof(where));

    char count[32];
    count[0] = '\0';
    if (gl_debug_state.repeat > 1) {
        snprintf(count, sizeof(count), " (repeated %llux)", (unsigned long long)gl_debug_state.repeat);
    }

    const char* src_str  = gl_debug_source_str(source);
    const char* type_str = gl_debug_type_str(type);
    const char* sev_str  = gl_debug_severity_str(severity);
    const char* msg      = message ? message : "(no message)";

    if (type == GL_DEBUG_TYPE_ERROR || severity == GL_DEBUG_SEVERITY_HIGH) {
        MD_LOG_DEBUG("GL %s [%s/%s] (id %u) during '%s'%s: %s",
                     type_str, src_str, sev_str, (unsigned int)id,
                     where[0] ? where : "no active debug group", count, msg);
    } else {
        MD_LOG_DEBUG("GL %s [%s/%s] (id %u) during '%s'%s: %s",
                    type_str, src_str, sev_str, (unsigned int)id,
                    where[0] ? where : "no active debug group", count, msg);
    }

    if (severity == GL_DEBUG_SEVERITY_HIGH) {
        ASSERT(false);
    }
}

bool initialize(Context* ctx, size_t width, size_t height, str_t title) {
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
            width  = (size_t)(dim_x * 0.9);
            height = (size_t)(dim_y * 0.8);
        }
    }

#if MD_PLATFORM_OSX
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    // Zero terminated
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    str_t ztitle = str_copy(title, md_temp_allocator(temp));
    GLFWwindow* window = glfwCreateWindow((int)width, (int)height, ztitle.ptr, NULL, NULL);
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

    // Log the context identity unconditionally: every later GL message in the log is
    // only actionable if we know which driver produced it, and this is the first thing
    // to ask for in a bug report.
    {
        auto gl_str = [](GLenum name) -> const char* {
            const char* s = (const char*)glGetString(name);
            return s ? s : "(unavailable)";
        };
        MD_LOG_INFO(
            "GL_VENDOR:   %s\n"
            "GL_RENDERER: %s\n"
            "GL_VERSION:  %s\n"
            "GL_SHADING_LANGUAGE_VERSION: %s",
            gl_str(GL_VENDOR),
            gl_str(GL_RENDERER),
            gl_str(GL_VERSION),
            gl_str(GL_SHADING_LANGUAGE_VERSION));

        if (glDispatchCompute) {
            GLint shared_mem = 0;
            GLint max_draw_buffers = 0;
            glGetIntegerv(GL_MAX_COMPUTE_SHARED_MEMORY_SIZE, &shared_mem);
            glGetIntegerv(GL_MAX_DRAW_BUFFERS, &max_draw_buffers);
            MD_LOG_INFO("GL_MAX_COMPUTE_SHARED_MEMORY_SIZE: %i", shared_mem);
            MD_LOG_INFO("GL_MAX_DRAW_BUFFERS: %i", max_draw_buffers);
        }
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
    
#if VIAMD_IMGUI_ENABLE_VIEWPORTS
    io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;  // Enable Multi-Viewport / Platform Windows
#endif
    // io.ConfigFlags |= ImGuiConfigFlags_ViewportsNoTaskBarIcons;
    // io.ConfigFlags |= ImGuiConfigFlags_ViewportsNoMerge;
    // io.ConfigDockingWithShift = true;

    io.ConfigWindowsMoveFromTitleBarOnly = true;
    //io.ConfigDpiScaleFonts = true;

    //float xscale, yscale;
    //glfwGetWindowContentScale(window, &xscale, &yscale);
    //const float dpi_scale = (xscale + yscale) * 0.5f;

    // Default range:               0x0020 - 0x00FF.
    // Greek and Coptik:            0x0370 - 0x03FF
    // General punctuation (dashes, quotes, ...): 0x2010 - 0x2027
    // Superscripts and Subscripts: 0x2070 - 0x209F
    const ImWchar ranges_characters[] = {0x0020, 0x00FF, 0x0370, 0x03FF, 0x2010, 0x2027, 0x2070, 0x209F, 0};
    const ImWchar ranges_icons[] = {ICON_MIN_FA, ICON_MAX_FA, 0};
    const float font_size = 18.0f;
    const char* font_name = "Dejavu Sans Mono";
    const float icons_size = font_size * 0.75f; // Scale to better fit buttons
    
    ImFontConfig config;
    snprintf(config.Name, sizeof(config.Name), "%s", font_name);
    config.OversampleV = 2;
    config.OversampleH = 3;
    config.PixelSnapH = true;

    // CHARACTERS
    ImGui::GetIO().Fonts->AddFontFromMemoryCompressedTTF((void*)dejavu_sans_mono_compressed_data, dejavu_sans_mono_compressed_size, font_size, &config, ranges_characters);

    // ICONS
    config.MergeMode = true;
    ImGui::GetIO().Fonts->AddFontFromMemoryCompressedTTF((void*)fa_solid_compressed_data, fa_solid_compressed_size, icons_size, &config, ranges_icons);

    io.FontDefault = io.Fonts->Fonts[0]; // Set default to 18px

    if (!ImGui_ImplGlfw_InitForOpenGL(window, false) ||
        !ImGui_ImplOpenGL3_Init("#version 150"))
    {
        MD_LOG_ERROR("Failed to initialize ImGui OpenGL");
        return false;
    }

    data.internal_ctx.window.ptr = window;
    data.internal_ctx.window.title = title;
    data.internal_ctx.window.width = (int)width;
    data.internal_ctx.window.height = (int)height;
    data.internal_ctx.window.vsync = true;

    
    int w, h;
    glfwGetFramebufferSize(window, &w, &h);
    data.internal_ctx.framebuffer.width = w;
    data.internal_ctx.framebuffer.height = h;
    glfwGetWindowContentScale(window, &data.internal_ctx.window.scale_factor, NULL);

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

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        str_t* str_paths = md_temp_alloc_array(temp, str_t, num_files);
        for (int i = 0; i < num_files; ++i) {
            str_paths[i] = {paths[i], strlen(paths[i])};
        }
        
        if (ctx->file_drop.callback) {
            ctx->file_drop.callback((size_t)num_files, str_paths, ctx->file_drop.user_data);
        }
    };
    glfwSetDropCallback(window, drop_cb);
    
#if MD_PLATFORM_WINDOWS
    HWND hwnd = glfwGetWin32Window(window);
    HINSTANCE hinst = GetModuleHandle(NULL);

    HICON hIcon = LoadIcon(hinst, "VIAMD_ICON");
    SendMessage(hwnd, WM_SETICON, ICON_SMALL, (LPARAM)hIcon);
    SendMessage(hwnd, WM_SETICON, ICON_BIG, (LPARAM)hIcon);
#endif

    if (NFD_Init() != NFD_OKAY) {
        MD_LOG_ERROR("Failed to initialize NativeFileDialog: %s", NFD_GetError());
    }

    MEMCPY(ctx, &data.internal_ctx, sizeof(Context));

    return true;
}

void shutdown(Context* ctx) {
    NFD_Quit();

    glfwDestroyWindow((GLFWwindow*)data.internal_ctx.window.ptr);
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    glfwTerminate();

    MEMSET(ctx, 0, sizeof(Context));
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

    MEMCPY(ctx, &data.internal_ctx, sizeof(Context));
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

bool file_dialog(char* str_buf, size_t str_cap, FileDialogFlag flags, str_t filter) {    
    md_temp_scope_t temp = md_temp_begin();
    nfdu8char_t* out_path = NULL;
    defer {
        md_temp_end(temp);
        if (out_path) NFD_FreePathU8(out_path);
    };

    nfdresult_t result = NFD_ERROR;

    // Zero terminated variant
    str_t zfilt = str_copy(filter, md_temp_allocator(temp));

    // NFD(e)'s filter spec is a comma-separated list of extensions without the leading dot
    // (e.g. "png,jpg,jpeg"), which is exactly the format callers already pass as `filter` here.
    const nfdu8filteritem_t filter_item = { "Files", zfilt.ptr };
    const nfdfiltersize_t filter_count = zfilt.len > 0 ? 1 : 0;
    const nfdu8filteritem_t* filter_list = filter_count ? &filter_item : NULL;

    // Associate the dialog with our window so the window manager (and the portal backend, if
    // enabled) parents and positions it correctly instead of it popping up as an unrelated
    // top-level window. Degrades to "no parent" (matching prior behavior) on platforms where
    // native window exposure isn't wired up above.
    nfdwindowhandle_t win_handle{};
    NFD_GetNativeWindowFromGLFWWindow((GLFWwindow*)data.internal_ctx.window.ptr, &win_handle);

    if (flags & FileDialogFlag_Open) {
        nfdopendialogu8args_t args{};
        args.filterList = filter_list;
        args.filterCount = filter_count;
        args.parentWindow = win_handle;
        args.title = "Open File";
        result = NFD_OpenDialogU8_With(&out_path, &args);
    } else if (flags & FileDialogFlag_Save) {
        nfdsavedialogu8args_t args{};
        args.filterList = filter_list;
        args.filterCount = filter_count;
        args.parentWindow = win_handle;
        args.title = "Save File";
        result = NFD_SaveDialogU8_With(&out_path, &args);
    }

    if (result == NFD_OKAY) {
        int len = snprintf(str_buf, str_cap, "%s", out_path);
        if (flags & FileDialogFlag_Save) {
            // If the user is saving through the dialogue and there is no extension
            // In such case we append the first extension found in the filter (if supplied)
            str_t ext;
            str_t path = str_from_cstr(out_path);
            if (!extract_ext(&ext, path) && filter) {
                // get ext from supplied filter (first match)
                ext = filter;
                str_find_char(&ext.len, ext, ',');
                len += snprintf(str_buf + len, str_cap - len, "." STR_FMT, STR_ARG(ext));
            }
        }
        
        if (0 < len && len < (int)str_cap) {
            replace_char(str_buf, len, '\\', '/');
            return true;
        }

        MD_LOG_ERROR("snprintf failed");
        return false;
    } else if (result == NFD_ERROR) {
        MD_LOG_ERROR("%s\n", NFD_GetError());
    }
    /* fallthrough for NFD_CANCEL */
    return false;
}

}  // namespace application
