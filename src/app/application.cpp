#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <app/application.h>

#include <core/md_log.h>
#include <core/md_common.h>
#include <core/md_platform.h>
#include <core/md_str.h>
#include <core/md_os.h>
#include <core/md_format.h>
#include <core/md_path.h>
#include <core/md_allocator.h>

#include <core/md_gpu.h>
#include <core/md_gpu_gl.h>

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#if MD_PLATFORM_WINDOWS
#define GLFW_EXPOSE_NATIVE_WIN32
#include <GLFW/glfw3native.h>
#elif MD_PLATFORM_OSX
#define GLFW_EXPOSE_NATIVE_COCOA
#include <GLFW/glfw3native.h>
#elif MD_PLATFORM_LINUX
#define GLFW_EXPOSE_NATIVE_X11
#include <GLFW/glfw3native.h>
#endif
#include <nfd.h>
#include <nfd_glfw3.h>

// NFD's GTK backend: on Linux, unless it is built for the xdg-desktop-portal, which defines NFD_PORTAL. GTK is linked
// in through nfd; the one function used here is declared rather than gtk.h included, which takes GTK's include paths.
#if MD_PLATFORM_LINUX && !defined(NFD_PORTAL)
#define VIAMD_NFD_GTK 1
extern "C" void gtk_disable_setlocale(void);
#else
#define VIAMD_NFD_GTK 0
#endif

#include <imgui.h>
#include <implot.h>
#include <implot3d.h>
#include <plot_input.h>

#include <app/imgui_impl_glfw.h>
#include <app/imgui_impl_md_gpu.h>

// Compressed fonts
#include <app/dejavu_sans_mono.inl>
#include <app/fa_solid.inl>
#include <app/IconsFontAwesome6.h>

#include <stdio.h> // snprintf
#include <stdint.h> // uint64_t
#include <stdlib.h> // free
#include <string.h> // strcmp
#include <locale.h> // setlocale

namespace application {

// Data
static struct {
    Context internal_ctx{};
} data;

// The decimal point is '.', always: in the workspace, in the script, in every number shown or typed in. What mdlib
// writes does not depend on the locale (md_snprintf), but the rest does: ImGui formats and parses its numeric fields
// through the C library, scanf included, and so do the labels and the .ini. Under a decimal comma locale 1.5 shows as
// "1,5" and a typed "1.5" reads as 1. So LC_NUMERIC stays "C", while the other categories follow the user (see
// initialize()).
//
// GTK, which NFD initializes on Linux, would switch it: gtk_init calls setlocale(LC_ALL, "") unless it is told not
// to, and it is. Nothing else is known to, and this check, once a frame, is for whatever does: loud, as that is a
// bug to find rather than something to put right in silence. It is put back as well, so numbers stay readable
// meanwhile. setlocale is process wide and not synchronized with the worker threads; it is only called here once
// something else has already switched the locale under them.
static void check_c_numeric_locale() {
    const char* cur = setlocale(LC_NUMERIC, NULL);
    if (cur && strcmp(cur, "C") == 0) {
        return;
    }
    MD_LOG_ERROR("Something switched the numeric locale to '%s', it is put back to \"C\": numbers would be shown and read with a decimal comma", cur ? cur : "unknown");
    setlocale(LC_NUMERIC, "C");
}

static void error_callback(int error, const char* description) {
    // The clipboard holding something other than text (an image, a file, nothing) is not a fault: GLFW reports it
    // when ImGui asks for the clipboard string on paste, and the paste simply does nothing.
    if (error == GLFW_FORMAT_UNAVAILABLE) {
        MD_LOG_DEBUG("%d: %s", error, description);
        return;
    }
    MD_LOG_ERROR("%d: %s\n", error, description);
}

// ## Presentation
//
// The window is shown through md_gpu: a surface (the window's swapchain) whose image ImGui's draw data is drawn into,
// the main view an image in that (main.cpp). Where md_gpu cannot present (an X session forwarded or virtualized without
// a Vulkan WSI, VirtualGL), or when asked to (VIAMD_PRESENT=gl, or gl-copy), through OpenGL instead: the same drawn into an md_gpu
// texture of the window's size and shown with one blit, the texture shared with GL where the two can (md_gpu_gl, the
// same GPU), else read back and uploaded. Only the blit is GL; its few entry points are loaded here.

#if MD_PLATFORM_WINDOWS
#define VIAMD_GLAPI __stdcall
#else
#define VIAMD_GLAPI
#endif

enum : uint32_t {
    VGL_VENDOR = 0x1F00, VGL_RENDERER = 0x1F01, VGL_VERSION = 0x1F02,
    VGL_TEXTURE_2D = 0x0DE1, VGL_RGBA8 = 0x8058, VGL_RGBA = 0x1908, VGL_UNSIGNED_BYTE = 0x1401,
    VGL_TEXTURE_MIN_FILTER = 0x2801, VGL_TEXTURE_MAG_FILTER = 0x2800, VGL_NEAREST = 0x2600,
    VGL_UNPACK_ALIGNMENT = 0x0CF5, VGL_UNPACK_ROW_LENGTH = 0x0CF2,
    VGL_READ_FRAMEBUFFER = 0x8CA8, VGL_DRAW_FRAMEBUFFER = 0x8CA9, VGL_COLOR_ATTACHMENT0 = 0x8CE0,
    VGL_FRAMEBUFFER_COMPLETE = 0x8CD5, VGL_COLOR_BUFFER_BIT = 0x4000,
};

struct GlFn {
    const unsigned char* (VIAMD_GLAPI* GetString)(uint32_t);
    void (VIAMD_GLAPI* GenTextures)(int, uint32_t*);
    void (VIAMD_GLAPI* DeleteTextures)(int, const uint32_t*);
    void (VIAMD_GLAPI* BindTexture)(uint32_t, uint32_t);
    void (VIAMD_GLAPI* TexImage2D)(uint32_t, int, int, int, int, int, uint32_t, uint32_t, const void*);
    void (VIAMD_GLAPI* TexSubImage2D)(uint32_t, int, int, int, int, int, uint32_t, uint32_t, const void*);
    void (VIAMD_GLAPI* TexParameteri)(uint32_t, uint32_t, int);
    void (VIAMD_GLAPI* PixelStorei)(uint32_t, int);
    void (VIAMD_GLAPI* GenFramebuffers)(int, uint32_t*);
    void (VIAMD_GLAPI* DeleteFramebuffers)(int, const uint32_t*);
    void (VIAMD_GLAPI* BindFramebuffer)(uint32_t, uint32_t);
    void (VIAMD_GLAPI* FramebufferTexture2D)(uint32_t, uint32_t, uint32_t, uint32_t, int);
    uint32_t (VIAMD_GLAPI* CheckFramebufferStatus)(uint32_t);
    void (VIAMD_GLAPI* BlitFramebuffer)(int, int, int, int, int, int, int, int, uint32_t, uint32_t);
    void (VIAMD_GLAPI* ClearColor)(float, float, float, float);
    void (VIAMD_GLAPI* Clear)(uint32_t);
};

static void* gl_loader(const char* name) { return (void*)glfwGetProcAddress(name); }

static bool gl_load(GlFn* f) {
    bool ok = true;
#define VGL_LOAD(x) f->x = (decltype(f->x))glfwGetProcAddress("gl" #x); if (!f->x) { MD_LOG_ERROR("OpenGL lacks gl" #x); ok = false; }
    VGL_LOAD(GetString) VGL_LOAD(GenTextures) VGL_LOAD(DeleteTextures) VGL_LOAD(BindTexture) VGL_LOAD(TexImage2D)
    VGL_LOAD(TexSubImage2D) VGL_LOAD(TexParameteri) VGL_LOAD(PixelStorei) VGL_LOAD(GenFramebuffers)
    VGL_LOAD(DeleteFramebuffers) VGL_LOAD(BindFramebuffer) VGL_LOAD(FramebufferTexture2D) VGL_LOAD(CheckFramebufferStatus)
    VGL_LOAD(BlitFramebuffer) VGL_LOAD(ClearColor) VGL_LOAD(Clear)
#undef VGL_LOAD
    return ok;
}

static struct {
    md_gpu_device_t  device = nullptr;
    md_gpu_stream_t  stream = nullptr;          // GRAPHICS: draws, presents
    md_gpu_surface_t surface = nullptr;         // the window's, unless through OpenGL
    bool             vsync = true;              // the surfaces' present mode
    md_gpu_sync_t    frame_done = {};           // the last frame's work: the next waits for it (update)
    // Through OpenGL
    bool                gl = false;
    GlFn                fn = {};
    md_gpu_texture_t    frame = nullptr;        // what is shown: drawn by md_gpu, of the window's size
    uint32_t            frame_w = 0, frame_h = 0;
    bool                share = false;          // the frame shared with GL, else read back and uploaded
    md_gpu_gl_texture_t shared = {};
    md_gpu_mem_t        readback = {};
    uint32_t            gl_tex = 0;             // the read back frame in GL
    uint32_t            fbo = 0;                // what the blit reads
} present;

static md_gpu_present_mode_t present_mode(bool vsync) {
    return vsync ? MD_GPU_PRESENT_VSYNC : MD_GPU_PRESENT_IMMEDIATE;
}

static md_gpu_surface_t make_surface(GLFWwindow* window, int width, int height, bool vsync, const char* label) {
    md_gpu_surface_desc_t d = {};
#if MD_PLATFORM_WINDOWS
    d.system = MD_GPU_WINDOW_WIN32;
    d.window = (void*)glfwGetWin32Window(window);
#elif MD_PLATFORM_OSX
    d.system = MD_GPU_WINDOW_COCOA;
    d.window = (void*)glfwGetCocoaView(window);
#else
    d.system  = MD_GPU_WINDOW_X11;
    d.window  = (void*)(uintptr_t)glfwGetX11Window(window);
    d.display = (void*)glfwGetX11Display();
#endif
    d.width  = (uint32_t)MD_MAX(width, 1);
    d.height = (uint32_t)MD_MAX(height, 1);
    d.format = MD_GPU_FORMAT_BGRA8_UNORM;
    d.present_mode = present_mode(vsync);
    d.label  = label;
    return md_gpu_surface_create(present.device, &d);
}

static void gl_frame_free() {
    if (present.share && present.shared.name) md_gpu_gl_release(&present.shared);
    present.shared = {};
    if (present.frame) md_gpu_texture_destroy(present.frame);
    present.frame = nullptr;
    if (present.readback.gpu) md_gpu_free(present.stream, present.readback.gpu);
    present.readback = {};
    if (present.gl_tex) present.fn.DeleteTextures(1, &present.gl_tex);
    present.gl_tex = 0;
    present.frame_w = present.frame_h = 0;
}

// The frame of the size, and GL's view of it
static bool gl_frame(uint32_t w, uint32_t h) {
    if (present.frame && present.frame_w == w && present.frame_h == h) return true;
    gl_frame_free();
    md_gpu_texture_desc_t d = {};
    d.type   = MD_GPU_TEX_2D;
    d.format = MD_GPU_FORMAT_RGBA8_UNORM;
    d.usage  = MD_GPU_TEX_RENDER_TARGET | (present.share ? MD_GPU_TEX_SHARED : 0);
    d.width  = w;
    d.height = h;
    d.label  = "window";
    present.frame = md_gpu_texture_create(present.stream, &d);
    if (present.share && !(present.frame && md_gpu_gl_import(&present.shared, present.frame, gl_loader))) {
        MD_LOG_INFO("The window's image cannot be shared with OpenGL (%s): it is read back and uploaded instead", md_gpu_last_error());
        gl_frame_free();
        present.share = false;
        d.usage = MD_GPU_TEX_RENDER_TARGET;
        present.frame = md_gpu_texture_create(present.stream, &d);
    }
    if (!present.frame) {
        MD_LOG_ERROR("The window's image (%u x %u) could not be made: %s", w, h, md_gpu_last_error());
        return false;
    }
    uint32_t name = 0, target = VGL_TEXTURE_2D;
    if (present.share) {
        name   = present.shared.name;
        target = present.shared.target;
    } else {
        present.readback = md_gpu_malloc(present.stream, MD_GPU_MEM_HOST_READ, (size_t)w * h * 4);
        if (!present.readback.cpu) {
            MD_LOG_ERROR("The window's image cannot be read back: %s", md_gpu_last_error());
            gl_frame_free();
            return false;
        }
        present.fn.GenTextures(1, &present.gl_tex);
        present.fn.BindTexture(VGL_TEXTURE_2D, present.gl_tex);
        present.fn.TexImage2D(VGL_TEXTURE_2D, 0, (int)VGL_RGBA8, (int)w, (int)h, 0, VGL_RGBA, VGL_UNSIGNED_BYTE, nullptr);
        present.fn.TexParameteri(VGL_TEXTURE_2D, VGL_TEXTURE_MIN_FILTER, (int)VGL_NEAREST);
        present.fn.TexParameteri(VGL_TEXTURE_2D, VGL_TEXTURE_MAG_FILTER, (int)VGL_NEAREST);
        present.fn.BindTexture(VGL_TEXTURE_2D, 0);
        name = present.gl_tex;
    }
    if (!present.fbo) present.fn.GenFramebuffers(1, &present.fbo);
    present.fn.BindFramebuffer(VGL_READ_FRAMEBUFFER, present.fbo);
    present.fn.FramebufferTexture2D(VGL_READ_FRAMEBUFFER, VGL_COLOR_ATTACHMENT0, target, name, 0);
    const bool complete = present.fn.CheckFramebufferStatus(VGL_READ_FRAMEBUFFER) == VGL_FRAMEBUFFER_COMPLETE;
    present.fn.BindFramebuffer(VGL_READ_FRAMEBUFFER, 0);
    if (!complete) {
        MD_LOG_ERROR("OpenGL cannot read the window's image");
        gl_frame_free();
        return false;
    }
    present.frame_w = w;
    present.frame_h = h;
    return true;
}

// The draw data into the frame, the frame into the window's back buffer (its rows from the bottom)
static void gl_present(ImDrawData* dd, uint32_t w, uint32_t h) {
    if (!w || !h || !gl_frame(w, h)) {
        present.fn.ClearColor(0, 0, 0, 1);
        present.fn.Clear(VGL_COLOR_BUFFER_BIT);
        return;
    }
    const float black[4] = {0, 0, 0, 1};
    if (present.share) md_gpu_texture_share_acquire(present.stream, present.frame);
    ImGui_ImplMdGpu_RenderDrawData(dd, present.frame, black);
    if (present.share) {
        md_gpu_texture_share_release(present.stream, present.frame);
        md_gpu_gl_begin_read(&present.shared);
    } else {
        md_gpu_copy_from_texture(present.stream, present.readback.gpu, present.frame, nullptr);
        md_gpu_stream_sync(present.stream);
        present.fn.BindTexture(VGL_TEXTURE_2D, present.gl_tex);
        present.fn.PixelStorei(VGL_UNPACK_ALIGNMENT, 4);
        present.fn.PixelStorei(VGL_UNPACK_ROW_LENGTH, 0);
        present.fn.TexSubImage2D(VGL_TEXTURE_2D, 0, 0, 0, (int)w, (int)h, VGL_RGBA, VGL_UNSIGNED_BYTE, present.readback.cpu);
        present.fn.BindTexture(VGL_TEXTURE_2D, 0);
    }
    present.fn.BindFramebuffer(VGL_READ_FRAMEBUFFER, present.fbo);
    present.fn.BindFramebuffer(VGL_DRAW_FRAMEBUFFER, 0);
    present.fn.BlitFramebuffer(0, 0, (int)w, (int)h, 0, (int)h, (int)w, 0, VGL_COLOR_BUFFER_BIT, VGL_NEAREST);
    present.fn.BindFramebuffer(VGL_READ_FRAMEBUFFER, 0);
    if (present.share) md_gpu_gl_end_read(&present.shared);
}

// The other windows ImGui makes (multiple viewports): a surface each
static void viewport_create(ImGuiViewport* vp) {
    int w = 0, h = 0;
    glfwGetFramebufferSize((GLFWwindow*)vp->PlatformHandle, &w, &h);
    vp->RendererUserData = make_surface((GLFWwindow*)vp->PlatformHandle, w, h, present.vsync, "imgui viewport");
    if (!vp->RendererUserData) MD_LOG_ERROR("A window's surface could not be made: %s", md_gpu_last_error());
}

static void viewport_destroy(ImGuiViewport* vp) {
    if (vp->RendererUserData) md_gpu_surface_destroy((md_gpu_surface_t)vp->RendererUserData);
    vp->RendererUserData = nullptr;
}

static void viewport_set_size(ImGuiViewport* vp, ImVec2) {
    int w = 0, h = 0;
    glfwGetFramebufferSize((GLFWwindow*)vp->PlatformHandle, &w, &h);
    if (vp->RendererUserData) md_gpu_surface_resize((md_gpu_surface_t)vp->RendererUserData, (uint32_t)MD_MAX(w, 1), (uint32_t)MD_MAX(h, 1));
}

static void viewport_render(ImGuiViewport* vp, void*) {
    md_gpu_surface_t s = (md_gpu_surface_t)vp->RendererUserData;
    if (!s) return;
    md_gpu_texture_t back = md_gpu_surface_acquire(present.stream, s);
    if (!back) return;
    const float black[4] = {0, 0, 0, 1};
    ImGui_ImplMdGpu_RenderDrawData(vp->DrawData, back, (vp->Flags & ImGuiViewportFlags_NoRendererClear) ? nullptr : black);
    md_gpu_surface_present(present.stream, s);
}

bool initialize(Context* ctx, size_t width, size_t height, md_str_t title) {
    glfwSetErrorCallback(error_callback);
    if (!glfwInit()) {
        MD_LOG_ERROR("Error while initializing glfw.");
        return false;
    }

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

    // The device everything is drawn with; VIAMD_GPU_VALIDATION set: with the backend's validation (Vulkan's layers)
    md_gpu_device_desc_t dd = {};
    dd.enable_validation = getenv("VIAMD_GPU_VALIDATION") != nullptr;
    present.device = md_gpu_device_create(&dd);
    if (!present.device) {
        MD_LOG_ERROR("No GPU device (md_gpu): %s", md_gpu_last_error());
        glfwTerminate();
        return false;
    }
    md_gpu_device_info_t info = {};
    md_gpu_device_info(present.device, &info);
    present.stream = md_gpu_stream_default(present.device, MD_GPU_STREAM_GRAPHICS);
    if (!info.supports_graphics || !present.stream) {
        MD_LOG_ERROR("The GPU device '%s' does not draw (md_gpu graphics)", info.name);
        md_gpu_device_destroy(present.device);
        present.device = nullptr;
        glfwTerminate();
        return false;
    }

    // Zero terminated
    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    md_str_t ztitle = md_alloc_str(md_temp_allocator(temp), title);

    // The window: presented by md_gpu, or through OpenGL where md_gpu cannot or is asked not to (VIAMD_PRESENT=gl;
    // gl-copy: the image read back and uploaded, not shared with GL)
    const char* env = getenv("VIAMD_PRESENT");
    const bool want_copy = env && strcmp(env, "gl-copy") == 0;
    const bool want_gl = want_copy || (env && strcmp(env, "gl") == 0);
    GLFWwindow* window = nullptr;
    if (!want_gl && info.supports_present) {
        glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
        window = glfwCreateWindow((int)width, (int)height, ztitle.ptr, NULL, NULL);
        if (window) {
            int w, h;
            glfwGetFramebufferSize(window, &w, &h);
            present.surface = make_surface(window, w, h, present.vsync, "window");
            if (!present.surface) {
                MD_LOG_INFO("md_gpu cannot present to the window (%s): it is shown through OpenGL instead", md_gpu_last_error());
                glfwDestroyWindow(window);
                window = nullptr;
            }
        }
    } else if (!want_gl) {
        MD_LOG_INFO("The GPU device '%s' cannot present: the window is shown through OpenGL", info.name);
    }
    if (!window) {
        present.gl = true;
        glfwDefaultWindowHints();
        glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_API);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
        window = glfwCreateWindow((int)width, (int)height, ztitle.ptr, NULL, NULL);
        if (!window) {
            MD_LOG_ERROR("Could not create a window, with or without OpenGL.");
            return false;
        }
        glfwMakeContextCurrent(window);
        glfwSwapInterval(present.vsync ? 1 : 0);
        if (!gl_load(&present.fn)) {
            MD_LOG_ERROR("Could not load the OpenGL functions the window is shown with.");
            return false;
        }
        present.share = info.supports_sharing && !want_copy;
        auto gl_str = [](uint32_t name) -> const char* {
            const char* s = (const char*)present.fn.GetString(name);
            return s ? s : "(unavailable)";
        };
        MD_LOG_INFO("Shown through OpenGL: %s, %s, %s", gl_str(VGL_VENDOR), gl_str(VGL_RENDERER), gl_str(VGL_VERSION));
    }
    MD_LOG_INFO("Drawn with '%s' (%s)%s", info.name, info.driver, present.gl ? "" : ", presented by it");

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImPlot3D::CreateContext();
    // The mouse as in every view (plot_input.h)
    plot_input_house(&ImPlot::GetInputMap());
    plot_input_house(&ImPlot3D::GetInputMap());
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard Controls
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;  // Enable Docking

#if VIAMD_IMGUI_ENABLE_VIEWPORTS
    // The other windows are presented by md_gpu: not through OpenGL
    if (!present.gl) io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;  // Enable Multi-Viewport / Platform Windows
#endif

    io.ConfigWindowsMoveFromTitleBarOnly = true;

    // Default range:               0x0020 - 0x00FF.
    // Greek and Coptik:            0x0370 - 0x03FF
    // General punctuation (dashes, quotes, ...): 0x2010 - 0x2027
    // Superscripts and Subscripts: 0x2070 - 0x209F
    // Nabla:                       0x2207
    const ImWchar ranges_characters[] = {0x0020, 0x00FF, 0x0370, 0x03FF, 0x2010, 0x2027, 0x2070, 0x209F, 0x2207, 0x2207, 0};
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

    const bool platform_ok = present.gl ? ImGui_ImplGlfw_InitForOpenGL(window, false) : ImGui_ImplGlfw_InitForOther(window, false);
    if (!platform_ok || !ImGui_ImplMdGpu_Init(present.device, present.stream)) {
        MD_LOG_ERROR("Failed to initialize ImGui's backends");
        return false;
    }
    if (!present.gl) {
        ImGuiPlatformIO& pio = ImGui::GetPlatformIO();
        pio.Renderer_CreateWindow  = viewport_create;
        pio.Renderer_DestroyWindow = viewport_destroy;
        pio.Renderer_SetWindowSize = viewport_set_size;
        pio.Renderer_RenderWindow  = viewport_render;
        io.BackendFlags |= ImGuiBackendFlags_RendererHasViewports;
    }

    data.internal_ctx.window.ptr = window;
    data.internal_ctx.window.title = title;
    data.internal_ctx.window.width = (int)width;
    data.internal_ctx.window.height = (int)height;
    data.internal_ctx.window.vsync = present.vsync;
    data.internal_ctx.gpu.device = present.device;
    data.internal_ctx.gpu.stream = present.stream;
    data.internal_ctx.gpu.gl = present.gl;

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
        MD_ASSERT(ctx);

        for (int i = 0; i < num_files; ++i) {
            MD_LOG_DEBUG("User dropped file: '%s'", paths[i]);
        }

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        md_str_t* str_paths = md_temp_alloc_array(temp, md_str_t, num_files);
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

#if VIAMD_NFD_GTK
    // gtk_init, in NFD_Init, would call setlocale(LC_ALL, "") and take LC_NUMERIC with it. GTK is kept from touching
    // the locale at all, and VIAMD sets it: the user's in every category but the numbers, so the dialog keeps its
    // language, its sorting of names and its dates. This holds whenever GTK gets initialized, not only in NFD_Init.
    gtk_disable_setlocale();
    setlocale(LC_ALL, "");
    setlocale(LC_NUMERIC, "C");
#endif
    if (NFD_Init() != NFD_OKAY) {
        MD_LOG_ERROR("Failed to initialize NativeFileDialog: %s", NFD_GetError());
    }

    md_memcpy(ctx, &data.internal_ctx, sizeof(Context));

    return true;
}

void shutdown(Context* ctx) {
    NFD_Quit();

    // What is in flight done, before what it uses goes
    if (present.stream) md_gpu_stream_sync(present.stream);
    if (ImGui::GetIO().ConfigFlags & ImGuiConfigFlags_ViewportsEnable) ImGui::DestroyPlatformWindows();
    ImGui_ImplMdGpu_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot3D::DestroyContext();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    if (present.gl) {
        gl_frame_free();
        if (present.fbo) present.fn.DeleteFramebuffers(1, &present.fbo);
        present.fbo = 0;
    }
    if (present.surface) md_gpu_surface_destroy(present.surface);
    present.surface = nullptr;
    if (present.device) md_gpu_device_destroy(present.device);    // after the window's last present, before the window
    present.device = nullptr;
    present.stream = nullptr;
    present.frame_done = {};
    glfwDestroyWindow((GLFWwindow*)data.internal_ctx.window.ptr);
    glfwTerminate();

    md_memset(ctx, 0, sizeof(Context));
}

void update(Context* ctx) {
    // The previous frame's work done before this one's input is taken: a frame ahead of the GPU at most. Acquire
    // alone does not hold the CPU back where the driver hands an image back early (NVIDIA): frames queued up to
    // md_gpu's limit of acquires, each a frame of the latency from input to the screen.
    md_gpu_sync_wait(present.frame_done);
    present.frame_done = {};
    glfwPollEvents();
    check_c_numeric_locale();

    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    GLFWwindow* window = (GLFWwindow*)data.internal_ctx.window.ptr;
    if (ctx->window.width != data.internal_ctx.window.width || ctx->window.height != data.internal_ctx.window.height) {
        glfwSetWindowSize(window, ctx->window.width, ctx->window.height);
        data.internal_ctx.window.width = ctx->window.width;
        data.internal_ctx.window.height = ctx->window.height;
    }
    int w, h;
    glfwGetFramebufferSize(window, &w, &h);
    data.internal_ctx.framebuffer.width = w;
    data.internal_ctx.framebuffer.height = h;
    if (present.surface && w > 0 && h > 0) md_gpu_surface_resize(present.surface, (uint32_t)w, (uint32_t)h);

    glfwGetWindowSize(window, &w, &h);
    data.internal_ctx.window.width = w;
    data.internal_ctx.window.height = h;

    // The vertical sync: GL's swap interval, or the surfaces' present mode, changed on them: a new surface for a
    // window whose old one is still being destroyed is refused by some drivers (VK_ERROR_NATIVE_WINDOW_IN_USE_KHR)
    if (ctx->window.vsync != data.internal_ctx.window.vsync) {
        data.internal_ctx.window.vsync = ctx->window.vsync;
        present.vsync = ctx->window.vsync;
        if (present.gl) {
            glfwSwapInterval((int)ctx->window.vsync);
        } else {
            if (present.surface) md_gpu_surface_set_present_mode(present.surface, present_mode(present.vsync));
            for (ImGuiViewport* vp : ImGui::GetPlatformIO().Viewports) {
                if (vp->RendererUserData) md_gpu_surface_set_present_mode((md_gpu_surface_t)vp->RendererUserData, present_mode(present.vsync));
            }
        }
    }

    if (ctx->file_drop.callback != data.internal_ctx.file_drop.callback) {
        data.internal_ctx.file_drop.callback = ctx->file_drop.callback;
    }

    if (ctx->file_drop.user_data != data.internal_ctx.file_drop.user_data) {
        data.internal_ctx.file_drop.user_data = ctx->file_drop.user_data;
    }

    data.internal_ctx.window.should_close = (bool)glfwWindowShouldClose(window);

    double t = glfwGetTime();
    data.internal_ctx.timing.delta_s = (t - data.internal_ctx.timing.total_s);
    data.internal_ctx.timing.total_s = t;

    md_memcpy(ctx, &data.internal_ctx, sizeof(Context));
}

void render_imgui(Context* ctx) {
    (void)ctx;
    ImGui::Render();
    ImDrawData* dd = ImGui::GetDrawData();
    const uint32_t w = data.internal_ctx.framebuffer.width, h = data.internal_ctx.framebuffer.height;
    if (present.gl) {
        gl_present(dd, w, h);
    } else if (present.surface && w > 0 && h > 0) {
        // Acquire waits for the image: what paces the frames with vertical sync. NULL: minimised
        md_gpu_texture_t back = md_gpu_surface_acquire(present.stream, present.surface);
        if (back) {
            const float black[4] = {0, 0, 0, 1};
            ImGui_ImplMdGpu_RenderDrawData(dd, back, black);
            md_gpu_surface_present(present.stream, present.surface);
        }
    }

    // Update and Render additional Platform Windows
    if (ImGui::GetIO().ConfigFlags & ImGuiConfigFlags_ViewportsEnable) {
        ImGui::UpdatePlatformWindows();
        ImGui::RenderPlatformWindowsDefault();
    }
    if (present.stream) present.frame_done = md_gpu_stream_record(present.stream);
}

void swap_buffers(Context* ctx) {
    if (present.gl) glfwSwapBuffers((GLFWwindow*)ctx->window.ptr);
}

bool read_frame(Context* ctx, uint32_t* rgba, uint32_t width, uint32_t height) {
    (void)ctx;
    ImDrawData* dd = ImGui::GetDrawData();
    if (!dd || !dd->Valid || !rgba || !width || !height || !present.stream) return false;
    md_gpu_texture_desc_t d = {};
    d.type   = MD_GPU_TEX_2D;
    d.format = MD_GPU_FORMAT_RGBA8_UNORM;
    d.usage  = MD_GPU_TEX_RENDER_TARGET;
    d.width  = width;
    d.height = height;
    d.label  = "window capture";
    md_gpu_texture_t tex = md_gpu_texture_create(present.stream, &d);
    const size_t bytes = (size_t)width * height * 4;
    md_gpu_mem_t rb = tex ? md_gpu_malloc(present.stream, MD_GPU_MEM_HOST_READ, bytes) : md_gpu_mem_t{};
    const float black[4] = {0, 0, 0, 1};
    bool ok = rb.cpu && ImGui_ImplMdGpu_RenderDrawData(dd, tex, black) && md_gpu_copy_from_texture(present.stream, rb.gpu, tex, nullptr);
    md_gpu_stream_sync(present.stream);
    if (ok) md_memcpy(rgba, rb.cpu, bytes);
    else MD_LOG_ERROR("The window could not be read: %s", md_gpu_last_error());
    if (rb.gpu) md_gpu_free(present.stream, rb.gpu);
    if (tex) md_gpu_texture_destroy(tex);
    return ok;
}

bool file_dialog(char* str_buf, size_t str_cap, FileDialogFlag flags, md_str_t filter) {    
    md_temp_scope_t temp = md_temp_begin();
    nfdu8char_t* out_path = NULL;
    MD_DEFER {
        md_temp_end(temp);
        if (out_path) NFD_FreePathU8(out_path);
    };

    nfdresult_t result = NFD_ERROR;

    // Zero terminated variant
    md_str_t zfilt = md_alloc_str(md_temp_allocator(temp), filter);

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
            md_str_t path = md_str_from_cstr(out_path);
            md_str_t ext = md_path_ext(path);
            if (md_str_empty(ext) && filter) {
                // get ext from supplied filter (first match)
                md_str_cut(&ext, NULL, filter, md_str_find_char(filter, ','));
                len += snprintf(str_buf + len, str_cap - len, "." MD_STR_FMT, MD_STR_ARG(ext));
            }
        }
        
        if (0 < len && len < (int)str_cap) {
            md_chars_replace(str_buf, len, '\\', '/');
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
