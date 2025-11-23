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
#include <SDL3/SDL.h>

#if MD_PLATFORM_WINDOWS
#include <SDL3/SDL_syswm.h>
#endif

#include <imgui.h>
#include <implot.h>

#include <app/imgui_impl_sdl3.h>
#include <app/imgui_impl_opengl3.h>

// Compressed fonts
#include <app/dejavu_sans_mono.inl>
#include <app/fa_solid.inl>
#include <app/IconsFontAwesome6.h>

#include <stdio.h> // snprintf
#include <stdlib.h> // free

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

bool initialize(Context* ctx, size_t width, size_t height, str_t title) {
    if (!SDL_Init(SDL_INIT_VIDEO | SDL_INIT_GAMEPAD)) {
        MD_LOG_ERROR("Error while initializing SDL: %s", SDL_GetError());
        return false;
    }

    if (width == 0 && height == 0) {
        int display_count = 0;
        SDL_DisplayID* displays = SDL_GetDisplays(&display_count);
        if (display_count > 0 && displays) {
            SDL_Rect usable_bounds;
            if (SDL_GetDisplayUsableBounds(displays[0], &usable_bounds)) {
                width  = (size_t)(usable_bounds.w * 0.9);
                height = (size_t)(usable_bounds.h * 0.8);
            }
            SDL_free(displays);
        }
    }

#if MD_PLATFORM_OSX
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG);
#endif

    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

    // Zero terminated
    str_t ztitle = str_copy(title, md_get_temp_allocator());
    SDL_Window* window = SDL_CreateWindow(ztitle.ptr, (int)width, (int)height, SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_HIGH_PIXEL_DENSITY);
    if (!window) {
        MD_LOG_ERROR("Could not create SDL window: %s", SDL_GetError());
        return false;
    }

    SDL_GLContext gl_context = SDL_GL_CreateContext(window);
    if (!gl_context) {
        MD_LOG_ERROR("Could not create OpenGL context: %s", SDL_GetError());
        SDL_DestroyWindow(window);
        return false;
    }
    
    SDL_GL_MakeCurrent(window, gl_context);
    SDL_GL_SetSwapInterval(1);
    
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

    int version = SDL_GetVersion();
    data.internal_ctx.gl_info.version.major = SDL_VERSIONNUM_MAJOR(version);
    data.internal_ctx.gl_info.version.minor = SDL_VERSIONNUM_MINOR(version);
    data.internal_ctx.gl_info.version.revision = SDL_VERSIONNUM_MICRO(version);

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
    // io.ConfigFlags |= ImGuiConfigFlags_DpiEnableScaleFonts;
    // io.ConfigDockingWithShift = true;
    io.ConfigWindowsMoveFromTitleBarOnly = true;

    //float xscale, yscale;
    //glfwGetWindowContentScale(window, &xscale, &yscale);
    //const float dpi_scale = (xscale + yscale) * 0.5f;

    // Default range:               0x0020 - 0x00FF.
    // Greek and Coptik:            0x0370 - 0x03FF
    // Superscripts and Subscripts: 0x2070 - 0x209F
    const ImWchar ranges_characters[] = {0x0020, 0x00FF, 0x0370, 0x03FF, 0x2070, 0x209F, 0};
    const ImWchar ranges_icons[] = {ICON_MIN_FA, ICON_MAX_FA, 0};
    const float font_size[] = {10, 12, 14, 16, 18, 20, 24, 32, 36, 40, 48};
    const char* font_names[] = {"10", "12", "14", "16", "18", "20", "24", "32", "36", "40", "48"};

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
    io.FontDefault = io.Fonts->Fonts[4]; // Set default to 18px

    if (!ImGui_ImplSDL3_InitForOpenGL(window, gl_context) ||
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
    SDL_GetWindowSizeInPixels(window, &w, &h);
    data.internal_ctx.framebuffer.width = w;
    data.internal_ctx.framebuffer.height = h;

    SDL_SetPointerProperty(SDL_GetWindowProperties(window), "context", &data.internal_ctx);
    
#if MD_PLATFORM_WINDOWS
    SDL_PropertiesID props = SDL_GetWindowProperties(window);
    HWND hwnd = (HWND)SDL_GetPointerProperty(props, SDL_PROP_WINDOW_WIN32_HWND_POINTER, NULL);
    if (hwnd) {
        HINSTANCE hinst = GetModuleHandle(NULL);
        HICON hIcon = LoadIcon(hinst, "VIAMD_ICON");
        SendMessage(hwnd, WM_SETICON, ICON_SMALL, (LPARAM)hIcon);
        SendMessage(hwnd, WM_SETICON, ICON_BIG, (LPARAM)hIcon);
    }
#endif

    MEMCPY(ctx, &data.internal_ctx, sizeof(Context));

    return true;
}

void shutdown(Context* ctx) {
    SDL_Window* window = (SDL_Window*)data.internal_ctx.window.ptr;
    SDL_GLContext gl_context = SDL_GL_GetCurrentContext();
    
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplSDL3_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    
    if (gl_context) {
        SDL_GL_DestroyContext(gl_context);
    }
    if (window) {
        SDL_DestroyWindow(window);
    }
    SDL_Quit();

    MEMSET(ctx, 0, sizeof(Context));
}

void update(Context* ctx) {
    SDL_Window* window = (SDL_Window*)data.internal_ctx.window.ptr;
    
    // Process SDL events
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        ImGui_ImplSDL3_ProcessEvent(&event);
        
        if (event.type == SDL_EVENT_QUIT) {
            data.internal_ctx.window.should_close = true;
        }
        
        if (event.type == SDL_EVENT_DROP_FILE) {
            Context* ctx = (Context*)SDL_GetPointerProperty(SDL_GetWindowProperties(window), "context", NULL);
            if (ctx && ctx->file_drop.callback && event.drop.data) {
                MD_LOG_DEBUG("User dropped file: '%s'", event.drop.data);
                str_t path = {event.drop.data, strlen(event.drop.data)};
                ctx->file_drop.callback(1, &path, ctx->file_drop.user_data);
            }
        }
    }

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplSDL3_NewFrame();
    ImGui::NewFrame();

    if (ctx->window.width != data.internal_ctx.window.width || ctx->window.height != data.internal_ctx.window.height) {
        SDL_SetWindowSize(window, ctx->window.width, ctx->window.height);
        data.internal_ctx.window.width = ctx->window.width;
        data.internal_ctx.window.height = ctx->window.height;
    }
    
    int w, h;
    SDL_GetWindowSizeInPixels(window, &w, &h);
    data.internal_ctx.framebuffer.width = w;
    data.internal_ctx.framebuffer.height = h;

    SDL_GetWindowSize(window, &w, &h);
    data.internal_ctx.window.width = w;
    data.internal_ctx.window.height = h;

    if (ctx->window.vsync != data.internal_ctx.window.vsync) {
        data.internal_ctx.window.vsync = ctx->window.vsync;
        SDL_GL_SetSwapInterval((int)ctx->window.vsync);
    }

    if (ctx->file_drop.callback != data.internal_ctx.file_drop.callback) {
        data.internal_ctx.file_drop.callback = ctx->file_drop.callback;
    }

    if (ctx->file_drop.user_data != data.internal_ctx.file_drop.user_data) {
        data.internal_ctx.file_drop.user_data = ctx->file_drop.user_data;
    }

    uint64_t t = SDL_GetTicks();
    double t_s = t / 1000.0;
    data.internal_ctx.timing.delta_s = (t_s - data.internal_ctx.timing.total_s);
    data.internal_ctx.timing.total_s = t_s;

    MEMCPY(ctx, &data.internal_ctx, sizeof(Context));
}

void render_imgui(Context* ctx) {
    (void)ctx;
    SDL_Window* window = (SDL_Window*)data.internal_ctx.window.ptr;
    SDL_GLContext gl_context = SDL_GL_GetCurrentContext();

    ImGui::Render();
    SDL_GL_MakeCurrent(window, gl_context);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    // Update and Render additional Platform Windows
    if (ImGui::GetIO().ConfigFlags & ImGuiConfigFlags_ViewportsEnable) {
        ImGui::UpdatePlatformWindows();
        ImGui::RenderPlatformWindowsDefault();
    }

    SDL_GL_MakeCurrent(window, gl_context);
}

void swap_buffers(Context* ctx) { 
    SDL_GL_SwapWindow((SDL_Window*)ctx->window.ptr); 
}

bool file_dialog(char* str_buf, size_t str_cap, FileDialogFlag flags, str_t filter) {
    struct DialogState {
        const char* result_path;
        bool completed;
    };
    
    static DialogState dialog_state = {NULL, false};
    dialog_state.result_path = NULL;
    dialog_state.completed = false;
    
    SDL_DialogFileFilter* sdl_filters = NULL;
    int filter_count = 0;
    
    // Parse filter string (e.g., "jpg,png,bmp")
    if (filter.ptr && filter.len > 0) {
        // Count comma-separated items
        filter_count = 1;
        for (size_t i = 0; i < filter.len; i++) {
            if (filter.ptr[i] == ',') filter_count++;
        }
        
        sdl_filters = (SDL_DialogFileFilter*)md_temp_push(filter_count * sizeof(SDL_DialogFileFilter));
        
        // Parse extensions
        const char* start = filter.ptr;
        int idx = 0;
        for (size_t i = 0; i <= filter.len; i++) {
            if (i == filter.len || filter.ptr[i] == ',') {
                size_t ext_len = (filter.ptr + i) - start;
                char* ext = (char*)md_temp_push(ext_len + 3); // "*.ext\0"
                snprintf(ext, ext_len + 3, "*.%.*s", (int)ext_len, start);
                
                sdl_filters[idx].name = ext;
                sdl_filters[idx].pattern = ext;
                idx++;
                start = filter.ptr + i + 1;
            }
        }
    }

    auto dialog_callback = [](void* userdata, const char* const* filelist, int filter_idx) {
        DialogState* state = (DialogState*)userdata;
        if (filelist && filelist[0]) {
            state->result_path = filelist[0];
        }
        state->completed = true;
    };

    if (flags & FileDialogFlag_Open) {
        SDL_ShowOpenFileDialog(
            dialog_callback,
            &dialog_state,
            NULL,  // parent window
            sdl_filters,
            filter_count,
            NULL,  // default location
            false  // allow_many
        );
    } else if (flags & FileDialogFlag_Save) {
        SDL_ShowSaveFileDialog(
            dialog_callback,
            &dialog_state,
            NULL,  // parent window
            sdl_filters,
            filter_count,
            NULL   // default location
        );
    } else {
        return false;
    }
    
    // Wait for dialog to complete (blocking event loop)
    while (!dialog_state.completed) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_EVENT_QUIT) {
                return false;
            }
        }
        SDL_Delay(10);
    }
    
    if (dialog_state.result_path && dialog_state.result_path[0] != '\0') {
        int len = snprintf(str_buf, str_cap, "%s", dialog_state.result_path);
        
        if (flags & FileDialogFlag_Save) {
            // If the user is saving and there is no extension, append the first filter extension
            str_t ext;
            str_t path = str_from_cstr(dialog_state.result_path);
            if (!extract_ext(&ext, path) && filter.ptr && filter.len > 0) {
                // get ext from supplied filter (first match)
                ext = filter;
                str_find_char(&ext.len, ext, ',');
                len += snprintf(str_buf + len, str_cap - len, "." STR_FMT, STR_ARG(ext));
            }
        }
        
        if (0 < len && (size_t)len < str_cap) {
            replace_char(str_buf, len, '\\', '/');
            return true;
        }
        
        MD_LOG_ERROR("snprintf failed or buffer too small");
        return false;
    }
    
    return false;
}

}  // namespace application
