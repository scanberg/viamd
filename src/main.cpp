#include <core/md_compiler.h>

#if MD_COMPILER_MSVC
#   ifndef _CRT_SECURE_NO_WARNINGS
#       define _CRT_SECURE_NO_WARNINGS
#   endif
#   pragma warning( disable : 26812 4244 )
#endif

#define IMGUI_DEFINE_MATH_OPERATORS

#include <md_util.h>
#include <md_gl.h>
#include <md_gfx.h>
#include <md_filter.h>
#include <md_script.h>
#include <md_system.h>
#include <md_xvg.h>
#include <md_csv.h>
#include <md_lammps.h>
#include <md_pdb.h>
#include <md_gto_int.h>

#include <core/md_log.h>
#include <core/md_str.h>
#include <core/md_array.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_ring_allocator.h>
#include <core/md_tracking_allocator.h>
#include <core/md_simd.h>
#include <core/md_os.h>
#include <core/md_unit.h>
#include <core/md_str_builder.h>
#include <core/md_parse.h>
#include <core/md_hash.h>

#include <gfx/gl.h>
#include <gfx/gl_utils.h>
#include <gfx/camera.h>
#include <gfx/camera_utils.h>
#include <gfx/immediate_draw_utils.h>
#include <gfx/postprocessing_utils.h>
#include <gfx/volumerender_utils.h>

#include <resource_path.h>
#include <app_settings.h>
#include <display_units.h>
#include <imgui_widgets.h>
#include <implot_widgets.h>
#include <task_system.h>
#include <color_utils.h>
#include <image.h>

#include <imgui.h>
#include <imgui_internal.h>
#include <app/imgui_impl_opengl3.h>

#include <implot.h>
#include <implot_internal.h>

#include <float.h>

#include <stdio.h>
#include <cmath>
#include <string>
#include <algorithm>
#include <chrono>

#include <viamd.h>
#include <viamd_icon.inl>
#include <script_reference.h>
#include <viamd_event.h>
#include <event.h>

#define EXPERIMENTAL_GFX_API 0
#define COMPILATION_TIME_DELAY_IN_SECONDS 1.0
#define NOTIFICATION_DISPLAY_TIME_IN_SECONDS 5.0
#define MEASURE_EVALUATION_TIME 1
#define VIAMD_RECOMPUTE_ORBITAL_PER_FRAME 0
#define VIS_FLAGS (MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY | MD_SCRIPT_VISUALIZE_TEXT)

// Global data for application
static md_allocator_i* frame_alloc = 0; // Linear allocator for scratch data which only is valid for the frame and then is reset
static md_allocator_i* persistent_alloc = 0;

static bool use_gfx = false;

constexpr str_t shader_output_snippet = STR_LIT(R"(
layout(location = 0) out vec4 out_color;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec4 out_velocity;
layout(location = 3) out vec4 out_atom_index;

vec2 encode_normal (vec3 n) {
   float p = sqrt(n.z * 8 + 8);
   return n.xy / p + 0.5;
}

vec4 encode_index(uint index) {
    return vec4(
        (index & 0x000000FFU) >> 0U,
        (index & 0x0000FF00U) >> 8U,
        (index & 0x00FF0000U) >> 16U,
        (index & 0xFF000000U) >> 24U) / 255.0;
}

vec2 compute_ss_vel(vec3 view_coord, vec3 view_vel) {
    vec3 prev_view_coord = view_coord - view_vel;
    vec4 prev_clip_coord = u_curr_view_to_prev_clip * vec4(prev_view_coord, 1);

    vec4 clip_coord = u_view_to_clip * vec4(view_coord, 1);
    vec2 curr_ndc = clip_coord.xy / clip_coord.w;
    vec2 prev_ndc = prev_clip_coord.xy / prev_clip_coord.w;
    vec2 ss_vel = (curr_ndc - prev_ndc) * 0.5 + (u_jitter_uv.xy - u_jitter_uv.zw);
    return ss_vel;
}

void write_fragment(vec3 view_coord, vec3 view_vel, vec3 view_normal, vec4 color, uint atom_index) {
   out_color  = color;
   out_normal = vec4(encode_normal(view_normal), 0, 0);
   out_velocity = vec4(compute_ss_vel(view_coord, view_vel), 0, 0);
   out_atom_index = encode_index(atom_index);
}
)");

constexpr str_t shader_output_snippet_lean_and_mean = STR_LIT(R"(
void write_fragment(vec3 view_coord, vec3 view_vel, vec3 view_normal, vec4 color, uint atom_index) {
}
)");


static double frame_to_time(double frame, const ApplicationState& data) {
    const int64_t num_frames = md_array_size(data.timeline.x_values);
    ASSERT(num_frames);
    const int64_t f0 = CLAMP((int64_t)frame + 0, 0, num_frames - 1);
    const int64_t f1 = CLAMP((int64_t)frame + 1, 0, num_frames - 1);
    return lerp(data.timeline.x_values[f0], data.timeline.x_values[f1], fract(frame));
}

// Try to map time t back into frame
static double time_to_frame(double time, const md_array(float) frame_times) {
    const int64_t num_frames = md_array_size(frame_times);
    if (!num_frames) return 0.0;

    const double beg = frame_times[0];
    const double end = frame_times[num_frames - 1];
    time = CLAMP(time, beg, end);

    // Estimate the frame
    const double frame_est = CLAMP(((time - beg) / (end-beg)) * (num_frames - 1), 0, num_frames - 1);

    int64_t prev_frame_idx = CLAMP((int64_t)frame_est,     0, num_frames - 1);
    int64_t next_frame_idx = CLAMP((int64_t)frame_est + 1, 0, num_frames - 1);

    if (time < (double)frame_times[prev_frame_idx]) {
        // Linear search down
        for (prev_frame_idx = prev_frame_idx - 1; prev_frame_idx >= 0; --prev_frame_idx) {
            next_frame_idx = prev_frame_idx + 1;
            if ((double)frame_times[prev_frame_idx] <= time && time <= (double)frame_times[next_frame_idx])
                break;
        }
    }
    else if (time > (double)frame_times[next_frame_idx]) {
        // Linear search up
        for (next_frame_idx = next_frame_idx + 1; next_frame_idx < num_frames; ++next_frame_idx) {
            prev_frame_idx = next_frame_idx - 1;
            if ((double)frame_times[prev_frame_idx] <= time && time <= (double)frame_times[next_frame_idx])
                break;
        }
    }

    // Compute true fraction between timestamps
    double t = (time - (double)frame_times[prev_frame_idx]) / ((double)frame_times[next_frame_idx] - (double)frame_times[prev_frame_idx]);
    t = CLAMP(t, 0.0, 1.0);

    // Compose frame value (base + fraction)
    return (double)prev_frame_idx + t;
}

static void update_view_param(ApplicationState* state);

//static void update_density_volume_texture(ApplicationState* state);

static void render(ApplicationState* state);
static void render_scene(ApplicationState* state, bool pip);
static void draw_representations_opaque(ApplicationState* state);
static void draw_representations_opaque_lean_and_mean(ApplicationState* state, uint32_t mask = 0xFFFFFFFFU);
// Returns true if anything was drawn into gbuffer.tex.transparency_hdr this frame
static bool draw_representations_transparent(ApplicationState* state);

static void draw_load_dataset_window(ApplicationState* state);
static void draw_main_menu(ApplicationState* state);
static void draw_context_popup(ApplicationState* state, const PickingHit& hit);
static void draw_selection_query_window(ApplicationState* state);
static void draw_selection_grow_window(ApplicationState* state);
static void draw_animation_window(ApplicationState* state);
static void draw_representations_window(ApplicationState* state);
static void draw_timeline_window(ApplicationState* state);
static void draw_distribution_window(ApplicationState* state);
static void draw_async_task_window(ApplicationState* state);
static void draw_script_editor_window(ApplicationState* state);
static void draw_script_reference_window(ApplicationState* state);
static void open_script_reference(ApplicationState* state, str_t topic, bool take_focus = true);
static void draw_coordinate_system_widget_window(ViewTransform* target, const ViewTransform& current);
static void draw_color_legend_windows(const ApplicationState& state);

static void draw_debug_window(ApplicationState* state);
static void draw_property_export_window(ApplicationState* state);
static void draw_structure_export_window(ApplicationState* state);
static void draw_notifications_window();

static void update_md_buffers(ApplicationState* state);

static bool export_xvg(const float* column_data[], const char* column_labels[], size_t num_columns, size_t num_rows, str_t filename);
static bool export_csv(const float* column_data[], const char* column_labels[], size_t num_columns, size_t num_rows, str_t filename);

static void create_screenshot(str_t path);

static void movie_recording_start(ApplicationState* state);
static void movie_recording_stop(ApplicationState* state);
static void movie_shutdown(ApplicationState* state);
static void movie_add_keyframe(ApplicationState* state);
static void movie_add_keyframe_from_view(ApplicationState* state);
static void movie_copy_keyframe_at_playhead(ApplicationState* state);
static void movie_paste_keyframe(ApplicationState* state);
static void update_movie_recording(ApplicationState* state);
static void update_movie_preview(ApplicationState* state);
static void update_movie_history(ApplicationState* state);
static void movie_undo(ApplicationState* state);
static void movie_redo(ApplicationState* state);
static void draw_movie_window(ApplicationState* state);
static void draw_movie_timeline_panel(ApplicationState* state);
static void draw_movie_recording_banner(ApplicationState* state);
static void movie_draw_camera_path(ApplicationState* state, ImDrawList* dl);
static float dof_focus_depth(const ApplicationState* state, const ViewTransform& view);
static bool movie_focus_active(const ApplicationState* state);
static bool movie_focus_owned(const ApplicationState* state);
static void movie_params_restore(ApplicationState* state);
static bool movie_target_center(const ApplicationState* state, const md_bitfield_t* mask, const md_system_state_t& system_state, vec3_t* out);
static bool movie_key_look_at_atom(ApplicationState* state, int key_idx, int32_t atom);
static bool movie_key_focus_on_atoms(ApplicationState* state, int key_idx, const std::vector<uint32_t>& atoms);
static bool movie_draw_timeline_markers(ApplicationState* state);
static void movie_blit_preview(ApplicationState* state);
static void movie_capture_frame(ApplicationState* state);
static void movie_overlays_draw(ImDrawList* dl, ImVec2 pos, ImVec2 size, double time, const ApplicationState* state);
static void movie_sort_keyframes(ApplicationState* state);
static bool movie_frame_guide(const ApplicationState* state, ImVec2* pos, ImVec2* size);
static bool movie_frame_guide_shift(const ApplicationState* state, float* sx, float* sy);
static void movie_frame_size(const ApplicationState* state, int* w, int* h);
static void movie_property_vis_apply(ApplicationState* state);
static void script_vis_text_draw(ImDrawList* dl, ImVec2 res, float scale, const ApplicationState& state);
static void movie_apply_time(ApplicationState* state, double time, bool apply_camera);
static void movie_follow_update(ApplicationState* state);
static double movie_duration(const ApplicationState* state);
static double movie_trajectory_frame(const ApplicationState* state, double time);
static int movie_num_frames(const ApplicationState* state);
static void movie_restore_state(ApplicationState* state);
static void movie_set_scene_view(ApplicationState* state, bool scene);
static void movie_play_mode(ApplicationState* state, bool on);
static void draw_movie_play_bar(ApplicationState* state);

// The sizes offered in the Settings menu. The stored setting is the size itself, not an
// index into this table, so the table can change without invalidating anyone's .ini.
static const float font_sizes[] = { 10.0f, 12.0f, 14.0f, 16.0f, 18.0f, 20.0f, 24.0f, 30.0f, 36.0f, 48.0f, 64.0f, 72.0f };
static const char* font_size_names[] = { "10", "12", "14", "16", "18", "20", "24", "30", "36", "48", "64", "72" };

static int nearest_font_size_index(float size) {
    int   best_idx  = 0;
    float best_dist = FLT_MAX;
    for (int i = 0; i < (int)ARRAY_SIZE(font_sizes); ++i) {
        const float dist = fabsf(font_sizes[i] - size);
        if (dist < best_dist) {
            best_dist = dist;
            best_idx  = i;
        }
    }
    return best_idx;
}

// Also the app_settings apply hook, so a size read from the .ini takes effect on the first frame.
static void apply_font_size(void* user_data) {
    ApplicationState* state = (ApplicationState*)user_data;
    ImGuiStyle& style = ImGui::GetStyle();
    style.FontSizeBase = state->settings.font_size;
    style._NextFrameFontSizeBase = style.FontSizeBase;  // From the demo, seems like a temporary fix
}

// timeline.x_values and view_range are stored in display units, so a change to the time preference
// has to be carried into them. The raw frame times are still with the trajectory, so this re-derives
// rather than compounding factors, and the visible range is moved by the ratio so the user keeps
// looking at the same stretch of trajectory.
static void update_timeline_time_unit(ApplicationState* state) {
    ASSERT(state);

    if (state->timeline.units_version == display_units::version()) {
        return;
    }
    state->timeline.units_version = display_units::version();

    const double prev_scl = state->timeline.time_scale;
    const double time_scl = display_units::factor(&state->timeline.time_unit, run_time_unit(state));
    state->timeline.time_scale = time_scl;

    if (time_scl == prev_scl) {
        return;
    }

    const double* frame_times = run_frame_times(state);
    const size_t  num_frames  = md_array_size(state->timeline.x_values);
    if (frame_times) {
        for (size_t i = 0; i < num_frames; ++i) {
            state->timeline.x_values[i] = (float)(frame_times[i] * time_scl);
        }
    }

    const double ratio = prev_scl != 0.0 ? time_scl / prev_scl : 1.0;
    state->timeline.view_range.beg_x *= ratio;
    state->timeline.view_range.end_x *= ratio;
}

static void modify_selection(ApplicationState* state, md_bitfield_t* atom_mask, SelectionOperator op = SelectionOperator::Set) {
    ASSERT(state);
    modify_field(&state->selection.selection_mask, atom_mask, op);
}

int main(int argc, char** argv) {
#if DEBUG
    persistent_alloc = md_tracking_allocator_create(md_get_heap_allocator());
#elif RELEASE
    persistent_alloc = md_get_heap_allocator();
#else
#error "Must define DEBUG or RELEASE"
#endif
    frame_alloc = md_vm_arena_create(GIGABYTES(4));

    struct NotificationState {
        md_mutex_t lock;
        uint64_t hash;
        md_tick_t time;
    };

    NotificationState notification_state = {
        .lock = md_mutex_create(),
        .hash = 0,
        .time = 0
    };

    md_logger_i notification_logger = {
        (md_logger_o*)&notification_state,
        [](struct md_logger_o* inst, enum md_log_type_t log_type, const char* msg) {
            NotificationState& state = *(NotificationState*)inst;            

            // Prevent spamming the logger with the same message by comparing its hash
            const md_tick_t time = md_tick_now();
            const uint64_t hash = md_hash64(msg, strlen(msg), 0);

            if (md_tick_to_seconds(time - state.time) < 1.0 && hash == state.hash) {
                return;
            }
            state.hash = hash;
            state.time = time;
            
            ImGuiToastType toast_type = ImGuiToastType_None;
            switch (log_type) {
            case MD_LOG_TYPE_INFO:
                toast_type = ImGuiToastType_Info;
                break;
            case MD_LOG_TYPE_ERROR:
                toast_type = ImGuiToastType_Error;
                break;
            case MD_LOG_TYPE_DEBUG:
            default:
                break;
            }
            if (toast_type != ImGuiToastType_None) {
                // @NOTE: This needs to be protected with a mutex as it pushes to an internal vector
                md_mutex_lock(&state.lock);
                ImGui::InsertNotification(ImGuiToast(toast_type, (uint64_t)(NOTIFICATION_DISPLAY_TIME_IN_SECONDS * 1000), msg));
                md_mutex_unlock(&state.lock);
            }
        }
    };

    md_log_register(&notification_logger);

    ApplicationState state;
    ViamdEventHandler event_handler(&state);

    state.allocator.persistent = persistent_alloc;
    state.allocator.frame = frame_alloc;
    state.file_queue.ring = md_ring_allocator_create(md_alloc(persistent_alloc, MEGABYTES(1)), MEGABYTES(1));
    // One arena per dataset, held by the system and the state that share it. Both are rewound
    // rather than destroyed between loads, so this handle is set once and never reassigned.
    state.mold.sys.alloc   = md_arena_allocator_create(persistent_alloc, MEGABYTES(1));
    state.mold.state.alloc = state.mold.sys.alloc;

    md_temp_arena_system_init();

    md_bitfield_init(&state.selection.selection_mask, persistent_alloc);
    md_bitfield_init(&state.selection.highlight_mask, persistent_alloc);
    md_bitfield_init(&state.selection.query.mask, persistent_alloc);
    md_bitfield_init(&state.selection.grow.mask, persistent_alloc);
    md_bitfield_init(&state.operations.selection_mask, persistent_alloc);
    md_bitfield_init(&state.movie.follow_mask, persistent_alloc);
    md_bitfield_init(&state.visuals.dof.target_mask, persistent_alloc);
    md_bitfield_init(&state.operations.recenter_query.mask, persistent_alloc);
    md_bitfield_init(&state.representation.visibility_mask, persistent_alloc);

    // Init platform
    VIAMD_LOG_DEBUG("Initializing GL...");
    if (!application::initialize(&state.app, 0, 0, STR_LIT("VIAMD"))) {
        VIAMD_LOG_ERROR("Could not initialize application...\n");
        return -1;
    }

    // Application settings live in the ImGui .ini. Bind first, then initialize: that call
    // reads the file, so the bound values are current by the time the loop starts.
    app_settings::bind(STR_LIT("keep_representations"), &state.settings.keep_representations);
    app_settings::bind(STR_LIT("exact_isosurfaces"), &state.settings.exact_isosurfaces);
    app_settings::bind(STR_LIT("font_size"), &state.settings.font_size);
    app_settings::on_apply(apply_font_size, &state);
    display_units::register_settings();
    app_settings::initialize();

#if MD_ENABLE_GPU
    VIAMD_LOG_DEBUG("Initializing GPU device...");
    state.gpu_device = md_gpu_device_create(nullptr);
    if (!state.gpu_device) {
        const char* reason = md_gpu_last_error();
        VIAMD_LOG_ERROR("Failed to create GPU device: %s", reason ? reason : "unknown");
    } else {
        // The scratch belongs to the device, not to whichever component happens to evaluate
        // first. Components borrow it; see the note on ApplicationState.
        state.gpu_stream = md_gpu_stream_default(state.gpu_device, MD_GPU_STREAM_COMPUTE);

        md_gpu_texture_desc_t vol_desc = {
            .type            = MD_GPU_TEX_3D,
            .format          = MD_GPU_FORMAT_R32_FLOAT,
            .usage           = MD_GPU_TEX_STORAGE,
            .width           = GPU_VOLUME_DIM,
            .height          = GPU_VOLUME_DIM,
            .depth_or_layers = GPU_VOLUME_DIM,
            .label           = "Evaluation volume",
        };
        state.gpu_volume = md_gpu_texture_create(state.gpu_stream, &vol_desc);
        if (!state.gpu_volume) {
            VIAMD_LOG_ERROR("Failed to create the GPU evaluation volume: %s", md_gpu_last_error());
        }

        // Potentials (surface fields). Kernels are made on first use, so this costs nothing until a
        // surface is coloured by one.
        md_gto_int_gpu_initialize(state.gpu_device);
    }
#endif

    state.app.window.vsync = true;
    state.app.file_drop.user_data = &state;
    state.app.file_drop.callback = [](size_t num_files, const str_t file_paths[], void* user_data) {
        ApplicationState* state = (ApplicationState*)user_data;
        ASSERT(state);

        for (size_t i = 0; i < num_files; ++i) {
            file_queue_push(&state->file_queue, file_paths[i], FileFlags_None);
        }
    };

    VIAMD_LOG_DEBUG("Initializing framebuffer...");
    gbuffer_init(&state.gbuffer, state.app.framebuffer.width, state.app.framebuffer.height);
    picking_surface_init(&state.picking_surface, interaction_surface_main);

    for (int i = 0; i < (int)ARRAY_SIZE(state.view.jitter.sequence); ++i) {
        state.view.jitter.sequence[i].x = md_halton(i + 1, 2);
        state.view.jitter.sequence[i].y = md_halton(i + 1, 3);
    }



    // Init subsystems
    VIAMD_LOG_DEBUG("Initializing immediate draw...");
    immediate::initialize();
    state.gfx.world   = immediate::queue_create("world");
	state.gfx.overlay = immediate::queue_create("overlay");

    VIAMD_LOG_DEBUG("Initializing post processing...");

    postprocess_pipeline::initialize(state.gbuffer.width, state.gbuffer.height);
    VIAMD_LOG_DEBUG("Initializing volume...");
    volume::initialize();
    VIAMD_LOG_DEBUG("Initializing task system...");
    md_os_sys_info_t sys_info = {0};
	md_os_sys_info_query(&sys_info);

    int num_threads = VIAMD_NUM_WORKER_THREADS == 0 ? sys_info.num_physical_cores : VIAMD_NUM_WORKER_THREADS;
	num_threads = CLAMP(num_threads, 2, sys_info.num_physical_cores);
    task_system::initialize((size_t)num_threads);

    md_gl_initialize();
    state.gl.shaders                = md_gl_shaders_create(shader_output_snippet);
    state.gl.shaders_lean_and_mean  = md_gl_shaders_create(shader_output_snippet_lean_and_mean);

    // Which of these are open is part of a workspace; the components register their own on initialize
    workspace_register_window("Timelines",       &state.timeline.show_window);
    workspace_register_window("Distributions",   &state.distributions.show_window);
    workspace_register_window("Representations", &state.representation.show_window);
    workspace_register_window("ScriptEditor",    &state.show_script_window);
    workspace_register_window("Animation",       &state.animation.show_window);
    workspace_register_window("Movie",           &state.movie.show_window);

    viamd::event_system_broadcast_event(viamd::EventType_ViamdInitialize, viamd::EventPayloadType_ApplicationState, &state);

#if EXPERIMENTAL_GFX_API
    md_gfx_initialize(state.gbuffer.width, state.gbuffer.height, 0);
#endif

    ImGui::init_theme();

    state.editor.SetLanguage(script_editor::language());
    state.editor.SetPalette(TextEditor::GetDarkPalette());
    state.editor.SetInsertSpacesOnTabs(true);
    script_editor::enable_autocomplete(state.editor, &state.mold.sys);

    {
#ifdef VIAMD_DEFAULT_DATASET
        {
            md_strb_t sb = md_strb_create(frame_alloc);
            str_t path = viamd::resource_path(&sb, STR_LIT(VIAMD_DEFAULT_DATASET));
            if (md_path_is_valid(path)) {
                // @NOTE: We want explicitly to disable writing of cache files for the default dataset
                // The motivation is that the dataset may reside in a shared folder on the system that has no write access.
                file_queue_push(&state.file_queue, path, FileFlags_DisableCacheWrite);
                state.editor.SetText("s1 = resname(\"ALA\")[2:8];\nd1 = distance(10,30);\na1 = angle(2,1,3) in resname(\"ALA\");\nr = rdf(element('C'), element('H'), 10.0);\nv = sdf(s1, element('H'), 10.0);\n{lin,plan,iso} = shape_weights(all);");
            }
        }
#endif
        if (argc > 1) {
            // Assume argv[1..] are files to load
            // Currently we do not support any command line flags
            // So anything here which is a file path is assumed to be a file to load
            for (int i = 1; i < argc; ++i) {
                str_t path = str_from_cstr(argv[i]);
                if (md_path_is_valid(path)) {
                    file_queue_push(&state.file_queue, path);
                }
            }
        }
    }

#if EXPERIMENTAL_SDF == 1
    draw::scan::test_scan();
#endif
    bool time_changed = true;
    bool time_stopped = true;

    //bool demo_window = true;

    // Main loop
    while (!state.app.window.should_close) {
        application::update(&state.app);
        
        // This needs to happen first (in imgui events) to enable docking of imgui windows
#if VIAMD_IMGUI_ENABLE_DOCKSPACE
        ImGui::CreateDockspace();
#endif

        const size_t num_frames  = run_num_frames(&state);
        const size_t last_frame  = num_frames > 0 ? num_frames - 1 : 0;
        const double   max_frame = (double)last_frame;

        file_queue_process(&state);

        state.movie.vis_fade = 1.0f;
        state.script.vis = {0};
        md_script_vis_init(&state.script.vis, state.allocator.frame);

        picking_handler_new_frame(&state.picking_handler);
        viamd::event_system_broadcast_event(viamd::EventType_ViamdPickingRangeReserve, viamd::EventPayloadType_PickingSpace, picking_handler_current_space(&state.picking_handler));

#if MD_ENABLE_GPU
        // Retires stream-ordered frees and fires completion callbacks for work
        // that finished since the last frame. Callbacks run on the calling
        // thread, so this must stay on the GL thread: they issue GL calls.
        // Not optional -- without it, freed device memory is never reclaimed.
        if (state.gpu_device) {
            md_gpu_device_poll(state.gpu_device);
        }
#endif

        viamd::event_system_broadcast_event(viamd::EventType_ViamdFrameTick, viamd::EventPayloadType_ApplicationState, &state);

        // GUI. In the movie's Preview only the bar that plays it is shown.
        if (state.movie.play_mode && (!state.movie.show_window || state.movie.state == MovieRecordingState::Recording)) movie_play_mode(&state, false);
        const bool play_mode = state.movie.play_mode;
        if (play_mode) draw_movie_play_bar(&state);
        const bool focus_owned = movie_focus_owned(&state);
        if (focus_owned && !state.movie.editor_focus_active) {
            state.movie.editor_focus_distance = state.view.camera.distance;
            state.movie.editor_saved_aperture = state.visuals.dof.aperture;
            state.visuals.dof.aperture = 0.0f;
        } else if (!focus_owned && state.movie.editor_focus_active) {
            movie_params_restore(&state);
            state.visuals.dof.aperture = state.movie.editor_saved_aperture;
        }
        state.movie.editor_focus_active = focus_owned;
        if (!play_mode) {
        if (state.show_script_window) draw_script_editor_window(&state);
        if (state.show_script_reference_window) draw_script_reference_window(&state);
        if (state.load_dataset.show_window) draw_load_dataset_window(&state);
        if (state.representation.show_window) draw_representations_window(&state);
        if (state.distributions.show_window) draw_distribution_window(&state);
        if (state.timeline.show_window) draw_timeline_window(&state);
        series_cache_gc(&state);
        if (state.selection.query.show_window) draw_selection_query_window(&state);
        if (state.selection.grow.show_window) draw_selection_grow_window(&state);
        if (state.show_property_export_window) draw_property_export_window(&state);
        if (state.structure_export.show_window) draw_structure_export_window(&state);
        if (state.show_debug_window) draw_debug_window(&state);
        if (state.animation.show_window) draw_animation_window(&state);
        if (state.movie.show_window) draw_movie_window(&state);
        }
        draw_movie_recording_banner(&state);

        draw_async_task_window(&state);
        if (!play_mode) draw_main_menu(&state);
        draw_notifications_window();

        //ImGui::ShowDemoWindow();

        if (!play_mode) draw_coordinate_system_widget_window(&state.view.target, state.view.camera);
        if (!play_mode) draw_color_legend_windows(state);
            
        ImGui::BeginCanvas("Main interaction window", true);
        ImVec2 view_size = ImGui::GetContentRegionAvail();
        InteractionSurfaceState surface_state = interaction_surface(interaction_surface_main, vec_cast(view_size));
        ImGui::EndCanvas();

        // The recording drives the frame and the camera, and the viewport shows the frames at the movie's
        // resolution rather than the scene under the mouse, so picking and navigation are off meanwhile
        const bool movie_recording = state.movie.state == MovieRecordingState::Recording;

        const mat4_t clip_to_world = camera_view_to_world_matrix(state.view.camera) * state.view.param.matrix.inv.proj;
        const mat4_t world_to_clip = state.view.param.matrix.curr.proj * state.view.param.matrix.curr.view;

        PickingHit hit = {};

        // A keyframe waits for an atom: the click picks it, so it must not select or rotate as well
        const bool look_pick = (state.movie.look_pick_key >= 0 || state.movie.focus_pick_key >= 0) && !movie_recording;
        // A handle of the camera path is hovered or dragged: the viewport must not select or rotate as well
        const bool surface_blocked = look_pick || ((state.movie.path_hot || state.movie.path_drag.active) && !movie_recording);
        if (look_pick && ImGui::IsKeyPressed(ImGuiKey_Escape)) state.movie.look_pick_key = state.movie.focus_pick_key = -1;

        if (surface_state.hovered && !movie_recording) {
            InteractionSurfaceHitArgs args = {
                .picking_surface = &state.picking_surface,
                .picking_handler = state.picking_handler,
                .fbo = state.gbuffer.fbo,
                .width = state.gbuffer.width,
                .height = state.gbuffer.height,
                .clip_to_world = clip_to_world,
            };

            interaction_surface_hit_extract(&hit, surface_state, args);

            if (look_pick) {
                ImGui::SetMouseCursor(ImGuiMouseCursor_Hand);
                if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && hit.domain == PickingDomain_Atom) {
                    if (state.movie.look_pick_key >= 0) movie_key_look_at_atom(&state, state.movie.look_pick_key, (int32_t)hit.local_idx);
                    else movie_key_focus_on_atoms(&state, state.movie.focus_pick_key, {(uint32_t)hit.local_idx});
                    state.movie.look_pick_key = state.movie.focus_pick_key = -1;
                }
            }

            InteractionSurfaceEvent event = {};
            interaction_surface_event_extract(&event, surface_state, hit);

            event.clip_to_world = clip_to_world;
            event.world_to_clip = world_to_clip;

            if (surface_blocked) {
                // handled above, or by the path
            } else if (event.kind == InteractionSurfaceEventKind::RegionSelect) {
                const md_bitfield_t* candidate_mask = &state.representation.visibility_mask;
                if (event.selection_mode == InteractionSelectionMode::Remove) {
                    // When removing, only consider currently selected atoms as candidates for region selection
                    candidate_mask = &state.selection.selection_mask;
                }
                point_set_region_mask_compute(&state.selection.highlight_mask,
                    state.mold.state.xyz, state.mold.state.num_atoms,
                    candidate_mask, world_to_clip, event.region_min, event.region_max, event.surface_size);

                grow_mask_by_selection_granularity(&state.selection.highlight_mask, state.selection.granularity, state.mold.sys);
                if (event.region_phase == InteractionSurfaceEventPhase::Commit) {
                    // Merge highlight into selection
                    if (event.selection_mode == InteractionSelectionMode::Append) {
                        md_bitfield_or_inplace(&state.selection.selection_mask, &state.selection.highlight_mask);
                    }
                    else if (event.selection_mode == InteractionSelectionMode::Remove) {
                        md_bitfield_andnot_inplace(&state.selection.selection_mask, &state.selection.highlight_mask);
                    }
                    md_bitfield_clear(&state.selection.highlight_mask);
                    // A region select is not a sequence of individual picks, so whatever order was
                    // recorded before no longer describes the selection.
                    single_selection_sequence_clear(&state.selection.single_selection_sequence);
                }
            } else if (event.kind == InteractionSurfaceEventKind::ContextMenu) {
                ImGui::OpenPopup("Context Popup");
            }
            
            // Since this is the main interaction view, we still broadcast all kinds of events, even if we have handled some explicitly here.
            if (!surface_blocked && event.kind != InteractionSurfaceEventKind::None) {
                viamd::event_system_broadcast_event(viamd::EventType_ViamdInteractionSurface, viamd::EventPayloadType_InteractionSurfaceEvent, &event);
            }
        }

        // In Scene view of the movie, the wheel zooms towards the mouse, around the middle of the frame. The viewport around the
        // frame has a wider field of view than the camera's, which is the frame's.
        float view_sx = 0.0f, view_sy = 0.0f;
        const bool scene_view = state.movie.show_window && !state.movie.show_frame && !state.movie.play_mode;
        if (scene_view) movie_frame_guide_shift(&state, &view_sx, &view_sy);
        Camera view_camera = state.view.camera;
        {
            ImVec2 gp, gs;
            if (movie_frame_guide(&state, &gp, &gs) && gs.y > 0.0f) view_camera.fov_y = movie_guide_fov_y(view_camera.fov_y, (float)state.app.window.height, gs.y);
        }
        InteractionSurfaceViewTransformArgs view_args = {
            .camera = view_camera,
            .trackball_param = state.view.trackball_param,
            .zoom_to_cursor = scene_view,
            .center_offset = {view_sx * 0.5f * surface_state.surface_size.x, -view_sy * 0.5f * surface_state.surface_size.y},
        };

        InteractionSurfaceViewTransformResult view_result = {};
        if (!movie_recording && !surface_blocked) {
            view_result = interaction_surface_view_transform_apply(&state.view.target, surface_state, view_args);
        }
        if (view_result.reset_requested) {
            ViewTransform reset_transform = {};
            if (hit.depth < 1.0f) {
                reset_transform.distance = state.view.target.distance;
                reset_transform.orientation = state.view.camera.orientation;
                reset_transform.position = hit.world_pos + state.view.camera.orientation * vec3_set(0, 0, state.view.target.distance);
            } else {
                reset_view(&reset_transform, state.mold.state, &state.representation.visibility_mask);
            }
            state.view.target = reset_transform;
        }

        draw_context_popup(&state, hit);

        camera_animate(&state.view.camera, state.view.target, state.app.timing.delta_s);

        ImGuiWindow* win = ImGui::GetCurrentContext()->HoveredWindow;
        if (win && strcmp(win->Name, "Main interaction window") == 0) {
            script_set_hovered_property(&state,  STR_LIT(""));
        }

        // Tab held for a moment: back to what it was before when it is let go
        if (state.movie.peek_active && !ImGui::IsKeyDown(ImGuiKey_Tab)) {
            state.movie.peek_active = false;
            if (ImGui::GetTime() - state.movie.peek_t0 > 0.3 && state.movie.show_window && !movie_recording) movie_set_scene_view(&state, state.movie.show_frame);
        }

        // Capture non-window specific keyboard events
        if (state.movie.show_window && !state.movie.play_mode && !movie_recording &&
            !ImGui::GetIO().WantTextInput) {
            if (ImGui::IsKeyChordPressed(ImGuiMod_Ctrl | ImGuiKey_Z)) {
                movie_undo(&state);
            } else if (ImGui::IsKeyChordPressed(ImGuiMod_Ctrl | ImGuiKey_Y) ||
                       ImGui::IsKeyChordPressed(ImGuiMod_Ctrl | ImGuiMod_Shift | ImGuiKey_Z)) {
                movie_redo(&state);
            }
        }
        if (state.movie.show_window && !state.movie.play_mode && !movie_recording &&
            !ImGui::GetIO().WantTextInput && ImGui::IsKeyPressed(ImGuiKey_Tab, false)) {
            movie_set_scene_view(&state, state.movie.show_frame);
            state.movie.peek_active = true;
            state.movie.peek_t0 = ImGui::GetTime();
        }
        if ((state.movie.show_window || state.movie.play_mode) && !movie_recording &&
            !ImGui::GetIO().WantTextInput && ImGui::GetIO().KeyShift &&
            ImGui::IsKeyPressed(KEY_ADD_MOVIE_KEYFRAME, false)) {
            movie_add_keyframe_from_view(&state);
        }
        if (!ImGui::GetIO().WantCaptureKeyboard) {
            if (ImGui::IsKeyPressed(KEY_SHOW_DEBUG_WINDOW)) {
                state.show_debug_window = true;
            }

            if (movie_recording && ImGui::IsKeyPressed(ImGuiKey_Escape, false)) {
                movie_recording_stop(&state);
            }

            // With the Movie window up, Space plays the movie (from the preview time), not the trajectory. Opening the window stops
            // the trajectory's playback (the movie's time drives the frame then); closing it stops the movie's.
            // Not while typing, nor while the keyboard focus outline is on a widget (Space presses that widget then).
            {
                static bool movie_window_was_open = false;
                if (state.movie.show_window != movie_window_was_open) {
                    if (state.movie.show_window) {
                        if (state.animation.mode == PlaybackMode::Playing) {
                            state.animation.mode = PlaybackMode::Stopped;
                            state.mold.dirty_gpu_buffers |= MolBit_ClearVelocity;
                        }
                    } else if (!state.movie.play_mode) {
                        state.movie.preview_playing = false;
                    }
                    movie_window_was_open = state.movie.show_window;
                }
            }
            if ((state.movie.play_mode || (state.movie.show_window && !movie_recording)) && !ImGui::GetIO().WantTextInput && !ImGui::GetIO().NavVisible && ImGui::IsKeyPressed(ImGuiKey_Space, false)) {
                state.movie.preview_playing = !state.movie.preview_playing;
                if (state.movie.preview_playing && state.movie.playhead >= (float)movie_duration(&state)) state.movie.playhead = 0.0f;
            }
            if (state.movie.play_mode) {
                if (ImGui::IsKeyPressed(ImGuiKey_Escape, false)) movie_play_mode(&state, false);
            }

            if (state.movie.show_window && !movie_recording && !state.movie.play_mode && !ImGui::GetIO().KeyShift && ImGui::IsKeyPressed(KEY_ADD_MOVIE_KEYFRAME, false)) {
                movie_add_keyframe_from_view(&state);
            }

            if (state.movie.show_window && !movie_recording && !state.movie.play_mode) {
                if (state.movie.shortcut_frame == ImGui::GetFrameCount()) {
                    // the lanes took a shortcut this frame
                } else if (ImGui::IsKeyChordPressed(ImGuiMod_Ctrl | ImGuiKey_C)) {
                    movie_copy_keyframe_at_playhead(&state);
                } else if (ImGui::IsKeyChordPressed(ImGuiMod_Ctrl | ImGuiKey_V)) {
                    movie_paste_keyframe(&state);
                }
            }

            if (ImGui::IsKeyPressed(KEY_RECOMPILE_SHADERS)) {
                VIAMD_LOG_INFO("Recompiling shaders and re-initializing volume");
                postprocess_pipeline::initialize(state.gbuffer.width, state.gbuffer.height);
                volume::initialize();
                md_gl_shaders_destroy(state.gl.shaders);
                state.gl.shaders = md_gl_shaders_create(shader_output_snippet);
            }

            if (!movie_recording && !state.movie.play_mode && !state.movie.show_window && ImGui::IsKeyPressed(KEY_PLAY_PAUSE)) {
                switch (state.animation.mode) {
                    case PlaybackMode::Playing:
                        state.animation.mode = PlaybackMode::Stopped;
                        state.mold.dirty_gpu_buffers |= MolBit_ClearVelocity;
                        break;
                    case PlaybackMode::Stopped:
                        state.animation.mode = PlaybackMode::Playing;
                        if (state.animation.frame == max_frame && state.animation.fps > 0) {
                            state.animation.frame = 0;
                        } else if (state.animation.frame == 0 && state.animation.fps < 0) {
                            state.animation.frame = max_frame;
                        }
                        break;
                    default:
                        ASSERT(false);
                }

            }

            if (!movie_recording && !state.movie.play_mode && state.movie.shortcut_frame != ImGui::GetFrameCount() && (ImGui::IsKeyPressed(KEY_SKIP_TO_PREV_FRAME) || ImGui::IsKeyPressed(KEY_SKIP_TO_NEXT_FRAME))) {
                double step = ImGui::IsKeyDown(ImGuiMod_Ctrl) ? 10.0 : 1.0;
                if (ImGui::IsKeyPressed(KEY_SKIP_TO_PREV_FRAME)) step = -step;
                state.animation.frame = CLAMP(state.animation.frame + step, 0.0, max_frame);
            }
        }

        if (state.animation.mode == PlaybackMode::Playing) {
            state.animation.frame += state.app.timing.delta_s * state.animation.fps;
            state.animation.frame = CLAMP(state.animation.frame, 0.0, max_frame);
            if (state.animation.frame >= max_frame) {
                state.animation.mode = PlaybackMode::Stopped;
                state.animation.frame = max_frame;
            } else if (state.animation.frame <= 0) {
                state.animation.mode = PlaybackMode::Stopped;
                state.animation.frame = 0;
            }
        }

        // Must run before the time_changed check below so a frame set here is interpolated
        // and uploaded to the GPU within this same loop iteration, before render() is called.
        update_movie_preview(&state);
        update_movie_recording(&state);
        update_movie_history(&state);

        {
            static auto prev_frame = state.animation.frame;
            if (state.animation.frame != prev_frame) {
                time_changed = true;
                prev_frame = state.animation.frame;
            }
            else {
                time_changed = false;
            }
        }

        if (state.timeline.filter.temporal_window.enabled) {
            const double pre_beg = state.timeline.filter.beg_frame;
            const double pre_end = state.timeline.filter.end_frame;
            const double half_window_ext = state.timeline.filter.temporal_window.extent_in_frames * 0.5;
            state.timeline.filter.beg_frame = CLAMP(round(state.animation.frame - half_window_ext), 0.0, max_frame);
            state.timeline.filter.end_frame = CLAMP(round(state.animation.frame + half_window_ext), 0.0, max_frame);
            if (state.script.ir && (state.timeline.filter.beg_frame != pre_beg || state.timeline.filter.end_frame != pre_end)) {
                state.script.evaluate_filt = true;
            }
        }

        if (state.timeline.filter.enabled) {
            static auto prev_filter_beg = state.timeline.filter.beg_frame;
            static auto prev_filter_end = state.timeline.filter.end_frame;
            if (state.timeline.filter.beg_frame != prev_filter_beg || state.timeline.filter.end_frame != prev_filter_end) {
                prev_filter_beg = state.timeline.filter.beg_frame;
                prev_filter_end = state.timeline.filter.end_frame;
                state.timeline.filter.fingerprint = generate_fingerprint();
            }
        }

        if (time_changed) {
            state.mold.interpolate_system_state = true;
            time_stopped = false;

            PUSH_CPU_SECTION("Flag dynamic representations for update")
            for (size_t i = 0; i < md_array_size(state.representation.reps); ++i) {
                auto& rep = state.representation.reps[i];
                if (rep.enabled && (rep.dynamic_evaluation || rep.color_mapping == ColorMapping::SecondaryStructure)) {
					flag_representation_as_dirty(&rep);
                }
            }
            POP_CPU_SECTION()
        } else {
            if (!time_stopped) {
                time_stopped = true;
                state.mold.dirty_gpu_buffers |= MolBit_DirtyPosition;
            }
        }

        if (state.mold.interpolate_system_state) {
			state.mold.interpolate_system_state = false;
            if (run_num_frames(&state) > 0) {
                PUSH_CPU_SECTION("Interpolate System State")
                interpolate_system_state(&state);
                POP_CPU_SECTION()
                viamd::event_system_enqueue_event(viamd::EventType_ViamdSystemStateChanged, viamd::EventPayloadType_ApplicationState, &state);
            }
        }

        {
            std::string text = state.editor.GetText();
            state.script.text_hash = md_hash64(text.c_str(), text.length(), 0);
        }

        if (state.script.compile_ir) {
            state.script.time_since_last_change += state.app.timing.delta_s;

            script_editor::markers_clear(&state.editor_markers);

            if (state.script.time_since_last_change > COMPILATION_TIME_DELAY_IN_SECONDS) {
                // We cannot recompile while it is evaluating.
                // Need to interrupt and wait for tasks to finish.
                if (state.script.full_eval) md_script_eval_interrupt(state.script.full_eval);
                if (state.script.filt_eval) md_script_eval_interrupt(state.script.filt_eval);

                // Try aquire all semaphores
                {
                    defer {
                        flag_all_representations_as_dirty(&state);
                    };

                    // Now we hold all semaphores for the script
                    state.script.compile_ir = false;
                    state.script.time_since_last_change = 0;

                    if (state.script.ir && state.script.ir != state.script.eval_ir) {
                        md_script_ir_free(state.script.ir);
                    }
                    
                    state.script.ir = md_script_ir_create(persistent_alloc);

                    std::string src = state.editor.GetText();
                    str_t src_str {src.data(), src.length()};

                    char buf[1024];
                    size_t len = md_path_write_cwd(buf, sizeof(buf));
                    str_t old_cwd = {buf, len};
                    defer {
                        md_path_set_cwd(old_cwd);
                    };
                    
                    str_t cwd = {};
                    if (state.files.workspace[0] != '\0') {
                        extract_folder_path(&cwd, str_from_cstr(state.files.workspace));
                    } else if (state.files.trajectory[0] != '\0') {
                        extract_folder_path(&cwd, str_from_cstr(state.files.trajectory));
                    } else if (state.files.molecule[0] != '\0') {
                        extract_folder_path(&cwd, str_from_cstr(state.files.molecule));
                    }
                    if (!str_empty(cwd)) {
                        md_path_set_cwd(cwd);
                    }

                    const size_t num_stored_selections = md_array_size(state.selection.stored_selections);                       
                    for (size_t i = 0; i < num_stored_selections; ++i) {
                        str_t name = str_from_cstr(state.selection.stored_selections[i].name);
                        const md_bitfield_t* bf = &state.selection.stored_selections[i].atom_mask;
                        md_script_ir_add_identifier_bitfield(state.script.ir, name, bf);
                    }
                    
                    if (src_str) {
                        md_script_ir_compile_from_source(state.script.ir, src_str, &state.mold.sys, NULL);

                        // The markers copy what they need from the IR (it is freed below if it did not compile)
                        script_editor::Markers* markers = &state.editor_markers;
                        script_editor::markers_set_source(markers, src_str);

                        const size_t num_errors = md_script_ir_num_errors(state.script.ir);
                        const md_log_token_t* errors = md_script_ir_errors(state.script.ir);
                        for (size_t i = 0; i < num_errors; ++i) {
                            // Errors win over everything they overlap
                            script_editor::markers_add(markers, script_editor::MarkerType_Error, INT32_MAX, errors[i].range, errors[i].text, errors[i].context);
                        }

                        const size_t num_warnings = md_script_ir_num_warnings(state.script.ir);
                        const md_log_token_t* warnings = md_script_ir_warnings(state.script.ir);
                        for (size_t i = 0; i < num_warnings; ++i) {
                            script_editor::markers_add(markers, script_editor::MarkerType_Warning, INT32_MAX - 1, warnings[i].range, warnings[i].text, warnings[i].context);
                        }

                        const size_t num_tokens = md_script_ir_num_vis_tokens(state.script.ir);
                        const md_script_vis_token_t* vis_tokens = md_script_ir_vis_tokens(state.script.ir);
                        for (size_t i = 0; i < num_tokens; ++i) {
                            const md_script_vis_token_t& tok = vis_tokens[i];
                            script_editor::markers_add(markers, script_editor::MarkerType_Visualization, tok.depth, tok.range, tok.text, nullptr, tok.payload);
                        }

                        if (md_script_ir_valid(state.script.ir)) {
                            uint64_t ir_figerprint = md_script_ir_fingerprint(state.script.ir);
                            if (state.script.ir_fingerprint != ir_figerprint) {
                                state.script.ir_fingerprint = ir_figerprint;
                            }
                        } else {
                            md_script_ir_free(state.script.ir);
                            state.script.ir = nullptr;
                        }
                    }
                }
            }
        }

        if (num_frames > 0) {
            if (state.script.eval_init) {
                if (task_system::task_is_running(state.tasks.evaluate_full)) md_script_eval_interrupt(state.script.full_eval);
                if (task_system::task_is_running(state.tasks.evaluate_filt)) md_script_eval_interrupt(state.script.filt_eval);
                    
                if (task_system::task_is_running(state.tasks.evaluate_full) == false &&
                    task_system::task_is_running(state.tasks.evaluate_filt) == false) {
                    state.script.eval_init = false;

                    if (state.script.full_eval) {
                        md_script_eval_free(state.script.full_eval);
                    }
                    if (state.script.filt_eval) {
                        md_script_eval_free(state.script.filt_eval);
                    }
                
                    if (md_script_ir_valid(state.script.ir)) {
                        if (state.script.ir != state.script.eval_ir) {
                            md_script_ir_free(state.script.eval_ir);
                            state.script.eval_ir = state.script.ir;
                        }
                        state.script.full_eval = md_script_eval_create(num_frames, state.script.eval_ir, state.allocator.persistent);
                        state.script.filt_eval = md_script_eval_create(num_frames, state.script.eval_ir, state.allocator.persistent);
                    }

                    // The evaluations' tables were freed and new ones allocated, possibly in the
                    // same place: nothing derived from the old ones may be taken for the new.
                    series_cache_free(&state);

                    state.script.evaluate_filt = true;
                    state.script.evaluate_full = true;
                }
            }

            if (state.script.full_eval && state.script.evaluate_full) {
                if (task_system::task_is_running(state.tasks.evaluate_full)) {
                    md_script_eval_interrupt(state.script.full_eval);
                } else {
                    if (md_script_ir_valid(state.script.eval_ir) &&
                        md_script_eval_ir_fingerprint(state.script.full_eval) == md_script_ir_fingerprint(state.script.eval_ir))
                    {
                        state.script.evaluate_full = false;
                        md_script_eval_clear_data(state.script.full_eval);

                        if (md_script_ir_property_count(state.script.eval_ir) > 0) {
                            state.tasks.evaluate_full = task_system::create_pool_task(STR_LIT("Eval Full"), (uint32_t)num_frames, [&state](uint32_t frame_beg, uint32_t frame_end, uint32_t thread_num) {
                                (void)thread_num;
                                md_script_eval_frame_range(state.script.full_eval, state.script.eval_ir, &state.mold.sys, str_from_cstr(state.mold.run), frame_beg, frame_end);
                            });
                            
#if MEASURE_EVALUATION_TIME
                            uint64_t time = (uint64_t)md_tick_now();
                            task_system::ID time_task = task_system::create_pool_task(STR_LIT("##Time Eval Full"), [t0 = time]() {
                                uint64_t t1 = md_tick_now();
                                double s = md_tick_to_seconds(t1 - t0);
                                VIAMD_LOG_INFO("Evaluation completed in: %.3fs", s);
                            });
#endif
                            task_system::set_task_dependency(time_task, state.tasks.evaluate_full);
                            task_system::enqueue_task(state.tasks.evaluate_full);
                        }
                    }
                }
            }

            if (state.script.filt_eval && state.script.evaluate_filt && state.timeline.filter.enabled) {
                if (task_system::task_is_running(state.tasks.evaluate_filt)) {
                    md_script_eval_interrupt(state.script.filt_eval);
                } else {
                    if (md_script_ir_valid(state.script.eval_ir) &&
                        md_script_eval_ir_fingerprint(state.script.filt_eval) == md_script_ir_fingerprint(state.script.eval_ir))
                    {
                        state.script.evaluate_filt = false;
                        md_script_eval_clear_data(state.script.filt_eval);

                        if (md_script_ir_property_count(state.script.eval_ir) > 0) {
                            const uint32_t traj_frames = (uint32_t)run_num_frames(&state);
                            const uint32_t beg_frame = CLAMP((uint32_t)state.timeline.filter.beg_frame, 0, traj_frames-1);
                            const uint32_t end_frame = CLAMP((uint32_t)state.timeline.filter.end_frame + 1, beg_frame + 1, traj_frames);
                            if (beg_frame != end_frame) {
                                state.tasks.evaluate_filt = task_system::create_pool_task(STR_LIT("Eval Filt"), end_frame - beg_frame, [offset = beg_frame, &state](uint32_t beg, uint32_t end, uint32_t thread_num) {
                                    (void)thread_num;
                                    md_script_eval_frame_range(state.script.filt_eval, state.script.eval_ir, &state.mold.sys, str_from_cstr(state.mold.run), offset + beg, offset + end);
                                });
                                task_system::enqueue_task(state.tasks.evaluate_filt);
                            }
                        }
                    }
                }
            }
        }

        recenter_update(&state);

		// Perform once per-frame updates of representations (if required)
		update_all_representations(&state);

        if (state.representation.atom_visibility_mask_dirty) {
            recompute_atom_visibility_mask(&state);
            state.representation.atom_visibility_mask_dirty = false;
        }

        // The visualization of the script property that a movie overlay shows
        movie_property_vis_apply(&state);

        if (state.script.vis.text) {
            PUSH_CPU_SECTION("Draw vis text");
            ImGuiWindow* window = ImGui::FindWindowByName("Main interaction window");
            if (window) {
                // The text is as large relative to the frame as it will be in the recording
                ImVec2 guide_pos, guide_size;
                const bool guided = movie_frame_guide(&state, &guide_pos, &guide_size);
                script_vis_text_draw(window->DrawList, ImVec2((float)state.app.window.width, (float)state.app.window.height),
                    guided ? guide_size.y / (float)state.app.window.height : 1.0f, state);
            }
            POP_CPU_SECTION();
        }

        {
            ImGuiWindow* window = ImGui::FindWindowByName("Main interaction window");
            ImVec2 guide_pos, guide_size;
            const bool guided = movie_frame_guide(&state, &guide_pos, &guide_size);
            if (window && guided) {
                // What is outside the frame is dimmed, and the frame is outlined
                const ImVec2 view(0, 0), end((float)state.app.window.width, (float)state.app.window.height);
                const ImVec2 g1(guide_pos.x + guide_size.x, guide_pos.y + guide_size.y);
                // In the Preview, what is outside the frame is black, as in the movie
                const ImU32 dim = state.movie.play_mode ? IM_COL32(0, 0, 0, 255) : IM_COL32(0, 0, 0, 120);
                ImDrawList* dl = window->DrawList;
                dl->AddRectFilled(view, ImVec2(end.x, guide_pos.y), dim);
                dl->AddRectFilled(ImVec2(view.x, g1.y), end, dim);
                dl->AddRectFilled(ImVec2(view.x, guide_pos.y), ImVec2(guide_pos.x, g1.y), dim);
                dl->AddRectFilled(ImVec2(g1.x, guide_pos.y), ImVec2(end.x, g1.y), dim);
                if (!state.movie.play_mode) {
                    dl->AddRect(guide_pos, g1, IM_COL32(255, 255, 255, 170), 0.0f, 0, 1.5f);
                    int fw, fh;
                    movie_frame_size(&state, &fw, &fh);
                    char label[48];
                    snprintf(label, sizeof(label), "%d x %d", fw, fh);
                    dl->AddText(ImVec2(guide_pos.x + 4.0f, guide_pos.y - ImGui::GetFontSize() - 2.0f), IM_COL32(255, 255, 255, 190), label);
                }
            }
            if (window && state.movie.show_overlay_preview && !state.movie.overlays.empty() && state.movie.state != MovieRecordingState::Recording &&
                state.movie.show_window) {
                // The movie's overlays at the preview time, laid out in the frame of the movie (the viewport while the frame is not shown)
                const ImVec2 pos = guided ? guide_pos : ImVec2(0, 0);
                const ImVec2 size = guided ? guide_size : ImVec2((float)state.app.window.width, (float)state.app.window.height);
                movie_overlays_draw(window->DrawList, pos, size, (double)state.movie.playhead, &state);
            }
            if (window && !state.movie.play_mode) movie_draw_camera_path(&state, window->DrawList);
        }

        if (ImGui::IsKeyPressed(KEY_RECENTER_ON_HIGHLIGHT) && !state.editor_focused) {
			ViewFitRequest fit_request = {
				.app = state,
				.surface_id = interaction_surface_main,
                .xyzw = nullptr,
                .alloc = state.allocator.frame,
				.round = ViewFitRound_Highlight,
            };
			viamd::event_system_broadcast_event(viamd::EventType_ViamdViewFit, viamd::EventPayloadType_ViewFitRequest, &fit_request);

			if (fit_request.xyzw == nullptr) {
                // If no highlight, fallback to selection
                fit_request.round = ViewFitRound_Selection;
                viamd::event_system_broadcast_event(viamd::EventType_ViamdViewFit, viamd::EventPayloadType_ViewFitRequest, &fit_request);
			}
			if (fit_request.xyzw == nullptr) {
                // If no selection, fallback to all
                fit_request.round = ViewFitRound_Visible;
                viamd::event_system_broadcast_event(viamd::EventType_ViamdViewFit, viamd::EventPayloadType_ViewFitRequest, &fit_request);
            }

            size_t count = md_array_size(fit_request.xyzw);
            if (count > 0) {
                // Extract OBB from returned pointcloud
                vec3_t com = md_util_com_compute_vec4(fit_request.xyzw, 0, count, 0);
                mat3_t cov = mat3_covariance_matrix_vec4(fit_request.xyzw, 0, count, com);
                mat3_eigen_t eigen = mat3_eigen(cov);
                mat3_t PCA = mat3_orthonormalize(mat3_extract_rotation(eigen.vectors));

                // Compute min and maximum extent along the PCA axes
                mat3_t basis = mat3_transpose(PCA);
                vec3_t min_obb = vec3_set1(FLT_MAX);
                vec3_t max_obb = vec3_set1(-FLT_MAX);

                // Transform the gto (x,y,z,cutoff) into the PCA frame to find the min and max extend within it
                for (size_t i = 0; i < count; ++i) {
                    vec3_t xyz = vec3_from_vec4(fit_request.xyzw[i]);
                    xyz = mat3_mul_vec3(PCA, xyz);
                    min_obb = vec3_min(min_obb, xyz);
                    max_obb = vec3_max(max_obb, xyz);
                }

				vec3_t ext = max_obb - min_obb;
				float max_ext = MAX(ext.x, MAX(ext.y, ext.z));

                if (max_ext > 3.0f) {
                    vec3_t half_ext = (max_obb - min_obb) * 0.5f;
                    state.view.target = compute_optimal_view(com, half_ext, basis);
                } else {
                    state.view.target.position = com + state.view.target.orientation * vec3_set(0, 0, state.view.target.distance);
                }
                state.view.target.position = mat4_mul_vec3(state.mold.unitcell_transform, state.view.target.position, 1.0f);
            }
            else {
                reset_view(&state.view.target, state.mold.state);
            }
        }

        // The motivation for doing this is to reduce the frequency at which we invalidate and upload the atom flag field to the GPU
        // For large systems, this can be a costly operation: Consider a system of 100'000'000 atoms
        // The size of each bitfield to represent the mask would be 12.5 MB
        // Given a throughput of 10GB/s results in a time of 1.25 ms per bitfield.

        uint64_t v_hash = state.representation.visibility_mask_hash;
        uint64_t h_hash = md_bitfield_hash64(&state.selection.highlight_mask, 0);
        uint64_t s_hash = md_bitfield_hash64(&state.selection.selection_mask, 0);
        uint64_t f_hash = md_hash64_combine(v_hash, md_hash64_combine(h_hash, s_hash));
        
        uint64_t r_hash = 0;
        for (size_t i = 0; i < md_array_size(state.representation.reps); ++i) {
			md_hash64(&state.representation.reps[i], sizeof(Representation), r_hash);
		}

        // These represent the 'current' state so we can compare against it to see if they were modified
        static uint64_t highlight_hash = 0;
        static uint64_t selection_hash = 0;
        static uint64_t representation_hash = 0;
        static uint64_t flag_hash = 0;

        if (h_hash != highlight_hash) {
            highlight_hash = h_hash;
            // enqueue the event here and do not broadcast (stall)
            viamd::event_system_enqueue_event(viamd::EventType_ViamdHighlightMaskChanged, viamd::EventPayloadType_ApplicationState, &state);
        }

        if (s_hash != selection_hash) {
            selection_hash = s_hash;
            // enqueue the event here and do not broadcast (stall)
            viamd::event_system_enqueue_event(viamd::EventType_ViamdSelectionMaskChanged, viamd::EventPayloadType_ApplicationState, &state);
        }

        if (f_hash != flag_hash) {
            flag_hash = f_hash;
            state.mold.dirty_gpu_buffers |= MolBit_DirtyFlags;
        }

        if (r_hash != representation_hash) {
            representation_hash = r_hash;
            // enqueue the event here and do not broadcast (stall)
            viamd::event_system_enqueue_event(viamd::EventType_ViamdRepresentationChanged, viamd::EventPayloadType_ApplicationState, &state);
        }

        // Process the queued events BEFORE uploading to the GPU and rendering.
        // The ViamdSystemStateChanged handler applies the toggled operations (recenter, apply
        // pbc, unwrap) directly to mold.state coordinates. update_md_buffers uploads those
        // coordinates and then clears dirty_gpu_buffers, so anything processed after it is not
        // merely rendered a frame late - it is discarded outright, because the next frame
        // re-interpolates mold.state from the trajectory before it is ever uploaded.
        viamd::event_system_process_event_queue();

        movie_follow_update(&state);
        update_md_buffers(&state);
        update_timeline_time_unit(&state);

        render(&state);

        task_system::execute_main_task_queue();

        // Reset frame allocator
        md_vm_arena_reset(frame_alloc);

        // Swap buffers
        application::swap_buffers(&state.app);
    }

    movie_shutdown(&state);
    interrupt_async_tasks(&state);
    series_cache_free(&state);

    viamd::event_system_broadcast_event(viamd::EventType_ViamdShutdown);

    // shutdown subsystems
    VIAMD_LOG_DEBUG("Shutting down immediate draw...");
    immediate::queue_destroy(state.gfx.world);
    immediate::queue_destroy(state.gfx.overlay);
    immediate::shutdown();

    state.gfx.world = nullptr;
    state.gfx.overlay = nullptr;

    VIAMD_LOG_DEBUG("Shutting down post processing...");
    postprocess_pipeline::shutdown();
    VIAMD_LOG_DEBUG("Shutting down volume...");
    volume::shutdown();
    VIAMD_LOG_DEBUG("Shutting down task system...");
    task_system::shutdown();

    gbuffer_free(&state.gbuffer);
#if MD_ENABLE_GPU
    if (state.gpu_device) {
        // A queued readback writes into a GL texture and frees its staging memory, so nothing
        // below may go until the queue has run out.
        gpu_volume_jobs_drain(&state);
        system_gpu_data_free(&state);
        md_gto_int_gpu_shutdown();

        md_gpu_free(state.gpu_stream, state.gpu_coeff);
        md_gpu_texture_destroy(state.gpu_volume);
        state.gpu_coeff = 0;
        state.gpu_coeff_capacity = 0;
        state.gpu_volume = nullptr;
        state.gpu_stream = nullptr;

        md_gpu_device_destroy(state.gpu_device);
        state.gpu_device = nullptr;
    }
#endif
    application::shutdown(&state.app);

    return 0;
}

// #misc
static void update_view_param(ApplicationState* data) {
    ViewParam& param = data->view.param;

    // While the frame of the movie is shown in the viewport, the camera looks at it through the frame: the field of view is
    // widened so that the frame, and not the whole viewport, has the vertical field of view that the recording will have
    Camera cam = data->view.camera;
    {
        ImVec2 guide_pos, guide_size;
        if (movie_frame_guide(data, &guide_pos, &guide_size) && guide_size.y > 0.0f) {
            cam.fov_y = movie_guide_fov_y(cam.fov_y, (float)data->app.window.height, guide_size.y);
        }
    }

    param.matrix.prev = param.matrix.curr;
    param.jitter.prev = param.jitter.curr;

    param.clip_planes.near = cam.near_plane;
    param.clip_planes.far = cam.far_plane;
    param.fov_y = cam.fov_y;
    param.resolution = {(float)data->gbuffer.width, (float)data->gbuffer.height};

    param.matrix.curr.view = camera_world_to_view_matrix(cam) * data->mold.unitcell_transform;
    param.matrix.inv.view  = mat4_inverse(data->mold.unitcell_transform) * camera_view_to_world_matrix(cam);

    const float n = cam.near_plane;
    const float f = cam.far_plane;
    const float aspect_ratio = (float)data->gbuffer.width / (float)data->gbuffer.height;

    if (data->visuals.temporal_aa.enabled && data->visuals.temporal_aa.jitter && !data->movie.pip_pass) {
        static uint32_t i = 0;
        i = (i+1) % (uint32_t)ARRAY_SIZE(data->view.jitter.sequence);
        param.jitter.curr = data->view.jitter.sequence[i] - 0.5f;
        if (data->view.mode == CameraMode::Perspective) {
            const vec2_t j = param.jitter.curr;
            param.matrix.curr.proj = camera_view_to_clip_matrix_persp(cam, data->gbuffer.width, data->gbuffer.height, j.x, j.y);
            param.matrix.inv.proj  = camera_clip_to_view_matrix_persp(cam, data->gbuffer.width, data->gbuffer.height, j.x, j.y);
            param.matrix.curr.proj_no_jitter = camera_view_to_clip_matrix_persp(cam, aspect_ratio);
        } else {
            const float h = cam.distance * tanf(cam.fov_y * 0.5f);
            const float w = aspect_ratio * h;
            const vec2_t scl = {w / data->gbuffer.width * 2.0f, h / data->gbuffer.height * 2.0f};
            const vec2_t j = param.jitter.curr * scl;
            param.matrix.curr.proj = camera_view_to_clip_matrix_ortho(-w + j.x, w + j.x, -h + j.y, h + j.y, n, f);
            param.matrix.inv.proj  = camera_clip_to_view_matrix_ortho(-w + j.x, w + j.x, -h + j.y, h + j.y, n, f);
            param.matrix.curr.proj_no_jitter = camera_view_to_clip_matrix_ortho(-w, w, -h, h, n, f);

        }
    } else {
        param.jitter.curr = {0,0};
        if (data->view.mode == CameraMode::Perspective) {
            param.matrix.curr.proj = camera_view_to_clip_matrix_persp(cam, aspect_ratio);
            param.matrix.inv.proj = camera_clip_to_view_matrix_persp(cam, (float)data->gbuffer.width / (float)data->gbuffer.height);
        } else {
            const float h = cam.distance * tanf(cam.fov_y * 0.5f);
            const float w = aspect_ratio * h;
            param.matrix.curr.proj = camera_view_to_clip_matrix_ortho(-w, w, -h, h, n, f);
            param.matrix.inv.proj = camera_clip_to_view_matrix_ortho(-w, w, -h, h, n, f);
        }
        param.matrix.curr.proj_no_jitter = param.matrix.curr.proj;
    }

    // The frame of the movie off the centre of the viewport (beside the Movie window): the picture is moved along with it
    float sx, sy;
    if (movie_frame_guide_shift(data, &sx, &sy)) {
        const mat4_t T = mat4_translate(sx, sy, 0.0f);
        param.matrix.curr.proj = T * param.matrix.curr.proj;
        param.matrix.curr.proj_no_jitter = T * param.matrix.curr.proj_no_jitter;
        param.matrix.inv.proj = param.matrix.inv.proj * mat4_translate(-sx, -sy, 0.0f);
    }

    param.matrix.curr.norm = mat4_transpose(param.matrix.inv.view);
}

// ### DRAW WINDOWS ###
static void draw_main_menu(ApplicationState* data) {
    ASSERT(data);
    bool new_clicked = false;
    char path_buf[2048] = "";

    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("Open File...", "CTRL+L")) {
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Open)) {
                    file_queue_push(&data->file_queue, str_from_cstr(path_buf), FileFlags_ShowDialogue);
                }
            }
            if (ImGui::MenuItem("Open Workspace...", "CTRL+O")) {
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Open, WORKSPACE_FILE_EXTENSION)) {
                    load_workspace(data, str_from_cstr(path_buf));
                }
            }
            if (ImGui::MenuItem("Save Workspace", "CTRL+S")) {
                if (strnlen(data->files.workspace, sizeof(data->files.workspace)) == 0) {
                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, WORKSPACE_FILE_EXTENSION)) {
                        save_workspace(data, {path_buf, strnlen(path_buf, sizeof(path_buf))});
                    }
                } else {
                    save_workspace(data, str_from_cstr(data->files.workspace));
                }
            }
            if (ImGui::MenuItem("Save Workspace As...")) {
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, WORKSPACE_FILE_EXTENSION)) {
                    save_workspace(data, {path_buf, strnlen(path_buf, sizeof(path_buf))});
                }
            }
            ImGui::Separator();
            if (ImGui::MenuItem("Quit", "ALT+F4")) {
                data->app.window.should_close = true;
            }
            ImGui::EndMenu();
        }
        /*
        if (ImGui::BeginMenu("Edit")) {
        if (ImGui::MenuItem("Undo", "CTRL+Z")) {
        }
        if (ImGui::MenuItem("Redo", "CTRL+Y", false, false)) {
        }  // Disabled item
        ImGui::Separator();
        if (ImGui::MenuItem("Cut", "CTRL+X")) {
        }
        if (ImGui::MenuItem("Copy", "CTRL+C")) {
        }
        if (ImGui::MenuItem("Paste", "CTRL+V")) {
        }
        ImGui::EndMenu();
        }
        */
        if (ImGui::BeginMenu("Visuals")) {
            if (ImGui::Button("Reset View")) {
                reset_view(&data->view.target, data->mold.state, &data->representation.visibility_mask);
            }
            ImGui::Separator();
            ImGui::Checkbox("VSync", &data->app.window.vsync);
            ImGui::SetItemTooltip("Limit the frame rate to the display's refresh rate");
            ImGui::Separator();

            ImGui::BeginGroup();
            ImGui::Text("Camera");
            {
                ImGui::Combo("Mode", (int*)(&data->view.mode), "Perspective\0Orthographic\0");
                if (data->view.mode == CameraMode::Perspective) {
                    float fov = RAD_TO_DEG(data->view.camera.fov_y);
                    if (ImGui::SliderFloat("Field of View", &fov, 12.5f, 80.0f, "%.1f deg")) {
                        data->view.camera.fov_y = DEG_TO_RAD(fov);
                    }
                }
            }
            ImGui::EndGroup();

            ImGui::BeginGroup();
            ImGui::Text("Background");
            ImGui::ColorEdit3Minimal("Color", data->visuals.background.color.elem);
            ImGui::SameLine();
            ImGui::SliderFloat("##Intensity", &data->visuals.background.intensity, 0.f, 100.f);
            ImGui::SetItemTooltip("Background brightness");
            ImGui::EndGroup();
            ImGui::Separator();
            ImGui::Checkbox("Anti-Aliasing (FXAA)", &data->visuals.fxaa.enabled);
            ImGui::SetItemTooltip("Smooth jagged edges in the rendered image");
            {
                int iso_mode = data->settings.exact_isosurfaces ? 1 : 0;
                if (ImGui::Combo("Isosurfaces", &iso_mode, "Fast\0Exact\0")) {
                    data->settings.exact_isosurfaces = iso_mode == 1;
                    app_settings::mark_dirty();
                }
                ImGui::SetItemTooltip("Fast: samples the volume once per voxel along each ray. Parts of a surface thinner than\n"
                                      "about a voxel can be missed, and up close or at grazing angles the sampling can show;\n"
                                      "a finer volume resolution makes both smaller.\n"
                                      "Exact: intersects every cell a ray passes through, so nothing is missed and the surfaces\n"
                                      "do not depend on any sampling. Roughly twice the GPU time of Fast.");
            }
            // Temporal
            ImGui::BeginGroup();
            {
                ImGui::Checkbox("Temporal Anti-Aliasing", &data->visuals.temporal_aa.enabled);
                ImGui::SetItemTooltip("Smooth edges by blending each frame with the previous ones");
                if (data->visuals.temporal_aa.enabled) {
                    // ImGui::Checkbox("Jitter Samples", &data->visuals.temporal_reprojection.jitter);
                    ImGui::SliderFloat("History Min", &data->visuals.temporal_aa.feedback_min, 0.5f, 1.0f);
                    ImGui::SetItemTooltip("How much of the previous frames is kept where the image changes: lower is sharper in motion, higher is smoother");
                    ImGui::SliderFloat("History Max", &data->visuals.temporal_aa.feedback_max, 0.5f, 1.0f);
                    ImGui::SetItemTooltip("How much of the previous frames is kept where the image is still");
                    ImGui::Checkbox("Motion Blur", &data->visuals.temporal_aa.motion_blur.enabled);
                    if (data->visuals.temporal_aa.motion_blur.enabled) {
                        ImGui::SliderFloat("Motion Blur Strength", &data->visuals.temporal_aa.motion_blur.motion_scale, 0.f, 2.0f);
                    }
                    ImGui::Checkbox("Sharpen", &data->visuals.sharpen.enabled);
                    if (data->visuals.sharpen.enabled) {
                        ImGui::SliderFloat("Sharpen Strength", &data->visuals.sharpen.weight, 0.0f, 4.0f);
                    }
                }
            }
            ImGui::EndGroup();
            ImGui::Separator();

            // SSAO
            ImGui::BeginGroup();
            ImGui::PushID("SSAO");
            ImGui::Checkbox("Ambient Occlusion", &data->visuals.ssao.enabled);
            ImGui::SetItemTooltip("Darken creases and cavities, where less light reaches (SSAO)");
            if (data->visuals.ssao.enabled) {
                const char* modes[postprocess_pipeline::SsaoMode_Count] = {"Performance", "Quality"};
                int mode = (int)data->visuals.ssao.mode;
                if (ImGui::Combo("Mode", &mode, modes, IM_ARRAYSIZE(modes))) {
                    data->visuals.ssao.mode = (postprocess_pipeline::SsaoMode)mode;
                }
                ImGui::SetItemTooltip("Performance: point samples.\nQuality: screen-space slices with a visibility bitmask, finds creases and cavities more reliably and avoids halos behind small atoms.");
                const postprocess_pipeline::SsaoTimings t = postprocess_pipeline::ssao_timings();
                if (t.total_ms > 0.0f) {
                    ImGui::SameLine();
                    ImGui::TextDisabled("%.2f ms", t.total_ms);
                    ImGui::SetItemTooltip("GPU time per frame\nEstimator: %.3f ms\nTotal (with filtering): %.3f ms", t.estimator_ms, t.total_ms);
                }
                ImGui::SliderFloat("Intensity", &data->visuals.ssao.intensity, 0.0f, 8.f);
                ImGui::Checkbox("Temporal Filter", &data->visuals.ssao.temporal);
                ImGui::SetItemTooltip("Accumulate the occlusion over frames: removes the flicker in crevasses and allows a sharper result.\nOff gives a fixed sample pattern instead.");
            }
            ImGui::PopID();
            ImGui::EndGroup();
            ImGui::Separator();

#if EXPERIMENTAL_CONE_TRACED_AO == 1
            // Cone Trace
            ImGui::BeginGroup();
            ImGui::PushID("Cone Trace");
            ImGui::Checkbox("Cone Traced AO", &data->visuals.cone_traced_ao.enabled);
            if (data->visuals.cone_traced_ao.enabled) {
                ImGui::SliderFloat("Intensity", &data->visuals.cone_traced_ao.intensity, 0.01f, 5.f);
                ImGui::SliderFloat("Step Scale", &data->visuals.cone_traced_ao.step_scale, 0.25f, 8.f);
            }
            ImGui::PopID();
            ImGui::EndGroup();
            ImGui::Separator();
#endif

            // DOF
            ImGui::BeginGroup();
            const bool driven = movie_focus_owned(data);
            ImGui::BeginDisabled(driven);
            ImGui::Checkbox("Depth of Field", &data->visuals.dof.enabled);
            ImGui::EndDisabled();
            if (data->visuals.dof.enabled || driven) {
                ImGui::TextDisabled(driven ? "Focus: driven by movie. Scene view keeps its own depth." : "Focus: Auto (look-at)");
                ImGui::SetItemTooltip("Ordinary viewport focus always follows the look-at/orbit point.\nDouble-click an object to reset that point. Movie targets are edited in Focus on.");
                ImGui::TextDisabled("Focus at %.2f from the camera", dof_focus_depth(data, data->view.camera));
                ImGui::BeginDisabled(driven);
                ImGui::SliderFloat("Blur Strength", &data->visuals.dof.aperture, 0.0f, 4.0f, "%.2f %%");
                ImGui::SetItemTooltip(driven ? "Movie blur is controlled in Looks > Depth of field blur." :
                    "Blur of distant objects, in percent of the view height. Independent of zoom level");
                ImGui::EndDisabled();
            }
            ImGui::EndGroup();
            ImGui::Separator();

            // Tonemapping
            ImGui::BeginGroup();
            ImGui::Checkbox("Tone Mapping", &data->visuals.tonemapping.enabled);
            ImGui::SetItemTooltip("Map the lighting into the colors a display can show");
            if (data->visuals.tonemapping.enabled) {
                // ImGui::Combo("Function", &data->visuals.tonemapping.tonemapper, "Passthrough\0Exposure Gamma\0Filmic\0\0");
                ImGui::SliderFloat("Exposure", &data->visuals.tonemapping.exposure, 0.01f, 10.f);
                ImGui::SliderFloat("Gamma", &data->visuals.tonemapping.gamma, 1.0f, 3.0f);
            }
            ImGui::EndGroup();
            ImGui::Separator();

            ImGui::BeginGroup();
            ImGui::Checkbox("Simulation Box", &data->simulation_box.enabled);
            ImGui::SetItemTooltip("Draw the outline of the periodic box");
            if (data->simulation_box.enabled) {
                ImGui::SameLine();
                ImGui::ColorEdit4Minimal("##Box-Color", data->simulation_box.color.elem);
            }
            ImGui::EndGroup();
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Windows")) {
            ImGui::Checkbox("Animation", &data->animation.show_window);
            ImGui::Checkbox("Movie", &data->movie.show_window);
            ImGui::Checkbox("Representations", &data->representation.show_window);
            ImGui::Checkbox("Script Editor", &data->show_script_window);
            ImGui::Checkbox("Script Reference", &data->show_script_reference_window);
            ImGui::Checkbox("Timelines", &data->timeline.show_window);
            ImGui::Checkbox("Distributions", &data->distributions.show_window);
            ImGui::Checkbox("Structure Export", &data->structure_export.show_window);

            viamd::event_system_broadcast_event(viamd::EventType_ViamdWindowDrawMenu);

            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Selection")) {
            ImGui::Combo("Select by", (int*)(&data->selection.granularity), selection_granularity_str, (int)SelectionGranularity::Count);
            ImGui::SetItemTooltip("What a click or a dragged region in the viewport selects: single atoms,\nwhole components (residues) or whole instances (chains, molecules)");
            size_t num_selected_atoms = md_bitfield_popcount(&data->selection.selection_mask);
            if (ImGui::MenuItem("Invert")) {
                md_bitfield_not_inplace(&data->selection.selection_mask, 0, data->mold.sys.atom.count);
            }
            if (ImGui::IsItemHovered()) {
                md_bitfield_not(&data->selection.highlight_mask, &data->selection.selection_mask, 0, data->mold.sys.atom.count);
            }
            if (ImGui::MenuItem("Query")) data->selection.query.show_window = true;
            if (num_selected_atoms == 0) ImGui::PushDisabled();
            if (ImGui::MenuItem("Grow"))  data->selection.grow.show_window = true;
            if (num_selected_atoms == 0) ImGui::PopDisabled();
            if (ImGui::MenuItem("Clear")) {
                md_bitfield_clear(&data->selection.selection_mask);
                single_selection_sequence_clear(&data->selection.single_selection_sequence);
            }
            ImGui::Spacing();
            ImGui::Separator();

            // STORED SELECTIONS
            {
                md_bitfield_clear(&data->selection.highlight_mask);
                // @NOTE(Robin): This ImGui ItemFlag can be used to force the menu to remain open after buttons are pressed.
                // Leave it here as a comment if we feel that it is needed in the future
                //ImGui::PushItemFlag(ImGuiItemFlags_SelectableDontClosePopup, true);

                ImGui::Text("Stored Selections");
                for (int i = 0; i < (int)md_array_size(data->selection.stored_selections); i++) {
                    auto& sel = data->selection.stored_selections[i];
                    const str_t name_str = str_from_cstr(sel.name);
                    bool is_valid = md_script_identifier_name_valid(name_str);
                    char error[64] = "";
                    if (!is_valid) {
                        snprintf(error, sizeof(error), "'%s' is not a valid identifier.", sel.name);
                    }

                    for (int j = 0; j < i; ++j) {
                        if (str_eq_cstr(name_str, data->selection.stored_selections[j].name)) {
                            is_valid = false;
                            snprintf(error, sizeof(error), "identifier '%s' is already taken.", sel.name);
                            break;
                        }
                    }

                    ImGui::PushID(i);
                    ImGui::InputQuery("##label", sel.name, sizeof(sel.name), is_valid, error);
                    ImGui::SameLine();
                    if (ImGui::Button("Load")) {
                        md_bitfield_copy(&data->selection.selection_mask, &sel.atom_mask);
                        flag_all_representations_as_dirty(data);
                    }
                    if (ImGui::IsItemHovered()) {
                        ImGui::SetTooltip("Load the stored selection as the active selection");
                        md_bitfield_copy(&data->selection.highlight_mask, &sel.atom_mask);
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Store")) {
                        ImGui::SetTooltip("Store the active selection into this selection");
                        md_bitfield_copy(&sel.atom_mask, &data->selection.selection_mask);
                        data->script.compile_ir = true;
                        flag_all_representations_as_dirty(data);
                    }
                    ImGui::SameLine();
                    if (ImGui::DeleteButton("Remove")) {
                        ImGui::SetTooltip("Remove this selection");
                        remove_selection(data, i);
                    }
                    ImGui::PopID();
                }

                if (ImGui::Button("Create New")) {
                    char name_buf[64];
                    snprintf(name_buf, sizeof(name_buf), "sel%i", (int)md_array_size(data->selection.stored_selections) + 1);
                    create_selection(data, str_from_cstr(name_buf), &data->selection.selection_mask);
                }

                //ImGui::PopItemFlag();
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Screenshot")) {
            ImGui::Checkbox("Hide GUI", &data->screenshot.hide_gui);
            if (data->screenshot.hide_gui) {
                data->screenshot.resolution = (ScreenshotResolution)MIN((int)data->screenshot.resolution, (int)ScreenshotResolution::Count - 1);
                if (ImGui::BeginCombo("Resolution", screenshot_resolution_str[(int)data->screenshot.resolution])) {
                    for (int i = 0; i < (int)ScreenshotResolution::Count; ++i) {
                        if (ImGui::Selectable(screenshot_resolution_str[i], (i == (int)data->screenshot.resolution))) {
                            data->screenshot.resolution = (ScreenshotResolution)i;
                        }
                    }
                    ImGui::EndCombo();
                }

                switch (data->screenshot.resolution) {
                case ScreenshotResolution::Window:
                    data->screenshot.res_x = data->gbuffer.width;
                    data->screenshot.res_y = data->gbuffer.height;
                    break;
                case ScreenshotResolution::FHD:
                    data->screenshot.res_x = 1920;
                    data->screenshot.res_y = 1080;
                    break;
                case ScreenshotResolution::QHD:
                    data->screenshot.res_x = 2560;
                    data->screenshot.res_y = 1440;
                    break;
                case ScreenshotResolution::UHD_4K:
                    data->screenshot.res_x = 3840;
                    data->screenshot.res_y = 2160;
                    break;
                case ScreenshotResolution::UHD_8K:
                    data->screenshot.res_x = 7680;
                    data->screenshot.res_y = 4320;
                    break;
                case ScreenshotResolution::Custom:
                    ImGui::InputInt("Res X", &data->screenshot.res_x);
                    ImGui::InputInt("Res Y", &data->screenshot.res_y);

                    data->screenshot.res_x = CLAMP(data->screenshot.res_x, 640, 16384);
                    data->screenshot.res_y = CLAMP(data->screenshot.res_y, 480, 16384);
                    break;
                default:
                    ASSERT(false);
                }
            } else {
                data->screenshot.res_x = data->gbuffer.width;
                data->screenshot.res_y = data->gbuffer.height;
            }
            if (ImGui::MenuItem("Take Screenshot", nullptr, false, data->movie.state != MovieRecordingState::Recording)) {
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, STR_LIT("jpg,png,bmp"))) {
                    size_t path_len = strnlen(path_buf, sizeof(path_buf));
                    str_t ext;
                    if (!extract_ext(&ext, {path_buf, path_len})) {
                        path_len += snprintf(path_buf + path_len, sizeof(path_buf) - path_len, ".jpg");
                        ext = STR_LIT("jpg");
                    }
                    if (str_eq_cstr_ignore_case(ext, "jpg") || str_eq_cstr_ignore_case(ext, "png") || str_eq_cstr_ignore_case(ext, "bmp")) {
                        data->screenshot.path_to_file = str_copy({path_buf, path_len}, persistent_alloc);
                        if (data->visuals.temporal_aa.enabled) {
                            data->screenshot.sample_target = JITTER_SEQUENCE_SIZE;
                        } else {
                            data->screenshot.sample_target = 1;
                        }
                    }
                    else {
                        VIAMD_LOG_ERROR("Supplied image extension is not supported");
                    }
                }
                ImGui::GetCurrentWindow()->Hidden = true;
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Operations")) {
            bool do_recenter = false;
            bool do_pbc = false;
            bool do_unwrap = false;
            bool do_bonds = false;
            bool redo_every_frame = false;

            // Each operation can be applied to the frame shown now, or to every frame as it is shown.
            // The labels say what happens to the atoms, not what the operation is called internally.
            ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg;
            if (ImGui::BeginTable("##table", 3, flags)) {
                if (ImGui::IsWindowHovered()) {
                    md_bitfield_clear(&data->selection.highlight_mask);
                }
                ImGui::TableSetupColumn("Now");
                ImGui::TableSetupColumn("Every frame");
                ImGui::TableSetupColumn("Options");
                ImGui::TableHeadersRow();

                const float button_width = (ImGui::GetFontSize() / 20.f) * 175.f;
                const md_bitfield_t& target_mask = recenter_get_active_target_mask(data);
                const bool recenter_available = !md_bitfield_empty(&target_mask);
                const char* no_target_tip =
                    "There is nothing to center on yet. Select atoms and choose 'Set as Centering Target' in the\n"
                    "right-click menu of the viewport, or tick " ICON_FA_COMMENT_DOTS " and pick the target with a query.";

                // ## Centering
                ImGui::TableNextRow();

                ImGui::TableSetColumnIndex(0);
                if (!recenter_available) {
                    ImGui::PushDisabled();
                }
                if (ImGui::Button("Center Target", ImVec2(button_width,0))) {
                    do_recenter = true;
                }
                if (!recenter_available) {
                    ImGui::PopDisabled();
                }
                if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
                    if (!recenter_available) {
                        ImGui::SetTooltip("%s", no_target_tip);
                    } else {
                        ImGui::SetTooltip("Move everything so the target (highlighted) is in the middle of the box,\nor at the origin when there is no box");
                        // Highlight the target atoms that the system will be recentered around
                        md_bitfield_copy (&data->selection.highlight_mask, &target_mask);
                    }
                }

                ImGui::TableSetColumnIndex(1);
                if (!recenter_available) {
                    ImGui::PushDisabled();
                }
                if (ImGui::Checkbox("##recenter", &data->operations.recenter) && data->operations.recenter) {
                    do_recenter = true;
                }
                if (!recenter_available) {
                    ImGui::PopDisabled();
                }
                if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
                    ImGui::SetTooltip("%s", recenter_available ? "Keep the target in the middle of the box in every frame" : no_target_tip);
                }

                ImGui::TableSetColumnIndex(2);
                if (ImGui::Checkbox(ICON_FA_ANCHOR_LOCK "##keep-orientation", &data->operations.fixate_orientation) && data->operations.recenter) {
                    // Turning it on or off while centering every frame takes effect on the frame shown, not the next one
                    redo_every_frame = true;
                }
                ImGui::SetItemTooltip("Keep orientation: when centering, also turn everything so the target\nkeeps the orientation it has in the first frame.\nEverything is wrapped into the box around the target before it is turned,\nand the box turns with it.");

                ImGui::SameLine();
                ImGui::Checkbox(ICON_FA_COMMENT_DOTS "##target-query", &data->operations.recenter_query.enabled);
                ImGui::SetItemTooltip("Target by query: center on the atoms a query picks (evaluated every frame\nwhen it depends on it) instead of the target set from a selection");

                if (data->operations.recenter_query.enabled) {
                    auto& recenter_query = data->operations.recenter_query;
                    ImGui::SameLine();
                    if (ImGui::InputQuery("##recenter-query", recenter_query.query, sizeof(recenter_query.query), recenter_query.valid, recenter_query.error)) {
                        recenter_mark_query_dirty(data);
                    }
                }

                // ## Periodic box
                ImGui::TableNextRow();

                ImGui::TableSetColumnIndex(0);
                if (ImGui::Button("Wrap into Box", ImVec2(button_width,0))) {
                    do_pbc = true;
                }
                ImGui::SetItemTooltip("Move every atom that is outside the periodic box back in through the opposite side\n(apply periodic boundary conditions). Molecules across the boundary are split.");

                ImGui::TableSetColumnIndex(1);
                if (ImGui::Checkbox("##pbc", &data->operations.apply_pbc) && data->operations.apply_pbc) {
                    do_pbc = true;
                }
                ImGui::SetItemTooltip("Wrap every atom into the periodic box in every frame");

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                if (ImGui::Button("Make Whole", ImVec2(button_width, 0))) {
                    do_unwrap = true;
                }
                ImGui::SetItemTooltip("Join molecules that are split across the periodic box boundary,\nso each one is drawn in one piece (unwrap)");

                ImGui::TableSetColumnIndex(1);
                if (ImGui::Checkbox("##unwrap", &data->operations.unwrap_structures) && data->operations.unwrap_structures) {
                    do_unwrap = true;
                }
                ImGui::SetItemTooltip("Keep every molecule whole in every frame");

                // ## Bonds
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                if (ImGui::Button("Recompute Bonds", ImVec2(button_width, 0))) {
                    do_bonds = true;
                }
                ImGui::SetItemTooltip("Guess the covalent bonds again from the distances between atoms in this frame,\nfor when bonds form or break, or the file's bonds are wrong");

                ImGui::TableSetColumnIndex(1);
                if (ImGui::Checkbox("##bonds", &data->operations.recalc_bonds) && data->operations.recalc_bonds) {
                    do_bonds = true;
                }
                ImGui::SetItemTooltip("Guess the covalent bonds from the distances in every frame (slower)");
                ImGui::EndTable();
            }

            if (redo_every_frame) {
                do_recenter = true;
                do_pbc     |= data->operations.apply_pbc;
                do_unwrap  |= data->operations.unwrap_structures;
            }

            // One pass, in the order that keeps them valid together, and in the lattice frame of the cell
            // even when the coordinates are turned (see apply_state_operations)
            if (apply_state_operations(data, do_recenter, do_pbc, do_unwrap)) {
                data->mold.dirty_gpu_buffers |= MolBit_DirtyPosition | MolBit_ClearVelocity;
            }

            if (do_bonds) {
                if (!task_system::task_is_running(data->tasks.evaluate_full) && !task_system::task_is_running(data->tasks.evaluate_filt)) {
                    // The whole frame nearest the animation time, or the coordinates shown when there is no run
                    recompute_covalent_bonds(data, (int64_t)(data->animation.frame + 0.5));
                } else {
                    MD_LOG_INFO("Cannot recalculate bonds while evaluation is occuring.");
                }
            }

            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Settings")) {
            if (ImGui::Checkbox("Keep Representations", &data->settings.keep_representations)) {
                app_settings::mark_dirty();
            }
            ImGui::SetItemTooltip("Keep representations when loading new topology (does not apply for workspaces)\n");

            // Font
            int font_size_idx = nearest_font_size_index(data->settings.font_size);
            if (ImGui::Combo("Font Size", &font_size_idx, font_size_names, (int)ARRAY_SIZE(font_size_names))) {
                data->settings.font_size = font_sizes[font_size_idx];
                apply_font_size(data);
                app_settings::mark_dirty();
            }

			if (ImGui::TreeNode("Units")) {
                display_units::draw_settings_menu_items();
				ImGui::TreePop();
			}

            ImGui::EndMenu();
        }
        {
            // Fps counter
            static int num_frames = 0;
            static double acc_ms  = 0;
            static double avg_ms  = 0;

            double ms = (data->app.timing.delta_s * 1000);
            acc_ms += ms;
            num_frames += 1;

            if (acc_ms > 500) {
                avg_ms = acc_ms / num_frames;
                acc_ms = 0;
                num_frames = 0;
            }

            char fps_buf[64];
            snprintf(fps_buf, ARRAY_SIZE(fps_buf), "%.2f ms (%.1f fps)", avg_ms, 1000.f / avg_ms);
            const float w = ImGui::CalcTextSize(fps_buf).x;
            ImGui::SetCursorPosX(ImGui::GetWindowContentRegionMax().x - w);
            ImGui::Text("%s", fps_buf);
            if (ImGui::IsItemHovered()) {
                // GPU time of the volume passes, all views together (main viewport and component windows)
                const volume::GpuTimings t = volume::timings_get();
                ImGui::BeginTooltip();
                ImGui::Text("Volume rendering, GPU per frame: %.3f ms", t.total_ms);
                for (int i = 0; i < volume::TimingStage_Count; ++i) {
                    ImGui::Text("  %-14s %.3f ms", volume::timing_stage_name((volume::TimingStage)i), t.ms[i]);
                }
                if (data->visuals.ssao.enabled) {
                    const postprocess_pipeline::SsaoTimings ao = postprocess_pipeline::ssao_timings();
                    ImGui::Text("Ambient occlusion, GPU per frame: %.3f ms", ao.total_ms);
                    ImGui::Text("  %-14s %.3f ms", "Estimator", ao.estimator_ms);
                }
                ImGui::EndTooltip();
            }
        }
        ImGui::EndMainMenuBar();
    }

    if (new_clicked) ImGui::OpenPopup("Warning New");
}

void draw_notifications_window() {
    // Render toasts on top of everything, at the end of your code!
    // You should push style vars here
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 5.f);
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(43.f / 255.f, 43.f / 255.f, 43.f / 255.f, 100.f / 255.f));
    ImGui::RenderNotifications();
    ImGui::PopStyleVar(1); // Don't forget to Pop()
    ImGui::PopStyleColor(1);
}

void draw_load_dataset_window(ApplicationState* data) {
    LoadDatasetWindowState& state = data->load_dataset;
    
    ImGui::SetNextWindowSize(ImVec2(300, 215), ImGuiCond_FirstUseEver);
    if (state.show_window) {
        ImGui::OpenPopup("Load Dataset");
    }

    if (ImGui::BeginPopupModal("Load Dataset", &state.show_window)) {
        const bool path_invalid = !state.path_is_valid && state.path_buf[0] != '\0';
        const int  loader_count = LoaderType_COUNT;

        if (path_invalid) ImGui::PushInvalid();
        if (ImGui::InputText("##path", state.path_buf, sizeof(state.path_buf))) {
            state.path_changed = true;
        }
        if (path_invalid) ImGui::PopInvalid();


        // @WORKAROUND(Robin): This show_file_dialog is only here to circumvent the issue that if you open a file dialog
        // Within the same frame as the button is clicked, the dialogue will open again after the closing the dialog.
        ImGui::SameLine();
        if (ImGui::Button("Browse") && !state.show_file_dialog) {
            state.show_file_dialog = true;
        }

        if (state.show_file_dialog) {
            state.show_file_dialog = false;
            if (application::file_dialog(state.path_buf, sizeof(state.path_buf), application::FileDialogFlag_Open)) {
                state.path_changed = true;
            }
        }

        str_t path = str_from_cstr(state.path_buf);

        if (state.path_changed) {
            state.path_changed = false;
            state.path_is_valid = md_path_is_valid(path) && !md_path_is_directory(path);

            // Try to assign loader_idx from extension
            state.loader_idx = 0;
            str_t ext;
            if (extract_ext(&ext, path)) {
                state.loader_idx = (int)loader::type_from_ext(ext);
            }
        }

        if (ImGui::BeginCombo("Loader", loader::type_name(state.loader_idx).ptr)) {
            for (int i = 1; i < loader_count; ++i) {
                const char* str = loader::type_name(i).ptr;
                if (!str) continue;
                if (ImGui::Selectable(str, state.loader_idx == i)) {
                    state.loader_idx = i;
                }
            }
            ImGui::EndCombo();
        }

        // True if the button should be enabled
        bool load_button_enabled = (state.path_is_valid && state.loader_idx > -1);

        LoaderType   type = (LoaderType)(state.loader_idx);
        LoaderFlags flags = loader::type_flags(type);

        // Draw Options
        bool show_cg = state.path_is_valid && (flags & LoaderFlag_System);
        if (show_cg) {
            ImGui::Checkbox("Coarse Grained", &state.coarse_grained);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Enable if the data should be interpreted as coarse grained particles");
            }
        }

        bool show_lammps_atom_format = state.path_is_valid && (type == LoaderType_LAMMPSDATA);
        if (show_lammps_atom_format) {
            const char** atom_format_names = md_lammps_atom_format_names();
            const char** atom_format_strings = md_lammps_atom_format_strings();

            if (state.atom_format_idx == -1) {
                // Try to determine format from file
                md_lammps_atom_format_t format = md_lammps_atom_format_from_file(path);
                state.atom_format_idx = format;
				const char* format_str = atom_format_strings[format];
                snprintf(state.atom_format_buf, sizeof(state.atom_format_buf), "%s", format_str);
            }
            state.atom_format_idx = CLAMP(state.atom_format_idx, 0, MD_LAMMPS_ATOM_FORMAT_COUNT - 1);

            if (ImGui::BeginCombo("Atom Format", state.atom_format_idx > 0 ? atom_format_names[state.atom_format_idx] : "user defined")) {
                for (int i = 0; i < MD_LAMMPS_ATOM_FORMAT_COUNT; ++i) {
                    if (ImGui::Selectable(i > 0 ? atom_format_names[i] : "user defined", state.atom_format_idx == i)) {
                        state.atom_format_idx = i;
                        int source_idx = i > 0 ? i : MD_LAMMPS_ATOM_FORMAT_FULL;
                        const char* format_str = atom_format_strings[source_idx];
                        snprintf(state.atom_format_buf, sizeof(state.atom_format_buf), "%s", format_str);
                    }
                }
                ImGui::EndCombo();
            }

            if (state.atom_format_idx == MD_LAMMPS_ATOM_FORMAT_UNKNOWN) {
                bool valid = state.atom_format_valid;
                if (!valid) ImGui::PushInvalid();
                if (ImGui::InputText("##atom_format", state.atom_format_buf, sizeof(state.atom_format_buf))) {
                    state.atom_format_valid = md_lammps_validate_atom_format(state.err_buf, sizeof(state.err_buf), state.atom_format_buf);
                }
                if (!valid) ImGui::PopInvalid();
                if (ImGui::IsItemHovered() && !valid) {
                    ImGui::SetTooltip("%s", state.err_buf);
                }
            }
            else {
                ImGui::PushDisabled();
                ImGui::InputText("##atom_format", state.atom_format_buf, sizeof(state.atom_format_buf), ImGuiInputTextFlags_ReadOnly);
                ImGui::PopDisabled();
            }

            if (state.atom_format_idx < 0 || (state.atom_format_idx == MD_LAMMPS_ATOM_FORMAT_UNKNOWN && !state.atom_format_valid)) {
                load_button_enabled = false;
            }
        }

        enum Action {
            Action_None,
            Action_Cancel,
            Action_Load,
        };
        Action action = Action_None;

        if (!load_button_enabled) ImGui::PushDisabled();
        if (ImGui::Button("Load")) {
            action = Action_Load;
        }
        if (!load_button_enabled) ImGui::PopDisabled();

        ImGui::SameLine();
        if (ImGui::Button("Cancel")) {
            action = Action_Cancel;
        }

        switch (action) {
        case Action_Load: {
            loader::LoaderState load_state = {};
            load_state.type = type;
            load_state.flags = flags;
            load_state.flags |= state.coarse_grained ? LoaderFlag_CoarseGrained : 0;
            if (type == LoaderType_LAMMPSDATA) {
                load_state.arg = state.atom_format_buf;
            }

            if (load_data_from_file(data, path, load_state)) {
                if ((flags & LoaderFlag_System) && !data->settings.keep_representations) {
                    remove_all_representations(data);
                    create_default_representations(data);
                }
                data->animation = {};
                reset_view(&data->view.target, data->mold.state, &data->representation.visibility_mask);
            }
            [[fallthrough]];
        }
        case Action_Cancel:
            // Reset state
            // @NOTE(Robin): Don't change this to {}, it won't work on GCC 9
            state = LoadDatasetWindowState();
            [[fallthrough]];
        case Action_None:
        default:
            break;
        }

        ImGui::EndPopup();
    }
}

/*
// TODO: Move these functions to dataset component
void clear_atom_elem_mappings(ApplicationState* data) {
    md_array_shrink(data->dataset.atom_element_remappings, 0);
}

AtomElementMapping* add_atom_elem_mapping(ApplicationState* data, str_t lbl, md_element_t elem) {
    // Check if we already have a mapping for the label -> overwrite
    size_t i = 0;
    for (; i < md_array_size(data->dataset.atom_element_remappings); ++i) {
        if (str_eq_cstr(lbl, data->dataset.atom_element_remappings[i].lbl)) break;
    }
    if (i == md_array_size(data->dataset.atom_element_remappings)) {
        AtomElementMapping mapping = {
            .elem = elem,
        };
        str_copy_to_char_buf(mapping.lbl, sizeof(mapping.lbl), lbl);
        md_array_push(data->dataset.atom_element_remappings, mapping, persistent_alloc);
        return md_array_last(data->dataset.atom_element_remappings);
    } else {
        data->dataset.atom_element_remappings[i].elem = elem;
        return &data->dataset.atom_element_remappings[i];
    }
}

void apply_atom_elem_mappings(ApplicationState* data) {
    if (data->mold.sys.atom.count == 0 || !data->mold.sys.atom.element) {
        return;
    }

    for (size_t j = 0; j < md_array_size(data->dataset.atom_element_remappings); ++j) {
        str_t lbl = str_from_cstr(data->dataset.atom_element_remappings[j].lbl);
        md_element_t elem = data->dataset.atom_element_remappings[j].elem;
        float radius = md_util_element_vdw_radius(elem);
        float mass = md_util_element_atomic_mass(elem);

        for (size_t i = 0; i < data->mold.sys.atom.count; ++i) {
            if (str_eq(lbl, data->mold.sys.atom.type[i])) {
                data->mold.sys.atom.element[i] = elem;
                data->mold.sys.atom.radius[i] = radius;
                data->mold.sys.atom.mass[i] = mass;
                data->mold.dirty_buffers |= MolBit_DirtyRadius;
            }
        }
    }
    md_system_t* mol = &data->mold.sys;
    
    
    md_array_free(mol->bond.pairs, data->mold.sys.alloc);
    md_array_free(mol->bond.order, data->mold.sys.alloc);

    md_array_free(mol->bond.conn.atom_idx, data->mold.sys.alloc);
    md_array_free(mol->bond.conn.bond_idx, data->mold.sys.alloc);

    md_index_data_free(&mol->structure);
    md_index_data_free(&mol->ring);
    
    md_system_state_t mol_state = {0};
    md_util_system_infer(mol, &mol_state, data->mold.sys.alloc, MD_UTIL_INFER_BOND_BIT | MD_UTIL_INFER_STRUCTURE_BIT);
    data->mold.dirty_buffers |= MolBit_DirtyBonds;

    flag_all_representations_as_dirty(data);
}
*/

// Create a textual script describing a selection from a bitfield with respect to some reference index
// @TODO(Robin): Clean this up, it is a mess. Just provide complete suggestions based on bitfield and molecule input.
static void write_script_range(md_strb_t& sb, const int* indices, size_t num_indices, int ref_idx = 0) {
    if (num_indices == 0) return;
    if (num_indices == 1) {
        md_strb_fmt(&sb, "%i", indices[0] - ref_idx + 1);
        return;
    }
    
    int range_beg = indices[0];
    int prev_idx  = -1;

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    md_array(md_irange_t) items = 0;
    
    for (size_t i = 0; i < num_indices; ++i) {
        int idx = indices[i];
        
        if (idx - prev_idx > 1) {
            if (prev_idx > range_beg) {
                md_irange_t item = {range_beg - ref_idx + 1, prev_idx - ref_idx + 1};
                md_array_push(items, item, temp_alloc);
            } else if (prev_idx != -1) {
                md_irange_t item = {prev_idx - ref_idx + 1, prev_idx - ref_idx + 1};
                md_array_push(items, item, temp_alloc);
            }
            range_beg = idx;
        }

        prev_idx = idx;
    }

    if (prev_idx - range_beg > 0) {
        md_irange_t item = {range_beg - ref_idx + 1, prev_idx - ref_idx + 1};
        md_array_push(items, item, temp_alloc);
    } else if (prev_idx != -1) {
        md_irange_t item = {prev_idx - ref_idx + 1, prev_idx - ref_idx + 1};
        md_array_push(items, item, temp_alloc);
    }

    const int64_t num_items = (int64_t)md_array_size(items);
    if (num_items > 1) sb += "{";
    for (int64_t i = 0; i < num_items; ++i) {
        md_irange_t item = items[i];
        if (item.beg == item.end) {
            md_strb_fmt(&sb, "%i", item.beg);
        }
        else {
            md_strb_fmt(&sb, "%i:%i", item.beg, item.end);
        }
        if (i < num_items - 1) {
            sb += ',';
        }
    }
    if (num_items > 1) sb += "}";
    md_temp_end(temp);
}

static md_array(str_t) generate_script_selection_suggestions(str_t ident, const md_bitfield_t* bf, const md_system_t* sys) {
    md_array(str_t) suggestions = 0;

    bool within_same_comp = true;
    bool within_same_inst = true;

    md_instance_idx_t inst_idx = -1;
    md_component_idx_t comp_idx = -1;

    md_bitfield_iter_t it = md_bitfield_iter_create(bf);
    while (md_bitfield_iter_next(&it)) {
        uint64_t a_idx = md_bitfield_iter_idx(&it);
        int32_t  i_idx = md_system_instance_find_by_atom_idx(sys, a_idx);
        int32_t  c_idx = md_system_component_find_by_atom_idx(sys, a_idx);

        if (inst_idx == -1 && comp_idx != -1) {
            inst_idx = i_idx;
        } else if (c_idx != -1 && inst_idx != i_idx) {
            within_same_inst = false;
        }

        if (comp_idx == -1 && c_idx != -1) {
            comp_idx = c_idx;
        } else if (c_idx != -1 && comp_idx != c_idx) {
            within_same_comp = false;                        
        }

        if (!within_same_comp && !within_same_inst) {
            break;
        }
    }
    
    const size_t popcount = md_bitfield_popcount(bf);
    
    md_strb_t sb = md_strb_create(frame_alloc);
    defer { md_strb_free(&sb); };

    auto write_atom_remainder = [](md_strb_t& sb, const md_bitfield_t* bf, int ref_idx = 0) {
        // Add any remainder
        const uint64_t remainder = md_bitfield_popcount(bf);
        if (remainder) {
            int* indices = (int*)md_alloc(frame_alloc, remainder * sizeof(int));
            defer { md_free(frame_alloc, indices, remainder * sizeof(int)); };

            md_bitfield_iter_extract_indices(indices, remainder, md_bitfield_iter_create(bf));
            sb += "atom(";
            write_script_range(sb, indices, remainder, ref_idx);
            sb += ")";
        }
    };

    if (comp_idx != -1 && within_same_comp) {
        const md_urange_t range = md_component_atom_range(&sys->component, comp_idx);
        if (popcount != range.end - range.beg) {
            md_strb_reset(&sb);
            sb += ident;
            sb += " = ";

            // Subset of residue is selected
            write_atom_remainder(sb, bf, range.beg);
            if (md_strb_len(sb) < 512) {
                str_t resname = md_component_name(&sys->component, comp_idx);
                md_strb_fmt(&sb, " in resname(\"" STR_FMT "\");", STR_ARG(resname));
                md_array_push(suggestions, str_copy(md_strb_to_str(sb), frame_alloc), frame_alloc);
            }
        }
    }

    else if (inst_idx != -1 && within_same_inst) {
        const md_urange_t range = md_system_instance_atom_range(sys, inst_idx);
        if (popcount != range.end - range.beg) {
            md_strb_reset(&sb);
            sb += ident;
            sb += " = ";

            // Subset of chain is selected
            write_atom_remainder(sb, bf, range.beg);
            if (md_strb_len(sb) < 512) {
				str_t inst_id = md_instance_id(&sys->instance, inst_idx);
				md_strb_fmt(&sb, " in chain(\"" STR_FMT "\");", STR_ARG(inst_id));
                md_array_push(suggestions, str_copy(md_strb_to_str(sb), frame_alloc), frame_alloc);
            }
        }
    }

    md_bitfield_t tmp_bf = md_bitfield_create(frame_alloc);
    md_bitfield_copy(&tmp_bf, bf);
    
    md_array(int) complete_chains = 0;
    md_array(int) complete_residues = 0;
    
    if (sys->instance.count) {
        for (size_t i = 0; i < sys->instance.count; ++i) {    
            const md_urange_t range = md_system_instance_atom_range(sys, i);
            if (md_bitfield_test_all_range(&tmp_bf, range.beg, range.end)) {
                md_array_push(complete_chains, (int)i, frame_alloc);
                md_bitfield_clear_range(&tmp_bf, range.beg, range.end);
            }
        }
    }

    if (sys->component.count) {
        for (size_t i = 0; i < sys->component.count; ++i) {    
            const md_urange_t range = md_component_atom_range(&sys->component, i);
            if (md_bitfield_test_all_range(&tmp_bf, range.beg, range.end)) {
                md_array_push(complete_residues, (int)i, frame_alloc);
                md_bitfield_clear_range(&tmp_bf, range.beg, range.end);
            }
        }
    }

    const uint64_t atom_remainder_count = md_bitfield_popcount(&tmp_bf);
    
    if (complete_chains) {
        md_strb_reset(&sb);
        sb += ident;
        sb += " = chain(";
        write_script_range(sb, complete_chains, md_array_size(complete_chains));
        sb += ")";

        if (complete_residues) {
            sb += " or residue(";
            write_script_range(sb, complete_residues, md_array_size(complete_residues));
            sb += ")";
        }

        if (atom_remainder_count) {
            sb += " or ";
            write_atom_remainder(sb, &tmp_bf);
        }
        
        if (md_strb_len(sb) < 512) {
            md_strb_push_char(&sb, ';');
            md_array_push(suggestions, str_copy(md_strb_to_str(sb), frame_alloc), frame_alloc);
        }
        
        md_strb_reset(&sb);
        sb += ident;
        sb += " = residue(";
        for (size_t i = 0; i < md_array_size(complete_chains); ++i) {
            md_urange_t range = md_instance_component_range(&sys->instance, complete_chains[i]);
            md_strb_fmt(&sb, "%i:%i,", range.beg + 1, range.end);
        }
        if (complete_residues) {
            write_script_range(sb, complete_residues, md_array_size(complete_residues));
        } else {
            md_strb_pop(&sb, 1);
        }
        sb += ")";

        if (atom_remainder_count) {
            sb += " or ";
            write_atom_remainder(sb, &tmp_bf);
        }
        
        if (md_strb_len(sb) < 512) {
            md_strb_push_char(&sb, ';');
            md_array_push(suggestions, str_copy(md_strb_to_str(sb), frame_alloc), frame_alloc);
        }
    }

    if (complete_residues) {
        md_strb_reset(&sb);
        sb += ident;
        sb += " = residue(";
        write_script_range(sb, complete_residues, md_array_size(complete_residues));
        sb += ")";

        if (atom_remainder_count) {
            sb += " or ";
            write_atom_remainder(sb, &tmp_bf);
        }

        if (md_strb_len(sb) < 512) {
            md_strb_push_char(&sb, ';');
            md_array_push(suggestions, str_copy(md_strb_to_str(sb), frame_alloc), frame_alloc);
        }
    }

    if (popcount) {
        md_strb_reset(&sb);
        sb += ident;
        sb += " = ";
        write_atom_remainder(sb, bf);
        if (md_strb_len(sb) < 512) {
            md_strb_push_char(&sb, ';');
            md_array_push(suggestions, str_copy(md_strb_to_str(sb), frame_alloc), frame_alloc);
        }
    }
    
    return suggestions;
}

static int64_t find_identifier(const md_script_ir_t* ir, str_t ident) {
    const int64_t num_ident = md_script_ir_num_identifiers(ir);
    const str_t* idents = md_script_ir_identifiers(ir);
    for (int64_t i = 0; i < num_ident; ++i) {
        if (str_eq(ident, idents[i])) return i;
    }
    return -1;
}

static str_t create_unique_identifier(const md_script_ir_t* ir, str_t base, md_allocator_i* alloc) {
    char buf[128];
    for (int64_t i = 1; i < 10; ++i) {
        int res = snprintf(buf, sizeof(buf), "%.*s%i", (int)base.len, base.ptr, (int)i);
        if (res > 0) {
            str_t ident = {buf, (size_t)res};
            if (find_identifier(ir, ident) == -1) {
                return str_copy(ident, alloc);
            }
        }
    }
    return str_t();
}

// # context_menu
void draw_context_popup(ApplicationState* state, const PickingHit& hit) {
    // @Robin: Will hit ever be used here? (Possibly remove argument)
    (void)hit;
    ASSERT(state);

    if (!state->mold.sys.atom.count) return;

    const size_t sss_count = single_selection_sequence_count(&state->selection.single_selection_sequence);
    const size_t num_atoms_selected = md_bitfield_popcount(&state->selection.selection_mask);

    if (ImGui::BeginPopup("Context Popup")) {
        if (num_atoms_selected == 2) {
            // Suggest construction of covalent bond
            int idx[2];
            md_bitfield_iter_extract_indices(idx, 2, md_bitfield_iter_create(&state->selection.selection_mask));

            md_bond_idx_t bond_idx = md_system_bond_find(&state->mold.sys, idx[0], idx[1]);
            if (bond_idx == -1) {
                char buf[256];
                snprintf(buf, sizeof(buf), "Create Bond (%i, %i)", idx[0] + 1, idx[1] + 1);
                if (ImGui::MenuItem(buf)) {
                    md_system_bond_insert(&state->mold.sys, idx[0], idx[1], md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_USER));
                    md_util_system_infer_coordination(&state->mold.sys);
                    md_system_bond_build_connectivity(&state->mold.sys);
                    state->mold.dirty_gpu_buffers |= MolBit_DirtyBonds;
                    ImGui::CloseCurrentPopup();
                }
            } else {
				md_bond_flags_t flags = md_system_bond_flags(&state->mold.sys, bond_idx);
                if (md_bond_origin(flags) == MD_BOND_ORIGIN_USER) {
                    char buf[256];
                    snprintf(buf, sizeof(buf), "Remove Bond (%i, %i)", idx[0] + 1, idx[1] + 1);
                    if (ImGui::MenuItem(buf)) {
                        md_system_bond_remove(&state->mold.sys, bond_idx);
                        md_system_bond_build_connectivity(&state->mold.sys);
                        state->mold.dirty_gpu_buffers |= MolBit_DirtyBonds;
                        ImGui::CloseCurrentPopup();
                    }
                }
            }
        }
        if (ImGui::BeginMenu("Script")) {
            bool any_suggestions = false;
            if (num_atoms_selected <= 4 && sss_count > 1) {
                int idx[4];
                MEMCPY(idx, state->selection.single_selection_sequence.idx, sizeof(idx));
                // Check if all selected atoms are within the same residue
                md_component_idx_t res_idx = md_component_find_by_atom_idx(&state->mold.sys.component, idx[0]);
                for (size_t i = 1; i < sss_count; ++i) {
                    md_component_idx_t ri = md_component_find_by_atom_idx(&state->mold.sys.component, idx[i]);
                    if (res_idx != ri) {
                        res_idx = -1;
                        break;
                    }
                }

                any_suggestions = true;
                char buf[256] = "";
                if (sss_count == 2) {
                    str_t ident = create_unique_identifier(state->script.ir, STR_LIT("dist"), frame_alloc);

                    snprintf(buf, sizeof(buf), STR_FMT " = distance(%i, %i);", STR_ARG(ident), idx[0]+1, idx[1]+1);
                    if (ImGui::MenuItem(buf)) {
                        script_editor::append_line(state->editor, str_from_cstr(buf));
                        ImGui::CloseCurrentPopup();
                    }
                    if (ImGui::IsItemHovered()) {
                        str_t str = str_from_cstr(buf);
                        script_visualize_str(state, str, VIS_FLAGS);
                    }

                    if (res_idx != -1) {
                        const md_urange_t range = md_component_atom_range(&state->mold.sys.component, res_idx);
                        idx[0] -= range.beg;
                        idx[1] -= range.beg;

                        snprintf(buf, sizeof(buf), STR_FMT " = distance(%i, %i) in residue(%i);", STR_ARG(ident), idx[0]+1, idx[1]+1, res_idx+1);
                        if (ImGui::MenuItem(buf)) {
                            script_editor::append_line(state->editor, str_from_cstr(buf));
                            ImGui::CloseCurrentPopup();
                        }
                        if (ImGui::IsItemHovered()) {
                            str_t str = str_from_cstr(buf);
                            script_visualize_str(state, str, VIS_FLAGS);
                        }

                        const int32_t resid = md_component_seq_id(&state->mold.sys.component, res_idx);
                        snprintf(buf, sizeof(buf), STR_FMT " = distance(%i, %i) in resid(%i);", STR_ARG(ident), idx[0]+1, idx[1]+1, resid);
                        if (ImGui::MenuItem(buf)) {
                            script_editor::append_line(state->editor, str_from_cstr(buf));
                            ImGui::CloseCurrentPopup();
                        }
                        if (ImGui::IsItemHovered()) {
                            str_t str = str_from_cstr(buf);
                            script_visualize_str(state, str, VIS_FLAGS);
                        }

                        str_t resname = md_component_name(&state->mold.sys.component, res_idx);
                        if (resname) {
                            snprintf(buf, sizeof(buf), STR_FMT " = distance(%i, %i) in resname(\"" STR_FMT "\");", STR_ARG(ident), idx[0]+1, idx[1]+1, STR_ARG(resname));
                            if (ImGui::MenuItem(buf)) {
                                script_editor::append_line(state->editor, str_from_cstr(buf));
                                ImGui::CloseCurrentPopup();
                            }
                            if (ImGui::IsItemHovered()) {
                                str_t str = str_from_cstr(buf);
                                script_visualize_str(state, str, VIS_FLAGS);
                            }
                        }
                    }
                }
                else if(sss_count == 3) {
                    str_t ident = create_unique_identifier(state->script.ir, STR_LIT("ang"), frame_alloc);

                    snprintf(buf, sizeof(buf), STR_FMT " = angle(%i, %i, %i);", STR_ARG(ident), idx[0]+1, idx[1]+1, idx[2]+1);
                    if (ImGui::MenuItem(buf)) {
                        script_editor::append_line(state->editor, str_from_cstr(buf));
                        ImGui::CloseCurrentPopup();
                    }
                    if (ImGui::IsItemHovered()) {
                        str_t str = str_from_cstr(buf);
                        script_visualize_str(state, str, VIS_FLAGS);
                    }

                    if (res_idx != -1) {
                        const md_urange_t range = md_component_atom_range(&state->mold.sys.component, res_idx);
                        idx[0] -= range.beg;
                        idx[1] -= range.beg;
                        idx[2] -= range.beg;

                        snprintf(buf, sizeof(buf), STR_FMT " = angle(%i, %i, %i) in residue(%i);", STR_ARG(ident), idx[0]+1, idx[1]+1, idx[2]+1, res_idx+1);
                        if (ImGui::MenuItem(buf)) {
                            script_editor::append_line(state->editor, str_from_cstr(buf));
                            ImGui::CloseCurrentPopup();
                        }
                        if (ImGui::IsItemHovered()) {
                            str_t str = str_from_cstr(buf);
                            script_visualize_str(state, str, VIS_FLAGS);
                        }

                        int32_t resid = md_component_seq_id(&state->mold.sys.component, res_idx);
                        snprintf(buf, sizeof(buf), STR_FMT " = angle(%i, %i, %i) in resid(%i);", STR_ARG(ident), idx[0]+1, idx[1]+1, idx[2]+1, resid);
                        if (ImGui::MenuItem(buf)) {
                            script_editor::append_line(state->editor, str_from_cstr(buf));
                            ImGui::CloseCurrentPopup();
                        }
                        if (ImGui::IsItemHovered()) {
                            str_t str = str_from_cstr(buf);
                            script_visualize_str(state, str, VIS_FLAGS);
                        }

                        str_t resname = md_component_name(&state->mold.sys.component, res_idx);
                        if (resname) {
                            snprintf(buf, sizeof(buf), STR_FMT " = angle(%i, %i, %i) in resname(\"" STR_FMT "\");", STR_ARG(ident), idx[0]+1, idx[1]+1, idx[2]+1, STR_ARG(resname));
                            if (ImGui::MenuItem(buf)) {
                                script_editor::append_line(state->editor, str_from_cstr(buf));
                                ImGui::CloseCurrentPopup();
                            }
                            if (ImGui::IsItemHovered()) {
                                str_t str = str_from_cstr(buf);
                                script_visualize_str(state, str, VIS_FLAGS);
                            }
                        }
                    }
                }
                else if(sss_count == 4) {
                    str_t ident = create_unique_identifier(state->script.ir, STR_LIT("dih"), frame_alloc);

                    snprintf(buf, sizeof(buf), "%.*s = dihedral(%i, %i, %i, %i);", STR_ARG(ident), idx[0]+1, idx[1]+1, idx[2]+1, idx[3]+1);
                    if (ImGui::MenuItem(buf)) {
                        script_editor::append_line(state->editor, str_from_cstr(buf));
                        ImGui::CloseCurrentPopup();
                    }
                    if (ImGui::IsItemHovered()) {
                        str_t str = str_from_cstr(buf);
                        script_visualize_str(state, str, VIS_FLAGS);
                    }

                    if (res_idx != -1) {
                        const md_urange_t range = md_component_atom_range(&state->mold.sys.component, res_idx);
                        idx[0] -= range.beg;
                        idx[1] -= range.beg;
                        idx[2] -= range.beg;
                        idx[3] -= range.beg;

                        snprintf(buf, sizeof(buf), STR_FMT " = dihedral(%i, %i, %i, %i) in residue(%i);", STR_ARG(ident), idx[0]+1, idx[1]+1, idx[2]+1, idx[3]+1, res_idx+1);
                        if (ImGui::MenuItem(buf)) {
                            script_editor::append_line(state->editor, str_from_cstr(buf));
                            ImGui::CloseCurrentPopup();
                        }
                        if (ImGui::IsItemHovered()) {
                            str_t str = str_from_cstr(buf);
                            script_visualize_str(state, str, VIS_FLAGS);
                        }

                        int32_t resid = md_component_seq_id(&state->mold.sys.component, res_idx);
                        snprintf(buf, sizeof(buf), STR_FMT " = dihedral(%i, %i, %i, %i) in resid(%i);", STR_ARG(ident), idx[0]+1, idx[1]+1, idx[2]+1, idx[3]+1, resid);
                        if (ImGui::MenuItem(buf)) {
                            script_editor::append_line(state->editor, str_from_cstr(buf));
                            ImGui::CloseCurrentPopup();
                        }
                        if (ImGui::IsItemHovered()) {
                            str_t str = str_from_cstr(buf);
                            script_visualize_str(state, str, VIS_FLAGS);
                        }

                        str_t resname = md_component_name(&state->mold.sys.component, res_idx);
                        if (resname) {
                            snprintf(buf, sizeof(buf), STR_FMT " = dihedral(%i, %i, %i, %i) in resname(\"" STR_FMT "\");", STR_ARG(ident), idx[0]+1, idx[1]+1, idx[2]+1, idx[3]+1, STR_ARG(resname));
                            if (ImGui::MenuItem(buf)) {
                                script_editor::append_line(state->editor, str_from_cstr(buf));
                                ImGui::CloseCurrentPopup();
                            }
                            if (ImGui::IsItemHovered()) {
                                str_t str = str_from_cstr(buf);
                                script_visualize_str(state, str, VIS_FLAGS);
                            }
                        }
                    }
                }
            }
            if (num_atoms_selected >= 1) {
                const md_bitfield_t* bf = &state->selection.selection_mask;
                str_t ident = create_unique_identifier(state->script.ir, STR_LIT("sel"), frame_alloc);
                
                md_array(str_t) suggestions = generate_script_selection_suggestions(ident, bf, &state->mold.sys);

			    char buf[128]; // Buffer for limiting menu item size
                for (size_t i = 0; i < md_array_size(suggestions); ++i) {
					const char* ptr = suggestions[i].ptr;
                    str_t s = suggestions[i];
                    if (str_len(s) > 120) {
						snprintf(buf, sizeof(buf), "%.*s...", 117, s.ptr);
						ptr = buf;
                    }
                    if (ImGui::MenuItem(ptr)) {
                        script_editor::append_line(state->editor, s);
                        ImGui::CloseCurrentPopup();
                    }
                }

                any_suggestions = any_suggestions || md_array_size(suggestions) > 0;
            }
            if (!any_suggestions) {
                ImGui::Text("No suggestions for current selection");
            }
            ImGui::EndMenu();
        }

        /*
        if (data->selection.atom_idx.right_click != -1 && data->mold.sys.atom.element) {
            int idx = data->selection.atom_idx.right_click;
            if (0 <= idx && idx < (int)data->mold.sys.atom.count) {
                char label[64] = "";
                str_t type = data->mold.sys.atom.type[idx];
                snprintf(label, sizeof(label), "Remap Element for '%.*s'", (int)type.len, type.ptr);
                if (ImGui::BeginMenu(label)) {
                    static char input_buf[32] = "";
                    md_element_t elem = data->mold.sys.atom.element[idx];
                    str_t name = md_util_element_name(elem);
                    str_t sym  = md_util_element_symbol(elem);

                    ImGui::Text("Current Element: %.*s (%.*s)", (int)name.len, name.ptr, (int)sym.len, sym.ptr);

                    str_t elem_str = {input_buf, strnlen(input_buf, sizeof(input_buf))};
                    md_element_t new_elem = md_util_element_lookup(elem_str);
                    const bool is_valid = new_elem != 0;

                    ImGui::InputQuery("##Symbol", input_buf, sizeof(input_buf), is_valid, "Cannot recognize Element symbol");
                    str_t new_name = md_util_element_name(new_elem);
                    str_t new_sym  = md_util_element_symbol(new_elem);
                    ImGui::Text("New Element: %.*s (%.*s)", (int)new_name.len, new_name.ptr, (int)new_sym.len, new_sym.ptr);
                    if (!is_valid) ImGui::PushDisabled();
                    if (ImGui::Button("Apply") && is_valid) {
                        add_atom_elem_mapping(data, type, new_elem);
                        apply_atom_elem_mappings(data);
                        ImGui::CloseCurrentPopup();
                    }
                    if (!is_valid) ImGui::PopDisabled();
                    ImGui::EndMenu();
                }
            }
        }
        */

        if (ImGui::BeginMenu("Selection")) {
            if (ImGui::MenuItem("Invert")) {
                md_bitfield_not_inplace(&state->selection.selection_mask, 0, state->mold.sys.atom.count);
                ImGui::CloseCurrentPopup();
            }
            if (ImGui::IsItemHovered()) {
                md_bitfield_not(&state->selection.highlight_mask, &state->selection.selection_mask, 0, state->mold.sys.atom.count);
            }
            if (ImGui::MenuItem("Query")) {
                state->selection.query.show_window = true;
                ImGui::CloseCurrentPopup();
            }
            if (num_atoms_selected > 0) {
                if (ImGui::MenuItem("Grow")) {
                    state->selection.grow.show_window = true;
                    ImGui::CloseCurrentPopup();
                }
                if (ImGui::MenuItem("Clear")) {
                    md_bitfield_clear(&state->selection.selection_mask);
                    single_selection_sequence_clear(&state->selection.single_selection_sequence);
                }
            }
            ImGui::EndMenu();
        }
        
        if (num_atoms_selected > 0) {
            if (ImGui::MenuItem("Set as Centering Target")) {
                md_bitfield_clear(&state->operations.selection_mask);
                md_bitfield_copy(&state->operations.selection_mask, &state->selection.selection_mask);
                state->operations.recenter_query.enabled = false;
                recenter_update_target_data(state);
                ImGui::CloseCurrentPopup();
            }
            if (ImGui::IsItemHovered()) {
                md_bitfield_clear(&state->selection.highlight_mask);
                md_bitfield_copy(&state->selection.highlight_mask, &state->selection.selection_mask);
            }
        }

        ImGui::EndPopup();
    }
}

static void draw_selection_grow_window(ApplicationState* data) {
    ImGui::SetNextWindowSize(ImVec2(300,150), ImGuiCond_Always);
    if (ImGui::Begin("Selection Grow", &data->selection.grow.show_window, ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoCollapse)) {
        ImGui::PushItemWidth(-1);
        static uint64_t sel_popcount = 0;
        const uint64_t popcount = md_bitfield_popcount(&data->selection.selection_mask);
        const bool mode_changed = ImGui::Combo("##Mode", (int*)(&data->selection.grow.mode), "Covalent Bond\0Radial\0\0");
        const char* fmt = (data->selection.grow.mode == SelectionGrowthMode::CovalentBond) ? "%.0f" : "%.2f";
        const bool extent_changed = ImGui::SliderFloat("##Extent", &data->selection.grow.extent, 1.0f, 20.f, fmt);
        const bool appearing = ImGui::IsWindowAppearing();
        const bool sel_changed = popcount != sel_popcount;
        ImGui::PopItemWidth();

        const bool apply = ImGui::Button("Apply");

        // Need to invalidate when selection changes
        data->selection.grow.mask_invalid |= (mode_changed || extent_changed || appearing || sel_changed);

        if (data->selection.grow.mask_invalid) {
            sel_popcount = popcount;
            data->selection.grow.mask_invalid = false;
            md_bitfield_copy(&data->selection.grow.mask, &data->selection.selection_mask);

            switch (data->selection.grow.mode) {
            case SelectionGrowthMode::CovalentBond:
                md_util_mask_grow_by_bonds(&data->selection.grow.mask, &data->mold.sys, (int)data->selection.grow.extent, &data->representation.visibility_mask);
                break;
            case SelectionGrowthMode::Radial: {
                md_util_mask_grow_by_radius(&data->selection.grow.mask, &data->mold.state, data->selection.grow.extent, &data->representation.visibility_mask);
                break;
            }
            default:
                ASSERT(false);
            }

            grow_mask_by_selection_granularity(&data->selection.grow.mask, data->selection.granularity, data->mold.sys);
        }

        const bool show_preview =   (ImGui::GetHoveredID() == ImGui::GetID("##Extent")) ||
                                    (ImGui::GetActiveID()  == ImGui::GetID("##Extent")) ||
                                    (ImGui::GetHoveredID() == ImGui::GetID("Apply"));

        if (show_preview) {
            md_bitfield_copy(&data->selection.highlight_mask, &data->selection.grow.mask);
        }
        if (apply) {
            md_bitfield_copy(&data->selection.selection_mask, &data->selection.grow.mask);
            data->selection.grow.mask_invalid = true;
        }
    }
    ImGui::End();
}

static void draw_selection_query_window(ApplicationState* data) {
    ImGui::SetNextWindowSize(ImVec2(300,100), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Selection Query", &data->selection.query.show_window, ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoCollapse)) {

        if (ImGui::IsKeyPressed(ImGuiKey_Escape, false)) {
            data->selection.query.show_window = false;
            return;
        }

        static double query_frame = 0.0;

        ImGui::PushItemWidth(-1);
        bool apply = ImGui::InputQuery("##query", data->selection.query.buf, sizeof(data->selection.query.buf), data->selection.query.query_ok, data->selection.query.error, ImGuiInputTextFlags_AutoSelectAll | ImGuiInputTextFlags_EnterReturnsTrue);
        ImGui::PopItemWidth();

        if (ImGui::IsItemEdited() || data->animation.frame != query_frame) {
            data->selection.query.query_invalid = true;
        }
        bool preview = ImGui::IsItemFocused() || ImGui::IsItemHovered();

        if (ImGui::IsWindowAppearing()) {
            ImGui::SetKeyboardFocusHere(-1);
        }

        if (!data->selection.query.query_ok) ImGui::PushDisabled();
        apply |= ImGui::Button("Apply");
        if (!data->selection.query.query_ok) ImGui::PopDisabled();

        preview |= ImGui::IsItemHovered();

        if (data->selection.query.query_invalid) {
            data->selection.query.query_invalid = false;
            data->selection.query.query_ok = md_filter(&data->selection.query.mask, str_from_cstr(data->selection.query.buf), &data->mold.sys, &data->mold.state, data->script.ir, NULL, data->selection.query.error, sizeof(data->selection.query.error));
            query_frame = data->animation.frame;

            if (data->selection.query.query_ok) {
                grow_mask_by_selection_granularity(&data->selection.query.mask, data->selection.granularity, data->mold.sys);
            } else {
                md_bitfield_clear(&data->selection.query.mask);
            }
        }

        if (preview) {
            md_bitfield_copy(&data->selection.highlight_mask, &data->selection.query.mask);
        }

        if (apply && data->selection.query.query_ok) {
            md_bitfield_copy(&data->selection.selection_mask, &data->selection.query.mask);
            data->selection.query.show_window = false;
        }
    }
    ImGui::End();
}

static void draw_animation_window(ApplicationState* data) {
    ASSERT(data);
    size_t num_frames = run_num_frames(data);
    if (num_frames == 0) return;

    ASSERT(data->timeline.x_values);
    ASSERT(md_array_size(data->timeline.x_values) == num_frames);

    ImGui::SetNextWindowSize({300,200}, ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Animation", &data->animation.show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
        // The recording sets the frame itself
        ImGui::BeginDisabled(data->movie.state == MovieRecordingState::Recording);
        ImGui::Text("Num Frames: %zu", num_frames);
        md_unit_t time_unit = data->timeline.time_unit;
        double t   = frame_to_time(data->animation.frame, *data);
        double min = data->timeline.x_values[0];
        double max = data->timeline.x_values[num_frames - 1];
        char time_label[64];
        if (md_unit_is_none(time_unit)) {
            snprintf(time_label, sizeof(time_label), "Time");
        } else {
            char unit_buf[32];
            md_unit_print(unit_buf, sizeof(unit_buf), time_unit);
            snprintf(time_label, sizeof(time_label), "Time (%s)", unit_buf);
        }
        const float w = ImGui::CalcTextSize("Time (ps)").x;
        const float item_width = MAX(ImGui::GetContentRegionAvail().x - w, 100.f);
        ImGui::PushItemWidth(item_width);
        if (ImGui::BeginCombo("Interp.", interpolation_mode_str[(int)data->animation.interpolation])) {
            for (int i = 0; i < (int)InterpolationMode::Count; ++i) {
                if (ImGui::Selectable(interpolation_mode_str[i], (int)data->animation.interpolation == i)) {
                    data->animation.interpolation = (InterpolationMode)i;
					data->mold.interpolate_system_state = true;
                    data->mold.dirty_gpu_buffers |= MolBit_ClearVelocity;
                }
            }
            ImGui::EndCombo();
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Interpolation Mode for Atom Positions");
        }
        if (ImGui::SliderScalar(time_label, ImGuiDataType_Double, &t, &min, &max, "%.2f")) {
            data->animation.frame = time_to_frame(t, data->timeline.x_values);
        }
        ImGui::SliderFloat("Speed", &data->animation.fps, -200.0f, 200.f, "%.2f", ImGuiSliderFlags_Logarithmic);
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Animation Speed in Frames Per Second");
        }
        if (data->animation.interpolation == InterpolationMode::CubicSpline) {
            ImGui::SliderFloat("Tension", &data->animation.tension, 0.0f, 1.0f, "%.2f");
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Tension of the Cubic Spline (0 = Catmull-Rom, 1 = ease in and out of every frame)");
            }
        }
        switch (data->animation.mode) {
            case PlaybackMode::Playing:
                if (ImGui::Button((const char*)ICON_FA_PAUSE)) {
                    data->animation.mode = PlaybackMode::Stopped;
                }
                break;
            case PlaybackMode::Stopped:
                if (ImGui::Button((const char*)ICON_FA_PLAY)) {
                    data->animation.mode = PlaybackMode::Playing;
                }
                break;
            default:
                ASSERT(false);
        }
        ImGui::SameLine();
        if (ImGui::Button((const char*)ICON_FA_STOP)) {
            data->animation.mode = PlaybackMode::Stopped;
            data->animation.frame = 0.0;
        }
        ImGui::PopItemWidth();
        ImGui::EndDisabled();
    }
    ImGui::End();
}

// How the isosurfaces of an electronic structure representation are coloured. The atoms' colours
// are not offered for a density property; a field is offered whenever the dataset has what one is
// made from.
static bool draw_surface_coloring(ElectronicStructureRepresentation& es, const md_system_t& sys, bool advanced, bool allow_atom_colors) {
    bool update_rep = false;

    bool any_field = false;
    for (int k = 0; k < (int)SurfaceFieldKind::Count; ++k) {
        any_field |= surface_field_available((SurfaceFieldKind)k, sys);
    }

    if (ImGui::BeginCombo("coloring", surface_coloring_str[(int)es.coloring])) {
        for (int i = 0; i < (int)SurfaceColoring::Count; ++i) {
            const SurfaceColoring c = (SurfaceColoring)i;
            const bool enabled = (c == SurfaceColoring::AtomColors) ? allow_atom_colors :
                                 (c == SurfaceColoring::Field)      ? any_field : true;
            ImGui::BeginDisabled(!enabled);
            if (ImGui::Selectable(surface_coloring_str[i], es.coloring == c)) {
                es.coloring = c;
                // A field this dataset can make, if the one the representation names is not
                if (c == SurfaceColoring::Field && !surface_field_available(es.field_kind, sys)) {
                    es.field_kind = surface_field_first_available(sys);
                }
                update_rep = true;
            }
            ImGui::EndDisabled();
            if (!enabled && c == SurfaceColoring::Field) {
                ImGui::SetItemTooltip("Nothing in this dataset to make a field from:\n"
                                      "the electrostatic potential needs a basis, a ground state density and the QM atoms' nuclear charges,\n"
                                      "the embedding potential classical charges (atom/charge) on atoms outside the QM region");
            }
        }
        ImGui::EndCombo();
    }

    if (es.coloring == SurfaceColoring::AtomColors && advanced) {
        const double min_power = 2.0;
        const double max_power = 20.0;
        update_rep |= ImGui::SliderScalar("gaussian power", ImGuiDataType_Double, &es.gaussian_splatting_power, &min_power, &max_power, "%.2f");
    }

    if (es.coloring != SurfaceColoring::Field) {
        return update_rep;
    }

    if (ImGui::BeginCombo("field", surface_field_kind_str[(int)es.field_kind])) {
        for (int k = 0; k < (int)SurfaceFieldKind::Count; ++k) {
            const SurfaceFieldKind kind = (SurfaceFieldKind)k;
            ImGui::BeginDisabled(!surface_field_available(kind, sys));
            if (ImGui::Selectable(surface_field_kind_str[k], es.field_kind == kind)) {
                es.field_kind = kind;
                update_rep = true;
            }
            ImGui::EndDisabled();
        }
        ImGui::EndCombo();
    }
    if (es.field_kind == SurfaceFieldKind::EmbeddingPotential) {
        ImGui::SetItemTooltip("The electrostatic potential of the environment's classical multipoles on the surface: the charges,\n"
                              "and the dipoles and quadrupoles where the potential has them.\n"
                              "The QM region is left out, whatever charges a topology gives its atoms.\n"
                              "For a polarizable embedding this is the potential of its PERMANENT multipoles:\n"
                              "the induced dipoles are solved for during the calculation and not stored.");
    } else if (es.field_kind == SurfaceFieldKind::ElectrostaticPotential) {
        ImGui::SetItemTooltip("The electrostatic potential of the QM region on the surface: its nuclei and its ground state\n"
                              "electron density, exact over the basis. Negative (red) where the electrons dominate - lone pairs,\n"
                              "pi systems, where an electrophile is drawn - and positive (blue) where the nuclei are less screened.\n"
                              "Mapped onto an isodensity surface (0.001-0.002 au) it is the molecular electrostatic potential.\n"
                              "An embedding's environment is not part of it.");
    }

    // The colours are applied per pixel as the surface is shaded, so a change to the scale needs
    // nothing re-evaluated: it is not an update of the representation
    const SurfaceFieldVolume& vol = es.field_vol;
    color_scale_draw_controls(&es.field_map, surface_field_span(vol), surface_field_unit(es.field_kind),
                              "on surface", "its 1st to 99th percentile on the surface");
    if (advanced && vol.num_surface_samples > 0) {
        ImGui::TextDisabled("%zu of %zu voxels (%dx%dx%d) evaluated on the %s", vol.num_evaluated, vol.num_voxels,
                            vol.grid.dim[0], vol.grid.dim[1], vol.grid.dim[2], vol.on_gpu ? "GPU" : "CPU");
    }

    return update_rep;
}

static bool draw_representations_window_electronic_structure(ApplicationState* state, Representation& rep) {
    ImGuiComboFlags flags = 0;
    ElectronicStructureRepresentation& es = rep.electronic_structure;

    // The key is an opaque hash now, so it is no longer worth showing beside the label: a property
    // with no label is named by its path's last segment already, which is the readable thing.
    auto density_property_combo_label = [](char* buf, size_t buf_cap, const DensityProperty* prop) {
        ASSERT(buf);
        ASSERT(buf_cap > 0);
        if (!prop || str_empty(prop->label)) {
            snprintf(buf, buf_cap, "None");
            return;
        }
        snprintf(buf, buf_cap, "%.*s", (int)prop->label.len, prop->label.ptr);
    };

    bool advanced = state->representation.advanced_mode;
    bool update_rep = false;

    const ElectronicStructureSourceFlags source_mask = es_source_mask(state->mold.sys);

    if (ImGui::BeginCombo("volume src", electronic_structure_source_str[(int)es.source], flags)) {
        for (int n = 0; n < (int)ElectronicStructureSource::Count; n++) {
            ElectronicStructureSource source = (ElectronicStructureSource)n;
            bool is_selected = (es.source == source);
            bool disabled = !electronic_structure_source_supported(source_mask, source);

            if (disabled) ImGui::PushDisabled();
            if (ImGui::Selectable(electronic_structure_source_str[n], is_selected)) {
                es.source = source;
                electronic_structure_set_source_defaults(&es);
                update_rep = true;
            }
            if (disabled) ImGui::PopDisabled();

            if (is_selected) {
                ImGui::SetItemDefaultFocus();
            }
        }
        ImGui::EndCombo();
    }
    if (advanced) {
        if (ImGui::Combo("volume res", (int*)&es.resolution, volume_resolution_str, IM_ARRAYSIZE(volume_resolution_str))) {
            update_rep = true;
        }
    }

    if (electronic_structure_is_density_property(es)) {
        // Gathered here rather than kept in a list somewhere: it is a handful of entries borrowed
        // from the attribute table, and a stored copy would be one more thing to invalidate on load.
        DensityProperty density_props[64];
        const int num_density_props = (int)MIN(density_properties_gather(density_props, ARRAY_SIZE(density_props), state->mold.sys), ARRAY_SIZE(density_props));
        if (num_density_props > 0) {
            // A key that names nothing in the current list - a fresh representation, or a workspace
            // saved against a file whose properties differ - falls back to the first one rather
            // than drawing nothing.
            const DensityProperty* selected = nullptr;
            for (int n = 0; n < num_density_props; ++n) {
                if (density_props[n].key == es.density_property_key) {
                    selected = &density_props[n];
                    break;
                }
            }
            if (!selected) {
                selected = &density_props[0];
                es.density_property_key = selected->key;
                update_rep = true;
            }

            char preview[256];
            density_property_combo_label(preview, sizeof(preview), selected);
            if (ImGui::BeginCombo("property", preview)) {
                for (int n = 0; n < num_density_props; ++n) {
                    char label[256];
                    density_property_combo_label(label, sizeof(label), &density_props[n]);
                    const bool is_selected = (es.density_property_key == density_props[n].key);
                    if (ImGui::Selectable(label, is_selected)) {
                        if (es.density_property_key != density_props[n].key) {
                            es.density_property_key = density_props[n].key;
                            update_rep = true;
                        }
                    }

                    if (is_selected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }
        } else {
            ImGui::LabelText("property", "No density properties available");
        }
    }

    const bool show_spin = electronic_structure_uses_spin(es);
    if (show_spin) {
        // Asked of the table: a second, DISTINCT set of coefficients. This used to be true for a
        // restricted file too, because beta is published there as an alias of alpha and the fan-in
        // only checked that the column existed - so the selector appeared and toggled between two
        // names for one datum.
        const bool has_beta = es_has_distinct_beta_orbitals(state->mold.sys);
        if (es.source == ElectronicStructureSource::MolecularOrbital) {
            if (es.spin != ElectronicStructureSpin::Alpha && es.spin != ElectronicStructureSpin::Beta) {
                es.spin = ElectronicStructureSpin::Alpha;
                update_rep = true;
            }
            if (has_beta) {
                const char* spin_options[] = { "Alpha", "Beta" };
                int spin_idx = es.spin == ElectronicStructureSpin::Beta ? 1 : 0;
                if (ImGui::Combo("spin", &spin_idx, spin_options, IM_ARRAYSIZE(spin_options))) {
                    es.spin = spin_idx == 1 ? ElectronicStructureSpin::Beta : ElectronicStructureSpin::Alpha;
                    update_rep = true;
                }
            }
        }
        else if (es.source == ElectronicStructureSource::ElectronDensity) {
            if (has_beta) {
                const ElectronicStructureSpin spin_options[] = {
                    ElectronicStructureSpin::Total,
                    ElectronicStructureSpin::Alpha,
                    ElectronicStructureSpin::Beta,
                    ElectronicStructureSpin::Difference,
                };
                const char* spin_labels[] = { "Total", "Alpha", "Beta", "Difference" };
                int spin_idx = 0;
                for (int i = 0; i < (int)IM_ARRAYSIZE(spin_options); ++i) {
                    if (es.spin == spin_options[i]) {
                        spin_idx = i;
                        break;
                    }
                }
                if (es.spin != spin_options[spin_idx]) {
                    es.spin = spin_options[spin_idx];
                    update_rep = true;
                }
                if (ImGui::Combo("spin", &spin_idx, spin_labels, IM_ARRAYSIZE(spin_labels))) {
                    es.spin = spin_options[spin_idx];
                    if (es.spin != ElectronicStructureSpin::Difference) {
                        es.use_magnitude = false;
                    }
                    update_rep = true;
                }
            } else if (es.spin != ElectronicStructureSpin::Total) {
                es.spin = ElectronicStructureSpin::Total;
                es.use_magnitude = false;
                update_rep = true;
            }
        }
    }

    if (es.source == ElectronicStructureSource::NaturalTransitionOrbital) {
        if (ImGui::Combo("component", (int*)&es.nto_component, electronic_structure_nto_component_str, IM_ARRAYSIZE(electronic_structure_nto_component_str))) {
            update_rep = true;
        }
    }

    if (es.source == ElectronicStructureSource::TransitionDensity) {
        if (ImGui::Combo("component", (int*)&es.transition_density_component, electronic_structure_transition_density_component_str, IM_ARRAYSIZE(electronic_structure_transition_density_component_str))) {
            update_rep = true;
        }
    }

    const bool show_molecular_orbitals = electronic_structure_uses_orbital_idx(es);
    const bool show_exited_states = electronic_structure_uses_excited_state_idx(es);
    const bool show_lambdas = electronic_structure_uses_nto_lambda_idx(es);

    if (show_molecular_orbitals) {
        // The (homo)/(lumo) markers sit at different indices for the two spins in an unrestricted
        // calculation, and the rendered volume already follows es.spin - so the frontier has to
        // follow it too, or the entry marked (homo) under spin=Beta is alpha's HOMO index.
        //
        // Formatted here rather than read out of a precomputed list: it is an index, an energy and
        // two markers, and a copy of it would be one more thing to keep in step with the table.
        const bool use_beta = (es.spin == ElectronicStructureSpin::Beta) && es_has_distinct_beta_orbitals(state->mold.sys);

        size_t num_orbitals = 0;
        OrbitalFrontier frontier = {};
        if (es_orbital_extent(state->mold.sys, &num_orbitals, nullptr) && num_orbitals > 0) {
            es_orbital_frontier(&frontier, state->mold.sys, use_beta ? es_path::beta_occupation : es_path::alpha_occupation);

            // Energies alongside the index, so the list can be read as an orbital diagram. A column
            // that does not span the orbital set is not one to index with, so it is dropped whole
            // rather than bounds checked per row - the list then falls back to bare indices.
            size_t num_energies = 0;
            const double* energy = es_orbital_energies(&num_energies, state->mold.sys, use_beta ? es_path::beta_energy : es_path::alpha_energy);
            if (num_energies != num_orbitals) energy = nullptr;

            // The UI font is monospaced (Dejavu Sans Mono), so a fixed field width per column is all
            // the alignment this needs: right aligned index, then a fixed number of decimals, which
            // puts every decimal point of the list on one column. Widths are passed in rather than
            // baked into the format so the preview - a single line, nothing to align against - can
            // ask for the same text unpadded.
            int idx_width = 1;
            for (size_t m = num_orbitals; m >= 10; m /= 10) idx_width++;

            auto orbital_label = [&](char* buf, size_t cap, int n, int idx_w, int ene_w) {
                const char* mark = (n == frontier.homo_idx) ? " (homo)" : (n == frontier.lumo_idx) ? " (lumo)" : "";
                if (energy) {
                    snprintf(buf, cap, "%*i  %*.4f%s", idx_w, n + 1, ene_w, energy[n], mark);
                } else {
                    snprintf(buf, cap, "%*i%s", idx_w, n + 1, mark);
                }
            };

            es.orbital_idx = CLAMP(es.orbital_idx, 0, (int)num_orbitals - 1);

            char preview[48];
            orbital_label(preview, sizeof(preview), es.orbital_idx, 0, 0);
            if (ImGui::BeginCombo("orbital idx", preview)) {
                // Highest index first, so the lowest lying orbital sits at the bottom - the
                // conventional orbital diagram, and the same order the Molecular Orbitals table
                // defaults to. Canonical MOs come out of the SCF in ascending energy, so descending
                // index IS descending energy; sorting on the energies instead would only differ for
                // a non-aufbau set, and would then disagree with that table.
                for (int n = (int)num_orbitals - 1; n >= 0; n--) {
                    bool is_selected = (es.orbital_idx == n);
                    char label[48];
                    orbital_label(label, sizeof(label), n, idx_width, 10);
                    if (ImGui::Selectable(label, is_selected)) {
                        if (es.orbital_idx != n) {
                            update_rep = true;
                        }
                        es.orbital_idx = n;
                    }

                    if (is_selected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }
        }
    }

    if (show_exited_states) {
        const int num_excited_states = (int)es_excited_state_count(state->mold.sys);
        if (num_excited_states > 0) {
            es.excited_state_idx = CLAMP(es.excited_state_idx, 0, num_excited_states - 1);
            char preview[16];
            snprintf(preview, sizeof(preview), "%i", es.excited_state_idx + 1);
            if (ImGui::BeginCombo("state idx", preview)) {
                for (int n = 0; n < num_excited_states; n++) {
                    const bool is_selected = (es.excited_state_idx == n);
                    char label[16];
                    snprintf(label, sizeof(label), "%i", n + 1);
                    if (ImGui::Selectable(label, is_selected)) {
                        if (es.excited_state_idx != n) {
                            update_rep = true;
                        }
                        es.excited_state_idx = n;
                    }

                    if (is_selected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }
            if (show_lambdas) {
                // The same 1e-3 the rest of the tree stops at, and the same place the zero padded
                // tail of the ragged lambda axis stops it.
                const double LAMBDA_CUTOFF = 1.0e-3;
                double lambdas[32];
                const int num_lambdas = (int)es_nto_lambdas(lambdas, ARRAY_SIZE(lambdas), state->mold.sys,
                                                            (size_t)es.excited_state_idx, LAMBDA_CUTOFF);
                if (num_lambdas > 0) {
                    es.nto_lambda_idx = CLAMP(es.nto_lambda_idx, 0, num_lambdas - 1);

                    auto lambda_label = [&lambdas](char* buf, size_t cap, int n) {
                        // The bytes rather than a u8 literal: u8"" is const char8_t* in C++20 and
                        // a format string has to be a plain char array. \xce\xbb is U+03BB, lambda.
                        snprintf(buf, cap, "\xce\xbb[%i] (%.3f)", n + 1, lambdas[n]);
                    };

                    char preview[32];
                    lambda_label(preview, sizeof(preview), es.nto_lambda_idx);
                    if (ImGui::BeginCombo("lambda idx", preview)) {
                        for (int n = 0; n < num_lambdas; n++) {
                            const bool is_selected = (es.nto_lambda_idx == n);
                            char label[32];
                            lambda_label(label, sizeof(label), n);
                            if (ImGui::Selectable(label, is_selected)) {
                                if (es.nto_lambda_idx != n) {
                                    update_rep = true;
                                }
                                es.nto_lambda_idx = n;
                            }

                            if (is_selected) {
                                ImGui::SetItemDefaultFocus();
                            }
                        }
                        ImGui::EndCombo();
                    }
                }
            }
        }
    }
    if (electronic_structure_uses_magnitude_toggle(es)) {
        const char* magnitude_label = es.source == ElectronicStructureSource::ElectronDensity ? (const char*)u8"magnitude |ρ|" : (const char*)u8"magnitude |Ψ|";
        if (ImGui::Checkbox(magnitude_label, &es.use_magnitude)) {
            update_rep = true;
        }
    }
    
    const double min_tau = 0.0;
    const double max_tau = 1.0;
    
    const double iso_min = 1.0e-2;
    const double iso_max = 5.0;

    if (electronic_structure_is_density_property(es)) {
        int num_isos = CLAMP(es.density_property.num_isos, 1, (int)ARRAY_SIZE(es.density_property.values));
        if (num_isos != es.density_property.num_isos) {
            es.density_property.num_isos = num_isos;
            update_rep = true;
        }

        if (ImGui::SliderInt("iso count", &num_isos, 1, (int)ARRAY_SIZE(es.density_property.values))) {
            es.density_property.num_isos = num_isos;
            update_rep = true;
        }

        for (int i = 0; i < es.density_property.num_isos; ++i) {
            char value_label[32];
            char color_label[32];
            snprintf(value_label, sizeof(value_label), "iso %d value", i + 1);
            snprintf(color_label, sizeof(color_label), "iso %d color", i + 1);
            update_rep |= ImGui::DragScalar(value_label, ImGuiDataType_Double, &es.density_property.values[i], 0.001f, NULL, NULL, "%.8f");
            update_rep |= ImGui::ColorEdit4(color_label, es.density_property.colors[i].elem);
        }

        if (advanced) {
            ImGui::SliderScalar((const char*)u8"iso τ", ImGuiDataType_Double, &es.iso_optical_density, &min_tau, &max_tau, "%.4f", ImGuiSliderFlags_Logarithmic);
            ImGui::SetItemTooltip("Optical density shared by all custom isosurfaces");
        }

        update_rep |= draw_surface_coloring(es, state->mold.sys, advanced, false);

        return update_rep;
    }
    
    const char* iso_label = electronic_structure_iso_value_label();
    
    if (electronic_structure_is_signed(es)) {
        // The band a field is evaluated in follows the isovalue
        update_rep |= ImGui::SliderScalar(iso_label, ImGuiDataType_Double, &rep.electronic_structure.iso_value, &iso_min, &iso_max, "%.8f", ImGuiSliderFlags_Logarithmic) && rep.electronic_structure.coloring == SurfaceColoring::Field;
        ImGui::SetItemTooltip("%s", electronic_structure_iso_value_tooltip(rep.electronic_structure));
        if (advanced) {
            ImGui::SliderScalar((const char*)u8"iso τ", ImGuiDataType_Double, &rep.electronic_structure.iso_optical_density, &min_tau, &max_tau, "%.4f", ImGuiSliderFlags_Logarithmic);
            ImGui::SetItemTooltip("Optical density of the isosurfaces");
        }
        update_rep |= draw_surface_coloring(rep.electronic_structure, state->mold.sys, advanced, true);
        if (rep.electronic_structure.coloring != SurfaceColoring::Uniform) {
            ImGui::ColorEdit4("tint positive", rep.electronic_structure.tint_psi_pos.elem);
            ImGui::ColorEdit4("tint negative", rep.electronic_structure.tint_psi_neg.elem);
        } else {
            ImGui::ColorEdit4("color positive", rep.electronic_structure.col_psi_pos.elem);
            ImGui::ColorEdit4("color negative", rep.electronic_structure.col_psi_neg.elem);
        }
    }
    else {
        // The band a field is evaluated in follows the isovalue
        update_rep |= ImGui::SliderScalar(iso_label, ImGuiDataType_Double, &rep.electronic_structure.iso_value, &iso_min, &iso_max, "%.8f", ImGuiSliderFlags_Logarithmic) && rep.electronic_structure.coloring == SurfaceColoring::Field;
        ImGui::SetItemTooltip("%s", electronic_structure_iso_value_tooltip(rep.electronic_structure));
        if (advanced) {
            ImGui::SliderScalar((const char*)u8"iso τ", ImGuiDataType_Double, &rep.electronic_structure.iso_optical_density, &min_tau, &max_tau, "%.4f", ImGuiSliderFlags_Logarithmic);
            ImGui::SetItemTooltip("Optical density of the isosurfaces");
        }
        update_rep |= draw_surface_coloring(rep.electronic_structure, state->mold.sys, advanced, true);
        if (rep.electronic_structure.coloring != SurfaceColoring::Uniform) {
            if (es.source == ElectronicStructureSource::TransitionDensity && es.transition_density_component == ElectronicStructureTransitionDensityComponent::Attachment) {
                ImGui::ColorEdit4("tint attachment", rep.electronic_structure.tint_att.elem);
            }
            else if (es.source == ElectronicStructureSource::TransitionDensity && es.transition_density_component == ElectronicStructureTransitionDensityComponent::Detachment) {
                ImGui::ColorEdit4("tint detachment", rep.electronic_structure.tint_det.elem);
            } else {
                ImGui::ColorEdit4("tint density", rep.electronic_structure.tint_den.elem);
            }
        } else {
            if (es.source == ElectronicStructureSource::TransitionDensity && es.transition_density_component == ElectronicStructureTransitionDensityComponent::Attachment) {
                ImGui::ColorEdit4("color attachment", rep.electronic_structure.col_att.elem);
            }
            else if (es.source == ElectronicStructureSource::TransitionDensity && es.transition_density_component == ElectronicStructureTransitionDensityComponent::Detachment) {
                ImGui::ColorEdit4("color detachment", rep.electronic_structure.col_det.elem);
            } else {
                ImGui::ColorEdit4("color density",  rep.electronic_structure.col_den.elem);
            }
        }
    }

    return update_rep;
}

static void draw_representations_window(ApplicationState* state) {
    if (!state->representation.show_window) return;

    ImGui::SetNextWindowSize({300,200}, ImGuiCond_FirstUseEver);
    ImGui::Begin("Representations", &state->representation.show_window, ImGuiWindowFlags_NoFocusOnAppearing);
    // The frames of a recording are made from the representations as they are, so they are not edited meanwhile
    const bool movie_locked = state->movie.state == MovieRecordingState::Recording;
    ImGui::BeginDisabled(movie_locked);
    if (ImGui::Button("create new")) {
        create_representation(state);
    }
    ImGui::SameLine();
    if (ImGui::DeleteButton("remove all")) {
        remove_all_representations(state);
    }
    ImGui::SameLine();

    const char* advanced_label = ICON_FA_SLIDERS;
    const float checkbox_width =
        ImGui::GetFrameHeight() +
        ImGui::GetStyle().ItemInnerSpacing.x +
        ImGui::CalcTextSize(advanced_label).x;
    ImGui::SetCursorPosX(ImGui::GetWindowContentRegionMax().x - checkbox_width);
    ImGui::Checkbox(advanced_label, &state->representation.advanced_mode);
    ImGui::SetItemTooltip("Advanced mode (show all representation settings)");

    bool advanced = state->representation.advanced_mode;

    ImGui::Spacing();
    ImGui::Separator();
    for (int rep_idx = 0; rep_idx < (int)md_array_size(state->representation.reps); rep_idx++) {
        bool update_rep = false;
        Representation& rep = state->representation.reps[rep_idx];

        char label[128];
        snprintf(label, sizeof(label), "%s###ID", rep.name);

        ImGui::PushID(rep_idx);
        
        const float pad = 3.0f;
        const float size = ImGui::GetFontSize() + pad * 2;
        const float spacing = 2.f;
        const float total_button_size = (size + 1) * 3;

        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0, pad));
        bool draw_content = ImGui::TreeNodeEx("##label", ImGuiTreeNodeFlags_FramePadding);
        ImGui::PopStyleVar();

        ImGui::SameLine();
        ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x - total_button_size);
        ImGui::InputText("##name", rep.name, sizeof(rep.name));

        ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - total_button_size, spacing);
        const char* eye_icon = rep.enabled ? ICON_FA_EYE : ICON_FA_EYE_SLASH;
        
        const ImVec2 btn_size = {size, size};
        if (ImGui::Button(eye_icon, btn_size)) {
            rep.enabled = !rep.enabled;
            state->representation.atom_visibility_mask_dirty = true;
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Show/Hide");
        }
        ImGui::SameLine(0, spacing);
        if (ImGui::Button(ICON_FA_COPY, btn_size)) {
            clone_representation(state, rep);
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Duplicate");
        }
        ImGui::SameLine(0, spacing);
        if (ImGui::DeleteButton(ICON_FA_XMARK, btn_size)) {
            remove_representation(state, rep_idx);
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Remove");
        }

        if (draw_content) {
            ImVec4 sub_group_color = ImGui::GetStyleColorVec4(ImGuiCol_Text);
            sub_group_color.w = 0.5f;

            const float inner_item_width = MAX(ImGui::GetContentRegionAvail().x - ImGui::GetStyle().IndentSpacing - ImGui::CalcTextSize("helix scale").x, 100.f);
            ImGui::PushItemWidth(inner_item_width);

            // @TODO: Only display the representations which can be used for the current dataset
            if (!rep.type_is_valid) ImGui::PushInvalid();
            if (ImGui::BeginCombo("type", representation_type_str[(int)rep.type])) {
                for (int i = 0; i < (int)RepresentationType::Count; ++i) {
                    if (i == (int)RepresentationType::ElectronicStructure) {
                        // Enlisted when the system carries ANY electronic structure source - orbitals are
                        // one, and a file without an SCF block can still carry density properties
                        if (es_source_mask(state->mold.sys) == 0) continue;
                    }
                    if (ImGui::Selectable(representation_type_str[(int)i], i == (int)rep.type)) {
                        rep.type = (RepresentationType)i;
                        if (rep.type == RepresentationType::ElectronicStructure) {
                            // The default source is orbitals, which this system may not have
                            electronic_structure_select_available_source(&rep.electronic_structure, state->mold.sys);
                        }
                        update_rep = true;
                    }
                }
                ImGui::EndCombo();
            }
            if (!rep.type_is_valid) ImGui::PopInvalid();

			const double dipole_min = 0.01;
			const double dipole_max = 100.0;

            switch (rep.type) {
            case RepresentationType::ElectronicStructure:
                update_rep |= draw_representations_window_electronic_structure(state, rep);
                break;
            case RepresentationType::DipoleMoment:
            {
                // The combo lists GROUPS, which is what the file publishes. A group holding one
                // moment per excited state is then picked apart by the index widget below, rather
                // than flattening every (group, element) pair into one list nobody can scan.
                DipoleGroup groups[16];
                const size_t num_groups = MIN(dipole_groups_gather(groups, ARRAY_SIZE(groups), state->mold.sys), ARRAY_SIZE(groups));

                // Asked of the key the representation holds, not searched for in the list: a key
                // that names nothing - a workspace pointed at a file without that group - simply
                // leaves the preview empty and the index widget away.
                DipoleGroup selected = {};
                const bool has_selection = dipole_group_from_key(&selected, state->mold.sys, rep.dipole.dipole_key);

                // A group can come back shorter than it was (fewer excited states), which would
                // otherwise leave the representation pointing past its end.
                if (has_selection && rep.dipole.dipole_index >= selected.count) {
                    rep.dipole.dipole_index = 0;
                    update_rep = true;
                }

                char preview[64] = "";
                if (has_selection) dipole_label_pretty(preview, sizeof(preview), selected.label);

                if (ImGui::BeginCombo("dipole", preview)) {
                    for (size_t i = 0; i < num_groups; ++i) {
                        char lbl[64];
                        dipole_label_pretty(lbl, sizeof(lbl), groups[i].label);
                        const bool is_selected = groups[i].key == rep.dipole.dipole_key;
                        if (ImGui::Selectable(lbl, is_selected)) {
                            if (!is_selected) {
                                rep.dipole.dipole_key   = groups[i].key;
                                rep.dipole.dipole_index = 0;   // a different group; the old element means nothing in it
                                update_rep = true;
                            }
                        }
                        if (is_selected) {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }

                // Only an array of moments has an index to choose. One moment is addressed by the
                // group alone, and a slider from 1 to 1 is a control that cannot do anything.
                if (has_selection && selected.count > 1) {
                    int index = (int)rep.dipole.dipole_index + 1;
                    if (ImGui::SliderInt("index", &index, 1, (int)selected.count)) {
                        rep.dipole.dipole_index = (uint32_t)CLAMP(index - 1, 0, (int)selected.count - 1);
                        update_rep = true;
                    }
                }
            }
            ImGui::ColorEdit4("color", rep.dipole.color.elem);
				ImGui::SliderScalar("scale", ImGuiDataType_Double, &rep.dipole.scale, &dipole_min, &dipole_max, "%.3f");
				ImGui::SliderFloat("radius", &rep.dipole.radius, 0.01f, 0.5f);
				ImGui::DragFloat3("offset", rep.dipole.offset.elem, 0.01f, -100.0f, 100.0f);
                break;
			default:
                if (ImGui::InputQuery("filter", rep.filt, sizeof(rep.filt), rep.filt_is_valid, rep.filt_error)) {
                    rep.filt_is_dirty = true;
                    update_rep = true;
                }
                break;
            }

            if (representation_uses_atom_colors(rep)) {
                if (ImGui::Combo("atom color", (int*)(&rep.color_mapping), color_mapping_str, IM_ARRAYSIZE(color_mapping_str))) {
                    update_rep = true;
                }

                if (rep.color_mapping == ColorMapping::Attribute) {
                    // The list of per atom fields IS the system's attribute table under atom/.
                    // Queried here rather than cached anywhere, so it cannot disagree with the data.
                    md_attribute_id_t prop_ids[64];
                    size_t num_props = MIN(atom_attribute_query(prop_ids, ARRAY_SIZE(prop_ids), state->mold.sys), ARRAY_SIZE(prop_ids));

                    const md_attributes_t& attributes = state->mold.sys.attributes;
                    const md_attribute_t* selected_prop = md_attributes_get(&attributes, rep.atom_attribute.key);

                    // A key the table no longer holds - a reload which dropped that field - falls
                    // back to the first available rather than leaving the representation blank.
                    if (!selected_prop && num_props > 0) {
                        atom_attribute_select(&rep.atom_attribute, prop_ids[0], state->mold.sys);
                        selected_prop = md_attributes_get(&attributes, rep.atom_attribute.key);
                        update_rep = true;
                    }

                    if (num_props > 0 && selected_prop) {
                        if (ImGui::BeginCombo("attribute", atom_attribute_label(selected_prop).ptr)) {
                            for (size_t i = 0; i < num_props; ++i) {
                                const md_attribute_t* attr = md_attributes_get(&attributes, prop_ids[i]);
                                if (!attr) continue;
                                bool selected = prop_ids[i] == rep.atom_attribute.key;
                                if (ImGui::Selectable(atom_attribute_label(attr).ptr, selected)) {
                                    atom_attribute_select(&rep.atom_attribute, prop_ids[i], state->mold.sys);
                                    update_rep = true;
                                }
                            }
                            ImGui::EndCombo();
                        }

                        const int num_variants = atom_attribute_variant_count(selected_prop);
                        if (num_variants > 1) {
                            int idx = rep.atom_attribute.variant_idx + 1;
                            const int min = 1;
                            const int max = num_variants;
                            if (ImGui::SliderInt("index", &idx, min, max)) {
                                update_rep = true;
                            }
                            rep.atom_attribute.variant_idx = CLAMP(idx - 1, 0, num_variants - 1);
                        }
                        
                        // The same scale an isosurface coloured by a field has. Its span is what
                        // the colours were last made from, the atoms of this representation.
                        const bool multiple = num_variants > 1;
                        update_rep |= color_scale_draw_controls(&rep.atom_attribute.scale, rep.atom_attribute.span, selected_prop->unit, "shown atoms",
                                                                multiple ? "the smallest to the largest value of the atoms shown, at every index"
                                                                         : "the smallest to the largest value of the atoms shown");
                    } else {
                        ImGui::TextDisabled("no per atom attributes in this dataset");
                    }
                }
                if (rep.filt_is_dynamic || rep.color_mapping == ColorMapping::Attribute) {
                    if (advanced) {
                        update_rep |= ImGui::Checkbox("auto-update", &rep.dynamic_evaluation);
                        if (!rep.dynamic_evaluation) {
                            ImGui::SameLine();
                            if (ImGui::Button("update")) {
                                rep.filt_is_dirty = true;
                                update_rep = true;
                            }
                        }
                    }
                } else {
                    rep.dynamic_evaluation = false;
                }

                if (rep.color_mapping == ColorMapping::Uniform) {
                    update_rep |= ImGui::ColorEdit3("##atom_base_color", rep.base_color.elem);
                }

                if (advanced) {
                    if (rep.type == RepresentationType::Licorice || rep.type == RepresentationType::BallAndStick) {
                        // Draw options for how bonds should be colored
					    update_rep |= ImGui::Combo("bond color", (int*)(&rep.bond_color), bond_color_mode_str, IM_ARRAYSIZE(bond_color_mode_str));
                        if (rep.bond_color == BondColorMode::SmoothAtom) {
                            ImGui::SliderFloat("sharpness", &rep.bond_sharpness, 0.0f, 1.0f);
                        } else if (rep.bond_color == BondColorMode::Uniform) {
                            ImGui::ColorEdit3("##bond_base_color", rep.bond_base_color.elem);
					    }
                    }
                }
            }

            if (advanced && representation_uses_atom_colors(rep)) {
                if (rep.type != RepresentationType::ElectronicStructure) {
                    ImGui::Spacing();
                    ImGui::TextColored(sub_group_color, "Geometric Scaling");
                    switch (rep.type) {
                    case RepresentationType::SpaceFill:
                        update_rep |= ImGui::SliderFloat("radius scale", &rep.scale[0], 0.1f, 4.f);
                        break;
                    case RepresentationType::Licorice:
                        update_rep |= ImGui::SliderFloat("radius scale", &rep.scale[0], 0.1f, 4.0f);
                        break;
                    case RepresentationType::BallAndStick:
                        update_rep |= ImGui::SliderFloat("ball scale",   &rep.scale[0], 0.1f, 4.f);
                        update_rep |= ImGui::SliderFloat("bond scale",   &rep.scale[1], 0.1f, 4.f);
                        break;
                    case RepresentationType::Ribbons:
                        update_rep |= ImGui::SliderFloat("width",        &rep.scale[0], 0.1f, 3.f);
                        update_rep |= ImGui::SliderFloat("thickness",    &rep.scale[1], 0.1f, 3.f);
                        break;
                    case RepresentationType::Cartoon:
                        update_rep |= ImGui::SliderFloat("coil",   &rep.scale[0], 0.1f, 3.f);
                        update_rep |= ImGui::SliderFloat("sheet",  &rep.scale[1], 0.1f, 3.f);
                        update_rep |= ImGui::SliderFloat("helix",  &rep.scale[2], 0.1f, 3.f);
                        break;
                    default:
                        break;
                    }
                }
                if (rep.color_mapping == ColorMapping::SecondaryStructure) {
                    ImGui::Spacing();
                    ImGui::TextColored(sub_group_color, "Secondary Structure Color");
                    update_rep |= ImGui::ColorEdit3("unknown", rep.secondary_structure.color_unknown.elem);
                    update_rep |= ImGui::ColorEdit3("coil",    rep.secondary_structure.color_coil.elem);
                    update_rep |= ImGui::ColorEdit3("helix",   rep.secondary_structure.color_helix.elem);
                    update_rep |= ImGui::ColorEdit3("sheet",   rep.secondary_structure.color_sheet.elem);
                }
                ImGui::Spacing();

                ImGui::TextColored(sub_group_color, "Post-Processing");
                update_rep |= ImGui::ColorEdit3("tint color", rep.tint_color.elem, ImGuiColorEditFlags_PickerHueWheel);
                update_rep |= ImGui::SliderFloat("tint scale", &rep.tint_scale, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_AlwaysClamp);
                update_rep |= ImGui::SliderFloat("saturation", &rep.saturation, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_AlwaysClamp);
            }

            ImGui::PopItemWidth();
            ImGui::Spacing();
            ImGui::TreePop();
        }

        ImGui::PopID();

        if (update_rep) {
            flag_representation_as_dirty(&rep);
        }
    }

    ImGui::EndDisabled();
    ImGui::End();
}



static void draw_async_task_window(ApplicationState* data) {
    constexpr float WIDTH = 300.f;
    constexpr float MARGIN = 10.f;

    task_system::ID tasks[256];
    size_t num_tasks = task_system::pool_running_tasks(tasks, ARRAY_SIZE(tasks));
    bool any_task_label_visible = false;
    for (size_t i = 0; i < num_tasks; i++) {
        str_t label = task_system::task_label(tasks[i]);
        if (!label || label[0] == '\0' || (label[0] == '#' && label[1] == '#')) continue;
        any_task_label_visible = true;
    }
    
    if (any_task_label_visible) {
        ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(viewport->Pos + ImVec2(data->app.window.width - WIDTH - MARGIN,
                                                       ImGui::GetCurrentContext()->FontSize + ImGui::GetStyle().FramePadding.y * 2.f + MARGIN));
        ImGui::SetNextWindowSize(ImVec2(WIDTH, 0));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
        ImGui::Begin("##Async Info", 0,
                     ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar |
                     ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing);

        const float pad = 3.0f;
        const float size = ImGui::GetFontSize() + pad * 2;

        char buf[64];
        for (size_t i = 0; i < MIN(num_tasks, 8); i++) {
            const auto id = tasks[i];
            str_t label = task_system::task_label(id);
            float fract = task_system::task_fraction_complete(id);

            if (!label || label[0] == '\0' || (label[0] == '#' && label[1] == '#')) continue;

            snprintf(buf, sizeof(buf), "%.*s %.1f%%", (int)label.len, label.ptr, fract * 100.f);
            ImGui::ProgressBar(fract, ImVec2(ImGui::GetContentRegionAvail().x - (size + pad),0), buf);
            ImGui::SameLine();
            if (ImGui::DeleteButton((const char*)ICON_FA_XMARK, ImVec2(size, size))) {
                task_system::task_interrupt(id);
                if (id == data->tasks.evaluate_full) {
                    md_script_eval_interrupt(data->script.full_eval);
                }
                else if(id == data->tasks.evaluate_filt) {
                    md_script_eval_interrupt(data->script.filt_eval);
                }
            }
        }

        ImGui::End();
        ImGui::PopStyleColor();
    }
}

struct TimelineArgs {
    const char* lbl;
    uint32_t col;
    int plot_height;

    struct {
        int count;
        int dim_y;

        const float* x;
        const float* y;
        const float* y_mean;
        const float* y_var;
        const float* y_min;
        const float* y_max;

        float min_y;
        float max_y;
        str_t unit;
    } values;

    struct {
        double* beg;
        double* end;
        double  min;
        double  max;
    } view_range;

    struct {
        bool* is_dragging;
        bool* is_selecting;
    } input;

    struct {
        bool show;
        bool enabled;

        double* beg;
        double* end;

        double min;
        double max;
    } filter;

    struct {
        bool enabled;
        double min;
        double max;
    } value_filter;

    double* time;
};

struct TimePayload {
    const TimelineArgs* args;
    int y_idx;
};

static ImPlotPoint get_time_point(int index, void* user_data) {
    const TimePayload* payload = (const TimePayload*)user_data;
    const TimelineArgs* args = payload->args;
    return ImPlotPoint(args->values.x[index], args->values.y[index * args->values.dim_y + payload->y_idx]);
}

bool draw_property_timeline(const ApplicationState& data, const TimelineArgs& args) {
    const ImPlotAxisFlags axis_flags = ImPlotAxisFlags_NoSideSwitch | ImPlotAxisFlags_NoHighlight;
    const ImPlotAxisFlags axis_flags_x = axis_flags;
    const ImPlotAxisFlags axis_flags_y = axis_flags | ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_RangeFit | ImPlotAxisFlags_NoLabel |ImPlotAxisFlags_NoTickLabels;
    
    const ImPlotFlags flags = ImPlotFlags_NoBoxSelect | ImPlotFlags_NoFrame;
    
    const float pad_x = ImPlot::GetStyle().PlotPadding.x;
    ImPlot::PushStyleVar(ImPlotStyleVar_PlotPadding, ImVec2(pad_x, 0));
    defer { ImPlot::PopStyleVar(1); };

    if (ImPlot::BeginPlot("##Timeline", ImVec2(-1,args.plot_height), flags)) {
        ImPlot::SetupAxisLinks(ImAxis_X1, args.view_range.beg, args.view_range.end);
        ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, args.view_range.min, args.view_range.max);
        ImPlot::SetupAxes(0, 0, axis_flags_x, axis_flags_y);
        ImPlot::SetupFinish();

        bool active = ImGui::IsItemActive();
        bool print_timeline_tooltip = false;

        if (args.value_filter.enabled) {
            float* y_vals = (float*)md_alloc(frame_alloc, args.values.count * sizeof(float));
            for (int i = 0; i < args.values.count; ++i) {
                float val = args.values.y[i];
                y_vals[i] = (args.value_filter.min < val && val < args.value_filter.max) ? args.values.max_y : -FLT_MAX;
            }

            const ImVec4 filter_frame_color = ImVec4(1, 1, 1, 1);
            ImPlot::SetNextFillStyle(filter_frame_color, 0.15f);
            ImPlot::PlotShaded("##value_filter", args.values.x, y_vals, args.values.count, -FLT_MAX);
        }

        if (args.filter.show) {
            if (!args.filter.enabled) ImGui::PushDisabled();
            ImPlot::DragRangeX("Time Filter", args.filter.beg, args.filter.end, args.filter.min, args.filter.max);
            if (!args.filter.enabled) ImGui::PopDisabled();
            *args.filter.beg = CLAMP(*args.filter.beg, args.filter.min, args.filter.max);
            *args.filter.end = CLAMP(*args.filter.end, args.filter.min, args.filter.max);
        }

        ImVec4 line_col = ImGui::ColorConvertU32ToFloat4(args.col);
        ImPlot::SetNextLineStyle(line_col);

        if (args.values.count > 0) {
            ASSERT(args.values.x);
            ASSERT(args.values.y);
            if (args.values.y_var) {
                ASSERT(args.values.y_min);
                ASSERT(args.values.y_max);
                char lbl[32];

                snprintf(lbl, sizeof(lbl), "%s", args.lbl);
                TimePayload payload = {
                    .args = &args,
                    .y_idx = 0,
                };
                for (int i = 0; i < args.values.dim_y; ++i) {
                    payload.y_idx = i;
                    ImPlot::PlotLineG(lbl, get_time_point, &payload, args.values.count);
                }

                ImPlot::SetNextLineStyle(line_col);
                snprintf(lbl, sizeof(lbl), "%s (mean)", args.lbl);
                ImPlot::PlotLine(lbl, args.values.x, args.values.y_mean, args.values.count);

                ImPlot::SetNextFillStyle(line_col, 0.4f);
                snprintf(lbl, sizeof(lbl), "%s (var)", args.lbl);
                ImPlot::PlotShadedG(lbl,
                    [](int idx, void* payload) -> ImPlotPoint {
                        TimelineArgs* args = (TimelineArgs*)payload;
                        return ImPlotPoint(args->values.x[idx], args->values.y_mean[idx] - args->values.y_var[idx]);
                    },
                    (void*)&args,
                    [](int idx, void* payload) -> ImPlotPoint {
                        TimelineArgs* args = (TimelineArgs*)payload;
                        return ImPlotPoint(args->values.x[idx], args->values.y_mean[idx] + args->values.y_var[idx]);
                    },
                    (void*)&args,
                    args.values.count
                );

                ImPlot::SetNextFillStyle(line_col, 0.2f);
                snprintf(lbl, sizeof(lbl), "%s (min,max)", args.lbl);
                ImPlot::PlotShadedG(lbl,
                    [](int idx, void* payload) -> ImPlotPoint {
                        TimelineArgs* args = (TimelineArgs*)payload;
                        return ImPlotPoint(args->values.x[idx], args->values.y_min[idx]);
                    },
                    (void*)&args,
                        [](int idx, void* payload) -> ImPlotPoint {
                        TimelineArgs* args = (TimelineArgs*)payload;
                        return ImPlotPoint(args->values.x[idx], args->values.y_max[idx]);
                    },
                    (void*)&args,
                    args.values.count
                );
            } 
            else {
                ImPlot::PlotLine(args.lbl, args.values.x, args.values.y, args.values.count);
            }
        }
        
        if (*args.input.is_dragging) {
            *args.time = ImPlot::GetPlotMousePos().x;
        }
        else if (*args.input.is_selecting) {
            *args.filter.end = MAX(ImPlot::GetPlotMousePos().x, *args.filter.beg);
        }
        else if (ImPlot::IsPlotHovered()) {
            if (active && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {           
                if (ImGui::GetIO().KeyMods == ImGuiMod_Shift) {
                    if (args.filter.show && args.filter.enabled) {
                        *args.input.is_selecting = true;
                        *args.filter.beg = ImPlot::GetPlotMousePos().x;
                    }
                } else {
                    *args.input.is_dragging = true;
                }
            }
        }

        if (ImPlot::IsPlotHovered()) {
            print_timeline_tooltip = true;
        }
        
        if (ImPlot::DragLineX(0, args.time, ImVec4(1,1,0,1))) {
            *args.time = CLAMP(*args.time, args.filter.min, args.filter.max);
        }
        if (ImGui::IsItemHovered()) {
            print_timeline_tooltip = true;
        }

        if (print_timeline_tooltip) {
            ImPlotPoint plot_pos = ImPlot::GetPlotMousePos();
            ImVec2 screen_pos = ImPlot::PlotToPixels(plot_pos);
            ImVec2 p0 = {screen_pos.x, ImPlot::GetPlotPos().y};
            ImVec2 p1 = {screen_pos.x, ImPlot::GetPlotPos().y + ImPlot::GetPlotSize().y};
            ImPlot::PushPlotClipRect();
            ImPlot::GetPlotDrawList()->AddLine(p0, p1, IM_COL32(255, 255, 255, 120));
            ImPlot::PopPlotClipRect();

            char buf[128] = "";
            int len = 0;

            double time = plot_pos.x;
            int frame_idx = CLAMP((int)(time_to_frame(time, data.timeline.x_values) + 0.5), 0, (int)md_array_size(data.timeline.x_values)-1);
            len += snprintf(buf + len, MAX(0, (int)sizeof(buf) - len), "time: %.2f", time);

            md_unit_t time_unit = data.timeline.time_unit;
            if (!md_unit_is_none(time_unit)) {
                char unit_buf[32];
                md_unit_print(unit_buf, sizeof(unit_buf), time_unit);
                len += snprintf(buf + len, MAX(0, (int)sizeof(buf) - len), " (%s)", unit_buf);
            }

            if (0 <= frame_idx && frame_idx < args.values.count) {
                const char* value_lbl = args.values.y_var ? "mean" : "value";
                if (args.values.y) {
                    len += snprintf(buf + len, MAX(0, (int)sizeof(buf) - len), ", %s: %.2f", value_lbl, args.values.y[frame_idx]);
                }
                if (args.values.y_var) {
                    ASSERT(args.values.y_min);
                    ASSERT(args.values.y_max);
                    len += snprintf(buf + len, MAX(0, (int)sizeof(buf) - len), ", var: %.2f, min: %.2f, max: %.2f",
                        args.values.y_var[frame_idx],
                        args.values.y_min[frame_idx],
                        args.values.y_max[frame_idx]);
                }

                if (!str_empty(args.values.unit)) {
                    len += snprintf(buf + len, MAX(0, (int)sizeof(buf) - len), " (%.*s)", (int)args.values.unit.len, args.values.unit.ptr);
                }
            }
            ImGui::SetTooltip("%.*s", len, buf);
        }

        ImPlot::EndPlot();
    }

    return true;
}

static double distance_to_linesegment(ImPlotPoint p0, ImPlotPoint p1, ImPlotPoint p) {
    double vx = p1.x - p0.x;
    double vy = p1.y - p0.y;

    double ux = p.x - p0.x;
    double uy = p.y - p0.y;

    double d_uv = ux*vx + uy*vy;
    double d_vv = vx*vx + vy*vy;

    if (d_vv < 1.0e-7) {
        double d_uu = ux*ux + uy*uy;
        return sqrt(d_uu);
    }

    double t = d_uv / d_vv;

    if (t < 0.0) {
        double d_uu = ux*ux + uy*uy;
        return sqrt(d_uu);
    } else if (t > 1.0) {
        double wx = p.x - p1.x;
        double wy = p.y - p1.y;
        double d_ww = wx*wx + wy*wy;
        return sqrt(d_ww);
    } else {
        double wx = p.x - (p0.x + vx * t);
        double wy = p.y - (p0.y + vy * t);
        double d_ww = wx*wx + wy*wy;
        return sqrt(d_ww);
    }
}

// #plots
//
// The Timelines and Distributions windows. Both draw series named by attribute path (see
// plot_series.h) into subplots that own them, and share the lists, the drag and drop and the
// styling of a legend entry.

// The colour of population member k of an entry drawn with a population of population_size
static ImVec4 plot_series_member_color(const PlotSeries& s, int k, int population_size) {
    ImVec4 color = s.color;
    if (s.use_colormap && population_size > 1) {
        if (ImPlot::ColormapQualitative(s.colormap)) {
            color = ImPlot::GetColormapColor(k, s.colormap);
        } else {
            color = ImPlot::SampleColormap((float)k / (float)(population_size - 1), s.colormap);
        }
        color.w *= s.colormap_alpha;
    }
    return color;
}

static ImVec4 plot_highlight(ImVec4 color) {
    const float scl = 1.5f;
    return ImVec4(ImSaturate(color.x * scl), ImSaturate(color.y * scl), ImSaturate(color.z * scl), color.w);
}

// The part of a legend entry's popup both windows share: colour or colormap, and which members of
// a population are drawn. A member hovered in the popup is written to hovered_pop_idx.
static void plot_series_style_popup(ApplicationState* data, PlotSeries& s, int dim, const char* script_ident, const md_script_vis_payload_o* vis_payload, int* hovered_pop_idx) {
    if (dim > 1) {
        const char* color_type_labels[] = {"Solid", "Colormap"};
        int color_type = s.use_colormap ? 1 : 0;
        if (ImGui::Combo("Color Type", &color_type, color_type_labels, IM_ARRAYSIZE(color_type_labels))) {
            s.use_colormap = color_type == 1;
        }
    }
    if (s.use_colormap && dim > 1) {
        ImPlot::ColormapSelection("##Colormap", &s.colormap);
        ImGui::SliderFloat("Alpha", &s.colormap_alpha, 0.0f, 1.0f);
    } else {
        ImGui::ColorEdit4("Color", &s.color.x);
    }

    if (dim > 1) {
        ImGui::Separator();
        if (ImGui::Button("Set All")) {
            s.population_mask.set();
        }
        ImGui::SameLine();
        if (ImGui::Button("Clear All")) {
            s.population_mask.reset();
        }

        const float sz = ImGui::GetFontSize() * 1.5f;
        ImGui::PushStyleVar(ImGuiStyleVar_SelectableTextAlign, ImVec2(0.5f, 0.5f));
        for (int k = 0; k < MIN(dim, PLOT_MAX_POPULATION); ++k) {
            char lbl[32];
            snprintf(lbl, sizeof(lbl), "%d", k+1);
            if (ImGui::Selectable(lbl, s.population_mask.test(k), ImGuiSelectableFlags_DontClosePopups, ImVec2(sz, sz))) {
                s.population_mask.flip(k);
            }
            if (ImGui::IsItemHovered()) {
                if (script_ident[0] != '\0') {
                    script_visualize_payload(data, vis_payload, k, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                    script_set_hovered_property(data, str_from_cstr(script_ident), k);
                }
                *hovered_pop_idx = k;
            }
            if (!k || ((k+1) % 10)) {
                ImGui::SameLine();
            }
        }
        ImGui::PopStyleVar();
    }
}

// One entry in a list of series: a click toggles it in the first subplot, a drag puts it in any
static void plot_series_list_item(ApplicationState* data, const char* dnd_type, PlotSubplot& first, const SeriesKey& key) {
    char label[96];
    series_label(label, sizeof(label), data, key);
    const ImVec4 color = series_default_color(data, key);
    const int idx = plot_find_series(first, key);

    ImGui::PushID(key.path);
    ImGui::PushID((int)key.variant * SeriesSource_Count + (int)key.source);
    ImPlot::ItemIcon(color);
    ImGui::SameLine();
    if (ImGui::Selectable(label, idx != -1, ImGuiSelectableFlags_DontClosePopups)) {
        if (idx != -1) {
            plot_remove_series(first, idx);
        } else {
            plot_add_series(data, first, key);
        }
    }
    if (ImGui::IsItemHovered()) {
        const str_t ident = series_script_ident(key);
        if (!str_empty(ident)) {
            const md_script_vis_payload_o* vis = data->script.eval_ir ? md_script_ir_property_vis_payload(data->script.eval_ir, ident) : nullptr;
            script_visualize_payload(data, vis, -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
            script_set_hovered_property(data, ident);
        } else {
            ImGui::SetTooltip("%s\nClick to show in the first subplot, drag into any", key.path);
        }
    }
    if (ImGui::BeginDragDropSource()) {
        series_set_drag_payload(dnd_type, key, -1, label, color);
        ImGui::EndDragDropSource();
    }
    ImGui::PopID();
    ImGui::PopID();
}

// The series loaded along the run, one submenu per file
static void plot_system_series_menu(ApplicationState* data, const char* dnd_type, PlotSubplot& first) {
    str_t groups[64];
    const size_t num_groups = MIN(system_series_groups(groups, ARRAY_SIZE(groups), data), ARRAY_SIZE(groups));
    size_t num_listed = 0;
    for (size_t g = 0; g < num_groups; ++g) {
        md_temp_scope_t temp = md_temp_begin();
        const size_t num = system_series_members(nullptr, 0, data, groups[g]);
        str_t* paths = md_temp_alloc_array(temp, str_t, num + 1);
        system_series_members(paths, num, data, groups[g]);
        if (num > 0) {
            char group_label[256];
            system_series_group_label(group_label, sizeof(group_label), data, groups[g]);
            if (ImGui::BeginMenu(group_label)) {
                for (size_t i = 0; i < num; ++i) {
                    plot_series_list_item(data, dnd_type, first, series_key(SeriesSource_System, paths[i]));
                }
                ImGui::EndMenu();
            }
            num_listed += 1;
        }
        md_temp_end(temp);
    }
    if (num_listed == 0) {
        ImGui::TextDisabled("Nothing loaded");
    }
}

// The members of a script property's population, if it has one: the extent of its value axis
static size_t script_property_population(const ApplicationState* data, const SeriesKey& key) {
    const md_attributes_t* table = series_table(data, key.source);
    const md_attribute_t* attr = table ? md_attributes_find(table, str_from_cstr(key.path)) : nullptr;
    if (!attr || attr->format.rank < 1 || attr->format.shape[0] == 0) return 0;
    return md_attribute_element_count(&attr->format) / attr->format.shape[0];
}

// #timeline

// What a getter is handed: the resolved series, and which member of its population
struct TemporalGetterPayload {
    const SeriesTemporalView* view;
    int k;
};

static ImPlotPoint temporal_getter_line(int i, void* payload) {
    const TemporalGetterPayload* p = (const TemporalGetterPayload*)payload;
    return series_temporal_point(*p->view, i, p->k);
}

static ImPlotPoint temporal_getter_band_lo(int i, void* payload) {
    const TemporalGetterPayload* p = (const TemporalGetterPayload*)payload;
    double lo, hi;
    series_temporal_band(*p->view, i, &lo, &hi);
    return ImPlotPoint(p->view->x[i], lo);
}

static ImPlotPoint temporal_getter_band_hi(int i, void* payload) {
    const TemporalGetterPayload* p = (const TemporalGetterPayload*)payload;
    double lo, hi;
    series_temporal_band(*p->view, i, &lo, &hi);
    return ImPlotPoint(p->view->x[i], hi);
}

static void draw_timeline_properties_menu(ApplicationState* data) {
    PlotSubplot& first = data->timeline.subplots[0];

    ImGui::SeparatorText("Script");
    int num_listed = 0;
    series_for_each_script_property(data, SeriesSource_Script, MD_SCRIPT_PROPERTY_FLAG_TEMPORAL, [&](const SeriesKey& key) {
        plot_series_list_item(data, TIMELINE_SERIES_DND, first, key);
        num_listed += 1;

        // A property with a population also has its summary over it
        char mean_path[SERIES_PATH_CAP + 8];
        const int len = snprintf(mean_path, sizeof(mean_path), "%s/mean", key.path);
        const md_attributes_t* table = series_table(data, key.source);
        if (table && len > 0 && md_attributes_find(table, str_t{mean_path, (size_t)len})) {
            ImGui::Indent();
            const SeriesVariant variants[] = { SeriesVariant_Mean, SeriesVariant_Sigma, SeriesVariant_Extent };
            for (SeriesVariant v : variants) {
                SeriesKey summary = key;
                summary.variant = v;
                plot_series_list_item(data, TIMELINE_SERIES_DND, first, summary);
            }
            ImGui::Unindent();
        }
    });
    if (num_listed == 0) {
        ImGui::TextDisabled("No temporal properties, define and evaluate them in the script editor");
    }

    ImGui::SeparatorText("Loaded along the run");
    plot_system_series_menu(data, TIMELINE_SERIES_DND, first);
}

static void draw_timeline_window(ApplicationState* data) {
    ASSERT(data);
    ImGui::SetNextWindowSize(ImVec2(600, 300), ImGuiCond_FirstUseEver);

    if (ImGui::Begin("Timelines", &data->timeline.show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_MenuBar)) {
        int& num_subplots = data->timeline.num_subplots;
        num_subplots = CLAMP(num_subplots, 1, PLOT_MAX_SUBPLOTS);

        double pre_filter_min = data->timeline.filter.beg_frame;
        double pre_filter_max = data->timeline.filter.end_frame;

        const float* x_values   = data->timeline.x_values;
        const int num_x_values  = (int)md_array_size(data->timeline.x_values);
        const float min_x_value = num_x_values > 0 ? x_values[0] : 0.0f;
        const float max_x_value = num_x_values > 0 ? x_values[num_x_values - 1] : 1.0f;

        ImPlot::PushStyleVar(ImPlotStyleVar_PlotPadding, ImVec2(ImPlot::GetStyle().PlotPadding.x, 2));
        defer { ImPlot::PopStyleVar(); };

        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("Properties")) {
                draw_timeline_properties_menu(data);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Filter")) {
                ImGui::Checkbox("Enabled", &data->timeline.filter.enabled);
                if (data->timeline.filter.enabled) {
                    ImGui::Checkbox("Temporal Window", &data->timeline.filter.temporal_window.enabled);
                    if (data->timeline.filter.temporal_window.enabled) {
                        const double extent_min = 1.0;
                        const double extent_max = num_x_values / 2.0;
                        ImGui::SliderScalar("Extent (frames)", ImGuiDataType_Double, &data->timeline.filter.temporal_window.extent_in_frames, &extent_min, &extent_max, "%1.0f");
                    }
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Subplots")) {
                ImGui::SliderInt("Num Subplots", &num_subplots, 1, PLOT_MAX_SUBPLOTS);
                if (ImGui::Button("Add Subplot")) {
                    num_subplots = CLAMP(num_subplots + 1, 1, PLOT_MAX_SUBPLOTS);
                }
                if (ImGui::Button("Clear All")) {
                    plot_clear(data->timeline.subplots, PLOT_MAX_SUBPLOTS);
                }
                ImGui::SeparatorText("Names");
                for (int s = 0; s < num_subplots; ++s) {
                    char hint[32];
                    snprintf(hint, sizeof(hint), "Subplot %d", s + 1);
                    ImGui::PushID(s);
                    ImGui::InputTextWithHint("##name", hint, data->timeline.subplots[s].name, sizeof(data->timeline.subplots[s].name));
                    ImGui::PopID();
                }
                ImGui::SetItemTooltip("Names the subplots, e.g. for the figures of the movie, which find a subplot by its identity and not by its position.");
                ImGui::EndMenu();
            }
            ImGui::EndMenuBar();
        }

        if (ImGui::IsWindowFocused() && !data->movie.show_window && ImGui::IsKeyPressed(KEY_PLAY_PAUSE, false)) {
            data->animation.mode = data->animation.mode == PlaybackMode::Playing ? PlaybackMode::Stopped : PlaybackMode::Playing;
        }

        if (num_x_values > 0) {
            ImPlotInputMap old_map = ImPlot::GetInputMap();

            static bool is_dragging  = false;
            static bool is_selecting = false;

            if (!ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
                is_dragging = false;
                is_selecting = false;
            }
            if (!ImGui::IsKeyDown(ImGuiKey_LeftCtrl)) {
                is_dragging = false;
            }
            if (!ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                is_selecting = false;
            }

            // Create a temporary 'time' representation of the filters min and max value
            // The visualization uses time units while we store 'frame' units
            double filter_beg = frame_to_time(data->timeline.filter.beg_frame, *data);
            double filter_end = frame_to_time(data->timeline.filter.end_frame, *data);
            double time = frame_to_time(data->animation.frame, *data);

            ImPlot::BeginSubplots("##Temporal", num_subplots, 1, ImVec2(-1,-1));

            const ImPlotFlags plot_flags = ImPlotFlags_NoBoxSelect | ImPlotFlags_NoFrame;
            const ImPlotAxisFlags axis_flags   = ImPlotAxisFlags_NoSideSwitch;
            const ImPlotAxisFlags axis_flags_y = axis_flags | ImPlotAxisFlags_Opposite | ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_RangeFit;

            char x_label[64] = "Frame";
            char x_unit_str[32] = "";
            md_unit_t x_unit = data->timeline.time_unit;
            if (!md_unit_is_none(x_unit)) {
                md_unit_print(x_unit_str, sizeof(x_unit_str), x_unit);
                snprintf(x_label, sizeof(x_label), "Time (%s)", x_unit_str);
            }

            for (int i = 0; i < num_subplots; ++i) {
                PlotSubplot& sp = data->timeline.subplots[i];

                // Every entry is resolved once per plot, and everything below reads the views
                SeriesTemporalView views[PLOT_MAX_SERIES_PER_SUBPLOT];
                bool resolved[PLOT_MAX_SERIES_PER_SUBPLOT] = {};
                for (int j = 0; j < sp.count; ++j) {
                    resolved[j] = series_resolve_temporal(&views[j], data, sp.series[j].key);
                }

                int remove_idx = -1;

                if (ImPlot::BeginPlot("", ImVec2(), plot_flags)) {
                    ImPlot::SetupAxisLinks(ImAxis_X1, &data->timeline.view_range.beg_x, &data->timeline.view_range.end_x);
                    ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, min_x_value, max_x_value);

                    ImPlotAxisFlags axis_flags_x = axis_flags | ImPlotAxisFlags_NoLabel;
                    if (i < num_subplots - 1) {
                        // Only show label and ticklabels for the last plot, since they are all synced on x-axis
                        axis_flags_x |= ImPlotAxisFlags_NoTickLabels;
                    }

                    // When every series in the subplot is shown in the same unit, that unit labels the y axis
                    char y_unit_str[32] = "";
                    for (int j = 0; j < sp.count; ++j) {
                        if (!resolved[j]) continue;
                        if (y_unit_str[0] == '\0') {
                            str_copy_to_char_buf(y_unit_str, sizeof(y_unit_str), str_from_cstr(views[j].unit_str));
                        } else if (strcmp(y_unit_str, views[j].unit_str) != 0) {
                            y_unit_str[0] = '\0';
                            break;
                        }
                    }

                    char y_label[64] = "";
                    if (y_unit_str[0] != '\0') {
                        snprintf(y_label, sizeof(y_label), "(%s)", y_unit_str);
                    }

                    ImPlot::SetupAxes(x_label, y_label, axis_flags_x, axis_flags_y);
                    ImPlot::SetupFinish();

                    if (data->timeline.filter.enabled) {
                        bool disabled = data->timeline.filter.temporal_window.enabled;
                        ImPlotDragRangeFlags flags = ImPlotDragToolFlags_NoFit;
                        if (i < num_subplots - 1) {
                            flags |= ImPlotDragRangeFlags_NoBar;
                        }
                        if (disabled) ImGui::PushDisabled();
                        ImPlot::DragRangeX("Time Filter", &filter_beg, &filter_end, min_x_value, max_x_value, flags);
                        if (disabled) ImGui::PopDisabled();
                    }

                    // The entry, and the member of its population, under the mouse
                    int  hovered_idx      = -1;
                    int  hovered_pop_idx  = -1;
                    char hovered_label[128] = "";

                    bool print_timeline_tooltip = false;

                    if (ImPlot::IsPlotHovered()) {
                        md_bitfield_clear(&data->selection.highlight_mask);
                        script_set_hovered_property(data,  STR_LIT(""));

                        print_timeline_tooltip = true;
                        const ImPlotPoint mouse_pos = ImPlot::GetPlotMousePos();
                        const ImVec2 mouse_coord = ImPlot::PlotToPixels(mouse_pos);
                        const float max_rad = 20; // 20 pixels
                        const float area_dist = max_rad * 0.2f;

                        float min_dist = max_rad;

                        for (int j = 0; j < sp.count; ++j) {
                            if (!resolved[j]) continue;
                            const PlotSeries& s = sp.series[j];
                            const SeriesTemporalView& v = views[j];

                            ImPlotItem* item = ImPlot::GetItem(v.plot_id);
                            if (!item || !item->Show) {
                                continue;
                            }

                            // Each series has an axis of its own, so the samples around the mouse are its own
                            const int n = v.num_samples;
                            const double fi = series_temporal_index_at(v, mouse_pos.x);
                            const int fn = CLAMP((int)(fi + 0.5), 0, n - 1);  // Nearest index
                            const int f[4] = {
                                CLAMP((int)fi - 1, 0, n - 1),
                                CLAMP((int)fi,     0, n - 1),
                                CLAMP((int)fi + 1, 0, n - 1),
                                CLAMP((int)fi + 2, 0, n - 1),
                            };

                            const int dim = CLAMP(v.dim, 1, PLOT_MAX_POPULATION);
                            for (int k = 0; k < dim; ++k) {
                                if (dim > 1 && !s.population_mask.test(k)) {
                                    continue;
                                }
                                TemporalGetterPayload payload = { &v, k };
                                double d = DBL_MAX;

                                if (v.band) {
                                    const ImVec2 p_min[4] = {
                                        ImPlot::PlotToPixels(temporal_getter_band_lo(f[0], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_band_lo(f[1], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_band_lo(f[2], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_band_lo(f[3], &payload)),
                                    };
                                    const ImVec2 p_max[4] = {
                                        ImPlot::PlotToPixels(temporal_getter_band_hi(f[0], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_band_hi(f[1], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_band_hi(f[2], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_band_hi(f[3], &payload)),
                                    };
                                    // Each segment forms a trapezoid with its left and right sides parallel to the y axis:
                                    // clamp the mouse to it and measure the distance to the clamped point
                                    for (int l = 0; l < 2; ++l) {
                                        const float x_min = MIN(p_min[l].x, p_min[l+1].x);
                                        const float x_max = MAX(p_min[l].x, p_min[l+1].x);
                                        const float t = (x_max > x_min) ? CLAMP((mouse_coord.x - x_min) / (x_max - x_min), 0.0f, 1.0f) : 0.0f;
                                        const float y[2] = {
                                            lerp(p_min[l].y, p_min[l+1].y, t),
                                            lerp(p_max[l].y, p_max[l+1].y, t)
                                        };
                                        const ImVec2 p = {
                                            CLAMP(mouse_coord.x, x_min, x_max),
                                            CLAMP(mouse_coord.y, MIN(y[0], y[1]), MAX(y[0], y[1]))
                                        };
                                        d = MIN(d, sqrt(ImLengthSqr(mouse_coord - p)));
                                    }
                                    d += area_dist;
                                } else if (s.plot_type == PlotType_Scatter) {
                                    ImVec2 p = ImPlot::PlotToPixels(temporal_getter_line(fn, &payload));
                                    d = sqrt(ImLengthSqr(mouse_coord - p));
                                } else {
                                    // Distance to the segments before, at and after the mouse:
                                    // the current one alone is not enough near a vertex
                                    ImVec2 p[4] = {
                                        ImPlot::PlotToPixels(temporal_getter_line(f[0], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_line(f[1], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_line(f[2], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_line(f[3], &payload)),
                                    };
                                    d = distance_to_linesegment(p[0], p[1], mouse_coord);
                                    d = MIN(distance_to_linesegment(p[1], p[2], mouse_coord), d);
                                    d = MIN(distance_to_linesegment(p[2], p[3], mouse_coord), d);
                                }

                                if (d < min_dist) {
                                    min_dist = (float)d;
                                    char value_buf[64] = "";
                                    series_temporal_print_value(value_buf, sizeof(value_buf), v, fn, k);

                                    hovered_idx = j;
                                    hovered_pop_idx = k;
                                    if (v.dim > 1) {
                                        snprintf(hovered_label, sizeof(hovered_label), "%s[%i]: %s %s", v.label, k + 1, value_buf, v.unit_str);
                                    } else {
                                        snprintf(hovered_label, sizeof(hovered_label), "%s: %s %s", v.label, value_buf, v.unit_str);
                                    }
                                }
                            }
                        }

                        if (hovered_idx != -1 && views[hovered_idx].script_ident[0] != '\0') {
                            script_set_hovered_property(data, str_from_cstr(views[hovered_idx].script_ident), hovered_pop_idx);
                        }
                    } else if (data->hovered_property_label[0] != '\0') {
                        // Hovered in another view of it (the other plot window, a list): light up its values here
                        for (int j = 0; j < sp.count; ++j) {
                            if (!resolved[j] || sp.series[j].key.variant != SeriesVariant_Values) continue;
                            if (strcmp(data->hovered_property_label, views[j].script_ident) == 0) {
                                hovered_idx = j;
                                hovered_pop_idx = data->hovered_property_pop_idx;
                                break;
                            }
                        }
                    }

                    if (is_dragging) {
                        time = ImPlot::GetPlotMousePos().x;
                    } else if (is_selecting) {
                        filter_end = MAX(ImPlot::GetPlotMousePos().x, filter_beg);
                        filter_beg = MIN(filter_beg, filter_end);
                    } else if (ImPlot::IsPlotHovered() && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
                        if (ImGui::IsKeyDown(ImGuiKey_LeftCtrl)) {
                            is_dragging = true;
                        } else if (ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                            filter_beg = ImPlot::GetPlotMousePos().x;
                            is_selecting = true;
                        }
                    }

                    for (int j = 0; j < sp.count; ++j) {
                        PlotSeries& s = sp.series[j];
                        const SeriesTemporalView& v = views[j];

                        if (ImPlot::IsLegendEntryHovered(v.plot_id)) {
                            if (v.script_ident[0] != '\0') {
                                script_visualize_payload(data, v.vis_payload, -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                                script_set_hovered_property(data, str_from_cstr(v.script_ident));
                            }
                            if (!resolved[j]) {
                                ImGui::SetTooltip("Nothing to plot at '%s' right now", s.key.path);
                            }
                            hovered_idx = j;
                            hovered_pop_idx = -1;
                        }

                        // legend context menu
                        if (ImPlot::BeginLegendPopup(v.plot_id)) {
                            if (ImGui::DeleteButton("Remove")) {
                                remove_idx = j;
                                ImGui::CloseCurrentPopup();
                            }

                            if (resolved[j] && !v.band) {
                                const char* plot_type_names[] = {"Line", "Scatter"};
                                int plot_type = s.plot_type == PlotType_Scatter ? 1 : 0;
                                if (ImGui::Combo("Plot Type", &plot_type, plot_type_names, IM_ARRAYSIZE(plot_type_names))) {
                                    s.plot_type = plot_type == 1 ? PlotType_Scatter : PlotType_Line;
                                }
                                if (s.plot_type == PlotType_Scatter) {
                                    if (ImGui::BeginCombo("Marker", ImPlot::GetMarkerName(s.marker))) {
                                        for (int k = 0; k < ImPlotMarker_COUNT; ++k) {
                                            if (ImGui::Selectable(ImPlot::GetMarkerName(k), s.marker == k)) {
                                                s.marker = (ImPlotMarker)k;
                                            }
                                        }
                                        ImGui::EndCombo();
                                    }
                                    ImGui::SliderFloat("Marker Size", &s.marker_size, 0.1f, 10.0f, "%.2f");
                                }
                            }

                            int popup_hovered_pop = -1;
                            plot_series_style_popup(data, s, resolved[j] ? v.dim : 1, v.script_ident, v.vis_payload, &popup_hovered_pop);
                            if (popup_hovered_pop != -1) {
                                hovered_idx = j;
                                hovered_pop_idx = popup_hovered_pop;
                            }
                            ImPlot::EndLegendPopup();
                        }

                        if (!resolved[j]) {
                            // Still in the legend, so it can be seen and removed
                            ImPlot::PlotDummy(v.plot_id);
                        } else {
                            const int population_size = CLAMP(v.dim, 1, PLOT_MAX_POPULATION);
                            auto plot = [&](int k) {
                                ImVec4 color = plot_series_member_color(s, k, population_size);
                                ImVec4 marker_line_color = {};
                                float  marker_line_weight = 0;
                                float  fill_alpha = 1.0f;
                                float  weight = 1.0f;

                                if (hovered_idx == j) {
                                    if (hovered_pop_idx == -1 || hovered_pop_idx == k) {
                                        color = plot_highlight(color);
                                        fill_alpha = 1.25f;
                                        marker_line_color = {1,1,1,1};
                                        marker_line_weight = 1.0f;
                                    }
                                    if (hovered_pop_idx == k) {
                                        weight = 2.0f;
                                    }
                                }

                                TemporalGetterPayload payload = { &v, k };
                                if (v.band) {
                                    ImPlot::SetNextFillStyle(color, fill_alpha);
                                    ImPlot::PlotShadedG(v.plot_id, temporal_getter_band_lo, &payload, temporal_getter_band_hi, &payload, v.num_samples);
                                } else if (s.plot_type == PlotType_Scatter) {
                                    ImPlot::SetNextMarkerStyle(s.marker, s.marker_size, color, marker_line_weight, marker_line_color);
                                    ImPlot::PlotScatterG(v.plot_id, temporal_getter_line, &payload, v.num_samples);
                                } else {
                                    ImPlot::SetNextLineStyle(color, weight);
                                    ImPlot::PlotLineG(v.plot_id, temporal_getter_line, &payload, v.num_samples);
                                }
                            };

                            for (int k = 0; k < population_size; ++k) {
                                if (population_size > 1 && !s.population_mask.test(k)) continue;
                                if (hovered_idx == j && hovered_pop_idx == k) continue;
                                plot(k);
                            }
                            // The hovered member last, on top
                            if (hovered_idx == j && hovered_pop_idx != -1 && hovered_pop_idx < population_size) {
                                plot(hovered_pop_idx);
                            }
                        }

                        if (ImPlot::BeginDragDropSourceItem(v.plot_id)) {
                            series_set_drag_payload(TIMELINE_SERIES_DND, s.key, i, v.label, s.color);
                            ImPlot::EndDragDropSource();
                        }
                    }

                    if (ImPlot::DragLineX(0, &time, ImVec4(1,1,0,1), 1.0f, ImPlotDragToolFlags_NoFit)) {
                        time = CLAMP(time, min_x_value, max_x_value);
                    }

                    if (movie_draw_timeline_markers(data)) {
                        print_timeline_tooltip = false;
                    }

                    if (ImPlot::IsPlotHovered() && hovered_idx != -1 && views[hovered_idx].script_ident[0] != '\0') {
                        const SeriesTemporalView& v = views[hovered_idx];
                        const int pop_idx = v.dim > 1 ? hovered_pop_idx : -1;
                        script_visualize_payload(data, v.vis_payload, pop_idx, ~MD_SCRIPT_VISUALIZE_SDF);
                        script_set_hovered_property(data, str_from_cstr(v.script_ident), hovered_pop_idx);
                    }

                    if (ImPlot::BeginDragDropTargetPlot()) {
                        if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(TIMELINE_SERIES_DND)) {
                            ASSERT(payload->DataSize == sizeof(SeriesDragPayload));
                            const SeriesDragPayload* dnd = (const SeriesDragPayload*)(payload->Data);
                            if (dnd->src_subplot >= 0 && dnd->src_subplot < PLOT_MAX_SUBPLOTS) {
                                // Out of one plot and into another moves it, keeping its style
                                plot_move_series(data, data->timeline.subplots[dnd->src_subplot], sp, dnd->key);
                            } else {
                                plot_add_series(data, sp, dnd->key);
                            }
                        }
                        ImPlot::EndDragDropTarget();
                    }

                    if (print_timeline_tooltip) {
                        ImPlotPoint plot_pos = ImPlot::GetPlotMousePos();
                        ImVec2 screen_pos = ImPlot::PlotToPixels(plot_pos);
                        ImVec2 p0 = {screen_pos.x, ImPlot::GetPlotPos().y};
                        ImVec2 p1 = {screen_pos.x, ImPlot::GetPlotPos().y + ImPlot::GetPlotSize().y};
                        ImPlot::PushPlotClipRect();
                        ImPlot::GetPlotDrawList()->AddLine(p0, p1, IM_COL32(255, 255, 255, 120));
                        ImPlot::PopPlotClipRect();

                        double t = plot_pos.x;
                        if (md_unit_is_none(x_unit)) {
                            int32_t frame_idx = CLAMP((int)(time_to_frame(t, x_values) + 0.5), 0, num_x_values-1);
                            ImGui::SetTooltip("Frame: %i\n%s", frame_idx, hovered_label);
                        } else {
                            ImGui::SetTooltip("Time: %.2f (%s)\n%s", t, x_unit_str, hovered_label);
                        }
                    }

                    ImPlot::EndPlot();
                }

                if (remove_idx != -1) {
                    plot_remove_series(sp, remove_idx);
                }
            }

            ImPlot::EndSubplots();

            time       = CLAMP(time, (double)min_x_value, (double)max_x_value);
            filter_beg = CLAMP(filter_beg, min_x_value, max_x_value);
            filter_end = CLAMP(filter_end, min_x_value, max_x_value);

            data->animation.frame = time_to_frame(time, data->timeline.x_values);
            data->timeline.filter.beg_frame = time_to_frame(filter_beg, data->timeline.x_values);
            data->timeline.filter.end_frame = time_to_frame(filter_end, data->timeline.x_values);

            ImPlot::GetInputMap() = old_map;
        }

        if (data->timeline.filter.enabled && (data->timeline.filter.beg_frame != pre_filter_min || data->timeline.filter.end_frame != pre_filter_max)) {
            data->script.evaluate_filt = true;
        }

        // A series dragged out of a plot and dropped outside any of them is taken out of that plot
        if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
            const ImGuiPayload* payload = ImGui::GetDragDropPayload();
            if (payload && payload->IsDataType(TIMELINE_SERIES_DND) && !ImGui::IsDragDropPayloadBeingAccepted()) {
                const SeriesDragPayload* dnd = (const SeriesDragPayload*)(payload->Data);
                if (dnd && dnd->src_subplot >= 0 && dnd->src_subplot < PLOT_MAX_SUBPLOTS) {
                    PlotSubplot& src = data->timeline.subplots[dnd->src_subplot];
                    plot_remove_series(src, plot_find_series(src, dnd->key));
                }
            }
        }
    }
    ImGui::End();
}

// #distribution_window

struct HistogramGetterPayload {
    const SeriesHistogramView* view;
    int k;
};

static ImPlotPoint histogram_getter(int i, void* payload) {
    const HistogramGetterPayload* p = (const HistogramGetterPayload*)payload;
    return series_histogram_point(*p->view, i, p->k);
}

static ImPlotPoint histogram_getter_zero(int i, void* payload) {
    const HistogramGetterPayload* p = (const HistogramGetterPayload*)payload;
    return ImPlotPoint(series_histogram_point(*p->view, i, p->k).x, 0.0);
}

static void draw_distribution_properties_menu(ApplicationState* data) {
    PlotSubplot& first = data->distributions.subplots[0];
    const uint32_t kinds = MD_SCRIPT_PROPERTY_FLAG_TEMPORAL | MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION;

    // A temporal property is binned over its frames, a distribution is shown as evaluated. Over the
    // timeline filter's frames too, while there is a filter.
    auto script_section = [&](SeriesSource source) -> int {
        int num_listed = 0;
        series_for_each_script_property(data, source, kinds, [&](const SeriesKey& key) {
            plot_series_list_item(data, DISTRIBUTION_SERIES_DND, first, key);
            num_listed += 1;
            const bool temporal = md_script_ir_property_flags(data->script.eval_ir, series_script_ident(key)) & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL;
            if (temporal && script_property_population(data, key) > 1) {
                SeriesKey agg = key;
                agg.variant = SeriesVariant_Aggregate;
                ImGui::Indent();
                plot_series_list_item(data, DISTRIBUTION_SERIES_DND, first, agg);
                ImGui::Unindent();
            }
        });
        return num_listed;
    };

    ImGui::SeparatorText("Script");
    if (script_section(SeriesSource_Script) == 0) {
        ImGui::TextDisabled("No temporal or distribution properties, define and evaluate them in the script editor");
    }
    if (data->timeline.filter.enabled) {
        ImGui::SeparatorText("Script, over the timeline filter");
        script_section(SeriesSource_ScriptFiltered);
    }

    ImGui::SeparatorText("Loaded along the run");
    plot_system_series_menu(data, DISTRIBUTION_SERIES_DND, first);
}

static void draw_distribution_window(ApplicationState* data) {
    ImGui::SetNextWindowSize(ImVec2(200, 300), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Distributions", &data->distributions.show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_MenuBar)) {
        int& num_subplots = data->distributions.num_subplots;
        num_subplots = CLAMP(num_subplots, 1, PLOT_MAX_SUBPLOTS);

        // Refit a subplot's axes to what it holds, after it gained a series or its binning changed
        static bool fit_subplot[PLOT_MAX_SUBPLOTS] = {};

        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("Properties")) {
                draw_distribution_properties_menu(data);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Subplots")) {
                ImGui::SliderInt("Num Subplots", &num_subplots, 1, PLOT_MAX_SUBPLOTS);
                if (ImGui::Button("Add Subplot")) {
                    num_subplots = CLAMP(num_subplots + 1, 1, PLOT_MAX_SUBPLOTS);
                }
                if (ImGui::Button("Clear All")) {
                    plot_clear(data->distributions.subplots, PLOT_MAX_SUBPLOTS);
                }
                ImGui::SeparatorText("Names");
                for (int s = 0; s < num_subplots; ++s) {
                    char hint[32];
                    snprintf(hint, sizeof(hint), "Subplot %d", s + 1);
                    ImGui::PushID(s);
                    ImGui::InputTextWithHint("##name", hint, data->distributions.subplots[s].name, sizeof(data->distributions.subplots[s].name));
                    ImGui::PopID();
                }
                ImGui::SetItemTooltip("Names the subplots, e.g. for the figures of the movie, which find a subplot by its identity and not by its position.");
                ImGui::EndMenu();
            }
            ImGui::EndMenuBar();
        }

        const ImPlotAxisFlags axis_flags   = ImPlotAxisFlags_NoSideSwitch | ImPlotAxisFlags_NoHighlight;
        const ImPlotAxisFlags axis_flags_x = axis_flags;
        const ImPlotAxisFlags axis_flags_y = axis_flags | ImPlotAxisFlags_AutoFit;
        const ImPlotFlags     plot_flags   = ImPlotFlags_NoBoxSelect | ImPlotFlags_NoFrame;

        if (ImPlot::BeginSubplots("##distribution_plots", num_subplots, 1, ImVec2(-1,-1))) {
            for (int i = 0; i < num_subplots; ++i) {
                PlotSubplot& sp = data->distributions.subplots[i];

                SeriesHistogramView views[PLOT_MAX_SERIES_PER_SUBPLOT];
                bool resolved[PLOT_MAX_SERIES_PER_SUBPLOT] = {};
                for (int j = 0; j < sp.count; ++j) {
                    resolved[j] = series_resolve_histogram(&views[j], data, sp.series[j].key, sp.series[j].num_bins);
                }

                int remove_idx = -1;

                if (fit_subplot[i]) {
                    ImPlot::SetNextAxesToFit();
                    fit_subplot[i] = false;
                }

                if (ImPlot::BeginPlot("", ImVec2(-1,0), plot_flags)) {
                    // Labelled with what the series are shown in; empty when there is no unit, or several
                    char x_label[64] = "";
                    char y_label[64] = "";
                    bool first_label = true;
                    for (int j = 0; j < sp.count; ++j) {
                        if (!resolved[j]) continue;
                        if (first_label) {
                            str_copy_to_char_buf(x_label, sizeof(x_label), str_from_cstr(views[j].x_unit_str));
                            str_copy_to_char_buf(y_label, sizeof(y_label), str_from_cstr(views[j].y_unit_str));
                            first_label = false;
                        } else {
                            if (strcmp(x_label, views[j].x_unit_str) != 0) x_label[0] = '\0';
                            if (strcmp(y_label, views[j].y_unit_str) != 0) y_label[0] = '\0';
                        }
                    }

                    ImPlot::SetupAxes(x_label, y_label, axis_flags_x, axis_flags_y);
                    ImPlot::SetupFinish();

                    int  hovered_idx       = -1;
                    int  hovered_pop_idx   = -1;
                    char hovered_label[128] = "";

                    if (ImPlot::IsPlotHovered()) {
                        script_set_hovered_property(data, STR_LIT(""));

                        const ImPlotPoint mouse_pos = ImPlot::GetPlotMousePos();
                        const ImVec2 mouse_coord = ImPlot::PlotToPixels(mouse_pos);

                        const double max_rad = 20; // 20 pixels
                        const double area_dist = max_rad * 0.2;
                        double min_dist = max_rad;

                        for (int j = 0; j < sp.count; ++j) {
                            if (!resolved[j]) continue;
                            const PlotSeries& s = sp.series[j];
                            const SeriesHistogramView& v = views[j];
                            if (v.num_bins <= 0 || v.x_max <= v.x_min) continue;

                            ImPlotItem* item = ImPlot::GetItem(v.plot_id);
                            if (!item || !item->Show) continue;

                            const int nb = v.num_bins;
                            const double x = series_histogram_index_at(v, mouse_pos.x);
                            const int xn = CLAMP((int)(x + 0.5), 0, nb - 1);
                            const int xi[4] = {
                                CLAMP((int)x - 1, 0, nb - 1),
                                CLAMP((int)x,     0, nb - 1),
                                CLAMP((int)x + 1, 0, nb - 1),
                                CLAMP((int)x + 2, 0, nb - 1),
                            };

                            const int dim = CLAMP(v.dim, 1, PLOT_MAX_POPULATION);
                            for (int k = 0; k < dim; ++k) {
                                if (dim > 1 && !s.population_mask.test(k)) continue;
                                HistogramGetterPayload payload = { &v, k };
                                // Later entries lie on top: prefer them a little when two are equally near
                                const double layer = j + k / (double)(PLOT_MAX_POPULATION - 1);
                                const double layer_dist = ((sp.count - layer) / sp.count) * (max_rad * 0.1);
                                double d = DBL_MAX;

                                if (s.plot_type == PlotType_Bars) {
                                    const double half_width = 0.5 * series_histogram_bin_width(v) * s.bar_width;
                                    for (int l = 0; l < 3; ++l) {
                                        const ImPlotPoint top = histogram_getter(xi[l], &payload);
                                        const ImVec2 p0 = ImPlot::PlotToPixels(ImPlotPoint(top.x - half_width, 0.0));
                                        const ImVec2 p1 = ImPlot::PlotToPixels(ImPlotPoint(top.x + half_width, top.y));
                                        const ImVec2 ll = ImMin(p0, p1);
                                        const ImVec2 ur = ImMax(p0, p1);
                                        d = MIN(d, sqrt(ImLengthSqr(mouse_coord - ImClamp(mouse_coord, ll, ur))));
                                    }
                                    d += area_dist + layer_dist;
                                } else if (s.plot_type == PlotType_Area) {
                                    for (int l = 0; l < 2; ++l) {
                                        const ImVec2 a = ImPlot::PlotToPixels(histogram_getter(xi[l], &payload));
                                        const ImVec2 b = ImPlot::PlotToPixels(histogram_getter(xi[l + 1], &payload));
                                        const ImVec2 z = ImPlot::PlotToPixels(histogram_getter_zero(xi[l], &payload));
                                        const float x_min = MIN(a.x, b.x);
                                        const float x_max = MAX(a.x, b.x);
                                        const float t = (x_max > x_min) ? CLAMP((mouse_coord.x - x_min) / (x_max - x_min), 0.0f, 1.0f) : 0.0f;
                                        const float y_top = lerp(a.y, b.y, t);
                                        const ImVec2 p = {
                                            CLAMP(mouse_coord.x, x_min, x_max),
                                            CLAMP(mouse_coord.y, MIN(y_top, z.y), MAX(y_top, z.y))
                                        };
                                        d = MIN(d, sqrt(ImLengthSqr(mouse_coord - p)));
                                    }
                                    d += area_dist + layer_dist;
                                } else {
                                    ImVec2 p[4] = {
                                        ImPlot::PlotToPixels(histogram_getter(xi[0], &payload)),
                                        ImPlot::PlotToPixels(histogram_getter(xi[1], &payload)),
                                        ImPlot::PlotToPixels(histogram_getter(xi[2], &payload)),
                                        ImPlot::PlotToPixels(histogram_getter(xi[3], &payload)),
                                    };
                                    d = distance_to_linesegment(p[0], p[1], mouse_coord);
                                    d = MIN(distance_to_linesegment(p[1], p[2], mouse_coord), d);
                                    d = MIN(distance_to_linesegment(p[2], p[3], mouse_coord), d);
                                    d += layer_dist;
                                }

                                if (d < min_dist) {
                                    min_dist = d;
                                    const ImPlotPoint pt = histogram_getter(xn, &payload);
                                    hovered_idx = j;
                                    hovered_pop_idx = dim > 1 ? k : -1;
                                    if (dim > 1) {
                                        snprintf(hovered_label, sizeof(hovered_label), "%s[%i]: %.4g at %.4g %s", v.label, k + 1, pt.y, pt.x, v.x_unit_str);
                                    } else {
                                        snprintf(hovered_label, sizeof(hovered_label), "%s: %.4g at %.4g %s", v.label, pt.y, pt.x, v.x_unit_str);
                                    }
                                }
                            }
                        }

                        if (hovered_idx != -1) {
                            const SeriesHistogramView& v = views[hovered_idx];
                            if (v.script_ident[0] != '\0') {
                                script_set_hovered_property(data, str_from_cstr(v.script_ident), hovered_pop_idx);
                                script_visualize_payload(data, v.vis_payload, hovered_pop_idx, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                            }
                            if (hovered_label[0] != '\0') {
                                ImGui::SetTooltip("%s", hovered_label);
                            }
                        }
                    } else if (data->hovered_property_label[0] != '\0') {
                        for (int j = 0; j < sp.count; ++j) {
                            if (!resolved[j]) continue;
                            if (strcmp(data->hovered_property_label, views[j].script_ident) == 0) {
                                hovered_idx = j;
                                hovered_pop_idx = data->hovered_property_pop_idx;
                                break;
                            }
                        }
                    }

                    for (int j = 0; j < sp.count; ++j) {
                        PlotSeries& s = sp.series[j];
                        const SeriesHistogramView& v = views[j];

                        if (ImPlot::IsLegendEntryHovered(v.plot_id)) {
                            if (v.script_ident[0] != '\0') {
                                script_visualize_payload(data, v.vis_payload, -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                                script_set_hovered_property(data, str_from_cstr(v.script_ident));
                            }
                            if (!resolved[j]) {
                                ImGui::SetTooltip("Nothing to show for '%s' right now", s.key.path);
                            }
                            hovered_idx = j;
                            hovered_pop_idx = -1;
                        }

                        if (ImPlot::BeginLegendPopup(v.plot_id)) {
                            if (ImGui::DeleteButton("Remove")) {
                                remove_idx = j;
                                ImGui::CloseCurrentPopup();
                            }

                            const char* plot_type_names[] = {"Line", "Area", "Bars"};
                            const PlotType plot_types[]   = {PlotType_Line, PlotType_Area, PlotType_Bars};
                            int plot_type = 0;
                            for (int k = 0; k < IM_ARRAYSIZE(plot_types); ++k) {
                                if (s.plot_type == plot_types[k]) plot_type = k;
                            }
                            if (ImGui::Combo("Plot Type", &plot_type, plot_type_names, IM_ARRAYSIZE(plot_type_names))) {
                                s.plot_type = plot_types[plot_type];
                            }
                            if (s.plot_type == PlotType_Bars) {
                                ImGui::SliderFloat("Bar Width", &s.bar_width, 0.01f, 1.0f, "%.3f");
                            }

                            // Powers of two: a script distribution only divides into those
                            int bins_exp = (int)log2((double)MAX(s.num_bins, 2));
                            char bins_fmt[32];
                            snprintf(bins_fmt, sizeof(bins_fmt), "%d", 1 << bins_exp);
                            if (ImGui::SliderInt("Num Bins", &bins_exp, 5, 10, bins_fmt)) {
                                s.num_bins = 1 << bins_exp;
                                fit_subplot[i] = true;
                            }
                            if (resolved[j] && v.num_bins != s.num_bins) {
                                ImGui::TextDisabled("Evaluated at %d bins", v.num_bins);
                            }

                            int popup_hovered_pop = -1;
                            plot_series_style_popup(data, s, resolved[j] ? v.dim : 1, v.script_ident, v.vis_payload, &popup_hovered_pop);
                            if (popup_hovered_pop != -1) {
                                hovered_idx = j;
                                hovered_pop_idx = popup_hovered_pop;
                            }
                            ImPlot::EndLegendPopup();
                        }

                        if (!resolved[j] || v.num_bins <= 0) {
                            ImPlot::PlotDummy(v.plot_id);
                        } else {
                            const int population_size = CLAMP(v.dim, 1, PLOT_MAX_POPULATION);
                            auto plot = [&](int k) {
                                ImVec4 color = plot_series_member_color(s, k, population_size);
                                float  fill_alpha = 1.0f;
                                float  weight = 1.0f;
                                if (hovered_idx == j) {
                                    if (hovered_pop_idx == -1 || hovered_pop_idx == k) {
                                        color = plot_highlight(color);
                                        fill_alpha = 1.25f;
                                    }
                                    if (hovered_pop_idx == k) {
                                        weight = 2.0f;
                                    }
                                }

                                HistogramGetterPayload payload = { &v, k };
                                switch (s.plot_type) {
                                case PlotType_Area:
                                    ImPlot::SetNextFillStyle(color, fill_alpha);
                                    ImPlot::PlotShadedG(v.plot_id, histogram_getter_zero, &payload, histogram_getter, &payload, v.num_bins);
                                    break;
                                case PlotType_Bars:
                                    ImPlot::SetNextFillStyle(color, fill_alpha);
                                    ImPlot::PlotBarsG(v.plot_id, histogram_getter, &payload, v.num_bins, series_histogram_bin_width(v) * s.bar_width);
                                    break;
                                default:
                                    ImPlot::SetNextLineStyle(color, weight);
                                    ImPlot::PlotLineG(v.plot_id, histogram_getter, &payload, v.num_bins);
                                    break;
                                }
                            };

                            for (int k = 0; k < population_size; ++k) {
                                if (population_size > 1 && !s.population_mask.test(k)) continue;
                                if (hovered_idx == j && hovered_pop_idx == k) continue;
                                plot(k);
                            }
                            if (hovered_idx == j && hovered_pop_idx != -1 && hovered_pop_idx < population_size) {
                                plot(hovered_pop_idx);
                            }
                        }

                        if (ImPlot::BeginDragDropSourceItem(v.plot_id)) {
                            series_set_drag_payload(DISTRIBUTION_SERIES_DND, s.key, i, v.label, s.color);
                            ImPlot::EndDragDropSource();
                        }
                    }

                    if (ImPlot::BeginDragDropTargetPlot()) {
                        if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(DISTRIBUTION_SERIES_DND)) {
                            ASSERT(payload->DataSize == sizeof(SeriesDragPayload));
                            const SeriesDragPayload* dnd = (const SeriesDragPayload*)(payload->Data);
                            if (dnd->src_subplot >= 0 && dnd->src_subplot < PLOT_MAX_SUBPLOTS) {
                                plot_move_series(data, data->distributions.subplots[dnd->src_subplot], sp, dnd->key);
                            } else {
                                plot_add_series(data, sp, dnd->key);
                            }
                            fit_subplot[i] = true;
                        } else if (const ImGuiPayload* tl_payload = ImGui::AcceptDragDropPayload(TIMELINE_SERIES_DND)) {
                            // Something plotted over time, dropped here: the distribution of its
                            // values, whichever view of them (a mean, a spread) was dragged
                            ASSERT(tl_payload->DataSize == sizeof(SeriesDragPayload));
                            SeriesKey key = ((const SeriesDragPayload*)tl_payload->Data)->key;
                            key.variant = SeriesVariant_Values;
                            plot_add_series(data, sp, key);
                            fit_subplot[i] = true;
                        }
                        ImPlot::EndDragDropTarget();
                    }

                    if (ImPlot::IsPlotHovered()) {
                        ImPlotPoint plot_pos = ImPlot::GetPlotMousePos();
                        ImVec2 screen_pos = ImPlot::PlotToPixels(plot_pos);
                        ImVec2 p0 = {screen_pos.x, ImPlot::GetPlotPos().y};
                        ImVec2 p1 = {screen_pos.x, ImPlot::GetPlotPos().y + ImPlot::GetPlotSize().y};
                        ImPlot::PushPlotClipRect();
                        ImPlot::GetPlotDrawList()->AddLine(p0, p1, IM_COL32(255, 255, 255, 120));
                        ImPlot::PopPlotClipRect();
                    }

                    ImPlot::EndPlot();
                }

                if (remove_idx != -1) {
                    plot_remove_series(sp, remove_idx);
                }
            }
            ImPlot::EndSubplots();
        }

        // A series dragged out of a plot and dropped outside any of them is taken out of that plot
        if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
            const ImGuiPayload* payload = ImGui::GetDragDropPayload();
            if (payload && payload->IsDataType(DISTRIBUTION_SERIES_DND) && !ImGui::IsDragDropPayloadBeingAccepted()) {
                const SeriesDragPayload* dnd = (const SeriesDragPayload*)(payload->Data);
                if (dnd && dnd->src_subplot >= 0 && dnd->src_subplot < PLOT_MAX_SUBPLOTS) {
                    PlotSubplot& src = data->distributions.subplots[dnd->src_subplot];
                    plot_remove_series(src, plot_find_series(src, dnd->key));
                }
            }
        }
    }
    ImGui::End();
}

static void draw_debug_window(ApplicationState* data) {
    ASSERT(data);

    ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Debug", &data->show_debug_window)) {       
        task_system::ID tasks[256]; 
        size_t num_tasks = task_system::pool_running_tasks(tasks, ARRAY_SIZE(tasks));
        if (num_tasks > 0) {
            ImGui::Text("Running Pool Tasks:");
            for (size_t i = 0; i < num_tasks; ++i) {
                str_t lbl = task_system::task_label(tasks[i]);
                ImGui::Text("[%i]: %.*s", (int)i, (int)lbl.len, lbl.ptr);
            }
        }

        ImGuiID active = ImGui::GetActiveID();
        ImGuiID hover  = ImGui::GetHoveredID();
        ImGui::Text("Active ID: %u, Hover ID: %u", active, hover);

        ImGui::Text("Mouse Pos: (%.3f, %.3f)", ImGui::GetMousePos().x, ImGui::GetMousePos().y);
        ImGui::Text("Camera Position: (%g, %g, %g)", data->view.camera.position.x, data->view.camera.position.y, data->view.camera.position.z);
        ImGui::Text("Camera Orientation: (%g, %g, %g, %g)", data->view.camera.orientation.x, data->view.camera.orientation.y, data->view.camera.orientation.z, data->view.camera.orientation.w);

        mat4_t P = data->view.param.matrix.curr.proj;
        ImGui::Text("proj_matrix:");
        ImGui::Text("[%g %g %g %g]", P.col[0].x, P.col[0].y, P.col[0].z, P.col[0].w);
        ImGui::Text("[%g %g %g %g]", P.col[1].x, P.col[1].y, P.col[1].z, P.col[1].w);
        ImGui::Text("[%g %g %g %g]", P.col[2].x, P.col[2].y, P.col[2].z, P.col[2].w);
        ImGui::Text("[%g %g %g %g]", P.col[3].x, P.col[3].y, P.col[3].z, P.col[3].w);
    }
    ImGui::End();
}

// Opens the script reference at topic (a procedure name, alias or heading anchor). Anything else is searched for,
// and an empty topic puts the focus in the search box (which always takes the focus). Without take_focus the keyboard
// focus stays where it is, so lookups from the script editor keep its caret.
static void open_script_reference(ApplicationState* state, str_t topic, bool take_focus) {
    ASSERT(state);
    if (str_empty(topic)) {
        script_reference::focus_search();
        take_focus = true;
    } else if (!script_reference::show(topic)) {
        script_reference::search(topic);
    }
    state->show_script_reference_window = true;
    script_reference::reveal(take_focus);
}

static void draw_script_reference_window(ApplicationState* state) {
    ASSERT(state);
    const script_reference::Action action = script_reference::draw_window(&state->show_script_reference_window);
    if (!str_empty(action.insert_code)) {
        // Examples go in as whole lines at the cursor of the script editor
        script_editor::insert_lines_at_cursor(state->editor, action.insert_code);
        state->show_script_window = true;
    }
}

static void draw_script_editor_window(ApplicationState* state) {
    ASSERT(state);

    ImGui::SetNextWindowSize({300,200}, ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Script Editor", &state->show_script_window, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
        ImGui::SetWindowSize(ImVec2(800, 600), ImGuiCond_FirstUseEver);
        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("File")) {
                char path_buf[1024] = "";
                if (ImGui::MenuItem("Load")) {
                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Open, STR_LIT("txt"))) {
                        str_t txt = load_textfile(str_from_cstr(path_buf), frame_alloc);
                        std::string str(txt.ptr, txt.len);
                        state->editor.SetText(str);
                    }
                }
                if (ImGui::MenuItem("Save")) {
                    auto textToSave = state->editor.GetText();
                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, STR_LIT("txt"))) {
                        str_t path = str_t{path_buf, strnlen(path_buf, sizeof(path_buf))};
                        md_file_t file = {0};
                        if (md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                            md_file_write(file, textToSave.c_str(), textToSave.length());
                            md_file_close(&file);
                        } else {
                            VIAMD_LOG_ERROR("Failed to open file '%s' for saving script", path_buf);
                        }
                    }
                }
                if (ImGui::MenuItem("Export")) {
                    state->show_property_export_window = true;
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Edit")) {
                bool ro = state->editor.IsReadOnlyEnabled();
                if (ImGui::MenuItem("Read-only mode", nullptr, &ro))
                    state->editor.SetReadOnlyEnabled(ro);
                ImGui::Separator();

                const bool sel = state->editor.AnyCursorHasSelection();
                if (ImGui::MenuItem("Undo", "Ctrl-Z", nullptr, !ro && state->editor.CanUndo()))
                    state->editor.Undo();
                if (ImGui::MenuItem("Redo", "Ctrl-Y", nullptr, !ro && state->editor.CanRedo()))
                    state->editor.Redo();

                ImGui::Separator();

                if (ImGui::MenuItem("Copy", "Ctrl-C", nullptr, sel))
                    state->editor.Copy();
                if (ImGui::MenuItem("Cut", "Ctrl-X", nullptr, !ro && sel))
                    state->editor.Cut();
                if (ImGui::MenuItem("Delete", "Del", nullptr, !ro && sel))
                    state->editor.ReplaceTextInAllCursors("");
                if (ImGui::MenuItem("Paste", "Ctrl-V", nullptr, !ro && ImGui::GetClipboardText() != nullptr))
                    state->editor.Paste();

                ImGui::Separator();

                if (ImGui::MenuItem("Select all", "Ctrl-A"))
                    state->editor.SelectAll();
                if (ImGui::MenuItem("Find / Replace", "Ctrl-F"))
                    state->editor.OpenFindReplaceWindow();

                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Settings")) {
                if (ImGui::MenuItem("Dark palette"))
                    state->editor.SetPalette(TextEditor::GetDarkPalette());
                if (ImGui::MenuItem("Light palette"))
                    state->editor.SetPalette(TextEditor::GetLightPalette());
                if (ImGui::MenuItem("Retro blue palette"))
                    state->editor.SetPalette(script_editor::retro_blue_palette());
                ImGui::Separator();
                ImGui::ColorEdit4("Point Color",    state->script.point_color.elem);
                ImGui::ColorEdit4("Line Color",     state->script.line_color.elem);
                ImGui::ColorEdit4("Triangle Color", state->script.triangle_color.elem);
                ImGui::ColorEdit4("Text Color",     state->script.text_color.elem);
                ImGui::ColorEdit4("Text Bg Color",  state->script.text_bg_color.elem);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Help")) {
                if (ImGui::MenuItem("Script reference", "F1")) {
                    open_script_reference(state, {});
                }
                const std::string word = script_editor::word_at_cursor(state->editor);
                char label[128];
                if (word.empty()) snprintf(label, sizeof(label), "Look up word under cursor");
                else snprintf(label, sizeof(label), "Look up '%s'", word.c_str());
                if (ImGui::MenuItem(label, "F1", nullptr, !word.empty())) {
                    open_script_reference(state, {word.data(), word.size()}, false);
                }
                ImGui::EndMenu();
            }

            ImGui::EndMenuBar();
        }

        const ImVec2 content_size = ImGui::GetContentRegionAvail();
        const char* btn_text = "Evaluate";
        const ImVec2 label_size = ImGui::CalcTextSize(btn_text, NULL, true) * ImVec2(1.4, 1.0);
        const ImVec2 btn_size = ImGui::CalcItemSize(ImVec2(0,0), label_size.x + ImGui::GetStyle().FramePadding.x * 2.0f, label_size.y + ImGui::GetStyle().FramePadding.y * 2.0f);
        const ImVec2 text_size(content_size - ImVec2(0, btn_size.y + ImGui::GetStyle().ItemSpacing.y));

        // Shift+Enter evaluates the script. The editor binds it to "insert line above", so claim it before the editor
        // sees it (shortcut routes are resolved from the previous frame, hence the focus from the last draw).
        bool eval = false;
        if (ImGui::Shortcut(KEY_SCRIPT_EVALUATE, ImGuiInputFlags_RouteGlobal | ImGuiInputFlags_RouteOverFocused)) {
            eval = true;
        }

        // While a visualization is hovered, the mouse wheel steps through its elements instead of scrolling
        const script_editor::Marker* prev_hovered = script_editor::markers_hovered(&state->editor_markers);
        ImGuiWindowFlags editor_flags = ImGuiWindowFlags_NoMove | ImGuiWindowFlags_HorizontalScrollbar;
        if (prev_hovered && prev_hovered->type == script_editor::MarkerType_Visualization) {
            editor_flags |= ImGuiWindowFlags_NoScrollWithMouse;
        }

        if (state->editor.Render("TextEditor", text_size, ImGuiChildFlags_None, editor_flags)) {
            state->script.compile_ir = true;
            state->script.time_since_last_change = 0;
        }
        const bool editor_hovered = ImGui::IsItemHovered();
        state->editor_focused = script_editor::has_focus_after_render();
        const script_editor::Marker* hovered_marker = script_editor::markers_update(&state->editor_markers, &state->editor, editor_hovered);

        // F1 looks up the identifier under the mouse in the script reference, or else the one at the text cursor.
        // Hovering does not require focus: opening the reference takes the focus, and a second F1 over another
        // identifier would otherwise be ignored until the editor is clicked again.
        if ((editor_hovered || state->editor_focused) && ImGui::IsKeyPressed(ImGuiKey_F1, false)) {
            std::string word;
            if (editor_hovered) {
                word = script_editor::word_at_mouse(state->editor, ImGui::GetMousePos());
            }
            if (word.empty() && state->editor_focused) {
                word = script_editor::word_at_cursor(state->editor);
            }
            open_script_reference(state, {word.data(), word.size()}, false);
        }

        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + content_size.x - btn_size.x);
        ImGui::SetCursorPosY(ImGui::GetCursorPosY() + ImGui::GetStyle().ItemSpacing.y);

        const bool valid = md_script_ir_valid(state->script.ir) && run_num_frames(state) > 0;  
        if (!valid) ImGui::PushDisabled();
        if (ImGui::Button(btn_text, btn_size)) {
            eval = true;
        }
        if (!valid) ImGui::PopDisabled();

        if (eval && valid) {
            state->script.eval_init = true;
        }

        if (editor_hovered) {
            md_bitfield_clear(&state->selection.highlight_mask);
        }

        if (hovered_marker) {
            if (!str_empty(hovered_marker->text)) {
                ImGui::BeginTooltip();
                ImGui::Text(STR_FMT, STR_ARG(hovered_marker->text));
                if (state->script.sub_idx != -1) {
                    ImGui::Text("Currently inspected idx: %i", state->script.sub_idx + 1);
                }
                ImGui::EndTooltip();
            }
            if (hovered_marker->atoms || hovered_marker->payload) {
                if (hovered_marker->atoms) {
                    // Errors and warnings: the atoms the message is about
                    md_bitfield_copy(&state->selection.highlight_mask, hovered_marker->atoms);
                }
                else if (hovered_marker->type == script_editor::MarkerType_Visualization) {
                    // Clear hovered property
                    script_set_hovered_property(state, STR_LIT(""));
                    if (md_script_ir_valid(state->script.ir)) {
                        const md_script_vis_payload_o* payload = hovered_marker->payload;
                        str_t payload_ident = md_script_payload_ident(payload);
                        int payload_dim  = md_script_payload_dim(payload); 

                        if (payload_dim > 1) {
                            int delta = (int)ImGui::GetIO().MouseWheel;
                            if (ImGui::IsKeyDown(ImGuiMod_Shift)) {
                                delta *= 10;
                            }
                            state->script.sub_idx += delta;
                            state->script.sub_idx = CLAMP(state->script.sub_idx, -1, (int)payload_dim - 1);
                        }

                        script_visualize_payload(state, payload, state->script.sub_idx, 0);
                        script_set_hovered_property(state, payload_ident, state->script.sub_idx);
                    }
                }

                bool lm_click = ImGui::IsMouseClicked(ImGuiMouseButton_Left);
                bool rm_click = ImGui::IsMouseClicked(ImGuiMouseButton_Right);
                if (ImGui::IsKeyDown(ImGuiMod_Shift) && (lm_click || rm_click)) {
                    SelectionOperator op = lm_click ? SelectionOperator::Or : SelectionOperator::AndNot;
                    modify_selection(state, &state->selection.highlight_mask, op);
                }
            }
        } else {
            state->script.sub_idx = -1;
        }
    }
    ImGui::End();
}

static bool export_xvg(const float* column_data[], const char* column_labels[], size_t num_columns, size_t num_rows, str_t filename) {
    ASSERT(column_data);
    ASSERT(column_labels);

    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        VIAMD_LOG_ERROR("Failed to open file '" STR_FMT "' to write data.", STR_ARG(filename));
        return false;
    }    

    time_t t;
    struct tm* info;
    time(&t);
    info = localtime(&t);

    // Print Header
    md_file_printf(file, "# This file was created %s", asctime(info));
    md_file_printf(file, "# Created by:\n");
    md_file_printf(file, "# VIAMD \n");

    // Print Legend Meta
    md_file_printf(file, "@    title \"VIAMD Properties\"\n");
    md_file_printf(file, "@    xaxis  label \"Time\"\n");
    md_file_printf(file, "@ TYPE xy\n");
    md_file_printf(file, "@ view 0.15, 0.15, 0.75, 0.85\n");
    md_file_printf(file, "@ legend on\n");
    md_file_printf(file, "@ legend box on\n");
    md_file_printf(file, "@ legend loctype view\n");
    md_file_printf(file, "@ legend 0.78, 0.8\n");
    md_file_printf(file, "@ legend length %i\n", num_columns);

    for (size_t j = 0; j < num_columns; ++j) {
        md_file_printf(file, "@ s%zu legend \"%s\"\n", j, column_labels[j]);
    }

    for (size_t i = 0; i < num_rows; ++i) {
        for (size_t j = 0; j < num_columns; ++j) {
            md_file_printf(file, "%12.6f ", column_data[j][i]);
        }
        md_file_printf(file, "\n");
    }

    md_file_close(&file);
    VIAMD_LOG_SUCCESS("Successfully exported XVG file to '%.*s'", (int)filename.len, filename.ptr);
    return true;
}

static bool export_csv(const float* column_data[], const char* column_labels[], size_t num_columns, size_t num_rows, str_t filename) {
    ASSERT(column_data);
    ASSERT(column_labels);

    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        VIAMD_LOG_ERROR("Failed to open file '%.*s' to write data.", (int)filename.len, filename.ptr);
        return false;
    }
    
    for (size_t i = 0; i < num_columns; ++i) {
        md_file_printf(file, "%s,", column_labels[i]);
    }
    md_file_printf(file, "\n");

    for (size_t i = 0; i < num_rows; ++i) {
        for (size_t j = 0; j < num_columns; ++j) {
            md_file_printf(file, "%.6g,", column_data[j][i]);
        }
        md_file_printf(file, "\n");
    }

    md_file_close(&file);
    VIAMD_LOG_SUCCESS("Successfully exported CSV file to '%.*s'", (int)filename.len, filename.ptr);
    return true;
}

static bool export_cube(const ApplicationState& data, const md_attribute_t* attr, const md_script_vis_payload_o* vis_payload, str_t filename) {
    // @NOTE: First we need to extract some meta data for the cube format, we need the atom indices/bits for any SDF
    // And the origin + extent of the volume in spatial coordinates (Ångström)

    const float* volume_values = (const float*)md_attribute_view(attr, MD_ATTRIBUTE_TYPE_F32, 1, 3);
    if (!volume_values) {
        VIAMD_LOG_ERROR("Export Cube: The property to be exported did not exist");
        return false;
    }

    if (!vis_payload) {
        VIAMD_LOG_ERROR("Export Cube: Missing input visualization data");
        return false;
    }

	md_temp_scope_t temp = md_temp_begin_in(frame_alloc);
	defer { md_temp_end(temp); };

    const md_system_t& sys = data.mold.sys;
    md_system_state_t state = { .alloc = temp.arena };
    if (!extract_frame(&data, 0, &state)) {
        return false;
    }

    bool result = false;
    md_script_vis_t vis = { 0 };
    md_script_vis_init(&vis, temp.arena);
        
    if (md_script_ir_valid(data.script.eval_ir)) {
        md_script_vis_ctx_t ctx = {
            .ir = data.script.eval_ir,
            .sys = &data.mold.sys,
            .state = &data.mold.state,
        };
        result = md_script_vis_eval_payload(&vis, vis_payload, 0, &ctx, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_SDF);
    }

    if (result == true) {
        md_file_t file = { 0 };
        if (!md_file_open(&file, filename, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
            VIAMD_LOG_ERROR("Failed to open file '%.*s' in order to write to it.", (int)filename.len, filename.ptr);
            return false;
        }

        // Two comment lines
        md_file_printf(file, "EXPORTED DENSITY VOLUME FROM VIAMD, UNITS IN BOHR\n");
        md_file_printf(file, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");

        if (md_array_size(vis.sdf.structures) > 0) {
            const float angstrom_to_bohr = (float)(1.0 / 0.529177210903);

            // transformation matrix from world to volume
            mat4_t M = vis.sdf.matrices[0];
            const md_bitfield_t* bf = &vis.sdf.structures[0];
            const int num_atoms = (int)md_bitfield_popcount(bf);
            const int vol_dim[3] = {(int)attr->format.shape[0], (int)attr->format.shape[1], (int)attr->format.shape[2]};
            const float* values = volume_values;
            const double extent = vis.sdf.extent * 2.0 * angstrom_to_bohr;
            const double voxel_ext[3] = {
                (double)extent / (double)vol_dim[0],
                (double)extent / (double)vol_dim[1],
                (double)extent / (double)vol_dim[2],
            };

            const double half_ext = extent * 0.5;

            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", -num_atoms, -half_ext, -half_ext, -half_ext);
            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol_dim[0], voxel_ext[0], 0.0, 0.0);
            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol_dim[1], 0.0, voxel_ext[1], 0.0);
            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol_dim[2], 0.0, 0.0, voxel_ext[2]);

            const float scl = angstrom_to_bohr;
            M = mat4_mul(mat4_scale(scl, scl, scl), M);

            size_t beg_bit = bf->beg_bit;
            size_t end_bit = bf->end_bit;
            while ((beg_bit = md_bitfield_scan(bf, beg_bit, end_bit)) != 0) {
                size_t i = beg_bit - 1;
                vec3_t coord = md_state_coord(&state, i);
                coord = mat4_mul_vec3(M, coord, 1.0f);
                int anum     = md_atom_atomic_number(&sys.atom, i);
                float charge = (float)anum;
                md_file_printf(file, "%5i %12.6f %12.6f %12.6f %12.6f\n", anum, charge, coord.x, coord.y, coord.z);
            }

            // This entry somehow relates to the number of densities
            md_file_printf(file, "%5i %5i\n", 1, 1);

            // Write density data
            int count = 0;
            for (int x = 0; x < vol_dim[0]; ++x) {
                for (int y = 0; y < vol_dim[1]; ++y) {
                    for (int z = 0; z < vol_dim[2]; ++z) {
                        int idx = z * vol_dim[0] * vol_dim[1] + y * vol_dim[0] + x;
                        float val = values[idx];
                        md_file_printf(file, " %12.6E", val);
                        if (++count % 6 == 0) md_file_printf(file, "\n");
                    }
                }
            }
        }

        md_file_close(&file);
    } else {
        VIAMD_LOG_ERROR("Failed to visualize volume for export.");
        return false;
    }

    return true;
}

#define APPEND_BUF(buf, len, fmt, ...) (len += snprintf(buf + len, MAX(0, (int)sizeof(buf) - len), fmt, ##__VA_ARGS__) + 1)

static void draw_property_export_window(ApplicationState* data) {
    ASSERT(data);

    struct ExportFormat {
        str_t lbl;
        str_t ext;
    };

    ExportFormat table_formats[] {
        {STR_INIT("XVG"), STR_INIT("xvg")},
        {STR_INIT("CSV"), STR_INIT("csv")}
    };

    ExportFormat volume_formats[] {
        {STR_INIT("Gaussian Cube"), STR_INIT("cube")},
    };

    enum ExportKind {
        ExportKind_Temporal = 0,
        ExportKind_Distribution,
        ExportKind_Volume,
    };

    if (ImGui::Begin("Property Export", &data->show_property_export_window)) {
        static int kind = ExportKind_Temporal;
        static SeriesKey selected = {};
        static int table_format  = 0;
        static int volume_format = 0;
        static int bins_exp      = 7;   // 128 bins

        if (task_system::task_is_running(data->tasks.evaluate_full)) {
            ImGui::Text("The properties are currently being evaluated, please wait...");
            ImGui::End();
            return;
        }

        md_temp_scope_t temp = md_temp_begin_in(frame_alloc);
        defer { md_temp_end(temp); };
        md_allocator_i* alloc = md_temp_allocator(temp);

        // What there is to export of the chosen kind: the same series the plots offer
        md_array(SeriesKey) candidates = 0;
        {
            const uint32_t script_kinds =
                kind == ExportKind_Temporal     ? (uint32_t)MD_SCRIPT_PROPERTY_FLAG_TEMPORAL :
                kind == ExportKind_Distribution ? (uint32_t)(MD_SCRIPT_PROPERTY_FLAG_TEMPORAL | MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) :
                                                  (uint32_t)MD_SCRIPT_PROPERTY_FLAG_VOLUME;
            auto push_script = [&](SeriesSource source) {
                series_for_each_script_property(data, source, script_kinds, [&](const SeriesKey& key) {
                    md_array_push(candidates, key, alloc);
                    const bool temporal = md_script_ir_property_flags(data->script.eval_ir, series_script_ident(key)) & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL;
                    if (kind == ExportKind_Distribution && temporal && script_property_population(data, key) > 1) {
                        SeriesKey agg = key;
                        agg.variant = SeriesVariant_Aggregate;
                        md_array_push(candidates, agg, alloc);
                    }
                });
            };
            push_script(SeriesSource_Script);
            if (data->timeline.filter.enabled && kind != ExportKind_Temporal) {
                push_script(SeriesSource_ScriptFiltered);
            }
            if (kind != ExportKind_Volume) {
                str_t groups[64];
                const size_t num_groups = MIN(system_series_groups(groups, ARRAY_SIZE(groups), data), ARRAY_SIZE(groups));
                for (size_t g = 0; g < num_groups; ++g) {
                    const size_t num = system_series_members(nullptr, 0, data, groups[g]);
                    str_t* paths = md_temp_alloc_array(temp, str_t, num + 1);
                    system_series_members(paths, num, data, groups[g]);
                    for (size_t i = 0; i < num; ++i) {
                        md_array_push(candidates, series_key(SeriesSource_System, paths[i]), alloc);
                    }
                }
            }
        }

        bool selected_valid = false;
        for (size_t i = 0; i < md_array_size(candidates); ++i) {
            if (series_key_equal(candidates[i], selected)) { selected_valid = true; break; }
        }
        if (!selected_valid) {
            selected = md_array_size(candidates) ? candidates[0] : SeriesKey{};
            selected_valid = md_array_size(candidates) > 0;
        }

        ImGui::PushItemWidth(200);
        ImGui::Combo("Data Type", &kind, "Temporal\0Distribution\0Density Volume\0");
        ImGui::Separator();

        char selected_label[96] = "";
        if (selected_valid) {
            series_label(selected_label, sizeof(selected_label), data, selected);
        }
        if (ImGui::BeginCombo("Property", selected_label)) {
            for (size_t i = 0; i < md_array_size(candidates); ++i) {
                char label[96];
                series_label(label, sizeof(label), data, candidates[i]);
                ImGui::PushID((int)i);
                if (ImGui::Selectable(label, series_key_equal(candidates[i], selected))) {
                    selected = candidates[i];
                }
                if (ImGui::IsItemHovered() && candidates[i].source == SeriesSource_System) {
                    ImGui::SetTooltip("%s", candidates[i].path);
                }
                ImGui::PopID();
            }
            ImGui::EndCombo();
        }
        if (md_array_size(candidates) == 0) {
            ImGui::TextDisabled("Nothing of this kind, evaluate the script or load an energy file");
        }

        if (kind == ExportKind_Distribution) {
            char bins_fmt[32];
            snprintf(bins_fmt, sizeof(bins_fmt), "%d", 1 << bins_exp);
            ImGui::SliderInt("Num Bins", &bins_exp, 5, 10, bins_fmt);
        }

        str_t file_extension = {};
        if (kind == ExportKind_Volume) {
            if (ImGui::BeginCombo("File Format", volume_formats[volume_format].lbl.ptr)) {
                for (int i = 0; i < (int)ARRAY_SIZE(volume_formats); ++i) {
                    if (ImGui::Selectable(volume_formats[i].lbl.ptr, volume_format == i)) {
                        volume_format = i;
                    }
                }
                ImGui::EndCombo();
            }
            file_extension = volume_formats[volume_format].ext;
        } else {
            if (ImGui::BeginCombo("File Format", table_formats[table_format].lbl.ptr)) {
                for (int i = 0; i < (int)ARRAY_SIZE(table_formats); ++i) {
                    if (ImGui::Selectable(table_formats[i].lbl.ptr, table_format == i)) {
                        table_format = i;
                    }
                }
                ImGui::EndCombo();
            }
            file_extension = table_formats[table_format].ext;
        }

        if (!selected_valid) ImGui::PushDisabled();
        const bool export_clicked = ImGui::Button("Export");
        if (!selected_valid) ImGui::PopDisabled();
        ImGui::PopItemWidth();

        char path_buf[1024];
        if (export_clicked && application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, file_extension)) {
            const str_t path = {path_buf, strnlen(path_buf, sizeof(path_buf))};
            char label[96];
            series_label(label, sizeof(label), data, selected);

            if (kind == ExportKind_Volume) {
                SeriesVolumeView view;
                if (!series_resolve_volume(&view, data, selected)) {
                    VIAMD_LOG_ERROR("'%s' has no volume to export", label);
                } else if (str_eq(file_extension, STR_LIT("cube")) && export_cube(*data, view.attr, view.vis_payload, path)) {
                    VIAMD_LOG_SUCCESS("Successfully exported property '%s' to '" STR_FMT "'", label, STR_ARG(path));
                }
            } else {
                // Columns in the units the plots show them in, so a column and its axis label agree
                md_array(const float*) column_data   = 0;
                md_array(str_t)        column_labels = 0;
                md_array(str_t)        legends       = 0;
                str_t x_label = {};
                str_t y_label = str_from_cstr(label);
                size_t num_rows = 0;
                bool ok = false;

                auto push_columns = [&](int dim, auto value) {
                    for (int k = 0; k < dim; ++k) {
                        float* col = md_temp_alloc_array(temp, float, num_rows);
                        for (size_t r = 0; r < num_rows; ++r) {
                            col[r] = value(r, k);
                        }
                        md_array_push(column_data, col, alloc);
                        if (dim > 1) {
                            const str_t legend = str_printf(alloc, "%s[%i]", label, k + 1);
                            md_array_push(legends, legend, alloc);
                            md_array_push(column_labels, legend, alloc);
                        } else {
                            md_array_push(column_labels, y_label, alloc);
                        }
                    }
                };

                if (kind == ExportKind_Temporal) {
                    SeriesTemporalView view;
                    if (series_resolve_temporal(&view, data, selected)) {
                        num_rows = (size_t)view.num_samples;
                        char time_unit_str[32] = "";
                        md_unit_print(time_unit_str, sizeof(time_unit_str), data->timeline.time_unit);
                        x_label = md_unit_is_none(data->timeline.time_unit) ? STR_LIT("Frame") : str_printf(alloc, "Time (%s)", time_unit_str);
                        if (view.unit_str[0]) y_label = str_printf(alloc, "%s (%s)", label, view.unit_str);

                        float* x = md_temp_alloc_array(temp, float, num_rows);
                        MEMCPY(x, view.x, num_rows * sizeof(float));
                        md_array_push(column_data, x, alloc);
                        md_array_push(column_labels, x_label, alloc);
                        push_columns(view.dim, [&](size_t r, int k) { return (float)series_temporal_point(view, (int)r, k).y; });
                        ok = true;
                    }
                } else {
                    SeriesHistogramView view;
                    if (series_resolve_histogram(&view, data, selected, 1 << bins_exp)) {
                        num_rows = (size_t)view.num_bins;
                        x_label = str_copy(str_from_cstr(view.x_unit_str), alloc);
                        if (view.y_unit_str[0]) y_label = str_printf(alloc, "%s (%s)", label, view.y_unit_str);

                        // The centre of each bin, where the plot draws it
                        float* x = md_temp_alloc_array(temp, float, num_rows);
                        for (size_t r = 0; r < num_rows; ++r) {
                            x[r] = (float)series_histogram_point(view, (int)r, 0).x;
                        }
                        md_array_push(column_data, x, alloc);
                        md_array_push(column_labels, x_label, alloc);
                        push_columns(view.dim, [&](size_t r, int k) { return (float)series_histogram_point(view, (int)r, k).y; });
                        ok = true;
                    }
                }

                if (!ok) {
                    VIAMD_LOG_ERROR("'%s' has nothing to export right now", label);
                } else {
                    str_t out_str = {};
                    if (str_eq(file_extension, STR_LIT("xvg"))) {
                        str_t header = md_xvg_format_header(str_from_cstr(label), x_label, y_label, md_array_size(legends), legends, alloc);
                        out_str = md_xvg_format(header, md_array_size(column_data), num_rows, column_data, alloc);
                    } else if (str_eq(file_extension, STR_LIT("csv"))) {
                        out_str = md_csv_write_to_str(column_data, column_labels, md_array_size(column_data), num_rows, alloc);
                    }
                    md_file_t file = {0};
                    if (!str_empty(out_str) && md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                        md_file_write(file, out_str.ptr, out_str.len);
                        md_file_close(&file);
                        VIAMD_LOG_SUCCESS("Successfully exported property '%s' to '" STR_FMT "'", label, STR_ARG(path));
                    } else {
                        VIAMD_LOG_ERROR("Failed to write '" STR_FMT "'", STR_ARG(path));
                    }
                }
            }
        }
    }
    ImGui::End();
}

md_array(int) extract_bitfield_indices(const md_bitfield_t* bf, md_allocator_i* alloc) {
    size_t count = md_bitfield_popcount(bf);
    md_array(int) indices = md_array_create(int, count, alloc);
    md_bitfield_iter_extract_indices(indices, count, md_bitfield_iter_create(bf));
    return indices;
}

static void xyz_write_frame(md_file_t file, const int* atomic_nr, const md_system_state_t* state, const int32_t* atom_indices, size_t num_atoms) {
    md_file_printf(file, "%zu\n", num_atoms);
    md_file_printf(file, "XYZ format exported from VIAMD\n");
    for (size_t i = 0; i < num_atoms; ++i) {
        int idx = atom_indices ? atom_indices[i] : (int)i;
        md_file_printf(file, "%-2d %12.6f %12.6f %12.6f\n", atomic_nr[idx], state->xyz[idx].x, state->xyz[idx].y, state->xyz[idx].z);
    }
}

void draw_structure_export_window(ApplicationState* data) {
    ASSERT(data);

    if (ImGui::Begin("Structure Export", &data->structure_export.show_window)) {
        md_temp_scope_t temp = md_temp_begin_in(frame_alloc);
        defer { md_temp_end(temp); };

        const md_system_t* sys = &data->mold.sys;
        const md_system_state_t* state = &data->mold.state;
        const bool traj = run_num_frames(data) > 0;
        auto& struct_exp = data->structure_export;

        static const char* atom_mask_options[] = {
            "Active Selection",
            "Visible Atoms",
            "Query",
            "All Atoms",
        };

        static const char* frame_mask_options[] = {
            "All Frames",
            "Current Frame",
            "Active Frame Range"
        };

        static const char* file_formats[] = {
            "xyz",
            "pdb"
        };

        //ImGui::PushItemWidth(200);

        struct_exp.selected_atom_filter = CLAMP(struct_exp.selected_atom_filter, 0, (int)ARRAY_SIZE(atom_mask_options) - 1);
        ImGui::Combo("Atoms to Export", &struct_exp.selected_atom_filter, atom_mask_options, (int)ARRAY_SIZE(atom_mask_options));

        if (struct_exp.selected_atom_filter == 2) { // Query
            auto& query = struct_exp.query;
            bool valid = query.is_valid;

            if (!md_bitfield_validate(&query.mask)) {
                md_bitfield_init(&query.mask, persistent_alloc);
            }
            
            if (!valid) ImGui::PushInvalid();
            ImGui::InputQuery("Query", query.buf, sizeof(query.buf), query.is_valid, query.error);
            if (!valid) ImGui::PopInvalid();
            
            query.requires_evaluation |= ImGui::IsItemEdited();
            bool preview = ImGui::IsItemFocused() || ImGui::IsItemHovered();
            
            if (query.requires_evaluation) {
                md_bitfield_clear(&query.mask);
                query.is_valid = md_filter(
                    &query.mask,
                    str_from_cstr(query.buf),
                    &data->mold.sys,
                    &data->mold.state,
                    data->script.ir,
                    &query.is_dynamic,
                    query.error,
                    sizeof(query.error)
                );
                query.requires_evaluation = false;
            }

            if (ImGui::IsWindowHovered()) {
                // Clear highlight mask
                md_bitfield_clear(&data->selection.highlight_mask);
            }

            if (query.is_valid && preview) {
                md_bitfield_copy(&data->selection.highlight_mask, &query.mask);
            }
        }

        if (traj) {
            struct_exp.selected_traj_filter = CLAMP(struct_exp.selected_traj_filter, 0, (int)ARRAY_SIZE(frame_mask_options) - 1);
            ImGui::Combo("Frames to Export", &struct_exp.selected_traj_filter, frame_mask_options, (int)ARRAY_SIZE(frame_mask_options));
        }

        struct_exp.selected_file_format = CLAMP(struct_exp.selected_file_format, 0, (int)ARRAY_SIZE(file_formats) - 1);
        ImGui::Combo("File Format", &struct_exp.selected_file_format, file_formats, (int)ARRAY_SIZE(file_formats));

        //ImGui::PopItemWidth();

        auto extract_atom_indices = [data](int* out_atom_indices, size_t atom_index_cap, const md_system_state_t* state) -> size_t {
            switch (data->structure_export.selected_atom_filter) {
                case 0: { // Active Selection
                    return md_bitfield_iter_extract_indices(out_atom_indices, atom_index_cap, md_bitfield_iter_create(&data->selection.selection_mask));
                } break;
                case 1: { // Visible Atoms
                    return md_bitfield_iter_extract_indices(out_atom_indices, atom_index_cap, md_bitfield_iter_create(&data->representation.visibility_mask));
                } break;
                case 2: { // Query
                    if (data->structure_export.query.is_dynamic) {
                        // IMPORTANT: Clear the previous mask before re-evaluating a dynamic query.
                        // Otherwise atoms that drop out in subsequent frames remain set from earlier frames.
                        md_bitfield_clear(&data->structure_export.query.mask);
                        data->structure_export.query.is_valid = md_filter(
                            &data->structure_export.query.mask,
                            str_from_cstr(data->structure_export.query.buf),
                            &data->mold.sys,
                            state,
                            data->script.ir,
                            &data->structure_export.query.is_dynamic,
                            data->structure_export.query.error,
                            sizeof(data->structure_export.query.error)
                        );
                    }
                    if (!data->structure_export.query.is_valid) {
                        VIAMD_LOG_ERROR("Cannot export structure, the query expression is invalid: " STR_FMT, STR_ARG(str_from_cstr(data->structure_export.query.error)));
                        return 0;
                    } else {
                        return md_bitfield_iter_extract_indices(out_atom_indices, atom_index_cap, md_bitfield_iter_create(&data->structure_export.query.mask));
                    }
                } break;
                case 3: { // All Atoms
                    for (size_t i = 0; i < data->mold.sys.atom.count && i < atom_index_cap; ++i) {
                        out_atom_indices[i] = (int)i;
                    }
                    return data->mold.sys.atom.count;
                } break;
                default: {
                    VIAMD_LOG_ERROR("Invalid atom filter selection.");
                } break;
            }
            return 0;
        };

        if (ImGui::Button("Export")) {
            md_array(int32_t) frame_indices = 0;
            if (traj) {
                switch (struct_exp.selected_traj_filter) {
                case 0: { // All Frames
                    size_t num_frames = run_num_frames(data);
                    for (size_t i = 0; i < num_frames; ++i) {
                        md_array_push(frame_indices, (int32_t)i, frame_alloc);
                    }
                } break;
                case 1: { // Current Frame
                    int32_t current_frame = (int)(data->animation.frame + 0.5); // Round to nearest
                    md_array_push(frame_indices, current_frame, frame_alloc);
                } break;
                case 2: { // Active Frame Range
                    int beg = data->timeline.filter.beg_frame;
                    int end = data->timeline.filter.end_frame;
                    for (int i = beg; i <= end; ++i) {
                        md_array_push(frame_indices, (int32_t)i, frame_alloc);
                    }
                } break;
                default: {
                    VIAMD_LOG_ERROR("Invalid trajectory filter selection.");
                } break;
                }
            }

            char path_buf[1024] = { 0 };
            str_t ext = str_from_cstr(file_formats[struct_exp.selected_file_format]);
            if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, ext)) {
                str_t path = str_from_cstrn(path_buf, sizeof(path_buf));
                md_file_t file = {0};
                if (md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                    // Write header for file format
                    switch (struct_exp.selected_file_format) {
                        case 0: { // XYZ
                            // No header needed
                        } break;
                        case 1: { // PDB
                            md_file_printf(file, "HEADER    EXPORTED FROM VIAMD\n");
                            md_file_printf(file, "REMARK    Exported from VIAMD\n");
                            
                            // Write CRYST1 record with unit cell if available
                            if (traj) {
                                // Coordinates are not needed, only the cell of the first exported
                                // frame, which is delivered on the state.
                                md_system_state_t frame_state = {};
                                int first_frame = md_array_size(frame_indices) > 0 ? frame_indices[0] : 0;
                                if (extract_frame(data, first_frame, &frame_state)) {
                                    if (frame_state.unitcell.flags != 0) {
                                        double a, b, c, alpha, beta, gamma;
                                        md_unitcell_extract_extent_angles(&a, &b, &c, &alpha, &beta, &gamma, &frame_state.unitcell);
                                        md_file_printf(file, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n", a, b, c, alpha, beta, gamma);
                                    }
                                }
                            } else if (state->unitcell.flags != 0) {
                                double a, b, c, alpha, beta, gamma;
                                md_unitcell_extract_extent_angles(&a, &b, &c, &alpha, &beta, &gamma, &state->unitcell);
                                md_file_printf(file, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n", a, b, c, alpha, beta, gamma);
                            }
                        } break;
                        default: {
                            VIAMD_LOG_ERROR("Unsupported export format: '" STR_FMT "'", STR_ARG(ext));
                        } break;
                    }

                    // Write body
                    md_array(int32_t) atom_indices = 0;
                    md_array_ensure(atom_indices, sys->atom.count, frame_alloc);

                    md_array(int) atomic_numbers = md_array_create(int, sys->atom.count, frame_alloc);
                    for (size_t i = 0; i < sys->atom.count; ++i) {
                        atomic_numbers[i] = md_atom_atomic_number(&sys->atom, i);
                    }

                    size_t num_frames = md_array_size(frame_indices);
                    if (num_frames > 0) {
						md_system_state_t frame_state = { .alloc = frame_alloc };
                        md_system_state_init(&frame_state, sys->atom.count);

                        // Many frames: one context for all of them, so the run's files stay open.
                        const str_t paths[] = { STR_INIT("atom/position"), STR_INIT("unitcell") };
                        md_system_extract_t* ex = md_system_extract_begin(sys, str_from_cstr(data->mold.run), paths, ARRAY_SIZE(paths), md_get_heap_allocator());
                        defer { md_system_extract_end(ex); };

                        for (size_t f = 0; f < num_frames; ++f) {
                            int frame_idx = frame_indices[f];
                            if (!ex || !md_system_extract_frame(ex, frame_idx, &frame_state)) {
                                VIAMD_LOG_ERROR("Failed to load frame %d from trajectory for structure export.", frame_idx);
                                continue;
                            }

                            size_t num_atoms = extract_atom_indices(atom_indices, md_array_capacity(atom_indices), &frame_state);
                            
                            if (struct_exp.selected_file_format == 0) { // XYZ
                                xyz_write_frame(file, atomic_numbers, &frame_state, atom_indices, num_atoms);
                            } else if (struct_exp.selected_file_format == 1) { // PDB
                                int model_num = (num_frames > 1) ? (int)(f + 1) : 0;
                                md_pdb_system_write_state_to_file(file, sys, &frame_state, atom_indices, num_atoms, model_num);
                            }
                        }
                    }
                    else {
                        // Single frame from system
                        size_t num_atoms = extract_atom_indices(atom_indices, md_array_capacity(atom_indices), &data->mold.state);
                        
                        if (struct_exp.selected_file_format == 0) { // XYZ
                            xyz_write_frame(file, atomic_numbers, state, atom_indices, num_atoms);
                        } else if (struct_exp.selected_file_format == 1) { // PDB
                            md_pdb_system_write_state_to_file(file, sys, state, atom_indices, num_atoms, 0);
                        }
                    }

                    VIAMD_LOG_INFO("Successfully exported structure to: '" STR_FMT "'", STR_ARG(path));
                    md_file_close(&file);
                } else {
                    VIAMD_LOG_ERROR("Failed to open file '" STR_FMT "' for writing.", STR_ARG(path));
                }
            }
        }
    }
    ImGui::End();
}

static void update_md_buffers(ApplicationState* data) {
    ASSERT(data);
    const md_system_t& sys = data->mold.sys;
	const md_system_state_t& state = data->mold.state;

    if (sys.atom.count == 0) return;

    if (data->mold.dirty_gpu_buffers) {
        data->mold.gpu_buffers_version += 1;
    }

    if (data->mold.dirty_gpu_buffers & MolBit_DirtyPosition) {
        vec3_t pbc_ext = { 0 };
        md_unitcell_diag_extract_float(pbc_ext.elem, &state.unitcell);
        md_gl_mol_set_atom_position(data->mold.gl_mol, 0, (uint32_t)state.num_atoms, state.xyz);
        if (!(data->mold.dirty_gpu_buffers & MolBit_ClearVelocity)) {
            md_gl_mol_compute_velocity(data->mold.gl_mol, pbc_ext.elem);
        }
#if EXPERIMENTAL_GFX_API
        md_gfx_structure_set_atom_position(data->mold.gfx_structure, state.xyz, (uint32_t)state.num_atoms, 0);
        md_gfx_structure_set_aabb(data->mold.gfx_structure, &data->mold.sys_aabb_min, &data->mold.sys_aabb_max);
#endif
    }

    if (data->mold.dirty_gpu_buffers & MolBit_ClearVelocity) {
        md_gl_mol_zero_velocity(data->mold.gl_mol);
    }

    if (data->mold.dirty_gpu_buffers & MolBit_ResetBackboneHistory) {
        md_gl_mol_reset_backbone_history(data->mold.gl_mol);
    }

    if (data->mold.dirty_gpu_buffers & MolBit_DirtyRadius) {
        md_temp_scope_t tmp = md_temp_begin_in(frame_alloc);
        defer { md_temp_end(tmp); };

        float* radii = (float*)md_vm_arena_push(frame_alloc, sys.atom.count * sizeof(float));
        md_atom_extract_radii(radii, 0, sys.atom.count, &sys.atom);

        md_gl_mol_set_atom_radius(data->mold.gl_mol, 0, (uint32_t)sys.atom.count, radii, 0);
#if EXPERIMENTAL_GFX_API
        md_gfx_structure_set_atom_radius(data->mold.gfx_structure, 0, (uint32_t)sys.atom.count, sys.atom.radius, 0);
#endif
    }

    if (data->mold.dirty_gpu_buffers & MolBit_DirtyFlags) {
        md_temp_scope_t tmp = md_temp_begin_in(frame_alloc);
        defer { md_temp_end(tmp); };

        uint8_t* flags = (uint8_t*)md_vm_arena_push(frame_alloc, sys.atom.count * sizeof(uint8_t));
        MEMSET(flags, 0, sys.atom.count * sizeof(uint8_t));

        {
            md_bitfield_iter_t it = md_bitfield_iter_create(&data->selection.highlight_mask);
            while (md_bitfield_iter_next(&it)) {
                uint64_t idx = md_bitfield_iter_idx(&it);
                flags[idx] |= AtomBit_Highlighted;
            }
        }
        {
            md_bitfield_iter_t it = md_bitfield_iter_create(&data->selection.selection_mask);
            while (md_bitfield_iter_next(&it)) {
                uint64_t idx = md_bitfield_iter_idx(&it);
                flags[idx] |= AtomBit_Selected;
            }
        }
        {
            md_bitfield_iter_t it = md_bitfield_iter_create(&data->representation.visibility_mask);
            while (md_bitfield_iter_next(&it)) {
                uint64_t idx = md_bitfield_iter_idx(&it);
                flags[idx] |= AtomBit_Visible;
            }
        }
        md_gl_mol_set_atom_flags(data->mold.gl_mol, 0, (uint32_t)sys.atom.count, flags, 0);
    }

    if (data->mold.dirty_gpu_buffers & MolBit_DirtyBonds) {
        md_gl_mol_set_bonds(data->mold.gl_mol, 0, (uint32_t)sys.bond.count, sys.bond.pairs, sizeof(md_atom_pair_t));
    }

    data->mold.dirty_gpu_buffers = 0;
}

void create_screenshot(str_t path) {
    md_temp_scope_t tmp = md_temp_begin_in(frame_alloc);
    defer { md_temp_end(tmp); };

    int viewport[4] = {};
    int fbo = 0;
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &fbo);

    int width  = viewport[2];
    int height = viewport[3];

    size_t bytes = width * height * sizeof(uint32_t);
    uint32_t* rgba = (uint32_t*)md_vm_arena_push(frame_alloc, bytes);
    defer { md_vm_arena_pop(frame_alloc, bytes); };

    GLenum read_buffer = fbo ? GL_COLOR_ATTACHMENT0 : GL_BACK;

    glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo);
    glReadBuffer(read_buffer);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, rgba);

    {
        // @NOTE: Swap Rows to flip image with respect to y-axis
        const uint32_t row_byte_size = width * sizeof(uint32_t);
        uint32_t* row_t = (uint32_t*)md_alloc(frame_alloc, row_byte_size);
        defer { md_free(frame_alloc, row_t, row_byte_size); };
        for (uint32_t i = 0; i < (uint32_t)height / 2; ++i) {
            uint32_t* row_a = rgba + i * width;
            uint32_t* row_b = rgba + (height - 1 - i) * width;
            if (row_a != row_b) {
                MEMCPY(row_t, row_a, row_byte_size);  // tmp = a;
                MEMCPY(row_a, row_b, row_byte_size);  // a = b;
                MEMCPY(row_b, row_t, row_byte_size);  // b = tmp;
            }
        }
    }

    str_t ext = {};
    extract_ext(&ext, path);

    if (str_eq_cstr_ignore_case(ext, "jpg")) {
        const int quality = 95;
        image_write_jpg(path, rgba, width, height, quality);
    } else if (str_eq_cstr_ignore_case(ext, "png")) {
        image_write_png(path, rgba, width, height);
    } else if (str_eq_cstr_ignore_case(ext, "bmp")) {
        image_write_bmp(path, rgba, width, height);
    } else {
        VIAMD_LOG_ERROR("Non supported file-extension '" STR_FMT "' when saving screenshot", STR_ARG(ext));
        return;
    }

    VIAMD_LOG_SUCCESS("Screenshot saved to: '" STR_FMT "'", STR_ARG(path));
}

// Length of the movie in seconds
static double movie_duration(const ApplicationState* state) {
    return (double)state->movie.duration;
}

// Sets the length of the movie and scales every time on its timeline with it, so the movie keeps its shape
static void movie_set_duration(ApplicationState* state, float new_len) {
    auto& m = state->movie;
    new_len = CLAMP(new_len, 0.01f, 3600.0f);
    const double old_len = (double)m.duration;
    m.duration = new_len;
    if (old_len <= 0.0 || (double)new_len == old_len) return;
    const double s = (double)new_len / old_len;
    MovieKeys keys = movie_keys_snapshot(state);
    movie_keys_scale_time(&keys, s);
    keys.duration = new_len;
    movie_keys_restore(state, keys);
    m.playhead = (float)(m.playhead * s);
    m.range_begin = (float)(m.range_begin * s);
    m.range_end = (float)(m.range_end * s);
    m.orbit_duration = (float)(m.orbit_duration * s);
}

// Keeps the trajectory anchors inside the movie, and moves an anchor out to a keyframe with a frame that is
// beyond it, as that keyframe takes its place
static void movie_clamp_anchors(ApplicationState* state) {
    auto& m = state->movie;
    for (size_t i = 0; i < md_array_size(m.keyframes); ++i) {
        const CameraKeyframe& k = m.keyframes[i];
        if (!k.use_frame) continue;
        if (k.time < m.traj_begin) m.traj_begin = (float)k.time;
        if (k.time > m.traj_end) m.traj_end = (float)k.time;
    }
    m.traj_begin = CLAMP(m.traj_begin, 0.0f, m.duration);
    m.traj_end = CLAMP(m.traj_end, m.traj_begin, m.duration);
}

static int movie_num_frames(const ApplicationState* state) {
    return (int)floor(movie_duration(state) * (double)state->movie.fps + 1.0e-6) + 1;
}

// The frames of the movie that are rendered, first to last. Frame i is at i / fps on the movie's timeline.
static void movie_render_range(const ApplicationState* state, int* first, int* last) {
    const auto& m = state->movie;
    movie_frame_range(movie_num_frames(state), (double)m.fps, m.range_enabled, (double)m.range_begin, (double)m.range_end, first, last);
}

static void movie_format_time(char* buf, size_t cap, double seconds) {
    const int s = (int)(seconds + 0.5);
    if (s >= 3600)    snprintf(buf, cap, "%d h %02d min", s / 3600, (s / 60) % 60);
    else if (s >= 60) snprintf(buf, cap, "%d min %02d s", s / 60, s % 60);
    else              snprintf(buf, cap, "%d s", s);
}

static bool movie_time_left(const ApplicationState* state, double* seconds) {
    const auto& m = state->movie;
    return movie_time_left(m.frame_index - m.rec_first, m.rec_last - m.rec_first + 1, m.rec_active_s, seconds);
}

// A time on the movie's timeline, moved to a frame when snapping is on
static double movie_snap_time(const ApplicationState* state, double time) {
    const auto& m = state->movie;
    if (m.snap_frames) return movie_snap_to_frame(time, (double)m.fps, (double)m.duration);
    return CLAMP(time, 0.0, (double)m.duration);
}

// Trajectory frame shown at a time on the movie timeline: from the keyframes that have one, otherwise
// the trajectory plays linearly between the times of the movie's timeline settings
static double movie_trajectory_frame(const ApplicationState* state, double time) {
    const auto& m = state->movie;
    const double f = camera_keyframes_frame_with_anchors(m.keyframes, md_array_size(m.keyframes),
        (double)m.traj_begin, m.start_frame, (double)m.traj_end, m.end_frame, time);
    const double last = (double)(run_num_frames(state) > 0 ? run_num_frames(state) - 1 : 0);
    return CLAMP(f, 0.0, last);
}

static void movie_sort_keyframes(ApplicationState* state) {
    std::stable_sort(state->movie.keyframes, state->movie.keyframes + md_array_size(state->movie.keyframes),
        [](const CameraKeyframe& a, const CameraKeyframe& b) { return a.time < b.time; });
}

// The look parameters that can be keyed. The id is what a key stores and what a workspace saves, so an
// id never changes meaning: new parameters go at the end.
enum MovieParamId : int {
    MovieParam_BackgroundColor = 0,
    MovieParam_BackgroundIntensity,
    MovieParam_SsaoIntensity,
    MovieParam_SsaoRadius, // Reserved for legacy scale-dependent SSAO.
    MovieParam_Exposure,
    MovieParam_DofStrength, // Reserved for legacy scale-dependent DOF.
    MovieParam_NearClip,
    MovieParam_FarClip,
    MovieParam_FocusDistance,
    MovieParam_DofAperture,
};
static_assert(MovieParam_DofAperture == MOVIE_PARAM_DOF_APERTURE, "the id of the blur is saved in workspaces");

struct MovieParamDesc {
    int         id;
    const char* label;
    int         comps;       // 3 for a colour
    bool        color;
    float       lo, hi;
    bool        log;         // Shown on a logarithmic axis, where the range spans orders of magnitude
    float*    (*ptr)(ApplicationState*);
    const char* tip;
};

static const MovieParamDesc movie_param_table[] = {
    { MovieParam_BackgroundColor,     "Background color",     3, true,  0.0f, 1.0f,     false, [](ApplicationState* s) { return s->visuals.background.color.elem; }, nullptr },
    { MovieParam_BackgroundIntensity, "Background intensity", 1, false, 0.0f, 100.0f,   false, [](ApplicationState* s) { return &s->visuals.background.intensity; }, nullptr },
    { MovieParam_SsaoIntensity,       "Ambient occlusion",    1, false, 0.0f, 10.0f,    false, [](ApplicationState* s) { return &s->visuals.ssao.intensity; }, "The intensity of the ambient occlusion. It has to be enabled." },
    { MovieParam_Exposure,            "Exposure",             1, false, 0.05f, 10.0f,   true,  [](ApplicationState* s) { return &s->visuals.tonemapping.exposure; }, "The exposure of the tonemapping. It has to be enabled." },
    { MovieParam_DofAperture,         "Depth of field blur",   1, false, 0.0f, 4.0f, false, [](ApplicationState* s) { return &s->visuals.dof.aperture; }, "Blur of distant objects, in percent of the view height. Depth of field has to be enabled." },
    { MovieParam_NearClip,            "Near clipping plane",  1, false, 0.01f, 5000.0f, true,  [](ApplicationState* s) { return &s->view.camera.near_plane; }, "Distance from the camera to where things start to show. Raise it to cut into the structure." },
    { MovieParam_FarClip,             "Far clipping plane",   1, false, 1.0f, 100000.0f,true,  [](ApplicationState* s) { return &s->view.camera.far_plane; }, "Distance from the camera to where things stop showing." },
};
static_assert(sizeof(movie_param_table) / sizeof(movie_param_table[0]) <= MOVIE_MAX_PARAMS, "more parameters than MOVIE_MAX_PARAMS");

static const MovieParamDesc* movie_param_desc(int id) {
    for (const MovieParamDesc& d : movie_param_table) {
        if (d.id == id) return &d;
    }
    return nullptr;
}

static const char* key_ease_str[(int)KeyEase::Count] = {
    "Smooth",
    "Ease in/out",
    "Linear",
    "Hold",
};

// Puts the keyed parameters at their values at a time. A parameter is taken hold of when it first has a key
// to follow, and what it was is kept, to be put back when the keys let go of it or the recording is over.
static void movie_params_apply(ApplicationState* state, double time) {
    auto& m = state->movie;
    for (const MovieParamDesc& d : movie_param_table) {
        float v[3] = {};
        const bool keyed = m.animate_params && param_keys_evaluate(v, d.comps, m.param_keys.data(), m.param_keys.size(), d.id, time);
        float* dst = d.ptr(state);
        if (keyed) {
            if (!m.param_saved_valid[d.id]) {
                for (int c = 0; c < d.comps; ++c) m.param_saved[d.id][c] = dst[c];
                m.param_saved_valid[d.id] = true;
            }
            for (int c = 0; c < d.comps; ++c) dst[c] = v[c];
        } else if (m.param_saved_valid[d.id]) {
            for (int c = 0; c < d.comps; ++c) dst[c] = m.param_saved[d.id][c];
            m.param_saved_valid[d.id] = false;
        }
    }
}

static void movie_params_restore(ApplicationState* state) {
    auto& m = state->movie;
    for (const MovieParamDesc& d : movie_param_table) {
        if (!m.param_saved_valid[d.id]) continue;
        float* dst = d.ptr(state);
        for (int c = 0; c < d.comps; ++c) dst[c] = m.param_saved[d.id][c];
        m.param_saved_valid[d.id] = false;
    }
}

static Representation* movie_find_rep(ApplicationState* state, uint32_t id) {
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        if (state->representation.reps[i].id == id) return &state->representation.reps[i];
    }
    return nullptr;
}

// What a property is called for a representation, or null when its type has no such thing. Also the range it is edited in.
static const char* movie_rep_prop_label(const Representation& rep, int prop, float* lo, float* hi) {
    *lo = 0.0f;
    *hi = 1.0f;
    switch ((RepProp)prop) {
    case RepProp::Visible:    return "Visible";
    case RepProp::TintScale:  return "Tint scale";
    case RepProp::Saturation: return "Saturation";
    case RepProp::BaseColor:  return "Base color";
    case RepProp::TintColor:  return "Tint color";
    case RepProp::Scale0:
    case RepProp::Scale1:
    case RepProp::Scale2: {
        const int c = prop - (int)RepProp::Scale0;
        *lo = 0.1f;
        *hi = 4.0f;
        switch (rep.type) {
        case RepresentationType::SpaceFill:
        case RepresentationType::Licorice:    return c == 0 ? "Radius scale" : nullptr;
        case RepresentationType::BallAndStick: return c == 0 ? "Ball scale" : c == 1 ? "Bond scale" : nullptr;
        case RepresentationType::Ribbons:      *hi = 3.0f; return c == 0 ? "Width" : c == 1 ? "Thickness" : nullptr;
        case RepresentationType::Cartoon:      *hi = 3.0f; return c == 0 ? "Coil" : c == 1 ? "Sheet" : "Helix";
        default: return nullptr;
        }
    }
    default: return nullptr;
    }
}

// out has rep_prop_comps(prop) numbers
static void movie_rep_prop_get(const Representation& rep, int prop, float* out) {
    switch ((RepProp)prop) {
    case RepProp::Visible:    out[0] = rep.enabled ? rep.presence : 0.0f; break;
    case RepProp::Scale0:     out[0] = rep.scale.elem[0]; break;
    case RepProp::Scale1:     out[0] = rep.scale.elem[1]; break;
    case RepProp::Scale2:     out[0] = rep.scale.elem[2]; break;
    case RepProp::TintScale:  out[0] = rep.tint_scale; break;
    case RepProp::Saturation: out[0] = rep.saturation; break;
    case RepProp::BaseColor:  for (int c = 0; c < 3; ++c) out[c] = rep.base_color.elem[c]; break;
    case RepProp::TintColor:  for (int c = 0; c < 3; ++c) out[c] = rep.tint_color.elem[c]; break;
    default: out[0] = 0.0f; break;
    }
}

// The scales are used as they are when drawing; what changes the colours of the atoms has to be updated
static void movie_rep_prop_set(ApplicationState* state, Representation* rep, int prop, const float* v) {
    float now[3] = {};
    movie_rep_prop_get(*rep, prop, now);
    bool same = true;
    for (int c = 0; c < rep_prop_comps(prop); ++c) same &= fabsf(now[c] - v[c]) < 1.0e-6f;
    if (same) return;
    switch ((RepProp)prop) {
    case RepProp::Visible: {
        // A value between 0 and 1 is a representation on its way in or out: shown, at that part of its size
        const bool was = rep->enabled;
        rep->enabled = v[0] > 0.002f;
        rep->presence = rep->enabled ? CLAMP(v[0], 0.0f, 1.0f) : 1.0f;
        if (was != rep->enabled) {
            state->representation.atom_visibility_mask_dirty = true;
            if (rep->enabled) flag_representation_as_dirty(rep);
        }
        break;
    }
    case RepProp::Scale0: rep->scale.elem[0] = v[0]; break;
    case RepProp::Scale1: rep->scale.elem[1] = v[0]; break;
    case RepProp::Scale2: rep->scale.elem[2] = v[0]; break;
    case RepProp::TintScale:  rep->tint_scale = v[0]; flag_representation_as_dirty(rep); break;
    case RepProp::Saturation: rep->saturation = v[0]; flag_representation_as_dirty(rep); break;
    case RepProp::BaseColor:
        for (int c = 0; c < 3; ++c) rep->base_color.elem[c] = v[c];
        if (rep->color_mapping == ColorMapping::Uniform) flag_representation_as_dirty(rep);
        break;
    case RepProp::TintColor:
        for (int c = 0; c < 3; ++c) rep->tint_color.elem[c] = v[c];
        if (rep->tint_scale > 0.0f) flag_representation_as_dirty(rep);
        break;
    default: break;
    }
}

static bool movie_rep_has_keys(const ApplicationState* state, uint32_t rep, int prop) {
    for (const RepKey& k : state->movie.rep_keys) {
        if (k.rep == rep && k.prop == prop) return true;
    }
    return false;
}

// Visible goes through its transition, the other properties are what their keys say
static bool movie_rep_eval(const ApplicationState* state, uint32_t rep, int prop, double time, float* out) {
    const auto& m = state->movie;
    if (prop == (int)RepProp::Visible) return rep_visible_factor(out, m.rep_keys.data(), m.rep_keys.size(), rep, time, (double)m.rep_transition);
    return rep_keys_evaluate(out, m.rep_keys.data(), m.rep_keys.size(), rep, prop, time);
}

// Puts the keyed properties of representations at their values at a time. A property is taken hold of when it first
// has a key to follow, and what it was is kept, to be put back when the keys let go of it or the recording is over.
static void movie_reps_apply(ApplicationState* state, double time) {
    auto& m = state->movie;
    if (m.animate_params) {
        for (size_t i = 0; i < m.rep_keys.size(); ++i) {
            const RepKey& k = m.rep_keys[i];
            bool seen = false;
            for (size_t j = 0; j < i && !seen; ++j) seen = m.rep_keys[j].rep == k.rep && m.rep_keys[j].prop == k.prop;
            Representation* rep = seen ? nullptr : movie_find_rep(state, k.rep);
            float v[3] = {};
            if (!rep || !movie_rep_eval(state, k.rep, k.prop, time, v)) continue;

            bool saved = false;
            for (const auto& s : m.rep_saved) saved |= s.rep == k.rep && s.prop == k.prop;
            if (!saved) {
                decltype(m.rep_saved)::value_type entry = {k.rep, k.prop, {}};
                movie_rep_prop_get(*rep, k.prop, entry.value);
                m.rep_saved.push_back(entry);
            }
            movie_rep_prop_set(state, rep, k.prop, v);
        }
    }
    for (size_t s = 0; s < m.rep_saved.size();) {
        const auto saved = m.rep_saved[s];
        Representation* rep = movie_find_rep(state, saved.rep);
        if (m.animate_params && rep && movie_rep_has_keys(state, saved.rep, saved.prop)) {
            ++s;
            continue;
        }
        if (rep) movie_rep_prop_set(state, rep, saved.prop, saved.value);
        m.rep_saved.erase(m.rep_saved.begin() + s);
    }
}

static void movie_reps_restore(ApplicationState* state) {
    auto& m = state->movie;
    for (const auto& s : m.rep_saved) {
        if (Representation* rep = movie_find_rep(state, s.rep)) movie_rep_prop_set(state, rep, s.prop, s.value);
    }
    m.rep_saved.clear();
}

// Where a selection is now: its centre, made whole across periodic boundaries, in the space the camera is in
static bool movie_target_center(const ApplicationState* state, const md_bitfield_t* mask, const md_system_state_t& system_state, vec3_t* out) {
    return movie_target_center(out, state->mold.sys, system_state, *mask, state->mold.unitcell_transform, state->allocator.frame);
}

// Where the movie's follow target is now (the global depth of field can follow it). False when there is none, or it
// holds atoms that are not in the system.
static bool movie_follow_center(const ApplicationState* state, vec3_t* out) {
    return movie_target_center(state, &state->movie.follow_mask, state->mold.state, out);
}

// What an expression selects in the system now, for a set or for the field that is typed. False with a message when it is not
// valid, or selects nothing.
static bool movie_expr_atoms(const ApplicationState* state, const std::string& expr, std::vector<uint32_t>* atoms, bool* dynamic, std::string* error) {
    atoms->clear();
    if (dynamic) *dynamic = false;
    const size_t num_atoms = state->mold.sys.atom.count;
    if (num_atoms == 0 || state->mold.state.num_atoms != num_atoms) {
        if (error) *error = "No system";
        return false;
    }
    md_temp_scope_t temp = md_temp_begin_in(state->allocator.frame);
    defer { md_temp_end(temp); };
    md_bitfield_t mask = {};
    md_bitfield_init(&mask, md_temp_allocator(temp));
    char message[256] = "";
    bool is_dynamic = false;
    const bool ok = md_filter(&mask, str_from_cstr(expr.c_str()), &state->mold.sys, &state->mold.state, state->script.ir, &is_dynamic, message, sizeof(message));
    if (dynamic) *dynamic = is_dynamic;
    if (!ok) {
        if (error) *error = message[0] ? message : "Not a valid selection";
        return false;
    }
    md_bitfield_iter_t it = md_bitfield_iter_create(&mask);
    while (md_bitfield_iter_next(&it)) atoms->push_back((uint32_t)md_bitfield_iter_idx(&it));
    if (atoms->empty()) {
        if (error) *error = "Selects no atoms";
        return false;
    }
    return true;
}

// The atoms of a set now: its list, or what its expression selects (again when the script or the system changed, or the frame if
// the expression depends on it). Null when the set is not there or its expression does not select anything.
static const std::vector<uint32_t>* movie_set_atoms(const ApplicationState* state, uint32_t id) {
    const MovieAtomSet* set = movie_atoms_find(state->movie.atom_sets, id);
    if (!set) return nullptr;
    if (set->expr.empty()) return &set->atoms;
    auto& cache = state->movie.set_eval;
    MovieSetEval* e = nullptr;
    for (auto& c : cache) if (c.id == id) e = &c;
    if (!e) {
        cache.emplace_back();
        e = &cache.back();
        e->id = id;
    }
    const uint64_t fingerprint = state->script.ir ? state->script.ir_fingerprint : 0;
    const size_t num_atoms = state->mold.sys.atom.count;
    const double frame = state->animation.frame;
    if (e->expr != set->expr || e->ir_fingerprint != fingerprint || e->num_atoms != num_atoms || (e->dynamic && e->frame != frame)) {
        e->expr = set->expr;
        e->ir_fingerprint = fingerprint;
        e->num_atoms = num_atoms;
        e->frame = frame;
        e->valid = movie_expr_atoms(state, set->expr, &e->atoms, &e->dynamic, &e->error);
    }
    return e->valid ? &e->atoms : nullptr;
}

// Where a set of atoms of the movie is now
static bool movie_set_center(const ApplicationState* state, uint32_t id, vec3_t* out) {
    const std::vector<uint32_t>* atoms = movie_set_atoms(state, id);
    if (!atoms || atoms->empty()) return false;
    md_temp_scope_t temp = md_temp_begin_in(state->allocator.frame);
    defer { md_temp_end(temp); };
    md_bitfield_t mask = {};
    md_bitfield_init(&mask, md_temp_allocator(temp));
    for (uint32_t atom : *atoms) md_bitfield_set_bit(&mask, atom);
    return movie_target_center(state, &mask, state->mold.state, out);
}

// The atoms of the selection
static std::vector<uint32_t> movie_selection_atoms(const ApplicationState* state) {
    std::vector<uint32_t> atoms;
    md_bitfield_iter_t it = md_bitfield_iter_create(&state->selection.selection_mask);
    while (md_bitfield_iter_next(&it)) atoms.push_back((uint32_t)md_bitfield_iter_idx(&it));
    return atoms;
}

// Where an atom is now, in the space the camera is in
static bool movie_atom_position(const ApplicationState* state, int32_t atom, vec3_t* out) {
    if (atom < 0 || (size_t)atom >= state->mold.sys.atom.count || state->mold.state.num_atoms != state->mold.sys.atom.count) return false;
    *out = mat4_mul_vec3(state->mold.unitcell_transform, state->mold.state.xyz[atom], 1.0f);
    return true;
}

static const vec3_t movie_up_axes[6] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}};
static const char* movie_up_axis_str[6] = {"+X", "+Y", "+Z", "-X", "-Y", "-Z"};

static const vec3_t* movie_up_vector(const ApplicationState* state) {
    return &movie_up_axes[CLAMP(state->movie.up_axis, 0, 5)];
}

// The world up the movie keeps the camera level about, or null when it does not
static const vec3_t* movie_upright(const ApplicationState* state) {
    return state->movie.keep_upright ? movie_up_vector(state) : nullptr;
}

// The keys, sorted by time, with where each of them looks at and focuses on now: the centre of its set of atoms where it has one
// (the saved point where the atoms are not in this system, which is reported once), otherwise its fixed point.
struct MovieResolved {
    std::vector<CameraKeyframe> keys;
    std::vector<vec3_t> look, focus;
};

static MovieResolved movie_resolve(const ApplicationState* state, const CameraKeyframe* keys, size_t n) {
    MovieResolved r;
    r.keys.assign(keys, keys + n);
    std::stable_sort(r.keys.begin(), r.keys.end(), [](const CameraKeyframe& a, const CameraKeyframe& b) { return a.time < b.time; });
    auto missing = [&]() {
        if (state->movie.invalid_target_reported) return;
        VIAMD_LOG_ERROR("A movie key looks at or focuses on atoms that are not in this system: its saved point is used. Set them again with Look at / Focus on in the key table.");
        state->movie.invalid_target_reported = true;
    };
    for (const CameraKeyframe& k : r.keys) {
        vec3_t look = camera_get_look_at(k.transform), focus = k.focus_point;
        if (k.look_set != 0 && !movie_set_center(state, k.look_set, &look)) {
            look = camera_get_look_at(k.transform);
            missing();
        }
        if (k.focus_on && k.focus_target == FocusTarget::Selection && !movie_set_center(state, k.focus_set, &focus)) {
            focus = k.focus_point;
            missing();
        }
        r.look.push_back(look);
        r.focus.push_back(focus);
    }
    return r;
}

// Whether any key tracks atoms: then the camera can only be placed once the trajectory frame is loaded
static bool movie_keys_track_atoms(const CameraKeyframe* keys, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        if (keys[i].look_set != 0 || (keys[i].focus_on && keys[i].focus_target == FocusTarget::Selection)) return true;
    }
    return false;
}

static void movie_pose(const ApplicationState* state, const MovieResolved& r, double time, ViewTransform* vt, float* fov_y) {
    movie_keys_pose(vt, fov_y, r.keys.data(), r.keys.size(), time, r.look.data(), movie_upright(state), *movie_up_vector(state));
}

// The camera the keys give at a time. 'keys' need not be sorted, and are not empty.
static void movie_camera_evaluate(const ApplicationState* state, double time, const CameraKeyframe* keys, size_t num_keys, ViewTransform* vt, float* fov_y) {
    movie_pose(state, movie_resolve(state, keys, num_keys), time, vt, fov_y);
}

// What is sharp for a camera at a time, from the keys that set the focus: before the first of them what the camera looks at is
// sharp.
static vec3_t movie_focus_evaluate(const MovieResolved& r, double time, const ViewTransform& camera) {
    return movie_focus_point(r.keys.data(), r.keys.size(), r.focus.data(), time, camera);
}

static bool movie_focus_owned(const ApplicationState* state) {
    const auto& m = state->movie;
    return m.show_window || m.play_mode || m.state == MovieRecordingState::Recording;
}

// The scene editor has its own fixed depth; only movie-camera surfaces evaluate movie focus.
static bool movie_focus_active(const ApplicationState* state) {
    const auto& m = state->movie;
    return movie_focus_owned(state) && (m.show_frame || m.pip_pass || m.play_mode || m.state == MovieRecordingState::Recording);
}

static double movie_focus_time(const ApplicationState* state) {
    return state->movie.state == MovieRecordingState::Recording ? state->movie.cur_time : (double)state->movie.playhead;
}

// Aims a key, from where its eye is, at the centre of a set of atoms that it then tracks through the trajectory
static bool movie_key_look_at_set(ApplicationState* state, int key_idx, uint32_t id) {
    auto& m = state->movie;
    if (key_idx < 0 || key_idx >= (int)md_array_size(m.keyframes) || id == 0) return false;
    vec3_t center;
    if (!movie_set_center(state, id, &center)) {
        VIAMD_LOG_ERROR("Nothing to look at: the atoms are not in the current system");
        return false;
    }
    CameraKeyframe& key = m.keyframes[key_idx];
    if (!camera_key_set_look(&key, center)) return false;
    key.look_set = id;
    m.invalid_target_reported = false;
    return true;
}

static bool movie_key_look_at_atoms(ApplicationState* state, int key_idx, const std::vector<uint32_t>& atoms) {
    if (atoms.empty()) return false;
    return movie_key_look_at_set(state, key_idx, movie_atoms_add(&state->movie.atom_sets, atoms));
}

// A selection in the language of the filters (protein, resname("AIN")), evaluated again wherever the movie is shown
static bool movie_key_look_at_expr(ApplicationState* state, int key_idx, const std::string& expr) {
    return movie_key_look_at_set(state, key_idx, movie_atoms_add_expr(&state->movie.atom_sets, expr));
}

// Makes a keyframe look at an atom and track it through the trajectory. The eye stays where it is.
static bool movie_key_look_at_atom(ApplicationState* state, int key_idx, int32_t atom) {
    vec3_t pos;
    if (atom < 0 || !movie_atom_position(state, atom, &pos)) return false;
    return movie_key_look_at_atoms(state, key_idx, {(uint32_t)atom});
}

// Makes a keyframe look at a fixed point. The eye stays where it is.
static bool movie_key_look_at_point(ApplicationState* state, int key_idx, vec3_t point) {
    auto& m = state->movie;
    if (key_idx < 0 || key_idx >= (int)md_array_size(m.keyframes)) return false;
    CameraKeyframe& key = m.keyframes[key_idx];
    if (!camera_key_set_look(&key, point)) return false;
    key.look_set = 0;
    return true;
}

// A key that sets the focus starts with what the camera looks at
static void movie_key_focus_begin(ApplicationState* state, int key_idx) {
    CameraKeyframe& key = state->movie.keyframes[key_idx];
    if (key.focus_on) return;
    key.focus_on = true;
    key.focus_target = FocusTarget::LookAt;
    key.focus_blur = -1.0f;
}

static bool movie_key_focus_on_set(ApplicationState* state, int key_idx, uint32_t id) {
    auto& m = state->movie;
    if (key_idx < 0 || key_idx >= (int)md_array_size(m.keyframes) || id == 0) return false;
    vec3_t center;
    if (!movie_set_center(state, id, &center)) {
        VIAMD_LOG_ERROR("Nothing to focus on: the atoms are not in the current system");
        return false;
    }
    movie_key_focus_begin(state, key_idx);
    CameraKeyframe& key = m.keyframes[key_idx];
    key.focus_target = FocusTarget::Selection;
    key.focus_set = id;
    key.focus_point = center;
    m.invalid_target_reported = false;
    return true;
}

static bool movie_key_focus_on_atoms(ApplicationState* state, int key_idx, const std::vector<uint32_t>& atoms) {
    if (atoms.empty()) return false;
    return movie_key_focus_on_set(state, key_idx, movie_atoms_add(&state->movie.atom_sets, atoms));
}

static bool movie_key_focus_on_expr(ApplicationState* state, int key_idx, const std::string& expr) {
    return movie_key_focus_on_set(state, key_idx, movie_atoms_add_expr(&state->movie.atom_sets, expr));
}

// How far from the camera depth of field is sharp, for a camera: from the movie's focus keys when it has some, otherwise from the
// focus mode of the Depth of Field settings
static float dof_focus_depth(const ApplicationState* state, const ViewTransform& view) {
    if (movie_focus_active(state)) {
        const auto& m = state->movie;
        const MovieResolved r = movie_resolve(state, m.keyframes, md_array_size(m.keyframes));
        return MAX(camera_depth_of_point(view, movie_focus_evaluate(r, movie_focus_time(state), view)), 1.0e-3f);
    }
    if (movie_focus_owned(state)) return MAX(state->movie.editor_focus_distance, 1.0e-3f);
    return MAX(view.distance, 1.0e-3f);
}

static void movie_camera_apply(ApplicationState* state, double time, const CameraKeyframe* keys, size_t num_keys) {
    ViewTransform vt;
    float fov_y;
    movie_camera_evaluate(state, time, keys, num_keys, &vt, &fov_y);
    state->view.target = vt;
    state->view.camera = vt;
    state->view.camera.fov_y = fov_y;
}

// Scene view: the viewport is the editor's own camera, for working on the path from outside. The movie camera is then only drawn in
// it (the green camera), and the movie moves it only while recording.
static bool movie_scene_view(const ApplicationState* state) {
    const auto& m = state->movie;
    return m.show_window && !m.show_frame && m.state != MovieRecordingState::Recording;
}

// Both view targets are set so the exponential smoothing in camera_animate does not lag behind.
// The keys must be sorted by time.
static void movie_apply_time_with_keys(ApplicationState* state, double time, bool apply_camera, const CameraKeyframe* keys, size_t num_keys) {
    auto& m = state->movie;
    state->animation.frame = movie_trajectory_frame(state, time);
    movie_params_apply(state, time);
    movie_reps_apply(state, time);
    m.follow_pending = false;
    if (apply_camera && num_keys > 0 && !movie_scene_view(state)) {
        if (movie_keys_track_atoms(keys, num_keys)) {
            // The atoms are not where they will be until the trajectory frame has been loaded, so the camera waits for that
            m.follow_keys.assign(keys, keys + num_keys);
            m.follow_time = time;
            m.follow_pending = true;
        } else {
            movie_camera_apply(state, time, keys, num_keys);
        }
    }
}

// Runs once the frame of the movie's time is in the system state: puts the camera where the keys say, around the atoms they track
static void movie_follow_update(ApplicationState* state) {
    auto& m = state->movie;
    if (!m.follow_pending) return;
    m.follow_pending = false;
    movie_camera_apply(state, m.follow_time, m.follow_keys.data(), m.follow_keys.size());
}

// Shows the movie at a time on its timeline: the trajectory frame and, if enabled, the camera.
static void movie_apply_time(ApplicationState* state, double time, bool apply_camera) {
    movie_apply_time_with_keys(state, time, apply_camera, state->movie.keyframes, md_array_size(state->movie.keyframes));
}

static void movie_restore_state(ApplicationState* state) {
    state->animation.mode = state->movie.prev_playback_mode;
    state->movie.follow_pending = false;
    movie_params_restore(state);
    movie_reps_restore(state);
    if (state->movie.camera_was_animated) {
        state->view.target = state->movie.prev_view_target;
        state->view.camera = state->movie.prev_view_target;
        state->view.camera.fov_y = state->movie.prev_fov_y;
        state->movie.camera_was_animated = false;
    }
}

static void movie_frame_size(const ApplicationState* state, int* w, int* h) {
    switch (state->movie.resolution) {
    case ScreenshotResolution::Window:
        *w = (int)state->app.framebuffer.width;
        *h = (int)state->app.framebuffer.height;
        break;
    case ScreenshotResolution::FHD:    *w = 1920; *h = 1080; break;
    case ScreenshotResolution::QHD:    *w = 2560; *h = 1440; break;
    case ScreenshotResolution::UHD_4K: *w = 3840; *h = 2160; break;
    case ScreenshotResolution::UHD_8K: *w = 7680; *h = 4320; break;
    default:
        *w = CLAMP(state->movie.res_x, 640, 16384);
        *h = CLAMP(state->movie.res_y, 480, 16384);
        break;
    }
    movie_scaled_size(w, h, state->movie.res_scale);
}

// The part of the viewport that the Movie window does not cover: the largest of the strips left, right, above and below it, judged
// by how large a picture of proportions 'aspect' fits in. The whole viewport when the window is not over it, in the Preview, or when
// it floats in the middle so that no strip beside it is worth it.
static void movie_viewport_free_rect(const ApplicationState* state, float aspect, float* x, float* y, float* w, float* h) {
    const auto& m = state->movie;
    const float vw = (float)state->app.window.width, vh = (float)state->app.window.height;
    *x = 0.0f; *y = 0.0f; *w = vw; *h = vh;
    if (!m.show_window || m.play_mode || ImGui::GetFrameCount() - m.window_rect_frame > 2 || vw <= 0.0f || vh <= 0.0f) return;
    const float x0 = MAX(m.window_rect[0], 0.0f), y0 = MAX(m.window_rect[1], 0.0f);
    const float x1 = MIN(m.window_rect[2], vw), y1 = MIN(m.window_rect[3], vh);
    if (x1 <= x0 || y1 <= y0) return;
    auto fit_h = [aspect](float rw, float rh) { return rw > 0.0f && rh > 0.0f ? MIN(rh, rw / MAX(aspect, 1e-3f)) : 0.0f; };
    const float cand[4][4] = {{0.0f, 0.0f, x0, vh}, {x1, 0.0f, vw - x1, vh}, {0.0f, 0.0f, vw, y0}, {0.0f, y1, vw, vh - y1}};
    int best = -1;
    float best_h = 0.30f * fit_h(vw, vh);
    for (int i = 0; i < 4; ++i) {
        const float fh = fit_h(cand[i][2], cand[i][3]);
        if (fh > best_h) { best_h = fh; best = i; }
    }
    if (best < 0) return;
    *x = cand[best][0]; *y = cand[best][1]; *w = cand[best][2]; *h = cand[best][3];
}

// Where the frame of the movie is in the viewport, fitted with a little room into the part the Movie window does not cover: the
// preview shows the movie as it will be recorded, in the proportions of the frame. Scene view shows the same frame, so that what is
// in it is what a key added from the view records, in both modes. False while there is nothing to show (not
// editing the movie, recording, a screenshot).
static bool movie_frame_guide(const ApplicationState* state, ImVec2* pos, ImVec2* size) {
    const auto& m = state->movie;
    if (m.state == MovieRecordingState::Recording || !m.show_window) return false;
    if (!str_empty(state->screenshot.path_to_file)) return false;
    int fw = 0, fh = 0;
    movie_frame_size(state, &fw, &fh);
    const float vw = (float)state->app.window.width, vh = (float)state->app.window.height;
    if (fw <= 0 || fh <= 0 || vw <= 0.0f || vh <= 0.0f) return false;
    float rx = 0.0f, ry = 0.0f, rw = vw, rh = vh;
    if (!m.pip_pass) movie_viewport_free_rect(state, (float)fw / (float)fh, &rx, &ry, &rw, &rh);
    movie_frame_fit(rw, rh, (float)fw, (float)fh, 0.92f, &pos->x, &pos->y, &size->x, &size->y);
    pos->x += rx;
    pos->y += ry;
    return true;
}

// How far the frame's centre is from the viewport's, in clip space (x right, y up). The projection is shifted by as much so that the
// frame shows the movie camera's view centred, as in the recording, wherever the frame is.
static bool movie_frame_guide_shift(const ApplicationState* state, float* sx, float* sy) {
    *sx = *sy = 0.0f;
    ImVec2 gp, gs;
    const float vw = (float)state->app.window.width, vh = (float)state->app.window.height;
    if (vw <= 0.0f || vh <= 0.0f || !movie_frame_guide(state, &gp, &gs)) return false;
    *sx = 2.0f * (gp.x + gs.x * 0.5f - vw * 0.5f) / vw;
    *sy = -2.0f * (gp.y + gs.y * 0.5f - vh * 0.5f) / vh;
    return *sx != 0.0f || *sy != 0.0f;
}

static void movie_pbo_free(ApplicationState* state) {
    auto& m = state->movie;
    for (int i = 0; i < MOVIE_RING_SIZE; ++i) {
        if (m.fence[i]) {
            glDeleteSync(m.fence[i]);
            m.fence[i] = nullptr;
        }
    }
    if (m.pbo[0]) {
        glDeleteBuffers(MOVIE_RING_SIZE, m.pbo);
        for (int i = 0; i < MOVIE_RING_SIZE; ++i) m.pbo[i] = 0;
    }
    m.pbo_head = 0;
    m.pbo_count = 0;
    m.pbo_bytes = 0;
}

static void movie_pbo_alloc(ApplicationState* state, size_t bytes) {
    movie_pbo_free(state);
    auto& m = state->movie;
    glGenBuffers(MOVIE_RING_SIZE, m.pbo);
    for (int i = 0; i < MOVIE_RING_SIZE; ++i) {
        glBindBuffer(GL_PIXEL_PACK_BUFFER, m.pbo[i]);
        glBufferData(GL_PIXEL_PACK_BUFFER, (GLsizeiptr)bytes, nullptr, GL_STREAM_READ);
    }
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    m.pbo_bytes = bytes;
}

// Hands the oldest read back frame to the sink. Without wait it only does so if the GPU is done with
// it and the sink has a buffer to spare, otherwise it leaves everything as it was and returns false.
static bool movie_pbo_pop(ApplicationState* state, bool wait) {
    auto& m = state->movie;
    if (m.pbo_count == 0 || !m.sink) return false;

    const int slot = m.pbo_head;
    if (!wait && m.fence[slot]) {
        const GLenum r = glClientWaitSync(m.fence[slot], 0, 0);
        if (r != GL_ALREADY_SIGNALED && r != GL_CONDITION_SATISFIED) return false;
    }

    uint8_t* dst = frame_sink::acquire(m.sink, wait);
    if (!dst) return false;

    glBindBuffer(GL_PIXEL_PACK_BUFFER, m.pbo[slot]);
    const void* src = glMapBufferRange(GL_PIXEL_PACK_BUFFER, 0, (GLsizeiptr)m.pbo_bytes, GL_MAP_READ_BIT);
    if (src) {
        memcpy(dst, src, m.pbo_bytes);
        glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
        frame_sink::submit(m.sink, dst, m.pbo_frame[slot]);
    } else {
        VIAMD_LOG_ERROR("Could not read back movie frame %d from the GPU", m.pbo_frame[slot]);
        frame_sink::release(m.sink, dst);
    }
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    if (m.fence[slot]) {
        glDeleteSync(m.fence[slot]);
        m.fence[slot] = nullptr;
    }
    m.pbo_head = (m.pbo_head + 1) % MOVIE_RING_SIZE;
    m.pbo_count -= 1;
    return true;
}

static void movie_pbo_flush(ApplicationState* state) {
    while (movie_pbo_pop(state, true)) {}
}

static const char* movie_overlay_type_str[(int)MovieOverlayType::Count] = {
    "Text",
    "Time stamp",
    "Scale bar",
    "Logo",
    "Image",
    "Time bar",
    "Timeline",
    "Distribution",
    "Property visualization",
    "Figure (old)",
};

static const char* movie_overlay_anchor_str[(int)MovieOverlayAnchor::Count] = {
    "Top left", "Top center", "Top right",
    "Middle left", "Center", "Middle right",
    "Bottom left", "Bottom center", "Bottom right",
};

// Makes a texture of rgba pixels, with mipmaps so that a small picture of a large image is smooth
static GLuint movie_upload_texture(const uint8_t* pixels, int w, int h) {
    GLint prev = 0;
    GLuint tex = 0;
    glActiveTexture(GL_TEXTURE0);
    glGetIntegerv(GL_TEXTURE_BINDING_2D, &prev);
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
    glGenerateMipmap(GL_TEXTURE_2D);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, (GLuint)prev);
    return tex;
}

// The logo as a texture, made the first time it is needed. 'aspect' is its width over its height. 0 if it cannot be made.
static GLuint movie_logo_texture(float* aspect) {
    static GLuint tex = 0;
    static float  logo_aspect = 4.0f;
    static bool   tried = false;
    if (!tried) {
        tried = true;
        int w = 0, h = 0;
        uint8_t* pixels = image_decode_rgba(viamd_png, viamd_png_size, &w, &h);
        if (pixels) {
            // The icon is a small drawing in the middle of a transparent square: the logo is only that drawing, with a little room
            int x0 = 0, y0 = 0, x1 = w, y1 = h;
            image_alpha_bounds(pixels, w, h, 8, &x0, &y0, &x1, &y1);
            const int pad = MAX((x1 - x0), (y1 - y0)) / 50;
            x0 = MAX(x0 - pad, 0); y0 = MAX(y0 - pad, 0); x1 = MIN(x1 + pad, w); y1 = MIN(y1 + pad, h);
            const int cw = x1 - x0, ch = y1 - y0;
            std::vector<uint8_t> cropped((size_t)cw * (size_t)ch * 4);
            for (int y = 0; y < ch; ++y) memcpy(&cropped[(size_t)y * (size_t)cw * 4], &pixels[((size_t)(y0 + y) * (size_t)w + (size_t)x0) * 4], (size_t)cw * 4);
            // The lettering of the icon is light grey, for a dark picture: black reads on a light one
            image_replace_light_grey(cropped.data(), cw, ch, 128, 24, 0, 0, 0);
            image_bleed_transparent(cropped.data(), cw, ch);
            tex = movie_upload_texture(cropped.data(), cw, ch);
            logo_aspect = (float)cw / (float)MAX(ch, 1);
            image_free(pixels);
        } else {
            VIAMD_LOG_ERROR("Could not read the logo for the movie overlays");
        }
    }
    *aspect = logo_aspect;
    return tex;
}

// The images of image overlays, read once per file. A file that cannot be read is remembered as such, so it is not tried every frame.
struct MovieImageTexture {
    std::string path;
    GLuint tex = 0;
    float  aspect = 1.0f;
};
static std::vector<MovieImageTexture> movie_image_cache;

static GLuint movie_image_texture(const char* path, float* aspect) {
    if (!path || !path[0]) return 0;
    for (const MovieImageTexture& c : movie_image_cache) {
        if (c.path == path) {
            *aspect = c.aspect;
            return c.tex;
        }
    }
    MovieImageTexture entry;
    entry.path = path;
    if (FILE* f = fopen(path, "rb")) {
        std::vector<uint8_t> bytes;
        uint8_t buf[65536];
        size_t n;
        while ((n = fread(buf, 1, sizeof(buf), f)) > 0) bytes.insert(bytes.end(), buf, buf + n);
        fclose(f);
        int w = 0, h = 0;
        uint8_t* pixels = image_decode_rgba(bytes.data(), bytes.size(), &w, &h);
        if (pixels) {
            image_bleed_transparent(pixels, w, h);
            entry.tex = movie_upload_texture(pixels, w, h);
            entry.aspect = (float)w / (float)MAX(h, 1);
            image_free(pixels);
        } else {
            VIAMD_LOG_ERROR("Movie overlay: '%s' is not a png or jpg image", path);
        }
    } else {
        VIAMD_LOG_ERROR("Movie overlay: could not open '%s'", path);
    }
    movie_image_cache.push_back(entry);
    *aspect = entry.aspect;
    return entry.tex;
}

// Reads the file again the next time it is needed
static void movie_image_forget(const char* path) {
    for (size_t i = 0; i < movie_image_cache.size(); ++i) {
        if (movie_image_cache[i].path == path) {
            if (movie_image_cache[i].tex) glDeleteTextures(1, &movie_image_cache[i].tex);
            movie_image_cache.erase(movie_image_cache.begin() + i);
            return;
        }
    }
}

// The trajectory time (or the frame, without times) at a movie time
static double movie_trajectory_quantity(const ApplicationState* state, double time) {
    const double frame = movie_trajectory_frame(state, time);
    return md_array_size(state->timeline.x_values) > 0 ? frame_to_time(frame, *state) : frame;
}

// How far the movie has taken the trajectory by each time, made again when what it depends on changes
static const MovieTimeBarProfile& movie_time_bar_profile_for(const ApplicationState* state) {
    static MovieTimeBarProfile profile;
    static uint64_t signature = 0;
    const auto& m = state->movie;
    uint64_t h = 7;
    for (size_t i = 0; i < md_array_size(m.keyframes); ++i) {
        const double v[] = {m.keyframes[i].time, m.keyframes[i].use_frame ? m.keyframes[i].frame : -1.0, (double)(int)m.keyframes[i].ease};
        h = md_hash64_combine(h, md_hash64(v, sizeof(v), 5));
    }
    const size_t nx = md_array_size(state->timeline.x_values);
    const double s[] = {(double)m.duration, (double)m.traj_begin, (double)m.traj_end, m.start_frame, m.end_frame, (double)run_num_frames(state), (double)nx,
        nx ? (double)state->timeline.x_values[0] : 0.0, nx ? (double)state->timeline.x_values[nx - 1] : 0.0};
    h = md_hash64_combine(h, md_hash64(s, sizeof(s), 6));
    if (h != signature || profile.distance.empty()) {
        signature = h;
        movie_time_bar_profile(&profile, (double)m.duration, 600, [state](double t) { return movie_trajectory_quantity(state, t); });
    }
    return profile;
}

// The labels of the script's visualization (distances, angles, ...), at the places they have on the screen. 'res' is the
// size of what is drawn into and 'scale' how much larger than the viewport's text the text is made, for a frame of another size.
static void script_vis_text_draw(ImDrawList* dl, ImVec2 res, float scale, const ApplicationState& state) {
    if (!state.script.vis.text) return;
    ImFont* font = ImGui::GetFont();
    const float font_px = ImGui::GetFontSize() * scale;
    const mat4_t mvp = state.view.param.matrix.curr.proj_no_jitter * state.view.param.matrix.curr.view;

    // The visualization of a movie overlay fades in and out with it
    vec4_t text_col = state.script.text_color, bg_col = state.script.text_bg_color;
    text_col.w *= state.movie.vis_fade;
    bg_col.w *= state.movie.vis_fade;
    const ImU32 text_color = convert_color(text_col);
    const ImU32 rect_color = convert_color(bg_col);
    const float rect_rounding = 5.f * scale;
    const ImVec2 rect_padding = ImVec2(4.f, 2.f) * scale;

    const size_t num_text = md_array_size(state.script.vis.text);
    for (size_t i = 0; i < num_text; ++i) {
        const md_script_vis_text_t& vis_text = state.script.vis.text[i];

        const vec4_t p = mat4_mul_vec4(mvp, vec4_from_vec3(vis_text.pos, 1.0f));
        const vec4_t c = p / p.w;

        // A label carrying a unit is a quantity, so it is shown in the same unit as the
        // plots and re-formatted here rather than taken as mdlib spelled it. Anything
        // else is drawn as it came.
        char buf[64];
        str_t str = vis_text.str;
        if (!md_unit_is_none(vis_text.unit)) {
            char unit_buf[32];
            const double scl = display_units::factor_print(unit_buf, sizeof(unit_buf), vis_text.unit);
            // The degree sign hugs the number the way the convention has it, everything else takes a space.
            const char* sep = strcmp(unit_buf, "\xC2\xB0") == 0 ? "" : " ";
            const int len = snprintf(buf, sizeof(buf), "%.2f%s%s", vis_text.value * scl, sep, unit_buf);
            str = {buf, (size_t)CLAMP(len, 0, (int)sizeof(buf) - 1)};
        }

        const ImVec2 text_size = font->CalcTextSizeA(font_px, FLT_MAX, 0.0f, str.beg(), str.end());

        if (-1 < c.x && c.x < 1 && -1 < c.y && c.y < 1 && -1 < c.z && c.z < 1) {
            ImVec2 tc = {(c.x * 0.5f + 0.5f) * res.x, (-c.y * 0.5f + 0.5f) * res.y};
            ImVec2 p0 = tc - text_size * 0.5f;
            ImVec2 p1 = tc + text_size * 0.5f;
            dl->AddRectFilled(p0 - rect_padding, p1 + rect_padding, rect_color, rect_rounding);
            dl->AddText(font, font_px, p0, text_color, str.beg(), str.end());
        }
    }
}

// Evaluates the visualization of the script property of each property overlay that is shown, in the viewport while the
// movie is previewed or recorded. What it makes (highlight, geometry, labels) is drawn as if the property were hovered.
static void movie_property_vis_apply(ApplicationState* state) {
    auto& m = state->movie;
    float fade = 0.0f;
    bool highlighted = false;
    const bool recording = m.state == MovieRecordingState::Recording;
    if (!m.overlays.empty() && state->script.eval_ir && (recording || (m.show_overlay_preview && m.show_window))) {
        const double time = recording ? m.cur_time : (double)m.playhead;
        for (const MovieOverlay& o : m.overlays) {
            const float alpha = movie_overlay_alpha(o, time);
            if (o.type != MovieOverlayType::PropertyVis || o.text[0] == '\0' || alpha <= 0.0f) continue;
            const md_script_vis_payload_o* payload = md_script_ir_property_vis_payload(state->script.eval_ir, str_from_cstr(o.text));
            if (payload) {
                script_visualize_payload(state, payload, -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                highlighted |= !md_bitfield_empty(&state->script.vis.atom_mask);
                fade = MAX(fade, alpha);
            }
        }
    }
    if (fade > 0.0f) m.vis_fade = fade;

    // The highlight is state that stays until something clears it: when the overlay that set it is not shown any more it goes too
    if (!highlighted && m.vis_highlight) md_bitfield_clear(&state->selection.highlight_mask);
    m.vis_highlight = highlighted;
}

// The value of a temporal series at a position on its axis, in the display unit
static double movie_series_value_at(const SeriesTemporalView& v, double x) {
    if (v.num_samples <= 0) return 0.0;
    const double idx = series_temporal_index_at(v, x);
    const int i0 = CLAMP((int)floor(idx), 0, v.num_samples - 1);
    const int i1 = MIN(i0 + 1, v.num_samples - 1);
    const double f = CLAMP(idx - (double)i0, 0.0, 1.0);
    const double a = (double)v.y[(size_t)i0 * v.stride] * v.y_scale;
    const double b = (double)v.y[(size_t)i1 * v.stride] * v.y_scale;
    return a + (b - a) * f;
}

// Decimals that a tick step needs
static int movie_tick_decimals(double step) {
    return step > 0.0 ? CLAMP((int)ceil(-log10(step) - 1e-9), 0, 6) : 0;
}

// What the panels of a plot overlay share
struct MoviePlotContext {
    ImDrawList* dl;
    ImFont* font;
    const MovieOverlay* o;
    ApplicationState* state;
    const MovieTimeBarProfile* prof;
    float fpx;
    float alpha;
    double time;
    double q;                        // The trajectory time (frame) that is shown
    double axis_begin, axis_end;     // The movie times that the horizontal axis covers: the whole movie, whenever the overlay is shown
    float thick;                     // Width of the lines
    int color_index = 0;             // How many series have been drawn, for the colours of a palette
};

// One subplot as a panel of 'w' by 'h' at 'p0': its legend with the value at the frame that is shown, the axes and the
// curves (or the bars). The labels of the horizontal axis are left out of a timeline panel that has another below it.
static void movie_plot_panel_draw(MoviePlotContext& c, const PlotSubplot& sp, const char* title_text, MoviePlotView view, ImVec2 p0, float w, float h, bool x_labels, bool marker_labels) {
    ImDrawList* dl = c.dl;
    ImFont* font = c.font;
    const MovieOverlay& o = *c.o;
    ApplicationState* state = c.state;
    const MovieTimeBarProfile& prof = *c.prof;
    const float fpx = c.fpx;
    const bool timeline = view == MoviePlotView::Timeline;
    const bool elapsed = timeline && o.plot_axis == MoviePlotAxis::Elapsed;
    const int ns = MIN(sp.count, 6);

    // Whenever a subplot comes in, it is drawn from the start of the movie, with what happened before it came in
    const double pt0 = c.axis_begin, pt1 = c.axis_end;
    const double tnow = CLAMP(c.time, pt0, pt1);
    double full_lo = 0.0, full_hi = 0.0, vis_lo = 0.0, vis_hi = 0.0;
    movie_time_bar_visited_between(prof, pt0, pt1, &full_lo, &full_hi);
    movie_time_bar_visited_between(prof, pt0, o.reveal ? tnow : pt1, &vis_lo, &vis_hi);
    const double d_axis0 = movie_time_bar_moved(prof, c.axis_begin);

    const int color_base = c.color_index;
    c.color_index += ns;
    auto series_color = [&](int s) -> ImVec4 {
        return o.palette > 0 ? ImPlot::GetColormapColor(color_base + s, o.palette - 1) : sp.series[s].color;
    };

    const char* title = title_text[0] != '\0' ? title_text : sp.name;
    const bool titled = o.show_titles && title[0] != '\0';
    const float title_h  = titled ? fpx * 1.35f : 0.0f;
    const float legend_h = title_h + (float)ns * fpx * 1.2f;
    const float axis_h   = x_labels ? fpx * 2.6f : fpx * 0.6f;
    const float label_w  = timeline ? fpx * 3.8f : fpx * 0.6f;
    const float pad      = fpx * 0.4f;
    const float ix0 = p0.x + label_w, ix1 = p0.x + w - pad;
    const float iy0 = p0.y + legend_h, iy1 = p0.y + h - axis_h;
    if (ix1 - ix0 < 8.0f || iy1 - iy0 < 8.0f) return;

    const ImVec4 oc(o.color[0], o.color[1], o.color[2], o.color[3]);
    auto col = [&](ImVec4 v, float a = 1.0f) { return ImGui::ColorConvertFloat4ToU32(ImVec4(v.x, v.y, v.z, v.w * a * c.alpha)); };
    const ImU32 shadow = IM_COL32(0, 0, 0, (int)(160.0f * o.color[3] * c.alpha));
    const float soff = MAX(fpx * 0.05f, 1.0f);
    auto text = [&](ImVec2 p, ImU32 colr, const char* s) {
        dl->AddText(font, fpx, ImVec2(p.x + soff, p.y + soff), shadow, s);
        dl->AddText(font, fpx, p, colr, s);
    };

    // What the series hold
    SeriesTemporalView tv[6];
    SeriesHistogramView hv[6];
    bool has_t[6] = {}, has_h[6] = {};
    for (int s = 0; s < ns; ++s) {
        has_t[s] = series_resolve_temporal(&tv[s], state, sp.series[s].key);
        if (!timeline) has_h[s] = series_resolve_histogram(&hv[s], state, sp.series[s].key, movie_distribution_bins(o, sp.series[s].num_bins));
    }

    // The axes
    double xa = 0.0, xb = 1.0, ya = 0.0, yb = 1.0;
    bool flip = false;   // A trajectory time axis for a movie that plays the trajectory backward
    if (timeline) {
        if (elapsed) {
            xb = movie_time_bar_moved(prof, c.axis_end) - d_axis0;
        } else {
            movie_time_bar_visited_between(prof, c.axis_begin, c.axis_end, &xa, &xb);
            flip = movie_trajectory_quantity(state, c.axis_end) < movie_trajectory_quantity(state, c.axis_begin);
            if (!(xb > xa)) {
                for (int s = 0; s < ns; ++s) {
                    if (has_t[s] && tv[s].num_samples > 1) { xa = tv[s].x[0]; xb = tv[s].x[tv[s].num_samples - 1]; break; }
                }
            }
        }
        bool any = false;
        double lo = 0.0, hi = 0.0;
        for (int s = 0; s < ns; ++s) {
            if (!has_t[s] || tv[s].num_samples < 1) continue;
            const int i0 = CLAMP((int)floor(series_temporal_index_at(tv[s], full_lo)), 0, tv[s].num_samples - 1);
            const int i1 = CLAMP((int)ceil(series_temporal_index_at(tv[s], full_hi)), 0, tv[s].num_samples - 1);
            for (int i = i0; i <= i1; ++i) {
                const double y = (double)tv[s].y[(size_t)i * tv[s].stride] * tv[s].y_scale;
                if (!any) { lo = hi = y; any = true; }
                else { lo = fmin(lo, y); hi = fmax(hi, y); }
            }
        }
        if (!any || !(xb > xa)) return;
        const double span = hi > lo ? hi - lo : MAX(fabs(hi), 1.0) * 0.1;
        ya = lo - 0.06 * span;
        yb = hi + 0.06 * span;
    } else {
        bool any = false;
        for (int s = 0; s < ns; ++s) {
            if (!has_h[s] || !(hv[s].x_max > hv[s].x_min)) continue;
            if (!any) { xa = hv[s].x_min; xb = hv[s].x_max; any = true; }
            else { xa = fmin(xa, hv[s].x_min); xb = fmax(xb, hv[s].x_max); }
        }
        if (!any) return;
        ya = 0.0;
        yb = 1.12;
    }
    auto sx = [&](double x) { const float t = (float)((x - xa) / (xb - xa)); return ix0 + (flip ? 1.0f - t : t) * (ix1 - ix0); };
    auto sy = [&](double y) { return iy1 - (float)((y - ya) / (yb - ya)) * (iy1 - iy0); };

    // Ticks and frame
    const ImU32 col_axis = col(oc, 0.9f);
    const ImU32 col_grid = col(oc, 0.15f);
    {
        double ticks[16], step = 0.0;
        const int nt = movie_nice_ticks(xa, xb, MAX((int)((ix1 - ix0) / (fpx * 5.0f)), 2), ticks, 16, &step);
        const int dec = movie_tick_decimals(step);
        for (int i = 0; i < nt; ++i) {
            const float x = sx(ticks[i]);
            dl->AddLine(ImVec2(x, iy0), ImVec2(x, iy1), col_grid);
            if (!x_labels) continue;
            char b[32];
            snprintf(b, sizeof(b), "%.*f", dec, ticks[i]);
            dl->AddLine(ImVec2(x, iy1), ImVec2(x, iy1 + fpx * 0.3f), col_axis);
            const float tw = font->CalcTextSizeA(fpx, FLT_MAX, 0.0f, b).x;
            text(ImVec2(CLAMP(x - tw * 0.5f, p0.x, p0.x + w - tw), iy1 + fpx * 0.35f), col_axis, b);
        }
        // What the axis is measured in
        char u[48] = "";
        if (timeline) {
            char unit[32] = "";
            if (md_array_size(state->timeline.x_values) > 0 && !md_unit_is_none(state->timeline.time_unit)) md_unit_print(unit, sizeof(unit), state->timeline.time_unit);
            else if (elapsed) snprintf(unit, sizeof(unit), "frames");
            if (elapsed) snprintf(u, sizeof(u), "elapsed %s", unit);
            else snprintf(u, sizeof(u), "%s", unit);
        } else {
            for (int s = 0; s < ns && !u[0]; ++s) if (has_h[s]) snprintf(u, sizeof(u), "%s", hv[s].x_unit_str);
        }
        if (u[0] && x_labels) {
            const float tw = font->CalcTextSizeA(fpx, FLT_MAX, 0.0f, u).x;
            text(ImVec2(ix1 - tw, iy1 + fpx * 1.5f), col_axis, u);
        }
    }
    if (timeline) {
        double ticks[16], step = 0.0;
        const int nt = movie_nice_ticks(ya, yb, MAX((int)((iy1 - iy0) / (fpx * 2.5f)), 2), ticks, 16, &step);
        const int dec = movie_tick_decimals(step);
        for (int i = 0; i < nt; ++i) {
            char b[32];
            snprintf(b, sizeof(b), "%.*f", dec, ticks[i]);
            const float y = sy(ticks[i]);
            dl->AddLine(ImVec2(ix0, y), ImVec2(ix1, y), col_grid);
            const float tw = font->CalcTextSizeA(fpx, FLT_MAX, 0.0f, b).x;
            text(ImVec2(ix0 - fpx * 0.3f - tw, y - fpx * 0.5f), col_axis, b);
        }
    }
    dl->AddLine(ImVec2(ix0, iy1), ImVec2(ix1, iy1), col_axis, MAX(fpx * 0.08f, 1.0f));
    dl->AddLine(ImVec2(ix0, iy0), ImVec2(ix0, iy1), col_axis, MAX(fpx * 0.08f, 1.0f));

    dl->PushClipRect(ImVec2(ix0, iy0 - 2.0f), ImVec2(ix1 + 2.0f, iy1 + 2.0f), true);
    std::vector<ImVec2> pts;
    std::vector<MovieCurvePoint> curve;
    std::vector<float> counts_full, counts_vis;
    const double s_now = movie_time_bar_moved(prof, tnow) - d_axis0;
    for (int s = 0; s < ns; ++s) {
        const PlotSeries& ps = sp.series[s];
        const ImVec4 sc = series_color(s);
        const float thick = c.thick;
        if (timeline) {
            if (!has_t[s] || tv[s].num_samples < 2) continue;
            const SeriesTemporalView& v = tv[s];
            pts.clear();
            double cursor_x = c.q;
            if (elapsed) {
                // The series along the path of the movie, to where it is (or to the end without 'reveal')
                movie_elapsed_curve_between(&curve, prof, pt0, o.reveal ? tnow : pt1, v.x, v.num_samples,
                    [&](int i) { return (double)v.y[(size_t)i * v.stride] * v.y_scale; }, [&](double x) { return movie_series_value_at(v, x); });
                const size_t step = MAX(curve.size() / (size_t)(4.0f * MAX(ix1 - ix0, 1.0f)), (size_t)1);
                for (size_t i = 0; i < curve.size(); i += step) pts.push_back(ImVec2(sx(curve[i].s), sy(curve[i].v)));
                if (step > 1 && !curve.empty()) pts.push_back(ImVec2(sx(curve.back().s), sy(curve.back().v)));
                cursor_x = s_now;
            } else {
                if (!(vis_hi >= vis_lo)) continue;
                const int j0 = CLAMP((int)ceil(series_temporal_index_at(v, vis_lo)), 0, v.num_samples - 1);
                const int j1 = CLAMP((int)floor(series_temporal_index_at(v, vis_hi)), 0, v.num_samples - 1);
                pts.push_back(ImVec2(sx(vis_lo), sy(movie_series_value_at(v, vis_lo))));
                const int step = MAX((j1 - j0) / MAX((int)(ix1 - ix0), 1), 1);
                for (int i = j0; i <= j1; i += step) pts.push_back(ImVec2(sx((double)v.x[i]), sy((double)v.y[(size_t)i * v.stride] * v.y_scale)));
                pts.push_back(ImVec2(sx(vis_hi), sy(movie_series_value_at(v, vis_hi))));
            }
            if (pts.size() < 2) continue;
            if (ps.plot_type == PlotType_Scatter) {
                for (const ImVec2& p : pts) dl->AddCircleFilled(p, thick, col(sc), 6);
            } else {
                dl->AddPolyline(pts.data(), (int)pts.size(), col(sc), ImDrawFlags_None, thick);
            }
            // Where the frame that is shown is
            const float mx = sx(cursor_x), my = sy(movie_series_value_at(v, c.q));
            dl->AddLine(ImVec2(mx, iy0), ImVec2(mx, iy1), col(sc, 0.35f));
            dl->AddCircleFilled(ImVec2(mx, my), thick * 2.0f, col(ImVec4(0, 0, 0, 1), 0.7f), 12);
            dl->AddCircleFilled(ImVec2(mx, my), thick * 1.5f, col(sc), 12);
        } else {
            if (!has_h[s]) continue;
            const SeriesHistogramView& hs = hv[s];
            const int nb = hs.num_bins;
            if (nb < 1 || !(hs.x_max > hs.x_min)) continue;
            const float* heights = nullptr;
            float peak = 0.0f;
            std::vector<float> own;
            if (o.reveal && has_t[s]) {
                // Counted over the frames that have been played, in the bins of the full distribution
                movie_histogram_counts(&counts_full, nb, hs.x_min, hs.x_max, tv[s].x, tv[s].y, tv[s].stride, tv[s].y_scale, tv[s].num_samples, full_lo, full_hi);
                movie_histogram_counts(&counts_vis, nb, hs.x_min, hs.x_max, tv[s].x, tv[s].y, tv[s].stride, tv[s].y_scale, tv[s].num_samples, vis_lo, vis_hi);
                for (float f : counts_full) peak = MAX(peak, f);
                heights = counts_vis.data();
            } else {
                own.assign(hs.bins, hs.bins + nb);
                for (float f : own) peak = MAX(peak, f);
                heights = own.data();
            }
            if (!(peak > 0.0f)) continue;
            const double bw = (hs.x_max - hs.x_min) / (double)nb;
            const ImU32 fill = col(sc, 0.45f), edge = col(sc);
            for (int b = 0; b < nb; ++b) {
                if (heights[b] <= 0.0f) continue;
                const float x0 = sx(hs.x_min + bw * b), x1 = sx(hs.x_min + bw * (b + 1));
                const float y = sy(1.0 * (double)heights[b] / (double)peak);
                dl->AddRectFilled(ImVec2(x0, y), ImVec2(MAX(x1, x0 + 1.0f), iy1), fill);
                dl->AddLine(ImVec2(x0, y), ImVec2(MAX(x1, x0 + 1.0f), y), edge, MAX(fpx * 0.1f, 1.0f));
            }
            if (has_t[s]) {
                const float mx = sx(movie_series_value_at(tv[s], c.q));
                dl->AddLine(ImVec2(mx, iy0), ImVec2(mx, iy1), edge, MAX(fpx * 0.12f, 1.0f));
            }
        }
    }
    // The markers of the movie, where the movie gets to them
    if (timeline && o.show_markers) {
        // Labels sit right of their marker in the first half of the plot and left of it (ending at it) in the second; one that
        // would overlap a label already written goes a line lower (in time order, so they do not jump as the movie plays)
        MovieLabelSpan placed[64];
        int num_placed = 0;
        const int max_rows = MAX(1, MIN(3, (int)((iy1 - iy0) / (fpx * 1.15f)) - 1));
        md_array(const MovieMarker*) order = 0;
        for (const MovieMarker& mk : state->movie.markers) md_array_push(order, &mk, frame_alloc);
        std::sort(order, order + md_array_size(order), [](const MovieMarker* a, const MovieMarker* b) { return a->time < b->time; });
        for (size_t mi = 0; mi < md_array_size(order); ++mi) {
            const MovieMarker& mk = *order[mi];
            if (!movie_marker_matches_subplot(mk, sp.id)) continue;
            if ((o.reveal && mk.time > c.time) || mk.time < c.axis_begin || mk.time > c.axis_end) continue;
            const float x = sx(elapsed ? movie_time_bar_moved(prof, mk.time) - d_axis0 : movie_trajectory_quantity(state, mk.time));
            if (x < ix0 || x > ix1) continue;
            float mc[4];
            movie_marker_color(state->movie.markers.data(), state->movie.markers.size(), (size_t)(&mk - state->movie.markers.data()), mc);
            const ImVec4 mcol(mc[0], mc[1], mc[2], mc[3]);
            dl->AddLine(ImVec2(x, iy0), ImVec2(x, iy1), col(mcol, 0.8f), MAX(fpx * 0.1f, 1.0f));
            if ((marker_labels || mk.subplot != 0) && mk.label[0]) {
                const float tw = font->CalcTextSizeA(fpx, FLT_MAX, 0.0f, mk.label).x;
                float lx = x;
                const int row = movie_label_place(placed, (size_t)num_placed, x, tw, ix0, ix1, fpx * 0.3f, fpx * 0.4f, max_rows, &lx);
                if (num_placed < (int)ARRAY_SIZE(placed)) placed[num_placed++] = { lx, lx + tw, row };
                text(ImVec2(lx, iy0 + fpx * 0.1f + (float)row * fpx * 1.15f), col(mcol), mk.label);
            }
        }
    }
    dl->PopClipRect();

    // The legend, with the value at the frame that is shown
    if (titled) text(ImVec2(p0.x, p0.y), col(oc), title);
    for (int s = 0; s < ns; ++s) {
        const float y = p0.y + title_h + (float)s * fpx * 1.2f;
        dl->AddRectFilled(ImVec2(p0.x, y + fpx * 0.2f), ImVec2(p0.x + fpx * 0.7f, y + fpx * 0.9f), col(series_color(s)));
        char b[160];
        const int n = snprintf(b, sizeof(b), "%s", tv[s].label[0] ? tv[s].label : hv[s].label);
        if (o.show_value && has_t[s] && tv[s].num_samples > 0) {
            const int i = CLAMP((int)(series_temporal_index_at(tv[s], c.q) + 0.5), 0, tv[s].num_samples - 1);
            char vb[48];
            series_temporal_print_value(vb, sizeof(vb), tv[s], i, 0);
            snprintf(b + n, sizeof(b) - (size_t)n, "   %s %s", vb, tv[s].unit_str);
        }
        text(ImVec2(p0.x + fpx * 1.0f, y), col(oc), b);
    }
}

// The panels of a figure whose subplot is there and has series in it, at most 2 * PLOT_MAX_SUBPLOTS
static size_t movie_figure_panels(const MovieOverlay& o, const ApplicationState* state, const PlotSubplot** shown, const MoviePlotPanel** shown_panel, MoviePlotView* views) {
    size_t k = 0;
    for (const MoviePlotPanel& panel : o.panels) {
        if ((panel.view == MoviePlotView::Timeline && o.type == MovieOverlayType::Distribution) || (panel.view == MoviePlotView::Distribution && o.type == MovieOverlayType::Timeline)) continue;
        const bool tl = panel.view == MoviePlotView::Timeline;
        const PlotSubplot* subs = tl ? state->timeline.subplots : state->distributions.subplots;
        const int idx = plot_find_subplot(subs, tl ? state->timeline.num_subplots : state->distributions.num_subplots, panel.subplot);
        if (idx < 0 || subs[idx].count == 0 || k >= 2 * PLOT_MAX_SUBPLOTS) continue;
        if (shown) shown[k] = &subs[idx];
        if (shown_panel) shown_panel[k] = &panel;
        if (views) views[k] = panel.view;
        k += 1;
    }
    return k;
}

// Where a figure of 'k' panels is on the frame: its top left corner, size, text size and plate padding
static void movie_figure_layout(const MovieOverlay& o, size_t k, ImVec2 pos, ImVec2 size, float margin, ImVec2* p0, float* plot_w, float* plot_h, float* fpx, float* pad) {
    *plot_h = MAX(movie_overlay_size_px(o, size.y), 40.0f);
    *plot_w = MAX(o.width * size.x, 80.0f);
    *fpx    = o.font_points > 0.0f ? CLAMP(o.font_points * size.y / MOVIE_OVERLAY_POINT_REFERENCE_HEIGHT, 4.0f, 400.0f)
                                   : CLAMP(*plot_h * 0.085f / sqrtf((float)MAX(k, (size_t)1)), 6.0f, 64.0f);
    *pad    = *fpx * 0.4f;
    const int ai = (int)o.anchor;
    *p0 = ImVec2(pos.x + margin + (size.x - 2.0f * margin - *plot_w) * 0.5f * (float)(ai % 3),
                 pos.y + margin + (size.y - 2.0f * margin - *plot_h) * 0.5f * (float)(ai / 3));
}

// A figure: the subplots of the Timelines and the Distributions windows that it holds, stacked, with one horizontal axis for
// timelines that are above each other. They fill a block of the overlay's width and height at its anchor, moved down by 'dy'.
static void movie_figure_draw(ImDrawList* dl, ImFont* font, const MovieOverlay& o, ImVec2 pos, ImVec2 size, float margin, float dy, float alpha, double time, const ApplicationState* cstate) {
    // Resolving a series fills the caches of the application, which is why it takes it mutable
    ApplicationState* state = const_cast<ApplicationState*>(cstate);

    const PlotSubplot* shown[2 * PLOT_MAX_SUBPLOTS];
    const MoviePlotPanel* shown_panel[2 * PLOT_MAX_SUBPLOTS];
    MoviePlotView views[2 * PLOT_MAX_SUBPLOTS];
    const size_t k = movie_figure_panels(o, state, shown, shown_panel, views);
    if (k == 0) return;

    ImVec2 p0;
    float plot_w, plot_h, fpx, pad;
    movie_figure_layout(o, k, pos, size, margin, &p0, &plot_w, &plot_h, &fpx, &pad);
    p0.y += dy;
    // The panels come in and go at times of their own; a panel that is not there yet keeps its place, so the others do not move
    const MovieTimeBarProfile& prof = movie_time_bar_profile_for(state);
    const float thick = o.line_points > 0.0f ? MAX(o.line_points * size.y / MOVIE_OVERLAY_POINT_REFERENCE_HEIGHT, 1.0f) : MAX(fpx * 0.16f, 1.5f);
    MoviePlotContext c = { dl, font, &o, state, &prof, fpx, alpha, time, movie_trajectory_quantity(state, time), 0.0, prof.duration, thick, 0 };

    float panel_alpha[2 * PLOT_MAX_SUBPLOTS];
    MoviePlotView seen_views[2 * PLOT_MAX_SUBPLOTS];
    size_t seen = 0;
    for (size_t i = 0; i < k; ++i) {
        panel_alpha[i] = movie_panel_alpha(*shown_panel[i], o, time);
        if (panel_alpha[i] > 0.0f) seen_views[seen++] = views[i];
    }

    const float panel_h = plot_h / (float)k;

    // The plate is behind the panels that are there, so it grows when another comes in (with its fade) instead of showing an empty slot
    if (o.background[3] > 0.0f) {
        const float round = fpx * 0.4f;
        for (size_t i = 0; i < k; ++i) {
            if (panel_alpha[i] <= 0.0f) continue;
            const bool above = i > 0 && panel_alpha[i - 1] > 0.0f, below = i + 1 < k && panel_alpha[i + 1] > 0.0f;
            const float y0 = p0.y + (float)i * panel_h - (above ? 0.0f : pad);
            const float y1 = p0.y + (float)(i + 1) * panel_h + (below ? 0.0f : pad);
            const ImDrawFlags flags = (above ? 0 : ImDrawFlags_RoundCornersTop) | (below ? 0 : ImDrawFlags_RoundCornersBottom);
            dl->AddRectFilled(ImVec2(p0.x - pad, y0), ImVec2(p0.x + plot_w + pad, y1),
                ImGui::ColorConvertFloat4ToU32(ImVec4(o.background[0], o.background[1], o.background[2], o.background[3] * panel_alpha[i])), flags ? round : 0.0f, flags ? flags : ImDrawFlags_RoundCornersNone);
        }
    }

    bool labelled_markers = false;
    size_t seen_index = 0;
    for (size_t i = 0; i < k; ++i) {
        if (panel_alpha[i] <= 0.0f) continue;
        const bool first_timeline = views[i] == MoviePlotView::Timeline && !labelled_markers;
        labelled_markers |= first_timeline;
        c.alpha = panel_alpha[i];   // Already inside the fades of the overlay
        movie_plot_panel_draw(c, *shown[i], shown_panel[i]->title, views[i], ImVec2(p0.x, p0.y + (float)i * panel_h), plot_w, panel_h, movie_figure_x_labels(seen_views, seen, seen_index), first_timeline);
        seen_index += 1;
    }
}

// The text of a text, time stamp or scale bar overlay at a movie time, and the length of the scale bar in pixels. False when
// there is nothing to draw.
static bool movie_overlay_text(const MovieOverlay& o, ImVec2 size, double time, const ApplicationState* state, char* buf, size_t cap, float* bar_px) {
    buf[0] = '\0';
    *bar_px = 0.0f;
    switch (o.type) {
    case MovieOverlayType::Text:
        snprintf(buf, cap, "%s", o.text);
        break;
    case MovieOverlayType::Timestamp: {
        // The time that has gone since the movie started, which only grows whichever way the trajectory plays
        const double elapsed = movie_time_bar_moved(movie_time_bar_profile_for(state), time);
        if (md_array_size(state->timeline.x_values) > 0) {
            char unit_buf[32] = "";
            if (!md_unit_is_none(state->timeline.time_unit)) md_unit_print(unit_buf, sizeof(unit_buf), state->timeline.time_unit);
            snprintf(buf, cap, "%.1f %s", elapsed, unit_buf);
        } else {
            snprintf(buf, cap, "%d frames", (int)(elapsed + 0.5));
        }
        break;
    }
    case MovieOverlayType::ScaleBar: {
        const double upp = movie_units_per_pixel(state->view.camera.distance, state->view.camera.fov_y, size.y);
        const float len = o.length > 0.0f ? o.length : movie_scale_bar_length(upp, size.x, 0.2);
        if (len <= 0.0f || upp <= 0.0) return false;
        *bar_px = (float)((double)len / upp);
        snprintf(buf, cap, "%g \xC3\x85", (double)len);
        break;
    }
    default: return false;
    }
    return true;
}

// The block (without plate) of a text, time stamp or scale bar overlay
static ImVec2 movie_overlay_text_block(ImFont* font, const MovieOverlay& o, float font_px, const char* buf, float bar_px, ImVec2* text_size) {
    *text_size = font->CalcTextSizeA(font_px, FLT_MAX, 0.0f, buf);
    const float bar_h = MAX(font_px * 0.18f, 2.0f);
    const float gap   = font_px * 0.2f;
    return o.type == MovieOverlayType::ScaleBar ? ImVec2(MAX(bar_px, text_size->x), text_size->y + gap + bar_h) : *text_size;
}

// The size of the time bar's block, and whether it has labels
static ImVec2 movie_time_bar_block(const MovieOverlay& o, ImVec2 size, float font_px) {
    const float bar_w = MAX(o.width * size.x, 8.0f);
    const float bar_h = MAX(font_px * 0.3f, 3.0f);
    const bool labels = o.show_elapsed || o.show_speed;
    return ImVec2(bar_w, (labels ? font_px + font_px * 0.25f : 0.0f) + bar_h);
}

static ImVec2 movie_anchor_pos(MovieOverlayAnchor anchor, ImVec2 pos, ImVec2 size, float margin, ImVec2 block) {
    const int a = (int)anchor;
    return ImVec2(pos.x + margin + (size.x - 2.0f * margin - block.x) * 0.5f * (float)(a % 3),
                  pos.y + margin + (size.y - 2.0f * margin - block.y) * 0.5f * (float)(a / 3));
}

// The rectangle an overlay takes on the frame, with its plate, at its anchor. False for one that is not on the frame.
static bool movie_overlay_rect(const MovieOverlay& o, ImFont* font, ImVec2 pos, ImVec2 size, float margin, double time, const ApplicationState* state, ImVec2* r0, ImVec2* r1) {
    const bool plate = o.background[3] > 0.0f;
    ImVec2 block, pad;
    switch (o.type) {
    case MovieOverlayType::TimeBar: {
        const float font_px = MAX(movie_overlay_size_px(o, size.y), 4.0f);
        block = movie_time_bar_block(o, size, font_px);
        pad = plate ? ImVec2(font_px * 0.3f, font_px * 0.18f) : ImVec2(0, 0);
        break;
    }
    case MovieOverlayType::Timeline:
    case MovieOverlayType::Distribution:
    case MovieOverlayType::Figure: {
        const size_t k = movie_figure_panels(o, state, nullptr, nullptr, nullptr);
        if (k == 0) return false;
        ImVec2 p0;
        float plot_w, plot_h, fpx, fpad;
        movie_figure_layout(o, k, pos, size, margin, &p0, &plot_w, &plot_h, &fpx, &fpad);
        const float e = plate ? fpad : 0.0f;
        *r0 = ImVec2(p0.x - e, p0.y - e);
        *r1 = ImVec2(p0.x + plot_w + e, p0.y + plot_h + e);
        return true;
    }
    case MovieOverlayType::Logo:
    case MovieOverlayType::Image: {
        float aspect = 4.0f;
        const GLuint tex = o.type == MovieOverlayType::Logo ? movie_logo_texture(&aspect) : movie_image_texture(o.path, &aspect);
        if (!tex) return false;
        const float logo_h = MAX(movie_overlay_size_px(o, size.y), 4.0f);
        block = ImVec2(logo_h * aspect, logo_h);
        pad = plate ? ImVec2(logo_h * 0.15f, logo_h * 0.15f) : ImVec2(0, 0);
        break;
    }
    case MovieOverlayType::Text:
    case MovieOverlayType::Timestamp:
    case MovieOverlayType::ScaleBar: {
        char buf[160];
        float bar_px = 0.0f;
        if (!movie_overlay_text(o, size, time, state, buf, sizeof(buf), &bar_px)) return false;
        const float font_px = MAX(movie_overlay_size_px(o, size.y), 4.0f);
        ImVec2 ts;
        block = movie_overlay_text_block(font, o, font_px, buf, bar_px, &ts);
        pad = plate ? ImVec2(font_px * 0.3f, font_px * 0.18f) : ImVec2(0, 0);
        break;
    }
    default: return false;
    }
    const ImVec2 p0 = movie_anchor_pos(o.anchor, pos, size, margin, block);
    *r0 = ImVec2(p0.x - pad.x, p0.y - pad.y);
    *r1 = ImVec2(p0.x + block.x + pad.x, p0.y + block.y + pad.y);
    return true;
}

// How far each overlay is moved down (up when negative) so that overlays shown at the same time do not overlap
static void movie_overlays_offsets(ImFont* font, ImVec2 pos, ImVec2 size, float margin, double time, const ApplicationState* state, float* dy) {
    const auto& ovs = state->movie.overlays;
    md_allocator_i* arena = frame_alloc;
    MovieOverlayBox* boxes = (MovieOverlayBox*)md_alloc(arena, sizeof(MovieOverlayBox) * MAX(ovs.size(), (size_t)1));
    size_t* index = (size_t*)md_alloc(arena, sizeof(size_t) * MAX(ovs.size(), (size_t)1));
    float* box_dy = (float*)md_alloc(arena, sizeof(float) * MAX(ovs.size(), (size_t)1));
    size_t n = 0;
    for (size_t i = 0; i < ovs.size(); ++i) {
        dy[i] = 0.0f;
        const MovieOverlay& o = ovs[i];
        if (!o.enabled || o.end <= o.begin) continue;
        ImVec2 r0, r1;
        if (!movie_overlay_rect(o, font, pos, size, margin, time, state, &r0, &r1)) continue;
        boxes[n] = MovieOverlayBox{ (int)o.anchor, r0.x, r0.y, r1.x, r1.y, o.begin, o.end };
        index[n++] = i;
    }
    movie_overlay_avoid(boxes, n, 0.012f * size.y, pos.y, pos.y + size.y, box_dy);
    for (size_t j = 0; j < n; ++j) dy[index[j]] = box_dy[j];
}

// The overlays that are visible at a movie time, drawn into the rectangle (pos, size) of a frame. Everything scales with the
// height of the frame, so a frame looks the same at any resolution. Overlays shown at the same time are moved apart.
static void movie_overlays_draw(ImDrawList* dl, ImVec2 pos, ImVec2 size, double time, const ApplicationState* state) {
    const auto& m = state->movie;
    ImFont* font = ImGui::GetFont();
    if (!font || size.x <= 0.0f || size.y <= 0.0f) return;

    const float margin = 0.035f * size.y;
    float* offsets = (float*)md_alloc(frame_alloc, sizeof(float) * MAX(m.overlays.size(), (size_t)1));
    movie_overlays_offsets(font, pos, size, margin, time, state, offsets);
    for (size_t oi = 0; oi < m.overlays.size(); ++oi) {
        const MovieOverlay& o = m.overlays[oi];
        const float dy = offsets[oi];
        const float alpha = movie_overlay_alpha(o, time);
        if (alpha <= 0.0f) continue;

        if (o.type == MovieOverlayType::TimeBar) {
            const MovieTimeBarProfile& profile = movie_time_bar_profile_for(state);
            const float font_px = MAX(movie_overlay_size_px(o, size.y), 4.0f);
            const float bar_w = MAX(o.width * size.x, 8.0f);
            const float bar_h = MAX(font_px * 0.3f, 3.0f);
            const float gap = font_px * 0.25f;
            const bool labels = o.show_elapsed || o.show_speed;
            const float label_h = labels ? font_px : 0.0f;
            const ImVec2 block = movie_time_bar_block(o, size, font_px);
            ImVec2 p0 = movie_anchor_pos(o.anchor, pos, size, margin, block);
            p0.y += dy;
            const float progress = (float)CLAMP(movie_time_bar_progress(profile, time), 0.0, 1.0);

            const ImU32 col_fill  = ImGui::ColorConvertFloat4ToU32(ImVec4(o.color[0], o.color[1], o.color[2], o.color[3] * alpha));
            const ImU32 col_track = ImGui::ColorConvertFloat4ToU32(ImVec4(0.4f, 0.4f, 0.4f, o.color[3] * alpha * 0.7f));
            const ImU32 col_shadow = IM_COL32(0, 0, 0, (int)(160.0f * o.color[3] * alpha));
            if (o.background[3] > 0.0f) {
                const float pad = font_px * 0.3f;
                const ImU32 plate = ImGui::ColorConvertFloat4ToU32(ImVec4(o.background[0], o.background[1], o.background[2], o.background[3] * alpha));
                dl->AddRectFilled(ImVec2(p0.x - pad, p0.y - pad * 0.6f), ImVec2(p0.x + block.x + pad, p0.y + block.y + pad * 0.6f), plate, font_px * 0.25f);
            }
            if (labels) {
                const float soff = MAX(font_px * 0.05f, 1.0f);
                if (o.show_elapsed) {
                    char text[96];
                    if (md_array_size(state->timeline.x_values) > 0) {
                        char unit_buf[32] = "";
                        if (!md_unit_is_none(state->timeline.time_unit)) md_unit_print(unit_buf, sizeof(unit_buf), state->timeline.time_unit);
                        snprintf(text, sizeof(text), "%.1f / %.1f %s", movie_time_bar_moved(profile, time), profile.total(), unit_buf);
                    } else {
                        snprintf(text, sizeof(text), "%.0f / %.0f frames", movie_time_bar_moved(profile, time), profile.total());
                    }
                    dl->AddText(font, font_px, ImVec2(p0.x + soff, p0.y + soff), col_shadow, text);
                    dl->AddText(font, font_px, p0, col_fill, text);
                }
                if (o.show_speed) {
                    // How fast the trajectory is played, as a multiple of the speed of the Animation panel
                    const double fps = fabs((double)state->animation.fps);
                    const double frames_per_s = movie_quantity_speed([state](double t) { return movie_trajectory_frame(state, t); }, time, 0.5 / MAX((double)m.fps, 1.0));
                    char text[32];
                    snprintf(text, sizeof(text), "x%.1f", fps > 0.0 ? frames_per_s / fps : 0.0);
                    const float tw = font->CalcTextSizeA(font_px, FLT_MAX, 0.0f, text).x;
                    const ImVec2 tp(p0.x + block.x - tw, p0.y);
                    dl->AddText(font, font_px, ImVec2(tp.x + soff, tp.y + soff), col_shadow, text);
                    dl->AddText(font, font_px, tp, col_fill, text);
                }
            }
            const float by = p0.y + label_h + (labels ? gap : 0.0f);
            dl->AddRectFilled(ImVec2(p0.x, by), ImVec2(p0.x + bar_w, by + bar_h), col_track, bar_h * 0.5f);
            if (progress > 0.0f) dl->AddRectFilled(ImVec2(p0.x, by), ImVec2(p0.x + MAX(bar_w * progress, bar_h), by + bar_h), col_fill, bar_h * 0.5f);
            continue;
        }

        if (o.type == MovieOverlayType::Timeline || o.type == MovieOverlayType::Distribution || o.type == MovieOverlayType::Figure) {
            movie_figure_draw(dl, font, o, pos, size, margin, dy, alpha, time, state);
            continue;
        }
        if (o.type == MovieOverlayType::PropertyVis) continue;   // In the viewport, not on the plane of the frame

        if (o.type == MovieOverlayType::Logo || o.type == MovieOverlayType::Image) {
            float aspect = 4.0f;
            const GLuint tex = o.type == MovieOverlayType::Logo ? movie_logo_texture(&aspect) : movie_image_texture(o.path, &aspect);
            if (!tex) continue;
            const float logo_h = MAX(movie_overlay_size_px(o, size.y), 4.0f);
            const ImVec2 logo_size(logo_h * aspect, logo_h);
            ImVec2 lp = movie_anchor_pos(o.anchor, pos, size, margin, logo_size);
            lp.y += dy;
            if (o.background[3] > 0.0f) {
                const float pad = logo_size.y * 0.15f;
                const ImU32 plate = ImGui::ColorConvertFloat4ToU32(ImVec4(o.background[0], o.background[1], o.background[2], o.background[3] * alpha));
                dl->AddRectFilled(ImVec2(lp.x - pad, lp.y - pad), ImVec2(lp.x + logo_size.x + pad, lp.y + logo_size.y + pad), plate, pad);
            }
            dl->AddImage((ImTextureID)(intptr_t)tex, lp, ImVec2(lp.x + logo_size.x, lp.y + logo_size.y), ImVec2(0, 0), ImVec2(1, 1),
                ImGui::ColorConvertFloat4ToU32(ImVec4(o.color[0], o.color[1], o.color[2], o.color[3] * alpha)));
            continue;
        }

        const float font_px = MAX(movie_overlay_size_px(o, size.y), 4.0f);
        const ImU32 col    = ImGui::ColorConvertFloat4ToU32(ImVec4(o.color[0], o.color[1], o.color[2], o.color[3] * alpha));
        const ImU32 shadow = IM_COL32(0, 0, 0, (int)(160.0f * o.color[3] * alpha));
        const float soff   = MAX(font_px * 0.05f, 1.0f);

        char buf[160];
        float bar_px = 0.0f;   // Scale bar only
        if (!movie_overlay_text(o, size, time, state, buf, sizeof(buf), &bar_px)) continue;
        ImVec2 text_size;
        const ImVec2 block = movie_overlay_text_block(font, o, font_px, buf, bar_px, &text_size);
        const float bar_h = MAX(font_px * 0.18f, 2.0f);
        const float gap   = font_px * 0.2f;
        ImVec2 p0 = movie_anchor_pos(o.anchor, pos, size, margin, block);
        p0.y += dy;

        if (o.background[3] > 0.0f) {
            const float pad = font_px * 0.3f;
            const ImU32 plate = ImGui::ColorConvertFloat4ToU32(ImVec4(o.background[0], o.background[1], o.background[2], o.background[3] * alpha));
            dl->AddRectFilled(ImVec2(p0.x - pad, p0.y - pad * 0.6f), ImVec2(p0.x + block.x + pad, p0.y + block.y + pad * 0.6f), plate, font_px * 0.25f);
        }

        if (o.type == MovieOverlayType::ScaleBar) {
            const float tx = p0.x + 0.5f * (block.x - text_size.x);
            dl->AddText(font, font_px, ImVec2(tx + soff, p0.y + soff), shadow, buf);
            dl->AddText(font, font_px, ImVec2(tx, p0.y), col, buf);
            const float bx = p0.x + 0.5f * (block.x - bar_px);
            const float by = p0.y + text_size.y + gap;
            dl->AddRectFilled(ImVec2(bx + soff, by + soff), ImVec2(bx + bar_px + soff, by + bar_h + soff), shadow);
            dl->AddRectFilled(ImVec2(bx, by), ImVec2(bx + bar_px, by + bar_h), col);
        } else {
            dl->AddText(font, font_px, ImVec2(p0.x + soff, p0.y + soff), shadow, buf);
            dl->AddText(font, font_px, p0, col, buf);
        }
    }
}

// Draws the overlays onto the frame in the G-buffer, after post processing and before it is read back
static void movie_overlays_render(ApplicationState* state, double time) {
    if (state->movie.overlays.empty()) return;

    const float w = (float)state->gbuffer.width;
    const float h = (float)state->gbuffer.height;

    ImDrawList dl(ImGui::GetDrawListSharedData());
    dl._ResetForNewFrame();
    dl.PushTexture(ImGui::GetIO().Fonts->TexRef);
    dl.PushClipRect(ImVec2(0, 0), ImVec2(w, h));
    movie_overlays_draw(&dl, ImVec2(0, 0), ImVec2(w, h), time, state);
    for (const MovieOverlay& o : state->movie.overlays) {
        // The labels of the visualization are drawn by the viewport, which the recording does not include
        if (o.type == MovieOverlayType::PropertyVis && movie_overlay_alpha(o, time) > 0.0f) {
            script_vis_text_draw(&dl, ImVec2(w, h), h / MAX((float)state->app.window.height, 1.0f), *state);
            break;
        }
    }
    if (dl.VtxBuffer.Size == 0) return;

    ImDrawData dd;
    dd.Valid = true;
    dd.CmdLists.push_back(&dl);
    dd.CmdListsCount = 1;
    dd.TotalVtxCount = dl.VtxBuffer.Size;
    dd.TotalIdxCount = dl.IdxBuffer.Size;
    dd.DisplayPos = ImVec2(0, 0);
    dd.DisplaySize = ImVec2(w, h);
    dd.FramebufferScale = ImVec2(1, 1);
    dd.Textures = &ImGui::GetPlatformIO().Textures;

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, state->gbuffer.fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    ImGui_ImplOpenGL3_RenderDrawData(&dd);
}

// Starts the read back of the frame that was just rendered. The pixels reach the sink some frames later,
// from movie_pbo_pop, so the GPU is not waited for here.
static void movie_capture_frame(ApplicationState* state) {
    auto& m = state->movie;
    ASSERT(m.pbo_count < MOVIE_RING_SIZE);

    if ((int)state->gbuffer.width != m.rec_w || (int)state->gbuffer.height != m.rec_h) {
        VIAMD_LOG_ERROR("Movie recording stopped: the render target is %dx%d, not the %dx%d of the movie",
            (int)state->gbuffer.width, (int)state->gbuffer.height, m.rec_w, m.rec_h);
        movie_recording_stop(state);
        return;
    }

    const int slot = (m.pbo_head + m.pbo_count) % MOVIE_RING_SIZE;
    glBindFramebuffer(GL_READ_FRAMEBUFFER, state->gbuffer.fbo);
    glReadBuffer(GL_COLOR_ATTACHMENT0);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, m.pbo[slot]);
    glReadPixels(0, 0, m.rec_w, m.rec_h, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    if (m.fence[slot]) glDeleteSync(m.fence[slot]);
    m.fence[slot] = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
    m.pbo_frame[slot] = m.frame_index;
    m.pbo_count += 1;
}

// Shows the frame that was just rendered in the viewport, scaled to fit. The GUI is drawn over it.
static void movie_blit_preview(ApplicationState* state) {
    const int fw = (int)state->app.framebuffer.width;
    const int fh = (int)state->app.framebuffer.height;
    const int gw = (int)state->gbuffer.width;
    const int gh = (int)state->gbuffer.height;

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    glDrawBuffer(GL_BACK);
    glViewport(0, 0, fw, fh);
    const float background[4] = {0.07f, 0.07f, 0.07f, 1.0f};
    glClearBufferfv(GL_COLOR, 0, background);

    if (fw <= 0 || fh <= 0 || gw <= 0 || gh <= 0) return;

    const float scale = MIN((float)fw / (float)gw, (float)fh / (float)gh);
    const int dw = MAX(1, (int)(gw * scale));
    const int dh = MAX(1, (int)(gh * scale));
    const int x0 = (fw - dw) / 2;
    const int y0 = (fh - dh) / 2;

    glBindFramebuffer(GL_READ_FRAMEBUFFER, state->gbuffer.fbo);
    glReadBuffer(GL_COLOR_ATTACHMENT0);
    glBlitFramebuffer(0, 0, gw, gh, x0, y0, x0 + dw, y0 + dh, GL_COLOR_BUFFER_BIT, GL_LINEAR);
}

// Begins a movie recording: the trajectory is stepped deterministically (independent of real
// render speed) from movie.start_frame to movie.end_frame. One frame is rendered per output frame
// at the movie's resolution and handed to a frame_sink, which writes it in the background.
static void movie_recording_start(ApplicationState* state) {
    ASSERT(state);
    auto& m = state->movie;

    if (m.state == MovieRecordingState::Recording) {
        return;
    }

    if (m.sink) {
        VIAMD_LOG_ERROR("Cannot start movie recording: the previous movie is still being written");
        return;
    }

    if (str_empty(m.output_dir)) {
        VIAMD_LOG_ERROR("Cannot start movie recording: no output directory selected");
        return;
    }

    const double max_frame = (double)(run_num_frames(state) > 0 ? run_num_frames(state) - 1 : 0);
    m.start_frame = CLAMP(m.start_frame, 0.0, max_frame);
    m.end_frame   = CLAMP(m.end_frame,   0.0, max_frame);
    if (movie_duration(state) <= 0.0) {
        VIAMD_LOG_ERROR("Cannot start movie recording: the movie has no duration");
        return;
    }

    if (m.filename_prefix[0] == '\0') {
        snprintf(m.filename_prefix, sizeof(m.filename_prefix), "frame");
    }

    int w = 0, h = 0;
    movie_frame_size(state, &w, &h);
    if (w <= 0 || h <= 0) {
        VIAMD_LOG_ERROR("Cannot start movie recording: the window has no size");
        return;
    }

    frame_sink::Desc desc;
    desc.kind   = movie_output_is_video(m.output) ? frame_sink::Kind::Ffmpeg : frame_sink::Kind::PngSequence;
    desc.codec  = m.output == MovieOutput::Mp4H265 ? frame_sink::Codec::H265 : m.output == MovieOutput::WebmVp9 ? frame_sink::Codec::Vp9 : frame_sink::Codec::H264;
    desc.dir    = m.output_dir;
    desc.prefix = str_from_cstr(m.filename_prefix);
    desc.width  = w;
    desc.height = h;
    desc.fps    = m.fps;
    desc.crf    = m.crf;
    desc.ffmpeg = str_from_cstr(m.ffmpeg_path);

    char err[256] = {};
    m.sink = frame_sink::create(desc, err, sizeof(err));
    if (!m.sink) {
        VIAMD_LOG_ERROR("Cannot start movie recording: %s%s", err,
            movie_output_is_video(m.output) ? ". The PNG sequence output does not need ffmpeg" : "");
        return;
    }

    movie_pbo_alloc(state, (size_t)w * (size_t)h * 4);

    if (movie_output_is_video(m.output)) {
        snprintf(m.rec_result, sizeof(m.rec_result), STR_FMT "/%s.%s", STR_ARG(m.output_dir), m.filename_prefix, frame_sink::file_extension(desc.codec));
    } else {
        snprintf(m.rec_result, sizeof(m.rec_result), STR_FMT, STR_ARG(m.output_dir));
    }

    // Nothing is highlighted in the frames but what an overlay of the movie shows
    md_bitfield_clear(&state->selection.highlight_mask);
    m.vis_highlight = false;

    m.rec_w = w;
    m.rec_h = h;
    m.rec_output = m.output;
    movie_render_range(state, &m.rec_first, &m.rec_last);
    m.cur_time = (double)m.rec_first / (double)m.fps;
    m.frame_index = m.rec_first;
    m.samples_done = 0;
    m.sample_target = state->visuals.temporal_aa.enabled ? (m.aa_samples > 0 ? CLAMP(m.aa_samples, 1, 256) : JITTER_SEQUENCE_SIZE) : 1;
    m.can_capture = true;
    m.paused = false;
    m.rec_active_s = 0.0;

    if (m.save_copy) {
        // Everything the movie was made from, so that it can be made again. Saving moves the workspace's name, which is put back.
        char copy[2048];
        snprintf(copy, sizeof(copy), STR_FMT "/%s." STR_FMT, STR_ARG(m.output_dir), m.filename_prefix, STR_ARG(WORKSPACE_FILE_EXTENSION));
        char prev[sizeof(state->files.workspace)];
        memcpy(prev, state->files.workspace, sizeof(prev));
        save_workspace(state, str_from_cstr(copy));
        memcpy(state->files.workspace, prev, sizeof(prev));
    }

    m.prev_playback_mode = state->animation.mode;
    state->animation.mode = PlaybackMode::Stopped;

    m.camera_was_animated = md_array_size(m.keyframes) > 0;
    m.prev_view_target = state->view.target;
    m.prev_fov_y = state->view.camera.fov_y;

    m.state = MovieRecordingState::Recording;

    VIAMD_LOG_INFO("Recording %d movie frame(s) of %d (%.2f s, %dx%d) to '%s'", m.rec_last - m.rec_first + 1, movie_num_frames(state), movie_duration(state), w, h, m.rec_result);
}

// Ends the recording. The frames captured so far are kept: the sink writes what it has and is
// closed, so a stopped movie is a shorter movie.
static void movie_recording_stop(ApplicationState* state) {
    ASSERT(state);
    auto& m = state->movie;

    if (m.state != MovieRecordingState::Recording) {
        return;
    }

    m.state = MovieRecordingState::Idle;
    m.paused = false;
    movie_restore_state(state);

    movie_pbo_flush(state);
    movie_pbo_free(state);
    if (m.sink) frame_sink::close(m.sink);

    VIAMD_LOG_INFO("Movie recording stopped after %d frame(s)", m.frame_index - m.rec_first);
}

static void movie_shutdown(ApplicationState* state) {
    movie_recording_stop(state);
    movie_pbo_free(state);
    if (state->movie.sink) {
        // Waits for whatever is still being written
        frame_sink::destroy(state->movie.sink);
        state->movie.sink = nullptr;
    }
}

// Plays the movie in the viewport at real time, with the trajectory and the camera as they will be recorded
static void update_movie_preview(ApplicationState* state) {
    auto& m = state->movie;
    if (!m.preview_playing) return;

    const float len = (float)movie_duration(state);
    if (m.state == MovieRecordingState::Recording || len <= 0.0f || run_num_frames(state) == 0 || (!m.show_window && !m.play_mode)) {
        m.preview_playing = false;
        return;
    }

    state->animation.mode = PlaybackMode::Stopped;
    m.playhead += (float)state->app.timing.delta_s;
    if (m.playhead >= len) {
        if (m.preview_loop) {
            m.playhead = fmodf(m.playhead, len);
        } else {
            m.playhead = len;
            m.preview_playing = false;
        }
    }
    movie_apply_time(state, (double)m.playhead, true);
}

// Drives the movie recording. Called once per application loop iteration, before render(): it sets the
// scene to the time of the frame being captured (and keeps it there while that frame's samples are
// accumulated), and passes finished read backs on to the sink.
static void update_movie_recording(ApplicationState* state) {
    ASSERT(state);
    auto& m = state->movie;

    if (m.state != MovieRecordingState::Recording) {
        // The recording is over but the sink may still be writing. Report when it is done.
        if (m.sink) {
            const frame_sink::Status st = frame_sink::status(m.sink);
            if (st.done) {
                if (st.ok) {
                    if (movie_output_is_video(m.rec_output)) {
                        VIAMD_LOG_SUCCESS("Movie saved to '%s' (%d frames)", m.rec_result, st.written);
                    } else {
                        VIAMD_LOG_SUCCESS("%d PNG frame(s) written to '%s'", st.written, m.rec_result);
                    }
                } else {
                    VIAMD_LOG_ERROR("Writing the movie failed: %s (%d of %d frame(s) written)", st.message, st.written, st.submitted);
                }
                frame_sink::destroy(m.sink);
                m.sink = nullptr;
            }
        }
        return;
    }

    // A writer that has failed (ffmpeg gone, disk full) ends the recording rather than rendering for nothing
    const frame_sink::Status st = frame_sink::status(m.sink);
    if (!st.ok) {
        VIAMD_LOG_ERROR("Movie recording stopped: %s", st.message);
        movie_recording_stop(state);
        return;
    }

    state->animation.mode = PlaybackMode::Stopped;

    while (movie_pbo_pop(state, false)) {}
    m.can_capture = m.pbo_count < MOVIE_RING_SIZE && !m.paused;
    if (!m.paused) m.rec_active_s += state->app.timing.delta_s;

    if (m.frame_index > m.rec_last) {
        const int num_frames = m.frame_index - m.rec_first;
        movie_pbo_flush(state);
        movie_pbo_free(state);
        frame_sink::close(m.sink);
        m.state = MovieRecordingState::Idle;
        movie_restore_state(state);
        VIAMD_LOG_INFO("Movie recording finished, %d frame(s) captured. Writing the rest in the background", num_frames);
        return;
    }

    // Computed from the index rather than accumulated, so the time does not drift. Applied on every
    // iteration so that nothing done in the meantime moves the scene away from this frame.
    m.cur_time = (double)m.frame_index / (double)m.fps;
    movie_apply_time(state, m.cur_time, m.camera_was_animated);
}

// ### MOVIE WINDOW ###

constexpr float MOVIE_DEG_TO_RAD = 3.14159265358979f / 180.0f;
constexpr float MOVIE_RAD_TO_DEG = 180.0f / 3.14159265358979f;

static const char* spin_axis_str[(int)SpinAxis::Count] = {
    "Camera up",
    "World Y",
    "World X",
    "World Z",
};

// A second key at the same time would leave the path undefined, so it replaces the first. With keep_extras
// the spin of the key it replaces stays, and so does its frame unless the new key has one.
static void movie_insert_key(ApplicationState* state, CameraKeyframe key, bool keep_extras) {
    auto& m = state->movie;
    for (size_t i = 0; i < md_array_size(m.keyframes); ++i) {
        const CameraKeyframe& old = m.keyframes[i];
        if (fabs(old.time - key.time) < 1.0e-3) {
            if (keep_extras) {
                key.spin_turns = old.spin_turns;
                key.spin_axis = old.spin_axis;
                key.spin_constant_speed = old.spin_constant_speed;
                if (!key.use_frame) {
                    key.use_frame = old.use_frame;
                    key.frame = old.frame;
                }
            }
            m.keyframes[i] = key;
            return;
        }
    }
    md_array_push(m.keyframes, key, state->allocator.persistent);
    movie_sort_keyframes(state);
}

// The pose a new key gets: what the viewport shows in the frame, in Scene view as in Movie preview
static void movie_camera_pose(const ApplicationState* state, ViewTransform* vt, float* fov_y) {
    *vt = state->view.target;
    *fov_y = state->view.camera.fov_y;
}

// A key at the preview time where the path already has the camera, so that the path keeps its shape (the view is not used)
static void movie_add_keyframe_on_path(ApplicationState* state) {
    auto& m = state->movie;
    const size_t n = md_array_size(m.keyframes);
    if (n == 0) {
        movie_add_keyframe(state);
        return;
    }
    const MovieResolved r = movie_resolve(state, m.keyframes, n);
    const double t = movie_snap_time(state, (double)m.playhead);
    ViewTransform pose;
    float fov_y;
    movie_pose(state, r, t, &pose, &fov_y);
    movie_insert_key(state, camera_key_on_path(pose, fov_y, r.keys, t, movie_upright(state)), false);
}

// The movie camera as it is now, at the playhead
static CameraKeyframe movie_current_key(ApplicationState* state) {
    auto& m = state->movie;
    CameraKeyframe key = {};
    movie_camera_pose(state, &key.transform, &key.fov_y);
    // With Keep upright the movie camera is levelled and tilted by the key's roll only: the view's own tilt goes into the roll,
    // so that the key plays back as it was seen
    if (const vec3_t* up = movie_upright(state)) key.roll = camera_roll(key.transform, *up);
    key.time = movie_snap_time(state, (double)m.playhead);
    if (m.key_includes_frame) {
        key.use_frame = true;
        key.frame = state->animation.frame;
    }
    return key;
}

static void movie_add_keyframe(ApplicationState* state) {
    movie_insert_key(state, movie_current_key(state), true);
}

// K / Add Keyframe: a key of the view at the preview time, or 2 s later when there is a key there already (and the preview time
// moves along), so that a key is never replaced by mistake
static void movie_add_keyframe_from_view(ApplicationState* state) {
    auto& m = state->movie;
    const double len = movie_duration(state);
    double t = movie_snap_time(state, (double)m.playhead);
    auto taken = [&](double time) {
        const double half_frame = 0.5 / MAX((double)m.fps, 1.0);
        for (size_t i = 0; i < md_array_size(m.keyframes); ++i) if (fabs(m.keyframes[i].time - time) <= half_frame) return true;
        return false;
    };
    while (taken(t)) {
        if (t + 2.0 > len + 1.0e-6) {
            VIAMD_LOG_ERROR("There is a key at the preview time and no room 2 s later: move the preview time or lengthen the movie");
            return;
        }
        t = movie_snap_time(state, t + 2.0);
    }
    m.playhead = (float)t;
    movie_add_keyframe(state);
}

// Remembers the key at the playhead (the nearest one within half a frame), to paste somewhere else
static void movie_copy_keyframe_at_playhead(ApplicationState* state) {
    auto& m = state->movie;
    const double half_frame = 0.5 / MAX((double)m.fps, 1.0);
    for (size_t i = 0; i < md_array_size(m.keyframes); ++i) {
        if (fabs(m.keyframes[i].time - (double)m.playhead) <= half_frame) {
            m.key_clipboard = m.keyframes[i];
            m.has_key_clipboard = true;
            return;
        }
    }
}

// The key that was copied from the table, at the playhead
static void movie_paste_keyframe(ApplicationState* state) {
    auto& m = state->movie;
    if (!m.has_key_clipboard) return;
    CameraKeyframe key = m.key_clipboard;
    key.time = movie_snap_time(state, (double)m.playhead);
    movie_insert_key(state, key, false);
}

// Adds a key that frames the selected atoms, seen from the direction the camera has now, and moves the view there
static void movie_add_selection_keyframe(ApplicationState* state) {
    const md_bitfield_t* mask = &state->selection.selection_mask;
    const size_t count = md_bitfield_popcount(mask);
    const size_t num_atoms = state->mold.sys.atom.count;
    uint64_t first = 0, last = 0;
    if (count == 0) {
        VIAMD_LOG_ERROR("Select some atoms to frame first");
        return;
    }
    if (state->mold.state.num_atoms != num_atoms || !md_bitfield_get_range(&first, &last, mask) || last >= num_atoms) {
        VIAMD_LOG_ERROR("The selection does not match the system");
        return;
    }

    // Placed together across periodic boundaries, so that a molecule split by the cell is framed as one
    md_temp_scope_t temp = md_temp_begin_in(state->allocator.frame);
    defer { md_temp_end(temp); };
    vec4_t* xyzw = md_temp_alloc_array(temp, vec4_t, count);
    vec3_t* xyz = md_temp_alloc_array(temp, vec3_t, count);
    md_util_system_extract_xyzw_from_mask(xyzw, mask, &state->mold.sys, &state->mold.state);
    vec3_t center = vec3_zero();
    md_util_deperiodize_self_vec4(xyzw, count, &state->mold.state.unitcell, &center);
    for (size_t i = 0; i < count; ++i) xyz[i] = vec3_from_vec4(xyzw[i]);

    CameraKeyframe key = movie_current_key(state);
    ViewTransform t = key.transform;
    t.distance = camera_fit_distance(xyz, nullptr, count, center, t.orientation, key.fov_y);
    t.position = camera_position_from_look_at(mat4_mul_vec3(state->mold.unitcell_transform, center, 1.0f), t.orientation, t.distance);
    key.transform = t;
    movie_insert_key(state, key, true);
    if (!movie_scene_view(state)) state->view.target = t;
}

// Two keys with the view as it is now, the second after the orbit's duration with whole turns around it.
// The camera leaves and comes back to the same pose.
static void movie_add_orbit(ApplicationState* state) {
    auto& m = state->movie;
    if (m.orbit_turns == 0) return;

    const double t0 = movie_snap_time(state, (double)m.playhead);
    const double t1 = movie_snap_time(state, t0 + (double)MAX(m.orbit_duration, 0.1f));
    if (t0 + (double)MAX(m.orbit_duration, 0.1f) > movie_duration(state) + 1.0e-6) {
        VIAMD_LOG_ERROR("The orbit ends at %.2f s, after the movie. Move the preview time earlier, shorten the orbit or increase the movie length", t0 + (double)MAX(m.orbit_duration, 0.1f));
        return;
    }

    movie_add_keyframe(state);

    CameraKeyframe end = {};
    for (size_t i = 0; i < md_array_size(m.keyframes); ++i) {
        if (fabs(m.keyframes[i].time - t0) < 1.0e-3) end = m.keyframes[i];
    }
    end.time = t1;
    end.spin_turns = m.orbit_turns;
    end.spin_axis = m.orbit_axis;
    end.spin_constant_speed = false;

    bool replaced = false;
    for (size_t i = 0; i < md_array_size(m.keyframes); ++i) {
        if (fabs(m.keyframes[i].time - t1) < 1.0e-3) {
            m.keyframes[i] = end;
            replaced = true;
        }
    }
    if (!replaced) {
        md_array_push(m.keyframes, end, state->allocator.persistent);
        movie_sort_keyframes(state);
    }
    m.playhead = (float)t1;
}

// Makes the camera path a loop: the movie ends in the pose it starts in, and the path is smooth across the seam
static void movie_close_loop(ApplicationState* state) {
    auto& m = state->movie;
    const size_t n = md_array_size(m.keyframes);
    if (n < 2) {
        VIAMD_LOG_ERROR("A loop needs at least two keyframes");
        return;
    }

    const CameraKeyframe first = m.keyframes[0];
    const double len = movie_duration(state);
    CameraKeyframe& last = m.keyframes[n - 1];
    if (last.time < len - 1.0e-3) {
        CameraKeyframe end = first;
        end.time = len;
        end.use_frame = false;
        end.spin_turns = 0;
        end.ease = KeyEase::Smooth;
        md_array_push(m.keyframes, end, state->allocator.persistent);
    } else {
        // The last key is at the end already: it becomes the first again
        last.transform = first.transform;
        last.fov_y = first.fov_y;
    }
    m.loop = true;
}

static void movie_param_sort(ApplicationState* state) {
    std::stable_sort(state->movie.param_keys.begin(), state->movie.param_keys.end(), [](const ParamKey& a, const ParamKey& b) {
        return a.param != b.param ? a.param < b.param : a.time < b.time;
    });
}

// Keys what the parameter is now, at the playhead
static void movie_key_param(ApplicationState* state, int id) {
    auto& m = state->movie;
    const MovieParamDesc* d = movie_param_desc(id);
    if (!d) return;

    ParamKey key;
    key.param = id;
    key.time = movie_snap_time(state, (double)m.playhead);
    const float* src = d->ptr(state);
    for (int c = 0; c < d->comps; ++c) key.value[c] = src[c];

    for (ParamKey& k : m.param_keys) {
        if (k.param == id && fabs(k.time - key.time) < 1.0e-3) {
            key.ease = k.ease;
            k = key;
            return;
        }
    }
    m.param_keys.push_back(key);
    movie_param_sort(state);
}

static void movie_undo(ApplicationState* state) {
    MovieKeys keys = movie_keys_snapshot(state);
    if (state->movie.history.undo(&keys)) movie_keys_restore(state, keys);
}

static void movie_redo(ApplicationState* state) {
    MovieKeys keys = movie_keys_snapshot(state);
    if (state->movie.history.redo(&keys)) movie_keys_restore(state, keys);
}

// Once per frame. What is edited is not reported, the history sees the keys change.
static void update_movie_history(ApplicationState* state) {
    if (state->movie.state == MovieRecordingState::Recording) return;
    const bool editing = ImGui::IsAnyItemActive() || ImGui::IsMouseDown(ImGuiMouseButton_Left);
    state->movie.history.update(movie_keys_snapshot(state), editing);
}

static void movie_goto_keyframe(ApplicationState* state, size_t idx) {
    auto& m = state->movie;
    if (idx >= md_array_size(m.keyframes)) return;
    const CameraKeyframe& key = m.keyframes[idx];
    if (!movie_scene_view(state)) {
        state->view.target = key.transform;
        if (const vec3_t* up = movie_upright(state)) camera_level(&state->view.target, *up, key.roll);
        state->view.camera.fov_y = key.fov_y;
    }
    m.playhead = CLAMP((float)key.time, 0.0f, (float)movie_duration(state));
    movie_apply_time(state, key.time, true);
}

// Moves the viewport back so that the whole camera path (the eye, what it looks at, and the keys) is in view, seen from the direction
// the viewport looks now
static void movie_scene_fit_path(ApplicationState* state) {
    auto& m = state->movie;
    const size_t n = md_array_size(m.keyframes);
    if (n == 0) return;
    const MovieResolved r = movie_resolve(state, m.keyframes, n);
    std::vector<vec3_t> pts;
    for (size_t i = 0; i < n; ++i) {
        pts.push_back(r.keys[i].transform.position);
        pts.push_back(r.look[i]);
        if (r.keys[i].focus_on) pts.push_back(r.focus[i]);
    }
    const double t0 = r.keys.front().time, t1 = r.keys.back().time;
    if (n >= 2 && t1 > t0) {
        const int samples = CLAMP((int)n * 8, 16, 256);
        for (int i = 0; i <= samples; ++i) {
            ViewTransform vt;
            float fov_y;
            const double t = t0 + (t1 - t0) * (double)i / (double)samples;
            movie_pose(state, r, t, &vt, &fov_y);
            pts.push_back(vt.position);
            pts.push_back(camera_get_look_at(vt));
        }
    }
    // Framed in the part of the viewport the Movie window leaves free (the picture is shifted there, see movie_frame_guide_shift),
    // as seen from the current direction: centred on the path's extent across the view, as close as its width and height allow
    ViewTransform t = state->view.target;
    const vec3_t right = quat_mul_vec3(t.orientation, vec3_set(1, 0, 0));
    const vec3_t up = quat_mul_vec3(t.orientation, vec3_set(0, 1, 0));
    const vec3_t back = quat_mul_vec3(t.orientation, vec3_set(0, 0, 1));
    float x0 = FLT_MAX, x1 = -FLT_MAX, y0 = FLT_MAX, y1 = -FLT_MAX, z0 = FLT_MAX, z1 = -FLT_MAX;
    for (const vec3_t& p : pts) {
        const float x = vec3_dot(p, right), y = vec3_dot(p, up), z = vec3_dot(p, back);
        x0 = MIN(x0, x); x1 = MAX(x1, x); y0 = MIN(y0, y); y1 = MAX(y1, y); z0 = MIN(z0, z); z1 = MAX(z1, z);
    }
    const vec3_t look = right * ((x0 + x1) * 0.5f) + up * ((y0 + y1) * 0.5f) + back * ((z0 + z1) * 0.5f);
    // The frame has the camera's field of view across its height (see movie_guide_fov_y)
    ImVec2 gp, gs;
    float aspect = (float)MAX(state->app.window.width, 1) / (float)MAX(state->app.window.height, 1);
    if (movie_frame_guide(state, &gp, &gs) && gs.y > 0.0f) aspect = gs.x / gs.y;
    const float tan_half = tanf(state->view.camera.fov_y * 0.5f) * 0.94f;
    const float tan_x = tan_half * aspect, tan_y = tan_half;
    float dist = 2.0f;
    for (const vec3_t& p : pts) {
        const vec3_t d = p - look;
        const float dz = vec3_dot(d, back) + 1.0f;
        dist = MAX(dist, dz + fabsf(vec3_dot(d, right)) / tan_x);
        dist = MAX(dist, dz + fabsf(vec3_dot(d, up)) / tan_y);
    }
    t.distance = dist;
    t.position = camera_position_from_look_at(look, t.orientation, t.distance);
    state->view.target = t;
}

// Changes between Scene view and Movie preview. Each keeps its own viewport pose: the first time Scene view is entered it frames the
// whole path, afterwards it is where it was left; Movie preview goes back to the movie camera.
static void movie_set_scene_view(ApplicationState* state, bool scene) {
    auto& m = state->movie;
    if (m.state == MovieRecordingState::Recording || scene == !m.show_frame) return;
    auto save = [&](decltype(m.pose_scene)& p) {
        p.valid = true;
        p.target = state->view.target;
        p.camera = state->view.camera;
        p.fov_y = state->view.camera.fov_y;
    };
    auto restore = [&](const decltype(m.pose_scene)& p) {
        state->view.target = p.target;
        state->view.camera.fov_y = p.fov_y;
    };
    if (scene) {
        save(m.pose_movie);
        m.show_frame = false;
        m.show_overlay_preview = false;
        if (m.pose_scene.valid) restore(m.pose_scene);
        else movie_scene_fit_path(state);
    } else {
        save(m.pose_scene);
        m.show_frame = true;
        m.show_overlay_preview = true;
        if (m.pose_movie.valid) restore(m.pose_movie);
        movie_apply_time(state, (double)m.playhead, true);
    }
}

// The Preview: the movie plays in the viewport through the movie camera, with the frame and the overlays, and every window is
// hidden but a small bar. Leaving it stops the playback and goes back to Scene view if that was shown.
static void movie_play_mode(ApplicationState* state, bool on) {
    auto& m = state->movie;
    if (on == m.play_mode) return;
    if (on) {
        if (m.state == MovieRecordingState::Recording || !m.show_window) return;
        m.play_mode_scene = !m.show_frame;
        movie_set_scene_view(state, false);
        // Watching the movie: from its start, whatever the preview time
        m.playhead = 0.0f;
        movie_apply_time(state, 0.0, true);
        m.preview_playing = true;
        m.play_mode = true;
        m.play_mode_moved_t = ImGui::GetTime();
    } else {
        m.play_mode = false;
        m.preview_playing = false;
        if (m.play_mode_scene) movie_set_scene_view(state, true);
    }
}

// The bar at the bottom of the viewport in the Preview: play and pause, back to the start, the time, Repeat and Exit. It hides while the
// mouse is still during playback.
static void draw_movie_play_bar(ApplicationState* state) {
    auto& m = state->movie;
    const float len = (float)movie_duration(state);
    const ImGuiIO& io = ImGui::GetIO();
    const double now = ImGui::GetTime();
    if (io.MouseDelta.x != 0.0f || io.MouseDelta.y != 0.0f || ImGui::IsAnyMouseDown() || !m.preview_playing) m.play_mode_moved_t = now;
    static bool hovered = false;
    if (!hovered && now - m.play_mode_moved_t > 2.5) return;

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    const float fs = ImGui::GetFontSize();
    const float w = CLAMP(vp->WorkSize.x * 0.5f, MIN(fs * 34.0f, vp->WorkSize.x), vp->WorkSize.x);
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x * 0.5f, vp->WorkPos.y + vp->WorkSize.y - fs * 1.2f), ImGuiCond_Always, ImVec2(0.5f, 1.0f));
    ImGui::SetNextWindowSize(ImVec2(w, 0.0f));
    ImGui::SetNextWindowBgAlpha(0.75f);
    const ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoDocking |
                                   ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_AlwaysAutoResize;
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, fs * 0.5f);
    if (ImGui::Begin("##movie_play_bar", nullptr, flags)) {
        hovered = ImGui::IsWindowHovered(ImGuiHoveredFlags_ChildWindows);
        const ImVec2 bs(ImGui::GetFrameHeight() * 1.3f, ImGui::GetFrameHeight());
        if (ImGui::Button(m.preview_playing ? (const char*)ICON_FA_PAUSE : (const char*)ICON_FA_PLAY, bs)) {
            m.preview_playing = !m.preview_playing;
            if (m.preview_playing && m.playhead >= len) m.playhead = 0.0f;
        }
        ImGui::SetItemTooltip(m.preview_playing ? "Pause (Space)" : "Play (Space)");
        ImGui::SameLine();
        if (ImGui::Button((const char*)ICON_FA_BACKWARD_STEP, bs)) {
            m.playhead = 0.0f;
            movie_apply_time(state, 0.0, true);
        }
        ImGui::SetItemTooltip("Back to the start");
        ImGui::SameLine();
        const float exit_w = ImGui::CalcTextSize("Exit preview").x + ImGui::GetStyle().FramePadding.x * 2.0f;
        const float repeat_w = ImGui::GetFrameHeight() + ImGui::CalcTextSize("Repeat").x + ImGui::GetStyle().ItemInnerSpacing.x;
        ImGui::SetNextItemWidth(MAX(ImGui::GetContentRegionAvail().x - exit_w - repeat_w - ImGui::GetStyle().ItemSpacing.x * 2.0f, fs * 6.0f));
        char fmt[48];
        snprintf(fmt, sizeof(fmt), "%%.2f / %.2f s", (double)len);
        if (ImGui::SliderFloat("##play_time", &m.playhead, 0.0f, len, fmt)) movie_apply_time(state, (double)m.playhead, true);
        if (ImGui::IsItemActivated()) {
            m.play_mode_scrub_resume = m.preview_playing;
            m.preview_playing = false;
        }
        if (ImGui::IsItemDeactivated() && m.play_mode_scrub_resume) {
            m.preview_playing = m.playhead < len || m.preview_loop;
            m.play_mode_scrub_resume = false;
        }
        ImGui::SameLine();
        ImGui::Checkbox("Repeat", &m.preview_loop);
        ImGui::SameLine();
        if (ImGui::Button("Exit preview")) movie_play_mode(state, false);
        ImGui::SetItemTooltip("Back to the editor (Esc)");
    } else {
        hovered = false;
    }
    ImGui::End();
    ImGui::PopStyleVar();
}

// The camera space (where the keys are) to the pixels of the viewport, for what is drawn over the viewport with the draw list
struct MoviePathView {
    mat4_t mvp = {};
    mat4_t inv = {};
    float  w = 1.0f, h = 1.0f;

    ImVec2 pixel(const vec4_t& c) const { return ImVec2((c.x / c.w * 0.5f + 0.5f) * w, (-c.y / c.w * 0.5f + 0.5f) * h); }

    // False when the point is behind the camera
    bool point(vec3_t p, ImVec2* out) const {
        const vec4_t c = mat4_mul_vec4(mvp, vec4_from_vec3(p, 1.0f));
        if (c.w <= 1.0e-3f) return false;
        *out = pixel(c);
        return true;
    }

    // A line between two points, cut where it goes behind the camera
    bool segment(vec3_t a, vec3_t b, ImVec2* pa, ImVec2* pb) const {
        vec4_t ca = mat4_mul_vec4(mvp, vec4_from_vec3(a, 1.0f));
        vec4_t cb = mat4_mul_vec4(mvp, vec4_from_vec3(b, 1.0f));
        if (!clip_segment_near(&ca, &cb)) return false;
        *pa = pixel(ca);
        *pb = pixel(cb);
        return true;
    }

    // The ray through a pixel
    void ray(ImVec2 px, vec3_t* origin, vec3_t* dir) const {
        const float nx = px.x / w * 2.0f - 1.0f, ny = -(px.y / h * 2.0f - 1.0f);
        const vec4_t n = mat4_mul_vec4(inv, vec4_set(nx, ny, -1.0f, 1.0f));
        const vec4_t f = mat4_mul_vec4(inv, vec4_set(nx, ny, 1.0f, 1.0f));
        const vec3_t a = vec3_set(n.x / n.w, n.y / n.w, n.z / n.w);
        const vec3_t b = vec3_set(f.x / f.w, f.y / f.w, f.z / f.w);
        *origin = a;
        *dir = vec3_normalize(b - a);
    }
};

static void movie_path_view(const ApplicationState* state, MoviePathView* v) {
    Camera cam = state->view.camera;
    ImVec2 guide_pos, guide_size;
    if (movie_frame_guide(state, &guide_pos, &guide_size) && guide_size.y > 0.0f) {
        cam.fov_y = movie_guide_fov_y(cam.fov_y, (float)state->app.window.height, guide_size.y);
    }
    const float aspect = (float)state->gbuffer.width / (float)MAX((int)state->gbuffer.height, 1);
    mat4_t P;
    if (state->view.mode == CameraMode::Perspective) {
        P = camera_view_to_clip_matrix_persp(cam, aspect);
    } else {
        const float h = cam.distance * tanf(cam.fov_y * 0.5f);
        const float w = aspect * h;
        P = camera_view_to_clip_matrix_ortho(-w, w, -h, h, cam.near_plane, cam.far_plane);
    }
    float sx, sy;
    if (movie_frame_guide_shift(state, &sx, &sy)) P = mat4_translate(sx, sy, 0.0f) * P;
    v->mvp = P * camera_world_to_view_matrix(cam);
    v->inv = mat4_inverse(v->mvp);
    v->w = (float)state->app.window.width;
    v->h = (float)state->app.window.height;
}

static bool movie_viewport_hovered() {
    const ImGuiWindow* win = ImGui::GetCurrentContext()->HoveredWindow;
    return win && strcmp(win->Name, "Main interaction window") == 0;
}

static ImU32 movie_path_fade(ImU32 c, float a) {
    return (c & 0x00FFFFFFu) | ((ImU32)((float)((c >> 24) & 0xFF) * a) << 24);
}

// A camera as a pyramid from the eye, with a roof on the top edge so that it can be told which way is up
static void movie_path_frustum(const MoviePathView& V, ImDrawList* dl, const ViewTransform& vt, float fov_y, float aspect, float length, ImU32 col, float thickness) {
    const vec3_t eye   = vt.position;
    const vec3_t fwd   = vt.orientation * vec3_t{0, 0, -1};
    const vec3_t right = vt.orientation * vec3_t{1, 0, 0};
    const vec3_t up    = vt.orientation * vec3_t{0, 1, 0};
    const float  hh = tanf(fov_y * 0.5f) * length;
    const float  hw = hh * aspect;
    const vec3_t c  = eye + fwd * length;
    const vec3_t p[4] = { c - right * hw - up * hh, c + right * hw - up * hh, c + right * hw + up * hh, c - right * hw + up * hh };
    auto line = [&](vec3_t a, vec3_t b) {
        ImVec2 pa, pb;
        if (V.segment(a, b, &pa, &pb)) dl->AddLine(pa, pb, col, thickness);
    };
    for (int i = 0; i < 4; ++i) {
        line(eye, p[i]);
        line(p[i], p[(i + 1) % 4]);
    }
    const vec3_t roof = c + up * (hh * 1.4f);
    line(p[2], roof);
    line(p[3], roof);
}

// The camera's paths in the viewport, drawn over it: the path of the eye (blue), of what it looks at (yellow) and of what is sharp
// (magenta, when a key sets the focus), brighter ahead of the preview time than behind it, with ticks at round times and chevrons
// that show the direction. Ticks that are close together are where the camera is slow. Each key has a handle on the eye, one on what
// it looks at and, when it sets the focus, one on what is sharp, with its number and name: click one to go to the key, drag it to
// edit that part of the key (the others stay; with Ctrl the key moves as a whole). Ctrl + click on a path adds a key there.
// The camera at the preview time is green, with the line it looks along. Where a key tracks atoms the paths are made with the atoms
// where they are at the trajectory frame shown.
static void movie_draw_camera_path(ApplicationState* state, ImDrawList* dl) {
    auto& m = state->movie;
    auto& d = m.path_drag;
    m.path_hot = false;
    m.path_hover_key = -1;
    const size_t n = md_array_size(m.keyframes);
    const bool shown = movie_path_visible(movie_scene_view(state), m.show_path, m.path_options) && n > 0 && m.show_window && m.state != MovieRecordingState::Recording && str_empty(state->screenshot.path_to_file);
    if (!shown) {
        d.active = false;
        m.lane_hover_key = -1;
        return;
    }

    MoviePathView V;
    movie_path_view(state, &V);
    const ImGuiIO& io = ImGui::GetIO();
    const ImVec2 mouse = io.MousePos;
    const bool over = movie_viewport_hovered() && ImGui::IsMousePosValid(&mouse);
    const bool picking = m.look_pick_key >= 0 || m.focus_pick_key >= 0;
    const float fs = ImGui::GetFontSize();
    const float handle_r = 9.0f;
    CameraKeyframe* keys = m.keyframes;
    const MovieResolved R = movie_resolve(state, keys, n);
    auto focus_world = [&](int i) {
        const CameraKeyframe& k = keys[i];
        if (k.focus_target == FocusTarget::LookAt) return R.look[(size_t)i];
        if (k.focus_target == FocusTarget::Distance) return k.transform.position + k.transform.orientation * vec3_t{0, 0, -k.focus_distance};
        return R.focus[(size_t)i];
    };
    auto handle_world = [&](int i, int kind) { return kind == 0 ? keys[i].transform.position : kind == 1 ? R.look[(size_t)i] : focus_world(i); };
    if (d.active && d.key >= (int)n) d.active = false;

    // ## Handles: hover, click and drag
    int hover_key = -1, hover_kind = 0;
    if (over && !d.active && !picking) {
        float best = (handle_r + 4.0f) * (handle_r + 4.0f);
        for (int i = 0; i < (int)n; ++i) {
            for (int kind = 0; kind < 3; ++kind) {
                if (kind == 2 && !keys[i].focus_on) continue;
                ImVec2 p;
                if (!V.point(handle_world(i, kind), &p)) continue;
                const float dx = p.x - mouse.x, dy = p.y - mouse.y;
                if (dx * dx + dy * dy < best) {
                    best = dx * dx + dy * dy;
                    hover_key = i;
                    hover_kind = kind;
                }
            }
        }
    }
    if (hover_key >= 0 && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
        d = decltype(m.path_drag){};
        if (!m.sel.contains(KeyKind::Camera, 0, keys[hover_key].time)) m.sel.set(KeyKind::Camera, 0, keys[hover_key].time);
        d.active = true;
        d.key = hover_key;
        d.kind = hover_kind;
        d.start = keys[hover_key];
        d.start_x = mouse.x;
        d.start_y = mouse.y;
        d.plane_point = handle_world(hover_key, hover_kind);
        d.plane_normal = state->view.camera.orientation * vec3_t{0, 0, -1};
        vec3_t o, dir, hit;
        V.ray(mouse, &o, &dir);
        d.grab = ray_plane_hit(o, dir, d.plane_point, d.plane_normal, &hit) ? hit - d.plane_point : vec3_t{0, 0, 0};
    }
    if (d.active) {
        if (ImGui::IsKeyPressed(ImGuiKey_Escape, false)) {
            if (d.moved) keys[d.key] = d.start;
            d.active = false;
        } else if (ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
            if (!d.moved && fabsf(mouse.x - d.start_x) + fabsf(mouse.y - d.start_y) > 4.0f) d.moved = true;
            if (d.moved) {
                vec3_t o, dir, hit;
                V.ray(mouse, &o, &dir);
                if (ray_plane_hit(o, dir, d.plane_point, d.plane_normal, &hit)) {
                    const vec3_t pos = hit - d.grab;
                    CameraKeyframe k = d.start;
                    if (io.KeyCtrl)         camera_key_translate(&k, pos - d.plane_point);
                    else if (d.kind == 0)   camera_key_set_eye(&k, pos);
                    else if (d.kind == 1) {
                        // A tracked look-at that is dragged becomes the fixed point it is dropped on
                        if (camera_key_set_look(&k, pos)) k.look_set = 0;
                    } else {
                        k.focus_target = FocusTarget::Point;
                        k.focus_set = 0;
                        k.focus_point = pos;
                    }
                    keys[d.key] = k;
                }
            }
        } else {
            if (!d.moved) movie_goto_keyframe(state, (size_t)d.key);
            d.active = false;
        }
    }
    if (hover_key >= 0 || d.active) {
        m.path_hot = true;
        ImGui::SetMouseCursor(ImGuiMouseCursor_Hand);
    }
    m.path_hover_key = d.active ? d.key : hover_key;

    // ## The paths themselves, made again when the keys or what they track change
    const std::vector<CameraKeyframe>& sorted = R.keys;
    const bool have_focus = movie_focus_keys_exist(sorted.data(), n);
    static uint64_t path_signature = 0;
    static CameraPathSamples path;
    static std::vector<vec3_t> path_focus;
    {
        uint64_t h = md_hash64(sorted.data(), n * sizeof(CameraKeyframe), 21);
        h = md_hash64_combine(h, md_hash64(R.look.data(), n * sizeof(vec3_t), 22));
        h = md_hash64_combine(h, md_hash64(R.focus.data(), n * sizeof(vec3_t), 23));
        const double scalars[] = {(double)m.up_axis, (double)m.keep_upright};
        h = md_hash64_combine(h, md_hash64(scalars, sizeof(scalars), 24));
        if (h != path_signature || path.time.empty()) {
            path_signature = h;
            path = {};
            path_focus.clear();
            const double t0 = sorted.front().time, t1 = sorted.back().time;
            const int samples = n >= 2 && t1 > t0 ? CLAMP((int)n * 24, 64, 400) : 0;
            for (int i = 0; i <= samples; ++i) {
                ViewTransform vt;
                float fov_y;
                const double t = samples > 0 ? t0 + (t1 - t0) * (double)i / (double)samples : t0;
                movie_pose(state, R, t, &vt, &fov_y);
                path.time.push_back(t);
                path.eye.push_back(vt.position);
                path.look.push_back(camera_get_look_at(vt));
                path_focus.push_back(have_focus ? movie_focus_evaluate(R, t, vt) : path.look.back());
            }
        }
    }

    const ImU32 col_eye  = IM_COL32(90, 200, 255, 255);
    const ImU32 col_look = IM_COL32(255, 200, 60, 255);
    const ImU32 col_focus = IM_COL32(255, 90, 220, 255);
    const ImU32 col_head = IM_COL32(80, 255, 120, 255);
    const ImU32 col_dark = IM_COL32(15, 15, 20, 255);
    const ImU32 col_spin = IM_COL32(200, 130, 250, 255);
    const double playhead = (double)m.playhead;
    const float aspect = [&] { int w = 0, h = 0; movie_frame_size(state, &w, &h); return h > 0 ? (float)w / (float)h : 1.0f; }();

    auto label = [&](ImVec2 p, const char* text, ImU32 col) {
        const ImVec2 size = ImGui::CalcTextSize(text);
        dl->AddRectFilled(ImVec2(p.x - 3.0f, p.y - 1.0f), ImVec2(p.x + size.x + 3.0f, p.y + size.y + 1.0f), IM_COL32(0, 0, 0, 150), 3.0f);
        dl->AddText(p, col, text);
    };

    const size_t samples = path.time.size();
    std::vector<vec2_t> screen_eye(samples), screen_look(samples), screen_focus(samples);
    std::vector<char> ok_eye(samples, 0), ok_look(samples, 0), ok_focus(samples, 0);
    for (size_t i = 0; i < samples; ++i) {
        ImVec2 p;
        if (V.point(path.eye[i], &p))  { screen_eye[i] = vec2_set(p.x, p.y);  ok_eye[i] = 1; }
        if (V.point(path.look[i], &p)) { screen_look[i] = vec2_set(p.x, p.y); ok_look[i] = 1; }
        if (have_focus && V.point(path_focus[i], &p)) { screen_focus[i] = vec2_set(p.x, p.y); ok_focus[i] = 1; }
    }

    // Round times along the path
    const double t_first = path.time.front(), t_last = path.time.back();
    const double step = camera_tick_step(t_last - t_first, 30);
    std::vector<double> ticks;
    if (t_last > t_first) {
        for (int64_t k = (int64_t)ceil(t_first / step - 1.0e-9); (double)k * step <= t_last + 1.0e-9; ++k) ticks.push_back((double)k * step);
    }

    // What the camera looks at from the eye, at each tick
    if ((m.path_options & 2) && !ticks.empty()) {
        for (double t : ticks) {
            vec3_t e, l;
            camera_path_at(path, t, &e, &l);
            ImVec2 a, b;
            if (V.segment(e, l, &a, &b)) dl->AddLine(a, b, IM_COL32(255, 255, 255, t < playhead ? 18 : 40), 1.0f);
        }
    }

    // The paths, dimmer before the preview time
    auto draw_path = [&](const std::vector<vec3_t>& pts, ImU32 col, float thickness) {
        for (size_t i = 1; i < samples; ++i) {
            ImVec2 a, b;
            if (!V.segment(pts[i - 1], pts[i], &a, &b)) continue;
            const double mid = 0.5 * (path.time[i - 1] + path.time[i]);
            dl->AddLine(a, b, mid < playhead ? movie_path_fade(col, 0.4f) : col, thickness);
        }
    };
    if (have_focus) draw_path(path_focus, col_focus, 2.0f);
    draw_path(path.look, col_look, 2.0f);
    draw_path(path.eye, col_eye, 3.0f);

    // Ticks, labels and chevrons
    if ((m.path_options & 1) && !ticks.empty()) {
        auto chevron = [&](ImVec2 p, ImVec2 q, ImU32 col) {
            ImVec2 dir(q.x - p.x, q.y - p.y);
            const float len = sqrtf(dir.x * dir.x + dir.y * dir.y);
            if (len < 1.0e-3f) return;
            dir = ImVec2(dir.x / len, dir.y / len);
            const ImVec2 nrm(-dir.y, dir.x);
            const ImVec2 back(p.x - dir.x * 7.0f, p.y - dir.y * 7.0f);
            dl->AddLine(p, ImVec2(back.x + nrm.x * 5.0f, back.y + nrm.y * 5.0f), col, 2.0f);
            dl->AddLine(p, ImVec2(back.x - nrm.x * 5.0f, back.y - nrm.y * 5.0f), col, 2.0f);
        };
        ImVec2 last_label(-1.0e9f, -1.0e9f), prev_tick(-1.0e9f, -1.0e9f);
        double prev_t = 0.0;
        bool have_prev = false;
        for (double t : ticks) {
            vec3_t e, l;
            camera_path_at(path, t, &e, &l);
            const ImU32 col = t < playhead ? movie_path_fade(col_eye, 0.5f) : col_eye;
            ImVec2 pe;
            if (V.point(e, &pe)) {
                dl->AddCircleFilled(pe, 4.0f, col_dark);
                dl->AddCircleFilled(pe, 3.0f, col);
                const float dx = pe.x - last_label.x, dy = pe.y - last_label.y;
                if (dx * dx + dy * dy > 60.0f * 60.0f) {
                    char buf[24];
                    snprintf(buf, sizeof(buf), "%g s", t);
                    label(ImVec2(pe.x + 7.0f, pe.y - fs - 3.0f), buf, col_eye);
                    last_label = pe;
                }
                // A chevron half way to the previous tick, pointing the way the camera goes
                if (have_prev) {
                    const double tm = 0.5 * (t + prev_t);
                    vec3_t em, lm, en, ln;
                    camera_path_at(path, tm, &em, &lm);
                    camera_path_at(path, tm + 0.02 * step, &en, &ln);
                    ImVec2 pm, pn;
                    const float dxp = pe.x - prev_tick.x, dyp = pe.y - prev_tick.y;
                    if (dxp * dxp + dyp * dyp > 26.0f * 26.0f && V.point(em, &pm) && V.point(en, &pn)) chevron(pm, pn, col);
                }
                prev_tick = pe;
                prev_t = t;
                have_prev = true;
            } else {
                have_prev = false;
            }
            ImVec2 pl;
            if (V.point(l, &pl)) {
                const ImU32 lc = t < playhead ? movie_path_fade(col_look, 0.5f) : col_look;
                dl->AddQuadFilled(ImVec2(pl.x, pl.y - 4.0f), ImVec2(pl.x + 4.0f, pl.y), ImVec2(pl.x, pl.y + 4.0f), ImVec2(pl.x - 4.0f, pl.y), lc);
            }
        }
    }

    // The ring that the camera goes round in a spin, turning the way the turns go
    if (m.path_options & 8) {
        for (const CameraBand& b : camera_bands(sorted)) {
            if (b.kind != CameraBandKind::Spin) continue;
            vec3_t c, u, v;
            float r;
            if (!camera_spin_ring(sorted[(size_t)b.first], sorted[(size_t)b.last], &c, &u, &v, &r)) continue;
            auto on_ring = [&](float a) { return c + (u * cosf(a) + v * sinf(a)) * r; };
            const int segments = 72;
            for (int s = 0; s < segments; ++s) {
                const float a0 = 6.2831853f * (float)s / (float)segments, a1 = 6.2831853f * (float)(s + 1) / (float)segments;
                ImVec2 pa, pb;
                if (V.segment(on_ring(a0), on_ring(a1), &pa, &pb)) dl->AddLine(pa, pb, movie_path_fade(col_spin, 0.55f), 2.0f);
            }
            const float dir = b.turns > 0 ? 1.0f : -1.0f;
            for (int k = 0; k < 4; ++k) {
                const float a = 0.7853982f + 1.5707963f * (float)k;
                ImVec2 pa, pb;
                if (V.point(on_ring(a), &pa) && V.point(on_ring(a + dir * 0.05f), &pb)) {
                    ImVec2 dd(pb.x - pa.x, pb.y - pa.y);
                    const float len = sqrtf(dd.x * dd.x + dd.y * dd.y);
                    if (len < 1.0e-3f) continue;
                    dd = ImVec2(dd.x / len, dd.y / len);
                    const ImVec2 nrm(-dd.y, dd.x);
                    const ImVec2 back(pa.x - dd.x * 8.0f, pa.y - dd.y * 8.0f);
                    dl->AddLine(pa, ImVec2(back.x + nrm.x * 5.0f, back.y + nrm.y * 5.0f), col_spin, 2.5f);
                    dl->AddLine(pa, ImVec2(back.x - nrm.x * 5.0f, back.y - nrm.y * 5.0f), col_spin, 2.5f);
                }
            }
            char buf[32];
            snprintf(buf, sizeof(buf), "%+d x", b.turns);
            ImVec2 p;
            if (V.point(on_ring(0.0f), &p)) label(ImVec2(p.x + 8.0f, p.y + 4.0f), buf, col_spin);
        }
    }

    // The cameras at the keys
    if (m.path_options & 4) {
        for (size_t i = 0; i < n; ++i) {
            const CameraKeyframe& k = keys[i];
            const bool at_head = fabs(k.time - playhead) < 1.0e-3;
            ViewTransform kt = k.transform;
            if (const vec3_t* up = movie_upright(state)) camera_level(&kt, *up, k.roll);
            movie_path_frustum(V, dl, kt, k.fov_y, aspect, k.transform.distance * 0.25f, at_head ? IM_COL32(255, 140, 50, 255) : IM_COL32(255, 255, 255, 190), at_head ? 2.0f : 1.0f);
        }
    }

    // The keys: a handle on the eye and one on what it looks at, joined, with the number and the name, and one on what is sharp
    for (int i = 0; i < (int)n; ++i) {
        const CameraKeyframe& k = keys[i];
        const bool at_head = fabs(k.time - playhead) < 1.0e-3;
        const bool hot = i == m.path_hover_key || i == m.lane_hover_key;
        ImU32 fill = IM_COL32(235, 235, 235, 255);
        if (k.look_set != 0) fill = IM_COL32(110, 150, 255, 255);
        if (at_head) fill = IM_COL32(255, 140, 50, 255);

        const vec3_t eye = k.transform.position, look = R.look[(size_t)i];
        ImVec2 a, b;
        if (V.segment(eye, look, &a, &b)) dl->AddLine(a, b, IM_COL32(255, 255, 255, hot ? 200 : 90), hot ? 2.0f : 1.0f);

        char num[16];
        snprintf(num, sizeof(num), "%d", i + 1);
        const ImVec2 num_size = ImGui::CalcTextSize(num);
        ImVec2 pe, pl;
        if (V.point(eye, &pe)) {
            if (hot) dl->AddCircle(pe, handle_r + 4.0f, IM_COL32(255, 235, 90, 255), 0, 2.5f);
            if (m.sel.contains(KeyKind::Camera, 0, k.time)) dl->AddCircle(pe, handle_r + 2.5f, IM_COL32(255, 255, 255, 255), 0, 2.0f);
            dl->AddCircleFilled(pe, handle_r, fill);
            dl->AddCircle(pe, handle_r, col_dark, 0, 1.5f);
            dl->AddText(ImVec2(pe.x - num_size.x * 0.5f, pe.y - num_size.y * 0.5f), col_dark, num);
            if (k.name[0] != '\0') label(ImVec2(pe.x + handle_r + 5.0f, pe.y - fs * 0.5f), k.name, IM_COL32(255, 255, 255, 255));
        }
        if (V.point(look, &pl)) {
            const float r = handle_r - 2.0f;
            if (hot) dl->AddCircle(pl, r + 4.0f, IM_COL32(255, 235, 90, 255), 0, 2.5f);
            dl->AddCircleFilled(pl, r, IM_COL32(20, 20, 20, 140));
            dl->AddCircle(pl, r, fill, 0, 2.5f);
            dl->AddLine(ImVec2(pl.x - r, pl.y), ImVec2(pl.x + r, pl.y), fill, 1.5f);
            dl->AddLine(ImVec2(pl.x, pl.y - r), ImVec2(pl.x, pl.y + r), fill, 1.5f);
            label(ImVec2(pl.x + r + 3.0f, pl.y - fs - 2.0f), num, col_look);
        }
        if (k.focus_on) {
            // What is sharp from here on: a diamond, joined to the eye; at the look-at point when it is the look-at that is sharp
            ImVec2 pf;
            if (V.point(focus_world(i), &pf)) {
                const float r = handle_r - 1.0f;
                if (V.point(eye, &pe)) dl->AddLine(pe, pf, movie_path_fade(col_focus, hot ? 0.8f : 0.3f), 1.0f);
                if (hot) dl->AddCircle(pf, r + 5.0f, IM_COL32(255, 235, 90, 255), 0, 2.5f);
                const ImVec2 dia[4] = {ImVec2(pf.x, pf.y - r), ImVec2(pf.x + r, pf.y), ImVec2(pf.x, pf.y + r), ImVec2(pf.x - r, pf.y)};
                dl->AddConvexPolyFilled(dia, 4, IM_COL32(20, 20, 20, 160));
                dl->AddPolyline(dia, 4, col_focus, ImDrawFlags_Closed, 2.5f);
                label(ImVec2(pf.x + r + 3.0f, pf.y + 2.0f), num, col_focus);
            }
        }
    }

    // The camera at the preview time: where it is, and where it looks
    ViewTransform head;
    float head_fov;
    movie_pose(state, R, playhead, &head, &head_fov);
    {
        const vec3_t eye = head.position, look = camera_get_look_at(head);
        movie_path_frustum(V, dl, head, head_fov, aspect, head.distance * 0.18f, col_head, 2.5f);
        ImVec2 a, b;
        if (V.segment(eye, look, &a, &b)) dl->AddLine(a, b, col_head, 3.0f);
        ImVec2 pe, pl;
        const bool have_eye = V.point(eye, &pe);
        if (have_eye) {
            dl->AddCircleFilled(pe, 10.0f, movie_path_fade(col_head, 0.85f));
            dl->AddCircle(pe, 10.0f, col_dark, 0, 1.5f);
            char buf[32];
            snprintf(buf, sizeof(buf), "camera %.2f s", m.playhead);
            label(ImVec2(pe.x + 14.0f, pe.y - fs * 0.5f), buf, col_head);
        }
        if (V.point(look, &pl)) {
            // The arrow head on the line of sight, at the target
            ImVec2 dir(0.0f, 1.0f);
            if (have_eye) dir = ImVec2(pl.x - pe.x, pl.y - pe.y);
            else dir = ImVec2(b.x - a.x, b.y - a.y);
            const float len = sqrtf(dir.x * dir.x + dir.y * dir.y);
            if (len > 1.0e-3f) {
                dir = ImVec2(dir.x / len, dir.y / len);
                const ImVec2 nrm(-dir.y, dir.x);
                const ImVec2 tip(pl.x - dir.x * 14.0f, pl.y - dir.y * 14.0f);
                dl->AddTriangleFilled(tip, ImVec2(tip.x - dir.x * 14.0f + nrm.x * 7.0f, tip.y - dir.y * 14.0f + nrm.y * 7.0f), ImVec2(tip.x - dir.x * 14.0f - nrm.x * 7.0f, tip.y - dir.y * 14.0f - nrm.y * 7.0f), col_head);
            }
            dl->AddCircle(pl, 13.0f, col_head, 0, 3.0f);
            dl->AddLine(ImVec2(pl.x - 20.0f, pl.y), ImVec2(pl.x + 20.0f, pl.y), col_head, 2.0f);
            dl->AddLine(ImVec2(pl.x, pl.y - 20.0f), ImVec2(pl.x, pl.y + 20.0f), col_head, 2.0f);
            label(ImVec2(pl.x + 17.0f, pl.y - fs * 0.5f - 14.0f), "looks at", col_head);
        }

        if (have_focus || state->visuals.dof.enabled) {
            // Where it is sharp: a frame the size of the view at that depth, and a cross where the focus is
            const float depth = dof_focus_depth(state, head);
            const vec3_t fwd   = head.orientation * vec3_t{0, 0, -1};
            const vec3_t right = head.orientation * vec3_t{1, 0, 0};
            const vec3_t up    = head.orientation * vec3_t{0, 1, 0};
            const vec3_t c = head.position + fwd * depth;
            const float hh = depth * tanf(head_fov * 0.5f);
            const float hw = hh * aspect;
            const vec3_t p[4] = { c - right * hw - up * hh, c + right * hw - up * hh, c + right * hw + up * hh, c - right * hw + up * hh };
            for (int i = 0; i < 4; ++i) {
                ImVec2 pa, pb;
                if (V.segment(p[i], p[(i + 1) % 4], &pa, &pb)) dl->AddLine(pa, pb, col_focus, 1.5f);
            }
            const float s = hh * 0.1f;
            ImVec2 pa, pb;
            if (V.segment(c - right * s, c + right * s, &pa, &pb)) dl->AddLine(pa, pb, col_focus, 1.5f);
            if (V.segment(c - up * s, c + up * s, &pa, &pb)) dl->AddLine(pa, pb, col_focus, 1.5f);
        }
    }

    // ## Ctrl + click on a path adds a key there
    double insert_at = -1.0;
    int insert_path = 0;
    if (over && io.KeyCtrl && hover_key < 0 && !d.active && !picking) {
        int seg[3] = {-1, -1, -1};
        float u[3] = {}, dist[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
        const vec2_t q = vec2_set(mouse.x, mouse.y);
        dist[0] = polyline_nearest(screen_eye, ok_eye, q, &seg[0], &u[0]);
        dist[1] = polyline_nearest(screen_look, ok_look, q, &seg[1], &u[1]);
        if (have_focus) dist[2] = polyline_nearest(screen_focus, ok_focus, q, &seg[2], &u[2]);
        int best = 0;
        for (int k = 1; k < 3; ++k) if (dist[k] < dist[best]) best = k;
        if (dist[best] < 10.0f && seg[best] >= 0) {
            const double t = movie_snap_time(state, path.time[(size_t)seg[best]] + (path.time[(size_t)seg[best] + 1] - path.time[(size_t)seg[best]]) * (double)u[best]);
            const std::vector<vec3_t>& pts = best == 0 ? path.eye : best == 1 ? path.look : path_focus;
            const vec3_t at = pts[(size_t)seg[best]] + (pts[(size_t)seg[best] + 1] - pts[(size_t)seg[best]]) * u[best];
            ImVec2 p;
            if (V.point(at, &p)) {
                m.path_hot = true;
                dl->AddCircle(p, 11.0f, IM_COL32(255, 255, 255, 255), 0, 2.5f);
                dl->AddLine(ImVec2(p.x - 6.0f, p.y), ImVec2(p.x + 6.0f, p.y), IM_COL32(255, 255, 255, 255), 2.0f);
                dl->AddLine(ImVec2(p.x, p.y - 6.0f), ImVec2(p.x, p.y + 6.0f), IM_COL32(255, 255, 255, 255), 2.0f);
                char buf[48];
                snprintf(buf, sizeof(buf), "add a key at %.2f s", t);
                label(ImVec2(p.x + 15.0f, p.y - fs * 0.5f), buf, IM_COL32(255, 255, 255, 255));
                if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) { insert_at = t; insert_path = best; }
            }
        }
    }

    if (hover_key >= 0 && !d.active) {
        const CameraKeyframe& k = keys[hover_key];
        ImGui::SetTooltip("Keyframe %s at %.2f s, field of view %.1f deg\n%s\nClick to go to it. Ctrl + drag moves the whole key. Esc cancels a drag.",
            camera_key_label(k, hover_key).c_str(), k.time, k.fov_y * MOVIE_RAD_TO_DEG,
            hover_kind == 0 ? "Drag to move the eye: the camera keeps looking at the same point." :
            hover_kind == 1 ? "Drag to move what it looks at: the eye stays. A tracked look-at becomes a fixed point." :
            "Drag to move what is sharp from this key on. It becomes a fixed point.");
    }

    if (insert_at >= 0.0) {
        bool taken = false;
        for (size_t i = 0; i < n; ++i) taken |= fabs(keys[i].time - insert_at) < 1.0e-3;
        if (!taken) {
            ViewTransform vt;
            float fov_y;
            movie_pose(state, R, insert_at, &vt, &fov_y);
            CameraKeyframe key = camera_key_on_path(vt, fov_y, sorted, insert_at, movie_upright(state));
            if (insert_path == 2) {
                // On the path of what is sharp: the key sets the focus where the path is, so that it can be dragged from there
                key.focus_on = true;
                key.focus_target = FocusTarget::Point;
                key.focus_point = movie_focus_evaluate(R, insert_at, vt);
                key.focus_transition = 0.0f;
            }
            movie_insert_key(state, key, false);
            m.playhead = (float)insert_at;
        }
    }
    m.lane_hover_key = -1;
}

// Keyframes in the Timelines plots, which are in trajectory time. A keyframe is at a time of the movie,
// so it is put where the trajectory is at that time; keys that fall on the same spot (the trajectory is
// held) are stacked. Read only: they are edited in the Movie window. Returns whether one is hovered.
static bool movie_draw_timeline_markers(ApplicationState* data) {
    auto& m = data->movie;
    const size_t n = md_array_size(m.keyframes);
    if ((n == 0 && !m.show_window) || md_array_size(data->timeline.x_values) == 0) return false;

    ImDrawList* dl = ImPlot::GetPlotDrawList();
    const ImVec2 plot_pos  = ImPlot::GetPlotPos();
    const ImVec2 plot_size = ImPlot::GetPlotSize();
    const float  r = ImGui::GetFontSize() * 0.4f;
    const ImU32  col = IM_COL32(255, 120, 40, 255);

    int hovered_key = -1;
    const ImVec2 mouse = ImGui::GetMousePos();
    const bool plot_hovered = ImPlot::IsPlotHovered();

    ImPlot::PushPlotClipRect();

    if (m.show_window) {
        // The part of the trajectory that the movie plays
        const double f0 = MIN(m.start_frame, m.end_frame);
        const double f1 = MAX(m.start_frame, m.end_frame);
        const float x0 = ImPlot::PlotToPixels(frame_to_time(f0, *data), 0.0).x;
        const float x1 = ImPlot::PlotToPixels(frame_to_time(f1, *data), 0.0).x;
        dl->AddRectFilled(ImVec2(x0, plot_pos.y), ImVec2(MAX(x1, x0 + 1.0f), plot_pos.y + plot_size.y), IM_COL32(255, 170, 40, 26));
    }

    std::vector<float> xs(n);
    for (size_t i = 0; i < n; ++i) {
        const double t = frame_to_time(movie_trajectory_frame(data, m.keyframes[i].time), *data);
        const float x = ImPlot::PlotToPixels(t, 0.0).x;
        xs[i] = x;

        int stack = 0;
        for (size_t j = 0; j < i; ++j) {
            if (fabsf(xs[j] - x) < 2.0f * r) stack += 1;
        }

        const ImVec2 c(x, plot_pos.y + r + 2.0f + (float)stack * (2.0f * r + 2.0f));
        dl->AddLine(ImVec2(x, plot_pos.y), ImVec2(x, plot_pos.y + plot_size.y), IM_COL32(255, 120, 40, 90));
        dl->AddQuadFilled(ImVec2(c.x, c.y - r), ImVec2(c.x + r, c.y), ImVec2(c.x, c.y + r), ImVec2(c.x - r, c.y), col);

        char label[16];
        snprintf(label, sizeof(label), "%d", (int)i + 1);
        dl->AddText(ImVec2(c.x + r + 2.0f, c.y - ImGui::GetFontSize() * 0.5f), IM_COL32(255, 200, 150, 255), label);

        if (plot_hovered) {
            const float dx = mouse.x - c.x;
            const float dy = mouse.y - c.y;
            if (dx * dx + dy * dy < (r * 1.6f) * (r * 1.6f)) hovered_key = (int)i;
        }
    }

    ImPlot::PopPlotClipRect();

    if (hovered_key >= 0) {
        const CameraKeyframe& k = m.keyframes[hovered_key];
        ImGui::SetTooltip("Camera keyframe %d\nMovie time: %.2f s\nField of view: %.1f deg\nClick to go to it",
            hovered_key + 1, k.time, k.fov_y * MOVIE_RAD_TO_DEG);
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && m.state != MovieRecordingState::Recording) {
            movie_goto_keyframe(data, (size_t)hovered_key);
        }
    }
    return hovered_key >= 0;
}

static void draw_movie_param_lane(ApplicationState* data, float movie_len, bool locked);
static void draw_movie_rep_lane(ApplicationState* data, float movie_len, bool locked);
static void draw_movie_overlay_lane(ApplicationState* data, float movie_len, bool locked);
static void draw_movie_rep_overview_lane(ApplicationState* data, float movie_len, bool locked);
// What the camera lane asks for, done by the caller once the lanes are drawn: the keys are read while they are
struct MovieCameraLaneEdit {
    int    remove = -1;
    double insert = -1.0;
    bool   held = false;
};
static void draw_movie_camera_lane(ApplicationState* data, float movie_len, bool locked, const std::vector<CameraKeyframe>& sorted, MovieCameraLaneEdit* edit);
static void draw_movie_ruler(ApplicationState* data, float movie_len, bool locked);
static void movie_lane_title(const char* text);
static void movie_reps_apply(ApplicationState* state, double time);
static std::vector<RepRow> movie_rep_overview_rows(const ApplicationState* data);
static std::vector<RepBlock> movie_rep_system_blocks(const ApplicationState* data, const std::string& group, double duration);
static float movie_rep_overview_height(const ApplicationState* data, double duration);

// ## Picking keys and moving them together

struct MovieKeyPointResult {
    bool clicked = false, hovered = false, held = false;
};

// A click on a key, a bar or a block: it is picked (Ctrl adds it to the picked ones or takes it away); a click on one that is picked
// leaves the group for when the mouse is released without a drag. A click is reported when the mouse is released over the item, which
// a drag ends with too, so after a drag there is nothing to do.
static void movie_pick_click(ApplicationState* data, KeyKind kind, int64_t subject, double time, double end, const KeyShift& lane) {
    auto& m = data->movie;
    auto& kd = m.key_drag;
    if (kd.active && kd.moved) return;
    if (ImGui::GetIO().KeyCtrl) {
        m.sel.toggle(kind, subject, time, end);
    } else if (!m.sel.contains(kind, subject, time)) {
        m.sel.set(kind, subject, time, end);
    } else {
        kd.collapse = true;
        kd.collapse_to = {kind, subject, time, end};
    }
    m.sel_lane = lane;
}

// A drag of a key, a bar or a block to 'new_time' (its start; for a key with a value to 'py'): every picked item moves by the same,
// made from the keys as they were when the drag began, so nothing is lost on the way. A key that is not picked is picked alone.
static void movie_group_drag(ApplicationState* data, KeyKind kind, int64_t subject, double time, double end, double y, double new_time, double py,
    const KeyShift& lane, bool toggles_value = false) {
    auto& m = data->movie;
    auto& kd = m.key_drag;
    if (!m.sel.contains(kind, subject, time)) m.sel.set(kind, subject, time, end);
    m.sel_lane = lane;
    if (!kd.active) {
        kd.active = true;
        kd.moved = false;
        kd.start = movie_keys_snapshot(data);
        if (lane.lane == KeyLane::Distance) {
            const auto resolved = movie_resolve(data, m.keyframes, md_array_size(m.keyframes));
            for (auto& key : kd.start.camera) {
                if (!m.sel.contains(KeyKind::Camera, 0, key.time)) continue;
                for (size_t i = 0; i < resolved.keys.size(); ++i) if (resolved.keys[i].time == key.time) {
                    camera_aim_at(&key.transform, resolved.look[i]);
                    break;
                }
            }
        }
        kd.start_sel = m.sel;
        kd.x0 = time;
        kd.y0 = y;
    }
    KeyShift s = lane;
    if (lane.lane != KeyLane::None) s.dy = lane.ratio ? (kd.y0 != 0.0 ? py / kd.y0 : 1.0) : py - kd.y0;
    MovieKeys out = kd.start;
    kd.dt = movie_keys_shift(&out, &m.sel, kd.start, kd.start_sel, new_time - kd.x0, s, (double)movie_duration(data));
    kd.dy = s.dy;
    kd.moved = true;
    movie_keys_restore(data, out);
    if (kind == KeyKind::Block) movie_reps_apply(data, (double)m.playhead);
    if (toggles_value && m.sel.size() == 1) {
        // A key of Visible has no value to move, it is shown or hidden by where it is dragged to
        for (RepKey& k : m.rep_keys) {
            if (kind == KeyKind::Rep && rep_key_subject(k.rep, k.prop) == subject && fabs(k.time - m.sel.ids[0].time) < 1.0e-9) k.value[0] = py >= 0.5 ? 1.0f : 0.0f;
        }
    }
}

// What a drag is doing, as a tooltip next to the mouse
static void movie_group_drag_tooltip(ApplicationState* data, const KeyShift& lane) {
    const auto& kd = data->movie.key_drag;
    char buf[96];
    const int count = (int)data->movie.sel.size();
    const int len = snprintf(buf, sizeof(buf), "%d item%s: %+.2f s", count, count == 1 ? "" : "s", kd.dt);
    if (lane.lane != KeyLane::None) snprintf(buf + len, sizeof(buf) - (size_t)len, lane.ratio ? "  x%.3g" : "  %+.3g", kd.dy);
    ImGui::SetTooltip("%s", buf);
}

// One key as a point of a lane. A click picks it, dragging moves every picked item together: sideways in time and, in a lane with a
// value, up and down. Where a key ends on another it replaces it when the mouse is released.
static MovieKeyPointResult movie_key_point(ApplicationState* data, bool locked, KeyKind kind, int64_t subject, double time, double y, int id,
    ImVec4 color, float size, const KeyShift& lane, bool* any_held, bool toggles_value = false) {
    auto& m = data->movie;
    auto& kd = m.key_drag;
    const double duration = (double)movie_duration(data);
    const ImPlotDragToolFlags flags = ImPlotDragToolFlags_NoFit | (locked ? ImPlotDragToolFlags_NoInputs : 0);
    MovieKeyPointResult r;
    double x = time, py = y;
    const bool changed = ImPlot::DragPoint(id, &x, &py, color, size, flags, &r.clicked, &r.hovered, &r.held);

    if (!locked) {
        if (r.clicked) movie_pick_click(data, kind, subject, time, 0.0, lane);
        if (changed) movie_group_drag(data, kind, subject, time, 0.0, y, movie_snap_time(data, x), py, lane, toggles_value);
        if (r.held) {
            *any_held = true;
            m.playhead = (float)CLAMP(kd.active && kd.moved ? kd.x0 + kd.dt : time, 0.0, duration);
        }
    }

    const bool dragged = r.held && kd.active && kd.moved;
    if (dragged) movie_group_drag_tooltip(data, lane);

    // A ring on the picked keys
    const double ring_time = dragged ? kd.x0 + kd.dt : time;
    if (m.sel.contains(kind, subject, ring_time)) {
        ImPlot::PushPlotClipRect();
        ImPlot::GetPlotDrawList()->AddCircle(ImPlot::PlotToPixels(ring_time, dragged ? py : y), size + 3.5f, IM_COL32(255, 255, 255, 255), 0, 2.0f);
        ImPlot::PopPlotClipRect();
    }
    return r;
}

// Dragging on the empty background of a lane draws a box and picks the keys in it (Ctrl or Shift: in addition to the picked ones);
// a click on the background puts them all down. 'vlines' are times of lines in the lane that are dragged themselves.
template <typename Pick>
static void movie_lane_box(ApplicationState* data, bool locked, int lane_id, bool over_key, const double* vlines, int num_vlines, Pick pick) {
    auto& b = data->movie.key_box;
    const ImGuiIO& io = ImGui::GetIO();
    if (locked) {
        b.active = false;
        return;
    }
    if (!b.active && !over_key && !data->movie.key_drag.active && ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
        bool on_line = false;
        for (int i = 0; i < num_vlines; ++i) on_line |= fabsf(ImPlot::PlotToPixels(vlines[i], 0.0).x - io.MousePos.x) < 6.0f;
        if (!on_line) {
            const ImPlotPoint p = ImPlot::GetPlotMousePos();
            b = {};
            b.active = true;
            b.lane = lane_id;
            b.x0 = p.x;
            b.y0 = p.y;
            b.px = io.MousePos.x;
            b.py = io.MousePos.y;
            b.add = io.KeyCtrl || io.KeyShift;
        }
    }
    if (!b.active || b.lane != lane_id) return;
    const ImPlotPoint p = ImPlot::GetPlotMousePos();
    if (io.MouseDown[0]) {
        if (fabsf(io.MousePos.x - b.px) + fabsf(io.MousePos.y - b.py) > 4.0f) b.moved = true;
        if (b.moved) {
            ImPlot::PushPlotClipRect();
            ImDrawList* dl = ImPlot::GetPlotDrawList();
            const ImVec2 a = ImPlot::PlotToPixels(b.x0, b.y0);
            dl->AddRectFilled(a, io.MousePos, IM_COL32(120, 170, 255, 45));
            dl->AddRect(a, io.MousePos, IM_COL32(150, 190, 255, 200), 0.0f, 0, 1.5f);
            ImPlot::PopPlotClipRect();
        }
    } else {
        b.active = false;
        if (!b.add) data->movie.sel.clear();
        if (b.moved) pick(MIN(b.x0, p.x), MAX(b.x0, p.x), MIN(b.y0, p.y), MAX(b.y0, p.y));
    }
}

// Once the lanes are drawn: a drag that ended puts the keys in order, a click on a picked key that was not a drag picks it alone
static void movie_key_drag_end(ApplicationState* data) {
    auto& m = data->movie;
    auto& kd = m.key_drag;
    if (ImGui::IsMouseDown(ImGuiMouseButton_Left)) return;
    if (kd.collapse && !(kd.active && kd.moved)) m.sel.set(kd.collapse_to.kind, kd.collapse_to.subject, kd.collapse_to.time);
    kd.collapse = false;
    if (kd.active) {
        MovieKeys cur = movie_keys_snapshot(data);
        movie_keys_resolve(&cur, m.sel);
        movie_keys_restore(data, cur);
        key_selection_prune(&m.sel, cur);
        kd.active = false;
        kd.moved = false;
        kd.start = MovieKeys();
    }
}

static void movie_params_apply(ApplicationState* state, double time);

// Edits of the picked keys that are not drags: arrow keys, delete, copy and paste. They are made from the keys as they are and end
// with the keys in order.
static void movie_selection_edit(ApplicationState* data, double dt, double steps, bool remove, bool paste) {
    auto& m = data->movie;
    MovieKeys keys = movie_keys_snapshot(data);
    const double duration = (double)movie_duration(data);
    if (remove) {
        movie_keys_delete(&keys, &m.sel);
    } else if (paste) {
        movie_keys_paste(&keys, &m.sel, m.key_clip, (double)m.playhead, duration);
    } else {
        KeyShift s = m.sel_lane;
        s.dy = s.lane == KeyLane::None ? 0.0 : s.step_ratio ? pow(1.05, steps) : steps * s.step;
        KeySelection moved;
        MovieKeys out = keys;
        movie_keys_shift(&out, &moved, keys, m.sel, dt, s, duration);
        keys = out;
        m.sel = moved;
        movie_keys_resolve(&keys, m.sel);
    }
    if (remove || paste) movie_keys_resolve(&keys, m.sel);
    movie_keys_restore(data, keys);
    key_selection_prune(&m.sel, keys);
    movie_params_apply(data, (double)m.playhead);
    movie_reps_apply(data, (double)m.playhead);
}

static void movie_selection_select_all(ApplicationState* data) {
    auto& m = data->movie;
    const int num_params = (int)(sizeof(movie_param_table) / sizeof(movie_param_table[0]));
    const int num_reps = (int)md_array_size(data->representation.reps);
    const int64_t param = m.timeline_param_lane ? (int64_t)movie_param_table[CLAMP(m.param_selected, 0, num_params - 1)].id : -1;
    const int64_t rep = m.timeline_rep_lane && num_reps > 0 ? rep_key_subject(data->representation.reps[CLAMP(m.rep_selected, 0, num_reps - 1)].id, m.rep_prop_selected) : -1;
    MovieKeys keys = movie_keys_snapshot(data);
    key_selection_all(&m.sel, keys, param, rep);
}

static void movie_selection_copy(ApplicationState* data) {
    auto& m = data->movie;
    m.key_clip = movie_keys_copy(movie_keys_snapshot(data), m.sel);
}

// Arrow keys move the picked keys (a frame, with Shift a second, with Ctrl ten frames; up and down change the value of the lane they
// were last picked in), Delete removes them, Esc puts them down, Ctrl + A picks all, Ctrl + C and Ctrl + V copy and paste at the preview
// time. They work with the mouse over the lanes.
static void movie_selection_shortcuts(ApplicationState* data, bool hovered) {
    auto& m = data->movie;
    const ImGuiIO& io = ImGui::GetIO();
    if (!hovered || io.WantTextInput || m.state == MovieRecordingState::Recording || ImGui::IsPopupOpen("", ImGuiPopupFlags_AnyPopupId | ImGuiPopupFlags_AnyPopupLevel)) return;
    bool took = false;
    const double frame = 1.0 / MAX((double)m.fps, 1.0);
    const double stride = io.KeyShift ? 1.0 : io.KeyCtrl ? 10.0 * frame : frame;
    const double vsteps = io.KeyShift ? 10.0 : 1.0;
    if (!m.sel.empty()) {
        if (ImGui::IsKeyPressed(ImGuiKey_LeftArrow))       { movie_selection_edit(data, -stride, 0.0, false, false); took = true; }
        else if (ImGui::IsKeyPressed(ImGuiKey_RightArrow)) { movie_selection_edit(data, stride, 0.0, false, false); took = true; }
        else if (ImGui::IsKeyPressed(ImGuiKey_UpArrow))    { movie_selection_edit(data, 0.0, vsteps, false, false); took = true; }
        else if (ImGui::IsKeyPressed(ImGuiKey_DownArrow))  { movie_selection_edit(data, 0.0, -vsteps, false, false); took = true; }
        else if (ImGui::IsKeyPressed(ImGuiKey_Delete) || ImGui::IsKeyPressed(ImGuiKey_Backspace)) { movie_selection_edit(data, 0.0, 0.0, true, false); took = true; }
        else if (ImGui::IsKeyPressed(ImGuiKey_Escape, false)) { m.sel.clear(); took = true; }
        else if (ImGui::IsKeyChordPressed(ImGuiMod_Ctrl | ImGuiKey_C)) { movie_selection_copy(data); took = true; }
    }
    if (ImGui::IsKeyChordPressed(ImGuiMod_Ctrl | ImGuiKey_A)) { movie_selection_select_all(data); took = true; }
    if (!m.key_clip.empty() && ImGui::IsKeyChordPressed(ImGuiMod_Ctrl | ImGuiKey_V)) { movie_selection_edit(data, 0.0, 0.0, false, true); took = true; }
    if (took) m.shortcut_frame = ImGui::GetFrameCount();
}

// A line of text in the corner of a lane, so that the lanes can be told apart
static void movie_lane_title(const char* text) {
    ImPlot::PushPlotClipRect();
    const ImVec2 p = ImPlot::GetPlotPos();
    const ImVec2 size = ImGui::CalcTextSize(text);
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    dl->AddRectFilled(ImVec2(p.x + 4.0f, p.y + 3.0f), ImVec2(p.x + size.x + 12.0f, p.y + size.y + 7.0f), IM_COL32(0, 0, 0, 110), 3.0f);
    dl->AddText(ImVec2(p.x + 8.0f, p.y + 5.0f), IM_COL32(230, 230, 230, 210), text);
    ImPlot::PopPlotClipRect();
}

// The time ruler: the seconds of the axis, a tick for every frame of the movie when they are far enough apart (numbered), the notes
// of the movie as triangles (click one to go there), and clicking or dragging anywhere in it moves the preview time.
// Numbers on the lanes' axes and in the mouse readout, short: two decimals at most (one from 100, none from 10000, three significant
// digits below 1), without trailing zeros
static int movie_lane_number_format(double value, char* buf, int size, void*) {
    const double a = fabs(value);
    int len;
    if (a < 1e-9) len = snprintf(buf, (size_t)size, "0");
    else if (a < 1.0) len = snprintf(buf, (size_t)size, "%.3g", value);
    else len = snprintf(buf, (size_t)size, "%.*f", a >= 10000.0 ? 0 : a >= 100.0 ? 1 : 2, value);
    if (len > 0 && len < size && strchr(buf, '.') && !strchr(buf, 'e')) {
        while (len > 0 && buf[len - 1] == '0') buf[--len] = '\0';
        if (len > 0 && buf[len - 1] == '.') buf[--len] = '\0';
    }
    if (len == 2 && buf[0] == '-' && buf[1] == '0') { buf[0] = '0'; buf[1] = '\0'; len = 1; }
    return len;
}

static void movie_lane_axes_format(bool y2 = false) {
    ImPlot::SetupAxisFormat(ImAxis_X1, movie_lane_number_format);
    ImPlot::SetupAxisFormat(ImAxis_Y1, movie_lane_number_format);
    if (y2) ImPlot::SetupAxisFormat(ImAxis_Y2, movie_lane_number_format);
}

static void draw_movie_ruler(ApplicationState* data, float movie_len, bool locked) {
    auto& m = data->movie;
    if (movie_len <= 0.0f) return;
    const ImPlotFlags plot_flags = ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMouseText | ImPlotFlags_NoTitle | ImPlotFlags_NoLegend;
    static bool scrubbing = false;
    if (!ImPlot::BeginPlot("##movie_ruler", ImVec2(-1, -1), plot_flags)) return;
    ImPlot::SetupAxes(nullptr, nullptr, 0, ImPlotAxisFlags_Lock | ImPlotAxisFlags_NoDecorations);
    movie_lane_axes_format();
    ImPlot::SetupAxisLinks(ImAxis_X1, &m.timeline_view_begin, &m.timeline_view_end);
    ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.0, ImPlotCond_Always);
    ImPlot::SetupFinish();

    const ImGuiIO& io = ImGui::GetIO();
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    const ImVec2 p0 = ImPlot::GetPlotPos(), size = ImPlot::GetPlotSize();
    const double fps = MAX((double)m.fps, 1.0);
    const double span = MAX(m.timeline_view_end - m.timeline_view_begin, 1.0e-9);
    const double px_per_frame = (double)size.x / span / fps;
    const int num_frames = movie_num_frames(data);
    const float fs = ImGui::GetFontSize();

    ImPlot::PushPlotClipRect();
    if (px_per_frame >= 4.0) {
        static const int steps[] = {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000};
        int step = steps[(sizeof(steps) / sizeof(steps[0])) - 1];
        for (int s : steps) {
            if (px_per_frame * s >= 44.0) { step = s; break; }
        }
        const int f0 = MAX((int)floor(m.timeline_view_begin * fps), 0), f1 = MIN((int)ceil(m.timeline_view_end * fps), num_frames - 1);
        for (int f = f0; f <= f1; ++f) {
            const float x = ImPlot::PlotToPixels((double)f / fps, 0.0).x;
            const bool major = f % step == 0;
            dl->AddLine(ImVec2(x, p0.y + size.y - (major ? 11.0f : 5.0f)), ImVec2(x, p0.y + size.y), IM_COL32(200, 200, 200, major ? 210 : 110));
            if (major) {
                char buf[16];
                snprintf(buf, sizeof(buf), "%d", f);
                dl->AddText(ImVec2(x + 3.0f, p0.y + size.y - fs - 3.0f), IM_COL32(170, 170, 170, 220), buf);
            }
        }
    }

    // The notes of the movie
    // Their labels go right of the note in the first half of the lane, left of it in the second, and a line lower
    // (above the frame numbers) when they would overlap one written before
    int hovered_marker = -1;
    md_array(int) note_order = 0;
    for (int i = 0; i < (int)m.markers.size(); ++i) md_array_push(note_order, i, frame_alloc);
    std::sort(note_order, note_order + md_array_size(note_order), [&m](int a, int b) { return m.markers[a].time < m.markers[b].time; });
    md_array(MovieLabelSpan) note_labels = 0;
    const float note_lo = ImPlot::GetPlotPos().x, note_hi = note_lo + ImPlot::GetPlotSize().x;
    const int note_rows = MAX(1, (int)((size.y - fs - 5.0f) / fs));
    for (size_t k = 0; k < md_array_size(note_order); ++k) {
        const int i = note_order[k];
        const MovieMarker& mk = m.markers[i];
        const float x = ImPlot::PlotToPixels(mk.time, 0.0).x;
        const float y = p0.y + 4.0f;
        float mc[4];
        movie_marker_color(m.markers.data(), m.markers.size(), (size_t)i, mc);
        const ImU32 mcol = ImGui::ColorConvertFloat4ToU32(ImVec4(mc[0], mc[1], mc[2], 0.95f));
        dl->AddTriangleFilled(ImVec2(x - 5.0f, y), ImVec2(x + 5.0f, y), ImVec2(x, y + 9.0f), mcol);
        if (mk.label[0] != '\0' && x >= note_lo && x <= note_hi) {
            const float tw = ImGui::CalcTextSize(mk.label).x;
            float lx = x;
            const int row = movie_label_place(note_labels, md_array_size(note_labels), x, tw, note_lo, note_hi, 7.0f, fs * 0.4f, note_rows, &lx);
            md_array_push(note_labels, (MovieLabelSpan{ lx, lx + tw, row }), frame_alloc);
            dl->AddText(ImVec2(lx, y - 2.0f + (float)row * fs), mcol, mk.label);
        }
        if (ImPlot::IsPlotHovered() && fabsf(io.MousePos.x - x) < 7.0f && io.MousePos.y > y - 3.0f && io.MousePos.y < y + 13.0f) hovered_marker = i;
    }

    const float ph = ImPlot::PlotToPixels((double)m.playhead, 0.0).x;
    dl->AddLine(ImVec2(ph, p0.y), ImVec2(ph, p0.y + size.y), IM_COL32(255, 255, 0, 255), 1.5f);
    if (locked) {
        const float cur = ImPlot::PlotToPixels(m.cur_time, 0.0).x;
        dl->AddLine(ImVec2(cur, p0.y), ImVec2(cur, p0.y + size.y), IM_COL32(255, 80, 80, 255), 1.5f);
    }
    ImPlot::PopPlotClipRect();

    if (hovered_marker >= 0) {
        const MovieMarker& mk = m.markers[hovered_marker];
        ImGui::SetTooltip("%s\n%.2f s. Click to go there.", mk.label[0] ? mk.label : "Note", mk.time);
    }
    if (!locked) {
        if (ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            if (hovered_marker >= 0) {
                m.playhead = (float)CLAMP(m.markers[hovered_marker].time, 0.0, (double)movie_len);
                movie_apply_time(data, (double)m.playhead, true);
            } else {
                scrubbing = true;
            }
        }
        if (scrubbing) {
            if (io.MouseDown[0]) {
                m.playhead = (float)movie_snap_time(data, ImPlot::GetPlotMousePos().x);
                movie_apply_time(data, (double)m.playhead, true);
            } else {
                scrubbing = false;
            }
        }
    } else {
        scrubbing = false;
    }
    movie_lane_title("Time");
    ImPlot::EndPlot();
}

// Subplots align the time axes and provide draggable row splitters.
static void draw_movie_strip(ApplicationState* data, float movie_len, bool locked, ImVec2 size) {
    auto& m = data->movie;
    if (movie_len <= 0.0f) return;
    const size_t n = md_array_size(m.keyframes);

    // Evaluated on a sorted copy, the keys themselves are only re-sorted when a drag ends
    std::vector<CameraKeyframe> sorted(m.keyframes, m.keyframes + n);
    std::stable_sort(sorted.begin(), sorted.end(), [](const CameraKeyframe& a, const CameraKeyframe& b) { return a.time < b.time; });

    constexpr int N = 200;
    float xs[N], dist[N], fov[N];
    const MovieResolved resolved = n > 0 ? movie_resolve(data, m.keyframes, n) : MovieResolved{};
    char distance_unit[32], distance_axis[64];
    const double distance_scale = display_units::factor_print(distance_unit, sizeof(distance_unit), md_unit_angstrom());
    snprintf(distance_axis, sizeof(distance_axis), "Distance (%s)", distance_unit);
    for (int i = 0; i < N; ++i) {
        xs[i] = movie_len * (float)i / (float)(N - 1);
        if (n > 0) {
            ViewTransform vt;
            float fov_y;
            const double t = xs[i];
            movie_pose(data, resolved, t, &vt, &fov_y);
            dist[i] = (float)(vt.distance * distance_scale);
            fov[i] = fov_y * MOVIE_RAD_TO_DEG;
        } else {
            dist[i] = (float)(data->view.camera.distance * distance_scale);
            fov[i] = data->view.camera.fov_y * MOVIE_RAD_TO_DEG;
        }
    }

    // The trajectory frame over the movie, against the whole trajectory (not scaled to fit, so that a point
    // can be dragged to a frame)
    const double last_frame = (double)(run_num_frames(data) > 0 ? run_num_frames(data) - 1 : 0);
    float frm[N];
    for (int i = 0; i < N; ++i) {
        const double t = (double)movie_len * (double)i / (double)(N - 1);
        frm[i] = (float)movie_trajectory_frame(data, t);
    }
    const ImPlotDragToolFlags drag_flags = ImPlotDragToolFlags_NoFit | (locked ? ImPlotDragToolFlags_NoInputs : 0);
    static bool resort_pending = false;

    const char* titles[3] = {"Trajectory", "Distance", "Field of view"};
    const char* axes[3] = {"Trajectory frame", distance_axis, "Field of view (deg)"};
    const ImVec4 colors[3] = {ImVec4(0.4f, 0.9f, 0.4f, 1), ImVec4(0.35f, 0.8f, 1, 1), ImVec4(1, 0.8f, 0.25f, 1)};
    int rows = 0;
    float ratios[9];
    int ratio_slot[9];  // Which of timeline_row_ratios each row is
    int overview_row = -1, overlay_row = -1, camera_row = -1;
    auto add_row = [&](int slot) { ratio_slot[rows] = slot; ratios[rows] = m.timeline_row_ratios[slot]; return rows++; };
    if (m.timeline_ruler) add_row(8);
    if (m.timeline_camera_lane) camera_row = add_row(7);
    if (m.timeline_tracks[0]) add_row(0);
    const bool lens_lane = m.timeline_tracks[1] || m.timeline_tracks[2];
    if (lens_lane) add_row(1);
    if (m.timeline_rep_overview) overview_row = add_row(6);
    if (m.timeline_overlay_lane) overlay_row = add_row(5);
    if (m.timeline_rep_lane) add_row(4);
    if (m.timeline_param_lane) add_row(3);
    const ImPlotFlags plot_flags = ImPlotFlags_NoBoxSelect | ImPlotFlags_NoLegend;
    if (rows == 0) {
        ImGui::TextDisabled("No lane is ticked.");
        return;
    }
    // With a fixed lane height, every lane is sized in pixels (the lane height times its ratio, or what its rows need
    // when that is more) and the window scrolls. Otherwise the lanes share the given height by their ratios.
    const bool fixed_height = size.y < 0.0f;
    float given[9];
    if (fixed_height) {
        const float fs = ImGui::GetFontSize();
        const int list_rows[3] = {overview_row, overlay_row, camera_row};
        const float list_needed[3] = {movie_rep_overview_height(data, movie_len), (float)m.overlays.size() * fs * 1.5f + fs * 4.0f, fs * 1.5f * 4.0f + fs * 4.0f};
        float total = 0.0f;
        for (int i = 0; i < rows; ++i) {
            ratios[i] = MAX(ratios[i], 0.05f) * m.timeline_lane_height;
            for (int k = 0; k < 3; ++k) {
                if (list_rows[k] == i) ratios[i] = MAX(ratios[i], list_needed[k]);
            }
            total += ratios[i];
        }
        size.y = total;
    }
    for (int i = 0; i < rows; ++i) given[i] = ratios[i];
    MovieCameraLaneEdit camera_edit;
    if (!m.sel.empty()) {
        MovieKeys cur;
        cur.camera.assign(m.keyframes, m.keyframes + n);
        cur.params = m.param_keys;
        cur.reps = m.rep_keys;
        cur.overlays = m.overlays;
        cur.duration = m.duration;
        key_selection_prune(&m.sel, cur);
    }
    // Dragging the background picks keys with a box, so the middle button pans (Shift + wheel and sideways scrolling pan too)
    const ImPlotInputMap old_input_map = ImPlot::GetInputMap();
    ImPlot::GetInputMap().Pan = ImGuiMouseButton_Middle;
    if (ImPlot::BeginSubplots("##movie_tracks", rows, 1, size, ImPlotSubplotFlags_NoTitle, ratios)) {
      if (m.timeline_ruler) draw_movie_ruler(data, movie_len, locked);
      if (m.timeline_camera_lane) draw_movie_camera_lane(data, movie_len, locked, sorted, &camera_edit);
      // The trajectory lane, then the camera lens lane: field of view on the left axis, distance on the right
      for (int lane_i = 0; lane_i < 2; ++lane_i) {
        if (lane_i == 0 ? !m.timeline_tracks[0] : !lens_lane) continue;
        int curves[2];
        int num_curves = 0;
        if (lane_i == 0) curves[num_curves++] = 0;
        else {
            curves[num_curves++] = 2;
            curves[num_curves++] = 1;
        }
        const bool lens = lane_i == 1;
        if (!ImPlot::BeginPlot(lens ? "Camera lens" : titles[0], ImVec2(-1, -1), lens && num_curves > 1 ? ImPlotFlags_NoBoxSelect : plot_flags)) continue;
        ImPlot::SetupAxes(nullptr, axes[curves[0]]);
        if (num_curves > 1) ImPlot::SetupAxis(ImAxis_Y2, axes[curves[1]], ImPlotAxisFlags_AuxDefault);
        movie_lane_axes_format(num_curves > 1);
        if (lens && num_curves > 1) ImPlot::SetupLegend(ImPlotLocation_NorthWest, ImPlotLegendFlags_Horizontal | ImPlotLegendFlags_NoButtons);
        ImPlot::SetupAxisLinks(ImAxis_X1, &m.timeline_view_begin, &m.timeline_view_end);
        for (int ci = 0; ci < num_curves; ++ci) {
            const int track = curves[ci];
            const float* values = track == 0 ? frm : track == 1 ? dist : fov;
            double lo = values[0], hi = values[0];
            for (int i = 1; i < N; ++i) {
                lo = MIN(lo, (double)values[i]);
                hi = MAX(hi, (double)values[i]);
            }
            const double padding = MAX((hi - lo) * 0.08, MAX(fabs(hi) * 0.05, 1.0e-3));
            // Refit when the curve's range changes (keys added, a workspace loaded), but not while something is dragged
            static double fitted[3][2] = {{1, 0}, {1, 0}, {1, 0}};
            const bool range_changed = fitted[track][0] != lo || fitted[track][1] != hi;
            const bool refit = range_changed && !ImGui::IsMouseDown(ImGuiMouseButton_Left);
            if (refit) {
                fitted[track][0] = lo;
                fitted[track][1] = hi;
            }
            ImPlot::SetupAxisLimits(ci == 0 ? ImAxis_Y1 : ImAxis_Y2, lo - padding, hi + padding, refit ? ImPlotCond_Always : ImPlotCond_Once);
        }

        bool anchors_held = false;
        if (lane_i == 0) {
            const double tb = (double)m.traj_begin;
            const double te = (double)m.traj_end;

            const ImVec4 anchor_col(0.3f, 0.6f, 1.0f, 1.0f);
            double ab = tb, ae = te;
            bool clicked = false, hov = false, hld = false;
            if (ImPlot::DragLineX(1100, &ab, anchor_col, 2.0f, drag_flags, &clicked, &hov, &hld)) {
                m.traj_begin = (float)CLAMP(movie_snap_time(data, ab), 0.0, (double)m.traj_end);
            }
            if (hov && !hld) ImGui::SetTooltip("Trajectory starts at frame %.0f, at %.2f s\nDrag to change when. Before it, the trajectory is held.", m.start_frame, m.traj_begin);
            anchors_held |= hld;
            if (ImPlot::DragLineX(1101, &ae, anchor_col, 2.0f, drag_flags, &clicked, &hov, &hld)) {
                m.traj_end = (float)CLAMP(movie_snap_time(data, ae), (double)m.traj_begin, (double)movie_len);
            }
            if (hov && !hld) ImGui::SetTooltip("Trajectory ends at frame %.0f, at %.2f s\nDrag to change when. After it, the trajectory is held.", m.end_frame, m.traj_end);
            anchors_held |= hld;
            ImPlot::PlotText("Start", m.traj_begin, m.start_frame, ImVec2(16, -12));
            ImPlot::PlotText("End", m.traj_end, m.end_frame, ImVec2(-16, -12));
        }

        bool any_held = false, over_key = false;
        for (int ci = 0; ci < num_curves; ++ci) {
            const int track = curves[ci];
            ImPlot::SetAxes(ImAxis_X1, ci == 0 ? ImAxis_Y1 : ImAxis_Y2);
            ImPlot::PushStyleColor(ImPlotCol_Line, colors[track]);
            ImPlot::PlotLine(titles[track], xs, track == 0 ? frm : track == 1 ? dist : fov, N);
            ImPlot::PopStyleColor();

            KeyShift lane;
            lane.lane = track == 0 ? KeyLane::Frame : track == 1 ? KeyLane::Distance : KeyLane::Fov;
            lane.lo = 0.0;
            lane.hi = last_frame;
            lane.unit = distance_scale;
            lane.step_ratio = track == 1;
            // On the lens lane each curve's keys have its colour, so it is clear which axis a key is read on
            const ImVec4 key_col = lens ? colors[track] : ImVec4(1.0f, 0.45f, 0.15f, 1.0f);
            for (size_t i = 0; i < n; ++i) {
                const CameraKeyframe key = m.keyframes[i];
                if (track == 0 && !key.use_frame) continue;
                ViewTransform live = key.transform;
                if (track == 1) {
                    float unused_fov;
                    movie_camera_evaluate(data, key.time, m.keyframes, n, &live, &unused_fov);
                }
                const double key_y = track == 0 ? key.frame : track == 1 ? live.distance * distance_scale : key.fov_y * MOVIE_RAD_TO_DEG;
                const MovieKeyPointResult r = movie_key_point(data, locked, KeyKind::Camera, 0, key.time, key_y, 2000 + 2 * (int)i + ci, key_col, 7.0f, lane, &any_held);
                over_key |= r.hovered;
                char label[16];
                snprintf(label, sizeof(label), "%d", (int)i + 1);
                ImPlot::PlotText(label, key.time, key_y, ImVec2(0, -14));
                if (r.hovered && !r.held) {
                    ImGui::SetTooltip("Keyframe %d: %.2f s\n%s: %.2f\nDrag to change the time and the value. Click picks it, Ctrl + click adds it to the picked keys:\nthey move together. Drag the background to pick with a box.", (int)i + 1, key.time, axes[track], key_y);
                }
            }
        }
        ImPlot::SetAxes(ImAxis_X1, ImAxis_Y1);
        {
            double vlines[3] = {(double)m.playhead, (double)m.traj_begin, (double)m.traj_end};
            movie_lane_box(data, locked, 10 + lane_i, over_key, vlines, lane_i == 0 ? 3 : 1, [&](double x0, double x1, double y0, double y1) {
                // The box is in the units of the left axis; a key of the right axis is in it when its pixel is
                const float py0 = ImPlot::PlotToPixels(x0, y0, ImAxis_X1, ImAxis_Y1).y, py1 = ImPlot::PlotToPixels(x0, y1, ImAxis_X1, ImAxis_Y1).y;
                for (size_t i = 0; i < md_array_size(m.keyframes); ++i) {
                    const CameraKeyframe& k = m.keyframes[i];
                    if (k.time < x0 || k.time > x1) continue;
                    for (int ci = 0; ci < num_curves; ++ci) {
                        const int track = curves[ci];
                        if (track == 0 && !k.use_frame) continue;
                        ViewTransform live = k.transform;
                        if (track == 1) {
                            float unused_fov;
                            movie_camera_evaluate(data, k.time, m.keyframes, md_array_size(m.keyframes), &live, &unused_fov);
                        }
                        const double ky = track == 0 ? k.frame : track == 1 ? live.distance * distance_scale : k.fov_y * MOVIE_RAD_TO_DEG;
                        const float py = ImPlot::PlotToPixels(k.time, ky, ImAxis_X1, ci == 0 ? ImAxis_Y1 : ImAxis_Y2).y;
                        if (py >= MIN(py0, py1) && py <= MAX(py0, py1)) m.sel.add(KeyKind::Camera, 0, k.time);
                    }
                }
            });
        }

        double playhead = (double)m.playhead;
        if (ImPlot::DragLineX(1000, &playhead, ImVec4(1, 1, 0, 1), 1.5f, drag_flags)) {
            m.playhead = (float)movie_snap_time(data, playhead);
            if (!locked) movie_apply_time_with_keys(data, (double)m.playhead, true, sorted.data(), n);
        }

        if (locked) {
            double cur = m.cur_time;
            ImPlot::DragLineX(1001, &cur, ImVec4(1.0f, 0.3f, 0.3f, 1), 1.5f, ImPlotDragToolFlags_NoInputs | ImPlotDragToolFlags_NoFit);
        }

        ImPlot::EndPlot();

        if (any_held || anchors_held) {
            sorted.assign(m.keyframes, m.keyframes + md_array_size(m.keyframes));
            std::stable_sort(sorted.begin(), sorted.end(), [](const CameraKeyframe& a, const CameraKeyframe& b) { return a.time < b.time; });
            movie_apply_time_with_keys(data, (double)m.playhead, true, sorted.data(), sorted.size());
        }
      }
      if (m.timeline_rep_overview) draw_movie_rep_overview_lane(data, movie_len, locked);
      if (m.timeline_overlay_lane) draw_movie_overlay_lane(data, movie_len, locked);
      if (m.timeline_rep_lane) draw_movie_rep_lane(data, movie_len, locked);
      if (m.timeline_param_lane) draw_movie_param_lane(data, movie_len, locked);
      if (ImPlot::IsSubplotsHovered() && !ImGui::GetIO().KeyCtrl) {
          const ImGuiIO& io = ImGui::GetIO();
          const double wheel = io.MouseWheelH != 0.0f ? (double)io.MouseWheelH : io.KeyShift ? -(double)io.MouseWheel : 0.0;
          if (wheel != 0.0) {
              const double shift = wheel * 0.1 * (m.timeline_view_end - m.timeline_view_begin);
              m.timeline_view_begin += shift;
              m.timeline_view_end += shift;
          }
      }
      ImPlot::EndSubplots();
      for (int i = 0; i < rows; ++i) {
          if (!fixed_height) {
              m.timeline_row_ratios[ratio_slot[i]] = ratios[i];
          } else if (fabsf(ratios[i] - given[i]) > 0.5f) {
              // A divider was dragged: the new height in pixels, as a multiple of the lane height
              m.timeline_row_ratios[ratio_slot[i]] = MAX(ratios[i] / m.timeline_lane_height, 0.05f);
          }
      }

      if (!locked) {
        if (camera_edit.remove >= 0 && camera_edit.remove < (int)md_array_size(m.keyframes)) {
            CameraKeyframe* keys = m.keyframes;
            const size_t count = md_array_size(keys);
            memmove(keys + camera_edit.remove, keys + camera_edit.remove + 1, (count - (size_t)camera_edit.remove - 1) * sizeof(CameraKeyframe));
            md_array_pop(keys);
        }
        if (camera_edit.insert >= 0.0) {
            bool taken = false;
            for (size_t i = 0; i < md_array_size(m.keyframes); ++i) taken |= fabs(m.keyframes[i].time - camera_edit.insert) < 1.0e-3;
            if (!taken) {
                CameraKeyframe key;
                if (md_array_size(m.keyframes) > 0) {
                    const MovieResolved r = movie_resolve(data, m.keyframes, md_array_size(m.keyframes));
                    ViewTransform pose;
                    float fov_y;
                    movie_pose(data, r, camera_edit.insert, &pose, &fov_y);
                    key = camera_key_on_path(pose, fov_y, r.keys, camera_edit.insert, movie_upright(data));
                } else {
                    key = movie_current_key(data);
                    key.time = camera_edit.insert;
                }
                movie_insert_key(data, key, false);
                m.playhead = (float)camera_edit.insert;
            }
        }
        if (camera_edit.held) movie_apply_time_with_keys(data, (double)m.playhead, true, sorted.data(), n);
      }
    }
    ImPlot::GetInputMap() = old_input_map;
    movie_key_drag_end(data);
    if (resort_pending && !ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
        movie_sort_keyframes(data);
        resort_pending = false;
    }
}

// The camera on one lane: its keys (numbered, with their names), where the look-at point follows the target or an atom, where
// the camera spins, and the keys that pin a trajectory frame. Drag a key sideways to change when, click it to go there, right
// click for its name and removal. A double click on an empty place adds a key on the path, without moving the camera.
// What a key's look-at is ("3 atoms", "resname(\"AIN\")"), and what its focus is, for the lane and the key table
static std::string movie_set_label(const ApplicationState* state, uint32_t id) {
    const MovieAtomSet* set = movie_atoms_find(state->movie.atom_sets, id);
    char buf[48];
    if (set && !set->expr.empty()) return set->expr.size() > 24 ? set->expr.substr(0, 23) + "..." : set->expr;
    if (!set || set->atoms.empty()) return "atoms";
    if (set->atoms.size() == 1) snprintf(buf, sizeof(buf), "atom %u", set->atoms[0] + 1);
    else snprintf(buf, sizeof(buf), "%zu atoms", set->atoms.size());
    return buf;
}

static std::string movie_focus_label(const ApplicationState* state, const CameraKeyframe& key) {
    char buf[96];
    switch (key.focus_target) {
    case FocusTarget::Selection: snprintf(buf, sizeof(buf), "%s", movie_set_label(state, key.focus_set).c_str()); break;
    case FocusTarget::LookAt:    snprintf(buf, sizeof(buf), "look-at"); break;
    case FocusTarget::Distance:  snprintf(buf, sizeof(buf), "%.1f A", key.focus_distance); break;
    default:                     snprintf(buf, sizeof(buf), "point"); break;
    }
    return buf;
}

static void draw_movie_camera_lane(ApplicationState* data, float movie_len, bool locked, const std::vector<CameraKeyframe>& sorted, MovieCameraLaneEdit* edit) {
    auto& m = data->movie;
    if (movie_len <= 0.0f) return;
    const size_t n = md_array_size(m.keyframes);
    const ImPlotDragToolFlags drag_flags = ImPlotDragToolFlags_NoFit | (locked ? ImPlotDragToolFlags_NoInputs : 0);
    const ImPlotFlags plot_flags = ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMouseText | ImPlotFlags_NoTitle | ImPlotFlags_NoLegend;
    static const char* row_labels[5] = {"Camera keys", "Look at", "Focus", "Spin", "Frame"};
    static const double row_pos[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    static int menu_key = -1;
    int open_menu = -1;
    bool any_hovered = false;

    if (ImPlot::BeginPlot("##camera_lane", ImVec2(-1, -1), plot_flags)) {
        ImPlot::SetupAxes("Movie time (s)", nullptr, 0, ImPlotAxisFlags_Lock | ImPlotAxisFlags_Invert);
        movie_lane_axes_format();
        ImPlot::SetupAxisLinks(ImAxis_X1, &m.timeline_view_begin, &m.timeline_view_end);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -0.6, 4.6, ImPlotCond_Always);
        ImPlot::SetupAxisTicks(ImAxis_Y1, row_pos, 5, row_labels);
        ImPlot::SetupFinish();
        movie_lane_title("Camera");

        ImDrawList* dl = ImPlot::GetPlotDrawList();
        const ImVec2 mouse = ImGui::GetMousePos();
        const bool plot_hovered = ImPlot::IsPlotHovered();
        auto px = [&](double t, double y) { return ImPlot::PlotToPixels(t, y); };

        // Spans between keys
        ImPlot::PushPlotClipRect();
        std::string band_tip;
        // After the last key (and before the first) the camera holds still: shaded, so that a Preview there is not
        // mistaken for a camera that does not move
        if (!sorted.empty()) {
            const ImU32 hold_fill = IM_COL32(128, 128, 128, 40), hold_line = IM_COL32(160, 160, 160, 70);
            auto hold_span = [&](double t0, double t1, const char* tip) {
                if (t1 - t0 <= 1.0e-6) return;
                const ImVec2 a = px(t0, -0.5), c = px(t1, 4.5);
                dl->AddRectFilled(a, c, hold_fill);
                const float step = ImGui::GetFontSize() * 0.8f;
                dl->PushClipRect(a, c, true);
                for (float x = a.x - (c.y - a.y); x < c.x; x += step) dl->AddLine(ImVec2(x, c.y), ImVec2(x + (c.y - a.y), a.y), hold_line);
                dl->PopClipRect();
                if (plot_hovered && mouse.x >= a.x && mouse.x <= c.x && mouse.y >= a.y && mouse.y <= c.y) band_tip = tip;
            };
            hold_span(0.0, sorted.front().time, "Before the first key the camera holds still on it");
            hold_span(sorted.back().time, (double)movie_len, "After the last key the camera holds still: add keys here for it to move,\nor shorten the movie (Timing > Movie length)");
        }
        for (const CameraBand& b : camera_bands(sorted)) {
            const int row = b.kind == CameraBandKind::Spin ? 3 : 1;
            ImVec2 a = px(b.begin, (double)row - 0.3), c = px(b.end, (double)row + 0.3);
            if (c.x - a.x < 8.0f) { a.x -= 4.0f; c.x += 4.0f; }
            ImVec4 col(0.2f, 0.75f, 0.7f, 1.0f);
            std::string text, tip;
            char buf[160];
            if (b.kind == CameraBandKind::LookAtSet) {
                col = ImVec4(0.45f, 0.6f, 1.0f, 1.0f);
                text = movie_set_label(data, b.set);
                snprintf(buf, sizeof(buf), "Looks at %s, tracked through the trajectory\nKeys %d to %d, %.2f s to %.2f s", text.c_str(), b.first + 1, b.last + 1, b.begin, b.end);
            } else {
                col = ImVec4(0.75f, 0.45f, 0.95f, 1.0f);
                snprintf(buf, sizeof(buf), "%+d x", b.turns);
                text = buf;
                snprintf(buf, sizeof(buf), "%+d turns around %s, from key %d to key %d\n%.2f s to %.2f s", b.turns, spin_axis_str[(int)sorted[(size_t)b.last].spin_axis], b.first + 1, b.last + 1, b.begin, b.end);
            }
            tip = buf;
            dl->AddRectFilled(a, c, ImGui::ColorConvertFloat4ToU32(ImVec4(col.x, col.y, col.z, 0.45f)), 3.0f);
            dl->AddRect(a, c, ImGui::ColorConvertFloat4ToU32(ImVec4(col.x, col.y, col.z, 0.9f)), 3.0f);
            dl->PushClipRect(a, c, true);
            dl->AddText(ImVec2(a.x + 4.0f, 0.5f * (a.y + c.y) - 0.5f * ImGui::GetFontSize()), IM_COL32(255, 255, 255, 235), text.c_str());
            dl->PopClipRect();
            if (plot_hovered && mouse.x >= a.x && mouse.x <= c.x && mouse.y >= a.y && mouse.y <= c.y) band_tip = tip;
        }

        // What is sharp: from each key that sets it to the next one that does, a band (its transition darker)
        {
            const ImVec4 col(1.0f, 0.4f, 0.85f, 1.0f);
            for (size_t i = 0; i < sorted.size(); ++i) {
                if (!sorted[i].focus_on) continue;
                double end = (double)movie_len;
                for (size_t j = i + 1; j < sorted.size(); ++j) {
                    if (sorted[j].focus_on) { end = sorted[j].time; break; }
                }
                ImVec2 a = px(sorted[i].time, 2.0 - 0.3), c = px(MAX(end, sorted[i].time), 2.0 + 0.3);
                if (c.x - a.x < 8.0f) c.x = a.x + 8.0f;
                const ImVec2 t = px(MIN(sorted[i].time + (double)sorted[i].focus_transition, end), 2.3);
                dl->AddRectFilled(a, c, ImGui::ColorConvertFloat4ToU32(ImVec4(col.x, col.y, col.z, 0.25f)), 3.0f);
                dl->AddRectFilled(a, ImVec2(MAX(t.x, a.x), c.y), ImGui::ColorConvertFloat4ToU32(ImVec4(col.x, col.y, col.z, 0.5f)), 3.0f);
                dl->AddRect(a, c, ImGui::ColorConvertFloat4ToU32(ImVec4(col.x, col.y, col.z, 0.9f)), 3.0f);
                const std::string text = movie_focus_label(data, sorted[i]);
                dl->PushClipRect(a, c, true);
                dl->AddText(ImVec2(a.x + 4.0f, 0.5f * (a.y + c.y) - 0.5f * ImGui::GetFontSize()), IM_COL32(255, 255, 255, 235), text.c_str());
                dl->PopClipRect();
                if (plot_hovered && mouse.x >= a.x && mouse.x <= c.x && mouse.y >= a.y && mouse.y <= c.y) {
                    char tip[200];
                    snprintf(tip, sizeof(tip), "Key %d sets what is sharp: %s\nChange over %.1f s, then it holds until the next key that sets it.\nHow blurred the rest is: Looks > Depth of field blur", (int)i + 1, text.c_str(), sorted[i].focus_transition);
                    band_tip = tip;
                }
            }
        }

        // The keys that pin a trajectory frame
        for (size_t i = 0; i < sorted.size(); ++i) {
            if (!sorted[i].use_frame) continue;
            const ImVec2 p = px(sorted[i].time, 4.0);
            const float r = ImGui::GetFontSize() * 0.4f;
            dl->AddQuadFilled(ImVec2(p.x, p.y - r), ImVec2(p.x + r, p.y), ImVec2(p.x, p.y + r), ImVec2(p.x - r, p.y), IM_COL32(110, 230, 110, 255));
            char buf[32];
            snprintf(buf, sizeof(buf), "%.0f", sorted[i].frame);
            dl->AddText(ImVec2(p.x + r + 3.0f, p.y - 0.5f * ImGui::GetFontSize()), IM_COL32(200, 255, 200, 255), buf);
            if (plot_hovered && fabsf(mouse.x - p.x) < r * 1.5f && fabsf(mouse.y - p.y) < r * 1.5f) {
                char tip[96];
                snprintf(tip, sizeof(tip), "Key %d shows trajectory frame %.0f at %.2f s", (int)i + 1, sorted[i].frame, sorted[i].time);
                band_tip = tip;
            }
        }
        ImPlot::PopPlotClipRect();

        // The keys, which are dragged sideways only
        for (size_t i = 0; i < n; ++i) {
            const CameraKeyframe key = m.keyframes[i];
            const bool at_playhead = fabs(key.time - (double)m.playhead) < 1.0e-3;
            const ImVec4 col = at_playhead ? ImVec4(1.0f, 0.85f, 0.3f, 1.0f) : ImVec4(1.0f, 0.45f, 0.15f, 1.0f);
            const MovieKeyPointResult r = movie_key_point(data, locked, KeyKind::Camera, 0, key.time, 0.0, 9500 + (int)i, col, m.path_hover_key == (int)i ? 11.0f : 7.0f, KeyShift{}, &edit->held);
            const bool hovered = r.hovered, held = r.held, clicked = r.clicked;
            any_hovered |= hovered;
            if (hovered) m.lane_hover_key = (int)i;
            if (clicked && !locked && !ImGui::GetIO().KeyCtrl && !(m.key_drag.active && m.key_drag.moved)) movie_goto_keyframe(data, i);
            if (hovered && !locked && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) open_menu = (int)i;
            const std::string label = camera_key_label(key, (int)i);
            ImPlot::PlotText(label.c_str(), key.time, 0.0, ImVec2(0.0f, -16.0f));
            if (hovered && !held) {
                ImGui::SetTooltip("Keyframe %s\n%.2f s, field of view %.1f deg%s\nDrag sideways to change when, click to go to it, right click for its name and removal.\nCtrl + click adds it to the picked keys, which move together. Drag the background to pick with a box.",
                    label.c_str(), key.time, key.fov_y * MOVIE_RAD_TO_DEG, key.use_frame ? ", pins a trajectory frame" : "");
            }
        }
        if (!any_hovered && !band_tip.empty()) ImGui::SetTooltip("%s", band_tip.c_str());
        {
            const double vline = (double)m.playhead;
            movie_lane_box(data, locked, 20, any_hovered, &vline, 1, [&](double x0, double x1, double y0, double y1) {
                for (size_t i = 0; i < md_array_size(m.keyframes); ++i) {
                    if (y0 <= 0 && y1 >= 0 && m.keyframes[i].time >= x0 && m.keyframes[i].time <= x1)
                        m.sel.add(KeyKind::Camera, 0, m.keyframes[i].time);
                }
            });
        }

        // A key is added on the path where it is double clicked
        if (!locked && !any_hovered && plot_hovered && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
            const auto mouse_point = ImPlot::GetPlotMousePos();
            const double time = movie_snap_time(data, mouse_point.x);
            edit->insert = time;
        }

        double playhead = (double)m.playhead;
        if (ImPlot::DragLineX(1000, &playhead, ImVec4(1, 1, 0, 1), 1.5f, drag_flags)) {
            m.playhead = (float)movie_snap_time(data, playhead);
            if (!locked) movie_apply_time_with_keys(data, (double)m.playhead, true, sorted.data(), sorted.size());
        }
        if (locked) {
            double cur = m.cur_time;
            ImPlot::DragLineX(1001, &cur, ImVec4(1.0f, 0.3f, 0.3f, 1), 1.5f, ImPlotDragToolFlags_NoInputs | ImPlotDragToolFlags_NoFit);
        }
        if (n == 0) ImPlot::PlotText("Double click to add a key", 0.5 * (double)movie_len, 0.0, ImVec2(0.0f, 0.0f));
        ImPlot::EndPlot();
    }

    if (open_menu >= 0) {
        menu_key = open_menu;
        ImGui::OpenPopup("##camera_key_menu");
    }
    if (ImGui::BeginPopup("##camera_key_menu")) {
        if (menu_key >= 0 && menu_key < (int)n) {
            CameraKeyframe& key = m.keyframes[menu_key];
            ImGui::TextDisabled("Keyframe %d", menu_key + 1);
            if (ImGui::IsWindowAppearing()) ImGui::SetKeyboardFocusHere();
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 12.0f);
            if (ImGui::InputTextWithHint("##camera_key_name", "Name", key.name, sizeof(key.name), ImGuiInputTextFlags_EnterReturnsTrue)) ImGui::CloseCurrentPopup();
            if (ImGui::Selectable("Go to")) movie_goto_keyframe(data, (size_t)menu_key);
            if (ImGui::Selectable("Remove")) {
                edit->remove = menu_key;
                ImGui::CloseCurrentPopup();
            }
        } else {
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
}

// A look parameter over the movie, as a curve whose keys can be dragged: sideways for their time, up and down
// for their value. A double click adds a key, a right click on one removes it. Fills the space that is left.
static void draw_movie_param_lane(ApplicationState* data, float movie_len, bool locked) {
    auto& m = data->movie;
    if (movie_len <= 0.0f) return;
    const int num_params = (int)(sizeof(movie_param_table) / sizeof(movie_param_table[0]));
    m.param_selected = CLAMP(m.param_selected, 0, num_params - 1);
    const MovieParamDesc& d = movie_param_table[m.param_selected];

    std::vector<int> mine;
    for (int i = 0; i < (int)m.param_keys.size(); ++i) {
        if (m.param_keys[i].param == d.id) mine.push_back(i);
    }

    constexpr int N = 200;
    float xs[N], ys[N];
    if (!d.color && !mine.empty()) {
        for (int i = 0; i < N; ++i) {
            const double t = (double)movie_len * (double)i / (double)(N - 1);
            float v = 0.0f;
            param_keys_evaluate(&v, 1, m.param_keys.data(), m.param_keys.size(), d.id, t);
            xs[i] = (float)t;
            ys[i] = v;
        }
    }

    const ImPlotDragToolFlags drag_flags = ImPlotDragToolFlags_NoFit | (locked ? ImPlotDragToolFlags_NoInputs : 0);
    bool any_held = false, any_hovered = false;
    int  remove_idx = -1;
    bool add_key = false;
    double add_time = 0.0, add_value = 0.0;

    const ImPlotFlags plot_flags = ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMouseText | ImPlotFlags_NoTitle | ImPlotFlags_NoLegend;
    if (ImPlot::BeginPlot("##param_lane", ImVec2(-1, -1), plot_flags)) {
        ImPlot::SetupAxes("Movie time (s)", d.color ? nullptr : d.label, 0, ImPlotAxisFlags_Lock | (d.color ? ImPlotAxisFlags_NoDecorations : 0));
        movie_lane_axes_format();
        ImPlot::SetupAxisLinks(ImAxis_X1, &m.timeline_view_begin, &m.timeline_view_end);
        if (d.color) {
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.0, ImPlotCond_Always);
        } else if (d.log) {
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            ImPlot::SetupAxisLimits(ImAxis_Y1, d.lo * 0.8, d.hi * 1.25, ImPlotCond_Always);
        } else {
            const double pad = 0.05 * (d.hi - d.lo);
            ImPlot::SetupAxisLimits(ImAxis_Y1, d.lo - pad, d.hi + pad, ImPlotCond_Always);
        }

        movie_lane_title((std::string("Look parameter: ") + d.label).c_str());
        if (!d.color && !mine.empty()) {
            ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0.75f, 0.5f, 1.0f, 1.0f));
            ImPlot::PlotLine(d.label, xs, ys, N);
            ImPlot::PopStyleColor();
        }

        KeyShift lane;
        if (!d.color) {
            lane.lane = KeyLane::Param;
            lane.subject = d.id;
            lane.lo = d.lo;
            lane.hi = d.hi;
            lane.ratio = d.log;
            lane.step = (d.hi - d.lo) / 100.0;
            lane.step_ratio = d.log;
        }
        for (int ki : mine) {
            const ParamKey key = m.param_keys[ki];
            const ImVec4 col = d.color ? ImVec4(key.value[0], key.value[1], key.value[2], 1.0f) : ImVec4(0.75f, 0.5f, 1.0f, 1.0f);
            const MovieKeyPointResult r = movie_key_point(data, locked, KeyKind::Param, d.id, key.time, d.color ? 0.5 : (double)key.value[0], 4000 + ki, col, 7.0f, lane, &any_held);
            if (r.hovered) {
                any_hovered = true;
                if (!r.held) {
                    if (d.color) ImGui::SetTooltip("%.2f s\nDrag to change its time. Its color is edited in the table of the Movie window.", key.time);
                    else         ImGui::SetTooltip("%.2f s, %.3g\nDrag to change it, right click to remove it.\nCtrl + click adds it to the picked keys, which move together.", key.time, key.value[0]);
                }
                if (!locked && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) remove_idx = ki;
            }
        }
        {
            double vline = (double)m.playhead;
            movie_lane_box(data, locked, 30, any_hovered, &vline, 1, [&](double x0, double x1, double y0, double y1) {
                for (int ki : mine) {
                    const ParamKey& k = m.param_keys[ki];
                    const double ky = d.color ? 0.5 : (double)k.value[0];
                    if (k.time >= x0 && k.time <= x1 && ky >= y0 && ky <= y1) m.sel.add(KeyKind::Param, d.id, k.time);
                }
            });
        }

        if (!locked && !any_hovered && ImPlot::IsPlotHovered() && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
            const ImPlotPoint p = ImPlot::GetPlotMousePos();
            add_key = true;
            add_time = CLAMP(p.x, 0.0, (double)movie_len);
            add_value = CLAMP(p.y, (double)d.lo, (double)d.hi);
        }

        double playhead = (double)m.playhead;
        if (ImPlot::DragLineX(1000, &playhead, ImVec4(1, 1, 0, 1), 1.5f, drag_flags)) {
            m.playhead = (float)movie_snap_time(data, playhead);
            if (!locked) movie_apply_time(data, (double)m.playhead, true);
        }

        if (locked) {
            double cur = m.cur_time;
            ImPlot::DragLineX(1001, &cur, ImVec4(1.0f, 0.3f, 0.3f, 1), 1.5f, ImPlotDragToolFlags_NoInputs | ImPlotDragToolFlags_NoFit);
        }

        ImPlot::EndPlot();
    }

    if (any_held) {
        movie_apply_time(data, (double)m.playhead, true);
    }

    if (remove_idx >= 0) {
        m.param_keys.erase(m.param_keys.begin() + remove_idx);
        // The parameter lets go of its values if this was its last key
        movie_params_apply(data, (double)m.playhead);
    }
    if (add_key) {
        ParamKey key;
        key.param = d.id;
        key.time = add_time;
        if (d.color) {
            // The colour it has at that time, or as it is now if there are no keys yet
            const float* now = d.ptr(data);
            float v[3] = {now[0], now[1], now[2]};
            param_keys_evaluate(v, 3, m.param_keys.data(), m.param_keys.size(), d.id, add_time);
            for (int c = 0; c < 3; ++c) key.value[c] = v[c];
        } else {
            key.value[0] = (float)add_value;
        }
        m.param_keys.push_back(key);
        movie_param_sort(data);
        movie_apply_time(data, (double)m.playhead, true);
    }
}

static void movie_rep_sort(ApplicationState* state);

// A property of a representation over the movie, with its keys to drag: sideways for their time, up and down for
// their value (Visible is shown or hidden). A double click adds a key, a right click on one removes it.
static void draw_movie_rep_lane(ApplicationState* data, float movie_len, bool locked) {
    auto& m = data->movie;
    if (movie_len <= 0.0f) return;
    const int num_reps = (int)md_array_size(data->representation.reps);
    m.rep_selected = CLAMP(m.rep_selected, 0, MAX(num_reps - 1, 0));
    m.rep_prop_selected = CLAMP(m.rep_prop_selected, 0, (int)RepProp::Count - 1);
    const Representation* rep = num_reps > 0 ? &data->representation.reps[m.rep_selected] : nullptr;
    float lo = 0.0f, hi = 1.0f;
    const char* label = rep ? movie_rep_prop_label(*rep, m.rep_prop_selected, &lo, &hi) : nullptr;
    if (rep && !label) {
        m.rep_prop_selected = (int)RepProp::Visible;
        label = movie_rep_prop_label(*rep, m.rep_prop_selected, &lo, &hi);
    }
    const int prop = m.rep_prop_selected;
    const bool visible_prop = prop == (int)RepProp::Visible;
    const bool color_prop = rep_prop_comps(prop) == 3;

    std::vector<int> mine;
    if (rep) {
        for (int i = 0; i < (int)m.rep_keys.size(); ++i) {
            if (m.rep_keys[i].rep == rep->id && m.rep_keys[i].prop == prop) mine.push_back(i);
        }
    }

    constexpr int N = 200;
    float xs[N], ys[N];
    if (!mine.empty() && !color_prop) {
        for (int i = 0; i < N; ++i) {
            const double t = (double)movie_len * (double)i / (double)(N - 1);
            float v[3] = {};
            movie_rep_eval(data, rep->id, prop, t, v);
            xs[i] = (float)t;
            ys[i] = v[0];
        }
    }

    const ImPlotDragToolFlags drag_flags = ImPlotDragToolFlags_NoFit | (locked ? ImPlotDragToolFlags_NoInputs : 0);
    bool any_held = false, any_hovered = false;
    int  remove_idx = -1;
    bool add_key = false;
    double add_time = 0.0, add_value = 0.0;

    char axis[96];
    if (rep && color_prop) axis[0] = '\0';
    else if (rep) snprintf(axis, sizeof(axis), "%s: %s", rep->name, label);
    else     snprintf(axis, sizeof(axis), "Representation");

    const ImPlotFlags plot_flags = ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMouseText | ImPlotFlags_NoTitle | ImPlotFlags_NoLegend;
    if (ImPlot::BeginPlot("##rep_lane", ImVec2(-1, -1), plot_flags)) {
        ImPlot::SetupAxes("Movie time (s)", color_prop ? nullptr : axis, 0, ImPlotAxisFlags_Lock | (color_prop ? ImPlotAxisFlags_NoDecorations : 0));
        movie_lane_axes_format();
        ImPlot::SetupAxisLinks(ImAxis_X1, &m.timeline_view_begin, &m.timeline_view_end);
        const double pad = 0.1 * (double)(hi - lo);
        ImPlot::SetupAxisLimits(ImAxis_Y1, (double)lo - pad, (double)hi + pad, ImPlotCond_Always);
        movie_lane_title(rep ? (std::string("Representation: ") + rep->name + ", " + label).c_str() : "Representation");
        if (rep && color_prop) ImPlot::PlotText(axis[0] ? axis : (std::string(rep->name) + ": " + label).c_str(), 0.02 * (double)movie_len, 0.95, ImVec2(0, 0));

        if (!rep) {
            ImPlot::PlotText("No representation yet", 0.5 * (double)movie_len, 0.5);
        }

        const ImVec4 col(0.95f, 0.6f, 0.3f, 1.0f);
        if (!mine.empty() && !color_prop) {
            ImPlot::PushStyleColor(ImPlotCol_Line, col);
            ImPlot::PlotLine(label, xs, ys, N);
            ImPlot::PopStyleColor();
        }

        KeyShift lane;
        const int64_t subject = rep ? rep_key_subject(rep->id, prop) : 0;
        if (rep && !color_prop && !visible_prop) {
            lane.lane = KeyLane::Rep;
            lane.subject = subject;
            lane.lo = lo;
            lane.hi = hi;
            lane.step = (double)(hi - lo) / 50.0;
        }
        for (int ki : mine) {
            const RepKey key = m.rep_keys[ki];
            const ImVec4 point_col = color_prop ? ImVec4(key.value[0], key.value[1], key.value[2], 1.0f) : col;
            const MovieKeyPointResult r = movie_key_point(data, locked, KeyKind::Rep, subject, key.time, color_prop ? 0.5 : (double)key.value[0], 5000 + ki, point_col, 7.0f, lane, &any_held, visible_prop);
            if (r.hovered) {
                any_hovered = true;
                if (!r.held) {
                    if (visible_prop) ImGui::SetTooltip("%.2f s: %s\nDrag sideways to change when. Right click to remove it", key.time, key.value[0] >= 0.5f ? "shown" : "hidden");
                    else if (color_prop) ImGui::SetTooltip("%.2f s\nDrag to change its time, right click to remove it. Its color is edited in the table of the Movie window.", key.time);
                    else              ImGui::SetTooltip("%.2f s, %.3g\nDrag to change it, right click to remove it", key.time, key.value[0]);
                }
                if (!locked && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) remove_idx = ki;
            }
        }

        if (rep) {
            double vline = (double)m.playhead;
            movie_lane_box(data, locked, 40, any_hovered, &vline, 1, [&](double x0, double x1, double y0, double y1) {
                for (int ki : mine) {
                    const RepKey& k = m.rep_keys[ki];
                    const double ky = color_prop ? 0.5 : (double)k.value[0];
                    if (k.time >= x0 && k.time <= x1 && ky >= y0 && ky <= y1) m.sel.add(KeyKind::Rep, subject, k.time);
                }
            });
        }

        if (rep && !locked && !any_hovered && ImPlot::IsPlotHovered() && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
            const ImPlotPoint p = ImPlot::GetPlotMousePos();
            add_key = true;
            add_time = movie_snap_time(data, p.x);
            add_value = CLAMP(p.y, (double)lo, (double)hi);
        }

        double playhead = (double)m.playhead;
        if (ImPlot::DragLineX(1000, &playhead, ImVec4(1, 1, 0, 1), 1.5f, drag_flags)) {
            m.playhead = (float)movie_snap_time(data, playhead);
            if (!locked) movie_apply_time(data, (double)m.playhead, true);
        }

        if (locked) {
            double cur = m.cur_time;
            ImPlot::DragLineX(1001, &cur, ImVec4(1.0f, 0.3f, 0.3f, 1), 1.5f, ImPlotDragToolFlags_NoInputs | ImPlotDragToolFlags_NoFit);
        }

        ImPlot::EndPlot();
    }

    if (any_held) movie_apply_time(data, (double)m.playhead, true);

    if (remove_idx >= 0) {
        m.rep_keys.erase(m.rep_keys.begin() + remove_idx);
        movie_reps_apply(data, (double)m.playhead);
    }
    if (add_key && rep) {
        RepKey key;
        key.rep = rep->id;
        key.prop = prop;
        key.time = add_time;
        if (visible_prop) {
            // A change of state at that time: what it is there now, turned around. With no keys yet, what it is now holds up to it.
            float cur[3] = {};
            movie_rep_prop_get(*rep, prop, cur);
            rep_keys_evaluate(cur, m.rep_keys.data(), m.rep_keys.size(), rep->id, prop, add_time);
            const float now = cur[0];
            if (mine.empty()) {
                RepKey first = key;
                first.time = 0.0;
                first.value[0] = now;
                first.ease = KeyEase::Hold;
                m.rep_keys.push_back(first);
            }
            key.value[0] = now >= 0.5f ? 0.0f : 1.0f;
            key.ease = KeyEase::Hold;
        } else if (color_prop) {
            // The colour it has at that time, or as it is now if there are no keys yet
            movie_rep_prop_get(*rep, prop, key.value);
            rep_keys_evaluate(key.value, m.rep_keys.data(), m.rep_keys.size(), rep->id, prop, add_time);
        } else {
            key.value[0] = (float)add_value;
        }
        m.rep_keys.push_back(key);
        movie_rep_sort(data);
        movie_apply_time(data, (double)m.playhead, true);
    }
}

// The overlays as bars, one row each: drag a bar for when it is shown, its ends for when it starts and stops. A bar that is
// moved takes the times of the subplots of a plot overlay with it. Fills the space that is left.
static void draw_movie_overlay_lane(ApplicationState* data, float movie_len, bool locked) {
    auto& m = data->movie;
    if (movie_len <= 0.0f) return;
    const int n = (int)m.overlays.size();
    const ImPlotDragToolFlags drag_flags = ImPlotDragToolFlags_NoFit | (locked ? ImPlotDragToolFlags_NoInputs : 0);
    const ImPlotFlags plot_flags = ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMouseText | ImPlotFlags_NoTitle | ImPlotFlags_NoLegend;
    static const ImVec4 type_colors[(int)MovieOverlayType::Count] = {
        ImVec4(0.90f, 0.90f, 0.90f, 1), ImVec4(0.65f, 0.65f, 1.00f, 1), ImVec4(0.55f, 0.90f, 0.55f, 1), ImVec4(1.00f, 0.80f, 0.40f, 1),
        ImVec4(1.00f, 0.60f, 0.40f, 1), ImVec4(0.80f, 0.50f, 1.00f, 1), ImVec4(0.40f, 0.80f, 1.00f, 1), ImVec4(0.40f, 1.00f, 0.80f, 1),
        ImVec4(1.00f, 0.50f, 0.70f, 1), ImVec4(0.80f, 0.80f, 0.80f, 1),
    };

    if (ImPlot::BeginPlot("##overlay_lane", ImVec2(-1, -1), plot_flags)) {
        ImPlot::SetupAxes("Movie time (s)", nullptr, 0, ImPlotAxisFlags_Lock | ImPlotAxisFlags_NoDecorations | ImPlotAxisFlags_Invert);
        movie_lane_axes_format();
        ImPlot::SetupAxisLinks(ImAxis_X1, &m.timeline_view_begin, &m.timeline_view_end);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -0.7, (double)MAX(n, 1) - 0.3, ImPlotCond_Always);
        movie_lane_title("Overlays");
        if (n == 0) ImPlot::PlotText("No overlay: add one in the Movie window", 0.5 * (double)movie_len, 0.0);
        bool any_hovered = false;

        for (int i = 0; i < n; ++i) {
            MovieOverlay& o = m.overlays[i];
            ImVec4 col = type_colors[CLAMP((int)o.type, 0, (int)MovieOverlayType::Count - 1)];
            if (!o.enabled) col.w = 0.35f;
            double x0 = o.begin, x1 = o.end, y0 = (double)i - 0.38, y1 = (double)i + 0.38;
            bool clicked = false, hovered = false, held = false;
            const bool changed = ImPlot::DragRect(7000 + i, &x0, &y0, &x1, &y1, col, drag_flags, &clicked, &hovered, &held);
            any_hovered |= hovered;
            if (changed && !locked) {
                double b = movie_snap_time(data, MIN(x0, x1)), e = movie_snap_time(data, MAX(x0, x1));
                const bool moved_only = fabs((x1 - x0) - (o.end - o.begin)) < 1.0e-6;
                b = CLAMP(b, 0.0, (double)movie_len);
                e = CLAMP(e, b + 0.05, (double)movie_len);
                if (moved_only) {
                    // The whole bar was moved, with the picked ones, and what is timed inside it moves too
                    movie_group_drag(data, KeyKind::Overlay, i, o.begin, o.end, 0.0, b, 0.0, KeyShift{});
                } else {
                    o.begin = b;
                    o.end = e;
                }
            }
            if (held && m.key_drag.active && m.key_drag.moved) movie_group_drag_tooltip(data, KeyShift{});
            if (m.sel.contains(KeyKind::Overlay, i, o.begin)) {
                ImPlot::PushPlotClipRect();
                ImPlot::GetPlotDrawList()->AddRect(ImPlot::PlotToPixels(o.begin, y0), ImPlot::PlotToPixels(o.end, y1), IM_COL32(255, 255, 255, 255), 0.0f, 0, 2.0f);
                ImPlot::PopPlotClipRect();
            }
            if (clicked && !locked) {
                movie_pick_click(data, KeyKind::Overlay, i, o.begin, o.end, KeyShift{});
                m.overlay_selected = i;
                m.editor_controls = true;
                m.editor_select_overlays = true;
            }
            if (hovered && !held) {
                ImGui::SetTooltip("%d  %s%s%s\n%.2f s to %.2f s\nDrag to move it, its ends to change when it starts and stops", i + 1,
                    movie_overlay_type_str[(int)o.type], o.type == MovieOverlayType::Text ? ": " : "", o.type == MovieOverlayType::Text ? o.text : "", o.begin, o.end);
            }

            char label[96];
            snprintf(label, sizeof(label), "%d %s%s%s", i + 1, movie_overlay_type_str[(int)o.type], o.type == MovieOverlayType::Text ? ": " : "", o.type == MovieOverlayType::Text ? o.text : "");
            ImPlot::PlotText(label, o.begin, (double)i, ImVec2(ImGui::CalcTextSize(label).x * 0.5f + 6.0f, 0));

            // When the subplots of a plot overlay come in
            if (o.type == MovieOverlayType::Timeline || o.type == MovieOverlayType::Distribution) {
                ImPlot::PushPlotClipRect();
                ImDrawList* dl = ImPlot::GetPlotDrawList();
                for (const MoviePlotPanel& panel : o.panels) {
                    if (panel.begin <= o.begin) continue;
                    dl->AddLine(ImPlot::PlotToPixels(panel.begin, (double)i - 0.38), ImPlot::PlotToPixels(panel.begin, (double)i + 0.38), ImGui::ColorConvertFloat4ToU32(col), 2.0f);
                }
                ImPlot::PopPlotClipRect();
            }
        }

        {
            const double vline = (double)m.playhead;
            movie_lane_box(data, locked, 50, any_hovered, &vline, 1, [&](double x0, double x1, double y0, double y1) {
                for (int i = 0; i < (int)m.overlays.size(); ++i) {
                    const MovieOverlay& o = m.overlays[i];
                    if (o.end >= x0 && o.begin <= x1 && (double)i + 0.38 >= y0 && (double)i - 0.38 <= y1) m.sel.add(KeyKind::Overlay, i, o.begin, o.end);
                }
            });
        }

        // The notes of the movie, at the bottom
        {
            ImPlot::PushPlotClipRect();
            ImDrawList* dl = ImPlot::GetPlotDrawList();
            for (size_t mi = 0; mi < m.markers.size(); ++mi) {
                const MovieMarker& mk = m.markers[mi];
                float mc[4];
                movie_marker_color(m.markers.data(), m.markers.size(), mi, mc);
                const ImVec2 p = ImPlot::PlotToPixels(mk.time, (double)n - 0.35);
                dl->AddTriangleFilled(ImVec2(p.x - 4.0f, p.y), ImVec2(p.x + 4.0f, p.y), ImVec2(p.x, p.y - 7.0f), ImGui::ColorConvertFloat4ToU32(ImVec4(mc[0], mc[1], mc[2], 0.9f)));
            }
            ImPlot::PopPlotClipRect();
        }

        double playhead = (double)m.playhead;
        if (ImPlot::DragLineX(1000, &playhead, ImVec4(1, 1, 0, 1), 1.5f, drag_flags)) {
            m.playhead = (float)movie_snap_time(data, playhead);
            if (!locked) movie_apply_time(data, (double)m.playhead, true);
        }
        if (locked) {
            double cur = m.cur_time;
            ImPlot::DragLineX(1001, &cur, ImVec4(1.0f, 0.3f, 0.3f, 1), 1.5f, ImPlotDragToolFlags_NoInputs | ImPlotDragToolFlags_NoFit);
        }
        ImPlot::EndPlot();
    }
}

// The representations that are in one group (named alike up to the first hyphen), by their place in the list
static std::vector<int> movie_rep_group_members(const ApplicationState* data, const std::string& group) {
    std::vector<int> out;
    for (int i = 0; i < (int)md_array_size(data->representation.reps); ++i) {
        std::string g, member;
        rep_name_split(data->representation.reps[i].name, &g, &member);
        if (g == group) out.push_back(i);
    }
    return out;
}

static std::vector<RepRow> movie_rep_overview_rows(const ApplicationState* data) {
    std::vector<std::string> names;
    for (size_t i = 0; i < md_array_size(data->representation.reps); ++i) names.push_back(data->representation.reps[i].name);
    return rep_system_rows(names);
}

static std::vector<RepBlock> movie_rep_system_blocks(const ApplicationState* data, const std::string& group, double duration) {
    std::vector<RepBlock> blocks;
    for (int member : movie_rep_group_members(data, group)) {
        const Representation& rep = data->representation.reps[member];
        const std::vector<RepInterval> intervals = rep_effective_intervals(data->movie.rep_keys, rep.id, duration, rep.enabled);
        for (int k = 0; k < (int)intervals.size(); ++k) blocks.push_back({member, k, intervals[k], 0});
    }
    return blocks;
}

// The bands of the systems, in plot units: one unit is a row of text
struct RepOverviewLayout {
    std::vector<RepRow> rows;
    std::vector<std::vector<RepBlock>> blocks;
    std::vector<int> slots;
    std::vector<double> band_begin, band_height;
    double total = 1.0;
};

static constexpr double REP_OVERVIEW_SYSTEM_GAP = 0.6;
static constexpr double REP_OVERVIEW_EQUAL_BAND = 1.5;

static RepOverviewLayout movie_rep_overview_layout(const ApplicationState* data, double duration) {
    RepOverviewLayout l;
    l.rows = movie_rep_overview_rows(data);
    double y = 0.0;
    for (const RepRow& row : l.rows) {
        l.blocks.push_back(movie_rep_system_blocks(data, row.group, duration));
        double h = 0.0;
        if (data->movie.timeline_rep_separate) {
            // A line for each representation of the system, in the order of the Representations window
            const std::vector<int> members = movie_rep_group_members(data, row.group);
            for (RepBlock& b : l.blocks.back()) {
                b.slot = (int)(std::find(members.begin(), members.end(), b.rep) - members.begin());
            }
            l.slots.push_back(MAX((int)members.size(), 1));
            h = (double)l.slots.back();
        } else {
            // Fades are not counted, so a block that starts where another ends cross-fades with it on the same line
            l.slots.push_back(rep_pack_blocks(&l.blocks.back(), 0.0));
            h = data->movie.timeline_rep_equal_rows ? REP_OVERVIEW_EQUAL_BAND : (double)l.slots.back();
        }
        h *= (double)CLAMP(data->movie.timeline_rep_block_height, 0.5f, 2.5f);
        l.band_begin.push_back(y);
        l.band_height.push_back(h);
        y += h + REP_OVERVIEW_SYSTEM_GAP;
    }
    if (!l.rows.empty()) l.total = y - REP_OVERVIEW_SYSTEM_GAP;
    return l;
}

static float movie_rep_overview_height(const ApplicationState* data, double duration) {
    return ImGui::GetFontSize() * ((float)movie_rep_overview_layout(data, duration).total * 1.8f + 4.0f);
}

// One row per system; overlapping representation blocks occupy separate slots within that row.
static void draw_movie_rep_overview_lane(ApplicationState* data, float movie_len, bool locked) {
    auto& m = data->movie;
    if (movie_len <= 0.0f) return;
    const double duration = (double)movie_len;
    const double ramp = MAX((double)m.rep_transition, 0.0);
    const RepOverviewLayout layout = movie_rep_overview_layout(data, duration);
    const std::vector<RepRow>& rows = layout.rows;
    const int n = (int)rows.size();
    const ImPlotDragToolFlags drag_flags = ImPlotDragToolFlags_NoFit | (locked ? ImPlotDragToolFlags_NoInputs : 0);
    const ImPlotFlags plot_flags = ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMouseText | ImPlotFlags_NoTitle | ImPlotFlags_NoLegend;

    std::vector<const char*> label_ptr((size_t)n);
    std::vector<double> label_pos((size_t)n);
    for (int r = 0; r < n; ++r) {
        label_ptr[(size_t)r] = rows[(size_t)r].label.c_str();
        label_pos[(size_t)r] = layout.band_begin[(size_t)r] + 0.5 * layout.band_height[(size_t)r];
    }
    const double half_gap = 0.5 * REP_OVERVIEW_SYSTEM_GAP;

    // What was dragged, done once the lane is drawn: the keys it edits are read while it is
    struct Edit { int kind = 0; uint32_t rep = 0, target = 0; RepInterval iv; double b = 0.0, e = 0.0; } edit;
    static std::string menu_group;
    static uint32_t menu_rep = 0;
    static double menu_begin = 0.0, menu_end = 0.0;
    static bool menu_add = false;
    bool any_hovered = false;

    if (ImPlot::BeginPlot("##rep_overview", ImVec2(-1, -1), plot_flags)) {
        ImPlot::SetupAxes("Movie time (s)", nullptr, 0, ImPlotAxisFlags_Lock | ImPlotAxisFlags_Invert);
        movie_lane_axes_format();
        ImPlot::SetupAxisLinks(ImAxis_X1, &m.timeline_view_begin, &m.timeline_view_end);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -half_gap, layout.total + half_gap, ImPlotCond_Always);
        if (n > 0) ImPlot::SetupAxisTicks(ImAxis_Y1, label_pos.data(), n, label_ptr.data());
        ImPlot::SetupFinish();
        movie_lane_title("Representations");
        if (n == 0) ImPlot::PlotText("No representation yet", 0.5 * duration, 0.5);

        ImDrawList* dl = ImPlot::GetPlotDrawList();
        for (int r = 0; r < n; ++r) {
            const RepRow& row = rows[(size_t)r];
            const std::vector<int> members = movie_rep_group_members(data, row.group);
            const std::vector<RepBlock>& blocks = layout.blocks[(size_t)r];
            const bool selected = std::find(members.begin(), members.end(), m.rep_selected) != members.end();
            const double band_top = layout.band_begin[(size_t)r];
            const double band_h = layout.band_height[(size_t)r];
            const double slot_h = band_h / (double)layout.slots[(size_t)r];
            const double pad = MIN(0.05, 0.1 * slot_h);

            // Each system has its band, shaded every other one and divided from the one above
            {
                const ImVec2 a = ImPlot::PlotToPixels(m.timeline_view_begin, band_top - half_gap * 0.5);
                const ImVec2 b = ImPlot::PlotToPixels(m.timeline_view_end, band_top + band_h + half_gap * 0.5);
                ImPlot::PushPlotClipRect();
                dl->AddRectFilled(a, b, selected ? IM_COL32(255, 255, 255, 26) : (r % 2 ? IM_COL32(255, 255, 255, 12) : IM_COL32(0, 0, 0, 30)));
                if (r > 0) {
                    const ImVec2 la = ImPlot::PlotToPixels(m.timeline_view_begin, band_top - half_gap);
                    const ImVec2 lb = ImPlot::PlotToPixels(m.timeline_view_end, band_top - half_gap);
                    dl->AddLine(la, lb, IM_COL32(255, 255, 255, 60));
                }
                ImPlot::PopPlotClipRect();
            }

            for (const RepBlock& block : blocks) {
                const Representation* rep = &data->representation.reps[block.rep];
                const RepInterval& iv = block.interval;
                std::string group, member;
                rep_name_split(rep->name, &group, &member);
                const char* label = member.empty() ? representation_type_str[(int)rep->type] : member.c_str();
                const float hue = (float)((int)rep->type * 137 % 360) / 360.0f;
                ImVec4 col(1, 1, 1, 1);
                ImGui::ColorConvertHSVtoRGB(hue, 0.6f, 0.95f, col.x, col.y, col.z);
                const double top = band_top + block.slot * slot_h + pad;
                const double bottom = top + slot_h - 2.0 * pad;
                double x0 = iv.begin, x1 = iv.end, y0 = top, y1 = bottom;
                bool clicked = false, hovered = false, held = false;
                ImGui::PushID((int)rep->id);
                const bool changed = ImPlot::DragRect(block.interval_index, &x0, &y0, &x1, &y1, col, drag_flags, &clicked, &hovered, &held);
                ImGui::PopID();
                any_hovered |= hovered;

                // The bar: it grows in over the transition after its start and shrinks away over it after its end
                const double head = MIN(ramp, iv.end - iv.begin);
                const double tail = MIN(ramp, MAX(duration - iv.end, 0.0));
                const ImU32 strong = ImGui::ColorConvertFloat4ToU32(ImVec4(col.x, col.y, col.z, 0.8f));
                const ImU32 faint = ImGui::ColorConvertFloat4ToU32(ImVec4(col.x, col.y, col.z, 0.12f));
                const ImU32 none = ImGui::ColorConvertFloat4ToU32(ImVec4(col.x, col.y, col.z, 0.0f));
                ImPlot::PushPlotClipRect();
                auto px = [&](double t, double y) { return ImPlot::PlotToPixels(t, y); };
                if (head > 0.0) dl->AddRectFilledMultiColor(px(iv.begin, top), px(iv.begin + head, bottom), faint, strong, strong, faint);
                dl->AddRectFilled(px(iv.begin + head, top), px(iv.end, bottom), strong);
                if (tail > 0.0) dl->AddRectFilledMultiColor(px(iv.end, top), px(iv.end + tail, bottom), strong, none, none, strong);
                const ImVec2 a = px(iv.begin, top), b = px(iv.end, bottom);
                dl->PushClipRect(a, b, true);
                dl->AddText(ImVec2(a.x + 4.0f, a.y + MAX((b.y - a.y - ImGui::GetFontSize()) * 0.5f, 0.0f)),
                    IM_COL32(255, 255, 255, 255), label);
                dl->PopClipRect();
                if (m.sel.contains(KeyKind::Block, rep->id, iv.begin)) dl->AddRect(a, b, IM_COL32(255, 255, 255, 255), 0.0f, 0, 2.0f);
                ImPlot::PopPlotClipRect();

                if (hovered && !held) {
                    ImGui::SetTooltip("%s (%s)\n%.2f s to %.2f s\nDrag to move or resize. Right click to change representation or remove.\nClick picks it, Ctrl + click adds it to the picked items, which move together.", rep->name,
                        representation_type_str[(int)rep->type], iv.begin, iv.end);
                }
                if (held && m.key_drag.active && m.key_drag.moved) movie_group_drag_tooltip(data, KeyShift{});
                if (clicked) {
                    m.rep_selected = block.rep;
                    m.rep_prop_selected = (int)RepProp::Visible;
                    if (!locked && iv.begin_key >= 0) movie_pick_click(data, KeyKind::Block, rep->id, iv.begin, iv.end, KeyShift{});
                }
                if (!locked && hovered && !held && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
                    menu_group = row.group;
                    menu_rep = rep->id;
                    menu_begin = iv.begin;
                    menu_end = iv.end;
                    menu_add = false;
                    ImGui::OpenPopup("System representation block");
                }
                if (changed && edit.kind == 0) {
                    const double length = iv.end - iv.begin;
                    const bool moved_only = fabs((x1 - x0) - length) < 1.0e-6;
                    double nb = movie_snap_time(data, MIN(x0, x1)), ne = movie_snap_time(data, MAX(x0, x1));
                    if (moved_only) ne = nb + length;
                    if (moved_only && iv.begin_key >= 0 && !locked) {
                        // A block that is moved goes with the other picked items
                        movie_group_drag(data, KeyKind::Block, rep->id, iv.begin, iv.end, 0.0, nb, 0.0, KeyShift{});
                    } else {
                        edit.kind = 1;
                        edit.iv = iv;
                        edit.rep = rep->id;
                        edit.b = nb;
                        edit.e = ne;
                    }
                }
                ImPlot::PushPlotClipRect();
                for (const RepKey& key : m.rep_keys) {
                    if (key.rep != rep->id || key.prop == (int)RepProp::Visible) continue;
                    if (key.time < iv.begin || key.time > iv.end) continue;
                    const ImVec2 p = ImPlot::PlotToPixels(key.time, bottom);
                    dl->AddTriangleFilled(ImVec2(p.x - 3.5f, p.y), ImVec2(p.x + 3.5f, p.y), ImVec2(p.x, p.y - 6.0f), IM_COL32(255, 255, 255, 220));
                }
                ImPlot::PopPlotClipRect();
            }
        }

        {
            const double vline = (double)m.playhead;
            movie_lane_box(data, locked, 60, any_hovered, &vline, 1, [&](double x0, double x1, double y0, double y1) {
                for (int r = 0; r < n; ++r) {
                    const double slot_h = layout.band_height[(size_t)r] / (double)layout.slots[(size_t)r];
                    const double pad = MIN(0.05, 0.1 * slot_h);
                    for (const RepBlock& block : layout.blocks[(size_t)r]) {
                        const double top = layout.band_begin[(size_t)r] + block.slot * slot_h + pad, bottom = top + slot_h - 2.0 * pad;
                        const RepInterval& iv = block.interval;
                        if (iv.begin_key >= 0 && iv.end >= x0 && iv.begin <= x1 && bottom >= y0 && top <= y1) {
                            m.sel.add(KeyKind::Block, data->representation.reps[block.rep].id, iv.begin, iv.end);
                        }
                    }
                }
            });
        }

        // A stretch is added where a representation is hidden
        if (!locked && !any_hovered && edit.kind == 0 && ImPlot::IsPlotHovered() && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
            const ImPlotPoint p = ImPlot::GetPlotMousePos();
            int r = -1;
            for (int i = 0; i < n; ++i) {
                const double b0 = layout.band_begin[(size_t)i] - half_gap, b1 = layout.band_begin[(size_t)i] + layout.band_height[(size_t)i] + half_gap;
                if (p.y >= b0 && p.y < b1) { r = i; break; }
            }
            if (r >= 0 && r < n) {
                menu_group = rows[r].group;
                menu_begin = movie_snap_time(data, p.x);
                menu_end = MIN(menu_begin + MAX(0.1 * duration, 1.0), duration);
                menu_add = true;
                ImGui::OpenPopup("System representation block");
            }
        }

        if (ImGui::BeginPopup("System representation block")) {
            ImGui::BeginDisabled(locked);
            RepInterval source;
            bool found = menu_add;
            bool enabled = false;
            for (int member : movie_rep_group_members(data, menu_group)) {
                if (data->representation.reps[member].id == menu_rep) enabled = data->representation.reps[member].enabled;
            }
            for (const RepInterval& iv : rep_effective_intervals(m.rep_keys, menu_rep, duration, enabled)) {
                if (fabs(iv.begin - menu_begin) < 1e-6 && fabs(iv.end - menu_end) < 1e-6) {
                    source = iv;
                    found = true;
                    break;
                }
            }
            ImGui::TextUnformatted(menu_add ? "Add representation block" : "Change representation");
            for (int member : movie_rep_group_members(data, menu_group)) {
                const Representation& target = data->representation.reps[member];
                ImGui::PushID((int)target.id);
                const auto spans = rep_effective_intervals(m.rep_keys, target.id, duration, target.enabled);
                bool available = menu_end - menu_begin >= 0.05;
                for (const RepInterval& span : spans) {
                    if (menu_begin >= span.begin - 0.05 && menu_begin < span.end + 0.05) available = false;
                }
                if (ImGui::MenuItem(target.name, nullptr, !menu_add && menu_rep == target.id, found && (!menu_add || available))) {
                    edit.kind = menu_add ? 4 : 5;
                    edit.rep = menu_add ? target.id : menu_rep;
                    edit.target = target.id;
                    edit.iv = source;
                    edit.b = menu_begin;
                    edit.e = menu_end;
                    m.rep_selected = member;
                    m.rep_prop_selected = (int)RepProp::Visible;
                }
                ImGui::PopID();
            }
            if (!menu_add && ImGui::MenuItem("Remove block", nullptr, false, found)) {
                edit.kind = 3;
                edit.iv = source;
                edit.rep = menu_rep;
            }
            ImGui::EndDisabled();
            ImGui::EndPopup();
        }
        double playhead = (double)m.playhead;
        if (ImPlot::DragLineX(1000, &playhead, ImVec4(1, 1, 0, 1), 1.5f, drag_flags)) {
            m.playhead = (float)movie_snap_time(data, playhead);
            if (!locked) movie_apply_time(data, (double)m.playhead, true);
        }
        if (locked) {
            double cur = m.cur_time;
            ImPlot::DragLineX(1001, &cur, ImVec4(1.0f, 0.3f, 0.3f, 1), 1.5f, ImPlotDragToolFlags_NoInputs | ImPlotDragToolFlags_NoFit);
        }
        ImPlot::EndPlot();
    }

    if (edit.kind != 0 && !locked) {
        if (edit.kind != 4 && edit.iv.begin_key < 0) {
            RepKey key;
            key.rep = edit.rep;
            key.prop = (int)RepProp::Visible;
            key.value[0] = 1.0f;
            key.ease = KeyEase::Hold;
            m.rep_keys.push_back(key);
            edit.iv = rep_shown_intervals(m.rep_keys, edit.rep, duration).front();
        }
        if (edit.kind == 5) {
            for (int member : movie_rep_group_members(data, menu_group)) {
                const Representation& rep = data->representation.reps[member];
                if (rep.id != edit.target || !rep.enabled) continue;
                const auto intervals = rep_effective_intervals(m.rep_keys, rep.id, duration, true);
                if (!intervals.empty() && intervals.front().begin_key < 0) {
                    RepKey key;
                    key.rep = rep.id;
                    key.prop = (int)RepProp::Visible;
                    key.value[0] = 1.0f;
                    key.ease = KeyEase::Hold;
                    m.rep_keys.push_back(key);
                }
            }
        }
        switch (edit.kind) {
        case 1: rep_move_interval(&m.rep_keys, edit.rep, edit.iv, edit.b, edit.e, duration); break;
        case 3: rep_remove_interval(&m.rep_keys, edit.iv); break;
        case 4:
            if (!rep_add_interval(&m.rep_keys, edit.rep, edit.b, edit.e, duration)) {
                MD_LOG_ERROR("Cannot add a representation block here: this representation is already shown or there is not enough room.");
            }
            break;
        case 5: rep_transfer_interval(&m.rep_keys, edit.rep, edit.target, edit.iv, duration); break;
        default: break;
        }
        movie_rep_sort(data);
        movie_reps_apply(data, (double)m.playhead);
    }
}

static void draw_movie_preview_controls(ApplicationState* data) {
    auto& m = data->movie;
    const bool recording = m.state == MovieRecordingState::Recording;

    const float movie_len = (float)movie_duration(data);
    m.playhead = CLAMP(m.playhead, 0.0f, movie_len);
    const float fs = ImGui::GetFontSize();

    ImGui::BeginDisabled(recording);
    if (ImGui::Button(m.preview_playing ? "Pause Preview" : "Play Preview")) {
        m.preview_playing = !m.preview_playing;
        if (m.preview_playing && m.playhead >= movie_len) m.playhead = 0.0f;
    }
    ImGui::SetItemTooltip("Plays the movie in the viewport at the speed it will have, from the preview time, without recording. Space");
    ImGui::SameLine();
    if (ImGui::Button("Play from start")) {
        m.playhead = 0.0f;
        movie_apply_time(data, 0.0, true);
        m.preview_playing = true;
    }
    ImGui::SetItemTooltip("Plays the movie from its beginning.");
    ImGui::SameLine();
    ImGui::Checkbox("Repeat", &m.preview_loop);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(CLAMP(ImGui::GetContentRegionAvail().x, fs * 8.0f, fs * 18.0f));
    if (ImGui::SliderFloat("##timeline_playhead", &m.playhead, 0.0f, movie_len, "%.2f s")) {
        movie_apply_time(data, (double)m.playhead, true);
    }
    ImGui::SetItemTooltip("Scrub the movie: shows the trajectory frame, the camera (in Movie preview) and the keyed looks.");
    ImGui::EndDisabled();
}

// Which lanes are shown, as a row of toggles with presets, and the layout settings in a popup
static void draw_movie_timeline_options(ApplicationState* data) {
    auto& m = data->movie;
    const bool recording = m.state == MovieRecordingState::Recording;
    const float movie_len = (float)movie_duration(data);
    const float fs = ImGui::GetFontSize();
    const int num_params = (int)(sizeof(movie_param_table) / sizeof(movie_param_table[0]));
    m.param_selected = CLAMP(m.param_selected, 0, num_params - 1);

    bool lens = m.timeline_tracks[1] || m.timeline_tracks[2];
    struct Chip { const char* label; bool* on; const char* tip; };
    const Chip chips[] = {
        {"Ruler",      &m.timeline_ruler,        "Time ruler: a tick per frame when zoomed in, the notes of the movie; click or drag to move the preview time."},
        {"Camera",     &m.timeline_camera_lane,  "Camera keys with their names, follow / look-at / spin bands and frame pins. Drag a key to change when,\nclick to go there, right click to rename or remove, double click to add a key on the path."},
        {"Trajectory", &m.timeline_tracks[0],    "The trajectory frame over the movie, with its blue Start and End and the keys that pin a frame."},
        {"Lens",       &lens,                    "Field of view (left axis) and camera distance (right axis) over the movie, with the camera keys on both."},
        {"Representations", &m.timeline_rep_overview, "When each representation is shown: a band per system with a block for each stretch. Drag to move or resize,\nright click to switch or remove, double click on empty space to add. Lines and block height are in Layout..."},
        {"Overlays",   &m.timeline_overlay_lane, "The overlays as bars, one row each: drag a bar to move it, its ends to change when it is shown."},
        {"Rep keys",   &m.timeline_rep_lane,     "The value of one property (color, radius, ...) of one representation over the movie, with its keys\n(the property and representation are picked next to the lanes)."},
        {"Look",       &m.timeline_param_lane,   "The keys of one look parameter (picked next to the lanes)."},
    };
    ImGui::AlignTextToFramePadding();
    ImGui::TextDisabled("Lanes");
    for (const Chip& c : chips) {
        ImGui::SameLine(0.0f, fs * 0.25f);
        const bool on = *c.on;
        if (on) {
            ImGui::PushStyleColor(ImGuiCol_Button, ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));
        }
        if (ImGui::SmallButton(c.label)) *c.on = !on;
        if (on) ImGui::PopStyleColor(2);
        ImGui::SetItemTooltip("%s", c.tip);
    }
    m.timeline_tracks[1] = m.timeline_tracks[2] = lens;

    auto preset = [&](bool ruler, bool camera, bool traj, bool lens_on, bool look, bool rep, bool overview, bool overlays) {
        m.timeline_ruler = ruler; m.timeline_camera_lane = camera; m.timeline_tracks[0] = traj;
        m.timeline_tracks[1] = m.timeline_tracks[2] = lens_on; m.timeline_param_lane = look;
        m.timeline_rep_lane = rep; m.timeline_rep_overview = overview; m.timeline_overlay_lane = overlays;
    };
    // The presets and the layout on a line of their own, under the lanes
    ImGui::Dummy(ImVec2(ImGui::CalcTextSize("Lanes").x, 0.0f));
    ImGui::SameLine(0.0f, fs * 0.25f);
    if (ImGui::SmallButton("All")) preset(true, true, true, true, true, false, true, false);
    ImGui::SetItemTooltip("Ruler, camera, trajectory, lens, representations and look parameter.");
    ImGui::SameLine(0.0f, fs * 0.25f);
    if (ImGui::SmallButton("Camera work")) preset(true, true, true, true, false, false, false, false);
    ImGui::SetItemTooltip("Ruler, camera, trajectory and lens.");
    ImGui::SameLine(0.0f, fs * 0.25f);
    if (ImGui::SmallButton("Scene work")) {
        preset(true, false, true, false, true, false, true, true);
        for (int i = 0; i < num_params; ++i)
            if (movie_param_table[i].id == MOVIE_PARAM_DOF_APERTURE) m.param_selected = i;
    }
    ImGui::SetItemTooltip("Ruler, trajectory, representations, overlays and depth of field blur.");
    ImGui::SameLine(0.0f, fs);
    if (ImGui::SmallButton("Layout...")) ImGui::OpenPopup("movie_lane_layout");
    ImGui::SetItemTooltip("Lane height, fit to window, snapping, the representations lane's lines and block height.");
    if (ImGui::BeginPopup("movie_lane_layout")) {
        ImGui::SetNextItemWidth(fs * 9.0f);
        ImGui::SliderFloat("Lane height", &m.timeline_lane_height, 60.0f, 400.0f, "%.0f px");
        ImGui::SetItemTooltip("The height of a lane at ratio 1, in pixels. Lanes keep their height and the window scrolls (the mouse wheel scrolls, Ctrl + wheel zooms the time). Lanes listing systems or overlays grow to fit them.");
        ImGui::Checkbox("Fit to window", &m.timeline_fit_window);
        ImGui::SetItemTooltip("The lanes share the height of the window, however small, with the mouse wheel zooming the time.");
        ImGui::Checkbox("Snap to frames", &m.snap_frames);
        ImGui::SetItemTooltip("Keys, the playhead and the trajectory's start and end that are dragged here land on a frame of the movie (at the Output FPS).");
        if (ImGui::Button("Show whole movie")) {
            m.timeline_view_begin = -0.03 * movie_len;
            m.timeline_view_end = 1.03 * movie_len;
        }
        ImGui::SeparatorText("Representations lane");
        {
            int mode = m.timeline_rep_separate ? 2 : (m.timeline_rep_equal_rows ? 0 : 1);
            ImGui::TextUnformatted("Lines");
            ImGui::SameLine();
            bool changed = ImGui::RadioButton("Compact", &mode, 0);
            ImGui::SetItemTooltip("One line per system: blocks that overlap in time share its height.");
            ImGui::SameLine();
            changed |= ImGui::RadioButton("Stacked", &mode, 1);
            ImGui::SetItemTooltip("A system gets another line only where its blocks overlap in time.");
            ImGui::SameLine();
            changed |= ImGui::RadioButton("Separate", &mode, 2);
            ImGui::SetItemTooltip("A line for each representation, grouped by system.");
            if (changed) {
                m.timeline_rep_separate = mode == 2;
                if (mode != 2) m.timeline_rep_equal_rows = mode == 0;
            }
            ImGui::BeginDisabled(m.timeline_fit_window);
            ImGui::SetNextItemWidth(fs * 9.0f);
            ImGui::SliderFloat("Block height", &m.timeline_rep_block_height, 0.5f, 2.5f, "x%.2f");
            ImGui::EndDisabled();
            ImGui::SetItemTooltip(m.timeline_fit_window
                ? "With Fit to window the lanes share the window's height; untick it to set the block height."
                : "How tall the blocks are, times the usual. The lane grows or shrinks with it.");
        }
        const int overview_reps = (int)md_array_size(data->representation.reps);
        std::vector<int> siblings;
        if (overview_reps > 0) {
            m.rep_selected = CLAMP(m.rep_selected, 0, overview_reps - 1);
            std::string g, member;
            rep_name_split(data->representation.reps[m.rep_selected].name, &g, &member);
            siblings = movie_rep_group_members(data, g);
        }
        ImGui::BeginDisabled(siblings.size() < 2 || recording);
        if (ImGui::SmallButton("Swap with the next of its group")) {
            const size_t at = (size_t)(std::find(siblings.begin(), siblings.end(), m.rep_selected) - siblings.begin());
            const int next = siblings[(at + 1) % siblings.size()];
            rep_swap_at(&m.rep_keys, data->representation.reps[m.rep_selected].id, data->representation.reps[next].id, movie_snap_time(data, (double)m.playhead));
            movie_rep_sort(data);
            movie_reps_apply(data, (double)m.playhead);
        }
        ImGui::EndDisabled();
        ImGui::SetItemTooltip("At the preview time the selected representation goes and the next one of its group (named alike up to the first hyphen) comes,\nwith the transition: e.g. protein-cartoon shrinks away while protein-cpk grows in.");
        ImGui::TextDisabled("Drag between lanes to resize them.");
        ImGui::EndPopup();
    }

    // What the look parameter and representation lanes show
    const int num_reps = (int)md_array_size(data->representation.reps);
    const bool rep_pickers = m.timeline_rep_lane && num_reps > 0;
    if (rep_pickers) {
        m.rep_selected = CLAMP(m.rep_selected, 0, num_reps - 1);
        const Representation& rep = data->representation.reps[m.rep_selected];
        float lo, hi;
        m.rep_prop_selected = CLAMP(m.rep_prop_selected, 0, (int)RepProp::Count - 1);
        if (!movie_rep_prop_label(rep, m.rep_prop_selected, &lo, &hi)) m.rep_prop_selected = (int)RepProp::Visible;
        ImGui::SetNextItemWidth(fs * 10.0f);
        if (ImGui::BeginCombo("##timeline_rep", rep.name)) {
            for (int i = 0; i < num_reps; ++i) {
                ImGui::PushID(i);
                if (ImGui::Selectable(data->representation.reps[i].name, i == m.rep_selected)) m.rep_selected = i;
                ImGui::PopID();
            }
            ImGui::EndCombo();
        }
        ImGui::SetItemTooltip("The representation in the Rep keys lane.");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(fs * 9.0f);
        if (ImGui::BeginCombo("##timeline_rep_prop", movie_rep_prop_label(rep, m.rep_prop_selected, &lo, &hi))) {
            for (int p = 0; p < (int)RepProp::Count; ++p) {
                if (const char* label = movie_rep_prop_label(rep, p, &lo, &hi)) {
                    if (ImGui::Selectable(label, p == m.rep_prop_selected)) m.rep_prop_selected = p;
                }
            }
            ImGui::EndCombo();
        }
        ImGui::SetItemTooltip("Its property in the Rep keys lane.");
    }
    if (m.timeline_param_lane) {
        if (rep_pickers) ImGui::SameLine(0.0f, fs);
        ImGui::SetNextItemWidth(fs * 12.0f);
        if (ImGui::BeginCombo("##look_param", movie_param_table[m.param_selected].label)) {
            for (int i = 0; i < num_params; ++i) {
                if (ImGui::Selectable(movie_param_table[i].label, i == m.param_selected)) m.param_selected = i;
            }
            ImGui::EndCombo();
        }
        ImGui::SetItemTooltip("The look parameter in the Look lane.");
    }
}

static void draw_movie_timeline_panel(ApplicationState* data) {
    auto& m = data->movie;
    const bool recording = m.state == MovieRecordingState::Recording;
    const float movie_len = (float)movie_duration(data);
    const float fs = ImGui::GetFontSize();
    if (movie_len <= 0.0f) {
        ImGui::TextDisabled("The movie has no duration yet.");
        return;
    }
    if (m.timeline_view_duration != movie_len) {
        m.timeline_view_begin = -0.03 * movie_len;
        m.timeline_view_end = 1.03 * movie_len;
        m.timeline_view_duration = movie_len;
    }
    draw_movie_timeline_options(data);
    {
        // The picked keys: what can be done with them is also on the keyboard, with the mouse over the lanes
        const size_t picked = m.sel.size();
        ImGui::BeginDisabled(recording);
        if (ImGui::SmallButton("Pick all")) movie_selection_select_all(data);
        ImGui::SetItemTooltip("Picks the camera keys and the keys of the look parameter and of the representation property that are shown. Ctrl + A.");
        ImGui::BeginDisabled(picked == 0);
        ImGui::SameLine();
        if (ImGui::SmallButton("Put down")) m.sel.clear();
        ImGui::SetItemTooltip("Esc");
        ImGui::SameLine();
        if (ImGui::SmallButton("Delete")) movie_selection_edit(data, 0.0, 0.0, true, false);
        ImGui::SetItemTooltip("Removes the picked keys. Delete.");
        ImGui::SameLine();
        if (ImGui::SmallButton("Copy")) movie_selection_copy(data);
        ImGui::SetItemTooltip("Remembers the picked keys, with the time between them. Ctrl + C.");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(fs * 5.5f);
        ImGui::DragFloat("##sel_scale", &m.sel_scale, 0.5f, 5.0f, 1000.0f, "%.0f %%");
        ImGui::SetItemTooltip("Stretches the time between the picked items (keys, bars, blocks, and what is timed in a bar): drag to the right to spread them out,\nto the left to bring them closer. 100 %% is as it is. It stops where something would leave the movie. Ctrl + click to type a value.");
        if (ImGui::IsItemActivated()) {
            auto& ks = m.key_scale;
            ks.active = true;
            ks.start = movie_keys_snapshot(data);
            ks.start_sel = m.sel;
            double t0 = 0.0, t1 = 0.0;
            key_selection_extent(ks.start, ks.start_sel, &t0, &t1);
            ks.anchor = m.sel_anchor == 1 ? (double)m.playhead : t0;
        }
        if (m.key_scale.active && ImGui::IsItemActive()) {
            auto& ks = m.key_scale;
            MovieKeys out = ks.start;
            movie_keys_scale(&out, &m.sel, ks.start, ks.start_sel, ks.anchor, (double)m.sel_scale / 100.0, (double)movie_duration(data));
            movie_keys_restore(data, out);
            movie_reps_apply(data, (double)m.playhead);
        }
        if (m.key_scale.active && ImGui::IsItemDeactivated()) {
            MovieKeys cur = movie_keys_snapshot(data);
            movie_keys_resolve(&cur, m.sel);
            movie_keys_restore(data, cur);
            key_selection_prune(&m.sel, cur);
            m.key_scale = decltype(m.key_scale){};
            m.sel_scale = 100.0f;
        }
        ImGui::SameLine();
        ImGui::SetNextItemWidth(fs * 8.0f);
        ImGui::Combo("##sel_anchor", &m.sel_anchor, "about first\0about preview time\0");
        ImGui::SetItemTooltip("What the stretch keeps in place: the first of the picked items, or the preview time.");
        ImGui::EndDisabled();
        ImGui::BeginDisabled(m.key_clip.empty());
        ImGui::SameLine();
        if (ImGui::SmallButton("Paste")) movie_selection_edit(data, 0.0, 0.0, false, true);
        ImGui::SetItemTooltip("Puts the copied keys in with the first one at the preview time. A key that lands on another replaces it. Ctrl + V.");
        ImGui::EndDisabled();
        ImGui::EndDisabled();
        ImGui::SameLine();
        if (picked > 0) ImGui::TextDisabled("%zu key%s picked: drag one to move them all, arrow keys nudge (Shift: 1 s, Ctrl: 10 frames; up and down change the value)", picked, picked == 1 ? "" : "s");
        else            ImGui::TextDisabled("Click a key to pick it, Ctrl + click adds, drag the background to pick with a box. Middle button or Shift + wheel pans.");
    }
    if (!m.sel.empty() && !recording) {
        // What the picked items have in common, to set for all of them
        double first = 1.0e30;
        int counts[5] = {};
        for (const KeyId& id : m.sel.ids) {
            first = MIN(first, id.time);
            counts[CLAMP((int)id.kind, 0, 4)] += 1;
        }
        static const char* kind_names[5] = {"camera", "look parameter", "representation", "overlay", "block"};
        char summary[160] = "";
        size_t len = 0;
        for (int k = 0; k < 5; ++k) {
            if (counts[k] == 0 || len + 40 >= sizeof(summary)) continue;
            len += (size_t)snprintf(summary + len, sizeof(summary) - len, "%s%d %s", len ? ", " : "", counts[k], kind_names[k]);
        }
        ImGui::AlignTextToFramePadding();
        ImGui::TextDisabled("%s: first at", summary);
        ImGui::SameLine();
        ImGui::SetNextItemWidth(fs * 6.0f);
        double start_time = first;
        if (ImGui::InputDouble("##sel_start", &start_time, 0.0, 0.0, "%.2f s", ImGuiInputTextFlags_EnterReturnsTrue)) movie_selection_edit(data, start_time - first, 0.0, false, false);
        ImGui::SetItemTooltip("When the first picked item is. Type a time (and Enter) to move all the picked items so that the first one is there.");

        KeyEase common = KeyEase::Smooth;
        bool mixed = false;
        int with_ease = key_selection_ease(m.keyframes, md_array_size(m.keyframes), m.param_keys, m.rep_keys, m.sel, &common, &mixed);
        if (with_ease > 0) {
            ImGui::SameLine();
            ImGui::SetNextItemWidth(fs * 8.0f);
            if (ImGui::BeginCombo("##sel_ease", mixed ? "(mixed)" : key_ease_str[(int)common])) {
                for (int e = 0; e < (int)KeyEase::Count; ++e) {
                    if (ImGui::Selectable(key_ease_str[e], !mixed && e == (int)common)) {
                        MovieKeys keys = movie_keys_snapshot(data);
                        movie_keys_set_ease(&keys, m.sel, (KeyEase)e);
                        movie_keys_restore(data, keys);
                        movie_params_apply(data, (double)m.playhead);
                        movie_reps_apply(data, (double)m.playhead);
                    }
                }
                ImGui::EndCombo();
            }
            ImGui::SetItemTooltip("How the movie gets to the picked keys (%d): they all take what you choose. The first camera key and the visibility keys have no easing to choose.", with_ease);
        }
        if (m.sel.size() == 1 && m.sel.ids[0].kind == KeyKind::Camera) {
            for (size_t i = 0; i < md_array_size(m.keyframes); ++i) {
                if (fabs(m.keyframes[i].time - m.sel.ids[0].time) >= 1.0e-9) continue;
                ImGui::SameLine();
                ImGui::SetNextItemWidth(fs * 9.0f);
                ImGui::InputTextWithHint("##sel_name", "name", m.keyframes[i].name, sizeof(m.keyframes[i].name));
                ImGui::SetItemTooltip("What the camera key is called, shown in the camera lane and the table.");
                break;
            }
        }
    }
    const ImVec2 strip_size(-1, MAX(ImGui::GetContentRegionAvail().y, fs * 12.0f));
    const ImVec2 strip_pos = ImGui::GetCursorScreenPos();
    const ImVec2 strip_avail = ImGui::GetContentRegionAvail();
    if (m.timeline_fit_window) {
        draw_movie_strip(data, movie_len, recording, strip_size);
    } else {
        // The lanes have their height and the window scrolls, so the wheel is for scrolling and Ctrl + wheel zooms the time
        const ImPlotInputMap old_map = ImPlot::GetInputMap();
        ImPlot::GetInputMap().ZoomMod = ImGuiMod_Ctrl;
        if (ImGui::BeginChild("##movie_strip_scroll", strip_size, ImGuiChildFlags_None, ImGuiWindowFlags_None)) {
            draw_movie_strip(data, movie_len, recording, ImVec2(-1, -1));
        }
        ImGui::EndChild();
        ImPlot::GetInputMap() = old_map;
    }
    movie_selection_shortcuts(data, ImGui::IsMouseHoveringRect(strip_pos, ImVec2(strip_pos.x + strip_avail.x, strip_pos.y + strip_avail.y)) && ImGui::IsWindowHovered(ImGuiHoveredFlags_ChildWindows));
}

static void draw_movie_overlay_section(ApplicationState* data, float movie_len) {
    auto& m = data->movie;

    auto add = [&](MovieOverlayType type) {
        MovieOverlay o;
        o.type = type;
        o.begin = (double)m.playhead;
        o.end = (double)MAX(movie_len, m.playhead);
        switch (type) {
        case MovieOverlayType::Text:      snprintf(o.text, sizeof(o.text), "Title"); o.anchor = MovieOverlayAnchor::BottomCenter; break;
        case MovieOverlayType::Timestamp: o.anchor = MovieOverlayAnchor::TopRight; break;
        case MovieOverlayType::ScaleBar:  o.anchor = MovieOverlayAnchor::BottomLeft; break;
        case MovieOverlayType::Image:     o.anchor = MovieOverlayAnchor::BottomRight; o.size = 0.1f; break;
        case MovieOverlayType::TimeBar:
            o.anchor = MovieOverlayAnchor::BottomCenter;
            o.size = 0.025f;
            o.background[0] = o.background[1] = o.background[2] = 0.15f;
            o.background[3] = 0.65f;
            break;
        case MovieOverlayType::Timeline:
        case MovieOverlayType::Distribution: {
            movie_overlay_plot_defaults(&o);
            if (type == MovieOverlayType::Distribution) o.num_bins = 128;
            // The subplots of its window that have series
            const bool tl = type == MovieOverlayType::Timeline;
            const PlotSubplot* subs = tl ? data->timeline.subplots : data->distributions.subplots;
            const int n = tl ? data->timeline.num_subplots : data->distributions.num_subplots;
            for (int s = 0; s < n && s < PLOT_MAX_SUBPLOTS; ++s) {
                if (subs[s].count > 0) o.panels.push_back({tl ? MoviePlotView::Timeline : MoviePlotView::Distribution, subs[s].id});
            }
            o.background[0] = 0.0f; o.background[1] = 0.0f; o.background[2] = 0.0f; o.background[3] = 0.5f;
            break;
        }
        case MovieOverlayType::Logo:      o = movie_overlay_default_logo(); o.begin = (double)m.playhead; o.end = (double)MAX(movie_len, m.playhead); break;
        default: break;
        }
        m.overlays.push_back(o);
        m.overlay_selected = (int)m.overlays.size() - 1;
    };

    if (ImGui::Button("Add overlay...")) ImGui::OpenPopup("Add overlay");
    if (ImGui::BeginPopup("Add overlay")) {
        for (int t = 0; t < (int)MovieOverlayType::Count; ++t) {
            if (t == (int)MovieOverlayType::Figure) continue;
            if (ImGui::MenuItem(movie_overlay_type_str[t])) {
                if (t == (int)MovieOverlayType::Image) {
                    char path_buf[2048] = "";
                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Open, STR_LIT("png,jpg,jpeg"))) {
                        add(MovieOverlayType::Image);
                        snprintf(m.overlays.back().path, sizeof(m.overlays.back().path), "%s", path_buf);
                    }
                } else {
                    add((MovieOverlayType)t);
                }
            }
        }
        ImGui::EndPopup();
    }
    m.overlay_selected = CLAMP(m.overlay_selected, 0, MAX((int)m.overlays.size() - 1, 0));
    int remove_idx = -1, duplicate_idx = -1;
    ImGui::SameLine();
    ImGui::BeginDisabled(m.overlays.empty());
    if (ImGui::Button("Duplicate")) duplicate_idx = m.overlay_selected;
    ImGui::SameLine();
    if (ImGui::Button("Remove")) remove_idx = m.overlay_selected;
    ImGui::EndDisabled();
    if (ImGui::BeginListBox("##overlays", ImVec2(-1, ImGui::GetTextLineHeightWithSpacing() * 5.5f))) {
        for (int i = 0; i < (int)m.overlays.size(); ++i) {
            MovieOverlay& o = m.overlays[i];
            ImGui::PushID(i);
            ImGui::Checkbox("##enabled", &o.enabled);
            ImGui::SameLine();
            char label[192];
            snprintf(label, sizeof(label), "%d  %s%s%s", i + 1, movie_overlay_type_str[(int)o.type],
                o.text[0] ? ": " : "", o.text);
            if (ImGui::Selectable(label, i == m.overlay_selected)) m.overlay_selected = i;
            ImGui::PopID();
        }
        ImGui::EndListBox();
    }
    if (m.overlays.empty()) ImGui::TextDisabled("Add an overlay to edit its settings.");
    for (int i = 0; i < (int)m.overlays.size(); ++i) {
        if (i != m.overlay_selected) continue;
        MovieOverlay& o = m.overlays[i];
        ImGui::PushID(i);
        if (!ImGui::BeginTabBar("Overlay inspector")) {
            ImGui::PopID();
            continue;
        }
        if (ImGui::BeginTabItem("Content")) {
            int type = (int)o.type;
            if (ImGui::BeginCombo("Type", movie_overlay_type_str[type])) {
                for (int t = 0; t < (int)MovieOverlayType::Count; ++t) {
                    // The first version of the plots had both kinds in one, which a workspace turns into the two when it is read
                    if (t == (int)MovieOverlayType::Figure) continue;
                    if (ImGui::Selectable(movie_overlay_type_str[t], t == type) && t != type) {
                        o.type = (MovieOverlayType)t;
                        if (o.type == MovieOverlayType::Timeline || o.type == MovieOverlayType::Distribution) {
                            // Only the subplots of its own window, and the place and size of its kind
                            const MoviePlotView own = o.type == MovieOverlayType::Timeline ? MoviePlotView::Timeline : MoviePlotView::Distribution;
                            o.panels.erase(std::remove_if(o.panels.begin(), o.panels.end(), [own](const MoviePlotPanel& panel) { return panel.view != own; }), o.panels.end());
                            movie_overlay_plot_defaults(&o);
                        }
                    }
                }
                ImGui::EndCombo();
            }
            if (o.type == MovieOverlayType::Text) {
                ImGui::InputText("Text", o.text, sizeof(o.text));
            }
            if (o.type == MovieOverlayType::Timeline || o.type == MovieOverlayType::Distribution) {
                const bool kind_tl = o.type == MovieOverlayType::Timeline;
                const MoviePlotView own = kind_tl ? MoviePlotView::Timeline : MoviePlotView::Distribution;
                if (!kind_tl) {
                    bool source_bins = o.num_bins == 0;
                    if (ImGui::Checkbox("Use source bin counts", &source_bins)) o.num_bins = source_bins ? 0 : 128;
                    if (!source_bins) {
                        ImGui::InputInt("Number of bins", &o.num_bins, 1, 16);
                        o.num_bins = CLAMP(o.num_bins, MOVIE_DISTRIBUTION_MIN_BINS, MOVIE_DISTRIBUTION_MAX_BINS);
                    }
                    ImGui::SetItemTooltip("Applies to all subplots in this overlay, without changing the Distributions window.\nScript distributions can only be coarsened to a divisor of their evaluated bin count.");
                }
                // The subplots in the stack, top first. A subplot is found by its id, so moving or renaming the subplots in the
                // windows does not change what is drawn.
                int remove_panel = -1, move_up = -1;
                for (int pi = 0; pi < (int)o.panels.size(); ++pi) {
                    const MoviePlotPanel& panel = o.panels[pi];
                    if (panel.view != own) continue;
                    const bool tl = panel.view == MoviePlotView::Timeline;
                    const PlotSubplot* subs = tl ? data->timeline.subplots : data->distributions.subplots;
                    const int idx = plot_find_subplot(subs, tl ? data->timeline.num_subplots : data->distributions.num_subplots, panel.subplot);
                    char name[48] = "";
                    if (idx >= 0) plot_subplot_label(name, sizeof(name), subs[idx], idx);
                    ImGui::PushID(pi);
                    const bool panel_open = ImGui::TreeNode("Subplot", "%s", idx >= 0 ? name : "Missing subplot");
                    if (panel_open) {
                        if (ImGui::SmallButton("Up") && pi > 0) move_up = pi;
                        ImGui::SameLine();
                        if (ImGui::SmallButton("Remove")) remove_panel = pi;
                        ImGui::SameLine();
                        if (idx < 0) ImGui::TextColored(ImVec4(1.0f, 0.7f, 0.3f, 1.0f), "%s: the subplot is not there any more", tl ? "Timelines" : "Distributions");
                        else if (subs[idx].count == 0) ImGui::TextColored(ImVec4(1.0f, 0.7f, 0.3f, 1.0f), "%s: %s (no series, not drawn)", tl ? "Timelines" : "Distributions", name);
                        else ImGui::Text("%s: %s (%d series)", tl ? "Timelines" : "Distributions", name, subs[idx].count);
                        // When it comes in and goes, inside the overlay's own range
                        float shown_from = (float)o.panels[pi].begin, shown_to = (float)o.panels[pi].end;
                        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 16.0f);
                        if (ImGui::DragFloatRange2("##panel_time", &shown_from, &shown_to, 0.05f, 0.0f, movie_len, "in at %.2f s", shown_to > shown_from ? "out at %.2f s" : "stays to the end")) {
                            o.panels[pi].begin = (double)shown_from;
                            o.panels[pi].end = (double)shown_to;
                        }
                        ImGui::SetItemTooltip("When this subplot comes in and goes, inside the range the overlay is shown in. A property that comes in later starts to be drawn there.\nIts place stays free until then, so the others do not move. Drag the right end back to the left end to let it stay to the end.");
                        ImGui::SameLine();
                        if (ImGui::SmallButton("In at preview time")) o.panels[pi].begin = MAX((double)m.playhead, o.begin);
                        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 16.0f);
                        ImGui::InputTextWithHint("##panel_title", idx >= 0 && subs[idx].name[0] != '\0' ? subs[idx].name : "Title (the name of the subplot)", o.panels[pi].title, sizeof(o.panels[pi].title));
                        ImGui::SetItemTooltip("The title written above this subplot when Titles is ticked. Empty takes the name of the subplot (set in the Subplots menu of its window).");
                        ImGui::TreePop();
                    }
                    ImGui::PopID();
                }
                if (remove_panel >= 0) o.panels.erase(o.panels.begin() + remove_panel);
                if (move_up > 0 && move_up < (int)o.panels.size()) std::swap(o.panels[move_up], o.panels[move_up - 1]);
                if (o.panels.empty()) ImGui::TextColored(ImVec4(1.0f, 0.7f, 0.3f, 1.0f), "No subplot: add one below (drag series into subplots in the %s window).", kind_tl ? "Timelines" : "Distributions");

                if (ImGui::BeginCombo("Add subplot", "Choose...")) {
                    for (int v = (int)own; v == (int)own; ++v) {
                        const bool tl = v == (int)MoviePlotView::Timeline;
                        const PlotSubplot* subs = tl ? data->timeline.subplots : data->distributions.subplots;
                        const int nsub = MAX(tl ? data->timeline.num_subplots : data->distributions.num_subplots, 1);
                        for (int s = 0; s < nsub; ++s) {
                            char name[48], item[96];
                            plot_subplot_label(name, sizeof(name), subs[s], s);
                            snprintf(item, sizeof(item), "%s: %s (%d series)###add%d_%d", tl ? "Timelines" : "Distributions", name, subs[s].count, v, s);
                            if (ImGui::Selectable(item)) {
                                // Added with the movie at a time inside the overlay, it comes in there
                                MoviePlotPanel added;
                                added.view = (MoviePlotView)v;
                                added.subplot = subs[s].id;
                                if ((double)m.playhead > o.begin + 1e-3) added.begin = (double)m.playhead;
                                o.panels.push_back(added);
                            }
                        }
                    }
                    ImGui::EndCombo();
                }
                ImGui::SetItemTooltip("Adds a subplot of one of the windows at the bottom of the stack. At most six series of each are drawn.");

                const char* axes[] = {"Elapsed time", "Trajectory time"};
                int axis = (int)o.plot_axis;
                if (kind_tl && ImGui::Combo("Horizontal axis", &axis, axes, (int)MoviePlotAxis::Count)) o.plot_axis = (MoviePlotAxis)axis;
                if (kind_tl) ImGui::SetItemTooltip("Elapsed time: the trajectory time that the movie has covered, like the time bar, so the curve always grows to the right,\nalso where the trajectory is played backward.\nTrajectory time: the time of the trajectory itself, turned around when the movie plays it backward.");
                ImGui::SliderFloat("Width", &o.width, 0.1f, 1.0f, "%.2f of the frame");
                ImGui::Checkbox("As the movie plays", &o.reveal);
                ImGui::SetItemTooltip("Only the part of the trajectory that the movie has played is drawn, and it grows (a timeline), and only the frames\nthat have been played are counted, so the bars grow (a distribution). A script distribution (not over frames) is drawn as it is.");
                ImGui::SameLine();
                ImGui::Checkbox("Titles", &o.show_titles);
                ImGui::SetItemTooltip("The name of each subplot above it, if it has one (name them in the Subplots menu of the window)");
                ImGui::SameLine();
                ImGui::Checkbox("Value", &o.show_value);
                ImGui::SetItemTooltip("The value at the frame that is shown, in the legend");
                if (kind_tl) {
                    ImGui::SameLine();
                    ImGui::Checkbox("Markers", &o.show_markers);
                    ImGui::SetItemTooltip("The markers of the movie (below), where the movie gets to them");
                }

                if (ImGui::TreeNode("Plot style")) {
                    ImGui::DragFloat("Text (points)", &o.font_points, 0.5f, 0.0f, 200.0f, o.font_points > 0.0f ? "%.0f pt" : "follows the height");
                    ImGui::SetItemTooltip("The size of the text. A point is a pixel of a frame that is 1080 pixels high, scaled with the frame.");
                    o.font_points = MAX(o.font_points, 0.0f);
                    ImGui::DragFloat("Lines (points)", &o.line_points, 0.1f, 0.0f, 20.0f, o.line_points > 0.0f ? "%.1f pt" : "follows the text");
                    o.line_points = MAX(o.line_points, 0.0f);
                    if (ImGui::BeginCombo("Colours", movie_plot_palette_name(o.palette))) {
                        for (int pal = 0; pal < MOVIE_PLOT_PALETTE_COUNT; ++pal) {
                            if (ImGui::Selectable(movie_plot_palette_name(pal), pal == o.palette)) o.palette = pal;
                        }
                        ImGui::EndCombo();
                    }
                    ImGui::SetItemTooltip("The colours of the series: the ones they have in the plots, or a set of its own, in the order of the stack.");
                    if (ImGui::SmallButton("Light text on a dark plate")) {
                        const float c[4] = {1, 1, 1, 1}, b[4] = {0, 0, 0, 0.5f};
                        memcpy(o.color, c, sizeof(c)); memcpy(o.background, b, sizeof(b));
                    }
                    ImGui::SameLine();
                    if (ImGui::SmallButton("Dark text on a light plate")) {
                        const float c[4] = {0.1f, 0.1f, 0.1f, 1}, b[4] = {1, 1, 1, 0.7f};
                        memcpy(o.color, c, sizeof(c)); memcpy(o.background, b, sizeof(b));
                    }
                    ImGui::TreePop();
                }
            }
            if (o.type == MovieOverlayType::PropertyVis) {
                ImGui::InputText("Property", o.text, sizeof(o.text));
                ImGui::SetItemTooltip("The identifier of a script property, as it is named in the script");
                if (ImGui::BeginCombo("##property_pick", "Pick from the script", ImGuiComboFlags_NoPreview)) {
                    const uint32_t kinds = MD_SCRIPT_PROPERTY_FLAG_TEMPORAL | MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION | MD_SCRIPT_PROPERTY_FLAG_VOLUME | MD_SCRIPT_PROPERTY_FLAG_SDF;
                    int listed = 0;
                    series_for_each_script_property(data, SeriesSource_Script, kinds, [&](const SeriesKey& key) {
                        const str_t ident = series_script_ident(key);
                        if (str_empty(ident)) return;
                        char name[128];
                        snprintf(name, sizeof(name), STR_FMT, STR_ARG(ident));
                        if (ImGui::Selectable(name, strcmp(name, o.text) == 0)) snprintf(o.text, sizeof(o.text), "%s", name);
                        listed += 1;
                    });
                    if (listed == 0) ImGui::TextDisabled("The script has no evaluated property");
                    ImGui::EndCombo();
                }
                ImGui::SameLine();
                ImGui::TextDisabled("Pick");
                if (o.text[0] != '\0' && (!data->script.eval_ir || !md_script_ir_property_vis_payload(data->script.eval_ir, str_from_cstr(o.text)))) {
                    ImGui::TextColored(ImVec4(1.0f, 0.7f, 0.3f, 1.0f), "The script has no such property (yet)");
                }
                ImGui::TextDisabled("Shown in the viewport at the preview time and in the recording, while this overlay is shown.");
            }
            if (o.type == MovieOverlayType::TimeBar) {
                ImGui::SliderFloat("Width", &o.width, 0.05f, 1.0f, "%.2f of the frame");
                ImGui::Checkbox("Time that has gone", &o.show_elapsed);
                ImGui::SetItemTooltip("The trajectory time that the movie has covered so far, over the whole (in the unit of the timeline).\nIt counts forward even when the trajectory is played backward.");
                ImGui::SameLine();
                ImGui::Checkbox("Speed", &o.show_speed);
                ImGui::SetItemTooltip("How fast the trajectory is played, as a multiple of the speed of the Animation panel (x1.0 is as there).");
            }
            if (o.type == MovieOverlayType::Image) {
                ImGui::InputText("File", o.path, sizeof(o.path));
                ImGui::SameLine();
                if (ImGui::SmallButton("Browse...")) {
                    char path_buf[2048] = "";
                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Open, STR_LIT("png,jpg,jpeg"))) {
                        snprintf(o.path, sizeof(o.path), "%s", path_buf);
                    }
                }
                ImGui::SameLine();
                if (ImGui::SmallButton("Reload")) movie_image_forget(o.path);
                ImGui::SetItemTooltip("Reads the file again, after it was changed");
                float aspect = 1.0f;
                if (o.path[0] == '\0')                            ImGui::TextDisabled("No file chosen");
                else if (movie_image_texture(o.path, &aspect) == 0) ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), "Could not read this file (a png or jpg image is needed)");
            }
            if (o.type == MovieOverlayType::ScaleBar) {
                ImGui::DragFloat("Length (\xC3\x85)", &o.length, 0.1f, 0.0f, 10000.0f, o.length > 0.0f ? "%.2f" : "automatic");
                ImGui::SetItemTooltip("0 chooses a length that suits the frame: 1, 2 or 5 times a power of ten");
                o.length = MAX(o.length, 0.0f);
            }
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Timing")) {
            float range[2] = {(float)o.begin, (float)o.end};
            if (ImGui::DragFloatRange2("Shown (s)", &range[0], &range[1], 0.05f, 0.0f, movie_len, "from %.2f", "to %.2f")) {
                o.begin = range[0];
                o.end = range[1];
            }
            if (ImGui::SmallButton("Start at preview time")) o.begin = MIN((double)m.playhead, o.end);
            ImGui::SameLine();
            if (ImGui::SmallButton("End at preview time")) o.end = MAX((double)m.playhead, o.begin);

            ImGui::DragFloat("Fade in (s)", &o.fade_in, 0.01f, 0.0f, 10.0f, "%.2f");
            ImGui::DragFloat("Fade out (s)", &o.fade_out, 0.01f, 0.0f, 10.0f, "%.2f");
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Appearance")) {
            int anchor = (int)o.anchor;
            if (ImGui::Combo("Position", &anchor, movie_overlay_anchor_str, (int)MovieOverlayAnchor::Count)) o.anchor = (MovieOverlayAnchor)anchor;
            {
                float lo, hi;
                int unit = (int)o.size_unit;
                const char* unit_names[] = {"% of frame height", "points"};
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 9.0f);
                if (ImGui::Combo("##size_unit", &unit, unit_names, (int)MovieOverlaySizeUnit::Count)) {
                    // The size is kept as it is, written in the other unit
                    o.size = movie_overlay_convert_size(o.size, o.size_unit, (MovieOverlaySizeUnit)unit);
                    o.size_unit = (MovieOverlaySizeUnit)unit;
                }
                ImGui::SetItemTooltip("Percent: the height is a part of the height of the frame.\nPoints: a point is a pixel of a frame that is 1080 pixels high, scaled with the frame, so a size looks the same at any resolution.");
                ImGui::SameLine();
                movie_overlay_size_range(o.size_unit, &lo, &hi);
                if (o.size_unit == MovieOverlaySizeUnit::Points) {
                    ImGui::SliderFloat("Size", &o.size, lo, hi, "%.0f pt", ImGuiSliderFlags_Logarithmic);
                } else {
                    float percent = o.size * 100.0f;
                    if (ImGui::SliderFloat("Size", &percent, lo * 100.0f, hi * 100.0f, "%.1f %%", ImGuiSliderFlags_Logarithmic)) o.size = percent * 0.01f;
                }
                o.size = CLAMP(o.size, lo, hi);
                ImGui::SetItemTooltip("The height of the text (of the logo)");
            }
            ImGui::ColorEdit4("Color", o.color, ImGuiColorEditFlags_AlphaBar | ImGuiColorEditFlags_NoInputs);
            ImGui::ColorEdit4("Background", o.background, ImGuiColorEditFlags_AlphaBar | ImGuiColorEditFlags_NoInputs);
            ImGui::SetItemTooltip("A plate behind it, to read it over a busy picture. None while its opacity is 0.");
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
        ImGui::PopID();
    }
    if (remove_idx >= 0) {
        m.overlays.erase(m.overlays.begin() + remove_idx);
        m.overlay_selected = MIN(remove_idx, MAX((int)m.overlays.size() - 1, 0));
    }
    if (duplicate_idx >= 0 && duplicate_idx < (int)m.overlays.size()) {
        const MovieOverlay copy = m.overlays[duplicate_idx];
        m.overlays.insert(m.overlays.begin() + duplicate_idx + 1, copy);
        m.overlay_selected = duplicate_idx + 1;
    }

    if (ImGui::TreeNode("Timeline markers")) {
        ImGui::TextWrapped("Notes placed on all timeline subplots or on one selected subplot. Each has its color, for its line and label in the overlays and its triangle in the lanes.");
        if (ImGui::Button("Add marker at the preview time")) {
            MovieMarker k;
            k.time = (double)m.playhead;
            snprintf(k.label, sizeof(k.label), "Marker");
            m.markers.push_back(k);
        }
        ImGui::SameLine();
        if (ImGui::Button("Add one at each camera key")) {
            for (size_t i = 0; i < md_array_size(data->movie.keyframes); ++i) {
                MovieMarker k;
                k.time = (double)data->movie.keyframes[i].time;
                snprintf(k.label, sizeof(k.label), "Key %d", (int)i + 1);
                m.markers.push_back(k);
            }
        }
        int remove_marker = -1;
        for (int i = 0; i < (int)m.markers.size(); ++i) {
            MovieMarker& k = m.markers[i];
            ImGui::PushID(i);
            float shown[4];
            movie_marker_color(m.markers.data(), m.markers.size(), (size_t)i, shown);
            if (ImGui::ColorEdit4("##color", shown, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoAlpha)) {
                for (int c = 0; c < 3; ++c) k.color[c] = shown[c];
                k.color[3] = 1.0f;
            }
            if (ImGui::BeginPopupContextItem("marker_color_menu")) {
                if (ImGui::MenuItem("Automatic color", nullptr, false, k.color[3] > 0.0f)) k.color[3] = 0.0f;
                ImGui::EndPopup();
            }
            ImGui::SetItemTooltip(k.color[3] > 0.0f ? "Its own color (right click: back to automatic)." : "Automatic: a color of its own by its place in time. Click to pick one.");
            ImGui::SameLine();
            float t = (float)k.time;
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 7.0f);
            if (ImGui::DragFloat("##time", &t, 0.05f, 0.0f, movie_len, "%.2f s")) k.time = (double)CLAMP(t, 0.0f, movie_len);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14.0f);
            ImGui::InputText("##label", k.label, sizeof(k.label));
            ImGui::SameLine();
            if (ImGui::SmallButton("Remove")) remove_marker = i;
            const int subplot_index = plot_find_subplot(data->timeline.subplots, data->timeline.num_subplots, k.subplot);
            char target[64] = "All timeline subplots";
            if (k.subplot != 0) {
                if (subplot_index >= 0) plot_subplot_label(target, sizeof(target), data->timeline.subplots[subplot_index], subplot_index);
                else snprintf(target, sizeof(target), "Missing subplot (%u)", k.subplot);
            }
            ImGui::SetNextItemWidth(-1);
            if (ImGui::BeginCombo("##marker_subplot", target)) {
                if (ImGui::Selectable("All timeline subplots", k.subplot == 0)) k.subplot = 0;
                for (int s = 0; s < data->timeline.num_subplots; ++s) {
                    const PlotSubplot& subplot = data->timeline.subplots[s];
                    char name[48];
                    plot_subplot_label(name, sizeof(name), subplot, s);
                    ImGui::PushID(s);
                    if (ImGui::Selectable(name, k.subplot == subplot.id)) k.subplot = subplot.id;
                    ImGui::PopID();
                }
                ImGui::EndCombo();
            }
            ImGui::SetItemTooltip("Draw this marker on every timeline subplot, or just the selected one. Saved by subplot identity.");
            if (k.subplot != 0 && subplot_index < 0) ImGui::TextColored(ImVec4(1.0f, 0.7f, 0.3f, 1.0f), "The target subplot is not available.");
            ImGui::PopID();
        }
        if (remove_marker >= 0) m.markers.erase(m.markers.begin() + remove_marker);
        ImGui::TreePop();
    }
}

// A field for a selection in the language of the filters; Enter applies a valid one. 'current' is shown when the popup opens.
static bool movie_expr_field(const ApplicationState* data, const char* label, const std::string& current, std::string* applied) {
    static char buf[256];
    static std::string checked, error;
    static bool valid = true;
    if (ImGui::IsWindowAppearing()) {
        snprintf(buf, sizeof(buf), "%s", current.c_str());
        checked = "\x01";
    }
    if (checked != buf) {
        checked = buf;
        std::vector<uint32_t> atoms;
        error.clear();
        valid = checked.empty() || movie_expr_atoms(data, checked, &atoms, nullptr, &error);
    }
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 16.0f);
    const bool entered = ImGui::InputQuery(label, buf, sizeof(buf), valid, error.c_str(), ImGuiInputTextFlags_EnterReturnsTrue);
    ImGui::SetItemTooltip("A selection: protein, resname(\"AIN\"), resid(12) ... Tracked through the trajectory. Enter applies it.");
    if (entered && valid && buf[0] != '\0') {
        *applied = buf;
        return true;
    }
    return false;
}

// Keys the blur of the depth of field at a time, as it is now
static void movie_key_blur_at(ApplicationState* state, double time) {
    auto& m = state->movie;
    const int num_params = (int)(sizeof(movie_param_table) / sizeof(movie_param_table[0]));
    for (int p = 0; p < num_params; ++p) {
        if (movie_param_table[p].id == MOVIE_PARAM_DOF_APERTURE) m.param_selected = p;
    }
    ParamKey key;
    key.param = MOVIE_PARAM_DOF_APERTURE;
    key.time = time;
    key.value[0] = state->visuals.dof.aperture;
    for (ParamKey& k : m.param_keys) {
        if (k.param == key.param && fabs(k.time - time) < 1.0e-3) return;
    }
    m.param_keys.push_back(key);
    movie_param_sort(state);
}

// What a key looks at: atoms that it tracks through the trajectory, or a fixed point
static void draw_movie_key_look_popup(ApplicationState* data, int i) {
    auto& m = data->movie;
    if (i < 0 || i >= (int)md_array_size(m.keyframes)) return;
    CameraKeyframe& key = m.keyframes[i];
    ImGui::TextDisabled("Keyframe %d looks at", i + 1);
    if (ImGui::RadioButton("Fixed point", key.look_set == 0)) {
        vec3_t point = camera_get_look_at(key.transform);
        if (key.look_set != 0 && !movie_set_center(data, key.look_set, &point))
            VIAMD_LOG_ERROR("Look-at target unavailable: freezing its saved point");
        movie_key_look_at_point(data, i, point);
    }
    ImGui::SameLine();
    if (ImGui::RadioButton("Track selection", key.look_set != 0)) {
        const auto atoms = movie_selection_atoms(data);
        if (!atoms.empty()) movie_key_look_at_atoms(data, i, atoms);
        else {
            m.look_pick_key = i;
            m.focus_pick_key = -1;
            ImGui::CloseCurrentPopup();
        }
    }
    ImGui::SetItemTooltip("Track selected atoms; with no selection, pick an atom in the viewport.\nYou can also type a selection expression below.");
    if (key.look_set == 0) {
        vec3_t point = camera_get_look_at(key.transform);
        if (ImGui::DragFloat3("Fixed point##look", point.elem, 0.1f)) movie_key_look_at_point(data, i, point);
        if (ImGui::Button("Use view aim")) movie_key_look_at_point(data, i, camera_get_look_at(data->view.target));
    } else ImGui::TextDisabled("Tracking: %s", movie_set_label(data, key.look_set).c_str());
    ImGui::Separator();
    ImGui::TextDisabled("Track selection");
    {
        std::string expr;
        const MovieAtomSet* set = movie_atoms_find(m.atom_sets, key.look_set);
        if (movie_expr_field(data, "##look_expr", set ? set->expr : std::string(), &expr) && movie_key_look_at_expr(data, i, expr)) ImGui::CloseCurrentPopup();
    }
    const std::vector<uint32_t> selected = movie_selection_atoms(data);
    char label[64];
    snprintf(label, sizeof(label), "The selection (%zu atoms)", selected.size());
    if (ImGui::MenuItem(label, nullptr, false, !selected.empty())) movie_key_look_at_atoms(data, i, selected);
    ImGui::SetItemTooltip("The centre of the selected atoms, tracked through the trajectory. The eye stays where it is.");
    if (ImGui::MenuItem("Pick an atom in the viewport")) {
        m.look_pick_key = i;
        m.focus_pick_key = -1;
    }
    ImGui::SetItemTooltip("Click an atom: it is tracked through the trajectory. Esc cancels.");
}

// Focus defaults to this key's look-at; an explicit target decouples it.
static void draw_movie_key_focus_popup(ApplicationState* data, int i) {
    auto& m = data->movie;
    if (i < 0 || i >= (int)md_array_size(m.keyframes)) return;
    CameraKeyframe& key = m.keyframes[i];
    ImGui::TextDisabled("Keyframe %d: what is sharp", i + 1);
    if (ImGui::RadioButton("Look-at", !key.focus_on || key.focus_target == FocusTarget::LookAt)) {
        key.focus_on = false;
        key.focus_target = FocusTarget::LookAt;
        key.focus_set = 0;
    }
    ImGui::SetItemTooltip("Focus follows the movie camera's aim at this key, even if an earlier key focused elsewhere.");
    ImGui::SameLine();
    if (ImGui::RadioButton("Track selection", key.focus_on && key.focus_target == FocusTarget::Selection)) {
        const auto atoms = movie_selection_atoms(data);
        if (!atoms.empty()) movie_key_focus_on_atoms(data, i, atoms);
        else {
            m.focus_pick_key = i;
            m.look_pick_key = -1;
            ImGui::CloseCurrentPopup();
        }
    }
    ImGui::SetItemTooltip("Track selected atoms; with no selection, pick an atom in the viewport.\nYou can also type a selection expression below.");
    if (key.focus_on && (key.focus_target == FocusTarget::Point || key.focus_target == FocusTarget::Distance))
        ImGui::TextDisabled("Saved legacy focus: %s (preserved until you choose another target)", movie_focus_label(data, key).c_str());
    ImGui::Separator();
    ImGui::TextDisabled("Track selection");
    {
        std::string expr;
        const MovieAtomSet* set = key.focus_target == FocusTarget::Selection ? movie_atoms_find(m.atom_sets, key.focus_set) : nullptr;
        if (movie_expr_field(data, "##focus_expr", set ? set->expr : std::string(), &expr)) movie_key_focus_on_expr(data, i, expr);
    }
    const std::vector<uint32_t> selected = movie_selection_atoms(data);
    char label[64];
    snprintf(label, sizeof(label), "The selection (%zu atoms)", selected.size());
    if (ImGui::Selectable(label, key.focus_target == FocusTarget::Selection, ImGuiSelectableFlags_NoAutoClosePopups | (selected.empty() ? ImGuiSelectableFlags_Disabled : 0))) movie_key_focus_on_atoms(data, i, selected);
    ImGui::SetItemTooltip("The centre of the selected atoms, tracked through the trajectory.");
    if (ImGui::Selectable("Pick an atom in the viewport")) {
        m.focus_pick_key = i;
        m.look_pick_key = -1;
        ImGui::CloseCurrentPopup();
    }
    ImGui::SetItemTooltip("Click an atom: it is tracked through the trajectory. Esc cancels.");
    ImGui::Separator();
    ImGui::BeginDisabled(!key.focus_on);
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8.0f);
    ImGui::DragFloat("Transition##focus_transition", &key.focus_transition, 0.05f, 0.0f, 60.0f, "%.1f s");
    ImGui::SetItemTooltip("How long the focus takes to change to this, from the key's time.");
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8.0f);
    int ease = (int)key.focus_ease;
    if (ImGui::Combo("Ease##focus_ease", &ease, key_ease_str, (int)KeyEase::Count)) key.focus_ease = (KeyEase)ease;
    ImGui::EndDisabled();
    ImGui::TextDisabled("How blurred the rest is: Looks > Depth of field blur");
    if (ImGui::SmallButton("Key the blur here")) movie_key_blur_at(data, key.time);
    ImGui::SetItemTooltip("Adds a key of Depth of field blur at this key's time, with the blur as it is now. Edit and smooth it in the Looks tab\nor the look parameter lane. The blur is on when it is above 0 and a key sets the focus.");
}

static void draw_movie_keyframe_table(ApplicationState* data, float movie_len, bool locked) {
    auto& m = data->movie;
    if (md_array_size(m.keyframes) == 0) {
        ImGui::TextDisabled("Move the view, then 'Add Keyframe' (K) at a time on the timeline.");
        return;
    }

    bool resort = false;
    int  remove_idx = -1;
    int  dup_idx = -1;
    int  move_from = -1;  // Dragged to the place of move_to
    int  move_to = -1;

    const double last_frame = (double)(run_num_frames(data) > 0 ? run_num_frames(data) - 1 : 0);
    double prev_frame = -1.0;

    // Columns can be resized, and the table scrolls sideways when they do not fit
    const float fs = ImGui::GetFontSize();
    if (ImGui::BeginTable("##keyframes_focus_v2", 15, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_Resizable | ImGuiTableFlags_ScrollX)) {
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, fs * 1.8f);
        ImGui::TableSetupColumn("Name", ImGuiTableColumnFlags_WidthFixed, fs * 7.0f);
        ImGui::TableSetupColumn("Time (s)", ImGuiTableColumnFlags_WidthFixed, fs * 5.0f);
        ImGui::TableSetupColumn("Frame", ImGuiTableColumnFlags_WidthFixed, fs * 5.0f);
        ImGui::TableSetupColumn("FOV (deg)", ImGuiTableColumnFlags_WidthFixed, fs * 5.0f);
        ImGui::TableSetupColumn("Go to", ImGuiTableColumnFlags_WidthFixed, fs * 2.5f);
        ImGui::TableSetupColumn("Update", ImGuiTableColumnFlags_WidthFixed, fs * 3.0f);
        ImGui::TableSetupColumn("Look at", ImGuiTableColumnFlags_WidthFixed, fs * 6.5f);
        ImGui::TableSetupColumn("Focus on", ImGuiTableColumnFlags_WidthFixed, fs * 6.5f);
        ImGui::TableSetupColumn("Ease", ImGuiTableColumnFlags_WidthFixed, fs * 6.0f);
        ImGui::TableSetupColumn("Spin", ImGuiTableColumnFlags_WidthFixed, fs * 7.0f);
        ImGui::TableSetupColumn("Roll (deg)", ImGuiTableColumnFlags_WidthFixed, fs * 6.0f);
        ImGui::TableSetupColumn("Copy", ImGuiTableColumnFlags_WidthFixed, fs * 2.5f);
        ImGui::TableSetupColumn("Dup", ImGuiTableColumnFlags_WidthFixed, fs * 2.5f);
        ImGui::TableSetupColumn("Remove", ImGuiTableColumnFlags_WidthFixed, fs * 3.0f);
        ImGui::TableHeadersRow();

        for (int i = 0; i < (int)md_array_size(m.keyframes); ++i) {
            CameraKeyframe& key = m.keyframes[i];
            ImGui::PushID(i);
            ImGui::TableNextRow();

            ImGui::TableNextColumn();
            char num[16];
            snprintf(num, sizeof(num), "%d", i + 1);
            ImGui::Selectable(num, false, ImGuiSelectableFlags_None);
            if (ImGui::BeginDragDropSource(ImGuiDragDropFlags_None)) {
                ImGui::SetDragDropPayload("MOVIE_KEYFRAME", &i, sizeof(int));
                ImGui::Text("Keyframe %d", i + 1);
                ImGui::EndDragDropSource();
            }
            if (ImGui::BeginDragDropTarget()) {
                if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload("MOVIE_KEYFRAME")) {
                    move_from = *(const int*)payload->Data;
                    move_to = i;
                }
                ImGui::EndDragDropTarget();
            }
            ImGui::SetItemTooltip("Drag to another row to move this keyframe, with its pose, easing and frame, to that place.\nThe times stay where they are in the list.");

            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::InputTextWithHint("##name", "name", key.name, sizeof(key.name));
            ImGui::SetItemTooltip("What this keyframe is called, shown in the camera lane. Saved with the workspace.");

            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            double t = key.time;
            if (ImGui::InputDouble("##time", &t, 0.0, 0.0, "%.2f")) {
                key.time = CLAMP(t, 0.0, (double)movie_len);
            }
            resort |= ImGui::IsItemDeactivatedAfterEdit();

            ImGui::TableNextColumn();
            bool use_frame = key.use_frame;
            if (ImGui::Checkbox("##useframe", &use_frame)) {
                // Taken from the movie as it is, so that switching it on does not change anything
                if (use_frame) key.frame = movie_trajectory_frame(data, key.time);
                key.use_frame = use_frame;
            }
            ImGui::SetItemTooltip("Key the trajectory frame shown at this time. With two or more, the trajectory plays\nfrom one to the next, so the speed can change between them.");
            if (key.use_frame) {
                const bool backward = i > 0 && prev_frame >= 0.0 && key.frame < prev_frame - 0.5;
                const double zero = 0.0;
                ImGui::SameLine();
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::DragScalar("##frame", ImGuiDataType_Double, &key.frame, 0.5f, &zero, &last_frame, "%.0f");
                if (backward) {
                    ImGui::SetItemTooltip("Behind the previous frame: the trajectory plays backward to here.");
                }
                prev_frame = key.frame;
            }

            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            float fov_deg = key.fov_y * MOVIE_RAD_TO_DEG;
            if (ImGui::DragFloat("##fov", &fov_deg, 0.1f, 1.0f, 170.0f, "%.1f")) {
                key.fov_y = fov_deg * MOVIE_DEG_TO_RAD;
            }

            ImGui::TableSetColumnIndex(7);
            {
                const bool picking = m.look_pick_key == i;
                if (picking) ImGui::PushStyleColor(ImGuiCol_Button, ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));
                const std::string text = std::string(ICON_FA_CROSSHAIRS " ") + (key.look_set != 0 ? "Track: " + movie_set_label(data, key.look_set) : "Fixed");
                if (ImGui::Button((text + "##look").c_str(), ImVec2(-FLT_MIN, 0.0f))) ImGui::OpenPopup("##look_popup");
                if (picking) ImGui::PopStyleColor();
                ImGui::SetItemTooltip(picking ? "Click an atom in the viewport. Esc cancels." :
                    "Look at: what this keyframe looks at, a set of atoms tracked through the trajectory or a fixed point.\nThe eye keeps its place.");
                if (ImGui::BeginPopup("##look_popup")) {
                    draw_movie_key_look_popup(data, i);
                    ImGui::EndPopup();
                }
            }

            ImGui::TableNextColumn();
            {
                const bool picking = m.focus_pick_key == i;
                if (picking) ImGui::PushStyleColor(ImGuiCol_Button, ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));
                const std::string text = std::string(ICON_FA_BULLSEYE " ") +
                    (key.focus_on && key.focus_target == FocusTarget::Selection ? "Track: " + movie_focus_label(data, key) :
                        key.focus_on ? movie_focus_label(data, key) : "Look-at");
                if (ImGui::Button((text + "##focus").c_str(), ImVec2(-FLT_MIN, 0.0f))) ImGui::OpenPopup("##focus_popup");
                if (picking) ImGui::PopStyleColor();
                ImGui::SetItemTooltip(picking ? "Click an atom in the viewport. Esc cancels." :
                    "Focus on: defaults to what the movie camera looks at at this key.\nChoose a separate target to track an object without changing aim.");
                if (ImGui::BeginPopup("##focus_popup")) {
                    draw_movie_key_focus_popup(data, i);
                    ImGui::EndPopup();
                }
            }

            ImGui::TableSetColumnIndex(5);
            if (ImGui::SmallButton(ICON_FA_ARROW_RIGHT_TO_BRACKET "##goto")) movie_goto_keyframe(data, (size_t)i);
            ImGui::SetItemTooltip("Go to: move the view to this keyframe");
            ImGui::TableNextColumn();
            if (ImGui::SmallButton(ICON_FA_CAMERA "##update")) {
                // The camera moves to the current view's eye, still looking at the point the keyframe looks at
                const vec3_t look = camera_get_look_at(key.transform);
                ViewTransform t = data->view.target;
                if (camera_aim_at(&t, look)) {
                    key.transform = t;
                    if (const vec3_t* up = movie_upright(data)) key.roll = camera_roll(t, *up);
                }
                key.fov_y = data->view.camera.fov_y;
                if (key.use_frame) key.frame = data->animation.frame;
            }
            ImGui::SetItemTooltip("Update position: move this keyframe's camera to the current view's position. What it looks at stays.");
            ImGui::TableSetColumnIndex(9);
            if (i == 0) {
                ImGui::TextDisabled("-");
                ImGui::SetItemTooltip("How the movie gets to a keyframe is set on that keyframe, so the first has none.");
            } else {
                ImGui::SetNextItemWidth(-FLT_MIN);
                int ease = (int)key.ease;
                if (ImGui::Combo("##ease", &ease, key_ease_str, (int)KeyEase::Count)) key.ease = (KeyEase)ease;
                ImGui::SetItemTooltip("How the camera, and the trajectory frame if it is keyed, move in the stretch leading to this keyframe.\n"
                    "Smooth: through the keys without stopping. Ease in/out: starts and ends slowly.\n"
                    "Linear: constant speed. Hold: stays as it was until the key, then jumps.");
            }

            ImGui::TableNextColumn();
            if (i == 0) {
                ImGui::TextDisabled("-");
                ImGui::SetItemTooltip("A spin is made in the stretch leading to a keyframe, so the first has none.");
            } else {
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 4.0f);
                ImGui::DragInt("##spin", &key.spin_turns, 0.05f, -16, 16, "%d x");
                ImGui::SetItemTooltip("Extra whole turns of the camera around what it looks at, made in the stretch leading to this\nkeyframe. Positive is counter-clockwise seen from the tip of the axis.");
                if (key.spin_turns != 0) {
                    ImGui::SameLine();
                    if (ImGui::SmallButton("...")) ImGui::OpenPopup("##spin_options");
                    if (ImGui::BeginPopup("##spin_options")) {
                        int axis = (int)key.spin_axis;
                        if (ImGui::Combo("Around", &axis, spin_axis_str, (int)SpinAxis::Count)) key.spin_axis = (SpinAxis)axis;
                        ImGui::Checkbox("Constant speed", &key.spin_constant_speed);
                        ImGui::SetItemTooltip("Otherwise it starts and ends slowly.");
                        ImGui::EndPopup();
                    }
                }
            }

            ImGui::TableNextColumn();
            if (m.keep_upright) {
                ImGui::SetNextItemWidth(fs * 3.5f);
                float roll_deg = key.roll * MOVIE_RAD_TO_DEG;
                if (ImGui::DragFloat("##roll", &roll_deg, 0.5f, -180.0f, 180.0f, "%.0f")) key.roll = CLAMP(roll_deg, -180.0f, 180.0f) * MOVIE_DEG_TO_RAD;
                ImGui::SetItemTooltip("The tilt of the camera at this keyframe, positive leaning it to the left.\nThe movie keeps the camera upright and goes smoothly from one keyframe's roll to the next.");
                ImGui::SameLine();
                if (ImGui::SmallButton("0")) key.roll = 0.0f;
                ImGui::SetItemTooltip("Level: no tilt");
            } else {
                const float tilt = camera_roll(key.transform, *movie_up_vector(data)) * MOVIE_RAD_TO_DEG;
                ImGui::AlignTextToFramePadding();
                ImGui::TextDisabled("%.0f", tilt);
                ImGui::SetItemTooltip("The tilt of this keyframe's camera about %s up. 'Keep upright' (in Camera) makes it a setting.", movie_up_axis_str[CLAMP(m.up_axis, 0, 5)]);
                ImGui::SameLine();
                if (ImGui::SmallButton("Level")) camera_level(&key.transform, *movie_up_vector(data));
                ImGui::SetItemTooltip("Turn this keyframe's camera about where it looks so that it is level, %s being up.", movie_up_axis_str[CLAMP(m.up_axis, 0, 5)]);
            }
            ImGui::TableNextColumn();
            if (ImGui::SmallButton(ICON_FA_COPY "##copy")) {
                m.key_clipboard = key;
                m.has_key_clipboard = true;
            }
            ImGui::SetItemTooltip("Copy this keyframe for Paste Keyframe");
            ImGui::TableNextColumn();
            if (ImGui::SmallButton(ICON_FA_CLONE "##duplicate")) dup_idx = i;
            ImGui::SetItemTooltip("Duplicate one second later");
            ImGui::TableNextColumn();
            if (ImGui::SmallButton(ICON_FA_TRASH_CAN "##remove")) remove_idx = i;
            ImGui::SetItemTooltip("Remove this keyframe");
            ImGui::PopID();
        }
        ImGui::EndTable();
    }

    if (md_array_size(m.keyframes) > 0) {
        // The camera in 3D at the preview time, which is what the distance alone does not say
        ViewTransform vt;
        float fov_y;
        movie_camera_evaluate(data, (double)m.playhead, m.keyframes, md_array_size(m.keyframes), &vt, &fov_y);
        const vec3_t look = camera_get_look_at(vt);
        ImGui::TextDisabled("At %.2f s: eye (%.2f, %.2f, %.2f)  looks at (%.2f, %.2f, %.2f)  distance %.2f",
            m.playhead, vt.position.x, vt.position.y, vt.position.z, look.x, look.y, look.z, vt.distance);
        if (data->visuals.dof.enabled || movie_focus_keys_exist(m.keyframes, md_array_size(m.keyframes))) {
            ImGui::TextDisabled("Focus at %.2f from the eye", dof_focus_depth(data, vt));
        }
    }

    if (locked) return;

    if (move_from >= 0 && move_from != move_to) {
        const int n = (int)md_array_size(m.keyframes);
        if (move_from < n && 0 <= move_to && move_to < n) {
            // The times stay with their places in the list: it is the keys that are moved along them
            std::vector<double> times(n);
            for (int k = 0; k < n; ++k) times[k] = m.keyframes[k].time;
            const CameraKeyframe moved = m.keyframes[move_from];
            if (move_from < move_to) memmove(m.keyframes + move_from, m.keyframes + move_from + 1, (size_t)(move_to - move_from) * sizeof(CameraKeyframe));
            else                     memmove(m.keyframes + move_to + 1, m.keyframes + move_to, (size_t)(move_from - move_to) * sizeof(CameraKeyframe));
            m.keyframes[move_to] = moved;
            for (int k = 0; k < n; ++k) m.keyframes[k].time = times[k];
        }
    }
    if (dup_idx >= 0) {
        CameraKeyframe copy = m.keyframes[dup_idx];
        copy.time = MIN(copy.time + 1.0, (double)movie_len);
        copy.spin_turns = 0;
        md_array_push(m.keyframes, copy, data->allocator.persistent);
        resort = true;
    }
    if (remove_idx >= 0) {
        // Keeps the order, unlike swap-and-pop
        CameraKeyframe* keys = m.keyframes;
        const size_t n = md_array_size(keys);
        memmove(keys + remove_idx, keys + remove_idx + 1, (n - remove_idx - 1) * sizeof(CameraKeyframe));
        md_array_pop(keys);
    }
    if (resort) {
        movie_sort_keyframes(data);
    }
}

// Look parameters keyed over the movie: you set one up in the Settings as it should look at some time, and key it.
static void draw_movie_param_section(ApplicationState* data, float movie_len) {
    auto& m = data->movie;
    const int num_params = (int)(sizeof(movie_param_table) / sizeof(movie_param_table[0]));

    if (ImGui::Checkbox("Animate parameters", &m.animate_params)) {
        // Turning it off lets go of them
        movie_params_apply(data, (double)m.playhead);
        movie_reps_apply(data, (double)m.playhead);
    }
    ImGui::SetItemTooltip("The keyed parameters, and the keyed properties of representations, follow their keys when the movie is\nscrubbed, previewed or recorded. Off, they stay as they are.");

    m.param_selected = CLAMP(m.param_selected, 0, num_params - 1);
    const MovieParamDesc& sel = movie_param_table[m.param_selected];
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14.0f);
    if (ImGui::BeginCombo("##param", sel.label)) {
        for (int i = 0; i < num_params; ++i) {
            if (ImGui::Selectable(movie_param_table[i].label, i == m.param_selected)) m.param_selected = i;
        }
        ImGui::EndCombo();
    }
    ImGui::SameLine();
    if (ImGui::Button("Key Now")) {
        movie_key_param(data, sel.id);
    }
    ImGui::SetItemTooltip("Keys the value the parameter has now at the preview time.\nSet it up in the Settings first, then key it.");
    if (sel.tip) ImGui::TextDisabled("%s", sel.tip);
    ImGui::TextDisabled("The keys can be dragged in the timeline, the left panel of this window.");

    if (m.param_keys.empty()) return;

    bool resort = false;
    int  remove_idx = -1;
    int  prev_param = -1;

    if (ImGui::BeginTable("##param_keys", 5, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("Parameter");
        ImGui::TableSetupColumn("Time (s)");
        ImGui::TableSetupColumn("Value");
        ImGui::TableSetupColumn("Ease");
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 4.5f);
        ImGui::TableHeadersRow();

        for (int i = 0; i < (int)m.param_keys.size(); ++i) {
            ParamKey& key = m.param_keys[i];
            const MovieParamDesc* d = movie_param_desc(key.param);
            const bool first_of_param = key.param != prev_param;
            prev_param = key.param;
            ImGui::PushID(i);
            ImGui::TableNextRow();

            ImGui::TableNextColumn();
            ImGui::TextUnformatted(d ? d->label : "(unknown)");

            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            double t = key.time;
            if (ImGui::InputDouble("##time", &t, 0.0, 0.0, "%.2f")) {
                key.time = CLAMP(t, 0.0, (double)movie_len);
            }
            resort |= ImGui::IsItemDeactivatedAfterEdit();

            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (d && d->color) {
                ImGui::ColorEdit3("##value", key.value, ImGuiColorEditFlags_NoInputs);
            } else if (d) {
                ImGui::DragFloat("##value", &key.value[0], (d->hi - d->lo) * 0.005f, d->lo, d->hi, "%.3g");
            }

            ImGui::TableNextColumn();
            if (first_of_param) {
                ImGui::TextDisabled("-");
            } else {
                ImGui::SetNextItemWidth(-FLT_MIN);
                int ease = (int)key.ease;
                if (ImGui::Combo("##ease", &ease, key_ease_str, (int)KeyEase::Count)) key.ease = (KeyEase)ease;
                ImGui::SetItemTooltip("How the value moves in the stretch leading to this key");
            }

            ImGui::TableNextColumn();
            if (ImGui::SmallButton("Remove")) remove_idx = i;
            ImGui::PopID();
        }
        ImGui::EndTable();
    }

    if (remove_idx >= 0) {
        m.param_keys.erase(m.param_keys.begin() + remove_idx);
        // The parameter lets go of its values if this was its last key
        movie_params_apply(data, (double)m.playhead);
    }
    if (resort) {
        movie_param_sort(data);
    }
}

static void movie_rep_sort(ApplicationState* state) {
    std::stable_sort(state->movie.rep_keys.begin(), state->movie.rep_keys.end(), [](const RepKey& a, const RepKey& b) {
        if (a.rep != b.rep) return a.rep < b.rep;
        return a.prop != b.prop ? a.prop < b.prop : a.time < b.time;
    });
}

// Keys what the property is now, at the preview time
static void movie_key_rep(ApplicationState* state, const Representation& rep, int prop) {
    auto& m = state->movie;
    RepKey key;
    key.rep = rep.id;
    key.prop = prop;
    key.time = movie_snap_time(state, (double)m.playhead);
    movie_rep_prop_get(rep, prop, key.value);
    if (prop == (int)RepProp::Visible) {
        key.value[0] = key.value[0] >= 0.5f ? 1.0f : 0.0f;
        key.ease = KeyEase::Hold;
    }

    for (RepKey& k : m.rep_keys) {
        if (k.rep == key.rep && k.prop == key.prop && fabs(k.time - key.time) < 1.0e-3) {
            key.ease = k.ease;
            k = key;
            return;
        }
    }
    m.rep_keys.push_back(key);
    movie_rep_sort(state);
}

// Properties of representations keyed over the movie: show or hide one, change its scale or tint. You set it up
// in the Representations window as it should look at some time, and key it.
static void draw_movie_rep_section(ApplicationState* data, float movie_len) {
    auto& m = data->movie;
    const int num_reps = (int)md_array_size(data->representation.reps);
    if (num_reps == 0 && m.rep_keys.empty()) {
        ImGui::TextDisabled("Create a representation in the Representations window first.");
        return;
    }

    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6.0f);
    ImGui::DragFloat("Transition (s)", &m.rep_transition, 0.05f, 0.0f, 60.0f, "%.2f");
    m.rep_transition = CLAMP(m.rep_transition, 0.0f, 60.0f);
    ImGui::SetItemTooltip("How long a representation takes to grow in or shrink away at a Visible key. It starts at the key.\n0 shows or hides it at once. Keys closer together than this make it turn around on the way.");

    if (num_reps > 0) {
        m.rep_selected = CLAMP(m.rep_selected, 0, num_reps - 1);
        const Representation& rep = data->representation.reps[m.rep_selected];
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 12.0f);
        if (ImGui::BeginCombo("##rep", rep.name)) {
            for (int i = 0; i < num_reps; ++i) {
                ImGui::PushID(i);
                if (ImGui::Selectable(data->representation.reps[i].name, i == m.rep_selected)) m.rep_selected = i;
                ImGui::PopID();
            }
            ImGui::EndCombo();
        }

        float lo, hi;
        m.rep_prop_selected = CLAMP(m.rep_prop_selected, 0, (int)RepProp::Count - 1);
        if (!movie_rep_prop_label(rep, m.rep_prop_selected, &lo, &hi)) m.rep_prop_selected = (int)RepProp::Visible;
        ImGui::SameLine();
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10.0f);
        if (ImGui::BeginCombo("##rep_prop", movie_rep_prop_label(rep, m.rep_prop_selected, &lo, &hi))) {
            for (int p = 0; p < (int)RepProp::Count; ++p) {
                if (const char* label = movie_rep_prop_label(rep, p, &lo, &hi)) {
                    if (ImGui::Selectable(label, p == m.rep_prop_selected)) m.rep_prop_selected = p;
                }
            }
            ImGui::EndCombo();
        }
        ImGui::SameLine();
        if (ImGui::Button("Key Now##rep")) {
            movie_key_rep(data, rep, m.rep_prop_selected);
        }
        ImGui::SetItemTooltip("Keys the value the property has now at the preview time. Set it up in the Representations window first\n(the eye shows or hides it). Visible changes at its keys, the others move smoothly between them.");
    }

    if (m.rep_keys.empty()) return;

    bool resort = false;
    int  remove_idx = -1;
    if (ImGui::BeginTable("##rep_keys", 6, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("Representation");
        ImGui::TableSetupColumn("Property");
        ImGui::TableSetupColumn("Time (s)");
        ImGui::TableSetupColumn("Value");
        ImGui::TableSetupColumn("Ease");
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 4.5f);
        ImGui::TableHeadersRow();

        uint32_t prev_rep = 0;
        int prev_prop = -1;
        for (int i = 0; i < (int)m.rep_keys.size(); ++i) {
            RepKey& key = m.rep_keys[i];
            const Representation* rep = movie_find_rep(data, key.rep);
            const bool first = key.rep != prev_rep || key.prop != prev_prop;
            prev_rep = key.rep;
            prev_prop = key.prop;
            float lo = 0.0f, hi = 1.0f;
            const char* label = rep ? movie_rep_prop_label(*rep, key.prop, &lo, &hi) : nullptr;
            ImGui::PushID(i);
            ImGui::TableNextRow();

            ImGui::TableNextColumn();
            ImGui::TextUnformatted(rep ? rep->name : "(removed)");

            ImGui::TableNextColumn();
            ImGui::TextUnformatted(label ? label : "(not for its type)");

            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            double t = key.time;
            if (ImGui::InputDouble("##time", &t, 0.0, 0.0, "%.2f")) {
                key.time = CLAMP(t, 0.0, (double)movie_len);
            }
            resort |= ImGui::IsItemDeactivatedAfterEdit();

            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (key.prop == (int)RepProp::Visible) {
                bool visible = key.value[0] >= 0.5f;
                if (ImGui::Checkbox("##visible", &visible)) key.value[0] = visible ? 1.0f : 0.0f;
            } else if (rep_prop_comps(key.prop) == 3) {
                ImGui::ColorEdit3("##value", key.value, ImGuiColorEditFlags_NoInputs);
            } else {
                ImGui::DragFloat("##value", &key.value[0], (hi - lo) * 0.005f, lo, hi, "%.3g");
            }

            ImGui::TableNextColumn();
            if (first || key.prop == (int)RepProp::Visible) {
                ImGui::TextDisabled(key.prop == (int)RepProp::Visible && !first ? "at the key" : "-");
            } else {
                ImGui::SetNextItemWidth(-FLT_MIN);
                int ease = (int)key.ease;
                if (ImGui::Combo("##ease", &ease, key_ease_str, (int)KeyEase::Count)) key.ease = (KeyEase)ease;
                ImGui::SetItemTooltip("How the value moves in the stretch leading to this key");
            }

            ImGui::TableNextColumn();
            if (ImGui::SmallButton("Remove")) remove_idx = i;
            ImGui::PopID();
        }
        ImGui::EndTable();
    }

    if (remove_idx >= 0) {
        m.rep_keys.erase(m.rep_keys.begin() + remove_idx);
        // The property lets go of its value if this was its last key
        movie_reps_apply(data, (double)m.playhead);
    }
    if (resort) movie_rep_sort(data);
}

static void draw_movie_settings_panel(ApplicationState* data) {
    ASSERT(data);
    auto& m = data->movie;
    const bool recording = m.state == MovieRecordingState::Recording;
    char path_buf[2048] = "";

    const double max_frame = (double)(run_num_frames(data) > 0 ? run_num_frames(data) - 1 : 0);
    int frame_w = 0, frame_h = 0;
    movie_frame_size(data, &frame_w, &frame_h);
    m.range_begin = CLAMP(m.range_begin, 0.0f, m.duration);
    m.range_end = CLAMP(m.range_end, m.range_begin, m.duration);
    int range_first = 0, range_last = 0;
    movie_render_range(data, &range_first, &range_last);

    // --- Recording (started from the Output tab) ---
    if (recording) {
        if (ImGui::Button("Stop Recording")) {
            movie_recording_stop(data);
        }
        ImGui::SameLine();
        if (ImGui::Button(m.paused ? "Resume" : "Pause")) {
            m.paused = !m.paused;
        }
        ImGui::SameLine();
        ImGui::Text("Recording frame %d / %d (%.2f s)", m.frame_index - m.rec_first, m.rec_last - m.rec_first + 1, m.cur_time);
        double left;
        if (movie_time_left(data, &left)) {
            char buf[48];
            movie_format_time(buf, sizeof(buf), left);
            ImGui::SameLine();
            ImGui::TextDisabled("about %s left", buf);
        }
    }

    ImGui::BeginDisabled(recording);
    ImGui::BeginDisabled(!m.history.can_undo());
    if (ImGui::Button("Undo")) movie_undo(data);
    ImGui::SetItemTooltip("Undo the last change to the keyframes or the look parameters. Ctrl+Z");
    ImGui::EndDisabled();
    ImGui::SameLine();
    ImGui::BeginDisabled(!m.history.can_redo());
    if (ImGui::Button("Redo")) movie_redo(data);
    ImGui::SetItemTooltip("Ctrl+Y or Ctrl+Shift+Z");
    ImGui::EndDisabled();
    ImGui::EndDisabled();

    // --- Keyframes: what is needed all the time, whatever the tab ---
    ImGui::SeparatorText("Keyframes");
    ImGui::BeginDisabled(recording);
    if (ImGui::Button("Add Keyframe")) {
        movie_add_keyframe_from_view(data);
    }
    if (ImGui::BeginPopupContextItem("add_key_menu")) {
        if (ImGui::MenuItem("Add a key on the path", nullptr, false, md_array_size(m.keyframes) > 0)) movie_add_keyframe_on_path(data);
        ImGui::SetItemTooltip("A key at the preview time where the path already has the camera: the path keeps its shape.\nDouble click in the Camera lane does the same at another time.");
        ImGui::EndPopup();
    }
    ImGui::SetItemTooltip("Adds a keyframe at the preview time with what the viewport shows in the frame, in Scene view as in Movie preview.\n"
        "When there is a key at the preview time already, the new one goes 2 s later and the preview time moves there.\nRight click: a key on the path instead. Shortcut: K");
    ImGui::SameLine();
    ImGui::Checkbox("with trajectory frame", &m.key_includes_frame);
    ImGui::SetItemTooltip("Also key the trajectory frame shown now. Keys with a frame decide how the trajectory plays,\nso the speed can change between them.");
    ImGui::SameLine();
    if (ImGui::Button("Key on Selection")) {
        movie_add_selection_keyframe(data);
    }
    ImGui::SetItemTooltip("Adds a keyframe at the preview time that frames the selected atoms, seen from the direction the camera\nhas now, and moves the view there.");
    ImGui::EndDisabled();
    ImGui::Spacing();

    if (!ImGui::BeginTabBar("Movie settings")) return;
    ImGui::BeginDisabled(recording);

    if (ImGui::BeginTabItem("Output")) {
        if (ImGui::Button("Select Output Folder...")) {
            if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Dir)) {
                size_t path_len = strnlen(path_buf, sizeof(path_buf));
                str_free(m.output_dir, data->allocator.persistent);
                m.output_dir = str_copy({path_buf, path_len}, data->allocator.persistent);
            }
        }
        if (!str_empty(m.output_dir)) {
            ImGui::TextWrapped(STR_FMT, STR_ARG(m.output_dir));
        } else {
            ImGui::TextDisabled("No output folder selected");
        }

        ImGui::InputText("Filename Prefix", m.filename_prefix, sizeof(m.filename_prefix));

        m.output = (MovieOutput)CLAMP((int)m.output, 0, (int)MovieOutput::Count - 1);
        if (ImGui::BeginCombo("Format", movie_output_str[(int)m.output])) {
            for (int i = 0; i < (int)MovieOutput::Count; ++i) {
                if (ImGui::Selectable(movie_output_str[i], i == (int)m.output)) {
                    m.output = (MovieOutput)i;
                }
            }
            ImGui::EndCombo();
        }
        if (movie_output_is_video(m.output)) {
            const int max_crf = m.output == MovieOutput::WebmVp9 ? 63 : 51;
            m.crf = CLAMP(m.crf, 0, max_crf);
            ImGui::SliderInt("Quality (CRF)", &m.crf, 0, max_crf);
            ImGui::SetItemTooltip("Constant rate factor. Lower is better and larger.\nH.264: 18 is close to lossless, 23 is the default. H.265 looks about the same at a CRF 4 to 6 higher, in a\nsmaller file (default 28). VP9: 15 to 35 is the usual range (default 32)."); 
            ImGui::InputText("ffmpeg", m.ffmpeg_path, sizeof(m.ffmpeg_path));
            ImGui::SetItemTooltip("The ffmpeg executable. Found on the PATH if it is only a name.\nFrames are piped straight into it, no image files are written.");
        } else {
            ImGui::TextWrapped("Writes a numbered PNG sequence ('%s_00000.png', ...), to be encoded into a video with ffmpeg.", m.filename_prefix);
            ImGui::BeginDisabled(str_empty(m.output_dir));
            if (ImGui::Button("Copy ffmpeg command")) {
                char cmd[2048];
                // Quoted for a shell; the folder is only copied, never executed from here. libx264 with yuv420p needs even dimensions, hence the scale.
                snprintf(cmd, sizeof(cmd), "ffmpeg -framerate %g -i \"" STR_FMT "/%s_%%05d.png\" -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" -c:v libx264 -pix_fmt yuv420p \"" STR_FMT "/%s.mp4\"",
                    m.fps, STR_ARG(m.output_dir), m.filename_prefix,
                    STR_ARG(m.output_dir), m.filename_prefix);
                ImGui::SetClipboardText(cmd);
            }
            ImGui::EndDisabled();
        }

        m.resolution = (ScreenshotResolution)MIN((int)m.resolution, (int)ScreenshotResolution::Count - 1);
        if (ImGui::BeginCombo("Resolution", screenshot_resolution_str[(int)m.resolution])) {
            for (int i = 0; i < (int)ScreenshotResolution::Count; ++i) {
                if (ImGui::Selectable(screenshot_resolution_str[i], (i == (int)m.resolution))) {
                    m.resolution = (ScreenshotResolution)i;
                }
            }
            ImGui::EndCombo();
        }
        if (m.resolution == ScreenshotResolution::Custom) {
            ImGui::InputInt("Res X", &m.res_x);
            ImGui::InputInt("Res Y", &m.res_y);
            m.res_x = CLAMP(m.res_x, 640, 16384);
            m.res_y = CLAMP(m.res_y, 480, 16384);
        } else if (m.resolution == ScreenshotResolution::Window) {
            ImGui::TextDisabled("%dx%d, the size of the window when recording starts", frame_w, frame_h);
        }

        ImGui::InputFloat("Output FPS", &m.fps, 1.0f, 5.0f, "%.1f");
        m.fps = CLAMP(m.fps, 1.0f, 240.0f);

        static const int scales[] = {100, 75, 50, 25};
        char scale_label[16];
        snprintf(scale_label, sizeof(scale_label), "%d%%", m.res_scale);
        if (ImGui::BeginCombo("Scale", scale_label)) {
            for (int s : scales) {
                char label[16];
                snprintf(label, sizeof(label), "%d%%", s);
                if (ImGui::Selectable(label, s == m.res_scale)) m.res_scale = s;
            }
            ImGui::EndCombo();
        }
        ImGui::SetItemTooltip("The size of the frames as a part of the size above. A smaller one renders much faster: use it to check\na movie before making it in full.");
        if (m.res_scale < 100) {
            ImGui::SameLine();
            ImGui::TextDisabled("%dx%d", frame_w, frame_h);
        }

        ImGui::BeginDisabled(!data->visuals.temporal_aa.enabled);
        ImGui::InputInt("Samples per frame", &m.aa_samples);
        m.aa_samples = CLAMP(m.aa_samples, 0, 256);
        ImGui::EndDisabled();
        ImGui::SetItemTooltip("How many rendered images are averaged into each frame, for smooth edges (temporal anti-aliasing,\nso it has to be on in the render settings). 0 uses the length of the jitter sequence. Fewer is faster.");

        ImGui::Checkbox("Render only a range", &m.range_enabled);
        ImGui::SetItemTooltip("Renders the frames between two times of the movie, e.g. to redo a part of it. With a PNG sequence\nthe files keep the numbers they have in the whole movie, so they can replace the old ones.");
        if (m.range_enabled) {
            float r[2] = {m.range_begin, m.range_end};
            if (ImGui::DragFloatRange2("Range (s)", &r[0], &r[1], 0.05f, 0.0f, m.duration, "from %.2f", "to %.2f")) {
                m.range_begin = r[0];
                m.range_end = r[1];
            }
            if (ImGui::SmallButton("Start at preview time")) m.range_begin = MIN(m.playhead, m.range_end);
            ImGui::SameLine();
            if (ImGui::SmallButton("End at preview time")) m.range_end = MAX(m.playhead, m.range_begin);
            ImGui::SameLine();
            ImGui::TextDisabled("frames %d to %d", range_first, range_last);
        }

        ImGui::Checkbox("Save a workspace copy with the movie", &m.save_copy);
        ImGui::SetItemTooltip("Writes '<prefix>.via' next to the movie, with the camera path, looks, overlays and settings it was made from.\nOpening it and recording again gives the same movie.");

        ImGui::Separator();
        {
            const bool can_start = !str_empty(m.output_dir) && !m.sink;
            ImGui::BeginDisabled(!can_start);
            if (ImGui::Button("Start Recording")) {
                movie_recording_start(data);
            }
            ImGui::EndDisabled();
            ImGui::SameLine();
            if (str_empty(m.output_dir)) {
                ImGui::TextDisabled("Select an output folder first");
            } else if (m.sink) {
                const frame_sink::Status st = frame_sink::status(m.sink);
                ImGui::TextDisabled("Still writing the previous movie, %d frame(s) left", st.queued);
            } else {
                ImGui::TextDisabled("%d of %d frames, %.2f s, %dx%d", range_last - range_first + 1, movie_num_frames(data), movie_duration(data), frame_w, frame_h);
            }
        }
        ImGui::EndTabItem();
    }

    if (ImGui::BeginTabItem("Timing")) {
        double frame_range[2] = { m.start_frame, m.end_frame };
        const double min_frame = 0.0;
        if (ImGui::SliderScalarN("Trajectory Frames", ImGuiDataType_Double, frame_range, 2, &min_frame, &max_frame, "%.0f")) {
            // Start after end plays the trajectory backwards
            m.start_frame = CLAMP(frame_range[0], 0.0, max_frame);
            m.end_frame   = CLAMP(frame_range[1], 0.0, max_frame);
        }

        float len = m.duration;
        if (ImGui::InputFloat("Movie length (s)", &len, 1.0f, 10.0f, "%.2f", ImGuiInputTextFlags_EnterReturnsTrue)) {
            movie_set_duration(data, len);
        }
        ImGui::SetItemTooltip("Everything on the timeline (keyframes, the trajectory, looks and overlays) is scaled with it,\nso the movie keeps its shape. Press Enter to apply.");
        ImGui::SameLine();
        ImGui::TextDisabled("%d frames", movie_num_frames(data));

        float span[2] = { m.traj_begin, m.traj_end };
        if (ImGui::SliderFloat2("Trajectory plays (s)", span, 0.0f, m.duration, "%.2f")) {
            m.traj_begin = CLAMP(span[0], 0.0f, m.duration);
            m.traj_end   = CLAMP(span[1], m.traj_begin, m.duration);
        }
        ImGui::SetItemTooltip("When the trajectory is at its first and at its last frame. Before and after, it is held,\ne.g. to fly over the structure first. They can also be dragged on the timeline.\nKeyframes with a frame in between change its speed.");

        const double span_s = (double)(m.traj_end - m.traj_begin);
        const double anim_fps = fabs((double)data->animation.fps);
        if (span_s > 0.0 && m.end_frame != m.start_frame && anim_fps > 0.0) {
            const double fps_traj = fabs(m.end_frame - m.start_frame) / span_s;
            ImGui::TextDisabled("On average %.1f trajectory frames per second, %.2fx the Animation speed", fps_traj, fps_traj / anim_fps);
            ImGui::SameLine();
            if (ImGui::SmallButton("Match Animation speed")) {
                // The whole movie is scaled so that the trajectory's part of it plays at the Animation panel's speed
                const double want = fabs(m.end_frame - m.start_frame) / anim_fps;
                movie_set_duration(data, (float)(m.duration * want / span_s));
            }
            ImGui::SetItemTooltip("Changes the length of the movie so that the trajectory plays at the Animation panel's speed.");
        }
        ImGui::EndTabItem();
    }
    movie_clamp_anchors(data);

    const float movie_len = (float)movie_duration(data);
    m.playhead = CLAMP(m.playhead, 0.0f, movie_len);

    ImGui::EndDisabled();

    if (ImGui::BeginTabItem("Camera")) {
        ImGui::BeginDisabled(recording);
        {
            bool ticks = (m.path_options & 1) != 0, sight = (m.path_options & 2) != 0, cams = (m.path_options & 4) != 0, rings = (m.path_options & 8) != 0;
            ImGui::TextDisabled("Path in the viewport (Show path at the top)");
            ImGui::SetItemTooltip("The path of the camera (blue) and of what it looks at (yellow), with the camera at each keyframe.\nThe green camera is where the preview time is, with the line it looks along.\n"
                "Each key has a handle on the eye (numbered) and one on what it looks at: click one to go to the key, drag it to edit the key\n"
                "(Ctrl + drag moves both, Esc cancels). Ctrl + click on the path adds a key there.");
            ImGui::BeginDisabled(!movie_path_visible(!m.show_frame, m.show_path, m.path_options));
            ImGui::Indent();
            if (ImGui::Checkbox("Time ticks", &ticks)) m.path_options = (m.path_options & ~1) | (ticks ? 1 : 0);
            ImGui::SetItemTooltip("Dots at round times along the path, with chevrons for the direction. Dots close together are where the camera is slow.");
            ImGui::SameLine();
            if (ImGui::Checkbox("Sight lines", &sight)) m.path_options = (m.path_options & ~2) | (sight ? 2 : 0);
            ImGui::SetItemTooltip("A thin line from the eye to what it looks at at each tick.");
            ImGui::SameLine();
            if (ImGui::Checkbox("Cameras", &cams)) m.path_options = (m.path_options & ~4) | (cams ? 4 : 0);
            ImGui::SetItemTooltip("A camera drawn at each keyframe.");
            ImGui::SameLine();
            if (ImGui::Checkbox("Spin rings", &rings)) m.path_options = (m.path_options & ~8) | (rings ? 8 : 0);
            ImGui::SetItemTooltip("The circle the camera goes round in a spin, with arrows the way it turns.");
            ImGui::Unindent();
            ImGui::EndDisabled();
        }

        ImGui::Checkbox("Seamless loop", &m.loop);
        ImGui::SetItemTooltip("The camera path is cyclic: it moves through the end into the start without a corner.\nFor that the movie has to end in the pose it starts in, 'Close Loop' sets that up.");
        ImGui::SameLine();
        if (ImGui::Button("Close Loop")) {
            movie_close_loop(data);
        }
        ImGui::SetItemTooltip("Ends the movie in the pose of the first keyframe and turns the loop on.");

        ImGui::Checkbox("Keep upright", &m.keep_upright);
        ImGui::SetItemTooltip("The camera stays level, so loops and spins cannot leave it tilted or upside down.\n"
            "The keyframes give where it looks from and at; its tilt is their Roll, 0 by default.\n"
            "Off: the camera goes through the keyframes' own tilts.");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 3.5f);
        ImGui::Combo("Up##movie_up", &m.up_axis, movie_up_axis_str, 6);
        ImGui::SetItemTooltip("The world direction that is up on the screen when the camera is level.");
        ImGui::SameLine();
        if (ImGui::Button("From view##movie_up")) {
            const vec3_t u = data->view.camera.orientation * vec3_t{0, 1, 0};
            int best = 1;
            for (int a = 0; a < 6; ++a) {
                if (vec3_dot(u, movie_up_axes[a]) > vec3_dot(u, movie_up_axes[best])) best = a;
            }
            m.up_axis = best;
        }
        ImGui::SetItemTooltip("Up is the world axis closest to the view's up now.");
        ImGui::BeginDisabled(!m.has_key_clipboard);
        if (ImGui::Button("Paste Keyframe")) {
            movie_paste_keyframe(data);
        }
        ImGui::EndDisabled();
        ImGui::SetItemTooltip("Puts the keyframe that was copied in the table at the preview time (replacing one that is there).\nCtrl+C copies the keyframe at the preview time, Ctrl+V pastes.");

        const float fs = ImGui::GetFontSize();
        ImGui::SetNextItemWidth(fs * 5.5f);
        ImGui::DragInt("##orbit_turns", &m.orbit_turns, 0.05f, -16, 16, "%d turn(s)");
        ImGui::SetItemTooltip("Whole turns around what the camera looks at. Positive is counter-clockwise seen from the tip of the axis.");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(fs * 5.0f);
        ImGui::InputFloat("##orbit_duration", &m.orbit_duration, 0.0f, 0.0f, "%.1f s");
        m.orbit_duration = CLAMP(m.orbit_duration, 0.1f, 3600.0f);
        ImGui::SetItemTooltip("How long the orbit takes");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(fs * 7.0f);
        int orbit_axis = (int)m.orbit_axis;
        if (ImGui::Combo("##orbit_axis", &orbit_axis, spin_axis_str, (int)SpinAxis::Count)) m.orbit_axis = (SpinAxis)orbit_axis;
        ImGui::SetItemTooltip("What the camera turns around: its own up direction, or a world axis");
        ImGui::SameLine();
        if (ImGui::Button("Add Orbit")) {
            movie_add_orbit(data);
        }
        ImGui::SetItemTooltip("Adds a keyframe of the current view at the preview time and another after the orbit's duration,\nwhere the camera is back in the same place after its turns.");
        ImGui::EndDisabled();

        ImGui::EndTabItem();
    }

    if (ImGui::BeginTabItem("Looks")) {
        ImGui::BeginDisabled(recording);
        draw_movie_param_section(data, movie_len);
        ImGui::EndDisabled();
        ImGui::EndTabItem();
    }

    if (ImGui::BeginTabItem("Representations")) {
        ImGui::BeginDisabled(recording);
        draw_movie_rep_section(data, movie_len);
        ImGui::EndDisabled();
        ImGui::EndTabItem();
    }

    if (ImGui::BeginTabItem("Overlays", nullptr, m.editor_select_overlays ? ImGuiTabItemFlags_SetSelected : 0)) {
        m.editor_select_overlays = false;
        ImGui::BeginDisabled(recording);
        draw_movie_overlay_section(data, movie_len);
        ImGui::EndDisabled();
        ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
}

static void draw_movie_window(ApplicationState* data) {
    auto& m = data->movie;
    const double max_frame = (double)(run_num_frames(data) > 0 ? run_num_frames(data) - 1 : 0);
    if (!m.duration_init && max_frame > 0.0) {
        m.start_frame = 0.0;
        m.end_frame = max_frame;
        m.duration = 60.0f;
        m.traj_begin = 0.0f;
        m.traj_end = m.duration;
        m.duration_init = true;
        movie_history_reset(data);
    }
    ImGui::SetNextWindowSize(ImVec2(1280, 800), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(720, 480), ImVec2(FLT_MAX, FLT_MAX));
    const bool movie_open = ImGui::Begin("Movie", &m.show_window, ImGuiWindowFlags_NoFocusOnAppearing);
    {
        // Only a window over the main viewport covers it, not one moved out into a window of its own
        const ImVec2 wp = ImGui::GetWindowPos(), ws = ImGui::GetWindowSize();
        if (ImGui::GetWindowViewport() == ImGui::GetMainViewport()) {
            const ImVec2 o = ImGui::GetMainViewport()->Pos;
            m.window_rect[0] = wp.x - o.x; m.window_rect[1] = wp.y - o.y;
            m.window_rect[2] = wp.x + ws.x - o.x; m.window_rect[3] = wp.y + ws.y - o.y;
            m.window_rect_frame = ImGui::GetFrameCount();
        }
    }
    if (movie_open) {
        ImGui::Checkbox("Timeline", &m.editor_timeline);
        ImGui::SameLine();
        ImGui::Checkbox("Controls", &m.editor_controls);
        if (!m.editor_timeline && !m.editor_controls) m.editor_controls = true;
        ImGui::SameLine();
        if (!m.show_frame) {
            ImGui::BeginDisabled(md_array_size(m.keyframes) == 0 || m.state == MovieRecordingState::Recording);
            if (ImGui::Button("Fit path")) movie_scene_fit_path(data);
            ImGui::EndDisabled();
            ImGui::SetItemTooltip("Move the view so that the whole camera path is in sight, in the part of the viewport this window leaves free.");
            ImGui::SameLine();
        }
        {
            bool shown = movie_path_visible(!m.show_frame, m.show_path, m.path_options);
            if (ImGui::Checkbox("Show path", &shown)) {
                if (!m.show_frame) m.show_path = shown;
                else {
                    m.path_options = (m.path_options & ~16) | (shown ? 16 : 0);
                }
            }
            ImGui::SetItemTooltip(m.show_frame
                ? "Draw the camera path in Movie preview too (Scene view has its own setting). How it is drawn is in the Camera tab."
                : "Draw the camera path in the viewport: blue for the camera, yellow for what it looks at, with a handle on each key.\nHow it is drawn is in the Camera tab.");
        }
        if (!m.show_frame) {
            ImGui::SameLine();
            ImGui::Checkbox("Live picture", &m.pip_enabled);
            ImGui::SetItemTooltip("Scene view: a small live picture of the movie camera in the lower right corner of this window (right click it for the size).\nIt costs one more render of the scene when something changes, and every half second otherwise.");
        }
        const float button_width = ImGui::GetFontSize() * 15.0f;
        const float play_width = ImGui::GetFontSize() * 7.0f;
        const float spacing = ImGui::GetStyle().ItemSpacing.x;
        ImGui::SameLine(MAX(ImGui::GetCursorPosX(), ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x - button_width - play_width - spacing));
        ImGui::BeginDisabled(m.state == MovieRecordingState::Recording || movie_duration(data) <= 0.0);
        ImGui::PushStyleColor(ImGuiCol_Button,        ImVec4(0.15f, 0.55f, 0.25f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.20f, 0.68f, 0.32f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive,  ImVec4(0.12f, 0.45f, 0.20f, 1.0f));
        if (ImGui::Button((const char*)ICON_FA_PLAY " Preview", ImVec2(play_width, ImGui::GetFrameHeight() * 1.4f))) movie_play_mode(data, true);
        ImGui::PopStyleColor(3);
        ImGui::SetItemTooltip("Plays the movie from its start in the viewport as it will be recorded, with every window hidden.\nA bar at the bottom plays, pauses and scrubs it (Space plays and pauses); Exit preview or Esc goes back.");
        ImGui::EndDisabled();
        ImGui::SameLine();
        ImGui::BeginDisabled(m.state == MovieRecordingState::Recording);
        if (ImGui::Button(m.show_frame ? "Switch to Scene view" : "Switch to Movie preview",
                          ImVec2(button_width, ImGui::GetFrameHeight() * 1.4f))) {
            movie_set_scene_view(data, m.show_frame);
        }
        ImGui::SetItemTooltip("Movie preview: the viewport is the movie camera, with the recording frame and the overlays.\nScene view: your own camera, framing the whole camera path (the movie camera is the green one), for editing the path.\nEach keeps its own view. Tab switches; hold Tab to look at the other one for a moment.");
        ImGui::EndDisabled();
        ImGui::Separator();
        draw_movie_preview_controls(data);
        const float available = MAX(1.0f, ImGui::GetContentRegionAvail().y);
        const float divider = 6.0f;
        const float min_pane = ImGui::GetFrameHeight() * 3.0f;
        const float bottom_height = CLAMP(available * m.editor_keyframe_fraction, MIN(min_pane, available * 0.25f), available * 0.65f);
        const float upper_height = MAX(1.0f, available - bottom_height - divider - ImGui::GetStyle().ItemSpacing.y * 2);
        if (ImGui::BeginChild("Movie upper panes", ImVec2(0, upper_height))) {
            const int columns = m.editor_timeline && m.editor_controls ? 2 : 1;
            if (ImGui::BeginTable("Movie editor split bottom keys", columns, ImGuiTableFlags_Resizable | ImGuiTableFlags_SizingStretchProp)) {
                if (m.editor_timeline) ImGui::TableSetupColumn("Timeline", ImGuiTableColumnFlags_WidthStretch, 0.68f);
                if (m.editor_controls) ImGui::TableSetupColumn("Controls", ImGuiTableColumnFlags_WidthStretch, 0.32f);
                if (m.editor_timeline) {
                    ImGui::TableNextColumn();
                    if (ImGui::BeginChild("Movie timeline", ImVec2(0, 0))) draw_movie_timeline_panel(data);
                    ImGui::EndChild();
                }
                if (m.editor_controls) {
                    ImGui::TableNextColumn();
                    if (ImGui::BeginChild("Movie controls", ImVec2(0, 0))) draw_movie_settings_panel(data);
                    ImGui::EndChild();
                }
                ImGui::EndTable();
            }
        }
        ImGui::EndChild();
        const ImVec2 divider_pos = ImGui::GetCursorScreenPos();
        ImGui::InvisibleButton("Movie keyframe divider", ImVec2(-FLT_MIN, divider));
        if (ImGui::IsItemHovered() || ImGui::IsItemActive()) ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeNS);
        if (ImGui::IsItemActive() && available > 0) m.editor_keyframe_fraction =
            CLAMP(m.editor_keyframe_fraction - ImGui::GetIO().MouseDelta.y / available, 0.15f, 0.65f);
        ImGui::GetWindowDrawList()->AddLine(ImVec2(divider_pos.x, divider_pos.y + divider * 0.5f),
            ImVec2(divider_pos.x + ImGui::GetItemRectSize().x, divider_pos.y + divider * 0.5f),
            ImGui::GetColorU32(ImGuiCol_Separator), 2);
        if (ImGui::BeginChild("Movie keyframes bottom", ImVec2(0, 0))) {
            ImGui::SeparatorText("Camera keyframes");
            ImGui::BeginDisabled(m.state == MovieRecordingState::Recording);
            draw_movie_keyframe_table(data, (float)movie_duration(data), m.state == MovieRecordingState::Recording);
            ImGui::EndDisabled();
        }
        ImGui::EndChild();

        // Keep the live picture above the full-width keyframe pane.
        if (m.pip_enabled && movie_scene_view(data)) {
            m.pip_requested_frame = ImGui::GetFrameCount();
            {
                static const float parts[4] = {0.18f, 0.26f, 0.36f, 0.50f};
                static const float least[4] = {200.0f, 280.0f, 360.0f, 480.0f};
                int fw = 0, fh = 0;
                movie_frame_size(data, &fw, &fh);
                const float aspect = (float)MAX(fw, 1) / (float)MAX(fh, 1);
                const ImVec2 ws = ImGui::GetWindowSize();
                const int s = CLAMP(m.pip_size, 0, 3);
                float w = MAX(parts[s] * ws.x, least[s]);
                w = MIN(w, MIN(0.65f * upper_height * aspect, 0.7f * ws.x));
                m.pip_target_w = MAX(64, ((int)w + 8) / 16 * 16);
            }
            if (m.pip_valid && m.pip_tex) {
                const ImVec2 size((float)m.pip_w, (float)m.pip_h);
                const ImVec2 wp = ImGui::GetWindowPos(), ws = ImGui::GetWindowSize();
                const float margin = ImGui::GetStyle().ScrollbarSize + 8.0f;
                ImGui::SetCursorScreenPos(ImVec2(wp.x + ws.x - size.x - margin, divider_pos.y - size.y - margin));
                ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
                ImGui::PushStyleVar(ImGuiStyleVar_ChildBorderSize, 0.0f);
                if (ImGui::BeginChild("##movie_pip", size, ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse | ImGuiWindowFlags_NoBackground)) {
                    ImDrawList* dl = ImGui::GetWindowDrawList();
                    const ImVec2 p0 = ImGui::GetCursorScreenPos();
                    const ImVec2 p1(p0.x + size.x, p0.y + size.y);
                    dl->AddImage((ImTextureID)(intptr_t)m.pip_tex, p0, p1, ImVec2(0, 1), ImVec2(1, 0));
                    if (!m.overlays.empty()) movie_overlays_draw(dl, p0, size, (double)m.playhead, data);
                    dl->AddRect(p0, p1, IM_COL32(255, 255, 255, 160));
                    char label[48];
                    snprintf(label, sizeof(label), "Movie camera  %.2f s", m.playhead);
                    dl->AddRectFilled(ImVec2(p0.x + 3.0f, p0.y + 3.0f), ImVec2(p0.x + ImGui::CalcTextSize(label).x + 9.0f, p0.y + ImGui::GetFontSize() + 7.0f), IM_COL32(0, 0, 0, 120), 3.0f);
                    dl->AddText(ImVec2(p0.x + 6.0f, p0.y + 5.0f), IM_COL32(255, 255, 255, 230), label);
                    ImGui::InvisibleButton("##movie_pip_hit", size);
                    ImGui::SetItemTooltip("What the movie camera sees at the preview time, as it will be recorded (without the look of the final render).\nRight click for the size.");
                    if (ImGui::BeginPopupContextItem("##movie_pip_menu")) {
                        ImGui::TextDisabled("Preview size");
                        if (ImGui::RadioButton("Small", m.pip_size == 0)) m.pip_size = 0;
                        if (ImGui::RadioButton("Medium", m.pip_size == 1)) m.pip_size = 1;
                        if (ImGui::RadioButton("Large", m.pip_size == 2)) m.pip_size = 2;
                        if (ImGui::RadioButton("Extra large", m.pip_size == 3)) m.pip_size = 3;
                        if (ImGui::MenuItem("Hide")) m.pip_enabled = false;
                        ImGui::EndPopup();
                    }
                }
                ImGui::EndChild();
                ImGui::PopStyleVar(2);
            }
        }
    }
    ImGui::End();
}

// Over the viewport while recording, and while the last frames are still being written afterwards
static void draw_movie_recording_banner(ApplicationState* state) {
    auto& m = state->movie;
    const bool recording = m.state == MovieRecordingState::Recording;
    if (!recording && !m.sink) return;

    const frame_sink::Status st = frame_sink::status(m.sink);

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x * 0.5f, vp->WorkPos.y + ImGui::GetFontSize()), ImGuiCond_Always, ImVec2(0.5f, 0.0f));
    ImGui::SetNextWindowBgAlpha(0.88f);
    const ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings |
        ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoNav;
    if (ImGui::Begin("##movie_recording_banner", nullptr, flags)) {
        const float bar_width = ImGui::GetFontSize() * 22.0f;
        if (recording) {
            const int n = m.rec_last - m.rec_first + 1;
            const int done = MIN(m.frame_index - m.rec_first, n);
            ImGui::Text("%s movie, frame %d / %d  (%.2f s)", m.paused ? "Paused" : "Recording", done, n, m.cur_time);
            ImGui::ProgressBar(n > 0 ? (float)done / (float)n : 0.0f, ImVec2(bar_width, 0));
            double left;
            char eta[48] = "";
            if (movie_time_left(state, &left)) {
                char t[40];
                movie_format_time(t, sizeof(t), left);
                snprintf(eta, sizeof(eta), ", about %s left", t);
            }
            ImGui::TextDisabled("%d written, %d waiting to be written%s", st.written, st.queued, eta);
            if (ImGui::Button(m.paused ? "Resume" : "Pause")) {
                m.paused = !m.paused;
            }
            ImGui::SameLine();
            if (ImGui::Button("Stop (Esc)")) {
                movie_recording_stop(state);
            }
        } else {
            ImGui::Text("Writing the movie, %d frame(s) left", st.queued);
            ImGui::ProgressBar(st.submitted > 0 ? (float)(st.written + st.failed) / (float)st.submitted : 1.0f, ImVec2(bar_width, 0));
        }
    }
    ImGui::End();
}

static void draw_coordinate_system_widget_window(ViewTransform* target, const ViewTransform& current) {
    const int def_size = 150;
    const int min_size = 100;
    const int max_size = 500;

    ImGui::SetNextWindowSize(ImVec2(def_size, def_size), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(min_size, min_size), ImVec2(max_size, max_size), [](ImGuiSizeCallbackData* data) {
        // Enforce a square aspect ratio
        data->DesiredSize = ImVec2(ImMax(data->DesiredSize.x, data->DesiredSize.y), ImMax(data->DesiredSize.x, data->DesiredSize.y));
    });

    // Scale down backround color to only show a faint outline
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImGui::GetStyleColorVec4(ImGuiCol_WindowBg) * ImVec4(1.0f, 1.0f ,1.0f ,0.1f));
    defer { ImGui::PopStyleColor(); };

    const bool editable = ImGui::IsKeyDown(ImGuiMod_Alt);

    int window_flags = ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoCollapse;
    if (!editable) {
        window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoBackground | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoNavFocus;
    }

    if (ImGui::Begin("Coordinate Widget", nullptr, window_flags)) {
        quat_t out_orientation = target->orientation;
        if (ImGui::CoordinateSystemWidget(&out_orientation, current.orientation, ImGui::GetContentRegionAvail())) {
            const vec3_t look_at = camera_get_look_at(*target);
            target->orientation = quat_normalize(out_orientation);
            target->position = camera_position_from_look_at(look_at, target->orientation, target->distance);
        }
    }
    ImGui::End();
}

// ## The live preview of the movie camera
//
// While Scene view is shown the Movie window has a small picture of what the movie camera sees at the preview time. It is the scene
// rendered once more through the same pipeline, before the main view so that the main view's picking, history and view matrices are
// the ones that last (those are saved and put back), into the G-buffer and from there scaled down into a small texture. It is only
// rendered when something it depends on changed (at most about 20 times a second) or twice a second otherwise, and only while the
// Movie window is open to show it.

static bool movie_pip_wanted(const ApplicationState* state) {
    const auto& m = state->movie;
    return m.pip_enabled && movie_scene_view(state) && str_empty(state->screenshot.path_to_file) && state->app.window.width > 0 && state->app.window.height > 0 &&
        ImGui::GetFrameCount() - m.pip_requested_frame <= 2;
}

static void movie_pip_free(ApplicationState* state) {
    auto& m = state->movie;
    if (m.pip_fbo) glDeleteFramebuffers(1, &m.pip_fbo);
    if (m.pip_tex) glDeleteTextures(1, &m.pip_tex);
    m.pip_fbo = 0;
    m.pip_tex = 0;
    m.pip_w = m.pip_h = 0;
    m.pip_valid = false;
}

static void movie_render_pip(ApplicationState* state) {
    auto& m = state->movie;
    if (!movie_pip_wanted(state)) return;

    int frame_w = 0, frame_h = 0;
    movie_frame_size(state, &frame_w, &frame_h);
    static const int widths[4] = {240, 360, 480, 720};
    const int w = m.pip_target_w > 0 ? m.pip_target_w : widths[CLAMP(m.pip_size, 0, 3)];
    const int h = MAX(1, (int)((double)w * (double)frame_h / (double)MAX(frame_w, 1) + 0.5));
    if (m.pip_w != w || m.pip_h != h || !m.pip_tex) {
        movie_pip_free(state);
        glGenTextures(1, &m.pip_tex);
        glBindTexture(GL_TEXTURE_2D, m.pip_tex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glBindTexture(GL_TEXTURE_2D, 0);
        glGenFramebuffers(1, &m.pip_fbo);
        glBindFramebuffer(GL_FRAMEBUFFER, m.pip_fbo);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m.pip_tex, 0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        m.pip_w = w;
        m.pip_h = h;
    }

    // What the picture depends on, apart from what the movie's looks do (those make it at most half a second old)
    const size_t n = md_array_size(m.keyframes);
    const double scalars[] = {(double)m.playhead, state->animation.frame, (double)m.loop, (double)m.keep_upright, (double)m.up_axis, (double)w, (double)h,
        (double)state->app.framebuffer.width, (double)state->app.framebuffer.height, (double)n};
    uint64_t hash = md_hash64(scalars, sizeof(scalars), 11);
    if (n > 0) hash = md_hash64_combine(hash, md_hash64(m.keyframes, n * sizeof(CameraKeyframe), 12));
    for (const MovieAtomSet& set : m.atom_sets) {
        hash = md_hash64_combine(hash, md_hash64(set.atoms.data(), set.atoms.size() * sizeof(uint32_t), 13 + set.id));
        hash = md_hash64_combine(hash, md_hash64(set.expr.data(), set.expr.size(), 14 + set.id));
    }
    const double now = ImGui::GetTime();
    const bool changed = hash != m.pip_hash || !m.pip_valid;
    if (!(changed && now - m.pip_time > 0.05) && now - m.pip_time < 0.5) return;

    // The movie camera at the preview time, where the trajectory is now
    ViewTransform vt = state->view.target;
    float fov_y = state->view.camera.fov_y;
    if (n > 0) {
        std::vector<CameraKeyframe> sorted(m.keyframes, m.keyframes + n);
        std::stable_sort(sorted.begin(), sorted.end(), [](const CameraKeyframe& a, const CameraKeyframe& b) { return a.time < b.time; });
        movie_camera_evaluate(state, (double)m.playhead, sorted.data(), n, &vt, &fov_y);
    } else if (m.pose_movie.valid) {
        vt = m.pose_movie.target;
        fov_y = m.pose_movie.fov_y;
    }

    const Camera saved_camera = state->view.camera;
    const ViewParam saved_param = state->view.param;
    state->view.camera = vt;
    state->view.camera.fov_y = fov_y;
    m.pip_pass = true;
    render_scene(state, true);

    // The frame of the movie is the part of the G-buffer that the widened view was made to match
    ImVec2 guide_pos, guide_size;
    if (movie_frame_guide(state, &guide_pos, &guide_size) && guide_size.x > 0.0f && guide_size.y > 0.0f) {
        const float sx = (float)state->gbuffer.width / (float)state->app.window.width;
        const float sy = (float)state->gbuffer.height / (float)state->app.window.height;
        const int x0 = (int)(guide_pos.x * sx), x1 = (int)((guide_pos.x + guide_size.x) * sx);
        const int y0 = (int)((float)state->gbuffer.height - (guide_pos.y + guide_size.y) * sy), y1 = (int)((float)state->gbuffer.height - guide_pos.y * sy);
        glBindFramebuffer(GL_READ_FRAMEBUFFER, state->gbuffer.fbo);
        glReadBuffer(GL_COLOR_ATTACHMENT0);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, m.pip_fbo);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        glBlitFramebuffer(x0, y0, x1, y1, 0, 0, m.pip_w, m.pip_h, GL_COLOR_BUFFER_BIT, GL_LINEAR);
        m.pip_valid = true;
    }
    m.pip_pass = false;
    state->view.camera = saved_camera;
    state->view.param = saved_param;
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    m.pip_hash = hash;
    m.pip_time = now;
}

// The legends of the representations that colour by a value: the field of an isosurface, or the
// attribute of its atoms. One per representation at most, since one representation colours by one
// thing, and only while what it shows exists - a legend of nothing would mislead.
static void draw_color_legend_windows(const ApplicationState& state) {
    const size_t num_reps = md_array_size(state.representation.reps);
    int slot = 0;
    for (size_t i = 0; i < num_reps; ++i) {
        const Representation& rep = state.representation.reps[i];
        if (!rep.enabled) continue;

        if (rep.type == RepresentationType::ElectronicStructure) {
            const ElectronicStructureRepresentation& es = rep.electronic_structure;
            if (es.coloring == SurfaceColoring::Field && es.field_map.show_legend && es.field_vol.tex_id) {
                color_scale_draw_legend(es.field_map, surface_field_unit(es.field_kind), surface_field_kind_str[(int)es.field_kind], rep.name, (int)i, slot++);
                continue;
            }
        }

        if (representation_uses_atom_colors(rep) && rep.color_mapping == ColorMapping::Attribute && rep.atom_attribute.scale.show_legend) {
            const md_attribute_t* attr = md_attributes_get(&state.mold.sys.attributes, rep.atom_attribute.key);
            if (!attr) continue;
            char label[128];
            atom_attribute_legend_label(label, sizeof(label), attr, rep.atom_attribute.variant_idx);
            color_scale_draw_legend(rep.atom_attribute.scale, attr->unit, label, rep.name, (int)i, slot++);
        }
    }
}

static void render(ApplicationState* state) {
    movie_render_pip(state);
    render_scene(state, false);
}

// Renders the scene. With 'pip' it is the live preview of the movie camera: no screenshot, recording, immediate drawing or temporal
// effects, and the result is left in the G-buffer's color attachment for the caller, which renders the main view after it.
static void render_scene(ApplicationState* state, bool pip) {
    // Frames of a recording are rendered at the movie's size into the G-buffer for as long as it lasts
    const bool movie_capture = !pip && state->movie.state == MovieRecordingState::Recording;
    bool do_screenshot = !pip && !str_empty(state->screenshot.path_to_file) && !movie_capture;

    uint32_t gbuffer_target_width  = state->app.framebuffer.width;
    uint32_t gbuffer_target_height = state->app.framebuffer.height;
    if (movie_capture) {
        gbuffer_target_width  = (uint32_t)state->movie.rec_w;
        gbuffer_target_height = (uint32_t)state->movie.rec_h;
    } else if (do_screenshot) {
        gbuffer_target_width  = state->screenshot.res_x;
        gbuffer_target_height = state->screenshot.res_y;
    }

    // Resize Framebuffer
    if ((state->gbuffer.width != gbuffer_target_width || state->gbuffer.height != gbuffer_target_height) &&
        (gbuffer_target_width != 0 && gbuffer_target_height != 0)) {
        gbuffer_init(&state->gbuffer, gbuffer_target_width, gbuffer_target_height);
        postprocess_pipeline::initialize(state->gbuffer.width, state->gbuffer.height);
    }

    update_view_param(state);

    volume::timings_new_frame();

    gbuffer_clear(&state->gbuffer);

    const GLenum draw_buffers_opaque[] = {GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
        GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_TRANSPARENCY };

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    // Enable all draw buffers
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, state->gbuffer.fbo);
    glDrawBuffers((int)ARRAY_SIZE(draw_buffers_opaque), draw_buffers_opaque);

    PUSH_GPU_SECTION("G-Buffer fill")

    if (!pip && state->simulation_box.enabled && state->mold.state.unitcell.flags != 0) {
        mat3_t A = { 0 };
		md_unitcell_A_extract_float(A.elem, &state->mold.state.unitcell);
        immediate::Scope scope(state->gfx.world, "simulation box");
        // Through the turn the coordinates carry, if the orientation is kept: the cell itself cannot hold one
        immediate::box_wireframe(scope, {0,0,0}, {1,1,1}, state->operations.state_rotation * mat4_from_mat3(A), convert_color(state->simulation_box.color));
    }

    if (!pip) {
        immediate::Scope vis_scope(state->gfx.overlay, "visualization");
        immediate::Scope vis_scope_depth(state->gfx.world, "visualization with depth");

        const md_script_vis_t& vis = state->script.vis;

        if (vis.points) {
            vec4_t c = state->script.point_color;
            c.w *= state->movie.vis_fade;
            immediate::points(vis_scope, (immediate::Vertex*)vis.points, md_array_size(vis.points), c);
        }

        if (vis.triangles) {
            vec4_t c = state->script.triangle_color;
            c.w *= state->movie.vis_fade;
            immediate::triangles(vis_scope, (immediate::Vertex*)vis.triangles, md_array_size(vis.triangles), c);
        }

        if (vis.lines) {
            vec4_t c = state->script.line_color;
            c.w *= state->movie.vis_fade;
            immediate::lines(vis_scope, (immediate::Vertex*)vis.lines, md_array_size(vis.lines), c);
        }

        const size_t num_matrices = md_array_size(vis.sdf.matrices);
        const size_t num_structures = md_array_size(vis.sdf.structures);
        const size_t num_sdf_items = MIN(MIN(num_matrices, num_structures), 100);

        if (num_matrices != num_structures) {
            VIAMD_LOG_DEBUG("SDF visualization returned mismatched matrix and structure counts (%zu, %zu)", num_matrices, num_structures);
        }
        if (num_sdf_items > 0) {
            const vec4_t col_x = { 1, 0, 0, 0.7f };
            const vec4_t col_y = { 0, 1, 0, 0.7f };
            const vec4_t col_z = { 0, 0, 1, 0.7f };
            const float ext = vis.sdf.extent * 0.25f;
            const vec3_t box_ext = vec3_set1(vis.sdf.extent);
            for (size_t i = 0; i < num_sdf_items; ++i) {
                const mat4_t model_matrix = mat4_inverse(vis.sdf.matrices[i]);
                immediate::basis(vis_scope, model_matrix, ext, col_x, col_y, col_z);
                immediate::box_wireframe(vis_scope_depth, -box_ext, box_ext, model_matrix);
            }
        }
    }

    if (!use_gfx) {
        // DRAW VELOCITY OF STATIC OBJECTS
        PUSH_GPU_SECTION("Blit Static Velocity")
        glDrawBuffer(GL_COLOR_ATTACHMENT_VELOCITY);
        glDepthMask(0);
        postprocessing::blit_static_velocity(state->gbuffer.tex.depth, state->view.param);
        glDepthMask(1);
        POP_GPU_SECTION()
    }
    glDepthMask(1);
    glColorMask(1, 1, 1, 1);

    // DRAW REPRESENTATIONS
    PUSH_GPU_SECTION("Draw Opaque")
    glDrawBuffers((int)ARRAY_SIZE(draw_buffers_opaque), draw_buffers_opaque);
    draw_representations_opaque(state);
    viamd::event_system_broadcast_event(viamd::EventType_ViamdRenderOpaque, viamd::EventPayloadType_ApplicationState, state);
    POP_GPU_SECTION()

    immediate::RenderParams params = {};
    params.view = state->view.param.matrix.curr.view;
    params.proj = state->view.param.matrix.curr.proj;

    // standard alpha blending
    PUSH_GPU_SECTION("Immediate world")
    //glEnable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    if (!pip) immediate::render(state->gfx.world, params);
    //glDisable(GL_BLEND);
    POP_GPU_SECTION()

    glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);

    if (!use_gfx) {
        PUSH_GPU_SECTION("Selection")
        const bool atom_selection_empty = md_bitfield_empty(&state->selection.selection_mask);
        const bool atom_highlight_empty = md_bitfield_empty(&state->selection.highlight_mask);

        glDepthMask(0);

        // @NOTE(Robin): This is a b*tch to get right, What we want is to separate in a single pass, the visible selected from the
        // non visible selected. In order to achieve this, we start with a cleared stencil of value 1 then either set it to zero selected and not visible
        // and to two if it is selected and visible. But the visible atoms should always be able to write over a non visible 0, but not the other way around.
        // Hence the GL_GREATER stencil test against the reference value of 2.

        if (!atom_selection_empty) {
            glColorMask(0, 0, 0, 0);

            glEnable(GL_DEPTH_TEST);
            glDepthFunc(GL_EQUAL);

            glEnable(GL_STENCIL_TEST);
            glStencilMask(0xFF);

            glClearStencil(1);
            glClear(GL_STENCIL_BUFFER_BIT);

            glStencilFunc(GL_GREATER, 0x02, 0xFF);
            glStencilOp(GL_KEEP, GL_ZERO, GL_REPLACE);
            draw_representations_opaque_lean_and_mean(state, AtomBit_Selected | AtomBit_Visible);

            glDisable(GL_DEPTH_TEST);

            glStencilMask(0x0);
            glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
            glColorMask(1, 1, 1, 1);

            glStencilFunc(GL_EQUAL, 2, 0xFF);
            postprocessing::blit_color(state->selection.color.selection.visible);

            glStencilFunc(GL_EQUAL, 0, 0xFF);
            postprocessing::blit_color(state->selection.color.selection.hidden);
        }

        if (!atom_highlight_empty) {
            glColorMask(0, 0, 0, 0);

            glEnable(GL_DEPTH_TEST);
            glDepthFunc(GL_EQUAL);

            glEnable(GL_STENCIL_TEST);
            glStencilMask(0xFF);

            glClearStencil(1);
            glClear(GL_STENCIL_BUFFER_BIT);

            glStencilFunc(GL_GREATER, 0x02, 0xFF);
            glStencilOp(GL_KEEP, GL_ZERO, GL_REPLACE);
            draw_representations_opaque_lean_and_mean(state, AtomBit_Highlighted | AtomBit_Visible);

            glDisable(GL_DEPTH_TEST);

            glStencilMask(0x0);
            glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
            glColorMask(1, 1, 1, 1);

            glStencilFunc(GL_EQUAL, 2, 0xFF);
            vec4_t col_vis = state->selection.color.highlight.visible;
            col_vis.w += sin(ImGui::GetTime() * HIGHLIGHT_PULSE_TIME_SCALE) * HIGHLIGHT_PULSE_ALPHA_SCALE;
            col_vis.w *= state->movie.vis_fade;   // A property of the movie fades in and out

            postprocessing::blit_color(col_vis);

            glStencilFunc(GL_EQUAL, 0, 0xFF);
            vec4_t col_hidden = state->selection.color.highlight.hidden;
            col_hidden.w *= state->movie.vis_fade;
            postprocessing::blit_color(col_hidden);
        }

        glDisable(GL_STENCIL_TEST);

        if (!atom_selection_empty) {
            PUSH_GPU_SECTION("Desaturate")
            const float saturation = state->selection.color.saturation;
            glDrawBuffer(GL_COLOR_ATTACHMENT_COLOR);
            postprocessing::scale_hsv(state->gbuffer.tex.color, vec3_t{1, saturation, 1});
            POP_GPU_SECTION()
        }

        glDepthFunc(GL_LESS);
        glDepthMask(0);
        glColorMask(1,1,1,1);
        POP_GPU_SECTION()
    }

    glDisable(GL_STENCIL_TEST);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glColorMask(1, 1, 1, 1);

    glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);

    const bool transparency_hdr_written = draw_representations_transparent(state);
    if (transparency_hdr_written) {
        // What the transparency buffer holds so far (the selection and highlight tints) lies under the
        // isosurfaces, which are composited before it, in HDR: cover it by their coverage here. The
        // overlays drawn after this stay on top of everything.
        PUSH_GPU_SECTION("Isosurface coverage")
        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_ZERO, GL_ONE_MINUS_SRC_ALPHA);
        postprocessing::blit_texture(state->gbuffer.tex.transparency_hdr);
        glDisable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
        POP_GPU_SECTION()
    }
    viamd::event_system_broadcast_event(viamd::EventType_ViamdRenderTransparent, viamd::EventPayloadType_ApplicationState, state);

    const GLenum draw_buffers_transparent[] = { GL_COLOR_ATTACHMENT_TRANSPARENCY, 0, GL_COLOR_ATTACHMENT_VELOCITY, GL_COLOR_ATTACHMENT_PICKING };
	glDrawBuffers(ARRAY_SIZE(draw_buffers_transparent), draw_buffers_transparent);

    glDisable(GL_CULL_FACE);

    // Blend colour only. glEnable(GL_BLEND) applies to every attachment, and
    // this pass also writes velocity and picking: alpha-blending those means a
    // translucent overlay fragment scales the velocity underneath it by
    // (1 - alpha) and mixes the picking index with whatever was already there,
    // neither of which is a meaningful operation on that data.
    glEnablei(GL_BLEND, 0);
    // The transparency buffer holds premultiplied colour: straight alpha sources, accumulated premultiplied
    glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    glDisablei(GL_BLEND, 1);
    glDisablei(GL_BLEND, 2);
    glDisablei(GL_BLEND, 3);

    PUSH_GPU_SECTION("Immediate overlay")
    glDisable(GL_DEPTH_TEST);
    if (!pip) immediate::render(state->gfx.overlay, params);
    POP_GPU_SECTION()

    glDisable(GL_BLEND);
    glEnable(GL_CULL_FACE);

    POP_GPU_SECTION()  // G-buffer

    if (movie_capture || pip || (do_screenshot && state->screenshot.hide_gui)) {
        // Activate gbuffer to store the frame
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, state->gbuffer.fbo);
        glViewport(0, 0, state->gbuffer.width, state->gbuffer.height);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
    } else {
        // Activate backbuffer
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glViewport(0, 0, state->app.framebuffer.width, state->app.framebuffer.height);
        glDrawBuffer(GL_BACK);
        glClear(GL_COLOR_BUFFER_BIT);
    }

    PUSH_GPU_SECTION("Postprocessing")
    postprocess_pipeline::Settings settings = {};
    postprocess_pipeline::Inputs inputs = {};

    settings.background_color = state->visuals.background.color * state->visuals.background.intensity;

    settings.ssao.enabled = state->visuals.ssao.enabled;
    settings.ssao.intensity = state->visuals.ssao.intensity;
    settings.ssao.temporal = state->visuals.ssao.temporal;
    settings.ssao.mode = state->visuals.ssao.mode;

    settings.tonemap.enabled = state->visuals.tonemapping.enabled;
    settings.tonemap.mode = state->visuals.tonemapping.tonemapper;
    settings.tonemap.exposure = state->visuals.tonemapping.exposure;
    settings.tonemap.gamma = state->visuals.tonemapping.gamma;

    settings.dof.enabled = state->visuals.dof.enabled;
    state->visuals.dof.focus_depth = dof_focus_depth(state, state->view.camera);
    settings.dof.focus_depth = state->visuals.dof.focus_depth;
    settings.dof.aperture = state->visuals.dof.aperture * 0.01f;
    if (movie_focus_active(state)) {
        // The keys that set the focus decide what is sharp; the blur is the look parameter
        const auto& m = state->movie;
        const MovieResolved r = movie_resolve(state, m.keyframes, md_array_size(m.keyframes));
        settings.dof.focus_depth = MAX(camera_depth_of_point(state->view.camera, movie_focus_evaluate(r, movie_focus_time(state), state->view.camera)), 1.0e-3f);
        settings.dof.enabled = state->visuals.dof.aperture > 0.0f;
    } else if (movie_focus_owned(state)) {
        settings.dof.aperture = state->movie.editor_saved_aperture * 0.01f;
    }

    settings.fxaa.enabled = state->visuals.fxaa.enabled;

    constexpr float MOTION_BLUR_REFERENCE_DT = 1.0f / 60.0f;
    const float dt_compensation = MOTION_BLUR_REFERENCE_DT / (float)state->app.timing.delta_s;
    const float motion_scale = state->visuals.temporal_aa.motion_blur.motion_scale * dt_compensation;
    settings.taa.enabled = state->visuals.temporal_aa.enabled;
    settings.taa.feedback_min = state->visuals.temporal_aa.feedback_min;
    settings.taa.feedback_max = state->visuals.temporal_aa.feedback_max;
    settings.taa.motion_blur.enabled = state->visuals.temporal_aa.motion_blur.enabled;
    settings.taa.motion_blur.motion_scale = motion_scale;

    settings.sharpen.enabled = state->visuals.temporal_aa.enabled && state->visuals.sharpen.enabled;
    settings.sharpen.weight = state->visuals.sharpen.weight;

    if (pip) {
        settings.taa.enabled = false;
        settings.taa.motion_blur.enabled = false;
        settings.sharpen.enabled = false;
    }

    inputs.depth = state->gbuffer.tex.depth;
    inputs.color = state->gbuffer.tex.color;
    inputs.normal = state->gbuffer.tex.normal;
    inputs.velocity = state->gbuffer.tex.velocity;
    inputs.transparency = state->gbuffer.tex.transparency;
    inputs.transparency_hdr = transparency_hdr_written ? state->gbuffer.tex.transparency_hdr : 0;
    inputs.history = settings.taa.enabled ? state->gbuffer.tex.history : 0;
    inputs.history_prev = settings.taa.enabled ? state->gbuffer.tex.history_prev : 0;

    postprocess_pipeline::execute(inputs, settings, state->view.param);
    POP_GPU_SECTION()

    if (settings.taa.enabled) {
        // Ping-pong: what was written this frame is read as the previous history next frame.
        uint32_t tmp = state->gbuffer.tex.history;
        state->gbuffer.tex.history = state->gbuffer.tex.history_prev;
        state->gbuffer.tex.history_prev = tmp;
    }

    if (pip) return;

    if (movie_capture) {
        // Samples of one output frame are accumulated by temporal AA, the last of them is the frame
        auto& m = state->movie;
        if (m.can_capture) {
            m.samples_done += 1;
            if (m.samples_done >= m.sample_target) {
                movie_overlays_render(state, m.cur_time);
                movie_capture_frame(state);
                m.frame_index += 1;
                m.samples_done = 0;
            }
        }
        movie_blit_preview(state);
    }

    if (do_screenshot && state->screenshot.hide_gui) {
        state->screenshot.sample_count += 1;

        ImGui::OpenPopup("Capturing Frame");
        ImVec2 size = ImGui::CalcTextSize("Capturing Frame");
        ImGui::SetNextWindowSize(size * ImVec2(2, 5), ImGuiCond_Always);
        if (ImGui::BeginPopupModal("Capturing Frame", 0, ImGuiWindowFlags_NoResize)) {
            float fraction = (float)state->screenshot.sample_count / (float)state->screenshot.sample_target;
            ImGui::ProgressBar(fraction);
            ImGui::EndPopup();
        }

        if (state->screenshot.sample_count == state->screenshot.sample_target) {
            create_screenshot(state->screenshot.path_to_file);
            str_free(state->screenshot.path_to_file, state->allocator.persistent);
            state->screenshot.path_to_file = {};
            state->screenshot.sample_count  = 0;
            state->screenshot.sample_target = 0;
            ImGui::CloseCurrentPopup();
        }

        // Activate backbuffer
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glViewport(0, 0, state->app.framebuffer.width, state->app.framebuffer.height);
        glDrawBuffer(GL_BACK);
        glClear(GL_COLOR_BUFFER_BIT);
    }

    PUSH_GPU_SECTION("Imgui render")
    application::render_imgui(&state->app);
    POP_GPU_SECTION()

    if (do_screenshot && !state->screenshot.hide_gui) {
        create_screenshot(state->screenshot.path_to_file);
        str_free(state->screenshot.path_to_file, state->allocator.persistent);
        state->screenshot.path_to_file = {};
    }
}

static void draw_representations_opaque(ApplicationState* state) {
    ASSERT(state);

    if (state->mold.sys.atom.count == 0) {
        return;
    }

#if EXPERIMENTAL_GFX_API
    if (use_gfx) {
        const uint32_t instance_count = 10000;
        static mat4_t* transforms = 0;

        if (transforms == 0) {
            auto rnd = []() -> float {
                return (float)rand() / RAND_MAX;
            };
            for (uint32_t i = 0; i < instance_count; ++i) {
                vec3_t axis = {rnd(), rnd(), rnd()};
                quat_t ori = quat_angle_axis(rnd() * TWO_PI, vec3_normalize(axis));
                mat4_t R = mat4_from_quat(ori);
                mat4_t T = mat4_translate(rnd() * 4000, rnd() * 4000, rnd() * 4000);
                mat4_t M = T * R;
                md_urange_t range = {0, (int32_t)state->mold.sys.atom.count};
                md_array_push(transforms, M, persistent_alloc);
            }
        }

        md_gfx_draw_op_t* draw_ops = 0;
        for (int64_t i = 0; i < md_array_size(state->representation.reps); ++i) {
            if (state->representation.reps[i].enabled) {
                md_gfx_draw_op_t op;
                op.structure = state->mold.gfx_structure;
                op.representation = state->representation.reps[i].gfx_rep;
                op.model_mat = NULL;
                md_array_push(draw_ops, op, frame_alloc);
                
                for (uint32_t j = 0; j < instance_count; ++j) {
                    md_gfx_draw_op_t op;
                    op.structure = state->mold.gfx_structure;
                    op.representation = state->representation.reps[i].gfx_rep;
                    op.model_mat = &transforms[j];
                    md_array_push(draw_ops, op, frame_alloc);
                }
                
            }
        }

        md_gfx_draw((uint32_t)md_array_size(draw_ops), draw_ops, &state->view.param.matrix.curr.proj, &state->view.param.matrix.curr.view, &state->view.param.matrix.inv.proj, &state->view.param.matrix.inv.view);
    } else {
#endif
        const size_t num_representations = md_array_size(state->representation.reps);
        if (num_representations == 0) return;

		glEnable(GL_DEPTH_TEST);

        md_array(md_gl_draw_op_t) draw_ops = 0;
        for (size_t i = 0; i < num_representations; ++i) {
            const Representation& rep = state->representation.reps[i];

            if (RepresentationType::SpaceFill <= rep.type && rep.type <= RepresentationType::Cartoon) {
                if (rep.enabled && rep.type_is_valid) {
                    md_gl_draw_op_t op = {
                        .type = (md_gl_rep_type_t)rep.type,
                        .args = {},
                        .rep = state->representation.reps[i].md_rep,
                        .model_matrix = NULL,
                    };
                    // While it is shown or hidden by the movie it grows in or shrinks away
                    const vec4_t sc = {rep.scale.x * rep.presence, rep.scale.y * rep.presence, rep.scale.z * rep.presence, rep.scale.w * rep.presence};
                    switch (rep.type) {
                    case RepresentationType::SpaceFill:
                        op.args.space_fill.radius_scale = sc.x;
                        break;
                    case RepresentationType::Licorice:
                        op.args.licorice.radius = sc.x;
                        op.args.licorice.color_mode = (md_gl_bond_mode_t)rep.bond_color;
                        op.args.licorice.sharpness = rep.bond_sharpness;
                        op.args.licorice.uniform_color = convert_color(rep.bond_base_color);
                        if (rep.tint_scale > 0.0f || rep.saturation < 1.0f) {
                            tint_colors(&op.args.licorice.uniform_color, 1, convert_color(rep.tint_color), rep.tint_scale, rep.saturation);
                        }
                        break;
                    case RepresentationType::BallAndStick:
                        op.args.ball_and_stick.ball_scale = sc.x;
                        op.args.ball_and_stick.stick_radius = sc.y;
                        op.args.ball_and_stick.color_mode = (md_gl_bond_mode_t)rep.bond_color;
                        op.args.ball_and_stick.sharpness = rep.bond_sharpness;
                        op.args.ball_and_stick.uniform_color = convert_color(rep.bond_base_color);
                        if (rep.tint_scale > 0.0f || rep.saturation < 1.0f) {
                            tint_colors(&op.args.ball_and_stick.uniform_color, 1, convert_color(rep.tint_color), rep.tint_scale, rep.saturation);
                        }
                        break;
                    case RepresentationType::Ribbons:
                        op.args.ribbons.width_scale = sc.x;
                        op.args.ribbons.thickness_scale = sc.y;
                        break;
                    case RepresentationType::Cartoon:
                        op.args.cartoon.coil_scale = sc.x;
                        op.args.cartoon.sheet_scale = sc.y;
                        op.args.cartoon.helix_scale = sc.z;
                        break;
                    default:
                        break;
                    }
                    md_array_push(draw_ops, op, frame_alloc);
                }
            } else if (rep.type == RepresentationType::DipoleMoment) {
                if (rep.enabled) {
                    // immediate draw of dipole moment as arrow
                    vec3_t dipole_vec = { 0, 0, 0 };
                    vec3_t dipole_org = { 0, 0, 0 };
                    if (dipole_moment_read(&dipole_vec, &dipole_org, state->mold.sys, rep.dipole.dipole_key, rep.dipole.dipole_index)) {
                        // The representation's own (key, index) IS the picking address: the group's
                        // range was reserved under that key this frame, and the shader adds the base to
                        // the per primitive index, so the element index is what goes on the primitive.
                        // A group that did not fit in the picking space draws with INVALID_PICKING_IDX,
                        // which the shader passes through untouched - visible, just not pickable.
                        const PickingSpace* space = picking_handler_current_space(&state->picking_handler);
                        const PickingRange* range = space ? picking_space_find_range(*space, PickingDomain_Dipole, rep.dipole.dipole_key) : nullptr;

                        immediate::Scope scope(state->gfx.world, "debug_dipole_moment");
                        immediate::set_picking_base_idx(scope, range ? range->beg : 0);

                        const vec3_t org = dipole_org;
                        const vec3_t vec = dipole_vec * (float)rep.dipole.scale;

                        // cylinder body
                        const float body_scale = 0.8f;

                        const float body_radius = rep.dipole.radius;
                        const float head_radius = body_radius * 1.5f;

                        vec3_t cyl_beg = rep.dipole.offset + org;
                        vec3_t cyl_end = rep.dipole.offset + org + vec * body_scale;
                        vec3_t arrow_end = rep.dipole.offset + org + vec;

                        uint32_t color_u32 = convert_color(rep.dipole.color);

                        uint32_t picking_idx = range ? rep.dipole.dipole_index : INVALID_PICKING_IDX;

                        immediate::cylinder(scope, cyl_beg, cyl_end, body_radius, color_u32, picking_idx);
                        immediate::cone(scope, cyl_end, arrow_end, head_radius, color_u32, picking_idx);
                    }
                }
            }
        }

		// MIN HALF BOX EXTENT AS MAX BOND LENGTH APPROXIMATION
        float max_bond_length = 10.0f;
        if (md_unitcell_flags(&state->mold.state.unitcell) == 0) {
            vec3_t aabb_ext = vec3_sub(state->mold.sys_aabb_max, state->mold.sys_aabb_min);
            max_bond_length = MAX(3.0f, vec3_length(aabb_ext) * 0.5f); // Max bond length should not exceed half the diagonal of the bounding box
        }
        else {
            vec3_t ext = { 0 };
            md_unitcell_diag_extract_float(ext.elem, &state->mold.state.unitcell);
            float min_half_box_ext = vec3_reduce_min(vec3_mul1(ext, 0.5f));
            max_bond_length = MAX(3.0, min_half_box_ext);
        }

        md_gl_draw_args_t args = {
            .shaders = state->gl.shaders,
            .draw_operations = {
                .count = (uint32_t)md_array_size(draw_ops),
                .ops = draw_ops,
            },
            .view_transform = {
                .view_matrix = &state->view.param.matrix.curr.view.elem[0][0],
                .proj_matrix = &state->view.param.matrix.curr.proj.elem[0][0],
                // These two are for temporal anti-aliasing reprojection (optional)
                .prev_view_matrix = &state->view.param.matrix.prev.view.elem[0][0],
                .prev_proj_matrix = &state->view.param.matrix.prev.proj.elem[0][0],
            },
            .picking_offset = gl_picking_offset(*state),
            .max_bond_length = max_bond_length,
        };

        md_gl_draw(&args);
#if EXPERIMENTAL_GFX_API
    }
#endif
}

static bool draw_representations_transparent(ApplicationState* state) {
    ASSERT(state);
    if (state->mold.sys.atom.count == 0) return false;

    const size_t num_representations = md_array_size(state->representation.reps);
    if (num_representations == 0) return false;

    bool written = false;

    for (size_t i = 0; i < num_representations; ++i) {
        Representation& rep = state->representation.reps[i];
        if (!rep.enabled) continue;
        if (rep.type == RepresentationType::ElectronicStructure) {
            const bool use_field = rep.electronic_structure.coloring == SurfaceColoring::Field && rep.electronic_structure.field_vol.tex_id;
            IsoDesc iso;
            electronic_structure_iso_desc_init(&iso, rep.electronic_structure);

#if VIAMD_RECOMPUTE_ORBITAL_PER_FRAME
            flag_representation_as_dirty(&state->representation.reps[i]);
#endif

            volume::IsoRenderDesc desc = {
                .render_target = {
                    .depth = state->gbuffer.tex.depth,
                    .color = state->gbuffer.tex.transparency_hdr,
                    .width = state->gbuffer.width,
                    .height = state->gbuffer.height,
                    .clear_color = !written,
                },
                .texture = {
                    .density_volume = rep.electronic_structure.density_vol.tex_id,
                    .color_volume = rep.electronic_structure.color_vol.tex_id,
                    .field_volume = use_field ? rep.electronic_structure.field_vol.tex_id : 0,
                    .field_colormap = use_field ? surface_field_colormap_texture(&rep.electronic_structure.field_vol, rep.electronic_structure.field_map.colormap) : 0,
                },
                .matrix = {
                    .model = rep.electronic_structure.density_vol.texture_to_world,
                    .view = state->view.param.matrix.curr.view,
                    .proj = state->view.param.matrix.curr.proj,
                    .inv_proj = state->view.param.matrix.inv.proj,
                },
                .clip_volume = {
                    .min = {0,0,0},
                    .max = {1,1,1},
                },
                .iso = {
                    .count = iso.count,
                    .values = iso.values,
                    .colors = iso.colors,
                    .optical_densities = iso.optical_densities,
                    .use_color_volume = rep.electronic_structure.coloring == SurfaceColoring::AtomColors,
                    .use_field = use_field,
                    .exact = state->settings.exact_isosurfaces,
                },
                .field = {
                    .range_beg = rep.electronic_structure.field_map.range_beg,
                    .range_end = rep.electronic_structure.field_map.range_end,
                },
                .shading = {
                    .env_radiance = state->visuals.background.color * state->visuals.background.intensity * 0.25,
                    .roughness = 0.3f,
                    .dir_radiance = {10,10,10},
                    .ior = 1.5f,
                },
            };

            if (volume::render_isosurfaces(desc)) {
                written = true;
            }

#if DEBUG
            {
				immediate::Scope scope(state->gfx.world, "debug_electronic_structure");
                immediate::box_wireframe(scope, { 0,0,0 }, { 1,1,1 }, rep.electronic_structure.density_vol.texture_to_world, immediate::COLOR_BLACK);
            }
#endif
        }
    }
    return written;
}

static void draw_representations_opaque_lean_and_mean(ApplicationState* data, uint32_t mask) {
    md_gl_draw_op_t* draw_ops = 0;
    for (size_t i = 0; i < md_array_size(data->representation.reps); ++i) {
        const Representation& rep = data->representation.reps[i];

        if (rep.type > RepresentationType::Cartoon) continue;

        if (rep.enabled && rep.type_is_valid) {
            md_gl_draw_op_t op = {
                .type = (md_gl_rep_type_t)rep.type,
                .args = {},
                .rep = data->representation.reps[i].md_rep,
                .model_matrix = NULL,
            };
            const vec4_t sc = {rep.scale.x * rep.presence, rep.scale.y * rep.presence, rep.scale.z * rep.presence, rep.scale.w * rep.presence};
            MEMCPY(&op.args, &sc, sizeof(op.args));
            md_array_push(draw_ops, op, frame_alloc);
        }
    }

    // MIN HALF BOX EXTENT AS MAX BOND LENGTH APPROXIMATION
    float max_bond_length = 10.0f;
    if (md_unitcell_flags(&data->mold.state.unitcell) == 0) {
        vec3_t aabb_ext = vec3_sub(data->mold.sys_aabb_max, data->mold.sys_aabb_min);
        max_bond_length = MAX(3.0f, vec3_length(aabb_ext) * 0.5f); // Max bond length should not exceed half the diagonal of the bounding box
    }
    else {
        vec3_t ext = { 0 };
        md_unitcell_diag_extract_float(ext.elem, &data->mold.state.unitcell);
        float min_half_box_ext = vec3_reduce_min(vec3_mul1(ext, 0.5f));
        max_bond_length = MAX(3.0, min_half_box_ext);
    }

    md_gl_draw_args_t args = {
        .shaders = data->gl.shaders_lean_and_mean,
        .draw_operations = {
            .count = md_array_size(draw_ops),
            .ops = draw_ops,
        },
        .view_transform = {
            .view_matrix = &data->view.param.matrix.curr.view.elem[0][0],
            .proj_matrix = &data->view.param.matrix.curr.proj.elem[0][0],
            // These two are for temporal anti-aliasing reprojection
            //.prev_model_view_matrix = &data->view.param.matrix.previous.view[0][0],
            //.prev_projection_matrix = &data->view.param.matrix.previous.proj[0][0],
        },
        .picking_offset = gl_picking_offset(*data),
        .atom_mask = mask,
		.max_bond_length = max_bond_length,
    };

    md_gl_draw(&args);
}
