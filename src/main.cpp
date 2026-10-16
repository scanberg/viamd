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
#include <core/md_writer.h>
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

#include <implot.h>
#include <implot_internal.h>

#include <stdio.h>
#include <cmath>
#include <string>

#include <viamd.h>
#include <representation.h>
#include <timeline_periods_ui.h>
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


// What Evaluate would do with the script as it is in the editor (md_script_eval_plan): how each of its properties
// stands against what is evaluated, and so how what is shown stands against what the script now says
struct ScriptPlan {
    bool valid = false;         // the editor's script compiles, and there are frames to evaluate
    bool running = false;       // an evaluation is running
    bool same = false;          // the editor's script is the one evaluated
    bool changes = false;       // Evaluate would change which properties there are, or what they are called
    const md_script_eval_plan_t* props = nullptr;   // one per property of the editor's script
    const str_t* names = nullptr;
    size_t   num_props = 0;
    uint32_t num_frames = 0;
    uint32_t num_pending = 0;   // the frames Evaluate would evaluate
};

static ScriptPlan script_plan(const ApplicationState& state, md_allocator_i* alloc) {
    ScriptPlan plan;
    const md_script_ir_t* ir = state.script.ir;
    plan.num_frames = (uint32_t)run_num_frames(&state);
    plan.valid = md_script_ir_valid(ir) && plan.num_frames > 0;
    plan.running = task_system::task_is_running(state.tasks.evaluate);
    if (!plan.valid) return plan;

    plan.same = state.script.eval && ir == state.script.eval_ir;
    plan.num_props = md_script_ir_property_count(ir);
    plan.names = md_script_ir_property_names(ir);
    md_script_eval_plan_t* props = (md_script_eval_plan_t*)md_alloc(alloc, MAX(plan.num_props, 1) * sizeof(md_script_eval_plan_t));
    md_bitset_t pending = md_bitset_empty_set();
    if (!md_script_eval_plan(props, plan.num_props, &pending, state.script.eval, ir, plan.num_frames, alloc)) {
        plan.valid = false;
        return plan;
    }
    plan.props = props;
    plan.num_pending = (uint32_t)md_bitset_count(pending);

    // A property gone, or renamed, changes what is shown without anything to evaluate
    const size_t num_evaluated = state.script.eval ? md_script_eval_property_count(state.script.eval) : 0;
    plan.changes = num_evaluated != plan.num_props;
    for (size_t i = 0; i < plan.num_props; ++i) {
        plan.changes |= props[i].status != MD_SCRIPT_EVAL_STATUS_HELD || !str_eq(props[i].from, plan.names[i]);
    }
    return plan;
}

// How property i of the plan stands, as a line of text
static void script_plan_describe(char* buf, size_t cap, const ScriptPlan& plan, size_t i) {
    const md_script_eval_plan_t& p = plan.props[i];
    const str_t name = plan.names[i];
    switch (p.status) {
    case MD_SCRIPT_EVAL_STATUS_NEW:
        snprintf(buf, cap, STR_FMT ": not evaluated yet, %u frames to evaluate", STR_ARG(name), p.num_pending);
        break;
    case MD_SCRIPT_EVAL_STATUS_CHANGED:
        snprintf(buf, cap, STR_FMT ": changed since it was evaluated, %u frames to evaluate", STR_ARG(name), p.num_pending);
        break;
    case MD_SCRIPT_EVAL_STATUS_HELD: {
        char from[96] = "";
        if (!str_eq(p.from, name)) snprintf(from, sizeof(from), ", computed as " STR_FMT " was", STR_ARG(p.from));
        if (p.num_pending == 0) {
            snprintf(buf, cap, STR_FMT ": up to date%s", STR_ARG(name), from);
        } else if (p.num_frames < plan.num_frames) {
            snprintf(buf, cap, STR_FMT ": %u of %u frames evaluated%s", STR_ARG(name), p.num_frames, plan.num_frames, from);
        } else {
            snprintf(buf, cap, STR_FMT ": %u frames to add or take out over periods%s", STR_ARG(name), p.num_pending, from);
        }
        break;
    }
    }
}

// A mark beside the line a property is named on: how it stands. It is also a handle on the property, which is dragged
// from it into the window which shows its kind.
struct ScriptGutterMark {
    int line;
    md_script_eval_status_t status;
    float done;                 // the fraction of what it needs which it holds
    char text[160];
    char ident[64];
    uint32_t kinds;             // MD_SCRIPT_PROPERTY_FLAG_*
};

static void draw_script_gutter_mark(ApplicationState* app, const TextEditor::Decorator& d, const ScriptGutterMark& m) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    const ImVec2 c  = p0 + ImVec2(d.width * 0.5f, d.height * 0.5f);
    const float  r  = d.glyphSize.y * 0.26f;
    const ImU32 green = IM_COL32(90, 200, 110, 255);
    const ImU32 amber = IM_COL32(240, 170, 60, 255);
    const ImU32 gray  = IM_COL32(150, 150, 150, 255);
    switch (m.status) {
    case MD_SCRIPT_EVAL_STATUS_NEW:
        dl->AddCircle(c, r, gray, 0, 1.5f);
        break;
    case MD_SCRIPT_EVAL_STATUS_CHANGED:
        dl->AddCircleFilled(c, r, amber);
        break;
    case MD_SCRIPT_EVAL_STATUS_HELD:
        if (m.done >= 1.0f) {
            dl->AddCircleFilled(c, r, green);
        } else {
            // How much of it there is, as a slice
            dl->AddCircle(c, r, green, 0, 1.5f);
            if (m.done > 0.0f) {
                dl->PathLineTo(c);
                dl->PathArcTo(c, r, -IM_PI * 0.5f, -IM_PI * 0.5f + 2.0f * IM_PI * m.done);
                dl->PathFillConvex(green);
            }
        }
        break;
    }

    // A temporal property goes to a timeline (which a distribution takes as well), a distribution to a distribution,
    // a volume to the density volume
    const char* dnd_type = nullptr;
    const char* hint = "";
    if (m.kinds & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) {
        dnd_type = TIMELINE_SERIES_DND;
        hint = "Drag into a Timelines or Distributions subplot";
    } else if (m.kinds & MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) {
        dnd_type = DISTRIBUTION_SERIES_DND;
        hint = "Drag into a Distributions subplot";
    } else if (m.kinds & MD_SCRIPT_PROPERTY_FLAG_VOLUME) {
        dnd_type = VOLUME_SERIES_DND;
        hint = "Drag into the Density Volume";
    }
    ImGui::InvisibleButton("##mark", ImVec2(d.width, d.height));
    if (dnd_type && ImGui::BeginDragDropSource()) {
        char path[SERIES_PATH_CAP];
        snprintf(path, sizeof(path), "script/%s", m.ident);
        const SeriesKey key = series_key(SeriesSource_Script, str_from_cstr(path));
        char label[96];
        series_label(label, sizeof(label), app, key);
        series_set_drag_payload(dnd_type, key, -1, label, series_default_color(app, key));
        ImGui::EndDragDropSource();
    } else if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("%s%s%s", m.text, hint[0] ? "\n" : "", hint);
    }
}

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
static ViewRegion main_view_uncovered_region(const ApplicationState& state, const ImGuiWindow* canvas);

//static void update_density_volume_texture(ApplicationState* state);

static void render(ApplicationState* state);
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
static void draw_timeline_window(ApplicationState* state);
static void draw_distribution_window(ApplicationState* state);
static void draw_async_task_window(ApplicationState* state);
static void draw_script_editor_window(ApplicationState* state);
static void draw_script_reference_window(ApplicationState* state);
static void open_script_reference(ApplicationState* state, str_t topic, bool take_focus = true);
static void draw_coordinate_system_widget_window(ViewTransform* target, const ViewTransform& current);

static void draw_debug_window(ApplicationState* state);
static void draw_property_export_window(ApplicationState* state);
static void draw_structure_export_window(ApplicationState* state);
static void draw_notifications_window();

static void update_md_buffers(ApplicationState* state);

static bool export_xvg(const float* column_data[], const char* column_labels[], size_t num_columns, size_t num_rows, str_t filename);
static bool export_csv(const float* column_data[], const char* column_labels[], size_t num_columns, size_t num_rows, str_t filename);

static void create_screenshot(str_t path);

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

static void modify_selection(ApplicationState* state, md_bitset_t atom_mask, SelectionOperator op = SelectionOperator::Set) {
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

    // The masks are edited in place. Over any atom index: they hold only the words they write, and a mask over a
    // system's atoms needs no resizing when another is loaded.
    md_bitmask_init(&state.selection.selection_mask,       0, MD_BITSET_MAX_INDEX, persistent_alloc);
    md_bitmask_init(&state.selection.highlight_mask,       0, MD_BITSET_MAX_INDEX, persistent_alloc);
    md_bitmask_init(&state.selection.query.mask,           0, MD_BITSET_MAX_INDEX, persistent_alloc);
    md_bitmask_init(&state.selection.grow.mask,            0, MD_BITSET_MAX_INDEX, persistent_alloc);
    md_bitmask_init(&state.representation.visibility_mask, 0, MD_BITSET_MAX_INDEX, persistent_alloc);
    // Sets, replaced whole
    state.operations.selection_mask = md_bitset_empty_set();
    state.operations.recenter_query.mask = md_bitset_empty_set();

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
    app_settings::bind(STR_LIT("focus_in_free_space"), &state.settings.focus_in_free_space);
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

    viamd::event_system_broadcast_event(viamd::EventType_ViamdInitialize, viamd::EventPayloadType_ApplicationState, &state);

#if EXPERIMENTAL_GFX_API
    md_gfx_initialize(state.gbuffer.width, state.gbuffer.height, 0);
#endif

    ImGui::init_theme();

    state.editor.SetLanguage(script_editor::language());
    state.editor.SetPalette(TextEditor::GetDarkPalette());
    state.editor.SetInsertSpacesOnTabs(true);
    // The marks of how each property stands, left of the line numbers (see draw_script_gutter_mark)
    state.editor.SetDecorationsBeforeLineNumbers(true);
    script_editor::enable_autocomplete(state.editor, &state.mold.sys);

    {
#ifdef VIAMD_DEFAULT_DATASET
        {
            md_str_writer_t sb = md_str_writer_create(frame_alloc);
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

        script_vis_begin_frame(&state);

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

        // GUI
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

        draw_async_task_window(&state);
        draw_main_menu(&state);
        draw_notifications_window();

        //ImGui::ShowDemoWindow();

        draw_coordinate_system_widget_window(&state.view.target, state.view.camera);
        draw_representation_legends(state);
            
        ImGui::BeginCanvas("Main interaction window", true);
        const ImGuiWindow* canvas_window = ImGui::GetCurrentWindow();
        ImVec2 view_size = ImGui::GetContentRegionAvail();
        InteractionSurfaceState surface_state = interaction_surface(interaction_surface_main, vec_cast(view_size));
        ImGui::EndCanvas();

        const mat4_t clip_to_world = camera_view_to_world_matrix(state.view.camera) * state.view.param.matrix.inv.proj;
        const mat4_t world_to_clip = state.view.param.matrix.curr.proj * state.view.param.matrix.curr.view;

        PickingHit hit = {};

        if (surface_state.hovered) {
            InteractionSurfaceHitArgs args = {
                .picking_surface = &state.picking_surface,
                .picking_handler = state.picking_handler,
                .fbo = state.gbuffer.fbo,
                .width = state.gbuffer.width,
                .height = state.gbuffer.height,
                .clip_to_world = clip_to_world,
            };

            interaction_surface_hit_extract(&hit, surface_state, args);

            InteractionSurfaceEvent event = {};
            interaction_surface_event_extract(&event, surface_state, hit);

            event.clip_to_world = clip_to_world;
            event.world_to_clip = world_to_clip;

            if (event.kind == InteractionSurfaceEventKind::RegionSelect) {
                md_bitset_t candidate_mask = md_bitmask_view(&state.representation.visibility_mask);
                if (event.selection_mode == InteractionSelectionMode::Remove) {
                    // When removing, only consider currently selected atoms as candidates for region selection
                    candidate_mask = md_bitmask_view(&state.selection.selection_mask);
                }
                point_set_region_mask_compute(&state.selection.highlight_mask,
                    state.mold.state.xyz, state.mold.state.num_atoms,
                    &candidate_mask, world_to_clip, event.region_min, event.region_max, event.surface_size);

                grow_mask_by_selection_granularity(&state.selection.highlight_mask, state.selection.granularity, state.mold.sys);
                if (event.region_phase == InteractionSurfaceEventPhase::Commit) {
                    // Merge highlight into selection
                    if (event.selection_mode == InteractionSelectionMode::Append) {
                        md_bitmask_or(&state.selection.selection_mask, md_bitmask_view(&state.selection.highlight_mask));
                    }
                    else if (event.selection_mode == InteractionSelectionMode::Remove) {
                        md_bitmask_andnot(&state.selection.selection_mask, md_bitmask_view(&state.selection.highlight_mask));
                    }
                    md_bitmask_clear(&state.selection.highlight_mask);
                    // A region select is not a sequence of individual picks, so whatever order was
                    // recorded before no longer describes the selection.
                    single_selection_sequence_clear(&state.selection.single_selection_sequence);
                }
            } else if (event.kind == InteractionSurfaceEventKind::ContextMenu) {
                ImGui::OpenPopup("Context Popup");
            }
            
            // Since this is the main interaction view, we still broadcast all kinds of events, even if we have handled some explicitly here.
            if (event.kind != InteractionSurfaceEventKind::None) {
                viamd::event_system_broadcast_event(viamd::EventType_ViamdInteractionSurface, viamd::EventPayloadType_InteractionSurfaceEvent, &event);
            }
        }

        InteractionSurfaceViewTransformArgs view_args = {
            .camera = state.view.camera,
            .trackball_param = state.view.trackball_param,
            .view_shift = state.view.shift,
        };

        InteractionSurfaceViewTransformResult view_result = interaction_surface_view_transform_apply(&state.view.target, surface_state, view_args);
        if (view_result.reset_requested) {
            // What is focused on is seen in the middle of what the windows leave in view, rather than behind them
            state.view.region = state.settings.focus_in_free_space ? main_view_uncovered_region(state, canvas_window) : ViewRegion{};
            ViewTransform reset_transform = {};
            if (hit.depth < 1.0f) {
                reset_transform.distance = state.view.target.distance;
                reset_transform.orientation = state.view.camera.orientation;
                reset_transform.position = hit.world_pos + state.view.camera.orientation * vec3_set(0, 0, state.view.target.distance);
            } else {
                reset_view(&reset_transform, state.mold.state, &state.representation.visibility_mask, state.view.camera.fov_y, state.view.region);
            }
            state.view.target = reset_transform;
        }

        draw_context_popup(&state, hit);

        camera_animate(&state.view.camera, state.view.target, state.app.timing.delta_s);
        camera_animate_shift(&state.view.shift, view_region_center(state.view.region), state.app.timing.delta_s);
        state.visuals.dof.focus_depth = state.view.camera.distance;

        ImGuiWindow* win = ImGui::GetCurrentContext()->HoveredWindow;
        if (win && strcmp(win->Name, "Main interaction window") == 0) {
            script_set_hovered_property(&state,  STR_LIT(""));
        }

        // Capture non-window specific keyboard events
        if (!ImGui::GetIO().WantCaptureKeyboard) {
            if (ImGui::IsKeyPressed(KEY_SHOW_DEBUG_WINDOW)) {
                state.show_debug_window = true;
            }

            if (ImGui::IsKeyPressed(KEY_RECOMPILE_SHADERS)) {
                VIAMD_LOG_INFO("Recompiling shaders and re-initializing volume");
                postprocess_pipeline::initialize(state.gbuffer.width, state.gbuffer.height);
                volume::initialize();
                md_gl_shaders_destroy(state.gl.shaders);
                state.gl.shaders = md_gl_shaders_create(shader_output_snippet);
            }

            if (ImGui::IsKeyPressed(KEY_PLAY_PAUSE)) {
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

            if (ImGui::IsKeyPressed(KEY_SKIP_TO_PREV_FRAME) || ImGui::IsKeyPressed(KEY_SKIP_TO_NEXT_FRAME)) {
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

        if (time_changed) {
            state.mold.interpolate_system_state = true;
            time_stopped = false;

            PUSH_CPU_SECTION("Flag dynamic representations for update")
            flag_dynamic_representations_as_dirty(&state);
            POP_CPU_SECTION()
        } else {
            if (!time_stopped) {
                time_stopped = true;
                state.mold.dirty_gpu_buffers |= MolBit_DirtyPosition;
            }
        }

        // The stride only applies while playing: starting or stopping playback, like changing the stride, changes
        // the frames the displayed state is made from without changing the time. Stopped, it is the trajectory's own.
        if (playback_stride(&state) != state.mold.last_interpolated_stride) {
            state.mold.interpolate_system_state = true;
            state.mold.dirty_gpu_buffers |= MolBit_ClearVelocity;
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
                // The evaluations are of what is about to be replaced, so they are stopped. Their IR is not
                // affected: the compilation goes into a new one, and the old is freed once nothing holds it.
                if (state.script.eval) md_script_eval_interrupt(state.script.eval);

                // Try aquire all semaphores
                {
                    defer {
                        flag_all_representations_as_dirty(&state);
                    };

                    // Now we hold all semaphores for the script
                    state.script.compile_ir = false;
                    state.script.time_since_last_change = 0;

                    // The previous IR is freed once nothing uses it (script_ir_collect)
                    state.script.ir = script_ir_create(&state);

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
                    } else if (const FileEntry* run_file = file_list_run_file(&state.files.list)) {
                        // The trajectory's folder, else the structure's
                        extract_folder_path(&cwd, str_from_cstr(run_file->path));
                    }
                    if (!str_empty(cwd)) {
                        md_path_set_cwd(cwd);
                    }

                    const size_t num_stored_selections = md_array_size(state.selection.stored_selections);                       
                    for (size_t i = 0; i < num_stored_selections; ++i) {
                        str_t name = str_from_cstr(state.selection.stored_selections[i].name);
                        const md_bitset_t* bf = &state.selection.stored_selections[i].atom_mask;
                        md_script_ir_add_identifier_bitset(state.script.ir, name, *bf);
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
                            script_editor::markers_add(markers, script_editor::MarkerType_Visualization, tok.depth, tok.range, tok.text, nullptr, tok.ref);
                        }

                        if (md_script_ir_valid(state.script.ir)) {
                            uint64_t ir_figerprint = md_script_ir_fingerprint(state.script.ir);
                            if (state.script.ir_fingerprint != ir_figerprint) {
                                state.script.ir_fingerprint = ir_figerprint;
                            }
                        } else {
                            state.script.ir = nullptr;
                        }
                    }
                }
            }
        }

        // IRs replaced, once no task may be on them. Before the evaluation is started again below, which could
        // otherwise keep a task running every frame
        script_ir_collect(&state);

        if (num_frames > 0) {
            if (state.script.eval_init) {
                if (task_system::task_is_running(state.tasks.evaluate)) {
                    md_script_eval_interrupt(state.script.eval);
                } else {
                    state.script.eval_init = false;

                    // The new evaluation takes over from the previous one: every property computed from the
                    // same thing keeps its values and is not evaluated again (md_script_eval_create_from).
                    // Not remade if the last compilation failed (in the same frame as the request).
                    md_script_eval_t* prev = state.script.eval;
                    state.script.eval = nullptr;

                    if (md_script_ir_valid(state.script.ir)) {
                        state.script.eval_ir = state.script.ir;
                        // The slots over periods come along, and what they hold where their property is taken over
                        const md_script_eval_desc_t desc = { .num_frames = num_frames };
                        state.script.eval = md_script_eval_create_from(prev, state.script.eval_ir, &desc, state.allocator.persistent);
                    }

                    if (prev) {
                        // The versions of the new evaluation carry on from the previous one's, so nothing derived
                        // from those can be taken for these: the caches stay
                        md_script_eval_free(prev);
                    } else {
                        // Its versions start over, in tables possibly allocated where earlier ones were
                        series_cache_free(&state);
                    }

                    state.script.evaluate = true;
                }
            }

            if (state.script.eval && md_script_ir_valid(state.script.eval_ir) &&
                md_script_eval_ir_fingerprint(state.script.eval) == md_script_ir_fingerprint(state.script.eval_ir))
            {
                // The slots over periods follow their periods at once, also while evaluating (series_resolve_histogram).
                // Nothing is evaluated for them until Evaluate is pressed: what needs evaluation only starts when
                // asked to. A distribution needs none, it is kept per frame.
                if (state.script.evaluate) {
                    if (task_system::task_is_running(state.tasks.evaluate)) {
                        // Started again once it has stopped, over the frames which remain then. Nothing it has
                        // done is lost.
                        md_script_eval_interrupt(state.script.eval);
                    } else {
                        state.script.evaluate = false;
                        md_script_eval_reset_interrupt(state.script.eval);

                        // What is left to do, for the properties and the slots
                        const md_bitset_t pending = md_script_eval_pending_frames(state.script.eval, frame_alloc);
                        const size_t count = md_bitset_count(pending);
                        if (count > 0) {
                            md_array_resize(state.script.eval_frames, count, state.allocator.persistent);
                            size_t n = 0;
                            md_bitset_iter_t it = md_bitset_iter(pending);
                            while (md_bitset_iter_next(&it) && n < count) {
                                state.script.eval_frames[n++] = (uint32_t)it.idx;
                            }

                            // The IR is taken now: eval_ir may be replaced while it runs, and the IR it runs on is
                            // not freed until it has finished. Nor are the frames touched.
                            const md_script_ir_t* ir = state.script.eval_ir;
                            md_script_eval_t* eval = state.script.eval;
                            const uint32_t* frames = state.script.eval_frames;
                            state.tasks.evaluate = task_system::create_pool_task(STR_LIT("Evaluate script"), (uint32_t)n, [&state, ir, eval, frames](uint32_t beg, uint32_t end, uint32_t thread_num) {
                                (void)thread_num;
                                md_script_eval_frames(eval, ir, &state.mold.sys, str_from_cstr(state.mold.run), frames + beg, end - beg);
                            });

#if MEASURE_EVALUATION_TIME
                            uint64_t time = (uint64_t)md_tick_now();
                            task_system::ID time_task = task_system::create_pool_task(STR_LIT("##Time Eval"), [t0 = time, n]() {
                                uint64_t t1 = md_tick_now();
                                double s = md_tick_to_seconds(t1 - t0);
                                VIAMD_LOG_INFO("Evaluation of %zu frames completed in: %.3fs", n, s);
                            });
                            task_system::set_task_dependency(time_task, state.tasks.evaluate);
#endif
                            task_system::enqueue_task(state.tasks.evaluate);
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

        // Everything which sets what the script visualizes has run
        script_vis_update(&state);

        if (state.script.vis_shown && state.script.vis_shown->text) {
            PUSH_CPU_SECTION("Draw vis text");
            ImGuiWindow* window = ImGui::FindWindowByName("Main interaction window");
            if (window) {
                ImDrawList* dl = window->DrawList;
                ASSERT(dl);

                const vec2_t res = { (float)state.app.window.width, (float)state.app.window.height };
                const mat4_t mvp = state.view.param.matrix.curr.proj_no_jitter * state.view.param.matrix.curr.view;

                // Script text
                const ImU32 text_color = convert_color(state.script.text_color);
                const ImU32 rect_color = convert_color(state.script.text_bg_color);
                const float rect_rounding = 5.f;
                const ImVec2 rect_padding = ImVec2(4.f, 2.f);

                const md_script_vis_text_t* texts = state.script.vis_shown->text;
                size_t num_text = md_array_size(texts);
                for (size_t i = 0; i < num_text; ++i) {
                    const md_script_vis_text_t& vis_text = texts[i];

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

                    const ImVec2 text_size = ImGui::CalcTextSize(str.beg(), str.end());

                    if (-1 < c.x && c.x < 1 && -1 < c.y && c.y < 1 && -1 < c.z && c.z < 1) {
                        ImVec2 tc = {(c.x * 0.5f + 0.5f) * res.x, (-c.y * 0.5f + 0.5f) * res.y};
                        ImVec2 p0 = tc - text_size * 0.5f;
                        ImVec2 p1 = tc + text_size * 0.5f;
                        dl->AddRectFilled(p0 - rect_padding, p1 + rect_padding, rect_color, rect_rounding);
                        dl->AddText(p0, text_color, str.beg(), str.end());
                    }
                }
            }
            POP_CPU_SECTION();
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
                reset_view(&state.view.target, state.mold.state, (const md_bitset_t*)nullptr, state.view.camera.fov_y, state.view.region);
            }
        }

        // The motivation for doing this is to reduce the frequency at which we invalidate and upload the atom flag field to the GPU
        // For large systems, this can be a costly operation: Consider a system of 100'000'000 atoms
        // The size of each bitfield to represent the mask would be 12.5 MB
        // Given a throughput of 10GB/s results in a time of 1.25 ms per bitfield.

        uint64_t v_hash = state.representation.visibility_mask_hash;
        uint64_t h_hash = md_bitset_hash64(md_bitmask_view(&state.selection.highlight_mask), 0);
        uint64_t s_hash = md_bitset_hash64(md_bitmask_view(&state.selection.selection_mask), 0);
        uint64_t f_hash = md_hash64_combine(v_hash, md_hash64_combine(h_hash, s_hash));
        
        // What each representation is: not the bytes of the struct, which hold handles, pointers and flags that
        // change without the representation changing (its atom mask, for one, is replaced whenever it is evaluated)
        uint64_t r_hash = md_array_size(state.representation.reps);
        for (size_t i = 0; i < md_array_size(state.representation.reps); ++i) {
            const Representation& rep = state.representation.reps[i];
            const uint32_t kind[3] = { (uint32_t)rep.type, (uint32_t)rep.color_mapping, (uint32_t)representation_shown(state, rep) };
            r_hash = md_hash64(kind, sizeof(kind), r_hash);
            const str_t filt = str_from_cstrn(rep.filter.text, sizeof(rep.filter.text));
            r_hash = md_hash64(filt.ptr, filt.len, r_hash);
            r_hash = md_bitset_hash64(rep.atom_mask, r_hash);
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

        update_md_buffers(&state);
        update_timeline_time_unit(&state);

        render(&state);

        task_system::execute_main_task_queue();

        // Reset frame allocator
        md_vm_arena_reset(frame_alloc);

        // Swap buffers
        application::swap_buffers(&state.app);
    }

    interrupt_async_tasks(&state);
    script_vis_reset(&state);
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
    param.matrix.prev = param.matrix.curr;
    param.jitter.prev = param.jitter.curr;

    param.clip_planes.near = data->view.camera.near_plane;
    param.clip_planes.far = data->view.camera.far_plane;
    param.fov_y = data->view.camera.fov_y;
    param.resolution = {(float)data->gbuffer.width, (float)data->gbuffer.height};

    param.matrix.curr.view = camera_world_to_view_matrix(data->view.camera) * data->mold.unitcell_transform;
    param.matrix.inv.view  = mat4_inverse(data->mold.unitcell_transform) * camera_view_to_world_matrix(data->view.camera);

    const float n = data->view.camera.near_plane;
    const float f = data->view.camera.far_plane;
    const float aspect_ratio = (float)data->gbuffer.width / (float)data->gbuffer.height;

    if (data->visuals.temporal_aa.enabled && data->visuals.temporal_aa.jitter) {
        static uint32_t i = 0;
        i = (i+1) % (uint32_t)ARRAY_SIZE(data->view.jitter.sequence);
        param.jitter.curr = data->view.jitter.sequence[i] - 0.5f;
        if (data->view.mode == CameraMode::Perspective) {
            const vec2_t j = param.jitter.curr;
            param.matrix.curr.proj = camera_view_to_clip_matrix_persp(data->view.camera, data->gbuffer.width, data->gbuffer.height, j.x, j.y);
            param.matrix.inv.proj  = camera_clip_to_view_matrix_persp(data->view.camera, data->gbuffer.width, data->gbuffer.height, j.x, j.y);
            param.matrix.curr.proj_no_jitter = camera_view_to_clip_matrix_persp(data->view.camera, aspect_ratio);
        } else {
            const float h = data->view.camera.distance * tanf(data->view.camera.fov_y * 0.5f);
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
            param.matrix.curr.proj = camera_view_to_clip_matrix_persp(data->view.camera, aspect_ratio);
            param.matrix.inv.proj = camera_clip_to_view_matrix_persp(data->view.camera, (float)data->gbuffer.width / (float)data->gbuffer.height);
        } else {
            const float h = data->view.camera.distance * tanf(data->view.camera.fov_y * 0.5f);
            const float w = aspect_ratio * h;
            param.matrix.curr.proj = camera_view_to_clip_matrix_ortho(-w, w, -h, h, n, f);
            param.matrix.inv.proj = camera_clip_to_view_matrix_ortho(-w, w, -h, h, n, f);
        }
        param.matrix.curr.proj_no_jitter = param.matrix.curr.proj;
    }

    // What the view is focused on is seen where the windows leave it uncovered (see ViewRegion): the picture is moved
    // over in clip space. Not in a screenshot without the windows, which shows it in the middle.
    const bool windows_shown = str_empty(data->screenshot.path_to_file) || !data->screenshot.hide_gui;
    const vec2_t shift = windows_shown ? data->view.shift : vec2_t{0, 0};
    if (shift != vec2_t{0, 0}) {
        const mat4_t to_shifted   = mat4_translate( shift.x,  shift.y, 0.0f);
        const mat4_t from_shifted = mat4_translate(-shift.x, -shift.y, 0.0f);
        param.matrix.curr.proj           = to_shifted * param.matrix.curr.proj;
        param.matrix.curr.proj_no_jitter = to_shifted * param.matrix.curr.proj_no_jitter;
        param.matrix.inv.proj            = param.matrix.inv.proj * from_shifted;
    }

    param.matrix.curr.norm = mat4_transpose(param.matrix.inv.view);
}

// The part of the main view that the windows leave in view, as a region of it: the main view is drawn over all of
// the main viewport, behind the windows. All of it when they leave too little to show anything in.
static ViewRegion main_view_uncovered_region(const ApplicationState& state, const ImGuiWindow* canvas) {
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    const ImRect bounds = ImRect(viewport->Pos, viewport->Pos + viewport->Size);
    const ImRect rect = ImGui::LargestUncoveredRect(bounds, canvas);

    ViewRegion region = {.aspect_ratio = (float)state.gbuffer.width / (float)state.gbuffer.height};

    // Smallest side, relative to that of the viewport, that is still worth fitting into
    const float MIN_SIDE = 0.25f;
    if (ImMin(rect.GetWidth(), rect.GetHeight()) >= MIN_SIDE * ImMin(bounds.GetWidth(), bounds.GetHeight())) {
        // Normalized device coordinates have y up
        const ImVec2 p0 = (rect.Min - bounds.Min) / bounds.GetSize();
        const ImVec2 p1 = (rect.Max - bounds.Min) / bounds.GetSize();
        region.min = {p0.x * 2.0f - 1.0f, 1.0f - p1.y * 2.0f};
        region.max = {p1.x * 2.0f - 1.0f, 1.0f - p0.y * 2.0f};
    }
    return region;
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
                reset_view(&data->view.target, data->mold.state, &data->representation.visibility_mask, data->view.camera.fov_y, data->view.region);
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

            const bool basic_shading = data->visuals.shading.mode == postprocess_pipeline::ShadingMode_Basic;

            ImGui::BeginGroup();
            ImGui::Text("Background");
            ImGui::ColorEdit3Minimal("Color", data->visuals.background.color.elem);
            if (!basic_shading) {
                // Basic shading shows the background as the colour itself
                ImGui::SameLine();
                ImGui::SliderFloat("##Intensity", &data->visuals.background.intensity, 0.f, 100.f);
                ImGui::SetItemTooltip("Background brightness");
            }
            ImGui::EndGroup();
            ImGui::Separator();

            // Shading
            ImGui::BeginGroup();
            ImGui::PushID("Shading");
            {
                const char* modes[postprocess_pipeline::ShadingMode_Count] = {"Physically Based", "Basic"};
                int mode = (int)data->visuals.shading.mode;
                if (ImGui::Combo("Shading", &mode, modes, IM_ARRAYSIZE(modes))) {
                    data->visuals.shading.mode = (postprocess_pipeline::ShadingMode)mode;
                }
                ImGui::SetItemTooltip("Physically Based: lit as a plastic-like material by a light and its surroundings, then tone mapped.\n"
                                      "Basic: the colors as they are picked where a surface faces the light, darker away from it. No tone mapping.");
                if (data->visuals.shading.mode == postprocess_pipeline::ShadingMode_Basic) {
                    ImGui::SliderFloat("Ambient", &data->visuals.shading.ambient, 0.0f, 1.0f);
                    ImGui::SetItemTooltip("How much of the color is left on the side facing away from the light");
                    ImGui::SliderFloat("Specular", &data->visuals.shading.specular, 0.0f, 1.0f);
                    ImGui::SetItemTooltip("Strength of the white highlight, 0 for none");
                }
            }
            ImGui::PopID();
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
            ImGui::Checkbox("Depth of Field", &data->visuals.dof.enabled);
            if (data->visuals.dof.enabled) {
                // ImGui::SliderFloat("Focus Point", &data->visuals.dof.focus_depth, 0.001f, 200.f);
                ImGui::SliderFloat("Blur Strength", &data->visuals.dof.aperture, 0.0f, 4.0f, "%.2f %%");
                ImGui::SetItemTooltip("Blur of distant objects, in percent of the view height. Independent of zoom level");
            }
            ImGui::EndGroup();
            ImGui::Separator();

            // Tonemapping
            ImGui::BeginGroup();
            ImGui::BeginDisabled(basic_shading);
            ImGui::Checkbox("Tone Mapping", &data->visuals.tonemapping.enabled);
            ImGui::SetItemTooltip(basic_shading ? "Not used by Basic shading, which shows the colors as they are picked"
                                                : "Map the lighting into the colors a display can show");
            if (data->visuals.tonemapping.enabled && !basic_shading) {
                // Neutral first: it is the default, and the one which keeps the colors
                static const postprocess_pipeline::Tonemapper tonemappers[] = {
                    postprocess_pipeline::Tonemapper_Neutral,
                    postprocess_pipeline::Tonemapper_ACES,
                    postprocess_pipeline::Tonemapper_Filmic,
                    postprocess_pipeline::Tonemapper_ExposureGamma,
                    postprocess_pipeline::Tonemapper_Passthrough,
                };
                static const char* tonemapper_names[] = {"Neutral", "ACES", "Filmic", "Exposure Gamma", "Passthrough"};
                STATIC_ASSERT(ARRAY_SIZE(tonemappers) == postprocess_pipeline::Tonemapper_Count, "every tone mapper is listed");
                int idx = 0;
                for (int i = 0; i < (int)ARRAY_SIZE(tonemappers); ++i) {
                    if (tonemappers[i] == data->visuals.tonemapping.tonemapper) idx = i;
                }
                if (ImGui::Combo("Function", &idx, tonemapper_names, IM_ARRAYSIZE(tonemapper_names))) {
                    data->visuals.tonemapping.tonemapper = tonemappers[idx];
                }
                ImGui::SetItemTooltip("Neutral: the colors as they are, only the highlights are compressed (Khronos PBR Neutral).\n"
                                      "ACES: a filmic curve with more contrast, which shifts and desaturates the colors.\n"
                                      "Filmic: the Uncharted 2 curve.\n"
                                      "Exposure Gamma: exposure and gamma only, no compression of the highlights.\n"
                                      "Passthrough: the values as they are, no gamma.");
                ImGui::SliderFloat("Exposure", &data->visuals.tonemapping.exposure, 0.01f, 10.f);
                ImGui::SliderFloat("Gamma", &data->visuals.tonemapping.gamma, 1.0f, 3.0f);
            }
            ImGui::EndDisabled();
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
            size_t num_selected_atoms = md_bitmask_count(&data->selection.selection_mask);
            if (ImGui::MenuItem("Invert")) {
                md_bitmask_invert_range(&data->selection.selection_mask, 0, data->mold.sys.atom.count);
            }
            if (ImGui::IsItemHovered()) {
                md_bitmask_assign(&data->selection.highlight_mask, md_bitmask_view(&data->selection.selection_mask));
                md_bitmask_invert_range(&data->selection.highlight_mask, 0, data->mold.sys.atom.count);
            }
            if (ImGui::MenuItem("Query")) data->selection.query.show_window = true;
            if (num_selected_atoms == 0) ImGui::PushDisabled();
            if (ImGui::MenuItem("Grow"))  data->selection.grow.show_window = true;
            if (num_selected_atoms == 0) ImGui::PopDisabled();
            if (ImGui::MenuItem("Clear")) {
                md_bitmask_clear(&data->selection.selection_mask);
                single_selection_sequence_clear(&data->selection.single_selection_sequence);
            }
            ImGui::Spacing();
            ImGui::Separator();

            // STORED SELECTIONS
            {
                md_bitmask_clear(&data->selection.highlight_mask);
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
                        md_bitmask_assign(&data->selection.selection_mask, sel.atom_mask);
                        flag_all_representations_as_dirty(data);
                    }
                    if (ImGui::IsItemHovered()) {
                        ImGui::SetTooltip("Load the stored selection as the active selection");
                        md_bitmask_assign(&data->selection.highlight_mask, sel.atom_mask);
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Store")) {
                        ImGui::SetTooltip("Store the active selection into this selection");
                        md_bitset_assign(&sel.atom_mask, md_bitmask_view(&data->selection.selection_mask), persistent_alloc);
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
                    create_selection(data, str_from_cstr(name_buf), md_bitmask_view(&data->selection.selection_mask));
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
            if (ImGui::MenuItem("Take Screenshot")) {
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
                    md_bitmask_clear(&data->selection.highlight_mask);
                }
                ImGui::TableSetupColumn("Now");
                ImGui::TableSetupColumn("Every frame");
                ImGui::TableSetupColumn("Options");
                ImGui::TableHeadersRow();

                const float button_width = (ImGui::GetFontSize() / 20.f) * 175.f;
                const md_bitset_t& target_mask = recenter_get_active_target_mask(data);
                const bool recenter_available = !md_bitset_empty(target_mask);
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
                        md_bitmask_assign(&data->selection.highlight_mask, target_mask);
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
                    if (filter_input("##recenter-query", &recenter_query.filter, filter_context(data))) {
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
                if (!task_system::task_is_running(data->tasks.evaluate) && !task_system::task_is_running(data->script.vis_task)) {
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

            if (ImGui::Checkbox("Focus in Free Space", &data->settings.focus_in_free_space)) {
                if (!data->settings.focus_in_free_space) {
                    // Back to the middle of the whole view
                    data->view.region = {};
                }
                app_settings::mark_dirty();
            }
            ImGui::SetItemTooltip("Double clicking the view frames what it focuses on in the largest part of the view\nthat no window covers, rather than in the middle of the whole view");

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
            if (type == LoaderType_LAMMPSDATA) {
                load_state.arg = state.atom_format_buf;
            }
            // A loader other than the extension would pick has to be stored with the file in a workspace
            str_t ext = {};
            extract_ext(&ext, path);
            if (type != loader::type_from_ext(ext)) {
                load_state.flags |= LoaderFlag_ExplicitType;
            }

            const LoadResult result = load_data_from_file(data, path, load_state);
            if (result == LoadResult_System || result == LoadResult_Trajectory) {
                if (result == LoadResult_System && !data->settings.keep_representations) {
                    remove_all_representations(data);
                    create_default_representations(data);
                }
                data->animation = {};
                reset_view(&data->view.target, data->mold.state, &data->representation.visibility_mask, data->view.camera.fov_y, data->view.region);
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
static void write_script_range(md_str_writer_t& sb, const int* indices, size_t num_indices, int ref_idx = 0) {
    if (num_indices == 0) return;
    if (num_indices == 1) {
        md_str_writer_fmt(&sb, "%i", indices[0] - ref_idx + 1);
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
            md_str_writer_fmt(&sb, "%i", item.beg);
        }
        else {
            md_str_writer_fmt(&sb, "%i:%i", item.beg, item.end);
        }
        if (i < num_items - 1) {
            sb += ',';
        }
    }
    if (num_items > 1) sb += "}";
    md_temp_end(temp);
}

static md_array(str_t) generate_script_selection_suggestions(str_t ident, md_bitset_t sel, const md_system_t* sys) {
    const md_bitset_t* bf = &sel;
    md_array(str_t) suggestions = 0;

    bool within_same_comp = true;
    bool within_same_inst = true;

    md_instance_idx_t inst_idx = -1;
    md_component_idx_t comp_idx = -1;

    md_bitset_iter_t it = md_bitset_iter(*bf);
    while (md_bitset_iter_next(&it)) {
        uint64_t a_idx = it.idx;
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
    
    const size_t popcount = md_bitset_count(*bf);
    
    md_str_writer_t sb = md_str_writer_create(frame_alloc);
    defer { md_str_writer_free(&sb); };

    auto write_atom_remainder = [](md_str_writer_t& sb, const md_bitset_t* bf, int ref_idx = 0) {
        // Add any remainder
        const uint64_t remainder = md_bitset_count(*bf);
        if (remainder) {
            int* indices = (int*)md_alloc(frame_alloc, remainder * sizeof(int));
            defer { md_free(frame_alloc, indices, remainder * sizeof(int)); };

            md_bitset_iter_extract_indices(indices, remainder, md_bitset_iter(*bf));
            sb += "atom(";
            write_script_range(sb, indices, remainder, ref_idx);
            sb += ")";
        }
    };

    if (comp_idx != -1 && within_same_comp) {
        const md_urange_t range = md_component_atom_range(&sys->component, comp_idx);
        if (popcount != range.end - range.beg) {
            md_str_writer_reset(&sb);
            sb += ident;
            sb += " = ";

            // Subset of residue is selected
            write_atom_remainder(sb, bf, range.beg);
            if (md_str_writer_len(&sb) < 512) {
                str_t resname = md_component_name(&sys->component, comp_idx);
                md_str_writer_fmt(&sb, " in resname(\"" STR_FMT "\");", STR_ARG(resname));
                md_array_push(suggestions, str_copy(md_str_writer_to_str(&sb), frame_alloc), frame_alloc);
            }
        }
    }

    else if (inst_idx != -1 && within_same_inst) {
        const md_urange_t range = md_system_instance_atom_range(sys, inst_idx);
        if (popcount != range.end - range.beg) {
            md_str_writer_reset(&sb);
            sb += ident;
            sb += " = ";

            // Subset of chain is selected
            write_atom_remainder(sb, bf, range.beg);
            if (md_str_writer_len(&sb) < 512) {
				str_t inst_id = md_instance_id(&sys->instance, inst_idx);
				md_str_writer_fmt(&sb, " in chain(\"" STR_FMT "\");", STR_ARG(inst_id));
                md_array_push(suggestions, str_copy(md_str_writer_to_str(&sb), frame_alloc), frame_alloc);
            }
        }
    }

    // What is left once complete chains and residues are taken out
    md_bitmask_t remainder = md_bitmask_create(0, (uint32_t)sys->atom.count, frame_alloc);
    md_bitmask_assign(&remainder, *bf);
    auto remainder_has_all = [&remainder](md_urange_t range) {
        return range.beg < range.end && md_bitmask_test_all_range(&remainder, range.beg, range.end);
    };

    md_array(int) complete_chains = 0;
    md_array(int) complete_residues = 0;
    
    if (sys->instance.count) {
        for (size_t i = 0; i < sys->instance.count; ++i) {    
            const md_urange_t range = md_system_instance_atom_range(sys, i);
            if (remainder_has_all(range)) {
                md_array_push(complete_chains, (int)i, frame_alloc);
                md_bitmask_clear_range(&remainder, range.beg, range.end);
            }
        }
    }

    if (sys->component.count) {
        for (size_t i = 0; i < sys->component.count; ++i) {    
            const md_urange_t range = md_component_atom_range(&sys->component, i);
            if (remainder_has_all(range)) {
                md_array_push(complete_residues, (int)i, frame_alloc);
                md_bitmask_clear_range(&remainder, range.beg, range.end);
            }
        }
    }

    const uint64_t atom_remainder_count = md_bitmask_count(&remainder);
    const md_bitset_t tmp_bf = md_bitmask_view(&remainder);
    
    if (complete_chains) {
        md_str_writer_reset(&sb);
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
        
        if (md_str_writer_len(&sb) < 512) {
            md_str_writer_push_char(&sb, ';');
            md_array_push(suggestions, str_copy(md_str_writer_to_str(&sb), frame_alloc), frame_alloc);
        }
        
        md_str_writer_reset(&sb);
        sb += ident;
        sb += " = residue(";
        for (size_t i = 0; i < md_array_size(complete_chains); ++i) {
            md_urange_t range = md_instance_component_range(&sys->instance, complete_chains[i]);
            md_str_writer_fmt(&sb, "%i:%i,", range.beg + 1, range.end);
        }
        if (complete_residues) {
            write_script_range(sb, complete_residues, md_array_size(complete_residues));
        } else {
            md_str_writer_pop(&sb, 1);
        }
        sb += ")";

        if (atom_remainder_count) {
            sb += " or ";
            write_atom_remainder(sb, &tmp_bf);
        }
        
        if (md_str_writer_len(&sb) < 512) {
            md_str_writer_push_char(&sb, ';');
            md_array_push(suggestions, str_copy(md_str_writer_to_str(&sb), frame_alloc), frame_alloc);
        }
    }

    if (complete_residues) {
        md_str_writer_reset(&sb);
        sb += ident;
        sb += " = residue(";
        write_script_range(sb, complete_residues, md_array_size(complete_residues));
        sb += ")";

        if (atom_remainder_count) {
            sb += " or ";
            write_atom_remainder(sb, &tmp_bf);
        }

        if (md_str_writer_len(&sb) < 512) {
            md_str_writer_push_char(&sb, ';');
            md_array_push(suggestions, str_copy(md_str_writer_to_str(&sb), frame_alloc), frame_alloc);
        }
    }

    if (popcount) {
        md_str_writer_reset(&sb);
        sb += ident;
        sb += " = ";
        write_atom_remainder(sb, bf);
        if (md_str_writer_len(&sb) < 512) {
            md_str_writer_push_char(&sb, ';');
            md_array_push(suggestions, str_copy(md_str_writer_to_str(&sb), frame_alloc), frame_alloc);
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
    const size_t num_atoms_selected = md_bitmask_count(&state->selection.selection_mask);

    if (ImGui::BeginPopup("Context Popup")) {
        if (num_atoms_selected == 2) {
            // Suggest construction of covalent bond
            int idx[2];
            md_bitset_iter_extract_indices(idx, 2, md_bitmask_iter(&state->selection.selection_mask));

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
                const md_bitset_t sel = md_bitmask_view(&state->selection.selection_mask);
                str_t ident = create_unique_identifier(state->script.ir, STR_LIT("sel"), frame_alloc);
                
                md_array(str_t) suggestions = generate_script_selection_suggestions(ident, sel, &state->mold.sys);

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
                md_bitmask_invert_range(&state->selection.selection_mask, 0, state->mold.sys.atom.count);
                ImGui::CloseCurrentPopup();
            }
            if (ImGui::IsItemHovered()) {
                md_bitmask_assign(&state->selection.highlight_mask, md_bitmask_view(&state->selection.selection_mask));
                md_bitmask_invert_range(&state->selection.highlight_mask, 0, state->mold.sys.atom.count);
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
                    md_bitmask_clear(&state->selection.selection_mask);
                    single_selection_sequence_clear(&state->selection.single_selection_sequence);
                }
            }
            ImGui::EndMenu();
        }
        
        if (num_atoms_selected > 0) {
            if (ImGui::MenuItem("Set as Centering Target")) {
                md_bitset_assign(&state->operations.selection_mask, md_bitmask_view(&state->selection.selection_mask), persistent_alloc);
                state->operations.recenter_query.enabled = false;
                recenter_update_target_data(state);
                ImGui::CloseCurrentPopup();
            }
            if (ImGui::IsItemHovered()) {
                md_bitmask_assign(&state->selection.highlight_mask, md_bitmask_view(&state->selection.selection_mask));
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
        const uint64_t popcount = md_bitmask_count(&data->selection.selection_mask);
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
            md_bitmask_t* grow = &data->selection.grow.mask;
            md_bitmask_assign(grow, md_bitmask_view(&data->selection.selection_mask));
            const md_bitset_t viable = md_bitmask_view(&data->representation.visibility_mask);

            switch (data->selection.grow.mode) {
            case SelectionGrowthMode::CovalentBond:
                md_util_mask_grow_by_bonds(grow, &data->mold.sys, (int)data->selection.grow.extent, &viable);
                break;
            case SelectionGrowthMode::Radial: {
                md_util_mask_grow_by_radius(grow, &data->mold.state, data->selection.grow.extent, &viable);
                break;
            }
            default:
                ASSERT(false);
            }

            grow_mask_by_selection_granularity(grow, data->selection.granularity, data->mold.sys);
        }

        const bool show_preview =   (ImGui::GetHoveredID() == ImGui::GetID("##Extent")) ||
                                    (ImGui::GetActiveID()  == ImGui::GetID("##Extent")) ||
                                    (ImGui::GetHoveredID() == ImGui::GetID("Apply"));

        if (show_preview) {
            md_bitmask_assign(&data->selection.highlight_mask, md_bitmask_view(&data->selection.grow.mask));
        }
        if (apply) {
            md_bitmask_assign(&data->selection.selection_mask, md_bitmask_view(&data->selection.grow.mask));
            data->selection.grow.mask_invalid = true;
        }
    }
    ImGui::End();
}

static void draw_selection_query_window(ApplicationState* data) {
    ImGui::SetNextWindowSize(ImVec2(300,100), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Selection Query", &data->selection.query.show_window, ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoCollapse)) {

        auto& query = data->selection.query;
        // Escape closes the window, unless it closes the suggestions of the query
        if (ImGui::IsKeyPressed(ImGuiKey_Escape, false) && !filter_completions_open()) {
            query.show_window = false;
        }

        static double query_frame = 0.0;

        ImGui::PushItemWidth(-1);
        bool apply = false;
        if (filter_input("##query", &query.filter, filter_context(data), ImGuiInputTextFlags_AutoSelectAll, &apply)) {
            query.query_invalid = true;
        }
        ImGui::PopItemWidth();

        // What depends on the coordinates is picked again in another frame
        if (query.filter.result.is_dynamic && data->animation.frame != query_frame) {
            query.query_invalid = true;
        }
        bool preview = ImGui::IsItemFocused() || ImGui::IsItemHovered();

        if (ImGui::IsWindowAppearing()) {
            ImGui::SetKeyboardFocusHere(-1);
        }

        if (!query.filter.result.success) ImGui::PushDisabled();
        apply |= ImGui::Button("Apply");
        if (!query.filter.result.success) ImGui::PopDisabled();

        preview |= ImGui::IsItemHovered();

        if (query.query_invalid) {
            query.query_invalid = false;
            md_bitset_t result = md_bitset_empty_set();
            filter_evaluate(&query.filter, &result, filter_context(data), frame_alloc);
            query_frame = data->animation.frame;

            md_bitmask_assign(&query.mask, query.filter.result.success ? result : md_bitset_empty_set());
            grow_mask_by_selection_granularity(&query.mask, data->selection.granularity, data->mold.sys);
        }

        if (preview) {
            md_bitmask_assign(&data->selection.highlight_mask, md_bitmask_view(&query.mask));
        }

        if (apply && query.filter.result.success) {
            md_bitmask_assign(&data->selection.selection_mask, md_bitmask_view(&data->selection.query.mask));
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
            static const char* mode_tooltip[(int)InterpolationMode::Count] = {
                "The nearest frame, as it is",
                "Straight from frame to frame, through each of them",
                "Smooth (C1) from frame to frame, through each of them",
                "Smoother (C2), and damps the jitter from frame to frame,\n"
                "but passes through none of the frames but the first and the last:\n"
                "the coordinates shown are never a frame's own, not even when stopped on one",
            };
            for (int i = 0; i < (int)InterpolationMode::Count; ++i) {
                if (ImGui::Selectable(interpolation_mode_str[i], (int)data->animation.interpolation == i)) {
                    data->animation.interpolation = (InterpolationMode)i;
					data->mold.interpolate_system_state = true;
                    data->mold.dirty_gpu_buffers |= MolBit_ClearVelocity;
                }
                ImGui::SetItemTooltip("%s", mode_tooltip[i]);
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
        {
            // animation.stride, named for what is seen: the frames kept are the keys, the motion between them is interpolated.
            // The slider covers up to 10 (or the last frame, beyond which a larger stride changes nothing), a larger one is
            // typed in (Ctrl+click): only the lower bound is enforced.
            const int max_stride = (int)CLAMP((int64_t)num_frames - 1, (int64_t)2, (int64_t)10);
            const char* fmt = data->animation.stride == 1 ? "every frame" : "every %d frames";
            ImGui::SliderInt("Keyframes", &data->animation.stride, 1, max_stride, fmt);
            data->animation.stride = MAX(1, data->animation.stride);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("While playing, interpolate between every Nth frame only: smoother motion\n"
                                  "at the expense of the frames in between (and the usual interpolation artifacts).\n"
                                  "The speed is unchanged. A stopped animation shows the trajectory's own frames.\n"
                                  "Ctrl+click to type a value beyond the slider's range.");
            }
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
    }
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
                if (id == data->tasks.evaluate) {
                    md_script_eval_interrupt(data->script.eval);
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
static void plot_series_style_popup(ApplicationState* data, PlotSeries& s, int dim, const char* script_ident, md_script_vis_ref_t vis_ref, int* hovered_pop_idx) {
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
                    script_visualize_ref(data, vis_ref, k, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
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
            script_visualize_ref(data, md_script_ir_property_vis_ref(data->script.eval_ir, ident), -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
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
    const md_attributes_t* table = series_table(data, key.source, str_from_cstr(key.path));
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
            if (ImGui::BeginMenu("Periods")) {
                timeline_periods_menu(&data->timeline.periods);
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
                ImGui::EndMenu();
            }
            ImGui::EndMenuBar();
        }

        if (ImGui::IsWindowFocused() && ImGui::IsKeyPressed(KEY_PLAY_PAUSE, false)) {
            data->animation.mode = data->animation.mode == PlaybackMode::Playing ? PlaybackMode::Stopped : PlaybackMode::Playing;
        }

        if (num_x_values > 0) {
            ImPlotInputMap old_map = ImPlot::GetInputMap();

            static bool is_dragging  = false;

            if (!ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
                is_dragging = false;
            }
            if (!ImGui::IsKeyDown(ImGuiKey_LeftCtrl)) {
                is_dragging = false;
            }

            // Shift-drag marks a period and ctrl-drag moves the current frame: neither may pan the plots
            // meanwhile, and nor may a period being dragged in the lane
            if (ImGui::GetIO().KeyShift || ImGui::GetIO().KeyCtrl || timeline_periods_dragging()) {
                ImPlot::GetInputMap().Pan = ImGuiMouseButton_Middle;
            }

            char x_label[64] = "Frame";
            char x_unit_str[32] = "";
            md_unit_t x_unit = data->timeline.time_unit;
            if (!md_unit_is_none(x_unit)) {
                md_unit_print(x_unit_str, sizeof(x_unit_str), x_unit);
                snprintf(x_label, sizeof(x_label), "Time (%s)", x_unit_str);
            }

            // The periods, above the plots and on their time axis
            TimelinePeriodAxis range_axis;
            range_axis.frame_times = x_values;
            range_axis.num_frames  = num_x_values;
            range_axis.view_beg    = &data->timeline.view_range.beg_x;
            range_axis.view_end    = &data->timeline.view_range.end_x;
            range_axis.time_unit   = x_unit_str;
            timeline_periods_lane(&data->timeline.periods, range_axis, &data->animation.frame);

            // After the lane, which ctrl-drag in it may have moved
            double time = frame_to_time(data->animation.frame, *data);

            ImPlot::BeginSubplots("##Temporal", num_subplots, 1, ImVec2(-1,-1));

            const ImPlotFlags plot_flags = ImPlotFlags_NoBoxSelect | ImPlotFlags_NoFrame;
            const ImPlotAxisFlags axis_flags   = ImPlotAxisFlags_NoSideSwitch;
            const ImPlotAxisFlags axis_flags_y = axis_flags | ImPlotAxisFlags_Opposite | ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_RangeFit;

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

                    // The periods shaded below the series, and shift-press marking a new one
                    timeline_periods_plot(&data->timeline.periods, range_axis, data->animation.frame);

                    // The entry, and the member of its population, under the mouse
                    int  hovered_idx      = -1;
                    int  hovered_pop_idx  = -1;
                    char hovered_label[128] = "";

                    bool print_timeline_tooltip = false;

                    if (ImPlot::IsPlotHovered()) {
                        md_bitmask_clear(&data->selection.highlight_mask);
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
                    } else if (ImPlot::IsPlotHovered() && ImGui::IsMouseDown(ImGuiMouseButton_Left) && ImGui::IsKeyDown(ImGuiKey_LeftCtrl)) {
                        is_dragging = true;
                    }

                    for (int j = 0; j < sp.count; ++j) {
                        PlotSeries& s = sp.series[j];
                        const SeriesTemporalView& v = views[j];

                        if (ImPlot::IsLegendEntryHovered(v.plot_id)) {
                            if (v.script_ident[0] != '\0') {
                                script_visualize_ref(data, v.vis_ref, -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                                script_set_hovered_property(data, str_from_cstr(v.script_ident));
                            }
                            if (!resolved[j]) {
                                ImGui::SetTooltip("Nothing to plot at '%s' right now", s.key.path);
                            } else if (v.state.status != SeriesStatus_Current) {
                                char text[160];
                                series_state_describe(text, sizeof(text), v.state);
                                ImGui::SetTooltip("%s", text);
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
                            plot_series_style_popup(data, s, resolved[j] ? v.dim : 1, v.script_ident, v.vis_ref, &popup_hovered_pop);
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
                                if (v.state.status >= SeriesStatus_Outdated) {
                                    // Not what the script says any more: shown, faintly, until it is evaluated
                                    color.w *= 0.35f;
                                }
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

                    if (ImPlot::IsPlotHovered() && hovered_idx != -1 && views[hovered_idx].script_ident[0] != '\0') {
                        const SeriesTemporalView& v = views[hovered_idx];
                        const int pop_idx = v.dim > 1 ? hovered_pop_idx : -1;
                        script_visualize_ref(data, v.vis_ref, pop_idx, ~MD_SCRIPT_VISUALIZE_SDF);
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

                    // A period being marked or dragged has a tooltip of its own
                    if (print_timeline_tooltip && !timeline_periods_dragging()) {
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

            time = CLAMP(time, (double)min_x_value, (double)max_x_value);
            data->animation.frame = time_to_frame(time, data->timeline.x_values);

            ImPlot::GetInputMap() = old_map;
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

// A polyline drawn dashed, carrying the dashes on across its corners and from one call to the next through phase
static void dashed_polyline(ImDrawList* dl, const ImVec2* pts, int count, ImU32 col, float* phase) {
    const float dash = 6.0f, gap = 4.0f;
    for (int i = 1; i < count; ++i) {
        const ImVec2 prev = pts[i - 1];
        const ImVec2 d = pts[i] - prev;
        const float len = sqrtf(d.x * d.x + d.y * d.y);
        if (len <= 0.0f) continue;
        const ImVec2 u = d / len;
        float t = 0.0f;
        while (t < len) {
            const float in_cycle = fmodf(*phase + t, dash + gap);
            const float step = in_cycle < dash ? MIN(dash - in_cycle, len - t) : MIN(dash + gap - in_cycle, len - t);
            if (in_cycle < dash) dl->AddLine(prev + u * t, prev + u * (t + step), col, 1.5f);
            t += MAX(step, 0.01f);
        }
        *phase = fmodf(*phase + len, dash + gap);
    }
}

// A histogram member over the plot as a preview, not a curve of the plot's own: drawn as the series is (a line, an area
// or bars), with a dashed outline and, for an area or bars, a faint fill
static void plot_histogram_preview(const SeriesHistogramView& v, int k, ImVec4 color, const PlotSeries& s) {
    if (v.num_bins <= 0) return;
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    ImPlot::PushPlotClipRect();
    const ImU32 line_col = ImGui::ColorConvertFloat4ToU32(color);
    const ImU32 fill_col = ImGui::ColorConvertFloat4ToU32(ImVec4(color.x, color.y, color.z, color.w * 0.25f));
    // Fills which meet edge to edge, without the seams anti-aliasing would leave between them
    const ImDrawListFlags flags = dl->Flags;
    float phase = 0.0f;

    if (s.plot_type == PlotType_Bars) {
        const double half = series_histogram_bin_width(v) * s.bar_width * 0.5;
        const bool joined = s.bar_width >= 0.999f;
        md_temp_scope_t temp = md_temp_begin();
        // Joined, the outline is one step along the tops of the bars; apart, each bar has its own
        ImVec2* outline = md_temp_alloc_array(temp, ImVec2, (size_t)v.num_bins * 2 + 2);
        int n = 0;
        dl->Flags &= ~ImDrawListFlags_AntiAliasedFill;
        for (int i = 0; i < v.num_bins; ++i) {
            const ImPlotPoint p = series_histogram_point(v, i, k);
            const ImVec2 a = ImPlot::PlotToPixels(p.x - half, 0.0);
            const ImVec2 b = ImPlot::PlotToPixels(p.x + half, p.y);
            dl->AddRectFilled(ImVec2(MIN(a.x, b.x), MIN(a.y, b.y)), ImVec2(MAX(a.x, b.x), MAX(a.y, b.y)), fill_col);
            if (joined) {
                if (n == 0) outline[n++] = a;
                outline[n++] = ImVec2(a.x, b.y);
                outline[n++] = b;
            } else {
                const ImVec2 bar[4] = { a, ImVec2(a.x, b.y), b, ImVec2(b.x, a.y) };
                dashed_polyline(dl, bar, 4, line_col, &phase);
            }
        }
        dl->Flags = flags;
        if (joined && n > 0) {
            outline[n] = ImVec2(outline[n - 1].x, outline[0].y);
            dashed_polyline(dl, outline, n + 1, line_col, &phase);
        }
        md_temp_end(temp);
    } else {
        md_temp_scope_t temp = md_temp_begin();
        ImVec2* pts = md_temp_alloc_array(temp, ImVec2, (size_t)v.num_bins);
        for (int i = 0; i < v.num_bins; ++i) {
            pts[i] = ImPlot::PlotToPixels(series_histogram_point(v, i, k));
        }
        if (s.plot_type == PlotType_Area) {
            // Down to zero, one quad per stretch between two bins
            dl->Flags &= ~ImDrawListFlags_AntiAliasedFill;
            for (int i = 1; i < v.num_bins; ++i) {
                const ImVec2 base0 = ImPlot::PlotToPixels(series_histogram_point(v, i - 1, k).x, 0.0);
                const ImVec2 base1 = ImPlot::PlotToPixels(series_histogram_point(v, i, k).x, 0.0);
                dl->AddQuadFilled(pts[i - 1], pts[i], base1, base0, fill_col);
            }
            dl->Flags = flags;
        }
        dashed_polyline(dl, pts, v.num_bins, line_col, &phase);
        md_temp_end(temp);
    }
    ImPlot::PopPlotClipRect();
}

// Over periods: every period's value at the bin, beside its frames, the hovered one first among equals
static void plot_periods_tooltip(const SeriesHistogramView& v, const PlotSeries& s, int bin, int hovered) {
    ImGui::BeginTooltip();
    const ImPlotPoint at = series_histogram_point(v, bin, 0);
    ImGui::Text("%s at %.4g %s", v.label, at.x, v.x_unit_str);
    if (ImGui::BeginTable("##periods", 3, ImGuiTableFlags_SizingFixedFit)) {
        for (int k = 0; k < v.num_periods; ++k) {
            const SeriesPeriodInfo& p = v.periods[k];
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            const ImVec4 color = p.own_color ? p.color : s.color;
            ImGui::ColorButton("##c", ImVec4(color.x, color.y, color.z, 1.0f), ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoBorder,
                ImVec2(ImGui::GetTextLineHeight(), ImGui::GetTextLineHeight()));
            ImGui::SameLine();
            if (k == hovered) ImGui::TextUnformatted(p.name);
            else              ImGui::TextDisabled("%s", p.name);
            ImGui::TableNextColumn();
            if (p.valid) ImGui::Text("%.4g", series_histogram_point(v, bin, k).y);
            else         ImGui::TextDisabled("-");
            ImGui::TableNextColumn();
            if (p.num_frames < p.num_wanted) ImGui::TextDisabled("%u of %u frames", p.num_frames, p.num_wanted);
            else                             ImGui::TextDisabled("%u frames", p.num_frames);
        }
        ImGui::EndTable();
    }
    ImGui::EndTooltip();
}

// The part of a script series' legend popup which says over which periods it is shown, one curve each
static void plot_periods_popup(ApplicationState* data, PlotSeries& s, bool* refit) {
    ImGui::SeparatorText("Periods");
    const TimelinePeriods& rs = data->timeline.periods;
    auto item = [&](uint32_t id, const char* name, const ImVec4* color, const char* frames) {
        bool on = plot_series_has_period(s, id);
        const bool last = on && s.num_periods == 1;
        const bool full = !on && s.num_periods >= PLOT_MAX_PERIODS;
        ImGui::PushID((int)id);
        if (last || full) ImGui::BeginDisabled();
        if (ImGui::Checkbox("##on", &on)) {
            if (on) plot_series_add_period(s, id);
            else    plot_series_remove_period(s, id);
            *refit = true;
        }
        if (last || full) ImGui::EndDisabled();
        ImGui::SameLine();
        if (color) {
            ImGui::ColorButton("##c", *color, ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoBorder,
                ImVec2(ImGui::GetTextLineHeight(), ImGui::GetTextLineHeight()));
            ImGui::SameLine();
        }
        ImGui::TextUnformatted(name);
        if (frames[0]) {
            ImGui::SameLine();
            ImGui::TextDisabled("%s", frames);
        }
        const uint32_t period = id != 0 ? id : s.key.period;
        if (period != 0 && ImGui::IsItemHovered()) timeline_periods_highlight(period);
        ImGui::PopID();
    };

    // The key's own frames first, in the series' colour: every frame, or its period's
    const int num_frames = (int)run_num_frames(data);
    auto frames_of = [&](char* buf, size_t cap, const TimelinePeriod& r) {
        int beg, end;
        if (timeline_period_frames(&beg, &end, r, data->animation.frame, num_frames)) snprintf(buf, cap, "%d to %d", beg, end);
    };
    {
        char frames[48] = "";
        const TimelinePeriod* own = timeline_periods_find(&rs, s.key.period);
        if (own) frames_of(frames, sizeof(frames), *own);
        item(0, s.key.period != 0 ? series_period_name(data, s.key.period) : "All frames", nullptr, frames);
    }
    for (int i = 0; i < rs.count; ++i) {
        const TimelinePeriod& r = rs.items[i];
        if (r.id == s.key.period) continue;
        char frames[48] = "";
        frames_of(frames, sizeof(frames), r);
        const ImVec4 color(r.color.x, r.color.y, r.color.z, 1.0f);
        item(r.id, r.name, &color, frames);
    }
    if (rs.count == 0) {
        ImGui::TextDisabled("Mark periods in the timeline to compare it over them");
    } else {
        ImGui::TextDisabled("Or pull a period out of the timeline onto the plot");
    }
}

static void draw_distribution_properties_menu(ApplicationState* data) {
    PlotSubplot& first = data->distributions.subplots[0];
    const uint32_t kinds = MD_SCRIPT_PROPERTY_FLAG_TEMPORAL | MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION;

    // A temporal property is binned over its frames, a distribution is shown as evaluated. Over the
    // frames of a timeline period too: "r (Production)".
    auto script_section = [&](uint32_t period) -> int {
        int num_listed = 0;
        series_for_each_script_property(data, SeriesSource_Script, kinds, [&](const SeriesKey& over_all) {
            SeriesKey key = over_all;
            key.period = period;
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
    if (script_section(0) == 0) {
        ImGui::TextDisabled("No temporal or distribution properties, define and evaluate them in the script editor");
    } else if (data->timeline.periods.count > 0) {
        ImGui::SeparatorText("Script, over a period");
        const TimelinePeriods& rs = data->timeline.periods;
        for (int i = 0; i < rs.count; ++i) {
            const TimelinePeriod& r = rs.items[i];
            char label[TIMELINE_PERIOD_NAME_CAP + 16];
            snprintf(label, sizeof(label), "%s###period%u", r.name, r.id);
            const bool open = ImGui::BeginMenu(label);
            if (ImGui::IsItemHovered()) timeline_periods_highlight(r.id);
            if (open) {
                script_section(r.id);
                ImGui::EndMenu();
            }
        }
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
                    PlotSeries& s = sp.series[j];
                    // A period removed takes its curve along
                    for (int p = s.num_periods - 1; p >= 0; --p) {
                        if (s.periods[p] != 0 && !timeline_periods_find(&data->timeline.periods, s.periods[p])) {
                            if (s.num_periods > 1) plot_series_remove_period(s, s.periods[p]);
                            else s.periods[0] = 0;
                        }
                    }
                    resolved[j] = series_resolve_histogram_periods(&views[j], data, s.key, s.num_bins, s.periods, s.num_periods);
                }

                // More to show than before (a workspace loaded before anything was evaluated, an evaluation which
                // came through): fitted to it
                static int num_resolved[PLOT_MAX_SUBPLOTS] = {};
                int resolved_now = 0;
                for (int j = 0; j < sp.count; ++j) resolved_now += resolved[j] ? 1 : 0;
                if (resolved_now > num_resolved[i]) fit_subplot[i] = true;
                num_resolved[i] = resolved_now;

                // A period under the pointer in the timeline: each series over it too, as a preview
                const uint32_t preview_id = timeline_periods_hovered();
                const TimelinePeriod* preview_range = timeline_periods_find(&data->timeline.periods, preview_id);

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
                    int  hovered_bin       = -1;
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
                                if (v.num_periods > 0 ? !v.periods[k].valid : (dim > 1 && !s.population_mask.test(k))) continue;
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
                                    hovered_pop_idx = (dim > 1 || v.num_periods > 0) ? k : -1;
                                    hovered_bin = xn;
                                    if (v.num_periods > 0) {
                                        snprintf(hovered_label, sizeof(hovered_label), "%s, %s: %.4g at %.4g %s", v.label, v.periods[k].name, pt.y, pt.x, v.x_unit_str);
                                    } else if (dim > 1) {
                                        snprintf(hovered_label, sizeof(hovered_label), "%s[%i]: %.4g at %.4g %s", v.label, k + 1, pt.y, pt.x, v.x_unit_str);
                                    } else {
                                        snprintf(hovered_label, sizeof(hovered_label), "%s: %.4g at %.4g %s", v.label, pt.y, pt.x, v.x_unit_str);
                                    }
                                }
                            }
                        }

                        if (hovered_idx != -1) {
                            const SeriesHistogramView& v = views[hovered_idx];
                            // A period is not a member of a population: the property as a whole is what is shown
                            const int member = v.num_periods > 0 ? -1 : hovered_pop_idx;
                            if (v.script_ident[0] != '\0') {
                                script_set_hovered_property(data, str_from_cstr(v.script_ident), member);
                                script_visualize_ref(data, v.vis_ref, member, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                            }
                            if (v.num_periods > 0 && hovered_bin >= 0) {
                                if (hovered_pop_idx >= 0 && v.periods[hovered_pop_idx].id != 0) {
                                    timeline_periods_highlight(v.periods[hovered_pop_idx].id);
                                }
                                plot_periods_tooltip(v, sp.series[hovered_idx], hovered_bin, hovered_pop_idx);
                            } else if (hovered_label[0] != '\0') {
                                ImGui::SetTooltip("%s", hovered_label);
                            }
                        }
                    } else if (data->hovered_property_label[0] != '\0') {
                        for (int j = 0; j < sp.count; ++j) {
                            if (!resolved[j]) continue;
                            if (strcmp(data->hovered_property_label, views[j].script_ident) == 0) {
                                hovered_idx = j;
                                // A member of a population elsewhere is not a period here
                                hovered_pop_idx = views[j].num_periods > 0 ? -1 : data->hovered_property_pop_idx;
                                break;
                            }
                        }
                    }

                    for (int j = 0; j < sp.count; ++j) {
                        PlotSeries& s = sp.series[j];
                        const SeriesHistogramView& v = views[j];

                        if (ImPlot::IsLegendEntryHovered(v.plot_id)) {
                            if (v.script_ident[0] != '\0') {
                                script_visualize_ref(data, v.vis_ref, -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                                script_set_hovered_property(data, str_from_cstr(v.script_ident));
                            }
                            if (!resolved[j]) {
                                ImGui::SetTooltip("Nothing to show for '%s' right now", s.key.path);
                            } else if (v.state.status != SeriesStatus_Current) {
                                char text[160];
                                series_state_describe(text, sizeof(text), v.state);
                                ImGui::SetTooltip("%s", text);
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
                            plot_series_style_popup(data, s, resolved[j] && v.num_periods == 0 ? v.dim : 1, v.script_ident, v.vis_ref, &popup_hovered_pop);
                            if (popup_hovered_pop != -1) {
                                hovered_idx = j;
                                hovered_pop_idx = popup_hovered_pop;
                            }
                            if (v.script_ident[0] != '\0') {
                                plot_periods_popup(data, s, &fit_subplot[i]);
                            }
                            ImPlot::EndLegendPopup();
                        }

                        if (!resolved[j] || v.num_bins <= 0) {
                            ImPlot::PlotDummy(v.plot_id);
                        } else {
                            const int population_size = CLAMP(v.dim, 1, PLOT_MAX_POPULATION);
                            auto plot = [&](int k) {
                                // A period in its own colour, so that it reads with the timeline
                                ImVec4 color = v.num_periods > 0 && v.periods[k].own_color ? v.periods[k].color : plot_series_member_color(s, k, population_size);
                                if (v.state.status >= SeriesStatus_Outdated) {
                                    // Not what the script says any more: shown, faintly, until it is evaluated
                                    color.w *= 0.35f;
                                }
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
                                if (v.num_periods > 0 ? !v.periods[k].valid : (population_size > 1 && !s.population_mask.test(k))) continue;
                                if (hovered_idx == j && hovered_pop_idx == k) continue;
                                plot(k);
                            }
                            if (hovered_idx == j && hovered_pop_idx != -1 && hovered_pop_idx < population_size) {
                                plot(hovered_pop_idx);
                            }
                        }

                        // The period under the pointer in the timeline, drawn as the series is but dashed, unless the series is
                        // shown over it already
                        if (preview_range && v.script_ident[0] != '\0' && !plot_series_has_period(s, preview_id)) {
                            SeriesHistogramView pv;
                            if (series_resolve_histogram_periods(&pv, data, s.key, s.num_bins, &preview_id, 1) && pv.num_bins > 0 && pv.periods[0].valid) {
                                const vec4_t rc = preview_range->color;
                                plot_histogram_preview(pv, 0, ImVec4(rc.x, rc.y, rc.z, 1.0f), s);
                            }
                        }

                        if (ImPlot::BeginDragDropSourceItem(v.plot_id)) {
                            series_set_drag_payload(DISTRIBUTION_SERIES_DND, s.key, i, v.label, s.color);
                            ImPlot::EndDragDropSource();
                        }
                    }

                    // A period pulled out of the timeline, dropped on a legend entry: that series over it too
                    if (ImPlot::BeginDragDropTargetLegend()) {
                        if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(TIMELINE_PERIOD_DND)) {
                            ASSERT(payload->DataSize == sizeof(uint32_t));
                            const uint32_t id = *(const uint32_t*)payload->Data;
                            for (int j = 0; j < sp.count; ++j) {
                                if (ImPlot::IsLegendEntryHovered(views[j].plot_id) && views[j].script_ident[0] != '\0') {
                                    plot_series_add_period(sp.series[j], id);
                                }
                            }
                        }
                        ImPlot::EndDragDropTarget();
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
                        } else if (const ImGuiPayload* range_payload = ImGui::AcceptDragDropPayload(TIMELINE_PERIOD_DND)) {
                            // A period pulled out of the timeline: every script series here over it too
                            ASSERT(range_payload->DataSize == sizeof(uint32_t));
                            const uint32_t id = *(const uint32_t*)range_payload->Data;
                            for (int j = 0; j < sp.count; ++j) {
                                if (views[j].script_ident[0] != '\0') {
                                    plot_series_add_period(sp.series[j], id);
                                }
                            }
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

        // What Evaluate would do. It says so, and does nothing when there is nothing to do.
        const ScriptPlan plan = script_plan(*state, frame_alloc);
        char btn_text[96];
        bool can_eval = false;
        if (!plan.valid) {
            snprintf(btn_text, sizeof(btn_text), "Evaluate###evaluate");
        } else if (plan.running && plan.same) {
            snprintf(btn_text, sizeof(btn_text), "Evaluating, %u frames left###evaluate", plan.num_pending);
        } else if (plan.num_pending > 0) {
            snprintf(btn_text, sizeof(btn_text), "Evaluate %u frames###evaluate", plan.num_pending);
            can_eval = true;
        } else if (plan.changes) {
            snprintf(btn_text, sizeof(btn_text), "Evaluate###evaluate");
            can_eval = true;
        } else {
            snprintf(btn_text, sizeof(btn_text), "Up to date###evaluate");
        }

        // How each property stands, beside the line it is named on
        md_array(ScriptGutterMark) marks = nullptr;
        if (plan.valid) {
            const uint32_t* line_offsets = state->editor_markers.line_offsets;
            const size_t num_lines = md_array_size(line_offsets);
            for (size_t i = 0; i < plan.num_props; ++i) {
                const md_script_range_marker_t range = md_script_ir_property_range(state->script.ir, plan.names[i]);
                if (range.end <= range.beg || num_lines == 0) continue;
                // The last line which begins at or before it
                size_t lo = 0, hi = num_lines;
                while (hi - lo > 1) {
                    const size_t mid = (lo + hi) / 2;
                    if (line_offsets[mid] <= (uint32_t)range.beg) lo = mid; else hi = mid;
                }
                const md_script_eval_plan_t& p = plan.props[i];
                ScriptGutterMark mark = { (int)lo, p.status, 0.0f, "", "", 0 };
                str_copy_to_char_buf(mark.ident, sizeof(mark.ident), plan.names[i]);
                mark.kinds = (uint32_t)md_script_ir_property_flags(state->script.ir, plan.names[i]);
                const uint32_t needed = p.num_frames + p.num_pending;
                mark.done = p.num_pending == 0 ? 1.0f : needed ? (float)p.num_frames / (float)needed : 0.0f;
                script_plan_describe(mark.text, sizeof(mark.text), plan, i);
                md_array_push(marks, mark, frame_alloc);
            }
        }
        const size_t num_marks = md_array_size(marks);
        state->editor.SetLineDecorator(1, [state, marks, num_marks](TextEditor::Decorator& d) {
            for (size_t i = 0; i < num_marks; ++i) {
                if ((size_t)marks[i].line == d.line) {
                    draw_script_gutter_mark(state, d, marks[i]);
                    break;
                }
            }
        });

        const ImVec2 content_size = ImGui::GetContentRegionAvail();
        const ImVec2 label_size = ImGui::CalcTextSize(btn_text, NULL, true) * ImVec2(1.2, 1.0);
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

        if (!can_eval) ImGui::PushDisabled();
        if (ImGui::Button(btn_text, btn_size)) {
            eval = true;
        }
        if (!can_eval) ImGui::PopDisabled();
        if (plan.valid && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled | ImGuiHoveredFlags_DelayShort)) {
            ImGui::BeginTooltip();
            char line[160];
            for (size_t i = 0; i < plan.num_props; ++i) {
                script_plan_describe(line, sizeof(line), plan, i);
                ImGui::TextUnformatted(line);
            }
            // What is shown which the script no longer has
            if (state->script.eval_ir && !plan.same) {
                const str_t* evaluated = md_script_ir_property_names(state->script.eval_ir);
                for (size_t j = 0; j < md_script_ir_property_count(state->script.eval_ir); ++j) {
                    bool kept = false;
                    for (size_t i = 0; i < plan.num_props && !kept; ++i) {
                        kept = str_eq(plan.names[i], evaluated[j]) || (plan.props[i].status == MD_SCRIPT_EVAL_STATUS_HELD && str_eq(plan.props[i].from, evaluated[j]));
                    }
                    if (!kept) ImGui::TextDisabled(STR_FMT ": no longer in the script", STR_ARG(evaluated[j]));
                }
            }
            ImGui::EndTooltip();
        }

        if (eval && can_eval) {
            // The script as it was evaluated: what remains of it. Otherwise a new evaluation, which takes over what
            // it can from this one.
            if (plan.same) {
                state->script.evaluate = true;
            } else {
                state->script.eval_init = true;
            }
        }

        if (editor_hovered) {
            md_bitmask_clear(&state->selection.highlight_mask);
        }

        if (hovered_marker) {
            if (hovered_marker->atoms || !md_script_vis_ref_empty(hovered_marker->vis_ref)) {
                if (hovered_marker->atoms) {
                    // Errors and warnings: the atoms the message is about
                    md_bitmask_assign(&state->selection.highlight_mask, *hovered_marker->atoms);
                }
                else if (hovered_marker->type == script_editor::MarkerType_Visualization) {
                    // Clear hovered property
                    script_set_hovered_property(state, STR_LIT(""));
                    // The markers come from the IR last compiled. One which did not compile is gone, and its
                    // markers resolve in nothing
                    const md_script_vis_ref_t ref = hovered_marker->vis_ref;
                    if (const md_script_ir_t* ir = script_ir_of(state, ref)) {
                        str_t ref_ident = md_script_vis_ref_ident(ir, ref);
                        int ref_dim  = md_script_vis_ref_dim(ir, ref);

                        if (ref_dim > 1) {
                            int delta = (int)ImGui::GetIO().MouseWheel;
                            if (ImGui::IsKeyDown(ImGuiMod_Shift)) {
                                delta *= 10;
                            }
                            state->script.sub_idx += delta;
                            state->script.sub_idx = CLAMP(state->script.sub_idx, -1, (int)ref_dim - 1);
                        }

                        script_visualize_ref(state, ref, state->script.sub_idx, 0);
                        script_set_hovered_property(state, ref_ident, state->script.sub_idx);
                    }
                }

                bool lm_click = ImGui::IsMouseClicked(ImGuiMouseButton_Left);
                bool rm_click = ImGui::IsMouseClicked(ImGuiMouseButton_Right);
                if (ImGui::IsKeyDown(ImGuiMod_Shift) && (lm_click || rm_click)) {
                    SelectionOperator op = lm_click ? SelectionOperator::Or : SelectionOperator::AndNot;
                    modify_selection(state, md_bitmask_view(&state->selection.highlight_mask), op);
                }
            }

            // After the visualization is set: the value is computed with it
            if (!str_empty(hovered_marker->text)) {
                ImGui::BeginTooltip();
                ImGui::Text(STR_FMT, STR_ARG(hovered_marker->text));
                if (hovered_marker->type == script_editor::MarkerType_Visualization) {
                    const str_t value = script_vis_value_text(state);
                    // A value known since the compilation is in the text already, unless it was too long to keep there
                    const bool dynamic = str_find_str(NULL, hovered_marker->text, STR_LIT("[d]"));
                    const bool has_value = str_find_char(NULL, hovered_marker->text, '\n');
                    if (!str_empty(value) && (dynamic || !has_value || state->script.sub_idx != -1)) {
                        ImGui::Separator();
                        if (dynamic) ImGui::TextDisabled("At the frame shown");
                        ImGui::Text(STR_FMT, STR_ARG(value));
                    }
                }
                if (state->script.sub_idx != -1) {
                    ImGui::TextDisabled("Element %i, the mouse wheel steps through them", state->script.sub_idx + 1);
                }
                ImGui::EndTooltip();
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

static bool export_cube(const ApplicationState& data, const SeriesVolumeView& volume, str_t filename) {
    // @NOTE: First we need to extract some meta data for the cube format, we need the atom indices/bits for any SDF
    // And the origin + extent of the volume in spatial coordinates (Ångström)

    const float* volume_values = volume.data;
    const md_script_vis_ref_t vis_ref = volume.vis_ref;
    if (!volume_values) {
        VIAMD_LOG_ERROR("Export Cube: The property to be exported did not exist");
        return false;
    }

    if (md_script_vis_ref_empty(vis_ref)) {
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
        result = md_script_vis_eval_ref(&vis, vis_ref, 0, &ctx, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_SDF);
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
            const md_bitset_t* bf = &vis.sdf.structures[0];
            const int num_atoms = (int)md_bitset_count(*bf);
            const int vol_dim[3] = {(int)volume.dim[0], (int)volume.dim[1], (int)volume.dim[2]};
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

            md_bitset_iter_t it = md_bitset_iter(*bf);
            while (md_bitset_iter_next(&it)) {
                size_t i = it.idx;
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

        if (task_system::task_is_running(data->tasks.evaluate)) {
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
            auto push_script = [&](uint32_t period) {
                series_for_each_script_property(data, SeriesSource_Script, script_kinds, [&](const SeriesKey& over_all) {
                    SeriesKey key = over_all;
                    key.period = period;
                    md_array_push(candidates, key, alloc);
                    const bool temporal = md_script_ir_property_flags(data->script.eval_ir, series_script_ident(key)) & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL;
                    if (kind == ExportKind_Distribution && temporal && script_property_population(data, key) > 1) {
                        SeriesKey agg = key;
                        agg.variant = SeriesVariant_Aggregate;
                        md_array_push(candidates, agg, alloc);
                    }
                });
            };
            // Over each period too, "r (Production)": what is over time is the same over a period
            push_script(0);
            if (kind != ExportKind_Temporal) {
                for (int i = 0; i < data->timeline.periods.count; ++i) {
                    push_script(data->timeline.periods.items[i].id);
                }
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
                } else if (str_eq(file_extension, STR_LIT("cube")) && export_cube(*data, view, path)) {
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

md_array(int) extract_bitfield_indices(const md_bitset_t* bf, md_allocator_i* alloc) {
    size_t count = md_bitset_count(*bf);
    md_array(int) indices = md_array_create(int, count, alloc);
    md_bitset_iter_extract_indices(indices, count, md_bitset_iter(*bf));
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
            query.requires_evaluation |= filter_input("Query", &query.filter, filter_context(data));
            bool preview = ImGui::IsItemFocused() || ImGui::IsItemHovered();

            if (query.requires_evaluation) {
                filter_evaluate(&query.filter, &query.mask, filter_context(data), persistent_alloc);
                query.requires_evaluation = false;
            }

            if (ImGui::IsWindowHovered()) {
                // Clear highlight mask
                md_bitmask_clear(&data->selection.highlight_mask);
            }

            if (query.filter.result.success && preview) {
                md_bitmask_assign(&data->selection.highlight_mask, query.mask);
            }
        }

        // Every frame, the current one, or those of a timeline period
        const TimelinePeriod* traj_period = timeline_periods_find(&data->timeline.periods, struct_exp.selected_traj_period);
        struct_exp.selected_traj_filter = CLAMP(struct_exp.selected_traj_filter, 0, 2);
        if (struct_exp.selected_traj_filter == 2 && !traj_period) struct_exp.selected_traj_filter = 0;
        if (traj) {
            const char* preview = struct_exp.selected_traj_filter == 2 ? traj_period->name : frame_mask_options[struct_exp.selected_traj_filter];
            if (ImGui::BeginCombo("Frames to Export", preview)) {
                for (int i = 0; i < (int)ARRAY_SIZE(frame_mask_options); ++i) {
                    if (ImGui::Selectable(frame_mask_options[i], struct_exp.selected_traj_filter == i)) struct_exp.selected_traj_filter = i;
                }
                const TimelinePeriods& rs = data->timeline.periods;
                if (rs.count > 0) ImGui::SeparatorText("Period");
                for (int i = 0; i < rs.count; ++i) {
                    const TimelinePeriod& r = rs.items[i];
                    char label[TIMELINE_PERIOD_NAME_CAP + 16];
                    snprintf(label, sizeof(label), "%s###period%u", r.name, r.id);
                    if (ImGui::Selectable(label, struct_exp.selected_traj_filter == 2 && r.id == struct_exp.selected_traj_period)) {
                        struct_exp.selected_traj_filter = 2;
                        struct_exp.selected_traj_period = r.id;
                    }
                    if (ImGui::IsItemHovered()) timeline_periods_highlight(r.id);
                }
                ImGui::EndCombo();
            }
        }

        struct_exp.selected_file_format = CLAMP(struct_exp.selected_file_format, 0, (int)ARRAY_SIZE(file_formats) - 1);
        ImGui::Combo("File Format", &struct_exp.selected_file_format, file_formats, (int)ARRAY_SIZE(file_formats));

        //ImGui::PopItemWidth();

        auto extract_atom_indices = [data](int* out_atom_indices, size_t atom_index_cap, const md_system_state_t* state) -> size_t {
            switch (data->structure_export.selected_atom_filter) {
                case 0: { // Active Selection
                    return md_bitset_iter_extract_indices(out_atom_indices, atom_index_cap, md_bitmask_iter(&data->selection.selection_mask));
                } break;
                case 1: { // Visible Atoms
                    return md_bitset_iter_extract_indices(out_atom_indices, atom_index_cap, md_bitmask_iter(&data->representation.visibility_mask));
                } break;
                case 2: { // Query
                    auto& query = data->structure_export.query;
                    if (query.filter.result.is_dynamic) {
                        // What it picks depends on the coordinates: picked again in each frame exported (md_filter
                        // replaces what the mask held, so nothing of the previous frame stays)
                        md_filter_context_t ctx = filter_context(data);
                        ctx.state = state;
                        filter_evaluate(&query.filter, &query.mask, ctx, persistent_alloc);
                    }
                    if (!query.filter.result.success) {
                        VIAMD_LOG_ERROR("Cannot export structure, the query expression is invalid: %s", query.filter.result.error);
                        return 0;
                    } else {
                        return md_bitset_iter_extract_indices(out_atom_indices, atom_index_cap, md_bitset_iter(query.mask));
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
                case 2: { // A timeline period
                    int beg = 0, end = -1;
                    if (!traj_period || !timeline_period_frames(&beg, &end, *traj_period, data->animation.frame, (int)run_num_frames(data))) {
                        beg = 0;
                        end = -1;
                    }
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
            md_bitset_iter_t it = md_bitmask_iter(&data->selection.highlight_mask);
            while (md_bitset_iter_next(&it)) {
                uint64_t idx = it.idx;
                flags[idx] |= AtomBit_Highlighted;
            }
        }
        {
            md_bitset_iter_t it = md_bitmask_iter(&data->selection.selection_mask);
            while (md_bitset_iter_next(&it)) {
                uint64_t idx = it.idx;
                flags[idx] |= AtomBit_Selected;
            }
        }
        {
            md_bitset_iter_t it = md_bitmask_iter(&data->representation.visibility_mask);
            while (md_bitset_iter_next(&it)) {
                uint64_t idx = it.idx;
                flags[idx] |= AtomBit_Visible;
            }
        }
        md_gl_mol_set_atom_flags(data->mold.gl_mol, 0, (uint32_t)sys.atom.count, flags, 0);
    }

    if (data->mold.dirty_gpu_buffers & MolBit_DirtyBonds) {
        md_gl_mol_set_bonds(data->mold.gl_mol, 0, (uint32_t)sys.bond.count, sys.bond.pairs, sizeof(md_atom_pair_t));
        md_gl_mol_set_bond_orders(data->mold.gl_mol, &sys);
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

static void render(ApplicationState* state) {
    bool do_screenshot = !str_empty(state->screenshot.path_to_file);

    uint32_t gbuffer_target_width  = state->app.framebuffer.width;
    uint32_t gbuffer_target_height = state->app.framebuffer.height;
    if (do_screenshot) {
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

    if (state->simulation_box.enabled && state->mold.state.unitcell.flags != 0) {
        mat3_t A = { 0 };
		md_unitcell_A_extract_float(A.elem, &state->mold.state.unitcell);
        immediate::Scope scope(state->gfx.world, "simulation box");
        // Through the turn the coordinates carry, if the orientation is kept: the cell itself cannot hold one
        immediate::box_wireframe(scope, {0,0,0}, {1,1,1}, state->operations.state_rotation * mat4_from_mat3(A), convert_color(state->simulation_box.color));
    }

    {
        immediate::Scope vis_scope(state->gfx.overlay, "visualization");
        immediate::Scope vis_scope_depth(state->gfx.world, "visualization with depth");

        static const md_script_vis_t nothing = {};
        const md_script_vis_t& vis = state->script.vis_shown ? *state->script.vis_shown : nothing;

        if (vis.points) {
            immediate::points(vis_scope, (immediate::Vertex*)vis.points, md_array_size(vis.points), state->script.point_color);
        }

        if (vis.triangles) {
            immediate::triangles(vis_scope, (immediate::Vertex*)vis.triangles, md_array_size(vis.triangles), state->script.triangle_color);
        }

        if (vis.lines) {
            immediate::lines(vis_scope, (immediate::Vertex*)vis.lines, md_array_size(vis.lines), state->script.line_color);
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
    immediate::render(state->gfx.world, params);
    //glDisable(GL_BLEND);
    POP_GPU_SECTION()

    glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);

    if (!use_gfx) {
        PUSH_GPU_SECTION("Selection")
        const bool atom_selection_empty = md_bitmask_empty(&state->selection.selection_mask);
        const bool atom_highlight_empty = md_bitmask_empty(&state->selection.highlight_mask);

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

            postprocessing::blit_color(col_vis);

            glStencilFunc(GL_EQUAL, 0, 0xFF);
            postprocessing::blit_color(state->selection.color.highlight.hidden);
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
    immediate::render(state->gfx.overlay, params);
    POP_GPU_SECTION()

    glDisable(GL_BLEND);
    glEnable(GL_CULL_FACE);

    POP_GPU_SECTION()  // G-buffer

    if (do_screenshot && state->screenshot.hide_gui) {
        // Activate gbuffer to store screenshot
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

    settings.background_color = state->visuals.background.color;
    settings.background_intensity = state->visuals.background.intensity;

    settings.shading.mode = state->visuals.shading.mode;
    settings.shading.ambient = state->visuals.shading.ambient;
    settings.shading.specular = state->visuals.shading.specular;

    settings.ssao.enabled = state->visuals.ssao.enabled;
    settings.ssao.intensity = state->visuals.ssao.intensity;
    settings.ssao.temporal = state->visuals.ssao.temporal;
    settings.ssao.mode = state->visuals.ssao.mode;

    settings.tonemap.enabled = state->visuals.tonemapping.enabled;
    settings.tonemap.mode = state->visuals.tonemapping.tonemapper;
    settings.tonemap.exposure = state->visuals.tonemapping.exposure;
    settings.tonemap.gamma = state->visuals.tonemapping.gamma;

    settings.dof.enabled = state->visuals.dof.enabled;
    settings.dof.focus_depth = state->visuals.dof.focus_depth;
    settings.dof.aperture = state->visuals.dof.aperture * 0.01f;

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
        // Ping-pong: what was written this frame is read as the previous history next frame (no copy needed)
        uint32_t tmp = state->gbuffer.tex.history;
        state->gbuffer.tex.history = state->gbuffer.tex.history_prev;
        state->gbuffer.tex.history_prev = tmp;
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
            if (representation_shown(*state, state->representation.reps[i])) {
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
                if (representation_shown(*state, rep) && rep.type_is_valid) {
                    md_gl_draw_op_t op = {
                        .type = (md_gl_rep_type_t)rep.type,
                        .args = {},
                        .rep = state->representation.reps[i].md_rep,
                        .model_matrix = NULL,
                    };
                    switch (rep.type) {
                    case RepresentationType::SpaceFill:
                        op.args.space_fill.radius_scale = rep.scale.x;
                        break;
                    case RepresentationType::Licorice:
                        op.args.licorice.radius = rep.scale.x;
                        op.args.licorice.color_mode = (md_gl_bond_mode_t)rep.bond_color;
                        op.args.licorice.sharpness = rep.bond_sharpness;
                        op.args.licorice.uniform_color = convert_color(rep.bond_base_color);
                        op.args.licorice.bond_orders = rep.bond_orders;
                        if (rep.tint_scale > 0.0f || rep.saturation < 1.0f) {
                            tint_colors(&op.args.licorice.uniform_color, 1, convert_color(rep.tint_color), rep.tint_scale, rep.saturation);
                        }
                        break;
                    case RepresentationType::BallAndStick:
                        op.args.ball_and_stick.ball_scale = rep.scale.x;
                        op.args.ball_and_stick.stick_radius = rep.scale.y;
                        op.args.ball_and_stick.color_mode = (md_gl_bond_mode_t)rep.bond_color;
                        op.args.ball_and_stick.sharpness = rep.bond_sharpness;
                        op.args.ball_and_stick.uniform_color = convert_color(rep.bond_base_color);
                        op.args.ball_and_stick.bond_orders = rep.bond_orders;
                        if (rep.tint_scale > 0.0f || rep.saturation < 1.0f) {
                            tint_colors(&op.args.ball_and_stick.uniform_color, 1, convert_color(rep.tint_color), rep.tint_scale, rep.saturation);
                        }
                        break;
                    case RepresentationType::Ribbons:
                        op.args.ribbons.width_scale = rep.scale.x;
                        op.args.ribbons.thickness_scale = rep.scale.y;
                        break;
                    case RepresentationType::Cartoon:
                        op.args.cartoon.coil_scale = rep.scale.x;
                        op.args.cartoon.sheet_scale = rep.scale.y;
                        op.args.cartoon.helix_scale = rep.scale.z;
                        break;
                    default:
                        break;
                    }
                    md_array_push(draw_ops, op, frame_alloc);
                }
            } else if (rep.type == RepresentationType::DipoleMoment) {
                if (representation_shown(*state, rep)) {
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
        if (!representation_shown(*state, rep)) continue;
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
                    .roughness = 0.3f,
                    .ior = 1.5f,
                    .mode = state->visuals.shading.mode,
                    .ambient = state->visuals.shading.ambient,
                    .specular = state->visuals.shading.specular,
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

        if (representation_shown(*data, rep) && rep.type_is_valid) {
            md_gl_draw_op_t op = {
                .type = (md_gl_rep_type_t)rep.type,
                .args = {},
                .rep = data->representation.reps[i].md_rep,
                .model_matrix = NULL,
            };
            // The scales lead every struct of the args. The geometry has to be that of the opaque pass, which this
            // draws again with GL_EQUAL depth: the bond orders too.
            MEMCPY(&op.args, &rep.scale, MIN(sizeof(op.args), sizeof(rep.scale)));
            if (rep.type == RepresentationType::Licorice) {
                op.args.licorice.bond_orders = rep.bond_orders;
            } else if (rep.type == RepresentationType::BallAndStick) {
                op.args.ball_and_stick.bond_orders = rep.bond_orders;
            }
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
