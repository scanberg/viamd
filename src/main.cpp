#include <core/md_compiler.h>

#if MD_COMPILER_MSVC
#   ifndef _CRT_SECURE_NO_WARNINGS
#       define _CRT_SECURE_NO_WARNINGS
#   endif
#   pragma warning( disable : 26812 4244 )
#endif

#define IMGUI_DEFINE_MATH_OPERATORS

#define VIAMD_RECOMPUTE_ORBITAL_PER_FRAME 0

#include <md_util.h>
#include <md_gl.h>
#include <md_gfx.h>
#include <md_filter.h>
#include <md_script.h>
#include <md_molecule.h>
#include <md_trajectory.h>
#include <md_xvg.h>
#include <md_csv.h>
#include <md_lammps.h>

#include <core/md_platform.h>
#include <core/md_log.h>
#include <core/md_str.h>
#include <core/md_array.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_ring_allocator.h>
#include <core/md_tracking_allocator.h>
#include <core/md_simd.h>
#include <core/md_os.h>
#include <core/md_base64.h>
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

#include <imgui_widgets.h>
#include <implot_widgets.h>
#include <task_system.h>
#include <color_utils.h>
#include <loader.h>
#include <image.h>
#include <app/application.h>
#include <app/IconsFontAwesome6.h>
#include <event.h>

#include <imgui.h>
#include <imgui_internal.h>

#include <implot.h>
#include <implot_internal.h>
#include <TextEditor.h>
#include <imgui_notify.h>

#include <stdio.h>
#include <bitset>

#include <viamd.h>
#include <serialization_utils.h>

#define MAX_POPULATION_SIZE 256
#define MAX_TEMPORAL_SUBPLOTS 10
#define MAX_DISTRIBUTION_SUBPLOTS 10
#define EXPERIMENTAL_GFX_API 0
#define PICKING_JITTER_HACK 0
#define COMPILATION_TIME_DELAY_IN_SECONDS 1.0
#define NOTIFICATION_DISPLAY_TIME_IN_SECONDS 5.0
#define IR_SEMAPHORE_MAX_COUNT 3
#define MEASURE_EVALUATION_TIME 1
#define FRAME_ALLOCATOR_BYTES MEGABYTES(256)

#define LOG_INFO  MD_LOG_INFO
#define LOG_DEBUG MD_LOG_DEBUG
#define LOG_ERROR MD_LOG_ERROR
#define LOG_SUCCESS(...) ImGui::InsertNotification(ImGuiToast(ImGuiToastType_Success, 6000, __VA_ARGS__))

// Global data for application
static md_allocator_i* frame_alloc = 0; // Linear allocator for scratch data which only is valid for the frame and then is reset
static md_allocator_i* persistent_alloc = 0;

static TextEditor editor {};
static bool use_gfx = false;

constexpr str_t header_snippet = STR_LIT(
R"(#01010110#01001001#01000001#01001101#01000100#01001101#01000001#01001001#01010110#
#                                                                                #
#            VIAMD — Visual Interactive Analysis of Molecular Dynamics           #
#                                                                                #
#                    github: https://github.com/scanberg/viamd                   #
#                 manual: https://github.com/scanberg/viamd/wiki                 #
#                    youtube playlist: https://bit.ly/4aRsPrh                    #
#                                twitter: @VIAMD_                                #
#                                                                                #
#                If you use VIAMD in your research, please cite:                 #
#   "VIAMD: a Software for Visual Interactive Analysis of Molecular Dynamics"    #
#       Robin Skånberg, Ingrid Hotz, Anders Ynnerman, and Mathieu Linares        #
#                 J. Chem. Inf. Model. 2023, 63, 23, 7382–7391                   #
#                   https://doi.org/10.1021/acs.jcim.3c01033                     #
#                                                                                #
#01010110#01001001#01000001#01001101#01000100#01001101#01000001#01001001#01010110#

)");

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

//constexpr ImGuiKey KEY_CONSOLE = ImGuiKey_GraveAccent;
constexpr ImGuiKey KEY_PLAY_PAUSE          = ImGuiKey_Space;
constexpr ImGuiKey KEY_SKIP_TO_PREV_FRAME  = ImGuiKey_LeftArrow;
constexpr ImGuiKey KEY_SKIP_TO_NEXT_FRAME  = ImGuiKey_RightArrow;
constexpr ImGuiKey KEY_RECOMPILE_SHADERS   = ImGuiKey_F5;
constexpr ImGuiKey KEY_SHOW_DEBUG_WINDOW   = ImGuiKey_F11;
constexpr ImGuiKey KEY_SCRIPT_EVALUATE     = ImGuiKey_Enter;
constexpr ImGuiKey KEY_SCRIPT_EVALUATE_MOD = ImGuiMod_Shift;
constexpr ImGuiKey KEY_RECENTER_ON_HIGHLIGHT = ImGuiKey_F1;

constexpr str_t WORKSPACE_FILE_EXTENSION = STR_LIT("via");

constexpr uint32_t PROPERTY_COLORS[] = {4293119554, 4290017311, 4287291314, 4281114675, 4288256763, 4280031971, 4285513725, 4278222847, 4292260554, 4288298346, 4288282623, 4280834481};
constexpr str_t SCRIPT_IMPORT_FILE_EXTENSIONS[] = { STR_LIT("edr"), STR_LIT("xvg"), STR_LIT("csv") };

inline const ImVec4& vec_cast(const vec4_t& v) { return *(const ImVec4*)(&v); }
inline const vec4_t& vec_cast(const ImVec4& v) { return *(const vec4_t*)(&v); }
inline const ImVec2& vec_cast(const vec2_t& v) { return *(const ImVec2*)(&v); }
inline const vec2_t& vec_cast(const ImVec2& v) { return *(const vec2_t*)(&v); }

inline ImVec4& vec_cast(vec4_t& v) { return *(ImVec4*)(&v); }
inline vec4_t& vec_cast(ImVec4& v) { return *(vec4_t*)(&v); }
inline ImVec2& vec_cast(vec2_t& v) { return *(ImVec2*)(&v); }
inline vec2_t& vec_cast(ImVec2& v) { return *(vec2_t*)(&v); }

enum LegendColorMapMode_ {
    LegendColorMapMode_Opaque,
    LegendColorMapMode_Transparent,
    LegendColorMapMode_Split,
};

enum RepBit_ {
    RepBit_DirtyColor   = 0x1,
    RepBit_DirtyFilter  = 0x2
};

enum MarkerType_{
    MarkerType_None,
    MarkerType_Error,
    MarkerType_Warning,
    MarkerType_Visualization,
};

// #struct Structure Declarations

// This is viamd's representation of a property
struct DisplayProperty {
    enum Type {
        Type_Temporal,
        Type_Distribution,
        Type_Volume,
        Type_Count
    };

    enum PlotType {
        PlotType_Line,      // Single line
        PlotType_Area,      // Shaded area
        PlotType_Bars,      // Bar chart
        PlotType_Scatter,   // Scatter plot
        PlotType_Count
    };

    enum ColorType {
        ColorType_Solid,
        ColorType_Colormap,
        ColorType_Count
    };

    // This is the payload passed to getters for display properties
    struct Payload {
        DisplayProperty* display_prop;
        int dim_idx;
    };

    // Callback signature for printing out the value (when hovering with mouse for example)
    typedef int (*PrintValue)(char* buf, size_t buf_cap, int sample_idx, Payload* data);

    struct Histogram {
        int num_bins;
        int dim = 0;
        // Can be multidimensional
        // Total number of entries will be dim * num_bins
        md_array(float) bins = 0;
        double x_min;
        double x_max;
        double y_min;
        double y_max;
        md_allocator_i* alloc;
    };

    Type type = Type_Temporal;

    char label[32] = "";

    ColorType color_type = ColorType_Solid;
    ImVec4 color = {1,1,1,1};
    ImPlotColormap colormap = ImPlotColormap_Plasma;
    float colormap_alpha = 1.0f;

    PlotType plot_type = PlotType_Line;
    ImPlotMarker marker_type = ImPlotMarker_Square;
    float marker_size = 1.0f;
    double bar_width_scale = 1.0;

    // We need two getters to support areas (min / max)
    ImPlotGetter getter[2] = {0,0};
    PrintValue   print_value = 0;

    bool aggregate_histogram = false;

    int dim = 1;                // Number of values per sample
    int num_samples = 0;        // Number of samples (length of x)
    const float* y_values = 0;  // Values (y)
    const float* x_values = 0;  // Corresponding x values

    int num_bins = 128;         // Requested number of bins for histogram

    md_unit_t unit[2] = {md_unit_none(), md_unit_none()};
    char unit_str[2][32] = {"",""};

    const md_script_eval_t* eval = NULL;

    md_script_property_flags_t prop_flags = MD_SCRIPT_PROPERTY_FLAG_NONE;
    const md_script_property_data_t* prop_data = NULL;
    const md_script_vis_payload_o* vis_payload = NULL;

    uint64_t prop_fingerprint = 0;

    // Encodes which temporal subplots this property is visible in
    uint32_t temporal_subplot_mask = 0;

    // Encodes which distribution subplots this property is visible in
    uint32_t distribution_subplot_mask = 0;

    bool show_in_volume = false;
    bool partial_evaluation = false;

    // Encodes which indices of the population to show (if applicable, i.e. dim > 1)
    std::bitset<MAX_POPULATION_SIZE> population_mask = {};

    STATIC_ASSERT(MAX_TEMPORAL_SUBPLOTS     <= sizeof(temporal_subplot_mask) * 8,     "Cannot fit temporal subplot mask");
    STATIC_ASSERT(MAX_DISTRIBUTION_SUBPLOTS <= sizeof(distribution_subplot_mask) * 8, "Cannot fit distribution subplot mask");

    Histogram hist = {};
};

struct LoadParam {
    md_molecule_loader_i*   mol_loader  = NULL;
    md_trajectory_loader_i* traj_loader = NULL;
    str_t file_path = STR_LIT("");
    bool coarse_grained = false;
    const void* mol_loader_arg = NULL;
    LoadTrajectoryFlags traj_loader_flags = 0;
};

const str_t* find_in_arr(str_t str, const str_t arr[], size_t len) {
    for (size_t i = 0; i < len; ++i) {
        if (str_eq(arr[i], str)) {
            return &arr[i];
        }
    }
    return NULL;
}

static inline bool file_queue_empty(const FileQueue* queue) {
    return queue->head == queue->tail;
}

static inline bool file_queue_full(const FileQueue* queue) {
    return (queue->head + 1) % ARRAY_SIZE(queue->arr) == queue->tail;
}

static inline void file_queue_push(FileQueue* queue, str_t path, FileFlags flags = FileFlags_None) {
    ASSERT(queue);
    ASSERT(!file_queue_full(queue));
    int prio = 0;

    str_t ext;
    if (extract_ext(&ext, path) && str_eq(ext, WORKSPACE_FILE_EXTENSION)) {
        prio = 1;
    } else if (load::mol::loader_from_ext(ext)) {
        prio = 2;
    } else if (load::traj::loader_from_ext(ext)) {
        prio = 3;
    } else if (find_in_arr(ext, SCRIPT_IMPORT_FILE_EXTENSIONS, ARRAY_SIZE(SCRIPT_IMPORT_FILE_EXTENSIONS))) {
        prio = 4;
    } else {
        // Unknown extension
        prio = 5;
        flags |= FileFlags_ShowDialogue;
    }

    uint64_t i = queue->head;
    queue->arr[queue->head] = {str_copy(path, queue->ring), flags, prio};
    queue->head = (queue->head + 1) % ARRAY_SIZE(queue->arr);

    // Sort queue based on prio
     while (i != queue->tail && queue->arr[i].prio < queue->arr[(i - 1) % ARRAY_SIZE(queue->arr)].prio) {
        FileQueue::Entry tmp = queue->arr[i];
        queue->arr[i] = queue->arr[(i - 1) % ARRAY_SIZE(queue->arr)];
        queue->arr[(i - 1) % ARRAY_SIZE(queue->arr)] = tmp;
        i = (i - 1) % ARRAY_SIZE(queue->arr);
     }
}

static inline FileQueue::Entry file_queue_front(const FileQueue* queue) {
    ASSERT(!file_queue_empty(queue));
    return queue->arr[queue->tail];
}

static inline FileQueue::Entry file_queue_pop(FileQueue* queue) {
    ASSERT(queue);
    ASSERT(!file_queue_empty(queue));
    FileQueue::Entry front = file_queue_front(queue);
    queue->tail = (queue->tail + 1) % ARRAY_SIZE(queue->arr);
    return front;
}

static inline uint64_t generate_fingerprint() {
    return (uint64_t)md_time_current();
}

static void free_histogram(DisplayProperty::Histogram* hist) {
    ASSERT(hist);
    ASSERT(hist->alloc);
    md_array_free(hist->bins, hist->alloc);
    hist->bins = 0;
}

static void compute_histogram(float* bins, int num_bins, float bin_range_min, float bin_range_max, const float* values, int num_values, float* bin_val_min, float* bin_val_max) {
    MEMSET(bins, 0, sizeof(float) * num_bins);

    const float range_ext = bin_range_max - bin_range_min;
    const float inv_range = 1.0f / range_ext;
    int count = 0;
    for (int i = 0; i < num_values; ++i) {
        if (values[i] < bin_range_min || bin_range_max < values[i]) continue;
        int idx = CLAMP((int)(((values[i] - bin_range_min) * inv_range) * num_bins), 0, num_bins - 1);
        bins[idx] += 1.0f;
        count += 1;
    }

    if (count == 0) {
        if (bin_val_min) *bin_val_min = 0;
        if (bin_val_max) *bin_val_max = 0;
        return;
    }
    
    float min_val = FLT_MAX;
    float max_val = -FLT_MAX;
    const float width = range_ext / num_bins;
    const float scl = 1.0f / (width * count);
    for (int i = 0; i < num_bins; ++i) {
        bins[i] *= scl;
        min_val = MIN(min_val, bins[i]);
        max_val = MAX(max_val, bins[i]);
    }
    
    if (bin_val_min) *bin_val_min = min_val;
    if (bin_val_max) *bin_val_max = max_val;
}

static void compute_histogram_masked(DisplayProperty::Histogram* hist, int num_bins, float value_range_min, float value_range_max, const float* values, int dim, const md_bitfield_t* mask, bool aggregate = false) {
    ASSERT(hist);
    ASSERT(values);
    ASSERT(mask);
    ASSERT(dim > 0);

    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(frame_alloc);
    defer { md_vm_arena_temp_end(temp); };

    hist->dim = aggregate ? 1 : dim;
    md_array_resize(hist->bins, (size_t)(hist->dim * num_bins), hist->alloc);
    MEMSET(hist->bins, 0, md_array_bytes(hist->bins));

    const size_t num_samples = md_bitfield_popcount(mask) * dim;
    if (num_samples == 0) return;

    const float range_ext = value_range_max - value_range_min;
    const float inv_range = range_ext > 0.0f ? 1.0f / range_ext : 0.0f;

    int* count = (int*)md_vm_arena_push_zero_array(frame_alloc, int, hist->dim);

    // We evaluate each frame, one at a time
    md_bitfield_iter_t it = md_bitfield_iter_create(mask);
    while (md_bitfield_iter_next(&it)) {
        const int val_idx = dim * (int)md_bitfield_iter_idx(&it);
        for (int i = 0; i < dim; ++i) {
            const float val = values[val_idx + i];
            if (val < value_range_min || value_range_max < val) continue;
            const int bin_idx = CLAMP((int)(((val - value_range_min) * inv_range) * num_bins), 0, num_bins - 1);

            if (aggregate) {
                hist->bins[bin_idx] += 1.0f;
                count[0] += 1;
            } else {
                hist->bins[num_bins * i + bin_idx] += 1.0f;
                count[i] += 1;
            }
        }
    }

    float min_bin =  FLT_MAX;
    float max_bin = -FLT_MAX;
    const float width = range_ext / num_bins;
    for (int i = 0; i < hist->dim; ++i) {
        const float scl = 1.0f / (width * count[i]);
        for (int j = 0; j < num_bins; ++j) {
            float& val = hist->bins[num_bins * i + j];
            val *= scl;
            min_bin = MIN(min_bin, val);
            max_bin = MAX(max_bin, val);
        }
    }

    hist->num_bins = num_bins;
    hist->x_min = value_range_min;
    hist->x_max = value_range_max;
    hist->y_min = min_bin;
    hist->y_max = max_bin;
}

static void downsample_histogram(float* dst_bins, int num_dst_bins, const float* src_bins, const float* src_weights, int num_src_bins) {
    ASSERT(dst_bins);
    ASSERT(src_bins);
    ASSERT(num_dst_bins <= num_src_bins);

    MEMSET(dst_bins, 0, sizeof(float) * num_dst_bins);

    const int factor = MAX(1, num_src_bins / num_dst_bins);
    for (int dst_idx = 0; dst_idx < num_dst_bins; ++dst_idx) {
        double bin = 0.0;
        double weight = 0.0;
        for (int i = 0; i < factor; ++i) {
            int src_idx = dst_idx * factor + i;
            bin += src_bins[src_idx];
            weight += src_weights ? src_weights[src_idx] : 1.0;
        }
        dst_bins[dst_idx] = (float)(bin / weight);
    }
}

static void scale_histogram(float* bins, const float* weights, int num_bins) {
    ASSERT(bins);
    ASSERT(weights);

    for (int i = 0; i < num_bins; ++i) {
        if (weights[i]) {
            bins[i] /= weights[i];
        }
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

//static void launch_prefetch_job(ApplicationState* data);

static void init_dataset_items(ApplicationState* data);
static void clear_dataset_items(ApplicationState* data);

static void init_display_properties(ApplicationState* data);
static void update_display_properties(ApplicationState* data);

static void update_density_volume(ApplicationState* data);
static void clear_density_volume(ApplicationState* data);

static void interpolate_atomic_properties(ApplicationState* data);
static void update_view_param(ApplicationState* data);
static void reset_view(ApplicationState* data, const md_bitfield_t* target, bool move_camera = false, bool smooth_transition = false);

static void handle_camera_interaction(ApplicationState* data);

//static void init_display_properties(ApplicationState* data);
//static void update_density_volume_texture(ApplicationState* data);
static void handle_picking(ApplicationState* data);

static void fill_gbuffer(ApplicationState* data);
static void apply_postprocessing(const ApplicationState& data);

static void draw_representations_opaque(ApplicationState*);
static void draw_representations_opaque_lean_and_mean(ApplicationState*, uint32_t mask = 0xFFFFFFFFU);
static void draw_representations_transparent(ApplicationState*);

static void draw_load_dataset_window(ApplicationState* data);
static void draw_main_menu(ApplicationState* data);
static void draw_context_popup(ApplicationState* data);
static void draw_selection_query_window(ApplicationState* data);
static void draw_selection_grow_window(ApplicationState* data);
static void draw_animation_window(ApplicationState* data);
static void draw_representations_window(ApplicationState* data);
static void draw_timeline_window(ApplicationState* data);
static void draw_distribution_window(ApplicationState* data);
static void draw_async_task_window(ApplicationState* data);
static void draw_density_volume_window(ApplicationState* data);
static void draw_script_editor_window(ApplicationState* data);
static void draw_dataset_window(ApplicationState* data);
static void draw_debug_window(ApplicationState* data);
static void draw_property_export_window(ApplicationState* data);
static void draw_notifications_window();

static void update_md_buffers(ApplicationState* data);

static void init_molecule_data(ApplicationState* data);
static void init_trajectory_data(ApplicationState* data);

static void interrupt_async_tasks(ApplicationState* data);

static bool load_dataset_from_file(ApplicationState* data, const LoadParam& param);

static void load_workspace(ApplicationState* data, str_t file);
static void save_workspace(ApplicationState* data, str_t file);

static bool export_xvg(const float* column_data[], const char* column_labels[], size_t num_columns, size_t num_rows, str_t filename);
static bool export_csv(const float* column_data[], const char* column_labels[], size_t num_columns, size_t num_rows, str_t filename);

static void create_screenshot(ApplicationState* data);

// Representations
static Representation* create_representation(ApplicationState* data, RepresentationType type = RepresentationType::SpaceFill,
                                             ColorMapping color_mapping = ColorMapping::Cpk, str_t filter = STR_LIT("all"));
static Representation* clone_representation(ApplicationState* data, const Representation& rep);
static void remove_representation(ApplicationState*, int idx);
static void update_representation(ApplicationState*, Representation* rep);
static void update_representation_info(ApplicationState*);
static void update_all_representations(ApplicationState*);
static void init_representation(ApplicationState*, Representation* rep);
static void init_all_representations(ApplicationState*);
static void clear_representations(ApplicationState*);
static void create_default_representations(ApplicationState*);
static void recompute_atom_visibility_mask(ApplicationState*);

// Selections
static Selection* create_selection(ApplicationState* data, str_t name, md_bitfield_t* bf = 0);
static void remove_selection(ApplicationState* data, int idx);
static void clear_selections(ApplicationState* data);

static bool filter_expression(ApplicationState* data, str_t expr, md_bitfield_t* mask, bool* is_dynamic, char* error_str, int error_cap);

static void modify_selection(ApplicationState* data, md_bitfield_t* atom_mask, SelectionOperator op = SelectionOperator::Set) {
    ASSERT(data);
    modify_field(&data->selection.selection_mask, atom_mask, op);
}

static void modify_selection(ApplicationState* data, md_range_t range, SelectionOperator op = SelectionOperator::Set) {
    ASSERT(data);
    modify_field(&data->selection.selection_mask, range, op);
}

static void set_hovered_property(ApplicationState* data, str_t label, int population_idx = -1) {
    data->hovered_display_property_label = label;
    data->hovered_display_property_pop_idx = population_idx;
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
        md_timestamp_t time;
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
            const md_timestamp_t time = md_time_current();
            const uint64_t hash = md_hash64(msg, strlen(msg), 0);

            if (md_time_as_seconds(time - state.time) < 1.0 && hash == state.hash) {
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
                md_mutex_lock(&state.lock);
                defer { md_mutex_unlock(&state.lock); };
                // @NOTE: This needs to be protected with a mutex as it pushes to an internal vector
                ImGui::InsertNotification(ImGuiToast(toast_type, (uint64_t)(NOTIFICATION_DISPLAY_TIME_IN_SECONDS * 1000), msg));
            }
        }
    };

    md_logger_add(&notification_logger);

    ApplicationState data;
    data.allocator.persistent = persistent_alloc;
    data.allocator.frame = frame_alloc;
    data.representation.info.alloc = md_arena_allocator_create(persistent_alloc, MEGABYTES(1));
    data.file_queue.ring = md_ring_allocator_create(md_alloc(persistent_alloc, MEGABYTES(1)), MEGABYTES(1));
    data.mold.mol_alloc  = md_arena_allocator_create(persistent_alloc, MEGABYTES(1));

    md_bitfield_init(&data.selection.selection_mask, persistent_alloc);
    md_bitfield_init(&data.selection.highlight_mask, persistent_alloc);
    md_bitfield_init(&data.selection.query.mask, persistent_alloc);
    md_bitfield_init(&data.selection.grow.mask, persistent_alloc);

    md_bitfield_init(&data.representation.visibility_mask, persistent_alloc);

    md_semaphore_init(&data.script.ir_semaphore, IR_SEMAPHORE_MAX_COUNT);

    // Init platform
    LOG_DEBUG("Initializing GL...");
    if (!application::initialize(&data.app, 0, 0, STR_LIT("VIAMD"))) {
        LOG_ERROR("Could not initialize application...\n");
        return -1;
    }

    data.app.window.vsync = true;
    data.app.file_drop.user_data = &data;
    data.app.file_drop.callback = [](size_t num_files, const str_t file_paths[], void* user_data) {
        ApplicationState* data = (ApplicationState*)user_data;
        ASSERT(data);

        for (int i = 0; i < num_files; ++i) {
            file_queue_push(&data->file_queue, file_paths[i], FileFlags_None);
        }
    };

    LOG_DEBUG("Initializing framebuffer...");
    init_gbuffer(&data.gbuffer, data.app.framebuffer.width, data.app.framebuffer.height);

    for (int i = 0; i < (int)ARRAY_SIZE(data.view.jitter.sequence); ++i) {
        data.view.jitter.sequence[i].x = md_halton(i + 1, 2);
        data.view.jitter.sequence[i].y = md_halton(i + 1, 3);
    }

    // Init subsystems
    LOG_DEBUG("Initializing immediate draw...");
    immediate::initialize();
    LOG_DEBUG("Initializing post processing...");
    postprocessing::initialize(data.gbuffer.width, data.gbuffer.height);
    LOG_DEBUG("Initializing volume...");
    volume::initialize();
    LOG_DEBUG("Initializing task system...");
    const size_t num_threads = VIAMD_NUM_WORKER_THREADS == 0 ? md_os_num_processors() : VIAMD_NUM_WORKER_THREADS;
    task_system::initialize(CLAMP(num_threads, 2, (uint32_t)md_os_num_processors()));

    md_gl_initialize();
    data.mold.gl_shaders                = md_gl_shaders_create(shader_output_snippet);
    data.mold.gl_shaders_lean_and_mean  = md_gl_shaders_create(shader_output_snippet_lean_and_mean);

    viamd::event_system_broadcast_event(viamd::EventType_ViamdInitialize, viamd::EventPayloadType_ApplicationState, &data);

#if EXPERIMENTAL_GFX_API
    md_gfx_initialize(data.gbuffer.width, data.gbuffer.height, 0);
#endif

    ImGui::init_theme();

    editor.SetLanguageDefinition(TextEditor::LanguageDefinition::VIAMD());
    editor.SetPalette(TextEditor::GetDarkPalette());

    {
#ifdef VIAMD_DATASET_DIR
        char exe[1024];
        size_t len = md_path_write_exe(exe, sizeof(exe));
        if (len) {
            md_strb_t sb = md_strb_create(frame_alloc);
            str_t folder;
            extract_folder_path(&folder, {exe, len});
            sb += folder;
            sb += VIAMD_DATASET_DIR "/1ALA-500.pdb";
            replace_char(md_strb_ptr(sb), md_strb_len(sb), '\\', '/');
            str_t path = md_strb_to_str(sb);
            if (md_path_is_valid(path)) {
                // @NOTE: We want explicitly to disable writing of cache files for the default dataset
                // The motivation is that the dataset may reside in a shared folder on the system that has no write access.
                file_queue_push(&data.file_queue, path, FileFlags_DisableCacheWrite);
                editor.SetText("s1 = resname(\"ALA\")[2:8];\nd1 = distance(10,30);\na1 = angle(2,1,3) in resname(\"ALA\");\nr = rdf(element('C'), element('H'), 10.0);\nv = sdf(s1, element('H'), 10.0);\n{lin,plan,iso} = shape_weights(all);");
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
                    file_queue_push(&data.file_queue, path);
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
    while (!data.app.window.should_close) {
        application::update(&data.app);
        
        // This needs to happen first (in imgui events) to enable docking of imgui windows
#if VIAMD_IMGUI_ENABLE_DOCKSPACE
        ImGui::CreateDockspace();
#endif

        const size_t num_frames  = md_trajectory_num_frames(data.mold.traj);
        const size_t last_frame  = num_frames > 0 ? num_frames - 1 : 0;
        const double   max_frame = (double)last_frame;

        if (!file_queue_empty(&data.file_queue) && !data.load_dataset.show_window) {
            FileQueue::Entry e = file_queue_pop(&data.file_queue);

            str_t ext;
            extract_ext(&ext, e.path);
            const str_t* res = 0;

            if (str_eq_ignore_case(ext, WORKSPACE_FILE_EXTENSION)) {
                load_workspace(&data, e.path);
            } else if ((res = find_in_arr(ext, SCRIPT_IMPORT_FILE_EXTENSIONS, ARRAY_SIZE(SCRIPT_IMPORT_FILE_EXTENSIONS)))) {
                char buf[1024];
                str_t base_path = {};
                if (data.files.workspace[0] != '\0') {
                    base_path = str_from_cstr(data.files.workspace);
                } else if (data.files.trajectory[0] != '\0') {
                    base_path = str_from_cstr(data.files.trajectory);
                } else if (data.files.molecule[0] != '\0') {
                    base_path = str_from_cstr(data.files.molecule);
                } else {
                    md_path_write_cwd(buf, sizeof(buf));
                    base_path = str_from_cstr(buf);
                }

                str_t rel_path = md_path_make_relative(base_path, e.path, frame_alloc);
                MD_LOG_DEBUG("Attempting to make relative path from '" STR_FMT "' to '" STR_FMT "'", STR_ARG(base_path), STR_ARG(e.path));
                MD_LOG_DEBUG("Relative path: '" STR_FMT "'", STR_ARG(rel_path));
                if (!str_empty(rel_path)) {
                    snprintf(buf, sizeof(buf), "table = import(\"%.*s\");\n", STR_ARG(rel_path));
                    TextEditor::Coordinates pos = editor.GetCursorPosition();
                    pos.mLine += 1;
                    editor.SetCursorPosition({0,0});
                    editor.InsertText(buf);
                    editor.SetCursorPosition(pos);
                }
            } else {
                load::LoaderState state = {};
                bool success = load::init_loader_state(&state, e.path, frame_alloc);

                if (!success || (e.flags & FileFlags_ShowDialogue) || (state.flags & LoaderStateFlag_RequiresDialogue)) {
                    data.load_dataset = LoadDatasetWindowState();
                    str_copy_to_char_buf(data.load_dataset.path_buf, sizeof(data.load_dataset.path_buf), e.path);
                    data.load_dataset.path_changed = true;
                    data.load_dataset.show_window = true;
                    data.load_dataset.coarse_grained = e.flags & FileFlags_CoarseGrained;
                } else if (success) {
                    LoadParam param = {};
                    param.mol_loader  = state.mol_loader;
                    param.traj_loader = state.traj_loader;
                    param.file_path      = e.path;
                    param.coarse_grained = e.flags & FileFlags_CoarseGrained;
                    param.mol_loader_arg = state.mol_loader_arg;
                    param.traj_loader_flags = (e.flags & FileFlags_DisableCacheWrite) ? LoadTrajectoryFlag_DisableCacheWrite : 0;
                    if (load_dataset_from_file(&data, param)) {
                        data.animation = {};
                        if (param.mol_loader) {
                            if (!data.settings.keep_representations) {
                                clear_representations(&data);
                                create_default_representations(&data);
                            }
                            interpolate_atomic_properties(&data);
                            data.mold.dirty_buffers |= MolBit_ClearVelocity;
                            reset_view(&data, &data.representation.visibility_mask, true, false);
                        }
                    }
                }
            }
        }

        viamd::event_system_broadcast_event(viamd::EventType_ViamdFrameTick, viamd::EventPayloadType_ApplicationState, &data);

        // GUI
        if (data.show_script_window) draw_script_editor_window(&data);
        if (data.load_dataset.show_window) draw_load_dataset_window(&data);
        if (data.representation.show_window) draw_representations_window(&data);
        if (data.density_volume.show_window) draw_density_volume_window(&data);
        if (data.distributions.show_window) draw_distribution_window(&data);
        if (data.timeline.show_window) draw_timeline_window(&data);
        if (data.dataset.show_window) draw_dataset_window(&data);
        if (data.selection.query.show_window) draw_selection_query_window(&data);
        if (data.selection.grow.show_window) draw_selection_grow_window(&data);
        if (data.show_property_export_window) draw_property_export_window(&data);
        if (data.show_debug_window) draw_debug_window(&data);

        data.selection.selecting = false;

        //ImGui::ShowDemoWindow();

        handle_camera_interaction(&data);
        camera_animate(&data.view.camera, data.view.animation.target_orientation, data.view.animation.target_position, data.view.animation.target_distance, data.app.timing.delta_s);
        data.visuals.dof.focus_depth = data.view.camera.focus_distance;

        update_view_param(&data);

        ImGuiWindow* win = ImGui::GetCurrentContext()->HoveredWindow;
        if (win && strcmp(win->Name, "Main interaction window") == 0) {
            set_hovered_property(&data,  STR_LIT(""));
        }

        draw_context_popup(&data);
        draw_async_task_window(&data);
        if (data.animation.show_window) draw_animation_window(&data);
        draw_main_menu(&data);
        draw_notifications_window();

        // Capture non-window specific keyboard events
        if (!ImGui::GetIO().WantCaptureKeyboard) {
#if EXPERIMENTAL_GFX_API
            if (ImGui::IsKeyPressed(ImGuiKey_F1)) {
                use_gfx = !use_gfx;
            }
#endif
            if (ImGui::IsKeyDown(KEY_SCRIPT_EVALUATE_MOD) && ImGui::IsKeyPressed(KEY_SCRIPT_EVALUATE)) {
                data.script.eval_init = true;
            }

            if (ImGui::IsKeyPressed(KEY_SHOW_DEBUG_WINDOW)) {
                data.show_debug_window = true;
            }

            if (ImGui::IsKeyPressed(KEY_RECOMPILE_SHADERS)) {
                LOG_INFO("Recompiling shaders and re-initializing volume");
                postprocessing::initialize(data.gbuffer.width, data.gbuffer.height);
                volume::initialize();
                md_gl_shaders_destroy(data.mold.gl_shaders);
                data.mold.gl_shaders = md_gl_shaders_create(shader_output_snippet);
            }

            if (ImGui::IsKeyPressed(KEY_PLAY_PAUSE)) {
                if (data.animation.mode == PlaybackMode::Stopped) {
                    if (data.animation.frame == max_frame && data.animation.fps > 0) {
                        data.animation.frame = 0;
                    } else if (data.animation.frame == 0 && data.animation.fps < 0) {
                        data.animation.frame = max_frame;
                    }
                }

                switch (data.animation.mode) {
                    case PlaybackMode::Playing:
                        data.animation.mode = PlaybackMode::Stopped;
                        break;
                    case PlaybackMode::Stopped:
                        data.animation.mode = PlaybackMode::Playing;
                        break;
                    default:
                        ASSERT(false);
                }

                data.mold.dirty_buffers |= MolBit_DirtyPosition;   // Update previous position to not get motion trail when paused
            }

            if (ImGui::IsKeyPressed(KEY_SKIP_TO_PREV_FRAME) || ImGui::IsKeyPressed(KEY_SKIP_TO_NEXT_FRAME)) {
                double step = ImGui::IsKeyDown(ImGuiMod_Ctrl) ? 10.0 : 1.0;
                if (ImGui::IsKeyPressed(KEY_SKIP_TO_PREV_FRAME)) step = -step;
                data.animation.frame = CLAMP(data.animation.frame + step, 0.0, max_frame);
            }
        }

        if (data.representation.atom_visibility_mask_dirty) {
            recompute_atom_visibility_mask(&data);
            data.representation.atom_visibility_mask_dirty = false;
        }

        if (data.animation.mode == PlaybackMode::Playing) {
            data.animation.frame += data.app.timing.delta_s * data.animation.fps;
            data.animation.frame = CLAMP(data.animation.frame, 0.0, max_frame);
            if (data.animation.frame >= max_frame) {
                data.animation.mode = PlaybackMode::Stopped;
                data.animation.frame = max_frame;
            } else if (data.animation.frame <= 0) {
                data.animation.mode = PlaybackMode::Stopped;
                data.animation.frame = 0;
            }

            if (data.settings.prefetch_frames) {
                if (!task_system::task_is_running(data.tasks.prefetch_frames)) {
                    uint32_t traj_frames = (uint32_t)md_trajectory_num_frames(data.mold.traj);
                    if (traj_frames > 0 && load::traj::num_cache_frames(data.mold.traj) < traj_frames) {
                        uint32_t frame_beg = 0;
                        uint32_t frame_end = 0;
                        // @NOTE: This is certainly something which can be improved upon.
                        // It prefetches frames in the direction of the animation.
                        // In a more optimal case, it should never yield until it catches up with the number of frames it expects to have as a buffer.

                        int look_ahead = CLAMP((int)load::traj::num_cache_frames(data.mold.traj) / 2, 1, 10);
                    
                        if (data.animation.fps > 0) {
                            frame_beg = (uint32_t)CLAMP((int)data.animation.frame             , 0, (int)traj_frames);
                            frame_end = (uint32_t)CLAMP((int)data.animation.frame + look_ahead, 0, (int)traj_frames);
                        } else {
                            frame_beg = (uint32_t)CLAMP((int)data.animation.frame - look_ahead, 0, (int)traj_frames);
                            frame_end = (uint32_t)CLAMP((int)data.animation.frame             , 0, (int)traj_frames);
                        }
                        if (frame_beg != frame_end) {
                            uint32_t frame_count = frame_end - frame_beg;
                            data.tasks.prefetch_frames = task_system::create_pool_task(STR_LIT("##Prefetch Frames"), frame_count, [&data, frame_offset = frame_beg](uint32_t frame_beg, uint32_t frame_end, uint32_t thread_num) {
                                (void)thread_num;                                
                                for (uint32_t i = frame_offset + frame_beg; i < frame_offset + frame_end; ++i) {
                                    md_trajectory_load_frame(data.mold.traj, i, 0, 0, 0, 0);
                                }
                            });
                            task_system::enqueue_task(data.tasks.prefetch_frames);
                        }
                    }
                }
            }
        }

        {
            static auto prev_frame = data.animation.frame;
            if (data.animation.frame != prev_frame) {
                time_changed = true;
                prev_frame = data.animation.frame;
            }
            else {
                time_changed = false;
            }
        }

        if (data.timeline.filter.temporal_window.enabled) {
            const double pre_beg = data.timeline.filter.beg_frame;
            const double pre_end = data.timeline.filter.end_frame;
            const double half_window_ext = data.timeline.filter.temporal_window.extent_in_frames * 0.5;
            data.timeline.filter.beg_frame = CLAMP(round(data.animation.frame - half_window_ext), 0.0, max_frame);
            data.timeline.filter.end_frame = CLAMP(round(data.animation.frame + half_window_ext), 0.0, max_frame);
            if (data.script.ir && (data.timeline.filter.beg_frame != pre_beg || data.timeline.filter.end_frame != pre_end)) {
                data.script.evaluate_filt = true;
            }
        }

        if (data.timeline.filter.enabled) {
            static auto prev_filter_beg = data.timeline.filter.beg_frame;
            static auto prev_filter_end = data.timeline.filter.end_frame;
            if (data.timeline.filter.beg_frame != prev_filter_beg || data.timeline.filter.end_frame != prev_filter_end) {
                prev_filter_beg = data.timeline.filter.beg_frame;
                prev_filter_end = data.timeline.filter.end_frame;
                data.timeline.filter.fingerprint = generate_fingerprint();
            }
        }

        if (time_changed) {
            time_stopped = false;

            PUSH_CPU_SECTION("Interpolate Position")
            if (data.mold.traj) {
                interpolate_atomic_properties(&data);
            }
            POP_CPU_SECTION()

            PUSH_CPU_SECTION("Update dynamic representations")
            for (size_t i = 0; i < md_array_size(data.representation.reps); ++i) {
                auto& rep = data.representation.reps[i];
                if (!rep.enabled) continue;
                if (rep.dynamic_evaluation || rep.color_mapping == ColorMapping::SecondaryStructure) {
                    update_representation(&data, &rep);
                }
            }
            POP_CPU_SECTION()
        } else {
            if (!time_stopped) {
                time_stopped = true;
                data.mold.dirty_buffers |= MolBit_DirtyPosition;
            }
        }

        {
            std::string text = editor.GetText();
            data.script.text = {text.c_str(), text.length()};
            data.script.text_hash = md_hash64(text.c_str(), text.length(), 0);
        }

        if (data.script.compile_ir) {
            data.script.time_since_last_change += data.app.timing.delta_s;

            editor.ClearMarkers();
            editor.ClearErrorMarkers();

            if (data.script.time_since_last_change > COMPILATION_TIME_DELAY_IN_SECONDS) {
                // We cannot recompile while it is evaluating.
                // Need to interrupt and wait for tasks to finish.
                if (data.script.full_eval) md_script_eval_interrupt(data.script.full_eval);
                if (data.script.filt_eval) md_script_eval_interrupt(data.script.filt_eval);

                // Try aquire all semaphores
                if (md_semaphore_try_aquire_n(&data.script.ir_semaphore, IR_SEMAPHORE_MAX_COUNT)) {
                    defer {
                        md_semaphore_release_n(&data.script.ir_semaphore, IR_SEMAPHORE_MAX_COUNT);
                        update_all_representations(&data);
                    };

                    // Now we hold all semaphores for the script
                    data.script.compile_ir = false;
                    data.script.time_since_last_change = 0;
                    
                    data.script.ir = md_script_ir_create(persistent_alloc);

                    std::string src = editor.GetText();
                    str_t src_str {src.data(), src.length()};

                    char buf[1024];
                    size_t len = md_path_write_cwd(buf, sizeof(buf));
                    str_t old_cwd = {buf, len};
                    defer {
                        md_path_set_cwd(old_cwd);
                    };
                    
                    str_t cwd = {};
                    if (data.files.workspace[0] != '\0') {
                        extract_folder_path(&cwd, str_from_cstr(data.files.workspace));
                    } else if (data.files.trajectory[0] != '\0') {
                        extract_folder_path(&cwd, str_from_cstr(data.files.trajectory));
                    } else if (data.files.molecule[0] != '\0') {
                        extract_folder_path(&cwd, str_from_cstr(data.files.molecule));
                    }
                    if (!str_empty(cwd)) {
                        md_path_set_cwd(cwd);
                    }
                    
                    if (src_str) {
                        const size_t num_stored_selections = md_array_size(data.selection.stored_selections);                       
                        for (size_t i = 0; i < num_stored_selections; ++i) {
                            str_t name = str_from_cstr(data.selection.stored_selections[i].name);
                            const md_bitfield_t* bf = &data.selection.stored_selections[i].atom_mask;
                            md_script_ir_add_identifier_bitfield(data.script.ir, name, bf);
                        }
                        md_script_ir_compile_from_source(data.script.ir, src_str, &data.mold.mol, data.mold.traj, NULL);

                        const size_t num_errors = md_script_ir_num_errors(data.script.ir);
                        const md_log_token_t* errors = md_script_ir_errors(data.script.ir);
                        
                        for (size_t i = 0; i < num_errors; ++i) {
                            TextEditor::Marker marker = {0};
                            auto first = editor.GetCharacterCoordinates(errors[i].range.beg);
                            auto last  = editor.GetCharacterCoordinates(errors[i].range.end);
                            marker.type = MarkerType_Error;
                            marker.begCol = first.mColumn;
                            marker.endCol = last.mColumn;
                            marker.prio = INT32_MAX;   // Ensures marker is rendered on top
                            marker.bgColor = IM_COL32(255, 0, 0, 128);
                            marker.hoverBgColor = 0;
                            marker.text = std::string(errors[i].text.ptr, errors[i].text.len);
                            marker.payload = errors[i].context;
                            marker.line = first.mLine + 1;
                            editor.AddMarker(marker);
                        }

                        const size_t num_warnings = md_script_ir_num_warnings(data.script.ir);
                        const md_log_token_t* warnings = md_script_ir_warnings(data.script.ir);
                        for (size_t i = 0; i < num_warnings; ++i) {
                            TextEditor::Marker marker = {0};
                            auto first = editor.GetCharacterCoordinates(warnings[i].range.beg);
                            auto last  = editor.GetCharacterCoordinates(warnings[i].range.end);
                            marker.type = MarkerType_Warning;
                            marker.begCol = first.mColumn;
                            marker.endCol = last.mColumn;
                            marker.prio = INT32_MAX - 1;   // Ensures marker is rendered on top (but bellow an error)
                            marker.bgColor = IM_COL32(255, 255, 0, 128);
                            marker.hoverBgColor = 0;
                            marker.text = std::string(warnings[i].text.ptr, warnings[i].text.len);
                            marker.payload = errors[i].context;
                            marker.line = first.mLine + 1;
                            editor.AddMarker(marker);
                        }

                        const size_t num_tokens = md_script_ir_num_vis_tokens(data.script.ir);
                        const md_script_vis_token_t* vis_tokens = md_script_ir_vis_tokens(data.script.ir);
                        for (size_t i = 0; i < num_tokens; ++i) {
                            const md_script_vis_token_t& vis_tok = vis_tokens[i];
                            TextEditor::Marker marker = {0};
                            auto first = editor.GetCharacterCoordinates(vis_tok.range.beg);
                            auto last  = editor.GetCharacterCoordinates(vis_tok.range.end);
                            marker.type = MarkerType_Visualization;
                            marker.begCol = first.mColumn;
                            marker.endCol = last.mColumn;
                            marker.prio = vis_tok.depth;
                            marker.bgColor = 0;
                            marker.hoverBgColor = IM_COL32(255, 255, 255, 128);
                            marker.text = std::string(vis_tok.text.ptr, vis_tok.text.len);
                            marker.payload = (void*)vis_tok.payload;
                            marker.line = first.mLine + 1;
                            editor.AddMarker(marker);
                        }

                        if (md_script_ir_valid(data.script.ir)) {
                            uint64_t ir_figerprint = md_script_ir_fingerprint(data.script.ir);
                            if (data.script.ir_fingerprint != ir_figerprint) {
                                data.script.ir_fingerprint = ir_figerprint;
                            }
                        } else {
                            md_script_ir_free(data.script.ir);
                            data.script.ir = nullptr;
                        }
                    }
                }
            }
        }

        if (num_frames > 0) {
            if (data.script.eval_init) {
                if (task_system::task_is_running(data.tasks.evaluate_full)) md_script_eval_interrupt(data.script.full_eval);
                if (task_system::task_is_running(data.tasks.evaluate_filt)) md_script_eval_interrupt(data.script.filt_eval);
                    
                if (task_system::task_is_running(data.tasks.evaluate_full) == false &&
                    task_system::task_is_running(data.tasks.evaluate_filt) == false) {
                    data.script.eval_init = false;

                    if (data.script.full_eval) {
                        md_script_eval_free(data.script.full_eval);
                    }
                    if (data.script.filt_eval) {
                        md_script_eval_free(data.script.filt_eval);
                    }
                
                    if (md_script_ir_valid(data.script.ir)) {
                        if (data.script.ir != data.script.eval_ir) {
                            md_script_ir_free(data.script.eval_ir);
                            data.script.eval_ir = data.script.ir;
                        }
                        data.script.full_eval = md_script_eval_create(num_frames, data.script.eval_ir, persistent_alloc);
                        data.script.filt_eval = md_script_eval_create(num_frames, data.script.eval_ir, persistent_alloc);
                    }

                    init_display_properties(&data);

                    data.script.evaluate_filt = true;
                    data.script.evaluate_full = true;
                }
            }

            if (data.script.full_eval && data.script.evaluate_full) {
                if (task_system::task_is_running(data.tasks.evaluate_full)) {
                    md_script_eval_interrupt(data.script.full_eval);
                } else {
                    if (md_script_ir_valid(data.script.eval_ir) &&
                        md_script_eval_ir_fingerprint(data.script.full_eval) == md_script_ir_fingerprint(data.script.eval_ir))
                    {
                        data.script.evaluate_full = false;
                        md_script_eval_clear_data(data.script.full_eval);

                        if (md_script_ir_property_count(data.script.eval_ir) > 0) {
                            data.tasks.evaluate_full = task_system::create_pool_task(STR_LIT("Eval Full"), (uint32_t)num_frames, [&data](uint32_t frame_beg, uint32_t frame_end, uint32_t thread_num) {
                                (void)thread_num;
                                md_script_eval_frame_range(data.script.full_eval, data.script.eval_ir, &data.mold.mol, data.mold.traj, frame_beg, frame_end);
                            });
                            
#if MEASURE_EVALUATION_TIME
                            uint64_t time = (uint64_t)md_time_current();
                            task_system::ID time_task = task_system::create_pool_task(STR_LIT("##Time Eval Full"), [t0 = time]() {
                                uint64_t t1 = md_time_current();
                                double s = md_time_as_seconds(t1 - t0);
                                LOG_INFO("Evaluation completed in: %.3fs", s);
                            });
#endif
                            task_system::set_task_dependency(time_task, data.tasks.evaluate_full);
                            task_system::enqueue_task(data.tasks.evaluate_full);
                        }
                    }
                }
            }

            if (data.script.filt_eval && data.script.evaluate_filt && data.timeline.filter.enabled) {
                if (task_system::task_is_running(data.tasks.evaluate_filt)) {
                    md_script_eval_interrupt(data.script.filt_eval);
                } else {
                    //if (md_semaphore_try_aquire(&data.script.ir_semaphore)) {
                        if (md_script_ir_valid(data.script.eval_ir) &&
                            md_script_eval_ir_fingerprint(data.script.filt_eval) == md_script_ir_fingerprint(data.script.eval_ir))
                        {
                            data.script.evaluate_filt = false;
                            md_script_eval_clear_data(data.script.filt_eval);

                            if (md_script_ir_property_count(data.script.eval_ir) > 0) {
                                const uint32_t traj_frames = (uint32_t)md_trajectory_num_frames(data.mold.traj);
                                const uint32_t beg_frame = CLAMP((uint32_t)data.timeline.filter.beg_frame, 0, traj_frames-1);
                                const uint32_t end_frame = CLAMP((uint32_t)data.timeline.filter.end_frame + 1, beg_frame + 1, traj_frames);
                                if (beg_frame != end_frame) {
                                    data.tasks.evaluate_filt = task_system::create_pool_task(STR_LIT("Eval Filt"), end_frame - beg_frame, [offset = beg_frame, &data](uint32_t beg, uint32_t end, uint32_t thread_num) {
                                        (void)thread_num;
                                        md_script_eval_frame_range(data.script.filt_eval, data.script.eval_ir, &data.mold.mol, data.mold.traj, offset + beg, offset + end);
                                    });
                                    task_system::enqueue_task(data.tasks.evaluate_filt);
                                }
                            }
                            
                            /*
                            task_system::pool_enqueue(STR_LIT("##Release IR Semaphore"), [](void* user_data)
                                {
                                    ApplicationState* data = (ApplicationState*)user_data;
                                    md_semaphore_release(&data->script.ir_semaphore);
                                }, &data, data.tasks.evaluate_filt);
                                */
                        }
                        //else {
                        //    md_semaphore_release(&data.script.ir_semaphore);
                        //}
                    //}
                }
            }
        }

        if (ImGui::IsKeyPressed(KEY_RECENTER_ON_HIGHLIGHT)) {
            reset_view(&data, &data.selection.highlight_mask, true, true);
        }

        // Resize Framebuffer
        if ((data.gbuffer.width != data.app.framebuffer.width || data.gbuffer.height != data.app.framebuffer.height) &&
            (data.app.framebuffer.width != 0 && data.app.framebuffer.height != 0)) {
            init_gbuffer(&data.gbuffer, data.app.framebuffer.width, data.app.framebuffer.height);
            postprocessing::initialize(data.gbuffer.width, data.gbuffer.height);
        }

        // The motivation for doing this is to reduce the frequency at which we invalidate and upload the atom flag field to the GPU
        // For large systems, this can be a costly operation: Consider a system of 100'000'000 atoms
        // The the size of each bitfield to represent the mask would be 12.5 MB
        // Given a throughput of 10GB/s for hash64 results in a time of 1.25 ms per bitfield.

        uint64_t v_hash = data.representation.visibility_mask_hash;
        uint64_t h_hash = md_bitfield_hash64(&data.selection.highlight_mask, 0);
        uint64_t s_hash = md_bitfield_hash64(&data.selection.selection_mask, 0);
        uint64_t f_hash = v_hash ^ h_hash ^ s_hash;

        // These represent the 'current' state so we can compare against it to see if they were modified
        static uint64_t highlight_hash = 0;
        static uint64_t selection_hash = 0;
        static uint64_t flag_hash = 0;

        if (h_hash != highlight_hash) {
            highlight_hash = h_hash;
            // enqueue the event here and do not broadcast (stall)
            viamd::event_system_enqueue_event(viamd::EventType_ViamdHoverMaskChanged, viamd::EventPayloadType_ApplicationState, &data);
        }

        if (s_hash != selection_hash) {
            selection_hash = s_hash;
            // enqueue the event here and do not broadcast (stall)
            viamd::event_system_enqueue_event(viamd::EventType_ViamdSelectionMaskChanged, viamd::EventPayloadType_ApplicationState, &data);
        }

        if (f_hash != flag_hash) {
            flag_hash = f_hash;
            data.mold.dirty_buffers |= MolBit_DirtyFlags;
        }

        update_md_buffers(&data);
        update_display_properties(&data);

        handle_picking(&data);
        clear_gbuffer(&data.gbuffer);
        fill_gbuffer(&data);

        // Activate backbuffer
        glDisable(GL_DEPTH_TEST);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glViewport(0, 0, data.app.framebuffer.width, data.app.framebuffer.height);
        glDrawBuffer(GL_BACK);
        glClear(GL_COLOR_BUFFER_BIT);

        apply_postprocessing(data);

        // Render Screenshot of backbuffer without GUI here
        if (data.screenshot.hide_gui && !str_empty(data.screenshot.path_to_file)) {
            create_screenshot(&data);
            str_free(data.screenshot.path_to_file, persistent_alloc);
            data.screenshot.path_to_file = {};
        }

        PUSH_GPU_SECTION("Imgui render")
        application::render_imgui(&data.app);
        POP_GPU_SECTION()

        if (!data.screenshot.hide_gui && !str_empty(data.screenshot.path_to_file)) {
            create_screenshot(&data);
            data.screenshot.path_to_file = {};
            str_free(data.screenshot.path_to_file, persistent_alloc);
        }

        viamd::event_system_process_event_queue();
        task_system::execute_main_task_queue();

        // Reset frame allocator
        md_vm_arena_reset(frame_alloc);

        // Swap buffers
        application::swap_buffers(&data.app);
    }

    interrupt_async_tasks(&data);

    viamd::event_system_broadcast_event(viamd::EventType_ViamdShutdown);

    // shutdown subsystems
    LOG_DEBUG("Shutting down immediate draw...");
    immediate::shutdown();
    LOG_DEBUG("Shutting down post processing...");
    postprocessing::shutdown();
    LOG_DEBUG("Shutting down volume...");
    volume::shutdown();
    LOG_DEBUG("Shutting down task system...");
    task_system::shutdown();

    destroy_gbuffer(&data.gbuffer);
    application::shutdown(&data.app);

    return 0;
}

static void init_dataset_items(ApplicationState* data) {
    clear_dataset_items(data);
    if (data->mold.mol.atom.count == 0) return;

    for (size_t i = 0; i < data->mold.mol.chain.count; ++i) {
        DatasetItem item = {};
        str_t str = LBL_TO_STR(data->mold.mol.chain.id[i]);
        snprintf(item.label, sizeof(item.label), "%.*s", (int)str.len, str.ptr);
        snprintf(item.query, sizeof(item.query), "chain(%d)", (int)(i+1));
        size_t atom_count = md_chain_atom_count(data->mold.mol.chain, i);
        item.count = 1;
        item.fraction = atom_count / (float)data->mold.mol.atom.count;
        md_array_push(data->dataset.chains, item, persistent_alloc);
    }

    for (size_t i = 0; i < data->mold.mol.residue.count; ++i) {
        const float fraction_size = md_residue_atom_count(data->mold.mol.residue, i) / (float)data->mold.mol.atom.count;
        {
            // Do resname
            str_t resname = LBL_TO_STR(data->mold.mol.residue.name[i]);
            DatasetItem* item = 0;
            for (size_t j = 0; j < md_array_size(data->dataset.residue_names); ++j) {
                if (strcmp(data->dataset.residue_names[j].label, resname.ptr) == 0) {
                    item = &data->dataset.residue_names[j];
                    break;
                }
            }
            if (!item) {
                DatasetItem it = {};
                snprintf(it.label, sizeof(it.label), "%.*s", (int)resname.len, resname.ptr);
                snprintf(it.query, sizeof(it.query), "resname('%.*s')", (int)resname.len, resname.ptr);
                it.count = 0;
                it.fraction = 0;
                md_array_push(data->dataset.residue_names, it, persistent_alloc);
                item = md_array_last(data->dataset.residue_names);
            }
            item->count += 1;
            item->fraction += fraction_size;
        }
    }

    for (size_t i = 0; i < data->mold.mol.atom.count; ++i) {
        {
            // Do atom label
            str_t label = LBL_TO_STR(data->mold.mol.atom.type[i]);
            DatasetItem* item = 0;
            for (size_t j = 0; j < md_array_size(data->dataset.atom_types); ++j) {
                if (strcmp(data->dataset.atom_types[j].label, label.ptr) == 0) {
                    item = &data->dataset.atom_types[j];
                    break;
                }
            }
            if (!item) {
                DatasetItem it = {};
                snprintf(it.label, sizeof(it.label), "%.*s", (int)label.len, label.ptr);
                snprintf(it.query, sizeof(it.query), "type('%.*s')", (int)label.len, label.ptr);
                it.count = 0;
                it.fraction = 0;
                md_array_push(data->dataset.atom_types, it, persistent_alloc);
                item = md_array_last(data->dataset.atom_types);
            }
            item->count += 1;
        }
    }

    for (size_t i = 0; i < md_array_size(data->dataset.atom_types); ++i) {
        data->dataset.atom_types[i].fraction = data->dataset.atom_types[i].count / (float)data->mold.mol.atom.count;
    }
}

static void clear_dataset_items(ApplicationState* data) {
    md_array_free(data->dataset.chains, persistent_alloc);
    md_array_free(data->dataset.residue_names, persistent_alloc);
    md_array_free(data->dataset.atom_types, persistent_alloc);
}

static void display_property_copy_param_from_old(DisplayProperty& item, const DisplayProperty* old_items, int64_t num_old_items) {
    // See if we have a matching item in the old list
    for (int64_t i = 0; i < num_old_items; ++i) {
        if (strcmp(item.label, old_items[i].label) == 0 && item.type == old_items[i].type) {
            // Copy relevant parameters from existing item which we want to be persistent
            item.temporal_subplot_mask      = old_items[i].temporal_subplot_mask;
            item.distribution_subplot_mask  = old_items[i].distribution_subplot_mask;
            item.color                      = old_items[i].color;
            item.temporal_subplot_mask      = old_items[i].temporal_subplot_mask;
            item.distribution_subplot_mask  = old_items[i].distribution_subplot_mask;
            item.show_in_volume             = old_items[i].show_in_volume;
            item.plot_type                  = old_items[i].plot_type;
            item.colormap_alpha			    = old_items[i].colormap_alpha;
            item.colormap                   = old_items[i].colormap;
            item.color_type                 = old_items[i].color_type;
            item.marker_size                = old_items[i].marker_size;
            item.marker_type                = old_items[i].marker_type;
            break;
        }
    }
}

static void init_display_properties(ApplicationState* data) {
    DisplayProperty* new_items = 0;
    DisplayProperty* old_items = data->display_properties;

    const md_script_ir_t* ir = data->script.eval_ir;

    const md_script_eval_t* evals[2] = {
        data->script.full_eval,
        data->script.filt_eval
    };

    const str_t eval_labels[2] = {
        {},
        STR_LIT("filt"),
    };

    for (size_t eval_idx = 0; eval_idx < ARRAY_SIZE(evals); ++eval_idx) {
        const md_script_eval_t* eval = evals[eval_idx];
        const size_t num_props = md_script_ir_property_count(ir);
        const str_t* prop_names = md_script_ir_property_names(ir);
        str_t eval_label = eval_labels[eval_idx];

        const bool partial_evaluation = (eval_idx > 0);

        for (size_t i = 0; i < num_props; ++i) {
            str_t prop_name = prop_names[i];
            md_script_property_flags_t prop_flags = md_script_ir_property_flags(ir, prop_name);
            const md_script_property_data_t* prop_data = md_script_eval_property_data(eval, prop_name);

            if (!prop_data) {
                MD_LOG_DEBUG("Failed to extract property data from property!");
                continue;
            }

            DisplayProperty item;
            if (!str_empty(eval_label)) {
                snprintf(item.label, sizeof(item.label), STR_FMT " " STR_FMT, STR_ARG(prop_name), STR_ARG(eval_label));
            } else {
                snprintf(item.label, sizeof(item.label), STR_FMT, STR_ARG(prop_name));
            }
            item.color = ImGui::ColorConvertU32ToFloat4(PROPERTY_COLORS[i % ARRAY_SIZE(PROPERTY_COLORS)]);
            item.unit[0] = prop_data->unit[0];
            item.unit[1] = prop_data->unit[1];
            item.prop_flags = prop_flags;
            item.prop_data = prop_data;
            item.vis_payload = md_script_ir_property_vis_payload(ir, prop_name);
            item.eval = eval;
            item.prop_fingerprint = 0;
            item.population_mask.set();
            item.temporal_subplot_mask = 0;
            item.distribution_subplot_mask = 0;
            item.hist = {};
            item.hist.alloc = persistent_alloc;
            item.partial_evaluation = partial_evaluation;

            md_unit_print(item.unit_str[0], sizeof(item.unit_str), item.unit[0]);
            md_unit_print(item.unit_str[1], sizeof(item.unit_str), item.unit[1]);

            if (prop_flags & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) {
                // Create a special distribution from the temporal (since we can)
                {
                    DisplayProperty item_dist_raw = item;
                    item_dist_raw.type = DisplayProperty::Type_Distribution;
                    item_dist_raw.plot_type = DisplayProperty::PlotType_Line;
                    item_dist_raw.unit[0] = item.unit[1];
                    item_dist_raw.unit[1] = md_unit_none();
                    MEMCPY(item.unit_str[0], item.unit_str[1], sizeof(item.unit_str[0]));

                    item_dist_raw.getter[0] = [](int sample_idx, void* payload) -> ImPlotPoint {
                        DisplayProperty::Payload* data = (DisplayProperty::Payload*)payload;

                        int num_bins = data->display_prop->hist.num_bins;
                        const double x_min = data->display_prop->hist.x_min;
                        const double x_max = data->display_prop->hist.x_max;
                        const double x_scl = (x_max - x_min) / (num_bins);
                        const double x_off = x_min + 0.5 * x_scl;
                        const double x = x_off + sample_idx * x_scl;
                        return ImPlotPoint(x, 0);
                    };
                    item_dist_raw.getter[1] = [](int sample_idx, void* payload) -> ImPlotPoint {
                        DisplayProperty::Payload* data = (DisplayProperty::Payload*)payload;

                        int num_bins = data->display_prop->hist.num_bins;
                        const double x_min = data->display_prop->hist.x_min;
                        const double x_max = data->display_prop->hist.x_max;
                        const double x_scl = (x_max - x_min) / (num_bins);
                        const double x_off = x_min + 0.5 * x_scl;
                        const double x = x_off + sample_idx * x_scl;
                        const double y = data->display_prop->hist.bins[data->dim_idx * num_bins + sample_idx];
                        return ImPlotPoint(x, y);
                    };
                    display_property_copy_param_from_old(item_dist_raw, old_items, md_array_size(old_items));
                    md_array_push(new_items, item_dist_raw, frame_alloc);

                    if (prop_data->dim[1] > 1) {
                        DisplayProperty item_dist_agg = item_dist_raw;
                        item_dist_agg.type = DisplayProperty::Type_Distribution;
                        snprintf(item_dist_agg.label, sizeof(item_dist_agg.label), "%s (agg)", item.label);
                        item_dist_agg.aggregate_histogram = true;
                        display_property_copy_param_from_old(item_dist_agg, old_items, md_array_size(old_items));
                        md_array_push(new_items, item_dist_agg, frame_alloc);
                    }
                }

                // Now do all of the real temporal ones
                if (!partial_evaluation) {
                    item.num_samples = (int)md_array_size(data->timeline.x_values);
                    item.x_values = data->timeline.x_values;
                    item.y_values = item.prop_data->values;

                    DisplayProperty item_raw = item;
                    item_raw.dim        = prop_data->dim[1];
                    item_raw.plot_type  = DisplayProperty::PlotType_Line;
                    item_raw.getter[0]  = [](int sample_idx, void* payload) -> ImPlotPoint {
                        DisplayProperty::Payload* data = (DisplayProperty::Payload*)payload;
                        int dim_idx = data->dim_idx;
                        int dim = data->display_prop->dim;
                        const float* y_values = data->display_prop->prop_data->values;
                        const float* x_values = data->display_prop->x_values;
                        return ImPlotPoint(x_values[sample_idx], y_values[sample_idx * dim + dim_idx]);
                    };
                    display_property_copy_param_from_old(item_raw, old_items, md_array_size(old_items));
                    md_array_push(new_items, item_raw, frame_alloc);

                    if (prop_data->aggregate) {
                        // Create 'pseudo' display properties which maps to the aggregate data
                        DisplayProperty item_mean = item;
                        snprintf(item_mean.label, sizeof(item_mean.label), "%s (mean)", item.label);
                        item_mean.dim = 1;
                        item_mean.y_values = item_mean.prop_data->aggregate->population_mean;
                        item_mean.plot_type = DisplayProperty::PlotType_Line;
                        item_mean.getter[0] = [](int sample_idx, void* payload) -> ImPlotPoint {
                            DisplayProperty* data = ((DisplayProperty::Payload*)payload)->display_prop;
                            const float* y_values = data->prop_data->aggregate->population_mean;
                            const float* x_values = data->x_values;
                            return ImPlotPoint(x_values[sample_idx], y_values[sample_idx]);
                        };
                        item_mean.print_value = [](char* buf, size_t cap, int sample_idx, DisplayProperty::Payload* payload) -> int {
                            const float* y_mean = payload->display_prop->prop_data->aggregate->population_mean;
                            return snprintf(buf, cap, "%.2f", y_mean[sample_idx]);
                        };
                        display_property_copy_param_from_old(item_mean, old_items, md_array_size(old_items));
                        md_array_push(new_items, item_mean, frame_alloc);

                        DisplayProperty item_var = item;
                        snprintf(item_var.label, sizeof(item_var.label), "%s (var)", item.label);
                        item_var.dim = 1;
                        item_var.y_values = item_mean.prop_data->aggregate->population_var;
                        item_var.color.w *= 0.4f;
                        item_var.plot_type = DisplayProperty::PlotType_Area;
                        item_var.getter[0] = [](int sample_idx, void* payload) -> ImPlotPoint {
                            DisplayProperty* data = ((DisplayProperty::Payload*)payload)->display_prop;
                            const float* y_mean = data->prop_data->aggregate->population_mean;
                            const float* y_var  = data->prop_data->aggregate->population_var;
                            const float* x_values = data->x_values;
                            return ImPlotPoint(x_values[sample_idx], y_mean[sample_idx] - y_var[sample_idx]);
                            };
                        item_var.getter[1] = [](int sample_idx, void* payload) -> ImPlotPoint {
                            DisplayProperty* data = ((DisplayProperty::Payload*)payload)->display_prop;
                            const float* y_mean = data->prop_data->aggregate->population_mean;
                            const float* y_var  = data->prop_data->aggregate->population_var;
                            const float* x_values = data->x_values;
                            return ImPlotPoint(x_values[sample_idx], y_mean[sample_idx] + y_var[sample_idx]);
                            };
                        item_var.print_value = [](char* buf, size_t cap, int sample_idx, DisplayProperty::Payload* payload) -> int {
                            const float* y_var = payload->display_prop->prop_data->aggregate->population_var;
                            return snprintf(buf, cap, "%.2f", y_var[sample_idx]);
                            };
                        display_property_copy_param_from_old(item_var, old_items, md_array_size(old_items));
                        md_array_push(new_items, item_var, frame_alloc);

                        DisplayProperty item_ext = item;
                        snprintf(item_ext.label, sizeof(item_ext.label), "%s (min/max)", item.label);
                        item_ext.dim = 2;
                        item_ext.y_values = (const float*)item_mean.prop_data->aggregate->population_ext;
                        item_ext.color.w *= 0.2f;
                        item_ext.plot_type = DisplayProperty::PlotType_Area;
                        item_ext.getter[0] = [](int sample_idx, void* payload) -> ImPlotPoint {
                            DisplayProperty* data = ((DisplayProperty::Payload*)payload)->display_prop;
                            const vec2_t* y_ext = data->prop_data->aggregate->population_ext;
                            const float* x_values = data->x_values;
                            return ImPlotPoint(x_values[sample_idx], y_ext[sample_idx].x);
                            };
                        item_ext.getter[1] = [](int sample_idx, void* payload) -> ImPlotPoint {
                            DisplayProperty* data = ((DisplayProperty::Payload*)payload)->display_prop;
                            const vec2_t* y_ext = data->prop_data->aggregate->population_ext;
                            const float* x_values = data->x_values;
                            return ImPlotPoint(x_values[sample_idx], y_ext[sample_idx].y);
                            };
                        item_ext.print_value = [](char* buf, size_t cap, int sample_idx, DisplayProperty::Payload* payload) -> int {
                            const vec2_t* y_ext = payload->display_prop->prop_data->aggregate->population_ext;
                            return snprintf(buf, cap, "%.2f, %.2f", y_ext[sample_idx].x, y_ext[sample_idx].y);
                            };
                        display_property_copy_param_from_old(item_ext, old_items, md_array_size(old_items));
                        md_array_push(new_items, item_ext, frame_alloc);
                    }
                }
            } else if (prop_flags & MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) {
                DisplayProperty item_dist = item;
                item_dist.type = DisplayProperty::Type_Distribution;
                item_dist.plot_type = DisplayProperty::PlotType_Line;
                item_dist.getter[0] = [](int sample_idx, void* payload) -> ImPlotPoint {
                    DisplayProperty::Payload* data = (DisplayProperty::Payload*)payload;

                    int num_bins = data->display_prop->hist.num_bins;
                    const double x_min = data->display_prop->hist.x_min;
                    const double x_max = data->display_prop->hist.x_max;
                    const double x_scl = (x_max - x_min) / (num_bins);
                    const double x_off = x_min + 0.5 * x_scl;
                    const double x = x_off + sample_idx * x_scl;
                    return ImPlotPoint(x, 0);
                    };
                item_dist.getter[1] = [](int sample_idx, void* payload) -> ImPlotPoint {
                    DisplayProperty::Payload* data = (DisplayProperty::Payload*)payload;

                    int num_bins = data->display_prop->hist.num_bins;
                    const double x_min = data->display_prop->hist.x_min;
                    const double x_max = data->display_prop->hist.x_max;
                    const double x_scl = (x_max - x_min) / (num_bins);
                    const double x_off = x_min + 0.5 * x_scl;
                    const double x = x_off + sample_idx * x_scl;
                    const double y = data->display_prop->hist.bins[data->dim_idx * num_bins + sample_idx];
                    return ImPlotPoint(x, y);
                    };
                display_property_copy_param_from_old(item_dist, old_items, md_array_size(old_items));
                md_array_push(new_items, item_dist, frame_alloc);
            } else if (prop_flags & MD_SCRIPT_PROPERTY_FLAG_VOLUME) {
                item.type = DisplayProperty::Type_Volume;
                item.show_in_volume = false;
                display_property_copy_param_from_old(item, old_items, md_array_size(old_items));
                md_array_push(new_items, item, frame_alloc);
            }
        }
    }

    for (size_t i = 0; i < md_array_size(old_items); ++i) {
        free_histogram(&old_items[i].hist);
    }

    md_array_resize(data->display_properties, md_array_size(new_items), persistent_alloc);
    MEMCPY(data->display_properties, new_items, md_array_size(new_items) * sizeof(DisplayProperty));
}

static void update_display_properties(ApplicationState* data) {
    ASSERT(data);

    for (size_t i = 0; i < md_array_size(data->display_properties); ++i) {
        DisplayProperty& dp = data->display_properties[i];
        if (dp.type == DisplayProperty::Type_Distribution) {
            if (dp.prop_fingerprint != dp.prop_data->fingerprint || dp.num_bins != dp.hist.num_bins) {
                dp.prop_fingerprint = dp.prop_data->fingerprint;
        
                if (dp.prop_flags & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) {
                    DisplayProperty::Histogram& hist = dp.hist;
                    compute_histogram_masked(&hist, dp.num_bins, dp.prop_data->min_range[0], dp.prop_data->max_range[0], dp.prop_data->values, dp.prop_data->dim[1], md_script_eval_frame_mask(dp.eval), dp.aggregate_histogram);
                }
                else if (dp.prop_flags & MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) {
                    DisplayProperty::Histogram& hist = dp.hist;
                    md_array_resize(hist.bins, (size_t)dp.num_bins, hist.alloc);
                    hist.num_bins = dp.num_bins;
                    hist.x_min = dp.prop_data->min_range[0];
                    hist.x_max = dp.prop_data->max_range[0];
                    hist.y_min = dp.prop_data->min_range[1];
                    hist.y_max = dp.prop_data->max_range[1];
                    hist.dim = 1;
                    downsample_histogram(hist.bins, hist.num_bins, dp.prop_data->values, dp.prop_data->weights, dp.prop_data->dim[2]);
                }
            }
        }
    }
}

static void update_density_volume(ApplicationState* data) {
    if (data->density_volume.dvr.tf.dirty) {
        data->density_volume.dvr.tf.dirty = false;
        // Update colormap texture
        volume::compute_transfer_function_texture_simple_ramp(&data->density_volume.dvr.tf.id, data->density_volume.dvr.tf.colormap, data->density_volume.dvr.tf.alpha_scale);
    }

    int64_t selected_property = -1;
    for (size_t i = 0; i < md_array_size(data->display_properties); ++i) {
        const DisplayProperty& dp = data->display_properties[i];
        if (dp.type == DisplayProperty::Type_Volume && dp.show_in_volume) {
            selected_property = i;
            break;
        }
    }

    const md_script_property_data_t* prop_data = 0;
    const md_script_vis_payload_o* vis_payload = 0;
    uint64_t data_fingerprint = 0;

    static int64_t s_selected_property = 0;
    if (s_selected_property != selected_property) {
        s_selected_property = selected_property;
        data->density_volume.dirty_vol = true;
        data->density_volume.dirty_rep = true;
    }

    if (selected_property != -1) {
        prop_data = data->display_properties[selected_property].prop_data;
        vis_payload = data->display_properties[selected_property].vis_payload;
        data_fingerprint = data->display_properties[selected_property].prop_data->fingerprint;
    }
    data->density_volume.show_density_volume = selected_property != -1;

    static uint64_t s_script_fingerprint = 0;
    if (s_script_fingerprint != md_script_ir_fingerprint(data->script.eval_ir)) {
        s_script_fingerprint = md_script_ir_fingerprint(data->script.eval_ir);
        data->density_volume.dirty_vol = true;
        data->density_volume.dirty_rep = true;
    }

    static uint64_t s_data_fingerprint = 0;
    if (s_data_fingerprint != data_fingerprint) {
        s_data_fingerprint = data_fingerprint;
        data->density_volume.dirty_vol = true;
    }

    static double s_frame = 0;
    if (s_frame != data->animation.frame) {
        s_frame = data->animation.frame;
        data->density_volume.dirty_rep = true;
    }

    if (data->density_volume.dirty_rep) {
        if (prop_data && vis_payload) {
            data->density_volume.dirty_rep = false;
            size_t num_reps = 0;
            bool result = false;
            md_script_vis_t vis = {};

            if (md_script_ir_valid(data->script.eval_ir)) {
                md_script_vis_init(&vis, frame_alloc);
                md_script_vis_ctx_t ctx = {
                    .ir = data->script.eval_ir,
                    .mol = &data->mold.mol,
                    .traj = data->mold.traj,
                };
                result = md_script_vis_eval_payload(&vis, vis_payload, 0, &ctx, MD_SCRIPT_VISUALIZE_SDF);
            }

            if (result) {
                if (vis.sdf.extent) {
                    const float s = vis.sdf.extent;
                    vec3_t min_aabb = { -s, -s, -s };
                    vec3_t max_aabb = { s, s, s };
                    data->density_volume.model_mat = volume::compute_model_to_world_matrix(min_aabb, max_aabb);
                    data->density_volume.voxel_spacing = vec3_t{2*s / prop_data->dim[1], 2*s / prop_data->dim[2], 2*s / prop_data->dim[3]};
                }
                num_reps = md_array_size(vis.sdf.structures);
            }

            // We need to limit this for performance reasons
            num_reps = MIN(num_reps, 100);

            const size_t old_size = md_array_size(data->density_volume.gl_reps);
            if (data->density_volume.gl_reps) {
                // Only free superflous entries
                for (size_t i = num_reps; i < old_size; ++i) {
                    md_gl_rep_destroy(data->density_volume.gl_reps[i]);
                }
            }
            md_array_resize(data->density_volume.gl_reps, num_reps, persistent_alloc);
            md_array_resize(data->density_volume.rep_model_mats, num_reps, persistent_alloc);

            for (size_t i = old_size; i < num_reps; ++i) {
                // Only init new entries
                data->density_volume.gl_reps[i] = md_gl_rep_create(data->mold.gl_mol);
            }

            const auto& mol = data->mold.mol;
            auto& rep = data->density_volume.rep;
            const size_t num_colors = data->mold.mol.atom.count;
            const size_t num_bytes = sizeof(uint32_t) * num_colors;
            uint32_t* colors = (uint32_t*)md_vm_arena_push(frame_alloc, num_bytes);
            defer { md_vm_arena_pop(frame_alloc, num_bytes); };

            switch (rep.colormap) {
            case ColorMapping::Uniform:
                color_atoms_uniform(colors, mol.atom.count, rep.color);
                break;
            case ColorMapping::Cpk:
                color_atoms_cpk(colors, mol.atom.count, mol);
                break;
            case ColorMapping::AtomLabel:
                color_atoms_type(colors, mol.atom.count, mol);
                break;
            case ColorMapping::AtomIndex:
                color_atoms_idx(colors, mol.atom.count, mol);
                break;
            case ColorMapping::ResName:
                color_atoms_res_name(colors, mol.atom.count, mol);
                break;
            case ColorMapping::ResIndex:
                color_atoms_res_id(colors, mol.atom.count, mol);
                break;
            case ColorMapping::ChainId:
                color_atoms_chain_id(colors, mol.atom.count, mol);
                break;
            case ColorMapping::ChainIndex:
                color_atoms_chain_idx(colors, mol.atom.count, mol);
                break;
            case ColorMapping::SecondaryStructure:
                color_atoms_sec_str(colors, mol.atom.count, mol);
                break;
            default:
                ASSERT(false);
                break;
            }

            for (size_t i = 0; i < num_reps; ++i) {
                filter_colors(colors, num_colors, &vis.sdf.structures[i]);
                md_gl_rep_set_color(data->density_volume.gl_reps[i], 0, (uint32_t)num_colors, colors, 0);
                data->density_volume.rep_model_mats[i] = vis.sdf.matrices[i];
            }
        }
    }

    if (data->density_volume.dirty_vol) {
        if (prop_data) {
            data->density_volume.dirty_vol = false;
            if (!data->density_volume.volume_texture.id) {
                int dim[3] = { prop_data->dim[1], prop_data->dim[2], prop_data->dim[3] };
                gl::init_texture_3D(&data->density_volume.volume_texture.id, dim[0], dim[1], dim[2], GL_R16F);
                MEMCPY(data->density_volume.volume_texture.dim, dim, sizeof(dim));
                data->density_volume.volume_texture.max_value = prop_data->max_value;
            }
            gl::set_texture_3D_data(data->density_volume.volume_texture.id, prop_data->values, GL_R32F);
        }
    }
}

static void clear_density_volume(ApplicationState* state) {
    md_array_shrink(state->density_volume.gl_reps, 0);
    md_array_shrink(state->density_volume.rep_model_mats, 0);
    state->density_volume.model_mat = {0};
}

static void interpolate_atomic_properties(ApplicationState* state) {
    ASSERT(state);
    auto& mol = state->mold.mol;
    const auto& traj = state->mold.traj;

    if (!mol.atom.count || !md_trajectory_num_frames(traj)) return;

    const int64_t last_frame = MAX(0LL, (int64_t)md_trajectory_num_frames(traj) - 1);
    // This is not actually time, but the fractional frame representation
    const double time = CLAMP(state->animation.frame, 0.0, double(last_frame));

    // Scaling factor for cubic spline
    const int64_t frame = (int64_t)time;
    const int64_t nearest_frame = CLAMP((int64_t)(time + 0.5), 0LL, last_frame);

    const int64_t frames[4] = {
        MAX(0LL, frame - 1),
        MAX(0LL, frame),
        MIN(frame + 1, last_frame),
        MIN(frame + 2, last_frame)
    };

    const size_t num_threads = task_system::pool_num_threads();

    const size_t stride = ALIGN_TO(mol.atom.count, 16);    // The interploation uses SIMD vectorization without bounds, so we make sure there is no overlap between the data segments
    const size_t bytes = stride * sizeof(float) * 3 * 4;
    
    // The number of atoms to be processed per thread when divided into chunks
    const uint32_t grain_size = 1024;

    md_vm_arena_temp_t tmp = md_vm_arena_temp_begin(frame_alloc);
    defer { md_vm_arena_temp_end(tmp); };

    void* mem = md_vm_arena_push(frame_alloc, bytes);

    struct Payload {
        ApplicationState* state;
        float s;
        float t;
        InterpolationMode mode;

        int64_t nearest_frame;
        int64_t frames[4];
        md_trajectory_frame_header_t headers[4];
        md_unit_cell_t unit_cell;

        size_t count;

        float* src_x[4];
        float* src_y[4];
        float* src_z[4];

        float* dst_x;
        float* dst_y;
        float* dst_z;

        vec3_t* aabb_min;
        vec3_t* aabb_max;
    };

    const InterpolationMode mode = (frames[1] != frames[2]) ? state->animation.interpolation : InterpolationMode::Nearest;

    Payload payload = {
        .state = state,
        .s = 1.0f - CLAMP(state->animation.tension, 0.0f, 1.0f),
        .t = (float)fractf(time),
        .mode = mode,
        .nearest_frame = nearest_frame,
        .frames = { frames[0], frames[1], frames[2], frames[3]},
        .count = mol.atom.count,
        .src_x = { (float*)mem + stride * 0, (float*)mem + stride * 1, (float*)mem + stride * 2,  (float*)mem + stride * 3 },
        .src_y = { (float*)mem + stride * 4, (float*)mem + stride * 5, (float*)mem + stride * 6,  (float*)mem + stride * 7 },
        .src_z = { (float*)mem + stride * 8, (float*)mem + stride * 9, (float*)mem + stride * 10, (float*)mem + stride * 11},
        .dst_x = mol.atom.x,
        .dst_y = mol.atom.y,
        .dst_z = mol.atom.z,
        .aabb_min = (vec3_t*)md_vm_arena_push(frame_alloc, num_threads * sizeof(vec3_t)),
        .aabb_max = (vec3_t*)md_vm_arena_push(frame_alloc, num_threads * sizeof(vec3_t)),
    };

    // This holds the chain of tasks we are about to submit
    task_system::ID tasks[16] = {0};
    int num_tasks = 0;

    switch (mode) {
        case InterpolationMode::Nearest: {
            task_system::ID load_task = task_system::create_pool_task(STR_LIT("## Load Frame"),[data = &payload]() {
                md_trajectory_frame_header_t header;
                md_trajectory_load_frame(data->state->mold.traj, data->nearest_frame, &header, data->dst_x, data->dst_y, data->dst_z);
                MEMCPY(&data->unit_cell, &header.unit_cell, sizeof(md_unit_cell_t));
            });

            tasks[num_tasks++] = load_task;
            break;
        }
        case InterpolationMode::Linear: {
            task_system::ID load_task = task_system::create_pool_task(STR_LIT("## Load Frame"), 2, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                for (uint32_t i = range_beg; i < range_end; ++i) {
                    md_trajectory_load_frame(data->state->mold.traj, data->frames[i+1], &data->headers[i], data->src_x[i], data->src_y[i], data->src_z[i]);
                }
            });

            task_system::ID interp_unit_cell_task = task_system::create_pool_task(STR_LIT("## Interp Unit Cell Data"), [data = &payload]() {
                if ((data->headers[0].unit_cell.flags & MD_UNIT_CELL_FLAG_ORTHO) && (data->headers[1].unit_cell.flags & MD_UNIT_CELL_FLAG_ORTHO)) {
                    double ext_x = lerp(data->headers[0].unit_cell.basis[0][0], data->headers[1].unit_cell.basis[0][0], data->t);
                    double ext_y = lerp(data->headers[0].unit_cell.basis[1][1], data->headers[1].unit_cell.basis[1][1], data->t);
                    double ext_z = lerp(data->headers[0].unit_cell.basis[2][2], data->headers[1].unit_cell.basis[2][2], data->t);
                    data->unit_cell = md_util_unit_cell_from_extent(ext_x, ext_y, ext_z);
                } else if ( (data->headers[0].unit_cell.flags & MD_UNIT_CELL_FLAG_TRICLINIC) || (data->headers[1].unit_cell.flags & MD_UNIT_CELL_FLAG_TRICLINIC)) {
                    data->unit_cell.basis = lerp(data->headers[0].unit_cell.basis, data->headers[1].unit_cell.basis, data->t);
                    data->unit_cell.inv_basis = mat3_inverse(data->state->mold.mol.unit_cell.basis);
                }
            });

            task_system::ID interp_coord_task = task_system::create_pool_task(STR_LIT("## Interp Coord Data"), (uint32_t)mol.atom.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                size_t count = range_end - range_beg;
                float* dst_x = data->dst_x + range_beg;
                float* dst_y = data->dst_y + range_beg;
                float* dst_z = data->dst_z + range_beg;
                const float* src_x[2] = { data->src_x[0] + range_beg, data->src_x[1] + range_beg};
                const float* src_y[2] = { data->src_y[0] + range_beg, data->src_y[1] + range_beg};
                const float* src_z[2] = { data->src_z[0] + range_beg, data->src_z[1] + range_beg};

                md_util_interpolate_linear(dst_x, dst_y, dst_z, src_x, src_y, src_z, count, &data->unit_cell, data->t);
            }, grain_size);

            tasks[num_tasks++] = load_task;
            tasks[num_tasks++] = interp_unit_cell_task;
            tasks[num_tasks++] = interp_coord_task;

            break;
        }
        case InterpolationMode::CubicSpline: {
            task_system::ID load_task = task_system::create_pool_task(STR_LIT("## Load Frame"), 4, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                for (uint32_t i = range_beg; i < range_end; ++i) {
                    md_trajectory_load_frame(data->state->mold.traj, data->frames[i], &data->headers[i], data->src_x[i], data->src_y[i], data->src_z[i]);
                }
            });

            task_system::ID interp_unit_cell_task = task_system::create_pool_task(STR_LIT("## Interp Unit Cell Data"), [data = &payload]() {
                if ((data->headers[0].unit_cell.flags & MD_UNIT_CELL_FLAG_ORTHO) &&
                    (data->headers[1].unit_cell.flags & MD_UNIT_CELL_FLAG_ORTHO) &&
                    (data->headers[2].unit_cell.flags & MD_UNIT_CELL_FLAG_ORTHO) &&
                    (data->headers[3].unit_cell.flags & MD_UNIT_CELL_FLAG_ORTHO))
                {
                    double ext_x = cubic_spline(data->headers[0].unit_cell.basis[0][0], data->headers[1].unit_cell.basis[0][0], data->headers[2].unit_cell.basis[0][0], data->headers[3].unit_cell.basis[0][0], data->t);
                    double ext_y = cubic_spline(data->headers[0].unit_cell.basis[1][1], data->headers[1].unit_cell.basis[1][1], data->headers[2].unit_cell.basis[1][1], data->headers[3].unit_cell.basis[1][1], data->t);
                    double ext_z = cubic_spline(data->headers[0].unit_cell.basis[2][2], data->headers[1].unit_cell.basis[2][2], data->headers[2].unit_cell.basis[2][2], data->headers[3].unit_cell.basis[2][2], data->t);
                    data->unit_cell = md_util_unit_cell_from_extent(ext_x, ext_y, ext_z);
                } else if ( (data->headers[0].unit_cell.flags & MD_UNIT_CELL_FLAG_TRICLINIC) ||
                            (data->headers[1].unit_cell.flags & MD_UNIT_CELL_FLAG_TRICLINIC) ||
                            (data->headers[2].unit_cell.flags & MD_UNIT_CELL_FLAG_TRICLINIC) ||
                            (data->headers[3].unit_cell.flags & MD_UNIT_CELL_FLAG_TRICLINIC))
                {
                    data->unit_cell.flags = data->headers[0].unit_cell.flags;
                    data->unit_cell.basis = cubic_spline(data->headers[0].unit_cell.basis, data->headers[1].unit_cell.basis, data->headers[2].unit_cell.basis, data->headers[3].unit_cell.basis, data->t, data->s);
                    data->unit_cell.inv_basis = mat3_inverse(data->state->mold.mol.unit_cell.basis);
                }
            });

            task_system::ID interp_coord_task = task_system::create_pool_task(STR_LIT("## Interp Coord Data"), (uint32_t)mol.atom.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                size_t count = range_end - range_beg;
                float* dst_x = data->dst_x + range_beg;
                float* dst_y = data->dst_y + range_beg;
                float* dst_z = data->dst_z + range_beg;
                const float* src_x[4] = { data->src_x[0] + range_beg, data->src_x[1] + range_beg, data->src_x[2] + range_beg, data->src_x[3] + range_beg};
                const float* src_y[4] = { data->src_y[0] + range_beg, data->src_y[1] + range_beg, data->src_y[2] + range_beg, data->src_y[3] + range_beg};
                const float* src_z[4] = { data->src_z[0] + range_beg, data->src_z[1] + range_beg, data->src_z[2] + range_beg, data->src_z[3] + range_beg};

                md_util_interpolate_cubic_spline(dst_x, dst_y, dst_z, src_x, src_y, src_z, count, &data->unit_cell, data->t, data->s);
            }, grain_size);

            tasks[num_tasks++] = load_task;
            tasks[num_tasks++] = interp_unit_cell_task;
            tasks[num_tasks++] = interp_coord_task;
            
            break;
        }
        default:
            ASSERT(false);
            break;
    }

    if (state->operations.apply_pbc) {
        task_system::ID pbc_task = task_system::create_pool_task(STR_LIT("## Apply PBC"), (uint32_t)mol.atom.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            (void)thread_num;
            size_t count = range_end - range_beg;
            float* x = data->dst_x + range_beg;
            float* y = data->dst_y + range_beg;
            float* z = data->dst_z + range_beg;
            md_util_pbc(x, y, z, 0, count, &data->unit_cell);
        });
        tasks[num_tasks++] = pbc_task;
    } 
    if (state->operations.unwrap_structures) {
        size_t num_structures = md_index_data_num_ranges(mol.structure);
        task_system::ID unwrap_task = task_system::create_pool_task(STR_LIT("## Unwrap Structures"), (uint32_t)num_structures, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            (void)thread_num;
            for (uint32_t i = range_beg; i < range_end; ++i) {
                int32_t* s_idx = md_index_range_beg(data->state->mold.mol.structure, i);
                size_t   s_len = md_index_range_size(data->state->mold.mol.structure, i);
                md_util_unwrap(data->dst_x, data->dst_y, data->dst_z, s_idx, s_len, &data->unit_cell);
            }
        });
        tasks[num_tasks++] = unwrap_task;
    }

    {
        task_system::ID aabb_task = task_system::create_pool_task(STR_LIT("## Compute AABB"), (uint32_t)mol.atom.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            size_t count = range_end - range_beg;
            const float* x = data->state->mold.mol.atom.x + range_beg;
            const float* y = data->state->mold.mol.atom.y + range_beg;
            const float* z = data->state->mold.mol.atom.z + range_beg;
            const float* r = data->state->mold.mol.atom.radius + range_beg;

            vec3_t aabb_min = vec3_set1(FLT_MAX);
            vec3_t aabb_max = vec3_set1(-FLT_MAX);
            md_util_aabb_compute(aabb_min.elem, aabb_max.elem, x, y, z, r, 0, count);

            data->aabb_min[thread_num] = aabb_min;
            data->aabb_max[thread_num] = aabb_max;
        });
        tasks[num_tasks++] = aabb_task;
    }

    if (mol.protein_backbone.angle) {
        switch (mode) {
            case InterpolationMode::Nearest: {
                task_system::ID angle_task = task_system::create_pool_task(STR_LIT("## Compute Backbone Angles"), [data = &payload]() {
                    const md_backbone_angles_t* src_angles[2] = {
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[1],
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[2],
                    };
                    const md_backbone_angles_t* src_angle = data->t < 0.5f ? src_angles[0] : src_angles[1];
                    MEMCPY(data->state->mold.mol.protein_backbone.angle, src_angle, data->state->mold.mol.protein_backbone.count * sizeof(md_backbone_angles_t));
                });

                tasks[num_tasks++] = angle_task;
                break;
            }
            case InterpolationMode::Linear: {
                task_system::ID angle_task = task_system::create_pool_task(STR_LIT("## Compute Backbone Angles"), (uint32_t)mol.protein_backbone.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                    (void)thread_num;
                    const md_backbone_angles_t* src_angles[2] = {
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[1],
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[2],
                    };
                    md_molecule_t& mol = data->state->mold.mol;
                    for (size_t i = range_beg; i < range_end; ++i) {
                        float phi[2] = {src_angles[0][i].phi, src_angles[1][i].phi};
                        float psi[2] = {src_angles[0][i].psi, src_angles[1][i].psi};

                        phi[1] = deperiodizef(phi[1], phi[0], (float)TWO_PI);
                        psi[1] = deperiodizef(psi[1], psi[0], (float)TWO_PI);

                        float final_phi = lerp(phi[0], phi[1], data->t);
                        float final_psi = lerp(psi[0], psi[1], data->t);
                        mol.protein_backbone.angle[i] = {deperiodizef(final_phi, 0, (float)TWO_PI), deperiodizef(final_psi, 0, (float)TWO_PI)};
                    }
                });

                tasks[num_tasks++] = angle_task;
                break;
            }
            case InterpolationMode::CubicSpline: {
                task_system::ID angle_task = task_system::create_pool_task(STR_LIT("## Interpolate Backbone Angles"), (uint32_t)mol.protein_backbone.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                    (void)thread_num;
                    const md_backbone_angles_t* src_angles[4] = {
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[0],
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[1],
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[2],
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[3],
                    };
                    md_molecule_t& mol = data->state->mold.mol;
                    for (size_t i = range_beg; i < range_end; ++i) {
                        float phi[4] = {src_angles[0][i].phi, src_angles[1][i].phi, src_angles[2][i].phi, src_angles[3][i].phi};
                        float psi[4] = {src_angles[0][i].psi, src_angles[1][i].psi, src_angles[2][i].psi, src_angles[3][i].psi};

                        phi[0] = deperiodizef(phi[0], phi[1], (float)TWO_PI);
                        phi[2] = deperiodizef(phi[2], phi[1], (float)TWO_PI);
                        phi[3] = deperiodizef(phi[3], phi[2], (float)TWO_PI);

                        psi[0] = deperiodizef(psi[0], psi[1], (float)TWO_PI);
                        psi[2] = deperiodizef(psi[2], psi[1], (float)TWO_PI);
                        psi[3] = deperiodizef(psi[3], psi[2], (float)TWO_PI);

                        float final_phi = cubic_spline(phi[0], phi[1], phi[2], phi[3], data->t, data->s);
                        float final_psi = cubic_spline(psi[0], psi[1], psi[2], psi[3], data->t, data->s);
                        mol.protein_backbone.angle[i] = {deperiodizef(final_phi, 0, (float)TWO_PI), deperiodizef(final_psi, 0, (float)TWO_PI)};
                    }
                });

                tasks[num_tasks++] = angle_task;
                break;
            }
            default:
                ASSERT(false);
                break;
        }
    }

    if (mol.protein_backbone.secondary_structure) {
        switch (mode) {
            case InterpolationMode::Nearest: {
                task_system::ID ss_task = task_system::create_pool_task(STR_LIT("## Interpolate Secondary Structures"), [data = &payload]() {
                    const md_secondary_structure_t* src_ss[2] = {
                        (md_secondary_structure_t*)data->state->trajectory_data.secondary_structure.data + data->state->trajectory_data.secondary_structure.stride * data->frames[1],
                        (md_secondary_structure_t*)data->state->trajectory_data.secondary_structure.data + data->state->trajectory_data.secondary_structure.stride * data->frames[2],
                    };
                    const md_secondary_structure_t* ss = data->t < 0.5f ? src_ss[0] : src_ss[1];
                    MEMCPY(data->state->mold.mol.protein_backbone.secondary_structure, ss, data->state->mold.mol.protein_backbone.count * sizeof(md_secondary_structure_t));
                });

                tasks[num_tasks++] = ss_task;
                break;
            }
            case InterpolationMode::Linear: {
                task_system::ID ss_task = task_system::create_pool_task(STR_LIT("## Interpolate Secondary Structures"), (uint32_t)mol.protein_backbone.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                    (void)thread_num;
                    const md_secondary_structure_t* src_ss[2] = {
                        (md_secondary_structure_t*)data->state->trajectory_data.secondary_structure.data + data->state->trajectory_data.secondary_structure.stride * data->frames[1],
                        (md_secondary_structure_t*)data->state->trajectory_data.secondary_structure.data + data->state->trajectory_data.secondary_structure.stride * data->frames[2],
                    };
                    for (size_t i = range_beg; i < range_end; ++i) {
                        const vec4_t ss_f[2] = {
                            convert_color((uint32_t)src_ss[0][i]),
                            convert_color((uint32_t)src_ss[1][i]),
                        };
                        const vec4_t ss_res = vec4_lerp(ss_f[0], ss_f[1], data->t);
                        data->state->mold.mol.protein_backbone.secondary_structure[i] = (md_secondary_structure_t)convert_color(ss_res);
                    }
                });

                tasks[num_tasks++] = ss_task;
                break;
            }
            case InterpolationMode::CubicSpline: {
                task_system::ID ss_task = task_system::create_pool_task(STR_LIT("## Interpolate Secondary Structures"), (uint32_t)mol.protein_backbone.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                    (void)thread_num;
                    const md_secondary_structure_t* src_ss[4] = {
                        (md_secondary_structure_t*)data->state->trajectory_data.secondary_structure.data + data->state->trajectory_data.secondary_structure.stride * data->frames[0],
                        (md_secondary_structure_t*)data->state->trajectory_data.secondary_structure.data + data->state->trajectory_data.secondary_structure.stride * data->frames[1],
                        (md_secondary_structure_t*)data->state->trajectory_data.secondary_structure.data + data->state->trajectory_data.secondary_structure.stride * data->frames[2],
                        (md_secondary_structure_t*)data->state->trajectory_data.secondary_structure.data + data->state->trajectory_data.secondary_structure.stride * data->frames[3],
                    };
                    for (size_t i = range_beg; i < range_end; ++i) {
                        const vec4_t ss_f[4] = {
                            convert_color((uint32_t)src_ss[0][i]),
                            convert_color((uint32_t)src_ss[1][i]),
                            convert_color((uint32_t)src_ss[2][i]),
                            convert_color((uint32_t)src_ss[3][i]),
                        };
                        const vec4_t ss_res = cubic_spline(ss_f[0], ss_f[1], ss_f[2], ss_f[3], data->t, data->s);
                        data->state->mold.mol.protein_backbone.secondary_structure[i] = (md_secondary_structure_t)convert_color(ss_res);
                    }
                });

                tasks[num_tasks++] = ss_task;
                break;
            }
            default:
                ASSERT(false);
                break;
        }
    }

    if (num_tasks > 0) {
        for (int i = 1; i < num_tasks; ++i) {
            task_system::set_task_dependency(tasks[i], tasks[i-1]);
        }
        task_system::enqueue_task(tasks[0]);
        task_system::task_wait_for(tasks[num_tasks - 1]);
    }

    vec3_t aabb_min = payload.aabb_min[0];
    vec3_t aabb_max = payload.aabb_max[0];
    for (size_t i = 1; i < task_system::pool_num_threads(); ++i) {
        aabb_min = vec3_min(aabb_min, payload.aabb_min[i]);
        aabb_max = vec3_max(aabb_max, payload.aabb_max[i]);
    }
    state->mold.mol_aabb_min = aabb_min;
    state->mold.mol_aabb_max = aabb_max;
    mol.unit_cell = payload.unit_cell;

    state->mold.dirty_buffers |= MolBit_DirtyPosition;
    state->mold.dirty_buffers |= MolBit_DirtySecondaryStructure;
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

    param.matrix.curr.view = camera_world_to_view_matrix(data->view.camera);
    param.matrix.inv.view  = camera_view_to_world_matrix(data->view.camera);

    if (data->visuals.temporal_aa.enabled && data->visuals.temporal_aa.jitter) {
        static uint32_t i = 0;
        i = (i+1) % (uint32_t)ARRAY_SIZE(data->view.jitter.sequence);
        param.jitter.curr = data->view.jitter.sequence[i] - 0.5f;
        if (data->view.mode == CameraMode::Perspective) {
            const vec2_t j = param.jitter.curr;
            const int w = data->gbuffer.width;
            const int h = data->gbuffer.height;
            param.matrix.curr.proj = camera_perspective_projection_matrix(data->view.camera, w, h, j.x, j.y);
            param.matrix.inv.proj  = camera_inverse_perspective_projection_matrix(data->view.camera, w, h, j.x, j.y);
        } else {
            const float aspect_ratio = (float)data->gbuffer.width / (float)data->gbuffer.height;
            const float h = data->view.camera.focus_distance * tanf(data->view.camera.fov_y * 0.5f);
            const float w = aspect_ratio * h;
            const vec2_t scl = {w / data->gbuffer.width * 2.0f, h / data->gbuffer.height * 2.0f};
            const vec2_t j = param.jitter.curr * scl;
            param.matrix.curr.proj = camera_orthographic_projection_matrix(-w + j.x, w + j.x, -h + j.y, h + j.y, data->view.camera.near_plane, data->view.camera.far_plane);
            param.matrix.inv.proj  = camera_inverse_orthographic_projection_matrix(-w + j.x, w + j.x, -h + j.y, h + j.y, data->view.camera.near_plane, data->view.camera.far_plane);
        }
    } else {
        param.jitter.curr = {0,0};
        if (data->view.mode == CameraMode::Perspective) {
            param.matrix.curr.proj = camera_perspective_projection_matrix(data->view.camera, (float)data->gbuffer.width / (float)data->gbuffer.height);
            param.matrix.inv.proj = camera_inverse_perspective_projection_matrix(data->view.camera, (float)data->gbuffer.width / (float)data->gbuffer.height);
        } else {
            const float aspect_ratio = (float)data->gbuffer.width / (float)data->gbuffer.height;
            const float h = data->view.camera.focus_distance * tanf(data->view.camera.fov_y * 0.5f);
            const float w = aspect_ratio * h;
            const float n = data->view.camera.near_plane;
            const float f = data->view.camera.far_plane;
            param.matrix.curr.proj = camera_orthographic_projection_matrix(-w, w, -h, h, n, f);
            param.matrix.inv.proj = camera_inverse_orthographic_projection_matrix(-w, w, -h, h, n, f);
        }
    }

    param.matrix.curr.norm = mat4_transpose(param.matrix.inv.view);
}

static void reset_view(ApplicationState* data, const md_bitfield_t* target, bool move_camera, bool smooth_transition) {
    ASSERT(data);
    if (!data->mold.mol.atom.count) return;
    const auto& mol = data->mold.mol;

    size_t popcount = 0;
    if (target) {
        popcount = md_bitfield_popcount(target);
    }
    vec3_t aabb_min = {};
    vec3_t aabb_max = {};
    
    if (0 < popcount && popcount < mol.atom.count) {
        md_vm_arena_temp_t tmp = md_vm_arena_temp_begin(frame_alloc);
        int32_t* indices = (int32_t*)md_vm_arena_push_array(frame_alloc, int32_t, popcount);
        size_t len = md_bitfield_iter_extract_indices(indices, popcount, md_bitfield_iter_create(target));
        if (len > popcount || len > mol.atom.count) {
            MD_LOG_DEBUG("Error: Invalid number of indices");
            len = MIN(popcount, mol.atom.count);
        }
        md_util_aabb_compute(aabb_min.elem, aabb_max.elem, mol.atom.x, mol.atom.y, mol.atom.z, nullptr, indices, len);
        md_vm_arena_temp_end(tmp);
    } else {
        md_util_aabb_compute(aabb_min.elem, aabb_max.elem, mol.atom.x, mol.atom.y, mol.atom.z, nullptr, nullptr, mol.atom.count);
    }

    vec3_t optimal_pos;
    quat_t optimal_ori;
    float  optimal_dist;
    camera_compute_optimal_view(&optimal_pos, &optimal_ori, &optimal_dist, aabb_min, aabb_max);

    if (move_camera) {
        data->view.animation.target_position    = optimal_pos;
        data->view.animation.target_orientation = optimal_ori;
        data->view.animation.target_distance    = optimal_dist;

        if (!smooth_transition) {
            data->view.camera.position       = optimal_pos;
            data->view.camera.orientation    = optimal_ori;
            data->view.camera.focus_distance = optimal_dist;
        }
    }

    float max_cell_ext = vec3_reduce_max(mat3_diag(data->mold.mol.unit_cell.basis));
    float max_aabb_ext = vec3_reduce_max(vec3_sub(aabb_max, aabb_min));

    data->view.camera.near_plane = 1.0f;
    data->view.camera.far_plane = 10000.0f;
    data->view.trackball_param.max_distance = MAX(max_cell_ext, max_aabb_ext) * 10.0f;
}

static bool filter_expression(ApplicationState* data, str_t expr, md_bitfield_t* mask, bool* is_dynamic = NULL, char* error_buf = NULL, int error_cap = 0) {
    if (data->mold.mol.atom.count == 0) return false;    
    
    bool success = false;

    if (md_semaphore_aquire(&data->script.ir_semaphore)) {
        defer { md_semaphore_release(&data->script.ir_semaphore); };
        return md_filter(mask, expr, &data->mold.mol, data->script.ir, is_dynamic, error_buf, error_cap);
    }

    return success;
}

// ### DRAW WINDOWS ###
static void draw_main_menu(ApplicationState* data) {
    ASSERT(data);
    bool new_clicked = false;
    char path_buf[2048] = "";

    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("Load File", "CTRL+L")) {
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Open)) {
                    file_queue_push(&data->file_queue, str_from_cstr(path_buf), FileFlags_ShowDialogue);
                }
            }
            if (ImGui::MenuItem("Open Workspace", "CTRL+O")) {
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
            if (ImGui::MenuItem("Save As")) {
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
                reset_view(data, &data->representation.visibility_mask, true, true);
            }
            ImGui::Separator();
            ImGui::Checkbox("Vsync", &data->app.window.vsync);
            ImGui::Separator();

            ImGui::BeginGroup();
            ImGui::Text("Camera");
            {
                ImGui::Combo("Mode", (int*)(&data->view.mode), "Perspective\0Orthographic\0");
                if (data->view.mode == CameraMode::Perspective) {
                    float fov = RAD_TO_DEG(data->view.camera.fov_y);
                    if (ImGui::SliderFloat("field of view", &fov, 12.5f, 80.0f)) {
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
            ImGui::EndGroup();
            ImGui::Separator();
            ImGui::Checkbox("FXAA", &data->visuals.fxaa.enabled);
            // Temporal
            ImGui::BeginGroup();
            {
                ImGui::Checkbox("Temporal AA", &data->visuals.temporal_aa.enabled);
                if (data->visuals.temporal_aa.enabled) {
                    // ImGui::Checkbox("Jitter Samples", &data->visuals.temporal_reprojection.jitter);
                    ImGui::SliderFloat("Feedback Min", &data->visuals.temporal_aa.feedback_min, 0.5f, 1.0f);
                    ImGui::SliderFloat("Feedback Max", &data->visuals.temporal_aa.feedback_max, 0.5f, 1.0f);
                    ImGui::Checkbox("Motion Blur", &data->visuals.temporal_aa.motion_blur.enabled);
                    if (data->visuals.temporal_aa.motion_blur.enabled) {
                        ImGui::SliderFloat("Motion Scale", &data->visuals.temporal_aa.motion_blur.motion_scale, 0.f, 2.0f);
                    }
                    ImGui::Checkbox("Sharpen", &data->visuals.sharpen.enabled);
                    if (data->visuals.sharpen.enabled) {
                        ImGui::SliderFloat("Weight", &data->visuals.sharpen.weight, 0.0f, 4.0f);
                    }
                }
            }
            ImGui::EndGroup();
            ImGui::Separator();

            // SSAO
            ImGui::BeginGroup();
            ImGui::PushID("SSAO");
            ImGui::Checkbox("SSAO", &data->visuals.ssao.enabled);
            if (data->visuals.ssao.enabled) {
                ImGui::SliderFloat("Intensity", &data->visuals.ssao.intensity, 0.5f, 12.f);
                ImGui::SliderFloat("Radius", &data->visuals.ssao.radius, 1.f, 30.f);
                ImGui::SliderFloat("Bias", &data->visuals.ssao.bias, 0.0f, 1.0f);
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
                ImGui::SliderFloat("Focus Scale", &data->visuals.dof.focus_scale, 0.001f, 100.f);
            }
            ImGui::EndGroup();
            ImGui::Separator();

            // Tonemapping
            ImGui::BeginGroup();
            ImGui::Checkbox("Tonemapping", &data->visuals.tonemapping.enabled);
            if (data->visuals.tonemapping.enabled) {
                // ImGui::Combo("Function", &data->visuals.tonemapping.tonemapper, "Passthrough\0Exposure Gamma\0Filmic\0\0");
                ImGui::SliderFloat("Exposure", &data->visuals.tonemapping.exposure, 0.01f, 10.f);
                ImGui::SliderFloat("Gamma", &data->visuals.tonemapping.gamma, 1.0f, 3.0f);
            }
            ImGui::EndGroup();
            ImGui::Separator();

            ImGui::BeginGroup();
            ImGui::Checkbox("Simulation Box", &data->simulation_box.enabled);
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
            ImGui::Checkbox("Timelines", &data->timeline.show_window);
            ImGui::Checkbox("Distributions", &data->distributions.show_window);
            ImGui::Checkbox("Density Volumes", &data->density_volume.show_window);
            ImGui::Checkbox("Dataset", &data->dataset.show_window);

            viamd::event_system_broadcast_event(viamd::EventType_ViamdDrawMenu);

            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Selection")) {
            ImGui::Combo("Granularity", (int*)(&data->selection.granularity), "Atom\0Residue\0Chain\0\0");
            int64_t num_selected_atoms = md_bitfield_popcount(&data->selection.selection_mask);
            if (ImGui::MenuItem("Invert")) {
                md_bitfield_not_inplace(&data->selection.selection_mask, 0, data->mold.mol.atom.count);
            }
            if (ImGui::IsItemHovered()) {
                md_bitfield_not(&data->selection.highlight_mask, &data->selection.selection_mask, 0, data->mold.mol.atom.count);
            }
            if (ImGui::MenuItem("Query")) data->selection.query.show_window = true;
            if (num_selected_atoms == 0) ImGui::PushDisabled();
            if (ImGui::MenuItem("Grow"))  data->selection.grow.show_window = true;
            if (num_selected_atoms == 0) ImGui::PopDisabled();
            if (ImGui::MenuItem("Clear")) {
                md_bitfield_clear(&data->selection.selection_mask);
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
                        update_all_representations(data);
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
                        update_all_representations(data);
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
                    }
                    else {
                        LOG_ERROR("Supplied image extension is not supported");
                    }
                }
                ImGui::GetCurrentWindow()->Hidden = true;
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Operations")) {
            ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg;
            if (ImGui::BeginTable("##table", 2, flags)) {
                /*
                ImGui::TableSetupColumn("Once", 0);
                ImGui::TableSetupColumn("Always", 0);
                ImGui::TableHeadersRow();
                */
                const float button_width = (ImGui::GetFontSize() / 20.f) * 150.f;

                bool do_pbc = false;
                bool do_unwrap = false;

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                if (ImGui::Button("PBC", ImVec2(button_width,0))) {
                    do_pbc = true;
                }
                ImGui::SetItemTooltip("Enforce Periodic Boundary Conditions (Once)");

                ImGui::TableSetColumnIndex(1);
                if (ImGui::Checkbox("##pbc", &data->operations.apply_pbc) && data->operations.apply_pbc) {
                    do_pbc = true;
                }
                ImGui::SetItemTooltip("Enforce Periodic Boundary Conditions (Always)");

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                if (ImGui::Button("Unwrap", ImVec2(button_width, 0))) {
                    do_unwrap = true;
                }
                ImGui::SetItemTooltip("Unwrap structures present in the system (Once)");

                ImGui::TableSetColumnIndex(1);
                if (ImGui::Checkbox("##unwrap", &data->operations.unwrap_structures) && data->operations.unwrap_structures) {
                    do_unwrap = true;
                }
                ImGui::SetItemTooltip("Unwrap structures present in the system (Always)");

                if (do_pbc) {
                    md_molecule_t& mol = data->mold.mol;
                    md_util_pbc(mol.atom.x, mol.atom.y, mol.atom.z, 0, mol.atom.count, &mol.unit_cell);
                    data->mold.dirty_buffers |= MolBit_DirtyPosition;
                }

                if (do_unwrap) {
                    md_molecule_t& mol = data->mold.mol;
                    size_t num_structures = md_index_data_num_ranges(mol.structure);
                    for (size_t i = 0; i < num_structures; ++i) {
                        const int32_t* s_idx = md_index_range_beg(mol.structure, i);
                        const size_t   s_len = md_index_range_size(mol.structure, i);
                        md_util_unwrap(mol.atom.x, mol.atom.y, mol.atom.z, s_idx, s_len, &mol.unit_cell);
                    }
                    data->mold.dirty_buffers |= MolBit_DirtyPosition;
                }

                ImGui::EndTable();
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Settings")) {
            ImGui::Checkbox("Prefetch Frames", &data->settings.prefetch_frames);
            ImGui::SetItemTooltip("Prefetch frames during animation\n");
            ImGui::Checkbox("Keep Representations", &data->settings.keep_representations);
            ImGui::SetItemTooltip("Keep representations when loading new topology (Does not apply for workspaces)\n");

            // Font
            ImFont* font_current = ImGui::GetFont();
            if (ImGui::BeginCombo("Font", font_current->GetDebugName()))
            {
                ImGuiIO& io = ImGui::GetIO();
                for (int n = 0; n < io.Fonts->Fonts.Size; n++) {
                    ImFont* font = io.Fonts->Fonts[n];
                    ImGui::PushID((void*)font);
                    if (ImGui::Selectable(font->GetDebugName(), font == font_current))
                        io.FontDefault = font;
                    ImGui::PopID();
                }
                ImGui::EndCombo();
            }

            /*
            ImGui::Text("Units");
            char buf[64];

            unit_print_long(buf, sizeof(buf), {.base = { .length = data->mold.unit_base.length }, .dim = { .length = 1 }});
            if (ImGui::BeginCombo("length", buf)) {
                for (uint32_t i = 0; i < UNIT_LENGTH_COUNT; ++i) {
                    unit_print_long(buf, sizeof(buf), {.base = { .length = i }, .dim = { .length = 1 }});
                    ImGui::PushID((int)i);
                    if (ImGui::Selectable(buf, i == data->mold.unit_base.length))
                        data->mold.unit_base.length = i;
                    ImGui::PopID();
                }
                ImGui::EndCombo();
            }

            unit_print_long(buf, sizeof(buf), {.base = { .time = data->mold.unit_base.time }, .dim = {.time = 1 }});
            if (ImGui::BeginCombo("time", buf)) {
                for (uint32_t i = 0; i < UNIT_TIME_COUNT; ++i) {
                    unit_print_long(buf, sizeof(buf), {.base = { .time = i }, .dim = {.time = 1}});
                    ImGui::PushID((int)i);
                    if (ImGui::Selectable(buf, i == data->mold.unit_base.time))
                        data->mold.unit_base.time = i;
                    ImGui::PopID();
                }
                ImGui::EndCombo();
            }

            unit_print_long(buf, sizeof(buf), {.base = { .angle = data->mold.unit_base.angle }, .dim = {.angle = 1 }});
            if (ImGui::BeginCombo("angle", buf)) {
                for (uint32_t i = 0; i < UNIT_ANGLE_COUNT; ++i) {
                    unit_print_long(buf, sizeof(buf), {.base = { .angle = i }, .dim = {.angle = 1 }});
                    ImGui::PushID((int)i);
                    if (ImGui::Selectable(buf, i == data->mold.unit_base.angle))
                        data->mold.unit_base.angle = i;
                    ImGui::PopID();
                }
                ImGui::EndCombo();
            }
            */

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
        bool path_invalid = !state.path_is_valid && state.path_buf[0] != '\0';
        const int    loader_count = (int)load::loader_count();
        const str_t* loader_ext_str = load::loader_extensions();
        const str_t* loader_name_str = load::loader_names();

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
            state.loader_idx = -1;
            str_t ext;
            if (extract_ext(&ext, path)) {
                for (int i = 0; i < loader_count; ++i) {
                    if (str_eq_ignore_case(ext, loader_ext_str[i])) {
                        state.loader_idx = i;
                        break;
                    }
                }
            }
        }


        if (ImGui::BeginCombo("Loader", state.loader_idx > -1 ? loader_name_str[state.loader_idx].ptr : "")) {
            for (int i = 0; i < loader_count; ++i) {
                const char* str = loader_name_str[i].ptr;
                if (!str) continue;
                if (ImGui::Selectable(str, state.loader_idx == i)) {
                    state.loader_idx = i;
                }
            }
            ImGui::EndCombo();
        }

        str_t cur_ext = {};
        if (state.loader_idx > -1) {
            cur_ext = loader_ext_str[state.loader_idx];
        }

        // True if the button should be enabled
        bool load_enabled = (state.path_is_valid && state.loader_idx > -1);

        // Draw Options
        md_molecule_loader_i* mol_loader = load::mol::loader_from_ext(cur_ext);
        bool show_cg = state.path_is_valid && mol_loader;
        if (show_cg) {
            ImGui::Checkbox("Coarse Grained", &state.coarse_grained);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Enable if the data should be interpreted as coarse grained particles");
            }
        }

        bool show_lammps_atom_format = state.path_is_valid && mol_loader && (mol_loader == md_lammps_molecule_api());
        if (show_lammps_atom_format) {
            const char** atom_format_names = md_lammps_atom_format_names();
            const char** atom_format_strings = md_lammps_atom_format_strings();

            if (state.atom_format_idx == -1) {
                // Try to determine format from file
                md_lammps_atom_format_t format = md_lammps_atom_format_from_file(path);
                state.atom_format_idx = format;
                strncpy(state.atom_format_buf, atom_format_strings[format], sizeof(state.atom_format_buf));
            }
            ASSERT(state.atom_format_idx > -1);

            if (ImGui::BeginCombo("Atom Format", state.atom_format_idx > 0 ? atom_format_names[state.atom_format_idx] : "user defined")) {
                for (int i = 0; i < MD_LAMMPS_ATOM_FORMAT_COUNT; ++i) {
                    if (ImGui::Selectable(i > 0 ? atom_format_names[i] : "user defined", state.atom_format_idx == i)) {
                        state.atom_format_idx = i;
                        int source_idx = i > 0 ? i : MD_LAMMPS_ATOM_FORMAT_FULL;
                        strncpy(state.atom_format_buf, atom_format_strings[source_idx], sizeof(state.atom_format_buf));
                    }
                }
                ImGui::EndCombo();
            }

            if (state.atom_format_idx == MD_LAMMPS_ATOM_FORMAT_UNKNOWN) {
                bool valid = state.atom_format_valid;
                if (!valid) ImGui::PushInvalid();
                if (ImGui::InputText("##atom_format", state.atom_format_buf, sizeof(state.atom_format_buf))) {
                    state.atom_format_valid = md_lammps_validate_atom_format(state.atom_format_buf, state.err_buf, sizeof(state.err_buf));
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
                load_enabled = false;
            }
        }

        md_trajectory_loader_i* traj_loader = load::traj::loader_from_ext(cur_ext);

        enum Action {
            Action_None,
            Action_Cancel,
            Action_Load,
        };
        Action action = Action_None;

        if (!load_enabled) ImGui::PushDisabled();
        if (ImGui::Button("Load")) {
            action = Action_Load;
        }
        if (!load_enabled) ImGui::PopDisabled();

        ImGui::SameLine();
        if (ImGui::Button("Cancel")) {
            action = Action_Cancel;
        }

        switch (action) {
        case Action_Load:
        {
            LoadParam param = {};
            param.file_path = path;
            param.mol_loader = mol_loader;
            param.traj_loader = traj_loader;
            param.coarse_grained = state.coarse_grained;

            md_lammps_molecule_loader_arg_t lammps_arg = {};
            if (mol_loader == md_lammps_molecule_api()) {
                lammps_arg = md_lammps_molecule_loader_arg(state.atom_format_buf);
                param.mol_loader_arg = &lammps_arg;
            }

            if (load_dataset_from_file(data, param)) {
                if (mol_loader && !data->settings.keep_representations) {
                    clear_representations(data);
                    create_default_representations(data);
                }
                data->animation = {};
                reset_view(data, &data->representation.visibility_mask, true, true);
            }
        }
            [[fallthrough]];
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
    if (data->mold.mol.atom.count == 0 || !data->mold.mol.atom.element) {
        return;
    }

    for (size_t j = 0; j < md_array_size(data->dataset.atom_element_remappings); ++j) {
        str_t lbl = str_from_cstr(data->dataset.atom_element_remappings[j].lbl);
        md_element_t elem = data->dataset.atom_element_remappings[j].elem;
        float radius = md_util_element_vdw_radius(elem);
        float mass = md_util_element_atomic_mass(elem);

        for (size_t i = 0; i < data->mold.mol.atom.count; ++i) {
            if (str_eq(lbl, data->mold.mol.atom.type[i])) {
                data->mold.mol.atom.element[i] = elem;
                data->mold.mol.atom.radius[i] = radius;
                data->mold.mol.atom.mass[i] = mass;
                data->mold.dirty_buffers |= MolBit_DirtyRadius;
            }
        }
    }
    md_molecule_t* mol = &data->mold.mol;
    
    
    md_array_free(mol->bond.pairs, data->mold.mol_alloc);
    md_array_free(mol->bond.order, data->mold.mol_alloc);

    md_array_free(mol->bond.conn.atom_idx, data->mold.mol_alloc);
    md_array_free(mol->bond.conn.bond_idx, data->mold.mol_alloc);

    md_index_data_free(&mol->structure);
    md_index_data_free(&mol->ring);
    
    md_util_molecule_postprocess(mol, data->mold.mol_alloc, MD_UTIL_POSTPROCESS_BOND_BIT | MD_UTIL_POSTPROCESS_STRUCTURE_BIT);
    data->mold.dirty_buffers |= MolBit_DirtyBonds;

    update_all_representations(data);
}

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

    size_t temp_pos = md_temp_get_pos();
    defer { md_temp_set_pos_back(temp_pos); };

    md_array(md_range_t) items = 0;
    
    for (size_t i = 0; i < num_indices; ++i) {
        int idx = indices[i];
        
        if (idx - prev_idx > 1) {
            if (prev_idx > range_beg) {
                md_range_t item = {range_beg - ref_idx + 1, prev_idx - ref_idx + 1};
                md_array_push(items, item, md_get_temp_allocator());
            } else if (prev_idx != -1) {
                md_range_t item = {prev_idx - ref_idx + 1, prev_idx - ref_idx + 1};
                md_array_push(items, item, md_get_temp_allocator());
            }
            range_beg = idx;
        }

        prev_idx = idx;
    }

    if (prev_idx - range_beg > 0) {
        md_range_t item = {range_beg - ref_idx + 1, prev_idx - ref_idx + 1};
        md_array_push(items, item, md_get_temp_allocator());
    } else if (prev_idx != -1) {
        md_range_t item = {prev_idx - ref_idx + 1, prev_idx - ref_idx + 1};
        md_array_push(items, item, md_get_temp_allocator());
    }

    const int64_t num_items = (int64_t)md_array_size(items);
    if (num_items > 1) sb += "{";
    for (int64_t i = 0; i < num_items; ++i) {
        md_range_t item = items[i];
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
}

static void write_script_atom_ranges(md_strb_t* sb, const md_bitfield_t* bf, int ref_idx = 0) {
    ASSERT(sb);
    ASSERT(bf);

    const int64_t popcount = md_bitfield_popcount(bf);

    md_strb_fmt(sb, "atom(");

    if (popcount > 1) {
        md_strb_push_char(sb, '{');
    }

    int range_beg = -1;
    int prev_idx  = -1;

    uint64_t beg_idx = bf->beg_bit;
    uint64_t end_idx = bf->end_bit;
    while ((beg_idx = md_bitfield_scan(bf, beg_idx, end_idx)) != 0) {
        int idx = (int)beg_idx - 1;
        if (range_beg == -1) range_beg = idx;

        if (idx - prev_idx > 1) {
            if (prev_idx > range_beg) {
                md_strb_fmt(sb, "%i:%i,", range_beg - ref_idx + 1, prev_idx - ref_idx + 1);
            } else if (prev_idx != -1) {
                md_strb_fmt(sb, "%i,", prev_idx - ref_idx + 1);
            }
            range_beg = idx;
        }

        prev_idx = idx;
    }

    md_strb_pop(sb, 1);
    if (prev_idx - range_beg > 0) {
        md_strb_fmt(sb, "%i:%i", range_beg - ref_idx + 1, prev_idx - ref_idx + 1);
    } else if (prev_idx != -1) {
        md_strb_fmt(sb, "%i", prev_idx - ref_idx + 1);
    }

    if (popcount > 1) {
        md_strb_push_char(sb, '}');
    }
    md_strb_push_char(sb, ')');
}

static md_array(str_t) generate_script_selection_suggestions(str_t ident, const md_bitfield_t* bf, const md_molecule_t* mol) {
    md_array(str_t) suggestions = 0;

    bool within_same_residue = true;
    bool within_same_chain   = true;

    int32_t     chain_idx = -1;
    int32_t     res_idx   = -1;

    md_bitfield_iter_t it = md_bitfield_iter_create(bf);
    while (md_bitfield_iter_next(&it)) {
        uint64_t i = md_bitfield_iter_idx(&it);
        if (mol->atom.chainid) {
            if (chain_idx == -1) {
                chain_idx = mol->atom.chain_idx[i];
            } else if (chain_idx != mol->atom.chain_idx[i]) {
                within_same_chain = false;
            }
        }

        if (mol->atom.res_idx) {
            if (res_idx == -1) {
                res_idx = mol->atom.res_idx[i];
            } else if (mol->atom.res_idx[i] != res_idx) {
                within_same_residue = false;                        
            }
        } 

        if (!within_same_residue && !within_same_chain) {
            break;
        }
    }
    
    const int64_t popcount = md_bitfield_popcount(bf);
    
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

    if (res_idx != -1 && within_same_residue) {
        const md_range_t range = md_residue_atom_range(mol->residue, res_idx);
        if (popcount != range.end - range.beg) {
            md_strb_reset(&sb);
            sb += ident;
            sb += " = ";

            // Subset of residue is selected
            write_atom_remainder(sb, bf, range.beg);
            if (md_strb_len(sb) < 512) {
                str_t resname = LBL_TO_STR(mol->residue.name[res_idx]);
                md_strb_fmt(&sb, " in resname(\"%.*s\");", (int)resname.len, resname.ptr);
                md_array_push(suggestions, str_copy(md_strb_to_str(sb), frame_alloc), frame_alloc);
            }
        }
    }

    else if (chain_idx != -1 && within_same_chain) {
        const md_range_t range = md_chain_atom_range(mol->chain, chain_idx);
        if (popcount != range.end - range.beg) {
            md_strb_reset(&sb);
            sb += ident;
            sb += " = ";

            // Subset of chain is selected
            write_atom_remainder(sb, bf, range.beg);
            if (md_strb_len(sb) < 512) {
                str_t chain_id = LBL_TO_STR(mol->chain.id[chain_idx]);
                md_strb_fmt(&sb, " in chain(\"%.*s\");", (int)chain_id.len, chain_id.ptr);
                md_array_push(suggestions, str_copy(md_strb_to_str(sb), frame_alloc), frame_alloc);
            }
        }
    }

    md_bitfield_t tmp_bf = md_bitfield_create(frame_alloc);
    md_bitfield_copy(&tmp_bf, bf);
    
    md_array(int) complete_chains = 0;
    md_array(int) complete_residues = 0;
    
    if (mol->chain.count) {
        for (size_t i = 0; i < mol->chain.count; ++i) {    
            const md_range_t range = md_chain_atom_range(mol->chain, i);
            if (md_bitfield_test_all_range(&tmp_bf, range.beg, range.end)) {
                md_array_push(complete_chains, (int)i, frame_alloc);
                md_bitfield_clear_range(&tmp_bf, range.beg, range.end);
            }
        }
    }

    if (mol->residue.count) {
        for (size_t i = 0; i < mol->residue.count; ++i) {    
            const md_range_t range = md_residue_atom_range(mol->residue, i);
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
            md_range_t range = md_chain_residue_range(mol->chain, complete_chains[i]);
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
void draw_context_popup(ApplicationState* data) {
    ASSERT(data);

    if (!data->mold.mol.atom.count) return;

    const int64_t sss_count = single_selection_sequence_count(&data->selection.single_selection_sequence);
    const int64_t num_frames = md_trajectory_num_frames(data->mold.traj);
    const int64_t num_atoms_selected = md_bitfield_popcount(&data->selection.selection_mask);

#if 0
    // FOR DEBUGGING
    if (ImGui::Begin("SSS windows")) {
        ImGui::Text("sel_seq = %i,%i,%i,%i",
            data->selection.single_selection_sequence.idx[0],
            data->selection.single_selection_sequence.idx[1],
            data->selection.single_selection_sequence.idx[2],
            data->selection.single_selection_sequence.idx[3]);
        ImGui::End();
    }
#endif

    if (ImGui::BeginPopup("AtomContextPopup")) {
        if (ImGui::BeginMenu("Script")) {
            bool any_suggestions = false;
            if (num_atoms_selected <= 4 && sss_count > 1) {
                any_suggestions = true;
                char buf[256] = "";
                if (sss_count == 2) {
                    int32_t idx[2] = {data->selection.single_selection_sequence.idx[0], data->selection.single_selection_sequence.idx[1]};
                    str_t ident = create_unique_identifier(data->script.ir, STR_LIT("dist"), frame_alloc);

                    snprintf(buf, sizeof(buf), "%.*s = distance(%i, %i);", (int)ident.len, ident.ptr, idx[0]+1, idx[1]+1);
                    if (ImGui::MenuItem(buf)) {
                        editor.AppendText("\n");
                        editor.AppendText(buf);
                        ImGui::CloseCurrentPopup();
                    }

                    if (data->mold.mol.residue.count) {
                        if (data->mold.mol.atom.res_idx[idx[0]] == data->mold.mol.atom.res_idx[idx[1]]) {
                            const int32_t res_idx = data->mold.mol.atom.res_idx[idx[0]];
                            const md_range_t range = md_residue_atom_range(data->mold.mol.residue, res_idx);
                            idx[0] -= range.beg;
                            idx[1] -= range.beg;

                            snprintf(buf, sizeof(buf), "%.*s = distance(%i, %i) in residue(%i);", (int)ident.len, ident.ptr, idx[0]+1, idx[1]+1, res_idx+1);
                            if (ImGui::MenuItem(buf)) {
                                editor.AppendText("\n");
                                editor.AppendText(buf);
                                ImGui::CloseCurrentPopup();
                            }

                            int32_t resid = data->mold.mol.residue.id[idx[0]];
                            snprintf(buf, sizeof(buf), "%.*s = distance(%i, %i) in resid(%i);", (int)ident.len, ident.ptr, idx[0]+1, idx[1]+1, resid);
                            if (ImGui::MenuItem(buf)) {
                                editor.AppendText("\n");
                                editor.AppendText(buf);
                                ImGui::CloseCurrentPopup();
                            }

                            str_t resname = data->mold.mol.residue.name[res_idx];
                            if (resname) {
                                snprintf(buf, sizeof(buf), "%.*s = distance(%i, %i) in resname(\"%s\");", (int)ident.len, ident.ptr, idx[0] + 1, idx[1] + 1, resname.ptr);
                                if (ImGui::MenuItem(buf)) {
                                    editor.AppendText("\n");
                                    editor.AppendText(buf);
                                    ImGui::CloseCurrentPopup();
                                }
                            }
                        }
                    }
                }
                else if(sss_count == 3) {
                    int32_t idx[3] = {data->selection.single_selection_sequence.idx[0], data->selection.single_selection_sequence.idx[1], data->selection.single_selection_sequence.idx[2]};
                    str_t ident = create_unique_identifier(data->script.ir, STR_LIT("ang"), frame_alloc);

                    snprintf(buf, sizeof(buf), "%.*s = angle(%i, %i, %i);", (int)ident.len, ident.ptr, idx[0]+1, idx[1]+1, idx[2]+1);
                    if (ImGui::MenuItem(buf)) {
                        editor.AppendText("\n");
                        editor.AppendText(buf);
                        ImGui::CloseCurrentPopup();
                    }

                    if (data->mold.mol.residue.count) {
                        if (data->mold.mol.atom.res_idx[idx[0]] == data->mold.mol.atom.res_idx[idx[1]] &&
                            data->mold.mol.atom.res_idx[idx[0]] == data->mold.mol.atom.res_idx[idx[2]]) {

                            const int32_t res_idx = data->mold.mol.atom.res_idx[idx[0]];
                            const md_range_t range = md_residue_atom_range(data->mold.mol.residue, res_idx);
                            idx[0] -= range.beg;
                            idx[1] -= range.beg;
                            idx[2] -= range.beg;

                            snprintf(buf, sizeof(buf), "%.*s = angle(%i, %i, %i) in residue(%i);", (int)ident.len, ident.ptr, idx[0]+1, idx[1]+1, idx[2]+1, res_idx+1);
                            if (ImGui::MenuItem(buf)) {
                                editor.AppendText("\n");
                                editor.AppendText(buf);
                                ImGui::CloseCurrentPopup();
                            }

                            int32_t resid = data->mold.mol.residue.id[idx[0]];
                            snprintf(buf, sizeof(buf), "%.*s = angle(%i, %i, %i) in resid(%i);", (int)ident.len, ident.ptr, idx[0]+1, idx[1]+1, idx[2]+1, resid);
                            if (ImGui::MenuItem(buf)) {
                                editor.AppendText("\n");
                                editor.AppendText(buf);
                                ImGui::CloseCurrentPopup();
                            }

                            str_t resname = data->mold.mol.residue.name[res_idx];
                            if (resname) {
                                snprintf(buf, sizeof(buf), "%.*s = angle(%i, %i, %i) in resname(\"%.*s\");", (int)ident.len, ident.ptr, idx[0]+1, idx[1]+1, idx[2]+1, (int)resname.len, resname.ptr);
                                if (ImGui::MenuItem(buf)) {
                                    editor.AppendText("\n");
                                    editor.AppendText(buf);
                                    ImGui::CloseCurrentPopup();
                                }
                            }
                        }
                    }
                }
                else if(sss_count == 4) {
                    int32_t idx[4] = {data->selection.single_selection_sequence.idx[0], data->selection.single_selection_sequence.idx[1], data->selection.single_selection_sequence.idx[2], data->selection.single_selection_sequence.idx[3]};
                    str_t ident = create_unique_identifier(data->script.ir, STR_LIT("dih"), frame_alloc);

                    snprintf(buf, sizeof(buf), "%.*s = dihedral(%i, %i, %i, %i);", (int)ident.len, ident.ptr, idx[0]+1, idx[1]+1, idx[2]+1, idx[3]+1);
                    if (ImGui::MenuItem(buf)) {
                        editor.AppendText("\n");
                        editor.AppendText(buf);
                        ImGui::CloseCurrentPopup();
                    }

                    if (data->mold.mol.residue.count) {
                        if (data->mold.mol.atom.res_idx[idx[0]] == data->mold.mol.atom.res_idx[idx[1]] &&
                            data->mold.mol.atom.res_idx[idx[0]] == data->mold.mol.atom.res_idx[idx[2]] &&
                            data->mold.mol.atom.res_idx[idx[0]] == data->mold.mol.atom.res_idx[idx[3]]) {

                            const int32_t res_idx = data->mold.mol.atom.res_idx[idx[0]];
                            const md_range_t range = md_residue_atom_range(data->mold.mol.residue, res_idx);
                            idx[0] -= range.beg;
                            idx[1] -= range.beg;
                            idx[2] -= range.beg;
                            idx[3] -= range.beg;

                            snprintf(buf, sizeof(buf), "%.*s = dihedral(%i, %i, %i, %i) in residue(%i);", (int)ident.len, ident.ptr, idx[0]+1, idx[1]+1, idx[2]+1, idx[3]+1, res_idx+1);
                            if (ImGui::MenuItem(buf)) {
                                editor.AppendText("\n");
                                editor.AppendText(buf);
                                ImGui::CloseCurrentPopup();
                            }

                            int32_t resid = data->mold.mol.residue.id[idx[0]];
                            snprintf(buf, sizeof(buf), "%.*s = dihedral(%i, %i, %i, %i) in resid(%i);", (int)ident.len, ident.ptr, idx[0]+1, idx[1]+1, idx[2]+1, idx[3]+1, resid);
                            if (ImGui::MenuItem(buf)) {
                                editor.AppendText("\n");
                                editor.AppendText(buf);
                                ImGui::CloseCurrentPopup();
                            }

                            str_t resname = data->mold.mol.residue.name[res_idx];
                            if (resname) {
                                snprintf(buf, sizeof(buf), "%.*s = dihedral(%i, %i, %i, %i) in resname(\"%.*s\");", (int)ident.len, ident.ptr, idx[0]+1, idx[1]+1, idx[2]+1, idx[3]+1, (int)resname.len, resname.ptr);
                                if (ImGui::MenuItem(buf)) {
                                    editor.AppendText("\n");
                                    editor.AppendText(buf);
                                    ImGui::CloseCurrentPopup();
                                }
                            }
                        }
                    }
                }
            }
            if (num_atoms_selected >= 1) {
                const md_bitfield_t* bf = &data->selection.selection_mask;
                str_t ident = create_unique_identifier(data->script.ir, STR_LIT("sel"), frame_alloc);
                
                md_array(str_t) suggestions = generate_script_selection_suggestions(ident, bf, &data->mold.mol);

                for (size_t i = 0; i < md_array_size(suggestions); ++i) {
                    str_t s = suggestions[i];
                    if (ImGui::MenuItem(s.ptr)) {
                        editor.AppendText("\n");
                        editor.AppendText(s.ptr);
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

        if (data->selection.atom_idx.right_click != -1 && data->mold.mol.atom.element) {
            int idx = data->selection.atom_idx.right_click;
            if (0 <= idx && idx < (int)data->mold.mol.atom.count) {
                char label[64] = "";
                str_t type = data->mold.mol.atom.type[idx];
                snprintf(label, sizeof(label), "Remap Element for '%.*s'", (int)type.len, type.ptr);
                if (ImGui::BeginMenu(label)) {
                    static char input_buf[32] = "";
                    md_element_t elem = data->mold.mol.atom.element[idx];
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
        if (data->selection.atom_idx.right_click != -1 && num_frames > 0) {
            if (ImGui::BeginMenu("Recenter Trajectory...")) {
                const int idx = data->selection.atom_idx.right_click;

                md_bitfield_t mask = {0};
                md_bitfield_init(&mask, frame_alloc);
                bool apply = false;

                apply |= ImGui::MenuItem("on Atom");
                if (ImGui::IsItemHovered()) {
                    md_bitfield_set_bit(&mask, idx);
                }

                if (data->mold.mol.residue.count > 0 && data->mold.mol.atom.res_idx && data->mold.mol.atom.res_idx[idx] != -1) {
                    apply |= ImGui::MenuItem("on Residue");
                    if (ImGui::IsItemHovered()) {
                        const int64_t res_idx = data->mold.mol.atom.res_idx[idx];
                        const md_range_t range = md_residue_atom_range(data->mold.mol.residue, res_idx);
                        md_bitfield_set_range(&mask, range.beg, range.end);
                    }
                }

                if (data->mold.mol.chain.count > 0 && data->mold.mol.atom.chain_idx && data->mold.mol.atom.chain_idx[idx] != -1) {
                    apply |= ImGui::MenuItem("on Chain");
                    if (ImGui::IsItemHovered()) {
                        const auto chain_idx = data->mold.mol.atom.chain_idx[idx];
                        const auto range = md_chain_atom_range(data->mold.mol.chain, chain_idx);
                        md_bitfield_set_range(&mask, range.beg, range.end);
                    }
                }

                if (num_atoms_selected > 0) {
                    apply |= ImGui::MenuItem("on Selection");
                    if (ImGui::IsItemHovered()) {
                        md_bitfield_copy(&mask, &data->selection.selection_mask);
                    }
                }

                if (!md_bitfield_empty(&mask)) {
                    md_bitfield_copy(&data->selection.highlight_mask, &mask);

                    if (apply) {
                        load::traj::set_recenter_target(data->mold.traj, &mask);
                        load::traj::clear_cache(data->mold.traj);
                        interpolate_atomic_properties(data);
                        data->mold.dirty_buffers |= MolBit_ClearVelocity;
                        ImGui::CloseCurrentPopup();
                    }
                }
                ImGui::EndMenu();
            }
        }
        if (ImGui::BeginMenu("Selection")) {
            if (ImGui::MenuItem("Invert")) {
                md_bitfield_not_inplace(&data->selection.selection_mask, 0, data->mold.mol.atom.count);
                ImGui::CloseCurrentPopup();
            }
            if (ImGui::IsItemHovered()) {
                md_bitfield_not(&data->selection.highlight_mask, &data->selection.selection_mask, 0, data->mold.mol.atom.count);
            }
            if (ImGui::MenuItem("Query")) {
                data->selection.query.show_window = true;
                ImGui::CloseCurrentPopup();
            }
            if (num_atoms_selected > 0) {
                if (ImGui::MenuItem("Grow")) {
                    data->selection.grow.show_window = true;
                    ImGui::CloseCurrentPopup();
                }
                if (ImGui::MenuItem("Clear")) {
                    md_bitfield_clear(&data->selection.selection_mask);
                }
            }
            ImGui::EndMenu();
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
        const char* fmt = (data->selection.grow.mode == SelectionGrowth::CovalentBond) ? "%.0f" : "%.2f";
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
            case SelectionGrowth::CovalentBond:
                md_util_mask_grow_by_bonds(&data->selection.grow.mask, &data->mold.mol, (int)data->selection.grow.extent, &data->representation.visibility_mask);
                break;
            case SelectionGrowth::Radial: {
                md_util_mask_grow_by_radius(&data->selection.grow.mask, &data->mold.mol, data->selection.grow.extent, &data->representation.visibility_mask);
                break;
            }
            default:
                ASSERT(false);
            }

            grow_mask_by_selection_granularity(&data->selection.grow.mask, data->selection.granularity, data->mold.mol);
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
            data->selection.query.query_ok = filter_expression(data, str_from_cstr(data->selection.query.buf), &data->selection.query.mask, NULL, data->selection.query.error, sizeof(data->selection.query.error));
            query_frame = data->animation.frame;

            if (data->selection.query.query_ok) {
                grow_mask_by_selection_granularity(&data->selection.query.mask, data->selection.granularity, data->mold.mol);
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
    int num_frames = (int)md_trajectory_num_frames(data->mold.traj);
    if (num_frames == 0) return;

    ASSERT(data->timeline.x_values);
    ASSERT(md_array_size(data->timeline.x_values) == num_frames);

    ImGui::SetNextWindowSize({300,200}, ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Animation", &data->animation.show_window)) {
        const float item_width = MAX(ImGui::GetContentRegionAvail().x - 80.f, 100.f);
        ImGui::PushItemWidth(item_width);

        ImGui::Text("Num Frames: %i", num_frames);

        md_unit_t time_unit = md_trajectory_time_unit(data->mold.traj);
        double t   = frame_to_time(data->animation.frame, *data);
        double min = data->timeline.x_values[0];
        double max = data->timeline.x_values[num_frames - 1];
        char time_label[64];
        if (md_unit_empty(time_unit)) {
            snprintf(time_label, sizeof(time_label), "Time");
        } else {
            char unit_buf[32];
            md_unit_print(unit_buf, sizeof(unit_buf), time_unit);
            snprintf(time_label, sizeof(time_label), "Time (%s)", unit_buf);
        }
        if (ImGui::BeginCombo("Interp.", interpolation_mode_str[(int)data->animation.interpolation])) {
            for (int i = 0; i < (int)InterpolationMode::Count; ++i) {
                if (ImGui::Selectable(interpolation_mode_str[i], (int)data->animation.interpolation == i)) {
                    data->animation.interpolation = (InterpolationMode)i;
                    data->mold.dirty_buffers |= MolBit_ClearVelocity;
                    interpolate_atomic_properties(data);
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
                ImGui::SetTooltip("Tension of the Cubic Spline");
            }
        }
        switch (data->animation.mode) {
            case PlaybackMode::Playing:
                if (ImGui::Button((const char*)ICON_FA_PAUSE)) data->animation.mode = PlaybackMode::Stopped;
                break;
            case PlaybackMode::Stopped:
                if (ImGui::Button((const char*)ICON_FA_PLAY)) data->animation.mode = PlaybackMode::Playing;
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

static void draw_representations_window(ApplicationState* state) {
    if (!state->representation.show_window) return;

    ImGui::SetNextWindowSize({300,200}, ImGuiCond_FirstUseEver);
    ImGui::Begin("Representations", &state->representation.show_window, ImGuiWindowFlags_NoFocusOnAppearing);
    if (ImGui::Button("create new")) {
        create_representation(state);
    }
    ImGui::SameLine();
    if (ImGui::DeleteButton("remove all")) {
        clear_representations(state);
    }
    ImGui::Spacing();
    ImGui::Separator();
    for (int rep_idx = 0; rep_idx < (int)md_array_size(state->representation.reps); rep_idx++) {
        bool update_rep = false;
        Representation& rep = state->representation.reps[rep_idx];
        const float item_width = MAX(ImGui::GetContentRegionAvail().x - 125.f, 100.f);
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
            ImGui::PushItemWidth(item_width);

            // @TODO: Only display the representations which can be used for the current dataset
            if (!rep.type_is_valid) ImGui::PushInvalid();
            if (ImGui::BeginCombo("type", representation_type_str[(int)rep.type])) {
                for (int i = 0; i < (int)RepresentationType::Count; ++i) {
                    switch (i) {
                    case (int)RepresentationType::Orbital:
                        if (md_array_size(state->representation.info.molecular_orbitals) == 0) continue;
                        break;
                    case (int)RepresentationType::DipoleMoment:
                        if (md_array_size(state->representation.info.dipole_moments) == 0) continue;
                        break;
                    default:
                        break;
                    }
                    if (ImGui::Selectable(representation_type_str[(int)i], i == (int)rep.type)) {
                        rep.type = (RepresentationType)i;
                        update_rep = true;
                    }
                }
                ImGui::EndCombo();
            }
            if (!rep.type_is_valid) ImGui::PopInvalid();

            if (rep.type <= RepresentationType::Cartoon) {
                if (ImGui::InputQuery("filter", rep.filt, sizeof(rep.filt), rep.filt_is_valid, rep.filt_error)) {
                    rep.filt_is_dirty = true;
                    update_rep = true;
                }
                if (ImGui::Combo("color", (int*)(&rep.color_mapping), color_mapping_str, IM_ARRAYSIZE(color_mapping_str))) {
                    update_rep = true;
                }
#if 0
                if (rep.color_mapping == ColorMapping::Property) {
                    // @TODO: Update this with something more proper, and use a filter window to allow the user to specify a property
                    /*
                    if (!rep.prop_is_valid) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
                    if (ImGui::InputText("property", rep.prop.cstr(), rep.prop.capacity())) {
                        update_color = true;
                    }
                    if (ImGui::IsItemHovered() && !rep.prop_is_valid) {
                        ImGui::SetTooltip("%s", rep.prop_error.cstr());
                    }
                    if (!rep.prop_is_valid) ImGui::PopStyleColor();
                    */

                    static int prop_idx = 0;
                    const md_script_property_t* props[32] = {0};
                    size_t num_props = 0;
                    for (size_t j = 0; j < md_array_size(data->display_properties); ++j) {
                        if (data->display_properties[j].type == DisplayProperty::Type_Temporal) {
                            props[num_props++] = data->display_properties[j].prop;
                        }
                        if (num_props == ARRAY_SIZE(props)) break;
                    }

                    rep.prop = NULL;
                    if (num_props > 0) {
                        prop_idx = CLAMP(prop_idx, 0, (int)num_props-1);
                        if (ImGui::BeginCombo("Prop", props[prop_idx]->ident.ptr)) {
                            for (size_t j = 0; j < num_props; ++j) {
                                if (ImGui::Selectable(props[j]->ident.ptr, prop_idx == i)) {
                                    prop_idx = (int)j;
                                    rep.map_beg = props[j]->data.min_value;
                                    rep.map_end = props[j]->data.max_value;
                                    update_rep = true;
                                }
                            }
                            ImGui::EndCombo();
                        }
                        rep.prop = props[prop_idx];

                        if (ImPlot::ColormapButton(ImPlot::GetColormapName(rep.color_map), ImVec2(item_width,0), rep.color_map)) {
                            ImGui::OpenPopup("Color Map Selector");
                        }
                        ImGui::DragFloatRange2("Min / Max",&rep.map_beg, &rep.map_end, 0.01f, rep.prop->data.min_value, rep.prop->data.max_value);
                        if (ImGui::BeginPopup("Color Map Selector")) {
                            for (int map = 0; map < ImPlot::GetColormapCount(); ++map) {
                                if (ImPlot::ColormapButton(ImPlot::GetColormapName(map), ImVec2(item_width,0), map)) {
                                    rep.color_map = map;
                                    ImGui::CloseCurrentPopup();
                                }
                            }
                            ImGui::EndPopup();
                        }
                    }
                }
#endif
                if (rep.filt_is_dynamic || rep.color_mapping == ColorMapping::Property) {
                    update_rep |= ImGui::Checkbox("auto-update", &rep.dynamic_evaluation);
                    if (!rep.dynamic_evaluation) {
                        ImGui::SameLine();
                        if (ImGui::Button("update")) {
                            rep.filt_is_dirty = true;
                            update_rep = true;
                        }
                    }
                } else {
                    rep.dynamic_evaluation = false;
                }
                ImGui::PopItemWidth();
                if (rep.color_mapping == ColorMapping::Uniform) {
                    update_rep |= ImGui::ColorEdit4("color", (float*)&rep.uniform_color, ImGuiColorEditFlags_NoInputs);
                }
                ImGui::PushItemWidth(item_width);
                update_rep |= ImGui::SliderFloat("saturation", &rep.saturation, 0.0f, 1.0f);
                switch (rep.type) {
                case RepresentationType::SpaceFill:
                    update_rep |= ImGui::SliderFloat("scale", &rep.scale[0], 0.1f, 4.f);
                    break;
                case RepresentationType::Licorice:
                    update_rep |= ImGui::SliderFloat("radius", &rep.scale[0], 0.1f, 4.0f);
                    break;
                case RepresentationType::BallAndStick:
                    update_rep |= ImGui::SliderFloat("ball scale", &rep.scale[0], 0.1f, 4.f);
                    update_rep |= ImGui::SliderFloat("bond scale", &rep.scale[1], 0.1f, 4.f);
                    break;
                case RepresentationType::Ribbons:
                    update_rep |= ImGui::SliderFloat("width",       &rep.scale[0], 0.1f, 3.f);
                    update_rep |= ImGui::SliderFloat("thickness",   &rep.scale[1], 0.1f, 3.f);
                    break;
                case RepresentationType::Cartoon:
                    update_rep |= ImGui::SliderFloat("coil scale",  &rep.scale[0], 0.1f, 3.f);
                    update_rep |= ImGui::SliderFloat("sheet scale", &rep.scale[1], 0.1f, 3.f);
                    update_rep |= ImGui::SliderFloat("helix scale", &rep.scale[2], 0.1f, 3.f);
                    break;
                default:
                    ASSERT(false);
                }

                ImGui::PopItemWidth();
                ImGui::Spacing();
                ImGui::Separator();
            }

            if (rep.type == RepresentationType::Orbital) {
                ImGuiComboFlags flags = 0;
                if (ImGui::BeginCombo("Orbital Type", orbital_type_str[(int)rep.orbital.type], flags)) {
                    for (int n = 0; n < (int)OrbitalType::Count; n++) {
                        const bool is_selected = ((int)rep.orbital.type == n);
                        if (ImGui::Selectable(orbital_type_str[n], is_selected)) {
                            rep.orbital.type = (OrbitalType)n;
                            update_rep = true;
                        }

                        if (is_selected) {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }

                if (rep.orbital.type == OrbitalType::MolecularOrbitalPsi || rep.orbital.type == OrbitalType::MolecularOrbitalPsiSquared) {
                    char lbl[32];
                    auto write_lbl = [&rep_info = state->representation.info, &lbl](int idx) -> const char* {
                        const char* suffix = "";
                        if (idx == rep_info.mo_homo_idx) {
                            suffix = "(HOMO)";
                        } else if (idx == rep_info.mo_lumo_idx) {
                            suffix = "(LUMO)";
                        }
                        snprintf(lbl, sizeof(lbl), "%i %s", idx + 1, suffix);
                        return lbl;
                    };

                    write_lbl(rep.orbital.orbital_idx);
                    if (ImGui::BeginCombo("Orbital Idx", lbl)) {
                        for (int n = 0; n < (int)md_array_size(state->representation.info.molecular_orbitals); n++) {
                            int idx = state->representation.info.molecular_orbitals[n].idx;
                            const bool is_selected = (rep.orbital.orbital_idx == idx);
                        
                            write_lbl(idx);
                            if (ImGui::Selectable(lbl, is_selected)) {
                                if (rep.orbital.orbital_idx != idx) {
                                    update_rep = true;
                                }
                                rep.orbital.orbital_idx = idx;
                            }

                            if (is_selected) {
                                ImGui::SetItemDefaultFocus();
                            }
                        }
                        ImGui::EndCombo();
                    }
                }

                if (ImGui::Combo("Volume Resolution", (int*)&rep.orbital.resolution, volume_resolution_str, IM_ARRAYSIZE(volume_resolution_str))) {
                    update_rep = true;
                }
#if 0
                // Currently we do not expose DVR, since we do not have a good way of exposing the alpha ramp for the transfer function...
                ImGui::Checkbox("Enable DVR", &rep.orbital.dvr.enabled);
                if (rep.orbital.dvr.enabled) {
                    const ImVec2 button_size = {160, 0};
                    if (ImPlot::ColormapButton(ImPlot::GetColormapName(rep.orbital.dvr.colormap), button_size, rep.orbital.dvr.colormap)) {
                        ImGui::OpenPopup("Colormap Selector");
                    }
                    if (ImGui::BeginPopup("Colormap Selector")) {
                        for (int map = 4; map < ImPlot::GetColormapCount(); ++map) {
                            if (ImPlot::ColormapButton(ImPlot::GetColormapName(map), button_size, map)) {
                                rep.orbital.dvr.colormap = map;
                                update_rep = true;
                                ImGui::CloseCurrentPopup();
                            }
                        }
                        ImGui::EndPopup();
                    }
                }
#endif

                //ImGui::Checkbox("Enable Iso-Surface", &rep.orbital.vol.iso.enabled);

                if (rep.orbital.type == OrbitalType::MolecularOrbitalPsi) {
                    const double iso_min = 1.0e-4;
                    const double iso_max = 5.0;
                    double iso_val = rep.orbital.iso_psi.values[0];
                    if (ImGui::SliderScalar("Iso Value", ImGuiDataType_Double, &iso_val, &iso_min, &iso_max, "%.6f", ImGuiSliderFlags_Logarithmic)) {
                        rep.orbital.iso_psi.values[0] =  (float)iso_val;
                        rep.orbital.iso_psi.values[1] = -(float)iso_val;
                        rep.orbital.iso_den.values[0] =  (float)(iso_val * iso_val);
                    }
                } else if (rep.orbital.type == OrbitalType::MolecularOrbitalPsiSquared ||
                           rep.orbital.type == OrbitalType::ElectronDensity ||
                           rep.orbital.type == OrbitalType::AverageLocalIonizationEnergy) {
                    const double iso_min = 1.0e-8;
                    const double iso_max = 5.0;
                    double iso_val = rep.orbital.iso_den.values[0];

                    if (ImGui::SliderScalar("Iso Value", ImGuiDataType_Double, &iso_val, &iso_min, &iso_max, "%.6f", ImGuiSliderFlags_Logarithmic)) {
                        rep.orbital.iso_psi.values[0] =  (float)sqrt(iso_val);
                        rep.orbital.iso_psi.values[1] = -(float)sqrt(iso_val);
                        rep.orbital.iso_den.values[0] =  (float)iso_val;
                        rep.orbital.iso_den.values[1] =  (float)iso_val;
                    }
                }

                if (rep.orbital.type == OrbitalType::MolecularOrbitalPsi) {
                    ImGui::ColorEdit4("Color Positive", rep.orbital.iso_psi.colors[0].elem);
                    ImGui::ColorEdit4("Color Negative", rep.orbital.iso_psi.colors[1].elem);
                } else if (rep.orbital.type == OrbitalType::AverageLocalIonizationEnergy) {
                    update_rep |= ImPlot::ColormapSelection("Colormap", &rep.orbital.alie.colormap);
                    update_rep |= ImGui::SliderFloat("Alpha", &rep.orbital.alie.tf_alpha, 0.0f, 1.0f);
                    //ImGui::RangeSliderFloat("Value Range", &rep.orbital.alie.min_tf_val, &rep.orbital.alie.max_tf_val, -0.01f, 0.01f);
                    ImGui::SliderFloat("Value Min", &rep.orbital.alie.min_tf_val, -0.1f, 0.1f, "%.8f");
                    ImGui::SliderFloat("Value Max", &rep.orbital.alie.max_tf_val, -0.1f, 0.1f, "%.8f");
                    rep.orbital.alie.max_tf_val = MAX(rep.orbital.alie.min_tf_val, rep.orbital.alie.max_tf_val);
                } else {
                    ImGui::ColorEdit4("Color Density",  rep.orbital.iso_den.colors[0].elem);
                    rep.orbital.iso_den.colors[1] = rep.orbital.iso_den.colors[0];
                }
            }
            ImGui::TreePop();
        }

        ImGui::PopID();

        if (update_rep) {
            update_representation(state, &rep);
        }
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
                                                       ImGui::GetCurrentContext()->FontBaseSize + ImGui::GetStyle().FramePadding.y * 2.f + MARGIN));
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
                y_vals[i] = (args.value_filter.min < val && val < args.value_filter.max) ? args.values.max_y : -INFINITY;
            }

            const ImVec4 filter_frame_color = ImVec4(1,1,1,1);
            ImPlot::SetNextFillStyle(filter_frame_color, 0.15f);
            ImPlot::PlotShaded("##value_filter", args.values.x, y_vals, args.values.count, -INFINITY);
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

            md_unit_t time_unit = md_trajectory_time_unit(data.mold.traj);
            if (!md_unit_empty(time_unit)) {
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

struct DisplayPropertyDragDropPayload {
    int prop_idx = 0;
    int src_plot_idx = -1;
};

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

static float distance_to_linesegment(vec2_t line_beg, vec2_t line_end, vec2_t point) {

    vec2_t v = vec2_sub(line_end, line_beg);
    vec2_t u = vec2_sub(point, line_beg);
    float dot = vec2_dot(u,v);
    float len2 = vec2_dot(v,v);

    if (len2 < 1.0e-5f) {
        return vec2_dist(point, line_beg);
    }

    float t = dot / len2;
    if (t < 0.0f) {
        return vec2_dist(point, line_beg);
    } else if (t > 1.0f) {
        return vec2_dist(point, line_end);
    } else {
        return vec2_dist(point, vec2_add(line_beg, vec2_mul_f(v, t)));
    }
}

static void visualize_payload(ApplicationState* data, const md_script_vis_payload_o* payload, int subidx, md_script_vis_flags_t flags) {
    md_script_vis_ctx_t ctx = {
        .ir   = data->script.eval_ir,
        .mol  = &data->mold.mol,
        .traj = data->mold.traj,
    };
    data->script.vis = {0};
    md_script_vis_init(&data->script.vis, frame_alloc);

    if (md_script_vis_eval_payload(&data->script.vis, payload, subidx, &ctx, flags)) {
        if (!md_bitfield_empty(&data->script.vis.atom_mask)) {
            md_bitfield_copy(&data->selection.highlight_mask, &data->script.vis.atom_mask);
        }
    }
}

// #timeline
static void draw_timeline_window(ApplicationState* data) {
    ASSERT(data);
    ImGui::SetNextWindowSize(ImVec2(600, 300), ImGuiCond_FirstUseEver);

    if (ImGui::Begin("Timelines", &data->timeline.show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_MenuBar)) {
        static int num_subplots = 1;

        double pre_filter_min = data->timeline.filter.beg_frame;
        double pre_filter_max = data->timeline.filter.end_frame;

        const float* x_values   = data->timeline.x_values;
        const int num_x_values  = (int)md_array_size(data->timeline.x_values);
        const float min_x_value = num_x_values > 0 ? x_values[0] : 0.0f;
        const float max_x_value = num_x_values > 0 ? x_values[num_x_values - 1] : 1.0f;

        ImPlot::PushStyleVar(ImPlotStyleVar_PlotPadding, ImVec2(ImPlot::GetStyle().PlotPadding.x, 2));
        defer { ImPlot::PopStyleVar(); };

        // Filter out temporal display properties
        int num_temp_props = 0;
        for (int i = 0; i < md_array_size(data->display_properties); ++i) {
            if (data->display_properties[i].type == DisplayProperty::Type_Temporal) {
                num_temp_props += 1;
            }
        }

        const int num_props = (int)md_array_size(data->display_properties);

        if (ImGui::BeginMenuBar()) {            
            if (ImGui::BeginMenu("Properties")) {
                if (num_temp_props) {
                    for (int i = 0; i < num_props; ++i) {
                        DisplayProperty& dp = data->display_properties[i];
                        if (dp.type != DisplayProperty::Type_Temporal) continue;

                        ImPlot::ItemIcon(dp.color);
                        ImGui::SameLine();
                        ImGui::Selectable(dp.label);

                        if (ImGui::IsItemHovered()) {
                            if ((dp.dim > MAX_POPULATION_SIZE)) {
                                ImGui::SetTooltip("The property has a large population, only the first %i items will be shown", MAX_POPULATION_SIZE);
                            }
                            visualize_payload(data, dp.vis_payload, -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                            set_hovered_property(data, str_from_cstr(dp.label));
                        }

                        if (ImGui::BeginDragDropSource()) {
                            DisplayPropertyDragDropPayload payload = {i};
                            ImGui::SetDragDropPayload("TEMPORAL_DND", &payload, sizeof(payload));
                            ImPlot::ItemIcon(dp.color); ImGui::SameLine();
                            ImGui::TextUnformatted(dp.label);
                            ImGui::EndDragDropSource();
                        }
                    }
                } else {
                    ImGui::Text("No temporal properties available, define and evaluate properties in the script editor");
                }
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
                ImGui::SliderInt("Num Subplots", &num_subplots, 1, MAX_TEMPORAL_SUBPLOTS);
                if (ImGui::Button("Add Subplot")) {
                    num_subplots = CLAMP(num_subplots + 1, 1, MAX_TEMPORAL_SUBPLOTS);
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
            md_unit_t x_unit = md_trajectory_time_unit(data->mold.traj);
            if (!md_unit_empty(x_unit)) {
                md_unit_print(x_unit_str, sizeof(x_unit_str), x_unit);
                snprintf(x_label, sizeof(x_label), "Time (%s)", x_unit_str);
            }

            for (int i = 0; i < num_subplots; ++i) {
                if (ImPlot::BeginPlot("", ImVec2(), plot_flags)) {
                    ImPlot::SetupAxisLinks(ImAxis_X1, &data->timeline.view_range.beg_x, &data->timeline.view_range.end_x);
                    ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, min_x_value, max_x_value);

                    ImPlotAxisFlags axis_flags_x = axis_flags | ImPlotAxisFlags_NoLabel;
                    if (i < num_subplots - 1) {
                        // Only show label and ticklabels for the last plot, since they are all synced on x-axis
                        axis_flags_x |= ImPlotAxisFlags_NoTickLabels;
                    }

                    // Check and see if every property within the current subplot share the same unit, if so, use it as y_label
                    md_unit_t y_unit = md_unit_none();
                    for (int j = 0; j < num_props; ++j) {
                        DisplayProperty& prop = data->display_properties[j];
                        if (prop.temporal_subplot_mask & (1 << i)) {
                            if (md_unit_equal(y_unit, md_unit_none())) {
                                y_unit = prop.unit[1];
                            } else if (!md_unit_equal(y_unit, prop.unit[1])) {
                                // unit conflict, drop it
                                y_unit = md_unit_none();
                                break;
                            }
                        }
                    }

                    char y_label[64] = "";
                    char y_unit_str[32] = "";
                    if (!md_unit_equal(y_unit, md_unit_none())) {
                        md_unit_print(y_unit_str, sizeof(y_unit_str), y_unit);
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

                    // Find the and set the index of hovered lines within the plot
                    int  hovered_prop_idx = -1;
                    int  hovered_pop_idx  = -1; // Population index (in the case that the property has a population of values (dim > 1))
                    char hovered_label[64] = "";
                    
                    bool print_timeline_tooltip = false;
                   
                    if (ImPlot::IsPlotHovered()) {
                        md_bitfield_clear(&data->selection.highlight_mask);
                        set_hovered_property(data,  STR_LIT(""));

                        print_timeline_tooltip = true;
                        const ImPlotPoint mouse_pos = ImPlot::GetPlotMousePos();
                        const ImVec2 mouse_coord = ImPlot::PlotToPixels(mouse_pos);
                        const double frame = time_to_frame(mouse_pos.x, data->timeline.x_values);
                        const int fn  = CLAMP((int)(frame + 0.5), 0, num_x_values - 1); // Nearest index
                        const int f[4] = {
                            CLAMP((int)frame - 1,   0, num_x_values - 1),
                            CLAMP((int)frame,       0, num_x_values - 1),
                            CLAMP((int)frame + 1,   0, num_x_values - 1),
                            CLAMP((int)frame + 2,   0, num_x_values - 1),
                        };
                        const float max_rad = 20; // 20 pixels
                        const float area_dist = max_rad * 0.2;

                        float min_dist = max_rad;
                    
                        for (int j = 0; j < num_props; ++j) {
                            DisplayProperty& prop = data->display_properties[j];
                            
                            ImPlotItem* item = ImPlot::GetItem(prop.label);
                            if (!item || !item->Show) {
                                continue;
                            }

                            DisplayProperty::Payload payload = {
                                .display_prop = &prop,
                            };
                            
                            if (prop.temporal_subplot_mask & (1 << i)) {
                                const int dim = CLAMP(1, prop.dim, MAX_POPULATION_SIZE);
                                for (int k = 0; k < dim; ++k) {
                                    if (dim > 1 && !(prop.population_mask.test(k))) {
                                        continue;
                                    }
                                    payload.dim_idx = k;
                                    double d = DBL_MAX;

                                    switch (prop.plot_type) {
                                    case DisplayProperty::PlotType_Line:
                                    {
                                        // Compute distance to line segments, prev, cur and next
                                        // It is not sufficient to only check the distance to the current line segment
                                        ImVec2 p[4] = {
                                            ImPlot::PlotToPixels(prop.getter[0](f[0], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[0](f[1], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[0](f[2], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[0](f[3], &payload)),
                                        };
                                        d = distance_to_linesegment(p[0], p[1], mouse_coord);
                                        d = MIN(distance_to_linesegment(p[1], p[2], mouse_coord), d);
                                        d = MIN(distance_to_linesegment(p[2], p[3], mouse_coord), d);
                                                
                                        break;
                                    }
                                    case DisplayProperty::PlotType_Area:
                                    {
                                        const ImVec2 p_min[4] = {
                                            ImPlot::PlotToPixels(prop.getter[0](f[0], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[0](f[1], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[0](f[2], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[0](f[3], &payload)),
                                        };
                                        const ImVec2 p_max[4] = {
                                            ImPlot::PlotToPixels(prop.getter[1](f[0], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[1](f[1], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[1](f[2], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[1](f[3], &payload)),
                                        };
                                        // Check if within area
                                        for (int l = 0; l < 2; ++l) {
                                            // Each segment forms a trapetzoid with left and right half parallel to the y axis
                                            // We want to clamp the mouse coordinate to the trapetzoid and compute the distance to the clamped point
                        
                                            const float x_min = MIN(p_min[l].x, p_min[l+1].x);
                                            const float x_max = MAX(p_min[l].x, p_min[l+1].x);

                                            // Bilinarly interpolate the y min/max
                                            const float t = CLAMP((mouse_coord.x - x_min) / (x_max - x_min), 0.0f, 1.0f);
                                            const float y[2] = {
                                                lerp(p_min[l].y, p_min[l+1].y, t),
                                                lerp(p_max[l].y, p_max[l+1].y, t)
                                            };
                                            const float y_min = MIN(y[0], y[1]);
                                            const float y_max = MAX(y[0], y[1]);

                                            ImVec2 p = {
                                                CLAMP(mouse_coord.x, x_min, x_max),
                                                CLAMP(mouse_coord.y, y_min, y_max)
                                            };

                                            d = MIN(d, sqrt(ImLengthSqr(mouse_coord - p)));
                                        }
                                        d += area_dist;
                                        break;
                                    }
                                    case DisplayProperty::PlotType_Scatter:
                                    {
                                        ImVec2 p = ImPlot::PlotToPixels(prop.getter[0](fn, &payload));
                                        d = sqrt(ImLengthSqr(mouse_coord - p));
                                        break;
                                    }
                                    default:
                                        // Should not end up here
                                        ASSERT(false);
                                        break;
                                    }

                                    if (d < min_dist) {
                                        min_dist = d;
                                        char value_buf[64] = "";
                                        if (prop.print_value) {
                                            prop.print_value(value_buf, sizeof(value_buf), fn, &payload);
                                        } else {
                                            ImPlotPoint p = prop.getter[0](fn, &payload);
                                            snprintf(value_buf, sizeof(value_buf), "%.2f", p.y);
                                        }

                                        hovered_prop_idx = j;
                                        hovered_pop_idx = k;
                                        if (prop.dim > 1) {
                                            snprintf(hovered_label, sizeof(hovered_label), "%s[%i]: %s", prop.label, k + 1, value_buf);
                                        } else {
                                            snprintf(hovered_label, sizeof(hovered_label), "%s: %s", prop.label, value_buf);
                                        }
                                    }
                                }
                            }
                        }

                        if (hovered_prop_idx != -1) {
                            set_hovered_property(data, str_from_cstr(data->display_properties[hovered_prop_idx].label), hovered_pop_idx);
                        }

                        if (int len = (int)strnlen(hovered_label, sizeof(hovered_label))) {
                            // Concat the hovered_label with the y-unit
                            snprintf(hovered_label + len, (int)sizeof(hovered_label) - len, " %s", y_label);
                        }
                    } else {
                        if (!str_empty(data->hovered_display_property_label)) {
                            for (int j = 0; j < num_props; ++j) {
                                DisplayProperty& dp = data->display_properties[j];
                                if (dp.type != DisplayProperty::Type_Temporal) continue;
                                if (str_eq_cstr(data->hovered_display_property_label, dp.label)) {
                                    hovered_prop_idx = j;
                                    hovered_pop_idx = data->hovered_display_property_pop_idx;
                                    break;
                                }
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
                    
                    for (int j = 0; j < num_props; ++j) {
                        DisplayProperty& dp = data->display_properties[j];
                        if ((dp.type != DisplayProperty::Type_Temporal)) continue;
                        if (!(dp.temporal_subplot_mask & (1 << i))) continue;

                        if (ImPlot::IsLegendEntryHovered(dp.label)) {
                            visualize_payload(data, dp.vis_payload, -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                            set_hovered_property(data, str_from_cstr(dp.label));
                            hovered_prop_idx = j;
                            hovered_pop_idx = -1;
                        }

                        // legend context menu
                        if (ImPlot::BeginLegendPopup(dp.label)) {
                            if (ImGui::DeleteButton("Remove")) {
                                dp.temporal_subplot_mask &= ~(1 << i);
                                ImGui::CloseCurrentPopup();
                            }

                            const char* plot_type_names[] = {"Line", "Area", "Bars", "Scatter"};
                            const bool  valid_plot_types[] = {true, false, false, true};
                            STATIC_ASSERT(ARRAY_SIZE(plot_type_names) == DisplayProperty::PlotType_Count);
                            STATIC_ASSERT(ARRAY_SIZE(valid_plot_types) == DisplayProperty::PlotType_Count);

                            // The user only has the option to choose between line and scatter if it is
                            // Line or scatter, which is initially determined by its type
                            if (valid_plot_types[dp.plot_type]) {
                                if (ImGui::BeginCombo("Plot Type", plot_type_names[dp.plot_type])) {
                                    for (int k = 0; k < DisplayProperty::PlotType_Count; ++k) {
                                        if (!valid_plot_types[k]) continue;
                                        if (ImGui::Selectable(plot_type_names[k], dp.plot_type == k)) {
                                            dp.plot_type = (DisplayProperty::PlotType)k;
                                        }
                                    }
                                    ImGui::EndCombo();
                                }
                            }

                            if (dp.plot_type == DisplayProperty::PlotType_Scatter) {
                                if (ImGui::BeginCombo("Marker", ImPlot::GetMarkerName(dp.marker_type))) {
                                    for (int k = 0; k < ImPlotMarker_COUNT; ++k) {
                                        if (ImGui::Selectable(ImPlot::GetMarkerName(k), dp.marker_type == k)) {
                                            dp.marker_type = (ImPlotMarker)k;
                                        }
                                    }
                                    ImGui::EndCombo();
                                }
                                ImGui::SliderFloat("Marker Size", &dp.marker_size, 0.1f, 10.0f, "%.2f");
                            }

                            if (dp.dim > 1) {
                                const char* color_type_labels[] = {"Solid", "Colormap"};
                                STATIC_ASSERT(ARRAY_SIZE(color_type_labels) == DisplayProperty::ColorType_Count);

                                if (ImGui::BeginCombo("Color Type", color_type_labels[dp.color_type])) {
                                    for (int k = 0; k < DisplayProperty::ColorType_Count; ++k) {
                                        if (ImGui::Selectable(color_type_labels[k], k == dp.color_type)) {
                                            dp.color_type = (DisplayProperty::ColorType)k;
                                        }
                                    }
                                    ImGui::EndCombo();
                                }
                            }
                            switch (dp.color_type) {
                            case DisplayProperty::ColorType_Solid:
                                ImGui::ColorEdit4("Color", &dp.color.x);
                                break;
                            case DisplayProperty::ColorType_Colormap:
                                ImPlot::ColormapSelection("##Colormap", &dp.colormap);
                                ImGui::SliderFloat("Alpha", &dp.colormap_alpha, 0.0f, 1.0f);
                                break;
                            default:
                                ASSERT(false);
                                break;
                            } 
                            if (dp.dim > 1) {
                                ImGui::Separator();
                                if (ImGui::Button("Set All")) {
                                    dp.population_mask.set();
                                }
                                ImGui::SameLine();
                                if (ImGui::Button("Clear All")) {
                                    dp.population_mask.reset();
                                }

                                const float sz = ImGui::GetFontSize() * 1.5f;
                                ImGui::PushStyleVar(ImGuiStyleVar_SelectableTextAlign, ImVec2(0.5f, 0.5f));
                                for (int k = 0; k < MIN(dp.dim, MAX_POPULATION_SIZE); ++k) {
                                    char lbl[32];
                                    snprintf(lbl, sizeof(lbl), "%d", k+1);
                                    if (ImGui::Selectable(lbl, dp.population_mask.test(k), ImGuiSelectableFlags_DontClosePopups, ImVec2(sz, sz))) {
                                        // Toggle bit for this population index
                                        dp.population_mask.flip(k);
                                    }
                                    if (ImGui::IsItemHovered()) {
                                        visualize_payload(data, dp.vis_payload, k, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                                        set_hovered_property(data, str_from_cstr(dp.label), k);
                                        hovered_prop_idx = j;
                                        hovered_pop_idx = k;
                                    }
                                    if (!k || ((k+1) % 10)) {
                                        ImGui::SameLine();
                                    }
                                }
                                ImGui::PopStyleVar();
                            }
                            ImPlot::EndLegendPopup();
                        }

                        auto plot = [j, &dp, hovered_prop_idx, hovered_pop_idx](int k) {
                            const float  hov_fill_alpha  = 1.25f;
                            const float  hov_line_weight = 2.0f;
                            const float  hov_col_scl = 1.5f;
                            const int    population_size = CLAMP(dp.dim, 1, MAX_POPULATION_SIZE);

                            ImVec4 color = {};
                            ImVec4 marker_line_color = {};
                            float  marker_line_weight = 0;
                            float  fill_alpha = 1.0f;
                            float  weight = 1.0f;

                            switch(dp.color_type) {
                            case DisplayProperty::ColorType_Solid:
                                color = dp.color;
                                break;
                            case DisplayProperty::ColorType_Colormap:
                                if (ImPlot::ColormapQualitative(dp.colormap)) {
                                    color = ImPlot::GetColormapColor(k, dp.colormap);
                                } else {
                                    color = ImPlot::SampleColormap( (float)k / (float)(population_size-1), dp.colormap);
                                }
                                color.w *= dp.colormap_alpha;
                                break;
                            default:
                                ASSERT(false);
                                break;
                            }

                            if (hovered_prop_idx == j) {
                                if (hovered_pop_idx == -1 || hovered_pop_idx == k) {
                                    color = ImVec4(ImSaturate(color.x * hov_col_scl), ImSaturate(color.y * hov_col_scl), ImSaturate(color.z * hov_col_scl), color.w);
                                    fill_alpha = hov_fill_alpha;
                                    marker_line_color = {1,1,1,1};
                                    marker_line_weight = 1.0f;
                                }
                                if (hovered_pop_idx == k) {
                                    weight = hov_line_weight;
                                }
                            }

                            DisplayProperty::Payload payload {
                                .display_prop = &dp,
                                .dim_idx = k,
                            };

                            switch (dp.plot_type) {
                            case DisplayProperty::PlotType_Line:
                                ImPlot::SetNextLineStyle(color, weight);
                                ImPlot::PlotLineG(dp.label, dp.getter[0], &payload, dp.num_samples);
                                break;
                            case DisplayProperty::PlotType_Area:
                                ImPlot::SetNextFillStyle(color, fill_alpha);
                                ImPlot::PlotShadedG(dp.label, dp.getter[0], &payload, dp.getter[1], &payload, dp.num_samples);
                                break;
                            case DisplayProperty::PlotType_Scatter:
                                ImPlot::SetNextMarkerStyle(dp.marker_type, dp.marker_size, color, marker_line_weight, marker_line_color);
                                ImPlot::PlotScatterG(dp.label, dp.getter[0], &payload, dp.num_samples);
                                break;
                            default:
                                // Should not end up here
                                ASSERT(false);
                                break;
                            }
                        };

                        // Draw regular lines
                        const int population_size = CLAMP(dp.dim, 1, MAX_POPULATION_SIZE);
                        for (int k = 0; k < population_size; ++k) {
                            if (population_size > 1 && !dp.population_mask.test(k)) {
                                continue;
                            }
                            if (hovered_prop_idx == j && hovered_pop_idx == k) {
                                continue;
                            }

                            plot(k);
                        }

                        // Draw hovered line
                        if (hovered_prop_idx == j && hovered_pop_idx != -1) {
                            plot(hovered_pop_idx);  
                        }

                        if (ImPlot::BeginDragDropSourceItem(dp.label)) {
                            DisplayPropertyDragDropPayload dnd_payload = {j, i};
                            ImGui::SetDragDropPayload("TEMPORAL_DND", &dnd_payload, sizeof(dnd_payload));
                            ImPlot::ItemIcon(dp.color); ImGui::SameLine();
                            ImGui::TextUnformatted(dp.label);
                            ImPlot::EndDragDropSource();
                        }
                    }

                    if (ImPlot::DragLineX(0, &time, ImVec4(1,1,0,1), 1.0f, ImPlotDragToolFlags_NoFit)) {
                        time = CLAMP(time, min_x_value, max_x_value);
                    }

                    if (ImPlot::IsPlotHovered()) {
                        if (hovered_prop_idx != -1) {
                            const int pop_idx = data->display_properties[hovered_prop_idx].dim > 1 ? hovered_pop_idx : -1;
                            visualize_payload(data, data->display_properties[hovered_prop_idx].vis_payload, pop_idx, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                            set_hovered_property(data, str_from_cstr(data->display_properties[hovered_prop_idx].label), hovered_pop_idx);
                        }
                    }

                    if (ImPlot::BeginDragDropTargetPlot()) {
                        if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload("TEMPORAL_DND")) {
                            ASSERT(payload->DataSize == sizeof(DisplayPropertyDragDropPayload));
                            DisplayPropertyDragDropPayload* dnd = (DisplayPropertyDragDropPayload*)(payload->Data);
                            data->display_properties[dnd->prop_idx].temporal_subplot_mask |= (1 << i);
                            if (dnd->src_plot_idx != -1 && dnd->src_plot_idx != i) {
                                // Clear bit from mask representing src plot index (only if it originated from another plot)
                                data->display_properties[dnd->prop_idx].temporal_subplot_mask &= ~(1 << dnd->src_plot_idx);
                            }
                        }
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
                        if (md_unit_empty(x_unit)) {
                            int32_t frame_idx = CLAMP((int)(time_to_frame(t, x_values) + 0.5), 0, num_x_values-1);
                            ImGui::SetTooltip("Frame: %i\n%s", frame_idx, hovered_label);
                        } else {
                            ImGui::SetTooltip("Time: %.2f (%s)\n%s", t, x_unit_str, hovered_label);
                        }
                    }
                    
                    ImPlot::EndPlot();
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

        // Try to handle the case when the user is dragging a payload and not dropping it within a valid target zone.
        // In such case if the property had a source plot index, remove the property from that plot
        if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
            const ImGuiPayload* payload = ImGui::GetDragDropPayload();
            if (payload && payload->IsDataType("TEMPORAL_DND") && !ImGui::IsDragDropPayloadBeingAccepted()) {
                DisplayPropertyDragDropPayload* dnd = (DisplayPropertyDragDropPayload*)(payload->Data);
                if (dnd && dnd->src_plot_idx != -1) {
                    data->display_properties[dnd->prop_idx].temporal_subplot_mask &= ~(1 << dnd->src_plot_idx);
                }
            }
        }
    }
    ImGui::End();
}

// #distribution_window
static void draw_distribution_window(ApplicationState* data) {
    ImGui::SetNextWindowSize(ImVec2(200, 300), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Distributions", &data->distributions.show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_MenuBar)) {
        static int num_subplots = 1;

        if (ImGui::BeginMenuBar())
        {
            DisplayProperty* props = data->display_properties;
            const int num_props = (int)md_array_size(props);

            if (ImGui::BeginMenu("Properties")) {
                if (num_props) {
                    for (int i = 0; i < num_props; ++i) {
                        DisplayProperty& dp = props[i];
                        if (dp.type == DisplayProperty::Type_Distribution) {
                            if (!data->timeline.filter.enabled && dp.partial_evaluation) {
                                // Hide the property as an option if the timeline filter is not enabled and the property is a partial evaluation
                                continue;
                            }
                            ImPlot::ItemIcon(dp.color);
                            ImGui::SameLine();
                            ImGui::Selectable(dp.label);
                            if (ImGui::IsItemHovered()) {
                                visualize_payload(data, dp.vis_payload, -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                                set_hovered_property(data, str_from_cstr(dp.label));
                            }
                            if (ImGui::BeginDragDropSource()) {
                                DisplayPropertyDragDropPayload payload = {i};
                                ImGui::SetDragDropPayload("DISTRIBUTION_DND", &payload, sizeof(payload));
                                ImPlot::ItemIcon(dp.color); ImGui::SameLine();
                                ImGui::TextUnformatted(dp.label);
                                ImGui::EndDragDropSource();
                            }
                        }
                    }
                } else {
                    ImGui::Text("No distribution properties available, try evaluating the script");
                }
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Subplots")) {
                ImGui::SliderInt("Num Subplots", &num_subplots, 1, MAX_DISTRIBUTION_SUBPLOTS);
                if (ImGui::Button("Add Subplot")) {
                    num_subplots = CLAMP(num_subplots + 1, 1, MAX_DISTRIBUTION_SUBPLOTS);
                }
                ImGui::EndMenu();
            }

            ImGui::EndMenuBar();
        }

        ImPlotAxisFlags axis_flags   = ImPlotAxisFlags_NoSideSwitch | ImPlotAxisFlags_NoHighlight;
        ImPlotAxisFlags axis_flags_x = axis_flags | 0;
        ImPlotAxisFlags axis_flags_y = axis_flags | ImPlotAxisFlags_AutoFit;

        ImPlotFlags     plot_flags   = ImPlotFlags_NoBoxSelect | ImPlotFlags_NoFrame;

        const int num_props = (int)md_array_size(data->display_properties);

        if (ImPlot::BeginSubplots("##distribution_plots", num_subplots, 1, ImVec2(-1,-1))) {
            for (int i = 0; i < num_subplots; ++i) {
                if (ImPlot::BeginPlot("", ImVec2(-1,0), plot_flags)) {

                    md_unit_t x_unit = {0};
                    md_unit_t y_unit = {0};
                    for (int j = 0; j < num_props; ++j) {
                        DisplayProperty& prop = data->display_properties[j];
                        if (prop.type != DisplayProperty::Type_Distribution) continue;
                        if (!(prop.distribution_subplot_mask & (1 << i))) continue;

                        if (md_unit_empty(x_unit)) {
                            x_unit = prop.unit[0];
                        } else if (!md_unit_equal(x_unit, prop.unit[0])) {
                            // Set to unitless
                            x_unit = md_unit_none();
                        }
                        if (md_unit_empty(y_unit)) {
                            y_unit = prop.unit[1];
                        } else if (!md_unit_equal(y_unit, prop.unit[1])) {
                            // Set to unitless
                            y_unit = md_unit_none();
                        }
                    }

                    char x_label[64] = "";
                    if (!md_unit_unitless(x_unit)) {
                        md_unit_print(x_label, sizeof(x_label), x_unit);
                    }
                    char y_label[64] = "";
                    if (!md_unit_unitless(y_unit)) {
                        md_unit_print(y_label, sizeof(y_label), y_unit);
                    }

                    ImPlot::SetupAxes(x_label, y_label, axis_flags_x, axis_flags_y);
                    ImPlot::SetupFinish();

                    int  hovered_prop_idx  = -1;
                    int  hovered_pop_idx   = -1;
                    char hovered_label[64] = "";

                    if (ImPlot::IsPlotHovered()) {
                        set_hovered_property(data, STR_LIT(""));
                        
                        const ImPlotPoint mouse_pos = ImPlot::GetPlotMousePos();
                        const ImVec2 mouse_coord = ImPlot::PlotToPixels(mouse_pos);
                        
                        const double max_rad = 20; // 20 pixels
                        const double area_dist = max_rad * 0.2;

                        double min_dist = max_rad;

                        for (int j = 0; j < num_props; ++j) {
                            DisplayProperty& prop = data->display_properties[j];

                            if (prop.hist.x_max <= prop.hist.x_min) continue;  // Collapsed x_axis (probably due to no data currently)
                            
                            // Do the reverse mapping that occurs within getters to go from x-coordinate to index
                            const double scl = (prop.hist.x_max - prop.hist.x_min) / (double)prop.hist.num_bins;
                            const double off = prop.hist.x_min + 0.5 * scl;
                            const double x   = ((mouse_pos.x - off) / scl);
                            const int xn = CLAMP(x + 0.5, 0, prop.hist.num_bins - 1);
                            const int xi[4] {
                                CLAMP((int)x - 1,   0, prop.hist.num_bins - 1),
                                CLAMP((int)x,       0, prop.hist.num_bins - 1),
                                CLAMP((int)x + 1,   0, prop.hist.num_bins - 1),
                                CLAMP((int)x + 2,   0, prop.hist.num_bins - 1),
                            };

                            ImPlotItem* item = ImPlot::GetItem(prop.label);
                            if (!item || !item->Show) {
                                continue;
                            }

                            DisplayProperty::Payload payload = {
                                .display_prop = &prop,
                            };

                            if (prop.distribution_subplot_mask & (1 << i)) {
                                for (int k = 0; k < MIN(prop.hist.dim, MAX_POPULATION_SIZE); ++k) {
                                    if (prop.hist.dim > 1 && !prop.population_mask.test(k)) {
                                        continue;
                                    }
                                    payload.dim_idx = k;
                                    double d = DBL_MAX;
                                    const double layer = j + k / (double)(MAX_POPULATION_SIZE - 1);
                                    const double layer_dist = ((num_props - layer) / num_props) * (max_rad * 0.1);

                                    switch (prop.plot_type) {
                                    case DisplayProperty::PlotType_Line:
                                    {
                                        // Compute distance to line segments, prev, cur and next
                                        // It is not sufficient to only check the distance to the current line segment
                                        ImVec2 p[4] = {
                                            ImPlot::PlotToPixels(prop.getter[1](xi[0], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[1](xi[1], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[1](xi[2], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[1](xi[3], &payload)),
                                        };
                                        d = distance_to_linesegment(p[0], p[1], mouse_coord);
                                        d = MIN(distance_to_linesegment(p[1], p[2], mouse_coord), d);
                                        d = MIN(distance_to_linesegment(p[2], p[3], mouse_coord), d);
                                        d += layer_dist;
                                        break;
                                    }
                                    case DisplayProperty::PlotType_Area:
                                    {
                                        const ImVec2 p_min[4] = {
                                            ImPlot::PlotToPixels(prop.getter[0](xi[0], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[0](xi[1], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[0](xi[2], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[0](xi[3], &payload)),
                                        };
                                        const ImVec2 p_max[4] = {
                                            ImPlot::PlotToPixels(prop.getter[1](xi[0], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[1](xi[1], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[1](xi[2], &payload)),
                                            ImPlot::PlotToPixels(prop.getter[1](xi[3], &payload)),
                                        };
                                        // Check if within area
                                        for (int l = 0; l < 2; ++l) {
                                            // Each segment forms a trapetzoid with left and right half parallel to the y axis
                                            // We want to clamp the mouse coordinate to the trapetzoid and compute the distance to the clamped point

                                            const float x_min = MIN(p_min[l].x, p_min[l+1].x);
                                            const float x_max = MAX(p_min[l].x, p_min[l+1].x);

                                            // Bilinarly interpolate the y min/max
                                            const float t = CLAMP((mouse_coord.x - x_min) / (x_max - x_min), 0.0f, 1.0f);
                                            const float y[2] = {
                                                lerp(p_min[l].y, p_min[l+1].y, t),
                                                lerp(p_max[l].y, p_max[l+1].y, t)
                                            };
                                            const float y_min = MIN(y[0], y[1]);
                                            const float y_max = MAX(y[0], y[1]);

                                            ImVec2 p = {
                                                CLAMP(mouse_coord.x, x_min, x_max),
                                                CLAMP(mouse_coord.y, y_min, y_max)
                                            };

                                            d = MIN(d, sqrt(ImLengthSqr(mouse_coord - p)));
                                        }
                                        d += area_dist + layer_dist;
                                        break;
                                    }
                                    case DisplayProperty::PlotType_Bars:
                                    {
                                        // @NOTE: There is a bug here in the computation of the width of the bars
                                        // And does not work properly atm at different zoom levels.
                                        const float bar_half_width = scl * prop.bar_width_scale * 0.5f;
                                        for (int l = 0; l < 3; ++l) {
                                            const ImVec2 p[2] = {
                                                ImPlot::PlotToPixels(ImPlotPoint(bar_half_width, 0)),
                                                ImPlot::PlotToPixels(prop.getter[1](xi[l], &payload))
                                            };
                                            const ImVec2 p_min = {p[1].x - p[0].x, ImMin(p[0].y, p[1].y)};
                                            const ImVec2 p_max = {p[1].x + p[0].x, ImMax(p[0].y, p[1].y)};

                                            const ImVec2 ll = {p_min.x, p_min.y};
                                            const ImVec2 ur = {p_max.x, p_max.y};
                                            const ImVec2 pos = ImClamp(mouse_coord, ll, ur);

                                            d = MIN(d, sqrt(ImLengthSqr(mouse_coord - pos)));
                                        }
                                        d += area_dist + layer_dist;
                                        break;
                                    }
                                    default:
                                        // Should not end up here
                                        break;
                                    }

                                    if (d < min_dist) {
                                        min_dist = d;
                                        char value_buf[64] = "";
                                        if (prop.print_value) {
                                            prop.print_value(value_buf, sizeof(value_buf), xi[1], &payload);
                                        } else {
                                            ImPlotPoint p = prop.getter[1](xn, &payload);
                                            snprintf(value_buf, sizeof(value_buf), "%.2f", p.y);
                                        }

                                        hovered_prop_idx = j;
                                        if (prop.hist.dim > 1) {
                                            hovered_pop_idx = k;
                                            snprintf(hovered_label, sizeof(hovered_label), "%s[%i]: %s", prop.label, k + 1, value_buf);
                                        } else {
                                            hovered_pop_idx = -1;
                                            snprintf(hovered_label, sizeof(hovered_label), "%s: %s", prop.label, value_buf);
                                        }
                                    }
                                }
                            }
                        }

                        if (hovered_prop_idx != -1) {
                            set_hovered_property(data, str_from_cstr(data->display_properties[hovered_prop_idx].label), hovered_pop_idx);
                            visualize_payload(data, data->display_properties[hovered_prop_idx].vis_payload, hovered_pop_idx, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);

                            if (strnlen(hovered_label, sizeof(hovered_label)) > 0) {
                                ImGui::SetTooltip("%s", hovered_label);
                            }
                        }
                    } else {
                        if (!str_empty(data->hovered_display_property_label)) {
                            for (size_t j = 0; j < md_array_size(data->display_properties); ++j) {
                                DisplayProperty& dp = data->display_properties[j];
                                if (dp.type != DisplayProperty::Type_Distribution) continue;
                                if (str_eq_cstr(data->hovered_display_property_label, dp.label)) {
                                    hovered_prop_idx = (int)j;
                                    hovered_pop_idx = data->hovered_display_property_pop_idx;
                                    break;
                                }
                            }
                        }
                    }
                    
                    for (int j = 0; j < num_props; ++j) {
                        DisplayProperty& dp = data->display_properties[j];
                        if (dp.type != DisplayProperty::Type_Distribution) continue;
                        if (!(dp.distribution_subplot_mask & (1 << i))) continue;

                        if (ImPlot::IsLegendEntryHovered(dp.label)) {
                            visualize_payload(data, dp.vis_payload, -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                            set_hovered_property(data, str_from_cstr(dp.label));
                            hovered_prop_idx = j;
                            hovered_pop_idx = -1;
                        }
                        
                        // legend context menu
                        if (ImPlot::BeginLegendPopup(dp.label)) {
                            if (ImGui::DeleteButton("Remove")) {
                                dp.distribution_subplot_mask &= ~(1 << i);
                                ImGui::CloseCurrentPopup();
                            }

                            const char* plot_type_names[] = {"Line", "Area", "Bars", "Scatter"};
                            const bool  valid_plot_types[] = {true, true, true, false};
                            STATIC_ASSERT(ARRAY_SIZE(plot_type_names) == DisplayProperty::PlotType_Count);
                            STATIC_ASSERT(ARRAY_SIZE(valid_plot_types) == DisplayProperty::PlotType_Count);

                            if (ImGui::BeginCombo("Plot Type", plot_type_names[dp.plot_type])) {
                                for (int k = 0; k < DisplayProperty::PlotType_Count; ++k) {
                                    if (!valid_plot_types[k]) continue;
                                    if (ImGui::Selectable(plot_type_names[k], dp.plot_type == k)) {
                                        dp.plot_type = (DisplayProperty::PlotType)k;
                                    }
                                }
                                ImGui::EndCombo();
                            }

                            if (dp.plot_type == DisplayProperty::PlotType_Bars) {
                                const double MIN_BAR_WIDTH = 0.01;
                                const double MAX_BAR_WIDTH = 1.00;
                                ImGui::SliderScalar("Bar Width", ImGuiDataType_Double, &dp.bar_width_scale, &MIN_BAR_WIDTH, &MAX_BAR_WIDTH, "%.3f");
                            }

                            if (dp.hist.dim > 1) {
                                const char* color_type_labels[] = {"Solid", "Colormap"};
                                STATIC_ASSERT(ARRAY_SIZE(color_type_labels) == DisplayProperty::ColorType_Count);

                                if (ImGui::BeginCombo("Color Type", color_type_labels[dp.color_type])) {
                                    for (int k = 0; k < DisplayProperty::ColorType_Count; ++k) {
                                        if (ImGui::Selectable(color_type_labels[k], k == dp.color_type)) {
                                            dp.color_type = (DisplayProperty::ColorType)k;
                                        }
                                    }
                                    ImGui::EndCombo();
                                }
                            }
                            switch (dp.color_type) {
                            case DisplayProperty::ColorType_Solid:
                                ImGui::ColorEdit4("Color", &dp.color.x);
                                break;
                            case DisplayProperty::ColorType_Colormap:
                                ImPlot::ColormapSelection("##Colormap", &dp.colormap);
                                ImGui::SliderFloat("Alpha", &dp.colormap_alpha, 0.0f, 1.0f);
                                break;
                            default:
                                ASSERT(false);
                                break;
                            } 

                            const int MIN_BINS = 32;
                            const int MAX_BINS = 1024;
                            if (ImGui::SliderInt("Num Bins", &dp.num_bins, MIN_BINS, MAX_BINS)) {
                                int next = next_power_of_two32(dp.num_bins);
                                int prev = next >> 1;
                                if (dp.num_bins - prev < next - dp.num_bins) {
                                    dp.num_bins = prev;
                                } else {
                                    dp.num_bins = next;
                                }
                                ImPlotPoint p_min = {dp.hist.x_min, dp.hist.y_min};
                                ImPlotPoint p_max = {dp.hist.x_max, dp.hist.y_max};
                                ImPlot::FitPoint(p_min);
                                ImPlot::FitPoint(p_max);
                            }

                            if (dp.hist.dim > 1) {
                                ImGui::Separator();
                                if (ImGui::Button("Set All")) {
                                    dp.population_mask = UINT64_MAX;
                                }
                                ImGui::SameLine();
                                if (ImGui::Button("Clear All")) {
                                    dp.population_mask = 0;
                                }

                                const float sz = ImGui::GetFontSize() * 1.5f;
                                ImGui::PushStyleVar(ImGuiStyleVar_SelectableTextAlign, ImVec2(0.5f, 0.5f));
                                for (int k = 0; k < MIN(dp.hist.dim, MAX_POPULATION_SIZE); ++k) {
                                    char lbl[32];
                                    snprintf(lbl, sizeof(lbl), "%d", k+1);
                                    if (ImGui::Selectable(lbl, dp.population_mask.test(k), ImGuiSelectableFlags_DontClosePopups, ImVec2(sz, sz))) {
                                        // Toggle bit for this population index
                                        dp.population_mask.flip(k);
                                    }
                                    if (ImGui::IsItemHovered()) {
                                        visualize_payload(data, dp.vis_payload, k, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                                        set_hovered_property(data, str_from_cstr(dp.label), k);
                                        hovered_prop_idx = j;
                                        hovered_pop_idx = k;
                                    }
                                    if (!k || ((k+1) % 10)) {
                                        ImGui::SameLine();
                                    }
                                }
                                ImGui::PopStyleVar();
                            }

                            ImPlot::EndLegendPopup();
                        }

                        auto plot = [j, &dp, hovered_prop_idx, hovered_pop_idx] (int k) {
                            const float  hov_fill_alpha  = 1.25f;
                            const float  hov_line_weight = 2.0f;
                            const float  hov_col_scl = 1.5f;
                            const int    population_size = CLAMP(dp.hist.dim, 1, MAX_POPULATION_SIZE);
                            const double bar_width = (dp.hist.x_max - dp.hist.x_min) / (dp.hist.num_bins);

                            ImVec4 color = {};
                            float  fill_alpha = 1.0f;
                            float  weight = 1.0f;

                            switch(dp.color_type) {
                            case DisplayProperty::ColorType_Solid:
                                color = dp.color;
                                break;
                            case DisplayProperty::ColorType_Colormap:
                                if (ImPlot::ColormapQualitative(dp.colormap)) {
                                    color = ImPlot::GetColormapColor(k, dp.colormap);
                                } else {
                                    color = ImPlot::SampleColormap( (float)k / (float)(population_size-1), dp.colormap);
                                }
                                color.w *= dp.colormap_alpha;
                                break;
                            default:
                                ASSERT(false);
                                break;
                            }

                            if (hovered_prop_idx == j) {
                                if (hovered_pop_idx == -1 || hovered_pop_idx == k) {
                                    color = ImVec4(ImSaturate(color.x * hov_col_scl), ImSaturate(color.y * hov_col_scl), ImSaturate(color.z * hov_col_scl), color.w);
                                    fill_alpha = hov_fill_alpha;
                                }
                                if (hovered_pop_idx == k) {
                                    weight = hov_line_weight;
                                }
                            }

                            DisplayProperty::Payload payload = {
                                .display_prop = &dp,
                                .dim_idx = k,
                            };

                            switch (dp.plot_type) {
                            case DisplayProperty::PlotType_Line:
                                ImPlot::SetNextLineStyle(color, weight);
                                ImPlot::PlotLineG(dp.label, dp.getter[1], &payload, dp.hist.num_bins);
                                break;
                            case DisplayProperty::PlotType_Area:
                                ImPlot::SetNextFillStyle(color, fill_alpha);
                                ImPlot::PlotShadedG(dp.label, dp.getter[0], &payload, dp.getter[1], &payload, dp.hist.num_bins);
                                break;
                            case DisplayProperty::PlotType_Bars:
                                ImPlot::SetNextFillStyle(color, fill_alpha);
                                ImPlot::PlotBarsG(dp.label, dp.getter[1], &payload, dp.hist.num_bins, bar_width * dp.bar_width_scale);
                                break;
                            default:
                                break;
                            }
                        };
                        
                        if (dp.hist.num_bins > 0) {
                            const int dim = CLAMP(dp.hist.dim, 1, MAX_POPULATION_SIZE);
                            for (int k = 0; k < dim; ++k) {
                                if (dp.hist.dim > 1 && !dp.population_mask.test(k)) {
                                    continue;
                                }
                                // Render this last, to make sure it is on top of the others
                                if (hovered_prop_idx == j && hovered_pop_idx == k) {
                                    continue;
                                }

                                plot(k);
                            }

                            if (hovered_prop_idx == j && hovered_pop_idx != -1) {
                                plot(hovered_pop_idx);
                            }
                        }

                        if (ImPlot::BeginDragDropSourceItem(dp.label)) {
                            DisplayPropertyDragDropPayload dnd_payload = {j, i};
                            ImGui::SetDragDropPayload("DISTRIBUTION_DND", &dnd_payload, sizeof(dnd_payload));
                            ImPlot::ItemIcon(dp.color); ImGui::SameLine();
                            ImGui::TextUnformatted(dp.label);
                            ImPlot::EndDragDropSource();
                        }
                    }

                    if (ImPlot::BeginDragDropTargetPlot()) {
                        if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload("DISTRIBUTION_DND")) {
                            ASSERT(payload->DataSize == sizeof(DisplayPropertyDragDropPayload));
                            DisplayPropertyDragDropPayload* dnd = (DisplayPropertyDragDropPayload*)(payload->Data);
                            data->display_properties[dnd->prop_idx].distribution_subplot_mask |= (1 << i);
                            if (dnd->src_plot_idx != -1 && dnd->src_plot_idx != i) {
                                // Clear bit from mask representing src plot index (only if it originated from another plot)
                                data->display_properties[dnd->prop_idx].distribution_subplot_mask &= ~(1 << dnd->src_plot_idx);
                            }

                            for (int j = 0; j < num_props; ++j) {
                                DisplayProperty& dp = data->display_properties[j];
                                if (dp.distribution_subplot_mask & (1 << i)) {
                                    ImPlotItem* item = ImPlot::GetItem(dp.label);
                                    bool just_dropped = dnd->prop_idx == j;
                                    bool previously_visible = item && item->Show;
                                    if (just_dropped || previously_visible) {
                                        ImPlot::GetCurrentPlot()->Axes[ImAxis_X1].ExtendFit(dp.hist.x_min);
                                        ImPlot::GetCurrentPlot()->Axes[ImAxis_X1].ExtendFit(dp.hist.x_max);
                                        ImPlot::GetCurrentPlot()->Axes[ImAxis_Y1].ExtendFit(dp.hist.y_min);
                                        ImPlot::GetCurrentPlot()->Axes[ImAxis_Y1].ExtendFit(dp.hist.y_max);
                                    }
                                }
                            }
                            
                            ImPlot::GetCurrentPlot()->Axes[ImAxis_X1].ApplyFit(ImPlot::GetStyle().FitPadding.x);
                            ImPlot::GetCurrentPlot()->Axes[ImAxis_Y1].ApplyFit(ImPlot::GetStyle().FitPadding.y);
                        }
                    }

                    if (ImPlot::IsPlotHovered()) {
                        if (hovered_prop_idx != -1) {
                            visualize_payload(data, data->display_properties[hovered_prop_idx].vis_payload, hovered_pop_idx, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                            set_hovered_property(data, str_from_cstr(data->display_properties[hovered_prop_idx].label), hovered_pop_idx);
                        }

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
            }
            ImPlot::EndSubplots();
        }

        // Try to handle the case when the user is dragging a payload and not dropping it within a valid target zone.
        // In such case if the property had a source plot index, remove the property from that plot
        if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
            const ImGuiPayload* payload = ImGui::GetDragDropPayload();
            if (payload && payload->IsDataType("DISTRIBUTION_DND") && !ImGui::IsDragDropPayloadBeingAccepted()) {
                DisplayPropertyDragDropPayload* dnd = (DisplayPropertyDragDropPayload*)(payload->Data);
                if (dnd && dnd->src_plot_idx != -1) {
                    data->display_properties[dnd->prop_idx].distribution_subplot_mask &= ~(1 << dnd->src_plot_idx);
                }
            }
        }
    }
    ImGui::End();
}

static void draw_density_volume_window(ApplicationState* data) {
    ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Density Volume", &data->density_volume.show_window, ImGuiWindowFlags_MenuBar)) {
        const ImVec2 button_size = {160, 0};
        bool volume_changed = false;

        if (ImGui::IsWindowFocused() && ImGui::IsKeyPressed(KEY_PLAY_PAUSE, false)) {
            data->animation.mode = data->animation.mode == PlaybackMode::Playing ? PlaybackMode::Stopped : PlaybackMode::Playing;
        }

        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("Property")) {
                int64_t selected_index = -1;
                int64_t candidate_count = 0;
                for (int64_t i = 0; i < (int64_t)md_array_size(data->display_properties); ++i) {
                    DisplayProperty& dp = data->display_properties[i];
                    if (dp.type != DisplayProperty::Type_Volume) continue;
                    if (!data->timeline.filter.enabled && dp.partial_evaluation) {
                        continue;
                    }
                    ImPlot::ItemIcon(dp.color); ImGui::SameLine();
                    if (ImGui::Selectable(dp.label, dp.show_in_volume)) {
                        selected_index = i;
                        volume_changed = true;
                    }
                    if (ImGui::IsItemHovered()) {
                        visualize_payload(data, dp.vis_payload, -1, MD_SCRIPT_VISUALIZE_DEFAULT);
                        set_hovered_property(data,  str_from_cstr(dp.label));
                    }
                    candidate_count += 1;
                }

                if (candidate_count == 0) {
                    ImGui::Text("No volume properties available.");
                }

                // Currently we only support viewing one volume at a time.
                // This will probably change over time but not now.
                if (selected_index != -1) {
                    for (int64_t i = 0; i < (int64_t)md_array_size(data->display_properties); ++i) {
                        if (selected_index == i) {
                            // Toggle bool
                            data->display_properties[i].show_in_volume = !data->display_properties[i].show_in_volume;
                        } else {
                            data->display_properties[i].show_in_volume = false;
                        }
                    }
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Render")) {
                ImGui::Checkbox("Direct Volume Rendering", &data->density_volume.dvr.enabled);
                if (data->density_volume.dvr.enabled) {
                    ImGui::Indent();
                    if (ImPlot::ColormapButton(ImPlot::GetColormapName(data->density_volume.dvr.tf.colormap), button_size, data->density_volume.dvr.tf.colormap)) {
                        ImGui::OpenPopup("Colormap Selector");
                    }
                    if (ImGui::BeginPopup("Colormap Selector")) {
                        for (int map = 4; map < ImPlot::GetColormapCount(); ++map) {
                            if (ImPlot::ColormapButton(ImPlot::GetColormapName(map), button_size, map)) {
                                data->density_volume.dvr.tf.colormap = map;
                                data->density_volume.dvr.tf.dirty = true;
                                ImGui::CloseCurrentPopup();
                            }
                        }
                        ImGui::EndPopup();
                    }
                    if (ImGui::SliderFloat("TF Alpha Scaling", &data->density_volume.dvr.tf.alpha_scale, 0.001f, 10.f, "%.3f", ImGuiSliderFlags_Logarithmic)) {
                        data->density_volume.dvr.tf.dirty = true;
                    }
                    ImGui::SliderFloat("TF Min Value", &data->density_volume.dvr.tf.min_val, 0.0f, 1000.f, "%.3f", ImGuiSliderFlags_Logarithmic);
                    ImGui::SameLine();
                    ImGui::SliderFloat("TF Max Value", &data->density_volume.dvr.tf.max_val, 0.0f, 1000.f, "%.3f", ImGuiSliderFlags_Logarithmic);
                    data->density_volume.dvr.tf.min_val = MIN(data->density_volume.dvr.tf.min_val, data->density_volume.dvr.tf.max_val);

                    ImGui::Unindent();
                }
                ImGui::Checkbox("Iso Surfaces", &data->density_volume.iso.enabled);
                if (data->density_volume.iso.enabled) {
                    ImGui::Indent();
                    for (int i = 0; i < data->density_volume.iso.count; ++i) {
                        ImGui::PushID(i);
                        ImGui::SliderFloat("##Isovalue", &data->density_volume.iso.values[i], 0.0f, 10.f, "%.3f", ImGuiSliderFlags_Logarithmic);
                        if (ImGui::IsItemDeactivatedAfterEdit()) {
                            // @TODO(Robin): Sort?
                        }
                        ImGui::SameLine();
                        ImGui::ColorEdit4Minimal("##Color", data->density_volume.iso.colors[i].elem);
                        ImGui::SameLine();
                        if (ImGui::DeleteButton(ICON_FA_XMARK)) {
                            for (int j = i; j < data->density_volume.iso.count - 1; ++j) {
                                data->density_volume.iso.colors[j] = data->density_volume.iso.colors[j+1];
                                data->density_volume.iso.values[j] = data->density_volume.iso.values[j+1];
                            }
                            data->density_volume.iso.count -= 1;
                        }
                        ImGui::PopID();
                    }
                    if ((data->density_volume.iso.count < (int)ARRAY_SIZE(data->density_volume.iso.values)) && ImGui::Button("Add", button_size)) {
                        size_t idx = data->density_volume.iso.count++;
                        data->density_volume.iso.values[idx] = 0.1f;
                        data->density_volume.iso.colors[idx] = { 0.2f, 0.1f, 0.9f, 1.0f };
                        // @TODO(Robin): Sort?
                    }
                        ImGui::SameLine();
                    if (ImGui::Button("Clear", button_size)) {
                        data->density_volume.iso.count = 0;
                    }
                    ImGui::Unindent();
                }
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Clip planes")) {
                ImGui::RangeSliderFloat("x", &data->density_volume.clip_volume.min.x, &data->density_volume.clip_volume.max.x, 0.0f, 1.0f);
                ImGui::RangeSliderFloat("y", &data->density_volume.clip_volume.min.y, &data->density_volume.clip_volume.max.y, 0.0f, 1.0f);
                ImGui::RangeSliderFloat("z", &data->density_volume.clip_volume.min.z, &data->density_volume.clip_volume.max.z, 0.0f, 1.0f);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Show")) {
                ImGui::Checkbox("Bounding Box", &data->density_volume.show_bounding_box);
                if (data->density_volume.show_bounding_box) {
                    ImGui::Indent();
                    ImGui::ColorEdit4("Color", data->density_volume.bounding_box_color.elem);
                    ImGui::Unindent();
                }
                ImGui::Checkbox("Reference Structure", &data->density_volume.show_reference_structures);
                if (data->density_volume.show_reference_structures) {
                    ImGui::Indent();
                    auto& rep = data->density_volume.rep;
                    ImGui::Checkbox("Show Superimposed Structures", &data->density_volume.show_reference_ensemble);

                    if (ImGui::BeginCombo("type", representation_type_str[(int)rep.type])) {
                        for (int i = 0; i < (int)RepresentationType::Orbital; ++i) {
                            if (ImGui::Selectable(representation_type_str[i], (int)rep.type == i)) {
                                rep.type = (RepresentationType)i;
                                data->density_volume.dirty_rep = true;
                            }
                        }
                        ImGui::EndCombo();
                    }

                    if (ImGui::BeginCombo("color", color_mapping_str[(int)rep.colormap])) {
                        for (int i = 0; i < (int)ColorMapping::Property; ++i) {
                            if (ImGui::Selectable(color_mapping_str[i], (int)rep.type == i)) {
                                rep.colormap = (ColorMapping)i;
                                data->density_volume.dirty_rep = true;
                            }
                        }
                        ImGui::EndCombo();
                    }

                    if (rep.colormap == ColorMapping::Uniform) {
                        data->density_volume.dirty_rep |= ImGui::ColorEdit4("color", rep.color.elem, ImGuiColorEditFlags_NoInputs);
                    }
                    if (rep.type == RepresentationType::SpaceFill || rep.type == RepresentationType::Licorice) {
                        data->density_volume.dirty_rep |= ImGui::SliderFloat("scale", &rep.param[0], 0.1f, 2.f);
                    }
                    if (rep.type == RepresentationType::Ribbons) {
                        data->density_volume.dirty_rep |= ImGui::SliderFloat("width", &rep.param[0], 0.1f, 2.f);
                        data->density_volume.dirty_rep |= ImGui::SliderFloat("thickness", &rep.param[1], 0.1f, 2.f);
                    }
                    if (rep.type == RepresentationType::Cartoon) {
                        data->density_volume.dirty_rep |= ImGui::SliderFloat("coil scale",  &rep.param[0], 0.1f, 3.f);
                        data->density_volume.dirty_rep |= ImGui::SliderFloat("sheet scale", &rep.param[1], 0.1f, 3.f);
                        data->density_volume.dirty_rep |= ImGui::SliderFloat("helix scale", &rep.param[2], 0.1f, 3.f);
                    }
                    ImGui::Unindent();
                }
                ImGui::Checkbox("Legend", &data->density_volume.legend.enabled);
                if (data->density_volume.legend.enabled) {
                    ImGui::Indent();
                    const char* colormap_modes[] = {"Opaque", "Transparent", "Split"};
                    if (ImGui::BeginCombo("Colormap", colormap_modes[data->density_volume.legend.colormap_mode])) {
                        for (int i = 0; i < IM_ARRAYSIZE(colormap_modes); ++i) {
                            if (ImGui::Selectable(colormap_modes[i])) {
                                data->density_volume.legend.colormap_mode = i;
                            }
                        }
                        ImGui::EndCombo();
                    }
                    ImGui::Checkbox("Use Checkerboard", &data->density_volume.legend.checkerboard);
                    if (ImGui::IsItemHovered()) {
                        ImGui::SetTooltip("Use a checkerboard background for transparent parts in the legend.");
                    }
                    ImGui::Unindent();
                }
                ImGui::Checkbox("Coordinate System Widget", &data->density_volume.show_coordinate_system_widget);
                ImGui::EndMenu();
            }

            ImGui::EndMenuBar();
        }

        update_density_volume(data);

        // Animate camera towards targets
        camera_animate(&data->density_volume.camera, data->density_volume.target_ori, data->density_volume.target_pos, data->density_volume.target_dist, data->app.timing.delta_s);

        /*
        // Canvas
        // Using InvisibleButton() as a convenience 1) it will advance the layout cursor and 2) allows us to use IsItemHovered()/IsItemActive()
        ImVec2 canvas_p0 = ImGui::GetCursorScreenPos();      // ImDrawList API uses screen coordinates!
        ImVec2 canvas_p1 = ImVec2(canvas_p0.x + canvas_sz.x, canvas_p0.y + canvas_sz.y);
        */
        ImVec2 canvas_sz = ImGui::GetContentRegionAvail();   // Resize canvas to what's available
        canvas_sz.x = MAX(canvas_sz.x, 50.0f);
        canvas_sz.y = MAX(canvas_sz.y, 50.0f);

        // This will catch our interactions
        ImGui::InvisibleButton("canvas", canvas_sz, ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_AllowOverlap);

        // Draw border and background color
        ImGuiIO& io = ImGui::GetIO();

        ImVec2 canvas_p0 = ImGui::GetItemRectMin();
        ImVec2 canvas_p1 = ImGui::GetItemRectMax();

        ImDrawList* draw_list = ImGui::GetWindowDrawList();
        draw_list->AddImage((ImTextureID)(intptr_t)data->density_volume.fbo.tex.transparency, canvas_p0, canvas_p1, { 0,1 }, { 1,0 });
        draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));

        if (data->density_volume.dvr.enabled && data->density_volume.legend.enabled) {
            ImVec2 canvas_ext = canvas_p1 - canvas_p0;
            ImVec2 cmap_ext = {MIN(canvas_ext.x * 0.5f, 250.0f), MIN(canvas_ext.y * 0.25f, 30.0f)};
            ImVec2 cmap_pad = {10, 10};
            ImVec2 cmap_pos = canvas_p1 - ImVec2(cmap_ext.x, cmap_ext.y) - cmap_pad;
            ImPlotColormap cmap = data->density_volume.dvr.tf.colormap;
            ImPlotContext& gp = *ImPlot::GetCurrentContext();
            ImU32 checker_bg = IM_COL32(255, 255, 255, 255);
            ImU32 checker_fg = IM_COL32(128, 128, 128, 255);
            float checker_size = 8.0f;
            ImVec2 checker_offset = ImVec2(0,0);

            int mode = data->density_volume.legend.colormap_mode;

            ImVec2 opaque_scl = ImVec2(1,1);
            ImVec2 transp_scl = ImVec2(0,0);

            if (mode == LegendColorMapMode_Split) {
                opaque_scl = ImVec2(1.0f, 0.5f);
                transp_scl = ImVec2(0.0f, 0.5f);
            }

            ImRect opaque_rect = ImRect(cmap_pos, cmap_pos + cmap_ext * opaque_scl);
            ImRect transp_rect = ImRect(cmap_pos + cmap_ext * transp_scl, cmap_pos + cmap_ext);
            
            // Opaque
            if (mode == LegendColorMapMode_Opaque || mode == LegendColorMapMode_Split) {
                ImPlot::RenderColorBar(gp.ColormapData.GetKeys(cmap),gp.ColormapData.GetKeyCount(cmap),*draw_list,opaque_rect,false,false,!gp.ColormapData.IsQual(cmap));
            }
            
            if (mode == LegendColorMapMode_Transparent || mode == LegendColorMapMode_Split) {
                if (data->density_volume.legend.checkerboard) {
                    // Checkerboard
                    ImGui::DrawCheckerboard(draw_list, transp_rect.Min, transp_rect.Max, checker_bg, checker_fg, checker_size, checker_offset);
                }
                // Transparent
                draw_list->AddImage((ImTextureID)(intptr_t)data->density_volume.dvr.tf.id, transp_rect.Min, transp_rect.Max);
            }
            
            // Boarder
            draw_list->AddRect(cmap_pos, cmap_pos + cmap_ext, IM_COL32(0, 0, 0, 255));
        }

        const bool is_hovered = ImGui::IsItemHovered();
        const bool is_active = ImGui::IsItemActive();
        const ImVec2 origin(canvas_p0.x, canvas_p0.y);  // Lock scrolled origin
        const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

        auto& gbuf = data->density_volume.fbo;
        int width  = MAX(1, (int)canvas_sz.x);
        int height = MAX(1, (int)canvas_sz.y);
        if ((int)gbuf.width != width || (int)gbuf.height != height) {
            init_gbuffer(&gbuf, width, height);
        }

        bool reset_hard = false;
        if (volume_changed) {
            static bool first_time = true;
            if (first_time) {
                reset_hard = true;
                first_time = false;
            }
        }
        bool reset_view = reset_hard;
        if (is_hovered) {
            if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                reset_view = true;
            }
        }

        if (reset_view) {
            vec3_t aabb_min = vec3_from_vec4(data->density_volume.model_mat * vec4_set(0,0,0,1));
            vec3_t aabb_max = vec3_from_vec4(data->density_volume.model_mat * vec4_set(1,1,1,1));
            camera_compute_optimal_view(&data->density_volume.target_pos, &data->density_volume.target_ori, &data->density_volume.target_dist, aabb_min, aabb_max);

            if (reset_hard) {
                data->density_volume.camera.position = data->density_volume.target_pos;
                data->density_volume.camera.orientation = data->density_volume.target_ori;
                data->density_volume.camera.focus_distance = data->density_volume.target_dist;
            }
        }

        if (is_active || is_hovered) {
            static const TrackballControllerParam param = {
                .min_distance = 1.0,
                .max_distance = 1000.0,
            };

            vec2_t delta = { io.MouseDelta.x, io.MouseDelta.y };
            vec2_t curr = {mouse_pos_in_canvas.x, mouse_pos_in_canvas.y};
            vec2_t prev = curr - delta;
            float  wheel_delta = io.MouseWheel;

            TrackballControllerInput input = {
                .rotate_button = is_active && ImGui::IsMouseDown(ImGuiMouseButton_Left),
                .pan_button    = is_active && ImGui::IsMouseDown(ImGuiMouseButton_Right),
                .dolly_button  = is_active && ImGui::IsMouseDown(ImGuiMouseButton_Middle),
                .dolly_delta   = is_hovered ? wheel_delta : 0.0f,
                .mouse_coord_prev = prev,
                .mouse_coord_curr = curr,
                .screen_size = {canvas_sz.x, canvas_sz.y},
                .fov_y = data->density_volume.camera.fov_y,
            };
            camera_controller_trackball(&data->density_volume.target_pos, &data->density_volume.target_ori, &data->density_volume.target_dist, input, param);
        }

        if (data->density_volume.show_coordinate_system_widget) {
            ImVec2 win_size = ImGui::GetWindowSize();
            float  ext = MIN(win_size.x, win_size.y) * 0.2f;
            float  pad = 0.1f * ext;

            CoordSystemWidgetParam params = {
                .pos = ImVec2(pad, win_size.y - ext - pad),
                .size = {ext, ext},
                .view_matrix = camera_world_to_view_matrix(data->density_volume.camera),
                .camera_ori  = data->density_volume.target_ori,
                .camera_pos  = data->density_volume.target_pos,
                .camera_dist = data->density_volume.target_dist,
            };

            ImGui::DrawCoordinateSystemWidget(params);
        }

        mat4_t view_mat = camera_world_to_view_matrix(data->density_volume.camera);
        mat4_t proj_mat = camera_perspective_projection_matrix(data->density_volume.camera, (float)canvas_sz.x / (float)canvas_sz.y);
        mat4_t inv_proj_mat = camera_inverse_perspective_projection_matrix(data->density_volume.camera, (float)canvas_sz.x / (float)canvas_sz.y);

        PUSH_GPU_SECTION("RENDER DENSITY VOLUME");
        clear_gbuffer(&gbuf);

        const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
            GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_TRANSPARENCY };

        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);

        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);
        glDepthFunc(GL_LESS);

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
        glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
        glViewport(0, 0, gbuf.width, gbuf.height);
        glScissor(0, 0,  gbuf.width, gbuf.height);

        int64_t selected_property = -1;
        for (size_t i = 0; i < md_array_size(data->display_properties); ++i) {
            const DisplayProperty& dp = data->display_properties[i];
            if (dp.type == DisplayProperty::Type_Volume && dp.show_in_volume) {
                selected_property = i;
                break;
            }
        }

        size_t num_reps = md_array_size(data->density_volume.gl_reps);
        if (selected_property > -1 && data->density_volume.show_reference_structures && num_reps > 0) {
            if (!data->density_volume.show_reference_ensemble) {
                num_reps = 1;
            }

            md_gl_draw_op_t* draw_ops = 0;

            md_gl_draw_op_t op = {};
            op.type = (md_gl_rep_type_t)data->density_volume.rep.type;
            MEMCPY(&op.args, data->density_volume.rep.param, sizeof(op.args));

            for (size_t i = 0; i < num_reps; ++i) {
                op.rep = data->density_volume.gl_reps[i];
                op.model_matrix = &data->density_volume.rep_model_mats[i].elem[0][0];
                md_array_push(draw_ops, op, frame_alloc);
            }

            md_gl_draw_args_t draw_args = {
                .shaders = data->mold.gl_shaders,
                .draw_operations = {
                    .count = (uint32_t)md_array_size(draw_ops),
                    .ops = draw_ops
                },
                .view_transform = {
                    .view_matrix = &view_mat.elem[0][0],
                    .proj_matrix = &proj_mat.elem[0][0],
                },
            };

            md_gl_draw(&draw_args);

            if (is_hovered) {
                const vec2_t coord = {mouse_pos_in_canvas.x, (float)gbuf.height - mouse_pos_in_canvas.y};
                uint32_t picking_idx = INVALID_PICKING_IDX;
                extract_picking_data(&picking_idx, NULL, &gbuf, (int)coord.x, (int)coord.y);
                if (picking_idx != INVALID_PICKING_IDX) {
                    draw_info_window(*data, picking_idx);
                }
            }
            glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
        }

        if (data->density_volume.show_density_volume) {
            if (data->density_volume.model_mat != mat4_t{ 0 }) {
                volume::RenderDesc vol_desc = {
                    .render_target = {
                        .depth  = gbuf.tex.depth,
                        .color  = gbuf.tex.transparency,
                        .width  = gbuf.width,
                        .height = gbuf.height,
                    },
                    .texture = {
                        .volume = data->density_volume.volume_texture.id,
                        .tf = data->density_volume.dvr.tf.id,
                    },
                    .matrix = {
                        .model = data->density_volume.model_mat,
                        .view = view_mat,
                        .proj = proj_mat,
                        .inv_proj = inv_proj_mat,
                    },
                    .clip_volume = {
                        .min = data->density_volume.clip_volume.min,
                        .max = data->density_volume.clip_volume.max,
                    },
                    .iso = {
                        .enabled = data->density_volume.iso.enabled,
                        .count = data->density_volume.iso.count,
                        .values = data->density_volume.iso.values,
                        .colors = data->density_volume.iso.colors,
                    },
                    .dvr = {
                        .enabled = data->density_volume.dvr.enabled,
                    },
                    .tf = {
                        .min_value = data->density_volume.dvr.tf.min_val,
                        .max_value = data->density_volume.dvr.tf.max_val,
                    },
                    .shading = {
                        .env_radiance = data->visuals.background.color * data->visuals.background.intensity * 0.25f,
                        .roughness = 0.3f,
                        .dir_radiance = {10,10,10},
                        .ior = 1.5f,
                    },
                    .voxel_spacing = data->density_volume.voxel_spacing
                    
                };
                volume::render_volume(vol_desc);
            }
        }

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
        glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
        glViewport(0, 0, gbuf.width, gbuf.height);
        glScissor(0, 0, gbuf.width, gbuf.height);

        if (data->density_volume.show_bounding_box) {
            glEnable(GL_DEPTH_TEST);
            glDepthMask(GL_TRUE);
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

            if (data->density_volume.model_mat != mat4_t{0}) {
                immediate::set_model_view_matrix(mat4_mul(view_mat, data->density_volume.model_mat));
                immediate::set_proj_matrix(proj_mat);

                uint32_t box_color = convert_color(data->density_volume.bounding_box_color);
                uint32_t clip_color = convert_color(data->density_volume.clip_volume_color);
                immediate::draw_box_wireframe({0,0,0}, {1,1,1}, box_color);
                immediate::draw_box_wireframe(data->density_volume.clip_volume.min, data->density_volume.clip_volume.max, clip_color);

                immediate::render();
            }
            glDisable(GL_BLEND);
        }

        PUSH_GPU_SECTION("Postprocessing")
        postprocessing::Descriptor postprocess_desc = {
            .background = {
                .color = data->visuals.background.color * data->visuals.background.intensity,
            },
            .tonemapping = {
                .enabled = data->visuals.tonemapping.enabled,
                .mode = data->visuals.tonemapping.tonemapper,
                .exposure = data->visuals.tonemapping.exposure,
                .gamma = data->visuals.tonemapping.gamma,
            },
            .ambient_occlusion = {
                .enabled = false,
            },
            .depth_of_field = {
                .enabled = false,
            },
            .fxaa = {
                .enabled = true,
            },
            .temporal_aa = {
                .enabled = false,
            },
            .sharpen = {
                .enabled = false,
            },
            .input_textures = {
                .depth = gbuf.tex.depth,
                .color = gbuf.tex.color,
                .normal = gbuf.tex.normal,
                .velocity = gbuf.tex.velocity,
                .transparency = gbuf.tex.transparency,
            }
        };

        ViewParam view_param = {
            .matrix = {
                .curr = {
                    .view = view_mat,
                    .proj = proj_mat,
                    .norm = view_mat,
                },
                .inv = {
                    .proj = inv_proj_mat,
                }
            },
            .clip_planes = {
                .near = data->density_volume.camera.near_plane,
                .far = data->density_volume.camera.far_plane,
            },
            .resolution = {canvas_sz.x, canvas_sz.y},
            .fov_y = data->density_volume.camera.fov_y,
        };

        postprocessing::shade_and_postprocess(postprocess_desc, view_param);
        POP_GPU_SECTION()

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glDrawBuffer(GL_BACK);

        POP_GPU_SECTION();
    }

    ImGui::End();
}

static void draw_dataset_window(ApplicationState* data) {
    ASSERT(data);

    ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Dataset", &data->dataset.show_window)) {
        ImGui::Text("Molecular data: %s", data->files.molecule);
        ImGui::Text("Num atoms:    %9i", (int)data->mold.mol.atom.count);
        ImGui::Text("Num residues: %9i", (int)data->mold.mol.residue.count);
        ImGui::Text("Num chains:   %9i", (int)data->mold.mol.chain.count);

        if (data->files.trajectory[0] != '\0') {
            ImGui::Separator();
            ImGui::Text("Trajectory data: %s", data->files.trajectory);
            ImGui::Text("Num frames:    %9i", (int)md_trajectory_num_frames(data->mold.traj));
            ImGui::Text("Num atoms:     %9i", (int)md_trajectory_num_atoms(data->mold.traj));
        }

        const char* lbls[] = {"Chains", "Residue Names", "Atom Labels"};
        const md_array(DatasetItem) items[] = {data->dataset.chains, data->dataset.residue_names, data->dataset.atom_types};
        STATIC_ASSERT(ARRAY_SIZE(lbls) == ARRAY_SIZE(items));

        const ImVec2 item_size = ImVec2(ImGui::GetFontSize() * 1.8f, ImGui::GetFontSize() * 1.1f);
        const float window_x_max = ImGui::GetWindowPos().x + ImGui::GetWindowContentRegionMax().x;
        for (size_t i = 0; i < ARRAY_SIZE(lbls); ++i) {
            const size_t count = md_array_size(items[i]);
            if (count) {
                if (ImGui::CollapsingHeader(lbls[i])) {
                    bool item_hovered = false;
                    for (size_t j = 0; j < count; ++j) {
                        const DatasetItem& item = items[i][j];
                        const float t = powf(item.fraction, 0.2f) * 0.5f;

                        ImGui::PushStyleColor(ImGuiCol_HeaderHovered, ImVec4(1, 1, 0.5, 0.3));
                        ImGui::PushStyleColor(ImGuiCol_Header, ImPlot::SampleColormap(t, ImPlotColormap_Plasma));
                        ImGui::Selectable(item.label, true, 0, item_size);
                        ImGui::PopStyleColor(2);
                        //We do not show an item as selected in the UI, as we don't keep track if the whole item group is selected
                        //Selecting will thus mark the atoms as selected, but not the item. It's a one way selection

                        if (ImGui::IsItemHovered()) {
                            ImGui::SetTooltip("%s: count %d (%.2f%%)", item.label, item.count, item.fraction * 100.f);
                            filter_expression(data, str_from_cstr(item.query), &data->selection.highlight_mask);
                            item_hovered = true;

                            //Select
                            if (ImGui::IsKeyDown(ImGuiKey_MouseLeft) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                                md_bitfield_or_inplace(&data->selection.selection_mask, &data->selection.highlight_mask);
                            }
                            //Deselect
                            else if (ImGui::IsKeyDown(ImGuiKey_MouseRight) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                                md_bitfield_andnot_inplace(&data->selection.selection_mask, &data->selection.highlight_mask);
                            }
                        }

                        float last_item_x = ImGui::GetItemRectMax().x;
                        float next_button_x = last_item_x + item_size.x;
                        if (j + 1 < count && next_button_x < window_x_max) {
                            ImGui::SameLine();
                        }

                        if (!item_hovered && ImGui::IsWindowHovered()) {
                            //Makes sure that we clear the highlight if we are in this window, but don't hover an item
                            md_bitfield_clear(&data->selection.highlight_mask);
                        }
                    }
                }
            }
        }

        const size_t num_mappings = md_array_size(data->dataset.atom_element_remappings);
        if (num_mappings) {
            if (ImGui::CollapsingHeader("Atom Element Mappings")) {
                for (size_t i = 0; i < num_mappings; ++i) {
                    const auto& mapping = data->dataset.atom_element_remappings[i];
                    ImGui::Text("%s -> %s (%s)", mapping.lbl, md_util_element_name(mapping.elem).ptr, md_util_element_symbol(mapping.elem).ptr);
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
        size_t sema_count = 0;
        if (md_semaphore_query_count(&data->script.ir_semaphore, &sema_count)) {
            ImGui::Text("Script IR semaphore count: %i", (int)sema_count);
        }
        
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

        if (ImGui::Button("Load Data")) {
            char buf[2048];
            if (application::file_dialog(buf, sizeof(buf), application::FileDialogFlag_Open, STR_LIT("out"))) {
                struct Payload {
                    ApplicationState* state;
                    str_t path;
                };
                Payload payload = {
                    .state = data,
                    .path  = str_from_cstr(buf),
                };
                viamd::event_system_broadcast_event(HASH_STR_LIT("Secret Sauce"), 0, &payload);
                update_representation_info(data);
                update_all_representations(data);
            }
        }
    }
    ImGui::End();
}

static void draw_script_editor_window(ApplicationState* data) {
    ASSERT(data);

    ImGui::SetNextWindowSize({300,200}, ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Script Editor", &data->show_script_window, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_MenuBar)) {
        ImGui::SetWindowSize(ImVec2(800, 600), ImGuiCond_FirstUseEver);
        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("File")) {
                char path_buf[1024] = "";
                if (ImGui::MenuItem("Load")) {
                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Open, STR_LIT("txt"))) {
                        str_t txt = load_textfile(str_from_cstr(path_buf), frame_alloc);
                        std::string str(txt.ptr, txt.len);
                        editor.SetText(str);
                    }
                }
                if (ImGui::MenuItem("Save")) {
                    auto textToSave = editor.GetText();
                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, STR_LIT("txt"))) {
                        str_t path = str_t{path_buf, strnlen(path_buf, sizeof(path_buf))};
                        md_file_o* file = md_file_open(path, MD_FILE_WRITE);
                        if (file) {
                            md_file_write(file, textToSave.c_str(), textToSave.length());
                            md_file_close(file);
                        } else {
                            LOG_ERROR("Failed to open file '%s' for saving script", path_buf);
                        }
                    }
                }
                if (ImGui::MenuItem("Export")) {
                    data->show_property_export_window = true;
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Edit")) {
                bool ro = editor.IsReadOnly();
                if (ImGui::MenuItem("Read-only mode", nullptr, &ro))
                    editor.SetReadOnly(ro);
                ImGui::Separator();

                if (ImGui::MenuItem("Undo", "ALT-Backspace", nullptr, !ro && editor.CanUndo()))
                    editor.Undo();
                if (ImGui::MenuItem("Redo", "Ctrl-Y", nullptr, !ro && editor.CanRedo()))
                    editor.Redo();

                ImGui::Separator();

                if (ImGui::MenuItem("Copy", "Ctrl-C", nullptr, editor.HasSelection()))
                    editor.Copy();
                if (ImGui::MenuItem("Cut", "Ctrl-X", nullptr, !ro && editor.HasSelection()))
                    editor.Cut();
                if (ImGui::MenuItem("Delete", "Del", nullptr, !ro && editor.HasSelection()))
                    editor.Delete();
                if (ImGui::MenuItem("Paste", "Ctrl-V", nullptr, !ro && ImGui::GetClipboardText() != nullptr))
                    editor.Paste();

                ImGui::Separator();

                if (ImGui::MenuItem("Select all", nullptr, nullptr))
                    editor.SetSelection(TextEditor::Coordinates(), TextEditor::Coordinates(editor.GetTotalLines(), 0));

                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Settings")) {
                if (ImGui::MenuItem("Dark palette"))
                    editor.SetPalette(TextEditor::GetDarkPalette());
                if (ImGui::MenuItem("Light palette"))
                    editor.SetPalette(TextEditor::GetLightPalette());
                if (ImGui::MenuItem("Retro blue palette"))
                    editor.SetPalette(TextEditor::GetRetroBluePalette());
                ImGui::Separator();
                ImGui::ColorEdit4("Point Color",      data->script.point_color.elem);
                ImGui::ColorEdit4("Line Color",       data->script.line_color.elem);
                ImGui::ColorEdit4("Triangle Color",   data->script.triangle_color.elem);

                ImGui::EndMenu();
            }

            ImGui::EndMenuBar();
        }

        if (editor.IsTextChanged()) {
            data->script.compile_ir = true;
            data->script.time_since_last_change = 0;
        }

        const ImVec2 content_size = ImGui::GetContentRegionAvail();
        const char* btn_text = "Evaluate";
        const ImVec2 label_size = ImGui::CalcTextSize(btn_text, NULL, true) * ImVec2(1.4, 1.0);
        const ImVec2 btn_size = ImGui::CalcItemSize(ImVec2(0,0), label_size.x + ImGui::GetStyle().FramePadding.x * 2.0f, label_size.y + ImGui::GetStyle().FramePadding.y * 2.0f);
        const ImVec2 text_size(content_size - ImVec2(0, btn_size.y + ImGui::GetStyle().ItemSpacing.y));

        editor.Render("TextEditor", text_size);
        bool editor_hovered = ImGui::IsItemHovered();
        bool eval = false;
        if (editor.IsFocused() && ImGui::IsKeyDown(KEY_SCRIPT_EVALUATE_MOD) && ImGui::IsKeyPressed(KEY_SCRIPT_EVALUATE)) {
            eval = true;
        }

        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + content_size.x - btn_size.x);
        ImGui::SetCursorPosY(ImGui::GetCursorPosY() + ImGui::GetStyle().ItemSpacing.y);

        const bool valid = md_script_ir_valid(data->script.ir);

        if (!valid) ImGui::PushDisabled();
        if (ImGui::Button(btn_text, btn_size)) {
            eval = true;
        }
        if (!valid) ImGui::PopDisabled();

        if (eval && valid) {
            data->script.eval_init = true;
        }

        if (editor_hovered) {
            md_bitfield_clear(&data->selection.highlight_mask);
        }

        const TextEditor::Marker* hovered_marker = editor.GetHoveredMarker();
        if (hovered_marker) {
            if (!hovered_marker->text.empty()) {
                ImGui::BeginTooltip();
                ImGui::Text("%.*s", (int)hovered_marker->text.length(), hovered_marker->text.c_str());
                if (data->script.sub_idx != -1) {
                    ImGui::Text("Currently inspected idx: %i", data->script.sub_idx + 1);
                }
                ImGui::EndTooltip();
            }
            if (hovered_marker->payload) {
                if (md_semaphore_try_aquire(&data->script.ir_semaphore)) {
                    defer { md_semaphore_release(&data->script.ir_semaphore); };

                    if (hovered_marker->type == MarkerType_Error || hovered_marker->type == MarkerType_Warning) {
                        const md_bitfield_t* bf = (const md_bitfield_t*)hovered_marker->payload;
                        md_bitfield_copy(&data->selection.highlight_mask, bf);
                    }
                    else if (hovered_marker->type == MarkerType_Visualization) {
                        // Clear hovered property
                        set_hovered_property(data, STR_LIT(""));

                        if (md_script_ir_valid(data->script.ir)) {
                            const md_script_vis_payload_o* payload = (const md_script_vis_payload_o*)hovered_marker->payload;
                            str_t payload_ident = md_script_payload_ident(payload);
                            size_t payload_dim  = md_script_payload_dim(payload); 

                            if (payload_dim > 1) {
                                int delta = (int)ImGui::GetIO().MouseWheel;
                                if (ImGui::IsKeyDown(ImGuiMod_Shift)) {
                                    delta *= 10;
                                }
                                data->script.sub_idx += delta;
                                data->script.sub_idx = CLAMP(data->script.sub_idx, -1, (int)payload_dim - 1);
                            }

                            visualize_payload(data, payload, data->script.sub_idx, 0);
                            set_hovered_property(data, payload_ident, data->script.sub_idx);
                        }
                    }

                    bool lm_click = ImGui::IsMouseClicked(ImGuiMouseButton_Left);
                    bool rm_click = ImGui::IsMouseClicked(ImGuiMouseButton_Right);
                    if (ImGui::IsKeyDown(ImGuiMod_Shift) && (lm_click || rm_click)) {
                        SelectionOperator op = lm_click ? SelectionOperator::Or : SelectionOperator::AndNot;
                        modify_selection(data, &data->selection.highlight_mask, op);
                    }
                }
            }
        } else {
            data->script.sub_idx = -1;
        }
    }
    ImGui::End();
}

static bool export_xvg(const float* column_data[], const char* column_labels[], size_t num_columns, size_t num_rows, str_t filename) {
    ASSERT(column_data);
    ASSERT(column_labels);
    ASSERT(num_columns >= 0);
    ASSERT(num_rows >= 0);

    md_file_o* file = md_file_open(filename, MD_FILE_WRITE);
    if (!file) {
        LOG_ERROR("Failed to open file '" STR_FMT "' to write data.", STR_ARG(filename));
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

    md_file_close(file);
    LOG_SUCCESS("Successfully exported XVG file to '%.*s'", (int)filename.len, filename.ptr);
    return true;
}

static bool export_csv(const float* column_data[], const char* column_labels[], size_t num_columns, size_t num_rows, str_t filename) {
    ASSERT(column_data);
    ASSERT(column_labels);
    ASSERT(num_columns >= 0);
    ASSERT(num_rows >= 0);

    md_file_o* file = md_file_open(filename, MD_FILE_WRITE);
    if (!file) {
        LOG_ERROR("Failed to open file '%.*s' to write data.", (int)filename.len, filename.ptr);
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

    md_file_close(file);
    LOG_SUCCESS("Successfully exported CSV file to '%.*s'", (int)filename.len, filename.ptr);
    return true;
}

static bool export_cube(const ApplicationState& data, const md_script_property_data_t* prop_data, const md_script_vis_payload_o* vis_payload, str_t filename) {
    // @NOTE: First we need to extract some meta data for the cube format, we need the atom indices/bits for any SDF
    // And the origin + extent of the volume in spatial coordinates (Ångström)

    if (!prop_data) {
        LOG_ERROR("Export Cube: The property to be exported did not exist");
        return false;
    }

    if (!vis_payload) {
        LOG_ERROR("Export Cube: Missing input visualization data");
        return false;
    }

    // Copy mol and replace with initial coords
    md_molecule_t mol = data.mold.mol;

    size_t stride = ALIGN_TO(data.mold.mol.atom.count, 8);
    float* coords = (float*)md_alloc(frame_alloc, stride * sizeof(float) * 3);
    mol.atom.x = coords + stride * 0;
    mol.atom.y = coords + stride * 1;
    mol.atom.z = coords + stride * 2;
    
    if (!md_trajectory_load_frame(data.mold.traj, 0, NULL, mol.atom.x, mol.atom.y, mol.atom.z)) {
        return false;
    }

    bool result = false;
    md_script_vis_t vis = { 0 };
    md_script_vis_init(&vis, frame_alloc);
    defer { md_script_vis_free(&vis); };
    
    //if (md_semaphore_aquire(&data.script.ir_semaphore)) {
    //    defer { md_semaphore_release(&data.script.ir_semaphore); };

        
    if (md_script_ir_valid(data.script.eval_ir)) {
        md_script_vis_ctx_t ctx = {
            .ir = data.script.eval_ir,
            .mol = &data.mold.mol,
            .traj = data.mold.traj,
        };
        result = md_script_vis_eval_payload(&vis, vis_payload, 0, &ctx, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_SDF);
    }

    if (result == true) {
        md_file_o* file = md_file_open(filename, MD_FILE_WRITE);
        if (!file) {
            LOG_ERROR("Failed to open file '%.*s' in order to write to it.", (int)filename.len, filename.ptr);
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
            const int vol_dim[3] = {prop_data->dim[1], prop_data->dim[2], prop_data->dim[3]};
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

            int64_t beg_bit = bf->beg_bit;
            int64_t end_bit = bf->end_bit;
            while ((beg_bit = md_bitfield_scan(bf, beg_bit, end_bit)) != 0) {
                int64_t i = beg_bit - 1;
                vec3_t coord = {mol.atom.x[i], mol.atom.y[i], mol.atom.z[i]};
                coord = mat4_mul_vec3(M, coord, 1.0f);
                // @NOTE(Robin): If we don't have any elements available for example in the case of coarse grained, we use a placeholder of 1 (Hydrogen).
                md_element_t elem = mol.atom.element ? mol.atom.element[i] : 1;
                md_file_printf(file, "%5i %12.6f %12.6f %12.6f %12.6f\n", elem, (float)elem, coord.x, coord.y, coord.z);
            }

            // This entry somehow relates to the number of densities
            md_file_printf(file, "%5i %5i\n", 1, 1);

            // Write density data
            int count = 0;
            for (int x = 0; x < vol_dim[0]; ++x) {
                for (int y = 0; y < vol_dim[1]; ++y) {
                    for (int z = 0; z < vol_dim[2]; ++z) {
                        int idx = z * vol_dim[0] * vol_dim[1] + y * vol_dim[0] + x;
                        float val = prop_data->values[idx];
                        md_file_printf(file, " %12.6E", val);
                        if (++count % 6 == 0) md_file_printf(file, "\n");
                    }
                }
            }
        }

        md_file_close(file);
    } else {
        LOG_ERROR("Failed to visualize volume for export.");
        return false;
    }

    return true;
}

#define APPEND_BUF(buf, len, fmt, ...) (len += snprintf(buf + len, MAX(0, (int)sizeof(buf) - len), fmt, ##__VA_ARGS__) + 1)

static md_array(float) sample_range(float beg, float end, int sample_count, md_allocator_i* alloc) {
    md_array(float) result = md_array_create(float, sample_count, alloc);
    double step = (end - beg) / (double)(sample_count - 1);
    for (int i = 0; i < sample_count; ++i) {
        result[i] = (float)((double)beg + step * (double)i);
    }
    return result;
}

static void draw_property_export_window(ApplicationState* data) {
    ASSERT(data);

    struct ExportFormat {
        str_t lbl;
        str_t ext;
    };

    ExportFormat table_formats[] {
        {STR_LIT("XVG"), STR_LIT("xvg")},
        {STR_LIT("CSV"), STR_LIT("csv")}
    };

    ExportFormat volume_formats[] {
        {STR_LIT("Gaussian Cube"), STR_LIT("cube")},
    };

    if (ImGui::Begin("Property Export", &data->show_property_export_window)) {
        static int type = DisplayProperty::Type_Temporal;
        static int property_idx  = 0;
        static int table_format  = 0;
        static int volume_format = 0;

        int num_properties = (int)md_array_size(data->display_properties);
        if (num_properties == 0) {
            ImGui::Text("No properties available for export, try evaluating the script.");
            ImGui::End();
            property_idx = 0;
            return;
        }

        if (task_system::task_is_running(data->tasks.evaluate_full)) {
            ImGui::Text("The properties are currently being evaluated, please wait...");
            ImGui::End();
            property_idx = 0;
            return;
        }

        ImGui::PushItemWidth(200);
        if (ImGui::Combo("Data Type", (int*)(&type), "Temporal\0Distribution\0Density Volume\0")) {
            if (property_idx != -1) {
                const char* cur_lbl = data->display_properties[property_idx].label;
                for (int i = 0; i < (int)md_array_size(data->display_properties); ++i) {
                    const DisplayProperty& dp = data->display_properties[i];
                    if (type == dp.type && strcmp(cur_lbl, dp.label) == 0) {
                        property_idx = i;
                        break;
                    }
                }
            }
        }
        ImGui::Separator();
        
        if (ImGui::BeginCombo("Property", property_idx != -1 ? data->display_properties[property_idx].label : "")) {
            for (int i = 0; i < (int)md_array_size(data->display_properties); ++i) {
                const DisplayProperty& dp = data->display_properties[i];
                if (dp.type == DisplayProperty::Type_Distribution && dp.hist.num_bins == 0) {
                    continue; 
                }
                if (type == dp.type) {
                    if (ImGui::Selectable(dp.label, property_idx == i)) {
                        property_idx = i;
                    }
                }
            }
            ImGui::EndCombo();
        }
        
        str_t file_extension = {};
        if (type == DisplayProperty::Type_Distribution || type == DisplayProperty::Type_Temporal) {
            if (ImGui::BeginCombo("File Format", table_formats[table_format].lbl.ptr)) {
                for (int i = 0; i < (int)ARRAY_SIZE(table_formats); ++i) {
                    if (ImGui::Selectable(table_formats[i].lbl.ptr, table_format == i)) {
                        table_format = i;
                    }
                }
                ImGui::EndCombo();
            }
            file_extension = table_formats[table_format].ext;
        } else if (type == DisplayProperty::Type_Volume) {
            if (ImGui::BeginCombo("File Format", volume_formats[volume_format].lbl.ptr)) {
                for (int i = 0; i < (int)ARRAY_SIZE(volume_formats); ++i) {
                    if (ImGui::Selectable(volume_formats[i].lbl.ptr, volume_format == i)) {
                        volume_format = i;
                    }
                }
                ImGui::EndCombo();
            }
            file_extension = volume_formats[volume_format].ext;
        }

        if (property_idx != -1) {
            if (type != data->display_properties[property_idx].type) {
                property_idx = -1;
            }
        }

        if (property_idx == -1) ImGui::PushDisabled();
        bool export_clicked = ImGui::Button("Export");
        if (property_idx == -1) ImGui::PopDisabled();

        if (export_clicked) {
            md_allocator_i* alloc = frame_alloc;
            md_vm_arena_temp_t temp = md_vm_arena_temp_begin(frame_alloc);
            defer { md_vm_arena_temp_end(temp); };

            ASSERT(property_idx != -1);
            char path_buf[1024];
            DisplayProperty& dp = data->display_properties[property_idx];
            md_array(const float*)  column_data = 0;
            md_array(str_t)         column_labels = 0;
            md_array(str_t)         legends = 0;

            if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, file_extension)) {
                str_t path = {path_buf, strnlen(path_buf, sizeof(path_buf))};
                if (dp.type == DisplayProperty::Type_Volume) {
                    if (str_eq(file_extension, STR_LIT("cube"))) {
                        if (export_cube(*data, dp.prop_data, dp.vis_payload, path)) {
                            LOG_SUCCESS("Successfully exported property '%s' to '" STR_FMT "'", dp.label, STR_ARG(path));
                        }
                    }
                } else {
                    md_file_o* file = md_file_open(path, MD_FILE_WRITE | MD_FILE_BINARY);
                    if (file) {
                        str_t out_str = {};
                        if (dp.type == DisplayProperty::Type_Temporal) {
                            const double* traj_times = md_trajectory_frame_times(data->mold.traj);
                            const size_t  num_frames = md_trajectory_num_frames(data->mold.traj);
                            md_array(float) time = md_array_create(float, num_frames, alloc);
                            for (size_t i = 0; i < num_frames; ++i) {
                                time[i] = (float)traj_times[i];
                            }

                            str_t x_label = STR_LIT("Frame");
                            str_t y_label = str_from_cstr(dp.label);

                            if (!md_unit_unitless(dp.unit[1])) {
                                y_label = str_printf(alloc, "%s (%s)", dp.label, dp.unit_str);
                            }

                            md_unit_t time_unit = md_trajectory_time_unit(data->mold.traj);
                            if (!md_unit_empty(time_unit)) {
                                char time_buf[64];
                                size_t len = md_unit_print(time_buf, sizeof(time_buf), time_unit);
                                x_label = str_printf(alloc, "Time (" STR_FMT ")", len, time_buf);
                            }

                            md_array_push(column_data, time, alloc);
                            md_array_push(column_labels, x_label, alloc);

                            if (dp.dim > 1) {
                                for (int i = 0; i < dp.dim; ++i) {
                                    str_t  legend = str_printf(alloc, "%s[%i]", dp.label, i + 1);
                                    float* values = (float*)md_alloc(alloc, sizeof(float) * num_frames);
                                    for (size_t j = 0; j < num_frames; ++j) {
                                        values[j] = dp.y_values[j * dp.dim + i];
                                    }
                                    md_array_push(column_data, values, alloc);
                                    md_array_push(legends, legend, alloc);
                                    md_array_push(column_labels, legend, alloc);
                                }
                            } else {
                                md_array_push(column_data, dp.y_values, alloc);
                                md_array_push(column_labels, y_label, alloc);
                            }

                            if (str_eq(file_extension, STR_LIT("xvg"))) {
                                str_t header = md_xvg_format_header(str_from_cstr(dp.label), x_label, y_label, md_array_size(legends), legends, alloc);
                                out_str = md_xvg_format(header, md_array_size(column_data), num_frames, column_data, alloc);
                            } else if (str_eq(file_extension, STR_LIT("csv"))) {
                                out_str = md_csv_write_to_str(column_data, column_labels, md_array_size(column_data), num_frames, alloc);
                            }


                        } else if (dp.type == DisplayProperty::Type_Distribution) {
                            md_array(float) x_values = sample_range(dp.hist.x_min, dp.hist.x_max, dp.hist.num_bins, alloc);

                            str_t x_label = str_from_cstr(dp.unit_str[0]);
                            str_t y_label = str_from_cstr(dp.label);
                            if (strlen(dp.unit_str[1]) > 0) {
                                y_label = str_printf(frame_alloc, "%s (%s)", dp.label, dp.unit_str);
                            }

                            md_array_push(column_data, x_values, alloc);
                            md_array_push(column_labels, x_label, alloc);

                            if (dp.hist.dim > 1) {
                                for (int i = 0; i < dp.hist.dim; ++i) {
                                    md_array_push(column_data, dp.hist.bins + i * dp.hist.num_bins, alloc);
                                    str_t legend = str_printf(alloc, "%s[%i]", dp.label, i + 1);
                                    md_array_push(legends, legend, alloc);
                                    md_array_push(column_labels, legend, alloc);
                                }
                            } else {
                                md_array_push(column_data, dp.hist.bins, alloc);
                                md_array_push(column_labels, y_label, alloc);
                            }

                            if (str_eq(file_extension, STR_LIT("xvg"))) {
                                str_t header = md_xvg_format_header(str_from_cstr(dp.label), x_label, y_label, md_array_size(legends), legends, alloc);
                                out_str = md_xvg_format(header, md_array_size(column_data), dp.hist.num_bins, column_data, alloc);
                            } else if (str_eq(file_extension, STR_LIT("csv"))) {
                                out_str = md_csv_write_to_str(column_data, column_labels, md_array_size(column_data), dp.hist.num_bins, alloc);
                            }
                        }
                        if (!str_empty(out_str)) {
                            md_file_write(file, out_str.ptr, out_str.len);
                            LOG_SUCCESS("Successfully exported property '%s' to '%.*s'", dp.label, (int)path.len, path.ptr);
                        }
                        md_file_close(file);
                    }
                }
            }
        }
        ImGui::PopItemWidth();
    }
    ImGui::End();
}

static void update_md_buffers(ApplicationState* data) {
    ASSERT(data);
    const auto& mol = data->mold.mol;

    if (mol.atom.count == 0) return;

    if (data->mold.dirty_buffers & MolBit_DirtyPosition) {
        const vec3_t pbc_ext = data->mold.mol.unit_cell.basis * vec3_t{1,1,1};
        md_gl_mol_set_atom_position(data->mold.gl_mol, 0, (uint32_t)mol.atom.count, mol.atom.x, mol.atom.y, mol.atom.z, 0);
        if (!(data->mold.dirty_buffers & MolBit_ClearVelocity)) {
            md_gl_mol_compute_velocity(data->mold.gl_mol, pbc_ext.elem);
        }
#if EXPERIMENTAL_GFX_API
        md_gfx_structure_set_atom_position(data->mold.gfx_structure, 0, (uint32_t)mol.atom.count, mol.atom.x, mol.atom.y, mol.atom.z, 0);
        md_gfx_structure_set_aabb(data->mold.gfx_structure, &data->mold.mol_aabb_min, &data->mold.mol_aabb_max);
#endif
    }

    if (data->mold.dirty_buffers & MolBit_ClearVelocity) {
        md_gl_mol_zero_velocity(data->mold.gl_mol);
    }

    if (data->mold.dirty_buffers & MolBit_DirtyRadius) {
        md_gl_mol_set_atom_radius(data->mold.gl_mol, 0, (uint32_t)mol.atom.count, mol.atom.radius, 0);
#if EXPERIMENTAL_GFX_API
        md_gfx_structure_set_atom_radius(data->mold.gfx_structure, 0, (uint32_t)mol.atom.count, mol.atom.radius, 0);
#endif
    }

    if (data->mold.dirty_buffers & MolBit_DirtyFlags) {
        md_vm_arena_temp_t tmp = md_vm_arena_temp_begin(frame_alloc);
        defer { md_vm_arena_temp_end(tmp); };

        uint8_t* flags = (uint8_t*)md_vm_arena_push(frame_alloc, mol.atom.count * sizeof(uint8_t));
        MEMSET(flags, 0, mol.atom.count * sizeof(uint8_t));

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
        md_gl_mol_set_atom_flags(data->mold.gl_mol, 0, (uint32_t)mol.atom.count, flags, 0);
    }

    if (data->mold.dirty_buffers & MolBit_DirtyBonds) {
        md_gl_mol_set_bonds(data->mold.gl_mol, 0, (uint32_t)mol.bond.count, mol.bond.pairs, sizeof(md_bond_pair_t));
    }

    if (data->mold.dirty_buffers & MolBit_DirtySecondaryStructure) {
        if (mol.protein_backbone.secondary_structure) {
            md_gl_mol_set_backbone_secondary_structure(data->mold.gl_mol, 0, (uint32_t)mol.protein_backbone.count, mol.protein_backbone.secondary_structure, 0);
        }
    }

    data->mold.dirty_buffers = 0;
}

static void interrupt_async_tasks(ApplicationState* data) {
    task_system::pool_interrupt_running_tasks();

    if (data->script.full_eval) md_script_eval_interrupt(data->script.full_eval);
    if (data->script.filt_eval) md_script_eval_interrupt(data->script.filt_eval);

    task_system::pool_wait_for_completion();
}

// #trajectorydata
static void free_trajectory_data(ApplicationState* data) {
    ASSERT(data);
    interrupt_async_tasks(data);

    if (data->mold.traj) {
        load::traj::close(data->mold.traj);
        data->mold.traj = nullptr;
    }
    data->files.trajectory[0] = '\0';
    
    data->mold.mol.unit_cell = {};
    md_array_free(data->timeline.x_values,  persistent_alloc);
    md_array_free(data->display_properties, persistent_alloc);

    md_array_free(data->trajectory_data.backbone_angles.data,     persistent_alloc);
    md_array_free(data->trajectory_data.secondary_structure.data, persistent_alloc);
}

static void init_trajectory_data(ApplicationState* data) {
    size_t num_frames = md_trajectory_num_frames(data->mold.traj);
    if (num_frames > 0) {
        size_t min_frame = 0;
        size_t max_frame = num_frames - 1;
        md_trajectory_header_t header;
        md_trajectory_get_header(data->mold.traj, &header);
        
        double min_time = header.frame_times[0];
        double max_time = header.frame_times[num_frames - 1];

        data->timeline.view_range = {min_time, max_time};
        data->timeline.filter.beg_frame = min_frame;
        data->timeline.filter.end_frame = max_frame;

        md_array_resize(data->timeline.x_values, num_frames, persistent_alloc);
        for (size_t i = 0; i < num_frames; ++i) {
            data->timeline.x_values[i] = header.frame_times[i];
        }

        data->animation.frame = CLAMP(data->animation.frame, (double)min_frame, (double)max_frame);
        int64_t frame_idx = CLAMP((int64_t)(data->animation.frame + 0.5), 0, (int64_t)max_frame);

        md_trajectory_frame_header_t frame_header;
        md_trajectory_load_frame(data->mold.traj, frame_idx, &frame_header, data->mold.mol.atom.x, data->mold.mol.atom.y, data->mold.mol.atom.z);
        data->mold.mol.unit_cell = frame_header.unit_cell;

        if (data->mold.mol.protein_backbone.count > 0) {
            data->trajectory_data.secondary_structure.stride = data->mold.mol.protein_backbone.count;
            data->trajectory_data.secondary_structure.count = data->mold.mol.protein_backbone.count * num_frames;
            md_array_resize(data->trajectory_data.secondary_structure.data, data->mold.mol.protein_backbone.count * num_frames, persistent_alloc);
            for (size_t i = 0; i < md_array_size(data->trajectory_data.secondary_structure.data); ++i) {
                data->trajectory_data.secondary_structure.data[i] = MD_SECONDARY_STRUCTURE_COIL;
            }
//            MEMSET(data->trajectory_data.secondary_structure.data, 0, md_array_size(data->trajectory_data.secondary_structure.data) * sizeof (md_secondary_structure_t));

            data->trajectory_data.backbone_angles.stride = data->mold.mol.protein_backbone.count;
            data->trajectory_data.backbone_angles.count = data->mold.mol.protein_backbone.count * num_frames;
            md_array_resize(data->trajectory_data.backbone_angles.data, data->mold.mol.protein_backbone.count * num_frames, persistent_alloc);
            MEMSET(data->trajectory_data.backbone_angles.data, 0, md_array_size(data->trajectory_data.backbone_angles.data) * sizeof (md_backbone_angles_t));

            // Launch work to compute the values
            task_system::task_interrupt_and_wait_for(data->tasks.backbone_computations);

            data->tasks.backbone_computations = task_system::create_pool_task(STR_LIT("Backbone Operations"), (uint32_t)num_frames, [data](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                // Create copy here of molecule since we use the full structure as input
                md_molecule_t mol = data->mold.mol;

                const size_t stride = ALIGN_TO(mol.atom.count, 8);
                const size_t bytes = stride * sizeof(float) * 3;
                float* coords = (float*)md_alloc(md_get_heap_allocator(), bytes);
                defer { md_free(md_get_heap_allocator(), coords, bytes); };
                // Overwrite the coordinate section, since we will load trajectory frame data into these
                mol.atom.x = coords + stride * 0;
                mol.atom.y = coords + stride * 1;
                mol.atom.z = coords + stride * 2;

                for (uint32_t frame_idx = range_beg; frame_idx < range_end; ++frame_idx) {
                    md_trajectory_load_frame(data->mold.traj, frame_idx, NULL, mol.atom.x, mol.atom.y, mol.atom.z);
                    md_util_backbone_angles_compute(data->trajectory_data.backbone_angles.data + data->trajectory_data.backbone_angles.stride * frame_idx, data->trajectory_data.backbone_angles.stride, &mol);
                    md_util_backbone_secondary_structure_compute(data->trajectory_data.secondary_structure.data + data->trajectory_data.secondary_structure.stride * frame_idx, data->trajectory_data.secondary_structure.stride, &mol);
                }
            });

            task_system::ID main_task = task_system::create_main_task(STR_LIT("Update Trajectory Data"), [data]() {
                data->trajectory_data.backbone_angles.fingerprint = generate_fingerprint();
                data->trajectory_data.secondary_structure.fingerprint = generate_fingerprint();
                
                interpolate_atomic_properties(data);
                data->mold.dirty_buffers |= MolBit_ClearVelocity;
                update_all_representations(data);

            });

            task_system::set_task_dependency(main_task, data->tasks.backbone_computations);
            task_system::enqueue_task(data->tasks.backbone_computations);
        }

        data->mold.dirty_buffers |= MolBit_DirtyPosition;
        data->mold.dirty_buffers |= MolBit_ClearVelocity;

        // Prefetch frames
        //launch_prefetch_job(data);
    }
}

static bool load_trajectory_data(ApplicationState* data, str_t filename, md_trajectory_loader_i* loader, LoadTrajectoryFlags flags) {
    md_trajectory_i* traj = load::traj::open_file(filename, loader, &data->mold.mol, persistent_alloc, flags);
    if (traj) {
        free_trajectory_data(data);
        data->mold.traj = traj;
        str_copy_to_char_buf(data->files.trajectory, sizeof(data->files.trajectory), filename);
        init_trajectory_data(data);
        data->animation.frame = 0;
        return true;
    }

    return false;
}

// #moleculedata
static void free_molecule_data(ApplicationState* data) {
    ASSERT(data);
    interrupt_async_tasks(data);

    clear_dataset_items(data);

    //md_molecule_free(&data->mold.mol, persistent_alloc);
    md_arena_allocator_reset(data->mold.mol_alloc);
    MEMSET(&data->mold.mol, 0, sizeof(data->mold.mol));

    md_gl_mol_destroy(data->mold.gl_mol);
    MEMSET(data->files.molecule, 0, sizeof(data->files.molecule));

    md_bitfield_clear(&data->selection.selection_mask);
    md_bitfield_clear(&data->selection.highlight_mask);
    if (data->script.ir) {
        md_script_ir_free(data->script.ir);
        data->script.ir = nullptr;
    }
    if (data->script.eval_ir && data->script.eval_ir != data->script.ir) {
        md_script_ir_free(data->script.eval_ir);
        data->script.eval_ir = nullptr;
    }
    if (data->script.full_eval) {
        md_script_eval_free(data->script.full_eval);
        data->script.full_eval = nullptr;
    }
    if (data->script.filt_eval) {
        md_script_eval_free(data->script.filt_eval);
        data->script.filt_eval = nullptr;
    }
    clear_density_volume(data);

    viamd::event_system_broadcast_event(viamd::EventType_ViamdTopologyFree, viamd::EventPayloadType_ApplicationState, data);
}

static void init_molecule_data(ApplicationState* data) {
    if (data->mold.mol.atom.count) {

        data->picking.idx = INVALID_PICKING_IDX;
        data->selection.atom_idx.hovered = -1;
        data->selection.atom_idx.right_click = -1;
        data->selection.bond_idx.hovered = -1;
        data->selection.bond_idx.right_click = -1;

        data->mold.gl_mol = md_gl_mol_create(&data->mold.mol);

#if EXPERIMENTAL_GFX_API
        const md_molecule_t& mol = data->mold.mol;
        vec3_t& aabb_min = data->mold.mol_aabb_min;
        vec3_t& aabb_max = data->mold.mol_aabb_max;
        md_util_compute_aabb_soa(&aabb_min, &aabb_max, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.radius, mol.atom.count);

        data->mold.gfx_structure = md_gfx_structure_create(mol.atom.count, mol.covalent.count, mol.backbone.count, mol.backbone.range_count, mol.residue.count, mol.instance.count);
        md_gfx_structure_set_atom_position(data->mold.gfx_structure, 0, mol.atom.count, mol.atom.x, mol.atom.y, mol.atom.z, 0);
        md_gfx_structure_set_atom_radius(data->mold.gfx_structure, 0, mol.atom.count, mol.atom.radius, 0);
        md_gfx_structure_set_aabb(data->mold.gfx_structure, &data->mold.mol_aabb_min, &data->mold.mol_aabb_max);
        if (mol.instance.count > 0) {
            md_gfx_structure_set_instance_atom_ranges(data->mold.gfx_structure, 0, mol.instance.count, (md_gfx_range_t*)mol.instance.atom_range, 0);
            md_gfx_structure_set_instance_transforms(data->mold.gfx_structure, 0, mol.instance.count, mol.instance.transform, 0);
        }
#endif
    }
    viamd::event_system_broadcast_event(viamd::EventType_ViamdTopologyInit, viamd::EventPayloadType_ApplicationState, data);

    update_representation_info(data);
    init_all_representations(data);
    update_all_representations(data);
    data->script.compile_ir = true;

    init_dataset_items(data);
}

static void launch_prefetch_job(ApplicationState* data) {
    const uint32_t num_frames = MIN((uint32_t)md_trajectory_num_frames(data->mold.traj), (uint32_t)load::traj::num_cache_frames(data->mold.traj));
    if (!num_frames) return;

    task_system::task_interrupt_and_wait_for(data->tasks.prefetch_frames);
    data->tasks.prefetch_frames = task_system::create_pool_task(STR_LIT("Prefetch Frames"), num_frames, [data](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
        (void)thread_num;
        for (uint32_t i = range_beg; i < range_end; ++i) {
            md_trajectory_frame_header_t header;
            md_trajectory_load_frame(data->mold.traj, i, &header, 0, 0, 0);
        }
    });

    task_system::ID main_task = task_system::create_main_task(STR_LIT("Prefetch Complete"), [data]() {
        interpolate_atomic_properties(data);
        data->mold.dirty_buffers |= MolBit_ClearVelocity;
    });

    task_system::set_task_dependency(main_task, data->tasks.prefetch_frames);
    task_system::enqueue_task(data->tasks.prefetch_frames);
}

static bool load_dataset_from_file(ApplicationState* data, const LoadParam& param) {
    ASSERT(data);

    str_t path_to_file = md_path_make_canonical(param.file_path, frame_alloc);
    if (path_to_file) {
        if (param.mol_loader) {
            interrupt_async_tasks(data);
            free_trajectory_data(data);
            free_molecule_data(data);

            if (!param.mol_loader->init_from_file(&data->mold.mol, path_to_file, param.mol_loader_arg, data->mold.mol_alloc)) {
                LOG_ERROR("Failed to load molecular data from file '" STR_FMT "'", STR_ARG(path_to_file));
                return false;
            }
            LOG_SUCCESS("Successfully loaded molecular data from file '" STR_FMT "'", STR_ARG(path_to_file));

            str_copy_to_char_buf(data->files.molecule, sizeof(data->files.molecule), path_to_file);
            data->files.coarse_grained = param.coarse_grained;
            // @NOTE: If the dataset is coarse-grained, then postprocessing must be aware
            md_util_postprocess_flags_t flags = param.coarse_grained ? MD_UTIL_POSTPROCESS_COARSE_GRAINED : MD_UTIL_POSTPROCESS_ALL;
            md_util_molecule_postprocess(&data->mold.mol, data->mold.mol_alloc, flags);
            init_molecule_data(data);

            // @NOTE: Some files contain both atomic coordinates and trajectory
            if (param.traj_loader) {
                LOG_INFO("File may also contain trajectory, attempting to load trajectory");
            } else {
                return true;
            }
        }

        if (param.traj_loader) {
            if (!data->mold.mol.atom.count) {
                LOG_ERROR("Before loading a trajectory, molecular data needs to be present");
                return false;
            }
            interrupt_async_tasks(data);

            bool success = load_trajectory_data(data, path_to_file, param.traj_loader, param.traj_loader_flags);
            if (success) {
                LOG_SUCCESS("Successfully opened trajectory from file '" STR_FMT "'", STR_ARG(path_to_file));
                return true;
            } else {
                if (param.mol_loader && param.traj_loader) {
                    // Don't record this as an error, as the trajectory may be optional (In case of PDB for example)
                    return true;
                }
                LOG_ERROR("Failed to opened trajectory from file '" STR_FMT "'", STR_ARG(path_to_file));
            }
        }
    }

    return false;
}

static void load_workspace(ApplicationState* data, str_t filename) {
    str_t txt = load_textfile(filename, frame_alloc);
    defer { str_free(txt, frame_alloc); };

    if (str_empty(txt)) {
        LOG_ERROR("Could not open workspace file: '" STR_FMT "'", STR_ARG(filename));
        return;
    }

    // Reset and clear things
    clear_selections(data);
    clear_representations(data);
    editor.SetText("");
    data->files.workspace[0]  = '\0';

    data->animation = {};
    reset_view(data, &data->representation.visibility_mask, false, true);

    str_t new_molecule_file   = {};
    str_t new_trajectory_file = {};
    bool  new_coarse_grained  = false;
    double new_frame = 0;

    str_t folder = {};
    extract_folder_path(&folder, filename);

    viamd::deserialization_state_t state = {
        .filename = filename,
        .text = txt,
    };
    str_t section;
    while (viamd::next_section_header(section, state)) {
        if (str_eq(section, STR_LIT("Files")) || str_eq(section, STR_LIT("File"))) {
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("MoleculeFile"))) {
                    str_t file;
                    viamd::extract_str(file, arg);
                    if (!str_empty(file)) {
                        md_strb_t path = md_strb_create(frame_alloc);
                        path += folder;
                        path += file;
                        new_molecule_file = md_path_make_canonical(path, frame_alloc);
                    }
                } else if (str_eq(ident, STR_LIT("TrajectoryFile"))) {
                    str_t file;
                    viamd::extract_str(file, arg);
                    if (!str_empty(file)) {
                        md_strb_t path = md_strb_create(frame_alloc);
                        path += folder;
                        path += file;
                        new_trajectory_file = md_path_make_canonical(path, frame_alloc);
                    }
                } else if (str_eq(ident, STR_LIT("CoarseGrained"))) {
                    viamd::extract_bool(new_coarse_grained, arg);
                }
            }
        } else if (str_eq(section, STR_LIT("Animation"))) {
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Frame"))) {
                    viamd::extract_dbl(new_frame, arg);
                } else if (str_eq(ident, STR_LIT("Fps"))) {
                    viamd::extract_flt(data->animation.fps, arg);
                } else if (str_eq(ident, STR_LIT("Interpolation"))) {
                    int mode;
                    viamd::extract_int(mode, arg);
                    data->animation.interpolation = (InterpolationMode)mode;
                }
            }
        } else if (str_eq(section, STR_LIT("RenderSettings"))) {
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("SsaoEnabled"))) {
                    viamd::extract_bool(data->visuals.ssao.enabled, arg);
                } else if (str_eq(ident, STR_LIT("SsaoIntensity"))) {
                    viamd::extract_flt(data->visuals.ssao.intensity, arg);
                } else if (str_eq(ident, STR_LIT("SsaoRadius"))) {
                    viamd::extract_flt(data->visuals.ssao.radius, arg);
                } else if (str_eq(ident, STR_LIT("SsaoBias"))) {
                    viamd::extract_flt(data->visuals.ssao.bias, arg);
                } else if (str_eq(ident, STR_LIT("DofEnabled"))) {
                    viamd::extract_bool(data->visuals.dof.enabled, arg);
                } else if (str_eq(ident, STR_LIT("DofFocusScale"))) {
                    viamd::extract_flt(data->visuals.dof.focus_scale, arg);
                }
            }
        } else if (str_eq(section, STR_LIT("Camera"))) {
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Position"))) {
                    viamd::extract_vec3(data->view.camera.position, arg);
                } else if (str_eq(ident, STR_LIT("Orientation"))) {
                    viamd::extract_quat(data->view.camera.orientation, arg);
                } else if (str_eq(ident, STR_LIT("Distance"))) {
                    viamd::extract_flt(data->view.camera.focus_distance, arg);
                } else if (str_eq(ident, STR_LIT("Rotation"))) {
                    // DEPRECATED
                    viamd::extract_quat(data->view.camera.orientation, arg);
                }
            }
        } else if (str_eq(section, STR_LIT("Representation"))) {
            Representation* rep = create_representation(data);
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Name"))) {
                    viamd::extract_to_char_buf(rep->name, sizeof(rep->name), arg);
                } else if (str_eq(ident, STR_LIT("Filter"))) {
                    viamd::extract_to_char_buf(rep->filt, sizeof(rep->filt), arg);
                } else if (str_eq(ident, STR_LIT("Enabled"))) {
                    viamd::extract_bool(rep->enabled, arg);
                } else if (str_eq(ident, STR_LIT("Type"))) {
                    int type;
                    viamd::extract_int(type, arg);
                    type = CLAMP(type, 0, (int)RepresentationType::Count);
                    rep->type = (RepresentationType)type;
                } else if (str_eq(ident, STR_LIT("ColorMapping"))) {
                    int mapping;
                    viamd::extract_int(mapping, arg);
                    mapping = CLAMP(mapping, 0, (int)ColorMapping::Count);
                    rep->color_mapping = (ColorMapping)mapping;
                } else if (str_eq(ident, STR_LIT("StaticColor"))) {
                    viamd::extract_flt_vec(rep->uniform_color.elem, 4, arg);
                } else if (str_eq(ident, STR_LIT("Saturation"))) {
                    viamd::extract_flt(rep->saturation, arg);
                } else if (str_eq(ident, STR_LIT("Radius"))) {
                    // DEPRECATED
                    viamd::extract_flt(rep->scale.x, arg);
                } else if (str_eq(ident, STR_LIT("Tension"))) {
                    // DEPRECATED
                } else if (str_eq(ident, STR_LIT("Width"))) {
                    viamd::extract_flt(rep->scale.x, arg);
                } else if (str_eq(ident, STR_LIT("Thickness"))) {
                    viamd::extract_flt(rep->scale.y, arg);
                } else if (str_eq(ident, STR_LIT("Param"))) {
                    viamd::extract_flt_vec(rep->scale.elem, 4, arg);
                } else if (str_eq(ident, STR_LIT("DynamicEval"))) {
                    viamd::extract_bool(rep->dynamic_evaluation, arg);
                } else if (str_eq(ident, STR_LIT("OrbIdx"))) {
                    viamd::extract_int(rep->orbital.orbital_idx, arg);
                } else if (str_eq(ident, STR_LIT("OrbRes"))) {
                    int res;
                    viamd::extract_int(res, arg);
                    rep->orbital.resolution = (VolumeResolution)res;
                } else if (str_eq(ident, STR_LIT("OrbType"))) {
                    int type;
                    viamd::extract_int(type, arg);
                    rep->orbital.type = (OrbitalType)type;
                } else if (str_eq(ident, STR_LIT("OrbIso"))) {
                    viamd::extract_flt(rep->orbital.iso_psi.values[0], arg);
                    rep->orbital.iso_psi.values[1] = - rep->orbital.iso_psi.values[0];
                    rep->orbital.iso_den.values[0] = rep->orbital.iso_psi.values[0] * rep->orbital.iso_psi.values[0];
                } else if (str_eq(ident, STR_LIT("OrbColPos"))) {
                    viamd::extract_vec4(rep->orbital.iso_psi.colors[0], arg);
                } else if (str_eq(ident, STR_LIT("OrbColNeg"))) {
                    viamd::extract_vec4(rep->orbital.iso_psi.colors[1], arg);
                } else if (str_eq(ident, STR_LIT("OrbColDen"))) {
                    viamd::extract_vec4(rep->orbital.iso_den.colors[0], arg);
                }
            }
        } else if (str_eq(section, STR_LIT("AtomElementMapping"))) {
            str_t lbl = {};
            int elem = 0;
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Label"))) {
                    viamd::extract_str(lbl, arg);
                } else if (str_eq(ident, STR_LIT("Element"))) {
                    viamd::extract_int(elem, arg);
                }
            }
            if (!str_empty(lbl) && elem) {
                add_atom_elem_mapping(data, lbl, (md_element_t)elem);
            }
        } else if (str_eq(section, STR_LIT("Script"))) {
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Text"))) {
                    str_t str;
                    viamd::extract_str(str, arg);
                    editor.SetText(std::string(str.ptr, str.len));
                }
            }
        } else if (str_eq(section, STR_LIT("Selection"))) {
            str_t ident, arg;
            str_t label = {};
            str_t mask_base64 = {};
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Label"))) {
                    viamd::extract_str(label, arg);
                } else if (str_eq(ident, STR_LIT("Mask"))) {
                    viamd::extract_str(mask_base64, arg);
                }
            }
            if (!str_empty(label) && !str_empty(mask_base64)) {
                Selection* sel = create_selection(data, label);
                md_bitfield_deserialize(&sel->atom_mask, mask_base64.ptr, mask_base64.len);
            }
        } else {
            viamd::event_system_broadcast_event(viamd::EventType_ViamdDeserialize, viamd::EventPayloadType_DeserializationState, &state);
        }
    }

    data->view.animation.target_position    = data->view.camera.position;
    data->view.animation.target_orientation = data->view.camera.orientation;
    data->view.animation.target_distance    = data->view.camera.focus_distance;
    
    str_copy_to_char_buf(data->files.workspace, sizeof(data->files.workspace), filename);
    
    data->files.coarse_grained  = new_coarse_grained;

    load::LoaderState loader_state = {};
    load::init_loader_state(&loader_state, new_molecule_file, frame_alloc);

    LoadParam param = {};
    param.file_path = new_molecule_file;
    param.mol_loader = loader_state.mol_loader;
    param.traj_loader = loader_state.traj_loader;
    param.coarse_grained = new_coarse_grained;
    param.mol_loader_arg = loader_state.mol_loader_arg;

    if (new_molecule_file && load_dataset_from_file(data, param)) {
        str_copy_to_char_buf(data->files.molecule, sizeof(data->files.molecule), new_molecule_file);
    } else {
        data->files.molecule[0]   = '\0';
    }

    if (new_trajectory_file) {
        load::init_loader_state(&loader_state, new_trajectory_file, frame_alloc);
        param.mol_loader = 0;
        param.traj_loader = loader_state.traj_loader;
        param.file_path = new_trajectory_file;
        if (load_dataset_from_file(data, param)) {
            str_copy_to_char_buf(data->files.trajectory, sizeof(data->files.trajectory), new_trajectory_file);
        }
        data->animation.frame = new_frame;
    } else {
        data->files.trajectory[0] = '\0';
    }

    apply_atom_elem_mappings(data);
}

static void save_workspace(ApplicationState* data, str_t filename) {
    md_file_o* file = md_file_open(filename, MD_FILE_WRITE);
    if (!file) {
        LOG_ERROR("Could not open workspace file for writing: '%.*s", (int)filename.len, filename.ptr);
        return;
    }
    defer { md_file_close(file); };

    viamd::serialization_state_t state {
        .filename = filename,
        .sb = md_strb_create(frame_alloc),
    };

    // Write big ass header
    state.sb += header_snippet;
    state.sb += '\n';

    str_t mol_file  = md_path_make_relative(filename, str_from_cstr(data->files.molecule),   frame_alloc);
    str_t traj_file = md_path_make_relative(filename, str_from_cstr(data->files.trajectory), frame_alloc);

    viamd::write_section_header(state, STR_LIT("Files"));
    viamd::write_str(state, STR_LIT("MoleculeFile"), mol_file);
    viamd::write_str(state, STR_LIT("TrajectoryFile"), traj_file);
    viamd::write_int(state, STR_LIT("CoarseGrained"), data->files.coarse_grained);

    viamd::write_section_header(state, STR_LIT("Animation"));
    viamd::write_dbl(state, STR_LIT("Frame"), data->animation.frame);
    viamd::write_flt(state, STR_LIT("Fps"), data->animation.fps);
    viamd::write_int(state, STR_LIT("Interpolation"), (int)data->animation.interpolation);

    viamd::write_section_header(state, STR_LIT("RenderSettings"));
    viamd::write_bool(state, STR_LIT("SsaoEnabled"), data->visuals.ssao.enabled);
    viamd::write_flt(state, STR_LIT("SsaoIntensity"), data->visuals.ssao.intensity);
    viamd::write_flt(state, STR_LIT("SsaoRadius"), data->visuals.ssao.radius);
    viamd::write_bool(state, STR_LIT("DofEnabled"), data->visuals.dof.enabled);
    viamd::write_flt(state, STR_LIT("DofFocusScale"), data->visuals.dof.focus_scale);

    viamd::write_section_header(state, STR_LIT("Camera"));
    viamd::write_vec3(state, STR_LIT("Position"), data->view.camera.position);
    viamd::write_quat(state, STR_LIT("Orientation"), data->view.camera.orientation);
    viamd::write_flt(state,  STR_LIT("Distance"), data->view.camera.focus_distance);
    viamd::write_int(state,  STR_LIT("Mode"), (int)data->view.mode);


    for (size_t i = 0; i < md_array_size(data->representation.reps); ++i) {
        const Representation& rep = data->representation.reps[i];
        viamd::write_section_header(state, STR_LIT("Representation"));
        viamd::write_str(state,  STR_LIT("Name"), str_from_cstr(rep.name));
        viamd::write_str(state,  STR_LIT("Filter"), str_from_cstr(rep.filt));
        viamd::write_bool(state, STR_LIT("Enabled"), rep.enabled);
        viamd::write_int(state,  STR_LIT("Type"), (int)rep.type);
        viamd::write_int(state,  STR_LIT("ColorMapping"), (int)rep.color_mapping);
        viamd::write_vec4(state, STR_LIT("StaticColor"), rep.uniform_color);
        viamd::write_flt(state,  STR_LIT("Saturation"), rep.saturation);
        viamd::write_vec4(state, STR_LIT("Param"), rep.scale);
        viamd::write_bool(state, STR_LIT("DynamicEval"), rep.dynamic_evaluation);

        if (rep.type == RepresentationType::Orbital) {
            viamd::write_int(state,  STR_LIT("OrbIdx"),      rep.orbital.orbital_idx);
            viamd::write_int(state,  STR_LIT("OrbType"),(int)rep.orbital.type);
            viamd::write_int(state,  STR_LIT("OrbRes"), (int)rep.orbital.resolution);
            viamd::write_flt(state,  STR_LIT("OrbIso"),      rep.orbital.iso_psi.values[0]);
            viamd::write_vec4(state, STR_LIT("OrbColPos"),   rep.orbital.iso_psi.colors[0]);
            viamd::write_vec4(state, STR_LIT("OrbColNeg"),   rep.orbital.iso_psi.colors[1]);
            viamd::write_vec4(state, STR_LIT("OrbColDen"),   rep.orbital.iso_den.colors[0]);
        }
    }

    for (size_t i = 0; i < md_array_size(data->dataset.atom_element_remappings); ++i) {
        const AtomElementMapping& mapping = data->dataset.atom_element_remappings[i];
        viamd::write_section_header(state, STR_LIT("AtomElementMapping"));
        viamd::write_str(state, STR_LIT("Label"), str_from_cstr(mapping.lbl));
        viamd::write_int(state, STR_LIT("Element"), mapping.elem);
    }

    {
        std::string text = editor.GetText();
        viamd::write_section_header(state, STR_LIT("Script"));
        viamd::write_str(state, STR_LIT("Text"), str_t{text.c_str(), text.size()});
    }

    for (size_t i = 0; i < md_array_size(data->selection.stored_selections); ++i) {
        const Selection& sel = data->selection.stored_selections[i];
        size_t cap = md_bitfield_serialize_size_in_bytes(&sel.atom_mask);
        char* buf = (char*)md_vm_arena_push(frame_alloc, cap);
        size_t len = md_bitfield_serialize(buf, &sel.atom_mask);
        str_t encoded_mask = {buf, len};
        viamd::write_section_header(state, STR_LIT("Selection"));
        viamd::write_str(state, STR_LIT("Label"), str_from_cstr(sel.name));
        viamd::write_str(state, STR_LIT("Mask"),  encoded_mask);
    }

    viamd::event_system_broadcast_event(viamd::EventType_ViamdSerialize, viamd::EventPayloadType_SerializationState, &state);

    str_t text = md_strb_to_str(state.sb);
    md_file_write(file, str_ptr(text), str_len(text));
}

void create_screenshot(ApplicationState* data) {
    ASSERT(data);

    md_vm_arena_temp_t tmp = md_vm_arena_temp_begin(frame_alloc);
    defer { md_vm_arena_temp_end(tmp); };

    str_t path = data->screenshot.path_to_file;

    int width  = data->gbuffer.width;
    int height = data->gbuffer.height;
    size_t bytes = width * height * sizeof(uint32_t);
    uint32_t* rgba = (uint32_t*)md_vm_arena_push(frame_alloc, bytes);
    defer { md_vm_arena_pop(frame_alloc, bytes); };

    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    glReadBuffer(GL_BACK);
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
        LOG_ERROR("Non supported file-extension '%.*s' when saving screenshot", (int)ext.len, ext.ptr);
        return;
    }

    LOG_SUCCESS("Screenshot saved to: '%.*s'", (int)path.len, path.ptr);
}

// #representation
static Representation* create_representation(ApplicationState* data, RepresentationType type, ColorMapping color_mapping, str_t filter) {
    ASSERT(data);
    md_array_push(data->representation.reps, Representation(), persistent_alloc);
    Representation* rep = md_array_last(data->representation.reps);
    rep->type = type;
    rep->color_mapping = color_mapping;
    if (!str_empty(filter)) {
        str_copy_to_char_buf(rep->filt, sizeof(rep->filt), filter);
    }
    rep->orbital.orbital_idx = data->representation.info.mo_homo_idx;
    init_representation(data, rep);
    update_representation(data, rep);
    return rep;
}

static Representation* clone_representation(ApplicationState* state, const Representation& rep) {
    ASSERT(state);
    md_array_push(state->representation.reps, rep, persistent_alloc);
    Representation* clone = md_array_last(state->representation.reps);
    clone->md_rep = {0};
    clone->atom_mask = {0};
    init_representation(state, clone);
    update_representation(state, clone);
    return clone;
}

static void remove_representation(ApplicationState* state, int idx) {
    ASSERT(state);
    ASSERT(idx < md_array_size(state->representation.reps));
    auto& rep = state->representation.reps[idx];
    md_bitfield_free(&rep.atom_mask);
    md_gl_rep_destroy(rep.md_rep);
    if (rep.orbital.vol.tex_id) gl::free_texture(&rep.orbital.vol.tex_id);
    if (rep.orbital.dvr.tf_tex) gl::free_texture(&rep.orbital.dvr.tf_tex);
    md_array_swap_back_and_pop(state->representation.reps, idx);
}

static void recompute_atom_visibility_mask(ApplicationState* state) {
    ASSERT(state);

    auto& mask = state->representation.visibility_mask;

    md_bitfield_clear(&mask);
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        auto& rep = state->representation.reps[i];
        if (!rep.enabled) continue;
        md_bitfield_or_inplace(&mask, &rep.atom_mask);
    }
    state->representation.visibility_mask_hash = md_bitfield_hash64(&mask, 0);
}

static void update_all_representations(ApplicationState* state) {
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        auto& rep = state->representation.reps[i];
        rep.filt_is_dirty = true;
        update_representation(state, &rep);
    }
}

static bool rep_type_uses_colors(RepresentationType type) {
    switch(type) {
        case RepresentationType::DipoleMoment:
        case RepresentationType::Orbital:
            return false;
        default:
            return true;
    }
}

static void update_representation(ApplicationState* state, Representation* rep) {
    ASSERT(state);
    ASSERT(rep);

    if (!rep->enabled) return;

    md_vm_arena_temp_t tmp = md_vm_arena_temp_begin(frame_alloc);
    defer { md_vm_arena_temp_end(tmp); };

    const size_t bytes = state->mold.mol.atom.count * sizeof(uint32_t);

    const auto& mol = state->mold.mol;

    //md_script_property_t prop = {0};
    //if (rep->color_mapping == ColorMapping::Property) {
        //rep->prop_is_valid = md_script_compile_and_eval_property(&prop, rep->prop, &data->mold.mol, frame_allocator, &data->script.ir, rep->prop_error.beg(), rep->prop_error.capacity());
    //}

    bool use_colors = rep_type_uses_colors(rep->type);
    uint32_t* colors = 0;

    if (use_colors) {
        colors = (uint32_t*)md_vm_arena_push(frame_alloc, sizeof(uint32_t) * mol.atom.count);

        switch (rep->color_mapping) {
            case ColorMapping::Uniform:
                color_atoms_uniform(colors, mol.atom.count, rep->uniform_color);
                break;
            case ColorMapping::Cpk:
                color_atoms_cpk(colors, mol.atom.count, mol);
                break;
            case ColorMapping::AtomLabel:
                color_atoms_type(colors, mol.atom.count, mol);
                break;
            case ColorMapping::AtomIndex:
                color_atoms_idx(colors, mol.atom.count, mol);
                break;
            case ColorMapping::ResName:
                color_atoms_res_name(colors, mol.atom.count, mol);
                break;
            case ColorMapping::ResIndex:
                color_atoms_res_idx(colors, mol.atom.count, mol);
                break;
            case ColorMapping::ChainId:
                color_atoms_chain_id(colors, mol.atom.count, mol);
                break;
            case ColorMapping::ChainIndex:
                color_atoms_chain_idx(colors, mol.atom.count, mol);
                break;
            case ColorMapping::SecondaryStructure:
                color_atoms_sec_str(colors, mol.atom.count, mol);
                break;
            case ColorMapping::Property:
                // @TODO: Map colors accordingly
                //color_atoms_uniform(colors, mol.atom.count, rep->uniform_color);
    #if 0
                if (rep->prop) {
                    MEMSET(colors, 0xFFFFFFFF, bytes);
                    md_script_pro
                    const float* values = rep->prop->data.values;
                    if (rep->prop->data.aggregate) {
                        const int dim = rep->prop->data.dim[0];
                        md_script_vis_t vis = {0};
                        bool result = false;
                    
                    //if (md_semaphore_aquire(&data->script.ir_semaphore)) {
                    //    defer { md_semaphore_release(&data->script.ir_semaphore); };
                        
                        if (md_script_ir_valid(state->script.eval_ir)) {
                            md_script_vis_init(&vis, frame_alloc);
                            md_script_vis_ctx_t ctx = {
                                .ir = state->script.eval_ir,
                                .mol = &state->mold.mol,
                                .traj = state->mold.traj,
                            };
                            result = md_script_vis_eval_payload(&vis, rep->prop->vis_payload, 0, &ctx, MD_SCRIPT_VISUALIZE_ATOMS);
                        }
                    //}
                    if (result) {
                        if (dim == (int)md_array_size(vis.structure)) {
                            int i0 = CLAMP((int)state->animation.frame + 0, 0, (int)rep->prop->data.num_values / dim - 1);
                            int i1 = CLAMP((int)state->animation.frame + 1, 0, (int)rep->prop->data.num_values / dim - 1);
                            float frame_fract = fractf((float)state->animation.frame);

                            md_bitfield_t mask = {0};
                            md_bitfield_init(&mask, frame_alloc);
                            for (int i = 0; i < dim; ++i) {
                                md_bitfield_and(&mask, &rep->atom_mask, &vis.structure[i]);
                                float value = lerpf(values[i0 * dim + i], values[i1 * dim + i], frame_fract);
                                float t = CLAMP((value - rep->map_beg) / (rep->map_end - rep->map_beg), 0, 1);
                                ImVec4 color = ImPlot::SampleColormap(t, rep->color_map);
                                color_atoms_uniform(colors, mol.atom.count, vec_cast(color), &mask);
                            }
                        }
                    }
                } else {
                    int i0 = CLAMP((int)state->animation.frame + 0, 0, (int)rep->prop->data.num_values - 1);
                    int i1 = CLAMP((int)state->animation.frame + 1, 0, (int)rep->prop->data.num_values - 1);
                    float value = lerpf(values[i0], values[i1], fractf((float)state->animation.frame));
                    float t = CLAMP((value - rep->map_beg) / (rep->map_end - rep->map_beg), 0, 1);
                    ImVec4 color = ImPlot::SampleColormap(t, rep->color_map);
                    color_atoms_uniform(colors, mol.atom.count, vec_cast(color));
                }
            } else {
                color_atoms_uniform(colors, mol.atom.count, rep->uniform_color);
            }
#endif
            break;
        default:
            ASSERT(false);
            break;
        }
    }

    if (rep->saturation != 1.0f) {
        scale_saturation(colors, mol.atom.count, rep->saturation);
    }

    switch (rep->type) {
    case RepresentationType::SpaceFill:
        rep->type_is_valid = mol.atom.count > 0;
        break;
    case RepresentationType::Licorice:
        rep->type_is_valid = mol.bond.count > 0;
        break;
    case RepresentationType::BallAndStick:
        rep->type_is_valid = mol.atom.count > 0;
        break;
    case RepresentationType::Ribbons:
    case RepresentationType::Cartoon:
        rep->type_is_valid = mol.protein_backbone.range.count > 0;
        break;
    case RepresentationType::Orbital: {
        size_t num_mol_orbitals = md_array_size(state->representation.info.molecular_orbitals);
        rep->type_is_valid = num_mol_orbitals > 0;
        if (num_mol_orbitals > 0) {
            uint64_t vol_hash = (uint64_t)rep->orbital.type | ((uint64_t)rep->orbital.resolution << 8);
            if (rep->orbital.type == OrbitalType::MolecularOrbitalPsi || rep->orbital.type == OrbitalType::MolecularOrbitalPsiSquared) {
                vol_hash |= rep->orbital.orbital_idx;
            }
            if (vol_hash != rep->orbital.vol_hash) {
                bool success = false;
                if (rep->orbital.type == OrbitalType::AverageLocalIonizationEnergy) {
                    EvaluateVolumePayload eval_alie = {
                        .id = orbital_type_id[(int)OrbitalType::AverageLocalIonizationEnergy],
                        .aux = 0,
                        .samples_per_angstrom = volume_resolution_samples_per_angstrom[(int)rep->orbital.resolution],
                        .dst_volume = &rep->orbital.tf_vol,
                    };
                    viamd::event_system_broadcast_event(viamd::EventType_RepresentationVolumeEval, viamd::EventPayloadType_EvaluateVolume, &eval_alie);

                    EvaluateVolumePayload eval_den = {
                        .id = orbital_type_id[(int)OrbitalType::ElectronDensity],
                        .aux = 0,
                        .samples_per_angstrom = volume_resolution_samples_per_angstrom[(int)rep->orbital.resolution],
                        .dst_volume = &rep->orbital.vol,
                    };
                    viamd::event_system_broadcast_event(viamd::EventType_RepresentationVolumeEval, viamd::EventPayloadType_EvaluateVolume, &eval_den);

                    success = eval_alie.output_written && eval_den.output_written;
                } else {
                    EvaluateVolumePayload data = {
                        .id = orbital_type_id[(int)rep->orbital.type],
                        .aux = rep->orbital.orbital_idx,
                        .samples_per_angstrom = volume_resolution_samples_per_angstrom[(int)rep->orbital.resolution],
                        .dst_volume = &rep->orbital.vol,
                    };
                    viamd::event_system_broadcast_event(viamd::EventType_RepresentationVolumeEval, viamd::EventPayloadType_EvaluateVolume, &data);
                    success = data.output_written;
                }

                if (success) {
#if !VIAMD_RECOMPUTE_ORBITAL_PER_FRAME
                    rep->orbital.vol_hash = vol_hash;
#endif
                }
            }
        }
        {
            uint64_t tf_hash = md_hash64(&rep->orbital.alie.colormap, sizeof(rep->orbital.alie.colormap), md_hash64(&rep->orbital.alie.tf_alpha, sizeof(rep->orbital.alie.tf_alpha), 0));
            if (tf_hash != rep->orbital.alie.tf_hash) {
                rep->orbital.alie.tf_hash = tf_hash;
                volume::compute_transfer_function_texture_simple(&rep->orbital.alie.tf_tex, rep->orbital.alie.colormap, rep->orbital.alie.tf_alpha);
            }
        }
        {
            uint64_t tf_hash = md_hash64(&rep->orbital.dvr.colormap, sizeof(rep->orbital.dvr.colormap), 0);
            if (tf_hash != rep->orbital.dvr.tf_hash) {
                rep->orbital.dvr.tf_hash = tf_hash;
                volume::compute_transfer_function_texture_simple(&rep->orbital.dvr.tf_tex, rep->orbital.dvr.colormap);
            }
        }
        break;
    }
    case RepresentationType::DipoleMoment:
        rep->type_is_valid = md_array_size(state->representation.info.dipole_moments) > 0;
        break;
    default:
        ASSERT(false);
        break;
    }

    if (use_colors) {
        if (rep->dynamic_evaluation) {
            rep->filt_is_dirty = true;
        }

        if (rep->filt_is_dirty) {
            rep->filt_is_valid = filter_expression(state, str_from_cstr(rep->filt), &rep->atom_mask, &rep->filt_is_dynamic, rep->filt_error, sizeof(rep->filt_error));
            rep->filt_is_dirty = false;
        }

        if (rep->filt_is_valid) {
            filter_colors(colors, mol.atom.count, &rep->atom_mask);
            state->representation.atom_visibility_mask_dirty = true;
            md_gl_rep_set_color(rep->md_rep, 0, (uint32_t)mol.atom.count, colors, 0);

 #if EXPERIMENTAL_GFX_API
            md_gfx_rep_attr_t attributes = {};
            attributes.spacefill.radius_scale = 1.0f;
            md_gfx_rep_set_type_and_attr(rep->gfx_rep, MD_GFX_REP_TYPE_SPACEFILL, &attributes);
            md_gfx_rep_set_color(rep->gfx_rep, 0, (uint32_t)mol.atom.count, (md_gfx_color_t*)colors, 0);
#endif
        }
    }
}

static void init_representation(ApplicationState* state, Representation* rep) {
#if EXPERIMENTAL_GFX_API
    rep->gfx_rep = md_gfx_rep_create(state->mold.mol.atom.count);
#endif
    rep->md_rep = md_gl_rep_create(state->mold.gl_mol);
    md_bitfield_init(&rep->atom_mask, persistent_alloc);
    rep->filt_is_dirty = true;
}

static void update_representation_info(ApplicationState* state) {
    md_arena_allocator_reset(state->representation.info.alloc);
    state->representation.info.dipole_moments = 0;
    state->representation.info.molecular_orbitals = 0;
    state->representation.info.mo_homo_idx = 0;
    state->representation.info.mo_lumo_idx = 0;
    viamd::event_system_broadcast_event(viamd::EventType_RepresentationInfoFill, viamd::EventPayloadType_RepresentationInfo, &state->representation.info);
}

static void init_all_representations(ApplicationState* state) {
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        auto& rep = state->representation.reps[i];
        init_representation(state, &rep);
    }
}

static void clear_representations(ApplicationState* state) {
    while (md_array_size(state->representation.reps) > 0) {
        remove_representation(state, (int32_t)md_array_size(state->representation.reps) - 1);
    }
}

static void create_default_representations(ApplicationState* state) {
    bool amino_acid_present = false;
    bool nucleic_present = false;
    bool ion_present = false;
    bool water_present = false;
    bool ligand_present = false;
    bool orbitals_present = md_array_size(state->representation.info.molecular_orbitals) > 0;

    if (state->mold.mol.atom.count > 4'000'000) {
        LOG_INFO("Large molecule detected, creating default representation for all atoms");
        Representation* rep = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Cpk, STR_LIT("all"));
        snprintf(rep->name, sizeof(rep->name), "default");
        goto done;
    }

    if (state->mold.mol.residue.count == 0) {
        // No residues present
        Representation* rep = create_representation(state, RepresentationType::BallAndStick, ColorMapping::Cpk, STR_LIT("all"));
        snprintf(rep->name, sizeof(rep->name), "default");
        goto done;
    }

    for (size_t i = 0; i < state->mold.mol.atom.count; ++i) {
        uint32_t flags = state->mold.mol.atom.flags[i];
        if (flags & MD_FLAG_AMINO_ACID) amino_acid_present = true;
        if (flags & MD_FLAG_NUCLEOTIDE) nucleic_present = true;
        if (flags & MD_FLAG_ION) ion_present = true;
        if (flags & MD_FLAG_WATER) water_present = true;

        if (!(flags & (MD_FLAG_AMINO_ACID | MD_FLAG_NUCLEOTIDE | MD_FLAG_ION | MD_FLAG_WATER))) {
            ligand_present = true;
        }
    }

    if (amino_acid_present) {
        RepresentationType type = RepresentationType::Cartoon;
        ColorMapping color = ColorMapping::SecondaryStructure;

        if (state->mold.mol.chain.count > 1) {
            color = ColorMapping::ChainId;
        } else {
            size_t res_count = md_chain_residue_count(state->mold.mol.chain, 0);
            if (res_count < 20) {
                type = RepresentationType::BallAndStick;
                color = ColorMapping::Cpk;
            }
        }

        Representation* prot = create_representation(state, type, color, STR_LIT("protein"));
        snprintf(prot->name, sizeof(prot->name), "protein");
    }
    if (nucleic_present) {
        Representation* nucl = create_representation(state, RepresentationType::BallAndStick, ColorMapping::Cpk, STR_LIT("nucleic"));
        snprintf(nucl->name, sizeof(nucl->name), "nucleic");
    }
    if (ion_present) {
        Representation* ion = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Cpk, STR_LIT("ion"));
        snprintf(ion->name, sizeof(ion->name), "ion");
    }
    if (ligand_present) {
        Representation* ligand = create_representation(state, RepresentationType::BallAndStick, ColorMapping::Cpk, STR_LIT("not (protein or nucleic or water or ion)"));
        snprintf(ligand->name, sizeof(ligand->name), "ligand");
    }
    if (water_present) {
        Representation* water = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Cpk, STR_LIT("water"));
        water->scale.x = 0.5f;
        snprintf(water->name, sizeof(water->name), "water");
        water->enabled = false;
        if (!amino_acid_present && !nucleic_present && !ligand_present) {
            water->enabled = true;
        }
    }

done:
    if (orbitals_present) {
        Representation* rep = create_representation(state, RepresentationType::Orbital);
        snprintf(rep->name, sizeof(rep->name), "orbital");
        rep->enabled = true;
    }

    recompute_atom_visibility_mask(state);
}

// #selection
static Selection* create_selection(ApplicationState* state, str_t name, md_bitfield_t* atom_mask) {
    ASSERT(state);
    Selection sel;
    str_copy_to_char_buf(sel.name, sizeof(sel.name), name);
    md_bitfield_init(&sel.atom_mask, persistent_alloc);
    if (atom_mask) {
        md_bitfield_copy(&sel.atom_mask, atom_mask);
    }
    md_array_push(state->selection.stored_selections, sel, persistent_alloc);
    return md_array_last(state->selection.stored_selections);
}

static void clear_selections(ApplicationState* state) {
    md_array_shrink(state->selection.stored_selections, 0);
}

#if 0
    static Selection* clone_selection(ApplicationState * data, const Selection& src) {
        ASSERT(data);
        Selection clone;
        clone.name = src.name;
        bitfield::init(&clone.atom_mask, data->selection.selection_mask);
        return &data->selection.stored_selections.push_back(clone);
    }
#endif

static void remove_selection(ApplicationState* state, int idx) {
    ASSERT(state);
    if (idx < 0 || (int)md_array_size(state->selection.stored_selections) <= idx) {
        LOG_ERROR("Index [%i] out of range when trying to remove selection", idx);
    }
    auto item = &state->selection.stored_selections[idx];
    md_bitfield_free(&item->atom_mask);
    
    state->selection.stored_selections[idx] = *md_array_last(state->selection.stored_selections);
    md_array_pop(state->selection.stored_selections);
}

// #camera-control
static void handle_camera_interaction(ApplicationState* data) {
    ASSERT(data);

    enum class RegionMode { Append, Remove };

#if 1
    // Coordinate system widget

    {
        static bool locked = true;
        const int def_size = 150;
        const int min_size = 20;
        const int max_size = 500;

        ImGui::SetNextWindowSize(ImVec2(def_size, def_size), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSizeConstraints(ImVec2(min_size, min_size), ImVec2(max_size, max_size), [](ImGuiSizeCallbackData* data) {
            // Enforce a square aspect ratio
            data->DesiredSize = ImVec2(ImMax(data->DesiredSize.x, data->DesiredSize.y), ImMax(data->DesiredSize.x, data->DesiredSize.y));
        });

        // Scale down backround color to only show a faint outline
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImGui::GetStyleColorVec4(ImGuiCol_WindowBg) * ImVec4(1,1,1,0.1f));
        defer { ImGui::PopStyleColor(); };

        int flags = ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoDocking;
        if (locked) {
            flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoBackground | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoNavFocus;
        }

        if (ImGui::Begin("Coordinate Widget", NULL, flags)) {
            if (ImGui::IsWindowHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
                ImGui::OpenPopup("Coordinate Widget Context");
            }

            if (ImGui::BeginPopup("Coordinate Widget Context")) {
                ImGui::Checkbox("Locked", &locked);
                ImGui::EndPopup();
            }

            CoordSystemWidgetParam param = {
                .pos = ImGui::GetWindowContentRegionMin(),
                .size = ImGui::GetContentRegionAvail(),
                .view_matrix = data->view.param.matrix.curr.view,
                .camera_ori = data->view.animation.target_orientation,
                .camera_pos = data->view.animation.target_position,
                .camera_dist = data->view.animation.target_distance,
            };

            ImGui::DrawCoordinateSystemWidget(param);
        }
        ImGui::End();
    }
#endif

    ImGui::BeginCanvas("Main interaction window", true);
    bool pressed = ImGui::InvisibleButton("canvas", ImGui::GetContentRegionAvail(), ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);
    bool open_atom_context = false;

    ImGuiWindow* window = ImGui::GetCurrentWindow();
    if (window) {
        ImDrawList* dl = window->DrawList;
        ASSERT(dl);

        if (pressed || ImGui::IsItemActive() || ImGui::IsItemDeactivated()) {
            if (ImGui::IsKeyPressed(ImGuiMod_Shift, false)) {
                ImGui::ResetMouseDragDelta();
            }

            if (ImGui::IsKeyDown(ImGuiMod_Shift)) {
                RegionMode mode = RegionMode::Append;
                if (ImGui::IsMouseDown(ImGuiMouseButton_Left) || ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                    mode = RegionMode::Append;
                }
                else if (ImGui::IsMouseDown(ImGuiMouseButton_Right) || ImGui::IsMouseReleased(ImGuiMouseButton_Right)) {
                    mode = RegionMode::Remove;
                }

                const ImVec2 ext = ImGui::GetMouseDragDelta(mode == RegionMode::Append ? ImGuiMouseButton_Left : ImGuiMouseButton_Right);
                const ImVec2 pos = ImGui::GetMousePos() - ext;
                const ImU32 fill_col = 0x22222222;
                const ImU32 line_col = 0x88888888;

                ASSERT(dl);
                dl->AddRectFilled(pos, pos + ext, fill_col);
                dl->AddRect(pos, pos + ext, line_col);

                const ImVec2 min_p = ImMin(pos, pos + ext) - window->Pos;
                const ImVec2 max_p = ImMax(pos, pos + ext) - window->Pos;

                md_bitfield_t mask = { 0 };
                md_bitfield_init(&mask, frame_alloc);

                if (min_p != max_p) {
                    md_bitfield_clear(&data->selection.highlight_mask);
                    data->selection.selecting = true;

                    const vec2_t res = { (float)data->app.window.width, (float)data->app.window.height };
                    const mat4_t mvp = data->view.param.matrix.curr.proj * data->view.param.matrix.curr.view;

                    md_bitfield_iter_t it = md_bitfield_iter_create(&data->representation.visibility_mask);
                    while (md_bitfield_iter_next(&it)) {
                        const uint64_t i = md_bitfield_iter_idx(&it);
                        const vec4_t p = mat4_mul_vec4(mvp, vec4_set(data->mold.mol.atom.x[i], data->mold.mol.atom.y[i], data->mold.mol.atom.z[i], 1.0f));
                        const vec2_t c = {
                            ( p.x / p.w * 0.5f + 0.5f) * res.x,
                            (-p.y / p.w * 0.5f + 0.5f) * res.y,
                        };

                        if (min_p.x <= c.x && c.x <= max_p.x && min_p.y <= c.y && c.y <= max_p.y) {
                            md_bitfield_set_bit(&mask, i);
                        }
                    }
                    grow_mask_by_selection_granularity(&mask, data->selection.granularity, data->mold.mol);

                    if (mode == RegionMode::Append) {
                        md_bitfield_or(&data->selection.highlight_mask, &data->selection.selection_mask, &mask);
                    }
                    else if (mode == RegionMode::Remove) {
                        md_bitfield_andnot(&data->selection.highlight_mask, &data->selection.selection_mask, &mask);
                    }
                    if (pressed || ImGui::IsMouseReleased(0)) {
                        md_bitfield_copy(&data->selection.selection_mask, &data->selection.highlight_mask);
                    }
                }
                else if (pressed) {
                    if (data->selection.atom_idx.hovered != -1 || data->selection.bond_idx.hovered != -1) {
                        if (mode == RegionMode::Append) {
                            if (data->selection.atom_idx.hovered != -1) {
                                single_selection_sequence_push_idx(&data->selection.single_selection_sequence, data->selection.atom_idx.hovered);
                                md_bitfield_set_bit(&mask, data->selection.atom_idx.hovered);
                            } else {
                                md_bond_pair_t pair = data->mold.mol.bond.pairs[data->selection.bond_idx.hovered];
                                md_bitfield_set_bit(&mask, pair.idx[0]);
                                md_bitfield_set_bit(&mask, pair.idx[1]);
                            }
                            grow_mask_by_selection_granularity(&mask, data->selection.granularity, data->mold.mol);
                            md_bitfield_or_inplace(&data->selection.selection_mask, &mask);
                        }
                        else if (mode == RegionMode::Remove) {
                            if (data->selection.atom_idx.hovered != -1) {
                                single_selection_sequence_pop_idx(&data->selection.single_selection_sequence, data->selection.atom_idx.hovered);
                                md_bitfield_set_bit(&mask, data->selection.atom_idx.hovered);
                            } else {
                                md_bond_pair_t pair = data->mold.mol.bond.pairs[data->selection.bond_idx.hovered];
                                md_bitfield_set_bit(&mask, pair.idx[0]);
                                md_bitfield_set_bit(&mask, pair.idx[1]);
                            }
                            grow_mask_by_selection_granularity(&mask, data->selection.granularity, data->mold.mol);
                            md_bitfield_andnot_inplace(&data->selection.selection_mask, &mask);
                        }
                    }
                    else {
                        single_selection_sequence_clear(&data->selection.single_selection_sequence);
                        md_bitfield_clear(&data->selection.selection_mask);
                        md_bitfield_clear(&data->selection.highlight_mask);
                    }
                }
            }
        }
        else if (ImGui::IsItemHovered() && !ImGui::IsAnyItemActive()) {
            md_bitfield_clear(&data->selection.highlight_mask);
            if (data->picking.idx != INVALID_PICKING_IDX) {
                if (data->selection.atom_idx.hovered != -1 && data->mold.mol.atom.count) {
                    md_bitfield_set_bit(&data->selection.highlight_mask, data->picking.idx);
                }
                else if (data->selection.bond_idx.hovered != -1 && data->selection.bond_idx.hovered < (int32_t)data->mold.mol.bond.count) {
                    md_bond_pair_t pair = data->mold.mol.bond.pairs[data->selection.bond_idx.hovered];
                    md_bitfield_set_bit(&data->selection.highlight_mask, pair.idx[0]);
                    md_bitfield_set_bit(&data->selection.highlight_mask, pair.idx[1]);
                }
                grow_mask_by_selection_granularity(&data->selection.highlight_mask, data->selection.granularity, data->mold.mol);

                draw_info_window(*data, data->picking.idx);
            }
        }

        if (ImGui::IsItemActive() || ImGui::IsItemHovered()) {
            if (!ImGui::IsKeyDown(ImGuiMod_Shift) && !data->selection.selecting) {
                const ImVec2 delta = ImGui::GetIO().MouseDelta;
                const ImVec2 coord = ImGui::GetMousePos() - ImGui::GetCurrentWindow()->Pos;
                const vec2_t mouse_delta = {delta.x, delta.y};
                const vec2_t mouse_coord = {coord.x, coord.y};
                const float  scroll_delta = ImGui::GetIO().MouseWheel;

                TrackballControllerInput input;
                input.rotate_button = ImGui::IsMouseDown(ImGuiMouseButton_Left);
                input.pan_button    = ImGui::IsMouseDown(ImGuiMouseButton_Right);
                input.dolly_button  = ImGui::IsMouseDown(ImGuiMouseButton_Middle);
                input.mouse_coord_curr = mouse_coord;
                input.mouse_coord_prev = mouse_coord - mouse_delta;
                input.screen_size = {(float)data->app.window.width, (float)data->app.window.height};
                input.dolly_delta = scroll_delta;
                input.fov_y = data->view.camera.fov_y;

                vec3_t pos = data->view.animation.target_position;
                quat_t ori = data->view.animation.target_orientation;
                float dist = data->view.animation.target_distance;
            
                TrackballFlags flags = TrackballFlags_AnyInteractionReturnsTrue;
                if (ImGui::IsItemActive()) {
                    flags |= TrackballFlags_EnableAllInteractions;
                } else {
                    flags |= TrackballFlags_DollyEnabled;
                }

                if (camera_controller_trackball(&pos, &ori, &dist, input, data->view.trackball_param, flags)) {
                    // @NOTE(Robin): We could make the camera interaction more snappy, by directly modifying camera[pos, ori, dist] here
                    // But for now its decent. I like smooth transitions rather than discontinous jumps
                }
                data->view.animation.target_position = pos;
                data->view.animation.target_orientation = ori;
                data->view.animation.target_distance = dist;

                if (ImGui::GetIO().MouseDoubleClicked[0]) {
                    if (data->picking.depth < 1.0f) {
                        const vec3_t forward = data->view.camera.orientation * vec3_t{0, 0, 1};
                        data->view.animation.target_position = data->picking.world_coord + forward * dist;
                    } else {
                        reset_view(data, &data->representation.visibility_mask, true, true);
                    }
                }

                data->visuals.dof.focus_depth = data->view.camera.focus_distance;
            
                if (ImGui::GetMouseDragDelta(ImGuiMouseButton_Right) == ImVec2(0,0) &&
                    ImGui::IsMouseReleased(ImGuiMouseButton_Right))
                {
                    open_atom_context = true;
                }
            }
        }
    }

    ImGui::EndCanvas();

    if (open_atom_context) {
        ImGui::OpenPopup("AtomContextPopup");
    }
}

static void fill_gbuffer(ApplicationState* data) {
    const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
        GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_TRANSPARENCY };

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    // Enable all draw buffers
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, data->gbuffer.fbo);
    glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);

    PUSH_GPU_SECTION("G-Buffer fill")

    // Immediate mode graphics

    if (data->simulation_box.enabled && data->mold.mol.unit_cell.basis != mat3_t{0}) {
        PUSH_GPU_SECTION("Draw Simulation Box")
        const mat4_t model_mat = mat4_from_mat3(data->mold.mol.unit_cell.basis);
        const mat4_t model_view_mat = data->view.param.matrix.curr.view * model_mat;
        immediate::set_model_view_matrix(model_view_mat);
        immediate::set_proj_matrix(data->view.param.matrix.curr.proj);
        immediate::draw_box_wireframe({0,0,0}, {1,1,1}, convert_color(data->simulation_box.color));
        immediate::render();
        POP_GPU_SECTION()
    }

#if 0
    // RENDER DEBUG INFORMATION (WITH DEPTH)
    PUSH_GPU_SECTION("Debug Draw") {
        glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
        immediate::set_model_view_matrix(data->view.param.matrix.curr.view);
        immediate::set_proj_matrix(data->view.param.matrix.curr.proj);
        immediate::flush();
    }
    POP_GPU_SECTION()

    PUSH_GPU_SECTION("Debug Draw Overlay") {
        glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);  // Post_Tonemap buffer
        glDisable(GL_DEPTH_TEST);
        glDepthMask(0);

        // immediate::set_model_view_matrix(data->view.param.matrix.current.view);
        // immediate::set_proj_matrix(data->view.param.matrix.current.proj);
        // immediate::flush();

        glEnable(GL_DEPTH_TEST);
        glDepthMask(1);
    }
    POP_GPU_SECTION()
#endif

    if (!use_gfx) {
        // DRAW VELOCITY OF STATIC OBJECTS
        PUSH_GPU_SECTION("Blit Static Velocity")
        glDrawBuffer(GL_COLOR_ATTACHMENT_VELOCITY);
        glDepthMask(0);
        postprocessing::blit_static_velocity(data->gbuffer.tex.depth, data->view.param);
        glDepthMask(1);
        POP_GPU_SECTION()
    }
    glDepthMask(1);
    glColorMask(1, 1, 1, 1);

    // DRAW REPRESENTATIONS
    PUSH_GPU_SECTION("Draw Opaque")
    glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
    draw_representations_opaque(data);
    viamd::event_system_broadcast_event(viamd::EventType_ViamdRenderOpaque, viamd::EventPayloadType_ApplicationState, data);
    POP_GPU_SECTION()

    glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);

    if (!use_gfx) {
        PUSH_GPU_SECTION("Selection")
        const bool atom_selection_empty = md_bitfield_popcount(&data->selection.selection_mask) == 0;
        const bool atom_highlight_empty = md_bitfield_popcount(&data->selection.highlight_mask) == 0;

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
            draw_representations_opaque_lean_and_mean(data, AtomBit_Selected | AtomBit_Visible);

            glDisable(GL_DEPTH_TEST);

            glStencilMask(0x0);
            glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
            glColorMask(1, 1, 1, 1);

            glStencilFunc(GL_EQUAL, 2, 0xFF);
            postprocessing::blit_color(data->selection.color.selection.visible);

            glStencilFunc(GL_EQUAL, 0, 0xFF);
            postprocessing::blit_color(data->selection.color.selection.hidden);
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
            draw_representations_opaque_lean_and_mean(data, AtomBit_Highlighted | AtomBit_Visible);

            glDisable(GL_DEPTH_TEST);

            glStencilMask(0x0);
            glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
            glColorMask(1, 1, 1, 1);

            glStencilFunc(GL_EQUAL, 2, 0xFF);
            postprocessing::blit_color(data->selection.color.highlight.visible);

            glStencilFunc(GL_EQUAL, 0, 0xFF);
            postprocessing::blit_color(data->selection.color.highlight.hidden);
        }

        glDisable(GL_STENCIL_TEST);

        if (!atom_selection_empty) {
            PUSH_GPU_SECTION("Desaturate") {
                const float saturation = data->selection.color.saturation;
                glDrawBuffer(GL_COLOR_ATTACHMENT_COLOR);
                postprocessing::scale_hsv(data->gbuffer.tex.color, vec3_t{1, saturation, 1});
            } POP_GPU_SECTION()
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

    PUSH_GPU_SECTION("Draw Transparent")
    draw_representations_transparent(data);
    viamd::event_system_broadcast_event(viamd::EventType_ViamdRenderTransparent, viamd::EventPayloadType_ApplicationState, data);
    POP_GPU_SECTION()

    PUSH_GPU_SECTION("Draw Visualization Geometry")
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    immediate::set_model_view_matrix(data->view.param.matrix.curr.view);
    immediate::set_proj_matrix(data->view.param.matrix.curr.proj);

    const md_script_vis_t& vis = data->script.vis;

    const uint32_t point_color      = convert_color(data->script.point_color);
    const uint32_t line_color       = convert_color(data->script.line_color);
    const uint32_t triangle_color   = convert_color(data->script.triangle_color);

    md_array(mat4_t) model_matrices = 0;
    if (md_array_size(vis.sdf.matrices) > 0) {
        model_matrices = md_array_create(mat4_t, md_array_size(vis.sdf.matrices), frame_alloc);
        for (size_t i = 0; i < md_array_size(vis.sdf.matrices); ++i) {
            model_matrices[i] = mat4_inverse(vis.sdf.matrices[i]);
        }
    }

    if (vis.points) {
        immediate::draw_points_v((immediate::Vertex*)vis.points, md_array_size(vis.points), data->script.point_color);
    }

    if (vis.triangles) {
        immediate::draw_triangles_v((immediate::Vertex*)vis.triangles, md_array_size(vis.triangles), data->script.triangle_color);
    }

    if (vis.lines) {
        immediate::draw_lines_v((immediate::Vertex*)vis.lines, md_array_size(vis.lines), data->script.line_color);
    }

    

    

    const vec4_t col_x = {1, 0, 0, 0.7f};
    const vec4_t col_y = {0, 1, 0, 0.7f};
    const vec4_t col_z = {0, 0, 1, 0.7f};
    const float ext = vis.sdf.extent * 0.25f;
    for (size_t i = 0; i < md_array_size(model_matrices); ++i) {
        immediate::draw_basis(model_matrices[i], ext, col_x, col_y, col_z);
    }
    
    immediate::render();

    immediate::set_model_view_matrix(data->view.param.matrix.curr.view);
    immediate::set_proj_matrix(data->view.param.matrix.curr.proj);
    
    glEnable(GL_DEPTH_TEST);
    
    const vec3_t box_ext = vec3_set1(vis.sdf.extent);
    for (size_t i = 0; i < md_array_size(model_matrices); ++i) {
        immediate::draw_box_wireframe(-box_ext, box_ext, model_matrices[i], data->density_volume.bounding_box_color);
    }

    immediate::render();

    md_array_free(model_matrices, frame_alloc);

    md_script_vis_free(&data->script.vis);

    glEnable(GL_CULL_FACE);
    POP_GPU_SECTION()

    POP_GPU_SECTION()  // G-buffer
}

static void handle_picking(ApplicationState* data) {
    ImGuiContext* ctx = ImGui::GetCurrentContext();
    ImGuiWindow* window = ImGui::FindWindowByName("Main interaction window");

    if (ctx && window && ctx->HoveredWindow == window) {
        PUSH_CPU_SECTION("PICKING")
        vec2_t mouse_pos = vec_cast(ImGui::GetMousePos() - ImGui::GetMainViewport()->Pos);
        vec2_t coord = {mouse_pos.x, ImGui::GetMainViewport()->Size.y - mouse_pos.y};
#if MD_PLATFORM_OSX
        coord = coord * vec_cast(ImGui::GetIO().DisplayFramebufferScale);
#endif
        if (coord.x < 0.f || coord.x >= (float)data->gbuffer.width || coord.y < 0.f || coord.y >= (float)data->gbuffer.height) {
            data->picking.idx = INVALID_PICKING_IDX;
            data->picking.depth = 1.f;
        } else {
#if PICKING_JITTER_HACK
            static uint32_t frame_idx = 0;
            static uint32_t ref_frame = 0;
            frame_idx = (frame_idx + 1) % JITTER_SEQUENCE_SIZE;
            // @NOTE: If we have jittering applied, we cannot? retreive the original pixel value (without the jitter)
            // Solution, pick one reference frame out of the jittering sequence and use that one...
            // Ugly hack but works...

            if (data->app.input.mouse.moving) {
                ref_frame = frame_idx;
            }

            if (ref_frame == frame_idx || data->view.param.jitter.curr == vec2_t{0, 0}) {
                data->picking = read_picking_data(data->gbuffer, (int32_t)round(coord.x), (int32_t)round(coord.y));
                if (data->picking.idx != INVALID_PICKING_IDX)
                    data->picking.idx = CLAMP(data->picking.idx, 0U, (uint32_t)data->mold.mol.atom.count - 1U);
                const vec4_t viewport = {0, 0, (float)data->gbuffer.width, (float)data->gbuffer.height};
                data->picking.world_coord = mat4_unproject({coord.x, coord.y, data->picking.depth}, data->view.param.matrix.inv.view_proj, viewport);
            }
#else
            data->picking = {};
            extract_picking_data(&data->picking.idx, &data->picking.depth, &data->gbuffer, (int)coord.x, (int)coord.y);
            const vec4_t viewport = {0, 0, (float)data->gbuffer.width, (float)data->gbuffer.height};
            const mat4_t inv_VP = data->view.param.matrix.inv.view * data->view.param.matrix.inv.proj;
            data->picking.world_coord = mat4_unproject({coord.x, coord.y, data->picking.depth}, inv_VP, viewport);
            data->picking.screen_coord = coord;
#endif
        }
        data->selection.atom_idx.hovered = -1;
        data->selection.bond_idx.hovered = -1;

        if (data->picking.idx != INVALID_PICKING_IDX) {
            // The index space is segmented into two parts, the first half is for atoms and the second half is for bonds
            if (data->picking.idx < 0x80000000) {
                data->selection.atom_idx.hovered = data->picking.idx;
            } else {
                data->selection.bond_idx.hovered = data->picking.idx & 0x7FFFFFFF;
            }
        }
        
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
            data->selection.atom_idx.right_click = data->selection.atom_idx.hovered;
            data->selection.bond_idx.right_click = data->selection.bond_idx.hovered;
        }
        POP_CPU_SECTION()
    }
}

static void apply_postprocessing(const ApplicationState& data) {
    PUSH_GPU_SECTION("Postprocessing")
    postprocessing::Descriptor desc;

    desc.background.color = data.visuals.background.color * data.visuals.background.intensity;

    desc.ambient_occlusion.enabled = data.visuals.ssao.enabled;
    desc.ambient_occlusion.intensity = data.visuals.ssao.intensity;
    desc.ambient_occlusion.radius = data.visuals.ssao.radius;
    desc.ambient_occlusion.bias = data.visuals.ssao.bias;

    desc.tonemapping.enabled = data.visuals.tonemapping.enabled;
    desc.tonemapping.mode = data.visuals.tonemapping.tonemapper;
    desc.tonemapping.exposure = data.visuals.tonemapping.exposure;
    desc.tonemapping.gamma = data.visuals.tonemapping.gamma;

    desc.depth_of_field.enabled = data.visuals.dof.enabled;
    desc.depth_of_field.focus_depth = data.visuals.dof.focus_depth;
    desc.depth_of_field.focus_scale = data.visuals.dof.focus_scale;

    desc.fxaa.enabled = data.visuals.fxaa.enabled;

    constexpr float MOTION_BLUR_REFERENCE_DT = 1.0f / 60.0f;
    const float dt_compensation = MOTION_BLUR_REFERENCE_DT / (float)data.app.timing.delta_s;
    const float motion_scale = data.visuals.temporal_aa.motion_blur.motion_scale * dt_compensation;
    desc.temporal_aa.enabled = data.visuals.temporal_aa.enabled;
    desc.temporal_aa.feedback_min = data.visuals.temporal_aa.feedback_min;
    desc.temporal_aa.feedback_max = data.visuals.temporal_aa.feedback_max;
    desc.temporal_aa.motion_blur.enabled = data.visuals.temporal_aa.motion_blur.enabled;
    desc.temporal_aa.motion_blur.motion_scale = motion_scale;

    desc.sharpen.enabled = data.visuals.temporal_aa.enabled && data.visuals.sharpen.enabled;
    desc.sharpen.weight  = data.visuals.sharpen.weight;

    desc.input_textures.depth = data.gbuffer.tex.depth;
    desc.input_textures.color = data.gbuffer.tex.color;
    desc.input_textures.normal = data.gbuffer.tex.normal;
    desc.input_textures.velocity = data.gbuffer.tex.velocity;
    desc.input_textures.transparency = data.gbuffer.tex.transparency;

    postprocessing::shade_and_postprocess(desc, data.view.param);
    POP_GPU_SECTION()
}

static void draw_representations_opaque(ApplicationState* data) {
    ASSERT(data);

    if (data->mold.mol.atom.count == 0) {
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
                md_range_t range = {0, (int32_t)data->mold.mol.atom.count};
                md_array_push(transforms, M, persistent_alloc);
            }
        }

        md_gfx_draw_op_t* draw_ops = 0;
        for (int64_t i = 0; i < md_array_size(data->representation.reps); ++i) {
            if (data->representation.reps[i].enabled) {
                md_gfx_draw_op_t op;
                op.structure = data->mold.gfx_structure;
                op.representation = data->representation.reps[i].gfx_rep;
                op.model_mat = NULL;
                md_array_push(draw_ops, op, frame_alloc);
                
                for (uint32_t j = 0; j < instance_count; ++j) {
                    md_gfx_draw_op_t op;
                    op.structure = data->mold.gfx_structure;
                    op.representation = data->representation.reps[i].gfx_rep;
                    op.model_mat = &transforms[j];
                    md_array_push(draw_ops, op, frame_alloc);
                }
                
            }
        }

        md_gfx_draw((uint32_t)md_array_size(draw_ops), draw_ops, &data->view.param.matrix.curr.proj, &data->view.param.matrix.curr.view, &data->view.param.matrix.inv.proj, &data->view.param.matrix.inv.view);
    } else {
#endif
        const size_t num_representations = md_array_size(data->representation.reps);
        if (num_representations == 0) return;

        md_array(md_gl_draw_op_t) draw_ops = 0;
        for (size_t i = 0; i < num_representations; ++i) {
            const Representation& rep = data->representation.reps[i];

            if (rep.type > RepresentationType::Cartoon) continue;

            if (rep.enabled && rep.type_is_valid) {
                md_gl_draw_op_t op = {
                    .type = (md_gl_rep_type_t)rep.type,
                    .args = {},
                    .rep = data->representation.reps[i].md_rep,
                    .model_matrix = NULL,
                };
                MEMCPY(&op.args, &rep.scale, sizeof(op.args));
                md_array_push(draw_ops, op, frame_alloc);
            }
        }

        md_gl_draw_args_t args = {
            .shaders = data->mold.gl_shaders,
            .draw_operations = {
                .count = (uint32_t)md_array_size(draw_ops),
                .ops = draw_ops,
            },
            .view_transform = {
                .view_matrix = &data->view.param.matrix.curr.view.elem[0][0],
                .proj_matrix = &data->view.param.matrix.curr.proj.elem[0][0],
                // These two are for temporal anti-aliasing reprojection (optional)
                .prev_view_matrix = &data->view.param.matrix.prev.view.elem[0][0],
                .prev_proj_matrix = &data->view.param.matrix.prev.proj.elem[0][0],
            },
        };

        md_gl_draw(&args);
#if EXPERIMENTAL_GFX_API
    }
#endif
}

static void draw_representations_transparent(ApplicationState* state) {
    ASSERT(state);
    if (state->mold.mol.atom.count == 0) return;

    const size_t num_representations = md_array_size(state->representation.reps);
    if (num_representations == 0) return;

    for (size_t i = 0; i < num_representations; ++i) {
        const Representation& rep = state->representation.reps[i];
        if (!rep.enabled) continue;
        if (rep.type != RepresentationType::Orbital) continue;

        const IsoDesc& iso = (rep.orbital.type == OrbitalType::MolecularOrbitalPsi) ? rep.orbital.iso_psi : rep.orbital.iso_den;

#if VIAMD_RECOMPUTE_ORBITAL_PER_FRAME
        update_representation(state, &state->representation.reps[i]);
#endif

        uint32_t tf_vol = 0;
        uint32_t tf_tex = 0;
        float tf_min = 0;
        float tf_max = 0;
        if (rep.orbital.type == OrbitalType::AverageLocalIonizationEnergy) {
            tf_vol = rep.orbital.tf_vol.tex_id;
            tf_tex = rep.orbital.alie.tf_tex;
            tf_min = rep.orbital.alie.min_tf_val;
            tf_max = rep.orbital.alie.max_tf_val;
        }

        volume::RenderDesc desc = {
            .render_target = {
                .depth = state->gbuffer.tex.depth,
                .color  = state->gbuffer.tex.transparency,
                .width  = state->gbuffer.width,
                .height = state->gbuffer.height,
            },
            .texture = {
                .volume = rep.orbital.vol.tex_id,
                .tf_volume = tf_vol,
                .tf = tf_tex,
            },
            .matrix = {
                .model = rep.orbital.vol.texture_to_world,
                .view  = state->view.param.matrix.curr.view,
                .proj  = state->view.param.matrix.curr.proj,
                .inv_proj = state->view.param.matrix.inv.proj,
            },
            .clip_volume = {
                .min = {0,0,0},
                .max = {1,1,1},
            },
            .temporal = {
                .enabled = state->visuals.temporal_aa.enabled,
            },
            .iso = {
                .enabled = iso.enabled,
                .count   = iso.count,
                .values  = iso.values,
                .colors  = iso.colors,
            },
            .dvr = {
                .enabled = false,
            },
            .tf = {
                .min_value = tf_min,
                .max_value = tf_max,
            },
            .shading = {
                .env_radiance = state->visuals.background.color * state->visuals.background.intensity * 0.25f,
                .roughness = 0.3f,
                .dir_radiance = {10,10,10},
                .ior = 1.5f,
        },
            .voxel_spacing = rep.orbital.vol.step_size,
        };

        volume::render_volume(desc);

#if DEBUG
        immediate::set_model_view_matrix(state->view.param.matrix.curr.view);
        immediate::set_proj_matrix(state->view.param.matrix.curr.proj);
        immediate::draw_box_wireframe({0,0,0}, {1,1,1}, rep.orbital.vol.texture_to_world, immediate::COLOR_BLACK);
        immediate::render();
#endif
    }
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
            MEMCPY(&op.args, &rep.scale, sizeof(op.args));
            md_array_push(draw_ops, op, frame_alloc);
        }
    }

    md_gl_draw_args_t args = {
        .shaders = data->mold.gl_shaders_lean_and_mean,
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
        .atom_mask = mask,
    };

    md_gl_draw(&args);
}
