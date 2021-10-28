#include <core/md_compiler.h>

#if MD_COMPILER_MSVC

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#pragma warning( disable : 26812 4244 )
#endif

#include <md_util.h>
#include <md_gl.h>
#include <md_filter.h>
#include <md_script.h>
#include <md_molecule.h>
#include <md_trajectory.h>

#include <core/md_sync.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_stack_allocator.h>
#include <core/md_tracking_allocator.h>
#include <core/md_log.h>
#include <core/md_simd.h>
#include <core/md_file.h>
#include <core/md_array.inl>
#include <core/md_os.h>
#include <core/md_spatial_hash.h>
#include <core/md_base64.h>

#include <imgui.h>
#define IMGUI_DEFINE_MATH_OPERATORS
#include <imgui_internal.h>

#include <implot.h>
#include "implot_internal.h"

#include <TextEditor.h>

#include "gfx/gl.h"
#include "gfx/gl_utils.h"
#include "gfx/camera.h"
#include "gfx/camera_utils.h"
#include "gfx/immediate_draw_utils.h"
#include "gfx/postprocessing_utils.h"
#include "gfx/volumerender_utils.h"
#include "gfx/conetracing_utils.h"

#include "string_util.h"
#include "random_util.h"
#include "imgui_widgets.h"
#include "implot_widgets.h"
#include "application/application.h"
#include "console.h"
#include "color_utils.h"
#include "isosurface.h"
#include "task_system.h"
#include "loader.h"

#include <atomic_queue.h>
#include <stdio.h>

#define PICKING_JITTER_HACK 0
#define SHOW_IMGUI_DEMO_WINDOW 0
#define EXPERIMENTAL_CONE_TRACED_AO 0
#define COMPILATION_TIME_DELAY_IN_SECONDS 1.0
#define IR_SEMAPHORE_MAX_COUNT 3

#if MD_PLATFORM_OSX
const Key::Key_t KEY_CONSOLE = Key::KEY_WORLD_1;
#else  // WIN32 and Linux
// @TODO: Make sure this is currect for Linux?
const Key::Key_t KEY_CONSOLE = Key::KEY_GRAVE_ACCENT;
#endif

const Key::Key_t KEY_PLAY_PAUSE = Key::KEY_SPACE;
const Key::Key_t KEY_SKIP_TO_PREV_FRAME = Key::KEY_LEFT;
const Key::Key_t KEY_SKIP_TO_NEXT_FRAME = Key::KEY_RIGHT;
const Key::Key_t KEY_TOGGLE_SCREENSHOT_MODE = Key::KEY_F10;
const Key::Key_t KEY_SHOW_DEBUG_WINDOW = Key::KEY_F11;

// For cpu profiling
#define PUSH_CPU_SECTION(lbl) {};
#define POP_CPU_SECTION() {};

// For gpu profiling
#define PUSH_GPU_SECTION(lbl)                                                                   \
{                                                                                               \
    if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); \
}
#define POP_GPU_SECTION()                   \
{                                           \
    if (glPopDebugGroup) glPopDebugGroup(); \
}

constexpr str_t FILE_EXTENSION = make_cstr("via"); 
constexpr uint32_t INVALID_PICKING_IDX = ~0U;

constexpr uint32_t TEXT_BG_ERROR_COLOR = 0xAA222299;
constexpr uint32_t PROPERTY_COLORS[] = {4293119554, 4290017311, 4287291314, 4281114675, 4288256763, 4280031971, 4285513725, 4278222847, 4292260554, 4288298346, 4288282623, 4280834481};

inline const ImVec4& vec_cast(const vec4_t& v) { return *(const ImVec4*)(&v); }
inline const vec4_t& vec_cast(const ImVec4& v) { return *(const vec4_t*)(&v); }
inline const ImVec2& vec_cast(const vec2_t& v) { return *(const ImVec2*)(&v); }
inline const vec2_t& vec_cast(const ImVec2& v) { return *(const vec2_t*)(&v); }

inline ImVec4& vec_cast(vec4_t& v) { return *(ImVec4*)(&v); }
inline vec4_t& vec_cast(ImVec4& v) { return *(vec4_t*)(&v); }
inline ImVec2& vec_cast(vec2_t& v) { return *(ImVec2*)(&v); }
inline vec2_t& vec_cast(ImVec2& v) { return *(vec2_t*)(&v); }

static inline bool operator == (const ImVec2& lhs, const ImVec2& rhs) { return lhs.x == rhs.x && lhs.y == rhs.y; }
static inline bool operator != (const ImVec2& lhs, const ImVec2& rhs) { return !(lhs == rhs); }

static inline bool operator == (const ImVec4& lhs, const ImVec4& rhs) { return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w; }
static inline bool operator != (const ImVec4& lhs, const ImVec4& rhs) { return !(lhs == rhs); }

static inline float rad_to_deg(float x) {
    return x * 180.0f / 3.1415926535f;
}

static inline float deg_to_rad(float x) {
    return x * 3.1415926535f / 180.0f;
}

enum class PlaybackMode { Stopped, Playing };
enum class InterpolationMode { Nearest, Linear, Cubic };
enum class SelectionLevel { Atom, Residue, Chain };
enum class SelectionOperator { Or, And, Replace, Clear };
enum class SelectionGrowth { CovalentBond, Radial };
enum class RepresentationType { SpaceFill, Licorice, Ribbons, Cartoon };
enum class TrackingMode { Absolute, Relative };
enum class CameraMode { Perspective, Orthographic };

#define GL_COLOR_ATTACHMENT_COLOR        GL_COLOR_ATTACHMENT0
#define GL_COLOR_ATTACHMENT_NORMAL       GL_COLOR_ATTACHMENT1
#define GL_COLOR_ATTACHMENT_VELOCITY     GL_COLOR_ATTACHMENT2
#define GL_COLOR_ATTACHMENT_POST_TONEMAP GL_COLOR_ATTACHMENT3
#define GL_COLOR_ATTACHMENT_PICKING      GL_COLOR_ATTACHMENT4

enum AtomBit_ {
    AtomBit_Highlighted = 0x1,
    AtomBit_Selected    = 0x2,
    AtomBit_Visible     = 0x4
};

enum MolBit_ {
    MolBit_DirtyPosition            = 0x01,
    MolBit_DirtyRadius              = 0x02,
    MolBit_DirtySecondaryStructure  = 0x04,
    MolBit_DirtyFlags               = 0x08,
    MolBit_DirtyBonds               = 0x10,
};

enum RepBit_ {
    RepBit_DirtyColor   = 0x1,
    RepBit_DirtyFilter  = 0x2
};

// #struct Structure Declarations

struct SingleSelectionSequence {
    int32_t idx[4] = {-1, -1, -1, -1};
};

struct PickingData {
    uint32_t idx = INVALID_PICKING_IDX;
    float depth = 1.0f;
    vec3_t world_coord = {0, 0, 0};
    vec2_t screen_coord = {0, 0};
};

struct GBuffer {
    struct {
        GLuint depth = 0;
        GLuint color = 0;
        GLuint normal = 0;
        GLuint velocity = 0;
        GLuint post_tonemap = 0;
        GLuint picking = 0;
        GLuint fbo = 0;
    } deferred;

    struct {
        // @NOTE: Two of each for ping-pong read / write (hopefully non blocking reads)
        // This means that we read with one frames latency
        GLuint color[2] = {0, 0};
        GLuint depth[2] = {0, 0};
        uint32_t frame = 0;
    } pbo_picking;

    int width = 0;
    int height = 0;
};

struct Representation {
    StrBuf<32> name = "rep";
    StrBuf<256> filt = "all";
    //StrBuf<256> prop = "";
    StrBuf<256> filt_error = "";
    //StrBuf<256> prop_error = "";

    RepresentationType type = RepresentationType::SpaceFill;
    ColorMapping color_mapping = ColorMapping::Cpk;
    md_exp_bitfield_t atom_mask{};
    md_gl_representation_t md_rep{};

    bool enabled = true;
    bool show_in_selection = true;
    bool filt_is_valid = false;
    bool filt_is_dynamic = false;
    bool dynamic_evaluation = false;
    //bool prop_is_valid = false;

    uint32_t flags = 0;

    // User defined color used in uniform mode
    vec4_t uniform_color = vec4_t{1,1,1,1};

    // For colormapping the property
    ImPlotColormap color_map = ImPlotColormap_Plasma;
    float map_beg = 0;
    float map_end = 1;
    float map_min = 0;
    float map_max = 1;

    // VDW and Ball & Stick
    float radius = 1.f;

    md_script_property_t* prop = NULL;

    // Ball & Stick and Licorice, Ribbons, Cartoon
    float thickness = 1.f;

    // Ribbons, Cartoon
    float tension = 0.5f;
    float width = 1.f;
};

struct Selection {
    StrBuf<64> name = "sel";
    md_exp_bitfield_t atom_mask{};
};

// This is viamd's representation of a property
struct DisplayProperty {
    enum ShowIn {
        ShowIn_Timeline = 1,
        ShowIn_Distribution = 2,
        ShowIn_Volume = 4
    };

    struct Range {
        float beg = 0;
        float end = 0;
    };

    struct Histogram {
        float bin[1024] = {};
        Range value_range = {};
    };

    StrBuf<32> lbl = "";
    uint32_t col = 0xFFFFFFFFU;

    md_script_property_t* full_prop = NULL;
    md_script_property_t* filt_prop = NULL;

    uint64_t full_prop_fingerprint = 0;
    uint64_t filt_prop_fingerprint = 0;

    uint32_t possible_display_mask = 0;
    uint32_t current_display_mask = 0;

    Range value_filter = {};
    Range frame_filter = {};

    Histogram full_hist = {};
    Histogram filt_hist = {};
};

struct AtomElementMapping {
    StrBuf<31> lbl = "";
    md_element_t elem = 0;
};

struct ApplicationData {
    // --- APPLICATION ---
    application::Context ctx {};

    // --- FILES ---
    // for keeping track of open files
    struct {
        StrBuf<512> molecule{};
        StrBuf<512> trajectory{};
        StrBuf<512> workspace{};
    } files;

    // --- CAMERA ---
    struct {
        Camera camera{};
        TrackballControllerParam trackball_param;
        ViewParam param{};
        CameraMode mode = CameraMode::Perspective;

        struct {
            vec2_t sequence[16] {};
        } jitter;

        struct {
            vec3_t target_position{};
            quat_t target_orientation{};
        } animation;
    } view;

    // --- MOLD DATA ---
    struct {
        md_allocator_i*     mol_alloc = NULL;

        md_gl_context_t     gl_ctx = {};
        md_gl_molecule_t    gl_mol = {};
        md_molecule_t       mol = {};
        md_trajectory_i     traj = {};

        struct {
            md_script_ir_t    ir = {};
            md_script_eval_t  full_eval = {};
            md_script_eval_t  filt_eval = {};

            // Semaphore to control access to IR
            md_semaphore_t ir_semaphore = {};

            md_exp_bitfield_t frame_mask = {};

            bool ir_is_valid = false;
            bool compile_ir = false;
            bool eval_init = false;
            bool evaluate_full = false;
            bool evaluate_filt = false;
            double time_since_last_change = 0.0;
        } script;
        uint32_t dirty_buffers = {0};
    } mold;

    DisplayProperty* display_properties = 0;

    // --- ASYNC TASKS HANDLES ---
    struct {
        task_system::ID backbone_computations = task_system::INVALID_ID;
        task_system::ID prefetch_frames = task_system::INVALID_ID;
        task_system::ID evaluate_full = task_system::INVALID_ID;
        task_system::ID evaluate_filt = task_system::INVALID_ID;
        task_system::ID evaluate_shape_space = task_system::INVALID_ID;
    } tasks;

    // --- ATOM SELECTION ---
    struct {
        SelectionLevel granularity = SelectionLevel::Atom;

        int32_t hovered = -1;
        int32_t right_clicked = -1;
        SingleSelectionSequence single_selection_sequence;

        md_exp_bitfield_t current_selection_mask{};
        md_exp_bitfield_t current_highlight_mask{};
        Selection* stored_selections = NULL;

        struct {
            struct {
                vec4_t fill_color = {1.0f, 1.0f, 1.0f, 0.5f};
                vec4_t outline_color = {1.0f, 0.5f, 0.0f, 0.5f};
                float outline_scale = 1.1f;
            } highlight;

            struct {
                vec4_t fill_color = {1.0f, 1.0f, 1.0f, 0.5f};
                vec4_t outline_color = {0.0f, 0.5f, 1.0f, 0.5f};
                float outline_scale = 1.2f;
            } selection;

            float selection_saturation = 0.3f;
        } color;

        bool selecting = false;

        struct {
            char buf[256] = "";
            char err_buf[256] = "";
            md_exp_bitfield_t mask = {0};
            bool query_ok = false;
            bool query_invalid = true;
            bool show_window = false;
        } query;

        struct {
            md_exp_bitfield_t mask = {0};
            SelectionGrowth mode = SelectionGrowth::CovalentBond;
            float extent = 0;
            bool mask_invalid = true;
            bool show_window = false;
        } grow;
    } selection;

    // --- FRAMEBUFFER ---
    GBuffer gbuffer {};

    PickingData picking {};

    // --- ANIMATION ---
    struct {
        double frame = 0.f;  // double precision for long trajectories
        float fps = 10.f;
        InterpolationMode interpolation = InterpolationMode::Cubic;
        PlaybackMode mode = PlaybackMode::Stopped;
        bool apply_pbc = false;
    } animation;

    // --- TIMELINE---
    struct {
        struct {
            bool enabled = true;
            double beg_frame = 0;
            double end_frame = 1;
            
            struct {
                bool enabled = false;
                double extent_in_frames = 10;
            } temporal_window;
        } filter;

        struct {
            double beg_x = 0;
            double end_x = 1;
        } view_range;

        // Holds the timestamps for each frame
        float* x_values = 0;

        bool show_window = false;
    } timeline;

    // --- DISTRIBUTIONS ---
    struct {
        struct {
            bool enabled = false;
        } filter;
        bool show_window = false;
    } distributions;

    struct {
        bool show_window = false;
        bool enabled = false;

        struct {
            bool enabled = true;
            float density_scale = 1.f;
            struct {
                GLuint id = 0;
                bool dirty = true;
                int width = 0;
                float alpha_scale = 1.f;
                StrBuf<512> path = VIAMD_IMAGE_DIR "/tf/tf.png";
            } tf;
        } dvr;

        struct {
            bool enabled = false;
            IsoSurfaces isosurfaces;
        } iso;

        struct {
            GLuint id = 0;
            bool dirty = false;
            int dim_x = 0;
            int dim_y = 0;
            int dim_z = 0;
            float max_value = 1.f;
        } volume_texture;

        GBuffer fbo = {0};

        struct {
            vec3_t min = {0, 0, 0};
            vec3_t max = {1, 1, 1};
        } clip_volume;

        // mat4 model_to_world_matrix{};
        // mat4 world_to_model_matrix{};
        vec3_t voxel_spacing{1.0f};
        float resolution_scale = 2.0f;

        vec4_t clip_volume_color = {1,0,0,1};
        vec4_t bounding_box_color = {0,0,0,1};

        bool show_bounding_box = true;
        bool show_reference_structures = true;
        bool show_reference_ensemble = false;
        bool show_target_atoms = false;

        md_gl_representation_t* gl_reps = 0;
        mat4_t* rep_model_mats = 0;
        mat4_t model_mat = {0};

        Camera camera = {};

    } density_volume;

    // --- VISUALS ---
    struct {
        struct {
            vec3_t color = {1, 1, 1};
            float intensity = 24.f;
        } background;

        struct { 
            bool enabled = true;
            float intensity = 6.0f;
            float radius = 6.0f;
            float bias = 0.1f;
        } ssao;

#if EXPERIMENTAL_CONE_TRACED_AO == 1
        struct {
            bool enabled = true;
            float intensity = 1.0f;
            float step_scale = 1.0f;
        } cone_traced_ao;
#endif

        struct {
            bool enabled = false;
            struct {
                float target = 5.0f;
                float current = 5.0f;
            } focus_depth;
            float focus_scale = 10.0f;
        } dof;

        struct {
            bool enabled = true;
            bool jitter = true;
            float feedback_min = 0.88f;
            float feedback_max = 0.97f;

            struct {
                bool enabled = true;
                float motion_scale = 1.0f;
            } motion_blur;
        } temporal_reprojection;

        struct {
            bool enabled = true;
            postprocessing::Tonemapping tonemapper = postprocessing::Tonemapping_Filmic;
            float exposure = 1.f;
            float gamma = 2.2;
        } tonemapping;

        struct {
            bool draw_control_points = false;
            bool draw_spline = false;
        } spline;
    } visuals;

    /*
    struct {
        bool enabled = false;
        bool dirty = true;
        bool overlay = false;
        vec4_t color = vec4_t(1, 0, 1, 1);
        float distance_cutoff = HYDROGEN_BOND_DISTANCE_CUTOFF_DEFAULT;  // In Ångström
        float angle_cutoff = HYDROGEN_BOND_ANGLE_CUTOFF_DEFAULT;        // In Degrees
        DynamicArray<HydrogenBond> bonds{};
    } hydrogen_bonds;
    */

    struct {
        bool enabled = false;
        vec4_t color = {0, 0, 0, 0.5f};
        mat3_t box = {0};
    } simulation_box;

#if EXPERIMENTAL_CONE_TRACED_AO == 1
    struct {
        cone_trace::GPUVolume vol;
        std::mutex vol_mutex{};
    } occupancy_volume;
#endif
    /*
    // --- RAMACHANDRAN ---
    struct {
        bool show_window = false;

        ramachandran::ColorMap color_map{};

        struct {
            bool enabled = false;
            float radius = 0.2f;
            vec4_t color = vec4_t(0, 0, 0, 1);
        } range;

        struct {
            bool enabled = true;
            vec4_t border_color = vec4_t(0, 0, 0, 1);
            struct {
                float radius = 2.5f;
                vec4_t fill_color = vec4_t(0.75f, 0.75f, 0.75f, 1);
            } base;
            struct {
                float radius = 3.5f;
                vec4_t selection_color = vec4_t(0.8f, 1.0f, 0.5f, 1.0f);
                vec4_t highlight_color = vec4_t(0.8f, 1.0f, 0.5f, 0.5f);
            } selection;
        } current;
    } ramachandran;
    */

    struct {
        StrBuf<256> input = "all";
        StrBuf<256> err_text = "";
        md_filter_result_t result = {0};
        bool evaluate    = true;
        bool show_window = false;
        bool input_valid = false;

        // coords and weights should be of size num_frames * num_structures
        int32_t num_frames;
        int32_t num_structures;
        vec3_t* weights = 0;
        vec2_t* coords = 0;

        float marker_size = 1.4f;
    } shape_space;

    // --- REPRESENTATIONS ---
    struct {
        Representation* buffer = 0;
        md_exp_bitfield_t atom_visibility_mask = {0};
        bool atom_visibility_mask_dirty = false;
        bool show_window = false;
    } representations;

    // --- CONSOLE ---
    Console console{};

    struct {
        TextEditor editor{};
        bool show_editor = true;
    } script;

    struct {
        bool show_window = false;
        AtomElementMapping* atom_element_remappings = 0;
    } dataset;

    struct {
        struct {
            // Data is packed, no need for stride
            int64_t count = 0;
            float* data = NULL;
        } timestamp;
        struct {
            int64_t stride = 0; // = mol.backbone.count. Multiply frame idx with this to get the data
            int64_t count = 0;  // = mol.backbone.count * num_frames. Defines the end of the data for assertions
            md_secondary_structure_t* data = NULL;
        } secondary_structure;
        struct {
            int64_t stride = 0; // = mol.backbone.count. Multiply frame idx with this to get the data
            int64_t count = 0;  // = mol.backbone.count * num_frames. Defines the end of the data for assertions
            md_backbone_angles_t* data = NULL;
        } backbone_angles;
    } trajectory_data;

    bool show_debug_window = false;
    bool show_property_export_window = false;
};

//static void postprocess_frame(md_frame_data_t* frame, void* user_data);

static void compute_histogram(float* bins, int num_bins, float min_bin_val, float max_bin_val, const float* values, int num_values) {
    memset(bins, 0, sizeof(float) * num_bins);

    const float bin_range = max_bin_val - min_bin_val;
    const float inv_range = 1.0f / bin_range;
    const float scl = 1.0f / num_values;
    for (int i = 0; i < num_values; ++i) {
        if (values[i] < min_bin_val || max_bin_val < values[i]) continue;
        int idx = CLAMP((int)(((values[i] - min_bin_val) * inv_range) * num_bins), 0, num_bins - 1);
        bins[idx] += scl;
    }
}

static void downsample_histogram(float* dst_bins, int num_dst_bins, const float* src_bins, int num_src_bins) {
    ASSERT(num_dst_bins <= num_src_bins);

    memset(dst_bins, 0, sizeof(float) * num_dst_bins);
    const int factor = MAX(1, num_src_bins / num_dst_bins);
    //const float scl = 1.0f / factor;
    for (int j = 0; j < num_src_bins; ++j) {
        dst_bins[j / factor] += src_bins[j];
    }
}

static double frame_to_time(double frame, const ApplicationData& data) {
    int64_t num_frames = md_array_size(data.timeline.x_values);
    ASSERT(num_frames);
    int64_t f0 = CLAMP((int64_t)frame, 0, num_frames - 1);
    int64_t f1 = CLAMP((int64_t)frame + 1, 0, num_frames - 1);
    return lerp(data.timeline.x_values[f0], data.timeline.x_values[f1], fract(frame));
}

static double time_to_frame(double time, const ApplicationData& data) {
    int64_t num_frames = md_array_size(data.timeline.x_values);
    ASSERT(num_frames);
    // Try to map time t back into frame

    double beg = data.timeline.x_values[0];
    double end = data.timeline.x_values[num_frames - 1];
    time = CLAMP(time, beg, end);

    // Estimate the frame
    double frame_est = ((time - beg) / (end-beg)) * (num_frames - 1);
    frame_est = CLAMP(frame_est, 0, num_frames - 1);

    int64_t prev_frame_idx = CLAMP((int64_t)frame_est,     0, num_frames - 1);
    int64_t next_frame_idx = CLAMP((int64_t)frame_est + 1, 0, num_frames - 1);

    if (time < (double)data.timeline.x_values[prev_frame_idx]) {
        // Linear search down
        for (prev_frame_idx = prev_frame_idx - 1; prev_frame_idx >= 0; --prev_frame_idx) {
            next_frame_idx = prev_frame_idx + 1;
            if ((double)data.timeline.x_values[prev_frame_idx] <= time && time <= (double)data.timeline.x_values[next_frame_idx])
                break;
        }
    }
    else if (time > (double)data.timeline.x_values[next_frame_idx]) {
        // Linear search up
        for (next_frame_idx = next_frame_idx + 1; next_frame_idx < num_frames; ++next_frame_idx) {
            prev_frame_idx = next_frame_idx - 1;
            if ((double)data.timeline.x_values[prev_frame_idx] <= time && time <= (double)data.timeline.x_values[next_frame_idx])
                break;
        }
    }

    // Compute true fraction between timestamps
    double t = (time - (double)data.timeline.x_values[prev_frame_idx]) / ((double)data.timeline.x_values[next_frame_idx] - (double)data.timeline.x_values[prev_frame_idx]);
    t = CLAMP(t, 0.0, 1.0);

    // Compose frame value (base + fraction)
    return (double)prev_frame_idx + t;
}

static void single_selection_sequence_clear(SingleSelectionSequence* seq) {
    ASSERT(seq);
    for (int64_t i = 0; i < ARRAY_SIZE(seq->idx); ++i) {
        seq->idx[i] = -1;
    }
}

static void single_selection_sequence_push_back(SingleSelectionSequence* seq, int32_t idx) {
    ASSERT(seq);
    for (int64_t i = 0; i < ARRAY_SIZE(seq->idx); ++i) {
        if (seq->idx[i] == -1) {
            seq->idx[i] = idx;
            break;
        }
    }
}

static void single_selection_sequence_pop_back(SingleSelectionSequence* seq) {
    ASSERT(seq);
    int64_t i = 0;
    for (; i < ARRAY_SIZE(seq->idx); ++i) {
        if (seq->idx[i] == -1) break;
    }
    if (i > 0) {
        seq->idx[i-1] = -1;
    }
}

static int32_t single_selection_sequence_last(const SingleSelectionSequence* seq) {
    ASSERT(seq);
    int64_t i = 0;
    for (; i < ARRAY_SIZE(seq->idx); ++i) {
        if (seq->idx[i] == -1) break;
    }
    if (i > 0) {
        return seq->idx[i-1];
    }
    return -1;
}

static int64_t single_selection_sequence_count(const SingleSelectionSequence* seq) {
    int64_t i = 0;
    for (; i < ARRAY_SIZE(seq->idx); ++i) {
        if (seq->idx[i] == -1) break;
    }
    return i;
}

static void launch_prefetch_job(ApplicationData* data);
static void clear_properties(DisplayProperty** prop_items);
static void init_display_properties(DisplayProperty** prop_items, md_script_property_t* full_props, md_script_property_t* filt_props, int64_t num_props);
static void update_display_properties(ApplicationData* data);

static void interpolate_atomic_properties(ApplicationData* data);
static void update_view_param(ApplicationData* data);
static void reset_view(ApplicationData* data, bool move_camera = false, bool smooth_transition = false);
static float compute_avg_ms(float dt);

static void compute_aabb(vec3_t* aabb_min, vec3_t* aabb_max, const float* x, const float* y, const float* z, int64_t count);

static bool handle_selection(ApplicationData* data);
// static void handle_animation(ApplicationData* data);
static void handle_camera_interaction(ApplicationData* data);
static void handle_camera_animation(ApplicationData* data);

//static void init_display_properties(ApplicationData* data);
//static void update_density_volume_texture(ApplicationData* data);
static void handle_picking(ApplicationData* data);

static void fill_gbuffer(ApplicationData* data);
static void apply_postprocessing(const ApplicationData& data);

#if EXPERIMENTAL_CONE_TRACED_AO
static void init_occupancy_volume(ApplicationData* data);
#endif

static void draw_representations(ApplicationData* data);
static void draw_representations_lean_and_mean(ApplicationData* data, uint32_t mask = 0xFFFFFFFFU);

static void draw_main_menu(ApplicationData* data);
static void draw_context_popup(ApplicationData* data);
static void draw_selection_query_window(ApplicationData* data);
static void draw_selection_grow_window(ApplicationData* data);
static void draw_animation_control_window(ApplicationData* data);
static void draw_representations_window(ApplicationData* data);
//static void draw_property_window(ApplicationData* data);
static void draw_timeline_window(ApplicationData* data);
static void draw_distribution_window(ApplicationData* data);
//static void draw_ramachandran_window(ApplicationData* data);
static void draw_atom_info_window(const ApplicationData& data, int atom_idx);
static void draw_molecule_dynamic_info_window(ApplicationData* data);
static void draw_async_task_window(ApplicationData* data);
//static void draw_reference_frame_window(ApplicationData* data);
static void draw_shape_space_window(ApplicationData* data);
static void draw_density_volume_window(ApplicationData* data);
static void draw_script_editor_window(ApplicationData* data);
static void draw_dataset_window(ApplicationData* data);
static void draw_debug_window(ApplicationData* data);
static void draw_property_export_window(ApplicationData* data);
// static void draw_density_volume_clip_plane_widgets(ApplicationData* data);
// static void draw_selection_window(ApplicationData* data);

static void clear_gbuffer(GBuffer* gbuf);
static void init_gbuffer(GBuffer* gbuf, int width, int height);
static void destroy_gbuffer(GBuffer* gbuf);
static PickingData read_picking_data(GBuffer* fbo, int32_t x, int32_t y);

static void update_md_buffers(ApplicationData* data);

static void init_molecule_data(ApplicationData* data);
static void init_trajectory_data(ApplicationData* data);

static bool load_dataset_from_file(ApplicationData* data, str_t file);
//static void free_dataset(ApplicationData* data);

static void load_workspace(ApplicationData* data, str_t file);
static void save_workspace(ApplicationData* data, str_t file);

static void create_screenshot(ApplicationData* data);

// Representations
static Representation* create_representation(ApplicationData* data, RepresentationType type = RepresentationType::SpaceFill,
                                             ColorMapping color_mapping = ColorMapping::Cpk, str_t filter = make_cstr("all"));
static Representation* clone_representation(ApplicationData* data, const Representation& rep);
static void remove_representation(ApplicationData* data, int idx);
static void update_representation(ApplicationData* data, Representation* rep);
static void update_all_representations(ApplicationData* data);
static void init_representation(ApplicationData* data, Representation* rep);
static void init_all_representations(ApplicationData* data);
static void clear_representations(ApplicationData* data);

static void recompute_atom_visibility_mask(ApplicationData* data);

// Selections
static Selection* create_selection(ApplicationData* data, str_t name, md_exp_bitfield_t* bf);
static Selection* clone_selection(ApplicationData* data, const Selection& sel);
static void remove_selection(ApplicationData* data, int idx);

static void reset_selections(ApplicationData* data);
static void clear_selections(ApplicationData* data);

static bool filter_expression(ApplicationData* data, str_t expr, md_exp_bitfield_t* mask, bool* is_dynamic, char* error_str, int64_t error_cap);

static void modify_field(md_exp_bitfield_t* bf, const md_exp_bitfield_t* mask, SelectionOperator op) {
    switch(op) {
    case SelectionOperator::Or:
        md_bitfield_or_inplace(bf, mask);
        break;
    case SelectionOperator::And:
        md_bitfield_and_inplace(bf, mask);
        break;
    case SelectionOperator::Replace:
        md_bitfield_copy(bf, mask);
        break;
    default:
        ASSERT(false);
    }
}

static void modify_field(md_exp_bitfield_t* bf, md_range_t range, SelectionOperator op) {
    switch(op) {
    case SelectionOperator::Or:
        md_bitfield_set_range(bf, range.beg, range.end);
        break;
    case SelectionOperator::And:
        md_bitfield_clear_range(bf, 0, range.beg);
        md_bitfield_clear_range(bf, range.end, bf->end_bit);
        break;
    case SelectionOperator::Replace:
        md_bitfield_clear(bf);
        md_bitfield_set_range(bf, range.beg, range.end);
        break;
    case SelectionOperator::Clear:
        md_bitfield_clear_range(bf, range.beg, range.end);
        break;
    default:
        ASSERT(false);
    }
}

static void clear_highlight(ApplicationData* data) {
    ASSERT(data);
    md_bitfield_clear(&data->selection.current_highlight_mask);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

static void modify_highlight(ApplicationData* data, md_exp_bitfield_t* mask, SelectionOperator op = SelectionOperator::Replace) {
    ASSERT(data);
    modify_field(&data->selection.current_highlight_mask, mask, op);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

static void modify_highlight(ApplicationData* data, md_range_t range, SelectionOperator op = SelectionOperator::Replace) {
    ASSERT(data);
    modify_field(&data->selection.current_highlight_mask, range, op);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

static void clear_selection(ApplicationData* data) {
    ASSERT(data);
    md_bitfield_clear(&data->selection.current_selection_mask);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

static void modify_selection(ApplicationData* data, md_exp_bitfield_t* atom_mask, SelectionOperator op = SelectionOperator::Replace) {
    ASSERT(data);
    modify_field(&data->selection.current_selection_mask, atom_mask, op);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

static void modify_selection(ApplicationData* data, md_range_t range, SelectionOperator op = SelectionOperator::Replace) {
    ASSERT(data);
    modify_field(&data->selection.current_selection_mask, range, op);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

static void on_trajectory_load_complete(ApplicationData* data);

// Global data for application
static md_allocator_i* frame_allocator = 0;
#if MD_DEBUG
static md_allocator_i* persistent_allocator = md_tracking_allocator_create(default_allocator);
#elif MD_RELEASE
static md_allocator_i* persistent_allocator = default_allocator;
#else
    
#endif

int main(int, char**) {
    const int64_t stack_size = MEGABYTES(64);
    void* stack_mem = md_alloc(default_allocator, stack_size);
    md_stack_allocator_t stack_alloc {};
    md_stack_allocator_init(&stack_alloc, stack_mem, stack_size);
    md_allocator_i stack_interface = md_stack_allocator_create_interface(&stack_alloc);
    frame_allocator = &stack_interface;

    ApplicationData data;

    md_bitfield_init(&data.selection.current_selection_mask, persistent_allocator);
    md_bitfield_init(&data.selection.current_highlight_mask, persistent_allocator);
    md_bitfield_init(&data.selection.query.mask, persistent_allocator);
    md_bitfield_init(&data.selection.grow.mask, persistent_allocator);

    md_bitfield_init(&data.representations.atom_visibility_mask, persistent_allocator);
    md_bitfield_init(&data.mold.script.frame_mask, persistent_allocator);

    md_semaphore_init(&data.mold.script.ir_semaphore, IR_SEMAPHORE_MAX_COUNT);

    md_logger_i logger = {
        .inst = (md_logger_o*)&data.console,
        .log = [](struct md_logger_o* inst, enum md_log_type log_type, const char* msg) {
            const char* modifier = "";
            switch (log_type) {
            case MD_LOG_TYPE_DEBUG:
                modifier = "[debug] ";
                break;
            case MD_LOG_TYPE_INFO:
                modifier = "[info] ";
                break;
            case MD_LOG_TYPE_ERROR:
                modifier = "[error] ";
                break;
            default:
                break;
            }
            ((Console*)inst)->AddLog("%s%s", modifier, msg);
        }
    };

      md_add_logger(&logger);

    // Init platform
    md_print(MD_LOG_TYPE_INFO, "Initializing GL...");
    if (!application::initialize(&data.ctx, 1920, 1080, "VIAMD")) {
        md_print(MD_LOG_TYPE_ERROR, "Could not initialize platform layer... terminating\n");
        return -1;
    }
    data.ctx.window.vsync = true;

    md_print(MD_LOG_TYPE_INFO, "Creating framebuffer...");
    init_gbuffer(&data.gbuffer, data.ctx.framebuffer.width, data.ctx.framebuffer.height);

    generate_halton_sequence(data.view.jitter.sequence, ARRAY_SIZE(data.view.jitter.sequence), 2, 3);

    // Init subsystems
    md_print(MD_LOG_TYPE_INFO, "Initializing immediate draw...");
    immediate::initialize();
    //md_print(MD_LOG_TYPE_INFO, "Initializing ramachandran...");
    //ramachandran::initialize();
    md_print(MD_LOG_TYPE_INFO, "Initializing post processing...");
    postprocessing::initialize(data.gbuffer.width, data.gbuffer.height);
    md_print(MD_LOG_TYPE_INFO, "Initializing volume...");
    volume::initialize();
#if EXPERIMENTAL_CONE_TRACED_AO == 1
    md_print(MD_LOG_TYPE_INFO, "Initializing cone tracing...");
    cone_trace::initialize();
#endif
    md_print(MD_LOG_TYPE_INFO, "Initializing task system...");
    task_system::initialize();

    md_gl_context_init(&data.mold.gl_ctx);

    ImGui::init_theme();

    const ImU32 dihedral_colors[] = {
        IM_COL32(255,  0,255,255),
        IM_COL32(0,    0,255,255),
        IM_COL32(255,255,255,255),
        IM_COL32(255,  0,  0,255),
        IM_COL32(255,  0,255,255)
    };

    ImPlot::AddColormap("Dihedral", dihedral_colors, ARRAY_SIZE(dihedral_colors), false);

    data.script.editor.SetLanguageDefinition(TextEditor::LanguageDefinition::VIAMD());
    data.script.editor.SetPalette(TextEditor::GetDarkPalette());

    load_dataset_from_file(&data, make_cstr(VIAMD_DATASET_DIR "/1ALA-500.pdb"));
    create_representation(&data, RepresentationType::SpaceFill, ColorMapping::Cpk, make_cstr("all"));
    data.script.editor.SetText("s1 = resname(\"ALA\")[2:8];\nd1 = distance(10,30);\na1 = angle(1,2,3) in resname(\"ALA\");\nv = sdf(s1, element('H'), 10.0);");

    reset_view(&data, true);
    recompute_atom_visibility_mask(&data);

    //init_density_volume(&data);
    interpolate_atomic_properties(&data);

#if EXPERIMENTAL_CONE_TRACED_AO == 1
    init_occupancy_volume(&data);
#endif

#if EXPERIMENTAL_SDF == 1
    draw::scan::test_scan();
#endif

    auto& mol = data.mold.mol;
    auto& traj = data.mold.traj;

    bool time_changed = true;
    bool time_stopped = true;

    // Main loop
    while (!data.ctx.window.should_close) {
        application::update(&data.ctx);
        
        // This needs to happen first (in imgui events) to enable docking of imgui windows
        ImGui::CreateDockspace();

#if SHOW_IMGUI_DEMO_WINDOW
        ImGui::ShowDemoWindow();
        ImPlot::ShowDemoWindow();
#endif

        const int64_t num_frames = md_trajectory_num_frames(&traj);
        const int64_t last_frame = MAX(0, num_frames - 1);
        const double max_frame = (double)MAX(0, last_frame);

        // #input
        if (data.ctx.input.key.hit[KEY_CONSOLE]) {
            if (data.console.Visible()) {
                data.console.Hide();
            } else if (!ImGui::GetIO().WantTextInput) {
                data.console.Show();
            }
        }

        if (data.ctx.input.key.hit[KEY_SHOW_DEBUG_WINDOW]) {
            data.show_debug_window = true;
        }

        if (!ImGui::GetIO().WantCaptureKeyboard) {
            if (data.ctx.input.key.hit[KEY_TOGGLE_SCREENSHOT_MODE]) {
                static bool screenshot_mode = false;
                screenshot_mode = !screenshot_mode;

                ImGuiStyle& style = ImGui::GetStyle();
                if (screenshot_mode) {
                    ImGui::StyleColorsClassic(&style);
                } else {
                    ImGui::StyleColorsLight(&style);
                    // style.Colors[ImGuiCol_WindowBg] = ImVec4(0.0f, 0.0f, 0.0f, 0.0f);
                    // style.Colors[ImGuiCol_ChildWindowBg] = ImVec4(0.0f, 0.0f, 0.0f, 0.0f);
                }
            }

            if (data.ctx.input.key.hit[Key::KEY_F5]) {
                md_print(MD_LOG_TYPE_INFO, "Recompiling shaders and re-initializing volume");
                postprocessing::initialize(data.gbuffer.width, data.gbuffer.height);
                volume::initialize();
#if EXPERIMENTAL_CONE_TRACED_AO == 1
                cone_trace::initialize();
#endif
                md_gl_context_free(&data.mold.gl_ctx);
                md_gl_context_init(&data.mold.gl_ctx);
            }

            if (data.ctx.input.key.hit[KEY_PLAY_PAUSE]) {
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

            if (data.ctx.input.key.hit[KEY_SKIP_TO_PREV_FRAME] || data.ctx.input.key.hit[KEY_SKIP_TO_NEXT_FRAME]) {
                double step = data.ctx.input.key.down[Key::KEY_LEFT_CONTROL] ? 10.0 : 1.0;
                if (data.ctx.input.key.hit[KEY_SKIP_TO_PREV_FRAME]) step = -step;
                data.animation.frame = CLAMP(data.animation.frame + step, 0.0, max_frame);
            }
        }

        if (data.representations.atom_visibility_mask_dirty) {
            recompute_atom_visibility_mask(&data);
            data.representations.atom_visibility_mask_dirty = false;
        }

        data.selection.selecting = false;
        if (!ImGui::GetIO().WantCaptureMouse) {
            data.selection.selecting |= handle_selection(&data);
        }
        // if (data.selection.selecting) update_all_representations(&data);

        handle_camera_interaction(&data);
        handle_camera_animation(&data);


        if (data.animation.mode == PlaybackMode::Playing) {
            data.animation.frame += data.ctx.timing.delta_s * data.animation.fps;
            data.animation.frame = CLAMP(data.animation.frame, 0.0, max_frame);
            if (data.animation.frame >= max_frame) {
                data.animation.mode = PlaybackMode::Stopped;
                data.animation.frame = max_frame;
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
            data.timeline.filter.beg_frame = CLAMP(data.animation.frame - half_window_ext, 0.0, max_frame);
            data.timeline.filter.end_frame = CLAMP(data.animation.frame + half_window_ext, 0.0, max_frame);
            if (data.mold.script.ir_is_valid && (data.timeline.filter.beg_frame != pre_beg || data.timeline.filter.end_frame != pre_end)) {
                data.mold.script.evaluate_filt = true;
            }
        }

        if (time_changed) {
            time_stopped = false;

            PUSH_CPU_SECTION("Interpolate Position")
            if (md_trajectory_num_frames(&traj)) {
                interpolate_atomic_properties(&data);

                if (data.animation.apply_pbc) {
                    md_util_apply_pbc_args_t args = {
                        .atom = {
                            .count = mol.atom.count,
                            .x = mol.atom.x,
                            .y = mol.atom.y,
                            .z = mol.atom.z,
                        },
                        .residue = {
                            .count = mol.residue.count,
                            .atom_range = mol.residue.atom_range,
                        },
                        .chain = {
                            .count = mol.chain.count,
                            .residue_range = mol.chain.residue_range,
                        },
                    };
                    memcpy(args.pbc.box, &data.simulation_box.box, sizeof(args.pbc.box));
                    md_util_apply_pbc(mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, args);
                }

#if EXPERIMENTAL_CONE_TRACED_AO
                if (data.visuals.cone_traced_ao.enabled) {
                    init_occupancy_volume(&data);
                }
#endif
            }
            POP_CPU_SECTION()

            PUSH_CPU_SECTION("Update dynamic representations")
            for (int64_t i = 0; i < md_array_size(data.representations.buffer); ++i) {
                auto& rep = data.representations.buffer[i];
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

        if (data.mold.script.compile_ir) {
            data.mold.script.time_since_last_change += data.ctx.timing.delta_s;

            if (data.mold.script.time_since_last_change > COMPILATION_TIME_DELAY_IN_SECONDS) {
                // We cannot recompile while it is evaluating.
                // Need to interrupt and wait for tasks to finish.
                md_script_eval_interrupt(&data.mold.script.full_eval);
                md_script_eval_interrupt(&data.mold.script.filt_eval);

                // Try aquire all semaphores
                if (md_semaphore_try_aquire_n(&data.mold.script.ir_semaphore, IR_SEMAPHORE_MAX_COUNT)) {
                    // Now we hold all semaphores for the script
                    data.mold.script.compile_ir = false;
                    data.mold.script.time_since_last_change = 0;

                    TextEditor& editor = data.script.editor;

                    std::string src = editor.GetText();
                    str_t src_str = {.ptr = src.data(), .len = (int64_t)src.length()};

                    editor.ClearMarkers();
                    editor.ClearErrorMarkers();

                    data.mold.script.ir_is_valid = md_script_ir_compile(&data.mold.script.ir, src_str, &data.mold.mol, persistent_allocator, NULL);

                    if (data.mold.script.ir_is_valid) {
                        data.mold.script.eval_init = true;
                        data.mold.script.evaluate_full = true;
                        data.mold.script.evaluate_filt = true;
                    } else {
                        TextEditor::ErrorMarkers markers;
                        const int64_t num_errors = data.mold.script.ir.num_errors;
                        if (num_errors) {
                            const md_script_error_t* errors = data.mold.script.ir.errors;
                            for (int64_t i = 0; i < num_errors; ++i) {
                                std::string err_str(errors[i].error.ptr, errors[i].error.len);
                                std::pair<int, std::string> pair = {errors[i].line, err_str};
                                markers.insert(pair);
                            }
                        }
                        editor.SetErrorMarkers(markers);
                    }

                    for (int64_t i = 0; i < data.mold.script.ir.num_tokens; ++i) {
                        const md_script_token_t& tok = data.mold.script.ir.tokens[i];
                        TextEditor::Marker marker = {0};
                        marker.begCol = tok.col_beg;
                        marker.endCol = tok.col_end;
                        marker.bgColor = ImVec4(1,1,1,0.5);
                        marker.depth = tok.depth;
                        marker.line = tok.line;
                        marker.onlyShowBgOnMouseOver = true;
                        marker.text = std::string(tok.text.ptr, tok.text.len);
                        marker.payload = (void*)tok.vis_token;
                        editor.AddMarker(marker);
                    }

                    md_semaphore_release_n(&data.mold.script.ir_semaphore, IR_SEMAPHORE_MAX_COUNT);
                }
            }
        }

        if (num_frames > 0) {
            if (data.mold.script.eval_init) {
                data.mold.script.eval_init = false;
                md_script_eval_init(&data.mold.script.full_eval, num_frames, &data.mold.script.ir, persistent_allocator);
                md_script_eval_init(&data.mold.script.filt_eval, num_frames, &data.mold.script.ir, persistent_allocator);
                const int64_t num_props = data.mold.script.full_eval.num_properties;
                ASSERT(data.mold.script.filt_eval.num_properties == num_props);
                init_display_properties(&data.display_properties, data.mold.script.full_eval.properties, data.mold.script.filt_eval.properties, num_props);
            }

            if (data.mold.script.evaluate_full) {
                if (data.tasks.evaluate_full.id != 0) {
                    md_script_eval_interrupt(&data.mold.script.full_eval);
                } else if (data.mold.script.ir_is_valid && md_semaphore_try_aquire(&data.mold.script.ir_semaphore)) {
                    data.mold.script.evaluate_full = false;
                    data.tasks.evaluate_full = task_system::enqueue_pool("Eval Full", 1, [&data](task_system::TaskSetRange) {
                        md_script_eval_compute(&data.mold.script.full_eval, &data.mold.script.ir, &data.mold.mol, &data.mold.traj, NULL);

                        data.tasks.evaluate_full.id = 0;
                        md_semaphore_release(&data.mold.script.ir_semaphore);
                    });
                }
            }

            if (data.timeline.filter.enabled && data.mold.script.evaluate_filt) {
                if (data.tasks.evaluate_filt.id != 0) {
                    md_script_eval_interrupt(&data.mold.script.filt_eval);
                } else if (data.mold.script.ir_is_valid && md_semaphore_try_aquire(&data.mold.script.ir_semaphore)) {
                    data.mold.script.evaluate_filt = false;
                    data.tasks.evaluate_filt = task_system::enqueue_pool("Eval Filt", 1, [&data](task_system::TaskSetRange) {
                        int64_t num_frames = md_trajectory_num_frames(&data.mold.traj);
                        int32_t beg_frame = CLAMP((int64_t)data.timeline.filter.beg_frame, 0, num_frames-1);
                        int32_t end_frame = CLAMP((int64_t)data.timeline.filter.end_frame, 0, num_frames);
                        end_frame = MAX(beg_frame + 1, end_frame);
                        md_bitfield_clear(&data.mold.script.frame_mask);
                        md_bitfield_set_range(&data.mold.script.frame_mask, beg_frame, end_frame);

                        md_script_eval_compute(&data.mold.script.filt_eval, &data.mold.script.ir, &data.mold.mol, &data.mold.traj, &data.mold.script.frame_mask);

                        data.tasks.evaluate_filt.id = 0;
                        md_semaphore_release(&data.mold.script.ir_semaphore);
                    });
                }
            }
        }

#if 0
        PUSH_CPU_SECTION("Hydrogen bonds")
        if (data.hydrogen_bonds.enabled && data.hydrogen_bonds.dirty) {
            data.hydrogen_bonds.bonds = hydrogen_bond::compute_bonds(
                {mol.hydrogen_bond.donor.data, mol.hydrogen_bond.donor.count}, {mol.hydrogen_bond.acceptor.data, mol.hydrogen_bond.acceptor.count},
                mol.atom.position, data.hydrogen_bonds.distance_cutoff, deg_to_rad(data.hydrogen_bonds.angle_cutoff));
            data.hydrogen_bonds.dirty = false;
        }
        POP_CPU_SECTION()
#endif

        // Resize Framebuffer
        if ((data.gbuffer.width != data.ctx.framebuffer.width || data.gbuffer.height != data.ctx.framebuffer.height) &&
            (data.ctx.framebuffer.width != 0 && data.ctx.framebuffer.height != 0)) {
            init_gbuffer(&data.gbuffer, data.ctx.framebuffer.width, data.ctx.framebuffer.height);
            postprocessing::initialize(data.gbuffer.width, data.gbuffer.height);
        }

        update_display_properties(&data);
        update_view_param(&data);
        update_md_buffers(&data);

        clear_gbuffer(&data.gbuffer);
        fill_gbuffer(&data);

#if EXPERIMENTAL_CONE_TRACED_AO == 1
        if (data.visuals.cone_traced_ao.enabled) {
            PUSH_GPU_SECTION("Cone-Trace AO")
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, data.gbuffer.deferred.gbuffer);
            glViewport(0, 0, data.gbuffer.width, data.gbuffer.height);
            glDrawBuffer(GL_COLOR_ATTACHMENT0);  // Modify to color buffer
            cone_trace::render_directional_occlusion(data.gbuffer.deferred.depth, data.gbuffer.deferred.normal, data.occupancy_volume.vol,
                                                     data.view.param.matrix.current.view, data.view.param.matrix.current.proj,
                                                     data.visuals.cone_traced_ao.intensity, data.visuals.cone_traced_ao.step_scale);
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
            POP_GPU_SECTION()
        }
#endif
        handle_picking(&data);

        // Activate backbuffer
        glDisable(GL_DEPTH_TEST);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glViewport(0, 0, data.ctx.framebuffer.width, data.ctx.framebuffer.height);
        glDrawBuffer(GL_BACK);
        glClear(GL_COLOR_BUFFER_BIT);

        apply_postprocessing(data);

        clear_highlight(&data);

        // GUI
        if (data.representations.show_window) draw_representations_window(&data);
        if (data.timeline.show_window) draw_timeline_window(&data);
        if (data.distributions.show_window) draw_distribution_window(&data);
        if (data.density_volume.show_window) draw_density_volume_window(&data);
        //if (data.ramachandran.show_window) draw_ramachandran_window(&data);
        if (data.shape_space.show_window) draw_shape_space_window(&data);
        if (data.script.show_editor) draw_script_editor_window(&data);
        if (data.dataset.show_window) draw_dataset_window(&data);
        if (data.selection.query.show_window) draw_selection_query_window(&data);
        if (data.selection.grow.show_window) draw_selection_grow_window(&data);
        if (data.show_property_export_window) draw_property_export_window(&data);
        if (data.show_debug_window) draw_debug_window(&data);
        // @NOTE: ImGui::GetIO().WantCaptureMouse does not work with Menu
        if (!ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow)) {
            if (data.picking.idx != INVALID_PICKING_IDX) {
                draw_atom_info_window(data, data.picking.idx);
            }
        }


        data.console.Draw("VIAMD", data.ctx.window.width, data.ctx.window.height, (float)data.ctx.timing.delta_s);
        draw_main_menu(&data);
        draw_context_popup(&data);
        draw_async_task_window(&data);
        draw_animation_control_window(&data);
        draw_molecule_dynamic_info_window(&data);

        PUSH_GPU_SECTION("Imgui render")
        application::render_imgui(&data.ctx);
        POP_GPU_SECTION()

        // Swap buffers
        application::swap_buffers(&data.ctx);

        task_system::run_main_tasks();

        // Reset frame allocator
        md_stack_allocator_reset(&stack_alloc);
    }

    // shutdown subsystems
    md_print(MD_LOG_TYPE_INFO, "Shutting down immediate draw...");
    immediate::shutdown();
    //md_print(MD_LOG_TYPE_INFO, "Shutting down ramachandran...");
    //ramachandran::shutdown();
    md_print(MD_LOG_TYPE_INFO, "Shutting down post processing...");
    postprocessing::shutdown();
    md_print(MD_LOG_TYPE_INFO, "Shutting down volume...");
    volume::shutdown();
#if EXPERIMENTAL_CONE_TRACED_AO == 1
    md_print(MD_LOG_TYPE_INFO, "Shutting down cone tracing...");
    cone_trace::shutdown();
#endif
    md_print(MD_LOG_TYPE_INFO, "Shutting down task system...");
    task_system::shutdown();

    destroy_gbuffer(&data.gbuffer);
    application::shutdown(&data.ctx);

    return 0;
}

static void init_display_properties(DisplayProperty** prop_items, md_script_property_t* full_props, md_script_property_t* filt_props, int64_t num_props) {
    ASSERT(prop_items);
    if (num_props > 0) {
        ASSERT(full_props);
        ASSERT(filt_props);
    }

    DisplayProperty* new_items = 0;
    DisplayProperty* old_items = *prop_items;

    for (int64_t i = 0; i < num_props; ++i) {
        uint32_t possible_display_mask = 0;
        switch (full_props[i].type) {
        case MD_SCRIPT_PROPERTY_TYPE_TEMPORAL:
            possible_display_mask = DisplayProperty::ShowIn_Timeline | DisplayProperty::ShowIn_Distribution;
            break;
        case MD_SCRIPT_PROPERTY_TYPE_DISTRIBUTION:
            possible_display_mask = DisplayProperty::ShowIn_Distribution;
            break;
        case MD_SCRIPT_PROPERTY_TYPE_VOLUME:
            possible_display_mask = DisplayProperty::ShowIn_Volume;
            break;
        case MD_SCRIPT_PROPERTY_TYPE_INVALID:
        default:
            ASSERT(false);
        }

        DisplayProperty item;
        item.lbl = full_props[i].ident;
        item.col = PROPERTY_COLORS[i % ARRAY_SIZE(PROPERTY_COLORS)];
        item.frame_filter = {0, 0};
        item.value_filter = {full_props[i].data.min_range[0], full_props[i].data.max_range[0]};
        item.full_prop = &full_props[i];
        item.filt_prop = &filt_props[i];
        item.full_prop_fingerprint = 0;
        item.filt_prop_fingerprint = 0;
        item.possible_display_mask = possible_display_mask;
        item.current_display_mask = possible_display_mask;

        for (int64_t j = 0; j < md_array_size(old_items); ++j) {
            if (compare_str(item.lbl, old_items[j].lbl)) {
                // Copy relevant parameters from existing item which we want to be persistent
                item.frame_filter = old_items[j].frame_filter;
                item.value_filter = old_items[j].value_filter;
                item.current_display_mask = old_items[j].current_display_mask & possible_display_mask;
                item.full_prop_fingerprint = old_items[j].full_prop_fingerprint;
                item.filt_prop_fingerprint = old_items[j].filt_prop_fingerprint;
            }
        }
        md_array_push(new_items, item, frame_allocator);
    }

    md_array_resize(*prop_items, md_array_size(new_items), persistent_allocator);
    memcpy(*prop_items, new_items, md_array_size(new_items) * sizeof(DisplayProperty));
}

static void update_display_properties(ApplicationData* data) {
    ASSERT(data);
    DisplayProperty* disp_props = data->display_properties;
    for (int64_t i = 0; i < md_array_size(disp_props); ++i) {
        if (disp_props[i].full_prop_fingerprint != disp_props[i].full_prop->data.fingerprint) {
            disp_props[i].full_prop_fingerprint = disp_props[i].full_prop->data.fingerprint;
            const md_script_property_t* p = disp_props[i].full_prop;
            if (p->type == MD_SCRIPT_PROPERTY_TYPE_TEMPORAL) {
                DisplayProperty::Histogram& hist = disp_props[i].full_hist;
                hist.value_range = {p->data.min_range[0], p->data.max_range[0]};
                compute_histogram(hist.bin, ARRAY_SIZE(hist.bin), hist.value_range.beg, hist.value_range.end, p->data.values, p->data.num_values);
            }
        }

        if (disp_props[i].filt_prop_fingerprint != disp_props[i].filt_prop->data.fingerprint) {
            disp_props[i].filt_prop_fingerprint = disp_props[i].filt_prop->data.fingerprint;
            const md_script_property_t* p = disp_props[i].full_prop;
            if (p->type == MD_SCRIPT_PROPERTY_TYPE_TEMPORAL) {
                DisplayProperty::Histogram& hist = disp_props[i].filt_hist;
                // We copy the full range here since we want the histograms to align on the x-axis
                // This also implys that we have a dependency to the full property, since we can only now the full range of values when full_prop has been completely evaluated.
                hist.value_range = disp_props[i].full_hist.value_range;
                int offset = (int)data->timeline.filter.beg_frame * p->data.dim[0];
                int length = ((int)data->timeline.filter.end_frame - (int)data->timeline.filter.beg_frame) * p->data.dim[0];
                compute_histogram(hist.bin, ARRAY_SIZE(hist.bin), hist.value_range.beg, hist.value_range.end, p->data.values + offset, length);
            }
        }
    }
}

static void interpolate_atomic_properties(ApplicationData* data) {
    ASSERT(data);
    const auto& mol = data->mold.mol;
    const auto& traj = data->mold.traj;

    if (!mol.atom.count || !md_trajectory_num_frames(&traj)) return;

    const int64_t last_frame = MAX(0LL, md_trajectory_num_frames(&traj) - 1);
    // This is not actually time, but the fractional frame representation
    const double time = CLAMP(data->animation.frame, 0.0, double(last_frame));

    const float t = (float)fractf(time);
    const int64_t frame = (int64_t)time;
    const int64_t nearest_frame = CLAMP((int64_t)(time + 0.5), 0LL, last_frame);

    const int64_t frames[4] = {
        MAX(0LL, frame - 1),
        MAX(0LL, frame),
        MIN(frame + 1, last_frame),
        MIN(frame + 2, last_frame)
    };

    mat3_t boxes[4] = {};

    int64_t stride = ROUND_UP(mol.atom.count, md_simd_widthf);    // The interploation uses SIMD vectorization without bounds, so we make sure there is no overlap between the data segments
    int64_t bytes = stride * sizeof(float) * 3 * 4;
    void* mem = md_alloc(frame_allocator, bytes);
    defer { md_free(frame_allocator, mem, bytes); };

    float* x[4] = {
        (float*)mem + stride * 0,
        (float*)mem + stride * 1,
        (float*)mem + stride * 2,
        (float*)mem + stride * 3,
    };
    float* y[4] = {
        (float*)mem + stride * 4,
        (float*)mem + stride * 5,
        (float*)mem + stride * 6,
        (float*)mem + stride * 7,
    };
    float* z[4] = {
        (float*)mem + stride * 8,
        (float*)mem + stride * 9,
        (float*)mem + stride * 10,
        (float*)mem + stride * 11,
    };

    const InterpolationMode mode = (frames[1] != frames[2]) ? data->animation.interpolation : InterpolationMode::Nearest;

    mat3_t box = {};
    switch (mode) {
        case InterpolationMode::Nearest:
        {
            md_trajectory_frame_header_t header = {0};
            md_trajectory_load_frame(&data->mold.traj, nearest_frame, &header, mol.atom.x, mol.atom.y, mol.atom.z);
            memcpy(&box, header.box, sizeof(box));
            break;
        }
        case InterpolationMode::Linear:
        {
            md_trajectory_frame_header_t header[2] = {0};

            md_trajectory_load_frame(&data->mold.traj, frames[1], &header[0], x[0], y[0], z[0]);
            md_trajectory_load_frame(&data->mold.traj, frames[2], &header[1], x[1], y[1], z[1]);

            memcpy(&boxes[0], header[0].box, sizeof(boxes[0]));
            memcpy(&boxes[1], header[1].box, sizeof(boxes[1]));

            box = lerp(boxes[0], boxes[1], t);

            md_util_linear_interpolation_args_t args = {
                .coord = {
                    .count = mol.atom.count,
                    .dst = {
                    .x = mol.atom.x,
                    .y = mol.atom.y,
                    .z = mol.atom.z,
                },
                .src = {
                    {
                        .x = x[0],
                        .y = y[0],
                        .z = z[0],
                    },
                    {
                        .x = x[1],
                        .y = y[1],
                        .z = z[1],
                    }
                },
                },
                .pbc = {0}, // memcpy this afterwards
                .t = t
            };
            memcpy(args.pbc.box, &box, sizeof(args.pbc.box));
            md_util_linear_interpolation(args);
        }
            break;
        case InterpolationMode::Cubic:
        {
            md_trajectory_frame_header_t header[4] = {0};

            md_trajectory_load_frame(&data->mold.traj, frames[0], &header[0], x[0], y[0], z[0]);
            md_trajectory_load_frame(&data->mold.traj, frames[1], &header[1], x[1], y[1], z[1]);
            md_trajectory_load_frame(&data->mold.traj, frames[2], &header[2], x[2], y[2], z[2]);
            md_trajectory_load_frame(&data->mold.traj, frames[3], &header[3], x[3], y[3], z[3]);

            memcpy(&boxes[0], header[0].box, sizeof(boxes[0]));
            memcpy(&boxes[1], header[1].box, sizeof(boxes[1]));
            memcpy(&boxes[2], header[2].box, sizeof(boxes[2]));
            memcpy(&boxes[3], header[3].box, sizeof(boxes[3]));

            box = cubic_spline(boxes[0], boxes[1], boxes[2], boxes[3], t);

            md_util_cubic_interpolation_args_t args = {
                .coord = {
                    .count = mol.atom.count,
                    .dst = {
                    .x = mol.atom.x,
                    .y = mol.atom.y,
                    .z = mol.atom.z,
                },
                .src = {
                    {
                        .x = x[0],
                        .y = y[0],
                        .z = z[0],
                    },
                    {
                        .x = x[1],
                        .y = y[1],
                        .z = z[1],
                    },
                    {
                        .x = x[2],
                        .y = y[2],
                        .z = z[2],
                    },
                    {
                        .x = x[3],
                        .y = y[3],
                        .z = z[3],
                    },
                },
                },
                .pbc = {0}, // memcpy this afterwards
                .t = t,
                .tension = 0.5f,
            };
            memcpy(args.pbc.box, &box, sizeof(args.pbc.box));
            md_util_cubic_interpolation(args);
        }
            break;
        default:
            ASSERT(false);
    }

    if (mol.backbone.angle) {
        const md_backbone_angles_t* src_angles[4] = {
            data->trajectory_data.backbone_angles.data + data->trajectory_data.backbone_angles.stride * frames[0],
            data->trajectory_data.backbone_angles.data + data->trajectory_data.backbone_angles.stride * frames[1],
            data->trajectory_data.backbone_angles.data + data->trajectory_data.backbone_angles.stride * frames[2],
            data->trajectory_data.backbone_angles.data + data->trajectory_data.backbone_angles.stride * frames[3]
        };

        switch (mode) {
        case InterpolationMode::Nearest: {
            const md_backbone_angles_t* src_angle = t < 0.5f ? src_angles[1] : src_angles[2];
            memcpy(mol.backbone.angle, src_angle, mol.backbone.count * sizeof(md_backbone_angles_t));
            break;
        }
        case InterpolationMode::Linear: {
            for (int64_t i = 0; i < mol.backbone.count; ++i) {
                float phi = lerp(src_angles[1][i].phi, src_angles[2][i].phi, t);
                float psi = lerp(src_angles[1][i].psi, src_angles[2][i].psi, t);
                mol.backbone.angle[i] = {phi, psi};
            }
            break;
        }
        case InterpolationMode::Cubic: {
            for (int64_t i = 0; i < mol.backbone.count; ++i) {
                float phi = cubic_spline(src_angles[0][i].phi, src_angles[1][i].phi, src_angles[2][i].phi, src_angles[3][i].phi, t);
                float psi = cubic_spline(src_angles[0][i].psi, src_angles[1][i].psi, src_angles[2][i].psi, src_angles[3][i].psi, t);
                mol.backbone.angle[i] = {phi, psi};
            }
            break;
        }
        default:
            ASSERT(false);
        }
    }

    if (mol.backbone.secondary_structure) {
        const md_secondary_structure_t* src_ss[4] = {
            (md_secondary_structure_t*)data->trajectory_data.secondary_structure.data + data->trajectory_data.secondary_structure.stride * frames[0],
            (md_secondary_structure_t*)data->trajectory_data.secondary_structure.data + data->trajectory_data.secondary_structure.stride * frames[1],
            (md_secondary_structure_t*)data->trajectory_data.secondary_structure.data + data->trajectory_data.secondary_structure.stride * frames[2],
            (md_secondary_structure_t*)data->trajectory_data.secondary_structure.data + data->trajectory_data.secondary_structure.stride * frames[3]
        };

        switch (mode) {
        case InterpolationMode::Nearest: {
            const md_secondary_structure_t* ss = t < 0.5f ? src_ss[1] : src_ss[2];
            memcpy(mol.backbone.secondary_structure, ss, mol.backbone.count * sizeof(md_secondary_structure_t));
            break;
        }
        case InterpolationMode::Linear: {
            for (int64_t i = 0; i < mol.backbone.count; ++i) {
                const vec4_t ss_f[2] = {
                    convert_color((uint32_t)src_ss[0][i]),
                    convert_color((uint32_t)src_ss[1][i]),
                };
                const vec4_t ss_res = vec4_lerp(ss_f[0], ss_f[1], t);
                mol.backbone.secondary_structure[i] = (md_secondary_structure_t)convert_color(ss_res);
            }
            break;
        }
        case InterpolationMode::Cubic: {
            for (int64_t i = 0; i < mol.backbone.count; ++i) {
                const vec4_t ss_f[4] = {
                    convert_color((uint32_t)src_ss[0][i]),
                    convert_color((uint32_t)src_ss[1][i]),
                    convert_color((uint32_t)src_ss[2][i]),
                    convert_color((uint32_t)src_ss[3][i]),
                };
                const vec4_t ss_res = cubic_spline(ss_f[0], ss_f[1], ss_f[2], ss_f[3], t);
                mol.backbone.secondary_structure[i] = (md_secondary_structure_t)convert_color(ss_res);
            }
            break;
        }
        default:
            ASSERT(false);
        }
    }

    data->mold.dirty_buffers |= MolBit_DirtyPosition;
    data->mold.dirty_buffers |= MolBit_DirtySecondaryStructure;
}

// #misc
static double compute_avg_ms(double dt) {
    constexpr double interval = 0.5f;
    static double avg = 0.f;
    static int num_frames = 0;
    static double t = 0;
    t += dt;
    num_frames++;

    if (t > interval) {
        avg = t / num_frames * 1000.0;
        t = 0;
        num_frames = 0;
    }

    return avg;
}

static void update_view_param(ApplicationData* data) {
    ViewParam& param = data->view.param;
    param.matrix.previous = param.matrix.current;
    param.jitter.previous = param.jitter.current;

    param.clip_planes.near = data->view.camera.near_plane;
    param.clip_planes.far = data->view.camera.far_plane;
    param.fov_y = data->view.camera.fov_y;
    param.resolution = {(float)data->gbuffer.width, (float)data->gbuffer.height};

    param.matrix.current.view = camera_world_to_view_matrix(data->view.camera);
    if (data->view.mode == CameraMode::Perspective) {
        param.matrix.current.proj = camera_perspective_projection_matrix(data->view.camera, data->gbuffer.width, data->gbuffer.height);
    } else {
        const float aspect_ratio = (float)data->gbuffer.width / (float)data->gbuffer.height;
        const float h = data->view.camera.focus_distance * tanf(data->view.camera.fov_y * 0.5f);
        const float w = aspect_ratio * h;
        param.matrix.current.proj = camera_orthographic_projection_matrix(-w, w, -h, h, data->view.camera.near_plane, data->view.camera.far_plane);
    }
    param.matrix.current.proj_jittered = param.matrix.current.proj;

    if (data->visuals.temporal_reprojection.enabled && data->visuals.temporal_reprojection.jitter) {
        static uint32_t i = 0;
        i = (++i) % ARRAY_SIZE(data->view.jitter.sequence);
        param.jitter.next    = data->view.jitter.sequence[(i + 1) % ARRAY_SIZE(data->view.jitter.sequence)] - 0.5f;
        param.jitter.current = data->view.jitter.sequence[i] - 0.5f;
        if (data->view.mode == CameraMode::Perspective) {
            param.matrix.current.proj_jittered = camera_perspective_projection_matrix(data->view.camera, data->gbuffer.width, data->gbuffer.height,
                param.jitter.current.x, param.jitter.current.y);
        } else {
            const float aspect_ratio = (float)data->gbuffer.width / (float)data->gbuffer.height;
            const float h = data->view.camera.focus_distance * tanf(data->view.camera.fov_y * 0.5f);
            const float w = aspect_ratio * h;
            const float scale_x = w / data->gbuffer.width * 2.0f;
            const float scale_y = h / data->gbuffer.height * 2.0f;
            const float j_x = param.jitter.current.x * scale_x;
            const float j_y = param.jitter.current.y * scale_y;
            param.matrix.current.proj_jittered = param.matrix.current.proj =
                camera_orthographic_projection_matrix(-w + j_x, w + j_x, -h + j_y, h + j_y, data->view.camera.near_plane, data->view.camera.far_plane);
        }
    }

    param.matrix.current.view_proj = param.matrix.current.proj * param.matrix.current.view;
    param.matrix.current.view_proj_jittered = param.matrix.current.proj_jittered * param.matrix.current.view;

    param.matrix.inverse.view = mat4_inverse(param.matrix.current.view);
    param.matrix.inverse.proj = mat4_inverse(param.matrix.current.proj);
    param.matrix.inverse.proj_jittered = mat4_inverse(param.matrix.current.proj_jittered);
    param.matrix.inverse.view_proj = mat4_inverse(param.matrix.current.view_proj);
    param.matrix.inverse.view_proj_jittered = mat4_inverse(param.matrix.current.view_proj_jittered);

    param.matrix.current.norm = mat4_transpose(param.matrix.inverse.view);
}

static void reset_view(ApplicationData* data, bool move_camera, bool smooth_transition) {
    ASSERT(data);
    if (!data->mold.mol.atom.count) return;
    const auto& mol = data->mold.mol;

    vec3_t aabb_min, aabb_max;
    compute_aabb(&aabb_min, &aabb_max, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count);
    const vec3_t ext = aabb_max - aabb_min;
    const vec3_t cen = (aabb_min + aabb_max) * 0.5f;
    const vec3_t pos = cen + ext * 3.f;

    if (move_camera) {
        if (!smooth_transition) data->view.camera.position = pos;
        data->view.animation.target_position = pos;
        data->view.camera.focus_distance = vec3_length(pos - cen);
        data->view.camera.orientation = (quat_from_mat4(look_at(data->view.animation.target_position, cen, {0, 1, 0})));
    }

    data->view.camera.near_plane = 1.f;
    data->view.camera.far_plane = vec3_length(ext) * 50.f;
    data->view.trackball_param.max_distance = vec3_length(ext) * 20.0f;
}

// #picking
static PickingData read_picking_data(GBuffer* gbuf, int32_t x, int32_t y) {
    ASSERT(gbuf);
    uint32_t frame = gbuf->pbo_picking.frame++;
    uint32_t curr = (frame + 0) % 2;
    uint32_t prev = (frame + 1) % 2;

    PickingData data{};

    PUSH_GPU_SECTION("READ PICKING DATA")
    glBindFramebuffer(GL_READ_FRAMEBUFFER, gbuf->deferred.fbo);
    glReadBuffer(GL_COLOR_ATTACHMENT_PICKING);

    // Queue async reads from current frame to pixel pack buffer
    glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.color[curr]);
    glReadPixels(x, y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.depth[curr]);
    glReadPixels(x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, 0);

    // Read values from previous frames pixel pack buffer
    glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.color[prev]);
    const GLubyte* color = (const GLubyte*)glMapBuffer(GL_PIXEL_PACK_BUFFER, GL_READ_ONLY);
    if (color) {
        data.idx = color[0] + (color[1] << 8) + (color[2] << 16) + (color[3] << 24);
        glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
    }

    glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.depth[prev]);
    const GLfloat* depth = (const GLfloat*)glMapBuffer(GL_PIXEL_PACK_BUFFER, GL_READ_ONLY);
    if (depth) {
        data.depth = depth[0];
        glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
    }

    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    POP_GPU_SECTION()

    return data;
}

static void compute_aabb(vec3_t* aabb_min, vec3_t* aabb_max, const float* x, const float* y, const float* z, int64_t count) {
    ASSERT(count >= 0);

    if (count < 1) {
        *aabb_min = *aabb_max = {0,0,0};
        return;
    }

    *aabb_min = {x[0], y[0], z[0]};
    *aabb_max = {x[0], y[0], z[0]};
    for (int64_t i = 1; i < count; ++i) {
        aabb_min->x = MIN(aabb_min->x, x[i]);
        aabb_max->x = MAX(aabb_max->x, x[i]);
        aabb_min->y = MIN(aabb_min->y, y[i]);
        aabb_max->y = MAX(aabb_max->y, y[i]);
        aabb_min->z = MIN(aabb_min->z, z[i]);
        aabb_max->z = MAX(aabb_max->z, z[i]);
    }
}

static void grow_mask_by_covalent_bond(md_exp_bitfield_t* mask, md_bond_t* bonds, int64_t num_bonds, int64_t extent) {
    md_exp_bitfield_t prev_mask;
    md_bitfield_init(&prev_mask, frame_allocator);
    defer { md_bitfield_free(&prev_mask); };

    for (int64_t i = 0; i < extent; i++) {
        for (int64_t j = 0; j < num_bonds; ++j) {
            const auto& bond = bonds[j];
            const int64_t idx[2] = {bond.idx[0], bond.idx[1]};
            if (md_bitfield_test_bit(&prev_mask, idx[0]) && !md_bitfield_test_bit(mask, idx[1])) {
                md_bitfield_set_bit(mask, idx[1]);
            } else if (md_bitfield_test_bit(&prev_mask, idx[1]) && !md_bitfield_test_bit(mask, idx[0])) {
                md_bitfield_set_bit(mask, idx[0]);
            }
        }
        md_bitfield_copy(&prev_mask, mask);
    }
}

static void grow_mask_by_radial_extent(md_exp_bitfield_t* dst_mask, const md_exp_bitfield_t* src_mask, const float* x, const float* y, const float* z, int64_t count, float ext) {
    if (ext > 0.0f) {
        const float cell_ext = CLAMP(ext / 3.0f, 3.0f, 12.0f);
        md_spatial_hash_t ctx = {0};
        md_spatial_hash_init(&ctx, x, y, z, count, cell_ext, frame_allocator);
        defer { md_spatial_hash_free(&ctx, frame_allocator); };

        md_bitfield_clear(dst_mask);

        int64_t beg_bit = src_mask->beg_bit;
        int64_t end_bit = src_mask->end_bit;
        while ((beg_bit = md_bitfield_scan(src_mask, beg_bit, end_bit)) != 0) {
            int64_t i = beg_bit - 1;
            vec4_t pos_rad = {x[i], y[i], z[i], ext};
            md_spatial_hash_query(&ctx, pos_rad, [](uint32_t idx, vec3_t, void* user_data) -> bool {
                md_exp_bitfield_t* mask = (md_exp_bitfield_t*)user_data;
                md_bitfield_set_bit(mask, (int64_t)idx);
                return true;
            }, dst_mask);
        }
    } else {
        md_bitfield_copy(dst_mask, src_mask);
    }
}

static void expand_mask(md_exp_bitfield_t* mask, const md_range_t ranges[], int64_t num_ranges) {
    for (int64_t i = 0; i < num_ranges; i++) {
        if (md_bitfield_popcount_range( mask, ranges[i].beg, ranges[i].end) != 0) {
            md_bitfield_set_range(mask, ranges[i].beg, ranges[i].end);
        }
    }
}

static bool filter_expression(ApplicationData* data, str_t expr, md_exp_bitfield_t* mask, bool* is_dynamic = NULL, char* error_str = NULL, int64_t error_cap = 0) {
    if (data->mold.mol.atom.count == 0) return false;    
    
    md_filter_result_t res = {0};
    bool success = false;

    md_semaphore_aquire(&data->mold.script.ir_semaphore);
    if (md_filter_evaluate(&res, expr, &data->mold.mol, &data->mold.script.ir, frame_allocator)) {
        if (is_dynamic) {
            *is_dynamic = res.is_dynamic;
        }
        if (mask) {
            md_bitfield_clear(mask);
            for (int64_t i = 0; i < res.num_bitfields; ++i) {
                md_bitfield_or_inplace(mask, &res.bitfields[i]);
            }
        }
        success = true;
    } else {
        if (error_str) {
            snprintf(error_str, error_cap, "%s", res.error_buf);
        }
    }
    md_semaphore_release(&data->mold.script.ir_semaphore);

    return success;
}

static bool valid_identifier(str_t str) {
    if (!str.ptr) return false;
    if (!str.len) return false;

    const char* beg = str.ptr;
    const char* end = str.ptr + str.len;

    if (!is_alpha(*beg) && *beg != '_') return false;
    for (const char* c = beg + 1; c < end; ++c) {
        if (!is_alpha(*c) && (*c != '_') && !is_digit(*c)) return false;
    }

    return true;
}

// ### DRAW WINDOWS ###
static void draw_main_menu(ApplicationData* data) {
    ASSERT(data);
    bool new_clicked = false;

    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("Load Data", "CTRL+L")) {
                auto res = application::file_dialog(application::FileDialogFlags_Open, {}, make_cstr("pdb,gro,xtc"));
                if (res.result == application::FileDialogResult::Ok) {
                    load_dataset_from_file(data, str_from_cstr(res.path));
                    if (!data->representations.buffer) {
                        create_representation(data); // Create default representation
                    }
                    data->animation = {};
                    reset_view(data, true);
                }
            }
            if (ImGui::MenuItem("Open Workspace", "CTRL+O")) {
                auto res = application::file_dialog(application::FileDialogFlags_Open, {}, FILE_EXTENSION);
                if (res.result == application::FileDialogResult::Ok) {
                    load_workspace(data, str_from_cstr(res.path));
                }
            }
            if (ImGui::MenuItem("Save Workspace", "CTRL+S")) {
                if (!data->files.workspace) {
                    auto res = application::file_dialog(application::FileDialogFlags_Save, {}, FILE_EXTENSION);
                    if (res.result == application::FileDialogResult::Ok) {
                        str_t ext = extract_ext(str_from_cstr(res.path));
                        if (!ext.len) {
                            snprintf(res.path + res.path_len, ARRAY_SIZE(res.path) - res.path_len, ".%.*s", (uint32_t)FILE_EXTENSION.len, FILE_EXTENSION.ptr);
                        }
                        save_workspace(data, str_from_cstr(res.path));
                    }
                } else {
                    save_workspace(data, data->files.workspace);
                }
            }
            if (ImGui::MenuItem("Save As")) {
                auto res = application::file_dialog(application::FileDialogFlags_Save, {}, FILE_EXTENSION);
                if (res.result == application::FileDialogResult::Ok) {
                    str_t ext = extract_ext(str_from_cstr(res.path));
                    if (!ext.len) {
                        snprintf(res.path + res.path_len, (size_t)(ARRAY_SIZE(res.path) - res.path_len), "%.*s", (uint32_t)FILE_EXTENSION.len, FILE_EXTENSION.ptr);
                    }
                    save_workspace(data, str_from_cstr(res.path));
                }
            }
            ImGui::Separator();
            if (ImGui::MenuItem("Quit", "ALT+F4")) {
                data->ctx.window.should_close = true;
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
                reset_view(data, true, true);
            }
            ImGui::Separator();
            ImGui::Checkbox("Vsync", &data->ctx.window.vsync);
            ImGui::Separator();

            ImGui::BeginGroup();
            ImGui::Text("Camera");
            {
                ImGui::Combo("Mode", (int*)(&data->view.mode), "Perspective\0Orthographic\0");
                if (data->view.mode == CameraMode::Perspective) {
                    float fov = rad_to_deg(data->view.camera.fov_y);
                    if (ImGui::SliderFloat("field of view", &fov, 12.5f, 80.0f)) {
                        data->view.camera.fov_y = deg_to_rad(fov);
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
            // Temporal
            ImGui::BeginGroup();
            {
                ImGui::Checkbox("Temporal Effects", &data->visuals.temporal_reprojection.enabled);
                if (data->visuals.temporal_reprojection.enabled) {
                    // ImGui::Checkbox("Jitter Samples", &data->visuals.temporal_reprojection.jitter);
                    // ImGui::SliderFloat("Feedback Min", &data->visuals.temporal_reprojection.feedback_min, 0.5f, 1.0f);
                    // ImGui::SliderFloat("Feedback Max", &data->visuals.temporal_reprojection.feedback_max, 0.5f, 1.0f);
                    ImGui::Checkbox("Motion Blur", &data->visuals.temporal_reprojection.motion_blur.enabled);
                    if (data->visuals.temporal_reprojection.motion_blur.enabled) {
                        ImGui::SliderFloat("Motion Scale", &data->visuals.temporal_reprojection.motion_blur.motion_scale, 0.f, 2.0f);
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
            /*
// DEBUG DRAW
ImGui::BeginGroup();
ImGui::Checkbox("Spline", &data->debug_draw.spline.enabled);
ImGui::Checkbox("Backbone", &data->debug_draw.backbone.enabled);
ImGui::EndGroup();

            */
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Windows")) {
            ImGui::Checkbox("Representations", &data->representations.show_window);
            ImGui::Checkbox("Timelines", &data->timeline.show_window);
            ImGui::Checkbox("Distributions", &data->distributions.show_window);
            ImGui::Checkbox("Density Volumes", &data->density_volume.show_window);
            //ImGui::Checkbox("Ramachandran", &data->ramachandran.show_window);
            ImGui::Checkbox("Shape Space", &data->shape_space.show_window);
            ImGui::Checkbox("Dataset", &data->dataset.show_window);

            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Selection")) {
            ImGui::Combo("Granularity", (int*)(&data->selection.granularity), "Atom\0Residue\0Chain\0\0");
            int64_t num_selected_atoms = md_bitfield_popcount(&data->selection.current_selection_mask);
            if (num_selected_atoms == 0) ImGui::PushDisabled();
            if (ImGui::MenuItem("Invert")) {
                md_bitfield_not_inplace(&data->selection.current_selection_mask, 0, data->mold.mol.atom.count);
                data->mold.dirty_buffers |= MolBit_DirtyFlags;
            }
            if (ImGui::MenuItem("Grow"))  data->selection.grow.show_window = true;
            if (num_selected_atoms == 0) ImGui::PopDisabled();

            if (ImGui::MenuItem("Query")) data->selection.query.show_window = true;

            ImGui::Spacing();
            ImGui::Separator();

            // STORED SELECTIONS
            {
                // @NOTE(Robin): This ImGui ItemFlag can be used to force the menu to remain open after buttons are pressed.
                // Leave it here as a comment if we feel that it is needed in the future
                //ImGui::PushItemFlag(ImGuiItemFlags_SelectableDontClosePopup, true);
                ImGui::Text("Stored Selections");
                for (int i = 0; i < (int)md_array_size(data->selection.stored_selections); i++) {
                    auto& sel = data->selection.stored_selections[i];
                    bool is_valid = valid_identifier(sel.name);
                    char err_buf[64] = "";
                    if (!is_valid) {
                        snprintf(err_buf, sizeof(err_buf), "'%s' is not a valid identifier.", sel.name.cstr());
                    }

                    for (int j = 0; j < i; ++j) {
                        if (compare_str(sel.name, data->selection.stored_selections[j].name)) {
                            is_valid = false;
                            snprintf(err_buf, sizeof(err_buf), "identifier '%s' is already taken.", sel.name.cstr());
                            break;
                        }
                    }

                    ImGui::PushID(i);
                    if (!is_valid) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
                    ImGui::InputText("##label", sel.name.beg(), sel.name.capacity());
                    if (!is_valid) ImGui::PopStyleColor();
                    if (ImGui::IsItemHovered() && !is_valid && err_buf[0] != '\0') {
                        ImGui::SetTooltip("%s", err_buf);
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Load")) {
                        md_bitfield_copy(&data->selection.current_selection_mask, &sel.atom_mask);
                        update_all_representations(data);
                    }
                    if (ImGui::IsItemHovered()) {
                        ImGui::SetTooltip("Load the stored selection into the active selection");
                        md_bitfield_copy(&data->selection.current_highlight_mask, &sel.atom_mask);
                        data->mold.dirty_buffers |= MolBit_DirtyFlags;
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Store")) {
                        ImGui::SetTooltip("Store the active selection into the stored selection");
                        md_bitfield_copy(&sel.atom_mask, &data->selection.current_selection_mask);
                        update_all_representations(data);
                    }
                    ImGui::SameLine();
                    if (ImGui::DeleteButton("Remove")) {
                        ImGui::SetTooltip("Remove the stored selection");
                        remove_selection(data, i);
                    }
                    ImGui::PopID();
                }

                if (ImGui::Button("Create New")) {
                    char name_buf[64];
                    snprintf(name_buf, sizeof(name_buf), "sel%lli", md_array_size(data->selection.stored_selections) + 1);
                    create_selection(data, str_from_cstr(name_buf), &data->selection.current_selection_mask);
                }

                //ImGui::PopItemFlag();
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Screenshot")) {
            if (ImGui::MenuItem("Take Screenshot")) {
                create_screenshot(data);
            }
            ImGui::EndMenu();
        }
        {
            // Fps counter
            const double ms = compute_avg_ms(data->ctx.timing.delta_s);
            char fps_buf[32] = {0};
            snprintf(fps_buf, ARRAY_SIZE(fps_buf), "%.2f ms (%.1f fps)", ms, 1000.f / ms);
            const float w = ImGui::CalcTextSize(fps_buf).x;
            ImGui::SetCursorPosX(ImGui::GetWindowContentRegionMax().x - w);
            ImGui::Text("%s", fps_buf);
        }
        ImGui::EndMainMenuBar();
    }

    if (new_clicked) ImGui::OpenPopup("Warning New");
}

void clear_atom_elem_mappings(ApplicationData* data) {
    md_array_shrink(data->dataset.atom_element_remappings, 0);
}

AtomElementMapping* add_atom_elem_mapping(ApplicationData* data, str_t lbl, md_element_t elem) {
    // Check if we already have a mapping for the label -> overwrite
    int64_t i = 0;
    for (; i < md_array_size(data->dataset.atom_element_remappings); ++i) {
        if (compare_str(lbl, data->dataset.atom_element_remappings[i].lbl)) break;
    }
    if (i == md_array_size(data->dataset.atom_element_remappings)) {
        AtomElementMapping mapping = {
            .lbl = lbl,
            .elem = elem,
        };
        return md_array_push(data->dataset.atom_element_remappings, mapping, persistent_allocator);
    } else {
        data->dataset.atom_element_remappings[i].elem = elem;
        return &data->dataset.atom_element_remappings[i];
    }
}

void apply_atom_elem_mappings(ApplicationData* data) {
    for (int64_t j = 0; j < md_array_size(data->dataset.atom_element_remappings); ++j) {
        str_t lbl = data->dataset.atom_element_remappings[j].lbl;
        md_element_t elem = data->dataset.atom_element_remappings[j].elem;
        float radius = md_util_element_vdw_radius(elem);

        for (int64_t i = 0; i < data->mold.mol.atom.count; ++i) {
            if (compare_str_cstr(lbl, data->mold.mol.atom.name[i])) {
                data->mold.mol.atom.element[i] = elem;
                data->mold.mol.atom.radius[i] = radius;
                data->mold.dirty_buffers |= MolBit_DirtyRadius;
            }
        }
    }

    auto& mol = data->mold.mol;
    md_util_covalent_bond_args_t args = {
        .atom = {
            .count = mol.atom.count,
            .x = mol.atom.x,
            .y = mol.atom.y,
            .z = mol.atom.z,
            .element = mol.atom.element,
        },
        .residue = {
            .count = mol.residue.count,
            .atom_range = mol.residue.atom_range,
            .internal_bond_range = mol.residue.internal_covalent_bond_range,
            .complete_bond_range = mol.residue.complete_covalent_bond_range,
        }
    };
    mol.covalent_bond.bond = md_util_extract_covalent_bonds(&args, data->mold.mol_alloc);
    mol.covalent_bond.count = md_array_size(mol.covalent_bond.bond);
    data->mold.dirty_buffers |= MolBit_DirtyBonds;

    update_all_representations(data);
}

// # context_menu
void draw_context_popup(ApplicationData* data) {
    ASSERT(data);

    if (!data->mold.mol.atom.count) return;

    const bool shift_down = ImGui::GetIO().KeyShift;
    const int64_t sss_count = single_selection_sequence_count(&data->selection.single_selection_sequence);
    const int64_t num_frames = md_trajectory_num_frames(&data->mold.traj);
    const int64_t num_atoms_selected = md_bitfield_popcount(&data->selection.current_selection_mask);

    if (data->ctx.input.mouse.clicked[1] && !shift_down && !ImGui::GetIO().WantTextInput) {
        ImGui::OpenPopup("AtomContextPopup");
    }

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
        if (sss_count > 1) {
            if (ImGui::BeginMenu("Create Property")) {
                char ident[128] = "prop";
                char buf[128] = "";
                ImGui::InputText("##identifier", ident, sizeof(ident));

                if (sss_count == 2) {
                    int32_t idx[2] = {data->selection.single_selection_sequence.idx[0], data->selection.single_selection_sequence.idx[1]};

                    snprintf(buf, sizeof(buf), "%s = distance(%i, %i);", ident, idx[0]+1, idx[1]+1);
                    if (ImGui::MenuItem(buf)) {
                        data->script.editor.AppendText("\n");
                        data->script.editor.AppendText(buf);
                        ImGui::CloseCurrentPopup();
                    }

                    if (data->mold.mol.atom.residue_idx[idx[0]] == data->mold.mol.atom.residue_idx[idx[1]]) {
                        int32_t res_idx = data->mold.mol.atom.residue_idx[idx[0]];
                        idx[0] -= data->mold.mol.residue.atom_range[res_idx].beg;
                        idx[1] -= data->mold.mol.residue.atom_range[res_idx].beg;

                        snprintf(buf, sizeof(buf), "%s = distance(%i, %i) in residue(%i);", ident, idx[0]+1, idx[1]+1, res_idx+1);
                        if (ImGui::MenuItem(buf)) {
                            data->script.editor.AppendText("\n");
                            data->script.editor.AppendText(buf);
                            ImGui::CloseCurrentPopup();
                        }

                        int32_t resid = data->mold.mol.residue.id[idx[0]];
                        snprintf(buf, sizeof(buf), "%s = distance(%i, %i) in resid(%i);", ident, idx[0]+1, idx[1]+1, resid);
                        if (ImGui::MenuItem(buf)) {
                            data->script.editor.AppendText("\n");
                            data->script.editor.AppendText(buf);
                            ImGui::CloseCurrentPopup();
                        }

                        const char* resname = data->mold.mol.residue.name[res_idx];
                        snprintf(buf, sizeof(buf), "%s = distance(%i, %i) in resname(%s);", ident, idx[0]+1, idx[1]+1, resname);
                        if (ImGui::MenuItem(buf)) {
                            data->script.editor.AppendText("\n");
                            data->script.editor.AppendText(buf);
                            ImGui::CloseCurrentPopup();
                        }
                    }
                }
                else if(sss_count == 3) {
                    int32_t idx[3] = {data->selection.single_selection_sequence.idx[0], data->selection.single_selection_sequence.idx[1], data->selection.single_selection_sequence.idx[2]};

                    snprintf(buf, sizeof(buf), "%s = angle(%i, %i, %i);", ident, idx[0]+1, idx[1]+1, idx[2]+1);
                    if (ImGui::MenuItem(buf)) {
                        data->script.editor.AppendText("\n");
                        data->script.editor.AppendText(buf);
                        ImGui::CloseCurrentPopup();
                    }

                    if (data->mold.mol.atom.residue_idx[idx[0]] == data->mold.mol.atom.residue_idx[idx[1]] &&
                        data->mold.mol.atom.residue_idx[idx[0]] == data->mold.mol.atom.residue_idx[idx[2]]) {

                        int32_t res_idx = data->mold.mol.atom.residue_idx[idx[0]];
                        idx[0] -= data->mold.mol.residue.atom_range[res_idx].beg;
                        idx[1] -= data->mold.mol.residue.atom_range[res_idx].beg;
                        idx[2] -= data->mold.mol.residue.atom_range[res_idx].beg;

                        snprintf(buf, sizeof(buf), "%s = angle(%i, %i, %i) in residue(%i);", ident, idx[0]+1, idx[1]+1, idx[2]+1, res_idx+1);
                        if (ImGui::MenuItem(buf)) {
                            data->script.editor.AppendText("\n");
                            data->script.editor.AppendText(buf);
                            ImGui::CloseCurrentPopup();
                        }

                        int32_t resid = data->mold.mol.residue.id[idx[0]];
                        snprintf(buf, sizeof(buf), "%s = angle(%i, %i, %i) in resid(%i);", ident, idx[0]+1, idx[1]+1, idx[2]+1, resid);
                        if (ImGui::MenuItem(buf)) {
                            data->script.editor.AppendText("\n");
                            data->script.editor.AppendText(buf);
                            ImGui::CloseCurrentPopup();
                        }

                        const char* resname = data->mold.mol.residue.name[res_idx];
                        snprintf(buf, sizeof(buf), "%s = angle(%i, %i, %i) in resname(%s);", ident, idx[0]+1, idx[1]+1, idx[2]+1, resname);
                        if (ImGui::MenuItem(buf)) {
                            data->script.editor.AppendText("\n");
                            data->script.editor.AppendText(buf);
                            ImGui::CloseCurrentPopup();
                        }
                    }
                }
                else if(sss_count == 4) {
                    int32_t idx[4] = {data->selection.single_selection_sequence.idx[0], data->selection.single_selection_sequence.idx[1], data->selection.single_selection_sequence.idx[2], data->selection.single_selection_sequence.idx[3]};

                    snprintf(buf, sizeof(buf), "%s = dihedral(%i, %i, %i, %i);", ident, idx[0]+1, idx[1]+1, idx[2]+1, idx[3]+1);
                    if (ImGui::MenuItem(buf)) {
                        data->script.editor.AppendText("\n");
                        data->script.editor.AppendText(buf);
                        ImGui::CloseCurrentPopup();
                    }

                    if (data->mold.mol.atom.residue_idx[idx[0]] == data->mold.mol.atom.residue_idx[idx[1]] &&
                        data->mold.mol.atom.residue_idx[idx[0]] == data->mold.mol.atom.residue_idx[idx[2]] &&
                        data->mold.mol.atom.residue_idx[idx[0]] == data->mold.mol.atom.residue_idx[idx[3]]) {

                        int32_t res_idx = data->mold.mol.atom.residue_idx[idx[0]];
                        idx[0] -= data->mold.mol.residue.atom_range[res_idx].beg;
                        idx[1] -= data->mold.mol.residue.atom_range[res_idx].beg;
                        idx[2] -= data->mold.mol.residue.atom_range[res_idx].beg;
                        idx[3] -= data->mold.mol.residue.atom_range[res_idx].beg;

                        snprintf(buf, sizeof(buf), "%s = dihedral(%i, %i, %i, %i) in residue(%i);", ident, idx[0]+1, idx[1]+1, idx[2]+1, idx[3]+1, res_idx+1);
                        if (ImGui::MenuItem(buf)) {
                            data->script.editor.AppendText("\n");
                            data->script.editor.AppendText(buf);
                            ImGui::CloseCurrentPopup();
                        }

                        int32_t resid = data->mold.mol.residue.id[idx[0]];
                        snprintf(buf, sizeof(buf), "%s = dihedral(%i, %i, %i, %i) in resid(%i);", ident, idx[0]+1, idx[1]+1, idx[2]+1, idx[3]+1, resid);
                        if (ImGui::MenuItem(buf)) {
                            data->script.editor.AppendText("\n");
                            data->script.editor.AppendText(buf);
                            ImGui::CloseCurrentPopup();
                        }

                        const char* resname = data->mold.mol.residue.name[res_idx];
                        snprintf(buf, sizeof(buf), "%s = dihedral(%i, %i, %i, %i) in resname(%s);", ident, idx[0]+1, idx[1]+1, idx[2]+1, idx[3]+1, resname);
                        if (ImGui::MenuItem(buf)) {
                            data->script.editor.AppendText("\n");
                            data->script.editor.AppendText(buf);
                            ImGui::CloseCurrentPopup();
                        }
                    }
                }
                ImGui::EndMenu();
            }
        }
        if (data->selection.right_clicked != -1 && num_atoms_selected == 0) {
            int idx = data->selection.right_clicked;
            if (0 <= idx && idx < data->mold.mol.atom.count) {
                char label[64] = "";
                snprintf(label, sizeof(label), "Remap Element for '%s'", data->mold.mol.atom.name[idx]);
                if (ImGui::BeginMenu(label)) {
                    static char input_buf[32] = "";
                    md_element_t elem = data->mold.mol.atom.element[idx];
                    str_t lbl = {data->mold.mol.atom.name[idx], (int64_t)strlen(data->mold.mol.atom.name[idx])};
                    str_t name = md_util_element_name(elem);
                    str_t sym  = md_util_element_symbol(elem);

                    ImGui::Text("Current Element: %.*s (%.*s)", (int)name.len, name.ptr, (int)sym.len, sym.ptr);

                    str_t elem_str = {input_buf, (int64_t)strnlen(input_buf, sizeof(input_buf))};
                    md_element_t new_elem = md_util_lookup_element(elem_str);
                    const bool is_valid = new_elem != 0;

                    if (!is_valid) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
                    ImGui::InputText("##Symbol", input_buf, sizeof(input_buf));
                    if (!is_valid) ImGui::PopStyleColor();
                    str_t new_name = md_util_element_name(new_elem);
                    str_t new_sym  = md_util_element_symbol(new_elem);
                    ImGui::Text("New Element: %.*s (%.*s)", (int)new_name.len, new_name.ptr, (int)new_sym.len, new_sym.ptr);
                    if (!is_valid) ImGui::PushDisabled();
                    if (ImGui::Button("Apply") && is_valid) {
                        add_atom_elem_mapping(data, lbl, new_elem);
                        apply_atom_elem_mappings(data);
                        ImGui::CloseCurrentPopup();
                    }
                    if (!is_valid) ImGui::PopDisabled();
                    ImGui::EndMenu();
                }
            }
        }
        if (data->selection.right_clicked != -1 && num_frames > 0) {
            if (ImGui::BeginMenu("Recenter Trajectory...")) {
                const int idx = data->selection.right_clicked;

                md_exp_bitfield_t mask = {0};
                md_bitfield_init(&mask, frame_allocator);
                bool apply = false;

                apply |= ImGui::MenuItem("on Atom");
                if (ImGui::IsItemHovered()) {
                    md_bitfield_set_bit(&mask, idx);
                }

                apply |= ImGui::MenuItem("on Residue");
                if (ImGui::IsItemHovered()) {
                    const auto res_idx = data->mold.mol.atom.residue_idx[idx];
                    const auto range = data->mold.mol.residue.atom_range[res_idx];
                    md_bitfield_set_range(&mask, range.beg, range.end);
                }

                if (data->mold.mol.atom.chain_idx && data->mold.mol.atom.chain_idx[idx] != -1) {
                    apply |= ImGui::MenuItem("on Chain");
                    if (ImGui::IsItemHovered()) {
                        const auto chain_idx = data->mold.mol.atom.chain_idx[idx];
                        const auto range = data->mold.mol.chain.atom_range[chain_idx];
                        md_bitfield_set_range(&mask, range.beg, range.end);
                    }
                }

                if (num_atoms_selected > 0) {
                    apply |= ImGui::MenuItem("on Selection");
                    if (ImGui::IsItemHovered()) {
                        md_bitfield_copy(&mask, &data->selection.current_selection_mask);
                    }
                }

                if (!md_bitfield_empty(&mask)) {
                    md_bitfield_copy(&data->selection.current_highlight_mask, &mask);
                    data->mold.dirty_buffers |= MolBit_DirtyFlags;

                    if (apply) {
                        load::traj::set_recenter_target(&data->mold.traj, &mask);
                        load::traj::clear_cache(&data->mold.traj);
                        launch_prefetch_job(data);
                        interpolate_atomic_properties(data);
                        data->mold.dirty_buffers |= MolBit_DirtyPosition;
                        update_md_buffers(data);
                        md_gl_molecule_update_atom_previous_position(&data->mold.gl_mol); // Do this explicitly to update the previous position to avoid motion blur trails
                        ImGui::CloseCurrentPopup();
                    }
                }
                ImGui::EndMenu();
            }
        }
        if (ImGui::BeginMenu("Selection")) {
            if (num_atoms_selected > 0) {
                if (ImGui::MenuItem("Invert")) {
                    md_bitfield_not_inplace(&data->selection.current_selection_mask, 0, data->mold.mol.atom.count);
                    ImGui::CloseCurrentPopup();
                }
                if (ImGui::MenuItem("Grow")) {
                    data->selection.grow.show_window = true;
                    ImGui::CloseCurrentPopup();
                }
            }
            if (ImGui::MenuItem("Query")) {
                data->selection.query.show_window = true;
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndMenu();
        }
        ImGui::EndPopup();
    }
}


static void draw_selection_grow_window(ApplicationData* data) {
    ImGui::SetNextWindowSize(ImVec2(300,120), ImGuiCond_Always);
    if (ImGui::Begin("Selection Grow", &data->selection.grow.show_window, ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse)) {
        ImGui::PushItemWidth(-1);
        const bool mode_changed = ImGui::Combo("##Mode", (int*)(&data->selection.grow.mode), "Covalent Bond\0Radial\0\0");
        const bool extent_changed = ImGui::SliderFloat("##Extent", &data->selection.grow.extent, 1.0f, 20.f);
        ImGui::PopItemWidth();
        const bool apply = ImGui::Button("Apply");

        data->selection.grow.mask_invalid |= (mode_changed || extent_changed);

        if (data->selection.grow.mask_invalid) {
            data->selection.grow.mask_invalid = false;
            switch (data->selection.grow.mode) {
            case SelectionGrowth::CovalentBond:
                md_bitfield_copy(&data->selection.grow.mask, &data->selection.current_selection_mask);
                grow_mask_by_covalent_bond(&data->selection.grow.mask, data->mold.mol.covalent_bond.bond, data->mold.mol.covalent_bond.count, (int64_t)data->selection.grow.extent);
                break;
            case SelectionGrowth::Radial: {
                const auto& mol = data->mold.mol;
                grow_mask_by_radial_extent(&data->selection.grow.mask, &data->selection.current_selection_mask, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, data->selection.grow.extent);
                break;
            }
            default:
                ASSERT(false);
            }

            switch (data->selection.granularity) {
            case SelectionLevel::Atom:
                // No need to expand the mask
                break;
            case SelectionLevel::Residue:
                expand_mask(&data->selection.grow.mask, data->mold.mol.residue.atom_range, data->mold.mol.residue.count);
                break;
            case SelectionLevel::Chain:
                expand_mask(&data->selection.grow.mask, data->mold.mol.chain.atom_range, data->mold.mol.chain.count);
                break;
            default:
                ASSERT(false);
            }
        }

        const bool show_preview =   (ImGui::GetHoveredID() == ImGui::GetID("##Extent")) ||
                                    (ImGui::GetActiveID()  == ImGui::GetID("##Extent")) ||
                                    (ImGui::GetHoveredID() == ImGui::GetID("Apply"));

        if (show_preview) {
            md_bitfield_copy(&data->selection.current_highlight_mask, &data->selection.grow.mask);
            data->mold.dirty_buffers |= MolBit_DirtyFlags;
        }
        if (apply) {
            md_bitfield_copy(&data->selection.current_selection_mask, &data->selection.grow.mask);
            data->selection.grow.show_window = false;
        }
    }
    ImGui::End();
}

static void draw_selection_query_window(ApplicationData* data) {
    ImGui::SetNextWindowSize(ImVec2(300,100), ImGuiCond_Always);
    if (ImGui::Begin("Selection Query", &data->selection.query.show_window, ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse)) {

        if (ImGui::IsKeyPressed(ImGuiKey_Escape, false)) {
            data->selection.query.show_window = false;
            return;
        }

        ImGui::PushItemWidth(-1);
        if (!data->selection.query.query_ok) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
        bool apply = ImGui::InputText("##query", data->selection.query.buf, sizeof(data->selection.query.buf), ImGuiInputTextFlags_AutoSelectAll | ImGuiInputTextFlags_EnterReturnsTrue);
        if (!data->selection.query.query_ok) ImGui::PopStyleColor();
        ImGui::PopItemWidth();

        data->selection.query.query_invalid |= ImGui::IsItemEdited();
        bool preview = ImGui::IsItemActive() || ImGui::IsItemHovered();

        if (ImGui::IsWindowAppearing()) {
            ImGui::SetKeyboardFocusHere(-1);
        }

        if (!data->selection.query.query_ok && data->selection.query.err_buf[0] != '\0' && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("%s", data->selection.query.err_buf);
        }

        if (!data->selection.query.query_ok) ImGui::PushDisabled();
        apply |= ImGui::Button("Apply");
        if (!data->selection.query.query_ok) ImGui::PopDisabled();

        preview |= ImGui::IsItemHovered();

        if (data->selection.query.query_invalid) {
            data->selection.query.query_invalid = false;
            data->selection.query.query_ok = filter_expression(data, str_from_cstr(data->selection.query.buf), &data->selection.query.mask, NULL, data->selection.query.err_buf, sizeof(data->selection.query.err_buf));

            if (data->selection.query.query_ok) {
                switch (data->selection.granularity) {
                case SelectionLevel::Atom:
                    // No need to expand the mask
                    break;
                case SelectionLevel::Residue:
                    expand_mask(&data->selection.query.mask, data->mold.mol.residue.atom_range, data->mold.mol.residue.count);
                    break;
                case SelectionLevel::Chain:
                    expand_mask(&data->selection.query.mask, data->mold.mol.chain.atom_range, data->mold.mol.chain.count);
                    break;
                default:
                    ASSERT(false);
                }
            } else {
                md_bitfield_clear(&data->selection.query.mask);
            }
        }

        if (preview) {
            md_bitfield_copy(&data->selection.current_highlight_mask, &data->selection.query.mask);
            data->mold.dirty_buffers |= MolBit_DirtyFlags;
        }

        if (apply && data->selection.query.query_ok) {
            md_bitfield_copy(&data->selection.current_selection_mask, &data->selection.query.mask);
            data->mold.dirty_buffers |= MolBit_DirtyFlags;
            data->selection.query.show_window = false;
        }
    }
    ImGui::End();
}

static void draw_animation_control_window(ApplicationData* data) {
    ASSERT(data);
    int num_frames = (int)md_trajectory_num_frames(&data->mold.traj);
    if (num_frames == 0) return;

    ASSERT(data->timeline.x_values);
    ASSERT(md_array_size(data->timeline.x_values) == num_frames);

    if (ImGui::Begin("Animation")) {
        ImGui::Text("Num Frames: %i", num_frames);

        double t   = frame_to_time(data->animation.frame, *data);
        double min = data->timeline.x_values[0];
        double max = data->timeline.x_values[num_frames - 1];
        if (ImGui::SliderScalar("Time", ImGuiDataType_Double, &t, &min, &max, "%.2f")) {
            data->animation.frame = time_to_frame(t, *data);
        }
        ImGui::SliderFloat("Speed", &data->animation.fps, -200.0f, 200.f, "%.2f", ImGuiSliderFlags_Logarithmic);
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Animation Speed in Frames Per Second");
        }
        if (ImGui::Combo("Interp.", (int*)(&data->animation.interpolation), "Nearest\0Linear\0Cubic\0\0")) {
            interpolate_atomic_properties(data);
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Interpolation Method for Atom Positions");
        }
        ImGui::Checkbox("Apply PBC", &data->animation.apply_pbc);
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
    }
    ImGui::End();
}

static void draw_representations_window(ApplicationData* data) {
    ImGui::Begin("Representations", &data->representations.show_window, ImGuiWindowFlags_NoFocusOnAppearing);
    if (ImGui::Button("create new")) {
        create_representation(data);
    }
    ImGui::SameLine();
    if (ImGui::DeleteButton("remove all")) {
        clear_representations(data);
    }
    ImGui::Spacing();
    ImGui::Separator();
    for (int i = 0; i < (int)md_array_size(data->representations.buffer); i++) {
        bool update_visual_rep = false;
        bool update_color = false;
        auto& rep = data->representations.buffer[i];
        const float item_width = CLAMP(ImGui::GetWindowContentRegionWidth() - 90.f, 100.f, 300.f);
        StrBuf<128> name;
        snprintf(name.cstr(), name.capacity(), "%s###ID", rep.name.cstr());

        ImGui::PushID(i);
        if (ImGui::CollapsingHeader(name.cstr())) {
            if (ImGui::Checkbox("enabled", &rep.enabled)) {
                data->representations.atom_visibility_mask_dirty = true;
            }
            ImGui::SameLine();
            ImGui::Checkbox("show in sel.", &rep.show_in_selection);
            ImGui::SameLine();
            if (ImGui::DeleteButton("remove")) {
                remove_representation(data, i);
            }
            ImGui::SameLine();
            if (ImGui::Button("clone")) {
                clone_representation(data, rep);
            }

            ImGui::PushItemWidth(item_width);
            ImGui::InputText("name", rep.name.cstr(), rep.name.capacity());
            if (!rep.filt_is_valid) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
            if (ImGui::InputText("filter", rep.filt.cstr(), rep.filt.capacity())) {
                update_color = true;
            }
            if (ImGui::IsItemHovered() && !rep.filt_is_valid) {
                ImGui::SetTooltip("%s", rep.filt_error.cstr());
            }
            if (!rep.filt_is_valid) ImGui::PopStyleColor();
            if (ImGui::Combo("type", (int*)(&rep.type), "VDW\0Licorice\0Ribbons\0Cartoon\0")) {
                update_visual_rep = true;
                
            }
            if (ImGui::Combo("color", (int*)(&rep.color_mapping),
                             "Uniform Color\0CPK\0Res Id\0Res Idx\0Chain Id\0Chain Idx\0Secondary Structure\0Property\0")) {
                update_color = true;
            }
            if (rep.color_mapping == ColorMapping::Property) {
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
                md_script_property_t* props[32] = {0};
                int num_props = 0;
                for (int64_t j = 0; j < md_array_size(data->display_properties); ++j) {
                    if (data->display_properties[j].full_prop->type == MD_SCRIPT_PROPERTY_TYPE_TEMPORAL) {
                        props[num_props++] = data->display_properties[j].full_prop;
                    }
                }

                rep.prop = NULL;
                if (num_props > 0) {
                    prop_idx = CLAMP(prop_idx, 0, num_props-1);
                    if (ImGui::BeginCombo("Prop", props[prop_idx]->ident.ptr)) {
                        for (int j = 0; j < num_props; ++j) {
                            if (ImGui::Selectable(props[j]->ident.ptr, prop_idx == i)) {
                                prop_idx = j;
                                rep.map_beg = props[j]->data.min_value;
                                rep.map_end = props[j]->data.max_value;
                                update_color = true;
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
            if (rep.filt_is_dynamic || rep.color_mapping == ColorMapping::Property) {
                ImGui::Checkbox("auto-update", &rep.dynamic_evaluation);
                if (!rep.dynamic_evaluation) {
                    ImGui::SameLine();
                    if (ImGui::Button("update")) update_color = true;
                }
            } else {
                rep.dynamic_evaluation = false;
            }
            ImGui::PopItemWidth();
            if (rep.color_mapping == ColorMapping::Uniform) {
                if (ImGui::ColorEdit4("color", (float*)&rep.uniform_color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel)) {
                    update_color = true;
                }
            }
            ImGui::PushItemWidth(item_width);
            if (rep.type == RepresentationType::SpaceFill || rep.type == RepresentationType::Licorice) {
                if (ImGui::SliderFloat("radii scale", &rep.radius, 0.1f, 2.f)) update_visual_rep = true;
            }
            if (rep.type == RepresentationType::Ribbons) {
                if (ImGui::SliderFloat("spline tension", &rep.tension, 0.f, 1.f)) update_visual_rep = true;
                if (ImGui::SliderFloat("spline width", &rep.width, 0.1f, 2.f)) update_visual_rep = true;
                if (ImGui::SliderFloat("spline thickness", &rep.thickness, 0.1f, 2.f)) update_visual_rep = true;
            }
            ImGui::PopItemWidth();
            ImGui::Spacing();
            ImGui::Separator();
        }

        ImGui::PopID();
        if (update_visual_rep) {
            md_gl_representation_type_t type = MD_GL_REP_DEFAULT;
            md_gl_representation_args_t args = {};
            switch(rep.type) {
            case RepresentationType::SpaceFill:
                type = MD_GL_REP_SPACE_FILL;
                args.space_fill.radius_scale = rep.radius;
                break;
            case RepresentationType::Licorice:
                type = MD_GL_REP_LICORICE;
                args.licorice.radius = rep.radius * 0.5f;
                break;
            case RepresentationType::Ribbons:
                type = MD_GL_REP_RIBBONS;
                args.ribbons.width_scale = rep.width;
                args.ribbons.thickness_scale = rep.thickness;
                break;
            case RepresentationType::Cartoon:
                type = MD_GL_REP_CARTOON;
                args.cartoon.width_scale = rep.width;
                args.cartoon.thickness_scale = rep.thickness;
                break;
            default: break;
            }

            md_gl_representation_set_type_and_args(&rep.md_rep, type, args);
        }

        if (update_color) {
            update_representation(data, &rep);
        }
    }

    ImGui::End();
}

static void draw_atom_info_window(const ApplicationData& data, int atom_idx) {
    const auto& mol = data.mold.mol;

    // @TODO: Assert things and make this failproof
    if (atom_idx < 0 || atom_idx >= mol.atom.count) return;

    int res_idx = mol.atom.residue_idx[atom_idx];
    const char* res_name = mol.residue.name[res_idx];
    const int res_id = mol.residue.id[res_idx];
    int local_idx = atom_idx - mol.residue.atom_range[res_idx].beg;
    const vec3_t pos = { mol.atom.x[atom_idx], mol.atom.y[atom_idx], mol.atom.z[atom_idx] };
    const char* label = mol.atom.name[atom_idx];
    str_t elem = md_util_element_name(mol.atom.element[atom_idx]);
    str_t symbol = md_util_element_symbol(mol.atom.element[atom_idx]);

    int chain_idx = -1;
    const char* chain_id = "\0";
    if (mol.atom.chain_idx) {
        chain_idx = mol.atom.chain_idx[atom_idx];
        if (chain_idx != -1 && mol.chain.count > 0) {
            chain_id = mol.chain.id[chain_idx];
        }
    }

    // External indices begin with 1 not 0
    res_idx += 1;
    chain_idx += 1;
    atom_idx += 1;
    local_idx += 1;

    char buff[256];
    int len = 0;
    len += snprintf(buff, ARRAY_SIZE(buff) - 1, "atom[%i][%i]: %s %.*s %.*s (%.2f, %.2f, %.2f)\n", atom_idx, local_idx, label, (int)elem.len, elem.ptr, (int)symbol.len, symbol.ptr, pos.x, pos.y, pos.z);
    len += snprintf(buff + len, ARRAY_SIZE(buff) - 1 - len, "res[%i]: %s %i\n", res_idx, res_name, res_id);
    if (chain_idx) {
        len += snprintf(buff + len, ARRAY_SIZE(buff) - 1 - len, "chain[%i]: %s\n", chain_idx, chain_id);
    }

    /*
    // @TODO: REIMPLEMENT THIS
    if (res_idx < mol.backbone.segment.angleangles.size() && res_idx < mol.backbone.segments.size() && valid_backbone_atoms(mol.backbone.segments[res_idx])) {
        const auto angles = rad_to_deg((vec2)mol.backbone.angles[res_idx]);
        len += snprintf(buff + len, 256 - len, u8"\u03C6: %.1f\u00b0, \u03C8: %.1f\u00b0\n", angles.x, angles.y);
    }
    */

    float x = data.ctx.input.mouse.win_coord.x;
    float y = data.ctx.input.mouse.win_coord.y;

    ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(x + 10.f, y + 18.f) + viewport->Pos);
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
    ImGui::Begin("##Atom Info", 0,
                 ImGuiWindowFlags_Tooltip | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoDocking);
    ImGui::Text("%s", buff);
    ImGui::End();
    ImGui::PopStyleColor();

    /*
ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
ImGui::BeginTooltip();
ImGui::Text("%s", buff);
ImGui::EndTooltip();
ImGui::PopStyleColor();
    */
}

static void draw_molecule_dynamic_info_window(ApplicationData* data) {
    ASSERT(data);
    if (!data->mold.mol.atom.count) return;

    const ImVec2 win_size = {200, 300};
    ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(viewport->Pos + viewport->Size - win_size);
    ImGui::SetNextWindowSize(win_size);
    ImGui::SetNextWindowViewport(viewport->ID);
    ImGui::SetNextWindowBgAlpha(0.0f);

    const ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse |
                                          ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoBringToFrontOnFocus |
                                          ImGuiWindowFlags_NoNavFocus | ImGuiWindowFlags_NoInputs | ImGuiWindowFlags_NoScrollbar;

    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, 0.5f);
    ImGui::Begin("##molecule_dynamic_info", NULL, window_flags);
    {
        const auto& mol = data->mold.mol;
        const auto& traj = data->mold.traj;

        if (mol.atom.count) {
            ImGui::Text("MOL");
            const auto file = extract_file(data->files.molecule);
            if (file.len) ImGui::Text("\"%.*s\"", (int32_t)file.len, file.ptr);
            ImGui::Text("# atoms: %lli", mol.atom.count);
            ImGui::Text("# residues: %lli", mol.residue.count);
            ImGui::Text("# chains: %lli", mol.chain.count);
        }
        if (md_trajectory_num_frames(&traj)) {
            ImGui::NewLine();
            ImGui::Text("TRAJ");
            const auto file = extract_file(data->files.trajectory);
            ImGui::Text("\"%.*s\"", (int32_t)file.len, file.ptr);
            ImGui::Text("# frames: %lli", md_trajectory_num_frames(&traj));
        }
        const int64_t selection_count = md_bitfield_popcount(&data->selection.current_selection_mask);
        const int64_t highlight_count = md_bitfield_popcount(&data->selection.current_highlight_mask);
        if (selection_count || highlight_count) {
            ImGui::NewLine();
            ImGui::Text("SELECTION");
            ImGui::Text("# atoms selected: %lli", selection_count);
            ImGui::Text("# atoms highlighted: %lli", highlight_count);
        }
    }
    ImGui::End();
    ImGui::PopStyleVar(4);
}

static void draw_async_task_window(ApplicationData* data) {
    constexpr float WIDTH = 300.f;
    constexpr float MARGIN = 10.f;
    constexpr float PROGRESSBAR_WIDTH_FRACT = 0.3f;

    //const float stats_fract = stats::fraction_done();

    task_system::ID* tasks = 0;
    const uint32_t num_tasks = task_system::get_tasks(&tasks);

    if (num_tasks) {

        ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(viewport->Pos + ImVec2(data->ctx.window.width - WIDTH - MARGIN,
                                                       ImGui::GetCurrentContext()->FontBaseSize + ImGui::GetStyle().FramePadding.y * 2.f + MARGIN));
        ImGui::SetNextWindowSize(ImVec2(WIDTH, 0));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
        ImGui::Begin("##Async Info", 0,
                     ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar |
                         ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoFocusOnAppearing);

        char buf[32];
        for (uint32_t i = 0; i < MIN(num_tasks, 8); i++) {
            const auto id = tasks[i];
            const char* label = task_system::get_task_label(id);
            const float fract = task_system::get_task_fraction_complete(id);

            if (!label || (label[0] == '#' && label[1] == '#')) continue;

            snprintf(buf, 32, "%.1f%%", fract * 100.f);
            ImGui::ProgressBar(fract, ImVec2(ImGui::GetWindowContentRegionWidth() * PROGRESSBAR_WIDTH_FRACT, 0), buf);
            ImGui::SameLine();
            ImGui::Text("%s", label);
            ImGui::SameLine();
            if (ImGui::Button("X")) {
                task_system::interrupt_task(id);
            }
        }

        ImGui::End();
        ImGui::PopStyleColor();
    }
}


bool draw_property_menu_widgets(DisplayProperty* props, int64_t num_props, uint32_t prop_mask, bool multi_selection = true) {
    bool changed = false;
    int64_t num_candidates = 0;
    int64_t selected_index = -1;

    for (int64_t i = 0; i < num_props; ++i) {
        DisplayProperty& prop = props[i];
        if (prop.possible_display_mask & prop_mask) {
            num_candidates += 1;
            ImPlot::ItemIcon(prop.col); ImGui::SameLine();
            if (ImGui::Selectable(prop.lbl.cstr(), prop.current_display_mask & prop_mask, 0, ImVec2(50,0))) {
                prop.current_display_mask ^= prop_mask;     // Toggle
                changed = true;
                if (prop.current_display_mask & prop_mask) {
                    selected_index = i;
                }
            }
        }
    }

    if (selected_index != -1 && !multi_selection) {
        for (int64_t i = 0; i < num_props; ++i) {
            if (i == selected_index) continue;
            props[i].current_display_mask &= ~prop_mask;    // Clear
        }
    }
    
    if (num_candidates == 0) {
        ImGui::Text("No properties to show.");
    }
    return changed;
    /*
    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    ASSERT(draw_list);
    ImVec2 pos = ImGui::GetWindowPos();

    ImGuiIO& IO = ImGui::GetIO();
    // vars
    const ImVec2 pad = {5,5};
    const ImVec2 spacing = {5,0};
    const float txt_ht      = ImGui::GetTextLineHeight();
    const float icon_size   = txt_ht;
    const float icon_shrink = 2;
    ImU32 col_txt           = ImGui::GetColorU32(ImGuiCol_Text);
    ImU32 col_txt_dis       = ImGui::GetColorU32(ImGuiCol_TextDisabled);
    // render each legend item
    float sum_label_width = 0;
    //bool any_item_hovered = false;
    for (int i = 0; i < num_items; ++i) {
        const char* label       = items[i].lbl.cstr();
        const float label_width = ImGui::CalcTextSize(label, NULL, true).x;
        const ImVec2 top_left   = pos + pad + ImVec2(0, i * (txt_ht + spacing.y));
        sum_label_width        += label_width;
        ImRect icon_bb;
        icon_bb.Min = top_left + ImVec2(icon_shrink,icon_shrink);
        icon_bb.Max = top_left + ImVec2(icon_size - icon_shrink, icon_size - icon_shrink);
        ImRect label_bb;
        label_bb.Min = top_left;
        label_bb.Max = top_left + ImVec2(label_width + icon_size, icon_size);
        ImU32 col_txt_hl;
        ImU32 col_item = ImAlphaU32(items[i].col,1);

        bool icon_hov = false;
        bool icon_hld = false;
        ImGui::PushID(i);
        uint32_t id = ImGui::GetID("btn");
        ImGui::ItemAdd(label_bb, id);
        bool icon_clk = ImGui::ButtonBehavior(icon_bb, id, &icon_hov, &icon_hld);
        ImGui::PopID();
        if (icon_clk)
            items[i].show = !items[i].show;

        if (icon_hov || label_bb.Contains(IO.MousePos)) {
            col_txt_hl = ImMixU32(col_txt, col_item, 64);
        }
        else {
            col_txt_hl = ImGui::GetColorU32(col_txt);
        }
        ImU32 col_icon;
        if (icon_hld)
            col_icon = items[i].show ? ImAlphaU32(col_item,0.5f) : ImGui::GetColorU32(ImGuiCol_TextDisabled, 0.5f);
        else if (icon_hov)
            col_icon = items[i].show ? ImAlphaU32(col_item,0.75f) : ImGui::GetColorU32(ImGuiCol_TextDisabled, 0.75f);
        else
            col_icon = items[i].show ? col_item : col_txt_dis;

        draw_list->AddRectFilled(icon_bb.Min, icon_bb.Max, col_icon, 1);
        const char* text_display_end = ImGui::FindRenderedTextEnd(label, NULL);
        if (label != text_display_end)
            draw_list->AddText(top_left + ImVec2(icon_size, 0), items[i].show ? col_txt_hl  : col_txt_dis, label, text_display_end);
    }
    */
}

struct TimelineArgs {
    const char* lbl;
    uint32_t col;
    int plot_height;

    struct {
        int count;
        const float* x;
        const float* y;
        const float* y_variance;

        float min_y;
        float max_y;
    } values;

    struct {
        double* min;
        double* max;
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

bool draw_property_timeline(const ApplicationData& data, const TimelineArgs& args) {
    const ImPlotAxisFlags axis_flags = 0;
    const ImPlotAxisFlags axis_flags_x = axis_flags;
    const ImPlotAxisFlags axis_flags_y = axis_flags | ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_RangeFit | ImPlotAxisFlags_NoLabel |ImPlotAxisFlags_NoTickLabels;
    const ImPlotFlags flags = ImPlotFlags_AntiAliased;

    ImPlot::LinkNextPlotLimits(args.view_range.min, args.view_range.max, 0, 0);
    if (ImPlot::BeginPlot("##Timeline", NULL, NULL, ImVec2(-1,args.plot_height), flags, axis_flags_x, axis_flags_y)) {
        bool active = ImGui::IsItemActive();

        if (args.value_filter.enabled) {
            float* y_vals = (float*)md_alloc(frame_allocator, args.values.count * sizeof(float));
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
            if (args.values.y_variance) {
                char lbl[32] = "";

                snprintf(lbl, sizeof(lbl), "%s (mean)", args.lbl);
                ImPlot::PlotLine(lbl, args.values.x, args.values.y, args.values.count);

                ImPlot::SetNextFillStyle(line_col, 0.2f);
                snprintf(lbl, sizeof(lbl), "%s (variance)", args.lbl);
                ImPlot::PlotShadedG(lbl,
                    [](void* payload, int idx) -> ImPlotPoint {
                        TimelineArgs* args = (TimelineArgs*)payload;
                        return ImPlotPoint(args->values.x[idx], args->values.y[idx] - args->values.y_variance[idx]);
                    },
                    (void*)&args,
                    [](void* payload, int idx) -> ImPlotPoint {
                        TimelineArgs* args = (TimelineArgs*)payload;
                        return ImPlotPoint(args->values.x[idx], args->values.y[idx] + args->values.y_variance[idx]);
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
            if (active && ImGui::IsMouseDown(0)) {           
                if (ImGui::GetIO().KeyMods == ImGuiKeyModFlags_Shift) {
                    if (args.filter.show && args.filter.enabled) {
                        *args.input.is_selecting = true;
                        *args.filter.beg = ImPlot::GetPlotMousePos().x;
                    }
                } else {
                    *args.input.is_dragging = true;
                }
            } else {
                ImPlotPoint plot_pos = ImPlot::GetPlotMousePos();
                ImVec2 screen_pos = ImPlot::PlotToPixels(plot_pos);
                ImVec2 p0 = {screen_pos.x, ImPlot::GetPlotPos().y};
                ImVec2 p1 = {screen_pos.x, ImPlot::GetPlotPos().y + ImPlot::GetPlotSize().y};
                ImPlot::PushPlotClipRect();
                ImPlot::GetPlotDrawList()->AddLine(p0, p1, IM_COL32(255, 255, 255, 60));
                ImPlot::PopPlotClipRect();

                char buf[128] = "";
                int len = 0;

                double time = plot_pos.x;
                int32_t frame_idx = CLAMP((int)(time_to_frame(time, data) + 0.5), 0, md_array_size(data.timeline.x_values)-1);
                len += snprintf(buf + len, sizeof(buf) - len, "frame: %i, time: %.2f", frame_idx, time);

                if (0 <= frame_idx && frame_idx < args.values.count) {
                    if (args.values.y) {
                        len += snprintf(buf + len, sizeof(buf) - len, ", value: %.2f", args.values.y[frame_idx]);
                    }
                    if (args.values.y_variance) {
                        len += snprintf(buf + len, sizeof(buf) - len, ", variance: %.2f", args.values.y_variance[frame_idx]);
                    }
                }
                ImGui::SetTooltip("%.*s", len, buf);
            }
        }

        if (ImPlot::DragLineX("Current Time", args.time, true, ImVec4(1,1,0,1))) {
            *args.time = CLAMP(*args.time, args.filter.min, args.filter.max);
        }
        ImPlot::EndPlot();
    }

    return true;
}

// #timeline
static void draw_timeline_window(ApplicationData* data) {
    ASSERT(data);
    ImGui::SetNextWindowSize(ImVec2(600, 300), ImGuiCond_FirstUseEver);

    if (ImGui::Begin("Temporal", &data->timeline.show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_MenuBar)) {

        constexpr int MIN_PLOT_HEIGHT = 10;
        constexpr int MAX_PLOT_HEIGHT = 1000;
        static int plot_height = 150;

        double pre_filter_min = data->timeline.filter.beg_frame;
        double pre_filter_max = data->timeline.filter.end_frame;

        const float* x_values = data->timeline.x_values;
        const int num_x_values = md_array_size(data->timeline.x_values);
        const float min_x_value = num_x_values > 0 ? x_values[0] : 0.0f;
        const float max_x_value = num_x_values > 0 ? x_values[num_x_values - 1] : 1.0f;

        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("Properties")) {
                draw_property_menu_widgets(data->display_properties, md_array_size(data->display_properties), DisplayProperty::ShowIn_Timeline);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Filter")) {
                ImGui::Checkbox("Enabled", &data->timeline.filter.enabled);
                if (data->timeline.filter.enabled) {
                    ImGui::Checkbox("Temporal Window", &data->timeline.filter.temporal_window.enabled);
                    if (data->timeline.filter.temporal_window.enabled) {
                        const double extent_min = 1.0;
                        const double extent_max = num_x_values / 2.0;
                        ImGui::SliderScalar("Extent", ImGuiDataType_Double, &data->timeline.filter.temporal_window.extent_in_frames, &extent_min, &extent_max, "%.1f");
                        //ImGui::Slider("Extent", &data->timeline.filter.temporal_window.extent_in_frames, 1, num_x_values / 2);
                    }
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Settings")) {
                ImGui::SliderInt("Plot Height", &plot_height, MIN_PLOT_HEIGHT, MAX_PLOT_HEIGHT);
                ImGui::EndMenu();
            }
            ImGui::EndMenuBar();
        }

        if (ImGui::IsWindowFocused() && ImGui::IsKeyPressed(KEY_PLAY_PAUSE, false)) {
            data->animation.mode = data->animation.mode == PlaybackMode::Playing ? PlaybackMode::Stopped : PlaybackMode::Playing;
        }


        if (num_x_values > 0) {
            ImPlotInputMap old_map = ImPlot::GetInputMap();

            ImPlotInputMap& map = ImPlot::GetInputMap();
            //map.BoxSelectButton = 0;
            //map.BoxSelectCancelButton = 0;
            //map.BoxSelectMod = 0;
            //map.PanMod
            map.PanButton = ImGuiMouseButton_Right;
            map.BoxSelectButton = ImGuiMouseButton_Right;
            map.BoxSelectCancelButton = ImGuiMouseButton_Left;
            map.BoxSelectMod = ImGuiKeyModFlags_Shift;
            map.QueryButton = -1;
            map.QueryMod = -1;
            map.QueryToggleMod = -1;
            map.ContextMenuButton = -1;
            map.FitButton = ImGuiMouseButton_Right;

            static bool is_dragging = false;
            static bool is_selecting = false;

            if (!ImGui::IsMouseDown(0)) {
                is_dragging = false;
                is_selecting = false;
            }

            int64_t num_props = md_array_size(data->display_properties);

            // Create a temporary 'time' representation of the filters min and max value
            // The visualization uses time units while we store 'frame' units
            double filter_beg = frame_to_time(data->timeline.filter.beg_frame, *data);
            double filter_end = frame_to_time(data->timeline.filter.end_frame, *data);
            double time = frame_to_time(data->animation.frame, *data);
            double view_beg = data->timeline.view_range.beg_x;
            double view_end = data->timeline.view_range.end_x;

            TimelineArgs args = {
                .lbl = "##Empty",
                .col = 0xFFFFFFFF,
                .plot_height = plot_height,
                .values = {
                    .count = 0,
                    .x = NULL,
                    .y = NULL,
                    .y_variance = NULL,
                },
                .view_range = {
                    &view_beg,
                    &view_end,
                },
                .input = {
                    .is_dragging = &is_dragging,
                    .is_selecting = &is_selecting,
                },
                .filter = {
                    .show = data->timeline.filter.enabled,
                    .enabled = data->timeline.filter.temporal_window.enabled == false,
                    .beg = &filter_beg,
                    .end = &filter_end,
                    .min = min_x_value,
                    .max = max_x_value,
                },
                .time = &time
            };

            int64_t num_shown_plots = 0;
            for (int64_t i = 0; i < num_props; ++i) {
                DisplayProperty& prop = data->display_properties[i];
                if ((prop.current_display_mask & DisplayProperty::ShowIn_Timeline) == 0) continue;

                ASSERT(prop.full_prop);
                ASSERT(prop.filt_prop);

                const float* y_values = NULL;
                const float* y_variance = NULL;
                int num_y_values = 0;
                    
                if (prop.full_prop->data.aggregate) {
                    y_values = prop.full_prop->data.aggregate->mean;
                    y_variance = prop.full_prop->data.aggregate->variance;
                    num_y_values = prop.full_prop->data.aggregate->num_values;
                } else {
                    y_values = prop.full_prop->data.values;
                    num_y_values = prop.full_prop->data.num_values;
                }
                    
                ASSERT(num_y_values == num_x_values);

                ImGui::PushID(i);
                args.lbl = prop.lbl.cstr();
                args.col = prop.col;
                args.values = {
                    .count = num_y_values,
                    .x = x_values,
                    .y = y_values,
                    .y_variance = y_variance,
                    .min_y = prop.full_prop->data.min_value,
                    .max_y = prop.full_prop->data.max_value
                };
                args.value_filter = {
                    .enabled = data->distributions.filter.enabled,
                    .min = prop.value_filter.beg,
                    .max = prop.value_filter.end
                };
                draw_property_timeline(*data, args);
                ImGui::PopID();

                num_shown_plots += 1;
            }
            if (num_shown_plots == 0) {
                draw_property_timeline(*data, args);
            }

            data->animation.frame = time_to_frame(time, *data);
            data->timeline.filter.beg_frame = time_to_frame(filter_beg, *data);
            data->timeline.filter.end_frame = time_to_frame(filter_end, *data);
            data->timeline.view_range.beg_x = CLAMP(view_beg, min_x_value, max_x_value);
            data->timeline.view_range.end_x = CLAMP(view_end, min_x_value, max_x_value);

            ImPlot::GetInputMap() = old_map;
        }
        
        if (data->timeline.filter.enabled && (data->timeline.filter.beg_frame != pre_filter_min || data->timeline.filter.end_frame != pre_filter_max)) {
            data->mold.script.evaluate_filt = true;
        }
    }
    ImGui::End();
}

// #distribution_window
static void draw_distribution_window(ApplicationData* data) {
    ImGui::SetNextWindowSize(ImVec2(200, 300), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Distributions", &data->distributions.show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_MenuBar)) {
        const int MIN_NUM_BINS = 8;
        const int MAX_NUM_BINS = 256;
        static int num_bins = 64;

        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("Properties")) {
                draw_property_menu_widgets(data->display_properties, md_array_size(data->display_properties), DisplayProperty::ShowIn_Distribution);
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Bins")) {
                if (ImGui::SliderInt("##bins", &num_bins, MIN_NUM_BINS, MAX_NUM_BINS, "%d", ImGuiSliderFlags_Logarithmic)) {
                    const int up   = next_power_of_two32(num_bins);
                    const int down = up / 2;
                    num_bins = abs(num_bins - down) < abs(num_bins - up) ? down : up;
                }
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Filter")) {
                ImGui::Checkbox("Enabled", &data->distributions.filter.enabled);
                ImGui::EndMenu();
            }

            ImGui::EndMenuBar();
        }

        if (ImGui::IsWindowFocused() && ImGui::IsKeyPressed(KEY_PLAY_PAUSE, false)) {
            data->animation.mode = data->animation.mode == PlaybackMode::Playing ? PlaybackMode::Stopped : PlaybackMode::Playing;
        }

        ImPlotAxisFlags axis_flags = 0;
        ImPlotAxisFlags axis_flags_x = axis_flags | 0;
        ImPlotAxisFlags axis_flags_y = axis_flags | ImPlotAxisFlags_NoTickLabels;
        ImPlotFlags flags = ImPlotFlags_AntiAliased;

        // The distribution properties are always computed as histograms with a resolution of 1024
        // If we have a different number of bins for our visualization, we need to recompute the bins
        float* bins = (float*)md_alloc(frame_allocator, num_bins * sizeof(float));
        float* filtered_bins = (float*)md_alloc(frame_allocator, num_bins * sizeof(float));

        int num_props = (int)md_array_size(data->display_properties);
        for (int i = 0; i < num_props; ++i) {
            DisplayProperty& prop = data->display_properties[i];
            if ((prop.current_display_mask & DisplayProperty::ShowIn_Distribution) == 0) continue;

            md_script_property_t& full_prop = *prop.full_prop;
            md_script_property_t& filt_prop = *prop.filt_prop;

            float min_x = full_prop.data.min_range[0];
            float max_x = full_prop.data.max_range[0];
            float min_y = full_prop.data.min_value;
            float max_y = full_prop.data.max_value;
            const float* full_src = 0;
            const float* filt_src = 0;
            int num_values_src = 0;

            if (full_prop.type == MD_SCRIPT_PROPERTY_TYPE_DISTRIBUTION) {
                full_src = full_prop.data.values;
                filt_src = filt_prop.data.values;
                num_values_src = full_prop.data.num_values;
            } else if (full_prop.type == MD_SCRIPT_PROPERTY_TYPE_TEMPORAL) {
                full_src = prop.full_hist.bin;
                filt_src = prop.filt_hist.bin;
                num_values_src = ARRAY_SIZE(prop.full_hist.bin);
            }

            // Downsample bins
            downsample_histogram(bins, num_bins, full_src, num_values_src);
            downsample_histogram(filtered_bins, num_bins, filt_src, num_values_src);

            min_y = 0;
            max_y = 0;
            for (int64_t j = 0; j < num_bins; ++j) {
                max_y = MAX(max_y, bins[j]);
            }
            
            ImGui::PushID(i);
            const double bar_width = (max_x - min_x) / (num_bins-1);
            const double bar_off = min_x;
            const double bar_scl = (max_x - min_x) / (num_bins-1);

            ImPlot::SetNextPlotLimits(min_x, max_x, 0, max_y * 1.1, ImGuiCond_Always);

            if (ImPlot::BeginPlot(prop.lbl.cstr(), 0, 0, ImVec2(-1,150), flags, axis_flags_x, axis_flags_y)) {
                //ImPlot::SetNextFillStyle(vec_cast(qualitative_color_scale(i)), 0.5f);
                //ImPlot::PlotBars(label, draw_bins, num_bins, bar_width, bar_off);

                bool is_active = ImGui::IsItemActive();
                bool is_hovered = ImGui::IsItemHovered();

                struct PlotData {
                    double offset;
                    double scale;
                    const float* bars;
                    int num_bars;
                } plot_data = {
                    .offset = bar_off,
                    .scale = bar_scl,
                    .bars = 0,
                    .num_bars = num_bins
                };

                auto getter = [](void* data, int idx) -> ImPlotPoint {
                    PlotData* pd = (PlotData*)data;
                    return ImPlotPoint(idx * pd->scale + pd->offset, (double)pd->bars[idx]);
                };

                plot_data.bars = bins;
                ImPlot::SetNextFillStyle(ImVec4(0,0,0,-1), 1.0f);
                ImPlot::SetNextLineStyle(ImVec4(0,0,0,0), 0);
                ImPlot::PlotBarsG("full", getter, &plot_data, num_bins, bar_width);

                if (data->timeline.filter.enabled) {
                    plot_data.bars = filtered_bins;
                    ImPlot::SetNextFillStyle(ImVec4(1,1,0,1), 0.3f);
                    ImPlot::SetNextLineStyle(ImVec4(0,0,0,0), 0);
                    ImPlot::PlotBarsG("filt", getter, &plot_data, num_bins, bar_width);
                }
                //ImPlot::PlotHistogram(label, prop.data.values, prop.data.num_values, bins, false, false, range);

                if (data->distributions.filter.enabled) {
                    if ((is_active || is_hovered) && ImGui::IsMouseClicked(0)) {
                        if (ImGui::GetIO().KeyMods == ImGuiKeyModFlags_Shift) {
                            prop.value_filter.beg = ImPlot::GetPlotMousePos().x;
                            prop.value_filter.end = ImPlot::GetPlotMousePos().x;
                        };
                    }

                    double beg = prop.value_filter.beg;
                    double end = prop.value_filter.end;
                    ImPlot::DragRangeX("filter", &beg, &end, min_x, max_x);
                    prop.value_filter.beg = beg;
                    prop.value_filter.end = end;
                }

                ImPlot::EndPlot();
            }
            ImGui::PopID();
        }
    }
    ImGui::End();
}

static void draw_shape_space_window(ApplicationData* data) {
    if (ImGui::Begin("Shape Space", &data->shape_space.show_window, ImGuiWindowFlags_MenuBar)) {
        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("Settings")) {
                const float marker_min_size =  0.01f;
                const float marker_max_size = 10.0f;
                ImGui::SliderFloat("Marker Size", &data->shape_space.marker_size, marker_min_size, marker_max_size);
                ImGui::EndMenu();
            }
            ImGui::EndMenuBar();
        }

        const float item_width = MAX(ImGui::GetWindowContentRegionWidth(), 150.f);
        ImGui::PushItemWidth(item_width);
        if (!data->shape_space.input_valid) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
        data->shape_space.evaluate |= ImGui::InputText("##input", data->shape_space.input.beg(), data->shape_space.input.capacity());
        if (!data->shape_space.input_valid) ImGui::PopStyleColor();
        ImGui::PopItemWidth();
        if (ImGui::IsItemHovered()) {
            if (data->shape_space.input_valid) {
                md_bitfield_clear(&data->selection.current_highlight_mask);
                for (int64_t i = 0; i < data->shape_space.result.num_bitfields; ++i) {
                    md_bitfield_or_inplace(&data->selection.current_highlight_mask, &data->shape_space.result.bitfields[i]);
                }
            } else if (data->shape_space.err_text[0] != '\0') {
                ImGui::SetTooltip("%s", data->shape_space.err_text.cstr());
            }
        }

        ImPlotFlags flags = ImPlotFlags_Equal | ImPlotFlags_AntiAliased | ImPlotFlags_NoMenus | ImPlotFlags_NoMousePos;
        ImPlotAxisFlags axis_flags = ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoTickLabels | ImPlotAxisFlags_NoTickMarks;

        const float x_reset[2] = {-0.025, 1.025};
        const float y_reset[2] = {-0.025, 0.891};
        ImPlot::SetNextPlotLimits(x_reset[0], x_reset[1], y_reset[0], y_reset[1], ImGuiCond_Once);

        if (ImPlot::BeginPlot("##Shape Space Plot", 0, 0, ImVec2(-1,-1), flags, axis_flags, axis_flags)) {
            ImVec2 p0 = ImPlot::PlotToPixels(ImPlotPoint(0.0f, 0.0f));
            ImVec2 p1 = ImPlot::PlotToPixels(ImPlotPoint(1.0f, 0.0f));
            ImVec2 p2 = ImPlot::PlotToPixels(ImPlotPoint(0.5f, 0.86602540378f));
            ImPlot::PushPlotClipRect();
            ImPlot::GetPlotDrawList()->AddTriangleFilled(p0, p1, p2, IM_COL32(255,255,255,20));
            ImPlot::GetPlotDrawList()->AddTriangle(p0, p1, p2, IM_COL32(255,255,255,50));
            ImPlot::PopPlotClipRect();

            ImPlot::PushStyleVar(ImPlotStyleVar_Marker, ImPlotMarker_Square);
            ImPlot::PushStyleColor(ImPlotCol_MarkerFill, ImVec4(0,0,0,0));
            ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, ImVec4(0,0,0,0));
            ImPlot::PlotScatter("##Hidden reset helper", x_reset, y_reset, 2);
            ImPlot::PopStyleColor(2);
            ImPlot::PopStyleVar();

            auto getter = [](void* user_data, int idx) -> ImPlotPoint {
                const vec2_t* coords = (vec2_t*)user_data;
                return {coords[idx].x, coords[idx].y};
            };

            int32_t hovered_structure_idx = -1;
            int32_t hovered_point_idx = -1;
            float   hovered_point_min_dist = FLT_MAX;

            vec2_t mouse_coord = {(float)ImPlot::GetPlotMousePos().x, (float)ImPlot::GetPlotMousePos().y};

            ImPlotLimits lim = ImPlot::GetPlotLimits();
            const float scl = (float)lim.X.Max - (float)lim.X.Min;
            const float MAX_D2 = 0.0001f * scl * scl;

            ImPlot::PushStyleVar(ImPlotStyleVar_MarkerSize, data->shape_space.marker_size);
            ImPlot::PushStyleVar(ImPlotStyleVar_Marker, ImPlotMarker_Square);
            for (int32_t i = 0; i < data->shape_space.num_structures; ++i) {
                int32_t offset = data->shape_space.num_frames * i;
                int32_t count  = data->shape_space.num_frames;
                if (data->timeline.filter.enabled) {
                    offset += data->timeline.filter.beg_frame;
                    count = MAX(data->timeline.filter.end_frame - data->timeline.filter.beg_frame, 0);
                }
                vec2_t* coordinates = data->shape_space.coords + offset;
                char buf[32] = "";
                if (data->shape_space.num_structures == 1) {
                    snprintf(buf, sizeof(buf), "##%i", i+1);
                } else {
                    snprintf(buf, sizeof(buf), "%i", i+1);
                }
                ImPlot::PlotScatterG(buf, getter, coordinates, count);

                auto item = ImPlot::GetItem(buf);
                if (item) {
                    if (item->LegendHovered) {
                        hovered_structure_idx = i;
                    }

                    if (item->Show) {
                        for (int32_t j = offset; j < offset + count; ++j ) {
                            vec2_t delta = mouse_coord - data->shape_space.coords[j];
                            float d2 = vec2_dot(delta, delta);
                            if (d2 < MAX_D2 && d2 < hovered_point_min_dist) {
                                hovered_point_min_dist = d2;
                                hovered_point_idx = j;
                            }
                        }
                    }
                }
            }
            ImPlot::PopStyleVar(2);

            // Redraw hovered index
            ImPlot::PushStyleVar(ImPlotStyleVar_Marker, ImPlotMarker_Square);
            ImPlot::PushStyleVar(ImPlotStyleVar_MarkerSize, data->shape_space.marker_size * 1.1f);
            ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, ImVec4(1,1,1,1));
            if (hovered_structure_idx != -1) {
                int32_t offset = data->shape_space.num_frames * hovered_structure_idx;
                int32_t count  = data->shape_space.num_frames;
                if (data->timeline.filter.enabled) {
                    offset += data->timeline.filter.beg_frame;
                    count = MAX(data->timeline.filter.end_frame - data->timeline.filter.beg_frame, 0);
                }
                vec2_t* coordinates = data->shape_space.coords + offset;
                ImPlot::PlotScatterG("##hovered structure", getter, coordinates, count);
                md_bitfield_copy(&data->selection.current_highlight_mask, &data->shape_space.result.bitfields[hovered_structure_idx]);
                data->mold.dirty_buffers |= MolBit_DirtyFlags;
            }
            if (hovered_point_idx != -1) {
                vec2_t* coordinates = data->shape_space.coords + hovered_point_idx;
                ImPlot::PlotScatterG("##hovered idx", getter, coordinates, 1);
                char buf[128] = "";
                int len = 0;
                int32_t structure_idx = hovered_point_idx / data->shape_space.num_frames;
                int32_t frame_idx = hovered_point_idx % data->shape_space.num_frames;
                vec3_t w = data->shape_space.weights[hovered_point_idx];
                if (data->shape_space.num_structures > 1) {
                    len += snprintf(buf, sizeof(buf), "Structure: %i, ", structure_idx + 1);
                    md_bitfield_copy(&data->selection.current_highlight_mask, &data->shape_space.result.bitfields[structure_idx]);
                    data->mold.dirty_buffers |= MolBit_DirtyFlags;
                }
                len += snprintf(buf + len, sizeof(buf) - len, "Frame: %i, Weight(l,p,s): (%.2f, %.2f, %.2f)", frame_idx, w.x, w.y, w.z);
                ImGui::SetTooltip("%s", buf);

                if (ImGui::IsWindowFocused() && ImPlot::IsPlotHovered() && ImGui::GetIO().MouseReleased[0]) {
                    data->animation.frame = frame_idx;
                }
            }
            ImPlot::PopStyleColor();
            ImPlot::PopStyleVar(2);

            ImPlot::EndPlot();
        }
    }
    ImGui::End();

    if (data->shape_space.evaluate) {
        data->shape_space.input_valid = false;
        const int64_t num_frames = md_trajectory_num_frames(&data->mold.traj);
        if (num_frames > 0) {
            if (data->tasks.evaluate_shape_space.id != 0) {
                task_system::interrupt_and_wait(data->tasks.evaluate_shape_space);
                data->tasks.evaluate_shape_space.id = 0;
            }
            if (md_semaphore_try_aquire(&data->mold.script.ir_semaphore)) {
                data->shape_space.evaluate = false;
                md_array_shrink(data->shape_space.coords, 0);
                md_array_shrink(data->shape_space.weights, 0);

                md_filter_free(&data->shape_space.result, default_allocator);
                if (md_filter_evaluate(&data->shape_space.result, data->shape_space.input, &data->mold.mol, &data->mold.script.ir, default_allocator)) {
                    data->shape_space.input_valid = true;
                    data->shape_space.num_structures = (int32_t)data->shape_space.result.num_bitfields;
                    if (data->shape_space.num_structures > 0) {
                        data->shape_space.num_frames = (int32_t)num_frames;
                        md_array_resize(data->shape_space.coords, num_frames * data->shape_space.num_structures, persistent_allocator);
                        md_array_resize(data->shape_space.weights, num_frames * data->shape_space.num_structures, persistent_allocator);

                        data->tasks.evaluate_shape_space = task_system::enqueue_pool("Eval Shape Space", (uint32_t)num_frames, [data](task_system::TaskSetRange frame_range) mutable {
                            int64_t stride = ROUND_UP(data->mold.mol.atom.count, md_simd_widthf);
                            float* coords = (float*)md_alloc(default_allocator, stride * 3 * sizeof(float));
                            float* x = coords + stride * 0;
                            float* y = coords + stride * 1;
                            float* z = coords + stride * 2;

                            const vec2_t p[3] = {{0.0f, 0.0f}, {1.0f, 0.0f}, {0.5f, 0.86602540378f}};

                            for (uint32_t frame_idx = frame_range.beg; frame_idx < frame_range.end; ++frame_idx) {
                                md_trajectory_load_frame(&data->mold.traj, frame_idx, NULL, x, y, z);
                                vec3_t* xyz = 0;
                                for (int64_t i = 0; i < data->shape_space.result.num_bitfields; ++i) {
                                    md_array_ensure(xyz, md_bitfield_popcount(&data->shape_space.result.bitfields[i]), default_allocator);
                                    int64_t beg_bit = data->shape_space.result.bitfields[i].beg_bit;
                                    int64_t end_bit = data->shape_space.result.bitfields[i].end_bit;
                                    int64_t count = 0;
                                    vec3_t com = {0,0,0};
                                    while ((beg_bit = md_bitfield_scan(&data->shape_space.result.bitfields[i], beg_bit, end_bit)) != 0) {
                                        int64_t src_idx = beg_bit - 1;
                                        vec3_t pos = {x[src_idx], y[src_idx], z[src_idx]};
                                        com = com + pos;
                                        xyz[count++] = pos;
                                    }
                                    com = com / (float)count;
                                    vec3_t eigen_vals;
                                    mat3_t eigen_vecs;
                                    mat3_eigen(mat3_covariance_matrix_vec3(xyz, com, count), eigen_vecs.col, eigen_vals.elem);
                                    float scl = 1.0f / (eigen_vals.x + eigen_vals.y + eigen_vals.z);
                                    int64_t dst_idx = data->shape_space.num_frames * i + frame_idx;
                                    vec3_t w = { (eigen_vals[0] - eigen_vals[1]) * scl, 2.0f * (eigen_vals[1] - eigen_vals[2]) * scl, 3.0f * eigen_vals[2] * scl };
                                    data->shape_space.weights[dst_idx] = w;
                                    data->shape_space.coords[dst_idx] = p[0] * w[0] + p[1] * w[1] + p[2] * w[2];
                                }
                                md_array_free(xyz, default_allocator);
                            }
                        });
                    } else {
                        snprintf(data->shape_space.err_text.beg(), data->shape_space.err_text.capacity(), "Expression did not evaluate into any bitfields");
                    }
                } else {
                    snprintf(data->shape_space.err_text.beg(), data->shape_space.err_text.capacity(), "%s", data->shape_space.result.error_buf);
                }
                md_semaphore_release(&data->mold.script.ir_semaphore);
            }
        } else {
            data->shape_space.evaluate = false;
            snprintf(data->shape_space.err_text.beg(), data->shape_space.err_text.capacity(), "Missing trajectory for evaluating expression");
        }
    }

}

#if 0
static void draw_ramachandran_window(ApplicationData* data) {
    // const int32 num_frames = data->mold.traj ? data->mold.md_trajectory_num_frames(&traj) : 0;
    // const int32 frame = (int32)data->time;
    const Range<int32_t> frame_range = data->timeline.range;
    const auto& mol = data->mold.mol;
    //Array<const BackboneAngle> backbone_angles = get_residue_backbone_angles(mol);
    //Array<const BackboneAtoms> backbone_atoms  = get_residue_backbone_atoms(mol);
    Bitfield& atom_selection = data->selection.current_selection_mask;
    Bitfield& atom_highlight = data->selection.current_highlight_mask;

    ImGui::SetNextWindowSizeConstraints(ImVec2(400, 200), ImVec2(10000, 10000));
    ImGui::Begin("Ramachandran", &data->ramachandran.show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoScrollbar);

    const float w = ImGui::GetContentRegionAvailWidth();

    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    ImGui::BeginChild("canvas", ImVec2(w, w), true, ImGuiWindowFlags_AlwaysHorizontalScrollbar | ImGuiWindowFlags_AlwaysVerticalScrollbar);
    ImGui::PopStyleVar(1);

    static float zoom_factor = 1.0f;
    const float max_c = ImGui::GetContentRegionAvail().x;
    const ImVec2 cont_min = ImGui::GetWindowContentRegionMin() + ImGui::GetWindowPos();
    const ImVec2 size = ImVec2(max_c, max_c) * zoom_factor;

    ImRect bb(ImGui::GetCurrentWindow()->DC.CursorPos, ImGui::GetCurrentWindow()->DC.CursorPos + size);
    ImGui::InvisibleButton("bg", bb.GetSize());

    ImDrawList* dl = ImGui::GetWindowDrawList();

    const ImVec2 mouse_pos = ImGui::GetMousePos();
    dl->ChannelsSplit(4);
    dl->ChannelsSetCurrent(0);
    dl->AddRectFilled(bb.Min, bb.Max, 0xffffffff);
    dl->ChannelsSetCurrent(1);
    dl->AddImage((ImTextureID)(intptr_t)ramachandran::get_gui_texture(), bb.Min, bb.Max);
    if (data->ramachandran.range.enabled) {
        dl->ChannelsSetCurrent(2);
        dl->AddImage((ImTextureID)(intptr_t)ramachandran::get_accumulation_texture(), bb.Min, bb.Max);
    }
    dl->ChannelsSetCurrent(3);

    constexpr float ONE_OVER_TWO_PI = 1.f / (2.f * math::PI);

    int64_t mouse_hover_idx = -1;

    if (data->ramachandran.current.enabled) {
        const uint32_t border_color = convert_color(data->ramachandran.current.border_color);
        const uint32_t base_color = convert_color(data->ramachandran.current.base.fill_color);
        const uint32_t selected_color = convert_color(data->ramachandran.current.selection.selection_color);
        const uint32_t highlight_color = convert_color(data->ramachandran.current.selection.highlight_color);
        const float base_radius = data->ramachandran.current.base.radius;
        const float selected_radius = data->ramachandran.current.selection.radius;

        for (int64_t ci = 0; ci < mol.chain.count; ++ci) {
            const auto range = mol.chain.backbone_range[ci];
            for (int64_t i = range.beg; i < range.end; i++) {
                const auto& angle = mol.backbone.angle[i];
                if (angle.phi == 0.f || angle.psi == 0.f) continue;

                Range<int> atom_range = { mol.residue.atom_range[i].beg, mol.residue.atom_range[i].end };

                const auto selected  = bitfield::any_bit_set_in_range(atom_selection, atom_range);
                const auto highlight = bitfield::any_bit_set_in_range(atom_highlight, atom_range);
                const auto radius = ((highlight || selected) ? selected_radius : base_radius) * zoom_factor;
                const auto fill_color = selected ? selected_color : (highlight ? highlight_color : base_color);

                const ImVec2 coord =
                    ImLerp(bb.Min, bb.Max, ImVec2(angle.phi * ONE_OVER_TWO_PI + 0.5f, -angle.psi * ONE_OVER_TWO_PI + 0.5f));  // [-PI, PI] -> [0, 1]
                const ImVec2 min_box(math::round(coord.x - radius), math::round(coord.y - radius));
                const ImVec2 max_box(math::round(coord.x + radius), math::round(coord.y + radius));
                if (radius > 1.f) {
                    dl->AddRectFilled(min_box, max_box, fill_color);
                    dl->AddRect(min_box, max_box, border_color, 0.f, 15, 1.f);
                } else {
                    dl->AddRectFilled(min_box, max_box, border_color);
                }
                if (min_box.x <= mouse_pos.x && mouse_pos.x <= max_box.x && min_box.y <= mouse_pos.y && mouse_pos.y <= max_box.y) {
                    mouse_hover_idx = i;
                }
            }
        }
    }

    const auto cx = math::round(math::mix(bb.Min.x, bb.Max.x, 0.5f));
    const auto cy = math::round(math::mix(bb.Min.y, bb.Max.y, 0.5f));
    dl->AddLine(ImVec2(cx, bb.Min.y), ImVec2(cx, bb.Max.y), 0xff000000, 0.5f);
    dl->AddLine(ImVec2(bb.Min.x, cy), ImVec2(bb.Max.x, cy), 0xff000000, 0.5f);
    dl->ChannelsMerge();
    dl->ChannelsSetCurrent(0);

    enum class Mode { Append, Remove };
    static ImVec2 region_x0 = {0, 0};
    static bool region_select = false;
    static Mode region_mode = Mode::Append;
    const ImVec2 region_x1 = ImGui::GetMousePos();

    const bool mouse_hit     = (data->ctx.input.mouse.hit[0] || data->ctx.input.mouse.hit[1]);
    const bool mouse_down    = (data->ctx.input.mouse.down[0] || data->ctx.input.mouse.down[1]);
    const bool mouse_clicked = (data->ctx.input.mouse.clicked[0] || data->ctx.input.mouse.clicked[1]);

    const bool shift_down = ImGui::GetIO().KeyShift;

    if (ImGui::IsItemHovered()) {
        clear_highlight(data);
        //bitfield::clear_all(data->selection.current_highlight_mask);
        //data->mold.dirty_buffers |= MolBit_DirtyFlags;

        const ImVec2 normalized_coord = ((ImGui::GetMousePos() - bb.Min) / (bb.Max - bb.Min) - ImVec2(0.5f, 0.5f)) * ImVec2(1, -1);
        const ImVec2 angles = normalized_coord * 2.f * 180.f;
        ImGui::BeginTooltip();
        ImGui::Text((const char*)u8"\u03C6: %.1f\u00b0, \u03C8: %.1f\u00b0", angles.x, angles.y);
        if (!region_select && mouse_hover_idx != -1) {
            const ResIdx res_idx = (ResIdx)mouse_hover_idx;
            const Range<int> atom_range = { mol.residue.atom_range[res_idx].beg, mol.residue.atom_range[res_idx].end };
            ImGui::Text("Residue[%i]: %s", res_idx, mol.residue.name[res_idx]);
            modify_highlight(data, atom_range, SelectionOperator::Replace);
        }
        ImGui::EndTooltip();

        if (shift_down && mouse_hit) {
            region_x0 = ImGui::GetIO().MousePos;
            region_mode = data->ctx.input.mouse.hit[0] ? Mode::Append : Mode::Remove;
        }

        if (shift_down && mouse_down && region_x1 != region_x0) {
            region_select = true;
        }

        if (shift_down && !region_select && mouse_clicked) {
            if (mouse_hover_idx != -1) {
                const ResIdx res_idx = (ResIdx)mouse_hover_idx;
                const Range<int> atom_range = { mol.residue.atom_range[res_idx].beg, mol.residue.atom_range[res_idx].end };
                const bool append = data->ctx.input.mouse.clicked[0];
                SelectionOperator op = append ? SelectionOperator::Or : SelectionOperator::Clear;                
                modify_selection(data, atom_range, op);
            } else if (data->ctx.input.mouse.clicked[1]) {
                bitfield::clear_all(data->selection.current_highlight_mask);
                bitfield::clear_all(data->selection.current_selection_mask);
            }
        }
    }

    if (region_select) {
        const ImVec2 x0 = ImMin(region_x0, region_x1);
        const ImVec2 x1 = ImMax(region_x0, region_x1);
        const ImU32 fill_col = 0x22222222;
        const ImU32 line_col = 0x88888888;
        dl->AddRectFilled(x0, x1, fill_col);
        dl->AddRect(x0, x1, line_col);

        //bitfield::clear_all(data->selection.current_highlight_mask);
        //clear_highlight(data);

        for (int64_t i = 0; i < mol.backbone.count; ++i) {
            const auto& angle = mol.backbone.angle[i];
            if (angle.phi == 0.f || angle.psi == 0.f) continue;
            const ImVec2 coord =
                ImLerp(bb.Min, bb.Max, ImVec2(angle.phi * ONE_OVER_TWO_PI + 0.5f, -angle.psi * ONE_OVER_TWO_PI + 0.5f));  // [-PI, PI] -> [0, 1]
            if (coord.x < x0.x || x1.x < coord.x) continue;
            if (coord.y < x0.y || x1.y < coord.y) continue;
            
            const Range<int> atom_range = { mol.residue.atom_range[i].beg, mol.residue.atom_range[i].end };
            modify_highlight(data, atom_range, SelectionOperator::Or);
            //bitfield::set_range(data->selection.current_highlight_mask, mol.residue.atom_range[i]);
        }
        
        if (region_mode == Mode::Remove) {
            bitfield::and_not_field(data->selection.current_highlight_mask, data->selection.current_selection_mask, data->selection.current_highlight_mask);
        }

        if (!shift_down || !(data->ctx.input.mouse.down[0] || data->ctx.input.mouse.down[1])) {
            // Commit range to selection
            if (region_mode == Mode::Append) {
                //bitfield::or_field(data->selection.current_selection_mask, data->selection.current_selection_mask, data->selection.current_highlight_mask);
                modify_selection(data, data->selection.current_highlight_mask, SelectionOperator::Or);
            } else {
                modify_selection(data, data->selection.current_highlight_mask, SelectionOperator::Replace);
                //bitfield::copy(data->selection.current_selection_mask, data->selection.current_highlight_mask);
            }
            //bitfield::clear_all(data->selection.current_highlight_mask);)
            region_select = false;
        }
        //data->mold.dirty_buffers |= MolBit_DirtyFlags;
    }

    if (ImGui::IsItemHovered()) {
        // ZOOM
        if (ImGui::GetIO().KeyCtrl && ImGui::GetIO().MouseWheel != 0.f) {
            const ImVec2 pos = (ImGui::GetIO().MousePos - cont_min);
            const float mouse_wheel_delta = ImGui::GetIO().MouseWheel;
            if (ImGui::GetIO().KeyCtrl && mouse_wheel_delta != 0.f) {
                constexpr float ZOOM_SCL = 0.1f;
                const float old_zoom = zoom_factor;
                const float new_zoom = CLAMP(zoom_factor + zoom_factor * ZOOM_SCL * mouse_wheel_delta, 1.f, 10.f);
                const ImVec2 delta = pos * (new_zoom / old_zoom) - pos;
                ImGui::SetScrollX(ImGui::GetScrollX() + delta.x);
                ImGui::SetScrollY(ImGui::GetScrollY() + delta.y);
                zoom_factor = new_zoom;
            }
        }

        // PAN
        if (!shift_down && zoom_factor > 1.0f && ImGui::GetIO().MouseDown[1] && ImGui::GetIO().MouseDelta != ImVec2(0,0)) {
            ImGui::SetScrollX(ImGui::GetScrollX() - ImGui::GetIO().MouseDelta.x);
            ImGui::SetScrollY(ImGui::GetScrollY() - ImGui::GetIO().MouseDelta.y);
        }
    }

    ImGui::EndChild();

    ImGui::End();
}
#endif

static void draw_density_volume_window(ApplicationData* data) {
    ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Density Volume", &data->density_volume.show_window, ImGuiWindowFlags_MenuBar)) {
        const ImVec2 button_size = {160, 20};
        bool property_selection_changed = false;

        if (ImGui::IsWindowFocused() && ImGui::IsKeyPressed(KEY_PLAY_PAUSE, false)) {
            data->animation.mode = data->animation.mode == PlaybackMode::Playing ? PlaybackMode::Stopped : PlaybackMode::Playing;
        }

        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("Property")) {
                property_selection_changed = draw_property_menu_widgets(data->display_properties, md_array_size(data->display_properties), DisplayProperty::ShowIn_Volume, false);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("DVR")) {
                ImGui::Checkbox("Enabled", &data->density_volume.dvr.enabled);
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
                ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0, 0, 0, 1.0f));
                ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
                ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(1.0f, 1.0f));
                if (ImGui::ImageButton((void*)(intptr_t)data->density_volume.dvr.tf.id, button_size)) {
                    auto res = application::file_dialog(application::FileDialogFlags_Open, {}, make_cstr("png,jpg"));
                    if (res.result == application::FileDialogResult::Ok) {
                        data->density_volume.dvr.tf.path = res.path;
                        data->density_volume.dvr.tf.dirty = true;
                    }
                }
                ImGui::PopStyleVar(2);
                ImGui::PopStyleColor(2);

                ImGui::SliderFloat("Density Scaling", &data->density_volume.dvr.density_scale, 0.001f, 10000.f, "%.3f", ImGuiSliderFlags_Logarithmic);
                ImGui::SliderFloat("Alpha Scaling", &data->density_volume.dvr.tf.alpha_scale, 0.001f, 10.f, "%.3f", ImGuiSliderFlags_Logarithmic);
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("ISO")) {
                ImGui::Checkbox("Enabled", &data->density_volume.iso.enabled);
                for (int i = 0; i < data->density_volume.iso.isosurfaces.count; ++i) {
                    ImGui::PushID(i);
                    ImGui::SliderFloat("##Isovalue", &data->density_volume.iso.isosurfaces.values[i], 0.0f, 10.f, "%.3f", ImGuiSliderFlags_Logarithmic);
                    if (ImGui::IsItemDeactivatedAfterEdit()) {
                        sort(data->density_volume.iso.isosurfaces);
                    }
                    ImGui::SameLine();
                    vec4_t col = convert_color(data->density_volume.iso.isosurfaces.colors[i]);
                    if (ImGui::ColorEdit4Minimal("##Color", col.elem)) {
                        data->density_volume.iso.isosurfaces.colors[i] = convert_color(col);
                    }
                    ImGui::PopID();
                }
                if ((data->density_volume.iso.isosurfaces.count < data->density_volume.iso.isosurfaces.MaxCount) &&
                    ImGui::Button("Add", button_size)) {
                    insert(data->density_volume.iso.isosurfaces, 0.1f, ImColor(0.2f, 0.1f, 0.9f, 1.0f));
                    sort(data->density_volume.iso.isosurfaces);
                }
                if (ImGui::Button("Clear", button_size)) {
                    clear(data->density_volume.iso.isosurfaces);
                }
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Clip planes")) {
                ImGui::RangeSliderFloat("x", &data->density_volume.clip_volume.min.x, &data->density_volume.clip_volume.max.x, 0.0f, 1.0f);
                ImGui::RangeSliderFloat("y", &data->density_volume.clip_volume.min.y, &data->density_volume.clip_volume.max.y, 0.0f, 1.0f);
                ImGui::RangeSliderFloat("z", &data->density_volume.clip_volume.min.z, &data->density_volume.clip_volume.max.z, 0.0f, 1.0f);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Visualize")) {
                ImGui::Checkbox("Bounding Box", &data->density_volume.show_bounding_box);
                if (data->density_volume.show_bounding_box) {
                    ImGui::SameLine();
                    ImGui::ColorEdit4Minimal("Color", data->density_volume.bounding_box_color.elem);
                }
                ImGui::Checkbox("Reference Structure", &data->density_volume.show_reference_structures);
                if (data->density_volume.show_reference_structures) {
                    ImGui::Checkbox("Show Superimposed Structures", &data->density_volume.show_reference_ensemble);
                }
                ImGui::EndMenu();
            }

            ImGui::EndMenuBar();
        }

        if (data->density_volume.dvr.tf.dirty) {
            image_t img = {0};
            if (read_image(&img, data->density_volume.dvr.tf.path.cstr(), frame_allocator)) {
                gl::init_texture_2D(&data->density_volume.dvr.tf.id, img.width, img.height, GL_RGBA8);
                gl::set_texture_2D_data(data->density_volume.dvr.tf.id, img.data, GL_RGBA8);
                data->density_volume.dvr.tf.dirty = false;
            }
        }

        bool update_volume = property_selection_changed;
        bool update_representations = property_selection_changed;
        md_script_property_t* prop = 0;
        uint64_t fingerprint = 0;

        for (int64_t i = 0; i < md_array_size(data->display_properties); ++i) {
            if (data->display_properties[i].current_display_mask & DisplayProperty::ShowIn_Volume) {
                if (data->timeline.filter.enabled) {
                    prop = data->display_properties[i].filt_prop;
                    fingerprint = data->display_properties[i].filt_prop->data.fingerprint;
                } else {
                    prop = data->display_properties[i].full_prop;
                    fingerprint = data->display_properties[i].full_prop->data.fingerprint;
                }
            }
        }

        static uint64_t s_fingerprint = 0;
        if (s_fingerprint != fingerprint) {
            s_fingerprint = fingerprint;
            update_volume = true;
            update_representations = true;
        }

        static double s_frame = 0;
        if (s_frame != data->animation.frame) {
            s_frame = data->animation.frame;
            update_representations = true;
        }

        // Canvas
        // Using InvisibleButton() as a convenience 1) it will advance the layout cursor and 2) allows us to use IsItemHovered()/IsItemActive()
        ImVec2 canvas_p0 = ImGui::GetCursorScreenPos();      // ImDrawList API uses screen coordinates!
        ImVec2 canvas_sz = ImGui::GetContentRegionAvail();   // Resize canvas to what's available
        if (canvas_sz.x < 50.0f) canvas_sz.x = 50.0f;
        if (canvas_sz.y < 50.0f) canvas_sz.y = 50.0f;
        ImVec2 canvas_p1 = ImVec2(canvas_p0.x + canvas_sz.x, canvas_p0.y + canvas_sz.y);
        
        // Draw border and background color
        ImGuiIO& io = ImGui::GetIO();
        ImDrawList* draw_list = ImGui::GetWindowDrawList();
        draw_list->AddImage((ImTextureID)(uint64_t)data->density_volume.fbo.deferred.post_tonemap, canvas_p0, canvas_p1, {0,1}, {1,0});
        draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));

        // This will catch our interactions
        ImGui::InvisibleButton("canvas", canvas_sz, ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);
        const bool is_hovered = ImGui::IsItemHovered(); // Hovered
        const bool is_active = ImGui::IsItemActive();   // Held
        const ImVec2 origin(canvas_p0.x, canvas_p0.y);  // Lock scrolled origin
        const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

        auto& gbuf = data->density_volume.fbo;
        if (gbuf.width != canvas_sz.x || gbuf.height != canvas_sz.y) {
            init_gbuffer(&gbuf, canvas_sz.x, canvas_sz.y);
        }

        static bool init = true;
        if (init) {
            const vec3_t ext = {10,10,10};
            const vec3_t cen = {0.5f,0.5f,0.5f};
            const vec3_t pos = cen + ext * vec3_t{1.0f,0.5f,0.7f} * 2.f;
            const vec3_t up = {0,1,0};

            data->density_volume.camera.position = pos;
            data->density_volume.camera.focus_distance = vec3_length(pos - cen);
            data->density_volume.camera.orientation = quat_from_mat4(look_at(data->density_volume.camera.position, cen, up));
            init = false;
        }

        static TrackballControllerParam param = {
            .min_distance = 1.0,
            .max_distance = 1000.0,
        };
        vec2_t delta = {data->ctx.input.mouse.win_delta.x, data->ctx.input.mouse.win_delta.y};
        vec2_t curr = {mouse_pos_in_canvas.x, mouse_pos_in_canvas.y};
        vec2_t prev = curr - delta;
        TrackballControllerInput input = {
            .rotate_button = is_active && data->ctx.input.mouse.down[0],
            .pan_button = is_active && data->ctx.input.mouse.down[1],
            .dolly_button = is_active && data->ctx.input.mouse.down[2],
            .dolly_delta = is_hovered ? data->ctx.input.mouse.scroll_delta : 0.0f,
            .mouse_coord_prev = prev,
            .mouse_coord_curr = curr,
            .screen_size = {canvas_sz.x, canvas_sz.y},
            .fov_y = data->density_volume.camera.fov_y,
        };
        camera_controller_trackball(&data->density_volume.camera.position, &data->density_volume.camera.orientation, &data->density_volume.camera.focus_distance, input, param);

        mat4_t view_mat = camera_world_to_view_matrix(data->density_volume.camera);
        mat4_t proj_mat = camera_perspective_projection_matrix(data->density_volume.camera, (int)canvas_sz.x, (int)canvas_sz.y);

        clear_gbuffer(&gbuf);
        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);

        if (update_volume) {
            data->density_volume.model_mat = volume::compute_model_to_world_matrix({0,0,0}, {1,1,1});
            if (prop) {
                if (!data->density_volume.volume_texture.id) {
                    gl::init_texture_3D(&data->density_volume.volume_texture.id, prop->data.dim[0], prop->data.dim[1], prop->data.dim[2], GL_R32F);
                    data->density_volume.volume_texture.dim_x = prop->data.dim[0];
                    data->density_volume.volume_texture.dim_y = prop->data.dim[1];
                    data->density_volume.volume_texture.dim_z = prop->data.dim[2];
                    data->density_volume.volume_texture.max_value = prop->data.max_value;
                }
                gl::set_texture_3D_data(data->density_volume.volume_texture.id, prop->data.values, GL_R32F);
            }
        }

        if (update_representations) {
            if (prop) {
                int64_t num_reps = 0;

                md_semaphore_aquire(&data->mold.script.ir_semaphore);
                md_script_visualization_t vis = {};
                md_script_visualization_args_t args = {
                    .token = prop->vis_token,
                    .ir = &data->mold.script.ir,
                    .mol = &data->mold.mol,
                    .traj = &data->mold.traj,
                    .alloc = frame_allocator,
                    .flags = MD_SCRIPT_VISUALIZE_SDF
                };
                if (md_script_visualization_init(&vis, args)) {
                    if (vis.sdf.extent) {
                        const float s = vis.sdf.extent;
                        vec3_t min_aabb = {-s, -s, -s};
                        vec3_t max_aabb = {s, s, s};
                        data->density_volume.model_mat = volume::compute_model_to_world_matrix(min_aabb, max_aabb);
                    }
                    num_reps = vis.sdf.count;
                }
                md_semaphore_release(&data->mold.script.ir_semaphore);

                if (data->density_volume.gl_reps) {
                    // Only free those required
                    for (int64_t i = num_reps; i < md_array_size(data->density_volume.gl_reps); ++i) {
                        md_gl_representation_free(&data->density_volume.gl_reps[i]);
                    }
                }
                md_array_resize(data->density_volume.gl_reps, num_reps, persistent_allocator);
                md_array_resize(data->density_volume.rep_model_mats, num_reps, persistent_allocator);

                const int64_t num_colors = data->mold.mol.atom.count;
                uint32_t* colors = (uint32_t*)md_alloc(frame_allocator, sizeof(uint32_t) * num_colors);
                for (int64_t i = 0; i < num_reps; ++i) {
                    color_atoms_cpk(colors, num_colors, data->mold.mol);
                    filter_colors(colors, num_colors, &vis.sdf.structures[i]);
                    md_gl_representation_init(&data->density_volume.gl_reps[i], &data->mold.gl_mol);
                    md_gl_representation_set_color(&data->density_volume.gl_reps[i], 0, num_colors, colors, 0);
                    md_gl_representation_set_type_and_args(&data->density_volume.gl_reps[i], MD_GL_REP_SPACE_FILL, {
                        .space_fill = {
                            .radius_scale = 1.0f
                        }
                    });
                    data->density_volume.rep_model_mats[i] = vis.sdf.matrices[i];
                }
            }
        }

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.deferred.fbo);
        glDrawBuffer(GL_COLOR_ATTACHMENT_POST_TONEMAP);
        glViewport(0, 0, gbuf.width, gbuf.height);
        glClearColor(1,1,1,1);
        glClear(GL_COLOR_BUFFER_BIT);

        int64_t num_reps = md_array_size(data->density_volume.gl_reps);
        if (prop && data->density_volume.show_reference_structures && num_reps > 0) {
            num_reps = data->density_volume.show_reference_ensemble ? num_reps : 1;

            md_gl_rendertarget_t render_target = {
                .width = (uint32_t)canvas_sz.x,
                .height = (uint32_t)canvas_sz.y,
                .texture_depth = gbuf.deferred.depth,
                .texture_color = gbuf.deferred.color,
                .texture_atom_index  = gbuf.deferred.picking,
                .texture_view_normal = gbuf.deferred.normal,
                .texture_view_velocity = gbuf.deferred.velocity,
            };

            const md_gl_representation_t** reps = (const md_gl_representation_t**)md_alloc(frame_allocator, num_reps * sizeof(void*));
            const float** mats = (const float**)md_alloc(frame_allocator, num_reps * sizeof(void*));

            for (int64_t i = 0; i < num_reps; ++i) {
                reps[i] = &data->density_volume.gl_reps[i];
                mats[i] = &data->density_volume.rep_model_mats[i].elem[0][0];
            }

            md_gl_draw_args_t draw_args = {
                .representation = {
                    .count = (uint32_t)num_reps,
                    .data = reps,
                    .model_matrix = mats,
                },
                .view_transform = {
                    .view_matrix = &view_mat.elem[0][0],
                    .projection_matrix = &proj_mat.elem[0][0],
                },
                .render_target = &render_target,
            };

            md_gl_draw(&data->mold.gl_ctx, &draw_args);

            if (is_hovered) {
                vec2_t coord = {mouse_pos_in_canvas.x, (float)gbuf.height - mouse_pos_in_canvas.y};
                PickingData pd = read_picking_data(&gbuf, (int)coord.x, (int)coord.y);
                if (pd.idx != INVALID_PICKING_IDX) {
                    draw_atom_info_window(*data, pd.idx);
                }
            }

            PUSH_GPU_SECTION("Postprocessing")
                postprocessing::Descriptor postprocess_desc = {
                .background = {
                        .intensity = data->visuals.background.color * data->visuals.background.intensity,
            },
            .bloom = {
                        .enabled = false,
            },
            .tonemapping = {
                        .enabled = data->visuals.tonemapping.enabled,
                        .mode = data->visuals.tonemapping.tonemapper,
                        .exposure = data->visuals.tonemapping.exposure,
                        .gamma = data->visuals.tonemapping.gamma,
            },
            .ambient_occlusion = {
                        .enabled = data->visuals.ssao.enabled,
                        .radius = data->visuals.ssao.radius,
                        .intensity = data->visuals.ssao.intensity,
                        .bias = data->visuals.ssao.bias,
            },
            .depth_of_field = {
                        .enabled = data->visuals.dof.enabled,
                        .focus_depth = data->visuals.dof.focus_depth.current,
                        .focus_scale = data->visuals.dof.focus_scale,
            },
            .temporal_reprojection = {
                        .enabled = false,
            },
            .input_textures = {
                        .depth = gbuf.deferred.depth,
                        .color = gbuf.deferred.color,
                        .normal = gbuf.deferred.normal,
                        .velocity = gbuf.deferred.velocity,
            }
            };

            ViewParam view_param = {
                .matrix = {
                    .current = {
                    .view = view_mat,
                    .proj = proj_mat,
                    .proj_jittered = proj_mat,
                    .view_proj = mat4_mul(proj_mat, view_mat),
                    .view_proj_jittered = mat4_mul(proj_mat, view_mat),
                    .norm = view_mat,
            },
            .inverse = {
                    .proj = mat4_inverse(proj_mat),
                    .proj_jittered = mat4_inverse(proj_mat),
            }
            },
                .clip_planes = {
                    .near = data->density_volume.camera.near_plane,
                    .far = data->density_volume.camera.far_plane,
            },
            .fov_y = data->density_volume.camera.fov_y,
            .resolution = {canvas_sz.x, canvas_sz.y}
            };

            postprocessing::shade_and_postprocess(postprocess_desc, view_param);
            POP_GPU_SECTION()
        }

        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);

        if (data->density_volume.show_bounding_box) {
            const vec3_t min_box = {0,0,0};
            const vec3_t max_box = {1,1,1};

            immediate::set_model_view_matrix(mat4_mul(view_mat, data->density_volume.model_mat));
            immediate::set_proj_matrix(proj_mat);

            uint32_t box_color = convert_color(data->density_volume.bounding_box_color);
            uint32_t clip_color = convert_color(data->density_volume.clip_volume_color);
            immediate::draw_box_wireframe(min_box, max_box, box_color);
            immediate::draw_box_wireframe(data->density_volume.clip_volume.min, data->density_volume.clip_volume.max, clip_color);

            immediate::flush();
        }

        if (prop) {
            volume::RenderDesc vol_desc = {
                .render_target = {
                    .texture = gbuf.deferred.post_tonemap,
                    .width   = gbuf.width,
                    .height  = gbuf.height,
                },
                .texture = {
                    .volume = data->density_volume.volume_texture.id,
                    .transfer_function = data->density_volume.dvr.tf.id,
                    .depth = gbuf.deferred.depth,
                },
                .matrix = {
                    .model = data->density_volume.model_mat,
                    .view = view_mat,
                    .proj = proj_mat
                },
                .clip_volume = {
                    .min = data->density_volume.clip_volume.min,
                    .max = data->density_volume.clip_volume.max
                },
                .bounding_box = {
                    .color = data->density_volume.bounding_box_color,
                    .enabled = data->density_volume.show_bounding_box,
                },
                .global_scaling = {
                    .density = data->density_volume.dvr.density_scale,
                    .alpha = data->density_volume.dvr.tf.alpha_scale
                },
                .isosurface = data->density_volume.iso.isosurfaces,
                .isosurface_enabled = data->density_volume.iso.enabled,
                .direct_volume_rendering_enabled = data->density_volume.dvr.enabled,

                .voxel_spacing = data->density_volume.voxel_spacing
            };
            volume::render_volume(vol_desc);
        }

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glDrawBuffer(GL_BACK);
    }

    ImGui::End();
}

static void draw_dataset_window(ApplicationData* data) {
    ASSERT(data);

    ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Dataset", &data->dataset.show_window)) {
        str_t mol_file = data->files.molecule;
        ImGui::Text("Molecular data: %.*s", (int)mol_file.len, mol_file.ptr);
        ImGui::Text("Num atoms:    %9lli", data->mold.mol.atom.count);
        ImGui::Text("Num residues: %9lli", data->mold.mol.residue.count);
        ImGui::Text("Num chains:   %9lli", data->mold.mol.chain.count);

        str_t traj_file = data->files.trajectory;
        if (traj_file.len) {
            ImGui::Separator();
            ImGui::Text("Trajectory data: %.*s", (int)traj_file.len, traj_file.ptr);
            ImGui::Text("Num frames:    %9lli", md_trajectory_num_frames(&data->mold.traj));
            ImGui::Text("Num atoms:     %9lli", md_trajectory_num_atoms(&data->mold.traj));
        }

        int64_t num_mappings = md_array_size(data->dataset.atom_element_remappings);
        if (num_mappings) {
            ImGui::Separator();
            ImGui::Text("Atom element mappings, label -> element");
            for (int64_t i = 0; i < num_mappings; ++i) {
                const auto& mapping = data->dataset.atom_element_remappings[i];
                ImGui::Text("%s -> %s (%s)", mapping.lbl.cstr(), md_util_element_name(mapping.elem).ptr, md_util_element_symbol(mapping.elem).ptr);
            }
        }
    }
    ImGui::End();
}

static void draw_debug_window(ApplicationData* data) {
    ASSERT(data);

    ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Debug", &data->show_debug_window)) {
        int32_t sema_count = 0;
        if (md_semaphore_query_count(&data->mold.script.ir_semaphore, &sema_count)) {
            ImGui::Text("Script IR semaphore count: %i", sema_count);
        }
    }
    ImGui::End();
}


static void draw_script_editor_window(ApplicationData* data) {
    ASSERT(data);

    TextEditor& editor = data->script.editor;

    if (ImGui::Begin("Script Editor", &data->script.show_editor, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_MenuBar)) {
        ImGui::SetWindowSize(ImVec2(800, 600), ImGuiCond_FirstUseEver);
        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("File")) {
                if (ImGui::MenuItem("Load")) {
                    application::FileDialogResult file_result = application::file_dialog(application::FileDialogFlags_Open, {}, make_cstr("txt"));
                    if (file_result.result == application::FileDialogResult::Ok) {
                        str_t txt = load_textfile({file_result.path, file_result.path_len}, default_allocator);
                        defer { free_str(txt, default_allocator); };
                        std::string str(txt.ptr, txt.len);
                        editor.SetText(str);
                    }
                }
                if (ImGui::MenuItem("Save")) {
                    auto textToSave = editor.GetText();
                    application::FileDialogResult file_result = application::file_dialog(application::FileDialogFlags_Open, {}, make_cstr("txt"));
                    if (file_result.result == application::FileDialogResult::Ok) {
                        StrBuf<1024> path = str_t{file_result.path, file_result.path_len};
                        if (str_empty(extract_ext(path))) {
                            path += ".txt";
                        }
                        md_file_o* file = md_file_open(path, MD_FILE_WRITE);
                        if (file) {
                            md_file_write(file, textToSave.c_str(), textToSave.length());
                            md_file_close(file);
                        } else {
                            md_printf(MD_LOG_TYPE_ERROR, "Failed to open file '%s' for saving script", path.cstr());
                        }
                    }
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

            if (ImGui::BeginMenu("View")) {
                if (ImGui::MenuItem("Dark palette"))
                    editor.SetPalette(TextEditor::GetDarkPalette());
                if (ImGui::MenuItem("Light palette"))
                    editor.SetPalette(TextEditor::GetLightPalette());
                if (ImGui::MenuItem("Retro blue palette"))
                    editor.SetPalette(TextEditor::GetRetroBluePalette());
                ImGui::EndMenu();
            }
            if (ImGui::MenuItem("Export")) {
                data->show_property_export_window = true;
            }

            ImGui::EndMenuBar();
        }

        if (editor.IsTextChanged()) {
            data->mold.script.compile_ir = true;
            data->mold.script.time_since_last_change = 0;
        }

        editor.Render("TextEditor");

        const TextEditor::Marker* hovered_marker = editor.GetHoveredMarker();
        if (hovered_marker) {
            md_script_visualization_t vis = {0};
            md_script_visualization_init(&vis, {
                .token = (const md_script_vis_token_t*)hovered_marker->payload,
                .ir = &data->mold.script.ir,
                .mol = &data->mold.mol,
                .traj = &data->mold.traj,
                .alloc = frame_allocator,
            });

            immediate::set_model_view_matrix(data->view.param.matrix.current.view);
            immediate::set_proj_matrix(data->view.param.matrix.current.proj_jittered);

            const vec3_t* vertices = (const vec3_t*)vis.vertex.pos;
            glDisable(GL_CULL_FACE);

            for (int64_t i = 0; i < vis.point.count; ++i) {
                ASSERT(vis.point.idx);
                uint16_t idx = vis.point.idx[i];
                immediate::draw_point(vertices[idx]);
            }

            for (int64_t i = 0; i < vis.line.count; ++i) {
                ASSERT(vis.line.idx);
                uint16_t idx[2] = { vis.line.idx[i * 2 + 0], vis.line.idx[i * 2 + 1] };
                immediate::draw_line(vertices[idx[0]], vertices[idx[1]], immediate::COLOR_CYAN);
            }

            for (int64_t i = 0; i < vis.triangle.count; ++i) {
                ASSERT(vis.triangle.idx);
                uint16_t idx[3] = { vis.triangle.idx[i * 3 + 0], vis.triangle.idx[i * 3 + 1], vis.triangle.idx[i * 3 + 2] };
                immediate::draw_triangle(vertices[idx[0]], vertices[idx[1]], vertices[idx[2]], 0x3300FFFF);
            }

            immediate::flush();

            glEnable(GL_CULL_FACE);

            if (!md_bitfield_empty(vis.atom_mask)) {
                md_bitfield_copy(&data->selection.current_highlight_mask, vis.atom_mask);
                data->mold.dirty_buffers |= MolBit_DirtyFlags;
            }
            md_script_visualization_free(&vis);
        }

    }
    ImGui::End();
}

static bool export_xvg(const float* column_data[], const char* column_labels[], int num_columns, int num_rows, str_t filename) {
    ASSERT(column_data);
    ASSERT(column_labels);
    ASSERT(num_columns >= 0);
    ASSERT(num_rows >= 0);

    md_file_o* file = md_file_open(filename, MD_FILE_WRITE);
    if (!file) {
        md_printf(MD_LOG_TYPE_ERROR, "Failed to open file '%.*s' to write data.", (int)filename.len, filename.ptr);
        return false;
    }    

    time_t t;
    struct tm* info;
    time(&t);
    info = localtime(&t);

    // Print Header
    md_file_printf(file, "# This file was created %s\n", asctime(info));
    md_file_printf(file, "# Created by:\n");
    md_file_printf(file, "# VIAMD \n");

    // Print Legend Meta
    md_file_printf(file, "@    title \"VIAMD Properties\"\n");
    md_file_printf(file, "@    xaxis  label \"Time\"\n");
    md_file_printf(file, "@TYPE xy\n");
    md_file_printf(file, "@ view 0.15, 0.15, 0.75, 0.85\n");
    md_file_printf(file, "@ legend on\n");
    md_file_printf(file, "@ legend box on\n");
    md_file_printf(file, "@ legend loctype view\n");
    md_file_printf(file, "@ legend 0.78, 0.8\n");
    md_file_printf(file, "@ legend length 2\n");

    for (int j = 0; j < num_columns; ++j) {
        md_file_printf(file, "@ s%i legend \"%s\"\n", j, column_labels[j]);
    }

    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_columns; ++j) {
            md_file_printf(file, "%12.6f ", column_data[j][i]);
        }
        md_file_printf(file, "\n");
    }

    md_file_close(file);
    return true;
}

static bool export_csv(const float* column_data[], const char* column_labels[], int num_columns, int num_rows, str_t filename) {
    ASSERT(column_data);
    ASSERT(column_labels);
    ASSERT(num_columns >= 0);
    ASSERT(num_rows >= 0);

    md_file_o* file = md_file_open(filename, MD_FILE_WRITE);
    if (!file) {
        md_printf(MD_LOG_TYPE_ERROR, "Failed to open file '%.*s' to write data.", (int)filename.len, filename.ptr);
        return false;
    }
    
    for (int i = 0; i < num_columns; ++i) {
        md_file_printf(file, "%s,", column_labels[i]);
    }
    md_file_printf(file, "\n");

    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_columns; ++j) {
            md_file_printf(file, "%.6g,", column_data[j][i]);
        }
        md_file_printf(file, "\n");
    }

    md_file_close(file);
    return true;
}

static bool export_cube(const ApplicationData& data, const md_script_property_t* prop, str_t filename) {
    // @NOTE: First we need to extract some meta data for the cube format, we need the atom indices/bits for any SDF
    // And the origin + extent of the volume in spatial coordinates (Ångström)

    // Copy mol and replace with initial coords
    md_molecule_t mol = data.mold.mol;

    int64_t stride = ROUND_UP(data.mold.mol.atom.count, md_simd_widthf);
    float* coords = (float*)md_alloc(frame_allocator, stride * sizeof(float) * 3);
    mol.atom.x = coords + stride * 0;
    mol.atom.y = coords + stride * 1;
    mol.atom.z = coords + stride * 2;

    md_trajectory_load_frame(&data.mold.traj, 0, NULL, mol.atom.x, mol.atom.y, mol.atom.z);

    md_script_visualization_args_t args = {
        .token = prop->vis_token,
        .ir = &data.mold.script.ir,
        .mol = &data.mold.mol,
        .traj = &data.mold.traj,
        .alloc = frame_allocator,
        .flags = MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_SDF,
    };
    md_script_visualization_t vis = {0};
    if (md_script_visualization_init(&vis, args)) {

        md_file_o* file = md_file_open(filename, MD_FILE_WRITE);
        if (!file) {
            md_printf(MD_LOG_TYPE_ERROR, "Failed to open file '%.*s' in order to write to it.", (int)filename.len, filename.ptr);
            return false;
        }

        // Two comment lines
        md_file_printf(file, "EXPORTED DENSITY VOLUME FROM VIAMD, UNITS IN BOHR\n");
        md_file_printf(file, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");

        if (vis.sdf.count > 0) {
            const float angstrom_to_bohr = (float)(1.0 / 0.529177210903);

            // transformation matrix from world to volume
            const mat4_t M = vis.sdf.matrices[0];
            const md_exp_bitfield_t* bf = &vis.sdf.structures[0];
            const int num_atoms = (int)md_bitfield_popcount(bf);
            const int vol_dim[3] = {prop->data.dim[0], prop->data.dim[1], prop->data.dim[2]};
            const double extent = vis.sdf.extent * angstrom_to_bohr;
            const double voxel_ext[3] = {
                (double)extent / (double)prop->data.dim[0],
                (double)extent / (double)prop->data.dim[1],
                (double)extent / (double)prop->data.dim[2],
            };

            const double half_ext = extent * 0.5;

            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", num_atoms, -half_ext, -half_ext, -half_ext);
            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol_dim[0], voxel_ext[0], 0, 0);
            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol_dim[1], 0, voxel_ext[1], 0);
            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol_dim[2], 0, 0, voxel_ext[2]);

            int64_t beg_bit = bf->beg_bit;
            int64_t end_bit = bf->end_bit;
            while ((beg_bit = md_bitfield_scan(bf, beg_bit, end_bit)) != 0) {
                int64_t i = beg_bit - 1;
                vec3_t coord = {mol.atom.x[i], mol.atom.y[i], mol.atom.z[i]};
                coord = mat4_mul_vec3(M, coord, 1.0f);
                md_element_t elem = mol.atom.element[i];
                md_file_printf(file, "%5i %12.6f %12.6f %12.6f %12.6f\n", elem, (float)elem, coord.x, coord.y, coord.z);
            }

            //md_file_printf(file, "%5i %5i\n", 1, 1);

            // Write density data
            int count = 0;
            for (int x = 0; x < vol_dim[0]; ++x) {
                for (int y = 0; y < vol_dim[1]; ++y) {
                    for (int z = 0; z < vol_dim[2]; ++z) {
                        int idx = z * vol_dim[0] * vol_dim[1] + y * vol_dim[0] + x;
                        float val = prop->data.values[idx];
                        md_file_printf(file, " %12.6E", val);
                        if (++count % 6 == 0) md_file_printf(file, "\n");
                    }
                }
            }
        }

        md_file_close(file);
    } else {
        md_print(MD_LOG_TYPE_ERROR, "Failed to visualize volume for export.");
        return false;
    }

    return true;
}

#define APPEND_BUF(buf, len, fmt, ...) (len += snprintf(buf + len, MAX(0, (int)sizeof(buf) - len), fmt, ##__VA_ARGS__) + 1)

static void draw_property_export_window(ApplicationData* data) {
    ASSERT(data);

    enum ExportType {
        Temporal = 0,
        Distribution,
        DensityVolume
    };

    struct ExportFormat {
        str_t label;
        str_t extension;
    };

    ExportFormat table_formats[] {
        {make_cstr("XVG"), make_cstr("xvg")},
        {make_cstr("CSV"), make_cstr("csv")}
    };

    ExportFormat volume_formats[] {
        {make_cstr("Gaussian Cube"), make_cstr("cube")},
        //{make_cstr("DAT + RAW"), make_cstr("dat")},
    };

    if (ImGui::Begin("Property Export", &data->show_property_export_window)) {
        static ExportType type = Temporal;
        ImGui::PushItemWidth(200);
        ImGui::Combo("Data Type", (int*)(&type), "Temporal\0Distribution\0Density Volume\0");
        ImGui::Separator();
        if (type == Temporal || type == Distribution) {
            static int format = 0;
            static int col_options[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
            static int num_columns = 4;
            static int values_option = 0;
            static int num_bins = 64;
            static bool temporal_filter = false;

            const int MIN_NUM_BINS = 8;
            const int MAX_NUM_BINS = 1024;

            struct ColData {
                const char* label;
                const float* values;
                int num_values;
                int dim;
            };

            ColData* col_data = 0;
            if (type == Temporal) {
                ColData time_data = {"Time", data->timeline.x_values, (int)md_array_size(data->timeline.x_values), 1};
                md_array_push(col_data, time_data, frame_allocator);
                for (int i = 0; i < (int)md_array_size(data->display_properties); ++i) {
                    const DisplayProperty& dp = data->display_properties[i];
                    if (dp.full_prop->type == MD_SCRIPT_PROPERTY_TYPE_TEMPORAL) {
                        ColData prop_data = {dp.lbl.cstr(), dp.full_prop->data.values, (int)dp.full_prop->data.num_values, dp.full_prop->data.dim[0]};
                        md_array_push(col_data, prop_data, frame_allocator);
                    }
                }
            } else if (type == Distribution) {
                for (int i = 0; i < (int)md_array_size(data->display_properties); ++i) {
                    const DisplayProperty& dp = data->display_properties[i];
                    if (dp.full_prop->type == MD_SCRIPT_PROPERTY_TYPE_TEMPORAL) {
                        ColData prop_data = {dp.lbl.cstr(), dp.full_hist.bin, ARRAY_SIZE(dp.full_hist.bin), 1};
                        md_array_push(col_data, prop_data, frame_allocator);
                        if (data->timeline.filter.enabled) {
                            str_t lbl = alloc_printf(frame_allocator, "%.*s(filt)", (int)dp.lbl.length(), dp.lbl.cstr());
                            ColData filt_data = {lbl.ptr, dp.filt_hist.bin, ARRAY_SIZE(dp.filt_hist.bin), 1};
                            md_array_push(col_data, filt_data, frame_allocator);
                        }
                    } else if (dp.full_prop->type == MD_SCRIPT_PROPERTY_TYPE_DISTRIBUTION) {
                        ColData prop_data = {dp.lbl.cstr(), dp.full_prop->data.values, (int)dp.full_prop->data.num_values, dp.full_prop->data.dim[0]};
                        md_array_push(col_data, prop_data, frame_allocator);
                        if (data->timeline.filter.enabled) {
                            str_t lbl = alloc_printf(frame_allocator, "%.*s(filt)", (int)dp.lbl.length(), dp.lbl.cstr());
                            ColData filt_data = {lbl.ptr, dp.filt_prop->data.values, (int)dp.filt_prop->data.num_values, dp.filt_prop->data.dim[0]};
                            md_array_push(col_data, filt_data, frame_allocator);
                        }
                    }
                }
            }

            ImGui::Text("Column Layout");
            ImGui::SliderInt("Num Columns", &num_columns, 1, ARRAY_SIZE(col_options));

            int num_col_data = (int)md_array_size(col_data);
            if (ImGui::BeginTable("Columns", num_columns, ImGuiTableFlags_Borders)) {
                for (int i = 0; i < num_columns; ++i) {
                    ImGui::TableNextColumn();
                    ImGui::PushItemWidth(-1);
                    ImGui::PushID(i);
                    const char* preview = (0 <= col_options[i] && col_options[i] < num_col_data) ? col_data[col_options[i]].label : "";
                    if (ImGui::BeginCombo("##col", preview)) {
                        if (ImGui::Selectable("", col_options[i] == -1)) {
                            col_options[i] = -1;
                        }
                        for (int j = 0; j < num_col_data; ++j) {
                            ImGui::PushID(j);
                            if (ImGui::Selectable(col_data[j].label, col_options[i] == j)) {
                                col_options[i] = j;
                            }
                            ImGui::PopID();
                        }
                        ImGui::EndCombo();
                    }
                    ImGui::PopID();
                }
                ImGui::EndTable();
            }

            if (type == Distribution) {
                if (ImGui::SliderInt("Num Bins", &num_bins, MIN_NUM_BINS, MAX_NUM_BINS, "%d", ImGuiSliderFlags_Logarithmic)) {
                    const int up   = next_power_of_two32(num_bins);
                    const int down = up / 2;
                    num_bins = abs(num_bins - down) < abs(num_bins - up) ? down : up;
                }
            }

            ImGui::Separator();
            ImGui::Text("File format");
            ImGui::Combo("##Export Format", (int*)(&format), "XVG\0CSV\0");
            if (ImGui::Button("Export")) {
                application::FileDialogResult file_dialog = application::file_dialog(application::FileDialogFlags_Save, {}, table_formats[format].extension);
                if (file_dialog.result == application::FileDialogResult::Ok) {
                    StrBuf<1024> path = str_t{file_dialog.path, file_dialog.path_len};
                    if (str_empty(extract_ext(path))) {
                        path += ".";
                        path += table_formats[format].extension;
                    }

                    const float** column_data = 0;
                    const char** column_labels = 0;
                    for (int i = 0; i < num_columns; ++i) {
                        int idx = col_options[i];
                        if (0 <= idx && idx < num_col_data) {
                            ASSERT(col_data[idx].dim >= 1);
                            if (type == Temporal) {
                                for (int j = 0; j < col_data[idx].dim; ++j) {
                                    int stride = col_data[idx].num_values / col_data[idx].dim;
                                    md_array_push(column_data, col_data[idx].values + stride * j, frame_allocator);
                                    if (col_data[idx].dim == 1) {
                                        md_array_push(column_labels, col_data[idx].label, frame_allocator);
                                    } else {
                                        str_t lbl = alloc_printf(frame_allocator, "%s[%i]", col_data[idx].label, j);
                                        md_array_push(column_labels, lbl.ptr, frame_allocator);
                                    }
                                }
                            } else if (type == Distribution) {
                                float* bins = (float*)md_alloc(frame_allocator, num_bins * sizeof(float));
                                downsample_histogram(bins, num_bins, col_data[idx].values, col_data[idx].num_values);
                                md_array_push(column_data, bins, frame_allocator);
                                md_array_push(column_labels, col_data[idx].label, frame_allocator);
                            }
                        }
                    }

                    int num_rows = (type == Distribution) ? num_bins : md_trajectory_num_frames(&data->mold.traj);

                    switch (format) {
                    case 0:
                        export_xvg(column_data, column_labels, md_array_size(column_data), num_rows, path);
                        break;
                    case 1:
                        export_csv(column_data, column_labels, md_array_size(column_data), num_rows, path);
                        break;
                    default:
                        ASSERT(false);
                    }
                }
            }
        }
        else if (type == DensityVolume) {
            static int format = 0;
            static int prop_idx = 0;
            static bool temporal_filter = false;
            const DisplayProperty** props = 0;

            for (int i = 0; i < (int)md_array_size(data->display_properties); ++i) {
                if (data->display_properties[i].full_prop->type == MD_SCRIPT_PROPERTY_TYPE_VOLUME) {
                    md_array_push(props, &data->display_properties[i], frame_allocator);
                }
            }

            const int num_props = (int)md_array_size(props);
            if (num_props > 0) {
                prop_idx = CLAMP(prop_idx, 0, num_props);
                if (ImGui::BeginCombo("Source", props[prop_idx]->lbl.cstr())) {
                    for (int i = 0; i < num_props; ++i) {
                        if (ImGui::Selectable(props[i]->lbl.cstr(), prop_idx == i)) {
                            prop_idx = i;
                        }
                    }
                    ImGui::EndCombo();
                }

                if (data->timeline.filter.enabled) {
                    ImGui::SameLine();
                    ImGui::Checkbox("Use Temporal Filter", &temporal_filter);
                } else {
                    temporal_filter = false;
                }

                format = CLAMP(format, 0, ARRAY_SIZE(volume_formats));
                if (ImGui::BeginCombo("Format", volume_formats[format].label.ptr)) {
                    for (int i = 0; i < (int)ARRAY_SIZE(volume_formats); ++i) {
                        if (ImGui::Selectable(volume_formats[i].label.ptr, format == i)) format = i;
                    }
                    ImGui::EndCombo();
                }

                if (ImGui::Button("Export")) {
                    application::FileDialogResult file_dialog = application::file_dialog(application::FileDialogFlags_Save, {}, volume_formats[format].extension);
                    if (file_dialog.result == application::FileDialogResult::Ok) {
                        StrBuf<1024> path = str_t{file_dialog.path, file_dialog.path_len};
                        if (str_empty(extract_ext(path))) {
                            path += ".";
                            path += volume_formats[format].extension;
                        }
                        switch (format) {
                        case 0:
                            export_cube(*data, temporal_filter ? props[prop_idx]->filt_prop : props[prop_idx]->filt_prop, path);
                            break;
                        case 1:
                            // @TODO: Export dat + raw
                            break;
                        default:
                            ASSERT(false);
                        }
                    }
                }
            }
        }
        else {
            ASSERT(false);
        }
        ImGui::PopItemWidth();
    }
    ImGui::End();
}

// static void draw_density_volume_clip_plane_widgets(ApplicationData* data) { ASSERT(data); }

// #gbuffer
static void init_gbuffer(GBuffer* gbuf, int width, int height) {
    ASSERT(gbuf);

    bool attach_textures_deferred = false;
    if (!gbuf->deferred.fbo) {
        glGenFramebuffers(1, &gbuf->deferred.fbo);
        attach_textures_deferred = true;
    }

    if (!gbuf->deferred.depth) glGenTextures(1, &gbuf->deferred.depth);
    if (!gbuf->deferred.color) glGenTextures(1, &gbuf->deferred.color);
    if (!gbuf->deferred.normal) glGenTextures(1, &gbuf->deferred.normal);
    if (!gbuf->deferred.velocity) glGenTextures(1, &gbuf->deferred.velocity);
    if (!gbuf->deferred.post_tonemap) glGenTextures(1, &gbuf->deferred.post_tonemap);
    if (!gbuf->deferred.picking) glGenTextures(1, &gbuf->deferred.picking);
    if (!gbuf->pbo_picking.color[0]) glGenBuffers(2, gbuf->pbo_picking.color);
    if (!gbuf->pbo_picking.depth[0]) glGenBuffers(2, gbuf->pbo_picking.depth);

    glBindTexture(GL_TEXTURE_2D, gbuf->deferred.depth);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH24_STENCIL8, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->deferred.color);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->deferred.normal);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16, width, height, 0, GL_RG, GL_UNSIGNED_SHORT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->deferred.velocity);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, width, height, 0, GL_RG, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->deferred.post_tonemap);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->deferred.picking);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.color[0]);
    glBufferData(GL_PIXEL_PACK_BUFFER, 4, NULL, GL_DYNAMIC_READ);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.color[1]);
    glBufferData(GL_PIXEL_PACK_BUFFER, 4, NULL, GL_DYNAMIC_READ);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.depth[0]);
    glBufferData(GL_PIXEL_PACK_BUFFER, 4, NULL, GL_DYNAMIC_READ);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.depth[1]);
    glBufferData(GL_PIXEL_PACK_BUFFER, 4, NULL, GL_DYNAMIC_READ);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    glBindTexture(GL_TEXTURE_2D, 0);

    gbuf->width = width;
    gbuf->height = height;

    const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
                                   GL_COLOR_ATTACHMENT_POST_TONEMAP, GL_COLOR_ATTACHMENT_PICKING};

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf->deferred.fbo);
    if (attach_textures_deferred) {
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, gbuf->deferred.depth, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_STENCIL_ATTACHMENT, GL_TEXTURE_2D, gbuf->deferred.depth, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_COLOR, GL_TEXTURE_2D, gbuf->deferred.color, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_NORMAL, GL_TEXTURE_2D, gbuf->deferred.normal, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_VELOCITY, GL_TEXTURE_2D, gbuf->deferred.velocity, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_POST_TONEMAP, GL_TEXTURE_2D, gbuf->deferred.post_tonemap, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_PICKING, GL_TEXTURE_2D, gbuf->deferred.picking, 0);
    }
    ASSERT(glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
    glDrawBuffers(ARRAY_SIZE(draw_buffers), draw_buffers);
    glClearColor(0, 0, 0, 0);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
}

static void destroy_gbuffer(GBuffer* gbuf) {
    ASSERT(gbuf);
    if (gbuf->deferred.fbo) glDeleteFramebuffers(1, &gbuf->deferred.fbo);
    if (gbuf->deferred.depth) glDeleteTextures(1, &gbuf->deferred.depth);
    if (gbuf->deferred.color) glDeleteTextures(1, &gbuf->deferred.color);
    if (gbuf->deferred.normal) glDeleteTextures(1, &gbuf->deferred.normal);
    if (gbuf->deferred.post_tonemap) glDeleteTextures(1, &gbuf->deferred.post_tonemap);
    if (gbuf->deferred.picking) glDeleteTextures(1, &gbuf->deferred.picking);

    if (gbuf->pbo_picking.color[0]) glDeleteBuffers(2, gbuf->pbo_picking.color);
    if (gbuf->pbo_picking.depth[0]) glDeleteBuffers(2, gbuf->pbo_picking.depth);
}

static void update_md_buffers(ApplicationData* data) {
    ASSERT(data);
    const auto& mol = data->mold.mol;

    if (data->mold.dirty_buffers & MolBit_DirtyPosition) {
        md_gl_molecule_update_atom_previous_position(&data->mold.gl_mol);
        md_gl_molecule_set_atom_position(&data->mold.gl_mol, 0, (uint32_t)mol.atom.count, mol.atom.x, mol.atom.y, mol.atom.z, 0);
    }

    if (data->mold.dirty_buffers & MolBit_DirtyRadius) {
        md_gl_molecule_set_atom_radius(&data->mold.gl_mol, 0, (uint32_t)mol.atom.count, mol.atom.radius, 0);
    }

    if (data->mold.dirty_buffers & MolBit_DirtyFlags) {
        if (data->mold.mol.atom.flags) {
            for (int64_t i = 0; i < mol.atom.count; i++) {
                uint8_t flags = 0;
                flags |= md_bitfield_test_bit(&data->selection.current_highlight_mask, i)     ? AtomBit_Highlighted : 0;
                flags |= md_bitfield_test_bit(&data->selection.current_selection_mask, i)     ? AtomBit_Selected : 0;
                flags |= md_bitfield_test_bit(&data->representations.atom_visibility_mask, i) ? AtomBit_Visible : 0;
                data->mold.mol.atom.flags[i] = flags;
            }
            md_gl_molecule_set_atom_flags(&data->mold.gl_mol, 0, (uint32_t)mol.atom.count, mol.atom.flags, 0);
        }
    }

    if (data->mold.dirty_buffers & MolBit_DirtyBonds) {
        md_gl_molecule_set_covalent_bonds(&data->mold.gl_mol, 0, (uint32_t)mol.covalent_bond.count, mol.covalent_bond.bond, 0);
    }

    if (data->mold.dirty_buffers & MolBit_DirtySecondaryStructure) {
        if (mol.backbone.secondary_structure) {
            md_gl_molecule_set_backbone_secondary_structure(&data->mold.gl_mol, 0, (uint32_t)mol.backbone.count, mol.backbone.secondary_structure, 0);
        }
    }

    data->mold.dirty_buffers = 0;
}

static void interrupt_async_tasks(ApplicationData* data) {
    if (data->tasks.backbone_computations.id != 0) {
        task_system::interrupt_task(data->tasks.backbone_computations);
    }

    if (data->tasks.evaluate_full.id != 0) {
        md_script_eval_interrupt(&data->mold.script.full_eval);
        task_system::interrupt_task(data->tasks.evaluate_full);
    }

    if (data->tasks.evaluate_filt.id != 0) {
        md_script_eval_interrupt(&data->mold.script.filt_eval);
        task_system::interrupt_task(data->tasks.evaluate_filt);
    }

    task_system::wait_for_task(data->tasks.backbone_computations);
    task_system::wait_for_task(data->tasks.evaluate_full);
    task_system::wait_for_task(data->tasks.evaluate_filt);

    data->tasks.backbone_computations.id = 0;
    data->tasks.evaluate_full.id = 0;
    data->tasks.evaluate_filt.id = 0;
}

// #trajectorydata
static void free_trajectory_data(ApplicationData* data) {
    ASSERT(data);
    interrupt_async_tasks(data);

    if (md_trajectory_num_frames(&data->mold.traj)) {
        load::traj::close(&data->mold.traj);
    }
    memset(data->simulation_box.box.elem, 0, sizeof(mat3_t));
    md_array_shrink(data->timeline.x_values,  0);
    md_array_shrink(data->display_properties, 0);

    data->shape_space.input_valid = false;
    data->shape_space.num_frames = 0;
    data->shape_space.num_structures = 0;
    md_array_shrink(data->shape_space.weights, 0);
    md_array_shrink(data->shape_space.coords, 0);
}


static void init_trajectory_data(ApplicationData* data) {
    int64_t num_frames = md_trajectory_num_frames(&data->mold.traj);
    if (num_frames > 0) {
        int64_t min_frame = 0;
        int64_t max_frame = num_frames - 1;
        double min_time = 0;
        double max_time = 0;
        {
            md_trajectory_frame_header_t header;
            md_trajectory_load_frame(&data->mold.traj, min_frame, &header, 0, 0, 0);
            min_time = header.timestamp;
        }

        {
            md_trajectory_frame_header_t header;
            md_trajectory_load_frame(&data->mold.traj, max_frame, &header, 0, 0, 0);
            max_time = header.timestamp;
        }

        data->timeline.view_range = {min_time, max_time};
        data->timeline.filter.beg_frame = min_frame;
        data->timeline.filter.end_frame = max_frame;

        data->shape_space.evaluate = true;

        md_array_resize(data->timeline.x_values, num_frames, persistent_allocator);
        for (int64_t i = 0; i < num_frames; ++i) {
            // We estimate the time values to grow linearly between min_time and max_time
            // This will be updated to an exact value later when prefetching the frames
            double t = i / (double)max_frame;
            data->timeline.x_values[i] = lerp(min_time, max_time, t);
        }

        data->animation.frame = CLAMP(data->animation.frame, (double)min_frame, (double)max_frame);
        int64_t frame_idx = CLAMP((int64_t)(data->animation.frame + 0.5), 0, max_frame);

        md_trajectory_frame_header_t header;
        md_trajectory_load_frame(&data->mold.traj, frame_idx, &header, data->mold.mol.atom.x, data->mold.mol.atom.y, data->mold.mol.atom.z);
        memcpy(&data->simulation_box.box, header.box, sizeof(header.box));
        data->mold.dirty_buffers |= MolBit_DirtyPosition;

        update_md_buffers(data);
        md_gl_molecule_update_atom_previous_position(&data->mold.gl_mol); // Do this explicitly to update the previous position to avoid motion blur trails

        if (data->mold.mol.backbone.count > 0) {
            data->trajectory_data.secondary_structure.stride = data->mold.mol.backbone.count;
            data->trajectory_data.secondary_structure.count = data->mold.mol.backbone.count * num_frames;
            md_array_resize(data->trajectory_data.secondary_structure.data, data->mold.mol.backbone.count * num_frames, persistent_allocator);

            data->trajectory_data.backbone_angles.stride = data->mold.mol.backbone.count,
                data->trajectory_data.backbone_angles.count = data->mold.mol.backbone.count * num_frames,
                md_array_resize(data->trajectory_data.backbone_angles.data, data->mold.mol.backbone.count * num_frames, persistent_allocator);

            // Launch work to prefetch frames
            launch_prefetch_job(data);

            // Launch work to compute the values
            if (data->tasks.backbone_computations.id != 0) {
                task_system::interrupt_and_wait(data->tasks.backbone_computations); // This should never happen
            }
            data->tasks.backbone_computations = task_system::enqueue_pool("Backbone Operations", (uint32_t)num_frames, [data](task_system::TaskSetRange range) {
                const auto& mol = data->mold.mol;
                const int64_t stride = ROUND_UP(mol.atom.count, md_simd_widthf);
                const int64_t bytes = stride * sizeof(float) * 3;
                float* coords = (float*)md_alloc(persistent_allocator, bytes);
                defer { md_free(persistent_allocator, coords, bytes); };
                float* x = coords + stride * 0;
                float* y = coords + stride * 1;
                float* z = coords + stride * 2;

                for (uint32_t frame_idx = range.beg; frame_idx < range.end; ++frame_idx) {
                    md_trajectory_load_frame(&data->mold.traj, frame_idx, NULL, x, y, z);

                    md_util_backbone_angle_args_t bb_args = {
                        .atom = {
                            .count = mol.atom.count,
                            .x = x,
                            .y = y,
                            .z = z,
                        },
                        .backbone = {
                            .count = mol.backbone.count,
                            .atoms = mol.backbone.atoms,
                        },
                        .chain = {
                            .count = mol.chain.count,
                            .backbone_range = mol.chain.backbone_range,
                        }
                    };
                    md_util_compute_backbone_angles(data->trajectory_data.backbone_angles.data + data->trajectory_data.backbone_angles.stride * frame_idx, data->trajectory_data.backbone_angles.stride, &bb_args);

                    md_util_secondary_structure_args_t ss_args = {
                        .atom = {
                            .count = data->mold.mol.atom.count,
                            .x = x,
                            .y = y,
                            .z = z,
                        },
                        .backbone = {
                            .count = mol.backbone.count,
                            .atoms = mol.backbone.atoms,
                        },
                        .chain = {
                            .count = data->mold.mol.chain.count,
                            .backbone_range = data->mold.mol.chain.backbone_range,
                        }
                    };
                    md_util_compute_secondary_structure(data->trajectory_data.secondary_structure.data + data->trajectory_data.secondary_structure.stride * frame_idx, data->trajectory_data.secondary_structure.stride, &ss_args);
                }
            });
        }
    }
}

static bool load_trajectory_data(ApplicationData* data, str_t filename) {
    interrupt_async_tasks(data);
    free_trajectory_data(data);
    data->files.trajectory = "";
    data->animation.frame = 0;

    if (load::traj::open_file(&data->mold.traj, filename, &data->mold.mol, persistent_allocator)) {
        data->files.trajectory = filename;
        init_trajectory_data(data);
        return true;
    }

    return false;
}

// #moleculedata
static void free_molecule_data(ApplicationData* data) {
    ASSERT(data);
    interrupt_async_tasks(data);

    if (data->mold.mol.atom.count) {
        data->files.molecule = "";
        load::mol::free(&data->mold.mol);
        data->mold.mol_alloc = NULL;
    }
    free_trajectory_data(data);

    md_bitfield_clear(&data->selection.current_selection_mask);
    md_bitfield_clear(&data->selection.current_highlight_mask);
    md_script_ir_free(&data->mold.script.ir);
    md_script_eval_free(&data->mold.script.full_eval);
}

static void init_molecule_data(ApplicationData* data) {
    if (data->mold.mol.atom.count) {
        data->mold.mol_alloc = load::mol::get_internal_allocator(&data->mold.mol);
        const auto& mol = data->mold.mol;
        data->picking.idx = INVALID_PICKING_IDX;
        data->selection.hovered = -1;
        data->selection.right_clicked = -1;

        md_gl_molecule_desc_t desc = {
            .atom = {
                .count = (uint32_t)mol.atom.count,
                .x = mol.atom.x,
                .y = mol.atom.y,
                .z = mol.atom.z,
                .radius = mol.atom.radius
            },
            .covalent_bond = {
                .count = (uint32_t)mol.covalent_bond.count,
                .atom_bond = mol.covalent_bond.bond
            },
            .residue = {
                .count = (uint32_t)mol.residue.count,
                .atom_range = mol.residue.atom_range,
            },
            .backbone = {
                .count = (uint32_t)mol.backbone.count,
                .atoms = mol.backbone.atoms,
                .secondary_structure = mol.backbone.secondary_structure
            },
            .chain = {
                .count = (uint32_t)mol.chain.count,
                .backbone_range = mol.chain.backbone_range
            }
        };
        md_gl_molecule_init(&data->mold.gl_mol, &desc);

        init_all_representations(data);
        update_all_representations(data);
        data->mold.script.compile_ir = true;
    }
}

static void launch_prefetch_job(ApplicationData* data) {
    uint32_t num_frames = (uint32_t)md_trajectory_num_frames(&data->mold.traj);
    if (!num_frames) return;

    if (data->tasks.prefetch_frames.id != 0) {
        task_system::interrupt_and_wait(data->tasks.prefetch_frames); // This should never happen
    }
    data->tasks.prefetch_frames = task_system::enqueue_pool("Prefetch Frames", num_frames, [data](task_system::TaskSetRange range) {
        for (uint32_t i = range.beg; i < range.end; ++i) {
            md_trajectory_frame_header_t header;
            md_trajectory_load_frame(&data->mold.traj, i, &header, 0, 0, 0);
            data->timeline.x_values[i] = header.timestamp;
        }
    });


#if 0
#define NUM_SLOTS 64

    // This is the number of slots we have to work with in parallel.
    // This number should ideally be more than the number of cores available.
    // We pre-allocate the number of slots * max frame data size, so don't go bananas here if you want to save some on memory.
    task_system::enqueue_pool("Preloading frames", 1, [data, num_frames](task_system::TaskSetRange)
        {
            timestamp_t t0 = md_os_time_current();
            const int64_t slot_size = md_trajectory_max_frame_data_size(&data->mold.traj);
            void* slot_mem = md_alloc(default_allocator, slot_size * NUM_SLOTS);

            void* slots[NUM_SLOTS];
            atomic_queue::AtomicQueue<uint32_t, NUM_SLOTS, 0xFFFFFFFF> slot_queue;

            for (uint32_t i = 0; i < NUM_SLOTS; ++i) {
                slots[i] = (char*)slot_mem + i * slot_size;
                slot_queue.push(i);
            }

            // Iterate over all frames and load the raw data, then spawn a task for each frame
            for (uint32_t i = 0; i < num_frames; ++i) {
                uint32_t slot_idx = slot_queue.pop();
                ASSERT(slot_idx < NUM_SLOTS);

                void* frame_mem = slots[slot_idx];
                const int64_t frame_size = data->mold.traj.fetch_frame_data(data->mold.traj.inst, i, NULL);
                data->mold.traj.fetch_frame_data(data->mold.traj.inst, i, frame_mem);

                // Spawn task: Load, Decode and postprocess
                //printf("Spawning task to decode frame %i\n", i);
                auto id = task_system::enqueue_pool("##Decode frame", 1, [slot_idx, frame_mem, frame_size, frame_idx = i, &slot_queue, data](task_system::TaskSetRange)
                    {
                        md_frame_data_t* frame_data;
                        md_frame_cache_lock_t* lock;
                        if (md_frame_cache_reserve_frame(&data->mold.frame_cache, frame_idx, &frame_data, &lock)) {
                            data->mold.traj.decode_frame_data(data->mold.traj.inst, frame_mem, frame_size, &frame_data->header, frame_data->x, frame_data->y, frame_data->z);

                            // Free the data slot directly here
                            slot_queue.push(slot_idx);

                            // deperiodize
                            const md_molecule_t& mol = data->mold.mol;
                            md_util_apply_pbc_args_t args = {
                                .atom = {
                                    .count = mol.atom.count,
                                    .x = frame_data->x,
                                    .y = frame_data->y,
                                    .z = frame_data->z,
                            },
                            .residue = {
                                    .count = mol.residue.count,
                                    .atom_range = mol.residue.atom_range,
                            },
                            .chain = {
                                    .count = mol.chain.count,
                                    .residue_range = mol.chain.residue_range,
                            }
                            };
                            memcpy(args.pbc.box, frame_data->header.box, sizeof(args.pbc.box));
                            md_util_apply_pbc(frame_data->x, frame_data->y, frame_data->z, mol.atom.count, args);

                            if (mol.backbone.count > 0) {
                                md_util_backbone_angle_args_t bb_args = {
                                    .atom = {
                                        .count = mol.atom.count,
                                        .x = frame_data->x,
                                        .y = frame_data->y,
                                        .z = frame_data->z,
                                },
                                .backbone = {
                                        .count = mol.backbone.count,
                                        .atoms = mol.backbone.atoms,
                                },
                                .chain = {
                                        .count = mol.chain.count,
                                        .backbone_range = mol.chain.backbone_range,
                                }
                                };
                                md_util_compute_backbone_angles(data->trajectory_data.backbone_angles.data + data->trajectory_data.backbone_angles.stride * frame_idx, data->trajectory_data.backbone_angles.stride, &bb_args);

                                md_util_secondary_structure_args_t ss_args = {
                                    .atom = {
                                        .count = data->mold.mol.atom.count,
                                        .x = frame_data->x,
                                        .y = frame_data->y,
                                        .z = frame_data->z,
                                },
                                .backbone = {
                                        .count = mol.backbone.count,
                                        .atoms = mol.backbone.atoms,
                                },
                                .chain = {
                                        .count = data->mold.mol.chain.count,
                                        .backbone_range = data->mold.mol.chain.backbone_range,
                                }
                                };
                                md_util_compute_secondary_structure(data->trajectory_data.secondary_structure.data + data->trajectory_data.secondary_structure.stride * frame_idx, data->trajectory_data.secondary_structure.stride, &ss_args);
                            }

                            md_frame_cache_release_frame_lock(lock);
                        } else {
                            slot_queue.push(slot_idx);
                        }
                        //printf("Finished decoding frame %i\n", frame_idx);
                    }
                );
            }

            //printf("Sitting back waiting for tasks to complete...\n");
            while (!slot_queue.was_full()) {
                _mm_pause(); // Back off for a bit
            }

            md_free(default_allocator, slot_mem, slot_size * NUM_SLOTS);

            timestamp_t t1 = md_os_time_current();
            md_printf(MD_LOG_TYPE_INFO, "Frame preload took %.2f seconds", md_os_time_delta_in_s(t0, t1));
        });
#undef NUM_SLOTS
#endif
}

static bool load_dataset_from_file(ApplicationData* data, str_t filename) {
    ASSERT(data);

    filename = md_os_path_make_canonical(filename, frame_allocator);

    if (filename.len) {
        if (load::mol::is_extension_supported(filename)) {
            if (compare_str(data->files.molecule, filename)) {
                // Same as loaded file
                return true;
            }

            interrupt_async_tasks(data);
            free_molecule_data(data);
            free_trajectory_data(data);
            data->files.molecule = filename;
            data->files.trajectory = "";

            md_printf(MD_LOG_TYPE_INFO, "Attempting to load molecular data from file '%.*s'", filename.len, filename.ptr);
            if (!load::mol::load_file(&data->mold.mol, filename, persistent_allocator)) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to load molecular data");
                data->files.molecule = "";
                return false;
            }
            init_molecule_data(data);

            // @NOTE: Some files contain both atomic coordinates and trajectory
            if (load::traj::is_extension_supported(filename)) {
                md_print(MD_LOG_TYPE_INFO, "File also contains trajectory, attempting to load trajectory");
                load_trajectory_data(data, filename);
            }
            return true;
        } else if (load::traj::is_extension_supported(filename)) {
            if (!data->mold.mol.atom.count) {
                md_print(MD_LOG_TYPE_ERROR, "Before loading a trajectory, molecular data needs to be present");
                return false;
            }

            if (compare_str(data->files.trajectory, filename)) {
                // Same as loaded file
                return true;
            }

            return load_trajectory_data(data, filename);
        } else {
            md_print(MD_LOG_TYPE_ERROR, "File extension not supported");
        }
    }
    return false;
}

// ### WORKSPACE ###
static RepresentationType get_rep_type(str_t str) {
    if (compare_str_cstr(str, "SPACE_FILL"))
        return RepresentationType::SpaceFill;
    else if (compare_str_cstr(str, "LICORICE"))
        return RepresentationType::Licorice;
    else if (compare_str_cstr(str, "BALL_AND_STICK"))    // Ball and stick is removed for now
        return RepresentationType::SpaceFill;
    else if (compare_str_cstr(str, "RIBBONS"))
        return RepresentationType::Ribbons;
    else if (compare_str_cstr(str, "CARTOON"))
        return RepresentationType::Cartoon;
    else
        return RepresentationType::SpaceFill;
}

static str_t get_rep_type_name(RepresentationType type) {
    switch (type) {
        case RepresentationType::SpaceFill:
            return make_cstr("SPACE_FILL");
        case RepresentationType::Licorice:
            return make_cstr("LICORICE");
        /*case RepresentationType::BallAndStick:
            return make_cstr("BALL_AND_STICK");*/
        case RepresentationType::Ribbons:
            return make_cstr("RIBBONS");
        case RepresentationType::Cartoon:
            return make_cstr("CARTOON");
        default:
            return make_cstr("UNKNOWN");
    }
}

static ColorMapping get_color_mapping(str_t str) {
    if (compare_str_cstr(str, "UNIFORM"))
        return ColorMapping::Uniform;
    else if (compare_str_cstr(str, "CPK"))
        return ColorMapping::Cpk;
    else if (compare_str_cstr(str, "RES_ID"))
        return ColorMapping::ResId;
    else if (compare_str_cstr(str, "RES_INDEX"))
        return ColorMapping::ResIndex;
    else if (compare_str_cstr(str, "CHAIN_ID"))
        return ColorMapping::ChainId;
    else if (compare_str_cstr(str, "CHAIN_INDEX"))
        return ColorMapping::ChainIndex;
    else if (compare_str_cstr(str, "SECONDARY_STRUCTURE"))
        return ColorMapping::SecondaryStructure;
    else
        return ColorMapping::Cpk;
}

static str_t get_color_mapping_name(ColorMapping mapping) {
    switch (mapping) {
        case ColorMapping::Uniform:
            return make_cstr("UNIFORM");
        case ColorMapping::Cpk:
            return make_cstr("CPK");
        case ColorMapping::ResId:
            return make_cstr("RES_ID");
        case ColorMapping::ResIndex:
            return make_cstr("RES_INDEX");
        case ColorMapping::ChainId:
            return make_cstr("CHAIN_ID");
        case ColorMapping::ChainIndex:
            return make_cstr("CHAIN_INDEX");
        case ColorMapping::SecondaryStructure:
            return make_cstr("SECONDARY_STRUCTURE");
        default:
            return make_cstr("UNDEFINED");
    }
}

static vec4_t parse_vec4(str_t txt, vec4_t default_val = {1,1,1,1}) {
    vec4_t res = default_val;
    str_t tok = {0};
    int i = 0;
    while (extract_next_token(&tok, &txt, ',') && i < 4) {
        res.elem[i] = parse_float(tok);
        ++i;
    }

    return res;
}

enum SerializationType {
    SerializationType_Invalid,
    SerializationType_Bool,
    SerializationType_String,
    SerializationType_Float,
    SerializationType_Double,
    SerializationType_Int8,
    SerializationType_Int16,
    SerializationType_Int32,
    SerializationType_Int64,
    SerializationType_Vec3,
    SerializationType_Vec4,
    SerializationType_Path,
    // Custom types
    SerializationType_Script,
    SerializationType_Bitfield,
};

struct SerializationObject {
    const char* group;
    const char* label;
    SerializationType type;
    size_t struct_byte_offset;
    size_t capacity = 0;
};

struct SerializationArray {
    const char* group;
    size_t array_byte_offset;
    size_t element_byte_size;
    void* (*create_item_func)(ApplicationData* data);
};

SerializationObject serialization_targets[] = {
    {"[Files]", "MoleculeFile",             SerializationType_Path,     offsetof(ApplicationData, files.molecule),    sizeof(ApplicationData::files.molecule)},
    {"[Files]", "TrajectoryFile",           SerializationType_Path,     offsetof(ApplicationData, files.trajectory),   sizeof(ApplicationData::files.trajectory)},
    
    {"[Animation]", "Frame",                SerializationType_Double,   offsetof(ApplicationData, animation.frame)},
    {"[Animation]", "Fps",                  SerializationType_Float,    offsetof(ApplicationData, animation.fps)},
    {"[Animation]", "Interpolation",        SerializationType_Int32,    offsetof(ApplicationData, animation.interpolation)},

    {"[RenderSettings]", "SsaoEnabled",     SerializationType_Bool,     offsetof(ApplicationData, visuals.ssao.enabled)},
    {"[RenderSettings]", "SsaoIntensity",   SerializationType_Float,    offsetof(ApplicationData, visuals.ssao.intensity)},
    {"[RenderSettings]", "SsaoRadius",      SerializationType_Float,    offsetof(ApplicationData, visuals.ssao.radius)},
    {"[RenderSettings]", "SsaoBias",        SerializationType_Float,    offsetof(ApplicationData, visuals.ssao.bias)},
    {"[RenderSettings]", "DofEnabled",      SerializationType_Bool,     offsetof(ApplicationData, visuals.dof.enabled)},
    {"[RenderSettings]", "DofFocusScale",   SerializationType_Bool,     offsetof(ApplicationData, visuals.dof.focus_scale)},

    {"[Camera]", "Position",                SerializationType_Vec3,     offsetof(ApplicationData, view.camera.position)},
    {"[Camera]", "Rotation",                SerializationType_Vec4,     offsetof(ApplicationData, view.camera.orientation)},
    {"[Camera]", "Distance",                SerializationType_Float,    offsetof(ApplicationData, view.camera.focus_distance)},

    {"[Representation]", "Name",            SerializationType_String,   offsetof(Representation, name),     sizeof(Representation::name)},
    {"[Representation]", "Filter",          SerializationType_String,   offsetof(Representation, filt),     sizeof(Representation::filt)},
    {"[Representation]", "Enabled",         SerializationType_Bool,     offsetof(Representation, enabled)},
    {"[Representation]", "Type",            SerializationType_Int32,    offsetof(Representation, type)},
    {"[Representation]", "ColorMapping",    SerializationType_Int32,    offsetof(Representation, color_mapping)},
    {"[Representation]", "StaticColor",     SerializationType_Vec4,     offsetof(Representation, uniform_color)},
    {"[Representation]", "Radius",          SerializationType_Float,    offsetof(Representation, radius)},
    {"[Representation]", "Tension",         SerializationType_Float,    offsetof(Representation, tension)},
    {"[Representation]", "Width",           SerializationType_Float,    offsetof(Representation, width)},
    {"[Representation]", "Thickness",       SerializationType_Float,    offsetof(Representation, thickness)},

    {"[AtomElementMapping]", "Label",       SerializationType_String,   offsetof(AtomElementMapping, lbl),  sizeof(AtomElementMapping::lbl)},
    {"[AtomElementMapping]", "Element",     SerializationType_Int8,     offsetof(AtomElementMapping, elem)},

    {"[Script]", "Text",                    SerializationType_Script,   0},

    {"[Selection]", "Label",                SerializationType_String,   offsetof(Selection, name), sizeof(Selection::name)},
    {"[Selection]", "Mask",                 SerializationType_Bitfield, offsetof(Selection, atom_mask)},
};

void* serialize_create_rep(ApplicationData* data) {
    return (void*)create_representation(data);
}

void* serialize_create_atom_elem_mapping(ApplicationData* data) {
    AtomElementMapping mapping = {};
    return md_array_push(data->dataset.atom_element_remappings, mapping, persistent_allocator);
}

void* serialize_create_selection(ApplicationData* data) {
    Selection sel = {};
    md_bitfield_init(&sel.atom_mask, persistent_allocator);
    return md_array_push(data->selection.stored_selections, sel, persistent_allocator);
}

SerializationArray serialization_array_groups[] = {
    {"[Representation]",        offsetof(ApplicationData, representations.buffer),          sizeof(Representation),     serialize_create_rep},
    {"[AtomElementMapping]",    offsetof(ApplicationData, dataset.atom_element_remappings), sizeof(AtomElementMapping), serialize_create_atom_elem_mapping},
    {"[Selection]",             offsetof(ApplicationData, selection.stored_selections),     sizeof(Selection),          serialize_create_selection},
};

#define COMPARE(str, ref) (compare_str_cstr_n(str, ref"", sizeof(ref) - 1))
#define EXTRACT(str, ref) (compare_str_cstr_n(str, ref"", sizeof(ref) - 1) && (line = trim_whitespace(substr(line, sizeof(ref) - 1))).len > 0)
#define EXTRACT_PARAM_LINE(line, txt) (c_txt.len && c_txt[0] != '[' && (extract_line(&line, &c_txt)))

static const SerializationObject* find_serialization_target(str_t group, str_t label) {
    for (int64_t i = 0; i < (int64_t)ARRAY_SIZE(serialization_targets); ++i) {
        if (compare_str_cstr(group, serialization_targets[i].group) && compare_str_cstr(label, serialization_targets[i].label)) {
            return &serialization_targets[i];
        }
    }
    return NULL;
}

static const SerializationArray* find_serialization_array_group(str_t group) {
    for (int64_t i = 0; i < ARRAY_SIZE(serialization_array_groups); ++i) {
        if (compare_str_cstr(group, serialization_array_groups[i].group)) {
            return &serialization_array_groups[i];
        }
    }
    return NULL;
}

static void deserialize_object(const SerializationObject* target, char* ptr, str_t* buf, str_t filename) {
    str_t line;
    if (extract_line(&line, buf)) {
        str_t arg = trim_whitespace(line);

        switch (target->type) {
        case SerializationType_Bool:
        {
            bool value = parse_int(arg) != 0;
            *(bool*)(ptr + target->struct_byte_offset) = value;
            break;
        }
        case SerializationType_String:
        {
            size_t copy_len = MIN(target->capacity - 1, (size_t)arg.len);
            memcpy(ptr + target->struct_byte_offset, arg.ptr, copy_len);
            (ptr + target->struct_byte_offset)[copy_len] = '\0';
            break;
        }
        case SerializationType_Float:
        {
            float value = (float)parse_float(arg);
            *(float*)(ptr + target->struct_byte_offset) = value;
            break;
        }
        case SerializationType_Double:
        {
            double value = parse_float(arg);
            *(double*)(ptr + target->struct_byte_offset) = value;
            break;
        }
        case SerializationType_Int8:
        {
            int8_t value = (int8_t)parse_int(arg);
            *(int8_t*)(ptr + target->struct_byte_offset) = value;
            break;
        }
        case SerializationType_Int16:
        {
            int16_t value = (int16_t)parse_int(arg);
            *(int16_t*)(ptr + target->struct_byte_offset) = value;
            break;
        }
        case SerializationType_Int32:
        {
            int32_t value = (int32_t)parse_int(arg);
            *(int32_t*)(ptr + target->struct_byte_offset) = value;
            break;
        }
        case SerializationType_Int64:
        {
            int64_t value = (int64_t)parse_int(arg);
            *(int64_t*)(ptr + target->struct_byte_offset) = value;
            break;
        }
        case SerializationType_Vec3:
        {
            vec3_t value = vec3_from_vec4(parse_vec4(arg));
            *(vec3_t*)(ptr + target->struct_byte_offset) = value;
            break;
        }
        case SerializationType_Vec4:
        {
            vec4_t value = parse_vec4(arg);
            *(vec4_t*)(ptr + target->struct_byte_offset) = value;
            break;
        }
        case SerializationType_Path:
        {
            StrBuf<512> path = extract_path_without_file(filename);
            path += arg;
            str_t can_path = md_os_path_make_canonical(path, frame_allocator);
            if (can_path.ptr && can_path.len > 0) {
                size_t copy_len = MIN(target->capacity - 1, (size_t)can_path.len);
                memcpy(ptr + target->struct_byte_offset, can_path.ptr, copy_len);
                (ptr + target->struct_byte_offset)[copy_len] = '\0';
            }
            break;
        }
        case SerializationType_Script:
        {
            // Script starts with """
            // and ends with """
            str_t token = make_cstr("\"\"\"");
            if (compare_str_n(arg, token, token.len)) {
                // Roll back buf to arg + 3
                const char* beg = arg.ptr + token.len;
                buf->len = buf->end() - beg;
                buf->ptr = beg;
                const char* end = str_find_str(*buf, token).ptr;
                if (end) {
                    std::string str(beg, end - beg);
                    ApplicationData* data = (ApplicationData*)ptr;
                    data->script.editor.SetText(str);
                    // Set buf pointer to after token
                    const char* pos = end + token.len;
                    buf->len = buf->end() - pos;
                    buf->ptr = pos;
                    break;
                } else {
                    md_print(MD_LOG_TYPE_ERROR, "Malformed end token for script");
                    return;
                }
            } else {
                md_print(MD_LOG_TYPE_ERROR, "Malformed start token for script");
                return;
            }
        }
        case SerializationType_Bitfield:
        {
            // Bitfield starts with ###
            // and ends with ###
            str_t token = make_cstr("###");
            if (compare_str_n(arg, token, token.len)) {
                // Roll back buf to arg + 3
                const char* beg = arg.ptr + token.len;
                buf->len = buf->end() - beg;
                buf->ptr = beg;
                const char* end = str_find_str(*buf, token).ptr;
                if (end) {
                    void* base64_data = md_alloc(frame_allocator, md_base64_decode_size_in_bytes(end-beg));
                    int64_t base64_size = md_base64_decode(base64_data, beg, end-beg);
                    md_exp_bitfield_t* bf = (md_exp_bitfield_t*)(ptr + target->struct_byte_offset);
                    if (!base64_size || !md_bitfield_deserialize(bf, base64_data, base64_size)) {
                        md_print(MD_LOG_TYPE_ERROR, "Failed to deserialize bitfield");
                        md_bitfield_clear(bf);
                        return;
                    }
                    // Set buf pointer to after token
                    const char* pos = end + token.len;
                    buf->len = buf->end() - pos;
                    buf->ptr = pos;
                    break;
                } else {
                    md_print(MD_LOG_TYPE_ERROR, "Malformed end token for bitfield");
                    return;
                }
            } else {
                md_print(MD_LOG_TYPE_ERROR, "Malformed start token for bitfield");
                return;
            }
        }
        case SerializationType_Invalid: // fallthrough
        default:
            ASSERT(false);
        }
    }
}

static void load_workspace(ApplicationData* data, str_t filename) {
    str_t txt = load_textfile(filename, frame_allocator);
    defer { free_str(txt, frame_allocator); };

    if (!txt.len) {
        md_printf(MD_LOG_TYPE_ERROR, "Could not open workspace file: '%.*s", (int)filename.len, filename.ptr);
        return;
    }

    // Reset and clear things
    clear_representations(data);
    data->script.editor.SetText("");

    data->animation = {};
    reset_view(data, false, true);

    str_t group = {0};
    str_t c_txt = txt;
    str_t line = {};

    str_t cur_molecule_file = copy_str(data->files.molecule, frame_allocator);
    str_t cur_trajectory_file = copy_str(data->files.trajectory, frame_allocator);

    const SerializationArray* arr_group = NULL;
    void* ptr = 0;

    while (extract_line(&line, &c_txt)) {
        line = trim_whitespace(line);
        if (line[0] == '[') {
            group = line;
            arr_group = find_serialization_array_group(group);
            if (arr_group) {
                ptr = arr_group->create_item_func(data);
            } else {
                ptr = data;
            }
        } else {
            int64_t loc = find_char(line, '=');
            if (loc != -1) {
                str_t label = trim_whitespace(substr(line, 0, loc));
                const SerializationObject* target = find_serialization_target(group, label);
                if (target) {
                    // Move read pointer back 'loc+1' => after '='.
                    const char* pos = line.ptr + loc + 1;
                    c_txt.len = c_txt.end() - pos;
                    c_txt.ptr = pos;
                    deserialize_object(target, (char*)ptr, &c_txt, filename);
                } else {
                    md_printf(MD_LOG_TYPE_ERROR, "Could not recognize serialization target '%.*s' in group '%.*s,", (int)label.len, label.ptr, (int)group.len, group.ptr);
                }
            }
        }
    }
    data->view.animation.target_position = data->view.camera.position;
    data->files.workspace = filename;

    str_t new_molecule_file = copy_str(data->files.molecule, frame_allocator);
    str_t new_trajectory_file = copy_str(data->files.trajectory, frame_allocator);

    // When we de-serialize we overwrite the two following paths, even though they are not loaded.
    // So we copy them back to their original values.
    data->files.molecule = cur_molecule_file;
    data->files.trajectory = cur_trajectory_file;

    if (new_molecule_file.len && load_dataset_from_file(data, new_molecule_file)) {
        init_all_representations(data);
        update_all_representations(data);
    }

    if (new_trajectory_file.len) {
        load_dataset_from_file(data, new_trajectory_file);
    }

    apply_atom_elem_mappings(data);
}

static void write_entry(FILE* file, SerializationObject target, const void* ptr, str_t filename) {
    fprintf(file, "%s=", target.label);

    switch (target.type) {
    case SerializationType_Bool:
    {
        bool value = *(bool*)((const char*)ptr + target.struct_byte_offset) != 0;
        fprintf(file, "%i\n", value ? 1 : 0);
        break;
    }
    case SerializationType_String:
    {
        const char* str = (const char*)((const char*)ptr + target.struct_byte_offset);
        int len = (int)strnlen(str, target.capacity);
        fprintf(file, "%.*s\n", len, str);
        break;
    }
    case SerializationType_Float:
    {
        float value = *(float*)((const char*)ptr + target.struct_byte_offset);
        fprintf(file, "%g\n", value);
        break;
    }
    case SerializationType_Double:
    {
        double value = *(double*)((const char*)ptr + target.struct_byte_offset);
        fprintf(file, "%g\n", value);
        break;
    }
    case SerializationType_Int8:
    {
        int32_t value = *(int8_t*)((const char*)ptr + target.struct_byte_offset);
        fprintf(file, "%i\n", value);
        break;
    }
    case SerializationType_Int16:
    {
        int32_t value = *(int16_t*)((const char*)ptr + target.struct_byte_offset);
        fprintf(file, "%i\n", value);
        break;
    }
    case SerializationType_Int32:
    {
        int32_t value = *(int32_t*)((const char*)ptr + target.struct_byte_offset);
        fprintf(file, "%i\n", value);
        break;
    }
    case SerializationType_Int64:
    {
        int32_t value = *(int64_t*)((const char*)ptr + target.struct_byte_offset);
        fprintf(file, "%i\n", value);
        break;
    }
    case SerializationType_Vec3:
    {
        vec3_t value = *(vec3_t*)((const char*)ptr + target.struct_byte_offset);
        fprintf(file, "%g,%g,%g\n", value.x, value.y, value.z);
        break;
    }
    case SerializationType_Vec4:
    {
        vec4_t value = *(vec4_t*)((const char*)ptr + target.struct_byte_offset);
        fprintf(file, "%g,%g,%g,%g\n", value.x, value.y, value.z, value.w);
        break;
    }
    case SerializationType_Path:
    {
        const char* str = (const char*)((const char*)ptr + target.struct_byte_offset);
        int len = (int)strnlen(str, target.capacity);

        // Make this sucker relative
        str_t rel_path = md_os_path_make_relative(filename, {str, len}, frame_allocator);
        if (rel_path.ptr && rel_path.len) {
            fprintf(file, "%.*s\n", (int)rel_path.len, rel_path.ptr);
        }     
        break;
    }
    case SerializationType_Script:
    {
        ApplicationData* data = (ApplicationData*)ptr;
        std::string str = data->script.editor.GetText();
        fprintf(file, "\"\"\"%s\"\"\"\n", str.c_str());
        break;
    }
    case SerializationType_Bitfield:
    {
        const md_exp_bitfield_t* bf = (const md_exp_bitfield_t*)((const char*)ptr + target.struct_byte_offset);
        void* serialized_data = md_alloc(frame_allocator, md_bitfield_serialize_size_in_bytes(bf));
        int64_t serialized_size = md_bitfield_serialize(serialized_data, bf);
        if (serialized_size) {
            char* base64_data = (char*)md_alloc(frame_allocator, md_base64_encode_size_in_bytes(serialized_size));
            int64_t base64_size = md_base64_encode(base64_data, serialized_data, serialized_size);
            if (base64_size) {
                fprintf(file, "###%.*s###\n", (int)base64_size, base64_data);
            }
        }
        break;
    }
    case SerializationType_Invalid: // fallthrough
    default:
        ASSERT(false);
    }
}

static void save_workspace(ApplicationData* data, str_t filename) {
    md_file_o* file = md_file_open(filename, MD_FILE_WRITE);

    if (!file) {
        md_printf(MD_LOG_TYPE_ERROR, "Could not open workspace file: '%.*s", (int)filename.len, filename.ptr);
        return;
    }
    defer { md_file_close(file); };

    const char* curr_group = "";

    for (int64_t i = 0; i < (int64_t)ARRAY_SIZE(serialization_targets); ++i) {
        const char* group = serialization_targets[i].group;

        const SerializationArray* arr_group = find_serialization_array_group({group, (int64_t)strlen(group)});
        if (arr_group) {
            // Special case for this since it is an array quantity, iterate over all array items then all serialization subfields marked with group
            const void* arr = *((const void**)((char*)data + arr_group->array_byte_offset));
            const int64_t arr_size = md_array_size(arr);

            if (arr_size) {
                int64_t beg_field_i = i;
                int64_t end_field_i = i + 1;
                while (end_field_i < (int64_t)ARRAY_SIZE(serialization_targets) && strcmp(serialization_targets[end_field_i].group, group) == 0) {
                    end_field_i += 1;
                }
                // Iterate over all array elements
                for (int64_t arr_idx = 0; arr_idx < arr_size; ++arr_idx) {
                    // Write group
                    fprintf((FILE*)file, "\n%s\n", group);

                    // Write all fields for item
                    const void* item_ptr = (const char*)arr + arr_idx * arr_group->element_byte_size;
                    for (int64_t j = beg_field_i; j < end_field_i; ++j) {
                        write_entry((FILE*)file, serialization_targets[j], item_ptr, filename);
                    }
                }
                i = end_field_i - 1;
            }

        }
        else {
            if (strcmp(group, curr_group) != 0) {
                fprintf((FILE*)file, "\n%s\n", group);
                curr_group = group;
            }

            write_entry((FILE*)file, serialization_targets[i], data, filename);
        }
    }
}

void create_screenshot(ApplicationData* data) {
    ASSERT(data);
    image_t img;
    init_image(&img, data->gbuffer.width, data->gbuffer.height, frame_allocator);
    defer { free_image(&img, frame_allocator); };

    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    glReadBuffer(GL_BACK);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    glReadPixels(0, 0, img.width, img.height, GL_RGBA, GL_UNSIGNED_BYTE, img.data);

    {
        // @NOTE: Swap Rows to flip image with respect to y-axis
        const uint32_t row_byte_size = img.width * sizeof(uint32_t);
        uint32_t* row_t = (uint32_t*)md_alloc(frame_allocator, row_byte_size);
        defer { md_free(frame_allocator, row_t, row_byte_size); };
        for (uint32_t i = 0; i < (uint32_t)img.height / 2; ++i) {
            uint32_t* row_a = img.data + i * img.width;
            uint32_t* row_b = img.data + (img.height - 1 - i) * img.width;
            if (row_a != row_b) {
                memcpy(row_t, row_a, row_byte_size);  // tmp = a;
                memcpy(row_a, row_b, row_byte_size);  // a = b;
                memcpy(row_b, row_t, row_byte_size);  // b = tmp;
            }
        }
    }

    application::FileDialogResult file_res = application::file_dialog(application::FileDialogFlags_Save, {}, make_cstr("jpg;png;bmp"));
    if (file_res.result == application::FileDialogResult::Ok) {
        str_t ext = extract_ext({file_res.path, file_res.path_len});
        if (ext.ptr == NULL) {
            snprintf(file_res.path + file_res.path_len, ARRAY_SIZE(file_res.path) - file_res.path_len, ".jpg");
            ext = make_cstr("jpg");
        }
        if (compare_str_cstr_ignore_case(ext, "jpg")) {
            const int quality = 95;
            write_image_jpg(img, file_res.path, quality);
        } else if (compare_str_cstr_ignore_case(ext, "png")) {
            write_image_png(img, file_res.path);
        } else if (compare_str_cstr_ignore_case(ext, "bmp")) {
            write_image_bmp(img, file_res.path);
        } else {
            md_print(MD_LOG_TYPE_ERROR, "Supplied image format is not supported");
        }
    }
}

// #representation
static Representation* create_representation(ApplicationData* data, RepresentationType type, ColorMapping color_mapping, str_t filter) {
    ASSERT(data);
    Representation rep;
    rep.type = type;
    rep.color_mapping = color_mapping;
    rep.filt = filter;   
    init_representation(data, &rep);
    update_representation(data, &rep);
    return md_array_push(data->representations.buffer, rep, persistent_allocator);
}

static Representation* clone_representation(ApplicationData* data, const Representation& rep) {
    ASSERT(data);
    Representation* clone = md_array_push(data->representations.buffer, rep, persistent_allocator);
    clone->md_rep = {0};
    clone->atom_mask = {0};
    init_representation(data, clone);
    update_representation(data, clone);
    return clone;
}

static void remove_representation(ApplicationData* data, int idx) {
    ASSERT(data);
    ASSERT(idx < md_array_size(data->representations.buffer));
    auto& rep = data->representations.buffer[idx];
    md_bitfield_free(&rep.atom_mask);
    md_gl_representation_free(&rep.md_rep);
    data->representations.buffer[idx] = *md_array_last(data->representations.buffer);
    md_array_pop(data->representations.buffer);
}

static void recompute_atom_visibility_mask(ApplicationData* data) {
    ASSERT(data);

    auto& mask = data->representations.atom_visibility_mask;

    md_bitfield_clear(&mask);
    for (int64_t i = 0; i < md_array_size(data->representations.buffer); ++i) {
        auto& rep = data->representations.buffer[i];
        if (!rep.enabled) continue;
        md_bitfield_or_inplace(&mask, &rep.atom_mask);
    }

    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

static void update_all_representations(ApplicationData* data) {
    for (int64_t i = 0; i < md_array_size(data->representations.buffer); ++i) {
        auto& rep = data->representations.buffer[i];
        update_representation(data, &rep);
    }
}

static void update_representation(ApplicationData* data, Representation* rep) {
    ASSERT(data);
    ASSERT(rep);

    const int64_t bytes = data->mold.mol.atom.count * sizeof(uint32_t);
    uint32_t* colors = (uint32_t*)md_alloc(frame_allocator, bytes);
    defer { md_free(frame_allocator, colors, bytes); };

    const auto& mol = data->mold.mol;

    //md_script_property_t prop = {0};
    //if (rep->color_mapping == ColorMapping::Property) {
        //rep->prop_is_valid = md_script_compile_and_eval_property(&prop, rep->prop, &data->mold.mol, frame_allocator, &data->mold.script.ir, rep->prop_error.beg(), rep->prop_error.capacity());
    //}

    switch (rep->color_mapping) {
        case ColorMapping::Uniform:
            color_atoms_uniform(colors, mol.atom.count, rep->uniform_color);
            break;
        case ColorMapping::Cpk:
            color_atoms_cpk(colors, mol.atom.count, mol);
            break;
        case ColorMapping::ResId:
            color_atoms_residue_id(colors, mol.atom.count, mol);
            break;
        case ColorMapping::ResIndex:
            color_atoms_residue_index(colors, mol.atom.count, mol);
            break;
        case ColorMapping::ChainId:
            color_atoms_chain_id(colors, mol.atom.count, mol);
            break;
        case ColorMapping::ChainIndex:
            color_atoms_chain_index(colors, mol.atom.count, mol);
            break;
        case ColorMapping::SecondaryStructure:
            color_atoms_secondary_structure(colors, mol.atom.count, mol);
            break;
        case ColorMapping::Property:
            // @TODO: Map colors accordingly
            //color_atoms_uniform(colors, mol.atom.count, rep->uniform_color);
            if (rep->prop) {
                color_atoms_uniform(colors, mol.atom.count, {1,1,1,1});
                const float* values = rep->prop->data.values;
                const int num_values = (int)rep->prop->data.num_values;
                if (rep->prop->data.aggregate) {
                    const int dim = rep->prop->data.dim[0];
                    md_script_visualization_args_t args = {
                        .token = rep->prop->vis_token,
                        .ir = &data->mold.script.ir,
                        .mol = &data->mold.mol,
                        .traj = NULL,
                        .alloc = frame_allocator,
                        .flags = MD_SCRIPT_VISUALIZE_ATOMS
                    };
                    md_script_visualization_t vis = {0};
                    if (md_script_visualization_init(&vis, args)) {
                        if (dim == (int)md_array_size(vis.structures.atom_masks)) {
                            int i0 = CLAMP((int)data->animation.frame + 0, 0, rep->prop->data.num_values / dim - 1);
                            int i1 = CLAMP((int)data->animation.frame + 1, 0, rep->prop->data.num_values / dim - 1);
                            float frame_fract = fractf((float)data->animation.frame);

                            md_exp_bitfield_t mask = {0};
                            md_bitfield_init(&mask, frame_allocator);
                            for (int i = 0; i < dim; ++i) {
                                md_bitfield_and(&mask, &rep->atom_mask, &vis.structures.atom_masks[i]);
                                float value = lerpf(values[i0 * dim + i], values[i1 * dim + i], frame_fract);
                                float t = CLAMP((value - rep->map_beg) / (rep->map_end - rep->map_beg), 0, 1);
                                ImVec4 color = ImPlot::SampleColormap(t, rep->color_map);
                                color_atoms_uniform(colors, mol.atom.count, vec_cast(color), &mask);
                            }
                        }
                    }
                } else {
                    int i0 = CLAMP((int)data->animation.frame + 0, 0, rep->prop->data.num_values - 1);
                    int i1 = CLAMP((int)data->animation.frame + 1, 0, rep->prop->data.num_values - 1);
                    float value = lerpf(values[i0], values[i1], fractf((float)data->animation.frame));
                    float t = CLAMP((value - rep->map_beg) / (rep->map_end - rep->map_beg), 0, 1);
                    ImVec4 color = ImPlot::SampleColormap(t, rep->color_map);
                    color_atoms_uniform(colors, mol.atom.count, vec_cast(color));
                }
            } else {
                color_atoms_uniform(colors, mol.atom.count, rep->uniform_color);
            }
            break;
        default:
            ASSERT(false);
            break;
    }

    rep->filt_is_valid = filter_expression(data, rep->filt, &rep->atom_mask, &rep->filt_is_dynamic, rep->filt_error.beg(), rep->filt_error.capacity());

    if (rep->filt_is_valid) {
        filter_colors(colors, mol.atom.count, &rep->atom_mask);
        data->representations.atom_visibility_mask_dirty = true;

        md_gl_representation_type_t type = MD_GL_REP_DEFAULT;
        md_gl_representation_args_t args = {};
        switch(rep->type) {
        case RepresentationType::SpaceFill:
            type = MD_GL_REP_SPACE_FILL;
            args.space_fill.radius_scale = rep->radius;
            break;
        case RepresentationType::Licorice:
            type = MD_GL_REP_LICORICE;
            args.licorice.radius = rep->radius * 0.5f;
            break;
        case RepresentationType::Ribbons:
            type = MD_GL_REP_RIBBONS;
            args.ribbons.width_scale = rep->width;
            args.ribbons.thickness_scale = rep->thickness;
            break;
        case RepresentationType::Cartoon:
            type = MD_GL_REP_CARTOON;
            args.cartoon.width_scale = rep->width;
            args.cartoon.thickness_scale = rep->thickness;
            break;
        default: break;
        }

        md_gl_representation_set_type_and_args(&rep->md_rep, type, args);
        md_gl_representation_set_color(&rep->md_rep, 0, (uint32_t)mol.atom.count, colors, 0);
    }
}

static void init_representation(ApplicationData* data, Representation* rep) {
    md_gl_representation_init(&rep->md_rep, &data->mold.gl_mol);
    md_bitfield_init(&rep->atom_mask, persistent_allocator);
}

static void init_all_representations(ApplicationData* data) {
    for (int64_t i = 0; i < md_array_size(data->representations.buffer); ++i) {
        auto& rep = data->representations.buffer[i];
        init_representation(data, &rep);
    }
}

static void clear_representations(ApplicationData* data) {
    ASSERT(data);
    while (md_array_size(data->representations.buffer) > 0) {
        remove_representation(data, (int32_t)md_array_size(data->representations.buffer) - 1);
    }
}

// #selection
static Selection* create_selection(ApplicationData* data, str_t name, md_exp_bitfield_t* atom_mask) {
    ASSERT(data);
    Selection sel;
    sel.name = name;
    md_bitfield_init(&sel.atom_mask, persistent_allocator);
    if (atom_mask) {
        md_bitfield_copy(&sel.atom_mask, atom_mask);
    }
    return md_array_push(data->selection.stored_selections, sel, persistent_allocator);
}

#if 0
    static Selection* clone_selection(ApplicationData * data, const Selection& src) {
        ASSERT(data);
        Selection clone;
        clone.name = src.name;
        bitfield::init(&clone.atom_mask, data->selection.current_selection_mask);
        return &data->selection.stored_selections.push_back(clone);
    }
#endif

static void remove_selection(ApplicationData* data, int idx) {
    ASSERT(data);
    if (idx < 0 || (int)md_array_size(data->selection.stored_selections) <= idx) {
        md_printf(MD_LOG_TYPE_ERROR, "Index [%i] out of range when trying to remove selection", idx);
    }
    auto item = &data->selection.stored_selections[idx];
    md_bitfield_free(&item->atom_mask);
    
    data->selection.stored_selections[idx] = *md_array_last(data->selection.stored_selections);
    md_array_pop(data->selection.stored_selections);
}

#if 0
static void reset_selections(ApplicationData* data) {
    UNUSED(data);
    // ASSERT(data);
    // @NOTE: What to do here?
}

static void clear_selections(ApplicationData* data) {
    ASSERT(data);
    while (data->selection.stored_selections.size() > 0) {
        remove_selection(data, (int32)data->selection.stored_selections.size() - 1);
    }
}
#endif

static bool handle_selection(ApplicationData* data) {
    ASSERT(data);
    enum class RegionMode { Append, Remove };

    const int64_t N = data->mold.mol.atom.count;
    if (N <= 0) return false;

    static RegionMode region_mode = RegionMode::Append;
    static bool region_select = false;
    static application::Coordinate x0;
    const application::Coordinate x1 = data->ctx.input.mouse.win_coord;
    const bool shift_down = ImGui::GetIO().KeyShift;
    const bool mouse_down = data->ctx.input.mouse.down[0] || data->ctx.input.mouse.down[1];

    md_exp_bitfield_t mask;
    md_bitfield_init(&mask, frame_allocator);
    defer { md_bitfield_free(&mask); };

    md_bitfield_clear(&data->selection.current_highlight_mask);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;

    if (data->picking.idx != INVALID_PICKING_IDX && !region_select) {
        ASSERT(0 <= data->picking.idx && data->picking.idx <= N);
        md_bitfield_set_bit(&mask, data->picking.idx);

        switch (data->selection.granularity) {
            case SelectionLevel::Atom:
                break;
            case SelectionLevel::Residue: {
                expand_mask(&mask, data->mold.mol.residue.atom_range, data->mold.mol.residue.count);
                break;
            }
            case SelectionLevel::Chain: {
                expand_mask(&mask, data->mold.mol.chain.atom_range, data->mold.mol.chain.count);
                break;
            }
            default:
                ASSERT(false);
                break;
        }
        md_bitfield_copy(&data->selection.current_highlight_mask, &mask);
        data->mold.dirty_buffers |= MolBit_DirtyFlags;
    }

    if (shift_down) {
        if (data->ctx.input.mouse.hit[0] || data->ctx.input.mouse.hit[1]) {
            x0 = data->ctx.input.mouse.win_coord;
            region_mode = data->ctx.input.mouse.hit[0] ? RegionMode::Append : RegionMode::Remove;
        }

        if (mouse_down && x1 != x0) {
            region_select = true;
        }

        const ImVec2 min_p = ImVec2(MIN(x0.x, x1.x), MIN(x0.y, x1.y));
        const ImVec2 max_p = ImVec2(MAX(x0.x, x1.x), MAX(x0.y, x1.y));

        if (region_select) {
            const vec2_t res = {(float)data->ctx.window.width, (float)data->ctx.window.height};
            const mat4_t mvp = data->view.param.matrix.current.view_proj;
            const md_exp_bitfield_t* vis_mask = &data->representations.atom_visibility_mask;

            int64_t beg_bit = vis_mask->beg_bit;
            int64_t end_bit = vis_mask->end_bit;
            while ((beg_bit = md_bitfield_scan(vis_mask, beg_bit, end_bit)) != 0) {
                int64_t i = beg_bit - 1;
                vec4_t p =  {data->mold.mol.atom.x[i], data->mold.mol.atom.y[i], data->mold.mol.atom.z[i], 1.0f};
                p = mat4_mul_vec4(mvp, p);

                vec2_t c = {
                    (p.x / p.w * 0.5f + 0.5f) * res.x,
                    (-p.y / p.w * 0.5f + 0.5f) * res.y
                };

                if (min_p.x <= c.x && c.x <= max_p.x && min_p.y <= c.y && c.y <= max_p.y) {
                    md_bitfield_set_bit(&mask, i);
                }
            }

            switch (data->selection.granularity) {
                case SelectionLevel::Atom:
                    break;
                case SelectionLevel::Residue:
                    expand_mask(&mask, data->mold.mol.residue.atom_range, data->mold.mol.residue.count);
                    break;
                case SelectionLevel::Chain:
                    expand_mask(&mask, data->mold.mol.chain.atom_range, data->mold.mol.residue.count);
                    break;
                default:
                    ASSERT(false);
            }

            md_exp_bitfield_t* dst_mask = mouse_down ? &data->selection.current_highlight_mask : &data->selection.current_selection_mask;
            md_exp_bitfield_t* src_mask = &data->selection.current_selection_mask;

            if (region_mode == RegionMode::Append) {
                if (dst_mask == src_mask) {
                    md_bitfield_or_inplace(dst_mask, &mask);
                } else {
                    md_bitfield_or(dst_mask, src_mask, &mask);
                }
            } else if (region_mode == RegionMode::Remove) {
                md_bitfield_not_inplace(&mask, 0, N);
                if (dst_mask == src_mask) {
                    md_bitfield_and_inplace(dst_mask, &mask);
                } else {
                    md_bitfield_and(dst_mask, src_mask, &mask);
                }
            }
            data->mold.dirty_buffers |= MolBit_DirtyFlags;

            if (!mouse_down) region_select = false;

            // Draw selection window
            const ImVec2 vp_pos = ImGui::GetMainViewport()->Pos;
            const ImVec2 pos0 = (ImVec2(min_p.x, min_p.y) + vp_pos);
            const ImVec2 pos1 = (ImVec2(max_p.x, max_p.y) + vp_pos);
            const ImU32 fill_col = 0x22222222;
            const ImU32 line_col = 0x88888888;

            ImGui::BeginCanvas("region select");
            auto dl = ImGui::GetCurrentWindow()->DrawList;
            dl->AddRectFilled(pos0, pos1, fill_col);
            dl->AddRect(pos0, pos1, line_col);
            ImGui::EndCanvas();
        } else if (data->ctx.input.mouse.clicked[0] || data->ctx.input.mouse.clicked[1]) {
            if (data->picking.idx != INVALID_PICKING_IDX) {
                const bool append = data->ctx.input.mouse.clicked[0];
                if (append) {
                    if (data->selection.granularity == SelectionLevel::Atom && data->picking.idx != ~0U) {
                        single_selection_sequence_push_back(&data->selection.single_selection_sequence, data->picking.idx);
                    }
                    md_bitfield_or_inplace(&data->selection.current_selection_mask, &mask);
                } else {
                    md_bitfield_not_inplace(&mask, 0, N);
                    md_bitfield_and_inplace(&data->selection.current_selection_mask, &mask);

                    if (data->selection.granularity == SelectionLevel::Atom && data->picking.idx != ~0U) {
                        if (single_selection_sequence_last(&data->selection.single_selection_sequence) == (int32_t)data->picking.idx) {
                            single_selection_sequence_pop_back(&data->selection.single_selection_sequence);
                        } else {
                            single_selection_sequence_clear(&data->selection.single_selection_sequence);
                        }
                    }
                }
            } else if (data->ctx.input.mouse.clicked[1]) {
                md_bitfield_clear(&data->selection.current_selection_mask);
                single_selection_sequence_clear(&data->selection.single_selection_sequence);
            }
        }
        return true;
    } else {
        region_select = false;
    }

    return false;
}

// #camera-control
static void handle_camera_interaction(ApplicationData* data) {
    if (!ImGui::GetIO().WantCaptureMouse && !data->selection.selecting) {
        const vec2_t mouse_delta = {data->ctx.input.mouse.win_delta.x, data->ctx.input.mouse.win_delta.y};
        TrackballControllerInput input;
        input.rotate_button = data->ctx.input.mouse.down[0];
        input.pan_button = data->ctx.input.mouse.down[1];
        input.dolly_button = data->ctx.input.mouse.down[2];
        input.mouse_coord_curr = {data->ctx.input.mouse.win_coord.x, data->ctx.input.mouse.win_coord.y};
        input.mouse_coord_prev = input.mouse_coord_curr - mouse_delta;
        input.screen_size = {(float)data->ctx.window.width, (float)data->ctx.window.height};
        input.dolly_delta = data->ctx.input.mouse.scroll_delta;
        input.fov_y = data->view.camera.fov_y;

        if (camera_controller_trackball(&data->view.camera.position, &data->view.camera.orientation, &data->view.camera.focus_distance, input, data->view.trackball_param)) {
            data->view.animation.target_position = data->view.camera.position;
        }

        if (ImGui::GetIO().MouseDoubleClicked[0]) {
            if (data->picking.depth < 1.0f) {
                const vec3_t forward = data->view.camera.orientation * vec3_t{0, 0, 1};
                const float dist = data->view.camera.focus_distance;
                data->view.animation.target_position = data->picking.world_coord + forward * dist;
            }
        }

        data->visuals.dof.focus_depth.target = data->view.camera.focus_distance;
    }
}

static void handle_camera_animation(ApplicationData* data) {
    const float dt = (float)MIN(data->ctx.timing.delta_s, 0.033);
    {
        // #camera-translation
        constexpr float speed = 10.0f;
        const vec3_t vel = (data->view.animation.target_position - data->view.camera.position) * speed;
        data->view.camera.position = data->view.camera.position + vel * dt;
    }
    {
        //data->view.camera.orientation = math::slerp(data->view.camera.orientation, data->view.animation.target_orientation, 0.99f);
    }
    {
        // #focus-depth
        constexpr float speed = 10.0f;
        const float vel = (data->visuals.dof.focus_depth.target - data->visuals.dof.focus_depth.current) * speed;
        data->visuals.dof.focus_depth.current += vel * dt;
    }
#if 0
        ImGui::Begin("Camera Debug Info");
        ImGui::Text("lin vel [%.2f %.2f %.2f]", vel.x, vel.y, vel.z);
        ImGui::Text("lin cur [%.2f %.2f %.2f]", data->view.camera.position.x, data->view.camera.position.y, data->view.camera.position.z);
        ImGui::Text("lin tar [%.2f %.2f %.2f]", data->view.animation.target_position.x, data->view.animation.target_position.y, data->view.animation.target_position.z);
        ImGui::Text("ang cur [%.2f %.2f %.2f %.2f]", data->view.camera.orientation.x, data->view.camera.orientation.y, data->view.camera.orientation.z, data->view.camera.orientation.w);
        ImGui::End();
#endif
}

#if EXPERIMENTAL_CONE_TRACED_AO == 1
static void init_occupancy_volume(ApplicationData* data) {
    const AABB box = compute_aabb(data->mold.mol.atom.position, data->mold.mol.atom.radius, data->mold.mol.atom.count);
    cone_trace::init_occlusion_volume(&data->occupancy_volume.vol, box.min, box.max, 8.0f);

    //data->occupancy_volume.vol_mutex.lock();
    if (data->occupancy_volume.vol.texture_id) {
        cone_trace::compute_occupancy_volume(data->occupancy_volume.vol, data->mold.mol.atom.position,
                                             data->mold.mol.atom.radius, data->representations.atom_visibility_mask);
    }
    //data->occupancy_volume.vol_mutex.unlock();
}
#endif

static void clear_gbuffer(GBuffer* gbuffer) {
    const vec4_t CLEAR_INDEX = vec4_t{1, 1, 1, 1};
    const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
        GL_COLOR_ATTACHMENT_POST_TONEMAP, GL_COLOR_ATTACHMENT_PICKING};

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuffer->deferred.fbo);
    glViewport(0, 0, gbuffer->width, gbuffer->height);

    glDepthMask(1);
    glColorMask(1, 1, 1, 1);

    // Setup gbuffer and clear textures
    PUSH_GPU_SECTION("Clear G-buffer") {
        // Clear color+alpha, normal, velocity, emissive, post_tonemap and depth
        glDrawBuffers(ARRAY_SIZE(draw_buffers), draw_buffers);
        glClearColor(0, 0, 0, 0);
        glClearDepthf(1.f);
        glStencilMask(0xFF);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

        // Clear picking buffer
        glDrawBuffer(GL_COLOR_ATTACHMENT_PICKING);
        glClearColor(CLEAR_INDEX.x, CLEAR_INDEX.y, CLEAR_INDEX.z, CLEAR_INDEX.w);
        glClear(GL_COLOR_BUFFER_BIT);
    }
    POP_GPU_SECTION()
}

static void fill_gbuffer(ApplicationData* data) {
    const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
        GL_COLOR_ATTACHMENT_POST_TONEMAP, GL_COLOR_ATTACHMENT_PICKING};

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    // Enable all draw buffers
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, data->gbuffer.deferred.fbo);
    glDrawBuffers(ARRAY_SIZE(draw_buffers), draw_buffers);

    PUSH_GPU_SECTION("G-Buffer fill")

    // SIMULATION BOX
    if (data->simulation_box.enabled && data->simulation_box.box != mat3_t{0}) {
        PUSH_GPU_SECTION("Draw Simulation Box")

        immediate::set_model_view_matrix(data->view.param.matrix.current.view);
        immediate::set_proj_matrix(data->view.param.matrix.current.proj_jittered);

        const mat3_t box = data->simulation_box.box;
        const vec3_t min_box = box * vec3_t{0, 0, 0};
        const vec3_t max_box = box * vec3_t{1, 1, 1};

        // Simulation box
        immediate::draw_box_wireframe(min_box, max_box, convert_color(data->simulation_box.color));

        immediate::flush();

        POP_GPU_SECTION()
    }

#if 0
    // RENDER DEBUG INFORMATION (WITH DEPTH)
    PUSH_GPU_SECTION("Debug Draw") {
        glDrawBuffer(GL_COLOR_ATTACHMENT_POST_TONEMAP);
        immediate::set_model_view_matrix(data->view.param.matrix.current.view);
        immediate::set_proj_matrix(data->view.param.matrix.current.proj);
        immediate::flush();
    }
    POP_GPU_SECTION()

    PUSH_GPU_SECTION("Debug Draw Overlay") {
        glDrawBuffer(GL_COLOR_ATTACHMENT_POST_TONEMAP);  // Post_Tonemap buffer
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

    // DRAW VELOCITY OF STATIC OBJECTS
    PUSH_GPU_SECTION("Blit Static Velocity")
    glDrawBuffer(GL_COLOR_ATTACHMENT_VELOCITY);
    glDepthMask(0);
    postprocessing::blit_static_velocity(data->gbuffer.deferred.depth, data->view.param);
    glDepthMask(1);
    POP_GPU_SECTION()

    glDepthMask(1);
    glColorMask(1, 1, 1, 1);

    // DRAW REPRESENTATIONS
    PUSH_GPU_SECTION("Representation")
    glDrawBuffers(ARRAY_SIZE(draw_buffers), draw_buffers);
    draw_representations(data);
    POP_GPU_SECTION()

    PUSH_GPU_SECTION("Selection")
    const bool atom_selection_empty = md_bitfield_popcount(&data->selection.current_selection_mask) == 0;
    const bool atom_highlight_empty = md_bitfield_popcount(&data->selection.current_highlight_mask) == 0;

    glDepthMask(0);
    glColorMask(0, 0, 0, 0);

    glDrawBuffer(GL_COLOR_ATTACHMENT_POST_TONEMAP);  // Post_Tonemap buffer
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_EQUAL);

    glEnable(GL_STENCIL_TEST);
    glStencilFunc(GL_ALWAYS, 0xFF, 0xFF);

    if (!atom_selection_empty)
    {
        glStencilOp(GL_KEEP, GL_REPLACE, GL_REPLACE);
        glStencilMask(0x02);
        draw_representations_lean_and_mean(data, AtomBit_Selected | AtomBit_Visible);
    }

    if (!atom_highlight_empty)
    {
        glStencilOp(GL_KEEP, GL_REPLACE, GL_REPLACE);
        glStencilMask(0x4);
        draw_representations_lean_and_mean(data, AtomBit_Highlighted | AtomBit_Visible);
    }
        
    if (!atom_selection_empty)
    {
        glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
        glStencilMask(0x1);
        draw_representations_lean_and_mean(data, AtomBit_Selected | AtomBit_Visible);
    }
        
    glDisable(GL_DEPTH_TEST);

    glStencilMask(0x00);
    glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
    glColorMask(1, 1, 1, 1);

    if (!atom_selection_empty) {
        glStencilFunc(GL_EQUAL, 2, 2);
        postprocessing::blit_color({0, 0, 1, 0.25});

        glStencilFunc(GL_EQUAL, 2, 3);
        postprocessing::blit_color({0, 0, 0.25, 0.4});
    }

    if (!atom_highlight_empty) {
        glStencilFunc(GL_EQUAL, 4, 4);
        postprocessing::blit_color({1, 1, 0, 0.25});
    }

    glDisable(GL_STENCIL_TEST);
    if (!atom_selection_empty) {
        PUSH_GPU_SECTION("DESATURATE")
        const float saturation = data->selection.color.selection_saturation;
        glDrawBuffer(GL_COLOR_ATTACHMENT_COLOR);
        postprocessing::scale_hsv(data->gbuffer.deferred.color, vec3_t{1, saturation, 1});
        POP_GPU_SECTION()
    }

    glDisable(GL_STENCIL_TEST);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glDepthMask(1);
    POP_GPU_SECTION()

    POP_GPU_SECTION()  // G-buffer

    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
}

static void handle_picking(ApplicationData* data) {
    PUSH_CPU_SECTION("PICKING") {
        vec2_t coord = {data->ctx.input.mouse.win_coord.x, (float)data->gbuffer.height - data->ctx.input.mouse.win_coord.y};
        if (coord.x < 0.f || coord.x >= (float)data->gbuffer.width || coord.y < 0.f || coord.y >= (float)data->gbuffer.height) {
            data->picking.idx = INVALID_PICKING_IDX;
            data->picking.depth = 1.f;
        } else {
#if PICKING_JITTER_HACK
            static uint32_t frame_idx = 0;
            static uint32_t ref_frame = 0;
            frame_idx = (frame_idx + 1) % 16;
            // @NOTE: If we have jittering applied, we cannot? retreive the original pixel value (without the jitter)
            // Solution, pick one reference frame out of the jittering sequence and use that one...
            // Ugly hack but works...

            if (data->ctx.input.mouse.moving) {
                ref_frame = frame_idx;
            }

            if (ref_frame == frame_idx || data->view.param.jitter.current == vec2_t{0, 0}) {
                data->picking = read_picking_data(data->gbuffer, (int32_t)round(coord.x), (int32_t)round(coord.y));
                if (data->picking.idx != INVALID_PICKING_IDX)
                    data->picking.idx = CLAMP(data->picking.idx, 0U, (uint32_t)data->mold.mol.atom.count - 1U);
                const vec4_t viewport = {0, 0, (float)data->gbuffer.width, (float)data->gbuffer.height};
                data->picking.world_coord = mat4_unproject({coord.x, coord.y, data->picking.depth}, data->view.param.matrix.inverse.view_proj_jittered, viewport);
            }
#else
            //coord -= data->view.param.jitter.current * 2.0f;
            //coord += 0.5f;
            data->picking = read_picking_data(&data->gbuffer, (int)coord.x, (int)coord.y);
            const vec4_t viewport = {0, 0, (float)data->gbuffer.width, (float)data->gbuffer.height};
            data->picking.world_coord = mat4_unproject({coord.x, coord.y, data->picking.depth}, data->view.param.matrix.inverse.view_proj_jittered, viewport);
#endif
        }
        data->selection.hovered = -1;
        if (data->picking.idx != INVALID_PICKING_IDX) {
            data->selection.hovered = data->picking.idx;
        }
        if (data->ctx.input.mouse.clicked[1]) {
            data->selection.right_clicked = data->selection.hovered;
        }
    }
    POP_CPU_SECTION()
}
static void apply_postprocessing(const ApplicationData& data) {
    PUSH_GPU_SECTION("Postprocessing")
    postprocessing::Descriptor desc;

    desc.background.intensity = data.visuals.background.color * data.visuals.background.intensity;

    desc.ambient_occlusion.enabled = data.visuals.ssao.enabled;
    desc.ambient_occlusion.intensity = data.visuals.ssao.intensity;
    desc.ambient_occlusion.radius = data.visuals.ssao.radius;
    desc.ambient_occlusion.bias = data.visuals.ssao.bias;

    desc.tonemapping.enabled = data.visuals.tonemapping.enabled;
    desc.tonemapping.mode = data.visuals.tonemapping.tonemapper;
    desc.tonemapping.exposure = data.visuals.tonemapping.exposure;
    desc.tonemapping.gamma = data.visuals.tonemapping.gamma;

    desc.depth_of_field.enabled = data.visuals.dof.enabled;
    desc.depth_of_field.focus_depth = data.visuals.dof.focus_depth.current;
    desc.depth_of_field.focus_scale = data.visuals.dof.focus_scale;

    constexpr float MOTION_BLUR_REFERENCE_DT = 1.0f / 60.0f;
    const float dt_compensation = MOTION_BLUR_REFERENCE_DT / (float)data.ctx.timing.delta_s;
    const float motion_scale = data.visuals.temporal_reprojection.motion_blur.motion_scale * dt_compensation;
    desc.temporal_reprojection.enabled = data.visuals.temporal_reprojection.enabled;
    desc.temporal_reprojection.feedback_min = data.visuals.temporal_reprojection.feedback_min;
    desc.temporal_reprojection.feedback_max = data.visuals.temporal_reprojection.feedback_max;
    desc.temporal_reprojection.motion_blur.enabled = data.visuals.temporal_reprojection.motion_blur.enabled;
    desc.temporal_reprojection.motion_blur.motion_scale = motion_scale;

    desc.input_textures.depth = data.gbuffer.deferred.depth;
    desc.input_textures.color = data.gbuffer.deferred.color;
    desc.input_textures.normal = data.gbuffer.deferred.normal;
    desc.input_textures.velocity = data.gbuffer.deferred.velocity;
    desc.input_textures.post_tonemap = data.gbuffer.deferred.post_tonemap;

    postprocessing::shade_and_postprocess(desc, data.view.param);
    POP_GPU_SECTION()
}

static void draw_representations(ApplicationData* data) {
    ASSERT(data);
    const md_gl_representation_t* rep_data[32] = { 0 };
    uint32_t rep_count = 0;
    for (int64_t i = 0; i < md_array_size(data->representations.buffer); i++) {
        const auto& rep = data->representations.buffer[i];
        if (rep.enabled) {
            rep_data[rep_count++] = &rep.md_rep;
        }
        if (rep_count == ARRAY_SIZE(rep_data)) break;
    }

    md_gl_rendertarget_t render_target = {
        .width = (uint32_t)data->gbuffer.width,
        .height = (uint32_t)data->gbuffer.height,
        .texture_depth = data->gbuffer.deferred.depth,
        .texture_color = data->gbuffer.deferred.color,
        .texture_atom_index = data->gbuffer.deferred.picking,
        .texture_view_normal = data->gbuffer.deferred.normal,
        .texture_view_velocity = data->gbuffer.deferred.velocity,
    };

    md_gl_draw_args_t desc = {
        .representation = {
            .count = rep_count,
            .data = rep_data,
        },
        .view_transform = {
            .view_matrix = &data->view.param.matrix.current.view.elem[0][0],
            .projection_matrix = &data->view.param.matrix.current.proj_jittered.elem[0][0],
            // These two are for temporal anti-aliasing reprojection (optional)
            .prev_view_matrix = &data->view.param.matrix.previous.view.elem[0][0],
            .prev_projection_matrix = &data->view.param.matrix.previous.proj_jittered.elem[0][0],
        },
        .render_target = &render_target,
    };

    md_gl_draw(&data->mold.gl_ctx, &desc);
}

static void draw_representations_lean_and_mean(ApplicationData* data, uint32_t mask) {
    const md_gl_representation_t* rep_data[32] = { 0 };
    uint32_t rep_count = 0;
    for (int64_t i = 0; i < md_array_size(data->representations.buffer); i++) {
        const auto& rep = data->representations.buffer[i];
        if (rep.enabled) {
            rep_data[rep_count++] = &data->representations.buffer[i].md_rep;
        }
        if (rep_count == ARRAY_SIZE(rep_data)) break;
    }

    md_gl_rendertarget_t render_target = {
        .width = (uint32_t)data->gbuffer.width,
        .height = (uint32_t)data->gbuffer.height,
        .texture_depth = data->gbuffer.deferred.depth,
    };

    md_gl_draw_args_t desc = {
        .representation = {
            .count = rep_count,
            .data = rep_data,
        },
        .view_transform = {
            .view_matrix = &data->view.param.matrix.current.view.elem[0][0],
            .projection_matrix = &data->view.param.matrix.current.proj_jittered.elem[0][0],
            // These two are for temporal anti-aliasing reprojection
            //.prev_model_view_matrix = &data->view.param.matrix.previous.view[0][0],
            //.prev_projection_matrix = &data->view.param.matrix.previous.proj_jittered[0][0],
        },
        .render_target = &render_target,
        .atom_mask = mask,
    };

    md_gl_draw(&data->mold.gl_ctx, &desc);
}
