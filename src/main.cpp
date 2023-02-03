#include <core/md_compiler.h>

#if MD_COMPILER_MSVC
#   ifndef _CRT_SECURE_NO_WARNINGS
#       define _CRT_SECURE_NO_WARNINGS
#   endif
#   pragma warning( disable : 26812 4244 )
#endif

#include <md_util.h>
#include <md_gl.h>
#include <md_gfx.h>
#include <md_filter.h>
#include <md_script.h>
#include <md_molecule.h>
#include <md_trajectory.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_linear_allocator.h>
#include <core/md_tracking_allocator.h>
#include <core/md_log.h>
#include <core/md_simd.h>
#include <core/md_array.h>
#include <core/md_os.h>
#include <core/md_spatial_hash.h>
#include <core/md_base64.h>
#include <core/md_unit.h>
#include <core/md_str_builder.h>

#include "gfx/gl.h"
#include "gfx/gl_utils.h"
#include "gfx/camera.h"
#include "gfx/camera_utils.h"
#include "gfx/immediate_draw_utils.h"
#include "gfx/postprocessing_utils.h"
#include "gfx/volumerender_utils.h"

#include "string_util.h"
#include "random_util.h"
#include "imgui_widgets.h"
#include "implot_widgets.h"
#include "task_system.h"
#include "color_utils.h"
#include "loader.h"
#include "ramachandran.h"
#include "image.h"
#include "application/application.h"
#include "application/IconsFontAwesome6.h"

#include <imgui.h>
#define IMGUI_DEFINE_MATH_OPERATORS
#include <imgui_internal.h>

#include <implot.h>
#include <implot_internal.h>
#include <TextEditor.h>
#include <imgui_notify.h>

#include <stdio.h>
#include <thread>
#include <functional>

#define EXPERIMENTAL_GFX_API 0
#define PICKING_JITTER_HACK 0
#define EXPERIMENTAL_CONE_TRACED_AO 0
#define COMPILATION_TIME_DELAY_IN_SECONDS 1.0
#define IR_SEMAPHORE_MAX_COUNT 3

#define GL_COLOR_ATTACHMENT_COLOR        GL_COLOR_ATTACHMENT0
#define GL_COLOR_ATTACHMENT_NORMAL       GL_COLOR_ATTACHMENT1
#define GL_COLOR_ATTACHMENT_VELOCITY     GL_COLOR_ATTACHMENT2
#define GL_COLOR_ATTACHMENT_PICKING      GL_COLOR_ATTACHMENT3
#define GL_COLOR_ATTACHMENT_POST_TONEMAP GL_COLOR_ATTACHMENT4

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

#define LOG_INFO  MD_LOG_INFO
#define LOG_DEBUG MD_LOG_DEBUG
#define LOG_ERROR MD_LOG_ERROR
#define LOG_SUCCESS(...) ImGui::InsertNotification(ImGuiToast(ImGuiToastType_Success, 6000, __VA_ARGS__))

constexpr const char* shader_output_snippet = R"(
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
)";

constexpr const char* shader_output_snippet_lean_and_mean = R"(
void write_fragment(vec3 view_coord, vec3 view_vel, vec3 view_normal, vec4 color, uint atom_index) {
}
)";

//constexpr ImGuiKey KEY_CONSOLE = ImGuiKey_GraveAccent;
constexpr ImGuiKey KEY_PLAY_PAUSE = ImGuiKey_Space;
constexpr ImGuiKey KEY_SKIP_TO_PREV_FRAME = ImGuiKey_LeftArrow;
constexpr ImGuiKey KEY_SKIP_TO_NEXT_FRAME = ImGuiKey_RightArrow;
constexpr ImGuiKey KEY_RECOMPILE_SHADERS = ImGuiKey_F5;
constexpr ImGuiKey KEY_TOGGLE_SCREENSHOT_MODE = ImGuiKey_F10;
constexpr ImGuiKey KEY_SHOW_DEBUG_WINDOW = ImGuiKey_F11;

constexpr const char* FILE_EXTENSION = "via"; 
constexpr uint32_t INVALID_PICKING_IDX = ~0U;

//constexpr uint32_t TEXT_BG_ERROR_COLOR = 0xAA222299;
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

enum class PlaybackMode { Stopped, Playing };
enum class InterpolationMode { Nearest, Linear, CubicSpline };
enum class SelectionLevel { Atom, Residue, Chain };
enum class SelectionOperator { Or, And, AndNot, Set, Clear };
enum class SelectionGrowth { CovalentBond, Radial };
enum class RepresentationType { SpaceFill, Licorice, Ribbons, Cartoon };
enum class TrackingMode { Absolute, Relative };
enum class CameraMode { Perspective, Orthographic };

enum class ColorMapping {
    Uniform,
    Cpk,
    AtomLabel,
    AtomIndex,
    ResId,
    ResIndex,
    ChainId,
    ChainIndex,
    SecondaryStructure,
    Property
};

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
        // @NOTE: Many of each, we submit the read and use it some frame(s) later
        // This means that we read with N-1 frames latency
        GLuint color[2] = {};
        GLuint depth[2] = {};
        uint32_t frame = 0;
    } pbo_picking;

    int width = 0;
    int height = 0;
};

struct Representation {
    struct PropertyColorMapping {
        char ident[32] = ""; // property identifier
        ImPlotColormap colormap = 0;
        float range_min = 0;
        float range_max = 0;
    };

    char name[32] = "rep";
    char filt[256] = "all";
    char filt_error[256] = "";

    RepresentationType type = RepresentationType::SpaceFill;
    ColorMapping color_mapping = ColorMapping::Cpk;
    md_bitfield_t atom_mask{};
    md_gl_representation_t md_rep{};
#if EXPERIMENTAL_GFX_API
    md_gfx_handle_t gfx_rep = {};
#endif

    bool enabled = true;
    bool type_is_valid = false;
    bool filt_is_dirty = true;
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

    // Ball & Stick and Licorice, Ribbons, Cartoon
    float thickness = 1.f;

    // Ribbons, Cartoon
    float tension = 0.5f;
    float width = 1.f;

    // Property color mapping
    PropertyColorMapping prop_mapping[8] = {};

    const md_script_property_t* prop = NULL;
};

struct Selection {
    StrBuf<64> name = "sel";
    md_bitfield_t atom_mask{};
};

// This is viamd's representation of a property
struct DisplayProperty {
    enum ShowIn {
        ShowIn_Timeline = 1,
        ShowIn_Distribution = 2,
        ShowIn_Volume = 4
    };

    struct Range {
        double beg = 0;
        double end = 0;
    };

    struct Histogram {
        float bin[1024] = {};
        Range value_range = {};
    };

    struct EvaluationResult {
        str_t prop_identifier = {};
        uint64_t prop_fingerprint = 0;
        const float* val    = 0;
        const float* mean   = 0;
        const float* min    = 0;
        const float* max    = 0;
        Histogram histogram;
    };

    struct EvaluationTarget {
        char label[64] = "";
        vec4_t color = {1,1,1,1};
        Range frame_filter = {};
        EvaluationResult* results = {0};
    };

    char label[32] = "";
    uint32_t col = 0xFFFFFFFFul;

    const float* values = 0;
    md_unit_t unit = {};

    char unit_str[32] = "";

    const md_script_property_t* full_prop = NULL;
    const md_script_property_t* filt_prop = NULL;

    uint64_t full_prop_fingerprint = 0;
    uint64_t filt_prop_fingerprint = 0;

    uint32_t possible_display_mask = 0;
    uint32_t current_display_mask = 0;

    Range value_filter = {};
    Range frame_filter = {};
    Range view_range   = {};

    Histogram full_hist = {};
    Histogram filt_hist = {};
};

struct AtomElementMapping {
    StrBuf<31> lbl = "";
    md_element_t elem = 0;
};

struct LoadDatasetWindowState {
    char path_buf[2048] = "";
    bool path_is_valid = false;
    int  loader_idx = -1;
    bool load_topology = false;
    bool load_trajectory = false;
    bool coarse_grained = false;
    bool deperiodize_on_load = true;
    bool show_window = false;
};

struct ApplicationData {
    // --- APPLICATION ---
    application::Context ctx {};

    LoadDatasetWindowState load_dataset;

    // --- FILES ---
    // for keeping track of open files
    struct {
        StrBuf<1024> molecule{};
        StrBuf<1024> trajectory{};
        StrBuf<1024> workspace{};

        bool coarse_grained = false;
        bool deperiodize    = false;
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
            quat_t target_orientation = {};
            float  target_distance = 0;
        } animation;
    } view;

    struct {
        bool  hide_gui = true;
        str_t path_to_file = {};
    } screenshot;

    // --- MOLD DATA ---
    struct {
        md_allocator_i*     mol_alloc = 0;
        md_gl_shaders_t     gl_shaders = {};
        md_gl_shaders_t     gl_shaders_lean_and_mean = {};
        md_gl_molecule_t    gl_mol = {};
#if EXPERIMENTAL_GFX_API
        md_gfx_handle_t     gfx_structure = {};
#endif
        md_molecule_t       mol = {};
        md_trajectory_i*    traj = 0;

        vec3_t              mol_aabb_min;
        vec3_t              mol_aabb_max;

        struct {
            md_script_ir_t*   ir = 0;
            md_script_eval_t* full_eval = 0;
            md_script_eval_t* filt_eval = 0;
            md_script_visualization_t vis = {};

            // Semaphore to control access to IR
            md_semaphore_t ir_semaphore = {};

            bool compile_ir = false;
            bool eval_init = false;
            bool evaluate_full = false;
            bool evaluate_filt = false;
            double time_since_last_change = 0.0;
            uint64_t ir_fingerprint = 0;
        } script;
        uint32_t dirty_buffers = {0};

        md_unit_base_t unit_base = {};
    } mold;

    DisplayProperty* display_properties = 0;

    // --- ASYNC TASKS HANDLES ---
    struct {
        task_system::ID backbone_computations = task_system::INVALID_ID;
        task_system::ID prefetch_frames = task_system::INVALID_ID;
        task_system::ID evaluate_full = task_system::INVALID_ID;
        task_system::ID evaluate_filt = task_system::INVALID_ID;
        task_system::ID shape_space_evaluate = task_system::INVALID_ID;
        task_system::ID ramachandran_compute_full_density = task_system::INVALID_ID;
        task_system::ID ramachandran_compute_filt_density = task_system::INVALID_ID;
    } tasks;

    // --- ATOM SELECTION ---
    struct {
        SelectionLevel granularity = SelectionLevel::Atom;

        int32_t hovered = -1;
        int32_t right_clicked = -1;
        SingleSelectionSequence single_selection_sequence;

        md_bitfield_t current_selection_mask{};
        md_bitfield_t current_highlight_mask{};
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
            char error[256] = "";
            md_bitfield_t mask = {0};
            bool query_ok = false;
            bool query_invalid = true;
            bool show_window = false;
        } query;

        struct {
            md_bitfield_t mask = {0};
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
        float tension = 0.0f;
        InterpolationMode interpolation = InterpolationMode::CubicSpline;
        PlaybackMode mode = PlaybackMode::Stopped;
        bool apply_pbc = false;
    } animation;

    // --- TIMELINE---
    struct {
        struct {
            bool enabled = false;
            double beg_frame = 0;
            double end_frame = 1;
            
            struct {
                bool enabled = false;
                double extent_in_frames = 10;
            } temporal_window;

            uint64_t fingerprint = 0;
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
                float alpha_scale = 1.f;
                ImPlotColormap colormap = ImPlotColormap_Plasma; // Corresponds to Plasma colormap
                bool dirty = true;
            } tf;
        } dvr;

        struct {
            bool enabled = false;
            float values[8] = {};
            vec4_t colors[8] = {};
            int count = 0;
            //IsoSurfaces isosurfaces;
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

        vec3_t voxel_spacing = {1.0f, 1.0f, 1.0f};
        float resolution_scale = 2.0f;

        vec4_t clip_volume_color = {1,0,0,1};
        vec4_t bounding_box_color = {0,0,0,1};

        bool show_bounding_box = true;
        bool show_reference_structures = true;
        bool show_reference_ensemble = false;
        bool show_density_volume = false;

        bool dirty_rep = false;
        bool dirty_vol = false;

        struct {
            RepresentationType type = RepresentationType::SpaceFill;
            ColorMapping colormap = ColorMapping::Cpk;
            float radius = 1.0f;
            float tension = 0.5f;
            float width = 1.0f;
            float thickness = 1.0f;
            vec4_t color = {1,1,1,1};
        } rep;

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
            float focus_depth = 10.0f;
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
            float gamma = 2.2f;
        } tonemapping;

        struct {
            bool draw_control_points = false;
            bool draw_spline = false;
        } spline;
    } visuals;

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
    
    // --- RAMACHANDRAN ---
    struct {
        StrBuf<256> input = "all";
        StrBuf<256> error = "";

        bool evaluate = true;
        bool input_valid = false;
        bool show_window = false;

        uint64_t backbone_fingerprint = 0;
        uint64_t full_fingerprint = 0;
        uint64_t filt_fingerprint = 0;

        rama_data_t data = {0};
        uint32_t* rama_type_indices[4] = {0};

        struct {
            vec4_t border_color     = {0.0f, 0.0f, 0.0f, 1.0f};
            vec4_t base_color       = {0.9f, 0.9f, 0.9f, 0.9f};
            vec4_t selection_color  = {0.5f, 0.5f, 0.8f, 0.9f};
            vec4_t highlight_color  = {0.8f, 0.8f, 0.5f, 0.9f};
            float base_radius       = 4.5f;
            float selected_radius   = 6.5f;
        } style;

        float blur_sigma = 5.0f;
    } ramachandran;

    struct {
        char input[256] = "all";
        char error[256] = "";
        
        bool evaluate    = true;
        bool input_valid = false;
        bool show_window = false;

        // coords and weights should be of size num_frames * num_structures
        int32_t num_frames;
        int32_t num_structures;

        vec3_t* weights = 0;
        vec2_t* coords = 0;

        md_array(md_bitfield_t) bitfields = 0;

        float marker_size = 1.4f;
    } shape_space;

    // --- REPRESENTATIONS ---
    struct {
        md_array(Representation) reps = 0;
        md_bitfield_t atom_visibility_mask = {0};
        bool atom_visibility_mask_dirty = false;
        bool show_window = false;
    } representation;

    struct {
        bool show_window = false;
        AtomElementMapping* atom_element_remappings = 0;
    } dataset;

    struct {
        struct {
            int64_t stride = 0; // = mol.backbone.count. Multiply frame idx with this to get the data
            int64_t count = 0;  // = mol.backbone.count * num_frames. Defines the end of the data for assertions
            md_secondary_structure_t* data = NULL;
            uint64_t fingerprint = 0;
        } secondary_structure;
        struct {
            int64_t stride = 0; // = mol.backbone.count. Multiply frame idx with this to get the data
            int64_t count = 0;  // = mol.backbone.count * num_frames. Defines the end of the data for assertions
            md_backbone_angles_t* data = NULL;
            uint64_t fingerprint = 0;
        } backbone_angles;
    } trajectory_data;

    struct {
        vec4_t point_color      = {1,1,1,1};
        vec4_t line_color       = {0,1,1,1};
        vec4_t triangle_color   = {1,1,0,0.7f};
    } script;

    bool show_script_window = true;
    bool show_debug_window = false;
    bool show_property_export_window = false;
};

static inline uint64_t generate_fingerprint() {
    return (uint64_t)md_time_current();
}

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

static void compute_histogram_masked(float* bins, int num_bins, float min_bin_val, float max_bin_val, const float* values, int dim, const md_bitfield_t* mask) {
    ASSERT(bins);
    ASSERT(values);
    ASSERT(mask);
    ASSERT(dim > 0);
    memset(bins, 0, sizeof(float) * num_bins);

    const int num_values = md_bitfield_popcount(mask);
    if (num_values == 0) return;

    const float bin_range = max_bin_val - min_bin_val;
    const float inv_range = 1.0f / bin_range;
    const float scl = 1.0f / num_values;

    uint64_t beg_bit = mask->beg_bit;
    uint64_t end_bit = mask->end_bit;

    // We evaluate each frame, one at a time
    while ((beg_bit = md_bitfield_scan(mask, beg_bit, end_bit)) != 0) {
        uint64_t i = dim * (beg_bit - 1);
        if (values[i] < min_bin_val || max_bin_val < values[i]) continue;
        int idx = CLAMP((int)(((values[i] - min_bin_val) * inv_range) * num_bins), 0, num_bins - 1);
        bins[idx] += scl;
    }
}

static void downsample_histogram(float* dst_bins, int num_dst_bins, const float* src_bins, int num_src_bins) {
    ASSERT(num_dst_bins <= num_src_bins);

    memset(dst_bins, 0, sizeof(float) * num_dst_bins);
    const int factor = MAX(1, num_src_bins / num_dst_bins);
    const float scl = 1.0f / (float)factor;
    for (int j = 0; j < num_src_bins; ++j) {
        dst_bins[j / factor] += src_bins[j] * scl;
    }
}

static double frame_to_time(double frame, const ApplicationData& data) {
    int64_t num_frames = md_array_size(data.timeline.x_values);
    ASSERT(num_frames);
    int64_t f0 = CLAMP((int64_t)frame, 0, num_frames - 1);
    int64_t f1 = CLAMP((int64_t)frame + 1, 0, num_frames - 1);
    return lerp(data.timeline.x_values[f0], data.timeline.x_values[f1], fract(frame));
}

// Try to map time t back into frame
static double time_to_frame(double time, const ApplicationData& data) {
    int64_t num_frames = md_array_size(data.timeline.x_values);
    ASSERT(num_frames);

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
    for (size_t i = 0; i < ARRAY_SIZE(seq->idx); ++i) {
        seq->idx[i] = -1;
    }
}

static void single_selection_sequence_push_idx(SingleSelectionSequence* seq, int32_t idx) {
    ASSERT(seq);
    for (size_t i = 0; i < ARRAY_SIZE(seq->idx); ++i) {
        if (seq->idx[i] == -1) {
            seq->idx[i] = idx;
            break;
        }
    }
}

static void single_selection_sequence_pop_idx(SingleSelectionSequence* seq, int32_t idx) {
    ASSERT(seq);
    for (size_t i = 0; i < ARRAY_SIZE(seq->idx); ++i) {
        if (seq->idx[i] == idx) {
            for (size_t j = i; j < ARRAY_SIZE(seq->idx) - 1; ++j) {
                seq->idx[j] = seq->idx[j+1];
            }
            seq->idx[ARRAY_SIZE(seq->idx)-1] = -1;
            break;
        }
    }
}

static void single_selection_sequence_pop_back(SingleSelectionSequence* seq) {
    ASSERT(seq);
    size_t i = 0;
    for (; i < ARRAY_SIZE(seq->idx); ++i) {
        if (seq->idx[i] == -1) break;
    }
    if (i > 0) {
        seq->idx[i-1] = -1;
    }
}

static int32_t single_selection_sequence_last(const SingleSelectionSequence* seq) {
    ASSERT(seq);
    size_t i = 0;
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
    for (; i < (int64_t)ARRAY_SIZE(seq->idx); ++i) {
        if (seq->idx[i] == -1) break;
    }
    return i;
}

static void launch_prefetch_job(ApplicationData* data);

static void init_display_properties(DisplayProperty** prop_items, const md_script_property_t* full_props, const md_script_property_t* filt_props, int64_t num_props);
static void update_display_properties(ApplicationData* data);

static void update_density_volume(ApplicationData* data);
static void clear_density_volume(ApplicationData* data);

static void interpolate_atomic_properties(ApplicationData* data);
static void update_view_param(ApplicationData* data);
static void reset_view(ApplicationData* data, bool move_camera = false, bool smooth_transition = false);

static void compute_aabb(vec3_t* aabb_min, vec3_t* aabb_max, const float* x, const float* y, const float* z, int64_t count);

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

static void draw_load_dataset_window(ApplicationData* data);
static void draw_main_menu(ApplicationData* data);
static void draw_context_popup(ApplicationData* data);
static void draw_selection_query_window(ApplicationData* data);
static void draw_selection_grow_window(ApplicationData* data);
static void draw_animation_window(ApplicationData* data);
static void draw_representations_window(ApplicationData* data);
static void draw_timeline_window(ApplicationData* data);
static void draw_distribution_window(ApplicationData* data);
static void draw_ramachandran_window(ApplicationData* data);
static void draw_atom_info_window(const ApplicationData& data, int atom_idx);
static void draw_async_task_window(ApplicationData* data);
static void draw_shape_space_window(ApplicationData* data);
static void draw_density_volume_window(ApplicationData* data);
static void draw_script_editor_window(ApplicationData* data);
static void draw_dataset_window(ApplicationData* data);
static void draw_debug_window(ApplicationData* data);
static void draw_property_export_window(ApplicationData* data);
static void draw_notifications_window();

static void clear_gbuffer(GBuffer* gbuf);
static void init_gbuffer(GBuffer* gbuf, int width, int height);
static void destroy_gbuffer(GBuffer* gbuf);
static PickingData read_picking_data(GBuffer* fbo, int32_t x, int32_t y);

static void update_md_buffers(ApplicationData* data);

static void init_molecule_data(ApplicationData* data);
static void init_trajectory_data(ApplicationData* data);

static void interrupt_async_tasks(ApplicationData* data);

static bool load_dataset_from_file(ApplicationData* data, str_t path_to_file, md_molecule_api* mol_api = NULL, md_trajectory_api* traj_api = NULL, bool coarse_grained = false, bool deperiodize_on_load = true);

static void load_workspace(ApplicationData* data, str_t file);
static void save_workspace(ApplicationData* data, str_t file);

static void create_screenshot(ApplicationData* data);

// Representations
static Representation* create_representation(ApplicationData* data, RepresentationType type = RepresentationType::SpaceFill,
                                             ColorMapping color_mapping = ColorMapping::Cpk, str_t filter = STR("all"));
static Representation* clone_representation(ApplicationData* data, const Representation& rep);
static void remove_representation(ApplicationData* data, int idx);
static void update_representation(ApplicationData* data, Representation* rep);
static void update_all_representations(ApplicationData* data);
static void init_representation(ApplicationData* data, Representation* rep);
static void init_all_representations(ApplicationData* data);
static void clear_representations(ApplicationData* data);

static void recompute_atom_visibility_mask(ApplicationData* data);

// Selections
static Selection* create_selection(ApplicationData* data, str_t name, md_bitfield_t* bf);
static Selection* clone_selection(ApplicationData* data, const Selection& sel);
static void remove_selection(ApplicationData* data, int idx);

static bool filter_expression(ApplicationData* data, str_t expr, md_bitfield_t* mask, bool* is_dynamic, char* error_str, int error_cap);

static void modify_field(md_bitfield_t* bf, const md_bitfield_t* mask, SelectionOperator op) {
    switch(op) {
    case SelectionOperator::Or:
        md_bitfield_or_inplace(bf, mask);
        break;
    case SelectionOperator::And:
        md_bitfield_and_inplace(bf, mask);
        break;
    case SelectionOperator::AndNot:
        md_bitfield_andnot_inplace(bf, mask);
        break;
    case SelectionOperator::Set:
        md_bitfield_copy(bf, mask);
        break;
    default:
        ASSERT(false);
    }
}

static void modify_field(md_bitfield_t* bf, md_range_t range, SelectionOperator op) {
    switch(op) {
    case SelectionOperator::Or:
        md_bitfield_set_range(bf, range.beg, range.end);
        break;
    case SelectionOperator::And:
        md_bitfield_clear_range(bf, 0, range.beg);
        md_bitfield_clear_range(bf, range.end, bf->end_bit);
        break;
    case SelectionOperator::AndNot:
        md_bitfield_clear_range(bf, range.beg, range.end);
        break;
    case SelectionOperator::Set:
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

static void modify_selection(ApplicationData* data, md_bitfield_t* atom_mask, SelectionOperator op = SelectionOperator::Set) {
    ASSERT(data);
    modify_field(&data->selection.current_selection_mask, atom_mask, op);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

static void modify_selection(ApplicationData* data, md_range_t range, SelectionOperator op = SelectionOperator::Set) {
    ASSERT(data);
    modify_field(&data->selection.current_selection_mask, range, op);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

// Global data for application
static md_allocator_i* frame_allocator = 0;
static TextEditor editor {};    // We do not want this within the application data since it messes up the layout and therefore the usage of offset_of
static bool use_gfx = false;
#if DEBUG
static md_allocator_i* persistent_allocator = md_tracking_allocator_create(default_allocator);
#elif RELEASE
static md_allocator_i* persistent_allocator = default_allocator;
#else
    
#endif

int main(int, char**) {
    const int64_t linear_size = MEGABYTES(256);
    void* linear_mem = md_alloc(default_allocator, linear_size);
    md_linear_allocator_t linear_alloc {};
    md_linear_allocator_init(&linear_alloc, linear_mem, linear_size);
    md_allocator_i linear_interface = md_linear_allocator_create_interface(&linear_alloc);
    frame_allocator = &linear_interface;

    md_logger_i notification_logger = {
        NULL,
        [](struct md_logger_o* inst, enum md_log_type_t log_type, const char* msg) {
            (void)inst;
    
            static char prev_msg[1024];
            if (strncmp(prev_msg, msg, sizeof(prev_msg)) == 0) {
                return;
            }
            strncpy(prev_msg, msg, sizeof(prev_msg));
            
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
                ImGui::InsertNotification(ImGuiToast(toast_type, 6000, msg));
            }
        }
    };

    md_logger_add(&notification_logger);

    ApplicationData data;

    data.mold.mol_alloc = md_arena_allocator_create(persistent_allocator, MEGABYTES(1));

    md_bitfield_init(&data.selection.current_selection_mask, persistent_allocator);
    md_bitfield_init(&data.selection.current_highlight_mask, persistent_allocator);
    md_bitfield_init(&data.selection.query.mask, persistent_allocator);
    md_bitfield_init(&data.selection.grow.mask, persistent_allocator);

    md_bitfield_init(&data.representation.atom_visibility_mask, persistent_allocator);

    md_semaphore_init(&data.mold.script.ir_semaphore, IR_SEMAPHORE_MAX_COUNT);

    data.mold.script.ir = md_script_ir_create(persistent_allocator);

    // Init platform
    LOG_DEBUG("Initializing GL...");
    if (!application::initialize(&data.ctx, 0, 0, "VIAMD")) {
        LOG_ERROR("Could not initialize platform layer... terminating\n");
        return -1;
    }
    data.ctx.window.vsync = true;

    LOG_DEBUG("Creating framebuffer...");
    init_gbuffer(&data.gbuffer, data.ctx.framebuffer.width, data.ctx.framebuffer.height);

    for (int i = 0; i < (int)ARRAY_SIZE(data.view.jitter.sequence); ++i) {
        data.view.jitter.sequence[i].x = halton(i + 1, 2);
        data.view.jitter.sequence[i].y = halton(i + 1, 3);
    }

    // Init subsystems
    LOG_DEBUG("Initializing immediate draw...");
    immediate::initialize();
    LOG_DEBUG("Initializing ramachandran...");
    ramachandran::initialize();
    LOG_DEBUG("Initializing post processing...");
    postprocessing::initialize(data.gbuffer.width, data.gbuffer.height);
    LOG_DEBUG("Initializing volume...");
    volume::initialize();
#if EXPERIMENTAL_CONE_TRACED_AO == 1
    LOG_DEBUG("Initializing cone tracing...");
    cone_trace::initialize();
#endif
    LOG_DEBUG("Initializing task system...");
    task_system::initialize(MIN(VIAMD_NUM_WORKER_THREADS, std::thread::hardware_concurrency()));

    md_gl_initialize();
    md_gl_shaders_init(&data.mold.gl_shaders, shader_output_snippet);
    md_gl_shaders_init(&data.mold.gl_shaders_lean_and_mean, shader_output_snippet_lean_and_mean);

#if EXPERIMENTAL_GFX_API
    md_gfx_initialize(data.gbuffer.width, data.gbuffer.height, 0);
#endif

    ImGui::init_theme();

    editor.SetLanguageDefinition(TextEditor::LanguageDefinition::VIAMD());
    editor.SetPalette(TextEditor::GetDarkPalette());

    str_t path = STR(VIAMD_DATASET_DIR "/1ALA-500.pdb");

    load_dataset_from_file(&data, path, load::mol::get_api(path), load::traj::get_api(path), false, true);
    create_representation(&data, RepresentationType::SpaceFill, ColorMapping::Cpk, STR("all"));
    editor.SetText("s1 = resname(\"ALA\")[2:8];\nd1 = distance(10,30);\na1 = angle(1,2,3) in resname(\"ALA\");\nr = rdf(element('C'), element('H'), 10.0);\nv = sdf(s1, element('H'), 10.0);");

    reset_view(&data, true);
    recompute_atom_visibility_mask(&data);
    interpolate_atomic_properties(&data);

    rama_init(&data.ramachandran.data);

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

    //bool demo_window = false;

    // Main loop
    while (!data.ctx.window.should_close) {
        application::update(&data.ctx);
        
        // This needs to happen first (in imgui events) to enable docking of imgui windows
#ifdef IMGUI_DOCKING
        ImGui::CreateDockspace();
#endif

        const int64_t num_frames = md_trajectory_num_frames(data.mold.traj);
        const int64_t last_frame = MAX(0, num_frames - 1);
        const double   max_frame = (double)last_frame;

        md_bitfield_clear(&data.selection.current_highlight_mask);

        // GUI
        if (data.load_dataset.show_window) draw_load_dataset_window(&data);
        if (data.representation.show_window) draw_representations_window(&data);
        if (data.density_volume.show_window) draw_density_volume_window(&data);
        if (data.timeline.show_window) draw_timeline_window(&data);
        if (data.distributions.show_window) draw_distribution_window(&data);
        if (data.ramachandran.show_window) draw_ramachandran_window(&data);
        if (data.shape_space.show_window) draw_shape_space_window(&data);
        if (data.show_script_window) draw_script_editor_window(&data);
        if (data.dataset.show_window) draw_dataset_window(&data);
        if (data.selection.query.show_window) draw_selection_query_window(&data);
        if (data.selection.grow.show_window) draw_selection_grow_window(&data);
        if (data.show_property_export_window) draw_property_export_window(&data);
        if (data.show_debug_window) draw_debug_window(&data);

        data.selection.selecting = false;

        //if (demo_window) ImGui::ShowDemoWindow(&demo_window);

        handle_camera_interaction(&data);
        handle_camera_animation(&data);
        update_view_param(&data);

        draw_context_popup(&data);
        draw_async_task_window(&data);
        draw_animation_window(&data);
        draw_main_menu(&data);
        draw_notifications_window();

        // Capture non-window specific keyboard events
        if (!ImGui::GetIO().WantCaptureKeyboard) {
#if EXPERIMENTAL_GFX_API
            if (ImGui::IsKeyPressed(ImGuiKey_F1)) {
                use_gfx = !use_gfx;
            }
#endif

            if (ImGui::IsKeyPressed(KEY_SHOW_DEBUG_WINDOW)) {
                data.show_debug_window = true;
            }

            if (ImGui::IsKeyPressed(KEY_TOGGLE_SCREENSHOT_MODE)) {
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

            if (ImGui::IsKeyPressed(KEY_RECOMPILE_SHADERS)) {
                LOG_INFO("Recompiling shaders and re-initializing volume");
                postprocessing::initialize(data.gbuffer.width, data.gbuffer.height);
                ramachandran::initialize();
                volume::initialize();
#if EXPERIMENTAL_CONE_TRACED_AO == 1
                cone_trace::initialize();
#endif
                md_gl_shaders_free(&data.mold.gl_shaders);
                md_gl_shaders_init(&data.mold.gl_shaders, shader_output_snippet);
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
            data.animation.frame += data.ctx.timing.delta_s * data.animation.fps;
            data.animation.frame = CLAMP(data.animation.frame, 0.0, max_frame);
            if (data.animation.frame >= max_frame) {
                data.animation.mode = PlaybackMode::Stopped;
                data.animation.frame = max_frame;
            }

            if (!task_system::task_is_running(data.tasks.prefetch_frames)) {
                uint32_t traj_frames = md_trajectory_num_frames(data.mold.traj);
                if (traj_frames > 0 && load::traj::num_cache_frames(data.mold.traj) < traj_frames) {
                    uint32_t frame_beg = 0;
                    uint32_t frame_end = 0;
                    // @NOTE: This is certainly something which can be improved upon.
                    // It prefetches frames in the direction of the animation.
                    // In a more optimal case, it should never yield until it catches up with the number of frames it expects to have as a buffer.
                    
                    if (data.animation.fps > 0) {
                        frame_beg = CLAMP((uint32_t)data.animation.frame, 0, traj_frames - 1);
                        frame_end = CLAMP((uint32_t)data.animation.frame + 20, 0, traj_frames);
                    } else {
                        frame_beg = CLAMP((uint32_t)data.animation.frame - 20, 0, traj_frames - 1);
                        frame_end = CLAMP((uint32_t)data.animation.frame, 0, traj_frames);
                    }
                    if (frame_beg != frame_end) {
                        data.tasks.prefetch_frames = task_system::pool_enqueue(STR("##Prefetch Frames"), frame_beg, frame_end, [](uint32_t frame_beg, uint32_t frame_end, void* user_data) {
                            ApplicationData* data = (ApplicationData*)user_data;
                            for (uint32_t i = frame_beg; i < frame_end; ++i) {
                                md_trajectory_load_frame(data->mold.traj, i, 0, 0, 0, 0);
                            }
                        }, &data);
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
            if (data.mold.script.ir && (data.timeline.filter.beg_frame != pre_beg || data.timeline.filter.end_frame != pre_end)) {
                data.mold.script.evaluate_filt = true;
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
            if (traj) {
                interpolate_atomic_properties(&data);

                if (data.animation.apply_pbc) {
                    const vec3_t pbc_ext = data.simulation_box.box * (vec3_t {1,1,1});
                    if (mol.covalent.count > 0) {
                        //md_util_apply_pbc_preserve_covalent(mol.atom.x, mol.atom.y, mol.atom.z, &mol.covalent, pbc_ext);
                        ASSERT(false);
                    } else {
                        md_util_pbc_ortho(mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, pbc_ext);
                    }
                }

#if EXPERIMENTAL_CONE_TRACED_AO
                if (data.visuals.cone_traced_ao.enabled) {
                    init_occupancy_volume(&data);
                }
#endif
            }
            POP_CPU_SECTION()

            PUSH_CPU_SECTION("Update dynamic representations")
            for (int64_t i = 0; i < md_array_size(data.representation.reps); ++i) {
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

        if (data.mold.script.compile_ir) {
            data.mold.script.time_since_last_change += data.ctx.timing.delta_s;

            if (data.mold.script.time_since_last_change > COMPILATION_TIME_DELAY_IN_SECONDS) {
                // We cannot recompile while it is evaluating.
                // Need to interrupt and wait for tasks to finish.
                if (data.mold.script.full_eval) md_script_eval_interrupt(data.mold.script.full_eval);
                if (data.mold.script.filt_eval) md_script_eval_interrupt(data.mold.script.filt_eval);

                // Try aquire all semaphores
                if (md_semaphore_try_aquire_n(&data.mold.script.ir_semaphore, IR_SEMAPHORE_MAX_COUNT)) {
                    defer { md_semaphore_release_n(&data.mold.script.ir_semaphore, IR_SEMAPHORE_MAX_COUNT); };

                    // Now we hold all semaphores for the script
                    data.mold.script.compile_ir = false;
                    data.mold.script.time_since_last_change = 0;
                    
                    md_script_ir_clear(data.mold.script.ir);

                    std::string src = editor.GetText();
                    str_t src_str {src.data(), (int64_t)src.length()};

                    editor.ClearMarkers();
                    editor.ClearErrorMarkers();
                    
                    if (src_str) {
                        const int64_t num_stored_selections = md_array_size(data.selection.stored_selections);                       
                        if (num_stored_selections > 0) {
                            md_script_bitfield_identifier_t* idents = 0;
                            for (int64_t i = 0; i < num_stored_selections; ++i) {
                                md_script_bitfield_identifier_t ident = {
                                    .identifier_name = data.selection.stored_selections[i].name,
                                    .bitfield = &data.selection.stored_selections[i].atom_mask,
                                };
                                md_array_push(idents, ident, frame_allocator);
                            }
                            md_script_ir_add_bitfield_identifiers(data.mold.script.ir, idents, md_array_size(idents));
                        }
                        md_script_ir_compile_from_source(data.mold.script.ir, src_str, &data.mold.mol, NULL);

                        if (md_script_ir_valid(data.mold.script.ir)) {
                            uint64_t ir_figerprint = md_script_ir_fingerprint(data.mold.script.ir);
                            if (data.mold.script.ir_fingerprint != ir_figerprint) {
                                data.mold.script.ir_fingerprint = ir_figerprint;
                            }
                        }

                        TextEditor::ErrorMarkers markers{};
                        const int64_t num_errors = md_script_ir_num_errors(data.mold.script.ir);
                        if (num_errors) {
                            const md_script_error_t* errors = md_script_ir_errors(data.mold.script.ir);
                            for (int64_t i = 0; i < num_errors; ++i) {
                                std::string err_str(errors[i].text.ptr, errors[i].text.len);
                                auto first = editor.GetCharacterCoordinates(errors[i].range.beg);
                                auto last  = editor.GetCharacterCoordinates(errors[i].range.end-1);
                                for (int j = first.mLine; j <= last.mLine; ++j) {
                                    std::pair<int, std::string> pair = {j + 1, err_str};
                                    markers.insert(pair);
                                }
                            }
                        }
                        editor.SetErrorMarkers(markers);

                        const int64_t num_tokens = md_script_ir_num_vis_tokens(data.mold.script.ir);
                        const md_script_vis_token_t* vis_tokens = md_script_ir_vis_tokens(data.mold.script.ir);
                        for (int64_t i = 0; i < num_tokens; ++i) {
                            const md_script_vis_token_t& vis_tok = vis_tokens[i];
                            TextEditor::Marker marker = {0};
                            auto first = editor.GetCharacterCoordinates(vis_tok.range.beg);
                            auto last  = editor.GetCharacterCoordinates(vis_tok.range.end);
                            marker.begCol = first.mColumn;
                            marker.endCol = last.mColumn;
                            marker.bgColor = ImVec4(1,1,1,0.5);
                            marker.depth = vis_tok.depth;
                            marker.onlyShowBgOnMouseOver = true;
                            marker.text = std::string(vis_tok.text.ptr, vis_tok.text.len);
                            marker.payload = (void*)vis_tok.payload;
                            marker.line = first.mLine + 1;
                            editor.AddMarker(marker);
                        }
                    }
                }
            }
        }

        if (num_frames > 0) {
            if (data.mold.script.eval_init) {
                data.mold.script.eval_init = false;

                if (data.mold.script.full_eval) {
                    md_script_eval_free(data.mold.script.full_eval);
                }
                if (data.mold.script.filt_eval) {
                    md_script_eval_free(data.mold.script.filt_eval);
                }
                
                if (md_script_ir_valid(data.mold.script.ir)) {
                    data.mold.script.full_eval = md_script_eval_create(num_frames, data.mold.script.ir, persistent_allocator);
                    data.mold.script.filt_eval = md_script_eval_create(num_frames, data.mold.script.ir, persistent_allocator);
                }

                const int64_t num_props = md_script_eval_num_properties(data.mold.script.full_eval);
                ASSERT(md_script_eval_num_properties(data.mold.script.filt_eval) == num_props);
                init_display_properties(&data.display_properties, md_script_eval_properties(data.mold.script.full_eval), md_script_eval_properties(data.mold.script.filt_eval), num_props);

                data.mold.script.evaluate_filt = true;
                data.mold.script.evaluate_full = true;
            }

            if (data.mold.script.full_eval && data.mold.script.evaluate_full) {
                if (task_system::task_is_running(data.tasks.evaluate_full)) {
                    md_script_eval_interrupt(data.mold.script.full_eval);
                } else {
                    if (md_semaphore_try_aquire(&data.mold.script.ir_semaphore)) {
                        if (md_script_ir_valid(data.mold.script.ir) &&
                            md_script_eval_ir_fingerprint(data.mold.script.full_eval) == md_script_ir_fingerprint(data.mold.script.ir))
                        {
                            data.mold.script.evaluate_full = false;
                            md_script_eval_clear(data.mold.script.full_eval);

                            data.tasks.evaluate_full = task_system::pool_enqueue(STR("Eval Full"), 0, (uint32_t)num_frames, [](uint32_t frame_beg, uint32_t frame_end, void* user_data) {
                                ApplicationData* data = (ApplicationData*)user_data;
                                md_script_eval_frame_range(data->mold.script.full_eval, data->mold.script.ir, &data->mold.mol, data->mold.traj, frame_beg, frame_end);
                            }, &data);
                            
                            task_system::pool_enqueue(STR("##Release IR Semaphore"), [](void* user_data) {
                                ApplicationData* data = (ApplicationData*)user_data;
                                md_semaphore_release(&data->mold.script.ir_semaphore);
                            }, &data, data.tasks.evaluate_full);
                        }
                        else {
                            md_semaphore_release(&data.mold.script.ir_semaphore);
                        }
                    }
                }
            }

            if (data.mold.script.filt_eval && data.mold.script.evaluate_filt && data.timeline.filter.enabled) {
                if (task_system::task_is_running(data.tasks.evaluate_filt)) {
                    md_script_eval_interrupt(data.mold.script.filt_eval);
                } else {
                    if (md_semaphore_try_aquire(&data.mold.script.ir_semaphore)) {
                        if (md_script_ir_valid(data.mold.script.ir) &&
                            md_script_eval_ir_fingerprint(data.mold.script.filt_eval) == md_script_ir_fingerprint(data.mold.script.ir))
                        {
                            data.mold.script.evaluate_filt = false;
                            md_script_eval_clear(data.mold.script.filt_eval);
                            
                            const uint32_t traj_frames = (uint32_t)md_trajectory_num_frames(data.mold.traj);
                            const uint32_t beg_frame = CLAMP((uint32_t)data.timeline.filter.beg_frame, 0, traj_frames-1);
                            const uint32_t end_frame = CLAMP((uint32_t)data.timeline.filter.end_frame + 1, beg_frame + 1, traj_frames);
                            data.tasks.evaluate_filt = task_system::pool_enqueue(STR("Eval Filt"), beg_frame, end_frame, [](uint32_t beg, uint32_t end, void* user_data)
                                {
                                    ApplicationData* data = (ApplicationData*)user_data;
                                    md_script_eval_frame_range(data->mold.script.filt_eval, data->mold.script.ir, &data->mold.mol, data->mold.traj, beg, end);
                                }, &data);
                            
                            task_system::pool_enqueue(STR("##Release IR Semaphore"), [](void* user_data)
                                {
                                    ApplicationData* data = (ApplicationData*)user_data;
                                    md_semaphore_release(&data->mold.script.ir_semaphore);
                                }, &data, data.tasks.evaluate_filt);
                        } else {
                            md_semaphore_release(&data.mold.script.ir_semaphore);
                        }
                    }
                }
            }
        }


#if 0
        PUSH_CPU_SECTION("Hydrogen bonds")
        if (data.hydrogen_bonds.enabled && data.hydrogen_bonds.dirty) {
            data.hydrogen_bonds.bonds = hydrogen::compute_bonds(
                {mol.hydrogen.donor.data, mol.hydrogen.donor.count}, {mol.hydrogen.acceptor.data, mol.hydrogen.acceptor.count},
                mol.atom.position, data.hydrogen_bonds.distance_cutoff, DEG_TO_RAD(data.hydrogen_bonds.angle_cutoff));
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

        update_md_buffers(&data);
        update_display_properties(&data);

        if (data.mold.mol.backbone.count > 0) {
            if (data.ramachandran.backbone_fingerprint != data.trajectory_data.backbone_angles.fingerprint) {
                data.ramachandran.backbone_fingerprint = data.trajectory_data.backbone_angles.fingerprint;

                md_array_shrink(data.ramachandran.rama_type_indices[0], 0);
                md_array_shrink(data.ramachandran.rama_type_indices[1], 0);
                md_array_shrink(data.ramachandran.rama_type_indices[2], 0);
                md_array_shrink(data.ramachandran.rama_type_indices[3], 0);

                for (uint32_t i = 0; i < (uint32_t)md_array_size(data.mold.mol.backbone.ramachandran_type); ++i) {
                    switch (data.mold.mol.backbone.ramachandran_type[i]) {
                    case MD_RAMACHANDRAN_TYPE_GENERAL: md_array_push(data.ramachandran.rama_type_indices[0], i, persistent_allocator); break;
                    case MD_RAMACHANDRAN_TYPE_GLYCINE: md_array_push(data.ramachandran.rama_type_indices[1], i, persistent_allocator); break;
                    case MD_RAMACHANDRAN_TYPE_PROLINE: md_array_push(data.ramachandran.rama_type_indices[2], i, persistent_allocator); break;
                    case MD_RAMACHANDRAN_TYPE_PREPROL: md_array_push(data.ramachandran.rama_type_indices[3], i, persistent_allocator); break;
                    default: break;
                    }
                }
            }

            if (data.ramachandran.full_fingerprint != data.trajectory_data.backbone_angles.fingerprint) {
                if (!task_system::task_is_running(data.tasks.ramachandran_compute_full_density)) {
                    data.ramachandran.full_fingerprint = data.trajectory_data.backbone_angles.fingerprint;
                    const uint32_t* indices[4] = {
                        data.ramachandran.rama_type_indices[0],
                        data.ramachandran.rama_type_indices[1],
                        data.ramachandran.rama_type_indices[2],
                        data.ramachandran.rama_type_indices[3],
                    };

                    const uint32_t frame_beg = 0;
                    const uint32_t frame_end = (uint32_t)num_frames;
                    const uint32_t frame_stride = (uint32_t)data.trajectory_data.backbone_angles.stride;

                    data.tasks.ramachandran_compute_full_density = rama_rep_compute_density(&data.ramachandran.data.full, data.trajectory_data.backbone_angles.data, indices, frame_beg, frame_end, frame_stride, data.ramachandran.blur_sigma);
                } else {
                    task_system::task_interrupt(data.tasks.ramachandran_compute_full_density);
                }
            }

            if (data.ramachandran.filt_fingerprint != data.timeline.filter.fingerprint) {
                if (!task_system::task_is_running(data.tasks.ramachandran_compute_filt_density)) {
                    data.ramachandran.filt_fingerprint = data.timeline.filter.fingerprint;

                    const uint32_t* indices[4] = {
                        data.ramachandran.rama_type_indices[0],
                        data.ramachandran.rama_type_indices[1],
                        data.ramachandran.rama_type_indices[2],
                        data.ramachandran.rama_type_indices[3],
                    };

                    const uint32_t frame_beg = (uint32_t)data.timeline.filter.beg_frame;
                    const uint32_t frame_end = (uint32_t)data.timeline.filter.end_frame;
                    const uint32_t frame_stride = (uint32_t)data.trajectory_data.backbone_angles.stride;

                    data.tasks.ramachandran_compute_filt_density = rama_rep_compute_density(&data.ramachandran.data.filt, data.trajectory_data.backbone_angles.data, indices, frame_beg, frame_end, frame_stride);
                }
                else {
                    task_system::task_interrupt(data.tasks.ramachandran_compute_filt_density);
                }
            }
        }

        handle_picking(&data);
        clear_gbuffer(&data.gbuffer);
        fill_gbuffer(&data);
        immediate::render();

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

        // Activate backbuffer
        glDisable(GL_DEPTH_TEST);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glViewport(0, 0, data.ctx.framebuffer.width, data.ctx.framebuffer.height);
        glDrawBuffer(GL_BACK);
        glClear(GL_COLOR_BUFFER_BIT);

        apply_postprocessing(data);

        // Render Screenshot of backbuffer without GUI here
        if (data.screenshot.hide_gui && !str_empty(data.screenshot.path_to_file)) {
            create_screenshot(&data);
            data.screenshot.path_to_file = {};
        }

        PUSH_GPU_SECTION("Imgui render")
        application::render_imgui(&data.ctx);
        POP_GPU_SECTION()

        // Render Screenshot of backbuffer with GUI here
        if (!data.screenshot.hide_gui && !str_empty(data.screenshot.path_to_file)) {
            create_screenshot(&data);
            data.screenshot.path_to_file = {};
        }

        // Swap buffers
        application::swap_buffers(&data.ctx);

        task_system::execute_queued_tasks();

        // Reset frame allocator
        md_linear_allocator_reset(&linear_alloc);
    }

    interrupt_async_tasks(&data);

    // shutdown subsystems
    LOG_DEBUG("Shutting down immediate draw...");
    immediate::shutdown();
    LOG_DEBUG("Shutting down ramachandran...");
    ramachandran::shutdown();
    LOG_DEBUG("Shutting down post processing...");
    postprocessing::shutdown();
    LOG_DEBUG("Shutting down volume...");
    volume::shutdown();
#if EXPERIMENTAL_CONE_TRACED_AO == 1
    LOG_INFO("Shutting down cone tracing...");
    cone_trace::shutdown();
#endif
    LOG_DEBUG("Shutting down task system...");
    task_system::shutdown();

    destroy_gbuffer(&data.gbuffer);
    application::shutdown(&data.ctx);

    return 0;
}

static void init_display_properties(DisplayProperty** prop_items, const md_script_property_t* full_props, const md_script_property_t* filt_props, int64_t num_props) {
    ASSERT(prop_items);
    if (num_props > 0) {
        ASSERT(full_props);
        ASSERT(filt_props);
    }

    DisplayProperty* new_items = 0;
    DisplayProperty* old_items = *prop_items;

    for (int64_t i = 0; i < num_props; ++i) {
        uint32_t possible_display_mask = 0;
        uint32_t flags = full_props[i].flags;

        if (flags & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) {
            possible_display_mask = DisplayProperty::ShowIn_Timeline | DisplayProperty::ShowIn_Distribution;
        }
        else if (flags & MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) {
            possible_display_mask = DisplayProperty::ShowIn_Distribution;
        }
        else if (flags & MD_SCRIPT_PROPERTY_FLAG_VOLUME) {
            possible_display_mask = DisplayProperty::ShowIn_Volume;
        }
        else {            
            ASSERT(false);
        }

        DisplayProperty item;
        strncpy(item.label, full_props[i].ident.ptr, MIN(sizeof(item.label), (size_t)full_props[i].ident.len));
        item.col = PROPERTY_COLORS[i % ARRAY_SIZE(PROPERTY_COLORS)];
        item.frame_filter = {0, 0};
        item.value_filter = {full_props[i].data.min_range[0], full_props[i].data.max_range[0]};
        item.view_range   = {full_props[i].data.min_range[0], full_props[i].data.max_range[0]};
        item.unit = full_props[i].data.unit;
        item.full_prop = &full_props[i];
        item.filt_prop = &filt_props[i];
        item.full_prop_fingerprint = 0;
        item.filt_prop_fingerprint = 0;
        item.possible_display_mask = possible_display_mask;
        item.current_display_mask  = possible_display_mask;

        unit_print(item.unit_str, sizeof(item.unit_str), item.unit);

        for (int64_t j = 0; j < md_array_size(old_items); ++j) {
            if (strcmp(item.label, old_items[j].label) == 0) {
                // Copy relevant parameters from existing item which we want to be persistent
                item.frame_filter = old_items[j].frame_filter;
                item.value_filter = old_items[j].value_filter;
                item.view_range   = old_items[j].view_range;
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
//    const int32_t num_frames = (int32_t)md_trajectory_num_frames(data->mold.traj);

    DisplayProperty* disp_props = data->display_properties;
    for (int64_t i = 0; i < md_array_size(disp_props); ++i) {
        if (disp_props[i].full_prop_fingerprint != disp_props[i].full_prop->data.fingerprint) {
            disp_props[i].full_prop_fingerprint = disp_props[i].full_prop->data.fingerprint;
            const md_script_property_t* p = disp_props[i].full_prop;
            if (p->flags & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) {
                DisplayProperty::Histogram& hist = disp_props[i].full_hist;
                hist.value_range = {p->data.min_range[0], p->data.max_range[0]};
                const md_bitfield_t* completed_frames = md_script_eval_completed_frames(data->mold.script.full_eval);
                //const int num_completed = md_bitfield_popcount(completed_frames);
                //const int num_values = p->data.dim[0] * num_completed;
                compute_histogram_masked(hist.bin, (int)ARRAY_SIZE(hist.bin), hist.value_range.beg, hist.value_range.end, p->data.values, p->data.dim[0], completed_frames);
            }
        }

        if (disp_props[i].filt_prop_fingerprint != disp_props[i].filt_prop->data.fingerprint) {
            disp_props[i].filt_prop_fingerprint = disp_props[i].filt_prop->data.fingerprint;
            const md_script_property_t* p = disp_props[i].full_prop;
            if (p->flags & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) {
                DisplayProperty::Histogram& hist = disp_props[i].filt_hist;
                // We copy the full range here since we want the histograms to align on the x-axis
                // This also implys that we have a dependency to the full property, since we can only now the full range of values when full_prop has been completely evaluated.
                hist.value_range = disp_props[i].full_hist.value_range;
                //int beg_frame = CLAMP((int)data->timeline.filter.beg_frame, 0, num_frames-1);
                //int end_frame = CLAMP((int)data->timeline.filter.end_frame+1, 0, num_frames);
                //end_frame = MAX(beg_frame + 1, end_frame);

                //int offset = beg_frame * p->data.dim[0];
                //int length = (end_frame - beg_frame) * p->data.dim[0];
                //ASSERT(offset + length <= p->data.num_values);
                const md_bitfield_t* completed_frames = md_script_eval_completed_frames(data->mold.script.filt_eval);

                compute_histogram_masked(hist.bin, (int)ARRAY_SIZE(hist.bin), hist.value_range.beg, hist.value_range.end, p->data.values, p->data.dim[0], completed_frames);
            }
        }
    }
}

static void update_density_volume(ApplicationData* data) {
    if (data->density_volume.dvr.tf.dirty) {
        data->density_volume.dvr.tf.dirty = false;

        uint32_t pixel_data[128];

        for (size_t i = 0; i < ARRAY_SIZE(pixel_data); ++i) {
            float t = (float)i / (float)(ARRAY_SIZE(pixel_data) - 1);
            ImVec4 col = ImPlot::SampleColormap(t, data->density_volume.dvr.tf.colormap);

            // This is a small alpha ramp in the start of the TF to avoid rendering low density values.
            col.w = MIN(160 * t*t, 0.341176);
            pixel_data[i] = ImGui::ColorConvertFloat4ToU32(col);
        }

        gl::init_texture_1D(&data->density_volume.dvr.tf.id, (int)ARRAY_SIZE(pixel_data), GL_RGBA8);
        gl::set_texture_1D_data(data->density_volume.dvr.tf.id, pixel_data, GL_RGBA8);
    }

    int64_t selected_property = -1;
    for (int64_t i = 0; i < md_array_size(data->display_properties); ++i) {
        if (data->display_properties[i].current_display_mask & DisplayProperty::ShowIn_Volume) {
            selected_property = i;
            break;
        }
    }

    const md_script_property_t* prop = 0;
    uint64_t data_fingerprint = 0;

    static int64_t s_selected_property = 0;
    if (s_selected_property != selected_property) {
        s_selected_property = selected_property;
        data->density_volume.dirty_vol = true;
        data->density_volume.dirty_rep = true;
    }

    if (selected_property != -1) {
        if (data->timeline.filter.enabled) {
            prop = data->display_properties[selected_property].filt_prop;
            data_fingerprint = data->display_properties[selected_property].filt_prop->data.fingerprint;
        }
        else {
            prop = data->display_properties[selected_property].full_prop;
            data_fingerprint = data->display_properties[selected_property].full_prop->data.fingerprint;
        }
    }
    data->density_volume.show_density_volume = selected_property != -1;

    static uint64_t s_script_fingerprint = 0;
    if (s_script_fingerprint != md_script_ir_fingerprint(data->mold.script.ir)) {
        s_script_fingerprint = md_script_ir_fingerprint(data->mold.script.ir);
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
        if (prop) {
            data->density_volume.dirty_rep = false;
            int64_t num_reps = 0;
            bool result = false;
            md_script_visualization_t vis = {};

            if (md_semaphore_aquire(&data->mold.script.ir_semaphore)) {
                defer { md_semaphore_release(&data->mold.script.ir_semaphore); };
                if (md_script_ir_valid(data->mold.script.ir)) {
                    md_script_visualization_args_t args = {
                        .payload = prop->vis_payload,
                        .ir = data->mold.script.ir,
                        .mol = &data->mold.mol,
                        .traj = data->mold.traj,
                        .alloc = frame_allocator,
                        .flags = MD_SCRIPT_VISUALIZE_SDF
                    };
                    result = md_script_visualization_init(&vis, args);
                }
            }

            if (result) {
                if (vis.sdf.extent) {
                    const float s = vis.sdf.extent;
                    vec3_t min_aabb = { -s, -s, -s };
                    vec3_t max_aabb = { s, s, s };
                    data->density_volume.model_mat = volume::compute_model_to_world_matrix(min_aabb, max_aabb);
                    data->density_volume.voxel_spacing = vec3_t{2*s / prop->data.dim[0], 2*s / prop->data.dim[1], 2*s / prop->data.dim[2]};
                }
                num_reps = vis.sdf.count;
            }

            // We need to limit this for performance reasons
            num_reps = MIN(num_reps, 100);

            const int64_t old_size = md_array_size(data->density_volume.gl_reps);
            if (data->density_volume.gl_reps) {
                // Only free superflous entries
                for (int64_t i = num_reps; i < old_size; ++i) {
                    md_gl_representation_free(&data->density_volume.gl_reps[i]);
                }
            }
            md_array_resize(data->density_volume.gl_reps, num_reps, persistent_allocator);
            md_array_resize(data->density_volume.rep_model_mats, num_reps, persistent_allocator);

            for (int64_t i = old_size; i < num_reps; ++i) {
                // Only init new entries
                md_gl_representation_init(&data->density_volume.gl_reps[i], &data->mold.gl_mol);
            }

            const auto& mol = data->mold.mol;
            auto& rep = data->density_volume.rep;
            const int64_t num_colors = data->mold.mol.atom.count;
            uint32_t* colors = (uint32_t*)md_alloc(frame_allocator, sizeof(uint32_t) * num_colors);
            md_gl_representation_type_t rep_type = MD_GL_REP_SPACE_FILL;
            md_gl_representation_args_t rep_args = {};

            switch (rep.colormap) {
            case ColorMapping::Uniform:
                color_atoms_uniform(colors, mol.atom.count, rep.color);
                break;
            case ColorMapping::Cpk:
                color_atoms_cpk(colors, mol.atom.count, mol);
                break;
            case ColorMapping::AtomLabel:
                color_atoms_label(colors, mol.atom.count, mol);
                break;
            case ColorMapping::AtomIndex:
                color_atoms_idx(colors, mol.atom.count, mol);
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
            default:
                ASSERT(false);
                break;
            }

            switch (rep.type) {
            case RepresentationType::SpaceFill:
                rep_type = MD_GL_REP_SPACE_FILL;
                rep_args.space_fill.radius_scale = rep.radius;
                break;
            case RepresentationType::Licorice:
                rep_type = MD_GL_REP_LICORICE;
                rep_args.licorice.radius = rep.radius * 0.5f;
                break;
            case RepresentationType::Ribbons:
                rep_type = MD_GL_REP_RIBBONS;
                rep_args.ribbons.width_scale = rep.width;
                rep_args.ribbons.thickness_scale = rep.thickness;
                break;
            case RepresentationType::Cartoon:
                rep_type = MD_GL_REP_CARTOON;
                rep_args.cartoon.width_scale = rep.width;
                rep_args.cartoon.thickness_scale = rep.thickness;
                break;
            default:
                ASSERT(false);
            }

            for (int64_t i = 0; i < num_reps; ++i) {
                filter_colors(colors, num_colors, &vis.sdf.structures[i]);
                md_gl_representation_set_color(&data->density_volume.gl_reps[i], 0, num_colors, colors, 0);
                md_gl_representation_set_type_and_args(&data->density_volume.gl_reps[i], rep_type, rep_args);
                data->density_volume.rep_model_mats[i] = vis.sdf.matrices[i];
            }
        }
    }

    if (data->density_volume.dirty_vol) {
        if (prop) {
            data->density_volume.dirty_vol = false;
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
}

static void clear_density_volume(ApplicationData* data) {
    md_array_shrink(data->density_volume.gl_reps, 0);
    md_array_shrink(data->density_volume.rep_model_mats, 0);
    data->density_volume.model_mat = {0};
}

static void interpolate_atomic_properties(ApplicationData* data) {
    ASSERT(data);
    const auto& mol = data->mold.mol;
    const auto& traj = data->mold.traj;

    if (!mol.atom.count || !md_trajectory_num_frames(traj)) return;

    const int64_t last_frame = MAX(0LL, md_trajectory_num_frames(traj) - 1);
    // This is not actually time, but the fractional frame representation
    const double time = CLAMP(data->animation.frame, 0.0, double(last_frame));

    // Scaling factor for cubic spline
    const float s = 1.0f - CLAMP(data->animation.tension, 0.0f, 1.0f);
    const float t = (float)fractf(time);
    const int64_t frame = (int64_t)time;
    const int64_t nearest_frame = CLAMP((int64_t)(time + 0.5), 0LL, last_frame);

    const int64_t frames[4] = {
        MAX(0LL, frame - 1),
        MAX(0LL, frame),
        MIN(frame + 1, last_frame),
        MIN(frame + 2, last_frame)
    };

    int64_t stride = ALIGN_TO(mol.atom.count, md_simd_f32_width);    // The interploation uses SIMD vectorization without bounds, so we make sure there is no overlap between the data segments
    int64_t bytes = stride * sizeof(float) * 3 * 4;
    void* mem = md_alloc(frame_allocator, bytes);
    defer { md_free(frame_allocator, mem, bytes); };

    md_vec3_soa_t src[4] = {
        {(float*)mem + stride * 0, (float*)mem + stride *  1, (float*)mem + stride *  2},
        {(float*)mem + stride * 3, (float*)mem + stride *  4, (float*)mem + stride *  5},
        {(float*)mem + stride * 6, (float*)mem + stride *  7, (float*)mem + stride *  8},
        {(float*)mem + stride * 9, (float*)mem + stride * 10, (float*)mem + stride * 11},
    };

    md_vec3_soa_t dst = {
        data->mold.mol.atom.x, data->mold.mol.atom.y, data->mold.mol.atom.z,
    };

    const InterpolationMode mode = (frames[1] != frames[2]) ? data->animation.interpolation : InterpolationMode::Nearest;
    switch (mode) {
        case InterpolationMode::Nearest:
        {
            md_trajectory_frame_header_t header = {0};
            md_trajectory_load_frame(data->mold.traj, nearest_frame, &header, mol.atom.x, mol.atom.y, mol.atom.z);
            data->simulation_box.box = header.box;
            break;
        }
        case InterpolationMode::Linear:
        {
            md_trajectory_frame_header_t header[2] = {0};
            md_trajectory_load_frame(data->mold.traj, frames[1], &header[0], src[0].x, src[0].y, src[0].z);
            md_trajectory_load_frame(data->mold.traj, frames[2], &header[1], src[1].x, src[1].y, src[1].z);
            data->simulation_box.box = lerp(header[0].box, header[1].box, t);
            const vec3_t pbc_ext = data->simulation_box.box * vec3_set1(1);

            md_util_linear_interpolation(dst, src, mol.atom.count, pbc_ext, t);
        }
            break;
        case InterpolationMode::CubicSpline:
        {
            md_trajectory_frame_header_t header[4] = {0};
            md_trajectory_load_frame(data->mold.traj, frames[0], &header[0], src[0].x, src[0].y, src[0].z);
            md_trajectory_load_frame(data->mold.traj, frames[1], &header[1], src[1].x, src[1].y, src[1].z);
            md_trajectory_load_frame(data->mold.traj, frames[2], &header[2], src[2].x, src[2].y, src[2].z);
            md_trajectory_load_frame(data->mold.traj, frames[3], &header[3], src[3].x, src[3].y, src[3].z);
            data->simulation_box.box = cubic_spline(header[0].box, header[1].box, header[2].box, header[3].box, t, s);
            const vec3_t pbc_ext = data->simulation_box.box * vec3_set1(1);

            md_util_cubic_spline_interpolation(dst, src, mol.atom.count, pbc_ext, t, s);
        }
            break;
        default:
            ASSERT(false);
    }

    md_util_compute_aabb_xyzr(&data->mold.mol_aabb_min, &data->mold.mol_aabb_max, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.radius, mol.atom.count);

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
                float phi[2] = {src_angles[1][i].phi, src_angles[2][i].phi};
                float psi[2] = {src_angles[1][i].psi, src_angles[2][i].psi};

                phi[1] = deperiodizef(phi[1], phi[0], (float)TWO_PI);
                psi[1] = deperiodizef(psi[1], psi[0], (float)TWO_PI);

                float final_phi = lerp(phi[0], phi[1], t);
                float final_psi = lerp(psi[0], psi[1], t);
                mol.backbone.angle[i] = {deperiodizef(final_phi, 0, (float)TWO_PI), deperiodizef(final_psi, 0, (float)TWO_PI)};
            }
            break;
        }
        case InterpolationMode::CubicSpline: {
            for (int64_t i = 0; i < mol.backbone.count; ++i) {
                float phi[4] = {src_angles[0][i].phi, src_angles[1][i].phi, src_angles[2][i].phi, src_angles[3][i].phi};
                float psi[4] = {src_angles[0][i].psi, src_angles[1][i].psi, src_angles[2][i].psi, src_angles[3][i].psi};

                phi[0] = deperiodizef(phi[0], phi[1], (float)TWO_PI);
                phi[2] = deperiodizef(phi[2], phi[1], (float)TWO_PI);
                phi[3] = deperiodizef(phi[3], phi[2], (float)TWO_PI);

                psi[0] = deperiodizef(psi[0], psi[1], (float)TWO_PI);
                psi[2] = deperiodizef(psi[2], psi[1], (float)TWO_PI);
                psi[3] = deperiodizef(psi[3], psi[2], (float)TWO_PI);

                float final_phi = cubic_spline(phi[0], phi[1], phi[2], phi[3], t, s);
                float final_psi = cubic_spline(psi[0], psi[1], psi[2], psi[3], t, s);
                mol.backbone.angle[i] = {deperiodizef(final_phi, 0, (float)TWO_PI), deperiodizef(final_psi, 0, (float)TWO_PI)};
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
        case InterpolationMode::CubicSpline: {
            for (int64_t i = 0; i < mol.backbone.count; ++i) {
                const vec4_t ss_f[4] = {
                    convert_color((uint32_t)src_ss[0][i]),
                    convert_color((uint32_t)src_ss[1][i]),
                    convert_color((uint32_t)src_ss[2][i]),
                    convert_color((uint32_t)src_ss[3][i]),
                };
                const vec4_t ss_res = cubic_spline(ss_f[0], ss_f[1], ss_f[2], ss_f[3], t, s);
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
static void update_view_param(ApplicationData* data) {
    ViewParam& param = data->view.param;
    param.matrix.previous = param.matrix.current;
    param.jitter.previous = param.jitter.current;

    param.clip_planes.near = data->view.camera.near_plane;
    param.clip_planes.far = data->view.camera.far_plane;
    param.fov_y = data->view.camera.fov_y;
    param.resolution = {(float)data->gbuffer.width, (float)data->gbuffer.height};

    param.matrix.current.view = camera_world_to_view_matrix(data->view.camera);
    param.matrix.inverse.view = camera_view_to_world_matrix(data->view.camera);

    if (data->view.mode == CameraMode::Perspective) {
        param.matrix.current.proj = camera_perspective_projection_matrix(data->view.camera, (float)data->gbuffer.width / (float)data->gbuffer.height);
        param.matrix.inverse.proj = camera_inverse_perspective_projection_matrix(data->view.camera, (float)data->gbuffer.width / (float)data->gbuffer.height);
    } else {
        const float aspect_ratio = (float)data->gbuffer.width / (float)data->gbuffer.height;
        const float h = data->view.camera.focus_distance * tanf(data->view.camera.fov_y * 0.5f);
        const float w = aspect_ratio * h;
        const float n = data->view.camera.near_plane;
        const float f = data->view.camera.far_plane;
        param.matrix.current.proj = camera_orthographic_projection_matrix(-w, w, -h, h, n, f);
        param.matrix.inverse.proj = camera_inverse_orthographic_projection_matrix(-w, w, -h, h, n, f);
    }
    param.matrix.current.proj_jittered = param.matrix.current.proj;

    if (data->visuals.temporal_reprojection.enabled && data->visuals.temporal_reprojection.jitter) {
        static uint32_t i = 0;
        i = (i+1) % ARRAY_SIZE(data->view.jitter.sequence);
        param.jitter.next    = data->view.jitter.sequence[(i + 1) % ARRAY_SIZE(data->view.jitter.sequence)] - 0.5f;
        param.jitter.current = data->view.jitter.sequence[i] - 0.5f;
        if (data->view.mode == CameraMode::Perspective) {
            const vec2_t j = param.jitter.current;
            const int w = data->gbuffer.width;
            const int h = data->gbuffer.height;
            param.matrix.current.proj_jittered = camera_perspective_projection_matrix(data->view.camera, w, h, j.x, j.y);
            param.matrix.inverse.proj_jittered = camera_inverse_perspective_projection_matrix(data->view.camera, w, h, j.x, j.y);
        } else {
            const float aspect_ratio = (float)data->gbuffer.width / (float)data->gbuffer.height;
            const float h = data->view.camera.focus_distance * tanf(data->view.camera.fov_y * 0.5f);
            const float w = aspect_ratio * h;
            const vec2_t scl = {w / data->gbuffer.width * 2.0f, h / data->gbuffer.height * 2.0f};
            const vec2_t j = param.jitter.current * scl;
            param.matrix.current.proj_jittered = camera_orthographic_projection_matrix(-w + j.x, w + j.x, -h + j.y, h + j.y, data->view.camera.near_plane, data->view.camera.far_plane);
            param.matrix.inverse.proj_jittered = camera_inverse_orthographic_projection_matrix(-w + j.x, w + j.x, -h + j.y, h + j.y, data->view.camera.near_plane, data->view.camera.far_plane);
        }
    }

    param.matrix.current.norm = mat4_transpose(param.matrix.inverse.view);
}

static void reset_view(ApplicationData* data, bool move_camera, bool smooth_transition) {
    ASSERT(data);
    if (!data->mold.mol.atom.count) return;
    const auto& mol = data->mold.mol;

    vec3_t aabb_min, aabb_max;
    if (data->simulation_box.box != mat3_t{0}) {
        aabb_min = {0,0,0};
        aabb_max = data->simulation_box.box * vec3_t{1,1,1};
    } else {
        compute_aabb(&aabb_min, &aabb_max, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count);
    }

    const vec3_t ext = aabb_max - aabb_min;
    const vec3_t cen = (aabb_min + aabb_max) * 0.5f;
    const vec3_t pos = cen + ext * 5.f;

    if (move_camera) {
        const quat_t ori = quat_from_mat4(look_at(pos, cen, {0, 1, 0}));
        const float dist = vec3_length(pos - cen);

        data->view.animation.target_position = pos;
        data->view.animation.target_orientation = ori;
        data->view.animation.target_distance = dist;

        if (!smooth_transition) {
            data->view.camera.position = pos;
            data->view.camera.orientation = ori;
            data->view.camera.focus_distance = dist;
        }
    }

    data->view.camera.near_plane = 1.0f;
    data->view.camera.far_plane = 10000.0f;
    data->view.trackball_param.max_distance = vec3_length(ext) * 10.0f;
}

// #picking
static PickingData read_picking_data(GBuffer* gbuf, int32_t x, int32_t y) {
    PickingData data{};
#if EXPERIMENTAL_GFX_API
    if (use_gfx) {
        data.idx = md_gfx_get_picking_idx();
        data.depth = md_gfx_get_picking_depth();
        md_gfx_query_picking((uint32_t)x, (uint32_t)y);
    }
    else {
#endif
    ASSERT(gbuf);
    uint32_t N = (uint32_t)ARRAY_SIZE(gbuf->pbo_picking.color);
    uint32_t frame = gbuf->pbo_picking.frame++;
    uint32_t queue = (frame) % N;
    uint32_t read  = (frame + N-1) % N;


    PUSH_GPU_SECTION("READ PICKING DATA")
    glBindFramebuffer(GL_READ_FRAMEBUFFER, gbuf->deferred.fbo);
    glReadBuffer(GL_COLOR_ATTACHMENT_PICKING);

    // Queue async reads from current frame to pixel pack buffer
    glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.color[queue]);
    glReadPixels(x, y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.depth[queue]);
    glReadPixels(x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, 0);

    // Read values from previous frames pixel pack buffer
    glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.color[read]);
    const GLubyte* color = (const GLubyte*)glMapBuffer(GL_PIXEL_PACK_BUFFER, GL_READ_ONLY);
    if (color) {
        data.idx = color[0] + (color[1] << 8) + (color[2] << 16) + (color[3] << 24);
        glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
    }

    glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.depth[read]);
    const GLfloat* depth = (const GLfloat*)glMapBuffer(GL_PIXEL_PACK_BUFFER, GL_READ_ONLY);
    if (depth) {
        data.depth = depth[0];
        glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
    }

    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    POP_GPU_SECTION()
#if EXPERIMENTAL_GFX_API
    }
#endif
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

static void grow_mask_by_covalent_bond(md_bitfield_t* mask, md_bond_t* bonds, int64_t num_bonds, int64_t extent) {
    md_bitfield_t prev_mask;
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

static void grow_mask_by_radial_extent(md_bitfield_t* dst_mask, const md_bitfield_t* src_mask, const md_molecule_t& mol, vec3_t pbc_ext, float rad_ext) {
    if (rad_ext > 0.0f) {
        md_spatial_hash_t ctx = {0};
        md_spatial_hash_init_soa(&ctx, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, pbc_ext, frame_allocator);
        defer { md_spatial_hash_free(&ctx); };

        md_bitfield_clear(dst_mask);

        int64_t beg_bit = src_mask->beg_bit;
        int64_t end_bit = src_mask->end_bit;
        while ((beg_bit = md_bitfield_scan(src_mask, beg_bit, end_bit)) != 0) {
            int64_t i = beg_bit - 1;
            vec3_t pos = {mol.atom.x[i], mol.atom.y[i], mol.atom.z[i]};
            float  rad = rad_ext;
            md_spatial_hash_query(&ctx, pos, rad, [](uint32_t idx, vec3_t, void* user_data) -> bool {
                md_bitfield_t* mask = (md_bitfield_t*)user_data;
                md_bitfield_set_bit(mask, (int64_t)idx);
                return true;
            }, dst_mask);
        }
    } else {
        md_bitfield_copy(dst_mask, src_mask);
    }
}

static void expand_mask(md_bitfield_t* mask, const md_range_t ranges[], int64_t num_ranges) {
    for (int64_t i = 0; i < num_ranges; i++) {
        if (md_bitfield_popcount_range( mask, ranges[i].beg, ranges[i].end) != 0) {
            md_bitfield_set_range(mask, ranges[i].beg, ranges[i].end);
        }
    }
}

static void grow_mask_by_current_selection_granularity(md_bitfield_t* mask, const ApplicationData& data) {
    ASSERT(mask);
    switch(data.selection.granularity) {
    case SelectionLevel::Atom:
        break;
    case SelectionLevel::Residue:
        expand_mask(mask, data.mold.mol.residue.atom_range, data.mold.mol.residue.count);
        break;
    case SelectionLevel::Chain:
        expand_mask(mask, data.mold.mol.chain.atom_range, data.mold.mol.chain.count);
        break;
    default:
        ASSERT(false);
    }
}

static bool filter_expression(ApplicationData* data, str_t expr, md_bitfield_t* mask, bool* is_dynamic = NULL, char* error_buf = NULL, int error_cap = 0) {
    if (data->mold.mol.atom.count == 0) return false;    
    
    bool success = false;

    if (md_semaphore_aquire(&data->mold.script.ir_semaphore)) {
        defer { md_semaphore_release(&data->mold.script.ir_semaphore); };
        return md_filter(mask, expr, &data->mold.mol, data->mold.script.ir, is_dynamic, error_buf, error_cap);
    }

    return success;
}

// ### DRAW WINDOWS ###
static void draw_main_menu(ApplicationData* data) {
    ASSERT(data);
    bool new_clicked = false;
    char path_buf[4096] = "";

    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("Load Data", "CTRL+L")) {
                data->load_dataset.show_window = true;
            }
            if (ImGui::MenuItem("Open Workspace", "CTRL+O")) {
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialog_Open, FILE_EXTENSION)) {
                    load_workspace(data, str_from_cstr(path_buf));
                }
            }
            if (ImGui::MenuItem("Save Workspace", "CTRL+S")) {
                if (!data->files.workspace) {
                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialog_Save, FILE_EXTENSION)) {
                        int path_len = (int)strnlen(path_buf, sizeof(path_buf));
                        str_t ext = extract_ext({path_buf, path_len});
                        if (str_empty(ext)) {
                            path_len += snprintf(path_buf + path_len, sizeof(path_buf) - path_len, ".%s", FILE_EXTENSION);
                        }
                        save_workspace(data, {path_buf, path_len});
                    }
                } else {
                    save_workspace(data, data->files.workspace);
                }
            }
            if (ImGui::MenuItem("Save As")) {
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialog_Save, FILE_EXTENSION)) {
                    int path_len = (int)strnlen(path_buf, sizeof(path_buf));
                    str_t ext = extract_ext({path_buf, path_len});
                    if (str_empty(ext)) {
                        path_len += snprintf(path_buf + path_len, sizeof(path_buf) - path_len, ".%s", FILE_EXTENSION);
                    }
                    save_workspace(data, {path_buf, path_len});
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
            ImGui::Checkbox("Representations", &data->representation.show_window);
            ImGui::Checkbox("Script", &data->show_script_window);
            ImGui::Checkbox("Timelines", &data->timeline.show_window);
            ImGui::Checkbox("Distributions", &data->distributions.show_window);
            ImGui::Checkbox("Density Volumes", &data->density_volume.show_window);
            ImGui::Checkbox("Ramachandran", &data->ramachandran.show_window);
            ImGui::Checkbox("Shape Space", &data->shape_space.show_window);
            ImGui::Checkbox("Dataset", &data->dataset.show_window);

            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Selection")) {
            ImGui::Combo("Granularity", (int*)(&data->selection.granularity), "Atom\0Residue\0Chain\0\0");
            int64_t num_selected_atoms = md_bitfield_popcount(&data->selection.current_selection_mask);
            if (ImGui::MenuItem("Invert")) {
                md_bitfield_not_inplace(&data->selection.current_selection_mask, 0, data->mold.mol.atom.count);
                data->mold.dirty_buffers |= MolBit_DirtyFlags;
            }
            if (ImGui::IsItemHovered()) {
                md_bitfield_not(&data->selection.current_highlight_mask, &data->selection.current_selection_mask, 0, data->mold.mol.atom.count);
                data->mold.dirty_buffers |= MolBit_DirtyFlags;
            }
            if (ImGui::MenuItem("Query")) data->selection.query.show_window = true;
            if (num_selected_atoms == 0) ImGui::PushDisabled();
            if (ImGui::MenuItem("Grow"))  data->selection.grow.show_window = true;
            if (num_selected_atoms == 0) ImGui::PopDisabled();

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
                    bool is_valid = md_script_identifier_name_valid(sel.name);
                    char error[64] = "";
                    if (!is_valid) {
                        snprintf(error, sizeof(error), "'%s' is not a valid identifier.", sel.name.cstr());
                    }

                    for (int j = 0; j < i; ++j) {
                        if (str_equal(sel.name, data->selection.stored_selections[j].name)) {
                            is_valid = false;
                            snprintf(error, sizeof(error), "identifier '%s' is already taken.", sel.name.cstr());
                            break;
                        }
                    }

                    ImGui::PushID(i);
                    ImGui::InputQuery("##label", sel.name.beg(), sel.name.capacity(), is_valid, error);
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
                        data->mold.script.compile_ir = true;
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
                    snprintf(name_buf, sizeof(name_buf), "sel%i", (int)md_array_size(data->selection.stored_selections) + 1);
                    create_selection(data, str_from_cstr(name_buf), &data->selection.current_selection_mask);
                }

                //ImGui::PopItemFlag();
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Screenshot")) {
            //ImGui::Checkbox("Hide GUI", &data->screenshot.hide_gui);
            if (ImGui::MenuItem("Take Screenshot")) {
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialog_Save, "jpg;png;bmp")) {
                    int path_len = (int)strnlen(path_buf, sizeof(path_buf));
                    str_t ext = extract_ext({path_buf, path_len});
                    if (str_empty(ext)) {
                        path_len += snprintf(path_buf + path_len, sizeof(path_buf) - path_len, ".jpg");
                        ext = STR("jpg");
                    }
                    if (str_equal_cstr_ignore_case(ext, "jpg") || str_equal_cstr_ignore_case(ext, "png") || str_equal_cstr_ignore_case(ext, "bmp")) {
                        data->screenshot.path_to_file = str_copy({path_buf, path_len}, frame_allocator);
                    }
                    else {
                        LOG_ERROR("Supplied image extension is not supported");
                    }
                }
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Settings")) {
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

            double ms = (data->ctx.timing.delta_s * 1000);
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

void draw_load_dataset_window(ApplicationData* data) {
    LoadDatasetWindowState& state = data->load_dataset;
    
    ImGui::SetNextWindowSize(ImVec2(300, 200), ImGuiCond_FirstUseEver);
    if (state.show_window) {
        ImGui::OpenPopup("Load Dataset");
    }

    if (ImGui::BeginPopupModal("Load Dataset", &state.show_window)) {
        bool path_changed = false;
        bool path_invalid = !state.path_is_valid && state.path_buf[0] != '\0';
        const str_t* loader_ext = load::get_supported_extensions();
        const int loader_ext_count = (int)load::get_supported_extension_count();

        if (path_invalid) ImGui::PushInvalid();
        if (ImGui::InputText("##path", state.path_buf, sizeof(state.path_buf))) {
            path_changed = true;
        }
        if (path_invalid) ImGui::PopInvalid();

        ImGui::SameLine();
        if (ImGui::Button("Browse")) {
            if (application::file_dialog(state.path_buf, sizeof(state.path_buf), application::FileDialog_Open)) {
                path_changed = true;
            }
        }

        str_t path = str_from_cstr(state.path_buf);
        if (path_changed) {
            state.path_is_valid = md_path_is_valid(path) && !md_path_is_directory(path);
            str_t ext = extract_ext(path);

            // Try to assign loader_idx from extension
            state.loader_idx = -1;
            for (int i = 0; i < loader_ext_count; ++i) {
                if (str_equal_ignore_case(ext, loader_ext[i])) {
                    state.loader_idx = i;
                    break;
                }
            }
        }

        if (ImGui::BeginCombo("Loader", state.loader_idx > -1 ? loader_ext[state.loader_idx].ptr : "")) {
            for (int i = 0; i < loader_ext_count; ++i) {
                if (ImGui::Selectable(loader_ext[i].ptr, state.loader_idx == i)) {
                    state.loader_idx = i;
                }
            }
            ImGui::EndCombo();
        }

        bool show_cg = state.path_is_valid && load::mol::get_api(path);
        if (show_cg) {
            ImGui::Checkbox("Coarse Grained", &state.coarse_grained);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Enable if the dataset is coarse grained");
            }
        }

        bool show_dp = state.path_is_valid && load::traj::get_api(path);
        if (show_dp) {
            ImGui::Checkbox("Deperiodize on Load", &state.deperiodize_on_load);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Enable if the loaded frames should be deperiodized");
            }
        }

        bool load_enabled = (state.path_is_valid && state.loader_idx > -1);
        if (!load_enabled) ImGui::PushDisabled();
        if (ImGui::Button("Load")) {
            // A bit cheeky, write a temp string with .[extension] to query for loaders
            char buf[32];
            int len = snprintf(buf, sizeof(buf), ".%s", loader_ext[state.loader_idx].ptr);
            str_t ext = {buf, len};
            if (load_dataset_from_file(data, path, load::mol::get_api(ext), load::traj::get_api(ext), show_cg && state.coarse_grained, show_dp && state.deperiodize_on_load)) {
                if (md_array_size(data->representation.reps) == 0) {
                    create_representation(data); // Create default representation
                }
                data->animation = {};
                reset_view(data, true, true);
            }
            state = LoadDatasetWindowState();
            state.show_window = false;
        }
        if (!load_enabled) ImGui::PopDisabled();

        ImGui::SameLine();
        if (ImGui::Button("Cancel")) {
            state.show_window = false;
        }
        ImGui::EndPopup();
    }
}

void clear_atom_elem_mappings(ApplicationData* data) {
    md_array_shrink(data->dataset.atom_element_remappings, 0);
}

AtomElementMapping* add_atom_elem_mapping(ApplicationData* data, str_t lbl, md_element_t elem) {
    // Check if we already have a mapping for the label -> overwrite
    int64_t i = 0;
    for (; i < md_array_size(data->dataset.atom_element_remappings); ++i) {
        if (str_equal(lbl, data->dataset.atom_element_remappings[i].lbl)) break;
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
    if (data->mold.mol.atom.count == 0 || !data->mold.mol.atom.element) {
        return;
    }

    for (int64_t j = 0; j < md_array_size(data->dataset.atom_element_remappings); ++j) {
        str_t lbl = data->dataset.atom_element_remappings[j].lbl;
        md_element_t elem = data->dataset.atom_element_remappings[j].elem;
        float radius = md_util_element_vdw_radius(elem);

        for (int64_t i = 0; i < data->mold.mol.atom.count; ++i) {
            if (str_equal(lbl, data->mold.mol.atom.name[i])) {
                data->mold.mol.atom.element[i] = elem;
                data->mold.mol.atom.radius[i] = radius;
                data->mold.dirty_buffers |= MolBit_DirtyRadius;
            }
        }
    }
    md_molecule_t* mol = &data->mold.mol;
    
    const vec3_t pbc_ext = md_util_compute_unit_cell_extent(mol->coord_frame);
    md_bond_data_clear(&mol->covalent);
    md_util_compute_covalent_bonds_and_connectivity(&mol->covalent, mol->atom.x, mol->atom.y, mol->atom.z, mol->atom.element, mol->atom.residue_idx, mol->atom.count, pbc_ext, data->mold.mol_alloc);
    data->mold.dirty_buffers |= MolBit_DirtyBonds;

    update_all_representations(data);
}

// Create a textual script describing a selection from a bitfield with respect to some reference index
static void create_script_ranges(md_strb_t* sb, const md_bitfield_t* bf, int ref_idx = 0) {
    ASSERT(sb);
    ASSERT(bf);

    md_strb_fmt(sb, "atom({");

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
        md_strb_fmt(sb, "%i:%i})", range_beg - ref_idx + 1, prev_idx - ref_idx + 1);
    } else if (prev_idx != -1) {
        md_strb_fmt(sb, "%i})", prev_idx - ref_idx + 1);
    }
}

static int64_t find_identifier(const md_script_ir_t* ir, str_t ident) {
    const int64_t num_ident = md_script_ir_num_identifiers(ir);
    const str_t* idents = md_script_ir_identifiers(ir);
    for (int64_t i = 0; i < num_ident; ++i) {
        if (str_equal(ident, idents[i])) return i;
    }
    return -1;
}

static str_t create_unique_identifier(const md_script_ir_t* ir, str_t base, md_allocator_i* alloc) {
    char buf[128];
    for (int64_t i = 1; i < 10; ++i) {
        str_t ident = {buf, snprintf(buf, sizeof(buf), "%.*s%i", (int)base.len, base.ptr, (int)i)};
        if (find_identifier(ir, ident) == -1) {
            return str_copy(ident, alloc);
        }
    }
    return str_t();
}

// # context_menu
void draw_context_popup(ApplicationData* data) {
    ASSERT(data);

    if (!data->mold.mol.atom.count) return;

    const int64_t sss_count = single_selection_sequence_count(&data->selection.single_selection_sequence);
    const int64_t num_frames = md_trajectory_num_frames(data->mold.traj);
    const int64_t num_atoms_selected = md_bitfield_popcount(&data->selection.current_selection_mask);

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
            if (num_atoms_selected <= 4 && sss_count > 1) {
                char buf[256] = "";

                if (sss_count == 2) {
                    int32_t idx[2] = {data->selection.single_selection_sequence.idx[0], data->selection.single_selection_sequence.idx[1]};
                    str_t ident = create_unique_identifier(data->mold.script.ir, STR("dist"), frame_allocator);

                    snprintf(buf, sizeof(buf), "%.*s = distance(%i, %i);", (int)ident.len, ident.ptr, idx[0]+1, idx[1]+1);
                    if (ImGui::MenuItem(buf)) {
                        editor.AppendText("\n");
                        editor.AppendText(buf);
                        ImGui::CloseCurrentPopup();
                    }

                    if (data->mold.mol.residue.count) {
                        if (data->mold.mol.atom.residue_idx[idx[0]] == data->mold.mol.atom.residue_idx[idx[1]]) {
                            int32_t res_idx = data->mold.mol.atom.residue_idx[idx[0]];
                            idx[0] -= data->mold.mol.residue.atom_range[res_idx].beg;
                            idx[1] -= data->mold.mol.residue.atom_range[res_idx].beg;

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
                    str_t ident = create_unique_identifier(data->mold.script.ir, STR("ang"), frame_allocator);

                    snprintf(buf, sizeof(buf), "%.*s = angle(%i, %i, %i);", (int)ident.len, ident.ptr, idx[0]+1, idx[1]+1, idx[2]+1);
                    if (ImGui::MenuItem(buf)) {
                        editor.AppendText("\n");
                        editor.AppendText(buf);
                        ImGui::CloseCurrentPopup();
                    }

                    if (data->mold.mol.residue.count) {
                        if (data->mold.mol.atom.residue_idx[idx[0]] == data->mold.mol.atom.residue_idx[idx[1]] &&
                            data->mold.mol.atom.residue_idx[idx[0]] == data->mold.mol.atom.residue_idx[idx[2]]) {

                            int32_t res_idx = data->mold.mol.atom.residue_idx[idx[0]];
                            idx[0] -= data->mold.mol.residue.atom_range[res_idx].beg;
                            idx[1] -= data->mold.mol.residue.atom_range[res_idx].beg;
                            idx[2] -= data->mold.mol.residue.atom_range[res_idx].beg;

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
                    str_t ident = create_unique_identifier(data->mold.script.ir, STR("dih"), frame_allocator);

                    snprintf(buf, sizeof(buf), "%.*s = dihedral(%i, %i, %i, %i);", (int)ident.len, ident.ptr, idx[0]+1, idx[1]+1, idx[2]+1, idx[3]+1);
                    if (ImGui::MenuItem(buf)) {
                        editor.AppendText("\n");
                        editor.AppendText(buf);
                        ImGui::CloseCurrentPopup();
                    }

                    if (data->mold.mol.residue.count) {
                        if (data->mold.mol.atom.residue_idx[idx[0]] == data->mold.mol.atom.residue_idx[idx[1]] &&
                            data->mold.mol.atom.residue_idx[idx[0]] == data->mold.mol.atom.residue_idx[idx[2]] &&
                            data->mold.mol.atom.residue_idx[idx[0]] == data->mold.mol.atom.residue_idx[idx[3]]) {

                            int32_t res_idx = data->mold.mol.atom.residue_idx[idx[0]];
                            idx[0] -= data->mold.mol.residue.atom_range[res_idx].beg;
                            idx[1] -= data->mold.mol.residue.atom_range[res_idx].beg;
                            idx[2] -= data->mold.mol.residue.atom_range[res_idx].beg;
                            idx[3] -= data->mold.mol.residue.atom_range[res_idx].beg;

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
            if (num_atoms_selected >= 3) {
                const md_bitfield_t* bf = &data->selection.current_selection_mask;
                str_t ident = create_unique_identifier(data->mold.script.ir, STR("sel"), frame_allocator);
                bool same_residue = true;
                bool same_chain   = true;

                md_chain_idx_t chain_idx = -1;
                md_residue_idx_t res_idx = -1;

                uint64_t beg_idx = md_bitfield_beg_bit(bf);
                uint64_t end_idx = md_bitfield_end_bit(bf);

                while ((beg_idx = md_bitfield_scan(bf, beg_idx, end_idx)) != 0) {
                    uint64_t i = beg_idx - 1;
                    if (data->mold.mol.atom.chain_idx) {
                        if (chain_idx == -1) {
                            chain_idx = data->mold.mol.atom.chain_idx[i];
                        } else if (data->mold.mol.atom.chain_idx[i] != chain_idx) {
                            same_chain = false;
                        }
                    }

                    if (data->mold.mol.atom.residue_idx) {
                        if (res_idx == -1) {
                            res_idx = data->mold.mol.atom.residue_idx[i];
                        } else if (data->mold.mol.atom.residue_idx[i] != res_idx) {
                            same_residue = false;                            
                        }
                    } 

                    if (!same_residue && !same_chain) {
                        break;
                    }
                }
                md_strb_t sb = {0};
                md_strb_init(&sb, frame_allocator);
                defer { md_strb_free(&sb); };

                if (res_idx != -1 && same_residue) {
                    md_strb_reset(&sb);
                    md_strb_str(&sb, ident);
                    md_strb_str(&sb, STR(" = "));
                    create_script_ranges(&sb, bf, data->mold.mol.residue.atom_range[res_idx].beg);
                    if (md_strb_len(&sb) < 512) {
                        str_t resname = LBL_TO_STR(data->mold.mol.residue.name[res_idx]);
                        md_strb_fmt(&sb, " in resname(\"%.*s\");", (int)resname.len, resname.ptr);
                        if (ImGui::MenuItem(sb.buf)) {
                            editor.AppendText("\n");
                            editor.AppendText(sb.buf);
                            ImGui::CloseCurrentPopup();
                        }
                    }
                }

                if (chain_idx != -1 && same_chain) {
                    md_strb_reset(&sb);
                    md_strb_str(&sb, ident);
                    md_strb_str(&sb, STR(" = "));
                    create_script_ranges(&sb, bf, data->mold.mol.chain.atom_range[chain_idx].beg);
                    if (md_strb_len(&sb) < 512) {
                        str_t chain_id = LBL_TO_STR(data->mold.mol.chain.id[chain_idx]);
                        md_strb_fmt(&sb, " in chain(\"%.*s\");", (int)chain_id.len, chain_id.ptr);
                        if (ImGui::MenuItem(sb.buf)) {
                            editor.AppendText("\n");
                            editor.AppendText(sb.buf);
                            ImGui::CloseCurrentPopup();
                        }
                    }
                }

                md_strb_reset(&sb);
                md_strb_str(&sb, ident);
                md_strb_str(&sb, STR(" = "));
                create_script_ranges(&sb, bf);
                if (md_strb_len(&sb) < 512) {
                    md_strb_char(&sb, ';');
                    if (ImGui::MenuItem(sb.buf)) {
                        editor.AppendText("\n");
                        editor.AppendText(sb.buf);
                        ImGui::CloseCurrentPopup();
                    }
                }
                md_strb_free(&sb);
            }
            ImGui::EndMenu();
        }

        if (data->selection.right_clicked != -1 && num_atoms_selected == 0 && data->mold.mol.atom.element) {
            int idx = data->selection.right_clicked;
            if (0 <= idx && idx < data->mold.mol.atom.count) {
                char label[64] = "";
                str_t atom_name = data->mold.mol.atom.name[idx];
                snprintf(label, sizeof(label), "Remap Element for '%.*s'", (int)atom_name.len, atom_name.ptr);
                if (ImGui::BeginMenu(label)) {
                    static char input_buf[32] = "";
                    md_element_t elem = data->mold.mol.atom.element[idx];
                    str_t name = md_util_element_name(elem);
                    str_t sym  = md_util_element_symbol(elem);

                    ImGui::Text("Current Element: %.*s (%.*s)", (int)name.len, name.ptr, (int)sym.len, sym.ptr);

                    str_t elem_str = {input_buf, (int64_t)strnlen(input_buf, sizeof(input_buf))};
                    md_element_t new_elem = md_util_element_lookup(elem_str);
                    const bool is_valid = new_elem != 0;

                    ImGui::InputQuery("##Symbol", input_buf, sizeof(input_buf), is_valid, "Cannot recognize Element symbol");
                    str_t new_name = md_util_element_name(new_elem);
                    str_t new_sym  = md_util_element_symbol(new_elem);
                    ImGui::Text("New Element: %.*s (%.*s)", (int)new_name.len, new_name.ptr, (int)new_sym.len, new_sym.ptr);
                    if (!is_valid) ImGui::PushDisabled();
                    if (ImGui::Button("Apply") && is_valid) {
                        add_atom_elem_mapping(data, atom_name, new_elem);
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

                md_bitfield_t mask = {0};
                md_bitfield_init(&mask, frame_allocator);
                bool apply = false;

                apply |= ImGui::MenuItem("on Atom");
                if (ImGui::IsItemHovered()) {
                    md_bitfield_set_bit(&mask, idx);
                }

                if (data->mold.mol.residue.count > 0 && data->mold.mol.atom.residue_idx && data->mold.mol.atom.residue_idx[idx] != -1) {
                    apply |= ImGui::MenuItem("on Residue");
                    if (ImGui::IsItemHovered()) {
                        const auto res_idx = data->mold.mol.atom.residue_idx[idx];
                        const auto range = data->mold.mol.residue.atom_range[res_idx];
                        md_bitfield_set_range(&mask, range.beg, range.end);
                    }
                }

                if (data->mold.mol.chain.count > 0 && data->mold.mol.atom.chain_idx && data->mold.mol.atom.chain_idx[idx] != -1) {
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
                        load::traj::set_recenter_target(data->mold.traj, &mask);
                        load::traj::clear_cache(data->mold.traj);
                        launch_prefetch_job(data);
                        interpolate_atomic_properties(data);
                        data->mold.dirty_buffers |= MolBit_DirtyPosition;
                        update_md_buffers(data);
                        md_gl_molecule_zero_velocity(&data->mold.gl_mol); // Do this explicitly to update the previous position to avoid motion blur trails
                        ImGui::CloseCurrentPopup();
                    }
                }
                ImGui::EndMenu();
            }
        }
        if (ImGui::BeginMenu("Selection")) {
            if (ImGui::MenuItem("Invert")) {
                md_bitfield_not_inplace(&data->selection.current_selection_mask, 0, data->mold.mol.atom.count);
                data->mold.dirty_buffers |= MolBit_DirtyFlags;
                ImGui::CloseCurrentPopup();
            }
            if (ImGui::IsItemHovered()) {
                md_bitfield_not(&data->selection.current_highlight_mask, &data->selection.current_selection_mask, 0, data->mold.mol.atom.count);
                data->mold.dirty_buffers |= MolBit_DirtyFlags;
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
                grow_mask_by_covalent_bond(&data->selection.grow.mask, data->mold.mol.covalent.bond, data->mold.mol.covalent.count, (int64_t)data->selection.grow.extent);
                break;
            case SelectionGrowth::Radial: {
                const auto& mol = data->mold.mol;
                vec3_t pbc_ext = mat3_mul_vec3(data->simulation_box.box, vec3_set1(1));
                grow_mask_by_radial_extent(&data->selection.grow.mask, &data->selection.current_selection_mask, mol, pbc_ext, data->selection.grow.extent);
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
        bool apply = ImGui::InputQuery("##query", data->selection.query.buf, sizeof(data->selection.query.buf), data->selection.query.query_ok, data->selection.query.error, ImGuiInputTextFlags_AutoSelectAll | ImGuiInputTextFlags_EnterReturnsTrue);
        ImGui::PopItemWidth();

        data->selection.query.query_invalid |= ImGui::IsItemEdited();
        bool preview = ImGui::IsItemActive() || ImGui::IsItemHovered();

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

static void draw_animation_window(ApplicationData* data) {
    ASSERT(data);
    int num_frames = (int)md_trajectory_num_frames(data->mold.traj);
    if (num_frames == 0) return;

    ASSERT(data->timeline.x_values);
    ASSERT(md_array_size(data->timeline.x_values) == num_frames);

    ImGui::SetNextWindowSize({300,200}, ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Animation")) {
        const float item_width = MAX(ImGui::GetContentRegionAvail().x - 80.f, 100.f);
        ImGui::PushItemWidth(item_width);

        ImGui::Text("Num Frames: %i", num_frames);

        md_unit_t time_unit = md_trajectory_time_unit(data->mold.traj);
        double t   = frame_to_time(data->animation.frame, *data);
        double min = data->timeline.x_values[0];
        double max = data->timeline.x_values[num_frames - 1];
        char time_label[64];
        if (unit_empty(time_unit)) {
            snprintf(time_label, sizeof(time_label), "Time");
        } else {
            char unit_buf[32];
            unit_print(unit_buf, sizeof(unit_buf), time_unit);
            snprintf(time_label, sizeof(time_label), "Time (%s)", unit_buf);
        }
        if (ImGui::Combo("Interp.", (int*)(&data->animation.interpolation), "Nearest\0Linear\0Cubic Spline\0\0")) {
            interpolate_atomic_properties(data);
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Interpolation Method for Atom Positions");
        }
        if (ImGui::SliderScalar(time_label, ImGuiDataType_Double, &t, &min, &max, "%.2f")) {
            data->animation.frame = time_to_frame(t, *data);
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
        ImGui::PopItemWidth();
    }
    ImGui::End();
}

static void draw_representations_window(ApplicationData* data) {

    ImGui::SetNextWindowSize({300,200}, ImGuiCond_FirstUseEver);
    ImGui::Begin("Representations", &data->representation.show_window, ImGuiWindowFlags_NoFocusOnAppearing);
    if (ImGui::Button("create new")) {
        create_representation(data);
    }
    ImGui::SameLine();
    if (ImGui::DeleteButton("remove all")) {
        clear_representations(data);
    }
    ImGui::Spacing();
    ImGui::Separator();
    for (int i = 0; i < (int)md_array_size(data->representation.reps); i++) {
        bool update_rep = false;
        auto& rep = data->representation.reps[i];
        const float item_width = MAX(ImGui::GetContentRegionAvail().x - 90.f, 100.f);
        char label[128];
        snprintf(label, sizeof(label), "%s###ID", rep.name);

        ImGui::PushID(i);
        
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
        const char* eye_icon = rep.enabled ? (const char*)ICON_FA_EYE : (const char*)ICON_FA_EYE_SLASH;
        
        const ImVec2 btn_size = {size, size};
        if (ImGui::Button(eye_icon, btn_size)) {
            rep.enabled = !rep.enabled;
            data->representation.atom_visibility_mask_dirty = true;
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Show/Hide");
        }
        ImGui::SameLine(0, spacing);
        if (ImGui::Button((const char*)ICON_FA_COPY, btn_size)) {
            clone_representation(data, rep);
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Duplicate");
        }
        ImGui::SameLine(0, spacing);
        if (ImGui::DeleteButton((const char*)ICON_FA_XMARK, btn_size)) {
            remove_representation(data, i);
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Remove");
        }

        if (draw_content) {
            ImGui::PushItemWidth(item_width);
            if (ImGui::InputQuery("filter", rep.filt, sizeof(rep.filt), rep.filt_is_valid, rep.filt_error)) {
                rep.filt_is_dirty = true;
                update_rep = true;
            }
            if (!rep.type_is_valid) ImGui::PushInvalid();
            if (ImGui::Combo("type", (int*)(&rep.type), "VDW\0Licorice\0Ribbons\0Cartoon\0")) {
                update_rep = true;
            }
            if (!rep.type_is_valid) ImGui::PopInvalid();

            if (ImGui::Combo("color", (int*)(&rep.color_mapping),
                             "Uniform Color\0CPK\0Atom Label\0Atom Idx\0Res Id\0Res Idx\0Chain Id\0Chain Idx\0Secondary Structure\0Property\0")) {
                update_rep = true;
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
                const md_script_property_t* props[32] = {0};
                int num_props = 0;
                for (int64_t j = 0; j < md_array_size(data->display_properties); ++j) {
                    if (data->display_properties[j].full_prop->flags & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) {
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
            if (rep.filt_is_dynamic || rep.color_mapping == ColorMapping::Property) {
                ImGui::Checkbox("auto-update", &rep.dynamic_evaluation);
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
            if (rep.type == RepresentationType::SpaceFill || rep.type == RepresentationType::Licorice) {
                update_rep |= ImGui::SliderFloat("scale", &rep.radius, 0.1f, 4.f);
            }
            if (rep.type == RepresentationType::Ribbons) {
                update_rep |= ImGui::SliderFloat("spline tension", &rep.tension, 0.f, 1.f);
                update_rep |= ImGui::SliderFloat("spline width", &rep.width, 0.1f, 2.f);
                update_rep |= ImGui::SliderFloat("spline thickness", &rep.thickness, 0.1f, 2.f);
            }
            ImGui::PopItemWidth();
            ImGui::Spacing();
            ImGui::Separator();

            ImGui::TreePop();
        }

        ImGui::PopID();

        if (update_rep) {
            update_representation(data, &rep);
        }
    }

    ImGui::End();
}

static void draw_atom_info_window(const ApplicationData& data, int atom_idx) {
    const auto& mol = data.mold.mol;

    // @TODO: Assert things and make this failproof
    if (atom_idx < 0 || atom_idx >= mol.atom.count) return;

    int local_idx = atom_idx;
    const vec3_t pos = { mol.atom.x[atom_idx], mol.atom.y[atom_idx], mol.atom.z[atom_idx] };
    str_t label = mol.atom.name ? mol.atom.name[atom_idx] : str_t{};
    str_t elem = mol.atom.element ? md_util_element_name(mol.atom.element[atom_idx]) : str_t{};
    str_t symbol = mol.atom.element ? md_util_element_symbol(mol.atom.element[atom_idx]) : str_t{};
    int valence = mol.atom.valence ? mol.atom.valence[atom_idx] : 0;

    int res_idx = -1;
    str_t res_name = {};
    int res_id = 0;
    if (mol.residue.count > 0 && mol.atom.residue_idx) {
        res_idx = mol.atom.residue_idx[atom_idx];
        res_name = mol.residue.name[res_idx];
        res_id = mol.residue.id[res_idx];
        local_idx = atom_idx - mol.residue.atom_range[res_idx].beg;
    }

    int chain_idx = -1;
    str_t chain_id = {};
    if (mol.chain.count > 0 && mol.atom.chain_idx) {
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

    char buf[256];
    int len = 0;
    len += snprintf(buf, sizeof(buf), "atom[%i][%i]: %.*s %.*s %.*s (%.2f, %.2f, %.2f)\n", atom_idx, local_idx, (int)label.len, label.ptr, (int)elem.len, elem.ptr, (int)symbol.len, symbol.ptr, pos.x, pos.y, pos.z);
    if (mol.atom.valence) {
        len += snprintf(buf + len, sizeof(buf) - len, "covalent-valence: %i\n", valence);
    }
    if (res_idx) {
        len += snprintf(buf + len, sizeof(buf) - len, "res[%i]: %.*s %i\n", res_idx, (int)res_name.len, res_name.ptr, res_id);
    }
    if (chain_idx) {
        len += snprintf(buf + len, sizeof(buf) - len, "chain[%i]: %.*s\n", chain_idx, (int)chain_id.len, chain_id.ptr);
    }

    /*
    // @TODO: REIMPLEMENT THIS
    if (res_idx < mol.backbone.segment.angleangles.size() && res_idx < mol.backbone.segments.size() && valid_backbone_atoms(mol.backbone.segments[res_idx])) {
        const auto angles = RAD_TO_DEG((vec2)mol.backbone.angles[res_idx]);
        len += snprintf(buff + len, 256 - len, u8"\u03C6: %.1f\u00b0, \u03C8: %.1f\u00b0\n", angles.x, angles.y);
    }
    */

    const ImVec2 offset = { 10.f, 18.f };
    ImGui::SetNextWindowPos(ImGui::GetMousePos() + offset);
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
    ImGui::Begin("##Atom Info", 0,
                 ImGuiWindowFlags_Tooltip | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoDocking);
    ImGui::Text("%s", buf);
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

static void draw_async_task_window(ApplicationData* data) {
    constexpr float WIDTH = 300.f;
    constexpr float MARGIN = 10.f;

    task_system::ID* tasks = task_system::pool_running_tasks(frame_allocator);
    uint32_t num_tasks = (uint32_t)md_array_size(tasks);
    bool any_task_label_visible = false;
    for (uint32_t i = 0; i < num_tasks; i++) {
        str_t label = task_system::task_label(tasks[i]);
        if (!label || label[0] == '\0' || (label[0] == '#' && label[1] == '#')) continue;
        any_task_label_visible = true;
    }
    
    if (any_task_label_visible) {
        ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(viewport->Pos + ImVec2(data->ctx.window.width - WIDTH - MARGIN,
                                                       ImGui::GetCurrentContext()->FontBaseSize + ImGui::GetStyle().FramePadding.y * 2.f + MARGIN));
        ImGui::SetNextWindowSize(ImVec2(WIDTH, 0));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
        ImGui::Begin("##Async Info", 0,
                     ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar |
                     ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing);

        const float pad = 3.0f;
        const float size = ImGui::GetFontSize() + pad * 2;

        char buf[64];
        for (uint32_t i = 0; i < MIN(num_tasks, 8); i++) {
            const auto id = tasks[i];
            str_t label = task_system::task_label(id);
            float fract = task_system::task_fraction_complete(id);

            /*
            if (id == data->tasks.evaluate_filt) {
                uint32_t completed = md_script_eval_num_frames_completed(data->mold.script.filt_eval);
                uint32_t total     = md_script_eval_num_frames_total(data->mold.script.filt_eval);
                if (total > 0) {
                    fract = (float)completed / (float)total;
                }
            }else if (id == data->tasks.evaluate_full) {
                uint32_t completed = md_script_eval_num_frames_completed(data->mold.script.full_eval);
                uint32_t total     = md_script_eval_num_frames_total(data->mold.script.full_eval);
                if (total > 0) {
                    fract = (float)completed / (float)total;
                }
            }
            */

            if (!label || label[0] == '\0' || (label[0] == '#' && label[1] == '#')) continue;

            snprintf(buf, sizeof(buf), "%.*s %.1f%%", (int)label.len, label.ptr, fract * 100.f);
            ImGui::ProgressBar(fract, ImVec2(ImGui::GetContentRegionAvail().x - (size + pad),0), buf);
            ImGui::SameLine();
            if (ImGui::DeleteButton((const char*)ICON_FA_XMARK, ImVec2(size, size))) {
                task_system::task_interrupt(id);
                if (id == data->tasks.evaluate_full) {
                    md_script_eval_interrupt(data->mold.script.full_eval);
                }
                else if(id == data->tasks.evaluate_filt) {
                    md_script_eval_interrupt(data->mold.script.filt_eval);
                }
            }
        }

        ImGui::End();
        ImGui::PopStyleColor();
    }
}


bool draw_property_menu_widgets(DisplayProperty* props, int64_t num_props, uint32_t prop_mask, bool multi_selection = true) {
    int64_t num_candidates = 0;
    int64_t cur_selected_index = -1;
    int64_t new_selected_index = -1;

    for (int64_t i = 0; i < num_props; ++i) {
        DisplayProperty& prop = props[i];
        if (prop.possible_display_mask & prop_mask) {
            num_candidates += 1;
            ImPlot::ItemIcon(prop.col); ImGui::SameLine();
            if (ImGui::Selectable(prop.label, prop.current_display_mask & prop_mask, 0, ImVec2(50,0))) {
                prop.current_display_mask ^= prop_mask;     // Toggle bit
                new_selected_index = i;
                cur_selected_index = i;
            }
            if (cur_selected_index == -1 && prop.current_display_mask & prop_mask) {
                cur_selected_index = i;
            }
        }
    }

    // Clear all other properties if we don't have multiple selection
    if (cur_selected_index != -1 && !multi_selection) {
        for (int64_t i = 0; i < num_props; ++i) {
            if (i == cur_selected_index) continue;
            props[i].current_display_mask &= ~prop_mask;    // Clear bit
        }
    }
    
    if (num_candidates == 0) {
        ImGui::Text("No properties to show.");
    }

    return new_selected_index != -1;
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

bool draw_property_timeline(const ApplicationData& data, const TimelineArgs& args) {
    const ImPlotAxisFlags axis_flags = 0;
    const ImPlotAxisFlags axis_flags_x = axis_flags;
    const ImPlotAxisFlags axis_flags_y = axis_flags | ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_RangeFit | ImPlotAxisFlags_NoLabel |ImPlotAxisFlags_NoTickLabels;
    
    //const ImPlotFlags flags = ImPlotFlags_AntiAliased;
    ImPlotFlags flags = 0;

    if (ImPlot::BeginPlot("##Timeline", ImVec2(-1,args.plot_height), flags)) {
        ImPlot::SetupAxisLinks(ImAxis_X1, args.view_range.beg, args.view_range.end);
        ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, args.view_range.min, args.view_range.max);
        ImPlot::SetupAxes(0, 0, axis_flags_x, axis_flags_y);
        ImPlot::SetupFinish();

        bool active = ImGui::IsItemActive();
        bool print_timeline_tooltip = false;

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
            if (args.values.y_var) {
                ASSERT(args.values.y_min);
                ASSERT(args.values.y_max);

                char lbl[32] = "";
                snprintf(lbl, sizeof(lbl), "%s (mean)", args.lbl);
                ImPlot::PlotLine(lbl, args.values.x, args.values.y, args.values.count);

                ImPlot::SetNextFillStyle(line_col, 0.4f);
                snprintf(lbl, sizeof(lbl), "%s (var)", args.lbl);
                ImPlot::PlotShadedG(lbl,
                    [](int idx, void* payload) -> ImPlotPoint {
                        TimelineArgs* args = (TimelineArgs*)payload;
                        return ImPlotPoint(args->values.x[idx], args->values.y[idx] - args->values.y_var[idx]);
                    },
                    (void*)&args,
                    [](int idx, void* payload) -> ImPlotPoint {
                        TimelineArgs* args = (TimelineArgs*)payload;
                        return ImPlotPoint(args->values.x[idx], args->values.y[idx] + args->values.y_var[idx]);
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
            int32_t frame_idx = CLAMP((int)(time_to_frame(time, data) + 0.5), 0, md_array_size(data.timeline.x_values)-1);
            len += snprintf(buf + len, MAX(0, (int)sizeof(buf) - len), "time: %.2f", time);

            md_unit_t time_unit = md_trajectory_time_unit(data.mold.traj);
            if (!unit_empty(time_unit)) {
                char unit_buf[32];
                unit_print(unit_buf, sizeof(unit_buf), time_unit);
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
                        ImGui::SliderScalar("Extent (frames)", ImGuiDataType_Double, &data->timeline.filter.temporal_window.extent_in_frames, &extent_min, &extent_max, "%1.0f");
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
            map.Pan = ImGuiMouseButton_Right;
            map.Select = ImGuiMouseButton_Right;
            map.SelectCancel = ImGuiMouseButton_Left;
            map.SelectMod = ImGuiMod_Shift;
            map.Menu = -1;
            map.Fit = ImGuiMouseButton_Right;

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
                    .y_var = NULL,
                    .y_min = NULL,
                    .y_max = NULL,
                },
                .view_range = {
                    .beg = &view_beg,
                    .end = &view_end,
                    .min = min_x_value,
                    .max = max_x_value,
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
                const float* y_var = NULL;
                const float* y_min = NULL;
                const float* y_max = NULL;

                int num_y_values = 0;
                    
                if (prop.full_prop->data.aggregate) {
                    y_values = prop.full_prop->data.aggregate->population_mean;
                    y_var = prop.full_prop->data.aggregate->population_var;
                    y_min = prop.full_prop->data.aggregate->population_min;
                    y_max = prop.full_prop->data.aggregate->population_max;
                    num_y_values = prop.full_prop->data.aggregate->num_values;
                } else {
                    y_values = prop.full_prop->data.values;
                    num_y_values = prop.full_prop->data.num_values;
                }
                    
                ASSERT(num_y_values == num_x_values);

                ImGui::PushID(i);
                args.lbl = prop.label;
                args.col = prop.col;
                args.values = {
                    .count = num_y_values,
                    .x = x_values,
                    .y = y_values,
                    .y_var = y_var,
                    .y_min = y_min,
                    .y_max = y_max,
                    .min_y = prop.full_prop->data.min_value,
                    .max_y = prop.full_prop->data.max_value,
                    .unit = {prop.unit_str, (int64_t)strnlen(prop.unit_str, sizeof(prop.unit_str))},
                };
                args.value_filter = {
                    .enabled = data->distributions.filter.enabled,
                    .min = prop.value_filter.beg,
                    .max = prop.value_filter.end,
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
        const int MAX_NUM_BINS = 1024;
        static int num_bins = 128;

        static constexpr int MIN_PLOT_HEIGHT = 10;
        static constexpr int MAX_PLOT_HEIGHT = 1000;
        static int plot_height = 180;

        static bool sync_axis = true;

        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("Properties")) {
                draw_property_menu_widgets(data->display_properties, md_array_size(data->display_properties), DisplayProperty::ShowIn_Distribution);
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Filter")) {
                ImGui::Checkbox("Enabled", &data->distributions.filter.enabled);
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Settings")) {
                ImGui::SliderInt("Plot Height", &plot_height, MIN_PLOT_HEIGHT, MAX_PLOT_HEIGHT);
                if (ImGui::SliderInt("Histogram Bins", &num_bins, MIN_NUM_BINS, MAX_NUM_BINS, "%d", ImGuiSliderFlags_Logarithmic)) {
                    const int up   = next_power_of_two32(num_bins);
                    const int down = up / 2;
                    num_bins = abs(num_bins - down) < abs(num_bins - up) ? down : up;
                }
                ImGui::Checkbox("Sync X-Axis per type", &sync_axis);
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
        //ImPlotFlags flags = ImPlotFlags_AntiAliased;
        ImPlotFlags flags = 0;

        // The distribution properties are always computed as histograms with a resolution of 1024
        // If we have a different number of bins for our visualization, we need to recompute the bins
        float* full_bins = (float*)md_alloc(frame_allocator, num_bins * sizeof(float));
        float* filt_bins = (float*)md_alloc(frame_allocator, num_bins * sizeof(float));

        const int num_props = (int)md_array_size(data->display_properties);
        
        md_array(uint64_t) prop_hash = 0;
        md_array(int)      prop_remap = 0;
        
        if (sync_axis) {
            md_array_resize(prop_hash,  num_props, frame_allocator);
            md_array_resize(prop_remap, num_props, frame_allocator);
            
            for (int i = 0; i < num_props; ++i) {
                DisplayProperty& prop = data->display_properties[i];
                if ((prop.current_display_mask & DisplayProperty::ShowIn_Distribution) == 0) continue;
                const uint64_t min_hash  = std::hash<double>()(prop.full_prop->data.min_range[0]);
                const uint64_t max_hash  = std::hash<double>()(prop.full_prop->data.max_range[0]);
                const uint64_t unit_hash = std::hash<std::string_view>()(prop.unit_str);
                const uint64_t hash = min_hash ^ (max_hash << 1) ^ (unit_hash << 2);
                prop_hash[i] = hash;
                prop_remap[i] = i;
                
                for (int j = 0; j < i; ++j) {
                    if (prop_hash[j] == hash) {
                        prop_remap[i] = j;
                        break;
                    }
                }
            }
        }

        for (int i = 0; i < num_props; ++i) {
            DisplayProperty& prop = data->display_properties[i];
            if ((prop.current_display_mask & DisplayProperty::ShowIn_Distribution) == 0) continue;

            const md_script_property_t* full_prop = prop.full_prop;
            const md_script_property_t* filt_prop = prop.filt_prop;

            double min_x = full_prop->data.min_range[0];
            double max_x = full_prop->data.max_range[0];
            
            double* view_beg = sync_axis ? &data->display_properties[prop_remap[i]].view_range.beg : &prop.view_range.beg;
            double* view_end = sync_axis ? &data->display_properties[prop_remap[i]].view_range.end : &prop.view_range.end;

            // We need to avoid the axis collapsing, otherwise Implot will get stuck in an infinite loop
            if (min_x == max_x) {
                max_x = min_x + 0.001;
            }

            const float* full_src = 0;
            const float* filt_src = 0;
            int num_values_src = 0;

            if (full_prop->flags & MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) {
                full_src = full_prop->data.values;
                filt_src = filt_prop->data.values;
                num_values_src = full_prop->data.num_values;
            } else if (full_prop->flags & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) {
                full_src = prop.full_hist.bin;
                filt_src = prop.filt_hist.bin;
                num_values_src = (int)ARRAY_SIZE(prop.full_hist.bin);
            }

            // Downsample bins
            downsample_histogram(full_bins, num_bins, full_src, num_values_src);
            downsample_histogram(filt_bins, num_bins, filt_src, num_values_src);

            double max_y = 0;
            for (int64_t j = 0; j < num_bins; ++j) {
                max_y = MAX(max_y, full_bins[j]);
            }
            max_y = MAX(max_y, 0.001);
            
            ImGui::PushID(i);
            const double bar_width = (max_x - min_x) / (num_bins-1);
            const double bar_off = min_x;
            const double bar_scl = (max_x - min_x) / (num_bins-1);

            char label[128];
            int len = snprintf(label, sizeof(label), "%s ", prop.label);
            if (!unit_empty(prop.unit)) {
                char unit_buf[64];
                int unit_len = unit_print(unit_buf, sizeof(unit_buf), prop.unit);
                snprintf(label + len, MAX(0, (int)sizeof(label) - len), "(%*s)", unit_len, unit_buf);
            }

            if (ImPlot::BeginPlot(label, ImVec2(-1,plot_height), flags)) {
                ImPlot::SetupAxisLimits(ImAxis_Y1, 0, max_y * 1.1, ImGuiCond_Always);
                ImPlot::SetupAxisLimits(ImAxis_X1, min_x, max_x, ImGuiCond_Once);
                ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, min_x, max_x);
                ImPlot::SetupAxisLinks(ImAxis_X1, view_beg, view_end);
                ImPlot::SetupAxes(0, 0, axis_flags_x, axis_flags_y);
                ImPlot::SetupFinish();

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

                auto getter = [](int idx, void* data) -> ImPlotPoint {
                    PlotData* pd = (PlotData*)data;
                    return ImPlotPoint(idx * pd->scale + pd->offset, (double)pd->bars[idx]);
                };

                plot_data.bars = full_bins;
                ImPlot::SetNextFillStyle(ImVec4(0,0,0,-1), 1.0f);
                ImPlot::SetNextLineStyle(ImVec4(0,0,0,0), 0);
                ImPlot::PlotBarsG("full", getter, &plot_data, num_bins, bar_width);

                if (data->timeline.filter.enabled) {
                    plot_data.bars = filt_bins;
                    ImPlot::SetNextFillStyle(ImVec4(1,1,0,1), 0.3f);
                    ImPlot::SetNextLineStyle(ImVec4(0,0,0,0), 0);
                    ImPlot::PlotBarsG("filt", getter, &plot_data, num_bins, bar_width);
                }

                if (data->distributions.filter.enabled) {
                    if ((is_active || is_hovered) && ImGui::GetIO().KeyMods == ImGuiMod_Shift) {
                        if (ImGui::IsMouseClicked(0)) {
                            prop.value_filter.beg = ImPlot::GetPlotMousePos().x;
                            prop.value_filter.end = ImPlot::GetPlotMousePos().x;
                        } else if (ImGui::IsMouseDragging(0)) {
                            prop.value_filter.end = ImPlot::GetPlotMousePos().x;
                            prop.value_filter.beg = MIN(prop.value_filter.beg, prop.value_filter.end);
                            prop.value_filter.end = MAX(prop.value_filter.beg, prop.value_filter.end);
                        }
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

    ImGui::SetNextWindowSize({300,350}, ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Shape Space", &data->shape_space.show_window, ImGuiWindowFlags_MenuBar)) {
        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("Settings")) {
                static constexpr float marker_min_size = 0.01f;
                static constexpr float marker_max_size = 10.0f;
                ImGui::SliderFloat("Marker Size", &data->shape_space.marker_size, marker_min_size, marker_max_size);
                ImGui::EndMenu();
            }
            ImGui::EndMenuBar();
        }

        const float item_width = MAX(ImGui::GetContentRegionAvail().x, 150.f);
        ImGui::PushItemWidth(item_width);
        data->shape_space.evaluate |= ImGui::InputQuery("##input", data->shape_space.input, sizeof(data->shape_space.input), data->shape_space.input_valid);
        ImGui::PopItemWidth();
        if (ImGui::IsItemHovered()) {
            if (data->shape_space.input_valid) {
                md_bitfield_clear(&data->selection.current_highlight_mask);
                for (int64_t i = 0; i < md_array_size(data->shape_space.bitfields); ++i) {
                    md_bitfield_or_inplace(&data->selection.current_highlight_mask, &data->shape_space.bitfields[i]);
                }
            } else if (data->shape_space.error[0] != '\0') {
                ImGui::SetTooltip("%s", data->shape_space.error);
            }
        }

        ImPlotFlags flags = ImPlotFlags_Equal | ImPlotFlags_NoMenus; // ImPlotFlags_AntiAliased;
        ImPlotAxisFlags axis_flags = ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoTickLabels | ImPlotAxisFlags_NoTickMarks;

        const float x_reset[2] = {-0.10f, 1.10f};
        const float y_reset[2] = {-0.10f, 0.98f};

        if (ImPlot::BeginPlot("##Shape Space Plot", ImVec2(-1,-1), flags)) {
            ImPlot::SetupAxesLimits(x_reset[0], x_reset[1], y_reset[0], y_reset[1], ImGuiCond_Appearing);
            ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, x_reset[0], x_reset[1]);
            ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, y_reset[0], y_reset[1]);
            ImPlot::SetupAxes(0, 0, axis_flags, axis_flags);
            ImPlot::SetupFinish();

			const ImPlotPoint lin(0.0, 0.0);
			const ImPlotPoint pla(1.0, 0.0);
			const ImPlotPoint iso(0.5, 0.86602540378);

            const ImVec2 p0 = ImPlot::PlotToPixels(lin);
            const ImVec2 p1 = ImPlot::PlotToPixels(pla);
            const ImVec2 p2 = ImPlot::PlotToPixels(iso);
            
            const ImVec2 pos_iso = ImPlot::PlotToPixels(iso);
            const ImVec2 pos_pla = ImPlot::PlotToPixels(pla);
			const ImVec2 pos_lin = ImPlot::PlotToPixels(lin);

            static const char* text_iso = "Isotropic";
            static const char* text_pla = "Planar";
            static const char* text_lin = "Linear";
            
            const ImVec2 text_iso_offset = pos_iso + ImGui::CalcTextSize(text_iso) * ImVec2(-0.5f, -1.0f);
            const ImVec2 text_pla_offset = pos_pla + ImGui::CalcTextSize(text_pla) * ImVec2(-0.5f,  0.0f);
            const ImVec2 text_lin_offset = pos_lin + ImGui::CalcTextSize(text_lin) * ImVec2(-0.5f,  0.0f);

            ImPlot::PushPlotClipRect();
            ImPlot::GetPlotDrawList()->AddTriangleFilled(p0, p1, p2, IM_COL32(255,255,255,20));
            ImPlot::GetPlotDrawList()->AddTriangle(p0, p1, p2, IM_COL32(255,255,255,50));
            
            ImPlot::GetPlotDrawList()->AddText(text_iso_offset, IM_COL32(255, 255, 255, 255), text_iso);
			ImPlot::GetPlotDrawList()->AddText(text_pla_offset, IM_COL32(255, 255, 255, 255), text_pla);
			ImPlot::GetPlotDrawList()->AddText(text_lin_offset, IM_COL32(255, 255, 255, 255), text_lin);
            ImPlot::PopPlotClipRect();

            ImPlot::PushStyleVar(ImPlotStyleVar_Marker, ImPlotMarker_Square);
            ImPlot::PushStyleColor(ImPlotCol_MarkerFill, ImVec4(0,0,0,0));
            ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, ImVec4(0,0,0,0));
            ImPlot::PlotScatter("##Hidden reset helper", x_reset, y_reset, 2);
            ImPlot::PopStyleColor(2);
            ImPlot::PopStyleVar();

            auto getter = [](int idx, void* user_data) -> ImPlotPoint {
                const vec2_t* coords = (vec2_t*)user_data;
                return {coords[idx].x, coords[idx].y};
            };

            int32_t hovered_structure_idx = -1;
            int32_t hovered_point_idx = -1;
            float   hovered_point_min_dist = FLT_MAX;

            vec2_t mouse_coord = {(float)ImPlot::GetPlotMousePos().x, (float)ImPlot::GetPlotMousePos().y};

            ImPlotRect lim = ImPlot::GetPlotLimits();
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
                md_bitfield_copy(&data->selection.current_highlight_mask, &data->shape_space.bitfields[hovered_structure_idx]);
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
                    md_bitfield_copy(&data->selection.current_highlight_mask, &data->shape_space.bitfields[structure_idx]);
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
        const int64_t num_frames = md_trajectory_num_frames(data->mold.traj);
        if (num_frames > 0) {
            if (task_system::task_is_running(data->tasks.shape_space_evaluate)) {
                task_system::task_interrupt(data->tasks.shape_space_evaluate);
            }
            else if (md_semaphore_try_aquire(&data->mold.script.ir_semaphore)) {
                defer { md_semaphore_release(&data->mold.script.ir_semaphore); };
                
                data->shape_space.evaluate = false;
                md_array_shrink(data->shape_space.coords, 0);
                md_array_shrink(data->shape_space.weights, 0);
                
                for (int64_t i = 0; i < md_array_size(data->shape_space.bitfields); ++i) {
                    md_bitfield_free(&data->shape_space.bitfields[i]);
                }
                md_array_shrink(data->shape_space.bitfields, 0);
                
                if (md_filter_evaluate(&data->shape_space.bitfields, str_from_cstr(data->shape_space.input), &data->mold.mol, data->mold.script.ir, NULL, data->shape_space.error, sizeof(data->shape_space.error), persistent_allocator)) {
                    data->shape_space.input_valid = true;
                    data->shape_space.num_structures = (int32_t)md_array_size(data->shape_space.bitfields);
                    
                    if (data->shape_space.num_structures > 0) {
                        data->shape_space.num_frames = (int32_t)num_frames;
                        md_array_resize(data->shape_space.coords,  num_frames * data->shape_space.num_structures, persistent_allocator);
                        MEMSET(data->shape_space.coords, 0, md_array_bytes(data->shape_space.coords));
                        md_array_resize(data->shape_space.weights, num_frames * data->shape_space.num_structures, persistent_allocator);
                        MEMSET(data->shape_space.weights, 0, md_array_bytes(data->shape_space.weights));

                        data->tasks.shape_space_evaluate = task_system::pool_enqueue(STR("Eval Shape Space"), 0, (uint32_t)num_frames, [](uint32_t range_beg, uint32_t range_end, void* user_data) {
                            ApplicationData* data = (ApplicationData*)user_data;
                            int64_t stride = ALIGN_TO(data->mold.mol.atom.count, md_simd_f32_width);
                            float* coords = (float*)md_alloc(default_allocator, stride * 3 * sizeof(float));
                            float* x = coords + stride * 0;
                            float* y = coords + stride * 1;
                            float* z = coords + stride * 2;

                            const vec2_t p[3] = {{0.0f, 0.0f}, {1.0f, 0.0f}, {0.5f, 0.86602540378f}};

                            vec3_t* xyz = 0;
                            for (uint32_t frame_idx = range_beg; frame_idx < range_end; ++frame_idx) {
                                md_trajectory_load_frame(data->mold.traj, frame_idx, NULL, x, y, z);
                                for (int64_t i = 0; i < md_array_size(data->shape_space.bitfields); ++i) {
                                    md_array_ensure(xyz, (int64_t)md_bitfield_popcount(&data->shape_space.bitfields[i]), default_allocator);
                                    int64_t beg_bit = data->shape_space.bitfields[i].beg_bit;
                                    int64_t end_bit = data->shape_space.bitfields[i].end_bit;
                                    int64_t count = 0;
                                    vec3_t com = {0,0,0};
                                    while ((beg_bit = md_bitfield_scan(&data->shape_space.bitfields[i], beg_bit, end_bit)) != 0) {
                                        int64_t src_idx = beg_bit - 1;
                                        vec3_t pos = {x[src_idx], y[src_idx], z[src_idx]};
                                        com = com + pos;
                                        xyz[count++] = pos;
                                    }
                                    com = com / (float)count;

                                    const mat3_eigen_t eigen = mat3_eigen(mat3_covariance_matrix_vec3(xyz, com, count));
                                    const float scl = 1.0f / (eigen.values[0] + eigen.values[1] + eigen.values[2]);
                                    const int64_t dst_idx = data->shape_space.num_frames * i + frame_idx;
                                    const vec3_t w = {(eigen.values[0] - eigen.values[1]) * scl, 2.0f * (eigen.values[1] - eigen.values[2]) * scl, 3.0f * eigen.values[2] * scl};

                                    data->shape_space.weights[dst_idx] = w;
                                    data->shape_space.coords[dst_idx] = p[0] * w[0] + p[1] * w[1] + p[2] * w[2];
                                }
                            }
                            md_array_free(xyz, default_allocator);
                        }, data);
                    } else {
                        snprintf(data->shape_space.error, sizeof(data->shape_space.error), "Expression did not evaluate into any bitfields");
                    }
                }
            }
        } else {
            data->shape_space.evaluate = false;
            snprintf(data->shape_space.error, sizeof(data->shape_space.error), "Missing trajectory for evaluating expression");
        }
    }
}

static void draw_ramachandran_window(ApplicationData* data) {

    ImGui::SetNextWindowSize({300,350}, ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Ramachandran", &data->ramachandran.show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_MenuBar)) {
        enum RamachandranDisplayMode {
            IsoLevels,
            IsoLines,
            Colormap,
        };

        constexpr const char* plot_labels[4] = {"##General", "##Glycine", "##Proline", "##Pre-Proline"};
        constexpr const char* layer_labels[3] = { "Reference", "Full Trajectory", "Filtered Trajectory" };
        constexpr const char* option_labels[3] = { "IsoLevels", "IsoLines", "Colormap" };

        constexpr const float min_ext = -180.0f;
        constexpr const float max_ext = 180.0f;
        constexpr const float reset_coords[2] = { min_ext, max_ext };
        constexpr const char* x_lbl = "\xc2\xb0\xcf\x86";   // utf8 Degree Phi
        constexpr const char* y_lbl = "\xc2\xb0\xcf\x88";   // utf8 Degree Psi

        constexpr const ImPlotFlags flags = ImPlotFlags_Equal | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect; // | ImPlotFlags_AntiAliased;
        constexpr const ImPlotFlags subplotflags = ImPlotSubplotFlags_NoResize | ImPlotSubplotFlags_NoMenus;
        constexpr const ImPlotAxisFlags axis_flags = ImPlotAxisFlags_Foreground | ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoTickLabels;

        static ImPlotRect selection_rect;
        static SelectionOperator op = SelectionOperator::Or;
        static bool is_selecting[4] = {false, false, false, false};
        
        static float ref_alpha  = 0.85f;
        static float full_alpha = 0.85f;
        static float filt_alpha = 0.85f;

        static RamachandranDisplayMode display_mode[3] = { IsoLevels, IsoLines, IsoLines };
        static ImPlotColormap colormap[3] = { ImPlotColormap_Hot, ImPlotColormap_Plasma, ImPlotColormap_Viridis };
        static vec4_t isoline_colors[3]   = { {1,1,1,1}, {1,1,1,1}, {1,1,1,1} };

        static int layout_mode = 0;
        static bool show_layer[4] = {true, true, true, true};

        static ImPlotRect viewrect = ImPlotRect(min_ext, max_ext, min_ext, max_ext);
        
        if (viewrect.Size().x < 1) {
            viewrect.X.Max = viewrect.X.Min + 1;
        }
        if (viewrect.Size().y < 1) {
            viewrect.Y.Max = viewrect.Y.Min + 1;
        }

        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("Layers")) {

                ImGui::Separator();
                ImGui::Text("Layers");
                for (int i = 0; i < 3; ++i) {
                    ImGui::Checkbox(layer_labels[i], &show_layer[i]);
                    if (show_layer[i]) {
                        ImGui::PushID(i);
                        if (ImGui::BeginCombo(layer_labels[i], option_labels[display_mode[i]])) {
                            if (ImGui::Selectable(option_labels[IsoLevels], display_mode[i] == IsoLevels)) display_mode[i] = IsoLevels;
                            if (ImGui::Selectable(option_labels[IsoLines],  display_mode[i] == IsoLines))  display_mode[i] = IsoLines;
                            if (ImGui::Selectable(option_labels[Colormap],  display_mode[i] == Colormap))  display_mode[i] = Colormap;
                            ImGui::EndCombo();
                        }
                        if (display_mode[i] == Colormap) {
                            ImPlot::ColorMapSelection("Color Map", &colormap[i]);
                        } else if (display_mode[i] == IsoLines) {
                            ImGui::ColorEdit4Minimal("Line Color", isoline_colors[i].elem);
                        }
                        ImGui::PopID();
                    }
                    ImGui::Separator();
                }
                ImGui::Checkbox("Current", &show_layer[3]);
                if (show_layer[3]) {
                    ImGui::SliderFloat("Point Size", &data->ramachandran.style.base_radius, 1.0f, 10.0f);
                    ImGui::ColorEdit4Minimal("Point Color", data->ramachandran.style.base_color.elem);
                }
                ImGui::Separator();
                if (ImGui::SliderFloat("Density Blur Sigma", &data->ramachandran.blur_sigma, 0.1f, 10.0f)) {
                    data->ramachandran.full_fingerprint = 0;
                    data->ramachandran.filt_fingerprint = 0;
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("View")) {
                if (ImGui::Selectable("Side-By-Side", layout_mode == 0)) layout_mode = 0;
                if (ImGui::Selectable("General",      layout_mode == 1)) layout_mode = 1;
                if (ImGui::Selectable("Glycine",      layout_mode == 2)) layout_mode = 2;
                if (ImGui::Selectable("Proline",      layout_mode == 3)) layout_mode = 3;
                if (ImGui::Selectable("Preproline",   layout_mode == 4)) layout_mode = 4;
                ImGui::EndMenu();
            }
            ImGui::EndMenuBar();
        }

        const auto& mol = data->mold.mol;
        md_bitfield_t* selection_mask = &data->selection.current_selection_mask;
        md_bitfield_t* highlight_mask = &data->selection.current_highlight_mask;

        const int plot_offset = MAX(0, layout_mode - 1);
        const int plot_cols = (layout_mode == 0) ? 2 : 1;
        const int plot_rows = (layout_mode == 0) ? 2 : 1;
        const int num_plots = plot_cols * plot_rows;

        ImPlotInputMap& map = ImPlot::GetInputMap();
        map.OverrideMod = ImGuiMod_Shift;

        auto formatter = [](double value, char* buff, int size, void* user_data) -> int {
            const char* suffix = (const char*)user_data;
            value = deperiodize(value, 0, 360.0);
            return snprintf(buff, size, "%g%s", value, suffix);
        };

        ImVec2 plot_size = {0,0};

        if (ImPlot::BeginSubplots("##Ramaplots", plot_rows, plot_cols, ImVec2(-1, -1), subplotflags | ImPlotSubplotFlags_ShareItems)) {
            for (int plot_idx = plot_offset; plot_idx < plot_offset + num_plots; ++plot_idx) {
                if (ImPlot::BeginPlot(plot_labels[plot_idx], ImVec2(), flags)) {
                    ImPlotPoint view_mid = { (viewrect.X.Min + viewrect.X.Max) * 0.5, (viewrect.Y.Min + viewrect.Y.Max) * 0.5 };

					ImPlot::SetupAxesLimits(min_ext, max_ext, min_ext, max_ext, ImPlotCond_Once);
                    ImPlot::SetupAxisLinks(ImAxis_X1, &viewrect.X.Min, &viewrect.X.Max);
                    ImPlot::SetupAxisLinks(ImAxis_Y1, &viewrect.Y.Min, &viewrect.Y.Max);
                    // @NOTE(Robin): This wont work out of the box due to the periodic domain.
                    //ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, -720, +720);
                    //ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, -720, +720);
                    ImPlot::SetupAxes(x_lbl, y_lbl, axis_flags, axis_flags);
                    
                    ImPlot::SetupAxisFormat(ImAxis_X1, formatter, (void*)x_lbl);
                    ImPlot::SetupAxisFormat(ImAxis_Y1, formatter, (void*)y_lbl);

                    ImPlot::SetupFinish();

                    viewrect = ImPlot::GetPlotLimits();

                    ImPlot::PushStyleVar(ImPlotStyleVar_Marker, ImPlotMarker_Square);
                    ImPlot::PushStyleColor(ImPlotCol_MarkerFill, ImVec4(0,0,0,0));
                    ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, ImVec4(0,0,0,0));
                    ImPlot::PlotScatter("##Hidden reset helper", reset_coords, reset_coords, 2);
                    ImPlot::PopStyleColor(2);
                    ImPlot::PopStyleVar();

                    //ImPlot::PlotDummy("Reference");
                    //ImPlot::PlotDummy("Full");
                    //ImPlot::PlotDummy("Filtered");
                    //ImPlot::PlotDummy("Current");

                    //ImPlot::GetCurrentContext()->CurrentSubplot->Items.GetLegendItem(3)->Color = 0xFFFFFFFF;

                    //bool show_ref  = ImPlot::GetCurrentContext()->CurrentSubplot->Items.GetLegendItem(0)->Show;
                    //bool show_full = ImPlot::GetCurrentContext()->CurrentSubplot->Items.GetLegendItem(1)->Show;
                    //bool show_filt = ImPlot::GetCurrentContext()->CurrentSubplot->Items.GetLegendItem(2)->Show;
                    //bool show_curr = ImPlot::GetCurrentContext()->CurrentSubplot->Items.GetLegendItem(3)->Show;
                    const bool show_ref  = show_layer[0];
                    const bool show_full = show_layer[1];
                    const bool show_filt = show_layer[2];
                    const bool show_curr = show_layer[3];

                    ImVec2 plot_min = ImPlot::PlotToPixels(viewrect.Min());
                    ImVec2 plot_max = ImPlot::PlotToPixels(viewrect.Max());
                    plot_size = { fabsf(plot_max.x - plot_min.x), fabsf(plot_max.y - plot_min.y) };

                    if (show_ref || show_full || show_filt) {
                        ImPlot::PushPlotClipRect();
                        ImDrawList* dl = ImPlot::GetPlotDrawList();

                        const uint32_t ref_tex  = display_mode[0] == Colormap ? data->ramachandran.data.ref.map_tex[plot_idx]  : data->ramachandran.data.ref.iso_tex[plot_idx];
                        const uint32_t full_tex = display_mode[1] == Colormap ? data->ramachandran.data.full.map_tex[plot_idx] : data->ramachandran.data.full.iso_tex[plot_idx];
                        const uint32_t filt_tex = display_mode[2] == Colormap ? data->ramachandran.data.filt.map_tex[plot_idx] : data->ramachandran.data.filt.iso_tex[plot_idx];

                        if (show_ref)  dl->AddImage((ImTextureID)(intptr_t)ref_tex,  plot_min, plot_max, {0,0}, {1,1}, ImColor(1.0f, 1.0f, 1.0f, ref_alpha));
                        if (show_full) dl->AddImage((ImTextureID)(intptr_t)full_tex, plot_min, plot_max, {0,0}, {1,1}, ImColor(1.0f, 1.0f, 1.0f, full_alpha));
                        if (show_filt) dl->AddImage((ImTextureID)(intptr_t)filt_tex, plot_min, plot_max, {0,0}, {1,1}, ImColor(1.0f, 1.0f, 1.0f, filt_alpha));

                        ImPlot::PopPlotClipRect();
                    }

                    const bool hovered = ImPlot::IsPlotHovered();
                    const bool shift   = ImGui::GetIO().KeyShift;
                    ImPlotPoint mouse_coord = ImPlot::GetPlotMousePos();

                    if (is_selecting[plot_idx]) {
                        data->mold.dirty_buffers |= MolBit_DirtyFlags;
                        selection_rect.X.Max = mouse_coord.x;
                        selection_rect.Y.Max = mouse_coord.y;
                        if (selection_rect.Size().x != 0 && selection_rect.Size().y != 0) {
                            md_bitfield_copy(highlight_mask, selection_mask);
                            ImPlot::DragRect(0, &selection_rect.X.Min, &selection_rect.Y.Min, &selection_rect.X.Max, &selection_rect.Y.Max, ImVec4(0,0,0,0.5f), ImPlotDragToolFlags_NoInputs);
                        }

                        if (!shift) {
                            is_selecting[plot_idx] = false;
                        }
                    } else if (hovered && shift && (ImGui::IsMouseClicked(0) || ImGui::IsMouseClicked(1))) {
                        is_selecting[plot_idx] = true;
                        selection_rect = ImPlotRect(mouse_coord.x, mouse_coord.x, mouse_coord.y, mouse_coord.y);
                        op = ImGui::IsMouseClicked(0) ? SelectionOperator::Or : SelectionOperator::AndNot;
                    }

                    const double cut_d2 = fabs(ImPlot::PixelsToPlot(9,0).x - ImPlot::PixelsToPlot(0,0).x);
                    double min_d2 = DBL_MAX;
                    int64_t mouse_hover_idx = -1;

                    if (show_curr && mol.backbone.angle) {
                        const uint32_t* indices = data->ramachandran.rama_type_indices[plot_idx];

                        double min_x = MIN(selection_rect.X.Min, selection_rect.X.Max);
                        double max_x = MAX(selection_rect.X.Min, selection_rect.X.Max);
                        double min_y = MIN(selection_rect.Y.Min, selection_rect.Y.Max);
                        double max_y = MAX(selection_rect.Y.Min, selection_rect.Y.Max);

                        double ref_x = (min_x + max_x) * 0.5;
                        double ref_y = (min_y + max_y) * 0.5;

                        mouse_coord.x = deperiodize(mouse_coord.x, ref_x, 360.0);
                        mouse_coord.y = deperiodize(mouse_coord.y, ref_y, 360.0);

                        for (int64_t i = 0; i < md_array_size(indices); ++i) {
                            uint32_t idx = indices[i];

                            if (mol.backbone.angle[idx].phi == 0 && mol.backbone.angle[idx].psi == 0) continue;

                            ImPlotPoint coord = ImPlotPoint(RAD_TO_DEG(mol.backbone.angle[idx].phi), RAD_TO_DEG(mol.backbone.angle[idx].psi));
                            coord.x = deperiodize(coord.x, ref_x, 360.0);
                            coord.y = deperiodize(coord.y, ref_y, 360.0);
                            
                            if (is_selecting[plot_idx]) {
                                if (min_x <= coord.x && coord.x <= max_x && min_y <= coord.y && coord.y <= max_y) {
                                    md_residue_idx_t res_idx = mol.backbone.residue_idx[idx];
                                    if (res_idx < mol.residue.count) {
                                        modify_field(highlight_mask, mol.residue.atom_range[res_idx], op);
                                    }
                                }
                            }
                            if (hovered) {
                                double d2 = (coord.x - mouse_coord.x) * (coord.x - mouse_coord.x) + (coord.y - mouse_coord.y) * (coord.y - mouse_coord.y);
                                if (d2 < cut_d2 && d2 < min_d2) {
                                    min_d2 = d2;
                                    mouse_hover_idx = idx;
                                }
                            }
                        }

                        if (is_selecting[plot_idx]) {
                            grow_mask_by_current_selection_granularity(highlight_mask, *data);
                            data->mold.dirty_buffers |= MolBit_DirtyFlags;
                        }

                        if (mouse_hover_idx != -1) {
                            if (mouse_hover_idx < mol.backbone.count) {
                                md_residue_idx_t res_idx = mol.backbone.residue_idx[mouse_hover_idx];
                                if (res_idx < mol.residue.count) {
                                    modify_field(highlight_mask, mol.residue.atom_range[res_idx], SelectionOperator::Or);
                                    grow_mask_by_current_selection_granularity(highlight_mask, *data);
                                    data->mold.dirty_buffers |= MolBit_DirtyFlags;
                                }
                            }
                        }

                        uint32_t* selection_indices = 0;
                        uint32_t* highlight_indices = 0;

                        for (uint32_t i = 0; i < (uint32_t)md_array_size(indices); ++i) {
                            uint32_t idx = indices[i];
                            if (mol.backbone.angle[idx].phi == 0 && mol.backbone.angle[idx].psi == 0) continue;

                            int64_t atom_idx = mol.backbone.atoms[idx].ca;
                            if (md_bitfield_test_bit(highlight_mask, atom_idx)) {
                                md_array_push(highlight_indices, idx, frame_allocator);
                            }

                            if (md_bitfield_test_bit(selection_mask, atom_idx)) {
                                md_array_push(selection_indices, idx, frame_allocator);
                            }
                        }

                        struct UserData {
                            const vec2_t* coords;
                            const uint32_t* indices;
                            ImPlotPoint     view_center; // We use this to deperiodize the coordinates
                        };

                        auto index_getter = [](int i, void* user_data) -> ImPlotPoint {
                            const UserData* data = (UserData*)user_data;
                            uint32_t idx = data->indices[i];

                            if (data->coords[idx].x == 0 && data->coords[idx].y == 0) {
                                return { INFINITY, INFINITY }; // Hide by INF!
                            }

                            float x = deperiodizef(RAD_TO_DEG(data->coords[idx].x), data->view_center.x, 360.0f);
                            float y = deperiodizef(RAD_TO_DEG(data->coords[idx].y), data->view_center.y, 360.0f);
                            return { x, y };
                        };

                        ImPlot::SetNextMarkerStyle(ImPlotMarker_Square, data->ramachandran.style.base_radius, vec_cast(data->ramachandran.style.base_color), 1.2f, vec_cast(data->ramachandran.style.border_color));
                        ImPlot::SetNextLineStyle(ImVec4(1, 1, 1, 1));
                        if (md_array_size(indices) > 0) {
                            UserData user_data = { (const vec2_t*)(mol.backbone.angle), indices, view_mid };
                            ImPlot::PlotScatterG("##Current", index_getter, &user_data, (int)md_array_size(indices));
                        }

                        if (md_array_size(selection_indices) > 0) {
                            UserData user_data = { (const vec2_t*)(mol.backbone.angle), selection_indices, view_mid };
                            ImPlot::SetNextMarkerStyle(ImPlotMarker_Square, data->ramachandran.style.base_radius + 1, vec_cast(data->ramachandran.style.selection_color), 2.0f, vec_cast(data->ramachandran.style.border_color));
                            ImPlot::PlotScatterG("##Selection", index_getter, &user_data, (int)md_array_size(selection_indices));
                        }

                        if (md_array_size(highlight_indices) > 0) {
                            UserData user_data = { (const vec2_t*)(mol.backbone.angle), highlight_indices, view_mid };
                            ImPlot::SetNextMarkerStyle(ImPlotMarker_Square, data->ramachandran.style.base_radius + 1, vec_cast(data->ramachandran.style.highlight_color), 2.0f, vec_cast(data->ramachandran.style.border_color));
                            ImPlot::PlotScatterG("##Highlight", index_getter, &user_data, (int)md_array_size(highlight_indices));
                        }
                    }

                    if (ImGui::IsMouseReleased(0) || ImGui::IsMouseReleased(1)) {
                        if (is_selecting[plot_idx]) {
                            is_selecting[plot_idx] = false;
                            auto size = selection_rect.Size();
                            if (size.x == 0 && size.y == 0) {
                                if (mouse_hover_idx != -1) {
                                    modify_field(selection_mask, highlight_mask, op);
                                } else {
                                    if (ImGui::IsMouseReleased(1)) {
                                        md_bitfield_clear(selection_mask);
                                    }
                                }
                            }
                            else {
                                // Commit whatever is in the highlight mask
                                modify_field(selection_mask, highlight_mask, SelectionOperator::Set);
                            }
                        } 
                    }

                    ImPlot::EndPlot();
                }
            }

            if (viewrect.X.Max < viewrect.X.Min) {
                ImSwap(viewrect.X.Min, viewrect.X.Max);
            }
            if (viewrect.Y.Max < viewrect.Y.Min) {
                ImSwap(viewrect.Y.Min, viewrect.Y.Max);
            }

            const double ratio = viewrect.Size().x / viewrect.Size().y;
            const double mid_x = (viewrect.X.Min + viewrect.X.Max) * 0.5;
            const double mid_y = (viewrect.Y.Min + viewrect.Y.Max) * 0.5;

            if (viewrect.Size().x > 360) {
                viewrect.X.Min = mid_x - 180;
                viewrect.X.Max = mid_x + 180;
                viewrect.Y.Min = mid_y - 0.5 * viewrect.Size().x / ratio;
                viewrect.Y.Max = mid_y + 0.5 * viewrect.Size().x / ratio;
            } else if (viewrect.Size().x < 10) {
                viewrect.X.Min = mid_x - 5;
                viewrect.X.Max = mid_x + 5;
                viewrect.Y.Min = mid_y - 0.5 * viewrect.Size().x / ratio;
                viewrect.Y.Max = mid_y + 0.5 * viewrect.Size().x / ratio;
            }

            if (viewrect.Size().y > 360) {
                viewrect.Y.Min = mid_y - 180;
                viewrect.Y.Max = mid_y + 180;
                viewrect.X.Min = mid_x - 0.5 * viewrect.Size().y * ratio;
                viewrect.X.Max = mid_x + 0.5 * viewrect.Size().y * ratio;
            } else if (viewrect.Size().y < 10) {
                viewrect.Y.Min = mid_y - 5;
                viewrect.Y.Max = mid_y + 5;
                viewrect.X.Min = mid_x - 0.5 * viewrect.Size().y * ratio;
                viewrect.X.Max = mid_x + 0.5 * viewrect.Size().y * ratio;
            }

            ImPlot::EndSubplots();
        }

        vec4_t viewport = { (float)viewrect.X.Min / 180.0f, (float)viewrect.Y.Min / 180.0f, (float)viewrect.X.Max / 180.0f, (float)viewrect.Y.Max / 180.0f };
        viewport = viewport * 0.5f + 0.5f;

        const float* ref_sum = data->ramachandran.data.ref.den_sum;
        const float* full_sum = data->ramachandran.data.full.den_sum;
        const float* filt_sum = data->ramachandran.data.filt.den_sum;

        const float ref_iso_values[4][3] = {
            {0, 0.0005f, 0.02f},  // 99.95%, 98% for General
            {0, 0.0020f, 0.02f},  // 99.80%, 98% for Others
            {0, 0.0020f, 0.02f},
            {0, 0.0020f, 0.02f},
        };

        const uint32_t ref_iso_level_colors[4][3] = {
            {0x00000000, 0xFFFFE8B3, 0xFFFFD97F},
            {0x00000000, 0xFFC5E8FF, 0xFF7FCCFF},
            {0x00000000, 0xFFC5FFD0, 0xFF8CFF7F},
            {0x00000000, 0xFFFFE8B3, 0xFFFFD97F}
        };

        PUSH_GPU_SECTION("RENDER RAMA");

        uint32_t tex_dim = (uint32_t)MAX(plot_size.x, plot_size.y);

        rama_rep_t* reps[3] = {
            &data->ramachandran.data.ref,
            &data->ramachandran.data.full,
            &data->ramachandran.data.filt
        };

        const float scl[3][4] = {
            { 1.0f, 1.0f, 1.0f, 1.0f },
            { full_sum[0] / ref_sum[0], full_sum[1] / ref_sum[1], full_sum[2] / ref_sum[2], full_sum[3] / ref_sum[3] },
            { filt_sum[0] / ref_sum[0], filt_sum[1] / ref_sum[1], filt_sum[2] / ref_sum[2], filt_sum[3] / ref_sum[3] },
        };

        for (uint32_t i = 0; i < 3; ++i) {
            if (display_mode[i] == IsoLevels || display_mode[i] == IsoLines) {
                const uint32_t count = 3;
                uint32_t colors[4][count] = {0};
                uint32_t lines[4][count] = {0};
                float values[4][count] = {0};

                memcpy(values, ref_iso_values, sizeof(ref_iso_values));
                for (uint32_t j = 0; j < count; ++j) {
                    values[0][j] *= MAX(scl[i][0], FLT_EPSILON);
                    values[1][j] *= MAX(scl[i][1], FLT_EPSILON);
                    values[2][j] *= MAX(scl[i][2], FLT_EPSILON);
                    values[3][j] *= MAX(scl[i][3], FLT_EPSILON);
                }

                if (display_mode[i] == IsoLevels) {
                    memcpy(colors, ref_iso_level_colors, sizeof(colors));
                    memcpy(lines,  ref_iso_level_colors, sizeof(lines));
                }
                else {
                    uint32_t line_colors[3] = {
                        0,
                        convert_color(vec4_from_vec3(hcl_to_rgb(rgb_to_hcl(vec3_from_vec4(isoline_colors[i])) * vec3_t { 1.0f, 0.5f, 1.0f }), isoline_colors[i].w)),
                        convert_color(vec4_from_vec3(hcl_to_rgb(rgb_to_hcl(vec3_from_vec4(isoline_colors[i])) * vec3_t { 1.0f, 1.0f, 1.0f }), isoline_colors[i].w))
                    };
                    for (int j = 0; j < 4; ++j) {
                        lines[j][0] = line_colors[0];
                        lines[j][1] = line_colors[1];
                        lines[j][2] = line_colors[2];
                    }
                }

                rama_isomap_t maps[4] = {
                    {
                        .values = values[0],
                        .level_colors = colors[0],
                        .contour_colors = lines[0],
                        .count = count,
                    },
                    {
                        .values = values[1],
                        .level_colors = colors[1],
                        .contour_colors = lines[1],
                        .count = count,
                    },
                    {
                        .values = values[2],
                        .level_colors = colors[2],
                        .contour_colors = lines[2],
                        .count = count,
                    },
                    {
                        .values = values[3],
                        .level_colors = colors[3],
                        .contour_colors = lines[3],
                        .count = count,
                    },
                };

                rama_rep_render_iso(reps[i], viewport.elem, maps, tex_dim);
            } else if (display_mode[i] == Colormap){
                // Fill in colors based on active color_map
                uint32_t colors[32] = {0};
                uint32_t num_colors = ImPlot::GetColormapSize(colormap[i]);
                for (uint32_t j = 0; j < num_colors; ++j) {
                    colors[j] = ImPlot::GetColormapColorU32(j, colormap[i]);
                }
                // Mask off the first color to be transparent
                colors[0] = colors[0] & 0x00FFFFFFU;

                rama_colormap_t maps[4] = {
                    {
                        .colors = colors,
                        .count = num_colors,
                        .min_value = 0,
                        .max_value = 0.5f * scl[i][0],
                    },
                    {
                        .colors = colors,
                        .count = num_colors,
                        .min_value = 0,
                        .max_value = 0.5f * scl[i][1],
                    },
                    {
                        .colors = colors,
                        .count = num_colors,
                        .min_value = 0,
                        .max_value = 0.5f * scl[i][2],
                    },
                    {
                        .colors = colors,
                        .count = num_colors,
                        .min_value = 0,
                        .max_value = 0.5f * scl[i][3],
                    },
                };
                rama_rep_render_map(reps[i], viewport.elem, maps, tex_dim);
            }
        }

        POP_GPU_SECTION();

        ImPlot::MapInputDefault();
    }
    ImGui::End();
}

static void draw_density_volume_window(ApplicationData* data) {
    ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Density Volume", &data->density_volume.show_window, ImGuiWindowFlags_MenuBar)) {
        const ImVec2 button_size = {160, 0};

        if (ImGui::IsWindowFocused() && ImGui::IsKeyPressed(KEY_PLAY_PAUSE, false)) {
            data->animation.mode = data->animation.mode == PlaybackMode::Playing ? PlaybackMode::Stopped : PlaybackMode::Playing;
        }

        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("Property")) {
                draw_property_menu_widgets(data->display_properties, md_array_size(data->display_properties), DisplayProperty::ShowIn_Volume, false);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("DVR")) {
                ImGui::Checkbox("Enabled", &data->density_volume.dvr.enabled);

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

                ImGui::SliderFloat("Density Scaling", &data->density_volume.dvr.density_scale, 0.001f, 10000.f, "%.3f", ImGuiSliderFlags_Logarithmic);
                ImGui::SliderFloat("Alpha Scaling", &data->density_volume.dvr.tf.alpha_scale, 0.001f, 10.f, "%.3f", ImGuiSliderFlags_Logarithmic);
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("ISO")) {
                ImGui::Checkbox("Enabled", &data->density_volume.iso.enabled);
                for (int i = 0; i < data->density_volume.iso.count; ++i) {
                    ImGui::PushID(i);
                    ImGui::SliderFloat("##Isovalue", &data->density_volume.iso.values[i], 0.0f, 10.f, "%.3f", ImGuiSliderFlags_Logarithmic);
                    if (ImGui::IsItemDeactivatedAfterEdit()) {
                        // @TODO(Robin): Sort?
                    }
                    ImGui::SameLine();
                    ImGui::ColorEdit4Minimal("##Color", data->density_volume.iso.colors[i].elem);
                    ImGui::PopID();
                }
                if ((data->density_volume.iso.count < (int)ARRAY_SIZE(data->density_volume.iso.values)) &&
                    ImGui::Button("Add", button_size)) {
                    int idx = data->density_volume.iso.count++;
                    data->density_volume.iso.values[idx] = 0.1f;
                    data->density_volume.iso.colors[idx] = { 0.2f, 0.1f, 0.9f, 1.0f };
                    // @TODO(Robin): Sort?
                }
                if (ImGui::Button("Clear", button_size)) {
                    data->density_volume.iso.count = 0;
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
                    auto& rep = data->density_volume.rep;
                    ImGui::Checkbox("Show Superimposed Structures", &data->density_volume.show_reference_ensemble);
                    if (ImGui::Combo("type", (int*)(&rep.type), "VDW\0Licorice\0Ribbons\0Cartoon\0")) {
                        data->density_volume.dirty_rep = true;
                    }
                    if (ImGui::Combo("color", (int*)(&rep.colormap), "Uniform Color\0CPK\0Atom Label\0Atom Idx\0Res Id\0Res Idx\0Chain Id\0Chain Idx\0Secondary Structure\0")) {
                        data->density_volume.dirty_rep = true;
                    }
                    if (rep.colormap == ColorMapping::Uniform) {
                        data->density_volume.dirty_rep |= ImGui::ColorEdit4("color", rep.color.elem, ImGuiColorEditFlags_NoInputs);
                    }
                    if (rep.type == RepresentationType::SpaceFill || rep.type == RepresentationType::Licorice) {
                        data->density_volume.dirty_rep |= ImGui::SliderFloat("scale", &rep.radius, 0.1f, 2.f);
                    }
                    if (rep.type == RepresentationType::Ribbons) {
                        data->density_volume.dirty_rep |= ImGui::SliderFloat("spline tension", &rep.tension, 0.f, 1.f);
                        data->density_volume.dirty_rep |= ImGui::SliderFloat("spline width", &rep.width, 0.1f, 2.f);
                        data->density_volume.dirty_rep |= ImGui::SliderFloat("spline thickness", &rep.thickness, 0.1f, 2.f);
                    }
                }
                ImGui::EndMenu();
            }

            ImGui::EndMenuBar();
        }

        update_density_volume(data);

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
        ImGui::InvisibleButton("canvas", canvas_sz, ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);

        // Draw border and background color
        ImGuiIO& io = ImGui::GetIO();

        ImVec2 canvas_p0 = ImGui::GetItemRectMin();
        ImVec2 canvas_p1 = ImGui::GetItemRectMax();

        ImDrawList* draw_list = ImGui::GetWindowDrawList();
        draw_list->AddImage((ImTextureID)(uint64_t)data->density_volume.fbo.deferred.post_tonemap, canvas_p0, canvas_p1, { 0,1 }, { 1,0 });
        draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));

        const bool is_hovered = ImGui::IsItemHovered();
        const bool is_active = ImGui::IsItemActive();
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
            camera_controller_trackball(&data->density_volume.camera.position, &data->density_volume.camera.orientation, &data->density_volume.camera.focus_distance, input, param);
        }

        mat4_t view_mat = camera_world_to_view_matrix(data->density_volume.camera);
        mat4_t proj_mat = camera_perspective_projection_matrix(data->density_volume.camera, (float)canvas_sz.x / (float)canvas_sz.y);
        mat4_t inv_proj_mat = camera_inverse_perspective_projection_matrix(data->density_volume.camera, (float)canvas_sz.x / (float)canvas_sz.y);

        clear_gbuffer(&gbuf);

        glDrawBuffer(GL_COLOR_ATTACHMENT_POST_TONEMAP);
        glClearColor(1, 1, 1, 1);
        glClear(GL_COLOR_BUFFER_BIT);

        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);

        const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
            GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_POST_TONEMAP };

        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);

        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.deferred.fbo);
        glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
        glViewport(0, 0, gbuf.width, gbuf.height);

        int64_t num_reps = md_array_size(data->density_volume.gl_reps);
        if (data->density_volume.show_reference_structures && num_reps > 0) {
            num_reps = data->density_volume.show_reference_ensemble ? num_reps : 1;

            md_gl_draw_op_t* draw_ops = 0;
            for (int64_t i = 0; i < num_reps; ++i) {
                md_gl_draw_op_t op = {};
                op.rep = &data->density_volume.gl_reps[i];
                op.model_matrix = &data->density_volume.rep_model_mats[i].elem[0][0];
                md_array_push(draw_ops, op, frame_allocator);
            }

            md_gl_draw_args_t draw_args = {
                .shaders = &data->mold.gl_shaders,
                .draw_operations = {
                    .count = (uint32_t)md_array_size(draw_ops),
                    .ops = draw_ops
                },
                .view_transform = {
                    .view_matrix = &view_mat.elem[0][0],
                    .projection_matrix = &proj_mat.elem[0][0],
                },
            };

            md_gl_draw(&draw_args);

            if (is_hovered) {
                vec2_t coord = {mouse_pos_in_canvas.x, (float)gbuf.height - mouse_pos_in_canvas.y};
                PickingData pd = read_picking_data(&gbuf, (int)coord.x, (int)coord.y);
                if (pd.idx != INVALID_PICKING_IDX) {
                    draw_atom_info_window(*data, pd.idx);
                }
            }

            glDrawBuffer(GL_COLOR_ATTACHMENT_POST_TONEMAP);

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
                            .focus_depth = data->visuals.dof.focus_depth,
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
                    .norm = view_mat,
                },
                .inverse = {
                    .proj = inv_proj_mat,
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

        glDrawBuffer(GL_COLOR_ATTACHMENT_POST_TONEMAP);
        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);

        if (data->density_volume.show_bounding_box) {
            if (data->density_volume.model_mat != mat4_t{0}) {
                const vec3_t min_box = {0,0,0};
                const vec3_t max_box = {1,1,1};

                immediate::set_model_view_matrix(mat4_mul(view_mat, data->density_volume.model_mat));
                immediate::set_proj_matrix(proj_mat);

                uint32_t box_color = convert_color(data->density_volume.bounding_box_color);
                uint32_t clip_color = convert_color(data->density_volume.clip_volume_color);
                immediate::draw_box_wireframe(min_box, max_box, box_color);
                immediate::draw_box_wireframe(data->density_volume.clip_volume.min, data->density_volume.clip_volume.max, clip_color);

                immediate::render();
            }
        }

        if (data->density_volume.show_density_volume) {
            if (data->density_volume.model_mat != mat4_t{ 0 }) {
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
                        .proj = proj_mat,
                    },
                    .clip_volume = {
                        .min = data->density_volume.clip_volume.min,
                        .max = data->density_volume.clip_volume.max,
                    },
                    .bounding_box = {
                        .color = data->density_volume.bounding_box_color,
                        .enabled = data->density_volume.show_bounding_box,
                    },
                    .global_scaling = {
                        .density = data->density_volume.dvr.density_scale,
                        .alpha = data->density_volume.dvr.tf.alpha_scale,
                    },
                    .iso_surface = {
                        .count = data->density_volume.iso.count,
                        .values = data->density_volume.iso.values,
                        .colors = data->density_volume.iso.colors,
                    },
                    .isosurface_enabled = data->density_volume.iso.enabled,
                    .direct_volume_rendering_enabled = data->density_volume.dvr.enabled,

                    .voxel_spacing = data->density_volume.voxel_spacing
                };
                volume::render_volume(vol_desc);
            }
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
        ImGui::Text("Num atoms:    %9i", (int)data->mold.mol.atom.count);
        ImGui::Text("Num residues: %9i", (int)data->mold.mol.residue.count);
        ImGui::Text("Num chains:   %9i", (int)data->mold.mol.chain.count);

        str_t traj_file = data->files.trajectory;
        if (traj_file.len) {
            ImGui::Separator();
            ImGui::Text("Trajectory data: %.*s", (int)traj_file.len, traj_file.ptr);
            ImGui::Text("Num frames:    %9i", (int)md_trajectory_num_frames(data->mold.traj));
            ImGui::Text("Num atoms:     %9i", (int)md_trajectory_num_atoms(data->mold.traj));
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

        
        task_system::ID* tasks = task_system::pool_running_tasks(default_allocator);
        int64_t num_tasks = md_array_size(tasks);
        if (num_tasks > 0) {
            ImGui::Text("Running Pool Tasks:");
            for (int64_t i = 0; i < md_array_size(tasks); ++i) {
                str_t lbl = task_system::task_label(tasks[i]);
                ImGui::Text("[%i]: %.*s", (int)i, (int)lbl.len, lbl.ptr);
            }
        }

        ImGuiID active = ImGui::GetActiveID();
        ImGuiID hover  = ImGui::GetHoveredID();
        ImGui::Text("Active ID: %u, Hover ID: %u", active, hover);
    }
    ImGui::End();
}

static void draw_script_editor_window(ApplicationData* data) {
    ASSERT(data);

    ImGui::SetNextWindowSize({300,200}, ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Script Editor", &data->show_script_window, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_MenuBar)) {
        ImGui::SetWindowSize(ImVec2(800, 600), ImGuiCond_FirstUseEver);
        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("File")) {
                char path_buf[1024] = "";
                if (ImGui::MenuItem("Load")) {
                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialog_Open, "txt")) {
                        str_t txt = load_textfile(str_from_cstr(path_buf), frame_allocator);
                        std::string str(txt.ptr, txt.len);
                        editor.SetText(str);
                    }
                }
                if (ImGui::MenuItem("Save")) {
                    auto textToSave = editor.GetText();
                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialog_Save, "txt")) {
                        int path_len = (int)strnlen(path_buf, sizeof(path_buf));
                        str_t path = str_t{path_buf, path_len};
                        if (str_empty(extract_ext(path))) {
                            path.len += snprintf(path_buf + path_len, sizeof(path_buf) - path_len, ".txt");
                        }
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
            data->mold.script.compile_ir = true;
            data->mold.script.time_since_last_change = 0;
        }

        const ImVec2 content_size = ImGui::GetContentRegionAvail();
        const char* btn_text = "Evaluate";
        const ImVec2 label_size = ImGui::CalcTextSize(btn_text, NULL, true) * ImVec2(1.4, 1.0);
        const ImVec2 btn_size = ImGui::CalcItemSize(ImVec2(0,0), label_size.x + ImGui::GetStyle().FramePadding.x * 2.0f, label_size.y + ImGui::GetStyle().FramePadding.y * 2.0f);
        const ImVec2 text_size(content_size - ImVec2(0, btn_size.y + ImGui::GetStyle().ItemSpacing.y));

        editor.Render("TextEditor", text_size);

        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + content_size.x - btn_size.x);
        ImGui::SetCursorPosY(ImGui::GetCursorPosY() + ImGui::GetStyle().ItemSpacing.y);

        const bool valid = md_script_ir_valid(data->mold.script.ir) && !task_system::task_is_running(data->tasks.evaluate_full);
        
        if (!valid) ImGui::PushDisabled();
        if (ImGui::Button(btn_text, btn_size)) {
            data->mold.script.eval_init = true;
        }
        if (!valid) ImGui::PopDisabled();

        data->mold.script.vis = {0};
        const TextEditor::Marker* hovered_marker = editor.GetHoveredMarker();
        if (hovered_marker) {
            if (md_semaphore_try_aquire(&data->mold.script.ir_semaphore)) {
                defer { md_semaphore_release(&data->mold.script.ir_semaphore); };
                
                if (md_script_ir_valid(data->mold.script.ir)) {
                    md_script_visualization_init(&data->mold.script.vis, {
                        .payload = (const md_script_vis_payload_t*)hovered_marker->payload,
                        .ir = data->mold.script.ir,
                        .mol = &data->mold.mol,
                        .traj = data->mold.traj,
                        .alloc = frame_allocator,
                    });
                    if (!md_bitfield_empty(data->mold.script.vis.atom_mask)) {
                        md_bitfield_copy(&data->selection.current_highlight_mask, data->mold.script.vis.atom_mask);
                        data->mold.dirty_buffers |= MolBit_DirtyFlags;
                    }
                }
            }
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
        LOG_ERROR("Failed to open file '%.*s' to write data.", (int)filename.len, filename.ptr);
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
        LOG_ERROR("Failed to open file '%.*s' to write data.", (int)filename.len, filename.ptr);
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

static bool export_cube(ApplicationData& data, const md_script_property_t* prop, str_t filename) {
    // @NOTE: First we need to extract some meta data for the cube format, we need the atom indices/bits for any SDF
    // And the origin + extent of the volume in spatial coordinates (Ångström)

    if (!prop) {
        LOG_ERROR("Export Cube: The property to be exported did not exist");
        return false;
    }

    // Copy mol and replace with initial coords
    md_molecule_t mol = data.mold.mol;

    int64_t stride = ALIGN_TO(data.mold.mol.atom.count, md_simd_f32_width);
    float* coords = (float*)md_alloc(frame_allocator, stride * sizeof(float) * 3);
    mol.atom.x = coords + stride * 0;
    mol.atom.y = coords + stride * 1;
    mol.atom.z = coords + stride * 2;
    
    if (!md_trajectory_load_frame(data.mold.traj, 0, NULL, mol.atom.x, mol.atom.y, mol.atom.z)) {
        return false;
    }

    bool result = false;
    md_script_visualization_t vis = { 0 };
    defer { md_script_visualization_free(&vis); };
    
    if (md_semaphore_aquire(&data.mold.script.ir_semaphore)) {
        defer { md_semaphore_release(&data.mold.script.ir_semaphore); };
        
        if (md_script_ir_valid(data.mold.script.ir)) {
            md_script_visualization_args_t args = {
                .payload = prop->vis_payload,
                .ir = data.mold.script.ir,
                .mol = &data.mold.mol,
                .traj = data.mold.traj,
                .alloc = frame_allocator,
                .flags = MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_SDF,
            };
            result = md_script_visualization_init(&vis, args);
        }
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

        if (vis.sdf.count > 0) {
            const float angstrom_to_bohr = (float)(1.0 / 0.529177210903);

            // transformation matrix from world to volume
            mat4_t M = vis.sdf.matrices[0];
            const md_bitfield_t* bf = &vis.sdf.structures[0];
            const int num_atoms = (int)md_bitfield_popcount(bf);
            const int vol_dim[3] = {prop->data.dim[0], prop->data.dim[1], prop->data.dim[2]};
            const double extent = vis.sdf.extent * 2.0 * angstrom_to_bohr;
            const double voxel_ext[3] = {
                (double)extent / (double)vol_dim[0],
                (double)extent / (double)vol_dim[1],
                (double)extent / (double)vol_dim[2],
            };

            const double half_ext = extent * 0.5;

            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", -num_atoms, -half_ext, -half_ext, -half_ext);
            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol_dim[0], voxel_ext[0], 0, 0);
            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol_dim[1], 0, voxel_ext[1], 0);
            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol_dim[2], 0, 0, voxel_ext[2]);

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
                        float val = prop->data.values[idx];
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

static void draw_property_export_window(ApplicationData* data) {
    ASSERT(data);

    enum ExportType {
        Temporal = 0,
        Distribution,
        DensityVolume
    };

    struct ExportFormat {
        const char* label;
        const char* extension;
    };

    ExportFormat table_formats[] {
        {"XVG", "xvg"},
        {"CSV", "csv"}
    };

    ExportFormat volume_formats[] {
        {"Gaussian Cube", "cube"},
    };

    if (ImGui::Begin("Property Export", &data->show_property_export_window)) {
        char path_buf[1024] = "";
        static ExportType type = Temporal;
        ImGui::PushItemWidth(200);
        ImGui::Combo("Data Type", (int*)(&type), "Temporal\0Distribution\0Density Volume\0");
        ImGui::Separator();
        if (type == Temporal || type == Distribution) {
            static int format = 0;
            static int col_options[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
            static int num_columns = 4;
            static int num_bins = 64;

            const int MIN_NUM_BINS = 8;
            const int MAX_NUM_BINS = 1024;

            struct ColData {
                const char* label;
                const float* values;
                int32_t num_values;
                int32_t dim;
            };

            ColData* col_data = 0;
            if (type == Temporal) {
                ColData time_data = {"Time", data->timeline.x_values, (int)md_array_size(data->timeline.x_values), 1};
                md_array_push(col_data, time_data, frame_allocator);
                for (int i = 0; i < (int)md_array_size(data->display_properties); ++i) {
                    const DisplayProperty& dp = data->display_properties[i];
                    if (dp.full_prop->flags & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) {
                        str_t lbl = {dp.label, (int64_t)strnlen(dp.label, sizeof(dp.label))};
                        if (dp.full_prop->data.dim[0] > 1) {
                            lbl = alloc_printf(frame_allocator, "%.*s[%d:%d]", (int)lbl.len, lbl.ptr, 1, dp.full_prop->data.dim[0]);
                        }
                        ColData prop_data = {lbl.ptr, dp.full_prop->data.values, (int)dp.full_prop->data.num_values, dp.full_prop->data.dim[0]};
                        md_array_push(col_data, prop_data, frame_allocator);
                    }
                }
            } else if (type == Distribution) {
                for (int i = 0; i < (int)md_array_size(data->display_properties); ++i) {
                    const DisplayProperty& dp = data->display_properties[i];
                    if (dp.full_prop->flags & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) {
                        ColData prop_data = {dp.label, dp.full_hist.bin, (int)ARRAY_SIZE(dp.full_hist.bin), 1};
                        md_array_push(col_data, prop_data, frame_allocator);
                        if (data->timeline.filter.enabled) {
                            str_t lbl = alloc_printf(frame_allocator, "%s(filt)", dp.label);
                            ColData filt_data = {lbl.ptr, dp.filt_hist.bin, (int)ARRAY_SIZE(dp.filt_hist.bin), 1};
                            md_array_push(col_data, filt_data, frame_allocator);
                        }
                    } else if (dp.full_prop->flags & MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) {
                        ColData prop_data = {dp.label, dp.full_prop->data.values, (int)dp.full_prop->data.num_values, dp.full_prop->data.dim[0]};
                        md_array_push(col_data, prop_data, frame_allocator);
                        if (data->timeline.filter.enabled) {
                            str_t lbl = alloc_printf(frame_allocator, "%s(filt)", dp.label);
                            ColData filt_data = {lbl.ptr, dp.filt_prop->data.values, (int)dp.filt_prop->data.num_values, dp.filt_prop->data.dim[0]};
                            md_array_push(col_data, filt_data, frame_allocator);
                        }
                    }
                }
            }

            ImGui::Text("Column Layout");
            ImGui::SliderInt("Num Columns", &num_columns, 1, (int)ARRAY_SIZE(col_options));

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

                //application::FileDialogResult file_dialog = application::file_dialog(application::FileDialog_Save, {}, table_formats[format].extension);
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialog_Save, table_formats[format].extension)) {
                    int path_len = (int)strnlen(path_buf, sizeof(path_buf));
                    if (str_empty(extract_ext({path_buf, path_len}))) {
                        path_len += snprintf(path_buf + path_len, sizeof(path_buf) - path_len, ".%s", table_formats[format].extension);
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
                                        str_t lbl = alloc_printf(frame_allocator, "%s[%i]", col_data[idx].label, j + 1);
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

                    int num_rows = (type == Distribution) ? num_bins : md_trajectory_num_frames(data->mold.traj);

                    switch (format) {
                    case 0:
                        export_xvg(column_data, column_labels, md_array_size(column_data), num_rows, {path_buf, path_len});
                        break;
                    case 1:
                        export_csv(column_data, column_labels, md_array_size(column_data), num_rows, {path_buf, path_len});
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
                if (data->display_properties[i].full_prop->flags & MD_SCRIPT_PROPERTY_FLAG_VOLUME) {
                    md_array_push(props, &data->display_properties[i], frame_allocator);
                }
            }

            const int num_props = (int)md_array_size(props);
            if (num_props > 0) {
                prop_idx = CLAMP(prop_idx, 0, num_props);
                if (ImGui::BeginCombo("Source", props[prop_idx]->label)) {
                    for (int i = 0; i < num_props; ++i) {
                        if (ImGui::Selectable(props[i]->label, prop_idx == i)) {
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

                format = CLAMP(format, 0, (int)ARRAY_SIZE(volume_formats));
                if (ImGui::BeginCombo("Format", volume_formats[format].label)) {
                    for (int i = 0; i < (int)ARRAY_SIZE(volume_formats); ++i) {
                        if (ImGui::Selectable(volume_formats[i].label, format == i)) format = i;
                    }
                    ImGui::EndCombo();
                }

                if (ImGui::Button("Export")) {
                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialog_Save, volume_formats[format].extension)) {
                        int path_len = (int)strnlen(path_buf, sizeof(path_buf));
                        if (str_empty(extract_ext({path_buf, path_len}))) {
                            path_len += snprintf(path_buf + path_len, sizeof(path_buf) - path_len, ".%s", volume_formats[format].extension);
                        }
                        switch (format) {
                        case 0:
                            export_cube(*data, temporal_filter ? props[prop_idx]->filt_prop : props[prop_idx]->full_prop, {path_buf, path_len});
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
    if (!gbuf->pbo_picking.color[0]) glGenBuffers((int)ARRAY_SIZE(gbuf->pbo_picking.color), gbuf->pbo_picking.color);
    if (!gbuf->pbo_picking.depth[0]) glGenBuffers((int)ARRAY_SIZE(gbuf->pbo_picking.depth), gbuf->pbo_picking.depth);

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
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG8, width, height, 0, GL_RG, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, gbuf->deferred.velocity);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16F, width, height, 0, GL_RG, GL_FLOAT, nullptr);
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

    for (uint32_t i = 0; i < ARRAY_SIZE(gbuf->pbo_picking.color); ++i) {
        glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.color[i]);
        glBufferData(GL_PIXEL_PACK_BUFFER, 4, NULL, GL_DYNAMIC_READ);
        glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    }

    for (uint32_t i = 0; i < ARRAY_SIZE(gbuf->pbo_picking.depth); ++i) {
        glBindBuffer(GL_PIXEL_PACK_BUFFER, gbuf->pbo_picking.depth[i]);
        glBufferData(GL_PIXEL_PACK_BUFFER, 4, NULL, GL_DYNAMIC_READ);
        glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    }

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
    glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
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
        const vec3_t pbc_ext = data->simulation_box.box * vec3_t{1,1,1};
        md_gl_molecule_set_atom_position(&data->mold.gl_mol, 0, (uint32_t)mol.atom.count, mol.atom.x, mol.atom.y, mol.atom.z, 0);
        md_gl_molecule_compute_velocity(&data->mold.gl_mol, pbc_ext.elem);
#if EXPERIMENTAL_GFX_API
        md_gfx_structure_set_atom_position(data->mold.gfx_structure, 0, (uint32_t)mol.atom.count, mol.atom.x, mol.atom.y, mol.atom.z, 0);
        md_gfx_structure_set_aabb(data->mold.gfx_structure, &data->mold.mol_aabb_min, &data->mold.mol_aabb_max);
#endif
    }

    if (data->mold.dirty_buffers & MolBit_DirtyRadius) {
        md_gl_molecule_set_atom_radius(&data->mold.gl_mol, 0, (uint32_t)mol.atom.count, mol.atom.radius, 0);
#if EXPERIMENTAL_GFX_API
        md_gfx_structure_set_atom_radius(data->mold.gfx_structure, 0, (uint32_t)mol.atom.count, mol.atom.radius, 0);
#endif
    }

    if (data->mold.dirty_buffers & MolBit_DirtyFlags) {
        if (data->mold.mol.atom.flags) {
            for (int64_t i = 0; i < mol.atom.count; i++) {
                md_flags_t flags = 0;
                flags |= md_bitfield_test_bit(&data->selection.current_highlight_mask, i)     ? AtomBit_Highlighted : 0;
                flags |= md_bitfield_test_bit(&data->selection.current_selection_mask, i)     ? AtomBit_Selected : 0;
                flags |= md_bitfield_test_bit(&data->representation.atom_visibility_mask, i) ? AtomBit_Visible : 0;
                data->mold.mol.atom.flags[i] = flags;
            }
            md_gl_molecule_set_atom_flags(&data->mold.gl_mol, 0, (uint32_t)mol.atom.count, mol.atom.flags, 0);
        }
    }

    if (data->mold.dirty_buffers & MolBit_DirtyBonds) {
        md_gl_molecule_set_covalent_bonds(&data->mold.gl_mol, 0, (uint32_t)mol.covalent.count, mol.covalent.bond, 0);
    }

    if (data->mold.dirty_buffers & MolBit_DirtySecondaryStructure) {
        if (mol.backbone.secondary_structure) {
            md_gl_molecule_set_backbone_secondary_structure(&data->mold.gl_mol, 0, (uint32_t)mol.backbone.count, mol.backbone.secondary_structure, 0);
        }
    }

    data->mold.dirty_buffers = 0;
}

static void interrupt_async_tasks(ApplicationData* data) {
    task_system::pool_interrupt_running_tasks();

    if (data->mold.script.full_eval) md_script_eval_interrupt(data->mold.script.full_eval);
    if (data->mold.script.filt_eval) md_script_eval_interrupt(data->mold.script.filt_eval);

    task_system::task_wait_for(data->tasks.backbone_computations);
    task_system::task_wait_for(data->tasks.evaluate_full);
    task_system::task_wait_for(data->tasks.evaluate_filt);
    task_system::task_wait_for(data->tasks.prefetch_frames);
    task_system::task_wait_for(data->tasks.ramachandran_compute_full_density);
    task_system::task_wait_for(data->tasks.ramachandran_compute_filt_density);
    task_system::task_wait_for(data->tasks.shape_space_evaluate);
}

// #trajectorydata
static void free_trajectory_data(ApplicationData* data) {
    ASSERT(data);
    interrupt_async_tasks(data);

    if (md_trajectory_num_frames(data->mold.traj)) {
        load::traj::close(data->mold.traj);
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
    int64_t num_frames = md_trajectory_num_frames(data->mold.traj);
    if (num_frames > 0) {
        int64_t min_frame = 0;
        int64_t max_frame = num_frames - 1;
        double min_time = 0;
        double max_time = 0;
        md_trajectory_frame_header_t header;

        {
            md_trajectory_load_frame(data->mold.traj, min_frame, &header, 0, 0, 0);
            min_time = header.timestamp;
        }

        {
            md_trajectory_load_frame(data->mold.traj, max_frame, &header, 0, 0, 0);
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

        md_trajectory_load_frame(data->mold.traj, frame_idx, &header, data->mold.mol.atom.x, data->mold.mol.atom.y, data->mold.mol.atom.z);
        data->simulation_box.box = header.box;

        if (data->mold.mol.backbone.count > 0) {
            data->trajectory_data.secondary_structure.stride = data->mold.mol.backbone.count;
            data->trajectory_data.secondary_structure.count = data->mold.mol.backbone.count * num_frames;
            md_array_resize(data->trajectory_data.secondary_structure.data, data->mold.mol.backbone.count * num_frames, persistent_allocator);
            memset(data->trajectory_data.secondary_structure.data, 0, md_array_size(data->trajectory_data.secondary_structure.data) * sizeof (md_secondary_structure_t));

            data->trajectory_data.backbone_angles.stride = data->mold.mol.backbone.count;
            data->trajectory_data.backbone_angles.count = data->mold.mol.backbone.count * num_frames;
            md_array_resize(data->trajectory_data.backbone_angles.data, data->mold.mol.backbone.count * num_frames, persistent_allocator);
            memset(data->trajectory_data.backbone_angles.data, 0, md_array_size(data->trajectory_data.backbone_angles.data) * sizeof (md_backbone_angles_t));

            // Launch work to compute the values
            task_system::task_interrupt_and_wait_for(data->tasks.backbone_computations);

            data->tasks.backbone_computations = task_system::pool_enqueue(STR("Backbone Operations"), 0, (uint32_t)num_frames, [](uint32_t range_beg, uint32_t range_end, void* user_data) {
                ApplicationData* data = (ApplicationData*)user_data;
                
                // Create copy here of molecule since we use the full structure as input
                md_molecule_t mol = data->mold.mol;

                const int64_t stride = ALIGN_TO(mol.atom.count, md_simd_f32_width);
                const int64_t bytes = stride * sizeof(float) * 3;
                float* coords = (float*)md_alloc(default_allocator, bytes);
                defer { md_free(default_allocator, coords, bytes); };
                // Overwrite the coordinate section, since we will load trajectory frame data into these
                mol.atom.x = coords + stride * 0;
                mol.atom.y = coords + stride * 1;
                mol.atom.z = coords + stride * 2;

                for (uint32_t frame_idx = range_beg; frame_idx < range_end; ++frame_idx) {
                    md_trajectory_load_frame(data->mold.traj, frame_idx, NULL, mol.atom.x, mol.atom.y, mol.atom.z);
                    md_util_backbone_angles_compute(data->trajectory_data.backbone_angles.data + data->trajectory_data.backbone_angles.stride * frame_idx, data->trajectory_data.backbone_angles.stride, &mol);
                    md_util_backbone_secondary_structure_compute(data->trajectory_data.secondary_structure.data + data->trajectory_data.secondary_structure.stride * frame_idx, data->trajectory_data.secondary_structure.stride, &mol);
                }
            }, data);

            task_system::main_enqueue(STR("Update Trajectory Data"), [](void* user_data) {
                ApplicationData* data = (ApplicationData*)user_data;
                data->trajectory_data.backbone_angles.fingerprint = generate_fingerprint();
                data->trajectory_data.secondary_structure.fingerprint = generate_fingerprint();
                interpolate_atomic_properties(data);
            }, data, data->tasks.backbone_computations);
        }

        data->mold.dirty_buffers |= MolBit_DirtyPosition;
        update_md_buffers(data);
        md_gl_molecule_zero_velocity(&data->mold.gl_mol); // Do this explicitly to update the previous position to avoid motion blur trails

                                                          // Prefetch frames
        launch_prefetch_job(data);
    }
}

static bool load_trajectory_data(ApplicationData* data, str_t filename, bool deperiodize_on_load) {
    md_trajectory_i* traj = load::traj::open_file(filename, &data->mold.mol, persistent_allocator, deperiodize_on_load);
    if (traj) {
        free_trajectory_data(data);
        data->mold.traj = traj;
        data->files.trajectory = filename;
        data->files.deperiodize = deperiodize_on_load;
        init_trajectory_data(data);
        data->animation.frame = 0;
        return true;
    }

    return false;
}

// #moleculedata
static void free_molecule_data(ApplicationData* data) {
    ASSERT(data);
    interrupt_async_tasks(data);

    //md_molecule_free(&data->mold.mol, persistent_allocator);
    md_arena_allocator_reset(data->mold.mol_alloc);
    memset(&data->mold.mol, 0, sizeof(data->mold.mol));

    md_gl_molecule_free(&data->mold.gl_mol);

    data->files.molecule = "";
    free_trajectory_data(data);

    md_bitfield_clear(&data->selection.current_selection_mask);
    md_bitfield_clear(&data->selection.current_highlight_mask);
    if (data->mold.script.ir) {
        md_script_ir_clear(data->mold.script.ir);
    }
    if (data->mold.script.full_eval) {
        md_script_eval_free(data->mold.script.full_eval);
        data->mold.script.full_eval = 0;
    }
    if (data->mold.script.filt_eval) {
        md_script_eval_free(data->mold.script.filt_eval);
        data->mold.script.filt_eval = 0;
    }
    clear_density_volume(data);
}

static void init_molecule_data(ApplicationData* data) {
    if (data->mold.mol.atom.count) {
        const md_molecule_t& mol = data->mold.mol;
        vec3_t& aabb_min = data->mold.mol_aabb_min;
        vec3_t& aabb_max = data->mold.mol_aabb_max;
        md_util_compute_aabb_xyzr(&aabb_min, &aabb_max, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.radius, mol.atom.count);

        data->picking.idx = INVALID_PICKING_IDX;
        data->selection.hovered = -1;
        data->selection.right_clicked = -1;

        md_gl_molecule_init(&data->mold.gl_mol, &data->mold.mol);

#if EXPERIMENTAL_GFX_API
        data->mold.gfx_structure = md_gfx_structure_create(mol.atom.count, mol.covalent.count, mol.backbone.count, mol.backbone.range_count, mol.residue.count, mol.instance.count);
        md_gfx_structure_set_atom_position(data->mold.gfx_structure, 0, mol.atom.count, mol.atom.x, mol.atom.y, mol.atom.z, 0);
        md_gfx_structure_set_atom_radius(data->mold.gfx_structure, 0, mol.atom.count, mol.atom.radius, 0);
        md_gfx_structure_set_aabb(data->mold.gfx_structure, &data->mold.mol_aabb_min, &data->mold.mol_aabb_max);
        if (mol.instance.count > 0) {
            md_gfx_structure_set_instance_atom_ranges(data->mold.gfx_structure, 0, mol.instance.count, (md_gfx_range_t*)mol.instance.atom_range, 0);
            md_gfx_structure_set_instance_transforms(data->mold.gfx_structure, 0, mol.instance.count, mol.instance.transform, 0);
        }
#endif
        init_all_representations(data);
        update_all_representations(data);
        data->mold.script.compile_ir = true;
    }
}

static void launch_prefetch_job(ApplicationData* data) {
    uint32_t num_frames = MIN((uint32_t)md_trajectory_num_frames(data->mold.traj), (uint32_t)load::traj::num_cache_frames(data->mold.traj));
    if (!num_frames) return;

    task_system::task_interrupt_and_wait_for(data->tasks.prefetch_frames);
    data->tasks.prefetch_frames = task_system::pool_enqueue(STR("Prefetch Frames"), 0, num_frames, [](uint32_t range_beg, uint32_t range_end, void* user_data) {
        ApplicationData* data = (ApplicationData*)user_data;
        for (uint32_t i = range_beg; i < range_end; ++i) {
            md_trajectory_frame_header_t header;
            md_trajectory_load_frame(data->mold.traj, i, &header, 0, 0, 0);
            data->timeline.x_values[i] = header.timestamp;
        }
    }, data);

    task_system::main_enqueue(STR("Prefetch Complete"), [](void* user_data) {
        ApplicationData* data = (ApplicationData*)user_data;
        interpolate_atomic_properties(data);
        update_md_buffers(data);
        md_gl_molecule_zero_velocity(&data->mold.gl_mol); // Do this explicitly to update the previous position to avoid motion blur trails
    }, data, data->tasks.prefetch_frames);
}

static bool load_dataset_from_file(ApplicationData* data, str_t path_to_file, md_molecule_api* mol_api, md_trajectory_api* traj_api, bool coarse_grained, bool deperiodize_on_load) {
    ASSERT(data);

    path_to_file = md_path_make_canonical(path_to_file, frame_allocator);
    if (path_to_file) {
        if (mol_api) {
            if (str_equal(data->files.molecule, path_to_file) && data->files.coarse_grained == coarse_grained) {
                // File already loaded as molecular data
                return true;
            }

            interrupt_async_tasks(data);
            free_molecule_data(data);
            free_trajectory_data(data);
            data->files.molecule = "";
            data->files.trajectory = "";

            if (!mol_api->init_from_file(&data->mold.mol, path_to_file, data->mold.mol_alloc)) {
                LOG_ERROR("Failed to load molecular data from file '%.*s'", path_to_file.len, path_to_file.ptr);
                return false;
            }
            LOG_SUCCESS("Successfully loaded molecular data from file '%.*s'", path_to_file.len, path_to_file.ptr);

            data->files.molecule = path_to_file;
            data->files.coarse_grained = coarse_grained;
            // @NOTE: If the dataset is coarse-grained, then postprocessing must be aware
            md_util_postprocess_flags_t flags = coarse_grained ? MD_UTIL_POSTPROCESS_COARSE_GRAINED : MD_UTIL_POSTPROCESS_ALL;
            md_util_postprocess_molecule(&data->mold.mol, data->mold.mol_alloc, flags);
            init_molecule_data(data);

            // @NOTE: Some files contain both atomic coordinates and trajectory
            if (traj_api) {
                LOG_INFO("File also contains trajectory, attempting to load trajectory");
            } else {
                return true;
            }
        }

        if (traj_api) {
            if (!data->mold.mol.atom.count) {
                LOG_ERROR("Before loading a trajectory, molecular data needs to be present");
                return false;
            }

            if (str_equal(data->files.trajectory, path_to_file) && data->files.deperiodize == deperiodize_on_load) {
                // Same as loaded file
                return true;
            }

            bool success = load_trajectory_data(data, path_to_file, deperiodize_on_load);
            if (success) {
                LOG_SUCCESS("Successfully opened trajectory from file '%.*s'", path_to_file.len, path_to_file.ptr);
                return true;
            } else {
                LOG_ERROR("Failed to opened trajectory from file '%.*s'", path_to_file.len, path_to_file.ptr);
            }
        }
    }

    return false;
}

// ### WORKSPACE ###
static RepresentationType get_rep_type(str_t str) {
    if (str_equal_cstr(str, "SPACE_FILL"))
        return RepresentationType::SpaceFill;
    else if (str_equal_cstr(str, "LICORICE"))
        return RepresentationType::Licorice;
    else if (str_equal_cstr(str, "BALL_AND_STICK"))    // Ball and stick is removed for now
        return RepresentationType::SpaceFill;
    else if (str_equal_cstr(str, "RIBBONS"))
        return RepresentationType::Ribbons;
    else if (str_equal_cstr(str, "CARTOON"))
        return RepresentationType::Cartoon;
    else
        return RepresentationType::SpaceFill;
}

static str_t get_rep_type_name(RepresentationType type) {
    switch (type) {
        case RepresentationType::SpaceFill:
            return STR("SPACE_FILL");
        case RepresentationType::Licorice:
            return STR("LICORICE");
        /*case RepresentationType::BallAndStick:
            return make_cstr("BALL_AND_STICK");*/
        case RepresentationType::Ribbons:
            return STR("RIBBONS");
        case RepresentationType::Cartoon:
            return STR("CARTOON");
        default:
            return STR("UNKNOWN");
    }
}

static ColorMapping get_color_mapping(str_t str) {
    if (str_equal_cstr(str, "UNIFORM"))
        return ColorMapping::Uniform;
    else if (str_equal_cstr(str, "CPK"))
        return ColorMapping::Cpk;
    else if (str_equal_cstr(str, "ATOM_LABEL"))
        return ColorMapping::AtomLabel;
    else if (str_equal_cstr(str, "ATOM_INDEX"))
        return ColorMapping::AtomIndex;
    else if (str_equal_cstr(str, "RES_ID"))
        return ColorMapping::ResId;
    else if (str_equal_cstr(str, "RES_INDEX"))
        return ColorMapping::ResIndex;
    else if (str_equal_cstr(str, "CHAIN_ID"))
        return ColorMapping::ChainId;
    else if (str_equal_cstr(str, "CHAIN_INDEX"))
        return ColorMapping::ChainIndex;
    else if (str_equal_cstr(str, "SECONDARY_STRUCTURE"))
        return ColorMapping::SecondaryStructure;
    else
        return ColorMapping::Cpk;
}

static str_t get_color_mapping_name(ColorMapping mapping) {
    switch (mapping) {
        case ColorMapping::Uniform:
            return STR("UNIFORM");
        case ColorMapping::Cpk:
            return STR("CPK");
        case ColorMapping::AtomLabel:
            return STR("ATOM_LABEL");
        case ColorMapping::AtomIndex:
            return STR("ATOM_INDEX");
        case ColorMapping::ResId:
            return STR("RES_ID");
        case ColorMapping::ResIndex:
            return STR("RES_INDEX");
        case ColorMapping::ChainId:
            return STR("CHAIN_ID");
        case ColorMapping::ChainIndex:
            return STR("CHAIN_INDEX");
        case ColorMapping::SecondaryStructure:
            return STR("SECONDARY_STRUCTURE");
        default:
            return STR("UNDEFINED");
    }
}

static vec4_t parse_vec4(str_t txt, vec4_t default_val = {1,1,1,1}) {
    vec4_t res = default_val;
    str_t tok = {};
    int i = 0;
    while (extract_token_delim(&tok, &txt, ',') && i < 4) {
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
    {"[Files]", "MoleculeFile",             SerializationType_Path,     offsetof(ApplicationData, files.molecule),     sizeof(ApplicationData::files.molecule)},
    {"[Files]", "TrajectoryFile",           SerializationType_Path,     offsetof(ApplicationData, files.trajectory),   sizeof(ApplicationData::files.trajectory)},
    {"[Files]", "CoarseGrained",            SerializationType_Bool,     offsetof(ApplicationData, files.coarse_grained)},
    {"[Files]", "Deperiodize",              SerializationType_Bool,     offsetof(ApplicationData, files.deperiodize)},
    
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
    {"[Representation]",        offsetof(ApplicationData, representation.reps),             sizeof(Representation),     serialize_create_rep},
    {"[AtomElementMapping]",    offsetof(ApplicationData, dataset.atom_element_remappings), sizeof(AtomElementMapping), serialize_create_atom_elem_mapping},
    {"[Selection]",             offsetof(ApplicationData, selection.stored_selections),     sizeof(Selection),          serialize_create_selection},
};

#define COMPARE(str, ref) (str_equal_cstr_n(str, ref"", sizeof(ref) - 1))
#define EXTRACT(str, ref) (str_equal_cstr_n(str, ref"", sizeof(ref) - 1) && (line = str_trim(str_substr(line, sizeof(ref) - 1))).len > 0)
#define EXTRACT_PARAM_LINE(line, txt) (c_txt.len && c_txt[0] != '[' && (str_extract_line(&line, &c_txt)))

static const SerializationObject* find_serialization_target(str_t group, str_t label) {
    for (size_t i = 0; i < ARRAY_SIZE(serialization_targets); ++i) {
        if (str_equal_cstr(group, serialization_targets[i].group) && str_equal_cstr(label, serialization_targets[i].label)) {
            return &serialization_targets[i];
        }
    }
    return NULL;
}

static const SerializationArray* find_serialization_array_group(str_t group) {
    for (size_t i = 0; i < ARRAY_SIZE(serialization_array_groups); ++i) {
        if (str_equal_cstr(group, serialization_array_groups[i].group)) {
            return &serialization_array_groups[i];
        }
    }
    return NULL;
}

static void deserialize_object(const SerializationObject* target, char* ptr, str_t* buf, str_t filename) {
    str_t line;
    if (str_extract_line(&line, buf)) {
        str_t arg = str_trim(line);

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
            md_strb_t path = md_strb_create(frame_allocator);
            path += extract_path_without_file(filename);
            path += arg;
            str_t can_path = md_path_make_canonical(path, frame_allocator);
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
            str_t token = STR("\"\"\"");
            if (str_equal_n(arg, token, token.len)) {
                // Roll back buf to arg + 3
                const char* beg = arg.ptr + token.len;
                buf->len = buf->end() - beg;
                buf->ptr = beg;
                const char* end = str_find_str(*buf, token).ptr;
                if (end) {
                    std::string str(beg, end - beg);
                    editor.SetText(str);
                    // Set buf pointer to after marker
                    const char* pos = end + token.len;
                    buf->len = buf->end() - pos;
                    buf->ptr = pos;
                    break;
                } else {
                    LOG_ERROR("Malformed end token for script");
                    return;
                }
            } else {
                LOG_ERROR("Malformed start token for script");
                return;
            }
        }
        case SerializationType_Bitfield:
        {
            // Bitfield starts with ###
            // and ends with ###
            str_t token = STR("###");
            if (str_equal_n(arg, token, token.len)) {
                // Roll back buf to arg + 3
                const char* beg = arg.ptr + token.len;
                buf->len = buf->end() - beg;
                buf->ptr = beg;
                const char* end = str_find_str(*buf, token).ptr;
                if (end) {
                    void* base64_data = md_alloc(frame_allocator, md_base64_decode_size_in_bytes(end-beg));
                    int64_t base64_size = md_base64_decode(base64_data, beg, end-beg);
                    md_bitfield_t* bf = (md_bitfield_t*)(ptr + target->struct_byte_offset);
                    if (!base64_size || !md_bitfield_deserialize(bf, base64_data, base64_size)) {
                        LOG_ERROR("Failed to deserialize bitfield");
                        md_bitfield_clear(bf);
                        return;
                    }
                    // Set buf pointer to after marker
                    const char* pos = end + token.len;
                    buf->len = buf->end() - pos;
                    buf->ptr = pos;
                    break;
                } else {
                    LOG_ERROR("Malformed end token for bitfield");
                    return;
                }
            } else {
                LOG_ERROR("Malformed start token for bitfield");
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
    defer { str_free(txt, frame_allocator); };

    if (!txt.len) {
        LOG_ERROR("Could not open workspace file: '%.*s", (int)filename.len, filename.ptr);
        return;
    }

    // Reset and clear things
    clear_representations(data);
    editor.SetText("");

    data->animation = {};
    reset_view(data, false, true);

    str_t group = {};
    str_t c_txt = txt;
    str_t line = {};

    str_t cur_molecule_file = str_copy(data->files.molecule, frame_allocator);
    str_t cur_trajectory_file = str_copy(data->files.trajectory, frame_allocator);
    bool  cur_coarse_grained = data->files.coarse_grained;
    bool  cur_deperiodize = data->files.deperiodize;

    const SerializationArray* arr_group = NULL;
    void* ptr = 0;

    while (str_extract_line(&line, &c_txt)) {
        line = str_trim(line);
        if (line[0] == '[') {
            group = line;
            arr_group = find_serialization_array_group(group);
            if (arr_group) {
                ptr = arr_group->create_item_func(data);
            } else {
                ptr = data;
            }
        } else {
            int64_t loc = str_find_char(line, '=');
            if (loc != -1) {
                str_t label = str_trim(str_substr(line, 0, loc));
                const SerializationObject* target = find_serialization_target(group, label);
                if (target) {
                    // Move read pointer back 'loc+1' => after '='.
                    const char* pos = line.ptr + loc + 1;
                    c_txt.len = c_txt.end() - pos;
                    c_txt.ptr = pos;
                    deserialize_object(target, (char*)ptr, &c_txt, filename);
                } else {
                    LOG_ERROR("Could not recognize serialization target '%.*s' in group '%.*s,", (int)label.len, label.ptr, (int)group.len, group.ptr);
                }
            }
        }
    }
    data->view.animation.target_position    = data->view.camera.position;
    data->view.animation.target_orientation = data->view.camera.orientation;
    data->view.animation.target_distance    = data->view.camera.focus_distance;
    
    data->files.workspace = filename;

    str_t new_molecule_file   = str_copy(data->files.molecule, frame_allocator);
    str_t new_trajectory_file = str_copy(data->files.trajectory, frame_allocator);
    bool  new_coarse_grained  = data->files.coarse_grained;
    bool  new_deperiodize     = data->files.deperiodize;

    // When we de-serialize we overwrite the two following paths, even though they are not loaded.
    // So we copy them back to their original values.
    data->files.molecule = cur_molecule_file;
    data->files.trajectory = cur_trajectory_file;
    data->files.coarse_grained = cur_coarse_grained;
    data->files.deperiodize = cur_deperiodize;

    md_molecule_api* mol_api = load::mol::get_api(new_molecule_file);
    md_trajectory_api* traj_api = load::traj::get_api(new_trajectory_file);

    if (new_molecule_file && load_dataset_from_file(data, new_molecule_file, mol_api, nullptr, new_coarse_grained, new_deperiodize)) {
        init_all_representations(data);
        update_all_representations(data);
    }

    if (new_trajectory_file) {
        load_dataset_from_file(data, new_trajectory_file, nullptr, traj_api);
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
        str_t rel_path = md_path_make_relative(filename, {str, len}, frame_allocator);
        if (rel_path.ptr && rel_path.len) {
            fprintf(file, "%.*s\n", (int)rel_path.len, rel_path.ptr);
        }     
        break;
    }
    case SerializationType_Script:
    {
        std::string str = editor.GetText();
        fprintf(file, "\"\"\"%s\"\"\"\n", str.c_str());
        break;
    }
    case SerializationType_Bitfield:
    {
        const md_bitfield_t* bf = (const md_bitfield_t*)((const char*)ptr + target.struct_byte_offset);
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
        LOG_ERROR("Could not open workspace file: '%.*s", (int)filename.len, filename.ptr);
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
    // @NOTE(Robin): Package this as a main(render thread) task to ensure that it is done when it has been 
    image_t img = {0};
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

    str_t ext = extract_ext(data->screenshot.path_to_file);
    if (str_equal_cstr_ignore_case(ext, "jpg")) {
        const int quality = 95;
        write_image_jpg(img, data->screenshot.path_to_file, quality);
    } else if (str_equal_cstr_ignore_case(ext, "png")) {
        write_image_png(img, data->screenshot.path_to_file);
    } else if (str_equal_cstr_ignore_case(ext, "bmp")) {
        write_image_bmp(img, data->screenshot.path_to_file);
    } else {
        LOG_ERROR("Could not match extension when saving screenshot");
    }
}

// #representation
static Representation* create_representation(ApplicationData* data, RepresentationType type, ColorMapping color_mapping, str_t filter) {
    ASSERT(data);
    Representation rep = {};
    rep.type = type;
    rep.color_mapping = color_mapping;
    if (!str_empty(filter)) {
        strncpy(rep.filt, filter.ptr, MIN((size_t)filter.len, sizeof(rep.filt)));
    }
    init_representation(data, &rep);
    update_representation(data, &rep);
    return md_array_push(data->representation.reps, rep, persistent_allocator);
}

static Representation* clone_representation(ApplicationData* data, const Representation& rep) {
    ASSERT(data);
    Representation* clone = md_array_push(data->representation.reps, rep, persistent_allocator);
    clone->md_rep = {0};
    clone->atom_mask = {0};
    init_representation(data, clone);
    update_representation(data, clone);
    return clone;
}

static void remove_representation(ApplicationData* data, int idx) {
    ASSERT(data);
    ASSERT(idx < md_array_size(data->representation.reps));
    auto& rep = data->representation.reps[idx];
    md_bitfield_free(&rep.atom_mask);
    md_gl_representation_free(&rep.md_rep);
    data->representation.reps[idx] = *md_array_last(data->representation.reps);
    md_array_pop(data->representation.reps);
}

static void recompute_atom_visibility_mask(ApplicationData* data) {
    ASSERT(data);

    auto& mask = data->representation.atom_visibility_mask;

    md_bitfield_clear(&mask);
    for (int64_t i = 0; i < md_array_size(data->representation.reps); ++i) {
        auto& rep = data->representation.reps[i];
        if (!rep.enabled) continue;
        md_bitfield_or_inplace(&mask, &rep.atom_mask);
    }

    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

static void update_all_representations(ApplicationData* data) {
    for (int64_t i = 0; i < md_array_size(data->representation.reps); ++i) {
        auto& rep = data->representation.reps[i];
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
        case ColorMapping::AtomLabel:
            color_atoms_label(colors, mol.atom.count, mol);
            break;
        case ColorMapping::AtomIndex:
            color_atoms_idx(colors, mol.atom.count, mol);
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
                memset(colors, 0xFFFFFFFF, bytes);
                const float* values = rep->prop->data.values;
                if (rep->prop->data.aggregate) {
                    const int dim = rep->prop->data.dim[0];
                    md_script_visualization_t vis = {0};
                    bool result = false;
                    
                    if (md_semaphore_aquire(&data->mold.script.ir_semaphore)) {
                        defer { md_semaphore_release(&data->mold.script.ir_semaphore); };
                        
                        if (md_script_ir_valid(data->mold.script.ir)) {
                            md_script_visualization_args_t args = {
                                .payload = rep->prop->vis_payload,
                                .ir = data->mold.script.ir,
                                .mol = &data->mold.mol,
                                .traj = NULL,
                                .alloc = frame_allocator,
                                .flags = MD_SCRIPT_VISUALIZE_ATOMS
                            };
                            result = md_script_visualization_init(&vis, args);
                        }
                    }
                    if (result) {
                        if (dim == (int)md_array_size(vis.structures.atom_masks)) {
                            int i0 = CLAMP((int)data->animation.frame + 0, 0, rep->prop->data.num_values / dim - 1);
                            int i1 = CLAMP((int)data->animation.frame + 1, 0, rep->prop->data.num_values / dim - 1);
                            float frame_fract = fractf((float)data->animation.frame);

                            md_bitfield_t mask = {0};
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

    switch (rep->type) {
    case RepresentationType::SpaceFill:
    case RepresentationType::Licorice:
        rep->type_is_valid = true;
        break;
    case RepresentationType::Ribbons:
    case RepresentationType::Cartoon:
        rep->type_is_valid = mol.backbone.range_count > 0;
        break;
    default:
        ASSERT(false);
    }

    if (rep->dynamic_evaluation) {
        rep->filt_is_dirty = true;
    }

    if (rep->filt_is_dirty) {
        rep->filt_is_valid = filter_expression(data, str_from_cstr(rep->filt), &rep->atom_mask, &rep->filt_is_dynamic, rep->filt_error, sizeof(rep->filt_error));
        rep->filt_is_dirty = false;
    }

    if (rep->filt_is_valid) {
        filter_colors(colors, mol.atom.count, &rep->atom_mask);
        data->representation.atom_visibility_mask_dirty = true;

        md_gl_representation_type_t type = MD_GL_REP_SPACE_FILL;
        md_gl_representation_args_t args = {};
        switch (rep->type) {
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

#if EXPERIMENTAL_GFX_API
        md_gfx_rep_attr_t attributes = {};
        attributes.spacefill.radius_scale = 1.0f;
        md_gfx_rep_set_type_and_attr(rep->gfx_rep, MD_GFX_REP_TYPE_SPACEFILL, &attributes);
        md_gfx_rep_set_color(rep->gfx_rep, 0, (uint32_t)mol.atom.count, (md_gfx_color_t*)colors, 0);
#endif
    }
}

static void init_representation(ApplicationData* data, Representation* rep) {
#if EXPERIMENTAL_GFX_API
    rep->gfx_rep = md_gfx_rep_create(data->mold.mol.atom.count);
#endif
    md_gl_representation_init(&rep->md_rep, &data->mold.gl_mol);
    md_bitfield_init(&rep->atom_mask, persistent_allocator);
    rep->filt_is_dirty = true;
}

static void init_all_representations(ApplicationData* data) {
    for (int64_t i = 0; i < md_array_size(data->representation.reps); ++i) {
        auto& rep = data->representation.reps[i];
        init_representation(data, &rep);
    }
}

static void clear_representations(ApplicationData* data) {
    ASSERT(data);
    while (md_array_size(data->representation.reps) > 0) {
        remove_representation(data, (int32_t)md_array_size(data->representation.reps) - 1);
    }
}

// #selection
static Selection* create_selection(ApplicationData* data, str_t name, md_bitfield_t* atom_mask) {
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
        LOG_ERROR("Index [%i] out of range when trying to remove selection", idx);
    }
    auto item = &data->selection.stored_selections[idx];
    md_bitfield_free(&item->atom_mask);
    
    data->selection.stored_selections[idx] = *md_array_last(data->selection.stored_selections);
    md_array_pop(data->selection.stored_selections);
}

// #camera-control
static void handle_camera_interaction(ApplicationData* data) {
    ASSERT(data);

    enum class RegionMode { Append, Remove };

    ImGui::BeginCanvas("Main interarction window", true);
    bool pressed = ImGui::InvisibleButton("canvas", ImGui::GetContentRegionAvail(), ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);

    if (pressed || ImGui::IsItemActive()) {
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

            ImGuiWindow* window = ImGui::GetCurrentWindow();
            ASSERT(window);
            ImDrawList* dl = window->DrawList;
            dl->AddRectFilled(pos, pos + ext, fill_col);
            dl->AddRect(pos, pos + ext, line_col);

            const ImVec2 min_p = ImMin(pos, pos + ext) - window->Pos;
            const ImVec2 max_p = ImMax(pos, pos + ext) - window->Pos;

            md_bitfield_t mask = { 0 };
            md_bitfield_init(&mask, frame_allocator);

            if (min_p != max_p) {
                md_bitfield_clear(&data->selection.current_highlight_mask);
                data->mold.dirty_buffers |= MolBit_DirtyFlags;
                data->selection.selecting = true;

                const vec2_t res = { (float)data->ctx.window.width, (float)data->ctx.window.height };
                const mat4_t mvp = data->view.param.matrix.current.proj * data->view.param.matrix.current.view;
                const md_bitfield_t* vis_mask = &data->representation.atom_visibility_mask;

                int64_t beg_bit = vis_mask->beg_bit;
                int64_t end_bit = vis_mask->end_bit;
                while ((beg_bit = md_bitfield_scan(vis_mask, beg_bit, end_bit)) != 0) {
                    int64_t i = beg_bit - 1;
                    vec4_t p = { data->mold.mol.atom.x[i], data->mold.mol.atom.y[i], data->mold.mol.atom.z[i], 1.0f };
                    p = mat4_mul_vec4(mvp, p);

                    vec2_t c = {
                        (p.x / p.w * 0.5f + 0.5f) * res.x,
                        (-p.y / p.w * 0.5f + 0.5f) * res.y
                    };

                    if (min_p.x <= c.x && c.x <= max_p.x && min_p.y <= c.y && c.y <= max_p.y) {
                        md_bitfield_set_bit(&mask, i);
                    }
                }
                grow_mask_by_current_selection_granularity(&mask, *data);

                if (mode == RegionMode::Append) {
                    md_bitfield_or(&data->selection.current_highlight_mask, &data->selection.current_selection_mask, &mask);
                }
                else if (mode == RegionMode::Remove) {
                    md_bitfield_andnot(&data->selection.current_highlight_mask, &data->selection.current_selection_mask, &mask);
                }
                if (pressed) {
                    md_bitfield_copy(&data->selection.current_selection_mask, &data->selection.current_highlight_mask);
                }
            }
            else if (pressed) {
                if (data->picking.idx != INVALID_PICKING_IDX) {
                    if (mode == RegionMode::Append) {
                        single_selection_sequence_push_idx(&data->selection.single_selection_sequence, data->picking.idx);
                        md_bitfield_set_bit(&mask, data->picking.idx);
                        grow_mask_by_current_selection_granularity(&mask, *data);
                        md_bitfield_or_inplace(&data->selection.current_selection_mask, &mask);
                    }
                    else if (mode == RegionMode::Remove) {
                        single_selection_sequence_pop_idx(&data->selection.single_selection_sequence, data->picking.idx);
                        md_bitfield_set_bit(&mask, data->picking.idx);
                        grow_mask_by_current_selection_granularity(&mask, *data);
                        md_bitfield_andnot_inplace(&data->selection.current_selection_mask, &mask);
                    }
                }
                else {
                    single_selection_sequence_clear(&data->selection.single_selection_sequence);
                    md_bitfield_clear(&data->selection.current_selection_mask);
                    md_bitfield_clear(&data->selection.current_highlight_mask);
                }
            }

            data->mold.dirty_buffers |= MolBit_DirtyFlags;
        }
    }
    else if (ImGui::IsItemHovered()) {
        if (data->picking.idx != INVALID_PICKING_IDX && data->picking.idx <= data->mold.mol.atom.count) {
            md_bitfield_clear(&data->selection.current_highlight_mask);
            md_bitfield_set_bit(&data->selection.current_highlight_mask, data->picking.idx);
            grow_mask_by_current_selection_granularity(&data->selection.current_highlight_mask, *data);
            data->mold.dirty_buffers |= MolBit_DirtyFlags;
            draw_atom_info_window(*data, data->picking.idx);
        }
    }

    bool open_atom_context = false;
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
            input.screen_size = {(float)data->ctx.window.width, (float)data->ctx.window.height};
            input.dolly_delta = scroll_delta;
            input.fov_y = data->view.camera.fov_y;

            vec3_t pos = data->view.animation.target_position;
            quat_t ori = data->view.animation.target_orientation;
            float dist = data->view.animation.target_distance;
            if (camera_controller_trackball(&pos, &ori, &dist, input, data->view.trackball_param, 0)) {
                // We could make the camera interaction more snappy, by directly modifying camera[pos, ori, dist] here
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
                    reset_view(data, true, true);
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

    ImGui::EndCanvas();

    if (open_atom_context) {
        ImGui::OpenPopup("AtomContextPopup");
    }
}

static void handle_camera_animation(ApplicationData* data) {
    // We use an exponential interpolation of the deltas with a common factor
    const float TARGET_DT = 1.0f / 100.0f;
    const float TARGET_FACTOR = 0.1f;
    const float dt = CLAMP(data->ctx.timing.delta_s, 1.0f / 1000.f, 1.0f / 20.f);
    const float interpolation_factor = TARGET_FACTOR * (dt / TARGET_DT);

    vec3_t& current_pos = data->view.camera.position;
    vec3_t& target_pos  = data->view.animation.target_position;

    quat_t& current_ori = data->view.camera.orientation;
    quat_t& target_ori  = data->view.animation.target_orientation;

    float&  current_dist = data->view.camera.focus_distance;
    float&  target_dist  = data->view.animation.target_distance;

    // We want to interpolate along an arc which is formed by maintaining a distance to the look_at position and smoothly interpolating the orientation,
    // This means that 
    // We linearly interpolate a look_at position which is implicitly defined by position, orientation and distance
    // There is some precision errors creeping into the posision because we transform back and forth to look at using the orientation

    current_dist = lerp(current_dist, target_dist, interpolation_factor);

    const vec3_t current_look_at = current_pos - current_ori * vec3_t{0, 0, current_dist};
    const vec3_t target_look_at  = target_pos  -  target_ori * vec3_t{0, 0,  target_dist};
    const vec3_t look_at = lerp(current_look_at, target_look_at, interpolation_factor);

    current_ori = quat_normalize(quat_slerp(current_ori, target_ori, interpolation_factor));
    current_pos = look_at + current_ori * vec3_t{0, 0, current_dist};
        
    data->visuals.dof.focus_depth = current_dist;
#if 0
    ImGui::Begin("Camera Debug Info");
    ImGui::Text("pos cur [%.4f %.4f %.4f]", current_pos.x, current_pos.y, current_pos.z);
    ImGui::Text("pos tar [%.4f %.4f %.4f]", target_pos.x, target_pos.y, target_pos.z);
    ImGui::Text("ori cur [%.4f %.4f %.4f %.4f]", current_ori.x, current_ori.y, current_ori.z, target_ori.w);
    ImGui::Text("ori tar [%.4f %.4f %.4f %.4f]", target_ori.x, target_ori.y, target_ori.z, target_ori.w);
    ImGui::Text("dis cur [%.4f]", current_dist);
    ImGui::Text("dis tar [%.4f]", target_dist);
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
                                             data->mold.mol.atom.radius, data->representation.atom_visibility_mask);
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
        glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
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
        GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_POST_TONEMAP };

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    // Enable all draw buffers
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, data->gbuffer.deferred.fbo);
    glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);

    PUSH_GPU_SECTION("G-Buffer fill")

    // Immediate mode graphics

    if (data->simulation_box.enabled && data->simulation_box.box != mat3_t{0}) {
        PUSH_GPU_SECTION("Draw Simulation Box")
        const mat3_t box = data->simulation_box.box;
        const vec3_t min_box = box * vec3_t{0, 0, 0};
        const vec3_t max_box = box * vec3_t{1, 1, 1};

        immediate::set_model_view_matrix(data->view.param.matrix.current.view);
        immediate::set_proj_matrix(data->view.param.matrix.current.proj_jittered);
        immediate::draw_box_wireframe(min_box, max_box, convert_color(data->simulation_box.color));
        immediate::render();
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

if (!use_gfx) {
    // DRAW VELOCITY OF STATIC OBJECTS
    PUSH_GPU_SECTION("Blit Static Velocity")
    glDrawBuffer(GL_COLOR_ATTACHMENT_VELOCITY);
    glDepthMask(0);
    postprocessing::blit_static_velocity(data->gbuffer.deferred.depth, data->view.param);
    glDepthMask(1);
    POP_GPU_SECTION()
}
    glDepthMask(1);
    glColorMask(1, 1, 1, 1);

    // DRAW REPRESENTATIONS
    PUSH_GPU_SECTION("Representation")
    glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
    draw_representations(data);
    POP_GPU_SECTION()

    glDrawBuffer(GL_COLOR_ATTACHMENT_POST_TONEMAP);  // Post_Tonemap buffer

if (!use_gfx) {
    PUSH_GPU_SECTION("Selection")
    const bool atom_selection_empty = md_bitfield_popcount(&data->selection.current_selection_mask) == 0;
    const bool atom_highlight_empty = md_bitfield_popcount(&data->selection.current_highlight_mask) == 0;

    glDepthMask(0);
    glColorMask(0, 0, 0, 0);


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
        postprocessing::blit_color({0, 0, 1, 0.25f});

        glStencilFunc(GL_EQUAL, 2, 3);
        postprocessing::blit_color({0, 0, 0.25f, 0.4f});
    }

    if (!atom_highlight_empty) {
        glStencilFunc(GL_EQUAL, 4, 4);
        postprocessing::blit_color({1, 1, 0, 0.25f});
    }

    glDisable(GL_STENCIL_TEST);
    if (!atom_selection_empty) {
        PUSH_GPU_SECTION("Desaturate")
        const float saturation = data->selection.color.selection_saturation;
        glDrawBuffer(GL_COLOR_ATTACHMENT_COLOR);
        postprocessing::scale_hsv(data->gbuffer.deferred.color, vec3_t{1, saturation, 1});
        POP_GPU_SECTION()
    }

    glDisable(GL_STENCIL_TEST);
    glDepthFunc(GL_LESS);
    glDepthMask(0);
    POP_GPU_SECTION()
}

    PUSH_GPU_SECTION("Draw Visualization Geometry")
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    immediate::set_model_view_matrix(data->view.param.matrix.current.view);
    immediate::set_proj_matrix(data->view.param.matrix.current.proj);

    const md_script_visualization_t& vis = data->mold.script.vis;
    const vec3_t* vertices = (const vec3_t*)vis.vertex.pos;

    const uint32_t point_color      = convert_color(data->script.point_color);
    const uint32_t line_color       = convert_color(data->script.line_color);
    const uint32_t triangle_color   = convert_color(data->script.triangle_color);

    for (int64_t i = 0; i < vis.triangle.count; ++i) {
        ASSERT(vis.triangle.idx);
        uint32_t idx[3] = { vis.triangle.idx[i * 3 + 0], vis.triangle.idx[i * 3 + 1], vis.triangle.idx[i * 3 + 2] };
        immediate::draw_triangle(vertices[idx[0]], vertices[idx[1]], vertices[idx[2]], triangle_color);
    }

    for (int64_t i = 0; i < vis.line.count; ++i) {
        ASSERT(vis.line.idx);
        uint32_t idx[2] = { vis.line.idx[i * 2 + 0], vis.line.idx[i * 2 + 1] };
        immediate::draw_line(vertices[idx[0]], vertices[idx[1]], line_color);
    }

    for (int64_t i = 0; i < vis.point.count; ++i) {
        ASSERT(vis.point.idx);
        uint32_t idx = vis.point.idx[i];
        immediate::draw_point(vertices[idx], point_color);
    }
    immediate::render();

    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    POP_GPU_SECTION()

    POP_GPU_SECTION()  // G-buffer
}

static void handle_picking(ApplicationData* data) {
    PUSH_CPU_SECTION("PICKING") {
        vec2_t mouse_pos = vec_cast(ImGui::GetMousePos() - ImGui::GetMainViewport()->Pos);
        vec2_t coord = {mouse_pos.x, (float)data->gbuffer.height - mouse_pos.y};
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
            data->picking = read_picking_data(&data->gbuffer, (int)coord.x, (int)coord.y);
            const vec4_t viewport = {0, 0, (float)data->gbuffer.width, (float)data->gbuffer.height};
            const mat4_t inv_VP = data->view.param.matrix.inverse.view * data->view.param.matrix.inverse.proj_jittered;
            data->picking.world_coord = mat4_unproject({coord.x, coord.y, data->picking.depth}, inv_VP, viewport);
#endif
        }
        data->selection.hovered = -1;
        if (data->picking.idx != INVALID_PICKING_IDX) {
            data->selection.hovered = data->picking.idx;
        }
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
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
    desc.depth_of_field.focus_depth = data.visuals.dof.focus_depth;
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
                md_array_push(transforms, M, persistent_allocator);
            }
        }

        md_gfx_draw_op_t* draw_ops = 0;
        for (int64_t i = 0; i < md_array_size(data->representation.reps); ++i) {
            if (data->representation.reps[i].enabled) {
                md_gfx_draw_op_t op;
                op.structure = data->mold.gfx_structure;
                op.representation = data->representation.reps[i].gfx_rep;
                op.model_mat = NULL;
                md_array_push(draw_ops, op, frame_allocator);
                
                for (uint32_t j = 0; j < instance_count; ++j) {
                    md_gfx_draw_op_t op;
                    op.structure = data->mold.gfx_structure;
                    op.representation = data->representation.reps[i].gfx_rep;
                    op.model_mat = &transforms[j];
                    md_array_push(draw_ops, op, frame_allocator);
                }
                
            }
        }

        md_gfx_draw((uint32_t)md_array_size(draw_ops), draw_ops, &data->view.param.matrix.current.proj, &data->view.param.matrix.current.view, &data->view.param.matrix.inverse.proj, &data->view.param.matrix.inverse.view);
    } else {
#endif
        md_gl_draw_op_t* draw_ops = 0;
        for (int64_t i = 0; i < md_array_size(data->representation.reps); ++i) {
            if (data->representation.reps[i].enabled && data->representation.reps[i].type_is_valid) {
                md_gl_draw_op_t op = {
                    &data->representation.reps[i].md_rep,
                    NULL,
                };
                md_array_push(draw_ops, op, frame_allocator);
            }
        }

        md_gl_draw_args_t args = {
            .shaders = &data->mold.gl_shaders,
            .draw_operations = {
                .count = (uint32_t)md_array_size(draw_ops),
                .ops = draw_ops,
        },
        .view_transform = {
                .view_matrix = &data->view.param.matrix.current.view.elem[0][0],
                .projection_matrix = &data->view.param.matrix.current.proj_jittered.elem[0][0],
                // These two are for temporal anti-aliasing reprojection (optional)
                .prev_view_matrix = &data->view.param.matrix.previous.view.elem[0][0],
                .prev_projection_matrix = &data->view.param.matrix.previous.proj_jittered.elem[0][0],
        },
        };

        md_gl_draw(&args);
#if EXPERIMENTAL_GFX_API
    }
#endif
}

static void draw_representations_lean_and_mean(ApplicationData* data, uint32_t mask) {
    md_gl_draw_op_t* draw_ops = 0;
    for (int64_t i = 0; i < md_array_size(data->representation.reps); ++i) {
        if (data->representation.reps[i].enabled && data->representation.reps[i].type_is_valid) {
            md_gl_draw_op_t op = { 
                &data->representation.reps[i].md_rep,
                NULL,
            };
            md_array_push(draw_ops, op, frame_allocator);
        }
    }

    md_gl_draw_args_t args = {
        .shaders = &data->mold.gl_shaders_lean_and_mean,
        .draw_operations = {
            .count = (uint32_t)md_array_size(draw_ops),
            .ops = draw_ops,
        },
        .view_transform = {
            .view_matrix = &data->view.param.matrix.current.view.elem[0][0],
            .projection_matrix = &data->view.param.matrix.current.proj_jittered.elem[0][0],
            // These two are for temporal anti-aliasing reprojection
            //.prev_model_view_matrix = &data->view.param.matrix.previous.view[0][0],
            //.prev_projection_matrix = &data->view.param.matrix.previous.proj_jittered[0][0],
        },
        .atom_mask = mask,
    };

    md_gl_draw(&args);
}
