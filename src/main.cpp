//#include <core/spatial_hash.h>

#include <core/md_compiler.h>

#define _CRT_SECURE_NO_WARNINGS
#if MD_COMPILER_MSVC
#pragma warning( disable : 26812 4244 )
#endif

#include <md_util.h>
#include <md_gl.h>
#include <md_filter.h>
#include <md_script.h>
#include <md_molecule.h>
#include <md_trajectory.h>
#include <md_frame_cache.h>

#include <core/md_sync.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_stack_allocator.h>
#include <core/md_pool_allocator.h>
#include <core/md_log.h>
#include <core/md_simd.h>
#include <core/md_file.h>
#include <core/md_array.inl>
#include <core/md_os.h>

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
//#include "gfx/molecule_draw.h"
#include "gfx/immediate_draw_utils.h"
#include "gfx/postprocessing_utils.h"
#include "gfx/volumerender_utils.h"
#include "gfx/conetracing_utils.h"

#include "string_util.h"
#include "random_util.h"
//#include "volume.h"
#include "imgui_widgets.h"
#include "implot_widgets.h"
//#include "plot_extended.h"
#include "application/application.h"
#include "console.h"
//#include "ramachandran.h"
#include "color_utils.h"
#include "isosurface.h"
#include "task_system.h"
#include "loader.h"

//#include "spatial_hash.h"

#include <atomic_queue.h>
#include <stdio.h>
//#include <filesystem>

#define PICKING_JITTER_HACK 0
#define DEPERIODIZE_ON_LOAD 1
#define SHOW_IMGUI_DEMO_WINDOW 0
#define VIAMD_RELEASE 0
#define EXPERIMENTAL_CONE_TRACED_AO 0
#define USE_MOLD 1
#define COMPILATION_TIME_DELAY_IN_SECONDS 1.5

#ifdef OS_MAC_OSX
const Key::Key_t KEY_CONSOLE = Key::KEY_WORLD_1;
#else  // WIN32 and Linux
// @TODO: Make sure this is currect for Linux?
const Key::Key_t KEY_CONSOLE = Key::KEY_GRAVE_ACCENT;
#endif

const Key::Key_t KEY_PLAY_PAUSE = Key::KEY_SPACE;
const Key::Key_t KEY_SKIP_TO_PREV_FRAME = Key::KEY_LEFT;
const Key::Key_t KEY_SKIP_TO_NEXT_FRAME = Key::KEY_RIGHT;
const Key::Key_t KEY_TOGGLE_SCREENSHOT_MODE = Key::KEY_F10;

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

#ifdef VIAMD_RELEASE
constexpr const char* CAFFINE_PDB = R"(
ATOM      1  N1  BENZ    1       5.040   1.944  -8.324                          
ATOM      2  C2  BENZ    1       6.469   2.092  -7.915                          
ATOM      3  C3  BENZ    1       7.431   0.865  -8.072                          
ATOM      4  C4  BENZ    1       6.916  -0.391  -8.544                          
ATOM      5  N5  BENZ    1       5.532  -0.541  -8.901                          
ATOM      6  C6  BENZ    1       4.590   0.523  -8.394                          
ATOM      7  C11 BENZ    1       4.045   3.041  -8.005                          
ATOM      8  H111BENZ    1       4.453   4.038  -8.264                          
ATOM      9  H112BENZ    1       3.101   2.907  -8.570                          
ATOM     10  H113BENZ    1       3.795   3.050  -6.926                          
ATOM     11  O21 BENZ    1       6.879   3.181  -7.503
ATOM     12  C51 BENZ    1       4.907  -1.659  -9.696
ATOM     13  H511BENZ    1       4.397  -1.273 -10.599                          
ATOM     14  H512BENZ    1       5.669  -2.391 -10.028                          
ATOM     15  H513BENZ    1       4.161  -2.209  -9.089
ATOM     16  O61 BENZ    1       3.470   0.208  -7.986                          
ATOM     17  N1  NSP3    1B      8.807   0.809  -7.799
ATOM     18  N1  NSP3    1C      7.982  -1.285  -8.604
ATOM     19  C1  CSP3    1D      9.015  -0.500  -8.152                          
ATOM     20  H1  CSP3    1D     10.007  -0.926  -8.079                          
ATOM     21  C1  CSP3    1E      9.756   1.835  -7.299                          
ATOM     22  H11 CSP3    1E     10.776   1.419  -7.199                          
ATOM     23  H12 CSP3    1E      9.437   2.207  -6.309                          
ATOM     24  H13 CSP3    1E      9.801   2.693  -7.994
)";
#endif

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
enum class RepresentationType { Vdw, Licorice, Ribbons, Cartoon };
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
    MolBit_DirtyPosition            = 0x1,
    MolBit_DirtySecondaryStructure  = 0x2,
    MolBit_DirtyFlags               = 0x4
};

enum RepBit_ {
    RepBit_DirtyColor   = 0x1,
    RepBit_DirtyFilter  = 0x2
};

// #struct Structure Declarations

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
    StrBuf<64> name = "rep";
    StrBuf<256> filter = "all";
    RepresentationType type = RepresentationType::Vdw;
    ColorMapping color_mapping = ColorMapping::Cpk;
    md_exp_bitfield_t atom_mask{};
    md_gl_representation_t md_rep{};

    bool enabled = true;
    bool show_in_selection = true;
    bool filter_is_valid = true;
    bool filter_is_dynamic = false;

    uint32_t flags = 0;

    // User defined color used in uniform mode
    vec4_t uniform_color = vec4_t{1,1,1,1};

    // VDW and Ball & Stick
    float radius = 1.f;

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

struct ApplicationData {
    // --- APPLICATION ---
    application::Context ctx {};

    // --- FILES ---
    // for keeping track of open files
    struct {
        StrBuf<256> molecule{};
        StrBuf<256> trajectory{};
        StrBuf<256> workspace{};
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
        md_gl_context_t     gl_ctx = {};
        md_gl_molecule_t    gl_mol = {};
        md_molecule_t       mol = {};
        md_trajectory_i     traj = {};
        md_frame_cache_t    frame_cache = {};

        struct {
            md_script_ir_t    ir = {};
            md_script_eval_t  full_eval = {};
            md_script_eval_t  filt_eval = {};

            md_exp_bitfield_t frame_mask = {};

            md_semaphore_t semaphore = {};

            bool compile_ir = false;
            bool evaluate_full = false;
            bool evaluate_filt = false;
            double time_since_last_change = 0.0;
        } script;
        uint32_t dirty_buffers = {0};
    } mold;

    // --- ASYNC TASKS HANDLES ---
    struct {
        task_system::ID prefetch_trajectory = task_system::INVALID_ID;
        task_system::ID evaluate_full = task_system::INVALID_ID;
        task_system::ID evaluate_filt = task_system::INVALID_ID;
    } tasks;

    // --- ATOM SELECTION ---
    struct {
        // bool show_window = false;
        SelectionLevel level_mode = SelectionLevel::Atom;
        SelectionOperator op_mode = SelectionOperator::Or;
        SelectionGrowth grow_mode = SelectionGrowth::CovalentBond;

        int32_t hovered = -1;
        int32_t right_clicked = -1;

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
    } selection;

    // --- STATISTICS ---
    struct {
        bool show_timeline_window = false;
        bool show_distribution_window = false;
        bool show_volume_window = false;

    } statistics;

    // --- FRAMEBUFFER ---
    GBuffer gbuffer {};

    PickingData picking {};

    // --- ANIMATION ---
    struct {
        double time = 0.f;  // double precision for long trajectories
        int32_t frame = 0;
        float fps = 10.f;
        InterpolationMode interpolation = InterpolationMode::Cubic;
        PlaybackMode mode = PlaybackMode::Stopped;
        bool apply_pbc = false;
    } animation;

    // --- TIMELINE---
    struct {
        struct {
            bool enabled = true;
            double min = 0;
            double max = 0;
            
            bool temporal_window = false;
            float window_extent = 10.f;
        } filter;

        struct {
            double min = 0;
            double max = 1;
        } view_range;
    } timeline;

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
        /*
        struct {
            int dim_x = 0;
            int dim_y = 0;
            GLuint gbuffer = 0;
            GLuint depth_tex = 0;
            GLuint color_tex = 0;
            GLuint normal_tex = 0;
            GLuint picking_tex = 0;
        } render_target;
        */

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

        Camera camera = {};

    } density_volume;

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

    // --- REPRESENTATIONS ---
    struct {
        Representation* buffer = {};
        md_exp_bitfield_t atom_visibility_mask{};
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
};

struct PropertyItem {
    StrBuf<32> lbl;
    uint32_t col;
    int idx;
    bool show;
};

//static void postprocess_frame(md_frame_data_t* frame, void* user_data);

static void interpolate_atomic_properties(ApplicationData* data);
static void update_view_param(ApplicationData* data);
static void reset_view(ApplicationData* data, bool move_camera = false, bool smooth_transition = false);
static float compute_avg_ms(float dt);

static void compute_aabb(vec3_t* aabb_min, vec3_t* aabb_max, const float* x, const float* y, const float* z, int64_t count);

static bool handle_selection(ApplicationData* data);
// static void handle_animation(ApplicationData* data);
static void handle_camera_interaction(ApplicationData* data);
static void handle_camera_animation(ApplicationData* data);

//static void update_properties(ApplicationData* data);
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
//static void draw_shape_space_window(ApplicationData* data);
static void draw_volume_window(ApplicationData* data);
static void draw_property_editor_window(ApplicationData* data);
static void draw_dataset_window(ApplicationData* data);
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
static Representation* create_representation(ApplicationData* data, RepresentationType type = RepresentationType::Vdw,
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

static bool filter_expression(const ApplicationData& data, str_t expr, md_exp_bitfield_t* mask);

static void modify_field(md_exp_bitfield_t* bf, const md_exp_bitfield_t* mask, SelectionOperator op) {
    switch(op) {
    case SelectionOperator::Or:
        md_bitfield_or(bf, bf, mask);
        break;
    case SelectionOperator::And:
        md_bitfield_and(bf, bf, mask);
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

// Global allocators for application
static md_allocator_i* frame_allocator = 0;
static md_allocator_i* persistent_allocator = default_allocator;

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
    md_bitfield_init(&data.representations.atom_visibility_mask, persistent_allocator);
    md_bitfield_init(&data.mold.script.frame_mask, persistent_allocator);

    md_semaphore_init(&data.mold.script.semaphore, 1);

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

    data.script.editor.SetLanguageDefinition(TextEditor::LanguageDefinition::VIAMD());
    data.script.editor.SetPalette(TextEditor::GetDarkPalette());

    load_dataset_from_file(&data, make_cstr(VIAMD_DATASET_DIR "/1ALA-500.pdb"));
    create_representation(&data, RepresentationType::Vdw, ColorMapping::Cpk, make_cstr("all"));
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

#if SHOW_IMGUI_DEMO_WINDOW
        ImGui::ShowDemoWindow();
        ImPlot::ShowDemoWindow();
#endif

        const int32_t num_frames = (int32_t)traj.num_frames;
        const int32_t last_frame = MAX(0, num_frames - 1);
        const double max_time = (double)MAX(0, last_frame);

        // #input
        if (data.ctx.input.key.hit[KEY_CONSOLE]) {
            if (data.console.Visible()) {
                data.console.Hide();
            } else if (!ImGui::GetIO().WantTextInput) {
                data.console.Show();
            }
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
                if (data.animation.mode == PlaybackMode::Stopped && data.animation.time == max_time) {
                    data.animation.time = 0;
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
                data.animation.time = CLAMP(data.animation.time + step, 0.0, max_time);
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

        // This needs to happen first (in imgui events) to enable docking of imgui windows
        ImGui::CreateDockspace();

        if (data.animation.mode == PlaybackMode::Playing) {
            data.animation.time += data.ctx.timing.delta_s * data.animation.fps;
            data.animation.time = CLAMP(data.animation.time, 0.0, max_time);
            if (data.animation.time >= max_time) {
                data.animation.mode = PlaybackMode::Stopped;
                data.animation.time = max_time;
            }
        }

        const int new_frame = CLAMP((int)(data.animation.time + 0.5f), 0, last_frame);
        const bool frame_changed = new_frame != data.animation.frame;
        data.animation.frame = new_frame;

        {
            static auto prev_time = data.animation.time;
            if (data.animation.time != prev_time) {
                time_changed = true;
                prev_time = data.animation.time;
            }
            else {
                time_changed = false;
            }
        }

        if (data.timeline.filter.temporal_window) {
            double pre_min = data.timeline.filter.min;
            double pre_max = data.timeline.filter.max;
            const auto half_window_ext = data.timeline.filter.window_extent * 0.5f;
            data.timeline.filter.min = CLAMP(data.animation.time - half_window_ext, 0.0, max_time);
            data.timeline.filter.max = CLAMP(data.animation.time + half_window_ext, 0.0, max_time);
            if (data.timeline.filter.min != pre_min || data.timeline.filter.max != pre_max) {
                data.mold.script.evaluate_filt = true;
            }
        }

        if (time_changed) {
            time_stopped = false;
            //data.hydrogen_bonds.dirty = true;

            PUSH_CPU_SECTION("Interpolate Position")
            if (traj.num_frames) {
                interpolate_atomic_properties(&data);
                //update_reference_frames(&data);

                /*
                if (data.ensemble_tracking.superimpose_structures) {
                    superimpose_ensemble(&data);
                }
                */

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
                if (rep.filter_is_dynamic || rep.color_mapping == ColorMapping::SecondaryStructure) {
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

                // Try aquire 2 semaphores
                if ( md_semaphore_try_aquire(&data.mold.script.semaphore)
                    // && md_semaphore_try_aquire(&data.mold.script.semaphore)
                    ) {
                    // Now we hold all resources for the script
                    data.mold.script.compile_ir = false;
                    data.mold.script.time_since_last_change = 0;

                    TextEditor& editor = data.script.editor;

                    std::string src = editor.GetText();
                    str_t str = {.ptr = src.data(), .len = (int64_t)src.length()};
                    md_script_ir_compile_args_t args = {
                        .src = str,
                        .mol = &data.mold.mol,
                        .alloc = default_allocator,
                    };

                    editor.ClearMarkers();
                    editor.ClearErrorMarkers();

                    if (md_script_ir_compile(&data.mold.script.ir, args)) {
                        data.mold.script.evaluate_full = true;
                        data.mold.script.evaluate_filt = true;
                    }

                    // Before we release the compute-dogs we want to allocate the data for the evaluations
                    md_script_eval_alloc(&data.mold.script.full_eval, data.mold.traj.num_frames, &data.mold.script.ir, persistent_allocator);
                    md_script_eval_alloc(&data.mold.script.filt_eval, data.mold.traj.num_frames, &data.mold.script.ir, persistent_allocator);

                    md_semaphore_release(&data.mold.script.semaphore);
                    //md_semaphore_release(&data.mold.script.semaphore);

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
                }
            }
        }

        if (data.mold.script.evaluate_full) {
            if (data.tasks.evaluate_full.id != 0) {
                md_script_eval_interrupt(&data.mold.script.full_eval);
            } else if (md_semaphore_try_aquire(&data.mold.script.semaphore)) {
                data.mold.script.evaluate_full = false;
                data.tasks.evaluate_full = task_system::enqueue_pool("Eval Full", 1, [&data](task_system::TaskSetRange) {
                    md_script_eval_args_t args = {
                        .ir = &data.mold.script.ir,
                        .mol = &data.mold.mol,
                        .traj = &data.mold.traj,
                        .frame_cache = &data.mold.frame_cache,
                        .filter_mask = NULL,
                    };

                    md_script_eval_compute(&data.mold.script.full_eval, args);

                    task_system::enqueue_main("Eval Complete", [&data]() {
                        data.tasks.evaluate_full.id = 0;
                        md_semaphore_release(&data.mold.script.semaphore);
                    });
                });
            }
        }

        if (data.mold.script.evaluate_filt) {
            if (data.tasks.evaluate_filt.id != 0) {
                md_script_eval_interrupt(&data.mold.script.filt_eval);
            } else if (md_semaphore_try_aquire(&data.mold.script.semaphore)) {
                data.mold.script.evaluate_filt = false;
                data.tasks.evaluate_filt = task_system::enqueue_pool("Eval Filt", 1, [&data](task_system::TaskSetRange) {
                    int64_t beg_frame = CLAMP((int64_t)data.timeline.filter.min, 0, data.mold.traj.num_frames);
                    int64_t end_frame = CLAMP((int64_t)data.timeline.filter.max, 0, data.mold.traj.num_frames);
                    end_frame = MAX(beg_frame + 1, end_frame);
                    md_bitfield_clear(&data.mold.script.frame_mask);
                    md_bitfield_set_range(&data.mold.script.frame_mask, beg_frame, end_frame);

                    md_script_eval_args_t args = {
                        .ir = &data.mold.script.ir,
                        .mol = &data.mold.mol,
                        .traj = &data.mold.traj,
                        .frame_cache = &data.mold.frame_cache,
                        .filter_mask = &data.mold.script.frame_mask
                    };

                    md_script_eval_compute(&data.mold.script.filt_eval, args);
                    task_system::enqueue_main("Eval Complete", [&data]() {
                        data.tasks.evaluate_filt.id = 0;
                        md_semaphore_release(&data.mold.script.semaphore);
                    });
                });
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

        update_view_param(&data);
        update_md_buffers(&data);

        //update_properties(&data);
        //update_density_volume_texture(&data);

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

        // GUI ELEMENTS
        data.console.Draw("VIAMD", data.ctx.window.width, data.ctx.window.height, (float)data.ctx.timing.delta_s);

        draw_main_menu(&data);
        draw_context_popup(&data);
        draw_async_task_window(&data);
        draw_animation_control_window(&data);
        draw_molecule_dynamic_info_window(&data);

        clear_highlight(&data);

        if (data.representations.show_window) draw_representations_window(&data);
        if (data.statistics.show_timeline_window) draw_timeline_window(&data);
        if (data.statistics.show_distribution_window) draw_distribution_window(&data);
        //if (data.ramachandran.show_window) draw_ramachandran_window(&data);
        //if (data.shape_space.show_window) draw_shape_space_window(&data);
        if (data.density_volume.show_window) draw_volume_window(&data);
        if (data.script.show_editor) draw_property_editor_window(&data);
        if (data.dataset.show_window) draw_dataset_window(&data);

        // @NOTE: ImGui::GetIO().WantCaptureMouse does not work with Menu
        if (!ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow)) {
            if (data.picking.idx != INVALID_PICKING_IDX) {
                draw_atom_info_window(data, data.picking.idx);
            }
        }

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

static void interpolate_atomic_properties(ApplicationData* data) {
    ASSERT(data);
    const auto& mol = data->mold.mol;
    const auto& traj = data->mold.traj;

    if (!mol.atom.count || !traj.num_frames) return;

    const int64_t last_frame = MAX(0LL, traj.num_frames - 1);
    const double time = CLAMP(data->animation.time, 0.0, double(last_frame));

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

    int64_t stride = ROUND_UP(mol.atom.count, md_simd_width);    // The interploation uses SIMD vectorization without bounds, so we make sure there is no overlap between the data segments
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
            md_frame_cache_load_frame_data(&data->mold.frame_cache, nearest_frame, mol.atom.x, mol.atom.y, mol.atom.z, box.elem, NULL, NULL, NULL);
            //load::traj::load_trajectory_frame_box(traj_ptr, (float*(*)[3])&box, nearest_frame);
            break;
        case InterpolationMode::Linear:
        {
            md_frame_cache_load_frame_data(&data->mold.frame_cache, frames[1], x[0], y[0], z[0], boxes[0].elem, NULL, NULL, NULL);
            md_frame_cache_load_frame_data(&data->mold.frame_cache, frames[2], x[1], y[1], z[1], boxes[1].elem, NULL, NULL, NULL);

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
            md_frame_cache_load_frame_data(&data->mold.frame_cache, frames[0], x[0], y[0], z[0], boxes[0].elem, NULL, NULL, NULL);
            md_frame_cache_load_frame_data(&data->mold.frame_cache, frames[1], x[1], y[1], z[1], boxes[1].elem, NULL, NULL, NULL);
            md_frame_cache_load_frame_data(&data->mold.frame_cache, frames[2], x[2], y[2], z[2], boxes[2].elem, NULL, NULL, NULL);
            md_frame_cache_load_frame_data(&data->mold.frame_cache, frames[3], x[3], y[3], z[3], boxes[3].elem, NULL, NULL, NULL);

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
    md_bitfield_copy(&prev_mask, mask);        
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

/*
static void grow_mask_by_radial_extent(Bitfield mask, const float* atom_x, const float* atom_y, const float* atom_z, int64_t atom_count, float extent) {
    Bitfield prev_mask;
    bitfield::init(&prev_mask, mask);
    defer { bitfield::free(&prev_mask); };

    spatialhash::Frame frame = spatialhash::compute_frame(atom_x, atom_y, atom_z, atom_count, vec3_t(extent));
    for (int64_t i = 0; i < atom_count; i++) {
        if (bitfield::get_bit(prev_mask, i)) {
            const vec3_t pos = {atom_x[i], atom_y[i], atom_z[i]};
            spatialhash::for_each_within(frame, pos, extent, [mask](int32_t idx, const vec3_t& pos) {
                UNUSED(pos);
                bitfield::set_bit(mask, idx);
            });
        }
    }
}
*/

static void expand_mask(md_exp_bitfield_t* mask, const md_range_t ranges[], int64_t num_ranges) {
    for (int64_t i = 0; i < num_ranges; i++) {
        if (md_bitfield_popcount_range( mask, ranges[i].beg, ranges[i].end) != 0) {
            md_bitfield_set_range(mask, ranges[i].beg, ranges[i].end);
        }
    }
}

static bool filter_expression(const ApplicationData& data, str_t expr, md_exp_bitfield_t* mask) {
    md_filter_stored_selection_t stored_sel[64] = {0};
    int64_t                      stored_sel_count = 0;
    /*
    * // @TODO: Reimplement selections with md_exp_bitfield type
    for (const auto& s : data.selection.stored_selections) {
        md_filter_stored_selection_t sel = {
            .ident = {.ptr = s.name.cstr(), .len = s.name.length() },
            .bitfield = {.bits = (uint64_t*)s.atom_mask.data(), .num_bits = s.atom_mask.size() }
        };
        stored_sel[stored_sel_count++] = sel;
        if (stored_sel_count == ARRAY_SIZE(stored_sel)) break;
    }
    */

    
    
    md_filter_context_t ctx {
        .mol = &data.mold.mol,
        .selection = {
            .count = stored_sel_count,
            .ptr = stored_sel,
        }
    };

    return md_filter_evaluate(expr, mask, ctx);
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
                            snprintf(res.path + res.path_len, ARRAY_SIZE(res.path) - res.path_len, "%.*s", (uint32_t)FILE_EXTENSION.len, FILE_EXTENSION.ptr);
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
                ImGui::Combo("Mode", (int*)(&data->view.mode), "Perspective\0Orthographic\0\0");
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
                ImGui::PushID("simulation_box");
                ImGui::ColorEdit4Minimal("Color", data->simulation_box.color.elem);
                ImGui::PopID();
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
            ImGui::Checkbox("Timelines", &data->statistics.show_timeline_window);
            ImGui::Checkbox("Distributions", &data->statistics.show_distribution_window);
            //ImGui::Checkbox("Ramachandran", &data->ramachandran.show_window);
            ImGui::Checkbox("Density Volumes", &data->density_volume.show_window);
            //ImGui::Checkbox("Shape Space", &data->shape_space.show_window);
            ImGui::Checkbox("Dataset", &data->dataset.show_window);

            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Selection")) {
            const auto atom_count = data->mold.mol.atom.count;
            md_exp_bitfield_t mask;
            md_bitfield_init(&mask, frame_allocator);
            defer { md_bitfield_free(&mask); };

            if (ImGui::MenuItem("Clear Selection")) {
                md_bitfield_clear(&data->selection.current_selection_mask);
                data->mold.dirty_buffers |= MolBit_DirtyFlags;
            }

            if (ImGui::MenuItem("Invert Selection")) {
                md_bitfield_clear(&data->selection.current_selection_mask);
                data->mold.dirty_buffers |= MolBit_DirtyFlags;
            }

            ImGui::Spacing();
            ImGui::Separator();

            // MODES
            ImGui::Combo("Level Mode", (int*)(&data->selection.level_mode), "Atom\0Residue\0Chain\0\0");

            ImGui::Spacing();
            ImGui::Separator();

            // QUERY
            {
                static char buf[256] = {0};
                static bool query_ok = false;

                ImGui::Text("Query");
                const auto TEXT_BG_DEFAULT_COLOR = ImGui::ColorConvertFloat4ToU32(ImGui::GetStyle().Colors[ImGuiCol_FrameBg]);
                ImGui::PushStyleColor(ImGuiCol_FrameBg, query_ok ? TEXT_BG_DEFAULT_COLOR : TEXT_BG_ERROR_COLOR);
                // ImGui::Combo("Mode", (int*)(&data->selection.op_mode), "Or\0And\0\0");
                const bool pressed_enter =
                    ImGui::InputText("##query", buf, ARRAY_SIZE(buf), ImGuiInputTextFlags_AutoSelectAll | ImGuiInputTextFlags_EnterReturnsTrue);
                ImGui::PopStyleColor();
                // ImGui::SameLine();
                if (!query_ok) ImGui::PushDisabled();
                const bool apply = ImGui::Button("Apply##query") || pressed_enter;
                if (!query_ok) ImGui::PopDisabled();

                if (pressed_enter) {
                    ImGui::SetKeyboardFocusHere(-1);
                }

                // if (ImGui::IsWindowAppearing()) {
                //    ImGui::SetKeyboardFocusHere(-1);
                //}

                {
                    query_ok = filter_expression(*data, str_from_cstr(buf), &mask);

                    if (query_ok) {
                        switch (data->selection.level_mode) {
                            case SelectionLevel::Atom:
                                // No need to expand the mask
                                break;
                            case SelectionLevel::Residue:
                                expand_mask(&mask, data->mold.mol.residue.atom_range, data->mold.mol.residue.count);
                                break;
                            case SelectionLevel::Chain:
                                expand_mask(&mask, data->mold.mol.chain.atom_range, data->mold.mol.chain.count);
                                break;
                            default:
                                ASSERT(false);
                        }
                    } else {
                        md_bitfield_clear(&mask);
                    }
                    // data->mold.dirty_buffers |= MolBit_DirtyFlags;
                }

                const bool show_preview = (ImGui::GetFocusID() == ImGui::GetID("##query")) || (ImGui::GetHoveredID() == ImGui::GetID("##query")) ||
                                          (ImGui::GetHoveredID() == ImGui::GetID("Apply##query"));

                if (show_preview) {
                    md_bitfield_copy(&data->selection.current_highlight_mask, &mask);
                    // if (query_ok) {
                    //memcpy(data->selection.current_highlight_mask.data(), mask.data(), mask.size_in_bytes());
                    /*
if (data->selection.op_mode == SelectionOperator::And) {
                            bitfield::and_field(data->selection.current_highlight_mask, data->selection.current_selection_mask, mask);
} else if (data->selection.op_mode == SelectionOperator::Or) {
                            bitfield::or_field(data->selection.current_highlight_mask, data->selection.current_selection_mask, mask);
}
                    */
                    data->mold.dirty_buffers |= MolBit_DirtyFlags;
                    //}
                }

                if (apply) {
                    md_bitfield_copy(&data->selection.current_selection_mask, &mask);

                    //memcpy(data->selection.current_selection_mask.data(), mask.data(), mask.size_in_bytes());
                    /*
if (data->selection.op_mode == SelectionOperator::And) {
                            bitfield::and_field(data->selection.current_selection_mask, data->selection.current_selection_mask, mask);
} else if (data->selection.op_mode == SelectionOperator::Or) {
                            bitfield::or_field(data->selection.current_selection_mask, data->selection.current_selection_mask, mask);
}
                    */
                    data->mold.dirty_buffers |= MolBit_DirtyFlags;
                }
                update_all_representations(data);
            }

            ImGui::Spacing();
            ImGui::Separator();

            // GROW
            {
                // static bool pos_dir = true;
                // static bool neg_dir = true;
                static float extent = 1.f;

                ImGui::Text("Grow");
                const bool mode_changed = ImGui::Combo("Mode", (int*)(&data->selection.grow_mode), "Covalent Bond\0Radial\0\0");
                const bool extent_changed = ImGui::SliderFloat("Extent", &extent, 1.0f, 16.f);
                const bool apply = ImGui::Button("Apply##grow");
                const bool show_preview = mode_changed || extent_changed || (ImGui::GetHoveredID() == ImGui::GetID("Extent")) ||
                                          (ImGui::GetHoveredID() == ImGui::GetID("Apply##grow")) || (ImGui::GetHoveredID() == ImGui::GetID("Mode"));

                if (show_preview) {
                    md_bitfield_copy(&mask, &data->selection.current_selection_mask);

                    switch (data->selection.grow_mode) {
                        case SelectionGrowth::CovalentBond:
                            grow_mask_by_covalent_bond(&mask, data->mold.mol.covalent_bond.bond, data->mold.mol.covalent_bond.count, (int64_t)extent);
                            break;
                        case SelectionGrowth::Radial: {
                            //const auto& mol = data->mold.mol;
                            //grow_mask_by_radial_extent(mask, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, extent);
                            break;
                        }
                        default:
                            ASSERT(false);
                    }

                    switch (data->selection.level_mode) {
                        case SelectionLevel::Atom:
                            // No need to expand the mask
                            break;
                        case SelectionLevel::Residue:
                            expand_mask(&mask, data->mold.mol.residue.atom_range, data->mold.mol.residue.count);
                            break;
                        case SelectionLevel::Chain:
                            expand_mask(&mask, data->mold.mol.chain.atom_range, data->mold.mol.chain.count);
                            break;
                        default:
                            ASSERT(false);
                    }

                    md_bitfield_copy(&data->selection.current_highlight_mask, &mask);
                    data->mold.dirty_buffers |= MolBit_DirtyFlags;

                    if (apply) {
                        md_bitfield_copy(&data->selection.current_selection_mask, &mask);
                    }
                    update_all_representations(data);
                }
            }

            ImGui::Spacing();
            ImGui::Separator();

            // STORED SELECTIONS
            {
                ImGui::Text("Stored Selections");
                const bool disable_new = md_bitfield_popcount(&data->selection.current_selection_mask) == 0;
                if (disable_new) ImGui::PushDisabled();
                if (ImGui::Button("Create New")) {
                    char name_buf[64];
                    snprintf(name_buf, 64, "sel%lli", md_array_size(data->selection.stored_selections) + 1);
                    create_selection(data, str_from_cstr(name_buf), &data->selection.current_selection_mask);
                }
                if (disable_new) ImGui::PopDisabled();
                for (int i = 0; i < (int)md_array_size(data->selection.stored_selections); i++) {
                    auto& sel = data->selection.stored_selections[i];
                    // const float32 item_width = CLAMP(ImGui::GetWindowContentRegionWidth() - 90.f, 100.f, 300.f);
                    StrBuf<128> name;
                    snprintf(name.cstr(), name.capacity(), "%s###ID", sel.name.cstr());

                    ImGui::PushID(i);
                    if (ImGui::CollapsingHeader(name.cstr())) {
                        if (ImGui::Button("Activate")) {
                            md_bitfield_copy(&data->selection.current_selection_mask, &sel.atom_mask);
                            data->mold.dirty_buffers |= MolBit_DirtyFlags;
                        }
                        ImGui::SameLine();
                        if (ImGui::DeleteButton("Remove")) {
                            remove_selection(data, i);
                        }
                        ImGui::SameLine();
                        if (ImGui::Button("Save")) {
                            md_bitfield_copy(&sel.atom_mask, &data->selection.current_selection_mask);
                            update_all_representations(data);
                        }
                    }

                    const auto h_id = ImGui::GetHoveredID();
                    bool show_preview = (h_id == ImGui::GetID(name.cstr()) || h_id == ImGui::GetID("Activate") || h_id == ImGui::GetID("Remove") ||
                                         h_id == ImGui::GetID("Clone"));

                    ImGui::PopID();

                    if (show_preview) {
                        md_bitfield_copy(&data->selection.current_highlight_mask, &sel.atom_mask);
                        data->mold.dirty_buffers |= MolBit_DirtyFlags;
                    }
                }
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

void draw_context_popup(ApplicationData* data) {
    ASSERT(data);

    bool valid_dynamic = data->mold.mol.atom.count && data->mold.traj.num_frames;

    const bool shift_down = data->ctx.input.key.down[Key::KEY_LEFT_SHIFT] || data->ctx.input.key.down[Key::KEY_RIGHT_SHIFT];
    if (data->ctx.input.mouse.clicked[1] && !shift_down && !ImGui::GetIO().WantTextInput) {
        if (data->selection.right_clicked != -1 && valid_dynamic) {
            ImGui::OpenPopup("AtomContextPopup");
        }
    }

    if (ImGui::BeginPopup("AtomContextPopup")) {
        /*
        if (data->selection.right_clicked != -1 && valid_dynamic) {
            if (ImGui::BeginMenu("Recenter Trajectory...")) {
                const int atom_idx = data->selection.right_clicked;
                AtomRange atom_range = {};

                if (ImGui::MenuItem("on Atom")) {
                    atom_range = {atom_idx, atom_idx + 1};
                }
                if (ImGui::IsItemHovered()) {
                }
                if (ImGui::MenuItem("on Residue")) {
                    const auto res_idx = data->mold.mol.atom.residue_idx[atom_idx];
                    atom_range = data->mold.mol.residue.atom_range[res_idx];
                }
                if (ImGui::MenuItem("on Chain")) {
                    const auto chain_idx = data->mold.mol.atom.chain_idx[atom_idx];
                    atom_range = data->mold.mol.chain.atom_range[chain_idx];
                }

                if (atom_range.ext() > 0) {
                    recenter_trajectory(&data->dynamic, atom_range);
                    interpolate_atomic_properties(data);
                    ImGui::CloseCurrentPopup();
                }

                ImGui::EndMenu();
            }
        }
        */
        ImGui::EndPopup();
    }
}

static void draw_animation_control_window(ApplicationData* data) {
    ASSERT(data);
    if (!data->mold.traj.num_frames) return;

    ImGui::Begin("Animation");
    const int32_t num_frames = (int32_t)data->mold.traj.num_frames;
    ImGui::Text("Num Frames: %i", num_frames);
    // ImGui::Checkbox("Apply post-interpolation pbc", &data->animation.apply_pbc);
    float t = (float)data->animation.time;
    if (ImGui::SliderFloat("Time", &t, 0, (float)(MAX(0, num_frames - 1)))) {
        data->animation.time = t;
    }
    //ImGui::SliderFloat("Speed", &data->animation.fps, 0.1f, 1000.f, "%.3f", 4.f);
    ImGui::SliderFloat("Speed", &data->animation.fps, 0.1f, 1000.f, "%.3f", ImGuiSliderFlags_Logarithmic);
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("Animation Speed in Frames Per Second");
    }
    if (ImGui::Combo("Interp.", (int*)(&data->animation.interpolation), "Nearest\0Linear\0Cubic\0\0")) {
        interpolate_atomic_properties(data);
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
        data->animation.time = 0.0;
    }
    ImGui::End();
}

static void draw_representations_window(ApplicationData* data) {
    // const auto old_hash = hash::crc64(data->representations.buffer);

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
        bool update_args = false;
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
            if (!rep.filter_is_valid) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
            if (ImGui::InputText("filter", rep.filter.cstr(), rep.filter.capacity(), ImGuiInputTextFlags_EnterReturnsTrue)) {
                update_color = true;
            }
            if (!rep.filter_is_valid) ImGui::PopStyleColor();
            ImGui::SameLine();
            ImGui::Checkbox("dynamic", &rep.filter_is_dynamic);
            if (ImGui::Combo("type", (int*)(&rep.type), "VDW\0Licorice\0Ribbons\0Cartoon\0\0")) {
                update_args = true;
                
            }
            if (ImGui::Combo("color mapping", (int*)(&rep.color_mapping),
                             "Uniform Color\0CPK\0Res Id\0Res Idx\0Chain Id\0Chain Idx\0Secondary Structure\0\0")) {
                update_color = true;
            }
            ImGui::PopItemWidth();
            if (rep.color_mapping == ColorMapping::Uniform) {
                ImGui::SameLine();
                if (ImGui::ColorEdit4("color", (float*)&rep.uniform_color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel)) {
                    update_color = true;
                }
            }
            ImGui::PushItemWidth(item_width);
            if (rep.type == RepresentationType::Vdw || rep.type == RepresentationType::Licorice) {
                if (ImGui::SliderFloat("radii scale", &rep.radius, 0.1f, 2.f)) update_args = true;
            }
            if (rep.type == RepresentationType::Ribbons) {
                if (ImGui::SliderFloat("spline tension", &rep.tension, 0.f, 1.f)) update_args = true;
                if (ImGui::SliderFloat("spline width", &rep.width, 0.1f, 2.f)) update_args = true;
                if (ImGui::SliderFloat("spline thickness", &rep.thickness, 0.1f, 2.f)) update_args = true;
            }
            ImGui::PopItemWidth();
            ImGui::Spacing();
            ImGui::Separator();
        }

        ImGui::PopID();
#if USE_MOLD
        if (update_args) {
            md_gl_representation_type_t type = MD_GL_REP_DEFAULT;
            md_gl_representation_args_t args = {};
            switch(rep.type) {
            case RepresentationType::Vdw:
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
#endif

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
        if (traj.num_frames) {
            ImGui::NewLine();
            ImGui::Text("TRAJ");
            const auto file = extract_file(data->files.trajectory);
            ImGui::Text("\"%.*s\"", (int32_t)file.len, file.ptr);
            ImGui::Text("# frames: %i", traj.num_frames);
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

void draw_property_menu_widgets(PropertyItem* items, int num_items) {
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
}

// #timeline
static void draw_timeline_window(ApplicationData* data) {
    ASSERT(data);
    ImGui::SetNextWindowSize(ImVec2(600, 300), ImGuiCond_FirstUseEver);

    if (ImGui::Begin("Temporal", &data->statistics.show_timeline_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_MenuBar)) {

        double pre_filter_min = data->timeline.filter.min;
        double pre_filter_max = data->timeline.filter.max;

        static PropertyItem s_props[32] = {};
        static int s_num_props = 0;
        static uint64_t s_fingerprint = 0;

        int num_time_values = (int)data->mold.traj.num_frames;
        float* time_values = (float*)md_alloc(frame_allocator, num_time_values * sizeof(float));
        defer {
            md_free(frame_allocator, time_values, num_time_values * sizeof(float));
        };
        for (int64_t i = 0; i < num_time_values; ++i) {
            //time_values[i] = (float)i / (float)(num_time_values - 1);
            time_values[i] = (float)i;
        }

        const double max_time_value = num_time_values > 0 ? time_values[num_time_values - 1] : 1.0;

        if (s_fingerprint != data->mold.script.full_eval.fingerprint) {
            PropertyItem new_props[ARRAY_SIZE(s_props)] = {};
            int num_new_props = 0;

            md_script_property_t* props = data->mold.script.full_eval.properties;
            for (int64_t i = 0; i < data->mold.script.full_eval.num_properties; ++i) {
                if (props[i].type != MD_SCRIPT_PROPERTY_TYPE_TEMPORAL) continue;

                bool show = true;
                for (int64_t j = 0; j < s_num_props; ++j) {
                    if (compare_str(props[i].ident, s_props[j].lbl)) {
                        show = s_props[j].show;
                    }
                }
                PropertyItem p = {
                    .lbl = props[i].ident,
                    .col = PROPERTY_COLORS[i % ARRAY_SIZE(PROPERTY_COLORS)],
                    .idx = (int)i,
                    .show = show
                };

                new_props[num_new_props++] = p;
            }
            
            memcpy(s_props, new_props, sizeof(new_props));
            s_num_props = num_new_props;
        }

        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("Properties")) {
                //draw_property_menu_widgets(s_props, s_num_props);
                
                if (s_num_props) {
                    for (int i = 0; i < s_num_props; ++i) {
                        ImPlot::ItemIcon(s_props[i].col); ImGui::SameLine();
                        ImGui::Selectable(s_props[i].lbl.cstr(), &s_props[i].show, 0, ImVec2(50, 0));
                    }
                } else {
                    ImGui::Text("No properties to show.");
                }
                
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Filter")) {
                ImGui::Checkbox("Enabled", &data->timeline.filter.enabled);
                if (data->timeline.filter.enabled) {
                    ImGui::Checkbox("Temporal Window", &data->timeline.filter.temporal_window);
                    if (data->timeline.filter.temporal_window) {
                        ImGui::SliderFloat("Extent", &data->timeline.filter.window_extent, 1.0f, (float)max_time_value);
                    }
                }

                ImGui::EndMenu();
            }

            ImGui::EndMenuBar();
        }

        if (ImGui::IsWindowFocused() && ImGui::IsKeyPressed(KEY_PLAY_PAUSE, false)) {
            data->animation.mode = data->animation.mode == PlaybackMode::Playing ? PlaybackMode::Stopped : PlaybackMode::Playing;
        }

        /*
        // convenience struct to manage DND items; do this however you like
        struct DndItem {
            char    label[16] = {};
            int     num_values = 0;
            float*  values = NULL;
            float*  variance = NULL;
            ImVec4  color = {};
            int     plot = 0;
        };

        static DndItem  dnd[32];
        int             num_dnd = 0;
        */

        // This is for the time stamps along the x-axis.
        // If we don't have any specific time points for the frames, we just use the indices as time points.



        if (num_time_values > 0) {

            /*
            for (int64_t i = 0; i < data->mold.script.full_eval.num_properties; i++) {
                auto& prop = data->mold.script.full_eval.properties[i];
                if (prop.type != MD_SCRIPT_PROPERTY_TYPE_TEMPORAL) continue;

                ASSERT(num_dnd < ARRAY_SIZE(dnd));
                int idx = num_dnd++;
                dnd[idx] = {
                    .label = {0},
                    .num_values = prop.data.aggregate ? (int)prop.data.aggregate->num_values : (int)prop.data.num_values,
                    .values = prop.data.aggregate ? prop.data.aggregate->mean : prop.data.values,
                    .variance = prop.data.aggregate ? prop.data.aggregate->variance : NULL,
                    .color = vec_cast(qualitative_color_scale(idx)),
                    .plot = dnd[idx].plot
                };
                const size_t cpy_size = ARRAY_SIZE(dnd[num_dnd-1].label) < prop.ident.len ? ARRAY_SIZE(dnd[num_dnd-1].label) : prop.ident.len;
                strncpy(dnd[num_dnd-1].label, prop.ident.ptr, cpy_size);
            }
            */

            /*
            // child window to serve as initial source for our DND items
            ImGui::BeginChild("DND_LEFT",ImVec2(100,-1));
            if (ImGui::Button("Reset", ImVec2(100, 0))) {
                for (int i = 0; i < num_dnd; ++i)
                    dnd[i].plot = 0;
            }
            for (int k = 0; k < num_dnd; ++k) {
                if (dnd[k].plot > 0)
                    continue;
                ImPlot::ItemIcon(dnd[k].color); ImGui::SameLine();
                ImGui::Selectable(dnd[k].label, false, 0, ImVec2(100, 0));
                if (ImGui::BeginDragDropSource(ImGuiDragDropFlags_None)) {
                    ImGui::SetDragDropPayload("MY_DND", &k, sizeof(int));
                    ImPlot::ItemIcon(dnd[k].color); ImGui::SameLine();
                    ImGui::TextUnformatted(dnd[k].label);
                    ImGui::EndDragDropSource();
                }
            }
            ImGui::EndChild();
            */
        
            /*
            if (ImGui::BeginDragDropTarget()) {
                if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload("MY_DND")) {
                    int i = *(int*)payload->Data;
                    dnd[i].plot = 0;
                }
                ImGui::EndDragDropTarget();
            }
            */

            //ImGui::SameLine();
            //ImGui::BeginChild("DND_RIGHT",ImVec2(-1,-1));
        
            ImPlotAxisFlags axis_flags = 0;
            ImPlotAxisFlags axis_flags_x = axis_flags;
            ImPlotAxisFlags axis_flags_y = axis_flags | ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_RangeFit;
            ImPlotFlags flags = ImPlotFlags_AntiAliased;
    /*
            static bool need_refit = false;
            if (need_refit) {
                ImPlot::FitNextPlotAxes();
                need_refit = false;
            }
            */

            static double x_min = 0;
            static double x_max = 1;

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

            const int64_t num_props = data->mold.script.full_eval.num_properties;

            static bool is_dragging = false;
            static bool is_selecting = false;

            if (!ImGui::IsMouseDown(0)) {
                is_dragging = false;
                is_selecting = false;
            }

            data->timeline.view_range.min = MAX(data->timeline.view_range.min, 0);
            data->timeline.view_range.max = MIN(data->timeline.view_range.max, num_time_values);

            if (num_props == 0) {
                ImPlot::LinkNextPlotLimits(&data->timeline.view_range.min, &data->timeline.view_range.max, 0, 0);
                if (ImPlot::BeginPlot("##Timeline", NULL, NULL, ImVec2(-1,150), flags, axis_flags_x, axis_flags_y)) {
                    if (data->timeline.filter.enabled) {
                        ImPlot::DragRangeX("Time Filter", &data->timeline.filter.min, &data->timeline.filter.max, 0, max_time_value);
                    }
                    ImPlot::DragLineX("Current Time", &data->animation.time, true, ImVec4(1,1,0,1));
                    ImPlot::EndPlot();
                }
            }
            else {
                for (int64_t i = 0; i < s_num_props; ++i) {
                    if (s_props[i].show == false) continue;
                    const md_script_property_t& prop = data->mold.script.full_eval.properties[s_props[i].idx];

                    ImGui::PushID(i);
                    ImPlot::LinkNextPlotLimits(&data->timeline.view_range.min, &data->timeline.view_range.max, 0, 0);
                    if (ImPlot::BeginPlot("##Timeline", NULL, NULL, ImVec2(-1,150), flags, axis_flags_x, axis_flags_y)) {


                        if (ImGui::IsItemHovered() && ImGui::IsItemActive()) {
                            double x = ImPlot::GetPlotMousePos().x;
                            if (is_dragging) {
                                data->animation.time = x;
                            }
                            else if (is_selecting) {
                                data->timeline.filter.max = CLAMP(x, data->timeline.filter.min, max_time_value);
                            }
                            else if (ImGui::IsMouseDown(0)) {
                                switch (ImGui::GetIO().KeyMods) {
                                case ImGuiKeyModFlags_Ctrl:
                                    is_dragging = true;
                                    break;
                                case ImGuiKeyModFlags_Shift:
                                    if (!data->timeline.filter.temporal_window) {
                                        is_selecting = true;
                                        data->timeline.filter.min = CLAMP(x, 0, max_time_value);
                                    }
                                    break;
                                default:
                                    break;
                                };
                            }
                        }

                        ImPlotPlot* plot = ImPlot::GetCurrentPlot();
                        plot->XAxis.Range.Min = MAX(plot->XAxis.Range.Min, 0);
                        plot->XAxis.Range.Max = MIN(plot->XAxis.Range.Max, max_time_value);

                        if (data->timeline.filter.enabled) {
                            const bool disabled = data->timeline.filter.temporal_window || is_selecting;
                            if (disabled) ImGui::PushDisabled();
                            ImPlot::DragRangeX("Time Filter", &data->timeline.filter.min, &data->timeline.filter.max, 0, max_time_value);
                            if (disabled) ImGui::PopDisabled();
                        }

                        ImVec4 prop_col = ImGui::ColorConvertU32ToFloat4(s_props[i].col);
                        ImPlot::SetNextLineStyle(prop_col);

                        if (prop.data.aggregate) {
                            ASSERT(num_time_values == prop.data.aggregate->num_values);

                            struct ShadedData {
                                float* x_vals;
                                float* y_mean;
                                float* y_variance;
                            } user_data = {
                                .x_vals = time_values,
                                .y_mean = prop.data.aggregate->mean,
                                .y_variance = prop.data.aggregate->variance
                            };
                            ImPlot::SetNextFillStyle(prop_col, 0.2f);
                            ImPlot::PlotShadedG(prop.ident.ptr,
                                [](void* payload, int idx) -> ImPlotPoint {
                                    ShadedData* data = (ShadedData*)payload;
                                    return ImPlotPoint(data->x_vals[idx], data->y_mean[idx] - data->y_variance[idx]);
                                },
                                &user_data,
                                [](void* payload, int idx) -> ImPlotPoint {
                                    ShadedData* data = (ShadedData*)payload;
                                    return ImPlotPoint(data->x_vals[idx], data->y_mean[idx] + data->y_variance[idx]);
                                },
                                &user_data,
                                prop.data.aggregate->num_values);

                            
                            ImPlot::PlotLine(prop.ident.ptr, time_values, prop.data.aggregate->mean, prop.data.aggregate->num_values);
                        } else {
                            ASSERT(num_time_values == prop.data.num_values);
                            ImPlot::PlotLine(prop.ident.ptr, time_values, prop.data.values, prop.data.num_values);
                        }

                        /*
                        for (int k = 0; k < num_dnd; ++k) {
                        if (dnd[k].plot == 1 && dnd[k].num_values > 0) {
                        ASSERT(dnd[k].num_values == num_time_values);
                        ImPlot::SetNextLineStyle(dnd[k].color);
                        ImPlot::PlotLine(dnd[k].label, time_values, dnd[k].values, dnd[k].num_values);
                        // allow legend item labels to be DND sources
                        if (ImPlot::BeginDragDropSourceItem(dnd[k].label)) {
                        ImGui::SetDragDropPayload("MY_DND", &k, sizeof(int));
                        ImPlot::ItemIcon(dnd[k].color); ImGui::SameLine();
                        ImGui::TextUnformatted(dnd[k].label);
                        ImPlot::EndDragDropSource();
                        }
                        }
                        }
                        */

                        ImPlot::DragLineX("Current Time", &data->animation.time, true, ImVec4(1,1,0,1));

                        /*
                        // allow the main plot area to be a DND target
                        if (ImPlot::BeginDragDropTarget()) {
                        if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload("MY_DND")) {
                        int i = *(int*)payload->Data;
                        dnd[i].plot = 1;
                        need_refit = true;
                        }
                        ImPlot::EndDragDropTarget();
                        }

                        // allow the legend to be a DND target
                        if (ImPlot::BeginDragDropTargetLegend()) {
                        if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload("MY_DND")) {
                        int i = *(int*)payload->Data;
                        dnd[i].plot = 1;
                        need_refit = true;
                        }
                        ImPlot::EndDragDropTarget();
                        }
                        */

                        ImPlot::EndPlot();
                    }
                    ImGui::PopID();
                }
            }

            ImPlot::GetInputMap() = old_map;

        }

        //ImGui::EndChild();

        
        if (data->timeline.filter.enabled && data->timeline.filter.min != pre_filter_min || data->timeline.filter.max != pre_filter_max) {
            data->mold.script.evaluate_filt = true;
        }
    }
    ImGui::End();
}

static void compute_histogram(float* bins, int num_bins, float min_bin_val, float max_bin_val, const float* values, int num_values) {
    memset(bins, 0, sizeof(float) * num_bins);

    float bin_range = max_bin_val - min_bin_val;
    float inv_range = 1.0f / bin_range;
    for (int i = 0; i < num_values; ++i) {
        int idx = CLAMP((int)(((values[i] - min_bin_val) * inv_range) * num_bins), 0, num_bins - 1);
        bins[idx] += 1;
    }

    const float scl = 1.0f / num_values;
    for (int i = 0; i < num_bins; ++i) {
        bins[i] *= scl;
    }
}

static void draw_distribution_window(ApplicationData* data) {
    ImGui::SetNextWindowSize(ImVec2(200, 300), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Distributions", &data->statistics.show_distribution_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_MenuBar)) {
        const int MIN_NUM_BINS = 8;
        const int MAX_NUM_BINS = 1024;

        static PropertyItem s_props[32] = {};
        static int s_num_props = 0;
        static uint64_t s_fingerprint = 0;
        static int num_bins = 128;

        if (s_fingerprint != data->mold.script.full_eval.fingerprint) {
            PropertyItem new_props[ARRAY_SIZE(s_props)] = {};
            int num_new_props = 0;

            md_script_property_t* props = data->mold.script.full_eval.properties;
            for (int64_t i = 0; i < data->mold.script.full_eval.num_properties; ++i) {
                if (props[i].type != MD_SCRIPT_PROPERTY_TYPE_TEMPORAL && props[i].type != MD_SCRIPT_PROPERTY_TYPE_DISTRIBUTION) continue;

                bool show = true;
                for (int64_t j = 0; j < s_num_props; ++j) {
                    if (compare_str(props[i].ident, s_props[j].lbl)) {
                        show = s_props[j].show;
                    }
                }
                PropertyItem p = {
                    .lbl = props[i].ident,
                    .col = PROPERTY_COLORS[i % ARRAY_SIZE(PROPERTY_COLORS)],
                    .idx = (int)i,
                    .show = show
                };

                new_props[num_new_props++] = p;
            }

            memcpy(s_props, new_props, sizeof(new_props));
            s_num_props = num_new_props;
        }

        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("Properties")) {
                //draw_property_menu_widgets(s_props, s_num_props);

                if (s_num_props) {
                    for (int i = 0; i < s_num_props; ++i) {
                        ImPlot::ItemIcon(s_props[i].col); ImGui::SameLine();
                        ImGui::Selectable(s_props[i].lbl.cstr(), &s_props[i].show, 0, ImVec2(50, 0));
                    }
                } else {
                    ImGui::Text("No properties to show.");
                }

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

            ImGui::EndMenuBar();
        }

        if (ImGui::IsWindowFocused() && ImGui::IsKeyPressed(KEY_PLAY_PAUSE, false)) {
            data->animation.mode = data->animation.mode == PlaybackMode::Playing ? PlaybackMode::Stopped : PlaybackMode::Playing;
        }

        ImPlotAxisFlags axis_flags = 0;
        ImPlotAxisFlags axis_flags_x = axis_flags | 0;
        ImPlotAxisFlags axis_flags_y = axis_flags | 0;
        ImPlotFlags flags = ImPlotFlags_AntiAliased;

        // The distribution properties are always computed as histograms with a resolution of 1024
        // If we have a different number of bins for our visualization, we need to recompute the bins
        float* bins = (float*)md_alloc(frame_allocator, num_bins * sizeof(float));
        float* filtered_bins = (float*)md_alloc(frame_allocator, num_bins * sizeof(float));

        for (int i = 0; i < s_num_props; ++i) {
            if (s_props[i].show == false) continue;
            md_script_property_t& prop = data->mold.script.full_eval.properties[s_props[i].idx];
            md_script_property_t& filt_prop = data->mold.script.filt_eval.properties[s_props[i].idx];

            if (prop.type != MD_SCRIPT_PROPERTY_TYPE_DISTRIBUTION && prop.type != MD_SCRIPT_PROPERTY_TYPE_TEMPORAL) continue;

            float min_x = prop.data.min_range[0];
            float max_x = prop.data.max_range[0];
            float min_y = prop.data.min_value;
            float max_y = prop.data.max_value;
            float range_x = max_x - min_x;
            float* draw_bins = 0;
            float* draw_filtered_bins = 0;

            if (prop.type == MD_SCRIPT_PROPERTY_TYPE_DISTRIBUTION) {
                if (num_bins != prop.data.num_values) {
                    memset(bins, 0, num_bins * sizeof(float));

                    // Downsample bins
                    const int factor = (int)prop.data.num_values / num_bins;
                    ASSERT(factor > 1);

                    for (int64_t j = 0; j < prop.data.num_values; ++j) {
                        int idx = j / factor;
                        bins[idx] += prop.data.values[j];
                    }

                    draw_bins = bins;
                }
                else {
                    draw_bins = prop.data.values;
                }
            }
            else if (prop.type == MD_SCRIPT_PROPERTY_TYPE_TEMPORAL) {
                compute_histogram(bins, num_bins, min_x, max_x, prop.data.values, prop.data.num_values);
                draw_bins = bins;

                if (data->timeline.filter.enabled) {
                    compute_histogram(filtered_bins, num_bins, min_x, max_x, filt_prop.data.values, filt_prop.data.num_values);
                    draw_filtered_bins = filtered_bins;
                }
            }

            if (draw_bins) {
                min_y = FLT_MAX;
                max_y = -FLT_MAX;
                for (int j = 0; j < num_bins; ++j) {
                    min_y = MIN(min_y, draw_bins[j]);
                    max_y = MAX(max_y, draw_bins[j]);
                }
            }

            char label[16] = {0};
            strncpy(label, prop.ident.ptr, MIN(ARRAY_SIZE(label) - 1, prop.ident.len));
            
            ImGui::PushID(i);
            const double bar_width = (max_x - min_x) / (num_bins-1);
            const double bar_off = min_x;
            const double bar_scl = (max_x - min_x) / (num_bins-1);

            ImPlot::SetNextPlotLimits(min_x, max_x, 0, max_y * 1.1, ImGuiCond_Always);

            if (ImPlot::BeginPlot(label, 0, 0, ImVec2(-1,150), flags, axis_flags_x, axis_flags_y)) {
                //ImPlot::SetNextFillStyle(vec_cast(qualitative_color_scale(i)), 0.5f);
                //ImPlot::PlotBars(label, draw_bins, num_bins, bar_width, bar_off);

                struct PlotData {
                    double offset;
                    double scale;
                    float* bars;
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

                plot_data.bars = draw_bins;
                ImPlot::SetNextFillStyle(ImVec4(0,0,0,-1), 1.0f);
                ImPlot::SetNextLineStyle(ImVec4(0,0,0,0), 0);
                ImPlot::PlotBarsG(label, getter, &plot_data, num_bins, bar_width);

                if (draw_filtered_bins) {
                    char buf[64];
                    snprintf(buf, sizeof(buf), "%s filt", label);
                    plot_data.bars = draw_filtered_bins;
                    ImPlot::SetNextFillStyle(ImVec4(1,1,0,1), 0.3f);
                    ImPlot::SetNextLineStyle(ImVec4(0,0,0,0), 0);
                    ImPlot::PlotBarsG(buf, getter, &plot_data, num_bins, bar_width);
                }
                //ImPlot::PlotHistogram(label, prop.data.values, prop.data.num_values, bins, false, false, range);
                ImPlot::EndPlot();
            }
            ImGui::PopID();
        }
    }
    ImGui::End();
}

static void draw_shapespace_window(ApplicationData* data) {
    
}

#if 0
static void draw_ramachandran_window(ApplicationData* data) {
    // const int32 num_frames = data->mold.traj ? data->mold.traj.num_frames : 0;
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

static void draw_volume_window(ApplicationData* data) {
    ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Density Volume", &data->density_volume.show_window, ImGuiWindowFlags_MenuBar)) {
        const ImVec2 button_size = {160, 20};

        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("File")) {
                if (ImGui::MenuItem("Export volume")) {
                    auto res = application::file_dialog(application::FileDialogFlags_Open, {}, make_cstr("raw"));
                    if (res.result == application::FileDialogResult::Ok) {
                        //volume::write_to_file(data->density_volume.volume, res.path);
                        md_print(MD_LOG_TYPE_INFO, "Wrote density volume");
                    }
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Property")) {
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

        char combo_buf[512] = {0};
        int  combo_buf_len = 0;
        int  combo_prop_idx[32] = {0};
        int  combo_size = 0;

        for (int64_t i = 0; i < data->mold.script.full_eval.num_properties; ++i) {
            md_script_property_t* prop = &data->mold.script.full_eval.properties[i];
            if (prop->type == MD_SCRIPT_PROPERTY_TYPE_VOLUME && prop->data.values) {
                combo_buf_len += snprintf(combo_buf + combo_buf_len, ARRAY_SIZE(combo_buf) - combo_buf_len, "%.*s\0", (int)prop->ident.len, prop->ident.ptr);
                combo_buf_len += 1; // Apparently snprintf does not print the final \0
                ASSERT(combo_size < ARRAY_SIZE(combo_prop_idx));
                combo_prop_idx[combo_size] = (int)i;
                combo_size += 1;
            }
        }

        static uint64_t s_fingerprint = 0;
        static md_script_property_t* curr_prop = NULL;

        md_script_property_t* prop = NULL;
        if (combo_size > 0) {
            static int curr_idx = 0;
            curr_idx = CLAMP(curr_idx, 0, combo_size - 1);
            if (combo_buf[0] != '\0') {
                ImGui::SetNextItemWidth(-1);
                ImGui::Combo("Volume", &curr_idx, combo_buf);
            }
            int prop_idx = combo_prop_idx[curr_idx];
            prop = &data->mold.script.full_eval.properties[prop_idx];
        }

        bool update_representations = false;
        if (curr_prop != prop) {
            curr_prop = prop;
            update_representations = true;
        }

        if (s_fingerprint != data->mold.script.full_eval.fingerprint) {
            s_fingerprint = data->mold.script.full_eval.fingerprint;
            update_representations = true;
        }

        // Canvas
        // Using InvisibleButton() as a convenience 1) it will advance the layout cursor and 2) allows us to use IsItemHovered()/IsItemActive()
        ImVec2 canvas_p0 = ImGui::GetCursorScreenPos();      // ImDrawList API uses screen coordinates!
        ImVec2 canvas_sz = ImGui::GetContentRegionAvail();   // Resize canvas to what's available
        if (canvas_sz.x < 50.0f) canvas_sz.x = 50.0f;
        if (canvas_sz.y < 50.0f) canvas_sz.y = 50.0f;
        ImVec2 canvas_p1 = ImVec2(canvas_p0.x + canvas_sz.x, canvas_p0.y + canvas_sz.y);
        // 
        // Draw border and background color
        ImGuiIO& io = ImGui::GetIO();
        ImDrawList* draw_list = ImGui::GetWindowDrawList();
        draw_list->AddImage((ImTextureID)data->density_volume.fbo.deferred.post_tonemap, canvas_p0, canvas_p1, {0,1}, {1,0});
        //draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));
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
            .max_distance = 100.0,
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

        mat4_t model_mat = volume::compute_model_to_world_matrix({0,0,0}, {1,1,1});
        mat4_t view_mat = camera_world_to_view_matrix(data->density_volume.camera);
        mat4_t proj_mat = camera_perspective_projection_matrix(data->density_volume.camera, (int)canvas_sz.x, (int)canvas_sz.y);

        clear_gbuffer(&gbuf);
        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);

        if (prop) {
            md_script_visualization_t vis = {};
            md_script_visualization_args_t args = {
                .token = prop->vis_token,
                .ir = &data->mold.script.ir,
                .mol = &data->mold.mol,
                .traj = &data->mold.traj,
                .frame_cache = &data->mold.frame_cache,
                .alloc = frame_allocator,
            };
            md_script_visualization_init(&vis, args);

            const float s = vis.sdf.extent;
            vec3_t min_aabb = {-s, -s, -s};
            vec3_t max_aabb = {s, s, s};
            model_mat = volume::compute_model_to_world_matrix(min_aabb, max_aabb);

            int64_t num_reps = vis.sdf.count;

            if (update_representations) {
                if (!data->density_volume.volume_texture.id) {
                    gl::init_texture_3D(&data->density_volume.volume_texture.id, prop->data.dim[0], prop->data.dim[1], prop->data.dim[2], GL_R32F);
                    data->density_volume.volume_texture.dim_x = prop->data.dim[0];
                    data->density_volume.volume_texture.dim_y = prop->data.dim[1];
                    data->density_volume.volume_texture.dim_z = prop->data.dim[2];
                    data->density_volume.volume_texture.max_value = prop->data.max_value;
                }
                gl::set_texture_3D_data(data->density_volume.volume_texture.id, prop->data.values, GL_R32F);

                if (data->density_volume.gl_reps) {
                    // Only free those required
                    for (int64_t i = num_reps; i < md_array_size(data->density_volume.gl_reps); ++i) {
                        md_gl_representation_free(&data->density_volume.gl_reps[i]);
                    }
                }
                md_array_resize(data->density_volume.gl_reps, num_reps, persistent_allocator);

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
                }
            }

            if (data->density_volume.show_reference_structures) {
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
                    mats[i] = &vis.sdf.matrices[i].elem[0][0];
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
            }
        }

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.deferred.fbo);
        glDrawBuffer(GL_COLOR_ATTACHMENT_POST_TONEMAP);
        glViewport(0, 0, gbuf.width, gbuf.height);

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

        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);

        if (data->density_volume.show_bounding_box) {
            const vec3_t min_box = {0,0,0};
            const vec3_t max_box = {1,1,1};

            immediate::set_model_view_matrix(mat4_mul(view_mat, model_mat));
            immediate::set_proj_matrix(proj_mat);

            uint32_t box_color = convert_color(data->density_volume.bounding_box_color);
            uint32_t clip_color = convert_color(data->density_volume.clip_volume_color);
            immediate::draw_box_wireframe(min_box, max_box, box_color);
            immediate::draw_box_wireframe(data->density_volume.clip_volume.min, data->density_volume.clip_volume.max, clip_color);

            immediate::flush();
        }

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
                .model = model_mat,
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

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glDrawBuffer(GL_BACK);
    }

    ImGui::End();
}

static void draw_dataset_window(ApplicationData* data) {
    ImGui::Begin("Dataset", &data->dataset.show_window);
    ImGui::SetWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);

    str_t mol_file = data->files.molecule;
    ImGui::Text("Molecular data: %.*s", (int)mol_file.len, mol_file.ptr);
    ImGui::Text("Num atoms:    %9lli", data->mold.mol.atom.count);
    ImGui::Text("Num residues: %9lli", data->mold.mol.residue.count);
    ImGui::Text("Num chains:   %9lli", data->mold.mol.chain.count);

    str_t traj_file = data->files.trajectory;
    if (traj_file.len) {
        ImGui::Separator();
        ImGui::Text("Trajectory data: %.*s", (int)traj_file.len, traj_file.ptr);
        ImGui::Text("Num frames:    %9lli", data->mold.traj.num_frames);
        ImGui::Text("Num atoms:     %9lli", data->mold.traj.num_atoms);
    }

    ImGui::End();
}

static void draw_property_editor_window(ApplicationData* data) {
    ASSERT(data);

    TextEditor& editor = data->script.editor;

    if (ImGui::Begin("Property Editor", &data->script.show_editor, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_MenuBar)) {
        ImGui::SetWindowSize(ImVec2(800, 600), ImGuiCond_FirstUseEver);
        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("File"))
            {
                if (ImGui::MenuItem("Save"))
                {
                    auto textToSave = editor.GetText();
                    /// save text....
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Edit"))
            {
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

            if (ImGui::BeginMenu("View"))
            {
                if (ImGui::MenuItem("Dark palette"))
                    editor.SetPalette(TextEditor::GetDarkPalette());
                if (ImGui::MenuItem("Light palette"))
                    editor.SetPalette(TextEditor::GetLightPalette());
                if (ImGui::MenuItem("Retro blue palette"))
                    editor.SetPalette(TextEditor::GetRetroBluePalette());
                ImGui::EndMenu();
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
                .frame_cache = &data->mold.frame_cache,
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
                immediate::draw_line(vertices[idx[0]], vertices[idx[1]]);
            }

            for (int64_t i = 0; i < vis.triangle.count; ++i) {
                ASSERT(vis.triangle.idx);
                uint16_t idx[3] = { vis.triangle.idx[i * 3 + 0], vis.triangle.idx[i * 3 + 1], vis.triangle.idx[i * 3 + 2] };
                immediate::draw_triangle(vertices[idx[0]], vertices[idx[1]], vertices[idx[2]], 0x3300FFFF);
            }

            immediate::flush();

            glEnable(GL_CULL_FACE);

            if (!md_bitfield_empty(&vis.atom_mask)) {
                md_bitfield_copy(&data->selection.current_highlight_mask, &vis.atom_mask);
                data->mold.dirty_buffers |= MolBit_DirtyFlags;
            }
            md_script_visualization_free(&vis);
        }

        ImGui::End();
    }

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

    if (data->mold.dirty_buffers & MolBit_DirtySecondaryStructure) {
        if (mol.backbone.secondary_structure) {
            md_gl_molecule_set_backbone_secondary_structure(&data->mold.gl_mol, 0, (uint32_t)mol.backbone.count, mol.backbone.secondary_structure, 0);
        }
    }

    data->mold.dirty_buffers = 0;
}

static void interrupt_async_tasks(ApplicationData* data) {
    if (data->tasks.prefetch_trajectory.id != 0) {
        task_system::interrupt_task(data->tasks.prefetch_trajectory);
    }

    if (data->tasks.evaluate_full.id != 0) {
        md_script_eval_interrupt(&data->mold.script.full_eval);
        task_system::interrupt_task(data->tasks.evaluate_full);
    }

    if (data->tasks.evaluate_filt.id != 0) {
        md_script_eval_interrupt(&data->mold.script.filt_eval);
        task_system::interrupt_task(data->tasks.evaluate_filt);
    }

    task_system::wait_for_task(data->tasks.evaluate_full);
    task_system::wait_for_task(data->tasks.evaluate_filt);

    data->tasks.evaluate_full.id = 0;
    data->tasks.evaluate_filt.id = 0;
}

// #moleculedata
static void free_trajectory_data(ApplicationData* data) {
    ASSERT(data);
    interrupt_async_tasks(data);

    if (data->mold.traj.num_frames) {
        md_frame_cache_free(&data->mold.frame_cache);
        load::traj::close(&data->mold.traj);
    }
}

static void free_molecule_data(ApplicationData* data) {
    ASSERT(data);
    interrupt_async_tasks(data);

    if (data->mold.mol.atom.count) {
        data->files.molecule = "";
        load::mol::free(&data->mold.mol);
    }
    free_trajectory_data(data);

    md_bitfield_clear(&data->selection.current_selection_mask);
    md_bitfield_clear(&data->selection.current_highlight_mask);
    md_script_ir_free(&data->mold.script.ir);
    md_script_eval_free(&data->mold.script.full_eval);
}

static void init_molecule_data(ApplicationData* data) {
    if (data->mold.mol.atom.count) {
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

static void init_trajectory_data(ApplicationData* data) {
    if (data->mold.traj.num_frames) {
        md_frame_cache_init(&data->mold.frame_cache, &data->mold.traj, persistent_allocator, 0);

        const double min_time = 0;
        const double max_time = (double)data->mold.traj.num_frames;

        data->timeline.view_range = {min_time, max_time};
        data->timeline.filter.min = min_time;
        data->timeline.filter.max = max_time;
        data->animation.time = CLAMP(data->animation.time, min_time, max_time);
        data->animation.frame = CLAMP((int32_t)data->animation.time, (int32_t)min_time, (int32_t)max_time - 1);
        int32_t frame = data->animation.frame;

        md_frame_cache_load_frame_data(&data->mold.frame_cache, frame, data->mold.mol.atom.x, data->mold.mol.atom.y, data->mold.mol.atom.z, data->simulation_box.box.elem, NULL, NULL, NULL);
        data->mold.dirty_buffers |= MolBit_DirtyPosition;

        update_md_buffers(data);
        md_gl_molecule_update_atom_previous_position(&data->mold.gl_mol); // Do this explicitly to update the previous position to avoid motion blur trails

        if (data->mold.mol.backbone.count > 0) {
            data->trajectory_data.secondary_structure.stride = data->mold.mol.backbone.count;
            data->trajectory_data.secondary_structure.count = data->mold.mol.backbone.count * data->mold.traj.num_frames;
            md_array_resize(data->trajectory_data.secondary_structure.data, data->mold.mol.backbone.count * data->mold.traj.num_frames, persistent_allocator);

            data->trajectory_data.backbone_angles.stride = data->mold.mol.backbone.count,
                data->trajectory_data.backbone_angles.count = data->mold.mol.backbone.count * data->mold.traj.num_frames,
                md_array_resize(data->trajectory_data.backbone_angles.data, data->mold.mol.backbone.count * data->mold.traj.num_frames, persistent_allocator);
        }
        
        // This is the number of slots we have to work with in parallel.
        // This number should ideally be more than the number of cores available.
        // We pre-allocate the number of slots * max frame data size, so don't go bananas here if you want to save some on memory.
        task_system::enqueue_pool("Preloading frames", 1, [data](task_system::TaskSetRange range)
            {
#define NUM_FRAME_SLOTS 64
                const int64_t slot_size = data->mold.traj.max_frame_data_size;
                void* slot_mem = md_alloc(default_allocator, slot_size * NUM_FRAME_SLOTS);

                void* slots[NUM_FRAME_SLOTS];
                atomic_queue::AtomicQueue<uint32_t, NUM_FRAME_SLOTS, 0xFFFFFFFF> slot_queue;

                for (uint32_t i = 0; i < NUM_FRAME_SLOTS; ++i) {
                    slots[i] = (char*)slot_mem + i * slot_size;
                    slot_queue.push(i);
                }

                // Iterate over all frames and load the raw data, then spawn a task for each frame
                for (uint32_t i = 0; i < data->mold.traj.num_frames; ++i) {
                    uint32_t slot_idx = slot_queue.pop();
                    ASSERT(slot_idx < NUM_FRAME_SLOTS);

                    void* frame_mem = slots[slot_idx];
                    const int64_t frame_size = data->mold.traj.extract_frame_data(data->mold.traj.inst, i, NULL);
                    data->mold.traj.extract_frame_data(data->mold.traj.inst, i, frame_mem);

                    // Spawn task: Load, Decode and postprocess
                    //printf("Spawning task to decode frame %i\n", i);
                    auto id = task_system::enqueue_pool("##Decode frame", 1, [slot_idx, frame_mem, frame_size, frame_idx = i, &slot_queue, slots, data](task_system::TaskSetRange)
                        {
                            md_frame_data_t* frame_data;
                            md_frame_cache_lock_t* lock;
                            if (md_frame_cache_reserve_frame(&data->mold.frame_cache, frame_idx, &frame_data, &lock)) {
                                md_trajectory_frame_header_t header;
                                data->mold.traj.decode_frame_header(data->mold.traj.inst, frame_mem, frame_size, &header);
                                data->mold.traj.decode_frame_coords(data->mold.traj.inst, frame_mem, frame_size, frame_data->x, frame_data->y, frame_data->z, frame_data->num_atoms);
                                // Free the data slot directly here
                                slot_queue.push(slot_idx);

                                memcpy(frame_data->box, header.box, sizeof(frame_data->box));
                                
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
                                memcpy(args.pbc.box, frame_data->box, sizeof(args.pbc.box));
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

                md_free(default_allocator, slot_mem, slot_size * NUM_FRAME_SLOTS);
#undef NUM_FRAME_SLOTS
            });
    }
}

static bool load_trajectory_data(ApplicationData* data, str_t filename) {
    interrupt_async_tasks(data);
    free_trajectory_data(data);
    data->files.trajectory = filename;
    data->animation.time = 0;
    data->animation.frame = 0;

    if (load::traj::open_file(&data->mold.traj, filename, &data->mold.mol, default_allocator)) {
        init_trajectory_data(data);
        return true;
    }

    return false;
}

static bool load_dataset_from_file(ApplicationData* data, str_t filename) {
    ASSERT(data);

    if (filename.len) {
        if (load::mol::is_extension_supported(filename)) {
            interrupt_async_tasks(data);
            free_molecule_data(data);
            free_trajectory_data(data);
            data->files.molecule = filename;
            data->files.trajectory = "";

            md_printf(MD_LOG_TYPE_INFO, "Attempting to load molecular data from file '%.*s'", filename.len, filename.ptr);
            if (!load::mol::load_file(&data->mold.mol, filename, default_allocator)) {
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
            return load_trajectory_data(data, filename);
        } else {
            md_print(MD_LOG_TYPE_ERROR, "File extension not supported");
        }
    }
    return false;
}

// ### WORKSPACE ###
static RepresentationType get_rep_type(str_t str) {
    if (compare_str_cstr(str, "VDW"))
        return RepresentationType::Vdw;
    else if (compare_str_cstr(str, "LICORICE"))
        return RepresentationType::Licorice;
    else if (compare_str_cstr(str, "BALL_AND_STICK"))    // Ball and stick is removed for now
        return RepresentationType::Vdw;
    else if (compare_str_cstr(str, "RIBBONS"))
        return RepresentationType::Ribbons;
    else if (compare_str_cstr(str, "CARTOON"))
        return RepresentationType::Cartoon;
    else
        return RepresentationType::Vdw;
}

static str_t get_rep_type_name(RepresentationType type) {
    switch (type) {
        case RepresentationType::Vdw:
            return make_cstr("VDW");
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

static void load_workspace(ApplicationData* data, str_t file) {
    ASSERT(data);
    clear_representations(data);
    //stats::remove_all_properties();
    //clear(data->density_volume.iso.isosurfaces);

    StrBuf<512> new_molecule_file;
    StrBuf<512> new_trajectory_file;

    str_t txt = load_textfile(file, frame_allocator);
    defer { free_str(txt, frame_allocator); };

    str_t c_txt = txt, line = {};
    while (extract_line(&line, &c_txt)) {
        if (compare_str_cstr_n(line, "[Files]", 7)) {
            while (c_txt.len && c_txt[0] != '[' && (extract_line(&line, &c_txt))) {
                if (compare_str_cstr_n(line, "MoleculeFile=", 13)) {
                    //std::filesystem::path path(file.ptr);
                    //path += trim_whitespace(substr(line, 13)).ptr;
                    //new_molecule_file = std::filesystem::canonical(path).string().c_str();
                    StrBuf<512> buf = file;
                    buf += trim_whitespace(substr(line, 13));
                    new_molecule_file = md_os_path_make_canonical(buf, default_temp_allocator);
                }
                if (compare_str_cstr_n(line, "TrajectoryFile=", 15)) {
                    //std::filesystem::path path(file.ptr);
                    //path += trim_whitespace(substr(line, 15)).ptr;
                    //new_trajectory_file = std::filesystem::canonical(path).string().c_str();
                    StrBuf<512> buf = file;
                    buf += trim_whitespace(substr(line, 15));
                    new_trajectory_file = md_os_path_make_canonical(buf, default_temp_allocator);
                }
            }
        } else if (compare_str_cstr_n(line, "[Representation]", 16)) {
            Representation* rep = create_representation(data);
            if (rep) {
                while (c_txt.len && c_txt[0] != '[' && (extract_line(&line, &c_txt))) {
                    if (compare_str_cstr_n(line, "Name=", 5)) rep->name = trim_whitespace(substr(line,5));
                    if (compare_str_cstr_n(line, "Filter=", 7)) rep->filter = trim_whitespace(substr(line,7));
                    if (compare_str_cstr_n(line, "Type=", 5)) rep->type = get_rep_type(trim_whitespace(substr(line,5)));
                    if (compare_str_cstr_n(line, "ColorMapping=", 13)) rep->color_mapping = get_color_mapping(trim_whitespace(substr(line,13)));
                    if (compare_str_cstr_n(line, "Enabled=", 8)) rep->enabled = parse_int(trim_whitespace(substr(line,8))) != 0;
                    if (compare_str_cstr_n(line, "StaticColor=", 12)) rep->uniform_color = parse_vec4(trim_whitespace(substr(line,12)));
                    if (compare_str_cstr_n(line, "Radius=", 7)) rep->radius = parse_float(trim_whitespace(substr(line,7)));
                    if (compare_str_cstr_n(line, "Tension=", 8)) rep->tension = parse_float(trim_whitespace(substr(line,8)));
                    if (compare_str_cstr_n(line, "Width=", 6)) rep->width = parse_float(trim_whitespace(substr(line,6)));
                    if (compare_str_cstr_n(line, "Thickness=", 10)) rep->thickness = parse_float(trim_whitespace(substr(line,10)));
                }
            }
        /*} else if (compare_str_cstr_n(line, "[Property]", 10)) {
            StrBuf<256> name, args;
            while (c_txt && c_txt[0] != '[' && (line = extract_line(c_txt))) {
                if (compare_str_cstr_n(line, "Name=", 5)) name = trim_whitespace(substr(line,5));
                if (compare_str_cstr_n(line, "Args=", 5)) args = trim_whitespace(substr(line,5));
            }
            //stats::create_property(name, args);
            */
        } else if (compare_str_cstr_n(line, "[RenderSettings]", 16)) {
            while (c_txt.len && c_txt[0] != '[' && (extract_line(&line, &c_txt))) {
                if (compare_str_cstr_n(line, "SsaoEnabled=", 12)) data->visuals.ssao.enabled = parse_int(trim_whitespace(substr(line,12))) != 0;
                if (compare_str_cstr_n(line, "SsaoIntensity=", 14)) data->visuals.ssao.intensity = parse_float(trim_whitespace(substr(line,14)));
                if (compare_str_cstr_n(line, "SsaoRadius=", 11)) data->visuals.ssao.radius = parse_float(trim_whitespace(substr(line,11)));
                if (compare_str_cstr_n(line, "SsaoBias=", 9)) data->visuals.ssao.bias = parse_float(trim_whitespace(substr(line,9)));
                if (compare_str_cstr_n(line, "DofEnabled=", 11)) data->visuals.dof.enabled = parse_int(trim_whitespace(substr(line,11))) != 0;
                if (compare_str_cstr_n(line, "DofFocusScale=", 14)) data->visuals.dof.focus_scale = parse_float(trim_whitespace(substr(line,14)));

                /*
                if (compare_str_cstr_n(line, "DensityVolumeEnabled=", 21)) data->density_volume.enabled = parse_int(trim_whitespace(substr(line,21))) != 0;
                if (compare_str_cstr_n(line, "DensityScale=", 13)) data->density_volume.dvr.density_scale = parse_float(trim_whitespace(substr(line,13)));
                if (compare_str_cstr_n(line, "AlphaScale=", 11)) data->density_volume.dvr.tf.alpha_scale = parse_float(trim_whitespace(substr(line,11)));
                if (compare_str_cstr_n(line, "TFFileName=", 11)) {
                    CStringView tfpath = trim_whitespace(substr(line,11));
                    if (!compare(tfpath, data->density_volume.dvr.tf.path)) {
                        data->density_volume.dvr.tf.path = tfpath;
                        data->density_volume.dvr.tf.dirty = true;
                    }
                }
                if (compare_str_cstr_n(line, "IsoSurfaceRenderingEnabled=", 27)) data->density_volume.iso.enabled = parse_int(trim_whitespace(substr(line,27))) != 0;
                */
            }
        /*} else if (compare_str_cstr_n(line, "[Isosurface]", 21)) {
            float isovalue = -1.0f;
            vec4_t isocolor{0.0f};
            while (c_txt && c_txt[0] != '[' && (line = extract_line(c_txt))) {
                if (compare_str_cstr_n(line, "Isovalue=", 9)) isovalue = parse_float(trim_whitespace(substr(line,9)));
                if (compare_str_cstr_n(line, "IsosurfaceColor=", 16)) isocolor = parse_vec4(trim_whitespace(substr(line,16)));
            }
            insert(data->density_volume.iso.isosurfaces, isovalue, isocolor);
            */
        } else if (compare_str_cstr_n(line, "[Camera]", 8)) {
            while (c_txt.len && c_txt[0] != '[' && (extract_line(&line, &c_txt))) {
                if (compare_str_cstr_n(line, "Position=", 9)) {
                    vec3_t pos = vec3_from_vec4(parse_vec4(trim_whitespace(substr(line,9))));
                    data->view.camera.position = pos;
                    data->view.animation.target_position = pos;
                }
                if (compare_str_cstr_n(line, "Rotation=", 9)) {
                    vec4_t v = parse_vec4(trim_whitespace(substr(line, 9,-1)));
                    data->view.camera.orientation.x = v.x;
                    data->view.camera.orientation.y = v.y;
                    data->view.camera.orientation.z = v.z;
                    data->view.camera.orientation.w = v.w;
                }
                if (compare_str_cstr_n(line, "Distance=", 9)) {
                    data->view.camera.focus_distance = parse_float(trim_whitespace(substr(line,9)));
                }
            }
        }
    }

    data->files.workspace = file;

    if (!compare_str(new_molecule_file, data->files.molecule) && new_molecule_file) {
        load_dataset_from_file(data, new_molecule_file);
    }

    if (!compare_str(new_trajectory_file, data->files.trajectory) && new_trajectory_file) {
        load_dataset_from_file(data, new_trajectory_file);
    }

    data->animation = {};
    reset_view(data, false, true);
    init_all_representations(data);
    update_all_representations(data);
}

static void save_workspace(ApplicationData* data, str_t filename) {
    md_file_o* file = md_file_open(filename, MD_FILE_WRITE);
    if (!file) {
        md_printf(MD_LOG_TYPE_ERROR, "ERROR! Could not save workspace to file '%.*s'\n", (int)filename.len, filename.ptr);
        return;
    }

    FILE* fptr = (FILE*)file;

    str_t rel_mol_path = {};
    str_t rel_traj_path = {};

    if (data->files.molecule) {
        //std::filesystem::path basep(filename.ptr);
        //std::filesystem::path filep(data->files.molecule.cstr());
        //std::filesystem::path relp = filep.lexically_relative(basep);
        //rel_mol_path = relp.string();
        rel_mol_path = md_os_path_make_relative(filename, data->files.molecule, default_temp_allocator);
    }

    if (data->files.trajectory) {
        //std::filesystem::path basep(filename.ptr);
        //std::filesystem::path filep(data->files.trajectory.cstr());
        //std::filesystem::path relp = filep.lexically_relative(basep);
        //rel_traj_path = relp.string();
        rel_traj_path = md_os_path_make_relative(filename, data->files.trajectory, default_temp_allocator);
    }

    fprintf(fptr, "[Files]\n");
    fprintf(fptr, "MoleculeFile=%.*s\n", (int)rel_mol_path.len, rel_mol_path.ptr);
    fprintf(fptr, "TrajectoryFile=%.*s\n", (int)rel_traj_path.len, rel_traj_path.ptr);
    fprintf(fptr, "\n");

    // REPRESENTATIONS
    for (int64_t i = 0; i < md_array_size(data->representations.buffer); ++i) {
        auto& rep = data->representations.buffer[i];
        fprintf(fptr, "[Representation]\n");
        fprintf(fptr, "Name=%s\n", rep.name.cstr());
        fprintf(fptr, "Filter=%s\n", rep.filter.cstr());
        fprintf(fptr, "Type=%s\n", get_rep_type_name(rep.type).ptr);
        fprintf(fptr, "ColorMapping=%s\n", get_color_mapping_name(rep.color_mapping).ptr);
        fprintf(fptr, "Enabled=%i\n", (int)rep.enabled);
        fprintf(fptr, "StaticColor=%g,%g,%g,%g\n", rep.uniform_color.x, rep.uniform_color.y, rep.uniform_color.z, rep.uniform_color.w);
        fprintf(fptr, "Radius=%g\n", rep.radius);
        fprintf(fptr, "Tension=%g\n", rep.tension);
        fprintf(fptr, "Width=%g\n", rep.width);
        fprintf(fptr, "Thickness=%g\n", rep.thickness);
        fprintf(fptr, "\n");
    }

    // PROPERTIES
    /*
    for (const auto prop : stats::get_properties()) {
        fprintf(fptr, "[Property]\n");
        fprintf(fptr, "Name=%s\n", prop->name_buf.cstr());
        fprintf(fptr, "Args=%s\n", prop->args_buf.cstr());
        fprintf(fptr, "\n");
    }
    */

    fprintf(fptr, "[RenderSettings]\n");
    fprintf(fptr, "SsaoEnabled=%i\n", data->visuals.ssao.enabled ? 1 : 0);
    fprintf(fptr, "SsaoIntensity=%g\n", data->visuals.ssao.intensity);
    fprintf(fptr, "SsaoRadius=%g\n", data->visuals.ssao.radius);
    fprintf(fptr, "SsaoBias=%g\n", data->visuals.ssao.bias);
    fprintf(fptr, "\n");

    fprintf(fptr, "DofEnabled=%i\n", data->visuals.dof.enabled ? 1 : 0);
    fprintf(fptr, "DofFocusScale=%g\n", data->visuals.dof.focus_scale);
    fprintf(fptr, "\n");

    /*
    fprintf(fptr, "DensityVolumeEnabled=%i\n", data->density_volume.enabled ? 1 : 0);
    fprintf(fptr, "DensityScale=%g\n", data->density_volume.dvr.density_scale);
    fprintf(fptr, "AlphaScale=%g\n", data->density_volume.dvr.tf.alpha_scale);
    fprintf(fptr, "TFFileName=%s\n", data->density_volume.dvr.tf.path.cstr());
    fprintf(fptr, "IsoSurfaceRenderingEnabled=%i\n", data->density_volume.iso.enabled ? 1 : 0);
    fprintf(fptr, "\n");

    for (int i = 0; i < data->density_volume.iso.isosurfaces.count; i++) {
        const auto value = data->density_volume.iso.isosurfaces.values[i];
        const auto color = data->density_volume.iso.isosurfaces.colors[i];
        fprintf(fptr, "[Isosurface]\n");
        fprintf(fptr, "Isovalue=%g\n", value);
        fprintf(fptr, "IsosurfaceColor=%g,%g,%g,%g\n", color.x, color.y, color.z, color.w);
    }
    */

    fprintf(fptr, "[Camera]\n");
    fprintf(fptr, "Position=%g,%g,%g\n", data->view.camera.position.x, data->view.camera.position.y, data->view.camera.position.z);
    fprintf(fptr, "Rotation=%g,%g,%g,%g\n", data->view.camera.orientation.x, data->view.camera.orientation.y, data->view.camera.orientation.z,
            data->view.camera.orientation.w);
    fprintf(fptr, "Distance=%g\n", data->view.camera.focus_distance);
    fprintf(fptr, "\n");

    fclose(fptr);

    data->files.workspace = filename;
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
    rep.filter = filter;   
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
        md_bitfield_or(&mask, &mask, &rep.atom_mask);
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
        default:
            ASSERT(false);
            break;
    }

    rep->filter_is_valid = filter_expression(*data, rep->filter, &rep->atom_mask);
    filter_colors(colors, mol.atom.count, &rep->atom_mask);

    data->representations.atom_visibility_mask_dirty = true;

    {
        md_gl_representation_type_t type = MD_GL_REP_DEFAULT;
        md_gl_representation_args_t args = {};
        switch(rep->type) {
        case RepresentationType::Vdw:
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
    md_bitfield_init(&rep->atom_mask, default_allocator);
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

    static RegionMode region_mode = RegionMode::Append;
    static bool region_select = false;
    static application::Coordinate x0;
    const application::Coordinate x1 = data->ctx.input.mouse.win_coord;
    const int64_t N = data->mold.mol.atom.count;
    const bool shift_down = data->ctx.input.key.down[Key::KEY_LEFT_SHIFT] || data->ctx.input.key.down[Key::KEY_RIGHT_SHIFT];
    const bool mouse_down = data->ctx.input.mouse.down[0] || data->ctx.input.mouse.down[1];

    md_exp_bitfield_t mask;
    md_bitfield_init(&mask, frame_allocator);
    defer { md_bitfield_free(&mask); };

    md_bitfield_clear(&data->selection.current_highlight_mask);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;

    if (data->picking.idx != INVALID_PICKING_IDX && !region_select) {
        ASSERT(0 <= data->picking.idx && data->picking.idx <= N);
        md_bitfield_set_bit(&mask, data->picking.idx);

        switch (data->selection.level_mode) {
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

            switch (data->selection.level_mode) {
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
                md_bitfield_or(dst_mask, src_mask, &mask);
            } else if (region_mode == RegionMode::Remove) {
                md_bitfield_not(&mask, &mask, 0, N);
                md_bitfield_and(dst_mask, src_mask, &mask);
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
                    md_bitfield_or(&data->selection.current_selection_mask, &data->selection.current_selection_mask, &mask);
                } else {
                    md_bitfield_not(&mask, &mask, 0, N);
                    md_bitfield_and(&data->selection.current_selection_mask, &data->selection.current_selection_mask, &mask);
                }
            } else if (data->ctx.input.mouse.clicked[1]) {
                md_bitfield_clear(&data->selection.current_selection_mask);
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
