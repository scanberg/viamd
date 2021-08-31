//#include <core/types.h>
//#include <core/file.h>
//#include <core/hash.h>
//#include <core/log.h>
//#include <core/bitfield.h>
//#include <core/math_utils.h>
//#include <core/string_utils.h>
#include <core/spatial_hash.h>
//#include <core/sync.h>

#include <md_util.h>
#include <md_gl.h>
#include <md_filter.h>
#include <md_script.h>
#include <md_molecule.h>
#include <md_trajectory.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_simd.h>

#include <imgui.h>
#define IMGUI_DEFINE_MATH_OPERATORS
#include <imgui_internal.h>
#include <implot.h>
#include <TextEditor.h>

#include "gfx/gl.h"
#include "gfx/camera.h"
#include "gfx/camera_utils.h"
//#include "gfx/molecule_draw.h"
#include "gfx/immediate_draw_utils.h"
#include "gfx/postprocessing_utils.h"
#include "gfx/volumerender_utils.h"
#include "gfx/conetracing_utils.h"

#include "volume.h"
#include "range_slider.h"
#include "plot_extended.h"
#include "platform/platform.h"
#include "console.h"
#include "ramachandran.h"
#include "color_utils.h"
#include "isosurface.h"
#include "task_system.h"
#include "loader.h"

#include <stdio.h>

#define PICKING_JITTER_HACK 1
#define DEPERIODIZE_ON_LOAD 1
#define SHOW_IMGUI_DEMO_WINDOW 1
#define VIAMD_RELEASE 0
#define EXPERIMENTAL_CONE_TRACED_AO 0
#define USE_MOLD 1

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

constexpr CStringView FILE_EXTENSION = "via"; 
constexpr u32 NO_PICKING_IDX = ~0U;

constexpr u32 DEL_BTN_COLOR = 0xFF1111CC;
constexpr u32 DEL_BTN_HOVER_COLOR = 0xFF3333DD;
constexpr u32 DEL_BTN_ACTIVE_COLOR = 0xFF5555FF;
constexpr u32 TEXT_BG_ERROR_COLOR = 0xAA222299;

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

inline const ImVec4& vec_cast(const vec4& v) { return *(const ImVec4*)(&v); }
inline const vec4& vec_cast(const ImVec4& v) { return *(const vec4*)(&v); }
inline const ImVec2& vec_cast(const vec2& v) { return *(const ImVec2*)(&v); }
inline const vec2& vec_cast(const ImVec2& v) { return *(const vec2*)(&v); }

inline ImVec4& vec_cast(vec4& v) { return *(ImVec4*)(&v); }
inline vec4& vec_cast(ImVec4& v) { return *(vec4*)(&v); }
inline ImVec2& vec_cast(vec2& v) { return *(ImVec2*)(&v); }
inline vec2& vec_cast(ImVec2& v) { return *(vec2*)(&v); }

inline str_t str_cast(CStringView view) { return {.ptr = view.cstr(), .len = view.length() }; }

static inline bool operator == (const ImVec2& lhs, const ImVec2& rhs) { return lhs.x == rhs.x && lhs.y == rhs.y; }
static inline bool operator != (const ImVec2& lhs, const ImVec2& rhs) { return !(lhs == rhs); }

static inline bool operator == (const ImVec4& lhs, const ImVec4& rhs) { return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w; }
static inline bool operator != (const ImVec4& lhs, const ImVec4& rhs) { return !(lhs == rhs); }

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
    AtomBit_Highlighted = BIT(0),
    AtomBit_Selected    = BIT(1),
    AtomBit_Visible     = BIT(2)
};

enum MolBit_ {
    MolBit_DirtyPosition            = BIT(0),
    MolBit_DirtySecondaryStructure  = BIT(1),
    MolBit_DirtyFlags               = BIT(2)
};

enum RepBit_ {
    RepBit_DirtyColor   = BIT(0),
    RepBit_DirtyFilter  = BIT(1)
};

// #struct Structure Declarations

struct PickingData {
    u32 idx = NO_PICKING_IDX;
    f32 depth = 1.0f;
    vec3 world_coord = {0, 0, 0};
    vec2 screen_coord = {0, 0};
};

struct MainFramebuffer {
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
    } pbo_picking;

    int width = 0;
    int height = 0;
};

struct Representation {
    StringBuffer<64> name = "rep";
    StringBuffer<256> filter = "all";
    RepresentationType type = RepresentationType::Vdw;
    ColorMapping color_mapping = ColorMapping::Cpk;
    Bitfield atom_mask{};
    md_gl_representation md_rep{};

    bool enabled = true;
    bool show_in_selection = true;
    bool filter_is_valid = true;
    bool filter_is_dynamic = false;

    uint32_t flags = 0;

    // User defined color used in uniform mode
    vec4 uniform_color = vec4(1);

    // VDW and Ball & Stick
    f32 radius = 1.f;

    // Ball & Stick and Licorice, Ribbons, Cartoon
    f32 thickness = 1.f;

    // Ribbons, Cartoon
    f32 tension = 0.5f;
    f32 width = 1.f;
};

struct Selection {
    StringBuffer<64> name = "sel";
    Bitfield atom_mask{};
};

struct ApplicationData {
    // --- PLATFORM ---
    platform::Context ctx;

    // --- FILES ---
    // for keeping track of open files
    struct {
        StringBuffer<512> molecule{};
        StringBuffer<512> trajectory{};
        StringBuffer<512> workspace{};
    } files;

    // --- CAMERA ---
    struct {
        Camera camera{};
        TrackballControllerParam trackball_param;
        ViewParam param{};
        CameraMode mode = CameraMode::Perspective;

        struct {
            vec2 sequence[16];
        } jitter;

        struct {
            vec3 target_position{};
            quat target_orientation{};
        } animation;
    } view;

    // --- MOLD DATA ---
    struct {
        md_gl_context       gl_ctx = {0};
        md_gl_molecule      gl_mol = {0};
        md_molecule_t       mol = {0};
        md_trajectory_i     traj = {0};
        struct {
            md_script_ir_t    ir = {0};
            md_script_eval_result_t eval = {0};
        } script;
        uint32_t dirty_buffers = {0};
    } mold;

    // --- ASYNC TASKS HANDLES ---
    struct {
        task_system::ID load_trajectory;
    } tasks;

    // --- ATOM SELECTION ---
    struct {
        // bool show_window = false;
        SelectionLevel level_mode = SelectionLevel::Atom;
        SelectionOperator op_mode = SelectionOperator::Or;
        SelectionGrowth grow_mode = SelectionGrowth::CovalentBond;

        i32 hovered = -1;
        i32 right_clicked = -1;

        Bitfield current_selection_mask{};
        Bitfield current_highlight_mask{};
        DynamicArray<Selection> stored_selections{};

        struct {
            struct {
                vec4 fill_color = vec4(1.0f, 1.0f, 1.0f, 0.5f);
                vec4 outline_color = vec4(1.0f, 0.5f, 0.0f, 0.5f);
                float outline_scale = 1.1f;
            } highlight;

            struct {
                vec4 fill_color = vec4(1.0f, 1.0f, 1.0f, 0.5f);
                vec4 outline_color = vec4(0.0f, 0.5f, 1.0f, 0.5f);
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
    MainFramebuffer fbo;

    PickingData picking;

    // --- ANIMATION ---
    struct {
        f64 time = 0.f;  // double precision for long trajectories
        i32 frame = 0;
        f32 fps = 10.f;
        InterpolationMode interpolation = InterpolationMode::Cubic;
        PlaybackMode mode = PlaybackMode::Stopped;
        bool apply_pbc = false;
    } animation;

    // --- TIME LINE FILTERING ---
    struct {
        bool enabled = true;
        Range<float> range{0, 0};
        bool dynamic_window = false;
        float window_extent = 10.f;
    } time_filter;

    // --- VISUALS ---
    struct {
        struct {
            vec3 color = {1, 1, 1};
            float intensity = 24.f;
        } background;

        struct {
            bool enabled = true;
            f32 intensity = 6.0f;
            f32 radius = 6.0f;
            f32 bias = 0.1f;
        } ssao;

#if EXPERIMENTAL_CONE_TRACED_AO == 1
        struct {
            bool enabled = true;
            f32 intensity = 1.0f;
            f32 step_scale = 1.0f;
        } cone_traced_ao;
#endif

        struct {
            bool enabled = false;
            struct {
                f32 target = 5.0f;
                f32 current = 5.0f;
            } focus_depth;
            f32 focus_scale = 10.0f;
        } dof;

        struct {
            bool enabled = true;
            bool jitter = true;
            f32 feedback_min = 0.88f;
            f32 feedback_max = 0.97f;

            struct {
                bool enabled = true;
                f32 motion_scale = 1.0f;
            } motion_blur;
        } temporal_reprojection;

        struct {
            bool enabled = true;
            postprocessing::Tonemapping tonemapper = postprocessing::Tonemapping_Filmic;
            f32 exposure = 1.f;
            f32 gamma = 2.2;
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
        vec4 color = vec4(1, 0, 1, 1);
        f32 distance_cutoff = HYDROGEN_BOND_DISTANCE_CUTOFF_DEFAULT;  // In Ångström
        f32 angle_cutoff = HYDROGEN_BOND_ANGLE_CUTOFF_DEFAULT;        // In Degrees
        DynamicArray<HydrogenBond> bonds{};
    } hydrogen_bonds;
    */

    struct {
        bool enabled = false;
        vec4 color = vec4(0, 0, 0, 0.5);
        mat3 box = mat3(0);
    } simulation_box;

    struct {
        bool show_window = false;
        bool enabled = false;

        struct {
            bool enabled = true;
            f32 density_scale = 1.f;
            struct {
                GLuint id = 0;
                bool dirty = true;
                int width = 0;
                f32 alpha_scale = 1.f;
                StringBuffer<512> path;
            } tf;
        } dvr;

        struct {
            bool enabled = false;
            IsoSurfaces isosurfaces;
        } iso;

        struct {
            GLuint id = 0;
            bool dirty = false;
            ivec3 dim = ivec3(0);
            f32 max_value = 1.f;
        } volume_texture;

        struct {
            GLuint id = 0;
            ivec2 dim = {0,0};
        } render_texture;

        struct {
            vec3 min = {0, 0, 0};
            vec3 max = {1, 1, 1};
        } clip_volume;

        // mat4 model_to_world_matrix{};
        // mat4 world_to_model_matrix{};
        vec3 voxel_spacing{1.0f};
        float resolution_scale = 2.0f;

    } density_volume;

#if EXPERIMENTAL_CONE_TRACED_AO == 1
    struct {
        cone_trace::GPUVolume vol;
        std::mutex vol_mutex{};
    } occupancy_volume;
#endif

    // --- RAMACHANDRAN ---
    struct {
        bool show_window = false;

        ramachandran::ColorMap color_map{};

        struct {
            bool enabled = false;
            f32 radius = 0.2f;
            vec4 color = vec4(0, 0, 0, 1);
        } range;

        struct {
            bool enabled = true;
            vec4 border_color = vec4(0, 0, 0, 1);
            struct {
                f32 radius = 2.5f;
                vec4 fill_color = vec4(0.75f, 0.75f, 0.75f, 1);
            } base;
            struct {
                f32 radius = 3.5f;
                vec4 selection_color = vec4(0.8f, 1.0f, 0.5f, 1.0f);
                vec4 highlight_color = vec4(0.8f, 1.0f, 0.5f, 0.5f);
            } selection;
        } current;
    } ramachandran;

    // --- REPRESENTATIONS ---
    struct {
        DynamicArray<Representation> buffer = {};
        Bitfield atom_visibility_mask{};
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

// custom ImGui procedures
namespace ImGui {
static bool DeleteButton(const char* label, const ImVec2& size = ImVec2(0, 0));
static void CreateDockspace();
static void BeginCanvas(const char* id);
static void EndCanvas();
static void PushDisabled() {
    ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
}
static void PopDisabled() {
    ImGui::PopItemFlag();
    ImGui::PopStyleVar();
}

void init_theme() {
    ImVec4* colors = ImGui::GetStyle().Colors;
    colors[ImGuiCol_Text]                   = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
    colors[ImGuiCol_TextDisabled]           = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
    colors[ImGuiCol_WindowBg]               = ImVec4(0.00f, 0.00f, 0.00f, 0.67f);
    colors[ImGuiCol_ChildBg]                = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    colors[ImGuiCol_PopupBg]                = ImVec4(0.11f, 0.11f, 0.14f, 0.92f);
    colors[ImGuiCol_Border]                 = ImVec4(0.50f, 0.50f, 0.50f, 0.50f);
    colors[ImGuiCol_BorderShadow]           = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    colors[ImGuiCol_FrameBg]                = ImVec4(0.43f, 0.43f, 0.43f, 0.39f);
    colors[ImGuiCol_FrameBgHovered]         = ImVec4(0.70f, 0.70f, 0.70f, 0.40f);
    colors[ImGuiCol_FrameBgActive]          = ImVec4(0.65f, 0.65f, 0.65f, 0.69f);
    colors[ImGuiCol_TitleBg]                = ImVec4(0.07f, 0.07f, 0.07f, 0.83f);
    colors[ImGuiCol_TitleBgActive]          = ImVec4(0.00f, 0.00f, 0.00f, 0.87f);
    colors[ImGuiCol_TitleBgCollapsed]       = ImVec4(0.00f, 0.00f, 0.00f, 0.20f);
    colors[ImGuiCol_MenuBarBg]              = ImVec4(0.55f, 0.55f, 0.55f, 0.80f);
    colors[ImGuiCol_ScrollbarBg]            = ImVec4(0.21f, 0.21f, 0.21f, 0.60f);
    colors[ImGuiCol_ScrollbarGrab]          = ImVec4(0.83f, 0.83f, 0.83f, 0.30f);
    colors[ImGuiCol_ScrollbarGrabHovered]   = ImVec4(0.83f, 0.83f, 0.83f, 0.40f);
    colors[ImGuiCol_ScrollbarGrabActive]    = ImVec4(0.81f, 0.81f, 0.81f, 0.60f);
    colors[ImGuiCol_CheckMark]              = ImVec4(0.90f, 0.90f, 0.90f, 0.50f);
    colors[ImGuiCol_SliderGrab]             = ImVec4(1.00f, 1.00f, 1.00f, 0.30f);
    colors[ImGuiCol_SliderGrabActive]       = ImVec4(0.82f, 0.82f, 0.82f, 0.60f);
    colors[ImGuiCol_Button]                 = ImVec4(0.64f, 0.64f, 0.64f, 0.62f);
    colors[ImGuiCol_ButtonHovered]          = ImVec4(0.72f, 0.72f, 0.72f, 0.79f);
    colors[ImGuiCol_ButtonActive]           = ImVec4(0.80f, 0.80f, 0.81f, 0.85f);
    colors[ImGuiCol_Header]                 = ImVec4(0.64f, 0.64f, 0.64f, 0.45f);
    colors[ImGuiCol_HeaderHovered]          = ImVec4(0.65f, 0.65f, 0.65f, 0.80f);
    colors[ImGuiCol_HeaderActive]           = ImVec4(0.85f, 0.85f, 0.85f, 0.80f);
    colors[ImGuiCol_Separator]              = ImVec4(0.50f, 0.50f, 0.50f, 1.00f);
    colors[ImGuiCol_SeparatorHovered]       = ImVec4(0.71f, 0.71f, 0.71f, 1.00f);
    colors[ImGuiCol_SeparatorActive]        = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
    colors[ImGuiCol_ResizeGrip]             = ImVec4(1.00f, 1.00f, 1.00f, 0.16f);
    colors[ImGuiCol_ResizeGripHovered]      = ImVec4(1.00f, 1.00f, 1.00f, 0.47f);
    colors[ImGuiCol_ResizeGripActive]       = ImVec4(1.00f, 0.96f, 1.00f, 0.63f);
    colors[ImGuiCol_Tab]                    = ImVec4(0.56f, 0.56f, 0.56f, 0.78f);
    colors[ImGuiCol_TabHovered]             = ImVec4(0.87f, 0.87f, 0.87f, 0.80f);
    colors[ImGuiCol_TabActive]              = ImVec4(0.73f, 0.73f, 0.73f, 0.84f);
    colors[ImGuiCol_TabUnfocused]           = ImVec4(0.57f, 0.57f, 0.57f, 0.82f);
    colors[ImGuiCol_TabUnfocusedActive]     = ImVec4(0.65f, 0.65f, 0.65f, 0.84f);
    colors[ImGuiCol_DockingPreview]         = ImVec4(0.90f, 0.90f, 0.90f, 0.31f);
    colors[ImGuiCol_DockingEmptyBg]         = ImVec4(0.20f, 0.20f, 0.20f, 1.00f);
    colors[ImGuiCol_PlotLines]              = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
    colors[ImGuiCol_PlotLinesHovered]       = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
    colors[ImGuiCol_PlotHistogram]          = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
    colors[ImGuiCol_PlotHistogramHovered]   = ImVec4(1.00f, 0.60f, 0.00f, 1.00f);
    colors[ImGuiCol_TextSelectedBg]         = ImVec4(0.00f, 0.00f, 1.00f, 0.35f);
    colors[ImGuiCol_DragDropTarget]         = ImVec4(1.00f, 1.00f, 0.00f, 0.90f);
    colors[ImGuiCol_NavHighlight]           = ImVec4(0.45f, 0.45f, 0.90f, 0.80f);
    colors[ImGuiCol_NavWindowingHighlight]  = ImVec4(1.00f, 1.00f, 1.00f, 0.70f);
    colors[ImGuiCol_NavWindowingDimBg]      = ImVec4(0.80f, 0.80f, 0.80f, 0.20f);
    colors[ImGuiCol_ModalWindowDimBg]       = ImVec4(0.20f, 0.20f, 0.20f, 0.35f);
}

}  // namespace ImGui

static void interpolate_atomic_properties(ApplicationData* data);
static void update_view_param(ApplicationData* data);
static void reset_view(ApplicationData* data, bool move_camera = false, bool smooth_transition = false);
static f32 compute_avg_ms(f32 dt);
static PickingData read_picking_data(const MainFramebuffer& fbo, i32 x, i32 y);

static void compute_aabb(vec3* aabb_min, vec3* aabb_max, const float* x, const float* y, const float* z, int64_t count);

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
static void draw_representations_lean_and_mean(ApplicationData* data, u32 mask = 0xFFFFFFFFU);

static void draw_main_menu(ApplicationData* data);
static void draw_context_popup(ApplicationData* data);
static void draw_animation_control_window(ApplicationData* data);
static void draw_representations_window(ApplicationData* data);
//static void draw_property_window(ApplicationData* data);
static void draw_timeline_window(ApplicationData* data);
static void draw_distribution_window(ApplicationData* data);
static void draw_ramachandran_window(ApplicationData* data);
static void draw_atom_info_window(const ApplicationData& data);
static void draw_molecule_dynamic_info_window(ApplicationData* data);
static void draw_async_task_window(ApplicationData* data);
//static void draw_reference_frame_window(ApplicationData* data);
//static void draw_shape_space_window(ApplicationData* data);
static void draw_density_volume_window(ApplicationData* data);
static void draw_script_editor_window(ApplicationData* data);
static void draw_dataset_window(ApplicationData* data);
// static void draw_density_volume_clip_plane_widgets(ApplicationData* data);
// static void draw_selection_window(ApplicationData* data);

static void init_framebuffer(MainFramebuffer* fbo, int width, int height);
static void destroy_framebuffer(MainFramebuffer* fbo);

static void update_md_buffers(ApplicationData* data);

static void init_molecule_data(ApplicationData* data);
static void init_trajectory_data(ApplicationData* data);

static bool load_dataset_from_file(ApplicationData* data, CStringView file);
//static void free_dataset(ApplicationData* data);

static void load_workspace(ApplicationData* data, CStringView file);
static void save_workspace(ApplicationData* data, CStringView file);

static void create_screenshot(ApplicationData* data);

// Representations
static Representation* create_representation(ApplicationData* data, RepresentationType type = RepresentationType::Vdw,
                                             ColorMapping color_mapping = ColorMapping::Cpk, CStringView filter = "all");
static Representation* clone_representation(ApplicationData* data, const Representation& rep);
static void remove_representation(ApplicationData* data, int idx);
static void update_representation(ApplicationData* data, Representation* rep);
static void update_all_representations(ApplicationData* data);
static void init_representation(ApplicationData* data, Representation* rep);
static void init_all_representations(ApplicationData* data);
static void reset_representations(ApplicationData* data);
static void clear_representations(ApplicationData* data);

static void recompute_atom_visibility_mask(ApplicationData* data);

// Selections
static Selection* create_selection(ApplicationData* data, CStringView name, const Bitfield);
static Selection* clone_selection(ApplicationData* data, const Selection& sel);
static void remove_selection(ApplicationData* data, int idx);

static void reset_selections(ApplicationData* data);
static void clear_selections(ApplicationData* data);

static bool filter_expression(const ApplicationData& data, CStringView expr, Bitfield mask);

static void modify_field(Bitfield field, const Bitfield mask, SelectionOperator op) {
    switch(op) {
    case SelectionOperator::Or:
        bitfield::or_field(field, field, mask);
        break;
    case SelectionOperator::And:
        bitfield::and_field(field, field, mask);
        break;
    case SelectionOperator::Replace:
        bitfield::copy(field, mask);
        break;
    default:
        ASSERT(false);
    }
}

template <typename Int>
static void modify_field(Bitfield field, Range<Int> range, SelectionOperator op) {
    switch(op) {
    case SelectionOperator::Or:
        bitfield::set_range(field, range);
        break;
    case SelectionOperator::And:
        bitfield::clear_range(field, Range<Int>(0, range.beg));
        bitfield::clear_range(field, Range<Int>(range.end, (int)field.size()));
        break;
    case SelectionOperator::Replace:
        bitfield::clear_all(field);
        bitfield::set_range(field, range);
        break;
    case SelectionOperator::Clear:
        bitfield::clear_range(field, range);
        break;
    default:
        ASSERT(false);
    }
}

static void clear_highlight(ApplicationData* data) {
    ASSERT(data);
    bitfield::clear_all(data->selection.current_highlight_mask);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

static void modify_highlight(ApplicationData* data, const Bitfield atom_mask, SelectionOperator op = SelectionOperator::Replace) {
    ASSERT(data);
    modify_field(data->selection.current_highlight_mask, atom_mask, op);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

template <typename Int>
static void modify_highlight(ApplicationData* data, const Range<Int> range, SelectionOperator op = SelectionOperator::Replace) {
    ASSERT(data);
    modify_field(data->selection.current_highlight_mask, range, op);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

static void clear_selection(ApplicationData* data) {
    ASSERT(data);
    bitfield::clear_all(data->selection.current_selection_mask);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

static void modify_selection(ApplicationData* data, const Bitfield atom_mask, SelectionOperator op = SelectionOperator::Replace) {
    ASSERT(data);
    modify_field(data->selection.current_selection_mask, atom_mask, op);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

template <typename Int>
static void modify_selection(ApplicationData* data, const Range<Int> range, SelectionOperator op = SelectionOperator::Replace) {
    ASSERT(data);
    modify_field(data->selection.current_selection_mask, range, op);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

static void on_trajectory_load_complete(ApplicationData* data);

// Global allocators for application
static md_allocator_i* frame_allocator = 0;
static md_allocator_i* persistent_allocator = default_allocator;

int main(int, char**) {
    frame_allocator = md_arena_allocator_create(default_allocator, MEGABYTES(8));

    ApplicationData data;

    // Init logging
    logging::initialize();
    logging::register_backend([](CStringView str, logging::Severity, void*) {
        print_string(str);
        printf("\n");
    });
    logging::register_backend(
        [](CStringView str, logging::Severity severity, void* usr_data) {
            const char* modifier = "";
            switch (severity) {
                case logging::Severity::Note:
                    modifier = "[note] ";
                    break;
                case logging::Severity::Warning:
                    modifier = "[warning] ";
                    break;
                case logging::Severity::Error:
                    modifier = "[error] ";
                    break;
                case logging::Severity::Fatal:
                    modifier = "[fatal] ";
                    break;
                default:
                    break;
            }
            ((Console*)usr_data)->AddLog("%s%s", modifier, str.cstr());
        },
        &data.console);

    // Init platform
    LOG_NOTE("Initializing GL...");
    if (!platform::initialize(&data.ctx, 1920, 1080, "VIAMD")) {
        LOG_ERROR("Could not initialize platform layer... terminating\n");
        return -1;
    }
    data.ctx.window.vsync = true;

    LOG_NOTE("Creating framebuffer...");
    init_framebuffer(&data.fbo, data.ctx.framebuffer.width, data.ctx.framebuffer.height);

    math::generate_halton_sequence(data.view.jitter.sequence, ARRAY_SIZE(data.view.jitter.sequence), 2, 3);

    // Init subsystems
    LOG_NOTE("Initializing immediate draw...");
    immediate::initialize();
    LOG_NOTE("Initializing ramachandran...");
    ramachandran::initialize();
    LOG_NOTE("Initializing post processing...");
    postprocessing::initialize(data.fbo.width, data.fbo.height);
    LOG_NOTE("Initializing volume...");
    volume::initialize();
#if EXPERIMENTAL_CONE_TRACED_AO == 1
    LOG_NOTE("Initializing cone tracing...");
    cone_trace::initialize();
#endif
    LOG_NOTE("Initializing task system...");
    task_system::initialize();

    md_gl_context_init(&data.mold.gl_ctx);

    ImGui::init_theme();

    data.script.editor.SetLanguageDefinition(TextEditor::LanguageDefinition::VIAMD());
    data.script.editor.SetPalette(TextEditor::GetRetroBluePalette());

#if VIAMD_RELEASE
    load::mol::load_string_pdb(&data.mold.mol, CAFFINE_PDB, default_allocator);
    init_molecule_data(&data);
    create_representation(&data);
#else
    //load_dataset_from_file(&data, "D:/data/1aon.pdb");
    load_dataset_from_file(&data, "D:/data/md/alanine/1ALA-560ns.pdb");
    //load_dataset_from_file(&data,  "D:/data/md/amyloid-PFT/centered.gro");
    // load_dataset_from_file(&data, VIAMD_DATA_DIR "/amyloid/centered.xtc");
    // load_dataset_from_file(&data, "D:/data/md/6T-water/16-6T-box.gro");
    // load_dataset_from_file(&data, "D:/data/md/6T-water/16-6T-box-md-npt.xtc");
    // load_dataset_from_file(&data, "D:/data/md/p-ftaa-water/p-FTAA-box-sol-ions-em.gro");
    // load_dataset_from_file(&data, "D:/data/md/p-ftaa-water/p-FTAA-box-sol-ions-md-npt.xtc");

    create_representation(&data, RepresentationType::Vdw, ColorMapping::Cpk, "protein");
    data.script.editor.SetText("d1 = distance(1,2);\na1 = angle(1,2,3) in resname(\"ALA\");");
#endif
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
        platform::update(&data.ctx);

#if SHOW_IMGUI_DEMO_WINDOW
        ImGui::ShowDemoWindow();
        ImPlot::ShowDemoWindow();
#endif

        const i32 num_frames = (i32)traj.num_frames;
        const i32 last_frame = math::max(0, num_frames - 1);
        const f64 max_time = (f64)math::max(0, last_frame);

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
                LOG_NOTE("Recompiling shaders and re-initializing volume");
                postprocessing::initialize(data.fbo.width, data.fbo.height);
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
                f64 step = data.ctx.input.key.down[Key::KEY_LEFT_CONTROL] ? 10.0 : 1.0;
                if (data.ctx.input.key.hit[KEY_SKIP_TO_PREV_FRAME]) step = -step;
                data.animation.time = math::clamp(data.animation.time + step, 0.0, max_time);
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
            data.animation.time = math::clamp(data.animation.time, 0.0, max_time);
            if (data.animation.time >= max_time) {
                data.animation.mode = PlaybackMode::Stopped;
                data.animation.time = max_time;
            }
        }

        const int new_frame = math::clamp((int)(data.animation.time + 0.5f), 0, last_frame);
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

        if (data.time_filter.dynamic_window) {
            const auto half_window_ext = data.time_filter.window_extent * 0.5f;
            data.time_filter.range.min = math::clamp((f32)data.animation.time - half_window_ext, 0.0f, (f32)max_time);
            data.time_filter.range.max = math::clamp((f32)data.animation.time + half_window_ext, 0.0f, (f32)max_time);
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
            for (auto& rep : data.representations.buffer) {
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

#if 0
        PUSH_CPU_SECTION("Hydrogen bonds")
        if (data.hydrogen_bonds.enabled && data.hydrogen_bonds.dirty) {
            data.hydrogen_bonds.bonds = hydrogen_bond::compute_bonds(
                {mol.hydrogen_bond.donor.data, mol.hydrogen_bond.donor.count}, {mol.hydrogen_bond.acceptor.data, mol.hydrogen_bond.acceptor.count},
                mol.atom.position, data.hydrogen_bonds.distance_cutoff, math::deg_to_rad(data.hydrogen_bonds.angle_cutoff));
            data.hydrogen_bonds.dirty = false;
        }
        POP_CPU_SECTION()
#endif

        // Resize Framebuffer
        if ((data.fbo.width != data.ctx.framebuffer.width || data.fbo.height != data.ctx.framebuffer.height) &&
            (data.ctx.framebuffer.width != 0 && data.ctx.framebuffer.height != 0)) {
            init_framebuffer(&data.fbo, data.ctx.framebuffer.width, data.ctx.framebuffer.height);
            postprocessing::initialize(data.fbo.width, data.fbo.height);
        }

        update_view_param(&data);
        update_md_buffers(&data);

        //update_properties(&data);
        //update_density_volume_texture(&data);

        fill_gbuffer(&data);

#if EXPERIMENTAL_CONE_TRACED_AO == 1
        if (data.visuals.cone_traced_ao.enabled) {
            PUSH_GPU_SECTION("Cone-Trace AO")
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, data.fbo.deferred.fbo);
            glViewport(0, 0, data.fbo.width, data.fbo.height);
            glDrawBuffer(GL_COLOR_ATTACHMENT0);  // Modify to color buffer
            cone_trace::render_directional_occlusion(data.fbo.deferred.depth, data.fbo.deferred.normal, data.occupancy_volume.vol,
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
        if (data.ramachandran.show_window) draw_ramachandran_window(&data);
        //if (data.shape_space.show_window) draw_shape_space_window(&data);
        if (data.density_volume.show_window) draw_density_volume_window(&data);
        if (data.script.show_editor) draw_script_editor_window(&data);
        if (data.dataset.show_window) draw_dataset_window(&data);

        // @NOTE: ImGui::GetIO().WantCaptureMouse does not work with Menu
        if (!ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow)) {
            if (data.picking.idx != NO_PICKING_IDX) {
                draw_atom_info_window(data);
            }
        }

        PUSH_GPU_SECTION("Imgui render")
        platform::render_imgui(&data.ctx);
        POP_GPU_SECTION()

        // Swap buffers
        platform::swap_buffers(&data.ctx);

        task_system::run_main_tasks();

        // Reset frame allocator
        md_arena_allocator_reset(frame_allocator);
    }

    // shutdown subsystems
    LOG_NOTE("Shutting down immediate draw...");
    immediate::shutdown();
    LOG_NOTE("Shutting down ramachandran...");
    ramachandran::shutdown();
    LOG_NOTE("Shutting down post processing...");
    postprocessing::shutdown();
    LOG_NOTE("Shutting down volume...");
    volume::shutdown();
#if EXPERIMENTAL_CONE_TRACED_AO == 1
    LOG_NOTE("Shutting down cone tracing...");
    cone_trace::shutdown();
#endif
    LOG_NOTE("Shutting down task system...");
    task_system::shutdown();

    destroy_framebuffer(&data.fbo);
    platform::shutdown(&data.ctx);

    return 0;
}

static void interpolate_atomic_properties(ApplicationData* data) {
    ASSERT(data);
    const auto& mol = data->mold.mol;
    const auto& traj = data->mold.traj;
    const auto traj_ptr = &data->mold.traj;

    if (!mol.atom.count || !traj.num_frames) return;

    const i64 last_frame = math::max(0LL, traj.num_frames - 1);
    const f64 time = math::clamp(data->animation.time, 0.0, f64(last_frame));

    const f32 t = (float)math::fract(time);
    const i64 frame = (i64)time;
    const i64 nearest_frame = math::clamp((i64)(time + 0.5), 0LL, last_frame);

    const i64 frames[4] = {
        math::max(0LL, frame - 1),
        math::max(0LL, frame),
        math::min(frame + 1, last_frame),
        math::min(frame + 2, last_frame)
    };

    mat3 boxes[4] = {};

    const InterpolationMode mode = (frames[1] != frames[2]) ? data->animation.interpolation : InterpolationMode::Nearest;

    mat3 box = {};
    switch (mode) {
        case InterpolationMode::Nearest:
            load::traj::load_trajectory_frame_box(traj_ptr, (float*(*)[3])&box, nearest_frame);
            break;
        case InterpolationMode::Linear:
            load::traj::load_trajectory_frame_box(traj_ptr, (float*(*)[3])&boxes[1], frames[1]);
            load::traj::load_trajectory_frame_box(traj_ptr, (float*(*)[3])&boxes[2], frames[2]);
            box = math::lerp(boxes[1], boxes[2], t);
            break;
        case InterpolationMode::Cubic:
            load::traj::load_trajectory_frame_box(traj_ptr, (float*(*)[3])&boxes[0], frames[0]);
            load::traj::load_trajectory_frame_box(traj_ptr, (float*(*)[3])&boxes[1], frames[1]);
            load::traj::load_trajectory_frame_box(traj_ptr, (float*(*)[3])&boxes[2], frames[2]);
            load::traj::load_trajectory_frame_box(traj_ptr, (float*(*)[3])&boxes[3], frames[3]);
            box = math::cubic_spline(boxes[0], boxes[1], boxes[2], boxes[3], t);
            break;
        default:
            ASSERT(false);
    }

    const bool pbc = (box != mat3(0));
    data->simulation_box.box = box;

    switch (mode) {
        case InterpolationMode::Nearest: {
            soa_vec3 in_pos = {0};
            load::traj::load_trajectory_frame_coords(traj_ptr, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, nearest_frame);
            break;
        }
        case InterpolationMode::Linear: {
            int64_t stride = ROUND_UP(mol.atom.count, md_simd_width);    // The interploation uses SIMD vectorization without bounds, so we make sure there is no overlap between the data segments
            int64_t bytes = stride * sizeof(float) * 3 * 2;
            void* mem = md_alloc(frame_allocator, bytes);
            defer { md_free(frame_allocator, mem, bytes); };

            float* x[2] = {
                (float*)mem + stride * 0,
                (float*)mem + stride * 1,
            };
            float* y[2] = {
                (float*)mem + stride * 2,
                (float*)mem + stride * 3,
            };
            float* z[2] = {
                (float*)mem + stride * 4,
                (float*)mem + stride * 5,
            };

            load::traj::load_trajectory_frame_coords(traj_ptr, x[0], y[0], z[0], mol.atom.count, frames[1]);
            load::traj::load_trajectory_frame_coords(traj_ptr, x[1], y[1], z[1], mol.atom.count, frames[2]);

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

            break;
        }
        case InterpolationMode::Cubic: {
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

            load::traj::load_trajectory_frame_coords(traj_ptr, x[0], y[0], z[0], mol.atom.count, frames[0]);
            load::traj::load_trajectory_frame_coords(traj_ptr, x[1], y[1], z[1], mol.atom.count, frames[1]);
            load::traj::load_trajectory_frame_coords(traj_ptr, x[2], y[2], z[2], mol.atom.count, frames[2]);
            load::traj::load_trajectory_frame_coords(traj_ptr, x[3], y[3], z[3], mol.atom.count, frames[3]);

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

            break;
        }
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
            const size_t bytes = mol.backbone.count * sizeof(md_backbone_angles_t);
            memcpy(mol.backbone.angle, src_angle, bytes);
            break;
        }
        case InterpolationMode::Linear: {
            for (int64_t i = 0; i < mol.backbone.count; ++i) {
                float phi = math::lerp(src_angles[1][i].phi, src_angles[2][i].phi, t);
                float psi = math::lerp(src_angles[1][i].psi, src_angles[2][i].psi, t);
                mol.backbone.angle[i] = {phi, psi};
            }
            break;
        }
        case InterpolationMode::Cubic: {
            for (int64_t i = 0; i < mol.backbone.count; ++i) {
                float phi = math::cubic_spline(src_angles[0][i].phi, src_angles[1][i].phi, src_angles[2][i].phi, src_angles[3][i].phi, t);
                float psi = math::cubic_spline(src_angles[0][i].psi, src_angles[1][i].psi, src_angles[2][i].psi, src_angles[3][i].psi, t);
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
            const size_t bytes = mol.backbone.count * sizeof(md_secondary_structure_t);
            memcpy(mol.backbone.secondary_structure, ss, bytes);
            break;
        }
        case InterpolationMode::Linear: {
            for (int64_t i = 0; i < mol.backbone.count; ++i) {
                const vec4 ss_f[2] = {
                    math::convert_color((u32)src_ss[0][i]),
                    math::convert_color((u32)src_ss[1][i]),
                };
                const vec4 ss_res = math::lerp(ss_f[0], ss_f[1], t);
                mol.backbone.secondary_structure[i] = (md_secondary_structure_t)math::convert_color(ss_res);
            }
            break;
        }
        case InterpolationMode::Cubic: {
            for (int64_t i = 0; i < mol.backbone.count; ++i) {
                const vec4 ss_f[4] = {
                    math::convert_color((u32)src_ss[0][i]),
                    math::convert_color((u32)src_ss[1][i]),
                    math::convert_color((u32)src_ss[2][i]),
                    math::convert_color((u32)src_ss[3][i]),
                };
                const vec4 ss_res = math::cubic_spline(ss_f[0], ss_f[1], ss_f[2], ss_f[3], t);
                mol.backbone.secondary_structure[i] = (md_secondary_structure_t)math::convert_color(ss_res);
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
static f64 compute_avg_ms(f64 dt) {
    constexpr f64 interval = 0.5f;
    static f64 avg = 0.f;
    static int num_frames = 0;
    static f64 t = 0;
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

    param.clip_volume.near = data->view.camera.near_plane;
    param.clip_volume.far = data->view.camera.far_plane;
    param.fov_y = data->view.camera.fov_y;
    param.resolution = vec2(data->fbo.width, data->fbo.height);

    param.matrix.current.view = camera_world_to_view_matrix(data->view.camera);
    if (data->view.mode == CameraMode::Perspective) {
        param.matrix.current.proj = camera_perspective_projection_matrix(data->view.camera, data->fbo.width, data->fbo.height);
    } else {
        const float aspect_ratio = (float)data->fbo.width / (float)data->fbo.height;
        const float h = data->view.camera.focus_distance * math::tan(data->view.camera.fov_y * 0.5f);
        const float w = aspect_ratio * h;
        param.matrix.current.proj = camera_orthographic_projection_matrix(-w, w, -h, h, data->view.camera.near_plane, data->view.camera.far_plane);
    }
    param.matrix.current.proj_jittered = param.matrix.current.proj;

    if (data->visuals.temporal_reprojection.enabled && data->visuals.temporal_reprojection.jitter) {
        static u32 i = 0;
        i = (++i) % ARRAY_SIZE(data->view.jitter.sequence);
        param.jitter.next    = data->view.jitter.sequence[(i + 1) % ARRAY_SIZE(data->view.jitter.sequence)] - 0.5f;
        param.jitter.current = data->view.jitter.sequence[i] - 0.5f;
        if (data->view.mode == CameraMode::Perspective) {
            param.matrix.current.proj_jittered = camera_perspective_projection_matrix(data->view.camera, data->fbo.width, data->fbo.height,
                param.jitter.current.x, param.jitter.current.y);
        } else {
            const float aspect_ratio = (float)data->fbo.width / (float)data->fbo.height;
            const float h = data->view.camera.focus_distance * math::tan(data->view.camera.fov_y * 0.5f);
            const float w = aspect_ratio * h;
            const float scale_x = w / data->fbo.width * 2.0f;
            const float scale_y = h / data->fbo.height * 2.0f;
            const float j_x = param.jitter.current.x * scale_x;
            const float j_y = param.jitter.current.y * scale_y;
            param.matrix.current.proj_jittered = param.matrix.current.proj =
                camera_orthographic_projection_matrix(-w + j_x, w + j_x, -h + j_y, h + j_y, data->view.camera.near_plane, data->view.camera.far_plane);
        }
    }

    param.matrix.current.view_proj = param.matrix.current.proj * param.matrix.current.view;
    param.matrix.current.view_proj_jittered = param.matrix.current.proj_jittered * param.matrix.current.view;

    param.matrix.inverse.view = math::inverse(param.matrix.current.view);
    param.matrix.inverse.proj = math::inverse(param.matrix.current.proj);
    param.matrix.inverse.proj_jittered = math::inverse(param.matrix.current.proj_jittered);
    param.matrix.inverse.view_proj = math::inverse(param.matrix.current.view_proj);
    param.matrix.inverse.view_proj_jittered = math::inverse(param.matrix.current.view_proj_jittered);

    param.matrix.current.norm = math::transpose(param.matrix.inverse.view);
}

static void reset_view(ApplicationData* data, bool move_camera, bool smooth_transition) {
    ASSERT(data);
    if (!data->mold.mol.atom.count) return;
    const auto& mol = data->mold.mol;

    vec3 aabb_min, aabb_max;
    compute_aabb(&aabb_min, &aabb_max, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count);
    const vec3 ext = aabb_max - aabb_min;
    const vec3 cen = (aabb_min + aabb_max) * 0.5f;
    const vec3 pos = cen + ext * 3.f;

    if (move_camera) {
        if (!smooth_transition) data->view.camera.position = pos;
        data->view.animation.target_position = pos;
        data->view.camera.focus_distance = math::length(pos - cen);
        data->view.camera.orientation = math::conjugate(math::quat_cast(look_at(data->view.animation.target_position, cen, vec3(0, 1, 0))));
    }

    data->view.camera.near_plane = 1.f;
    data->view.camera.far_plane = math::length(ext) * 50.f;
    data->view.trackball_param.max_distance = math::length(ext) * 20.0f;
}

// #picking
static PickingData read_picking_data(const MainFramebuffer& framebuffer, i32 x, i32 y) {    
    static u32 frame = 0;
    u32 curr = (frame + 0) % 2;
    u32 prev = (frame + 1) % 2;

    PickingData data{};

    PUSH_GPU_SECTION("READ PICKING DATA")
    glBindFramebuffer(GL_READ_FRAMEBUFFER, framebuffer.deferred.fbo);
    glReadBuffer(GL_COLOR_ATTACHMENT_PICKING);

    // Queue async reads from current frame to pixel pack buffer
    glBindBuffer(GL_PIXEL_PACK_BUFFER, framebuffer.pbo_picking.color[curr]);
    glReadPixels(x, y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, framebuffer.pbo_picking.depth[curr]);
    glReadPixels(x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, 0);

    // Read values from previous frames pixel pack buffer
    glBindBuffer(GL_PIXEL_PACK_BUFFER, framebuffer.pbo_picking.color[prev]);
    const GLubyte* color = (const GLubyte*)glMapBuffer(GL_PIXEL_PACK_BUFFER, GL_READ_ONLY);
    if (color) {
        data.idx = color[0] + (color[1] << 8) + (color[2] << 16) + (color[3] << 24);
        glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
    }

    glBindBuffer(GL_PIXEL_PACK_BUFFER, framebuffer.pbo_picking.depth[prev]);
    const GLfloat* depth = (const GLfloat*)glMapBuffer(GL_PIXEL_PACK_BUFFER, GL_READ_ONLY);
    if (depth) {
        data.depth = depth[0];
        glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
    }

    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    POP_GPU_SECTION()

    ++frame;
    return data;
}

static void compute_aabb(vec3* aabb_min, vec3* aabb_max, const float* x, const float* y, const float* z, int64_t count) {
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

// #imgui
namespace ImGui {

static void CreateDockspace() {
    // Invisible dockspace
    ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(viewport->Pos);
    ImGui::SetNextWindowSize(viewport->Size);
    ImGui::SetNextWindowViewport(viewport->ID);
    ImGui::SetNextWindowBgAlpha(0.0f);

    const ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar |
ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoBringToFrontOnFocus |
ImGuiWindowFlags_NoNavFocus;

    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
    ImGui::Begin("DockspaceWindow", NULL, window_flags);
    ImGui::PopStyleVar(3);

    // ImGuiDockNodeFlags dockspace_flags = ImGuiDockNodeFlags_PassthruDockspace;
    const ImGuiID id = ImGui::GetID("Dockspace");
    const ImGuiDockNodeFlags flags = ImGuiDockNodeFlags_PassthruCentralNode;
    ImGui::DockSpace(id, ImVec2(0.0f, 0.0f), flags);

    ImGui::End();
}

static void BeginCanvas(const char* id) {
    // Invisible Canvas
    ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(viewport->Pos);
    ImGui::SetNextWindowSize(viewport->Size);
    ImGui::SetNextWindowViewport(viewport->ID);
    ImGui::SetNextWindowBgAlpha(0.0f);

    const ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar |
                                          ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                                          ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus | ImGuiWindowFlags_NoInputs;

    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
    ImGui::Begin(id, NULL, window_flags);
    ImGui::PopStyleVar(3);
}

static void EndCanvas() { ImGui::End(); }

static bool DeleteButton(const char* label, const ImVec2& size) {
    PushStyleColor(ImGuiCol_Button, DEL_BTN_COLOR);
    PushStyleColor(ImGuiCol_ButtonHovered, DEL_BTN_HOVER_COLOR);
    PushStyleColor(ImGuiCol_ButtonActive, DEL_BTN_ACTIVE_COLOR);
    defer { PopStyleColor(3); };
    return ImGui::Button(label, size);
}

}  // namespace ImGui

static void grow_mask_by_covalent_bond(Bitfield mask, Array<const Bond> bonds, i64 extent) {
    Bitfield prev_mask;
    bitfield::init(&prev_mask, mask);
    defer { bitfield::free(&prev_mask); };

    for (i64 i = 0; i < extent; i++) {
        for (const auto& bond : bonds) {
            const i64 idx[2] = {bond.idx[0], bond.idx[1]};
            if (bitfield::get_bit(prev_mask, idx[0]) && !bitfield::get_bit(mask, idx[1])) {
                bitfield::set_bit(mask, idx[1]);
            } else if (bitfield::get_bit(prev_mask, idx[1]) && !bitfield::get_bit(mask, idx[0])) {
                bitfield::set_bit(mask, idx[0]);
            }
        }
        memcpy(prev_mask.data(), mask.data(), prev_mask.size_in_bytes());
    }
}

static void grow_mask_by_radial_extent(Bitfield mask, const float* atom_x, const float* atom_y, const float* atom_z, i64 atom_count, float extent) {
    Bitfield prev_mask;
    bitfield::init(&prev_mask, mask);
    defer { bitfield::free(&prev_mask); };

    spatialhash::Frame frame = spatialhash::compute_frame(atom_x, atom_y, atom_z, atom_count, vec3(extent));
    for (i64 i = 0; i < atom_count; i++) {
        if (bitfield::get_bit(prev_mask, i)) {
            const vec3 pos = {atom_x[i], atom_y[i], atom_z[i]};
            spatialhash::for_each_within(frame, pos, extent, [mask](i32 idx, const vec3& pos) {
                UNUSED(pos);
                bitfield::set_bit(mask, idx);
            });
        }
    }
}

static void expand_mask(Bitfield mask, const md_range_t in_range[], i64 count) {
    for (i64 i = 0; i < count; i++) {
        Range<int> range = {in_range[i].beg, in_range[i].end};
        if (bitfield::any_bit_set_in_range(mask, range)) {
            bitfield::set_range(mask, range);
        }
    }
}

static bool filter_expression(const ApplicationData& data, CStringView expr, Bitfield mask) {
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

    str_t expression = {
        .ptr = expr.data(),
        .len = expr.size()
    };
    
    md_filter_context_t ctx {
        .mol = &data.mold.mol,
        .selection = {
            .count = stored_sel_count,
            .ptr = stored_sel,
        }
    };

    md_exp_bitfield_t bf = {0};
    md_bitfield_init(&bf, frame_allocator);

    bool result = md_filter_evaluate(expression, &bf, ctx);
    if (result) {
        md_bitfield_extract_u64(mask.data(), mask.size(), &bf);
    }
    else {
        memset(mask.data(), 0, mask.size_in_bytes());
    }

    md_bitfield_free(&bf);
    return result;
}

// ### DRAW WINDOWS ###
static void draw_main_menu(ApplicationData* data) {
    ASSERT(data);
    bool new_clicked = false;

    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("Load Data", "CTRL+L")) {
                auto res = platform::file_dialog(platform::FileDialogFlags_Open, {}, "pdb,gro,xtc");
                if (res.result == platform::FileDialogResult::Ok) {
                    load_dataset_from_file(data, res.path);
                    if (data->representations.buffer.empty()) {
                        create_representation(data); // Create default representation
                    }
                    data->animation = {};
                    reset_view(data, true);
                }
            }
            if (ImGui::MenuItem("Open Workspace", "CTRL+O")) {
                auto res = platform::file_dialog(platform::FileDialogFlags_Open, {}, FILE_EXTENSION);
                if (res.result == platform::FileDialogResult::Ok) {
                    load_workspace(data, res.path);
                }
            }
            if (ImGui::MenuItem("Save Workspace", "CTRL+S")) {
                if (!data->files.workspace) {
                    auto res = platform::file_dialog(platform::FileDialogFlags_Save, {}, FILE_EXTENSION);
                    if (res.result == platform::FileDialogResult::Ok) {
                        if (!get_file_extension(res.path)) {
                            snprintf(res.path.cstr() + strnlen(res.path.cstr(), res.path.capacity()), res.path.capacity(), ".%s",
                                     FILE_EXTENSION.cstr());
                        }
                        save_workspace(data, res.path);
                    }
                } else {
                    save_workspace(data, data->files.workspace);
                }
            }
            if (ImGui::MenuItem("Save As")) {
                auto res = platform::file_dialog(platform::FileDialogFlags_Save, {}, FILE_EXTENSION);
                if (res.result == platform::FileDialogResult::Ok) {
                    if (!get_file_extension(res.path)) {
                        snprintf(res.path.cstr() + strnlen(res.path.cstr(), res.path.capacity()), res.path.capacity(), ".%s", FILE_EXTENSION.cstr());
                    }
                    save_workspace(data, res.path);
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
                    float fov = math::rad_to_deg(data->view.camera.fov_y);
                    if (ImGui::SliderFloat("field of view", &fov, 12.5f, 80.0f)) {
                        data->view.camera.fov_y = math::deg_to_rad(fov);
                    }
                }
            }
            ImGui::EndGroup();

            ImGui::BeginGroup();
            ImGui::Text("Background");
            ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
            ImGui::ColorEdit3("Color", &data->visuals.background.color[0], ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_Float);
            ImGui::PopStyleVar();
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
                ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
                ImGui::ColorEdit4("Color", (float*)&data->simulation_box.color, ImGuiColorEditFlags_NoInputs);
                ImGui::PopStyleVar();
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
            ImGui::Checkbox("Ramachandran", &data->ramachandran.show_window);
            ImGui::Checkbox("Density Volumes", &data->density_volume.show_window);
            //ImGui::Checkbox("Shape Space", &data->shape_space.show_window);
            ImGui::Checkbox("Dataset", &data->dataset.show_window);

            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Selection")) {
            const auto atom_count = data->mold.mol.atom.count;
            Bitfield mask;
            bitfield::init(&mask, atom_count);
            defer { bitfield::free(&mask); };

            if (ImGui::MenuItem("Clear Selection")) {
                bitfield::clear_all(data->selection.current_selection_mask);
                data->mold.dirty_buffers |= MolBit_DirtyFlags;
            }

            if (ImGui::MenuItem("Invert Selection")) {
                bitfield::invert_all(data->selection.current_selection_mask);
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
                    query_ok = filter_expression(*data, CStringView(buf), mask);

                    if (query_ok) {
                        switch (data->selection.level_mode) {
                            case SelectionLevel::Atom:
                                // No need to expand the mask
                                break;
                            case SelectionLevel::Residue:
                                expand_mask(mask, data->mold.mol.residue.atom_range, data->mold.mol.residue.count);
                                break;
                            case SelectionLevel::Chain:
                                expand_mask(mask, data->mold.mol.chain.atom_range, data->mold.mol.chain.count);
                                break;
                            default:
                                ASSERT(false);
                        }
                    } else {
                        bitfield::clear_all(mask);
                    }
                    // data->mold.dirty_buffers |= MolBit_DirtyFlags;
                }

                const bool show_preview = (ImGui::GetFocusID() == ImGui::GetID("##query")) || (ImGui::GetHoveredID() == ImGui::GetID("##query")) ||
                                          (ImGui::GetHoveredID() == ImGui::GetID("Apply##query"));

                if (show_preview) {
                    // if (query_ok) {
                    memcpy(data->selection.current_highlight_mask.data(), mask.data(), mask.size_in_bytes());
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
                    memcpy(data->selection.current_selection_mask.data(), mask.data(), mask.size_in_bytes());
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
                    memcpy(mask.data(), data->selection.current_selection_mask.data(), data->selection.current_selection_mask.size_in_bytes());

                    switch (data->selection.grow_mode) {
                        case SelectionGrowth::CovalentBond:
                            grow_mask_by_covalent_bond(mask, {(Bond*)data->mold.mol.covalent_bond.bond, data->mold.mol.covalent_bond.count}, (i64)extent);
                            break;
                        case SelectionGrowth::Radial: {
                            const auto& mol = data->mold.mol;
                            grow_mask_by_radial_extent(mask, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, extent);
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
                            expand_mask(mask, data->mold.mol.residue.atom_range, data->mold.mol.residue.count);
                            break;
                        case SelectionLevel::Chain:
                            expand_mask(mask, data->mold.mol.chain.atom_range, data->mold.mol.chain.count);
                            break;
                        default:
                            ASSERT(false);
                    }

                    memcpy(data->selection.current_highlight_mask.data(), mask.data(), mask.size_in_bytes());
                    data->mold.dirty_buffers |= MolBit_DirtyFlags;

                    if (apply) {
                        memcpy(data->selection.current_selection_mask.data(), mask.data(), mask.size_in_bytes());
                    }
                    update_all_representations(data);
                }
            }

            ImGui::Spacing();
            ImGui::Separator();

            // STORED SELECTIONS
            {
                ImGui::Text("Stored Selections");
                const bool disable_new = bitfield::any_bit_set(data->selection.current_selection_mask) == false;
                if (disable_new) ImGui::PushDisabled();
                if (ImGui::Button("Create New")) {
                    char name_buf[64];
                    snprintf(name_buf, 64, "sel%i", (int)data->selection.stored_selections.size() + 1);
                    create_selection(data, name_buf, data->selection.current_selection_mask);
                }
                if (disable_new) ImGui::PopDisabled();
                for (int i = 0; i < data->selection.stored_selections.size(); i++) {
                    auto& sel = data->selection.stored_selections[i];
                    // const float32 item_width = math::clamp(ImGui::GetWindowContentRegionWidth() - 90.f, 100.f, 300.f);
                    StringBuffer<128> name;
                    snprintf(name.cstr(), name.capacity(), "%s###ID", sel.name.cstr());

                    ImGui::PushID(i);
                    if (ImGui::CollapsingHeader(name.cstr())) {
                        if (ImGui::Button("Activate")) {
                            memcpy(data->selection.current_selection_mask.data(), sel.atom_mask.data(), sel.atom_mask.size_in_bytes());
                            data->mold.dirty_buffers |= MolBit_DirtyFlags;
                        }
                        ImGui::SameLine();
                        if (ImGui::DeleteButton("Remove")) {
                            remove_selection(data, i);
                        }
                        ImGui::SameLine();
                        if (ImGui::Button("Save")) {
                            memcpy(sel.atom_mask.data(), data->selection.current_selection_mask.data(), sel.atom_mask.size_in_bytes());
                            update_all_representations(data);
                        }
                    }

                    const auto h_id = ImGui::GetHoveredID();
                    bool show_preview = (h_id == ImGui::GetID(name.cstr()) || h_id == ImGui::GetID("Activate") || h_id == ImGui::GetID("Remove") ||
                                         h_id == ImGui::GetID("Clone"));

                    ImGui::PopID();

                    if (show_preview) {
                        memcpy(data->selection.current_highlight_mask.data(), sel.atom_mask.data(), sel.atom_mask.size_in_bytes());
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
            const f64 ms = compute_avg_ms(data->ctx.timing.delta_s);
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
    const i32 num_frames = (int32_t)data->mold.traj.num_frames;
    ImGui::Text("Num Frames: %i", num_frames);
    // ImGui::Checkbox("Apply post-interpolation pbc", &data->animation.apply_pbc);
    f32 t = (float)data->animation.time;
    if (ImGui::SliderFloat("Time", &t, 0, (float)(math::max(0, num_frames - 1)))) {
        data->animation.time = t;
    }
    //ImGui::SliderFloat("Speed", &data->animation.fps, 0.1f, 1000.f, "%.3f", 4.f);
    ImGui::SliderFloat("Speed", &data->animation.fps, 0.1f, 1000.f, "%.3f", ImGuiSliderFlags_Logarithmic);
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("Animation Speed in Frames Per Second");
    }
    if (ImGui::Combo("Interpolation", (int*)(&data->animation.interpolation), "Nearest\0Linear\0Cubic\0\0")) {
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
    for (int i = 0; i < data->representations.buffer.size(); i++) {
        bool update_args = false;
        bool update_color = false;
        auto& rep = data->representations.buffer[i];
        const f32 item_width = math::clamp(ImGui::GetWindowContentRegionWidth() - 90.f, 100.f, 300.f);
        StringBuffer<128> name;
        snprintf(name.cstr(), name.capacity(), "%s###ID", rep.name.buffer);

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
            md_gl_representation_type type = MD_GL_REP_DEFAULT;
            md_gl_representation_args args = {};
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

static void draw_atom_info_window(const ApplicationData& data) {

    const auto& mol = data.mold.mol;
    int atom_idx = data.picking.idx;

    // @TODO: Assert things and make this failproof
    if (atom_idx < 0 || atom_idx >= mol.atom.count) return;

    int res_idx = mol.atom.residue_idx[atom_idx];
    const char* res_name = mol.residue.name[res_idx];
    const int res_id = mol.residue.id[res_idx];
    int local_idx = atom_idx - mol.residue.atom_range[res_idx].beg;
    const vec3 pos = { mol.atom.x[atom_idx], mol.atom.y[atom_idx], mol.atom.z[atom_idx] };
    const char* label = mol.atom.name[atom_idx];
    str_t elem = md_util_element_name(mol.atom.element[atom_idx]);
    str_t symbol = md_util_element_symbol(mol.atom.element[atom_idx]);

    int chain_idx = -1;
    const char* chain_id = "\0";
    if (mol.atom.chain_idx) {
        chain_idx = mol.atom.chain_idx[atom_idx];
        if (chain_idx != INVALID_CHAIN_IDX && mol.chain.count > 0) {
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
        const auto angles = math::rad_to_deg((vec2)mol.backbone.angles[res_idx]);
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
            const auto file = get_file(data->files.molecule);
            if (file) ImGui::Text("\"%.*s\"", (i32)file.size(), file.data());
            ImGui::Text("# atoms: %lli", mol.atom.count);
            ImGui::Text("# residues: %lli", mol.residue.count);
            ImGui::Text("# chains: %lli", mol.chain.count);
        }
        if (traj.num_frames) {
            ImGui::NewLine();
            ImGui::Text("TRAJ");
            const auto file = get_file(data->files.trajectory);
            ImGui::Text("\"%.*s\"", (i32)file.size(), file.data());
            ImGui::Text("# frames: %i", traj.num_frames);
        }
        const i64 selection_count = bitfield::number_of_bits_set(data->selection.current_selection_mask);
        const i64 highlight_count = bitfield::number_of_bits_set(data->selection.current_highlight_mask);
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
    constexpr f32 WIDTH = 300.f;
    constexpr f32 MARGIN = 10.f;
    constexpr f32 PROGRESSBAR_WIDTH_FRACT = 0.3f;

    //const f32 stats_fract = stats::fraction_done();

    const u32 num_tasks = task_system::get_num_tasks();
    task_system::ID* tasks = task_system::get_tasks();

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
        for (u32 i = 0; i < num_tasks; i++) {
            const auto id = tasks[i];
            const f32 fract = task_system::get_task_fraction_complete(id);
            snprintf(buf, 32, "%.1f%%", fract * 100.f);
            ImGui::ProgressBar(fract, ImVec2(ImGui::GetWindowContentRegionWidth() * PROGRESSBAR_WIDTH_FRACT, 0), buf);
            ImGui::SameLine();
            ImGui::Text("%s", task_system::get_task_label(id));
            ImGui::SameLine();
            if (ImGui::Button("X")) {
                task_system::interrupt_task(id);
            }
        }

        ImGui::End();
        ImGui::PopStyleColor();
    }
}

static void draw_timeline_window(ApplicationData* data) {
    ASSERT(data);
    ImGui::SetNextWindowSize(ImVec2(600, 300), ImGuiCond_FirstUseEver);

    if (ImGui::Begin("Temporal", &data->statistics.show_timeline_window, ImGuiWindowFlags_NoFocusOnAppearing)) {

        // convenience struct to manage DND items; do this however you like
        struct DndItem {
            char                label[16];
            int                 num_values;
            float*              values;
            ImVec4              color;
            int                 plot;
        };

        static DndItem  dnd[32];
        int             num_dnd = 0;

        // This is for the time stamps along the x-axis.
        // If we don't have any specific time points for the frames, we just use the indices as time points.

        int num_time_values = (int)data->mold.traj.num_frames;
        float* time_values = (float*)md_alloc(frame_allocator, num_time_values * sizeof(float));
        defer {
            md_free(frame_allocator, time_values, num_time_values * sizeof(float));
        };
        for (int64_t i = 0; i < num_time_values; ++i) {
            //time_values[i] = (float)i / (float)(num_time_values - 1);
            time_values[i] = (float)i;
        }

        for (int64_t i = 0; i < data->mold.script.eval.num_properties; i++) {
            auto& prop = data->mold.script.eval.properties[i];
            if (prop.type != MD_SCRIPT_PROPERTY_TYPE_TEMPORAL) continue;

            ASSERT(num_dnd < ARRAY_SIZE(dnd));
            int idx = num_dnd++;
            dnd[idx] = {
                .label = {0},
                .num_values = (int)prop.data.num_values,
                .values = prop.data.values,
                .color = vec_cast(qualitative_color_scale(idx)),
                .plot = dnd[idx].plot
            };
            const size_t cpy_size = ARRAY_SIZE(dnd[num_dnd-1].label) < prop.ident.len ? ARRAY_SIZE(dnd[num_dnd-1].label) : prop.ident.len;
            strncpy(dnd[num_dnd-1].label, prop.ident.ptr, cpy_size);
        }

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
        
        if (ImGui::BeginDragDropTarget()) {
            if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload("MY_DND")) {
                int i = *(int*)payload->Data;
                dnd[i].plot = 0;
            }
            ImGui::EndDragDropTarget();
        }

        ImGui::SameLine();
        ImGui::BeginChild("DND_RIGHT",ImVec2(-1,-1));
        
        //ImPlot::SetNextPlotLimitsX(0.0, data->mold.traj.num_frames, ImGuiCond_Always);
        ImPlotAxisFlags axis_flags = 0;
        ImPlotAxisFlags axis_flags_x = axis_flags;
        ImPlotAxisFlags axis_flags_y = axis_flags | ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_RangeFit;
        ImPlotFlags flags = ImPlotFlags_AntiAliased;
        //ImPlot::SetNextPlotLimitsX(-10, 11.0, ImGuiCond_Always);

        //ImPlot::SetNextPlotLimitsX(0.0, (double)num_time_values, ImGuiCond_FirstUseEver);

        static ImPlotLimits plot_limits;

        static bool need_refit = false;
        if (need_refit) {
            ImPlot::FitNextPlotAxes();
            need_refit = false;
        }

        bool set_limits = false;
        if (plot_limits.X.Min < 0) {
            plot_limits.X.Min = 0;
            set_limits = true;
        }

        if (plot_limits.X.Max > time_values[num_time_values-1]) {
            plot_limits.X.Max = time_values[num_time_values-1];
            set_limits = true;
        }

        if (set_limits) ImPlot::SetNextPlotLimits(plot_limits.X.Min, plot_limits.X.Max, plot_limits.Y.Min, plot_limits.Y.Max, ImGuiCond_Always);

        if (ImPlot::BeginPlot("##DND1", NULL, NULL, ImVec2(-1,-1), flags, axis_flags_x, axis_flags_y)) {
            for (int k = 0; k < num_dnd; ++k) {
                if (dnd[k].plot == 1 && dnd[k].num_values > 0) {
                    ASSERT(dnd[k].num_values == num_time_values);
                    //ImPlot::SetPlotYAxis(0);
                    
                    ImPlot::DragLineX("Current Time", &data->animation.time, true, ImVec4(1,1,0,1));

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

            plot_limits = ImPlot::GetPlotLimits();

            ImPlot::EndPlot();
        }
        ImGui::EndChild();
    }
    ImGui::End();
}

static void draw_distribution_window(ApplicationData* data) {
    ImGui::SetNextWindowSize(ImVec2(200, 300), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Distributions", &data->statistics.show_distribution_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
        static const int MIN_NUM_BINS = 8;
        static const int MAX_NUM_BINS = 1024;

        static int bins = 32;
        static int bin_technique = 0;
        /*
        ImGui::Combo("Bins", &bin_technique, "N Bins\0Sqrt\0Sturges\0Rice\0Scott\0\0");

        if (bin_technique == 0) {
            // Explicit number of bins
            if (ImGui::DragInt("Num Bins", &bins, 1, MIN_NUM_BINS, MAX_NUM_BINS)) {
                // @TODO: Round this to some even multiple of 2???
            }
        }
        */

        ImPlotAxisFlags axis_flags = 0;
        ImPlotAxisFlags axis_flags_x = axis_flags | ImPlotAxisFlags_AutoFit;
        ImPlotAxisFlags axis_flags_y = axis_flags | ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_RangeFit;
        ImPlotFlags flags = ImPlotFlags_AntiAliased;

        for (int64_t i = 0; i < data->mold.script.eval.num_properties; ++i) {
            md_script_property_t& prop = data->mold.script.eval.properties[i];
            if (prop.type != MD_SCRIPT_PROPERTY_TYPE_TEMPORAL && prop.type != MD_SCRIPT_PROPERTY_TYPE_DISTRIBUTION) continue;

            char label[16] = {0};
            strncpy(label, prop.ident.ptr, MIN(ARRAY_SIZE(label) - 1, prop.ident.len));
            
            ImGui::PushID(i);
            if (ImPlot::BeginPlot("##Distributions", 0, 0, ImVec2(-1,200), flags, axis_flags_x, axis_flags_y)) {
                ImPlot::SetNextFillStyle(vec_cast(qualitative_color_scale(i)), 0.5f);
                ImPlotRange range;
                ImPlot::PlotHistogram(label, prop.data.values, prop.data.num_values, bins, false, true, range);
                ImPlot::EndPlot();
            }
            ImGui::PopID();
        }
    }
    ImGui::End();
}

static void draw_ramachandran_window(ApplicationData* data) {
    // const int32 num_frames = data->mold.traj ? data->mold.traj.num_frames : 0;
    // const int32 frame = (int32)data->time;
    const Range<i32> frame_range = data->time_filter.range;
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

    constexpr f32 ONE_OVER_TWO_PI = 1.f / (2.f * math::PI);

    i64 mouse_hover_idx = -1;

    if (data->ramachandran.current.enabled) {
        const u32 border_color = math::convert_color(data->ramachandran.current.border_color);
        const u32 base_color = math::convert_color(data->ramachandran.current.base.fill_color);
        const u32 selected_color = math::convert_color(data->ramachandran.current.selection.selection_color);
        const u32 highlight_color = math::convert_color(data->ramachandran.current.selection.highlight_color);
        const f32 base_radius = data->ramachandran.current.base.radius;
        const f32 selected_radius = data->ramachandran.current.selection.radius;

        for (i64 ci = 0; ci < mol.chain.count; ++ci) {
            const auto range = mol.chain.backbone_range[ci];
            for (i64 i = range.beg; i < range.end; i++) {
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

        for (i64 i = 0; i < mol.backbone.count; ++i) {
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
                const float new_zoom = math::clamp(zoom_factor + zoom_factor * ZOOM_SCL * mouse_wheel_delta, 1.f, 10.f);
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

static void draw_density_volume_window(ApplicationData* data) {
    ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Density Volume", &data->density_volume.show_window, ImGuiWindowFlags_MenuBar)) {
        const ImVec2 button_size = {160, 20};

        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("File")) {
                if (ImGui::MenuItem("Export volume")) {
                    auto res = platform::file_dialog(platform::FileDialogFlags_Open, {}, "raw");
                    if (res.result == platform::FileDialogResult::Ok) {
                        //volume::write_to_file(data->density_volume.volume, res.path);
                        LOG_NOTE("Wrote density volume");
                    }
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("DVR")) {
                ImGui::Checkbox("Enabled", &data->density_volume.dvr.enabled);
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
                ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0, 0, 0, 1.0f));
                ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
                ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(1.0f, 1.0f));
                if (ImGui::ImageButton((void*)(intptr_t)data->density_volume.dvr.tf.id, button_size)) {
                    auto res = platform::file_dialog(platform::FileDialogFlags_Open, {}, "png,jpg");
                    if (res.result == platform::FileDialogResult::Ok) {
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
                    ImGui::ColorEdit4("##Color", &data->density_volume.iso.isosurfaces.colors[i][0],
                        ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_Float);
                    ImGui::PopID();
                }
                if ((data->density_volume.iso.isosurfaces.count < data->density_volume.iso.isosurfaces.MaxCount) &&
                    ImGui::Button("Add", button_size)) {
                    insert(data->density_volume.iso.isosurfaces, 0.1f, {0.2f, 0.1f, 0.9f, 1.0f});
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

            if (ImGui::BeginMenu("Scale")) {
                static int res_scale = 1;
                if (ImGui::SliderInt("Scale", &res_scale, 1, 8)) {

                }
                ImGui::EndMenu();
            }

            ImGui::EndMenuBar();
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
        draw_list->AddImage((ImTextureID)data->density_volume.render_texture.id, canvas_p0, canvas_p1, {0,1}, {1,0});
        //draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));
        draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));

        // This will catch our interactions
        ImGui::InvisibleButton("canvas", canvas_sz, ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);
        const bool is_hovered = ImGui::IsItemHovered(); // Hovered
        const bool is_active = ImGui::IsItemActive();   // Held
        const ImVec2 origin(canvas_p0.x, canvas_p0.y);  // Lock scrolled origin
        const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

        static Camera cam;

        static bool init = true;
        if (init) {

            const vec3 ext = {1,1,1};
            const vec3 cen = {0.5,0.5,0.5};
            const vec3 pos = cen + ext * vec3(1,0.5,0.7) * 2.f;
            const vec3 up = {0,1,0};

            cam.position = pos;
            cam.focus_distance = math::length(pos - cen);
            cam.orientation = math::conjugate(math::quat_cast(look_at(cam.position, cen, up)));
            init = false;
        }

        static TrackballControllerParam param = {
            .min_distance = 1.0,
            .max_distance = 100.0,
        };
        vec2 delta = {data->ctx.input.mouse.win_delta.x, data->ctx.input.mouse.win_delta.y};
        vec2 curr = {mouse_pos_in_canvas.x, mouse_pos_in_canvas.y};
        vec2 prev = curr - delta;
        TrackballControllerInput input = {
            .rotate_button = is_active && data->ctx.input.mouse.down[0],
            .pan_button = is_active && data->ctx.input.mouse.down[1],
            .dolly_button = is_active && data->ctx.input.mouse.down[2],
            .dolly_delta = is_hovered ? data->ctx.input.mouse.scroll_delta : 0.0f,
            .mouse_coord_prev = prev,
            .mouse_coord_curr = curr,
            .screen_size = {canvas_sz.x, canvas_sz.y},
            .fov_y = cam.fov_y,
        };
        camera_controller_trackball(&cam.position, &cam.orientation, &cam.focus_distance, input, param);

        //cam.position = {5,5,5};
        //cam.orientation = math::conjugate(math::quat_cast(look_at(cam.position, {0,0,0}, {0,1,0})));

        if (data->density_volume.render_texture.id == 0||
            data->density_volume.render_texture.dim.x != canvas_sz.x ||
            data->density_volume.render_texture.dim.y != canvas_sz.y) {
            data->density_volume.render_texture.dim = {canvas_sz.x, canvas_sz.y};
            volume::init_texture_2D(&data->density_volume.render_texture.id, data->density_volume.render_texture.dim.x, data->density_volume.render_texture.dim.y, GL_RGBA8);
        }

        if (!data->density_volume.volume_texture.id) {
            volume::init_texture_3D(&data->density_volume.volume_texture.id, 128, 128, 128, GL_R32F);
        }

        for (int64_t i = 0; i < data->mold.script.eval.num_properties; ++i) {
            md_script_property_t* prop = &data->mold.script.eval.properties[i];
            if (prop->type == MD_SCRIPT_PROPERTY_TYPE_VOLUME && prop->data.values) {
                volume::set_texture_3D_data(data->density_volume.volume_texture.id, prop->data.values, GL_R32F);
            }
        }

        vec3 min_aabb = {0,0,0};
        vec3 max_aabb = {1,1,1};

        mat4 model_mat = volume::compute_model_to_world_matrix(min_aabb, max_aabb);
        mat4 view_mat = camera_world_to_view_matrix(cam);
        mat4 proj_mat = camera_perspective_projection_matrix(cam, canvas_sz.x, canvas_sz.y);

        volume::RenderDesc desc = {
            .render_target = {
                .texture = data->density_volume.render_texture.id,
                .width = data->density_volume.render_texture.dim.x,
                .height = data->density_volume.render_texture.dim.y,
            },
            .texture = {
                .volume = data->density_volume.volume_texture.id,
                .transfer_function = data->density_volume.dvr.tf.id,
                .depth = data->fbo.deferred.depth,
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
            .isosurface = data->density_volume.iso.isosurfaces,
            .isosurface_enabled = data->density_volume.iso.enabled,
            .direct_volume_rendering_enabled = data->density_volume.dvr.enabled,
            .bounding_box_enabled = true,
            .voxel_spacing = data->density_volume.voxel_spacing
        };
        volume::render_volume(desc);
    }

    ImGui::End();
}

static void draw_dataset_window(ApplicationData* data) {
    ImGui::Begin("Dataset", &data->dataset.show_window);
    ImGui::SetWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);

    CStringView mol_file = data->files.molecule;
    ImGui::Text("Molecular data: %.*s", (int)mol_file.length(), mol_file.cstr());
    ImGui::Text("Num atoms:    %9lli", data->mold.mol.atom.count);
    ImGui::Text("Num residues: %9lli", data->mold.mol.residue.count);
    ImGui::Text("Num chains:   %9lli", data->mold.mol.chain.count);

    CStringView traj_file = data->files.trajectory;
    if (traj_file) {
        ImGui::Separator();
        ImGui::Text("Trajectory data: %.*s", (int)traj_file.length(), traj_file.cstr());
        ImGui::Text("Num frames:    %9lli", data->mold.traj.num_frames);
        ImGui::Text("Num atoms:     %9lli", data->mold.traj.num_atoms);
    }

    ImGui::End();
}

static void compile_script_in_editor(ApplicationData* data) {
    TextEditor& editor = data->script.editor;

    std::string src = editor.GetText();
    str_t str = {.ptr = src.data(), .len = (int64_t)src.length()};
    md_script_ir_compile_args_t args = {
        .src = str,
        .mol = &data->mold.mol,
        .alloc = default_allocator,
    };
    md_script_ir_compile(&data->mold.script.ir, args);

    TextEditor::ErrorMarkers markers;
    const int64_t num_errors = data->mold.script.ir.num_errors;
    if (num_errors) {
        const md_script_error_t* errors = data->mold.script.ir.errors;
        for (int64_t i = 0; i < num_errors; ++i) {
            std::string err_str(errors[i].error.ptr, errors[i].error.len);
            std::pair<int, std::string> pair = {errors[i].line, err_str};
            markers.insert(pair);
        }
    }
    editor.SetErrorMarkers(markers);
}

static void draw_script_editor_window(ApplicationData* data) {
    ASSERT(data);

    TextEditor& editor = data->script.editor;

    ImGui::Begin("Script Editor", &data->script.show_editor, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_MenuBar);
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

        if (editor.IsTextChanged()) {
            compile_script_in_editor(data);
        }

        if (ImGui::MenuItem("Evaluate")) {
            md_script_eval_args_t args = {
                .ir = &data->mold.script.ir,
                .mol = &data->mold.mol,
                .traj = &data->mold.traj,
                .frame_mask = 0,
                .alloc = default_allocator
            };
            if (md_script_eval(&data->mold.script.eval, args)) {
                LOG_NOTE("Evaluation successful");
            } else {
                LOG_NOTE("Evaluation failed");
            }
        }
        ImGui::EndMenuBar();
    }

    md_script_tokens_t tokens = {0};
    md_script_tokens_init(&tokens, &data->mold.script.ir, frame_allocator);

    editor.ClearMarkers();
    if (tokens.num_tokens > 0) {
        for (int64_t i = 0; i < tokens.num_tokens; ++i) {
            const md_script_token_t& tok = tokens.tokens[i];
            TextEditor::Marker marker = {0};
            marker.begCol = tok.col_beg;
            marker.endCol = tok.col_end;
            marker.bgColor = ImVec4(1,1,1,0.5);
            marker.depth = tok.depth;
            marker.line = tok.line;
            marker.onlyShowBgOnMouseOver = true;
            marker.text = std::string(tok.text.ptr, tok.text.len);
            marker.payload = (void*)&tok;
            editor.AddMarker(marker);
        }
    }

    editor.Render("TextEditor");

    const TextEditor::Marker* hovered_marker = editor.GetHoveredMarker();
    if (hovered_marker) {
        md_script_visualization_t vis = {0};
        md_script_visualization_init(&vis, {
            .token = (md_script_token_t*)hovered_marker->payload,
            .ir = &data->mold.script.ir,
            .mol = &data->mold.mol,
            .alloc = frame_allocator,
        });

        immediate::set_model_view_matrix(data->view.param.matrix.current.view);
        immediate::set_proj_matrix(data->view.param.matrix.current.proj_jittered);

        const vec3* vertices = (const vec3*)vis.vertex.pos;
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
            md_bitfield_extract_u64(data->selection.current_highlight_mask.data(), data->selection.current_highlight_mask.size(), &vis.atom_mask);
            data->mold.dirty_buffers |= MolBit_DirtyFlags;
        }
        md_script_visualization_free(&vis);
    }

    ImGui::End();

}

// static void draw_density_volume_clip_plane_widgets(ApplicationData* data) { ASSERT(data); }

// #framebuffer
static void init_framebuffer(MainFramebuffer* fbo, int width, int height) {
    ASSERT(fbo);

    bool attach_textures_deferred = false;
    if (!fbo->deferred.fbo) {
        glGenFramebuffers(1, &fbo->deferred.fbo);
        attach_textures_deferred = true;
    }

    if (!fbo->deferred.depth) glGenTextures(1, &fbo->deferred.depth);
    if (!fbo->deferred.color) glGenTextures(1, &fbo->deferred.color);
    if (!fbo->deferred.normal) glGenTextures(1, &fbo->deferred.normal);
    if (!fbo->deferred.velocity) glGenTextures(1, &fbo->deferred.velocity);
    if (!fbo->deferred.post_tonemap) glGenTextures(1, &fbo->deferred.post_tonemap);
    if (!fbo->deferred.picking) glGenTextures(1, &fbo->deferred.picking);
    if (!fbo->pbo_picking.color[0]) glGenBuffers(2, fbo->pbo_picking.color);
    if (!fbo->pbo_picking.depth[0]) glGenBuffers(2, fbo->pbo_picking.depth);

    glBindTexture(GL_TEXTURE_2D, fbo->deferred.depth);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH24_STENCIL8, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, fbo->deferred.color);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, fbo->deferred.normal);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16, width, height, 0, GL_RG, GL_UNSIGNED_SHORT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, fbo->deferred.velocity);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, width, height, 0, GL_RG, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, fbo->deferred.post_tonemap);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, fbo->deferred.picking);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, fbo->pbo_picking.color[0]);
    glBufferData(GL_PIXEL_PACK_BUFFER, 4, NULL, GL_DYNAMIC_READ);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, fbo->pbo_picking.color[1]);
    glBufferData(GL_PIXEL_PACK_BUFFER, 4, NULL, GL_DYNAMIC_READ);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, fbo->pbo_picking.depth[0]);
    glBufferData(GL_PIXEL_PACK_BUFFER, 4, NULL, GL_DYNAMIC_READ);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, fbo->pbo_picking.depth[1]);
    glBufferData(GL_PIXEL_PACK_BUFFER, 4, NULL, GL_DYNAMIC_READ);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    glBindTexture(GL_TEXTURE_2D, 0);

    fbo->width = width;
    fbo->height = height;

    const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
                                   GL_COLOR_ATTACHMENT_POST_TONEMAP, GL_COLOR_ATTACHMENT_PICKING};

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo->deferred.fbo);
    if (attach_textures_deferred) {
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, fbo->deferred.depth, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_STENCIL_ATTACHMENT, GL_TEXTURE_2D, fbo->deferred.depth, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_COLOR, GL_TEXTURE_2D, fbo->deferred.color, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_NORMAL, GL_TEXTURE_2D, fbo->deferred.normal, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_VELOCITY, GL_TEXTURE_2D, fbo->deferred.velocity, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_POST_TONEMAP, GL_TEXTURE_2D, fbo->deferred.post_tonemap, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT_PICKING, GL_TEXTURE_2D, fbo->deferred.picking, 0);
    }
    ASSERT(glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
    glDrawBuffers(ARRAY_SIZE(draw_buffers), draw_buffers);
    glClearColor(0, 0, 0, 0);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
}

static void destroy_framebuffer(MainFramebuffer* fbo) {
    ASSERT(fbo);
    if (fbo->deferred.fbo) glDeleteFramebuffers(1, &fbo->deferred.fbo);
    if (fbo->deferred.depth) glDeleteTextures(1, &fbo->deferred.depth);
    if (fbo->deferred.color) glDeleteTextures(1, &fbo->deferred.color);
    if (fbo->deferred.normal) glDeleteTextures(1, &fbo->deferred.normal);
    if (fbo->deferred.post_tonemap) glDeleteTextures(1, &fbo->deferred.post_tonemap);
    if (fbo->deferred.picking) glDeleteTextures(1, &fbo->deferred.picking);

    if (fbo->pbo_picking.color[0]) glDeleteBuffers(2, fbo->pbo_picking.color);
    if (fbo->pbo_picking.depth[0]) glDeleteBuffers(2, fbo->pbo_picking.depth);
}

static void update_md_buffers(ApplicationData* data) {
    ASSERT(data);
    const auto& mol = data->mold.mol;

    if (data->mold.dirty_buffers & MolBit_DirtyPosition) {
        md_gl_molecule_update_atom_previous_position(&data->mold.gl_mol);
        md_gl_molecule_set_atom_position(&data->mold.gl_mol, 0, (u32)mol.atom.count, mol.atom.x, mol.atom.y, mol.atom.z, 0);
    }

    if (data->mold.dirty_buffers & MolBit_DirtyFlags) {
        if (data->mold.mol.atom.flags) {
            for (i64 i = 0; i < mol.atom.count; i++) {
                u8 flags = 0;
                flags |= bitfield::get_bit(data->selection.current_highlight_mask, i)     ? AtomBit_Highlighted : 0;
                flags |= bitfield::get_bit(data->selection.current_selection_mask, i)     ? AtomBit_Selected : 0;
                flags |= bitfield::get_bit(data->representations.atom_visibility_mask, i) ? AtomBit_Visible : 0;
                data->mold.mol.atom.flags[i] = flags;
            }
            md_gl_molecule_set_atom_flags(&data->mold.gl_mol, 0, (u32)mol.atom.count, mol.atom.flags, 0);
        }
    }

    if (data->mold.dirty_buffers & MolBit_DirtySecondaryStructure) {
        if (mol.backbone.secondary_structure) {
            md_gl_molecule_set_backbone_secondary_structure(&data->mold.gl_mol, 0, (u32)mol.backbone.count, mol.backbone.secondary_structure, 0);
        }
    }

    data->mold.dirty_buffers = 0;
}

static void interrupt_async_tasks(ApplicationData* data) {
    if (data->tasks.load_trajectory.id != 0) {
        task_system::interrupt_task(data->tasks.load_trajectory);
    }
    task_system::wait_for_task(data->tasks.load_trajectory);
    data->tasks.load_trajectory.id = 0;
}

// #moleculedata
static void free_trajectory_data(ApplicationData* data) {
    ASSERT(data);
    if (data->mold.traj.num_frames) {
        if (data->tasks.load_trajectory.id != 0) task_system::interrupt_and_wait(data->tasks.load_trajectory);
        load::traj::close(&data->mold.traj);
    }
}

static void free_molecule_data(ApplicationData* data) {
    ASSERT(data);
    if (data->mold.mol.atom.count) {
        data->files.molecule = "";
        load::mol::free(&data->mold.mol);
    }
    if (data->mold.traj.num_atoms) {
        data->files.trajectory = "";
        load::traj::close(&data->mold.traj);
    }
    bitfield::clear_all(data->selection.current_selection_mask);
    bitfield::clear_all(data->selection.current_highlight_mask);
}

static void init_molecule_data(ApplicationData* data) {
    if (data->mold.mol.atom.count) {
        const auto& mol = data->mold.mol;
        bitfield::init(&data->selection.current_selection_mask, mol.atom.count);
        bitfield::init(&data->selection.current_highlight_mask, mol.atom.count);
        data->picking.idx = NO_PICKING_IDX;
        data->selection.hovered = -1;
        data->selection.right_clicked = -1;

        md_gl_molecule_desc desc = {
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
        compile_script_in_editor(data);
    }
}

static void init_trajectory_data(ApplicationData* data) {
    if (data->mold.traj.num_frames) {
        data->time_filter.range = {0, (float)data->mold.traj.num_frames};
        data->animation.time = math::clamp(data->animation.time, (f64)data->time_filter.range.min, (f64)data->time_filter.range.max);
        data->animation.frame = math::clamp((i32)data->animation.time, 0, (int32_t)data->mold.traj.num_frames - 1);
        i32 frame = data->animation.frame;

        load::traj::load_trajectory_frame_box(&data->mold.traj, (float*(*)[3])&data->simulation_box.box, frame);
        load::traj::load_trajectory_frame_coords(&data->mold.traj, data->mold.mol.atom.x, data->mold.mol.atom.y, data->mold.mol.atom.z, data->mold.mol.atom.count, frame);
        data->mold.dirty_buffers |= MolBit_DirtyPosition;

        update_md_buffers(data);
        md_gl_molecule_update_atom_previous_position(&data->mold.gl_mol); // Do this explicitly to update the previous position to avoid motion blur trails

        if (data->mold.mol.backbone.count > 0) {
            data->trajectory_data.secondary_structure = {
                .stride = data->mold.mol.backbone.count,
                .count = data->mold.mol.backbone.count * data->mold.traj.num_frames,
                .data = (md_secondary_structure_t*)REALLOC(data->trajectory_data.secondary_structure.data, data->mold.mol.backbone.count * data->mold.traj.num_frames * sizeof(md_secondary_structure_t))
            };

            data->trajectory_data.backbone_angles = {
                .stride = data->mold.mol.backbone.count,
                .count = data->mold.mol.backbone.count * data->mold.traj.num_frames,
                .data = (md_backbone_angles_t*)REALLOC(data->trajectory_data.backbone_angles.data, data->mold.mol.backbone.count * data->mold.traj.num_frames * sizeof(md_backbone_angles_t))
            };

            // Compute angles
            task_system::enqueue_pool("Performing trajectory data computations",
                //(uint32_t)data->mold.traj.num_frames,
                1,
                [&traj_data = data->trajectory_data, data](task_system::TaskSetRange range) {
                    int64_t bytes = data->mold.mol.atom.count * 3 * sizeof(float);
                    void* mem = md_alloc(default_allocator, bytes);
                    defer { md_free(default_allocator, mem, bytes); };

                    float* x = (float*)mem + 0 * data->mold.mol.atom.count;
                    float* y = (float*)mem + 1 * data->mold.mol.atom.count;
                    float* z = (float*)mem + 2 * data->mold.mol.atom.count;

                    for (int64_t f_idx = 0; f_idx < data->mold.traj.num_frames; ++f_idx) {
                    //for (u32 f_idx = range.beg; f_idx < range.end; f_idx++) {
                        load::traj::load_trajectory_frame_coords(&data->mold.traj, x, y, z, data->mold.mol.atom.count, f_idx);

                        md_util_backbone_angle_args_t bb_args = {
                            .atom = {
                                .count = data->mold.mol.atom.count,
                                .x = x,
                                .y = y,
                                .z = z,
                            },
                            .backbone = {
                                .count = data->mold.mol.backbone.count,
                                .atoms = data->mold.mol.backbone.atoms,
                            },
                            .chain = {
                                .count = data->mold.mol.chain.count,
                                .backbone_range = data->mold.mol.chain.backbone_range,
                            }
                        };
                        md_util_compute_backbone_angles(traj_data.backbone_angles.data + traj_data.backbone_angles.stride * f_idx, traj_data.backbone_angles.stride, &bb_args);

                        md_util_secondary_structure_args_t ss_args = {
                            .atom = {
                                .count = data->mold.mol.atom.count,
                                .x = x,
                                .y = y,
                                .z = z,
                            },
                            .backbone = {
                                .count = data->mold.mol.backbone.count,
                                .atoms = data->mold.mol.backbone.atoms,
                            },
                            .chain = {
                                .count = data->mold.mol.chain.count,
                                .backbone_range = data->mold.mol.chain.backbone_range,
                            }
                        };
                        md_util_compute_secondary_structure(traj_data.secondary_structure.data + traj_data.secondary_structure.stride * f_idx, traj_data.secondary_structure.stride, &ss_args);
                    }
                });
        }
    }
}

static void on_trajectory_load_complete(ApplicationData* data) {
    if (data->mold.traj.num_frames) {
        init_trajectory_data(data);
        //update_all_representations(data);

        //stats::set_all_property_flags(true, true);

        Array<const BackboneAngle> angle_data = {(BackboneAngle*)data->trajectory_data.backbone_angles.data, data->trajectory_data.backbone_angles.count};
        int64_t angle_stride = data->trajectory_data.backbone_angles.stride;
        ramachandran::init_vbo(angle_data, angle_stride);
    }
}

static bool load_trajectory_data(ApplicationData* data, CStringView filename) {
    interrupt_async_tasks(data);
    free_trajectory_data(data);
    data->files.trajectory = filename;
    data->animation.time = 0;
    data->animation.frame = 0;

    if (load::traj::open_file(&data->mold.traj, {.ptr = filename.cstr(), .len = filename.length()}, &data->mold.mol, default_allocator)) {
        on_trajectory_load_complete(data);
        return true;
    }

    return false;
}

static bool load_dataset_from_file(ApplicationData* data, CStringView filename) {
    ASSERT(data);

    str_t file = str_cast(filename);

    if (filename) {
        if (load::mol::is_extension_supported(file)) {
            interrupt_async_tasks(data);
            free_molecule_data(data);
            free_trajectory_data(data);
            data->files.molecule = filename;
            data->files.trajectory = "";

            LOG_NOTE("Attempting to load molecular data from file '%.*s'", file.len, file.ptr);
            if (!load::mol::load_file(&data->mold.mol, file, default_allocator)) {
                LOG_ERROR("Failed to load molecular data");
                data->files.molecule = "";
                return false;
            }
            init_molecule_data(data);

            // @NOTE: Some files contain both atomic coordinates and trajectory
            if (load::traj::is_extension_supported(file)) {
                LOG_NOTE("File also contains trajectory, attempting to load trajectory");
                load_trajectory_data(data, filename);
            }
            return true;
        } else if (load::traj::is_extension_supported(file)) {
            if (!data->mold.mol.atom.count) {
                LOG_ERROR("Before loading a trajectory, molecular data needs to be present");
                return false;
            }
            return load_trajectory_data(data, filename);
        } else {
            LOG_ERROR("File extension not supported");
        }
    }
    return false;
}

// ### WORKSPACE ###
static RepresentationType get_rep_type(CStringView str) {
    if (compare(str, "VDW"))
        return RepresentationType::Vdw;
    else if (compare(str, "LICORICE"))
        return RepresentationType::Licorice;
    else if (compare(str, "BALL_AND_STICK"))    // Ball and stick is removed for now
        return RepresentationType::Vdw;
    else if (compare(str, "RIBBONS"))
        return RepresentationType::Ribbons;
    else if (compare(str, "CARTOON"))
        return RepresentationType::Cartoon;
    else
        return RepresentationType::Vdw;
}

static CStringView get_rep_type_name(RepresentationType type) {
    switch (type) {
        case RepresentationType::Vdw:
            return "VDW";
        case RepresentationType::Licorice:
            return "LICORICE";
        /*case RepresentationType::BallAndStick:
            return "BALL_AND_STICK";*/
        case RepresentationType::Ribbons:
            return "RIBBONS";
        case RepresentationType::Cartoon:
            return "CARTOON";
        default:
            return "UNKNOWN";
    }
}

static ColorMapping get_color_mapping(CStringView str) {
    if (compare(str, "UNIFORM"))
        return ColorMapping::Uniform;
    else if (compare(str, "CPK"))
        return ColorMapping::Cpk;
    else if (compare(str, "RES_ID"))
        return ColorMapping::ResId;
    else if (compare(str, "RES_INDEX"))
        return ColorMapping::ResIndex;
    else if (compare(str, "CHAIN_ID"))
        return ColorMapping::ChainId;
    else if (compare(str, "CHAIN_INDEX"))
        return ColorMapping::ChainIndex;
    else if (compare(str, "SECONDARY_STRUCTURE"))
        return ColorMapping::SecondaryStructure;
    else
        return ColorMapping::Cpk;
}

static CStringView get_color_mapping_name(ColorMapping mapping) {
    switch (mapping) {
        case ColorMapping::Uniform:
            return "UNIFORM";
        case ColorMapping::Cpk:
            return "CPK";
        case ColorMapping::ResId:
            return "RES_ID";
        case ColorMapping::ResIndex:
            return "RES_INDEX";
        case ColorMapping::ChainId:
            return "CHAIN_ID";
        case ColorMapping::ChainIndex:
            return "CHAIN_INDEX";
        case ColorMapping::SecondaryStructure:
            return "SECONDARY_STRUCTURE";
        default:
            return "UNDEFINED";
    }
}

static vec4 to_vec4(CStringView txt, const vec4& default_val = vec4(1)) {
    vec4 res = default_val;
    DynamicArray<CStringView> tokens = tokenize(txt, ",");
    i32 count = (i32)tokens.size() < 4 ? (i32)tokens.size() : 4;
    for (int i = 0; i < count; i++) {
        res[i] = to_float(tokens[i]);
    }
    return res;
}

static void load_workspace(ApplicationData* data, CStringView file) {
    ASSERT(data);
    clear_representations(data);
    //stats::remove_all_properties();
    //clear(data->density_volume.iso.isosurfaces);

    StringBuffer<256> new_molecule_file;
    StringBuffer<256> new_trajectory_file;

    StringView txt = allocate_and_read_textfile(file);
    defer { free_string(&txt); };

    CStringView c_txt = txt;
    CStringView line;
    while ((line = extract_line(c_txt))) {
        if (compare_n(line, "[Files]", 7)) {
            while (c_txt && c_txt[0] != '[' && (line = extract_line(c_txt))) {
                if (compare_n(line, "MoleculeFile=", 13)) {
                    new_molecule_file = get_absolute_path(file, trim(line.substr(13)));
                }
                if (compare_n(line, "TrajectoryFile=", 15)) {
                    new_trajectory_file = get_absolute_path(file, trim(line.substr(15)));
                }
            }
        } else if (compare_n(line, "[Representation]", 16)) {
            Representation* rep = create_representation(data);
            if (rep) {
                while (c_txt && c_txt[0] != '[' && (line = extract_line(c_txt))) {
                    if (compare_n(line, "Name=", 5)) rep->name = trim(line.substr(5));
                    if (compare_n(line, "Filter=", 7)) rep->filter = trim(line.substr(7));
                    if (compare_n(line, "Type=", 5)) rep->type = get_rep_type(trim(line.substr(5)));
                    if (compare_n(line, "ColorMapping=", 13)) rep->color_mapping = get_color_mapping(trim(line.substr(13)));
                    if (compare_n(line, "Enabled=", 8)) rep->enabled = to_int(trim(line.substr(8))) != 0;
                    if (compare_n(line, "StaticColor=", 12)) rep->uniform_color = to_vec4(trim(line.substr(12)));
                    if (compare_n(line, "Radius=", 7)) rep->radius = to_float(trim(line.substr(7)));
                    if (compare_n(line, "Tension=", 8)) rep->tension = to_float(trim(line.substr(8)));
                    if (compare_n(line, "Width=", 6)) rep->width = to_float(trim(line.substr(6)));
                    if (compare_n(line, "Thickness=", 10)) rep->thickness = to_float(trim(line.substr(10)));
                }
            }
        /*} else if (compare_n(line, "[Property]", 10)) {
            StringBuffer<256> name, args;
            while (c_txt && c_txt[0] != '[' && (line = extract_line(c_txt))) {
                if (compare_n(line, "Name=", 5)) name = trim(line.substr(5));
                if (compare_n(line, "Args=", 5)) args = trim(line.substr(5));
            }
            //stats::create_property(name, args);
            */
        } else if (compare_n(line, "[RenderSettings]", 16)) {
            while (c_txt && c_txt[0] != '[' && (line = extract_line(c_txt))) {
                if (compare_n(line, "SsaoEnabled=", 12)) data->visuals.ssao.enabled = to_int(trim(line.substr(12))) != 0;
                if (compare_n(line, "SsaoIntensity=", 14)) data->visuals.ssao.intensity = to_float(trim(line.substr(14)));
                if (compare_n(line, "SsaoRadius=", 11)) data->visuals.ssao.radius = to_float(trim(line.substr(11)));
                if (compare_n(line, "SsaoBias=", 9)) data->visuals.ssao.bias = to_float(trim(line.substr(9)));
                if (compare_n(line, "DofEnabled=", 11)) data->visuals.dof.enabled = to_int(trim(line.substr(11))) != 0;
                if (compare_n(line, "DofFocusScale=", 14)) data->visuals.dof.focus_scale = to_float(trim(line.substr(14)));

                /*
                if (compare_n(line, "DensityVolumeEnabled=", 21)) data->density_volume.enabled = to_int(trim(line.substr(21))) != 0;
                if (compare_n(line, "DensityScale=", 13)) data->density_volume.dvr.density_scale = to_float(trim(line.substr(13)));
                if (compare_n(line, "AlphaScale=", 11)) data->density_volume.dvr.tf.alpha_scale = to_float(trim(line.substr(11)));
                if (compare_n(line, "TFFileName=", 11)) {
                    CStringView tfpath = trim(line.substr(11));
                    if (!compare(tfpath, data->density_volume.dvr.tf.path)) {
                        data->density_volume.dvr.tf.path = tfpath;
                        data->density_volume.dvr.tf.dirty = true;
                    }
                }
                if (compare_n(line, "IsoSurfaceRenderingEnabled=", 27)) data->density_volume.iso.enabled = to_int(trim(line.substr(27))) != 0;
                */
            }
        /*} else if (compare_n(line, "[Isosurface]", 21)) {
            float isovalue = -1.0f;
            vec4 isocolor{0.0f};
            while (c_txt && c_txt[0] != '[' && (line = extract_line(c_txt))) {
                if (compare_n(line, "Isovalue=", 9)) isovalue = to_float(trim(line.substr(9)));
                if (compare_n(line, "IsosurfaceColor=", 16)) isocolor = to_vec4(trim(line.substr(16)));
            }
            insert(data->density_volume.iso.isosurfaces, isovalue, isocolor);
            */
        } else if (compare_n(line, "[Camera]", 8)) {
            while (c_txt && c_txt[0] != '[' && (line = extract_line(c_txt))) {
                if (compare_n(line, "Position=", 9)) {
                    vec3 pos = vec3(to_vec4(trim(line.substr(9))));
                    data->view.camera.position = pos;
                    data->view.animation.target_position = pos;
                }
                if (compare_n(line, "Rotation=", 9)) {
                    vec4 v = to_vec4(trim(line.substr(9)));
                    data->view.camera.orientation.x = v.x;
                    data->view.camera.orientation.y = v.y;
                    data->view.camera.orientation.z = v.z;
                    data->view.camera.orientation.w = v.w;
                }
                if (compare_n(line, "Distance=", 9)) {
                    data->view.camera.focus_distance = to_float(trim(line.substr(9)));
                }
            }
        }
    }

    data->files.workspace = file;

    if (!compare(new_molecule_file, data->files.molecule) && new_molecule_file) {
        load_dataset_from_file(data, new_molecule_file);
    }

    if (!compare(new_trajectory_file, data->files.trajectory) && new_trajectory_file) {
        load_dataset_from_file(data, new_trajectory_file);
    }

    data->animation = {};
    reset_view(data, false, true);
    init_all_representations(data);
    update_all_representations(data);
}

static void save_workspace(ApplicationData* data, CStringView file) {
    FILE* fptr = fopen(file, CStringView("w"));
    if (!fptr) {
        printf("ERROR! Could not save workspace to file '%s'\n", file.beg());
        return;
    }

    fprintf(fptr, "[Files]\n");
    fprintf(fptr, "MoleculeFile=%s\n", data->files.molecule ? get_relative_path(file, data->files.molecule).cstr() : "");
    fprintf(fptr, "TrajectoryFile=%s\n", data->files.trajectory ? get_relative_path(file, data->files.trajectory).cstr() : "");
    fprintf(fptr, "\n");

    // REPRESENTATIONS
    for (const auto& rep : data->representations.buffer) {
        fprintf(fptr, "[Representation]\n");
        fprintf(fptr, "Name=%s\n", rep.name.cstr());
        fprintf(fptr, "Filter=%s\n", rep.filter.cstr());
        fprintf(fptr, "Type=%s\n", get_rep_type_name(rep.type).cstr());
        fprintf(fptr, "ColorMapping=%s\n", get_color_mapping_name(rep.color_mapping).cstr());
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

    data->files.workspace = file;
}

void create_screenshot(ApplicationData* data) {
    ASSERT(data);
    Image img;
    init_image(&img, data->fbo.width, data->fbo.height);
    defer { free_image(&img); };

    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    glReadBuffer(GL_BACK);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    glReadPixels(0, 0, img.width, img.height, GL_RGBA, GL_UNSIGNED_BYTE, img.data);

    {
        // @NOTE: Swap Rows to flip image with respect to y-axis
        const u32 row_byte_size = img.width * sizeof(u32);
        u32* row_t = (u32*)md_alloc(frame_allocator, row_byte_size);
        defer { md_free(frame_allocator, row_t, row_byte_size); };
        for (u32 i = 0; i < (u32)img.height / 2; ++i) {
            u32* row_a = img.data + i * img.width;
            u32* row_b = img.data + (img.height - 1 - i) * img.width;
            if (row_a != row_b) {
                memcpy(row_t, row_a, row_byte_size);  // tmp = a;
                memcpy(row_a, row_b, row_byte_size);  // a = b;
                memcpy(row_b, row_t, row_byte_size);  // b = tmp;
            }
        }
    }

    platform::FileDialogResult file_res = platform::file_dialog(platform::FileDialogFlags_Save, {}, "jpg;png;bmp");
    if (file_res.result == platform::FileDialogResult::Ok) {
        CStringView ext = get_file_extension(file_res.path);
        if (!ext) {
            file_res.path += ".jpg";
            ext = "jpg";
        }
        if (compare_ignore_case(ext, "jpg")) {
            const int quality = 95;
            write_image_jpg(img, file_res.path, quality);
        } else if (compare_ignore_case(ext, "png")) {
            write_image_png(img, file_res.path);
        } else if (compare_ignore_case(ext, "bmp")) {
            write_image_bmp(img, file_res.path);
        } else {
            LOG_ERROR("Supplied image format is not supported");
        }
    }
}

// #representation
static Representation* create_representation(ApplicationData* data, RepresentationType type, ColorMapping color_mapping, CStringView filter) {
    ASSERT(data);
    Representation& rep = data->representations.buffer.push_back({});
    rep.type = type;
    rep.color_mapping = color_mapping;
    rep.filter = filter;   
    init_representation(data, &rep);
    update_representation(data, &rep);
    return &rep;
}

static Representation* clone_representation(ApplicationData* data, const Representation& rep) {
    ASSERT(data);
    Representation& clone = data->representations.buffer.push_back(rep);
    init_representation(data, &clone);
    update_representation(data, &clone);
    return &clone;
}

static void remove_representation(ApplicationData* data, int idx) {
    ASSERT(data);
    ASSERT(idx < data->representations.buffer.size());
    auto& rep = data->representations.buffer[idx];
    if (rep.atom_mask) bitfield::free(&rep.atom_mask);
    md_gl_representation_free(&rep.md_rep);
    data->representations.buffer.remove(&rep);
}

static void recompute_atom_visibility_mask(ApplicationData* data) {
    ASSERT(data);

    auto& atom_visibility_mask = data->representations.atom_visibility_mask;
    if (atom_visibility_mask.size() != data->mold.mol.atom.count) {
        bitfield::init(&atom_visibility_mask, data->mold.mol.atom.count);
    }

    bitfield::clear_all(atom_visibility_mask);
    for (const auto& rep : data->representations.buffer) {
        if (!rep.enabled) continue;
        bitfield::or_field(atom_visibility_mask, atom_visibility_mask, rep.atom_mask);
    }

    data->mold.dirty_buffers |= MolBit_DirtyFlags;
}

static void update_all_representations(ApplicationData* data) {
    for (auto& rep : data->representations.buffer) {
        update_representation(data, &rep);
    }
}

static void update_representation(ApplicationData* data, Representation* rep) {
    ASSERT(data);
    ASSERT(rep);

    const int64_t bytes = data->mold.mol.atom.count * sizeof(u32);
    void* mem = md_alloc(frame_allocator, bytes);
    defer { md_free(frame_allocator, mem, bytes); };

    Array<u32> colors((u32*)mem, data->mold.mol.atom.count);
    const auto& mol = data->mold.mol;

    switch (rep->color_mapping) {
        case ColorMapping::Uniform:
            color_atoms_uniform(colors, rep->uniform_color);
            break;
        case ColorMapping::Cpk:
            color_atoms_cpk(colors, mol);
            break;
        case ColorMapping::ResId:
            color_atoms_residue_id(colors, mol);
            break;
        case ColorMapping::ResIndex:
            color_atoms_residue_index(colors, mol);
            break;
        case ColorMapping::ChainId:
            color_atoms_chain_id(colors, mol);
            break;
        case ColorMapping::ChainIndex:
            color_atoms_chain_index(colors, mol);
            break;
        case ColorMapping::SecondaryStructure:
            color_atoms_secondary_structure(colors, mol);
            break;
        default:
            ASSERT(false);
            break;
    }

    rep->filter_is_valid = filter_expression(*data, CStringView(rep->filter.buffer), rep->atom_mask);

    filter_colors(colors, rep->atom_mask);
    data->representations.atom_visibility_mask_dirty = true;

    {
        md_gl_representation_type type = MD_GL_REP_DEFAULT;
        md_gl_representation_args args = {};
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
        md_gl_representation_set_color(&rep->md_rep, 0, (uint32_t)mol.atom.count, colors.data(), 0);
    }
}

static void init_representation(ApplicationData* data, Representation* rep) {
    md_gl_representation_init(&rep->md_rep, &data->mold.gl_mol);
    if (rep->atom_mask.size() != data->mold.mol.atom.count) {
        bitfield::init(&rep->atom_mask, data->mold.mol.atom.count);
    }
}

static void init_all_representations(ApplicationData* data) {
    for (auto& rep : data->representations.buffer) {
        init_representation(data, &rep);
    }
}

static void reset_representations(ApplicationData* data) {
    ASSERT(data);
    for (auto& rep : data->representations.buffer) {
        update_representation(data, &rep);
    }
}

static void clear_representations(ApplicationData* data) {
    ASSERT(data);
    while (data->representations.buffer.size() > 0) {
        remove_representation(data, (i32)data->representations.buffer.size() - 1);
    }
}

static void update_representation_buffers(Representation* rep) {
    ASSERT(rep);

    if (rep->flags & RepBit_DirtyFilter) {

    }

    if (rep->flags & RepBit_DirtyColor) {

    }

    rep->flags = 0;
}

// #selection
static Selection* create_selection(ApplicationData* data, CStringView name, Bitfield atom_mask) {
    ASSERT(data);
    Selection sel;
    sel.name = name;
    bitfield::init(&sel.atom_mask, atom_mask);
    return &data->selection.stored_selections.push_back(sel);
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
    if (idx < 0 || (int)data->selection.stored_selections.size() <= idx) {
        LOG_ERROR("Index [%i] out of range when trying to remove selection", idx);
    }
    auto item = &data->selection.stored_selections[idx];
    if (item->atom_mask) bitfield::free(&item->atom_mask);
    data->selection.stored_selections.remove(item);
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
    static platform::Coordinate x0;
    const platform::Coordinate x1 = data->ctx.input.mouse.win_coord;
    const i64 N = data->mold.mol.atom.count;
    const bool shift_down = data->ctx.input.key.down[Key::KEY_LEFT_SHIFT] || data->ctx.input.key.down[Key::KEY_RIGHT_SHIFT];
    const bool mouse_down = data->ctx.input.mouse.down[0] || data->ctx.input.mouse.down[1];

    Bitfield mask;
    bitfield::init(&mask, N);
    defer { bitfield::free(&mask); };

    bitfield::clear_all(data->selection.current_highlight_mask);
    data->mold.dirty_buffers |= MolBit_DirtyFlags;
        
    // Range<int32> picking_range = {0, 0};

    if (data->picking.idx != NO_PICKING_IDX && !region_select) {
        ASSERT(0 <= data->picking.idx && data->picking.idx <= N);
        bitfield::set_bit(mask, data->picking.idx);

        switch (data->selection.level_mode) {
            case SelectionLevel::Atom:
                break;
            case SelectionLevel::Residue: {
                expand_mask(mask, data->mold.mol.residue.atom_range, data->mold.mol.residue.count);
                break;
            }
            case SelectionLevel::Chain: {
                expand_mask(mask, data->mold.mol.chain.atom_range, data->mold.mol.chain.count);
                break;
            }
            default:
                ASSERT(false);
                break;
        }
        memcpy(data->selection.current_highlight_mask.data(), mask.data(), mask.size_in_bytes());
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

        const ImVec2 min_p = ImVec2(math::min(x0.x, x1.x), math::min(x0.y, x1.y));
        const ImVec2 max_p = ImVec2(math::max(x0.x, x1.x), math::max(x0.y, x1.y));

        if (region_select) {
            const vec2 res = {data->ctx.window.width, data->ctx.window.height};
            const mat4 mvp = data->view.param.matrix.current.view_proj;
            const Bitfield vis_mask = data->representations.atom_visibility_mask;

            for (i64 i = 0; i < data->mold.mol.atom.count; ++i) {
                if (vis_mask[i]) {
                    const float x = data->mold.mol.atom.x[i];
                    const float y = data->mold.mol.atom.y[i];
                    const float z = data->mold.mol.atom.z[i];

                    const float p_x = mvp[0][0] * x + mvp[1][0] * y + mvp[2][0] * z + mvp[3][0];
                    const float p_y = mvp[0][1] * x + mvp[1][1] * y + mvp[2][1] * z + mvp[3][1];
                    const float p_w = mvp[0][3] * x + mvp[1][3] * y + mvp[2][3] * z + mvp[3][3];

                    vec2 c = {
                        (p_x / p_w * 0.5f + 0.5f) * res.x,
                        (-p_y / p_w * 0.5f + 0.5f) * res.y
                    };


                    if (min_p.x <= c.x && c.x <= max_p.x && min_p.y <= c.y && c.y <= max_p.y) {
                        bitfield::set_bit(mask, i);
                    }
                }
            }

            switch (data->selection.level_mode) {
                case SelectionLevel::Atom:
                    break;
                case SelectionLevel::Residue:
                    expand_mask(mask, data->mold.mol.residue.atom_range, data->mold.mol.residue.count);
                    break;
                case SelectionLevel::Chain:
                    expand_mask(mask, data->mold.mol.chain.atom_range, data->mold.mol.residue.count);
                    break;
                default:
                    ASSERT(false);
            }

            Bitfield dst_mask = mouse_down ? data->selection.current_highlight_mask : data->selection.current_selection_mask;
            Bitfield src_mask = data->selection.current_selection_mask;

            if (region_mode == RegionMode::Append) {
                bitfield::or_field(dst_mask, src_mask, mask);
            } else if (region_mode == RegionMode::Remove) {
                bitfield::and_not_field(dst_mask, src_mask, mask);
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
            if (data->picking.idx != NO_PICKING_IDX) {
                const bool append = data->ctx.input.mouse.clicked[0];
                if (append) {
                    bitfield::or_field(data->selection.current_selection_mask, data->selection.current_selection_mask, mask);
                } else {
                    bitfield::and_not_field(data->selection.current_selection_mask, data->selection.current_selection_mask, mask);
                }
            } else if (data->ctx.input.mouse.clicked[1]) {
                bitfield::clear_all(data->selection.current_selection_mask);
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
        const vec2 mouse_delta = vec2(data->ctx.input.mouse.win_delta.x, data->ctx.input.mouse.win_delta.y);
        TrackballControllerInput input;
        input.rotate_button = data->ctx.input.mouse.down[0];
        input.pan_button = data->ctx.input.mouse.down[1];
        input.dolly_button = data->ctx.input.mouse.down[2];
        input.mouse_coord_curr = {data->ctx.input.mouse.win_coord.x, data->ctx.input.mouse.win_coord.y};
        input.mouse_coord_prev = input.mouse_coord_curr - mouse_delta;
        input.screen_size = vec2(data->ctx.window.width, data->ctx.window.height);
        input.dolly_delta = data->ctx.input.mouse.scroll_delta;
        input.fov_y = data->view.camera.fov_y;

        if (camera_controller_trackball(&data->view.camera.position, &data->view.camera.orientation, &data->view.camera.focus_distance, input, data->view.trackball_param)) {
            data->view.animation.target_position = data->view.camera.position;
        }

        if (ImGui::GetIO().MouseDoubleClicked[0]) {
            if (data->picking.depth < 1.0f) {
                const vec3 forward = data->view.camera.orientation * vec3(0, 0, 1);
                const f32 dist = data->view.camera.focus_distance;
                data->view.animation.target_position = data->picking.world_coord + forward * dist;
            }
        }

        data->visuals.dof.focus_depth.target = data->view.camera.focus_distance;
    }
}

static void handle_camera_animation(ApplicationData* data) {
    const f32 dt = (f32)math::min(data->ctx.timing.delta_s, 0.033);
    {
        // #camera-translation
        constexpr f32 speed = 10.0f;
        const vec3 vel = (data->view.animation.target_position - data->view.camera.position) * speed;
        data->view.camera.position += vel * dt;
    }
    {
        //data->view.camera.orientation = math::slerp(data->view.camera.orientation, data->view.animation.target_orientation, 0.99f);
    }
    {
        // #focus-depth
        constexpr f32 speed = 10.0f;
        const f32 vel = (data->visuals.dof.focus_depth.target - data->visuals.dof.focus_depth.current) * speed;
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

static void fill_gbuffer(ApplicationData* data) {
    const vec4 CLEAR_INDEX = vec4(1, 1, 1, 1);
    const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
                                   GL_COLOR_ATTACHMENT_POST_TONEMAP, GL_COLOR_ATTACHMENT_PICKING};

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, data->fbo.deferred.fbo);
    glViewport(0, 0, data->fbo.width, data->fbo.height);

    glDepthMask(1);
    glColorMask(1, 1, 1, 1);

    // Setup fbo and clear textures
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

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    // Enable all draw buffers
    glDrawBuffers(ARRAY_SIZE(draw_buffers), draw_buffers);

    PUSH_GPU_SECTION("G-Buffer fill")

    // SIMULATION BOX
    if (data->simulation_box.enabled && data->simulation_box.box != mat3(0)) {
        PUSH_GPU_SECTION("Draw Simulation Box")

        immediate::set_model_view_matrix(data->view.param.matrix.current.view);
        immediate::set_proj_matrix(data->view.param.matrix.current.proj_jittered);

        const mat3 box = data->simulation_box.box;
        const vec3 min_box = box * vec3(0.0f);
        const vec3 max_box = box * vec3(1.0f);

        // Simulation box
        immediate::draw_box_wireframe(min_box, max_box, math::convert_color(data->simulation_box.color));

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
    postprocessing::blit_static_velocity(data->fbo.deferred.depth, data->view.param);
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
    const bool atom_selection_empty = !bitfield::any_bit_set(data->selection.current_selection_mask);
    const bool atom_highlight_empty = !bitfield::any_bit_set(data->selection.current_highlight_mask);

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
        postprocessing::scale_hsv(data->fbo.deferred.color, vec3(1, saturation, 1));
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
        vec2 coord = {data->ctx.input.mouse.win_coord.x, (float)data->fbo.height - data->ctx.input.mouse.win_coord.y};
        if (coord.x < 0.f || coord.x >= (float)data->fbo.width || coord.y < 0.f || coord.y >= (float)data->fbo.height) {
            data->picking.idx = NO_PICKING_IDX;
            data->picking.depth = 1.f;
        } else {
#if PICKING_JITTER_HACK
            static u32 frame_idx = 0;
            static u32 ref_frame = 0;
            frame_idx = (frame_idx + 1) % 16;
            // @NOTE: If we have jittering applied, we cannot? retreive the original pixel value (without the jitter)
            // Solution, pick one reference frame out of the jittering sequence and use that one...
            // Ugly hack but works...

            if (data->ctx.input.mouse.moving) {
                ref_frame = frame_idx;
            }

            if (ref_frame == frame_idx || data->view.param.jitter.current == vec2(0, 0)) {
                data->picking = read_picking_data(data->fbo, (i32)math::round(coord.x), (i32)math::round(coord.y));
                if (data->picking.idx != NO_PICKING_IDX)
                    data->picking.idx = math::clamp(data->picking.idx, 0U, (u32)data->mold.mol.atom.count - 1U);
                const vec4 viewport(0, 0, data->fbo.width, data->fbo.height);
                data->picking.world_coord =
                    math::unproject(vec3(coord.x, coord.y, data->picking.depth), data->view.param.matrix.inverse.view_proj_jittered, viewport);
            }
#else
            //coord -= data->view.param.jitter.current * 2.0f;
            //coord += 0.5f;
            data->picking = read_picking_data(data->fbo, (int)coord.x, (int)coord.y);
            const vec4 viewport(0, 0, data->fbo.width, data->fbo.height);
            data->picking.world_coord =
                math::unproject(vec3(coord.x, coord.y, data->picking.depth), data->view.param.matrix.inverse.view_proj_jittered, viewport);
#endif
        }
        data->selection.hovered = -1;
        if (data->picking.idx != NO_PICKING_IDX) {
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

    desc.input_textures.depth = data.fbo.deferred.depth;
    desc.input_textures.color = data.fbo.deferred.color;
    desc.input_textures.normal = data.fbo.deferred.normal;
    desc.input_textures.velocity = data.fbo.deferred.velocity;
    desc.input_textures.post_tonemap = data.fbo.deferred.post_tonemap;

    postprocessing::shade_and_postprocess(desc, data.view.param);
    POP_GPU_SECTION()
}

static void draw_representations(ApplicationData* data) {
    ASSERT(data);
    const md_gl_representation* rep_data[32] = { 0 };
    u32 rep_count = 0;
    for (i64 i = 0; i < data->representations.buffer.size(); i++) {
        if (data->representations.buffer[i].enabled) {
            rep_data[rep_count++] = &data->representations.buffer[i].md_rep;
        }
        if (rep_count == ARRAY_SIZE(rep_data)) break;
    }

    md_gl_rendertarget render_target = {
        .width = (u32)data->fbo.width,
        .height = (u32)data->fbo.height,
        .texture_depth = data->fbo.deferred.depth,
        .texture_color = data->fbo.deferred.color,
        .texture_atom_index = data->fbo.deferred.picking,
        .texture_view_normal = data->fbo.deferred.normal,
        .texture_view_velocity = data->fbo.deferred.velocity,
    };

    md_gl_draw_args desc = {
        .representation = {
            .count = rep_count,
            .data = rep_data,
        },
        .view_transform = {
            .model_view_matrix = &data->view.param.matrix.current.view[0][0],
            .projection_matrix = &data->view.param.matrix.current.proj_jittered[0][0],
            // These two are for temporal anti-aliasing reprojection (optional)
            .prev_model_view_matrix = &data->view.param.matrix.previous.view[0][0],
            .prev_projection_matrix = &data->view.param.matrix.previous.proj_jittered[0][0],
        },
        .render_target = &render_target,
    };

    md_draw(&data->mold.gl_ctx, &desc);
}

static void draw_representations_lean_and_mean(ApplicationData* data, u32 mask) {
    const md_gl_representation* rep_data[32] = { 0 };
    u32 rep_count = 0;
    for (i64 i = 0; i < data->representations.buffer.size(); i++) {
        if (data->representations.buffer[i].enabled) {
            rep_data[rep_count++] = &data->representations.buffer[i].md_rep;
        }
        if (rep_count == ARRAY_SIZE(rep_data)) break;
    }

    md_gl_rendertarget render_target = {
        .width = (u32)data->fbo.width,
        .height = (u32)data->fbo.height,
        .texture_depth = data->fbo.deferred.depth,
    };

    md_gl_draw_args desc = {
        .representation = {
            .count = rep_count,
            .data = rep_data,
        },
        .view_transform = {
            .model_view_matrix = &data->view.param.matrix.current.view[0][0],
            .projection_matrix = &data->view.param.matrix.current.proj_jittered[0][0],
            // These two are for temporal anti-aliasing reprojection
            //.prev_model_view_matrix = &data->view.param.matrix.previous.view[0][0],
            //.prev_projection_matrix = &data->view.param.matrix.previous.proj_jittered[0][0],
        },
        .render_target = &render_target,
        .mol_mask = mask,
    };

    md_draw(&data->mold.gl_ctx, &desc);
}
