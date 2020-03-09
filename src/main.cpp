#include <core/types.h>
#include <core/file.h>
#include <core/hash.h>
#include <core/log.h>
#include <core/bitfield.h>
#include <core/math_utils.h>
#include <core/string_utils.h>
#include <core/spatial_hash.h>
#include <core/sync.h>

#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>
#include <mol/trajectory_utils.h>
#include <mol/molecule_utils.h>
#include <mol/hydrogen_bond.h>
#include <mol/filter.h>
#include <mol/pdb_utils.h>
#include <mol/gro_utils.h>
#include <mol/xtc_utils.h>

#include <imgui.h>
#define IMGUI_DEFINE_MATH_OPERATORS
#include <imgui_internal.h>

//#include <TaskScheduler.h>

#include "gfx/gl.h"
#include "gfx/camera.h"
#include "gfx/camera_utils.h"
#include "gfx/molecule_draw.h"
#include "gfx/immediate_draw_utils.h"
#include "gfx/postprocessing_utils.h"
#include "gfx/volumerender_utils.h"
#include "gfx/conetracing_utils.h"

#include "volume.h"
#include "range_slider.h"
#include "plot_extended.h"
#include "platform/platform.h"
#include "console.h"
#include "stats.h"
#include "ramachandran.h"
#include "color_utils.h"
#include "structure_tracking.h"
#include "isosurface.h"
#include "task_system.h"
#include "trajectory_loader.h"

#include <stdio.h>
#include <thread>
#include <mutex>

#define PICKING_JITTER_HACK 1
#define DEPERIODIZE_ON_FRAME_LOADED 1
#define SHOW_IMGUI_DEMO_WINDOW 0
#define VIAMD_RELEASE 0
#define EXPERIMENTAL_CULLING 0
#define EXPERIMENTAL_SDF 0

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
#define PUSH_GPU_SECTION(lbl)                                                                       \
    {                                                                                               \
        if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); \
    }
#define POP_GPU_SECTION()                       \
    {                                           \
        if (glPopDebugGroup) glPopDebugGroup(); \
    }

constexpr unsigned int NO_PICKING_IDX = 0xFFFFFFFFU;
constexpr CStringView FILE_EXTENSION = "via";

constexpr u32 DEL_BTN_COLOR = 0xFF1111CC;
constexpr u32 DEL_BTN_HOVER_COLOR = 0xFF3333DD;
constexpr u32 DEL_BTN_ACTIVE_COLOR = 0xFF5555FF;
constexpr u32 TEXT_BG_ERROR_COLOR = 0xAA222299;

constexpr f32 HYDROGEN_BOND_DISTANCE_CUTOFF_DEFAULT = 3.0f;
constexpr f32 HYDROGEN_BOND_DISTANCE_CUTOFF_MIN = 0.1f;
constexpr f32 HYDROGEN_BOND_DISTANCE_CUTOFF_MAX = 12.0f;

constexpr f32 HYDROGEN_BOND_ANGLE_CUTOFF_DEFAULT = 20.f;
constexpr f32 HYDROGEN_BOND_ANGLE_CUTOFF_MIN = 5.f;
constexpr f32 HYDROGEN_BOND_ANGLE_CUTOFF_MAX = 90.f;

constexpr f32 BALL_AND_STICK_VDW_SCALE = 0.30f;
constexpr f32 BALL_AND_STICK_LICORICE_SCALE = 0.5f;

constexpr i32 SPLINE_SUBDIVISION_COUNT = 16;

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

enum class PlaybackMode { Stopped, Playing };
enum class InterpolationMode { Nearest, Linear, Cubic };
enum class SelectionLevel { Atom, Residue, Chain };
enum class SelectionOperator { Or, And };
enum class SelectionGrowth { CovalentBond, Radial };
enum class RepresentationType { Vdw, Licorice, BallAndStick, Ribbons, Cartoon };
enum class TrackingMode { Absolute, Relative };
enum class CameraMode { Perspective, Orthographic };

enum AtomBit_ { AtomBit_Highlighted = BIT(0), AtomBit_Selected = BIT(1), AtomBit_Visible = BIT(2) };

// #struct Structure Declarations

struct PickingData {
    u32 idx = NO_PICKING_IDX;
    f32 depth = 1.0f;
    vec3 world_coord = {0, 0, 0};
};

struct MainFramebuffer {
    struct {
        GLuint depth = 0;
        GLuint color = 0;
        GLuint normal = 0;
        GLuint velocity = 0;
        GLuint emissive = 0;
        GLuint post_tonemap = 0;
        GLuint picking = 0;
        GLuint fbo = 0;
    } deferred;

    struct {
        // @NOTE: two of each for ping-pong read / write (hopefully non blocking reads)
        GLuint color[2] = {0, 0};
        GLuint depth[2] = {0, 0};
    } pbo_picking;

    int width = 0;
    int height = 0;
};

struct MoleculeBuffers {
    GLuint position = 0;
    GLuint old_position = 0;
    GLuint velocity = 0;
    GLuint radius = 0;
    GLuint selection = 0;
    GLuint bond = 0;

    struct {
        GLuint residue = 0;
        GLuint aabb = 0;
        GLuint visibility = 0;
    } experimental;

    struct {
        GLuint backbone_segment_index = 0;  // Stores indices to atoms which are needed for control points and support vectors
        GLuint control_point = 0;           // Stores extracted control points[vec3] and support vectors[vec3] before subdivision
        GLuint control_point_index =
            0;  // Stores draw element indices to compute splines with adjacent index information (each chain separated by restart-index 0xFFFFFFFFU).
        GLuint spline = 0;        // Stores subdivided spline data control points[vec3] + support vector[vec3]
        GLuint spline_index = 0;  // Stores draw element indices to render splines (each chain separated by restart-index 0xFFFFFFFFU).

        i32 num_backbone_segment_indices = 0;
        i32 num_control_point_indices = 0;
        i32 num_spline_indices = 0;
    } backbone;

    struct {
        bool position = false;
        bool selection = false;
        bool backbone = false;
    } dirty;
};

struct Representation {
    StringBuffer<64> name = "rep";
    StringBuffer<128> filter = "all";
    RepresentationType type = RepresentationType::Vdw;
    ColorMapping color_mapping = ColorMapping::Cpk;
    Bitfield atom_mask{};
    GLuint color_buffer = 0;

    bool enabled = true;
    bool show_in_selection = true;
    bool filter_is_ok = true;

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

struct ReferenceFrame {
    StringBuffer<64> name = "ref";
    StringBuffer<128> filter = "label CA";
    Bitfield atom_mask{};
    bool active = false;
    bool show = false;
    bool filter_is_ok = false;
    structure_tracking::ID id = 0;
    TrackingMode tracking_mode = TrackingMode::Relative;

    mat4 world_to_reference = {};
    mat4 reference_to_world = {};

    // For debugging
    mat4 basis = {};
    mat3 pca = {};
};

struct EnsembleStructure {
    int offset = 0;
    structure_tracking::ID id = 0;
    mat4 alignment_matrix = {};
    mat4 world_to_structure = {};
    bool enabled = false;
};

/*
struct ThreadSyncData {
    std::thread thread{};
    i32 running = 0;
    i32 stop_signal = 0;

    void signal_stop() { stop_signal = true; }

    void wait_until_finished() {
        while (running) {
            platform::sleep(1);
        }
    }

    void signal_stop_and_wait() {
        signal_stop();
        wait_until_finished();
    }
};
*/

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
        TrackballControllerState trackball_state{};
        ViewParam param;
        CameraMode mode = CameraMode::Perspective;

        struct {
            vec3 target_position{};
            // quat target_orientation{};
        } animation;
    } view;

    // --- MOLECULAR DATA ---
    MoleculeDynamic dynamic{};

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
        bool show_property_window = false;
        bool show_timeline_window = false;
        bool show_distribution_window = false;
    } statistics;

    // --- FRAMEBUFFER ---
    MainFramebuffer fbo;

    PickingData picking;

    // --- MOLECULE GPU BUFFERS ---
    MoleculeBuffers gpu_buffers;

    // --- PLAYBACK ---
    struct {
        f64 time = 0.f;  // double precision for long trajectories
        i32 frame = 0;
        f32 fps = 10.f;
        InterpolationMode interpolation = InterpolationMode::Cubic;
        PlaybackMode mode = PlaybackMode::Stopped;
        bool apply_pbc = false;
    } playback;

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

        struct {
            bool enabled = true;
            f32 intensity = 1.0f;
            f32 voxel_ext = 1.0f;
        } cone_traced_ao;

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

    struct {
        bool enabled = false;
        bool dirty = true;
        bool overlay = false;
        vec4 color = vec4(1, 0, 1, 1);
        f32 distance_cutoff = HYDROGEN_BOND_DISTANCE_CUTOFF_DEFAULT;  // In Ångström
        f32 angle_cutoff = HYDROGEN_BOND_ANGLE_CUTOFF_DEFAULT;        // In Degrees
        DynamicArray<HydrogenBond> bonds{};
    } hydrogen_bonds;

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
        } texture;

        struct {
            vec3 min = {0, 0, 0};
            vec3 max = {1, 1, 1};
        } clip_planes;

        struct {
            structure_tracking::ID reference_id = 0;
            bool enable_view = false;
        } reference_structure;

        Volume volume{};
        std::mutex volume_data_mutex{};

        // mat4 model_to_world_matrix{};
        // mat4 world_to_model_matrix{};
        vec3 voxel_spacing{1.0f};
        float resolution_scale = 2.0f;

    } density_volume;

    struct {
        cone_trace::GPUVolume vol;
        std::mutex vol_mutex{};
    } occupancy_volume;

    struct {
        Bitfield atom_mask;
        DynamicArray<EnsembleStructure> structures;
        bool superimpose_structures = false;
        TrackingMode tracking_mode = TrackingMode::Absolute;
    } ensemble_tracking;

    // --- RAMACHANDRAN ---
    struct {
        bool show_window = false;
        // int frame_range_min = 0;
        // int frame_range_max = 0;

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
                f32 radius = 1.5f;
                vec4 fill_color = vec4(1, 1, 1, 1);
            } base;
            struct {
                f32 radius = 3.0f;
                vec4 selection_color = vec4(0.8f, 1.0f, 0.5f, 1.0f);
                vec4 highlight_color = vec4(0.8f, 1.0f, 0.5f, 0.5f);
            } selection;
        } current;

        // BackboneAnglesTrajectory backbone_angles{};
    } ramachandran;

    // --- REPRESENTATIONS ---
    struct {
        DynamicArray<Representation> buffer = {};
        Bitfield atom_visibility_mask{};
        bool atom_visibility_mask_dirty = false;
        bool show_window = false;
        // bool changed = false;
    } representations;

    // --- REFERENCE FRAME ---
    struct {
        DynamicArray<ReferenceFrame> frames = {};
        bool show_window = false;
    } reference_frame;

    struct {
        i32 target = 0;
        bool show_window = false;
    } shape_space;

    // --- CONSOLE ---
    Console console{};
    bool show_console = false;
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

// From https://github.com/procedural/gpulib/blob/master/gpulib_imgui.h
struct ImVec3 {
    float x, y, z;
    ImVec3(float _x = 0.0f, float _y = 0.0f, float _z = 0.0f) {
        x = _x;
        y = _y;
        z = _z;
    }
};

void imgui_easy_theming(ImVec3 color_for_text, ImVec3 color_for_head, ImVec3 color_for_area, ImVec3 color_for_body, ImVec3 color_for_pops) {
    ImGuiStyle& style = ImGui::GetStyle();

    style.Colors[ImGuiCol_Text] = ImVec4(color_for_text.x, color_for_text.y, color_for_text.z, 1.00f);
    style.Colors[ImGuiCol_TextDisabled] = ImVec4(color_for_text.x, color_for_text.y, color_for_text.z, 0.58f);
    style.Colors[ImGuiCol_WindowBg] = ImVec4(color_for_body.x, color_for_body.y, color_for_body.z, 0.95f);
    style.Colors[ImGuiCol_ChildWindowBg] = ImVec4(color_for_area.x, color_for_area.y, color_for_area.z, 0.58f);
    style.Colors[ImGuiCol_Border] = ImVec4(color_for_body.x, color_for_body.y, color_for_body.z, 0.00f);
    style.Colors[ImGuiCol_BorderShadow] = ImVec4(color_for_body.x, color_for_body.y, color_for_body.z, 0.00f);
    style.Colors[ImGuiCol_FrameBg] = ImVec4(color_for_area.x, color_for_area.y, color_for_area.z, 1.00f);
    style.Colors[ImGuiCol_FrameBgHovered] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 0.78f);
    style.Colors[ImGuiCol_FrameBgActive] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 1.00f);
    style.Colors[ImGuiCol_TitleBg] = ImVec4(color_for_area.x, color_for_area.y, color_for_area.z, 1.00f);
    style.Colors[ImGuiCol_TitleBgCollapsed] = ImVec4(color_for_area.x, color_for_area.y, color_for_area.z, 0.75f);
    style.Colors[ImGuiCol_TitleBgActive] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 1.00f);
    style.Colors[ImGuiCol_MenuBarBg] = ImVec4(color_for_area.x, color_for_area.y, color_for_area.z, 0.47f);
    style.Colors[ImGuiCol_ScrollbarBg] = ImVec4(color_for_area.x, color_for_area.y, color_for_area.z, 1.00f);
    style.Colors[ImGuiCol_ScrollbarGrab] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 0.21f);
    style.Colors[ImGuiCol_ScrollbarGrabHovered] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 0.78f);
    style.Colors[ImGuiCol_ScrollbarGrabActive] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 1.00f);
    // style.Colors[ImGuiCol_ComboBg] = ImVec4(color_for_area.x, color_for_area.y, color_for_area.z, 1.00f);
    style.Colors[ImGuiCol_CheckMark] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 0.80f);
    style.Colors[ImGuiCol_SliderGrab] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 0.50f);
    style.Colors[ImGuiCol_SliderGrabActive] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 1.00f);
    style.Colors[ImGuiCol_Button] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 0.50f);
    style.Colors[ImGuiCol_ButtonHovered] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 0.86f);
    style.Colors[ImGuiCol_ButtonActive] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 1.00f);
    style.Colors[ImGuiCol_Header] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 0.76f);
    style.Colors[ImGuiCol_HeaderHovered] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 0.86f);
    style.Colors[ImGuiCol_HeaderActive] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 1.00f);
    style.Colors[ImGuiCol_Column] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 0.32f);
    style.Colors[ImGuiCol_ColumnHovered] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 0.78f);
    style.Colors[ImGuiCol_ColumnActive] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 1.00f);
    style.Colors[ImGuiCol_ResizeGrip] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 0.15f);
    style.Colors[ImGuiCol_ResizeGripHovered] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 0.78f);
    style.Colors[ImGuiCol_ResizeGripActive] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 1.00f);
    // style.Colors[ImGuiCol_CloseButton] = ImVec4(color_for_text.x, color_for_text.y, color_for_text.z, 0.16f);
    // style.Colors[ImGuiCol_CloseButtonHovered] = ImVec4(color_for_text.x, color_for_text.y, color_for_text.z, 0.39f);
    // style.Colors[ImGuiCol_CloseButtonActive] = ImVec4(color_for_text.x, color_for_text.y, color_for_text.z, 1.00f);
    style.Colors[ImGuiCol_PlotLines] = ImVec4(color_for_text.x, color_for_text.y, color_for_text.z, 0.63f);
    style.Colors[ImGuiCol_PlotLinesHovered] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 1.00f);
    style.Colors[ImGuiCol_PlotHistogram] = ImVec4(color_for_text.x, color_for_text.y, color_for_text.z, 0.63f);
    style.Colors[ImGuiCol_PlotHistogramHovered] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 1.00f);
    style.Colors[ImGuiCol_TextSelectedBg] = ImVec4(color_for_head.x, color_for_head.y, color_for_head.z, 0.43f);
    style.Colors[ImGuiCol_PopupBg] = ImVec4(color_for_pops.x, color_for_pops.y, color_for_pops.z, 0.92f);
    style.Colors[ImGuiCol_ModalWindowDarkening] = ImVec4(color_for_area.x, color_for_area.y, color_for_area.z, 0.73f);
}

void SetupImGuiStyle2() {
    static ImVec3 color_for_text = ImVec3(236.f / 255.f, 240.f / 255.f, 241.f / 255.f);
    static ImVec3 color_for_head = ImVec3(41.f / 255.f, 128.f / 255.f, 185.f / 255.f);
    static ImVec3 color_for_area = ImVec3(57.f / 255.f, 79.f / 255.f, 105.f / 255.f);
    static ImVec3 color_for_body = ImVec3(44.f / 255.f, 62.f / 255.f, 80.f / 255.f);
    static ImVec3 color_for_pops = ImVec3(33.f / 255.f, 46.f / 255.f, 60.f / 255.f);
    imgui_easy_theming(color_for_text, color_for_head, color_for_area, color_for_body, color_for_pops);
}

}  // namespace ImGui

static vec3 compute_simulation_box_ext(const ApplicationData& data);
static void interpolate_atomic_positions(ApplicationData* data);
static void reset_view(ApplicationData* data, bool move_camera = false, bool smooth_transition = false);
static f32 compute_avg_ms(f32 dt);
static PickingData read_picking_data(const MainFramebuffer& fbo, i32 x, i32 y);

static bool handle_selection(ApplicationData* data);
// static void handle_animation(ApplicationData* data);
static void handle_camera_interaction(ApplicationData* data);
static void handle_camera_animation(ApplicationData* data);

static void update_properties(ApplicationData* data);
static void update_density_volume_texture(ApplicationData* data);
static void handle_picking(ApplicationData* data);
static void compute_backbone_spline(ApplicationData* data);
static void compute_velocity(ApplicationData* data);
static void compute_residue_aabbs(ApplicationData* data);
static void cull_residue_aabbs(ApplicationData* data);
static void fill_gbuffer(const ApplicationData& data);
static void apply_postprocessing(const ApplicationData& data);
static void draw_residue_aabbs(const ApplicationData& data);
static void init_occupancy_volume(ApplicationData* data);

static void draw_representations(const ApplicationData& data);
static void draw_representations_lean_and_mean(const ApplicationData& data, vec4 color = vec4(1, 1, 1, 1), float scale = 1.0f,
                                               u32 mask = 0xFFFFFFFFU);

static void draw_main_menu(ApplicationData* data);
static void draw_context_popup(ApplicationData* data);
static void draw_animation_control_window(ApplicationData* data);
static void draw_representations_window(ApplicationData* data);
static void draw_property_window(ApplicationData* data);
static void draw_timeline_window(ApplicationData* data);
static void draw_distribution_window(ApplicationData* data);
static void draw_ramachandran_window(ApplicationData* data);
static void draw_atom_info_window(const MoleculeStructure& mol, int atom_idx, int x, int y);
static void draw_molecule_dynamic_info_window(ApplicationData* data);
static void draw_async_task_window(ApplicationData* data);
static void draw_reference_frame_window(ApplicationData* data);
static void draw_shape_space_window(ApplicationData* data);
static void draw_density_volume_window(ApplicationData* data);
// static void draw_density_volume_clip_plane_widgets(ApplicationData* data);
// static void draw_selection_window(ApplicationData* data);

static void init_framebuffer(MainFramebuffer* fbo, int width, int height);
static void destroy_framebuffer(MainFramebuffer* fbo);

static void init_molecule_buffers(ApplicationData* data);
static void free_molecule_buffers(ApplicationData* data);

static void clear_buffer(GLuint buffer);
static void copy_buffer(GLuint dst_buffer, GLuint src_buffer);
static void copy_molecule_data_to_buffers(ApplicationData* data);

static void init_molecule_data(ApplicationData* data);
static void init_trajectory_data(ApplicationData* data);

static bool load_dataset_from_file(ApplicationData* data, CStringView file);
static void free_dataset(ApplicationData* data);

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
static void reset_representations(ApplicationData* data);
static void clear_representations(ApplicationData* data);

static void recompute_atom_visibility_mask(ApplicationData* data);

// Selections
static Selection* create_selection(ApplicationData* data, CStringView name, const Bitfield);
static Selection* clone_selection(ApplicationData* data, const Selection& sel);
static void remove_selection(ApplicationData* data, int idx);

static void reset_selections(ApplicationData* data);
static void clear_selections(ApplicationData* data);

// Reference frames
static ReferenceFrame* create_reference_frame(ApplicationData* data, CStringView name = "ref", CStringView filter = "label CA");
static bool remove_reference_frame(ApplicationData* data, ReferenceFrame* ref);
static void reset_reference_frame(ApplicationData* data);
static void clear_reference_frames(ApplicationData* data);
static void compute_reference_frame(ApplicationData* data, ReferenceFrame* ref);
static ReferenceFrame* get_active_reference_frame(ApplicationData* data);
static void update_reference_frames(ApplicationData* data);

static void superimpose_ensemble(ApplicationData* data);

static void init_density_volume(ApplicationData* data);

// Async operations
// static void load_trajectory_async(ApplicationData* data);

// static void on_trajectory_frame_loaded(ApplicationData* data, int frame_idx);
static void on_trajectory_load_complete(ApplicationData* data);

// Temp functionality

int main(int, char**) {
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

    // Init subsystems
    LOG_NOTE("Initializing immediate draw...");
    immediate::initialize();
    LOG_NOTE("Initializing molecule draw...");
    draw::initialize();
    LOG_NOTE("Initializing ramachandran...");
    ramachandran::initialize(data.dynamic);
    LOG_NOTE("Initializing stats...");
    stats::initialize();
    LOG_NOTE("Initializing filter...");
    filter::initialize();
    LOG_NOTE("Initializing post processing...");
    postprocessing::initialize(data.fbo.width, data.fbo.height);
    LOG_NOTE("Initializing volume...");
    volume::initialize();
    LOG_NOTE("Initializing cone tracing...");
    cone_trace::initialize();
    LOG_NOTE("Initializing structure tracking...");
    structure_tracking::initialize();
    LOG_NOTE("Initializing task system...");
    task_system::initialize();

    // ImGui::SetupImGuiStyle2();
    ImGui::StyleColorsLight();
    // ImGui::StyleColorsClassic();

    vec2 halton_sequence[16];
    math::generate_halton_sequence(halton_sequence, ARRAY_SIZE(halton_sequence), 2, 3);

#if VIAMD_RELEASE
    pdb::load_molecule_from_string(&data.dynamic.molecule, CAFFINE_PDB);
    init_molecule_data(&data);
#else
    load_dataset_from_file(&data, "D:/data/3u75.pdb");
    // load_dataset_from_file(&data, VIAMD_DATA_DIR "/1ALA-250ns-2500frames.pdb");
    // load_dataset_from_file(&data, VIAMD_DATA_DIR "/amyloid/centered.gro");
    // load_dataset_from_file(&data, VIAMD_DATA_DIR "/amyloid/centered.xtc");
    // load_dataset_from_file(&data, "D:/data/md/6T-water/16-6T-box.gro");
    // load_dataset_from_file(&data, "D:/data/md/6T-water/16-6T-box-md-npt.xtc");
    // load_dataset_from_file(&data, "D:/data/md/p-ftaa-water/p-FTAA-box-sol-ions-em.gro");
    // load_dataset_from_file(&data, "D:/data/md/p-ftaa-water/p-FTAA-box-sol-ions-md-npt.xtc");

    // create_reference_frame(&data, "dna", "dna");
    // create_reference_frame(&data, "res1", "residue 1");
    // data.simulation_box.enabled = true;
    // stats::create_property("d1", "distance atom(*) com(atom(*))");
    // data.density_volume.enabled = true;
    // data.reference_frame.frames[0].active = true;
#endif
    reset_view(&data, true);
    create_representation(&data, RepresentationType::Vdw, ColorMapping::ResId, "not water");
    // create_representation(&data, RepresentationType::Vdw, ColorMapping::ResId, "residue 1");

    init_density_volume(&data);

    interpolate_atomic_positions(&data);
    copy_molecule_data_to_buffers(&data);
    copy_buffer(data.gpu_buffers.old_position, data.gpu_buffers.position);

    init_occupancy_volume(&data);

#if EXPERIMENTAL_SDF == 1
    draw::scan::test_scan();
#endif

    // Main loop
    while (!data.ctx.window.should_close) {
        platform::update(&data.ctx);

#if SHOW_IMGUI_DEMO_WINDOW
        ImGui::ShowDemoWindow();
#endif

        const i32 num_frames = data.dynamic.trajectory ? data.dynamic.trajectory.num_frames : 0;
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
                draw::initialize();
                postprocessing::initialize(data.fbo.width, data.fbo.height);
                volume::initialize();
                cone_trace::initialize();
            }

            if (data.ctx.input.key.hit[KEY_PLAY_PAUSE]) {
                if (data.playback.mode == PlaybackMode::Stopped && data.playback.time == max_time) {
                    data.playback.time = 0;
                }

                switch (data.playback.mode) {
                    case PlaybackMode::Playing:
                        data.playback.mode = PlaybackMode::Stopped;
                        break;
                    case PlaybackMode::Stopped:
                        data.playback.mode = PlaybackMode::Playing;
                        break;
                    default:
                        ASSERT(false);
                }
            }

            if (data.ctx.input.key.hit[KEY_SKIP_TO_PREV_FRAME] || data.ctx.input.key.hit[KEY_SKIP_TO_NEXT_FRAME]) {
                f64 step = data.ctx.input.key.down[Key::KEY_LEFT_CONTROL] ? 10.0 : 1.0;
                if (data.ctx.input.key.hit[KEY_SKIP_TO_PREV_FRAME]) step = -step;
                data.playback.time = math::clamp(data.playback.time + step, 0.0, max_time);
            }
        }

        recompute_atom_visibility_mask(&data);

        data.selection.selecting = false;
        if (!ImGui::GetIO().WantCaptureMouse) {
            data.selection.selecting |= handle_selection(&data);
        }
        // if (data.selection.selecting) update_all_representations(&data);

        handle_camera_interaction(&data);
        handle_camera_animation(&data);

        // This needs to happen first (in imgui events) to enable docking of imgui windows
        // ImGui::CreateDockspace();

        if (data.playback.mode == PlaybackMode::Playing) {
            data.playback.time += data.ctx.timing.delta_s * data.playback.fps;
            data.playback.time = math::clamp(data.playback.time, 0.0, max_time);
            if (data.playback.time >= max_time) {
                data.playback.mode = PlaybackMode::Stopped;
                data.playback.time = max_time;
            }
        }

        const int new_frame = math::clamp((int)(data.playback.time + 0.5f), 0, last_frame);
        const bool frame_changed = new_frame != data.playback.frame;
        data.playback.frame = new_frame;

        bool time_changed = false;
        {
            static auto prev_time = data.playback.time;
            if (data.playback.time != prev_time) {
                time_changed = true;
                prev_time = data.playback.time;
            }
        }

        if (data.time_filter.dynamic_window) {
            const auto half_window_ext = data.time_filter.window_extent * 0.5f;
            data.time_filter.range.min = math::clamp((f32)data.playback.time - half_window_ext, 0.0f, (f32)max_time);
            data.time_filter.range.max = math::clamp((f32)data.playback.time + half_window_ext, 0.0f, (f32)max_time);

            if (frame_changed && data.dynamic.trajectory) {
                stats::set_all_property_flags(false, true);
            }
        }

        if (time_changed) {
            auto& mol = data.dynamic.molecule;
            auto& traj = data.dynamic.trajectory;

            data.hydrogen_bonds.dirty = true;
            data.gpu_buffers.dirty.backbone = true;

            PUSH_CPU_SECTION("Interpolate Position")
            if (traj) {
                interpolate_atomic_positions(&data);
                update_reference_frames(&data);

                if (data.ensemble_tracking.superimpose_structures) {
                    superimpose_ensemble(&data);
                }

                if (data.playback.apply_pbc) {
                    apply_pbc(mol.atom.position.x, mol.atom.position.y, mol.atom.position.z, mol.atom.mass, mol.sequences, data.simulation_box.box);
                }

                init_occupancy_volume(&data);

#if 0
                if (const auto ref = get_active_reference_frame(&data)) {
                    if (const auto tracking_data = structure_tracking::get_tracking_data(ref->id)) {
                        const vec3 com = ref->basis[3];
                        int idx = 0;
                        for (int i = 0; i < ref->atom_mask.size(); i++) {
                            if (ref->atom_mask[i]) {
                                data.dynamic.molecule.atom.position.x[i] = com.x + tracking_data->average_structure.x[idx];
                                data.dynamic.molecule.atom.position.y[i] = com.y + tracking_data->average_structure.y[idx];
                                data.dynamic.molecule.atom.position.z[i] = com.z + tracking_data->average_structure.z[idx];
                                idx++;
                            }
                        }
                    }
                }
#endif

                data.gpu_buffers.dirty.position = true;
            }
            POP_CPU_SECTION()

            PUSH_CPU_SECTION("Compute backbone angles") {
                zero_array(mol.backbone.angles);
                for (const auto& bb_seq : mol.backbone.sequences) {
                    auto segments = get_backbone(mol, bb_seq);
                    compute_backbone_angles(mol.backbone.angles.data(), segments.data(), mol.atom.position.x, mol.atom.position.y,
                                            mol.atom.position.z, segments.size());
                }
            }
            POP_CPU_SECTION()

            PUSH_CPU_SECTION("Update dynamic representations")
            for (auto& rep : data.representations.buffer) {
                if (rep.color_mapping == ColorMapping::SecondaryStructure) {
                    update_representation(&data, &rep);
                }
            }
            POP_CPU_SECTION()
        } else {
            static auto prev_time = data.playback.time;
            if (data.playback.time != prev_time) {
                copy_buffer(data.gpu_buffers.old_position, data.gpu_buffers.position);
            }
            prev_time = data.playback.time;
        }

        PUSH_CPU_SECTION("Hydrogen bonds")
        if (data.hydrogen_bonds.enabled && data.hydrogen_bonds.dirty) {
            const auto& mol = data.dynamic.molecule;
            data.hydrogen_bonds.bonds = hydrogen_bond::compute_bonds(mol.hydrogen_bond.donors, mol.hydrogen_bond.acceptors, mol.atom.position.x,
                                                                     mol.atom.position.y, mol.atom.position.z, data.hydrogen_bonds.distance_cutoff,
                                                                     math::deg_to_rad(data.hydrogen_bonds.angle_cutoff));
            data.hydrogen_bonds.dirty = false;
        }
        POP_CPU_SECTION()

        // Resize Framebuffer
        if ((data.fbo.width != data.ctx.framebuffer.width || data.fbo.height != data.ctx.framebuffer.height) &&
            (data.ctx.framebuffer.width != 0 && data.ctx.framebuffer.height != 0)) {
            init_framebuffer(&data.fbo, data.ctx.framebuffer.width, data.ctx.framebuffer.height);
            postprocessing::initialize(data.fbo.width, data.fbo.height);
        }

        {
            ViewParam& param = data.view.param;
            param.matrix.previous = param.matrix.current;
            param.jitter.previous = param.jitter.current;

            param.clip_planes.near = data.view.camera.near_plane;
            param.clip_planes.far = data.view.camera.far_plane;
            param.fov_y = data.view.camera.fov_y;

            param.matrix.current.view = compute_world_to_view_matrix(data.view.camera);
            if (data.view.mode == CameraMode::Perspective) {
                param.matrix.current.proj = perspective_projection_matrix(data.view.camera, data.fbo.width, data.fbo.height);
            } else {
                const float aspect_ratio = (float)data.fbo.width / (float)data.fbo.height;
                const float h = data.view.trackball_state.distance * math::tan(data.view.camera.fov_y * 0.5f);
                const float w = aspect_ratio * h;
                param.matrix.current.proj = orthographic_projection_matrix(-w, w, -h, h, data.view.camera.near_plane, data.view.camera.far_plane);
            }
            param.matrix.current.proj_jittered = param.matrix.current.proj;

            param.resolution = vec2(data.fbo.width, data.fbo.height);
            if (data.visuals.temporal_reprojection.enabled && data.visuals.temporal_reprojection.jitter) {
                static u32 i = 0;
                i = (++i) % ARRAY_SIZE(halton_sequence);
                param.jitter.current = halton_sequence[i] - 0.5f;
                if (data.view.mode == CameraMode::Perspective) {
                    param.matrix.current.proj_jittered = perspective_projection_matrix(data.view.camera, data.fbo.width, data.fbo.height,
                                                                                       param.jitter.current.x, param.jitter.current.y);
                } else {
                    const float aspect_ratio = (float)data.fbo.width / (float)data.fbo.height;
                    const float h = data.view.trackball_state.distance * math::tan(data.view.camera.fov_y * 0.5f);
                    const float w = aspect_ratio * h;
                    const float scale_x = w / data.fbo.width * 2.0f;
                    const float scale_y = h / data.fbo.height * 2.0f;
                    const float j_x = param.jitter.current.x * scale_x;
                    const float j_y = param.jitter.current.y * scale_y;
                    param.matrix.current.proj_jittered = param.matrix.current.proj =
                        orthographic_projection_matrix(-w + j_x, w + j_x, -h + j_y, h + j_y, data.view.camera.near_plane, data.view.camera.far_plane);
                }
            }

            if (const ReferenceFrame* ref = get_active_reference_frame(&data)) {
                param.matrix.current.view = param.matrix.current.view * ref->world_to_reference;
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

        copy_molecule_data_to_buffers(&data);
#if EXPERIMENTAL_SDF == 1
        PUSH_GPU_SECTION("COMPUTE VDW SDF");
        const AABB aabb = compute_aabb(data.dynamic.molecule.atom.position.x, data.dynamic.molecule.atom.position.y,
                                       data.dynamic.molecule.atom.position.z, data.dynamic.molecule.atom.radius, data.dynamic.molecule.atom.count);
        float max_radius = 0.0f;
        for (int64 i = 0; i < data.dynamic.molecule.atom.count; i++) {
            max_radius = math::max(data.dynamic.molecule.atom.radius[i], max_radius);
        }
        draw::sdf::compute_vdw_sdf(data.gpu_buffers.position, data.gpu_buffers.radius, data.dynamic.molecule.atom.count, aabb, max_radius);
        // draw::sdf::compute_vdw_sdf(data.dynamic.molecule.atom.position.x, data.dynamic.molecule.atom.position.y,
        // data.dynamic.molecule.atom.position.z, data.dynamic.molecule.atom.radius, data.dynamic.molecule.atom.count, aabb);
        POP_GPU_SECTION();
#endif

        update_properties(&data);
        update_density_volume_texture(&data);
        compute_backbone_spline(&data);
        compute_velocity(&data);
#if EXPERIMENTAL_CULLING == 1
        compute_residue_aabbs(&data);
#endif
        fill_gbuffer(data);
#if EXPERIMENTAL_CULLING == 1
        cull_residue_aabbs(&data);
#endif
        if (data.visuals.cone_traced_ao.enabled) {
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, data.fbo.deferred.fbo);
            glViewport(0, 0, data.fbo.width, data.fbo.height);
            glDrawBuffer(GL_COLOR_ATTACHMENT0);  // Modify to color buffer
            cone_trace::render_directional_occlusion(data.fbo.deferred.depth, data.fbo.deferred.normal, data.occupancy_volume.vol,
                                                     data.view.param.matrix.current.view, data.view.param.matrix.current.proj,
                                                     data.visuals.cone_traced_ao.intensity);
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        }
#if EXPERIMENTAL_SDF == 1
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, data.fbo.deferred.fbo);
        glViewport(0, 0, data.fbo.width, data.fbo.height);
        glDrawBuffer(GL_COLOR_ATTACHMENT4);
        glEnable(GL_DEPTH_TEST);
        // draw::sdf::draw_sdf_debug(data.view.param);
        glDisable(GL_DEPTH_TEST);
#endif
        handle_picking(&data);

        // Activate backbuffer
        glDisable(GL_DEPTH_TEST);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glViewport(0, 0, data.ctx.framebuffer.width, data.ctx.framebuffer.height);
        glDrawBuffer(GL_BACK);
        glClear(GL_COLOR_BUFFER_BIT);

        apply_postprocessing(data);
#if EXPERIMENTAL_SDF == 1
        /*
        glDisable(GL_DEPTH_TEST);
        draw::sdf::draw_sdf(data.view.param);
        */
#endif

        // GUI ELEMENTS
        data.console.Draw("VIAMD", data.ctx.window.width, data.ctx.window.height, data.ctx.timing.delta_s);

        draw_main_menu(&data);
        draw_context_popup(&data);
        draw_async_task_window(&data);
        draw_animation_control_window(&data);
        draw_molecule_dynamic_info_window(&data);

        if (data.representations.show_window) draw_representations_window(&data);
        if (data.statistics.show_property_window) draw_property_window(&data);
        if (data.statistics.show_timeline_window) draw_timeline_window(&data);
        if (data.statistics.show_distribution_window) draw_distribution_window(&data);
        if (data.ramachandran.show_window) draw_ramachandran_window(&data);
        if (data.reference_frame.show_window) draw_reference_frame_window(&data);
        if (data.shape_space.show_window) draw_shape_space_window(&data);
        if (data.density_volume.show_window) draw_density_volume_window(&data);

        // @NOTE: ImGui::GetIO().WantCaptureMouse does not work with Menu
        if (!ImGui::IsMouseHoveringAnyWindow()) {
            if (data.picking.idx != NO_PICKING_IDX) {
                draw_atom_info_window(data.dynamic.molecule, data.picking.idx, (int)data.ctx.input.mouse.win_coord.x,
                                      (int)data.ctx.input.mouse.win_coord.y);
            }
        }

        PUSH_GPU_SECTION("Imgui render")
        platform::render_imgui(&data.ctx);
        POP_GPU_SECTION()

        // Swap buffers
        platform::swap_buffers(&data.ctx);

        task_system::clear_completed_tasks();
    }

    stats::signal_stop_and_wait();

    // shutdown subsystems
    LOG_NOTE("Shutting down immediate draw...");
    immediate::shutdown();
    LOG_NOTE("Shutting down molecule draw...");
    draw::shutdown();
    LOG_NOTE("Shutting down ramachandran...");
    ramachandran::shutdown();
    LOG_NOTE("Shutting down stats...");
    stats::shutdown();
    LOG_NOTE("Shutting down filter...");
    filter::shutdown();
    LOG_NOTE("Shutting down post processing...");
    postprocessing::shutdown();
    LOG_NOTE("Shutting down volume...");
    volume::shutdown();
    LOG_NOTE("Shutting down cone tracing...");
    cone_trace::shutdown();
    LOG_NOTE("Shutting down structure tracking...");
    structure_tracking::shutdown();
    LOG_NOTE("Shutting down task system...");
    task_system::shutdown();

    destroy_framebuffer(&data.fbo);
    platform::shutdown(&data.ctx);

    return 0;
}

static void interpolate_atomic_positions(ApplicationData* data) {
    ASSERT(data);
    const auto& mol = data->dynamic.molecule;
    const auto& traj = data->dynamic.trajectory;

    if (!mol || !traj) return;

    const int last_frame = math::max(0, data->dynamic.trajectory.num_frames - 1);
    const f64 time = math::clamp(data->playback.time, 0.0, f64(last_frame));

    const f32 t = (float)math::fract(time);
    const int frame = (int)time;
    const int nearest_frame = math::clamp((int)(time + 0.5), 0, last_frame);
    const int prev_frame_2 = math::max(0, frame - 1);
    const int prev_frame_1 = math::max(0, frame);
    const int next_frame_1 = math::min(frame + 1, last_frame);
    const int next_frame_2 = math::min(frame + 2, last_frame);

    const mat3 boxes[4] = {
        get_trajectory_frame(traj, prev_frame_2).box,
        get_trajectory_frame(traj, prev_frame_1).box,
        get_trajectory_frame(traj, next_frame_1).box,
        get_trajectory_frame(traj, next_frame_2).box,
    };

    const InterpolationMode mode = (prev_frame_1 != next_frame_1) ? data->playback.interpolation : InterpolationMode::Nearest;

    mat3 box;
    switch (mode) {
        case InterpolationMode::Nearest:
            box = boxes[1];
            break;
        case InterpolationMode::Linear:
            box = math::lerp(boxes[1], boxes[2], t);
            break;
        case InterpolationMode::Cubic:
            box = math::cubic_spline(boxes[0], boxes[1], boxes[2], boxes[3], t);
            break;
        default:
            ASSERT(false);
    }

    const bool pbc = (box != mat3(0));
    data->simulation_box.box = box;

    switch (mode) {
        case InterpolationMode::Nearest: {
            const auto x = get_trajectory_position_x(traj, nearest_frame);
            const auto y = get_trajectory_position_y(traj, nearest_frame);
            const auto z = get_trajectory_position_z(traj, nearest_frame);
            memcpy(mol.atom.position.x, x.data(), x.size_in_bytes());
            memcpy(mol.atom.position.y, y.data(), y.size_in_bytes());
            memcpy(mol.atom.position.z, z.data(), z.size_in_bytes());
            break;
        }
        case InterpolationMode::Linear: {
            const float* x[2] = {get_trajectory_position_x(traj, prev_frame_1).data(), get_trajectory_position_x(traj, next_frame_1).data()};
            const float* y[2] = {get_trajectory_position_y(traj, prev_frame_1).data(), get_trajectory_position_y(traj, next_frame_1).data()};
            const float* z[2] = {get_trajectory_position_z(traj, prev_frame_1).data(), get_trajectory_position_z(traj, next_frame_1).data()};

            if (pbc) {
                linear_interpolation_pbc(mol.atom.position.x, mol.atom.position.y, mol.atom.position.z, x[0], y[0], z[0], x[1], y[1], z[1],
                                         mol.atom.count, t, box);
            } else {
                linear_interpolation(mol.atom.position.x, mol.atom.position.y, mol.atom.position.z, x[0], y[0], z[0], x[1], y[1], z[1],
                                     mol.atom.count, t);
            }

            break;
        }
        case InterpolationMode::Cubic: {
            const float* x[4] = {get_trajectory_position_x(traj, prev_frame_2).data(), get_trajectory_position_x(traj, prev_frame_1).data(),
                                 get_trajectory_position_x(traj, next_frame_1).data(), get_trajectory_position_x(traj, next_frame_2).data()};
            const float* y[4] = {get_trajectory_position_y(traj, prev_frame_2).data(), get_trajectory_position_y(traj, prev_frame_1).data(),
                                 get_trajectory_position_y(traj, next_frame_1).data(), get_trajectory_position_y(traj, next_frame_2).data()};
            const float* z[4] = {get_trajectory_position_z(traj, prev_frame_2).data(), get_trajectory_position_z(traj, prev_frame_1).data(),
                                 get_trajectory_position_z(traj, next_frame_1).data(), get_trajectory_position_z(traj, next_frame_2).data()};
            if (pbc) {
                cubic_interpolation_pbc(mol.atom.position.x, mol.atom.position.y, mol.atom.position.z, x[0], y[0], z[0], x[1], y[1], z[1], x[2], y[2],
                                        z[2], x[3], y[3], z[3], mol.atom.count, t, box);
            } else {
                cubic_interpolation(mol.atom.position.x, mol.atom.position.y, mol.atom.position.z, x[0], y[0], z[0], x[1], y[1], z[1], x[2], y[2],
                                    z[2], x[3], y[3], z[3], mol.atom.count, t);
            }
            break;
        }
        default:
            ASSERT(false);
    }
}

// #misc
static f32 compute_avg_ms(f32 dt) {
    // @NOTE: Perhaps this can be done with a simple running mean?
    constexpr f32 interval = 0.5f;
    static f32 avg = 0.f;
    static int num_frames = 0;
    static f32 t = 0;
    t += dt;
    num_frames++;

    if (t > interval) {
        avg = t / num_frames * 1000.f;
        t = 0;
        num_frames = 0;
    }

    return avg;
}

static void reset_view(ApplicationData* data, bool move_camera, bool smooth_transition) {
    ASSERT(data);
    if (!data->dynamic.molecule) return;
    const auto& mol = data->dynamic.molecule;

    const AABB aabb = compute_aabb(mol.atom.position.x, mol.atom.position.y, mol.atom.position.z, mol.atom.count);
    const vec3 ext = aabb.max - aabb.min;
    const vec3 cen = (aabb.min + aabb.max) * 0.5f;
    const vec3 pos = cen + ext * 3.f;

    if (move_camera) {
        if (!smooth_transition) data->view.camera.position = pos;
        data->view.animation.target_position = pos;
        data->view.trackball_state.distance = math::length(pos - cen);
        data->view.camera.orientation = math::conjugate(math::quat_cast(look_at(data->view.animation.target_position, cen, vec3(0, 1, 0))));
    }

    data->view.camera.near_plane = 1.f;
    data->view.camera.far_plane = math::length(ext) * 50.f;
    data->view.trackball_state.params.max_distance = math::length(ext) * 20.0f;
}

// #picking
static PickingData read_picking_data(const MainFramebuffer& framebuffer, i32 x, i32 y) {
    static u32 curr = 0;
    u32 prev = (curr + 1) % 2;

    PickingData data{};

    PUSH_GPU_SECTION("READ PICKING DATA")
    glBindFramebuffer(GL_READ_FRAMEBUFFER, framebuffer.deferred.fbo);
    glReadBuffer(GL_COLOR_ATTACHMENT5);

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

    curr = prev;
    return data;
}

// #imgui
namespace ImGui {

#if 0
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
#endif

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

void grow_mask_by_covalent_bond(Bitfield mask, Array<const Bond> bonds, i32 extent) {
    Bitfield prev_mask;
    bitfield::init(&prev_mask, mask);
    defer { bitfield::free(&prev_mask); };

    for (i32 i = 0; i < extent; i++) {
        for (const auto& bond : bonds) {
            const i32 idx[2] = {bond.idx[0], bond.idx[1]};
            if (bitfield::get_bit(prev_mask, idx[0]) && !bitfield::get_bit(mask, idx[1])) {
                bitfield::set_bit(mask, idx[1]);
            } else if (bitfield::get_bit(prev_mask, idx[1]) && !bitfield::get_bit(mask, idx[0])) {
                bitfield::set_bit(mask, idx[0]);
            }
        }
        memcpy(prev_mask.data(), mask.data(), prev_mask.size_in_bytes());
    }
}

void grow_mask_by_radial_extent(Bitfield mask, const float* atom_x, const float* atom_y, const float* atom_z, i64 atom_count, float extent) {
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

template <typename AtomRangeOwner>
void expand_mask(Bitfield mask, Array<AtomRangeOwner> sequences) {
    for (const auto& seq : sequences) {
        const AtomRange range = seq.atom_range;
        if (bitfield::any_bit_set_in_range(mask, range)) {
            bitfield::set_range(mask, range);
        }
    }
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
                    if (data->representations.buffer.size() > 0) {
                        reset_representations(data);
                    } else {
                        create_representation(data);
                    }
                    data->playback = {};
                    stats::clear_all_properties();
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
                        ImGui::SliderFloat("Motion Scale", &data->visuals.temporal_reprojection.motion_blur.motion_scale, 0.f, 1.5f);
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

            // Cone Trace
            ImGui::BeginGroup();
            ImGui::PushID("Cone Trace");
            ImGui::Checkbox("Cone Traced AO", &data->visuals.cone_traced_ao.enabled);
            if (data->visuals.cone_traced_ao.enabled) {
                ImGui::SliderFloat("Intensity", &data->visuals.cone_traced_ao.intensity, 0.01f, 5.f);
            }
            ImGui::PopID();
            ImGui::EndGroup();
            ImGui::Separator();

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

            /*
            ImGui::BeginGroup();
            ImGui::Text("Spatial Selection");
            ImGui::PushID("spatial");
            ImGui::ColorEdit4("Selection Fill", &data->selection.color.selection.fill_color[0], ImGuiColorEditFlags_NoInputs |
            ImGuiColorEditFlags_Float); ImGui::ColorEdit4("Selection Outline", &data->selection.color.selection.outline_color[0],
            ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_Float); ImGui::SliderFloat("Selection Outline Scale",
            &data->selection.color.selection.outline_scale, 1.01f, 2.0f);

            ImGui::ColorEdit4("Highlight Fill", &data->selection.color.highlight.fill_color[0], ImGuiColorEditFlags_NoInputs |
            ImGuiColorEditFlags_Float); ImGui::ColorEdit4("Highlight Outline", &data->selection.color.highlight.outline_color[0],
            ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_Float); ImGui::SliderFloat("Highlight Outline Scale",
            &data->selection.color.highlight.outline_scale, 1.01f, 2.0f);

            ImGui::Spacing();
            ImGui::SliderFloat("Non selected saturation", &data->selection.color.selection_saturation, 0.0f, 1.0f);
            ImGui::PopID();
            ImGui::EndGroup();
            ImGui::Separator();
            */

            /*
            ImGui::BeginGroup();
            ImGui::Text("Ramachandran Selection");
            ImGui::PushID("rama");
            ImGui::ColorEdit4("Selection Fill", &data->ramachandran.current.selection.selection_color[0], ImGuiColorEditFlags_NoInputs |
            ImGuiColorEditFlags_Float); ImGui::ColorEdit4("Highlight Fill", &data->ramachandran.current.selection.highlight_color[0],
            ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_Float); ImGui::PopID(); ImGui::EndGroup(); ImGui::Separator();

            ImGui::BeginGroup();
            ImGui::Checkbox("Draw Control Points", &data->visuals.spline.draw_control_points);
            ImGui::Checkbox("Draw Spline", &data->visuals.spline.draw_spline);
            ImGui::EndGroup();
            ImGui::Separator();
            */

            /*
            // Property Overlay
            ImGui::BeginGroup();
            ImGui::Text("Property Style");
            auto style = stats::get_style();
            ImGui::Text("point colors ");
            for (int32 i = 0; i < style->NUM_COLORS; i++) {
                ImGui::SameLine();
                ImGui::PushID(i);
                ImVec4 color = ImColor(style->point_colors[i]);
                if (ImGui::ColorEdit4("PointColor", (float*)&color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel))
            style->point_colors[i] = ImColor(color); ImGui::PopID();
            }
            ImGui::Text("line color   ");
            ImGui::SameLine();
            ImVec4 color = ImColor(style->line_color);
            if (ImGui::ColorEdit4("LineColor", (float*)&color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel)) style->line_color =
            ImColor(color); ImGui::EndGroup(); ImGui::Separator();
            */

            ImGui::BeginGroup();
            ImGui::Checkbox("Hydrogen Bond", &data->hydrogen_bonds.enabled);
            if (data->hydrogen_bonds.enabled) {
                ImGui::PushID("hydrogen_bond");
                if (ImGui::SliderFloat("Distance Cutoff", &data->hydrogen_bonds.distance_cutoff, HYDROGEN_BOND_DISTANCE_CUTOFF_MIN,
                                       HYDROGEN_BOND_DISTANCE_CUTOFF_MAX)) {
                    data->hydrogen_bonds.dirty = true;
                }
                if (ImGui::SliderFloat("Angle Cutoff", &data->hydrogen_bonds.angle_cutoff, HYDROGEN_BOND_ANGLE_CUTOFF_MIN,
                                       HYDROGEN_BOND_ANGLE_CUTOFF_MAX)) {
                    data->hydrogen_bonds.dirty = true;
                }
                ImGui::Checkbox("Overlay", &data->hydrogen_bonds.overlay);
                ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
                ImGui::ColorEdit4("Color", (float*)&data->hydrogen_bonds.color, ImGuiColorEditFlags_NoInputs);
                ImGui::PopStyleVar();
                ImGui::PopID();
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
            ImGui::Checkbox("Properties", &data->statistics.show_property_window);
            ImGui::Checkbox("Timelines", &data->statistics.show_timeline_window);
            ImGui::Checkbox("Distributions", &data->statistics.show_distribution_window);
            ImGui::Checkbox("Ramachandran", &data->ramachandran.show_window);
            // ImGui::Checkbox("Selection", &data->selection.show_window);
            ImGui::Checkbox("Reference Frames", &data->reference_frame.show_window);
            ImGui::Checkbox("Density Volume", &data->density_volume.show_window);
            ImGui::Checkbox("Shape Space", &data->shape_space.show_window);

            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Selection")) {
            const auto atom_count = data->dynamic.molecule.atom.count;
            Bitfield mask;
            bitfield::init(&mask, atom_count);
            defer { bitfield::free(&mask); };

            if (ImGui::MenuItem("Clear Selection")) {
                bitfield::clear_all(data->selection.current_selection_mask);
                data->gpu_buffers.dirty.selection = true;
            }

            if (ImGui::MenuItem("Invert Selection")) {
                bitfield::invert_all(data->selection.current_selection_mask);
                data->gpu_buffers.dirty.selection = true;
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
                    DynamicArray<StoredSelection> sel;
                    sel.push_back({"current", data->selection.current_selection_mask});
                    for (const auto& s : data->selection.stored_selections) {
                        sel.push_back({s.name, s.atom_mask});
                    }

                    query_ok = filter::compute_filter_mask(mask, buf, data->dynamic.molecule, sel);

                    if (query_ok) {
                        switch (data->selection.level_mode) {
                            case SelectionLevel::Atom:
                                // No need to expand the mask
                                break;
                            case SelectionLevel::Residue:
                                expand_mask(mask, data->dynamic.molecule.residues);
                                break;
                            case SelectionLevel::Chain:
                                expand_mask(mask, data->dynamic.molecule.sequences);
                                break;
                            default:
                                ASSERT(false);
                        }
                    } else {
                        bitfield::clear_all(mask);
                    }
                    // data->gpu_buffers.dirty.selection = true;
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
                    data->gpu_buffers.dirty.selection = true;
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
                    data->gpu_buffers.dirty.selection = true;
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
                            grow_mask_by_covalent_bond(mask, get_covalent_bonds(data->dynamic.molecule), (i32)extent);
                            break;
                        case SelectionGrowth::Radial: {
                            const auto& mol = data->dynamic.molecule;
                            grow_mask_by_radial_extent(mask, mol.atom.position.x, mol.atom.position.y, mol.atom.position.z, mol.atom.count, extent);
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
                            expand_mask(mask, data->dynamic.molecule.residues);
                            break;
                        case SelectionLevel::Chain:
                            expand_mask(mask, data->dynamic.molecule.sequences);
                            break;
                        default:
                            ASSERT(false);
                    }

                    memcpy(data->selection.current_highlight_mask.data(), mask.data(), mask.size_in_bytes());
                    data->gpu_buffers.dirty.selection = true;

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
                            data->gpu_buffers.dirty.selection = true;
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
                        data->gpu_buffers.dirty.selection = true;
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
            const f32 ms = compute_avg_ms(data->ctx.timing.delta_s);
            char fps_buf[32];
            snprintf(fps_buf, 32, "%.2f ms (%.1f fps)", ms, 1000.f / ms);
            const float w = ImGui::CalcTextSize(fps_buf).x;
            ImGui::SetCursorPosX(ImGui::GetWindowContentRegionMax().x - w);
            ImGui::Text("%s", fps_buf);
        }
        ImGui::EndMainMenuBar();
    }

    if (new_clicked) ImGui::OpenPopup("Warning New");

    /*
            // This does not work atm, probably an Imgui Issue / Bug
if (ImGui::BeginPopupModal("Warning New")) {
    ImGui::Text("By creating a new workspace you will loose any unsaved progress.");
    ImGui::Text("Are you sure?");
    ImGui::Separator();
    if (ImGui::Button("OK", ImVec2(120, 0))) {
            ImGui::CloseCurrentPopup();
    }
    ImGui::SameLine();
    if (ImGui::Button("Cancel", ImVec2(120, 0))) {
            ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
}
    */
}

/*
void draw_selection_window(ApplicationData* data) {
    ASSERT(data);
    if (!data->selection.show_window) return;

    static DynamicArray<bool> mask;
    mask.resize(data->dynamic.molecule.atom.count);

    ImGui::Begin("Selection", &data->selection.show_window);

    // MODE
    ImGui::Combo("Selection Mode", (int*)(&data->selection.level_mode), "Atom\0Residue\0Chain\0\0");

    ImGui::Spacing();
    ImGui::Separator();

    // QUERY
    {
        enum CombineMode { AND, OR };
        static char buf[256] = {0};
        static CombineMode mode = OR;
        static bool query_ok = false;

        ImGui::Text("Query");
        const auto TEXT_BG_DEFAULT_COLOR = ImGui::ColorConvertFloat4ToU32(ImGui::GetStyle().Colors[ImGuiCol_FrameBg]);
        ImGui::PushStyleColor(ImGuiCol_FrameBg, query_ok ? TEXT_BG_DEFAULT_COLOR : TEXT_BG_ERROR_COLOR);
        bool query_modified = ImGui::InputText("##query", buf, ARRAY_SIZE(buf), ImGuiInputTextFlags_AutoSelectAll);
        bool pressed_enter = ImGui::IsItemActivePreviousFrame() && !ImGui::IsItemActive() &&
ImGui::IsKeyPressed(ImGui::GetIO().KeyMap[ImGuiKey_Enter]); ImGui::PopStyleColor();

        if (ImGui::IsWindowAppearing()) {
            ImGui::SetKeyboardFocusHere(-1);
        }

        if (ImGui::RadioButton("OR", mode == OR)) mode = OR;
        ImGui::SameLine();
        if (ImGui::RadioButton("AND", mode == AND)) mode = AND;

        if (!query_ok) {
            ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
            ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
        }
        bool apply = (ImGui::Button("Ok") || pressed_enter);
        bool show = (ImGui::GetFocusID() == ImGui::GetID("##query") || ImGui::GetHoveredID() == ImGui::GetID("Ok"));

        if (!query_ok) {
            ImGui::PopItemFlag();
            ImGui::PopStyleVar();
        }

        ImGui::SameLine();
        if (ImGui::Button("Cancel")) {
            data->selection.show_window = false;
        }

        if (show) {
            if (query_modified) {
                query_ok = filter::compute_filter_mask(mask, data->dynamic, buf);
                if (!query_ok) {
                    memset_array(mask, false);
                }
                data->gpu_buffers.dirty.selection = true;
            }

            if (query_ok) {
                for (int64 i = 0; i < data->selection.current_highlight.size(); i++) {
                    const bool mask_val = mask[i];
                    const bool curr_val = data->selection.current_selection[i];

                    if (mode == AND) {
                        data->selection.current_highlight[i] = curr_val & mask_val;
                    } else if (mode == OR) {
                        data->selection.current_highlight[i] = curr_val | mask_val;
                    }
                }
                data->gpu_buffers.dirty.selection = true;
            }
        }

        if (apply) {
            data->selection.show_window = false;

            for (int64 i = 0; i < data->selection.current_selection.size(); i++) {
                const bool mask_val = mask[i];
                const bool curr_val = data->selection.current_selection[i];

                if (mode == AND) {
                    data->selection.current_selection[i] = curr_val & mask_val;
                } else if (mode == OR) {
                    data->selection.current_selection[i] = curr_val | mask_val;
                }
            }
            data->gpu_buffers.dirty.selection = true;
        }
    }

    ImGui::Spacing();
    ImGui::Separator();

    // GROW
    {
        enum GrowMode { CovalentBond, Radial };
        static GrowMode grow_mode = CovalentBond;
        static bool pos_dir = true;
        static bool neg_dir = true;
        static float extent = 1.f;

        ImGui::Text("Grow");
        ImGui::Combo("Mode", (int*)(&grow_mode), "Covalent Bond\0Radial\0\0");
        if (grow_mode == CovalentBond) {
            // ImGui::Text("Direction: ");
            // ImGui::SameLine();
            // ImGui::Checkbox("Positive", &pos_dir);
            // ImGui::SameLine();
            // ImGui::Checkbox("Negative", &neg_dir);
        }
        ImGui::SliderFloat("Extent", &extent, 1.0f, 10.f);
        bool apply = ImGui::Button("Apply");
        bool show = (ImGui::GetHoveredID() == ImGui::GetID("Extent") || ImGui::GetActiveID() == ImGui::GetID("Extent") || ImGui::GetHoveredID() ==
ImGui::GetID("Apply"));

        if (show) {
            Array<bool> prev_mask = {(bool*)TMP_MALLOC(mask.size_in_bytes()), mask.size()};
            defer { TMP_FREE(prev_mask.data()); };

            memset_array(mask, false);
            memcpy(prev_mask.data(), data->selection.current_selection.data(), data->selection.current_selection.size_in_bytes());

            if (grow_mode == CovalentBond) {
                for (int64 i = 0; i < (int64)extent; i++) {
                    for (const auto& bond : data->dynamic.molecule.covalent_bonds) {
                        const int32 idx[2] = {bond.idx[0], bond.idx[1]};
                        if (prev_mask[idx[0]] && !mask[idx[1]]) {
                            mask[idx[1]] = true;
                        } else if (prev_mask[idx[1]] && !mask[idx[0]]) {
                            mask[idx[0]] = true;
                        }
                    }
                    memcpy(prev_mask.data(), mask.data(), prev_mask.size_in_bytes());
                }
            } else if (grow_mode == Radial) {
                const auto positions = get_positions(data->dynamic.molecule);
                spatialhash::Frame frame = spatialhash::compute_frame(positions, vec3(extent));
                for (int64 i = 0; i < positions.size(); i++) {
                    if (prev_mask[i]) {
                        spatialhash::for_each_within(frame, positions[i], extent, [mask = mask.data()](int32 idx, const vec3& pos) {
                            if (!mask[idx]) {
                                mask[idx] = true;
                            }
                        });
                    }
                }
            }

            switch (data->selection.level_mode) {
                case SelectionLevel::Atom:
                    break;
                case SelectionLevel::Residue:
                    expand_mask_to_residue(mask, data->dynamic.molecule.residues);
                    break;
                case SelectionLevel::Chain:
                    expand_mask_to_chain(mask, data->dynamic.molecule.chains);
                    break;
                default:
                    ASSERT(false);
            }

            memcpy(data->selection.current_highlight.data(), mask.data(), mask.size_in_bytes());
            data->gpu_buffers.dirty.selection = true;

            if (apply) {
                memcpy(data->selection.current_selection.data(), mask.data(), mask.size_in_bytes());
            }
        }
    }

    ImGui::Text("Stored Selections");
    if (ImGui::Button("Store Active Selection")) {
        char name_buf[64];
        snprintf(name_buf, 64, "selection%i", (int)data->selection.stored_selections.size());
        create_selection(data, name_buf, data->selection.current_selection);
    }
    for (int i = 0; i < data->selection.stored_selections.size(); i++) {
        auto& sel = data->selection.stored_selections[i];
        const float32 item_width = math::clamp(ImGui::GetWindowContentRegionWidth() - 90.f, 100.f, 300.f);
        StringBuffer<128> name;
        snprintf(name, name.size(), "%s###ID", sel.name.cstr());

        ImGui::PushID(i);
        if (ImGui::CollapsingHeader(name)) {
            if (ImGui::Button("Activate")) {
                memcpy(data->selection.current_selection.data(), sel.atom_mask.data(), sel.atom_mask.size_in_bytes());
                data->gpu_buffers.dirty.selection = true;
            }
            ImGui::SameLine();
            if (ImGui::DeleteButton("Remove")) {
                remove_selection(data, i);
            }
            ImGui::SameLine();
            if (ImGui::Button("Clone")) {
                clone_selection(data, sel);
            }
        }

        const auto h_id = ImGui::GetHoveredID();
        bool show = (h_id == ImGui::GetID(name) || h_id == ImGui::GetID("Activate") || h_id == ImGui::GetID("Remove") || h_id ==
ImGui::GetID("Clone"));

        ImGui::PopID();

        if (show) {
            memcpy(data->selection.current_highlight.data(), sel.atom_mask.data(), sel.atom_mask.size_in_bytes());
            data->gpu_buffers.dirty.selection = true;
        }
    }

    ImGui::End();
}
*/

// static void highlight_range(ApplicationData* data, AtomRange range) {}

void draw_context_popup(ApplicationData* data) {
    ASSERT(data);

    const bool shift_down = data->ctx.input.key.down[Key::KEY_LEFT_SHIFT] || data->ctx.input.key.down[Key::KEY_RIGHT_SHIFT];
    if (data->ctx.input.mouse.clicked[1] && !shift_down && !ImGui::GetIO().WantTextInput) {
        if (data->selection.right_clicked != -1 && data->dynamic) {
            ImGui::OpenPopup("AtomContextPopup");
        }
    }

    if (ImGui::BeginPopup("AtomContextPopup")) {
        if (data->selection.right_clicked != -1 && data->dynamic) {
            if (ImGui::BeginMenu("Recenter Trajectory...")) {
                const int atom_idx = data->selection.right_clicked;
                AtomRange atom_range = {};

                if (ImGui::MenuItem("on Atom")) {
                    atom_range = {atom_idx, atom_idx + 1};
                }
                if (ImGui::IsItemHovered()) {
                }
                if (ImGui::MenuItem("on Residue")) {
                    const auto res_idx = data->dynamic.molecule.atom.res_idx[atom_idx];
                    atom_range = data->dynamic.molecule.residues[res_idx].atom_range;
                }
                if (ImGui::MenuItem("on Sequence")) {
                    const auto seq_idx = data->dynamic.molecule.atom.seq_idx[atom_idx];
                    atom_range = data->dynamic.molecule.sequences[seq_idx].atom_range;
                }

                if (atom_range.ext() > 0) {
                    Bitfield atom_mask;
                    bitfield::init(&atom_mask, data->dynamic.molecule.atom.count);
                    defer { bitfield::free(&atom_mask); };
                    bitfield::set_range(atom_mask, atom_range);
                    recenter_trajectory(&data->dynamic, atom_mask);
                    interpolate_atomic_positions(data);
                    data->gpu_buffers.dirty.position = true;
                    ImGui::CloseCurrentPopup();
                }

                ImGui::EndMenu();
            }
        }
        ImGui::EndPopup();
    }
}

static void draw_animation_control_window(ApplicationData* data) {
    ASSERT(data);
    if (!data->dynamic.trajectory) return;

    ImGui::Begin("Animation");
    const i32 num_frames = data->dynamic.trajectory.num_frames;
    ImGui::Text("Num Frames: %i", num_frames);
    // ImGui::Checkbox("Apply post-interpolation pbc", &data->playback.apply_pbc);
    f32 t = (float)data->playback.time;
    if (ImGui::SliderFloat("Time", &t, 0, (float)(math::max(0, num_frames - 1)))) {
        data->playback.time = t;
    }
    ImGui::SliderFloat("fps", &data->playback.fps, 0.1f, 1000.f, "%.3f", 4.f);
    ImGui::Combo("type", (int*)(&data->playback.interpolation), "Nearest\0Linear\0Cubic\0\0");
    switch (data->playback.mode) {
        case PlaybackMode::Playing:
            if (ImGui::Button("Pause")) data->playback.mode = PlaybackMode::Stopped;
            break;
        case PlaybackMode::Stopped:
            if (ImGui::Button("Play")) data->playback.mode = PlaybackMode::Playing;
            break;
        default:
            ASSERT(false);
    }
    ImGui::SameLine();
    if (ImGui::Button("Stop")) {
        data->playback.mode = PlaybackMode::Stopped;
        data->playback.time = 0.0;
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
        bool recompute_colors = false;
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
            if (!rep.filter_is_ok) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
            if (ImGui::InputText("filter", rep.filter.cstr(), rep.filter.capacity(), ImGuiInputTextFlags_EnterReturnsTrue)) {
                recompute_colors = true;
            }
            if (!rep.filter_is_ok) ImGui::PopStyleColor();
            ImGui::Combo("type", (int*)(&rep.type), "VDW\0Licorice\0Ball & Stick\0Ribbons\0Cartoon\0\0");
            if (ImGui::Combo("color mapping", (int*)(&rep.color_mapping),
                             "Uniform Color\0CPK\0Res Id\0Res Idx\0Chain Id\0Chain Idx\0Secondary Structure\0\0")) {
                recompute_colors = true;
            }
            ImGui::PopItemWidth();
            if (rep.color_mapping == ColorMapping::Uniform) {
                ImGui::SameLine();
                if (ImGui::ColorEdit4("color", (float*)&rep.uniform_color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel)) {
                    recompute_colors = true;
                }
            }
            ImGui::PushItemWidth(item_width);
            if (rep.type == RepresentationType::Vdw || rep.type == RepresentationType::Licorice) {
                ImGui::SliderFloat("radii scale", &rep.radius, 0.1f, 2.f);
            }
            if (rep.type == RepresentationType::Ribbons) {
                ImGui::SliderFloat("spline tension", &rep.tension, 0.f, 1.f);
                ImGui::SliderFloat("spline width", &rep.width, 0.1f, 2.f);
                ImGui::SliderFloat("spline thickness", &rep.thickness, 0.1f, 2.f);
            }
            ImGui::PopItemWidth();
            ImGui::Spacing();
            ImGui::Separator();
        }

        ImGui::PopID();

        if (recompute_colors) {
            update_representation(data, &rep);
        }
    }

    ImGui::End();

    // const auto new_hash = hash::crc64(data->representations.buffer);
    // data->representations.changed = (new_hash != old_hash);
}

static void draw_reference_frame_window(ApplicationData* data) {
    ImGui::Begin("Reference Frames", &data->reference_frame.show_window, ImGuiWindowFlags_NoFocusOnAppearing);
    if (ImGui::Button("create new")) {
        create_reference_frame(data);
    }
    ImGui::SameLine();
    if (ImGui::DeleteButton("remove all")) {
        clear_reference_frames(data);
    }
    ImGui::Spacing();
    ImGui::Separator();
    for (int i = 0; i < data->reference_frame.frames.size(); i++) {
        bool recompute_frame = false;
        auto& ref = data->reference_frame.frames[i];
        const f32 item_width = math::clamp(ImGui::GetWindowContentRegionWidth() - 90.f, 100.f, 300.f);
        StringBuffer<128> name;
        snprintf(name.cstr(), name.capacity(), "%s###ID", ref.name.buffer);

        ImGui::PushID(i);
        if (ImGui::CollapsingHeader(name.cstr())) {
            if (ImGui::Checkbox("active", &ref.active)) {
                if (ref.active) {
                    for (auto& other : data->reference_frame.frames) {
                        if (&other == &ref) continue;
                        other.active = false;
                    }
                }
            }
            ImGui::SameLine();
            ImGui::Checkbox("show", &ref.show);
            ImGui::SameLine();
            if (ImGui::DeleteButton("remove")) {
                remove_reference_frame(data, &ref);
            }
            ImGui::SameLine();
            /*
            if (ImGui::Button("shape space")) {
                data->shape_space.show_window = true;
                data->shape_space.structure_list = {ref.id};
                draw_shape_space_window(data);
            }
            */

            ImGui::PushItemWidth(item_width);
            ImGui::InputText("name", ref.name.cstr(), ref.name.capacity());
            if (!ref.filter_is_ok) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
            if (ImGui::InputText("filter", ref.filter.cstr(), ref.filter.capacity(), ImGuiInputTextFlags_EnterReturnsTrue)) {
                recompute_frame = true;
            }
            // if (ImGui::SliderInt("reference frame idx", ))
            if (!ref.filter_is_ok) ImGui::PopStyleColor();
            ImGui::Combo("tracking mode", (int*)(&ref.tracking_mode), "Absolute\0Relative\0\0");

            if (ImGui::Button("Transform trajectory to internal reference frame") && ref.filter_is_ok) {
                Bitfield mask;
                bitfield::init(&mask, data->dynamic.molecule.atom.count);
                defer { bitfield::free(&mask); };
                filter::compute_filter_mask(mask, ref.filter, data->dynamic.molecule);
                structure_tracking::transform_to_internal_frame(data->dynamic, mask, 0);
            }

            ImGui::PopItemWidth();

            ImGui::PushItemWidth(-1);
            constexpr auto plot_height = 150.0f;
#if 0
            {
                const auto ev_data = structure_tracking::get_eigen_values(ref.id);
                float* tmp_data = (float*)TMP_MALLOC(ev_data.size_in_bytes());
                defer { TMP_FREE(tmp_data); };
                float* x = tmp_data + 0 * ev_data.size();
                float* y = tmp_data + 1 * ev_data.size();
                float* z = tmp_data + 2 * ev_data.size();
                for (int64 j = 0; j < ev_data.size(); j++) {
                    x[j] = ev_data[j].x;
                    y[j] = ev_data[j].y;
                    z[j] = ev_data[j].z;
                }

                const ImVec2 x_range = {0.0f, (float)ev_data.size()};
                const ImVec2 y_range = {-0.05f, 1.05f};
                ImGui::BeginPlot("Eigen Values", ImVec2(0, plot_height), x_range, y_range, ImGui::LinePlotFlags_AxisX | ImGui::LinePlotFlags_ShowXVal);
                ImGui::PlotValues("0", x, (int)ev_data.size(), 0xFF5555FF);
                ImGui::PlotValues("1", y, (int)ev_data.size(), 0xFF55FF55);
                ImGui::PlotValues("2", z, (int)ev_data.size(), 0xFFFF5555);
                float x_val;
                if (ImGui::ClickingAtPlot(&x_val)) {
                    data->playback.time = x_val;
                }
                ImGui::EndPlot();
            }
#endif
            const auto tracking_data = structure_tracking::get_tracking_data(ref.id);
            if (tracking_data && tracking_data->count > 0) {
                const i32 N = (i32)tracking_data->count;
                const quat* abs_data = tracking_data->transform.rotation.absolute;
                const quat* rel_data = tracking_data->transform.rotation.relative;
                // const quat* hyb_data = tracking_data->transform.rotation.hybrid;

                // Scratch data
                float* tmp_data = (float*)TMP_MALLOC(sizeof(float) * N * 3);
                defer { TMP_FREE(tmp_data); };

                float* abs_angle = tmp_data + 0 * N;
                float* rel_angle = tmp_data + 1 * N;
                float* del_angle = tmp_data + 2 * N;

                auto export_to_csv = [abs_angle, rel_angle, del_angle, N]() {
                    auto res = platform::file_dialog(platform::FileDialogFlags_Save, {}, "csv");
                    if (res.result == platform::FileDialogResult::Ok) {
                        FILE* file = fopen(res.path, "w");
                        if (file) {
                            fprintf(file, "time, ABS, REL, DEL\n");
                            for (int i = 0; i < N; i++) {
                                fprintf(file, "%i, %.4f, %.4f, %.4f\n", i, abs_angle[i], rel_angle[i], del_angle[i]);
                            }
                            fclose(file);
                        }
                    }
                };

                {
                    abs_angle[0] = 0.0f;
                    rel_angle[0] = 0.0f;
                    del_angle[0] = 0.0f;
                    const quat ref_q = {1, 0, 0, 0};
                    for (int j = 0; j < N; j++) {
                        abs_angle[j] = math::rad_to_deg(math::geodesic_distance(ref_q, abs_data[j]));
                        rel_angle[j] = math::rad_to_deg(math::geodesic_distance(ref_q, rel_data[j]));
                        del_angle[j] = math::rad_to_deg(math::geodesic_distance(abs_data[j], rel_data[j]));
                    }

                    if (ImGui::Button("Export to CSV##1")) export_to_csv();

                    const ImVec2 x_range = {0.0f, (float)N};
                    const ImVec2 y_range = {-0.05f, 180.5f};
                    ImGui::BeginPlot("Distance To Ref", ImVec2(0, plot_height), x_range, y_range,
                                     ImGui::LinePlotFlags_AxisX | ImGui::LinePlotFlags_ShowXVal);
                    ImGui::PlotValues("absolute", abs_angle, N, 0xFF5555FF);
                    ImGui::PlotValues("relative", rel_angle, N, 0xFF55FF55);
                    ImGui::PlotValues("delta", del_angle, N, 0xFFFF5555);
                    float x_val;
                    if (ImGui::ClickingAtPlot(&x_val)) {
                        data->playback.time = x_val;
                    }
                    ImGui::EndPlot();
                }
                {
                    abs_angle[0] = 0.0f;
                    rel_angle[0] = 0.0f;
                    for (int j = 1; j < N; j++) {
                        abs_angle[j] = math::rad_to_deg(math::angle(abs_data[j - 1], abs_data[j]));
                        rel_angle[j] = math::rad_to_deg(math::angle(rel_data[j - 1], rel_data[j]));
                    }

                    if (ImGui::Button("Export to CSV##2")) export_to_csv();

                    const ImVec2 x_range = {0.0f, (float)N};
                    const ImVec2 y_range = {-0.05f, 180.5f};
                    ImGui::BeginPlot("Preceding Angle Delta", ImVec2(0, plot_height), x_range, y_range,
                                     ImGui::LinePlotFlags_AxisX | ImGui::LinePlotFlags_ShowXVal);
                    ImGui::PlotValues("absolute", abs_angle, N, 0xFF5555FF);
                    ImGui::PlotValues("relative", rel_angle, N, 0xFF55FF55);
                    float x_val;
                    if (ImGui::ClickingAtPlot(&x_val)) {
                        data->playback.time = x_val;
                    }
                    ImGui::EndPlot();
                }
                {
                    abs_angle[0] = 0.0f;
                    rel_angle[0] = 0.0f;
                    for (int j = 1; j < N; j++) {
                        abs_angle[j] = abs_angle[j - 1] + math::rad_to_deg(math::angle(abs_data[j - 1], abs_data[j]));
                        rel_angle[j] = rel_angle[j - 1] + math::rad_to_deg(math::angle(rel_data[j - 1], rel_data[j]));
                    }

                    if (ImGui::Button("Export to CSV##3")) export_to_csv();

                    const float max_val = math::max(abs_angle[N - 1], rel_angle[N - 1]);
                    const ImVec2 x_range = {0.0f, (float)N};
                    const ImVec2 y_range = {0.0f, max_val};
                    ImGui::BeginPlot("Accumulated Angle Delta", ImVec2(0, plot_height), x_range, y_range,
                                     ImGui::LinePlotFlags_AxisX | ImGui::LinePlotFlags_ShowXVal);
                    ImGui::PlotValues("absolute", abs_angle, N, 0xFF5555FF);
                    ImGui::PlotValues("relative", rel_angle, N, 0xFF55FF55);
                    float x_val;
                    if (ImGui::ClickingAtPlot(&x_val)) {
                        data->playback.time = x_val;
                    }
                    ImGui::EndPlot();
                }
            }
            ImGui::PopItemWidth();

            ImGui::Spacing();
            ImGui::Separator();
        }

        ImGui::PopID();

        if (recompute_frame) {
            compute_reference_frame(data, &ref);
        }
    }

    ImGui::End();
}

static void draw_property_window(ApplicationData* data) {
    static bool first_time_shown = true;
    ImGui::Begin("Properties", &data->statistics.show_property_window, ImGuiWindowFlags_NoFocusOnAppearing);

    ImGui::PushID("PROPERTIES");
    ImGui::PushItemWidth(-1);
    // ImGui::Columns(4, "columns", true);
    const ImGuiColorEditFlags flags = 0;  // ImGuiColumnsFlags_NoPreserveWidths;
    ImGui::BeginColumns("columns", 4, flags);
    ImGui::Separator();

    if (first_time_shown) {
        first_time_shown = false;
        ImGui::SetColumnWidth(0, ImGui::GetWindowContentRegionWidth() * 0.15f);
        ImGui::SetColumnWidth(1, ImGui::GetWindowContentRegionWidth() * 0.65f);
        ImGui::SetColumnWidth(2, ImGui::GetWindowContentRegionWidth() * 0.1f);
        ImGui::SetColumnWidth(3, ImGui::GetWindowContentRegionWidth() * 0.1f);
    }

    ImGui::Text("name");
    ImGui::NextColumn();
    ImGui::Text("args");
    ImGui::NextColumn();
    ImGui::Text("S/T/D/V");
    ImGui::NextColumn();
    ImGui::NextColumn();

    // bool compute_stats = false;

    auto properties = stats::get_properties();
    for (i32 i = 0; i < (i32)properties.count; i++) {
        auto prop = properties[i];

        ImGui::Separator();
        ImGui::PushID(i);

        ImGui::PushItemWidth(-1);
        if (!prop->valid) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
        if (ImGui::InputText("##name", prop->name_buf.cstr(), prop->name_buf.capacity(), ImGuiInputTextFlags_EnterReturnsTrue)) {
            prop->data_dirty = true;
            // compute_stats = true;
        }
        if (!prop->valid) ImGui::PopStyleColor();
        ImGui::PopItemWidth();
        ImGui::NextColumn();

        constexpr auto buf_len = 128;
        char key_buf[buf_len] = {0};
        bool paste_buf = false;

        if (ImGui::BeginPopup("AtomContextMenu")) {
            if (data->dynamic.molecule) {
                const i32 atom_range = data->selection.right_clicked;
                if (atom_range != -1) {
                    ASSERT(atom_range < data->dynamic.molecule.atom.count);
                    const i32 residue_idx = data->dynamic.molecule.atom.res_idx[data->selection.right_clicked];
                    const i32 chain_idx = data->dynamic.molecule.atom.chain_idx[data->selection.right_clicked];

                    char buf[buf_len];

                    snprintf(buf, buf_len, "atom(%i) ", data->selection.right_clicked + 1);
                    if (ImGui::MenuItem(buf)) {
                        memcpy(key_buf, buf, buf_len);
                        paste_buf = true;
                    }

                    if (residue_idx > -1) {
                        const auto& residue = data->dynamic.molecule.residues[residue_idx];

                        snprintf(buf, buf_len, "residue(%i) ", residue_idx + 1);
                        if (ImGui::MenuItem(buf)) {
                            memcpy(key_buf, buf, buf_len);
                            paste_buf = true;
                        }
                        snprintf(buf, buf_len, "resid(%i) ", residue.id);
                        if (ImGui::MenuItem(buf)) {
                            memcpy(key_buf, buf, buf_len);
                            paste_buf = true;
                        }
                        snprintf(buf, buf_len, "resname(%s) ", residue.name.cstr());
                        if (ImGui::MenuItem(buf)) {
                            memcpy(key_buf, buf, buf_len);
                            paste_buf = true;
                        }

                        if (ImGui::BeginMenu("resatom...")) {
                            snprintf(buf, buf_len, "resatom(resid(%i), %i) ", residue.id, atom_range + 1);
                            if (ImGui::MenuItem(buf)) {
                                memcpy(key_buf, buf, buf_len);
                                paste_buf = true;
                            }
                            snprintf(buf, buf_len, "resatom(resname(%s), %i) ", residue.name.cstr(), atom_range + 1);
                            if (ImGui::MenuItem(buf)) {
                                memcpy(key_buf, buf, buf_len);
                                paste_buf = true;
                            }
                            ImGui::EndMenu();
                        }
                    }
                    if (chain_idx > -1) {
                        snprintf(buf, buf_len, "chain(%i) ", chain_idx + 1);
                        if (ImGui::MenuItem(buf)) {
                            memcpy(key_buf, buf, buf_len);
                            paste_buf = true;
                        }
                    }
                }
            }
            ImGui::EndPopup();
        }

        if (paste_buf) {
            ImGui::SetActiveID(ImGui::GetID("##args"), ImGui::GetCurrentWindow());
            ImGui::SetKeyboardFocusHere();
            ImGui::GetIO().AddInputCharactersUTF8(key_buf);
        }

        ImGui::PushItemWidth(-1);
        if (!prop->valid) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
        if (ImGui::InputText("##args", prop->args_buf.cstr(), prop->args_buf.capacity(), ImGuiInputTextFlags_EnterReturnsTrue)) {
            prop->data_dirty = true;
        }
        if (!prop->valid) ImGui::PopStyleColor();
        ImGui::PopItemWidth();

        if (ImGui::IsItemActive() && data->selection.hovered != -1 && data->ctx.input.mouse.release[1]) {
            ImGui::OpenPopup("AtomContextMenu");
        }

        if (!prop->valid) {
            if (!prop->valid && prop->error_msg_buf && ImGui::GetHoveredID() == ImGui::GetID("##args")) {
                ImGui::SetTooltip("%s", prop->error_msg_buf.cstr());
            }
        }
        ImGui::NextColumn();
        ImGui::Checkbox("##visualize", &prop->enable_visualization);
        ImGui::SameLine();
        ImGui::Checkbox("##timeline", &prop->enable_timeline);
        ImGui::SameLine();
        ImGui::Checkbox("##distribution", &prop->enable_distribution);
        ImGui::SameLine();
        if (ImGui::Checkbox("##volume", &prop->enable_volume)) {
            // Trigger update of volume
            prop->filter_dirty = true;
        }

        ImGui::NextColumn();
        if (ImGui::ArrowButton("up", ImGuiDir_Up)) {
            stats::move_property_up(prop);
        }
        ImGui::SameLine();
        if (ImGui::ArrowButton("down", ImGuiDir_Down)) {
            stats::move_property_down(prop);
        }
        ImGui::SameLine();
        if (ImGui::DeleteButton("remove")) {
            stats::remove_property(prop);
        }
        ImGui::NextColumn();
        ImGui::PopID();
    }
    ImGui::EndColumns();
    // ImGui::Columns(1);
    ImGui::Separator();
    ImGui::PopID();
    ImGui::PopItemWidth();

    if (ImGui::Button("create new")) {
        stats::create_property();
    }
    ImGui::SameLine();
    if (ImGui::DeleteButton("remove all")) {
        stats::remove_all_properties();
    }
    ImGui::End();

    // if (compute_stats) {
    // stats::compute_stats(data->dynamic);
    // compute_statistics_async(data);
    //}
}

static void draw_atom_info_window(const MoleculeStructure& mol, int atom_idx, int x, int y) {

    // @TODO: Assert things and make this failproof
    if (atom_idx < 0 || atom_idx >= mol.atom.count) return;

    int res_idx = mol.atom.res_idx[atom_idx];
    const Residue& res = mol.residues[res_idx];
    const char* res_name = res.name.cstr();
    const int res_id = res.id;
    int local_idx = atom_idx - res.atom_range.beg;
    const float pos_x = mol.atom.position.x[atom_idx];
    const float pos_y = mol.atom.position.y[atom_idx];
    const float pos_z = mol.atom.position.z[atom_idx];
    const char* label = mol.atom.label[atom_idx].cstr();
    const char* elem = element::name(mol.atom.element[atom_idx]);
    const char* symbol = element::symbol(mol.atom.element[atom_idx]);

    int chain_idx = mol.atom.chain_idx[atom_idx];
    const char* chain_id = "\0";
    if (chain_idx != -1 && mol.chains.size() > 0) {
        const Chain& chain = mol.chains[chain_idx];
        chain_id = chain.id.cstr();
    }

    int seq_idx = mol.atom.seq_idx[atom_idx];

    // External indices begin with 1 not 0
    res_idx += 1;
    chain_idx += 1;
    seq_idx += 1;
    atom_idx += 1;
    local_idx += 1;

    char buff[256];
    int len = 0;
    len += snprintf(buff, 256, "atom[%i][%i]: %s %s %s (%.2f, %.2f, %.2f)\n", atom_idx, local_idx, label, elem, symbol, pos_x, pos_y, pos_z);
    len += snprintf(buff + len, 256 - len, "res[%i]: %s %i\n", res_idx, res_name, res_id);
    if (chain_idx) {
        len += snprintf(buff + len, 256 - len, "chain[%i]: %s\n", chain_idx, chain_id);
    }
    if (seq_idx) {
        len += snprintf(buff + len, 256 - len, "seq[%i]\n", seq_idx);
    }

    if (res_idx < mol.backbone.angles.size() && res_idx < mol.backbone.segments.size() && valid_segment(mol.backbone.segments[res_idx])) {
        const auto angles = math::rad_to_deg((vec2)mol.backbone.angles[res_idx]);
        len += snprintf(buff + len, 256 - len, u8"\u03C6: %.1f\u00b0, \u03C8: %.1f\u00b0\n", angles.x, angles.y);
    }

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
    if (!data->dynamic.molecule) return;

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
        if (const auto& mol = data->dynamic.molecule) {
            ImGui::Text("MOL");
            const auto file = get_file(data->files.molecule);
            if (file) ImGui::Text("\"%.*s\"", file.size(), file.data());
            ImGui::Text("# atoms: %lli", mol.atom.count);
            ImGui::Text("# residues: %lli", mol.residues.count);
            ImGui::Text("# chains: %lli", mol.chains.count);
            ImGui::Text("# sequences: %lli", mol.sequences.count);
        }
        if (const auto& traj = data->dynamic.trajectory) {
            ImGui::NewLine();
            ImGui::Text("TRAJ");
            const auto file = get_file(data->files.trajectory);
            ImGui::Text("\"%.*s\"", file.size(), file.data());
            ImGui::Text("# frames: %lli", traj.num_frames);
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

    const f32 stats_fract = stats::fraction_done();

    const u32 num_tasks = task_system::get_num_tasks();
    task_system::ID* tasks = task_system::get_tasks();

    if ((0.f < stats_fract && stats_fract < 1.f || num_tasks > 0)) {

        ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(viewport->Pos + ImVec2(data->ctx.window.width - WIDTH - MARGIN,
                                                       ImGui::GetCurrentContext()->FontBaseSize + ImGui::GetStyle().FramePadding.y * 2.f + MARGIN));
        ImGui::SetNextWindowSize(ImVec2(WIDTH, 0));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
        ImGui::Begin("##Async Info", 0,
                     ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar |
                         ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoFocusOnAppearing);

        char buf[32];
        if (0.f < stats_fract && stats_fract < 1.f) {
            snprintf(buf, 32, "%.1f%%", stats_fract * 100.f);
            ImGui::ProgressBar(stats_fract, ImVec2(ImGui::GetWindowContentRegionWidth() * PROGRESSBAR_WIDTH_FRACT, 0), buf);
            ImGui::SameLine();
            ImGui::Text("Computing Statistics");
        }

        for (u32 i = 0; i < num_tasks; i++) {
            const auto id = tasks[i];
            const f32 fract = task_system::get_task_fraction_complete(id);
            snprintf(buf, 32, "%.1f%%", fract * 100.f);
            ImGui::ProgressBar(fract, ImVec2(ImGui::GetWindowContentRegionWidth() * PROGRESSBAR_WIDTH_FRACT, 0), buf);
            ImGui::SameLine();
            ImGui::Text(task_system::get_task_label(id));
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
    ImGui::SetNextWindowSize(ImVec2(400, 150), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Timelines", &data->statistics.show_timeline_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
        static float zoom = 1.f;
        const float num_frames = (float)data->dynamic.trajectory.num_frames;

        ImGui::Checkbox("Dynamic Framewindow", &data->time_filter.dynamic_window);
        if (data->time_filter.dynamic_window) {
            ImGui::SameLine();
            ImGui::SliderFloat("Window Extent", &data->time_filter.window_extent, 1.f, num_frames);
        }
        const float width = ImGui::GetContentRegionAvail().x - ImGui::GetStyle().WindowPadding.x * 2.0f;
        ImGui::SetNextWindowContentWidth(width * zoom);
        ImGui::BeginChild("Scroll Region", ImVec2(0.0f, 0.0f), true, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_HorizontalScrollbar);

        const Range<float> frame_range = {0, num_frames - 1};
        const Range<float> old_range = data->time_filter.range;

        const float cont_min_x = ImGui::GetWindowContentRegionMin().x + ImGui::GetWindowPos().x;

        ImGui::PushItemWidth(-1.0f);
        if (ImGui::RangeSliderFloat("###selection_range", &data->time_filter.range.beg, &data->time_filter.range.end, frame_range.beg,
                                    frame_range.end)) {
            if (data->time_filter.dynamic_window) {
                if (data->time_filter.range != old_range) {
                    data->playback.time = math::lerp(data->time_filter.range.beg, data->time_filter.range.end, 0.5f);
                } else {
                    if (data->time_filter.range.beg != old_range.beg) {
                        data->time_filter.window_extent = 2.f * math::abs((float)data->playback.time - data->time_filter.range.beg);
                    } else if (data->time_filter.range.end != old_range.end) {
                        data->time_filter.window_extent = 2.f * math::abs((float)data->playback.time - data->time_filter.range.end);
                    }
                }
            }
        }

        // const int32 prop_count = stats::get_property_count();
        // const float32 plot_height = ImGui::GetContentRegionAvail().y / (float)prop_count;
        const f32 plot_height = 100.f;
        const u32 bar_fill_color = ImColor(1.f, 1.f, 1.f, 0.25f);
        const u32 var_fill_color = ImColor(1.f, 1.f, 0.3f, 0.1f);
        const u32 var_line_color = ImColor(1.f, 1.f, 0.3f, 0.3f);
        const u32 var_text_color = ImColor(1.f, 1.f, 0.3f, 0.5f);

        static f32 selection_start;
        static bool is_selecting = false;

        auto style = ImGui::GetStyle();
        ImGui::GetStyleColorVec4(ImGuiCol_PlotLines);

        const auto properties = stats::get_properties();
        for (int i = 0; i < (i32)properties.count; i++) {
            auto prop = properties[i];

            if (!prop->enable_timeline) continue;
            Array<float> prop_data = prop->avg_data;
            CStringView prop_name = prop->name_buf;
            auto prop_range = prop->avg_data_range;
            if (!prop_data) continue;
            const float pad = math::abs(prop_range.end - prop_range.beg) * 0.75f;
            Range<float> display_range = {prop_range.beg - pad, prop_range.end + pad};
            if (display_range.beg == display_range.end) {
                display_range.beg -= 1.f;
                display_range.end += 1.f;
            }

            const ImGuiID id = ImGui::GetID(prop_name.cstr());
            ImGui::PushID(i);

            ImGui::BeginPlot(prop_name.cstr(), ImVec2(0, plot_height), ImVec2(frame_range.beg, frame_range.end),
                             ImVec2(display_range.beg, display_range.end), ImGui::LinePlotFlags_AxisX);
            const ImRect inner_bb(ImGui::GetItemRectMin() + ImGui::GetStyle().FramePadding, ImGui::GetItemRectMax() - ImGui::GetStyle().FramePadding);
            ImGui::PushClipRect(ImGui::GetItemRectMin(), ImGui::GetItemRectMax(), true);

            if (ImGui::IsItemHovered()) ImGui::SetHoveredID(id);

            ImGui::PlotVerticalBars(prop->filter_fraction.data(), (i32)prop->filter_fraction.size(), bar_fill_color);
            if (prop->std_dev_data.data()[0] > 0.f) {
                ImGui::PlotVariance(prop->avg_data.data(), prop->std_dev_data.data(), (i32)prop->std_dev_data.size(), 1.f, var_line_color,
                                    var_fill_color);
            }
            ImGui::PlotValues(prop->name_buf.cstr(), prop_data.data(), (i32)prop_data.size(),
                              ImGui::ColorConvertFloat4ToU32(style.Colors[ImGuiCol_PlotLines]));

            ImGui::PopClipRect();

            if (ImGui::IsItemHovered() && ImGui::GetIO().MouseClicked[0]) {
                ImGui::SetActiveID(id, ImGui::GetCurrentWindow());
            }

            if (ImGui::IsItemHovered() && ImGui::GetIO().MouseClicked[1] && ImGui::GetIO().KeyCtrl) {
                data->time_filter.range = frame_range;
            }

            if (ImGui::GetActiveID() == id) {
                const float t = ImClamp((ImGui::GetIO().MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x), 0.0f, 1.0f);
                if (ImGui::GetIO().MouseClicked[0] && ImGui::GetIO().KeyCtrl) {
                    selection_start = ImLerp(frame_range.beg, frame_range.end, t);
                    data->time_filter.range.beg = selection_start;
                    data->time_filter.range.end = selection_start;
                    is_selecting = true;
                } else if (is_selecting) {
                    const float v = ImLerp(frame_range.beg, frame_range.end, t);
                    if (v < data->time_filter.range.beg) {
                        data->time_filter.range.beg = v;
                    } else if (v > data->time_filter.range.beg && v < data->time_filter.range.end) {
                        if (selection_start < v) {
                            data->time_filter.range.end = v;
                        } else {
                            data->time_filter.range.beg = v;
                        }
                    } else if (v > data->time_filter.range.end) {
                        data->time_filter.range.end = v;
                    }
                } else if (ImGui::GetIO().MouseDown[0]) {
                    data->playback.time = ImLerp(frame_range.beg, frame_range.end, t);
                }

                if (!ImGui::GetIO().MouseDown[0] && !ImGui::IsItemHovered()) {
                    ImGui::ClearActiveID();
                    is_selecting = false;
                }
            }

            data->time_filter.range.beg = ImClamp(data->time_filter.range.beg, frame_range.beg, frame_range.end);
            data->time_filter.range.end = ImClamp(data->time_filter.range.end, frame_range.beg, frame_range.end);

            // SELECTION RANGE
            {
                constexpr ImU32 SELECTION_RANGE_COLOR = 0x55bbbbbb;
                const f32 t0 = (data->time_filter.range.beg - frame_range.beg) / frame_range.ext();
                const f32 t1 = (data->time_filter.range.end - frame_range.beg) / frame_range.ext();
                const ImVec2 pos0 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t0, 0));
                const ImVec2 pos1 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t1, 1));
                ImGui::GetCurrentWindow()->DrawList->AddRectFilled(pos0, pos1, SELECTION_RANGE_COLOR);
            }

            // CURRENT FRAME POSITION
            {
                constexpr ImU32 CURRENT_LINE_COLOR = 0xaa33ffff;
                const f32 t = ((float)data->playback.time - frame_range.beg) / frame_range.ext();
                const ImVec2 pos0 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t, 0));
                const ImVec2 pos1 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t, 1));
                ImGui::GetCurrentWindow()->DrawList->AddLine(pos0, pos1, CURRENT_LINE_COLOR);
            }

            // HOVERED CURSOR POSITION
            if (ImGui::IsItemHovered()) {
                constexpr ImU32 HOVER_LINE_COLOR = 0xaaffffff;
                const ImVec2 pos0(roundf(ImGui::GetIO().MousePos.x), inner_bb.Min.y);
                const ImVec2 pos1(roundf(ImGui::GetIO().MousePos.x), inner_bb.Max.y);
                ImGui::GetCurrentWindow()->DrawList->AddLine(pos0, pos1, HOVER_LINE_COLOR);
            }

            // TOOLTIP
            if (ImGui::GetActiveID() == id || ImGui::GetHoveredID() == id) {
                const f32 min_x = ImGui::GetItemRectMin().x;
                const f32 max_x = ImGui::GetItemRectMax().x;
                const f32 t = ImClamp((ImGui::GetIO().MousePos.x - min_x) / (max_x - min_x), 0.f, 1.f);
                i32 idx = ImClamp((i32)ImLerp(frame_range.beg, frame_range.end, t), 0, (i32)prop->avg_data.size() - 1);
                const f32 value = prop->std_dev_data[idx];

                ImGui::BeginTooltip();
                ImGui::Text("%i: %g ", idx, prop->avg_data[idx]);
                if (value > 0.f) {
                    ImGui::SameLine();
                    ImGui::TextColored(ImColor(var_text_color), "(%g)", value);
                }
                ImGui::EndTooltip();
            }

            ImGui::EndPlot();
            ImGui::PopID();
        }
        ImGui::PopItemWidth();

        if (data->time_filter.range != old_range) {
            stats::set_all_property_flags(false, true);
        }

        const float x = (ImGui::GetIO().MousePos.x - cont_min_x);
        if (ImGui::IsWindowHovered()) {
            const float mouse_wheel_delta = ImGui::GetIO().MouseWheel;
            if (ImGui::GetIO().KeyCtrl && mouse_wheel_delta != 0.f) {
                constexpr float ZOOM_SCL = 0.1f;
                const float old_zoom = zoom;
                const float new_zoom = math::clamp(zoom + zoom * ZOOM_SCL * mouse_wheel_delta, 1.f, 100.f);
                const float delta_x = (new_zoom / old_zoom) * x - x;
                ImGui::SetScrollX(ImGui::GetScrollX() + delta_x);
                zoom = new_zoom;
            }

            if (ImGui::GetIO().KeyShift && ImGui::GetIO().MouseDown[0] && ImGui::GetIO().MouseDelta.x != 0.0f) {
                ImGui::SetScrollX(ImGui::GetScrollX() - ImGui::GetIO().MouseDelta.x);
            }
        }
        ImGui::EndChild();
    }

    ImGui::End();
}

static void draw_distribution_window(ApplicationData* data) {
    ImGui::SetNextWindowSize(ImVec2(200, 300), ImGuiCond_FirstUseEver);
    ImGui::Begin("Distributions", &data->statistics.show_distribution_window, ImGuiWindowFlags_NoFocusOnAppearing);
    ImGuiWindow* window = ImGui::GetCurrentWindow();
    const ImGuiStyle& style = ImGui::GetStyle();
    ImGui::PushItemWidth(-1);

    const auto properties = stats::get_properties();
    // constexpr float32 RANGE_SLIDER_HEIGHT = 26.f;
    // const float32 plot_height = ImGui::GetContentRegionAvail().y / (float)properties.count - RANGE_SLIDER_HEIGHT;
    const f32 plot_height = 100.f;
    const ImVec2 frame_size{ImGui::CalcItemWidth(), plot_height};

    constexpr u32 FULL_FILL_COLOR = 0x99cc9e66;
    constexpr u32 FULL_LINE_COLOR = 0xffcc9e66;
    constexpr u32 FULL_TEXT_COLOR = 0xffcc9e66;

    constexpr u32 FILT_FILL_COLOR = 0x3333ffff;
    constexpr u32 FILT_LINE_COLOR = 0xaa33ffff;
    constexpr u32 FILT_TEXT_COLOR = 0xaa33ffff;
    constexpr ImU32 SELECTION_RANGE_COLOR = 0x55bbbbbb;

    for (i32 i = 0; i < (i32)properties.count; i++) {
        stats::Property* prop = properties[i];
        if (!prop->enable_distribution) continue;
        ImGui::PushID(i);

        const ImRect frame_bb(window->DC.CursorPos, window->DC.CursorPos + ImVec2(frame_size.x, frame_size.y));
        const ImRect inner_bb(frame_bb.Min + style.FramePadding, frame_bb.Max - style.FramePadding);
        const ImRect total_bb(frame_bb.Min, frame_bb.Max);
        ImGui::ItemSize(total_bb, style.FramePadding.y);
        if (ImGui::ItemAdd(total_bb, 0)) {
            ImGui::RenderFrame(frame_bb.Min, frame_bb.Max, ImGui::GetColorU32(ImGuiCol_FrameBg), true, style.FrameRounding);

            // const float32 max_val = math::max(prop->full_histogram.bin_range.y, prop->filt_histogram.bin_range.y);
            ImGui::PushClipRect(inner_bb.Min, inner_bb.Max, true);
            const f32 max_val = prop->full_histogram.bin_range.end * 1.5f;
            ImGui::DrawFilledLine(inner_bb.Min, inner_bb.Max, prop->full_histogram.bins.ptr, (i32)prop->full_histogram.bins.count, max_val,
                                  FULL_LINE_COLOR, FULL_FILL_COLOR);

            ImGui::DrawFilledLine(inner_bb.Min, inner_bb.Max, prop->filt_histogram.bins.ptr, (i32)prop->filt_histogram.bins.count, max_val,
                                  FILT_LINE_COLOR, FILT_FILL_COLOR);
            // ImGui::PopClipRect();

            // SELECTION RANGE
            {
                const f32 t0 = (prop->filter.beg - prop->total_data_range.beg) / (prop->total_data_range.ext());
                const f32 t1 = (prop->filter.end - prop->total_data_range.beg) / (prop->total_data_range.ext());
                const ImVec2 pos0 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t0, 0));
                const ImVec2 pos1 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t1, 1));
                ImGui::GetCurrentWindow()->DrawList->AddRectFilled(pos0, pos1, SELECTION_RANGE_COLOR);
            }

            ImGui::PopClipRect();

            if (ImGui::IsItemHovered()) {
                window->DrawList->AddLine(ImVec2(ImGui::GetIO().MousePos.x, inner_bb.Min.y), ImVec2(ImGui::GetIO().MousePos.x, inner_bb.Max.y),
                                          0xffffffff);
                const f32 t = (ImGui::GetIO().MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x);
                const i32 count = (i32)prop->full_histogram.bins.count;
                const i32 idx = ImClamp((i32)(t * (count - 1)), 0, count - 1);
                const f32 full_val = prop->full_histogram.bins.ptr[idx];
                const f32 filt_val = prop->filt_histogram.bins.ptr[idx];
                ImGui::BeginTooltip();
                ImGui::Text("%.3f:", ImLerp(prop->filt_histogram.value_range.beg, prop->filt_histogram.value_range.end, t));
                ImGui::TextColored(ImColor(FULL_TEXT_COLOR), "%g", full_val * 100.f);
                ImGui::TextColored(ImColor(FILT_TEXT_COLOR), "%g", filt_val * 100.f);
                ImGui::EndTooltip();
            }

            if (ImGui::RangeSliderFloat("##filter", &prop->filter.beg, &prop->filter.end, prop->total_data_range.beg, prop->total_data_range.end)) {
                prop->filter_dirty = true;
            }
        }
        ImGui::PopID();
    }
    ImGui::PopItemWidth();
    ImGui::End();
}

static void draw_ramachandran_window(ApplicationData* data) {
    // const int32 num_frames = data->dynamic.trajectory ? data->dynamic.trajectory.num_frames : 0;
    // const int32 frame = (int32)data->time;
    const Range<i32> frame_range = data->time_filter.range;
    const auto& mol = data->dynamic.molecule;
    // Array<const BackboneAngle> trajectory_angles = get_backbone_angles(traj_angles, frame_range.beg, frame_range.ext());
    Array<const BackboneAngle> current_angles = mol.backbone.angles;
    Array<const BackboneSegment> backbone_segments = mol.backbone.segments;
    Array<const Residue> residues = mol.residues;
    Bitfield& atom_selection = data->selection.current_selection_mask;
    Bitfield& atom_highlight = data->selection.current_highlight_mask;

    ImGui::SetNextWindowSizeConstraints(ImVec2(400, 200), ImVec2(10000, 10000));
    ImGui::Begin("Ramachandran", &data->ramachandran.show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoScrollbar);

    ImGui::BeginColumns("cols", 2, ImGuiColumnsFlags_NoResize);
    {
        ImGui::Checkbox("Show Current Frame", &data->ramachandran.current.enabled);
        if (data->ramachandran.current.enabled) {
            ImGui::PushID("current");
            ImGui::SliderFloat("##base_radius", &data->ramachandran.current.base.radius, 0.5f, 5.f, "base radius %1.1f");
            // ImGui::SameLine();
            // ImGui::ColorEdit4("base color", (float*)&data->ramachandran.current.fill_color, ImGuiColorEditFlags_NoInputs |
            // ImGuiColorEditFlags_NoLabel);

            ImGui::SliderFloat("##selected_radius", &data->ramachandran.current.selection.radius, 0.5, 5.f, "selected radius %1.1f");
            // ImGui::SameLine();
            // ImGui::ColorEdit4("selected color", (float*)&data->ramachandran.selected.fill_color, ImGuiColorEditFlags_NoInputs |
            // ImGuiColorEditFlags_NoLabel);

            data->ramachandran.current.base.radius = math::round(data->ramachandran.current.base.radius * 2.f) / 2.f;
            data->ramachandran.current.selection.radius = math::round(data->ramachandran.current.selection.radius * 2.f) / 2.f;
            // ImGui::ColorEdit4("border color", (float*)&data->ramachandran.current.border_color, ImGuiColorEditFlags_NoInputs |
            // ImGuiColorEditFlags_NoLabel);
            ImGui::PopID();
        }
    }
    ImGui::NextColumn();
    {
        ImGui::Checkbox("Show Range", &data->ramachandran.range.enabled);
        if (data->ramachandran.range.enabled) {
            ImGui::PushID("range");
            ImGui::SliderFloat("", &data->ramachandran.range.radius, 0.1f, 5.f, "radius %1.1f");
            ImGui::SameLine();
            ImGui::ColorEdit4("color", (float*)&data->ramachandran.range.color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Fill color for trajectory range");
            }
            // ImGui::RangeSliderFloat("range", &data->time_filter.range.min, &data->time_filter.range.max, 0.f,
            //                        (float)data->dynamic.trajectory.num_frames, "(%.1f, %.1f)");
            ImGui::PopID();

            ramachandran::clear_accumulation_texture();
            ramachandran::render_accumulation_texture(frame_range, data->ramachandran.range.color, data->ramachandran.range.radius);
        }
    }
    ImGui::EndColumns();

    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    ImGui::BeginChild("canvas", ImVec2(-1, -1), true);
    ImGui::PopStyleVar(1);

    static float zoom_factor = 1.0f;
    const float max_c = ImMax(ImGui::GetContentRegionAvail().x, ImGui::GetContentRegionAvail().y);
    const ImVec2 size = ImVec2(max_c, max_c) * zoom_factor;
    // const ImGuiID id = ImGui::GetID("bg");

    ImRect bb(ImGui::GetCurrentWindow()->DC.CursorPos, ImGui::GetCurrentWindow()->DC.CursorPos + size);
    // ImGui::ItemSize(bb);
    // ImGui::ItemAdd(bb, id);

    // bool hovered = false, held = false;
    // bool pressed = ImGui::ButtonBehavior(bb, id, &hovered, &held);
    ImGui::InvisibleButton("bg", bb.Max - bb.Min);

    ImDrawList* dl = ImGui::GetWindowDrawList();

    const ImVec2 mouse_pos = ImGui::GetIO().MousePos;
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

        for (i64 i = 0; i < backbone_segments.size(); i++) {
            const auto& angle = current_angles[i];
            // const auto& seg = backbone_segments[i];
            const auto& res = residues[i];
            if (angle.phi == 0.f || angle.psi == 0.f) continue;

            const auto selected = bitfield::any_bit_set_in_range(atom_selection, res.atom_range);
            const auto highlight = bitfield::any_bit_set_in_range(atom_highlight, res.atom_range);
            const auto radius = (highlight || selected) ? selected_radius : base_radius;
            const auto fill_color = highlight ? highlight_color : (selected ? selected_color : base_color);

            const ImVec2 coord =
                ImLerp(bb.Min, bb.Max, ImVec2(angle.phi * ONE_OVER_TWO_PI + 0.5f, -angle.psi * ONE_OVER_TWO_PI + 0.5f));  // [-PI, PI] -> [0, 1]
            const ImVec2 min_box(math::round(coord.x - radius), math::round(coord.y - radius));
            const ImVec2 max_box(math::round(coord.x + radius), math::round(coord.y + radius));
            if (radius > 1.f) {
                dl->AddRectFilled(min_box, max_box, fill_color);
                dl->AddRect(min_box, max_box, border_color, 0.f, 15, 2.f);  // @TODO: REVERT THIS LATER
            } else {
                dl->AddRectFilled(min_box, max_box, border_color);
            }
            if (min_box.x <= mouse_pos.x && mouse_pos.x <= max_box.x && min_box.y <= mouse_pos.y && mouse_pos.y <= max_box.y) {
                mouse_hover_idx = i;
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

    if (ImGui::IsItemHovered()) {

        if (ImGui::GetIO().KeyCtrl && ImGui::GetIO().MouseWheel != 0.f) {
            // const float old_zoom = zoom_factor;
            const float new_zoom = ImClamp(zoom_factor - ImGui::GetIO().MouseWheel * 0.1f, 1.0f, 10.0f);
            // ImGui::GetCurrentWindow()->ScrollTarget.x = ImGui::GetCurrentWindow()->Scroll.x * new_zoom / old_zoom;
            zoom_factor = new_zoom;
        }

        if (!region_select && ImGui::GetIO().KeyShift && (ImGui::GetIO().MouseClicked[0] || ImGui::GetIO().MouseClicked[1])) {
            region_select = true;
            region_mode = ImGui::GetIO().MouseClicked[0] ? Mode::Append : Mode::Remove;
            region_x0 = ImGui::GetMousePos();
        }

        const ImVec2 normalized_coord = ((ImGui::GetMousePos() - bb.Min) / (bb.Max - bb.Min) - ImVec2(0.5f, 0.5f)) * ImVec2(1, -1);
        const ImVec2 angles = normalized_coord * 2.f * 180.f;
        ImGui::BeginTooltip();
        ImGui::Text(u8"\u03C6: %.1f\u00b0, \u03C8: %.1f\u00b0", angles.x, angles.y);
        if (!region_select && mouse_hover_idx != -1) {
            const auto res_idx = mouse_hover_idx;
            const auto& res = get_residues(data->dynamic.molecule)[res_idx];
            ImGui::Text("Residue[%lli]: %s", res_idx, res.name.cstr());
            bitfield::clear_all(data->selection.current_highlight_mask);
            bitfield::set_range(data->selection.current_highlight_mask, res.atom_range);
            data->gpu_buffers.dirty.selection = true;
        }
        ImGui::EndTooltip();
    }

    if (region_select) {
        const ImVec2 region_x1 = ImGui::GetMousePos();
        const ImVec2 x0 = ImMin(region_x0, region_x1);
        const ImVec2 x1 = ImMax(region_x0, region_x1);
        const ImU32 fill_col = 0x22222222;
        const ImU32 line_col = 0x88888888;
        dl->AddRectFilled(x0, x1, fill_col);
        dl->AddRect(x0, x1, line_col);

        Bitfield mask;
        bitfield::init(&mask, mol.atom.count);
        defer { bitfield::free(&mask); };

        for (i64 i = 0; i < residues.size(); i++) {
            const auto& angle = current_angles[i];
            if (angle.phi == 0.f || angle.psi == 0.f) continue;
            const ImVec2 coord =
                ImLerp(bb.Min, bb.Max, ImVec2(angle.phi * ONE_OVER_TWO_PI + 0.5f, -angle.psi * ONE_OVER_TWO_PI + 0.5f));  // [-PI, PI] -> [0, 1]
            if (coord.x < x0.x || x1.x < coord.x) continue;
            if (coord.y < x0.y || x1.y < coord.y) continue;

            bitfield::set_bit(mask, backbone_segments[i].ca_idx);
        }

        switch (data->selection.level_mode) {
            case SelectionLevel::Atom:
                break;
            case SelectionLevel::Residue:
                expand_mask(mask, mol.residues);
                break;
            case SelectionLevel::Chain:
                expand_mask(mask, mol.sequences);
                break;
            default:
                ASSERT(false);
        }

        bitfield::clear_all(data->selection.current_highlight_mask);
        switch (region_mode) {
            case Mode::Append:
                bitfield::or_field(data->selection.current_highlight_mask, data->selection.current_selection_mask, mask);
                break;
            case Mode::Remove:
                bitfield::and_not_field(data->selection.current_highlight_mask, data->selection.current_selection_mask, mask);
                break;
            default:
                ASSERT(false);
        }

        if (!ImGui::GetIO().KeyShift || !(ImGui::GetIO().MouseDown[0] || ImGui::GetIO().MouseDown[1])) {
            switch (region_mode) {
                case Mode::Append:
                    bitfield::or_field(data->selection.current_selection_mask, data->selection.current_selection_mask, mask);
                    break;
                case Mode::Remove:
                    bitfield::and_not_field(data->selection.current_selection_mask, data->selection.current_selection_mask, mask);
                    break;
                default:
                    ASSERT(false);
            }
            region_select = false;
            data->gpu_buffers.dirty.selection = true;
            update_all_representations(data);
        }
    }
    ImGui::EndChild();
    ImGui::End();
}

static void draw_shape_space_window(ApplicationData* data) {
    static struct Style {
        vec4 window_bg = {1.0f, 1.0f, 1.0f, 0.8f};

        struct {
            vec4 fill_color = {0.0f, 0.0f, 0.0f, 0.03f};
            vec4 line_color = {0.0f, 0.0f, 0.0f, 0.3f};
            float line_thickness = 3.0f;
        } triangle;

        struct Case {
            vec3 line_color;
            float alpha;
            float line_thickness;
            float radius;
        };
        bool use_base_for_all = true;
        Case base = {{0.15f, 0.15f, 0.15f}, 1.0f, 0.15f, 4.5f};
        Case hover = {{1.f, 1.f, 1.f}, 3.0f, 1.0f, 7.5f};
        Case selected = {{1.f, 1.f, 0.f}, 3.0f, 0.9f, 6.5f};

        const int num_segments = 10;
    } style;

    const Range<i32> frame_range = data->time_filter.range;

    ImGui::PushStyleColor(ImGuiCol_WindowBg, math::convert_color(style.window_bg));
    // ImGui::SetNextWindowSizeConstraints(ImVec2(200, 200), ImVec2(10000, 10000));
    ImGui::Begin("Shape Space", &data->shape_space.show_window, ImGuiWindowFlags_NoFocusOnAppearing);
    ImGui::PopStyleColor();

    static bool show_style_window = false;
    if (ImGui::Button("Style")) show_style_window = !show_style_window;

    if (show_style_window) {
        ImGui::Begin("Shape Space Style", &show_style_window, ImGuiWindowFlags_NoFocusOnAppearing);
        ImGui::ColorEdit4("Window BG", &style.window_bg[0]);
        ImGui::ColorEdit4("Triangle fill color", &style.triangle.fill_color[0]);
        ImGui::ColorEdit4("Triangle line color", &style.triangle.line_color[0]);
        ImGui::SliderFloat("Triangle line thickness", &style.triangle.line_thickness, 0.0f, 5.0f);
        ImGui::Checkbox("Use base for all", &style.use_base_for_all);
        ImGui::Text("Base");
        ImGui::PushID(0);
        ImGui::ColorEdit3("line color", &style.base.line_color[0]);
        ImGui::SliderFloat("line thickness", &style.base.line_thickness, 0, 1);
        ImGui::SliderFloat("alpha", &style.base.alpha, 0, 1);
        ImGui::SliderFloat("radius", &style.base.radius, 0.1f, 10.f);
        ImGui::Separator();
        ImGui::PopID();
        if (!style.use_base_for_all) {
            ImGui::Text("Hover");
            ImGui::PushID(1);
            ImGui::ColorEdit3("line color", &style.hover.line_color[0]);
            ImGui::SliderFloat("line thickness", &style.hover.line_thickness, 0, 1);
            ImGui::SliderFloat("alpha", &style.hover.alpha, 0, 1);
            ImGui::SliderFloat("radius", &style.hover.radius, 0.1f, 10.f);
            ImGui::Separator();
            ImGui::PopID();
            ImGui::PushID(2);
            ImGui::Text("Selected");
            ImGui::ColorEdit3("line color", &style.selected.line_color[0]);
            ImGui::SliderFloat("line thickness", &style.selected.line_thickness, 0, 1);
            ImGui::SliderFloat("alpha", &style.selected.alpha, 0, 1);
            ImGui::SliderFloat("radius", &style.selected.radius, 0.1f, 10.f);
            ImGui::Separator();
            ImGui::PopID();
        }

        ImGui::End();
    }

    const i32 ensemble_count = data->ensemble_tracking.structures.empty() ? 0 : 1;
    const i32 reference_count = (i32)data->reference_frame.frames.size();
    const i32 count = ensemble_count + reference_count;

    structure_tracking::ID structure_ids[512];
    i32 num_structure_ids = 0;

    if (count > 0) {
        StringBuffer<512> buf = {};
        for (const auto ref_frame : data->reference_frame.frames) {
            buf += "Reference Frame: ";
            buf += ref_frame.name;
        }
        if (!data->ensemble_tracking.structures.empty()) {
            if (buf.length() == 0) {
                buf += "Ensemble\0";
            } else {
                buf += "\0Ensemble\0";
            }
        }

        ImGui::Combo("Target", (int*)(&data->shape_space.target), buf.cstr());

        int i = data->shape_space.target;
        if (i < reference_count) {
            structure_ids[0] = data->reference_frame.frames[i].id;
            num_structure_ids = 1;
        } else if (i < count) {
            for (const auto& s : data->ensemble_tracking.structures) {
                if (s.enabled) {
                    structure_ids[num_structure_ids] = s.id;
                    num_structure_ids++;
                }
            }
        }
    }

    constexpr float height_ratio = 0.86602540f;  // sqrt(3.0f) / 2.0f   base -> height ratio in equilateral triangle
    const ImVec2 size = ImVec2(ImGui::GetContentRegionAvail().x, ImGui::GetContentRegionAvail().x * height_ratio);
    const ImRect bb(ImGui::GetCurrentWindow()->DC.CursorPos, ImGui::GetCurrentWindow()->DC.CursorPos + size);

    const ImGuiID id = ImGui::GetCurrentWindow()->GetID("canvas");
    ImGui::ItemAdd(bb, id);

    // ImGui::RenderFrame(bb.Min, bb.Max, ImGui::GetColorU32(ImGuiCol_FrameBg), true, ImGui::GetStyle().FrameRounding);

    // ImGui::InvisibleButton("bg", bb.Max - bb.Min);

    const float cube_size_a = 30;
    const float cube_size_b = 10;
    const ImU32 cube_face_color_a = 0xCC111111;
    const ImU32 cube_face_color_b = 0xCC777777;
    const ImU32 cube_face_color_c = 0xCCCCCCCC;

    const ImVec2 tri_pad = {2 * cube_size_a, 2 * cube_size_a};
    const ImVec2 tri_a = ImLerp(bb.Min + tri_pad, bb.Max - tri_pad, ImVec2(0.0f, 1.0f));
    const ImVec2 tri_b = ImLerp(bb.Min + tri_pad, bb.Max - tri_pad, ImVec2(1.0f, 1.0f));
    const ImVec2 tri_c = ImLerp(bb.Min + tri_pad, bb.Max - tri_pad, ImVec2(0.5f, 0.0f));

    const ImVec2 cube_pad = tri_pad * 0.5f;
    const ImVec2 cube_a = ImLerp(bb.Min + cube_pad, bb.Max - cube_pad, ImVec2(0.0f, 1.0f));
    const ImVec2 cube_b = ImLerp(bb.Min + cube_pad, bb.Max - cube_pad, ImVec2(1.0f, 1.0f));
    const ImVec2 cube_c = ImLerp(bb.Min + cube_pad, bb.Max - cube_pad, ImVec2(0.5f, 0.0f));

    ImDrawList* dl = ImGui::GetWindowDrawList();

    const auto draw_ortho_cube = [dl](float x, float y, float w, float h, float d, u32 col_a, u32 col_b, u32 col_c) {
        const vec3 f_a[4] = {{0, 0, d}, {w, 0, d}, {w, h, d}, {0, h, d}};
        const vec3 f_b[4] = {{w, 0, d}, {w, 0, 0}, {w, h, 0}, {w, h, d}};
        const vec3 f_c[4] = {{0, h, d}, {w, h, d}, {w, h, 0}, {0, h, 0}};

        const mat4 T = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {-0.5f * w, -0.5f * h, -0.5f * d, 1.0f}};
        const mat4 view_mat = look_at(vec3(1.0f, 0.5f, 1.0f), vec3(0, 0, 0), vec3(0, 1, 0));
        const mat4 proj_mat = {{1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 1}};
        const mat4 mvp = proj_mat * view_mat * T;

        const vec4 v_a[4] = {mvp * vec4(f_a[0], 1.0f), mvp * vec4(f_a[1], 1.0f), mvp * vec4(f_a[2], 1.0f), mvp * vec4(f_a[3], 1.0f)};
        const vec4 v_b[4] = {mvp * vec4(f_b[0], 1.0f), mvp * vec4(f_b[1], 1.0f), mvp * vec4(f_b[2], 1.0f), mvp * vec4(f_b[3], 1.0f)};
        const vec4 v_c[4] = {mvp * vec4(f_c[0], 1.0f), mvp * vec4(f_c[1], 1.0f), mvp * vec4(f_c[2], 1.0f), mvp * vec4(f_c[3], 1.0f)};

        const ImVec2 o = {x, y};
        const ImVec2 p_a[4] = {o + ImVec2(v_a[0].x, v_a[0].y), o + ImVec2(v_a[1].x, v_a[1].y), o + ImVec2(v_a[2].x, v_a[2].y),
                               o + ImVec2(v_a[3].x, v_a[3].y)};
        const ImVec2 p_b[4] = {o + ImVec2(v_b[0].x, v_b[0].y), o + ImVec2(v_b[1].x, v_b[1].y), o + ImVec2(v_b[2].x, v_b[2].y),
                               o + ImVec2(v_b[3].x, v_b[3].y)};
        const ImVec2 p_c[4] = {o + ImVec2(v_c[0].x, v_c[0].y), o + ImVec2(v_c[1].x, v_c[1].y), o + ImVec2(v_c[2].x, v_c[2].y),
                               o + ImVec2(v_c[3].x, v_c[3].y)};

        dl->AddQuadFilled(p_a[0], p_a[1], p_a[2], p_a[3], col_a);
        dl->AddQuadFilled(p_b[0], p_b[1], p_b[2], p_b[3], col_b);
        dl->AddQuadFilled(p_c[0], p_c[1], p_c[2], p_c[3], col_c);
    };

    dl->ChannelsSplit(4);
    dl->ChannelsSetCurrent(0);
    dl->AddTriangleFilled(tri_a, tri_b, tri_c, math::convert_color(style.triangle.fill_color));
    dl->AddTriangle(tri_a, tri_b, tri_c, math::convert_color(style.triangle.line_color), style.triangle.line_thickness);

    draw_ortho_cube(cube_a.x, cube_a.y, cube_size_a * 1.5f, cube_size_b, cube_size_b, cube_face_color_a, cube_face_color_b, cube_face_color_c);
    draw_ortho_cube(cube_b.x, cube_b.y, cube_size_a * 1.2f, cube_size_b, cube_size_a * 1.2f, cube_face_color_a, cube_face_color_b, cube_face_color_c);
    draw_ortho_cube(cube_c.x, cube_c.y, cube_size_a, cube_size_a, cube_size_a, cube_face_color_a, cube_face_color_b, cube_face_color_c);

    float mouse_hover_d2 = FLT_MAX;
    int mouse_hover_idx = -1;
    vec3 mouse_hover_ev = {0, 0, 0};
    vec3 mouse_hover_w = {0, 0, 0};

    if (!data->shape_space.target) {
        for (int j = 0; j < num_structure_ids; j++) {
            const auto structure_id = structure_ids[j];
            if (const auto* tracking_data = structure_tracking::get_tracking_data(structure_id)) {
                const i32 N = (i32)tracking_data->count;
                const vec3* eigen_values = tracking_data->eigen.values;

                const auto compute_shape_space_weights = [](const vec3& ev) -> vec3 {
                    const float l_sum = ev[0] + ev[1] + ev[2];
                    const float one_over_denom = 1.0f / l_sum;
                    const float cl = (ev[0] - ev[1]) * one_over_denom;
                    const float cp = 2.0f * (ev[1] - ev[2]) * one_over_denom;
                    const float cs = 3.0f * ev[2] * one_over_denom;
                    return {cl, cp, cs};
                };

                const vec2 p[3] = {vec_cast(tri_a), vec_cast(tri_b), vec_cast(tri_c)};
                for (i32 i = 0; i < N; i++) {
                    const vec3 ev = eigen_values[i];
                    const vec3 w = compute_shape_space_weights(ev);
                    const ImVec2 coord = vec_cast(math::barycentric_to_cartesian(p[0], p[1], p[2], w));

                    const bool in_range = (frame_range.beg <= i && i <= frame_range.end);
                    const bool selected = (!style.use_base_for_all) && (i == data->playback.frame);

                    // Draw channel, higher number => later draw
                    int channel = 1 + (int)in_range + (int)selected;

                    const float t = (float)i / (float)(N - 1);
                    const vec3 color = green_color_scale(t);

                    Style::Case s = style.base;
                    if (!style.use_base_for_all && selected) s = style.selected;

                    dl->ChannelsSetCurrent(channel);
                    dl->AddCircleFilled(coord, s.radius, math::convert_color(vec4(color, s.alpha)), style.num_segments);
                    dl->AddCircle(coord, s.radius, math::convert_color(vec4(s.line_color, s.alpha)), style.num_segments, s.line_thickness);

                    if (ImGui::IsItemHovered()) {
                        const float d2 = ImLengthSqr(coord - ImGui::GetIO().MousePos);
                        if (d2 < style.hover.radius * style.hover.radius && d2 < mouse_hover_d2) {
                            mouse_hover_d2 = d2;
                            mouse_hover_idx = i;
                            mouse_hover_ev = ev;
                            mouse_hover_w = w;
                        }
                    }
                }

                if (mouse_hover_idx != -1) {
                    // Draw hovered explicity
                    const int i = mouse_hover_idx;
                    const vec3 ev = eigen_values[i];
                    const ImVec2 coord = vec_cast(math::barycentric_to_cartesian(p[0], p[1], p[2], compute_shape_space_weights(ev)));
                    const float t = (float)i / (float)(N - 1);
                    const vec3 color = green_color_scale(t);

                    dl->ChannelsSetCurrent(3);
                    dl->AddCircleFilled(coord, style.hover.radius, math::convert_color(vec4(color, style.hover.alpha)), style.num_segments);
                    dl->AddCircle(coord, style.hover.radius, math::convert_color(vec4(style.hover.line_color, style.hover.alpha)), style.num_segments,
                                  style.hover.line_thickness);
                }
            } else {
                LOG_ERROR("Could not find tracking data for '%lu'", structure_id);
            }
        }
    }

    dl->ChannelsMerge();
    dl->ChannelsSetCurrent(0);

    if (mouse_hover_idx != -1) {
        ImGui::BeginTooltip();
        ImGui::Text("Frame[%i]", (i32)mouse_hover_idx);
        if (ImGui::GetIO().MouseDown[0]) {
            data->playback.time = (float)mouse_hover_idx;
        }
        ImGui::Text(u8"\u03BB: (%.2f, %.2f, %.2f)", mouse_hover_ev.x, mouse_hover_ev.y, mouse_hover_ev.z);
        ImGui::Text(u8"\u03B1: %1.2f, \u03B2: %1.2f, \u03B3: %1.2f", mouse_hover_w.x, mouse_hover_w.y, mouse_hover_w.z);
        ImGui::EndTooltip();
    }

    ImGui::End();
}

static void append_trajectory_density(Volume* vol, Bitfield atom_mask, const MoleculeTrajectory& traj, const mat4& rotation,
                                      const structure_tracking::TrackingData& tracking_data, TrackingMode mode, float radial_cutoff) {
    ASSERT(vol);

    quat* rot_data = nullptr;
    switch (mode) {
        case TrackingMode::Absolute:
            rot_data = tracking_data.transform.rotation.absolute;
            break;
        case TrackingMode::Relative:
            rot_data = tracking_data.transform.rotation.relative;
            break;
        default:
            ASSERT(false);
    }

    const mat4 model_to_texture = volume::compute_model_to_texture_matrix(vol->dim);

    const float cutoff_squared = radial_cutoff * radial_cutoff;

    for (int frame_idx = 0; frame_idx < traj.num_frames; frame_idx++) {
        const auto& frame = get_trajectory_frame(traj, frame_idx);
        const mat3 box = frame.box;
        const vec3 half_box = box * vec3(0.5f);
        const vec3 min_aabb = box * vec3(0.0f);
        const vec3 max_aabb = box * vec3(1.0f);
        const vec3 com = tracking_data.transform.com[frame_idx];
        const mat4 R = math::mat4_cast(rot_data[frame_idx]);
        const mat4 T_com = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(-com, 1));
        const mat4 T_box = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(half_box, 1));
        const mat4 world_to_model = volume::compute_world_to_model_matrix(min_aabb, max_aabb);
        // const mat4 M = volume_world_to_texture * math::transpose(R) * T_com;
        const mat4 M = model_to_texture * world_to_model * T_box * rotation * R * T_com;

        // @TODO: Replace when proper iterators are implemented for bitfield.
        bitfield::for_each_bit_set(atom_mask, [&frame, &M, &com, cutoff_squared, &vol](i64 i) {
            const vec4 wc = {frame.atom_position.x[i], frame.atom_position.y[i], frame.atom_position.z[i], 1.0f};  // world coord
            if (math::distance2(com, vec3(wc)) > cutoff_squared) return;
            const vec4 vc = M * wc;  // volume coord [0,1]
            const vec3 tc = apply_pbc(vec3(vc));
            const ivec3 c = (ivec3)(tc * vec3(vol->dim - 1) + 0.5f);
            ASSERT(-1 < c.x && c.x < vol->dim.x);
            ASSERT(-1 < c.y && c.y < vol->dim.y);
            ASSERT(-1 < c.z && c.z < vol->dim.z);
            const i32 voxel_idx = c.z * vol->dim.x * vol->dim.y + c.y * vol->dim.x + c.x;
            atomic_fetch_and_add(vol->voxel_data.ptr + voxel_idx, 1);
        });
        /*
        for (i64 i = 0; i < atom_mask.size(); i++) {
            if (bitfield::get_bit(atom_mask, i)) {
                const vec4 wc = {frame.atom_position.x[i], frame.atom_position.y[i], frame.atom_position.z[i], 1.0f};  // world coord
                if (math::distance2(com, vec3(wc)) > cutoff_squared) continue;
                const vec4 vc = M * wc;  // volume coord [0,1]
                const vec3 tc = apply_pbc(vec3(vc));
                const ivec3 c = (ivec3)(tc * vec3(vol->dim - 1) + 0.5f);
                ASSERT(-1 < c.x && c.x < vol->dim.x);
                ASSERT(-1 < c.y && c.y < vol->dim.y);
                ASSERT(-1 < c.z && c.z < vol->dim.z);
                const i32 voxel_idx = c.z * vol->dim.x * vol->dim.y + c.y * vol->dim.x + c.x;
                platform::atomic_fetch_and_add(vol->voxel_data.ptr + voxel_idx, 1);
                // vol->voxel_data[voxel_idx]++;
            }
        }
        */
    }
}

static void append_trajectory_density(Volume* vol, Bitfield atom_mask, const MoleculeTrajectory& traj, const mat4& world_to_volume_texture) {
    ASSERT(vol);
    for (int frame_idx = 0; frame_idx < traj.num_frames; frame_idx++) {
        const auto& frame = get_trajectory_frame(traj, frame_idx);
        const mat4& M = world_to_volume_texture;

        // @TODO: Replace when proper iterators are implemented for bitfield.
        for (int i = 0; i < (int)atom_mask.size(); i++) {
            if (bitfield::get_bit(atom_mask, i)) {
                const vec4 wc = {frame.atom_position.x[i], frame.atom_position.y[i], frame.atom_position.z[i], 1.0f};  // world coord
                const vec4 vc = M * wc;                                                                                // volume coord [0,1]
                const vec4 tc = math::fract(vc);                                                                       // PBC
                const ivec3 c = vec3(tc) * (vec3)vol->dim;
                const i32 voxel_idx = c.z * vol->dim.x * vol->dim.y + c.y * vol->dim.x + c.x;
                vol->voxel_data[voxel_idx]++;
            }
        }
    }
}

static void draw_density_volume_window(ApplicationData* data) {
    const ImVec2 button_size = {250, 20};
    ImGui::Begin("Density Volume", &data->density_volume.show_window);
    ImGui::Checkbox("Enabled", &data->density_volume.enabled);
    if (data->density_volume.enabled) {
        ImGui::Separator();
        ImGui::Checkbox("Direct Volume Rendering", &data->density_volume.dvr.enabled);
        if (data->density_volume.dvr.enabled) {
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

            ImGui::SliderFloat("Density Scaling", &data->density_volume.dvr.density_scale, 0.001f, 10000.f, "%.3f", 5.f);
            ImGui::SliderFloat("Alpha Scaling", &data->density_volume.dvr.tf.alpha_scale, 0.001f, 10.f, "%.3f", 3.f);
        }

        ImGui::Separator();
        ImGui::Checkbox("Isosurface Rendering", &data->density_volume.iso.enabled);
        for (int i = 0; i < data->density_volume.iso.isosurfaces.count; ++i) {
            ImGui::PushID(i);
            ImGui::SliderFloat("##Isovalue", &data->density_volume.iso.isosurfaces.values[i], 0.0f, 10.f, "%.3f", 3.f);
            if (ImGui::IsItemDeactivatedAfterEdit()) {
                sort(data->density_volume.iso.isosurfaces);
            }
            ImGui::SameLine();
            ImGui::ColorEdit4("##Color", &data->density_volume.iso.isosurfaces.colors[i][0],
                              ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_Float);
            ImGui::PopID();
        }

        if ((data->density_volume.iso.isosurfaces.count < data->density_volume.iso.isosurfaces.MaxCount) &&
            ImGui::Button("Add Isosurface", button_size)) {
            insert(data->density_volume.iso.isosurfaces, 0.1f, {0.2f, 0.1f, 0.9f, 1.0f});
            sort(data->density_volume.iso.isosurfaces);
        }
        if (ImGui::Button("Clear Isosurfaces", button_size)) {
            clear(data->density_volume.iso.isosurfaces);
        }

        ImGui::Separator();
        ImGui::Text("Clip planes");
        ImGui::RangeSliderFloat("x", &data->density_volume.clip_planes.min.x, &data->density_volume.clip_planes.max.x, 0.0f, 1.0f);
        ImGui::RangeSliderFloat("y", &data->density_volume.clip_planes.min.y, &data->density_volume.clip_planes.max.y, 0.0f, 1.0f);
        ImGui::RangeSliderFloat("z", &data->density_volume.clip_planes.min.z, &data->density_volume.clip_planes.max.z, 0.0f, 1.0f);

        ImGui::Separator();
        ImGui::Text("Resolution Scale");
        {
            static f32 x = 1.0f;
            if (ImGui::SliderFloat("Scale", &x, 0.125f, 8.0f, "%.3f", 1.0f)) {
                f32 rounded_x;
                f32 rounded_cur;
                if (x < 1.0f) {
                    rounded_x = math::pow(2.0f, math::round(math::log(x) / math::log(2.0f)));
                    rounded_cur = math::pow(2.0f, math::round(math::log(data->density_volume.resolution_scale) / math::log(2.0f)));
                } else {
                    rounded_x = math::floor(x);
                    rounded_cur = math::floor(data->density_volume.resolution_scale);
                }

                if (rounded_x != rounded_cur) {
                    data->density_volume.resolution_scale = rounded_x;
                    init_density_volume(data);
                }
            }
        }

        ImGui::Separator();
        ImGui::Text("Ensemble");
        if (!data->ensemble_tracking.structures.empty()) {
            static char filter_buf[128] = {};
            ImGui::Checkbox("Superimpose ensemble", &data->ensemble_tracking.superimpose_structures);
            ImGui::Combo("Tracking mode", (int*)(&data->ensemble_tracking.tracking_mode), "Absolute\0Relative\0\0");

            ImGui::InputText("Density filter", filter_buf, 128, ImGuiInputTextFlags_EnterReturnsTrue);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Speficy a filter for what should be included in the density computation");
            }
            static float cutoff = 10.0f;
            ImGui::SliderFloat("Radial Cutoff", &cutoff, 1.0f, 100.0f, "%.3f", 2.0f);
            {
                ImGuiStyle& style = ImGui::GetStyle();
                const ImVec2 button_sz(32, 32);
                const int button_count = (int)data->ensemble_tracking.structures.size();
                const float window_visible_x2 = ImGui::GetWindowPos().x + ImGui::GetWindowContentRegionMax().x;
                const Bitfield ensemble_mask = data->ensemble_tracking.atom_mask;
                bool begin_popup = false;

                for (int i = 0; i < button_count; i++) {
                    auto& structure = data->ensemble_tracking.structures[i];
                    char lbl[4];
                    snprintf(lbl, 4, "%i", i);
                    ImGui::PushID(i);
                    const bool apply_style = structure.enabled;
                    if (apply_style) ImGui::PushStyleColor(ImGuiCol_Button, style.Colors[ImGuiCol_ButtonActive]);
                    if (ImGui::Button(lbl, button_sz)) {
                        structure.enabled = !structure.enabled;
                    }
                    if (ImGui::IsItemHovered()) {
                        auto mask = data->selection.current_highlight_mask;
                        bitfield::clear_all(mask);
                        for (i64 j = 0; j < ensemble_mask.size(); ++j) {
                            if (ensemble_mask[j]) {
                                bitfield::set_bit(data->selection.current_highlight_mask, structure.offset + j);
                            }
                        }
                        data->gpu_buffers.dirty.selection = true;
                    }

                    if (ImGui::IsItemClicked(1)) {
                        begin_popup = true;
                    }

                    if (apply_style) ImGui::PopStyleColor();

                    const float last_button_x2 = ImGui::GetItemRectMax().x;
                    const float next_button_x2 = last_button_x2 + style.ItemSpacing.x + button_sz.x;
                    if (i + 1 < button_count && next_button_x2 < window_visible_x2) ImGui::SameLine();
                    ImGui::PopID();
                }

                if (begin_popup) ImGui::OpenPopup("EnsemblePopup");
                if (ImGui::BeginPopup("EnsemblePopup")) {
                    if (ImGui::MenuItem("Select All")) {
                        for (auto& item : data->ensemble_tracking.structures) {
                            item.enabled = true;
                        }
                    }
                    if (ImGui::MenuItem("Deselect All")) {
                        for (auto& item : data->ensemble_tracking.structures) {
                            item.enabled = false;
                        }
                    }
                    ImGui::EndPopup();
                }
            }

            if (ImGui::Button("Compute Density")) {
                task_system::create_task("Big Density", [data](task_system::TaskSetRange, task_system::TaskData) {
                    Bitfield filter_mask;
                    bitfield::init(&filter_mask, data->dynamic.molecule.atom.count);
                    defer { bitfield::free(&filter_mask); };

                    DynamicArray<StoredSelection> sel;
                    sel.push_back({"current", data->selection.current_selection_mask});
                    for (const auto& s : data->selection.stored_selections) {
                        sel.push_back({s.name, s.atom_mask});
                    }
                    const bool filter_ok = filter::compute_filter_mask(filter_mask, filter_buf, data->dynamic.molecule, sel);

                    if (bitfield::number_of_bits_set(filter_mask) > 0) {
                        const auto traj = data->dynamic.trajectory;
                        const auto ensemble_structures = data->ensemble_tracking.structures;

                        data->density_volume.volume_data_mutex.lock();
                        clear_volume(&data->density_volume.volume);

                        auto t0 = platform::get_time();
                        task_system::ID task_system = task_system::create_task(
                            "Small Density", (u32)ensemble_structures.size(),
                            [data, filter_mask, &ensemble_structures](task_system::TaskSetRange range, task_system::TaskData) {
                                for (u32 i = range.beg; i < range.end; ++i) {
                                    auto& structure = ensemble_structures[i];
                                    if (!structure.enabled) continue;
                                    const structure_tracking::TrackingData* tracking_data = structure_tracking::get_tracking_data(structure.id);
                                    if (!tracking_data) {
                                        // LOG_ERROR("Could not find tracking data for structure: %u", structure.id);
                                        continue;
                                    }
                                    append_trajectory_density(&data->density_volume.volume, filter_mask, data->dynamic.trajectory,
                                                              structure.alignment_matrix, *tracking_data, data->ensemble_tracking.tracking_mode,
                                                              cutoff);
                                }
                            });
                        task_system::wait_for_task(task_system);
                        auto t1 = platform::get_time();

                        data->density_volume.volume.voxel_range = {data->density_volume.volume.voxel_data[0],
                                                                   data->density_volume.volume.voxel_data[0]};
                        for (auto& v : data->density_volume.volume.voxel_data) {
                            if (v < data->density_volume.volume.voxel_range.min) data->density_volume.volume.voxel_range.min = v;
                            if (v > data->density_volume.volume.voxel_range.max) data->density_volume.volume.voxel_range.max = v;
                        }

                        LOG_NOTE("Done!");
                        LOG_NOTE("Time taken: %.2f ms", platform::compute_delta_ms(t0, t1));
                        LOG_NOTE("Max density: %u", data->density_volume.volume.voxel_range.max);
                        data->density_volume.texture.dirty = true;
                        data->density_volume.volume_data_mutex.unlock();
                    }
                });
            }
        }

        if (ImGui::Button("Create Ensemble", button_size)) {
            ImGui::OpenPopup("Create Ensemble Popup");
        }

        if (ImGui::BeginPopup("Create Ensemble Popup")) {
            static char reference_buf[128] = {};

            bool compute = false;
            compute |= ImGui::InputText("Reference structure filter", reference_buf, 128, ImGuiInputTextFlags_EnterReturnsTrue);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip(
                    "Specify the filter of the reference structure for the ensamble. All structure which matches this reference will be included "
                    "in the ensamble");
            }
            if (ImGui::Button("Create ensemble")) {
                ImGui::CloseCurrentPopup();
                Bitfield ref_mask;
                bitfield::init(&ref_mask, data->dynamic.molecule.atom.count);
                defer { bitfield::free(&ref_mask); };

                DynamicArray<StoredSelection> sel;
                sel.push_back({"current", data->selection.current_selection_mask});
                for (const auto& s : data->selection.stored_selections) {
                    sel.push_back({s.name, s.atom_mask});
                }
                const bool ref_ok = filter::compute_filter_mask(ref_mask, reference_buf, data->dynamic.molecule, sel);
                if (ref_ok) {
                    const auto& dyn = data->dynamic;
                    const auto& mol = data->dynamic.molecule;
                    const auto& traj = data->dynamic.trajectory;

                    Bitfield& ensemble_mask = data->ensemble_tracking.atom_mask;
                    DynamicArray<EnsembleStructure>& ensemble_structures = data->ensemble_tracking.structures;

                    const auto first_bit = bitfield::find_first_bit_set(ref_mask);
                    const auto last_bit = bitfield::find_last_bit_set(ref_mask);
                    if (first_bit == -1 || last_bit == -1) {
                        LOG_ERROR("Bitfield was empty");
                        return;
                    }

                    const auto mask_offset = first_bit;
                    const auto mask_size = last_bit - first_bit + 1;

                    // @NOTE: Make a smaller 'ensemble' reference mask for the smaller range surrounding the structure of interest
                    bitfield::init(&ensemble_mask, mask_size);
                    bitfield::clear_all(ensemble_mask);
                    for (int i = 0; i < mask_size; i++) {
                        if (ref_mask[mask_offset + i]) bitfield::set_bit(ensemble_mask, i);
                    }

                    LOG_NOTE("Matching structures to reference...");
                    {
                        DynamicArray<int> ensemble_offsets = find_equivalent_structures(data->dynamic.molecule, ensemble_mask, (int)mask_offset);
                        if (!ensemble_structures.empty()) {
                            for (const auto& s : ensemble_structures) {
                                structure_tracking::remove_structure(s.id);
                            }
                        }
                        ensemble_structures.clear();
                        ensemble_structures.push_back(
                            {(int)mask_offset, structure_tracking::create_structure(), mat4(1), mat4(1), true});  // Reference structure at index 0
                        for (const auto& offset : ensemble_offsets) {
                            ensemble_structures.push_back({(int)offset, structure_tracking::create_structure(), mat4(1), mat4(1), true});
                        }
                    }
                    LOG_NOTE("Done!");

                    {
                        LOG_NOTE("Computing Internal Reference Frames...");
                        auto t0 = platform::get_time();
                        u32 completed_count = 0;
                        task_system::ID task_system = task_system::create_task(
                            "Computing Internal Reference Frames...", (u32)ensemble_structures.size(),
                            [dyn, ensemble_mask, &ensemble_structures, &completed_count](task_system::TaskSetRange range, task_system::TaskData) {
                                for (u32 i = range.beg; i < range.end; ++i) {
                                    structure_tracking::compute_trajectory_transform_data(ensemble_structures[i].id, dyn, ensemble_mask,
                                                                                          ensemble_structures[i].offset);
                                }
                            });
                        task_system::wait_for_task(task_system);
                        auto t1 = platform::get_time();

                        LOG_NOTE("Done!");
                        LOG_NOTE("Time taken: %.2f ms", platform::compute_delta_ms(t0, t1));
                    }

                    {
                        LOG_NOTE("Computing structure Alignment Matrices...");
                        auto t0 = platform::get_time();
                        const auto atom_count = bitfield::number_of_bits_set(ensemble_mask);
                        void* mem = TMP_MALLOC(sizeof(float) * 7 * atom_count);
                        defer { TMP_FREE(mem); };
                        ASSERT(mem);

                        float* ref_x = (float*)mem + 0 * atom_count;
                        float* ref_y = (float*)mem + 1 * atom_count;
                        float* ref_z = (float*)mem + 2 * atom_count;

                        float* x = (float*)mem + 3 * atom_count;
                        float* y = (float*)mem + 4 * atom_count;
                        float* z = (float*)mem + 5 * atom_count;

                        float* mass = (float*)mem + 6 * atom_count;

                        const auto frame0 = get_trajectory_frame(traj, 0);
                        bitfield::gather_masked(ref_x, frame0.atom_position.x, ensemble_mask, ensemble_structures[0].offset);
                        bitfield::gather_masked(ref_y, frame0.atom_position.y, ensemble_mask, ensemble_structures[0].offset);
                        bitfield::gather_masked(ref_z, frame0.atom_position.z, ensemble_mask, ensemble_structures[0].offset);
                        bitfield::gather_masked(mass, mol.atom.mass, ensemble_mask, ensemble_structures[0].offset);
                        const vec3 ref_com = compute_com(ref_x, ref_y, ref_z, mass, atom_count);
                        const mat3 PCA = structure_tracking::get_tracking_data(ensemble_structures[0].id)->simulation_box_aligned_pca;

                        ensemble_structures[0].alignment_matrix = PCA;
                        for (i64 i = 1; i < ensemble_structures.size(); i++) {
                            auto& structure = ensemble_structures[i];
                            bitfield::gather_masked(x, frame0.atom_position.x, ensemble_mask, structure.offset);
                            bitfield::gather_masked(y, frame0.atom_position.y, ensemble_mask, structure.offset);
                            bitfield::gather_masked(z, frame0.atom_position.z, ensemble_mask, structure.offset);
                            const vec3 com = compute_com(x, y, z, mass, atom_count);
                            const mat3 R = structure_tracking::compute_rotation(x, y, z, ref_x, ref_y, ref_z, mass, atom_count, com, ref_com);
                            structure.alignment_matrix = PCA * math::transpose(R);
                        }
                        auto t1 = platform::get_time();
                        LOG_NOTE("Done!");
                        LOG_NOTE("Time taken: %.2f ms", platform::compute_delta_ms(t0, t1));
                    }
                }
            }
            ImGui::EndPopup();
        }
        ImGui::Separator();
        if (ImGui::Button("Export Density...", ImVec2(250, 20))) {
            auto res = platform::file_dialog(platform::FileDialogFlags_Open, {}, "raw");
            if (res.result == platform::FileDialogResult::Ok) {
                volume::write_to_file(data->density_volume.volume, res.path);
                LOG_NOTE("Wrote density volume");
            }
        }
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
    if (!fbo->deferred.emissive) glGenTextures(1, &fbo->deferred.emissive);
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
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16F, width, height, 0, GL_RG, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, fbo->deferred.emissive);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R11F_G11F_B10F, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
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
    glBufferData(GL_PIXEL_PACK_BUFFER, 4, 0, GL_DYNAMIC_READ);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, fbo->pbo_picking.color[1]);
    glBufferData(GL_PIXEL_PACK_BUFFER, 4, 0, GL_DYNAMIC_READ);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, fbo->pbo_picking.depth[0]);
    glBufferData(GL_PIXEL_PACK_BUFFER, 4, 0, GL_DYNAMIC_READ);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, fbo->pbo_picking.depth[1]);
    glBufferData(GL_PIXEL_PACK_BUFFER, 4, 0, GL_DYNAMIC_READ);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    glBindTexture(GL_TEXTURE_2D, 0);

    fbo->width = width;
    fbo->height = height;

    const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2,
                                   GL_COLOR_ATTACHMENT3, GL_COLOR_ATTACHMENT4, GL_COLOR_ATTACHMENT5};

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo->deferred.fbo);
    if (attach_textures_deferred) {
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, fbo->deferred.depth, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_STENCIL_ATTACHMENT, GL_TEXTURE_2D, fbo->deferred.depth, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, fbo->deferred.color, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, fbo->deferred.normal, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, fbo->deferred.velocity, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, GL_TEXTURE_2D, fbo->deferred.emissive, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT4, GL_TEXTURE_2D, fbo->deferred.post_tonemap, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT5, GL_TEXTURE_2D, fbo->deferred.picking, 0);
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
    if (fbo->deferred.emissive) glDeleteTextures(1, &fbo->deferred.emissive);
    if (fbo->deferred.post_tonemap) glDeleteTextures(1, &fbo->deferred.post_tonemap);
    if (fbo->deferred.picking) glDeleteTextures(1, &fbo->deferred.picking);

    if (fbo->pbo_picking.color[0]) glDeleteBuffers(2, fbo->pbo_picking.color);
    if (fbo->pbo_picking.depth[0]) glDeleteBuffers(2, fbo->pbo_picking.depth);
}

static void init_molecule_buffers(ApplicationData* data) {
    ASSERT(data);
    const auto& mol = data->dynamic.molecule;

    DynamicArray<u32> backbone_index_data;
    DynamicArray<u32> control_point_index_data;
    DynamicArray<u32> spline_index_data;

    {
        i32 control_idx = 0;
        i32 spline_idx = 0;
        for (const auto& seq : mol.backbone.sequences) {
            const auto backbone = get_backbone(mol, seq);
            const i32 count = (i32)backbone.size();

            /*
                // These indices are needed to compute the backbone angles
                                phi = math::dihedral_angle(c_idx[i-1], n_idx[i], ca_idx[i], c_idx[i]);
                                psi = math::dihedral_angle(n_idx[i], ca_idx[i], c_idx[i], n_idx[i+1]);
            */

            for (i32 i = 0; i < count; i++) {
                const bool first = (i == 0);
                const bool last = (i == count - 1);

                const auto ca_i = backbone[i].ca_idx;
                const auto c_i = backbone[i].c_idx;
                const auto o_i = backbone[i].o_idx;
                const auto n_i = backbone[i].n_idx;
                const auto c_im1 = backbone[math::max(i - 1, 0)].c_idx;
                const auto n_ip1 = backbone[math::min(i + 1, count - 1)].n_idx;

                backbone_index_data.push_back(ca_i);
                backbone_index_data.push_back(c_i);
                backbone_index_data.push_back(o_i);
                backbone_index_data.push_back(n_i);
                backbone_index_data.push_back(c_im1);
                backbone_index_data.push_back(n_ip1);
                control_point_index_data.push_back(control_idx);

                // @NOTE: Pad with extra index on first and last to help cubic spline construction
                if (first || last) {
                    control_point_index_data.push_back(control_idx);
                }
                control_idx++;

                // @NOTE: For every 'base' control point we generate N spline control points (for subdivision)
                if (!last) {
                    for (i32 j = 0; j < SPLINE_SUBDIVISION_COUNT; j++) {
                        spline_index_data.push_back(spline_idx);
                        spline_idx++;
                    }
                } else {
                    spline_index_data.push_back(0xFFFFFFFFU);
                }
            }
            control_point_index_data.push_back(0xFFFFFFFFU);
        }
    }

    data->gpu_buffers.backbone.num_backbone_segment_indices = (i32)backbone_index_data.size();
    data->gpu_buffers.backbone.num_control_point_indices = (i32)control_point_index_data.size();
    data->gpu_buffers.backbone.num_spline_indices = (i32)spline_index_data.size();

    LOG_NOTE("num backbone segment indices: %i", (i32)backbone_index_data.size());
    LOG_NOTE("num control point indices: %i", (i32)control_point_index_data.size());
    LOG_NOTE("num spline indices: %i", (i32)spline_index_data.size());

    const i64 num_backbone_segments = backbone_index_data.size() / 6;
    const i64 position_buffer_size = mol.atom.count * 3 * sizeof(float);
    const i64 velocity_buffer_size = mol.atom.count * 3 * sizeof(float);
    const i64 radius_buffer_size = mol.atom.count * sizeof(float);
    const i64 bond_buffer_size = mol.covalent_bonds.size() * sizeof(u32) * 2;
    const i64 selection_buffer_size = mol.atom.count * sizeof(u8);
    const i64 control_point_buffer_size = num_backbone_segments * sizeof(draw::ControlPoint);
    const i64 spline_buffer_size = control_point_buffer_size * SPLINE_SUBDIVISION_COUNT;

    if (!data->gpu_buffers.position) glGenBuffers(1, &data->gpu_buffers.position);
    if (!data->gpu_buffers.old_position) glGenBuffers(1, &data->gpu_buffers.old_position);
    if (!data->gpu_buffers.velocity) glGenBuffers(1, &data->gpu_buffers.velocity);
    if (!data->gpu_buffers.radius) glGenBuffers(1, &data->gpu_buffers.radius);
    if (!data->gpu_buffers.bond) glGenBuffers(1, &data->gpu_buffers.bond);
    if (!data->gpu_buffers.selection) glGenBuffers(1, &data->gpu_buffers.selection);
    if (!data->gpu_buffers.backbone.backbone_segment_index) glGenBuffers(1, &data->gpu_buffers.backbone.backbone_segment_index);
    if (!data->gpu_buffers.backbone.control_point) glGenBuffers(1, &data->gpu_buffers.backbone.control_point);
    if (!data->gpu_buffers.backbone.control_point_index) glGenBuffers(1, &data->gpu_buffers.backbone.control_point_index);
    if (!data->gpu_buffers.backbone.spline) glGenBuffers(1, &data->gpu_buffers.backbone.spline);
    if (!data->gpu_buffers.backbone.spline_index) glGenBuffers(1, &data->gpu_buffers.backbone.spline_index);

    glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.position);
    glBufferData(GL_ARRAY_BUFFER, position_buffer_size, NULL, GL_DYNAMIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.old_position);
    glBufferData(GL_ARRAY_BUFFER, position_buffer_size, NULL, GL_DYNAMIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.velocity);
    glBufferData(GL_ARRAY_BUFFER, velocity_buffer_size, NULL, GL_DYNAMIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.radius);
    glBufferData(GL_ARRAY_BUFFER, radius_buffer_size, mol.atom.radius, GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.bond);
    glBufferData(GL_ARRAY_BUFFER, bond_buffer_size, mol.covalent_bonds.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.selection);
    glBufferData(GL_ARRAY_BUFFER, selection_buffer_size, NULL, GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.backbone.backbone_segment_index);
    glBufferData(GL_ARRAY_BUFFER, backbone_index_data.size_in_bytes(), backbone_index_data.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.backbone.control_point);
    glBufferData(GL_ARRAY_BUFFER, control_point_buffer_size, NULL, GL_DYNAMIC_COPY);

    glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.backbone.control_point_index);
    glBufferData(GL_ARRAY_BUFFER, control_point_index_data.size_in_bytes(), control_point_index_data.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.backbone.spline);
    glBufferData(GL_ARRAY_BUFFER, spline_buffer_size, NULL, GL_DYNAMIC_COPY);

    glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.backbone.spline_index);
    glBufferData(GL_ARRAY_BUFFER, spline_index_data.size_in_bytes(), spline_index_data.data(), GL_STATIC_DRAW);

    // Experimental
    draw::generate_aabb_buffer(&data->gpu_buffers.experimental.aabb, (i32)mol.residues.count);
    draw::generate_residue_buffer(&data->gpu_buffers.experimental.residue, (i32)mol.residues.count);
    draw::generate_visibility_buffer(&data->gpu_buffers.experimental.visibility, (i32)mol.residues.count);

    DynamicArray<ivec2> res_offset_count(mol.residues.count);
    for (i64 i = 0; i < mol.residues.count; i++) {
        res_offset_count[i] = {(i32)mol.residues[i].atom_range.beg, (i32)mol.residues[i].atom_range.ext()};
    }

    glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.experimental.residue);
    glBufferSubData(GL_ARRAY_BUFFER, 0, res_offset_count.size_in_bytes(), res_offset_count.data());

    data->gpu_buffers.dirty.position = true;
    data->gpu_buffers.dirty.selection = true;
    data->gpu_buffers.dirty.backbone = true;

    copy_molecule_data_to_buffers(data);
    data->gpu_buffers.dirty.position = true;

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

static void free_molecule_buffers(ApplicationData* data) {
    ASSERT(data);

    auto free_buffer = [](GLuint& id) {
        if (id) {
            glDeleteBuffers(1, &id);
            id = 0;
        }
    };

    free_buffer(data->gpu_buffers.position);
    free_buffer(data->gpu_buffers.old_position);
    free_buffer(data->gpu_buffers.velocity);
    free_buffer(data->gpu_buffers.radius);
    free_buffer(data->gpu_buffers.selection);
    free_buffer(data->gpu_buffers.backbone.backbone_segment_index);
    free_buffer(data->gpu_buffers.backbone.control_point);
    free_buffer(data->gpu_buffers.backbone.control_point_index);
    free_buffer(data->gpu_buffers.backbone.spline);
    free_buffer(data->gpu_buffers.backbone.spline_index);
    free_buffer(data->gpu_buffers.bond);

    // experimental
    free_buffer(data->gpu_buffers.experimental.residue);
    free_buffer(data->gpu_buffers.experimental.aabb);
    free_buffer(data->gpu_buffers.experimental.visibility);
}

static void clear_buffer(GLuint buffer) {
    GLint size;
    glBindBuffer(GL_COPY_WRITE_BUFFER, buffer);
    glGetBufferParameteriv(GL_COPY_WRITE_BUFFER, GL_BUFFER_SIZE, &size);
    glBufferSubData(GL_COPY_WRITE_BUFFER, 0, size, NULL);
    glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
}

static void copy_buffer(GLuint dst_buffer, GLuint src_buffer) {
    GLint size;
    glBindBuffer(GL_COPY_READ_BUFFER, src_buffer);
    glBindBuffer(GL_COPY_WRITE_BUFFER, dst_buffer);
    glGetBufferParameteriv(GL_COPY_READ_BUFFER, GL_BUFFER_SIZE, &size);
    glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0, size);
    glBindBuffer(GL_COPY_READ_BUFFER, 0);
    glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
}

static void copy_molecule_data_to_buffers(ApplicationData* data) {
    ASSERT(data);
    const auto N = data->dynamic.molecule.atom.count;

    if (data->gpu_buffers.dirty.position) {
        data->gpu_buffers.dirty.position = false;
        const float* pos_x = data->dynamic.molecule.atom.position.x;
        const float* pos_y = data->dynamic.molecule.atom.position.y;
        const float* pos_z = data->dynamic.molecule.atom.position.z;

        // Copy position to old_position
        copy_buffer(data->gpu_buffers.old_position, data->gpu_buffers.position);

        // Update data inside position buffer
        glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.position);

        float* pos_gpu = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        if (!pos_gpu) {
            LOG_ERROR("Could not map position buffer");
            return;
        }

        for (i64 i = 0; i < N; i++) {
            pos_gpu[i * 3 + 0] = pos_x[i];
            pos_gpu[i * 3 + 1] = pos_y[i];
            pos_gpu[i * 3 + 2] = pos_z[i];
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }

    if (data->gpu_buffers.dirty.selection) {
        data->gpu_buffers.dirty.selection = false;

        glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.selection);
        u8* sel_gpu = (u8*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

        if (!sel_gpu) {
            LOG_ERROR("Could not map selection buffer");
            return;
        }

        for (i64 i = 0; i < N; i++) {
            const bool selection = bitfield::get_bit(data->selection.current_selection_mask, i);
            const bool highlight = bitfield::get_bit(data->selection.current_highlight_mask, i);
            sel_gpu[i] = (selection ? 0x2 : 0) | (highlight ? 0x1 : 0);
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// #moleculedata
static void free_trajectory_data(ApplicationData* data) {
    ASSERT(data);
    if (data->dynamic.trajectory) {
        if (data->tasks.load_trajectory) task_system::interrupt_and_wait(data->tasks.load_trajectory);
        stats::signal_stop_and_wait();
        // close_file_handle(&data->dynamic.trajectory);
        free_trajectory(&data->dynamic.trajectory);
    }
}

static void free_molecule_data(ApplicationData* data) {
    ASSERT(data);
    if (data->dynamic.molecule) {
        data->files.molecule = "";
        free_molecule_structure(&data->dynamic.molecule);
    }
    if (data->dynamic.trajectory) {
        data->files.trajectory = "";
        free_trajectory_data(data);
    }
    free_molecule_buffers(data);
    data->hydrogen_bonds.bonds.clear();
    data->hydrogen_bonds.dirty = true;
    data->gpu_buffers.dirty.backbone = true;
    bitfield::clear_all(data->selection.current_selection_mask);
    bitfield::clear_all(data->selection.current_highlight_mask);
}

static void init_molecule_data(ApplicationData* data) {
    if (data->dynamic.molecule) {
        const auto atom_count = data->dynamic.molecule.atom.count;
        bitfield::init(&data->selection.current_selection_mask, atom_count);
        bitfield::init(&data->selection.current_highlight_mask, atom_count);
        init_molecule_buffers(data);
        data->picking.idx = NO_PICKING_IDX;
        data->selection.hovered = -1;
        data->selection.right_clicked = -1;
        /*
        const auto& mol = data->dynamic.molecule;
        LOG_NOTE("Computing VDW sdf...");
        auto t0 = platform::get_time();
        draw::sdf::compute_vdw_sdf(mol.atom.position.x, mol.atom.position.y, mol.atom.position.z, mol.atom.radius, mol.atom.count);
        auto t1 = platform::get_time();
        LOG_NOTE("Done! Time taken: %.2f ms", platform::compute_delta_ms(t0, t1));
        */
    }
}

static void init_trajectory_data(ApplicationData* data) {
    if (data->dynamic.trajectory) {
        data->time_filter.range = {0, (float)data->dynamic.trajectory.num_frames};
        data->playback.time = math::clamp(data->playback.time, (f64)data->time_filter.range.min, (f64)data->time_filter.range.max);
        data->playback.frame = math::clamp((i32)data->playback.time, 0, data->dynamic.trajectory.num_frames);
        i32 frame = data->playback.frame;

        data->simulation_box.box = get_trajectory_frame(data->dynamic.trajectory, data->playback.frame).box;
        const Array<float> pos_x = get_trajectory_position_x(data->dynamic.trajectory, frame);
        const Array<float> pos_y = get_trajectory_position_y(data->dynamic.trajectory, frame);
        const Array<float> pos_z = get_trajectory_position_z(data->dynamic.trajectory, frame);

        memcpy(data->dynamic.molecule.atom.position.x, pos_x.data(), pos_x.size_in_bytes());
        memcpy(data->dynamic.molecule.atom.position.y, pos_y.data(), pos_y.size_in_bytes());
        memcpy(data->dynamic.molecule.atom.position.z, pos_z.data(), pos_z.size_in_bytes());

        data->gpu_buffers.dirty.position = true;
        copy_molecule_data_to_buffers(data);
        copy_buffer(data->gpu_buffers.old_position, data->gpu_buffers.position);

        init_density_volume(data);
    }
}

static void interrupt_async_tasks(ApplicationData* data) {
    if (data->tasks.load_trajectory) {
        task_system::interrupt_task(data->tasks.load_trajectory);
    }
    task_system::wait_for_task(data->tasks.load_trajectory);
    data->tasks.load_trajectory = 0;
}

static bool load_trajectory_data(ApplicationData* data, CStringView filename) {
    interrupt_async_tasks(data);
    free_trajectory_data(data);
    data->files.trajectory = filename;
    data->playback.time = 0;
    data->playback.frame = 0;

    LOG_NOTE("Attempting to load trajectory data from file '%.*s", filename.count, filename.ptr);
    data->tasks.load_trajectory = task_system::create_task("Loading trajectory", [data](task_system::TaskSetRange, task_system::TaskData) {
        const auto t0 = platform::get_time();
        if (!load::traj::load_trajectory(&data->dynamic.trajectory, (i32)data->dynamic.molecule.atom.count, data->files.trajectory)) {
            LOG_ERROR("Failed to load trajectory");
            free_trajectory_data(data);
            data->files.trajectory = "";
            return;
        }
        const auto t1 = platform::get_time();
        LOG_NOTE("Trajectory loaded in %.2f seconds", platform::compute_delta_ms(t0, t1) / 1000.0f);
        init_trajectory_data(data);
        on_trajectory_load_complete(data);
    });
    return true;
}

static bool load_dataset_from_file(ApplicationData* data, CStringView filename) {
    ASSERT(data);

    if (filename) {
        if (load::mol::is_extension_supported(filename)) {
            interrupt_async_tasks(data);
            free_molecule_data(data);
            free_trajectory_data(data);
            data->files.molecule = filename;
            data->files.trajectory = "";

            LOG_NOTE("Attempting to load molecular data from file '%.*s'", filename.count, filename.ptr);
            if (!load::mol::load_molecule(&data->dynamic.molecule, filename)) {
                LOG_ERROR("Failed to load molecular data");
                data->files.molecule = "";
                return false;
            }
            init_molecule_data(data);

            // @NOTE: Some files contain both atomic coordinates and trajectory
            i32 num_frames = 0;
            if (load::traj::is_extension_supported(filename) && load::traj::read_num_frames(&num_frames, filename) && num_frames > 1) {
                // @TODO: Halt here and query user if this is actual frames or individual models of a larger assembly.
                return load_trajectory_data(data, filename);
            }
            return true;
        } else if (load::traj::is_extension_supported(filename)) {
            if (!data->dynamic.molecule) {
                LOG_ERROR("Before loading a trajectory a molecular data needs to be present");
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
    else if (compare(str, "BALL_AND_STICK"))
        return RepresentationType::BallAndStick;
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
        case RepresentationType::BallAndStick:
            return "BALL_AND_STICK";
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
    stats::remove_all_properties();
    clear(data->density_volume.iso.isosurfaces);

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
        } else if (compare_n(line, "[Property]", 10)) {
            StringBuffer<256> name, args;
            while (c_txt && c_txt[0] != '[' && (line = extract_line(c_txt))) {
                if (compare_n(line, "Name=", 5)) name = trim(line.substr(5));
                if (compare_n(line, "Args=", 5)) args = trim(line.substr(5));
            }
            stats::create_property(name, args);
        } else if (compare_n(line, "[RenderSettings]", 16)) {
            while (c_txt && c_txt[0] != '[' && (line = extract_line(c_txt))) {
                if (compare_n(line, "SsaoEnabled=", 12)) data->visuals.ssao.enabled = to_int(trim(line.substr(12))) != 0;
                if (compare_n(line, "SsaoIntensity=", 14)) data->visuals.ssao.intensity = to_float(trim(line.substr(14)));
                if (compare_n(line, "SsaoRadius=", 11)) data->visuals.ssao.radius = to_float(trim(line.substr(11)));
                if (compare_n(line, "SsaoBias=", 9)) data->visuals.ssao.bias = to_float(trim(line.substr(9)));
                if (compare_n(line, "DofEnabled=", 11)) data->visuals.dof.enabled = to_int(trim(line.substr(11))) != 0;
                if (compare_n(line, "DofFocusScale=", 14)) data->visuals.dof.focus_scale = to_float(trim(line.substr(14)));

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
            }
        } else if (compare_n(line, "[Isosurface]", 21)) {
            float isovalue = -1.0f;
            vec4 isocolor{0.0f};
            while (c_txt && c_txt[0] != '[' && (line = extract_line(c_txt))) {
                if (compare_n(line, "Isovalue=", 9)) isovalue = to_float(trim(line.substr(9)));
                if (compare_n(line, "IsosurfaceColor=", 16)) isocolor = to_vec4(trim(line.substr(16)));
            }
            insert(data->density_volume.iso.isosurfaces, isovalue, isocolor);
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
                    data->view.trackball_state.distance = to_float(trim(line.substr(9)));
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

    data->playback = {};
    reset_view(data, false, true);
    reset_representations(data);
}

static void save_workspace(ApplicationData* data, CStringView file) {
    FILE* fptr = fopen(file, "w");
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
    for (const auto prop : stats::get_properties()) {
        fprintf(fptr, "[Property]\n");
        fprintf(fptr, "Name=%s\n", prop->name_buf.cstr());
        fprintf(fptr, "Args=%s\n", prop->args_buf.cstr());
        fprintf(fptr, "\n");
    }

    fprintf(fptr, "[RenderSettings]\n");
    fprintf(fptr, "SsaoEnabled=%i\n", data->visuals.ssao.enabled ? 1 : 0);
    fprintf(fptr, "SsaoIntensity=%g\n", data->visuals.ssao.intensity);
    fprintf(fptr, "SsaoRadius=%g\n", data->visuals.ssao.radius);
    fprintf(fptr, "SsaoBias=%g\n", data->visuals.ssao.bias);
    fprintf(fptr, "\n");

    fprintf(fptr, "DofEnabled=%i\n", data->visuals.dof.enabled ? 1 : 0);
    fprintf(fptr, "DofFocusScale=%g\n", data->visuals.dof.focus_scale);
    fprintf(fptr, "\n");

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

    fprintf(fptr, "[Camera]\n");
    fprintf(fptr, "Position=%g,%g,%g\n", data->view.camera.position.x, data->view.camera.position.y, data->view.camera.position.z);
    fprintf(fptr, "Rotation=%g,%g,%g,%g\n", data->view.camera.orientation.x, data->view.camera.orientation.y, data->view.camera.orientation.z,
            data->view.camera.orientation.w);
    fprintf(fptr, "Distance=%g\n", data->view.trackball_state.distance);
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
        u32* row_t = (u32*)TMP_MALLOC(row_byte_size);
        defer { TMP_FREE(row_t); };
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
        auto ext = get_file_extension(file_res.path);
        if (!ext) {
            file_res.path += ".jpg";
            ext = "jpg";
        }
        if (ext == "jpg") {
            const int quality = 95;
            write_image_jpg(img, file_res.path, quality);
        } else if (ext == "png") {
            write_image_png(img, file_res.path);
        } else if (ext == "bmp") {
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
    update_representation(data, &rep);
    return &rep;
}

static Representation* clone_representation(ApplicationData* data, const Representation& rep) {
    ASSERT(data);
    Representation& clone = data->representations.buffer.push_back(rep);
    clone.color_buffer = 0;
    update_representation(data, &clone);
    return &clone;
}

static void remove_representation(ApplicationData* data, int idx) {
    ASSERT(data);
    ASSERT(idx < data->representations.buffer.size());

    auto& rep = data->representations.buffer[idx];
    if (rep.color_buffer) glDeleteBuffers(1, &rep.color_buffer);
    if (rep.atom_mask) bitfield::free(&rep.atom_mask);
    data->representations.buffer.remove(&rep);
}

static void recompute_atom_visibility_mask(ApplicationData* data) {
    ASSERT(data);
    if (!data->representations.atom_visibility_mask_dirty) return;

    auto& atom_visibility_mask = data->representations.atom_visibility_mask;
    if (atom_visibility_mask.size() != data->dynamic.molecule.atom.count) {
        bitfield::init(&atom_visibility_mask, data->dynamic.molecule.atom.count);
    }

    bitfield::clear_all(atom_visibility_mask);
    for (const auto& rep : data->representations.buffer) {
        if (!rep.enabled) continue;
        bitfield::or_field(atom_visibility_mask, atom_visibility_mask, rep.atom_mask);
    }

    data->representations.atom_visibility_mask_dirty = false;
}

static void update_all_representations(ApplicationData* data) {
    for (auto& rep : data->representations.buffer) {
        update_representation(data, &rep);
    }
}

static void update_representation(ApplicationData* data, Representation* rep) {
    ASSERT(data);
    ASSERT(rep);

    void* mem = TMP_MALLOC(data->dynamic.molecule.atom.count * sizeof(u32));
    defer { TMP_FREE(mem); };

    Array<u32> colors((u32*)mem, data->dynamic.molecule.atom.count);
    const auto& mol = data->dynamic.molecule;

    switch (rep->color_mapping) {
        case ColorMapping::Uniform:
            color_atoms_uniform(colors, rep->uniform_color);
            break;
        case ColorMapping::Cpk:
            color_atoms_cpk(colors, get_elements(mol));
            break;
        case ColorMapping::ResId:
            color_atoms_residue_id(colors, get_residues(mol));
            break;
        case ColorMapping::ResIndex:
            color_atoms_residue_index(colors, get_residues(mol));
            break;
        case ColorMapping::ChainId:
            color_atoms_chain_id(colors, get_chains(mol));
            break;
        case ColorMapping::ChainIndex:
            color_atoms_chain_index(colors, get_chains(mol));
            break;
        case ColorMapping::SecondaryStructure:
            color_atoms_backbone_angles(colors, get_residues(mol), mol.backbone.sequences, mol.backbone.angles, ramachandran::get_color_image());
            break;
        default:
            ASSERT(false);
            break;
    }

    if (rep->atom_mask.size() != mol.atom.count) {
        bitfield::init(&rep->atom_mask, mol.atom.count);
    }

    DynamicArray<StoredSelection> sel;
    sel.push_back({"current", data->selection.current_selection_mask});
    for (const auto& s : data->selection.stored_selections) {
        sel.push_back({s.name, s.atom_mask});
    }

    rep->filter_is_ok = filter::compute_filter_mask(rep->atom_mask, rep->filter.buffer, data->dynamic.molecule, sel);
    filter_colors(colors, rep->atom_mask);
    data->representations.atom_visibility_mask_dirty = true;

    if (!rep->color_buffer) glGenBuffers(1, &rep->color_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, rep->color_buffer);
    glBufferData(GL_ARRAY_BUFFER, colors.size_in_bytes(), colors.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
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
    const i64 N = data->dynamic.molecule.atom.count;
    const bool shift_down = data->ctx.input.key.down[Key::KEY_LEFT_SHIFT] || data->ctx.input.key.down[Key::KEY_RIGHT_SHIFT];
    const bool mouse_down = data->ctx.input.mouse.down[0] || data->ctx.input.mouse.down[1];

    Bitfield mask;
    bitfield::init(&mask, N);
    defer { bitfield::free(&mask); };

    bitfield::clear_all(mask);
    bitfield::clear_all(data->selection.current_highlight_mask);
    data->gpu_buffers.dirty.selection = true;

    // Range<int32> picking_range = {0, 0};

    if (data->picking.idx != NO_PICKING_IDX && !region_select) {
        ASSERT(0 <= data->picking.idx && data->picking.idx <= N);
        bitfield::set_bit(mask, data->picking.idx);

        switch (data->selection.level_mode) {
            case SelectionLevel::Atom:
                break;
            case SelectionLevel::Residue: {
                expand_mask(mask, data->dynamic.molecule.residues);
                break;
            }
            case SelectionLevel::Chain: {
                expand_mask(mask, data->dynamic.molecule.chains);
                break;
            }
            default:
                ASSERT(false);
                break;
        }
        memcpy(data->selection.current_highlight_mask.data(), mask.data(), mask.size_in_bytes());
        data->gpu_buffers.dirty.selection = true;
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

            for (i64 i = 0; i < N; i++) {
                if (!bitfield::get_bit(data->representations.atom_visibility_mask, i)) continue;

                // @PERF: Do the projection manually. GLM is super slow in doing the matrix vector multiplication on msvc for some reason...
                const float x = data->dynamic.molecule.atom.position.x[i];
                const float y = data->dynamic.molecule.atom.position.y[i];
                const float z = data->dynamic.molecule.atom.position.z[i];

                const float p_x = mvp[0][0] * x + mvp[1][0] * y + mvp[2][0] * z + mvp[3][0];
                const float p_y = mvp[0][1] * x + mvp[1][1] * y + mvp[2][1] * z + mvp[3][1];
                const float p_w = mvp[0][3] * x + mvp[1][3] * y + mvp[2][3] * z + mvp[3][3];

                vec2 c;
                c.x = (p_x / p_w * 0.5f + 0.5f) * res.x;
                c.y = (-p_y / p_w * 0.5f + 0.5f) * res.y;

                if (min_p.x <= c.x && c.x <= max_p.x && min_p.y <= c.y && c.y <= max_p.y) {
                    bitfield::set_bit(mask, i);
                }
            }

            switch (data->selection.level_mode) {
                case SelectionLevel::Atom:
                    break;
                case SelectionLevel::Residue:
                    expand_mask(mask, data->dynamic.molecule.residues);
                    break;
                case SelectionLevel::Chain:
                    expand_mask(mask, data->dynamic.molecule.chains);
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
            data->gpu_buffers.dirty.selection = true;

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

// static void handle_animation(ApplicationData* data) {}

// #camera-control
static void handle_camera_interaction(ApplicationData* data) {
    if (!ImGui::GetIO().WantCaptureMouse && !data->selection.selecting) {
        const vec2 mouse_delta = vec2(data->ctx.input.mouse.win_delta.x, data->ctx.input.mouse.win_delta.y);
        data->view.trackball_state.input.rotate_button = data->ctx.input.mouse.down[0];
        data->view.trackball_state.input.pan_button = data->ctx.input.mouse.down[1];
        data->view.trackball_state.input.dolly_button = data->ctx.input.mouse.down[2];
        data->view.trackball_state.input.mouse_coord_curr = {data->ctx.input.mouse.win_coord.x, data->ctx.input.mouse.win_coord.y};
        data->view.trackball_state.input.mouse_coord_prev = data->view.trackball_state.input.mouse_coord_curr - mouse_delta;
        data->view.trackball_state.input.screen_size = vec2(data->ctx.window.width, data->ctx.window.height);
        data->view.trackball_state.input.dolly_delta = data->ctx.input.mouse.scroll_delta;
        data->view.trackball_state.input.fov_y = data->view.camera.fov_y;

        {
            vec3 pos = data->view.animation.target_position;
            quat ori = data->view.camera.orientation;
            if (camera_controller_trackball(&pos, &ori, &data->view.trackball_state,
                                            TrackballFlags_RotateReturnsTrue | TrackballFlags_PanReturnsTrue)) {
                data->view.camera.position = pos;
                data->view.camera.orientation = ori;
            }
            data->view.animation.target_position = pos;
        }

        if (ImGui::GetIO().MouseDoubleClicked[0]) {
            if (data->picking.depth < 1.0f) {
                const vec3 forward = data->view.camera.orientation * vec3(0, 0, 1);
                const f32 dist = data->view.trackball_state.distance;
                data->view.animation.target_position = data->picking.world_coord + forward * dist;
            }
        }

        if (data->picking.depth < 1.0f) {
            data->visuals.dof.focus_depth.target = math::distance(data->view.camera.position, data->picking.world_coord);
        } else {
            data->visuals.dof.focus_depth.target = data->view.trackball_state.distance;
        }
    }
}

static void handle_camera_animation(ApplicationData* data) {
    const f32 dt = math::min(data->ctx.timing.delta_s, 0.033f);
    {
        // #camera-translation
        constexpr f32 speed = 10.0f;
        const vec3 vel = (data->view.animation.target_position - data->view.camera.position) * speed;
        data->view.camera.position += vel * dt;
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

static void compute_backbone_spline(ApplicationData* data) {
    PUSH_GPU_SECTION("Compute Backbone Spline") {
        bool has_spline_rep = false;
        for (const auto& rep : data->representations.buffer) {
            if (rep.type == RepresentationType::Ribbons || rep.type == RepresentationType::Cartoon) {
                has_spline_rep = true;
                break;
            }
        }
        has_spline_rep |= data->visuals.spline.draw_control_points || data->visuals.spline.draw_spline;

        data->gpu_buffers.dirty.backbone = true;
        if (has_spline_rep && data->gpu_buffers.dirty.backbone) {
            data->gpu_buffers.dirty.backbone = false;
            draw::compute_backbone_control_points(data->gpu_buffers.backbone.control_point, data->gpu_buffers.position,
                                                  data->gpu_buffers.backbone.backbone_segment_index,
                                                  data->gpu_buffers.backbone.num_backbone_segment_indices, ramachandran::get_segmentation_texture());
            draw::compute_backbone_spline(data->gpu_buffers.backbone.spline, data->gpu_buffers.backbone.control_point,
                                          data->gpu_buffers.backbone.control_point_index, data->gpu_buffers.backbone.num_control_point_indices);
        }
    }
    POP_GPU_SECTION()
}

static void compute_velocity(ApplicationData* data) {
    PUSH_GPU_SECTION("Compute View Velocity")
    draw::compute_pbc_view_velocity(data->gpu_buffers.velocity, data->gpu_buffers.position, data->gpu_buffers.old_position,
                                    (int)data->dynamic.molecule.atom.count, data->view.param, data->simulation_box.box * vec3(1, 1, 1));
    POP_GPU_SECTION()
}

static void compute_residue_aabbs(ApplicationData* data) {
    PUSH_GPU_SECTION("Compute Residue AABBs")
    draw::culling::compute_residue_aabbs(data->gpu_buffers.experimental.aabb, data->gpu_buffers.position, data->gpu_buffers.radius,
                                         data->gpu_buffers.experimental.residue, (i32)data->dynamic.molecule.residues.size());
    POP_GPU_SECTION()
}

static void cull_residue_aabbs(ApplicationData* data) {
    PUSH_GPU_SECTION("Cull Residue AABBs")
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glDepthMask(0);

    draw::culling::cull_aabbs(data->gpu_buffers.experimental.visibility, data->gpu_buffers.experimental.aabb, data->view.param,
                              (i32)data->dynamic.molecule.residues.size());

    glDepthMask(1);
    glDisable(GL_DEPTH_TEST);
    POP_GPU_SECTION()
}

static void draw_residue_aabbs(const ApplicationData& data) {
    PUSH_GPU_SECTION("Draw Residue AABBs")

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    // draw::culling::draw_aabbs(data.gpu_buffers.experimental.aabb, data.view.param, data.dynamic.molecule.residues.size());
    draw::culling::draw_culled_aabbs(data.gpu_buffers.experimental.visibility, data.gpu_buffers.experimental.aabb, data.view.param,
                                     (i32)data.dynamic.molecule.residues.size());

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    POP_GPU_SECTION()
}

static void init_occupancy_volume(ApplicationData* data) {

    const AABB box = compute_aabb(data->dynamic.molecule.atom.position.x, data->dynamic.molecule.atom.position.y, data->dynamic.molecule.atom.position.z,
                                  data->dynamic.molecule.atom.radius, data->dynamic.molecule.atom.count);
    cone_trace::init_occlusion_volume(&data->occupancy_volume.vol, box.min, box.max);

    data->occupancy_volume.vol_mutex.lock();
    if (data->occupancy_volume.vol.texture_id) {
        cone_trace::compute_occupancy_volume(data->occupancy_volume.vol, data->dynamic.molecule.atom.position.x,
                                             data->dynamic.molecule.atom.position.y, data->dynamic.molecule.atom.position.z,
                                             data->dynamic.molecule.atom.radius, data->representations.atom_visibility_mask);
    }
    data->occupancy_volume.vol_mutex.unlock();
}

static void fill_gbuffer(const ApplicationData& data) {
    const vec4 CLEAR_INDEX = vec4(1, 1, 1, 1);
    const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2,
                                   GL_COLOR_ATTACHMENT3, GL_COLOR_ATTACHMENT4, GL_COLOR_ATTACHMENT5};

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, data.fbo.deferred.fbo);
    glViewport(0, 0, data.fbo.width, data.fbo.height);

    glDepthMask(1);
    glColorMask(1, 1, 1, 1);

    // Setup fbo and clear textures
    PUSH_GPU_SECTION("Clear G-buffer") {
        // Clear color+alpha, normal, velocity, emissive, post_tonemap and depth
        glDrawBuffers(5, draw_buffers);
        glClearColor(0, 0, 0, 0);
        glClearDepthf(1.f);
        glStencilMask(0xFF);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

        // Clear picking buffer
        glDrawBuffer(GL_COLOR_ATTACHMENT5);
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
    if (data.simulation_box.enabled && data.dynamic.trajectory) {
        PUSH_GPU_SECTION("Draw Simulation Box")

        immediate::set_model_view_matrix(data.view.param.matrix.current.view);
        immediate::set_proj_matrix(data.view.param.matrix.current.proj_jittered);

        const mat3 box = get_trajectory_frame(data.dynamic.trajectory, data.playback.frame).box;
        const vec3 min_box = box * vec3(0.0f);
        const vec3 max_box = box * vec3(1.0f);

        // Simulation domain
        immediate::draw_box_wireframe(min_box, max_box, math::convert_color(data.simulation_box.color));

        for (const auto& ref : data.reference_frame.frames) {
            if (ref.show) {
                immediate::draw_box_wireframe(min_box, max_box, ref.reference_to_world, immediate::COLOR_GREEN);
            }
        }

        immediate::flush();

        POP_GPU_SECTION()
    }

    // RENDER DEBUG INFORMATION (WITH DEPTH)
    PUSH_GPU_SECTION("Debug Draw") {
        glDrawBuffer(GL_COLOR_ATTACHMENT4);  // Post_Tonemap buffer

        // draw_residue_aabbs(data);

        immediate::set_model_view_matrix(data.view.param.matrix.current.view);
        immediate::set_proj_matrix(data.view.param.matrix.current.proj);

        // HYDROGEN BONDS
        if (data.hydrogen_bonds.enabled && !data.hydrogen_bonds.overlay) {
            const auto& mol = data.dynamic.molecule;
            for (const auto& bond : data.hydrogen_bonds.bonds) {
                const vec3 pos0 = get_position_xyz(mol, bond.acc_idx);
                const vec3 pos1 = get_position_xyz(mol, bond.hyd_idx);
                immediate::draw_line(pos0, pos1, math::convert_color(data.hydrogen_bonds.color));
            }
        }

        immediate::flush();
    }
    POP_GPU_SECTION()

    PUSH_GPU_SECTION("Debug Draw Overlay") {
        glDrawBuffer(GL_COLOR_ATTACHMENT4);  // Post_Tonemap buffer
        glDisable(GL_DEPTH_TEST);
        glDepthMask(0);

        immediate::set_model_view_matrix(data.view.param.matrix.current.view);
        immediate::set_proj_matrix(data.view.param.matrix.current.proj);
        stats::visualize(data.dynamic);

        // HYDROGEN BONDS
        if (data.hydrogen_bonds.enabled && data.hydrogen_bonds.overlay) {
            for (const auto& bond : data.hydrogen_bonds.bonds) {
                const vec3 p0 = get_position_xyz(data.dynamic.molecule, bond.acc_idx);
                const vec3 p1 = get_position_xyz(data.dynamic.molecule, bond.hyd_idx);
                immediate::draw_line(p0, p1, math::convert_color(data.hydrogen_bonds.color));
            }
        }

        // REFERENCE FRAME
        for (const auto& ref : data.reference_frame.frames) {
            if (ref.show) {
                const mat4 B = ref.basis;

                const auto tracking_data = structure_tracking::get_tracking_data(ref.id);
                if (tracking_data) {
                    const mat3 pca = tracking_data->eigen.vectors[data.playback.frame];
                    immediate::draw_line(B[3], vec3(B[3]) + pca[0] * 2.0f, immediate::COLOR_RED);
                    immediate::draw_line(B[3], vec3(B[3]) + pca[1] * 2.0f, immediate::COLOR_GREEN);
                    immediate::draw_line(B[3], vec3(B[3]) + pca[2] * 2.0f, immediate::COLOR_BLUE);
                }

                // immediate::draw_basis(B, 4.0f);
                // immediate::draw_plane_wireframe(B[3], B[0], B[1], immediate::COLOR_BLACK);
            }
        }

        immediate::flush();

        PUSH_GPU_SECTION("Draw Control Points")
        if (data.visuals.spline.draw_control_points) {
            draw::draw_spline(data.gpu_buffers.backbone.control_point, data.gpu_buffers.backbone.control_point_index,
                              data.gpu_buffers.backbone.num_control_point_indices, data.view.param);
        }
        if (data.visuals.spline.draw_spline) {
            draw::draw_spline(data.gpu_buffers.backbone.spline, data.gpu_buffers.backbone.spline_index, data.gpu_buffers.backbone.num_spline_indices,
                              data.view.param);
        }
        POP_GPU_SECTION()

        glEnable(GL_DEPTH_TEST);
        glDepthMask(1);
    }
    POP_GPU_SECTION()

    // DRAW VELOCITY OF STATIC OBJECTS
    PUSH_GPU_SECTION("Blit Static Velocity")
    glDrawBuffer(GL_COLOR_ATTACHMENT2);
    glDepthMask(0);
    postprocessing::blit_static_velocity(data.fbo.deferred.depth, data.view.param);
    glDepthMask(1);
    POP_GPU_SECTION()

    // DRAW REPRESENTATIONS
    glDrawBuffers(ARRAY_SIZE(draw_buffers), draw_buffers);
    draw_representations(data);

    // VOLUME RENDERING
    if (data.density_volume.enabled) {
        PUSH_GPU_SECTION("Volume Rendering")

        glDrawBuffer(GL_COLOR_ATTACHMENT4);  // Post_Tonemap buffer

        glDepthMask(0);
        glDisable(GL_DEPTH_TEST);

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        volume::VolumeRenderDesc desc;
        const vec3 min_aabb = data.simulation_box.box * vec3(0, 0, 0);
        const vec3 max_aabb = data.simulation_box.box * vec3(1, 1, 1);
        desc.matrix.model = volume::compute_model_to_world_matrix(min_aabb, max_aabb);
        desc.matrix.view = data.view.param.matrix.current.view;
        desc.matrix.proj = data.view.param.matrix.current.proj_jittered;
        desc.texture.depth = data.fbo.deferred.depth;
        desc.texture.volume = data.density_volume.texture.id;
        desc.texture.transfer_function = data.density_volume.dvr.tf.id;
        desc.global_scaling.density = data.density_volume.dvr.density_scale;
        desc.global_scaling.alpha = data.density_volume.dvr.tf.alpha_scale;
        desc.isosurface = data.density_volume.iso.isosurfaces;
        desc.isosurface_enabled = data.density_volume.iso.enabled;
        desc.direct_volume_rendering_enabled = data.density_volume.dvr.enabled;
        desc.clip_volume.min = data.density_volume.clip_planes.min;
        desc.clip_volume.max = data.density_volume.clip_planes.max;
        desc.voxel_spacing = data.density_volume.voxel_spacing;

        volume::render_volume_texture(desc);

        glDisable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
        glDepthMask(1);

        POP_GPU_SECTION()
    }

    {
        PUSH_GPU_SECTION("Selection")
        const bool atom_selection_empty = !bitfield::any_bit_set(data.selection.current_selection_mask);
        const bool atom_highlight_empty = !bitfield::any_bit_set(data.selection.current_highlight_mask);

        glDrawBuffer(GL_COLOR_ATTACHMENT4);  // Post_Tonemap buffer
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_STENCIL_TEST);

        glDepthMask(0);
        glColorMask(0, 0, 0, 0);

        if (!atom_selection_empty) {
            glStencilOp(GL_KEEP, GL_REPLACE, GL_REPLACE);

            glStencilMask(0x2);
            glStencilFunc(GL_ALWAYS, 2, 0xFF);

            const vec4 color = data.selection.color.selection.fill_color;
            draw_representations_lean_and_mean(data, color, 1.0f, AtomBit_Selected);
        }

        if (!atom_highlight_empty) {
            glStencilOp(GL_KEEP, GL_REPLACE, GL_REPLACE);

            glStencilMask(0x4);
            glStencilFunc(GL_ALWAYS, 4, 0xFF);

            const vec4 color = data.selection.color.highlight.fill_color;
            draw_representations_lean_and_mean(data, color, 1.0f, AtomBit_Highlighted);
        }

        if (!atom_selection_empty || !atom_highlight_empty) {
            glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
            glDepthFunc(GL_LEQUAL);

            glStencilMask(0x1);
            glStencilFunc(GL_ALWAYS, 1, 0xFF);

            draw_representations_lean_and_mean(data);
        }

        glDisable(GL_DEPTH_TEST);

        glStencilMask(0x00);
        glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
        glColorMask(1, 1, 1, 1);

        if (!atom_selection_empty) {
            // Selection
            glStencilFunc(GL_NOTEQUAL, 2, 0x2);

            const vec4 color = data.selection.color.selection.outline_color;
            const float scale = data.selection.color.selection.outline_scale;
            draw_representations_lean_and_mean(data, color, scale, AtomBit_Selected);
        }

        if (!atom_highlight_empty) {
            // Highlight
            glStencilFunc(GL_NOTEQUAL, 4, 0x4);

            const vec4 color = data.selection.color.highlight.outline_color;
            const float scale = data.selection.color.highlight.outline_scale;
            draw_representations_lean_and_mean(data, color, scale, AtomBit_Highlighted);
        }

        if (!atom_selection_empty) {
            glStencilFunc(GL_NOTEQUAL, 1, 0x1);
            const float saturation = data.selection.color.selection_saturation;
            glDrawBuffer(GL_COLOR_ATTACHMENT0);
            postprocessing::scale_hsv(data.fbo.deferred.color, vec3(1, saturation, 1));
        }

        glDisable(GL_STENCIL_TEST);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        glDepthMask(1);
        POP_GPU_SECTION()
    }
    POP_GPU_SECTION()  // G-buffer

    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
}

static void update_properties(ApplicationData* data) {
    // @NOTE: This is stupid and should be updated
    stats::async_update(
        data->dynamic, {(i32)data->time_filter.range.beg, (i32)data->time_filter.range.end},
        [](void* usr_data) {
            ApplicationData* data = (ApplicationData*)usr_data;
            if (data->density_volume.enabled) {
                data->density_volume.volume_data_mutex.lock();
                const Range<i32> range = {(i32)data->time_filter.range.beg, (i32)data->time_filter.range.end};

                const ReferenceFrame* ref = get_active_reference_frame(data);
                if (ref) {
                    const structure_tracking::ID id = ref->id;
                    stats::compute_density_volume_with_basis(
                        &data->density_volume.volume, data->dynamic.trajectory, range, [ref, data](const vec4& world_pos, i32 frame_idx) -> vec4 {
                            const auto tracking_data = structure_tracking::get_tracking_data(ref->id);
                            if (!tracking_data) {
                                LOG_ERROR("Could not find tracking data for supplied id: '%lu'", ref->id);
                                return {};
                            }
                            const vec3* com_data = tracking_data->transform.com;
                            const quat* rot_data = nullptr;
                            switch (ref->tracking_mode) {
                                case TrackingMode::Absolute:
                                    rot_data = tracking_data->transform.rotation.absolute;
                                    break;
                                case TrackingMode::Relative:
                                    rot_data = tracking_data->transform.rotation.relative;
                                    break;
                                default:
                                    ASSERT(false);
                                    break;
                            }

                            if (com_data && rot_data) {
                                const mat3 box = data->dynamic.trajectory.frame_buffer[frame_idx].box;
                                const vec3 com = com_data[frame_idx];
                                const vec3 half_box = box * vec3(0.5f);
                                const mat4 R = math::mat4_cast(rot_data[frame_idx]);
                                const mat4 T_com = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(-com, 1));
                                const mat4 T_box = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(half_box, 1));
                                const mat4 M = T_box * math::transpose(R) * T_com;
                                const mat4 world_to_model = volume::compute_world_to_model_matrix(box * vec3(0, 0, 0), box * vec3(1, 1, 1));
                                return world_to_model * M * world_pos;
                            } else {
                                LOG_ERROR("Could not get com data or rotation data");
                                return {};
                            }
                        });
                } else {
                    const vec3 min_aabb = data->simulation_box.box * vec3(0, 0, 0);
                    const vec3 max_aabb = data->simulation_box.box * vec3(1, 1, 1);
                    const mat4 world_to_model = volume::compute_world_to_model_matrix(min_aabb, max_aabb);
                    stats::compute_density_volume(&data->density_volume.volume, data->dynamic.trajectory, range, world_to_model);
                }
                data->density_volume.volume_data_mutex.unlock();
                data->density_volume.texture.dirty = true;
            }
        },
        data);
}

static void update_density_volume_texture(ApplicationData* data) {
    // If gpu representation of volume is not up to date, upload data
    if (data->density_volume.texture.dirty) {
        // @NOTE: We need to lock the data here in case the data is changing while doing this.
        if (data->density_volume.volume_data_mutex.try_lock()) {
            volume::create_volume_texture(&data->density_volume.texture.id, data->density_volume.volume.dim);
            volume::set_volume_texture_data(data->density_volume.texture.id, data->density_volume.volume.dim,
                                            data->density_volume.volume.voxel_data.ptr, data->density_volume.volume.voxel_range.max);
            data->density_volume.volume_data_mutex.unlock();
            data->density_volume.texture.max_value = (float)data->density_volume.volume.voxel_range.end;
            data->density_volume.texture.dirty = false;
        }
    }

    if (data->density_volume.dvr.tf.dirty) {
        if (data->density_volume.dvr.tf.path.length() > 0) {
            volume::create_tf_texture(&data->density_volume.dvr.tf.id, &data->density_volume.dvr.tf.width, data->density_volume.dvr.tf.path);
        } else {
            volume::create_tf_texture(&data->density_volume.dvr.tf.id, &data->density_volume.dvr.tf.width, VIAMD_IMAGE_DIR "/tf/default.png");
        }
        data->density_volume.dvr.tf.dirty = false;
    }
}

static void handle_picking(ApplicationData* data) {
    PUSH_CPU_SECTION("PICKING") {
        const auto ndc = data->ctx.input.mouse.ndc_coord;
        vec2 coord = vec2(ndc.x, ndc.y) * 0.5f + 0.5f;
        coord *= vec2(data->fbo.width - 1, data->fbo.height - 1);
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
                    data->picking.idx = math::clamp(data->picking.idx, 0U, (u32)data->dynamic.molecule.atom.count - 1U);
                const vec4 viewport(0, 0, data->fbo.width, data->fbo.height);
                data->picking.world_coord =
                    math::unproject(vec3(coord.x, coord.y, data->picking.depth), data->view.param.matrix.inverse.view_proj_jittered, viewport);
            }
#else
            coord += data->view.param.jitter.current;
            // coord -= 0.5f;
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
    const float dt_compensation = MOTION_BLUR_REFERENCE_DT / data.ctx.timing.delta_s;
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
    desc.input_textures.emissive = data.fbo.deferred.emissive;
    desc.input_textures.post_tonemap = data.fbo.deferred.post_tonemap;

    postprocessing::shade_and_postprocess(desc, data.view.param);
    POP_GPU_SECTION()
}

// #reference-frame
static ReferenceFrame* create_reference_frame(ApplicationData* data, CStringView name, CStringView filter) {
    ASSERT(data);
    ReferenceFrame ref;
    ref.id = structure_tracking::create_structure();
    ref.name = name;
    ref.filter = filter;
    bitfield::init(&ref.atom_mask, data->dynamic.molecule.atom.count);
    compute_reference_frame(data, &ref);
    return &data->reference_frame.frames.push_back(ref);
}

static bool remove_reference_frame(ApplicationData* data, ReferenceFrame* ref) {
    ASSERT(data);
    ASSERT(ref);
    bitfield::free(&ref->atom_mask);
    data->reference_frame.frames.remove(ref);
    return true;
}

#if 0
static void reset_reference_frame(ApplicationData* data) {
    ASSERT(data);
    UNUSED(data);
    // @NOTE: What to do here?
}
#endif

static void clear_reference_frames(ApplicationData* data) {
    ASSERT(data);
    while (data->reference_frame.frames.size() > 0) {
        remove_reference_frame(data, &data->reference_frame.frames.back());
    }
}

static void compute_reference_frame(ApplicationData* data, ReferenceFrame* ref) {
    ASSERT(data);
    ASSERT(ref);
    const auto& mol = data->dynamic.molecule;

    if (ref->atom_mask.size() != mol.atom.count) {
        bitfield::init(&ref->atom_mask, mol.atom.count);
    }

    DynamicArray<StoredSelection> sel;
    for (const auto& s : data->selection.stored_selections) {
        sel.push_back({s.name, s.atom_mask});
    }
    sel.push_back({"current", data->selection.current_selection_mask});

    ref->filter_is_ok = filter::compute_filter_mask(ref->atom_mask, ref->filter, data->dynamic.molecule, sel);
    if (ref->filter_is_ok) {
        structure_tracking::compute_trajectory_transform_data(ref->id, data->dynamic, ref->atom_mask);
    }
}

static ReferenceFrame* get_active_reference_frame(ApplicationData* data) {
    ASSERT(data);
    for (auto& ref : data->reference_frame.frames) {
        if (ref.active) return &ref;
    }
    return nullptr;
}

static void update_reference_frames(ApplicationData* data) {
    ASSERT(data);

    const auto& mol = data->dynamic.molecule;
    const auto& traj = data->dynamic.trajectory;

    const f64 time = data->playback.time;
    const int last_frame = math::max(0, traj.num_frames - 1);
    const f32 t = (float)math::fract(data->playback.time);
    const int frame = math::clamp((int)time, 0, last_frame);
    const int p2 = math::max(0, frame - 1);
    const int p1 = math::max(0, frame);
    const int n1 = math::min(frame + 1, last_frame);
    const int n2 = math::min(frame + 2, last_frame);

    const InterpolationMode mode = (p1 != n1) ? data->playback.interpolation : InterpolationMode::Nearest;

    for (auto& ref : data->reference_frame.frames) {
        const structure_tracking::TrackingData* tracking_data = structure_tracking::get_tracking_data(ref.id);
        if (tracking_data) {
            const i64 masked_count = bitfield::number_of_bits_set(ref.atom_mask);
            if (masked_count == 0) break;

            void* w_mem = TMP_MALLOC(sizeof(float) * 1 * masked_count);
            void* f_mem = TMP_MALLOC(sizeof(float) * 3 * masked_count);
            void* c_mem = TMP_MALLOC(sizeof(float) * 3 * masked_count);
            defer {
                TMP_FREE(w_mem);
                TMP_FREE(f_mem);
                TMP_FREE(c_mem);
            };

            float* weight = (float*)w_mem;
            float* frame_x = (float*)f_mem;
            float* frame_y = (float*)f_mem + 1 * masked_count;
            float* frame_z = (float*)f_mem + 2 * masked_count;
            float* current_x = (float*)c_mem;
            float* current_y = (float*)c_mem + 1 * masked_count;
            float* current_z = (float*)c_mem + 2 * masked_count;

            bitfield::gather_masked(weight, mol.atom.mass, ref.atom_mask);
            bitfield::gather_masked(frame_x, get_trajectory_position_x(traj, 0).data(), ref.atom_mask);
            bitfield::gather_masked(frame_y, get_trajectory_position_y(traj, 0).data(), ref.atom_mask);
            bitfield::gather_masked(frame_z, get_trajectory_position_z(traj, 0).data(), ref.atom_mask);
            bitfield::gather_masked(current_x, mol.atom.position.x, ref.atom_mask);
            bitfield::gather_masked(current_y, mol.atom.position.y, ref.atom_mask);
            bitfield::gather_masked(current_z, mol.atom.position.z, ref.atom_mask);

            const vec3 frame_com = compute_com(frame_x, frame_y, frame_z, weight, masked_count);
            const vec3 current_com = compute_com(current_x, current_y, current_z, weight, masked_count);

            const quat* rot_data = nullptr;
            switch (ref.tracking_mode) {
                case TrackingMode::Absolute:
                    rot_data = tracking_data->transform.rotation.absolute;
                    break;
                case TrackingMode::Relative:
                    rot_data = tracking_data->transform.rotation.relative;
                    break;
                default:
                    ASSERT(false);
                    break;
            }

            // clang-format off
            const quat q[4] = { rot_data[p2],
                                rot_data[p1],
                                rot_data[n1],
                                rot_data[n2] };

            const vec3 box[4] = { get_trajectory_frame(traj, p2).box * vec3(0.5f),
                                  get_trajectory_frame(traj, p1).box * vec3(0.5f),
                                  get_trajectory_frame(traj, n1).box * vec3(0.5f),
                                  get_trajectory_frame(traj, n2).box * vec3(0.5f) };
            // clang-format on

            vec3 box_c;
            mat4 R;
            switch (mode) {
                case InterpolationMode::Nearest:
                    R = math::mat4_cast(t < 0.5f ? q[1] : q[2]);
                    box_c = t < 0.5f ? box[1] : box[2];
                    break;
                case InterpolationMode::Linear:
                    R = math::mat4_cast(math::slerp(q[1], q[2], t));
                    box_c = math::lerp(box[1], box[2], t);
                    break;
                case InterpolationMode::Cubic:
                    R = math::mat4_cast(math::cubic_slerp(q[0], q[1], q[2], q[3], t));
                    box_c = math::cubic_spline(box[0], box[1], box[2], box[3], t);
                    break;
                default:
                    ASSERT(false);
            }

            const mat4 T_com = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(-current_com, 1));
            const mat4 T_box = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(box_c, 1));
            const mat4 R_inv = math::transpose(R);

            const mat4 PCA = tracking_data->simulation_box_aligned_pca;

            ref.world_to_reference = T_box * PCA * R * T_com;
            ref.reference_to_world = math::inverse(ref.world_to_reference);

            // if (ref.active) {
            //    ref.basis = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {box_c, 1}};
            //} else {
            ref.basis = math::inverse(PCA * R * T_com);
            ref.pca = tracking_data->eigen.vectors[p1];
            //}
        } else {
            LOG_ERROR("Could not find tracking data for '%lu'", ref.id);
        }
    }
}

static void superimpose_ensemble(ApplicationData* data) {
    ASSERT(data);

    const auto& mol = data->dynamic.molecule;
    const auto& traj = data->dynamic.trajectory;

    const f64 time = data->playback.time;
    const int last_frame = math::max(0, traj.num_frames - 1);
    const f32 t = (float)math::fract(data->playback.time);
    const int frame = math::clamp((int)time, 0, last_frame);
    const int p2 = math::max(0, frame - 1);
    const int p1 = math::max(0, frame);
    const int n1 = math::min(frame + 1, last_frame);
    const int n2 = math::min(frame + 2, last_frame);

    const InterpolationMode mode = (p1 != n1) ? data->playback.interpolation : InterpolationMode::Nearest;
    const Bitfield atom_mask = data->ensemble_tracking.atom_mask;

    Bitfield ensemble_atom_mask;
    bitfield::init(&ensemble_atom_mask, mol.atom.count);
    defer { bitfield::free(&ensemble_atom_mask); };
    bitfield::clear_all(ensemble_atom_mask);

    for (auto& structure : data->ensemble_tracking.structures) {
        const int idx = (int)(&structure - data->ensemble_tracking.structures.beg());

        const structure_tracking::TrackingData* tracking_data = structure_tracking::get_tracking_data(structure.id);
        if (tracking_data) {

            const i64 masked_count = bitfield::number_of_bits_set(atom_mask);
            if (masked_count == 0) break;

            void* w_mem = TMP_MALLOC(sizeof(float) * 1 * masked_count);
            void* f_mem = TMP_MALLOC(sizeof(float) * 3 * masked_count);
            void* c_mem = TMP_MALLOC(sizeof(float) * 3 * masked_count);
            defer {
                TMP_FREE(w_mem);
                TMP_FREE(f_mem);
                TMP_FREE(c_mem);
            };

            float* weight = (float*)w_mem;
            float* frame_x = (float*)f_mem;
            float* frame_y = (float*)f_mem + 1 * masked_count;
            float* frame_z = (float*)f_mem + 2 * masked_count;
            float* current_x = (float*)c_mem;
            float* current_y = (float*)c_mem + 1 * masked_count;
            float* current_z = (float*)c_mem + 2 * masked_count;

            bitfield::gather_masked(weight, mol.atom.mass, atom_mask, structure.offset);
            bitfield::gather_masked(frame_x, get_trajectory_position_x(traj, 0).data(), atom_mask, structure.offset);
            bitfield::gather_masked(frame_y, get_trajectory_position_y(traj, 0).data(), atom_mask, structure.offset);
            bitfield::gather_masked(frame_z, get_trajectory_position_z(traj, 0).data(), atom_mask, structure.offset);
            bitfield::gather_masked(current_x, mol.atom.position.x, atom_mask, structure.offset);
            bitfield::gather_masked(current_y, mol.atom.position.y, atom_mask, structure.offset);
            bitfield::gather_masked(current_z, mol.atom.position.z, atom_mask, structure.offset);

            const vec3 frame_com = compute_com(frame_x, frame_y, frame_z, weight, masked_count);
            const vec3 current_com = compute_com(current_x, current_y, current_z, weight, masked_count);

            const quat* rot_data = nullptr;
            switch (data->ensemble_tracking.tracking_mode) {
                case TrackingMode::Absolute:
                    rot_data = tracking_data->transform.rotation.absolute;
                    break;
                case TrackingMode::Relative:
                    rot_data = tracking_data->transform.rotation.relative;
                    break;
                default:
                    ASSERT(false);
                    break;
            }

            // clang-format off
            const quat q[4] = { rot_data[p2],
                                rot_data[p1],
                                rot_data[n1],
                                rot_data[n2] };

            const vec3 box[4] = { get_trajectory_frame(traj, p2).box * vec3(0.5f),
                                  get_trajectory_frame(traj, p1).box * vec3(0.5f),
                                  get_trajectory_frame(traj, n1).box * vec3(0.5f),
                                  get_trajectory_frame(traj, n2).box * vec3(0.5f) };
            // clang-format on

            vec3 box_c;
            mat4 R;
            switch (mode) {
                case InterpolationMode::Nearest:
                    R = math::mat4_cast(t < 0.5f ? q[1] : q[2]);
                    box_c = t < 0.5f ? box[1] : box[2];
                    break;
                case InterpolationMode::Linear:
                    R = math::mat4_cast(math::slerp(q[1], q[2], t));
                    box_c = math::lerp(box[1], box[2], t);
                    break;
                case InterpolationMode::Cubic:
                    R = math::mat4_cast(math::cubic_slerp(q[0], q[1], q[2], q[3], t));
                    box_c = math::cubic_spline(box[0], box[1], box[2], box[3], t);
                    break;
                default:
                    ASSERT(false);
            }

            const mat4 T_com = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(-current_com, 1));
            const mat4 T_box = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(box_c, 1));

            // mat4 PCA = tracking_data->eigen.vector[0];
            // PCA = mat4(1);

            structure.world_to_structure = T_box * structure.alignment_matrix * R * T_com;
            for (i64 i = 0; i < atom_mask.size(); i++) {
                if (atom_mask[i]) {
                    bitfield::set_bit(ensemble_atom_mask, structure.offset + i);
                    float& x = mol.atom.position.x[structure.offset + i];
                    float& y = mol.atom.position.y[structure.offset + i];
                    float& z = mol.atom.position.z[structure.offset + i];
                    vec4 v = {x, y, z, 1.0f};
                    v = structure.world_to_structure * v;
                    x = v.x;
                    y = v.y;
                    z = v.z;
                }
            }
        } else {
            LOG_ERROR("Could not find tracking data for '%lu'", structure.id);
        }
    }

    // @NOTE: Modify all non ensemble included atoms with the reference structures transform
    for (i64 i = 0; i < ensemble_atom_mask.size(); i++) {
        if (!ensemble_atom_mask[i]) {
            float& x = mol.atom.position.x[i];
            float& y = mol.atom.position.y[i];
            float& z = mol.atom.position.z[i];
            vec4 v = {x, y, z, 1.0f};
            v = data->ensemble_tracking.structures[0].world_to_structure * v;
            x = v.x;
            y = v.y;
            z = v.z;
        }
    }
}

static void on_trajectory_load_complete(ApplicationData* data) {
    for (auto& ref : data->reference_frame.frames) {
        compute_reference_frame(data, &ref);
    }
    update_reference_frames(data);
    stats::set_all_property_flags(true, true);
    ramachandran::initialize(data->dynamic);
}

static void init_density_volume(ApplicationData* data) {
    data->density_volume.volume_data_mutex.lock();
    const vec3 min_box = data->simulation_box.box * vec3(0, 0, 0);
    const vec3 max_box = data->simulation_box.box * vec3(1, 1, 1);
    const ivec3 dim = math::max(ivec3(1), ivec3(max_box * data->density_volume.resolution_scale));
    init_volume(&data->density_volume.volume, dim);
    clear_volume(&data->density_volume.volume);
    data->density_volume.texture.dirty = true;
    data->density_volume.reference_structure = {};
    data->density_volume.voxel_spacing = (max_box - min_box) / vec3(data->density_volume.volume.dim);
    data->density_volume.volume_data_mutex.unlock();
}

static void draw_representations(const ApplicationData& data) {
    const i32 atom_count = (i32)data.dynamic.molecule.atom.count;
    const i32 bond_count = (i32)data.dynamic.molecule.covalent_bonds.count;
    const i32 res_count = (i32)data.dynamic.molecule.residues.count;

    PUSH_GPU_SECTION("Full Detail") for (const auto& rep : data.representations.buffer) {
        if (!rep.enabled) continue;
        switch (rep.type) {
            case RepresentationType::Vdw:
                PUSH_GPU_SECTION("Vdw")
#if EXPERIMENTAL_CULLING == 1
                draw::culling::draw_culled_vdw(data.gpu_buffers.position, data.gpu_buffers.radius, rep.color_buffer, data.gpu_buffers.velocity,
                                               data.gpu_buffers.experimental.visibility, data.gpu_buffers.experimental.residue, res_count,
                                               data.view.param, rep.radius);
#else
                draw::draw_vdw(data.gpu_buffers.position, data.gpu_buffers.radius, rep.color_buffer, data.gpu_buffers.velocity, atom_count,
                               data.view.param, rep.radius);
#endif
                POP_GPU_SECTION()
                break;
            case RepresentationType::Licorice:
                PUSH_GPU_SECTION("Licorice")
                draw::draw_licorice(data.gpu_buffers.position, rep.color_buffer, data.gpu_buffers.velocity, data.gpu_buffers.bond, bond_count,
                                    data.view.param, rep.radius);
                POP_GPU_SECTION()
                break;
            case RepresentationType::BallAndStick:
                PUSH_GPU_SECTION("Vdw")
                draw::draw_vdw(data.gpu_buffers.position, data.gpu_buffers.radius, rep.color_buffer, data.gpu_buffers.velocity, atom_count,
                               data.view.param, rep.radius * BALL_AND_STICK_VDW_SCALE);
                POP_GPU_SECTION()
                PUSH_GPU_SECTION("Licorice")
                draw::draw_licorice(data.gpu_buffers.position, rep.color_buffer, data.gpu_buffers.velocity, data.gpu_buffers.bond, bond_count,
                                    data.view.param, rep.radius * BALL_AND_STICK_LICORICE_SCALE);
                POP_GPU_SECTION()
                break;
            case RepresentationType::Ribbons:
                PUSH_GPU_SECTION("Ribbons")
                draw::draw_ribbons(data.gpu_buffers.backbone.spline, data.gpu_buffers.backbone.spline_index, rep.color_buffer,
                                   data.gpu_buffers.velocity, data.gpu_buffers.backbone.num_spline_indices, data.view.param);
                POP_GPU_SECTION()
                break;
            case RepresentationType::Cartoon:
                PUSH_GPU_SECTION("Cartoon")
                draw::draw_cartoon(data.gpu_buffers.backbone.spline, data.gpu_buffers.backbone.spline_index, rep.color_buffer,
                                   data.gpu_buffers.backbone.num_spline_indices, data.view.param);
                POP_GPU_SECTION()
                break;
        }
    }
    POP_GPU_SECTION()
}

static void draw_representations_lean_and_mean(const ApplicationData& data, vec4 color, float scale, u32 mask) {
    const i32 atom_count = (i32)data.dynamic.molecule.atom.count;
    const i32 bond_count = (i32)data.dynamic.molecule.covalent_bonds.size();

    PUSH_GPU_SECTION("Lean and Mean")
    for (const auto& rep : data.representations.buffer) {
        if (!rep.enabled) continue;
        // if (!rep.show_in_selection) continue;
        switch (rep.type) {
            case RepresentationType::Vdw:
                PUSH_GPU_SECTION("Vdw")
                draw::lean_and_mean::draw_vdw(data.gpu_buffers.position, data.gpu_buffers.radius, rep.color_buffer, data.gpu_buffers.selection,
                                              atom_count, data.view.param, rep.radius * scale, color, mask);
                POP_GPU_SECTION()
                break;
            case RepresentationType::Licorice:
                PUSH_GPU_SECTION("Licorice")
                draw::lean_and_mean::draw_licorice(data.gpu_buffers.position, rep.color_buffer, data.gpu_buffers.selection, data.gpu_buffers.bond,
                                                   bond_count, data.view.param, rep.radius * scale, color, mask);
                POP_GPU_SECTION()
                break;
            case RepresentationType::BallAndStick:
                PUSH_GPU_SECTION("Vdw")
                draw::lean_and_mean::draw_vdw(data.gpu_buffers.position, data.gpu_buffers.radius, rep.color_buffer, data.gpu_buffers.selection,
                                              atom_count, data.view.param, rep.radius * scale * BALL_AND_STICK_VDW_SCALE, color, mask);
                POP_GPU_SECTION()
                PUSH_GPU_SECTION("Licorice")
                draw::lean_and_mean::draw_licorice(data.gpu_buffers.position, rep.color_buffer, data.gpu_buffers.selection, data.gpu_buffers.bond,
                                                   bond_count, data.view.param, rep.radius * scale * BALL_AND_STICK_LICORICE_SCALE, color, mask);
                POP_GPU_SECTION()
                break;
            case RepresentationType::Ribbons:
                PUSH_GPU_SECTION("Ribbons")
                draw::lean_and_mean::draw_ribbons(data.gpu_buffers.backbone.spline, data.gpu_buffers.backbone.spline_index, rep.color_buffer,
                                                  data.gpu_buffers.selection, data.gpu_buffers.backbone.num_spline_indices, data.view.param, scale,
                                                  color, mask);
                POP_GPU_SECTION()
                break;
            default:
                break;
        }
    }
    POP_GPU_SECTION()
}
