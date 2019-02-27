#include <core/gl.h>
#include <core/types.h>
#include <core/hash.h>
#include <core/log.h>
#include <core/math_utils.h>
#include <core/camera.h>
#include <core/camera_utils.h>
#include <core/string_utils.h>
#include <core/volume.h>

#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>
#include <mol/trajectory_utils.h>
#include <mol/molecule_utils.h>
#include <mol/hydrogen_bond.h>
#include <mol/filter.h>
#include <mol/pdb_utils.h>
#include <mol/gro_utils.h>

#include <mol/spatial_hash.h>

#include <gfx/molecule_draw.h>
#include <gfx/immediate_draw_utils.h>
#include <gfx/postprocessing_utils.h>
#include <gfx/volume_utils.h>

#include <imgui.h>
#define IMGUI_DEFINE_MATH_OPERATORS
#include <imgui_internal.h>

#include <range_slider.h>
#include <plot_extended.h>

#include <glm/gtx/io.hpp>

#include <stdio.h>
#include <thread>
#include <atomic>
#include <mutex>
#include <iostream>

#include "platform/platform.h"
#include "console.h"
#include "stats.h"
#include "ramachandran.h"
#include "color_utils.h"

//#define VIAMD_RELEASE

#define ARRAY_SIZE(x) (sizeof(x) / sizeof(x[0]))

#ifdef OS_MAC_OSX
constexpr Key::Key_t CONSOLE_KEY = Key::KEY_WORLD_1;
#else  // WIN32 and Linux
// @TODO: Make sure this is currect for Linux?
constexpr Key::Key_t CONSOLE_KEY = Key::KEY_GRAVE_ACCENT;
#endif
constexpr Key::Key_t PLAY_PAUSE_KEY = Key::KEY_SPACE;

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
constexpr const char* FILE_EXTENSION = "via";

constexpr uint32 DEL_BTN_COLOR = 0xFF1111CC;
constexpr uint32 DEL_BTN_HOVER_COLOR = 0xFF3333DD;
constexpr uint32 DEL_BTN_ACTIVE_COLOR = 0xFF5555FF;
constexpr uint32 TEXT_BG_ERROR_COLOR = 0xAA222299;

constexpr float32 HYDROGEN_BOND_DISTANCE_CUTOFF_DEFAULT = 3.0f;
constexpr float32 HYDROGEN_BOND_DISTANCE_CUTOFF_MIN = 0.1f;
constexpr float32 HYDROGEN_BOND_DISTANCE_CUTOFF_MAX = 12.0f;

constexpr float32 HYDROGEN_BOND_ANGLE_CUTOFF_DEFAULT = 20.f;
constexpr float32 HYDROGEN_BOND_ANGLE_CUTOFF_MIN = 5.f;
constexpr float32 HYDROGEN_BOND_ANGLE_CUTOFF_MAX = 90.f;

constexpr float32 BALL_AND_STICK_VDW_SCALE = 0.25f;
constexpr float32 BALL_AND_STICK_LICORICE_SCALE = 0.5f;

constexpr int32 VOLUME_DOWNSAMPLE_FACTOR = 2;

constexpr int32 SPLINE_SUBDIVISION_COUNT = 16;

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

enum class PlaybackInterpolationMode { Nearest, Linear, LinearPbc, Cubic, CubicPbc };
enum class SelectionLevel { Atom, Residue, Chain };
enum class SelectionOperator { Or, And };
enum class SelectionGrowth { CovalentBond, Radial };
enum class RepresentationType { Vdw, Licorice, BallAndStick, Ribbons, Cartoon };

struct PickingData {
    uint32 idx = NO_PICKING_IDX;
    float32 depth = 1.0f;
    vec3 world_coord = {0, 0, 0};
};

struct MainFramebuffer {
    struct {
        GLuint depth = 0;
        GLuint color = 0;
        GLuint normal = 0;
        GLuint velocity = 0;
        GLuint emissive = 0;
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
    GLuint velocity = 0;
    GLuint radius = 0;
    GLuint selection = 0;
    GLuint bond = 0;

    struct {
        GLuint backbone_segment_index = 0;  // Stores indices to atoms which are needed for control points and support vectors
        GLuint control_point = 0;           // Stores extracted control points[vec3] and support vectors[vec3] before subdivision
        GLuint control_point_index = 0;     // Stores draw element indices to compute splines with adjacent index information (each chain separated by restart-index 0xFFFFFFFFU).
        GLuint spline = 0;                  // Stores subdivided spline data control points[vec3] + support vector[vec3]
        GLuint spline_index = 0;            // Stores draw element indices to render splines (each chain separated by restart-index 0xFFFFFFFFU).

        int32 num_backbone_segment_indices = 0;
        int32 num_control_point_indices = 0;
        int32 num_spline_indices = 0;
    } backbone;

    struct {
        bool position = false;
        bool velocity = false;
        // bool radius = false;
        bool selection = false;
        // bool bond = false;
        bool backbone = false;
    } dirty;
};

struct Representation {

    StringBuffer<128> name = "rep";
    StringBuffer<128> filter = "all";
    RepresentationType type = RepresentationType::Vdw;
    ColorMapping color_mapping = ColorMapping::Cpk;
    Array<bool> atom_mask{};
    GLuint color_buffer = 0;

    bool enabled = true;
    bool filter_is_ok = true;

    // For ColorMapping::StaticColor mode
    vec4 static_color = vec4(1);

    // VDW and Ball & Stick
    float32 radius = 1.f;

    // Ball & Stick and Licorice, Ribbons, Cartoon
    float32 thickness = 1.f;

    // Ribbons, Cartoon
    float32 tension = 0.5f;
    float32 width = 1.f;
};

struct Selection {
    StringBuffer<128> name = "sel";
    Array<bool> atom_mask{};
};

struct ThreadSyncData {
    std::thread thread{};
    std::atomic<bool> running{false};
    std::atomic<bool> stop_signal{false};

    void signal_stop() { stop_signal = true; }

    void wait_until_finished() {
        while (running) {
            platform::sleep(1);
        }
        // thread.join();
    }

    void signal_stop_and_wait() {
        signal_stop();
        wait_until_finished();
    }
};

struct ApplicationData {
    // --- PLATFORM ---
    platform::Context ctx;

    uint64 dirty_flag = 0;

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

        struct {
            vec3 target_position{};
            // quat target_orientation{};
        } animation;
    } view;

    // --- MOLECULAR DATA ---
    struct {
        MoleculeDynamic dynamic{};
        DynamicArray<float> atom_radii{};
        DynamicArray<vec3> atom_velocity{};
    } mol_data;

    // --- THREAD SYNCHRONIZATION ---
    struct {
        struct {
            ThreadSyncData sync{};
            float32 fraction = 0.f;
        } trajectory;

        struct {
            ThreadSyncData sync{};
            float32 fraction = 0.f;
            bool query_update = false;
        } backbone_angles;
    } async;

    // --- ATOM SELECTION ---
    struct {
        bool show_window = false;
        SelectionLevel level_mode = SelectionLevel::Atom;
        SelectionOperator op_mode = SelectionOperator::Or;
        SelectionGrowth grow_mode = SelectionGrowth::CovalentBond;

        int32 hovered = -1;
        int32 right_clicked = -1;

        DynamicArray<bool> current_selection_mask{};
        DynamicArray<bool> current_highlight_mask{};
        DynamicArray<Selection> stored_selections{};

        vec3 highlight_color = vec3(1, 1, 1);
        float highlight_scale = 2.0f;
        vec3 selection_color = vec3(1, 1, 0);
        float selection_scale = 2.0f;
        vec3 outline_color = vec3(1, 1, 1);
        float outline_scale = 10.0f;

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
    uint64 frame = 0;
    float64 time = 0.f;  // needs to be double precision for long trajectories
    float32 frames_per_second = 10.f;
    bool is_playing = false;
    PlaybackInterpolationMode interpolation = PlaybackInterpolationMode::CubicPbc;

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
            bool enabled = false;
            float32 intensity = 3.0f;
            float32 radius = 6.0f;
            float32 bias = 0.1f;
        } ssao;

        struct {
            bool enabled = false;
            float32 focus_depth = 0.5f;
            float32 focus_scale = 1.0f;
        } dof;

        struct {
            bool enabled = true;
            bool jitter = true;
            float32 feedback_min = 0.88f;
            float32 feedback_max = 0.97f;

            struct {
                bool enabled = true;
                float32 motion_scale = 0.5f;
            } motion_blur;
        } temporal_reprojection;

        struct {
            bool enabled = true;
            postprocessing::Tonemapping tonemapper = postprocessing::Tonemapping_Filmic;
            float32 exposure = 1.f;
            float32 gamma = 2.2;
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
        float32 distance_cutoff = HYDROGEN_BOND_DISTANCE_CUTOFF_DEFAULT;  // In Ångström
        float32 angle_cutoff = HYDROGEN_BOND_ANGLE_CUTOFF_DEFAULT;        // In Degrees
        DynamicArray<HydrogenBond> bonds{};
    } hydrogen_bonds;

    struct {
        bool enabled = false;
        vec4 color = vec4(0, 0, 0, 0.5);
    } simulation_box;

    struct {
        bool enabled = false;
        vec3 color = vec3(1, 0, 0);
        float32 density_scale = 1.f;

        struct {
            GLuint id = 0;
            bool dirty = false;
            ivec3 dim = ivec3(0);
            float32 max_value = 1.f;
        } texture;

        Volume volume{};
        std::mutex volume_data_mutex{};

        mat4 model_to_world_matrix{};
        mat4 texture_to_model_matrix{};
        mat4 world_to_texture_matrix{};
    } density_volume;

    // --- RAMACHANDRAN ---
    struct {
        bool show_window = false;
        int frame_range_min = 0;
        int frame_range_max = 0;

        ramachandran::ColorMap color_map{};

        struct {
            bool enabled = false;
            float32 radius = 0.2f;
            vec4 color = vec4(0, 0, 0, 1);
        } range;

        struct {
            bool enabled = true;
            float32 radius = 1.5f;
            vec4 border_color = vec4(0, 0, 0, 1);
            vec4 fill_color = vec4(1, 1, 1, 1);
        } current;

        struct {
            bool enabled = true;
            float32 radius = 2.5f;
            vec4 border_color = vec4(0, 0, 0, 1);
            vec4 fill_color = vec4(1, 1, 0, 1);
        } selected;

        BackboneAnglesTrajectory backbone_angles{};
    } ramachandran;

    // --- REPRESENTATIONS ---
    struct {
        DynamicArray<Representation> buffer{};
        DynamicArray<bool> atom_visibility_mask{};
        bool show_window = false;
        bool changed = false;
    } representations;

    /*
struct {
    spatialhash::Frame frame = {};
    vec3 cell_ext = vec3(4.0f);  // in Ångström
    bool dirty = true;
} spatial_hash;
    */

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

static bool IsItemActivePreviousFrame() {
    ImGuiContext& g = *GImGui;
    if (g.ActiveIdPreviousFrame) return g.ActiveIdPreviousFrame == GImGui->CurrentWindow->DC.LastItemId;
    return false;
}
}  // namespace ImGui

static void interpolate_atomic_positions(Array<vec3> dst_pos, const MoleculeTrajectory& traj, float64 time, PlaybackInterpolationMode interpolation_mode);
static void compute_atomic_velocities(Array<vec3> dst_vel, Array<const vec3> pos, Array<const vec3> old_pos, const vec3& box_ext = vec3(0));
static void reset_view(ApplicationData* data, bool move_camera = false, bool smooth_transition = false);
static float32 compute_avg_ms(float32 dt);
static PickingData read_picking_data(const MainFramebuffer& fbo, int32 x, int32 y);
static bool handle_selection(ApplicationData* data);
static void draw_representations_lean_and_mean(ApplicationData* data, vec4 color = vec4(1,1,1,1), float scale = 1.0f, uint32 mask = 0xFFFFFFFFU);

static void draw_main_menu(ApplicationData* data);
static void draw_context_popup(ApplicationData* data);
static void draw_animation_control_window(ApplicationData* data);
static void draw_representations_window(ApplicationData* data);
static void draw_property_window(ApplicationData* data);
static void draw_timeline_window(ApplicationData* data);
static void draw_distribution_window(ApplicationData* data);
static void draw_ramachandran_window(ApplicationData* data);
static void draw_atom_info_window(const MoleculeStructure& mol, int atom_idx, int x, int y);
static void draw_async_info(ApplicationData* data);
// static void draw_selection_window(ApplicationData* data);

static void init_framebuffer(MainFramebuffer* fbo, int width, int height);
static void destroy_framebuffer(MainFramebuffer* fbo);

static void init_molecule_buffers(ApplicationData* data);
static void free_molecule_buffers(ApplicationData* data);

static void copy_molecule_data_to_buffers(ApplicationData* data);

static void init_molecule_data(ApplicationData* data);
static void init_trajectory_data(ApplicationData* data);

static void load_molecule_data(ApplicationData* data, CString file);
static void free_molecule_data(ApplicationData* data);

static void load_workspace(ApplicationData* data, CString file);
static void save_workspace(ApplicationData* data, CString file);

static void create_screenshot(ApplicationData* data);

// Representations
static Representation* create_representation(ApplicationData* data, RepresentationType type = RepresentationType::Vdw, ColorMapping color_mapping = ColorMapping::Cpk, CString filter = "all");
static Representation* clone_representation(ApplicationData* data, const Representation& rep);
static void remove_representation(ApplicationData* data, int idx);
static void update_representation(ApplicationData* data, Representation* rep);
static void reset_representations(ApplicationData* data);
static void clear_representations(ApplicationData* data);

// Selections
static Selection* create_selection(ApplicationData* data, CString name, Array<const bool> atom_mask);
static Selection* clone_selection(ApplicationData* data, const Selection& sel);
static void remove_selection(ApplicationData* data, int idx);

static void reset_selections(ApplicationData* data);
static void clear_selections(ApplicationData* data);

static void create_volume(ApplicationData* data);

// Async operations
static void load_trajectory_async(ApplicationData* data);
static void compute_backbone_angles_async(ApplicationData* data);

int main(int, char**) {
    ApplicationData data;

    // Init logging
    logging::initialize();
    logging::register_backend([](CString str, logging::Severity, void*) {
        print_string(str);
        printf("\n");
    });
    logging::register_backend(
        [](CString str, logging::Severity severity, void* usr_data) {
            const char* modifier = "";
            switch (severity) {
                case logging::Note:
                    modifier = "[note] ";
                    break;
                case logging::Warning:
                    modifier = "[warning] ";
                    break;
                case logging::Error:
                    modifier = "[error] ";
                    break;
                case logging::Fatal:
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
    ramachandran::initialize();
    LOG_NOTE("Initializing stats...");
    stats::initialize();
    LOG_NOTE("Initializing filter...");
    filter::initialize();
    LOG_NOTE("Initializing post processing...");
    postprocessing::initialize(data.fbo.width, data.fbo.height);
    LOG_NOTE("Initializing volume...");
    volume::initialize();

    // Setup IMGUI style
    {
        ImGui::StyleColorsClassic();
        // ImGui::StyleColorsLight();
        auto& style = ImGui::GetStyle();
        style.WindowRounding = 0.0f;
        style.Colors[ImGuiCol_TitleBgCollapsed] = ImVec4(0.40f, 0.40f, 0.80f, 0.30f);

        /*
        style.Colors[ImGuiCol_Text] = ImVec4(0.00f, 0.00f, 0.00f, 1.00f);
        style.Colors[ImGuiCol_TextDisabled] = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
        //style.Colors[ImGuiCol_TextHovered] = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
        //style.Colors[ImGuiCol_TextActive] = ImVec4(1.00f, 1.00f, 0.00f, 1.00f);
        style.Colors[ImGuiCol_WindowBg] = ImVec4(0.94f, 0.94f, 0.94f, 1.00f);
        style.Colors[ImGuiCol_ChildWindowBg] = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
        style.Colors[ImGuiCol_Border] = ImVec4(0.00f, 0.00f, 0.00f, 0.39f);
        style.Colors[ImGuiCol_BorderShadow] = ImVec4(1.00f, 1.00f, 1.00f, 0.10f);
        style.Colors[ImGuiCol_FrameBg] = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
        style.Colors[ImGuiCol_FrameBgHovered] = ImVec4(0.26f, 0.59f, 0.98f, 0.40f);
        style.Colors[ImGuiCol_FrameBgActive] = ImVec4(0.26f, 0.59f, 0.98f, 0.67f);
        style.Colors[ImGuiCol_TitleBg] = ImVec4(0.96f, 0.96f, 0.96f, 1.00f);
        style.Colors[ImGuiCol_TitleBgCollapsed] = ImVec4(1.00f, 1.00f, 1.00f, 0.51f);
        style.Colors[ImGuiCol_TitleBgActive] = ImVec4(0.82f, 0.82f, 0.82f, 1.00f);
        style.Colors[ImGuiCol_MenuBarBg] = ImVec4(0.86f, 0.86f, 0.86f, 1.00f);
        style.Colors[ImGuiCol_ScrollbarBg] = ImVec4(0.98f, 0.98f, 0.98f, 0.53f);
        style.Colors[ImGuiCol_ScrollbarGrab] = ImVec4(0.69f, 0.69f, 0.69f, 0.80f);
        style.Colors[ImGuiCol_ScrollbarGrabHovered] = ImVec4(0.49f, 0.49f, 0.49f, 0.80f);
        style.Colors[ImGuiCol_ScrollbarGrabActive] = ImVec4(0.49f, 0.49f, 0.49f, 1.00f);
        //style.Colors[ImGuiCol_ComboBg] = ImVec4(0.86f, 0.86f, 0.86f, 0.99f);
        style.Colors[ImGuiCol_CheckMark] = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
        style.Colors[ImGuiCol_SliderGrab] = ImVec4(0.26f, 0.59f, 0.98f, 0.78f);
        style.Colors[ImGuiCol_SliderGrabActive] = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
        style.Colors[ImGuiCol_Button] = ImVec4(0.26f, 0.59f, 0.98f, 0.40f);
        style.Colors[ImGuiCol_ButtonHovered] = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
        style.Colors[ImGuiCol_ButtonActive] = ImVec4(0.06f, 0.53f, 0.98f, 1.00f);
        style.Colors[ImGuiCol_Header] = ImVec4(0.26f, 0.59f, 0.98f, 0.31f);
        style.Colors[ImGuiCol_HeaderHovered] = ImVec4(0.26f, 0.59f, 0.98f, 0.80f);
        style.Colors[ImGuiCol_HeaderActive] = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
        style.Colors[ImGuiCol_Column] = ImVec4(0.39f, 0.39f, 0.39f, 1.00f);
        style.Colors[ImGuiCol_ColumnHovered] = ImVec4(0.26f, 0.59f, 0.98f, 0.78f);
        style.Colors[ImGuiCol_ColumnActive] = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
        style.Colors[ImGuiCol_ResizeGrip] = ImVec4(1.00f, 1.00f, 1.00f, 0.00f);
        style.Colors[ImGuiCol_ResizeGripHovered] = ImVec4(0.26f, 0.59f, 0.98f, 0.67f);
        style.Colors[ImGuiCol_ResizeGripActive] = ImVec4(0.26f, 0.59f, 0.98f, 0.95f);
        //style.Colors[ImGuiCol_CloseButton] = ImVec4(0.59f, 0.59f, 0.59f, 0.50f);
        //style.Colors[ImGuiCol_CloseButtonHovered] = ImVec4(0.98f, 0.39f, 0.36f, 1.00f);
        //style.Colors[ImGuiCol_CloseButtonActive] = ImVec4(0.98f, 0.39f, 0.36f, 1.00f);
        style.Colors[ImGuiCol_PlotLines] = ImVec4(0.39f, 0.39f, 0.39f, 1.00f);
        style.Colors[ImGuiCol_PlotLinesHovered] = ImVec4(1.00f, 0.43f, 0.35f, 1.00f);
        style.Colors[ImGuiCol_PlotHistogram] = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
        style.Colors[ImGuiCol_PlotHistogramHovered] = ImVec4(1.00f, 0.60f, 0.00f, 1.00f);
        style.Colors[ImGuiCol_TextSelectedBg] = ImVec4(0.26f, 0.59f, 0.98f, 0.35f);
        //style.Colors[ImGuiCol_TooltipBg] = ImVec4(1.00f, 1.00f, 1.00f, 0.94f);
        style.Colors[ImGuiCol_ModalWindowDarkening] = ImVec4(0.20f, 0.20f, 0.20f, 0.35f);
        */
    }

    // const vec4 CLEAR_COLOR = vec4(0, 0, 0, 0);
    const vec4 CLEAR_INDEX = vec4(1, 1, 1, 1);

    vec2 halton_23[16];
    math::generate_halton_sequence(halton_23, ARRAY_SIZE(halton_23), 2, 3);

#ifdef VIAMD_RELEASE
    allocate_and_parse_pdb_from_string(&data.mol_data.dynamic, CAFFINE_PDB);
    init_molecule_data(&data);
#else
    load_molecule_data(&data, VIAMD_DATA_DIR "/1af6.pdb");
#endif
    reset_view(&data, true);
    create_representation(&data, RepresentationType::Vdw, ColorMapping::ResId);
    create_volume(&data);

    // Main loop
    while (!data.ctx.window.should_close) {
        platform::Coordinate previous_mouse_coord = data.ctx.input.mouse.win_coord;
        platform::update(&data.ctx);

        // Try to fix false move on touch
        if (data.ctx.input.mouse.hit[0]) {
            previous_mouse_coord = data.ctx.input.mouse.win_coord;
        }

        // #input
        if (data.ctx.input.key.hit[CONSOLE_KEY]) {
            if (data.console.Visible()) {
                data.console.Hide();
            } else if (!ImGui::GetIO().WantTextInput) {
                data.console.Show();
            }
        }

        if (!ImGui::GetIO().WantTextInput) {
            if (data.ctx.input.key.hit[Key::KEY_F5]) {
                draw::initialize();
                postprocessing::initialize(data.fbo.width, data.fbo.height);
            }

            if (data.ctx.input.key.hit[PLAY_PAUSE_KEY]) {
                const int32 num_frames = data.mol_data.dynamic.trajectory ? data.mol_data.dynamic.trajectory.num_frames : 0;
                const float64 max_time = (float64)math::max(0, num_frames - 1);
                if (!data.is_playing && data.time == max_time) {
                    data.time = 0;
                }
                data.is_playing = !data.is_playing;
            }
        }

        data.selection.selecting = false;

        if (!ImGui::GetIO().WantCaptureMouse) {
            // #selection
            bool selecting = handle_selection(&data);
            data.selection.selecting |= selecting;

            // #camera-control
            if (!selecting) {
                // CAMERA CONTROLS
                data.view.trackball_state.input.rotate_button = data.ctx.input.mouse.down[0];
                data.view.trackball_state.input.pan_button = data.ctx.input.mouse.down[1];
                data.view.trackball_state.input.dolly_button = data.ctx.input.mouse.down[2];
                data.view.trackball_state.input.mouse_coord_prev = {previous_mouse_coord.x, previous_mouse_coord.y};
                data.view.trackball_state.input.mouse_coord_curr = {data.ctx.input.mouse.win_coord.x, data.ctx.input.mouse.win_coord.y};
                data.view.trackball_state.input.screen_size = vec2(data.ctx.window.width, data.ctx.window.height);
                data.view.trackball_state.input.dolly_delta = data.ctx.input.mouse.scroll_delta;

                {
                    vec3 pos = data.view.animation.target_position;
                    quat ori = data.view.camera.orientation;
                    if (camera_controller_trackball(&pos, &ori, &data.view.trackball_state, TrackballFlags_RotateReturnsTrue | TrackballFlags_PanReturnsTrue)) {
                        data.view.camera.position = pos;
                        data.view.camera.orientation = ori;
                    }
                    data.view.animation.target_position = pos;
                }

                if (ImGui::GetIO().MouseDoubleClicked[0]) {
                    if (data.picking.depth < 1.f) {
                        const vec3 forward = data.view.camera.orientation * vec3(0, 0, 1);
                        const float32 dist = data.view.trackball_state.distance;
                        const vec3 camera_target_pos = data.picking.world_coord + forward * dist;

                        data.view.animation.target_position = camera_target_pos;
                    }
                }
            }
        }
        // #animate-camera
        {
            const float32 dt = math::min(data.ctx.timing.delta_s, 0.033f);
            const float32 speed = 10.0f;

            const vec3 vel = (data.view.animation.target_position - data.view.camera.position) * speed;
            data.view.camera.position += vel * dt;
#if 0
            ImGui::Begin("Camera Debug Info");
            ImGui::Text("lin vel [%.2f %.2f %.2f]", vel.x, vel.y, vel.z);
            ImGui::Text("lin cur [%.2f %.2f %.2f]", data.camera.camera.position.x, data.camera.camera.position.y, data.camera.camera.position.z);
            ImGui::Text("lin tar [%.2f %.2f %.2f]", data.camera.animation.target_position.x, data.camera.animation.target_position.y,
                        data.camera.animation.target_position.z);
            ImGui::Text("ang cur [%.2f %.2f %.2f %.2f]", data.camera.camera.orientation.x, data.camera.camera.orientation.y,
                        data.camera.camera.orientation.z, data.camera.camera.orientation.w);
            ImGui::Text("ang tar [%.2f %.2f %.2f %.2f]", data.camera.animation.target_orientation.x, data.camera.animation.target_orientation.y,
                        data.camera.animation.target_orientation.z, data.camera.animation.target_orientation.w);
            ImGui::End();
#endif
        }

        // This needs to happen first (in imgui events) to enable docking of imgui windows
        // ImGui::CreateDockspace();

        {
            if (data.density_volume.enabled) {
                stats::async_update(
                    data.mol_data.dynamic, {(int32)data.time_filter.range.beg, (int32)data.time_filter.range.end},
                    [](void* usr_data) {
                        ApplicationData* data = (ApplicationData*)usr_data;
                        data->density_volume.volume_data_mutex.lock();

                        stats::compute_density_volume(&data->density_volume.volume, data->density_volume.world_to_texture_matrix, data->mol_data.dynamic.trajectory,
                                                      {(int32)data->time_filter.range.beg, (int32)data->time_filter.range.end});

                        data->density_volume.volume_data_mutex.unlock();
                        data->density_volume.texture.dirty = true;
                    },
                    &data);
            } else {
                stats::async_update(data.mol_data.dynamic, {(int32)data.time_filter.range.beg, (int32)data.time_filter.range.end});
            }
        }

        // If gpu representation of volume is not up to date, upload data
        if (data.density_volume.texture.dirty) {
            if (data.density_volume.volume_data_mutex.try_lock()) {
                if (data.density_volume.texture.dim != data.density_volume.volume.dim) {
                    data.density_volume.texture.dim = data.density_volume.volume.dim;
                    volume::create_volume_texture(&data.density_volume.texture.id, data.density_volume.texture.dim);
                }

                volume::set_volume_texture_data(data.density_volume.texture.id, data.density_volume.texture.dim, data.density_volume.volume.voxel_data.ptr);
                data.density_volume.volume_data_mutex.unlock();
                data.density_volume.texture.max_value = data.density_volume.volume.voxel_range.y;
                data.density_volume.texture.dirty = false;
            }
        }

        bool time_changed = false;
        bool frame_changed = false;

        if (data.is_playing) {
            const int32 num_frames = data.mol_data.dynamic.trajectory ? data.mol_data.dynamic.trajectory.num_frames : 0;
            const float64 max_time = (float64)math::max(0, num_frames - 1);
            data.time += data.ctx.timing.delta_s * data.frames_per_second;
            data.time = math::clamp(data.time, 0.0, max_time);
            if (data.time >= max_time) {
                data.is_playing = false;
                data.time = max_time;
            }
        }

        {
            static float64 prev_time = data.time;
            static int32 prev_frame = (int32)data.time;

            if (data.time != prev_time) {
                time_changed = true;
            }
            prev_time = data.time;

            int32 frame = (int32)data.time;
            if (frame != prev_frame) {
                frame_changed = true;
            }
            prev_frame = frame;
        }

        if (data.time_filter.dynamic_window) {
            const float32 max_frame = data.mol_data.dynamic.trajectory ? (float32)data.mol_data.dynamic.trajectory.num_frames : 1.f;
            data.time_filter.range.min = math::clamp((int32)data.time - data.time_filter.window_extent / 2.f, 0.f, max_frame);
            data.time_filter.range.max = math::clamp((int32)data.time + data.time_filter.window_extent / 2.f, 0.f, max_frame);
        }

        if (frame_changed) {
            if (data.mol_data.dynamic.trajectory) {
                if (data.time_filter.dynamic_window) {
                    stats::set_all_property_flags(false, true);
                }
            }
        }

        if (time_changed) {
            auto& mol = data.mol_data.dynamic.molecule;
            auto& traj = data.mol_data.dynamic.trajectory;

            data.hydrogen_bonds.dirty = true;
            data.gpu_buffers.dirty.backbone = true;

            PUSH_CPU_SECTION("Interpolate Position")
            if (traj) {
                auto pos = get_positions(mol);
                vec3* old_pos = (vec3*)TMP_MALLOC(pos.size_in_bytes());
                defer { TMP_FREE(old_pos); };
                memcpy(old_pos, pos.data(), pos.size_in_bytes());

                const int current_frame = math::clamp((int)data.time, 0, math::max(0, data.mol_data.dynamic.trajectory.num_frames - 1));
                const vec3 box_ext = get_trajectory_frame(data.mol_data.dynamic.trajectory, current_frame).box * vec3(1.0f);

                interpolate_atomic_positions(pos, traj, data.time, data.interpolation);
                // compute_atomic_velocities(data.mol_data.atom_velocity, pos, Array<const vec3>(old_pos, pos.size()), box_ext);
#if 0
                if (data.interpolation != PlaybackInterpolationMode::Nearest) {
                    const auto& box = get_trajectory_frame(data.mol_data.dynamic.trajectory, (int)data.time).box;
                    apply_pbc_residues(get_positions(data.mol_data.dynamic.molecule), data.mol_data.dynamic.molecule.residues, box);
                    // apply_pbc_chains(get_positions(data.mol_data.dynamic.molecule), data.mol_data.dynamic.molecule.chains, data.mol_data.dynamic.molecule.residues, box);
                }
#endif

                data.gpu_buffers.dirty.position = true;
                data.gpu_buffers.dirty.velocity = true;
            }
            POP_CPU_SECTION()

            PUSH_CPU_SECTION("Compute backbone angles")
            zero_array(mol.backbone.angles);
            compute_backbone_angles(mol.backbone.angles, get_positions(mol), mol.backbone.segments, mol.backbone.sequences);
            POP_CPU_SECTION()

            PUSH_CPU_SECTION("Update dynamic representations")
            for (auto& rep : data.representations.buffer) {
                if (rep.color_mapping == ColorMapping::SecondaryStructure) {
                    update_representation(&data, &rep);
                }
            }
            POP_CPU_SECTION()

#if 0
			PUSH_CPU_SECTION("Spatial Hash") {
				if (data.spatial_hash.dirty_flag) {
					data.spatial_hash.dirty_flag = false;
					spatialhash::compute_frame(&data.spatial_hash.frame, get_positions(data.mol_data.dynamic.molecule), data.spatial_hash.cell_ext);
				}
			}
			POP_CPU_SECTION()
#endif
        } else {  // Not time_changed
            // @FIXME: DO THIS ONLY ONCE
            // Clear velocity buffer (once)
            // zero_array(data.mol_data.atom_velocity);
            // copy_molecule_data_to_buffers(&data);
        }

        PUSH_CPU_SECTION("Hydrogen bonds")
        if (data.hydrogen_bonds.enabled && data.hydrogen_bonds.dirty) {
            data.hydrogen_bonds.bonds.clear();
            hydrogen_bond::compute_bonds(&data.hydrogen_bonds.bonds, data.mol_data.dynamic.molecule.hydrogen_bond.donors, data.mol_data.dynamic.molecule.hydrogen_bond.acceptors,
                                         get_positions(data.mol_data.dynamic.molecule), data.hydrogen_bonds.distance_cutoff, data.hydrogen_bonds.angle_cutoff * math::DEG_TO_RAD);
            data.hydrogen_bonds.dirty = false;
        }
        POP_CPU_SECTION()

        if (data.async.trajectory.sync.running) {
            constexpr float32 TICK_INTERVAL_SEC = 3.f;
            static float32 time = 0.f;
            time += data.ctx.timing.delta_s;
            if (time > TICK_INTERVAL_SEC) {
                time = 0.f;
                stats::set_all_property_flags(true, true);
                compute_backbone_angles_async(&data);
            }
        }

        bool visuals_changed = false;
        {
            static auto old_hash = hash::crc64(&data.visuals, sizeof(data.visuals));
            const auto new_hash = hash::crc64(&data.visuals, sizeof(data.visuals));
            visuals_changed = (new_hash != old_hash);
            old_hash = new_hash;
        }

        // Resize Framebuffer
        if ((data.fbo.width != data.ctx.framebuffer.width || data.fbo.height != data.ctx.framebuffer.height) && (data.ctx.framebuffer.width != 0 && data.ctx.framebuffer.height != 0)) {
            init_framebuffer(&data.fbo, data.ctx.framebuffer.width, data.ctx.framebuffer.height);
            postprocessing::initialize(data.fbo.width, data.fbo.height);
        }

        PUSH_GPU_SECTION("Compute Backbone Spline") {
            bool has_spline_rep = false;
            for (const auto& rep : data.representations.buffer) {
                if (rep.type == RepresentationType::Ribbons || rep.type == RepresentationType::Cartoon) {
                    has_spline_rep = true;
                    break;
                }
            }
            has_spline_rep |= data.visuals.spline.draw_control_points || data.visuals.spline.draw_spline;

            data.gpu_buffers.dirty.backbone = true;
            if (has_spline_rep && data.gpu_buffers.dirty.backbone) {
                data.gpu_buffers.dirty.backbone = false;
                draw::compute_backbone_control_points(data.gpu_buffers.backbone.control_point, data.gpu_buffers.position, data.gpu_buffers.backbone.backbone_segment_index,
                                                      data.gpu_buffers.backbone.num_backbone_segment_indices, ramachandran::get_segmentation_texture());
                draw::compute_backbone_spline(data.gpu_buffers.backbone.spline, data.gpu_buffers.backbone.control_point, data.gpu_buffers.backbone.control_point_index,
                                              data.gpu_buffers.backbone.num_control_point_indices);
            }
        }
        POP_GPU_SECTION()

        {
            mat4 view_mat = compute_world_to_view_matrix(data.view.camera);
            mat4 proj_mat = compute_perspective_projection_matrix(data.view.camera, data.fbo.width, data.fbo.height);

            const vec2 res = vec2(data.fbo.width, data.fbo.height);
            vec2 jitter = vec2(0, 0);
            if (data.visuals.temporal_reprojection.enabled && data.visuals.temporal_reprojection.jitter) {
                static uint32 i = 0;
                i = (++i) % ARRAY_SIZE(halton_23);
                jitter = halton_23[i] - 0.5f;
                proj_mat = compute_perspective_projection_matrix(data.view.camera, data.fbo.width, data.fbo.height, jitter.x, jitter.y);
            }

            auto& param = data.view.param;
            param.previous.matrix.proj = param.matrix.proj;
            param.previous.matrix.view = param.matrix.view;
            param.previous.matrix.view_proj = param.matrix.view_proj;
            param.previous.jitter = param.jitter;

            param.matrix.view = view_mat;
            param.matrix.proj = proj_mat;
            param.matrix.view_proj = proj_mat * view_mat;

            param.matrix.inverse.view = math::inverse(view_mat);
            param.matrix.inverse.proj = math::inverse(proj_mat);
            param.matrix.inverse.view_proj = math::inverse(param.matrix.view_proj);
            param.matrix.norm = math::transpose(param.matrix.inverse.view);

            param.jitter = jitter;
            param.resolution = res;
        }

        const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3, GL_COLOR_ATTACHMENT4};

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, data.fbo.deferred.fbo);
        glViewport(0, 0, data.fbo.width, data.fbo.height);

        glDepthMask(GL_TRUE);
        glColorMask(1, 1, 1, 1);

        // Setup fbo and clear textures
        PUSH_GPU_SECTION("Clear G-buffer") {
            // Clear color+alpha, normal, velocity, emissive and depth
            glDrawBuffers(4, draw_buffers);
            glClearColor(0, 0, 0, 0);
            glClearDepthf(1.f);
            glStencilMask(0xFF);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

            // Clear picking buffer
            glDrawBuffer(GL_COLOR_ATTACHMENT4);
            glClearColor(CLEAR_INDEX.x, CLEAR_INDEX.y, CLEAR_INDEX.z, CLEAR_INDEX.w);
            glClear(GL_COLOR_BUFFER_BIT);
        }
        POP_GPU_SECTION()

        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);

        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        glDepthMask(GL_TRUE);

        // Enable all draw buffers
        glDrawBuffers(ARRAY_SIZE(draw_buffers), draw_buffers);

        PUSH_GPU_SECTION("G-Buffer fill") {
            for (const auto& rep : data.representations.buffer) {
                if (!rep.enabled) continue;
                switch (rep.type) {
                    case RepresentationType::Vdw:
                        PUSH_GPU_SECTION("Vdw")
                        draw::draw_vdw(data.gpu_buffers.position, data.gpu_buffers.radius, rep.color_buffer, data.gpu_buffers.velocity, (int)data.mol_data.dynamic.molecule.atom.count, data.view.param,
                                       rep.radius);
                        POP_GPU_SECTION()
                        break;
                    case RepresentationType::Licorice:
                        PUSH_GPU_SECTION("Licorice")
                        draw::draw_licorice(data.gpu_buffers.position, rep.color_buffer, data.gpu_buffers.velocity, data.gpu_buffers.bond, (int)data.mol_data.dynamic.molecule.covalent_bonds.size(),
                                            data.view.param, rep.radius);
                        POP_GPU_SECTION()
                        break;
                    case RepresentationType::BallAndStick:
                        PUSH_GPU_SECTION("Vdw")
                        draw::draw_vdw(data.gpu_buffers.position, data.gpu_buffers.radius, rep.color_buffer, data.gpu_buffers.velocity, (int)data.mol_data.dynamic.molecule.atom.count, data.view.param,
                                       rep.radius * BALL_AND_STICK_VDW_SCALE);
                        POP_GPU_SECTION()
                        PUSH_GPU_SECTION("Licorice")
                        draw::draw_licorice(data.gpu_buffers.position, rep.color_buffer, data.gpu_buffers.velocity, data.gpu_buffers.bond, (int)data.mol_data.dynamic.molecule.covalent_bonds.size(),
                                            data.view.param, rep.radius * BALL_AND_STICK_LICORICE_SCALE);
                        POP_GPU_SECTION()
                        break;
                    case RepresentationType::Ribbons:
                        PUSH_GPU_SECTION("Ribbons")
                        draw::draw_ribbons(data.gpu_buffers.backbone.spline, data.gpu_buffers.backbone.spline_index, rep.color_buffer, data.gpu_buffers.velocity,
                                           data.gpu_buffers.backbone.num_spline_indices, data.view.param);
                        POP_GPU_SECTION()
                        break;
                    case RepresentationType::Cartoon:
                        PUSH_GPU_SECTION("Cartoon")
                        draw::draw_cartoon(data.gpu_buffers.backbone.spline, data.gpu_buffers.backbone.spline_index, rep.color_buffer, data.gpu_buffers.backbone.num_spline_indices, data.view.param);
                        POP_GPU_SECTION()
                        break;
                }
            }

            // RENDER DEBUG INFORMATION (WITH DEPTH)
            PUSH_GPU_SECTION("Debug Draw") {
                immediate::set_view_matrix(data.view.param.matrix.view);
                immediate::set_proj_matrix(data.view.param.matrix.proj);
                // immediate::Material plane_mat = immediate::MATERIAL_GLOSSY_WHITE;
                // plane_mat.f0 = {data.immediate_gfx.material.f0, data.immediate_gfx.material.f0, data.immediate_gfx.material.f0};
                // plane_mat.smoothness = data.immediate_gfx.material.smoothness;
                // immediate::set_material(plane_mat);
                // immediate::draw_plane({-30, -30, -50}, {100, 0, 0}, {0, 0, 100});

                // HYDROGEN BONDS
                if (data.hydrogen_bonds.enabled && !data.hydrogen_bonds.overlay) {
                    for (const auto& bond : data.hydrogen_bonds.bonds) {
                        immediate::draw_line(data.mol_data.dynamic.molecule.atom.positions[bond.acc_idx], data.mol_data.dynamic.molecule.atom.positions[bond.hyd_idx],
                                             math::convert_color(data.hydrogen_bonds.color));
                    }
                }

                // SIMULATION BOX
                if (data.simulation_box.enabled && data.mol_data.dynamic.trajectory.num_frames > 0) {
                    auto frame_idx = math::clamp((int)data.time, 0, data.mol_data.dynamic.trajectory.num_frames - 1);
                    auto frame = get_trajectory_frame(data.mol_data.dynamic.trajectory, frame_idx);
                    immediate::draw_aabb_lines(vec3(0), frame.box * vec3(1), math::convert_color(data.simulation_box.color));
                }

                immediate::flush();
            }
            POP_GPU_SECTION()
        }
        POP_GPU_SECTION()  // G-buffer

		const bool atom_selection_empty = is_array_zero(data.selection.current_selection_mask);
		const bool atom_highlight_empty = is_array_zero(data.selection.current_highlight_mask);

        PUSH_GPU_SECTION("Selection") {
            const uint32 atom_count = (uint32)data.mol_data.dynamic.molecule.atom.count;
            glDepthMask(0);

            glDrawBuffer(GL_COLOR_ATTACHMENT3);  // Emission buffer
            glEnable(GL_DEPTH_TEST);
            glEnable(GL_STENCIL_TEST);

            {
                glStencilOp(GL_KEEP, GL_REPLACE, GL_REPLACE);

                glColorMask(0, 0, 0, 0);
                glStencilMask(0xFF);
                glStencilFunc(GL_ALWAYS, 2, 0xFF);

                // const vec4 fill_color = vec4(0, 0, 10, 0);
				draw_representations_lean_and_mean(&data);
            }

            {
                glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
                glDepthFunc(GL_LEQUAL);

                glColorMask(0, 0, 0, 0);
                glStencilMask(0xFF);
                glStencilFunc(GL_ALWAYS, 4, 0xFF);

                //const vec4 visible_color = vec4(1, 1, 1, 0);
				draw_representations_lean_and_mean(&data);

            }

            glDisable(GL_DEPTH_TEST);

            {
				// Highlight
                glStencilFunc(GL_GREATER, 1, 0xFF);
                glStencilOp(GL_KEEP, GL_REPLACE, GL_REPLACE);
                glStencilMask(0x00);
                glColorMask(1, 1, 1, 0);

                const vec4 color = vec4(0, 0.5, 1.0, 0) * 5.0f;
				draw_representations_lean_and_mean(&data, color, 1.2f, 0x02);
            }

			{
				glStencilFunc(GL_GREATER, 1, 0xFF);
				glStencilOp(GL_KEEP, GL_REPLACE, GL_REPLACE);
				glStencilMask(0x00);
				glColorMask(1, 1, 1, 0);

				const vec4 color = vec4(1.0, 0.5, 0, 0) * 5.0f;
				draw_representations_lean_and_mean(&data, color, 1.2f, 0x01);
			}


            if (!atom_selection_empty) {
                glStencilFunc(GL_NOTEQUAL, 4, 0xFF);
                glStencilMask(0x00);
                glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
                glDrawBuffer(GL_COLOR_ATTACHMENT0);
                postprocessing::scale_hsv(data.fbo.deferred.color, vec3(1, 0.15, 1));
            }

            glDisable(GL_STENCIL_TEST);
            glEnable(GL_DEPTH_TEST);
			glDepthFunc(GL_LESS);
            glDepthMask(1);
        }
        POP_GPU_SECTION()

        // PICKING
        PUSH_GPU_SECTION("Picking") {
            vec2 coord = {data.ctx.input.mouse.win_coord.x, data.fbo.height - data.ctx.input.mouse.win_coord.y};
            if (coord.x < 0.f || coord.x >= (float)data.fbo.width || coord.y < 0.f || coord.y >= (float)data.fbo.height) {
                data.picking.idx = NO_PICKING_IDX;
                data.picking.depth = 1.f;
            } else {
                static uint32 frame_idx = 0;
                static uint32 ref_frame = 0;
                frame_idx = (frame_idx + 1) % 16;
                // @NOTE: If we have jittering applied, we cannot? retreive the original pixel value (without the jitter)
                // Solution, pick one reference frame out of the jittering sequence and use that one...
                // Ugly hack but works...

                if (data.ctx.input.mouse.moving) {
                    ref_frame = frame_idx;
                }

                if (ref_frame == frame_idx || data.view.param.jitter == vec2(0, 0)) {
                    data.picking = read_picking_data(data.fbo, (int32)math::round(coord.x), (int32)math::round(coord.y));
                    if (data.picking.idx != NO_PICKING_IDX) data.picking.idx = math::clamp(data.picking.idx, 0U, (uint32)data.mol_data.dynamic.molecule.atom.count - 1U);
                    const vec4 viewport(0, 0, data.fbo.width, data.fbo.height);
                    data.picking.world_coord = math::unproject(vec3(coord.x, coord.y, data.picking.depth), data.view.param.matrix.inverse.view_proj, viewport);
                }
            }

            data.selection.hovered = -1;
            if (data.picking.idx != NO_PICKING_IDX) {
                data.selection.hovered = data.picking.idx;
            }
            if (data.ctx.input.mouse.clicked[1]) {
                data.selection.right_clicked = data.selection.hovered;
            }
        }
        POP_GPU_SECTION()

        glDisable(GL_DEPTH_TEST);
        glDepthMask(GL_FALSE);

#if 0
        PUSH_GPU_SECTION("Highlight Selection") {

            const bool selection_empty = is_array_zero(data.selection.current_selection_mask);

            // glDrawBuffer(GL_COLOR_ATTACHMENT3);  // emission as intermediate target
            // postprocessing::desaturate_selection(data.fbo.deferred.color, data.fbo.deferred.picking, data.gpu_buffers.selection, data.selection.selecting || !selection_empty);

            // glDrawBuffer(GL_COLOR_ATTACHMENT0);
            // postprocessing::blit_texture(data.fbo.deferred.emissive);

            /*
glDrawBuffer(GL_COLOR_ATTACHMENT3);  // Emission buffer
glEnable(GL_BLEND);
glBlendFunc(GL_ONE, GL_ONE);
const vec3 highlight = data.selection.highlight_color * data.selection.highlight_scale;
const vec3 selection = data.selection.selection_color * data.selection.selection_scale;
const vec3 outline = data.selection.outline_color * data.selection.outline_scale;
postprocessing::highlight_selection(data.fbo.deferred.picking, data.gpu_buffers.selection, highlight, selection, outline);
glDisable(GL_BLEND);
            */
        }
        POP_GPU_SECTION();
#endif

        // Activate backbuffer
        glViewport(0, 0, data.ctx.framebuffer.width, data.ctx.framebuffer.height);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glDrawBuffer(GL_BACK);

        PUSH_GPU_SECTION("Postprocessing") {
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

            data.visuals.dof.focus_depth = data.view.trackball_state.distance;

            desc.depth_of_field.enabled = data.visuals.dof.enabled;
            desc.depth_of_field.focus_depth = data.visuals.dof.focus_depth;
            desc.depth_of_field.focus_scale = data.visuals.dof.focus_scale;

            desc.temporal_reprojection.enabled = data.visuals.temporal_reprojection.enabled;
            desc.temporal_reprojection.feedback_min = data.visuals.temporal_reprojection.feedback_min;
            desc.temporal_reprojection.feedback_max = data.visuals.temporal_reprojection.feedback_max;
            desc.temporal_reprojection.motion_blur.enabled = data.visuals.temporal_reprojection.motion_blur.enabled;
            desc.temporal_reprojection.motion_blur.motion_scale = data.visuals.temporal_reprojection.motion_blur.motion_scale;

            desc.input_textures.depth = data.fbo.deferred.depth;
            desc.input_textures.color = data.fbo.deferred.color;
            desc.input_textures.normal = data.fbo.deferred.normal;
            desc.input_textures.velocity = data.fbo.deferred.velocity;
            desc.input_textures.emissive = data.fbo.deferred.emissive;

            postprocessing::shade_and_postprocess(desc, data.view.param);
        }
        POP_GPU_SECTION()

        // DRAW DEBUG GRAPHICS W/O DEPTH
        PUSH_GPU_SECTION("Debug Draw Overlay") {
            immediate::set_view_matrix(data.view.param.matrix.view);
            immediate::set_proj_matrix(data.view.param.matrix.proj);
            stats::visualize(data.mol_data.dynamic);

            // HYDROGEN BONDS
            if (data.hydrogen_bonds.enabled && data.hydrogen_bonds.overlay) {
                for (const auto& bond : data.hydrogen_bonds.bonds) {
                    immediate::draw_line(data.mol_data.dynamic.molecule.atom.positions[bond.acc_idx], data.mol_data.dynamic.molecule.atom.positions[bond.hyd_idx],
                                         math::convert_color(data.hydrogen_bonds.color));
                }
            }
            immediate::flush();

            PUSH_GPU_SECTION("Draw Control Points")
            if (data.visuals.spline.draw_control_points) {
                draw::draw_spline(data.gpu_buffers.backbone.control_point, data.gpu_buffers.backbone.control_point_index, data.gpu_buffers.backbone.num_control_point_indices, data.view.param);
            }
            if (data.visuals.spline.draw_spline) {
                draw::draw_spline(data.gpu_buffers.backbone.spline, data.gpu_buffers.backbone.spline_index, data.gpu_buffers.backbone.num_spline_indices, data.view.param);
            }
            POP_GPU_SECTION()
        }
        POP_GPU_SECTION()

        // Activate backbuffer
        glViewport(0, 0, data.ctx.framebuffer.width, data.ctx.framebuffer.height);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glDrawBuffer(GL_BACK);

        if (data.density_volume.enabled) {
            PUSH_GPU_SECTION("Volume Rendering")
            const float32 scl = 1.f * data.density_volume.density_scale / data.density_volume.texture.max_value;
            volume::render_volume_texture(data.density_volume.texture.id, data.fbo.deferred.depth, data.density_volume.texture_to_model_matrix, data.density_volume.model_to_world_matrix,
                                          data.view.param.matrix.view, data.view.param.matrix.proj, data.density_volume.color, scl);
            POP_GPU_SECTION()
        }

        // GUI ELEMENTS
        data.console.Draw("VIAMD", data.ctx.window.width, data.ctx.window.height, data.ctx.timing.delta_s);

        draw_main_menu(&data);
        draw_context_popup(&data);

        if (data.representations.show_window) draw_representations_window(&data);
        if (data.statistics.show_property_window) draw_property_window(&data);
        if (data.statistics.show_timeline_window) draw_timeline_window(&data);
        if (data.statistics.show_distribution_window) draw_distribution_window(&data);
        // if (data.selection.show_window) draw_selection_window(&data);

        if (data.ramachandran.show_window) {
            draw_ramachandran_window(&data);
        }

        // ImGui::GetIO().WantCaptureMouse does not work with Menu
        if (!ImGui::IsMouseHoveringAnyWindow()) {
            if (data.picking.idx != NO_PICKING_IDX) {
                draw_atom_info_window(data.mol_data.dynamic.molecule, data.picking.idx, (int)data.ctx.input.mouse.win_coord.x, (int)data.ctx.input.mouse.win_coord.y);
            }
        }

        draw_async_info(&data);
        draw_animation_control_window(&data);

        PUSH_GPU_SECTION("Imgui render")
        platform::render_imgui(&data.ctx);
        POP_GPU_SECTION()

        // Swap buffers
        platform::swap_buffers(&data.ctx);
        data.frame++;

        PUSH_GPU_SECTION("Update Buffers")
        copy_molecule_data_to_buffers(&data);
        POP_GPU_SECTION()
    }

    data.async.trajectory.sync.signal_stop_and_wait();
    stats::signal_stop_and_wait();
    data.async.backbone_angles.sync.signal_stop_and_wait();

    destroy_framebuffer(&data.fbo);

    platform::shutdown(&data.ctx);

    return 0;
}

static void interpolate_atomic_positions(Array<vec3> dst_pos, const MoleculeTrajectory& traj, float64 time, PlaybackInterpolationMode interpolation_mode) {
    const int last_frame = traj.num_frames - 1;
    time = math::clamp(time, 0.0, float64(last_frame));

    const int frame = (int)time;
    const int prev_frame_2 = math::max(0, frame - 1);
    const int prev_frame_1 = math::max(0, frame);
    const int next_frame_1 = math::min(frame + 1, last_frame);
    const int next_frame_2 = math::min(frame + 2, last_frame);
    const mat3& box = get_trajectory_frame(traj, prev_frame_1).box;

    if (prev_frame_1 == next_frame_1) {
        copy_trajectory_positions(dst_pos, traj, prev_frame_1);
    } else {
        const float32 t = (float)math::fract(time);

        // INTERPOLATE
        switch (interpolation_mode) {
            case PlaybackInterpolationMode::Nearest: {
                const int nearest_frame = math::clamp((int)(time + 0.5), 0, last_frame);
                copy_trajectory_positions(dst_pos, traj, nearest_frame);
                break;
            }
            case PlaybackInterpolationMode::Linear: {
                const auto prev = get_trajectory_positions(traj, prev_frame_1);
                const auto next = get_trajectory_positions(traj, next_frame_1);
                linear_interpolation(dst_pos, prev, next, t);
                break;
            }
            case PlaybackInterpolationMode::LinearPbc: {
                const auto prev = get_trajectory_positions(traj, prev_frame_1);
                const auto next = get_trajectory_positions(traj, next_frame_1);
                linear_interpolation_periodic(dst_pos, prev, next, t, box);
                break;
            }
            case PlaybackInterpolationMode::Cubic: {
                const Array<const vec3> pos[4] = {get_trajectory_positions(traj, prev_frame_2), get_trajectory_positions(traj, prev_frame_1), get_trajectory_positions(traj, next_frame_1),
                                                  get_trajectory_positions(traj, next_frame_2)};
                cubic_interpolation(dst_pos, pos[0], pos[1], pos[2], pos[3], t);
                break;
            }
            case PlaybackInterpolationMode::CubicPbc: {
                const Array<const vec3> pos[4] = {get_trajectory_positions(traj, prev_frame_2), get_trajectory_positions(traj, prev_frame_1), get_trajectory_positions(traj, next_frame_1),
                                                  get_trajectory_positions(traj, next_frame_2)};
                cubic_interpolation_periodic(dst_pos, pos[0], pos[1], pos[2], pos[3], t, box);
                break;
            }

            default:
                break;
        }
    }
}

static void compute_atomic_velocities(Array<vec3> dst_vel, Array<const vec3> pos, Array<const vec3> old_pos, const vec3& box_ext) {
    ASSERT(dst_vel.size() == pos.size());
    ASSERT(dst_vel.size() == old_pos.size());

    for (int64 i = 0; i < dst_vel.size(); i++) {
        // De-periodize previous position
        const vec3 p1 = pos[i];
        const vec3 p0 = old_pos[i];

        const vec3 delta = p1 - p0;
        const vec3 signed_mask = sign(delta) * step(box_ext * 0.5f, abs(delta));
        const vec3 dp_p0 = p0 + box_ext * signed_mask;

        dst_vel[i] = p1 - dp_p0;
    }
}

// #misc
static float32 compute_avg_ms(float32 dt) {
    // @NOTE: Perhaps this can be done with a simple running mean?
    constexpr float32 interval = 0.5f;
    static float32 avg = 0.f;
    static int num_frames = 0;
    static float32 t = 0;
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
    if (!data->mol_data.dynamic.molecule) return;

    vec3 min_box, max_box;
    compute_bounding_box(&min_box, &max_box, get_positions(data->mol_data.dynamic.molecule));
    vec3 size = max_box - min_box;
    vec3 cent = (min_box + max_box) * 0.5f;
    vec3 pos = cent + size * 3.f;

    if (move_camera) {
        if (!smooth_transition) data->view.camera.position = pos;
        data->view.animation.target_position = pos;
        data->view.trackball_state.distance = math::length(pos - cent);
        look_at(&data->view.animation.target_position, &data->view.camera.orientation, cent, vec3(0, 1, 0));
    }

    data->view.camera.near_plane = 1.f;
    data->view.camera.far_plane = math::length(size) * 50.f;
}

// #picking
static PickingData read_picking_data(const MainFramebuffer& framebuffer, int32 x, int32 y) {
    static uint32 frame = 0;
    uint32 next = (frame + 1) % 2;

    PickingData data{};

    glBindFramebuffer(GL_READ_FRAMEBUFFER, framebuffer.deferred.fbo);
    glReadBuffer(GL_COLOR_ATTACHMENT4);

    // Queue async reads from current frame to pixel pack buffer
    glBindBuffer(GL_PIXEL_PACK_BUFFER, framebuffer.pbo_picking.color[frame]);
    glReadPixels(x, y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, framebuffer.pbo_picking.depth[frame]);
    glReadPixels(x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, 0);

    // Read values from previous frames pixel pack buffer
    glBindBuffer(GL_PIXEL_PACK_BUFFER, framebuffer.pbo_picking.color[next]);
    GLubyte* color = (GLubyte*)glMapBuffer(GL_PIXEL_PACK_BUFFER, GL_READ_ONLY);
    if (color) {
        data.idx = color[0] + (color[1] << 8) + (color[2] << 16) + (color[3] << 24);
        glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
    }

    glBindBuffer(GL_PIXEL_PACK_BUFFER, framebuffer.pbo_picking.depth[next]);
    GLfloat* depth = (GLfloat*)glMapBuffer(GL_PIXEL_PACK_BUFFER, GL_READ_ONLY);
    if (depth) {
        data.depth = depth[0];
        glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
    }

    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);

    frame = next;
    return data;
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

    ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking;
    window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
    window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus;

    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
    ImGui::Begin("DockspaceWindow", NULL, window_flags);
    ImGui::PopStyleVar(3);

    ImGuiID dockspace_id = ImGui::GetID("Dockspace");
    ImGuiDockNodeFlags dockspace_flags = ImGuiDockNodeFlags_PassthruDockspace;
    ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), dockspace_flags);

    ImGui::End();
}

static void BeginCanvas(const char* id) {
    // Invisible Canvas
    ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(viewport->Pos);
    ImGui::SetNextWindowSize(viewport->Size);
    ImGui::SetNextWindowViewport(viewport->ID);
    ImGui::SetNextWindowBgAlpha(0.0f);

    ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking;
    window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
    window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus | ImGuiWindowFlags_NoInputs;

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

void grow_mask_by_covalent_bond(Array<bool> mask, Array<const Bond> bonds, int32 extent) {
    Array<bool> prev_mask = {(bool*)TMP_MALLOC(mask.size_in_bytes()), mask.size()};
    defer { TMP_FREE(prev_mask.data()); };
    memcpy(prev_mask.data(), mask.data(), mask.size_in_bytes());

    for (int32 i = 0; i < extent; i++) {
        for (const auto& bond : bonds) {
            const int32 idx[2] = {bond.idx[0], bond.idx[1]};
            if (prev_mask[idx[0]] && !mask[idx[1]]) {
                mask[idx[1]] = true;
            } else if (prev_mask[idx[1]] && !mask[idx[0]]) {
                mask[idx[0]] = true;
            }
        }
        memcpy(prev_mask.data(), mask.data(), prev_mask.size_in_bytes());
    }
}

void grow_mask_by_radial_extent(Array<bool> mask, Array<const vec3> positions, float extent) {
    Array<bool> prev_mask = {(bool*)TMP_MALLOC(mask.size_in_bytes()), mask.size()};
    defer { TMP_FREE(prev_mask.data()); };
    memcpy(prev_mask.data(), mask.data(), mask.size_in_bytes());

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

void expand_mask_to_residue(Array<bool> mask, Array<const Residue> residues) {
    for (const auto& res : residues) {
        for (int64 i = res.atom_idx.beg; i != res.atom_idx.end; i++) {
            if (mask[i]) {
                memset_array(mask.subarray(res.atom_idx.beg, res.atom_idx.end - res.atom_idx.beg), true);
                break;
            }
        }
    }
}

void expand_mask_to_chain(Array<bool> mask, Array<const Chain> chains) {
    for (const auto& chain : chains) {
        for (int64 i = chain.atom_idx.beg; i != chain.atom_idx.end; i++) {
            if (mask[i]) {
                memset_array(mask.subarray(chain.atom_idx.beg, chain.atom_idx.end - chain.atom_idx.beg), true);
                break;
            }
        }
    }
}

void apply_mask(Array<bool> dst_mask, Array<const bool> mask_a, Array<const bool> mask_b, SelectionOperator op) {
    ASSERT(dst_mask.size() == mask_a.size());
    ASSERT(dst_mask.size() == mask_b.size());

    switch (op) {
        case SelectionOperator::Or:
            for (int64 i = 0; i < dst_mask.size(); i++) {
                dst_mask[i] = mask_a[i] || mask_b[i];
            }
            break;
        case SelectionOperator::And:
            for (int64 i = 0; i < dst_mask.size(); i++) {
                dst_mask[i] = mask_a[i] && mask_b[i];
            }
            break;
        default:
            ASSERT(false);
    }
}

void invert_mask(Array<bool> mask) {
    for (auto& v : mask) {
        v = !v;
    }
}

// ### DRAW WINDOWS ###
static void draw_main_menu(ApplicationData* data) {
    ASSERT(data);
    bool new_clicked = false;

    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("New", "CTRL+N")) {
                new_clicked = true;
            }
            if (ImGui::MenuItem("Load Data", "CTRL+L")) {
                auto res = platform::file_dialog(platform::FileDialogFlags_Open, {}, "pdb,gro,xtc");
                if (res.result == platform::FileDialogResult::FILE_OK) {
                    load_molecule_data(data, res.path);
                    if (data->representations.buffer.size() > 0) {
                        reset_representations(data);
                    } else {
                        create_representation(data);
                    }
                    stats::clear_all_properties();
                    reset_view(data, true);
                }
            }
            if (ImGui::MenuItem("Open", "CTRL+O")) {
                auto res = platform::file_dialog(platform::FileDialogFlags_Open, {}, FILE_EXTENSION);
                if (res.result == platform::FileDialogResult::FILE_OK) {
                    load_workspace(data, res.path);
                }
            }
            if (ImGui::MenuItem("Save", "CTRL+S")) {
                if (!data->files.workspace) {
                    auto res = platform::file_dialog(platform::FileDialogFlags_Save, {}, FILE_EXTENSION);
                    if (res.result == platform::FileDialogResult::FILE_OK) {
                        if (!get_file_extension(res.path)) {
                            snprintf(res.path.cstr() + strnlen(res.path.cstr(), res.path.capacity()), res.path.capacity(), ".%s", FILE_EXTENSION);
                        }
                        save_workspace(data, res.path);
                    }
                } else {
                    save_workspace(data, data->files.workspace);
                }
            }
            if (ImGui::MenuItem("Save As")) {
                auto res = platform::file_dialog(platform::FileDialogFlags_Save, {}, FILE_EXTENSION);
                if (res.result == platform::FileDialogResult::FILE_OK) {
                    if (!get_file_extension(res.path)) {
                        snprintf(res.path.cstr() + strnlen(res.path.cstr(), res.path.capacity()), res.path.capacity(), ".%s", FILE_EXTENSION);
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
            ImGui::Text("Background");
            ImGui::ColorEdit3("Color", &data->visuals.background.color[0], ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_Float);
            ImGui::SliderFloat("##Intensity", &data->visuals.background.intensity, 0.f, 100.f);
            ImGui::EndGroup();
            ImGui::Separator();
            // Temporal
            ImGui::BeginGroup();
            {
                ImGui::Checkbox("Temporal Effects", &data->visuals.temporal_reprojection.enabled);
                if (data->visuals.temporal_reprojection.enabled) {
                    ImGui::Checkbox("Jitter Samples", &data->visuals.temporal_reprojection.jitter);
                    ImGui::SliderFloat("Feedback Min", &data->visuals.temporal_reprojection.feedback_min, 0.5f, 1.0f);
                    ImGui::SliderFloat("Feedback Max", &data->visuals.temporal_reprojection.feedback_max, 0.5f, 1.0f);
                    ImGui::Checkbox("Motion Blur", &data->visuals.temporal_reprojection.motion_blur.enabled);
                    if (data->visuals.temporal_reprojection.motion_blur.enabled) {
                        ImGui::SliderFloat("Motion Scale", &data->visuals.temporal_reprojection.motion_blur.motion_scale, 0.f, 1.0f);
                    }
                }
            }
            ImGui::EndGroup();
            ImGui::Separator();

            // SSAO
            ImGui::BeginGroup();
            ImGui::Checkbox("SSAO", &data->visuals.ssao.enabled);
            if (data->visuals.ssao.enabled) {
                ImGui::SliderFloat("Intensity", &data->visuals.ssao.intensity, 0.5f, 12.f);
                ImGui::SliderFloat("Radius", &data->visuals.ssao.radius, 1.f, 30.f);
                ImGui::SliderFloat("Bias", &data->visuals.ssao.bias, 0.0f, 1.0f);
            }
            ImGui::EndGroup();
            ImGui::Separator();

            // DOF
            ImGui::BeginGroup();
            ImGui::Checkbox("Depth of Field", &data->visuals.dof.enabled);
            if (data->visuals.dof.enabled) {
                ImGui::SliderFloat("Focus Point", &data->visuals.dof.focus_depth, 0.001f, 200.f);
                ImGui::SliderFloat("Focus Scale", &data->visuals.dof.focus_scale, 0.001f, 100.f);
            }
            ImGui::EndGroup();
            ImGui::Separator();

            // Tonemapping
            ImGui::BeginGroup();
            ImGui::Checkbox("Tonemapping", &data->visuals.tonemapping.enabled);
            if (data->visuals.tonemapping.enabled) {
                ImGui::Combo("Function", &data->visuals.tonemapping.tonemapper, "Passthrough\0Exposure Gamma\0Filmic\0\0");
                ImGui::SliderFloat("Exposure", &data->visuals.tonemapping.exposure, 0.01f, 10.f);
                ImGui::SliderFloat("Gamma", &data->visuals.tonemapping.gamma, 1.0f, 3.0f);
            }
            ImGui::EndGroup();
            ImGui::Separator();

            ImGui::BeginGroup();
            ImGui::PushID("Selection");
            ImGui::Text("Selection");
            ImGui::ColorEdit3("Color", &data->selection.selection_color[0], ImGuiColorEditFlags_NoInputs);
            ImGui::SameLine();
            ImGui::SliderFloat("Scale", &data->selection.selection_scale, 0.f, 20.f);
            ImGui::PopID();
            ImGui::PushID("Highlight");
            ImGui::Text("Highlight");
            ImGui::ColorEdit3("Color", &data->selection.highlight_color[0], ImGuiColorEditFlags_NoInputs);
            ImGui::SameLine();
            ImGui::SliderFloat("Scale", &data->selection.highlight_scale, 0.f, 20.f);
            ImGui::PopID();
            ImGui::PushID("Outline");
            ImGui::Text("Outline");
            ImGui::ColorEdit3("Color", &data->selection.outline_color[0], ImGuiColorEditFlags_NoInputs);
            ImGui::SameLine();
            ImGui::SliderFloat("Scale", &data->selection.outline_scale, 0.f, 20.f);
            ImGui::PopID();
            ImGui::EndGroup();
            ImGui::Separator();

            ImGui::BeginGroup();
            ImGui::Checkbox("Draw Control Points", &data->visuals.spline.draw_control_points);
            ImGui::Checkbox("Draw Spline", &data->visuals.spline.draw_spline);
            ImGui::EndGroup();
            ImGui::Separator();

            // Property Overlay
            ImGui::BeginGroup();
            ImGui::Text("Property Style");
            auto style = stats::get_style();
            ImGui::Text("point colors ");
            for (int32 i = 0; i < style->NUM_COLORS; i++) {
                ImGui::SameLine();
                ImGui::PushID(i);
                ImVec4 color = ImColor(style->point_colors[i]);
                if (ImGui::ColorEdit4("PointColor", (float*)&color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel)) style->point_colors[i] = ImColor(color);
                ImGui::PopID();
            }
            ImGui::Text("line color   ");
            ImGui::SameLine();
            ImVec4 color = ImColor(style->line_color);
            if (ImGui::ColorEdit4("LineColor", (float*)&color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel)) style->line_color = ImColor(color);
            ImGui::EndGroup();
            ImGui::Separator();

            ImGui::BeginGroup();
            ImGui::Checkbox("Hydrogen Bond", &data->hydrogen_bonds.enabled);
            if (data->hydrogen_bonds.enabled) {
                ImGui::PushID("hydrogen_bond");
                if (ImGui::SliderFloat("Distance Cutoff", &data->hydrogen_bonds.distance_cutoff, HYDROGEN_BOND_DISTANCE_CUTOFF_MIN, HYDROGEN_BOND_DISTANCE_CUTOFF_MAX)) {
                    data->hydrogen_bonds.dirty = true;
                }
                if (ImGui::SliderFloat("Angle Cutoff", &data->hydrogen_bonds.angle_cutoff, HYDROGEN_BOND_ANGLE_CUTOFF_MIN, HYDROGEN_BOND_ANGLE_CUTOFF_MAX)) {
                    data->hydrogen_bonds.dirty = true;
                }
                ImGui::Checkbox("Overlay", &data->hydrogen_bonds.overlay);
                ImGui::ColorEdit4("Color", (float*)&data->hydrogen_bonds.color, ImGuiColorEditFlags_NoInputs);
                ImGui::PopID();
            }
            ImGui::EndGroup();
            ImGui::Separator();

            ImGui::BeginGroup();
            ImGui::Checkbox("Simulation Box", &data->simulation_box.enabled);
            if (data->simulation_box.enabled) {
                ImGui::PushID("simulation_box");
                ImGui::ColorEdit4("Color", (float*)&data->simulation_box.color, ImGuiColorEditFlags_NoInputs);
                ImGui::PopID();
            }
            ImGui::EndGroup();
            ImGui::Separator();

            ImGui::BeginGroup();
            if (ImGui::Checkbox("Density Volume", &data->density_volume.enabled)) {
                // if (data->density_volume.enabled) data->density_volume.texture.dirty = true;
            }
            if (data->density_volume.enabled) {
                ImGui::PushID("density_volume");
                ImGui::ColorEdit3("Color", (float*)&data->density_volume.color, ImGuiColorEditFlags_NoInputs);
                ImGui::SliderFloat("Scale", &data->density_volume.density_scale, 0.001f, 10.f, "%.3f", 3.f);
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
            ImGui::Checkbox("Selection", &data->selection.show_window);

            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Selection")) {
            static DynamicArray<bool> mask;
            mask.resize(data->mol_data.dynamic.molecule.atom.count);

            if (ImGui::MenuItem("Clear Selection")) {
                memset_array(data->selection.current_selection_mask, false);
                data->gpu_buffers.dirty.selection = true;
            }

            if (ImGui::MenuItem("Invert Selection")) {
                for (auto& v : data->selection.current_selection_mask) {
                    v = !v;
                }
                data->gpu_buffers.dirty.selection = true;
            }

            ImGui::Spacing();
            ImGui::Separator();

            // MODES
            ImGui::Combo("Level Mode", (int*)(&data->selection.level_mode), "Atom\0Residue\0Chain\0\0");
            ImGui::Combo("Operator Mode", (int*)(&data->selection.op_mode), "Or\0And\0\0");

            ImGui::Spacing();
            ImGui::Separator();

            // QUERY
            {
                static char buf[256] = {0};
                static bool query_ok = false;

                ImGui::Text("Query");
                const auto TEXT_BG_DEFAULT_COLOR = ImGui::ColorConvertFloat4ToU32(ImGui::GetStyle().Colors[ImGuiCol_FrameBg]);
                ImGui::PushStyleColor(ImGuiCol_FrameBg, query_ok ? TEXT_BG_DEFAULT_COLOR : TEXT_BG_ERROR_COLOR);
                const bool query_modified = ImGui::InputText("##query", buf, ARRAY_SIZE(buf), ImGuiInputTextFlags_AutoSelectAll);
                const bool pressed_enter = ImGui::IsItemActivePreviousFrame() && !ImGui::IsItemActive() && ImGui::IsKeyPressed(ImGui::GetIO().KeyMap[ImGuiKey_Enter]);
                ImGui::PopStyleColor();
                // ImGui::SameLine();
                if (!query_ok) ImGui::PushDisabled();
                const bool apply = ImGui::Button("Apply##query") || pressed_enter;
                if (!query_ok) ImGui::PopDisabled();

                // if (ImGui::IsWindowAppearing()) {
                //    ImGui::SetKeyboardFocusHere(-1);
                //}

                const bool show_preview =
                    (ImGui::GetFocusID() == ImGui::GetID("##query")) || (ImGui::GetHoveredID() == ImGui::GetID("##query")) || (ImGui::GetHoveredID() == ImGui::GetID("Apply##query"));

                if (show_preview) {
                    if (query_modified) {
                        query_ok = filter::compute_filter_mask(mask, data->mol_data.dynamic, buf);
                        if (!query_ok) {
                            memset_array(mask, false);
                        }
                        data->gpu_buffers.dirty.selection = true;
                    }

                    if (query_ok) {
                        for (int64 i = 0; i < mask.size(); i++) {
                            const bool mask_val = mask[i];
                            const bool curr_val = data->selection.current_selection_mask[i];

                            if (data->selection.op_mode == SelectionOperator::And) {
                                data->selection.current_highlight_mask[i] = curr_val & mask_val;
                            } else if (data->selection.op_mode == SelectionOperator::Or) {
                                data->selection.current_highlight_mask[i] = curr_val | mask_val;
                            }
                        }
                        data->gpu_buffers.dirty.selection = true;
                    }
                }

                if (apply) {
                    for (int64 i = 0; i < mask.size(); i++) {
                        const bool mask_val = mask[i];
                        const bool curr_val = data->selection.current_selection_mask[i];

                        if (data->selection.op_mode == SelectionOperator::And) {
                            data->selection.current_selection_mask[i] = curr_val & mask_val;
                        } else if (data->selection.op_mode == SelectionOperator::Or) {
                            data->selection.current_selection_mask[i] = curr_val | mask_val;
                        }
                    }
                    data->gpu_buffers.dirty.selection = true;
                }
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
                const bool extent_changed = ImGui::SliderFloat("Extent", &extent, 1.0f, 10.f);
                const bool apply = ImGui::Button("Apply##grow");
                const bool show_preview = mode_changed || extent_changed || (ImGui::GetHoveredID() == ImGui::GetID("Extent")) || (ImGui::GetHoveredID() == ImGui::GetID("Apply##grow")) ||
                                          (ImGui::GetHoveredID() == ImGui::GetID("Mode"));

                if (show_preview) {
                    memcpy(mask.data(), data->selection.current_selection_mask.data(), data->selection.current_selection_mask.size_in_bytes());

                    switch (data->selection.grow_mode) {
                        case SelectionGrowth::CovalentBond:
                            grow_mask_by_covalent_bond(mask, get_covalent_bonds(data->mol_data.dynamic.molecule), (int32)extent);
                            break;
                        case SelectionGrowth::Radial:
                            grow_mask_by_radial_extent(mask, get_positions(data->mol_data.dynamic.molecule), extent);
                            break;
                        default:
                            ASSERT(false);
                    }

                    switch (data->selection.level_mode) {
                        case SelectionLevel::Atom:
                            // No need to expand the mask
                            break;
                        case SelectionLevel::Residue:
                            expand_mask_to_residue(mask, data->mol_data.dynamic.molecule.residues);
                            break;
                        case SelectionLevel::Chain:
                            expand_mask_to_chain(mask, data->mol_data.dynamic.molecule.chains);
                            break;
                        default:
                            ASSERT(false);
                    }

                    memcpy(data->selection.current_highlight_mask.data(), mask.data(), mask.size_in_bytes());
                    data->gpu_buffers.dirty.selection = true;

                    if (apply) {
                        memcpy(data->selection.current_selection_mask.data(), mask.data(), mask.size_in_bytes());
                    }
                }
            }

            ImGui::Spacing();
            ImGui::Separator();

            // STORED SELECTIONS
            {
                ImGui::Text("Stored Selections");
                const bool disable_new = is_array_zero((Array<const bool>)(data->selection.current_selection_mask));
                if (disable_new) ImGui::PushDisabled();
                if (ImGui::Button("Create New")) {
                    char name_buf[64];
                    snprintf(name_buf, 64, "selection%i", (int)data->selection.stored_selections.size());
                    create_selection(data, name_buf, data->selection.current_selection_mask);
                }
                if (disable_new) ImGui::PopDisabled();
                for (int i = 0; i < data->selection.stored_selections.size(); i++) {
                    auto& sel = data->selection.stored_selections[i];
                    const float32 item_width = math::clamp(ImGui::GetWindowContentRegionWidth() - 90.f, 100.f, 300.f);
                    StringBuffer<128> name;
                    snprintf(name, name.size(), "%s###ID", sel.name.cstr());

                    ImGui::PushID(i);
                    if (ImGui::CollapsingHeader(name)) {
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
                        }
                    }

                    const auto h_id = ImGui::GetHoveredID();
                    bool show_preview = (h_id == ImGui::GetID(name) || h_id == ImGui::GetID("Activate") || h_id == ImGui::GetID("Remove") || h_id == ImGui::GetID("Clone"));

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
            const float32 ms = compute_avg_ms(data->ctx.timing.delta_s);
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
    mask.resize(data->mol_data.dynamic.molecule.atom.count);

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
        bool pressed_enter = ImGui::IsItemActivePreviousFrame() && !ImGui::IsItemActive() && ImGui::IsKeyPressed(ImGui::GetIO().KeyMap[ImGuiKey_Enter]);
        ImGui::PopStyleColor();

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
                query_ok = filter::compute_filter_mask(mask, data->mol_data.dynamic, buf);
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
        bool show = (ImGui::GetHoveredID() == ImGui::GetID("Extent") || ImGui::GetActiveID() == ImGui::GetID("Extent") || ImGui::GetHoveredID() == ImGui::GetID("Apply"));

        if (show) {
            Array<bool> prev_mask = {(bool*)TMP_MALLOC(mask.size_in_bytes()), mask.size()};
            defer { TMP_FREE(prev_mask.data()); };

            memset_array(mask, false);
            memcpy(prev_mask.data(), data->selection.current_selection.data(), data->selection.current_selection.size_in_bytes());

            if (grow_mode == CovalentBond) {
                for (int64 i = 0; i < (int64)extent; i++) {
                    for (const auto& bond : data->mol_data.dynamic.molecule.covalent_bonds) {
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
                const auto positions = get_positions(data->mol_data.dynamic.molecule);
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
                    expand_mask_to_residue(mask, data->mol_data.dynamic.molecule.residues);
                    break;
                case SelectionLevel::Chain:
                    expand_mask_to_chain(mask, data->mol_data.dynamic.molecule.chains);
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
        bool show = (h_id == ImGui::GetID(name) || h_id == ImGui::GetID("Activate") || h_id == ImGui::GetID("Remove") || h_id == ImGui::GetID("Clone"));

        ImGui::PopID();

        if (show) {
            memcpy(data->selection.current_highlight.data(), sel.atom_mask.data(), sel.atom_mask.size_in_bytes());
            data->gpu_buffers.dirty.selection = true;
        }
    }

    ImGui::End();
}
*/

void draw_context_popup(ApplicationData* data) {
    ASSERT(data);

    const bool shift_down = data->ctx.input.key.down[Key::KEY_LEFT_SHIFT] || data->ctx.input.key.down[Key::KEY_RIGHT_SHIFT];
    if (data->ctx.input.mouse.clicked[1] && !shift_down && !ImGui::GetIO().WantTextInput) {
        if (data->selection.right_clicked != -1 && data->mol_data.dynamic) {
            ImGui::OpenPopup("AtomContextPopup");
        }
    }

    if (ImGui::BeginPopup("AtomContextPopup")) {
        if (data->selection.right_clicked != -1 && data->mol_data.dynamic) {
            if (ImGui::MenuItem("Recenter Trajectory")) {
                recenter_trajectory(&data->mol_data.dynamic, data->mol_data.dynamic.molecule.atom.residue_indices[data->selection.right_clicked]);
                interpolate_atomic_positions(get_positions(data->mol_data.dynamic.molecule), data->mol_data.dynamic.trajectory, data->time, data->interpolation);
                data->gpu_buffers.dirty.position = true;
                ImGui::CloseCurrentPopup();
            }
        }
        ImGui::EndPopup();
    }
}

static void draw_animation_control_window(ApplicationData* data) {
    ASSERT(data);
    if (!data->mol_data.dynamic.trajectory) return;

    ImGui::Begin("Control");
    int32 num_frames = data->mol_data.dynamic.trajectory.num_frames;
    ImGui::Text("Num Frames: %i", num_frames);
    float32 t = (float)data->time;
    if (ImGui::SliderFloat("Time", &t, 0, (float)(num_frames - 1))) {
        data->time = t;
    }
    ImGui::SliderFloat("fps", &data->frames_per_second, 0.1f, 100.f, "%.3f", 4.f);
    ImGui::Combo("type", (int*)(&data->interpolation), "Nearest\0Linear\0Linear Periodic\0Cubic\0Cubic Periodic\0\0");
    if (data->is_playing) {
        if (ImGui::Button("Pause")) data->is_playing = false;
    } else {
        if (ImGui::Button("Play")) data->is_playing = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("Stop")) {
        data->is_playing = false;
        data->time = 0.0;
    }
    ImGui::End();
}

static void draw_representations_window(ApplicationData* data) {
    const auto old_hash = hash::crc64(data->representations.buffer);

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
        const float32 item_width = math::clamp(ImGui::GetWindowContentRegionWidth() - 90.f, 100.f, 300.f);
        StringBuffer<128> name;
        snprintf(name, name.size(), "%s###ID", rep.name.buffer);

        ImGui::PushID(i);
        if (ImGui::CollapsingHeader(name)) {
            ImGui::Checkbox("enabled", &rep.enabled);
            ImGui::SameLine();
            if (ImGui::DeleteButton("remove")) {
                remove_representation(data, i);
            }
            ImGui::SameLine();
            if (ImGui::Button("clone")) {
                clone_representation(data, rep);
            }

            ImGui::PushItemWidth(item_width);
            ImGui::InputText("name", rep.name, rep.name.capacity());
            if (!rep.filter_is_ok) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
            if (ImGui::InputText("filter", rep.filter, rep.filter.capacity(), ImGuiInputTextFlags_EnterReturnsTrue)) {
                recompute_colors = true;
            }
            if (!rep.filter_is_ok) ImGui::PopStyleColor();
            ImGui::Combo("type", (int*)(&rep.type), "VDW\0Licorice\0Ball & Stick\0Ribbons\0Cartoon\0\0");
            if (ImGui::Combo("color mapping", (int*)(&rep.color_mapping), "Static Color\0CPK\0Res Id\0Res Idx\0Chain Id\0Chain Idx\0Secondary Structure\0\0")) {
                recompute_colors = true;
            }
            ImGui::PopItemWidth();
            if (rep.color_mapping == ColorMapping::Static) {
                ImGui::SameLine();
                if (ImGui::ColorEdit4("color", (float*)&rep.static_color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel)) {
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

    const auto new_hash = hash::crc64(data->representations.buffer);
    data->representations.changed = (new_hash != old_hash);
}

static void draw_property_window(ApplicationData* data) {
    static bool first_time_shown = true;
    ImGui::Begin("Properties", &data->statistics.show_property_window, ImGuiWindowFlags_NoFocusOnAppearing);

    ImGui::PushID("PROPERTIES");
    ImGui::PushItemWidth(-1);
    // ImGui::Columns(4, "columns", true);
    ImGui::BeginColumns("columns", 4, ImGuiColumnsFlags_NoPreserveWidths);
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
    for (int32 i = 0; i < (int32)properties.count; i++) {
        auto prop = properties[i];

        ImGui::Separator();
        ImGui::PushID(i);

        ImGui::PushItemWidth(-1);
        if (!prop->valid) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
        if (ImGui::InputText("##name", prop->name_buf, prop->name_buf.size(), ImGuiInputTextFlags_EnterReturnsTrue)) {
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
            if (data->mol_data.dynamic.molecule) {
                const int32 atom_idx = data->selection.right_clicked;
                if (atom_idx != -1) {
                    ASSERT(atom_idx < data->mol_data.dynamic.molecule.atom.count);
                    const int32 residue_idx = data->mol_data.dynamic.molecule.atom.residue_indices[data->selection.right_clicked];
                    const int32 chain_idx = data->mol_data.dynamic.molecule.residues[residue_idx].chain_idx;

                    char buf[buf_len];

                    snprintf(buf, buf_len, "atom(%i) ", data->selection.right_clicked + 1);
                    if (ImGui::MenuItem(buf)) {
                        memcpy(key_buf, buf, buf_len);
                        paste_buf = true;
                    }

                    if (residue_idx > -1) {
                        const auto& residue = data->mol_data.dynamic.molecule.residues[residue_idx];

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
                            snprintf(buf, buf_len, "resatom(resid(%i), %i) ", residue.id, atom_idx + 1);
                            if (ImGui::MenuItem(buf)) {
                                memcpy(key_buf, buf, buf_len);
                                paste_buf = true;
                            }
                            snprintf(buf, buf_len, "resatom(resname(%s), %i) ", residue.name.cstr(), atom_idx + 1);
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
        if (ImGui::InputText("##args", prop->args_buf, prop->args_buf.size(), ImGuiInputTextFlags_EnterReturnsTrue)) {
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
    // stats::compute_stats(data->mol_data.dynamic);
    // compute_statistics_async(data);
    //}
}

static void draw_atom_info_window(const MoleculeStructure& mol, int atom_idx, int x, int y) {

    // @TODO: Assert things and make this failproof
    if (atom_idx < 0 || atom_idx >= mol.atom.count) return;

    int res_idx = mol.atom.residue_indices[atom_idx];
    const Residue& res = mol.residues[res_idx];
    const char* res_id = res.name;
    int local_idx = atom_idx - res.atom_idx.beg;
    const vec3& pos = mol.atom.positions[atom_idx];
    const char* label = mol.atom.labels[atom_idx];
    const char* elem = element::name(mol.atom.elements[atom_idx]);
    const char* symbol = element::symbol(mol.atom.elements[atom_idx]);

    int chain_idx = res.chain_idx;
    const char* chain_id = "\0";
    if (chain_idx != -1 && mol.chains.size() > 0) {
        const Chain& chain = mol.chains[chain_idx];
        chain_id = chain.id;
        chain_idx = res.chain_idx;
    }

    // External indices begin with 1 not 0
    res_idx += 1;
    chain_idx += 1;
    atom_idx += 1;
    local_idx += 1;

    char buff[256];
    int len = 0;
    len += snprintf(buff, 256, "atom[%i][%i]: %s %s %s (%.2f, %.2f, %.2f)\n", atom_idx, local_idx, label, elem, symbol, pos.x, pos.y, pos.z);
    len += snprintf(buff + len, 256 - len, "res[%i]: %s\n", res_idx, res_id);
    if (chain_idx) {
        len += snprintf(buff + len, 256 - len, "chain[%i]: %s\n", chain_idx, chain_id);
    }

    if (res_idx < mol.backbone.angles.size() && res_idx < mol.backbone.segments.size() && valid_segment(mol.backbone.segments[res_idx])) {
        const auto angles = mol.backbone.angles[res_idx] * math::RAD_TO_DEG;
        len += snprintf(buff + len, 256 - len, u8"\u03C6: %.1f\u00b0, \u03C8: %.1f\u00b0\n", angles.x, angles.y);
    }

    ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(x + 10.f, y + 18.f) + viewport->Pos);
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
    ImGui::Begin("##Atom Info", 0, ImGuiWindowFlags_Tooltip | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoDocking);
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

static void draw_async_info(ApplicationData* data) {
    constexpr float32 WIDTH = 300.f;
    constexpr float32 MARGIN = 10.f;
    constexpr float32 PROGRESS_FRACT = 0.3f;

    float32 traj_fract = data->async.trajectory.fraction;
    float32 angle_fract = data->async.backbone_angles.fraction;
    float32 stats_fract = stats::fraction_done();

    if ((0.f < traj_fract && traj_fract < 1.f) || (0.f < angle_fract && angle_fract < 1.f) || (0.f < stats_fract && stats_fract < 1.f)) {

        ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(viewport->Pos + ImVec2(data->ctx.window.width - WIDTH - MARGIN, ImGui::GetCurrentContext()->FontBaseSize + ImGui::GetStyle().FramePadding.y * 2.f + MARGIN));
        ImGui::SetNextWindowSize(ImVec2(WIDTH, 0));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
        ImGui::Begin("##Async Info", 0,
                     ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoSavedSettings |
                         ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoFocusOnAppearing);

        char buf[32];
        if (0.f < traj_fract && traj_fract < 1.f) {
            snprintf(buf, 32, "%.1f%%", traj_fract * 100.f);
            ImGui::ProgressBar(traj_fract, ImVec2(ImGui::GetWindowContentRegionWidth() * PROGRESS_FRACT, 0), buf);
            ImGui::SameLine();
            ImGui::Text("Reading Trajectory");
            ImGui::SameLine();
            if (ImGui::Button("X")) {
                data->async.trajectory.sync.signal_stop_and_wait();
                compute_backbone_angles_async(data);
                // compute_statistics_async(data);
                data->async.trajectory.fraction = 0.f;
            }
        }
        if (0.f < angle_fract && angle_fract < 1.f) {
            snprintf(buf, 32, "%.1f%%", angle_fract * 100.f);
            ImGui::ProgressBar(angle_fract, ImVec2(ImGui::GetWindowContentRegionWidth() * PROGRESS_FRACT, 0), buf);
            ImGui::SameLine();
            ImGui::Text("Computing Backbone Angles");
        }
        if (0.f < stats_fract && stats_fract < 1.f) {
            snprintf(buf, 32, "%.1f%%", stats_fract * 100.f);
            ImGui::ProgressBar(stats_fract, ImVec2(ImGui::GetWindowContentRegionWidth() * PROGRESS_FRACT, 0), buf);
            ImGui::SameLine();
            ImGui::Text("Computing Statistics");
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
        const float num_frames = (float)data->mol_data.dynamic.trajectory.num_frames;

        ImGui::Checkbox("Dynamic Framewindow", &data->time_filter.dynamic_window);
        if (data->time_filter.dynamic_window) {
            ImGui::SameLine();
            ImGui::SliderFloat("Window Extent", &data->time_filter.window_extent, 1.f, num_frames);
        }
        ImGui::BeginChild("Scroll Region", ImVec2(0, 0), true, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_HorizontalScrollbar);

        const Range<float> frame_range = {0, num_frames};
        auto old_range = data->time_filter.range;

        ImGui::PushItemWidth(ImGui::GetWindowContentRegionWidth() * zoom);
        if (ImGui::RangeSliderFloat("###selection_range", &data->time_filter.range.beg, &data->time_filter.range.end, 0.f, num_frames)) {
            if (data->time_filter.dynamic_window) {
                if (data->time_filter.range.x != old_range.x && data->time_filter.range.y != old_range.y) {
                    data->time = math::lerp(data->time_filter.range.x, data->time_filter.range.y, 0.5f);
                } else {
                    if (data->time_filter.range.x != old_range.x) {
                        data->time_filter.window_extent = 2.f * math::abs((float)data->time - data->time_filter.range.x);
                    } else if (data->time_filter.range.y != old_range.y) {
                        data->time_filter.window_extent = 2.f * math::abs((float)data->time - data->time_filter.range.y);
                    }
                }
            }
        }

        // const int32 prop_count = stats::get_property_count();
        // const float32 plot_height = ImGui::GetContentRegionAvail().y / (float)prop_count;
        const float32 plot_height = 100.f;
        const uint32 bar_fill_color = ImColor(1.f, 1.f, 1.f, 0.25f);
        const uint32 var_fill_color = ImColor(1.f, 1.f, 0.3f, 0.1f);
        const uint32 var_line_color = ImColor(1.f, 1.f, 0.3f, 0.3f);
        const uint32 var_text_color = ImColor(1.f, 1.f, 0.3f, 0.5f);

        static float32 selection_start;
        static bool is_selecting = false;

        const auto properties = stats::get_properties();
        for (int i = 0; i < (int32)properties.count; i++) {
            auto prop = properties[i];

            if (!prop->enable_timeline) continue;
            Array<float> prop_data = prop->avg_data;
            CString prop_name = prop->name_buf;
            auto prop_range = prop->avg_data_range;
            if (!prop_data) continue;
            float32 pad = math::abs(prop_range.y - prop_range.x) * 0.75f;
            Range<float> display_range = {prop_range.x - pad, prop_range.y + pad};
            if (display_range.x == display_range.y) {
                display_range.x -= 1.f;
                display_range.y += 1.f;
            }
            // float32 val = (float)time;
            ImGuiID id = ImGui::GetID(prop_name);

            ImGui::PushID(i);

            ImGui::BeginPlot(prop_name, ImVec2(0, plot_height), ImVec2(frame_range.x, frame_range.y), ImVec2(display_range.x, display_range.y), ImGui::LinePlotFlags_AxisX);
            const ImRect inner_bb(ImGui::GetItemRectMin() + ImGui::GetStyle().FramePadding, ImGui::GetItemRectMax() - ImGui::GetStyle().FramePadding);
            ImGui::PushClipRect(ImGui::GetItemRectMin(), ImGui::GetItemRectMax(), true);

            if (ImGui::IsItemHovered()) ImGui::SetHoveredID(id);

            ImGui::PlotVerticalBars(prop->filter_fraction.ptr, (int32)prop->filter_fraction.count, bar_fill_color);
            if (prop->std_dev_data.ptr[0] > 0.f) {
                ImGui::PlotVariance(prop->avg_data.ptr, prop->std_dev_data.ptr, (int32)prop->std_dev_data.count, 1.f, var_line_color, var_fill_color);
            }
            ImGui::PlotValues(prop->name_buf.cstr(), prop_data.ptr, (int32)prop_data.count);

            ImGui::PopClipRect();

            if (ImGui::IsItemHovered() && ImGui::GetIO().MouseClicked[0]) {
                ImGui::SetActiveID(id, ImGui::GetCurrentWindow());
            }

            if (ImGui::IsItemHovered() && ImGui::GetIO().MouseClicked[1] && ImGui::GetIO().KeyCtrl) {
                data->time_filter.range = frame_range;
            }

            if (ImGui::GetActiveID() == id) {
                if (ImGui::GetIO().MouseClicked[0] && ImGui::GetIO().KeyCtrl) {
                    float32 t = (ImGui::GetIO().MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x);
                    selection_start = ImLerp(frame_range.x, frame_range.y, t);
                    data->time_filter.range.x = selection_start;
                    data->time_filter.range.y = selection_start;
                    is_selecting = true;
                } else if (is_selecting) {
                    float32 t = (ImGui::GetIO().MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x);
                    float32 v = ImLerp(frame_range.x, frame_range.y, t);
                    if (v < data->time_filter.range.x) {
                        data->time_filter.range.x = v;
                    } else if (v > data->time_filter.range.x && v < data->time_filter.range.y) {
                        if (selection_start < v) {
                            data->time_filter.range.y = v;
                        } else {
                            data->time_filter.range.x = v;
                        }
                    } else if (v > data->time_filter.range.y) {
                        data->time_filter.range.y = v;
                    }
                } else if (ImGui::GetIO().MouseDown[0]) {
                    float32 t = ImClamp((ImGui::GetIO().MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x), 0.f, 1.f);
                    data->time = ImLerp(frame_range.x, frame_range.y, t);
                }

                if (!ImGui::GetIO().MouseDown[0] && !ImGui::IsItemHovered()) {
                    ImGui::ClearActiveID();
                    is_selecting = false;
                }
            }

            data->time_filter.range.x = ImClamp(data->time_filter.range.x, frame_range.x, frame_range.y);
            data->time_filter.range.y = ImClamp(data->time_filter.range.y, frame_range.x, frame_range.y);

            // SELECTION RANGE
            {
                constexpr ImU32 SELECTION_RANGE_COLOR = 0x55bbbbbb;
                const float32 t0 = (data->time_filter.range.x - frame_range.x) / (frame_range.y - frame_range.x);
                const float32 t1 = (data->time_filter.range.y - frame_range.x) / (frame_range.y - frame_range.x);
                const ImVec2 pos0 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t0, 0));
                const ImVec2 pos1 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t1, 1));
                ImGui::GetCurrentWindow()->DrawList->AddRectFilled(pos0, pos1, SELECTION_RANGE_COLOR);
            }

            // CURRENT FRAME POSITION
            {
                constexpr ImU32 CURRENT_LINE_COLOR = 0xaa33ffff;
                const float32 t = ((float)data->time - frame_range.x) / (frame_range.y - frame_range.x);
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
                const float32 min_x = ImGui::GetItemRectMin().x;
                const float32 max_x = ImGui::GetItemRectMax().x;
                float32 t = ImClamp((ImGui::GetIO().MousePos.x - min_x) / (max_x - min_x), 0.f, 1.f);
                int idx = ImClamp((int32)ImLerp(frame_range.x, frame_range.y, t), 0, (int32)prop->avg_data.count - 1);

                ImGui::BeginTooltip();
                ImGui::Text("%i: %g ", idx, prop->avg_data[idx]);
                ImGui::SameLine();
                ImGui::TextColored(ImColor(var_text_color), "(%g)", prop->std_dev_data[idx]);
                ImGui::EndTooltip();
            }

            ImGui::EndPlot();
            ImGui::PopID();
        }
        ImGui::PopItemWidth();

        if (data->time_filter.range != old_range) {
            stats::set_all_property_flags(false, true);
        }

        if (ImGui::IsWindowHovered() && ImGui::GetIO().MouseWheel != 0.f && ImGui::GetIO().KeyCtrl) {
            constexpr float32 ZOOM_SCL = 0.24f;
            float32 pre_coord = ImGui::GetScrollX() + (ImGui::GetIO().MousePos.x - ImGui::GetWindowPos().x) * zoom;
            zoom = math::clamp(zoom + ZOOM_SCL * ImGui::GetIO().MouseWheel, 1.f, 100.f);
            float32 post_coord = ImGui::GetScrollX() + (ImGui::GetIO().MousePos.x - ImGui::GetWindowPos().x) * zoom;
            float32 delta = pre_coord - post_coord;
            ImGui::SetScrollX(ImGui::GetScrollX() - delta);
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
    const float32 plot_height = 100.f;

    ImVec2 frame_size{ImGui::CalcItemWidth(), plot_height};

    constexpr uint32 FULL_FILL_COLOR = 0x99cc9e66;
    constexpr uint32 FULL_LINE_COLOR = 0xffcc9e66;
    constexpr uint32 FULL_TEXT_COLOR = 0xffcc9e66;

    constexpr uint32 FILT_FILL_COLOR = 0x3333ffff;
    constexpr uint32 FILT_LINE_COLOR = 0xaa33ffff;
    constexpr uint32 FILT_TEXT_COLOR = 0xaa33ffff;
    constexpr ImU32 SELECTION_RANGE_COLOR = 0x55bbbbbb;

    for (int32 i = 0; i < (int32)properties.count; i++) {
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
            const float32 max_val = prop->full_histogram.bin_range.y * 1.5f;
            ImGui::DrawFilledLine(inner_bb.Min, inner_bb.Max, prop->full_histogram.bins.ptr, (int32)prop->full_histogram.bins.count, max_val, FULL_LINE_COLOR, FULL_FILL_COLOR);

            ImGui::DrawFilledLine(inner_bb.Min, inner_bb.Max, prop->filt_histogram.bins.ptr, (int32)prop->filt_histogram.bins.count, max_val, FILT_LINE_COLOR, FILT_FILL_COLOR);
            // ImGui::PopClipRect();

            // SELECTION RANGE
            {
                const float32 t0 = (prop->filter.x - prop->total_data_range.x) / (prop->total_data_range.y - prop->total_data_range.x);
                const float32 t1 = (prop->filter.y - prop->total_data_range.x) / (prop->total_data_range.y - prop->total_data_range.x);
                const ImVec2 pos0 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t0, 0));
                const ImVec2 pos1 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t1, 1));
                ImGui::GetCurrentWindow()->DrawList->AddRectFilled(pos0, pos1, SELECTION_RANGE_COLOR);
            }

            ImGui::PopClipRect();

            if (ImGui::IsItemHovered()) {
                window->DrawList->AddLine(ImVec2(ImGui::GetIO().MousePos.x, inner_bb.Min.y), ImVec2(ImGui::GetIO().MousePos.x, inner_bb.Max.y), 0xffffffff);
                float32 t = (ImGui::GetIO().MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x);
                int32 count = (int32)prop->full_histogram.bins.count;
                int32 idx = ImClamp((int32)(t * (count - 1)), 0, count - 1);
                float32 full_val = prop->full_histogram.bins.ptr[idx];
                float32 filt_val = prop->filt_histogram.bins.ptr[idx];
                ImGui::BeginTooltip();
                ImGui::Text("%.3f:", ImLerp(prop->filt_histogram.value_range.x, prop->filt_histogram.value_range.y, t));
                ImGui::TextColored(ImColor(FULL_TEXT_COLOR), "%g", full_val * 100.f);
                ImGui::TextColored(ImColor(FILT_TEXT_COLOR), "%g", filt_val * 100.f);
                ImGui::EndTooltip();
            }

            if (ImGui::RangeSliderFloat("##filter", &prop->filter.x, &prop->filter.y, prop->total_data_range.x, prop->total_data_range.y)) {
                prop->filter_dirty = true;
            }
        }
        ImGui::PopID();
    }
    ImGui::PopItemWidth();
    ImGui::End();
}

static void draw_ramachandran_window(ApplicationData* data) {
    // const int32 num_frames = data->mol_data.dynamic.trajectory ? data->mol_data.dynamic.trajectory.num_frames : 0;
    // const int32 frame = (int32)data->time;
    const IntRange frame_range = {(int32)data->time_filter.range.x, (int32)data->time_filter.range.y};
    const auto& mol = data->mol_data.dynamic.molecule;
    Array<const BackboneAngle> trajectory_angles = get_backbone_angles(data->ramachandran.backbone_angles, frame_range.x, frame_range.y - frame_range.x);
    Array<const BackboneAngle> current_angles = mol.backbone.angles;
    Array<const BackboneSegment> backbone_segments = mol.backbone.segments;
    Array<const Residue> residues = mol.residues;
    Array<const bool> atom_selection = data->selection.current_selection_mask;
    Array<const bool> atom_highlight = data->selection.current_highlight_mask;

    ImGui::SetNextWindowSizeConstraints(ImVec2(400, 200), ImVec2(10000, 10000));
    ImGui::Begin("Ramachandran", &data->ramachandran.show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoScrollbar);

    ImGui::BeginColumns("cols", 2, ImGuiColumnsFlags_NoResize);
    {
        ImGui::Checkbox("Show Current Frame", &data->ramachandran.current.enabled);
        if (data->ramachandran.current.enabled) {
            ImGui::PushID("current");
            ImGui::SliderFloat("##base_radius", &data->ramachandran.current.radius, 0.5f, 5.f, "base radius %1.1f");
            // ImGui::SameLine();
            // ImGui::ColorEdit4("base color", (float*)&data->ramachandran.current.fill_color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);

            ImGui::SliderFloat("##selected_radius", &data->ramachandran.selected.radius, 0.5, 5.f, "selected radius %1.1f");
            // ImGui::SameLine();
            // ImGui::ColorEdit4("selected color", (float*)&data->ramachandran.selected.fill_color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);

            data->ramachandran.current.radius = math::round(data->ramachandran.current.radius * 2.f) / 2.f;
            data->ramachandran.selected.radius = math::round(data->ramachandran.selected.radius * 2.f) / 2.f;
            // ImGui::ColorEdit4("border color", (float*)&data->ramachandran.current.border_color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);
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
            ImGui::RangeSliderFloat("range", &data->time_filter.range.min, &data->time_filter.range.max, 0.f, (float)data->mol_data.dynamic.trajectory.num_frames, "(%.1f, %.1f)");
            ImGui::PopID();

            ramachandran::clear_accumulation_texture();
            ramachandran::compute_accumulation_texture(trajectory_angles, data->ramachandran.range.color, data->ramachandran.range.radius);
        }
    }
    ImGui::EndColumns();

    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    ImGui::BeginChild("canvas", ImVec2(-1, -1), true);
    ImGui::PopStyleVar(1);

    static float zoom_factor = 1.0f;
    const float max_c = ImMax(ImGui::GetContentRegionAvail().x, ImGui::GetContentRegionAvail().y);
    const ImVec2 size = ImVec2(max_c, max_c) * zoom_factor;
    const ImGuiID id = ImGui::GetID("bg");

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
    dl->ChannelsSetCurrent(2);
    dl->AddImage((ImTextureID)(intptr_t)ramachandran::get_accumulation_texture(), bb.Min, bb.Max);
    dl->ChannelsSetCurrent(3);

    constexpr float32 ONE_OVER_TWO_PI = 1.f / (2.f * math::PI);

    int64 mouse_hover_idx = -1;

    if (data->ramachandran.current.enabled) {
        const uint32 border_color = math::convert_color(data->ramachandran.current.border_color);
        const uint32 base_color = math::convert_color(data->ramachandran.current.fill_color);
        const uint32 selected_color = math::convert_color(vec4(data->selection.selection_color, 1));
        const uint32 highlight_color = math::convert_color(vec4(data->selection.highlight_color, 1));
        const float32 base_radius = data->ramachandran.current.radius;
        const float32 selected_radius = data->ramachandran.selected.radius;

        for (int64 i = 0; i < backbone_segments.size(); i++) {
            const auto& angle = current_angles[i];
            const auto& seg = backbone_segments[i];
            const auto& res = residues[i];
            if (angle.x == 0.f || angle.y == 0.f) continue;

            const auto selected = !is_array_zero(atom_selection.subarray(res.atom_idx));
            const auto highlight = !is_array_zero(atom_highlight.subarray(res.atom_idx));
            const auto radius = (highlight || selected) ? selected_radius : base_radius;
            const auto fill_color = highlight ? highlight_color : (selected ? selected_color : base_color);

            const ImVec2 coord = ImLerp(bb.Min, bb.Max, ImVec2(angle.x * ONE_OVER_TWO_PI + 0.5f, -angle.y * ONE_OVER_TWO_PI + 0.5f));  // [-PI, PI] -> [0, 1]
            const ImVec2 min_box(math::round(coord.x - radius), math::round(coord.y - radius));
            const ImVec2 max_box(math::round(coord.x + radius), math::round(coord.y + radius));
            if (radius > 1.f) {
                dl->AddRectFilled(min_box, max_box, fill_color);
                dl->AddRect(min_box, max_box, border_color);
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

    enum Mode { Append, Remove };
    static ImVec2 region_x0 = {0, 0};
    static bool region_select = false;
    static Mode region_mode = Append;

    if (ImGui::IsItemHovered()) {

        if (ImGui::GetIO().KeyCtrl && ImGui::GetIO().MouseWheel != 0.f) {
            const float old_zoom = zoom_factor;
            const float new_zoom = ImClamp(zoom_factor - ImGui::GetIO().MouseWheel * 0.1f, 1.0f, 10.0f);
            // ImGui::GetCurrentWindow()->ScrollTarget.x = ImGui::GetCurrentWindow()->Scroll.x * new_zoom / old_zoom;
            zoom_factor = new_zoom;
        }

        if (!region_select && ImGui::GetIO().KeyShift && (ImGui::GetIO().MouseClicked[0] || ImGui::GetIO().MouseClicked[1])) {
            region_select = true;
            region_mode = ImGui::GetIO().MouseClicked[0] ? Append : Remove;
            region_x0 = ImGui::GetMousePos();
        }

        const ImVec2 normalized_coord = ((ImGui::GetMousePos() - bb.Min) / (bb.Max - bb.Min) - ImVec2(0.5f, 0.5f)) * ImVec2(1, -1);
        const ImVec2 angles = normalized_coord * 2.f * 180.f;
        ImGui::BeginTooltip();
        ImGui::Text(u8"\u03C6: %.1f\u00b0, \u03C8: %.1f\u00b0", angles.x, angles.y);
        if (!region_select && mouse_hover_idx != -1) {
            const auto res_idx = mouse_hover_idx;
            const auto& res = get_residues(data->mol_data.dynamic.molecule)[res_idx];
            ImGui::Text("Residue[%lli]: %s", res_idx, res.name.cstr());
            memset_array(data->selection.current_highlight_mask, false);
            memset_array(data->selection.current_highlight_mask.subarray(res.atom_idx), true);
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
        ImDrawList* dl = ImGui::GetWindowDrawList();
        dl->AddRectFilled(x0, x1, fill_col);
        dl->AddRect(x0, x1, line_col);

        static DynamicArray<bool> mask;
        mask.resize(mol.atom.count);
        memset_array(mask, false);

        for (int64 i = 0; i < residues.size(); i++) {
            const auto& angle = current_angles[i];
            if (angle.x == 0.f || angle.y == 0.f) continue;
            const ImVec2 coord = ImLerp(bb.Min, bb.Max, ImVec2(angle.x * ONE_OVER_TWO_PI + 0.5f, -angle.y * ONE_OVER_TWO_PI + 0.5f));  // [-PI, PI] -> [0, 1]
            if (coord.x < x0.x || x1.x < coord.x) continue;
            if (coord.y < x0.y || x1.y < coord.y) continue;

            const auto& seg = backbone_segments[i];
            mask[seg.ca_idx] = true;
        }

        switch (data->selection.level_mode) {
            case SelectionLevel::Atom:
                break;
            case SelectionLevel::Residue:
                expand_mask_to_residue(mask, mol.residues);
                break;
            case SelectionLevel::Chain:
                expand_mask_to_chain(mask, mol.chains);
                break;
            default:
                ASSERT(false);
        }

        memset_array(data->selection.current_highlight_mask, false);

        if (region_mode == Append) {
            for (int64 i = 0; i < mask.size(); i++) {
                data->selection.current_highlight_mask[i] = data->selection.current_selection_mask[i] || mask[i];
            }
        } else if (region_mode == Remove) {
            for (int64 i = 0; i < mask.size(); i++) {
                data->selection.current_highlight_mask[i] = data->selection.current_selection_mask[i] && !mask[i];
            }
        }

        if (!ImGui::GetIO().KeyShift || !(ImGui::GetIO().MouseDown[0] || ImGui::GetIO().MouseDown[1])) {
            if (region_mode == Append) {
                for (int64 i = 0; i < mask.size(); i++) {
                    data->selection.current_selection_mask[i] |= mask[i];
                }
            } else if (region_mode == Remove) {
                for (int64 i = 0; i < mask.size(); i++) {
                    data->selection.current_selection_mask[i] &= !mask[i];
                }
            }
            region_select = false;
        }
        data->gpu_buffers.dirty.selection = true;
    }

    ImGui::EndChild();

    ImGui::End();
}

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
    if (!fbo->deferred.picking) glGenTextures(1, &fbo->deferred.picking);
    if (!fbo->pbo_picking.color[0]) glGenBuffers(2, fbo->pbo_picking.color);
    if (!fbo->pbo_picking.depth[0]) glGenBuffers(2, fbo->pbo_picking.depth);

    glBindTexture(GL_TEXTURE_2D, fbo->deferred.depth);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH24_STENCIL8, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, fbo->deferred.color);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
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
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, fbo->deferred.emissive);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R11F_G11F_B10F, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
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

    const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3, GL_COLOR_ATTACHMENT4};

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo->deferred.fbo);
    if (attach_textures_deferred) {
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, fbo->deferred.depth, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_STENCIL_ATTACHMENT, GL_TEXTURE_2D, fbo->deferred.depth, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, fbo->deferred.color, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, fbo->deferred.normal, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, fbo->deferred.velocity, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, GL_TEXTURE_2D, fbo->deferred.emissive, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT4, GL_TEXTURE_2D, fbo->deferred.picking, 0);
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
    if (fbo->deferred.picking) glDeleteTextures(1, &fbo->deferred.picking);

    if (fbo->pbo_picking.color[0]) glDeleteBuffers(2, fbo->pbo_picking.color);
    if (fbo->pbo_picking.depth[0]) glDeleteBuffers(2, fbo->pbo_picking.depth);
}

static void init_molecule_buffers(ApplicationData* data) {
    ASSERT(data);
    const auto& mol = data->mol_data.dynamic.molecule;

    DynamicArray<uint32> backbone_index_data;
    DynamicArray<uint32> control_point_index_data;
    DynamicArray<uint32> spline_index_data;

    {
        int32 control_idx = 0;
        int32 spline_idx = 0;
        for (const auto& seq : mol.backbone.sequences) {
            const auto backbone = get_backbone(mol, seq);
            const int32 count = (int32)backbone.size();

            /*
                // These indices are needed to compute the backbone angles
                                phi = math::dihedral_angle(c_idx[i-1], n_idx[i], ca_idx[i], c_idx[i]);
                                psi = math::dihedral_angle(n_idx[i], ca_idx[i], c_idx[i], n_idx[i+1]);
            */

            for (int32 i = 0; i < count; i++) {
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

                // @NOTE: For every control point we generate N spline control points
                if (!last) {
                    for (int32 j = 0; j < SPLINE_SUBDIVISION_COUNT; j++) {
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

    data->gpu_buffers.backbone.num_backbone_segment_indices = (int32)backbone_index_data.size();
    data->gpu_buffers.backbone.num_control_point_indices = (int32)control_point_index_data.size();
    data->gpu_buffers.backbone.num_spline_indices = (int32)spline_index_data.size();

    LOG_NOTE("num backbone segment indices: %i", (int32)backbone_index_data.size());
    LOG_NOTE("num control point indices: %i", (int32)control_point_index_data.size());
    LOG_NOTE("num spline indices: %i", (int32)spline_index_data.size());

    const int64 num_backbone_segments = backbone_index_data.size() / 6;
    const int64 position_buffer_size = mol.atom.count * 3 * sizeof(float);
    const int64 velocity_buffer_size = mol.atom.count * 3 * sizeof(float);
    const int64 bond_buffer_size = mol.covalent_bonds.size() * sizeof(uint32) * 2;
    const int64 selection_buffer_size = mol.atom.count * sizeof(uint8);
    const int64 control_point_buffer_size = num_backbone_segments * sizeof(draw::ControlPoint);
    const int64 spline_buffer_size = control_point_buffer_size * SPLINE_SUBDIVISION_COUNT;

    if (!data->gpu_buffers.position) glGenBuffers(1, &data->gpu_buffers.position);
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

    glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.velocity);
    glBufferData(GL_ARRAY_BUFFER, velocity_buffer_size, NULL, GL_DYNAMIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.radius);
    glBufferData(GL_ARRAY_BUFFER, data->mol_data.atom_radii.size_in_bytes(), data->mol_data.atom_radii.data(), GL_STATIC_DRAW);

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

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    data->gpu_buffers.dirty.position = true;
    data->gpu_buffers.dirty.velocity = true;
    data->gpu_buffers.dirty.selection = true;
    data->gpu_buffers.dirty.backbone = true;
}

static void free_molecule_buffers(ApplicationData* data) {
    ASSERT(data);
    if (data->gpu_buffers.position) {
        glDeleteBuffers(1, &data->gpu_buffers.position);
        data->gpu_buffers.position = 0;
    }
    if (data->gpu_buffers.velocity) {
        glDeleteBuffers(1, &data->gpu_buffers.velocity);
        data->gpu_buffers.velocity = 0;
    }
    if (data->gpu_buffers.radius) {
        glDeleteBuffers(1, &data->gpu_buffers.radius);
        data->gpu_buffers.radius = 0;
    }
    if (data->gpu_buffers.selection) {
        glDeleteBuffers(1, &data->gpu_buffers.selection);
        data->gpu_buffers.selection = 0;
    }
    if (data->gpu_buffers.backbone.backbone_segment_index) {
        glDeleteBuffers(1, &data->gpu_buffers.backbone.backbone_segment_index);
        data->gpu_buffers.backbone.backbone_segment_index = 0;
    }
    if (data->gpu_buffers.backbone.control_point) {
        glDeleteBuffers(1, &data->gpu_buffers.backbone.control_point);
        data->gpu_buffers.backbone.control_point = 0;
    }
    if (data->gpu_buffers.backbone.control_point_index) {
        glDeleteBuffers(1, &data->gpu_buffers.backbone.control_point_index);
        data->gpu_buffers.backbone.control_point_index = 0;
    }
    if (data->gpu_buffers.backbone.spline) {
        glDeleteBuffers(1, &data->gpu_buffers.backbone.spline);
        data->gpu_buffers.backbone.spline = 0;
    }
    if (data->gpu_buffers.backbone.spline_index) {
        glDeleteBuffers(1, &data->gpu_buffers.backbone.spline_index);
        data->gpu_buffers.backbone.spline_index = 0;
    }
    if (data->gpu_buffers.bond) {
        glDeleteBuffers(1, &data->gpu_buffers.bond);
        data->gpu_buffers.bond = 0;
    }
}

void copy_molecule_data_to_buffers(ApplicationData* data) {
    ASSERT(data);
    const auto N = data->mol_data.dynamic.molecule.atom.count;

    if (data->gpu_buffers.dirty.position) {
        data->gpu_buffers.dirty.position = false;
        const vec3* position = data->mol_data.dynamic.molecule.atom.positions;

        // Update data inside position buffer
        glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.position);

        float* pos_gpu = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        if (!pos_gpu) {
            LOG_ERROR("Could not map position buffer");
            return;
        }

        // @NOTE: we cannot do a pure memcpy since the CPU buffer using vec3 can be a m128 with 16-byte alignment.
        for (int64 i = 0; i < N; i++) {
            pos_gpu[i * 3 + 0] = position[i][0];
            pos_gpu[i * 3 + 1] = position[i][1];
            pos_gpu[i * 3 + 2] = position[i][2];
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }

    if (data->gpu_buffers.dirty.velocity) {
        data->gpu_buffers.dirty.velocity = false;
        const vec3* velocity = data->mol_data.atom_velocity.data();

        // Update data inside position buffer
        glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.velocity);

        float* vel_gpu = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        if (!vel_gpu) {
            LOG_ERROR("Could not map velocity buffer");
            return;
        }

        // @NOTE: we cannot do a pure memcpy since the CPU buffer using vec3 can be a m128 with 16-byte alignment.
        for (int64 i = 0; i < N; i++) {
            vel_gpu[i * 3 + 0] = velocity[i][0];
            vel_gpu[i * 3 + 1] = velocity[i][1];
            vel_gpu[i * 3 + 2] = velocity[i][2];
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }

    if (data->gpu_buffers.dirty.selection) {
        data->gpu_buffers.dirty.selection = false;

        glBindBuffer(GL_ARRAY_BUFFER, data->gpu_buffers.selection);
        uint8* sel_gpu = (uint8*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

        if (!sel_gpu) {
            LOG_ERROR("Could not map selection buffer");
            return;
        }

        for (int64 i = 0; i < N; i++) {
            sel_gpu[i] = (data->selection.current_selection_mask[i] ? 2 : 0) | (data->selection.current_highlight_mask[i] ? 1 : 0);
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// #moleculedata
static void free_trajectory_data(ApplicationData* data) {
    ASSERT(data);
    if (data->mol_data.dynamic.trajectory) {
        data->async.trajectory.sync.signal_stop_and_wait();
        stats::signal_stop_and_wait();
        data->async.backbone_angles.sync.signal_stop_and_wait();
        close_file_handle(&data->mol_data.dynamic.trajectory);
        free_trajectory(&data->mol_data.dynamic.trajectory);
    }
}

static void free_molecule_data(ApplicationData* data) {
    ASSERT(data);
    if (data->mol_data.dynamic.molecule) {
        data->files.molecule = "";
        free_molecule_structure(&data->mol_data.dynamic.molecule);
    }
    if (data->mol_data.dynamic.trajectory) {
        data->files.trajectory = "";
        free_trajectory_data(data);
    }
    data->mol_data.atom_radii.clear();
    data->mol_data.atom_velocity.clear();

    free_molecule_buffers(data);
    free_backbone_angles_trajectory(&data->ramachandran.backbone_angles);
    data->ramachandran.backbone_angles = {};
    data->hydrogen_bonds.bonds.clear();
    data->hydrogen_bonds.dirty = true;
    data->gpu_buffers.dirty.backbone = true;
    data->selection.current_selection_mask.clear();
    data->selection.current_highlight_mask.clear();
}

static void init_molecule_data(ApplicationData* data) {
    if (data->mol_data.dynamic.molecule) {
        const auto atom_count = data->mol_data.dynamic.molecule.atom.count;
        data->mol_data.atom_radii = compute_atom_radii(get_elements(data->mol_data.dynamic.molecule));
        data->mol_data.atom_velocity.resize(atom_count);
        zero_array(data->mol_data.atom_velocity);
        data->selection.current_selection_mask.resize(atom_count);
        data->selection.current_highlight_mask.resize(atom_count);
        zero_array(data->selection.current_selection_mask);
        zero_array(data->selection.current_highlight_mask);
        init_molecule_buffers(data);
        data->picking.idx = NO_PICKING_IDX;
        data->selection.hovered = -1;
        data->selection.right_clicked = -1;
    }
}

static void init_trajectory_data(ApplicationData* data) {
    if (data->mol_data.dynamic.trajectory) {
        if (data->mol_data.dynamic.trajectory.num_atoms != data->mol_data.dynamic.molecule.atom.count) {
            LOG_ERROR("ERROR! The number of atoms in the molecule does not match the number of atoms in the trajectory.");
            free_trajectory_data(data);
            return;
        }

        read_next_trajectory_frame(&data->mol_data.dynamic.trajectory);  // read first frame explicitly
        auto frame_0_pos = get_trajectory_positions(data->mol_data.dynamic.trajectory, 0);
        memcpy(data->mol_data.dynamic.molecule.atom.positions, frame_0_pos.data(), frame_0_pos.size_in_bytes());
        data->gpu_buffers.dirty.position = true;

        load_trajectory_async(data);

        create_volume(data);
#if 1
        if (data->mol_data.dynamic.trajectory.num_frames > 0) {
            vec3 box_ext = data->mol_data.dynamic.trajectory.frame_buffer[0].box * vec3(1);
            init_volume(&data->density_volume.volume, math::max(ivec3(1), ivec3(box_ext) / VOLUME_DOWNSAMPLE_FACTOR));
            data->density_volume.model_to_world_matrix = volume::compute_model_to_world_matrix(vec3(0), box_ext);
            data->density_volume.texture_to_model_matrix = volume::compute_texture_to_model_matrix(data->density_volume.volume.dim);
        }
#endif

        init_backbone_angles_trajectory(&data->ramachandran.backbone_angles, data->mol_data.dynamic);
        compute_backbone_angles_trajectory(&data->ramachandran.backbone_angles, data->mol_data.dynamic);
    }
}

static void load_molecule_data(ApplicationData* data, CString file) {
    ASSERT(data);
    if (file.count > 0) {
        data->is_playing = false;
        CString ext = get_file_extension(file);
        LOG_NOTE("Loading molecular data from file '%.*s'...", file.count, file.ptr);
        auto t0 = platform::get_time();
        if (compare_ignore_case(ext, "pdb")) {
            free_molecule_data(data);
            if (!allocate_and_load_pdb_from_file(&data->mol_data.dynamic, file)) {
                LOG_ERROR("ERROR! Failed to load pdb file.");
            }
            data->files.molecule = file;
            init_molecule_data(data);
            init_trajectory_data(data);
        } else if (compare_ignore_case(ext, "gro")) {
            free_molecule_data(data);
            if (!allocate_and_load_gro_from_file(&data->mol_data.dynamic.molecule, file)) {
                LOG_ERROR("ERROR! Failed to load gro file.");
                return;
            }
            data->files.molecule = file;
            init_molecule_data(data);
        } else if (compare_ignore_case(ext, "xtc")) {
            if (!data->mol_data.dynamic.molecule) {
                LOG_ERROR("ERROR! Must have molecule structure before trajectory can be loaded.");
                return;
            }
            free_trajectory_data(data);
            if (!load_and_allocate_trajectory(&data->mol_data.dynamic.trajectory, file)) {
                LOG_ERROR("ERROR! Problem loading trajectory.");
                return;
            }
            init_trajectory_data(data);
        } else {
            LOG_ERROR("ERROR! file extension is not supported!\n");
            return;
        }
        auto t1 = platform::get_time();
        LOG_NOTE("Success! operation took %.3fs.", platform::compute_delta_ms(t0, t1) / 1000.f);
        LOG_NOTE("Number of chains: %i", (int32)data->mol_data.dynamic.molecule.chains.size());
        LOG_NOTE("Number of residues: %i", (int32)data->mol_data.dynamic.molecule.residues.size());
        LOG_NOTE("Number of atoms: %i", (int32)data->mol_data.dynamic.molecule.atom.count);
    }
}

// ### WORKSPACE ###
static RepresentationType get_rep_type(CString str) {
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

static CString get_rep_type_name(RepresentationType type) {
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

static ColorMapping get_color_mapping(CString str) {
    if (compare(str, "STATIC_COLOR"))
        return ColorMapping::Static;
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

static CString get_color_mapping_name(ColorMapping mapping) {
    switch (mapping) {
        case ColorMapping::Static:
            return "STATIC_COLOR";
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

static vec4 to_vec4(CString txt, const vec4& default_val = vec4(1)) {
    vec4 res = default_val;
    auto tokens = ctokenize(txt, ",");
    int32 count = (int32)tokens.count < 4 ? (int32)tokens.count : 4;
    for (int i = 0; i < count; i++) {
        res[i] = to_float(tokens[i]);
    }
    return res;
}

static void load_workspace(ApplicationData* data, CString file) {
    ASSERT(data);
    clear_representations(data);
    stats::remove_all_properties();

    StringBuffer<256> new_molecule_file;
    StringBuffer<256> new_trajectory_file;

    String txt = allocate_and_read_textfile(file);
    defer { FREE(txt); };

    CString c_txt = txt;
    CString line;
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
                    if (compare_n(line, "StaticColor=", 12)) rep->static_color = to_vec4(trim(line.substr(12)));
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
            }
        } else if (compare_n(line, "[Camera]", 8)) {
            while (c_txt && c_txt[0] != '[' && (line = extract_line(c_txt))) {
                if (compare_n(line, "Position=", 9)) {
                    vec3 pos = vec3(to_vec4(trim(line.substr(9))));
                    data->view.camera.position = pos;
                    data->view.animation.target_position = pos;
                }
                if (compare_n(line, "Rotation=", 9)) {
                    quat rot = quat(to_vec4(trim(line.substr(9))));
                    data->view.camera.orientation = rot;
                }
                if (compare_n(line, "Distance=", 9)) data->view.trackball_state.distance = to_float(trim(line.substr(9)));
            }
        }
    }

    data->files.workspace = file;

    if (!compare(new_molecule_file, data->files.molecule) && new_molecule_file) {
        load_molecule_data(data, new_molecule_file);
    }

    if (!compare(new_trajectory_file, data->files.trajectory) && new_trajectory_file) {
        load_molecule_data(data, new_trajectory_file);
    }

    reset_view(data, false, true);
    reset_representations(data);
}

static void save_workspace(ApplicationData* data, CString file) {
    FILE* fptr = fopen(file, "w");
    if (!fptr) {
        printf("ERROR! Could not save workspace to file '%s'\n", file.beg());
        return;
    }

    // @TODO: Make relative paths
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
        fprintf(fptr, "StaticColor=%g,%g,%g,%g\n", rep.static_color.r, rep.static_color.g, rep.static_color.b, rep.static_color.a);
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

    fprintf(fptr, "[Camera]\n");
    fprintf(fptr, "Position=%g,%g,%g\n", data->view.camera.position.x, data->view.camera.position.y, data->view.camera.position.z);
    fprintf(fptr, "Rotation=%g,%g,%g,%g\n", data->view.camera.orientation.x, data->view.camera.orientation.y, data->view.camera.orientation.z, data->view.camera.orientation.w);
    fprintf(fptr, "Distance=%g\n", data->view.trackball_state.distance);
    fprintf(fptr, "\n");

    fclose(fptr);

    data->files.workspace = file;
}

void create_screenshot(ApplicationData* data) {
    ASSERT(data);
    Image img;
    init_image(&img, data->fbo.width, data->fbo.height);

    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    glReadBuffer(GL_BACK);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    glReadPixels(0, 0, img.width, img.height, GL_RGBA, GL_UNSIGNED_BYTE, img.data);

    time_t now = time(0);
    struct tm tstruct;
    tstruct = *localtime(&now);
    char time_str[80];
    strftime(time_str, sizeof(time_str), "%Y-%m-%d_%X", &tstruct);

    StringBuffer<256> path;
    path += VIAMD_SCREENSHOT_DIR;
    path += "/screenshot";
    path += "";
    path += ".bmp";

    bool res = write_image(img, path);
}

// #representation
static Representation* create_representation(ApplicationData* data, RepresentationType type, ColorMapping color_mapping, CString filter) {
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
    if (rep.atom_mask) free(&rep.atom_mask);
    data->representations.buffer.remove(&rep);
}

static void recompute_atom_visibility_mask(ApplicationData* data) {
    ASSERT(data);
    const auto N = data->mol_data.dynamic.molecule.atom.count;
    data->representations.atom_visibility_mask.resize(N);
    memset_array(data->representations.atom_visibility_mask, false);

    for (const auto& rep : data->representations.buffer) {
        for (int64 i = 0; i < N; i++) {
            data->representations.atom_visibility_mask[i] |= rep.atom_mask[i];
        }
    }
}

static void update_representation(ApplicationData* data, Representation* rep) {
    ASSERT(data);
    ASSERT(rep);
    uint32 static_color = ImGui::ColorConvertFloat4ToU32(vec_cast(rep->static_color));
    DynamicArray<uint32> colors(data->mol_data.dynamic.molecule.atom.count);
    const auto& mol = data->mol_data.dynamic.molecule;

    switch (rep->color_mapping) {
        case ColorMapping::Static:
            memset_array(colors, static_color);
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
        free_array(&rep->atom_mask);
        rep->atom_mask = allocate_array<bool>(mol.atom.count);
    }
    rep->filter_is_ok = filter::compute_filter_mask(rep->atom_mask, data->mol_data.dynamic, rep->filter.buffer);
    filter::filter_colors(colors, rep->atom_mask);
    recompute_atom_visibility_mask(data);

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
        remove_representation(data, (int32)data->representations.buffer.size() - 1);
    }
}

// #selection
static Selection* create_selection(ApplicationData* data, CString name, Array<const bool> atom_mask) {
    ASSERT(data);
    Selection sel;
    sel.name = name;
    sel.atom_mask = allocate_array<bool>(data->selection.current_selection_mask.size());
    memcpy(sel.atom_mask.data(), atom_mask.data(), math::min(sel.atom_mask.size_in_bytes(), atom_mask.size_in_bytes()));
    return &data->selection.stored_selections.push_back(sel);
}

static Selection* clone_selection(ApplicationData* data, const Selection& src) {
    ASSERT(data);
    Selection clone;
    clone.name = src.name;
    clone.atom_mask = allocate_array<bool>(data->selection.current_selection_mask.size());
    memcpy(clone.atom_mask.data(), src.atom_mask.data(), math::min(clone.atom_mask.size_in_bytes(), src.atom_mask.size_in_bytes()));
    return &data->selection.stored_selections.push_back(clone);
}

static void remove_selection(ApplicationData* data, int idx) {
    ASSERT(data);
    if (idx < 0 || (int)data->selection.stored_selections.size() <= idx) {
        LOG_ERROR("Index [%i] out of range when trying to remove selection", idx);
    }
    auto item = &data->selection.stored_selections[idx];
    free_array(&item->atom_mask);
    data->selection.stored_selections.remove(item);
}

static void reset_selections(ApplicationData* data) {
    ASSERT(data);
    // @NOTE: What to do here?
}

static void clear_selections(ApplicationData* data) {
    ASSERT(data);
    while (data->selection.stored_selections.size() > 0) {
        remove_selection(data, (int32)data->selection.stored_selections.size() - 1);
    }
}

static bool handle_selection(ApplicationData* data) {
    ASSERT(data);
    enum class RegionMode { Append, Remove };

    static RegionMode region_mode = RegionMode::Append;
    static DynamicArray<bool> mask = {};
    static bool region_select = false;
    static platform::Coordinate x0;
    const platform::Coordinate x1 = data->ctx.input.mouse.win_coord;
    const int64 N = data->mol_data.dynamic.molecule.atom.count;
    const bool shift_down = data->ctx.input.key.down[Key::KEY_LEFT_SHIFT] || data->ctx.input.key.down[Key::KEY_RIGHT_SHIFT];
    const bool mouse_down = data->ctx.input.mouse.down[0] || data->ctx.input.mouse.down[1];

    mask.resize(data->mol_data.dynamic.molecule.atom.count);
    memset_array(mask, false);

    memset_array(data->selection.current_highlight_mask, false);

    Range<int32> picking_range = {0, 0};

    if (data->picking.idx != NO_PICKING_IDX && !region_select) {
        switch (data->selection.level_mode) {
            case SelectionLevel::Atom:
                picking_range = {(int32)data->picking.idx, (int32)data->picking.idx + 1};
                break;
            case SelectionLevel::Residue: {
                const auto res_idx = data->mol_data.dynamic.molecule.atom.residue_indices[data->picking.idx];
                if (0 <= res_idx && res_idx < data->mol_data.dynamic.molecule.residues.size()) {
                    const auto& res = data->mol_data.dynamic.molecule.residues[res_idx];
                    picking_range = {res.atom_idx.beg, res.atom_idx.end};
                }
                break;
            }
            case SelectionLevel::Chain: {
                const auto res_idx = data->mol_data.dynamic.molecule.atom.residue_indices[data->picking.idx];
                if (0 <= res_idx && res_idx < data->mol_data.dynamic.molecule.residues.size()) {
                    const auto chain_idx = data->mol_data.dynamic.molecule.residues[res_idx].chain_idx;
                    if (0 <= chain_idx && chain_idx < data->mol_data.dynamic.molecule.chains.size()) {
                        const auto& chain = data->mol_data.dynamic.molecule.chains[chain_idx];
                        picking_range = {chain.atom_idx.beg, chain.atom_idx.end};
                    }
                }
                break;
            }
            default:
                assert(false);
                break;
        }

        for (int64 i = picking_range.x; i < picking_range.y; i++) {
            mask[i] = true;
        }
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
            const mat4 mvp = compute_perspective_projection_matrix(data->view.camera, data->ctx.window.width, data->ctx.window.height) * data->view.param.matrix.view;
            const Array<const vec3> positions = get_positions(data->mol_data.dynamic.molecule);

            for (int64 i = 0; i < N; i++) {
                if (!data->representations.atom_visibility_mask[i]) continue;
                vec4 p = mvp * vec4(positions[i], 1);
                p /= p.w;
                const vec2 c = (vec2(p.x, -p.y) * 0.5f + 0.5f) * res;
                mask[i] = (min_p.x <= c.x && c.x <= max_p.x && min_p.y <= c.y && c.y <= max_p.y);
            }

            switch (data->selection.level_mode) {
                case SelectionLevel::Atom:
                    break;
                case SelectionLevel::Residue:
                    expand_mask_to_residue(mask, data->mol_data.dynamic.molecule.residues);
                    break;
                case SelectionLevel::Chain:
                    expand_mask_to_chain(mask, data->mol_data.dynamic.molecule.chains);
                    break;
                default:
                    ASSERT(false);
            }

            /*
switch (region_mode) {
    case RegionMode::Append:
        for (int64 i = 0; i < N; i++) {
            mask[i] = data->selection.current_selection_mask[i] && mask[i];
        }
        break;
    case RegionMode::Remove:
        for (int64 i = 0; i < N; i++) {
            mask[i] = data->selection.current_selection_mask[i] && !mask[i];
        }
        break;
    default:
        ASSERT(false);
}
            */

            if (!mouse_down) {
                for (int64 i = 0; i < data->selection.current_selection_mask.size(); i++) {
                    if (region_mode == RegionMode::Append) {
                        data->selection.current_selection_mask[i] |= mask[i];
                    } else if (region_mode == RegionMode::Remove) {
                        data->selection.current_selection_mask[i] &= !mask[i];
                    }
                }
                region_select = false;
            }
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
                const bool val = data->ctx.input.mouse.clicked[0];
                for (int64 i = picking_range.x; i < picking_range.y; i++) {
                    data->selection.current_selection_mask[i] = val;
                }
            } else if (data->ctx.input.mouse.clicked[1]) {
                // Clear selection
                memset_array(data->selection.current_selection_mask, false);
            }
        }
    } else {
        region_select = false;
    }

    for (int64 i = 0; i < N; i++) {
        data->selection.current_highlight_mask[i] = mask[i];
    }
    data->gpu_buffers.dirty.selection = true;

    return region_select;
}

// #async
static void load_trajectory_async(ApplicationData* data) {
    ASSERT(data);
    // Wait for thread to finish if already running
    if (data->async.trajectory.sync.running) {
        data->async.trajectory.sync.signal_stop_and_wait();
    }

    if (data->mol_data.dynamic.trajectory.file_handle) {
        data->async.trajectory.sync.stop_signal = false;
        data->async.trajectory.sync.running = true;
        data->async.trajectory.sync.thread = std::thread([data]() {
            while (read_next_trajectory_frame(&data->mol_data.dynamic.trajectory)) {
                data->async.trajectory.fraction = data->mol_data.dynamic.trajectory.num_frames / (float)data->mol_data.dynamic.trajectory.frame_offsets.count;
                if (data->async.trajectory.sync.stop_signal) break;
            }
            data->async.trajectory.sync.running = false;
            data->async.trajectory.sync.stop_signal = false;

            // compute_statistics_async(data);
            stats::set_all_property_flags(true, true);
            compute_backbone_angles_async(data);
        });
        data->async.trajectory.sync.thread.detach();
    }
}

static void compute_backbone_angles_async(ApplicationData* data) {
    ASSERT(data);
    data->async.backbone_angles.query_update = true;
    if (data->async.backbone_angles.sync.running == false) {
        data->async.backbone_angles.sync.running = true;

        data->async.backbone_angles.sync.thread = std::thread([data]() {
            data->async.backbone_angles.fraction = 0.0f;
            while (data->async.backbone_angles.query_update) {
                data->async.backbone_angles.query_update = false;
                data->async.backbone_angles.fraction = 0.5f;
                compute_backbone_angles_trajectory(&data->ramachandran.backbone_angles, data->mol_data.dynamic);
                if (data->async.backbone_angles.sync.stop_signal) break;
            }
            data->async.backbone_angles.fraction = 1.f;
            data->async.backbone_angles.sync.running = false;
            data->async.backbone_angles.sync.stop_signal = false;
        });
        data->async.backbone_angles.sync.thread.detach();
    }
}

static void create_volume(ApplicationData* data) {
    const vec3 min_box = vec3(0);
    const vec3 max_box = data->mol_data.dynamic.trajectory.num_frames > 0 ? data->mol_data.dynamic.trajectory.frame_buffer[0].box * vec3(1) : vec3(1);
    const ivec3 dim = math::max(ivec3(1), ivec3(max_box) / VOLUME_DOWNSAMPLE_FACTOR);
    init_volume(&data->density_volume.volume, dim);
    // volume::create_volume_texture(&data->density_volume.texture, dim);
    data->density_volume.model_to_world_matrix = volume::compute_model_to_world_matrix(min_box, max_box);
    data->density_volume.texture_to_model_matrix = volume::compute_texture_to_model_matrix(dim);
    data->density_volume.world_to_texture_matrix = math::inverse(data->density_volume.model_to_world_matrix * data->density_volume.texture_to_model_matrix);
}

static void draw_representations_lean_and_mean(ApplicationData* data, vec4 color, float scale, uint32 mask) {
	ASSERT(data);
	const int32 atom_count = (int32)data->mol_data.dynamic.molecule.atom.count;
	const int32 bond_count = (int32)data->mol_data.dynamic.molecule.covalent_bonds.size();

	PUSH_GPU_SECTION("Lean and Mean")
	for (const auto& rep : data->representations.buffer) {
		if (!rep.enabled) continue;
		switch (rep.type) {
		case RepresentationType::Vdw:
			PUSH_GPU_SECTION("Vdw")
			draw::lean_and_mean::draw_vdw(data->gpu_buffers.position, data->gpu_buffers.radius, data->gpu_buffers.selection, atom_count, data->view.param, rep.radius * scale, color, mask);
			POP_GPU_SECTION()
			break;
		case RepresentationType::Licorice:
			PUSH_GPU_SECTION("Licorice")
			draw::lean_and_mean::draw_licorice(data->gpu_buffers.position, data->gpu_buffers.selection, data->gpu_buffers.bond, bond_count, data->view.param, rep.radius * scale, color, mask);
			POP_GPU_SECTION()
			break;
		case RepresentationType::BallAndStick:
			PUSH_GPU_SECTION("Vdw")
			draw::lean_and_mean::draw_vdw(data->gpu_buffers.position, data->gpu_buffers.radius, data->gpu_buffers.selection, atom_count, data->view.param, rep.radius * scale * BALL_AND_STICK_VDW_SCALE, color, mask);
			POP_GPU_SECTION()
			PUSH_GPU_SECTION("Licorice")
			draw::lean_and_mean::draw_licorice(data->gpu_buffers.position, data->gpu_buffers.selection, data->gpu_buffers.bond, bond_count, data->view.param, rep.radius * scale * BALL_AND_STICK_LICORICE_SCALE, color, mask);
			POP_GPU_SECTION()
			break;
		case RepresentationType::Ribbons:
			PUSH_GPU_SECTION("Ribbons")
			draw::lean_and_mean::draw_ribbons(data->gpu_buffers.backbone.spline, data->gpu_buffers.backbone.spline_index, data->gpu_buffers.selection, data->gpu_buffers.backbone.num_spline_indices, data->view.param, scale, color, mask);
			POP_GPU_SECTION()
			break;
		default:
			break;
		}
	}
	POP_GPU_SECTION()

}