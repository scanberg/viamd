#include <imgui.h>
#include <imgui_internal.h>
#include <core/platform.h>
#include <core/gl.h>
#include <core/types.h>
#include <core/hash.h>
#include <core/log.h>
#include <core/math_utils.h>
#include <core/camera.h>
#include <core/camera_utils.h>
#include <core/console.h>
#include <core/string_utils.h>
#include <core/volume.h>
#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>
#include <mol/trajectory_utils.h>
#include <mol/molecule_utils.h>
#include <mol/pdb_utils.h>
#include <mol/gro_utils.h>
#include <mol/stats.h>
#include <mol/spatial_hash.h>
#include <gfx/immediate_draw_utils.h>
#include <gfx/postprocessing_utils.h>
#include <gfx/volume_utils.h>
#include <gfx/raytracing_utils.h>

#include <mol/radial_basis.h>

#include <stdio.h>
#include <thread>
#include <mutex>

#include <iostream>
#include <glm/gtx/io.hpp>

//#define VIAMD_RELEASE

#ifdef OS_MAC_OSX
constexpr Key::Key_t CONSOLE_KEY = Key::KEY_WORLD_1;
#else  // WIN32 and Linux
// @TODO: Make sure this is currect for Linux?
constexpr Key::Key_t CONSOLE_KEY = Key::KEY_GRAVE_ACCENT;
#endif

constexpr unsigned int NO_PICKING_IDX = 0xffffffff;
constexpr const char* FILE_EXTENSION = "via";

constexpr uint32 DEL_BTN_COLOR = 0xff1111cc;
constexpr uint32 DEL_BTN_HOVER_COLOR = 0xff3333dd;
constexpr uint32 DEL_BTN_ACTIVE_COLOR = 0xff5555ff;
constexpr uint32 TEXT_BG_ERROR_COLOR = 0xaa222299;

constexpr float HYDROGEN_BOND_DISTANCE_CUTOFF_DEFAULT = 3.0f;
constexpr float HYDROGEN_BOND_DISTANCE_CUTOFF_MIN = 0.1f;
constexpr float HYDROGEN_BOND_DISTANCE_CUTOFF_MAX = 12.0f;

constexpr float HYDROGEN_BOND_ANGLE_CUTOFF_DEFAULT = 20.f;
constexpr float HYDROGEN_BOND_ANGLE_CUTOFF_MIN = 5.f;
constexpr float HYDROGEN_BOND_ANGLE_CUTOFF_MAX = 90.f;

constexpr int32 VOLUME_DOWNSAMPLE_FACTOR = 2;

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

inline ImVec4& vec_cast(vec4& v) { return *(ImVec4*)(&v); }
inline vec4& vec_cast(ImVec4& v) { return *(vec4*)(&v); }
inline ImVec2& vec_cast(vec2& v) { return *(ImVec2*)(&v); }
inline vec2& vec_cast(ImVec2& v) { return *(vec2*)(&v); }

inline ImVec2 operator+(const ImVec2& a, const ImVec2& b) { return {a.x + b.x, a.y + b.y}; }
inline ImVec2 operator-(const ImVec2& a, const ImVec2& b) { return {a.x - b.x, a.y - b.y}; }
inline ImVec2 operator*(const ImVec2& a, const ImVec2& b) { return {a.x * b.x, a.y * b.y}; }
inline ImVec2 operator/(const ImVec2& a, const ImVec2& b) { return {a.x / b.x, a.y / b.y}; }

enum PlaybackInterpolationMode { NEAREST, LINEAR, LINEAR_PERIODIC, CUBIC, CUBIC_PERIODIC };

struct CameraTransformation {
    mat4 world_to_view;
    mat4 view_to_world;
    // Projection
    mat4 view_to_clip;
    mat4 clip_to_view;
};

struct MainFramebuffer {
    GLuint id = 0;
    GLuint tex_depth = 0;
    GLuint tex_color = 0;
    GLuint tex_normal = 0;
    GLuint tex_picking = 0;
    int width = 0;
    int height = 0;
};

struct Representation {
    enum Type { VDW, LICORICE, RIBBONS };

    StringBuffer<128> name = "rep";
    StringBuffer<128> filter = "all";
    Type type = VDW;
    ColorMapping color_mapping = ColorMapping::CPK;
    Array<uint32> colors{};

    bool enabled = true;
    bool filter_is_ok = true;

    // Static color mode
    vec4 static_color = vec4(1);

    // VDW and Licorice
    float radius = 1.f;

    // Ribbons and other spline primitives
    int num_subdivisions = 8;
    float tension = 0.5f;
    float width = 1.f;
    float thickness = 1.f;
};

struct AtomSelection {
    int32 atom_idx = -1;
    int32 residue_idx = -1;
    int32 chain_idx = -1;
};

struct ThreadSyncData {
    std::thread thread{};
    volatile bool running = false;
    volatile bool stop_signal = false;

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

    // --- FILES ---
    // for keeping track of open files
    struct {
        String molecule{};
        String trajectory{};
        String workspace{};
    } files;

    // --- CAMERA ---
    struct {
        Camera camera;
        TrackballControllerState trackball_state{};
        CameraTransformation matrices{};
    } camera;

    // --- MOLECULAR DATA ---
    struct {
        MoleculeDynamic dynamic{};
        DynamicArray<float> atom_radii{};
    } mol_data;

    // --- THREAD SYNCHRONIZATION ---
    struct {
        struct {
            ThreadSyncData sync{};
            float fraction = 0.f;
        } trajectory;

        struct {
            ThreadSyncData sync{};
            float fraction = 0.f;
            bool query_update = false;
        } backbone_angles;
    } async;

    struct {
        bool show_window = false;
        DynamicArray<Representation> data{};
    } representations;

    // --- ATOM SELECTION ---
    AtomSelection hovered{};
    AtomSelection selected{};
    AtomSelection right_clicked{};
    ImGuiID text_field_target = 0;

    // --- STATISTICS ---
    struct {
        bool show_property_window = false;
        bool show_timeline_window = false;
        bool show_distribution_window = false;
    } statistics;

    // --- FRAMEBUFFER ---
    MainFramebuffer fbo;
    unsigned int picking_idx = NO_PICKING_IDX;

    // --- PLAYBACK ---
    float64 time = 0.f;  // needs to be double precision for long trajectories
    float frames_per_second = 10.f;
    bool is_playing = false;
    PlaybackInterpolationMode interpolation = PlaybackInterpolationMode::CUBIC_PERIODIC;

    // --- TIME LINE FILTERING ---
    struct {
        bool enabled = true;
        vec2 range{0, 0};
        bool dynamic_window = false;
        float window_extent = 10.f;
    } time_filter;

    // --- VISUALS ---
    // SSAO
    struct {
        bool enabled = false;
        float intensity = 2.0f;
        float radius = 6.0f;
    } ssao;

    // HYDROGEN BONDS
    struct {
        bool enabled = false;
        bool dirty = true;
        vec4 color = vec4(1, 0, 1, 1);
        float distance_cutoff = HYDROGEN_BOND_DISTANCE_CUTOFF_DEFAULT;  // In ?gstr?m
        float angle_cutoff = HYDROGEN_BOND_ANGLE_CUTOFF_DEFAULT;        // In Degrees
        DynamicArray<HydrogenBond> bonds{};
    } hydrogen_bonds;

    // SIMULATION BOX
    struct {
        bool enabled = false;
        vec4 color = vec4(1, 1, 0, 0.5);
    } simulation_box;

    // VOLUME
    struct {
        bool enabled = false;
        vec3 color = vec3(1, 0, 0);
        float density_scale = 1.f;

        struct {
            GLuint id = 0;
            bool dirty = false;
            ivec3 dim = ivec3(0);
            float max_value = 1.f;
        } texture;

        Volume volume{};
        std::mutex volume_data_mutex{};

        mat4 model_to_world_matrix{};
        mat4 texture_to_model_matrix{};
        mat4 world_to_texture_matrix{};
    } density_volume;

    // RAMACHANDRAN
    struct {
        bool show_window = false;
        float radius = 1.f;
        float opacity = 1.f;
        int frame_range_min = 0;
        int frame_range_max = 0;

        BackboneAnglesTrajectory backbone_angles{};
        Array<BackboneAngles> current_backbone_angles{};
    } ramachandran;

    // DYNAMIC REFERENCE FRAME
    struct {
        Range atom_range{0, 0};
        DynamicArray<bool> atom_mask;

        int32 frame_index = 10;
        mat4 world_to_reference{1};
        mat4 reference_to_world{1};

        bool use_rbf_refinement = false;
        bool view_in_reference = false;
        RadialBasis refinement_basis;

        bool show_grid = false;
        bool show_grid_points = false;
        bool show_error = false;

        float beta = 0.5f;
        bool dirty_flag = true;
    } dynamic_frame;

    struct {
        spatialhash::Frame frame = {};
        vec3 cell_ext = vec3(4.f);  // in Ångström
        bool dirty_flag = true;
    } spatial_hash;

    // --- CONSOLE ---
    Console console{};
    bool show_console = false;

    bool high_res_font = false;
};

namespace dynamic_structure {

enum class StructureType { Undefined = 0, Vector, Plane, CoordinateSystem };

struct DynamicStructure {
    StringBuffer<32> identifier = {};
    StructureType type = StructureType::Undefined;
};

struct VectorStructure : DynamicStructure {
    Array<vec3> frame_data;
};

struct PlaneStructure : DynamicStructure {
    Array<vec4> frame_data;
};

struct CoordinateSystemStructure : DynamicStructure {
    Array<mat4> frame_data;
};

typedef bool (*StructureComputeFunc)(DynamicStructure* structure, CString argument, const MoleculeDynamic& dynamic);

DynamicArray<DynamicStructure> structures;

}  // namespace dynamic_structure

static void reset_view(ApplicationData* data, bool reposition_camera = true);
static float compute_avg_ms(float dt);
static uint32 get_picking_id(uint32 fbo, int32 x, int32 y);

static void draw_main_menu(ApplicationData* data);
static void draw_representations_window(ApplicationData* data);
static void draw_property_window(ApplicationData* data);
static void draw_timeline_window(ApplicationData* data);
static void draw_distribution_window(ApplicationData* data);
static void draw_ramachandran_window(ApplicationData* data);
static void draw_atom_info_window(const MoleculeStructure& mol, int atom_idx, int x, int y);
static void draw_async_info(ApplicationData* data);

static void init_main_framebuffer(MainFramebuffer* fbo, int width, int height);
static void destroy_main_framebuffer(MainFramebuffer* fbo);

static void free_molecule_data(ApplicationData* data);
static void load_molecule_data(ApplicationData* data, CString file);

static void load_workspace(ApplicationData* data, CString file);
static void save_workspace(ApplicationData* data, CString file);

static void create_default_representation(ApplicationData* data);
static void remove_representation(ApplicationData* data, int idx);
static void reset_representations(ApplicationData* data);
static void clear_representations(ApplicationData* data);

static void create_volume(ApplicationData* data);

static mat4 compute_volume_basis(const mat4& world_to_reference, const mat3& box) {
    mat4 world_to_volume;
    world_to_volume[0] = world_to_reference[0] / box[0][0];
    world_to_volume[1] = world_to_reference[1] / box[1][1];
    world_to_volume[2] = world_to_reference[2] / box[2][2];
    world_to_volume[3] = vec4(-vec3(0.5f), 0);

    return world_to_volume;
}

// Async operations
static void load_trajectory_async(ApplicationData* data);
static void compute_backbone_angles_async(ApplicationData* data);

int main(int, char**) {

    defer { printf("Good bye!"); };

    ApplicationData data;

    //__attribute__((unused)) int pelle = 0;

    // Init logging
    logging::initialize();
    logging::register_backend([](CString str, logging::Severity, void*) { printf("%s\n", str.cstr()); });
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
    if (!platform::initialize(&data.ctx, 1920, 1080, "VIAMD")) {
        printf("Could not initialize platform layer... terminating\n");
        return 1;
    }
    data.ctx.window.vsync = false;

    init_main_framebuffer(&data.fbo, data.ctx.framebuffer.width, data.ctx.framebuffer.height);

    // Init subsystems
    immediate::initialize();
    draw::initialize();
    ramachandran::initialize();
    stats::initialize();
    filter::initialize();
    postprocessing::initialize(data.fbo.width, data.fbo.height);
    volume::initialize();

    // Setup IMGUI style
    ImGui::StyleColorsClassic();

    bool show_demo_window = false;
    const vec4 CLEAR_COLOR = vec4(1, 1, 1, 1);
    const vec4 CLEAR_NORMAL = vec4(0.5f, 0, 0, 0);
    const vec4 CLEAR_INDEX = vec4(1, 1, 1, 1);

#ifdef VIAMD_RELEASE
    allocate_and_parse_pdb_from_string(&data.mol_data.dynamic, CAFFINE_PDB);
    data.mol_data.atom_radii = compute_atom_radii(data.mol_data.dynamic.molecule.atom_elements);
#else
    
stats::create_property("b1", "distance resatom(resname(ALA), 1) com(resname(ALA))");
load_molecule_data(&data, PROJECT_SOURCE_DIR "/data/1ALA-250ns-2500frames.pdb");
data.dynamic_frame.atom_range = {0, 152};
    /*
    stats::create_property("b1", "distance resname(DE3) com(resname(DE3))");
    load_molecule_data(&data, PROJECT_SOURCE_DIR "/data/haofan/for_VIAMD.pdb");
    data.dynamic_frame.atom_range = {0, 1277};
    */
#endif
    reset_view(&data);
    create_default_representation(&data);
    create_volume(&data);

    // Main loop
    while (!data.ctx.window.should_close) {
        platform::Coordinate previous_mouse_coord = data.ctx.input.mouse.coord;
        platform::update(&data.ctx);

        // Try to fix false move on touch
        if (data.ctx.input.mouse.hit[0]) {
            previous_mouse_coord = data.ctx.input.mouse.coord;
        }

        if (data.density_volume.enabled) {
            stats::async_update(
                data.mol_data.dynamic, data.time_filter.range,
                [](void* usr_data) {
                    ApplicationData* data = (ApplicationData*)usr_data;
                    data->density_volume.volume_data_mutex.lock();

                    stats::compute_density_volume_with_basis(
                        &data->density_volume.volume, data->mol_data.dynamic.trajectory, data->time_filter.range, [data](int32 frame_idx) -> mat4 {
                            /*
        const auto p_ref = get_trajectory_positions(data->mol_data.dynamic.trajectory, data->dynamic_frame.frame_index)
                               .sub_array(data->dynamic_frame.atom_range.x, data->dynamic_frame.atom_range.y);
        const auto p_cur = get_trajectory_positions(data->mol_data.dynamic.trajectory, frame_idx)
                               .sub_array(data->dynamic_frame.atom_range.x, data->dynamic_frame.atom_range.y);
        mat4 world_to_volume = compute_transform(p_cur, p_ref);
        const auto box = get_trajectory_frame(data->mol_data.dynamic.trajectory, frame_idx).box;
        world_to_volume[3] += vec4(box * vec3(0.5f), 0);
        world_to_volume[0] /= box[0][0];
        world_to_volume[1] /= box[1][1];
        world_to_volume[2] /= box[2][2];
                                    */

                            const auto box = get_trajectory_frame(data->mol_data.dynamic.trajectory, frame_idx).box;
                            const auto p_ref = get_trajectory_positions(data->mol_data.dynamic.trajectory, data->dynamic_frame.frame_index)
                                                   .sub_array(data->dynamic_frame.atom_range.x, data->dynamic_frame.atom_range.y);
                            const auto p_cur = get_trajectory_positions(data->mol_data.dynamic.trajectory, frame_idx)
                                                   .sub_array(data->dynamic_frame.atom_range.x, data->dynamic_frame.atom_range.y);

                            auto cur_to_ref = compute_linear_transform(p_cur, p_ref);
                            return compute_volume_basis(cur_to_ref, box);
                        });

                    // stats::compute_density_volume(&data->density_volume.volume, data->density_volume.world_to_texture_matrix,
                    //                              data->mol_data.dynamic.trajectory, data->time_filter.range);

                    data->density_volume.volume_data_mutex.unlock();
                    data->density_volume.texture.dirty = true;
                },
                &data);
        } else {
            stats::async_update(data.mol_data.dynamic, data.time_filter.range);
        }

        // If gpu representation of volume is not up to date, upload data
        if (data.density_volume.texture.dirty) {
            if (data.density_volume.volume_data_mutex.try_lock()) {
                if (data.density_volume.texture.dim != data.density_volume.volume.dim) {
                    data.density_volume.texture.dim = data.density_volume.volume.dim;
                    volume::create_volume_texture(&data.density_volume.texture.id, data.density_volume.texture.dim);
                }

                volume::set_volume_texture_data(data.density_volume.texture.id, data.density_volume.texture.dim,
                                                data.density_volume.volume.voxel_data.data);
                data.density_volume.volume_data_mutex.unlock();
                data.density_volume.texture.max_value = data.density_volume.volume.voxel_range.y;
                data.density_volume.texture.dirty = false;
            }
        }

        // RESIZE FRAMEBUFFER?
        if ((data.fbo.width != data.ctx.framebuffer.width || data.fbo.height != data.ctx.framebuffer.height) &&
            (data.ctx.framebuffer.width != 0 && data.ctx.framebuffer.height != 0)) {
            init_main_framebuffer(&data.fbo, data.ctx.framebuffer.width, data.ctx.framebuffer.height);
            postprocessing::initialize(data.fbo.width, data.fbo.height);
        }

        // Setup fbo and clear textures
        glViewport(0, 0, data.fbo.width, data.fbo.height);

        const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2};
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, data.fbo.id);

        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        glDepthMask(GL_TRUE);

        // Clear color, normal and depth buffer
        glDrawBuffers(2, draw_buffers);
        glClearColor(CLEAR_COLOR.x, CLEAR_COLOR.y, CLEAR_COLOR.z, CLEAR_COLOR.w);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Clear picking buffer
        glDrawBuffer(GL_COLOR_ATTACHMENT2);
        glClearColor(CLEAR_INDEX.x, CLEAR_INDEX.y, CLEAR_INDEX.z, CLEAR_INDEX.w);
        glClear(GL_COLOR_BUFFER_BIT);

        // Enable all draw buffers
        glDrawBuffers(3, draw_buffers);

        if (data.ctx.input.key.hit[CONSOLE_KEY]) {
            data.console.visible = !data.console.visible;
        }

        if (data.async.trajectory.sync.running) {
            constexpr float TICK_INTERVAL_SEC = 3.f;
            static float time = 0.f;
            time += data.ctx.timing.delta_s;
            if (time > TICK_INTERVAL_SEC) {
                time = 0.f;
                stats::set_all_property_flags(true, true);
                compute_backbone_angles_async(&data);
            }
        }

        float ms = compute_avg_ms(data.ctx.timing.delta_s);
        bool time_changed = false;
        bool frame_changed = false;

        if (data.is_playing) {
            data.time += data.ctx.timing.delta_s * data.frames_per_second;
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
            float max_frame = data.mol_data.dynamic.trajectory ? data.mol_data.dynamic.trajectory.num_frames : 1.f;
            data.time_filter.range.x = math::max((float)data.time - data.time_filter.window_extent * 0.5f, 0.f);
            data.time_filter.range.y = math::min((float)data.time + data.time_filter.window_extent * 0.5f, max_frame);
        }

        if (data.mol_data.dynamic.trajectory && time_changed) {
            data.dynamic_frame.dirty_flag = true;
            data.spatial_hash.dirty_flag = true;

            int last_frame = data.mol_data.dynamic.trajectory.num_frames - 1;
            data.time = math::clamp(data.time, 0.0, float64(last_frame));
            if (data.time == float64(last_frame)) data.is_playing = false;

            int frame = (int)data.time;
            int prev_frame_2 = math::max(0, frame - 1);
            int prev_frame_1 = math::max(0, frame);
            int next_frame_1 = math::min(frame + 1, last_frame);
            int next_frame_2 = math::min(frame + 2, last_frame);

            if (prev_frame_1 == next_frame_1) {
                copy_trajectory_positions(data.mol_data.dynamic.molecule.atom_positions, data.mol_data.dynamic.trajectory, prev_frame_1);
            } else {
                float t = (float)math::fract(data.time);

                // INTERPOLATE
                switch (data.interpolation) {
                    case PlaybackInterpolationMode::NEAREST:
                        // @NOTE THIS IS ACTUALLY FLOORING
                        copy_trajectory_positions(data.mol_data.dynamic.molecule.atom_positions, data.mol_data.dynamic.trajectory, prev_frame_1);
                        break;
                    case PlaybackInterpolationMode::LINEAR: {
                        auto prev_frame = get_trajectory_frame(data.mol_data.dynamic.trajectory, prev_frame_1);
                        auto next_frame = get_trajectory_frame(data.mol_data.dynamic.trajectory, next_frame_1);
                        linear_interpolation(data.mol_data.dynamic.molecule.atom_positions, prev_frame.atom_positions, next_frame.atom_positions, t);
                        break;
                    }
                    case PlaybackInterpolationMode::LINEAR_PERIODIC: {
                        auto prev_frame = get_trajectory_frame(data.mol_data.dynamic.trajectory, prev_frame_1);
                        auto next_frame = get_trajectory_frame(data.mol_data.dynamic.trajectory, next_frame_1);
                        linear_interpolation_periodic(data.mol_data.dynamic.molecule.atom_positions, prev_frame.atom_positions,
                                                      next_frame.atom_positions, t, prev_frame.box);
                        break;
                    }
                    case PlaybackInterpolationMode::CUBIC: {
                        auto prev_2 = get_trajectory_frame(data.mol_data.dynamic.trajectory, prev_frame_2);
                        auto prev_1 = get_trajectory_frame(data.mol_data.dynamic.trajectory, prev_frame_1);
                        auto next_1 = get_trajectory_frame(data.mol_data.dynamic.trajectory, next_frame_1);
                        auto next_2 = get_trajectory_frame(data.mol_data.dynamic.trajectory, next_frame_2);
                        spline_interpolation(data.mol_data.dynamic.molecule.atom_positions, prev_2.atom_positions, prev_1.atom_positions,
                                             next_1.atom_positions, next_2.atom_positions, t);
                        break;
                    }
                    case PlaybackInterpolationMode::CUBIC_PERIODIC: {
                        auto prev_2 = get_trajectory_frame(data.mol_data.dynamic.trajectory, prev_frame_2);
                        auto prev_1 = get_trajectory_frame(data.mol_data.dynamic.trajectory, prev_frame_1);
                        auto next_1 = get_trajectory_frame(data.mol_data.dynamic.trajectory, next_frame_1);
                        auto next_2 = get_trajectory_frame(data.mol_data.dynamic.trajectory, next_frame_2);
                        spline_interpolation_periodic(data.mol_data.dynamic.molecule.atom_positions, prev_2.atom_positions, prev_1.atom_positions,
                                                      next_1.atom_positions, next_2.atom_positions, t, prev_1.box);
                        break;
                    }

                    default:
                        break;
                }
            }
        }

        if (data.spatial_hash.dirty_flag) {
            data.spatial_hash.dirty_flag = false;

            spatialhash::compute_frame(&data.spatial_hash.frame, data.mol_data.dynamic.molecule.atom_positions, data.spatial_hash.cell_ext);
        }

        if (data.dynamic_frame.dirty_flag) {
            data.dynamic_frame.dirty_flag = false;

            int frame = (int)data.time;
            const auto box = get_trajectory_frame(data.mol_data.dynamic.trajectory, frame).box;
            const auto p_ref = get_trajectory_positions(data.mol_data.dynamic.trajectory, data.dynamic_frame.frame_index)
                                   .sub_array(data.dynamic_frame.atom_range.x, data.dynamic_frame.atom_range.y);
            auto p_cur = data.mol_data.dynamic.molecule.atom_positions.sub_array(data.dynamic_frame.atom_range.x, data.dynamic_frame.atom_range.y);

            mat4 M = compute_linear_transform(p_ref, p_cur);
            mat3 A = M;
            mat3 R, S;
            decompose(A, &R, &S);
            float det = math::abs(math::determinant(A));
            // LOG_NOTE("Determinant is: %.5f", det);
            A = A / math::pow(det, 1.f / 3.f);

            mat4 M2 = data.dynamic_frame.beta * A + (1.f - data.dynamic_frame.beta) * R;
            M2[3] = M[3];
            data.dynamic_frame.reference_to_world = M2;

            // Transform origin to corner instead of center
            // data.dynamic_frame.reference_to_world[3] -= data.dynamic_frame.reference_to_world * vec4(box * vec3(0.5f), 0);
            data.dynamic_frame.world_to_reference = math::inverse(data.dynamic_frame.reference_to_world);
            data.density_volume.model_to_world_matrix = math::inverse(compute_volume_basis(data.dynamic_frame.world_to_reference, box));

            if (data.dynamic_frame.use_rbf_refinement) {
                DynamicArray<vec3> basis_points;
                DynamicArray<vec3> basis_values;
                const auto elements = data.mol_data.dynamic.molecule.atom_elements;
                for (int32 i = 0; i < p_cur.count; i++) {
                    // if (elements[i] == Element::C || elements[i] == Element::O || elements[i] == Element::P) {
                    if (elements[i] == Element::C) {
                        vec3 delta = vec3(data.dynamic_frame.reference_to_world * vec4(p_ref[i], 1.f)) - p_cur[i];
                        basis_points.push_back(p_cur[i]);
                        basis_values.push_back(delta);
                    }
                }

                // DynamicArray<vec3> p_ref_mod(p_ref.count);
                // for (int32 i = 0; i < p_ref_mod.count; i++) {
                //    p_ref_mod[i] = vec3(data.dynamic_frame.reference_to_world * vec4(p_ref[i], 1.f));
                //}

                data.dynamic_frame.refinement_basis = compute_radial_basis(basis_points, basis_values);

                for (auto& p : data.mol_data.dynamic.molecule.atom_positions) {
                    auto v = evaluate_radial_basis(data.dynamic_frame.refinement_basis, p);
                    p = p + v;
                }
            }
        }

        if (frame_changed) {
            data.hydrogen_bonds.dirty = true;
            if (data.mol_data.dynamic.trajectory) {
                if (data.time_filter.dynamic_window) {
                    stats::set_all_property_flags(false, true);
                }
            }
        }

        if (data.hydrogen_bonds.enabled && data.hydrogen_bonds.dirty) {
            data.hydrogen_bonds.bonds.clear();
            hydrogen_bond::compute_bonds(&data.hydrogen_bonds.bonds, data.mol_data.dynamic.molecule.hydrogen_bond.donors,
                                         data.mol_data.dynamic.molecule.hydrogen_bond.acceptors, data.mol_data.dynamic.molecule.atom_positions,
                                         data.hydrogen_bonds.distance_cutoff, data.hydrogen_bonds.angle_cutoff * math::DEG_TO_RAD);
            data.hydrogen_bonds.dirty = false;
        }

        // CAMERA CONTROLS
        if (!ImGui::GetIO().WantCaptureMouse) {
            data.camera.trackball_state.input.rotate_button = data.ctx.input.mouse.down[0];
            data.camera.trackball_state.input.pan_button = data.ctx.input.mouse.down[1];
            data.camera.trackball_state.input.dolly_button = data.ctx.input.mouse.down[2];
            data.camera.trackball_state.input.mouse_coord_prev = {previous_mouse_coord.x, previous_mouse_coord.y};
            data.camera.trackball_state.input.mouse_coord_curr = {data.ctx.input.mouse.coord.x, data.ctx.input.mouse.coord.y};
            data.camera.trackball_state.input.screen_size = vec2(data.ctx.window.width, data.ctx.window.height);
            data.camera.trackball_state.input.dolly_delta = data.ctx.input.mouse.scroll_delta;

            camera_controller_trackball(&data.camera.camera, &data.camera.trackball_state);

            if (data.ctx.input.mouse.release[0] && data.ctx.input.mouse.coord == previous_mouse_coord) {
                data.selected = data.hovered;
            }
        }
        if (!ImGui::GetIO().WantCaptureKeyboard) {
            if (data.ctx.input.key.hit[Key::KEY_SPACE]) data.is_playing = !data.is_playing;
        }

        // RENDER TO FBO
        mat4 view_mat = compute_world_to_view_matrix(data.camera.camera);
        if (data.dynamic_frame.view_in_reference) {
            view_mat = view_mat * data.dynamic_frame.world_to_reference;
        }
        mat4 proj_mat = compute_perspective_projection_matrix(data.camera.camera, data.fbo.width, data.fbo.height);
        mat4 inv_proj_mat = math::inverse(proj_mat);

        for (const auto& rep : data.representations.data) {
            if (!rep.enabled) continue;
            switch (rep.type) {
                case Representation::VDW:
                    draw::draw_vdw(data.mol_data.dynamic.molecule.atom_positions, data.mol_data.atom_radii, rep.colors, view_mat, proj_mat,
                                   rep.radius);
                    break;
                case Representation::LICORICE:
                    draw::draw_licorice(data.mol_data.dynamic.molecule.atom_positions, data.mol_data.dynamic.molecule.covalent_bonds, rep.colors,
                                        view_mat, proj_mat, rep.radius);
                    break;
                case Representation::RIBBONS:
                    draw::draw_ribbons(data.mol_data.dynamic.molecule.backbone_segments, data.mol_data.dynamic.molecule.chains,
                                       data.mol_data.dynamic.molecule.atom_positions, rep.colors, view_mat, proj_mat, rep.num_subdivisions,
                                       rep.tension, rep.width, rep.thickness);
                    break;
            }
        }

        // RENDER DEBUG INFORMATION (WITH DEPTH)
        {
            immediate::set_view_matrix(view_mat);
            immediate::set_proj_matrix(proj_mat);

            if (data.hydrogen_bonds.enabled) {
                for (const auto& bond : data.hydrogen_bonds.bonds) {
                    immediate::draw_line(data.mol_data.dynamic.molecule.atom_positions[bond.acc_idx],
                                         data.mol_data.dynamic.molecule.atom_positions[bond.hyd_idx], ImColor(vec_cast(data.hydrogen_bonds.color)));
                }
            }

            if (data.simulation_box.enabled && data.mol_data.dynamic.trajectory.num_frames > 0) {
                int32 frame_idx = math::clamp((int)data.time, 0, data.mol_data.dynamic.trajectory.num_frames - 1);
                TrajectoryFrame frame = get_trajectory_frame(data.mol_data.dynamic.trajectory, frame_idx);
                immediate::draw_aabb(vec3(0), frame.box * vec3(1), ImColor(vec_cast(data.simulation_box.color)));
            }

            if (data.mol_data.dynamic.trajectory) {
                const mat4& mat = data.dynamic_frame.reference_to_world;
                const TrajectoryFrame frame = get_trajectory_frame(data.mol_data.dynamic.trajectory, 0);
                const vec3 ext = frame.box * vec3(1.0f);

                const vec3 min_val = vec3(0);
                const vec3 max_val = ext;
                const vec3 step = ext / 4.f;

                /*

for (float x = min_val.x; x <= max_val.x; x += step.x) {
    for (float y = min_val.y; y <= max_val.y; y += step.y) {
        immediate::draw_line(mat * vec4(x, y, min_val.z, 1), mat * vec4(x, y, max_val.z, 1), 0xff000000);
    }
}
for (float x = min_val.x; x <= max_val.x; x += step.x) {
    for (float z = min_val.z; z <= max_val.z; z += step.z) {
        immediate::draw_line(mat * vec4(x, min_val.y, z, 1), mat * vec4(x, max_val.y, z, 1), 0xff000000);
    }
}
for (float y = min_val.y; y <= max_val.y; y += step.y) {
    for (float z = min_val.z; z <= max_val.z; z += step.z) {
        immediate::draw_line(mat * vec4(min_val.x, y, z, 1), mat * vec4(max_val.x, y, z, 1), 0xff000000);
    }
}
                */

                const ivec3 RES = {15, 15, 15};
                const vec3 STEP = {ext.x / RES.x, ext.y / RES.y, ext.z / RES.z};

                DynamicArray<vec3> grid_points(RES.x * RES.y * RES.z);
                for (int32 z = 0; z < RES.z; z++) {
                    for (int32 y = 0; y < RES.y; y++) {
                        for (int32 x = 0; x < RES.x; x++) {
                            int32 idx = z * RES.x * RES.y + y * RES.y + x;
                            vec3 p = mat * vec4(min_val + STEP * vec3(x, y, z), 1.f);
                            grid_points[idx] = p;
                            if (data.dynamic_frame.use_rbf_refinement) {
                                auto v = evaluate_radial_basis(data.dynamic_frame.refinement_basis, p);
                                grid_points[idx] += v;
                            }
                            if (data.dynamic_frame.show_grid_points) {
                                immediate::draw_point(grid_points[idx], immediate::COLOR_CYAN);
                            }
                        }
                    }
                }

                if (data.dynamic_frame.show_grid) {
                    for (int32 x = 0; x < RES.x; x++) {
                        for (int32 y = 0; y < RES.y; y++) {
                            for (int32 z = 0; z < RES.z - 1; z++) {
                                int32 i = z * RES.x * RES.y + y * RES.y + x;
                                int32 j = (z + 1) * RES.x * RES.y + y * RES.y + x;
                                immediate::draw_line(grid_points[i], grid_points[j], 0xff000000);
                            }
                        }
                    }

                    for (int32 x = 0; x < RES.x; x++) {
                        for (int32 z = 0; z < RES.z; z++) {
                            for (int32 y = 0; y < RES.y - 1; y++) {
                                int32 i = z * RES.x * RES.y + y * RES.y + x;
                                int32 j = z * RES.x * RES.y + (y + 1) * RES.y + x;
                                immediate::draw_line(grid_points[i], grid_points[j], 0xff000000);
                            }
                        }
                    }

                    for (int32 y = 0; y < RES.y; y++) {
                        for (int32 z = 0; z < RES.z; z++) {
                            for (int32 x = 0; x < RES.x - 1; x++) {
                                int32 i = z * RES.x * RES.y + y * RES.y + x;
                                int32 j = z * RES.x * RES.y + y * RES.y + x + 1;
                                immediate::draw_line(grid_points[i], grid_points[j], 0xff000000);
                            }
                        }
                    }
                }
            }

            if (data.dynamic_frame.show_error) {
                const auto p_ref = get_trajectory_positions(data.mol_data.dynamic.trajectory, data.dynamic_frame.frame_index)
                                       .sub_array(data.dynamic_frame.atom_range.x, data.dynamic_frame.atom_range.y);
                const auto p_cur =
                    data.mol_data.dynamic.molecule.atom_positions.sub_array(data.dynamic_frame.atom_range.x, data.dynamic_frame.atom_range.y);

                for (int32 i = 0; i < p_cur.count; i++) {
                    vec3 v0 = p_cur[i];
                    vec3 v1 = vec3(data.dynamic_frame.reference_to_world * vec4(p_ref[i], 1.f));
                    immediate::draw_line(v0, v1, immediate::COLOR_MAGENTA);
                }
            }

            immediate::flush();
        }

        // PICKING
        {
            ivec2 coord = {data.ctx.input.mouse.coord.x, data.ctx.framebuffer.height - data.ctx.input.mouse.coord.y};
            if (coord.x < 0 || coord.y < 0 || coord.x >= data.ctx.framebuffer.width || coord.y >= data.ctx.framebuffer.height) {
                data.picking_idx = NO_PICKING_IDX;
            } else {
                // data.picking_idx = get_picking_id(data.fbo.id, coord.x, coord.y);
            }

            data.hovered = {};
            if (data.picking_idx != NO_PICKING_IDX) {
                data.hovered.atom_idx = data.picking_idx;
                if (-1 < data.hovered.atom_idx && data.hovered.atom_idx < data.mol_data.dynamic.molecule.atom_residue_indices.count) {
                    data.hovered.residue_idx = data.mol_data.dynamic.molecule.atom_residue_indices[data.hovered.atom_idx];
                }
                if (-1 < data.hovered.residue_idx && data.hovered.residue_idx < data.mol_data.dynamic.molecule.residues.count) {
                    data.hovered.chain_idx = data.mol_data.dynamic.molecule.residues[data.hovered.residue_idx].chain_idx;
                }

                if (data.ctx.input.mouse.hit[1]) {
                    data.right_clicked = data.hovered;
                }
            }
        }

        // Activate backbuffer
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glDrawBuffer(GL_BACK);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glDisable(GL_DEPTH_TEST);
        glDepthMask(GL_FALSE);

        // Render deferred
        postprocessing::render_deferred(data.fbo.tex_depth, data.fbo.tex_color, data.fbo.tex_normal, inv_proj_mat);

        // Apply post processing
        // postprocessing::apply_tonemapping(data.fbo.tex_color);
        if (data.ssao.enabled) {
            postprocessing::apply_ssao(data.fbo.tex_depth, data.fbo.tex_normal, proj_mat, data.ssao.intensity, data.ssao.radius);
        }

        if (data.density_volume.enabled) {
            const float scl = 1.f * data.density_volume.density_scale / data.density_volume.texture.max_value;
            volume::render_volume_texture(data.density_volume.texture.id, data.fbo.tex_depth, data.density_volume.texture_to_model_matrix,
                                          data.density_volume.model_to_world_matrix, view_mat, proj_mat, data.density_volume.color, scl);
        }

        // DRAW DEBUG GRAPHICS W/O DEPTH
        {
            immediate::set_view_matrix(view_mat);
            immediate::set_proj_matrix(proj_mat);
            stats::visualize(data.mol_data.dynamic);

            immediate::draw_basis(data.dynamic_frame.reference_to_world, 5.f);

            /*
for (const auto& p : ps.particle_positions) {
    immediate::draw_point(p, immediate::COLOR_RED);
}
for (const auto& d_c : ps.constraints.distance) {
    const vec3& p0 = ps.particle_positions[d_c.idx[0]];
    const vec3& p1 = ps.particle_positions[d_c.idx[1]];
    immediate::draw_line(p0, p1, immediate::COLOR_GREEN);
}
            */

            immediate::flush();

            render::draw_spatial_hash_cells(data.spatial_hash.frame, view_mat, proj_mat);
        }

        // GUI ELEMENTS
        data.console.Draw("VIAMD", data.ctx.window.width, data.ctx.window.height, data.ctx.timing.delta_s);

        draw_main_menu(&data);

        if (data.representations.show_window) draw_representations_window(&data);
        if (data.statistics.show_property_window) draw_property_window(&data);
        if (data.statistics.show_timeline_window) draw_timeline_window(&data);
        if (data.statistics.show_distribution_window) draw_distribution_window(&data);

        if (data.ramachandran.show_window) {
            draw_ramachandran_window(&data);
        }

        if (!ImGui::GetIO().WantCaptureMouse) {
            if (data.picking_idx != NO_PICKING_IDX) {
                draw_atom_info_window(data.mol_data.dynamic.molecule, data.picking_idx, data.ctx.input.mouse.coord.x, data.ctx.input.mouse.coord.y);
            }
        }

        draw_async_info(&data);

        // MISC WINDOW
        ImGui::Begin("Misc");
        ImGui::Text("%.2f ms (%.1f fps)", ms, 1000.f / (ms));
        ImGui::Checkbox("Show Demo Window", &show_demo_window);
        if (ImGui::Button("Reset View")) {
            reset_view(&data);
        }
        if (data.mol_data.dynamic.trajectory) {
            int32 num_frames = data.mol_data.dynamic.trajectory.num_frames;
            ImGui::Text("Num Frames: %i", num_frames);
            float t = (float)data.time;
            if (ImGui::SliderFloat("Time", &t, 0, (float)(num_frames - 1))) {
                time_changed = true;
                data.time = t;
            }
            ImGui::SliderFloat("Frames Per Second", &data.frames_per_second, 0.1f, 100.f, "%.3f", 4.f);
            if (data.is_playing) {
                if (ImGui::Button("Pause")) data.is_playing = false;
            } else {
                if (ImGui::Button("Play")) data.is_playing = true;
            }
            ImGui::SameLine();
            if (ImGui::Button("Stop")) {
                data.is_playing = false;
                data.time = 0.0;
                time_changed = true;
            }
            ImGui::Combo("type", (int*)(&data.interpolation), "Nearest\0Linear\0Linear Periodic\0Cubic\0Cubic Periodic\0\0");
            ImGui::Checkbox("Dynamic Framewindow", &data.time_filter.dynamic_window);
            if (data.time_filter.dynamic_window) {
                ImGui::SliderFloat("Window Extent", &data.time_filter.window_extent, 1.f, (float)num_frames);
            }
        }
        ImGui::End();

        /*
if (data.right_clicked.atom_idx != NO_PICKING_IDX && data.ctx.input.mouse.hit[1]) {
    // data.text_field_target = 0;
    // if (ImGui::GetIO().WantTextInput) {
    //    data.text_field_target = ImGui::GetActiveID();
    //}
    ImGui::OpenPopup("AtomContextMenu");
}
        */

        // Show the ImGui demo window. Most of the sample code is in ImGui::ShowDemoWindow().
        if (show_demo_window) {
            ImGui::SetNextWindowPos(ImVec2(650, 20),
                                    ImGuiCond_FirstUseEver);  // Normally user code doesn't need/want to call this because positions are saved in .ini
                                                              // file anyway. Here we just want to make the demo initial state a bit more friendly!
            ImGui::ShowDemoWindow(&show_demo_window);
        }
        ImGui::Render();

        // Swap buffers
        platform::swap_buffers(&data.ctx);
    }

    data.async.trajectory.sync.signal_stop_and_wait();
    stats::signal_stop_and_wait();
    data.async.backbone_angles.sync.signal_stop_and_wait();

    destroy_main_framebuffer(&data.fbo);

    platform::shutdown(&data.ctx);

    return 0;
}

// ### MISC FUNCTIONS ###
static float compute_avg_ms(float dt) {
    // @NOTE: Perhaps this can be done with a simple running mean?
    constexpr float interval = 0.5f;
    static float avg = 0.f;
    static int num_frames = 0;
    static float t = 0;
    t += dt;
    num_frames++;

    if (t > interval) {
        avg = t / num_frames * 1000.f;
        t = 0;
        num_frames = 0;
    }

    return avg;
}

static void reset_view(ApplicationData* data, bool reposition_camera) {
    ASSERT(data);
    if (!data->mol_data.dynamic.molecule) return;

    vec3 min_box, max_box;
    compute_bounding_box(&min_box, &max_box, data->mol_data.dynamic.molecule.atom_positions);
    vec3 size = max_box - min_box;
    vec3 cent = (min_box + max_box) * 0.5f;
    vec3 pos = cent + size * 4.f;

    if (reposition_camera) {
        data->camera.camera.position = pos;
        data->camera.trackball_state.distance = math::length(pos - cent);
        look_at(&data->camera.camera.position, &data->camera.camera.orientation, cent, vec3(0, 1, 0));
    }
    data->camera.camera.near_plane = 1.f;
    data->camera.camera.far_plane = math::length(size) * 50.f;
}

static uint32 get_picking_id(uint32 fbo_id, int32 x, int32 y) {
    unsigned char color[4];
    glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo_id);
    glReadBuffer(GL_COLOR_ATTACHMENT2);
    glReadPixels(x, y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, color);
    glReadBuffer(GL_NONE);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    return color[0] + (color[1] << 8) + (color[2] << 16) + (color[3] << 24);
}

// ### DRAW WINDOWS ###
static void draw_main_menu(ApplicationData* data) {
    ASSERT(data);
    bool new_clicked = false;
    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("New", "CTRL+N")) new_clicked = true;
            if (ImGui::MenuItem("Load Data", "CTRL+L")) {
                auto res = platform::file_dialog(platform::FileDialogFlags_Open, {}, "pdb;gro;xtc");
                if (res.result == platform::FileDialogResult::FILE_OK) {
                    load_molecule_data(data, res.path);
                    if (data->representations.data.count > 0) {
                        reset_representations(data);
                    } else {
                        create_default_representation(data);
                    }
                    stats::clear_all_properties();
                    reset_view(data);
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
                            snprintf(res.path.buffer + strnlen(res.path.buffer, res.path.MAX_LENGTH), res.path.MAX_LENGTH, ".%s", FILE_EXTENSION);
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
                        snprintf(res.path.buffer + strnlen(res.path.buffer, res.path.MAX_LENGTH), res.path.MAX_LENGTH, ".%s", FILE_EXTENSION);
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

            // SSAO
            ImGui::BeginGroup();
            ImGui::Checkbox("SSAO", &data->ssao.enabled);
            if (data->ssao.enabled) {
                ImGui::SliderFloat("Intensity", &data->ssao.intensity, 0.5f, 6.f);
                ImGui::SliderFloat("Radius", &data->ssao.radius, 1.f, 30.f);
            }
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
                if (ImGui::ColorEdit4("PointColor", (float*)&color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel))
                    style->point_colors[i] = ImColor(color);
                ImGui::PopID();
            }
            ImGui::Text("line color   ");
            ImGui::SameLine();
            ImVec4 color = ImColor(style->line_color);
            if (ImGui::ColorEdit4("LineColor", (float*)&color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel))
                style->line_color = ImColor(color);
            ImGui::EndGroup();
            ImGui::Separator();

            ImGui::BeginGroup();
            if (ImGui::Checkbox("Use high-res font", &data->high_res_font)) {
                if (ImGui::GetIO().Fonts->Fonts.size() > 1 && data->high_res_font) {
                    ImGui::GetIO().FontDefault = ImGui::GetIO().Fonts->Fonts[1];
                    ImGui::GetIO().FontGlobalScale = 0.75f;
                } else {
                    ImGui::GetIO().FontDefault = ImGui::GetIO().Fonts->Fonts[0];
                    ImGui::GetIO().FontGlobalScale = 1.0f;
                }
            }
            ImGui::EndGroup();
            ImGui::Separator();

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
            ImGui::Checkbox("Ramachandran", &data->ramachandran.show_window);
            ImGui::Checkbox("Properties", &data->statistics.show_property_window);
            ImGui::Checkbox("Timelines", &data->statistics.show_timeline_window);
            ImGui::Checkbox("Distributions", &data->statistics.show_distribution_window);

            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("Test")) {
            if (ImGui::Checkbox("RBF-refinement", &data->dynamic_frame.use_rbf_refinement)) {
                data->dynamic_frame.dirty_flag = true;
            }
            ImGui::Checkbox("View in Reference", &data->dynamic_frame.view_in_reference);
            ImGui::Checkbox("Show error", &data->dynamic_frame.show_error);
            ImGui::Checkbox("Show grid", &data->dynamic_frame.show_grid);
            ImGui::Checkbox("Show grid points", &data->dynamic_frame.show_grid_points);
            if (ImGui::SliderFloat("Beta", &data->dynamic_frame.beta, 0.f, 1.f)) {
                data->dynamic_frame.dirty_flag = true;
            }

            ImGui::EndMenu();
        }

        ImGui::EndMainMenuBar();
    }

    if (new_clicked) ImGui::OpenPopup("Warning New");
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
}

static void draw_representations_window(ApplicationData* data) {
    ImGui::Begin("Representations", &data->representations.show_window, ImGuiWindowFlags_NoFocusOnAppearing);

    if (ImGui::Button("create new")) {
        create_default_representation(data);
    }
    ImGui::SameLine();
    ImGui::PushStyleColor(ImGuiCol_Button, DEL_BTN_COLOR);
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, DEL_BTN_HOVER_COLOR);
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, DEL_BTN_ACTIVE_COLOR);
    if (ImGui::Button("remove all")) {
        data->representations.data.clear();
    }
    ImGui::PopStyleColor(3);
    ImGui::Spacing();
    ImGui::Separator();
    for (int i = 0; i < data->representations.data.count; i++) {
        bool recompute_colors = false;
        auto& rep = data->representations.data[i];
        const float item_width = math::clamp(ImGui::GetWindowContentRegionWidth() - 90.f, 100.f, 300.f);
        StringBuffer<128> name;
        snprintf(name.buffer, name.MAX_LENGTH, "%s###ID", rep.name.buffer);

        ImGui::PushID(i);
        if (ImGui::CollapsingHeader(name.buffer)) {
            ImGui::Checkbox("enabled", &rep.enabled);
            ImGui::SameLine();
            ImGui::PushStyleColor(ImGuiCol_Button, DEL_BTN_COLOR);
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, DEL_BTN_HOVER_COLOR);
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, DEL_BTN_ACTIVE_COLOR);
            if (ImGui::Button("remove")) {
                remove_representation(data, i);
            }
            ImGui::PopStyleColor(3);
            ImGui::SameLine();
            if (ImGui::Button("clone")) {
                Representation clone = rep;
                clone.colors = {(uint32*)MALLOC(rep.colors.size_in_bytes()), rep.colors.count};
                memcpy(clone.colors.data, rep.colors.data, rep.colors.size_in_bytes());
                data->representations.data.insert(&rep, clone);
            }

            ImGui::PushItemWidth(item_width);
            ImGui::InputText("name", rep.name.buffer, rep.name.MAX_LENGTH);
            if (!rep.filter_is_ok) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
            if (ImGui::InputText("filter", rep.filter.buffer, rep.filter.MAX_LENGTH, ImGuiInputTextFlags_EnterReturnsTrue)) {
                recompute_colors = true;
            }
            if (!rep.filter_is_ok) ImGui::PopStyleColor();
            ImGui::Combo("type", (int*)(&rep.type), "VDW\0Licorice\0Ribbons\0\0");
            if (ImGui::Combo("color mapping", (int*)(&rep.color_mapping), "Static Color\0CPK\0Res Id\0Res Idx\0Chain Id\0Chain Idx\0\0")) {
                recompute_colors = true;
            }
            ImGui::PopItemWidth();
            if (rep.color_mapping == ColorMapping::STATIC_COLOR) {
                ImGui::SameLine();
                if (ImGui::ColorEdit4("color", (float*)&rep.static_color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel)) {
                    recompute_colors = true;
                }
            }
            ImGui::PushItemWidth(item_width);
            if (rep.type == Representation::VDW || rep.type == Representation::LICORICE) {
                ImGui::SliderFloat("radii scale", &rep.radius, 0.1f, 2.f);
            }
            if (rep.type == Representation::RIBBONS) {
                ImGui::SliderInt("spline subdivisions", &rep.num_subdivisions, 1, 16);
                ImGui::SliderFloat("spline tension", &rep.tension, 0.f, 1.f);
                ImGui::SliderFloat("spline width", &rep.width, 0.1f, 2.f);
                ImGui::SliderFloat("spline thickness", &rep.thickness, 0.1f, 2.f);
            }
            ImGui::PopItemWidth();
            ImGui::Spacing();
            ImGui::Separator();
        }

        // ENABLE DRAGGING TO REORDER ELEMENTS
        // THIS IS BROKEN BECAUSE OF THE SPACE BETWEEN ELEMENTS

        /*
if (ImGui::GetActiveID() == ImGui::GetID(name.buffer) && !ImGui::IsItemHovered()) {
    float drag_dy = ImGui::GetMouseDragDelta(0).y;
    if (drag_dy < 0.0f && i > 0) {
        // Swap
        Representation tmp = data->representations.data[i];
        data->representations.data[i] = data->representations.data[i - 1];
        data->representations.data[i - 1] = tmp;
        ImGui::ResetMouseDragDelta();
    } else if (drag_dy > 0.0f && i < data->representations.data.count - 1) {
        Representation tmp = data->representations.data[i];
        data->representations.data[i] = data->representations.data[i + 1];
        data->representations.data[i + 1] = tmp;
        ImGui::ResetMouseDragDelta();
    }
}
        */
        ImGui::PopID();

        if (recompute_colors) {
            compute_atom_colors(rep.colors, data->mol_data.dynamic.molecule, rep.color_mapping,
                                ImGui::ColorConvertFloat4ToU32(vec_cast(rep.static_color)));
            DynamicArray<bool> mask(data->mol_data.dynamic.molecule.atom_elements.count, false);
            rep.filter_is_ok = filter::compute_filter_mask(mask, data->mol_data.dynamic, rep.filter.buffer);
            filter::filter_colors(rep.colors, mask);
        }
    }

    ImGui::End();
}

static void draw_property_window(ApplicationData* data) {
    static bool first_time_shown = true;
    ImGui::Begin("Properties", &data->statistics.show_property_window, ImGuiWindowFlags_NoFocusOnAppearing);

    ImGui::PushID("PROPERTIES");
    ImGui::PushItemWidth(-1);
    ImGui::Columns(4, "columns", true);
    ImGui::Separator();

    if (first_time_shown) {
        first_time_shown = false;
        ImGui::SetColumnWidth(0, ImGui::GetWindowContentRegionWidth() * 0.1f);
        ImGui::SetColumnWidth(1, ImGui::GetWindowContentRegionWidth() * 0.7f);
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
        if (ImGui::InputText("##name", prop->name_buf.buffer, prop->name_buf.MAX_LENGTH, ImGuiInputTextFlags_EnterReturnsTrue)) {
            prop->data_dirty = true;
            // compute_stats = true;
        }
        if (!prop->valid) ImGui::PopStyleColor();
        ImGui::PopItemWidth();
        ImGui::NextColumn();

        char insert_buf[32];
        bool insert_buf_set = false;
        if (ImGui::BeginPopup("AtomContextMenu")) {
            char buf[32];
            if (data->right_clicked.atom_idx != -1) {
                snprintf(buf, 32, "atom(%i) ", data->right_clicked.atom_idx + 1);
                if (ImGui::MenuItem(buf)) {
                    memcpy(insert_buf, buf, 32);
                    insert_buf_set = true;
                }
            }
            if (data->right_clicked.residue_idx != -1) {
                snprintf(buf, 32, "residue(%i) ", data->right_clicked.residue_idx + 1);
                if (ImGui::MenuItem(buf)) {
                    memcpy(insert_buf, buf, 32);
                    insert_buf_set = true;
                }
                snprintf(buf, 32, "resid(%i) ", data->mol_data.dynamic.molecule.residues[data->right_clicked.residue_idx].id);
                if (ImGui::MenuItem(buf)) {
                    memcpy(insert_buf, buf, 32);
                    insert_buf_set = true;
                }
                snprintf(buf, 32, "resname(%s) ", data->mol_data.dynamic.molecule.residues[data->right_clicked.residue_idx].name.cstr());
                if (ImGui::MenuItem(buf)) {
                    memcpy(insert_buf, buf, 32);
                    insert_buf_set = true;
                }
                /*
                                // This does not work as the internal buffer of imgui has a length of 17 characters.
                                if (ImGui::BeginMenu("resatom...")) {
                                        snprintf(buf, 32, "resatom(resid(%i), %i) ",
data->mol_data.dynamic.molecule.residues[data->right_clicked.residue_idx].id, data->right_clicked.atom_idx + 1); if (ImGui::MenuItem(buf)) {
                                                memcpy(insert_buf, buf, 32);
                                                insert_buf_set = true;
                                        }
                                        snprintf(buf, 32, "resatom(resname(%s), %i) ",
                                                         data->mol_data.dynamic.molecule.residues[data->right_clicked.residue_idx].name.cstr(),
data->right_clicked.atom_idx + 1); if (ImGui::MenuItem(buf)) { memcpy(insert_buf, buf, 32); insert_buf_set = true;
                                        }
                                        ImGui::EndMenu();
}
                */
            }
            if (data->right_clicked.chain_idx != -1) {
                snprintf(buf, 32, "chain(%i) ", data->right_clicked.chain_idx + 1);
                if (ImGui::MenuItem(buf)) {
                    memcpy(insert_buf, buf, 32);
                    insert_buf_set = true;
                }
            }
            ImGui::EndPopup();
        }

        if (insert_buf_set) {
            ImGui::SetActiveID(ImGui::GetID("##args"), ImGui::GetCurrentWindow());
            ImGui::SetKeyboardFocusHere();
            ImGui::GetIO().AddInputCharactersUTF8(insert_buf);
        }

        ImGui::PushItemWidth(-1);
        if (!prop->valid) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
        if (ImGui::InputText("##args", prop->args_buf.buffer, prop->args_buf.MAX_LENGTH, ImGuiInputTextFlags_EnterReturnsTrue)) {
            prop->data_dirty = true;
        }
        if (!prop->valid) ImGui::PopStyleColor();
        ImGui::PopItemWidth();

        if (ImGui::IsItemActive() && data->hovered.atom_idx != -1 && data->ctx.input.mouse.release[1]) {
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
        ImGui::PushStyleColor(ImGuiCol_Button, DEL_BTN_COLOR);
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, DEL_BTN_HOVER_COLOR);
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, DEL_BTN_ACTIVE_COLOR);
        if (ImGui::Button("remove")) {
            stats::remove_property(prop);
        }
        ImGui::PopStyleColor(3);
        ImGui::NextColumn();
        ImGui::PopID();
    }
    ImGui::Columns(1);
    ImGui::Separator();
    ImGui::PopID();
    ImGui::PopItemWidth();

    if (ImGui::Button("create new")) {
        stats::create_property();
    }
    ImGui::SameLine();
    ImGui::PushStyleColor(ImGuiCol_Button, DEL_BTN_COLOR);
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, DEL_BTN_HOVER_COLOR);
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, DEL_BTN_ACTIVE_COLOR);
    if (ImGui::Button("remove all")) {
        stats::remove_all_properties();
    }
    ImGui::PopStyleColor(3);
    ImGui::End();

    // if (compute_stats) {
    // stats::compute_stats(data->mol_data.dynamic);
    // compute_statistics_async(data);
    //}
}

static void draw_atom_info_window(const MoleculeStructure& mol, int atom_idx, int x, int y) {

    // @TODO: Assert things and make this failproof
    if (atom_idx < 0 || atom_idx >= mol.atom_positions.count) return;

    int res_idx = mol.atom_residue_indices[atom_idx];
    const Residue& res = mol.residues[res_idx];
    const char* res_id = res.name;
    int local_idx = atom_idx - res.beg_atom_idx;
    const char* label = mol.atom_labels[atom_idx];
    const char* elem = element::name(mol.atom_elements[atom_idx]);
    const char* symbol = element::symbol(mol.atom_elements[atom_idx]);

    int chain_idx = res.chain_idx;
    const char* chain_id = 0;
    if (chain_idx != -1) {
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
    int len = snprintf(buff, 256, "atom[%i][%i]: %s %s %s\nres[%i]: %s", atom_idx, local_idx, label, elem, symbol, res_idx, res_id);
    if (chain_idx) {
        snprintf(buff + len, 256 - len, "\nchain[%i]: %s\n", chain_idx, chain_id);
    }

    ImVec2 text_size = ImGui::CalcTextSize(buff);
    ImGui::SetNextWindowPos(ImVec2(x + 10.f, y + 10.f));
    ImGui::SetNextWindowSize(ImVec2(text_size.x + 20.f, text_size.y + 15.f));
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
    ImGui::Begin("##Atom Info", 0,
                 ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar |
                     ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoInputs | ImGuiWindowFlags_NoBringToFrontOnFocus |
                     ImGuiWindowFlags_NoFocusOnAppearing);
    ImGui::Text("%s", buff);
    ImGui::End();
    ImGui::PopStyleColor();
}

static void draw_async_info(ApplicationData* data) {
    constexpr float WIDTH = 300.f;
    constexpr float MARGIN = 10.f;
    constexpr float PROGRESS_FRACT = 0.3f;

    float traj_fract = data->async.trajectory.fraction;
    float angle_fract = data->async.backbone_angles.fraction;
    float stats_fract = stats::fraction_done();

    if ((0.f < traj_fract && traj_fract < 1.f) || (0.f < angle_fract && angle_fract < 1.f) || (0.f < stats_fract && stats_fract < 1.f)) {

        ImGui::SetNextWindowPos(ImVec2(data->ctx.window.width - WIDTH - MARGIN,
                                       ImGui::GetCurrentContext()->FontBaseSize + ImGui::GetStyle().FramePadding.y * 2.f + MARGIN));
        ImGui::SetNextWindowSize(ImVec2(WIDTH, 0));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
        ImGui::Begin("##Async Info", 0,
                     ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar |
                         ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoFocusOnAppearing);

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
    ImGui::SetNextWindowSize(ImVec2(400, 150), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Timelines", &data->statistics.show_timeline_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
        static float zoom = 1.f;
        ImGui::BeginChild("Scroll Region", ImVec2(0, 0), true, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_HorizontalScrollbar);

        const int max_frame = data->mol_data.dynamic.trajectory.num_frames;
        const Range frame_range(0, (float)max_frame);
        auto old_range = data->time_filter.range;

        ImGui::PushItemWidth(ImGui::GetWindowContentRegionWidth() * zoom);
        if (ImGui::RangeSliderFloat("###selection_range", &data->time_filter.range.x, &data->time_filter.range.y, 0.f, (float)max_frame)) {
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
        // const float plot_height = ImGui::GetContentRegionAvail().y / (float)prop_count;
        const float plot_height = 100.f;
        const uint32 bar_fill_color = ImColor(1.f, 1.f, 1.f, 0.25f);
        const uint32 var_fill_color = ImColor(1.f, 1.f, 0.3f, 0.1f);
        const uint32 var_line_color = ImColor(1.f, 1.f, 0.3f, 0.3f);
        const uint32 var_text_color = ImColor(1.f, 1.f, 0.3f, 0.5f);

        static float selection_start;
        static bool is_selecting = false;

        const auto properties = stats::get_properties();
        for (int i = 0; i < (int32)properties.count; i++) {
            auto prop = properties[i];

            if (!prop->enable_timeline) continue;
            Array<float> prop_data = prop->avg_data;
            CString prop_name = prop->name_buf;
            Range prop_range = prop->avg_data_range;
            if (!prop_data) continue;
            float pad = math::abs(prop_range.y - prop_range.x) * 0.75f;
            vec2 display_range = prop_range + vec2(-pad, pad);
            if (display_range.x == display_range.y) {
                display_range.x -= 1.f;
                display_range.y += 1.f;
            }
            // float val = (float)data->time;
            ImGuiID id = ImGui::GetID(prop_name);

            ImGui::PushID(i);

            ImGui::BeginPlot(prop_name, ImVec2(0, plot_height), ImVec2(frame_range.x, frame_range.y), ImVec2(display_range.x, display_range.y),
                             ImGui::LinePlotFlags_AxisX);
            const ImRect inner_bb(ImGui::GetItemRectMin() + ImGui::GetStyle().FramePadding, ImGui::GetItemRectMax() - ImGui::GetStyle().FramePadding);
            ImGui::PushClipRect(ImGui::GetItemRectMin(), ImGui::GetItemRectMax(), true);

            if (ImGui::IsItemHovered()) ImGui::SetHoveredID(id);

            ImGui::PlotVerticalBars(prop->filter_fraction.data, (int32)prop->filter_fraction.count, bar_fill_color);
            if (prop->std_dev_data.data[0] > 0.f) {
                ImGui::PlotVariance(prop->avg_data.data, prop->std_dev_data.data, (int32)prop->std_dev_data.count, 1.f, var_line_color,
                                    var_fill_color);
            }
            ImGui::PlotValues(prop->name_buf.cstr(), prop_data.data, (int32)prop_data.count);

            ImGui::PopClipRect();

            if (ImGui::IsItemHovered() && ImGui::GetIO().MouseClicked[0]) {
                ImGui::SetActiveID(id, ImGui::GetCurrentWindow());
            }

            if (ImGui::IsItemHovered() && ImGui::GetIO().MouseClicked[1] && ImGui::GetIO().KeyCtrl) {
                data->time_filter.range = frame_range;
            }

            if (ImGui::GetActiveID() == id) {
                if (ImGui::GetIO().MouseClicked[0] && ImGui::GetIO().KeyCtrl) {
                    float t = (ImGui::GetIO().MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x);
                    selection_start = ImLerp(frame_range.x, frame_range.y, t);
                    data->time_filter.range.x = selection_start;
                    data->time_filter.range.y = selection_start;
                    is_selecting = true;
                } else if (is_selecting) {
                    float t = (ImGui::GetIO().MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x);
                    float v = ImLerp(frame_range.x, frame_range.y, t);
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
                    float t = ImClamp((ImGui::GetIO().MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x), 0.f, 1.f);
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
                const float t0 = (data->time_filter.range.x - frame_range.x) / (frame_range.y - frame_range.x);
                const float t1 = (data->time_filter.range.y - frame_range.x) / (frame_range.y - frame_range.x);
                const ImVec2 pos0 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t0, 0));
                const ImVec2 pos1 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t1, 1));
                ImGui::GetCurrentWindow()->DrawList->AddRectFilled(pos0, pos1, SELECTION_RANGE_COLOR);
            }

            // CURRENT FRAME POSITION
            {
                constexpr ImU32 CURRENT_LINE_COLOR = 0xaa33ffff;
                const float t = ((float)data->time - frame_range.x) / (frame_range.y - frame_range.x);
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
                const float min_x = ImGui::GetItemRectMin().x;
                const float max_x = ImGui::GetItemRectMax().x;
                float t = ImClamp((ImGui::GetIO().MousePos.x - min_x) / (max_x - min_x), 0.f, 1.f);
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
            constexpr float ZOOM_SCL = 0.24f;
            float pre_coord = ImGui::GetScrollX() + (ImGui::GetIO().MousePos.x - ImGui::GetWindowPos().x) * zoom;
            zoom = math::clamp(zoom + ZOOM_SCL * ImGui::GetIO().MouseWheel, 1.f, 100.f);
            float post_coord = ImGui::GetScrollX() + (ImGui::GetIO().MousePos.x - ImGui::GetWindowPos().x) * zoom;
            float delta = pre_coord - post_coord;
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
    constexpr float RANGE_SLIDER_HEIGHT = 26.f;
    const float plot_height = ImGui::GetContentRegionAvail().y / (float)properties.count - RANGE_SLIDER_HEIGHT;

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
        if (ImGui::ItemAdd(total_bb, NULL)) {
            ImGui::RenderFrame(frame_bb.Min, frame_bb.Max, ImGui::GetColorU32(ImGuiCol_FrameBg), true, style.FrameRounding);

            // const float max_val = math::max(prop->full_histogram.bin_range.y, prop->filt_histogram.bin_range.y);
            ImGui::PushClipRect(inner_bb.Min, inner_bb.Max, true);
            const float max_val = prop->full_histogram.bin_range.y * 1.5f;
            ImGui::DrawFilledLine(inner_bb.Min, inner_bb.Max, prop->full_histogram.bins.data, (int32)prop->full_histogram.bins.count, max_val,
                                  FULL_LINE_COLOR, FULL_FILL_COLOR);

            ImGui::DrawFilledLine(inner_bb.Min, inner_bb.Max, prop->filt_histogram.bins.data, (int32)prop->filt_histogram.bins.count, max_val,
                                  FILT_LINE_COLOR, FILT_FILL_COLOR);
            // ImGui::PopClipRect();

            // SELECTION RANGE
            {
                const float t0 = (prop->filter.x - prop->total_data_range.x) / (prop->total_data_range.y - prop->total_data_range.x);
                const float t1 = (prop->filter.y - prop->total_data_range.x) / (prop->total_data_range.y - prop->total_data_range.x);
                const ImVec2 pos0 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t0, 0));
                const ImVec2 pos1 = ImLerp(inner_bb.Min, inner_bb.Max, ImVec2(t1, 1));
                ImGui::GetCurrentWindow()->DrawList->AddRectFilled(pos0, pos1, SELECTION_RANGE_COLOR);
            }

            ImGui::PopClipRect();

            if (ImGui::IsItemHovered()) {
                window->DrawList->AddLine(ImVec2(ImGui::GetIO().MousePos.x, inner_bb.Min.y), ImVec2(ImGui::GetIO().MousePos.x, inner_bb.Max.y),
                                          0xffffffff);
                float t = (ImGui::GetIO().MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x);
                int32 count = (int32)prop->full_histogram.bins.count;
                int32 idx = ImClamp((int32)(t * (count - 1)), 0, count - 1);
                float full_val = prop->full_histogram.bins.data[idx];
                float filt_val = prop->filt_histogram.bins.data[idx];
                ImVec2 val_range = vec_cast(prop->filt_histogram.value_range);
                ImGui::BeginTooltip();
                ImGui::Text("%.3f:", ImLerp(val_range.x, val_range.y, t));
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
    constexpr vec2 res(512, 512);
    ImGui::SetNextWindowContentSize(ImVec2(res.x, res.y));
    ImGui::Begin("Ramachandran", &data->ramachandran.show_window, ImGuiWindowFlags_NoFocusOnAppearing);

    int32 num_frames = data->mol_data.dynamic.trajectory ? data->mol_data.dynamic.trajectory.num_frames : 0;

    ImGui::SliderFloat("opacity", &data->ramachandran.opacity, 0.f, 10.f);
    ImGui::SliderFloat("radius", &data->ramachandran.radius, 0.1f, 2.f);
    ImGui::RangeSliderFloat("framerange", &data->time_filter.range.x, &data->time_filter.range.y, 0, (float)math::max(0, num_frames));

    int32 frame = (int32)data->time;
    Array<BackboneAngles> accumulated_angles = get_backbone_angles(data->ramachandran.backbone_angles, (int32)data->time_filter.range.x,
                                                                   (int32)data->time_filter.range.y - (int32)data->time_filter.range.x);
    Array<BackboneAngles> current_angles = get_backbone_angles(data->ramachandran.backbone_angles, frame);

    ramachandran::clear_accumulation_texture();

    const vec4 ordinary_color(1.f, 1.f, 1.f, 0.1f * data->ramachandran.opacity);
    ramachandran::compute_accumulation_texture(accumulated_angles, ordinary_color, data->ramachandran.radius);

    float dim = math::min(ImGui::GetWindowWidth(), ImGui::GetWindowHeight());
    ImVec2 win_pos = ImGui::GetCursorScreenPos();
    ImVec2 canvas_size(dim, dim);
    ImDrawList* dl = ImGui::GetWindowDrawList();

    ImVec2 x0 = win_pos;
    ImVec2 x1(win_pos.x + canvas_size.x, win_pos.y + canvas_size.y);

    dl->ChannelsSplit(3);
    dl->ChannelsSetCurrent(0);
    dl->AddImage((ImTextureID)(intptr_t)ramachandran::get_segmentation_texture(), x0, x1);
    dl->ChannelsSetCurrent(1);
    dl->AddImage((ImTextureID)(intptr_t)ramachandran::get_accumulation_texture(), x0, x1);
    dl->ChannelsSetCurrent(2);
    constexpr float ONE_OVER_TWO_PI = 1.f / (2.f * math::PI);
    for (const auto& angle : current_angles) {
        if (angle.phi == 0 || angle.psi == 0) continue;
        ImVec2 coord(angle.phi * ONE_OVER_TWO_PI + 0.5f, angle.psi * ONE_OVER_TWO_PI + 0.5f);  // [-PI, PI] -> [0, 1]
        coord.y = 1.f - coord.y;
        coord = ImLerp(x0, x1, coord);
        float radius = data->ramachandran.radius * 5.f;
        ImVec2 min_box(math::round(coord.x - radius), math::round(coord.y - radius));
        ImVec2 max_box(math::round(coord.x + radius), math::round(coord.y + radius));
        dl->AddRectFilled(min_box, max_box, 0xff00ffff);
        dl->AddRect(min_box, max_box, 0xff000000);
    }
    dl->ChannelsMerge();

    dl->ChannelsSetCurrent(0);

    ImGui::End();
}

// ### FRAMEBUFFER ###
static void init_main_framebuffer(MainFramebuffer* fbo, int width, int height) {
    ASSERT(fbo);

    bool attach_textures = false;
    if (!fbo->id) {
        glGenFramebuffers(1, &fbo->id);
        attach_textures = true;
    }

    if (!fbo->tex_depth) glGenTextures(1, &fbo->tex_depth);
    if (!fbo->tex_color) glGenTextures(1, &fbo->tex_color);
    if (!fbo->tex_normal) glGenTextures(1, &fbo->tex_normal);
    if (!fbo->tex_picking) glGenTextures(1, &fbo->tex_picking);

    glBindTexture(GL_TEXTURE_2D, fbo->tex_depth);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, fbo->tex_color);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, fbo->tex_normal);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16, width, height, 0, GL_RG, GL_UNSIGNED_SHORT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, fbo->tex_picking);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, 0);

    fbo->width = width;
    fbo->height = height;

    glBindFramebuffer(GL_FRAMEBUFFER, fbo->id);
    if (attach_textures) {
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, fbo->tex_depth, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, fbo->tex_color, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, fbo->tex_normal, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, fbo->tex_picking, 0);
    }

    ASSERT(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

static void destroy_main_framebuffer(MainFramebuffer* fbo) {
    ASSERT(fbo);
    if (fbo->id) glDeleteFramebuffers(1, &fbo->id);
    if (fbo->tex_depth) glDeleteTextures(1, &fbo->tex_depth);
    if (fbo->tex_color) glDeleteTextures(1, &fbo->tex_color);
    if (fbo->tex_normal) glDeleteTextures(1, &fbo->tex_normal);
    if (fbo->tex_picking) glDeleteTextures(1, &fbo->tex_picking);
}

// ### MOLECULE DATA ###
static void free_trajectory_data(ApplicationData* data) {
    if (data->mol_data.dynamic.trajectory) {
        data->async.trajectory.sync.signal_stop_and_wait();
        stats::signal_stop_and_wait();
        data->async.backbone_angles.sync.signal_stop_and_wait();
        close_file_handle(&data->mol_data.dynamic.trajectory);
        free_trajectory(&data->mol_data.dynamic.trajectory);
    }
}

static void free_molecule_data(ApplicationData* data) {
    if (data->mol_data.dynamic.molecule) {
        free_molecule_structure(&data->mol_data.dynamic.molecule);
    }
    if (data->mol_data.dynamic.trajectory) {
        free_trajectory_data(data);
    }
    data->mol_data.atom_radii.clear();
    free_backbone_angles_trajectory(&data->ramachandran.backbone_angles);
    data->ramachandran.backbone_angles = {};
    data->ramachandran.current_backbone_angles = {};
}

static void load_molecule_data(ApplicationData* data, CString file) {
    ASSERT(data);
    if (file.count > 0) {
        data->is_playing = false;
        CString ext = get_file_extension(file);
        printf("'%s'\n", ext.beg());
        if (compare_n(ext, "pdb", 3, true)) {
            free_molecule_data(data);
            free_string(&data->files.molecule);
            free_string(&data->files.trajectory);
            allocate_and_load_pdb_from_file(&data->mol_data.dynamic, file);

            if (!data->mol_data.dynamic.molecule) {
                printf("ERROR! Failed to load pdb file.\n");
                return;
            }

            data->files.molecule = allocate_string(file);
            data->mol_data.atom_radii = compute_atom_radii(data->mol_data.dynamic.molecule.atom_elements);
            if (data->mol_data.dynamic.trajectory) {
                create_volume(data);
                /*
if (data->mol_data.dynamic.trajectory.num_frames > 0) {
    vec3 box_ext = data->mol_data.dynamic.trajectory.frame_buffer[0].box * vec3(1);
    init_volume(&data->density_volume.volume, math::max(ivec3(1), ivec3(box_ext) / VOLUME_DOWNSAMPLE_FACTOR));
    data->density_volume.model_to_world_matrix = volume::compute_model_to_world_matrix(vec3(0), box_ext);
    data->density_volume.texture_to_model_matrix = volume::compute_texture_to_model_matrix(data->density_volume.volume.dim);
}
                */

                init_backbone_angles_trajectory(&data->ramachandran.backbone_angles, data->mol_data.dynamic);
                compute_backbone_angles_trajectory(&data->ramachandran.backbone_angles, data->mol_data.dynamic);
                // stats::update(data->mol_data.dynamic, );
                // compute_statistics_async(data);
            }
        } else if (compare_n(ext, "gro", 3, true)) {
            free_molecule_data(data);
            free_string(&data->files.molecule);
            free_string(&data->files.trajectory);
            allocate_and_load_gro_from_file(&data->mol_data.dynamic.molecule, file);

            if (!data->mol_data.dynamic.molecule) {
                printf("ERROR! Failed to load gro file.\n");
                return;
            }

            data->files.molecule = allocate_string(file);
            data->mol_data.atom_radii = compute_atom_radii(data->mol_data.dynamic.molecule.atom_elements);
        } else if (compare_n(ext, "xtc", 3, true)) {
            if (!data->mol_data.dynamic.molecule) {
                printf("ERROR! Must have molecule loaded before trajectory can be loaded!\n");
            } else {
                if (data->mol_data.dynamic.trajectory) {
                    free_trajectory_data(data);
                }
                if (!load_and_allocate_trajectory(&data->mol_data.dynamic.trajectory, file)) {
                    printf("ERROR! Problem loading trajectory\n");
                    return;
                }
                if (data->mol_data.dynamic.trajectory) {
                    if (data->mol_data.dynamic.trajectory.num_atoms != data->mol_data.dynamic.molecule.atom_positions.count) {
                        printf("ERROR! The number of atoms in the molecule does not match the number of atoms in the trajectory\n");
                        free_trajectory_data(data);
                        return;
                    }
                    data->files.trajectory = allocate_string(file);

                    init_backbone_angles_trajectory(&data->ramachandran.backbone_angles, data->mol_data.dynamic);
                    load_trajectory_async(data);
                }
            }
        } else {
            printf("ERROR! file extension not supported!\n");
        }
    }
}

// ### WORKSPACE ###
static Representation::Type get_rep_type(CString str) {
    if (compare(str, "VDW"))
        return Representation::VDW;
    else if (compare(str, "LICORICE"))
        return Representation::LICORICE;
    else if (compare(str, "RIBBONS"))
        return Representation::RIBBONS;
    else
        return Representation::VDW;
}

static CString get_rep_type_name(Representation::Type type) {
    switch (type) {
        case Representation::VDW:
            return "VDW";
        case Representation::LICORICE:
            return "LICORICE";
        case Representation::RIBBONS:
            return "RIBBONS";
        default:
            return "UNKNOWN";
    }
}

static ColorMapping get_color_mapping(CString str) {
    if (compare(str, "STATIC_COLOR"))
        return ColorMapping::STATIC_COLOR;
    else if (compare(str, "CPK"))
        return ColorMapping::CPK;
    else if (compare(str, "RES_ID"))
        return ColorMapping::RES_ID;
    else if (compare(str, "RES_INDEX"))
        return ColorMapping::RES_INDEX;
    else if (compare(str, "CHAIN_ID"))
        return ColorMapping::CHAIN_ID;
    else if (compare(str, "CHAIN_INDEX"))
        return ColorMapping::CHAIN_INDEX;
    else
        return ColorMapping::CPK;
}

static CString get_color_mapping_name(ColorMapping mapping) {
    switch (mapping) {
        case ColorMapping::STATIC_COLOR:
            return "STATIC_COLOR";
        case ColorMapping::CPK:
            return "CPK";
        case ColorMapping::RES_ID:
            return "RES_ID";
        case ColorMapping::RES_INDEX:
            return "RES_INDEX";
        case ColorMapping::CHAIN_ID:
            return "CHAIN_ID";
        case ColorMapping::CHAIN_INDEX:
            return "CHAIN_INDEX";
        default:
            return "UNKNOWN";
    }
}

static vec4 to_vec4(CString txt, vec4 default_val = vec4(1)) {
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
    CString c_txt = txt;
    CString line;
    while (extract_line(line, c_txt)) {
        if (compare(line, "[Files]")) {
            while (c_txt.beg() != c_txt.end() && c_txt[0] != '[') {
                extract_line(line, c_txt);
                if (compare_n(line, "MoleculeFile=", 13)) {
                    new_molecule_file = get_absolute_path(file, trim(line.substr(13)));
                }
                if (compare_n(line, "TrajectoryFile=", 15)) {
                    new_trajectory_file = get_absolute_path(file, trim(line.substr(15)));
                }
            }
        } else if (compare(line, "[Representation]")) {
            Representation rep{};
            while (c_txt.beg() != c_txt.end() && c_txt[0] != '[') {
                extract_line(line, c_txt);
                if (compare_n(line, "Name=", 5)) rep.name = trim(line.substr(5));
                if (compare_n(line, "Filter=", 7)) rep.filter = trim(line.substr(7));
                if (compare_n(line, "Type=", 5)) rep.type = get_rep_type(trim(line.substr(5)));
                if (compare_n(line, "ColorMapping=", 13)) rep.color_mapping = get_color_mapping(trim(line.substr(13)));
                if (compare_n(line, "Enabled=", 8)) rep.enabled = to_int(trim(line.substr(8))) != 0;
                if (compare_n(line, "StaticColor=", 12)) rep.static_color = to_vec4(trim(line.substr(12)));
                if (compare_n(line, "Radius=", 7)) rep.radius = to_float(trim(line.substr(7)));
                if (compare_n(line, "NumSubdivisions=", 16)) rep.num_subdivisions = to_int(trim(line.substr(16)));
                if (compare_n(line, "Tension=", 8)) rep.tension = to_float(trim(line.substr(8)));
                if (compare_n(line, "Width=", 6)) rep.width = to_float(trim(line.substr(6)));
                if (compare_n(line, "Thickness=", 10)) rep.thickness = to_float(trim(line.substr(10)));
            }
            data->representations.data.push_back(rep);
        } else if (compare(line, "[Property]")) {
            StringBuffer<256> name, args;
            while (c_txt.beg() != c_txt.end() && c_txt[0] != '[') {
                extract_line(line, c_txt);
                if (compare_n(line, "Name=", 5)) name = trim(line.substr(5));
                if (compare_n(line, "Args=", 5)) args = trim(line.substr(5));
            }
            stats::create_property(name, args);
        } else if (compare(line, "[RenderSettings]")) {
            while (c_txt.beg() != c_txt.end() && c_txt[0] != '[') {
                extract_line(line, c_txt);
                if (compare_n(line, "SsaoEnabled=", 12)) data->ssao.enabled = to_int(trim(line.substr(12))) != 0;
                if (compare_n(line, "SsaoIntensity=", 14)) data->ssao.intensity = to_float(trim(line.substr(14)));
                if (compare_n(line, "SsaoRadius=", 11)) data->ssao.radius = to_float(trim(line.substr(11)));
            }
        } else if (compare(line, "[Camera]")) {
            while (c_txt.beg() != c_txt.end() && c_txt[0] != '[') {
                extract_line(line, c_txt);
                if (compare_n(line, "Position=", 9)) {
                    vec3 pos = vec3(to_vec4(trim(line.substr(9))));
                    data->camera.camera.position = pos;
                }
                if (compare_n(line, "Rotation=", 9)) {
                    quat rot = quat(to_vec4(trim(line.substr(9))));
                    data->camera.camera.orientation = rot;
                }
                if (compare_n(line, "Distance=", 9)) data->camera.trackball_state.distance = to_float(trim(line.substr(9)));
            }
        }
    }

    // Store Loaded Molecule File Relative Path
    // (Store Loaded Trajectory File Relative Path)
    // Store Representations
    // Store Groups and Properties
    // Store Rendersettings
    // ...

    if (txt) FREE(txt.data);

    if (data->files.workspace) free_string(&data->files.workspace);
    data->files.workspace = allocate_string(file);

    if (!compare(new_molecule_file, data->files.molecule) && new_molecule_file) {
        load_molecule_data(data, new_molecule_file);
    }

    if (!compare(new_trajectory_file, data->files.trajectory) && new_trajectory_file) {
        load_molecule_data(data, new_trajectory_file);
    }

    reset_view(data, false);
    reset_representations(data);
}

static void save_workspace(ApplicationData* data, CString file) {
    FILE* fptr = fopen(file.beg(), "w");
    if (!fptr) {
        printf("ERROR! Could not save workspace to file '%s'\n", file.beg());
        return;
    }

    // @TODO: Make relative paths
    fprintf(fptr, "[Files]\n");
    fprintf(fptr, "MoleculeFile=%s\n", data->files.molecule ? get_relative_path(file, data->files.molecule).beg() : "");
    fprintf(fptr, "TrajectoryFile=%s\n", data->files.trajectory ? get_relative_path(file, data->files.trajectory).beg() : "");
    fprintf(fptr, "\n");

    // REPRESENTATIONS
    for (const auto& rep : data->representations.data) {
        fprintf(fptr, "[Representation]\n");
        fprintf(fptr, "Name=%s\n", rep.name.beg());
        fprintf(fptr, "Filter=%s\n", rep.filter.beg());
        fprintf(fptr, "Type=%s\n", get_rep_type_name(rep.type).beg());
        fprintf(fptr, "ColorMapping=%s\n", get_color_mapping_name(rep.color_mapping).beg());
        fprintf(fptr, "Enabled=%i\n", rep.enabled);
        if (rep.color_mapping == ColorMapping::STATIC_COLOR) {
            fprintf(fptr, "StaticColor=%i\n", rep.enabled);
        }
        if (rep.type == Representation::VDW || Representation::LICORICE)
            fprintf(fptr, "Radius=%g\n", rep.radius);
        else if (rep.type == Representation::RIBBONS) {
            fprintf(fptr, "NumSubdivisions=%i\n", rep.num_subdivisions);
            fprintf(fptr, "Tension=%g\n", rep.tension);
            fprintf(fptr, "Width=%g\n", rep.width);
            fprintf(fptr, "Thickness=%g\n", rep.thickness);
        }
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
    fprintf(fptr, "SsaoEnabled=%i\n", data->ssao.enabled ? 1 : 0);
    fprintf(fptr, "SsaoIntensity=%g\n", data->ssao.intensity);
    fprintf(fptr, "SsaoRadius=%g\n", data->ssao.radius);
    fprintf(fptr, "\n");

    fprintf(fptr, "[Camera]\n");
    fprintf(fptr, "Position=%g,%g,%g\n", data->camera.camera.position.x, data->camera.camera.position.y, data->camera.camera.position.z);
    fprintf(fptr, "Rotation=%g,%g,%g,%g\n", data->camera.camera.orientation.x, data->camera.camera.orientation.y, data->camera.camera.orientation.z,
            data->camera.camera.orientation.w);
    fprintf(fptr, "Distance=%g\n", data->camera.trackball_state.distance);
    fprintf(fptr, "\n");

    fclose(fptr);

    if (data->files.workspace) free_string(&data->files.workspace);
    data->files.workspace = allocate_string(file);
}

// ### REPRESENTATIONS ###
static void create_default_representation(ApplicationData* data) {
    ASSERT(data);
    if (!data->mol_data.dynamic.molecule) return;

    auto& rep = data->representations.data.push_back({});
    rep.colors.count = data->mol_data.dynamic.molecule.atom_positions.count;
    rep.colors.data = (uint32*)MALLOC(rep.colors.count * sizeof(uint32));
    compute_atom_colors(rep.colors, data->mol_data.dynamic.molecule, rep.color_mapping);
}

static void remove_representation(ApplicationData* data, int idx) {
    ASSERT(idx < data->representations.data.count);
    auto& rep = data->representations.data[idx];
    if (rep.colors) {
        FREE(rep.colors.data);
    }
    data->representations.data.remove(&rep);
}

static void reset_representations(ApplicationData* data) {
    ASSERT(data);

    for (auto& rep : data->representations.data) {
        if (rep.colors) {
            FREE(rep.colors.data);
        }
        rep.colors.count = data->mol_data.dynamic.molecule.atom_positions.count;
        if (rep.colors.count == 0) {
            rep.colors.data = nullptr;
            continue;
        }
        rep.colors.data = (uint32*)MALLOC(rep.colors.count * sizeof(uint32));

        compute_atom_colors(rep.colors, data->mol_data.dynamic.molecule, rep.color_mapping,
                            ImGui::ColorConvertFloat4ToU32(vec_cast(rep.static_color)));
        DynamicArray<bool> mask(data->mol_data.dynamic.molecule.atom_elements.count, false);
        filter::compute_filter_mask(mask, data->mol_data.dynamic, rep.filter.buffer);
        filter::filter_colors(rep.colors, mask);
    }
}

static void clear_representations(ApplicationData* data) {
    ASSERT(data);
    for (auto& rep : data->representations.data) {
        if (rep.colors) {
            FREE(rep.colors.data);
        }
    }
    data->representations.data.clear();
}

// ### ASYNC OPERATIONS ON DATA ###
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
            // read first frame expicitly
            read_next_trajectory_frame(&data->mol_data.dynamic.trajectory);

            create_volume(data);

            /*
if (data->mol_data.dynamic.trajectory.num_frames > 0) {
    vec3 box_ext = data->mol_data.dynamic.trajectory.frame_buffer[0].box * vec3(1);
    init_volume(&data->density_volume.volume, math::max(ivec3(1), ivec3(box_ext) / VOLUME_DOWNSAMPLE_FACTOR));
    data->density_volume.model_to_world_matrix = volume::compute_model_to_world_matrix(vec3(0), box_ext);
    data->density_volume.texture_to_model_matrix = volume::compute_texture_to_model_matrix(data->density_volume.volume.dim);
}
            */

            while (read_next_trajectory_frame(&data->mol_data.dynamic.trajectory)) {
                data->async.trajectory.fraction =
                    data->mol_data.dynamic.trajectory.num_frames / (float)data->mol_data.dynamic.trajectory.frame_offsets.count;
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
    data->density_volume.world_to_texture_matrix =
        math::inverse(data->density_volume.model_to_world_matrix * data->density_volume.texture_to_model_matrix);
}
