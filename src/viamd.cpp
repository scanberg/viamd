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
#include <mol/molecule.h>
#include <mol/trajectory.h>
#include <mol/trajectory_utils.h>
#include <mol/molecule_utils.h>
#include <mol/pdb_utils.h>
#include <mol/gro_utils.h>
#include <mol/stats.h>
#include <gfx/immediate_draw_utils.h>
#include <gfx/postprocessing_utils.h>

#include <stdio.h>
#include <thread>

//#define VIAMD_RELEASE 1

#ifdef _WIN32
constexpr Key::Key_t CONSOLE_KEY = Key::KEY_GRAVE_ACCENT;
#elif __APPLE__
constexpr Key::Key_t CONSOLE_KEY = Key::KEY_WORLD_1;
#else
// @TODO: Make sure this is right for Linux?
constexpr Key::Key_t CONSOLE_KEY = Key::KEY_GRAVE_ACCENT;
#endif
constexpr unsigned int NO_PICKING_IDX = 0xffffffff;
constexpr const char* FILE_EXTENSION = "via";

constexpr uint32 DEL_BTN_COLOR = 0xff1111cc;
constexpr uint32 DEL_BTN_HOVER_COLOR = 0xff3333dd;
constexpr uint32 DEL_BTN_ACTIVE_COLOR = 0xff5555ff;
constexpr uint32 TEXT_BG_ERROR_COLOR = 0xaa222299;

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

inline ImVec4& vec_cast(vec4& v) { return *(ImVec4*)(&v); }
inline vec4& vec_cast(ImVec4& v) { return *(vec4*)(&v); }
inline ImVec2& vec_cast(vec2& v) { return *(ImVec2*)(&v); }
inline vec2& vec_cast(ImVec2& v) { return *(vec2*)(&v); }

inline ImVec2 operator+(const ImVec2& a, const ImVec2& b) { return {a.x + b.x, a.y + b.y}; }
inline ImVec2 operator-(const ImVec2& a, const ImVec2& b) { return {a.x - b.x, a.y - b.y}; }
inline ImVec2 operator*(const ImVec2& a, const ImVec2& b) { return {a.x * b.x, a.y * b.y}; }
inline ImVec2 operator/(const ImVec2& a, const ImVec2& b) { return {a.x / b.x, a.y / b.y}; }

enum PlaybackInterpolationMode { NEAREST, LINEAR, LINEAR_PERIODIC, CUBIC, CUBIC_PERIODIC };

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

struct MoleculeData {
    MoleculeDynamic dynamic{};
    DynamicArray<float> atom_radii{};
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

    struct {
        String molecule{};
        String trajectory{};
        String workspace{};
    } files;

    // --- CAMERA ---
    Camera camera;
    TrackballController controller;

    // --- MOL DATA ---
    MoleculeData mol_data;

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
        } statistics;

        struct {
            ThreadSyncData sync{};
            float fraction = 0.f;
            bool query_update = false;
        } backbone_angles;
    } async;

    struct {
        bool show_window;
        DynamicArray<Representation> data;
    } representations;

    // --- ATOM SELECTION ---
    AtomSelection hovered;
    AtomSelection selected;

    // --- STATISTICAL DATA ---
    struct {
        bool show_property_window = false;
        bool show_timeline_window = false;
        bool show_distribution_window = false;
    } statistics;

    // --- FRAMEBUFFER ---
    MainFramebuffer fbo;
    unsigned int picking_idx = NO_PICKING_IDX;

    // --- PLAYBACK ---
    float64 time = 0.f;  // needs to be double precision
    float frames_per_second = 10.f;
    bool is_playing = false;
    PlaybackInterpolationMode interpolation = PlaybackInterpolationMode::LINEAR_PERIODIC;

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
        float intensity = 1.5f;
        float radius = 6.0f;
    } ssao;

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

    // DEBUG DRAW
    struct {
        struct {
            bool enabled = false;
        } spline;

        struct {
            bool enabled = false;
        } backbone;
    } debug_draw;

    // --- CONSOLE ---
    Console console;
    bool show_console;

	bool high_res_font = false;
};

static void reset_view(ApplicationData* data, bool reposition_camera = true);
static float compute_avg_ms(float dt);
static uint32 get_picking_id(uint32 fbo, int32 x, int32 y);

static void draw_main_menu(ApplicationData* data);
static void draw_representations_window(ApplicationData* data);
static void draw_property_window(ApplicationData* data);
static void draw_timeline_window(ApplicationData* data);
static void draw_distribution_window(ApplicationData* data);
static void draw_ramachandran_window(ApplicationData* data);
static void draw_atom_info(const MoleculeStructure& mol, int atom_idx, int x, int y);
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

// Async operations
static void load_trajectory_async(ApplicationData* data);
//static void remove_properties(ApplicationData* data);
//static void clear_properties(ApplicationData* data);
//static void compute_statistics_async(ApplicationData* data);
static void compute_backbone_angles_async(ApplicationData* data);

int main(int, char**) {
    ApplicationData data;

    // Init platform
    platform::initialize(&data.ctx, 1920, 1080, "VIAMD");
    data.ctx.window.vsync = false;

    init_main_framebuffer(&data.fbo, data.ctx.framebuffer.width, data.ctx.framebuffer.height);

    // Init subsystems
    logging::initialize();
    immediate::initialize();
    draw::initialize();
    ramachandran::initialize();
    stats::initialize();
    filter::initialize();
    postprocessing::initialize(data.fbo.width, data.fbo.height);

    // Standard output
    logging::register_backend([](CString str, logging::Severity, void*) { printf("%s\n", str.cstr()); });

    // App Console output
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

    // Setup style
    ImGui::StyleColorsClassic();

    bool show_demo_window = false;
    const vec4 CLEAR_COLOR = vec4(1, 1, 1, 1);
    const vec4 CLEAR_INDEX = vec4(1, 1, 1, 1);

#ifdef VIAMD_RELEASE
    allocate_and_parse_pdb_from_string(&data.mol_data.dynamic, CAFFINE_PDB);
    data.mol_data.atom_radii = compute_atom_radii(data.mol_data.dynamic.molecule.atom_elements);
#else
    stats::create_property("b1", "distance resatom(resname(ALA), 1) com(resname(ALA))");
    load_molecule_data(&data, PROJECT_SOURCE_DIR "/data/1ALA-250ns-2500frames.pdb");
#endif
    reset_view(&data);
    create_default_representation(&data);

    // Main loop
    while (!data.ctx.window.should_close) {
        platform::update(&data.ctx);
		stats::async_update(data.mol_data.dynamic, data.time_filter.range);

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
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);

        if (data.ctx.input.key.hit[CONSOLE_KEY]) {
            data.console.visible = !data.console.visible;
        }

        if (data.async.trajectory.sync.running) {
            constexpr float TICK_INTERVAL_SEC = 3.f;
            static float time = 0.f;
            time += data.ctx.timing.dt;
            if (time > TICK_INTERVAL_SEC) {
                time = 0.f;
				stats::set_all_property_flags(true, true, true);
                compute_backbone_angles_async(&data);
            }
        }

        float ms = compute_avg_ms(data.ctx.timing.dt);
        bool time_changed = false;

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

        if (data.is_playing) {
            data.time += data.ctx.timing.dt * data.frames_per_second;
        }

        {
            static float64 prev_time = data.time;
            if (data.time != prev_time) {
                time_changed = true;
            }
            prev_time = data.time;
        }

        if (data.time_filter.dynamic_window) {
            float max_frame = data.mol_data.dynamic.trajectory ? data.mol_data.dynamic.trajectory.num_frames : 1.f;
            data.time_filter.range.x = math::max((float)data.time - data.time_filter.window_extent * 0.5f, 0.f);
            data.time_filter.range.y = math::min((float)data.time + data.time_filter.window_extent * 0.5f, max_frame);
        }

        if (data.mol_data.dynamic.trajectory && time_changed) {
            int last_frame = data.mol_data.dynamic.trajectory.num_frames - 1;
            data.time = math::clamp(data.time, 0.0, float64(last_frame));
            if (data.time == float64(last_frame)) data.is_playing = false;

            if (data.time_filter.dynamic_window) {
                stats::set_all_property_flags(false, false, true);
                //compute_statistics_async(&data);
            }

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

		// CAMERA CONTROLS
        if (!ImGui::GetIO().WantCaptureMouse) {
            data.controller.input.rotate_button = data.ctx.input.mouse.down[0];
            data.controller.input.pan_button = data.ctx.input.mouse.down[1];
            data.controller.input.dolly_button = data.ctx.input.mouse.down[2];
            data.controller.input.mouse_coord_prev = data.ctx.input.mouse.coord_prev;
            data.controller.input.mouse_coord_curr = data.ctx.input.mouse.coord_curr;
            data.controller.input.screen_size = vec2(data.ctx.window.width, data.ctx.window.height);
            data.controller.input.dolly_delta = data.ctx.input.mouse.scroll.y;
            data.controller.update();
            data.camera.position = data.controller.position;
            data.camera.orientation = data.controller.orientation;

            if (data.ctx.input.mouse.release[0] && data.ctx.input.mouse.velocity == vec2(0)) {
                data.selected = data.hovered;
            }
        }
        if (!ImGui::GetIO().WantCaptureKeyboard) {
            if (data.ctx.input.key.hit[Key::KEY_SPACE]) data.is_playing = !data.is_playing;
        }

        // RENDER TO FBO
        mat4 view_mat = compute_world_to_view_matrix(data.camera);
        mat4 proj_mat = compute_perspective_projection_matrix(data.camera, data.fbo.width, data.fbo.height);
        mat4 inv_proj_mat = math::inverse(proj_mat);

        immediate::set_view_matrix(view_mat);
        immediate::set_proj_matrix(proj_mat);

        for (const auto& rep : data.representations.data) {
            if (!rep.enabled) continue;
            switch (rep.type) {
                case Representation::VDW:
                    draw::draw_vdw(data.mol_data.dynamic.molecule.atom_positions, data.mol_data.atom_radii, rep.colors, view_mat, proj_mat,
                                   rep.radius);
                    break;
                case Representation::LICORICE:
                    draw::draw_licorice(data.mol_data.dynamic.molecule.atom_positions, data.mol_data.dynamic.molecule.bonds, rep.colors, view_mat,
                                        proj_mat, rep.radius);
                    break;
                case Representation::RIBBONS:
                    draw::draw_ribbons(data.mol_data.dynamic.molecule.backbone_segments, data.mol_data.dynamic.molecule.chains,
                                       data.mol_data.dynamic.molecule.atom_positions, rep.colors, view_mat, proj_mat, rep.num_subdivisions,
                                       rep.tension, rep.width, rep.thickness);
                    break;
            }
        }

        // PICKING
        {
            ivec2 coord = {data.ctx.input.mouse.coord_curr.x, data.ctx.framebuffer.height - data.ctx.input.mouse.coord_curr.y};
            data.picking_idx = get_picking_id(data.fbo.id, coord.x, coord.y);

            data.hovered = {};
            if (data.picking_idx != NO_PICKING_IDX) {
                data.hovered.atom_idx = data.picking_idx;
                if (-1 < data.hovered.atom_idx && data.hovered.atom_idx < data.mol_data.dynamic.molecule.atom_residue_indices.count) {
                    data.hovered.residue_idx = data.mol_data.dynamic.molecule.atom_residue_indices[data.hovered.atom_idx];
                }
                if (-1 < data.hovered.residue_idx && data.hovered.residue_idx < data.mol_data.dynamic.molecule.residues.count) {
                    data.hovered.chain_idx = data.mol_data.dynamic.molecule.residues[data.hovered.residue_idx].chain_idx;
                }
            }
        }

        // Activate backbuffer
        glDisable(GL_DEPTH_TEST);
        glDepthFunc(GL_ALWAYS);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glDrawBuffer(GL_BACK);
        glClear(GL_COLOR_BUFFER_BIT);

        // Render deferred
        postprocessing::render_deferred(data.fbo.tex_depth, data.fbo.tex_color, data.fbo.tex_normal, inv_proj_mat);

        // Apply post processing
        // postprocessing::apply_tonemapping(data.fbo.tex_color);
        if (data.ssao.enabled) {
            postprocessing::apply_ssao(data.fbo.tex_depth, data.fbo.tex_normal, proj_mat, data.ssao.intensity, data.ssao.radius);
        }

        stats::visualize(data.mol_data.dynamic);
        immediate::flush();

        // GUI ELEMENTS
        data.console.Draw("VIAMD", data.ctx.window.width, data.ctx.window.height, data.ctx.timing.dt);

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
                ivec2 pos = data.ctx.input.mouse.coord_curr;
                draw_atom_info(data.mol_data.dynamic.molecule, data.picking_idx, pos.x, pos.y);
            }
        }

        draw_async_info(&data);

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

    data.async.trajectory.sync.signal_stop();
    data.async.backbone_angles.sync.signal_stop();
    data.async.statistics.sync.signal_stop();

    data.async.trajectory.sync.wait_until_finished();
    data.async.backbone_angles.sync.wait_until_finished();
    data.async.statistics.sync.wait_until_finished();

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

    if (reposition_camera) {
        data->controller.look_at(cent, cent + size * 2.f);
        data->camera.position = data->controller.position;
        data->camera.orientation = data->controller.orientation;
    }
    data->camera.near_plane = 1.f;
    data->camera.far_plane = math::length(size) * 30.f;
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
                auto res = platform::open_file_dialog("pdb,gro,xtc");
                if (res.action == platform::FileDialogResult::FILE_OK) {
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
                auto res = platform::open_file_dialog(FILE_EXTENSION);
                if (res.action == platform::FileDialogResult::FILE_OK) {
                    load_workspace(data, res.path);
                }
            }
            if (ImGui::MenuItem("Save", "CTRL+S")) {
                if (!data->files.workspace) {
                    auto res = platform::save_file_dialog({}, FILE_EXTENSION);
                    if (res.action == platform::FileDialogResult::FILE_OK) {
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
                auto res = platform::save_file_dialog({}, FILE_EXTENSION);
                if (res.action == platform::FileDialogResult::FILE_OK) {
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

			ImGui::BeginGroup();
			if (ImGui::Checkbox("Use high-res font", &data->high_res_font)) {
				if (ImGui::GetIO().Fonts->Fonts.size() > 1 && data->high_res_font) {
					ImGui::GetIO().FontDefault = ImGui::GetIO().Fonts->Fonts[1];
					ImGui::GetIO().FontGlobalScale = 0.75f;
				}
				else {
					ImGui::GetIO().FontDefault = ImGui::GetIO().Fonts->Fonts[0];
					ImGui::GetIO().FontGlobalScale = 1.0f;
				}
			}
			ImGui::EndGroup();

            // DEBUG DRAW
            ImGui::BeginGroup();
            ImGui::Checkbox("Spline", &data->debug_draw.spline.enabled);
            ImGui::Checkbox("Backbone", &data->debug_draw.backbone.enabled);
            ImGui::EndGroup();

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
    ImGui::Spacing();

	ImGui::PushItemWidth(-1);
    ImGui::Columns(3, "columns", true);
    ImGui::Separator();

    if (first_time_shown) {
        first_time_shown = false;
        ImGui::SetColumnWidth(0, ImGui::GetWindowContentRegionWidth() * 0.1f);
        ImGui::SetColumnWidth(1, ImGui::GetWindowContentRegionWidth() * 0.7f);
        ImGui::SetColumnWidth(2, ImGui::GetWindowContentRegionWidth() * 0.2f);
    }

    ImGui::Text("name");
    ImGui::NextColumn();
    ImGui::Text("args");
    ImGui::NextColumn();
    ImGui::NextColumn();

    //bool compute_stats = false;

    for (int i = 0; i < stats::get_property_count(); i++) {
        auto prop = stats::get_property(i);

        ImGui::Separator();
        ImGui::PushID(i);

        if (!prop->valid) ImGui::PushStyleColor(ImGuiCol_FrameBg, TEXT_BG_ERROR_COLOR);
        ImGui::PushItemWidth(-1);
        if (ImGui::InputText("##name", prop->name.buffer, prop->name.MAX_LENGTH, ImGuiInputTextFlags_EnterReturnsTrue)) {
            prop->data_dirty = true;
            //compute_stats = true;
        }
        ImGui::PopItemWidth();
        ImGui::NextColumn();
        ImGui::PushItemWidth(-1);
        if (ImGui::InputText("##args", prop->args.buffer, prop->args.MAX_LENGTH, ImGuiInputTextFlags_EnterReturnsTrue)) {
            prop->data_dirty = true;
            //compute_stats = true;
        }
        ImGui::PopItemWidth();
        if (!prop->valid) {
            ImGui::PopStyleColor();
            if (!prop->valid && prop->error_msg && ImGui::GetHoveredID() == ImGui::GetID("##args")) {
                ImGui::SetTooltip("%s", prop->error_msg.cstr());
            }
        }
        ImGui::NextColumn();
        if (ImGui::ArrowButton("up", ImGuiDir_Up)) {
            stats::move_property_up(prop);
            //compute_stats = true;
        }
        ImGui::SameLine();
        if (ImGui::ArrowButton("down", ImGuiDir_Down)) {
            stats::move_property_down(prop);
            //compute_stats = true;
        }
        ImGui::SameLine();
        ImGui::PushStyleColor(ImGuiCol_Button, DEL_BTN_COLOR);
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, DEL_BTN_HOVER_COLOR);
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, DEL_BTN_ACTIVE_COLOR);
        if (ImGui::Button("remove")) {
            stats::remove_property(prop);
        }
        ImGui::SameLine();
        ImGui::Checkbox("show", &prop->visualize);
        ImGui::PopStyleColor(3);
        ImGui::NextColumn();
        ImGui::PopID();
    }
    ImGui::Columns(1);
    ImGui::Separator();
    ImGui::PopID();
	ImGui::PopItemWidth();
    ImGui::End();

    //if (compute_stats) {
        // stats::compute_stats(data->mol_data.dynamic);
        //compute_statistics_async(data);
    //}
}

static void draw_atom_info(const MoleculeStructure& mol, int atom_idx, int x, int y) {

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
    float stats_fract = data->async.statistics.fraction;

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
                //compute_statistics_async(data);
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
    if (ImGui::Begin("Timelines", &data->statistics.show_timeline_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
        static float zoom = 1.f;
        ImGui::BeginChild("Scroll Region", ImVec2(0, 0), true, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_HorizontalScrollbar);

        const int max_frame = data->mol_data.dynamic.trajectory.num_frames;
        const ImVec2 frame_range(0, (float)max_frame);
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

        for (int i = 0; i < stats::get_property_count(); i++) {
            auto prop = stats::get_property(i);
            Array<float> prop_data = prop->data;
            CString prop_name = prop->name;
            Range prop_range = prop->data_range;
            if (!prop_data) continue;
            float pad = math::max((prop_range.y - prop_range.x) * 0.1f, 1.f);
            vec2 display_range = prop_range + vec2(-pad, pad);
            float val = (float)data->time;

            ImGui::PushID(i);
            if (ImGui::BeginPlot(prop_name, ImVec2(0, 100), frame_range, ImVec2(display_range.x, display_range.y), &val,
                                 &vec_cast(data->time_filter.range), ImGui::LinePlotFlags_AxisX | ImGui::LinePlotFlags_ShowXVal)) {
                data->time = val;
            }
            ImGui::PlotValues("Najs", prop_data.data, (int)prop_data.count);
            ImGui::EndPlot();
            ImGui::PopID();
        }
        ImGui::PopItemWidth();

        if (data->time_filter.range != old_range) {
            stats::set_all_property_flags(false, false, true);
            //compute_statistics_async(data);
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
    ImGui::Begin("Distributions", &data->statistics.show_distribution_window, ImGuiWindowFlags_NoFocusOnAppearing);
    ImGuiWindow* window = ImGui::GetCurrentWindow();
    const ImGuiStyle& style = ImGui::GetStyle();
    ImGui::PushItemWidth(-1);
    ImVec2 frame_size{ImGui::CalcItemWidth(), 100.f};

	constexpr uint32 FULL_BAR_COLOR = 0x77cc9e66;
    constexpr uint32 FULL_TXT_COLOR = 0xffcc9e66;
    constexpr uint32 FILT_BAR_COLOR = 0x5533ffff;
    constexpr uint32 FILT_TXT_COLOR = 0xff33ffff;

    for (int i = 0; i < stats::get_property_count(); i++) {
        stats::Property* prop = stats::get_property(i);
        ImGui::PushID(i);

        const ImRect frame_bb(window->DC.CursorPos, window->DC.CursorPos + ImVec2(frame_size.x, frame_size.y));
        const ImRect inner_bb(frame_bb.Min + style.FramePadding, frame_bb.Max - style.FramePadding);
        const ImRect total_bb(frame_bb.Min, frame_bb.Max);
        ImGui::ItemSize(total_bb, style.FramePadding.y);
        if (ImGui::ItemAdd(total_bb, NULL)) {
            ImGui::RenderFrame(frame_bb.Min, frame_bb.Max, ImGui::GetColorU32(ImGuiCol_FrameBg), true, style.FrameRounding);

            const float max_val = math::max(prop->filt_histogram.bin_range.y, prop->filt_histogram.bin_range.y);
            ImGui::DrawHistogram(inner_bb.Min, inner_bb.Max, prop->full_histogram.bins.data, (int32)prop->full_histogram.bins.count, max_val, FULL_BAR_COLOR);
            ImGui::DrawHistogram(inner_bb.Min, inner_bb.Max, prop->filt_histogram.bins.data, (int32)prop->filt_histogram.bins.count, max_val, FILT_BAR_COLOR);

            if (ImGui::IsItemHovered()) {
                window->DrawList->AddLine(ImVec2(ImGui::GetIO().MousePos.x, inner_bb.Min.y), ImVec2(ImGui::GetIO().MousePos.x, inner_bb.Max.y), 0xffffffff);
                float t = (ImGui::GetIO().MousePos.x - inner_bb.Min.x) / (inner_bb.Max.x - inner_bb.Min.x);
                int32 count = (int32)prop->full_histogram.bins.count;
                int32 idx = ImClamp((int32)(t * (count - 1)), 0, count - 1);
                float full_val = prop->full_histogram.bins.data[idx];
                float filt_val = prop->filt_histogram.bins.data[idx];
                ImVec2 val_range = vec_cast(prop->filt_histogram.value_range);
                ImGui::BeginTooltip();
                ImGui::Text("%.3f:", ImLerp(val_range.x, val_range.y, t));
                ImGui::TextColored(ImColor(FULL_TXT_COLOR), "%g", full_val * 100.f);
                ImGui::TextColored(ImColor(FILT_TXT_COLOR), "%g", filt_val * 100.f);
                ImGui::EndTooltip();
            }

            if (ImGui::RangeSliderFloat("##filter", &prop->filter.x, &prop->filter.y, prop->data_range.x, prop->data_range.y)) {
                prop->filt_hist_dirty = true;
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
        data->async.statistics.sync.signal_stop_and_wait();
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
                init_backbone_angles_trajectory(&data->ramachandran.backbone_angles, data->mol_data.dynamic);
                compute_backbone_angles_trajectory(&data->ramachandran.backbone_angles, data->mol_data.dynamic);
                // stats::update(data->mol_data.dynamic, );
                //compute_statistics_async(data);
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
                    /*
read_trajectory_async(&data->mol_data.dynamic.trajectory, [data]() {
    stats::compute_stats(data->mol_data.dynamic);
    compute_backbone_angles_trajectory(&data->ramachandran.backbone_angles, data->mol_data.dynamic);
});
                    */
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
                    data->camera.position = pos;
                    data->controller.position = pos;
                }
                if (compare_n(line, "Rotation=", 9)) {
                    quat rot = quat(to_vec4(trim(line.substr(9))));
                    data->camera.orientation = rot;
                    data->controller.orientation = rot;
                }
                if (compare_n(line, "Distance=", 9)) data->controller.distance = to_float(trim(line.substr(9)));
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
        if (rep.color_mapping == ColorMapping::STATIC_COLOR) fprintf(fptr, "StaticColor=%i\n", rep.enabled);
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
    for (int i = 0; i < stats::get_property_count(); i++) {
        auto prop = stats::get_property(i);
        fprintf(fptr, "[Property]\n");
        fprintf(fptr, "Name=%s\n", prop->name.cstr());
        fprintf(fptr, "Args=%s\n", prop->args.cstr());
        fprintf(fptr, "\n");
    }

    fprintf(fptr, "[RenderSettings]\n");
    fprintf(fptr, "SsaoEnabled=%i\n", data->ssao.enabled ? 1 : 0);
    fprintf(fptr, "SsaoIntensity=%g\n", data->ssao.intensity);
    fprintf(fptr, "SsaoRadius=%g\n", data->ssao.radius);
    fprintf(fptr, "\n");

    fprintf(fptr, "[Camera]\n");
    fprintf(fptr, "Position=%g,%g,%g\n", data->camera.position.x, data->camera.position.y, data->camera.position.z);
    fprintf(fptr, "Rotation=%g,%g,%g,%g\n", data->camera.orientation.x, data->camera.orientation.y, data->camera.orientation.z,
            data->camera.orientation.w);
    fprintf(fptr, "Distance=%g\n", data->controller.distance);
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
            while (read_next_trajectory_frame(&data->mol_data.dynamic.trajectory)) {
                data->async.trajectory.fraction =
                    data->mol_data.dynamic.trajectory.num_frames / (float)data->mol_data.dynamic.trajectory.frame_offsets.count;
                if (data->async.trajectory.sync.stop_signal) break;
            }
            data->async.trajectory.sync.running = false;
            data->async.trajectory.sync.stop_signal = false;

            //compute_statistics_async(data);
			stats::set_all_property_flags(true, true, true);
            compute_backbone_angles_async(data);
        });
        data->async.trajectory.sync.thread.detach();
    }
}

//static void remove_properties(ApplicationData* data) {
    //data->async.statistics.sync.signal_stop_and_wait();
    //stats::remove_all_properties();
//}

//static void clear_properties(ApplicationData* data) {
    //data->async.statistics.sync.signal_stop_and_wait();
    //stats::clear_all_properties();
//}

/*
static void compute_statistics_async(ApplicationData* data) {
    ASSERT(data);
    data->async.statistics.query_update = true;
    if (data->async.statistics.sync.running == false) {
        data->async.statistics.sync.running = true;

        data->async.statistics.sync.thread = std::thread([data]() {
            data->async.statistics.fraction = 0.0f;
            while (data->async.statistics.query_update) {
                data->async.statistics.query_update = false;
                data->async.statistics.fraction = 0.5f;
                if (data->time_filter.enabled) {
                    stats::update(data->mol_data.dynamic, &data->time_filter.range);
                } else {
                    stats::update(data->mol_data.dynamic);
                }
                if (data->async.statistics.sync.stop_signal) break;
            }
            data->async.statistics.fraction = 1.f;
            data->async.statistics.sync.running = false;
            data->async.statistics.sync.stop_signal = false;
        });
        data->async.statistics.sync.thread.detach();
    }
}
*/

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
