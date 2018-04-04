#include <imgui.h>
#include <core/platform.h>
#include <core/gl.h>
#include <core/types.h>
#include <core/hash.h>
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
#include <stats/stats.h>
#include <gfx/immediate_draw_utils.h>
#include <gfx/postprocessing_utils.h>

#include <stdio.h>

#ifdef _WIN32
constexpr Key::Key_t CONSOLE_KEY = Key::KEY_GRAVE_ACCENT;
#elif __APPLE__
constexpr Key::Key_t CONSOLE_KEY = Key::KEY_WORLD_1;
#else
// @TODO: Make sure this is right for Linux?
constexpr Key::Key_t CONSOLE_KEY = Key::KEY_GRAVE_ACCENT;
#endif
constexpr unsigned int NO_PICKING_IDX = 0xffffffff;

constexpr const char* testosterone_pdb = R"(
ATOM      1  C           1      26.334  10.561 -22.799  1.00  0.00
ATOM      2  C           1      26.570  10.310 -21.300  1.00  0.00
ATOM      3  C           1      25.355  10.680 -20.409  1.00  0.00
ATOM      4  C           1      24.998  12.157 -20.713  1.00  0.00
ATOM      5  C           1      24.767  12.453 -22.204  1.00  0.00
ATOM      6  C           1      25.970  12.029 -23.053  1.00  0.00
ATOM      7  C           1      25.204  10.780 -16.421  1.00  0.00
ATOM      8  C           1      25.656   9.334 -16.165  1.00  0.00
ATOM      9  C           1      26.790   9.024 -17.158  1.00  0.00
ATOM     10  C           1      26.396   9.103 -14.833  1.00  0.00
ATOM     11  H           1      25.945   8.688 -12.826  1.00  0.00
ATOM     12  C           1      24.749  10.968 -17.880  1.00  0.00
ATOM     13  H           1      26.643  11.188 -18.796  1.00  0.00
ATOM     14  C           1      24.457   8.379 -16.309  1.00  0.00
ATOM     15  H           1      27.554   9.796 -17.059  1.00  0.00
ATOM     17  C           1      24.153   9.767 -20.712  1.00  0.00
ATOM     18  H           1      26.415   9.800 -23.575  1.00  0.00
ATOM     19  C           1      26.333   9.093 -18.621  1.00  0.00
ATOM     20  C           1      25.795  10.514 -18.923  1.00  0.00
ATOM     21  H           1      24.087  12.408 -20.170  1.00  0.00
ATOM     22  H           1      25.816  12.786 -20.362  1.00  0.00
ATOM     23  C           1      27.071   8.880 -21.042  1.00  0.00
ATOM     24  C           1      27.477   8.712 -19.574  1.00  0.00
ATOM     25  H           1      23.886  11.907 -22.541  1.00  0.00
ATOM     26  H           1      24.602  13.522 -22.331  1.00  0.00
ATOM     27  O           1      26.607  12.821 -23.890  1.00  0.00
ATOM     29  H           1      24.374  11.016 -15.755  1.00  0.00
ATOM     30  H           1      26.036  11.454 -16.217  1.00  0.00
ATOM     31  C           1      27.397   7.722 -16.623  1.00  0.00
ATOM     32  C           1      27.273   7.862 -15.092  1.00  0.00
ATOM     33  O           1      25.421   8.842 -13.723  1.00  0.00
ATOM     34  H           1      27.012   9.966 -14.580  1.00  0.00
ATOM     35  H           1      23.838  10.389 -18.033  1.00  0.00
ATOM     36  H           1      24.539  12.025 -18.041  1.00  0.00
ATOM     37  H           1      24.808   7.348 -16.280  1.00  0.00
ATOM     38  H           1      23.758   8.548 -15.489  1.00  0.00
ATOM     39  H           1      23.955   8.566 -17.258  1.00  0.00
ATOM     40  H           1      24.432   8.728 -20.536  1.00  0.00
ATOM     41  H           1      23.320  10.035 -20.061  1.00  0.00
ATOM     42  H           1      23.855   9.892 -21.753  1.00  0.00
ATOM     43  H           1      25.521   8.379 -18.764  1.00  0.00
ATOM     44  H           1      26.275   8.174 -21.277  1.00  0.00
ATOM     45  H           1      27.933   8.682 -21.679  1.00  0.00
ATOM     46  H           1      27.751   7.672 -19.399  1.00  0.00
ATOM     47  H           1      28.335   9.353 -19.369  1.00  0.00
ATOM     48  H           1      26.838   6.857 -16.981  1.00  0.00
ATOM     49  H           1      28.441   7.629 -16.922  1.00  0.00
ATOM     50  H           1      26.809   6.965 -14.682  1.00  0.00
ATOM     51  H           1      28.265   7.989 -14.657  1.00  0.00
)";

inline ImVec4 vec_cast(vec4 v) { return ImVec4(v.x, v.y, v.z, v.w); }
inline vec4 vec_cast(ImVec4 v) { return vec4(v.x, v.y, v.z, v.w); }
inline ImVec2 vec_cast(vec2 v) { return ImVec2(v.x, v.y); }
inline vec2 vec_cast(ImVec2 v) { return vec2(v.x, v.y); }

enum struct PlaybackInterpolationMode { NEAREST_FRAME, LINEAR, LINEAR_PERIODIC, CUBIC, CUBIC_PERIODIC };

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
    Array<uint32> colors;

    // @TODO Fill in options heres
    vec4 static_color = vec4(1);
    float radii_scale = 1.f;
    int num_subdivisions = 8;
};

struct AtomSelection {
    int32 atom_idx = -1;
    int32 residue_idx = -1;
    int32 chain_idx = -1;
};

struct MoleculeData {
    MoleculeDynamic dynamic{};
    DynamicArray<float> atom_radii{};

    BackboneAnglesTrajectory backbone_angles{};
    Array<BackboneAngles> current_backbone_angles{};
};

struct ApplicationData {
    // --- PLATFORM ---
    platform::Context ctx;

    // --- CAMERA ---
    Camera camera;
    TrackballController controller;

    // --- MOL DATA ---
    MoleculeData mol_data;

    struct {
        bool show_window;
        DynamicArray<Representation> data;
    } representations;

    // --- ATOM SELECTION ---
    AtomSelection hovered;
    AtomSelection selected;

    // --- STATISTICAL DATA ---
    // stats::StatisticsContext stats;

    // Framebuffer
    MainFramebuffer fbo;
    unsigned int picking_idx = NO_PICKING_IDX;

    // --- PLAYBACK ---
    float64 time = 0.f;  // needs to be double precision
    float frames_per_second = 10.f;
    bool is_playing = false;
    PlaybackInterpolationMode interpolation = PlaybackInterpolationMode::LINEAR_PERIODIC;

    // --- VISUALS ---
    // SSAO
    struct {
        bool enabled = false;
        float intensity = 1.5f;
        float radius = 6.0f;
    } ssao;

    struct {
        bool enabled = false;
        float radius = 1.f;
        float opacity = 1.f;
        int frame_range_min = 0;
        int frame_range_max = 0;
    } ramachandran;

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
};

static void draw_main_menu(ApplicationData* data);
static void draw_representations_window(ApplicationData* data);
static void draw_atom_info(const MoleculeStructure& mol, int atom_idx, int x, int y);
static void draw_statistics(ApplicationData* data);
static void draw_ramachandran(ApplicationData* data);

static void init_main_framebuffer(MainFramebuffer* fbo, int width, int height);
static void destroy_main_framebuffer(MainFramebuffer* fbo);

static void reset_view(ApplicationData* data);
static float compute_avg_ms(float dt);
static uint32 get_picking_id(uint32 fbo, int32 x, int32 y);

static void load_molecule_data(MoleculeData* mol_data, CString file);

static void create_default_representation(ApplicationData* data);
static void remove_representation(ApplicationData* data, int idx);
static void reset_representations(ApplicationData* data);

int main(int, char**) {
    ApplicationData data;

    auto dir_list = platform::list_directory(".");
    for (const auto& d : dir_list) {
        printf("%s\n", d.name.beg());
    }
    printf("%s\n", platform::get_cwd().beg());

    // Init platform
    platform::initialize(&data.ctx, 1920, 1080, "VIAMD");
    data.ctx.window.vsync = false;

    init_main_framebuffer(&data.fbo, data.ctx.framebuffer.width, data.ctx.framebuffer.height);

    // Init subsystems
    immediate::initialize();
    draw::initialize();
	ramachandran::initialize();
    stats::initialize();
    filter::initialize();
    postprocessing::initialize(data.fbo.width, data.fbo.height);

    // Setup style
    ImGui::StyleColorsClassic();

    bool show_demo_window = false;
    vec4 clear_color = vec4(1, 1, 1, 1);
    vec4 clear_index = vec4(1, 1, 1, 1);

    data.mol_data.dynamic = allocate_and_parse_pdb_from_string(testosterone_pdb);
    data.mol_data.atom_radii = compute_atom_radii(data.mol_data.dynamic.molecule->atom_elements);
    create_default_representation(&data);
    reset_view(&data);

    // load_molecule_data(&data.mol_data, PROJECT_SOURCE_DIR "/data/1ALA-250ns-2500frames.pdb");

    // data.dynamic.molecule = allocate_and_load_gro_from_file(PROJECT_SOURCE_DIR "/data/bta-gro/20-mol-p.gro");
    // data.dynamic.molecule= allocate_and_load_gro_from_file(PROJECT_SOURCE_DIR "/data/peptides/box_2.gro");
    // data.dynamic.molecule = allocate_and_load_gro_from_file(PROJECT_SOURCE_DIR "/data/shaoqi/md-nowater.gro");
    // auto res = load_gro_from_file(PROJECT_SOURCE_DIR "/data/peptides/box_2.gro");
    // data.dynamic.molecule = allocate_and_load_gro_from_file(PROJECT_SOURCE_DIR "/data/amyloid/centered.gro");
    // data.dynamic.molecule = allocate_and_load_gro_from_file(PROJECT_SOURCE_DIR "/data/water/water.gro");
    // data.dynamic.molecule = allocate_and_load_gro_from_file(PROJECT_SOURCE_DIR "/data/amyloid-6T/conf-60-6T.gro");
    // auto res = load_gro_from_file(PROJECT_SOURCE_DIR "/data/yuya/nowat_npt.gro");
    // data.dynamic = allocate_and_load_pdb_from_file(PROJECT_SOURCE_DIR "/data/5ulj.pdb");
    // data.dynamic = allocate_and_load_pdb_from_file(PROJECT_SOURCE_DIR "/data/1ALA-250ns-2500frames.pdb");

    // data.dynamic.trajectory = allocate_trajectory(PROJECT_SOURCE_DIR "/data/bta-gro/traj-centered.xtc");
    // data.dynamic.trajectory = allocate_trajectory(PROJECT_SOURCE_DIR "/data/peptides/md_0_1_noPBC_2.xtc");
    // data.dynamic.trajectory = allocate_trajectory(PROJECT_SOURCE_DIR "/data/shaoqi/md-centered.xtc");
    // Trajectory* traj = allocate_trajectory(PROJECT_SOURCE_DIR "/data/peptides/md_0_1_noPBC_2.xtc");
    // data.dynamic.trajectory = allocate_trajectory(PROJECT_SOURCE_DIR "/data/amyloid/centered.xtc");
    // data.dynamic.trajectory = allocate_trajectory(PROJECT_SOURCE_DIR "/data/amyloid-6T/prod-centered.xtc");
    // Trajectory* traj = allocate_trajectory(PROJECT_SOURCE_DIR "/data/yuya/traj-centered.xtc");
    // if (data.dynamic.trajectory)
    //	read_trajectory_async(data.dynamic.trajectory);

    auto g1 = stats::create_group("group1", "resid ALA");
    stats::create_property(g1, "b1", "dist 1 2");
    stats::create_property(g1, "a1", "angle 1 2 3");
    stats::create_property(g1, "d1", "dist 1 2 3 4");

    // stats::compute_stats(&data.dynamic);

    /*
    if (data.mol_data.dynamic.molecule->chains.count > 0) {
            backbone = compute_backbone_segments(get_residues(*data.mol_data.dynamic.molecule, get_chain(*data.mol_data.dynamic.molecule, 0)),
    data.mol_data.dynamic.molecule->atom_labels); if (data.mol_data.dynamic.trajectory) backbone_angles =
    compute_backbone_angles_trajectory(*data.mol_data.dynamic.trajectory, data.mol_data.dynamic.molecule->backbone_segments); current_backbone_angles
    = compute_backbone_angles(data.mol_data.dynamic.molecule->atom_positions, backbone); current_spline =
    compute_spline(data.mol_data.dynamic.molecule->atom_positions, backbone, 8);
    }
    */

    // Main loop
    while (!data.ctx.window.should_close) {
        platform::update(&data.ctx);

        // RESIZE FRAMEBUFFER?
        if (data.fbo.width != data.ctx.framebuffer.width || data.fbo.height != data.ctx.framebuffer.height) {
            init_main_framebuffer(&data.fbo, data.ctx.framebuffer.width, data.ctx.framebuffer.height);
            postprocessing::initialize(data.fbo.width, data.fbo.height);
        }

        // Setup fbo and clear textures
        glViewport(0, 0, data.fbo.width, data.fbo.height);

        const GLenum draw_buffers[] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2};
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, data.fbo.id);

        // Clear color, normal and depth buffer
        glDrawBuffers(2, draw_buffers);
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Clear picking buffer
        glDrawBuffer(GL_COLOR_ATTACHMENT2);
        glClearColor(clear_index.x, clear_index.y, clear_index.z, clear_index.w);
        glClear(GL_COLOR_BUFFER_BIT);

        // Enable all draw buffers
        glDrawBuffers(3, draw_buffers);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);

        if (data.ctx.input.key.hit[CONSOLE_KEY]) {
            data.console.visible = !data.console.visible;
        }

        float ms = compute_avg_ms(data.ctx.timing.dt);
        bool time_changed = false;

        ImGui::Begin("Misc");
        ImGui::Text("%.2f ms (%.1f fps)", ms, 1000.f / (ms));
        ImGui::Text("MouseVel: %g, %g", data.ctx.input.mouse.velocity.x, data.ctx.input.mouse.velocity.y);
        ImGui::Text("Camera Pos: %g, %g, %g", data.camera.position.x, data.camera.position.y, data.camera.position.z);
        ImGui::Checkbox("Show Demo Window", &show_demo_window);
        if (ImGui::Button("Reset View")) {
            reset_view(&data);
        }
        if (data.mol_data.dynamic.trajectory) {
            ImGui::Text("Num Frames: %i", data.mol_data.dynamic.trajectory->num_frames);
            float t = (float)data.time;
            if (ImGui::SliderFloat("Time", &t, 0, (float)(data.mol_data.dynamic.trajectory->num_frames - 1))) {
                time_changed = true;
                data.time = t;
            }
            ImGui::SliderFloat("Frames Per Second", &data.frames_per_second, 0.1f, 100.f);
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
        }
        ImGui::End();

        if (data.is_playing) {
            data.time += data.ctx.timing.dt * data.frames_per_second;
            time_changed = true;
        }

        if (data.mol_data.dynamic.trajectory && time_changed) {
            data.time = math::clamp(data.time, 0.0, float64(data.mol_data.dynamic.trajectory->num_frames - 1));
            if (data.time == float64(data.mol_data.dynamic.trajectory->num_frames - 1)) data.is_playing = false;

            int prev_frame_idx = math::max(0, (int)data.time);
            int next_frame_idx = math::min(prev_frame_idx + 1, data.mol_data.dynamic.trajectory->num_frames - 1);
            if (prev_frame_idx == next_frame_idx) {
                copy_trajectory_positions(data.mol_data.dynamic.molecule->atom_positions, *data.mol_data.dynamic.trajectory, prev_frame_idx);
            } else {
                // INTERPOLATE
                auto prev_frame = get_trajectory_frame(*data.mol_data.dynamic.trajectory, prev_frame_idx);
                auto next_frame = get_trajectory_frame(*data.mol_data.dynamic.trajectory, next_frame_idx);

                float t = (float)math::fract(data.time);
                linear_interpolation_periodic(data.mol_data.dynamic.molecule->atom_positions, prev_frame.atom_positions, next_frame.atom_positions, t,
                                              prev_frame.box);
            }

            for (auto& r : data.representations.data) {
                if (r.type == Representation::RIBBONS) {
                }
            }
            // data.mol_data.current_spline = compute_spline(data.mol_data.dynamic.molecule->atom_positions, backbone, 8);
        }

        if (!ImGui::GetIO().WantCaptureMouse) {
            data.controller.input.rotate_button = data.ctx.input.mouse.down[0];
            data.controller.input.pan_button = data.ctx.input.mouse.down[1];
            data.controller.input.dolly_button = data.ctx.input.mouse.down[2];
            data.controller.input.mouse_coord_prev = data.ctx.input.mouse.coord_prev;
            data.controller.input.mouse_coord_curr = data.ctx.input.mouse.coord_curr;
            data.controller.input.screen_size = vec2(data.ctx.window.width, data.ctx.window.height);

            // data.controller.input.prev_mouse_ndc = data.ctx.input.mouse.ndc_prev;
            // data.controller.input.curr_mouse_ndc = data.ctx.input.mouse.ndc_curr;
            data.controller.input.dolly_delta = data.ctx.input.mouse.scroll.y;
            data.controller.update();
            data.camera.position = data.controller.position;
            data.camera.orientation = data.controller.orientation;

            if (data.ctx.input.mouse.release[0] && data.ctx.input.mouse.velocity == vec2(0)) {
                data.selected = data.hovered;
            }
        }

        // RENDER TO FBO
        mat4 view_mat = compute_world_to_view_matrix(data.camera);
        mat4 proj_mat = compute_perspective_projection_matrix(data.camera, data.fbo.width, data.fbo.height);
        mat4 inv_proj_mat = math::inverse(proj_mat);

        for (const auto& rep : data.representations.data) {
            switch (rep.type) {
                case Representation::VDW:
                    draw::draw_vdw(data.mol_data.dynamic.molecule->atom_positions, data.mol_data.atom_radii, rep.colors, view_mat, proj_mat,
                                   rep.radii_scale);
                    break;
                case Representation::LICORICE:
                    draw::draw_licorice(data.mol_data.dynamic.molecule->atom_positions, data.mol_data.dynamic.molecule->bonds, rep.colors, view_mat,
                                        proj_mat, rep.radii_scale);
                    break;
                case Representation::RIBBONS:
                    draw::draw_ribbons(data.mol_data.dynamic.molecule->backbone_segments, data.mol_data.dynamic.molecule->chains,
                                       data.mol_data.dynamic.molecule->atom_positions, rep.colors, view_mat, proj_mat, rep.num_subdivisions);
                    break;
            }
        }
        // draw::draw_vdw(data.dynamic.molecule->atom_positions, data.atom_radii, data.atom_colors, view_mat, proj_mat, radii_scale);
        // draw::draw_licorice(data.mol_struct->atom_positions, data.mol_struct->bonds, data.atom_colors, view_mat, proj_mat, radii_scale);
        // draw::draw_ribbons(current_spline, view_mat, proj_mat);

        // PICKING
        {
            ivec2 coord = {data.ctx.input.mouse.coord_curr.x, data.ctx.framebuffer.height - data.ctx.input.mouse.coord_curr.y};
            data.picking_idx = get_picking_id(data.fbo.id, coord.x, coord.y);

            data.hovered = {};
            if (data.picking_idx != NO_PICKING_IDX) {
                data.hovered.atom_idx = data.picking_idx;
                if (-1 < data.hovered.atom_idx && data.hovered.atom_idx < data.mol_data.dynamic.molecule->atom_residue_indices.count) {
                    data.hovered.residue_idx = data.mol_data.dynamic.molecule->atom_residue_indices[data.hovered.atom_idx];
                }
                if (-1 < data.hovered.residue_idx && data.hovered.residue_idx < data.mol_data.dynamic.molecule->residues.count) {
                    data.hovered.chain_idx = data.mol_data.dynamic.molecule->residues[data.hovered.residue_idx].chain_idx;
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

        // Apply tone mapping
        // postprocessing::apply_tonemapping(data.fbo.tex_color);
        if (data.ssao.enabled) {
            postprocessing::apply_ssao(data.fbo.tex_depth, data.fbo.tex_normal, proj_mat, data.ssao.intensity, data.ssao.radius);
        }

        /*
if (data.debug_draw.backbone.enabled) {
    draw::draw_backbone(backbone, data.mol_data.dynamic.molecule->atom_positions, view_mat, proj_mat);
}
if (data.debug_draw.spline.enabled) {
    draw::draw_spline(current_spline, view_mat, proj_mat);
}
        */

        // GUI ELEMENTS
        if (data.ramachandran.enabled) {
			if (data.mol_data.dynamic.trajectory && data.mol_data.dynamic.trajectory->is_loading) {
				static int32 prev_frame = 0;
				if (get_backbone_angles_trajectory_current_frame_count(data.mol_data.backbone_angles) - prev_frame > 100) {
					compute_backbone_angles_trajectory(&data.mol_data.backbone_angles, data.mol_data.dynamic);
					prev_frame = data.mol_data.backbone_angles.num_frames;
					data.ramachandran.frame_range_max = data.mol_data.backbone_angles.num_frames;
				}
			}
            draw_ramachandran(&data);
        }

        data.console.Draw("VIAMD", data.ctx.window.width, data.ctx.window.height, data.ctx.timing.dt);

        draw_main_menu(&data);
        if (data.representations.show_window) {
            draw_representations_window(&data);
        }

        if (!ImGui::GetIO().WantCaptureMouse) {
            if (data.picking_idx != NO_PICKING_IDX) {
                ivec2 pos = data.ctx.input.mouse.coord_curr;
                draw_atom_info(*data.mol_data.dynamic.molecule, data.picking_idx, pos.x, pos.y);
            }
        }

        // draw_statistics(&data);

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

    destroy_main_framebuffer(&data.fbo);

    platform::shutdown(&data.ctx);

    return 0;
}

void draw_random_triangles(const mat4& mvp) {
    immediate::set_view_matrix(mvp);
    immediate::set_proj_matrix(mat4(1));
    math::set_rnd_seed(0);
    for (int i = 0; i < 500; i++) {
        vec3 v0 = vec3(math::rnd(), math::rnd(), math::rnd()) * 50.f - 50.f;
        vec3 v1 = vec3(math::rnd(), math::rnd(), math::rnd()) * 50.f - 50.f;
        vec3 v2 = vec3(math::rnd(), math::rnd(), math::rnd()) * 50.f - 50.f;
        immediate::draw_triangle(&v0[0], &v1[0], &v2[0], immediate::COLOR_RED);
    }
    immediate::flush();
}

// @NOTE: Perhaps this can be done with a simple running mean?
static float compute_avg_ms(float dt) {
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

    /*
    constexpr int num_frames = 100;
    static float ms_buffer[num_frames] = {};
    static int next_idx = 0;
    next_idx = (next_idx + 1) % num_frames;
    ms_buffer[next_idx] = dt * 1000.f; // seconds to milliseconds

    float avg = 0.f;
    for (int i = 0; i < num_frames; i++) {
            avg += ms_buffer[i];
    }
    return avg / (float)num_frames;
    */
}

uint32 get_picking_id(uint32 fbo_id, int32 x, int32 y) {
    unsigned char color[4];
    glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo_id);
    glReadBuffer(GL_COLOR_ATTACHMENT2);
    glReadPixels(x, y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, color);
    glReadBuffer(GL_NONE);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    return color[0] + (color[1] << 8) + (color[2] << 16) + (color[3] << 24);
}

static void draw_main_menu(ApplicationData* data) {
    ASSERT(data);
    bool new_clicked = false;
    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("New", "CTRL+N")) new_clicked = true;
            if (ImGui::MenuItem("Load Data", "CTRL+L")) {
                platform::Path path = platform::open_file_dialog("pdb,gro,xtc");
                load_molecule_data(&data->mol_data, path);
                if (data->representations.data.count > 0) {
                    reset_representations(data);
                } else {
                    create_default_representation(data);
                }
                reset_view(data);
            }
            if (ImGui::MenuItem("Open", "CTRL+O")) {
                platform::open_file_dialog();
            }
            if (ImGui::MenuItem("Save", "CTRL+S")) {
            }
            if (ImGui::MenuItem("Save As")) {
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

            // DEBUG DRAW
            ImGui::BeginGroup();
            ImGui::Checkbox("Spline", &data->debug_draw.spline.enabled);
            ImGui::Checkbox("Backbone", &data->debug_draw.backbone.enabled);
            ImGui::EndGroup();

            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Windows")) {
            ImGui::Checkbox("Representations", &data->representations.show_window);
            ImGui::Checkbox("Ramachandran", &data->ramachandran.enabled);

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

    if (ImGui::Button("Create New Representation")) {
        create_default_representation(data);
    }

    for (int i = 0; i < data->representations.data.count; i++) {
        auto& rep = data->representations.data[i];
        ImGui::Separator();
        ImGui::BeginGroup();

        bool recompute_colors = false;
        ImGui::PushID(i);
        ImGui::InputText("name", rep.name.buffer, rep.name.MAX_LENGTH);
        ImGui::SameLine();
        if (ImGui::Button("remove")) {
            remove_representation(data, i);
        }
        if (ImGui::InputText("filter", rep.filter.buffer, rep.filter.MAX_LENGTH, ImGuiInputTextFlags_EnterReturnsTrue)) {
            recompute_colors = true;
        }
        ImGui::Combo("type", (int*)(&rep.type), "VDW\0Licorice\0Ribbons\0\0");
        if (ImGui::Combo("color mapping", (int*)(&rep.color_mapping), "Static Color\0CPK\0Res Id\0Res Idx\0Chain Id\0Chain Idx\0\0")) {
            recompute_colors = true;
        }
        if (rep.type == Representation::VDW || rep.type == Representation::LICORICE) {
            ImGui::SliderFloat("radii scale", &rep.radii_scale, 0.1f, 2.f);
        }
        if (rep.type == Representation::RIBBONS) {
            ImGui::SliderInt("spline subdivisions", &rep.num_subdivisions, 1, 16);
        }
        if (rep.color_mapping == ColorMapping::STATIC_COLOR) {
            if (ImGui::ColorEdit4("color", (float*)&rep.static_color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel)) {
                recompute_colors = true;
            }
        }
        ImGui::PopID();
        ImGui::EndGroup();
        ImGui::Spacing();

        if (recompute_colors) {
            compute_atom_colors(rep.colors, *data->mol_data.dynamic.molecule, rep.color_mapping,
                                ImGui::ColorConvertFloat4ToU32(vec_cast(rep.static_color)));
            DynamicArray<bool> mask(data->mol_data.dynamic.molecule->atom_elements.count, false);
            filter::compute_filter_mask(mask, data->mol_data.dynamic, rep.filter.buffer);
            filter::filter_colors(rep.colors, mask);
        }
    }

    ImGui::End();
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

static void draw_statistics(ApplicationData* data) {
    // if (!stats) return;

    ImGui::Begin("Timelines");
    auto group_id = stats::get_group("group1");
    stats::get_property_count(group_id);

    int32 frame_idx = (int32)data->time;

    for (int i = 0; i < stats::get_property_count(group_id); i++) {
        auto prop_id = stats::get_property(group_id, i);
        auto avg_data = stats::get_property_data(prop_id, 0);
        auto count = stats::get_property_data_count(prop_id);
        auto frame = ImGui::BeginPlotFrame(stats::get_property_name(prop_id), ImVec2(0, 100), 0, count, -2.f, 2.f);
        ImGui::PlotFrameLine(frame, "group1", (float*)avg_data, ImGui::FrameLineStyle(), frame_idx);
        int32 new_frame_idx = ImGui::EndPlotFrame(frame, frame_idx);
        if (new_frame_idx != -1) {
            frame_idx = new_frame_idx;
            data->time = (float64)frame_idx;
        }
    }

    // stats::get_group_properties();
    /*
    if (stats->properties.size() > 0) {
            const auto& g = stats->groups.front();
            const auto& i = stats->instances[g.instance_avg_idx];
            const auto& p = stats->properties[i.property_beg_idx];
            auto frame = ImGui::BeginPlotFrame("Prop0", ImVec2(0, 100), 0, p.count, -2.f, 2.f);
            ImGui::PlotFrameLine(frame, "najs", (float*)p.data);
            ImGui::EndPlotFrame(frame);
    }
    */
    ImGui::End();
}

static void draw_ramachandran(ApplicationData* data) {
	constexpr vec2 res(512, 512);
	ImGui::SetNextWindowContentSize(ImVec2(res.x, res.y));
	ImGui::Begin("Ramachandran", &data->ramachandran.enabled, ImGuiWindowFlags_NoFocusOnAppearing);

	int32 num_frames = data->mol_data.dynamic.trajectory ? data->mol_data.dynamic.trajectory->num_frames : 0;
	float range_min = (float)data->ramachandran.frame_range_min;
	float range_max = (float)data->ramachandran.frame_range_max;

	ImGui::SliderFloat("opacity", &data->ramachandran.opacity, 0.f, 2.f);
	ImGui::SliderFloat("radius", &data->ramachandran.radius, 0.1f, 2.f);
	ImGui::RangeSliderFloat("framerange", &range_min, &range_max, 0, math::max(0, num_frames - 1));
	ImGui::SameLine();
	if (ImGui::Button("reset")) {
		range_min = 0;
		range_max = num_frames - 1;
	}
	data->ramachandran.frame_range_min = (int32)range_min;
	data->ramachandran.frame_range_max = (int32)range_max;

	int32 frame = (int32)data->time;
	Array<BackboneAngles> accumulated_angles = get_backbone_angles(data->mol_data.backbone_angles, data->ramachandran.frame_range_min, data->ramachandran.frame_range_max - data->ramachandran.frame_range_min);
	Array<BackboneAngles> current_angles = get_backbone_angles(data->mol_data.backbone_angles, frame);

	ramachandran::compute_accumulation_texture(accumulated_angles, current_angles, data->ramachandran.radius, data->ramachandran.opacity);

	float dim = math::min(ImGui::GetWindowWidth(), ImGui::GetWindowHeight());
	ImVec2 win_pos = ImGui::GetCursorScreenPos();
	ImVec2 canvas_size(dim, dim);
	ImDrawList* dl = ImGui::GetWindowDrawList();

	ImVec2 x0 = win_pos;
	ImVec2 x1(win_pos.x + canvas_size.x, win_pos.y + canvas_size.y);

	dl->ChannelsSplit(2);
	dl->ChannelsSetCurrent(0);
	// ImGui::Image((ImTextureID)ramachandran::segmentation_tex, canvas_size);
	dl->AddImage((ImTextureID)(intptr_t)ramachandran::get_segmentation_texture(), x0, x1);
	dl->ChannelsSetCurrent(1);
	// ImGui::Image((ImTextureID)ramachandran::accumulation_tex, canvas_size);
	dl->AddImage((ImTextureID)(intptr_t)ramachandran::get_accumulation_texture(), x0, x1);
	dl->ChannelsMerge();

	dl->ChannelsSetCurrent(0);

	ImGui::End();
}

static void reset_view(ApplicationData* data) {
    ASSERT(data);
    if (!data->mol_data.dynamic.molecule) return;

    vec3 min_box, max_box;
    compute_bounding_box(&min_box, &max_box, data->mol_data.dynamic.molecule->atom_positions);
    vec3 size = max_box - min_box;
    vec3 cent = (min_box + max_box) * 0.5f;

    data->controller.look_at(cent, cent + size * 2.f);
    data->camera.position = data->controller.position;
    data->camera.orientation = data->controller.orientation;
    data->camera.near_plane = 1.f;
    data->camera.far_plane = math::length(size) * 10.f;
}

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

static void free_mol_data(MoleculeData* data) {
    if (data->dynamic.molecule) {
        free_molecule_structure(data->dynamic.molecule);
        data->dynamic.molecule = nullptr;
    }
    if (data->dynamic.trajectory) {
        FREE(data->dynamic.trajectory);
        data->dynamic.trajectory = nullptr;
    }
    free_backbone_angles_trajectory(&data->backbone_angles);
}

static void load_molecule_data(MoleculeData* mol_data, CString file) {
    ASSERT(mol_data);
    if (file.count > 0) {
        CString ext = get_file_extension(file);
        printf("'%s'\n", ext.beg());
        if (compare_n(ext, "pdb", 3, true)) {
            free_mol_data(mol_data);
            mol_data->dynamic = allocate_and_load_pdb_from_file(file);
            if (!mol_data->dynamic.molecule) {
                printf("ERROR! Failed to load pdb file.\n");
                return;
            }
            mol_data->atom_radii = compute_atom_radii(mol_data->dynamic.molecule->atom_elements);
            if (mol_data->dynamic.trajectory) {
                init_backbone_angles_trajectory(&mol_data->backbone_angles, mol_data->dynamic);
                compute_backbone_angles_trajectory(&mol_data->backbone_angles, mol_data->dynamic);
            }
        } else if (compare_n(ext, "gro", 3, true)) {
            free_mol_data(mol_data);
            mol_data->dynamic.molecule = allocate_and_load_gro_from_file(file);
            mol_data->atom_radii = compute_atom_radii(mol_data->dynamic.molecule->atom_elements);
            if (!mol_data->dynamic.molecule) {
                printf("ERROR! Failed to load gro file.\n");
            }
        } else if (compare_n(ext, "xtc", 3, true)) {
            if (!mol_data->dynamic.molecule) {
                printf("ERROR! Must have molecule loaded before trajectory can be loaded!\n");
            } else {
                mol_data->dynamic.trajectory = allocate_trajectory(file);
                init_backbone_angles_trajectory(&mol_data->backbone_angles, mol_data->dynamic);
                if (mol_data->dynamic.trajectory) {
                    read_trajectory_async(mol_data->dynamic.trajectory, [mol_data]() {
						compute_backbone_angles_trajectory(&mol_data->backbone_angles, mol_data->dynamic);
					});
                }
            }
        } else {
            printf("ERROR! file extension not supported!\n");
        }
    }
}

static void create_default_representation(ApplicationData* data) {
    ASSERT(data);
    if (!data->mol_data.dynamic.molecule) return;

    auto& rep = data->representations.data.push_back({});
    rep.colors.count = data->mol_data.dynamic.molecule->atom_positions.count;
    rep.colors.data = (uint32*)MALLOC(rep.colors.count * sizeof(uint32));
    compute_atom_colors(rep.colors, *data->mol_data.dynamic.molecule, rep.color_mapping);
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
        rep.colors.count = data->mol_data.dynamic.molecule->atom_positions.count;
        rep.colors.data = (uint32*)MALLOC(rep.colors.count * sizeof(uint32));

        compute_atom_colors(rep.colors, *data->mol_data.dynamic.molecule, rep.color_mapping,
                            ImGui::ColorConvertFloat4ToU32(vec_cast(rep.static_color)));
        DynamicArray<bool> mask(data->mol_data.dynamic.molecule->atom_elements.count, false);
        filter::compute_filter_mask(mask, data->mol_data.dynamic, rep.filter.buffer);
        filter::filter_colors(rep.colors, mask);
    }
}

