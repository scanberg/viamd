#include <imgui.h>
#include <core/platform.h>
#include <core/gl.h>
#include <core/types.h>
#include <core/math_utils.h>
#include <core/camera.h>
#include <core/camera_utils.h>
#include <core/string_utils.h>
#include <mol/molecule.h>
#include <mol/trajectory.h>
#include <mol/trajectory_utils.h>
#include <mol/molecule_utils.h>
#include <mol/pdb_utils.h>
#include <mol/gro_utils.h>
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

struct MainFramebuffer {
	GLuint id = 0;
	GLuint tex_depth = 0;
	GLuint tex_color = 0;
	GLuint tex_picking = 0;
	int width = 0;
	int height = 0;
};

constexpr unsigned int NO_PICKING_IDX = 0xffffffff;

enum struct PlaybackInterpolationMode {
	NEAREST_FRAME,
	LINEAR,
	LINEAR_PERIODIC,
	CUBIC,
	CUBIC_PERIODIC
};

struct ApplicationData {
	// --- PLATFORM ---
    platform::Window* main_window;

	// --- CAMERA ---
    Camera camera;
    TrackballController controller;

    // --- MOL DATA ---
    MoleculeStructure* mol_struct;
    Trajectory* trajectory;
    DynamicArray<float> atom_radii;
    DynamicArray<uint32> atom_colors;

	// Framebuffer
	MainFramebuffer fbo;
    unsigned int picking_idx = NO_PICKING_IDX;

	// --- PLAYBACK ---
	float64 time = 0.f; 	// needs to be double precision
	float frames_per_second = 10.f;
    bool is_playing = false;
	PlaybackInterpolationMode interpolation = PlaybackInterpolationMode::LINEAR_PERIODIC;

	// --- VISUALS ---
	// SSAO
    bool use_ssao = false;
	float ssao_intensity = 1.5f;
	float ssao_radius = 6.f;

    bool show_console = false;
};

static float compute_avg_ms(float dt);
static void draw_main_menu(ApplicationData* data);
static void draw_console(ApplicationData* data, int width, int height, float dt);
static void draw_atom_info(const MoleculeStructure& mol, int atom_idx, int x, int y);
static void init_main_framebuffer(MainFramebuffer* fbo, int width, int height);
static void destroy_main_framebuffer(MainFramebuffer* fbo);
static void reset_view(ApplicationData* data);

int main(int, char**) {
	ApplicationData data;

    int display_w = 1920;
    int display_h = 1080;

	float radii_scale = 1.0f;

    //auto gro_res = load_gro_from_file(PROJECT_SOURCE_DIR "/data/bta-gro/20-mol-p.gro");
    //auto gro_res = load_gro_from_file(PROJECT_SOURCE_DIR "/data/peptides/box_2.gro");
    //auto gro_res = load_gro_from_file(PROJECT_SOURCE_DIR "/data/shaoqi/md-nowater.gro");
	//auto gro_res = load_gro_from_file(PROJECT_SOURCE_DIR "/data/peptides/box_2.gro");
	//auto gro_res = load_gro_from_file(PROJECT_SOURCE_DIR "/data/amyloid/centered.gro");
	//auto gro_res = load_gro_from_file(PROJECT_SOURCE_DIR "/data/water/water.gro");
	auto gro_res = load_gro_from_file(PROJECT_SOURCE_DIR "/data/amyloid-6T/conf-60-6T.gro");
    //auto pdb_res = load_pdb_from_file(PROJECT_SOURCE_DIR "/data/5ulj.pdb");

    //Trajectory* traj = allocate_trajectory(PROJECT_SOURCE_DIR "/data/bta-gro/traj-centered.xtc");
    //Trajectory* traj = allocate_trajectory(PROJECT_SOURCE_DIR "/data/peptides/md_0_1_noPBC_2.xtc");
	//Trajectory* traj = allocate_trajectory(PROJECT_SOURCE_DIR "/data/shaoqi/md-centered.xtc");
	//Trajectory* traj = allocate_trajectory(PROJECT_SOURCE_DIR "/data/peptides/md_0_1_noPBC_2.xtc");
	//Trajectory* traj = allocate_trajectory(PROJECT_SOURCE_DIR "/data/amyloid/centered.xtc");
	Trajectory* traj = allocate_trajectory(PROJECT_SOURCE_DIR "/data/amyloid-6T/prod-centered.xtc");

    if (traj) {
        MoleculeStructure* mol = &gro_res.gro;
        read_trajectory_async(traj, [traj, mol]() {
            auto backbone = compute_backbone(mol->chains[0], mol->residues, mol->atom_labels);
            auto traj_angles = compute_backbone_angles_trajectory(*traj, backbone);
            for (int i = 0; i < traj_angles.num_frames; i++) {
                auto angles = get_backbone_angles(traj_angles, i);
                printf("omega  phi   psi\n");
                for (const auto& ba : angles) {
                    printf("% 6.1f % 6.1f % 6.1f\n", ba.omega * math::RAD_TO_DEG, ba.phi * math::RAD_TO_DEG, ba.psi * math::RAD_TO_DEG);
                }
            }
        });
    }
	
    data.mol_struct = &gro_res.gro;
    data.trajectory = traj;

	DynamicArray<BackboneSegment> backbone;
	DynamicArray<SplineSegment> spline;

    if (data.trajectory && data.trajectory->num_frames > 0)
        copy_trajectory_positions(data.mol_struct->atom_positions, *data.trajectory, 0);
    reset_view(&data);

	if (data.mol_struct->chains.count > 0) {
		backbone = compute_backbone(data.mol_struct->chains[0], data.mol_struct->residues, data.mol_struct->atom_labels);
		spline = compute_spline(data.mol_struct->atom_positions, backbone, 8);
		//auto backbone_angles = compute_backbone_angles(data.mol_struct->atom_positions, backbone);
		//printf("omega  phi   psi\n");
		//for (const auto& ba : backbone_angles) {
		//	printf("% 6.1f % 6.1f % 6.1f\n", ba.omega * math::RAD_TO_DEG, ba.phi * math::RAD_TO_DEG, ba.psi * math::RAD_TO_DEG);
		//}
	}

    platform::initialize();
    data.main_window = platform::create_window(display_w, display_h, "VIAMD");

    platform::set_vsync(false);

    data.atom_radii = compute_atom_radii(data.mol_struct->atom_elements);
    data.atom_colors = compute_atom_colors(*data.mol_struct, ColorMapping::CPK);

    immediate::initialize();
    postprocessing::initialize(display_w, display_h);
    draw::initialize();
	init_main_framebuffer(&data.fbo, display_w, display_h);

    // Setup style
    ImGui::StyleColorsClassic();

    bool show_demo_window = false;
    vec4 clear_color = vec4(0.6, 0.6, 0.6, 1);
    vec4 clear_index = vec4(1, 1, 1, 1);

    // Main loop
    while (!(platform::window_should_close(data.main_window))) {
        platform::update();
        platform::InputState* input = platform::get_input_state();
		platform::get_framebuffer_size(data.main_window, &display_w, &display_h);
		float dt = (float)platform::get_delta_time();

		if (data.fbo.width != display_w || data.fbo.height != display_h) {
			init_main_framebuffer(&data.fbo, display_w, display_h);
			postprocessing::initialize(display_w, display_h);
		}

        {
            int x = input->mouse_screen_coords.x;
            int y = input->mouse_screen_coords.y;
            unsigned char color[4];
            glBindFramebuffer(GL_READ_FRAMEBUFFER, data.fbo.id);
            glReadBuffer(GL_COLOR_ATTACHMENT1);
            glReadPixels(x, display_h - y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, color);
            glReadBuffer(GL_NONE);
            glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
            data.picking_idx = color[0] + (color[1] << 8) + (color[2] << 16) + (color[3] << 24);
        }

		if (input->key_hit[CONSOLE_KEY]) data.show_console = !data.show_console;

        bool time_changed = false;

		float ms = compute_avg_ms(dt);

		ImGui::Begin("Misc");
		ImGui::Text("%.2f ms (%.1f fps)", ms, 1000.f / (ms));
		ImGui::Text("MouseVel: %g, %g", input->mouse_velocity.x, input->mouse_velocity.y);
		ImGui::Text("Camera Pos: %g, %g, %g", data.camera.position.x, data.camera.position.y, data.camera.position.z);
		ImGui::Text("Mouse Buttons: [%i, %i, %i, %i, %i]", input->mouse_down[0], input->mouse_down[1], input->mouse_down[2], input->mouse_down[3], input->mouse_down[4]);
		ImGui::SliderFloat("Radius Scale", &radii_scale, 0.1f, 2.f);
        if (ImGui::Button("Reset View")) {
            reset_view(&data);
        }
		if (data.trajectory) {
			ImGui::Text("Num Frames: %i", data.trajectory->num_frames);
			float t = (float)data.time;
			if (ImGui::SliderFloat("Time", &t, 0, (float)(data.trajectory->num_frames - 1))) {
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
            data.time += dt * data.frames_per_second;
			time_changed = true;
        }

        if (data.trajectory && time_changed) {
            data.time = math::clamp(data.time, 0.0, float64(data.trajectory->num_frames - 1));
			if (data.time == float64(data.trajectory->num_frames - 1)) data.is_playing = false;

			int prev_frame_idx = math::max(0, (int)data.time);
			int next_frame_idx = math::min(prev_frame_idx + 1, data.trajectory->num_frames - 1);
			if (prev_frame_idx == next_frame_idx) {
				copy_trajectory_positions(data.mol_struct->atom_positions, *data.trajectory, prev_frame_idx);
			}
			else {
				// INTERPOLATE
				auto prev_frame = get_trajectory_frame(*data.trajectory, prev_frame_idx);
				auto next_frame = get_trajectory_frame(*data.trajectory, next_frame_idx);

				float t = (float)math::fract(data.time);
				linear_interpolation_periodic(data.mol_struct->atom_positions, prev_frame.atom_positions, next_frame.atom_positions, t, prev_frame.box);
			}
        }

		if (!ImGui::GetIO().WantCaptureMouse) {
            data.controller.input.rotate_button = input->mouse_down[0];
            data.controller.input.pan_button = input->mouse_down[1];
            data.controller.input.dolly_button = input->mouse_down[2];
            data.controller.input.prev_mouse_ndc = input->prev_mouse_ndc_coords;
            data.controller.input.curr_mouse_ndc = input->mouse_ndc_coords;
            data.controller.input.dolly_delta = input->mouse_scroll.y;
			data.controller.update();
			data.camera.position = data.controller.position;
			data.camera.orientation = data.controller.orientation;

            if (data.picking_idx != NO_PICKING_IDX) {
                ivec2 pos = input->mouse_screen_coords;
                draw_atom_info(*data.mol_struct, data.picking_idx, pos.x, pos.y);
            }
		}

        // GUI ELEMENTS
		draw_console(&data, display_w, display_h, dt);
        draw_main_menu(&data);

        // 3. Show the ImGui demo window. Most of the sample code is in ImGui::ShowDemoWindow().
        if (show_demo_window) {
            ImGui::SetNextWindowPos(ImVec2(650, 20),
                                    ImGuiCond_FirstUseEver);  // Normally user code doesn't need/want to call this because positions are saved in .ini
                                                              // file anyway. Here we just want to make the demo initial state a bit more friendly!
            ImGui::ShowDemoWindow(&show_demo_window);
        }

		mat4 model_mat = mat4(1);
		mat4 view_mat = compute_world_to_view_matrix(data.camera);
		mat4 proj_mat = compute_perspective_projection_matrix(data.camera, display_w, display_h);

        // Rendering
        glViewport(0, 0, display_w, display_h);

		const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1 };
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, data.fbo.id);

		// Clear picking buffer
		glDrawBuffer(GL_COLOR_ATTACHMENT1);
		glClearColor(clear_index.x, clear_index.y, clear_index.z, clear_index.w);
		glClear(GL_COLOR_BUFFER_BIT);

		// Clear color buffer
		glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
		glDrawBuffer(GL_COLOR_ATTACHMENT0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Enable both color and picking
		glDrawBuffers(2, draw_buffers);
		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LESS);

        draw::draw_vdw(data.mol_struct->atom_positions, data.atom_radii, data.atom_colors, view_mat, proj_mat, radii_scale);
		//draw::draw_licorice(data.mol_struct->atom_positions, data.mol_struct->bonds, data.atom_colors, view_mat, proj_mat, radii_scale);
		//draw::draw_ribbons(spline, view_mat, proj_mat);

        // Activate backbuffer
        glDisable(GL_DEPTH_TEST);
		glDepthFunc(GL_ALWAYS);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
		glDrawBuffer(GL_BACK);
        glClear(GL_COLOR_BUFFER_BIT);

        // Apply tone mapping
        postprocessing::apply_tonemapping(data.fbo.tex_color);

		draw::draw_backbone(backbone, data.mol_struct->atom_positions, view_mat, proj_mat);
		//draw::draw_spline(spline, view_mat, proj_mat);

		if (data.use_ssao) {
			postprocessing::apply_ssao(data.fbo.tex_depth, proj_mat, data.ssao_intensity, data.ssao_radius);
		}

        // Render Imgui
        ImGui::Render();

        // Swap buffers
        platform::swap_buffers(data.main_window);
    }

    platform::shutdown();

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
}

static void draw_main_menu(ApplicationData* data) {
    ASSERT(data);
    bool new_clicked = false;
    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("New", "CTRL+N")) new_clicked = true;

            if (ImGui::MenuItem("Open", "CTRL+O")) {
            }
            if (ImGui::BeginMenu("Open Recent")) {
                ImGui::MenuItem("fish_hat.c");
                ImGui::MenuItem("fish_hat.inl");
                ImGui::MenuItem("fish_hat.h");
                ImGui::EndMenu();
            }
            if (ImGui::MenuItem("Save", "CTRL+S")) {
            }
            if (ImGui::MenuItem("Save As..")) {
            }
            if (ImGui::MenuItem("Export", "CTRL+E")) {
            }
            ImGui::Separator();
            if (ImGui::MenuItem("Quit", "ALT+F4")) {
                platform::set_window_should_close(data->main_window, true);
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
            ImGui::Checkbox("SSAO", &data->use_ssao);
			ImGui::SliderFloat("Intensity", &data->ssao_intensity, 0.5f, 6.f);
			ImGui::SliderFloat("Radius", &data->ssao_radius, 1.f, 30.f);
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

// Demonstrating creating a simple console window, with scrolling, filtering, completion and history.
// For the console example, here we are using a more C++ like approach of declaring a class to hold the data and the functions.
struct Console
{
    char                  InputBuf[256];
    ImVector<char*>       Items;
    bool                  ScrollToBottom;
    ImVector<char*>       History;
    int                   HistoryPos;    // -1: new line, 0..History.Size-1 browsing history.
    ImVector<const char*> Commands;
	float				  YPos = -10000;

    Console()
    {
        ClearLog();
        memset(InputBuf, 0, sizeof(InputBuf));
        HistoryPos = -1;
        Commands.push_back("HELP");
        Commands.push_back("HISTORY");
        Commands.push_back("CLEAR");
        Commands.push_back("CLASSIFY");  // "classify" is here to provide an example of "C"+[tab] completing to "CL" and displaying matches.
        AddLog("Welcome to ImGui!");
    }
    ~Console()
    {
        ClearLog();
        for (int i = 0; i < History.Size; i++)
            free(History[i]);
    }

    // Portable helpers
    static int   Stricmp(const char* str1, const char* str2)         { int d; while ((d = toupper(*str2) - toupper(*str1)) == 0 && *str1) { str1++; str2++; } return d; }
    static int   Strnicmp(const char* str1, const char* str2, int n) { int d = 0; while (n > 0 && (d = toupper(*str2) - toupper(*str1)) == 0 && *str1) { str1++; str2++; n--; } return d; }
    static char* Strdup(const char *str)                             { size_t len = strlen(str) + 1; void* buff = malloc(len); return (char*)memcpy(buff, (const void*)str, len); }

    void ClearLog()
    {
        for (int i = 0; i < Items.Size; i++)
            free(Items[i]);
        Items.clear();
        ScrollToBottom = true;
    }

    void AddLog(const char* fmt, ...) IM_FMTARGS(2)
    {
        // FIXME-OPT
        char buf[1024];
        va_list args;
        va_start(args, fmt);
        vsnprintf(buf, IM_ARRAYSIZE(buf), fmt, args);
        buf[IM_ARRAYSIZE(buf)-1] = 0;
        va_end(args);
        Items.push_back(Strdup(buf));
        ScrollToBottom = true;
    }

    void Draw(const char* title, ApplicationData* data, int width, int height, float dt)
    {
		constexpr int WINDOW_FLAGS = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoScrollbar |
			ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_Modal;
		
		float console_width = width;
		float console_height = height * 0.75f;

		// Quarter of a second to hide/show
		float speed = console_height * dt * 4.f;

		float target_hide_y = -console_height;
		float target_show_y = 0;
		float target_y = data->show_console ? target_show_y : target_hide_y;

		if (YPos != target_y) {
			float delta = target_y < YPos ? -speed : speed;
			YPos = math::clamp(YPos + delta, target_hide_y, target_show_y);
		}

		bool console_fully_shown = (YPos == target_show_y);
		bool console_fully_hidden = (YPos == target_hide_y);

		if (console_fully_hidden) return;

		ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.f);

		ImGui::SetNextWindowSize(ImVec2(console_width, console_height));
        ImGui::SetNextWindowPos(ImVec2(0.f, YPos));
		ImGui::Begin(title, &data->show_console, WINDOW_FLAGS);


        // As a specific feature guaranteed by the library, after calling Begin() the last Item represent the title bar. So e.g. IsItemHovered() will return true when hovering the title bar.
        // Here we create a context menu only available from the title bar.
        //if (ImGui::BeginPopupContextItem())
        //{
        //    if (ImGui::MenuItem("Close"))
        //        *p_open = false;
        // ImGui::EndPopup();
        //}

        //ImGui::TextWrapped("This example implements a console with basic coloring, completion and history. A more elaborate implementation may want to store entries along with extra data such as timestamp, emitter, etc.");
        //ImGui::TextWrapped("Enter 'HELP' for help, press TAB to use text completion.");

        // TODO: display items starting from the bottom

        //if (ImGui::SmallButton("Add Dummy Text")) { AddLog("%d some text", Items.Size); AddLog("some more text"); AddLog("display very important message here!"); } ImGui::SameLine();
        //if (ImGui::SmallButton("Add Dummy Error")) { AddLog("[error] something went wrong"); } ImGui::SameLine();
        //if (ImGui::SmallButton("Clear")) { ClearLog(); } ImGui::SameLine();
        //bool copy_to_clipboard = ImGui::SmallButton("Copy"); ImGui::SameLine();
        //if (ImGui::SmallButton("Scroll to bottom")) ScrollToBottom = true;
        //static float t = 0.0f; if (ImGui::GetTime() - t > 0.02f) { t = ImGui::GetTime(); AddLog("Spam %f", t); }

        //ImGui::Separator();

        //ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0,0));
        //static ImGuiTextFilter filter;
        //filter.Draw("Filter (\"incl,-excl\") (\"error\")", 180);
        //ImGui::PopStyleVar();
        //ImGui::Separator();

        const float footer_height_to_reserve = ImGui::GetStyle().ItemSpacing.y + ImGui::GetFrameHeightWithSpacing(); // 1 separator, 1 input text
        ImGui::BeginChild("ScrollingRegion", ImVec2(0, -footer_height_to_reserve), false, ImGuiWindowFlags_HorizontalScrollbar); // Leave room for 1 separator + 1 InputText
        //if (ImGui::BeginPopupContextWindow())
        //{
        //    if (ImGui::Selectable("Clear")) ClearLog();
        //    ImGui::EndPopup();
        //}

        // Display every line as a separate entry so we can change their color or add custom widgets. If you only want raw text you can use ImGui::TextUnformatted(log.begin(), log.end());
        // NB- if you have thousands of entries this approach may be too inefficient and may require user-side clipping to only process visible items.
        // You can seek and display only the lines that are visible using the ImGuiListClipper helper, if your elements are evenly spaced and you have cheap random access to the elements.
        // To use the clipper we could replace the 'for (int i = 0; i < Items.Size; i++)' loop with:
        //     ImGuiListClipper clipper(Items.Size);
        //     while (clipper.Step())
        //         for (int i = clipper.DisplayStart; i < clipper.DisplayEnd; i++)
        // However take note that you can not use this code as is if a filter is active because it breaks the 'cheap random-access' property. We would need random-access on the post-filtered list.
        // A typical application wanting coarse clipping and filtering may want to pre-compute an array of indices that passed the filtering test, recomputing this array when user changes the filter,
        // and appending newly elements as they are inserted. This is left as a task to the user until we can manage to improve this example code!
        // If your items are of variable size you may want to implement code similar to what ImGuiListClipper does. Or split your data into fixed height items to allow random-seeking into your list.
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(4,1)); // Tighten spacing
        //if (copy_to_clipboard)
        //    ImGui::LogToClipboard();
        ImVec4 col_default_text = ImGui::GetStyleColorVec4(ImGuiCol_Text);
        for (int i = 0; i < Items.Size; i++) {
            const char* item = Items[i];
            //if (!filter.PassFilter(item))
            //    continue;
            ImVec4 col = col_default_text;
            if (strstr(item, "[error]")) col = ImColor(1.0f,0.4f,0.4f,1.0f);
            else if (strncmp(item, "# ", 2) == 0) col = ImColor(1.0f,0.78f,0.58f,1.0f);
            ImGui::PushStyleColor(ImGuiCol_Text, col);
            ImGui::TextUnformatted(item);
            ImGui::PopStyleColor();
        }
        //if (copy_to_clipboard)
        //    ImGui::LogFinish();
        if (ScrollToBottom)
            ImGui::SetScrollHere();
        ScrollToBottom = false;
        ImGui::PopStyleVar();
        ImGui::EndChild();
        ImGui::Separator();

		if (console_fully_shown) {
			// Command-line
			ImGui::PushItemWidth(-1);
			if (ImGui::InputText("##Input", InputBuf, IM_ARRAYSIZE(InputBuf), ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CallbackCompletion | ImGuiInputTextFlags_CallbackHistory, &TextEditCallbackStub, (void*)this))
			{
				char* input_end = InputBuf + strlen(InputBuf);
				while (input_end > InputBuf && input_end[-1] == ' ') { input_end--; } *input_end = 0;
				if (InputBuf[0])
					ExecCommand(InputBuf);
				strcpy(InputBuf, "");
			}
			ImGui::PopItemWidth();
			ImGui::SetKeyboardFocusHere(-1);
		}

        ImGui::End();
		ImGui::PopStyleVar();
    }

    void ExecCommand(const char* command_line)
    {
        AddLog("# %s\n", command_line);

        // Insert into history. First find match and delete it so it can be pushed to the back. This isn't trying to be smart or optimal.
        HistoryPos = -1;
        for (int i = History.Size-1; i >= 0; i--)
            if (Stricmp(History[i], command_line) == 0)
            {
                free(History[i]);
                History.erase(History.begin() + i);
                break;
            }
        History.push_back(Strdup(command_line));

        // Process command
        if (Stricmp(command_line, "CLEAR") == 0)
        {
            ClearLog();
        }
        else if (Stricmp(command_line, "HELP") == 0)
        {
            AddLog("Commands:");
            for (int i = 0; i < Commands.Size; i++)
                AddLog("- %s", Commands[i]);
        }
        else if (Stricmp(command_line, "HISTORY") == 0)
        {
            int first = History.Size - 10;
            for (int i = first > 0 ? first : 0; i < History.Size; i++)
                AddLog("%3d: %s\n", i, History[i]);
        }
        else
        {
            AddLog("Unknown command: '%s'\n", command_line);
        }
    }

    static int TextEditCallbackStub(ImGuiTextEditCallbackData* data) // In C++11 you are better off using lambdas for this sort of forwarding callbacks
    {
        Console* console = (Console*)data->UserData;
        return console->TextEditCallback(data);
    }

    int TextEditCallback(ImGuiTextEditCallbackData* data)
    {
        //AddLog("cursor: %d, selection: %d-%d", data->CursorPos, data->SelectionStart, data->SelectionEnd);
        switch (data->EventFlag)
        {
        case ImGuiInputTextFlags_CallbackCompletion:
            {
                // Example of TEXT COMPLETION

                // Locate beginning of current word
                const char* word_end = data->Buf + data->CursorPos;
                const char* word_start = word_end;
                while (word_start > data->Buf)
                {
                    const char c = word_start[-1];
                    if (c == ' ' || c == '\t' || c == ',' || c == ';')
                        break;
                    word_start--;
                }

                // Build a list of candidates
                ImVector<const char*> candidates;
                for (int i = 0; i < Commands.Size; i++)
                    if (Strnicmp(Commands[i], word_start, (int)(word_end-word_start)) == 0)
                        candidates.push_back(Commands[i]);

                if (candidates.Size == 0)
                {
                    // No match
                    AddLog("No match for \"%.*s\"!\n", (int)(word_end-word_start), word_start);
                }
                else if (candidates.Size == 1)
                {
                    // Single match. Delete the beginning of the word and replace it entirely so we've got nice casing
                    data->DeleteChars((int)(word_start-data->Buf), (int)(word_end-word_start));
                    data->InsertChars(data->CursorPos, candidates[0]);
                    data->InsertChars(data->CursorPos, " ");
                }
                else
                {
                    // Multiple matches. Complete as much as we can, so inputing "C" will complete to "CL" and display "CLEAR" and "CLASSIFY"
                    int match_len = (int)(word_end - word_start);
                    for (;;)
                    {
                        int c = 0;
                        bool all_candidates_matches = true;
                        for (int i = 0; i < candidates.Size && all_candidates_matches; i++)
                            if (i == 0)
                                c = toupper(candidates[i][match_len]);
                            else if (c == 0 || c != toupper(candidates[i][match_len]))
                                all_candidates_matches = false;
                        if (!all_candidates_matches)
                            break;
                        match_len++;
                    }

                    if (match_len > 0)
                    {
                        data->DeleteChars((int)(word_start - data->Buf), (int)(word_end-word_start));
                        data->InsertChars(data->CursorPos, candidates[0], candidates[0] + match_len);
                    }

                    // List matches
                    AddLog("Possible matches:\n");
                    for (int i = 0; i < candidates.Size; i++)
                        AddLog("- %s\n", candidates[i]);
                }

                break;
            }
        case ImGuiInputTextFlags_CallbackHistory:
            {
                // Example of HISTORY
                const int prev_history_pos = HistoryPos;
                if (data->EventKey == ImGuiKey_UpArrow)
                {
                    if (HistoryPos == -1)
                        HistoryPos = History.Size - 1;
                    else if (HistoryPos > 0)
                        HistoryPos--;
                }
                else if (data->EventKey == ImGuiKey_DownArrow)
                {
                    if (HistoryPos != -1)
                        if (++HistoryPos >= History.Size)
                            HistoryPos = -1;
                }

                // A better implementation would preserve the data on the current input line along with cursor position.
                if (prev_history_pos != HistoryPos)
                {
                    data->CursorPos = data->SelectionStart = data->SelectionEnd = data->BufTextLen = (int)snprintf(data->Buf, (size_t)data->BufSize, "%s", (HistoryPos >= 0) ? History[HistoryPos] : "");
                    data->BufDirty = true;
                }
            }
        }
        return 0;
    }
};

static void draw_console(ApplicationData* data, int width, int height, float dt) {
    static Console console;
    console.Draw("Console", data, width, height, dt);
}

static void draw_atom_info(const MoleculeStructure& mol, int atom_idx, int x, int y) {
    
    // @TODO: Assert things and make this failproof
	if (atom_idx < 0 || atom_idx >= mol.atom_positions.count) return;

    int res_idx = mol.atom_residue_indices[atom_idx];
    const Residue& res = mol.residues[res_idx];
    const char* res_id = res.id;
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
        ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
        ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar |
        ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoInputs |
        ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoFocusOnAppearing);
    ImGui::Text("%s", buff);
    ImGui::End();
    ImGui::PopStyleColor();
}


static void reset_view(ApplicationData* data) {
    ASSERT(data);
    ASSERT(data->mol_struct);

    vec3 min_box, max_box;
    compute_bounding_box(&min_box, &max_box, data->mol_struct->atom_positions);
    vec3 size = max_box - min_box;
    vec3 cent = (min_box + max_box) * 0.5f;

	data->controller.look_at(cent, cent + size * 2.f);
    data->camera.position = data->controller.position;
    data->camera.orientation = data->controller.orientation;
	data->camera.near_plane = 1.f;
	data->camera.far_plane = math::length(size) * 3.f;
}

static void init_main_framebuffer(MainFramebuffer* fbo, int width, int height) {
	ASSERT(fbo);

	bool attach_textures = false;
	if (!fbo->id) {
		glGenFramebuffers(1, &fbo->id);
		attach_textures = true;
	}

	if (!fbo->tex_depth)
		glGenTextures(1, &fbo->tex_depth);
	if (!fbo->tex_color)
		glGenTextures(1, &fbo->tex_color);
	if (!fbo->tex_picking)
		glGenTextures(1, &fbo->tex_picking);

	glBindTexture(GL_TEXTURE_2D, fbo->tex_depth);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	glBindTexture(GL_TEXTURE_2D, fbo->tex_color);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, width, height, 0, GL_RGBA, GL_FLOAT, nullptr);
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

	if (attach_textures) {
		glBindFramebuffer(GL_FRAMEBUFFER, fbo->id);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, fbo->tex_depth, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, fbo->tex_color, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, fbo->tex_picking, 0);

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}
}

static void destroy_main_framebuffer(MainFramebuffer* fbo) {
	ASSERT(fbo);
	if (fbo->id) glDeleteFramebuffers(1, &fbo->id);
	if (fbo->tex_depth) glDeleteTextures(1, &fbo->tex_depth);
	if (fbo->tex_color) glDeleteTextures(1, &fbo->tex_color);
	if (fbo->tex_picking) glDeleteTextures(1, &fbo->tex_picking);
}