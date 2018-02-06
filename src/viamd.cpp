#include <imgui.h>
#include <core/platform.h>
#include <core/gl.h>
#include <core/types.h>
#include <core/math_utils.h>
#include <core/camera.h>
#include <core/camera_utils.h>
#include <mol/molecule.h>
#include <mol/trajectory.h>
#include <mol/molecule_utils.h>
#include <mol/pdb_utils.h>
#include <immediate_draw_utils.h>

#include <stdio.h>

struct ApplicationData {
    platform::Window* main_window;
    Camera camera;

    // Perhaps move these into a struct
    MoleculeStructure* mol_struct;
    DynamicArray<float> atom_radii;
    DynamicArray<uint32> atom_colors;
    Trajectory* trajectory;
};

void draw_main_menu(platform::Window* window);

int main(int, char**) {
    ApplicationData data;

	Trackball trackball;

    data.camera.position = vec3(0, 0, 30);

    platform::initialize();
    data.main_window = platform::create_window(1920, 1080, "VIAMD");

    platform::set_vsync(true);

    CString str(PROJECT_SOURCE_DIR "/data/5ulj.pdb");
    auto pdb_res = load_pdb_from_file(str);
    data.mol_struct = &pdb_res.pdb;
    data.atom_radii = molecule::compute_atom_radii(data.mol_struct->atom_elements);
    data.atom_colors = molecule::compute_atom_colors(*data.mol_struct, ColorMapping::CPK);

    immediate::initialize();
    molecule::draw::initialize();

    // Setup style
    ImGui::StyleColorsClassic();

    bool show_demo_window = true;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    // Main loop
    while (!(platform::window_should_close(data.main_window))) {
        platform::update();
        platform::InputState* input = platform::get_input_state();
        float dt = platform::get_delta_time();
		int display_w, display_h;
		platform::get_framebuffer_size(data.main_window, &display_w, &display_h);

		ImGui::Begin("Input");
		ImGui::Text("MouseVel: %g, %g", input->mouse_velocity.x, input->mouse_velocity.y);
		ImGui::Text("Look at: %g, %g, %g", trackball.target.x, trackball.target.y, trackball.target.z);
		ImGui::Text("Camera Pos: %g, %g, %g", data.camera.position.x, data.camera.position.y, data.camera.position.z);

		ImGui::Text("MOUSE_BUTTONS [%i, %i, %i, %i, %i]", input->mouse_down[0], input->mouse_down[1], input->mouse_down[2], input->mouse_down[3], input->mouse_down[4]);

		ImGui::End();

		if (!ImGui::GetIO().WantCaptureKeyboard) {
			//camera_controller_trackball(&data.camera, input->mouse_down[0], input->mouse_down[1], input->key_down[2],
			//	input->mouse_velocity, display_w, display_h);
			trackball.update(&data.camera, input->mouse_down[0], input->mouse_down[2], input->mouse_down[1],
				input->mouse_velocity, input->mouse_scroll.y, display_w, display_h);
		}

        // MAIN MENU BAR
        draw_main_menu(data.main_window);

        // 3. Show the ImGui demo window. Most of the sample code is in ImGui::ShowDemoWindow().
        if (show_demo_window) {
            ImGui::SetNextWindowPos(ImVec2(650, 20),
                                    ImGuiCond_FirstUseEver);  // Normally user code doesn't need/want to call this because positions are saved in .ini
                                                              // file anyway. Here we just want to make the demo initial state a bit more friendly!
            ImGui::ShowDemoWindow(&show_demo_window);
        }

        // Rendering
        glViewport(0, 0, display_w, display_h);
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glDepthMask(GL_TRUE);
		glDepthFunc(GL_LESS);

        mat4 model_mat = mat4(1);
        mat4 view_mat = compute_world_to_view_matrix(data.camera);
        mat4 proj_mat = compute_perspective_projection_matrix(data.camera, display_w, display_h);
        mat4 mvp = proj_mat * view_mat;

		srand(0);
		for (int i = 0; i < 500; i++) {
			vec3 c = math::ballRand(100.f);
			vec3 v0 = c + mat3(view_mat) * vec3(-1, 0, 0);
			vec3 v1 = c + mat3(view_mat) * vec3(1, 0, 0);
			vec3 v2 = c + mat3(view_mat) * vec3(0.5f, 1, 0);
			unsigned char color[4] = { 255, 255, 255, 255 };
			immediate::draw_triangle(&v0[0], &v1[0], &v2[0], color);
		}
        float p0[3] = {-1, 0, 0};
        float p1[3] = {1, 0, 0};
        float p2[3] = {0, 2, 0};

        unsigned char c[4] = {255, 0, 0, 255};
		//immediate::draw_point(p0, c);
		//immediate::draw_point(p1, c);
		//immediate::draw_point(p2, c);
        immediate::draw_triangle(p0, p1, p2, c);
        immediate::flush(&mvp[0][0]);

        //molecule::draw::draw_vdw(data.mol_struct->atom_positions, data.atom_radii, data.atom_colors, model_mat, view_mat, proj_mat);

        ImGui::Render();
        platform::swap_buffers(data.main_window);
    }

    platform::shutdown();

    return 0;
}

void draw_main_menu(platform::Window* main_window) {
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
                platform::set_window_should_close(main_window, true);
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