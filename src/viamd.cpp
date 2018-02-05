#include <imgui.h>
#include <core/platform.h>
#include <core/gl.h>
#include <core/types.h>
#include <core/math_utils.h>
#include <core/camera.h>
#include <core/camera_utils.h>
#include <mol/molecule.h>
#include <mol/trajectory.h>
#include <immediate_draw_utils.h>

#include <stdio.h>

struct ApplicationData {
    platform::Window*   main_window;
    Camera              camera;
    MoleculeStructure*  mol_struct;
    Trajectory*         trajectory;
};

void draw_main_menu(platform::Window* window);

int main(int, char**)
{
    ApplicationData data;

    data.camera.position = vec3(0,0,3);

    platform::initialize();
    data.main_window = platform::create_window(1024, 768, "VIAMD");

    immediate::initialize();

    // Setup style
    ImGui::StyleColorsClassic();

    bool show_demo_window = true;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    // Main loop
    while (!(platform::window_should_close(data.main_window))) {
        platform::update();

		// MAIN MENU BAR
		draw_main_menu(data.main_window);

        // 3. Show the ImGui demo window. Most of the sample code is in ImGui::ShowDemoWindow().
        if (show_demo_window)
        {
            ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiCond_FirstUseEver); // Normally user code doesn't need/want to call this because positions are saved in .ini file anyway. Here we just want to make the demo initial state a bit more friendly!
            ImGui::ShowDemoWindow(&show_demo_window);
        }

        // Rendering
        int display_w, display_h;
        platform::get_framebuffer_size(data.main_window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        mat4 view_mat = glm::translate(mat4(), -data.camera.position);
        mat4 proj_mat = compute_perspective_projection_matrix(data.camera, display_w, display_h);
        mat4 mvp = proj_mat * view_mat;

        float p0[3] = {0,0,0};
        float p1[3] = {1,0,0};
        float p2[3] = {0.5,1,0};

        unsigned char c[4] = {255,255,255,255};
        immediate::draw_triangle(p0, p1, p2, c);
        immediate::flush(&mvp[0][0]);

        ImGui::Render();
        platform::swap_buffers(data.main_window);
    }

    platform::shutdown();

    return 0;
}


void draw_main_menu(platform::Window* main_window) {
    bool new_clicked = false;
	if (ImGui::BeginMainMenuBar())
	{
		if (ImGui::BeginMenu("File"))
		{
			if (ImGui::MenuItem("New", "CTRL+N"))
				new_clicked = true;

			if (ImGui::MenuItem("Open", "CTRL+O")) {}
			if (ImGui::BeginMenu("Open Recent"))
			{
				ImGui::MenuItem("fish_hat.c");
				ImGui::MenuItem("fish_hat.inl");
				ImGui::MenuItem("fish_hat.h");
				ImGui::EndMenu();
			}
			if (ImGui::MenuItem("Save", "CTRL+S")) {}
			if (ImGui::MenuItem("Save As..")) {}
            if (ImGui::MenuItem("Export", "CTRL+E")) {}
			ImGui::Separator();
			if (ImGui::MenuItem("Quit", "ALT+F4")) {
				platform::set_window_should_close(main_window, true);
			}
			ImGui::EndMenu();
		}
		if (ImGui::BeginMenu("Edit"))
		{
			if (ImGui::MenuItem("Undo", "CTRL+Z")) {}
			if (ImGui::MenuItem("Redo", "CTRL+Y", false, false)) {}  // Disabled item
			ImGui::Separator();
			if (ImGui::MenuItem("Cut", "CTRL+X")) {}
			if (ImGui::MenuItem("Copy", "CTRL+C")) {}
			if (ImGui::MenuItem("Paste", "CTRL+V")) {}
			ImGui::EndMenu();
		}
		ImGui::EndMainMenuBar();
	}

    if (new_clicked)
        ImGui::OpenPopup("Warning New");
    if (ImGui::BeginPopupModal("Warning New")) {
        ImGui::Text("By creating a new workspace you will loose any unsaved progress.");
        ImGui::Text("Are you sure?");
        ImGui::Separator();
        if (ImGui::Button("OK", ImVec2(120,0))) { ImGui::CloseCurrentPopup(); }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120,0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
    }
}