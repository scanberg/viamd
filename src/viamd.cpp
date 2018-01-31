#include <imgui.h>
#include <core/platform.h>
#include <core/types.h>
#include <core/math_utils.h>
#include <mol/molecule.h>

#include <stdio.h>

void draw_main_menu(platform::Window* window);



int main(int, char**)
{


    platform::initialize();
    platform::Window* window = platform::create_window(1024, 768, "VIAMD");

    // Setup style
    ImGui::StyleColorsClassic();

    bool show_demo_window = true;
    bool show_another_window = false;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    // Main loop
    while (!(platform::window_should_close(window) && !platform::get_input_state()->key_hit[Key::KEY_ESCAPE])) {
        platform::update();


		// MAIN MENU BAR
		draw_main_menu(window);

        // 1. Show a simple window.
        // Tip: if we don't call ImGui::Begin()/ImGui::End() the widgets automatically appears in a window called "Debug".
        {
            static float f = 0.0f;
            ImGui::Text("Hello, world!");                           // Some text (you can use a format string too)
            ImGui::SliderFloat("float", &f, 0.0f, 1.0f);            // Edit 1 float as a slider from 0.0f to 1.0f
            ImGui::ColorEdit3("clear color", (float*)&clear_color); // Edit 3 floats as a color
            if (ImGui::Button("Demo Window"))                       // Use buttons to toggle our bools. We could use Checkbox() as well.
                show_demo_window ^= 1;
            if (ImGui::Button("Another Window")) 
                show_another_window ^= 1;
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
        }

        // 2. Show another simple window. In most cases you will use an explicit Begin/End pair to name the window.
        if (show_another_window)
        {
            ImGui::Begin("Another Window", &show_another_window);
            ImGui::Text("Hello from another window!");
            ImGui::End();
        }

        // 3. Show the ImGui demo window. Most of the sample code is in ImGui::ShowDemoWindow().
        if (show_demo_window)
        {
            ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiCond_FirstUseEver); // Normally user code doesn't need/want to call this because positions are saved in .ini file anyway. Here we just want to make the demo initial state a bit more friendly!
            ImGui::ShowDemoWindow(&show_demo_window);
        }

        // Rendering
        int display_w, display_h;
        platform::get_framebuffer_size(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui::Render();
        platform::swap_buffers(window);
    }

    platform::shutdown();

    return 0;
}






void draw_main_menu(platform::Window* window) {
	if (ImGui::BeginMainMenuBar())
	{
		if (ImGui::BeginMenu("File"))
		{
			if (ImGui::MenuItem("New")) {
				ImGui::OpenPopup("Are you sure?");
			}

			bool open;
			if (ImGui::BeginPopupModal("Are you sure?", &open))
			{
				ImGui::Text("Hello from Stacked The First\nUsing style.Colors[ImGuiCol_ModalWindowDarkening] for darkening.");
				static int item = 1;
				ImGui::Combo("Combo", &item, "aaaa\0bbbb\0cccc\0dddd\0eeee\0\0");
				static float color[4] = { 0.4f,0.7f,0.0f,0.5f };
				ImGui::ColorEdit4("color", color);  // This is to test behavior of stacked regular popups over a modal

				if (ImGui::Button("Add another modal.."))
					ImGui::OpenPopup("Stacked 2");
				if (ImGui::BeginPopupModal("Stacked 2"))
				{
					ImGui::Text("Hello from Stacked The Second!");
					if (ImGui::Button("Close"))
						ImGui::CloseCurrentPopup();
					ImGui::EndPopup();
				}

				if (ImGui::Button("Close"))
					ImGui::CloseCurrentPopup();
				ImGui::EndPopup();
			}

			if (ImGui::MenuItem("Open", "Ctrl+O")) {}
			if (ImGui::BeginMenu("Open Recent"))
			{
				ImGui::MenuItem("fish_hat.c");
				ImGui::MenuItem("fish_hat.inl");
				ImGui::MenuItem("fish_hat.h");
				ImGui::EndMenu();
			}
			if (ImGui::MenuItem("Save", "Ctrl+S")) {}
			if (ImGui::MenuItem("Save As..")) {}
			ImGui::Separator();
			if (ImGui::MenuItem("Quit", "Alt+F4")) {
				platform::set_window_should_close(window, true);
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
}