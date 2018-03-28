#pragma once

#include <core/common.h>
#include <imgui.h>

struct Console
{
	char                  input_buffer[256];
	ImVector<char*>       items;
	bool                  scroll_to_bottom;
	ImVector<char*>       history;
	int                   history_pos;    // -1: new line, 0..History.Size-1 browsing history.
	ImVector<const char*> commands;
	float				  y_pos = -10000;
	bool				  visible = false;

	Console()
	{
		ClearLog();
		memset(input_buffer, 0, sizeof(input_buffer));
		history_pos = -1;
		commands.push_back("HELP");
		commands.push_back("HISTORY");
		commands.push_back("CLEAR");
		commands.push_back("CLASSIFY");  // "classify" is here to provide an example of "C"+[tab] completing to "CL" and displaying matches.
		AddLog("Welcome to ImGui!");
	}
	~Console()
	{
		ClearLog();
		for (int i = 0; i < history.Size; i++)
			free(history[i]);
	}

	// Portable helpers
	static int   Stricmp(const char* str1, const char* str2) { int d; while ((d = toupper(*str2) - toupper(*str1)) == 0 && *str1) { str1++; str2++; } return d; }
	static int   Strnicmp(const char* str1, const char* str2, int n) { int d = 0; while (n > 0 && (d = toupper(*str2) - toupper(*str1)) == 0 && *str1) { str1++; str2++; n--; } return d; }
	static char* Strdup(const char *str) { size_t len = strlen(str) + 1; void* buff = malloc(len); return (char*)memcpy(buff, (const void*)str, len); }

	void ClearLog()
	{
		for (int i = 0; i < items.Size; i++)
			free(items[i]);
		items.clear();
		scroll_to_bottom = true;
	}

	void AddLog(const char* fmt, ...) IM_FMTARGS(2)
	{
		// FIXME-OPT
		char buf[1024];
		va_list args;
		va_start(args, fmt);
		vsnprintf(buf, IM_ARRAYSIZE(buf), fmt, args);
		buf[IM_ARRAYSIZE(buf) - 1] = 0;
		va_end(args);
		items.push_back(Strdup(buf));
		scroll_to_bottom = true;
	}

	void Draw(const char* title, int width, int height, float dt)
	{
		constexpr int WINDOW_FLAGS = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoScrollbar |
			ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_Modal;

		float console_width = (float)width;
		float console_height = (float)height * 0.75f;

		// Quarter of a second to hide/show
		float speed = console_height * dt * 4.f;

		float target_hide_y = -console_height;
		float target_show_y = 0;
		float target_y = visible ? target_show_y : target_hide_y;

		if (y_pos != target_y) {
			float delta = target_y < y_pos ? -speed : speed;
			y_pos += delta;
			if (y_pos < target_hide_y) y_pos = target_hide_y;
			if (y_pos > target_show_y) y_pos = target_show_y;
		}

		bool console_fully_shown = (y_pos == target_show_y);
		bool console_fully_hidden = (y_pos == target_hide_y);

		if (console_fully_hidden) return;

		ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.f);

		ImGui::SetNextWindowSize(ImVec2(console_width, console_height));
		ImGui::SetNextWindowPos(ImVec2(0.f, y_pos));
		ImGui::Begin(title, &visible, WINDOW_FLAGS);

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
		ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(4, 1)); // Tighten spacing
																	  //if (copy_to_clipboard)
																	  //    ImGui::LogToClipboard();
		ImVec4 col_default_text = ImGui::GetStyleColorVec4(ImGuiCol_Text);
		for (int i = 0; i < items.Size; i++) {
			const char* item = items[i];
			//if (!filter.PassFilter(item))
			//    continue;
			ImVec4 col = col_default_text;
			if (strstr(item, "[error]")) col = ImColor(1.0f, 0.4f, 0.4f, 1.0f);
			else if (strncmp(item, "# ", 2) == 0) col = ImColor(1.0f, 0.78f, 0.58f, 1.0f);
			ImGui::PushStyleColor(ImGuiCol_Text, col);
			ImGui::TextUnformatted(item);
			ImGui::PopStyleColor();
		}
		//if (copy_to_clipboard)
		//    ImGui::LogFinish();
		if (scroll_to_bottom)
			ImGui::SetScrollHere();
		scroll_to_bottom = false;
		ImGui::PopStyleVar();
		ImGui::EndChild();
		ImGui::Separator();

		if (console_fully_shown) {
			// Command-line
			ImGui::PushItemWidth(-1);
			if (ImGui::InputText("##Input", input_buffer, IM_ARRAYSIZE(input_buffer), ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CallbackCompletion | ImGuiInputTextFlags_CallbackHistory, &TextEditCallbackStub, (void*)this))
			{
				char* input_end = input_buffer + strlen(input_buffer);
				while (input_end > input_buffer && input_end[-1] == ' ') { input_end--; } *input_end = 0;
				if (input_buffer[0])
					ExecCommand(input_buffer);
				strcpy(input_buffer, "");
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
		history_pos = -1;
		for (int i = history.Size - 1; i >= 0; i--)
			if (Stricmp(history[i], command_line) == 0)
			{
				free(history[i]);
				history.erase(history.begin() + i);
				break;
			}
		history.push_back(Strdup(command_line));

		// Process command
		if (Stricmp(command_line, "CLEAR") == 0)
		{
			ClearLog();
		}
		else if (Stricmp(command_line, "HELP") == 0)
		{
			AddLog("Commands:");
			for (int i = 0; i < commands.Size; i++)
				AddLog("- %s", commands[i]);
		}
		else if (Stricmp(command_line, "HISTORY") == 0)
		{
			int first = history.Size - 10;
			for (int i = first > 0 ? first : 0; i < history.Size; i++)
				AddLog("%3d: %s\n", i, history[i]);
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
			for (int i = 0; i < commands.Size; i++)
				if (Strnicmp(commands[i], word_start, (int)(word_end - word_start)) == 0)
					candidates.push_back(commands[i]);

			if (candidates.Size == 0)
			{
				// No match
				AddLog("No match for \"%.*s\"!\n", (int)(word_end - word_start), word_start);
			}
			else if (candidates.Size == 1)
			{
				// Single match. Delete the beginning of the word and replace it entirely so we've got nice casing
				data->DeleteChars((int)(word_start - data->Buf), (int)(word_end - word_start));
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
					data->DeleteChars((int)(word_start - data->Buf), (int)(word_end - word_start));
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
			const int prev_history_pos = history_pos;
			if (data->EventKey == ImGuiKey_UpArrow)
			{
				if (history_pos == -1)
					history_pos = history.Size - 1;
				else if (history_pos > 0)
					history_pos--;
			}
			else if (data->EventKey == ImGuiKey_DownArrow)
			{
				if (history_pos != -1)
					if (++history_pos >= history.Size)
						history_pos = -1;
			}

			// A better implementation would preserve the data on the current input line along with cursor position.
			if (prev_history_pos != history_pos)
			{
				data->CursorPos = data->SelectionStart = data->SelectionEnd = data->BufTextLen = (int)snprintf(data->Buf, (size_t)data->BufSize, "%s", (history_pos >= 0) ? history[history_pos] : "");
				data->BufDirty = true;
			}
		}
		}
		return 0;
	}
};