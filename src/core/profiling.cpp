#include "profiling.h"
#include <core/hash.h>
#include <core/platform.h>
#include <core/gl.h>
#include <imgui.h>

namespace cpu_profiling {

typedef unsigned int Key;

struct Section {
    Key key = 0;
    float32 total_ms = 0;
    platform::Timestamp t0 = 0;

    int32 parent_idx = -1;
    StringBuffer<32> name{};
};

constexpr int MAX_SECTIONS = 1024;
static Section sections[1024];
static int section_count = 0;
static int current_section_idx = -1;

void initialize() {}

void shutdown() {}

int32 find_section_idx(Key key) {
    for (int32 i = 0; i < section_count; i++) {
        if (key == sections[i].key) return i;
    }
    return -1;
}

void push_section(CString sec) {
    Key key = hash::crc32(sec);
    int idx = find_section_idx(key);
    if (idx == -1) {
        idx = section_count++;
        sections[idx] = {key, 0, 0, current_section_idx, sec};
    }
    current_section_idx = idx;
    sections[current_section_idx].t0 = platform::get_time();
}

void pop_section(CString sec) {
#ifdef DEBUG
    ASSERT(current_section_idx != -1);
    Key key = hash::crc32(sec);
    int idx = find_section_idx(key);
    ASSERT(current_section_idx == idx);
#else
    (void)sec;
#endif
    auto& section = sections[current_section_idx];
    section.total_ms += platform::compute_delta_ms(section.t0, platform::get_time());
    current_section_idx = sections[current_section_idx].parent_idx;
}

void clear() {
    section_count = 0;
    current_section_idx = -1;
}
void finish() {}

}  // namespace cpu_profiling

namespace gpu_profiling {

typedef unsigned int Key;

struct Section {
    Key key = 0;
	float beg_ms = 0;
	float end_ms = 0;
    unsigned int query_ids[2] = {0, 0};

    int parent_idx = -1;
    StringBuffer<32> name{};
};

constexpr int MAX_SECTIONS = 1024;
static Section sections[MAX_SECTIONS];
static int section_count = 0;
static int current_section_idx = -1;

void initialize() {}
void shutdown() {
	for (int i = 0; i < section_count; i++) {
		glDeleteQueries(2, sections[i].query_ids);
	}
}

int32 find_section_idx(Key key) {
    for (int32 i = 0; i < section_count; i++) {
        if (key == sections[i].key) return i;
    }
    return -1;
}

void push_section(CString sec) {
    Key key = hash::crc32(sec);
    int idx = find_section_idx(key);
    if (idx == -1) {
        idx = section_count++;
        sections[idx] = {key, 0, 0, {0, 0}, current_section_idx, sec};
        glGenQueries(2, sections[idx].query_ids);
    }
    current_section_idx = idx;
    glQueryCounter(sections[idx].query_ids[0], GL_TIMESTAMP);
}

void pop_section(CString sec) {
#ifdef DEBUG
    ASSERT(current_section_idx != -1);
    Key key = hash::crc32(sec);
    int idx = find_section_idx(key);
    ASSERT(current_section_idx == idx);
#else
    (void)sec;
#endif
    auto& section = sections[current_section_idx];
    glQueryCounter(section.query_ids[1], GL_TIMESTAMP);
    current_section_idx = section.parent_idx;
}

void clear() {
	//for (int i = 0; i < section_count; i++) {
	//	sections[i].beg_ms = 0;
	//	sections[i].end_ms = 0;
	//}
}

void finish() {
	static unsigned int frame_count = 0;
	frame_count++;
    for (int i = 0; i < section_count; i++) {
        int done = 0;
        while (!done) {
            glGetQueryObjectiv(sections[i].query_ids[1], GL_QUERY_RESULT_AVAILABLE, &done);
        }
        GLuint64 t0, t1;
        glGetQueryObjectui64v(sections[i].query_ids[0], GL_QUERY_RESULT, &t0);
        glGetQueryObjectui64v(sections[i].query_ids[1], GL_QUERY_RESULT, &t1);
		sections[i].beg_ms = t0 / 1000000.f;
		sections[i].end_ms = t1 / 1000000.f;
    }
}

int get_depth(int section_idx) {
    int count = 0;
    Section* sec = &sections[section_idx];
    while (sec->parent_idx != -1) {
        count++;
        sec = &sections[sec->parent_idx];
    }
    return count;
}

void print() {
    StringBuffer<32> buff;
    for (int i = 0; i < section_count; i++) {
        int count = 0;
        count += snprintf(buff.beg(), 32, "%*c [%s]:", get_depth(i) * 4, ' ', sections[i].name.cstr());
        printf("%s%8.4f\n", buff.cstr(), (sections[i].end_ms - sections[i].beg_ms));
    }
    printf("\n");
}

void draw_window() {
	ImGui::SetNextWindowSize(ImVec2(520, 600), ImGuiSetCond_FirstUseEver);
	if (!ImGui::Begin("GPU Profiler"))
	{
		ImGui::End();
		return;
	}

	if (section_count == 0) {
		ImGui::End();
		return;
	}

	static int current_section_idx = 1;
	char tmps[1024];
	int len = 0;
	for (int i = 0; i < section_count; i++) {
		len += snprintf(tmps + len, 1024 - len, "%s\0", sections[i].name.cstr()) + 1;
	}
	snprintf(tmps + len, 1024 - len, "\0");
	static int timeRange = 0;
	//ImGui::PushItemWidth(-100);

	/*
	if (ImGui::Button("Copy output to clipboard"))
	{
		ImGui::LogToClipboard();
		YAP::LogDump((int(*)(const char *szFormat, ...))&ImGui::LogText);
		ImGui::LogFinish();
	}
	ImGui::SameLine();
	if (ImGui::Button(gProfilerPaused ? "Resume" : "Pause"))
		gProfilerPaused = !gProfilerPaused;
	*/

	//ImGui::PopItemWidth();

	const float total_min_time = sections[0].beg_ms;
	const float total_max_time = sections[0].end_ms;
	const float width = ImGui::GetWindowContentRegionWidth();
	const float scale = width / (total_max_time - total_min_time);
	const float height = 20.f;
	const float offset = 20.f;

	// GUI
	ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 1)); // Tighten spacing
	for (int section_idx = 0; section_idx < section_count; section_idx++)
	{
		int depth = get_depth(section_idx);

		ImGui::SetCursorPos(ImVec2((sections[section_idx].beg_ms - total_min_time) * scale, depth * height + offset));
		//ImGui::SetCursorPosY();
		//ImGui::BeginGroup();
		
		/*
		float preTime = (float)(sections[section_idx].beg_ms - total_min_time);
		float invisibleButWidth = preTime * scale;
		if (invisibleButWidth >= 1.f)
		{
			ImGui::InvisibleButton("", ImVec2(invisibleButWidth, height));
			ImGui::SameLine();
		}
		*/
		float duration = (float)(sections[section_idx].end_ms - sections[section_idx].beg_ms);
		ImVec4 buttonColor(0.1,0.6,0.1,1);
		ImGui::PushStyleColor(ImGuiCol_Button, buttonColor);
		ImGui::PushStyleColor(ImGuiCol_ButtonHovered, buttonColor);
		ImGui::PushStyleColor(ImGuiCol_ButtonActive, buttonColor);
		float butWidth = duration * scale;
		if (butWidth >= 1.f)
		{
			ImGui::Button(sections[section_idx].name.cstr(), ImVec2(butWidth, height));
			if (ImGui::IsItemHovered())
			{
				ImGui::SetTooltip("%s\n\nStart %5.3f ms\nEnd %5.3f ms\nDuration %5.3f ms",
					sections[section_idx].name.cstr(),
					(float)(sections[section_idx].beg_ms - total_min_time),
					(float)(sections[section_idx].end_ms - sections[section_idx].beg_ms),
					duration
				);
			}
		}
		ImGui::PopStyleColor(3);
		//ImGui::SameLine();
		//ImGui::EndGroup();
	}
	ImGui::PopStyleVar();
	ImGui::End();
}

};  // namespace gpu_profiling
