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

constexpr int NUM_SAMPLES = 256;
constexpr int MAX_SECTIONS = 1024;

constexpr int NUM_COLORS = 13;
constexpr uint32 color_table[NUM_COLORS] = {0xffffffff, 0xffc7d38d, 0xffb3ffff, 0xffdababe, 0xff7280fb, 0xffd3b180, 0xff62b4fd,
                                            0xff69deb3, 0xffe5cdfc, 0xffd9d9d9, 0xffbd80bc, 0xffc5ebcc, 0xff6fedff};

struct Section {
    unsigned int query_ids[2] = {0, 0};
    int parent_idx = -1;
    StringBuffer<32> name{};
    float beg_avg = 0;
    float end_avg = 0;
    float beg_ms[NUM_SAMPLES] = {0};
    float end_ms[NUM_SAMPLES] = {0};
};

static Section sections[MAX_SECTIONS];
static int section_count = 0;
static int current_section_idx = -1;

void initialize() {
    for (int i = 0; i < MAX_SECTIONS; i++) {
        sections[i] = {};
        glGenQueries(2, sections[i].query_ids);
    }
}
void shutdown() {
    for (int i = 0; i < MAX_SECTIONS; i++) {
        glDeleteQueries(2, sections[i].query_ids);
    }
}

void push_section(CString sec) {
    int i = section_count++;
    sections[i].parent_idx = current_section_idx;
    sections[i].name = sec;
    current_section_idx = i;
    glQueryCounter(sections[i].query_ids[0], GL_TIMESTAMP);
}

void pop_section() {
    ASSERT(current_section_idx != -1);
    glQueryCounter(sections[current_section_idx].query_ids[1], GL_TIMESTAMP);
    current_section_idx = sections[current_section_idx].parent_idx;
}

void begin_frame() {
    section_count = 0;
    push_section("FRAME");
}

void end_frame() {
    // glFlush();
    // glFinish();
    pop_section();  // "FRAME";
    static unsigned int frame_count = 0;
    frame_count++;
    int sample = frame_count % NUM_SAMPLES;
    for (int i = 0; i < section_count; i++) {
        int done = 0;
        while (!done) {
            glGetQueryObjectiv(sections[i].query_ids[1], GL_QUERY_RESULT_AVAILABLE, &done);
        }
        GLuint64 t0, t1;
        glGetQueryObjectui64v(sections[i].query_ids[0], GL_QUERY_RESULT, &t0);
        glGetQueryObjectui64v(sections[i].query_ids[1], GL_QUERY_RESULT, &t1);
        sections[i].beg_ms[sample] = t0 / 1000000.f;
        sections[i].end_ms[sample] = t1 / 1000000.f;

        sections[i].beg_avg = 0;
        sections[i].end_avg = 0;
        for (int u = 0; u < NUM_SAMPLES; u++) {
            sections[i].beg_avg += sections[i].beg_ms[u];
            sections[i].end_avg += sections[i].end_ms[u];
        }
        sections[i].beg_avg /= (float)NUM_SAMPLES;
        sections[i].end_avg /= (float)NUM_SAMPLES;
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

void draw_window(bool* window_open) {
    ImGui::SetNextWindowSize(ImVec2(600, 200), ImGuiSetCond_FirstUseEver);
    if (!ImGui::Begin("GPU Profiler", window_open, ImGuiWindowFlags_NoFocusOnAppearing)) {
        ImGui::End();
        return;
    }

    if (section_count == 0) {
        ImGui::End();
        return;
    }

    const float total_min_time = sections[0].beg_avg;
    const float total_max_time = sections[0].end_avg;
    const float width = ImGui::GetWindowContentRegionWidth();
    const float scale = width / (total_max_time - total_min_time);
    const float height = 20.f;
    const float offset = 20.f;

    // GUI
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 1));  // Tighten spacing
    for (int section_idx = 0; section_idx < section_count; section_idx++) {
        int depth = get_depth(section_idx);

        ImGui::SetCursorPos(ImVec2((sections[section_idx].beg_avg - total_min_time) * scale, depth * height + offset));

        float duration = (float)(sections[section_idx].end_avg - sections[section_idx].beg_avg);
        float butWidth = duration * scale;
        if (butWidth >= 1.f) {
            uint32 c = color_table[section_idx % NUM_COLORS];
            const float scl = 1.f / 255.f;
            ImVec4 color = {(c >> 0 & 0xff) * scl, (c >> 8 & 0xff) * scl, (c >> 16 & 0xff) * scl, (c >> 24 & 0xff) * scl};
            ImGui::PushStyleColor(ImGuiCol_Text, {0, 0, 0, 1});
            ImGui::PushStyleColor(ImGuiCol_Button, color);
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, color);
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, color);
            ImGui::Button(sections[section_idx].name.cstr(), ImVec2(butWidth, height));
            ImGui::PopStyleColor(4);

            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("%s\n\nStart %5.3f ms\nEnd %5.3f ms\nDuration %5.3f ms", sections[section_idx].name.cstr(),
                                  (float)(sections[section_idx].beg_avg - total_min_time),
                                  (float)(sections[section_idx].end_avg - sections[section_idx].beg_avg), duration);
            }
        }
    }
    ImGui::PopStyleVar();
    ImGui::End();
}
};  // namespace gpu_profiling
