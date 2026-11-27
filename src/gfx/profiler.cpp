#include <gfx/profiler.h>

#include <core/md_common.h>
#include <core/md_gpu.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_writer.h>

#include <imgui.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

namespace profiler {

static constexpr int MAX_NODES      = 1024;     // sections, by their path; beyond: not timed
static constexpr int MAX_DEPTH      = 64;       // sections open at once; deeper: not timed
static constexpr int PUBLISH_FRAMES = 30;       // frames a published average is over
static constexpr int NAME_LEN       = 48;
static constexpr int MAX_STREAMS    = 4;        // streams timed on the GPU

// A section, by its path (its name and its parent's path): what the times are gathered into
struct Node {
    uint64_t key = 0;
    int      parent = -1;
    int      depth = 0;
    int      first_child = -1;     // the children in the order they were first seen
    int      last_child = -1;
    int      next_sibling = -1;
    int      seen_frame = -1;      // the last frame it was pushed in
    char     name[NAME_LEN] = {};

    bool     cpu = false;          // pushed as a section, not only made by the GPU's timings
    int      gpu_frame = -1;       // the last frame a GPU timing landed in it

    // Over the frames of the current window
    double acc_cpu_ns = 0.0;       // from its pushes to its pops
    double acc_calls = 0.0;
    double acc_gpu_s = 0.0;        // the GPU timings that landed in it, itself (not those within it)
    double acc_gpu_calls = 0.0;

    // Published, per frame
    float cpu_ms = 0.0f;
    float gpu_ms = 0.0f;           // of it and all within it; 0 where nothing is timed
    float calls = 0.0f;
    bool  gpu = false;             // GPU time within it, in the published window
    bool  shown = false;           // pushed, or timed on the GPU, within the published window
};

// What a push left on the stack: the section it timed, or none
struct Entry {
    int       node;                // -1: not timed
    md_tick_t cpu_start;
};

static struct {
    Node   nodes[MAX_NODES];
    int    node_count = 0;
    int    table[2 * MAX_NODES];   // open addressing: node index + 1 by key, 0 empty

    Entry  stack[MAX_DEPTH];
    int    depth = 0;              // entries on the stack
    int    overflow = 0;           // pushes beyond MAX_DEPTH, not timed

    bool   in_frame = false;
    int    frame = 0;
    int    window_start = 0;       // the frame the current window began in
    int    cpu_frames = 0;         // frames in the current window
    int    published_cpu_frames = 0;   // the same of the window published
    bool   warned = false;

    md_gpu_stream_t streams[MAX_STREAMS] = {};
    int    stream_count = 0;
    bool   gpu_on = false;
    md_gpu_timing_t timings[512];
} g;

static const char* ROOT_NAME = "Frame";

static uint64_t hash_name(uint64_t seed, const char* name) {
    uint64_t h = seed ^ 0xcbf29ce484222325ULL;
    for (const char* c = name; *c; ++c) {
        h ^= (uint8_t)*c;
        h *= 0x100000001b3ULL;
    }
    return h ? h : 1;
}

static void warn_once(const char* what) {
    if (!g.warned) {
        MD_LOG_DEBUG("profiler: %s", what);
        g.warned = true;
    }
}

// The section named so under 'parent' (-1: the root), made the first time it is seen; -1 when there is no room
static int node_get(int parent, const char* name) {
    const uint64_t seed = parent >= 0 ? g.nodes[parent].key : 0;
    const uint64_t key  = hash_name(seed, name);
    const int mask = (int)MD_COUNTOF(g.table) - 1;
    int slot = (int)(key & (uint64_t)mask);
    while (g.table[slot]) {
        const int i = g.table[slot] - 1;
        if (g.nodes[i].key == key) return i;
        slot = (slot + 1) & mask;
    }
    if (g.node_count >= MAX_NODES) {
        warn_once("more sections than it keeps, the rest are not timed");
        return -1;
    }
    const int i = g.node_count++;
    Node& n = g.nodes[i];
    n = Node{};
    n.key = key;
    n.parent = parent;
    n.depth = parent >= 0 ? g.nodes[parent].depth + 1 : 0;
    strncpy(n.name, name, NAME_LEN - 1);
    n.name[NAME_LEN - 1] = '\0';
    if (parent >= 0) {
        Node& p = g.nodes[parent];
        if (p.last_child >= 0) g.nodes[p.last_child].next_sibling = i;
        else                   p.first_child = i;
        p.last_child = i;
    }
    g.table[slot] = i + 1;
    return i;
}

static int stack_top_node() {
    for (int i = g.depth - 1; i >= 0; --i) {
        if (g.stack[i].node >= 0) return g.stack[i].node;
    }
    return -1;
}

static void publish() {
    const double cpu_frames = (double)MD_MAX(g.cpu_frames, 1);
    for (int i = 0; i < g.node_count; ++i) {
        Node& n = g.nodes[i];
        n.gpu    = n.gpu_frame >= g.window_start;
        n.shown  = n.seen_frame >= g.window_start || n.gpu;
        n.cpu_ms = (float)(n.acc_cpu_ns / cpu_frames * 1.0e-6);
        n.gpu_ms = (float)(n.acc_gpu_s / cpu_frames * 1.0e3);
        n.calls  = (float)((n.cpu ? n.acc_calls : n.acc_gpu_calls) / cpu_frames);
        n.acc_cpu_ns = n.acc_calls = n.acc_gpu_s = n.acc_gpu_calls = 0.0;
    }
    // A section's GPU time: all within it. A node is made after its parent, so later nodes are added up first.
    for (int i = g.node_count - 1; i > 0; --i) {
        const Node& n = g.nodes[i];
        if (n.parent >= 0) g.nodes[n.parent].gpu_ms += n.gpu_ms;
    }
    g.published_cpu_frames = g.cpu_frames;
    g.window_start = g.frame;
    g.cpu_frames = 0;
}

// GPU

// The section the streams' next operations are in: tag node + 1, 0 none
static void gpu_tag(int node) {
    if (!g.gpu_on) return;
    for (int i = 0; i < g.stream_count; ++i) md_gpu_stream_profile_tag(g.streams[i], node >= 0 ? (uint32_t)(node + 1) : 0u);
}

// A timing into the node its name makes under its section: a dotted name a path of nodes
static void gpu_land(const md_gpu_timing_t& t) {
    int node = (int)t.tag - 1;
    if (node < 0 || node >= g.node_count) node = 0;                 // outside every section: the frame's
    char part[NAME_LEN];
    for (const char* c = t.name; *c && node >= 0;) {
        const char* dot = strchr(c, '.');
        const size_t n = MD_MIN(dot ? (size_t)(dot - c) : strlen(c), (size_t)NAME_LEN - 1);
        memcpy(part, c, n);
        part[n] = '\0';
        if (n) node = node_get(node, part);
        c += n + (dot ? 1 : 0);
    }
    if (node < 0) return;
    Node& n = g.nodes[node];
    n.acc_gpu_s += t.seconds;
    n.acc_gpu_calls += 1.0;
    for (int p = node; p >= 0; p = g.nodes[p].parent) g.nodes[p].gpu_frame = g.frame;
}

static void gpu_take() {
    for (int i = 0; i < g.stream_count; ++i) {
        size_t n;
        while ((n = md_gpu_stream_profile_take(g.streams[i], g.timings, MD_COUNTOF(g.timings))) > 0) {
            for (size_t k = 0; k < n; ++k) gpu_land(g.timings[k]);
        }
    }
}

void set_gpu_streams(md_gpu_stream* const* streams, int count) {
    const bool on = g.gpu_on;
    set_gpu(false);
    g.stream_count = 0;
    for (int i = 0; i < count && g.stream_count < MAX_STREAMS; ++i) {
        if (streams[i]) g.streams[g.stream_count++] = streams[i];
    }
    set_gpu(on);
}

void set_gpu(bool on) {
    if (on == g.gpu_on) return;
    for (int i = 0; i < g.stream_count; ++i) md_gpu_stream_profile(g.streams[i], on);
    if (!on) gpu_take();                                         // what is done; the rest waits for the next time
    g.gpu_on = on;
    gpu_tag(stack_top_node());
}

void frame_begin() {
    if (g.in_frame) frame_end();
    if (g.gpu_on) gpu_take();

    if (g.cpu_frames >= PUBLISH_FRAMES) publish();

    g.in_frame = true;
    g.depth = 0;
    g.overflow = 0;
    push(ROOT_NAME);
}

void frame_end() {
    if (!g.in_frame) return;
    if (g.depth != 1 || g.overflow) {
        warn_once("a section was pushed and not popped within the frame");
        while (g.depth > 1 || g.overflow) pop();
    }
    pop();
    g.in_frame = false;
    g.cpu_frames += 1;
    g.frame += 1;
}

void shutdown() {}

void push(const char* name) {
    if (!g.in_frame) return;
    if (g.depth >= MAX_DEPTH) {
        g.overflow += 1;
        return;
    }
    const int parent = stack_top_node();
    const int node = (g.depth == 0 || parent >= 0) ? node_get(parent, name) : -1;
    g.stack[g.depth++] = {node, md_tick_now()};
    if (node < 0) return;

    g.nodes[node].seen_frame = g.frame;
    g.nodes[node].acc_calls += 1.0;
    g.nodes[node].cpu = true;
    gpu_tag(node);
}

void pop() {
    if (!g.in_frame) return;
    if (g.overflow) {
        g.overflow -= 1;
        return;
    }
    if (g.depth == 0) {
        warn_once("a section was popped that was not pushed");
        return;
    }
    const Entry e = g.stack[--g.depth];
    if (e.node < 0) return;

    g.nodes[e.node].acc_cpu_ns += md_tick_to_nanoseconds(md_tick_now() - e.cpu_start);
    gpu_tag(stack_top_node());
}

// RESULTS

static bool within_named(const Node& n, const char* within) {
    if (!within) return true;
    for (int p = n.parent; p >= 0; p = g.nodes[p].parent) {
        if (strcmp(g.nodes[p].name, within) == 0) return true;
    }
    return false;
}

float cpu_ms(const char* name, const char* within) {
    float ms = 0.0f;
    for (int i = 0; i < g.node_count; ++i) {
        const Node& n = g.nodes[i];
        if (n.shown && strcmp(n.name, name) == 0 && within_named(n, within)) ms += n.cpu_ms;
    }
    return ms;
}

float frame_cpu_ms() { return g.node_count ? g.nodes[0].cpu_ms : 0.0f; }

float gpu_ms(const char* name, const char* within) {
    float ms = 0.0f;
    for (int i = 0; i < g.node_count; ++i) {
        const Node& n = g.nodes[i];
        if (n.gpu && strcmp(n.name, name) == 0 && within_named(n, within)) ms += n.gpu_ms;
    }
    return ms;
}

float frame_gpu_ms() { return g.node_count && g.nodes[0].gpu ? g.nodes[0].gpu_ms : 0.0f; }

// A bar under the number: the share of the frame's time
static void time_cell(float ms, float frame_ms) {
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const float  w = ImGui::GetContentRegionAvail().x;
    const float  h = ImGui::GetTextLineHeight();
    if (frame_ms > 0.0f && ms > 0.0f) {
        const float f = MD_MIN(ms / frame_ms, 1.0f);
        ImVec4 c = ImGui::GetStyleColorVec4(ImGuiCol_PlotHistogram);
        c.w *= 0.35f;
        ImGui::GetWindowDrawList()->AddRectFilled(p, ImVec2(p.x + w * f, p.y + h), ImGui::GetColorU32(c));
    }
    ImGui::Text("%7.3f", ms);
}

static bool has_shown_child(const Node& n) {
    for (int c = n.first_child; c >= 0; c = g.nodes[c].next_sibling) {
        if (g.nodes[c].shown) return true;
    }
    return false;
}

static void draw_row(int i, float frame_ms) {
    const Node& n = g.nodes[i];
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::PushID(i);
    ImGuiTreeNodeFlags flags = ImGuiTreeNodeFlags_SpanAllColumns;
    const bool leaf = !has_shown_child(n);
    if (leaf)          flags |= ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;
    if (n.depth < 3)   flags |= ImGuiTreeNodeFlags_DefaultOpen;
    if (!n.cpu) ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));   // GPU work only
    const bool open = ImGui::TreeNodeEx("##node", flags, "%s", n.name);
    if (!n.cpu) ImGui::PopStyleColor();
    ImGui::TableNextColumn();
    if (n.cpu) time_cell(n.cpu_ms, frame_ms);
    ImGui::TableNextColumn();
    if (n.gpu) time_cell(n.gpu_ms, frame_gpu_ms());
    ImGui::TableNextColumn();
    if (n.calls != 1.0f && n.calls > 0.0f) ImGui::TextDisabled(n.calls < 10.0f && n.calls != (float)(int)n.calls ? "%.1f" : "%.0f", n.calls);
    ImGui::PopID();
    if (open && !leaf) {
        for (int c = n.first_child; c >= 0; c = g.nodes[c].next_sibling) {
            if (g.nodes[c].shown) draw_row(c, frame_ms);
        }
        ImGui::TreePop();
    }
}

void draw_table() {
    if (!g.node_count || !g.nodes[0].shown) {
        ImGui::TextDisabled("Measuring...");
        return;
    }
    const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_Resizable |
                                  ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit;
    if (ImGui::BeginTable("##profiler", 4, flags)) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("Section", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("CPU ms", ImGuiTableColumnFlags_WidthFixed, ImGui::CalcTextSize("0000.000000").x);
        ImGui::TableSetupColumn("GPU ms", ImGuiTableColumnFlags_WidthFixed, ImGui::CalcTextSize("0000.000000").x);
        ImGui::TableSetupColumn("Calls",  ImGuiTableColumnFlags_WidthFixed, ImGui::CalcTextSize("00000").x);
        ImGui::TableHeadersRow();
        draw_row(0, g.nodes[0].cpu_ms);
        ImGui::EndTable();
    }
}

void draw_summary() {
    if (!g.node_count || !g.nodes[0].shown) {
        ImGui::TextDisabled("Measuring...");
        return;
    }
    const Node& root = g.nodes[0];
    if (root.gpu) ImGui::Text("Per frame: CPU %.2f ms, GPU %.2f ms", root.cpu_ms, root.gpu_ms);
    else          ImGui::Text("Per frame: CPU %.2f ms", root.cpu_ms);
    for (int c = root.first_child; c >= 0; c = g.nodes[c].next_sibling) {
        const Node& n = g.nodes[c];
        if (!n.shown) continue;
        if (n.gpu) ImGui::Text("  %-22s CPU %7.3f ms  GPU %7.3f ms", n.name, n.cpu_ms, n.gpu_ms);
        else       ImGui::Text("  %-22s CPU %7.3f ms", n.name, n.cpu_ms);
    }
    ImGui::TextDisabled("Every section: Windows > Profiler");
}

// EXPORT

// A field of a CSV row, quoted when it has to be
static void write_csv_field(md_writer_i* w, const char* s) {
    if (!strpbrk(s, ",\"\r\n")) {
        md_write_cstr(w, s);
        return;
    }
    md_write_char(w, '"');
    for (const char* c = s; *c; ++c) {
        if (*c == '"') md_write_char(w, '"');
        md_write_char(w, *c);
    }
    md_write_char(w, '"');
}

// The section and those shown within it, in the order of the table
static void write_csv_rows(md_writer_i* w, int i) {
    const Node& n = g.nodes[i];
    // Its path: the names from the root, written before the rows within it are (so one buffer does for all of them)
    static char path[MAX_DEPTH * NAME_LEN];
    int chain[MAX_DEPTH + 1];
    int count = 0;
    for (int p = i; p >= 0 && count < (int)MD_COUNTOF(chain); p = g.nodes[p].parent) chain[count++] = p;
    size_t len = 0;
    path[0] = '\0';
    for (int k = count - 1; k >= 0 && len + 1 < sizeof(path); --k) {
        const int m = snprintf(path + len, sizeof(path) - len, k == count - 1 ? "%s" : "/%s", g.nodes[chain[k]].name);
        if (m < 0) break;
        len = MD_MIN(len + (size_t)m, sizeof(path) - 1);
    }
    write_csv_field(w, path);
    if (n.cpu) md_writef(w, ",%.4f", n.cpu_ms); else md_write_char(w, ',');
    if (n.gpu) md_writef(w, ",%.4f", n.gpu_ms); else md_write_char(w, ',');
    md_writef(w, ",%.2f\n", n.calls);
    for (int c = n.first_child; c >= 0; c = g.nodes[c].next_sibling) {
        if (g.nodes[c].shown) write_csv_rows(w, c);
    }
}

void write_csv(md_writer_i* w, const char* info) {
    char when[64] = "";
    const time_t t = time(nullptr);
    if (const struct tm* tm = localtime(&t)) strftime(when, sizeof(when), "%Y-%m-%d %H:%M:%S", tm);

    md_writef(w, "# VIAMD profiler, %s\n", when);
#ifdef NDEBUG
    md_write_cstr(w, "# build: NDEBUG\n");
#else
    md_write_cstr(w, "# build: debug (asserts on)\n");
#endif
    md_writef(w, "# per frame, averaged over %d frames\n", g.published_cpu_frames);
    for (const char* line = info; line && *line;) {
        const char* end = strchr(line, '\n');
        const size_t n = end ? (size_t)(end - line) : strlen(line);
        md_write_cstr(w, "# ");
        md_write(w, line, n);
        md_write_char(w, '\n');
        line = end ? end + 1 : line + n;
    }
    md_writef(w, "# GPU times %s\n", g.gpu_on ? "on (each operation's, by its kernel or pass)" : "off");
    md_write_cstr(w, "section,cpu_ms,gpu_ms,calls\n");
    if (g.node_count && g.nodes[0].shown) write_csv_rows(w, 0);
}

}  // namespace profiler
