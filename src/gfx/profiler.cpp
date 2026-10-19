#include <gfx/profiler.h>

#include <gfx/gl.h>

#include <core/md_common.h>
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

static constexpr int MAX_NODES      = 1024;     // sections, by their path; beyond: debug groups only
static constexpr int MAX_DEPTH      = 64;       // sections open at once; deeper: debug groups only
static constexpr int FRAMES         = 4;        // frames whose queries can be in flight; a fifth waits for the oldest
static constexpr int PUBLISH_FRAMES = 30;       // frames a published average is over
static constexpr int NAME_LEN       = 48;

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

    // Over the frames of the current window
    double acc_gpu_self_ns = 0.0;  // while it was the innermost
    double acc_cpu_ns = 0.0;       // from its pushes to its pops
    double acc_calls = 0.0;

    // Published, per frame
    float gpu_self_ms = 0.0f;
    float gpu_ms = 0.0f;           // its own and its sections'
    float cpu_ms = 0.0f;
    float calls = 0.0f;
    bool  shown = false;           // pushed within the published window
};

// The GPU time of a frame: a query per stretch of it, and the section it went to
struct Segment {
    GLuint query;
    int    node;
};

struct Record {
    GLuint*  queries = nullptr;    // owned, made as more are needed
    int      query_cap = 0;
    Segment* segments = nullptr;
    int      count = 0;
    bool     pending = false;      // its results are yet to be read
};

// What a push left on the stack: the section it timed, or none (a debug group only)
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
    int    overflow = 0;           // pushes beyond MAX_DEPTH, debug groups only

    Record records[FRAMES];
    Record* rec = nullptr;         // this frame's, or none (the oldest still in flight: no GPU timing this frame)
    bool   query_running = false;
    int    paused = 0;

    bool   in_frame = false;
    int    frame = 0;
    int    window_start = 0;       // the frame the current window began in
    int    cpu_frames = 0;         // frames in the current window
    int    gpu_frames = 0;         // frames whose GPU times went into it
    int    published_cpu_frames = 0;   // the same of the window published
    int    published_gpu_frames = 0;
    int    publish_count = 0;
    bool   warned = false;
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
    const int mask = (int)ARRAY_SIZE(g.table) - 1;
    int slot = (int)(key & (uint64_t)mask);
    while (g.table[slot]) {
        const int i = g.table[slot] - 1;
        if (g.nodes[i].key == key) return i;
        slot = (slot + 1) & mask;
    }
    if (g.node_count >= MAX_NODES) {
        warn_once("more sections than it keeps, the rest are debug groups only");
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

// GPU TIMING
// One query runs at a time, for the innermost section

static void segment_begin(int node) {
    if (!g.rec || g.paused || node < 0) return;
    Record& r = *g.rec;
    if (r.count == r.query_cap) {
        const int cap = r.query_cap ? 2 * r.query_cap : 256;
        r.queries  = (GLuint*)realloc(r.queries, sizeof(GLuint) * cap);
        r.segments = (Segment*)realloc(r.segments, sizeof(Segment) * cap);
        glGenQueries(cap - r.query_cap, r.queries + r.query_cap);
        r.query_cap = cap;
    }
    const GLuint q = r.queries[r.count];
    glBeginQuery(GL_TIME_ELAPSED, q);
    r.segments[r.count] = {q, node};
    g.query_running = true;
}

static void segment_end() {
    if (!g.query_running) return;
    glEndQuery(GL_TIME_ELAPSED);
    g.rec->count += 1;
    g.query_running = false;
}

static int stack_top_node() {
    for (int i = g.depth - 1; i >= 0; --i) {
        if (g.stack[i].node >= 0) return g.stack[i].node;
    }
    return -1;
}

// Reads the results of the frames that are done; a frame is read whole, once its last query is
static void collect() {
    for (int f = 0; f < FRAMES; ++f) {
        Record& r = g.records[f];
        if (!r.pending) continue;
        GLint available = 0;
        glGetQueryObjectiv(r.segments[r.count - 1].query, GL_QUERY_RESULT_AVAILABLE, &available);
        if (!available) continue;
        for (int i = 0; i < r.count; ++i) {
            GLuint64 ns = 0;
            glGetQueryObjectui64v(r.segments[i].query, GL_QUERY_RESULT, &ns);
            g.nodes[r.segments[i].node].acc_gpu_self_ns += (double)ns;
        }
        r.pending = false;
        g.gpu_frames += 1;
    }
}

static void publish() {
    const double cpu_frames = (double)MAX(g.cpu_frames, 1);
    const double gpu_frames = (double)g.gpu_frames;
    for (int i = 0; i < g.node_count; ++i) {
        Node& n = g.nodes[i];
        n.shown       = n.seen_frame >= g.window_start;
        n.gpu_self_ms = gpu_frames > 0 ? (float)(n.acc_gpu_self_ns / gpu_frames * 1.0e-6) : 0.0f;
        n.gpu_ms      = n.gpu_self_ms;
        n.cpu_ms      = (float)(n.acc_cpu_ns / cpu_frames * 1.0e-6);
        n.calls       = (float)(n.acc_calls / cpu_frames);
        n.acc_gpu_self_ns = n.acc_cpu_ns = n.acc_calls = 0.0;
    }
    // A child is made after its parent: from the last, every section's time is complete when it is added to its parent's
    for (int i = g.node_count - 1; i >= 0; --i) {
        const Node& n = g.nodes[i];
        if (n.parent >= 0) g.nodes[n.parent].gpu_ms += n.gpu_ms;
    }
    g.published_cpu_frames = g.cpu_frames;
    g.published_gpu_frames = g.gpu_frames;
    g.publish_count += 1;
    g.window_start = g.frame;
    g.cpu_frames = 0;
    g.gpu_frames = 0;
}

void frame_begin() {
    if (g.in_frame) frame_end();

    collect();
    if (g.cpu_frames >= PUBLISH_FRAMES) publish();

    Record& r = g.records[g.frame % FRAMES];
    g.rec = r.pending ? nullptr : &r;
    if (g.rec) g.rec->count = 0;

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
    if (g.rec) g.rec->pending = g.rec->count > 0;
    g.rec = nullptr;
    g.in_frame = false;
    g.cpu_frames += 1;
    g.frame += 1;
}

void shutdown() {
    if (g.query_running) segment_end();
    for (Record& r : g.records) {
        if (r.query_cap) glDeleteQueries(r.query_cap, r.queries);
        free(r.queries);
        free(r.segments);
        r = Record{};
    }
    g.rec = nullptr;
}

void push(const char* name) {
    if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, 0, -1, name);
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
    segment_end();
    segment_begin(node);
}

void pop() {
    if (glPopDebugGroup) glPopDebugGroup();
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
    segment_end();
    segment_begin(stack_top_node());
}

Pause::Pause() {
    if (g.paused++ == 0) segment_end();
}

Pause::~Pause() {
    if (--g.paused == 0 && g.in_frame) segment_begin(stack_top_node());
}

// RESULTS

static bool within_named(const Node& n, const char* within) {
    if (!within) return true;
    for (int p = n.parent; p >= 0; p = g.nodes[p].parent) {
        if (strcmp(g.nodes[p].name, within) == 0) return true;
    }
    return false;
}

float gpu_ms(const char* name, const char* within) {
    float ms = 0.0f;
    for (int i = 0; i < g.node_count; ++i) {
        const Node& n = g.nodes[i];
        if (n.shown && strcmp(n.name, name) == 0 && within_named(n, within)) ms += n.gpu_ms;
    }
    return ms;
}

float cpu_ms(const char* name, const char* within) {
    float ms = 0.0f;
    for (int i = 0; i < g.node_count; ++i) {
        const Node& n = g.nodes[i];
        if (n.shown && strcmp(n.name, name) == 0 && within_named(n, within)) ms += n.cpu_ms;
    }
    return ms;
}

float frame_gpu_ms() { return g.node_count ? g.nodes[0].gpu_ms : 0.0f; }
float frame_cpu_ms() { return g.node_count ? g.nodes[0].cpu_ms : 0.0f; }

// A bar under the number: the share of the frame's GPU time
static void time_cell(float ms, float frame_ms) {
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const float  w = ImGui::GetContentRegionAvail().x;
    const float  h = ImGui::GetTextLineHeight();
    if (frame_ms > 0.0f && ms > 0.0f) {
        const float f = MIN(ms / frame_ms, 1.0f);
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
    const bool open = ImGui::TreeNodeEx("##node", flags, "%s", n.name);
    if (ImGui::IsItemHovered() && n.first_child >= 0) {
        ImGui::SetItemTooltip("Its own GPU time: %.3f ms", n.gpu_self_ms);
    }
    ImGui::TableNextColumn();
    time_cell(n.gpu_ms, frame_ms);
    ImGui::TableNextColumn();
    ImGui::Text("%7.3f", n.cpu_ms);
    ImGui::TableNextColumn();
    if (n.calls != 1.0f) ImGui::TextDisabled("%.0f", n.calls);
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
        ImGui::TableSetupColumn("GPU ms", ImGuiTableColumnFlags_WidthFixed, ImGui::CalcTextSize("0000.000000").x);
        ImGui::TableSetupColumn("CPU ms", ImGuiTableColumnFlags_WidthFixed, ImGui::CalcTextSize("0000.0000").x);
        ImGui::TableSetupColumn("Calls",  ImGuiTableColumnFlags_WidthFixed, ImGui::CalcTextSize("00000").x);
        ImGui::TableHeadersRow();
        draw_row(0, g.nodes[0].gpu_ms);
        ImGui::EndTable();
    }
}

void draw_summary() {
    if (!g.node_count || !g.nodes[0].shown) {
        ImGui::TextDisabled("Measuring...");
        return;
    }
    const Node& root = g.nodes[0];
    ImGui::Text("Per frame: GPU %.2f ms, CPU %.2f ms", root.gpu_ms, root.cpu_ms);
    for (int c = root.first_child; c >= 0; c = g.nodes[c].next_sibling) {
        const Node& n = g.nodes[c];
        if (n.shown) ImGui::Text("  %-22s GPU %7.3f ms   CPU %7.3f ms", n.name, n.gpu_ms, n.cpu_ms);
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

// The section and those shown within it, in the order of the table: all of them, or with 'subtree' those named so,
// what is within them and the frame
static void write_csv_node(md_writer_i* w, int i, const char* subtree, const char* first, bool inside) {
    const Node& n = g.nodes[i];
    inside = inside || !subtree || strcmp(n.name, subtree) == 0;
    if (inside || i == 0) {
        // Its path: the names from the root, written before the rows within it are (so one buffer does for all of them)
        static char path[MAX_DEPTH * NAME_LEN];
        int chain[MAX_DEPTH + 1];
        int count = 0;
        for (int p = i; p >= 0 && count < (int)ARRAY_SIZE(chain); p = g.nodes[p].parent) chain[count++] = p;
        size_t len = 0;
        path[0] = '\0';
        for (int k = count - 1; k >= 0 && len + 1 < sizeof(path); --k) {
            const int m = snprintf(path + len, sizeof(path) - len, k == count - 1 ? "%s" : "/%s", g.nodes[chain[k]].name);
            if (m < 0) break;
            len = MIN(len + (size_t)m, sizeof(path) - 1);
        }
        if (first) {
            write_csv_field(w, first);
            md_write_char(w, ',');
        }
        write_csv_field(w, path);
        md_writef(w, ",%.4f,%.4f,%.4f,%.2f\n", n.gpu_ms, n.gpu_self_ms, n.cpu_ms, n.calls);
    }
    for (int c = n.first_child; c >= 0; c = g.nodes[c].next_sibling) {
        if (g.nodes[c].shown) write_csv_node(w, c, subtree, first, inside);
    }
}

void write_csv_header(md_writer_i* w, const char* info, const char* first_column) {
    char when[64] = "";
    const time_t t = time(nullptr);
    if (const struct tm* tm = localtime(&t)) strftime(when, sizeof(when), "%Y-%m-%d %H:%M:%S", tm);
    const char* vendor   = (const char*)glGetString(GL_VENDOR);
    const char* renderer = (const char*)glGetString(GL_RENDERER);
    const char* version  = (const char*)glGetString(GL_VERSION);

    md_writef(w, "# VIAMD profiler, %s\n", when);
#ifdef NDEBUG
    md_write_cstr(w, "# build: NDEBUG\n");
#else
    md_write_cstr(w, "# build: debug (asserts on)\n");
#endif
    md_writef(w, "# GL: %s | %s | %s\n", vendor ? vendor : "?", renderer ? renderer : "?", version ? version : "?");
    md_writef(w, "# per frame, averaged over %d frames (the GPU's over %d of them)\n", g.published_cpu_frames, g.published_gpu_frames);
    for (const char* line = info; line && *line;) {
        const char* end = strchr(line, '\n');
        const size_t n = end ? (size_t)(end - line) : strlen(line);
        md_write_cstr(w, "# ");
        md_write(w, line, n);
        md_write_char(w, '\n');
        line = end ? end + 1 : line + n;
    }
    if (first_column) {
        write_csv_field(w, first_column);
        md_write_char(w, ',');
    }
    md_write_cstr(w, "section,gpu_ms,gpu_self_ms,cpu_ms,calls\n");
}

void write_csv_rows(md_writer_i* w, const char* subtree, const char* first) {
    if (g.node_count && g.nodes[0].shown) write_csv_node(w, 0, subtree, first, false);
}

void write_csv(md_writer_i* w, const char* info) {
    write_csv_header(w, info);
    write_csv_rows(w);
}

int publish_count() { return g.publish_count; }

}  // namespace profiler
