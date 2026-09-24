#include "script_reference.h"

#include <imgui.h>
#include <imgui_md.h>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>

// The document is baked into the executable at configure time (see CMakeLists.txt, create_resources).
#include <script_reference.inl>

////////////////////////////////////////////////////////////////////////////////
// Outline parsing
////////////////////////////////////////////////////////////////////////////////

namespace {

std::string to_lower(const std::string& s) {
    std::string r = s;
    for (char& c : r) c = (char)tolower((unsigned char)c);
    return r;
}

// Turns the source of a heading into the text a reader sees (no backticks, emphasis markers or link targets).
std::string plain_heading(const std::string& s) {
    std::string out;
    for (size_t i = 0; i < s.size(); ++i) {
        const char c = s[i];
        if (c == '`' || c == '*') continue;
        if (c == '\\' && i + 1 < s.size()) { out += s[++i]; continue; }
        if (c == '[') continue;
        if (c == ']') {
            if (i + 1 < s.size() && s[i + 1] == '(') {
                const size_t close = s.find(')', i);
                if (close != std::string::npos) i = close;
            }
            continue;
        }
        out += c;
    }
    return out;
}

// GitHub style anchor: lower case, spaces become '-', punctuation is dropped.
std::string make_slug(const std::string& plain) {
    std::string slug;
    for (unsigned char c : plain) {
        if (c == ' ') slug += '-';
        else if (c == '-' || c == '_') slug += (char)c;
        else if (c >= 0x80) slug += (char)c;
        else if (isalnum(c)) slug += (char)tolower(c);
    }
    return slug;
}

std::string trim(const std::string& s) {
    size_t b = 0, e = s.size();
    while (b < e && isspace((unsigned char)s[b])) ++b;
    while (e > b && isspace((unsigned char)s[e - 1])) --e;
    return s.substr(b, e - b);
}

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    size_t b = 0;
    while (b <= s.size()) {
        size_t e = s.find(sep, b);
        if (e == std::string::npos) e = s.size();
        std::string t = trim(s.substr(b, e - b));
        if (!t.empty()) out.push_back(t);
        b = e + 1;
    }
    return out;
}

// Value of "key=value" in the text of a <!-- proc ... --> comment.
std::string meta_value(const std::string& line, const char* key) {
    const std::string k = std::string(key) + "=";
    size_t p = 0;
    while ((p = line.find(k, p)) != std::string::npos) {
        if (p == 0 || isspace((unsigned char)line[p - 1])) {
            const size_t b = p + k.size();
            size_t e = b;
            while (e < line.size() && !isspace((unsigned char)line[e])) ++e;
            return line.substr(b, e - b);
        }
        p += k.size();
    }
    return "";
}

// A fence opens with 3+ backticks or tildes and closes with at least as many of the same character.
struct FenceState {
    char ch = 0;
    int  len = 0;
    bool update(const std::string& line) {
        size_t i = 0;
        while (i < line.size() && i < 4 && line[i] == ' ') ++i;
        if (i >= 4 || i >= line.size()) return ch != 0;
        const char c = line[i];
        if (c != '`' && c != '~') return ch != 0;
        int n = 0;
        while (i + n < line.size() && line[i + n] == c) ++n;
        if (n < 3) return ch != 0;
        if (ch == 0) {
            ch = c; len = n;
        } else if (c == ch && n >= len && trim(line.substr(i + n)).empty()) {
            ch = 0; len = 0;
        }
        return true; // this line belongs to the fence (opening or closing)
    }
    bool inside() const { return ch != 0; }
};

} // namespace

ScriptReference::ScriptReference() {}

void ScriptReference::set_document(const char* markdown, size_t size) {
    m_doc.assign(markdown, size);
    m_pages.clear();
    m_sections.clear();
    m_anchors.clear();
    m_anchors_lower.clear();
    m_page = 0;
    m_current = 0;
    m_history.clear();
    m_history_pos = -1;

    Page first;
    first.title = "Overview";
    m_pages.push_back(first);

    std::unordered_map<std::string, int> slug_count;
    std::vector<size_t> section_begin;
    FenceState fence;
    bool seen_heading = false;
    int cur = -1; // section that metadata comments attach to

    size_t pos = 0;
    while (pos < m_doc.size()) {
        size_t eol = m_doc.find('\n', pos);
        if (eol == std::string::npos) eol = m_doc.size();
        const std::string line = m_doc.substr(pos, eol - pos);
        const size_t line_begin = pos;
        pos = eol + 1;

        const bool was_inside = fence.inside();
        fence.update(line);
        if (was_inside || fence.inside()) continue;

        // ATX heading: up to 3 spaces, 1-6 '#', then a space or the end of the line
        size_t i = 0;
        while (i < line.size() && i < 4 && line[i] == ' ') ++i;
        int level = 0;
        while (i + level < line.size() && line[i + level] == '#') ++level;
        const bool is_heading = i < 4 && level >= 1 && level <= 6 && (i + level == line.size() || line[i + level] == ' ');

        if (is_heading) {
            std::string raw = trim(line.substr(i + level));
            while (!raw.empty() && raw.back() == '#') raw.pop_back();
            raw = trim(raw);

            Section s;
            s.title = plain_heading(raw);
            s.slug = make_slug(s.title);
            const int n = slug_count[s.slug]++;
            if (n > 0) s.slug += "-" + std::to_string(n);
            s.level = level;

            if (level <= 2 && !(level == 1 && !seen_heading)) {
                m_pages.back().end = line_begin;
                Page p;
                p.title = s.title;
                p.begin = line_begin;
                m_pages.push_back(p);
            } else if (level == 1 && !seen_heading) {
                m_pages.back().title = s.title;
            }
            seen_heading = true;

            s.page = (int)m_pages.size() - 1;
            s.heading_index = (int)m_pages.back().sections.size();
            m_pages.back().sections.push_back((int)m_sections.size());
            cur = (int)m_sections.size();
            m_sections.push_back(s);
            section_begin.push_back(line_begin);
            continue;
        }

        if (cur < 0) continue;
        Section& s = m_sections[cur];

        // <!-- proc name=x aliases=a,b category=c -->
        if (line.compare(0, 9, "<!-- proc") == 0) {
            const std::string name = meta_value(line, "name");
            if (!name.empty()) s.proc = name;
            s.category = meta_value(line, "category");
            for (const std::string& a : split(meta_value(line, "aliases"), ',')) s.aliases.push_back(a);
        }
        // <a id="x"></a>
        if (line.compare(0, 5, "<a id") == 0 || line.compare(0, 7, "<a name") == 0) {
            size_t p = 0;
            while ((p = line.find("id=\"", p)) != std::string::npos) {
                p += 4;
                const size_t e = line.find('"', p);
                if (e == std::string::npos) break;
                const std::string id = line.substr(p, e - p);
                if (!id.empty() && std::find(s.aliases.begin(), s.aliases.end(), id) == s.aliases.end()) s.aliases.push_back(id);
                p = e;
            }
        }
    }
    m_pages.back().end = m_doc.size();

    // Section bodies for searching, and anchors
    for (size_t i = 0; i < m_sections.size(); ++i) {
        const size_t b = section_begin[i];
        const size_t e = (i + 1 < m_sections.size()) ? section_begin[i + 1] : m_doc.size();
        m_sections[i].body_lower = to_lower(m_doc.substr(b, e - b));
    }
    for (size_t i = 0; i < m_sections.size(); ++i) m_anchors.emplace(m_sections[i].slug, (int)i);
    for (size_t i = 0; i < m_sections.size(); ++i) {
        if (!m_sections[i].proc.empty()) m_anchors.emplace(m_sections[i].proc, (int)i);
        for (const std::string& a : m_sections[i].aliases) m_anchors.emplace(a, (int)i);
    }
    for (const auto& kv : m_anchors) m_anchors_lower.emplace(to_lower(kv.first), kv.second);

    m_history.push_back(0);
    m_history_pos = 0;
}

int ScriptReference::find_anchor(const std::string& name) const {
    auto it = m_anchors.find(name);
    if (it != m_anchors.end()) return it->second;
    it = m_anchors_lower.find(to_lower(name));
    if (it != m_anchors_lower.end()) return it->second;
    return -1;
}

////////////////////////////////////////////////////////////////////////////////
// Navigation
////////////////////////////////////////////////////////////////////////////////

void ScriptReference::navigate(int section, bool push_history) {
    if (section < 0 || section >= (int)m_sections.size()) return;
    const Section& s = m_sections[section];
    m_page = s.page;
    m_current = section;
    m_target_heading = s.level >= 3 ? s.heading_index : -1;
    m_scroll_top = (m_target_heading < 0);
    m_scroll_settle = 4;
    m_reveal_nav = true;
    if (push_history && !(m_history_pos >= 0 && m_history[m_history_pos] == section)) {
        m_history.resize((size_t)m_history_pos + 1);
        m_history.push_back(section);
        m_history_pos = (int)m_history.size() - 1;
    }
}

void ScriptReference::history_step(int dir) {
    const int p = m_history_pos + dir;
    if (p < 0 || p >= (int)m_history.size()) return;
    m_history_pos = p;
    navigate(m_history[p], false);
}

bool ScriptReference::show(const char* topic) {
    if (!topic || !*topic) return false;
    const int idx = find_anchor(topic);
    if (idx < 0) return false;
    m_query[0] = 0;
    navigate(idx);
    return true;
}

void ScriptReference::search(const char* query) {
    snprintf(m_query, sizeof(m_query), "%s", query ? query : "");
    m_query_applied.clear();
    m_focus_search = false;
}

void ScriptReference::follow_link(const std::string& url) {
    if (!url.empty() && url[0] == '#') {
        const std::string name = url.substr(1);
        const int idx = find_anchor(name);
        if (idx >= 0) {
            m_query[0] = 0;
            navigate(idx);
        } else {
            search(name.c_str());
        }
    } else if (on_open_url) {
        on_open_url(url);
    } else {
        ImGui::SetClipboardText(url.c_str());
    }
}

void ScriptReference::update_search_results() {
    const std::string q = to_lower(trim(m_query));
    if (q == m_query_applied) return;
    m_query_applied = q;
    m_hits.clear();
    if (q.empty()) return;

    for (size_t i = 0; i < m_sections.size(); ++i) {
        const Section& s = m_sections[i];
        int score = 0;
        const std::string title = to_lower(s.title);
        const std::string proc = to_lower(s.proc);
        auto consider = [&](const std::string& name, int exact, int prefix, int contains) {
            if (name.empty()) return;
            if (name == q) score = std::max(score, exact);
            else if (name.compare(0, q.size(), q) == 0) score = std::max(score, prefix);
            else if (name.find(q) != std::string::npos) score = std::max(score, contains);
        };
        consider(proc, 1000, 800, 500);
        for (const std::string& a : s.aliases) consider(to_lower(a), 950, 750, 450);
        consider(title, 900, 700, 600);
        if (score == 0) {
            size_t count = 0, p = 0;
            while (count < 20 && (p = s.body_lower.find(q, p)) != std::string::npos) { ++count; p += q.size(); }
            if (count) score = 100 + (int)count;
        }
        if (score > 0) m_hits.push_back({(int)i, score});
    }
    std::stable_sort(m_hits.begin(), m_hits.end(), [](const Hit& a, const Hit& b) { return a.score > b.score; });
    if (m_hits.size() > 100) m_hits.resize(100);
}

////////////////////////////////////////////////////////////////////////////////
// Rendering
////////////////////////////////////////////////////////////////////////////////

namespace {

struct Palette {
    ImU32 comment, number, function, keyword, string, param, type, text, punct;
};

Palette make_palette() {
    const ImVec4 bg = ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
    const bool dark = (bg.x + bg.y + bg.z) / 3.0f < 0.5f;
    auto c = [&](float r, float g, float b) {
        const float k = dark ? 1.0f : 0.55f;
        return ImGui::GetColorU32(ImVec4(r * k, g * k, b * k, 1.0f));
    };
    Palette p;
    p.comment  = c(0.47f, 0.62f, 0.47f);
    p.number   = c(0.74f, 0.63f, 0.96f);
    p.function = c(0.55f, 0.80f, 1.00f);
    p.keyword  = c(0.96f, 0.56f, 0.66f);
    p.string   = c(0.92f, 0.77f, 0.46f);
    p.param    = c(0.72f, 0.86f, 0.72f);
    p.type     = c(0.46f, 0.86f, 0.76f);
    p.text     = ImGui::GetColorU32(ImGuiCol_Text);
    ImVec4 punct = ImGui::GetStyle().Colors[ImGuiCol_Text];
    punct.w *= 0.75f;
    p.punct    = ImGui::GetColorU32(punct);
    return p;
}

bool is_ident_start(char c) { return isalpha((unsigned char)c) || c == '_'; }
bool is_ident(char c) { return isalnum((unsigned char)c) || c == '_'; }

bool is_keyword(const std::string& w) {
    static const char* kw[] = {"in", "out", "and", "or", "not", "xor", "true", "false"};
    for (const char* k : kw) if (w == k) return true;
    return false;
}

bool is_type_name(const std::string& w) {
    static const char* types[] = {"int", "float", "bool", "string", "irange", "frange", "bitfield", "position",
                                  "distribution", "volume", "vec2", "vec3", "vec4", "mat3", "mat4", "any"};
    for (const char* t : types) if (w == t) return true;
    return false;
}

// Draws one line of code with simple colouring and returns its width.
// 'script' selects mdscript rules ('#' comments, named arguments); otherwise the line is a signature.
float draw_code_line(ImDrawList* dl, ImVec2 pos, const char* line, size_t n, bool script, const Palette& pal) {
    float x = pos.x;
    auto emit = [&](size_t b, size_t e, ImU32 col) {
        if (e <= b) return;
        dl->AddText(ImVec2(x, pos.y), col, line + b, line + e);
        x += ImGui::CalcTextSize(line + b, line + e).x;
    };

    int depth = 0;
    bool after_arrow = false;
    size_t i = 0;
    while (i < n) {
        const char c = line[i];
        if (script && c == '#') {
            emit(i, n, pal.comment);
            break;
        }
        if (c == '"' || c == '\'') {
            size_t j = i + 1;
            while (j < n && line[j] != c) ++j;
            if (j < n) ++j;
            emit(i, j, pal.string);
            i = j;
        } else if (isdigit((unsigned char)c) || (c == '.' && i + 1 < n && isdigit((unsigned char)line[i + 1]))) {
            size_t j = i + 1;
            while (j < n && (isalnum((unsigned char)line[j]) || line[j] == '.')) ++j;
            emit(i, j, pal.number);
            i = j;
        } else if (is_ident_start(c)) {
            size_t j = i + 1;
            while (j < n && is_ident(line[j])) ++j;
            const std::string w(line + i, j - i);
            size_t k = j;
            while (k < n && line[k] == ' ') ++k;
            ImU32 col = pal.text;
            if (script) {
                if (is_keyword(w)) col = pal.keyword;
                else if (k < n && line[k] == '(') col = pal.function;
                else if (depth > 0 && k < n && line[k] == '=' && !(k + 1 < n && line[k + 1] == '=')) col = pal.param;
            } else {
                if (after_arrow || is_type_name(w)) col = pal.type;
                else if (i == 0 && k < n && line[k] == '(') col = pal.function;
                else if (k < n && line[k] == ':') col = pal.param;
            }
            emit(i, j, col);
            i = j;
        } else {
            if (c == '(' || c == '[') ++depth;
            if (c == ')' || c == ']') depth = std::max(0, depth - 1);
            size_t j = i + 1;
            if (!script && c == '-' && j < n && line[j] == '>') { after_arrow = true; ++j; }
            emit(i, j, script ? pal.text : (after_arrow ? pal.punct : pal.punct));
            i = j;
        }
    }
    return x - pos.x;
}

} // namespace

struct ScriptReferenceRenderer : imgui_md {
    ScriptReference& ref;
    explicit ScriptReferenceRenderer(ScriptReference& r) : ref(r) {}

    // Blocks of code that come from the reference get Copy / Insert buttons and colouring.
    void render_code_block(const std::string& lang, const std::string& code) override {
        const ImGuiStyle& style = ImGui::GetStyle();
        const float fs = ImGui::GetFontSize();
        const float lh = ImGui::GetTextLineHeight();
        const float pad = fs * 0.5f;
        const float header_h = fs + style.FramePadding.y * 2.0f + 4.0f;
        const bool script = (lang == "mdscript");

        // Split into lines and measure
        std::vector<std::pair<size_t, size_t>> lines;
        float max_w = 0.0f;
        for (size_t b = 0; b <= code.size();) {
            size_t e = code.find('\n', b);
            if (e == std::string::npos) e = code.size();
            lines.emplace_back(b, e);
            max_w = std::max(max_w, ImGui::CalcTextSize(code.c_str() + b, code.c_str() + e).x);
            b = e + 1;
        }

        const float avail = ImGui::GetContentRegionAvail().x;
        const bool need_hscroll = max_w + pad * 2.0f > avail;
        const float body_h = (float)lines.size() * lh + pad + (need_hscroll ? style.ScrollbarSize : 0.0f);
        const float total_h = header_h + body_h;

        const ImVec2 p = ImGui::GetCursorScreenPos();
        if (!ImGui::IsRectVisible(ImVec2(avail, total_h))) {
            ImGui::Dummy(ImVec2(avail, total_h));
            return;
        }

        ImDrawList* dl = ImGui::GetWindowDrawList();
        ImVec4 bg = style.Colors[ImGuiCol_FrameBg];
        ImVec4 head = style.Colors[ImGuiCol_TableHeaderBg];
        dl->AddRectFilled(p, ImVec2(p.x + avail, p.y + total_h), ImGui::GetColorU32(bg), 4.0f);
        dl->AddRectFilled(p, ImVec2(p.x + avail, p.y + header_h), ImGui::GetColorU32(head), 4.0f, ImDrawFlags_RoundCornersTop);

        // Header: language on the left, buttons on the right
        ImGui::Dummy(ImVec2(avail, header_h));
        const char* label = lang.empty() ? "code" : lang.c_str();
        ImVec4 dim = style.Colors[ImGuiCol_TextDisabled];
        dl->AddText(ImVec2(p.x + pad, p.y + (header_h - fs) * 0.5f), ImGui::GetColorU32(dim), label);

        const int id = m_code_index;
        const bool feedback = (ref.m_feedback_id == id) && ImGui::GetTime() < ref.m_feedback_until;
        const bool can_insert = script && (bool)ref.on_insert_code;
        const char* copy_label = (feedback && !strcmp(ref.m_feedback_text, "Copied")) ? "Copied" : "Copy";
        const char* ins_label = (feedback && !strcmp(ref.m_feedback_text, "Inserted")) ? "Inserted" : "Insert";
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(style.FramePadding.x, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
        const float bh = fs + 2.0f;
        float bw = ImGui::CalcTextSize(copy_label).x + style.FramePadding.x * 2.0f;
        if (can_insert) bw += ImGui::CalcTextSize(ins_label).x + style.FramePadding.x * 2.0f + style.ItemSpacing.x;
        ImGui::SetCursorScreenPos(ImVec2(p.x + avail - bw - pad * 0.5f, p.y + (header_h - bh) * 0.5f));
        if (ImGui::Button(copy_label)) {
            ImGui::SetClipboardText(code.c_str());
            ref.m_feedback_id = id;
            ref.m_feedback_until = ImGui::GetTime() + 1.2;
            ref.m_feedback_text = "Copied";
        }
        if (can_insert) {
            ImGui::SameLine();
            if (ImGui::Button(ins_label)) {
                ref.on_insert_code(code);
                ref.m_feedback_id = id;
                ref.m_feedback_until = ImGui::GetTime() + 1.2;
                ref.m_feedback_text = "Inserted";
            }
        }
        ImGui::PopStyleColor();
        ImGui::PopStyleVar();

        // Body: scrolls horizontally when a line is too long
        ImGui::SetCursorScreenPos(ImVec2(p.x, p.y + header_h));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(pad, pad * 0.5f));
        ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0, 0, 0, 0));
        if (ImGui::BeginChild("##code", ImVec2(avail, body_h), ImGuiChildFlags_AlwaysUseWindowPadding, ImGuiWindowFlags_HorizontalScrollbar)) {
            static thread_local Palette pal;
            pal = make_palette();
            ImDrawList* cdl = ImGui::GetWindowDrawList();
            const ImVec2 origin = ImGui::GetCursorScreenPos();
            for (size_t l = 0; l < lines.size(); ++l) {
                const ImVec2 lp(origin.x, origin.y + (float)l * lh);
                draw_code_line(cdl, lp, code.c_str() + lines[l].first, lines[l].second - lines[l].first, script, pal);
            }
            ImGui::Dummy(ImVec2(max_w, (float)lines.size() * lh));
        }
        ImGui::EndChild();
        ImGui::PopStyleColor();
        ImGui::PopStyleVar();
    }

    void on_heading(int, const std::string&, float screen_y) override {
        if (ref.m_target_heading == ref.m_headings_rendered) {
            const ImGuiStyle& style = ImGui::GetStyle();
            // The target stays armed for a few frames: text wraps differently once the scrollbar has appeared, so the
            // position of the heading can still move after the first frame.
            ImGui::SetScrollY(screen_y - ImGui::GetWindowPos().y + ImGui::GetScrollY() - style.WindowPadding.y);
        }
        ++ref.m_headings_rendered;
    }

    void open_url(const std::string& url) override { ref.follow_link(url); }
};

void ScriptReference::draw_toolbar() {
    const float fs = ImGui::GetFontSize();

    ImGui::BeginDisabled(m_history_pos <= 0);
    if (ImGui::ArrowButton("##back", ImGuiDir_Left)) history_step(-1);
    ImGui::EndDisabled();
    ImGui::SameLine();
    ImGui::BeginDisabled(m_history_pos + 1 >= (int)m_history.size());
    if (ImGui::ArrowButton("##forward", ImGuiDir_Right)) history_step(+1);
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (ImGui::Button("Contents")) {
        m_query[0] = 0;
        navigate(0);
    }
    ImGui::SameLine();

    if (m_focus_search) {
        ImGui::SetKeyboardFocusHere();
        m_focus_search = false;
    }
    ImGui::SetNextItemWidth(std::max(fs * 8.0f, ImGui::GetContentRegionAvail().x));
    const bool enter = ImGui::InputTextWithHint("##search", "Search the reference (procedure, keyword...)", m_query, sizeof(m_query),
                                                ImGuiInputTextFlags_EscapeClearsAll | ImGuiInputTextFlags_EnterReturnsTrue);
    if (enter) {
        update_search_results();
        if (!m_hits.empty()) navigate(m_hits.front().section);
    }

    if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows) && ImGui::IsKeyPressed(ImGuiKey_F) && ImGui::GetIO().KeyCtrl) {
        m_focus_search = true;
    }
}

void ScriptReference::draw_sidebar() {
    const ImGuiStyle& style = ImGui::GetStyle();
    const float fs = ImGui::GetFontSize();
    const ImVec4 dim = style.Colors[ImGuiCol_TextDisabled];

    if (m_query[0] != 0) {
        update_search_results();
        if (m_hits.empty()) {
            ImGui::TextDisabled("No matches.");
            return;
        }
        ImGui::TextDisabled("%d match%s", (int)m_hits.size(), m_hits.size() == 1 ? "" : "es");
        for (const Hit& h : m_hits) {
            const Section& s = m_sections[h.section];
            ImGui::PushID(h.section);
            const std::string label = s.proc.empty() ? s.title : s.proc;
            if (ImGui::Selectable(label.c_str(), h.section == m_current)) navigate(h.section);
            if (m_pages[s.page].title != s.title) {
                ImGui::SameLine();
                ImGui::TextColored(dim, "%s", m_pages[s.page].title.c_str());
            }
            ImGui::PopID();
        }
        return;
    }

    for (size_t pi = 0; pi < m_pages.size(); ++pi) {
        const Page& page = m_pages[pi];
        const bool open = ((int)pi == m_page);
        const int head = page.sections.empty() ? 0 : page.sections.front();
        ImGui::PushID((int)pi);
        ImVec4 page_col = ImGui::GetStyleColorVec4(open ? ImGuiCol_Text : ImGuiCol_TextDisabled);
        if (!open) page_col.w = std::min(1.0f, page_col.w * 1.6f);
        ImGui::PushStyleColor(ImGuiCol_Text, page_col);
        const bool clicked = ImGui::Selectable(pi == 0 ? "Overview" : page.title.c_str(), open && m_current == head);
        ImGui::PopStyleColor();
        if (open && m_reveal_nav && m_current == head) ImGui::SetScrollHereY(0.3f);
        if (clicked) navigate(head);
        if (ImGui::IsItemHovered() && ImGui::CalcTextSize(page.title.c_str()).x > ImGui::GetContentRegionAvail().x) ImGui::SetTooltip("%s", page.title.c_str());

        if (open) {
            ImGui::Indent(fs * 0.8f);
            for (int si : page.sections) {
                const Section& s = m_sections[si];
                if (s.level < 3) continue;
                ImGui::PushID(si);
                const bool is_proc = !s.proc.empty();
                if (!is_proc) ImGui::PushStyleColor(ImGuiCol_Text, dim);
                if (ImGui::Selectable(s.title.c_str(), si == m_current)) navigate(si);
                if (!is_proc) ImGui::PopStyleColor();
                if (si == m_current && m_reveal_nav) ImGui::SetScrollHereY(0.5f);
                ImGui::PopID();
            }
            ImGui::Unindent(fs * 0.8f);
        }
        ImGui::PopID();
    }
    m_reveal_nav = false;
}

void ScriptReference::draw_content() {
    if (m_page < 0 || m_page >= (int)m_pages.size()) return;
    const Page& page = m_pages[m_page];

    if (m_scroll_top) {
        ImGui::SetScrollY(0.0f);
        m_scroll_top = false;
    }

    m_headings_rendered = 0;
    ScriptReferenceRenderer md(*this);
    md.print(m_doc.data() + page.begin, m_doc.data() + page.end);

    if (m_target_heading >= 0 && --m_scroll_settle <= 0) m_target_heading = -1;

    // Previous / next page
    ImGui::Dummy(ImVec2(0.0f, ImGui::GetFontSize() * 0.5f));
    ImGui::Separator();
    if (m_page > 0) {
        if (ImGui::Button(("<  " + (m_page == 1 ? std::string("Overview") : m_pages[m_page - 1].title)).c_str())) {
            navigate(m_pages[m_page - 1].sections.front());
        }
        ImGui::SameLine();
    }
    if (m_page + 1 < (int)m_pages.size()) {
        if (ImGui::Button((m_pages[m_page + 1].title + "  >").c_str())) {
            navigate(m_pages[m_page + 1].sections.front());
        }
    }
    ImGui::Dummy(ImVec2(0.0f, ImGui::GetFontSize()));
}

void ScriptReference::draw() {
    if (!has_document()) {
        ImGui::TextDisabled("The script reference is not available in this build.");
        return;
    }
    const float fs = ImGui::GetFontSize();
    if (m_nav_width <= 0.0f) m_nav_width = fs * 17.0f;

    draw_toolbar();
    ImGui::Spacing();

    const ImGuiStyle& style = ImGui::GetStyle();
    const float total_w = ImGui::GetContentRegionAvail().x;
    const float splitter_w = 6.0f;
    m_nav_width = std::clamp(m_nav_width, fs * 8.0f, std::max(fs * 8.0f, total_w * 0.6f));

    if (ImGui::BeginChild("##ref_nav", ImVec2(m_nav_width, 0.0f), ImGuiChildFlags_Borders)) {
        draw_sidebar();
    }
    ImGui::EndChild();

    ImGui::SameLine(0.0f, 0.0f);
    ImGui::InvisibleButton("##ref_split", ImVec2(splitter_w, ImGui::GetContentRegionAvail().y));
    if (ImGui::IsItemActive()) m_nav_width += ImGui::GetIO().MouseDelta.x;
    if (ImGui::IsItemHovered() || ImGui::IsItemActive()) ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeEW);
    ImGui::SameLine(0.0f, 0.0f);

    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(fs * 0.9f, fs * 0.6f));
    if (ImGui::BeginChild("##ref_content", ImVec2(0.0f, 0.0f), ImGuiChildFlags_Borders | ImGuiChildFlags_AlwaysUseWindowPadding)) {
        draw_content();
    }
    ImGui::EndChild();
    ImGui::PopStyleVar();
    (void)style;
}

void ScriptReference::draw_window(bool* p_open, const char* title) {
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 60.0f, ImGui::GetFontSize() * 40.0f), ImGuiCond_FirstUseEver);
    if (ImGui::Begin(title, p_open)) {
        draw();
    }
    ImGui::End();
}

////////////////////////////////////////////////////////////////////////////////

ScriptReference& script_reference() {
    static ScriptReference ref;
    static bool init = false;
    if (!init) {
        init = true;
        ref.set_document((const char*)script_reference_md, script_reference_md_size - 1); // drop the terminating zero
    }
    return ref;
}
