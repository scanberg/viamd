#include "script_reference.h"

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_str.h>

#include <imgui.h>
#include <imgui_md.h>

#include <stdio.h>
#include <string.h>
#include <algorithm>

// The document is baked into the executable at configure time (see CMakeLists.txt, create_resources).
#include <script_reference.inl>

// Everything the viewer knows about the document is a view into the embedded text, except for the few strings that
// have to be generated (heading titles without markup, anchors). Those are made once, when the document is loaded,
// in an arena that lives as long as the application. Drawing does not allocate.

namespace script_reference {
namespace {

struct Section {
    str_t title;            // plain text of the heading (zero terminated)
    str_t slug;             // GitHub style anchor of the heading (zero terminated)
    str_t proc;             // procedure name from <!-- proc name= --> (zero terminated), empty if not a procedure
    uint32_t beg, end;      // byte range in the document, from the heading to the next one
    int level;              // heading level, 1..6
    int page;
    int heading_index;      // index among the headings of its page, in document order
};

struct Page {
    str_t title;            // zero terminated
    uint32_t beg, end;      // byte range in the document
    int first_section;      // the sections of a page are consecutive
    int num_sections;
};

enum AnchorKind { AnchorKind_Slug, AnchorKind_Proc, AnchorKind_Alias };

// Every name that leads to a section. There are a few hundred, and they are only looked up when the user navigates,
// so a linear scan in priority order (slugs, then procedure names, then aliases) is all that is needed.
struct Anchor {
    str_t name;
    int section;
    AnchorKind kind;
};

struct Hit {
    int section;
    int score;
};

enum Feedback { Feedback_None, Feedback_Copied, Feedback_Inserted };

constexpr int MAX_HISTORY = 128;
constexpr int MAX_HITS    = 100;
constexpr int MAX_QUERY   = 128;

struct Palette {
    ImU32 comment, number, function, keyword, string, param, type, text, punct;
};

struct Renderer : imgui_md {
    Palette pal = {};
    void render_code_block(const char* lang, const char* lang_end, const char* code, const char* code_end) override;
    void on_heading(int level, float screen_y) override;
    void open_url(const char* url, const char* url_end) override;
};

struct State {
    bool loaded = false;
    md_allocator_i* arena = nullptr;

    str_t    doc = {};
    Page*    pages = nullptr;       // md_array
    Section* sections = nullptr;    // md_array
    Anchor*  anchors = nullptr;     // md_array

    // View state
    int  page = 0;
    int  current = 0;               // section that was navigated to last
    int  target_heading = -1;       // heading (on the current page) to scroll to, -1 for none
    int  scroll_settle = 0;         // frames the scroll target stays armed
    int  headings_rendered = 0;     // headings seen while drawing the current page (for scroll targets)
    bool scroll_top = false;
    bool reveal_nav = false;        // scroll the sidebar to the current entry
    bool focus_search = false;
    float nav_width = 0.0f;         // sidebar width in pixels (0 -> derive from font size)

    int history[MAX_HISTORY] = {};
    int history_len = 0;
    int history_pos = -1;

    char query[MAX_QUERY] = "";
    char query_applied[MAX_QUERY] = "";
    bool query_dirty = false;
    Hit* hits = nullptr;            // md_array, capacity reserved on load
    int* scores = nullptr;          // md_array, one per section

    // Feedback for the copy / insert buttons of code blocks
    int      feedback_id = -1;
    double   feedback_until = 0.0;
    Feedback feedback = Feedback_None;

    // Example code the user asked to insert (md_array on the heap, reused)
    char* insert_buf = nullptr;
    bool  insert_pending = false;

    Renderer renderer;
};

State g;

////////////////////////////////////////////////////////////////////////////////
// String helpers
////////////////////////////////////////////////////////////////////////////////

bool is_space(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; }
bool is_alnum(char c) { return is_alpha(c) || is_digit(c); }

// Case-insensitive str_find_str
bool find_ignore_case(size_t* loc, str_t hay, str_t needle) {
    if (needle.len == 0 || needle.len > hay.len) return false;
    for (size_t i = 0; i + needle.len <= hay.len; ++i) {
        if (str_eq_ignore_case({hay.ptr + i, needle.len}, needle)) {
            if (loc) *loc = i;
            return true;
        }
    }
    return false;
}

bool begins_with_ignore_case(str_t str, str_t prefix) {
    return str.len >= prefix.len && str_eq_ignore_case(str_substr(str, 0, prefix.len), prefix);
}

// Turns the source of a heading into the text a reader sees (no backticks, emphasis markers or link targets).
str_t plain_heading(str_t s, md_allocator_i* alloc) {
    char* out = (char*)md_alloc(alloc, s.len + 1);
    size_t n = 0;
    for (size_t i = 0; i < s.len; ++i) {
        const char c = s.ptr[i];
        if (c == '`' || c == '*' || c == '[') continue;
        if (c == '\\' && i + 1 < s.len) { out[n++] = s.ptr[++i]; continue; }
        if (c == ']') {
            // Skip the target of a link: "[text](target)"
            size_t close;
            if (i + 1 < s.len && s.ptr[i + 1] == '(' && str_find_char(&close, str_substr(s, i), ')')) i += close;
            continue;
        }
        out[n++] = c;
    }
    out[n] = '\0';
    return {out, n};
}

// GitHub style anchor: lower case, spaces become '-', punctuation is dropped. 'suffix' > 0 appends "-suffix", which is
// how repeated headings are told apart.
str_t make_slug(str_t plain, int suffix, md_allocator_i* alloc) {
    char* out = (char*)md_alloc(alloc, plain.len + 16);
    size_t n = 0;
    for (size_t i = 0; i < plain.len; ++i) {
        const unsigned char c = (unsigned char)plain.ptr[i];
        if (c == ' ') out[n++] = '-';
        else if (c == '-' || c == '_' || c >= 0x80) out[n++] = (char)c;
        else if (is_alnum((char)c)) out[n++] = (char)to_lower(c);
    }
    if (suffix > 0) n += (size_t)snprintf(out + n, 16, "-%d", suffix);
    out[n] = '\0';
    return {out, n};
}

// Value of "key=value" in the text of a <!-- proc ... --> comment.
str_t meta_value(str_t line, str_t key) {
    for (size_t p = 0; p + key.len < line.len; ++p) {
        if ((p == 0 || is_space(line.ptr[p - 1])) && str_eq(str_substr(line, p, key.len), key) && line.ptr[p + key.len] == '=') {
            const size_t b = p + key.len + 1;
            size_t e = b;
            while (e < line.len && !is_space(line.ptr[e])) ++e;
            return str_substr(line, b, e - b);
        }
    }
    return {};
}

// A fence opens with 3+ backticks or tildes and closes with at least as many of the same character.
struct Fence {
    char ch = 0;
    size_t len = 0;

    void update(str_t line) {
        size_t i = 0;
        while (i < line.len && i < 4 && line.ptr[i] == ' ') ++i;
        if (i >= 4 || i >= line.len) return;
        const char c = line.ptr[i];
        if (c != '`' && c != '~') return;
        size_t n = 0;
        while (i + n < line.len && line.ptr[i + n] == c) ++n;
        if (n < 3) return;
        if (ch == 0) {
            ch = c;
            len = n;
        } else if (c == ch && n >= len && str_empty(str_trim(str_substr(line, i + n)))) {
            ch = 0;
            len = 0;
        }
    }
    bool inside() const { return ch != 0; }
};

////////////////////////////////////////////////////////////////////////////////
// Outline parsing
////////////////////////////////////////////////////////////////////////////////

void load(str_t doc) {
    State& s = g;
    ASSERT(!s.loaded);
    s.loaded = true;
    s.arena = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    md_allocator_i* arena = s.arena;
    s.doc = doc;

    // Aliases in document order. They are appended to the anchors after all slugs and procedure names, which is the
    // priority order for lookups.
    Anchor* aliases = nullptr;
    str_t*  base_slugs = nullptr;  // slug of every heading before a "-n" suffix, for numbering repeats

    Page overview = {};
    overview.title = STR_LIT("Overview");
    md_array_push(s.pages, overview, arena);

    Fence fence;
    bool seen_heading = false;

    str_t rest = doc;
    str_t line;
    while (true) {
        const uint32_t line_beg = (uint32_t)(rest.ptr - doc.ptr);
        if (!str_extract_line(&line, &rest)) break;

        const bool was_inside = fence.inside();
        fence.update(line);
        if (was_inside || fence.inside()) continue;

        // ATX heading: up to 3 spaces, 1-6 '#', then a space or the end of the line
        size_t i = 0;
        while (i < line.len && i < 4 && line.ptr[i] == ' ') ++i;
        int level = 0;
        while (i + level < line.len && line.ptr[i + level] == '#') ++level;
        const bool is_heading = i < 4 && 1 <= level && level <= 6 && (i + level == line.len || line.ptr[i + level] == ' ');

        if (is_heading) {
            str_t raw = str_trim(str_substr(line, i + level));
            while (raw.len > 0 && raw.ptr[raw.len - 1] == '#') --raw.len;
            raw = str_trim(raw);

            Section sec = {};
            sec.title = plain_heading(raw, arena);
            sec.level = level;

            str_t base = make_slug(sec.title, 0, arena);
            int repeats = 0;
            for (size_t j = 0; j < md_array_size(base_slugs); ++j) {
                if (str_eq(base_slugs[j], base)) ++repeats;
            }
            md_array_push(base_slugs, base, arena);
            sec.slug = repeats ? make_slug(sec.title, repeats, arena) : base;

            if (level == 1 && !seen_heading) {
                // The first top level heading names the document, it is the title of the overview page
                s.pages[0].title = sec.title;
            } else if (level <= 2) {
                md_array_back(s.pages).end = line_beg;
                Page p = {};
                p.title = sec.title;
                p.beg = line_beg;
                p.first_section = (int)md_array_size(s.sections);
                md_array_push(s.pages, p, arena);
            }
            seen_heading = true;

            Page& page = md_array_back(s.pages);
            sec.page = (int)md_array_size(s.pages) - 1;
            sec.heading_index = page.num_sections++;
            sec.beg = line_beg;
            if (md_array_size(s.sections) > 0) md_array_back(s.sections).end = line_beg;
            md_array_push(s.sections, sec, arena);
            continue;
        }

        if (md_array_size(s.sections) == 0) continue;
        const int cur = (int)md_array_size(s.sections) - 1;  // metadata attaches to the last heading
        Section& sec = s.sections[cur];

        // <!-- proc name=x aliases=a,b category=c -->
        if (str_begins_with(line, STR_LIT("<!-- proc"))) {
            const str_t name = meta_value(line, STR_LIT("name"));
            if (!str_empty(name)) sec.proc = str_copy(name, arena);
            str_t list = meta_value(line, STR_LIT("aliases"));
            while (!str_empty(list)) {
                size_t comma = list.len;
                str_find_char(&comma, list, ',');
                const str_t alias = str_trim(str_substr(list, 0, comma));
                if (!str_empty(alias)) md_array_push(aliases, (Anchor{alias, cur, AnchorKind_Alias}), arena);
                list = str_substr(list, MIN(comma + 1, list.len));
            }
        }
        // <a id="x"></a>
        if (str_begins_with(line, STR_LIT("<a id")) || str_begins_with(line, STR_LIT("<a name"))) {
            str_t rem = line;
            size_t p;
            while (str_find_str(&p, rem, STR_LIT("id=\""))) {
                rem = str_substr(rem, p + 4);
                size_t q;
                if (!str_find_char(&q, rem, '"')) break;
                const str_t id = str_substr(rem, 0, q);
                bool dup = str_empty(id);
                for (size_t j = 0; j < md_array_size(aliases) && !dup; ++j) {
                    dup = aliases[j].section == cur && str_eq(aliases[j].name, id);
                }
                if (!dup) md_array_push(aliases, (Anchor{id, cur, AnchorKind_Alias}), arena);
                rem = str_substr(rem, q);
            }
        }
    }
    md_array_back(s.pages).end = (uint32_t)doc.len;
    if (md_array_size(s.sections) > 0) md_array_back(s.sections).end = (uint32_t)doc.len;

    const int num_sections = (int)md_array_size(s.sections);
    for (int i = 0; i < num_sections; ++i) {
        md_array_push(s.anchors, (Anchor{s.sections[i].slug, i, AnchorKind_Slug}), arena);
    }
    size_t a = 0;
    for (int i = 0; i < num_sections; ++i) {
        if (!str_empty(s.sections[i].proc)) md_array_push(s.anchors, (Anchor{s.sections[i].proc, i, AnchorKind_Proc}), arena);
        for (; a < md_array_size(aliases) && aliases[a].section == i; ++a) md_array_push(s.anchors, aliases[a], arena);
    }

    md_array_resize(s.scores, (size_t)num_sections, arena);
    md_array_ensure(s.hits, (size_t)num_sections, arena);

    s.history[0] = 0;
    s.history_len = 1;
    s.history_pos = 0;
}

void ensure_loaded() {
    if (!g.loaded) {
        load({(const char*)script_reference_md, script_reference_md_size - 1});  // drop the terminating zero
    }
}

int find_anchor(str_t name) {
    for (size_t i = 0; i < md_array_size(g.anchors); ++i) {
        if (str_eq(g.anchors[i].name, name)) return g.anchors[i].section;
    }
    for (size_t i = 0; i < md_array_size(g.anchors); ++i) {
        if (str_eq_ignore_case(g.anchors[i].name, name)) return g.anchors[i].section;
    }
    return -1;
}

////////////////////////////////////////////////////////////////////////////////
// Navigation
////////////////////////////////////////////////////////////////////////////////

void navigate(int section, bool push_history = true) {
    State& s = g;
    if (section < 0 || section >= (int)md_array_size(s.sections)) return;
    const Section& sec = s.sections[section];
    s.page = sec.page;
    s.current = section;
    s.target_heading = sec.level >= 3 ? sec.heading_index : -1;
    s.scroll_top = (s.target_heading < 0);
    s.scroll_settle = 4;
    s.reveal_nav = true;
    if (push_history && s.history[s.history_pos] != section) {
        s.history_len = s.history_pos + 1;
        if (s.history_len == MAX_HISTORY) {
            memmove(s.history, s.history + 1, sizeof(int) * (MAX_HISTORY - 1));
            --s.history_len;
        }
        s.history[s.history_len++] = section;
        s.history_pos = s.history_len - 1;
    }
}

void history_step(int dir) {
    const int p = g.history_pos + dir;
    if (p < 0 || p >= g.history_len) return;
    g.history_pos = p;
    navigate(g.history[p], false);
}

void set_query(str_t q) {
    str_copy_to_char_buf(g.query, sizeof(g.query), q);
    g.query_dirty = true;
}

void follow_link(const char* url, const char* url_end) {
    const str_t u = {url, (size_t)(url_end - url)};
    if (u.len > 0 && u.ptr[0] == '#') {
        const str_t name = str_substr(u, 1);
        const int idx = find_anchor(name);
        if (idx >= 0) {
            g.query[0] = '\0';
            navigate(idx);
        } else {
            set_query(name);
        }
    } else {
        ImGui::SetClipboardText(url);  // zero terminated, see imgui_md
    }
}

void update_search_results() {
    State& s = g;
    const str_t q = str_trim(str_from_cstr(s.query));
    if (!s.query_dirty && str_eq_cstr(q, s.query_applied)) return;
    s.query_dirty = false;
    str_copy_to_char_buf(s.query_applied, sizeof(s.query_applied), q);
    md_array_shrink(s.hits, 0);
    if (str_empty(q)) return;

    const int num_sections = (int)md_array_size(s.sections);
    auto consider = [&q](int* score, str_t name, int exact, int prefix, int contains) {
        if (str_empty(name)) return;
        if (str_eq_ignore_case(name, q))            *score = MAX(*score, exact);
        else if (begins_with_ignore_case(name, q))  *score = MAX(*score, prefix);
        else if (find_ignore_case(nullptr, name, q)) *score = MAX(*score, contains);
    };

    for (int i = 0; i < num_sections; ++i) s.scores[i] = 0;
    for (size_t i = 0; i < md_array_size(s.anchors); ++i) {
        const Anchor& a = s.anchors[i];
        if (a.kind == AnchorKind_Alias) consider(&s.scores[a.section], a.name, 950, 750, 450);
    }
    for (int i = 0; i < num_sections; ++i) {
        const Section& sec = s.sections[i];
        int score = s.scores[i];
        consider(&score, sec.proc, 1000, 800, 500);
        consider(&score, sec.title, 900, 700, 600);
        if (score == 0) {
            // Mentioned in the text: the more often, the better (up to a point)
            str_t body = str_substr(s.doc, sec.beg, sec.end - sec.beg);
            int count = 0;
            size_t p;
            while (count < 20 && find_ignore_case(&p, body, q)) {
                ++count;
                body = str_substr(body, p + q.len);
            }
            if (count) score = 100 + count;
        }
        if (score > 0) md_array_push_no_grow(s.hits, (Hit{i, score}));
    }
    // Best first, document order among equals
    std::sort(s.hits, s.hits + md_array_size(s.hits), [](const Hit& a, const Hit& b) {
        return a.score != b.score ? a.score > b.score : a.section < b.section;
    });
    md_array_shrink(s.hits, MIN((int)md_array_size(s.hits), MAX_HITS));
}

////////////////////////////////////////////////////////////////////////////////
// Rendering
////////////////////////////////////////////////////////////////////////////////

Palette make_palette() {
    const ImVec4 bg = ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
    const bool dark = (bg.x + bg.y + bg.z) / 3.0f < 0.5f;
    auto c = [dark](float r, float g, float b) {
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

bool is_ident_start(char c) { return is_alpha(c) || c == '_'; }
bool is_ident(char c) { return is_alnum(c) || c == '_'; }

bool is_one_of(str_t word, const char* const* list, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        if (str_eq_cstr(word, list[i])) return true;
    }
    return false;
}

bool is_keyword(str_t w) {
    static const char* const kw[] = {"in", "out", "and", "or", "not", "xor", "true", "false"};
    return is_one_of(w, kw, ARRAY_SIZE(kw));
}

bool is_type_name(str_t w) {
    static const char* const types[] = {"int", "float", "bool", "string", "irange", "frange", "bitfield", "position",
                                        "distribution", "volume", "vec2", "vec3", "vec4", "mat3", "mat4", "any"};
    return is_one_of(w, types, ARRAY_SIZE(types));
}

// Draws one line of code with simple colouring.
// 'script' selects mdscript rules ('#' comments, named arguments); otherwise the line is a signature.
void draw_code_line(ImDrawList* dl, ImVec2 pos, str_t line, bool script, const Palette& pal) {
    const char* s = line.ptr;
    const size_t n = line.len;
    float x = pos.x;
    auto emit = [&](size_t b, size_t e, ImU32 col) {
        if (e <= b) return;
        dl->AddText(ImVec2(x, pos.y), col, s + b, s + e);
        x += ImGui::CalcTextSize(s + b, s + e).x;
    };

    int depth = 0;
    bool after_arrow = false;
    size_t i = 0;
    while (i < n) {
        const char c = s[i];
        if (script && c == '#') {
            emit(i, n, pal.comment);
            break;
        }
        if (c == '"' || c == '\'') {
            size_t j = i + 1;
            while (j < n && s[j] != c) ++j;
            if (j < n) ++j;
            emit(i, j, pal.string);
            i = j;
        } else if (is_digit(c) || (c == '.' && i + 1 < n && is_digit(s[i + 1]))) {
            size_t j = i + 1;
            while (j < n && (is_alnum(s[j]) || s[j] == '.')) ++j;
            emit(i, j, pal.number);
            i = j;
        } else if (is_ident_start(c)) {
            size_t j = i + 1;
            while (j < n && is_ident(s[j])) ++j;
            const str_t w = {s + i, j - i};
            size_t k = j;
            while (k < n && s[k] == ' ') ++k;
            ImU32 col = pal.text;
            if (script) {
                if (is_keyword(w)) col = pal.keyword;
                else if (k < n && s[k] == '(') col = pal.function;
                else if (depth > 0 && k < n && s[k] == '=' && !(k + 1 < n && s[k + 1] == '=')) col = pal.param;
            } else {
                if (after_arrow || is_type_name(w)) col = pal.type;
                else if (i == 0 && k < n && s[k] == '(') col = pal.function;
                else if (k < n && s[k] == ':') col = pal.param;
            }
            emit(i, j, col);
            i = j;
        } else {
            if (c == '(' || c == '[') ++depth;
            if (c == ')' || c == ']') depth = MAX(0, depth - 1);
            size_t j = i + 1;
            if (!script && c == '-' && j < n && s[j] == '>') { after_arrow = true; ++j; }
            emit(i, j, script ? pal.text : pal.punct);
            i = j;
        }
    }
}

void set_feedback(int id, Feedback f) {
    g.feedback_id = id;
    g.feedback = f;
    g.feedback_until = ImGui::GetTime() + 1.2;
}

// Blocks of code that come from the reference get Copy / Insert buttons and colouring.
void Renderer::render_code_block(const char* lang_beg, const char* lang_end, const char* code_beg, const char* code_end) {
    const ImGuiStyle& style = ImGui::GetStyle();
    const float fs = ImGui::GetFontSize();
    const float lh = ImGui::GetTextLineHeight();
    const float pad = fs * 0.5f;
    const float header_h = fs + style.FramePadding.y * 2.0f + 4.0f;
    const str_t lang = {lang_beg, (size_t)(lang_end - lang_beg)};
    const str_t code = {code_beg, (size_t)(code_end - code_beg)};
    const bool script = str_eq(lang, STR_LIT("mdscript"));

    // Measure
    int num_lines = 0;
    float max_w = 0.0f;
    {
        str_t rest = code, line;
        while (str_extract_line(&line, &rest)) {
            max_w = MAX(max_w, ImGui::CalcTextSize(line.ptr, line.ptr + line.len).x);
            ++num_lines;
        }
        num_lines = MAX(num_lines, 1);
    }

    const float avail = ImGui::GetContentRegionAvail().x;
    const bool need_hscroll = max_w + pad * 2.0f > avail;
    const float body_h = (float)num_lines * lh + pad + (need_hscroll ? style.ScrollbarSize : 0.0f);
    const float total_h = header_h + body_h;

    const ImVec2 p = ImGui::GetCursorScreenPos();
    if (!ImGui::IsRectVisible(ImVec2(avail, total_h))) {
        ImGui::Dummy(ImVec2(avail, total_h));
        return;
    }

    ImDrawList* dl = ImGui::GetWindowDrawList();
    dl->AddRectFilled(p, ImVec2(p.x + avail, p.y + total_h), ImGui::GetColorU32(ImGuiCol_FrameBg), 4.0f);
    dl->AddRectFilled(p, ImVec2(p.x + avail, p.y + header_h), ImGui::GetColorU32(ImGuiCol_TableHeaderBg), 4.0f, ImDrawFlags_RoundCornersTop);

    // Header: language on the left, buttons on the right
    ImGui::Dummy(ImVec2(avail, header_h));
    dl->AddText(ImVec2(p.x + pad, p.y + (header_h - fs) * 0.5f), ImGui::GetColorU32(ImGuiCol_TextDisabled),
                str_empty(lang) ? "code" : lang_beg, str_empty(lang) ? nullptr : lang_end);

    const int id = m_code_index;
    const bool feedback = (g.feedback_id == id) && ImGui::GetTime() < g.feedback_until;
    const char* copy_label = (feedback && g.feedback == Feedback_Copied)   ? "Copied"   : "Copy";
    const char* ins_label  = (feedback && g.feedback == Feedback_Inserted) ? "Inserted" : "Insert";
    ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(style.FramePadding.x, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
    const float bh = fs + 2.0f;
    float bw = ImGui::CalcTextSize(copy_label).x + style.FramePadding.x * 2.0f;
    if (script) bw += ImGui::CalcTextSize(ins_label).x + style.FramePadding.x * 2.0f + style.ItemSpacing.x;
    ImGui::SetCursorScreenPos(ImVec2(p.x + avail - bw - pad * 0.5f, p.y + (header_h - bh) * 0.5f));
    if (ImGui::Button(copy_label)) {
        ImGui::SetClipboardText(code_beg);  // zero terminated, see imgui_md
        set_feedback(id, Feedback_Copied);
    }
    if (script) {
        ImGui::SameLine();
        if (ImGui::Button(ins_label)) {
            // Handed to the host through draw_window(), as one or more complete lines
            md_allocator_i* heap = md_get_heap_allocator();
            md_array_shrink(g.insert_buf, 0);
            md_array_push_array(g.insert_buf, code.ptr, code.len, heap);
            if (code.len == 0 || code.ptr[code.len - 1] != '\n') md_array_push(g.insert_buf, '\n', heap);
            md_array_push(g.insert_buf, '\0', heap);
            g.insert_pending = true;
            set_feedback(id, Feedback_Inserted);
        }
    }
    ImGui::PopStyleColor();
    ImGui::PopStyleVar();

    // Body: scrolls horizontally when a line is too long
    ImGui::SetCursorScreenPos(ImVec2(p.x, p.y + header_h));
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(pad, pad * 0.5f));
    ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0, 0, 0, 0));
    if (ImGui::BeginChild("##code", ImVec2(avail, body_h), ImGuiChildFlags_AlwaysUseWindowPadding, ImGuiWindowFlags_HorizontalScrollbar)) {
        ImDrawList* cdl = ImGui::GetWindowDrawList();
        const ImVec2 origin = ImGui::GetCursorScreenPos();
        str_t rest = code, line;
        for (int l = 0; str_extract_line(&line, &rest); ++l) {
            draw_code_line(cdl, ImVec2(origin.x, origin.y + (float)l * lh), line, script, pal);
        }
        ImGui::Dummy(ImVec2(max_w, (float)num_lines * lh));
    }
    ImGui::EndChild();
    ImGui::PopStyleColor();
    ImGui::PopStyleVar();
}

void Renderer::on_heading(int, float screen_y) {
    if (g.target_heading == g.headings_rendered) {
        // The target stays armed for a few frames: text wraps differently once the scrollbar has appeared, so the
        // position of the heading can still move after the first frame.
        ImGui::SetScrollY(screen_y - ImGui::GetWindowPos().y + ImGui::GetScrollY() - ImGui::GetStyle().WindowPadding.y);
    }
    ++g.headings_rendered;
}

void Renderer::open_url(const char* url, const char* url_end) {
    follow_link(url, url_end);
}

void draw_toolbar() {
    State& s = g;
    const float fs = ImGui::GetFontSize();

    ImGui::BeginDisabled(s.history_pos <= 0);
    if (ImGui::ArrowButton("##back", ImGuiDir_Left)) history_step(-1);
    ImGui::EndDisabled();
    ImGui::SameLine();
    ImGui::BeginDisabled(s.history_pos + 1 >= s.history_len);
    if (ImGui::ArrowButton("##forward", ImGuiDir_Right)) history_step(+1);
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (ImGui::Button("Contents")) {
        s.query[0] = '\0';
        navigate(0);
    }
    ImGui::SameLine();

    if (s.focus_search) {
        ImGui::SetKeyboardFocusHere();
        s.focus_search = false;
    }
    ImGui::SetNextItemWidth(MAX(fs * 8.0f, ImGui::GetContentRegionAvail().x));
    const bool enter = ImGui::InputTextWithHint("##search", "Search the reference (procedure, keyword...)", s.query, sizeof(s.query),
                                                ImGuiInputTextFlags_EscapeClearsAll | ImGuiInputTextFlags_EnterReturnsTrue);
    if (enter) {
        update_search_results();
        if (md_array_size(s.hits) > 0) navigate(s.hits[0].section);
    }

    if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows) && ImGui::IsKeyPressed(ImGuiKey_F) && ImGui::GetIO().KeyCtrl) {
        s.focus_search = true;
    }
}

void draw_sidebar() {
    State& s = g;
    const float fs = ImGui::GetFontSize();
    const ImVec4 dim = ImGui::GetStyle().Colors[ImGuiCol_TextDisabled];

    if (s.query[0] != '\0') {
        update_search_results();
        const int num_hits = (int)md_array_size(s.hits);
        if (num_hits == 0) {
            ImGui::TextDisabled("No matches.");
            return;
        }
        ImGui::TextDisabled("%d match%s", num_hits, num_hits == 1 ? "" : "es");
        for (int i = 0; i < num_hits; ++i) {
            const int si = s.hits[i].section;
            const Section& sec = s.sections[si];
            const Page& page = s.pages[sec.page];
            ImGui::PushID(si);
            if (ImGui::Selectable(str_empty(sec.proc) ? sec.title.ptr : sec.proc.ptr, si == s.current)) navigate(si);
            if (!str_eq(page.title, sec.title)) {
                ImGui::SameLine();
                ImGui::TextColored(dim, "%s", page.title.ptr);
            }
            ImGui::PopID();
        }
        return;
    }

    for (int pi = 0; pi < (int)md_array_size(s.pages); ++pi) {
        const Page& page = s.pages[pi];
        const bool open = (pi == s.page);
        const int head = page.first_section;
        const char* label = pi == 0 ? "Overview" : page.title.ptr;
        ImGui::PushID(pi);
        ImVec4 page_col = ImGui::GetStyleColorVec4(open ? ImGuiCol_Text : ImGuiCol_TextDisabled);
        if (!open) page_col.w = MIN(1.0f, page_col.w * 1.6f);
        ImGui::PushStyleColor(ImGuiCol_Text, page_col);
        const bool clicked = ImGui::Selectable(label, open && s.current == head);
        ImGui::PopStyleColor();
        if (open && s.reveal_nav && s.current == head) ImGui::SetScrollHereY(0.3f);
        if (clicked) navigate(head);
        if (ImGui::IsItemHovered() && ImGui::CalcTextSize(label).x > ImGui::GetContentRegionAvail().x) ImGui::SetTooltip("%s", label);

        if (open) {
            ImGui::Indent(fs * 0.8f);
            for (int si = page.first_section; si < page.first_section + page.num_sections; ++si) {
                const Section& sec = s.sections[si];
                if (sec.level < 3) continue;
                ImGui::PushID(si);
                const bool is_proc = !str_empty(sec.proc);
                if (!is_proc) ImGui::PushStyleColor(ImGuiCol_Text, dim);
                if (ImGui::Selectable(sec.title.ptr, si == s.current)) navigate(si);
                if (!is_proc) ImGui::PopStyleColor();
                if (si == s.current && s.reveal_nav) ImGui::SetScrollHereY(0.5f);
                ImGui::PopID();
            }
            ImGui::Unindent(fs * 0.8f);
        }
        ImGui::PopID();
    }
    s.reveal_nav = false;
}

void draw_content() {
    State& s = g;
    if (s.page < 0 || s.page >= (int)md_array_size(s.pages)) return;
    const Page& page = s.pages[s.page];

    if (s.scroll_top) {
        ImGui::SetScrollY(0.0f);
        s.scroll_top = false;
    }

    s.headings_rendered = 0;
    s.renderer.pal = make_palette();
    s.renderer.print(s.doc.ptr + page.beg, s.doc.ptr + page.end);

    if (s.target_heading >= 0 && --s.scroll_settle <= 0) s.target_heading = -1;

    // Previous / next page
    ImGui::Dummy(ImVec2(0.0f, ImGui::GetFontSize() * 0.5f));
    ImGui::Separator();
    char label[256];
    if (s.page > 0) {
        const Page& prev = s.pages[s.page - 1];
        snprintf(label, sizeof(label), "<  %s", s.page == 1 ? "Overview" : prev.title.ptr);
        if (ImGui::Button(label)) navigate(prev.first_section);
        ImGui::SameLine();
    }
    if (s.page + 1 < (int)md_array_size(s.pages)) {
        const Page& next = s.pages[s.page + 1];
        snprintf(label, sizeof(label), "%s  >", next.title.ptr);
        if (ImGui::Button(label)) navigate(next.first_section);
    }
    ImGui::Dummy(ImVec2(0.0f, ImGui::GetFontSize()));
}

void draw() {
    State& s = g;
    if (md_array_size(s.sections) == 0) {
        ImGui::TextDisabled("The script reference is not available in this build.");
        return;
    }
    const float fs = ImGui::GetFontSize();
    if (s.nav_width <= 0.0f) s.nav_width = fs * 17.0f;

    draw_toolbar();
    ImGui::Spacing();

    const float total_w = ImGui::GetContentRegionAvail().x;
    const float splitter_w = 6.0f;
    s.nav_width = CLAMP(s.nav_width, fs * 8.0f, MAX(fs * 8.0f, total_w * 0.6f));

    if (ImGui::BeginChild("##ref_nav", ImVec2(s.nav_width, 0.0f), ImGuiChildFlags_Borders)) {
        draw_sidebar();
    }
    ImGui::EndChild();

    ImGui::SameLine(0.0f, 0.0f);
    ImGui::InvisibleButton("##ref_split", ImVec2(splitter_w, ImGui::GetContentRegionAvail().y));
    if (ImGui::IsItemActive()) s.nav_width += ImGui::GetIO().MouseDelta.x;
    if (ImGui::IsItemHovered() || ImGui::IsItemActive()) ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeEW);
    ImGui::SameLine(0.0f, 0.0f);

    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(fs * 0.9f, fs * 0.6f));
    if (ImGui::BeginChild("##ref_content", ImVec2(0.0f, 0.0f), ImGuiChildFlags_Borders | ImGuiChildFlags_AlwaysUseWindowPadding)) {
        draw_content();
    }
    ImGui::EndChild();
    ImGui::PopStyleVar();
}

}  // namespace

////////////////////////////////////////////////////////////////////////////////
// Interface
////////////////////////////////////////////////////////////////////////////////

bool show(str_t topic) {
    ensure_loaded();
    if (str_empty(topic)) return false;
    const int idx = find_anchor(topic);
    if (idx < 0) return false;
    g.query[0] = '\0';
    navigate(idx);
    return true;
}

void search(str_t query) {
    ensure_loaded();
    set_query(query);
    g.focus_search = false;
}

void focus_search() {
    g.focus_search = true;
}

Action draw_window(bool* p_open) {
    ensure_loaded();
    g.insert_pending = false;

    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 60.0f, ImGui::GetFontSize() * 40.0f), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Script Reference", p_open)) {
        draw();
    }
    ImGui::End();

    Action action;
    if (g.insert_pending) action.insert_code = {g.insert_buf, md_array_size(g.insert_buf) - 1};
    return action;
}

}  // namespace script_reference
