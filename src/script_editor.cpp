#include "script_editor.h"

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>

#include <imgui.h>
#include <imgui_internal.h>

#include <algorithm>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace script_editor {

////////////////////////////////////////////////////////////////////////////////
// Language and palette
////////////////////////////////////////////////////////////////////////////////

const TextEditor::Language* language() {
    static TextEditor::Language lang;
    static bool initialized = false;
    if (!initialized) {
        initialized = true;
        lang.name = "mdscript";
        lang.caseSensitive = true;
        lang.singleLineComment = "#";
        lang.hasSingleQuotedStrings = true;
        lang.hasDoubleQuotedStrings = true;
        lang.stringEscape = '\\';

        // Keywords and built-in names come from mdlib, so the colouring follows the language
        const str_t* keywords = md_script_keywords();
        for (size_t i = 0; i < md_script_num_keywords(); ++i) {
            lang.keywords.emplace(keywords[i].ptr, keywords[i].len);
        }
        const size_t num_identifiers = md_script_builtin_identifiers(nullptr, 0);
        std::vector<str_t> identifiers(num_identifiers);
        md_script_builtin_identifiers(identifiers.data(), identifiers.size());
        for (str_t id : identifiers) {
            lang.identifiers.emplace(id.ptr, id.len);
        }

        // mdscript identifiers, numbers and punctuation follow C, so reuse the editor's own C tokenizers
        const TextEditor::Language* c = TextEditor::Language::C();
        lang.isPunctuation = c->isPunctuation;
        lang.getIdentifier = c->getIdentifier;
        lang.getNumber     = c->getNumber;
    }
    return &lang;
}

const TextEditor::Palette& retro_blue_palette() {
    static const TextEditor::Palette palette = {{
        0xff00ffff,  // text
        0xffffff00,  // keyword
        0xffffff00,  // declaration
        0xff00ff00,  // number
        0xff808000,  // string
        0xffffffff,  // punctuation
        0xff008000,  // preprocessor
        0xff00ffff,  // identifier
        0xffffffff,  // known identifier
        0xff808080,  // comment
        0xff800000,  // background
        0xff0080ff,  // cursor
        0x80ffff00,  // selection
        0x40ffffff,  // whitespace
        0x40ffffff,  // matchingBracketBackground
        0x80ffffff,  // matchingBracketActive
        0xff00ffff,  // matchingBracketLevel1
        0xffffff00,  // matchingBracketLevel2
        0xffff00ff,  // matchingBracketLevel3
        0xa00000ff,  // matchingBracketError
        0xff808000,  // line number
        0xffffff00,  // current line number
        0x40000000,  // current line highlight
        0x40000000,  // current line highlight border
    }};
    return palette;
}

////////////////////////////////////////////////////////////////////////////////
// Helpers
////////////////////////////////////////////////////////////////////////////////

static bool is_identifier(const std::string& s) {
    if (s.empty() || !(is_alpha(s[0]) || s[0] == '_')) return false;
    for (char c : s) {
        if (!(is_alpha(c) || is_digit(c) || c == '_')) return false;
    }
    return true;
}

std::string word_at_cursor(const TextEditor& editor) {
    const TextEditor::DocPos pos = editor.GetMainCursorPosition();
    auto word_at = [&editor](TextEditor::DocPos p) {
        return editor.GetSectionText(editor.FindWordStart(p), editor.FindWordEnd(p));
    };
    std::string word = word_at(pos);
    if (!is_identifier(word) && pos.index > 0) {
        // The cursor is typically right after the word, e.g. "distance|("
        word = word_at(TextEditor::DocPos(pos.line, pos.index - 1));
    }
    return is_identifier(word) ? word : std::string();
}

std::string word_at_mouse(const TextEditor& editor, ImVec2 mouse_pos) {
    if (!editor.IsMousePosOverGlyph(mouse_pos)) return {};
    const TextEditor::DocPos pos = editor.GetDocPosAtMousePos(mouse_pos);
    std::string word = editor.GetSectionText(editor.FindWordStart(pos), editor.FindWordEnd(pos));
    return is_identifier(word) ? word : std::string();
}

void append_line(TextEditor& editor, str_t line) {
    const TextEditor::DocPos end(editor.GetLineCount() - 1, SIZE_MAX);  // clamped to the end of the last line
    std::string text;
    if (!editor.IsEmpty()) text += '\n';
    text.append(line.ptr, line.len);
    editor.ReplaceSectionText(end, end, text);
}

void insert_lines_at_cursor(TextEditor& editor, str_t code) {
    std::string text;
    if (editor.GetCurrentCursorPosition().index > 0) text += '\n';
    text.append(code.ptr, code.len);
    editor.ReplaceTextInCurrentCursor(text);
}

bool has_focus_after_render() {
    const ImGuiContext& g = *ImGui::GetCurrentContext();
    const ImGuiID child_id = g.LastItemData.ID;
    // The focused window is the editor's child window, or a window inside it (such as its find / replace box)
    for (const ImGuiWindow* w = g.NavWindow; w; w = w->ParentWindow) {
        if (w->ChildId == child_id) return true;
        if (!(w->Flags & ImGuiWindowFlags_ChildWindow)) break;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////
// Autocomplete
////////////////////////////////////////////////////////////////////////////////

static inline bool ascii_lower(char c) { return 'a' <= c && c <= 'z'; }
static inline bool ascii_upper(char c) { return 'A' <= c && c <= 'Z'; }
static inline char to_lower(char c)    { return ascii_upper(c) ? (char)(c - 'A' + 'a') : c; }

// Whether word[j] starts the word or a part of it: after an underscore or at a lower to upper case change
static bool starts_part(std::string_view word, size_t j) {
    if (j == 0) return true;
    const char prev = word[j - 1], cur = word[j];
    return (prev == '_' && cur != '_') || (ascii_lower(prev) && ascii_upper(cur));
}

// Fuzzy match in the style of VS Code. The pattern has to appear in order in the word, ignoring case, and its first
// character has to start the word or a part of it. Every matched character scores, more if the case is the same, and
// gets a bonus for following the previous match directly or for starting a part. Characters skipped between matches
// cost a little. The best scoring alignment is found with dynamic programming (O(pattern * word), both are short).
// Returns false if there is no match.
static bool fuzzy_match(std::string_view pat, std::string_view word, int* out_score) {
    enum : int { MATCH = 1, SAME_CASE = 1, START = 8, PART = 6, CONSECUTIVE = 5, GAP = 1, NONE = -1000000 };
    const size_t n = pat.size(), m = word.size();
    if (n == 0) { *out_score = 0; return true; }
    if (n > m) return false;

    // Cheap rejection before the real thing: is the pattern a subsequence at all?
    size_t k = 0;
    for (size_t j = 0; j < m && k < n; ++j) {
        if (to_lower(word[j]) == to_lower(pat[k])) ++k;
    }
    if (k < n) return false;

    // prev[j] / cur[j]: best score with the pattern matched up to the previous / current character, which is at word[j]
    static thread_local std::vector<int> prev, cur;
    prev.assign(m, NONE);
    cur.assign(m, NONE);

    for (size_t i = 0; i < n; ++i) {
        int run = NONE;  // best prev[k] - GAP * (j - 1 - k) over k <= j - 2, i.e. reaching word[j] across a gap
        for (size_t j = 0; j < m; ++j) {
            if (j >= 2) run = std::max(run, prev[j - 2]) - GAP;
            cur[j] = NONE;
            if (to_lower(word[j]) != to_lower(pat[i])) continue;

            const bool part = starts_part(word, j);
            const int  s    = MATCH + (word[j] == pat[i] ? SAME_CASE : 0);
            if (i == 0) {
                if (part) cur[j] = s + (j == 0 ? START : PART);
            } else {
                const int after_prev = (j >= 1 && prev[j - 1] > NONE / 2) ? prev[j - 1] + CONSECUTIVE : NONE;
                const int after_gap  = run > NONE / 2 ? run + (part ? PART : 0) : NONE;
                const int best = std::max(after_prev, after_gap);
                if (best > NONE / 2) cur[j] = best + s;
            }
        }
        std::swap(prev, cur);
    }

    int best = NONE;
    for (size_t j = 0; j < m; ++j) best = std::max(best, prev[j]);
    if (best <= NONE / 2) return false;
    *out_score = best;
    return true;
}

// Where a suggestion comes from, which is also the order among equally good matches
enum Origin {
    Origin_Document,  // variables (and anything else) named in the script
    Origin_Builtin,   // built-in procedures and constants
    Origin_Keyword,
};

struct Vocabulary {
    std::vector<std::string> builtins;
    std::vector<std::string> keywords;
    std::unordered_set<std::string> known;  // all of the above
};

static const Vocabulary& vocabulary() {
    static const Vocabulary vocab = [] {
        Vocabulary v;
        const size_t num = md_script_builtin_identifiers(nullptr, 0);
        std::vector<str_t> names(num);
        md_script_builtin_identifiers(names.data(), names.size());
        for (str_t name : names) {
            std::string s(name.ptr, name.len);
            if (is_identifier(s) && v.known.insert(s).second) v.builtins.push_back(s);  // skip operators
        }
        const str_t* keywords = md_script_keywords();
        for (size_t i = 0; i < md_script_num_keywords(); ++i) {
            std::string s(keywords[i].ptr, keywords[i].len);
            if (v.known.insert(s).second) v.keywords.push_back(s);
        }
        return v;
    }();
    return vocab;
}

static void suggest(const TextEditor& editor, TextEditor::AutoCompleteState& state) {
    const size_t MAX_SUGGESTIONS = 100;
    state.suggestions.clear();

    const std::string& term = state.searchTerm;
    if (state.inNumber || (!term.empty() && is_digit(term[0]))) return;

    const Vocabulary& vocab = vocabulary();

    // Identifiers in the document, counted, so that the word being completed can be left out unless it also appears
    // elsewhere. The word is what a suggestion would replace, which extends past the cursor.
    std::unordered_map<std::string, int> in_document;
    editor.IterateIdentifiers([&](const std::string& id) {
        if (is_identifier(id) && !vocab.known.count(id)) in_document[id] += 1;
    });
    const std::string word = editor.GetSectionText(state.searchTermStart, editor.FindWordEnd(state.searchTermStart, true));
    if (auto it = in_document.find(word); it != in_document.end() && --it->second == 0) {
        in_document.erase(it);
    }

    struct Match {
        std::string_view name;
        Origin origin;
        int score;
    };
    std::vector<Match> matches;
    auto consider = [&](std::string_view name, Origin origin) {
        int score;
        if (fuzzy_match(term, name, &score)) matches.push_back({name, origin, score});
    };
    for (const auto& [name, count] : in_document) consider(name, Origin_Document);
    for (const std::string& name : vocab.builtins) consider(name, Origin_Builtin);
    for (const std::string& name : vocab.keywords) consider(name, Origin_Keyword);

    // What has been typed already is complete, nothing to offer
    if (matches.size() == 1 && matches[0].name == term) return;

    // Best score first. Among equals: an exact match, then by origin, then the shortest (the one needing the least
    // typing to tell apart), then alphabetically. Without a term everything scores 0 and is listed alphabetically.
    std::sort(matches.begin(), matches.end(), [&term](const Match& a, const Match& b) {
        if (a.score != b.score) return a.score > b.score;
        const bool a_exact = a.name == term, b_exact = b.name == term;
        if (a_exact != b_exact) return a_exact;
        if (a.origin != b.origin) return a.origin < b.origin;
        if (!term.empty() && a.name.size() != b.name.size()) return a.name.size() < b.name.size();
        return a.name < b.name;
    });

    const size_t count = std::min(matches.size(), MAX_SUGGESTIONS);
    state.suggestions.reserve(count);
    for (size_t i = 0; i < count; ++i) state.suggestions.emplace_back(matches[i].name);
}

void enable_autocomplete(TextEditor& editor) {
    TextEditor::AutoCompleteConfig config;
    config.callback = [ed = &editor](TextEditor::AutoCompleteState& state) { suggest(*ed, state); };
    editor.SetAutoCompleteConfig(&config);
}

////////////////////////////////////////////////////////////////////////////////
// Markers
////////////////////////////////////////////////////////////////////////////////

// Squiggle types, so that the editor can tell ours apart (we only ever clear all of them)
enum SquiggleType : size_t {
    SquiggleType_Error = 1,
    SquiggleType_Warning,
    SquiggleType_Hovered,
};

static const ImU32 SQUIGGLE_COLOR_ERROR   = IM_COL32(255,  64,  64, 255);
static const ImU32 SQUIGGLE_COLOR_WARNING = IM_COL32(255, 210,   0, 255);
static const ImU32 SQUIGGLE_COLOR_HOVERED = IM_COL32(255, 255, 255, 255);

void markers_clear(Markers* m) {
    ASSERT(m);
    if (md_array_size(m->list) > 0) m->dirty = true;
    if (m->arena) md_arena_allocator_reset(m->arena);
    m->list = nullptr;
    m->line_offsets = nullptr;
    m->source = {};
    m->hovered = -1;
}

void markers_set_source(Markers* m, str_t source) {
    ASSERT(m);
    markers_clear(m);
    if (!m->arena) m->arena = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    m->source = source;
    md_array_push(m->line_offsets, 0u, m->arena);
    for (size_t i = 0; i < source.len; ++i) {
        if (source.ptr[i] == '\n') md_array_push(m->line_offsets, (uint32_t)(i + 1), m->arena);
    }
}

// Document position of a byte offset into the source: its line, and the number of code points before it on that line
static TextEditor::DocPos doc_pos(const Markers* m, int offset) {
    const uint32_t* lines = m->line_offsets;
    const size_t num_lines = md_array_size(lines);
    const uint32_t off = (uint32_t)CLAMP(offset, 0, (int)m->source.len);
    // Last line that starts at or before the offset
    size_t lo = 0, hi = num_lines;
    while (hi - lo > 1) {
        const size_t mid = (lo + hi) / 2;
        if (lines[mid] <= off) lo = mid; else hi = mid;
    }
    size_t index = 0;
    for (uint32_t i = lines[lo]; i < off; ++i) {
        if (((unsigned char)m->source.ptr[i] & 0xC0) != 0x80) ++index;  // count UTF-8 lead bytes
    }
    return TextEditor::DocPos(lo, index);
}

void markers_add(Markers* m, MarkerType type, int prio, md_script_range_marker_t range, str_t text,
                 const md_bitfield_t* atoms, const md_script_vis_payload_o* payload) {
    ASSERT(m);
    ASSERT(m->arena && "markers_set_source() must be called before markers are added");
    Marker marker = {};
    marker.type = type;
    marker.prio = prio;
    marker.beg  = doc_pos(m, range.beg);
    marker.end  = doc_pos(m, range.end);
    marker.text = str_copy(text, m->arena);
    if (atoms) {
        md_bitfield_t* copy = (md_bitfield_t*)md_alloc(m->arena, sizeof(md_bitfield_t));
        md_bitfield_init(copy, m->arena);
        md_bitfield_copy(copy, atoms);
        marker.atoms = copy;
    }
    marker.payload = payload;
    md_array_push(m->list, marker, m->arena);
    m->dirty = true;
}

const Marker* markers_update(Markers* m, TextEditor* editor, bool editor_hovered) {
    ASSERT(m);
    ASSERT(editor);
    const int num_markers = (int)md_array_size(m->list);

    m->hovered = -1;
    const ImVec2 mouse = ImGui::GetMousePos();
    if (num_markers > 0 && editor_hovered && !ImGui::IsAnyItemActive() && editor->IsMousePosOverGlyph(mouse)) {
        const TextEditor::DocPos pos = editor->GetDocPosAtMousePos(mouse);
        for (int i = 0; i < num_markers; ++i) {
            const Marker& marker = m->list[i];
            if (marker.beg <= pos && pos < marker.end && (m->hovered == -1 || marker.prio > m->list[m->hovered].prio)) {
                m->hovered = i;
            }
        }
    }

    // A glyph carries at most one squiggle, so all of them are laid down again whenever anything changes. The hovered
    // one goes last, on top.
    if (m->dirty || m->hovered != m->shown_hovered) {
        editor->ClearSquiggles();
        for (int i = 0; i < num_markers; ++i) {
            const Marker& marker = m->list[i];
            if (marker.type == MarkerType_Error) {
                editor->AddSquiggle(marker.beg, marker.end, SquiggleType_Error, SQUIGGLE_COLOR_ERROR);
            } else if (marker.type == MarkerType_Warning) {
                editor->AddSquiggle(marker.beg, marker.end, SquiggleType_Warning, SQUIGGLE_COLOR_WARNING);
            }
        }
        if (m->hovered != -1 && m->list[m->hovered].type == MarkerType_Visualization) {
            editor->AddSquiggle(m->list[m->hovered].beg, m->list[m->hovered].end, SquiggleType_Hovered, SQUIGGLE_COLOR_HOVERED);
        }
        m->dirty = false;
        m->shown_hovered = m->hovered;
    }

    return markers_hovered(m);
}

const Marker* markers_hovered(const Markers* m) {
    ASSERT(m);
    return m->hovered != -1 ? &m->list[m->hovered] : nullptr;
}

}  // namespace script_editor
