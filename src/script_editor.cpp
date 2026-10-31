#include "script_editor.h"

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>

#include <imgui.h>
#include <imgui_internal.h>

#include <algorithm>
#include <string_view>
#include <vector>

#include <stdio.h>

namespace script_editor {

////////////////////////////////////////////////////////////////////////////////
// Language and palette
////////////////////////////////////////////////////////////////////////////////

const TextEditor::Language* language(script::Language script_lang) {
    // One for each script language, made the first time it is asked for
    static TextEditor::Language langs[3];
    static bool initialized[3] = {};
    const size_t idx = (size_t)script_lang < 3 ? (size_t)script_lang : 0;
    TextEditor::Language& lang = langs[idx];
    if (!initialized[idx]) {
        initialized[idx] = true;
        lang.name = script_lang == script::Language::V2 ? "VIAL" : "mdscript";
        lang.caseSensitive = true;
        lang.singleLineComment = "#";
        lang.hasSingleQuotedStrings = true;
        lang.hasDoubleQuotedStrings = true;
        lang.stringEscape = '\\';

        // Keywords and built-in names come from the script, so the colouring follows the language
        const str_t* keywords = script::keywords(script_lang);
        for (size_t i = 0; i < script::num_keywords(script_lang); ++i) {
            lang.keywords.emplace(keywords[i].ptr, keywords[i].len);
        }
        const size_t num_identifiers = script::builtin_names(script_lang, nullptr, 0);
        std::vector<str_t> identifiers(num_identifiers);
        script::builtin_names(script_lang, identifiers.data(), identifiers.size());
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
// Dropping text in
////////////////////////////////////////////////////////////////////////////////

bool attr_expression(char* buf, size_t cap, script::Language lang, const md_attributes_t* attributes, const md_attribute_t* attr) {
    ASSERT(buf && cap > 0);
    buf[0] = '\0';
    const str_t path = script::attr_path(lang, attributes, attr);
    if (str_empty(path)) return false;
    // Nothing in it that the script's lexer ends (or gives up on) a string at
    for (size_t i = 0; i < path.len; ++i) {
        const char c = path.ptr[i];
        if (c == '"' || c == ';' || c == '\n' || c == '\r') return false;
    }
    const int len = snprintf(buf, cap, "attr(\"" STR_FMT "\")", STR_ARG(path));
    if (len < 0 || (size_t)len >= cap) {
        buf[0] = '\0';
        return false;
    }
    return true;
}

bool accept_drop(TextEditor& editor, const char* payload_type, str_t text) {
    if (editor.IsReadOnlyEnabled() || str_empty(text)) return false;

    const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(payload_type, ImGuiDragDropFlags_AcceptBeforeDelivery);
    if (!payload) return false;

    const TextEditor::DocPos pos = editor.GetCursorPosAtMousePos(ImGui::GetMousePos());
    const std::string_view str(text.ptr, text.len);

    if (!payload->IsDelivery()) {
        // A caret where it goes, over the editor (its child window is drawn over the window that the target is in)
        const ImVec2 p = editor.GetScreenPosOfDocPos(pos);
        const float w = std::max(1.0f, ImGui::GetFontSize() * (1.0f / 16.0f));  // half its width
        ImDrawList* dl = ImGui::GetForegroundDrawList();
        dl->PushClipRect(ImGui::GetItemRectMin(), ImGui::GetItemRectMax(), true);
        dl->AddRectFilled(ImVec2(p.x - w, p.y), ImVec2(p.x + w, p.y + editor.GetLineHeight()), editor.GetPalette().get(TextEditor::Color::cursor));
        dl->PopClipRect();

        // And what goes in, in place of the source's preview
        ImGui::BeginTooltip();
        ImGui::TextUnformatted(str.data(), str.data() + str.size());
        ImGui::EndTooltip();
        return false;
    }

    // One edit, after which the cursor is at its end
    editor.ReplaceSectionText(pos, pos, str);
    editor.SetFocus();
    return true;
}

////////////////////////////////////////////////////////////////////////////////
// Autocomplete
////////////////////////////////////////////////////////////////////////////////

static inline bool ascii_lower(char c) { return 'a' <= c && c <= 'z'; }
static inline bool ascii_upper(char c) { return 'A' <= c && c <= 'Z'; }
static inline char to_lower(char c)    { return ascii_upper(c) ? (char)(c - 'A' + 'a') : c; }

static inline bool ascii_alnum(char c) { return ascii_lower(c) || ascii_upper(c) || ('0' <= c && c <= '9'); }

// Whether word[j] starts the word or a part of it: after a separator (shape_weights, edr/potential) or at a lower to
// upper case change (loadFile)
static bool starts_part(std::string_view word, size_t j) {
    if (j == 0) return true;
    const char prev = word[j - 1], cur = word[j];
    return (!ascii_alnum(prev) && ascii_alnum(cur)) || (ascii_lower(prev) && ascii_upper(cur));
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

// The order among equally good matches
static int kind_order(script::CompletionKind kind) {
    switch (kind) {
    case script::CompletionKind::Value:     return 0;
    case script::CompletionKind::Parameter: return 1;
    case script::CompletionKind::Variable:  return 2;
    case script::CompletionKind::Procedure: return 3;
    case script::CompletionKind::Constant:  return 4;
    default:                                return 5;
    }
}

bool completion_match(std::string_view typed, const script::Completion& item, int* score) {
    return fuzzy_match(typed, std::string_view(item.label.ptr, item.label.len), score);
}

bool completion_before(const script::Completion& a, int a_score, const script::Completion& b, int b_score, std::string_view typed) {
    const std::string_view a_label(a.label.ptr, a.label.len), b_label(b.label.ptr, b.label.len);
    if (a.fits != b.fits) return a.fits;
    const bool a_exact = a_label == typed, b_exact = b_label == typed;
    if (a_exact != b_exact) return a_exact;
    if (a_score != b_score) return a_score > b_score;
    const int a_kind = kind_order(a.kind), b_kind = kind_order(b.kind);
    if (a_kind != b_kind) return a_kind < b_kind;
    if (!typed.empty() && a_label.size() != b_label.size()) return a_label.size() < b_label.size();
    return a_label < b_label;
}

// What can be written at the cursor comes from the script (script::complete), which knows the language and the system.
// What matches what has been typed, and in which order, is decided here.
static void suggest(const TextEditor& editor, script::Language lang, const md_system_t* sys, TextEditor::AutoCompleteState& state) {
    const size_t MAX_SUGGESTIONS = 100;
    state.suggestions.clear();

    const std::string src = editor.GetText();
    auto byte_offset = [&editor](TextEditor::DocPos pos) { return editor.GetSectionText(TextEditor::DocPos(0, 0), pos).size(); };
    const int cursor = (int)byte_offset(state.searchTermEnd);

    md_temp_scope_t temp = md_temp_begin();
    const script::Completions completions = script::complete(lang, {src.data(), src.size()}, cursor, sys, md_temp_allocator(temp));

    // mdlib says which bytes a completion replaces: an identifier, or the contents of a string, such as an attribute
    // path. The editor replaces the word the cursor is in, which can be less (the 'potential' of 'edr/potential'), so
    // a completion is cut down to what goes in its place. The rest of mdlib's range has to be in the completion
    // already, as it is when that is what has been typed.
    const TextEditor::DocPos word_end = editor.FindWordEnd(state.searchTermStart, true);
    const size_t word_beg_byte = byte_offset(state.searchTermStart);
    const size_t word_end_byte = byte_offset(word_end);
    const size_t range_beg = (size_t)completions.range.beg, range_end = (size_t)completions.range.end;
    const bool word_in_range = range_beg <= word_beg_byte && word_end_byte <= range_end && range_end <= src.size();
    const std::string_view lead = word_in_range ? std::string_view(src).substr(range_beg, word_beg_byte - range_beg) : std::string_view();
    const std::string_view tail = word_in_range ? std::string_view(src).substr(word_end_byte, range_end - word_end_byte) : std::string_view();
    const std::string_view word = std::string_view(src).substr(std::min(word_beg_byte, src.size()), word_end_byte - std::min(word_beg_byte, word_end_byte));

    // Matched against all that has been typed of the value, not just the word
    const std::string_view typed(completions.prefix.ptr, completions.prefix.len);

    struct Match {
        const script::Completion* item;
        std::string_view insert;    // what replaces the word
        int score;
    };
    std::vector<Match> matches;
    for (size_t i = 0; word_in_range && i < completions.count; ++i) {
        const script::Completion& item = completions.items[i];
        const std::string_view text(item.text.ptr, item.text.len);
        if (text.size() < lead.size() + tail.size() || text.substr(0, lead.size()) != lead || text.substr(text.size() - tail.size()) != tail) continue;
        int score;
        if (completion_match(typed, item, &score)) {
            matches.push_back({&item, text.substr(lead.size(), text.size() - lead.size() - tail.size()), score});
        }
    }

    // Nothing to offer if what has been typed is already complete
    if (!(matches.size() == 1 && matches[0].insert == word)) {
        // Without a term everything scores 0, and what fits the argument comes first, then by kind and name
        std::sort(matches.begin(), matches.end(), [&](const Match& a, const Match& b) {
            return completion_before(*a.item, a.score, *b.item, b.score, typed);
        });

        const size_t count = std::min(matches.size(), MAX_SUGGESTIONS);
        state.suggestions.reserve(count);
        for (size_t i = 0; i < count; ++i) state.suggestions.emplace_back(matches[i].insert);
    }

    md_temp_end(temp);
}

void enable_autocomplete(TextEditor& editor, script::Language lang, const md_system_t* sys) {
    TextEditor::AutoCompleteConfig config;
    // Strings are for mdlib to judge: it offers the names a selector such as resname() takes, and nothing elsewhere
    config.triggerInStrings = true;
    config.callback = [ed = &editor, lang, sys](TextEditor::AutoCompleteState& state) { suggest(*ed, lang, sys, state); };
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

void markers_add(Markers* m, MarkerType type, int prio, script::Range range, str_t text,
                 const md_bitset_t* atoms, script::Ref vis_ref) {
    ASSERT(m);
    ASSERT(m->arena && "markers_set_source() must be called before markers are added");
    Marker marker = {};
    marker.type = type;
    marker.prio = prio;
    marker.beg  = doc_pos(m, range.beg);
    marker.end  = doc_pos(m, range.end);
    marker.text = str_copy(text, m->arena);
    if (atoms) {
        md_bitset_t* copy = (md_bitset_t*)md_alloc(m->arena, sizeof(md_bitset_t));
        *copy = md_bitset_copy(*atoms, m->arena);
        marker.atoms = copy;
    }
    marker.vis_ref = vis_ref;
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
