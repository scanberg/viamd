#include "script_editor.h"

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>

#include <imgui.h>
#include <imgui_internal.h>

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
