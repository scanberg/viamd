#pragma once

// viamd's side of the script editor.
//
// The editor widget (ext/ImGuiColorTextEdit, https://github.com/goossens/ImGuiColorTextEdit) is vendored unmodified so
// that it can be updated by copying in a newer version. Everything viamd needs on top of it lives here: the mdscript
// language definition, autocompletion, markers for compiler errors, warnings and visualization tokens, and a few helpers.
//
// Markers are kept here rather than in the editor. The editor is only asked two things: which glyph the mouse is over
// (TextEditor::GetDocPosAtMousePos) and to draw squiggles. Errors and warnings are underlined in red and yellow, and
// the visualization token under the mouse in white.

#include <TextEditor.h>

#include <core/md_str.h>
#include <core/md_bitfield.h>
#include <md_script.h>

#include <string>

struct md_allocator_i;

namespace script_editor {

// Syntax colouring for mdscript
const TextEditor::Language* language();

// The palette of the previous editor version, which the fork does not provide (it has dark and light)
const TextEditor::Palette& retro_blue_palette();

// Identifier under (or just before) the main cursor, empty if there is none
std::string word_at_cursor(const TextEditor& editor);

// Identifier under the mouse position (in screen coordinates), empty if there is none
std::string word_at_mouse(const TextEditor& editor, ImVec2 mouse_pos);

// Appends text as a new line at the end of the document. The cursor moves to the end of it.
void append_line(TextEditor& editor, str_t line);

// Inserts code at the current cursor as whole lines: a line break is added in front if the cursor is not at the start
// of a line. Does nothing in read-only mode.
void insert_lines_at_cursor(TextEditor& editor, str_t code);

// True if the editor has keyboard focus. Call right after TextEditor::Render(), while its child window is the last item.
bool has_focus_after_render();

// Turns on autocompletion for mdscript. Suggestions are the keywords, the built-in procedures and constants, and the
// identifiers used in the document, ranked by how well they match what has been typed. The characters typed have to
// appear in order in a suggestion (ignoring case) and the first one has to start it or one of its parts, so 'sw' finds
// shape_weights. Suggestions pop up while typing an identifier, Ctrl+Space asks for them (also on macOS); Tab or Enter
// inserts the selected one, Escape closes the list.
// The editor keeps a reference to itself in its configuration, so it must not move while autocomplete is on.
void enable_autocomplete(TextEditor& editor);

enum MarkerType {
    MarkerType_Error,
    MarkerType_Warning,
    MarkerType_Visualization,
};

struct Marker {
    MarkerType type;
    int prio;                               // where markers overlap, the one with the highest priority is hovered
    TextEditor::DocPos beg, end;            // [beg, end)
    str_t text;                             // tooltip, zero terminated
    const md_bitfield_t* atoms;             // errors and warnings: the atoms the message is about (a copy), may be NULL
    const md_script_vis_payload_o* payload; // visualization tokens (owned by the script IR the markers came from)
};

struct Markers {
    md_allocator_i* arena = nullptr;    // everything below, reset when the markers are cleared
    Marker*   list = nullptr;           // md_array
    uint32_t* line_offsets = nullptr;   // md_array, byte offset of every line in the source
    str_t     source = {};              // the text that the byte ranges given to markers_add() refer to
    int  hovered = -1;
    int  shown_hovered = -1;            // hovered marker that the squiggles currently show
    bool dirty = false;                 // squiggles do not match the list
};

// Removes all markers. Their squiggles disappear the next time markers_update() runs.
void markers_clear(Markers* markers);

// Clears the markers and sets the source text that the byte ranges of the next markers_add() calls refer to. This is
// the text that was compiled, i.e. what TextEditor::GetText() returned. It is only read by markers_add(), so it has to
// stay alive until the markers have been added.
void markers_set_source(Markers* markers, str_t source);

// Adds a marker for the byte range [beg, end) of the source. text and atoms are copied.
void markers_add(Markers* markers, MarkerType type, int prio, md_script_range_marker_t range, str_t text,
                 const md_bitfield_t* atoms = nullptr, const md_script_vis_payload_o* payload = nullptr);

// Call after TextEditor::Render(). Finds the marker under the mouse and updates the squiggles.
// editor_hovered: whether the editor itself is hovered (ImGui::IsItemHovered() right after Render()).
// Returns the hovered marker, or NULL.
const Marker* markers_update(Markers* markers, TextEditor* editor, bool editor_hovered);

// The marker found by the last markers_update(), or NULL
const Marker* markers_hovered(const Markers* markers);

}  // namespace script_editor
