#pragma once

// viamd's side of the script editor.
//
// The editor widget (ext/ImGuiColorTextEdit, https://github.com/goossens/ImGuiColorTextEdit) is vendored so that it can
// be updated by copying in a newer version, with two additions to carry over when that is done: line decorations can go
// left of the line numbers (SetDecorationsBeforeLineNumbers), where the marks of the properties are, and what a drop
// needs to know of where the mouse is (GetCursorPosAtMousePos, GetScreenPosOfDocPos). Everything viamd needs on top of
// it lives here: the mdscript language definition, autocompletion, markers for compiler errors, warnings and
// visualization tokens, dropping text in, and a few helpers.
//
// Markers are kept here rather than in the editor. The editor is only asked two things: which glyph the mouse is over
// (TextEditor::GetDocPosAtMousePos) and to draw squiggles. Errors and warnings are underlined in red and yellow, and
// the visualization token under the mouse in white.

#include <TextEditor.h>

#include <core/md_str.h>
#include <core/md_bitset.h>
#include <md_script.h>

#include <string>
#include <string_view>

struct md_allocator_i;
struct md_attribute_t;
struct md_attributes_t;

// ImGui drag and drop payload type of text for the script, zero terminated, which the editor takes as it is. What an
// attribute that the script reads is dragged as where it carries no series (the Timelines and Distributions windows
// cannot show it, such as a frame axis).
#define SCRIPT_TEXT_DND "SCRIPT_TEXT"

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

// The expression that reads an attribute of the system (attributes: the table attr is in) in the script: attr("<path>"),
// with the path in the form that completion offers it (md_script_attr_path), "edr/potential" rather than
// "run/a/edr/potential" where that names it alone. Writes a terminated string of at most cap - 1 characters. False, with
// buf empty, if attr() cannot read the attribute or the expression does not fit.
bool attr_expression(char* buf, size_t cap, const md_attributes_t* attributes, const md_attribute_t* attr);

// Accepts the payload being dragged over the editor as text: it goes where a left click at the mouse would put the
// cursor. Call between ImGui::BeginDragDropTarget() and EndDragDropTarget(), right after TextEditor::Render() so that
// its child window is the target, with the type of the payload and the text it stands for. While it is dragged, a
// caret shows where the text would go and the tooltip what it is. Once dropped, the text is inserted as one edit
// (undone as one) with the cursor after it, and the editor takes the focus: true on that frame. Not in read-only mode.
bool accept_drop(TextEditor& editor, const char* payload_type, str_t text);

// Turns on autocompletion for mdscript. What can be written at the cursor comes from mdlib (md_script_complete):
// variables, procedures, constants and keywords, the parameters of the procedure being called, and the values that
// selectors take, such as the residue names of the system in resname("|") or the attribute paths in attr("|").
// Suggestions are ranked by how well they match what has been typed: the characters typed have to appear in order in
// a suggestion (ignoring case) and the first one has to start it or one of its parts, so 'sw' finds shape_weights and
// 'pot' edr/potential. They pop up while typing and Ctrl+Space asks for them (also on macOS); Tab or Enter inserts the
// selected one, Escape closes the list.
// sys is where the values come from, it is read when suggestions are made and may be NULL. The editor keeps a
// reference to itself in its configuration, so neither it nor sys may move while autocomplete is on.
void enable_autocomplete(TextEditor& editor, const md_system_t* sys);

// How completions (md_script_complete, md_filter_complete) are matched against what has been typed of them and ranked,
// here and in the inputs of filters (filter_input.h). Whether the label of item matches typed, and how well: the
// characters typed have to appear in it in order (ignoring case), the first one starting it or one of its parts.
bool completion_match(std::string_view typed, const md_script_completion_t& item, int* score);
// Whether a goes before b, both matching typed with those scores: what has the type expected at the cursor (fits)
// first, then an exact match, the better match, by kind (values, parameters, variables, procedures, constants,
// keywords), the shorter (the one needing the least typing to tell apart) and by name.
bool completion_before(const md_script_completion_t& a, int a_score, const md_script_completion_t& b, int b_score, std::string_view typed);

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
    const md_bitset_t* atoms;             // errors and warnings: the atoms the message is about (a copy), may be NULL
    md_script_vis_ref_t vis_ref;            // visualization tokens: resolves in the script IR the markers came from
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
                 const md_bitset_t* atoms = nullptr, md_script_vis_ref_t vis_ref = {});

// Call after TextEditor::Render(). Finds the marker under the mouse and updates the squiggles.
// editor_hovered: whether the editor itself is hovered (ImGui::IsItemHovered() right after Render()).
// Returns the hovered marker, or NULL.
const Marker* markers_update(Markers* markers, TextEditor* editor, bool editor_hovered);

// The marker found by the last markers_update(), or NULL
const Marker* markers_hovered(const Markers* markers);

}  // namespace script_editor
