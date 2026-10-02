#pragma once

// In-app viewer for the mdlib script reference (ext/mdlib/docs/script_reference.md).
//
// The reference is plain Markdown that lives next to the implementation and is baked into the executable. The viewer
// parses its outline once (headings, procedure metadata comments, alias anchors), shows one section ("page") at a
// time with a navigable sidebar, history and search, and renders the Markdown with imgui_md.
//
// The document is expected to follow the conventions listed at the top of script_reference.md:
//   - `##` headings are the pages, `###` headings are entries within a page.
//   - Every procedure has a `### name` heading followed by `<!-- proc name=... aliases=a,b category=... -->`
//     and optionally `<a id="alias"></a>` anchors.
//   - Code blocks tagged `mdscript` are examples that can be copied or inserted into the script editor.

#include <core/md_str.h>

namespace script_reference {

// Shows the entry for a procedure name, alias or heading anchor (case-insensitive).
// Returns false if there is no such entry, in which case nothing changes.
bool show(str_t topic);

// Puts the query in the search box and lists the matches in the sidebar.
void search(str_t query);

// Moves keyboard focus to the search box the next time the window is drawn.
void focus_search();

// Makes the window visible the next time it is drawn: brought to the front, or its tab selected if it is docked.
// With take_focus it also gets keyboard focus, otherwise the focus stays where it is (e.g. in the script editor).
void reveal(bool take_focus);

// What the user asked the host to do while the window was drawn.
struct Action {
    str_t insert_code = {};  // "Insert" was pressed on an example. Zero terminated, ends with a newline, valid until the next draw_window().
};

// Draws the viewer as a window named "Script Reference".
Action draw_window(bool* p_open);

}  // namespace script_reference
