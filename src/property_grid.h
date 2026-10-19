#pragma once

// PROPERTY GRID
// How the settings of a thing are laid out (the properties of a representation): rows of a label in a column on the
// left and its control filling the rest, gathered in sections whose headers fold.
//
// The label column is one for a whole grid, however many panels draw into it (every representation and group of the
// Representations window): it ends where the furthest reaching label any of them drew the frame before ended, in the
// window, so the controls line up across sections, panels and depths of indentation, and no label is cut. It is
// kept in the window's state storage under the grid's name.
//
//   prop::begin("representation");                  // where the grid starts, this panel of it
//   prop::row("type");                               // the label; the control after it fills the rest of the row
//   ImGui::Combo("##type", ...);
//   if (prop::section("Color")) {                    // a folding header, rows indented under it
//       prop::row("coloring");
//       ImGui::Combo("##coloring", ...);
//       prop::hint("%zu atoms", n);                  // dimmed text in the control column
//       prop::section_end();
//   }
//   prop::end();
//
// Controls are ordinary ImGui calls with hidden labels ("##name"); a row sets the width of the next item to what is
// left of the row (SetNextItemWidth(-FLT_MIN)). Several controls on a row share it with SameLine and their own
// widths. A row with an empty label continues the one above it (the entries of a list).

#include <imgui.h>

namespace prop {

void begin(const char* grid_name);
void end();

// A folding header and the rows under it, indented. Open: returns true and section_end() is to follow; closed:
// returns false and there is nothing to end. The state is ImGui's, by the name, within the current ID scope.
bool section(const char* name, bool default_open = true);
void section_end();

// A section about something that can be turned off (an effect): the header's name is a checkbox for 'enabled', and
// the rows are there while it is on - returns true, and section_end() is to follow, when it is open and on.
// 'toggled' (optional) says whether the checkbox was clicked this frame; 'tooltip' is the checkbox's.
bool section_toggle(const char* name, bool* enabled, bool* toggled = nullptr, bool default_open = true, const char* tooltip = nullptr);

// A section about several things that are each turned on and off in its rows (effects that combine): the header's name
// is a checkbox for all of them, checked when all are on, mixed when some are and clear when none is. A click turns
// all of them off when any is on, else all on. The rows are there while it is open - returns true, and section_end()
// is to follow, when it is open. 'toggled' (optional) says whether the checkbox was clicked this frame.
bool section_toggle_all(const char* name, bool* const* enabled, int count, bool* toggled = nullptr, bool default_open = true, const char* tooltip = nullptr);

// A row: the label, and the control column for what follows. The tooltip, if any, is the label's.
void row(const char* label, const char* tooltip = nullptr);

// The width to give a control so that 'text' (dimmed, after it on the same row: a timing, a unit) fits beside it
float control_width_before(const char* text);

// Dimmed text in the control column, wrapped: what a setting is measured against, why there is nothing to set
void hint(const char* fmt, ...) IM_FMTARGS(1);

// A colour as a swatch across the control column; a click opens the picker. With alpha, rgba is 4 floats.
bool color(const char* id, float* rgba, bool alpha = false);

// The x (window local) where the controls start, and the width of the control column from there
float control_x();
float control_width();

}  // namespace prop
