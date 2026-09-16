#pragma once

#include <stddef.h>
#include <stdint.h>

#include <core/md_str.h>
#include <core/md_unit.h>

// The units quantities are shown in.
//
// One preferred unit per physical quantity, persisted in the ImGui .ini through app_settings.
// These are display defaults, not a property of the data: a value keeps whatever unit its source
// gave it and the conversion happens on the way to the screen, so a view with a context of its own
// -- an orbital energy in Hartree, say -- is free to ignore all of this and print its own unit.
//
// A display site converts by asking for the factor and the unit the value lands in:
//
//     char unit_buf[32];
//     const double scl = display_units::factor_print(unit_buf, sizeof(unit_buf), prop.unit);
//     // plot 'value * scl', label the axis with 'unit_buf'
//
// Anything the preferences do not cover comes back unchanged with a factor of one, so a site can
// route every value through this without first asking whether it is convertible.

namespace display_units {

enum Quantity {
    Quantity_Length,
    Quantity_Angle,
    Quantity_Time,
    Quantity_Energy,
    Quantity_Mass,
    Quantity_Charge,
    Quantity_Temperature,
    Quantity_Pressure,
    Quantity_Count,
};

// Binds the preferences to the .ini. Call before app_settings::initialize().
void register_settings();

md_unit_t get(Quantity quantity);
void      set(Quantity quantity, md_unit_t unit);

// Bumps whenever a preference changes, including when the .ini is read. Anything holding values it
// has already converted compares this against the generation it converted at, and redoes the work
// when the two differ. Never zero, so a zero initialised generation always reads as stale.
uint64_t version();

// The unit 'src' is displayed in. Returns 'src' itself when no preference covers it.
md_unit_t target(md_unit_t src);

// Factor taking a value expressed in 'src' into its display unit, which is written to 'dst'.
// Returns 1 and leaves 'dst' equal to 'src' when there is nothing to convert.
double factor(md_unit_t* dst, md_unit_t src);

// Same, but prints the display unit into 'buf' the way md_unit_print does.
double factor_print(char* buf, size_t cap, md_unit_t src);

// The combos, for the Settings menu.
void draw_settings_menu_items();

}  // namespace display_units
