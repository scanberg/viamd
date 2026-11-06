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
//
// Two kinds of quantity. The base ones (length, time, mass, temperature, angle) also build the
// units nothing names: Å²/ps, count/nm³. The composite ones (energy, force, ...) are recognised
// whole by their dimensions, because their usual units are not products of the base preferences:
// a force in kJ/(mol·nm) is not u·Å/ps². Velocity and diffusion may instead be Derived, built from
// length and time like the units nothing names; velocity is by default, diffusion is in 10⁻⁵ cm²/s.

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
    Quantity_Velocity,
    Quantity_Force,
    Quantity_Density,
    Quantity_Diffusion,
    Quantity_Dipole,
    Quantity_Count,
};

// Binds the preferences to the .ini. Call before app_settings::initialize().
void register_settings();

// The preferred unit; none for a quantity which is Derived
md_unit_t get(Quantity quantity);

// Bumps whenever a preference changes, including when the .ini is read. Anything holding values it
// has already converted compares this against the generation it converted at, and redoes the work
// when the two differ. Never zero, so a zero initialised generation always reads as stale.
uint64_t version();

// The unit 'src' is displayed in. Returns 'src' itself when no preference covers it.
md_unit_t target(md_unit_t src);

// Factor taking a value expressed in 'src' into its display unit, which is written to 'dst'.
// Returns 1 and leaves 'dst' equal to 'src' when there is nothing to convert.
double factor(md_unit_t* dst, md_unit_t src);

// Prints a unit the way it is shown: a preference by its own spelling (kJ/(mol·nm), not TN/mol),
// a unit built from the base preferences from theirs (nm/ns, not m/s), anything else as md_unit
// spells it, with small powers as superscripts.
size_t print(char* buf, size_t cap, md_unit_t unit);

// factor(), and the display unit printed into 'buf' by print()
double factor_print(char* buf, size_t cap, md_unit_t src);

// The combos, for the Settings menu.
void draw_settings_menu_items();

// ### For the settings and the tests ###

const char* quantity_label(Quantity quantity);
size_t      num_options(Quantity quantity);
const char* option_label(Quantity quantity, size_t option);
size_t      current_option(Quantity quantity);   // num_options() for a unit only a hand edited .ini names
void        select(Quantity quantity, size_t option);

// The text the .ini holds for each quantity, and the call that takes it after the .ini is read
const char* ini_key(Quantity quantity);
char* ini_buffer(Quantity quantity, size_t* cap);
void  adopt_ini();

// Every preference back to its default
void reset();

}  // namespace display_units
