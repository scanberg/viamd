#include <display_units.h>
#include <app_settings.h>

#include <core/md_common.h>
#include <core/md_log.h>

#include <imgui.h>

#include <string.h>

namespace display_units {

// The units offered per quantity. They are parsed with md_unit_parse, so these strings are also
// exactly what ends up in the .ini, and a hand written one only has to match a dimension rather
// than appear in this list.
static const char* const length_options[]      = { "Å", "nm", "pm", "bohr", "m" };
static const char* const angle_options[]       = { "°", "rad" };   // "deg" also parses, from a hand edited .ini
static const char* const time_options[]        = { "fs", "ps", "ns", "μs", "s" };
static const char* const energy_options[]      = { "kJ/mol", "kcal/mol", "eV", "Ha", "kJ", "J" };
static const char* const mass_options[]        = { "u", "kg" };
static const char* const charge_options[]      = { "e", "C" };
// md_unit_t is a scale without an offset, so the affine scales (°C, °F) have no representation.
static const char* const temperature_options[] = { "K" };
static const char* const pressure_options[]    = { "bar", "Pa", "kPa", "MPa" };

struct QuantityDesc {
    const char*        key;      // .ini key
    const char*        label;    // shown in the Settings menu
    const char* const* options;
    size_t             num_options;
    size_t             default_option;
};

#define DESC(key, label, opts, def) { key, label, opts, ARRAY_SIZE(opts), def }
static const QuantityDesc s_desc[Quantity_Count] = {
    DESC("unit_length",      "Length",      length_options,      0),  // Å
    DESC("unit_angle",       "Angle",       angle_options,       0),  // deg
    DESC("unit_time",        "Time",        time_options,        1),  // ps
    DESC("unit_energy",      "Energy",      energy_options,      0),  // kJ/mol
    DESC("unit_mass",        "Mass",        mass_options,        0),  // u
    DESC("unit_charge",      "Charge",      charge_options,      0),  // e
    DESC("unit_temperature", "Temperature", temperature_options, 0),  // K
    DESC("unit_pressure",    "Pressure",    pressure_options,    0),  // bar
};
#undef DESC

// The .ini stores the unit as text and that is what app_settings binds; the parsed form is derived
// from it. select() and adopt_str() are the only two places which write both.
static char      s_str[Quantity_Count][32];
static md_unit_t s_unit[Quantity_Count];

static uint64_t s_version = 1;

static md_unit_t parse_option(Quantity quantity, size_t option) {
    md_unit_t unit;
    const char* str = s_desc[quantity].options[option];
    if (!md_unit_parse(&unit, str_from_cstr(str))) {
        // An option this file spells wrong, which no .ini can cause and no user can fix.
        MD_LOG_ERROR("display_units: '%s' is not a parsable unit", str);
        return md_unit_none();
    }
    return unit;
}

static void select(Quantity quantity, size_t option) {
    ASSERT(option < s_desc[quantity].num_options);
    str_copy_to_char_buf(s_str[quantity], sizeof(s_str[quantity]), str_from_cstr(s_desc[quantity].options[option]));
    s_unit[quantity] = parse_option(quantity, option);
    s_version += 1;
}

// A unit is accepted for a quantity if it shares dimensions with one of the offered options, so a
// hand edited .ini can say 'Mm' or 'μs' without those having to be listed.
static bool dimension_is_valid(Quantity quantity, md_unit_t unit) {
    for (size_t i = 0; i < s_desc[quantity].num_options; ++i) {
        if (md_unit_base_equal(unit, parse_option(quantity, i))) {
            return true;
        }
    }
    return false;
}

// Takes whatever the .ini left in s_str and turns it into the parsed preference.
static void adopt_str(void*) {
    for (int i = 0; i < Quantity_Count; ++i) {
        const Quantity quantity = (Quantity)i;

        md_unit_t unit;
        if (!md_unit_parse(&unit, str_from_cstr(s_str[i])) || md_unit_is_none(unit)) {
            MD_LOG_ERROR("display_units: '%s' is not a unit, showing %s in %s instead", s_str[i], s_desc[i].label,
                         s_desc[i].options[s_desc[i].default_option]);
            select(quantity, s_desc[i].default_option);
            continue;
        }
        if (!dimension_is_valid(quantity, unit)) {
            MD_LOG_ERROR("display_units: '%s' is not a %s, showing %s instead", s_str[i], s_desc[i].label,
                         s_desc[i].options[s_desc[i].default_option]);
            select(quantity, s_desc[i].default_option);
            continue;
        }
        s_unit[i] = unit;
        s_version += 1;
    }
}

void register_settings() {
    for (int i = 0; i < Quantity_Count; ++i) {
        // The default has to stand before the .ini is read, since a file which says nothing about a
        // quantity leaves its string untouched.
        select((Quantity)i, s_desc[i].default_option);
        app_settings::bind(str_from_cstr(s_desc[i].key), s_str[i], sizeof(s_str[i]));
    }
    app_settings::on_apply(adopt_str, nullptr);
}

md_unit_t get(Quantity quantity) {
    ASSERT(quantity < Quantity_Count);
    return s_unit[quantity];
}

void set(Quantity quantity, md_unit_t unit) {
    ASSERT(quantity < Quantity_Count);
    s_unit[quantity] = unit;
    md_unit_print(s_str[quantity], sizeof(s_str[quantity]), unit);
    s_version += 1;
    app_settings::mark_dirty();
}

uint64_t version() {
    return s_version;
}

// True if 'src' is the preferred unit's quantity, allowing the two to disagree about the per-mole
// factor: a kJ/mol preference should still serve a plain joule, and it is the preference which is
// then adjusted rather than the value left unconverted.
static bool match_composite(md_unit_t* dst, md_unit_t src, md_unit_t pref) {
    if (md_unit_is_none(pref)) {
        return false;
    }

    md_unit_base_t a = src.base;
    md_unit_base_t b = pref.base;
    const int mole_diff = a.dim.mole - b.dim.mole;
    a.dim.mole = 0;
    b.dim.mole = 0;
    if (a.raw_bits != b.raw_bits) {
        return false;
    }

    *dst = mole_diff != 0 ? md_unit_mul(pref, md_unit_pow(md_unit_mole(), mole_diff)) : pref;
    return true;
}

static md_unit_t pref_or(Quantity quantity, md_unit_t fallback) {
    const md_unit_t unit = s_unit[quantity];
    return md_unit_is_none(unit) ? fallback : unit;
}

md_unit_t target(md_unit_t src) {
    if (md_unit_is_none(src)) {
        return src;
    }

    // Energy, charge and pressure are matched as a whole. Their preferred unit can carry a
    // dimension which does not survive being taken apart -- kJ/mol is an energy with a mole in it,
    // and the elementary charge is an ampere-second nobody wants spelled out that way -- so they
    // are recognised by their dimension signature rather than rebuilt from the base preferences.
    md_unit_t dst;
    if (match_composite(&dst, src, s_unit[Quantity_Energy]))   return dst;
    if (match_composite(&dst, src, s_unit[Quantity_Charge]))   return dst;
    if (match_composite(&dst, src, s_unit[Quantity_Pressure])) return dst;

    // Everything else is rebuilt one base dimension at a time, which is what carries Å² to nm²,
    // Å/ps to nm/ns and count/Å³ to count/nm³ without any of those being spelled out anywhere.
    // A unit reaching into current or mole without having matched above -- a dipole moment, say --
    // has no sensible decomposition and is left exactly as it came.
    if (src.base.dim.current != 0 || src.base.dim.mole != 0) {
        return src;
    }

    dst = md_unit_none();
    dst = md_unit_mul(dst, md_unit_pow(pref_or(Quantity_Length,      md_unit_meter()),    src.base.dim.length));
    dst = md_unit_mul(dst, md_unit_pow(pref_or(Quantity_Mass,        md_unit_kilogram()), src.base.dim.mass));
    dst = md_unit_mul(dst, md_unit_pow(pref_or(Quantity_Time,        md_unit_second()),   src.base.dim.time));
    dst = md_unit_mul(dst, md_unit_pow(pref_or(Quantity_Temperature, md_unit_kelvin()),   src.base.dim.temp));
    dst = md_unit_mul(dst, md_unit_pow(pref_or(Quantity_Angle,       md_unit_radian()),   src.base.dim.angle));
    // No preference for a bare count, but it has to be carried so the dimensions still line up.
    dst = md_unit_mul(dst, md_unit_pow(md_unit_count(), src.base.dim.count));
    return dst;
}

double factor(md_unit_t* dst, md_unit_t src) {
    const md_unit_t unit = target(src);

    double scl = 1.0;
    if (!md_unit_conversion_factor(&scl, src, unit)) {
        // target() only ever returns a unit of the same dimensions, so this is unreachable short of
        // a bug in here. Showing the value in the unit it arrived in is the harmless way to fail.
        ASSERT(false);
        if (dst) *dst = src;
        return 1.0;
    }

    if (dst) *dst = unit;
    return scl;
}

double factor_print(char* buf, size_t cap, md_unit_t src) {
    md_unit_t unit;
    const double scl = factor(&unit, src);
    md_unit_print(buf, cap, unit);
    return scl;
}

void draw_settings_menu_items() {
    for (int i = 0; i < Quantity_Count; ++i) {
        const QuantityDesc& desc = s_desc[i];
        if (ImGui::BeginCombo(desc.label, s_str[i])) {
            for (size_t j = 0; j < desc.num_options; ++j) {
                if (ImGui::Selectable(desc.options[j], strcmp(s_str[i], desc.options[j]) == 0)) {
                    select((Quantity)i, j);
                    app_settings::mark_dirty();
                }
            }
            ImGui::EndCombo();
        }
    }
}

}  // namespace display_units
