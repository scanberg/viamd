#include <display_units.h>

#include <core/md_common.h>
#include <core/md_log.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

namespace display_units {

// The units offered per quantity. The first string is parsed with md_unit_parse and is what the .ini holds, so a hand
// written one only has to match a dimension rather than appear in this list; the second is how it is shown, which
// md_unit cannot always spell (it prints kJ/(mol*nm) as TN/mol, by its scale).
struct Option {
    const char* ini;
    const char* label;
};

// Velocity and diffusion may be built from length and time like the units nothing names, which is what they did
// before they could be chosen. Velocity is by default, so it reads in the length shown everywhere else; a diffusion
// coefficient has a unit of its own by convention, 10⁻⁵ cm²/s, in gmx msd and in the literature alike, and is shown in it.
static const char DERIVED[] = "derived";

static const Option length_options[]      = { {"Å", "Å"}, {"nm", "nm"}, {"pm", "pm"}, {"bohr", "bohr"}, {"m", "m"} };
static const Option angle_options[]       = { {"°", "°"}, {"rad", "rad"} };   // "deg" also parses, from a hand edited .ini
static const Option time_options[]        = { {"fs", "fs"}, {"ps", "ps"}, {"ns", "ns"}, {"μs", "μs"}, {"s", "s"} };
static const Option energy_options[]      = { {"kJ/mol", "kJ/mol"}, {"kcal/mol", "kcal/mol"}, {"eV", "eV"}, {"Ha", "Ha"}, {"kJ", "kJ"}, {"J", "J"} };
static const Option mass_options[]        = { {"u", "u"}, {"kg", "kg"} };
static const Option charge_options[]      = { {"e", "e"}, {"C", "C"} };
// md_unit_t is a scale without an offset, so the affine scales (°C, °F) have no representation.
static const Option temperature_options[] = { {"K", "K"} };
static const Option pressure_options[]    = { {"bar", "bar"}, {"Pa", "Pa"}, {"kPa", "kPa"}, {"MPa", "MPa"} };
static const Option velocity_options[]    = { {DERIVED, "Derived (length/time)"}, {"Å/ps", "Å/ps"}, {"nm/ps", "nm/ps"}, {"Å/fs", "Å/fs"}, {"m/s", "m/s"}, {"km/s", "km/s"} };
// GROMACS writes forces in kJ/(mol·nm), so that is where a .trr's forces are shown as they came
static const Option force_options[]       = { {"kJ/(mol*nm)", "kJ/(mol·nm)"}, {"kJ/(mol*Å)", "kJ/(mol·Å)"}, {"kcal/(mol*Å)", "kcal/(mol·Å)"}, {"eV/Å", "eV/Å"}, {"pN", "pN"} };
// md_unit has no gram, so g/cm³ is spelled by its scale
static const Option density_options[]     = { {"1e3*kg/m^3", "g/cm³"}, {"kg/m^3", "kg/m³"}, {"u/Å^3", "u/Å³"} };
static const Option diffusion_options[]   = { {DERIVED, "Derived (length²/time)"}, {"1e-5*cm^2/s", "10⁻⁵ cm²/s"}, {"nm^2/ps", "nm²/ps"}, {"Å^2/ps", "Å²/ps"}, {"μm^2/s", "μm²/s"}, {"cm^2/s", "cm²/s"} };
static const Option dipole_options[]      = { {"D", "D"}, {"e*Å", "e·Å"}, {"e*nm", "e·nm"} };

struct QuantityDesc {
    const char*   key;      // .ini key
    const char*   label;    // shown in the Settings menu
    const Option* options;
    size_t        num_options;
    size_t        default_option;
};

#define DESC(key, label, opts, def) { key, label, opts, MD_COUNTOF(opts), def }
static const QuantityDesc s_desc[Quantity_Count] = {
    DESC("unit_length",      "Length",      length_options,      0),  // Å
    DESC("unit_angle",       "Angle",       angle_options,       0),  // deg
    DESC("unit_time",        "Time",        time_options,        1),  // ps
    DESC("unit_energy",      "Energy",      energy_options,      0),  // kJ/mol
    DESC("unit_mass",        "Mass",        mass_options,        0),  // u
    DESC("unit_charge",      "Charge",      charge_options,      0),  // e
    DESC("unit_temperature", "Temperature", temperature_options, 0),  // K
    DESC("unit_pressure",    "Pressure",    pressure_options,    0),  // bar
    DESC("unit_velocity",    "Velocity",    velocity_options,    0),  // derived
    DESC("unit_force",       "Force",       force_options,       0),  // kJ/(mol·nm)
    DESC("unit_density",     "Density",     density_options,     0),  // g/cm³
    DESC("unit_diffusion",   "Diffusion",   diffusion_options,   1),  // 10⁻⁵ cm²/s, as gmx msd reports it
    DESC("unit_dipole",      "Dipole",      dipole_options,      0),  // D
};
#undef DESC

// The .ini stores the unit as text and that is what app_settings binds; the parsed form is derived
// from it. select() and adopt_ini() are the only two places which write both.
static char      s_str[Quantity_Count][32];
static md_unit_t s_unit[Quantity_Count];

static uint64_t s_version = 1;

static bool is_derived(const char* str) { return strcmp(str, DERIVED) == 0; }

static md_unit_t parse_option(Quantity quantity, size_t option) {
    const char* str = s_desc[quantity].options[option].ini;
    if (is_derived(str)) return md_unit_none();
    md_unit_t unit;
    if (!md_unit_parse(&unit, md_str_from_cstr(str))) {
        // An option this file spells wrong, which no .ini can cause and no user can fix.
        MD_LOG_ERROR("display_units: '%s' is not a parsable unit", str);
        return md_unit_none();
    }
    return unit;
}

void select(Quantity quantity, size_t option) {
    MD_ASSERT(quantity < Quantity_Count && option < s_desc[quantity].num_options);
    md_str_copy_to_char_buf(s_str[quantity], sizeof(s_str[quantity]), md_str_from_cstr(s_desc[quantity].options[option].ini));
    s_unit[quantity] = parse_option(quantity, option);
    s_version += 1;
}

void reset() {
    for (int i = 0; i < Quantity_Count; ++i) select((Quantity)i, s_desc[i].default_option);
}

// A unit is accepted for a quantity if it shares dimensions with one of the offered options, so a
// hand edited .ini can say 'Mm' or 'μs' without those having to be listed.
static bool dimension_is_valid(Quantity quantity, md_unit_t unit) {
    for (size_t i = 0; i < s_desc[quantity].num_options; ++i) {
        const md_unit_t opt = parse_option(quantity, i);
        if (!md_unit_is_none(opt) && md_unit_base_equal(unit, opt)) {
            return true;
        }
    }
    return false;
}

static bool has_derived(Quantity quantity) {
    for (size_t i = 0; i < s_desc[quantity].num_options; ++i) {
        if (is_derived(s_desc[quantity].options[i].ini)) return true;
    }
    return false;
}

// Takes whatever the .ini left in s_str and turns it into the parsed preference.
void adopt_ini() {
    for (int i = 0; i < Quantity_Count; ++i) {
        const Quantity quantity = (Quantity)i;
        if (is_derived(s_str[i]) && has_derived(quantity)) {
            s_unit[i] = md_unit_none();
            s_version += 1;
            continue;
        }

        md_unit_t unit;
        if (!md_unit_parse(&unit, md_str_from_cstr(s_str[i])) || md_unit_is_none(unit)) {
            MD_LOG_ERROR("display_units: '%s' is not a unit, showing %s in %s instead", s_str[i], s_desc[i].label,
                         s_desc[i].options[s_desc[i].default_option].label);
            select(quantity, s_desc[i].default_option);
            continue;
        }
        if (!dimension_is_valid(quantity, unit)) {
            MD_LOG_ERROR("display_units: '%s' is not a %s, showing %s instead", s_str[i], s_desc[i].label,
                         s_desc[i].options[s_desc[i].default_option].label);
            select(quantity, s_desc[i].default_option);
            continue;
        }
        s_unit[i] = unit;
        s_version += 1;
    }
}

char* ini_buffer(Quantity quantity, size_t* cap) {
    MD_ASSERT(quantity < Quantity_Count);
    if (cap) *cap = sizeof(s_str[quantity]);
    return s_str[quantity];
}

const char* ini_key(Quantity quantity) { return s_desc[quantity].key; }
const char* quantity_label(Quantity quantity) { return s_desc[quantity].label; }
size_t num_options(Quantity quantity) { return s_desc[quantity].num_options; }
const char* option_label(Quantity quantity, size_t option) { return s_desc[quantity].options[option].label; }

size_t current_option(Quantity quantity) {
    for (size_t i = 0; i < s_desc[quantity].num_options; ++i) {
        if (strcmp(s_str[quantity], s_desc[quantity].options[i].ini) == 0) return i;
    }
    return s_desc[quantity].num_options;
}

// How a preference is shown: its option's label, or what a hand edited .ini wrote
static const char* current_label(Quantity quantity) {
    const size_t i = current_option(quantity);
    return i < s_desc[quantity].num_options ? s_desc[quantity].options[i].label : s_str[quantity];
}

md_unit_t get(Quantity quantity) {
    MD_ASSERT(quantity < Quantity_Count);
    return s_unit[quantity];
}

uint64_t version() {
    return s_version;
}

// The quantities recognised whole by their dimensions, before anything is built from the base ones
static const Quantity s_composites[] = {
    Quantity_Energy, Quantity_Charge, Quantity_Pressure, Quantity_Force, Quantity_Dipole, Quantity_Density,
    Quantity_Velocity, Quantity_Diffusion,
};

// True if 'src' is the preferred unit's quantity, allowing the two to disagree about the per-mole
// factor: a kJ/mol preference should still serve a plain joule, and it is the preference which is
// then adjusted rather than the value left unconverted. The difference goes to 'mole_diff'.
static bool match_composite(md_unit_t* dst, int* mole_diff, md_unit_t src, md_unit_t pref) {
    if (md_unit_is_none(pref)) {
        return false;
    }

    md_unit_base_t a = src.base;
    md_unit_base_t b = pref.base;
    const int diff = a.dim.mole - b.dim.mole;
    a.dim.mole = 0;
    b.dim.mole = 0;
    if (a.raw_bits != b.raw_bits) {
        return false;
    }

    *dst = diff != 0 ? md_unit_mul(pref, md_unit_pow(md_unit_mole(), diff)) : pref;
    if (mole_diff) *mole_diff = diff;
    return true;
}

// The base quantities a unit is built from, in the order its label names them
struct BaseDim {
    Quantity    quantity;   // Quantity_Count for count, which has no preference
    md_unit_t (*fallback)();
    const char* fallback_label;
};

static md_unit_t unit_count() { return md_unit_count(); }

static const BaseDim s_base[] = {
    {Quantity_Mass,        md_unit_kilogram, "kg"},
    {Quantity_Length,      md_unit_meter,    "m"},
    {Quantity_Time,        md_unit_second,   "s"},
    {Quantity_Temperature, md_unit_kelvin,   "K"},
    {Quantity_Angle,       md_unit_radian,   "rad"},
    {Quantity_Count,       unit_count,       "count"},
};

static int base_power(md_unit_t u, size_t i) {
    switch (s_base[i].quantity) {
    case Quantity_Mass:        return u.base.dim.mass;
    case Quantity_Length:      return u.base.dim.length;
    case Quantity_Time:        return u.base.dim.time;
    case Quantity_Temperature: return u.base.dim.temp;
    case Quantity_Angle:       return u.base.dim.angle;
    default:                   return u.base.dim.count;
    }
}

static md_unit_t base_unit(size_t i) {
    const Quantity q = s_base[i].quantity;
    if (q != Quantity_Count && !md_unit_is_none(s_unit[q])) return s_unit[q];
    return s_base[i].fallback();
}

static const char* base_label(size_t i) {
    const Quantity q = s_base[i].quantity;
    if (q != Quantity_Count && !md_unit_is_none(s_unit[q])) return current_label(q);
    return s_base[i].fallback_label;
}

static bool can_rebuild(md_unit_t src) {
    // A unit reaching into current or mole without having matched as a whole -- a dipole moment in
    // C·m with no dipole preference, say -- has no sensible decomposition.
    return src.base.dim.current == 0 && src.base.dim.mole == 0;
}

static md_unit_t rebuild(md_unit_t src) {
    md_unit_t dst = md_unit_none();
    for (size_t i = 0; i < MD_COUNTOF(s_base); ++i) {
        const int p = base_power(src, i);
        if (p) dst = md_unit_mul(dst, md_unit_pow(base_unit(i), p));
    }
    return dst;
}

md_unit_t target(md_unit_t src) {
    if (md_unit_is_none(src)) {
        return src;
    }

    // The composites are matched as a whole. Their preferred unit can carry a dimension which does
    // not survive being taken apart -- kJ/mol is an energy with a mole in it, and the elementary
    // charge is an ampere-second nobody wants spelled out that way -- so they are recognised by
    // their dimension signature rather than rebuilt from the base preferences.
    md_unit_t dst;
    for (size_t i = 0; i < MD_COUNTOF(s_composites); ++i) {
        if (match_composite(&dst, nullptr, src, s_unit[s_composites[i]])) return dst;
    }

    // Everything else is rebuilt one base dimension at a time, which is what carries Å² to nm²,
    // Å/ps to nm/ns and count/Å³ to count/nm³ without any of those being spelled out anywhere.
    if (!can_rebuild(src)) {
        return src;
    }
    return rebuild(src);
}

double factor(md_unit_t* dst, md_unit_t src) {
    const md_unit_t unit = target(src);

    double scl = 1.0;
    if (!md_unit_conversion_factor(&scl, src, unit)) {
        // target() only ever returns a unit of the same dimensions, so this is unreachable short of
        // a bug in here. Showing the value in the unit it arrived in is the harmless way to fail.
        MD_ASSERT(false);
        if (dst) *dst = src;
        return 1.0;
    }

    if (dst) *dst = unit;
    return scl;
}

// ### PRINTING ###

struct Out {
    char*  buf;
    size_t cap;
    size_t len;
};

static void put(Out* o, const char* s) {
    for (; *s; ++s) {
        if (o->len + 1 < o->cap) o->buf[o->len] = *s;
        o->len += 1;
    }
    if (o->cap) o->buf[MD_MIN(o->len, o->cap - 1)] = '\0';
}

// ² to ⁹ as UTF-8 bytes; a power past 9 is written ^n
static void put_power(Out* o, int p) {
    static const char* sup[] = {"", "", "\xC2\xB2", "\xC2\xB3", "\xE2\x81\xB4", "\xE2\x81\xB5", "\xE2\x81\xB6", "\xE2\x81\xB7", "\xE2\x81\xB8", "\xE2\x81\xB9"};
    if (p <= 1) return;
    if (p <= 9) { put(o, sup[p]); return; }
    char tmp[16];
    snprintf(tmp, sizeof(tmp), "^%d", p);
    put(o, tmp);
}

static bool same_unit(md_unit_t a, md_unit_t b) {
    if (a.base.raw_bits != b.base.raw_bits) return false;
    const double sa = md_unit_scale(a), sb = md_unit_scale(b);
    return fabs(sa - sb) <= 1e-9 * MD_MAX(fabs(sa), fabs(sb));
}

// The label of a unit built from the base preferences: those raised to a power above, those to a power below
static size_t print_rebuilt(char* buf, size_t cap, md_unit_t unit) {
    Out o = {buf, cap, 0};
    if (cap) buf[0] = '\0';
    int num = 0, den = 0;
    for (size_t i = 0; i < MD_COUNTOF(s_base); ++i) {
        const int p = base_power(unit, i);
        num += p > 0;
        den += p < 0;
    }
    bool first = true;
    for (size_t i = 0; i < MD_COUNTOF(s_base); ++i) {
        const int p = base_power(unit, i);
        if (p <= 0) continue;
        if (!first) put(&o, "\xC2\xB7");   // ·
        put(&o, base_label(i));
        put_power(&o, p);
        first = false;
    }
    if (den) {
        if (!num) put(&o, "1");
        put(&o, "/");
        if (den > 1) put(&o, "(");
        first = true;
        for (size_t i = 0; i < MD_COUNTOF(s_base); ++i) {
            const int p = base_power(unit, i);
            if (p >= 0) continue;
            if (!first) put(&o, "\xC2\xB7");
            put(&o, base_label(i));
            put_power(&o, -p);
            first = false;
        }
        if (den > 1) put(&o, ")");
    }
    return o.len;
}

// A composite's label with one mole fewer or more: kJ/mol for a plain joule is kJ, kJ/(mol·nm) is kJ/nm, eV for a
// molar energy is eV/mol. False when the label is not one of those shapes.
static bool print_mole_adjusted(char* buf, size_t cap, const char* label, int mole_diff) {
    char tmp[64];
    const size_t n = strlen(label);
    if (mole_diff == 1) {
        if (n > 4 && strcmp(label + n - 4, "/mol") == 0) {
            snprintf(buf, cap, "%.*s", (int)(n - 4), label);
            return true;
        }
        const char* m = strstr(label, "/(mol\xC2\xB7");
        if (m) {
            const char* rest = m + 7;   // past "/(mol·"
            const size_t rn = strlen(rest);
            if (rn > 0 && rest[rn - 1] == ')' && !strchr(rest, '\xC2')) {
                snprintf(tmp, sizeof(tmp), "%.*s/%.*s", (int)(m - label), label, (int)(rn - 1), rest);
                snprintf(buf, cap, "%s", tmp);
                return true;
            }
        }
        return false;
    }
    if (mole_diff == -1 && !strchr(label, '/')) {
        snprintf(buf, cap, "%s/mol", label);
        return true;
    }
    return false;
}

size_t print(char* buf, size_t cap, md_unit_t unit) {
    if (!buf || cap == 0) return 0;
    buf[0] = '\0';
    if (md_unit_is_none(unit)) return 0;

    // A preference, or a composite's preference a mole apart
    for (int i = 0; i < Quantity_Count; ++i) {
        const md_unit_t pref = s_unit[i];
        if (md_unit_is_none(pref)) continue;
        if (same_unit(unit, pref)) {
            snprintf(buf, cap, "%s", current_label((Quantity)i));
            return strlen(buf);
        }
    }
    for (size_t i = 0; i < MD_COUNTOF(s_composites); ++i) {
        const Quantity q = s_composites[i];
        md_unit_t m;
        int diff = 0;
        if (match_composite(&m, &diff, unit, s_unit[q]) && diff != 0 && same_unit(unit, m) &&
            print_mole_adjusted(buf, cap, current_label(q), diff)) {
            return strlen(buf);
        }
    }

    // Built from the base preferences
    if (can_rebuild(unit) && same_unit(unit, rebuild(unit))) {
        return print_rebuilt(buf, cap, unit);
    }

    // Anything else as md_unit spells it, small powers raised
    char raw[64];
    const size_t n = md_unit_print(raw, sizeof(raw), unit);
    Out o = {buf, cap, 0};
    for (size_t i = 0; i < n && i < sizeof(raw) - 1; ++i) {
        if (raw[i] == '^' && i + 1 < n && raw[i + 1] >= '2' && raw[i + 1] <= '9' && (i + 2 >= n || raw[i + 2] < '0' || raw[i + 2] > '9')) {
            put_power(&o, raw[i + 1] - '0');
            i += 1;
            continue;
        }
        const char c[2] = {raw[i], 0};
        put(&o, c);
    }
    return o.len;
}

double factor_print(char* buf, size_t cap, md_unit_t src) {
    md_unit_t unit;
    const double scl = factor(&unit, src);
    print(buf, cap, unit);
    return scl;
}

}  // namespace display_units
