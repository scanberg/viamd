#include "utest.h"

#include <display_units.h>
#include <script_engine.h>

#include <string.h>
#include <string>

using namespace display_units;

static md_unit_t U(const char* s) {
    md_unit_t u;
    if (!md_unit_parse(&u, str_from_cstr(s))) return md_unit_none();
    return u;
}

// The display unit of a value in src, printed, and the factor taking the value into it
static std::string shown(const char* src, double* scl = nullptr) {
    char buf[64];
    const double f = factor_print(buf, sizeof(buf), U(src));
    if (scl) *scl = f;
    return buf;
}

static size_t option_of(Quantity q, const char* label) {
    for (size_t i = 0; i < num_options(q); ++i) {
        if (strcmp(option_label(q, i), label) == 0) return i;
    }
    return num_options(q);
}

static void choose(Quantity q, const char* label) {
    const size_t i = option_of(q, label);
    MD_ASSERT(i < num_options(q));
    select(q, i);
}

UTEST(display_units, every_option_parses) {
    for (int q = 0; q < Quantity_Count; ++q) {
        for (size_t i = 0; i < num_options((Quantity)q); ++i) {
            select((Quantity)q, i);
            const bool derived = strncmp(option_label((Quantity)q, i), "Derived", 7) == 0;
            EXPECT_TRUE(derived == md_unit_is_none(get((Quantity)q)));
        }
    }
    reset();
}

// Built from length and time, and printed from their labels: nm/ns, never m/s
UTEST(display_units, velocity_derived) {
    reset();
    double f = 0;
    { const std::string got_ = shown("nm/ps", &f); EXPECT_STREQ("Å/ps", got_.c_str()); }
    EXPECT_NEAR(10.0, f, 1e-9);
    choose(Quantity_Length, "nm");
    choose(Quantity_Time, "ns");
    { const std::string got_ = shown("Å/ps", &f); EXPECT_STREQ("nm/ns", got_.c_str()); }
    EXPECT_NEAR(100.0, f, 1e-9);
    { const std::string got_ = shown("Å^3/ps", &f); EXPECT_STREQ("nm³/ns", got_.c_str()); }   // no quantity of its own
    { const std::string got_ = shown("count/Å^3"); EXPECT_STREQ("count/nm³", got_.c_str()); }
    reset();
}

UTEST(display_units, velocity_chosen) {
    reset();
    choose(Quantity_Velocity, "m/s");
    double f = 0;
    { const std::string got_ = shown("nm/ps", &f); EXPECT_STREQ("m/s", got_.c_str()); }
    EXPECT_NEAR(1000.0, f, 1e-9);
    // Not an angular velocity
    { const std::string got_ = shown("deg/ps"); EXPECT_STREQ("°/ps", got_.c_str()); }
    reset();
}

UTEST(display_units, force) {
    reset();
    double f = 0;
    { const std::string got_ = shown("kJ/(mol*nm)", &f); EXPECT_STREQ("kJ/(mol·nm)", got_.c_str()); }
    EXPECT_NEAR(1.0, f, 1e-12);
    { const std::string got_ = shown("kcal/(mol*Å)", &f); EXPECT_STREQ("kJ/(mol·nm)", got_.c_str()); }
    EXPECT_NEAR(41.84, f, 1e-9);
    // A force on one particle, a mole fewer
    { const std::string got_ = shown("eV/Å", &f); EXPECT_STREQ("kJ/nm", got_.c_str()); }
    choose(Quantity_Force, "pN");
    { const std::string got_ = shown("nN", &f); EXPECT_STREQ("pN", got_.c_str()); }
    EXPECT_NEAR(1000.0, f, 1e-9);
    reset();
}

UTEST(display_units, density_diffusion_dipole) {
    reset();
    double f = 0;
    { const std::string got_ = shown("u/Å^3", &f); EXPECT_STREQ("g/cm³", got_.c_str()); }
    EXPECT_NEAR(1.66053906892, f, 1e-6);
    // Diffusion as gmx msd reports it, unless Derived
    { const std::string got_ = shown("nm^2/ps", &f); EXPECT_STREQ("10⁻⁵ cm²/s", got_.c_str()); }
    EXPECT_NEAR(1000.0, f, 1e-6);
    choose(Quantity_Diffusion, "Derived (length²/time)");
    { const std::string got_ = shown("nm^2/ps", &f); EXPECT_STREQ("Å²/ps", got_.c_str()); }
    EXPECT_NEAR(100.0, f, 1e-9);
    { const std::string got_ = shown("e*Å", &f); EXPECT_STREQ("D", got_.c_str()); }
    EXPECT_NEAR(4.80320, f, 1e-4);
    reset();
}

UTEST(display_units, energy_per_particle) {
    reset();
    { const std::string got_ = shown("kcal/mol"); EXPECT_STREQ("kJ/mol", got_.c_str()); }
    { const std::string got_ = shown("eV"); EXPECT_STREQ("kJ", got_.c_str()); }
    choose(Quantity_Energy, "eV");
    { const std::string got_ = shown("kJ/mol"); EXPECT_STREQ("eV/mol", got_.c_str()); }
    reset();
}

UTEST(display_units, ini) {
    reset();
    size_t cap = 0;
    // A hand edited unit of the right dimension stands; Derived only where it is offered
    snprintf(ini_buffer(Quantity_Length, &cap), 32, "Mm");
    snprintf(ini_buffer(Quantity_Force, &cap), 32, "derived");
    snprintf(ini_buffer(Quantity_Velocity, &cap), 32, "derived");
    snprintf(ini_buffer(Quantity_Density, &cap), 32, "ps");
    adopt_ini();
    EXPECT_TRUE(md_unit_base_equal(get(Quantity_Length), md_unit_meter()));
    EXPECT_EQ(num_options(Quantity_Length), current_option(Quantity_Length));
    EXPECT_STREQ("kJ/(mol·nm)", option_label(Quantity_Force, current_option(Quantity_Force)));
    EXPECT_TRUE(md_unit_is_none(get(Quantity_Velocity)));
    EXPECT_STREQ("g/cm³", option_label(Quantity_Density, current_option(Quantity_Density)));
    reset();
}

// A hover's text: its unit and its value in the display units, the rest as it came
UTEST(display_units, mark_text) {
    reset();
    script::Mark m = {};
    m.text = STR_LIT("int km/s -> 10");
    m.unit = U("nm/ps");
    m.unit_at = 4;
    m.unit_len = 4;
    m.type_len = 8;
    m.has_number = true;
    m.number = 10.0;
    char buf[128];
    str_t t = script::mark_text(buf, sizeof(buf), m);
    { const std::string got_(t.ptr, t.len); EXPECT_STREQ("int Å/ps -> 100", got_.c_str()); }

    // Known per frame, no value
    m.text = STR_LIT("float[3] km/s, frame");
    m.unit_at = 9;
    m.type_len = (uint32_t)m.text.len;
    m.has_number = false;
    t = script::mark_text(buf, sizeof(buf), m);
    { const std::string got_(t.ptr, t.len); EXPECT_STREQ("float[3] Å/ps, frame", got_.c_str()); }

    // No unit: as it came
    script::Mark n = {};
    n.text = STR_LIT("int -> 3");
    t = script::mark_text(buf, sizeof(buf), n);
    { const std::string got_(t.ptr, t.len); EXPECT_STREQ("int -> 3", got_.c_str()); }
}
