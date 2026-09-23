#include "neutron_core.h"

namespace neutron {

namespace {

// Sears, Neutron News 3 (1992) 26, natural isotopic abundance. b_coh (fm), sigma_inc (barn).
struct Entry {
    unsigned char z;
    double b_coh;
    double sigma_inc;
};

const Entry table[] = {
    { 1, -3.7390, 80.26},   // H
    { 2,  3.26,    0.0},    // He
    { 3, -1.90,    0.92},   // Li
    { 4,  7.79,    0.0018}, // Be
    { 5,  5.30,    1.70},   // B (absorbing)
    { 6,  6.6460,  0.001},  // C
    { 7,  9.36,    0.50},   // N
    { 8,  5.803,   0.0008}, // O
    { 9,  5.654,   0.0008}, // F
    {10,  4.566,   0.008},  // Ne
    {11,  3.63,    1.62},   // Na
    {12,  5.375,   0.08},   // Mg
    {13,  3.449,   0.0082}, // Al
    {14,  4.1491,  0.004},  // Si
    {15,  5.13,    0.005},  // P
    {16,  2.847,   0.007},  // S
    {17,  9.5770,  5.3},    // Cl
    {18,  1.909,   0.225},  // Ar
    {19,  3.67,    0.27},   // K
    {20,  4.70,    0.05},   // Ca
    {21, 12.29,    4.5},    // Sc
    {22, -3.438,   2.87},   // Ti
    {23, -0.3824,  5.08},   // V
    {24,  3.635,   1.83},   // Cr
    {25, -3.73,    0.40},   // Mn
    {26,  9.45,    0.40},   // Fe
    {27,  2.49,    4.8},    // Co
    {28, 10.3,     5.2},    // Ni
    {29,  7.718,   0.55},   // Cu
    {30,  5.680,   0.077},  // Zn
    {31,  7.288,   0.16},   // Ga
    {32,  8.185,   0.18},   // Ge
    {33,  6.58,    0.06},   // As
    {34,  7.970,   0.32},   // Se
    {35,  6.795,   0.10},   // Br
    {36,  7.81,    0.01},   // Kr
    {37,  7.09,    0.5},    // Rb
    {38,  7.02,    0.06},   // Sr
    {40,  7.16,    0.02},   // Zr
    {42,  6.715,   0.04},   // Mo
    {47,  5.922,   0.58},   // Ag
    {48,  4.87,    3.46},   // Cd (absorbing)
    {50,  6.225,   0.022},  // Sn
    {53,  5.28,    0.31},   // I
    {55,  5.42,    0.21},   // Cs
    {56,  5.07,    0.15},   // Ba
    {78,  9.60,    0.13},   // Pt
    {79,  7.63,    0.43},   // Au
    {82,  9.405,   0.003},  // Pb
};

inline double clamp01(double x) { return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x); }

}  // namespace

bool element(int z, Element* out) {
    for (const Entry& e : table) {
        if (e.z == z) {
            out->b_coh = e.b_coh;
            out->sigma_inc = e.sigma_inc;
            return true;
        }
    }
    out->b_coh = 0.0;
    out->sigma_inc = 0.0;
    return false;
}

double water_sld(double d2o_fraction) {
    const double x = clamp01(d2o_fraction);
    return (1.0 - x) * SLD_H2O + x * SLD_D2O;
}

double hydrogen_b(double x) {
    x = clamp01(x);
    return (1.0 - x) * B_H + x * B_D;
}

// The incoherent cross section of an H/D mixture on one site is not the linear average in general (the spin
// incoherence of each isotope adds, but so does the isotope incoherence (b_H - b_D)^2 x (1 - x) 4 pi). The isotope
// term is included, since for a partially deuterated site it is of the same order as sigma_inc(D).
double hydrogen_sigma_inc(double x) {
    x = clamp01(x);
    const double db_fm = B_H - B_D;
    const double isotope = 4.0 * 3.14159265358979323846 * x * (1.0 - x) * db_fm * db_fm * 0.01;   // fm^2 -> barn
    return (1.0 - x) * SIGMA_INC_H + x * SIGMA_INC_D + isotope;
}

bool atom(Scattering* out, int z, double mass, bool exchangeable, const HydrogenModel& h, double mass_density) {
    const double m = (z == 1) ? H_MASS : mass;
    out->volume = mass_density > 0.0 ? m / (mass_density * DA_PER_A3_PER_GCM3) : 0.0;
    if (z == 1) {
        double x;
        if (exchangeable && h.exchange) {
            x = h.exchange_d;
        } else if (mass > 1.5) {
            x = 1.0;    // Deuterium in the topology
        } else {
            x = h.deuteration;
        }
        out->b = hydrogen_b(x);
        out->sigma_inc = hydrogen_sigma_inc(x);
        return true;
    }
    Element e;
    const bool known = element(z, &e);
    out->b = e.b_coh;
    out->sigma_inc = e.sigma_inc;
    return known;
}

Scattering particle(const Particle& p, const HydrogenModel& h) {
    const double x_ex = h.exchange ? h.exchange_d : h.deuteration;
    const double dbd = B_D - B_H;
    Scattering s;
    s.b = p.b_protiated + (p.n_h * clamp01(h.deuteration) + p.n_ex * clamp01(x_ex)) * dbd;
    s.sigma_inc = p.n_h * hydrogen_sigma_inc(h.deuteration) + p.n_ex * hydrogen_sigma_inc(x_ex);
    s.volume = p.volume;
    return s;
}

double absorption_sld(double number_density, double sigma_abs_thermal) {
    return number_density * sigma_abs_thermal * BARN_TO_A2 / (2.0 * 1.798) * 1.0e6;
}

}  // namespace neutron
