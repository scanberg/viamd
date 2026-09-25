#pragma once

// Neutron scattering lengths and H/D contrast for the scattering component.
//
// For neutrons the scattering weight of a particle is its excess coherent scattering length relative to the
// background medium,
//     w_j = b_j - SLD_bg * V_j
// where b_j is the coherent scattering length and V_j the volume the particle displaces. Unlike the X-ray path, where
// a uniform contrast factor is applied to the whole intensity, this is exact for mixtures, partial deuteration and
// (partially) contrast matched components, which is what contrast variation experiments rely on.
//
// Hydrogen is the reason neutrons are used: b_H = -3.74 fm, b_D = +6.67 fm. Hydrogen sites are either
//   - non-exchangeable: D fraction = deuteration (synthetic labelling), or
//   - exchangeable (bound to N, O or S): D fraction = D fraction of the exchange reservoir, typically the H2O/D2O
//     solvent or vapour.
// A hydrogen whose topology mass is that of deuterium (> 1.5 Da) is taken as fully deuterated (unless exchangeable).
//
// Units: b in fm, cross sections in barn, volumes in Å^3, SLD in 10^-6 Å^-2 unless stated otherwise.
// 1 fm/Å^3 = 10 x 10^-6 Å^-2 = 1e-5 Å^-2.

namespace neutron {

constexpr double B_H = -3.7390;             // fm (Sears 1992)
constexpr double B_D =  6.671;              // fm
constexpr double SIGMA_INC_H = 80.26;       // barn
constexpr double SIGMA_INC_D = 2.05;        // barn
constexpr double SLD_H2O = -0.56;           // 10^-6 Å^-2 (0.997 g/cm^3)
constexpr double SLD_D2O =  6.36;           // 10^-6 Å^-2 (1.104 g/cm^3)
constexpr double FM_TO_ANGSTROM = 1.0e-5;   // scattering length fm -> Å
constexpr double SLD_E6_TO_FM_PER_A3 = 0.1; // 10^-6 Å^-2 -> fm/Å^3
constexpr double BARN_TO_A2 = 1.0e-8;
constexpr double DA_PER_A3_PER_GCM3 = 0.602214076;  // V [Å^3] = m [Da] / (rho [g/cm^3] * 0.6022)
constexpr double H_MASS = 1.008;            // Da, used for the volume of any hydrogen isotope

// Natural abundance bound coherent scattering length (fm) and incoherent cross section (barn) of an element.
// Returns false for elements not in the table (out is zeroed). Imaginary parts (B, Cd, Gd, ...) are ignored.
struct Element {
    double b_coh;
    double sigma_inc;
};
bool element(int z, Element* out);

// Linear H2O/D2O mixture (volume fraction of D2O), 10^-6 Å^-2
double water_sld(double d2o_fraction);

// Hydrogen site with D fraction x in [0, 1]
double hydrogen_b(double x);
double hydrogen_sigma_inc(double x);

struct HydrogenModel {
    double deuteration = 0.0;   // D fraction of non-exchangeable hydrogen
    double exchange_d = 0.0;    // D fraction on exchangeable sites
    bool   exchange = true;     // false: exchangeable sites are treated as non-exchangeable
};

struct Scattering {
    double b = 0.0;             // fm
    double sigma_inc = 0.0;     // barn
    double volume = 0.0;        // Å^3
};

// Atom with atomic number z and topology mass (Da). exchangeable: hydrogen bound to N, O or S.
// The volume is mass / mass_density with the mass of any hydrogen isotope taken as H_MASS.
// Returns false if the element is unknown (b = sigma_inc = 0, the volume is still set).
bool atom(Scattering* out, int z, double mass, bool exchangeable, const HydrogenModel& h, double mass_density);

// Coarse grained particle: b_protiated is the coherent scattering length with every hydrogen as 1H, n_h the number of
// non-exchangeable and n_ex the number of exchangeable hydrogens (both included in b_protiated), volume in Å^3.
// sigma_inc only accounts for the hydrogens (the other common elements are negligible in comparison).
struct Particle {
    double b_protiated = 0.0;
    double n_h = 0.0;
    double n_ex = 0.0;
    double volume = 0.0;
};
Scattering particle(const Particle& p, const HydrogenModel& h);

// Excess scattering length (fm) relative to a background of the given SLD (10^-6 Å^-2)
inline double excess(const Scattering& s, double sld_bg) { return s.b - sld_bg * SLD_E6_TO_FM_PER_A3 * s.volume; }

// SLD (10^-6 Å^-2) of a material with total scattering length b (fm) in volume v (Å^3)
inline double sld(double b, double v) { return v > 0.0 ? b / (v * SLD_E6_TO_FM_PER_A3) : 0.0; }

// Absorption SLD (imaginary part, 10^-6 Å^-2). With sigma_a proportional to lambda this is wavelength independent:
// Im SLD = N sigma_a(1.798 Å) / (2 * 1.798 Å), N in 1/Å^3, sigma_a in barn.
double absorption_sld(double number_density, double sigma_abs_thermal);

}  // namespace neutron
