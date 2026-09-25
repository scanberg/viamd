#include "utest.h"

#include "../neutron_core.h"

#include <math.h>

// Neutron contrast for the scattering component. The numbers checked here are the ones a SANS / GISANS user would
// compare against: the scattering length of a known molecule, the SLD of water and of the common substrates, and the
// match point of cellulose, which depends on hydrogen exchange.

namespace {

constexpr double H_DA = 1.008, C_DA = 12.011, O_DA = 15.999;

// Anhydroglucose unit C6H10O5 (cellulose), 3 exchangeable (hydroxyl) hydrogens, 7 non-exchangeable
neutron::Scattering anhydroglucose(const neutron::HydrogenModel& h, double rho) {
    neutron::Scattering sum;
    auto add = [&](int z, double m, bool ex) {
        neutron::Scattering s;
        neutron::atom(&s, z, m, ex, h, rho);
        sum.b += s.b; sum.sigma_inc += s.sigma_inc; sum.volume += s.volume;
    };
    for (int i = 0; i < 6; ++i) add(6, C_DA, false);
    for (int i = 0; i < 5; ++i) add(8, O_DA, false);
    for (int i = 0; i < 7; ++i) add(1, H_DA, false);
    for (int i = 0; i < 3; ++i) add(1, H_DA, true);
    return sum;
}

}  // namespace

UTEST(neutron, element_table) {
    neutron::Element e;
    ASSERT_TRUE(neutron::element(6, &e));
    EXPECT_NEAR(e.b_coh, 6.646, 1e-3);
    ASSERT_TRUE(neutron::element(1, &e));
    EXPECT_NEAR(e.b_coh, neutron::B_H, 1e-9);
    EXPECT_FALSE(neutron::element(0, &e));
    EXPECT_EQ(e.b_coh, 0.0);
    EXPECT_FALSE(neutron::element(118, &e));
}

UTEST(neutron, water_sld) {
    EXPECT_NEAR(neutron::water_sld(0.0), -0.56, 1e-9);
    EXPECT_NEAR(neutron::water_sld(1.0), 6.36, 1e-9);
    // Clamped outside [0, 1]
    EXPECT_NEAR(neutron::water_sld(2.0), 6.36, 1e-9);
    // Independent check from composition: D2O at 1.104 g/cm^3
    const double b = 2 * neutron::B_D + 5.803;
    const double n = 1.104 * neutron::DA_PER_A3_PER_GCM3 / 20.028;    // molecules / Å^3
    EXPECT_NEAR(b * n / neutron::SLD_E6_TO_FM_PER_A3, 6.36, 0.02);
}

UTEST(neutron, substrate_slds) {
    // Si: 2.329 g/cm^3, 28.086 Da
    neutron::Element si;
    neutron::element(14, &si);
    const double n_si = 2.329 * neutron::DA_PER_A3_PER_GCM3 / 28.086;
    EXPECT_NEAR(si.b_coh * n_si / neutron::SLD_E6_TO_FM_PER_A3, 2.07, 0.01);
    // Amorphous SiO2: 2.2 g/cm^3
    const double b_sio2 = si.b_coh + 2 * 5.803;
    const double n_sio2 = 2.2 * neutron::DA_PER_A3_PER_GCM3 / 60.084;
    EXPECT_NEAR(b_sio2 * n_sio2 / neutron::SLD_E6_TO_FM_PER_A3, 3.47, 0.01);
    // Absorption of Si is negligible (sigma_a = 0.171 b at 1.798 Å)
    EXPECT_LT(neutron::absorption_sld(n_si, 0.171), 1.0e-4);
}

UTEST(neutron, cellulose_contrast) {
    neutron::HydrogenModel h;
    h.exchange = true;
    const double rho = 1.5;

    // Fully protiated: b = 6 * 6.646 + 5 * 5.803 + 10 * -3.739 = 31.50 fm, SLD ~1.75
    h.exchange_d = 0.0;
    const neutron::Scattering s_h = anhydroglucose(h, rho);
    EXPECT_NEAR(s_h.b, 31.50, 0.01);
    const double sld_h = neutron::sld(s_h.b, s_h.volume);
    EXPECT_NEAR(sld_h, 1.75, 0.02);

    // In D2O the three hydroxyl H exchange: +3 * (b_D - b_H) = +31.23 fm, SLD ~3.5
    h.exchange_d = 1.0;
    const neutron::Scattering s_d = anhydroglucose(h, rho);
    EXPECT_NEAR(s_d.b - s_h.b, 3 * (neutron::B_D - neutron::B_H), 1e-9);
    EXPECT_NEAR(neutron::sld(s_d.b, s_d.volume), 3.49, 0.03);
    EXPECT_NEAR(s_d.volume, s_h.volume, 1e-9);   // Isotopes do not change the volume

    // Match point: the D2O fraction x where the exchanged cellulose SLD equals the water SLD. Both are linear in x.
    auto excess_at = [&](double x) {
        neutron::HydrogenModel hx = h;
        hx.exchange_d = x;
        return neutron::excess(anhydroglucose(hx, rho), neutron::water_sld(x));
    };
    const double e0 = excess_at(0.0), e1 = excess_at(1.0);
    const double x_match = e0 / (e0 - e1);
    EXPECT_GT(x_match, 0.35);
    EXPECT_LT(x_match, 0.45);
    EXPECT_NEAR(excess_at(x_match), 0.0, 1e-9);

    // Without exchange the hydroxyl H stay 1H in D2O
    h.exchange = false;
    h.exchange_d = 1.0;
    EXPECT_NEAR(anhydroglucose(h, rho).b, s_h.b, 1e-9);
}

UTEST(neutron, particle_matches_atoms) {
    // 22 anhydroglucose units as one coarse grained bead
    neutron::HydrogenModel h;
    h.deuteration = 0.25;
    h.exchange_d = 0.6;
    neutron::Scattering ref = anhydroglucose(h, 1.5);
    neutron::Particle p;
    neutron::HydrogenModel h0;
    h0.exchange_d = 0.0;
    p.b_protiated = 22 * anhydroglucose(h0, 1.5).b;
    p.n_h = 22 * 7;
    p.n_ex = 22 * 3;
    p.volume = 22 * ref.volume;
    const neutron::Scattering s = neutron::particle(p, h);
    EXPECT_NEAR(s.b, 22 * ref.b, 1e-6);
    // Hydrogens only in the particle model; C and O contribute ~0.005 b per unit
    EXPECT_NEAR(s.sigma_inc, 22 * ref.sigma_inc, 22 * 0.01);
}

UTEST(neutron, hydrogen_incoherent) {
    EXPECT_NEAR(neutron::hydrogen_sigma_inc(0.0), neutron::SIGMA_INC_H, 1e-9);
    EXPECT_NEAR(neutron::hydrogen_sigma_inc(1.0), neutron::SIGMA_INC_D, 1e-9);
    // Isotope disorder: 4 pi x (1 - x) (b_H - b_D)^2, 3.4 b at x = 0.5
    const double mid = 0.5 * (neutron::SIGMA_INC_H + neutron::SIGMA_INC_D);
    EXPECT_NEAR(neutron::hydrogen_sigma_inc(0.5) - mid, 3.40, 0.02);
}

UTEST(neutron, deuterium_from_mass) {
    neutron::HydrogenModel h;
    neutron::Scattering s;
    neutron::atom(&s, 1, 2.014, false, h, 1.0);
    EXPECT_NEAR(s.b, neutron::B_D, 1e-9);
    // Volume of a D is that of an H
    neutron::Scattering sh;
    neutron::atom(&sh, 1, 1.008, false, h, 1.0);
    EXPECT_NEAR(s.volume, sh.volume, 1e-12);
    // Exchangeable D follows the reservoir
    h.exchange_d = 0.0;
    neutron::atom(&s, 1, 2.014, true, h, 1.0);
    EXPECT_NEAR(s.b, neutron::B_H, 1e-9);
}
