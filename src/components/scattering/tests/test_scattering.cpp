#include "utest.h"

#include "../scattering_core.h"

#include <math.h>
#include <vector>

// GISAXS support of the scattering component: optics, geometry, the q-map operations (lookup on the non-uniform
// rings, instrument resolution, cuts) and the detector projection. The physics of md_gisaxs itself is validated in
// mdlib (unittest/test_gisaxs.c and tools/gisaxs_validation against BornAgain, sasmodels and refnx).

using namespace gisaxs;

namespace {

// Map with rings at non-uniform q (like lattice shells) and a linear intensity in q_par and q_z
QMap make_map(size_t rows = 11) {
    QMap m;
    m.ring_q    = {1.0, 1.4142, 2.0, 2.2361, 3.0, 4.0};
    m.ring_edge = {0.5, 1.2071, 1.7071, 2.1180, 2.6180, 3.5, 4.5};
    m.ring_count = {4, 4, 4, 8, 12, 20};
    m.qz0 = 0.0;
    m.dqz = 0.1;
    m.rows = rows;
    m.I.resize(rows * m.cols());
    for (size_t r = 0; r < rows; ++r)
        for (size_t c = 0; c < m.cols(); ++c) m.I[r * m.cols() + c] = (float)(2.0 * m.ring_q[c] + 10.0 * m.qz(r));
    return m;
}

}  // namespace

UTEST(scattering, critical_angle_silicon) {
    // Si (0.6991 e/Å^3) at 12 keV: alpha_c ~ 0.148 deg, delta ~ 3.4e-6
    const double lambda = wavelength_from_energy(12.0);
    EXPECT_NEAR(lambda, 1.0332, 1e-4);
    const double ac = critical_angle(sld_from_electron_density(0.6991), lambda) * kRadToDeg;
    EXPECT_NEAR(ac, 0.148, 0.002);
    // alpha_c = sqrt(2 delta) for small angles
    const double delta = delta_from_electron_density(0.6991, lambda);
    EXPECT_NEAR(ac * kDegToRad, sqrt(2.0 * delta), 1e-6);
    EXPECT_NEAR(electron_density_from_delta(delta, lambda), 0.6991, 1e-9);
    EXPECT_EQ(critical_angle(-1e-6, lambda), 0.0);
}

UTEST(scattering, geometry_round_trip) {
    Beam b;
    b.wavelength = 1.0;
    b.alpha_i = 0.2 * kDegToRad;
    for (double af_deg : {0.0, 0.05, 0.3, 1.5}) {
        const double qz = qz_from_alpha_f(b, af_deg * kDegToRad);
        EXPECT_NEAR(alpha_f_from_qz(b, qz) * kRadToDeg, af_deg, 1e-9);
    }
    EXPECT_NEAR(horizon_qz(b), b.p(), 1e-15);
    const double d_sld = sld_from_electron_density(0.6991);
    EXPECT_NEAR(yoneda_qz(b, d_sld), qz_from_alpha_f(b, critical_angle(d_sld, 1.0)), 1e-15);
}

UTEST(scattering, setup_model) {
    Setup s;
    s.wavelength = 1.0;
    s.rho_ambient = 0.3342;     // water
    s.rho_material = 0.478;     // cellulose
    s.rho_substrate = 0.6991;
    s.beta_substrate = 3.6e-8;
    const double c = 1.0 - 0.3342 / 0.478;
    EXPECT_NEAR(contrast_amplitude(s), c, 1e-12);
    const md_gisaxs_model_t m = make_model(s);
    EXPECT_NEAR(m.intensity_scale, MD_GISAXS_R_E * MD_GISAXS_R_E * c * c, 1e-24);
    EXPECT_NEAR(m.sld_substrate - m.sld_ambient, substrate_sld(s), 1e-18);
    EXPECT_NEAR(m.sld_substrate_abs, 2.0 * kPi * 3.6e-8, 1e-18);
    // Contrast matched: zero intensity, not md_gisaxs' default scale of 1
    s.rho_ambient = s.rho_material;
    EXPECT_TRUE(make_model(s).intensity_scale > 0.0);
    EXPECT_TRUE(make_model(s).intensity_scale < 1e-300);
}

UTEST(scattering, film_sld) {
    // Slices at or above half the maximum: 0.4, 0.4, 0.4, 0.2
    const double prof[] = {0.0, 0.1, 0.4, 0.4, 0.4, 0.2, 0.0};
    EXPECT_NEAR(film_sld(prof, 7, 2.0), 2.0 * 0.35, 1e-12);
    EXPECT_EQ(film_sld(prof, 0, 2.0), 0.0);
}

UTEST(scattering, ring_lookup) {
    const QMap m = make_map();
    EXPECT_EQ(ring_at(m, 0.49), -1L);
    EXPECT_EQ(ring_at(m, 0.5), 0L);
    EXPECT_EQ(ring_at(m, 1.2070), 0L);
    EXPECT_EQ(ring_at(m, 1.2071), 1L);
    EXPECT_EQ(ring_at(m, 4.5), 5L);     // upper edge belongs to the last ring
    EXPECT_EQ(ring_at(m, 4.51), -1L);

    size_t c0, c1;
    double t;
    ASSERT_TRUE(ring_interp(m, 2.1, &c0, &c1, &t));
    EXPECT_EQ(c0, (size_t)2);
    EXPECT_EQ(c1, (size_t)3);
    EXPECT_NEAR(t, (2.1 - 2.0) / (2.2361 - 2.0), 1e-12);
    ASSERT_TRUE(ring_interp(m, 0.7, &c0, &c1, &t));     // below the first ring center: nearest ring
    EXPECT_EQ(c0, (size_t)0);
    EXPECT_EQ(c1, (size_t)0);
    EXPECT_FALSE(ring_interp(m, 0.2, &c0, &c1, &t));    // specular rod, no data

    // The map is linear in q_par between ring centers and in q_z: bilinear lookup is exact there
    float v;
    ASSERT_TRUE(qmap_sample(m, 2.6, 0.37, &v));
    EXPECT_NEAR(v, 2.0 * 2.6 + 10.0 * 0.37, 1e-5);
    ASSERT_TRUE(qmap_sample(m, 3.0, 1.0, &v));          // last row
    EXPECT_NEAR(v, 2.0 * 3.0 + 10.0, 1e-5);
    EXPECT_FALSE(qmap_sample(m, 3.0, 1.01, &v));
    EXPECT_FALSE(qmap_sample(m, 3.0, -0.01, &v));
}

UTEST(scattering, resolution) {
    QMap m = make_map();
    // A constant map is unchanged (normalized kernels), also with the non-uniform rings
    for (float& x : m.I) x = 3.0f;
    QMap s;
    qmap_apply_resolution(&s, m, 0.8, 0.25, -1.0);
    for (float x : s.I) EXPECT_NEAR(x, 3.0f, 1e-5f);
    // No resolution: identity
    m = make_map();
    qmap_apply_resolution(&s, m, 0.0, 0.0, -1.0);
    for (size_t i = 0; i < m.I.size(); ++i) EXPECT_EQ(s.I[i], m.I[i]);
    // Rows below the horizon are zero and do not leak into the rows above
    qmap_apply_resolution(&s, m, 0.0, 0.3, 0.45);
    for (size_t c = 0; c < m.cols(); ++c) {
        EXPECT_EQ(s.at(4, c), 0.0f);
        EXPECT_TRUE(s.at(5, c) > m.at(5, c));       // only rows above (larger I) contribute at the edge
    }
    // Smoothing along q_par conserves a linear profile away from the ends
    qmap_apply_resolution(&s, m, 0.3, 0.0, -1.0);
    EXPECT_NEAR(s.at(3, 2), m.at(3, 2), 0.2f);
}

UTEST(scattering, cuts) {
    const QMap m = make_map();
    std::vector<double> h(m.cols()), v(m.rows);
    EXPECT_EQ(qmap_cut_horizontal(m, 0.5, 0.0, h.data()), (size_t)1);
    EXPECT_NEAR(h[4], 2.0 * 3.0 + 5.0, 1e-5);
    EXPECT_EQ(qmap_cut_horizontal(m, 0.5, 0.25, h.data()), (size_t)3);     // rows 0.4, 0.5, 0.6
    EXPECT_NEAR(h[0], 2.0 + 5.0, 1e-5);
    EXPECT_EQ(qmap_cut_vertical(m, 2.1, 0.5, v.data()), (size_t)2);        // rings 2.0 and 2.2361
    EXPECT_NEAR(v[0], (4.0 + 4.4722) / 2.0, 1e-4);
    // Narrow band between two rings: interpolation instead of an empty cut
    EXPECT_EQ(qmap_cut_vertical(m, 2.1, 0.0, v.data()), (size_t)2);
    EXPECT_NEAR(v[2], 2.0 * 2.1 + 2.0, 1e-5);
    EXPECT_EQ(qmap_cut_vertical(m, 0.1, 0.0, v.data()), (size_t)0);
}

UTEST(scattering, detector_geometry) {
    Detector d;
    d.gaps = false;
    Beam b;
    b.wavelength = 0.97;
    b.alpha_i = 0.42 * kDegToRad;
    DetectorQ q;
    // Direct beam: q = 0 but below the horizon
    EXPECT_FALSE(detector_q(d, b, 0.0, 0.0, &q));
    EXPECT_NEAR(q.qpar, 0.0, 1e-12);
    EXPECT_NEAR(q.qz, 0.0, 1e-12);
    // Specular spot: q_par = 0, q_z = 2 k sin(alpha_i), alpha_f = alpha_i
    const double ys = tan(2.0 * b.alpha_i) * d.sdd_mm;
    EXPECT_TRUE(detector_q(d, b, 0.0, ys, &q));
    EXPECT_NEAR(q.qpar, 0.0, 1e-9);
    EXPECT_NEAR(q.qz, 2.0 * b.p(), 1e-9);
    EXPECT_NEAR(q.alpha_f, b.alpha_i, 1e-9);
    // Horizontal offset at the horizon height: q_y = k sin(2 theta_h)
    const double yh = tan(b.alpha_i) * d.sdd_mm;
    EXPECT_TRUE(detector_q(d, b, 50.0, yh, &q));
    EXPECT_NEAR(q.alpha_f, 0.0, 1e-4);
    EXPECT_NEAR(q.qy, b.k0() * 50.0 / sqrt(50.0 * 50.0 + d.sdd_mm * d.sdd_mm / (cos(b.alpha_i) * cos(b.alpha_i))), 1e-5);
    // Masks
    EXPECT_TRUE(detector_masked(d, b, (int)d.beam_x_px, (int)d.beam_y_px));
    EXPECT_TRUE(detector_masked(d, b, (int)d.beam_x_px, (int)(d.beam_y_px + ys / d.pixel_mm)));
    EXPECT_FALSE(detector_masked(d, b, 10, 500));
    d.gaps = true;
    EXPECT_TRUE(detector_masked(d, b, d.module_w + 1, 500));
}

UTEST(scattering, detector_render) {
    // A map that is 1 everywhere above the horizon renders 1 on every unmasked pixel inside its q-range
    Beam b;
    b.wavelength = 1.0;
    b.alpha_i = 0.3 * kDegToRad;
    QMap m;
    m.ring_q = {0.01, 0.02, 0.04};
    m.ring_edge = {0.005, 0.015, 0.03, 0.05};
    m.ring_count = {4, 8, 12};
    m.qz0 = 0.0;
    m.dqz = 0.001;
    m.rows = 101;
    m.I.assign(m.rows * m.cols(), 1.0f);
    Detector d;
    d.npx_h = 200; d.npx_v = 200; d.beam_x_px = 100; d.beam_y_px = 20; d.binning = 2; d.gaps = false;
    d.sdd_mm = 3000; d.pixel_mm = 0.172;
    DetectorImage img;
    detector_render(&img, d, b, m);
    EXPECT_EQ(img.width, 100);
    EXPECT_EQ(img.height, 100);
    size_t ones = 0, masked = 0;
    for (float v : img.I) { if (v < 0.0f) masked += 1; else { EXPECT_NEAR(v, 1.0f, 1e-6f); ones += 1; } }
    EXPECT_TRUE(ones > 0);
    EXPECT_TRUE(masked > 0);    // below the horizon, the specular rod (q_par < edge[0]) and the beamstops
    EXPECT_TRUE(img.qy_min < 0.0 && img.qy_max > 0.0);
}
