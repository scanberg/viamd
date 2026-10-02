#pragma once

// GISAXS support for the scattering component: X-ray optics, the scattering geometry, the evaluated q-map (with the
// instrument resolution and the cuts) and the projection onto an area detector. Free of ImGui and ApplicationState, so
// it is linked into viamd_test (see test/CMakeLists.txt) and tested in tests/test_scattering.cpp.
//
// Conventions
// - Lengths in Å, wave vectors in Å^-1, angles in radians. The UI converts to nm^-1 and degrees for display only.
// - Angles are measured in the ambient medium with the vacuum wave number k0 = 2 pi / lambda (as md_gisaxs and
//   BornAgain). The vertical wave vector of the incident beam in the ambient is p = k0 sin(alpha_i), and
//   q_z = k0 (sin alpha_i + sin alpha_f). The sample horizon is alpha_f = 0, i.e. q_z = p.
// - SLDs of the substrate and film are relative to the ambient medium.
// - Particle weights are electrons. The intensity of make_model() is dsigma/dOmega per unit sample area (sr^-1).

#include <md_gisaxs.h>

#include <stddef.h>
#include <stdint.h>
#include <vector>

namespace gisaxs {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDegToRad = kPi / 180.0;
constexpr double kRadToDeg = 180.0 / kPi;
constexpr double kHcKevAngstrom = 12.398419843320026;  // lambda [Å] = hc / E [keV]

// -------------------------------------------------------------------------------------------------
// X-ray optics
// -------------------------------------------------------------------------------------------------

inline double wavelength_from_energy(double energy_kev) { return kHcKevAngstrom / energy_kev; }
inline double energy_from_wavelength(double lambda)     { return kHcKevAngstrom / lambda; }

// Scattering length density (Å^-2) of a medium with electron density rho_e (e/Å^3)
inline double sld_from_electron_density(double rho_e) { return MD_GISAXS_R_E * rho_e; }

// Absorption SLD (Å^-2, imaginary part) for the refractive index n = 1 - delta + i beta
inline double absorption_sld(double beta, double lambda) { return 2.0 * kPi * beta / (lambda * lambda); }

// delta = lambda^2 r_e rho_e / (2 pi)
inline double delta_from_electron_density(double rho_e, double lambda) { return lambda * lambda * MD_GISAXS_R_E * rho_e / (2.0 * kPi); }
inline double electron_density_from_delta(double delta, double lambda) { return 2.0 * kPi * delta / (lambda * lambda * MD_GISAXS_R_E); }

// Critical angle (rad) of total external reflection for an SLD step d_sld (Å^-2, medium below minus the ambient).
// 0 when there is no total reflection (d_sld <= 0).
double critical_angle(double d_sld, double lambda);

// -------------------------------------------------------------------------------------------------
// Scattering geometry
// -------------------------------------------------------------------------------------------------

struct Beam {
    double wavelength = 1.0;    // Å
    double alpha_i = 0.0;       // rad

    double k0() const;          // 2 pi / lambda
    double p() const;           // k0 sin(alpha_i)
};

double qz_from_alpha_f(const Beam& beam, double alpha_f);
double alpha_f_from_qz(const Beam& beam, double qz);       // asin clamped to [-pi/2, pi/2]
inline double horizon_qz(const Beam& beam) { return qz_from_alpha_f(beam, 0.0); }
// Exit angle of the Yoneda peak is the critical angle of the SLD step d_sld (relative to the ambient)
double yoneda_qz(const Beam& beam, double d_sld);

// -------------------------------------------------------------------------------------------------
// Experimental setup -> md_gisaxs model
// -------------------------------------------------------------------------------------------------

struct Setup {
    double wavelength = 1.0;        // Å
    double alpha_i = 0.0;           // rad
    bool   dwba = true;             // false: Born approximation without substrate
    bool   graded = true;           // laterally averaged film in the DWBA reference medium
    double z_substrate = 0.0;       // Å
    double rho_ambient = 0.0;       // e/Å^3
    double rho_substrate = 0.6991;  // e/Å^3
    double beta_substrate = 0.0;
    double roughness = 0.0;         // Å RMS, substrate interface (Nevot-Croce)
    double rho_material = 0.0;      // e/Å^3 of the particle material, <= 0: no displaced ambient (contrast 1)
    double beta_material = 0.0;     // absorption of the particle material at rho_material (graded film only)
};

// Every particle displaces (electrons / rho_material) of the ambient, which scales all amplitudes by
// 1 - rho_ambient / rho_material (a common material density for all particles).
double contrast_amplitude(const Setup& setup);
inline double contrast_factor(const Setup& setup) { const double c = contrast_amplitude(setup); return c * c; }

// SLD steps relative to the ambient (Å^-2)
double substrate_sld(const Setup& setup);
double substrate_absorption_sld(const Setup& setup);

// Model for particle weights in electrons: I = r_e^2 contrast^2 <|F|^2> / A (sr^-1). The film SLD is
// r_e contrast * profile (profile in e/Å^3, see md_gisaxs_slice_profile).
md_gisaxs_model_t make_model(const Setup& setup);

// Mean film SLD (Å^-2, relative to the ambient) over the slices with a density above half the maximum, or 0
double film_sld(const double* profile, size_t num_slices, double profile_sld_scale);

// -------------------------------------------------------------------------------------------------
// q-map
// -------------------------------------------------------------------------------------------------

// Intensity on md_gisaxs rings (non-uniform in q_par, see md_gisaxs_ring_edges) times uniform q_z rows.
struct QMap {
    std::vector<double>   ring_q;       // mean |q_par| per ring, increasing
    std::vector<double>   ring_edge;    // num_rings + 1 ring boundaries
    std::vector<unsigned> ring_count;   // lattice points (full plane) per ring
    double qz0 = 0.0;                   // q_z of row 0
    double dqz = 0.0;                   // row spacing
    size_t rows = 0;
    std::vector<float> I;               // [row][ring]

    size_t cols() const { return ring_q.size(); }
    bool   empty() const { return rows == 0 || ring_q.empty(); }
    double qz(size_t row) const { return qz0 + (double)row * dqz; }
    float  at(size_t row, size_t ring) const { return I[row * ring_q.size() + ring]; }
    double qpar_min() const { return ring_edge.empty() ? 0.0 : ring_edge.front(); }
    double qpar_max() const { return ring_edge.empty() ? 0.0 : ring_edge.back(); }
};

// Rings, edges and counts of a computed context, rows = 0 and no intensity
void qmap_init_rings(QMap* map, const md_gisaxs_t* ctx);

// Ring with edge[r] <= q_par < edge[r + 1] (the last ring includes its upper edge), -1 outside [edge[0], edge[R]]
long ring_at(const QMap& map, double qpar);

// Linear interpolation in q_par between ring centers: value = (1 - t) I[c0] + t I[c1]. Between the outer edges and the
// first / last ring center the nearest ring is used. False outside [edge[0], edge[R]].
bool ring_interp(const QMap& map, double qpar, size_t* c0, size_t* c1, double* t);

// Bilinear lookup at (q_par, q_z). False outside the map.
bool qmap_sample(const QMap& map, double qpar, double qz, float* out);

// Gaussian instrument resolution (FWHM in Å^-1, <= 0 disables a direction). Along q_par the convolution uses the
// non-uniform rings as quadrature nodes (weight = ring width). Rows with q_z < qz_valid_min (below the sample
// horizon) carry no intensity: they are excluded from the q_z convolution and set to zero.
void qmap_apply_resolution(QMap* dst, const QMap& src, double fwhm_qpar, double fwhm_qz, double qz_valid_min);

// Horizontal cut I(q_par): mean over the rows within |q_z - qz| <= max(width / 2, dqz / 2). out: cols() values.
// Returns the number of rows averaged.
size_t qmap_cut_horizontal(const QMap& map, double qz, double width, double* out);

// Vertical cut I(q_z): mean over the rings with |ring_q - qpar| <= width / 2, or the interpolation between the two
// neighbouring rings when no ring is inside the band. out: rows values. Returns the number of rings used (0: none).
size_t qmap_cut_vertical(const QMap& map, double qpar, double width, double* out);

// -------------------------------------------------------------------------------------------------
// Area detector (flat, normal to the direct beam)
// -------------------------------------------------------------------------------------------------

struct Detector {
    double sdd_mm = 4976.0;         // sample-detector distance
    double pixel_mm = 0.172;
    int    npx_h = 981;
    int    npx_v = 1043;
    double beam_x_px = 490.5;       // direct beam column, from the left edge
    double beam_y_px = 23.0;        // direct beam row, from the bottom edge
    int    binning = 4;
    bool   gaps = true;
    int    module_w = 487, module_h = 195;
    int    gap_w = 7, gap_h = 17;
    int    bs_direct_px = 20;       // square beamstops, full width
    int    bs_specular_px = 12;
};

struct DetectorQ {
    double qpar;    // |q_par| = sqrt(q_x^2 + q_y^2)
    double qy;
    double qz;
    double alpha_f; // rad
};

// Scattering vector of a detector position (mm, relative to the direct beam, y upwards).
// Returns false below the sample horizon (alpha_f < 0); out is filled regardless.
bool detector_q(const Detector& det, const Beam& beam, double x_mm, double y_mm, DetectorQ* out);

// Pixel (px, py from the bottom-left corner) in a module gap or behind a beamstop
bool detector_masked(const Detector& det, const Beam& beam, int px, int py);

struct DetectorImage {
    int width = 0, height = 0;      // binned pixels
    std::vector<float> I;           // [row 0 = bottom][column], < 0: masked or outside the q-map
    double qy_min = 0, qy_max = 0;  // linear axis bounds along the central row / column, for display
    double qz_min = 0, qz_max = 0;
};

void detector_render(DetectorImage* img, const Detector& det, const Beam& beam, const QMap& map);

// Binned pixel -> the detector position (mm) of its center relative to the direct beam
void detector_bin_center_mm(const Detector& det, int bx, int by, double* x_mm, double* y_mm);

}  // namespace gisaxs
