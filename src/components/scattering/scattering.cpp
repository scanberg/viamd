// Scattering component: grazing incidence small-angle X-ray scattering (GISAXS)
//
// Computes the rotationally (azimuthally) averaged grazing incidence intensity I(q_par, q_z) of the current system
// state. The electron density is approximated with one Gaussian per particle (bead / atom), or with a cross-section
// swept along coarse grained fibrils (fibril_core.h). The heavy lifting is done in mdlib (md_gisaxs, see md_gisaxs.h
// for the method and tools/gisaxs_validation for its validation against BornAgain, sasmodels and refnx); the X-ray
// optics, the scattering geometry, the q-map operations and the detector projection are in scattering_core.h.
//
// The computation has two stages:
//   1. Structure stage (heavy): slices, FFTs and ring averaged cross spectral matrices. Depends on the particle
//      coordinates, the selection, the particle model (electrons, Gaussian width, density model) and the sampling.
//   2. Model stage (cheap): DWBA / Born evaluation on top of stage 1. Depends on the beam, the substrate and the
//      ambient medium, and is re-evaluated automatically when any of those change. The instrument resolution, the
//      detector image and the cuts are derived from its result without touching md_gisaxs.
//
// The ambient medium enters as the reference medium of the DWBA and as a uniform contrast factor of the particles:
// every particle displaces (electrons / material electron density) of the ambient, which scales the intensity by
// (1 - rho_ambient / rho_material)^2.
//
// Neutrons (GISANS): the neutron path (neutron_core.h: excess scattering lengths, H/D exchange, incoherent background)
// is compiled with VIAMD_SCATTERING_NEUTRON = 1. It is off for now while the X-ray path is being consolidated.

#ifndef VIAMD_SCATTERING_NEUTRON
#define VIAMD_SCATTERING_NEUTRON 0
#endif

#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <core/md_hash.h>
#include <core/md_fft.h>

#include <md_gisaxs.h>
#include <md_filter.h>
#include <md_system.h>
#include <md_util.h>
#include <md_csv.h>

#include "scattering_core.h"
#include "fibril_core.h"
#if VIAMD_SCATTERING_NEUTRON
#include "neutron_core.h"
#endif

#include <viamd_event.h>
#include <viamd.h>
#include <event.h>
#include <task_system.h>
#include <serialization_utils.h>

#include <imgui_widgets.h>
#include <implot_internal.h>

#include <atomic>
#include <algorithm>
#include <utility>
#include <vector>
#include <math.h>
#include <float.h>

namespace {

constexpr double DEG_TO_RAD = gisaxs::kDegToRad;   // PI comes from md_common.h
constexpr double RAD_TO_DEG = gisaxs::kRadToDeg;

// Electron densities in e/Å^3
struct MediumPreset {
    const char* name;
    double electron_density;
    double beta;    // Absorption (imaginary part of the refractive index), substrates only
};

// Ambient media
const MediumPreset ambient_presets[] = {
    {"Vacuum",  0.0,      0.0},
    {"Air",     3.62e-4,  0.0},     // 1.205 kg/m^3, Z/A ~ 0.499
    {"Water",   0.3342,   0.0},
    {"Custom",  0.0,      0.0},
};
constexpr int AMBIENT_CUSTOM = 3;

// Substrates, beta given for ~12 keV
const MediumPreset substrate_presets[] = {
    {"Silicon",     0.6991, 3.6e-8},  // 2.329 g/cm^3, mu/rho ~19 cm^2/g at 12 keV
    {"SiO2 (glass)",0.6620, 2.4e-8},  // 2.2 g/cm^3, mu/rho ~13 cm^2/g at 12 keV
    {"Custom",      0.6991, 0.0},
};
constexpr int SUB_CUSTOM = 2;

// Experimental setups. Applying a preset sets the beam, the substrate optical constants, the q-range covered by the
// detector, the detector and the default cut positions. Substrate optical constants are given as delta / beta at the
// preset's wavelength.
struct BeamPreset {
    const char* name;
    const char* description;
    double wavelength_nm;
    double alpha_i_deg;
    // Detector (flat, normal to the direct beam)
    double sdd_mm;
    double pixel_mm;
    int    det_px_h;            // horizontal pixels
    double alpha_f_max_deg;     // vertical range used, from the sample horizon
    // Substrate
    double sub_delta;
    double sub_beta;
    // Cuts
    double cut_alpha_f_deg;     // horizontal cut
    double cut_qpar_nm;         // vertical cut
    // Detector image
    int    det_px_v;            // vertical pixels
    double beam_y_px;           // direct beam row, counted from the bottom edge
    double beam_x_px;           // direct beam column, counted from the left edge
    int    module_w, module_h;  // detector module size (px), 0 = no gaps
    int    gap_w, gap_h;        // gaps between modules (px)
    int    bs_direct_px;        // direct beamstop (square, px)
    int    bs_specular_px;      // specular beamstop (square, px)
    // Film material
    double material_beta;
};

const BeamPreset beam_presets[] = {
    {
        "P03 / PETRA III (Ohm 2018)",
        "Ohm et al., J. Coat. Technol. Res. 15, 759 (2018), as used in cnf_to_bornagain_simple.py:\n"
        "lambda = 0.097 nm, alpha_i = 0.42 deg, Pilatus 1M (981 x 1043 px, 172 um) at SDD = 4976 mm,\n"
        "alpha_f up to 1.6 deg, Si substrate (delta 2.979e-6, beta 2.78e-8).\n"
        "Horizontal cut at the cellulose Yoneda (alpha_f = 0.117 deg), vertical cut at q_y = 0.05 nm^-1.",
        0.097, 0.42,
        4976.0, 0.172, 981, 1.6,
        2.979e-6, 2.78e-8,
        0.117, 0.05,
        // Pilatus 1M: 2 x 5 modules of 487 x 195 px, gaps 7 px (horizontal) and 17 px (vertical).
        // The direct beam row is chosen such that the top edge is at alpha_f ~ 1.6 deg (as in the script).
        // Beamstop sizes are estimates, the paper does not state them.
        1043, 23.0, 490.5,
        487, 195, 7, 17,
        20, 12,
        2.0e-9,
    },
    {
        "P03 / PETRA III (Brett 2019)",
        "Brett et al., Macromolecules 52, 4721 (2019), spray-deposited CNF films:\n"
        "13 keV (lambda = 0.09537 nm), alpha_i = 0.41 deg, Pilatus 1M (981 x 1043 px, 172 um) at SDD = 5007 mm,\n"
        "alpha_f up to 1.61 deg, Si substrate (delta 2.88e-6, beta 2.6e-8, critical angle 0.137 deg).\n"
        "Horizontal cut at the CNF Yoneda (alpha_f = 0.13 deg). The vertical cut of the paper is at q_y = 0, on the\n"
        "specular rod, which the model does not contain: it is placed at the lowest q_par of the map instead.\n"
        "Direct beam position and beamstop sizes are estimated from the Pilatus 1M patterns (SI, Fig. S8).",
        0.09537, 0.41,
        5007.0, 0.172, 981, 1.614,
        // Si at 13 keV. delta reproduces the critical angle given in the paper (0.137 deg); delta and beta are
        // scaled from the values of the Ohm 2018 setup (delta ~ lambda^2, beta ~ lambda^4). The native oxide is not
        // modelled, the paper finds the substrate contribution negligible.
        2.88e-6, 2.6e-8,
        // The paper cuts at qy = 0 nm^-1 (Fig. 4b); 0 is clamped to the lowest q_par of the map
        0.13, 0.0,
        // Not stated in the paper, estimated from Fig. S8, calibrated on the two horizontal module gaps visible
        // there (212 px pitch): the specular beamstop sits on the gap at rows 407-423 and the CNF Yoneda is ~135 px
        // below that gap, which puts the direct beam at row ~14 (+-5) and the top edge at alpha_f ~1.61 deg, as in
        // Ohm 2018. The specular rod runs ~41 px right of the vertical module gap (column ~532). The specular
        // beamstop is a disc of ~33 px (~5.7 mm), modelled as a square. The direct beamstop is outside the figure.
        1043, 14.0, 532.0,
        487, 195, 7, 17,
        20, 33,
        // Cellulose at 13 keV, scaled from 2e-9 at 12.8 keV
        1.9e-9,
    },
};

enum CutMode : int {
    CutMode_Yoneda = 0,     // Horizontal cut follows the substrate Yoneda position
    CutMode_AlphaF = 1,     // Horizontal cut at a fixed exit angle
    CutMode_Qz     = 2,     // Horizontal cut at a fixed q_z
    CutMode_FilmYoneda = 3, // Horizontal cut follows the film Yoneda position (graded DWBA)
};
const char* cut_mode_lbl[] = {"Substrate Yoneda", "Fixed exit angle", "Fixed q_z", "Film Yoneda"};

enum ElectronMode : int { ElectronMode_Constant = 0, ElectronMode_AtomicNumber = 1, ElectronMode_Mass = 2 };
const char* electron_mode_lbl[] = {"Constant", "Atomic number (Z)", "From mass"};

enum SigmaMode : int { SigmaMode_Constant = 0, SigmaMode_Radius = 1 };
const char* sigma_mode_lbl[] = {"Constant", "From radius (R / sqrt(5))"};

enum DensityModel : int { DensityModel_Beads = 0, DensityModel_Fibril = 1 };
const char* density_model_lbl[] = {"Beads (Gaussian per particle)", "Fibrils (swept cross-section)"};

enum FibrilTemplateMode : int { FibrilTemplate_Beads = 0, FibrilTemplate_Disk = 1, FibrilTemplate_File = 2 };
const char* fibril_template_lbl[] = {"Bead layout", "Uniform disk", "File"};

#if VIAMD_SCATTERING_NEUTRON
enum Radiation : int { Radiation_Xray = 0, Radiation_Neutron = 1 };
const char* radiation_lbl[] = {"X-rays (GISAXS)", "Neutrons (GISANS)"};

// Neutron media, SLD and absorption SLD (imaginary part) in 10^-6 Å^-2. The absorption SLD is wavelength
// independent (sigma_a ~ lambda), see neutron::absorption_sld. Incoherent attenuation is not included.
struct SldPreset {
    const char* name;
    double sld;
    double abs;
};

const SldPreset neutron_background_presets[] = {
    {"Vacuum / air",     0.0, 0.0},
    {"Water (H2O/D2O)",  0.0, 0.0},    // SLD from the D2O fraction
    {"Custom",           0.0, 0.0},
};
constexpr int NBG_WATER = 1;
constexpr int NBG_CUSTOM = 2;

const SldPreset neutron_substrate_presets[] = {
    {"Silicon",             2.07, 2.4e-5},  // 2.329 g/cm^3
    {"SiO2 (fused silica)", 3.47, 1.1e-5},  // 2.2 g/cm^3
    {"Quartz",              4.18, 1.3e-5},  // 2.65 g/cm^3
    {"Sapphire (Al2O3)",    5.71, 3.0e-5},  // 3.98 g/cm^3
    {"Gold",                4.50, 1.6e-2},  // 19.32 g/cm^3
    {"Custom",              2.07, 0.0},
};
constexpr int NSUB_CUSTOM = 5;

enum NeutronWeight : int { NeutronWeight_Constant = 0, NeutronWeight_Element = 1 };
const char* neutron_weight_lbl[] = {"Constant per particle", "Per element (Sears 1992)"};
#endif

enum ComputeState : int {
    ComputeState_Idle = 0,
    ComputeState_Running,
    ComputeState_Done,
    ComputeState_Failed,
};

// Combo over a fixed label array. Clamps the value into range. Returns true when the selection changed.
template <size_t N>
bool combo(const char* label, int* value, const char* const (&items)[N]) {
    *value = CLAMP(*value, 0, (int)N - 1);
    bool changed = false;
    if (ImGui::BeginCombo(label, items[*value])) {
        for (int i = 0; i < (int)N; ++i) {
            if (ImGui::Selectable(items[i], *value == i)) {
                changed = *value != i;
                *value = i;
            }
        }
        ImGui::EndCombo();
    }
    return changed;
}

// Index range [beg, end) of strictly positive values, for log scale plotting
void positive_range(const double* v, size_t n, size_t* beg, size_t* end) {
    size_t b = 0;
    while (b < n && !(v[b] > 0.0)) ++b;
    size_t e = n;
    while (e > b && !(v[e - 1] > 0.0)) --e;
    *beg = b;
    *end = e;
}

}  // namespace

struct ScatteringComponent : viamd::EventHandler {
    ApplicationState* app_state = nullptr;
    bool show_window = false;

    // --- Sample: selection and particle model ---
    char  filter[256] = "all";
    char  filter_err[256] = "";
    bool  filter_valid = true;
    int   electron_mode = ElectronMode_Constant;
    float electrons_per_particle = 1900.0f;     // e.g. ~22 glucose units (86 e each) per sCG bead
    float electrons_per_dalton = 0.530f;         // cellulose C6H10O5: 86 e / 162.14 Da
    int   sigma_mode = SigmaMode_Constant;
    float sigma_constant = 6.0f;                 // Å
    float material_density = 0.478f;             // e/Å^3 (cellulose ~1.5 g/cm^3)
    double material_beta = 2.0e-9;               // absorption of the particle material at material_density

    // --- Sample: density model ---
    int   density_model = DensityModel_Beads;
    char  fib_center_name[16] = "CC";
    int   fib_stride = 7;
    int   fib_template = FibrilTemplate_Beads;
    float fib_bead_sigma = 6.0f;                 // Å, cross-section Gaussians of the bead layout template
    float fib_disk_radius = 19.1f;               // Å
    float fib_disk_spacing = 8.0f;               // Å
    float fib_ds = 0.0f;                         // Å, sample spacing along the fibril, 0 = smallest template sigma
    char  fib_template_path[512] = "";
    char  fib_info[256] = "";

    // Self scattering of the (original) beads, groups of (sigma, sum w^2), for the diagnostic in the horizontal cut
    std::vector<std::pair<float, double>> bead_self;
    double bead_self_area = 1.0;
    bool   show_bead_self = false;

    // --- Ambient medium and substrate ---
    int    ambient_preset = 0;
    double ambient_density = 0.0;                // e/Å^3
    bool   dwba = true;
    int    sub_preset = 0;
    double sub_density = substrate_presets[0].electron_density;
    double sub_beta = substrate_presets[0].beta;
    bool   sub_auto_z = true;
    float  sub_z = 0.0f;                         // Å
    float  sub_z_offset = 0.0f;                  // Å, offset relative to the lowest particle when auto
    float  sub_roughness = 3.0f;                 // Å RMS
    bool   film_graded = true;                   // Laterally averaged film as part of the DWBA reference medium

    // --- Beam ---
    float energy_kev = 12.0f;
    float alpha_i_deg = 0.20f;
    int   beam_preset = -1;                      // Last applied preset, -1 = none

    // --- Sampling (q in nm^-1 in the UI) ---
    float q_par_max_nm = 2.0f;
    float q_z_max_nm = 2.0f;
    float oversampling = 3.0f;
    float ring_rel_width = 0.1f;                 // Lattice shells closer than this (relative to q_par) share a ring
    int   num_qz = 512;
    int   max_slices = 1024;

    // --- Detector and resolution ---
    float res_fwhm_qpar = 0.0f;                  // Resolution (FWHM, nm^-1)
    float res_fwhm_qz = 0.0f;
    bool  view_detector = false;                 // Show the result as a detector image
    gisaxs::Detector det;

    // --- Display ---
    bool  log_scale = true;
    float decades = 6.0f;
    ImPlotColormap colormap = ImPlotColormap_Jet;
    bool  show_profile = false;

    // --- Cuts (positions in nm^-1) ---
    bool   show_cuts = true;
    double cut_qz = 0.0;          // Horizontal cut, I(q_par) at fixed q_z
    double cut_qpar = 0.0;        // Vertical cut, I(q_z) at fixed q_par
    float  cut_width_qz = 0.0f;   // Integration band width (nm^-1), 0 -> single row
    float  cut_width_qpar = 0.0f; // Integration band width (nm^-1), 0 -> single ring
    int    cut_mode = CutMode_Yoneda;
    float  cut_alpha_f_deg = 0.117f;  // Used when cut_mode == CutMode_AlphaF
    bool   cut_log_x = false;
    bool   vcut_vs_alpha_f = false;   // Plot the vertical cut against alpha_f instead of q_z
    bool   cut_init = false;

#if VIAMD_SCATTERING_NEUTRON
    // --- Neutrons ---
    int   radiation = Radiation_Xray;
    int   ctx_radiation = Radiation_Xray;        // Radiation the current ctx (particle weights) was computed for
    // Particle model (defaults: one sCG bead of 22 anhydroglucose units, cellulose at 1.5 g/cm^3)
    int   neu_weight_mode = NeutronWeight_Constant;
    float neu_b_protiated = 693.0f;              // fm per particle with every H as 1H (22 x 31.50 fm)
    float neu_n_h = 154.0f;                      // non-exchangeable H per particle (22 x 7)
    float neu_n_ex = 66.0f;                      // exchangeable H per particle (22 x 3 hydroxyl)
    float neu_volume = 3949.0f;                  // Å^3 per particle (22 x 179.5 Å^3)
    float neu_mass_density = 1.5f;               // g/cm^3, per element: particle volume = mass / density
    float neu_deuteration = 0.0f;                // D fraction of the non-exchangeable H
    bool  neu_exchange = true;                   // H bound to N, O, S exchange with the reservoir
    bool  neu_exchange_follow = true;            // Reservoir = background water (D2O fraction)
    float neu_exchange_d = 0.0f;                 // D fraction of the reservoir when not following the background
    bool  neu_incoherent = true;                 // Add the flat incoherent background
    // Media (SLD in 10^-6 Å^-2) and beam
    int    neu_bg_preset = 0;
    double neu_bg_sld = 0.0;                     // Custom background
    float  neu_d2o = 1.0f;                       // D2O volume fraction of a water background
    int    neu_sub_preset = 0;
    double neu_sub_sld = neutron_substrate_presets[0].sld;
    double neu_sub_abs = neutron_substrate_presets[0].abs;
    float  neu_wavelength = 6.0f;                // Å
    // Summary of the selection at the last compute
    double neu_sum_b = 0.0;                      // fm, coherent scattering length (not excess)
    double neu_sum_v = 0.0;                      // Å^3
    double neu_sum_inc = 0.0;                    // barn
    double neu_bg_sld_used = 0.0;                // 10^-6 Å^-2, background SLD the weights were computed with
    size_t neu_num_unknown = 0;                  // particles with an element not in the table (b = 0)
    size_t neu_num_exchangeable = 0;
    double inc_level = 0.0;                      // Flat incoherent background (sr^-1) of the current ctx
    double eval_inc = 0.0;                       // Incoherent background of the evaluation in flight
#endif

    // --- Structure stage ---
    md_gisaxs_t* ctx = nullptr;              // Owned. Valid for evaluation when compute_state == Done
    std::atomic<int> compute_state{ComputeState_Idle};
    std::atomic<bool> compute_cancel{false};
    task_system::ID task_slices = task_system::INVALID_ID;
    task_system::ID task_rings  = task_system::INVALID_ID;
    task_system::ID task_finish = task_system::INVALID_ID;
    md_array(void*) scratch = nullptr;       // per worker thread slice scratch
    size_t scratch_bytes = 0;
    uint64_t structure_hash = 0;             // Hash of the inputs used for the current ctx
    uint64_t ctx_gen = 0;                    // Incremented for every created ctx, identifies it in the model hash
    bool structure_stale = false;            // The system changed since the current ctx was computed
    double ctx_q_z_max = 0.0;                // Å^-1, as used when ctx was created
    double particle_z_min = 0.0;
    char status[256] = "";
    double compute_time_start = 0.0;
    double compute_time = 0.0;

    // --- Model stage ---
    task_system::ID task_eval_range = task_system::INVALID_ID;   // The evaluation itself
    task_system::ID task_eval = task_system::INVALID_ID;         // Completion of the evaluation (sets eval_ready)
    std::atomic<bool> eval_ready{false};
    uint64_t eval_hash_pending = 0;          // Hash of the model currently being evaluated
    uint64_t eval_hash_shown = 0;            // Hash of the model currently shown
    md_gisaxs_model_t eval_model = {};       // Model of the evaluation in flight, then of the shown result
    std::vector<double> eval_qz;             // Å^-1
    std::vector<float>  eval_out;            // num_qz * num_rings, written by the evaluation task

    // --- Shown result ---
    gisaxs::QMap raw;                        // As evaluated
    gisaxs::QMap smooth;                     // With the instrument resolution
    bool   res_dirty = false;                // smooth / display need to be rebuilt
    uint64_t res_version = 0;
    double lattice_dq = 0.0;                 // Å^-1, max(2 pi / Lx, 2 pi / Ly), the widest ring
    double qpar_first = 0.0;                 // Å^-1, smallest sampled q_par, min(2 pi / Lx, 2 pi / Ly)
    // Display (nm^-1)
    std::vector<double> map_display;         // [ring][row 0 = highest q_z], log10 if log_scale
    std::vector<double> ring_q_nm;           // Mean q_par per ring, increasing but non-uniform
    std::vector<double> ring_edge_nm;        // num_rings + 1
    std::vector<double> qz_nm;               // Row q_z, ascending
    double res_min = 0.0, res_max = 1.0;     // Color range
    double res_qpar_min = 0.0, res_qpar_max = 1.0;
    double res_qz_min = 0.0, res_qz_max = 1.0;
    std::vector<double> prof_z;              // Slice profile (laterally averaged density)
    std::vector<double> prof_rho;
    // Cuts
    std::vector<double> hcut, vcut, vcut_af;
    size_t hcut_rows = 0;
    size_t vcut_cols = 0;
    double link_qpar_min = 0.0, link_qpar_max = 1.0;
    double link_qz_min = 0.0, link_qz_max = 1.0;
    // Detector image
    gisaxs::DetectorImage det_img;
    std::vector<double> det_display;         // [row 0 = top], display values
    uint64_t det_hash = 0;

    md_allocator_i* alloc = nullptr;

    ScatteringComponent() { viamd::event_system_register_handler(*this); }

    // ---------------------------------------------------------------------------------------------
    // Events
    // ---------------------------------------------------------------------------------------------

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize:
                app_state = (ApplicationState*)e.payload;
                alloc = md_get_heap_allocator();
                break;
            case viamd::EventType_ViamdShutdown:
                cancel_all();
                destroy_ctx();
                clear_result();
                break;
            case viamd::EventType_ViamdFrameTick:
                update();
                draw_window();
                break;
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox(VIAMD_SCATTERING_NEUTRON ? "Scattering (GISAXS / GISANS)" : "Scattering (GISAXS)", &show_window);
                break;
            case viamd::EventType_ViamdSystemFree:
                cancel_all();
                destroy_ctx();
                clear_result();
                status[0] = '\0';
                break;
            case viamd::EventType_ViamdSystemStateChanged:
            case viamd::EventType_ViamdTrajectoryInit:
                if (ctx || compute_state.load() == ComputeState_Running) {
                    structure_stale = true;
                }
                break;
            case viamd::EventType_ViamdSerialize:
                serialize(*(viamd::serialization_state_t*)e.payload);
                break;
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t& state = *(viamd::deserialization_state_t*)e.payload;
                // "GISAXS" is the section name before the component was generalized
                const str_t header = viamd::section_header(state);
                if (str_eq(header, STR_LIT("Scattering")) || str_eq(header, STR_LIT("GISAXS"))) {
                    deserialize(state);
                }
                break;
            }
            default:
                break;
            }
        }
    }

    // ---------------------------------------------------------------------------------------------
    // Derived quantities
    // ---------------------------------------------------------------------------------------------

#if VIAMD_SCATTERING_NEUTRON
    bool is_neutron() const { return radiation == Radiation_Neutron; }
#else
    static constexpr bool is_neutron() { return false; }
#endif

    double wavelength() const {
#if VIAMD_SCATTERING_NEUTRON
        if (is_neutron()) return (double)MAX(neu_wavelength, 0.1f);
#endif
        return gisaxs::wavelength_from_energy(MAX(energy_kev, 1.0e-3f));
    }

    double substrate_z() const {
        return sub_auto_z ? particle_z_min + sub_z_offset : sub_z;
    }

    gisaxs::Setup xray_setup() const {
        gisaxs::Setup s;
        s.wavelength = wavelength();
        s.alpha_i = alpha_i_deg * DEG_TO_RAD;
        s.dwba = dwba;
        s.graded = film_graded;
        s.z_substrate = substrate_z();
        s.rho_ambient = ambient_density;
        s.rho_substrate = sub_density;
        s.beta_substrate = sub_beta;
        s.roughness = sub_roughness;
        s.rho_material = material_density;
        s.beta_material = material_beta;
        return s;
    }

#if VIAMD_SCATTERING_NEUTRON
    // Neutron background SLD (10^-6 Å^-2)
    double neutron_bg_sld() const {
        return neu_bg_preset == NBG_WATER ? neutron::water_sld(neu_d2o) : neu_bg_sld;
    }

    // D fraction of the exchangeable hydrogen reservoir
    double neutron_exchange_d() const {
        return (neu_exchange_follow && neu_bg_preset == NBG_WATER) ? (double)neu_d2o : (double)neu_exchange_d;
    }

    // Weights are excess scattering lengths in fm: the contrast is in the weights, the profile is the excess SLD
    md_gisaxs_model_t neutron_model() const {
        md_gisaxs_model_t m = {};
        m.wavelength = wavelength();
        m.alpha_i = alpha_i_deg * DEG_TO_RAD;
        m.dwba = dwba;
        m.graded = film_graded;
        m.z_substrate = substrate_z();
        m.sld_ambient = neutron_bg_sld() * 1.0e-6;
        m.sld_substrate = neu_sub_sld * 1.0e-6;
        m.sld_substrate_abs = neu_sub_abs * 1.0e-6;
        m.substrate_roughness = sub_roughness;
        m.profile_sld_scale = neutron::FM_TO_ANGSTROM;
        m.profile_abs_scale = 0.0;
        m.intensity_scale = neutron::FM_TO_ANGSTROM * neutron::FM_TO_ANGSTROM;
        return m;
    }

    // Flat incoherent background (sr^-1) added to the result
    double incoherent_level() const {
        return (ctx_radiation == Radiation_Neutron && neu_incoherent) ? inc_level : 0.0;
    }
#endif

    md_gisaxs_model_t build_model() const {
#if VIAMD_SCATTERING_NEUTRON
        if (is_neutron()) return neutron_model();
#endif
        return gisaxs::make_model(xray_setup());
    }

    // Critical angle (deg) of the substrate for the current settings, 0 without total reflection
    double critical_angle_deg() const {
        const md_gisaxs_model_t m = build_model();
        return gisaxs::critical_angle(m.sld_substrate - m.sld_ambient, m.wavelength) * RAD_TO_DEG;
    }

    uint64_t hash_structure_inputs() const {
        uint64_t h = md_hash64(filter, strnlen(filter, sizeof(filter)), 0x6153u);
#if VIAMD_SCATTERING_NEUTRON
        h = md_hash64(&radiation, sizeof(radiation), h);
        if (is_neutron()) {
            const double v[] = { (double)neu_weight_mode, neu_b_protiated, neu_n_h, neu_n_ex, neu_volume, neu_mass_density,
                                 neu_deuteration, (double)neu_exchange, neutron_exchange_d(), neutron_bg_sld() };
            h = md_hash64(v, sizeof(v), h);
        } else
#endif
        {
            h = md_hash64(&electron_mode, sizeof(electron_mode), h);
            h = md_hash64(&electrons_per_particle, sizeof(electrons_per_particle), h);
            h = md_hash64(&electrons_per_dalton, sizeof(electrons_per_dalton), h);
        }
        h = md_hash64(&sigma_mode, sizeof(sigma_mode), h);
        h = md_hash64(&sigma_constant, sizeof(sigma_constant), h);
        h = md_hash64(&q_par_max_nm, sizeof(q_par_max_nm), h);
        h = md_hash64(&q_z_max_nm, sizeof(q_z_max_nm), h);
        h = md_hash64(&oversampling, sizeof(oversampling), h);
        h = md_hash64(&ring_rel_width, sizeof(ring_rel_width), h);
        h = md_hash64(&max_slices, sizeof(max_slices), h);
        h = md_hash64(&density_model, sizeof(density_model), h);
        if (density_model == DensityModel_Fibril) {
            h = md_hash64(fib_center_name, strnlen(fib_center_name, sizeof(fib_center_name)), h);
            h = md_hash64(&fib_stride, sizeof(fib_stride), h);
            h = md_hash64(&fib_template, sizeof(fib_template), h);
            h = md_hash64(&fib_bead_sigma, sizeof(fib_bead_sigma), h);
            h = md_hash64(&fib_disk_radius, sizeof(fib_disk_radius), h);
            h = md_hash64(&fib_disk_spacing, sizeof(fib_disk_spacing), h);
            h = md_hash64(&fib_ds, sizeof(fib_ds), h);
            h = md_hash64(fib_template_path, strnlen(fib_template_path, sizeof(fib_template_path)), h);
        }
        return h;
    }

    uint64_t hash_model(const md_gisaxs_model_t& m) const {
        const double vals[] = { m.wavelength, m.alpha_i, m.dwba ? 1.0 : 0.0, m.graded ? 1.0 : 0.0, m.z_substrate, m.sld_ambient,
                                m.sld_substrate, m.sld_substrate_abs, m.substrate_roughness, m.profile_sld_scale,
                                m.profile_abs_scale, m.intensity_scale };
        uint64_t h = md_hash64(vals, sizeof(vals), 0x4d6fu);
#if VIAMD_SCATTERING_NEUTRON
        const double inc = incoherent_level();
        h = md_hash64(&inc, sizeof(inc), h);
#endif
        h = md_hash64(&num_qz, sizeof(num_qz), h);
        // Not the ctx pointer: a recompute frees and reallocates ctx, which usually lands at the same address, and the
        // result of the new structure (e.g. another trajectory frame) would then never be evaluated
        h = md_hash64(&ctx_gen, sizeof(ctx_gen), h);
        h = md_hash64(&structure_hash, sizeof(structure_hash), h);
        return h;
    }

    // Geometry of the shown result (nm^-1 and degrees, for the plots)
    gisaxs::Beam eval_beam() const {
        gisaxs::Beam b;
        b.wavelength = eval_model.wavelength > 0.0 ? eval_model.wavelength : wavelength();
        b.alpha_i = eval_model.alpha_i;
        return b;
    }
    double horizon_qz_nm() const { return gisaxs::horizon_qz(eval_beam()) * 10.0; }
    double yoneda_qz_nm() const { return gisaxs::yoneda_qz(eval_beam(), eval_model.sld_substrate - eval_model.sld_ambient) * 10.0; }
    double film_sld_rel() const { return gisaxs::film_sld(prof_rho.data(), prof_rho.size(), eval_model.profile_sld_scale); }
    double film_yoneda_qz_nm() const { return gisaxs::yoneda_qz(eval_beam(), film_sld_rel()) * 10.0; }
    double qz_nm_from_alpha_f(double af_deg) const { return gisaxs::qz_from_alpha_f(eval_beam(), af_deg * DEG_TO_RAD) * 10.0; }
    double alpha_f_deg(double qz) const { return gisaxs::alpha_f_from_qz(eval_beam(), qz * 0.1) * RAD_TO_DEG; }
    bool   has_result() const { return !smooth.empty() && !map_display.empty(); }

    // ---------------------------------------------------------------------------------------------
    // Beam presets
    // ---------------------------------------------------------------------------------------------

    struct PresetValues {
        float  energy_kev;
        float  alpha_i_deg;
        float  q_par_max_nm;
        float  q_z_max_nm;
        double sub_density;
        double sub_beta;
    };

    static PresetValues preset_values(const BeamPreset& p) {
        const double lambda_a = p.wavelength_nm * 10.0;
        const double k0 = 2.0 * PI / p.wavelength_nm;   // nm^-1
        // Largest in-plane angle: the detector edge farthest from the direct beam, so the map covers every pixel
        const double edge_px = MAX(p.beam_x_px, p.det_px_h - p.beam_x_px);
        const double two_theta_max = atan(edge_px * p.pixel_mm / p.sdd_mm);
        const double ai = p.alpha_i_deg * DEG_TO_RAD;
        const double af = p.alpha_f_max_deg * DEG_TO_RAD;
        PresetValues v;
        v.energy_kev   = (float)gisaxs::energy_from_wavelength(lambda_a);
        v.alpha_i_deg  = (float)p.alpha_i_deg;
        v.q_par_max_nm = (float)(k0 * sin(two_theta_max));
        v.q_z_max_nm   = (float)(k0 * (sin(af) + sin(ai)));
        v.sub_density  = gisaxs::electron_density_from_delta(p.sub_delta, lambda_a);
        v.sub_beta     = p.sub_beta;
        return v;
    }

    bool preset_matches(const BeamPreset& p) const {
        const PresetValues v = preset_values(p);
        auto eq = [](double a, double b) { return fabs(a - b) <= 1.0e-4 * MAX(fabs(a), fabs(b)); };
        return eq(energy_kev, v.energy_kev) && eq(alpha_i_deg, v.alpha_i_deg) &&
               eq(q_par_max_nm, v.q_par_max_nm) && eq(q_z_max_nm, v.q_z_max_nm) &&
               dwba && eq(sub_density, v.sub_density) && eq(sub_beta, v.sub_beta);
    }

    void apply_preset(int idx) {
        if (idx < 0 || idx >= (int)ARRAY_SIZE(beam_presets)) return;
        const BeamPreset& p = beam_presets[idx];
        const PresetValues v = preset_values(p);
#if VIAMD_SCATTERING_NEUTRON
        radiation    = Radiation_Xray;   // Presets carry X-ray optical constants (delta, beta)
#endif
        energy_kev   = v.energy_kev;
        alpha_i_deg  = v.alpha_i_deg;
        q_par_max_nm = v.q_par_max_nm;
        q_z_max_nm   = v.q_z_max_nm;
        dwba         = true;
        sub_preset   = SUB_CUSTOM;
        sub_density  = v.sub_density;
        sub_beta     = v.sub_beta;
        cut_mode        = CutMode_AlphaF;
        cut_alpha_f_deg = (float)p.cut_alpha_f_deg;
        cut_qpar        = p.cut_qpar_nm;
        cut_init        = true;
        vcut_vs_alpha_f = true;     // As in the papers (Ohm 2018 Fig. 4c, Brett 2019 Fig. 4b)
        film_graded     = true;
        material_beta   = p.material_beta;
        det.sdd_mm      = p.sdd_mm;
        det.pixel_mm    = p.pixel_mm;
        det.npx_h       = p.det_px_h;
        det.npx_v       = p.det_px_v;
        det.beam_x_px   = p.beam_x_px;
        det.beam_y_px   = p.beam_y_px;
        det.module_w    = p.module_w;
        det.module_h    = p.module_h;
        det.gap_w       = p.gap_w;
        det.gap_h       = p.gap_h;
        det.gaps        = p.module_w > 0;
        det.bs_direct_px   = p.bs_direct_px;
        det.bs_specular_px = p.bs_specular_px;
        // Resolution: one detector pixel (the angular beam divergence at P03 is of the same order or smaller)
        const float dq = (float)(2.0 * PI / p.wavelength_nm * p.pixel_mm / p.sdd_mm);
        res_fwhm_qpar = dq;
        res_fwhm_qz = dq;
        res_dirty = true;
        beam_preset = idx;
    }

    // ---------------------------------------------------------------------------------------------
    // Lifetime helpers
    // ---------------------------------------------------------------------------------------------

    void cancel_all() {
        compute_cancel = true;
        task_system::task_interrupt_and_wait_for(task_slices);
        task_system::task_interrupt_and_wait_for(task_rings);
        task_system::task_wait_for(task_finish);
        // A partial evaluation still completes (and sets eval_ready), it is discarded below
        task_system::task_interrupt(task_eval_range);
        task_system::task_wait_for(task_eval);
        task_slices = task_rings = task_finish = task_eval_range = task_eval = task_system::INVALID_ID;
        if (compute_state.load() == ComputeState_Running) {
            compute_state = ComputeState_Idle;
        }
        eval_ready = false;
        compute_cancel = false;
    }

    void free_scratch() {
        for (size_t i = 0; i < md_array_size(scratch); ++i) {
            if (scratch[i]) md_fft_free(scratch[i]);
        }
        md_array_free(scratch, alloc);
        scratch = nullptr;
        scratch_bytes = 0;
    }

    void destroy_ctx() {
        free_scratch();
        if (ctx) {
            md_gisaxs_destroy(ctx);
            ctx = nullptr;
        }
        compute_state = ComputeState_Idle;
    }

    void clear_result() {
        raw = gisaxs::QMap();
        smooth = gisaxs::QMap();
        map_display.clear();
        ring_q_nm.clear();
        ring_edge_nm.clear();
        qz_nm.clear();
        prof_z.clear();
        prof_rho.clear();
        hcut.clear();
        vcut.clear();
        vcut_af.clear();
        det_img = gisaxs::DetectorImage();
        det_display.clear();
        det_hash = 0;
        eval_hash_shown = 0;
        eval_hash_pending = 0;
        res_dirty = false;
    }

    // ---------------------------------------------------------------------------------------------
    // Structure stage
    // ---------------------------------------------------------------------------------------------

    // Detects slices (groups of fib_stride consecutive beads containing one center bead) and chains of consecutive
    // slices, and sweeps the cross-section template along them.
    // x, y, z, w are the gathered (selected) particles, atom_idx their atom indices (ascending). w must be positive, it
    // is used for the slice centroids and orientations as well as for the distribution of the scattering weight.
    bool build_fibrils(FibrilOutput* out, const std::vector<uint32_t>& atom_idx, const float* x, const float* y, const float* z, const float* w, const md_unitcell_t& cell) {
        const md_system_t& sys = app_state->mold.sys;
        const int S = fib_stride;
        const size_t n = atom_idx.size();
        if (S < 2 || S > 64) {
            snprintf(status, sizeof(status), "Fibril model: beads per slice must be in [2, 64]");
            return false;
        }
        const str_t center = str_from_cstr(fib_center_name);
        const double box[3] = {cell.x, cell.y, cell.z};
        auto min_img = [&](double d, int k) { return box[k] > 0.0 ? d - box[k] * floor(d / box[k] + 0.5) : d; };

        // Selected center beads (positions in the gathered arrays)
        std::vector<size_t> centers;
        for (size_t i = 0; i < n; ++i) {
            if (str_eq(md_atom_name(&sys.atom, atom_idx[i]), center)) centers.push_back(i);
        }
        if (centers.empty()) {
            snprintf(status, sizeof(status), "Fibril model: no '%s' beads in the selection", fib_center_name);
            return false;
        }

        // A group [g, g + S) is valid if it lies within the gathered set, is contiguous in atom index and contains
        // exactly one center bead.
        auto group_valid = [&](long g) -> bool {
            if (g < 0 || g + S > (long)n) return false;
            if (atom_idx[g + S - 1] - atom_idx[g] != (uint32_t)(S - 1)) return false;
            int nc = 0;
            for (int b = 0; b < S; ++b) nc += str_eq(md_atom_name(&sys.atom, atom_idx[g + b]), center) ? 1 : 0;
            return nc == 1;
        };

        // Position of the center bead within its slice: choose the offset that keeps the beads closest to the center
        int best_off = 0;
        double best_cost = DBL_MAX;
        const size_t sample = MIN(centers.size(), (size_t)2000);
        for (int off = 0; off < S; ++off) {
            double cost = 0.0;
            size_t valid = 0;
            for (size_t k = 0; k < sample; ++k) {
                const size_t ci = centers[k * centers.size() / sample];
                const long g = (long)ci - off;
                if (!group_valid(g)) continue;
                for (int b = 0; b < S; ++b) {
                    const double dx = min_img(x[g + b] - x[ci], 0), dy = min_img(y[g + b] - y[ci], 1), dz = min_img(z[g + b] - z[ci], 2);
                    cost += sqrt(dx * dx + dy * dy + dz * dz);
                }
                valid += 1;
            }
            if (valid == 0) continue;
            cost /= (double)valid;
            // Prefer offsets that are valid for (almost) all slices
            cost *= (double)sample / (double)valid;
            if (cost < best_cost) { best_cost = cost; best_off = off; }
        }

        // Slices
        std::vector<long> slice_start;
        for (size_t ci : centers) {
            const long g = (long)ci - best_off;
            if (group_valid(g)) slice_start.push_back(g);
        }
        if (slice_start.size() < 2) {
            snprintf(status, sizeof(status), "Fibril model: could not identify slices of %i beads around '%s'", S, fib_center_name);
            return false;
        }

        // Chains: consecutive slices in atom order whose center beads are bonded (or close, without bonds)
        std::vector<double> gaps;
        for (size_t k = 0; k + 1 < slice_start.size(); ++k) {
            const size_t a = slice_start[k] + best_off, b = slice_start[k + 1] + best_off;
            const double dx = min_img(x[b] - x[a], 0), dy = min_img(y[b] - y[a], 1), dz = min_img(z[b] - z[a], 2);
            gaps.push_back(sqrt(dx * dx + dy * dy + dz * dz));
        }
        std::vector<double> sorted = gaps;
        std::sort(sorted.begin(), sorted.end());
        const double median_gap = sorted[sorted.size() / 2];

        size_t bonded = 0, tested = 0;
        const bool has_bonds = sys.bond.count > 0 && sys.bond.conn.offset != nullptr;
        if (has_bonds) {
            for (size_t k = 0; k + 1 < slice_start.size() && tested < 2000; ++k) {
                const uint32_t a = atom_idx[slice_start[k] + best_off], b = atom_idx[slice_start[k + 1] + best_off];
                tested += 1;
                if (md_bond_find(&sys.bond, a, b) != (md_bond_idx_t)-1) bonded += 1;
            }
        }
        const bool use_bonds = has_bonds && bonded * 2 > tested;

        std::vector<uint32_t> chain_offset;
        chain_offset.push_back(0);
        for (size_t k = 0; k + 1 < slice_start.size(); ++k) {
            bool connected = slice_start[k + 1] == slice_start[k] + S &&
                             atom_idx[slice_start[k + 1]] == atom_idx[slice_start[k]] + (uint32_t)S;
            if (connected) {
                if (use_bonds) {
                    const uint32_t a = atom_idx[slice_start[k] + best_off], b = atom_idx[slice_start[k + 1] + best_off];
                    connected = md_bond_find(&sys.bond, a, b) != (md_bond_idx_t)-1;
                } else {
                    connected = gaps[k] < 1.6 * median_gap;
                }
            }
            if (!connected) chain_offset.push_back((uint32_t)(k + 1));
        }
        chain_offset.push_back((uint32_t)slice_start.size());

        // Bead data in slice order
        const size_t num_slices = slice_start.size();
        std::vector<float> bx(num_slices * S), by(num_slices * S), bz(num_slices * S), be(num_slices * S);
        for (size_t k = 0; k < num_slices; ++k) {
            for (int b = 0; b < S; ++b) {
                const size_t src = slice_start[k] + b;
                bx[k * S + b] = x[src];
                by[k * S + b] = y[src];
                bz[k * S + b] = z[src];
                be[k * S + b] = w[src];
            }
        }

        FibrilInput in;
        in.num_slices = num_slices;
        in.stride = S;
        in.x = bx.data(); in.y = by.data(); in.z = bz.data(); in.e = be.data();
        in.chain_offset = chain_offset.data();
        in.num_chains = chain_offset.size() - 1;
        in.box[0] = box[0]; in.box[1] = box[1]; in.box[2] = box[2];
        in.ref_bead = best_off == 0 ? 1 : 0;

        std::vector<float> mu, mv, mf;
        if (!fibril_mean_layout(mu, mv, mf, in)) {
            snprintf(status, sizeof(status), "Fibril model: no chains with at least two slices");
            return false;
        }

        FibrilTemplate tmpl;
        char err[256] = "";
        switch (fib_template) {
        case FibrilTemplate_Disk:
            tmpl = fibril_template_disk(mu, mv, fib_disk_radius, fib_disk_spacing);
            break;
        case FibrilTemplate_File:
            if (!fibril_template_load(&tmpl, fib_template_path, err, sizeof(err))) {
                snprintf(status, sizeof(status), "Fibril template: %s", err);
                return false;
            }
            if ((int)tmpl.anchor_u.size() != S) {
                // No (or incompatible) anchors in the file, orient with the mean bead layout
                if (!tmpl.anchor_u.empty()) {
                    MD_LOG_INFO("Scattering fibrils: template has %zu anchors, expected %i, orienting with the mean MD layout instead", tmpl.anchor_u.size(), S);
                }
                tmpl.anchor_u = mu;
                tmpl.anchor_v = mv;
            } else {
                double rms = 0.0;
                const bool mirrored = fibril_template_match_handedness(tmpl, mu, mv, &rms);
                MD_LOG_INFO("Scattering fibrils: template anchors match the mean MD slice layout to %.2f Å RMS%s", rms, mirrored ? " (template mirrored)" : "");
            }
            break;
        default:
            tmpl = fibril_template_beads(mu, mv, mf, fib_bead_sigma);
            break;
        }

        if (!fibril_sweep(out, in, tmpl, fib_ds, err, sizeof(err))) {
            snprintf(status, sizeof(status), "Fibril model: %s", err);
            return false;
        }
        const FibrilStats& st = out->stats;
        snprintf(fib_info, sizeof(fib_info), "%zu slices in %zu chains (%s), spacing %.1f Å, bead radius %.1f Å, %zu template Gaussians, %zu points",
            num_slices, in.num_chains, use_bonds ? "bonds" : "distance", st.mean_spacing, st.mean_radius, tmpl.gauss.size(), st.num_points);
        MD_LOG_INFO("Scattering fibrils: %s (weight %.6g -> %.6g, RMS axial bead offset %.2f Å)", fib_info, st.electrons_in, st.electrons_out, st.mean_axial_offset);
        return true;
    }

    // Gathered particles of the selection. w: scattering weight (X-rays: electrons, neutrons: excess scattering length
    // in fm), gw: positive weight for the fibril geometry (X-rays: electrons, neutrons: volume), s: Gaussian width.
    struct Particles {
        std::vector<float> x, y, z, w, gw, s;
        std::vector<uint32_t> atom_idx;
        size_t size() const { return x.size(); }
    };

    bool gather_particles(Particles* p, const md_bitfield_t& mask) {
        const md_system_t& sys = app_state->mold.sys;
        const md_system_state_t& state = app_state->mold.state;
        const size_t num_atoms = sys.atom.count;
        const bool has_types = sys.atom.type_idx != nullptr;
        const size_t count = md_bitfield_popcount(&mask);
        p->x.reserve(count); p->y.reserve(count); p->z.reserve(count);
        p->w.reserve(count); p->gw.reserve(count); p->s.reserve(count);
        p->atom_idx.reserve(count);

#if VIAMD_SCATTERING_NEUTRON
        const bool neutron_mode = is_neutron();
        neutron::HydrogenModel hmodel;
        hmodel.deuteration = neu_deuteration;
        hmodel.exchange = neu_exchange;
        hmodel.exchange_d = neutron_exchange_d();
        const double bg_sld = neutron_bg_sld();
        neutron::Particle neu_particle;
        neu_particle.b_protiated = neu_b_protiated;
        neu_particle.n_h = neu_n_h;
        neu_particle.n_ex = neu_n_ex;
        neu_particle.volume = neu_volume;
        const neutron::Scattering neu_const = neutron::particle(neu_particle, hmodel);
        const bool per_element = neutron_mode && neu_weight_mode == NeutronWeight_Element;
        const bool has_bonds = sys.bond.count > 0 && sys.bond.conn.offset != nullptr;
        if (per_element && neu_exchange && !has_bonds) {
            MD_LOG_INFO("Scattering: the system has no bonds, no hydrogen is treated as exchangeable");
        }
        neu_sum_b = neu_sum_v = neu_sum_inc = 0.0;
        neu_num_unknown = neu_num_exchangeable = 0;
#endif

        size_t num_zero_weight = 0;
        md_bitfield_iter_t it = md_bitfield_iter_create(&mask);
        while (md_bitfield_iter_next(&it)) {
            const size_t idx = md_bitfield_iter_idx(&it);
            if (idx >= num_atoms) continue;
            float w = 0.0f, gw = 0.0f;
#if VIAMD_SCATTERING_NEUTRON
            if (neutron_mode) {
                neutron::Scattering sc = neu_const;
                if (per_element) {
                    const int zn = (int)md_atom_atomic_number(&sys.atom, idx);
                    const double mass = has_types ? md_atom_mass(&sys.atom, idx) : 0.0;
                    // Hydrogen bound to N, O or S is exchangeable
                    bool exch = false;
                    if (zn == 1 && neu_exchange && has_bonds) {
                        for (md_bond_iter_t bit = md_bond_iter(&sys.bond, idx); md_bond_iter_has_next(&bit); md_bond_iter_next(&bit)) {
                            const int zo = (int)md_atom_atomic_number(&sys.atom, md_bond_iter_atom_index(&bit));
                            if (zo == 7 || zo == 8 || zo == 16) { exch = true; break; }
                        }
                    }
                    if (!neutron::atom(&sc, zn, mass, exch, hmodel, neu_mass_density)) neu_num_unknown += 1;
                    if (exch) neu_num_exchangeable += 1;
                }
                w  = (float)neutron::excess(sc, bg_sld);
                gw = (float)sc.volume;
                neu_sum_b   += sc.b;
                neu_sum_v   += sc.volume;
                neu_sum_inc += sc.sigma_inc;
            } else
#endif
            {
                if (electron_mode == ElectronMode_AtomicNumber) {
                    w = (float)md_atom_atomic_number(&sys.atom, idx);
                } else if (electron_mode == ElectronMode_Mass) {
                    w = has_types ? md_atom_mass(&sys.atom, idx) * electrons_per_dalton : 0.0f;
                } else {
                    w = electrons_per_particle;
                }
                gw = w;
            }
            if (w == 0.0f) num_zero_weight += 1;
            p->x.push_back(state.xyz[idx].x);
            p->y.push_back(state.xyz[idx].y);
            p->z.push_back(state.xyz[idx].z);
            p->w.push_back(w);
            p->gw.push_back(gw);
            p->s.push_back(sigma_mode == SigmaMode_Radius ? md_atom_radius(&sys.atom, idx) / sqrtf(5.0f) : sigma_constant);
            p->atom_idx.push_back((uint32_t)idx);
        }
        const size_t n = p->size();

#if VIAMD_SCATTERING_NEUTRON
        if (neutron_mode) {
            neu_bg_sld_used = bg_sld;
            if (per_element && neu_num_unknown == n) {
                snprintf(status, sizeof(status), "No particle has an element with a tabulated scattering length, use a constant scattering length per particle");
                return false;
            }
            if (neu_sum_v <= 0.0) {
                snprintf(status, sizeof(status), per_element ? "Particle volumes are zero (masses missing?), use a constant scattering length per particle"
                                                             : "The particle volume must be positive");
                return false;
            }
            if (num_zero_weight == n) {
                snprintf(status, sizeof(status), "The selection is contrast matched to the background (zero excess scattering length)");
                return false;
            }
            if (neu_num_unknown) {
                MD_LOG_INFO("Scattering: %zu particles have an element without a tabulated neutron scattering length (b = 0)", neu_num_unknown);
            }
            return true;
        }
#endif
        if (num_zero_weight == n) {
            snprintf(status, sizeof(status), "All particles have zero electrons (atomic numbers or masses missing?), use a constant electron count");
            return false;
        }
        return true;
    }

    bool start_compute() {
        cancel_all();
        destroy_ctx();
        structure_stale = false;

        const md_system_t& sys = app_state->mold.sys;
        const md_system_state_t& state = app_state->mold.state;
        if (sys.atom.count == 0 || !state.xyz) {
            snprintf(status, sizeof(status), "No system loaded");
            return false;
        }

        const md_unitcell_t& cell = state.unitcell;
        const bool ortho = (cell.flags & MD_UNITCELL_ORTHO) ||
            ((cell.flags & MD_UNITCELL_TRICLINIC) && cell.xy == 0.0 && cell.xz == 0.0 && cell.yz == 0.0);
        if (!ortho || cell.x <= 0.0 || cell.y <= 0.0) {
            snprintf(status, sizeof(status), "An orthorhombic periodic box is required (periodic in XY)");
            return false;
        }

        // Selection
        md_bitfield_t mask = {0};
        md_bitfield_init(&mask, alloc);
        defer { md_bitfield_free(&mask); };
        filter_valid = md_filter(&mask, str_from_cstr(filter), &sys, &state, app_state->script.ir, NULL, filter_err, sizeof(filter_err));
        if (!filter_valid) {
            snprintf(status, sizeof(status), "Invalid selection: %s", filter_err);
            return false;
        }
        if (md_bitfield_popcount(&mask) == 0) {
            snprintf(status, sizeof(status), "Selection is empty");
            return false;
        }

        Particles p;
        if (!gather_particles(&p, mask)) return false;
        const size_t n = p.size();

        // Self scattering of the beads (diagnostic)
        bead_self.clear();
        bead_self_area = cell.x * cell.y;
        for (size_t i = 0; i < n; ++i) {
            size_t g = 0;
            for (; g < bead_self.size(); ++g) if (bead_self[g].first == p.s[i]) break;
            if (g == bead_self.size()) {
                if (bead_self.size() >= 16) g = bead_self.size() - 1;   // Enough for a diagnostic
                else bead_self.push_back({p.s[i], 0.0});
            }
            bead_self[g].second += (double)p.w[i] * p.w[i];
        }

        // Fibril model: replace the beads by a continuous swept cross-section
        FibrilOutput fib;
        fib_info[0] = '\0';
        if (density_model == DensityModel_Fibril) {
            if (!build_fibrils(&fib, p.atom_idx, p.x.data(), p.y.data(), p.z.data(), p.gw.data(), cell)) {
                return false;
            }
#if VIAMD_SCATTERING_NEUTRON
            if (is_neutron()) {
                // The geometry was swept with the volumes; carry the excess scattering length with the mean excess SLD
                // of the selection, consistent with the template being a homogeneous cross-section.
                double sum_w = 0.0;
                for (size_t i = 0; i < n; ++i) sum_w += p.w[i];
                const float ratio = (float)(sum_w / neu_sum_v);
                for (float& fw : fib.w) fw *= ratio;
            }
#endif
        }

        md_gisaxs_input_t input = {};
        if (density_model == DensityModel_Fibril) {
            input.count = fib.x.size();
            input.x = fib.x.data();
            input.y = fib.y.data();
            input.z = fib.z.data();
            input.weight = fib.w.data();
            input.sigma = fib.s.data();
        } else {
            input.count = n;
            input.x = p.x.data();
            input.y = p.y.data();
            input.z = p.z.data();
            input.weight = p.w.data();
            input.sigma = (sigma_mode == SigmaMode_Constant) ? nullptr : p.s.data();
            input.sigma_uniform = sigma_constant;
        }
        input.box_x = cell.x;
        input.box_y = cell.y;

        md_gisaxs_params_t params = {};
        params.q_par_max = q_par_max_nm * 0.1;
        params.q_z_max = q_z_max_nm * 0.1;
        params.oversampling = oversampling;
        params.ring_rel_width = ring_rel_width;
        params.max_slices = (size_t)MAX(max_slices, 8);

        ctx = md_gisaxs_create(&input, &params, alloc);
        if (!ctx) {
            snprintf(status, sizeof(status), "Failed to initialize the scattering computation (see log)");
            return false;
        }
        ctx_gen += 1;
#if VIAMD_SCATTERING_NEUTRON
        ctx_radiation = radiation;
        // Flat incoherent background: sum sigma_inc / (4 pi A), per unit area and solid angle like the coherent part
        inc_level = is_neutron() ? neu_sum_inc * neutron::BARN_TO_A2 / (4.0 * PI * cell.x * cell.y) : 0.0;
#endif
        ctx_q_z_max = params.q_z_max;
        particle_z_min = md_gisaxs_particle_z_min(ctx);
        structure_hash = hash_structure_inputs();

        md_gisaxs_info_t info;
        md_gisaxs_get_info(ctx, &info);
        MD_LOG_INFO("Scattering: %zu particles, grid %i x %i (dx %.2f Å), %zu slices (dz %.2f Å), %zu rings, %zu classes, %.1f MB spectra, %.1f MB matrices",
            info.num_particles, info.nx, info.ny, info.dx, info.num_slices, info.dz, info.num_rings, info.num_classes,
            info.spectra_bytes / (1024.0 * 1024.0), info.matrix_bytes / (1024.0 * 1024.0));

        // Per thread scratch, allocated lazily by each worker
        const size_t num_threads = MAX(task_system::pool_num_threads(), 1);
        md_array_resize(scratch, num_threads, alloc);
        for (size_t i = 0; i < num_threads; ++i) scratch[i] = nullptr;
        scratch_bytes = md_gisaxs_slice_scratch_bytes(ctx);

        compute_state = ComputeState_Running;
        compute_cancel = false;
        compute_time_start = ImGui::GetTime();

        task_slices = task_system::create_pool_task(STR_LIT("Scattering slices"), (uint32_t)info.num_slices, [this](uint32_t beg, uint32_t end, uint32_t thread_num) {
            if (compute_cancel) return;
            if (thread_num >= md_array_size(scratch)) return;   // Should not happen
            if (!scratch[thread_num]) {
                scratch[thread_num] = md_fft_alloc(scratch_bytes);
                if (!scratch[thread_num]) { compute_cancel = true; return; }
            }
            md_gisaxs_compute_slices(ctx, beg, end, scratch[thread_num]);
        });

        task_rings = task_system::create_pool_task(STR_LIT("Scattering rings"), (uint32_t)info.num_rings, [this](uint32_t beg, uint32_t end, uint32_t) {
            if (compute_cancel) return;
            md_gisaxs_compute_rings(ctx, beg, end);
        });

        task_finish = task_system::create_pool_task(STR_LIT("##Scattering finalize"), [this]() {
            md_gisaxs_release_spectra(ctx);
            compute_state = compute_cancel ? ComputeState_Failed : ComputeState_Done;
        });

        task_system::set_task_dependency(task_rings, task_slices);
        task_system::set_task_dependency(task_finish, task_rings);
        task_system::enqueue_task(task_slices);

        snprintf(status, sizeof(status), "Computing: %zu particles, %i x %i grid, %zu slices, %zu rings",
            info.num_particles, info.nx, info.ny, info.num_slices, info.num_rings);
        return true;
    }

    // ---------------------------------------------------------------------------------------------
    // Model stage
    // ---------------------------------------------------------------------------------------------

    void start_eval(const md_gisaxs_model_t& model, uint64_t hash) {
        const size_t R = md_gisaxs_num_rings(ctx);
        const size_t Nq = (size_t)MAX(num_qz, 2);
        eval_qz.resize(Nq);
        eval_out.assign(Nq * R, 0.0f);
        for (size_t i = 0; i < Nq; ++i) {
            eval_qz[i] = ctx_q_z_max * (double)i / (double)(Nq - 1);
        }
        eval_model = model;
#if VIAMD_SCATTERING_NEUTRON
        eval_inc = incoherent_level();
#endif
        eval_hash_pending = hash;
        eval_ready = false;

        task_eval_range = task_system::create_pool_task(STR_LIT("Scattering evaluate"), (uint32_t)Nq, [this](uint32_t beg, uint32_t end, uint32_t) {
            md_gisaxs_evaluate_range(ctx, &eval_model, eval_qz.data(), beg, end, eval_out.data());
        }, 4);
        // Bookkeeping only, hidden from the async task overlay ("##")
        task_eval = task_system::create_pool_task(STR_LIT("##Scattering evaluate done"), [this]() {
            eval_ready = true;
        });
        task_system::set_task_dependency(task_eval, task_eval_range);
        task_system::enqueue_task(task_eval_range);
    }

    void accept_eval() {
        const size_t Nq = eval_qz.size();
        gisaxs::qmap_init_rings(&raw, ctx);
        const size_t R = raw.cols();
        raw.rows = Nq;
        raw.qz0 = eval_qz[0];
        raw.dqz = Nq > 1 ? eval_qz[1] - eval_qz[0] : 0.0;
        raw.I = eval_out;
#if VIAMD_SCATTERING_NEUTRON
        if (eval_inc > 0.0) {
            // Flat incoherent background, above the sample horizon only when there is a substrate
            const double horizon = eval_model.dwba ? gisaxs::horizon_qz(eval_beam()) : -DBL_MAX;
            for (size_t i = 0; i < Nq; ++i) {
                if (raw.qz(i) < horizon) continue;
                for (size_t c = 0; c < R; ++c) raw.I[i * R + c] += (float)eval_inc;
            }
        }
#endif

        md_gisaxs_info_t info;
        md_gisaxs_get_info(ctx, &info);
        lattice_dq = info.dq_ring;
        qpar_first = info.q_par_min;

        // Display axes (nm^-1). The rings are non-uniform in q_par: ring r has mean q_par ring_q[r] and spans
        // [edge[r], edge[r + 1]).
        ring_q_nm.resize(R);
        ring_edge_nm.resize(R + 1);
        qz_nm.resize(Nq);
        for (size_t i = 0; i < R; ++i) ring_q_nm[i] = raw.ring_q[i] * 10.0;
        for (size_t i = 0; i <= R; ++i) ring_edge_nm[i] = raw.ring_edge[i] * 10.0;
        for (size_t i = 0; i < Nq; ++i) qz_nm[i] = raw.qz(i) * 10.0;

        const double prev_bounds[4] = {res_qpar_min, res_qpar_max, res_qz_min, res_qz_max};
        res_qpar_min = ring_edge_nm[0];
        res_qpar_max = ring_edge_nm[R];
        res_qz_min = (raw.qz0 - 0.5 * raw.dqz) * 10.0;
        res_qz_max = (raw.qz(Nq - 1) + 0.5 * raw.dqz) * 10.0;
        // The map and the cut plots share linked axes. Reset the view when the q-range changed.
        if (prev_bounds[0] != res_qpar_min || prev_bounds[1] != res_qpar_max || prev_bounds[2] != res_qz_min || prev_bounds[3] != res_qz_max) {
            link_qpar_min = res_qpar_min;
            link_qpar_max = res_qpar_max;
            link_qz_min = res_qz_min;
            link_qz_max = res_qz_max;
        }

        const size_t S = md_gisaxs_num_slices(ctx);
        const double* sz = md_gisaxs_slice_z(ctx);
        const double* sp = md_gisaxs_slice_profile(ctx);
        prof_z.assign(sz, sz + S);
        prof_rho.assign(sp, sp + S);

        eval_hash_shown = eval_hash_pending;
        res_dirty = true;
        eval_ready = false;
        task_eval_range = task_eval = task_system::INVALID_ID;
    }

    // ---------------------------------------------------------------------------------------------
    // Per frame update
    // ---------------------------------------------------------------------------------------------

    void update() {
        const int cs = compute_state.load();
        if (cs == ComputeState_Done && scratch) {
            // Structure stage finished
            free_scratch();
            compute_time = ImGui::GetTime() - compute_time_start;
            task_slices = task_rings = task_finish = task_system::INVALID_ID;
            md_gisaxs_info_t info;
            md_gisaxs_get_info(ctx, &info);
            snprintf(status, sizeof(status), "%zu particles, %i x %i grid (%.1f Å), %zu slices (%.1f Å), %zu rings, %zu class(es) - %.2f s",
                info.num_particles, info.nx, info.ny, info.dx, info.num_slices, info.dz, info.num_rings, info.num_classes, compute_time);
        } else if (cs == ComputeState_Failed) {
            free_scratch();
            if (ctx) { md_gisaxs_destroy(ctx); ctx = nullptr; }
            compute_state = ComputeState_Idle;
            snprintf(status, sizeof(status), "Computation was cancelled or failed");
        }

        if (eval_ready.load()) {
            accept_eval();
        }

        bool ctx_matches = compute_state.load() == ComputeState_Done && ctx;
#if VIAMD_SCATTERING_NEUTRON
        // The particle weights depend on the radiation, a result computed for the other one is not re-evaluated
        ctx_matches = ctx_matches && ctx_radiation == radiation;
#endif
        if (ctx_matches && !task_system::task_is_running(task_eval) && !eval_ready.load()) {
            const md_gisaxs_model_t model = build_model();
            const uint64_t h = hash_model(model);
            if (h != eval_hash_shown) {
                start_eval(model, h);
            }
        }

        if (res_dirty && !raw.empty()) {
            refresh_display();
        }

        if (view_detector && has_result()) {
            const uint64_t h = hash_detector();
            if (h != det_hash) {
                compute_detector();
                det_hash = h;
            }
        }
    }

    // Instrument resolution, color range and the heat map (one column per ring)
    void refresh_display() {
        const double horizon = eval_model.dwba ? gisaxs::horizon_qz(eval_beam()) : -DBL_MAX;
        gisaxs::qmap_apply_resolution(&smooth, raw, res_fwhm_qpar * 0.1, res_fwhm_qz * 0.1, horizon);

        double vmax = -DBL_MAX, vmin = DBL_MAX;
        for (float v : smooth.I) {
            if (v > 0.0f) { vmax = MAX(vmax, (double)v); vmin = MIN(vmin, (double)v); }
        }
        if (vmax <= 0.0 || vmax == -DBL_MAX) { vmax = 1.0; vmin = 0.0; }
        if (log_scale) {
            res_max = log10(vmax);
            res_min = res_max - MAX(decades, 0.5f);
        } else {
            res_max = vmax;
            res_min = 0.0;
        }

        const size_t R = smooth.cols(), Q = smooth.rows;
        map_display.resize(R * Q);
        for (size_t c = 0; c < R; ++c) {
            double* dst = map_display.data() + c * Q;
            for (size_t row = 0; row < Q; ++row) {
                double v = smooth.at(Q - 1 - row, c);   // Row 0 of the heat map is drawn at the top (highest q_z)
                if (log_scale) v = v > 0.0 ? MAX(log10(v), res_min) : res_min;
                dst[row] = v;
            }
        }
        res_dirty = false;
        res_version += 1;
    }

    // ---------------------------------------------------------------------------------------------
    // Detector image
    // ---------------------------------------------------------------------------------------------

    uint64_t hash_detector() const {
        const double v[] = { eval_model.wavelength, eval_model.alpha_i, (double)log_scale, det.sdd_mm, det.pixel_mm,
            (double)det.npx_h, (double)det.npx_v, det.beam_x_px, det.beam_y_px, (double)det.binning, (double)det.gaps,
            (double)det.module_w, (double)det.module_h, (double)det.gap_w, (double)det.gap_h,
            (double)det.bs_direct_px, (double)det.bs_specular_px };
        return md_hash64(v, sizeof(v), res_version);
    }

    void compute_detector() {
        gisaxs::detector_render(&det_img, det, eval_beam(), smooth);
        const int W = det_img.width, H = det_img.height;
        det_display.resize((size_t)W * H);
        const double masked = res_min - 1.0e3 * (res_max - res_min + 1.0);   // below the color range
        for (int by = 0; by < H; ++by) {
            const float* src = det_img.I.data() + (size_t)by * W;
            double* dst = det_display.data() + (size_t)(H - 1 - by) * W;   // row 0 at the top
            for (int bx = 0; bx < W; ++bx) {
                const double v = src[bx];
                if (v < 0.0)        dst[bx] = masked;
                else if (log_scale) dst[bx] = v > 0.0 ? MAX(log10(v), res_min) : res_min;
                else                dst[bx] = v;
            }
        }
    }

    // ---------------------------------------------------------------------------------------------
    // Cuts
    // ---------------------------------------------------------------------------------------------

    // Horizontal cut: I(q_par) averaged over the q_z band centered at cut_qz
    // Vertical cut:   I(q_z)   averaged over the q_par band centered at cut_qpar
    void compute_cuts() {
        hcut.resize(smooth.cols());
        vcut.resize(smooth.rows);
        vcut_af.resize(smooth.rows);
        hcut_rows = gisaxs::qmap_cut_horizontal(smooth, cut_qz * 0.1, cut_width_qz * 0.1, hcut.data());
        vcut_cols = gisaxs::qmap_cut_vertical(smooth, cut_qpar * 0.1, cut_width_qpar * 0.1, vcut.data());
        for (size_t r = 0; r < smooth.rows; ++r) vcut_af[r] = alpha_f_deg(qz_nm[r]);
    }

    // Positions the horizontal cut according to its mode and keeps both cuts inside the map
    void update_cut_positions() {
        if (!cut_init) {
            cut_qpar = res_qpar_min + 0.25 * (res_qpar_max - res_qpar_min);
            cut_qz = res_qz_min + 0.5 * (res_qz_max - res_qz_min);
            cut_init = true;
        }
        if (cut_mode == CutMode_Yoneda && eval_model.dwba) {
            cut_qz = yoneda_qz_nm();
        } else if (cut_mode == CutMode_AlphaF) {
            cut_qz = qz_nm_from_alpha_f(cut_alpha_f_deg);
        } else if (cut_mode == CutMode_FilmYoneda && eval_model.dwba && eval_model.graded) {
            cut_qz = film_yoneda_qz_nm();
        }
        cut_qz   = CLAMP(cut_qz,   res_qz_min,   res_qz_max);
        cut_qpar = CLAMP(cut_qpar, res_qpar_min, res_qpar_max);
    }

    // ---------------------------------------------------------------------------------------------
    // Export
    // ---------------------------------------------------------------------------------------------

    static bool export_columns(const float* const* cols, const str_t* names, size_t num_cols, size_t n) {
        char path_buf[2048];
        if (!application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, STR_LIT("csv"))) return false;
        const str_t path = {path_buf, strnlen(path_buf, sizeof(path_buf))};
        if (!md_csv_write_to_file(cols, names, num_cols, n, path)) {
            MD_LOG_ERROR("Scattering: failed to write '%.*s'", (int)path.len, path.ptr);
            return false;
        }
        return true;
    }

    void export_cut(bool horizontal) {
        const size_t n = horizontal ? smooth.cols() : smooth.rows;
        std::vector<float> x(n), y(n), a(n);
        for (size_t i = 0; i < n; ++i) {
            x[i] = (float)(horizontal ? ring_q_nm[i] : qz_nm[i]);
            y[i] = (float)(horizontal ? hcut[i] : vcut[i]);
            a[i] = horizontal ? (float)alpha_f_deg(cut_qz) : (float)alpha_f_deg(qz_nm[i]);
        }
        const float* cols[3] = {x.data(), y.data(), a.data()};
        const str_t names[3] = {
            horizontal ? STR_LIT("q_par [nm^-1]") : STR_LIT("q_z [nm^-1]"),
            STR_LIT("I [sr^-1]"),
            STR_LIT("alpha_f [deg]"),
        };
        export_columns(cols, names, 3, n);
    }

    // Full map in long format: one line per (ring, row) with the instrument resolution applied
    void export_map() {
        const size_t R = smooth.cols(), Q = smooth.rows, n = R * Q;
        std::vector<float> qp(n), qz(n), af(n), I(n), cnt(n);
        for (size_t r = 0; r < Q; ++r) {
            for (size_t c = 0; c < R; ++c) {
                const size_t i = r * R + c;
                qp[i]  = (float)ring_q_nm[c];
                qz[i]  = (float)qz_nm[r];
                af[i]  = (float)alpha_f_deg(qz_nm[r]);
                I[i]   = smooth.at(r, c);
                cnt[i] = (float)smooth.ring_count[c];
            }
        }
        const float* cols[5] = {qp.data(), qz.data(), af.data(), I.data(), cnt.data()};
        const str_t names[5] = {STR_LIT("q_par [nm^-1]"), STR_LIT("q_z [nm^-1]"), STR_LIT("alpha_f [deg]"), STR_LIT("I [sr^-1]"), STR_LIT("ring points")};
        export_columns(cols, names, 5, n);
    }

    // ---------------------------------------------------------------------------------------------
    // UI: settings
    // ---------------------------------------------------------------------------------------------

    void draw_window() {
        if (!show_window) return;

        ImGui::SetNextWindowSize({1000, 650}, ImGuiCond_FirstUseEver);
        char title[64];
        snprintf(title, sizeof(title), "Scattering (%s)###Scattering", is_neutron() ? "GISANS" : "GISAXS");
        if (!ImGui::Begin(title, &show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            ImGui::End();
            return;
        }

        const float settings_w = 340.0f;
        ImGui::BeginChild("##scat_settings", ImVec2(settings_w, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
        draw_settings();
        ImGui::EndChild();

        ImGui::SameLine();
        ImGui::BeginChild("##scat_plot", ImVec2(0, 0));
        draw_plot();
        ImGui::EndChild();

        ImGui::End();
    }

    void draw_compute_bar() {
        const bool running = compute_state.load() == ComputeState_Running;
        if (running) {
            float frac = 1.0f;
            if (task_system::task_is_running(task_slices)) {
                frac = 0.8f * task_system::task_fraction_complete(task_slices);
            } else if (task_system::task_is_running(task_rings)) {
                frac = 0.8f + 0.2f * task_system::task_fraction_complete(task_rings);
            }
            const float bw = ImGui::GetContentRegionAvail().x;
            ImGui::ProgressBar(frac, ImVec2(bw * 0.7f, 0));
            ImGui::SameLine();
            if (ImGui::Button("Cancel", ImVec2(-1, 0))) {
                compute_cancel = true;
                task_system::task_interrupt(task_slices);
                task_system::task_interrupt(task_rings);
            }
        } else {
            const bool inputs_changed = ctx && (hash_structure_inputs() != structure_hash);
            const bool stale = ctx && (structure_stale || inputs_changed);
            if (stale) ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.75f, 0.45f, 0.1f, 1.0f));
            if (ImGui::Button(ctx ? "Recompute" : "Compute", ImVec2(-1, 0))) {
                start_compute();
            }
            if (stale) {
                ImGui::PopStyleColor();
                ImGui::SetItemTooltip("%s", structure_stale ? "The system changed since the last computation"
                                                            : "Sample or sampling settings changed since the last computation");
            }
        }
        if (status[0]) {
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextDisabled("%s", status);
            ImGui::PopTextWrapPos();
        }
        ImGui::Separator();
    }

    void draw_settings() {
        draw_compute_bar();
        ImGui::PushItemWidth(150.0f);
#if VIAMD_SCATTERING_NEUTRON
        combo("Radiation", &radiation, radiation_lbl);
        ImGui::SetItemTooltip("X-rays scatter from the electron density, neutrons from the nuclear scattering length\n"
                              "density. Each has its own particle model, media and beam settings.");
#endif
        draw_beam_settings();
        draw_sample_settings();
        draw_environment_settings();
        draw_detector_settings();
        draw_sampling_settings();
        draw_view_settings();
        ImGui::PopItemWidth();
    }

    void draw_beam_settings() {
        if (!ImGui::CollapsingHeader("Beam", ImGuiTreeNodeFlags_DefaultOpen)) return;
        if (!is_neutron()) {
            const bool preset_valid = beam_preset >= 0 && beam_preset < (int)ARRAY_SIZE(beam_presets);
            char preview[128];
            if (preset_valid) {
                snprintf(preview, sizeof(preview), "%s%s", beam_presets[beam_preset].name, preset_matches(beam_presets[beam_preset]) ? "" : " (modified)");
            } else {
                snprintf(preview, sizeof(preview), "Custom");
            }
            if (ImGui::BeginCombo("Setup", preview)) {
                for (int i = 0; i < (int)ARRAY_SIZE(beam_presets); ++i) {
                    if (ImGui::Selectable(beam_presets[i].name, beam_preset == i)) apply_preset(i);
                    ImGui::SetItemTooltip("%s", beam_presets[i].description);
                }
                ImGui::EndCombo();
            }
            ImGui::SetItemTooltip("Experimental setup: beam, substrate, detector, resolution and cut positions");
            ImGui::InputFloat("Energy (keV)", &energy_kev, 0, 0, "%.3f");
            energy_kev = CLAMP(energy_kev, 0.1f, 1000.0f);
            ImGui::SameLine();
            ImGui::TextDisabled("%.4f Å", wavelength());
        }
#if VIAMD_SCATTERING_NEUTRON
        else {
            ImGui::InputFloat("Wavelength (Å)", &neu_wavelength, 0, 0, "%.3f");
            neu_wavelength = CLAMP(neu_wavelength, 0.5f, 50.0f);
            // E = h^2 / (2 m lambda^2) = 81.8042 meV Å^2 / lambda^2, v = h / (m lambda) = 3956 m/s Å / lambda
            ImGui::TextDisabled("Energy: %.3f meV, %.0f m/s", 81.8042 / (neu_wavelength * neu_wavelength), 3956.03 / neu_wavelength);
        }
#endif
        ImGui::SliderFloat("Incidence (°)", &alpha_i_deg, 0.01f, 2.0f, "%.3f");
        if (dwba) {
            const double ac = critical_angle_deg();
            if (ac > 0.0) {
                ImGui::TextDisabled("Substrate critical angle %.3f° (%s)", ac, alpha_i_deg < ac ? "below: evanescent" : "above");
            }
        }
    }

    void draw_sample_settings() {
        if (!ImGui::CollapsingHeader("Sample", ImGuiTreeNodeFlags_DefaultOpen)) return;
        ImGui::InputQuery("Selection", filter, sizeof(filter), filter_valid, filter_err);
#if VIAMD_SCATTERING_NEUTRON
        if (is_neutron()) {
            draw_neutron_particles();
        } else
#endif
        {
            combo("Electrons", &electron_mode, electron_mode_lbl);
            ImGui::SetItemTooltip("Number of electrons of each particle (the scattering weight)");
            if (electron_mode == ElectronMode_Constant) {
                ImGui::InputFloat("e / particle", &electrons_per_particle, 0, 0, "%.1f");
                electrons_per_particle = MAX(electrons_per_particle, 0.0f);
            } else if (electron_mode == ElectronMode_Mass) {
                ImGui::InputFloat("e / Da", &electrons_per_dalton, 0, 0, "%.4f");
                ImGui::SetItemTooltip("Electrons per dalton of mass. Cellulose 0.530, water 0.555, pure C/N/O 0.5");
            }
            ImGui::InputFloat("Material density (e/Å³)", &material_density, 0, 0, "%.4f");
            material_density = MAX(material_density, 0.0f);
            ImGui::SetItemTooltip("Electron density of the particle material. Each particle displaces (electrons /\n"
                                  "material density) of the ambient medium, which sets the contrast. 0: no displacement.\n"
                                  "Cellulose (1.5 g/cm³): 0.478");
            ImGui::InputDouble("Material beta", &material_beta, 0, 0, "%.3e");
            material_beta = MAX(material_beta, 0.0);
            ImGui::SetItemTooltip("Absorption (imaginary part of the refractive index) of the particle material at the\n"
                                  "material density. Only used for the film in the graded DWBA. Cellulose ~2e-9 at 12.8 keV.");
        }
        combo("Gaussian width", &sigma_mode, sigma_mode_lbl);
        if (sigma_mode == SigmaMode_Constant) {
            ImGui::InputFloat("Sigma (Å)", &sigma_constant, 0, 0, "%.2f");
            sigma_constant = MAX(sigma_constant, 0.0f);
        }

        combo("Density model", &density_model, density_model_lbl);
        ImGui::SetItemTooltip("Beads: one isotropic Gaussian per particle.\n"
                              "Fibrils: slices of beads (one center bead + surrounding beads) define a centerline and an\n"
                              "orientation, and a cross-section template is swept continuously along it. Removes the\n"
                              "artificial scattering of discrete beads (high-q plateau, peak at 2 pi / bead spacing).");
        if (density_model == DensityModel_Fibril) {
            ImGui::Indent();
            ImGui::InputText("Center bead", fib_center_name, sizeof(fib_center_name));
            ImGui::InputInt("Beads per slice", &fib_stride);
            fib_stride = CLAMP(fib_stride, 2, 64);
            combo("Cross-section", &fib_template, fibril_template_lbl);
            if (fib_template == FibrilTemplate_Beads) {
                ImGui::InputFloat("Bead sigma (Å)", &fib_bead_sigma, 0, 0, "%.2f");
                fib_bead_sigma = MAX(fib_bead_sigma, 0.5f);
                ImGui::SetItemTooltip("One Gaussian per bead at the mean bead layout of the slices, swept along the fibril");
            } else if (fib_template == FibrilTemplate_Disk) {
                ImGui::InputFloat("Radius (Å)", &fib_disk_radius, 0, 0, "%.2f");
                ImGui::InputFloat("Grid spacing (Å)", &fib_disk_spacing, 0, 0, "%.2f");
                fib_disk_radius = MAX(fib_disk_radius, 1.0f);
                fib_disk_spacing = MAX(fib_disk_spacing, 1.0f);
                ImGui::SetItemTooltip("Uniform disk (as the cylinders in the BornAgain script) built from Gaussians with\n"
                                      "sigma = spacing / 2 on a hexagonal grid");
            } else {
                ImGui::InputText("##fib_path", fib_template_path, sizeof(fib_template_path));
                ImGui::SameLine();
                if (ImGui::Button("Browse")) {
                    char buf[512];
                    if (application::file_dialog(buf, sizeof(buf), application::FileDialogFlag_Open, STR_LIT("txt"))) {
                        snprintf(fib_template_path, sizeof(fib_template_path), "%s", buf);
                    }
                }
                ImGui::SetItemTooltip("Text file with lines 'gauss u v weight sigma' (and optionally 'anchor u v' per bead), Å");
            }
            ImGui::InputFloat("Sample spacing (Å)", &fib_ds, 0, 0, "%.2f");
            fib_ds = MAX(fib_ds, 0.0f);
            ImGui::SetItemTooltip("Spacing of the cross-section samples along the fibril. 0 = smallest template sigma,\n"
                                  "which makes the density continuous along the fibril.");
            if (fib_info[0]) {
                ImGui::PushTextWrapPos(0.0f);
                ImGui::TextDisabled("%s", fib_info);
                ImGui::PopTextWrapPos();
            }
            ImGui::Unindent();
        }
    }

    void draw_environment_settings() {
        if (!ImGui::CollapsingHeader("Ambient and substrate", ImGuiTreeNodeFlags_DefaultOpen)) return;

        // --- Ambient medium ---
#if VIAMD_SCATTERING_NEUTRON
        if (is_neutron()) {
            if (ImGui::BeginCombo("Ambient", neutron_background_presets[neu_bg_preset].name)) {
                for (int i = 0; i < (int)ARRAY_SIZE(neutron_background_presets); ++i) {
                    if (ImGui::Selectable(neutron_background_presets[i].name, neu_bg_preset == i)) {
                        neu_bg_preset = i;
                        if (i != NBG_CUSTOM && i != NBG_WATER) neu_bg_sld = neutron_background_presets[i].sld;
                    }
                }
                ImGui::EndCombo();
            }
            if (neu_bg_preset == NBG_WATER) {
                ImGui::SliderFloat("D2O fraction", &neu_d2o, 0.0f, 1.0f, "%.3f");
                ImGui::SetItemTooltip("Volume fraction of D2O in H2O/D2O. The SLD is linear from %.2f (H2O) to %.2f (D2O) x 10⁻⁶ Å⁻².",
                                      neutron::SLD_H2O, neutron::SLD_D2O);
            }
            ImGui::BeginDisabled(neu_bg_preset != NBG_CUSTOM);
            double shown = neutron_bg_sld();
            if (ImGui::InputDouble("SLD (10⁻⁶ Å⁻²)", &shown, 0, 0, "%.4f") && neu_bg_preset == NBG_CUSTOM) neu_bg_sld = shown;
            ImGui::EndDisabled();
            ImGui::TextDisabled("Part of the particle contrast (recompute)");
        } else
#endif
        {
            if (ImGui::BeginCombo("Ambient", ambient_presets[ambient_preset].name)) {
                for (int i = 0; i < (int)ARRAY_SIZE(ambient_presets); ++i) {
                    if (ImGui::Selectable(ambient_presets[i].name, ambient_preset == i)) {
                        ambient_preset = i;
                        if (i != AMBIENT_CUSTOM) ambient_density = ambient_presets[i].electron_density;
                    }
                }
                ImGui::EndCombo();
            }
            ImGui::SetItemTooltip("Medium the beam travels in and the particles are embedded in");
            ImGui::BeginDisabled(ambient_preset != AMBIENT_CUSTOM);
            ImGui::InputDouble("Density (e/Å³)", &ambient_density, 0, 0, "%.5f");
            ImGui::EndDisabled();
            ambient_density = MAX(ambient_density, 0.0);
            const double cf = gisaxs::contrast_factor(xray_setup());
            if (cf < 1.0) ImGui::TextDisabled("Contrast factor %.4f", cf);
            if (cf <= 0.0) ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.3f, 1.0f), "Contrast matched: no particle scattering");
        }

        // --- Substrate ---
        ImGui::Checkbox("Substrate (DWBA)", &dwba);
        ImGui::SetItemTooltip("Unchecked: Born approximation without a substrate");
        ImGui::BeginDisabled(!dwba);
#if VIAMD_SCATTERING_NEUTRON
        if (is_neutron()) {
            if (ImGui::BeginCombo("Material", neutron_substrate_presets[neu_sub_preset].name)) {
                for (int i = 0; i < (int)ARRAY_SIZE(neutron_substrate_presets); ++i) {
                    if (ImGui::Selectable(neutron_substrate_presets[i].name, neu_sub_preset == i)) {
                        neu_sub_preset = i;
                        if (i != NSUB_CUSTOM) {
                            neu_sub_sld = neutron_substrate_presets[i].sld;
                            neu_sub_abs = neutron_substrate_presets[i].abs;
                        }
                    }
                }
                ImGui::EndCombo();
            }
            ImGui::BeginDisabled(neu_sub_preset != NSUB_CUSTOM);
            ImGui::InputDouble("SLD (10⁻⁶ Å⁻²)##sub", &neu_sub_sld, 0, 0, "%.4f");
            ImGui::InputDouble("Absorption (10⁻⁶ Å⁻²)", &neu_sub_abs, 0, 0, "%.3e");
            ImGui::SetItemTooltip("Imaginary part of the SLD, N sigma_a / (2 lambda), which is wavelength independent\n"
                                  "for 1/v absorbers. Negligible for Si and SiO2.");
            ImGui::EndDisabled();
            neu_sub_abs = MAX(neu_sub_abs, 0.0);
        } else
#endif
        {
            if (ImGui::BeginCombo("Material", substrate_presets[sub_preset].name)) {
                for (int i = 0; i < (int)ARRAY_SIZE(substrate_presets); ++i) {
                    if (ImGui::Selectable(substrate_presets[i].name, sub_preset == i)) {
                        sub_preset = i;
                        if (i != SUB_CUSTOM) {
                            sub_density = substrate_presets[i].electron_density;
                            sub_beta = substrate_presets[i].beta;
                        }
                    }
                }
                ImGui::EndCombo();
            }
            ImGui::BeginDisabled(sub_preset != SUB_CUSTOM);
            ImGui::InputDouble("Density (e/Å³)##sub", &sub_density, 0, 0, "%.4f");
            ImGui::InputDouble("Beta", &sub_beta, 0, 0, "%.3e");
            ImGui::SetItemTooltip("Imaginary part of the refractive index. Presets are for ~12 keV.");
            ImGui::EndDisabled();
        }
        ImGui::Checkbox("Below lowest particle", &sub_auto_z);
        ImGui::SetItemTooltip("Place the substrate interface relative to the lowest selected particle");
        if (sub_auto_z) {
            ImGui::InputFloat("Offset (Å)", &sub_z_offset, 0, 0, "%.2f");
            ImGui::SetItemTooltip("Substrate z relative to the lowest selected particle (negative = below)");
            if (ctx) {
                ImGui::SameLine();
                ImGui::TextDisabled("z = %.1f Å", substrate_z());
            }
        } else {
            ImGui::InputFloat("z (Å)", &sub_z, 0, 0, "%.2f");
        }
        ImGui::InputFloat("Roughness (Å)", &sub_roughness, 0, 0, "%.2f");
        sub_roughness = MAX(sub_roughness, 0.0f);
        ImGui::SetItemTooltip("RMS roughness of the substrate interface (Nevot-Croce factor)");
        const md_gisaxs_model_t m = build_model();
        if (!(m.sld_substrate > m.sld_ambient)) {
            ImGui::TextDisabled("No total reflection (substrate SLD below the ambient)");
            ImGui::SetItemTooltip("The beam enters from the ambient side. Geometries where the beam enters through the\n"
                                  "substrate (e.g. a solid/liquid interface) are not modelled.");
        }
        ImGui::Checkbox("Film in reference medium (graded DWBA)", &film_graded);
        ImGui::SetItemTooltip("Include the laterally averaged density of the particles (the density profile)\n"
                              "in the DWBA reference medium: ambient / graded film / substrate.\n"
                              "Gives refraction inside the film and the film Yoneda peak.");
        if (film_graded && !prof_rho.empty()) {
            const double d = film_sld_rel();
            if (d > 0.0) ImGui::TextDisabled("Film critical angle %.3f°", gisaxs::critical_angle(d, wavelength()) * RAD_TO_DEG);
        }
        ImGui::EndDisabled();
    }

    void draw_detector_settings() {
        if (!ImGui::CollapsingHeader("Detector and resolution")) return;
        bool changed = false;
        changed |= ImGui::InputFloat("Resolution q_par (nm⁻¹)", &res_fwhm_qpar, 0, 0, "%.4f");
        ImGui::SetItemTooltip("Instrument resolution (FWHM) along q_par, applied as a Gaussian convolution");
        changed |= ImGui::InputFloat("Resolution q_z (nm⁻¹)", &res_fwhm_qz, 0, 0, "%.4f");
        ImGui::SetItemTooltip("Instrument resolution (FWHM) along q_z, applied as a Gaussian convolution");
        res_fwhm_qpar = MAX(res_fwhm_qpar, 0.0f);
        res_fwhm_qz = MAX(res_fwhm_qz, 0.0f);
        if (changed) res_dirty = true;

        ImGui::Checkbox("Show as detector image", &view_detector);
        ImGui::SetItemTooltip("Project the (azimuthally averaged) result onto a flat area detector,\n"
                              "with module gaps and beamstops");
        if (!view_detector) return;
        ImGui::InputDouble("Distance (mm)", &det.sdd_mm, 0, 0, "%.1f");
        ImGui::InputDouble("Pixel size (mm)", &det.pixel_mm, 0, 0, "%.4f");
        ImGui::InputInt("Pixels horizontal", &det.npx_h);
        ImGui::InputInt("Pixels vertical", &det.npx_v);
        ImGui::InputDouble("Direct beam x (px)", &det.beam_x_px, 0, 0, "%.1f");
        ImGui::InputDouble("Direct beam y (px)", &det.beam_y_px, 0, 0, "%.1f");
        ImGui::SetItemTooltip("Row of the direct beam, counted from the bottom edge of the detector");
        ImGui::SliderInt("Binning", &det.binning, 1, 16);
        ImGui::Checkbox("Module gaps", &det.gaps);
        if (det.gaps) {
            ImGui::InputInt("Module width (px)", &det.module_w);
            ImGui::InputInt("Module height (px)", &det.module_h);
            ImGui::InputInt("Gap horizontal (px)", &det.gap_w);
            ImGui::InputInt("Gap vertical (px)", &det.gap_h);
        }
        ImGui::InputInt("Direct beamstop (px)", &det.bs_direct_px);
        ImGui::InputInt("Specular beamstop (px)", &det.bs_specular_px);
        det.sdd_mm = MAX(det.sdd_mm, 1.0);
        det.pixel_mm = MAX(det.pixel_mm, 1.0e-3);
        det.npx_h = CLAMP(det.npx_h, 1, 8192);
        det.npx_v = CLAMP(det.npx_v, 1, 8192);
        det.module_w = MAX(det.module_w, 0); det.module_h = MAX(det.module_h, 0);
        det.gap_w = MAX(det.gap_w, 0); det.gap_h = MAX(det.gap_h, 0);
        det.bs_direct_px = MAX(det.bs_direct_px, 0); det.bs_specular_px = MAX(det.bs_specular_px, 0);
    }

    void draw_sampling_settings() {
        if (!ImGui::CollapsingHeader("Sampling")) return;
        ImGui::InputFloat("q_par max (nm⁻¹)", &q_par_max_nm, 0, 0, "%.3f");
        ImGui::InputFloat("q_z max (nm⁻¹)", &q_z_max_nm, 0, 0, "%.3f");
        q_par_max_nm = MAX(q_par_max_nm, 0.01f);
        q_z_max_nm = MAX(q_z_max_nm, 0.01f);
        ImGui::SliderInt("q_z rows", &num_qz, 16, 1024);
        ImGui::SliderFloat("Oversampling", &oversampling, 1.5f, 4.0f, "%.2f");
        ImGui::SetItemTooltip("Grid oversampling relative to q max (in-plane and z). The B-spline aliasing error at\n"
                              "q max is ~4%% at 2, ~0.4%% at 3 and ~0.1%% at 4, at the cost of memory and time.");
        ImGui::SliderFloat("Ring width (Δq/q)", &ring_rel_width, 0.01f, 1.0f, "%.3f", ImGuiSliderFlags_Logarithmic);
        ring_rel_width = CLAMP(ring_rel_width, 0.001f, 10.0f);
        ImGui::SetItemTooltip("The periodic box restricts q_par to the reciprocal lattice (2π i / Lx, 2π j / Ly).\n"
                              "Lattice shells closer than this fraction of q_par are averaged into one ring.\n"
                              "At low q_par every shell is its own ring (exact q_par, few points), rings are never wider than\n"
                              "the lattice spacing max(2π/Lx, 2π/Ly). Smaller: finer q_par sampling, noisier rings, more memory.");
        ImGui::InputInt("Max slices", &max_slices);
        max_slices = CLAMP(max_slices, 8, 8192);
        if (ctx) {
            md_gisaxs_info_t info;
            md_gisaxs_get_info(ctx, &info);
            ImGui::TextDisabled("Grid %i x %i, dx = %.2f Å", info.nx, info.ny, info.dx);
            ImGui::TextDisabled("%zu slices, dz = %.2f Å", info.num_slices, info.dz);
            ImGui::TextDisabled("%zu rings, lowest q_par %.4f nm⁻¹", info.num_rings, info.q_par_min * 10.0);
            ImGui::SetItemTooltip("The lowest q_par is set by the box size, 2π / max(Lx, Ly). Lower q needs a larger box.");
            ImGui::TextDisabled("Ring matrices: %.1f MB", info.matrix_bytes / (1024.0 * 1024.0));
        }
    }

    void draw_view_settings() {
        if (!ImGui::CollapsingHeader("Display and cuts", ImGuiTreeNodeFlags_DefaultOpen)) return;
        if (ImGui::Checkbox("Log scale", &log_scale)) res_dirty = true;
        if (log_scale) {
            ImGui::SameLine();
            ImGui::SetNextItemWidth(-1);
            if (ImGui::SliderFloat("##decades", &decades, 1.0f, 12.0f, "%.1f decades")) res_dirty = true;
        }
        if (ImGui::BeginCombo("Colormap", ImPlot::GetColormapName(colormap))) {
            for (int i = 0; i < ImPlot::GetColormapCount(); ++i) {
                if (ImGui::Selectable(ImPlot::GetColormapName(i), colormap == i)) colormap = i;
            }
            ImGui::EndCombo();
        }
        ImGui::Checkbox("Density profile", &show_profile);
        ImGui::SameLine();
        ImGui::Checkbox("Cuts", &show_cuts);
        if (!show_cuts) return;

        combo("Horizontal cut", &cut_mode, cut_mode_lbl);
        ImGui::SetItemTooltip("How the horizontal cut (I vs q_par) is positioned.\nDragging the line on the map switches to a fixed q_z.");
        if (cut_mode == CutMode_AlphaF) {
            ImGui::InputFloat("alpha_f (°)", &cut_alpha_f_deg, 0, 0, "%.4f");
        } else if (cut_mode == CutMode_Qz) {
            ImGui::InputDouble("q_z (nm⁻¹)", &cut_qz, 0, 0, "%.4f");
        }
        ImGui::InputFloat("q_z band (nm⁻¹)", &cut_width_qz, 0, 0, "%.4f");
        ImGui::SetItemTooltip("Width of the integration band. 0 uses the nearest row.");
        ImGui::InputDouble("Vertical cut q_par (nm⁻¹)", &cut_qpar, 0, 0, "%.4f");
        ImGui::InputFloat("q_par band (nm⁻¹)", &cut_width_qpar, 0, 0, "%.4f");
        ImGui::SetItemTooltip("Width of the integration band. 0 interpolates between the two nearest rings.");
        cut_width_qz = MAX(cut_width_qz, 0.0f);
        cut_width_qpar = MAX(cut_width_qpar, 0.0f);
        ImGui::Checkbox("Log q_par axis", &cut_log_x);
        ImGui::SameLine();
        ImGui::Checkbox("Vertical vs alpha_f", &vcut_vs_alpha_f);
        ImGui::Checkbox("Bead self scattering", &show_bead_self);
        ImGui::SetItemTooltip("Show the (Born) self scattering of the discrete beads in the horizontal cut.\n"
                              "Where the computed intensity follows this line, bead discreteness dominates.");
        if (has_result()) {
            ImGui::TextDisabled("Averaged over %zu rows / %zu rings", hcut_rows, vcut_cols);
        }
        ImGui::BeginDisabled(!has_result());
        if (ImGui::Button("Export cuts...")) ImGui::OpenPopup("##scat_export");
        if (ImGui::BeginPopup("##scat_export")) {
            if (ImGui::Selectable("Horizontal cut I(q_par) (CSV)")) export_cut(true);
            if (ImGui::Selectable("Vertical cut I(q_z) (CSV)")) export_cut(false);
            if (ImGui::Selectable("Full map I(q_par, q_z) (CSV)")) export_map();
            ImGui::EndPopup();
        }
        ImGui::EndDisabled();
    }

#if VIAMD_SCATTERING_NEUTRON
    void draw_neutron_particles() {
        if (ImGui::BeginCombo("Scattering length", neutron_weight_lbl[neu_weight_mode])) {
            for (int i = 0; i < (int)ARRAY_SIZE(neutron_weight_lbl); ++i) {
                if (ImGui::Selectable(neutron_weight_lbl[i], neu_weight_mode == i)) neu_weight_mode = i;
            }
            ImGui::EndCombo();
        }
        ImGui::SetItemTooltip("Coherent scattering length b of each particle. The scattering weight is the excess\n"
                              "b - SLD_background * V, exact per particle (mixtures, labelling, contrast matching).");
        if (neu_weight_mode == NeutronWeight_Constant) {
            ImGui::InputFloat("b, all ¹H (fm)", &neu_b_protiated, 0, 0, "%.2f");
            ImGui::SetItemTooltip("Coherent scattering length per particle with every hydrogen as ¹H (including the\n"
                                  "exchangeable ones). Cellulose: 31.50 fm per anhydroglucose unit (C6H10O5).");
            ImGui::InputFloat("Non-exchangeable H", &neu_n_h, 0, 0, "%.1f");
            ImGui::SetItemTooltip("Hydrogens bound to C per particle, labelled by the deuteration. Cellulose: 7 per unit.");
            ImGui::InputFloat("Exchangeable H", &neu_n_ex, 0, 0, "%.1f");
            ImGui::SetItemTooltip("Hydrogens bound to N, O or S per particle, which exchange with the reservoir.\n"
                                  "Cellulose: 3 (hydroxyl) per anhydroglucose unit.");
            ImGui::InputFloat("Volume (Å³)", &neu_volume, 0, 0, "%.1f");
            ImGui::SetItemTooltip("Volume displaced by one particle. Cellulose: 179.5 Å³ per anhydroglucose unit at 1.5 g/cm³.");
            neu_n_h = MAX(neu_n_h, 0.0f);
            neu_n_ex = MAX(neu_n_ex, 0.0f);
            neu_volume = MAX(neu_volume, 0.0f);
        } else {
            ImGui::InputFloat("Mass density (g/cm³)", &neu_mass_density, 0, 0, "%.3f");
            ImGui::SetItemTooltip("Volume of each atom = mass / density (hydrogen isotopes with the ¹H mass).\n"
                                  "Hydrogen with a deuterium mass in the topology is taken as D. Hydrogen bound\n"
                                  "to N, O or S (from the bonds) is exchangeable.");
            neu_mass_density = MAX(neu_mass_density, 0.01f);
        }
        ImGui::SliderFloat("Deuteration", &neu_deuteration, 0.0f, 1.0f, "%.3f");
        ImGui::SetItemTooltip("D fraction of the non-exchangeable hydrogen (synthetic labelling)");
        ImGui::Checkbox("H/D exchange", &neu_exchange);
        ImGui::SetItemTooltip("Exchangeable hydrogen (bound to N, O or S) takes the D fraction of the reservoir:\n"
                              "the background water, or e.g. D2O vapour for films measured in a humidity cell.");
        if (neu_exchange) {
            const bool water = neu_bg_preset == NBG_WATER;
            if (water) ImGui::Checkbox("Exchange with the background", &neu_exchange_follow);
            if (!water || !neu_exchange_follow) {
                ImGui::SliderFloat("Reservoir D fraction", &neu_exchange_d, 0.0f, 1.0f, "%.3f");
            } else {
                ImGui::TextDisabled("Reservoir D fraction: %.3f", neu_d2o);
            }
        }
        ImGui::Checkbox("Incoherent background", &neu_incoherent);
        ImGui::SetItemTooltip("Flat (Born) incoherent background sum(sigma_inc) / (4 pi A), dominated by ¹H (80 b).\n"
                              "Added above the sample horizon. Transmission (DWBA) factors are not applied to it.");
        if (ctx && ctx_radiation == Radiation_Neutron && neu_sum_v > 0.0) {
            const double sld_sel = neutron::sld(neu_sum_b, neu_sum_v);
            ImGui::TextDisabled("Selection SLD: %.3f x 10⁻⁶ Å⁻²", sld_sel);
            ImGui::TextDisabled("Contrast: %.3f x 10⁻⁶ Å⁻²", sld_sel - neu_bg_sld_used);
            ImGui::SetItemTooltip("Mean SLD of the selection minus the background SLD, at the last compute");
            if (inc_level > 0.0) ImGui::TextDisabled("Incoherent: %.3e sr⁻¹", inc_level);
            if (neu_weight_mode == NeutronWeight_Element) {
                ImGui::TextDisabled("%zu exchangeable H, %zu unknown elements", neu_num_exchangeable, neu_num_unknown);
            }
        }
    }
#endif

    // ---------------------------------------------------------------------------------------------
    // UI: plots
    // ---------------------------------------------------------------------------------------------

    void draw_detector(float w, float h, ImVec4 col_h, ImVec4 col_v) {
        if (!ImPlot::BeginPlot("##scat_detector", ImVec2(w, h), ImPlotFlags_NoLegend | ImPlotFlags_Equal)) return;
        const double qy0 = det_img.qy_min * 10.0, qy1 = det_img.qy_max * 10.0;
        const double qz0 = det_img.qz_min * 10.0, qz1 = det_img.qz_max * 10.0;
        ImPlot::SetupAxis(ImAxis_X1, "q_y [nm⁻¹]");
        ImPlot::SetupAxis(ImAxis_Y1, "q_z [nm⁻¹]");
        ImPlot::SetupAxesLimits(qy0, qy1, qz0, qz1, ImPlotCond_Once);
        ImPlot::SetupFinish();

        ImPlot::PlotHeatmap("##det", det_display.data(), det_img.height, det_img.width, res_min, res_max, nullptr,
                            ImPlotPoint(qy0, qz0), ImPlotPoint(qy1, qz1));

        if (eval_model.dwba) {
            const double horizon = horizon_qz_nm();
            ImPlot::SetNextLineStyle(ImVec4(1, 1, 1, 0.4f), 1.0f);
            ImPlot::PlotInfLines("Horizon", &horizon, 1, ImPlotInfLinesFlags_Horizontal);
        }
        if (show_cuts) {
            ImPlot::SetNextLineStyle(ImVec4(col_h.x, col_h.y, col_h.z, 0.8f), 1.0f);
            ImPlot::PlotInfLines("##hcut", &cut_qz, 1, ImPlotInfLinesFlags_Horizontal);
            const double v[2] = {-cut_qpar, cut_qpar};
            ImPlot::SetNextLineStyle(ImVec4(col_v.x, col_v.y, col_v.z, 0.8f), 1.0f);
            ImPlot::PlotInfLines("##vcut", v, 2);
        }

        if (ImPlot::IsPlotHovered() && det_img.width > 0 && det_img.height > 0) {
            // Map the mouse back to a detector pixel through the same linear axes used for display
            const ImPlotPoint mp = ImPlot::GetPlotMousePos();
            const double fx = (mp.x - qy0) / (qy1 - qy0);
            const double fy = (mp.y - qz0) / (qz1 - qz0);
            if (fx >= 0 && fx < 1 && fy >= 0 && fy < 1) {
                const int bx = (int)(fx * det_img.width);
                const int by = (int)(fy * det_img.height);
                double x_mm, y_mm;
                gisaxs::detector_bin_center_mm(det, bx, by, &x_mm, &y_mm);
                gisaxs::DetectorQ q;
                gisaxs::detector_q(det, eval_beam(), x_mm, y_mm, &q);
                const double px = det.beam_x_px + x_mm / det.pixel_mm, py = det.beam_y_px + y_mm / det.pixel_mm;
                const float I = det_img.I[(size_t)by * det_img.width + bx];
                if (I >= 0.0f) {
                    ImGui::SetTooltip("Pixel: %.0f, %.0f\nq_y: %.4f nm⁻¹\nq_z: %.4f nm⁻¹\nq_par: %.4f nm⁻¹\nalpha_f: %.3f°\nI: %.4e sr⁻¹",
                        px, py, q.qy * 10.0, q.qz * 10.0, q.qpar * 10.0, q.alpha_f * RAD_TO_DEG, I);
                } else {
                    ImGui::SetTooltip("Pixel: %.0f, %.0f\nMasked (gap, beamstop, shadow or outside the computed q-range)", px, py);
                }
            }
        }
        ImPlot::EndPlot();
    }

    void draw_plot() {
        if (!has_result()) {
            const bool running = compute_state.load() == ComputeState_Running;
            ImGui::TextDisabled("%s", running ? "Computing..." : "No result yet. Press Compute.");
            return;
        }

        update_cut_positions();
        if (show_cuts) compute_cuts();

        const ImGuiStyle& style = ImGui::GetStyle();
        const float colorbar_w = 70.0f;
        const float avail_h = ImGui::GetContentRegionAvail().y;
        const bool  lower_row = show_cuts || show_profile;
        const float lower_h = lower_row ? MAX(160.0f, avail_h * 0.38f) : 0.0f;
        const float map_h = MAX(100.0f, avail_h - lower_h - (lower_row ? style.ItemSpacing.y : 0.0f));
        const float map_w = ImGui::GetContentRegionAvail().x - colorbar_w - style.ItemSpacing.x;

        const ImVec4 col_h = ImVec4(1.0f, 0.35f, 0.35f, 1.0f);   // horizontal cut (fixed q_z)
        const ImVec4 col_v = ImVec4(0.35f, 0.75f, 1.0f, 1.0f);   // vertical cut (fixed q_par)
        const size_t R = smooth.cols(), Q = smooth.rows;

        ImPlot::PushColormap(colormap);
        ImPlot::PushStyleColor(ImPlotCol_PlotBg, ImPlot::SampleColormap(0.0f, colormap));

        if (view_detector && !det_display.empty()) {
            draw_detector(map_w, map_h, col_h, col_v);
        } else if (ImPlot::BeginPlot("##scat_map", ImVec2(map_w, map_h), ImPlotFlags_NoLegend)) {
            ImPlot::SetupAxis(ImAxis_X1, "q_par [nm⁻¹]");
            ImPlot::SetupAxis(ImAxis_Y1, "q_z [nm⁻¹]");
            ImPlot::SetupAxesLimits(res_qpar_min, res_qpar_max, res_qz_min, res_qz_max, ImPlotCond_Once);
            ImPlot::SetupAxisLinks(ImAxis_X1, &link_qpar_min, &link_qpar_max);
            ImPlot::SetupAxisLinks(ImAxis_Y1, &link_qz_min, &link_qz_max);
            ImPlot::SetupFinish();

            // The rings are non-uniform in q_par: one single column heat map per ring, spanning its edges
            for (size_t c = 0; c < R; ++c) {
                ImGui::PushID((int)c);
                ImPlot::PlotHeatmap("##I", map_display.data() + c * Q, (int)Q, 1, res_min, res_max, nullptr,
                                    ImPlotPoint(ring_edge_nm[c], res_qz_min), ImPlotPoint(ring_edge_nm[c + 1], res_qz_max));
                ImGui::PopID();
            }

            if (eval_model.dwba) {
                const double horizon = horizon_qz_nm();
                ImPlot::SetNextLineStyle(ImVec4(1, 1, 1, 0.5f), 1.0f);
                ImPlot::PlotInfLines("Horizon", &horizon, 1, ImPlotInfLinesFlags_Horizontal);
                if (eval_model.sld_substrate > eval_model.sld_ambient) {
                    const double yoneda = yoneda_qz_nm();
                    ImPlot::SetNextLineStyle(ImVec4(1, 0.6f, 0.2f, 0.5f), 1.0f);
                    ImPlot::PlotInfLines("Yoneda", &yoneda, 1, ImPlotInfLinesFlags_Horizontal);
                }
                if (eval_model.graded && film_sld_rel() > 0.0) {
                    const double yf = film_yoneda_qz_nm();
                    ImPlot::SetNextLineStyle(ImVec4(0.5f, 1.0f, 0.5f, 0.5f), 1.0f);
                    ImPlot::PlotInfLines("Film Yoneda", &yf, 1, ImPlotInfLinesFlags_Horizontal);
                }
            }

            if (show_cuts) {
                // Band edges
                if (cut_width_qz > 0.0f) {
                    const double e[2] = {cut_qz - 0.5 * cut_width_qz, cut_qz + 0.5 * cut_width_qz};
                    ImPlot::SetNextLineStyle(ImVec4(col_h.x, col_h.y, col_h.z, 0.4f), 1.0f);
                    ImPlot::PlotInfLines("##hband", e, 2, ImPlotInfLinesFlags_Horizontal);
                }
                if (cut_width_qpar > 0.0f) {
                    const double e[2] = {cut_qpar - 0.5 * cut_width_qpar, cut_qpar + 0.5 * cut_width_qpar};
                    ImPlot::SetNextLineStyle(ImVec4(col_v.x, col_v.y, col_v.z, 0.4f), 1.0f);
                    ImPlot::PlotInfLines("##vband", e, 2);
                }
                if (ImPlot::DragLineY(0, &cut_qz, col_h, 1.5f)) {
                    cut_mode = CutMode_Qz;
                }
                ImPlot::DragLineX(1, &cut_qpar, col_v, 1.5f);
                ImPlot::TagY(cut_qz, col_h, "%.3f", cut_qz);
                ImPlot::TagX(cut_qpar, col_v, "%.3f", cut_qpar);
            }

            if (ImPlot::IsPlotHovered()) {
                const ImPlotPoint mp = ImPlot::GetPlotMousePos();
                const long col = gisaxs::ring_at(smooth, mp.x * 0.1);
                const long row = smooth.dqz > 0 ? (long)floor((mp.y * 0.1 - smooth.qz0) / smooth.dqz + 0.5) : -1;
                if (col >= 0 && col < (long)R && row >= 0 && row < (long)Q) {
                    ImGui::SetTooltip("q_par: %.4f nm⁻¹ (ring %.4f nm⁻¹, %u points)\nq_z: %.4f nm⁻¹\nalpha_f: %.3f°\nI: %.4e sr⁻¹",
                        mp.x, ring_q_nm[col], smooth.ring_count[col], mp.y, alpha_f_deg(mp.y), smooth.at((size_t)row, (size_t)col));
                }
            }
            ImPlot::EndPlot();
        }
        ImPlot::PopStyleColor();

        ImGui::SameLine();
        ImPlot::ColormapScale(log_scale ? "log10 I [sr⁻¹]" : "I [sr⁻¹]", res_min, res_max, ImVec2(colorbar_w, map_h), log_scale ? "%.1f" : "%.1e");
        ImPlot::PopColormap();

        if (!lower_row) return;

        const int num_plots = (show_cuts ? 2 : 0) + (show_profile ? 1 : 0);
        const float plot_w = (ImGui::GetContentRegionAvail().x - style.ItemSpacing.x * (num_plots - 1)) / num_plots;
        const ImPlotScale y_scale = log_scale ? ImPlotScale_Log10 : ImPlotScale_Linear;

        if (show_cuts) {
            // Horizontal cut: I(q_par)
            char title[128];
            snprintf(title, sizeof(title), "q_z = %.3f nm⁻¹ (α_f = %.3f°)###hcut", cut_qz, alpha_f_deg(cut_qz));
            if (ImPlot::BeginPlot(title, ImVec2(plot_w, lower_h), ImPlotFlags_NoLegend)) {
                ImPlot::SetupAxis(ImAxis_X1, "q_par [nm⁻¹]");
                ImPlot::SetupAxis(ImAxis_Y1, "I [sr⁻¹]", ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisScale(ImAxis_Y1, y_scale);
                if (cut_log_x) ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                ImPlot::SetupAxisLinks(ImAxis_X1, &link_qpar_min, &link_qpar_max);
                ImPlot::SetupFinish();
                size_t b = 0, e = R;
                if (log_scale) positive_range(hcut.data(), R, &b, &e);
                if (e > b) {
                    ImPlot::SetNextLineStyle(col_h, 1.5f);
                    ImPlot::PlotLine("I", ring_q_nm.data() + b, hcut.data() + b, (int)(e - b));
                }
                if (show_bead_self && !bead_self.empty()) {
                    // Born self scattering of the discrete beads: sum_j w_j^2 exp(-(q_par^2 + q_z^2) sigma_j^2) / A.
                    // Where the computed curve approaches it, the signal is dominated by bead discreteness.
                    const double qz_a = cut_qz * 0.1;
                    const double scl = eval_model.intensity_scale / bead_self_area;   // as the result
                    std::vector<double> self_I(R);
                    for (size_t c = 0; c < R; ++c) {
                        const double qp = smooth.ring_q[c];
                        double sum = 0.0;
                        for (const auto& g : bead_self) sum += g.second * exp(-(qp * qp + qz_a * qz_a) * g.first * g.first);
                        self_I[c] = sum * scl;
                    }
                    ImPlot::SetNextLineStyle(ImVec4(0.8f, 0.8f, 0.8f, 0.7f), 1.0f);
                    ImPlot::PlotLine("Bead self scattering", ring_q_nm.data(), self_I.data(), (int)R);
                }
                ImPlot::SetNextLineStyle(ImVec4(col_v.x, col_v.y, col_v.z, 0.6f), 1.0f);
                ImPlot::PlotInfLines("##qpar", &cut_qpar, 1);
                ImPlot::EndPlot();
            }

            // Vertical cut: I(q_z), against q_z (linked with the map) or the exit angle (as in most papers)
            ImGui::SameLine();
            snprintf(title, sizeof(title), "q_par = %.3f nm⁻¹###vcut", cut_qpar);
            if (ImPlot::BeginPlot(title, ImVec2(plot_w, lower_h), ImPlotFlags_NoLegend)) {
                const bool use_af = vcut_vs_alpha_f;
                const double* xs = use_af ? vcut_af.data() : qz_nm.data();
                auto to_x = [&](double qz) { return use_af ? alpha_f_deg(qz) : qz; };
                ImPlot::SetupAxis(ImAxis_X1, use_af ? "alpha_f [°]" : "q_z [nm⁻¹]", use_af ? ImPlotAxisFlags_AutoFit : 0);
                ImPlot::SetupAxis(ImAxis_Y1, "I [sr⁻¹]", ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisScale(ImAxis_Y1, y_scale);
                if (!use_af) ImPlot::SetupAxisLinks(ImAxis_X1, &link_qz_min, &link_qz_max);
                ImPlot::SetupFinish();
                size_t b = 0, e = Q;
                if (log_scale) positive_range(vcut.data(), Q, &b, &e);
                if (use_af && eval_model.dwba) {
                    // Only above the sample horizon
                    const double horizon = horizon_qz_nm();
                    while (b < e && qz_nm[b] < horizon) ++b;
                }
                if (e > b) {
                    ImPlot::SetNextLineStyle(col_v, 1.5f);
                    ImPlot::PlotLine("I", xs + b, vcut.data() + b, (int)(e - b));
                }
                const double xcut = to_x(cut_qz);
                ImPlot::SetNextLineStyle(ImVec4(col_h.x, col_h.y, col_h.z, 0.6f), 1.0f);
                ImPlot::PlotInfLines("##qz", &xcut, 1);
                if (eval_model.dwba) {
                    const double horizon = to_x(horizon_qz_nm());
                    ImPlot::SetNextLineStyle(ImVec4(1, 1, 1, 0.4f), 1.0f);
                    ImPlot::PlotInfLines("##horizon", &horizon, 1);
                    if (eval_model.sld_substrate > eval_model.sld_ambient) {
                        const double yoneda = to_x(yoneda_qz_nm());
                        ImPlot::SetNextLineStyle(ImVec4(1, 0.6f, 0.2f, 0.4f), 1.0f);
                        ImPlot::PlotInfLines("##yoneda", &yoneda, 1);
                    }
                }
                ImPlot::EndPlot();
            }
        }

        if (show_profile && !prof_z.empty()) {
            if (show_cuts) ImGui::SameLine();
            if (ImPlot::BeginPlot("Density profile##scat_profile", ImVec2(plot_w, lower_h), ImPlotFlags_NoLegend)) {
                const size_t np = prof_z.size();
                bool neu = false;
#if VIAMD_SCATTERING_NEUTRON
                neu = ctx_radiation == Radiation_Neutron;
#endif
                ImPlot::SetupAxis(ImAxis_X1, "z [Å]");
                ImPlot::SetupAxis(ImAxis_Y1, neu ? "ΔSLD [10⁻⁶ Å⁻²]" : "rho_e [e/Å³]", ImPlotAxisFlags_AutoFit);
                ImPlot::SetupFinish();
                std::vector<double> y(prof_rho);
#if VIAMD_SCATTERING_NEUTRON
                if (neu) for (double& v : y) v /= neutron::SLD_E6_TO_FM_PER_A3;   // fm/Å^3 -> 10^-6 Å^-2
#endif
                ImPlot::PlotLine("Profile", prof_z.data(), y.data(), (int)np);
                if (eval_model.dwba) {
                    const double zs = eval_model.z_substrate;
                    ImPlot::SetNextLineStyle(ImVec4(1, 0.6f, 0.2f, 0.8f), 1.0f);
                    ImPlot::PlotInfLines("Substrate", &zs, 1);
                }
                ImPlot::EndPlot();
            }
        }
    }

    // ---------------------------------------------------------------------------------------------
    // Serialization
    // ---------------------------------------------------------------------------------------------

    void serialize(viamd::serialization_state_t& state) {
        viamd::write_section_header(state, STR_LIT("Scattering"));
        viamd::write_str(state, STR_LIT("Filter"), str_from_cstr(filter));
#if VIAMD_SCATTERING_NEUTRON
        viamd::write_int(state, STR_LIT("Radiation"), radiation);
        viamd::write_int(state, STR_LIT("NeuWeightMode"), neu_weight_mode);
        viamd::write_flt(state, STR_LIT("NeuBProtiated"), neu_b_protiated);
        viamd::write_flt(state, STR_LIT("NeuNH"), neu_n_h);
        viamd::write_flt(state, STR_LIT("NeuNEx"), neu_n_ex);
        viamd::write_flt(state, STR_LIT("NeuVolume"), neu_volume);
        viamd::write_flt(state, STR_LIT("NeuMassDensity"), neu_mass_density);
        viamd::write_flt(state, STR_LIT("NeuDeuteration"), neu_deuteration);
        viamd::write_bool(state, STR_LIT("NeuExchange"), neu_exchange);
        viamd::write_bool(state, STR_LIT("NeuExchangeFollow"), neu_exchange_follow);
        viamd::write_flt(state, STR_LIT("NeuExchangeD"), neu_exchange_d);
        viamd::write_bool(state, STR_LIT("NeuIncoherent"), neu_incoherent);
        viamd::write_int(state, STR_LIT("NeuBackgroundPreset"), neu_bg_preset);
        viamd::write_dbl(state, STR_LIT("NeuBackgroundSld"), neu_bg_sld);
        viamd::write_flt(state, STR_LIT("NeuD2O"), neu_d2o);
        viamd::write_int(state, STR_LIT("NeuSubstratePreset"), neu_sub_preset);
        viamd::write_dbl(state, STR_LIT("NeuSubstrateSld"), neu_sub_sld);
        viamd::write_dbl(state, STR_LIT("NeuSubstrateAbs"), neu_sub_abs);
        viamd::write_flt(state, STR_LIT("NeuWavelength"), neu_wavelength);
#endif
        viamd::write_int(state, STR_LIT("ElectronMode"), electron_mode);
        viamd::write_flt(state, STR_LIT("Electrons"), electrons_per_particle);
        viamd::write_flt(state, STR_LIT("ElectronsPerDalton"), electrons_per_dalton);
        viamd::write_int(state, STR_LIT("SigmaMode"), sigma_mode);
        viamd::write_flt(state, STR_LIT("Sigma"), sigma_constant);
        viamd::write_flt(state, STR_LIT("MaterialDensity"), material_density);
        viamd::write_dbl(state, STR_LIT("MaterialBeta"), material_beta);
        viamd::write_int(state, STR_LIT("DensityModel"), density_model);
        viamd::write_str(state, STR_LIT("FibrilCenter"), str_from_cstr(fib_center_name));
        viamd::write_int(state, STR_LIT("FibrilStride"), fib_stride);
        viamd::write_int(state, STR_LIT("FibrilTemplate"), fib_template);
        viamd::write_flt(state, STR_LIT("FibrilBeadSigma"), fib_bead_sigma);
        viamd::write_flt(state, STR_LIT("FibrilDiskRadius"), fib_disk_radius);
        viamd::write_flt(state, STR_LIT("FibrilDiskSpacing"), fib_disk_spacing);
        viamd::write_flt(state, STR_LIT("FibrilSampleSpacing"), fib_ds);
        viamd::write_str(state, STR_LIT("FibrilTemplatePath"), str_from_cstr(fib_template_path));
        viamd::write_int(state, STR_LIT("BackgroundPreset"), ambient_preset);
        viamd::write_dbl(state, STR_LIT("BackgroundDensity"), ambient_density);
        viamd::write_bool(state, STR_LIT("DWBA"), dwba);
        viamd::write_int(state, STR_LIT("SubstratePreset"), sub_preset);
        viamd::write_dbl(state, STR_LIT("SubstrateDensity"), sub_density);
        viamd::write_dbl(state, STR_LIT("SubstrateBeta"), sub_beta);
        viamd::write_bool(state, STR_LIT("SubstrateAutoZ"), sub_auto_z);
        viamd::write_flt(state, STR_LIT("SubstrateZ"), sub_z);
        viamd::write_flt(state, STR_LIT("SubstrateZOffset"), sub_z_offset);
        viamd::write_flt(state, STR_LIT("SubstrateRoughness"), sub_roughness);
        viamd::write_bool(state, STR_LIT("FilmGraded"), film_graded);
        viamd::write_flt(state, STR_LIT("EnergyKeV"), energy_kev);
        viamd::write_flt(state, STR_LIT("AlphaIDeg"), alpha_i_deg);
        viamd::write_int(state, STR_LIT("BeamPreset"), beam_preset);
        viamd::write_flt(state, STR_LIT("QParMax"), q_par_max_nm);
        viamd::write_flt(state, STR_LIT("QZMax"), q_z_max_nm);
        viamd::write_flt(state, STR_LIT("Oversampling"), oversampling);
        viamd::write_flt(state, STR_LIT("RingRelWidth"), ring_rel_width);
        viamd::write_int(state, STR_LIT("NumQz"), num_qz);
        viamd::write_int(state, STR_LIT("MaxSlices"), max_slices);
        viamd::write_flt(state, STR_LIT("ResFwhmQpar"), res_fwhm_qpar);
        viamd::write_flt(state, STR_LIT("ResFwhmQz"), res_fwhm_qz);
        viamd::write_bool(state, STR_LIT("ViewDetector"), view_detector);
        viamd::write_dbl(state, STR_LIT("DetSDD"), det.sdd_mm);
        viamd::write_dbl(state, STR_LIT("DetPixel"), det.pixel_mm);
        viamd::write_int(state, STR_LIT("DetPxH"), det.npx_h);
        viamd::write_int(state, STR_LIT("DetPxV"), det.npx_v);
        viamd::write_dbl(state, STR_LIT("DetBeamX"), det.beam_x_px);
        viamd::write_dbl(state, STR_LIT("DetBeamY"), det.beam_y_px);
        viamd::write_int(state, STR_LIT("DetBinning"), det.binning);
        viamd::write_bool(state, STR_LIT("DetGaps"), det.gaps);
        viamd::write_int(state, STR_LIT("DetModuleW"), det.module_w);
        viamd::write_int(state, STR_LIT("DetModuleH"), det.module_h);
        viamd::write_int(state, STR_LIT("DetGapW"), det.gap_w);
        viamd::write_int(state, STR_LIT("DetGapH"), det.gap_h);
        viamd::write_int(state, STR_LIT("DetBsDirect"), det.bs_direct_px);
        viamd::write_int(state, STR_LIT("DetBsSpecular"), det.bs_specular_px);
        viamd::write_bool(state, STR_LIT("LogScale"), log_scale);
        viamd::write_flt(state, STR_LIT("Decades"), decades);
        viamd::write_bool(state, STR_LIT("ShowBeadSelf"), show_bead_self);
        viamd::write_bool(state, STR_LIT("ShowCuts"), show_cuts);
        viamd::write_int(state, STR_LIT("CutMode"), cut_mode);
        viamd::write_bool(state, STR_LIT("VCutVsAlphaF"), vcut_vs_alpha_f);
        viamd::write_flt(state, STR_LIT("CutAlphaF"), cut_alpha_f_deg);
        viamd::write_dbl(state, STR_LIT("CutQz"), cut_qz);
        viamd::write_dbl(state, STR_LIT("CutQpar"), cut_qpar);
        viamd::write_flt(state, STR_LIT("CutWidthQz"), cut_width_qz);
        viamd::write_flt(state, STR_LIT("CutWidthQpar"), cut_width_qpar);
    }

    void deserialize(viamd::deserialization_state_t& state) {
        str_t ident, arg;
        while (viamd::next_entry(ident, arg, state)) {
            if      (str_eq(ident, STR_LIT("Filter")))            viamd::extract_to_char_buf(filter, sizeof(filter), arg);
#if VIAMD_SCATTERING_NEUTRON
            else if (str_eq(ident, STR_LIT("Radiation")))         viamd::extract_int(radiation, arg);
            else if (str_eq(ident, STR_LIT("NeuWeightMode")))     viamd::extract_int(neu_weight_mode, arg);
            else if (str_eq(ident, STR_LIT("NeuBProtiated")))     viamd::extract_flt(neu_b_protiated, arg);
            else if (str_eq(ident, STR_LIT("NeuNH")))             viamd::extract_flt(neu_n_h, arg);
            else if (str_eq(ident, STR_LIT("NeuNEx")))            viamd::extract_flt(neu_n_ex, arg);
            else if (str_eq(ident, STR_LIT("NeuVolume")))         viamd::extract_flt(neu_volume, arg);
            else if (str_eq(ident, STR_LIT("NeuMassDensity")))    viamd::extract_flt(neu_mass_density, arg);
            else if (str_eq(ident, STR_LIT("NeuDeuteration")))    viamd::extract_flt(neu_deuteration, arg);
            else if (str_eq(ident, STR_LIT("NeuExchange")))       viamd::extract_bool(neu_exchange, arg);
            else if (str_eq(ident, STR_LIT("NeuExchangeFollow"))) viamd::extract_bool(neu_exchange_follow, arg);
            else if (str_eq(ident, STR_LIT("NeuExchangeD")))      viamd::extract_flt(neu_exchange_d, arg);
            else if (str_eq(ident, STR_LIT("NeuIncoherent")))     viamd::extract_bool(neu_incoherent, arg);
            else if (str_eq(ident, STR_LIT("NeuBackgroundPreset"))) viamd::extract_int(neu_bg_preset, arg);
            else if (str_eq(ident, STR_LIT("NeuBackgroundSld")))  viamd::extract_dbl(neu_bg_sld, arg);
            else if (str_eq(ident, STR_LIT("NeuD2O")))            viamd::extract_flt(neu_d2o, arg);
            else if (str_eq(ident, STR_LIT("NeuSubstratePreset"))) viamd::extract_int(neu_sub_preset, arg);
            else if (str_eq(ident, STR_LIT("NeuSubstrateSld")))   viamd::extract_dbl(neu_sub_sld, arg);
            else if (str_eq(ident, STR_LIT("NeuSubstrateAbs")))   viamd::extract_dbl(neu_sub_abs, arg);
            else if (str_eq(ident, STR_LIT("NeuWavelength")))     viamd::extract_flt(neu_wavelength, arg);
#endif
            else if (str_eq(ident, STR_LIT("ElectronMode")))      viamd::extract_int(electron_mode, arg);
            else if (str_eq(ident, STR_LIT("Electrons")))         viamd::extract_flt(electrons_per_particle, arg);
            else if (str_eq(ident, STR_LIT("ElectronsPerDalton"))) viamd::extract_flt(electrons_per_dalton, arg);
            else if (str_eq(ident, STR_LIT("SigmaMode")))         viamd::extract_int(sigma_mode, arg);
            else if (str_eq(ident, STR_LIT("Sigma")))             viamd::extract_flt(sigma_constant, arg);
            else if (str_eq(ident, STR_LIT("MaterialDensity")))   viamd::extract_flt(material_density, arg);
            else if (str_eq(ident, STR_LIT("MaterialBeta")))      viamd::extract_dbl(material_beta, arg);
            else if (str_eq(ident, STR_LIT("DensityModel")))      viamd::extract_int(density_model, arg);
            else if (str_eq(ident, STR_LIT("FibrilCenter")))      viamd::extract_to_char_buf(fib_center_name, sizeof(fib_center_name), arg);
            else if (str_eq(ident, STR_LIT("FibrilStride")))      viamd::extract_int(fib_stride, arg);
            else if (str_eq(ident, STR_LIT("FibrilTemplate")))    viamd::extract_int(fib_template, arg);
            else if (str_eq(ident, STR_LIT("FibrilBeadSigma")))   viamd::extract_flt(fib_bead_sigma, arg);
            else if (str_eq(ident, STR_LIT("FibrilDiskRadius")))  viamd::extract_flt(fib_disk_radius, arg);
            else if (str_eq(ident, STR_LIT("FibrilDiskSpacing"))) viamd::extract_flt(fib_disk_spacing, arg);
            else if (str_eq(ident, STR_LIT("FibrilSampleSpacing"))) viamd::extract_flt(fib_ds, arg);
            else if (str_eq(ident, STR_LIT("FibrilTemplatePath"))) viamd::extract_to_char_buf(fib_template_path, sizeof(fib_template_path), arg);
            else if (str_eq(ident, STR_LIT("BackgroundPreset")))  viamd::extract_int(ambient_preset, arg);
            else if (str_eq(ident, STR_LIT("BackgroundDensity"))) viamd::extract_dbl(ambient_density, arg);
            else if (str_eq(ident, STR_LIT("DWBA")))              viamd::extract_bool(dwba, arg);
            else if (str_eq(ident, STR_LIT("SubstratePreset")))   viamd::extract_int(sub_preset, arg);
            else if (str_eq(ident, STR_LIT("SubstrateDensity")))  viamd::extract_dbl(sub_density, arg);
            else if (str_eq(ident, STR_LIT("SubstrateBeta")))     viamd::extract_dbl(sub_beta, arg);
            else if (str_eq(ident, STR_LIT("SubstrateAutoZ")))    viamd::extract_bool(sub_auto_z, arg);
            else if (str_eq(ident, STR_LIT("SubstrateZ")))        viamd::extract_flt(sub_z, arg);
            else if (str_eq(ident, STR_LIT("SubstrateZOffset")))  viamd::extract_flt(sub_z_offset, arg);
            else if (str_eq(ident, STR_LIT("SubstrateRoughness"))) viamd::extract_flt(sub_roughness, arg);
            else if (str_eq(ident, STR_LIT("FilmGraded")))        viamd::extract_bool(film_graded, arg);
            else if (str_eq(ident, STR_LIT("EnergyKeV")))         viamd::extract_flt(energy_kev, arg);
            else if (str_eq(ident, STR_LIT("AlphaIDeg")))         viamd::extract_flt(alpha_i_deg, arg);
            else if (str_eq(ident, STR_LIT("BeamPreset")))        viamd::extract_int(beam_preset, arg);
            else if (str_eq(ident, STR_LIT("QParMax")))           viamd::extract_flt(q_par_max_nm, arg);
            else if (str_eq(ident, STR_LIT("QZMax")))             viamd::extract_flt(q_z_max_nm, arg);
            else if (str_eq(ident, STR_LIT("Oversampling")))      viamd::extract_flt(oversampling, arg);
            else if (str_eq(ident, STR_LIT("RingRelWidth")))      viamd::extract_flt(ring_rel_width, arg);
            else if (str_eq(ident, STR_LIT("NumQz")))             viamd::extract_int(num_qz, arg);
            else if (str_eq(ident, STR_LIT("MaxSlices")))         viamd::extract_int(max_slices, arg);
            else if (str_eq(ident, STR_LIT("ResFwhmQpar")))       viamd::extract_flt(res_fwhm_qpar, arg);
            else if (str_eq(ident, STR_LIT("ResFwhmQz")))         viamd::extract_flt(res_fwhm_qz, arg);
            else if (str_eq(ident, STR_LIT("ViewDetector")))      viamd::extract_bool(view_detector, arg);
            else if (str_eq(ident, STR_LIT("DetSDD")))            viamd::extract_dbl(det.sdd_mm, arg);
            else if (str_eq(ident, STR_LIT("DetPixel")))          viamd::extract_dbl(det.pixel_mm, arg);
            else if (str_eq(ident, STR_LIT("DetPxH")))            viamd::extract_int(det.npx_h, arg);
            else if (str_eq(ident, STR_LIT("DetPxV")))            viamd::extract_int(det.npx_v, arg);
            else if (str_eq(ident, STR_LIT("DetBeamX")))          viamd::extract_dbl(det.beam_x_px, arg);
            else if (str_eq(ident, STR_LIT("DetBeamY")))          viamd::extract_dbl(det.beam_y_px, arg);
            else if (str_eq(ident, STR_LIT("DetBinning")))        viamd::extract_int(det.binning, arg);
            else if (str_eq(ident, STR_LIT("DetGaps")))           viamd::extract_bool(det.gaps, arg);
            else if (str_eq(ident, STR_LIT("DetModuleW")))        viamd::extract_int(det.module_w, arg);
            else if (str_eq(ident, STR_LIT("DetModuleH")))        viamd::extract_int(det.module_h, arg);
            else if (str_eq(ident, STR_LIT("DetGapW")))           viamd::extract_int(det.gap_w, arg);
            else if (str_eq(ident, STR_LIT("DetGapH")))           viamd::extract_int(det.gap_h, arg);
            else if (str_eq(ident, STR_LIT("DetBsDirect")))       viamd::extract_int(det.bs_direct_px, arg);
            else if (str_eq(ident, STR_LIT("DetBsSpecular")))     viamd::extract_int(det.bs_specular_px, arg);
            else if (str_eq(ident, STR_LIT("LogScale")))          viamd::extract_bool(log_scale, arg);
            else if (str_eq(ident, STR_LIT("Decades")))           viamd::extract_flt(decades, arg);
            else if (str_eq(ident, STR_LIT("ShowBeadSelf")))      viamd::extract_bool(show_bead_self, arg);
            else if (str_eq(ident, STR_LIT("ShowCuts")))          viamd::extract_bool(show_cuts, arg);
            else if (str_eq(ident, STR_LIT("CutMode")))           viamd::extract_int(cut_mode, arg);
            else if (str_eq(ident, STR_LIT("VCutVsAlphaF")))      viamd::extract_bool(vcut_vs_alpha_f, arg);
            else if (str_eq(ident, STR_LIT("CutAlphaF")))         viamd::extract_flt(cut_alpha_f_deg, arg);
            else if (str_eq(ident, STR_LIT("CutQz")))             { viamd::extract_dbl(cut_qz, arg); cut_init = true; }
            else if (str_eq(ident, STR_LIT("CutQpar")))           { viamd::extract_dbl(cut_qpar, arg); cut_init = true; }
            else if (str_eq(ident, STR_LIT("CutWidthQz")))        viamd::extract_flt(cut_width_qz, arg);
            else if (str_eq(ident, STR_LIT("CutWidthQpar")))      viamd::extract_flt(cut_width_qpar, arg);
        }
        electron_mode = CLAMP(electron_mode, 0, (int)ARRAY_SIZE(electron_mode_lbl) - 1);
        sigma_mode = CLAMP(sigma_mode, 0, (int)ARRAY_SIZE(sigma_mode_lbl) - 1);
        density_model = CLAMP(density_model, 0, (int)ARRAY_SIZE(density_model_lbl) - 1);
        fib_template = CLAMP(fib_template, 0, (int)ARRAY_SIZE(fibril_template_lbl) - 1);
        ambient_preset = CLAMP(ambient_preset, 0, (int)ARRAY_SIZE(ambient_presets) - 1);
        sub_preset = CLAMP(sub_preset, 0, (int)ARRAY_SIZE(substrate_presets) - 1);
        cut_mode = CLAMP(cut_mode, 0, (int)ARRAY_SIZE(cut_mode_lbl) - 1);
        beam_preset = CLAMP(beam_preset, -1, (int)ARRAY_SIZE(beam_presets) - 1);
#if VIAMD_SCATTERING_NEUTRON
        radiation = CLAMP(radiation, 0, (int)ARRAY_SIZE(radiation_lbl) - 1);
        neu_weight_mode = CLAMP(neu_weight_mode, 0, (int)ARRAY_SIZE(neutron_weight_lbl) - 1);
        neu_bg_preset = CLAMP(neu_bg_preset, 0, (int)ARRAY_SIZE(neutron_background_presets) - 1);
        neu_sub_preset = CLAMP(neu_sub_preset, 0, (int)ARRAY_SIZE(neutron_substrate_presets) - 1);
#endif
        res_dirty = true;
    }
};

static ScatteringComponent instance;
