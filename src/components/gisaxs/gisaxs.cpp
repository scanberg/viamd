// GISAXS component
//
// Computes the rotationally (azimuthally) averaged GISAXS intensity I(q_par, q_z) of the current system state.
// The electron density is approximated with one Gaussian per particle (bead / atom). The heavy lifting is done
// in mdlib (md_gisaxs), see md_gisaxs.h for a description of the method.
//
// The computation is split in two stages:
//   1. Structure stage (heavy): slices, FFTs and ring averaged cross spectral matrices. Depends on the particle
//      coordinates, selection, particle model (electrons, Gaussian width) and the q-grid.
//   2. Model stage (cheap): DWBA / Born evaluation on top of stage 1. Depends on the beam, substrate and
//      background (ambient) medium. Re-evaluated automatically when any of those change.
//
// The background medium enters in two places: as the ambient medium in the DWBA (refraction/reflection at the
// substrate) and as a contrast factor for the particles, where every particle is assumed to displace a volume
// of background given by electrons / material electron density. With a common material density for all particles,
// this is a uniform scale (1 - rho_bg / rho_mat)^2 of the intensity, which is why it lives in the cheap stage.

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

#include "fibril_core.h"

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

constexpr double HC_KEV_ANGSTROM = 12.398419843320026;  // lambda [Å] = hc / E [keV]
constexpr double DEG_TO_RAD = 3.14159265358979323846 / 180.0;
constexpr double RAD_TO_DEG = 180.0 / 3.14159265358979323846;

// Electron densities in e/Å^3
struct MediumPreset {
    const char* name;
    double electron_density;
    double beta;    // Absorption (imaginary part of refractive index), only used for substrates
};

// Background (ambient) media
const MediumPreset background_presets[] = {
    {"Vacuum",  0.0,      0.0},
    {"Air",     3.62e-4,  0.0},     // 1.205 kg/m^3, Z/A ~ 0.499
    {"Water",   0.3342,   0.0},
    {"Custom",  0.0,      0.0},
};
constexpr int BG_CUSTOM = 3;

// Substrates, beta given for ~12 keV
const MediumPreset substrate_presets[] = {
    {"Silicon",     0.6991, 3.6e-8},  // 2.329 g/cm^3, mu/rho ~19 cm^2/g at 12 keV
    {"SiO2 (glass)",0.6620, 2.4e-8},  // 2.2 g/cm^3, mu/rho ~13 cm^2/g at 12 keV
    {"Custom",      0.6991, 0.0},
};
constexpr int SUB_CUSTOM = 2;

// Experimental setups. Applying a preset sets the beam, the substrate optical constants, the q-range covered by the
// detector and the default cut positions. Substrate optical constants are given as delta / beta at the preset's
// wavelength, delta is converted to an electron density through delta = lambda^2 r_e rho_e / (2 pi).
struct BeamPreset {
    const char* name;
    const char* description;
    double wavelength_nm;
    double alpha_i_deg;
    // Detector (flat, normal to the direct beam)
    double sdd_mm;
    double pixel_mm;
    int    det_px_h;            // horizontal pixels (beam centered horizontally)
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
        1043, 23.0,
        487, 195, 7, 17,
        20, 12,
        2.0e-9,
    },
};

enum CutMode : int {
    CutMode_Yoneda = 0,     // Horizontal cut follows the substrate Yoneda position
    CutMode_AlphaF = 1,     // Horizontal cut at a fixed exit angle
    CutMode_Qz     = 2,     // Horizontal cut at a fixed q_z
    CutMode_FilmYoneda = 3, // Horizontal cut follows the film Yoneda position (graded DWBA)
};
const char* cut_mode_lbl[] = {"Substrate Yoneda", "Fixed exit angle", "Fixed q_z", "Film Yoneda"};

// Flat area detector, normal to the direct beam
struct Detector {
    double sdd_mm   = 4976.0;
    double pixel_mm = 0.172;
    int    npx_h    = 981;
    int    npx_v    = 1043;
    double beam_x_px = 490.5;       // direct beam column (from the left edge)
    double beam_y_px = 23.0;        // direct beam row (from the bottom edge)
    int    binning  = 4;
    int    module_w = 487, module_h = 195;
    int    gap_w = 7, gap_h = 17;
    bool   gaps = true;
    int    bs_direct_px = 20;
    int    bs_specular_px = 12;
};

enum ElectronMode : int { ElectronMode_Constant = 0, ElectronMode_AtomicNumber = 1, ElectronMode_Mass = 2 };
enum DensityModel : int { DensityModel_Beads = 0, DensityModel_Fibril = 1 };
enum FibrilTemplateMode : int { FibrilTemplate_Beads = 0, FibrilTemplate_Disk = 1, FibrilTemplate_File = 2 };
const char* density_model_lbl[] = {"Beads (Gaussian per particle)", "Fibrils (swept cross-section)"};
const char* fibril_template_lbl[] = {"Bead layout", "Uniform disk", "File"};
enum SigmaMode    : int { SigmaMode_Constant = 0, SigmaMode_Radius = 1 };

const char* electron_mode_lbl[] = {"Constant", "Atomic number (Z)", "From mass"};
const char* sigma_mode_lbl[]    = {"Constant", "From radius (R / sqrt(5))"};

enum ComputeState : int {
    ComputeState_Idle = 0,
    ComputeState_Running,
    ComputeState_Done,
    ComputeState_Failed,
};

}  // namespace

struct Gisaxs : viamd::EventHandler {
    ApplicationState* app_state = nullptr;
    bool show_window = false;

    // --- Selection ---
    char filter[256] = "all";
    char filter_err[256] = "";
    bool filter_valid = true;

    // --- Particle model ---
    int   electron_mode = ElectronMode_Constant;
    float electrons_per_particle = 1900.0f;     // e.g. ~22 glucose units (86 e each) per sCG bead
    int   sigma_mode = SigmaMode_Constant;
    float sigma_constant = 6.0f;                 // Å
    float material_density = 0.478f;             // e/Å^3 (cellulose ~1.5 g/cm^3)
    float electrons_per_dalton = 0.530f;         // cellulose C6H10O5: 86 e / 162.14 Da

    // --- Density model ---
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
    bool   show_bead_self = true;

    // --- Background / ambient ---
    int    bg_preset = 0;
    double bg_density = 0.0;                     // e/Å^3

    // --- Substrate ---
    bool   dwba = true;
    int    sub_preset = 0;
    double sub_density = substrate_presets[0].electron_density;
    double sub_beta = substrate_presets[0].beta;
    bool   sub_auto_z = true;
    float  sub_z = 0.0f;                         // Å
    float  sub_z_offset = 0.0f;                  // Å, offset relative to lowest particle when auto
    float  sub_roughness = 3.0f;                 // Å RMS
    bool   film_graded = true;                   // Laterally averaged film as part of the DWBA reference medium
    double material_beta = 2.0e-9;               // Absorption of the particle material at material_density

    // --- Beam ---
    float energy_kev = 12.0f;
    float alpha_i_deg = 0.20f;

    // --- q-grid (nm^-1 in the UI) ---
    float q_par_max_nm = 2.0f;
    float q_z_max_nm = 2.0f;
    float oversampling = 2.0f;
    int   num_qz = 256;
    int   max_slices = 1024;

    // --- Display ---
    bool  log_scale = true;

    // --- Instrument ---
    float res_fwhm_qpar = 0.0f;       // Resolution (FWHM, nm^-1)
    float res_fwhm_qz = 0.0f;
    bool  view_detector = false;      // Show the result as a detector image
    Detector det;
    uint64_t det_hash = 0;
    size_t det_rows = 0, det_cols = 0;
    md_array(double) det_map = nullptr;     // [row 0 = top], display values (log10 if log_scale)
    md_array(float)  det_raw = nullptr;     // [row 0 = bottom], intensity, < 0 masked
    double det_qy_min = 0, det_qy_max = 1, det_qz_min = 0, det_qz_max = 1;   // nm^-1, approximate axis bounds
    uint64_t res_version = 0;
    float decades = 6.0f;
    ImPlotColormap colormap = ImPlotColormap_Viridis;
    bool  show_profile = false;

    // --- Cuts (positions in nm^-1) ---
    bool   show_cuts = true;
    double cut_qz = 0.0;          // Horizontal cut, I(q_par) at fixed q_z
    double cut_qpar = 0.0;        // Vertical cut, I(q_z) at fixed q_par
    float  cut_width_qz = 0.0f;   // Integration band width (nm^-1), 0 -> single row
    float  cut_width_qpar = 0.0f; // Integration band width (nm^-1), 0 -> single ring
    int    cut_mode = CutMode_Yoneda;
    float  cut_alpha_f_deg = 0.117f;  // Used when cut_mode == CutMode_AlphaF
    int    beam_preset = -1;          // Last applied preset, -1 = none
    bool   cut_log_x = false;
    bool   vcut_vs_alpha_f = false;   // Plot the vertical cut against alpha_f instead of q_z
    bool   cut_init = false;
    size_t hcut_rows = 0;
    size_t vcut_cols = 0;
    double link_qpar_min = 0.0, link_qpar_max = 1.0;
    double link_qz_min = 0.0, link_qz_max = 1.0;

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
    bool structure_stale = false;
    double ctx_q_z_max = 0.0;                // Å^-1, as used when ctx was created
    double particle_z_min = 0.0;
    char status[256] = "";
    double compute_time_start = 0.0;
    double compute_time = 0.0;

    // --- Model stage ---
    task_system::ID task_eval = task_system::INVALID_ID;
    std::atomic<bool> eval_ready{false};
    uint64_t eval_hash_pending = 0;          // Hash of the model currently being evaluated
    uint64_t eval_hash_shown = 0;            // Hash of the model currently shown
    md_gisaxs_model_t eval_model = {};
    md_array(double) eval_qz = nullptr;
    md_array(float)  eval_out = nullptr;     // num_qz * num_rings, written by the evaluation task

    // Shown result
    size_t res_rows = 0;
    size_t res_cols = 0;
    md_array(float)  res_raw = nullptr;      // [row = qz ascending][ring]
    md_array(double) res_map = nullptr;      // [row 0 = highest qz][ring], log10 if log_scale
    md_array(float)  res_smooth = nullptr;   // res_raw with the instrument resolution applied, same layout
    double res_min = 0.0, res_max = 1.0;
    double res_qpar_min = 0.0, res_qpar_max = 1.0;   // nm^-1 bounds
    double res_qz_min = 0.0, res_qz_max = 1.0;       // nm^-1 bounds
    double res_dq = 0.0;                             // Å^-1
    double res_dqz = 0.0;                            // Å^-1
    bool   res_dirty_map = false;
    md_array(double) prof_z = nullptr;
    md_array(double) prof_rho = nullptr;
    md_array(double) res_ring_q_nm = nullptr;   // Ring centers (nm^-1)
    md_array(unsigned) res_ring_count = nullptr;
    md_array(double) res_qz_nm = nullptr;       // Row q_z (nm^-1), ascending
    md_array(double) hcut = nullptr;
    md_array(double) vcut_af = nullptr;         // alpha_f (deg) per row, for the vertical cut
    md_array(double) vcut = nullptr;

    md_allocator_i* alloc = nullptr;

    Gisaxs() { viamd::event_system_register_handler(*this); }

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
                free_results();
                break;
            case viamd::EventType_ViamdFrameTick:
                update();
                draw_window();
                break;
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("GISAXS", &show_window);
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
                if (str_eq(viamd::section_header(state), STR_LIT("GISAXS"))) {
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

    double wavelength() const { return HC_KEV_ANGSTROM / MAX(energy_kev, 1.0e-3f); }

    // X-ray SLD (1/Å^2) from electron density (e/Å^3)
    static double sld(double rho_e) { return MD_GISAXS_R_E * rho_e; }

    // Critical angle (rad) of the substrate relative to the ambient
    double critical_angle() const {
        const double lambda = wavelength();
        const double d_sld = sld(sub_density) - sld(bg_density);
        if (d_sld <= 0.0) return 0.0;
        const double k0 = 2.0 * 3.14159265358979323846 / lambda;
        return asin(MIN(1.0, sqrt(4.0 * 3.14159265358979323846 * d_sld) / k0));
    }

    double contrast_factor() const {
        if (material_density <= 0.0f) return 1.0;
        const double f = 1.0 - bg_density / (double)material_density;
        return f * f;
    }

    double substrate_z() const {
        return sub_auto_z ? particle_z_min + sub_z_offset : sub_z;
    }

    uint64_t hash_structure_inputs() const {
        uint64_t h = md_hash64(filter, strnlen(filter, sizeof(filter)), 0x6153u);
        h = md_hash64(&electron_mode, sizeof(electron_mode), h);
        h = md_hash64(&electrons_per_particle, sizeof(electrons_per_particle), h);
        h = md_hash64(&sigma_mode, sizeof(sigma_mode), h);
        h = md_hash64(&sigma_constant, sizeof(sigma_constant), h);
        h = md_hash64(&q_par_max_nm, sizeof(q_par_max_nm), h);
        h = md_hash64(&q_z_max_nm, sizeof(q_z_max_nm), h);
        h = md_hash64(&oversampling, sizeof(oversampling), h);
        h = md_hash64(&max_slices, sizeof(max_slices), h);
        h = md_hash64(&electrons_per_dalton, sizeof(electrons_per_dalton), h);
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

    md_gisaxs_model_t build_model() const {
        md_gisaxs_model_t m = {};
        const double lambda = wavelength();
        m.wavelength = lambda;
        m.alpha_i = alpha_i_deg * DEG_TO_RAD;
        m.dwba = dwba;
        m.z_substrate = substrate_z();
        m.sld_ambient = sld(bg_density);
        m.sld_substrate = sld(sub_density);
        m.sld_substrate_abs = 2.0 * 3.14159265358979323846 * sub_beta / (lambda * lambda);
        m.substrate_roughness = sub_roughness;
        // Graded film: the laterally averaged particle density, contrast corrected like the particles themselves
        m.graded = film_graded;
        const double contrast_amp = material_density > 0.0f ? 1.0 - bg_density / (double)material_density : 1.0;
        m.profile_sld_scale = MD_GISAXS_R_E * contrast_amp;
        m.profile_abs_scale = material_density > 0.0f ? 2.0 * 3.14159265358979323846 * material_beta / (lambda * lambda * material_density) : 0.0;
        // dsigma/dOmega per unit area (sr^-1): r_e^2 * <|F|^2> / A, with the displaced background as contrast
        m.intensity_scale = MD_GISAXS_R_E * MD_GISAXS_R_E * contrast_factor();
        return m;
    }

    uint64_t hash_model(const md_gisaxs_model_t& m) const {
        const double vals[] = { m.wavelength, m.alpha_i, m.dwba ? 1.0 : 0.0, m.graded ? 1.0 : 0.0, m.z_substrate, m.sld_ambient,
                                m.sld_substrate, m.sld_substrate_abs, m.substrate_roughness, m.profile_sld_scale,
                                m.profile_abs_scale, m.intensity_scale };
        uint64_t h = md_hash64(vals, sizeof(vals), 0x4d6fu);
        h = md_hash64(&num_qz, sizeof(num_qz), h);
        h = md_hash64(&ctx, sizeof(ctx), h);
        h = md_hash64(&structure_hash, sizeof(structure_hash), h);
        return h;
    }

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
        const double lambda_nm = p.wavelength_nm;
        const double lambda_a = lambda_nm * 10.0;
        const double k0 = 2.0 * 3.14159265358979323846 / lambda_nm;   // nm^-1
        // Largest in-plane angle: half the detector width (beam centered horizontally)
        const double two_theta_max = atan(0.5 * p.det_px_h * p.pixel_mm / p.sdd_mm);
        const double ai = p.alpha_i_deg * DEG_TO_RAD;
        const double af = p.alpha_f_max_deg * DEG_TO_RAD;
        PresetValues v;
        v.energy_kev   = (float)(HC_KEV_ANGSTROM / lambda_a);
        v.alpha_i_deg  = (float)p.alpha_i_deg;
        v.q_par_max_nm = (float)(k0 * sin(two_theta_max));
        v.q_z_max_nm   = (float)(k0 * (sin(af) + sin(ai)));
        v.sub_density  = 2.0 * 3.14159265358979323846 * p.sub_delta / (lambda_a * lambda_a * MD_GISAXS_R_E);
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
        vcut_vs_alpha_f = true;     // As in the paper (Fig. 4c)
        film_graded     = true;
        material_beta   = p.material_beta;
        det.sdd_mm      = p.sdd_mm;
        det.pixel_mm    = p.pixel_mm;
        det.npx_h       = p.det_px_h;
        det.npx_v       = p.det_px_v;
        det.beam_x_px   = 0.5 * p.det_px_h;
        det.beam_y_px   = p.beam_y_px;
        det.module_w    = p.module_w;
        det.module_h    = p.module_h;
        det.gap_w       = p.gap_w;
        det.gap_h       = p.gap_h;
        det.gaps        = p.module_w > 0;
        det.bs_direct_px   = p.bs_direct_px;
        det.bs_specular_px = p.bs_specular_px;
        // Resolution: one detector pixel (the angular beam divergence at P03 is of the same order or smaller)
        {
            const double k = 2.0 * 3.14159265358979323846 / p.wavelength_nm;
            const float dq = (float)(k * p.pixel_mm / p.sdd_mm);
            res_fwhm_qpar = dq;
            res_fwhm_qz = dq;
        }
        beam_preset  = idx;
    }

    // ---------------------------------------------------------------------------------------------
    // Lifetime helpers
    // ---------------------------------------------------------------------------------------------

    void cancel_all() {
        compute_cancel = true;
        task_system::task_interrupt_and_wait_for(task_slices);
        task_system::task_interrupt_and_wait_for(task_rings);
        task_system::task_wait_for(task_finish);
        task_system::task_wait_for(task_eval);
        task_slices = task_rings = task_finish = task_eval = task_system::INVALID_ID;
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
        res_rows = res_cols = 0;
        md_array_shrink(res_raw, 0);
        md_array_shrink(res_map, 0);
        md_array_shrink(res_smooth, 0);
        md_array_shrink(det_map, 0);
        md_array_shrink(det_raw, 0);
        det_rows = det_cols = 0;
        md_array_shrink(prof_z, 0);
        md_array_shrink(prof_rho, 0);
        eval_hash_shown = 0;
        eval_hash_pending = 0;
    }

    void free_results() {
        md_array_free(res_raw, alloc);
        md_array_free(res_map, alloc);
        md_array_free(res_smooth, alloc);
        md_array_free(det_map, alloc);
        md_array_free(det_raw, alloc);
        res_smooth = nullptr; det_map = nullptr; det_raw = nullptr;
        md_array_free(eval_qz, alloc);
        md_array_free(eval_out, alloc);
        md_array_free(prof_z, alloc);
        md_array_free(prof_rho, alloc);
        md_array_free(res_ring_q_nm, alloc);
        md_array_free(res_ring_count, alloc);
        md_array_free(res_qz_nm, alloc);
        md_array_free(hcut, alloc);
        md_array_free(vcut_af, alloc);
        vcut_af = nullptr;
        md_array_free(vcut, alloc);
        res_ring_q_nm = nullptr; res_ring_count = nullptr; res_qz_nm = nullptr; hcut = nullptr; vcut = nullptr;
        res_raw = nullptr; res_map = nullptr; eval_qz = nullptr; eval_out = nullptr; prof_z = nullptr; prof_rho = nullptr;
    }

    // ---------------------------------------------------------------------------------------------
    // Structure stage
    // ---------------------------------------------------------------------------------------------

    // Detects slices (groups of fib_stride consecutive beads containing one center bead) and chains of consecutive
    // slices, and sweeps the cross-section template along them.
    // x, y, z, w are the gathered (selected) particles, atom_idx their atom indices (ascending).
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
                    MD_LOG_INFO("GISAXS fibrils: template has %zu anchors, expected %i, orienting with the mean MD layout instead", tmpl.anchor_u.size(), S);
                }
                tmpl.anchor_u = mu;
                tmpl.anchor_v = mv;
            } else {
                double rms = 0.0;
                const bool mirrored = fibril_template_match_handedness(tmpl, mu, mv, &rms);
                MD_LOG_INFO("GISAXS fibrils: template anchors match the mean MD slice layout to %.2f Å RMS%s", rms, mirrored ? " (template mirrored)" : "");
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
        MD_LOG_INFO("GISAXS fibrils: %s (electrons %.6g -> %.6g, RMS axial bead offset %.2f Å)", fib_info, st.electrons_in, st.electrons_out, st.mean_axial_offset);
        return true;
    }

    bool start_compute() {
        cancel_all();
        destroy_ctx();

        const md_system_t& sys = app_state->mold.sys;
        const md_system_state_t& state = app_state->mold.state;
        const size_t num_atoms = sys.atom.count;
        if (num_atoms == 0 || !state.x) {
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
        const size_t count = md_bitfield_popcount(&mask);
        if (count == 0) {
            snprintf(status, sizeof(status), "Selection is empty");
            return false;
        }

        // Gather particles (temporary copies, md_gisaxs_create keeps its own sorted copy)
        md_array(float) x = nullptr; md_array(float) y = nullptr; md_array(float) z = nullptr;
        md_array(float) w = nullptr; md_array(float) s = nullptr;
        std::vector<uint32_t> atom_idx;
        atom_idx.reserve(count);
        defer {
            md_array_free(x, alloc); md_array_free(y, alloc); md_array_free(z, alloc);
            md_array_free(w, alloc); md_array_free(s, alloc);
        };
        md_array_resize(x, count, alloc);
        md_array_resize(y, count, alloc);
        md_array_resize(z, count, alloc);
        md_array_resize(w, count, alloc);
        md_array_resize(s, count, alloc);

        size_t n = 0;
        size_t num_zero_weight = 0;
        md_bitfield_iter_t it = md_bitfield_iter_create(&mask);
        while (md_bitfield_iter_next(&it)) {
            const size_t idx = md_bitfield_iter_idx(&it);
            if (idx >= num_atoms) continue;
            x[n] = state.x[idx];
            y[n] = state.y[idx];
            z[n] = state.z[idx];
            atom_idx.push_back((uint32_t)idx);
            if (electron_mode == ElectronMode_AtomicNumber) {
                w[n] = (float)md_atom_atomic_number(&sys.atom, idx);
            } else if (electron_mode == ElectronMode_Mass) {
                w[n] = md_atom_mass(&sys.atom, idx) * electrons_per_dalton;
            } else {
                w[n] = electrons_per_particle;
            }
            if (w[n] == 0.0f) num_zero_weight += 1;
            if (sigma_mode == SigmaMode_Radius) {
                s[n] = md_atom_radius(&sys.atom, idx) / sqrtf(5.0f);
            } else {
                s[n] = sigma_constant;
            }
            ++n;
        }
        if (num_zero_weight == n) {
            snprintf(status, sizeof(status), "All particles have zero electrons (atomic numbers or masses missing?), use a constant electron count");
            return false;
        }

        // Self scattering of the beads (diagnostic)
        bead_self.clear();
        bead_self_area = cell.x * cell.y;
        for (size_t i = 0; i < n; ++i) {
            const float sig = s[i];
            size_t g = 0;
            for (; g < bead_self.size(); ++g) if (bead_self[g].first == sig) break;
            if (g == bead_self.size()) {
                if (bead_self.size() >= 16) g = bead_self.size() - 1;   // Enough for a diagnostic
                else bead_self.push_back({sig, 0.0});
            }
            bead_self[g].second += (double)w[i] * w[i];
        }

        // Fibril model: replace the beads by a continuous swept cross-section
        FibrilOutput fib;
        fib_info[0] = '\0';
        if (density_model == DensityModel_Fibril) {
            if (!build_fibrils(&fib, atom_idx, x, y, z, w, cell)) {
                return false;
            }
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
            input.x = x;
            input.y = y;
            input.z = z;
            input.weight = w;
            input.sigma = (sigma_mode == SigmaMode_Constant) ? nullptr : s;
            input.sigma_uniform = sigma_constant;
        }
        input.box_x = cell.x;
        input.box_y = cell.y;

        md_gisaxs_params_t params = {};
        params.q_par_max = q_par_max_nm * 0.1;
        params.q_z_max = q_z_max_nm * 0.1;
        params.oversampling = oversampling;
        params.max_slices = (size_t)MAX(max_slices, 8);

        ctx = md_gisaxs_create(&input, &params, alloc);
        if (!ctx) {
            snprintf(status, sizeof(status), "Failed to initialize GISAXS computation (see log)");
            return false;
        }
        ctx_q_z_max = params.q_z_max;
        particle_z_min = md_gisaxs_particle_z_min(ctx);
        structure_hash = hash_structure_inputs();
        structure_stale = false;

        md_gisaxs_info_t info;
        md_gisaxs_get_info(ctx, &info);
        MD_LOG_INFO("GISAXS: %zu particles, grid %i x %i (dx %.2f Å), %zu slices (dz %.2f Å), %zu rings, %zu classes, %.1f MB spectra, %.1f MB matrices",
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

        const uint32_t num_slices = (uint32_t)info.num_slices;
        const uint32_t num_rings  = (uint32_t)info.num_rings;

        task_slices = task_system::create_pool_task(STR_LIT("GISAXS slices"), num_slices, [this](uint32_t beg, uint32_t end, uint32_t thread_num) {
            if (compute_cancel) return;
            if (thread_num >= md_array_size(scratch)) return;   // Should not happen
            if (!scratch[thread_num]) {
                scratch[thread_num] = md_fft_alloc(scratch_bytes);
                if (!scratch[thread_num]) { compute_cancel = true; return; }
            }
            md_gisaxs_compute_slices(ctx, beg, end, scratch[thread_num]);
        });

        task_rings = task_system::create_pool_task(STR_LIT("GISAXS rings"), num_rings, [this](uint32_t beg, uint32_t end, uint32_t) {
            if (compute_cancel) return;
            md_gisaxs_compute_rings(ctx, beg, end);
        });

        task_finish = task_system::create_pool_task(STR_LIT("GISAXS finalize"), [this]() {
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
        md_array_resize(eval_qz, Nq, alloc);
        md_array_resize(eval_out, Nq * R, alloc);
        for (size_t i = 0; i < Nq; ++i) {
            eval_qz[i] = ctx_q_z_max * (double)i / (double)(Nq - 1);
        }
        eval_model = model;
        eval_hash_pending = hash;
        eval_ready = false;

        task_eval = task_system::create_pool_task(STR_LIT("GISAXS evaluate"), (uint32_t)Nq, [this](uint32_t beg, uint32_t end, uint32_t) {
            md_gisaxs_evaluate_range(ctx, &eval_model, eval_qz, beg, end, eval_out);
        }, 4);
        const task_system::ID done = task_system::create_pool_task(STR_LIT("GISAXS evaluate done"), [this]() {
            eval_ready = true;
        });
        task_system::set_task_dependency(done, task_eval);
        task_system::enqueue_task(task_eval);
        task_eval = done;
    }

    void accept_eval() {
        const size_t R = md_gisaxs_num_rings(ctx);
        const size_t Nq = md_array_size(eval_qz);
        res_rows = Nq;
        res_cols = R;
        md_array_resize(res_raw, Nq * R, alloc);
        MEMCPY(res_raw, eval_out, sizeof(float) * Nq * R);

        const double dq = R > 0 ? md_gisaxs_ring_q(ctx)[0] : 0.0;
        md_gisaxs_info_t info;
        md_gisaxs_get_info(ctx, &info);
        res_dq = info.dq_ring;
        res_dqz = Nq > 1 ? eval_qz[1] - eval_qz[0] : 0.0;
        (void)dq;
        // Ring r is centered at (r + 1) * dq
        const double prev_bounds[4] = {res_qpar_min, res_qpar_max, res_qz_min, res_qz_max};
        res_qpar_min = 0.5 * res_dq * 10.0;
        res_qpar_max = (R + 0.5) * res_dq * 10.0;
        res_qz_min = (eval_qz[0] - 0.5 * res_dqz) * 10.0;
        res_qz_max = (eval_qz[Nq - 1] + 0.5 * res_dqz) * 10.0;

        // The map and the cut plots share linked axes. Reset the view when the q-range changed.
        if (prev_bounds[0] != res_qpar_min || prev_bounds[1] != res_qpar_max || prev_bounds[2] != res_qz_min || prev_bounds[3] != res_qz_max) {
            link_qpar_min = res_qpar_min;
            link_qpar_max = res_qpar_max;
            link_qz_min = res_qz_min;
            link_qz_max = res_qz_max;
        }

        const size_t S = md_gisaxs_num_slices(ctx);
        md_array_resize(prof_z, S, alloc);
        md_array_resize(prof_rho, S, alloc);
        const double* sz = md_gisaxs_slice_z(ctx);
        const double* sp = md_gisaxs_slice_profile(ctx);
        for (size_t i = 0; i < S; ++i) {
            prof_z[i] = sz[i];
            prof_rho[i] = sp[i];
        }

        md_array_resize(res_ring_q_nm, R, alloc);
        md_array_resize(res_ring_count, R, alloc);
        const double* rq = md_gisaxs_ring_q(ctx);
        const unsigned* rc = md_gisaxs_ring_count(ctx);
        for (size_t i = 0; i < R; ++i) {
            res_ring_q_nm[i] = rq[i] * 10.0;
            res_ring_count[i] = rc[i];
        }
        md_array_resize(res_qz_nm, Nq, alloc);
        for (size_t i = 0; i < Nq; ++i) res_qz_nm[i] = eval_qz[i] * 10.0;

        eval_hash_shown = eval_hash_pending;
        res_dirty_map = true;
        eval_ready = false;
        task_eval = task_system::INVALID_ID;
    }

    void rebuild_map() {
        const size_t N = res_rows * res_cols;
        md_array_resize(res_map, N, alloc);
        double vmax = -DBL_MAX;
        double vmin = DBL_MAX;
        for (size_t i = 0; i < N; ++i) {
            if (res_smooth[i] > 0.0f) {
                vmax = MAX(vmax, (double)res_smooth[i]);
                vmin = MIN(vmin, (double)res_smooth[i]);
            }
        }
        if (vmax <= 0.0 || vmax == -DBL_MAX) {
            vmax = 1.0;
            vmin = 0.0;
        }
        if (log_scale) {
            res_max = log10(vmax);
            res_min = res_max - MAX(decades, 0.5f);
        } else {
            res_max = vmax;
            res_min = 0.0;
        }
        for (size_t row = 0; row < res_rows; ++row) {
            // Row 0 of the heat map is drawn at the top, i.e. highest q_z
            const float* src = res_smooth + (res_rows - 1 - row) * res_cols;
            double* dst = res_map + row * res_cols;
            for (size_t c = 0; c < res_cols; ++c) {
                double v = src[c];
                if (log_scale) {
                    v = v > 0.0 ? log10(v) : res_min;
                    v = MAX(v, res_min);
                }
                dst[c] = v;
            }
        }
        res_dirty_map = false;
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

        if (compute_state.load() == ComputeState_Done && ctx) {
            if (!task_system::task_is_running(task_eval) && !eval_ready.load()) {
                const md_gisaxs_model_t model = build_model();
                const uint64_t h = hash_model(model);
                if (h != eval_hash_shown) {
                    start_eval(model, h);
                }
            }
        }

        if (res_dirty_map && res_rows && res_cols) {
            apply_resolution();
            rebuild_map();
            res_version += 1;
        }

        if (view_detector && res_rows && res_cols) {
            const uint64_t h = hash_detector();
            if (h != det_hash) {
                compute_detector();
                det_hash = h;
            }
        }
    }

    // ---------------------------------------------------------------------------------------------
    // Instrument resolution
    // ---------------------------------------------------------------------------------------------

    // Gaussian smoothing of res_raw into res_smooth along q_par (rings) and q_z (rows).
    // Rows below the sample horizon (DWBA) carry no intensity and are excluded from the q_z smoothing.
    void apply_resolution() {
        const size_t R = res_cols, Q = res_rows;
        md_array_resize(res_smooth, R * Q, alloc);
        MEMCPY(res_smooth, res_raw, sizeof(float) * R * Q);

        const double fwhm_to_sigma = 1.0 / 2.354820045;
        const double sc = res_dq > 0.0 ? res_fwhm_qpar * fwhm_to_sigma / (res_dq * 10.0) : 0.0;
        const double sr = res_dqz > 0.0 ? res_fwhm_qz * fwhm_to_sigma / (res_dqz * 10.0) : 0.0;

        md_array(float) tmp = nullptr;
        md_array(double) kern = nullptr;
        defer { md_array_free(tmp, alloc); md_array_free(kern, alloc); };

        auto make_kernel = [&](double sigma) -> int {
            const int rad = (int)ceil(3.0 * sigma);
            md_array_resize(kern, (size_t)(2 * rad + 1), alloc);
            for (int i = -rad; i <= rad; ++i) kern[i + rad] = exp(-0.5 * (i * i) / (sigma * sigma));
            return rad;
        };

        if (sc > 0.05) {
            const int rad = make_kernel(sc);
            md_array_resize(tmp, R, alloc);
            for (size_t r = 0; r < Q; ++r) {
                float* row = res_smooth + r * R;
                for (size_t c = 0; c < R; ++c) {
                    double sum = 0.0, wsum = 0.0;
                    for (int i = -rad; i <= rad; ++i) {
                        const long cc = (long)c + i;
                        if (cc < 0 || cc >= (long)R) continue;
                        sum += kern[i + rad] * row[cc];
                        wsum += kern[i + rad];
                    }
                    tmp[c] = (float)(wsum > 0.0 ? sum / wsum : row[c]);
                }
                MEMCPY(row, tmp, sizeof(float) * R);
            }
        }

        if (sr > 0.05) {
            const int rad = make_kernel(sr);
            md_array_resize(tmp, Q, alloc);
            const double horizon = eval_model.dwba ? horizon_qz_nm() : -DBL_MAX;
            for (size_t c = 0; c < R; ++c) {
                for (size_t r = 0; r < Q; ++r) {
                    if (res_qz_nm[r] < horizon) { tmp[r] = 0.0f; continue; }
                    double sum = 0.0, wsum = 0.0;
                    for (int i = -rad; i <= rad; ++i) {
                        const long rr = (long)r + i;
                        if (rr < 0 || rr >= (long)Q || res_qz_nm[rr] < horizon) continue;
                        sum += kern[i + rad] * res_smooth[rr * R + c];
                        wsum += kern[i + rad];
                    }
                    tmp[r] = (float)(wsum > 0.0 ? sum / wsum : 0.0);
                }
                for (size_t r = 0; r < Q; ++r) res_smooth[r * R + c] = tmp[r];
            }
        }
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

    bool in_gap(int px, int py) const {
        if (!det.gaps || det.module_w <= 0 || det.module_h <= 0) return false;
        const int pw = det.module_w + det.gap_w;
        const int ph = det.module_h + det.gap_h;
        return (px % pw) >= det.module_w || (py % ph) >= det.module_h;
    }

    // Bilinear lookup of the (smoothed) q-map at q_par, q_z (nm^-1). Returns false outside the data.
    bool lookup(double qpar_nm, double qz_nm, float* out) const {
        const double cf = qpar_nm / (res_dq * 10.0) - 1.0;     // ring r is centered at (r + 1) dq
        const double rf = qz_nm / (res_dqz * 10.0);            // row r at r * dqz
        if (cf < -0.5 || cf > (double)res_cols - 1.0) return false;   // specular rod / beyond q_par max
        if (rf < 0.0 || rf > (double)res_rows - 1.0) return false;
        const double c = MAX(cf, 0.0);
        const size_t c0 = MIN((size_t)c, res_cols - 1);
        const size_t c1 = MIN(c0 + 1, res_cols - 1);
        const size_t r0 = MIN((size_t)rf, res_rows - 1);
        const size_t r1 = MIN(r0 + 1, res_rows - 1);
        const double tc = c - c0, tr = rf - r0;
        const float* m = res_smooth;
        const double v = (1 - tr) * ((1 - tc) * m[r0 * res_cols + c0] + tc * m[r0 * res_cols + c1]) +
                         tr       * ((1 - tc) * m[r1 * res_cols + c0] + tc * m[r1 * res_cols + c1]);
        *out = (float)v;
        return true;
    }

    // Scattering geometry of a detector position (mm, relative to the direct beam, y upwards).
    // Returns false below the sample horizon.
    bool pixel_q(double x_mm, double y_mm, double* qpar_nm, double* qy_nm, double* qz_nm, double* af_deg) const {
        const double lambda_nm = eval_model.wavelength * 0.1;
        const double k = 2.0 * 3.14159265358979323846 / lambda_nm;
        const double ai = eval_model.alpha_i;
        // Direct beam direction and detector axes (sample frame, z normal to the surface)
        const double kin[3] = {cos(ai), 0.0, -sin(ai)};
        const double eup[3] = {sin(ai), 0.0, cos(ai)};
        double P[3] = {
            det.sdd_mm * kin[0] + y_mm * eup[0],
            x_mm,
            det.sdd_mm * kin[2] + y_mm * eup[2],
        };
        const double len = sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2]);
        P[0] /= len; P[1] /= len; P[2] /= len;
        const double qx = k * (P[0] - kin[0]);
        const double qy = k * P[1];
        const double qz = k * (P[2] - kin[2]);
        *qpar_nm = sqrt(qx * qx + qy * qy);
        *qy_nm = qy;
        *qz_nm = qz;
        *af_deg = asin(CLAMP(P[2], -1.0, 1.0)) * RAD_TO_DEG;
        return P[2] >= 0.0;
    }

    void compute_detector() {
        const int b = MAX(det.binning, 1);
        const int W = MAX(det.npx_h / b, 1);
        const int H = MAX(det.npx_v / b, 1);
        det_cols = (size_t)W;
        det_rows = (size_t)H;
        md_array_resize(det_raw, det_cols * det_rows, alloc);
        md_array_resize(det_map, det_cols * det_rows, alloc);

        const double ai = eval_model.alpha_i;
        const double spec_dy_px = tan(2.0 * ai) * det.sdd_mm / det.pixel_mm;   // specular spot above the direct beam
        const double bd = 0.5 * det.bs_direct_px;
        const double bs = 0.5 * det.bs_specular_px;

        for (int by = 0; by < H; ++by) {
            for (int bx = 0; bx < W; ++bx) {
                double sum = 0.0;
                int n = 0;
                for (int sy = 0; sy < b; ++sy) {
                    for (int sx = 0; sx < b; ++sx) {
                        const int px = bx * b + sx;
                        const int py = by * b + sy;
                        if (in_gap(px, py)) continue;
                        const double dxp = px + 0.5 - det.beam_x_px;
                        const double dyp = py + 0.5 - det.beam_y_px;
                        if (fabs(dxp) < bd && fabs(dyp) < bd) continue;                       // direct beamstop
                        if (fabs(dxp) < bs && fabs(dyp - spec_dy_px) < bs) continue;          // specular beamstop
                        double qpar, qy, qz, af;
                        if (!pixel_q(dxp * det.pixel_mm, dyp * det.pixel_mm, &qpar, &qy, &qz, &af)) continue;   // sample shadow
                        float v;
                        if (!lookup(qpar, qz, &v)) continue;
                        sum += v;
                        n += 1;
                    }
                }
                det_raw[(size_t)by * W + bx] = n ? (float)(sum / n) : -1.0f;
            }
        }

        // Display values, row 0 at the top
        for (int by = 0; by < H; ++by) {
            const float* src = det_raw + (size_t)by * W;
            double* dst = det_map + (size_t)(H - 1 - by) * W;
            for (int bx = 0; bx < W; ++bx) {
                const double v = src[bx];
                if (v < 0.0) {
                    dst[bx] = res_min - 1.0e3 * (res_max - res_min + 1.0);   // masked: below the color range
                } else if (log_scale) {
                    dst[bx] = v > 0.0 ? MAX(log10(v), res_min) : res_min;
                } else {
                    dst[bx] = v;
                }
            }
        }

        // Approximate (linear) axis bounds along the central row / column, used for display only
        double qpar, qy, qz, af;
        const double x0 = (0.0 - det.beam_x_px) * det.pixel_mm;
        const double x1 = (W * b - det.beam_x_px) * det.pixel_mm;
        const double y0 = (0.0 - det.beam_y_px) * det.pixel_mm;
        const double y1 = (H * b - det.beam_y_px) * det.pixel_mm;
        pixel_q(x0, 0.0, &qpar, &qy, &qz, &af); det_qy_min = qy;
        pixel_q(x1, 0.0, &qpar, &qy, &qz, &af); det_qy_max = qy;
        pixel_q(0.0, y0, &qpar, &qy, &qz, &af); det_qz_min = qz;
        pixel_q(0.0, y1, &qpar, &qy, &qz, &af); det_qz_max = qz;
    }

    // ---------------------------------------------------------------------------------------------
    // UI
    // ---------------------------------------------------------------------------------------------

    void draw_window() {
        if (!show_window) return;

        ImGui::SetNextWindowSize({900, 600}, ImGuiCond_FirstUseEver);
        if (!ImGui::Begin("GISAXS", &show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            ImGui::End();
            return;
        }

        const float settings_w = 330.0f;
        ImGui::BeginChild("##gisaxs_settings", ImVec2(settings_w, 0), ImGuiChildFlags_Borders);
        draw_settings();
        ImGui::EndChild();

        ImGui::SameLine();
        ImGui::BeginChild("##gisaxs_plot", ImVec2(0, 0));
        draw_plot();
        ImGui::EndChild();

        ImGui::End();
    }

    void draw_settings() {
        const bool running = compute_state.load() == ComputeState_Running;
        const float item_w = 150.0f;

        // --- Compute ---
        ImGui::BeginDisabled(running);
        if (ImGui::Button("Compute", ImVec2(-1, 0))) {
            start_compute();
        }
        ImGui::EndDisabled();
        if (running) {
            float frac = 0.0f;
            if (task_system::task_is_running(task_slices)) {
                frac = 0.8f * task_system::task_fraction_complete(task_slices);
            } else if (task_system::task_is_running(task_rings)) {
                frac = 0.8f + 0.2f * task_system::task_fraction_complete(task_rings);
            } else {
                frac = 1.0f;
            }
            ImGui::ProgressBar(frac, ImVec2(-1, 0));
            if (ImGui::Button("Cancel", ImVec2(-1, 0))) {
                compute_cancel = true;
                task_system::task_interrupt(task_slices);
                task_system::task_interrupt(task_rings);
            }
        }
        if (status[0]) {
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextDisabled("%s", status);
            ImGui::PopTextWrapPos();
        }
        const bool inputs_changed = ctx && (hash_structure_inputs() != structure_hash);
        if (!running && (structure_stale || inputs_changed)) {
            ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.3f, 1.0f), "Structure or grid changed, recompute");
        }

        ImGui::PushItemWidth(item_w);

        // --- Particles ---
        if (ImGui::CollapsingHeader("Particles", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::InputQuery("Selection", filter, sizeof(filter), filter_valid, filter_err);
            if (ImGui::BeginCombo("Electrons", electron_mode_lbl[electron_mode])) {
                for (int i = 0; i < (int)ARRAY_SIZE(electron_mode_lbl); ++i) {
                    if (ImGui::Selectable(electron_mode_lbl[i], electron_mode == i)) electron_mode = i;
                }
                ImGui::EndCombo();
            }
            if (electron_mode == ElectronMode_Constant) {
                ImGui::InputFloat("e / particle", &electrons_per_particle, 0, 0, "%.1f");
            } else if (electron_mode == ElectronMode_Mass) {
                ImGui::InputFloat("e / Da", &electrons_per_dalton, 0, 0, "%.4f");
                ImGui::SetItemTooltip("Electrons per dalton of mass. Cellulose 0.530, water 0.555, pure C/N/O 0.5");
            }
            if (ImGui::BeginCombo("Gaussian width", sigma_mode_lbl[sigma_mode])) {
                for (int i = 0; i < 2; ++i) {
                    if (ImGui::Selectable(sigma_mode_lbl[i], sigma_mode == i)) sigma_mode = i;
                }
                ImGui::EndCombo();
            }
            if (sigma_mode == SigmaMode_Constant) {
                ImGui::InputFloat("Sigma (Å)", &sigma_constant, 0, 0, "%.2f");
                sigma_constant = MAX(sigma_constant, 0.0f);
            }
            ImGui::InputFloat("Material density (e/Å³)", &material_density, 0, 0, "%.4f");
            ImGui::SetItemTooltip("Electron density of the particle material. Each particle displaces\n"
                                  "a volume (electrons / material density) of the background medium.");
            ImGui::InputDouble("Material beta", &material_beta, 0, 0, "%.3e");

            ImGui::Separator();
            if (ImGui::BeginCombo("Density model", density_model_lbl[density_model])) {
                for (int i = 0; i < (int)ARRAY_SIZE(density_model_lbl); ++i) {
                    if (ImGui::Selectable(density_model_lbl[i], density_model == i)) density_model = i;
                }
                ImGui::EndCombo();
            }
            ImGui::SetItemTooltip("Beads: one isotropic Gaussian per particle.\n"
                                  "Fibrils: slices of beads (one center bead + surrounding beads) define a centerline and an\n"
                                  "orientation, and a cross-section template is swept continuously along it. Removes the\n"
                                  "artificial scattering of discrete beads (high-q plateau, peak at 2 pi / bead spacing).");
            if (density_model == DensityModel_Fibril) {
                ImGui::InputText("Center bead", fib_center_name, sizeof(fib_center_name));
                ImGui::InputInt("Beads per slice", &fib_stride);
                fib_stride = CLAMP(fib_stride, 2, 64);
                if (ImGui::BeginCombo("Cross-section", fibril_template_lbl[fib_template])) {
                    for (int i = 0; i < (int)ARRAY_SIZE(fibril_template_lbl); ++i) {
                        if (ImGui::Selectable(fibril_template_lbl[i], fib_template == i)) fib_template = i;
                    }
                    ImGui::EndCombo();
                }
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
            }
            ImGui::SetItemTooltip("Absorption (imaginary part of the refractive index) of the particle material at the\n"
                                  "material density. Only used for the film in the graded DWBA. Cellulose ~2e-9 at 12.8 keV.");
        }

        // --- Background ---
        if (ImGui::CollapsingHeader("Background medium", ImGuiTreeNodeFlags_DefaultOpen)) {
            if (ImGui::BeginCombo("Medium", background_presets[bg_preset].name)) {
                for (int i = 0; i < (int)ARRAY_SIZE(background_presets); ++i) {
                    if (ImGui::Selectable(background_presets[i].name, bg_preset == i)) {
                        bg_preset = i;
                        if (i != BG_CUSTOM) bg_density = background_presets[i].electron_density;
                    }
                }
                ImGui::EndCombo();
            }
            ImGui::BeginDisabled(bg_preset != BG_CUSTOM);
            ImGui::InputDouble("Density (e/Å³)", &bg_density, 0, 0, "%.5f");
            ImGui::EndDisabled();
            bg_density = MAX(bg_density, 0.0);
            ImGui::TextDisabled("Contrast factor: %.4f", contrast_factor());
        }

        // --- Substrate ---
        if (ImGui::CollapsingHeader("Substrate (DWBA)", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Checkbox("Enable substrate (DWBA)", &dwba);
            ImGui::SetItemTooltip("Unchecked: Born approximation without a substrate");
            ImGui::BeginDisabled(!dwba);
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
            ImGui::Checkbox("Place below lowest particle", &sub_auto_z);
            if (sub_auto_z) {
                ImGui::InputFloat("Offset (Å)", &sub_z_offset, 0, 0, "%.2f");
                ImGui::SetItemTooltip("Substrate z relative to the lowest selected particle (negative = below)");
                if (ctx) ImGui::TextDisabled("z substrate: %.2f Å", substrate_z());
            } else {
                ImGui::InputFloat("z (Å)", &sub_z, 0, 0, "%.2f");
            }
            ImGui::InputFloat("Roughness (Å)", &sub_roughness, 0, 0, "%.2f");
            sub_roughness = MAX(sub_roughness, 0.0f);
            ImGui::SetItemTooltip("RMS roughness of the substrate interface (Nevot-Croce factor)");
            ImGui::TextDisabled("Critical angle: %.3f°", critical_angle() * RAD_TO_DEG);
            ImGui::Checkbox("Film in reference medium (graded DWBA)", &film_graded);
            ImGui::SetItemTooltip("Include the laterally averaged density of the particles (the density profile)\n"
                                  "in the DWBA reference medium: vacuum / graded film / substrate.\n"
                                  "Gives refraction inside the film and the film Yoneda peak.");
            if (film_graded && md_array_size(prof_rho)) {
                const double d = film_sld_rel();
                const double k0 = 2.0 * 3.14159265358979323846 / wavelength();
                if (d > 0.0) ImGui::TextDisabled("Film critical angle: %.3f°", asin(MIN(1.0, sqrt(4.0 * 3.14159265358979323846 * d) / k0)) * RAD_TO_DEG);
            }
            ImGui::EndDisabled();
        }

        // --- Beam ---
        if (ImGui::CollapsingHeader("Beam", ImGuiTreeNodeFlags_DefaultOpen)) {
            const bool preset_valid = beam_preset >= 0 && beam_preset < (int)ARRAY_SIZE(beam_presets);
            char preview[128];
            if (preset_valid) {
                snprintf(preview, sizeof(preview), "%s%s", beam_presets[beam_preset].name, preset_matches(beam_presets[beam_preset]) ? "" : " (modified)");
            } else {
                snprintf(preview, sizeof(preview), "Custom");
            }
            if (ImGui::BeginCombo("Preset", preview)) {
                for (int i = 0; i < (int)ARRAY_SIZE(beam_presets); ++i) {
                    if (ImGui::Selectable(beam_presets[i].name, beam_preset == i)) {
                        apply_preset(i);
                    }
                    ImGui::SetItemTooltip("%s", beam_presets[i].description);
                }
                ImGui::EndCombo();
            }
            if (preset_valid) {
                ImGui::SetItemTooltip("%s", beam_presets[beam_preset].description);
            }
            ImGui::InputFloat("Energy (keV)", &energy_kev, 0, 0, "%.3f");
            energy_kev = CLAMP(energy_kev, 0.1f, 1000.0f);
            ImGui::TextDisabled("Wavelength: %.4f Å", wavelength());
            ImGui::SliderFloat("Incidence (°)", &alpha_i_deg, 0.01f, 2.0f, "%.3f");
        }

        // --- Grid ---
        if (ImGui::CollapsingHeader("q-range and sampling")) {
            ImGui::InputFloat("q_par max (nm⁻¹)", &q_par_max_nm, 0, 0, "%.3f");
            ImGui::InputFloat("q_z max (nm⁻¹)", &q_z_max_nm, 0, 0, "%.3f");
            q_par_max_nm = MAX(q_par_max_nm, 0.01f);
            q_z_max_nm = MAX(q_z_max_nm, 0.01f);
            ImGui::SliderInt("q_z rows", &num_qz, 16, 1024);
            ImGui::SliderFloat("Oversampling", &oversampling, 1.5f, 4.0f, "%.2f");
            ImGui::SetItemTooltip("Grid oversampling relative to q max (in-plane and z).\nHigher reduces aliasing at the cost of memory and time.");
            ImGui::InputInt("Max slices", &max_slices);
            max_slices = CLAMP(max_slices, 8, 8192);
            if (ctx) {
                md_gisaxs_info_t info;
                md_gisaxs_get_info(ctx, &info);
                ImGui::TextDisabled("Grid %i x %i, dx = %.2f Å", info.nx, info.ny, info.dx);
                ImGui::TextDisabled("%zu slices, dz = %.2f Å", info.num_slices, info.dz);
                ImGui::TextDisabled("%zu rings, dq = %.4f nm⁻¹", info.num_rings, info.dq_ring * 10.0);
                ImGui::TextDisabled("Ring matrices: %.1f MB", info.matrix_bytes / (1024.0 * 1024.0));
            }
        }

        // --- Instrument ---
        if (ImGui::CollapsingHeader("Instrument")) {
            bool changed = false;
            changed |= ImGui::InputFloat("Resolution q_par (FWHM, nm⁻¹)", &res_fwhm_qpar, 0, 0, "%.4f");
            changed |= ImGui::InputFloat("Resolution q_z (FWHM, nm⁻¹)", &res_fwhm_qz, 0, 0, "%.4f");
            res_fwhm_qpar = MAX(res_fwhm_qpar, 0.0f);
            res_fwhm_qz = MAX(res_fwhm_qz, 0.0f);
            if (changed) res_dirty_map = true;
            if (res_dq > 0.0) {
                ImGui::TextDisabled("Ring spacing (box limit): %.4f nm⁻¹", res_dq * 10.0);
            }
            ImGui::Checkbox("Show detector image", &view_detector);
            ImGui::BeginDisabled(!view_detector);
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
            ImGui::EndDisabled();
        }

        // --- Display ---
        if (ImGui::CollapsingHeader("Display")) {
            if (ImGui::Checkbox("Log scale", &log_scale)) res_dirty_map = true;
            if (log_scale) {
                if (ImGui::SliderFloat("Decades", &decades, 1.0f, 12.0f, "%.1f")) res_dirty_map = true;
            }
            if (ImGui::BeginCombo("Colormap", ImPlot::GetColormapName(colormap))) {
                for (int i = 0; i < ImPlot::GetColormapCount(); ++i) {
                    if (ImGui::Selectable(ImPlot::GetColormapName(i), colormap == i)) colormap = i;
                }
                ImGui::EndCombo();
            }
            ImGui::Checkbox("Show density profile", &show_profile);
        }

        // --- Cuts ---
        if (ImGui::CollapsingHeader("Cuts", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Checkbox("Show cuts", &show_cuts);
            ImGui::BeginDisabled(!show_cuts);
            if (ImGui::BeginCombo("Horizontal cut", cut_mode_lbl[cut_mode])) {
                for (int i = 0; i < (int)ARRAY_SIZE(cut_mode_lbl); ++i) {
                    if (ImGui::Selectable(cut_mode_lbl[i], cut_mode == i)) cut_mode = i;
                }
                ImGui::EndCombo();
            }
            ImGui::SetItemTooltip("How the horizontal cut (I vs q_par) is positioned.\nDragging the line on the map switches to a fixed q_z.");
            if (cut_mode == CutMode_AlphaF) {
                ImGui::InputFloat("alpha_f (°)", &cut_alpha_f_deg, 0, 0, "%.4f");
            }
            ImGui::BeginDisabled(cut_mode != CutMode_Qz);
            ImGui::InputDouble("q_z (nm⁻¹)", &cut_qz, 0, 0, "%.4f");
            ImGui::EndDisabled();
            ImGui::InputFloat("q_z band (nm⁻¹)", &cut_width_qz, 0, 0, "%.4f");
            ImGui::InputDouble("q_par (nm⁻¹)", &cut_qpar, 0, 0, "%.4f");
            ImGui::InputFloat("q_par band (nm⁻¹)", &cut_width_qpar, 0, 0, "%.4f");
            cut_width_qz = MAX(cut_width_qz, 0.0f);
            cut_width_qpar = MAX(cut_width_qpar, 0.0f);
            ImGui::SetItemTooltip("Width of the integration band. 0 uses the nearest row / ring.");
            ImGui::Checkbox("Log q_par axis", &cut_log_x);
            ImGui::Checkbox("Vertical cut vs alpha_f", &vcut_vs_alpha_f);
            ImGui::Checkbox("Bead self scattering", &show_bead_self);
            ImGui::SetItemTooltip("Show the (Born) self scattering of the discrete beads in the horizontal cut.\n"
                                  "Where the computed intensity follows this line, bead discreteness dominates.");
            if (res_rows) {
                ImGui::TextDisabled("Averaged over %zu rows / %zu rings", hcut_rows, vcut_cols);
            }
            ImGui::BeginDisabled(!res_rows);
            if (ImGui::Button("Export horizontal cut (CSV)")) export_cut(true);
            if (ImGui::Button("Export vertical cut (CSV)")) export_cut(false);
            ImGui::EndDisabled();
            ImGui::EndDisabled();
        }

        ImGui::PopItemWidth();
    }

    // ---------------------------------------------------------------------------------------------
    // Cuts
    // ---------------------------------------------------------------------------------------------

    // Horizontal cut: I(q_par) averaged over the q_z band centered at cut_qz
    // Vertical cut:   I(q_z)   averaged over the q_par band centered at cut_qpar
    void compute_cuts() {
        const size_t R = res_cols;
        const size_t Q = res_rows;
        md_array_resize(hcut, R, alloc);
        md_array_resize(vcut, Q, alloc);

        const double half_qz = MAX(0.5 * cut_width_qz, 0.5 * res_dqz * 10.0);
        const double half_qp = MAX(0.5 * cut_width_qpar, 0.5 * res_dq * 10.0);

        for (size_t c = 0; c < R; ++c) hcut[c] = 0.0;
        size_t nrows = 0;
        for (size_t r = 0; r < Q; ++r) {
            if (fabs(res_qz_nm[r] - cut_qz) > half_qz) continue;
            for (size_t c = 0; c < R; ++c) hcut[c] += res_smooth[r * R + c];
            nrows += 1;
        }
        if (nrows) for (size_t c = 0; c < R; ++c) hcut[c] /= (double)nrows;
        hcut_rows = nrows;

        size_t ncols = 0;
        size_t c_beg = R, c_end = 0;
        for (size_t c = 0; c < R; ++c) {
            if (fabs(res_ring_q_nm[c] - cut_qpar) > half_qp) continue;
            c_beg = MIN(c_beg, c);
            c_end = MAX(c_end, c + 1);
            ncols += 1;
        }
        for (size_t r = 0; r < Q; ++r) {
            double sum = 0.0;
            for (size_t c = c_beg; c < c_end; ++c) sum += res_smooth[r * R + c];
            vcut[r] = ncols ? sum / (double)ncols : 0.0;
        }
        vcut_cols = ncols;

        md_array_resize(vcut_af, Q, alloc);
        for (size_t r = 0; r < Q; ++r) vcut_af[r] = alpha_f_deg(res_qz_nm[r]);
    }

    // Mean SLD (relative to the ambient) of the film from the density profile, over the slices above half the maximum
    double film_sld_rel() const {
        const size_t n = md_array_size(prof_rho);
        if (!n) return 0.0;
        double vmax = 0.0;
        for (size_t i = 0; i < n; ++i) vmax = MAX(vmax, prof_rho[i]);
        double sum = 0.0; size_t cnt = 0;
        for (size_t i = 0; i < n; ++i) {
            if (prof_rho[i] >= 0.5 * vmax) { sum += prof_rho[i]; cnt += 1; }
        }
        return cnt ? eval_model.profile_sld_scale * sum / cnt : 0.0;
    }

    double film_yoneda_qz_nm() const {
        const double k0 = 2.0 * 3.14159265358979323846 / eval_model.wavelength;
        const double p = k0 * sin(eval_model.alpha_i);
        const double d = film_sld_rel();
        const double sin_ac = d > 0.0 ? sqrt(4.0 * 3.14159265358979323846 * d) / k0 : 0.0;
        return (p + k0 * sin_ac) * 10.0;
    }

    double yoneda_qz_nm() const {
        const double k0 = 2.0 * 3.14159265358979323846 / eval_model.wavelength;
        const double p = k0 * sin(eval_model.alpha_i);
        const double d_sld = eval_model.sld_substrate - eval_model.sld_ambient;
        const double sin_ac = d_sld > 0.0 ? sqrt(4.0 * 3.14159265358979323846 * d_sld) / k0 : 0.0;
        return (p + k0 * sin_ac) * 10.0;
    }

    double horizon_qz_nm() const {
        const double k0 = 2.0 * 3.14159265358979323846 / eval_model.wavelength;
        return k0 * sin(eval_model.alpha_i) * 10.0;
    }

    double qz_nm_from_alpha_f(double af_deg) const {
        const double k0 = 2.0 * 3.14159265358979323846 / eval_model.wavelength;
        return k0 * (sin(af_deg * DEG_TO_RAD) + sin(eval_model.alpha_i)) * 10.0;
    }

    double alpha_f_deg(double qz_nm) const {
        const double k0 = 2.0 * 3.14159265358979323846 / eval_model.wavelength;
        const double s = qz_nm * 0.1 / k0 - sin(eval_model.alpha_i);
        return fabs(s) <= 1.0 ? asin(s) * RAD_TO_DEG : 0.0;
    }

    // Index range [beg, end) of strictly positive values, for log scale plotting
    static void positive_range(const double* v, size_t n, size_t* beg, size_t* end) {
        size_t b = 0;
        while (b < n && !(v[b] > 0.0)) ++b;
        size_t e = n;
        while (e > b && !(v[e - 1] > 0.0)) --e;
        *beg = b;
        *end = e;
    }

    void export_cut(bool horizontal) {
        char path_buf[2048];
        if (!application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, STR_LIT("csv"))) return;
        str_t path = {path_buf, strnlen(path_buf, sizeof(path_buf))};

        const size_t n = horizontal ? res_cols : res_rows;
        md_array(float) x = nullptr;
        md_array(float) y = nullptr;
        md_array(float) a = nullptr;
        defer { md_array_free(x, alloc); md_array_free(y, alloc); md_array_free(a, alloc); };
        md_array_resize(x, n, alloc);
        md_array_resize(y, n, alloc);
        md_array_resize(a, n, alloc);
        for (size_t i = 0; i < n; ++i) {
            x[i] = (float)(horizontal ? res_ring_q_nm[i] : res_qz_nm[i]);
            y[i] = (float)(horizontal ? hcut[i] : vcut[i]);
            a[i] = horizontal ? (float)alpha_f_deg(cut_qz) : (float)alpha_f_deg(res_qz_nm[i]);
        }
        const float* cols[3] = {x, y, a};
        str_t names[3] = {
            horizontal ? str_from_cstr("q_par [nm^-1]") : str_from_cstr("q_z [nm^-1]"),
            STR_LIT("I [sr^-1]"),
            STR_LIT("alpha_f [deg]"),
        };
        if (!md_csv_write_to_file(cols, names, 3, n, path)) {
            MD_LOG_ERROR("GISAXS: failed to write '%.*s'", (int)path.len, path.ptr);
        }
    }

    void draw_detector(float w, float h, ImVec4 col_h, ImVec4 col_v) {
        if (!ImPlot::BeginPlot("##gisaxs_detector", ImVec2(w, h), ImPlotFlags_NoLegend | ImPlotFlags_Equal)) return;
        ImPlot::SetupAxis(ImAxis_X1, "q_y [nm⁻¹]");
        ImPlot::SetupAxis(ImAxis_Y1, "q_z [nm⁻¹]");
        ImPlot::SetupAxesLimits(det_qy_min, det_qy_max, det_qz_min, det_qz_max, ImPlotCond_Once);
        ImPlot::SetupFinish();

        ImPlot::PlotHeatmap("##det", det_map, (int)det_rows, (int)det_cols, res_min, res_max, nullptr,
                            ImPlotPoint(det_qy_min, det_qz_min), ImPlotPoint(det_qy_max, det_qz_max));

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

        if (ImPlot::IsPlotHovered()) {
            // Map the mouse back to a detector pixel through the same linear axes used for display
            const ImPlotPoint mp = ImPlot::GetPlotMousePos();
            const int b = MAX(det.binning, 1);
            const double fx = (mp.x - det_qy_min) / (det_qy_max - det_qy_min);
            const double fy = (mp.y - det_qz_min) / (det_qz_max - det_qz_min);
            if (fx >= 0 && fx < 1 && fy >= 0 && fy < 1) {
                const int bx = (int)(fx * det_cols);
                const int by = (int)(fy * det_rows);
                const double px = (bx + 0.5) * b, py = (by + 0.5) * b;
                double qpar, qy, qz, af;
                pixel_q((px - det.beam_x_px) * det.pixel_mm, (py - det.beam_y_px) * det.pixel_mm, &qpar, &qy, &qz, &af);
                const float I = det_raw[(size_t)by * det_cols + bx];
                if (I >= 0.0f) {
                    ImGui::SetTooltip("Pixel: %.0f, %.0f\nq_y: %.4f nm⁻¹\nq_z: %.4f nm⁻¹\nq_par: %.4f nm⁻¹\nalpha_f: %.3f°\nI: %.4e sr⁻¹",
                        px, py, qy, qz, qpar, af, I);
                } else {
                    ImGui::SetTooltip("Pixel: %.0f, %.0f\nMasked (gap, beamstop, shadow or outside the computed q-range)", px, py);
                }
            }
        }
        ImPlot::EndPlot();
    }

    void draw_plot() {
        if (!res_rows || !res_cols || !res_map) {
            ImGui::TextDisabled("No result yet. Press Compute.");
            return;
        }

        // Keep cut positions valid
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

        ImPlot::PushColormap(colormap);
        ImPlot::PushStyleColor(ImPlotCol_PlotBg, ImPlot::SampleColormap(0.0f, colormap));

        if (view_detector && det_rows && det_cols) {
            draw_detector(map_w, map_h, col_h, col_v);
        } else if (ImPlot::BeginPlot("##gisaxs_map", ImVec2(map_w, map_h), ImPlotFlags_NoLegend)) {
            ImPlot::SetupAxis(ImAxis_X1, "q_par [nm⁻¹]");
            ImPlot::SetupAxis(ImAxis_Y1, "q_z [nm⁻¹]");
            ImPlot::SetupAxesLimits(res_qpar_min, res_qpar_max, res_qz_min, res_qz_max, ImPlotCond_Once);
            ImPlot::SetupAxisLinks(ImAxis_X1, &link_qpar_min, &link_qpar_max);
            ImPlot::SetupAxisLinks(ImAxis_Y1, &link_qz_min, &link_qz_max);
            ImPlot::SetupFinish();

            ImPlot::PlotHeatmap("##I", res_map, (int)res_rows, (int)res_cols, res_min, res_max, nullptr,
                                ImPlotPoint(res_qpar_min, res_qz_min), ImPlotPoint(res_qpar_max, res_qz_max));

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
                const double qpar = mp.x * 0.1;   // Å^-1
                const double qz = mp.y * 0.1;
                const long col = (long)floor(qpar / res_dq + 0.5) - 1;
                const long row = res_dqz > 0 ? (long)floor(qz / res_dqz + 0.5) : -1;
                if (col >= 0 && col < (long)res_cols && row >= 0 && row < (long)res_rows) {
                    const float I = res_smooth[row * res_cols + col];
                    ImGui::SetTooltip("q_par: %.4f nm⁻¹\nq_z: %.4f nm⁻¹\nalpha_f: %.3f°\nI: %.4e sr⁻¹\nRing points: %u",
                        mp.x, mp.y, alpha_f_deg(mp.y), I, (col < (long)md_array_size(res_ring_count)) ? res_ring_count[col] : 0u);
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
        bool first = true;

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
                size_t b, e;
                positive_range(hcut, res_cols, &b, &e);
                if (!log_scale) { b = 0; e = res_cols; }
                if (e > b) {
                    ImPlot::SetNextLineStyle(col_h, 1.5f);
                    ImPlot::PlotLine("I", res_ring_q_nm + b, hcut + b, (int)(e - b));
                }
                if (show_bead_self && !bead_self.empty()) {
                    // Born self scattering of the discrete beads: sum_j w_j^2 exp(-(q_par^2 + q_z^2) sigma_j^2) / A.
                    // Where the computed curve approaches it, the signal is dominated by bead discreteness.
                    const double qz_a = cut_qz * 0.1;
                    const double scl = MD_GISAXS_R_E * MD_GISAXS_R_E * contrast_factor() / bead_self_area;
                    md_array(double) selfI = nullptr;
                    md_array_resize(selfI, res_cols, alloc);
                    for (size_t c = 0; c < res_cols; ++c) {
                        const double qp = res_ring_q_nm[c] * 0.1;
                        double sum = 0.0;
                        for (const auto& g : bead_self) sum += g.second * exp(-(qp * qp + qz_a * qz_a) * g.first * g.first);
                        selfI[c] = sum * scl;
                    }
                    ImPlot::SetNextLineStyle(ImVec4(0.8f, 0.8f, 0.8f, 0.7f), 1.0f);
                    ImPlot::PlotLine("Bead self scattering", res_ring_q_nm, selfI, (int)res_cols);
                    md_array_free(selfI, alloc);
                }
                ImPlot::SetNextLineStyle(ImVec4(col_v.x, col_v.y, col_v.z, 0.6f), 1.0f);
                ImPlot::PlotInfLines("##qpar", &cut_qpar, 1);
                ImPlot::EndPlot();
            }
            first = false;

            // Vertical cut: I(q_z)
            ImGui::SameLine();
            snprintf(title, sizeof(title), "q_par = %.3f nm⁻¹###vcut", cut_qpar);
            if (ImPlot::BeginPlot(title, ImVec2(plot_w, lower_h), ImPlotFlags_NoLegend)) {
                // Either q_z (linked with the map) or the exit angle alpha_f (as in most experimental papers)
                const bool   use_af = vcut_vs_alpha_f;
                const double* xs = use_af ? vcut_af : res_qz_nm;
                auto to_x = [&](double qz_nm) { return use_af ? alpha_f_deg(qz_nm) : qz_nm; };
                ImPlot::SetupAxis(ImAxis_X1, use_af ? "alpha_f [°]" : "q_z [nm⁻¹]", use_af ? ImPlotAxisFlags_AutoFit : 0);
                ImPlot::SetupAxis(ImAxis_Y1, "I [sr⁻¹]", ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisScale(ImAxis_Y1, y_scale);
                if (!use_af) ImPlot::SetupAxisLinks(ImAxis_X1, &link_qz_min, &link_qz_max);
                ImPlot::SetupFinish();
                size_t b, e;
                positive_range(vcut, res_rows, &b, &e);
                if (!log_scale) { b = 0; e = res_rows; }
                if (use_af && eval_model.dwba) {
                    // Only above the sample horizon
                    const double horizon = horizon_qz_nm();
                    while (b < e && res_qz_nm[b] < horizon) ++b;
                }
                if (e > b) {
                    ImPlot::SetNextLineStyle(col_v, 1.5f);
                    ImPlot::PlotLine("I", xs + b, vcut + b, (int)(e - b));
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

        if (show_profile && prof_z && md_array_size(prof_z)) {
            if (!first) ImGui::SameLine();
            if (ImPlot::BeginPlot("Density profile##gisaxs_profile", ImVec2(plot_w, lower_h), ImPlotFlags_NoLegend)) {
                ImPlot::SetupAxis(ImAxis_X1, "z [Å]");
                ImPlot::SetupAxis(ImAxis_Y1, "rho_e [e/Å³]", ImPlotAxisFlags_AutoFit);
                ImPlot::SetupFinish();
                ImPlot::PlotLine("Profile", prof_z, prof_rho, (int)md_array_size(prof_z));
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
        viamd::write_section_header(state, STR_LIT("GISAXS"));
        viamd::write_str(state, STR_LIT("Filter"), str_from_cstr(filter));
        viamd::write_int(state, STR_LIT("ElectronMode"), electron_mode);
        viamd::write_flt(state, STR_LIT("Electrons"), electrons_per_particle);
        viamd::write_int(state, STR_LIT("SigmaMode"), sigma_mode);
        viamd::write_flt(state, STR_LIT("Sigma"), sigma_constant);
        viamd::write_flt(state, STR_LIT("MaterialDensity"), material_density);
        viamd::write_int(state, STR_LIT("BackgroundPreset"), bg_preset);
        viamd::write_dbl(state, STR_LIT("BackgroundDensity"), bg_density);
        viamd::write_bool(state, STR_LIT("DWBA"), dwba);
        viamd::write_int(state, STR_LIT("SubstratePreset"), sub_preset);
        viamd::write_dbl(state, STR_LIT("SubstrateDensity"), sub_density);
        viamd::write_dbl(state, STR_LIT("SubstrateBeta"), sub_beta);
        viamd::write_bool(state, STR_LIT("SubstrateAutoZ"), sub_auto_z);
        viamd::write_flt(state, STR_LIT("SubstrateZ"), sub_z);
        viamd::write_flt(state, STR_LIT("SubstrateZOffset"), sub_z_offset);
        viamd::write_flt(state, STR_LIT("EnergyKeV"), energy_kev);
        viamd::write_flt(state, STR_LIT("AlphaIDeg"), alpha_i_deg);
        viamd::write_flt(state, STR_LIT("QParMax"), q_par_max_nm);
        viamd::write_flt(state, STR_LIT("QZMax"), q_z_max_nm);
        viamd::write_flt(state, STR_LIT("Oversampling"), oversampling);
        viamd::write_int(state, STR_LIT("NumQz"), num_qz);
        viamd::write_int(state, STR_LIT("MaxSlices"), max_slices);
        viamd::write_bool(state, STR_LIT("LogScale"), log_scale);
        viamd::write_flt(state, STR_LIT("Decades"), decades);
        viamd::write_flt(state, STR_LIT("SubstrateRoughness"), sub_roughness);
        viamd::write_flt(state, STR_LIT("ElectronsPerDalton"), electrons_per_dalton);
        viamd::write_int(state, STR_LIT("DensityModel"), density_model);
        viamd::write_str(state, STR_LIT("FibrilCenter"), str_from_cstr(fib_center_name));
        viamd::write_int(state, STR_LIT("FibrilStride"), fib_stride);
        viamd::write_int(state, STR_LIT("FibrilTemplate"), fib_template);
        viamd::write_flt(state, STR_LIT("FibrilBeadSigma"), fib_bead_sigma);
        viamd::write_flt(state, STR_LIT("FibrilDiskRadius"), fib_disk_radius);
        viamd::write_flt(state, STR_LIT("FibrilDiskSpacing"), fib_disk_spacing);
        viamd::write_flt(state, STR_LIT("FibrilSampleSpacing"), fib_ds);
        viamd::write_str(state, STR_LIT("FibrilTemplatePath"), str_from_cstr(fib_template_path));
        viamd::write_bool(state, STR_LIT("ShowBeadSelf"), show_bead_self);
        viamd::write_bool(state, STR_LIT("FilmGraded"), film_graded);
        viamd::write_dbl(state, STR_LIT("MaterialBeta"), material_beta);
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
        viamd::write_bool(state, STR_LIT("ShowCuts"), show_cuts);
        viamd::write_int(state, STR_LIT("CutMode"), cut_mode);
        viamd::write_bool(state, STR_LIT("VCutVsAlphaF"), vcut_vs_alpha_f);
        viamd::write_flt(state, STR_LIT("CutAlphaF"), cut_alpha_f_deg);
        viamd::write_int(state, STR_LIT("BeamPreset"), beam_preset);
        viamd::write_dbl(state, STR_LIT("CutQz"), cut_qz);
        viamd::write_dbl(state, STR_LIT("CutQpar"), cut_qpar);
        viamd::write_flt(state, STR_LIT("CutWidthQz"), cut_width_qz);
        viamd::write_flt(state, STR_LIT("CutWidthQpar"), cut_width_qpar);
    }

    void deserialize(viamd::deserialization_state_t& state) {
        str_t ident, arg;
        while (viamd::next_entry(ident, arg, state)) {
            if      (str_eq(ident, STR_LIT("Filter")))            viamd::extract_to_char_buf(filter, sizeof(filter), arg);
            else if (str_eq(ident, STR_LIT("ElectronMode")))      viamd::extract_int(electron_mode, arg);
            else if (str_eq(ident, STR_LIT("Electrons")))         viamd::extract_flt(electrons_per_particle, arg);
            else if (str_eq(ident, STR_LIT("SigmaMode")))         viamd::extract_int(sigma_mode, arg);
            else if (str_eq(ident, STR_LIT("Sigma")))             viamd::extract_flt(sigma_constant, arg);
            else if (str_eq(ident, STR_LIT("MaterialDensity")))   viamd::extract_flt(material_density, arg);
            else if (str_eq(ident, STR_LIT("BackgroundPreset")))  viamd::extract_int(bg_preset, arg);
            else if (str_eq(ident, STR_LIT("BackgroundDensity"))) viamd::extract_dbl(bg_density, arg);
            else if (str_eq(ident, STR_LIT("DWBA")))              viamd::extract_bool(dwba, arg);
            else if (str_eq(ident, STR_LIT("SubstratePreset")))   viamd::extract_int(sub_preset, arg);
            else if (str_eq(ident, STR_LIT("SubstrateDensity")))  viamd::extract_dbl(sub_density, arg);
            else if (str_eq(ident, STR_LIT("SubstrateBeta")))     viamd::extract_dbl(sub_beta, arg);
            else if (str_eq(ident, STR_LIT("SubstrateAutoZ")))    viamd::extract_bool(sub_auto_z, arg);
            else if (str_eq(ident, STR_LIT("SubstrateZ")))        viamd::extract_flt(sub_z, arg);
            else if (str_eq(ident, STR_LIT("SubstrateZOffset")))  viamd::extract_flt(sub_z_offset, arg);
            else if (str_eq(ident, STR_LIT("EnergyKeV")))         viamd::extract_flt(energy_kev, arg);
            else if (str_eq(ident, STR_LIT("AlphaIDeg")))         viamd::extract_flt(alpha_i_deg, arg);
            else if (str_eq(ident, STR_LIT("QParMax")))           viamd::extract_flt(q_par_max_nm, arg);
            else if (str_eq(ident, STR_LIT("QZMax")))             viamd::extract_flt(q_z_max_nm, arg);
            else if (str_eq(ident, STR_LIT("Oversampling")))      viamd::extract_flt(oversampling, arg);
            else if (str_eq(ident, STR_LIT("NumQz")))             viamd::extract_int(num_qz, arg);
            else if (str_eq(ident, STR_LIT("MaxSlices")))         viamd::extract_int(max_slices, arg);
            else if (str_eq(ident, STR_LIT("LogScale")))          viamd::extract_bool(log_scale, arg);
            else if (str_eq(ident, STR_LIT("Decades")))           viamd::extract_flt(decades, arg);
            else if (str_eq(ident, STR_LIT("SubstrateRoughness"))) viamd::extract_flt(sub_roughness, arg);
            else if (str_eq(ident, STR_LIT("ElectronsPerDalton"))) viamd::extract_flt(electrons_per_dalton, arg);
            else if (str_eq(ident, STR_LIT("DensityModel")))      viamd::extract_int(density_model, arg);
            else if (str_eq(ident, STR_LIT("FibrilCenter")))      viamd::extract_to_char_buf(fib_center_name, sizeof(fib_center_name), arg);
            else if (str_eq(ident, STR_LIT("FibrilStride")))      viamd::extract_int(fib_stride, arg);
            else if (str_eq(ident, STR_LIT("FibrilTemplate")))    viamd::extract_int(fib_template, arg);
            else if (str_eq(ident, STR_LIT("FibrilBeadSigma")))   viamd::extract_flt(fib_bead_sigma, arg);
            else if (str_eq(ident, STR_LIT("FibrilDiskRadius")))  viamd::extract_flt(fib_disk_radius, arg);
            else if (str_eq(ident, STR_LIT("FibrilDiskSpacing"))) viamd::extract_flt(fib_disk_spacing, arg);
            else if (str_eq(ident, STR_LIT("FibrilSampleSpacing"))) viamd::extract_flt(fib_ds, arg);
            else if (str_eq(ident, STR_LIT("FibrilTemplatePath"))) viamd::extract_to_char_buf(fib_template_path, sizeof(fib_template_path), arg);
            else if (str_eq(ident, STR_LIT("ShowBeadSelf")))      viamd::extract_bool(show_bead_self, arg);
            else if (str_eq(ident, STR_LIT("FilmGraded")))        viamd::extract_bool(film_graded, arg);
            else if (str_eq(ident, STR_LIT("MaterialBeta")))      viamd::extract_dbl(material_beta, arg);
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
            else if (str_eq(ident, STR_LIT("ShowCuts")))          viamd::extract_bool(show_cuts, arg);
            else if (str_eq(ident, STR_LIT("CutMode")))           viamd::extract_int(cut_mode, arg);
            else if (str_eq(ident, STR_LIT("VCutVsAlphaF")))      viamd::extract_bool(vcut_vs_alpha_f, arg);
            else if (str_eq(ident, STR_LIT("CutAlphaF")))         viamd::extract_flt(cut_alpha_f_deg, arg);
            else if (str_eq(ident, STR_LIT("BeamPreset")))        viamd::extract_int(beam_preset, arg);
            else if (str_eq(ident, STR_LIT("CutQz")))             { viamd::extract_dbl(cut_qz, arg); cut_init = true; }
            else if (str_eq(ident, STR_LIT("CutQpar")))           { viamd::extract_dbl(cut_qpar, arg); cut_init = true; }
            else if (str_eq(ident, STR_LIT("CutWidthQz")))        viamd::extract_flt(cut_width_qz, arg);
            else if (str_eq(ident, STR_LIT("CutWidthQpar")))      viamd::extract_flt(cut_width_qpar, arg);
        }
        electron_mode = CLAMP(electron_mode, 0, (int)ARRAY_SIZE(electron_mode_lbl) - 1);
        density_model = CLAMP(density_model, 0, (int)ARRAY_SIZE(density_model_lbl) - 1);
        fib_template = CLAMP(fib_template, 0, (int)ARRAY_SIZE(fibril_template_lbl) - 1);
        sigma_mode = CLAMP(sigma_mode, 0, 1);
        bg_preset = CLAMP(bg_preset, 0, (int)ARRAY_SIZE(background_presets) - 1);
        sub_preset = CLAMP(sub_preset, 0, (int)ARRAY_SIZE(substrate_presets) - 1);
        cut_mode = CLAMP(cut_mode, 0, (int)ARRAY_SIZE(cut_mode_lbl) - 1);
        beam_preset = CLAMP(beam_preset, -1, (int)ARRAY_SIZE(beam_presets) - 1);
    }
};

static Gisaxs instance;
