#define IMGUI_DEFINE_MATH_OPERATORS

#include <event.h>
#include <viamd_event.h>
#include <viamd.h>
#include <task_system.h>
#include <color_utils.h>

#include <gfx/volumerender_utils.h>
#include <gfx/gl_utils.h>
#include <gfx/immediate_draw_utils.h>

#include <md_gto.h>
#include <md_vlx.h>
#include <md_util.h>
#include <md_topo.h>
#include <md_vector_graphics.h>
#include <md_flow.h>
#include <md_csv.h>
#include <md_xvg.h>
#include <core/md_vec_math.h>
#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>

#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#endif

#include <serialization_utils.h>

#include <imgui_widgets.h>
#include <imgui_internal.h>

#include "flow_draw.h"

#include <implot_widgets.h>
#include <implot_internal.h>

#include <algorithm>
#include <app/IconsFontAwesome6.h>

#define BLK_DIM 8
#define ROTATORY_STRENGTH_IN_CGS 471.443648175
#define EXTINCTION_COEFFICIENT_FROM_BETA 19.603697575813566
#define ROTATORY_STRENGTH_TO_DELTA_EPSILON 0.01386075702557538652

#define HARTREE_TO_KJ_PER_MOL 2625.4996394799
#define HARTREE_TO_EV 27.2114079527
#define EV_TO_HARTREE (1.0 / HARTREE_TO_EV)

#define AU_TO_GM 40479.02797814119
#define OSCILLATOR_STRENGTH_TO_EPSILON 1054.9516366171872

#define DEFAULT_SAMPLES_PER_ANGSTROM 8
#define DEFAULT_GTO_CUTOFF_VALUE 1.0e-6
#define MAX_NTO_GROUPS 16
#define MAX_NTO_LAMBDAS 3
#define NTO_LAMBDA_CUTOFF_VALUE 0.1

#define FORCE_CPU_PATH 0

#define USE_AABB_GRID 0

// Resolution for broadened plots
#define NUM_SAMPLES 1024

#define IM_GREEN ImVec4(0, 1, 0, 1)
#define IM_RED ImVec4(1, 0, 0, 1)
#define IM_YELLOW ImVec4(1, 1, 0.5, 0.3)
#define IM_BLUE ImVec4(0.5, 0.5, 1, 0.3)

#define COLOR_PEAK_SELECTED ImVec4(0.5, 0.5, 1.0, 1.0)
#define COLOR_PEAK_HOVER    ImVec4(1.0, 1.0, 0.5, 1.0)
// COLOR_PEAK_HOVER as a table row background. Kept translucent: the row text and the selectable's
// own fill both draw on top of it, and an opaque tint makes either unreadable.
#define COLOR_ROW_HIGHLIGHT ImVec4(1.0, 1.0, 0.5, 0.25)

#define COLOR_TABLE_SELECTED IM_BLUE
#define COLOR_TABLE_HOVERED  IM_YELLOW

#define U32_MAGENTA IM_COL32(255, 0, 255, 255)

#define U32_VELOXCHEM_GREEN IM_COL32(0, 162, 135, 191)
#define VEC4_VELOXCHEM_GREEN {0, 162.0f/255.0f, 135.0f/255.0f, 0.75f}

// Complement to veloxchem green (magentaish)
#define VEC4_VELOXCHEM_MAGENTA {162.0f/255.0f, 35.0f/255.0f, 135.0f/255.0f, 0.75f}

constexpr uint64_t interaction_surface_nto = HASH_STR_LIT64("interaction surface nto");
constexpr uint64_t interaction_surface_orb = HASH_STR_LIT64("interaction surface orb");

constexpr PickingDomainID PickingDomain_CriticalPoints = HASH_STR_LIT64("Picking Domain Critical Points");

// This is the internal storage order of volumes and textures related to NTOs
enum NTO {
    NTO_Attachment,
    NTO_Detachment,
};

// The reference the SCF was solved against. It belongs with the method and the basis set because
// it is part of how the calculation is named - RKS/B3LYP/def2-SVP - and because it is what decides
// whether the beta orbitals are their own set or a second name for alpha's.
static const char* vlx_scf_type_str(md_vlx_scf_type_t type) {
    switch (type) {
    case MD_VLX_SCF_RESTRICTED:           return "Restricted";
    case MD_VLX_SCF_RESTRICTED_OPENSHELL: return "Restricted Open-Shell";
    case MD_VLX_SCF_UNRESTRICTED:         return "Unrestricted";
    case MD_VLX_SCF_UNKNOWN:              // fallthrough
    default:                              return "Unknown";
    }
}

static const char* electronic_structure_value_mode_str(ElectronicStructureSource source, ElectronicStructureSpin spin, bool use_magnitude) {
    switch (source) {
    case ElectronicStructureSource::MolecularOrbital:
    case ElectronicStructureSource::NaturalTransitionOrbital:
        return use_magnitude ? "Magnitude" : "Signed";
    case ElectronicStructureSource::TransitionDensity:
        return "Density";
    case ElectronicStructureSource::ElectronDensity:
        return spin == ElectronicStructureSpin::Difference ? (use_magnitude ? "Magnitude" : "Signed") : "Density";
    default:
        return "";
    }
}

// What the Charge Transfer diagram routes its flow through.
//
// Groups is the M1 view and reproduces the old two-column diagram exactly; it is kept as the
// comparison the NTO view can be checked against, not as a feature.
enum flow_mode_t {
    FLOW_MODE_GROUPS = 0,   // group -> group, via compute_transition_matrix
    FLOW_MODE_NTO,          // atom/group -> NTO pair -> atom/group
    FLOW_MODE_ORBITALS,     // occupied MO -> NTO pair -> virtual MO
    FLOW_MODE_COUNT,
};

static const char* flow_mode_str[] = { "Groups only", "Atoms via NTOs", "Orbitals" };

enum x_unit_t {
    X_UNIT_EV,
    X_UNIT_NM,
    X_UNIT_CM_INVERSE,
    X_UNIT_HARTREE,
    X_UNIT_COUNT,
};

static const char* x_unit_full_str[] = {"Energy (eV)", "Wavelength (nm)", (const char*)u8"Wavenumber (cm⁻¹)", "Energy (au)"};
static const char* x_unit_label_str[] = {"Energy", "Wavelength", "Wavenumber", "Energy"};
static const char* x_unit_short_str[] = {"eV", "nm", (const char*)u8"cm⁻¹", "au"};

enum broadening_mode_t {
    BROADENING_MODE_GAUSSIAN,
    BROADENING_MODE_LORENTZIAN,
    BROADENING_MODE_COUNT,
};

static const char* broadening_mode_str[] = { "Gaussian", "Lorentzian" };

// Attribute charge contributions of an AO-basis density matrix D to atoms or groups
// using a Mulliken-style partitioning: q_g = sum_{mu in g} sum_nu D[mu,nu] * S[mu,nu].
// Both D and S are symmetric (num_ao x num_ao) row-major matrices.
// The total tr(D*S) is exactly preserved: sum_g out_charge[g] == tr(D*S).
static void attribute_charge_density(double out_charge[], const int* ao_to_idx, const double* D, const double* S, size_t num_ao) {
	ASSERT(D);
	ASSERT(S);
	ASSERT(ao_to_idx);

	for (size_t mu = 0; mu < num_ao; ++mu) {
		const int idx = ao_to_idx[mu];
		const double* D_row = &D[mu * num_ao];
		const double* S_row = &S[mu * num_ao];
		double s = 0.0;
		for (size_t nu = 0; nu < num_ao; ++nu) {
			s += D_row[nu] * S_row[nu];
		}
		out_charge[idx] += s;
	}
}

// Mulliken partitioning of a SINGLE orbital vector, which is the same quantity as above for the
// rank-one density D = C C^T. Written out rather than routed through attribute_charge_density so
// it does not have to materialize an N x N matrix per orbital - with 32 NTOs and a few hundred
// AOs, forming those densities is most of the work and none of the answer.
//
//   q_g = sum_{mu in g} C[mu] * (S C)[mu],  and  sum_g q_g = C^T S C = 1 for a normalized orbital.
static void attribute_orbital_density(double out_charge[], const int* ao_to_idx, const double* C, const double* S, size_t num_ao) {
	ASSERT(C);
	ASSERT(S);
	ASSERT(ao_to_idx);

	for (size_t mu = 0; mu < num_ao; ++mu) {
		const double* S_row = &S[mu * num_ao];
		double sc = 0.0;
		for (size_t nu = 0; nu < num_ao; ++nu) {
			sc += S_row[nu] * C[nu];
		}
		out_charge[ao_to_idx[mu]] += C[mu] * sc;
	}
}

struct VeloxChem : viamd::EventHandler {
    VeloxChem() { viamd::event_system_register_handler(*this); }
    md_vlx_t* vlx = nullptr;

    bool use_gpu_path = false;
#if MD_ENABLE_GPU
    // No GPU handles are kept here any more. The device, its stream, its pools and the evaluation
    // scratch belong to the application; the uploaded basis and the atom buffer belong to the
    // dataset. Everything this component does with them goes through the application's evaluators,
    // which read them from the ApplicationState it is already handed - so a stale copy taken at
    // load time is not something that can exist.
    //
    // The one exception is the critical point extraction below, which records its own kernel over
    // the scratch volume. It reads state.gpu_* at the point of use, for the same reason.
#endif

    // GL representations
    md_gl_rep_t gl_rep = {};

    int homo_idx[2] = {};
    int lumo_idx[2] = {};

    AABB aabb = {};
    OABB oabb = {};

    dvec3_t center_of_charge = {};

    // Up to date packed atom xyz in BOHR ready to be supplied for evaluation
    // The w component is padding for now. 16-bytes are used for better alignment (gpu buffer 1:1)
    vec4_t* atom_xyzw = nullptr;

    // This is the default view which is used as a reset view target
    ViewTransform default_view = {};

    // If this is not NULL, then it means the qm data represented in the vlx object is a subset of the actual system
    // This array then maps the local qm indices into system-wide atom indices.
    const int32_t* qm_to_atom_idx = nullptr;
    md_bitfield_t* qm_mask = nullptr;

    struct Summary {
        bool show_window = false;
    } summary;

    struct CriticalPoints {
        bool enabled = false;

        uint64_t raw_hash = 0;
        uint64_t simp_hash = 0;

#if MD_ENABLE_GPU
        md_topo_gpu_context_t* topo_ctx      = nullptr;
        uint32_t               topo_dims[3]  = {};
#endif
        md_topo_extremum_graph_t raw_graph  = { 0 };
        md_topo_extremum_graph_t simp_graph = { 0 };

        PickingRange picking_range = {};

        md_bitfield_t selection_mask = {};
        md_bitfield_t highlight_mask = {};
    } critical_points;

    struct Opt {
        int hovered = -1;
        int selected = -1;

        bool coord_modified = false;

        uint64_t energy_hash = 0;
    } opt;

    struct Orb {
        bool show_window = false;
        Volume   vol[16] = {};
        int      vol_mo_idx[16] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
        md_vlx_spin_t vol_mo_type[16] = {};
        uint32_t iso_tex[16] = {};
        task_system::ID vol_task[16] = {};
        int num_x  = 3;
        int num_y  = 3;
        int mo_idx = -1;
        int scroll_to_idx = -1;

        IsoDesc iso = {
            .count = 2,
            .values = {0.05f, -0.05},
            .colors = {{0.f/255.f,75.f/255.f,135.f/255.f,0.75f}, {255.f/255.f,205.f/255.f,0.f/255.f,0.75f}},
        };

        GBuffer gbuf = {};
        PickingSurface picking_surface = {};
        Camera camera = {};
        ViewTransform target;

        float distance_scale = 2.0f;
        bool show_coordinate_system_widget = true;
    } orb;

    struct Nto {
        bool show_window = false;
        Volume   vol[2] = {};
        uint32_t iso_tex[2] = {};

        int sel_nto_idx = -1;

        md_grid_t grid = {0};

        md_array(uint32_t) atom_group_idx = nullptr;

        md_gl_rep_t gl_rep = {0};

        size_t transition_matrix_dim = 0;
        float* transition_matrix = nullptr;
        float* transition_density_hole = nullptr;
        float* transition_density_part = nullptr;

        struct {
            size_t count = 0;

            char   label[MAX_NTO_GROUPS][64] = {};
            vec4_t color[MAX_NTO_GROUPS] = {};
            int8_t hovered_index = -1;
        } group;

        double iso_val = 0.05;
        vec4_t col_pos = { 0.0f, 0.294f, 0.529f, 0.75f };
        vec4_t col_neg = { 1.0f, 0.804f, 0.0f,   0.75f };
        vec4_t col_den = { 1.0f, 1.0f,   1.0f,   0.75f };
        vec4_t col_att = VEC4_VELOXCHEM_GREEN;
        vec4_t col_det = VEC4_VELOXCHEM_MAGENTA;

        struct {
            bool enabled = false;
            vec4_t colors[3] = {
                {0.0f, 1.0f, 1.0f, 1.0f},
                {1.0f, 0.0f, 1.0f, 1.0f}
            };
            float vector_scale = 1.0f;
            bool show_angle = false;
        } dipole;

        GBuffer gbuf = {};
        PickingSurface picking_surface = {};
        Camera camera = {};
		ViewTransform target = {};

        float distance_scale = 2.0f;
        bool link_attachment_detachment_density = false;
        bool show_coordinate_system_widget = true;
    } nto;

    // Charge Transfer Analysis: the layered flow diagram (docs/transition_flow_design.md).
    //
    // Separate window, separate state, but the SAME group definitions as the NTO window - asking
    // someone to define their groups twice in two windows would be worse than any coupling this
    // introduces.
    struct Flow {
        bool show_window = false;

        flow_mode_t mode = FLOW_MODE_NTO;

        md_flow_graph_t  graph  = {};
        md_flow_cut_t    cut    = {};
        md_flow_layout_t layout = {};
        md_flow_layout_params_t params = {};

        // Node labels are borrowed by the graph rather than owned, so they have to outlive it.
        // Group labels already live in nto.group.label; these two are ours.
        char nto_label[MD_VLX_NTO_MAX_LAMBDAS][40] = {};

        // Second line: the dominant MO character of each NTO, e.g. "HOMO-1 -> LUMO". Empty when
        // no single pair dominates enough to name honestly.
        char nto_character[MD_VLX_NTO_MAX_LAMBDAS][40] = {};

        // Fixed-stride blocks, allocated once per file. The graph BORROWS these strings, so they
        // have to outlive it; rebuilding them per graph rebuild would leak into the arena on every
        // state change. Atom labels never change at all; MO labels are rewritten in place.
        char*  atom_labels = nullptr;
        size_t num_atom_labels = 0;
        char*  mo_labels = nullptr;
        size_t num_mo_labels = 0;

        // Most atoms contribute almost nothing to a given transition, so the diagram opens with a
        // threshold rather than with a hundred hairline ribbons. The slider is in the toolbar.
        float initial_threshold = 0.01f;

        // How hard to flatten the ORBITAL columns' heights. 0 = strictly proportional, 1 = nearly
        // uniform. Applied only to orbital columns: an atom column's proportions ARE the charge
        // transfer story, while an NTO column is a channel picker, and a channel drawn two pixels
        // tall cannot be read, labelled or clicked.
        float size_normalize = 0.75f;

        // NTO pairs below this share of the transition are dropped and the rest renormalized.
        // Distinct from the cut threshold, which hides nodes that are still in the graph.
        float lambda_cutoff = 1.0e-3f;

        FlowDrawStyle style = {};

        // View transform, in graph units. Kept here rather than in the layout so that panning and
        // zooming never trigger a relayout.
        vec2_t pan  = {0, 0};
        float  zoom = 1.0f;

        // Last frame's hit, so the dimming does not flicker while the cursor crosses a seam.
        FlowHit hover = {};

        // Deferred because fitting needs a laid-out graph and a canvas size, neither of which the
        // rebuild has. Consumed on the next draw.
        bool fit_requested = true;

        // Rubber-band region select. The button is captured at press time rather than read live,
        // so a drag started with the left button stays an ADD even if the right one gets pressed
        // part way through.
        bool   region_active = false;
        bool   region_add    = true;
        ImVec2 region_start  = {0, 0};

        // ATOM selection is NOT held here. It lives in state.selection.selection_mask, and this
        // window reads it every frame - so a selection made in the 3D view, by a script, or in any
        // other window shows up here, and one made here shows up everywhere else. A private copy
        // would be a second source of truth that silently drifts from the first.
        //
        // ORBITAL selection - NTO pairs and canonical MOs - is a different quantity: it picks a
        // transition channel rather than atoms, so it has nowhere else to live and is held here,
        // by key, like expansion.
        md_array(uint64_t) selected_orbital_keys = nullptr;

        uint64_t data_hash = 0;
        bool built = false;
    } flow;

    // Workspace group assignments, waiting for a molecule.
    //
    // load_workspace parses the ENTIRE workspace text before it loads any data, so when the
    // [VeloxChem] section is read there is no vlx object and no atom count yet. Groups are
    // per-atom, so they have to be buffered and stamped on at the end of init_from_file.
    //
    // The buffer lives on the heap, NOT on 'arena': reset_data() resets that arena, and a
    // workspace load goes through reset_data before it gets anywhere near init_from_file.
    struct PendingGroups {
        bool   pending = false;
        size_t count = 0;
        char   label[MAX_NTO_GROUPS][64] = {};
        vec4_t color[MAX_NTO_GROUPS] = {};
        md_bitfield_t atoms[MAX_NTO_GROUPS] = {};
        bool   atoms_initialized = false;
    } pending_groups;

    struct Rsp {
        bool show_window = false;
        bool show_export_window = false;
        int hovered = -1;
        int selected = -1;

        // Broadening parameter
        double broadening_fwhm = 0.123; // FWHM in eV, default value is 0.123 eV
        broadening_mode_t broadening_mode = BROADENING_MODE_LORENTZIAN;
        x_unit_t x_unit = X_UNIT_EV;

        uint64_t hash = 0; // Hash for tracking the unit of the x-axis, if it changes we need to recalculate the spectra (and shown peaks)
    } rsp;

    // Settings + cached results for the RIXS map.
    //
    // The struct owns its own arena ('alloc'), created at ViamdInitialize and released at
    // ViamdShutdown. The cached arrays below are the *only* thing allocated from it, and
    // compute_rixs_map() rebuilds all of them from scratch, so the arena is simply reset at the top
    // of that function: memory never exceeds the footprint of a single map, and no array can be
    // leaked by forgetting to free it individually.
    //
    // A vm arena is used rather than md_arena_allocator because reset() is a pointer rewind that
    // keeps the pages committed. The map is rebuilt on every settings change (including while a
    // slider is being dragged), and md_arena_allocator_reset() would release every page back to the
    // backing allocator only to re-acquire it on the next frame.
    //
    // NOTE: 'alloc' must survive the value-initialization of this struct in reset_data().
    struct Rixs {
        // --- Settings (mirrors the keyword arguments of spectrumplot.plot_rixs_map) ---

        // Include the elastic (Rayleigh) line as an additional stick per photon energy.
        bool plot_elastic_line = false;

        // Normalize the map by its global maximum ('global_max' in the python version).
        bool normalize = true;

        broadening_mode_t broadening_mode = BROADENING_MODE_LORENTZIAN;
        double broadening_fwhm_ev = 0.24;   // FWHM in eV applied along the energy loss axis
        double x_step_ev = 0.01;            // Grid spacing along the energy loss axis, in eV
        ImPlotColormap colormap = ImPlotColormap_Viridis;

        // Arena backing every pointer below. Owned by the module, see the comment above the struct.
        md_allocator_i* alloc = nullptr;

        // --- Cached results, rebuilt whenever 'hash' changes ---
        uint64_t hash  = 0;
        bool     valid = false;

        // The map, row-major [num_rows][num_cols].
        // NOTE: row 0 holds the HIGHEST photon energy. ImPlot::PlotHeatmap renders the first row at
        // bounds_max.y (top), so reversing the row order here reproduces matplotlib's origin='lower'.
        double* map = nullptr;

        double* grid      = nullptr;  // [num_cols] energy loss in eV
        double* photon_ev = nullptr;  // [num_rows] incoming photon energies in eV, ascending
        double* xas_x     = nullptr;  // [num_xas]  broadened XAS energies in eV
        double* xas_y     = nullptr;  // [num_xas]  broadened XAS intensities
        double* core_ev   = nullptr;  // [num_core] core excitation energies in eV, ascending
        double* core_f    = nullptr;  // [num_core] core oscillator strengths, sorted with core_ev

        size_t num_rows = 0;  // == number of incoming photon energies (P)
        size_t num_cols = 0;  // == number of grid samples along the energy loss axis (N)
        size_t num_xas  = 0;
        size_t num_core = 0;

        double map_min = 0.0, map_max = 0.0;
        double x_min_ev = 0.0, x_max_ev = 0.0;  // extent of the energy loss axis
        double y_min_ev = 0.0, y_max_ev = 0.0;  // extent of the photon energy axis (pixel edges)
        double xas_max  = 0.0;

        // Linked photon-energy axis, shared between the map and the XAS panel
        double y_link_min = 0.0, y_link_max = 1.0;

        int hovered_core  = -1;
        int selected_core = -1;

        // --- RIXS spectrum (a single incoming photon energy, i.e. one row of the map) ---
        // The spectrum is cheap enough to rebuild every frame, so nothing is cached here beyond the
        // selection state. The broadening settings above are shared with the map.
        int selected_photon = 0;   // index into the ascending photon energy list
        int hovered_final   = -1;
        int selected_final  = -1;
        uint64_t spectrum_hash = 0;  // only used to refit the axes when the plotted data changes

        // Set once the non-uniform photon spacing has been reported, to avoid spamming the log while
        // a settings slider is being dragged.
        bool warned_non_uniform = false;
    } rixs;

    // Settings + interaction state for the XPS plot.
    // No cached arrays: the broadened curve is generated on the fly by an ImPlotGetter and the peak
    // arrays are a few tens of doubles built into the frame temp arena, so there is nothing worth
    // keeping between frames.
    struct Xps {
        broadening_mode_t broadening_mode = BROADENING_MODE_LORENTZIAN;
        double broadening_fwhm_ev = 0.5;    // Matches the default of plot_xps_spectrum()

        bool show_sticks = true;
        bool invert_x    = false;   // XPS binding energy axes are often drawn decreasing to the right

        // false -> every core hole contributes a stick of height 1 (what VeloxChem's plot does)
        // true  -> stick height is the atom's share of the core MO
        // These differ only for delocalized holes, see the comment in draw_xps_plot().
        bool weight_by_contribution = false;

        // Which element is on screen. One at a time: a peak is now an ATOM rather than an element,
        // and core binding energies of different elements are hundreds of eV apart, so putting two
        // elements on one axis leaves each of them an unreadable spike. Zero means nothing has been
        // picked yet - the first group with entries is adopted on the first draw.
        md_element_t element = 0;

        // Index into the FLAT entry array, not into a group. -1 when the cursor is not on a peak.
        // There is deliberately no 'selected' counterpart: a peak draws as selected when its atom is
        // in the application selection mask, so the plot is a view of that mask and cannot drift
        // out of sync with the viewport or with a selection made anywhere else.
        int hovered = -1;

        // MO of the hovered peak, handed to the orbital grid window so it can highlight the matching
        // table row, or -1 when no peak is hovered. Written by draw_xps_plot, consumed and cleared by
        // draw_orb_window - see the comment at the top of that function for the ordering.
        int32_t highlight_mo_idx = -1;

        uint64_t hash = 0;  // Refit the axes when the settings change
    } xps;

    struct Vib {
        int hovered = -1;
        int selected = -1;

        // In case of resonance raman, this is the index of the selected electronic transition
        int external_frequency_index = 0;

        double broadening_fwhm = 20.0; // FWHM in cm^-1, default value is 20 cm^-1
        broadening_mode_t broadening_mode = BROADENING_MODE_LORENTZIAN;
        bool coord_modified = false;
        float displacement_amp_scl  = 1.0f;
        float displacement_freq_scl = 1.0f;
        bool invert_x = true;
        bool invert_y = false;
        bool displace_aos = false;

        // Time accumulation for vibrational mode visualization (pertubation of atoms)
        double t = 0;

        // Scaling applied to input frequencies (to account for inaccuracies of the chosen basis set)
        double freq_scaling_factor = 1.0;

        // Hash to track if the underlying peaks have changed
        uint64_t hash = 0;
    } vib;

    struct Export {
        char export_path[1024] = {};
        ElectronicStructureSource source = ElectronicStructureSource::MolecularOrbital;
        bool use_magnitude = false;
        ElectronicStructureSpin spin = ElectronicStructureSpin::Total;
        ElectronicStructureNtoComponent nto_component = ElectronicStructureNtoComponent::Particle;
        ElectronicStructureTransitionDensityComponent transition_density_component = ElectronicStructureTransitionDensityComponent::Attachment;
        VolumeResolution resolution = VolumeResolution::Mid;

        struct {
            md_vlx_spin_t type = MD_VLX_SPIN_ALPHA;
            int idx = 0;
        } mo;

        struct {
            md_vlx_nto_type_t type = MD_VLX_NTO_PARTICLE;
            int lambda_idx = 0;
            int idx = 1;
        } nto;

        bool use_obb = true;
        bool show_window = false;
    md_system_state_t export_state = {0};
    } export_state;

    // Arena for persistent allocations for the veloxchem module (tied to the lifetime of the VLX object)
    md_allocator_i* arena = 0;

    // Extracted GTO basis from the VLX file. One shell per contracted radial shell,
    // pure radial normalization coefficients (sph→cart factors baked in by gto.c).
    // Allocated from arena; zeroed when reset_data() is called.
    md_gto_basis_t basis = {};

    size_t num_molecular_orbitals() const {
        return md_vlx_scf_number_of_molecular_orbitals(vlx);
    }

    size_t num_natural_transition_orbitals() const {
        return md_vlx_rsp_number_of_excited_states(vlx);
    }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t event_idx = 0; event_idx < num_events; ++event_idx) {
            const viamd::Event& e = events[event_idx];

            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState& state = *(ApplicationState*)e.payload;

                arena = md_arena_allocator_create(state.allocator.persistent, MEGABYTES(4));

                // Reservation only, pages are committed on demand. The high-water mark is a single
                // RIXS map (num_photon_energies * num_grid_samples doubles, the latter capped at
                // 4096) plus the XAS panel, so this is generous by a wide margin.
                rixs.alloc = md_vm_arena_create(MEGABYTES(512));

                int gl_major, gl_minor;
                glGetIntegerv(GL_MAJOR_VERSION, &gl_major);
                glGetIntegerv(GL_MINOR_VERSION, &gl_minor);
                use_gpu_path = (!FORCE_CPU_PATH && gl_major >= 4 && gl_minor >= 3);
#if MD_ENABLE_GPU
                md_gto_gpu_initialize(state.gpu_device);
                md_topo_gpu_initialize(state.gpu_device);
#endif
                break;
            }
            case viamd::EventType_ViamdShutdown:
#if MD_ENABLE_GPU
                gl::pbo_upload_shutdown();
                md_gto_gpu_shutdown();
                md_topo_gpu_shutdown();
                // Nothing is freed here: the device scratch and the dataset's uploaded basis are
                // both the application's, and it destroys them.
#endif
                md_arena_allocator_destroy(arena);
                if (rixs.alloc) {
                    md_vm_arena_destroy(rixs.alloc);
                    rixs.alloc = nullptr;
                }
                break;
            case viamd::EventType_ViamdFrameTick: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState& state = *(ApplicationState*)e.payload;

                if (vlx) {
                    // Before any window draws, so both transition windows see the same numbers and
                    // neither has to be open for the other to have data.
                    update_nto_derived_data();

                    draw_orb_window(state);
                    draw_summary_window(state);
                    draw_rsp_window(state);
                    if (md_vlx_rsp_has_nto(vlx)) {
                        draw_nto_window(state);
                        draw_flow_window(state);
                    }
                    draw_export_window(state);
                }
                break;
            }
            case viamd::EventType_ViamdWindowDrawMenu:
                if (vlx) {
                    if (ImGui::BeginMenu("VeloxChem")) {
                        ImGui::Checkbox("Summary", &summary.show_window);
                        ImGui::Checkbox("Response", &rsp.show_window);
                        ImGui::Checkbox("Orbital Grid", &orb.show_window);
                        if (md_vlx_rsp_has_nto(vlx)) {
                            ImGui::Checkbox("Transition Analysis", &nto.show_window);
                            ImGui::Checkbox("Charge Transfer Analysis", &flow.show_window);
                        }
                        ImGui::Checkbox("Export", &export_state.show_window);
                        ImGui::EndMenu();
                    }
                }
                break;
            case viamd::EventType_ViamdRenderTransparent: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
				const ApplicationState& state = *(ApplicationState*)e.payload;

                md_temp_scope_t temp = md_temp_begin_in(state.allocator.frame);
                defer { md_temp_end(temp); };

                if (critical_points.enabled && critical_points.simp_graph.num_vertices > 0) {
                    immediate::Scope scope(state.gfx.overlay, "veloxchem_critical_points");
                    // Render topology as points if available
                    immediate::set_picking_base_idx(scope, critical_points.picking_range.beg);

                    glDisable(GL_DEPTH_TEST);

                    for (size_t i = 0; i < critical_points.simp_graph.num_edges; ++i) {
						uint32_t i0 = critical_points.simp_graph.edges[i].from;
						uint32_t i1 = critical_points.simp_graph.edges[i].to;
                        vec3_t p0 = { critical_points.simp_graph.vertices[i0].x,
                                      critical_points.simp_graph.vertices[i0].y,
									  critical_points.simp_graph.vertices[i0].z };
                        vec3_t p1 = { critical_points.simp_graph.vertices[i1].x,
									  critical_points.simp_graph.vertices[i1].y,
                                      critical_points.simp_graph.vertices[i1].z };
						p0 *= BOHR_TO_ANGSTROM;
						p1 *= BOHR_TO_ANGSTROM;
                        immediate::line(scope, p0, p1, immediate::COLOR_BLACK);
                    }

                    /*
                        MD_TOPO_UNDEFINED = 0,
                        MD_TOPO_MAXIMUM = 1,
                        MD_TOPO_SPLIT_SADDLE = 2,
                        MD_TOPO_MINIMUM = 3,
                        MD_TOPO_JOIN_SADDLE = 4,
                    */
                    const uint32_t type_colors[5] = {
                        0xFF7A7A7A, // unknown
                        0xFF5B5BD9, // maximum
                        0xFF5FAE5C, // split saddle
                        0xFFD49A4C, // minimum
                        0xFFB9A14A  // join saddle
                    };

                    size_t num_verts = critical_points.simp_graph.num_vertices;
                    immediate::Vertex* vertices = md_temp_alloc_array(temp, immediate::Vertex, num_verts);

                    vec4_t selection_color = state.selection.color.selection.visible;
                    vec4_t highlight_color = state.selection.color.highlight.visible;
                    highlight_color.w += (float)(sin(ImGui::GetTime() * HIGHLIGHT_PULSE_TIME_SCALE) * HIGHLIGHT_PULSE_ALPHA_SCALE);
                    highlight_color.w = CLAMP(highlight_color.w, 0.0f, 1.0f);

                    float saturation_scale = md_bitfield_popcount(&critical_points.selection_mask) > 0 ? state.selection.color.saturation : 1.0f;

                    for (size_t i = 0; i < num_verts; ++i) {
                        md_topo_critical_point_type_t type = md_topo_vertex_type(&critical_points.simp_graph, i);
                        vec4_t base_color = vec4_from_u32(type_colors[type]);

                        bool selected = md_bitfield_test_bit(&critical_points.selection_mask, i);
                        bool highlighted = md_bitfield_test_bit(&critical_points.highlight_mask, i);

                        vec4_t color = selected ? selection_color : tint_color(base_color, base_color, 0.0f, saturation_scale);
                        color = highlighted ? vec4_lerp(color, highlight_color, highlight_color.w) : color;

                        vertices[i].coord = {critical_points.simp_graph.vertices[i].x, critical_points.simp_graph.vertices[i].y, critical_points.simp_graph.vertices[i].z};
                        vertices[i].coord *= BOHR_TO_ANGSTROM;
                        vertices[i].color = u32_from_vec4(color);
                        vertices[i].picking_idx = (uint32_t)i;
                    }

                    // Sort on view Z to improve picking, since we cannot use depth testing
                    const mat4_t view = camera_world_to_view_matrix(state.view.camera);
                    std::sort(vertices, vertices + num_verts, [&view](const immediate::Vertex& a, const immediate::Vertex& b) {
                        float za = view[0][2] * a.coord.x + view[1][2] * a.coord.y + view[2][2] * a.coord.z + view[3][2];
                        float zb = view[0][2] * b.coord.x + view[1][2] * b.coord.y + view[2][2] * b.coord.z + view[3][2];
                        return za < zb;
                    });

                    immediate::points(scope, vertices, num_verts);
                }

                //ApplicationState& state = *(ApplicationState*)e.payload;
                //draw_orb_volume(state);
                break;
            }
            case viamd::EventType_ViamdInteractionSurface: {
                ASSERT(e.payload_type == viamd::EventPayloadType_InteractionSurfaceEvent);
                ASSERT(e.payload);
                InteractionSurfaceEvent* event = (InteractionSurfaceEvent*)e.payload;
                if (event->surface_id == interaction_surface_main) {
                    if (critical_points.enabled && md_bitfield_validate(&critical_points.highlight_mask)) {
                        switch (event->kind) {
                        case InteractionSurfaceEventKind::Hover: [[fallthrough]];
                        case InteractionSurfaceEventKind::Click: {
                            md_bitfield_clear(&critical_points.highlight_mask);
                            if (event->hit.domain == PickingDomain_CriticalPoints) {
                                md_bitfield_set_bit(&critical_points.highlight_mask, event->hit.local_idx);
                            }
                            if (event->selection_mode == InteractionSelectionMode::Append) {
                                md_bitfield_or_inplace(&critical_points.selection_mask, &critical_points.highlight_mask);
                            } else if (event->selection_mode == InteractionSelectionMode::Remove) {
                                if (event->hit.domain == PickingDomain_CriticalPoints) {
                                    md_bitfield_andnot_inplace(&critical_points.selection_mask, &critical_points.highlight_mask);
                                } else if (event->hit.domain == 0) {
                                    md_bitfield_clear(&critical_points.selection_mask);
                                }
                            }
                            break;
                        }
                        case InteractionSurfaceEventKind::RegionSelect: {
                            md_bitfield_clear(&critical_points.highlight_mask);
                            if (event->selection_mode == InteractionSelectionMode::Remove) {
                                // When removing, only consider currently selected items as candidates for region selection
                                md_bitfield_iter_t it = md_bitfield_iter_create(&critical_points.selection_mask);
                                while (md_bitfield_iter_next(&it)) {
                                    size_t idx = md_bitfield_iter_idx(&it);
                                    vec3_t xyz = { critical_points.simp_graph.vertices[idx].x,
                                                   critical_points.simp_graph.vertices[idx].y,
                                                   critical_points.simp_graph.vertices[idx].z };
                                    xyz *= BOHR_TO_ANGSTROM;
                                    vec2_t coord = world_to_surface_project(xyz, event->world_to_clip, event->surface_size);
                                    if (coord.x >= event->region_min.x && coord.x <= event->region_max.x &&
                                        coord.y >= event->region_min.y && coord.y <= event->region_max.y)
                                    {
                                        md_bitfield_set_bit(&critical_points.highlight_mask, idx);
                                    }
                                }
                            } else {
                                for (size_t i = 0; i < critical_points.simp_graph.num_vertices; ++i) {
                                    vec3_t xyz = { critical_points.simp_graph.vertices[i].x,
                                                   critical_points.simp_graph.vertices[i].y,
                                                   critical_points.simp_graph.vertices[i].z };
                                    xyz *= BOHR_TO_ANGSTROM;
                                    vec2_t coord = world_to_surface_project(xyz, event->world_to_clip, event->surface_size);
                                    if (coord.x >= event->region_min.x && coord.x <= event->region_max.x &&
                                        coord.y >= event->region_min.y && coord.y <= event->region_max.y)
                                    {
                                        md_bitfield_set_bit(&critical_points.highlight_mask, i);
                                    }
                                }
                            }

                            if (event->region_phase == InteractionSurfaceEventPhase::Commit) {
                                // Merge highlight into selection
                                if (event->selection_mode == InteractionSelectionMode::Append) {
                                    md_bitfield_or_inplace(&critical_points.selection_mask, &critical_points.highlight_mask);
                                }
                                else if (event->selection_mode == InteractionSelectionMode::Remove) {
                                    md_bitfield_andnot_inplace(&critical_points.selection_mask, &critical_points.highlight_mask);
                                }
                                md_bitfield_clear(&critical_points.highlight_mask);
                            }
                            break;
                        }
                        default:
                            break;
                        }
                    }
                }
                break;
            }
            case viamd::EventType_ViamdPickingTooltipTextRequest: {
                ASSERT(e.payload_type == viamd::EventPayloadType_PickingTooltipTextRequest);
                ASSERT(e.payload);
                PickingTooltipTextRequest* req = (PickingTooltipTextRequest*)e.payload;
                if (req->hit.domain == PickingDomain_CriticalPoints) {
                    uint32_t cp_idx = req->hit.local_idx;
                    if (cp_idx < critical_points.simp_graph.num_vertices) {
                        md_topo_critical_point_type_t type = md_topo_vertex_type(&critical_points.simp_graph, cp_idx);
                        const char* str = md_topo_critical_point_type_str(type);
                        float value = critical_points.simp_graph.vertices[cp_idx].value;
                        md_strb_fmt(&req->sb, "Type: %s\nValue: %.7f", str, value);
                    }
                }
                break;
            }
#if 0
            case viamd::EventType_ViamdSystemInit: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState& state = *(ApplicationState*)e.payload;
                init_from_file(str_from_cstr(state.files.molecule), state);
                break;
            }
#endif
            case viamd::EventType_ViamdLoadData: {
                 ASSERT(e.payload_type == viamd::EventPayloadType_LoadData);
                 LoadDataPayload& payload = *(LoadDataPayload*)e.payload;
                 if (payload.loader_state.type == LoaderType_VLX_H5) {
                     init_from_file(payload.path_to_file, *payload.app_state);
                 }
                 break;
            }
            case viamd::EventType_ViamdSystemFree: {
                reset_data();
                break;
            }
            case viamd::EventType_ViamdSerialize: {
                ASSERT(e.payload_type == viamd::EventPayloadType_SerializationState);
                serialize_workspace(*(viamd::serialization_state_t*)e.payload);
                break;
            }
            case viamd::EventType_ViamdDeserialize: {
                ASSERT(e.payload_type == viamd::EventPayloadType_DeserializationState);
                viamd::deserialization_state_t& deser = *(viamd::deserialization_state_t*)e.payload;
                if (str_eq(viamd::section_header(deser), STR_LIT("VeloxChem"))) {
                    deserialize_workspace(deser);
                }
                break;
            }
            case viamd::EventType_ViamdSystemStateChanged: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState& state = *(ApplicationState*)e.payload;

                if (vlx) {
                    // Update atom_xyz
					// Calculate nuclei charge weighted center of charge for later use in orbital centering
					dvec3_t nucl_dipole = {0};

					size_t count = qm_to_atom_idx ? md_vlx_number_of_atoms(vlx) : state.mold.sys.atom.count;
                    ASSERT(md_array_size(atom_xyzw) == count);
                    for (size_t i = 0; i < count; ++i) {
                        int idx = qm_to_atom_idx ? qm_to_atom_idx[i] : (int)i;
                        atom_xyzw[i].x = (float)(state.mold.state.x[idx] * ANGSTROM_TO_BOHR);
                        atom_xyzw[i].y = (float)(state.mold.state.y[idx] * ANGSTROM_TO_BOHR);
                        atom_xyzw[i].z = (float)(state.mold.state.z[idx] * ANGSTROM_TO_BOHR);

						int z = md_atom_atomic_number(&state.mold.sys.atom, idx);
						nucl_dipole.x += atom_xyzw[i].x * z;
						nucl_dipole.y += atom_xyzw[i].y * z;
						nucl_dipole.z += atom_xyzw[i].z * z;
                    }

                    size_t num_electrons = md_vlx_number_of_electrons(vlx, MD_VLX_SPIN_ALPHA) + md_vlx_number_of_electrons(vlx, MD_VLX_SPIN_BETA);
                    nucl_dipole /= num_electrons > 0 ? (double)num_electrons : 1.0;
                    center_of_charge = nucl_dipole - md_vlx_scf_ground_state_dipole_moment(vlx);
                    
                    // Recalculate OABB and AABB
                    oabb.orientation = mat3_PCA(atom_xyzw, md_array_size(atom_xyzw));
                    calculate_bounds(oabb.min_ext.elem, oabb.max_ext.elem, atom_xyzw, md_array_size(atom_xyzw), oabb.orientation);
                    calculate_bounds(aabb.min_ext.elem, aabb.max_ext.elem, atom_xyzw, md_array_size(atom_xyzw));
                }

                for (size_t i = 0; i < md_array_size(state.representation.reps); ++i) {
                    if (state.representation.reps[i].type == RepresentationType::ElectronicStructure) {
                        flag_representation_as_dirty(&state.representation.reps[i]);
                    }
                }
                
                break;
            }
            case viamd::EventType_ViamdRepresentationInfoFill: {
                ASSERT(e.payload_type == viamd::EventPayloadType_RepresentationInfo);
                RepresentationInfo& info = *(RepresentationInfo*)e.payload;

                info.alpha.homo_idx = homo_idx[0];
                info.alpha.lumo_idx = lumo_idx[0];

                info.beta.homo_idx  = homo_idx[1];
                info.beta.lumo_idx  = lumo_idx[1];

                size_t num_mos = num_molecular_orbitals();
                if (num_mos) {
                    const double* occ = md_vlx_scf_mo_occupancy(vlx, MD_VLX_SPIN_ALPHA);
                    const double* ene = md_vlx_scf_mo_energy   (vlx, MD_VLX_SPIN_ALPHA);

                    if (occ && ene) {
                        info.alpha.num_orbitals = num_mos;
                        md_array_resize(info.alpha.label,       num_mos, info.alloc);
                        md_array_resize(info.alpha.occupation,  num_mos, info.alloc);
                        md_array_resize(info.alpha.energy,      num_mos, info.alloc);

                        for (size_t i = 0; i < num_mos; ++i) {
                            const char* lbl = "";
                            if (i == info.alpha.homo_idx) {
                                lbl = " (homo)";
                            } else if (i == info.alpha.lumo_idx) {
                                lbl = " (lumo)";
                            }
                            info.alpha.label[i]      = str_printf(info.alloc, "%zu%s", i + 1, lbl);
                            info.alpha.energy[i]     = ene[i];
                            info.alpha.occupation[i] = occ[i];
                        }
                    }
                }

                // @TODO: Check condition of wheter or not to include beta orbitals
                if (true) {
                    const double* occ = md_vlx_scf_mo_occupancy(vlx, MD_VLX_SPIN_BETA);
                    const double* ene = md_vlx_scf_mo_energy   (vlx, MD_VLX_SPIN_BETA);

                    if (occ && ene) {
                        info.beta.num_orbitals = num_mos;
                        md_array_resize(info.beta.label,       num_mos, info.alloc);
                        md_array_resize(info.beta.occupation,  num_mos, info.alloc);
                        md_array_resize(info.beta.energy,      num_mos, info.alloc);

                        for (size_t i = 0; i < num_mos; ++i) {
                            const char* lbl = "";
                            if (i == info.beta.homo_idx) {
                                lbl = " (homo)";
                            } else if (i == info.beta.lumo_idx) {
                                lbl = " (lumo)";
                            }
                            info.beta.label[i]      = str_printf(info.alloc, "%zu%s", i + 1, lbl);
                            info.beta.energy[i]     = ene[i];
                            info.beta.occupation[i] = occ[i];
                        }
                    }
                }
                
                size_t num_excited_states = md_vlx_rsp_number_of_excited_states(vlx);
                if (md_vlx_rsp_has_nto(vlx) && num_excited_states > 0) {
                    info.nto.num_orbitals = num_excited_states;
                    md_array_resize(info.nto.label, num_excited_states, info.alloc);
                    for (size_t i = 0; i < num_excited_states; ++i) {
                        info.nto.label[i] = str_printf(info.alloc, "%zu", i + 1);

                        const double LAMBDA_CUTOFF = 1.0e-3;
                        double lambdas[16] = {0};
                        size_t lambda_count = md_vlx_rsp_nto_lambdas_extract(lambdas, vlx, i, ARRAY_SIZE(lambdas));
                        NaturalTransitionOrbitalLambda lambda_info = {};
                        for (size_t j = 0; j < lambda_count; ++j) {
                            if (lambdas[j] < LAMBDA_CUTOFF) break;
                            str_t lbl = str_printf(info.alloc, (const char*)u8"λ[%zu] (%.3f)", j + 1, lambdas[j]);
                            md_array_push(lambda_info.label, lbl, info.alloc);
                            md_array_push(lambda_info.value, lambdas[j], info.alloc);
                            lambda_info.num_lambdas += 1;
                        }
                        md_array_push(info.nto.lambda, lambda_info, info.alloc);
                    }
                }
                
                // Neither the source mask nor the density property list is filled here any more.
                // Both are answers about what the SYSTEM holds, and update_representation_info reads
                // them straight off its attribute table - which is the exact test, and one a second
                // reader satisfies without touching this file.
                //
                // Atomic properties and dipoles went the same way earlier: published into
                // sys->attributes at load, gathered from there, no round trip back into here.
                break;
            }
            case viamd::EventType_ViamdPickingRangeReserve: {
                ASSERT(e.payload_type == viamd::EventPayloadType_PickingSpace);
                PickingSpace* space = (PickingSpace*)e.payload;
                if (critical_points.enabled) {
                    size_t num_cp = md_topo_num_critical_points(&critical_points.simp_graph);
                    picking_range_reserve(&critical_points.picking_range, space, PickingDomain_CriticalPoints, num_cp);
                }
                break;
            }
            case viamd::EventType_ViamdViewFit: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ViewFitRequest);
                ViewFitRequest* req = (ViewFitRequest*)e.payload;
                if (req && req->surface_id == interaction_surface_main) {
                    md_bitfield_t* bf = nullptr;
                    switch (req->round) {
                    case ViewFitRound_Highlight:
                        bf = &critical_points.highlight_mask; break;
                    case ViewFitRound_Selection:
                        bf = &critical_points.selection_mask; break;
					case ViewFitRound_Visible: [[fallthrough]];
                    default:
                        break;
                    }

                    // Add highlighted atoms to the mask used for view fitting, so that they are included in the computed optimal view
                    if (critical_points.enabled) {
                        size_t popcount = bf ? md_bitfield_popcount(bf) : 0;
                        if (popcount > 0) {
                            vec4_t* dst_xyzw = md_array_extend(req->xyzw, popcount, req->alloc);
                            if (dst_xyzw) {
                                md_bitfield_iter_t it = md_bitfield_iter_create(bf);
                                while (md_bitfield_iter_next(&it)) {
                                    size_t idx = md_bitfield_iter_idx(&it);
                                    dst_xyzw->x = (float)(critical_points.simp_graph.vertices[idx].x * BOHR_TO_ANGSTROM);
                                    dst_xyzw->y = (float)(critical_points.simp_graph.vertices[idx].y * BOHR_TO_ANGSTROM);
                                    dst_xyzw->z = (float)(critical_points.simp_graph.vertices[idx].z * BOHR_TO_ANGSTROM);
                                    dst_xyzw->w = 1.0f;
                                    dst_xyzw++;
                                }
                            }
                        }
                    }
                }
                break;
            }
            default:
                break;
            }
        }
    }

    void reset_data() {
#if MD_ENABLE_GPU
        // The topo context is the only GPU object this component owns. The uploaded basis goes with
        // the system it was built from, released by free_system_data.
        md_topo_gpu_context_destroy(critical_points.topo_ctx);
#endif
        //md_gl_mol_destroy(gl_mol);
        md_gl_rep_destroy(gl_rep);
        md_vlx_destroy(vlx);
        md_arena_allocator_reset(arena);
        md_topo_extremum_graph_free(&critical_points.raw_graph);
        md_topo_extremum_graph_free(&critical_points.simp_graph);
        critical_points = {};
        critical_points.raw_graph.alloc = arena;
        critical_points.simp_graph.alloc = arena;
        basis = {};  // arena reset above invalidates allocations; zero the struct
        vlx = nullptr;
        atom_xyzw = nullptr;
        orb = VeloxChem::Orb{};
        nto = VeloxChem::Nto{};
        // The arena reset above already released everything these held.
        flow = VeloxChem::Flow{};
        rsp = VeloxChem::Rsp{};
        // The rixs arena outlives the VLX object, so carry it across the value-initialization.
        // Rewinding it first drops the cached arrays that the zeroed struct would otherwise orphan.
        md_allocator_i* rixs_alloc = rixs.alloc;
        reset_rixs_cache(rixs);
        rixs = VeloxChem::Rixs{};
        rixs.alloc = rixs_alloc;
        vib = VeloxChem::Vib{};
        opt = VeloxChem::Opt{};
    }

    void update_nto_group_colors() {
        md_temp_scope_t temp = md_temp_begin();
		defer { md_temp_end(temp); };
        size_t num_atoms = md_vlx_number_of_atoms(vlx);
        uint32_t* colors = (uint32_t*)md_temp_alloc_array(temp, uint32_t, num_atoms);
        // Group 0 is "not in any group", so it keeps the element colour the molecule would have
        // anyway. Stamping the unassigned grey over everything would make an ungrouped system look
        // like a system whose groups all happen to be grey.
        const md_element_t* atom_z = md_vlx_atomic_numbers(vlx);
        for (size_t i = 0; i < num_atoms; ++i) {
            const size_t group_idx = nto.atom_group_idx[i];
            uint32_t color = U32_MAGENTA;
            if (group_idx == 0) {
                color = atom_z ? md_util_element_cpk_color(atom_z[i]) : U32_MAGENTA;
            } else if (group_idx < nto.group.count) {
                color = convert_color(nto.group.color[group_idx]);
            }
            colors[i] = color;
        }

        md_gl_rep_set_atom_colors(nto.gl_rep, 0, (uint32_t)num_atoms, colors, 0);
    }

    void init_from_file(str_t filename, ApplicationState& state) {
        str_t ext;
        if (extract_ext(&ext, filename)) {
            if (str_eq_ignore_case(ext, STR_LIT("h5"))) {
                MD_LOG_INFO("Attempting to load VeloxChem data from file '" STR_FMT "'", STR_ARG(filename));
                
                if (!vlx) {
                    vlx = md_vlx_create(arena);
                } else {
                    md_vlx_reset(vlx);
                }
                
                if (md_vlx_parse_file(vlx, filename)) {
                    // Publish first, then read back. Everything this file carries which fits the
                    // attribute model goes into the system's table here - the dipole groups, the
                    // SCF history, orbital energies, response and vibrational quantities, normal
                    // modes, the GTO basis and the MO coefficients. mdlib owns the paths, formats
                    // and units; this component is a consumer of them from here on, like anything
                    // else which can read a system.
                    md_vlx_publish_attributes(&state.mold.sys, vlx);

                    // And the basis comes back out of that table rather than out of the vlx object,
                    // which is the whole point: an evaluator needs the system, not the reader.
                    md_gto_basis_free(&basis, arena);
                    if (!md_gto_basis_extract_attributes(&basis, &state.mold.sys.attributes, arena)) {
                        MD_LOG_ERROR("Failed to build a GTO basis from the system's attributes");
                        return;
                    }
                    size_t num_vlx_atoms = md_vlx_number_of_atoms(vlx);

                    const int* local_to_global = md_vlx_local_to_global_atom_idx(vlx);
                    if (local_to_global) {
                        // This means that the vlx file is derived from a larger system and only contains data for a subset.
                        // We can support this by remapping atom indices.
                        // *HOWEVER* The VLX h5 file can be loaded as a stand alone system and in such case, the mapping should not be applied.
                        if (num_vlx_atoms < state.mold.sys.atom.count) {
                            // @TODO: Check atomic numbers consistency as well
                            const uint8_t* atomic_numbers = md_vlx_atomic_numbers(vlx);
                            bool match = true;
                            for (size_t i = 0; i < num_vlx_atoms; ++i) {
                                int idx = local_to_global[i];
                                if (md_atom_atomic_number(&state.mold.sys.atom, idx) != atomic_numbers[i]) {
                                   match = false;
                                   break;
                                }
                            }
                            if (match) {
                                qm_to_atom_idx = local_to_global;
                            } else {
                                MD_LOG_ERROR("VLX file '%s' contains a subset of atoms, but their atomic numbers do not match with the system. Atom index remapping will not be applied.", STR_ARG(filename));
                                md_vlx_reset(vlx);
                                return;
                            }
                        }
                    }

                    // Which system atom each basis atom is, published beside the basis so that an
                    // evaluation can gather positions with nothing but the system in reach. NULL is
                    // the standalone case and clears any map a previous load left behind.
                    basis_atom_map_publish(&state.mold.sys, qm_to_atom_idx, num_vlx_atoms);

                    homo_idx[0] = (int)md_vlx_scf_homo_idx(vlx, MD_VLX_SPIN_ALPHA);
                    homo_idx[1] = (int)md_vlx_scf_homo_idx(vlx, MD_VLX_SPIN_BETA);

                    lumo_idx[0] = (int)md_vlx_scf_lumo_idx(vlx, MD_VLX_SPIN_ALPHA);
                    lumo_idx[1] = (int)md_vlx_scf_lumo_idx(vlx, MD_VLX_SPIN_BETA);

                    size_t num_colors = state.mold.sys.atom.count;
                    uint32_t* colors = (uint32_t*)md_vm_arena_push(state.allocator.frame, num_colors * sizeof(uint32_t));
                    color_atoms_type(colors, num_colors, state.mold.sys);

                    if (qm_to_atom_idx) {
                        md_bitfield_t mask = md_bitfield_create(state.allocator.frame);
                        md_bitfield_set_indices_u32(&mask, (uint32_t*)qm_to_atom_idx, md_vlx_number_of_atoms(vlx));
                        filter_colors(colors, num_colors, &mask);
                        reset_view(&default_view, state.mold.state, &mask);
                    } else {
                        reset_view(&default_view, state.mold.state);
                    }

                    gl_rep = md_gl_rep_create(state.mold.gl_mol);
                    md_gl_rep_set_atom_colors(gl_rep, 0, (uint32_t)num_colors, colors, 0);

					const md_element_t* atom_z = md_vlx_atomic_numbers(vlx);
					dvec3_t nucl_dipole = { 0, 0, 0 };

                    md_array_resize(atom_xyzw, num_vlx_atoms, arena);
                    const dvec3_t* atom_vlx = md_vlx_atom_coordinates(vlx);
                    for (size_t i = 0; i < num_vlx_atoms; ++i) {
                        dvec3_t xyz = atom_vlx[i] * ANGSTROM_TO_BOHR;
                        atom_xyzw[i] = vec4_set((float)xyz.x, (float)xyz.y, (float)xyz.z, 1.0f);
						nucl_dipole += xyz * atom_z[i];
                    }

					dvec3_t ground_state_dipole = md_vlx_scf_ground_state_dipole_moment(vlx);

                    size_t num_electrons = md_vlx_number_of_electrons(vlx, MD_VLX_SPIN_ALPHA) + md_vlx_number_of_electrons(vlx, MD_VLX_SPIN_BETA);
                    const double inv_ne = num_electrons > 0 ? 1.0 / (double)num_electrons : 1.0;
                    center_of_charge = (nucl_dipole - ground_state_dipole) * inv_ne;

                    oabb.orientation = mat3_PCA(atom_xyzw, md_array_size(atom_xyzw));
                    calculate_bounds(oabb.min_ext.elem, oabb.max_ext.elem, atom_xyzw, md_array_size(atom_xyzw), oabb.orientation);
                    calculate_bounds(aabb.min_ext.elem, aabb.max_ext.elem, atom_xyzw, md_array_size(atom_xyzw));

#if MD_ENABLE_GPU
                    // The dataset's uploaded basis is built from the basis/ attributes just
                    // published, by the application rather than here: it is derived from the
                    // system, every consumer of that system wants the same one, and this component
                    // is only the first to ask.
                    system_gpu_data_update(&state, DEFAULT_GTO_CUTOFF_VALUE);
#endif

                    // NTO

                    size_t num_excited_states = md_vlx_rsp_number_of_excited_states(vlx);
                    if (num_excited_states > 0) {
                        //nto.show_window = true;

                        picking_surface_init(&nto.picking_surface, interaction_surface_nto);
						vec3_t center = (oabb.min_ext + oabb.max_ext) * 0.5f * BOHR_TO_ANGSTROM;
						vec3_t half_ext = (oabb.max_ext - oabb.min_ext) * 0.5f * BOHR_TO_ANGSTROM;
                        nto.target = compute_optimal_view(center, half_ext, oabb.orientation, nto.distance_scale);
                        nto.camera = nto.target;

						nto.atom_group_idx = md_array_create(uint32_t, num_vlx_atoms, arena);
                        MEMSET(nto.atom_group_idx, 0, sizeof(uint32_t) * num_vlx_atoms);

                        snprintf(nto.group.label[0], sizeof(nto.group.label[0]), "Unassigned");
                        nto.group.color[0] = vec4_t{ 0.25f, 0.25f, 0.25f, 1.0f };

                        for (int i = 1; i < (int)ARRAY_SIZE(nto.group.color); ++i) {
                            ImVec4 color = ImPlot::GetColormapColor(i - 1, ImPlotColormap_Deep);
                            nto.group.color[i] = vec_cast(color);
                            snprintf(nto.group.label[i], sizeof(nto.group.label[i]), "Group %i", i);
                        }

                        // No groups until the user makes one. Slot 0 is the "everything else"
                        // bucket every atom starts in; it is not a group the user chose and is
                        // never drawn as one. Splitting the molecule in half, or stamping a
                        // hardcoded mapping onto a file with a particular name, put invented
                        // structure in front of the user and then invited them to read charge
                        // transfer off it.
                        nto.group.count = 1;
                        nto.gl_rep = md_gl_rep_create(state.mold.gl_mol);
                        update_nto_group_colors();

                        // Callculate ballpark scaling factor for dipole vectors
                        vec3_t extent = aabb.max_ext - aabb.min_ext;
                        float max_ext = MAX(extent.x, MAX(extent.y, extent.z));
                        float max_len = 0;
                        const dvec3_t* electric_dp = md_vlx_rsp_electric_transition_dipole_moments(vlx);
                        const dvec3_t* magnetic_dp = md_vlx_rsp_magnetic_transition_dipole_moments(vlx);
                        if (electric_dp && magnetic_dp) {
                            for (size_t i = 0; i < num_excited_states; ++i) {
                                max_len = MAX(max_len, (float)dvec3_length(electric_dp[i]));
                                max_len = MAX(max_len, (float)dvec3_length(magnetic_dp[i]));
                            }
                            nto.dipole.vector_scale = CLAMP((max_ext * 0.75f) / max_len, 0.1f, 10.0f);
                        }

#if USE_AABB_GRID
                        init_grid(&nto.grid, mat3_ident(), aabb.min_ext, aabb.max_ext, DEFAULT_SAMPLES_PER_ANGSTROM * BOHR_TO_ANGSTROM);
#else
                        init_grid(&nto.grid, oabb.orientation, oabb.min_ext, oabb.max_ext, DEFAULT_SAMPLES_PER_ANGSTROM * BOHR_TO_ANGSTROM);
#endif
                        nto.target = default_view;
                        nto.camera = nto.target;
                    }

                    // RSP
					md_vlx_rsp_type_t rsp_type = md_vlx_rsp_type(vlx);
                    if (rsp_type != MD_VLX_RSP_UNKNOWN) {
                        if (rsp_type == MD_VLX_RSP_LINEAR) {
                            rsp.hovered = -1;
                            rsp.selected = 0;
						}
                    }

                    // OPT
                    if (md_vlx_opt_number_of_steps(vlx) > 0) {
                        opt.selected = (int)(md_vlx_opt_number_of_steps(vlx) - 1);
                    }

                    // ORB
                    //orb.show_window = true;
                    picking_surface_init(&orb.picking_surface, interaction_surface_orb);
                    orb.target = default_view;
                    orb.camera = orb.target;
                    orb.mo_idx = homo_idx[0];
                    orb.scroll_to_idx = homo_idx[0];

                    // Export
                    export_state.mo.idx = homo_idx[0];

                    // CriPoAl
                    md_bitfield_init(&critical_points.selection_mask, arena);
                    md_bitfield_init(&critical_points.highlight_mask, arena);

                    // Last, so it overrides the placeholder grouping that init sets up.
                    apply_pending_groups();

                    MD_LOG_INFO("Successfully initialized VeloxChem data");
                }
                else {
                    MD_LOG_INFO("Failed to initialize VeloxChem data");
                    reset_data();
                }
            }
        }
    }

    // Every one of these forwards to the application's evaluator. What used to live here - packing
    // coefficients for the device, launching the kernel, queueing the readback, and the GL fallback
    // beside it - was never specific to VeloxChem; it is the same work for any system that
    // publishes a basis, and it now lives once in viamd.cpp where the GPU resources do.
    //
    // What is left is the translation from this component's vocabulary to an attribute path and a
    // slice, which is the only part that was ever ours.

    bool evaluate_nto(ApplicationState& state, uint32_t vol_tex, const md_grid_t& grid, size_t nto_idx, size_t lambda_idx, md_vlx_nto_type_t type, md_gto_op_t op) {
        const str_t path = (type == MD_VLX_NTO_PARTICLE) ? es_path::nto_particle : es_path::nto_hole;
        const md_attribute_slice_t slice = md_attribute_slice_2((uint32_t)nto_idx, (uint32_t)lambda_idx);
        return orbital_evaluate(&state, vol_tex, grid, path, &slice, MD_GTO_EVAL_MODE_PSI, op, DEFAULT_GTO_CUTOFF_VALUE);
    }

    bool evaluate_mo(ApplicationState& state, uint32_t vol_tex, const md_grid_t& grid, md_vlx_spin_t mo_type, size_t mo_idx, md_gto_op_t op) {
        const str_t path = (mo_type == MD_VLX_SPIN_BETA) ? es_path::beta_coefficient : es_path::alpha_coefficient;
        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)mo_idx);
        return orbital_evaluate(&state, vol_tex, grid, path, &slice, MD_GTO_EVAL_MODE_PSI, op, DEFAULT_GTO_CUTOFF_VALUE);
    }

    bool evaluate_transition_density(ApplicationState& state, uint32_t vol_tex, const md_grid_t& grid, size_t state_idx, ElectronicStructureTransitionDensityComponent component) {
        str_t path = es_path::attachment_density;
        switch (component) {
        case ElectronicStructureTransitionDensityComponent::Detachment: path = es_path::detachment_density; break;
        case ElectronicStructureTransitionDensityComponent::Difference: path = es_path::transition_diff;    break;
        case ElectronicStructureTransitionDensityComponent::Attachment:
        default: break;
        }
        // The attribute is VIRTUAL and {S,A,A}: the slice is what tells the provider to reconstruct
        // this one state instead of every state, which is the only reason asking is affordable.
        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)state_idx);
        return density_evaluate(&state, vol_tex, grid, path, &slice, MD_GTO_OP_SET);
    }

    bool evaluate_electron_density(ApplicationState& state, uint32_t vol_tex, const md_grid_t& grid, ElectronicStructureSpin spin, md_gto_op_t op) {
        return density_evaluate(&state, vol_tex, grid, es_electron_density_path(spin), nullptr, op);
    }

    static inline ImVec4 make_highlight_color(const ImVec4& color, float factor = 0.2f) {
        // Ensure the factor is not too high, to avoid over-brightening
        factor = CLAMP(factor, 0.0f, 1.0f);

        // Increase the brightness of each component by the factor
        float r = color.x + factor;
        float g = color.y + factor;
        float b = color.z + factor;

        // Clamp the values to the range [0, 1]
        r = ImClamp(r, 0.0f, 1.0f);
        g = ImClamp(g, 0.0f, 1.0f);
        b = ImClamp(b, 0.0f, 1.0f);

        // Return the new highlighted color with the same alpha
        return ImVec4(r, g, b, color.w);
    }

    static inline vec2_t vec2_from_imvec2(ImVec2 v) {
        return { v.x, v.y };
    }

    static inline md_vg_color_t vg_color_from_imvec4(const ImVec4& color) {
        return md_vg_color_from_vec4({ color.x, color.y, color.z, color.w });
    }

    static inline ImVec2 aligned_text_pos(const char* text, ImVec2 pos, ImVec2 alignment = { 0,0 }) {
        return pos - ImGui::CalcTextSize(text) * alignment;
    }

    static inline void compute_vertical_sankey_bezier(ImVec2* p1, ImVec2* p2, ImVec2* p3, ImVec2* p4, ImVec2 source_pos, ImVec2 dest_pos, float thickness) {
        ASSERT(p1);
        ASSERT(p2);
        ASSERT(p3);
        ASSERT(p4);

        *p1 = ImVec2(source_pos.x + thickness * 0.5f, source_pos.y);
        *p4 = ImVec2(dest_pos.x + thickness * 0.5f, dest_pos.y);

        const float dist = fabsf(p1->y - p4->y);
        const float curve_offset = dist * 0.25f;

        *p2 = ImVec2(source_pos.x + thickness * 0.5f, source_pos.y - curve_offset);
        *p3 = ImVec2(dest_pos.x + thickness * 0.5f, dest_pos.y + curve_offset);
    }

    // Draws a vertical flow based on cubic bezier and returns the abs delta from the mouse cursor to this line
    static inline ImVec2 draw_vertical_sankey_flow(ImDrawList* draw_list, ImVec2 source_pos, ImVec2 dest_pos, float thickness, ImU32 flow_color) {
        ImVec2 p1 = {};
        ImVec2 p2 = {};
        ImVec2 p3 = {};
        ImVec2 p4 = {};
        compute_vertical_sankey_bezier(&p1, &p2, &p3, &p4, source_pos, dest_pos, thickness);

        const int num_segments = 50;

        // Define the color and thickness for the flow
        //ImU32 flow_color = IM_COL32(100, 149, 237, 255);  // Cornflower Blue
        // ImBezierCubicCalc Use this for calculating mouse distance to curve
        // Draw the Bezier curve representing the flow
        draw_list->AddBezierCubic(p1, p2, p3, p4, flow_color, thickness, num_segments);

        ImVec2 mouse_pos = ImGui::GetMousePos();
        float clamped_y = ImClamp(mouse_pos.y, ImMin(p1.y, p4.y), ImMax(p1.y, p4.y));

        ImVec2 point_on_bezier = ImBezierCubicClosestPointCasteljau(p1, p2, p3, p4, mouse_pos, 1.0e-3f);
        ImVec2 delta = mouse_pos - point_on_bezier;
        return {ImAbs(delta.x), ImAbs(clamped_y - mouse_pos.y)};
    }

    static inline void draw_aligned_text(ImDrawList* draw_list, const char* text, ImVec2 pos, ImVec2 alignment = { 0,0 }) {
        draw_list->AddText(aligned_text_pos(text, pos, alignment), ImGui::ColorConvertFloat4ToU32({ 0,0,0,1 }), text);
    }

    static inline void vg_sankey_diagram(md_vg_scene_t* scene, ImRect area, Nto* nto, bool hide_text_overlaps) {
        if (!scene) return;
        if (!nto->transition_density_hole) return;
        if (!nto->transition_density_part) return;

        md_temp_scope_t temp = md_temp_begin_avoid(scene->alloc);
        defer { md_temp_end(temp); };

        md_vg_add_rect_filled(scene, vec2_from_imvec2(area.Min), vec2_from_imvec2(area.Max), MD_VG_COLOR_RGB(255, 255, 255), 0.0f);

        ImRect plot_area = area;
        const float plot_percent = 0.8f;
        plot_area.Expand({ area.GetWidth() * -(1 - plot_percent), area.GetHeight() * -(1 - plot_percent) });

        const float bar_height = plot_area.GetHeight() * 0.05f;
        int num_bars = 0;
        for (size_t i = 0; i < nto->group.count; ++i) {
            if (nto->transition_density_hole[i] != 0.0f) {
                num_bars += 1;
            }
        }
        const int num_gaps = MAX(num_bars - 1, 0);
        const float gap_size = plot_area.GetWidth() * 0.05f;
        const float bars_avail_width = plot_area.GetWidth() - gap_size * num_gaps;

        float hole_sum = 0.0f;
        float part_sum = 0.0f;
        for (size_t i = 0; i < nto->group.count; ++i) {
            hole_sum += nto->transition_density_hole[i];
            part_sum += nto->transition_density_part[i];
        }
        hole_sum = MAX(hole_sum, 1.0e-6f);
        part_sum = MAX(part_sum, 1.0e-6f);

        float* hole_percentages = md_temp_alloc_array(temp, float, nto->group.count);
        float* part_percentages = md_temp_alloc_array(temp, float, nto->group.count);
        for (size_t i = 0; i < nto->group.count; ++i) {
            hole_percentages[i] = nto->transition_density_hole[i] / hole_sum;
            part_percentages[i] = nto->transition_density_part[i] / part_sum;
        }

        float* start_positions = md_temp_alloc_array(temp, float, nto->group.count);
        float cur_bottom_pos = plot_area.Min.x;
        for (size_t i = 0; i < nto->group.count; ++i) {
            start_positions[i] = cur_bottom_pos;
            cur_bottom_pos += bars_avail_width * hole_percentages[i];
            if (hole_percentages[i] != 0.0f) {
                cur_bottom_pos += gap_size;
            }
        }

        float* end_positions = md_temp_alloc_array(temp, float, nto->group.count);
        float* sub_end_positions = md_temp_alloc_array(temp, float, nto->group.count);
        float cur_pos = plot_area.Min.x;
        for (size_t end_i = 0; end_i < nto->group.count; ++end_i) {
            end_positions[end_i] = cur_pos;
            sub_end_positions[end_i] = cur_pos;
            cur_pos += bars_avail_width * part_percentages[end_i];
            if (hole_percentages[end_i] != 0.0f) {
                cur_pos += gap_size;
            }
        }

        ImVec2* curve_midpoints  = md_temp_alloc_array(temp, ImVec2, nto->group.count * nto->group.count);
        ImVec2* curve_directions = md_temp_alloc_array(temp, ImVec2, nto->group.count * nto->group.count);
        float* curve_percentages = md_temp_alloc_array(temp, float,  nto->group.count * nto->group.count);

        const float font_size = ImGui::GetFontSize();
        const md_vg_color_t text_color = MD_VG_COLOR_RGB(0, 0, 0);
        const md_vg_color_t outline_color = MD_VG_COLOR_RGBA(0, 0, 0, 128);

        for (size_t beg_i = 0; beg_i < nto->group.count; ++beg_i) {
            if (hole_percentages[beg_i] == 0.0f) continue;

            ImVec2 start_pos = { start_positions[beg_i], plot_area.Max.y - bar_height + 0.1f * bar_height };
            ImVec4 start_col = vec_cast(nto->group.color[beg_i]);
            start_col.w = 0.5f;

            for (size_t end_i = 0; end_i < nto->group.count; ++end_i) {
                const float percentage = nto->transition_matrix[end_i * nto->group.count + beg_i];
                if (percentage <= 0.0f) continue;

                const float width = bars_avail_width * percentage;
                ImVec2 end_pos = { sub_end_positions[end_i], plot_area.Min.y + bar_height - 0.1f * bar_height };
                ImVec4 end_col = vec_cast(nto->group.color[end_i]);

                ImVec2 p1 = {};
                ImVec2 p2 = {};
                ImVec2 p3 = {};
                ImVec2 p4 = {};
                compute_vertical_sankey_bezier(&p1, &p2, &p3, &p4, start_pos, end_pos, width);

                md_vg_style_t flow_style = md_vg_style_stroke(vg_color_from_imvec4(start_col), width);
                flow_style.line_cap = MD_VG_LINE_CAP_BUTT;
                flow_style.line_join = MD_VG_LINE_JOIN_ROUND;
                if (end_i != beg_i) {
                    flow_style.stroke = md_vg_paint_linear_gradient(vec2_from_imvec2({ start_pos.x, start_pos.y }), vec2_from_imvec2({ start_pos.x, end_pos.y }), vg_color_from_imvec4(start_col), vg_color_from_imvec4(end_col));
                }
                md_vg_add_bezier_cubic_styled(scene, vec2_from_imvec2(p1), vec2_from_imvec2(p2), vec2_from_imvec2(p3), vec2_from_imvec2(p4), &flow_style);

                char label[32];
                snprintf(label, sizeof(label), "%3.2f%%", percentage * 100.0f);
                const ImVec2 label_size = ImGui::CalcTextSize(label);
                if (width > label_size.x) {
                    const ImVec2 midpoint = (start_pos + end_pos) * 0.5f + ImVec2{ width * 0.5f, 0.0f };
                    const vec2_t direction = vec2_normalize({ start_pos.x - end_pos.x, start_pos.y - end_pos.y });
                    const size_t idx = beg_i * nto->group.count + end_i;
                    curve_midpoints[idx] = midpoint;
                    curve_directions[idx] = { direction.x, direction.y };
                    curve_percentages[idx] = percentage;
                }

                start_pos.x += width;
                sub_end_positions[end_i] += width;
            }
        }

        bool text_overlap = true;
        ImVec2 curve_text_size = ImGui::CalcTextSize("99.99%");
        const float text_spacing = 4.0f;

        while (text_overlap) {
            text_overlap = false;
            for (size_t beg_i = 0; beg_i < nto->group.count; ++beg_i) {
                for (size_t end_i = 0; end_i < nto->group.count; ++end_i) {
                    const size_t index1 = beg_i * nto->group.count + end_i;
                    ImRect rect1 = { curve_midpoints[index1] - curve_text_size * 0.5f, curve_midpoints[index1] + curve_text_size * 0.5f };

                    for (size_t beg_j = 0; beg_j < nto->group.count; ++beg_j) {
                        for (size_t end_j = 0; end_j < nto->group.count; ++end_j) {
                            const size_t index2 = beg_j * nto->group.count + end_j;
                            if (index1 == index2) continue;

                            ImRect rect2 = { curve_midpoints[index2] - curve_text_size * 0.5f, curve_midpoints[index2] + curve_text_size * 0.5f };
                            if (curve_midpoints[index1].y == 0.0f || curve_midpoints[index2].y == 0.0f) continue;

                            while (rect1.Overlaps(rect2)) {
                                text_overlap = true;
                                curve_midpoints[index1] += curve_directions[index1] * curve_text_size.y * text_spacing;
                                curve_midpoints[index2] -= curve_directions[index2] * curve_text_size.y * text_spacing;
                                rect1 = { curve_midpoints[index1] - curve_text_size * 0.5f, curve_midpoints[index1] + curve_text_size * 0.5f };
                                rect2 = { curve_midpoints[index2] - curve_text_size * 0.5f, curve_midpoints[index2] + curve_text_size * 0.5f };
                            }
                        }
                    }
                }
            }
        }

        for (size_t beg_i = 0; beg_i < nto->group.count; ++beg_i) {
            for (size_t end_i = 0; end_i < nto->group.count; ++end_i) {
                const size_t index = beg_i * nto->group.count + end_i;
                const float percentage = curve_percentages[index];
                if (percentage <= 0.0f) continue;

                char label[16];
                snprintf(label, sizeof(label), "%3.2f%%", percentage * 100.0f);
                const ImVec2 pos = aligned_text_pos(label, curve_midpoints[index], { 0.5f, 0.5f });
                md_vg_add_text(scene, vec2_from_imvec2(pos), font_size, text_color, str_from_cstr(label));
            }
        }

        bool* show_start_text = md_temp_alloc_array(temp, bool, nto->group.count);
        bool* show_end_text = md_temp_alloc_array(temp, bool, nto->group.count);
        int* size_order = md_temp_alloc_array(temp, int, nto->group.count);
        float* text_sizes = md_temp_alloc_array(temp, float, nto->group.count);
        ImRect* rectangles = md_temp_alloc_array(temp, ImRect, nto->group.count);

        for (size_t i = 0; i < nto->group.count; ++i) {
            show_start_text[i] = true;
            show_end_text[i] = true;
        }

        if (hide_text_overlaps) {
            for (int i = 0; i < (int)nto->group.count; ++i) {
                size_order[i] = i;
                text_sizes[i] = MAX(ImGui::CalcTextSize(nto->group.label[i]).x, ImGui::CalcTextSize("99.99%").x);
                const float midpoint = start_positions[i] + bars_avail_width * hole_percentages[i] * 0.5f;
                rectangles[i] = { midpoint - text_sizes[i] * 0.5f, 0.0f, midpoint + text_sizes[i] * 0.5f, 1.0f };
            }

            std::sort(size_order, size_order + nto->group.count, [values = hole_percentages](int a, int b) {
                return values[a] > values[b];
            });

            for (int i = 1; i < (int)nto->group.count; ++i) {
                for (int j = i - 1; j >= 0; --j) {
                    if (hole_percentages[i] == 0.0f) {
                        show_start_text[size_order[i]] = false;
                        break;
                    }
                    if (show_start_text[size_order[j]] && rectangles[size_order[i]].Overlaps(rectangles[size_order[j]])) {
                        show_start_text[size_order[i]] = false;
                        break;
                    }
                }
            }

            for (int i = 0; i < (int)nto->group.count; ++i) {
                size_order[i] = i;
                const float midpoint = end_positions[i] + bars_avail_width * part_percentages[i] * 0.5f;
                rectangles[i] = { midpoint - text_sizes[i] * 0.5f, 0.0f, midpoint + text_sizes[i] * 0.5f, 1.0f };
            }

            std::sort(size_order, size_order + nto->group.count, [values = part_percentages](int a, int b) {
                return values[a] > values[b];
            });

            for (int i = 1; i < (int)nto->group.count; ++i) {
                for (int j = i - 1; j >= 0; --j) {
                    if (hole_percentages[i] == 0.0f) {
                        show_end_text[size_order[i]] = false;
                        break;
                    }
                    if (show_end_text[size_order[j]] && rectangles[size_order[i]].Overlaps(rectangles[size_order[j]])) {
                        show_end_text[size_order[i]] = false;
                        break;
                    }
                }
            }
        }

        for (size_t i = 0; i < nto->group.count; ++i) {
            if (hole_percentages[i] == 0.0f) continue;

            const ImVec4 bar_color = vec_cast(nto->group.color[i]);

            char start_label[16];
            char end_label[16];
            snprintf(start_label, sizeof(start_label), "%3.2f%%", hole_percentages[i] * 100.0f);
            snprintf(end_label, sizeof(end_label), "%3.2f%%", part_percentages[i] * 100.0f);

            ImVec2 start_p0 = { start_positions[i], plot_area.Max.y - bar_height };
            ImVec2 start_p1 = { start_positions[i] + bars_avail_width * hole_percentages[i], plot_area.Max.y };
            ImVec2 start_midpoint = { (start_p0.x + start_p1.x) * 0.5f, start_p1.y };

            ImVec2 end_p0 = { end_positions[i], plot_area.Min.y };
            ImVec2 end_p1 = { end_positions[i] + bars_avail_width * part_percentages[i], plot_area.Min.y + bar_height };
            ImVec2 end_midpoint = { (end_p0.x + end_p1.x) * 0.5f, end_p0.y };

            md_vg_add_rect_filled(scene, vec2_from_imvec2(start_p0), vec2_from_imvec2(start_p1), vg_color_from_imvec4(bar_color), 0.0f);
            md_vg_add_rect(scene, vec2_from_imvec2(start_p0), vec2_from_imvec2(start_p1), outline_color, 0.0f, 1.0f);
            if (show_start_text[i]) {
                md_vg_add_text(scene, vec2_from_imvec2(aligned_text_pos(nto->group.label[i], start_midpoint, { 0.5f, -0.2f })), font_size, text_color, str_from_cstr(nto->group.label[i]));
                md_vg_add_text(scene, vec2_from_imvec2(aligned_text_pos(start_label, start_midpoint, { 0.5f, -1.2f })), font_size, text_color, str_from_cstr(start_label));
            }

            md_vg_add_rect_filled(scene, vec2_from_imvec2(end_p0), vec2_from_imvec2(end_p1), vg_color_from_imvec4(bar_color), 0.0f);
            md_vg_add_rect(scene, vec2_from_imvec2(end_p0), vec2_from_imvec2(end_p1), outline_color, 0.0f, 1.0f);
            if (show_end_text[i]) {
                md_vg_add_text(scene, vec2_from_imvec2(aligned_text_pos(nto->group.label[i], end_midpoint, { 0.5f, 1.2f })), font_size, text_color, str_from_cstr(nto->group.label[i]));
                md_vg_add_text(scene, vec2_from_imvec2(aligned_text_pos(end_label, end_midpoint, { 0.5f, 2.2f })), font_size, text_color, str_from_cstr(end_label));
            }
        }
    }

    static inline bool export_transition_diagram(ImVec2 size, Nto* nto, bool hide_text_overlaps) {
        if (!nto || !nto->transition_density_hole || !nto->transition_density_part || !nto->transition_matrix) {
            MD_LOG_ERROR("Transition diagram export failed: transition diagram data is not initialized");
            return false;
        }

        char path_buf[1024] = {};
        if (!application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, STR_LIT("svg"))) {
            return false;
        }

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };
        md_vg_scene_t scene = {};
        md_vg_scene_init(&scene, { size.x, size.y }, md_temp_allocator(temp));
        vg_sankey_diagram(&scene, { 0.0f, 0.0f, size.x, size.y }, nto, hide_text_overlaps);

        const bool success = md_vg_scene_write_svg_file(&scene, str_from_cstr(path_buf));
        if (success) {
            MD_LOG_INFO("Exported transition diagram to '%s'", path_buf);
        } else {
            MD_LOG_ERROR("Failed to export transition diagram to '%s'", path_buf);
        }
        return success;
    }

    static inline void im_sankey_diagram(ApplicationState* state, ImRect area, Nto* nto, bool hide_text_overlaps) {
        if (!nto->transition_density_hole) return;
        if (!nto->transition_density_part) return;

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        /*
        * A sankey diagram needs to implement bezier curves, that are connecting two points
        */
        ImDrawList* draw_list = ImGui::GetWindowDrawList();

        //Draw background
        //draw_list->AddRectFilled(area.Min, area.Max, ImGui::ColorConvertFloat4ToU32({ 1,1,1,1 }));
        //draw_list->AddRect(area.Min, area.Max, ImGui::ColorConvertFloat4ToU32({ 0,0,0,1 }));

        ImRect plot_area = area;
        const float plot_percent = 0.8f;
        plot_area.Expand({ area.GetWidth() * -(1 - plot_percent), area.GetHeight() * -(1 - plot_percent) });
        //draw_list->AddRect(plot_area.Min, plot_area.Max, ImGui::ColorConvertFloat4ToU32({ 1,0,0,1 })); //Use this to draw debug of plot area

        ////The given data
        //const char* names[] = { "THIO", "QUIN" };
        //float initial_percentages[2] = { 0.3, 0.7 };
        //float transitions[2][2] = {
        //    /*
        //    *   A B   col -->
        //    * A
        //    * B
        //    * |
        //    * |
        //    * v
        //    * row
        //    */
        //    {0.22, 0.78}, //Read as "Starting from A, 22% of A_Start goes to A_End, and 78% goes to B_End
        //    {0.0, 1.0}
        //};

        //Bar definitions
        const float bar_height = plot_area.GetHeight() * 0.05f;
        int num_bars = 0;
        for (size_t i = 0; i < nto->group.count; i++) {
            if (nto->transition_density_hole[i] != 0.0f) {
                num_bars++;
            }
        }
        int num_gaps = num_bars - 1;
        float gap_size = plot_area.GetWidth() * 0.05f;
        float bars_avail_width = plot_area.GetWidth() - gap_size * num_gaps;

        //Calculate bar percentages
        float hole_sum = 0;
        float part_sum = 0;
        for (size_t i = 0; i < nto->group.count; i++) {
            hole_sum += nto->transition_density_hole[i];
            part_sum += nto->transition_density_part[i];
        }

        // Prevent division by zero
        hole_sum = MAX(hole_sum, 1.0e-6f);
        part_sum = MAX(part_sum, 1.0e-6f);

        float* hole_percentages = md_temp_alloc_array(temp, float, nto->group.count);
        float* part_percentages = md_temp_alloc_array(temp, float, nto->group.count);
        for (size_t i = 0; i < nto->group.count; i++) {
            hole_percentages[i] = nto->transition_density_hole[i] / hole_sum;
            part_percentages[i] = nto->transition_density_part[i] / part_sum;
        }

        //Calculate start positions
        float* start_positions = md_temp_alloc_array(temp, float, nto->group.count);
        float cur_bottom_pos = plot_area.Min.x;
        for (size_t i = 0; i < nto->group.count; i++) {
            start_positions[i] = cur_bottom_pos;
            cur_bottom_pos += bars_avail_width * hole_percentages[i];
            if (hole_percentages[i] != 0.0) {
                cur_bottom_pos += gap_size;
            }
        }

        //Calculate end positions
        float* end_positions = md_temp_alloc_array(temp, float, nto->group.count);
        float* sub_end_positions = md_temp_alloc_array(temp, float, nto->group.count);
        float cur_pos = plot_area.Min.x;
        for (size_t end_i = 0; end_i < nto->group.count; end_i++) {
            end_positions[end_i] = cur_pos;
            sub_end_positions[end_i] = cur_pos;
            cur_pos += bars_avail_width * part_percentages[end_i];
            if (hole_percentages[end_i] != 0.0) {
                cur_pos += gap_size;
            }
        }

        float  hovered_flow_width = 0.0f;
        ImVec2 hovered_flow_beg = {};
        ImVec2 hovered_flow_end = {};
        ImVec2 mouse_pos = ImGui::GetMousePos();

        char mouse_label[32] = { 0 };

        const float min_test_width = 1.0f;

        //Draw curves
        float* curve_widths      = md_temp_alloc_array(temp, float, nto->group.count * nto->group.count);
        ImVec2* curve_midpoints  = md_temp_alloc_array(temp, ImVec2, nto->group.count * nto->group.count);
        ImVec2* curve_directions = md_temp_alloc_array(temp, ImVec2, nto->group.count * nto->group.count);
        float* curve_percentages = md_temp_alloc_array(temp, float, nto->group.count * nto->group.count);

        for (size_t beg_i = 0; beg_i < nto->group.count; beg_i++) {
            if (hole_percentages[beg_i] != 0.0) {

                ImVec2 start_pos = { start_positions[beg_i], plot_area.Max.y - bar_height + 0.1f * bar_height };
                ImVec4 start_col = vec_cast(nto->group.color[beg_i]);

                start_col.w = 0.5f;
                for (size_t end_i = 0; end_i < nto->group.count; end_i++) {
                    float percentage = nto->transition_matrix[end_i * nto->group.count + beg_i];
                    ImVec4 end_col = vec_cast(nto->group.color[end_i]);

                    if (percentage > 0.0f) {
                        float width = bars_avail_width * percentage;
                        ImVec2 end_pos = { sub_end_positions[end_i], plot_area.Min.y + bar_height - 0.1f * bar_height };

                        int vert_beg = draw_list->VtxBuffer.Size;
                        ImVec2 mouse_delta = draw_vertical_sankey_flow(draw_list, start_pos, end_pos, width, ImGui::ColorConvertFloat4ToU32(start_col));
                        int vert_end = draw_list->VtxBuffer.Size;

                        // The width test here is clamped to be atleast some pixel in width, otherwise it would be impossible to hit thin flows
                        bool flow_hovered = mouse_delta.x < MAX(width * 0.5f, min_test_width) && mouse_delta.y < 1.0e-4;

                        // Apply a gradient based on y-value from start color to end color if they belong to different groups
                        if (end_i != beg_i) {
                            ImVec2 grad_p0 = { start_pos.x, start_pos.y };
                            ImVec2 grad_p1 = { start_pos.x, end_pos.y };
                            ImGui::ShadeVertsLinearColorGradientKeepAlpha(draw_list, vert_beg, vert_end, grad_p0, grad_p1, ImGui::ColorConvertFloat4ToU32(start_col), ImGui::ColorConvertFloat4ToU32(end_col));
                        }

                        char label[32];
                        snprintf(label, sizeof(label), "%3.2f%%", percentage * 100);
                        const ImVec2 label_size = ImGui::CalcTextSize(label);

                        // @NOTE: The comparison for width here is to ensure that the narrowest flow is always the one hovered (when multiple flows are overlapping)
                        if (flow_hovered && (hovered_flow_width == 0.0f || width < hovered_flow_width)) {
                            hovered_flow_width = width;
                            hovered_flow_beg = start_pos;
                            hovered_flow_end = end_pos;
                            MEMCPY(mouse_label, label, sizeof(label));
                        }

                        if (width > label_size.x) {
                            ImVec2 midpoint = (start_pos + end_pos) * 0.5 + ImVec2{ width / 2, 0 };
                            ImVec2 direction = vec_cast(vec2_normalize(vec_cast(start_pos - end_pos)));
                            curve_midpoints[beg_i * nto->group.count + end_i] = midpoint;
                            curve_directions[beg_i * nto->group.count + end_i] = direction;
                            curve_widths[beg_i * nto->group.count + end_i] = width;
                            curve_percentages[beg_i * nto->group.count + end_i] = percentage;
                            //draw_aligned_text(draw_list, label, midpoint, { 0.5, 0.5 });
                        }

                        start_pos.x += width;
                        sub_end_positions[end_i] += width;
                    }
                }
            }
        }

        //Draw Curve Text
        //For every midpoint, we check if another midpoint is to close to that midpoint
        bool text_overlap = true;
        ImVec2 text_size = ImGui::CalcTextSize("99.99%");
        float text_spacing = 4;

        while (text_overlap) {
            text_overlap = false;
            for (size_t beg_i = 0; beg_i < nto->group.count; beg_i++) {
                for (size_t end_i = 0; end_i < nto->group.count; end_i++) {
                    size_t index1 = beg_i * nto->group.count + end_i;
                    ImVec2 midpoint1 = curve_midpoints[index1];
                    ImVec2 direction1 = curve_directions[index1];
                    ImRect rect1 = { midpoint1 - (text_size / 2), midpoint1 + (text_size / 2) };

                    for (size_t beg_j = 0; beg_j < nto->group.count; beg_j++) {
                        for (size_t end_j = 0; end_j < nto->group.count; end_j++) {
                            size_t index2 = (beg_j * nto->group.count + end_j);
                            if (index1 != index2) {
                                ImVec2 midpoint2 = curve_midpoints[index2];
                                ImVec2 direction2 = curve_directions[index2] * -1.0;
                                ImRect rect2 = { midpoint2 - (text_size / 2), midpoint2 + (text_size / 2) };
                    
                                if ((midpoint1.y != 0 && midpoint2.y != 0)) {
                                    while (rect1.Overlaps(rect2)){
                                        text_overlap = true;
                                        curve_midpoints[index1] += direction1 * text_size.y * text_spacing;
                                        curve_midpoints[index2] += direction2 * text_size.y * text_spacing;
                                        rect1 = { curve_midpoints[index1] - (text_size / 2), curve_midpoints[index1] + (text_size / 2) };
                                        rect2 = { curve_midpoints[index2] - (text_size / 2), curve_midpoints[index2] + (text_size / 2) };
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (size_t beg_i = 0; beg_i < nto->group.count; beg_i++) {
            for (size_t end_i = 0; end_i < nto->group.count; end_i++) {
                size_t index = beg_i * nto->group.count + end_i;
                ImVec2 midpoint = curve_midpoints[index];
                float percentage = curve_percentages[index];
                if (percentage > 0) {
                    char label[16];
                    snprintf(label, sizeof(label), "%3.2f%%", curve_percentages[beg_i * nto->group.count + end_i] * 100);
                    draw_aligned_text(draw_list, label, midpoint, {0.5, 0.5});
                }
            }
        }

        // Some flow is hovered, So we render a new flow on top of the existing ones to ensure that this is rendered on top
        if (hovered_flow_width > 0.0f) {
            draw_vertical_sankey_flow(draw_list, hovered_flow_beg, hovered_flow_end, hovered_flow_width, ImColor(1.0f, 1.0f, 1.0f, 0.2f));
        }

        //Calculate text visibility
        bool* show_start_text  = md_temp_alloc_array(temp, bool,   nto->group.count);
        bool* show_end_text    = md_temp_alloc_array(temp, bool,   nto->group.count);
        int*  size_order       = md_temp_alloc_array(temp, int,    nto->group.count);
        float* text_sizes      = md_temp_alloc_array(temp, float,  nto->group.count);
        ImRect* rectangles     = md_temp_alloc_array(temp, ImRect, nto->group.count);

        for (size_t i = 0; i < nto->group.count; i++) {
            show_start_text[i] = true;
            show_end_text[i] = true;
        }
        
        if (hide_text_overlaps) {
            //Start texts
            for (int i = 0; i < (int)nto->group.count; i++) {
                size_order[i] = i;
                text_sizes[i] = MAX(ImGui::CalcTextSize(nto->group.label[i]).x, ImGui::CalcTextSize("99.99%").x);

                float midpoint = start_positions[i] + bars_avail_width * hole_percentages[i] * 0.5f;
                rectangles[i] = { midpoint - text_sizes[i] * 0.5f, 0, midpoint + text_sizes[i] * 0.5f, 1};
            }
        
            std::sort(size_order, size_order + nto->group.count, [values = hole_percentages](int a, int b) {
                return values[a] > values[b];  // Sort in descending order
            });

            for (int i = 1; i < (int)nto->group.count; i++) { //First one is always drawn
                for (int j = i - 1; j >= 0 ; j--) {//The bigger ones
                    if (hole_percentages[i] == 0.0f) {
                        show_start_text[size_order[i]] = false;
                        break;
                    }
                    else if (show_start_text[size_order[j]] && rectangles[size_order[i]].Overlaps(rectangles[size_order[j]])) {
                        show_start_text[size_order[i]] = false;
                        break;
                    }
                }
            }

            //End texts
            for (int i = 0; i < (int)nto->group.count; i++) {
                size_order[i] = i;

                float midpoint = end_positions[i] + bars_avail_width * part_percentages[i] * 0.5f;
                rectangles[i] = { midpoint - text_sizes[i] * 0.5f, 0, midpoint + text_sizes[i] * 0.5f, 1 };
            }

            std::sort(size_order, size_order + nto->group.count, [values = part_percentages](int a, int b) {
                return values[a] > values[b];  // Sort in descending order
            });

            for (int i = 1; i < (int)nto->group.count; i++) { //First one is always drawn
                for (int j = i - 1; j >= 0; j--) {//The bigger ones
                    if (hole_percentages[i] == 0.0f) {
                        show_end_text[size_order[i]] = false;
                        break;
                    }
                    else if (show_end_text[size_order[j]] && rectangles[size_order[i]].Overlaps(rectangles[size_order[j]])) {
                        show_end_text[size_order[i]] = false;
                        break;
                    }
                }
            }
        }

        //Draw bars
        for (size_t i = 0; i < nto->group.count; i++) {
            if (hole_percentages[i] != 0.0f) {
                ImVec4 bar_color = vec_cast(nto->group.color[i]);

                char start_label1[16];
                char end_label1[16]; 
                snprintf(start_label1, sizeof(start_label1), "%3.2f%%", hole_percentages[i] * 100);
                snprintf(end_label1, sizeof(end_label1), "%3.2f%%", part_percentages[i] * 100);

                char start_label2[16];
                char end_label2[16];
                snprintf(start_label2, sizeof(start_label2), "%3.2f%%", hole_percentages[i + 1] * 100);
                snprintf(end_label2, sizeof(end_label2), "%3.2f%%", part_percentages[i + 1] * 100);

                //Calculate start
                ImVec2 start_p0 = { start_positions[i], plot_area.Max.y - bar_height };
                ImVec2 start_p1 = { start_positions[i] + bars_avail_width * hole_percentages[i], plot_area.Max.y };
                ImVec2 start_midpoint = { (start_p0.x + start_p1.x) * 0.5f, start_p1.y };
                ImRect start_bar = ImRect{ start_p0, start_p1 };

                //Calculate end
                ImVec2 end_p0 = { end_positions[i], plot_area.Min.y };
                ImVec2 end_p1 = { end_positions[i] + bars_avail_width * part_percentages[i], plot_area.Min.y + bar_height };
                ImRect end_bar = ImRect{ end_p0, end_p1 };
                ImVec2 end_midpoint = { (end_p0.x + end_p1.x) * 0.5f, end_p0.y };

                bool index_hovered = false;
                if (start_bar.Contains(mouse_pos)) {
                    MEMCPY(mouse_label, start_label1, sizeof(start_label1));
                    index_hovered = true;
                }
                else if (end_bar.Contains(mouse_pos)) {
                    MEMCPY(mouse_label, end_label1, sizeof(end_label1));
                    index_hovered = true;
                }

                if (index_hovered) {
                    nto->group.hovered_index = (int8_t)i;
                    bar_color = make_highlight_color(bar_color);
                    for (size_t atom_i = 0; atom_i < md_array_size(nto->atom_group_idx); atom_i++) {
                        if (nto->atom_group_idx[atom_i] == (uint32_t)i) {
                            md_bitfield_set_bit(&state->selection.highlight_mask, atom_i);
                        }
                    }
                }

                //Draw start
                draw_list->AddRectFilled(start_p0, start_p1, ImGui::ColorConvertFloat4ToU32(bar_color));
                draw_list->AddRect(start_p0, start_p1, ImGui::ColorConvertFloat4ToU32({ 0,0,0,0.5 }));
                if (show_start_text[i]) {
                    draw_aligned_text(draw_list, nto->group.label[i], start_midpoint, { 0.5, -0.2 });
                    draw_aligned_text(draw_list, start_label1, start_midpoint, { 0.5, -1.2 });
                }

                //Draw end
                draw_list->AddRectFilled(end_p0, end_p1, ImGui::ColorConvertFloat4ToU32(bar_color));
                draw_list->AddRect(end_p0, end_p1, ImGui::ColorConvertFloat4ToU32({ 0,0,0,0.5 }));
                if (show_end_text[i]) {
                    draw_aligned_text(draw_list, nto->group.label[i], end_midpoint, { 0.5, 1.2 });
                    draw_aligned_text(draw_list, end_label1, end_midpoint, { 0.5, 2.2 });
                }
            }
        }

        size_t mouse_label_len = strnlen(mouse_label, sizeof(mouse_label));
        if (mouse_label_len > 0) {
            ImVec2 offset = { 15, 15 };
            ImVec2 pos = ImGui::GetMousePos() + offset;
            draw_list->AddText(pos, IM_COL32_BLACK, mouse_label, mouse_label + mouse_label_len);
        }
    }

    static inline double convert_value_from_au(double value, x_unit_t unit) {
        switch (unit) {
        case X_UNIT_EV:
            return value * HARTREE_TO_EV;
        case X_UNIT_NM:
            return 45.56331418628337 / value;
        case X_UNIT_CM_INVERSE:
            return value * 219479.86946633097;
        case X_UNIT_HARTREE:
            return value;
        default:
            ASSERT(false); // Should not happen
            return value;
        }
    }

    static inline double convert_value_to_au(double value, x_unit_t unit) {
        switch (unit) {
        case X_UNIT_EV:
            return value * EV_TO_HARTREE;
        case X_UNIT_NM:
            return 45.56331418628337 / value;
        case X_UNIT_CM_INVERSE:
            return value / 219479.86946633097;
        case X_UNIT_HARTREE:
            return value;
        default:
            ASSERT(false); // Should not happen
            return value;
        }
    }

    static inline void convert_values_from_au(double* out_values, const double* in_values, size_t num_values, x_unit_t unit) {
        if (!out_values || !in_values || num_values == 0) return;

        switch (unit) {
        case X_UNIT_EV:
            for (size_t i = 0; i < num_values; ++i) {
                out_values[i] = in_values[i] * HARTREE_TO_EV;
            }
            break;
        case X_UNIT_NM:
            for (size_t i = 0; i < num_values; ++i) {
                out_values[i] = 45.56331418628337 / in_values[i];
            }
            break;
        case X_UNIT_CM_INVERSE:
            for (size_t i = 0; i < num_values; ++i) {
                out_values[i] = in_values[i] * 219479.86946633097;
            }
            break;
        case X_UNIT_HARTREE:
            for (size_t i = 0; i < num_values; ++i) {
                out_values[i] = in_values[i];
            }
            break;
        default:
            ASSERT(false); // Should not happen
            break;
        }
    }

    static inline double sigma_from_fwhm(double fwhm) {
        return fwhm * 0.21233045007200476068;
    }

    // Evaluates lorentzian broadening around a point x, given a set of peaks and their strengths.
    // This is a NON normalized version.
    // To normalize this, just multiply with the normalization factor
	static inline double lorentzian_base(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
		const double gamma = fwhm * 0.5;
		const double gamma2 = gamma * gamma;

		double y = 0;
		for (size_t k = 0; k < num_peaks; k++) {
			double dx = x - peaks_x[k];
			y += peaks_y[k] * gamma / fma(dx, dx, gamma2);
		}
		return y;
	}

    static inline double lorentzian_abs(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        const double gamma = fwhm * 0.5;
		const double gamma2 = gamma * gamma;
        double y = 0.0;
		for (size_t k = 0; k < num_peaks; k++) {
			double dx = x - peaks_x[k];
			y += peaks_y[k] / peaks_x[k] * gamma / fma(dx, dx, gamma2);
		}
        return y / PI;
    }

    static inline double lorentzian_ecd(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        return lorentzian_base(x, peaks_x, peaks_y, num_peaks, fwhm);
    }

    static inline double lorentzian_vib(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        return (1.0 / PI) * lorentzian_base(x, peaks_x, peaks_y, num_peaks, fwhm);
    }

    static inline double lorentzian_tpa(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        return x * x * lorentzian_base(x, peaks_x, peaks_y, num_peaks, fwhm);
    }

    static inline double lorentzian_xps(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        return (1.0 / PI) * lorentzian_base(x, peaks_x, peaks_y, num_peaks, fwhm);
    }

    // Evaluates gaussian broadening around a point x, given a set of peaks and their strengths.
    // This is a NON normalized version.
    // To normalize this, just multiply with the normalization factor
	static inline double gaussian_base(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
		const double sigma = sigma_from_fwhm(fwhm);
		const double factor = -0.5 / (sigma * sigma);

        double y = 0.0;
        for (size_t k = 0; k < num_peaks; k++) {
            double dx = x - peaks_x[k];
            y += peaks_y[k] * exp(factor * dx * dx);
        }
		return y;
	}

    static inline double gaussian_abs(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        const double sigma = sigma_from_fwhm(fwhm);
		const double factor = -0.5 / (sigma * sigma);
        double y = 0.0;
        for (size_t k = 0; k < num_peaks; k++) {
            double dx = x - peaks_x[k];
            y += peaks_y[k] / peaks_x[k] * exp(factor * dx * dx);
        }
		return y / (sigma * sqrt(2.0 * PI));
    }

    static inline double gaussian_ecd(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        const double sigma = sigma_from_fwhm(fwhm);
		const double y = gaussian_base(x, peaks_x, peaks_y, num_peaks, fwhm);
		return y * PI / (sigma * sqrt(2.0 * PI));
    }

    static inline double gaussian_vib(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        const double sigma = sigma_from_fwhm(fwhm);
		const double y = gaussian_base(x, peaks_x, peaks_y, num_peaks, fwhm);
		return y * sqrt(2.0) / (sigma * sqrt(PI));
    }

    static inline double gaussian_xps(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        const double sigma = sigma_from_fwhm(fwhm);
		const double y = gaussian_base(x, peaks_x, peaks_y, num_peaks, fwhm);
		return y / (sigma * sqrt(2.0 * PI));
    }


    //Constructs plot limits from peaks
    static inline ImPlotRect get_plot_limits(const double* x_samples, const double* y_peaks, size_t num_peaks, size_t num_samples, double ext_fac = 0.1) {
        ImPlotRect lim = { MIN(x_samples[0],x_samples[num_samples - 1]), MAX(x_samples[0],x_samples[num_samples - 1]),0,0};
        for (size_t i = 0; i < num_peaks; i++) {
            //Use Contains to check if values are within the limits, or if they should extend the limits
            if (lim.Y.Max < y_peaks[i]) {
                lim.Y.Max = y_peaks[i];
            }
            else if (lim.Y.Min > y_peaks[i]) {
                lim.Y.Min = y_peaks[i];
            }
        }

        double height = lim.Y.Max - lim.Y.Min;
        double width = lim.X.Max - lim.X.Min;

        lim.Y.Max += height * ext_fac;
        lim.Y.Min -= height * ext_fac;
        lim.X.Max += width * ext_fac;
        lim.X.Min -= width * ext_fac;

        //The y limits needs to be symmetric so that the spectra is not clipped. For example, if the y peak max is a positive value but the y spectra value is negative, we get issues otherwise. This ensures space for all the data.
        //This is needed because Y2 scaling is based on the maximum peak value
        double abs_max_y = MAX(fabs(lim.Y.Min), fabs(lim.Y.Max));
        lim.Y.Min = -abs_max_y;
        lim.Y.Max = abs_max_y;
        return lim;
    }

    //converts x and y peaks into pixel points in context of the current plot. Use between BeginPlot() and EndPlot()
    static inline void peaks_to_pixels(ImVec2* pixel_peaks, const double* x_peaks, const double* y_peaks, size_t num_peaks) {
        for (size_t i = 0; i < num_peaks; i++) {
            pixel_peaks[i] = ImPlot::PlotToPixels(ImPlotPoint{ x_peaks[i], y_peaks[i] });
        }
    }

    static inline int get_hovered_plot_peak(const double* xs, const double* ys, size_t count, double distance_threshold = 10.0) {
        ImPlotPoint mouse_pos = ImPlot::GetPlotMousePos();
        int closest_idx = -1;
        double closest_dist = DBL_MAX;
        for (size_t i = 0; i < count; ++i) {
            double dx = mouse_pos.x - xs[i];
            double dy = CLAMP(mouse_pos.y, 0.0, ys[i]) - mouse_pos.y;
            double dist = sqrt(dx * dx + dy * dy);
            if (dist < distance_threshold && dist < closest_dist) {
                closest_idx = (int)i;
                closest_dist = dist;
            }
        }
        return closest_idx;
    }

    static inline int get_hovered_plot_peak(const double* values, size_t count, double x_base = 0.0, double x_step = 1.0, double distance_threshold = 10.0) {
        ImPlotPoint mouse_pos = ImGui::GetMousePos();
        int closest_idx = -1;
        double closest_dist = DBL_MAX;
        ImVec2 zero = ImPlot::PlotToPixels(0, 0);
        for (size_t i = 0; i < count; ++i) {
            ImVec2 point = ImPlot::PlotToPixels(x_base + i * x_step, values[i]);
            double dx = mouse_pos.x - point.x;
            double dy = CLAMP(mouse_pos.y, MIN(zero.y, point.y), MAX(zero.y, point.y)) - mouse_pos.y;
            double dist = sqrt(dx * dx + dy * dy);
            if (dist < distance_threshold && dist < closest_dist) {
                closest_idx = (int)i;
                closest_dist = dist;
            }
        }
        return closest_idx;
    }

    // Returns peak index closest to mouse pixel position, assumes that x-values are sorted.
    static inline int get_hovered_peak(const ImVec2 mouse_pos, const ImVec2* pixel_peaks, const ImVec2* pixel_points, size_t num_peaks, bool y_flipped = false, double proxy_distance = 10.0) {
        int closest_idx = 0;
        double x = mouse_pos.x;
        double y = mouse_pos.y;
        double y_max = 0;
        double y_min = 0;
        double distance_x = 0;
        double distance_y = 0;
        double closest_distance = 0;
        double pixel_y0 = ImPlot::PlotToPixels(0, 0).y;

        //Keep in mind that pixel y is 0 at the top, so you flip the comparison compared to plot y. The code below still seems to work as intended though.

        for (size_t i = 0; i < num_peaks; i++) {
            y_max = MAX(pixel_peaks[i].y, pixel_y0);
            y_min = MIN(pixel_peaks[i].y, pixel_y0);

            //Check if the y location is within the range of y_min,ymax
            if (y > y_max) {
                distance_y = fabs(y - y_max);
                if (y_flipped) {
                    distance_y = fabs(y - pixel_points[i].y) < distance_y ? fabs(y - pixel_points[i].y) : distance_y;
                }
            }
            else if (y < y_min) {
                distance_y = fabs(y - y_min);
                if (!y_flipped) {
                    distance_y = fabs(y - pixel_points[i].y) < distance_y ? fabs(y - pixel_points[i].y) : distance_y;
                }
            }
            else {
                distance_y = 0;
            }

            distance_x = fabs(pixel_peaks[i].x - x);

            // We need a special case for i == 0 to set a reference for comparison, so that closest_distance does not start on 0;
            if (i == 0 && distance_y == 0 ){
                closest_distance = distance_x;
                closest_idx = 0;

            }
            else if (i == 0) {
                closest_distance = sqrt(pow(distance_x, 2) + pow(distance_y, 2)); // Is there a better function for doing this? Is this expensive?
                closest_idx = 0;
            }
            else if (distance_y == 0 && distance_x < closest_distance) {
                closest_distance = fabs(pixel_peaks[i].x - x);
                closest_idx = (int)i;
            }
            else if (sqrt(pow(distance_x, 2) + pow(distance_y, 2)) < closest_distance) {
                closest_distance = sqrt(pow(distance_x, 2) + pow(distance_y, 2));
                closest_idx = (int)i;
                //ImPlot::Annotation()
            }
            else if (distance_x > closest_distance){
                // We are now so far away that a closer bar will not occur, no matter the y value.
                break;
            }
        }
        return closest_distance < proxy_distance ? closest_idx : -1;
    }

    // Draws a bar within a plot at x = x, and from y = 0 to y = y, with specified width and color.
    // It also draws a point at the top with a radius
    static inline void draw_bar(double x, double y, double bar_width_in_pixels, double point_radius_in_pixels, ImVec4 color) {
        ImVec2 rect_min = ImPlot::PlotToPixels(x, 0, IMPLOT_AUTO, IMPLOT_AUTO);
        ImVec2 rect_max = ImPlot::PlotToPixels(x, y, IMPLOT_AUTO, IMPLOT_AUTO);

        rect_min = ImVec2{rect_min.x - (float)bar_width_in_pixels * 0.5f, rect_min.y};
        rect_max = ImVec2{rect_max.x + (float)bar_width_in_pixels * 0.5f, rect_max.y};
        ImVec2 pos = ImPlot::PlotToPixels(x, y, IMPLOT_AUTO, IMPLOT_AUTO);
        float  rad = (float)point_radius_in_pixels;

        ImPlot::PushPlotClipRect();
        ImDrawList& DrawList = *ImPlot::GetPlotDrawList();

        ImU32 col32 = ImGui::ColorConvertFloat4ToU32(color);
        DrawList.AddRectFilled(rect_min, rect_max, col32);
        DrawList.AddCircleFilled(pos, rad, col32);

        ImPlot::PopPlotClipRect();
    }

    // Compute a simple display multiplier to map oscillator-strength values (f) to eps for axis-syncing.
    // Strategy: Match peak magnitudes: multiplier = max(abs(eps_samples))/max(abs(f_peaks)).
    double compute_display_multiplier_from_peaks(const double* eps_samples, size_t eps_count, const double* peak_values, size_t peak_count) {
        double eps_max = 0.0;
        for (size_t i = 0; i < eps_count; ++i) {
            eps_max = MAX(eps_max, fabs(eps_samples[i]));
        }
        double peak_max = 0.0;
        for (size_t i = 0; i < peak_count; ++i) {
            peak_max = MAX(peak_max, fabs(peak_values[i]));
        }
        if (peak_max <= 0.0 || !isfinite(peak_max) || !isfinite(eps_max)) return 1.0;

        return eps_max / peak_max;
    }

    enum PlotPeaksFlags_ {
        PlotPeaksFlags_None = 0,
        PlotPeaksFlags_Bars = 1 << 0,
        PlotPeaksFlags_Points = 1 << 1,
        PlotPeaksFlags_All = PlotPeaksFlags_Bars | PlotPeaksFlags_Points
    };

    typedef uint32_t PlotPeaksFlags;

    // This is a helper function to plot peaks and handle hovering
    // Only valid to call within a plot context
    void plot_peaks(const char* title, const double x_values[], const double y_values[], size_t num_peaks, int& selected, int& hovered, PlotPeaksFlags flags = PlotPeaksFlags_All) {
        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        const double bar_width_in_pixels = 2.0;
        const double point_radius_in_pixels = 3.0;

        // If x values are not given we implicitly assume that they correspond to 1..num_peaks
        if (!x_values) {
            double* x_vals = md_temp_alloc_array(temp, double, num_peaks);
            for (size_t i = 0; i < num_peaks; i++) {
                x_vals[i] = (double)(i + 1);
            }
            x_values = x_vals;
        }

        // An item hidden through the legend draws nothing, so it must not claim a hover or take a
        // click either. ImPlot only resolves that inside BeginItem, which happens further down -
        // after both would already have run - so the item is looked up directly instead. It is
        // absent on the very first frame it is drawn, which BeginItem will treat as shown.
        const ImPlotItem* item = ImPlot::GetItem(title);
        const bool item_visible = (item == nullptr) || item->Show;

        // Dropped whenever the cursor is inside the plot, hidden or not, so that a hover taken
        // before the item was toggled off cannot keep driving a tooltip in the caller.
        if (ImPlot::IsPlotHovered()) {
            hovered = -1;
        }

        if (item_visible && ImPlot::IsPlotHovered()) {
            // Calculate the hovered index as the closest peak or 'bar' to the mouse position
            // This comparison needs to be done in pixel space

            const ImVec2 mouse_pos = ImGui::GetMousePos();
            const double pixel_tolerance = 7.0;
            const double min_line_d2  = pow(bar_width_in_pixels * 0.5 + pixel_tolerance, 2);  // 7 pixels tolerance
            const double min_point_d2 = pow(point_radius_in_pixels + pixel_tolerance, 2);     // 7 pixels tolerance

            double min_dist = DBL_MAX;
            for (size_t i = 0; i < num_peaks; i++) {
                ImVec2 min_raw = ImPlot::PlotToPixels(x_values[i], 0.0, IMPLOT_AUTO, IMPLOT_AUTO);
                ImVec2 max_raw = ImPlot::PlotToPixels(x_values[i], y_values[i], IMPLOT_AUTO, IMPLOT_AUTO);

                // The axis may be flipped, so we need to determine min and max correctly
                ImVec2 min = ImMin(min_raw, max_raw);
                ImVec2 max = ImMax(min_raw, max_raw);

                if (flags & PlotPeaksFlags_Bars) {
                    // Calculate distance to the line which forms the bar (1D clamp trick)
                    ImVec2 line_dx = mouse_pos - ImClamp(mouse_pos, min, max);
                    double line_d2 = line_dx.x * line_dx.x + line_dx.y * line_dx.y;

                    if (line_d2 < min_line_d2 && line_d2 < min_dist) {
                        hovered = (int)i;
                        min_dist = line_d2;
                    }
                }

                if (flags & PlotPeaksFlags_Points) {
                    // Calculate distance to the point at the top of the bar
                    ImVec2 point_dx = mouse_pos - max_raw;
                    double point_d2 = point_dx.x * point_dx.x + point_dx.y * point_dx.y;

                    if (point_d2 < min_point_d2 && point_d2 < min_dist) {
                        hovered = (int)i;
                        min_dist = point_d2;
                    }
                }
            }
        }

        // Update selected peak on click
        if (item_visible && ImGui::IsMouseReleased(ImGuiMouseButton_Left) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Left) &&
            ImPlot::IsPlotHovered()) {
            selected = hovered == selected ? -1 : hovered;
        }

        if (ImPlot::BeginItem(title)) {
            ImPlot::PushPlotClipRect();
            ImDrawList& DrawList = *ImPlot::GetPlotDrawList();
            const float rad = (float)point_radius_in_pixels;

            for (int i = 0; i < (int)num_peaks; i++) {
                const ImPlotNextItemData& s = ImPlot::GetItemData();
                ImVec4 color = (i == hovered) ? COLOR_PEAK_HOVER : (i == selected) ? COLOR_PEAK_SELECTED : s.Colors[ImPlotCol_Fill];

                const double x = x_values[i];
                const double y = y_values[i];

                ImVec2 rect_min = ImPlot::PlotToPixels(x, 0, IMPLOT_AUTO, IMPLOT_AUTO);
                ImVec2 rect_max = ImPlot::PlotToPixels(x, y, IMPLOT_AUTO, IMPLOT_AUTO);

                rect_min = ImVec2{rect_min.x - (float)bar_width_in_pixels * 0.5f, rect_min.y};
                rect_max = ImVec2{rect_max.x + (float)bar_width_in_pixels * 0.5f, rect_max.y};
                ImVec2 pos = ImPlot::PlotToPixels(x, y, IMPLOT_AUTO, IMPLOT_AUTO);

                ImU32 col32 = ImGui::ColorConvertFloat4ToU32(color);
                if (flags & PlotPeaksFlags_Bars) {
                    DrawList.AddRectFilled(rect_min, rect_max, col32);
                }

                if (flags & PlotPeaksFlags_Points) {
                    DrawList.AddCircleFilled(pos, rad, col32);
                }

                if (ImPlot::FitThisFrame()) {
                    ImPlot::FitPoint(ImPlotPoint(x, 0.0));
                    ImPlot::FitPoint(ImPlotPoint(x, y));
                }
            }
            ImPlot::PopPlotClipRect();

            ImPlot::EndItem();
        }
    }

    // =================================================================================================
    // RIXS map & spectrum
    //
    // ImGui/ImPlot port of VeloxChem's spectrumplot.plot_rixs_map(), plus a 1D spectrum widget that
    // shows a single row of that map (one incoming photon energy).
    // The implementation is split into a pure data pass (compute_rixs_map) and a drawing pass
    // (draw_rixs_map) so that the map can be rebuilt only when the settings or the input actually
    // change. Both take the input data through RixsMapInput and are therefore independent of the
    // md_vlx object; rixs_map_input_from_vlx() provides the glue.
    // =================================================================================================

    // Raw inputs, mirroring the 'rixs_results' dict consumed by plot_rixs_map().
    // Dimensions:
    //   P = num_photon_energies (incoming photons)
    //   F = num_final_states    (valence / final states)
    //   C = num_core_states     (core-excited states, only used for the XAS panel)
    // The 2D arrays are row-major [F][P], i.e. element (f, p) is at [f * P + p].
    struct RixsMapInput {
        const double* photon_energies_au     = nullptr;  // [P], a.u.
        const double* elastic_cross_sections = nullptr;  // [P], a.u.    (optional)
        const double* cross_sections         = nullptr;  // [F][P], a.u.
        const double* energy_losses_au       = nullptr;  // [F][P], a.u.
        const double* core_eigenvalues_au    = nullptr;  // [C], a.u.    (optional)
        const double* core_osc_strengths     = nullptr;  // [C]          (optional)
        double gamma_fwhm_ev = 0.0;                      // core-hole lifetime broadening, eV
        size_t num_photon_energies = 0;
        size_t num_final_states    = 0;
        size_t num_core_states     = 0;
    };

    // Collects the RIXS input from the currently loaded vlx object.
    // Returns false if the loaded data is not a RIXS calculation or is incomplete.
    bool rixs_map_input_from_vlx(RixsMapInput& out) {
        if (!vlx || md_vlx_rsp_type(vlx) != MD_VLX_RSP_RIXS) return false;

        out.photon_energies_au     = md_vlx_rsp_rixs_photon_energies(vlx);
        out.elastic_cross_sections = md_vlx_rsp_rixs_elastic_cross_sections(vlx);
        out.cross_sections         = md_vlx_rsp_rixs_cross_sections(vlx);
        out.energy_losses_au       = md_vlx_rsp_rixs_energy_losses(vlx);
        out.core_eigenvalues_au    = md_vlx_rsp_rixs_core_eigenvalues(vlx);
        out.core_osc_strengths     = md_vlx_rsp_rixs_core_osc_strengths(vlx);
        out.gamma_fwhm_ev          = md_vlx_rsp_rixs_gamma_fwhm_ev(vlx);
        out.num_photon_energies    = md_vlx_rsp_rixs_number_of_photon_energies(vlx);
        out.num_final_states       = md_vlx_rsp_rixs_number_of_final_states(vlx);
        out.num_core_states        = md_vlx_rsp_rixs_number_of_core_states(vlx);

        // The map itself needs the photon energies, the cross-sections and the energy losses.
        // The XAS panel is optional and is skipped if the core data is missing.
        if (out.num_photon_energies == 0 || out.num_final_states == 0) return false;
        if (!out.photon_energies_au || !out.cross_sections) return false;
        if (!out.energy_losses_au) return false;

        return true;
    }

    // Drops every cached array by rewinding the arena. The pointers are nulled rather than freed
    // individually: the arena owns nothing else, so a reset is sufficient and exhaustive.
    static void reset_rixs_cache(Rixs& r) {
        if (r.alloc) {
            md_vm_arena_reset(r.alloc);
        }
        r.map = r.grid = r.photon_ev = nullptr;
        r.xas_x = r.xas_y = r.core_ev = r.core_f = nullptr;
        r.num_rows = r.num_cols = r.num_xas = r.num_core = 0;
        r.valid = false;
    }

    // Evaluates the configured broadening kernel at x for a set of sticks.
    // Both kernels use exactly the same normalization as lorentzian_xps / gaussian_xps in
    // VeloxChem's spectrumplot.py, so the resulting intensities are directly comparable.
    static inline double rixs_broaden(broadening_mode_t mode, double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        return (mode == BROADENING_MODE_LORENTZIAN)
            ? lorentzian_xps(x, peaks_x, peaks_y, num_peaks, fwhm)
            : gaussian_xps  (x, peaks_x, peaks_y, num_peaks, fwhm);
    }

    // Rebuilds the cached RIXS map, XAS curve and axis extents from 'in' using the settings in 'r'.
    // This is the expensive part; it is O(P * N * F) and is therefore guarded by a hash in draw_rixs_map().
    void compute_rixs_map(Rixs& r, const RixsMapInput& in) {
        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        // Every array below is rebuilt from scratch, so nothing from the previous map needs to
        // survive. Rewinding here keeps the footprint at exactly one map regardless of how many
        // times the settings are changed. Also nulls the pointers, which matters for the early
        // returns further down: they leave the cache invalid but must not leave it dangling.
        reset_rixs_cache(r);

        md_allocator_i* alloc = r.alloc;
        ASSERT(alloc);

        const size_t P = in.num_photon_energies;
        const size_t F = in.num_final_states;

        const double* x_src = in.energy_losses_au;
        if (P == 0 || F == 0 || !in.photon_energies_au || !in.cross_sections || !x_src) return;

        // Include the elastic line as one extra stick per photon energy, if requested and available.
        const bool use_elastic = r.plot_elastic_line && in.elastic_cross_sections != nullptr;
        const size_t num_sticks = F + (use_elastic ? 1 : 0);

        // ---- Sort the incoming photon energies ascending (np.argsort) -------------------------------
        int* sort_idx = md_temp_alloc_array(temp, int, P);
        for (size_t i = 0; i < P; ++i) sort_idx[i] = (int)i;
        std::sort(sort_idx, sort_idx + P, [&](int a, int b) {
            return in.photon_energies_au[a] < in.photon_energies_au[b];
        });

        md_array_resize(r.photon_ev, P, alloc);
        for (size_t i = 0; i < P; ++i) {
            r.photon_ev[i] = in.photon_energies_au[sort_idx[i]] * HARTREE_TO_EV;
        }

        // ---- Gather the sticks per photon energy, column-contiguous ---------------------------------
        // stick_x[i * num_sticks + s], stick_y[i * num_sticks + s]
        double* stick_x = md_temp_alloc_array(temp, double, P * num_sticks);
        double* stick_y = md_temp_alloc_array(temp, double, P * num_sticks);

        double x_lo =  DBL_MAX;
        double x_hi = -DBL_MAX;

        for (size_t i = 0; i < P; ++i) {
            const size_t p = (size_t)sort_idx[i];
            double* sx = stick_x + i * num_sticks;
            double* sy = stick_y + i * num_sticks;

            for (size_t f = 0; f < F; ++f) {
                const double x = x_src[f * P + p] * HARTREE_TO_EV;
                sx[f] = x;
                sy[f] = in.cross_sections[f * P + p];
                x_lo = MIN(x_lo, x);
                x_hi = MAX(x_hi, x);
            }

            if (use_elastic) {
                // Elastic scattering: zero energy loss.
                const double x = 0.0;
                sx[F] = x;
                sy[F] = in.elastic_cross_sections[p];
                x_lo = MIN(x_lo, x);
                x_hi = MAX(x_hi, x);
            }
        }

        if (!(x_lo <= x_hi)) return;  // catches NaN / empty input

        // ---- Build the grid along the energy loss axis -------------------------------------------
        // The python version uses np.arange(val_min, val_max, xstep). Here the sample count is derived
        // from the same step but clamped, so that a small step cannot blow up the allocation.
        const double pad = 0.5;
        const double val_min = x_lo - pad;
        const double val_max = x_hi + pad;

        const double step = (r.x_step_ev > 0.0) ? r.x_step_ev : 0.01;
        size_t N = (size_t)((val_max - val_min) / step + 0.5);
        N = CLAMP(N, (size_t)2, (size_t)4096);

        md_array_resize(r.grid, N, alloc);
        for (size_t j = 0; j < N; ++j) {
            r.grid[j] = val_min + (val_max - val_min) * (double)j / (double)(N - 1);
        }

        // ---- Broaden each photon-energy column into a row of the map ---------------------------------
        md_array_resize(r.map, P * N, alloc);

        double map_lo =  DBL_MAX;
        double map_hi = -DBL_MAX;

        for (size_t i = 0; i < P; ++i) {
            const double* sx = stick_x + i * num_sticks;
            const double* sy = stick_y + i * num_sticks;

            // Reverse the row order so that the lowest photon energy ends up at the bottom of the
            // heatmap, matching matplotlib's origin='lower'.
            double* row = r.map + (P - 1 - i) * N;

            for (size_t j = 0; j < N; ++j) {
                const double y = rixs_broaden(r.broadening_mode, r.grid[j], sx, sy, num_sticks, r.broadening_fwhm_ev);
                row[j] = y;
                map_lo = MIN(map_lo, y);
                map_hi = MAX(map_hi, y);
            }
        }

        // ---- Normalize to the global maximum ---------------------------------------------------------
        if (r.normalize && map_hi > 0.0) {
            const double scl = 1.0 / map_hi;
            for (size_t k = 0; k < P * N; ++k) {
                r.map[k] *= scl;
            }
            map_lo *= scl;
            map_hi  = 1.0;
        }

        r.num_rows = P;
        r.num_cols = N;
        r.map_min  = map_lo;
        r.map_max  = map_hi;
        r.x_min_ev = val_min;
        r.x_max_ev = val_max;

        // ---- Photon energy axis extent (pixel edges, half a step beyond the outermost centers) -------
        double photon_step = 1.0;
        if (P > 1) {
            photon_step = (r.photon_ev[P - 1] - r.photon_ev[0]) / (double)(P - 1);

            // The heatmap assumes a uniform spacing; warn instead of asserting like the python version.
            for (size_t i = 1; i < P; ++i) {
                const double d = r.photon_ev[i] - r.photon_ev[i - 1];
                if (fabs(d - photon_step) > 1.0e-5 * MAX(1.0, fabs(photon_step))) {
                    if (!r.warned_non_uniform) {
                        r.warned_non_uniform = true;
                        MD_LOG_INFO("RIXS map: incoming photon energies are not uniformly spaced, the map will be distorted.");
                    }
                    break;
                }
            }
        }
        r.y_min_ev = r.photon_ev[0]     - 0.5 * photon_step;
        r.y_max_ev = r.photon_ev[P - 1] + 0.5 * photon_step;

        // ---- XAS side panel --------------------------------------------------------------------------
        const size_t C = in.num_core_states;
        if (C > 0 && in.core_eigenvalues_au && in.core_osc_strengths) {
            int* xas_sort = md_temp_alloc_array(temp, int, C);
            for (size_t i = 0; i < C; ++i) xas_sort[i] = (int)i;
            std::sort(xas_sort, xas_sort + C, [&](int a, int b) {
                return in.core_eigenvalues_au[a] < in.core_eigenvalues_au[b];
            });

            md_array_resize(r.core_ev, C, alloc);
            md_array_resize(r.core_f,  C, alloc);
            for (size_t i = 0; i < C; ++i) {
                r.core_ev[i] = in.core_eigenvalues_au[xas_sort[i]] * HARTREE_TO_EV;
                r.core_f[i]  = in.core_osc_strengths[xas_sort[i]];
            }
            r.num_core = C;

            const double xas_min = r.core_ev[0]     - 1.0;
            const double xas_max_e = r.core_ev[C - 1] + 1.0;

            size_t M = (size_t)((xas_max_e - xas_min) / 0.01 + 0.5);
            M = CLAMP(M, (size_t)2, (size_t)4096);

            md_array_resize(r.xas_x, M, alloc);
            md_array_resize(r.xas_y, M, alloc);

            // The XAS is broadened with the core-hole lifetime (gamma), not with the map broadening.
            const double xas_fwhm = (in.gamma_fwhm_ev > 0.0) ? in.gamma_fwhm_ev : r.broadening_fwhm_ev;

            double y_hi = 0.0;
            for (size_t j = 0; j < M; ++j) {
                const double x = xas_min + (xas_max_e - xas_min) * (double)j / (double)(M - 1);
                const double y = rixs_broaden(r.broadening_mode, x, r.core_ev, r.core_f, C, xas_fwhm);
                r.xas_x[j] = x;
                r.xas_y[j] = y;
                y_hi = MAX(y_hi, y);
            }

            // The sticks are drawn on the same axis, so make sure they are not clipped either.
            for (size_t i = 0; i < C; ++i) {
                y_hi = MAX(y_hi, r.core_f[i]);
            }
            r.num_xas = M;
            r.xas_max = y_hi;
        }

        r.valid = true;
    }

    // Draws the XAS side panel content: a filled broadened curve plotted as x = f(y), plus the
    // core-excitation sticks as horizontal bars. Only valid inside a plot context.
    void plot_xas_sidebar(const Rixs& r, int& hovered, ImVec4 curve_col, ImVec4 stick_col) {
        ImDrawList& draw_list = *ImPlot::GetPlotDrawList();
        ImPlot::PushPlotClipRect();
        defer { ImPlot::PopPlotClipRect(); };

        // Fill between x = 0 and the curve. ImPlot::PlotShaded only shades along y, so the horizontal
        // equivalent of matplotlib's fill_betweenx is drawn manually, one quad per segment.
        if (r.num_xas > 1) {
            const ImU32 fill_col = ImGui::ColorConvertFloat4ToU32(ImVec4(curve_col.x, curve_col.y, curve_col.z, curve_col.w * 0.25f));
            for (size_t j = 0; j + 1 < r.num_xas; ++j) {
                const ImVec2 a = ImPlot::PlotToPixels(0.0,          r.xas_x[j],     IMPLOT_AUTO, IMPLOT_AUTO);
                const ImVec2 b = ImPlot::PlotToPixels(r.xas_y[j],   r.xas_x[j],     IMPLOT_AUTO, IMPLOT_AUTO);
                const ImVec2 c = ImPlot::PlotToPixels(r.xas_y[j+1], r.xas_x[j + 1], IMPLOT_AUTO, IMPLOT_AUTO);
                const ImVec2 d = ImPlot::PlotToPixels(0.0,          r.xas_x[j + 1], IMPLOT_AUTO, IMPLOT_AUTO);
                draw_list.AddQuadFilled(a, b, c, d, fill_col);
            }
        }

        // Sticks (matplotlib barh), drawn as horizontal bars of a fixed pixel thickness.
        const float bar_thickness_in_pixels = 2.0f;
        hovered = -1;

        if (ImPlot::IsPlotHovered() && r.num_core > 0) {
            const ImVec2 mouse_pos = ImGui::GetMousePos();
            const float  tolerance = 7.0f;
            double min_dist = DBL_MAX;

            for (size_t i = 0; i < r.num_core; ++i) {
                const ImVec2 p0 = ImPlot::PlotToPixels(0.0,          r.core_ev[i], IMPLOT_AUTO, IMPLOT_AUTO);
                const ImVec2 p1 = ImPlot::PlotToPixels(r.core_f[i],  r.core_ev[i], IMPLOT_AUTO, IMPLOT_AUTO);
                const ImVec2 lo = ImMin(p0, p1);
                const ImVec2 hi = ImMax(p0, p1);

                const ImVec2 dx = mouse_pos - ImClamp(mouse_pos, lo, hi);
                const double d2 = dx.x * dx.x + dx.y * dx.y;
                if (d2 < (tolerance * tolerance) && d2 < min_dist) {
                    hovered = (int)i;
                    min_dist = d2;
                }
            }
        }

        for (size_t i = 0; i < r.num_core; ++i) {
            const ImVec4 col = ((int)i == hovered)            ? COLOR_PEAK_HOVER
                             : ((int)i == r.selected_core)   ? COLOR_PEAK_SELECTED
                                                             : stick_col;
            ImVec2 p0 = ImPlot::PlotToPixels(0.0,         r.core_ev[i], IMPLOT_AUTO, IMPLOT_AUTO);
            ImVec2 p1 = ImPlot::PlotToPixels(r.core_f[i], r.core_ev[i], IMPLOT_AUTO, IMPLOT_AUTO);
            p0.y -= bar_thickness_in_pixels * 0.5f;
            p1.y += bar_thickness_in_pixels * 0.5f;
            draw_list.AddRectFilled(ImMin(p0, p1), ImMax(p0, p1), ImGui::ColorConvertFloat4ToU32(col));
        }
    }

    // The settings for the whole RIXS section. Drawn once, above both plots: the spectrum is a single
    // row of the map, so they share the broadening and the elastic line and must not drift apart.
    void draw_rixs_settings(const RixsMapInput& in, Rixs& r) {
        const float avail_width = ImGui::GetContentRegionAvail().x;
        ImGui::PushItemWidth(MIN(avail_width, 200.0f));
        defer { ImGui::PopItemWidth(); };

        static const double broadening_min = 0.01;
        static const double broadening_max = 2.0;
        ImGui::SliderScalar((const char*)u8"Broadening FWHM (eV)", ImGuiDataType_Double, &r.broadening_fwhm_ev, &broadening_min, &broadening_max);
        ImGui::Combo("Broadening mode", (int*)(&r.broadening_mode), broadening_mode_str, BROADENING_MODE_COUNT);

        if (in.elastic_cross_sections) {
            ImGui::Checkbox("Elastic line", &r.plot_elastic_line);
            ImGui::SameLine();
        }
        // Only affects the map, the spectrum is drawn on its own auto fitted axes.
        ImGui::Checkbox("Normalize map", &r.normalize);
    }

    // The RIXS section: one set of settings, then the spectrum and the map, each in its own tree node
    // so either can be collapsed away.
    void draw_rixs_section(const RixsMapInput& in, Rixs& r) {
        // The enclosing tree node uses NoTreePushOnOpen, so without this everything drawn here would
        // share the ID scope of the parent window - and collide with the identically labelled
        // broadening widgets of the XPS section, which sits under the same kind of node.
        ImGui::PushID("rixs");
        defer { ImGui::PopID(); };

        const ImGuiTreeNodeFlags flags = ImGuiTreeNodeFlags_DefaultOpen | ImGuiTreeNodeFlags_NoTreePushOnOpen;

        draw_rixs_settings(in, r);

        ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
        if (ImGui::TreeNodeEx("Spectrum", flags)) {
            draw_rixs_spectrum(in, r);
        }

        ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
        if (ImGui::TreeNodeEx("Map", flags)) {
            draw_rixs_map(in, r);
        }
    }

    // Standalone RIXS spectrum widget: the energy loss spectrum for one incoming photon energy, which
    // is exactly one row of the RIXS map. The photon energy is picked from a combo listing every
    // incoming energy of the calculation, the broadening is shared with the map.
    // Call from inside an ImGui window.
    void draw_rixs_spectrum(const RixsMapInput& in, Rixs& r, ImVec2 size = ImVec2(-1.0f, 350.0f)) {
        const size_t P = in.num_photon_energies;
        const size_t F = in.num_final_states;
        if (P == 0 || F == 0 || !in.photon_energies_au || !in.cross_sections || !in.energy_losses_au) return;

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        constexpr int num_broadened_samples = 2048;

        // Photon energies ascending, so the combo lists them in order and the indices agree with the
        // row order of the map.
        int* sort_idx = md_temp_alloc_array(temp, int, P);
        for (size_t i = 0; i < P; ++i) sort_idx[i] = (int)i;
        std::sort(sort_idx, sort_idx + P, [&](int a, int b) {
            return in.photon_energies_au[a] < in.photon_energies_au[b];
        });

        double* photon_ev = md_temp_alloc_array(temp, double, P);
        for (size_t i = 0; i < P; ++i) {
            photon_ev[i] = in.photon_energies_au[sort_idx[i]] * HARTREE_TO_EV;
        }

        // The selection is kept across data sets, so it has to be validated against the current one.
        r.selected_photon = CLAMP(r.selected_photon, 0, (int)P - 1);

        // ---- Which row of the map to show ------------------------------------------------------------
        // Not a setting of the section: it selects the data, the broadening settings above apply to
        // both plots and are drawn once by draw_rixs_settings().
        {
            const float avail_width = ImGui::GetContentRegionAvail().x;
            ImGui::PushItemWidth(MIN(avail_width, 200.0f));
            defer { ImGui::PopItemWidth(); };

            char preview[32];
            snprintf(preview, sizeof(preview), "%.4f eV", photon_ev[r.selected_photon]);
            if (ImGui::BeginCombo("Photon energy", preview)) {
                for (int i = 0; i < (int)P; ++i) {
                    char item[32];
                    snprintf(item, sizeof(item), "%.4f eV", photon_ev[i]);
                    if (ImGui::Selectable(item, i == r.selected_photon)) {
                        r.selected_photon = i;
                        r.hovered_final   = -1;
                        r.selected_final  = -1;
                    }
                    if (i == r.selected_photon) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }
        }

        // ---- Sticks for the selected photon energy ---------------------------------------------------
        // Column p of the [F][P] arrays, gathered into the tightly packed arrays the broadening kernels
        // and plot_peaks expect. A handful of doubles, from the frame temp arena.
        const size_t p = (size_t)sort_idx[r.selected_photon];
        const bool   use_elastic = r.plot_elastic_line && in.elastic_cross_sections != nullptr;
        const size_t num_peaks   = F + (use_elastic ? 1 : 0);

        double* peaks_x = md_temp_alloc_array(temp, double, num_peaks);
        double* peaks_y = md_temp_alloc_array(temp, double, num_peaks);

        double x_lo =  DBL_MAX;
        double x_hi = -DBL_MAX;
        for (size_t f = 0; f < F; ++f) {
            peaks_x[f] = in.energy_losses_au[f * P + p] * HARTREE_TO_EV;
            peaks_y[f] = in.cross_sections[f * P + p];
            x_lo = MIN(x_lo, peaks_x[f]);
            x_hi = MAX(x_hi, peaks_x[f]);
        }
        if (use_elastic) {
            // Elastic scattering: zero energy loss.
            peaks_x[F] = 0.0;
            peaks_y[F] = in.elastic_cross_sections[p];
            x_lo = MIN(x_lo, peaks_x[F]);
            x_hi = MAX(x_hi, peaks_x[F]);
        }

        if (!(x_lo <= x_hi)) return;  // catches NaN / empty input

        // Same padding as the map, so both cover the same range for a given photon energy.
        const double pad = 0.5;
        const double x_min = x_lo - pad;
        const double x_max = x_hi + pad;

        // ---- Refit whenever the plotted data changes -------------------------------------------------
        uint64_t hash = md_hash64(&r.selected_photon,   sizeof(r.selected_photon),   0);
        hash = md_hash64(&r.plot_elastic_line,  sizeof(r.plot_elastic_line),  hash);
        hash = md_hash64(&r.broadening_mode,    sizeof(r.broadening_mode),    hash);
        hash = md_hash64(&r.broadening_fwhm_ev, sizeof(r.broadening_fwhm_ev), hash);
        hash = md_hash64(&in,                   sizeof(in),                   hash);

        if (hash != r.spectrum_hash) {
            r.spectrum_hash = hash;
            ImPlot::SetNextAxesToFit();
        }

        if (ImPlot::BeginPlot("RIXS spectrum", size)) {
            // The energy-loss convention draws the axis reversed (high loss on the left), as in the map.
            ImPlot::SetupAxis(ImAxis_X1, "Energy loss [eV]", ImPlotAxisFlags_Invert);
            ImPlot::SetupAxis(ImAxis_Y1, "Intensity [a.u.]");
            ImPlot::SetupAxis(ImAxis_Y2, "Cross section [a.u.]", ImPlotAxisFlags_AuxDefault);
            ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
            ImPlot::SetupFinish();

            BroadenedCurve curve = {
                peaks_x, peaks_y, num_peaks, x_min, x_max,
                r.broadening_fwhm_ev, num_broadened_samples, r.broadening_mode,
            };
            ImPlot::SetAxis(ImAxis_Y1);
            ImPlot::PlotLineG("Broadened Spectrum", broadened_curve_getter, &curve, num_broadened_samples);

            ImPlot::SetAxis(ImAxis_Y2);
            plot_peaks("Cross Section", peaks_x, peaks_y, num_peaks, r.selected_final, r.hovered_final);

            if (ImPlot::IsPlotHovered() && r.hovered_final != -1) {
                const int i = r.hovered_final;
                if (ImGui::BeginTooltip()) {
                    if (use_elastic && (size_t)i == F) {
                        ImGui::TextUnformatted("Elastic line");
                    } else {
                        ImGui::Text("Final state %i", i + 1);
                    }
                    ImGui::Text("Energy loss: %.4f eV", peaks_x[i]);
                    ImGui::Text("Cross section: %.4g", peaks_y[i]);
                    ImGui::EndTooltip();
                }
            }

            ImPlot::EndPlot();
        }
    }

    // RIXS map widget: 2D map, XAS side panel and colorbar. The settings it reads are drawn by
    // draw_rixs_settings().
    // Call from inside an ImGui window. 'size' is the total size of the plot row; pass -1 for the
    // width to use the available content region.
    void draw_rixs_map(const RixsMapInput& in, Rixs& r, ImVec2 size = ImVec2(-1.0f, 350.0f)) {
        // ---- Rebuild the cache only when something actually changed -----------------------------------
        uint64_t hash = md_hash64(&r.plot_elastic_line, sizeof(r.plot_elastic_line), 0);
        hash = md_hash64(&r.normalize,          sizeof(r.normalize),          hash);
        hash = md_hash64(&r.broadening_mode,    sizeof(r.broadening_mode),    hash);
        hash = md_hash64(&r.broadening_fwhm_ev, sizeof(r.broadening_fwhm_ev), hash);
        hash = md_hash64(&r.x_step_ev,          sizeof(r.x_step_ev),          hash);
        hash = md_hash64(&in,                   sizeof(in),                   hash);

        const bool refit = (hash != r.hash);
        if (refit) {
            r.hash = hash;
            compute_rixs_map(r, in);
        }

        if (!r.valid) {
            ImGui::TextUnformatted("No RIXS map could be computed from the loaded data.");
            return;
        }

        // ---- Layout: map | XAS | colorbar --------------------------------------------------------------
        const float spacing       = ImGui::GetStyle().ItemSpacing.x;
        const float colorbar_w    = 90.0f;
        const float height        = (size.y > 0.0f) ? size.y : 350.0f;
        const float total_w       = ((size.x > 0.0f) ? size.x : ImGui::GetContentRegionAvail().x) - colorbar_w - 2.0f * spacing;
        const bool  has_xas       = r.num_xas > 1;
        // Same 4 : 1.25 width ratio as the matplotlib gridspec.
        const float map_w         = has_xas ? total_w * (4.0f / 5.25f) : total_w;
        const float xas_w         = total_w - map_w;

        const char* x_label = "Energy loss [eV]";

        // The energy-loss convention draws the axis reversed (high loss on the left).
        const ImPlotAxisFlags x_flags = ImPlotAxisFlags_Invert;

        // Only the heatmap and the colorbar are drawn with the map colormap. The XAS panel is a
        // regular line plot and picks up the default item colors, like every other spectrum.
        ImPlot::PushColormap(r.colormap);

        // PlotHeatmap only rasterizes the cells it was given, so anything outside the data would fall
        // through to the window background. Painting the canvas in the color a zero sample maps to
        // makes the map read as an infinite field of zeros instead of a floating tile.
        const double map_range = MAX(r.map_max - r.map_min, DBL_EPSILON);
        const float  zero_t    = ImClamp((float)((0.0 - r.map_min) / map_range), 0.0f, 1.0f);
        ImPlot::PushStyleColor(ImPlotCol_PlotBg, ImPlot::SampleColormap(zero_t, r.colormap));

        ImPlot::SetNextAxisLinks(ImAxis_Y1, &r.y_link_min, &r.y_link_max);
        if (ImPlot::BeginPlot("RIXS map", ImVec2(map_w, height), ImPlotFlags_NoLegend)) {
            ImPlot::SetupAxis(ImAxis_X1, x_label, x_flags);
            ImPlot::SetupAxis(ImAxis_Y1, "Photon energy [eV]");
            ImPlot::SetupAxesLimits(r.x_min_ev, r.x_max_ev, r.y_min_ev, r.y_max_ev, refit ? ImPlotCond_Always : ImPlotCond_Once);

            // Keep panning and zooming inside the data. Combined with the background above, there is
            // no way to end up looking at a region the map says nothing about.
            ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, r.x_min_ev, r.x_max_ev);
            ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, r.y_min_ev, r.y_max_ev);
            ImPlot::SetupFinish();

            ImPlot::PlotHeatmap("##rixs_map", r.map, (int)r.num_rows, (int)r.num_cols,
                                r.map_min, r.map_max, nullptr,
                                ImPlotPoint(r.x_min_ev, r.y_min_ev),
                                ImPlotPoint(r.x_max_ev, r.y_max_ev));

            // The row currently shown in the spectrum. White reads against the whole colormap, which
            // a colormap color would not.
            if (r.selected_photon >= 0 && (size_t)r.selected_photon < r.num_rows) {
                const double y = r.photon_ev[r.selected_photon];
                ImPlot::SetNextLineStyle(ImVec4(1.0f, 1.0f, 1.0f, 0.85f), 1.5f);
                ImPlot::PlotInfLines("##selected_photon", &y, 1, ImPlotInfLinesFlags_Horizontal);
            }

            // Clicking anywhere in the map moves the spectrum to the closest incoming photon energy.
            // Released rather than clicked, and dragging excluded, so that panning the plot does not
            // change the selection.
            if (ImPlot::IsPlotHovered() && ImGui::IsMouseReleased(ImGuiMouseButton_Left) &&
                !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Left) && r.num_rows > 0) {
                const double y = ImPlot::GetPlotMousePos().y;

                int    closest   = 0;
                double closest_d = DBL_MAX;
                for (size_t i = 0; i < r.num_rows; ++i) {
                    const double d = fabs(r.photon_ev[i] - y);
                    if (d < closest_d) {
                        closest_d = d;
                        closest   = (int)i;
                    }
                }

                if (closest != r.selected_photon) {
                    r.selected_photon = closest;
                    r.hovered_final   = -1;
                    r.selected_final  = -1;
                }
            }

            if (ImPlot::IsPlotHovered()) {
                const ImPlotPoint p = ImPlot::GetPlotMousePos();

                // Map the cursor back onto a cell to read out the intensity.
                const double u = (p.x - r.x_min_ev) / (r.x_max_ev - r.x_min_ev);
                const double v = (p.y - r.y_min_ev) / (r.y_max_ev - r.y_min_ev);
                if (u >= 0.0 && u < 1.0 && v >= 0.0 && v < 1.0) {
                    const size_t col = (size_t)(u * (double)r.num_cols);
                    const size_t row = r.num_rows - 1 - (size_t)(v * (double)r.num_rows);
                    ImGui::BeginTooltip();
                    ImGui::Text("Energy loss: %.3f eV", p.x);
                    ImGui::Text("Photon energy: %.3f eV", p.y);
                    ImGui::Text("Intensity: %.4g", r.map[row * r.num_cols + col]);
                    ImGui::EndTooltip();
                }
            }
            ImPlot::EndPlot();
        }
        ImPlot::PopStyleColor();
        ImPlot::PopColormap();

        if (has_xas) {
            ImGui::SameLine();
            ImPlot::SetNextAxisLinks(ImAxis_Y1, &r.y_link_min, &r.y_link_max);
            const ImPlotFlags xas_flags = ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText;
            if (ImPlot::BeginPlot("XAS", ImVec2(xas_w, height), xas_flags)) {
                ImPlot::SetupAxis(ImAxis_X1, "Intensity",  ImPlotAxisFlags_NoGridLines);
                // The photon energy axis is shared with the map, so it is not labelled again here.
                ImPlot::SetupAxis(ImAxis_Y1, nullptr, ImPlotAxisFlags_NoTickLabels | ImPlotAxisFlags_NoLabel);
                ImPlot::SetupAxesLimits(0.0, MAX(r.xas_max, 1.0e-12) * 1.1, r.y_min_ev, r.y_max_ev, refit ? ImPlotCond_Always : ImPlotCond_Once);
                // The photon energy axis is linked to the map, so it carries the same constraint.
                ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, r.y_min_ev, r.y_max_ev);
                ImPlot::SetupFinish();

                // The first two colors of the current colormap, as two ordinary plot items would get.
                // NOTE: GetColormapColor, not NextColormapColor: the latter advances a counter that
                // lives in the plot and is never reset per frame, so calling it for items that are not
                // registered with ImPlot walks through the colormap one step per frame (flickering).
                const ImVec4 curve_col = ImPlot::GetColormapColor(0);
                const ImVec4 stick_col = ImPlot::GetColormapColor(1);

                plot_xas_sidebar(r, r.hovered_core, curve_col, stick_col);

                ImPlot::SetNextLineStyle(curve_col);
                ImPlot::PlotLine("##xas", r.xas_y, r.xas_x, (int)r.num_xas);

                if (r.hovered_core != -1) {
                    ImGui::BeginTooltip();
                    ImGui::Text("Core state %i", r.hovered_core + 1);
                    ImGui::Text("Energy: %.3f eV", r.core_ev[r.hovered_core]);
                    ImGui::Text("Osc. strength: %.4g", r.core_f[r.hovered_core]);
                    ImGui::EndTooltip();
                }

                if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Left) && ImPlot::IsPlotHovered()) {
                    r.selected_core = (r.hovered_core == r.selected_core) ? -1 : r.hovered_core;
                }

                ImPlot::EndPlot();
            }
        }

        ImGui::SameLine();
        ImPlot::PushColormap(r.colormap);
        ImPlot::ColormapScale("Intensity [arb. u.]", r.map_min, r.map_max, ImVec2(colorbar_w, height),
                              r.normalize ? "%.2f" : "%.1e");
        ImPlot::PopColormap();
    }

    // =================================================================================================
    // XPS
    //
    // ImGui/ImPlot port of VeloxChem's spectrumplot.plot_xps_spectrum(). Like the python version this
    // draws one element at a time, picked from a combo box listing only the elements that actually
    // have core hole entries.
    //
    // Within that element every peak is one ATOM, so the plot doubles as a way to find the atom a
    // binding energy belongs to: hovering a peak highlights its atom in the viewport and clicking
    // adds it to the selection. That is also why there is no CPK coloring any more - the colour used
    // to encode the element, which the combo box now carries, and per atom CPK would encode nothing
    // at all since every peak on screen shares one element.
    //
    // Two ImPlot items, matching the other spectra in the Response window: "Broadened Spectrum" on Y1
    // and "Peaks" on Y2, each with its own legend entry and show/hide state.
    // =================================================================================================

    // Per-element input to the broadening getter. Lives on the stack for the duration of one
    // PlotLineG call, which consumes it immediately.
    // A broadened curve sampled on the fly by ImPlot, shared by the XPS and RIXS spectra.
    struct BroadenedCurve {
        const double*     peaks_x;
        const double*     peaks_y;
        size_t            num_peaks;
        double            x_min;
        double            x_max;
        double            fwhm;
        int               num_samples;
        broadening_mode_t mode;
    };

    static ImPlotPoint broadened_curve_getter(int idx, void* user_data) {
        const BroadenedCurve* c = (const BroadenedCurve*)user_data;
        const double x = c->x_min + (c->x_max - c->x_min) * idx / (double)(c->num_samples - 1);
        const double y = (c->mode == BROADENING_MODE_LORENTZIAN)
            ? lorentzian_xps(x, c->peaks_x, c->peaks_y, c->num_peaks, c->fwhm)
            : gaussian_xps  (x, c->peaks_x, c->peaks_y, c->num_peaks, c->fwhm);
        return ImPlotPoint{ x, y };
    }

    // Draws the stick spectrum of one element as an ImPlot item under 'label'.
    // 'hovered' is an index into the flat entry array and 'base' is the offset of this group within
    // it, so a hover survives being compared against entries from any group.
    // 'stick_selected' is per entry and comes from the application selection mask - the plot does not
    // own a selection of its own.
    // Returns false if the item is hidden through the legend, in which case nothing was drawn.
    static bool plot_xps_sticks(const char* label, const md_vlx_xps_entry_t* entries, const double* stick_y,
                                const bool* stick_selected, size_t count, size_t base, int& hovered) {
        const float bar_width_in_pixels    = 2.0f;
        const float point_radius_in_pixels = 3.0f;

        if (!ImPlot::BeginItem(label)) {
            // Hidden via the legend. Clear any hover that belonged to this element so a stale
            // tooltip does not survive the toggle.
            if (hovered >= (int)base && hovered < (int)(base + count)) hovered = -1;
            return false;
        }
        defer { ImPlot::EndItem(); };

        ImPlot::PushPlotClipRect();
        defer { ImPlot::PopPlotClipRect(); };

        ImDrawList& draw_list = *ImPlot::GetPlotDrawList();
        const ImPlotNextItemData& s = ImPlot::GetItemData();
        const ImVec4 item_color = s.Colors[ImPlotCol_Fill];

        const bool  plot_hovered = ImPlot::IsPlotHovered();
        const ImVec2 mouse_pos   = ImGui::GetMousePos();
        const double min_line_d2  = pow(bar_width_in_pixels * 0.5 + 7.0, 2);
        const double min_point_d2 = pow(point_radius_in_pixels + 7.0, 2);
        double min_dist = DBL_MAX;

        for (size_t i = 0; i < count; ++i) {
            const double x = entries[i].ionization_energy;
            const double y = stick_y[i];

            const ImVec2 p0_raw = ImPlot::PlotToPixels(x, 0.0, IMPLOT_AUTO, IMPLOT_AUTO);
            const ImVec2 p1_raw = ImPlot::PlotToPixels(x, y,   IMPLOT_AUTO, IMPLOT_AUTO);
            const ImVec2 p_min  = ImMin(p0_raw, p1_raw);
            const ImVec2 p_max  = ImMax(p0_raw, p1_raw);

            if (plot_hovered) {
                const ImVec2 line_dx = mouse_pos - ImClamp(mouse_pos, p_min, p_max);
                const double line_d2 = line_dx.x * line_dx.x + line_dx.y * line_dx.y;
                if (line_d2 < min_line_d2 && line_d2 < min_dist) {
                    hovered  = (int)(base + i);
                    min_dist = line_d2;
                }
                const ImVec2 point_dx = mouse_pos - p1_raw;
                const double point_d2 = point_dx.x * point_dx.x + point_dx.y * point_dx.y;
                if (point_d2 < min_point_d2 && point_d2 < min_dist) {
                    hovered  = (int)(base + i);
                    min_dist = point_d2;
                }
            }

            const int    flat_idx = (int)(base + i);
            const ImVec4 color    = (flat_idx == hovered)   ? COLOR_PEAK_HOVER
                                  : (stick_selected[i])     ? COLOR_PEAK_SELECTED
                                                            : item_color;
            const ImU32  col32    = ImGui::ColorConvertFloat4ToU32(color);

            draw_list.AddRectFilled(ImVec2{ p_min.x - bar_width_in_pixels * 0.5f, p_min.y },
                                    ImVec2{ p_max.x + bar_width_in_pixels * 0.5f, p_max.y }, col32);
            draw_list.AddCircleFilled(p1_raw, point_radius_in_pixels, col32);

            // A delocalized core hole is not attributable to a single atom, so mark it.
            if (entries[i].is_delocalized) {
                draw_list.AddCircle(p1_raw, point_radius_in_pixels + 2.5f, col32, 0, 1.5f);
            }

            if (ImPlot::FitThisFrame()) {
                ImPlot::FitPoint(ImPlotPoint(x, 0.0));
                ImPlot::FitPoint(ImPlotPoint(x, y));
            }
        }

        return true;
    }

    // Maps an atom index as it appears in the vlx object onto an index into the loaded system. The
    // two differ when the qm data covers only a subset of the system. Returns -1 when the entry has
    // no atom at all - a delta-SCF entry can arrive without one - or when the mapping lands outside
    // the system, in which case there is nothing to highlight.
    int32_t xps_atom_index(int32_t vlx_atom_index, const ApplicationState& state) const {
        if (vlx_atom_index < 0) return -1;
        if ((size_t)vlx_atom_index >= md_vlx_number_of_atoms(vlx)) return -1;
        const int32_t idx = qm_to_atom_idx ? qm_to_atom_idx[vlx_atom_index] : vlx_atom_index;
        return (idx >= 0 && (size_t)idx < state.mold.sys.atom.count) ? idx : -1;
    }

    void draw_xps_plot(ApplicationState& state, ImVec2 size = ImVec2(-1.0f, 350.0f)) {
        const size_t num_groups  = md_vlx_xps_group_count(vlx);
        const size_t num_entries = md_vlx_xps_count(vlx);
        if (num_groups == 0 || num_entries == 0) return;

        // The enclosing tree node uses NoTreePushOnOpen, so this section would otherwise draw into the
        // parent window's ID scope. Two collisions follow from that: ImPlot derives the plot ID from
        // Window->GetID(title), so BeginPlot("XPS") lands on the exact ID of the TreeNodeEx("XPS")
        // header - which makes ImPlot report the plot as held while the header is dragged, panning the
        // axes - and the broadening widgets below collide with the identically labelled ones in the
        // RIXS section when a file carries both.
        ImGui::PushID("xps");
        defer { ImGui::PopID(); };

        const md_vlx_xps_entry_t* all_entries = md_vlx_xps_entries(vlx);

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        constexpr int num_broadened_samples = 2048;

        // Resolves the selected element, adopting the first non empty group when the selection is
        // unset or no longer present. Called again after the combo so a change lands on the same
        // frame it is made rather than one frame later.
        auto resolve_group = [&]() -> const md_vlx_xps_group_t* {
            const md_vlx_xps_group_t* g = md_vlx_xps_group_by_element(vlx, xps.element);
            if (g && g->count > 0) return g;
            for (size_t i = 0; i < num_groups; ++i) {
                g = md_vlx_xps_group_by_index(vlx, i);
                if (g && g->count > 0) {
                    xps.element = g->element;
                    return g;
                }
            }
            return nullptr;
        };

        const md_vlx_xps_group_t* grp = resolve_group();
        if (!grp) return;

        // ---- Settings ---------------------------------------------------------------------------
        const float avail_width = ImGui::GetContentRegionAvail().x;
        ImGui::PushItemWidth(MIN(avail_width, 200.0f));
        defer { ImGui::PopItemWidth(); };

        // The group list is exactly the set of elements that were computed, so listing it offers no
        // element that would draw an empty plot.
        {
            const str_t cur_sym = md_util_element_symbol(grp->element);
            char preview[8];
            snprintf(preview, sizeof(preview), "%.*s", (int)cur_sym.len, cur_sym.ptr);
            if (ImGui::BeginCombo("Element", preview)) {
                for (size_t g = 0; g < num_groups; ++g) {
                    const md_vlx_xps_group_t* it = md_vlx_xps_group_by_index(vlx, g);
                    if (!it || it->count == 0) continue;
                    const str_t sym = md_util_element_symbol(it->element);
                    char item[32];
                    snprintf(item, sizeof(item), "%.*s (%i)", (int)sym.len, sym.ptr, (int)it->count);
                    if (ImGui::Selectable(item, it->element == xps.element)) {
                        xps.element = it->element;
                        xps.hovered = -1;
                    }
                }
                ImGui::EndCombo();
            }
            // NOTE: no tooltip on the combo. IsItemHovered() after EndCombo() reports on whatever the
            // popup left in LastItemData, not on the combo frame. The bracketed count is the number
            // of core hole states computed for that element, one per atom.
            grp = resolve_group();
            if (!grp) return;
        }

        static const double fwhm_min = 0.05;
        static const double fwhm_max = 5.0;
        ImGui::SliderScalar((const char*)u8"Broadening FWHM (eV)", ImGuiDataType_Double, &xps.broadening_fwhm_ev, &fwhm_min, &fwhm_max);
        ImGui::Combo("Broadening mode", (int*)(&xps.broadening_mode), broadening_mode_str, BROADENING_MODE_COUNT);

        ImGui::Checkbox("Sticks", &xps.show_sticks);
        ImGui::SameLine();
        ImGui::Checkbox("Invert energy axis", &xps.invert_x);
        ImGui::SameLine();
        ImGui::Checkbox("Weight by contribution", &xps.weight_by_contribution);
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Off: every core hole contributes a unit peak (what VeloxChem plots).\n"
                              "On:  each peak is scaled by the atom's share of the core orbital, so a\n"
                              "     hole delocalized over N atoms sums to one instead of N.");
        }

        // ---- Per peak arrays ---------------------------------------------------------------------
        // The broadening kernels take tightly packed arrays, so the strided entry fields are gathered
        // here. A handful of doubles for one element, from the frame temp arena.
        double* peaks_x        = md_temp_alloc_array(temp, double, grp->count);
        double* peaks_y        = md_temp_alloc_array(temp, double, grp->count);
        bool*   peaks_selected = md_temp_alloc_array(temp, bool,   grp->count);
        for (size_t i = 0; i < grp->count; ++i) {
            peaks_x[i] = grp->entries[i].ionization_energy;
            peaks_y[i] = xps.weight_by_contribution ? grp->entries[i].contribution : 1.0;
            const int32_t atom = xps_atom_index(grp->entries[i].atom_index, state);
            peaks_selected[i]  = (atom >= 0) && md_bitfield_test_bit(&state.selection.selection_mask, (uint64_t)atom);
        }

        // ---- Energy range ------------------------------------------------------------------------
        // Over the selected element only. Core binding energies of different elements sit hundreds of
        // eV apart, so a range spanning all of them would compress the peaks actually on screen into
        // a single unreadable spike.
        double x_min = DBL_MAX;
        double x_max = -DBL_MAX;
        for (size_t i = 0; i < grp->count; ++i) {
            x_min = MIN(x_min, peaks_x[i]);
            x_max = MAX(x_max, peaks_x[i]);
        }
        const double pad = 5.0;
        x_min = MAX(0.0, x_min - pad);
        x_max = x_max + pad;

        uint64_t hash = md_hash64(&xps.broadening_fwhm_ev, sizeof(xps.broadening_fwhm_ev), xps.broadening_mode);
        hash = md_hash64(&xps.weight_by_contribution, sizeof(xps.weight_by_contribution), hash);
        hash = md_hash64(&xps.element, sizeof(xps.element), hash);
        hash = md_hash64(grp->entries, grp->count * sizeof(md_vlx_xps_entry_t), hash);
        const bool refit = (hash != xps.hash);
        xps.hash = hash;

        if (refit) {
            ImPlot::SetNextAxesToFit();
        }

        const ImPlotAxisFlags x_flags = xps.invert_x ? ImPlotAxisFlags_Invert : ImPlotAxisFlags_None;

        // NoMenus: right click removes an atom from the selection here, and ImPlot would otherwise
        // open its context menu on the very same event. The ramachandran plot gives up its context
        // menu for the same reason.
        if (ImPlot::BeginPlot("XPS", size, ImPlotFlags_NoMenus)) {
            ImPlot::SetupAxis(ImAxis_X1, "Binding energy (eV)", x_flags);
            ImPlot::SetupAxis(ImAxis_Y1, "Intensity (a.u.)");
            if (xps.show_sticks) {
                ImPlot::SetupAxis(ImAxis_Y2, "Peak intensity", ImPlotAxisFlags_AuxDefault);
            }
            ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
            ImPlot::SetupFinish();

            // Recomputed from scratch every frame by plot_xps_sticks. Cleared unconditionally so that
            // hiding the sticks - through the checkbox or the legend - cannot leave a stale hover
            // driving the tooltip and the highlight.
            xps.hovered = -1;

            BroadenedCurve curve = {
                peaks_x, peaks_y, grp->count, x_min, x_max,
                xps.broadening_fwhm_ev, num_broadened_samples, xps.broadening_mode,
            };
            ImPlot::SetAxis(ImAxis_Y1);
            ImPlot::PlotLineG("Broadened Spectrum", broadened_curve_getter, &curve, num_broadened_samples);

            if (xps.show_sticks) {
                ImPlot::SetAxis(ImAxis_Y2);
                const size_t base = (size_t)(grp->entries - all_entries);
                plot_xps_sticks("Peaks", grp->entries, peaks_y, peaks_selected, grp->count, base, xps.hovered);
            }

            // ---- Interaction ---------------------------------------------------------------------
            const bool plot_hovered = ImPlot::IsPlotHovered();
            const md_vlx_xps_entry_t* hov_entry =
                (xps.hovered >= 0 && xps.hovered < (int)num_entries) ? &all_entries[xps.hovered] : nullptr;
            const int32_t hov_atom = hov_entry ? xps_atom_index(hov_entry->atom_index, state) : -1;

            // While the cursor is inside the plot this owns the highlight mask, which is the contract
            // every other hover provider in the application follows. Cleared unconditionally, so
            // moving off a peak but staying inside the plot drops the highlight rather than leaving
            // the last atom lit up.
            if (plot_hovered) {
                md_bitfield_clear(&state.selection.highlight_mask);
                if (hov_atom >= 0) {
                    md_bitfield_set_bit(&state.selection.highlight_mask, (uint64_t)hov_atom);
                }
            }

            // Same contract towards the orbital grid window: only written while the cursor is inside
            // the plot, so moving off a peak drops the row highlight exactly as it drops the atom
            // highlight. A delocalized hole names one shared core MO, so every entry of that hole
            // lights the same row - which is the truth of it.
            // md_vlx.c range checks atom_index on load but not mo_index, so it is checked here.
            if (plot_hovered && hov_entry && hov_entry->mo_index >= 0 &&
                (size_t)hov_entry->mo_index < num_molecular_orbitals()) {
                xps.highlight_mo_idx = hov_entry->mo_index;
            }

            if (plot_hovered && hov_entry) {
                const str_t sym = md_util_element_symbol(hov_entry->element);
                if (ImGui::BeginTooltip()) {
                    // Atom index is shown 1-based to match the labels VeloxChem prints.
                    if (hov_entry->atom_index >= 0) {
                        ImGui::Text("%.*s%i", (int)sym.len, sym.ptr, hov_entry->atom_index + 1);
                    } else {
                        ImGui::Text("%.*s (unassigned atom)", (int)sym.len, sym.ptr);
                    }
                    ImGui::Text("Binding energy: %.3f eV", hov_entry->ionization_energy);
                    ImGui::Text("Contribution: %.4g", hov_entry->contribution);
                    ImGui::Text("MO index: %i", hov_entry->mo_index);
                    if (hov_entry->is_delocalized) {
                        // The hole is shared with the other atoms carrying this same mo_index; only
                        // the atom this entry names is highlighted.
                        ImGui::TextUnformatted("Delocalized core hole");
                    }
                    ImGui::EndTooltip();
                }
            }

            // Left adds the atom to the selection, right removes it - the same convention the
            // ramachandran plot and the dataset window use. Drags are excluded so that panning or
            // box selecting does not also grab whatever sat under the cursor when the drag began.
            if (plot_hovered && hov_atom >= 0) {
                if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Left)) {
                    md_bitfield_set_bit(&state.selection.selection_mask, (uint64_t)hov_atom);
                } else if (ImGui::IsMouseReleased(ImGuiMouseButton_Right) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Right)) {
                    md_bitfield_clear_bit(&state.selection.selection_mask, (uint64_t)hov_atom);
                }
            }

            ImPlot::EndPlot();
        }
    }

    void set_atom_coordinates(ApplicationState& state, const dvec3_t* atom_coords) {
        const dvec3_t* coords = atom_coords ? atom_coords : md_vlx_atom_coordinates(vlx);
        size_t num_atoms = md_vlx_number_of_atoms(vlx);
        for (size_t i = 0; i < num_atoms; i++) {
            state.mold.state.x[i] = (float)coords[i].x;
            state.mold.state.y[i] = (float)coords[i].y;
            state.mold.state.z[i] = (float)coords[i].z;
        }
		viamd::event_system_broadcast_event(viamd::EventType_ViamdSystemStateChanged, viamd::EventPayloadType_ApplicationState, &state);
        state.mold.dirty_gpu_buffers |= MolBit_ClearVelocity | MolBit_DirtyPosition;
    }

    // Reads one rank 1 series of doubles out of the system's attribute table.
    //
    // Hands back the stored buffer rather than extracting to float, deliberately: the vlx/ series
    // are resident F64 and this panel prints total energies to ten decimals, which a float does not
    // carry. That is the opt-in fast path - it is only valid because the attribute is resident, and
    // the null return covers every case where it is not.
    static const double* vlx_series(size_t* out_count, const md_system_t& sys, str_t path) {
        const md_attribute_t* attr = md_attributes_find(&sys.attributes, path);
        if (!attr || !attr->data) return nullptr;
        if (attr->format.type != MD_ATTRIBUTE_TYPE_F64) return nullptr;
        if (attr->format.rank != 1 || md_attribute_components(&attr->format) != 1) return nullptr;

        if (out_count) *out_count = attr->format.shape[0];
        return (const double*)attr->data;
    }

    void draw_summary_window(ApplicationState& state) {
        if (!summary.show_window) {
            opt.selected = (int)md_vlx_opt_number_of_steps(vlx) - 1;
            return;
        }

        static ImPlotRect lims{ 0,1,0,1 };

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        md_vlx_rsp_type_t rsp_type = md_vlx_rsp_type(vlx);

        // The actual plot
        ImGui::SetNextWindowSize({ 300, 350 }, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Summary", &summary.show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::TreeNode("Level of Calculation")) {
                str_t basis_set = md_vlx_basis_set_ident(vlx);
                str_t dft_func  = md_vlx_dft_func_label(vlx);
                ImGui::Text("Method: %s", str_ptr(dft_func));
                ImGui::Text("SCF Type: %s", vlx_scf_type_str(md_vlx_scf_type(vlx)));
                ImGui::Text("Basis Set: %s", str_ptr(basis_set));
                ImGui::Spacing();
                
                ImGui::TreePop();
            }
            if (ImGui::TreeNode("System Information")) {
                ImGui::Text("Num Atoms:           %-6zu", md_vlx_number_of_atoms(vlx));
                ImGui::Text("Num Alpha Electrons: %-6zu", md_vlx_number_of_electrons(vlx, MD_VLX_SPIN_ALPHA));
                ImGui::Text("Num Beta Electrons:  %-6zu", md_vlx_number_of_electrons(vlx, MD_VLX_SPIN_BETA));
                ImGui::Text("Molecular Charge:    %-6f",  md_vlx_molecular_charge(vlx));
                ImGui::Text("Spin Multiplicity:   %-6zu", md_vlx_spin_multiplicity(vlx));
                if (rsp_type == MD_VLX_RSP_C6) {
                    ImGui::Text("C6 Value:            %-12.6f (au)", md_vlx_c6_value(vlx));
                }
                ImGui::Spacing();
                ImGui::TreePop();
            }

            if (ImGui::TreeNode("SCF")) {
                // The convergence history is read from the system's attribute table rather than
                // from the vlx object. md_vlx_publish_attributes put it there when the file was
                // loaded, so this panel is one consumer among others and a file which carried no
                // history simply has no such path - there is no separate "was it parsed" question.
                size_t num_iter = 0;
                size_t num_grad = 0;
                const double* energy    = vlx_series(&num_iter, state.mold.sys, STR_LIT("vlx/scf/history/energy"));
                const double* grad_norm = vlx_series(&num_grad, state.mold.sys, STR_LIT("vlx/scf/history/gradient_norm"));

                // Siblings in one group share a shape; a file which somehow disagrees is not a
                // history anyone can plot.
                if (!energy || !grad_norm || num_grad != num_iter) {
                    energy = nullptr;
                    grad_norm = nullptr;
                    num_iter = 0;
                }

                double* energy_offsets = nullptr;
                double ref_energy = 0.0;

                // We set up iterations as doubles for easier use
                if (num_iter > 0) {
                    energy_offsets = md_temp_alloc_array(temp, double, num_iter);
                    ref_energy = energy[num_iter - 1];
                    for (size_t i = 0; i < num_iter; i++) {
                        energy_offsets[i] = fabs(energy[i] - ref_energy);
                    }
                }

                if (num_iter > 1) {
                    if (ImPlot::BeginPlot("SCF")) {
                        ImPlot::SetupAxisLimits(ImAxis_X1, 1.0, (int)num_iter);
                        ImPlot::SetupLegend(ImPlotLocation_NorthEast);
                        ImPlot::SetupAxes("Iteration", "Gradient Norm (au)");
                        // We draw 2 y axis as "Energy total" has values in a different range then the rest of the data
                        ImPlot::SetupAxis(ImAxis_Y2, "Energy (au)", ImPlotAxisFlags_AuxDefault);
                        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

                        ImPlot::PlotLine("Gradient", grad_norm, (int)num_iter, 1.0, 1.0);
                        ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
                        ImPlot::PlotLine("Energy", energy_offsets, (int)num_iter, 1.0, 1.0);
                        lims = ImPlot::GetPlotLimits(ImAxis_X1, ImAxis_Y1);
                        //ImPlot::PlotLine("Density Change", iter, vlx.scf.iter.density_change, (int)vlx.scf.iter.count);
                        //ImPlot::PlotLine("Energy Change", iter, vlx.scf.iter.energy_change, (int)vlx.scf.iter.count);
                        //ImPlot::PlotLine("Max Gradient", iter, vlx.scf.iter.max_gradient, (int)vlx.scf.iter.count);
                        ImPlot::EndPlot();
                    }
                } else {
                    ImGui::Text("There is no history in the supplied veloxchem data");
                }
                ImGui::Spacing();
                if (num_iter > 0) {
                    ImGui::Text("Total energy:              %16.10f (au)", energy[num_iter - 1]);
                    ImGui::Text("Gradient norm:             %16.10f (au)", grad_norm[num_iter - 1]);
                }
                ImGui::Text("Nuclear repulsion energy:  %16.10f (au)", md_vlx_nuclear_repulsion_energy(vlx));
                ImGui::Spacing();
                ImGui::TreePop();
            }

            {
                size_t num_steps = md_vlx_opt_number_of_steps(vlx);
                if (num_steps > 0) {
                    if (ImGui::TreeNode("Optimization")) {
                        ImPlot::PushStyleVar(ImPlotStyleVar_FitPadding, ImVec2(0.05f, 0.1f));
                        defer { ImPlot::PopStyleVar(); };

                        md_vlx_opt_type_t opt_type = md_vlx_opt_type(vlx);
                        size_t ts_index = md_vlx_opt_irc_ts_index(vlx);
                        const double* energies = md_vlx_opt_energies(vlx);
                        const char* y_axis_label = "Relative Energy [kJ/mol]";
                        const char* plot_label = "OPT";

                        uint64_t energy_hash = md_hash64((const uint8_t*)energies, sizeof(double) * num_steps, 0);

                        if (energy_hash != opt.energy_hash) {
                            opt.energy_hash = energy_hash;
                            ImPlot::SetNextAxesToFit();
                        }

                        // Find minima as the reference energy
                        double ref_energy = energies[0];
                        if (opt_type == MD_VLX_OPT_IRC) {
                            if (ts_index < num_steps) {
                                ref_energy = energies[ts_index];
                            }
                            plot_label = "Intrinsic Reaction Coordinate (IRC)";
                        } else {
                            for (size_t i = 1; i < num_steps; i++) {
                                ref_energy = MIN(ref_energy, energies[i]);
                            }
                        }

                        double* x_vals = md_temp_alloc_array(temp, double, num_steps);
                        double* y_vals = md_temp_alloc_array(temp, double, num_steps);

                        for (size_t i = 0; i < num_steps; ++i) {
                            y_vals[i] = (energies[i] - ref_energy) * HARTREE_TO_KJ_PER_MOL;
                            x_vals[i] = (double)(i + 1);
                        }

                        const ImPlotFlags plot_flags = ImPlotFlags_NoMouseText;
                        if (ImPlot::BeginPlot(plot_label, ImVec2(-1, 0), plot_flags)) {
                            
                            ImPlot::SetupAxes("Step", y_axis_label);
                            ImPlot::SetupLegend(ImPlotLocation_NorthEast);

                            if (opt_type == MD_VLX_OPT_CONSTRAINED) {
                                ImPlot::PlotLineSpline("Energy", x_vals, y_vals, (int)num_steps);
                            } else {
                                ImPlot::PlotLine("Energy", x_vals, y_vals, (int)num_steps);
                            }

                            if (opt.selected >= 0 && opt.selected < (int)num_steps) {
                                double val = x_vals[opt.selected];
                                ImPlot::PlotInfLines("##Selected", &val, 1);
                            }

                            plot_peaks("##opt_peaks", x_vals, y_vals, num_steps, opt.selected, opt.hovered, PlotPeaksFlags_Points);

                            if (opt_type == MD_VLX_OPT_IRC) {
                                // Plot a scatter point for the ts_index
                                const double x = x_vals[ts_index];
                                const double y = y_vals[ts_index];

                                ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 6.0f, ImVec4(0,0,0,-1), 2.0f);
                                ImPlot::SetNextFillStyle(ImVec4(0, 0, 0, -1), 0.5f);

                                ImPlot::PlotScatter("Transition State", &x, &y, 1);
                            }

                            if (ImPlot::IsPlotHovered() && opt.hovered != -1) {
                                ImGui::SetTooltip("Step %d\nEnergy: %.8f kJ/mol", opt.hovered + 1, y_vals[opt.hovered]);
                            }

                            ImPlot::EndPlot();
                        }
                        ImGui::Spacing();
                        int value = opt.selected + 1;
                        //ImGui::SetNextItemWidth(-FLT_MIN);
						const char* step_label = "Step";
						const ImVec2 label_size = ImGui::CalcTextSize(step_label);
						const float width = label_size.x + ImGui::GetStyle().ItemInnerSpacing.x;
                        ImGui::SetNextItemWidth(-width);
                        ImGui::SliderInt("Step", &value, 0, (int)num_steps);
                        opt.selected = value - 1;

                        ImGui::TreePop();
                    }
                    else {
                        opt.selected = (int)num_steps - 1;
                    }

                    static int prev_idx = -1;
                    if (opt.selected != prev_idx) {
                        prev_idx = opt.selected;
                        set_atom_coordinates(state, md_vlx_opt_coordinates(vlx, opt.selected));
                    }
                }
            }

            size_t num_atoms = md_vlx_number_of_atoms(vlx);
            if (num_atoms > 0 && ImGui::TreeNode("Geometry")) {
                static const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollX |
                                                        ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit;

                static const ImGuiTableColumnFlags columns_base_flags = ImGuiTableColumnFlags_NoSort;

                if (ImGui::BeginTable("Geometry Table", 5, flags, ImVec2(500, -1), 0)) {
                    const dvec3_t* atom_coord = md_vlx_opt_coordinates(vlx, opt.selected) ? md_vlx_opt_coordinates(vlx, opt.selected) : md_vlx_atom_coordinates(vlx);
                    const uint8_t* atom_nr    = md_vlx_atomic_numbers(vlx);

                    ImGui::TableSetupColumn("Atom", columns_base_flags, 0.0f);
                    ImGui::TableSetupColumn("Symbol", columns_base_flags, 0.0f);
                    ImGui::TableSetupColumn("Coord X", columns_base_flags, 0.0f);
                    ImGui::TableSetupColumn("Coord Y", columns_base_flags, 0.0f);
                    ImGui::TableSetupColumn("Coord Z", columns_base_flags | ImGuiTableColumnFlags_WidthFixed, 0.0f);
                    ImGui::TableSetupScrollFreeze(0, 1);
                    ImGui::TableHeadersRow();

                    ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);
                    ImGui::PushStyleColor(ImGuiCol_Header, IM_BLUE);
                    bool item_hovered = false;
                    for (size_t row_n = 0; row_n < num_atoms; row_n++) {

                        ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                        bool is_sel = md_bitfield_test_bit(&state.selection.selection_mask, row_n); //If atom is selected, mark it as such
                        bool is_hov = md_bitfield_test_bit(&state.selection.highlight_mask, row_n); //If atom is hovered,  mark it as such
                        ImGui::TableNextRow(ImGuiTableRowFlags_None, 0);
                        ImGui::TableNextColumn();

                        if (is_hov) {
                            ImGui::PushStyleColor(ImGuiCol_Header, IM_YELLOW);
                        }
                        else {
                            ImGui::PushStyleColor(ImGuiCol_Header, IM_BLUE);
                        }

                        char label[24];
                        snprintf(label, sizeof(label), "%zu", row_n + 1);
                        ImGui::Selectable(label, is_sel || is_hov, selectable_flags);
                        if (ImGui::TableGetHoveredRow() == (int)(row_n + 1)) {
                            if (state.mold.sys.atom.count > row_n) {
                                md_bitfield_clear(&state.selection.highlight_mask);
                                md_bitfield_set_bit(&state.selection.highlight_mask, row_n);
                                item_hovered = true;

                                //Selection
                                if (ImGui::IsKeyDown(ImGuiKey_MouseLeft) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                                    md_bitfield_set_bit(&state.selection.selection_mask, row_n);
                                }
                                //Deselect
                                else if (ImGui::IsKeyDown(ImGuiKey_MouseRight) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                                    md_bitfield_clear_bit(&state.selection.selection_mask, row_n);
                                }
                            }
                        }

                        ImGui::TableNextColumn();
                        str_t sym = md_util_element_symbol(atom_nr[row_n]);
                        ImGui::Text(STR_FMT, STR_ARG(sym));
                        ImGui::TableNextColumn();
                        ImGui::Text("%12.6f", atom_coord[row_n].x);
                        ImGui::TableNextColumn();
                        ImGui::Text("%12.6f", atom_coord[row_n].y);
                        ImGui::TableNextColumn();
                        ImGui::Text("%12.6f", atom_coord[row_n].z);

                        ImGui::PopStyleColor(1);
                                
                    }
                    if (!item_hovered && ImGui::IsWindowHovered()) {
                        //Makes sure that we clear the highlight if we are in this window, but don't hover an item
                        md_bitfield_clear(&state.selection.highlight_mask);
                    }

                    ImGui::PopStyleColor(2);
                    ImGui::EndTable();
                }
                ImGui::TreePop();
            }
            if (ImGui::TreeNode("Critical Points")) {
                ImGui::Checkbox("Enable Critical Points Analysis", &critical_points.enabled);

                static VolumeResolution vol_res = VolumeResolution::Mid;
                ImGui::Combo("Volume Resolution", (int*)&vol_res, volume_resolution_str, (int)VolumeResolution::Count);

                if (critical_points.enabled) {
                    static bool  reevaluate_graph_per_frame = false;
                    static bool  enable_simplification = true;
                    static float simp_threshold = 1.0E-4f;

                    ImGui::Checkbox("Extract graph per frame", &reevaluate_graph_per_frame);
                    if (reevaluate_graph_per_frame) {
                        critical_points.raw_hash = 0;
                    }

                    // Update volume if required
                    const double samples_per_unit_length = (float)(volume_resolution_samples_per_angstrom[(int)vol_res] * BOHR_TO_ANGSTROM);
                    md_grid_t grid = { 0 };
#if USE_OABB_GRID
                    init_grid(&grid, oabb.orientation, oabb.center, oabb.extents, samples_per_unit_length);
#else
                    init_grid(&grid, mat3_ident(), aabb.min_ext, aabb.max_ext, samples_per_unit_length);
#endif

                    // Compute a hash of the current state
                    size_t frame_idx = (size_t)(state.animation.frame + 0.5);
                    uint64_t raw_hash = md_hash64(&grid, sizeof(grid), frame_idx);

                    if (critical_points.raw_hash != raw_hash) {
                        critical_points.raw_hash = raw_hash;
                        md_topo_extremum_graph_free(&critical_points.raw_graph);

                        // The density comes out of the attribute table, evaluated by the same code
                        // path a representation uses. The only thing special here is the tail: the
                        // topology pass reads what the density pass wrote, on the same stream, so
                        // the volume is never read back to the host at all.
#if MD_ENABLE_GPU
                        // Ensure the topo context exists with the correct volume dimensions
                        uint32_t new_dims[3] = { (uint32_t)grid.dim[0], (uint32_t)grid.dim[1], (uint32_t)grid.dim[2] };
                        if (!critical_points.topo_ctx ||
                            critical_points.topo_dims[0] != new_dims[0] ||
                            critical_points.topo_dims[1] != new_dims[1] ||
                            critical_points.topo_dims[2] != new_dims[2])
                        {
                            md_topo_gpu_context_destroy(critical_points.topo_ctx);
                            critical_points.topo_ctx = md_topo_gpu_context_create(state.gpu_device, new_dims[0], new_dims[1], new_dims[2]);
                            critical_points.topo_dims[0] = new_dims[0];
                            critical_points.topo_dims[1] = new_dims[1];
                            critical_points.topo_dims[2] = new_dims[2];
                        }

                        if (!critical_points.topo_ctx) {
                            MD_LOG_ERROR("Failed to create topo GPU context for electron density");
                        } else if (density_evaluate_to_gpu_volume(&state, grid, es_path::alpha_density, nullptr, MD_GTO_OP_SET)) {
                            // Program order on the stream is the only dependency needed.
                            md_topo_gpu_record(state.gpu_stream, critical_points.topo_ctx, state.gpu_volume, &grid, 0.0f);
                            md_gpu_stream_sync(state.gpu_stream);
                            md_topo_gpu_context_extract(&critical_points.raw_graph, critical_points.topo_ctx);
                        } else {
                            MD_LOG_ERROR("Failed to evaluate the electron density for critical points analysis");
                        }
#else
                        uint32_t tex_id = 0;
                        if (gl::init_texture_3D(&tex_id, grid.dim[0], grid.dim[1], grid.dim[2], GL_R32F)) {
                            if (density_evaluate_gl(tex_id, grid, state.mold.sys, state.mold.state, es_path::alpha_density, nullptr, MD_GTO_OP_SET)) {
                                if (!md_topo_compute_extremum_graph_GPU(&critical_points.raw_graph, tex_id, &grid, 0.0f)) {
                                    MD_LOG_ERROR("Failed to compute extremum graph for electron density");
                                }
                            }
                            glDeleteTextures(1, &tex_id);
                        }
#endif
                    }

                    ImGui::Checkbox("Enable Graph Simplification", &enable_simplification);

                    uint64_t simp_hash = raw_hash;
                    if (enable_simplification) {
                        ImGui::SliderFloat("CP threshold", &simp_threshold, 0.0f, 1.0E-2f, "%.6f", ImGuiSliderFlags_Logarithmic);
                        simp_hash = md_hash64(&simp_threshold, sizeof(simp_threshold), raw_hash);
                    }

                    if (critical_points.simp_hash != simp_hash) {
                        critical_points.simp_hash = simp_hash;

                        md_topo_extremum_graph_free(&critical_points.simp_graph);
                        if (enable_simplification && critical_points.raw_graph.num_vertices > 0) {
                            md_topo_simplify(&critical_points.simp_graph, &critical_points.raw_graph, simp_threshold, true);
                        } else {
                            md_topo_extremum_graph_copy(&critical_points.simp_graph, &critical_points.raw_graph);
                        }
                    }

                    uint32_t vertex_count[MD_TOPO_NUM_TYPES] = { 0 };
                    md_topo_count_vertex_types(vertex_count, &critical_points.simp_graph);

                    const md_topo_extremum_graph_t& graph = critical_points.simp_graph;
                    ImGui::Text("Number of Edges:           %u", graph.num_edges);
                    ImGui::Text("Number of Critical Points: %u", graph.num_vertices);
                    ImGui::Text("\tNumber of Maxima:        %u", vertex_count[MD_TOPO_MAXIMUM]);
                    ImGui::Text("\tNumber of Split Saddles: %u", vertex_count[MD_TOPO_SPLIT_SADDLE]);
                    ImGui::Text("\tNumber of Minima:        %u", vertex_count[MD_TOPO_MINIMUM]);
                    ImGui::Text("\tNumber of Join Saddles:  %u", vertex_count[MD_TOPO_JOIN_SADDLE]);
                    if (ImGui::Button("Print Critical Points Info")) {
                        for (size_t i = 0; i < graph.num_vertices; ++i) {
                            const md_topo_vert_t& v = graph.vertices[i];
                            const md_topo_critical_point_type_t type = md_topo_vertex_type(&graph, i);
                            const char* type_cstr = md_topo_critical_point_type_str(type);
                            printf("Vertex %zu: Type=%s, Value=%.6f, Pos=(%.4f, %.4f, %.4f)\n", i, type_cstr, v.value, v.x, v.y, v.z);
                        }
                    }

                    if (graph.num_vertices > 0) {
                        enum CriticalPointColumn {
                            CriticalPointColumn_Idx,
                            CriticalPointColumn_X,
                            CriticalPointColumn_Y,
                            CriticalPointColumn_Z,
                            CriticalPointColumn_Type,
                            CriticalPointColumn_Value,
							CriticalPointColumn_Count
                        };

                        uint32_t* row_indices = md_temp_alloc_array(temp, uint32_t, graph.num_vertices);
                        for (uint32_t i = 0; i < graph.num_vertices; ++i) {
                            row_indices[i] = i;
                        }

                        static const ImGuiTableFlags table_flags = ImGuiTableFlags_RowBg |
                                                                   ImGuiTableFlags_Borders |
                                                                   ImGuiTableFlags_ScrollY |
                                                                   ImGuiTableFlags_Sortable |
                                                                   ImGuiTableFlags_SortMulti |
                                                                   ImGuiTableFlags_SizingStretchProp;

                        if (ImGui::IsWindowHovered()) {
                            md_bitfield_clear(&critical_points.highlight_mask);
                        }

                        if (ImGui::BeginTable("Critical Points Table", CriticalPointColumn_Count, table_flags)) {
							ImGui::TableSetupColumn("Idx",   ImGuiTableColumnFlags_DefaultSort, 0.0f, CriticalPointColumn_Idx);
                            ImGui::TableSetupColumn("X",     ImGuiTableColumnFlags_None,        0.0f, CriticalPointColumn_X);
                            ImGui::TableSetupColumn("Y",     ImGuiTableColumnFlags_None,        0.0f, CriticalPointColumn_Y);
                            ImGui::TableSetupColumn("Z",     ImGuiTableColumnFlags_None,        0.0f, CriticalPointColumn_Z);
                            ImGui::TableSetupColumn("Type",  ImGuiTableColumnFlags_None,        0.0f, CriticalPointColumn_Type);
                            ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_None,        0.0f, CriticalPointColumn_Value);
                            ImGui::TableSetupScrollFreeze(0, 1);
                            ImGui::TableHeadersRow();

                            if (ImGui::IsWindowHovered()) {
                                md_bitfield_clear(&critical_points.highlight_mask);
                            }

                            if (ImGuiTableSortSpecs* sort_specs = ImGui::TableGetSortSpecs()) {
                                if (sort_specs->SpecsCount > 0) {
                                    auto compare_rows = [&graph, sort_specs](uint32_t lhs_idx, uint32_t rhs_idx) {
                                        const md_topo_vert_t& lhs = graph.vertices[lhs_idx];
                                        const md_topo_vert_t& rhs = graph.vertices[rhs_idx];

                                        for (int spec_idx = 0; spec_idx < sort_specs->SpecsCount; ++spec_idx) {
                                            const ImGuiTableColumnSortSpecs& spec = sort_specs->Specs[spec_idx];
                                            int cmp = 0;
                                            switch ((CriticalPointColumn)spec.ColumnUserID) {
                                            case CriticalPointColumn_Idx:
                                                cmp = (lhs_idx < rhs_idx) ? -1 : (lhs_idx > rhs_idx ? 1 : 0);
												break;
                                            case CriticalPointColumn_X:
                                                cmp = (lhs.x < rhs.x) ? -1 : (lhs.x > rhs.x ? 1 : 0);
                                                break;
                                            case CriticalPointColumn_Y:
                                                cmp = (lhs.y < rhs.y) ? -1 : (lhs.y > rhs.y ? 1 : 0);
                                                break;
                                            case CriticalPointColumn_Z:
                                                cmp = (lhs.z < rhs.z) ? -1 : (lhs.z > rhs.z ? 1 : 0);
                                                break;
                                            case CriticalPointColumn_Type:
                                                cmp = (graph.types[lhs_idx] < graph.types[rhs_idx]) ? -1 : (graph.types[lhs_idx] > graph.types[rhs_idx] ? 1 : 0);
                                                break;
                                            case CriticalPointColumn_Value:
                                                cmp = (lhs.value < rhs.value) ? -1 : (lhs.value > rhs.value ? 1 : 0);
                                                break;
                                            default:
                                                break;
                                            }

                                            if (cmp != 0) {
                                                return spec.SortDirection == ImGuiSortDirection_Ascending ? cmp < 0 : cmp > 0;
                                            }
                                        }

                                        return lhs_idx < rhs_idx;
                                    };

                                    std::sort(row_indices, row_indices + graph.num_vertices, compare_rows);
                                    sort_specs->SpecsDirty = false;
                                }
                            }

                            ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);
                            ImGui::PushStyleColor(ImGuiCol_Header, IM_BLUE);

                            for (size_t row = 0; row < graph.num_vertices; ++row) {
                                const uint32_t idx = row_indices[row];
                                const md_topo_vert_t& v = graph.vertices[idx];
                                const int type = graph.types[idx];

                                ImGui::TableNextRow();
                                ImGui::TableNextColumn();

                                ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                                bool is_sel = md_bitfield_test_bit(&critical_points.selection_mask, idx);
                                bool is_hov = md_bitfield_test_bit(&critical_points.highlight_mask, idx);

                                if (is_hov) {
                                    ImGui::PushStyleColor(ImGuiCol_Header, IM_YELLOW);
                                }
                                else {
                                    ImGui::PushStyleColor(ImGuiCol_Header, IM_BLUE);
                                }

                                char lbl[16];
                                snprintf(lbl, sizeof(lbl), "%u", idx);
                                ImGui::Selectable(lbl, is_sel || is_hov, selectable_flags);

                                if (ImGui::TableGetHoveredRow() == (int)(row + 1)) {
                                    if (idx < graph.num_vertices) {
                                        md_bitfield_clear(&critical_points.highlight_mask);
                                        md_bitfield_set_bit(&critical_points.highlight_mask, idx);

                                        //Selection
                                        if (ImGui::IsKeyDown(ImGuiMod_Shift)) {
                                            if (ImGui::IsKeyDown(ImGuiKey_MouseLeft)) {
                                                md_bitfield_set_bit(&critical_points.selection_mask, idx);
                                            } else if (ImGui::IsKeyDown(ImGuiKey_MouseRight)) {
                                                md_bitfield_clear_bit(&critical_points.selection_mask, idx);
                                            }
                                        }
                                    }
                                }
                                ImGui::PopStyleColor(1);
                                ImGui::TableNextColumn();
                                ImGui::Text("%.6f", v.x);
                                ImGui::TableNextColumn();
                                ImGui::Text("%.6f", v.y);
                                ImGui::TableNextColumn();
                                ImGui::Text("%.6f", v.z);
                                ImGui::TableNextColumn();
                                ImGui::TextUnformatted(md_topo_critical_point_type_str(type));
                                ImGui::TableNextColumn();
                                ImGui::Text("%.8f", v.value);                                
                            }

                            ImGui::PopStyleColor(2);

                            ImGui::EndTable();
                        }
                    }


                }

                ImGui::TreePop();
            }
        }
        ImGui::End();
    }

    typedef struct {
        double har_freq;
        double redu_mass;
        double force_const;
        double ir_intens;
        double* x;
        double* y;
        double* z;
    } vibration_mode;

    void draw_rsp_spectra_export_window(ApplicationState& state) {
#if 0
        (void)state;
        ASSERT(&state);

        struct ExportFormat {
            str_t lbl;
            str_t ext;
        };

        ExportFormat table_formats[]{
            {STR_LIT("XVG"), STR_LIT("xvg")},
            {STR_LIT("CSV"), STR_LIT("csv")}
        };

        struct ExportProperty {
            double* x = 0;
            double* y = 0;
            str_t lable;
            str_t y_unit;
        };

        ExportProperty properties[]{
            {rsp.x_spectra, rsp.y_spectra_abs,   STR_LIT("Absorption"), STR_LIT((const char*)u8"ε (L mol⁻¹ cm⁻¹)")},
            {rsp.x_spectra, rsp.y_spectra_ecd,   STR_LIT("ECD"),        STR_LIT((const char*)u8"Δε(ω) (L mol⁻¹ cm⁻¹)")},
            {rsp.vib_x,          rsp.vib_y, STR_LIT("Vibration"),  STR_LIT("IR Intensity (km/mol)")}
        };
        
        if (ImGui::Begin("Spectra Export", &rsp.show_export_window)) {
            static int table_format = 0;
            static int property_idx = 0;
            const char* x_unit = "DEBUG";
            
            //TODO: Add sanity checks
            
            ImGui::PushItemWidth(200);

            str_t file_extension = {};
            if (ImGui::BeginCombo("File Format", table_formats[table_format].lbl.ptr)) {
                for (int i = 0; i < (int)ARRAY_SIZE(table_formats); ++i) {
                    if (ImGui::Selectable(table_formats[i].lbl.ptr, table_format == i)) {
                        table_format = i;
                    }
                }
                ImGui::EndCombo();
            }
            file_extension = table_formats[table_format].ext;

            if (ImGui::BeginCombo("Property", properties[property_idx].lable.ptr)) {
                for (int i = 0; i < (int)ARRAY_SIZE(properties); ++i) {
                    if (ImGui::Selectable(properties[i].lable.ptr, property_idx == i)) {
                        property_idx = i;
                    }
                }
                ImGui::EndCombo();
            }

            if (property_idx == 0 || property_idx == 1) {
                x_unit = x_unit_full_str[rsp.x_unit];
            }
            else if (property_idx == 2) {
                x_unit = (const char*)u8"Harmonic Frequency (cm⁻¹)";
            }

            static bool export_valid = true;
            bool export_clicked = ImGui::Button("Export");
            if (export_clicked) {
                export_valid = true;

                if (export_valid) {
                    char path_buf[1024];
                    md_array(const float*)  column_data = 0;
                    md_array(str_t)   column_labels = 0;
                    md_array(str_t)         legends = 0;

                    md_array(float) x_values = md_array_create(float, 1024, arena);
                    md_array(float) y_values = md_array_create(float, 1024, arena);

                    for (size_t i = 0; i < 1024; i++) {
                        x_values[i] = (float)properties[property_idx].x[i];
                        y_values[i] = (float)properties[property_idx].y[i];
                    }

                    //str_t x_label = str_from_cstr(x_unit);
                    md_array_push(column_labels, str_from_cstr(x_unit), arena);
                    //str_t y_label = str_from_cstr("Y");
                    md_array_push(column_labels, properties[property_idx].y_unit, arena);

                    md_array_push(column_data, x_values, arena);
                    md_array_push(column_data, y_values, arena);

                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, file_extension)) {
                        str_t path = { path_buf, strnlen(path_buf, sizeof(path_buf)) };
                        if (table_format == 0) {
                            //md_xvg_write
                            //TODO: Implement md_xvg_write_to_file
                            str_t header = md_xvg_format_header(properties[property_idx].lable, str_from_cstr(x_unit), properties[property_idx].y_unit, 0, legends, arena);
                            str_t xvg = md_xvg_format(header, 2, 1024, column_data, arena);
                            md_file_t file = { 0 };
                            if (md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                                const size_t written_bytes = md_file_write(file, xvg.ptr, xvg.len);
                                if (written_bytes != xvg.len) {
                                    MD_LOG_ERROR("CSV: Unexpected error, some bytes were not written");
                                }
                                md_file_close(&file);
                            }
                            else {
                                MD_LOG_ERROR("CSV: File could not be opened for writing: '%.*s'", (int)path.len, path.ptr);
                            }

                        }
                        else if (table_format == 1) {
                            //export_csv(column_data, column_labels, 0, 0, path);
                            md_csv_write_to_file(column_data, column_labels, 2, 1024, path);
                        }
                    }
                }
            }
            ImGui::PopItemWidth();
            if (!export_valid) { ImGui::TextWrapped("Values are not valid, make sure that you have opened plot windows once to trigger calculations"); }
        }
        ImGui::End();
#endif
    }

    template<typename T>
    static inline void find_min_max(int* min_idx, int* max_idx, const T* arr, size_t n) {
        if (min_idx) *min_idx = 0;
        if (max_idx) *max_idx = 0;
        for (size_t i = 1; i < n; i++) {
            if (min_idx) {
                if (arr[i] < arr[*min_idx]) {
                    *min_idx = (int)i;
                }
            }
            if (max_idx) {
                if (arr[i] > arr[*max_idx]) {
                    *max_idx = (int)i;
                }
            }
        }
    }

    void draw_rsp_window(ApplicationState& state) {
        md_temp_scope_t temp = md_temp_begin_in(state.allocator.frame);
        defer { md_temp_end(temp); };

        if (!rsp.show_window) return;
		md_vlx_rsp_type_t rsp_type = md_vlx_rsp_type(vlx);
        const size_t num_frequencies  = md_vlx_rsp_number_of_frequencies(vlx);
        const size_t num_normal_modes = md_vlx_vib_number_of_normal_modes(vlx);
        // If the spetrum is broadened from peaks, this is the number of samples in the visual space (x_min -> x_max) used to sample that signal.
        constexpr int num_broadened_samples = 2048;

        // XPS is a delta-SCF property and can be the only thing in the file, so it has to be able to
        // keep this window alive on its own.
        const bool has_xps = md_vlx_has_xps(vlx);

        if (num_frequencies == 0 && num_normal_modes == 0 && !has_xps) return;

        // The frequencies are given in a.u. but can either be peaks to broaden if LINEAR, or spectrum samples if CPP
        const double* x_freq_au   = md_vlx_rsp_frequencies(vlx);

        double x_freq_au_min = DBL_MAX;
        double x_freq_au_max = -DBL_MAX;
        for (size_t i = 0; i < num_frequencies; i++) {
            x_freq_au_min = MIN(x_freq_au_min, x_freq_au[i]);
            x_freq_au_max = MAX(x_freq_au_max, x_freq_au[i]);
        }

        ImPlot::PushStyleVar(ImPlotStyleVar_FitPadding, ImVec2(0.1f, 0.1f));
        defer { ImPlot::PopStyleVar(); };

        ImGui::SetNextWindowSize({ 300, 350 }, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Response", &rsp.show_window, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("File")) {
                    if (ImGui::MenuItem("Export")) {
                        rsp.show_export_window = true;
                    }
                    ImGui::EndMenu();
                }
                if (num_normal_modes > 0) {
                    if (ImGui::BeginMenu("Settings")) {
                        ImGui::Checkbox("VIB: Displace Atomic Orbitals", &vib.displace_aos);
                        ImGui::EndMenu();
                    }
                }
                ImGui::EndMenuBar();
            }

            struct SpectrumGetterData {
                const double* x_peaks;
                const double* y_peaks;
                size_t num_peaks;
                double x_min;
                double x_max;
                double broadening_fwhm;
                size_t num_samples;
                x_unit_t x_unit = X_UNIT_EV;
            };

            const ImGuiTreeNodeFlags tree_flags = ImGuiTreeNodeFlags_DefaultOpen | ImGuiTreeNodeFlags_NoTreePushOnOpen;

            // RIXS is its own thing: the 2D map has no meaningful representation in the shared
            // 1D spectrum plots below, so it gets its own section.
            if (rsp_type == MD_VLX_RSP_RIXS) {
                RixsMapInput rixs_input = {};
                if (rixs_map_input_from_vlx(rixs_input)) {
                    ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
                    if (ImGui::TreeNodeEx("RIXS", tree_flags)) {
                        draw_rixs_section(rixs_input, rixs);
                    }
                }
            }

            // XPS is independent of the response type, so it sits outside the electronic spectroscopy
            // section below (which is driven by the rsp eigenvalues).
            if (has_xps) {
                ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
                if (ImGui::TreeNodeEx("XPS", tree_flags)) {
                    draw_xps_plot(state);
                }
            }

            // Hide the ECD plot if the RSP is a VIB calculation, since ECD is not relevant for VIB, but can be involved as an intermediate result holder.
            if (num_normal_modes == 0 && num_frequencies > 0 && ImGui::TreeNode("Electronic Spectroscopy")) {
                // Explicit x axis limits for linking
                static double x_min = 0.0;
                static double x_max = 10.0;

                const float avail_width = ImGui::GetContentRegionAvail().x;
                ImGui::PushItemWidth(MIN(avail_width, 200));
                // Broadening min max values
                static const double broadening_min = 0.05;
                static const double broadening_max = 1.0;

                ImGui::SliderScalar((const char*)u8"Broadening FWHM (eV)", ImGuiDataType_Double, &rsp.broadening_fwhm, &broadening_min, &broadening_max);
                ImGui::Combo("Broadening mode", (int*)(&rsp.broadening_mode), broadening_mode_str, BROADENING_MODE_COUNT);
                ImGui::Combo("X unit", (int*)(&rsp.x_unit), x_unit_full_str, X_UNIT_COUNT);
                ImGui::PopItemWidth();

                // Samples are only used in CPP RSP, these hold intermediate results, i.e. converted to the selected x unit.
                size_t num_samples = 0;
                double* x_samples = NULL;
                const double* y_samples_sigma = md_vlx_rsp_sigma(vlx);
                const double* y_samples_delta_epsilons = md_vlx_rsp_delta_epsilons(vlx);
                const double* y_samples_ord = md_vlx_rsp_optical_rotations(vlx);
                const double* y_samples_cs = md_vlx_rsp_tpa_cross_sections(vlx);

                // Peaks are only used in LINEAR RSP, not in CPP
				size_t num_peaks = 0;
				double* x_peaks = NULL;
                const double* y_peaks_osc = md_vlx_rsp_oscillator_strengths(vlx);
                const double* y_peaks_cgs = md_vlx_rsp_rotatory_strengths(vlx);
                const double* y_peaks_tpa_trans_linear = md_vlx_rsp_tpa_trans_linear(vlx);
                const double* y_peaks_tpa_trans_circular = md_vlx_rsp_tpa_trans_circular(vlx);

                bool refit = false;

                // Construct a hash to compare against if we need to refit.
                uint64_t hash = md_hash64(&rsp.broadening_fwhm, sizeof(rsp.broadening_fwhm), rsp.broadening_mode);
                hash = md_hash64(x_freq_au, num_frequencies * sizeof(double), hash);
                hash = md_hash64(&rsp.x_unit, sizeof(rsp.x_unit), hash);

                if (hash != rsp.hash) {
                    rsp.hash = hash;
                    refit = true;
                }

				if (rsp_type == MD_VLX_RSP_LINEAR || rsp_type == MD_VLX_RSP_TPA_TRANSITION) {
					num_peaks = num_frequencies;
                    x_peaks = md_temp_alloc_array(temp, double, num_peaks);
                    convert_values_from_au(x_peaks, x_freq_au, num_peaks, rsp.x_unit);
				} else {
                    num_samples = num_frequencies;
                    x_samples = md_temp_alloc_array(temp, double, num_samples);
                    convert_values_from_au(x_samples, x_freq_au, num_samples, rsp.x_unit);
                }

                if ((num_samples > 0 && y_samples_sigma) || (num_peaks > 0 && y_peaks_osc)) {
                    ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
                    if (ImGui::TreeNodeEx("Absorption Plot", tree_flags)) {
                        ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

                        if (refit) {
                            ImPlot::SetNextAxesToFit();
                        }

                        auto getter_lorentzian = [](int idx, void* user_data) -> ImPlotPoint {
                            const SpectrumGetterData* data = (const SpectrumGetterData*)user_data;

                            // Construct an x value for idx as a linear space between x_min and x_max.
                            double x = data->x_min + (data->x_max - data->x_min) * idx / (data->num_samples - 1);

                            // Convert x to a.u. for the broadening function
                            double x_au = convert_value_to_au(x, data->x_unit);

                            // Compute the corresponding y value by broadening the peaks at this x value.
                            double y = lorentzian_abs(x_au, data->x_peaks, data->y_peaks, data->num_peaks, data->broadening_fwhm);
                            double eps = x_au * y * OSCILLATOR_STRENGTH_TO_EPSILON;

                            return ImPlotPoint{x, eps};
                        };

                        auto getter_gaussian = [](int idx, void* user_data) -> ImPlotPoint {
                            const SpectrumGetterData* data = (const SpectrumGetterData*)user_data;

                            // Construct an x value for idx as a linear space between x_min and x_max.
                            double x = data->x_min + (data->x_max - data->x_min) * idx / (data->num_samples - 1);

                            // Convert x to a.u. for the broadening function
                            double x_au = convert_value_to_au(x, data->x_unit);

                            // Compute the corresponding y value by broadening the peaks at this x value.
                            double y = gaussian_abs(x_au, data->x_peaks, data->y_peaks, data->num_peaks, data->broadening_fwhm);
                            double eps = x_au * y * OSCILLATOR_STRENGTH_TO_EPSILON;

                            return ImPlotPoint{ x, eps };
                        };

                        ImPlotGetter getter = (rsp.broadening_mode == BROADENING_MODE_LORENTZIAN) ? getter_lorentzian : getter_gaussian;

                        double* y_samples = md_temp_alloc_zero_array(temp, double, num_samples);

                        if (ImPlot::BeginPlot("Absorption")) {
                            ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                            ImPlot::SetupAxis(ImAxis_X1, x_unit_full_str[rsp.x_unit]);
                            ImPlot::SetupAxis(ImAxis_Y1, (const char*)u8"ε (L mol⁻¹ cm⁻¹)");
                            if (num_peaks > 0) {
                                ImPlot::SetupAxis(ImAxis_Y2, "f", ImPlotAxisFlags_AuxDefault);
                            }
                            ImPlot::SetupFinish();

                            ImPlot::SetAxis(ImAxis_Y1);
                            if (rsp_type == MD_VLX_RSP_CPP) {
                                // Convert y_samples to eps
                                constexpr double sigma_to_epsilon = 7323.816924863764;
                                for (size_t i = 0; i < num_samples; ++i) {
                                    y_samples[i] = y_samples_sigma[i] * sigma_to_epsilon;
                                }
                                ImPlot::PlotLineSpline("Interpolated Spectrum", x_samples, y_samples, (int)num_samples);
                                int selected = -1;
                                int hovered = -1;
                                plot_peaks("Samples", x_samples, y_samples, num_samples, selected, hovered, PlotPeaksFlags_Points);

                                if (ImPlot::IsPlotHovered() && hovered != -1) {
                                    double x_val = x_samples[hovered];
                                    double y_val = y_samples[hovered];
                                    const char* x_unit = x_unit_short_str[rsp.x_unit];
                                    const char* y_unit = (const char*)u8"ε (L mol⁻¹ cm⁻¹)";
                                    if (ImGui::BeginTooltip()) {
                                        ImGui::Text("%f %s, %f %s", x_val, x_unit, y_val, y_unit);
                                        ImGui::EndTooltip();
                                    }
                                }
                            } else {
                                SpectrumGetterData data = { x_freq_au, y_peaks_osc, num_peaks, x_min, x_max, rsp.broadening_fwhm * EV_TO_HARTREE, num_broadened_samples, rsp.x_unit};
                                ImPlot::PlotLineG("Broadened Spectrum", getter, &data, num_broadened_samples);
                            }

                            if (num_peaks > 0 && x_peaks && y_peaks_osc) {
                                ImPlot::SetAxis(ImAxis_Y2);
                                plot_peaks("Oscillator Strength", x_peaks, y_peaks_osc, num_peaks, rsp.selected, rsp.hovered);

                                if (ImPlot::IsPlotHovered() && rsp.hovered != -1) {
                                    double x_val = x_peaks[rsp.hovered];
                                    double y_val = y_peaks_osc[rsp.hovered];
                                    int state_idx = rsp.hovered + 1;
                                    const char* x_label = x_unit_label_str[rsp.x_unit];
                                    const char* x_unit = x_unit_short_str[rsp.x_unit];
                                    if (ImGui::BeginTooltip()) {
                                        ImGui::Text("State %i: %s = %f (%s), f = %f", state_idx, x_label, x_val, x_unit, y_val);
                                        ImGui::EndTooltip();
                                    }
                                }
                            }

                            if (ImPlot::FitThisFrame()) {
                                ImPlotPlot* plot = ImPlot::GetCurrentPlot();
                                if (plot) {
                                    const double min_x = convert_value_from_au(x_freq_au_min, rsp.x_unit);
                                    const double max_x = convert_value_from_au(x_freq_au_max, rsp.x_unit);
                                    plot->Axes[ImAxis_X1].FitExtents = ImPlotRange(MIN(min_x, max_x), MAX(min_x, max_x));

                                    // FIND MAX Y AXES
                                    if (rsp_type == MD_VLX_RSP_CPP) {
                                        int max_idx;
                                        find_min_max(NULL, &max_idx, y_samples, num_samples);
                                        plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(0.0, y_samples[max_idx]);
                                    } else {
                                        // Find the maximum peak and calculate the maximum value over the spectrum.
                                        int max_peak_idx = 0;
                                        find_min_max(NULL, &max_peak_idx, y_peaks_osc, num_peaks);
                                        double max_y = y_peaks_osc[max_peak_idx];

                                        double max_eps = 0.0;
                                        SpectrumGetterData data = { x_freq_au, y_peaks_osc, num_peaks, min_x, max_x, rsp.broadening_fwhm * EV_TO_HARTREE, num_broadened_samples, rsp.x_unit};
                                        for (int i = 0; i < num_broadened_samples; ++i) {
                                            ImPlotPoint point = getter(i, &data);
                                            max_eps = MAX(max_eps, point.y);
                                        }
                                        
                                        plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(0.0, max_eps);
                                        plot->Axes[ImAxis_Y2].FitExtents = ImPlotRange(0.0, max_y);
                                    }
                                }
                            }
                            ImPlot::EndPlot();
                        }
                    }
                }
                
                if ((num_samples > 0 && y_samples_delta_epsilons) || (num_peaks > 0 && y_peaks_cgs)) {
                    ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
                    if (ImGui::TreeNodeEx("ECD Plot", tree_flags)) {
                        ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

                        if (refit) {
                            ImPlot::SetNextAxesToFit();
                        }

                        auto getter_lorentzian = [](int idx, void* user_data) -> ImPlotPoint {
                            const SpectrumGetterData* data = (const SpectrumGetterData*)user_data;

                            // Construct an x value for idx as a linear space between x_min and x_max.
                            double x = data->x_min + (data->x_max - data->x_min) * idx / (data->num_samples - 1);

                            // Convert x to a.u. for the broadening function
                            double x_au = convert_value_to_au(x, data->x_unit);

                            // Compute the corresponding y value by broadening the peaks at this x value.
                            double y = lorentzian_ecd(x_au, data->x_peaks, data->y_peaks, data->num_peaks, data->broadening_fwhm);
                            double delta_eps = x_au * y * ROTATORY_STRENGTH_TO_DELTA_EPSILON;

                            return ImPlotPoint{x, delta_eps};
                        };

                        auto getter_gaussian = [](int idx, void* user_data) -> ImPlotPoint {
                            const SpectrumGetterData* data = (const SpectrumGetterData*)user_data;

                            // Construct an x value for idx as a linear space between x_min and x_max.
                            double x = data->x_min + (data->x_max - data->x_min) * idx / (data->num_samples - 1);

                            // Convert x to a.u. for the broadening function
                            double x_au = convert_value_to_au(x, data->x_unit);

                            // Compute the corresponding y value by broadening the peaks at this x value.
                            double y = gaussian_ecd(x_au, data->x_peaks, data->y_peaks, data->num_peaks, data->broadening_fwhm);
                            double delta_eps = x_au * y * ROTATORY_STRENGTH_TO_DELTA_EPSILON;

                            return ImPlotPoint{ x, delta_eps };
                        };

                        ImPlotGetter getter = (rsp.broadening_mode == BROADENING_MODE_LORENTZIAN) ? getter_lorentzian : getter_gaussian;

                        if (ImPlot::BeginPlot("ECD")) {
                            ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                            ImPlot::SetupAxis(ImAxis_X1, x_unit_full_str[rsp.x_unit]);
                            ImPlot::SetupAxis(ImAxis_Y1, (const char*)u8"Δε(ω) (L mol⁻¹ cm⁻¹)");
                            if (num_peaks > 0) {
                                ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"R (10⁻⁴⁰ cgs)", ImPlotAxisFlags_AuxDefault);
                            }
                            ImPlot::SetupFinish();

                            ImPlot::SetAxis(ImAxis_Y1);
                            if (rsp_type == MD_VLX_RSP_CPP) {
                                ImPlot::PlotLineSpline("Interpolated Spectrum", x_samples, y_samples_delta_epsilons, (int)num_samples);

                                int selected = -1;
                                int hovered = -1;
                                plot_peaks("Samples", x_samples, y_samples_delta_epsilons, num_samples, selected, hovered, PlotPeaksFlags_Points);

                                if (ImPlot::IsPlotHovered() && hovered != -1) {
                                    double x_val = x_samples[hovered];
                                    double y_val = y_samples_delta_epsilons[hovered];
                                    const char* x_unit = x_unit_short_str[rsp.x_unit];
                                    const char* y_unit = (const char*)u8"Δε(ω) (L mol⁻¹ cm⁻¹)";
                                    if (ImGui::BeginTooltip()) {
                                        ImGui::Text("%f %s, %f %s", x_val, x_unit, y_val, y_unit);
                                        ImGui::EndTooltip();
                                    }
                                }
                            } else {
                                SpectrumGetterData data = { x_freq_au, y_peaks_cgs, num_peaks, x_min, x_max, rsp.broadening_fwhm * EV_TO_HARTREE, num_broadened_samples, rsp.x_unit};
                                ImPlot::PlotLineG("Broadened Spectrum", getter, &data, num_broadened_samples);
                            }

                            if (num_peaks > 0 && x_peaks && y_peaks_cgs) {
                                ImPlot::SetAxis(ImAxis_Y2);
                                plot_peaks("Rotatory Strength", x_peaks, y_peaks_cgs, num_peaks, rsp.selected, rsp.hovered);

                                if (ImPlot::IsPlotHovered() && rsp.hovered != -1) {
                                    double x_val = x_peaks[rsp.hovered];
                                    double y_val = y_peaks_cgs[rsp.hovered];
                                    int state_idx = rsp.hovered + 1;
                                    const char* x_label = x_unit_label_str[rsp.x_unit];
                                    const char* x_unit  = x_unit_short_str[rsp.x_unit];
                                    const char* y_unit  = (const char*)u8"10⁻⁴⁰ cgs";
                                    if (ImGui::BeginTooltip()) {
                                        ImGui::Text((const char*)u8"State %i: %s = %f (%s), R = %f (%s)", state_idx, x_label, x_val, x_unit, y_val, y_unit);
                                        ImGui::EndTooltip();
                                    }
                                }
                            }

                            if (ImPlot::FitThisFrame()) {
                                ImPlotPlot* plot = ImPlot::GetCurrentPlot();
                                if (plot) {
                                    const double min_x = convert_value_from_au(x_freq_au_min, rsp.x_unit);
                                    const double max_x = convert_value_from_au(x_freq_au_max, rsp.x_unit);
                                    plot->Axes[ImAxis_X1].FitExtents = ImPlotRange(MIN(min_x, max_x), MAX(min_x, max_x));

                                    // FIND MAX Y AXES
                                    if (rsp_type == MD_VLX_RSP_CPP) {
                                        int min_idx, max_idx;
                                        find_min_max(&min_idx, &max_idx, y_samples_delta_epsilons, num_samples);
                                        double max_delta_eps = MAX(fabs(y_samples_delta_epsilons[min_idx]), fabs(y_samples_delta_epsilons[max_idx]));
                                        plot->Axes[ImAxis_Y2].FitExtents = ImPlotRange(-max_delta_eps, max_delta_eps);
                                    } else {
                                        // Find the maximum peak and calculate the maximum value over the spectrum.
                                        int min_peak_idx, max_peak_idx;
                                        find_min_max(&min_peak_idx, &max_peak_idx, y_peaks_cgs, num_peaks);
                                        int peak_idx = (fabs(y_peaks_cgs[min_peak_idx]) > fabs(y_peaks_cgs[max_peak_idx])) ? min_peak_idx : max_peak_idx;
                                        double max_y = fabs(y_peaks_cgs[peak_idx]);

                                        double max_delta_eps = 0.0;
                                        SpectrumGetterData data = { x_freq_au, y_peaks_cgs, num_peaks, min_x, max_x, rsp.broadening_fwhm * EV_TO_HARTREE, num_broadened_samples, rsp.x_unit};
                                        for (int i = 0; i < num_broadened_samples; ++i) {
                                            ImPlotPoint point = getter(i, &data);
                                            max_delta_eps = MAX(max_delta_eps, fabs(point.y));
                                        }
                                        
                                        plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(-max_delta_eps, max_delta_eps);
                                        plot->Axes[ImAxis_Y2].FitExtents = ImPlotRange(-max_y, max_y);
                                    }
                                }
                            }
                            ImPlot::EndPlot();
                        }
                    }
                }

                if (num_samples > 0 && y_samples_ord) {
                    ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);

                    if (ImGui::TreeNodeEx("ORD Plot", tree_flags)) {
                        ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

                        if (refit) {
                            ImPlot::SetNextAxesToFit();
                        }

                        if (ImPlot::BeginPlot("Optical Rotatory Dispersion")) {
                            const char* y_unit = (const char*)u8"ORD (10³ ° cm^3 g⁻¹ dm⁻¹)";
                            ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                            ImPlot::SetupAxis(ImAxis_X1, x_unit_full_str[rsp.x_unit]);
                            ImPlot::SetupAxis(ImAxis_Y1, y_unit);
                            ImPlot::SetupFinish();

                            ImPlot::SetAxis(ImAxis_Y1);
                            ImPlot::PlotLineSpline("Interpolated Spectrum", x_samples, y_samples_ord, (int)num_samples);
                            int selected = -1;
                            int hovered = -1;
                            plot_peaks("Samples", x_samples, y_samples_ord, num_samples, selected, hovered, PlotPeaksFlags_Points);

                            if (ImPlot::IsPlotHovered() && hovered != -1) {
                                double x_val = x_samples[hovered];
                                double y_val = y_samples_ord[hovered];
                                const char* x_unit = x_unit_short_str[rsp.x_unit];                                    
                                if (ImGui::BeginTooltip()) {
                                    ImGui::Text("%f %s, %f %s", x_val, x_unit, y_val, y_unit);
                                    ImGui::EndTooltip();
                                }
                            }

                            if (ImPlot::FitThisFrame()) {
                                ImPlotPlot* plot = ImPlot::GetCurrentPlot();
                                if (plot) {
                                    const double min_x = convert_value_from_au(x_freq_au_min, rsp.x_unit);
                                    const double max_x = convert_value_from_au(x_freq_au_max, rsp.x_unit);
                                    plot->Axes[ImAxis_X1].FitExtents = ImPlotRange(MIN(min_x, max_x), MAX(min_x, max_x));

                                    int min_idx, max_idx;
                                    find_min_max(&min_idx, &max_idx, y_samples_ord, num_samples);
                                    double max_ord = MAX(fabs(y_samples_ord[min_idx]), fabs(y_samples_ord[max_idx]));
                                    plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(-max_ord, max_ord);
                                }
                            }
                            ImPlot::EndPlot();
                        }
                    }
                }

                if ((num_samples > 0 && y_samples_cs) || (num_peaks > 0 && y_peaks_tpa_trans_linear)) {
                    ImGui::SetNextItemOpen(true, ImGuiCond_Once);
                    if (ImGui::TreeNodeEx("Two Photon Absorption", tree_flags)) {

                        char x_label[64];
                        snprintf(x_label, sizeof(x_label), "Photon Energy [%s]", x_unit_short_str[rsp.x_unit]);
                        const char* title = (rsp_type == MD_VLX_RSP_TPA_TRANSITION) ? "TPA Transition" : "TPA";

                        auto getter = [](int idx, void* user_data) -> ImPlotPoint {
                            const SpectrumGetterData* data = (const SpectrumGetterData*)user_data;

                            // Construct an x value for idx as a linear space between x_min and x_max.
                            double x = data->x_min + (data->x_max - data->x_min) * idx / (data->num_samples - 1);

                            const double x_au = convert_value_to_au(x, data->x_unit);

                            // Compute the corresponding y value by broadening the peaks at this x value.
                            // Perform broadening in a.u.
                            // Lorentzian broadening
                            const double factor = AU_TO_GM * x_au * x_au;
                            double y = factor * lorentzian_tpa(x_au, data->x_peaks, data->y_peaks, data->num_peaks, data->broadening_fwhm);

                            return ImPlotPoint{x, y};
                        };

                        ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

                        if (refit) {
                            ImPlot::SetNextAxesToFit();
                        }

                        if (ImPlot::BeginPlot(title)) {
                            ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                            ImPlot::SetupAxis(ImAxis_X1, x_label);
                            ImPlot::SetupAxis(ImAxis_Y1, (const char*)u8"TPA Cross-Section [GM]");
                            if (rsp_type == MD_VLX_RSP_TPA_TRANSITION) {
                                ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"TPA Strengths (linear) [a.u.]", ImPlotAxisFlags_AuxDefault);
                            }
                            ImPlot::SetupFinish();

                            ImPlot::SetAxis(ImAxis_Y1);
                            if (rsp_type == MD_VLX_RSP_TPA_TRANSITION) {
                                SpectrumGetterData getter_data = { x_freq_au, y_peaks_tpa_trans_linear, num_peaks, x_min, x_max, rsp.broadening_fwhm * EV_TO_HARTREE, num_broadened_samples, rsp.x_unit};
                                ImPlot::PlotLineG("Broadened Spectrum", getter, &getter_data, num_broadened_samples);

                                ImPlot::SetAxis(ImAxis_Y2);
                                plot_peaks("Linear TPA Strength", x_peaks, y_peaks_tpa_trans_linear, num_peaks, rsp.selected, rsp.hovered);
                            } else {
                                ImPlot::PlotLineSpline("Interpolated Spectrum", x_samples, y_samples_cs, (int)num_samples);
                                plot_peaks("Computed Cross-Sections", x_samples, y_samples_cs, num_samples, rsp.selected, rsp.hovered, PlotPeaksFlags_Points);
                            }

                            if (ImPlot::FitThisFrame() || refit) {
                                ImPlotPlot* plot = ImPlot::GetCurrentPlot();
                                if (plot) {
                                    const double min_x = convert_value_from_au(x_freq_au_min, rsp.x_unit);
                                    const double max_x = convert_value_from_au(x_freq_au_max, rsp.x_unit);
                                    plot->Axes[ImAxis_X1].FitExtents = ImPlotRange(MIN(min_x, max_x), MAX(min_x, max_x));

                                    // FIND MAX Y AXES
                                    if (rsp_type == MD_VLX_RSP_TPA) {
                                        int max_idx;
                                        find_min_max(NULL, &max_idx, y_samples_cs, num_samples);
                                        double max_y = y_samples_cs[max_idx];
                                        plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(0.0, max_y);
                                    } else {
                                        int max_idx;
                                        find_min_max(NULL, &max_idx, y_peaks_tpa_trans_linear, num_peaks);
                                        double max_y = y_peaks_tpa_trans_linear[max_idx];
                                        double max_cs = 0.0;
                                        SpectrumGetterData data = { x_freq_au, y_peaks_tpa_trans_linear, num_peaks, min_x, max_x, rsp.broadening_fwhm * EV_TO_HARTREE, num_broadened_samples, rsp.x_unit};
                                        for (int i = 0; i < num_broadened_samples; ++i) {
                                            ImPlotPoint point = getter(i, &data);
                                            max_cs = MAX(max_cs, point.y);
                                        }
                                        
                                        plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(0.0, max_cs);
                                        plot->Axes[ImAxis_Y2].FitExtents = ImPlotRange(0.0, max_y);
                                    }
                                }
                            }

                            ImPlot::EndPlot();
                        }
                    }
                }

                if (num_peaks > 0 && ImGui::TreeNodeEx("State Table", tree_flags)) {
                    // Table
                    static const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY |
                        ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_Sortable;

                    static const ImGuiTableColumnFlags columns_base_flags = ImGuiTableColumnFlags_DefaultSort;

                    ImVec2 table_size = { 0, 0 };

                    struct Column {
                        const char* label;
                        ImGuiTableColumnFlags flags;
                        const double* data;
                    };

                    Column cols[8];
                    int num_cols = 0;

                    // Populate colums with existing data
                    if (x_peaks) {
                        cols[num_cols++] = { x_unit_full_str[rsp.x_unit], columns_base_flags, x_peaks };
                    }
                    if (y_peaks_osc) {
                        cols[num_cols++] = { "Osc. Str.", columns_base_flags, y_peaks_osc };
                    }
                    if (y_peaks_cgs) {
                        cols[num_cols++] = { (const char*)u8"R (10⁻⁴⁰ cgs)", columns_base_flags, y_peaks_cgs };
                    }
                    if (y_peaks_tpa_trans_linear) {
                        cols[num_cols++] = { (const char*)u8"TPA Str. (lin) [a.u.]", columns_base_flags, y_peaks_tpa_trans_linear };
                    }
                    if (y_peaks_tpa_trans_circular) {
                        cols[num_cols++] = { (const char*)u8"TPA Str. (cir) [a.u.]", columns_base_flags, y_peaks_tpa_trans_circular };
                    }

                    if (ImGui::BeginTable("rsp table", num_cols + 1, flags, table_size, 0)) {
                        ImGui::TableSetupColumn("Index", columns_base_flags);
                        for (int col = 0; col < num_cols; ++col) {
                            ImGui::TableSetupColumn(cols[col].label, cols[col].flags);
                        }
                        ImGui::TableSetupScrollFreeze(0, 1);
                        ImGui::TableHeadersRow();

                        // TODO: Add sorting to the table once the vibs data structure is properly defined

                        ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);
                        ImGui::PushStyleColor(ImGuiCol_Header, IM_BLUE);

                        for (int row_n = 0; row_n < (int)num_frequencies; row_n++) {
                            ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                            bool is_sel = row_n == rsp.selected;
                            bool is_hov = row_n == rsp.hovered;

                            ImGui::TableNextRow(ImGuiTableRowFlags_None, 0);
                            ImGui::TableNextColumn();

                            if (is_sel) {
                                ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_BLUE);
                            } else {
                                ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);
                            }

                            char label[16];
                            snprintf(label, sizeof(label), "%i", row_n + 1);
                            if (ImGui::Selectable(label, is_sel || is_hov, selectable_flags)) {
                                rsp.selected = (rsp.selected == row_n) ? -1 : row_n;
                            }
                            ImGui::TableNextColumn();
                            for (int col = 0; col < num_cols; ++col) {
                                ImGui::Text("%12.6f", *(cols[col].data + row_n));
                                ImGui::TableNextColumn();
                            }

                            ImGui::PopStyleColor(1);
                        }
                        if (ImGui::IsWindowHovered() && ImGui::TableGetHoveredRow() > 0) {
                            rsp.hovered = ImGui::TableGetHoveredRow() - 1;
                        }

                        ImGui::PopStyleColor(2);
                        ImGui::EndTable();
                    }
                }
                ImGui::TreePop();
            }

            if (num_normal_modes > 0) {
                ImGui::SetNextItemOpen(true, ImGuiCond_Once);
                if (ImGui::TreeNode("Vibrational Spectroscopy")) {
                    size_t num_atoms = md_vlx_number_of_atoms(vlx);

                    // Broadening gamma limits
                    static const double gamma_min = 1.0;
                    static const double gamma_max = 100.0;

                    // Frequency scaling factor limits
                    static const double scale_min = 0.85;
                    static const double scale_max = 1.05;

                    // Explicit x axis limits for linking
                    static double x_min = 0.0;
                    static double x_max = 10.0;

                    const float avail_width = ImGui::GetContentRegionAvail().x;
                    ImGui::PushItemWidth(MIN(avail_width, 200));

                    bool refit = false;

                    refit |= ImGui::SliderScalar((const char*)u8"Broadening FWHM (cm⁻¹)", ImGuiDataType_Double, &vib.broadening_fwhm, &gamma_min, &gamma_max, "%.3f", ImGuiSliderFlags_Logarithmic);
                    refit |= ImGui::Combo("Broadening mode", (int*)(&vib.broadening_mode), broadening_mode_str, IM_ARRAYSIZE(broadening_mode_str));
                    refit |= ImGui::SliderScalar("Frequency scaling factor", ImGuiDataType_Double, &vib.freq_scaling_factor, &scale_min, &scale_max, "%.4f");
                    ImGui::PopItemWidth();

                    ImGui::SetItemTooltip("Frequency scaling factor");

                    auto spectrum_getter_lorentzian = [](int idx, void* user_data) -> ImPlotPoint {
                        const SpectrumGetterData* data = (const SpectrumGetterData*)user_data;

                        // Construct an x value for idx as a linear space between x_min and x_max.
                        double x = data->x_min + (data->x_max - data->x_min) * idx / (data->num_samples - 1);

                        // Compute the corresponding y value by broadening the peaks at this x value.
                        double y = lorentzian_vib(x, data->x_peaks, data->y_peaks, data->num_peaks, data->broadening_fwhm);
                        return ImPlotPoint{x, y};
                    };

                    auto spectrum_getter_gaussian = [](int idx, void* user_data) -> ImPlotPoint {
                        const SpectrumGetterData* data = (const SpectrumGetterData*)user_data;

                        // Construct an x value for idx as a linear space between x_min and x_max.
                        double x = data->x_min + (data->x_max - data->x_min) * idx / (data->num_samples - 1);

                        // Compute the corresponding y value by broadening the peaks at this x value.
                        double y = gaussian_vib(x, data->x_peaks, data->y_peaks, data->num_peaks, data->broadening_fwhm);
                        return ImPlotPoint{ x, y };
                    };

					ImPlotGetter getter = nullptr;
					if (vib.broadening_mode == BROADENING_MODE_LORENTZIAN) {
						getter = spectrum_getter_lorentzian;
					} else {
						getter = spectrum_getter_gaussian;
					}

                    const double* x_peaks_raw = md_vlx_vib_frequencies(vlx);
                    const double* y_peaks_ir = md_vlx_vib_ir_intensities(vlx);
                        
                    size_t num_external_frequencies = md_vlx_vib_number_of_external_frequencies(vlx);
                    const double* external_frequencies = md_vlx_vib_external_frequencies(vlx);
                    const double** y_raman_activities = NULL;
                    if (num_external_frequencies > 0) {
                        y_raman_activities = md_temp_alloc_array(temp, const double*, num_external_frequencies);
                        for (size_t i = 0; i < num_external_frequencies; ++i) {
                            y_raman_activities[i] = md_vlx_vib_raman_activities(vlx, i);
                        }
                    }

                    double* x_peaks = md_temp_alloc_array(temp, double, num_normal_modes);
                    for (size_t i = 0; i < num_normal_modes; ++i) {
                        x_peaks[i] = x_peaks_raw[i] * vib.freq_scaling_factor;
                    }

                    uint64_t x_hash = md_hash64(x_peaks_raw, sizeof(double) * num_normal_modes, 0);
                    if (vib.hash != x_hash) {
                        vib.hash = x_hash;
                        refit = true;
                    }

				    const size_t num_samples = 512;

                    /*
                    // These do not need to be exposed as the user can right click the plot and invert the axes
                    ImGui::Checkbox("Invert X", &vib.invert_x);
                    ImGui::SameLine();
                    ImGui::Checkbox("Invert Y", &vib.invert_y);
                    */

                    ImPlotAxisFlags x_flags = vib.invert_x ? ImPlotAxisFlags_Invert : 0;
                    ImPlotAxisFlags y_flags = vib.invert_y ? ImPlotAxisFlags_Invert : 0;

                    if (x_peaks && y_peaks_ir && num_normal_modes) {
                        if (ImGui::TreeNodeEx("IR", tree_flags)) {

                            ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

                            if (refit) {
                                ImPlot::SetNextAxesToFit();
                            }

                            if (ImPlot::BeginPlot("IR Spectrum")) {
                                ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);

                                ImPlot::SetupAxis(ImAxis_X1, (const char*)u8"Harmonic Frequency (cm⁻¹)", x_flags);
                                ImPlot::SetupAxis(ImAxis_Y1, (const char*)u8"Absorption Cross Section", y_flags);
                                ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"IR Intensity (km/mol)", y_flags | ImPlotAxisFlags_AuxDefault);

                                ImPlot::SetupFinish();

                                ImPlot::SetAxis(ImAxis_Y1);
                                SpectrumGetterData getter_data = { x_peaks, y_peaks_ir, num_normal_modes, x_min, x_max, vib.broadening_fwhm, num_samples };
                                ImPlot::PlotLineG("Broadening", getter, &getter_data, num_samples);

                                ImPlot::SetAxis(ImAxis_Y2);
                                plot_peaks("IR Intensity", x_peaks, y_peaks_ir, num_normal_modes, vib.selected, vib.hovered);

                                if (ImPlot::FitThisFrame()) {
                                    ImPlotPlot* plot = ImPlot::GetCurrentPlot();
                                    if (plot) {
                                        double min_x = DBL_MAX, max_x = -DBL_MAX, max_y = -DBL_MAX;
                                        for (size_t i = 0; i < num_normal_modes; ++i) {
                                            min_x = MIN(min_x, x_peaks[i]);
                                            max_x = MAX(max_x, x_peaks[i]);
                                            max_y = MAX(max_y, y_peaks_ir[i]);
                                        }

                                        double normalization_factor = 1.0;
                                        if (vib.broadening_mode == BROADENING_MODE_LORENTZIAN) {
                                            normalization_factor = 1.0 / PI;
                                        }
                                        else {
                                            normalization_factor = sqrt(2.0) / (sigma_from_fwhm(vib.broadening_fwhm) * sqrt(PI));
                                        }

                                        plot->Axes[ImAxis_X1].FitExtents = ImPlotRange(min_x, max_x);
                                        plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(0.0, 1.2 * max_y * normalization_factor);
                                        plot->Axes[ImAxis_Y2].FitExtents = ImPlotRange(0.0, 1.1 * max_y);
                                    }
                                }

                                ImPlot::EndPlot();
                            }
                        }
                    }

                    if (num_external_frequencies > 0) {
                        if (ImGui::TreeNodeEx("Raman Activity", tree_flags)) {
                            
                            if (num_external_frequencies > 1) {
                                ImGui::PushItemWidth(MIN(avail_width, 200));

                                char label[32];
                                snprintf(label, sizeof(label), "%.4f", external_frequencies[vib.external_frequency_index]);
                                if (ImGui::BeginCombo("External Frequency", label)) {
                                    for (size_t i = 0; i < num_external_frequencies; ++i) {
                                        bool is_selected = ((size_t)vib.external_frequency_index == i);

                                        snprintf(label, sizeof(label), "%.4f", external_frequencies[i]);
                                        if (ImGui::Selectable(label, is_selected)) {
                                            vib.external_frequency_index = (int)i;
                                        }
                                        if (is_selected) {
                                            ImGui::SetItemDefaultFocus();
                                        }
                                    }
                                    ImGui::EndCombo();
                                }

                                ImGui::PopItemWidth();
                            }

                            ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

                            if (refit) {
                                ImPlot::SetNextAxesToFit();
                            }

                            const double* y_peaks_raman = y_raman_activities[vib.external_frequency_index];

                            char label[32];
                            double freq = external_frequencies[vib.external_frequency_index];
                            if (num_external_frequencies == 1 && freq == 0.0) {
                                snprintf(label, sizeof(label), "Raman Spectrum (static)");
                            } else {
                                snprintf(label, sizeof(label), "Raman Spectrum (%.4f)", freq);
                            }

                            if (ImPlot::BeginPlot(label)) {
                                ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                                ImPlot::SetupAxis(ImAxis_X1, (const char*)u8"Harmonic Frequency (cm⁻¹)", x_flags);
                                ImPlot::SetupAxis(ImAxis_Y1, (const char*)u8"Scattering Cross Section", y_flags);
                                ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"Raman Activity (Å⁴/amu)", y_flags | ImPlotAxisFlags_AuxDefault);
                                ImPlot::SetupFinish();

                                ImPlot::SetAxis(ImAxis_Y1);
								SpectrumGetterData data = { x_peaks, y_peaks_raman, num_normal_modes, x_min, x_max, vib.broadening_fwhm, num_broadened_samples };
                                ImPlot::PlotLineG("Broadening", getter, &data, num_broadened_samples);

                                ImPlot::SetAxis(ImAxis_Y2);
                                plot_peaks("Raman Activity", x_peaks, y_peaks_raman, num_normal_modes, vib.selected, vib.hovered);

                                if (ImPlot::FitThisFrame() || refit) {
                                    ImPlotPlot* plot = ImPlot::GetCurrentPlot();
                                    if (plot) {
                                        double min_x = DBL_MAX, max_x = -DBL_MAX, max_y = -DBL_MAX;
                                        for (size_t i = 0; i < num_normal_modes; ++i) {
                                            min_x = MIN(min_x, x_peaks[i]);
                                            max_x = MAX(max_x, x_peaks[i]);
                                            max_y = MAX(max_y, y_peaks_raman[i]);
                                        }

                                        plot->Axes[ImAxis_X1].FitExtents = ImPlotRange(min_x, max_x);
                                        plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(0.0, 1.2 * max_y / (PI * vib.broadening_fwhm));
                                        plot->Axes[ImAxis_Y2].FitExtents = ImPlotRange(0.0, 1.1 * max_y);
                                    }
                                }

                                ImPlot::EndPlot();
                            }
                        }
                    }

                    ImGui::PushItemWidth(MIN(avail_width, 200));

                    int idx = vib.selected + 1;
                    ImGui::SliderInt("Index", &idx, 0, (int)num_normal_modes);
                    vib.selected = CLAMP(idx - 1, -1, (int)num_normal_modes - 1);
                    ImGui::SliderFloat((const char*)"Amplitude", &vib.displacement_amp_scl, 0.25f, 5.0f);
                    ImGui::SetItemTooltip("Displacement Amplitude Scale");
                    ImGui::SliderFloat((const char*)"Frequency", &vib.displacement_freq_scl, 0.25f, 5.0f);
                    ImGui::SetItemTooltip("Displacement Frequency Scale");
                    ImGui::PopItemWidth();

                    // Table
                    static const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY |
                                                         ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_Sortable;

                    static const ImGuiTableColumnFlags columns_base_flags = ImGuiTableColumnFlags_DefaultSort;

                    ImVec2 table_size = {0, 0};

                    int num_cols = 1 + (int)num_external_frequencies;
					if (x_peaks) {
						num_cols++;
					}
					if (y_peaks_ir) {
						num_cols++;
					}
                    if (ImGui::BeginTable("table", num_cols, flags, table_size, 0)) {
                        ImGui::TableSetupColumn("Vibration mode", columns_base_flags, 0.0f);
                        if (x_peaks) {
                            ImGui::TableSetupColumn("Harmonic Frequency", columns_base_flags, 0.0f);
                        }
                        if (y_peaks_ir) {
                            ImGui::TableSetupColumn("IR Intensity", columns_base_flags, 0.0f);
                        }
                        for (size_t i = 0; i < num_external_frequencies; ++i) {
                            char label[32];
                            double freq = external_frequencies[i];
                            snprintf(label, sizeof(label), "Raman Activity (%.4f)", freq);
                            ImGui::TableSetupColumn(label, columns_base_flags, 0.0f);
                        }
                        ImGui::TableSetupScrollFreeze(0, 1);
                        ImGui::TableHeadersRow();

                        // TODO: Add sorting to the table once the vibs data structure is properly defined

                        ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);
                        ImGui::PushStyleColor(ImGuiCol_Header, IM_BLUE);
                        
                        for (int row_n = 0; row_n < (int)num_normal_modes; row_n++) {

                            ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                            bool is_sel = row_n == vib.selected;
                            bool is_hov = row_n == vib.hovered;
                            
                            ImGui::TableNextRow(ImGuiTableRowFlags_None, 0);
                            ImGui::TableNextColumn();

                            if (is_sel) {
                                ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_BLUE);
                            } else {
                                ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);
                            }

                            char label[16];
                            snprintf(label, sizeof(label), "%i", row_n + 1);
                            if (ImGui::Selectable(label, is_sel || is_hov, selectable_flags)) {
                                vib.selected = (vib.selected == row_n) ? -1 : row_n;
                            }
                            if (x_peaks) {
                                ImGui::TableNextColumn();
                                ImGui::Text("%12.6f", x_peaks[row_n]);
                            }
							if (y_peaks_ir) {
                                ImGui::TableNextColumn();
                                ImGui::Text("%12.6f", y_peaks_ir[row_n]);
                            }
                            for (size_t i = 0; i < num_external_frequencies; ++i) {
                                ImGui::TableNextColumn();
                                ImGui::Text("%12.6f", y_raman_activities[i][row_n]);
                            }

                            ImGui::PopStyleColor(1);
                        }
                        if (ImGui::IsWindowHovered() && ImGui::TableGetHoveredRow() > 0) {
                            vib.hovered = ImGui::TableGetHoveredRow() - 1;
                        } else {
                            vib.hovered = -1;
                        }

                        ImGui::PopStyleColor(2);
                        ImGui::EndTable();
                    }

                    const dvec3_t* atom_coord = md_vlx_atom_coordinates(vlx);

                    // Check selected state
                    if (vib.selected != -1) {
                        // Animate
                        vib.t += state.app.timing.delta_s * vib.displacement_freq_scl * 8.0;
                        const dvec3_t* norm_modes = md_vlx_vib_normal_mode(vlx, vib.selected);

                        if (norm_modes) {

#if 0
                            md_temp_scope_t temp = md_temp_begin();
                            defer { md_temp_end(temp); };

                            const double scl = vib.displacement_amp_scl * 1.0;

                            float* vx = (float*)md_temp_alloc(temp, sizeof(float) * num_atoms);
                            float* vy = (float*)md_temp_alloc(temp, sizeof(float) * num_atoms);
                            float* vz = (float*)md_temp_alloc(temp, sizeof(float) * num_atoms);

                            for (size_t i = 0; i < num_atoms; i++) {
                                /*
                                state.mold.sys.atom.x[i] = (float)(atom_coord[i].x + norm_modes[i].x * scl);
                                state.mold.sys.atom.y[i] = (float)(atom_coord[i].y + norm_modes[i].y * scl);
                                state.mold.sys.atom.z[i] = (float)(atom_coord[i].z + norm_modes[i].z * scl);
                                */
                                vx[i] = (float)(norm_modes[i].x * scl);
                                vy[i] = (float)(norm_modes[i].y * scl);
                                vz[i] = (float)(norm_modes[i].z * scl);
                            }

                            //md_gl_mol_set_atom_position(state.mold.gl_mol, 0, num_atoms, state.mold.sys.atom.x, state.mold.sys.atom.y, state.mold.sys.atom.z, 0);
                            md_gl_mol_set_atom_velocity(state.mold.gl_mol, 0, num_atoms, vx, vy, vz, 0);
                            

#else
                            const double scl = vib.displacement_amp_scl * 0.25 * sin(vib.t);
                            for (size_t i = 0; i < num_atoms; i++) {
                                state.mold.state.x[i] = (float)(atom_coord[i].x + norm_modes[i].x * scl);
                                state.mold.state.y[i] = (float)(atom_coord[i].y + norm_modes[i].y * scl);
                                state.mold.state.z[i] = (float)(atom_coord[i].z + norm_modes[i].z * scl);
                            }
                            if (vib.displace_aos) {
                            }
                            state.mold.dirty_gpu_buffers |= MolBit_DirtyPosition;
#endif
                            vib.coord_modified = true;

                            if (vib.displace_aos) {
                                viamd::event_system_broadcast_event(viamd::EventType_ViamdSystemStateChanged, viamd::EventPayloadType_ApplicationState, &state);
                                size_t num_reps = md_array_size(state.representation.reps);
                                for (size_t i = 0; i < num_reps; ++i) {
                                    if (state.representation.reps[i].type == RepresentationType::ElectronicStructure) {
                                        flag_representation_as_dirty(&state.representation.reps[i]);
                                    }
                                }
                            }
                        }
                    }
                    // If all is deselected, reset coords once
                    else if (vib.coord_modified) {
                        for (size_t i = 0; i < num_atoms; i++) {
                            state.mold.state.x[i] = (float)atom_coord[i].x;
                            state.mold.state.y[i] = (float)atom_coord[i].y;
                            state.mold.state.z[i] = (float)atom_coord[i].z;
                        }
                        viamd::event_system_broadcast_event(viamd::EventType_ViamdSystemStateChanged, viamd::EventPayloadType_ApplicationState, &state);
                        state.mold.dirty_gpu_buffers |= MolBit_DirtyPosition | MolBit_ClearVelocity;
                        vib.coord_modified = false;
                    }

                    ImGui::TreePop();
                }
            }
        }
        if (!ImGui::IsWindowHovered()) {
            rsp.hovered = -1;
        }
        ImGui::End();

        if (rsp.show_export_window) { draw_rsp_spectra_export_window(state); }
    }

    void draw_orb_window(ApplicationState& state) {
        // The XPS plot hands over the MO of the peak under the cursor. draw_rsp_window runs after
        // this one in the frame tick, so what is read here is one frame old - invisible for a hover
        // highlight, and cheaper than reordering the windows. Consumed unconditionally, ahead of the
        // early returns, so a hover that has gone away - or a Response window that stopped drawing
        // altogether - cannot leave a row lit up forever.
        const int32_t xps_highlight_mo = xps.highlight_mo_idx;
        xps.highlight_mo_idx = -1;

        if (!orb.show_window) return;
        if (num_molecular_orbitals() == 0) return;
        ImGui::SetNextWindowSize({600, 300}, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Orbital Grid", &orb.show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {

            const double* occ_alpha = md_vlx_scf_mo_occupancy(vlx, MD_VLX_SPIN_ALPHA);
            const double* occ_beta = md_vlx_scf_mo_occupancy(vlx, MD_VLX_SPIN_BETA);
            const double* ene_alpha = md_vlx_scf_mo_energy(vlx, MD_VLX_SPIN_ALPHA);
            const double* ene_beta = md_vlx_scf_mo_energy(vlx, MD_VLX_SPIN_BETA);

            const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();
            md_vlx_scf_type_t type = md_vlx_scf_type(vlx);

            int num_x = (type == MD_VLX_SCF_UNRESTRICTED) ? 2 : orb.num_x;
            int num_y = orb.num_y;

            int orb_homo_idx = MAX(homo_idx[0], homo_idx[1]);
            int orb_lumo_idx = MIN(lumo_idx[0], lumo_idx[1]);
            if (type == MD_VLX_SCF_RESTRICTED_OPENSHELL) {
                orb_lumo_idx = MAX(lumo_idx[0], lumo_idx[1]);
            }

            int num_mos = num_x * num_y;
            int beg_mo_idx = orb.mo_idx - num_mos / 2 + (num_mos % 2 == 0 ? 1 : 0);
            int window_size = num_mos;
            if (type == MD_VLX_SCF_UNRESTRICTED) {
                beg_mo_idx = orb.mo_idx - num_y / 2 + (num_y % 2 == 0 ? 1 : 0);
                window_size = num_y;
            }

            // LEFT PANE
            {
                ImGui::BeginChild("left pane", ImVec2(300, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
                ImGui::SameLine();

                ImGui::PushItemWidth(-1);
                ImGui::BeginGroup();

                ImGui::SliderInt("##Rows", &orb.num_y, 1, 4);
                if (type == MD_VLX_SCF_UNRESTRICTED) ImGui::PushDisabled();
                ImGui::SliderInt("##Cols", &orb.num_x, 1, 4);
                if (type == MD_VLX_SCF_UNRESTRICTED) ImGui::PopDisabled();

                const double iso_min = 1.0e-4;
                const double iso_max = 5.0;
                double iso_val = orb.iso.values[0];
                ImGui::TextUnformatted(electronic_structure_iso_value_label());
                ImGui::SliderScalar("##Iso Value", ImGuiDataType_Double, &iso_val, &iso_min, &iso_max, "%.6f", ImGuiSliderFlags_Logarithmic);
                ImGui::SetItemTooltip("Visual surface threshold. Orbital grids render both +isovalue and -isovalue surfaces.");

                orb.iso.values[0] = (float)iso_val;
                orb.iso.values[1] = -(float)iso_val;
                orb.iso.count = 2;

                ImGui::ColorEdit4("##Color Positive", orb.iso.colors[0].elem);
                ImGui::SetItemTooltip("Color Positive");
                ImGui::ColorEdit4("##Color Negative", orb.iso.colors[1].elem);
                ImGui::SetItemTooltip("Color Negative");

                enum {
                    Col_Idx,
                    Col_Occ_Alpha,
                    Col_Occ_Beta,
                    Col_Ene_Alpha,
                    Col_Ene_Beta,
                };

                const char* btn_text = "Goto HOMO";
                if (type == MD_VLX_SCF_RESTRICTED_OPENSHELL && (homo_idx[0] != homo_idx[1])) {
                    btn_text = "Goto SOMO";
                }

                if (ImGui::IsWindowAppearing()) {
                    orb.scroll_to_idx = orb.mo_idx;
                }
                if (ImGui::Button(btn_text, ImVec2(-1, 0))) {
                    orb.scroll_to_idx = homo_idx[0];
                }

                int num_cols = (type == MD_VLX_SCF_UNRESTRICTED) ? 5 : 3;

                const ImGuiTableFlags flags =
                    ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable | ImGuiTableFlags_RowBg |
                    ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersV | ImGuiTableFlags_NoBordersInBody | ImGuiTableFlags_ScrollY;
                if (ImGui::BeginTable("Molecular Orbitals", num_cols, flags))//, ImVec2(0.0f, TEXT_BASE_HEIGHT * 15), 0.0f))
                {
                    // Declare columns
                    // We use the "user_id" parameter of TableSetupColumn() to specify a user id that will be stored in the sort specifications.
                    // This is so our sort function can identify a column given our own identifier. We could also identify them based on their index!
                    // Demonstrate using a mixture of flags among available sort-related flags:
                    // - ImGuiTableColumnFlags_DefaultSort
                    // - ImGuiTableColumnFlags_NoSort / ImGuiTableColumnFlags_NoSortAscending / ImGuiTableColumnFlags_NoSortDescending
                    // - ImGuiTableColumnFlags_PreferSortAscending / ImGuiTableColumnFlags_PreferSortDescending
                    if (type == MD_VLX_SCF_UNRESTRICTED) {
                        ImGui::TableSetupColumn("MO", ImGuiTableColumnFlags_DefaultSort | ImGuiTableColumnFlags_WidthFixed, 0.0f, Col_Idx);
                        ImGui::TableSetupColumn((const char*)u8"Occ. α", ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed, 0.0f, Col_Occ_Alpha);
                        ImGui::TableSetupColumn((const char*)u8"Occ. β", ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed, 0.0f, Col_Occ_Beta);
                        ImGui::TableSetupColumn((const char*)u8"Ene. α", ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed, 0.0f, Col_Ene_Alpha);
                        ImGui::TableSetupColumn((const char*)u8"Ene. β", ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed, 0.0f, Col_Ene_Beta);
                    }
                    else {
                        ImGui::TableSetupColumn("MO", ImGuiTableColumnFlags_DefaultSort | ImGuiTableColumnFlags_WidthFixed, 0.0f, Col_Idx);
                        ImGui::TableSetupColumn("Occupancy", ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed, 0.0f, Col_Occ_Alpha);
                        ImGui::TableSetupColumn("Energy", ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed, 0.0f, Col_Ene_Alpha);
                    }
                    ImGui::TableSetupScrollFreeze(0, 1); // Make row always visible
                    ImGui::TableHeadersRow();

                    for (int n = (int)num_molecular_orbitals() - 1; n >= 0; n--) {
                        ImGui::PushID(n + 1);
                        ImGui::TableNextRow();
                        if (n == xps_highlight_mo) {
                            ImGui::TableSetBgColor(ImGuiTableBgTarget_RowBg1, ImGui::GetColorU32(COLOR_ROW_HIGHLIGHT));
                        }
                        bool is_selected = (beg_mo_idx <= n && n < beg_mo_idx + window_size);
                        ImGui::TableNextColumn();
                        if (orb.scroll_to_idx != -1 && n == orb.scroll_to_idx) {
                            orb.scroll_to_idx = -1;
                            ImGui::SetScrollHereY();
                        }

                        char buf[32];
                        snprintf(buf, sizeof(buf), "%i", n + 1);

                        ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                        if (ImGui::Selectable(buf, is_selected, selectable_flags)) {
                            if (orb.mo_idx != n) {
                                orb.mo_idx = n;
                            }
                        }
                        ImGui::TableNextColumn();
                        if (occ_alpha) {
                            const char* lbl = "";
                            double occ = occ_alpha[n];
                            if (type == MD_VLX_SCF_UNRESTRICTED) {
                                if (n == homo_idx[0] && n == orb_homo_idx) {
                                    lbl = "HOMO";
                                }
                                else if (n == lumo_idx[0] && n == orb_lumo_idx) {
                                    lbl = "LUMO";
                                }
                            }
                            else {
                                if (occ_beta) {
                                    occ += occ_beta[n];
                                }
                                if (n == orb_homo_idx) {
                                    lbl = (homo_idx[0] == homo_idx[1]) ? "HOMO" : "SOMO";
                                }
                                else if (n == orb_lumo_idx) {
                                    lbl = "LUMO";
                                }
                            }
                            ImGui::Text("%.1f %s", occ, lbl);
                        }
                        else {
                            ImGui::Text("-");
                        }
                        if (type == MD_VLX_SCF_UNRESTRICTED) {
                            ImGui::TableNextColumn();
                            if (occ_beta) {
                                const char* lbl = "";
                                if (n == homo_idx[1] && n == orb_homo_idx) {
                                    lbl = "HOMO";
                                }
                                else if (n == lumo_idx[1] && n == orb_lumo_idx) {
                                    lbl = "LUMO";
                                }
                                ImGui::Text("%.1f %s", occ_beta[n], lbl);
                            }
                            else {
                                ImGui::Text("-");
                            }
                        }
                        ImGui::TableNextColumn();
                        if (ene_alpha) {
                            ImGui::Text("%.4f", ene_alpha[n]);
                        }
                        else {
                            ImGui::Text("-");
                        }
                        if (type == MD_VLX_SCF_UNRESTRICTED) {
                            ImGui::TableNextColumn();
                            if (ene_beta) {
                                ImGui::Text("%.4f", ene_beta[n]);
                            }
                            else {
                                ImGui::Text("-");
                            }
                        }
                        ImGui::PopID();
                    }
                    ImGui::EndTable();
                }

                ImGui::EndGroup();
                ImGui::PopItemWidth();

                ImGui::EndChild();
            }

            ImGui::SameLine();

            // These represent the new mo_idx we want to have in each slot
            int vol_mo_idx[16] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };
            md_vlx_spin_t vol_mo_type[16] = {};
            for (int i = 0; i < num_mos; ++i) {
                if (type == MD_VLX_SCF_UNRESTRICTED) {
                    vol_mo_idx[i] = beg_mo_idx + i / 2;
                    vol_mo_type[i] = (i & 1) ? MD_VLX_SPIN_ALPHA : MD_VLX_SPIN_BETA;
                }
                else {
                    vol_mo_idx[i] = beg_mo_idx + i;
                    vol_mo_type[i] = MD_VLX_SPIN_ALPHA;
                }
            }

            int job_queue[16];
            int num_jobs = 0;
            // Find and reuse volume data from existing slots (if applicable)
            // If there is no existing volume, we queue up a new job
            for (int i = 0; i < num_mos; ++i) {
                // Check if we already have that entry in the correct slot
                if (orb.vol_mo_idx[i] == vol_mo_idx[i] && orb.vol_mo_type[i] == vol_mo_type[i]) continue;

                // Try to find the entry in the existing list
                bool found = false;
                for (int j = 0; j < num_mos; ++j) {
                    if (i == j) continue;
                    if (vol_mo_idx[i] == orb.vol_mo_idx[j]) {
                        // Swap to correct location
                        ImSwap(orb.vol[i], orb.vol[j]);
                        ImSwap(orb.vol_mo_idx[i], orb.vol_mo_idx[j]);
                        ImSwap(orb.vol_mo_type[i], orb.vol_mo_type[j]);
                        found = true;
                        break;
                    }
                }

                // If not found, put in job queue to compute the volume
                if (!found) {
                    job_queue[num_jobs++] = i;
                }
            }

            if (num_jobs > 0) {
                const float samples_per_unit_length = DEFAULT_SAMPLES_PER_ANGSTROM * BOHR_TO_ANGSTROM;
                md_grid_t grid { 0 };
#if USE_AABB_GRID
                init_grid(&grid, mat3_ident(), aabb.min_ext, aabb.max_ext, samples_per_unit_length);
#else
                init_grid(&grid, oabb.orientation, oabb.min_ext, oabb.max_ext, samples_per_unit_length);
#endif
                const int num_tot_mos = (int)num_molecular_orbitals();

                for (int i = 0; i < num_jobs; ++i) {
                    int slot_idx = job_queue[i];
                    int mo_idx = vol_mo_idx[slot_idx];
                    md_vlx_spin_t mo_type = vol_mo_type[slot_idx];
                    orb.vol_mo_idx[slot_idx] = mo_idx;
                    orb.vol_mo_type[slot_idx] = mo_type;

                    if (-1 < mo_idx && mo_idx < num_tot_mos) {
                        init_volume(&orb.vol[slot_idx], grid);
                        evaluate_mo(state, orb.vol[slot_idx].tex_id, grid, mo_type, mo_idx, MD_GTO_OP_SET);
                    }
                }
            }

            ImVec2 canvas_sz = ImMax(ImVec2(50.0f, 50.0f), ImGui::GetContentRegionAvail());   // Resize canvas to what's available

            // Animate camera towards targets
            const double dt = state.app.timing.delta_s;
            camera_animate(&orb.camera, orb.target, dt);

            ImGui::Dummy(canvas_sz);

            // Draw border and background color
            ImGuiIO& io = ImGui::GetIO();

            ImVec2 canvas_min = ImGui::GetItemRectMin();
            ImVec2 canvas_max = ImGui::GetItemRectMax();
            ImVec2 orb_win_sz = ImFloor(canvas_sz / ImVec2((float)num_x, (float)num_y));

            const float aspect = orb_win_sz.x / orb_win_sz.y;

            //const mat4_t world_to_view = camera_world_to_view_matrix(orb.camera);
            //const mat4_t view_to_clip  = camera_view_to_clip_matrix_persp(orb.camera, aspect);

            const mat4_t view_to_world = camera_view_to_world_matrix(orb.camera);
            const mat4_t clip_to_view  = camera_clip_to_view_matrix_persp(orb.camera, aspect);
            const mat4_t clip_to_world = mat4_mul(view_to_world, clip_to_view);

            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->AddRectFilled(canvas_min, canvas_max, IM_COL32(255, 255, 255, 255));

            const int num_orbs = (int)num_molecular_orbitals();

            for (int i = 0; i < num_mos; ++i) {
                ImGui::PushID(i);
                defer { ImGui::PopID(); };

                int mo_idx = beg_mo_idx + i;
                if (type == MD_VLX_SCF_UNRESTRICTED) {
                    mo_idx = beg_mo_idx + i / 2;
                }
                int x = num_x - i % num_x - 1;
                int y = num_y - i / num_x - 1;
                if (-1 < mo_idx && mo_idx < num_orbs) {
                    const ImVec2 p0 = canvas_min + orb_win_sz * ImVec2((float)(x + 0), (float)(y + 0));
                    const ImVec2 p1 = canvas_min + orb_win_sz * ImVec2((float)(x + 1), (float)(y + 1));
                    const ImVec2 text_pos_bl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p1.y - TEXT_BASE_HEIGHT);
                    const ImVec2 text_pos_tl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p0.y + TEXT_BASE_HEIGHT * 0.25f);
                    const ImVec2 text_pos_br = ImVec2(p1.x - TEXT_BASE_HEIGHT * 0.5f, p1.y - TEXT_BASE_HEIGHT);

                    const ImVec2 sz = p1 - p0;

                    ImGui::SetCursorScreenPos(p0);
                    InteractionSurfaceState surface_state = interaction_surface(interaction_surface_orb, vec_cast(sz), InteractionSurfaceFlags_NoRegionSelect);

                    PickingHit hit = {};
                    if (surface_state.hovered) {
                        InteractionSurfaceHitArgs hit_args = {
                            .picking_surface = &orb.picking_surface,
                            .picking_handler = state.picking_handler,
                            .fbo = orb.gbuf.fbo,
                            .width = orb.gbuf.width,
                            .height = orb.gbuf.height,
                            .clip_to_world = clip_to_world,
                        };
                        interaction_surface_hit_extract(&hit, surface_state, hit_args);

                        InteractionSurfaceEvent event = {};
                        interaction_surface_event_extract(&event, surface_state, hit);

                        viamd::event_system_broadcast_event(viamd::EventType_ViamdInteractionSurface, viamd::EventPayloadType_InteractionSurfaceEvent, &event);
                    }

                    InteractionSurfaceViewTransformArgs view_args = {
                        .camera = orb.camera,
                        .trackball_param = state.view.trackball_param,
                    };

                    InteractionSurfaceViewTransformResult view_result = interaction_surface_view_transform_apply(&orb.target, surface_state, view_args);
                    if (view_result.reset_requested) {
                        ViewTransform reset_transform = default_view;
                        if (hit.depth < 1.0f) {
                            reset_transform.distance = orb.target.distance;
                            reset_transform.orientation = orb.camera.orientation;
                            reset_transform.position = hit.world_pos + orb.camera.orientation * vec3_set(0, 0, orb.target.distance);
                        }
                        orb.target = reset_transform;
                    }


                    const char* lbl = "";
                    if (type == MD_VLX_SCF_UNRESTRICTED) {
                        int j = (i & 1) ? 0 : 1;
                        if (mo_idx == orb_homo_idx && orb_homo_idx == homo_idx[j]) {
                            lbl = "(HOMO)";
                        }
                        else if (mo_idx == orb_lumo_idx && orb_lumo_idx == lumo_idx[j]) {
                            lbl = "(LUMO)";
                        }
                    }
                    else {
                        if (mo_idx == orb_homo_idx) {
                            lbl = (homo_idx[0] == homo_idx[1]) ? "(HOMO)" : "(SOMO)";
                        }
                        else if (mo_idx == orb_lumo_idx) {
                            lbl = "(LUMO)";
                        }
                    }

                    char buf[32];
                    snprintf(buf, sizeof(buf), "%i %s", mo_idx + 1, lbl);
                    draw_list->AddImage((ImTextureID)(intptr_t)orb.gbuf.tex.transparency, p0, p1, { 0,1 }, { 1,0 });
                    draw_list->AddImage((ImTextureID)(intptr_t)orb.iso_tex[i], p0, p1, { 0,1 }, { 1,0 });
                    draw_list->AddText(text_pos_bl, ImColor(0, 0, 0), buf);

                    if (type == MD_VLX_SCF_UNRESTRICTED) {
                        draw_list->AddText(text_pos_tl, ImColor(0, 0, 0), (i & 1) ? (const char*)u8"α" : (const char*)u8"β");
                        snprintf(buf, sizeof(buf), "%.4f", (i & 1) ? ene_alpha[mo_idx] : ene_beta[mo_idx]);
                    }
                    else {
                        snprintf(buf, sizeof(buf), "%.4f", ene_alpha[mo_idx]);
                    }
                    float width = ImGui::CalcTextSize(buf).x;
                    draw_list->AddText(text_pos_br - ImVec2(width, 0), ImColor(0, 0, 0), buf);
                }
            }

            for (int x = 1; x < orb.num_x; ++x) {
                ImVec2 p0 = { canvas_min.x + orb_win_sz.x * x, canvas_min.y };
                ImVec2 p1 = { canvas_min.x + orb_win_sz.x * x, canvas_max.y };
                draw_list->AddLine(p0, p1, IM_COL32(0, 0, 0, 255));
            }

            for (int y = 1; y < orb.num_y; ++y) {
                ImVec2 p0 = { canvas_min.x, canvas_min.y + orb_win_sz.y * y };
                ImVec2 p1 = { canvas_max.x, canvas_min.y + orb_win_sz.y * y };
                draw_list->AddLine(p0, p1, IM_COL32(0, 0, 0, 255));
            }

            const ImVec2 origin(canvas_min.x, canvas_min.y);  // Lock scrolled origin
            const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

            int width  = MAX(1, (int)(orb_win_sz.x * io.DisplayFramebufferScale.x));
            int height = MAX(1, (int)(orb_win_sz.y * io.DisplayFramebufferScale.y));

            auto& gbuf = orb.gbuf;
            if ((int)gbuf.width != width || (int)gbuf.height != height) {
                gbuffer_init(&gbuf, width, height);
                for (int i = 0; i < num_mos; ++i) {
                    gl::init_texture_2D(orb.iso_tex + i, width, height, GL_RGBA8);
                }
            }

            if (orb.show_coordinate_system_widget) {
                float ext = MIN(orb_win_sz.x, orb_win_sz.y) * 0.4f;
                float pad = 20.0f;
                ImVec2 size = ImVec2(ext, ext);
                ImGui::SetCursorScreenPos(ImVec2(canvas_min.x + pad, canvas_max.y - ext - pad));

                quat_t out_orientation = orb.target.orientation;
                if (ImGui::CoordinateSystemWidget(&out_orientation, orb.camera.orientation, size)) {
                    const vec3_t look_at = camera_get_look_at(orb.target);
                    orb.target.orientation = quat_normalize(out_orientation);
                    orb.target.position = camera_position_from_look_at(look_at, orb.target.orientation, orb.target.distance);
                }
            }

            if (gl_rep.id) {
                const float aspect_ratio = orb_win_sz.x / orb_win_sz.y;
                mat4_t view_mat = camera_world_to_view_matrix(orb.camera);
                mat4_t proj_mat = camera_view_to_clip_matrix_persp(orb.camera, aspect_ratio);
                mat4_t inv_proj_mat = camera_clip_to_view_matrix_persp(orb.camera, aspect_ratio);

                gbuffer_clear(&gbuf);

                const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
                    GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_TRANSPARENCY };

                glEnable(GL_CULL_FACE);
                glCullFace(GL_BACK);

                glEnable(GL_DEPTH_TEST);
                glDepthFunc(GL_LESS);
                glDepthMask(GL_TRUE);
                glEnable(GL_SCISSOR_TEST);

                glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
                glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
                glViewport(0, 0, gbuf.width, gbuf.height);
                glScissor(0, 0, gbuf.width, gbuf.height);

                md_gl_draw_op_t draw_op = {};
                draw_op.type = MD_GL_REP_BALL_AND_STICK;
                draw_op.args.ball_and_stick.ball_scale   = 1.0f;
                draw_op.args.ball_and_stick.stick_radius = 1.0f;
                draw_op.rep = gl_rep;

                md_gl_draw_args_t draw_args = {
                    .shaders = state.gl.shaders,
                    .draw_operations = {
                        .count = 1,
                        .ops = &draw_op
                    },
                    .view_transform = {
                        .view_matrix = (const float*)view_mat.elem,
                        .proj_matrix = (const float*)proj_mat.elem,
                    },
                    .picking_offset = {
                       .atom_base = state.picking_range_atom.beg,
                       .bond_base = state.picking_range_bond.beg,
                    },
                };

                md_gl_draw(&draw_args);

                glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
                glClearColor(1, 1, 1, 0);
                glClear(GL_COLOR_BUFFER_BIT);

                PUSH_GPU_SECTION("Postprocessing")
                postprocess_pipeline::Settings postprocess_settings = {};
                postprocess_pipeline::Inputs postprocess_inputs = {};

                postprocess_settings.background_color = {24.f, 24.f, 24.f};
                postprocess_settings.tonemap.enabled = state.visuals.tonemapping.enabled;
                postprocess_settings.tonemap.mode = state.visuals.tonemapping.tonemapper;
                postprocess_settings.tonemap.exposure = state.visuals.tonemapping.exposure;
                postprocess_settings.tonemap.gamma = state.visuals.tonemapping.gamma;
                postprocess_settings.ssao.enabled = false;
                postprocess_settings.dof.enabled = false;
                postprocess_settings.fxaa.enabled = true;
                postprocess_settings.taa.enabled = false;
                postprocess_settings.sharpen.enabled = false;

                postprocess_inputs.depth = orb.gbuf.tex.depth;
                postprocess_inputs.color = orb.gbuf.tex.color;
                postprocess_inputs.normal = orb.gbuf.tex.normal;
                postprocess_inputs.velocity = orb.gbuf.tex.velocity;

                ViewParam view_param = {
                    .matrix = {
                        .curr = {
                        .view = view_mat,
                        .proj = proj_mat,
                        .norm = view_mat,
                    },
                    .inv = {
                        .proj = inv_proj_mat,
                    }
                    },
                    .clip_planes = {
                        .near = orb.camera.near_plane,
                        .far  = orb.camera.far_plane,
                    },
                    .resolution = {orb_win_sz.x, orb_win_sz.y},
                    .fov_y = orb.camera.fov_y,
                };

                postprocess_pipeline::execute(postprocess_inputs, postprocess_settings, view_param);
                POP_GPU_SECTION()

                glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
                glDrawBuffer(GL_BACK);
                glDisable(GL_SCISSOR_TEST);

                PUSH_GPU_SECTION("ORB GRID RAYCAST")
                    for (int i = 0; i < num_mos; ++i) {
                        volume::RenderDesc vol_desc = {
                            .render_target = {
                                .depth  = orb.gbuf.tex.depth,
                                .color  = orb.iso_tex[i],
                                .width  = orb.gbuf.width,
                                .height = orb.gbuf.height,
                                .clear_color = true,
                            },
                            .texture = {
                                .density_volume = orb.vol[i].tex_id,
                            },
                            .matrix = {
                                .model = orb.vol[i].texture_to_world,
                                .view  = view_mat,
                                .proj  = proj_mat,
                                .inv_proj = inv_proj_mat,
                            },
                            .iso = {
                                .enabled = true,
                                .count  = (size_t)orb.iso.count,
                                .values = orb.iso.values,
                                .colors = orb.iso.colors,
                            },
                            .shading = {
                                .env_radiance = state.visuals.background.color * state.visuals.background.intensity * 0.25f,
                                .roughness = 0.3f,
                                .dir_radiance = {10,10,10},
                                .ior = 1.5f,
                            },
                            .voxel_spacing = orb.vol[i].voxel_size,
                        };
                        volume::render_volume(vol_desc);
                    }
                POP_GPU_SECTION();
            }
        }
        ImGui::End();
    }

    void draw_export_window(ApplicationState& state) {
        if (!export_state.show_window) return;
        if (ImGui::Begin("VeloxChem Export", &export_state.show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::BeginCombo("Source", electronic_structure_source_str[(int)export_state.source])) {
                for (int i = 0; i < (int)ElectronicStructureSource::Count; i++) {
                    ElectronicStructureSource source = (ElectronicStructureSource)i;
                    bool is_selected = (export_state.source == source);
                    bool disabled = source == ElectronicStructureSource::DensityProperty ||
                                    ((source == ElectronicStructureSource::NaturalTransitionOrbital ||
                                     source == ElectronicStructureSource::TransitionDensity) &&
                                    (md_vlx_rsp_number_of_excited_states(vlx) == 0));

                    if (disabled) ImGui::PushDisabled();
                    if (ImGui::Selectable(electronic_structure_source_str[i], is_selected)) {
                        export_state.source = source;
                        export_state.use_magnitude = false;
                        if (export_state.source == ElectronicStructureSource::ElectronDensity) {
                            export_state.spin = ElectronicStructureSpin::Total;
                        }
                    }
                    if (disabled) ImGui::PopDisabled();

                    if (is_selected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }

            if (export_state.source == ElectronicStructureSource::MolecularOrbital ||
                export_state.source == ElectronicStructureSource::NaturalTransitionOrbital ||
                (export_state.source == ElectronicStructureSource::ElectronDensity && export_state.spin == ElectronicStructureSpin::Difference)) {
                const char* magnitude_label = "Magnitude";
                ImGui::Checkbox(magnitude_label, &export_state.use_magnitude);
            }

            if (export_state.source == ElectronicStructureSource::NaturalTransitionOrbital) {
                ImGui::Combo("Component", (int*)&export_state.nto_component, electronic_structure_nto_component_str, IM_ARRAYSIZE(electronic_structure_nto_component_str));
            }

            if (export_state.source == ElectronicStructureSource::TransitionDensity) {
                ImGui::Combo("Component", (int*)&export_state.transition_density_component, electronic_structure_transition_density_component_str, IM_ARRAYSIZE(electronic_structure_transition_density_component_str));
            }

            if (export_state.source == ElectronicStructureSource::ElectronDensity) {
                const bool has_beta_density = md_vlx_scf_density_matrix_data(vlx, MD_VLX_SPIN_BETA) != nullptr;
                if (has_beta_density) {
                    const ElectronicStructureSpin spin_options[] = {
                        ElectronicStructureSpin::Total,
                        ElectronicStructureSpin::Alpha,
                        ElectronicStructureSpin::Beta,
                        ElectronicStructureSpin::Difference,
                    };
                    const char* spin_labels[] = { "Total", "Alpha", "Beta", "Difference" };
                    int spin_idx = 0;
                    for (int i = 0; i < (int)IM_ARRAYSIZE(spin_options); ++i) {
                        if (export_state.spin == spin_options[i]) {
                            spin_idx = i;
                            break;
                        }
                    }
                    if (export_state.spin != spin_options[spin_idx]) {
                        export_state.spin = spin_options[spin_idx];
                    }
                    if (ImGui::Combo("Spin", &spin_idx, spin_labels, IM_ARRAYSIZE(spin_labels))) {
                        export_state.spin = spin_options[spin_idx];
                        if (export_state.spin != ElectronicStructureSpin::Difference) {
                            export_state.use_magnitude = false;
                        }
                    }
                } else {
                    export_state.spin = ElectronicStructureSpin::Total;
                    export_state.use_magnitude = false;
                }
            }

            int orb_idx = 0;
            bool show_mo_combo = false;
            bool show_lambda_combo = false;
            bool show_excited_state_combo = false;

            switch (export_state.source) {
            case ElectronicStructureSource::MolecularOrbital:
                show_mo_combo = true;
                break;
            case ElectronicStructureSource::NaturalTransitionOrbital:
                show_excited_state_combo = true;
                show_lambda_combo = true;
                break;
            case ElectronicStructureSource::TransitionDensity:
                show_excited_state_combo = true;
                break;
            case ElectronicStructureSource::ElectronDensity:
                break;
            default:
                break;
            }

            if (show_mo_combo) {
                md_vlx_scf_type_t type = md_vlx_scf_type(vlx);
                if (type == MD_VLX_SCF_RESTRICTED_OPENSHELL) {
                    const char* options[2] = {"Alpha", "Beta"};
                    if (ImGui::BeginCombo("##MO_TYPE", options[export_state.mo.type])) {
                        for (size_t i = 0; i < ARRAY_SIZE(options); ++i) {
                            if (ImGui::Selectable(options[i])) {
                                export_state.mo.type = (md_vlx_spin_t)i;
                            }
                        }
                        ImGui::EndCombo();
                    }
                }

                int orb_homo_idx = homo_idx[export_state.mo.type];
                int orb_lumo_idx = lumo_idx[export_state.mo.type];

                auto write_lbl = [orb_homo_idx, orb_lumo_idx](char* buf, size_t cap, int idx) {
                    const char* suffix = "";
                    if (idx == orb_homo_idx) {
                        suffix = "(HOMO)";
                    } else if (idx == orb_lumo_idx) {
                        suffix = "(LUMO)";
                    }
                    snprintf(buf, cap, "%i %s", idx + 1, suffix);
                };

                char lbl[32];
                write_lbl(lbl, sizeof(lbl), export_state.mo.idx);
                int num_orbs = (int)num_molecular_orbitals();
                if (ImGui::BeginCombo("MO Idx", lbl)) {
                    for (int i = 0; i < num_orbs; i++) {                        
                        const bool is_selected = (export_state.mo.idx == i);

                        write_lbl(lbl, sizeof(lbl), i);
                        if (ImGui::Selectable(lbl, is_selected)) {
                            export_state.mo.idx = i;
                        }

                        if (is_selected) {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }
                orb_idx = export_state.mo.idx;
            }

            if (show_excited_state_combo) {
                int num_orbs = (int)num_natural_transition_orbitals();
                char lbl[32];
                snprintf(lbl, sizeof(lbl), "%i", export_state.nto.idx + 1);
                if (ImGui::BeginCombo("Excited State Idx", lbl)) {
                    for (int i = 0; i < num_orbs; i++) {
                        const bool is_selected = (export_state.nto.idx == i);

                        snprintf(lbl, sizeof(lbl), "%i", i + 1);
                        if (ImGui::Selectable(lbl, is_selected)) {
                            export_state.nto.idx = i;
                        }

                        if (is_selected) {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }
                orb_idx = export_state.nto.idx;
            }

            if (show_lambda_combo) {
                double lambdas[4] = {0};
                size_t lambda_count = md_vlx_rsp_nto_lambdas_extract(lambdas, vlx, (size_t)export_state.nto.idx, ARRAY_SIZE(lambdas));
                if (lambda_count > 0) {
                    if ((size_t)export_state.nto.lambda_idx >= lambda_count) {
                        export_state.nto.lambda_idx = (int)lambda_count - 1;
                    }
                    char lbl[32];
                    int i = export_state.nto.lambda_idx;
                    snprintf(lbl, sizeof(lbl), "%i (%.4f)", i + 1, lambdas[i]);
                    if (ImGui::BeginCombo("Lambda Idx", lbl)) {
                        for (i = 0; i < (int)lambda_count; i++) {
                            const bool is_selected = (export_state.nto.lambda_idx == i);
                            snprintf(lbl, sizeof(lbl), "%i (%.4f)", i + 1, lambdas[i]);
                            if (ImGui::Selectable(lbl, is_selected)) {
                                export_state.nto.lambda_idx = i;
                            }

                            if (is_selected) {
                                ImGui::SetItemDefaultFocus();
                            }
                        }
                        ImGui::EndCombo();
                    }
                }
            }

            vec3_t extent;
            if (export_state.use_obb) {
                extent = oabb.max_ext - oabb.min_ext;
            } else {
                extent = aabb.max_ext - aabb.min_ext;
            }

            int dim[3];
            compute_dim(dim, extent, volume_resolution_samples_per_angstrom[(int)export_state.resolution] * BOHR_TO_ANGSTROM);

            char lbl[32];
            static const char* fmt = (const char*)u8"%s (%i×%i×%i)";
            snprintf(lbl, sizeof(lbl), fmt, volume_resolution_str[(int)export_state.resolution], dim[0], dim[1], dim[2]);
            if (ImGui::BeginCombo("Volume Resolution (XYZ)", lbl)) {
                for (int i = 0; i < IM_ARRAYSIZE(volume_resolution_str); ++i) {
                    const bool is_selected = ((int)export_state.resolution == i);
                    compute_dim(dim, extent, volume_resolution_samples_per_angstrom[i] * BOHR_TO_ANGSTROM);
                    snprintf(lbl, sizeof(lbl), fmt, volume_resolution_str[i], dim[0], dim[1], dim[2]);
                    if (ImGui::Selectable(lbl, is_selected)) {
                        export_state.resolution = (VolumeResolution)i;
                    }

                    if (is_selected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }

            enum {
                EXPORT_FILE_FORMAT_CUBE = 0,
                EXPORT_FILE_FORMAT_RAW_MHD_XYZ,
                EXPORT_FILE_FORMAT_COUNT,
            };
            const char* file_format_str[] = {"cube", "raw + mhd + xyz"};
            str_t file_format_ext[] = {STR_LIT("cube"), STR_LIT("raw")};

            STATIC_ASSERT(ARRAY_SIZE(file_format_str) == EXPORT_FILE_FORMAT_COUNT);
            STATIC_ASSERT(ARRAY_SIZE(file_format_ext) == EXPORT_FILE_FORMAT_COUNT);

            static int file_format = 0;
            ImGui::Combo("File format", &file_format, file_format_str, IM_ARRAYSIZE(file_format_str));

            ImGui::Checkbox("Use OBB", &export_state.use_obb);
            ImGui::SetItemTooltip("Use Oriented Bounding Box for volume (rotate and fit the volume to the data)");

            if (ImGui::Button("Export")) {
                str_t ext = file_format_ext[file_format];
                char path_buf[2048];
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, ext)) {
                    str_t path = {path_buf, strnlen(path_buf, sizeof(path_buf))};

                    md_file_t file = { 0 };
                    if (!md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                        MD_LOG_ERROR("Failed to open file for writing: '" STR_FMT "'", STR_ARG(path));
                        return;
                    }

                    const double samples_per_unit_length = volume_resolution_samples_per_angstrom[(int)export_state.resolution] * BOHR_TO_ANGSTROM;
                    md_grid_t grid = {};
                    if (export_state.use_obb) {
                        init_grid(&grid, oabb.orientation, oabb.min_ext, oabb.max_ext, samples_per_unit_length);
                    } else {
                        init_grid(&grid, mat3_ident(), aabb.min_ext, aabb.max_ext, samples_per_unit_length);
                    }

                    md_temp_scope_t temp = md_temp_begin();

                    defer {
                        md_file_close(&file); 
                        md_temp_end(temp);
                    };

                    Volume vol = {};
                    init_volume(&vol, grid, GL_R32F);

                    switch (export_state.source) {
                    case ElectronicStructureSource::MolecularOrbital:
                    {
                        md_gto_op_t op = es_gto_op(export_state.use_magnitude);
                        evaluate_mo(state, vol.tex_id, grid, export_state.mo.type, export_state.mo.idx, op);
                        break;
                    }
                    case ElectronicStructureSource::NaturalTransitionOrbital:
                    {
                        md_vlx_nto_type_t type = export_state.nto_component == ElectronicStructureNtoComponent::Particle ? MD_VLX_NTO_PARTICLE : MD_VLX_NTO_HOLE;
                        md_gto_op_t op = es_gto_op(export_state.use_magnitude);
                        
                        evaluate_nto(state, vol.tex_id, grid, export_state.nto.idx, export_state.nto.lambda_idx, type, op);
                        break;
                    }
                    case ElectronicStructureSource::TransitionDensity:
                    {
                        evaluate_transition_density(state, vol.tex_id, grid, export_state.nto.idx, export_state.transition_density_component);
                        break;
                    }
                    case ElectronicStructureSource::ElectronDensity:
                        evaluate_electron_density(state, vol.tex_id, grid, export_state.spin, es_gto_op(export_state.spin == ElectronicStructureSpin::Difference && export_state.use_magnitude));
                        break;
                    default:
                        MD_LOG_ERROR("Unsupported export type");
                        goto done;
                    }

                    int num_samples = vol.dim[0] * vol.dim[1] * vol.dim[2];
                    size_t              natoms = md_vlx_number_of_atoms(vlx);
                    const dvec3_t* vlx_coords  = md_vlx_atom_coordinates(vlx);
                    const uint8_t* vlx_numbers = md_vlx_atomic_numbers(vlx);
                    float* data = md_temp_alloc_array(temp, float, num_samples);

                    vec3_t origin = grid.origin;
                    vec3_t step_x = grid.orientation[0] * grid.spacing[0];
                    vec3_t step_y = grid.orientation[1] * grid.spacing[1];
                    vec3_t step_z = grid.orientation[2] * grid.spacing[2];

                    // Readbacks are queued, not completed, by the evaluate_* calls above.
                    // Exporting reads the texture directly, so it is one of the few places that
                    // genuinely needs the data now rather than next frame.
                    gpu_volume_jobs_drain(&state);

                    // Extract data from OpenGL Texture
                    glBindTexture(GL_TEXTURE_3D, vol.tex_id);
                    glGetTexImage(GL_TEXTURE_3D, 0, GL_RED, GL_FLOAT, data);
                    glBindTexture(GL_TEXTURE_3D, 0);

                    switch (file_format) {
                    case EXPORT_FILE_FORMAT_CUBE: {
                        md_file_printf(file, "Cube file generated by VIAMD\n");
                        if (export_state.source == ElectronicStructureSource::ElectronDensity) {
                            md_file_printf(file, "%s %s %s\n", electronic_structure_source_str[(int)export_state.source], electronic_structure_value_mode_str(export_state.source, export_state.spin, export_state.use_magnitude), electronic_structure_spin_str[(int)export_state.spin]);
                        } else {
                            md_file_printf(file, "%s %s\n", electronic_structure_source_str[(int)export_state.source], electronic_structure_value_mode_str(export_state.source, export_state.spin, export_state.use_magnitude));
                        }

                        md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", -(int)natoms, origin.x, origin.y, origin.z);
                        md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol.dim[0], step_x.x, step_x.y, step_x.z);
                        md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol.dim[1], step_y.x, step_y.y, step_y.z);
                        md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol.dim[2], step_z.x, step_z.y, step_z.z);

                        for (size_t i = 0; i < natoms; ++i) {
                            md_file_printf(file, "%5i %12.6f %12.6f %12.6f %12.6f\n", (int)vlx_numbers[i], (float)vlx_numbers[i], vlx_coords[i].x * ANGSTROM_TO_BOHR, vlx_coords[i].y * ANGSTROM_TO_BOHR, vlx_coords[i].z * ANGSTROM_TO_BOHR);
                        }

                        // This entry somehow relates to the number of densities
                        md_file_printf(file, "%5i %5i\n", 1, orb_idx + 1);

                        // Write density data
                        int count = 0;
                        for (int x = 0; x < vol.dim[0]; ++x) {
                            for (int y = 0; y < vol.dim[1]; ++y) {
                                for (int z = 0; z < vol.dim[2]; ++z) {
                                    int idx = z * vol.dim[0] * vol.dim[1] + y * vol.dim[0] + x;
                                    float val = data[idx];
                                    md_file_printf(file, " %12.6E", val);
                                    if (++count % 6 == 0) md_file_printf(file, "\n");
                                }
                            }
                        }
                        MD_LOG_INFO("Successfully exported electronic structure to '" STR_FMT "'", STR_ARG(path));
                        break;
                    }
                    case EXPORT_FILE_FORMAT_RAW_MHD_XYZ: {
                        /*
                        * # EXAMPLE
                        ObjectType = Image
                        NDims = 3
                        DimSize = 256 256 64
                        ElementType = MET_USHORT
                        HeaderSize = -1
                        ElementSpacing = 1 1 1
                        ElementByteOrderMSB = False
                        ElementDataFile = image.raw
                        */

                        str_t basepath;
                        if (!extract_file_path_without_ext(&basepath, path)) {
                            MD_LOG_ERROR("Failed to extract base path");
                            goto done;
                        }

                        str_t filename_raw;
                        if (!extract_file(&filename_raw, path)) {
                            MD_LOG_ERROR("Failed to extract filename raw");
                            goto done;
                        }

                        str_t filename;
                        if (!extract_file_path_without_ext(&filename, filename_raw)) {
                            MD_LOG_ERROR("Failed to extract filename");
                            goto done;
                        }

                        md_file_write(file, data, num_samples * sizeof(float));

                        md_allocator_i* temp_arena = md_temp_allocator(temp);
                        str_t mhd_path = str_printf(temp_arena, STR_FMT ".mhd", STR_ARG(basepath));
                        str_t xyz_path = str_printf(temp_arena, STR_FMT ".xyz", STR_ARG(basepath));

                        md_file_t mhd_file = { 0 };
                        if (!md_file_open(&mhd_file, mhd_path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                            MD_LOG_ERROR("Failed to open .mhd file");
                            goto done;
                        }

                        vec3_t pos   = origin * BOHR_TO_ANGSTROM;
                        vec3_t step  = vol.voxel_size * BOHR_TO_ANGSTROM;

                        md_file_printf(mhd_file, "ObjectType = Image\n");
                        md_file_printf(mhd_file, "NDims = 3\n");
                        md_file_printf(mhd_file, "DimSize = %i %i %i\n", vol.dim[0], vol.dim[1], vol.dim[2]);
                        md_file_printf(mhd_file, "ElementType = MET_FLOAT\n");
                        md_file_printf(mhd_file, "HeaderSize = -1\n");
                        md_file_printf(mhd_file, "ElementSpacing = %12.6f %12.6f %12.6f\n", step.x, step.y, step.z);
                        md_file_printf(mhd_file, "ElementByteOrderMSB = False\n");
                        md_file_printf(mhd_file, "ElementDataFile = " STR_FMT "\n", STR_ARG(filename_raw));
                        md_file_printf(mhd_file, "Position = %12.6f %12.6f %12.6f\n", pos.x, pos.y, pos.z);
                        md_file_printf(mhd_file, "Orientation =\n");
                        md_file_printf(mhd_file, "%12.6f %12.6f %12.6f\n", grid.orientation[0].x, grid.orientation[0].y, grid.orientation[0].z);
                        md_file_printf(mhd_file, "%12.6f %12.6f %12.6f\n", grid.orientation[1].x, grid.orientation[1].y, grid.orientation[1].z);
                        md_file_printf(mhd_file, "%12.6f %12.6f %12.6f\n", grid.orientation[2].x, grid.orientation[2].y, grid.orientation[2].z);

                        md_file_close(&mhd_file);

                        md_file_t xyz_file = { 0 };
                        if (!md_file_open(&xyz_file, xyz_path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                            MD_LOG_ERROR("Failed to open .xyz file");
                            goto done;
                        }
                        defer { md_file_close(&xyz_file); };
                        
                        // XYZ
                        md_file_printf(xyz_file, "%zu\n", natoms);
                        md_file_printf(xyz_file, "Geometry extracted from VeloxChem dataset: '%s'\n", state.files.molecule);
                        for (size_t i = 0; i < natoms; ++i) {
                            str_t sym = md_util_element_symbol(vlx_numbers[i]);
                            md_file_printf(xyz_file, "%-2s %12.6f %12.6f %12.6f\n", sym.ptr, vlx_coords[i].x, vlx_coords[i].y, vlx_coords[i].z);
                        }

                        MD_LOG_INFO("Successfully exported electronic structure to '" STR_FMT "'", STR_ARG(path));

                        break;
                    }
                    default:
                        ASSERT(false);
                        break;
                    }
                }
            }
        }
        done:
        ImGui::End();
    }

    // Calculates the transition matrix heuristic
    static inline void compute_transition_matrix(float* out_matrix, const size_t num_groups, const float* hole_charges, const float* part_charges) {
        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        if (num_groups == 0) return;

        int* donors         = md_temp_alloc_zero_array(temp, int, num_groups);
        int* acceptors      = md_temp_alloc_zero_array(temp, int, num_groups);
        double* charge_diff = md_temp_alloc_zero_array(temp, double, num_groups);

        int num_donors = 0;
        int num_acceptors = 0;

        // Sanitize inputs: clamp small negative numerical noise to zero before normalization.
        double hole_sum = 0;
        double part_sum = 0;

        for (size_t i = 0; i < num_groups; i++) {
            const double h = MAX(0.0, (double)hole_charges[i]);
            const double p = MAX(0.0, (double)part_charges[i]);
            hole_sum += h;
            part_sum += p;
        }

        // Guard against degenerate inputs (no charge to attribute).
        if (hole_sum <= 0.0 || part_sum <= 0.0) {
            MEMSET(out_matrix, 0, sizeof(float) * num_groups * num_groups);
            return;
        }

        for (size_t i = 0; i < num_groups; i++) {
            // percentages (clamped, normalized so each sums to 1)
            double gsCharge = MAX(0.0, (double)hole_charges[i]) / hole_sum;
            double esCharge = MAX(0.0, (double)part_charges[i]) / part_sum;

            if (gsCharge > esCharge) {
                donors[num_donors++] = (int)i;
            } else {
                acceptors[num_acceptors++] = (int)i;
            }

            out_matrix[i * num_groups + i] = (float)MIN(gsCharge, esCharge);
            charge_diff[i] = esCharge - gsCharge;
        }

        double total_acceptor_charge = 0;
        for (int i = 0; i < num_acceptors; i++) {
            total_acceptor_charge += MAX(0.0, charge_diff[acceptors[i]]);
        }

        if (total_acceptor_charge <= 0.0) {
            return;
        }

        for (int i = 0; i < num_donors; i++) {
            int donor = donors[i];
            double charge_deficit = MAX(0.0, -charge_diff[donor]);
            for (int j = 0; j < num_acceptors; j++) {
                int acceptor = acceptors[j];
                double acc = MAX(0.0, charge_diff[acceptor]);
                double contrib = charge_deficit * acc / total_acceptor_charge;
                if (contrib < 0.0) contrib = 0.0;
                out_matrix[acceptor * num_groups + donor] = (float)contrib;
            }
        }
    }

    void delete_group(int index) {
        vec4_t deleted_color = nto.group.color[index];
        for (size_t i = 0; i < md_array_size(nto.atom_group_idx); i++) {
            if (nto.atom_group_idx[i] == (uint32_t)index) {
                nto.atom_group_idx[i] = 0;
            }
            else if (nto.atom_group_idx[i] > (uint32_t)index) {
                nto.atom_group_idx[i] -= 1;
            }
        }
        for (int i = index; i < (int)nto.group.count - 1; ++i) {
            nto.group.color[i] = nto.group.color[i + 1];
            MEMCPY(nto.group.label[i], nto.group.label[i + 1], sizeof(nto.group.label[i]));
        }
        nto.group.count--;
        snprintf(nto.group.label[nto.group.count], sizeof(nto.group.label[nto.group.count]), "Group %i", (int)nto.group.count);
        nto.group.color[nto.group.count] = deleted_color;
    }

    void copy_group(int from_idx, int to_idx) {
        nto.group.color[to_idx] = nto.group.color[from_idx];
        MEMCPY(nto.group.label[to_idx], nto.group.label[from_idx], sizeof(nto.group.label[from_idx]));
        for (size_t i = 0; i < md_array_size(nto.atom_group_idx); i++) {
            if ((int)nto.atom_group_idx[i] == to_idx) {
                nto.atom_group_idx[i] = (uint32_t)from_idx;
            }
        }
    }

    void swap_groups(int a, int b) {
        vec4_t color = nto.group.color[a];
        nto.group.color[a] = nto.group.color[b];
        nto.group.color[b] = color;

        char label_buf[sizeof(nto.group.label[a])];
        //sprintf(label_buf, nto.group.label[row_n]);
        MEMCPY(label_buf, nto.group.label[a], sizeof(nto.group.label[a])); //Current to buf
        MEMCPY(nto.group.label[a], nto.group.label[b], sizeof(nto.group.label[b])); //Next to current
        MEMCPY(nto.group.label[b], label_buf, sizeof(label_buf)); //Buf to next

        for (size_t i = 0; i < md_array_size(nto.atom_group_idx); i++) {
            if ((int)nto.atom_group_idx[i] == a) {
                nto.atom_group_idx[i] = (uint32_t)b;
            }
            else if ((int)nto.atom_group_idx[i] == b) {
                nto.atom_group_idx[i] = (uint32_t)a;
            }
        }
    }

    static inline void imgui_delayed_hover_tooltip(const char* text) {
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayNormal)) {
            ImGui::SetTooltip("%s", text);
        }
    }

    void draw_nto_window(ApplicationState& state) {
        if (!nto.show_window) return;

        size_t num_excited_states = md_vlx_rsp_number_of_excited_states(vlx);

        if (num_excited_states == 0) return;

        bool open_context_menu = false;
        bool export_transition_diagram_requested = false;
        static bool edit_mode = false;

        ImGui::SetNextWindowSize(ImVec2(500, 300), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Transition Analysis", &nto.show_window, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
            nto.group.hovered_index = -1;
            bool viewport_hovered = false;

            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("File")) {
                    export_transition_diagram_requested = ImGui::MenuItem("Export Transition Diagram");
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Settings")) {
                    ImGui::Text("Iso Colors");
                    ImGui::Checkbox("Link Attachment/Detachment Density Colors", &nto.link_attachment_detachment_density);
                    if (nto.link_attachment_detachment_density) {
                        ImGui::ColorEdit4("##Color Density", nto.col_den.elem);
                        ImGui::SetItemTooltip("Color Density");
                    } else {
                        ImGui::ColorEdit4("##Color Attachment", nto.col_att.elem);
                        ImGui::SetItemTooltip("Color Attachment");
                        ImGui::ColorEdit4("##Color Detachment", nto.col_det.elem);
                        ImGui::SetItemTooltip("Color Detachment");
                    }

                    ImGui::ColorEdit4("##Color Positive", nto.col_pos.elem);
                    ImGui::SetItemTooltip("Color Positive");
                    ImGui::ColorEdit4("##Color Negative", nto.col_neg.elem);
                    ImGui::SetItemTooltip("Color Negative");

                    ImGui::Text("Transition Dipole Moments");
                    ImGui::Checkbox("Enabled", &nto.dipole.enabled);
                    ImGui::SliderFloat("Scale", &nto.dipole.vector_scale, 1.0f, 10.0f);
                    ImGui::Checkbox("Show Angle", &nto.dipole.show_angle);
                    ImGui::ColorEdit4("Color Electric", nto.dipole.colors[0].elem);
                    ImGui::SetItemTooltip("Color Electric");
                    ImGui::ColorEdit4("Color Magnetic", nto.dipole.colors[1].elem);
                    ImGui::SetItemTooltip("Color Magnetic");
                    ImGui::EndMenu();
                }
                ImGui::EndMenuBar();
            }

            const ImVec2 outer_size = {300.f + edit_mode * 80, 0.f};
            ImGui::PushItemWidth(outer_size.x);
            ImGui::BeginGroup();

            ImGui::Text("Transition State Index");
            ImGui::Spacing();
            if (ImGui::BeginListBox("##NTO Index", outer_size)) {
                if (ImGui::IsWindowHovered()) {
                    rsp.hovered = -1;
                }
                for (int i = 0; i < (int)num_excited_states; ++i) {
                    bool is_selected = rsp.selected == i;
                    bool is_hovered  = rsp.hovered  == i;
                    char buf[32];
                    snprintf(buf, sizeof(buf), "%i", i + 1);
                    if (is_hovered) {
                        ImGui::PushStyleColor(ImGuiCol_Header, ImGui::GetColorU32(ImGuiCol_HeaderHovered));
                    }
                    if (ImGui::Selectable(buf, is_selected || is_hovered)) {
                        rsp.selected = i;
                    }
                    if (is_hovered) {
                        ImGui::PopStyleColor();
                    }
                    if (ImGui::IsItemHovered()) {
                        rsp.hovered = i;
                    }
                }
                ImGui::EndListBox();
            }
            ImGui::Spacing();
            const double iso_min = 1.0e-8;
            const double iso_max = 5.0;
            
            ImGui::Spacing();
            ImGui::TextUnformatted(electronic_structure_iso_value_label()); 
            ImGui::SliderScalar("##Iso Value", ImGuiDataType_Double, &nto.iso_val, &iso_min, &iso_max, "%.6f", ImGuiSliderFlags_Logarithmic);
            ImGui::SetItemTooltip("Visual surface threshold. Attachment and detachment density views render at isovalue^2 for perceptual consistency with NTO amplitudes.");

            ImGui::Spacing();

            static const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersH | ImGuiTableFlags_SizingFixedFit;
            static const ImGuiTableColumnFlags columns_base_flags = ImGuiTableColumnFlags_NoSort;
            
            static bool hide_overlap_text = true;
            ImVec2 button_size(ImGui::GetFontSize() + ImGui::GetStyle().FramePadding.x, 0.f);

            ImGui::Checkbox("Edit mode", &edit_mode);
            ImGui::SameLine();
            ImGui::Checkbox("Hide overlapping text", &hide_overlap_text);

            int group_counts[MAX_NTO_GROUPS] = {0};
            for (size_t i = 0; i < md_array_size(nto.atom_group_idx); ++i) {
                int group_idx = nto.atom_group_idx[i] < nto.group.count ? nto.atom_group_idx[i] : 0;
                group_counts[group_idx] += 1;
            }

            if (ImGui::BeginTable("Group Table", 3 + edit_mode * 4, flags, outer_size, 0)) {
                ImGui::TableSetupColumn("Group", columns_base_flags, 150.f);
                ImGui::TableSetupColumn("Count", columns_base_flags, 0.0f);
                ImGui::TableSetupColumn("Color", columns_base_flags, 0.0f);
                if (edit_mode) {
                    ImGui::TableSetupColumn("##Delete", columns_base_flags, button_size.x);
                    ImGui::TableSetupColumn("##Move down", columns_base_flags, button_size.x);
                    ImGui::TableSetupColumn("##Move up", columns_base_flags, button_size.x);
                    ImGui::TableSetupColumn("##Merge up", columns_base_flags, button_size.x);
                }

                //ImGui::TableSetupColumn("Coord Y", columns_base_flags, 0.0f);
                //ImGui::TableSetupColumn("Coord Z", columns_base_flags | ImGuiTableColumnFlags_WidthFixed, 0.0f);
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableHeadersRow();

                ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);
                ImGui::PushStyleColor(ImGuiCol_Header, IM_BLUE);
                bool show_unassigned = group_counts[0] > 0;
                int row_n = show_unassigned ? 0 : 1;
                for (; row_n < (int)nto.group.count; row_n++) {
                    ImGui::PushID((int)row_n);
                    ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                    ImGui::TableNextRow(ImGuiTableRowFlags_None, 0);

                    ImGui::TableNextColumn();
                    ImGui::AlignTextToFramePadding();
                    if (edit_mode) {
                        ImGui::PushItemWidth(150.f);
                        if (row_n > 0) {
                            ImGui::InputText("##label", nto.group.label[row_n], sizeof(nto.group.label[row_n]));
                        } else {
                            ImGui::Text("%s", nto.group.label[0]);
                        }
                        ImGui::PopItemWidth();
                    }
                    else {
                        ImGui::Selectable(nto.group.label[row_n], false, selectable_flags);
                        if (ImGui::TableGetHoveredRow() == (int)(row_n + show_unassigned)) {
                            nto.group.hovered_index = (int8_t)row_n;
                            md_bitfield_clear(&state.selection.highlight_mask);
                            for (size_t j = 0; j < md_array_size(nto.atom_group_idx); j++) {
                                if (nto.atom_group_idx[j] == (uint32_t)row_n) {
                                    //Add j to highlight mask
                                    md_bitfield_set_bit(&state.selection.highlight_mask, j); //TODO: Hover only works if you hold leftMouseBtn
                                }
                            }

                            //Selection
                            if (ImGui::IsKeyDown(ImGuiKey_MouseLeft) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                                md_bitfield_or_inplace(&state.selection.selection_mask, &state.selection.highlight_mask);
                            }
                            //Deselect
                            else if (ImGui::IsKeyDown(ImGuiKey_MouseRight) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                                md_bitfield_andnot_inplace(&state.selection.selection_mask, &state.selection.highlight_mask);
                            }
                        }
                    }


                    ImGui::TableNextColumn();
                    ImGui::Text("%i", group_counts[row_n]);

                    ImGui::TableNextColumn();

                    //Center the color picker
                    ImVec2 cell_size = ImGui::GetContentRegionAvail();
                    float color_size = ImGui::GetFrameHeight();
                    ImVec2 padding((cell_size.x - color_size) * 0.5f, 0.0f);
                    ImGui::SetCursorPos(ImGui::GetCursorPos() + padding);

                    if (edit_mode && row_n > 0) { //You cannot edit the Unassigned color
                        ImGui::ColorEdit4Minimal("##color", nto.group.color[row_n].elem);
                    } else {
                        ImGui::ColorButton("##color", vec_cast(nto.group.color[row_n]));
                    }

                    if (edit_mode) {
                        if (row_n > 0) {
                            ImGui::TableNextColumn();
                            ImGui::SetCursorPosX(ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x - button_size.x);
                            if (ImGui::DeleteButton(ICON_FA_XMARK, button_size)) {
                                delete_group(row_n);
                            }
                            imgui_delayed_hover_tooltip("Delete this group (related atoms will be unassigned)");

                            //Down Arrow
                            ImGui::TableNextColumn();
                            if (row_n == (int)nto.group.count - 1) {
                                ImGui::Dummy(button_size);
                            }
                            else {
                                if (ImGui::IsKeyDown(ImGuiKey_LeftCtrl)) {
                                    if (ImGui::Button(ICON_FA_ANGLES_DOWN, button_size)) {
                                        for (int i = row_n; i < (int)nto.group.count - 1; i++) {
                                            swap_groups(i, i + 1);
                                        }
                                    }
                                    imgui_delayed_hover_tooltip("Move to bottom");

                                }
                                else {
                                    if (ImGui::Button(ICON_FA_ANGLE_DOWN, button_size)) {
                                        swap_groups(row_n, row_n + 1);
                                    }
                                    imgui_delayed_hover_tooltip("Move down (Ctrl-click to move to bottom)");
                                }
                            }

                            //Up Arrow
                            ImGui::TableNextColumn();
                            if (row_n == 1) {
                                ImGui::Dummy(button_size);
                            }
                            else {
                                if (ImGui::IsKeyDown(ImGuiKey_LeftCtrl)) {
                                    if (ImGui::Button(ICON_FA_ANGLES_UP, button_size)) {
                                        int last = (int)nto.group.count - 1;
                                        for (int i = last; i >= (int)row_n - 1; i--) {
                                            swap_groups(i, i - 1);
                                        }
                                    }
                                    imgui_delayed_hover_tooltip("Move to top");
                                }
                                else {
                                    if (ImGui::Button(ICON_FA_ANGLE_UP, button_size)) {
                                        swap_groups(row_n, row_n - 1);
                                    }
                                    imgui_delayed_hover_tooltip("Move up (Ctrl-click to move to top)");
                                }
                            }

                            //Merge with above arrow
                            ImGui::TableNextColumn();
                            if (row_n == 1) {
                                ImGui::Dummy(button_size);
                            }
                            else {
                                if (ImGui::Button(ICON_FA_TRASH_ARROW_UP, button_size)) {
                                    for (size_t i = 0; i < md_array_size(nto.atom_group_idx); i++) {
                                        if (nto.atom_group_idx[i] == (uint32_t)row_n) {
                                            nto.atom_group_idx[i] = (uint32_t)(row_n - 1);
                                        }
                                    }
                                    delete_group(row_n);
                                }
                                imgui_delayed_hover_tooltip("Merge with above");
                            }
                        }
                    }
                    ImGui::PopID();
                }

                ImGui::PopStyleColor(2);
                ImGui::EndTable();
            }

            if (ImGui::Button("Add new group")) {
                nto.group.count++;
            }

            ImGui::EndGroup();
            ImGui::PopItemWidth();

            ImGui::SameLine();

            // Animate camera towards targets
            const double dt = state.app.timing.delta_s;
            camera_animate(&nto.camera, nto.target, dt);

            ImVec2 canvas_sz = ImGui::GetContentRegionAvail();   // Resize canvas to what's available
            canvas_sz.x = MAX(canvas_sz.x, 50.0f);
            canvas_sz.y = MAX(canvas_sz.y, 50.0f);

            if (export_transition_diagram_requested) {
                ImVec2 export_sz = { MAX(canvas_sz.x * 0.5f, 1.0f), MAX(canvas_sz.y, 1.0f) };
                export_sz *= 2.0f;
                export_transition_diagram(export_sz, &nto, hide_overlap_text);
            }

            ImGui::Dummy(canvas_sz);

            // Draw border and background color
            ImGuiIO& io = ImGui::GetIO();

            ImVec2 canvas_p0 = ImGui::GetItemRectMin();
            ImVec2 canvas_p1 = ImGui::GetItemRectMax();

            bool reset_view = false;

            if (rsp.selected != -1) {
                if (nto.sel_nto_idx != rsp.selected) {
                    nto.sel_nto_idx = rsp.selected;
                    size_t nto_idx = (size_t)rsp.selected;

                    init_volume(&nto.vol[NTO_Attachment], nto.grid, GL_R32F);
                    init_volume(&nto.vol[NTO_Detachment], nto.grid, GL_R32F);

                    evaluate_transition_density(state, nto.vol[NTO_Attachment].tex_id, nto.grid, nto_idx, ElectronicStructureTransitionDensityComponent::Attachment);
                    evaluate_transition_density(state, nto.vol[NTO_Detachment].tex_id, nto.grid, nto_idx, ElectronicStructureTransitionDensityComponent::Detachment);
                }
            }

            const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();

            ImVec2 grid_p0 = canvas_p0;
            ImVec2 grid_p1 = canvas_p0 + canvas_sz * ImVec2(0.5f, 1.0f);
            ImVec2 win_sz = (grid_p1 - grid_p0) / ImVec2(1.0f, 2.0f);
            win_sz.x = floorf(win_sz.x);
            win_sz.y = floorf(win_sz.y);

            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->ChannelsSplit(2);
            draw_list->ChannelsSetCurrent(0);
            draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(255, 255, 255, 255));

            const float aspect_ratio = win_sz.x / win_sz.y;

            mat4_t view_mat     = camera_world_to_view_matrix(nto.camera);
            mat4_t proj_mat     = camera_view_to_clip_matrix_persp(nto.camera, aspect_ratio);
            mat4_t inv_view_mat = camera_view_to_world_matrix(nto.camera);
            mat4_t inv_proj_mat = camera_clip_to_view_matrix_persp(nto.camera, aspect_ratio);

            mat4_t world_to_clip = proj_mat * view_mat;
            mat4_t clip_to_world = inv_view_mat * inv_proj_mat;

            if (rsp.selected != -1) {
                const int nto_target_idx[2] = {
                    NTO_Attachment,
                    NTO_Detachment,
                };

                // Draw Attachment / Detachment
                for (int i = 0; i < 2; ++i) {
					ImGui::PushID(i);
					defer { ImGui::PopID(); };

                    int idx = nto_target_idx[i];

                    ImVec2 p0 = grid_p0 + win_sz * ImVec2(0.0f, (float)(i+0));
                    ImVec2 p1 = grid_p0 + win_sz * ImVec2(1.0f, (float)(i+1));
                    ImRect rect = {p0, p1};
                    ImGui::SetCursorScreenPos(p0);

                    draw_list->PushClipRect(p0, p1);
                    draw_list->ChannelsSetCurrent(1);
                    ImGui::PushID(i);
                    defer {
                        draw_list->PopClipRect();
                        ImGui::PopID();
                    };

                    InteractionSurfaceState surface = interaction_surface(interaction_surface_nto, vec_cast(rect.GetSize()));

                    InteractionSurfaceViewTransformArgs view_args = {
                        .camera = nto.camera,
                        .trackball_param = state.view.trackball_param,
                    };
                    InteractionSurfaceViewTransformResult view_result = interaction_surface_view_transform_apply(&nto.target, surface, view_args);
                    if (view_result.reset_requested) {
                        nto.target = default_view;
                    }

                    if (surface.hovered) {
                        InteractionSurfaceHitArgs hit_args = {
                            .picking_surface = &nto.picking_surface,
                            .picking_handler = state.picking_handler,
                            .fbo = nto.gbuf.fbo,
                            .width = nto.gbuf.width,
                            .height = nto.gbuf.height,
                            .clip_to_world = clip_to_world,
                        };

                        PickingHit hit = {};
                        interaction_surface_hit_extract(&hit, surface, hit_args);

                        InteractionSurfaceEvent event = {};
                        interaction_surface_event_extract(&event, surface, hit);

                        event.clip_to_world = clip_to_world;
                        event.world_to_clip = world_to_clip;

                        if (event.kind == InteractionSurfaceEventKind::ContextMenu) {
                            open_context_menu = true;
                        } else if (event.kind == InteractionSurfaceEventKind::RegionSelect) {
                            const md_bitfield_t* candidate_mask = nullptr;
                            if (event.selection_mode == InteractionSelectionMode::Remove) {
                                // When removing, only consider currently selected atoms as candidates for region selection
                                candidate_mask = &state.selection.selection_mask;
                            }
                            point_set_region_mask_compute(&state.selection.highlight_mask,
                                state.mold.state.x, state.mold.state.y, state.mold.state.z, state.mold.state.num_atoms,
                                candidate_mask, event.world_to_clip, event.region_min, event.region_max, event.surface_size);

                            grow_mask_by_selection_granularity(&state.selection.highlight_mask, state.selection.granularity, state.mold.sys);
                            if (event.region_phase == InteractionSurfaceEventPhase::Commit) {
                                // Merge highlight into selection
                                if (event.selection_mode == InteractionSelectionMode::Append) {
                                    md_bitfield_or_inplace(&state.selection.selection_mask, &state.selection.highlight_mask);
                                }
                                else if (event.selection_mode == InteractionSelectionMode::Remove) {
                                    md_bitfield_andnot_inplace(&state.selection.selection_mask, &state.selection.highlight_mask);
                                }
                                md_bitfield_clear(&state.selection.highlight_mask);
                            }
                        } else {
                            viamd::event_system_broadcast_event(viamd::EventType_ViamdInteractionSurface, viamd::EventPayloadType_InteractionSurfaceEvent, &event);
                        }
                    }

                    draw_list->ChannelsSetCurrent(0);

                    //ImVec2 text_pos_bl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p1.y - TEXT_BASE_HEIGHT);
                    ImVec2 text_pos_tl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p0.y + TEXT_BASE_HEIGHT * 0.5f);
                    const char* lbl = ((i & 1) == 0) ? "Attachment" : "Detachment";

                    //char buf[32];
                    //snprintf(buf, sizeof(buf), (const char*)u8"λ: %.3f", nto_lambda[i / 2]);
                    draw_list->AddImage((ImTextureID)(intptr_t)nto.gbuf.tex.transparency, p0, p1, { 0,1 }, { 1,0 });
                    draw_list->AddImage((ImTextureID)(intptr_t)nto.iso_tex[idx], p0, p1, { 0,1 }, { 1,0 });
                    //draw_list->AddText(text_pos_bl, ImColor(0,0,0), buf);
                    draw_list->AddText(text_pos_tl, ImColor(0,0,0), lbl);

                    if (nto.dipole.enabled) {
                        const vec3_t mid = vec3_lerp(aabb.min_ext, aabb.max_ext, 0.5f);

                        const dvec3_t* magnetic_dp = md_vlx_rsp_magnetic_transition_dipole_moments(vlx);
                        const dvec3_t* electric_dp = md_vlx_rsp_electric_transition_dipole_moments(vlx);

                        vec3_t magn_vec = {(float)magnetic_dp[rsp.selected].x, (float)magnetic_dp[rsp.selected].y, (float)magnetic_dp[rsp.selected].z};
                        vec3_t elec_vec = {(float)electric_dp[rsp.selected].x, (float)electric_dp[rsp.selected].y, (float)electric_dp[rsp.selected].z};

                        magn_vec *= (float)(nto.dipole.vector_scale * BOHR_TO_ANGSTROM);
                        elec_vec *= (float)(nto.dipole.vector_scale * BOHR_TO_ANGSTROM);

                        auto proj_point = [world_to_clip, win_sz, p0](vec3_t point) -> ImVec2 {
                            const vec4_t p = mat4_mul_vec4(world_to_clip, vec4_from_vec3(point, 1.0));
                            const ImVec2 c = {
                                p0.x + ( p.x / p.w * 0.5f + 0.5f) * win_sz.x,
                                p0.y + (-p.y / p.w * 0.5f + 0.5f) * win_sz.y,
                            };
                            return c;
                        };

                        auto draw_arrow = [draw_list](ImVec2 beg, ImVec2 end, ImU32 color, float thickness) {
                            float len2 = ImLengthSqr(end-beg);
                            if (len2 > 1.0e-3f) {
                                float len = sqrtf(len2);
                                ImVec2 d = (end - beg) / len;
                                ImVec2 o = {-d.y, d.x};
                                float d_scl = ImMin(thickness * 3, len);
                                float o_scl = 0.57735026918962576451f * d_scl;
                                draw_list->AddLine(beg, end - d * thickness, color, thickness);
                                draw_list->AddTriangleFilled(end, end - d * d_scl + o * o_scl, end - d * d_scl - o * o_scl, color);
                            }
                        };

                        const ImU32 elec_color = convert_color(nto.dipole.colors[0]);
                        const ImU32 magn_color = convert_color(nto.dipole.colors[1]);
                        const ImU32 text_color = IM_COL32_BLACK;
                        const float arrow_thickness = 3.0f;

                        ImVec2 c      = proj_point(mid);
                        ImVec2 c_elec = proj_point(mid + elec_vec);
                        ImVec2 c_magn = proj_point(mid + magn_vec);

                        draw_arrow(c, c_elec, elec_color, arrow_thickness);
                        draw_list->AddText(c_elec, text_color, (const char*)u8"μe");

                        draw_arrow(c, c_magn, magn_color, arrow_thickness);
                        draw_list->AddText(c_magn, text_color, (const char*)u8"μm");

                        if (nto.dipole.show_angle) {
                            const vec3_t magn_dir = vec3_normalize(magn_vec);
                            const vec3_t elec_dir = vec3_normalize(elec_vec);
                            const float angle = acosf(vec3_dot(magn_dir, elec_dir));
                            const ImU32 angle_color = IM_COL32(0, 0, 0, 128);
                            const float angle_thickness = 1.0f;
                            const float angle_vec_scale = 0.1f;

                            const size_t num_seg = 20;
                            const vec3_t axis = vec3_normalize(vec3_cross(elec_dir, magn_dir));

                            for (size_t seg = 0; seg <= num_seg; ++seg) {
                                float  t = (float)(seg) / (float)num_seg;
                                vec3_t v = quat_mul_vec3(quat_axis_angle(axis, angle * t), elec_dir);
                                ImVec2 p = proj_point(mid + v * angle_vec_scale);
                                draw_list->PathLineTo(p);
                            }
                            draw_list->PathStroke(angle_color, 0, angle_thickness);

                            char buf[32];
                            snprintf(buf, sizeof(buf), (const char*)u8"%.2f°", RAD_TO_DEG(angle));
                            draw_list->AddText(c, text_color, buf);
                        }
                    }
                }
                // Draw Sankey Diagram of Transition Matrix
                {
                    //compute_transition_matrix(nto.transition_matrix, nto.group.count, nto.transition_density_hole, nto.transition_density_part);
                    ImVec2 p0 = canvas_p0 + canvas_sz * ImVec2(0.5f, 0.0f);
                    ImVec2 p1 = canvas_p1;
                    im_sankey_diagram(&state, {p0.x, p0.y, p1.x, p1.y}, &nto, hide_overlap_text);
                    ImVec2 text_pos_bl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p1.y - TEXT_BASE_HEIGHT);
                    draw_list->AddText(text_pos_bl, ImColor(0, 0, 0, 255), "Transition Diagram");
                }
                // Draw grid
                {
                    ImVec2 p0 = {floorf(canvas_p0.x + canvas_sz.x * 0.5f), canvas_p0.y};
                    ImVec2 p1 = {floorf(canvas_p0.x + canvas_sz.x * 0.5f), canvas_p1.y};
                    draw_list->AddLine(p0, p1, IM_COL32(0, 0, 0, 255));
                }

                int num_subwindows = 2;
                for (int i = 1; i < num_subwindows; ++i) {
                    float y = floorf(canvas_p0.y + canvas_sz.y / ((float)num_subwindows) * i);
                    float x0 = canvas_p0.x;
                    float x1 = floorf(canvas_p0.x + canvas_sz.x * (i & 1 ? 0.5f : 1.0f));
                    draw_list->AddLine({x0, y}, {x1, y}, IM_COL32(0, 0, 0, 255));
                }

                // Draw stuff

                const bool is_hovered = ImGui::IsItemHovered();
                const ImVec2 origin(canvas_p0.x, canvas_p0.y);  // Lock scrolled origin
                const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

                int width  = MAX(1, (int)(win_sz.x * io.DisplayFramebufferScale.x));
                int height = MAX(1, (int)(win_sz.y * io.DisplayFramebufferScale.y));

                auto& gbuf = nto.gbuf;
                if ((int)gbuf.width != width || (int)gbuf.height != height) {
                    gbuffer_init(&gbuf, width, height);
                    for (int i : nto_target_idx) {
                        gl::init_texture_2D(nto.iso_tex + i, width, height, GL_RGBA8);
                    }
                }

                if (is_hovered) {
                    if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                        reset_view = true;
                    }
                }

                if (reset_view) {
					vec3_t center = vec3_lerp(aabb.min_ext, aabb.max_ext, 0.5f) * BOHR_TO_ANGSTROM;
					vec3_t half_ext = (aabb.max_ext - aabb.min_ext) * 0.5f * BOHR_TO_ANGSTROM;
                    nto.target = compute_optimal_view(center, half_ext, oabb.orientation, nto.distance_scale);
                }

                if (nto.show_coordinate_system_widget) {
                    float ext = MIN(win_sz.x, win_sz.y) * 0.3f;
                    float pad = 20.0f;
                    ImVec2 size = ImVec2(ext, ext);
                    ImGui::SetCursorScreenPos(canvas_p0 + ImVec2(pad, 2.0f * win_sz.y - ext - pad));

                    quat_t out_orientation = nto.target.orientation;
                    if (ImGui::CoordinateSystemWidget(&out_orientation, nto.camera.orientation, size)) {
                        const vec3_t look_at = camera_get_look_at(nto.target);
                        nto.target.orientation = quat_normalize(out_orientation);
                        nto.target.position = camera_position_from_look_at(look_at, nto.target.orientation, nto.target.distance);
                    }
                }

                gbuffer_clear(&gbuf);

                const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
                    GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_TRANSPARENCY };

                glEnable(GL_CULL_FACE);
                glCullFace(GL_BACK);

                glEnable(GL_DEPTH_TEST);
                glDepthFunc(GL_LESS);
                glDepthMask(GL_TRUE);
                glEnable(GL_SCISSOR_TEST);

                glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
                glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
                glViewport(0, 0, gbuf.width, gbuf.height);
                glScissor(0, 0, gbuf.width, gbuf.height);

                auto draw_rep = [&state](md_gl_rep_t& rep, md_gl_shaders_t& shaders, mat4_t& view_mat, mat4_t& proj_mat, uint32_t atom_mask = 0) {
                    md_gl_draw_op_t draw_op = {};
                    draw_op.type = MD_GL_REP_BALL_AND_STICK;
                    draw_op.args.ball_and_stick.ball_scale   = 1.0f;
                    draw_op.args.ball_and_stick.stick_radius = 1.0f;
                    draw_op.rep = rep;

                    md_gl_draw_args_t draw_args = {
                        .shaders = shaders,
                        .draw_operations = {
                            .count = 1,
                            .ops = &draw_op
                        },
                        .view_transform = {
                            .view_matrix = (const float*)view_mat.elem,
                            .proj_matrix = (const float*)proj_mat.elem,
                        },

                        .picking_offset = {
                            .atom_base = state.picking_range_atom.beg,
                            .bond_base = state.picking_range_bond.beg,
                        },

                        .atom_mask = atom_mask,
                    };

                    md_gl_draw(&draw_args);
                };

                draw_rep(nto.gl_rep, state.gl.shaders, view_mat, proj_mat);

                glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
                glClearColor(1, 1, 1, 0);
                glClear(GL_COLOR_BUFFER_BIT);

                if (true) {
                    PUSH_GPU_SECTION("Selection")
                    const bool atom_selection_empty = md_bitfield_popcount(&state.selection.selection_mask) == 0;
                    const bool atom_highlight_empty = md_bitfield_popcount(&state.selection.highlight_mask) == 0;

                    glDepthMask(0);

                    // @NOTE(Robin): This is a b*tch to get right, What we want is to separate in a single pass, the visible selected from the
                    // non visible selected. In order to achieve this, we start with a cleared stencil of value 1 then either set it to zero selected and not visible
                    // and to two if it is selected and visible. But the visible atoms should always be able to write over a non visible 0, but not the other way around.
                    // Hence the GL_GREATER stencil test against the reference value of 2.

                    if (!atom_selection_empty) {
                        glColorMask(0, 0, 0, 0);

                        glEnable(GL_DEPTH_TEST);
                        glDepthFunc(GL_EQUAL);

                        glEnable(GL_STENCIL_TEST);
                        glStencilMask(0xFF);

                        glClearStencil(1);
                        glClear(GL_STENCIL_BUFFER_BIT);

                        glStencilFunc(GL_GREATER, 0x02, 0xFF);
                        glStencilOp(GL_KEEP, GL_ZERO, GL_REPLACE);
                        draw_rep(nto.gl_rep, state.gl.shaders_lean_and_mean, view_mat, proj_mat, AtomBit_Selected);

                        glDisable(GL_DEPTH_TEST);

                        glStencilMask(0x0);
                        glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
                        glColorMask(1, 1, 1, 1);

                        glStencilFunc(GL_EQUAL, 2, 0xFF);
                        postprocessing::blit_color(state.selection.color.selection.visible);

                        glStencilFunc(GL_EQUAL, 0, 0xFF);
                        postprocessing::blit_color(state.selection.color.selection.hidden);
                    }

                    if (!atom_highlight_empty) {
                        glColorMask(0, 0, 0, 0);

                        glEnable(GL_DEPTH_TEST);
                        glDepthFunc(GL_EQUAL);

                        glEnable(GL_STENCIL_TEST);
                        glStencilMask(0xFF);

                        glClearStencil(1);
                        glClear(GL_STENCIL_BUFFER_BIT);

                        glStencilFunc(GL_GREATER, 0x02, 0xFF);
                        glStencilOp(GL_KEEP, GL_ZERO, GL_REPLACE);
                        draw_rep(nto.gl_rep, state.gl.shaders_lean_and_mean, view_mat, proj_mat, AtomBit_Highlighted);

                        glDisable(GL_DEPTH_TEST);

                        glStencilMask(0x0);
                        glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
                        glColorMask(1, 1, 1, 1);

                        glStencilFunc(GL_EQUAL, 2, 0xFF);
                        postprocessing::blit_color(state.selection.color.highlight.visible);

                        glStencilFunc(GL_EQUAL, 0, 0xFF);
                        postprocessing::blit_color(state.selection.color.highlight.hidden);
                    }

                    glDisable(GL_STENCIL_TEST);

                    /*
                    if (!atom_selection_empty) {
                    PUSH_GPU_SECTION("Desaturate") {
                    glColorMask(1, 1, 1, 1);
                    glDrawBuffer(GL_COLOR_ATTACHMENT_COLOR);
                    postprocessing::scale_hsv(nto.gbuf.tex.color, vec3_t{1, state.selection.color.saturation, 1});
                    } POP_GPU_SECTION()
                    }
                    */

                    glDepthFunc(GL_LESS);
                    glDepthMask(0);
                    glColorMask(1,1,1,1);
                    POP_GPU_SECTION()
                }

                PUSH_GPU_SECTION("Postprocessing")
                postprocess_pipeline::Settings postprocess_settings = {};
                postprocess_pipeline::Inputs postprocess_inputs = {};

                postprocess_settings.background_color = state.visuals.background.color * state.visuals.background.intensity;
                postprocess_settings.tonemap.enabled = state.visuals.tonemapping.enabled;
                postprocess_settings.tonemap.mode = state.visuals.tonemapping.tonemapper;
                postprocess_settings.tonemap.exposure = state.visuals.tonemapping.exposure;
                postprocess_settings.tonemap.gamma = state.visuals.tonemapping.gamma;
                postprocess_settings.ssao.enabled = false;
                postprocess_settings.dof.enabled = false;
                postprocess_settings.fxaa.enabled = true;
                postprocess_settings.taa.enabled = false;
                postprocess_settings.sharpen.enabled = false;

                postprocess_inputs.depth = nto.gbuf.tex.depth;
                postprocess_inputs.color = nto.gbuf.tex.color;
                postprocess_inputs.normal = nto.gbuf.tex.normal;
                postprocess_inputs.velocity = nto.gbuf.tex.velocity;
                postprocess_inputs.transparency = nto.gbuf.tex.transparency;

                ViewParam view_param = {
                    .matrix = {
                        .curr = {
                        .view = view_mat,
                        .proj = proj_mat,
                        .norm = view_mat,
                    },
                    .inv = {
                        .proj = inv_proj_mat,
                    }
                    },
                    .clip_planes = {
                        .near = nto.camera.near_plane,
                        .far  = nto.camera.far_plane,
                    },
                    .resolution = {win_sz.x, win_sz.y},
                    .fov_y = nto.camera.fov_y,
                };

                postprocess_pipeline::execute(postprocess_inputs, postprocess_settings, view_param);
                POP_GPU_SECTION()

                PUSH_GPU_SECTION("NTO RAYCAST")
                for (int i : nto_target_idx) {
                    bool is_density = (i == NTO_Attachment || i == NTO_Detachment);
                    
                    bool enabled = true;
                    size_t count = is_density ? 1 : 2;
                    float  values[2];
                    vec4_t colors[2];
                    
                    if (is_density) {
                        if (nto.link_attachment_detachment_density) {
                            colors[0] = nto.col_den;
                        } else {
                            if (i == NTO_Attachment) {
                                colors[0] = nto.col_att;
                            } else {
                                colors[0] = nto.col_det;
                            }
                        }
                        values[0] = (float)(nto.iso_val * nto.iso_val);
                    } else {
                        colors[0] = nto.col_pos;
                        colors[1] = nto.col_neg;
                        values[0] =  (float)nto.iso_val;
                        values[1] = -(float)nto.iso_val;
                    }

                    if (enabled) {
                        volume::RenderDesc vol_desc = {
                            .render_target = {
                                .depth  = nto.gbuf.tex.depth,
                                .color  = nto.iso_tex[i],
                                .width  = nto.gbuf.width,
                                .height = nto.gbuf.height,
                                .clear_color = true,
                            },
                            .texture = {
                                .density_volume = nto.vol[i].tex_id,
                            },
                            .matrix = {
                                .model = nto.vol[i].texture_to_world,
                                .view  = view_mat,
                                .proj  = proj_mat,
                                .inv_proj = inv_proj_mat,
                            },
                            .iso = {
                                .enabled = enabled,
                                .count   = count,
                                .values  = values,
                                .colors  = colors,
                            },
                            .shading = {
                                .env_radiance = state.visuals.background.color * state.visuals.background.intensity * 0.25f,
                                .roughness = 0.3f,
                                .dir_radiance = {10,10,10},
                                .ior = 1.5f,
                            },
                            .voxel_spacing = nto.vol[i].voxel_size,
                        };
                        volume::render_volume(vol_desc);
                    }
                }
                POP_GPU_SECTION();

                // Reset state
                glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
                glDrawBuffer(GL_BACK);
                glDisable(GL_SCISSOR_TEST);
            }
            if (ImGui::IsWindowHovered() && nto.group.hovered_index == -1 && !viewport_hovered) {
                md_bitfield_clear(&state.selection.highlight_mask);
            }
        }
        ImGui::End();

        if (open_context_menu) {
            ImGui::OpenPopup("Context Menu");
        }

        if (ImGui::BeginPopup("Context Menu")) {
			size_t count = md_bitfield_popcount(&state.selection.selection_mask);
			if (count == 0) {
                ImGui::Text("No atoms selected");
            } else {
                if (ImGui::BeginMenu("Assign to group")) {
                    for (size_t i = 0; i < nto.group.count; i++) {
                        if (ImGui::MenuItem(nto.group.label[i])) {
                            for (size_t j = 0; j < md_array_size(nto.atom_group_idx); j++) {
                                //Is the atom selected?
                                if (md_bitfield_test_bit(&state.selection.selection_mask, j)) {
                                    //If so, set its group to the selected group index
                                    nto.atom_group_idx[j] = (uint32_t)i;
                                }
                            }
                        }
                    }
                    ImGui::EndMenu();
                }
                if (nto.group.count < MAX_NTO_GROUPS) {
                    if (ImGui::MenuItem("Assign to new group")) {
                        //Create a new group and add an item to it
                        uint32_t group_idx = (uint32_t)(nto.group.count++);
                        for (size_t i = 0; i < md_array_size(nto.atom_group_idx); i++) {
                            if (md_bitfield_test_bit(&state.selection.selection_mask, i)) {
                                //If so, set its group to the selected group index
                                nto.atom_group_idx[i] = group_idx;
                            }
                        }
                    }
                }
            }
            ImGui::EndPopup();
        }

    }
    // Everything derived from the group assignment: the representation colours and the transition
    // matrix. Lifted out of draw_nto_window so it runs whether or not that window is open - the
    // Charge Transfer window needs the same numbers, and a window should not have to be visible
    // for another one to have data.
    //
    // Driven by rsp.selected rather than nto.sel_nto_idx: the latter tracks which state's VOLUMES
    // are resident, which is a rendering concern and is only updated while the NTO window draws.
    // ---------------------------------------------------------------------------
    // Workspace persistence
    // ---------------------------------------------------------------------------

    void ensure_pending_group_storage() {
        if (pending_groups.atoms_initialized) return;
        for (size_t i = 0; i < MAX_NTO_GROUPS; ++i) {
            md_bitfield_init(&pending_groups.atoms[i], md_get_heap_allocator());
        }
        pending_groups.atoms_initialized = true;
    }

    void serialize_workspace(viamd::serialization_state_t& ser) {
        if (!vlx || !nto.atom_group_idx) return;

        viamd::write_section_header(ser, STR_LIT("VeloxChem"));

        viamd::write_bool(ser, STR_LIT("ShowTransitionAnalysis"), nto.show_window);
        viamd::write_bool(ser, STR_LIT("ShowChargeTransfer"),     flow.show_window);
        viamd::write_int (ser, STR_LIT("FlowMode"),               (int64_t)flow.mode);
        viamd::write_flt (ser, STR_LIT("FlowThreshold"),          flow.cut.threshold);
        viamd::write_flt (ser, STR_LIT("FlowOthersCap"),          flow.cut.other_max);
        viamd::write_flt (ser, STR_LIT("FlowSizeNormalize"),      flow.size_normalize);
        viamd::write_int (ser, STR_LIT("FlowOthersExpanded"),     (int64_t)flow.cut.other_expanded);
        viamd::write_flt (ser, STR_LIT("FlowNodeThickness"),      flow.params.node_thickness);
        viamd::write_flt (ser, STR_LIT("FlowNodeGap"),            flow.params.node_gap);
        viamd::write_flt (ser, STR_LIT("FlowRibbonAlpha"),        flow.style.ribbon_alpha);
        viamd::write_flt (ser, STR_LIT("FlowRibbonCurvature"),    flow.style.ribbon_curvature);
        viamd::write_int (ser, STR_LIT("FlowOrder"),              (int64_t)flow.params.order);
        viamd::write_int (ser, STR_LIT("FlowSweeps"),             (int64_t)flow.params.barycentre_sweeps);
        viamd::write_flt (ser, STR_LIT("FlowRibbonMinPx"),        flow.style.ribbon_min_px);
        viamd::write_flt (ser, STR_LIT("FlowNodeFillTint"),       flow.style.node_fill_tint);
        viamd::write_flt (ser, STR_LIT("FlowBorderWidth"),        flow.style.border_width);
        viamd::write_bool(ser, STR_LIT("FlowColorLabels"),        flow.style.color_labels);

        viamd::write_int(ser, STR_LIT("GroupCount"), (int64_t)nto.group.count);

        // One bitfield per group rather than one index per atom: the same information, but it
        // survives being looked at in a text editor.
        md_bitfield_t mask = md_bitfield_create(md_get_heap_allocator());
        defer { md_bitfield_free(&mask); };

        for (size_t g = 0; g < nto.group.count; ++g) {
            md_bitfield_clear(&mask);
            for (size_t i = 0; i < md_array_size(nto.atom_group_idx); ++i) {
                if (nto.atom_group_idx[i] == (uint32_t)g) {
                    md_bitfield_set_bit(&mask, i);
                }
            }
            viamd::write_str(ser, STR_LIT("GroupLabel"), str_from_cstr(nto.group.label[g]));
            viamd::write_vec4(ser, STR_LIT("GroupColor"), nto.group.color[g]);
            viamd::write_bitfield(ser, STR_LIT("GroupAtoms"), &mask);
        }
    }

    void deserialize_workspace(viamd::deserialization_state_t& deser) {
        ensure_pending_group_storage();

        pending_groups.pending = true;
        pending_groups.count = 0;
        for (size_t i = 0; i < MAX_NTO_GROUPS; ++i) {
            md_bitfield_clear(&pending_groups.atoms[i]);
            pending_groups.label[i][0] = '\0';
            pending_groups.color[i] = vec4_t{0.5f, 0.5f, 0.5f, 1.0f};
        }

        // 'GroupLabel' opens a group; the color and atom entries that follow belong to it. The
        // index only advances on a label, so a malformed section cannot walk off the end.
        int64_t group_idx = -1;

        str_t ident = {};
        str_t arg   = {};
        while (viamd::next_entry(ident, arg, deser)) {
            if (str_eq_cstr(ident, "ShowTransitionAnalysis")) {
                viamd::extract_bool(nto.show_window, arg);
            } else if (str_eq_cstr(ident, "ShowChargeTransfer")) {
                viamd::extract_bool(flow.show_window, arg);
            } else if (str_eq_cstr(ident, "FlowMode")) {
                int mode = FLOW_MODE_NTO;
                if (viamd::extract_int(mode, arg)) {
                    flow.mode = (flow_mode_t)CLAMP(mode, 0, FLOW_MODE_COUNT - 1);
                }
            } else if (str_eq_cstr(ident, "FlowThreshold")) {
                viamd::extract_flt(flow.cut.threshold, arg);
            } else if (str_eq_cstr(ident, "FlowOthersCap")) {
                viamd::extract_flt(flow.cut.other_max, arg);
            } else if (str_eq_cstr(ident, "FlowSizeNormalize")) {
                viamd::extract_flt(flow.size_normalize, arg);
            } else if (str_eq_cstr(ident, "FlowOthersExpanded")) {
                int mask = 0;
                if (viamd::extract_int(mask, arg)) flow.cut.other_expanded = (uint32_t)MAX(0, mask);
            } else if (str_eq_cstr(ident, "FlowNodeThickness")) {
                viamd::extract_flt(flow.params.node_thickness, arg);
            } else if (str_eq_cstr(ident, "FlowNodeGap")) {
                viamd::extract_flt(flow.params.node_gap, arg);
            } else if (str_eq_cstr(ident, "FlowRibbonAlpha")) {
                viamd::extract_flt(flow.style.ribbon_alpha, arg);
            } else if (str_eq_cstr(ident, "FlowRibbonCurvature")) {
                viamd::extract_flt(flow.style.ribbon_curvature, arg);
            } else if (str_eq_cstr(ident, "FlowOrder")) {
                int order = MD_FLOW_ORDER_WEIGHT;
                if (viamd::extract_int(order, arg)) {
                    flow.params.order = (md_flow_order_t)CLAMP(order, 0, 1);
                }
            } else if (str_eq_cstr(ident, "FlowRibbonMinPx")) {
                viamd::extract_flt(flow.style.ribbon_min_px, arg);
            } else if (str_eq_cstr(ident, "FlowNodeFillTint")) {
                viamd::extract_flt(flow.style.node_fill_tint, arg);
            } else if (str_eq_cstr(ident, "FlowBorderWidth")) {
                viamd::extract_flt(flow.style.border_width, arg);
            } else if (str_eq_cstr(ident, "FlowColorLabels")) {
                viamd::extract_bool(flow.style.color_labels, arg);
            } else if (str_eq_cstr(ident, "FlowSweeps")) {
                int sweeps = 2;
                if (viamd::extract_int(sweeps, arg)) {
                    flow.params.barycentre_sweeps = (uint32_t)CLAMP(sweeps, 0, 8);
                }
            } else if (str_eq_cstr(ident, "GroupCount")) {
                int count = 0;
                if (viamd::extract_int(count, arg)) {
                    pending_groups.count = (size_t)CLAMP(count, 0, MAX_NTO_GROUPS);
                }
            } else if (str_eq_cstr(ident, "GroupLabel")) {
                group_idx += 1;
                if (group_idx < MAX_NTO_GROUPS) {
                    viamd::extract_to_char_buf(pending_groups.label[group_idx], sizeof(pending_groups.label[0]), arg);
                }
            } else if (str_eq_cstr(ident, "GroupColor")) {
                if (0 <= group_idx && group_idx < MAX_NTO_GROUPS) {
                    viamd::extract_vec4(pending_groups.color[group_idx], arg);
                }
            } else if (str_eq_cstr(ident, "GroupAtoms")) {
                if (0 <= group_idx && group_idx < MAX_NTO_GROUPS) {
                    viamd::extract_bitfield(&pending_groups.atoms[group_idx], arg);
                }
            }
        }

        // Trust the labels actually present over a GroupCount that disagrees with them.
        pending_groups.count = MAX(pending_groups.count, (size_t)(group_idx + 1));
        pending_groups.count = MIN(pending_groups.count, (size_t)MAX_NTO_GROUPS);
    }

    // Consume-once: a buffered grouping belongs to the workspace that was just opened, and must
    // not be stamped onto whatever file is loaded next.
    void apply_pending_groups() {
        if (!pending_groups.pending) return;
        pending_groups.pending = false;

        if (!nto.atom_group_idx || pending_groups.count == 0) return;

        const size_t num_atoms = md_array_size(nto.atom_group_idx);
        MEMSET(nto.atom_group_idx, 0, sizeof(uint32_t) * num_atoms);

        nto.group.count = pending_groups.count;
        for (size_t g = 0; g < pending_groups.count; ++g) {
            if (pending_groups.label[g][0]) {
                snprintf(nto.group.label[g], sizeof(nto.group.label[g]), "%s", pending_groups.label[g]);
            }
            nto.group.color[g] = pending_groups.color[g];

            for (size_t i = 0; i < num_atoms; ++i) {
                if (md_bitfield_test_bit(&pending_groups.atoms[g], i)) {
                    nto.atom_group_idx[i] = (uint32_t)g;
                }
            }
        }

        update_nto_group_colors();
        MD_LOG_INFO("Restored %zu VeloxChem groups from the workspace", pending_groups.count);
    }

    // ---------------------------------------------------------------------------
    // Charge Transfer Analysis
    // ---------------------------------------------------------------------------

    // Node keys are (column, level, index) so that expansion, selection and hover survive a
    // rebuild. Anything keyed by array index would reset every time the state or the grouping
    // changes, which is exactly when the user least wants to lose what they had open.
    enum flow_level_t {
        FLOW_LEVEL_GROUP = 0,
        FLOW_LEVEL_ATOM  = 1,
        FLOW_LEVEL_MO    = 2,   // a canonical molecular orbital; owns no atoms
    };

    // Whether a node stands for a set of atoms. Orbital nodes - NTOs and MOs alike - do not, so
    // they take no part in atom selection, region sweeps or group assignment.
    static inline bool flow_level_has_atoms(uint32_t level) {
        return level == FLOW_LEVEL_GROUP || level == FLOW_LEVEL_ATOM;
    }

    static inline uint64_t flow_node_key(uint32_t column, uint32_t level, uint32_t index) {
        return ((uint64_t)column << 40) | ((uint64_t)level << 32) | (uint64_t)index;
    }
    static inline uint32_t flow_key_index(uint64_t key) { return (uint32_t)(key & 0xFFFFFFFFu); }
    static inline uint32_t flow_key_level(uint64_t key) { return (uint32_t)((key >> 32) & 0xFFu); }

    static const size_t FLOW_ATOM_LABEL_STRIDE = 12;   // "Xx" + up to 8 digits + nul

    const char* flow_atom_label(uint32_t atom_idx) const {
        if (!flow.atom_labels || atom_idx >= flow.num_atom_labels) return "?";
        return flow.atom_labels + (size_t)atom_idx * FLOW_ATOM_LABEL_STRIDE;
    }

    // Colour says what an atom belongs to. Inside a group that is the group; outside one there is
    // no group to name, so it falls back to the atom's ROLE - donating or accepting.
    //
    // Not CPK here, deliberately. Now that colour lives on the border rather than the fill, a
    // white element (hydrogen) would draw a white outline on a near-white box and vanish, and a
    // near-black one (carbon) would say nothing about the transition. Role colour is both always
    // visible and the thing the diagram is actually about; the 3D viewport still uses CPK, because
    // there a molecule should look like a molecule.
    vec4_t flow_atom_color(uint32_t atom_idx, uint32_t column) const {
        const uint32_t group_idx = (atom_idx < md_array_size(nto.atom_group_idx)) ? nto.atom_group_idx[atom_idx] : 0;
        if (group_idx > 0 && group_idx < nto.group.count) {
            return nto.group.color[group_idx];
        }
        return (column == 0) ? vec4_t{0.78f, 0.24f, 0.22f, 1.0f}    // donating
                             : vec4_t{0.13f, 0.56f, 0.35f, 1.0f};   // accepting
    }

    char* flow_mo_label_storage(uint32_t mo_idx) {
        if (!flow.mo_labels || mo_idx >= flow.num_mo_labels) return nullptr;
        return flow.mo_labels + (size_t)mo_idx * FLOW_ATOM_LABEL_STRIDE;
    }

    void flow_build_label_storage() {
        const size_t num_mo = md_vlx_scf_number_of_molecular_orbitals(vlx);
        if (num_mo > 0 && !flow.mo_labels) {
            flow.mo_labels = (char*)md_alloc(arena, num_mo * FLOW_ATOM_LABEL_STRIDE);
            MEMSET(flow.mo_labels, 0, num_mo * FLOW_ATOM_LABEL_STRIDE);
            flow.num_mo_labels = num_mo;
        }
        flow_build_atom_labels();
    }

    void flow_build_atom_labels() {
        const size_t num_atoms = md_vlx_number_of_atoms(vlx);
        if (num_atoms == 0 || flow.atom_labels) return;

        const md_element_t* z = md_vlx_atomic_numbers(vlx);
        flow.atom_labels = (char*)md_alloc(arena, num_atoms * FLOW_ATOM_LABEL_STRIDE);
        flow.num_atom_labels = num_atoms;

        for (size_t i = 0; i < num_atoms; ++i) {
            char* dst = flow.atom_labels + i * FLOW_ATOM_LABEL_STRIDE;
            const str_t sym = z ? md_util_element_symbol(z[i]) : str_from_cstr("X");
            // One-based, matching how the rest of VIAMD numbers atoms to the user.
            snprintf(dst, FLOW_ATOM_LABEL_STRIDE, "%.*s%zu", (int)sym.len, sym.ptr, i + 1);
        }
    }

    // The M1 view: two columns of groups, exactly what the old diagram showed, but through the
    // md_flow model. Kept so the NTO view below has something to be checked against.
    bool build_flow_graph_groups() {
        md_flow_graph_init(&flow.graph, 2, arena);

        if (!nto.transition_matrix || nto.transition_matrix_dim != nto.group.count) return false;

        const size_t num_groups = nto.group.count;
        if (num_groups == 0) return false;

        double hole_sum = 0.0;
        double part_sum = 0.0;
        for (size_t i = 0; i < num_groups; ++i) {
            hole_sum += MAX(0.0f, nto.transition_density_hole[i]);
            part_sum += MAX(0.0f, nto.transition_density_part[i]);
        }
        if (hole_sum <= 0.0 || part_sum <= 0.0) return false;

        uint32_t src_node[MAX_NTO_GROUPS];
        uint32_t dst_node[MAX_NTO_GROUPS];
        for (size_t i = 0; i < MAX_NTO_GROUPS; ++i) {
            src_node[i] = MD_FLOW_INVALID_INDEX;
            dst_node[i] = MD_FLOW_INVALID_INDEX;
        }

        // A group can donate without accepting (or the reverse), so the two columns hold
        // different node sets rather than one set drawn twice.
        for (size_t i = 0; i < num_groups; ++i) {
            const float w = (float)(MAX(0.0f, nto.transition_density_hole[i]) / hole_sum);
            if (w <= 0.0f) continue;
            md_flow_node_t node = {
                .column = 0,
                .parent = MD_FLOW_INVALID_INDEX,
                .level  = 0,
                .weight = w,
                .color  = nto.group.color[i],
                .label  = str_from_cstr(nto.group.label[i]),
                .key    = flow_node_key(0, FLOW_LEVEL_GROUP, (uint32_t)i),
            };
            src_node[i] = md_flow_graph_add_node(&flow.graph, &node);
        }
        for (size_t i = 0; i < num_groups; ++i) {
            const float w = (float)(MAX(0.0f, nto.transition_density_part[i]) / part_sum);
            if (w <= 0.0f) continue;
            md_flow_node_t node = {
                .column = 1,
                .parent = MD_FLOW_INVALID_INDEX,
                .level  = 0,
                .weight = w,
                .color  = nto.group.color[i],
                .label  = str_from_cstr(nto.group.label[i]),
                .key    = flow_node_key(1, FLOW_LEVEL_GROUP, (uint32_t)i),
            };
            dst_node[i] = md_flow_graph_add_node(&flow.graph, &node);
        }

        // compute_transition_matrix() already normalizes so that a donor's row sums to its hole
        // share and an acceptor's column to its particle share, which is exactly invariant (4).
        // If that ever stops being true the validation below says so instead of drawing nonsense.
        for (size_t src = 0; src < num_groups; ++src) {
            if (src_node[src] == MD_FLOW_INVALID_INDEX) continue;
            for (size_t dst = 0; dst < num_groups; ++dst) {
                if (dst_node[dst] == MD_FLOW_INVALID_INDEX) continue;
                const float w = nto.transition_matrix[dst * num_groups + src];
                if (w > 0.0f) {
                    md_flow_graph_add_link(&flow.graph, src_node[src], dst_node[dst], w);
                }
            }
        }

        return true;
    }

    // The M2 view: group -> NTO pair -> group.
    //
    // This is the first version that measures the group-to-group flow instead of guessing it.
    // compute_transition_matrix() only ever had the two marginals and split the difference
    // proportionally; routing through the NTOs gives
    //
    //     flow(g1 -> g2) = sum_k lambda_k * q_hole[k][g1] * q_part[k][g2]
    //
    // which is a real joint distribution. It is not the FULL one - the product form asserts that
    // hole and particle are independent once you know which NTO pair you are in, which is exactly
    // what the NTO decomposition claims and no more. The occupied x virtual amplitude matrix (M4)
    // drops that assumption too.
    bool build_flow_graph_nto() {
        md_flow_graph_init(&flow.graph, 3, arena);

        if (rsp.selected < 0 || !nto.atom_group_idx) return false;

        const size_t num_groups = nto.group.count;
        const size_t num_ao     = md_vlx_scf_number_of_atomic_orbitals(vlx);
        const double* S         = md_vlx_scf_overlap_matrix_data(vlx);
        const int* ao_to_atom   = md_vlx_ao_to_atom_idx(vlx);

        if (num_groups == 0 || num_ao == 0 || !S || !ao_to_atom) return false;

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        const size_t max_lambdas = MD_VLX_NTO_MAX_LAMBDAS;

        double* lambdas = md_temp_alloc_zero_array(temp, double, max_lambdas);
        double* C_hole  = md_temp_alloc_zero_array(temp, double, max_lambdas * num_ao);
        double* C_part  = md_temp_alloc_zero_array(temp, double, max_lambdas * num_ao);

        const size_t num_hole = md_vlx_rsp_nto_coefficients_extract(C_hole, lambdas, vlx, (size_t)rsp.selected, MD_VLX_NTO_HOLE, max_lambdas);
        const size_t num_part = md_vlx_rsp_nto_coefficients_extract(C_part, nullptr, vlx, (size_t)rsp.selected, MD_VLX_NTO_PARTICLE, max_lambdas);
        const size_t num_nto  = MIN(num_hole, num_part);
        if (num_nto == 0) return false;

        // Attribute to ATOMS, not to groups. The group numbers are then sums over their members,
        // which is what makes the collapsed and expanded views incapable of disagreeing (invariant
        // 3). Attributing per group and separately per atom would be two computations of the same
        // thing, and the moment they diverge the diagram lies in one of its two states.
        const size_t num_atoms = md_array_size(nto.atom_group_idx);
        if (num_atoms == 0) return false;

        int* ao_to_atom_clamped = md_temp_alloc_array(temp, int, num_ao);
        for (size_t i = 0; i < num_ao; ++i) {
            const int atom_idx = ao_to_atom[i];
            ao_to_atom_clamped[i] = (atom_idx >= 0 && (size_t)atom_idx < num_atoms) ? atom_idx : 0;
        }

        double* q_hole = md_temp_alloc_zero_array(temp, double, num_nto * num_atoms);
        double* q_part = md_temp_alloc_zero_array(temp, double, num_nto * num_atoms);

        // Mulliken populations of a single orbital are not guaranteed non-negative - that is a
        // known artifact of the partitioning, not a bug here. A negative share of a flow means
        // nothing on a diagram, so clamp and renormalize, and accept that the renormalization is
        // an approximation rather than pretending the raw numbers were clean.
        double* lambda_keep = md_temp_alloc_zero_array(temp, double, num_nto);
        double lambda_sum = 0.0;
        size_t num_keep = 0;

        for (size_t k = 0; k < num_nto; ++k) {
            if (lambdas[k] < (double)flow.lambda_cutoff) continue;

            attribute_orbital_density(q_hole + num_keep * num_atoms, ao_to_atom_clamped, C_hole + k * num_ao, S, num_ao);
            attribute_orbital_density(q_part + num_keep * num_atoms, ao_to_atom_clamped, C_part + k * num_ao, S, num_ao);

            double hsum = 0.0, psum = 0.0;
            for (size_t a = 0; a < num_atoms; ++a) {
                double* h = q_hole + num_keep * num_atoms + a;
                double* p = q_part + num_keep * num_atoms + a;
                *h = MAX(0.0, *h);
                *p = MAX(0.0, *p);
                hsum += *h;
                psum += *p;
            }
            if (hsum <= 0.0 || psum <= 0.0) continue;
            for (size_t a = 0; a < num_atoms; ++a) {
                q_hole[num_keep * num_atoms + a] /= hsum;
                q_part[num_keep * num_atoms + a] /= psum;
            }

            lambda_keep[num_keep] = lambdas[k];
            lambda_sum += lambdas[k];
            snprintf(flow.nto_label[num_keep], sizeof(flow.nto_label[0]), "NTO %zu", k + 1);
            flow_nto_character(flow.nto_character[num_keep], sizeof(flow.nto_character[0]),
                               C_hole + k * num_ao, C_part + k * num_ao, 0.35f);
            num_keep++;
        }

        if (num_keep == 0 || lambda_sum <= 0.0) return false;

        // Lambdas are truncated at MD_VLX_NTO_MAX_LAMBDAS and again by the cutoff, so they do not
        // sum to one on their own. Renormalizing over what is kept is what makes the diagram add
        // up; it also means the percentages are shares of the RETAINED transition, which the
        // window says out loud when anything was dropped.
        for (size_t k = 0; k < num_keep; ++k) {
            lambda_keep[k] /= lambda_sum;
        }

        // Per-atom column weights, and the group totals they add up to.
        double* atom_hole = md_temp_alloc_zero_array(temp, double, num_atoms);
        double* atom_part = md_temp_alloc_zero_array(temp, double, num_atoms);
        for (size_t a = 0; a < num_atoms; ++a) {
            for (size_t k = 0; k < num_keep; ++k) {
                atom_hole[a] += lambda_keep[k] * q_hole[k * num_atoms + a];
                atom_part[a] += lambda_keep[k] * q_part[k * num_atoms + a];
            }
        }

        uint32_t* atom_src = md_temp_alloc_array(temp, uint32_t, num_atoms);
        uint32_t* atom_dst = md_temp_alloc_array(temp, uint32_t, num_atoms);
        uint32_t* mid_node = md_temp_alloc_array(temp, uint32_t, num_keep);
        for (size_t a = 0; a < num_atoms; ++a) {
            atom_src[a] = MD_FLOW_INVALID_INDEX;
            atom_dst[a] = MD_FLOW_INVALID_INDEX;
        }

        // md_flow requires a parent to be an EARLIER node in the same column, so each group node
        // goes in before its atoms - which also means the group's weight has to be known first.
        // Two passes per column, and no post-hoc fixups.
        const int columns[2] = { 0, 2 };
        const double* column_weight[2] = { atom_hole, atom_part };
        uint32_t* column_atom_node[2] = { atom_src, atom_dst };

        for (int ci = 0; ci < 2; ++ci) {
            const uint32_t column = (uint32_t)columns[ci];
            const double* weight = column_weight[ci];
            uint32_t* atom_node = column_atom_node[ci];

            // Atoms are the fundamental level; a group is an OPTIONAL layer above them. Group 0 is
            // the "not in any group" bucket every atom starts in, so its atoms are roots rather
            // than children of a node standing for a group the user never made.
            for (size_t a = 0; a < num_atoms; ++a) {
                if (nto.atom_group_idx[a] != 0 || weight[a] <= 0.0) continue;
                md_flow_node_t node = {
                    .column = column, .parent = MD_FLOW_INVALID_INDEX, .level = FLOW_LEVEL_ATOM,
                    .weight = (float)weight[a],
                    .color  = flow_atom_color((uint32_t)a, column),
                    .label  = str_from_cstr(flow_atom_label((uint32_t)a)),
                    .key    = flow_node_key(column, FLOW_LEVEL_ATOM, (uint32_t)a),
                };
                atom_node[a] = md_flow_graph_add_node(&flow.graph, &node);
            }

            // md_flow requires a parent to be an EARLIER node in the same column, so each group
            // node goes in before its atoms - which also means the group's weight has to be known
            // first. Two passes per group, and no post-hoc fixups.
            for (size_t g = 1; g < num_groups; ++g) {
                double group_weight = 0.0;
                for (size_t a = 0; a < num_atoms; ++a) {
                    if (nto.atom_group_idx[a] == (uint32_t)g) group_weight += weight[a];
                }
                if (group_weight <= 0.0) continue;

                md_flow_node_t group_node = {
                    .column = column, .parent = MD_FLOW_INVALID_INDEX, .level = FLOW_LEVEL_GROUP,
                    .weight = (float)group_weight, .color = nto.group.color[g],
                    .label = str_from_cstr(nto.group.label[g]),
                    .key = flow_node_key(column, FLOW_LEVEL_GROUP, (uint32_t)g),
                };
                const uint32_t parent = md_flow_graph_add_node(&flow.graph, &group_node);
                if (parent == MD_FLOW_INVALID_INDEX) continue;

                for (size_t a = 0; a < num_atoms; ++a) {
                    if (nto.atom_group_idx[a] != (uint32_t)g || weight[a] <= 0.0) continue;

                    // An atom inside a group wears the group's colour: the expanded view has to
                    // read as a closer look at the collapsed one, not as a different quantity.
                    md_flow_node_t node = {
                        .column = column, .parent = parent, .level = FLOW_LEVEL_ATOM,
                        .weight = (float)weight[a], .color = nto.group.color[g],
                        .label = str_from_cstr(flow_atom_label((uint32_t)a)),
                        .key = flow_node_key(column, FLOW_LEVEL_ATOM, (uint32_t)a),
                    };
                    atom_node[a] = md_flow_graph_add_node(&flow.graph, &node);
                }
            }
        }

        for (size_t k = 0; k < num_keep; ++k) {
            // Middle-column nodes take their colour from the colormap rather than from a group:
            // an NTO belongs to the transition, not to either side of it.
            const ImVec4 c = ImPlot::GetColormapColor((int)k, ImPlotColormap_Deep);
            md_flow_node_t node = {
                .column = 1, .parent = MD_FLOW_INVALID_INDEX, .level = FLOW_LEVEL_GROUP,
                .weight = (float)lambda_keep[k], .color = vec_cast(c),
                .label = str_from_cstr(flow.nto_label[k]),
                .sublabel = str_from_cstr(flow.nto_character[k]),
                .key = flow_node_key(1, FLOW_LEVEL_GROUP, (uint32_t)k),
            };
            mid_node[k] = md_flow_graph_add_node(&flow.graph, &node);
        }

        // Links live between LEAVES - atoms here. Everything coarser is a sum over these, which is
        // why collapsing a group can never disagree with expanding it.
        for (size_t k = 0; k < num_keep; ++k) {
            for (size_t a = 0; a < num_atoms; ++a) {
                if (atom_src[a] != MD_FLOW_INVALID_INDEX) {
                    const float w = (float)(lambda_keep[k] * q_hole[k * num_atoms + a]);
                    if (w > 0.0f) md_flow_graph_add_link(&flow.graph, atom_src[a], mid_node[k], w);
                }
                if (atom_dst[a] != MD_FLOW_INVALID_INDEX) {
                    const float w = (float)(lambda_keep[k] * q_part[k * num_atoms + a]);
                    if (w > 0.0f) md_flow_graph_add_link(&flow.graph, mid_node[k], atom_dst[a], w);
                }
            }
        }

        return true;
    }

    // "HOMO-2", "LUMO+3" - how a chemist names an orbital, which is the entire point of showing
    // this column. Falls back to a bare index when the HOMO is unknown.
    void flow_mo_label(char* buf, size_t cap, size_t mo_idx, size_t homo_idx) {
        if (mo_idx <= homo_idx) {
            const size_t below = homo_idx - mo_idx;
            if (below == 0) snprintf(buf, cap, "HOMO");
            else            snprintf(buf, cap, "HOMO-%zu", below);
        } else {
            const size_t above = mo_idx - homo_idx - 1;
            if (above == 0) snprintf(buf, cap, "LUMO");
            else            snprintf(buf, cap, "LUMO+%zu", above);
        }
    }

    // The MO pair an NTO is mostly made of, as a label like "HOMO-1 -> LUMO".
    //
    // This is what the sketch's "pi* / sigma* / Rydberg" is reaching for, but those are symmetry
    // assignments and we have no classifier - inventing one would put a confident chemical claim on
    // screen that nothing computed. The dominant MO pair is the honest version: it is measured, it
    // is what people write in papers, and it says something the NTO index alone never can.
    //
    // Left blank unless one MO leads on each side by a clear margin. A 40/35 split has no dominant
    // character and saying otherwise would be worse than saying nothing.
    void flow_nto_character(char* buf, size_t cap, const double* C_hole_k, const double* C_part_k, float min_share) {
        buf[0] = '\0';

        const size_t num_ao = md_vlx_scf_number_of_atomic_orbitals(vlx);
        const size_t num_mo = md_vlx_scf_number_of_molecular_orbitals(vlx);
        const double* S     = md_vlx_scf_overlap_matrix_data(vlx);
        const size_t homo   = md_vlx_scf_homo_idx(vlx, MD_VLX_SPIN_ALPHA);
        if (!S || num_ao == 0 || num_mo == 0 || homo + 1 >= num_mo) return;

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };
        double* sc = md_temp_alloc_array(temp, double, num_ao);

        size_t best[2] = { num_mo, num_mo };
        double share[2] = { 0.0, 0.0 };

        for (int side = 0; side < 2; ++side) {
            const double* C_nto = (side == 0) ? C_hole_k : C_part_k;
            const size_t beg = (side == 0) ? 0    : homo + 1;
            const size_t end = (side == 0) ? homo + 1 : num_mo;

            for (size_t mu = 0; mu < num_ao; ++mu) {
                double sum = 0.0;
                for (size_t nu = 0; nu < num_ao; ++nu) sum += S[mu * num_ao + nu] * C_nto[nu];
                sc[mu] = sum;
            }

            double norm = 0.0, top = 0.0;
            size_t top_mo = num_mo;
            for (size_t mo = beg; mo < end; ++mo) {
                const double* C_mo = md_vlx_scf_mo_coefficients(vlx, mo, MD_VLX_SPIN_ALPHA);
                if (!C_mo) continue;
                double ov = 0.0;
                for (size_t mu = 0; mu < num_ao; ++mu) ov += C_mo[mu] * sc[mu];
                const double w = ov * ov;
                norm += w;
                if (w > top) { top = w; top_mo = mo; }
            }
            if (norm <= 1.0e-9 || top_mo >= num_mo) return;
            best[side] = top_mo;
            share[side] = top / norm;
        }

        if (share[0] < (double)min_share || share[1] < (double)min_share) return;

        char lo[24], hi[24];
        flow_mo_label(lo, sizeof(lo), best[0], homo);
        flow_mo_label(hi, sizeof(hi), best[1], homo);
        snprintf(buf, cap, "%s " "\xe2\x86\x92" " %s", lo, hi);
    }

    // occupied MO -> NTO pair -> virtual MO.
    //
    // NTOs ARE a rotation of the MOs: hole NTO k = sum_i U_ik phi_i over occupied MOs, particle
    // NTO k = sum_a V_ak phi_a over virtual ones. U and V are orthogonal, so |U_ik|^2 is an exact
    // normalized share and every band here is exact - no Mulliken partitioning anywhere.
    //
    // This is deliberately a route WITHOUT atoms rather than an extra column inside the atom
    // route. Chaining atom <- MO <- NTO would keep only the diagonal of
    //     q_k(A) = sum_ij U_ik U_jk P_ij(A)
    // and drop the MO-MO interference. Measured on test_data/vlx/amide.h5, that moves 33.6% of the
    // population on the dominant transition: a diagram that adds up perfectly and is wrong by a
    // third about which atoms donated. See docs/transition_flow_design.md.
    bool build_flow_graph_orbitals() {
        md_flow_graph_init(&flow.graph, 3, arena);

        if (rsp.selected < 0) return false;

        const size_t num_ao  = md_vlx_scf_number_of_atomic_orbitals(vlx);
        const size_t num_mo  = md_vlx_scf_number_of_molecular_orbitals(vlx);
        const double* S      = md_vlx_scf_overlap_matrix_data(vlx);
        const size_t homo    = md_vlx_scf_homo_idx(vlx, MD_VLX_SPIN_ALPHA);

        if (num_ao == 0 || num_mo == 0 || !S) return false;
        const size_t nocc = homo + 1;
        if (nocc == 0 || nocc >= num_mo) return false;

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        const size_t max_lambdas = MD_VLX_NTO_MAX_LAMBDAS;
        double* lambdas = md_temp_alloc_zero_array(temp, double, max_lambdas);
        double* C_hole  = md_temp_alloc_zero_array(temp, double, max_lambdas * num_ao);
        double* C_part  = md_temp_alloc_zero_array(temp, double, max_lambdas * num_ao);

        const size_t num_hole = md_vlx_rsp_nto_coefficients_extract(C_hole, lambdas, vlx, (size_t)rsp.selected, MD_VLX_NTO_HOLE, max_lambdas);
        const size_t num_part = md_vlx_rsp_nto_coefficients_extract(C_part, nullptr, vlx, (size_t)rsp.selected, MD_VLX_NTO_PARTICLE, max_lambdas);
        const size_t num_nto  = MIN(num_hole, num_part);
        if (num_nto == 0) return false;

        // Projection onto the MO basis: U_ik = <phi_i | NTO_k> = C_i^T S C_k. S*C_k is formed once
        // per NTO and then dotted with each MO, which keeps this at O(k * nao^2) rather than
        // O(k * nmo * nao^2).
        double* sc = md_temp_alloc_array(temp, double, num_ao);
        double* w_occ = md_temp_alloc_zero_array(temp, double, num_mo);
        double* w_vir = md_temp_alloc_zero_array(temp, double, num_mo);
        double* proj_occ = md_temp_alloc_zero_array(temp, double, num_nto * num_mo);
        double* proj_vir = md_temp_alloc_zero_array(temp, double, num_nto * num_mo);

        double* lambda_keep = md_temp_alloc_zero_array(temp, double, num_nto);
        double lambda_sum = 0.0;
        size_t num_keep = 0;

        for (size_t k = 0; k < num_nto; ++k) {
            if (lambdas[k] < (double)flow.lambda_cutoff) continue;

            for (int side = 0; side < 2; ++side) {
                const double* C_nto = (side == 0 ? C_hole : C_part) + k * num_ao;
                double* proj = (side == 0 ? proj_occ : proj_vir) + num_keep * num_mo;
                const size_t mo_beg = (side == 0) ? 0    : nocc;
                const size_t mo_end = (side == 0) ? nocc : num_mo;

                for (size_t mu = 0; mu < num_ao; ++mu) {
                    double sum = 0.0;
                    for (size_t nu = 0; nu < num_ao; ++nu) sum += S[mu * num_ao + nu] * C_nto[nu];
                    sc[mu] = sum;
                }

                double norm = 0.0;
                for (size_t mo = mo_beg; mo < mo_end; ++mo) {
                    const double* C_mo = md_vlx_scf_mo_coefficients(vlx, mo, MD_VLX_SPIN_ALPHA);
                    if (!C_mo) continue;
                    double overlap = 0.0;
                    for (size_t mu = 0; mu < num_ao; ++mu) overlap += C_mo[mu] * sc[mu];
                    proj[mo] = overlap * overlap;   // |U_ik|^2
                    norm += proj[mo];
                }
                // Should already be 1 by orthogonality; renormalizing absorbs the numerical drift
                // rather than letting it leak into the conservation check.
                if (norm > 1.0e-9) {
                    for (size_t mo = mo_beg; mo < mo_end; ++mo) proj[mo] /= norm;
                } else {
                    continue;
                }
            }

            snprintf(flow.nto_label[num_keep], sizeof(flow.nto_label[0]), "NTO %zu", k + 1);
            flow_nto_character(flow.nto_character[num_keep], sizeof(flow.nto_character[0]),
                               C_hole + k * num_ao, C_part + k * num_ao, 0.35f);

            lambda_keep[num_keep] = lambdas[k];
            lambda_sum += lambdas[k];
            num_keep++;
        }

        if (num_keep == 0 || lambda_sum <= 0.0) return false;
        for (size_t k = 0; k < num_keep; ++k) lambda_keep[k] /= lambda_sum;

        for (size_t k = 0; k < num_keep; ++k) {
            for (size_t mo = 0; mo < num_mo; ++mo) {
                w_occ[mo] += lambda_keep[k] * proj_occ[k * num_mo + mo];
                w_vir[mo] += lambda_keep[k] * proj_vir[k * num_mo + mo];
            }
        }

        uint32_t* occ_node = md_temp_alloc_array(temp, uint32_t, num_mo);
        uint32_t* vir_node = md_temp_alloc_array(temp, uint32_t, num_mo);
        uint32_t* mid_node = md_temp_alloc_array(temp, uint32_t, num_keep);
        for (size_t mo = 0; mo < num_mo; ++mo) { occ_node[mo] = vir_node[mo] = MD_FLOW_INVALID_INDEX; }

        // Every orbital with any weight becomes a node. There is no truncation here on purpose:
        // dropping the tail would break conservation, and the cut's threshold already folds the
        // small ones into "Others" at view time, where it can be undone.
        char label[32];
        for (size_t mo = 0; mo < num_mo; ++mo) {
            const bool occupied = mo <= homo;
            const double w = occupied ? w_occ[mo] : w_vir[mo];
            if (w <= 0.0) continue;

            flow_mo_label(label, sizeof(label), mo, homo);
            char* dst = flow_mo_label_storage((uint32_t)mo);
            if (!dst) continue;
            snprintf(dst, FLOW_ATOM_LABEL_STRIDE, "%s", label);

            md_flow_node_t node = {
                .column = occupied ? 0u : 2u,
                .parent = MD_FLOW_INVALID_INDEX,
                .level  = FLOW_LEVEL_MO,
                .weight = (float)w,
                .color  = occupied ? vec4_t{0.78f, 0.24f, 0.22f, 1.0f} : vec4_t{0.13f, 0.56f, 0.35f, 1.0f},
                .label  = str_from_cstr(dst),
                .key    = flow_node_key(occupied ? 0u : 2u, FLOW_LEVEL_MO, (uint32_t)mo),
            };
            (occupied ? occ_node : vir_node)[mo] = md_flow_graph_add_node(&flow.graph, &node);
        }

        for (size_t k = 0; k < num_keep; ++k) {
            const ImVec4 c = ImPlot::GetColormapColor((int)k, ImPlotColormap_Deep);
            md_flow_node_t node = {
                .column = 1, .parent = MD_FLOW_INVALID_INDEX, .level = FLOW_LEVEL_GROUP,
                .weight = (float)lambda_keep[k], .color = vec_cast(c),
                .label = str_from_cstr(flow.nto_label[k]),
                .sublabel = str_from_cstr(flow.nto_character[k]),
                .key = flow_node_key(1, FLOW_LEVEL_GROUP, (uint32_t)k),
            };
            mid_node[k] = md_flow_graph_add_node(&flow.graph, &node);
        }

        for (size_t k = 0; k < num_keep; ++k) {
            for (size_t mo = 0; mo < num_mo; ++mo) {
                if (occ_node[mo] != MD_FLOW_INVALID_INDEX) {
                    const float w = (float)(lambda_keep[k] * proj_occ[k * num_mo + mo]);
                    if (w > 0.0f) md_flow_graph_add_link(&flow.graph, occ_node[mo], mid_node[k], w);
                }
                if (vir_node[mo] != MD_FLOW_INVALID_INDEX) {
                    const float w = (float)(lambda_keep[k] * proj_vir[k * num_mo + mo]);
                    if (w > 0.0f) md_flow_graph_add_link(&flow.graph, mid_node[k], vir_node[mo], w);
                }
            }
        }

        return true;
    }

    void rebuild_flow_graph() {
        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        // The cut stores expansion by node INDEX, and a rebuild renumbers every node - a new
        // excited state, a regroup or a route change would otherwise leave the user's open groups
        // pointing at whatever now happens to sit at those indices. Carry the state across by KEY,
        // which is what the key field exists for.
        md_array(uint64_t) expanded_keys = 0;
        md_array(uint64_t) known_keys = 0;
        for (size_t i = 0; i < md_array_size(flow.graph.nodes); ++i) {
            md_array_push(known_keys, flow.graph.nodes[i].key, md_temp_allocator(temp));
            if (md_flow_cut_is_expanded(&flow.cut, (uint32_t)i)) {
                md_array_push(expanded_keys, flow.graph.nodes[i].key, md_temp_allocator(temp));
            }
        }

        md_flow_graph_free(&flow.graph);
        flow.built = false;
        flow.hover = FlowHit{};     // a cut index, and stale now (selection is by key, so it lives)

        const bool ok = (flow.mode == FLOW_MODE_ORBITALS) ? build_flow_graph_orbitals()
                      : (flow.mode == FLOW_MODE_NTO)      ? build_flow_graph_nto()
                                                          : build_flow_graph_groups();
        if (!ok) {
            md_flow_graph_free(&flow.graph);
            return;
        }

        // Validating every build is cheap and it is the only thing standing between a arithmetic
        // slip upstream and a diagram that looks plausible while being wrong.
        const md_flow_result_t result = md_flow_graph_validate(&flow.graph, 1.0e-3f);
        if (result != MD_FLOW_OK) {
            MD_LOG_ERROR("Charge transfer graph is inconsistent: %s", md_flow_result_str(result));
            md_flow_graph_clear(&flow.graph);
            return;
        }

        // Atoms are the level the diagram is about, so a group opens by default and its band
        // brackets its members. A group the user has deliberately closed stays closed - which is
        // why "was this key expanded" and "did this key exist at all" are different questions.
        md_flow_cut_collapse_all(&flow.cut);
        for (size_t i = 0; i < md_array_size(flow.graph.nodes); ++i) {
            const uint64_t key = flow.graph.nodes[i].key;

            bool was_expanded = false;
            for (size_t k = 0; k < md_array_size(expanded_keys); ++k) {
                if (key == expanded_keys[k]) { was_expanded = true; break; }
            }
            bool was_known = false;
            for (size_t k = 0; k < md_array_size(known_keys); ++k) {
                if (key == known_keys[k]) { was_known = true; break; }
            }

            if (was_expanded || !was_known) {
                md_flow_cut_set_expanded(&flow.cut, (uint32_t)i, true);
            }
        }

        flow.built = true;
        flow.fit_requested = true;
        md_flow_cut_resolve(&flow.cut, &flow.graph);
    }

    // Frame the diagram with room to breathe.
    //
    // Resetting to pan=0, zoom=1 makes the unit square exactly fill the canvas, which sounds right
    // and looks wrong: the layout's unit square bounds the NODES, while bands hang off the sides
    // and captions sit above the top row. Those are sized in pixels by the drawing backend, so the
    // fit has to reserve pixels for them and only then convert to a zoom.
    void flow_fit_view(ImVec2 canvas_size) {
        if (!flow.built || md_array_size(flow.layout.nodes) == 0) {
            flow.pan = {0, 0};
            flow.zoom = 1.0f;
            return;
        }

        vec2_t bb_min = { 1.0e30f,  1.0e30f };
        vec2_t bb_max = {-1.0e30f, -1.0e30f };
        for (size_t i = 0; i < md_array_size(flow.layout.nodes); ++i) {
            const md_flow_layout_node_t* n = flow.layout.nodes + i;
            bb_min.x = MIN(bb_min.x, n->min.x);  bb_min.y = MIN(bb_min.y, n->min.y);
            bb_max.x = MAX(bb_max.x, n->max.x);  bb_max.y = MAX(bb_max.y, n->max.y);
        }

        const float bb_w = MAX(bb_max.x - bb_min.x, 1.0e-4f);
        const float bb_h = MAX(bb_max.y - bb_min.y, 1.0e-4f);

        // Room for the decorations the layout knows nothing about, plus a plain margin so the
        // diagram is not welded to the window edge.
        const float margin   = 12.0f;
        const float caption  = ImGui::GetTextLineHeight() + 8.0f;
        const float band_w   = flow.style.band_offset_px + flow.style.band_width_px;
        const float band_lbl = md_array_size(flow.layout.bands) > 0 ? 64.0f : 0.0f;

        const float pad_x = margin + band_w + band_lbl;
        const float pad_y_top = margin + caption;
        const float pad_y_bot = margin;

        const float avail_w = MAX(canvas_size.x - pad_x * 2.0f, 16.0f);
        const float avail_h = MAX(canvas_size.y - pad_y_top - pad_y_bot, 16.0f);

        // Graph units are fractions of the canvas, so a graph-space extent of bb_w occupies
        // bb_w * zoom * canvas_size.x pixels.
        const float zoom_x = avail_w / (bb_w * canvas_size.x);
        const float zoom_y = avail_h / (bb_h * canvas_size.y);
        flow.zoom = CLAMP(MIN(zoom_x, zoom_y), 0.25f, 12.0f);

        // Centre the content in what is left after the padding. The vertical target is offset by
        // half the difference between the top and bottom pads, so the caption band does not push
        // the diagram off-centre.
        const vec2_t bb_center = { (bb_min.x + bb_max.x) * 0.5f, (bb_min.y + bb_max.y) * 0.5f };
        const float target_x = 0.5f;
        const float target_y = 0.5f + (pad_y_top - pad_y_bot) * 0.5f / canvas_size.y;

        flow.pan.x = target_x - bb_center.x * flow.zoom;
        flow.pan.y = target_y - bb_center.y * flow.zoom;
    }

    // Visible nodes whose on-screen rect overlaps a region, in cut index space.
    //
    // Orbital nodes are skipped: they own no atoms, so sweeping over one and calling that
    // "selecting atoms" would either do nothing or, worse, quietly select the whole molecule.
    void flow_nodes_in_region(md_array(uint32_t)* out_nodes, const FlowView& view, ImVec2 rect_min, ImVec2 rect_max,
                              md_allocator_i* alloc)
    {
        for (size_t i = 0; i < md_array_size(flow.layout.nodes); ++i) {
            const md_flow_layout_node_t* ln = flow.layout.nodes + i;
            const md_flow_node_t* node = md_flow_cut_node(&flow.cut, &flow.graph, ln->cut_idx);
            if (!node || !flow_level_has_atoms(node->level)) continue;

            const ImVec2 min = flow_view_to_screen(view, ln->min);
            const ImVec2 max = flow_view_to_screen(view, ln->max);

            // Overlap, not containment: a sweep across a column should catch the rows it crosses,
            // not only the ones it happens to swallow whole.
            if (max.x < rect_min.x || min.x > rect_max.x) continue;
            if (max.y < rect_min.y || min.y > rect_max.y) continue;

            md_array_push(*out_nodes, ln->cut_idx, alloc);
        }
    }

    // Group membership from whatever is currently selected - wherever that selection was made.
    void flow_assign_selection_to_group(const ApplicationState& state, uint32_t group_idx) {
        if (group_idx >= nto.group.count) return;
        for (size_t i = 0; i < md_array_size(nto.atom_group_idx); ++i) {
            if (md_bitfield_test_bit(&state.selection.selection_mask, i)) {
                nto.atom_group_idx[i] = group_idx;
            }
        }
    }

    // Returns the new group index, or 0 (the ungrouped bucket) when there is no room.
    uint32_t flow_create_group_from_selection(const ApplicationState& state) {
        if (nto.group.count >= MAX_NTO_GROUPS) {
            MD_LOG_ERROR("Cannot create another group: the limit of %d is reached", (int)MAX_NTO_GROUPS);
            return 0;
        }
        const uint32_t group_idx = (uint32_t)nto.group.count++;
        flow_assign_selection_to_group(state, group_idx);
        return group_idx;
    }

    bool flow_orbital_is_selected(uint64_t key) const {
        for (size_t i = 0; i < md_array_size(flow.selected_orbital_keys); ++i) {
            if (flow.selected_orbital_keys[i] == key) return true;
        }
        return false;
    }

    void flow_orbital_set_selected(uint64_t key, bool selected) {
        for (size_t i = 0; i < md_array_size(flow.selected_orbital_keys); ++i) {
            if (flow.selected_orbital_keys[i] != key) continue;
            if (!selected) {
                flow.selected_orbital_keys[i] = flow.selected_orbital_keys[md_array_size(flow.selected_orbital_keys) - 1];
                md_array_shrink(flow.selected_orbital_keys, md_array_size(flow.selected_orbital_keys) - 1);
            }
            return;
        }
        if (selected) md_array_push(flow.selected_orbital_keys, key, arena);
    }

    // One node's contribution to a selection change. Atom-bearing nodes edit the application mask
    // directly - there is no diagram-local copy to keep in step - while NTO nodes are a separate
    // axis and never touch it.
    void flow_set_selected(ApplicationState& state, uint32_t cut_idx, bool selected) {
        const md_flow_node_t* node = md_flow_cut_node(&flow.cut, &flow.graph, cut_idx);
        if (!node) return;

        if (!flow_level_has_atoms(node->level)) {
            flow_orbital_set_selected(node->key, selected);
            return;
        }

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        md_bitfield_t mask = md_bitfield_create(md_temp_allocator(temp));
        flow_node_atom_mask(&mask, cut_idx);
        if (selected) {
            md_bitfield_or_inplace(&state.selection.selection_mask, &mask);
        } else {
            md_bitfield_andnot_inplace(&state.selection.selection_mask, &mask);
        }
    }

    void flow_clear_selection(ApplicationState& state) {
        md_bitfield_clear(&state.selection.selection_mask);
        md_array_shrink(flow.selected_orbital_keys, 0);
    }

    // Atoms belonging to a visible node, for the two-way link with the viewport.
    void flow_node_atom_mask(md_bitfield_t* out_mask, uint32_t cut_idx) {
        md_bitfield_clear(out_mask);
        const md_flow_node_t* node = md_flow_cut_node(&flow.cut, &flow.graph, cut_idx);
        if (!node || (node->flags & MD_FLOW_NODE_FLAG_OTHER)) return;

        // Orbitals own no atoms; selecting "all the atoms an NTO touches" would be the whole
        // molecule and mean nothing. Gated on the LEVEL rather than the column, because in the
        // Orbitals route the outer columns are orbitals too.
        const md_flow_node_t* n = node;
        if (!flow_level_has_atoms(n->level)) return;

        const uint32_t index = flow_key_index(n->key);
        if (flow_key_level(n->key) == FLOW_LEVEL_ATOM) {
            md_bitfield_set_bit(out_mask, index);
            return;
        }
        for (size_t i = 0; i < md_array_size(nto.atom_group_idx); ++i) {
            if (nto.atom_group_idx[i] == index) {
                md_bitfield_set_bit(out_mask, i);
            }
        }
    }

    void draw_flow_window(ApplicationState& state) {
        if (!flow.show_window) return;

        const size_t num_excited_states = md_vlx_rsp_number_of_excited_states(vlx);
        if (num_excited_states == 0) return;

        if (!flow.cut.alloc) {
            md_flow_cut_init(&flow.cut, arena);
            md_flow_layout_init(&flow.layout, arena);
            flow.params = md_flow_layout_params_default();
            flow.cut.threshold = flow.initial_threshold;
            // Uncapped on purpose. A cap makes "Show above 1%" stop being true - folding halts
            // early and the LARGEST sub-threshold rows stay on screen, with nothing to tell the
            // reader why. A big "Others" is fine now that it opens.
            flow.cut.other_max = 1.0f;
            flow_build_label_storage();
        }

        // Rebuild only when the numbers behind the picture actually changed. Hashing the matrix
        // itself rather than a bookkeeping counter keeps this window independent of whatever the
        // NTO window happens to be doing.
        uint64_t hash = md_hash64(nto.atom_group_idx, md_array_bytes(nto.atom_group_idx), 0);
        hash = md_hash64(nto.group.color, sizeof(nto.group.color), hash);
        hash ^= (uint64_t)nto.group.count;
        hash ^= (uint64_t)(rsp.selected + 1) << 8;
        hash ^= (uint64_t)flow.mode << 24;
        if (nto.transition_matrix) {
            const size_t dim = nto.transition_matrix_dim;
            hash = md_hash64(nto.transition_matrix, sizeof(float) * dim * (dim + 2), hash);
        }
        if (hash != flow.data_hash) {
            flow.data_hash = hash;
            rebuild_flow_graph();
        }

        bool export_requested = false;

        ImGui::SetNextWindowSize(ImVec2(720, 420), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Charge Transfer Analysis", &flow.show_window, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("File")) {
                    export_requested = ImGui::MenuItem("Export Diagram (SVG)", nullptr, false, flow.built);
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Settings")) {
                    ImGui::SliderFloat("Box tint", &flow.style.node_fill_tint, 0.0f, 0.5f);
                    ImGui::SetItemTooltip("How much of a node's colour survives in its interior.\n"
                                          "0 is white; the colour itself lives on the border.");
                    ImGui::SliderFloat("Border width", &flow.style.border_width, 0.5f, 4.0f, "%.1f px");
                    ImGui::Checkbox("Colour the labels", &flow.style.color_labels);
                    ImGui::Separator();
                    ImGui::SliderFloat("Ribbon opacity", &flow.style.ribbon_alpha, 0.05f, 1.0f);
                    ImGui::SliderFloat("Ribbon curvature", &flow.style.ribbon_curvature, 0.0f, 1.0f);
                    ImGui::SliderInt("Ribbon segments", &flow.style.ribbon_segments, 4, 64);
                    ImGui::Separator();
                    ImGui::SliderFloat("Node thickness", &flow.params.node_thickness, 0.01f, 0.25f);
                    ImGui::SliderFloat("Node gap", &flow.params.node_gap, 0.0f, 0.1f);
                    ImGui::BeginDisabled(flow.params.order != MD_FLOW_ORDER_BARYCENTRE);
                    int sweeps = (int)flow.params.barycentre_sweeps;
                    if (ImGui::SliderInt("Crossing sweeps", &sweeps, 0, 4)) {
                        flow.params.barycentre_sweeps = (uint32_t)MAX(0, sweeps);
                    }
                    ImGui::EndDisabled();
                    ImGui::SetItemTooltip("Only used by the 'Fewest crossings' row order.");
                    ImGui::Separator();
                    ImGui::SliderFloat("Hide ribbons under", &flow.style.ribbon_min_px, 0.0f, 6.0f, "%.1f px");
                    ImGui::SetItemTooltip("A ribbon thinner than this is not drawn - it reads as noise rather\n"
                                          "than as a quantity. The flow is still in the model and still counted.");
                    ImGui::EndMenu();
                }
                ImGui::EndMenuBar();
            }

            ImGui::PushItemWidth(160.0f);
            int state_idx = rsp.selected + 1;
            if (ImGui::SliderInt("State", &state_idx, 1, (int)num_excited_states)) {
                rsp.selected = state_idx - 1;
            }
            ImGui::PopItemWidth();

            ImGui::SameLine();
            ImGui::PushItemWidth(130.0f);
            int mode = (int)flow.mode;
            if (ImGui::Combo("Route", &mode, flow_mode_str, FLOW_MODE_COUNT)) {
                flow.mode = (flow_mode_t)mode;
            }
            ImGui::SetItemTooltip("'Groups only' reproduces the old two-column diagram, where the\n"
                                  "group-to-group flow is inferred from the two marginals.\n"
                                  "'Atoms via NTOs' measures it instead.\n"
                                  "'Orbitals' drops atoms entirely and shows which canonical MOs the\n"
                                  "NTOs are built from - every band there is exact.");
            ImGui::PopItemWidth();

            ImGui::SameLine();
            ImGui::PushItemWidth(160.0f);
            float threshold_pct = flow.cut.threshold * 100.0f;
            if (ImGui::SliderFloat("Show above", &threshold_pct, 0.0f, 25.0f, "%.1f%%")) {
                flow.cut.threshold = threshold_pct * 0.01f;
                md_flow_cut_resolve(&flow.cut, &flow.graph);
            }
            ImGui::SetItemTooltip("Share of a COLUMN, not of the whole diagram. Everything below folds\n"
                                  "into 'Others'; double click that to draw its members individually.");
            ImGui::PopItemWidth();

            ImGui::SameLine();
            ImGui::PushItemWidth(140.0f);
            float other_pct = flow.cut.other_max * 100.0f;
            if (ImGui::SliderFloat("Others cap", &other_pct, 1.0f, 100.0f, "%.0f%%")) {
                flow.cut.other_max = other_pct * 0.01f;
                md_flow_cut_resolve(&flow.cut, &flow.graph);
            }
            ImGui::SetItemTooltip("The most a single 'Others' may hold. Below 100%% folding stops early,\n"
                                  "which means sub-threshold rows STAY on screen - so 'Show above' no\n"
                                  "longer means what it says. Leave it at 100%% unless you want that.");
            ImGui::PopItemWidth();

            ImGui::SameLine();
            ImGui::PushItemWidth(140.0f);
            ImGui::SliderFloat("Even sizes", &flow.size_normalize, 0.0f, 1.0f, "%.2f");
            ImGui::SetItemTooltip("Flattens the ORBITAL columns' heights so a 1%% channel is still big\n"
                                  "enough to read and click. Atom columns stay strictly proportional -\n"
                                  "their sizes are the charge transfer itself.\n"
                                  "Ribbons keep their true widths and simply taper into the column.");
            ImGui::PopItemWidth();

            ImGui::SameLine();
            // 'Groups only' has no atom level to open. Its links come from the two marginals and
            // carry no information about which ATOM within a group donated to which - splitting
            // them proportionally would invent detail one level finer than the data supports,
            // which is the exact failure this window exists to stop doing.
            const bool expandable = (flow.mode == FLOW_MODE_NTO);   // only that route has a hierarchy
            ImGui::BeginDisabled(!expandable);
            if (ImGui::Button("Expand all")) {
                md_flow_cut_expand_all(&flow.cut, &flow.graph);
                md_flow_cut_resolve(&flow.cut, &flow.graph);
            }
            ImGui::SameLine();
            if (ImGui::Button("Collapse all")) {
                md_flow_cut_collapse_all(&flow.cut);
                md_flow_cut_resolve(&flow.cut, &flow.graph);
            }
            ImGui::EndDisabled();
            if (!expandable) {
                ImGui::SetItemTooltip("Only 'Atoms via NTOs' has a hierarchy to open.\n"
                                      "'Groups only' infers its flow from the two marginals and says\n"
                                      "nothing about which atom went where; 'Orbitals' has no atoms.");
            }

            ImGui::SameLine();
            ImGui::PushItemWidth(150.0f);
            int order = (int)flow.params.order;
            static const char* order_str[] = { "By contribution", "Fewest crossings" };
            if (ImGui::Combo("Order", &order, order_str, IM_ARRAYSIZE(order_str))) {
                flow.params.order = (md_flow_order_t)order;
            }
            ImGui::SetItemTooltip("Rows within a column. 'By contribution' puts the largest at the top;\n"
                                  "either way a group's atoms stay together inside the group's slot.");
            ImGui::PopItemWidth();

            ImGui::SameLine();
            ImGui::TextDisabled("(?)");
            ImGui::SetItemTooltip("Left drag pans, wheel zooms at the cursor.\n"
                                  "Shift+left selects, shift+right deselects - click one item,\n"
                                  "or drag a box over several. Shift+right on empty space clears.\n"
                                  "Double click opens a node; double click on empty space refits.\n"
                                  "Right click for group actions on the selection.");

            const ImVec2 canvas_p0 = ImGui::GetCursorScreenPos();
            ImVec2 canvas_sz = ImGui::GetContentRegionAvail();
            canvas_sz.x = MAX(canvas_sz.x, 64.0f);
            canvas_sz.y = MAX(canvas_sz.y, 64.0f);
            const ImVec2 canvas_p1 = canvas_p0 + canvas_sz;

            ImGui::InvisibleButton("##flow_canvas", canvas_sz,
                ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_MouseButtonMiddle);
            const bool canvas_hovered = ImGui::IsItemHovered();
            const bool canvas_active  = ImGui::IsItemActive();

            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(255, 255, 255, 255));
            draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(0, 0, 0, 60));

            if (!flow.built) {
                const char* msg = "No transition data for this state";
                draw_list->AddText(canvas_p0 + (canvas_sz - ImGui::CalcTextSize(msg)) * 0.5f, IM_COL32(0, 0, 0, 160), msg);
                ImGui::End();
                return;
            }

            ImGuiIO& io = ImGui::GetIO();

            // Zoom about the cursor. Screen = area + (p * zoom + pan) * size, so the graph point
            // under the pointer is p = (rel - pan) / zoom, and holding it fixed across a zoom
            // change means pan' = rel - p * zoom'. That is an ASSIGNMENT, not an increment - the
            // increment form leaves an extra pan term in and the view creeps toward the corner.
            if (canvas_hovered && io.MouseWheel != 0.0f) {
                const ImVec2 rel = (io.MousePos - canvas_p0) / canvas_sz;
                const vec2_t p = { (rel.x - flow.pan.x) / flow.zoom, (rel.y - flow.pan.y) / flow.zoom };
                flow.zoom = CLAMP(flow.zoom * (1.0f + io.MouseWheel * 0.12f), 0.25f, 12.0f);
                flow.pan.x = rel.x - p.x * flow.zoom;
                flow.pan.y = rel.y - p.y * flow.zoom;
            }

            // Left drag pans. Selection moved onto shift, so the plain drag is free for the thing
            // people reach for first - and a shift drag is a region, so panning has to stand down
            // for it or the diagram would slide out from under the rubber band.
            if (canvas_active && !flow.region_active && !io.KeyShift &&
                                 (ImGui::IsMouseDragging(ImGuiMouseButton_Left) ||
                                  ImGui::IsMouseDragging(ImGuiMouseButton_Middle) ||
                                  ImGui::IsMouseDragging(ImGuiMouseButton_Right))) {
                flow.pan.x += io.MouseDelta.x / canvas_sz.x;
                flow.pan.y += io.MouseDelta.y / canvas_sz.y;
            }

            // Map the slider onto a per-column exponent. Orbital columns get it; atom columns stay
            // proportional. In the Orbitals route every column is orbitals.
            {
                const float e = 1.0f - flow.size_normalize * 0.85f;   // 1.0 .. 0.15
                for (uint32_t c = 0; c < MD_FLOW_MAX_COLUMNS; ++c) flow.params.weight_exponent[c] = 1.0f;
                flow.params.weight_exponent[1] = e;
                if (flow.mode == FLOW_MODE_ORBITALS) {
                    flow.params.weight_exponent[0] = e;
                    flow.params.weight_exponent[2] = e;
                }
            }

            md_flow_layout_compute(&flow.layout, &flow.graph, &flow.cut, &flow.params);

            // The layout lives in graph space and knows nothing about pan or zoom, so fitting
            // after it is computed costs nothing and needs no second pass.
            if (flow.fit_requested) {
                flow.fit_requested = false;
                flow_fit_view(canvas_sz);
            }

            FlowView view = {};
            view.area_min = canvas_p0;
            view.area_max = canvas_p1;
            view.pan  = { flow.pan.x, flow.pan.y };
            view.zoom = flow.zoom;

            // Node state is DERIVED from the application masks every frame rather than stored.
            // That is what makes a selection made in the 3D view, or by a script, light up here
            // without any event plumbing - and it is why there is no local copy to fall out of
            // step with the rest of VIAMD.
            md_temp_scope_t draw_temp = md_temp_begin();
            defer { md_temp_end(draw_temp); };

            const size_t cut_space = md_array_size(flow.graph.nodes) + md_array_size(flow.cut.other);
            uint8_t* node_state = md_temp_alloc_zero_array(draw_temp, uint8_t, cut_space ? cut_space : 1);

            const size_t num_atoms = md_array_size(nto.atom_group_idx);

            // One pass over the atoms gives both the per-atom answer and the per-group tallies, so
            // a group node is a lookup rather than a rescan.
            uint32_t* group_total = md_temp_alloc_zero_array(draw_temp, uint32_t, MAX_NTO_GROUPS);
            uint32_t* group_sel   = md_temp_alloc_zero_array(draw_temp, uint32_t, MAX_NTO_GROUPS);
            uint32_t* group_hi    = md_temp_alloc_zero_array(draw_temp, uint32_t, MAX_NTO_GROUPS);
            bool selection_touches_atoms = false;

            for (size_t a = 0; a < num_atoms; ++a) {
                const uint32_t g = nto.atom_group_idx[a] < MAX_NTO_GROUPS ? nto.atom_group_idx[a] : 0;
                const bool sel = md_bitfield_test_bit(&state.selection.selection_mask, a);
                const bool hi  = md_bitfield_test_bit(&state.selection.highlight_mask, a);
                group_total[g] += 1;
                group_sel[g]   += sel ? 1 : 0;
                group_hi[g]    += hi  ? 1 : 0;
                selection_touches_atoms |= sel;
            }

            for (size_t i = 0; i < md_array_size(flow.graph.nodes); ++i) {
                const md_flow_node_t* n = flow.graph.nodes + i;
                uint8_t st = 0;

                if (!flow_level_has_atoms(n->level)) {
                    if (flow_orbital_is_selected(n->key)) st |= FlowNodeState_Selected;
                } else {
                    const uint32_t index = flow_key_index(n->key);
                    if (flow_key_level(n->key) == FLOW_LEVEL_ATOM) {
                        if (index < num_atoms) {
                            if (md_bitfield_test_bit(&state.selection.selection_mask, index)) st |= FlowNodeState_Selected;
                            if (md_bitfield_test_bit(&state.selection.highlight_mask, index)) st |= FlowNodeState_Highlight;
                        }
                    } else if (index < MAX_NTO_GROUPS && group_total[index] > 0) {
                        // A group is only "selected" when all of it is. Anything less is partial -
                        // otherwise one selected atom would claim the whole group, and there would
                        // be no way to see the difference.
                        if (group_sel[index] == group_total[index])  st |= FlowNodeState_Selected;
                        else if (group_sel[index] > 0)               st |= FlowNodeState_Partial;
                        if (group_hi[index] > 0)                     st |= FlowNodeState_Highlight;
                    }
                }
                node_state[i] = st;
            }

            FlowEmphasis emphasis = {};
            emphasis.node_state = node_state;
            emphasis.num_state_entries = cut_space;
            emphasis.any_selected = selection_touches_atoms || md_array_size(flow.selected_orbital_keys) > 0;
            emphasis.hover_node = flow.hover.node;
            emphasis.hover_link = flow.hover.link;

            draw_list->PushClipRect(canvas_p0, canvas_p1, true);
            const FlowHit hit = flow_draw_imgui(draw_list, view, &flow.graph, &flow.cut, &flow.layout,
                                                flow.style, emphasis);
            draw_list->PopClipRect();

            flow.hover = canvas_hovered ? hit : FlowHit{};

            // Column captions. They belong to the window rather than to the backend: the backend
            // draws a graph and has no idea what its columns mean.
            // Captions belong to the COLUMNS, so they live in graph space and pan and zoom with
            // them. Pinned to the window corners they would drift away from what they name the
            // moment anyone moved the view.
            {
                const ImU32 caption_col = IM_COL32(0, 0, 0, 200);
                const char* caption_atoms[3]    = { "DONATING (hole)", "NTO PAIRS", "ACCEPTING (particle)" };
                const char* caption_orbitals[3] = { "OCCUPIED MOs",    "NTO PAIRS", "VIRTUAL MOs" };
                const char** caption = (flow.mode == FLOW_MODE_ORBITALS) ? caption_orbitals : caption_atoms;

                for (uint32_t c = 0; c < flow.graph.num_columns; ++c) {
                    float x0 = 0.0f, x1 = 0.0f, top = 1.0e30f;
                    bool found = false;
                    for (size_t i = 0; i < md_array_size(flow.layout.nodes); ++i) {
                        const md_flow_node_t* n = md_flow_cut_node(&flow.cut, &flow.graph, flow.layout.nodes[i].cut_idx);
                        if (!n || n->column != c) continue;
                        x0 = flow.layout.nodes[i].min.x;
                        x1 = flow.layout.nodes[i].max.x;
                        top = MIN(top, flow.layout.nodes[i].min.y);
                        found = true;
                    }
                    if (!found) continue;

                    const char* lbl = (flow.graph.num_columns == 2 && c == 1) ? caption[2] : caption[MIN(c, 2u)];
                    const ImVec2 mid = flow_view_to_screen(view, vec2_t{ (x0 + x1) * 0.5f, top });
                    const ImVec2 size = ImGui::CalcTextSize(lbl);
                    draw_list->AddText({ mid.x - size.x * 0.5f, mid.y - size.y - 6.0f }, caption_col, lbl);
                }
            }

            // ---- selection: rubber band, or a click when the band never opened ----
            //
            // Press, drag and release are one gesture, so they are handled in one place rather than
            // split across the "is something under the cursor" branches - a drag that starts on a
            // node routinely ends over empty space, and the release is the half that matters.
            bool region_preview = false;
            {
                // Every button here also pans, so a "click" means pressed and released without
                // travelling. Deciding that at RELEASE, from the distance actually covered, is what
                // lets one gesture be either a click or a band without the user declaring which.
                const float slop = 4.0f;

                if (canvas_hovered && io.KeyShift && !flow.region_active &&
                    (ImGui::IsMouseClicked(ImGuiMouseButton_Left) || ImGui::IsMouseClicked(ImGuiMouseButton_Right)))
                {
                    flow.region_active = true;
                    flow.region_add    = ImGui::IsMouseClicked(ImGuiMouseButton_Left);
                    flow.region_start  = io.MousePos;
                }

                if (flow.region_active) {
                    const ImGuiMouseButton button = flow.region_add ? ImGuiMouseButton_Left : ImGuiMouseButton_Right;
                    const ImVec2 rect_min = { ImMin(flow.region_start.x, io.MousePos.x), ImMin(flow.region_start.y, io.MousePos.y) };
                    const ImVec2 rect_max = { ImMax(flow.region_start.x, io.MousePos.x), ImMax(flow.region_start.y, io.MousePos.y) };
                    const bool dragged = ImLengthSqr(io.MousePos - flow.region_start) >= slop * slop;

                    md_temp_scope_t rtemp = md_temp_begin();
                    defer { md_temp_end(rtemp); };
                    md_allocator_i* ralloc = md_temp_allocator(rtemp);

                    if (ImGui::IsMouseDown(button)) {
                        if (dragged) {
                            region_preview = true;
                            const ImU32 fill = flow.region_add ? IM_COL32(60, 130, 220, 48) : IM_COL32(210, 70, 70, 48);
                            const ImU32 line = flow.region_add ? IM_COL32(60, 130, 220, 200) : IM_COL32(210, 70, 70, 200);
                            draw_list->PushClipRect(canvas_p0, canvas_p1, true);
                            draw_list->AddRectFilled(rect_min, rect_max, fill);
                            draw_list->AddRect(rect_min, rect_max, line);
                            draw_list->PopClipRect();

                            // Preview into the highlight mask so the 3D view shows what the release
                            // is about to do, for both add and remove.
                            md_array(uint32_t) nodes = 0;
                            flow_nodes_in_region(&nodes, view, rect_min, rect_max, ralloc);

                            md_bitfield_t node_mask = md_bitfield_create(ralloc);
                            md_bitfield_clear(&state.selection.highlight_mask);
                            for (size_t i = 0; i < md_array_size(nodes); ++i) {
                                flow_node_atom_mask(&node_mask, nodes[i]);
                                md_bitfield_or_inplace(&state.selection.highlight_mask, &node_mask);
                            }
                        }
                    } else {
                        if (dragged) {
                            md_array(uint32_t) nodes = 0;
                            flow_nodes_in_region(&nodes, view, rect_min, rect_max, ralloc);
                            for (size_t i = 0; i < md_array_size(nodes); ++i) {
                                flow_set_selected(state, nodes[i], flow.region_add);
                            }
                        } else if (hit.node >= 0) {
                            // Never opened into a band: the same gesture on a single item.
                            flow_set_selected(state, (uint32_t)hit.node, flow.region_add);
                        } else if (!flow.region_add) {
                            // Shift+right on nothing clears everything - the gesture that removes
                            // one item, aimed at no item.
                            flow_clear_selection(state);
                        }
                        flow.region_active = false;
                    }
                }
            }

            // Plain right click (no shift, no drag) opens the group menu. Shift+right is the
            // deselect gesture and a right DRAG pans, so this is the one right-button meaning left
            // over - and it is the one that needs somewhere to live, since the new window had no
            // way to make a group at all.
            if (canvas_hovered && !io.KeyShift && !flow.region_active &&
                ImGui::IsMouseReleased(ImGuiMouseButton_Right) &&
                io.MouseDragMaxDistanceSqr[ImGuiMouseButton_Right] < 16.0f)
            {
                ImGui::OpenPopup("##flow_context");
            }

            if (ImGui::BeginPopup("##flow_context")) {
                const size_t num_selected = md_bitfield_popcount(&state.selection.selection_mask);
                if (num_selected == 0) {
                    ImGui::TextDisabled("No atoms selected");
                    ImGui::Separator();
                    ImGui::TextDisabled("Shift+click or shift+drag to select");
                } else {
                    ImGui::Text("%zu atom%s selected", num_selected, num_selected == 1 ? "" : "s");
                    ImGui::Separator();

                    if (ImGui::MenuItem("Create new group", nullptr, false, nto.group.count < MAX_NTO_GROUPS)) {
                        flow_create_group_from_selection(state);
                    }
                    if (ImGui::BeginMenu("Add to group", nto.group.count > 1)) {
                        // Slot 0 is the ungrouped bucket, not a destination anyone means to pick.
                        for (size_t g = 1; g < nto.group.count; ++g) {
                            if (ImGui::MenuItem(nto.group.label[g])) {
                                flow_assign_selection_to_group(state, (uint32_t)g);
                            }
                        }
                        ImGui::EndMenu();
                    }
                    if (ImGui::MenuItem("Remove from groups")) {
                        flow_assign_selection_to_group(state, 0);
                    }
                    ImGui::Separator();
                    if (ImGui::MenuItem("Clear selection")) {
                        flow_clear_selection(state);
                    }
                }
                ImGui::EndPopup();
            }

            if (canvas_hovered && hit.node >= 0) {
                const md_flow_node_t* node = md_flow_cut_node(&flow.cut, &flow.graph, (uint32_t)hit.node);
                if (node && (node->flags & MD_FLOW_NODE_FLAG_OTHER)) {
                    // An "Others" is a sum, so its tooltip has to say what it is a sum OF - an
                    // unexplained grey box holding a chunk of the transfer is worse than no box.
                    const uint32_t other_idx = (uint32_t)hit.node - (uint32_t)md_array_size(flow.graph.nodes);
                    const uint32_t count = other_idx < md_array_size(flow.cut.other_count)
                                         ? flow.cut.other_count[other_idx] : 0;
                    ImGui::SetTooltip("Others\n%.1f%% across %u entries, each below %.1f%%\n\n"
                                      "Double click to draw them individually",
                                      node->weight * 100.0f, count, flow.cut.threshold * 100.0f);
                } else if (node) {
                    const char* what = (node->column == 0) ? "detachment"
                                     : (node->column == 1) ? "transition" : "attachment";
                    const bool has_children = (uint32_t)hit.node < md_array_size(flow.graph.nodes)
                                              && !md_flow_node_is_leaf(&flow.graph, (uint32_t)hit.node);
                    if (has_children) {
                        const bool open = md_flow_cut_is_expanded(&flow.cut, (uint32_t)hit.node);
                        ImGui::SetTooltip("%.*s\n%.1f%% of the %s\n\nDouble click to %s",
                            (int)node->label.len, node->label.ptr, node->weight * 100.0f, what,
                            open ? "collapse" : "show its atoms");
                    } else {
                        ImGui::SetTooltip("%.*s\n%.1f%% of the %s",
                            (int)node->label.len, node->label.ptr, node->weight * 100.0f, what);
                    }
                }
                // Hovering previews in the viewport; committing to a selection takes shift, so a
                // plain drag can pan without the selection changing under the cursor. While a band
                // is open it owns the preview, or the node under the cursor would keep overwriting
                // what the band is about to do.
                if (!region_preview) {
                    flow_node_atom_mask(&state.selection.highlight_mask, (uint32_t)hit.node);
                }

                // Double click opens or closes a node - including an "Others", which opens by
                // suspending the threshold for its column. That is what lets the threshold hide
                // things honestly: everything below it goes, and you can still look.
                if (!io.KeyShift && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                    const uint32_t idx = (uint32_t)hit.node;
                    const md_flow_node_t* n = md_flow_cut_node(&flow.cut, &flow.graph, idx);
                    if (n && (n->flags & MD_FLOW_NODE_FLAG_OTHER)) {
                        md_flow_cut_set_other_expanded(&flow.cut, n->column, true);
                        md_flow_cut_resolve(&flow.cut, &flow.graph);
                    } else if (idx < md_array_size(flow.graph.nodes) && !md_flow_node_is_leaf(&flow.graph, idx)) {
                        md_flow_cut_set_expanded(&flow.cut, idx, !md_flow_cut_is_expanded(&flow.cut, idx));
                        md_flow_cut_resolve(&flow.cut, &flow.graph);
                    } else if (n && md_flow_cut_is_other_expanded(&flow.cut, n->column)) {
                        // Double clicking a node in an opened column closes it again - the pair of
                        // gestures is symmetric, since the "Others" itself is gone while open.
                        md_flow_cut_set_other_expanded(&flow.cut, n->column, false);
                        md_flow_cut_resolve(&flow.cut, &flow.graph);
                    }
                }
            } else if (canvas_hovered && hit.node < 0 && hit.link < 0) {
                // Cursor is in this canvas and over nothing, so whatever this window last put in
                // the highlight mask is stale. Cleared only under those two conditions: clearing it
                // unconditionally would stamp on the 3D viewport's own hover.
                if (!region_preview) {
                    md_bitfield_clear(&state.selection.highlight_mask);
                }

                // Empty space: double click reframes. A button for this would occupy toolbar room
                // permanently to undo something that only occasionally needs undoing.
                if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                    flow_fit_view(canvas_sz);
                }
            }
            if (canvas_hovered && hit.link >= 0) {
                const md_flow_link_t* link = flow.cut.links + hit.link;
                const md_flow_node_t* src = md_flow_cut_node(&flow.cut, &flow.graph, link->src);
                const md_flow_node_t* dst = md_flow_cut_node(&flow.cut, &flow.graph, link->dst);
                if (src && dst) {
                    ImGui::SetTooltip("%.*s  ->  %.*s\n%.2f%%",
                        (int)src->label.len, src->label.ptr,
                        (int)dst->label.len, dst->label.ptr,
                        link->weight * 100.0f);
                }
            }
        }
        ImGui::End();

        if (export_requested) {
            export_flow_diagram(ImVec2(1024, 640));
        }
    }

    bool export_flow_diagram(ImVec2 size) {
        if (!flow.built) {
            MD_LOG_ERROR("Charge transfer export failed: no diagram to export");
            return false;
        }

        char path_buf[1024] = {};
        if (!application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, STR_LIT("svg"))) {
            return false;
        }

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        md_vg_scene_t scene = {};
        md_vg_scene_init(&scene, { size.x, size.y }, md_temp_allocator(temp));

        // The export draws the SAME layout the screen does, at a different size. That is the
        // point of keeping layout separate from drawing: there is no second implementation to
        // drift out of agreement with the first.
        md_flow_layout_t layout = {};
        md_flow_layout_init(&layout, md_temp_allocator(temp));
        md_flow_layout_compute(&layout, &flow.graph, &flow.cut, &flow.params);

        FlowView view = {};
        view.area_min = { 0.0f, 0.0f };
        view.area_max = size;
        view.zoom = 1.0f;

        flow_draw_vg(&scene, view, &flow.graph, &flow.cut, &layout, flow.style);

        const bool success = md_vg_scene_write_svg_file(&scene, str_from_cstr(path_buf));
        if (success) {
            MD_LOG_INFO("Exported charge transfer diagram to '%s'", path_buf);
        } else {
            MD_LOG_ERROR("Failed to export charge transfer diagram to '%s'", path_buf);
        }
        return success;
    }

    void update_nto_derived_data() {
        if (!vlx || !md_vlx_rsp_has_nto(vlx) || !nto.atom_group_idx) return;

        // Create hash to check for changes to trigger update of colors in gl representation
        uint64_t atom_idx_hash = md_hash64(nto.atom_group_idx, md_array_bytes(nto.atom_group_idx), 0);
        uint64_t color_hash    = md_hash64(nto.group.color, sizeof(nto.group.color), atom_idx_hash);

        static uint64_t cur_color_hash = 0;
        if (color_hash != cur_color_hash) {
            cur_color_hash = color_hash;
            update_nto_group_colors();
        }

        // Create hash to check for changes to trigger recomputation of transition matrix
        uint64_t matrix_hash = atom_idx_hash ^ (uint64_t)rsp.selected ^ nto.group.count;
        static uint64_t cur_matrix_hash = 0;

        if (matrix_hash != cur_matrix_hash) {
            cur_matrix_hash = matrix_hash;

            // Resize transition matrix to the correct size
            if (nto.transition_matrix_dim != nto.group.count) {
                if (nto.transition_matrix) {
                    // The allocated size contains matrix N*N + 2*N for hole/part arrays
                    const size_t cur_mem = sizeof(float) * nto.transition_matrix_dim * (nto.transition_matrix_dim + 2);
                    md_free(arena, nto.transition_matrix, cur_mem);
                }

                nto.transition_matrix_dim = nto.group.count;
                const size_t new_mem = sizeof(float) * nto.transition_matrix_dim * (nto.transition_matrix_dim + 2);
                nto.transition_matrix = (float*)md_alloc(arena, new_mem);
                nto.transition_density_hole = nto.transition_matrix + (nto.transition_matrix_dim * nto.transition_matrix_dim);
                nto.transition_density_part = nto.transition_density_hole + nto.transition_matrix_dim;
            }

            // Clear all of the values to zero (matrix, group values)
            MEMSET(nto.transition_matrix, 0, sizeof(float) * nto.transition_matrix_dim * (nto.transition_matrix_dim + 2));

			if (rsp.selected != -1) {
				const size_t nto_idx = (size_t)rsp.selected;
#if 1
				{
					md_temp_scope_t temp = md_temp_begin();
                    defer { md_temp_end(temp); };

					const size_t num_aos = md_vlx_rsp_transition_density_matrix_size(vlx, nto_idx);
					const double* S = md_vlx_scf_overlap_matrix_data(vlx);

					if (num_aos > 0 && S != NULL) {
						// Full-density Mulliken attribution from attachment / detachment AO matrices.
						// hole charges  <- detachment density (D-)
						// part charges  <- attachment density (D+)
                        double* D_attach = md_temp_alloc_array(temp, double, num_aos * num_aos);
                        double* D_detach = md_temp_alloc_array(temp, double, num_aos * num_aos);

						const size_t a_dim = md_vlx_rsp_transition_density_matrix_extract(D_attach, vlx, nto_idx, MD_VLX_TRANSITION_ATTACHMENT);
						const size_t d_dim = md_vlx_rsp_transition_density_matrix_extract(D_detach, vlx, nto_idx, MD_VLX_TRANSITION_DETACHMENT);

						if (a_dim == num_aos && d_dim == num_aos) {
							double group_density_part[MAX_NTO_GROUPS] = {0};
							double group_density_hole[MAX_NTO_GROUPS] = {0};

							const int* ao_to_atom = md_vlx_ao_to_atom_idx(vlx);
                            int* ao_to_group = md_temp_alloc_array(temp, int, num_aos);

							for (size_t i = 0; i < num_aos; i++) {
								int atom_idx = ao_to_atom[i];
								ASSERT(0 <= atom_idx && (size_t)atom_idx < md_vlx_number_of_atoms(vlx));
								int g = (int)nto.atom_group_idx[atom_idx];
								if (g < 0 || g >= MAX_NTO_GROUPS) g = 0;
								ao_to_group[i] = g;
							}

							attribute_charge_density(group_density_part, ao_to_group, D_attach, S, num_aos);
							attribute_charge_density(group_density_hole, ao_to_group, D_detach, S, num_aos);

							for (size_t g1 = 0; g1 < nto.group.count; g1++) {
								// Mulliken populations of the attachment / detachment density are non-negative
								// in exact arithmetic. Clamp tiny negative numerical noise to zero so that
								// the downstream transition matrix cannot acquire spurious negative entries.
								nto.transition_density_part[g1] = (float)MAX(0.0, group_density_part[g1]);
								nto.transition_density_hole[g1] = (float)MAX(0.0, group_density_hole[g1]);
							}
							compute_transition_matrix(nto.transition_matrix, nto.group.count, nto.transition_density_hole, nto.transition_density_part);
						}
					}
				}

#else
                const float samples_per_unit_length = DEFAULT_SAMPLES_PER_ANGSTROM * BOHR_TO_ANGSTROM;

                if (use_gpu_path) {
                    md_gto_segment_and_attribute_to_groups_GPU(nto.transition_density_part, nto.group.count, nto.vol[NTO_Attachment].tex_id, &nto.grid, (const float*)nto.atom_xyzr, nto.atom_group_idx, nto.num_atoms);
                    md_gto_segment_and_attribute_to_groups_GPU(nto.transition_density_hole, nto.group.count, nto.vol[NTO_Detachment].tex_id, &nto.grid, (const float*)nto.atom_xyzr, nto.atom_group_idx, nto.num_atoms);
                    compute_transition_matrix(nto.transition_matrix, nto.group.count, nto.transition_density_hole, nto.transition_density_part);
                } else {
                    task_system::ID eval_attach = 0;
                    task_system::ID seg_attach  = 0;
                    task_system::ID eval_detach = 0;
                    task_system::ID seg_detach  = 0;

                    if (compute_transition_group_values_async(&eval_attach, &seg_attach, nto.transition_density_part, nto.group.count, nto.grid, nto.atom_group_idx, nto.atom_xyzr, nto.num_atoms, nto_idx, MD_VLX_NTO_TYPE_PARTICLE, MD_GTO_EVAL_MODE_PSI, samples_per_unit_length) &&
                        compute_transition_group_values_async(&eval_detach, &seg_detach, nto.transition_density_hole, nto.group.count, nto.grid, nto.atom_group_idx, nto.atom_xyzr, nto.num_atoms, nto_idx, MD_VLX_NTO_HOLE,     MD_GTO_EVAL_MODE_PSI, samples_per_unit_length))
                    {
                        task_system::ID compute_matrix_task = task_system::create_main_task(STR_LIT("##Compute Transition Matrix"), [nto = &nto]() {
                            compute_transition_matrix(nto->transition_matrix, nto->group.count, nto->transition_density_hole, nto->transition_density_part);
                            // Print out the groups and the values for debugging
#if 1
                            printf("Computed Transition Matrix for NTO index %d\n", nto->sel_nto_idx);
                            printf("Number of Groups: %zu\n", nto->group.count);
                            printf("Atom to Group Mapping:\n");
                            for (size_t i = 0; i < nto->num_atoms; i++) {
                                printf("    Atom %zu: Group %u\n", i, nto->atom_group_idx[i]);
                            }
                            printf("Transition Matrix:\n");
                            for (size_t i = 0; i < nto->group.count; i++) {
                                printf("  From Group %zu:\n", i);
                                for (size_t j = 0; j < nto->group.count; j++) {
                                    printf("    To Group %zu: %f\n", j, nto->transition_matrix[i * nto->group.count + j]);
                                }
                            }
                            printf("Group Densities (Particle):\n");
                            for (size_t i = 0; i < nto->group.count; i++) {
                                printf("    Group %zu: %f\n", i, nto->transition_density_part[i]);
                            }
                            printf("Group Densities (Hole):\n");
                            for (size_t i = 0; i < nto->group.count; i++) {
                                printf("    Group %zu: %f\n", i, nto->transition_density_hole[i]);
                            }
#endif          
                        });

                        task_system::set_task_dependency(compute_matrix_task, seg_attach);
                        task_system::set_task_dependency(compute_matrix_task, seg_detach);

                        task_system::enqueue_task(eval_attach);
                        task_system::enqueue_task(eval_detach);

                        nto.seg_task[0] = seg_attach;
                        nto.seg_task[1] = seg_detach;
                    } else {
                        MD_LOG_DEBUG("An error occured when computing nto group values");
                    }
                }
#endif
            }
        }
    }

};

static VeloxChem instance = {};
