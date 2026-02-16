#define IMGUI_DEFINE_MATH_OPERATORS

#include <event.h>
#include <viamd.h>
#include <task_system.h>
#include <color_utils.h>

#include <md_gto.h>
#include <md_vlx.h>
#include <md_util.h>
#include <md_topo.h>
#include <core/md_vec_math.h>
#include <core/md_log.h>
#include <core/md_arena_allocator.h>

#include <gfx/volumerender_utils.h>
#include <gfx/gl_utils.h>
#include <gfx/immediate_draw_utils.h>

#include <imgui_widgets.h>
#include <imgui_internal.h>

#include <implot_widgets.h>
#include <implot_internal.h>

#include <md_csv.h>
#include <md_xvg.h>

#include <algorithm>
#include <app/IconsFontAwesome6.h>

#define BLK_DIM 8
#define ANGSTROM_TO_BOHR 1.8897261246257702
#define BOHR_TO_ANGSTROM 0.529177210903

#define HARTREE_TO_KJ_PER_MOL 2625.4996394799

#define DEFAULT_SAMPLES_PER_ANGSTROM 8
#define DEFAULT_GTO_CUTOFF_VALUE 1.0e-6
#define MAX_NTO_GROUPS 16
#define MAX_NTO_LAMBDAS 3
#define NTO_LAMBDA_CUTOFF_VALUE 0.1

#define FORCE_CPU_PATH 0

// Resolution for broadened plots
#define NUM_SAMPLES 1024

#define IM_GREEN ImVec4(0, 1, 0, 1)
#define IM_RED ImVec4(1, 0, 0, 1)
#define IM_YELLOW ImVec4(1, 1, 0.5, 0.3)
#define IM_BLUE ImVec4(0.5, 0.5, 1, 0.3)

#define COLOR_PEAK_SELECTED ImVec4(0.5, 0.5, 1.0, 1.0)
#define COLOR_PEAK_HOVER    ImVec4(1.0, 1.0, 0.5, 1.0)

#define COLOR_TABLE_SELECTED IM_BLUE
#define COLOR_TABLE_HOVERED  IM_YELLOW

#define U32_MAGENTA IM_COL32(255, 0, 255, 255)

#define U32_VELOXCHEM_GREEN IM_COL32(0, 162, 135, 191)
#define F32_VELOXCHEM_GREEN {0, 162.0f/255.0f, 135.0f/255.0f, 0.75f}

// Complement to veloxchem green (magentaish)
#define F32_VELOXCHEM_MAGENTA {162.0f/255.0f, 35.0f/255.0f, 135.0f/255.0f, 0.75f}

static const uint64_t ATOM_PROPERTY_RESP_CHARGE = HASH_STR_LIT("ATOM PROPERTY RESP CHARGE");

// Broadening min max values
static const double gamma_min = 0.1;
static const double gamma_max = 1.0;

inline const ImVec4& vec_cast(const vec4_t& v) { return *(const ImVec4*)(&v); }
inline const vec4_t& vec_cast(const ImVec4& v) { return *(const vec4_t*)(&v); }
inline const ImVec2& vec_cast(const vec2_t& v) { return *(const ImVec2*)(&v); }
inline const vec2_t& vec_cast(const ImVec2& v) { return *(const vec2_t*)(&v); }

inline ImVec4& vec_cast(vec4_t& v) { return *(ImVec4*)(&v); }
inline vec4_t& vec_cast(ImVec4& v) { return *(vec4_t*)(&v); }
inline ImVec2& vec_cast(vec2_t& v) { return *(ImVec2*)(&v); }
inline vec2_t& vec_cast(ImVec2& v) { return *(vec2_t*)(&v); }

// This is the internal storage order of volumes and textures related to NTOs
enum NTO {
    NTO_Attachment,
    NTO_Part_0,
    NTO_Part_1,
    NTO_Part_2,
    NTO_Detachment,
    NTO_Hole_0,
    NTO_Hole_1,
    NTO_Hole_2,
};

enum class AttachmentDetachmentType {
    Attachment,
    Detachment,
};

enum x_unit_t {
    X_UNIT_EV,
    X_UNIT_NM,
    X_UNIT_CM_INVERSE,
    X_UNIT_HARTREE,
    X_UNIT_COUNT,
};

// Predefined samples per Ångström for corresponding VolumeRes (Low, Mid, High)
static const float volume_resolution_samples_per_angstrom[3] = {
    4.0f,
    8.0f,
    16.0f,
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

// Construct texture to world transformation matrix for Volume
// extent is the extent of the volume (dim * voxel_size)
static inline mat4_t compute_texture_to_world_mat(const mat3_t& orientation, const vec3_t& origin, const vec3_t& extent) {
    mat4_t T = mat4_translate_vec3(origin);
    mat4_t R = mat4_from_mat3(orientation);
    mat4_t S = mat4_scale_vec3(extent);
    return T * R * S;
}

static inline mat4_t compute_world_to_model_mat(const mat3_t& orientation, const vec3_t& origin) {
    mat4_t world_to_model = mat4_from_mat3(mat3_transpose(orientation)) * mat4_translate_vec3(-origin);
    return world_to_model;
}

static inline mat4_t compute_index_to_world_mat(const mat3_t& orientation, const vec3_t& in_origin, const vec3_t& stepsize) {
    vec3_t step_x = orientation.col[0] * stepsize.x;
    vec3_t step_y = orientation.col[1] * stepsize.y;
    vec3_t step_z = orientation.col[2] * stepsize.z;
    // Shift origin by half voxel
    vec3_t origin = in_origin + orientation * (stepsize * 0.5f);

    mat4_t index_to_world = {
        step_x.x, step_x.y, step_x.z, 0.0f,
        step_y.x, step_y.y, step_y.z, 0.0f,
        step_z.x, step_z.y, step_z.z, 0.0f,
        origin.x, origin.y, origin.z, 1.0f,
    };

    return index_to_world;
}

// Attempts to compute fitting volume dimensions given an input extent and a suggested number of samples per length unit
static inline void compute_dim(int out_dim[3], const vec3_t& in_ext, float samples_per_unit_length) {
    out_dim[0] = CLAMP(ALIGN_TO((int)(in_ext.x * samples_per_unit_length), 8), 8, 512);
    out_dim[1] = CLAMP(ALIGN_TO((int)(in_ext.y * samples_per_unit_length), 8), 8, 512);
    out_dim[2] = CLAMP(ALIGN_TO((int)(in_ext.z * samples_per_unit_length), 8), 8, 512);
}

// Attribute charge contributions encoded in the ao coefficients and overlap matrix to atoms or groups
// a: ao coefficients of NTO (num_lambdas x num_ao) array
// S: overlap matrix (num_ao x num_ao) array (row-major)
static void attribute_charge(double out_charge[], const int* ao_to_idx, const double* lambdas, const double** a, const double* S, size_t num_lambdas, size_t num_ao) {

    // Temporary buffer for r[mu] = sum_nu S[mu,nu] * a_k[nu]
	size_t temp_pos = md_temp_get_pos();
    double *r = (double*)md_temp_push(num_ao * sizeof(double));
    ASSERT(r);

    for (int k = 0; k < num_lambdas; ++k) {
        double lam = lambdas[k];

        const double* ak = a[k];

		// Initialize r to zero
		MEMSET(r, 0, num_ao * sizeof(double));

        // r = S * ak  (row-major S)
        for (int mu = 0; mu < num_ao; ++mu) {
            double a_mu = ak[mu];
            const double *S_row = &S[mu * num_ao];
            // Diagonal term
            r[mu] += S_row[mu] * a_mu;

            // Upper-triangular terms
            for (int nu = mu + 1; nu < num_ao; ++nu) {
                double s = S_row[nu];
                r[mu] += s * ak[nu];
                r[nu] += s * a_mu;
            }
        }

        // accumulate group contributions
        for (int mu = 0; mu < num_ao; ++mu) {
            int idx = ao_to_idx[mu];
            out_charge[idx] += lam * ak[mu] * r[mu];
        }
    }

	md_temp_set_pos_back(temp_pos);
}

// Voronoi segmentation
static void grid_segment_and_attribute_to_point(double out_point_value[], const vec4_t point_xyzr[], size_t num_points, const float* grid_values, const md_grid_t& grid) {
    float step_x[3] = {
        grid.orientation[0][0] * grid.spacing[0],
        grid.orientation[0][1] * grid.spacing[0],
        grid.orientation[0][2] * grid.spacing[0],
    };

    float step_y[3] = {
        grid.orientation[1][0] * grid.spacing[1],
        grid.orientation[1][1] * grid.spacing[1],
        grid.orientation[1][2] * grid.spacing[1],
    };

    float step_z[3] = {
        grid.orientation[2][0] * grid.spacing[2],
        grid.orientation[2][1] * grid.spacing[2],
        grid.orientation[2][2] * grid.spacing[2],
    };

    mat4_t index_to_world = compute_index_to_world_mat(grid.orientation, grid.origin, grid.spacing);

    for (int iz = 0; iz < grid.dim[2]; ++iz) {
        for (int iy = 0; iy < grid.dim[1]; ++iy) {
            for (int ix = 0; ix < grid.dim[0]; ++ix) {
                int index = ix + iy * grid.dim[0] + iz * grid.dim[0] * grid.dim[1];
                float value = grid_values[index];

                // Skip if its does not contribute
                if (value == 0.0f) continue;

                vec4_t coord = index_to_world * vec4_set((float)ix, (float)iy, (float)iz, 1.0f);
                coord.w = 0.0f;

                float  min_dist = FLT_MAX;
                size_t point_idx = 0;

                // find closest point to grid point
                for (size_t i = 0; i < num_points; ++i) {
                    vec4_t point = point_xyzr[i];
                    float r = point.w;
                    point.w = 0.0f;

                    float dist = vec4_distance_squared(coord, point) - r * r;
                    if (dist < min_dist) {
                        min_dist = dist;
                        point_idx = i;
                    }
                }

                out_point_value[point_idx] += value;
            }
        }
    }
}

struct VeloxChem : viamd::EventHandler {
    VeloxChem() { viamd::event_system_register_handler(*this); }
    md_vlx_t* vlx = nullptr;

    bool use_gpu_path = false;

    // Used for clearing volumes
    uint32_t vol_fbo = 0;

    // GL representations
    //md_gl_mol_t gl_mol = {};
    md_gl_rep_t gl_rep = {};

    int homo_idx[2] = {};
    int lumo_idx[2] = {};

    // Use same OBB to align all volumes
    struct OBB {
        mat3_t orientation = mat3_ident();
        vec3_t min_ext = {0};
        vec3_t max_ext = {0};
    } obb;

    struct AABB {
        vec3_t min_ext = {0};
        vec3_t max_ext = {0};
    } aabb;

    struct Summary {
        bool show_window = false;
    } summary;

    struct CriticalPoints {
        bool enabled = false;
        uint64_t vol_hash = 0;
        uint64_t raw_hash = 0;
        uint64_t simp_hash = 0;
        md_topo_extremum_graph_t raw_graph = { 0 };
        md_topo_extremum_graph_t simp_graph = { 0 };
        Volume density_vol = {0};
        md_grid_t grid = {0};
    } critical_points;

    struct Opt {
        int hovered = -1;
        int selected = -1;

        bool coord_modified = false;
    } opt;

    struct Orb {
        bool show_window = false;
        Volume   vol[16] = {};
        int      vol_mo_idx[16] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
        md_vlx_mo_type_t vol_mo_type[16] = {};
        uint32_t iso_tex[16] = {};
        task_system::ID vol_task[16] = {};
        int num_x  = 3;
        int num_y  = 3;
        int mo_idx = -1;
        int scroll_to_idx = -1;

        IsoDesc iso = {
            .enabled = true,
            .count = 2,
            .values = {0.05f, -0.05},
            .colors = {{0.f/255.f,75.f/255.f,135.f/255.f,0.75f}, {255.f/255.f,205.f/255.f,0.f/255.f,0.75f}},
        };

        GBuffer gbuf = {};
        Camera camera = {};

        struct {
            quat_t ori = {};
            vec3_t pos = {};
            float dist = {};
        } target;

        float distance_scale = 2.0f;

        bool show_coordinate_system_widget = true;
    } orb;

    struct Nto {
        bool show_window = false;
        // We have a maximum of 4 orbital slots for each particle and hole
        // A Maximum of 8 will be used in practice, 2 for attachment / detachment + 3*2 for individual lambdas (particle + hole)
        Volume   vol[8] = {};
        uint32_t iso_tex[8] = {};
        task_system::ID vol_task[8] = {};
        task_system::ID seg_task[2] = {};
        int sel_nto_idx = -1;

        md_grid_t grid = {0};

        size_t num_atoms = 0;
        vec4_t*   atom_xyzr = nullptr;      // in Atomic Units (Bohr)
        uint32_t* atom_group_idx = nullptr;

        md_gl_rep_t gl_rep = {0};

        // Square transition matrix, should have dim == num_groups
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

        vec4_t col_den = { 1.0f, 1.0f, 1.0f, 0.75f };
        vec4_t col_att = F32_VELOXCHEM_GREEN;
        vec4_t col_det = F32_VELOXCHEM_MAGENTA;

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
        Camera camera = {};
        PickingData picking = {};

        struct {
            quat_t ori = {};
            vec3_t pos = {};
            float dist = {};
        } target;

        float distance_scale = 2.0f;

        bool link_attachment_detachment_density = false;
        bool show_coordinate_system_widget = true;
    } nto;

    struct Rsp {
        bool show_window = false;
        bool show_export_window = false;
        int hovered = -1;
        int selected = -1;

        double x_ev_samples[NUM_SAMPLES] = {};
        double x_unit_samples[NUM_SAMPLES] = {};
        double* x_unit_peaks = nullptr;

        //Spectra y values, calculated from osc
        double eps[NUM_SAMPLES] = {};
        //Spectra y values, calculated from rot
        double ecd[NUM_SAMPLES] = {};

        double y_axis_sync_factor_eps = 1.0;
        double y_axis_sync_factor_ecd = 1.0;

        // Broadening gamma parameter
        double broadening_gamma = 0.123;
        broadening_mode_t broadening_mode = BROADENING_MODE_LORENTZIAN;
        x_unit_t x_unit = X_UNIT_EV;

        bool first_plot_rot_ecd = true;

    } rsp;

    struct Vib {
        int hovered = -1;
        int selected = -1;

        float gamma = 5.0f;
        broadening_mode_t broadening_mode = BROADENING_MODE_LORENTZIAN;
        bool coord_modified = false;
        float displacement_amp_scl  = 1.0f;
        float displacement_freq_scl = 1.0f;
        bool invert_x = true;
        bool invert_y = false;

        double x_samples[NUM_SAMPLES] = {}; // in cm^-1
        double y_samples[NUM_SAMPLES] = {};

        // Time accumulation for vibrational mode visualization (pertubation of atoms)
        double t = 0;

        // Scaling applied to input frequencies (to account for inaccuracies of the chosen basis set)
        double freq_scaling_factor = 1.0;

        // x_peaks and y_peaks are fetched directly from vlx
        bool first_plot = true;
    } vib;

    struct Export {
        char export_path[1024] = {};
        ElectronicStructureType orb_type = ElectronicStructureType::MolecularOrbital;
        VolumeResolution resolution = VolumeResolution::Mid;

        struct {
            md_vlx_mo_type_t type = MD_VLX_MO_TYPE_ALPHA;
            int idx = 0;
        } mo;

        struct {
            md_vlx_nto_type_t type = MD_VLX_NTO_TYPE_PARTICLE;
            int lambda_idx = 0;
            int idx = 1;
        } nto;

        bool use_obb = true;
        bool show_window = false;
    } export_state;

    // Arena for persistent allocations for the veloxchem module (tied to the lifetime of the VLX object)
    md_allocator_i* arena = 0;

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
                int gl_major, gl_minor;
                glGetIntegerv(GL_MAJOR_VERSION, &gl_major);
                glGetIntegerv(GL_MINOR_VERSION, &gl_minor);
                use_gpu_path = (!FORCE_CPU_PATH && gl_major >= 4 && gl_minor >= 3);

                break;
            }
            case viamd::EventType_ViamdShutdown:
                md_arena_allocator_destroy(arena);
                break;
            case viamd::EventType_ViamdFrameTick: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState& state = *(ApplicationState*)e.payload;

                if (vlx) {
                    draw_orb_window(state);
                    draw_summary_window(state);
                    draw_rsp_window(state);
                    if (md_vlx_rsp_has_nto(vlx)) {
                        draw_nto_window(state);
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
                        }
                        ImGui::Checkbox("Export", &export_state.show_window);
                        ImGui::EndMenu();
                    }
                }
                break;
            case viamd::EventType_ViamdRenderTransparent: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
				const ApplicationState& state = *(ApplicationState*)e.payload;

                if (critical_points.enabled && critical_points.simp_graph.num_vertices > 0) {

				    glBindFramebuffer(GL_FRAMEBUFFER, state.gbuffer.fbo);
                    // Render topology as points if available
				    immediate::set_model_view_matrix(state.view.param.matrix.curr.view);
				    immediate::set_proj_matrix(state.view.param.matrix.curr.proj);

                    glDisable(GL_DEPTH_TEST);

					uint32_t maxima_offset = md_topo_offset_maxima(&critical_points.simp_graph);
					uint32_t maxima_count  = critical_points.simp_graph.num_maxima;
                    for (size_t i = maxima_offset; i < maxima_offset + maxima_count; ++i) {
                        vec3_t pos = {critical_points.simp_graph.vertices[i].x, critical_points.simp_graph.vertices[i].y, critical_points.simp_graph.vertices[i].z};
                        pos *= BOHR_TO_ANGSTROM;
                        immediate::draw_point(pos, immediate::COLOR_RED);
                    }

					uint32_t splits_offset = md_topo_offset_split_saddles(&critical_points.simp_graph);
					uint32_t splits_count = critical_points.simp_graph.num_split_saddles;
                    for (size_t i = splits_offset; i < splits_offset + splits_count; ++i) {
                        vec3_t pos = {critical_points.simp_graph.vertices[i].x, critical_points.simp_graph.vertices[i].y, critical_points.simp_graph.vertices[i].z};
                        pos *= BOHR_TO_ANGSTROM;
                        immediate::draw_point(pos, immediate::COLOR_GREEN);
					}

					uint32_t minima_offset = md_topo_offset_minima(&critical_points.simp_graph);
					uint32_t minima_count = critical_points.simp_graph.num_minima;
                    for (size_t i = minima_offset; i < minima_offset + minima_count; ++i) {
                        vec3_t pos = {critical_points.simp_graph.vertices[i].x, critical_points.simp_graph.vertices[i].y, critical_points.simp_graph.vertices[i].z};
                        pos *= BOHR_TO_ANGSTROM;
						immediate::draw_point(pos, immediate::COLOR_BLUE);
					}

					uint32_t join_saddles_offset = md_topo_offset_join_saddles(&critical_points.simp_graph);
					uint32_t join_saddles_count = critical_points.simp_graph.num_join_saddles;
                    for (size_t i = join_saddles_offset; i < join_saddles_offset + join_saddles_count; ++i) {
                        vec3_t pos = {critical_points.simp_graph.vertices[i].x, critical_points.simp_graph.vertices[i].y, critical_points.simp_graph.vertices[i].z};
                        pos *= BOHR_TO_ANGSTROM;
                        immediate::draw_point(pos, immediate::COLOR_CYAN);
					}

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
						immediate::draw_line(p0, p1, immediate::COLOR_BLACK);
                    }

                    immediate::render();
                }

                //ApplicationState& state = *(ApplicationState*)e.payload;
                //draw_orb_volume(state);
                break;
            }
            case viamd::EventType_ViamdTopologyInit: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState& state = *(ApplicationState*)e.payload;
                init_from_file(str_from_cstr(state.files.molecule), state);
                break;
            }
            case viamd::EventType_ViamdTopologyFree:
                reset_data();
                break;

            case viamd::EventType_ViamdRepresentationInfoFill: {
                ASSERT(e.payload_type == viamd::EventPayloadType_RepresentationInfo);
                RepresentationInfo& info = *(RepresentationInfo*)e.payload;

                info.alpha.homo_idx = homo_idx[0];
                info.alpha.lumo_idx = lumo_idx[0];

                info.beta.homo_idx  = homo_idx[1];
                info.beta.lumo_idx  = lumo_idx[1];

                size_t num_mos = num_molecular_orbitals();
                if (num_mos) {
                    const double* occ = md_vlx_scf_mo_occupancy(vlx, MD_VLX_MO_TYPE_ALPHA);
                    const double* ene = md_vlx_scf_mo_energy   (vlx, MD_VLX_MO_TYPE_ALPHA);

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
                if (false) {
                    const double* occ = md_vlx_scf_mo_occupancy(vlx, MD_VLX_MO_TYPE_BETA);
                    const double* ene = md_vlx_scf_mo_energy   (vlx, MD_VLX_MO_TYPE_BETA);

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
                        const double* lambdas = md_vlx_rsp_nto_lambdas(vlx, i);
                        NaturalTransitionOrbitalLambda lambda_info = {};
                        for (size_t j = 0; j < 16; ++j) {
                            if (lambdas[j] < LAMBDA_CUTOFF) break;
                            str_t lbl = str_printf(info.alloc, (const char*)u8"λ[%zu] (%.3f)", j + 1, lambdas[j]);
                            md_array_push(lambda_info.label, lbl, info.alloc);
                            md_array_push(lambda_info.value, lambdas[j], info.alloc);
                            lambda_info.num_lambdas += 1;
                        }
                        md_array_push(info.nto.lambda, lambda_info, info.alloc);
                    }
                }
                
                if (md_vlx_scf_number_of_molecular_orbitals(vlx) > 0) {
                    info.electronic_structure_type_mask |= (1 << (int)ElectronicStructureType::MolecularOrbital);
                    info.electronic_structure_type_mask |= (1 << (int)ElectronicStructureType::MolecularOrbitalDensity);
                    info.electronic_structure_type_mask |= (1 << (int)ElectronicStructureType::ElectronDensity);
                }

                if (md_vlx_rsp_number_of_excited_states(vlx) > 0) {
                    info.electronic_structure_type_mask |= (1 << (int)ElectronicStructureType::NaturalTransitionOrbitalParticle);
                    info.electronic_structure_type_mask |= (1 << (int)ElectronicStructureType::NaturalTransitionOrbitalHole);
                    info.electronic_structure_type_mask |= (1 << (int)ElectronicStructureType::NaturalTransitionOrbitalDensityParticle);
                    info.electronic_structure_type_mask |= (1 << (int)ElectronicStructureType::NaturalTransitionOrbitalDensityHole);
                    info.electronic_structure_type_mask |= (1 << (int)ElectronicStructureType::AttachmentDensity);
                    info.electronic_structure_type_mask |= (1 << (int)ElectronicStructureType::DetachmentDensity);
                }

                const double* resp_charges = md_vlx_scf_resp_charges(vlx);
                if (resp_charges) {
                    double value_min =  DBL_MAX;
                    double value_max = -DBL_MAX;

                    for (size_t i = 0; i < md_vlx_number_of_atoms(vlx); ++i) {
                        value_min = MIN(value_min, resp_charges[i]);
                        value_max = MAX(value_max, resp_charges[i]);
                    }

                    AtomProperty prop = {
                        .id = ATOM_PROPERTY_RESP_CHARGE,
                        .label = STR_LIT("Resp Charge"),
                        .num_idx = 0,
                        .value_min = (float)value_min,
                        .value_max = (float)value_max,
                    };

                    md_array_push(info.atom_properties, prop, info.alloc);
                }

                // @TODO: Fill in dipole information
                break;
            }
            case viamd::EventType_ViamdRepresentationEvalElectronicStructure: {
                ASSERT(e.payload_type == viamd::EventPayloadType_EvalElectronicStructure);
                EvalElectronicStructure& data = *(EvalElectronicStructure*)e.payload;

                if (!data.output_written) {
                    const float samples_per_unit_length = data.samples_per_angstrom * BOHR_TO_ANGSTROM;
                    md_grid_t grid = {};
                    init_grid(&grid, obb.orientation, obb.min_ext, obb.max_ext, samples_per_unit_length);
                    init_volume(data.dst_volume, grid, GL_R32F);

                    switch (data.type) {
                    case ElectronicStructureType::MolecularOrbital:
                    case ElectronicStructureType::MolecularOrbitalDensity:
                    {
                        md_gto_eval_mode_t mode = (data.type == ElectronicStructureType::MolecularOrbital) ? MD_GTO_EVAL_MODE_PSI : MD_GTO_EVAL_MODE_PSI_SQUARED;
                        if (use_gpu_path) {
                            data.output_written = compute_mo_GPU(data.dst_volume->tex_id, grid, MD_VLX_MO_TYPE_ALPHA, data.major_idx, mode);
                        }
                        else {
                            data.output_written = (compute_mo_async(data.dst_volume->tex_id, grid, MD_VLX_MO_TYPE_ALPHA, data.major_idx, mode) != task_system::INVALID_ID);
                        }
                        break;
                    }
                    case ElectronicStructureType::NaturalTransitionOrbitalParticle:
                    case ElectronicStructureType::NaturalTransitionOrbitalHole:
                    case ElectronicStructureType::NaturalTransitionOrbitalDensityParticle:
                    case ElectronicStructureType::NaturalTransitionOrbitalDensityHole:
                    {
                        md_vlx_nto_type_t type  = (data.type == ElectronicStructureType::NaturalTransitionOrbitalParticle ||
                                                   data.type == ElectronicStructureType::NaturalTransitionOrbitalDensityParticle)
                                                   ? MD_VLX_NTO_TYPE_PARTICLE : MD_VLX_NTO_TYPE_HOLE;
                        md_gto_eval_mode_t mode = (data.type == ElectronicStructureType::NaturalTransitionOrbitalParticle ||
                                                   data.type == ElectronicStructureType::NaturalTransitionOrbitalHole)
                                                   ? MD_GTO_EVAL_MODE_PSI : MD_GTO_EVAL_MODE_PSI_SQUARED;

                        if (use_gpu_path) {
                            data.output_written = compute_nto_GPU(data.dst_volume->tex_id, grid, data.major_idx, data.minor_idx, type, mode);
                        } else {
                            data.output_written = (compute_nto_async(data.dst_volume->tex_id, grid, data.major_idx, data.minor_idx, type, mode) != task_system::INVALID_ID);
                        }
                        break;
                    }
                    case ElectronicStructureType::AttachmentDensity:
                    case ElectronicStructureType::DetachmentDensity:
                    {
                        AttachmentDetachmentType nto_type = (data.type == ElectronicStructureType::AttachmentDensity) ? AttachmentDetachmentType::Attachment : AttachmentDetachmentType::Detachment;
                        if (use_gpu_path) {
                            data.output_written = compute_attachment_detachment_density_GPU(data.dst_volume->tex_id, grid, data.major_idx, nto_type);
                        } else {
                            data.output_written = (compute_attachment_detachment_density_async(data.dst_volume->tex_id, grid, data.major_idx, nto_type) != task_system::INVALID_ID);
                        }
                        break;
                    }
                    case ElectronicStructureType::ElectronDensity:
                    {
                        if (use_gpu_path) {
                            data.output_written = compute_electron_density_GPU(data.dst_volume->tex_id, grid, MD_VLX_MO_TYPE_ALPHA);
                        } else {
                            data.output_written = (compute_electron_density_async(data.dst_volume->tex_id, grid, MD_VLX_MO_TYPE_ALPHA) != task_system::INVALID_ID);
                        }
                        break;
                    }
                    default:
                        MD_LOG_ERROR("Invalid Orbital Type supplied to Compute Orbital Event");
                        break;
                    }
                }

                break;
            }
            case viamd::EventType_ViamdRepresentationEvalAtomProperty: {
                ASSERT(e.payload_type == viamd::EventPayloadType_EvalAtomProperty);
                EvalAtomProperty& data = *(EvalAtomProperty*)e.payload;

                switch (data.property_id) {
                case ATOM_PROPERTY_RESP_CHARGE: {
                    size_t num_atoms = md_vlx_number_of_atoms(vlx);
                    if (data.num_values != num_atoms) {
                        MD_LOG_ERROR("Invalid number of values, did not match number of atoms");
                        break;
                    }

                    if (data.dst_values == nullptr) {
                        MD_LOG_ERROR("Dst values not supplied");
                        break;
                    }

                    const double* resp_charges = md_vlx_scf_resp_charges(vlx);
                    if (resp_charges) {
                        for (size_t i = 0; i < num_atoms; ++i) {
                            data.dst_values[i] = (float)resp_charges[i];
                        }
                        data.output_written = true;
                    }
                    break;
                }
                default:
                    // Do not fail here, as conceptually there may be other components that can calculate this property
                    break;
                }

                break;
            }
            default:
                break;
            }
        }
    }

    void reset_data() {
        //md_gl_mol_destroy(gl_mol);
        md_gl_rep_destroy(gl_rep);
        md_vlx_destroy(vlx);
        md_arena_allocator_reset(arena);
        gl::free_texture(&critical_points.density_vol.tex_id);
        md_topo_extremum_graph_free(&critical_points.raw_graph);
        md_topo_extremum_graph_free(&critical_points.simp_graph);
        critical_points = {};
        critical_points.raw_graph.alloc = arena;
        critical_points.simp_graph.alloc = arena;
        vlx = nullptr;
        orb = VeloxChem::Orb{};
        nto = VeloxChem::Nto{};
        rsp = VeloxChem::Rsp{};
        vib = VeloxChem::Vib{};
        opt = VeloxChem::Opt{};
    }

    void update_nto_group_colors() {
        ScopedTemp temp_reset;
        uint32_t* colors = (uint32_t*)md_temp_push(sizeof(uint32_t) * nto.num_atoms);
        for (size_t i = 0; i < nto.num_atoms; ++i) {
            const size_t group_idx = nto.atom_group_idx[i];
            const uint32_t color = group_idx < nto.group.count ? convert_color(nto.group.color[group_idx]) : U32_MAGENTA;
            colors[i] = color;
        }

        md_gl_rep_set_color(nto.gl_rep, 0, (uint32_t)nto.num_atoms, colors, 0);
    }

    void init_from_file(str_t filename, ApplicationState& state) {
        str_t ext;
        if (extract_ext(&ext, filename)) {
            if (str_eq_ignore_case(ext, STR_LIT("out")) || str_eq_ignore_case(ext, STR_LIT("h5"))) {
                MD_LOG_INFO("Attempting to load VeloxChem data from file '" STR_FMT "'", STR_ARG(filename));
                
                if (!vlx) {
                    vlx = md_vlx_create(arena);
                } else {
                    md_vlx_reset(vlx);
                }
                
                if (md_vlx_parse_file(vlx, filename)) {
                    MD_LOG_INFO("Successfully loaded VeloxChem data");

                    if (!vol_fbo) glGenFramebuffers(1, &vol_fbo);

                    // Scf
                    //scf.show_window = true;

                    homo_idx[0] = (int)md_vlx_scf_homo_idx(vlx, MD_VLX_MO_TYPE_ALPHA);
                    homo_idx[1] = (int)md_vlx_scf_homo_idx(vlx, MD_VLX_MO_TYPE_BETA);

                    lumo_idx[0] = (int)md_vlx_scf_lumo_idx(vlx, MD_VLX_MO_TYPE_ALPHA);
                    lumo_idx[1] = (int)md_vlx_scf_lumo_idx(vlx, MD_VLX_MO_TYPE_BETA);

                    size_t num_atoms = md_vlx_number_of_atoms(vlx);
                    const dvec3_t* coords = md_vlx_atom_coordinates(vlx);
                    const uint8_t* atomic_numbers = md_vlx_atomic_numbers(vlx);

                    nto.atom_xyzr = (vec4_t*)md_arena_allocator_push(arena, sizeof(vec4_t) * num_atoms);
                    nto.num_atoms = num_atoms;

                    // Compute the PCA of the provided geometry
                    // This is used in determining a better fitting volume for the orbitals
                    vec4_t* xyzw = (vec4_t*)md_vm_arena_push(state.allocator.frame, sizeof(vec4_t) * num_atoms);
                    for (size_t i = 0; i < num_atoms; ++i) {
                        nto.atom_xyzr[i] = vec4_set((float)coords[i].x, (float)coords[i].y, (float)coords[i].z, md_util_element_vdw_radius(atomic_numbers[i])) * ANGSTROM_TO_BOHR;
                        xyzw[i] = vec4_set((float)coords[i].x, (float)coords[i].y, (float)coords[i].z, 1.0f);
                    }

                    md_system_t mol = { 0 };
                    md_vlx_system_init(&mol, vlx, state.allocator.frame);
                    md_util_molecule_postprocess(&mol, state.allocator.frame, MD_UTIL_POSTPROCESS_BOND_BIT | MD_UTIL_POSTPROCESS_STRUCTURE_BIT);
                    //gl_mol = md_gl_mol_create(&mol);

                    uint32_t* colors = (uint32_t*)md_vm_arena_push(state.allocator.frame, mol.atom.count * sizeof(uint32_t));
                    color_atoms_type(colors, mol.atom.count, mol);

                    gl_rep = md_gl_rep_create(state.mold.gl_mol);
                    md_gl_rep_set_color(gl_rep, 0, (uint32_t)mol.atom.count, colors, 0);

                    vec3_t com = md_util_com_compute_vec4(xyzw, 0, num_atoms, 0);
                    mat3_t cov = mat3_covariance_matrix_vec4(xyzw, 0, num_atoms, com);
                    mat3_eigen_t eigen = mat3_eigen(cov);
                    mat3_t PCA = mat3_orthonormalize(mat3_extract_rotation(eigen.vectors));

                    // Compute min and maximum extent along the PCA axes
                    obb.orientation = mat3_transpose(PCA);
                    obb.min_ext  = vec3_set1( FLT_MAX);
                    obb.max_ext  = vec3_set1(-FLT_MAX);
                    aabb.min_ext = vec3_set1( FLT_MAX);
                    aabb.max_ext = vec3_set1(-FLT_MAX);

                    // Transform the gto (x,y,z,cutoff) into the PCA frame to find the min and max extend within it
                    for (size_t i = 0; i < num_atoms; ++i) {
                        vec3_t xyz = vec3_from_vec4(xyzw[i]) * ANGSTROM_TO_BOHR;
                        aabb.min_ext = vec3_min(aabb.min_ext, xyz);
                        aabb.max_ext = vec3_max(aabb.max_ext, xyz);

                        xyz = mat3_mul_vec3(PCA, xyz);
                        obb.min_ext = vec3_min(obb.min_ext, xyz);
                        obb.max_ext = vec3_max(obb.max_ext, xyz);
                    }

                    // This is the extra padding we apply to the 'bounding volumes'
                    const float pad = 6.0f;
                    aabb.min_ext -= pad;
                    aabb.max_ext += pad;

                    obb.min_ext -= pad;
                    obb.max_ext += pad;

                    // NTO
                    size_t num_excited_states = md_vlx_rsp_number_of_excited_states(vlx);
                    if (num_excited_states > 0) {
                        //nto.show_window = true;
                        camera_compute_optimal_view(&nto.target.pos, &nto.target.ori, &nto.target.dist, obb.orientation, obb.min_ext * BOHR_TO_ANGSTROM, obb.max_ext * BOHR_TO_ANGSTROM, nto.distance_scale);
                        nto.atom_group_idx = (uint32_t*)md_alloc(arena, sizeof(uint32_t) * mol.atom.count);
                        MEMSET(nto.atom_group_idx, 0, sizeof(uint32_t) * mol.atom.count);

                        snprintf(nto.group.label[0], sizeof(nto.group.label[0]), "Unassigned");
                        nto.group.color[0] = vec4_t{ 0, 0, 0, 1 };

                        for (int i = 1; i < (int)ARRAY_SIZE(nto.group.color); ++i) {
                            ImVec4 color = ImPlot::GetColormapColor(i - 1, ImPlotColormap_Deep);
                            nto.group.color[i] = vec_cast(color);
                            snprintf(nto.group.label[i], sizeof(nto.group.label[i]), "Group %i", i);
                        }

                        str_t file = {};
                        extract_file(&file, filename);

                        if (str_eq_cstr(file, "tq.out")) {
                            uint32_t index_from_text[23] = { 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 };
                            //nto.atom_group_idx = index_from_text;
                            nto.group.count = 3;
                            MEMCPY(nto.atom_group_idx, index_from_text, sizeof(index_from_text));
                            snprintf(nto.group.label[1], sizeof(nto.group.label[1]), "Thio");
                            snprintf(nto.group.label[2], sizeof(nto.group.label[2]), "Quin");
                        }
                        else {

                            // @TODO: Remove once proper interface is there
                            nto.group.count = 3;
                            // Assign half of the atoms to group 1
                            for (size_t i = 0; i < mol.atom.count; ++i) {
                                nto.atom_group_idx[i] = i < mol.atom.count / 2 ? 1 : 2;
                            }
                        }
                        nto.gl_rep = md_gl_rep_create(state.mold.gl_mol);
                        update_nto_group_colors();

                        // Callculate ballpark scaling factor for dipole vectors
                        vec3_t extent = aabb.max_ext - aabb.min_ext;
                        float max_ext = MAX(extent.x, MAX(extent.y, extent.z));
                        float max_len = 0;
                        const dvec3_t* electric_dp = md_vlx_rsp_electric_transition_dipole_moments(vlx);
                        const dvec3_t* magnetic_dp = md_vlx_rsp_magnetic_transition_dipole_moments(vlx);
                        ASSERT(electric_dp);
                        ASSERT(magnetic_dp);
                        for (int i = 0; i < num_excited_states; ++i) {
                            max_len = MAX(max_len, (float)dvec3_length(electric_dp[i]));
                            max_len = MAX(max_len, (float)dvec3_length(magnetic_dp[i]));
                        }
                        nto.dipole.vector_scale = CLAMP((max_ext * 0.75f) / max_len, 0.1f, 10.0f);

                        init_grid(&nto.grid, obb.orientation, obb.min_ext, obb.max_ext, DEFAULT_SAMPLES_PER_ANGSTROM * BOHR_TO_ANGSTROM);
                    }

                    // RSP
                    if (num_excited_states > 0) {
                        //rsp.show_window = true;
                        rsp.hovered = -1;
                        rsp.selected = -1;

                        md_array_resize(rsp.x_unit_peaks, num_excited_states, arena);

                        // Populate x values
                        const double* abs_ev = md_vlx_rsp_absorption_ev(vlx);
                        if (abs_ev) {
                            const double x_min = abs_ev[0] - 1.0;
                            const double x_max = abs_ev[num_excited_states - 1] + 1.0;
                            for (int i = 0; i < NUM_SAMPLES; ++i) {
                                double t = (double)i / (double)(NUM_SAMPLES - 1);
                                double value = lerp(x_min, x_max, t);
                                rsp.x_ev_samples[i] = value;
                            }
                        }
                    }

                    // OPT
                    if (md_vlx_opt_number_of_steps(vlx) > 0) {
                        opt.selected = (int)(md_vlx_opt_number_of_steps(vlx) - 1);
                    }

                    // ORB
                    //orb.show_window = true;
                    camera_compute_optimal_view(&orb.target.pos, &orb.target.ori, &orb.target.dist, obb.orientation, obb.min_ext * BOHR_TO_ANGSTROM, obb.max_ext * BOHR_TO_ANGSTROM, orb.distance_scale);
                    orb.mo_idx = homo_idx[0];
                    orb.scroll_to_idx = homo_idx[0];

                    // Export
                    export_state.mo.idx = homo_idx[0];
                }
                else {
                    MD_LOG_INFO("Failed to load VeloxChem data");
                    reset_data();
                }
            }
        }
    }

    void init_grid(md_grid_t* grid, const mat3_t& orientation, const vec3_t& min_ext, const vec3_t& max_ext, float samples_per_unit_length) {
        vec3_t extent = max_ext - min_ext;
        compute_dim(grid->dim, extent, samples_per_unit_length);
        vec3_t voxel_size = vec3_div(extent, vec3_set((float)grid->dim[0], (float)grid->dim[1], (float)grid->dim[2]));
        grid->orientation = orientation;
        grid->origin = orientation * min_ext;
        grid->spacing = voxel_size;
    }

    void init_volume(Volume* vol, const md_grid_t& grid, GLenum format = GL_R16F) {
        ASSERT(vol);
        MEMCPY(vol->dim, grid.dim, sizeof(vol->dim));
        const float scl = BOHR_TO_ANGSTROM;

        vec3_t extent = md_grid_extent(&grid);
        vol->texture_to_world = compute_texture_to_world_mat(grid.orientation, grid.origin * scl, extent * scl);
        vol->voxel_size = grid.spacing;
        gl::init_texture_3D(&vol->tex_id, vol->dim[0], vol->dim[1], vol->dim[2], format);
    }

    bool compute_nto_GPU(uint32_t vol_tex, const md_grid_t& grid, size_t nto_idx, size_t lambda_idx, md_vlx_nto_type_t type, md_gto_eval_mode_t mode, double cutoff_value = DEFAULT_GTO_CUTOFF_VALUE) {
        ScopedTemp temp_reset;

        size_t num_gtos = md_vlx_nto_gto_count(vlx);
        md_gto_t* gtos  = (md_gto_t*)md_temp_push(sizeof(md_gto_t) * num_gtos);

        if (!md_vlx_nto_gto_extract(gtos, vlx, nto_idx, lambda_idx, type)) {
            MD_LOG_ERROR("Failed to extract NTO gto for nto index: %zu and lambda: %zu", nto_idx, lambda_idx);
            return task_system::INVALID_ID;
        }
        num_gtos = md_gto_cutoff_compute_and_filter(gtos, num_gtos, cutoff_value);

        // Dispatch Compute shader evaluation
        md_gto_grid_evaluate_GPU(vol_tex, &grid, gtos, num_gtos, mode);

        return true;
    }

    bool extract_attachment_detachment_orb_data(md_orbital_data_t* orb_data, double cutoff, AttachmentDetachmentType type, size_t nto_idx, md_allocator_i* alloc) {
        ASSERT(orb_data);

        md_vlx_nto_type_t nto_type = (type == AttachmentDetachmentType::Attachment) ? MD_VLX_NTO_TYPE_PARTICLE : MD_VLX_NTO_TYPE_HOLE;
        size_t num_gtos_per_lambda = md_vlx_nto_gto_count(vlx);

        size_t num_lambdas = 0;
        const double* lambda = md_vlx_rsp_nto_lambdas(vlx, nto_idx);
        if (lambda) {
            for (size_t i = 0; i < MAX_NTO_LAMBDAS; ++i) {
                if (lambda[i] < NTO_LAMBDA_CUTOFF_VALUE) {
                    break;
                }
                num_lambdas += 1;
            }
        }

        md_array_ensure(orb_data->gtos, num_gtos_per_lambda * num_lambdas, alloc);
        md_array_ensure(orb_data->orb_offsets, num_lambdas + 1, alloc);
        md_array_ensure(orb_data->orb_scaling, num_lambdas, alloc);

        size_t temp_pos = md_temp_get_pos();
        defer { md_temp_set_pos_back(temp_pos); };

        md_gto_t* temp_gtos = (md_gto_t*)md_temp_push(num_gtos_per_lambda * sizeof(md_gto_t));
        size_t num_temp_gtos = num_gtos_per_lambda;

        md_array_push(orb_data->orb_offsets, (uint32_t)md_array_size(orb_data->gtos), alloc);

        for (size_t i = 0; i < num_lambdas; ++i) {
            if (!md_vlx_nto_gto_extract(temp_gtos, vlx, nto_idx, i, nto_type)) {
                MD_LOG_ERROR("Failed to extract NTO gto for nto index: %zu and lambda: %zu", nto_idx, i);
                return false;
            }
            size_t num_pruned = md_gto_cutoff_compute_and_filter(temp_gtos, num_temp_gtos, cutoff);
            md_array_push_array(orb_data->gtos, temp_gtos, num_pruned, alloc);
            md_array_push(orb_data->orb_offsets, (uint32_t)md_array_size(orb_data->gtos), alloc);
            md_array_push(orb_data->orb_scaling, (float)lambda[i], alloc);
        }

        orb_data->num_gtos = md_array_size(orb_data->gtos);
        orb_data->num_orbs = md_array_size(orb_data->orb_scaling);

        return true;
    }

    bool extract_electron_density_orb_data(md_orbital_data_t* orb_data, double cutoff, md_allocator_i* alloc) {
        size_t num_mo = MAX(md_vlx_scf_lumo_idx(vlx, MD_VLX_MO_TYPE_ALPHA), md_vlx_scf_lumo_idx(vlx, MD_VLX_MO_TYPE_BETA));
        size_t num_gtos_per_mo = md_vlx_mo_gto_count(vlx);

        const double* occ_a = md_vlx_scf_mo_occupancy(vlx, MD_VLX_MO_TYPE_ALPHA);
        const double* occ_b = md_vlx_scf_mo_occupancy(vlx, MD_VLX_MO_TYPE_BETA);
        md_vlx_scf_type_t scf_type = md_vlx_scf_type(vlx);

        md_array_ensure(orb_data->gtos, num_gtos_per_mo * num_mo, alloc);
        md_array_ensure(orb_data->orb_offsets, num_mo + 1, alloc);
        md_array_ensure(orb_data->orb_scaling, num_mo, alloc);

        size_t temp_pos = md_temp_get_pos();
        defer { md_temp_set_pos_back(temp_pos); };

        md_gto_t* temp_gtos = (md_gto_t*)md_temp_push(num_gtos_per_mo * sizeof(md_gto_t));
        size_t num_temp_gtos = num_gtos_per_mo;

        md_array_push(orb_data->orb_offsets, (uint32_t)md_array_size(orb_data->gtos), alloc);

        if (scf_type == MD_VLX_SCF_TYPE_RESTRICTED || scf_type == MD_VLX_SCF_TYPE_RESTRICTED_OPENSHELL) {
            for (size_t mo_idx = 0; mo_idx < num_mo; ++mo_idx) {
                double occ = occ_a[mo_idx] + occ_b[mo_idx];
                size_t num_gtos = md_vlx_mo_gto_extract(temp_gtos, vlx, mo_idx, MD_VLX_MO_TYPE_ALPHA, cutoff);
                md_array_push_array(orb_data->gtos, temp_gtos, num_gtos, alloc);
                md_array_push(orb_data->orb_offsets, (uint32_t)md_array_size(orb_data->gtos), alloc);
                md_array_push(orb_data->orb_scaling, (float)occ, alloc);
            }
        }
        else if (scf_type == MD_VLX_SCF_TYPE_UNRESTRICTED) {
            for (size_t mo_idx = 0; mo_idx < num_mo; ++mo_idx) {
                if (occ_a[mo_idx] == 0.0 && occ_b[mo_idx] == 0.0) {
                    break;
                }

                if (occ_a[mo_idx] > 0.0) {
                    double occ = occ_a[mo_idx];
                    size_t num_gtos = md_vlx_mo_gto_extract(temp_gtos, vlx, mo_idx, MD_VLX_MO_TYPE_ALPHA, cutoff);
                    md_array_push_array(orb_data->gtos, temp_gtos, num_gtos, alloc);
                    md_array_push(orb_data->orb_offsets, (uint32_t)md_array_size(orb_data->gtos), alloc);
                    md_array_push(orb_data->orb_scaling, (float)occ, alloc);
                }

                if (occ_b[mo_idx] > 0.0) {
                    double occ = occ_b[mo_idx];
                    size_t num_gtos = md_vlx_mo_gto_extract(temp_gtos, vlx, mo_idx, MD_VLX_MO_TYPE_BETA, cutoff);
                    md_array_push_array(orb_data->gtos, temp_gtos, num_gtos, alloc);
                    md_array_push(orb_data->orb_offsets, (uint32_t)md_array_size(orb_data->gtos), alloc);
                    md_array_push(orb_data->orb_scaling, (float)occ, alloc);
                }
            }
        }

        orb_data->num_gtos = md_array_size(orb_data->gtos);
        orb_data->num_orbs = md_array_size(orb_data->orb_scaling);

        return true;
    }

    // Compute attachment / detachment density
    bool compute_attachment_detachment_density_GPU(uint32_t vol_tex, const md_grid_t& grid, size_t nto_idx, AttachmentDetachmentType type, double cutoff_value = DEFAULT_GTO_CUTOFF_VALUE) {
        md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(64));
        defer { md_vm_arena_destroy(alloc); };

        md_orbital_data_t orb_data = {0};
        if (!extract_attachment_detachment_orb_data(&orb_data, cutoff_value, type, nto_idx, alloc)) {
            return false;
        }

        // Dispatch Compute shader evaluation
        md_gto_grid_evaluate_orb_GPU(vol_tex, &grid, &orb_data, MD_GTO_EVAL_MODE_PSI_SQUARED);

        return true;
    }

    // Compute attachment / detachment density
    task_system::ID compute_attachment_detachment_density_async(uint32_t vol_tex, const md_grid_t& grid, size_t nto_idx, AttachmentDetachmentType type, double cutoff_value = DEFAULT_GTO_CUTOFF_VALUE) {
        md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

        md_orbital_data_t orb_data = {0};
        if (!extract_attachment_detachment_orb_data(&orb_data, cutoff_value, type, nto_idx, alloc)) {
            md_vm_arena_destroy(alloc);
            return task_system::INVALID_ID;
        }

        struct Payload {
            AsyncGridEvalArgs args;
            md_allocator_i* arena;
            uint32_t tex;
        };

        // Allocate data for payload struct + volume
        Payload* payload = (Payload*)md_vm_arena_push(alloc, sizeof(Payload));
        *payload = {
            .args = {
                .grid = grid,
                .grid_data = (float*)md_vm_arena_push_zero(alloc, sizeof(float) * grid.dim[0] * grid.dim[1] * grid.dim[2]),
                .orb  = orb_data,
                .mode = MD_GTO_EVAL_MODE_PSI_SQUARED,
            },
            .arena = alloc,
            .tex = vol_tex,
        };

        task_system::ID async_task = evaluate_gto_on_grid_async(&payload->args);

        // Launch task for main (render) thread to update the volume texture
        task_system::ID main_task = task_system::create_main_task(STR_LIT("##Update Volume"), [data = payload]() {
            // Ensure that the dimensions of the texture have not changed during evaluation
            int dim[3];
            if (gl::get_texture_dim(dim, data->tex) && MEMCMP(dim, data->args.grid.dim, sizeof(dim)) == 0) {
                gl::set_texture_3D_data(data->tex, data->args.grid_data, GL_R32F);
            }

            md_vm_arena_destroy(data->arena);
        });

        task_system::set_task_dependency(main_task, async_task);
        task_system::enqueue_task(async_task);

        return async_task;
    }

    task_system::ID compute_nto_async(uint32_t vol_tex, const md_grid_t& grid, size_t nto_idx, size_t lambda_idx, md_vlx_nto_type_t type, md_gto_eval_mode_t mode, double cutoff_value = DEFAULT_GTO_CUTOFF_VALUE) {
        md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

        md_orbital_data_t orb_data = {0};
        orb_data.num_gtos = md_vlx_nto_gto_count(vlx);
        orb_data.gtos = (md_gto_t*)md_vm_arena_push(alloc, sizeof(md_gto_t) * orb_data.num_gtos);
        if (!md_vlx_nto_gto_extract(orb_data.gtos, vlx, nto_idx, lambda_idx, type)) {
            MD_LOG_ERROR("Failed to extract NTO gto for nto index: %zu and lambda: %zu", nto_idx, lambda_idx);
            md_vm_arena_destroy(alloc);
            return task_system::INVALID_ID;
        }
        orb_data.num_gtos = md_gto_cutoff_compute_and_filter(orb_data.gtos, orb_data.num_gtos, cutoff_value);

        struct Payload {
            AsyncGridEvalArgs args;
            md_allocator_i* alloc;
            uint32_t tex;
        };

        Payload* payload = (Payload*)md_vm_arena_push(alloc, sizeof(Payload));
        *payload = {
            .args = {
                .grid = grid,
                .grid_data = (float*)md_vm_arena_push_zero(alloc, sizeof(float) * grid.dim[0] * grid.dim[1] * grid.dim[2]),
                .orb = orb_data,
                .mode = mode,
            },
            .alloc = alloc,
            .tex = vol_tex,
        };

        task_system::ID async_task = evaluate_gto_on_grid_async(&payload->args);

        // Launch task for main (render) thread to update the volume texture
        task_system::ID main_task = task_system::create_main_task(STR_LIT("##Update Volume"), [data = payload]() {
            // Ensure that the dimensions of the texture have not changed during evaluation
            int dim[3];
            if (gl::get_texture_dim(dim, data->tex) && MEMCMP(dim, data->args.grid.dim, sizeof(dim)) == 0) {
                gl::set_texture_3D_data(data->tex, data->args.grid_data, GL_R32F);
            }

            md_vm_arena_destroy(data->alloc);
        });

        task_system::set_task_dependency(main_task, async_task);
        task_system::enqueue_task(async_task);

        return async_task;
    }

    struct AsyncGridEvalArgs {
        md_grid_t  grid;
        float*     grid_data;
        md_orbital_data_t orb;
        md_gto_eval_mode_t mode;
    };

    bool compute_mo_GPU(uint32_t vol_tex, const md_grid_t& grid, md_vlx_mo_type_t mo_type, size_t mo_idx, md_gto_eval_mode_t mode, double cutoff_value = DEFAULT_GTO_CUTOFF_VALUE) {
        ScopedTemp reset_temp;

        size_t num_gtos = md_vlx_mo_gto_count(vlx);
        md_gto_t* gtos  = (md_gto_t*)md_temp_push(sizeof(md_gto_t) * num_gtos);

        num_gtos = md_vlx_mo_gto_extract(gtos, vlx, mo_idx, mo_type, cutoff_value);
        if (num_gtos == 0) {
            MD_LOG_ERROR("Failed to extract molecular gto for orbital index: %zu", mo_idx);
            return false;
        }

        // Dispatch Compute shader evaluation
        md_gto_grid_evaluate_GPU(vol_tex, &grid, gtos, num_gtos, mode);

        return true;
    }

    // The full electron density that includes all orbitals weighted by their occupancy
    bool compute_electron_density_GPU(uint32_t vol_tex, const md_grid_t& grid, md_vlx_mo_type_t mo_type, double cutoff_value = DEFAULT_GTO_CUTOFF_VALUE) {
        md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(2));
        defer { md_vm_arena_destroy(alloc); };

        #if 0
        md_orbital_data_t orb_data = {0};
        if (!extract_electron_density_orb_data(&orb_data, cutoff_value, alloc)) {
            return false;
        }

        // Dispatch Compute shader evaluation
        md_gto_grid_evaluate_orb_GPU(vol_tex, &grid, &orb_data, MD_GTO_EVAL_MODE_PSI_SQUARED);
        #else
        
        md_gto_data_t gto_data = {0};
        if (!md_vlx_scf_extract_gto_data(&gto_data, vlx, cutoff_value, alloc)) {
            MD_LOG_ERROR("Failed to extract gto data for electron density");
            return false;
        }
        
        //size_t matrix_size = md_vlx_scf_density_matrix_size(vlx);
        size_t matrix_size = md_vlx_scf_upper_triangular_density_matrix_size(vlx);
        size_t matrix_dim  = md_vlx_scf_density_matrix_size(vlx); 
        float* density_matrix_data = (float*)md_vm_arena_push(alloc, sizeof(float) * matrix_size);
        if (!md_vlx_scf_extract_upper_triangular_density_matrix_data(density_matrix_data, vlx, mo_type)) {
            MD_LOG_ERROR("Failed to extract density matrix data");
            return false;
        }
        
        // Use the density matrix and calculate density directly
        md_gto_grid_evaluate_matrix_GPU(vol_tex, &grid, &gto_data, density_matrix_data, matrix_dim, false);

        #endif

        return true;
    }

    task_system::ID compute_electron_density_async(uint32_t vol_tex, const md_grid_t& grid, md_vlx_mo_type_t mo_type, double cutoff_value = DEFAULT_GTO_CUTOFF_VALUE) {
        md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

        md_orbital_data_t orb_data = {0};
        if (!extract_electron_density_orb_data(&orb_data, cutoff_value, alloc)) {
            md_vm_arena_destroy(alloc);
            return task_system::INVALID_ID;
        }

        // Clear volume texture
        const float zero[4] = {0,0,0,0};
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, vol_fbo);
        glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, vol_tex, 0);
        glClearBufferfv(GL_COLOR, 0, zero);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

        struct Payload {
            AsyncGridEvalArgs args;
            uint32_t tex;
            md_allocator_i* arena;
        };

        Payload* payload = (Payload*)md_vm_arena_push(alloc, sizeof(Payload));
        *payload = {
            .args = {
                .grid = grid,
                .grid_data = (float*)md_vm_arena_push_zero(alloc, sizeof(float) * grid.dim[0] * grid.dim[1] * grid.dim[2]),
                .orb = orb_data,
                .mode = MD_GTO_EVAL_MODE_PSI_SQUARED,
             },
            .tex = vol_tex,
            .arena = alloc,
        };

        task_system::ID async_task = evaluate_gto_on_grid_async(&payload->args);

        // Launch task for main (render) thread to update the volume texture
        task_system::ID main_task = task_system::create_main_task(STR_LIT("##Update Volume"), [data = payload]() {
            // Ensure that the dimensions of the texture have not changed during evaluation
            int dim[3];
            if (gl::get_texture_dim(dim, data->tex) && MEMCMP(dim, data->args.grid.dim, sizeof(dim)) == 0) {
                gl::set_texture_3D_data(data->tex, data->args.grid_data, GL_R32F);
            }

            md_vm_arena_destroy(data->arena);
        });

        task_system::set_task_dependency(main_task, async_task);
        task_system::enqueue_task(async_task);

        return async_task;
    }

    task_system::ID compute_mo_async(uint32_t vol_tex, const md_grid_t& grid, md_vlx_mo_type_t mo_type, size_t mo_idx, md_gto_eval_mode_t mode, double cutoff_value = DEFAULT_GTO_CUTOFF_VALUE) {
        md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
        size_t num_gtos = md_vlx_mo_gto_count(vlx);
        md_gto_t* gtos  = (md_gto_t*)md_vm_arena_push(alloc, sizeof(md_gto_t)* num_gtos);

        num_gtos = md_vlx_mo_gto_extract(gtos, vlx, mo_idx, mo_type, cutoff_value);
        if (num_gtos == 0) {
            MD_LOG_ERROR("Failed to extract molecular gto for orbital index: %zu", mo_idx);
            md_vm_arena_destroy(alloc);
            return task_system::INVALID_ID;
        }

        md_orbital_data_t orb_data = {
            .num_gtos = num_gtos,
            .gtos = gtos,
        };

        struct Payload {
            AsyncGridEvalArgs args;
            uint32_t tex;
            md_allocator_i* arena;
        };

        Payload* payload = (Payload*)md_vm_arena_push(alloc, sizeof(Payload));
        *payload = {
            .args = {
                .grid = grid,
                .grid_data = (float*)md_vm_arena_push_zero(alloc, sizeof(float) * md_grid_num_points(&grid)),
                .orb = orb_data,
                .mode = mode,
            },
            .tex = vol_tex,
            .arena = alloc,
        };

        task_system::ID async_task = evaluate_gto_on_grid_async(&payload->args);

        // Launch task for main (render) thread to update the volume texture
        task_system::ID main_task = task_system::create_main_task(STR_LIT("##Update Volume"), [data = payload]() {
            // Ensure that the dimensions of the texture have not changed during evaluation
            int dim[3];
            if (gl::get_texture_dim(dim, data->tex) && MEMCMP(dim, data->args.grid.dim, sizeof(dim)) == 0) {
                gl::set_texture_3D_data(data->tex, data->args.grid_data, GL_R32F);
            }

            md_vm_arena_destroy(data->arena);
        });

        task_system::set_task_dependency(main_task, async_task);
        task_system::enqueue_task(async_task);

        return async_task;
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

    // Draws a vertical flow based on cubic bezier and returns the abs delta from the mouse cursor to this line
    static inline ImVec2 draw_vertical_sankey_flow(ImDrawList* draw_list, ImVec2 source_pos, ImVec2 dest_pos, float thickness, ImU32 flow_color) {
        // Define the control points for the Bezier curve
        ImVec2 p1 = ImVec2(source_pos.x + thickness / 2, source_pos.y);  // Start point (bottom-center of source node)
        ImVec2 p4 = ImVec2(dest_pos.x + thickness / 2, dest_pos.y);  // End point (top-center of destination node)

        float dist = fabsf(p1.y - p4.y);
        float curve_offset = dist / 4;

        ImVec2 p2 = ImVec2(source_pos.x + thickness / 2, source_pos.y - curve_offset);  // Control point 1
        ImVec2 p3 = ImVec2(dest_pos.x + thickness / 2, dest_pos.y + curve_offset);  // Control point 2

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
        ImVec2 text_size = ImGui::CalcTextSize(text);
        ImVec2 text_pos = pos - text_size * alignment;
        draw_list->AddText(text_pos, ImGui::ColorConvertFloat4ToU32({ 0,0,0,1 }), text);
    }

    static inline void im_sankey_diagram(ApplicationState* state, ImRect area, Nto* nto, bool hide_text_overlaps) {
        ScopedTemp temp_reset;
        md_allocator_i* temp_alloc = md_get_temp_allocator();
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

        md_array(float) hole_percentages = md_array_create(float, nto->group.count, temp_alloc);
        md_array(float) part_percentages = md_array_create(float, nto->group.count, temp_alloc);
        for (size_t i = 0; i < nto->group.count; i++) {
            hole_percentages[i] = nto->transition_density_hole[i] / hole_sum;
            part_percentages[i] = nto->transition_density_part[i] / part_sum;
        }

        //Calculate start positions
        md_array(float) start_positions = md_array_create(float, nto->group.count, temp_alloc);
        float cur_bottom_pos = plot_area.Min.x;
        for (int i = 0; i < nto->group.count; i++) {
            start_positions[i] = cur_bottom_pos;
            cur_bottom_pos += bars_avail_width * hole_percentages[i];
            if (hole_percentages[i] != 0.0) {
                cur_bottom_pos += gap_size;
            }
        }

        //Calculate end positions
        md_array(float) end_positions = md_array_create(float, nto->group.count, temp_alloc);
        md_array(float) sub_end_positions = md_array_create(float, nto->group.count, temp_alloc);
        float cur_pos = plot_area.Min.x;
        for (int end_i = 0; end_i < nto->group.count; end_i++) {
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
        md_array(float) curve_widths = md_array_create(float, nto->group.count * nto->group.count, temp_alloc);
        md_array(ImVec2) curve_midpoints = md_array_create(ImVec2, nto->group.count * nto->group.count, temp_alloc);
        md_array(ImVec2) curve_directions = md_array_create(ImVec2, nto->group.count * nto->group.count, temp_alloc);
        md_array(float) curve_percentages = md_array_create(float, nto->group.count * nto->group.count, temp_alloc);
        for (size_t i = 0; i < nto->group.count * nto->group.count; i++) {
            curve_widths[i] = 0;
            curve_midpoints[i] = { 0, 0 };
            curve_directions[i] = { 0, 0 };
            curve_percentages[i] = 0;
        }
        /*for (size_t start_i = 0; start_i < nto->group.count; start_i++) {
            for (size_t end_i = 0; end_i < nto->group.count; end_i++) {
                curve_widths[start_i * nto->group.count + end_i] = 0.0;
                curve_midpoints[start_i * nto->group.count + end_i] = ImVec2{};
                curve_directions[start_i * nto->group.count + end_i] = ImVec2{};
            }
        }*/

        for (int start_i = 0; start_i < nto->group.count; start_i++) {
            if (hole_percentages[start_i] != 0.0) {

                ImVec2 start_pos = { start_positions[start_i], plot_area.Max.y - bar_height + 0.1f * bar_height };
                ImVec4 start_col = vec_cast(nto->group.color[start_i]);

                start_col.w = 0.5;
                for (int end_i = 0; end_i < nto->group.count; end_i++) {
                    float percentage = nto->transition_matrix[end_i * nto->group.count + start_i];
                    ImVec4 end_col = vec_cast(nto->group.color[end_i]);

                    if (percentage != 0) {
                        float width = bars_avail_width * percentage;
                        ImVec2 end_pos = { sub_end_positions[end_i], plot_area.Min.y + bar_height - 0.1f * bar_height };

                        int vert_beg = draw_list->VtxBuffer.Size;
                        ImVec2 mouse_delta = draw_vertical_sankey_flow(draw_list, start_pos, end_pos, width, ImGui::ColorConvertFloat4ToU32(start_col));
                        int vert_end = draw_list->VtxBuffer.Size;

                        // The width test here is clamped to be atleast some pixel in width, otherwise it would be impossible to hit thin flows
                        bool flow_hovered = mouse_delta.x < MAX(width * 0.5f, min_test_width) && mouse_delta.y < 1.0e-4;

                        // Apply a gradient based on y-value from start color to end color if they belong to different groups
                        if (end_i != start_i) {
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
                            curve_midpoints[start_i * nto->group.count + end_i] = midpoint;
                            curve_directions[start_i * nto->group.count + end_i] = direction;
                            curve_widths[start_i * nto->group.count + end_i] = width;
                            curve_percentages[start_i * nto->group.count + end_i] = percentage;
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
            for (size_t start_i = 0; start_i < nto->group.count; start_i++) {
                for (size_t end_i = 0; end_i < nto->group.count; end_i++) {
                    size_t index1 = start_i * nto->group.count + end_i;
                    ImVec2 midpoint1 = curve_midpoints[index1];
                    ImVec2 direction1 = curve_directions[index1];
                    ImRect rect1 = { midpoint1 - (text_size / 2), midpoint1 + (text_size / 2) };

                    for (size_t start_j = 0; start_j < nto->group.count; start_j++) {
                        for (size_t end_j = 0; end_j < nto->group.count; end_j++) {
                            size_t index2 = (start_j * nto->group.count + end_j);
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

        for (size_t start_i = 0; start_i < nto->group.count; start_i++) {
            for (size_t end_i = 0; end_i < nto->group.count; end_i++) {
                size_t index = start_i * nto->group.count + end_i;
                ImVec2 midpoint = curve_midpoints[index];
                float percentage = curve_percentages[index];
                if (percentage > 0) {
                    char label[16];
                    sprintf(label, "%3.2f%%", curve_percentages[start_i * nto->group.count + end_i] * 100);
                    const ImVec2 label_size = ImGui::CalcTextSize(label);
                    draw_aligned_text(draw_list, label, midpoint, {0.5, 0.5});
                }
            }
        }

        // Some flow is hovered, So we render a new flow on top of the existing ones to ensure that this is rendered on top
        if (hovered_flow_width > 0.0f) {
            draw_vertical_sankey_flow(draw_list, hovered_flow_beg, hovered_flow_end, hovered_flow_width, ImColor(1.0f, 1.0f, 1.0f, 0.2f));
        }
        //bool next_start_text_visible = true;
        //bool next_end_text_visible = true;

        //Calculate text visibility
        md_array(bool) show_start_text  = md_array_create(bool,   nto->group.count, temp_alloc);
        md_array(bool) show_end_text    = md_array_create(bool,   nto->group.count, temp_alloc);
        md_array(int)  size_order       = md_array_create(int,    nto->group.count, temp_alloc);
        md_array(float) text_sizes      = md_array_create(float,  nto->group.count, temp_alloc);
        md_array(ImRect) rectangles     = md_array_create(ImRect, nto->group.count, temp_alloc);
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

            for (int i = 1; i < nto->group.count; i++) { //First one is always drawn
                for (int j = i - 1; j >= 0 ; j--) {//The bigger ones
                    if (hole_percentages[i] == 0.0) {
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
                    if (hole_percentages[i] == 0.0) {
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
        for (int i = 0; i < nto->group.count; i++) {
            if (hole_percentages[i] != 0.0) {
                ImVec4 bar_color = vec_cast(nto->group.color[i]);

                char start_label1[16];
                char end_label1[16]; 

                snprintf(start_label1, sizeof(start_label1), "%3.2f%%", hole_percentages[i] * 100);
                snprintf(end_label1, sizeof(end_label1), "%3.2f%%", part_percentages[i] * 100);

                char start_label2[16];
                char end_label2[16];
                snprintf(start_label2, sizeof(start_label1), "%3.2f%%", hole_percentages[i + 1] * 100);
                snprintf(end_label2, sizeof(end_label1), "%3.2f%%", part_percentages[i + 1] * 100);

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
                    for (size_t atom_i = 0; atom_i < nto->num_atoms; atom_i++) {
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

        if (mouse_label[0] != '\0') {
            ImVec2 offset = { 15, 15 };
            ImVec2 pos = ImGui::GetMousePos() + offset;
            draw_list->AddText(pos, IM_COL32_BLACK, mouse_label);
        }
    }
    
    bool compute_transition_group_values_async(task_system::ID* out_eval_task, task_system::ID* out_segment_task, float* out_group_values, size_t num_groups, const md_grid_t& grid, const uint32_t* point_group_idx, const vec4_t* point_xyzr, size_t num_points, size_t nto_idx, md_vlx_nto_type_t type, md_gto_eval_mode_t mode, float samples_per_angstrom = DEFAULT_SAMPLES_PER_ANGSTROM) {       
        md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

        md_orbital_data_t orb_data = {0};
        AttachmentDetachmentType attachment_detachment_type = (type == MD_VLX_NTO_TYPE_PARTICLE) ? AttachmentDetachmentType::Attachment : AttachmentDetachmentType::Detachment;
        if (!extract_attachment_detachment_orb_data(&orb_data, DEFAULT_GTO_CUTOFF_VALUE, attachment_detachment_type, nto_idx, alloc)) {
            MD_LOG_ERROR("Failed to extract attachment/detachment orbital data for NTO index: %zu", nto_idx);
            md_vm_arena_destroy(alloc);
            return false;
        }

        struct Payload {
            AsyncGridEvalArgs args;

            float* dst_group_values;
            size_t num_groups;

            const uint32_t* point_group_idx;
            const vec4_t* point_xyzr;
            size_t num_points;

            md_allocator_i* arena;
        };

        Payload* payload = (Payload*)md_vm_arena_push(alloc, sizeof(Payload));
        *payload = {
            .args = {
                .grid = grid,
                .grid_data = (float*)md_vm_arena_push_zero(alloc, sizeof(float) * md_grid_num_points(&grid)),
                .orb = orb_data,
                .mode = mode,
            },
            .dst_group_values = out_group_values,
            .num_groups = num_groups,
            
            .point_group_idx = point_group_idx,
            .point_xyzr = point_xyzr,
            .num_points = num_points,
            .arena = alloc,
        };

        task_system::ID eval_task = evaluate_gto_on_grid_async(&payload->args);

        // @TODO: This should be performed as a range task in parallel
        task_system::ID segment_task = task_system::create_pool_task(STR_LIT("##Segment Volume"), [data = payload]() {

            MD_LOG_DEBUG("Starting segmentation of volume");
            double* point_values = (double*)md_vm_arena_push_zero_array(data->arena, double, data->num_points);
            grid_segment_and_attribute_to_point(point_values, data->point_xyzr, data->num_points, data->args.grid_data, data->args.grid);

            // Accumulate values into groups
            for (size_t i = 0; i < data->num_points; ++i) {
                uint32_t group_idx = data->point_group_idx[i];
                if (group_idx >= data->num_groups) {
                    continue;
                }
                data->dst_group_values[group_idx] += (float)point_values[i];
            }

            // Print out point values for debugging
            for (size_t i = 0; i < data->num_points; ++i) {
                printf(" %zu: Value = %f\n", i, point_values[i]);
            }

            MD_LOG_DEBUG("Finished segmentation of volume");

            md_vm_arena_destroy(data->arena);
        });

        task_system::set_task_dependency(segment_task, eval_task);

        if (out_eval_task) {
            *out_eval_task = eval_task;
        }
        if (out_segment_task) {
            *out_segment_task = segment_task;
        }

        return true;
    }

    // This sets up and returns a async task which evaluates and orbital data on a grid in parallel
    task_system::ID evaluate_gto_on_grid_async(AsyncGridEvalArgs* args) {
        ASSERT(args);

        // We evaluate the in parallel over smaller NxNxN blocks
        const uint32_t num_blocks = (args->grid.dim[0] / BLK_DIM) * (args->grid.dim[1] / BLK_DIM) * (args->grid.dim[2] / BLK_DIM);
        task_system::ID async_task = task_system::create_pool_task(STR_LIT("Evaluate Orbital"), num_blocks, [data = args](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            (void)thread_num;
            MD_LOG_DEBUG("Starting async eval of orbital grid [%i][%i][%i]", data->grid.dim[0], data->grid.dim[1], data->grid.dim[2]);

            // Number of NxNxN blocks in each dimension
            int num_blk[3] = {
                data->grid.dim[0] / BLK_DIM,
                data->grid.dim[1] / BLK_DIM,
                data->grid.dim[2] / BLK_DIM,
            };

            md_grid_t* grid = &data->grid;
            mat4_t world_to_model = compute_world_to_model_mat(grid->orientation, grid->origin);

            size_t temp_pos = md_temp_get_pos();
            md_gto_t* sub_gtos = (md_gto_t*)md_temp_push(sizeof(md_gto_t) * data->orb.num_gtos);

            for (int blk_idx = (int)range_beg; blk_idx < (int)range_end; ++blk_idx) {
                // Determine block index from linear input index i
                int blk[3] = {
                    (blk_idx % num_blk[0]),
                    (blk_idx / num_blk[0]) % num_blk[1],
                    (blk_idx / (num_blk[0] * num_blk[1])),
                };

                int off_idx[3] = { blk[0] * BLK_DIM, blk[1] * BLK_DIM, blk[2] * BLK_DIM };
                int len_idx[3] = { BLK_DIM, BLK_DIM, BLK_DIM };

                int beg_idx[3] = {off_idx[0], off_idx[1], off_idx[2]};
                int end_idx[3] = {off_idx[0] + len_idx[0], off_idx[1] + len_idx[1], off_idx[2] + len_idx[2]};

                // The aabb is in model space
                vec4_t aabb_min = {
                    beg_idx[0] * data->grid.spacing.x,
                    beg_idx[1] * data->grid.spacing.y,
                    beg_idx[2] * data->grid.spacing.z,
                    0
                };
                vec4_t aabb_max = {
                    end_idx[0] * data->grid.spacing.x,
                    end_idx[1] * data->grid.spacing.y,
                    end_idx[2] * data->grid.spacing.z,
                    0
                };

                if (data->orb.num_orbs <= 1) {
                    size_t num_sub_gtos = 0;
                    // Transform GTO xyz into model space of the volume and perform AABB overlap test
                    for (size_t i = 0; i < data->orb.num_gtos; ++i) {
                        float cutoff = data->orb.gtos[i].cutoff;
                        if (cutoff == 0.0f) continue;

                        vec4_t coord = world_to_model * vec4_set(data->orb.gtos[i].x, data->orb.gtos[i].y, data->orb.gtos[i].z, 1.0f);
                        vec4_t clamped = vec4_clamp(coord, aabb_min, aabb_max);
                        if (vec4_distance_squared(coord, clamped) < cutoff * cutoff) {
                            sub_gtos[num_sub_gtos++] = data->orb.gtos[i];
                        }
                    }

                    md_gto_grid_evaluate_sub(data->grid_data, &data->grid, off_idx, len_idx, sub_gtos, num_sub_gtos, data->mode);
                } else {
                    for (size_t orb_idx = 0; orb_idx < data->orb.num_orbs; ++orb_idx) {
                        size_t num_sub_gtos = 0;
                        // Transform GTO xyz into model space of the volume and perform AABB overlap test
                        size_t beg = data->orb.orb_offsets[orb_idx];
                        size_t end = data->orb.orb_offsets[orb_idx + 1];
                        for (size_t i = beg; i < end; ++i) {
                            float cutoff = data->orb.gtos[i].cutoff;
                            if (cutoff == 0.0f) continue;

                            vec4_t coord = world_to_model * vec4_set(data->orb.gtos[i].x, data->orb.gtos[i].y, data->orb.gtos[i].z, 1.0f);
                            vec4_t clamped = vec4_clamp(coord, aabb_min, aabb_max);
                            if (vec4_distance_squared(coord, clamped) < cutoff * cutoff) {
                                sub_gtos[num_sub_gtos++] = data->orb.gtos[i];
                            }
                        }

                        md_gto_grid_evaluate_sub(data->grid_data, &data->grid, off_idx, len_idx, sub_gtos, num_sub_gtos, data->mode);
                    }
                }
            }

            md_temp_set_pos_back(temp_pos);
        });

        return async_task;
    }

    static inline double axis_conversion_multiplier(const double* y1_array, const double* y2_array, size_t y1_array_size, size_t y2_array_size) {
        double y1_max = 0;
        double y2_max = 0;
        for (size_t i = 0; i < y1_array_size; i++) {
            y1_max = MAX(y1_max, fabs(y1_array[i]));
        }
        for (size_t i = 0; i < y2_array_size; i++) {
            y2_max = MAX(y2_max, fabs(y2_array[i]));
        }

        return y2_max / y1_max;
    }

    static inline void convert_values(double* out_values, const double* in_values, size_t num_values, x_unit_t unit) {
        switch (unit) {
        case X_UNIT_EV:
            for (size_t i = 0; i < num_values; ++i) {
                out_values[i] = in_values[i];
            }
            break;
        case X_UNIT_NM:
            for (size_t i = 0; i < num_values; ++i) {
                out_values[i] = 1239.84193 / in_values[i];
            }
            break;
        case X_UNIT_CM_INVERSE:
            for (size_t i = 0; i < num_values; ++i) {
                out_values[i] = in_values[i] * 8065.73;
            }
            break;
        case X_UNIT_HARTREE:
            for (size_t i = 0; i < num_values; ++i) {
                out_values[i] = in_values[i] * 0.0367502;
            }
            break;
        default:
            ASSERT(false); // Should not happen
            break;
        }
    }

#if 0
    static inline double lorentzian(double x, double x_0, double gamma) {
        double sigma = gamma / 2;
        double res = (1 / PI) * sigma / (pow((x - x_0), 2) + pow(sigma, 2));
        return res;
    }

    static inline double gaussian(double x, double x_0, double gamma) {
        double sigma = gamma / 2.3548;
        return (1 / (sigma * sqrt(2 * PI))) * exp(-(pow(x - x_0, 2) / (2 * pow(sigma, 2)))); 
    }

#endif
    static inline double phys_lorentzian(double x, double x_0, double gamma, double intensity) {
        double sigma = gamma / 2;
        double res = intensity * pow(sigma, 2) / (pow((x - x_0), 2) + pow(sigma, 2));
        return res;
    }

    //TODO: Check that phys_gaussian implementation is actually correct
    static inline double phys_gaussian(double p, double p_0, double gamma, double intensity) {
        double sigma = gamma / 2;
        double x = (p - p_0) / sigma;
        return intensity * exp(-log(2) * pow(x,2));
    }

    static inline double broaden_normalized(double x, double x_0, double gamma, broadening_mode_t mode) {
        constexpr double inv_pi = 0.31830988618379067154;
        constexpr double inv_sqrt_2pi = 0.39894228040143267794;
        constexpr double inv_fwhm_to_sigma = 0.42466090014400953; // 1 / 2.3548

        const double dx = x - x_0;

        switch (mode) {
        case BROADENING_MODE_GAUSSIAN: {
            const double sigma = gamma * inv_fwhm_to_sigma;
            const double d = dx / sigma;
            return (inv_sqrt_2pi / sigma) * exp(-0.5 * d * d);
        }
        case BROADENING_MODE_LORENTZIAN: {
            const double sigma = 0.5 * gamma;
            const double denom = dx * dx + sigma * sigma;
            return (inv_pi * sigma) / denom;
        }
        default:
            ASSERT(false); // Should not happen
            return 0.0;
        }
    }

    static inline double broaden_physical(double x, double x0, double gamma, broadening_mode_t mode) {
        constexpr double inv_fwhm_to_sigma = 0.42466090014400953; // 1 / 2.3548
        const double dx = x - x0;
        switch (mode) {
        case BROADENING_MODE_GAUSSIAN: {
            const double sigma = gamma * inv_fwhm_to_sigma;
            const double d = dx / sigma;
            return exp(-0.5 * d * d);
        }
        case BROADENING_MODE_LORENTZIAN: {
            const double sigma = 0.5 * gamma;
            const double denom = dx * dx + sigma * sigma;
            return (sigma * sigma) / denom;
        }
        default:
            ASSERT(false); // Should not happen
            return 0.0;
		}
    }
    
    /*
    * We return to this at a later stage
    void save_absorption(str_t filename, md_array(double)* x_values, const char* x_lable, md_array(double)* y_values_osc, md_array(double)* y_values_cgs, int step) {
        md_file_o* file = md_file_open(filename, MD_FILE_WRITE);
        if (!file) {
            MD_LOG_ERROR("Could not open workspace file for writing: '%.*s", (int)filename.len, filename.ptr);
            return;
        }
        defer{ md_file_close(file); };

        // md_array(char*) rows = md_array_create(char*, md_array_size(*x_values), )
        //double* x_values = (double*)md_temp_push(sizeof(double) * num_samples);
        // I want to create an array with the rows to print to the file. What type should this be, and how do I create it?

        
        for (int i = 0; i < md_array_size(*x_values); i += 5) {

        }
    }


    */

    static inline void osc_to_eps(double* eps_out, const double* x, size_t num_samples, const double* osc_peaks, const double* x_peaks, size_t num_peaks, broadening_mode_t broadening_mode, double gamma) {
        // Physical constants
        const double c_au = 137.035999;             // speed of light (atomic units used in original expression)
        const double a0_m = 5.29177210903e-11;      // Bohr radius in meters
        const double NA = 6.02214076e23;            // Avogadro's number
        const double eV_to_au = 1.0 / 27.211396;    // electronvolt -> atomic units
        const double ln10 = 2.30258509299404568401799;

        // Precompute conversions used repeatedly
        const double gamma_au = gamma * eV_to_au;   // gamma in atomic units

        // conversion factor from cross section (m^2) to molar absorptivity ε (L mol^-1 cm^-1)
        const double to_molar_eps = (20.0 * PI * PI * a0_m * a0_m * NA) / (c_au * ln10);

        for (size_t si = 0; si < num_samples; ++si) {
            // Convert sample x to atomic units
            const double x_au = x[si] * eV_to_au;

            // Sum contributions from each oscillator
            double sum = 0.0;
            for (size_t pi = 0; pi < num_peaks; ++pi) {
                // Convert peak center to atomic units and compute the oscillator factor
                const double x0_au = x_peaks[pi] * eV_to_au;
                const double peak_factor = osc_peaks[pi] / (x0_au);
                sum += broaden_normalized(x_au, x0_au, gamma_au, broadening_mode) * peak_factor;
            }

            // Convert to molar absorptivity ε (L mol^-1 cm^-1)
            eps_out[si] = x_au * sum * to_molar_eps;
        }
    }

    static inline void rot_to_eps_delta(double* eps_out, const double* x, size_t num_samples, const double* rot_peaks, const double* x_peaks, size_t num_peaks, broadening_mode_t broadening_mode, double gamma) {
        static const double scl = 1.0 / (22.94);

        for (size_t si = 0; si < num_samples; si++) {
            double sum = 0;
            for (size_t pi = 0; pi < num_peaks; pi++) {
                sum += broaden_normalized(x[si], x_peaks[pi], gamma, broadening_mode) * rot_peaks[pi] * x_peaks[pi];
            }
            
            eps_out[si] = sum * scl;
        }
    }

    static inline void general_broadening(double* y_out, const double* x, size_t num_samples, const double* y_peaks, const double* x_peaks,
                                          size_t num_peaks, broadening_mode_t broadening_mode, double gamma) {
        double integral = 0;
        double dist = x[1] - x[0];
        for (size_t si = 0; si < num_samples; si++) {
            double sum = 0;
            for (size_t pi = 0; pi < num_peaks; pi++) {
                sum += y_peaks[pi] * broaden_physical(x[si], x_peaks[pi], gamma, broadening_mode);
            }
            y_out[si] = sum;
            integral += y_out[si] * dist;
        }

        double i_sum = integral;
    }

    // Compute a simple display multiplier to map oscillator-strength values (f) to eps for axis-syncing.
    // Strategy: use peak magnitudes: multiplier = max_abs(eps_samples)/max_abs(f_peaks).
    double compute_display_multiplier_from_peaks(const double* eps_samples, size_t eps_count,
        const double* peak_values, size_t peak_count) {
        double eps_max = 0.0;
        for (size_t i = 0; i < eps_count; ++i) eps_max = MAX(eps_max, fabs(eps_samples[i]));
        double peak_max = 0.0;
        for (size_t i = 0; i < peak_count; ++i) peak_max = MAX(peak_max, fabs(peak_values[i]));
        if (peak_max <= 0.0 || !isfinite(peak_max) || !isfinite(eps_max)) return 1.0;
        return eps_max / peak_max;
    }

    // This is a helper function to plot peaks and handle hovering
    // Only valid to call within a plot context
    void plot_peaks(const char* title, const double x_values[], const double y_values[], size_t num_peaks, int& selected, int& hovered) {
        ScopedTemp temp_reset;

        const double bar_width_in_pixels = 2.0;
        const double point_radius_in_pixels = 3.0;

        // If x values are not given we implicitly assume that they correspond to 1..num_peaks
        if (!x_values) {
            double* x_vals = (double*)md_temp_push(sizeof(double) * num_peaks);
            for (size_t i = 0; i < num_peaks; i++) {
                x_vals[i] = (double)(i + 1);
            }
            x_values = x_vals;
        }

        if (ImPlot::IsPlotHovered()) {
            // Calculate the hovered index as the closest peak or 'bar' to the mouse position
            // This comparison needs to be done in pixel space
            hovered = -1;

            const ImVec2 mouse_pos = ImGui::GetMousePos();
            const double pixel_tolerance = 7.0;
            const double min_line_d = (bar_width_in_pixels + 0.5 * pixel_tolerance);
			const double min_line_d2 = min_line_d * min_line_d;
            const double min_point_d = (point_radius_in_pixels + pixel_tolerance);
			const double min_point_d2 = min_point_d * min_point_d;

            double min_dist = DBL_MAX;
            for (size_t i = 0; i < num_peaks; i++) {
                ImVec2 min_raw = ImPlot::PlotToPixels(x_values[i], 0.0, IMPLOT_AUTO, IMPLOT_AUTO);
                ImVec2 max_raw = ImPlot::PlotToPixels(x_values[i], y_values[i], IMPLOT_AUTO, IMPLOT_AUTO);

                // The axis may be flipped, so we need to determine min and max correctly
                ImVec2 min = ImMin(min_raw, max_raw);
                ImVec2 max = ImMax(min_raw, max_raw);

                // Calculate distance to the line which forms the bar (1D clamp trick)
                ImVec2 line_dx = mouse_pos - ImClamp(mouse_pos, min, max);
                double line_d2 = line_dx.x * line_dx.x + line_dx.y * line_dx.y;

                ImVec2 point_dx = mouse_pos - max_raw;
                double point_d2 = point_dx.x * point_dx.x + point_dx.y * point_dx.y;

                if (line_d2  < min_line_d2  && line_d2  < min_dist ||
                    point_d2 < min_point_d2 && point_d2 < min_dist) {
                    hovered = (int)i;
                    min_dist = MIN(line_d2, point_d2);
                }
            }
        }

        // Update selected peak on click
        if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Left) &&
            ImPlot::IsPlotHovered()) {
            selected = hovered == selected ? -1 : hovered;
        }

        if (ImPlot::BeginItem(title)) {
            ImPlot::PushPlotClipRect();
            ImDrawList& DrawList = *ImPlot::GetPlotDrawList();
			const float rad = (float)point_radius_in_pixels;
			const float half_bar = (float)bar_width_in_pixels * 0.5f;
            for (int i = 0; i < (int)num_peaks; i++) {
                const ImPlotNextItemData& s = ImPlot::GetItemData();
                ImVec4 color = (i == hovered) ? COLOR_PEAK_HOVER : (i == selected) ? COLOR_PEAK_SELECTED : s.Colors[ImPlotCol_Fill];
                ImU32  col = ImGui::ColorConvertFloat4ToU32(color);
				ImVec2 btm = ImPlot::PlotToPixels(x_values[i], 0.0, IMPLOT_AUTO, IMPLOT_AUTO);
				ImVec2 top = ImPlot::PlotToPixels(x_values[i], y_values[i], IMPLOT_AUTO, IMPLOT_AUTO);

                if (half_bar > 0.0f) {
                    ImVec2 rect_min = ImVec2{ btm.x - half_bar, btm.y };
                    ImVec2 rect_max = ImVec2{ top.x + half_bar, top.y };
                    DrawList.AddRectFilled(rect_min, rect_max, col);
                }

                if (rad > 0.0f) {
                    DrawList.AddCircleFilled(top, rad, col);
                }

                if (ImPlot::FitThisFrame()) {
                    ImPlot::FitPoint(ImPlotPoint(x_values[i], 0.0));
                    ImPlot::FitPoint(ImPlotPoint(x_values[i], y_values[i]));
                }
            }
            ImPlot::PopPlotClipRect();
            ImPlot::EndItem();
        }
    }

    void set_atom_coordinates(ApplicationState& state, const dvec3_t* atom_coords, uint32_t flags = MolBit_DirtyPosition | MolBit_ClearVelocity) {
        const dvec3_t* coords = atom_coords ? atom_coords : md_vlx_atom_coordinates(vlx);
        size_t num_atoms = md_vlx_number_of_atoms(vlx);
        for (size_t i = 0; i < num_atoms; i++) {
            state.mold.sys.atom.x[i] = (float)coords[i].x;
            state.mold.sys.atom.y[i] = (float)coords[i].y;
            state.mold.sys.atom.z[i] = (float)coords[i].z;
        }
        state.mold.dirty_gpu_buffers |= flags;
    }

    void draw_summary_window(ApplicationState& state) {
        if (!summary.show_window) {
            opt.selected = (int)md_vlx_opt_number_of_steps(vlx) - 1;
            return;
        }

        static ImPlotRect lims{ 0,1,0,1 };

        size_t temp_pos = md_temp_get_pos();
        defer {  md_temp_set_pos_back(temp_pos); };

        // The actual plot
        ImGui::SetNextWindowSize({ 300, 350 }, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Summary", &summary.show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::TreeNode("Level of Calculation")) {
                str_t basis_set = md_vlx_basis_set_ident(vlx);
                str_t dft_func  = md_vlx_dft_func_label(vlx);
                ImGui::Text("Method: %s", str_ptr(dft_func));
                ImGui::Text("Basis Set: %s", str_ptr(basis_set));
                ImGui::Spacing();
                
                ImGui::TreePop();
            }
            if (ImGui::TreeNode("System Information")) {
                ImGui::Text("Num Atoms:           %-6zu", md_vlx_number_of_atoms(vlx));
                ImGui::Text("Num Alpha Electrons: %-6zu", md_vlx_number_of_alpha_electrons(vlx));
                ImGui::Text("Num Beta Electrons:  %-6zu", md_vlx_number_of_beta_electrons(vlx));
                ImGui::Text("Molecular Charge:    %-6f",  md_vlx_molecular_charge(vlx));
                ImGui::Text("Spin Multiplicity:   %-6zu", md_vlx_spin_multiplicity(vlx));
                ImGui::Spacing();
                ImGui::TreePop();
            }

            if (ImGui::TreeNode("SCF")) {
                size_t num_iter = md_vlx_scf_history_size(vlx);
                const double* grad_norm = md_vlx_scf_history_gradient_norm(vlx);
                const double* energy = md_vlx_scf_history_energy(vlx);
                double* energy_offsets = nullptr;
                double ref_energy = 0.0;
                double y1_to_y2_mult = 1.0;

                // We set up iterations as doubles for easier use
                if (num_iter > 0) {
                    ASSERT(energy);
                    ASSERT(grad_norm);

                    energy_offsets = (double*)md_temp_push(sizeof(double) * num_iter);
                    ref_energy = energy[num_iter - 1];
                    for (size_t i = 0; i < num_iter; i++) {
                        energy_offsets[i] = fabs(energy[i] - ref_energy);
                    }
                    y1_to_y2_mult = axis_conversion_multiplier(grad_norm, energy_offsets, num_iter, num_iter);
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
                        const double* energies = md_vlx_opt_energies(vlx);
                        double ref_energy = energies[num_steps - 1];
                        double* energy_offsets = (double*)md_temp_push(sizeof(double) * num_steps);

                        for (size_t i = 0; i < num_steps; ++i) {
                            energy_offsets[i] = fabs(energies[i] - ref_energy) * HARTREE_TO_KJ_PER_MOL;
                        }

                        if (ImPlot::BeginPlot("OPT")) {
                            ImPlot::SetupAxes("Step", "Energy (kJ/mol)");
                            ImPlot::SetupAxisLimits(ImAxis_X1, 1.0, (double)num_steps);
                            ImPlot::SetupLegend(ImPlotLocation_NorthEast);
                            ImPlot::PlotLine("Energy", energy_offsets, (int)num_steps, 1.0, 1.0);

                            plot_peaks("##opt_peaks", nullptr, energy_offsets, num_steps, opt.selected, opt.hovered);

                            ImPlot::EndPlot();
                        }
                        ImGui::Spacing();
                        int value = opt.selected + 1;
                        ImGui::SliderInt("Step", &value, 1, (int)num_steps);
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
                    const dvec3_t* atom_coord = md_vlx_atom_coordinates(vlx);
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
                    for (int row_n = 0; row_n < num_atoms; row_n++) {

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

                        char lable[16];
                        sprintf(lable, "%i", row_n + 1);
                        ImGui::Selectable(lable, is_sel || is_hov, selectable_flags);
                        if (ImGui::TableGetHoveredRow() == row_n + 1) {
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
                ImGui::Combo("Volume Resolution", (int*)&vol_res, "Low\0Mid\0High\0");

                if (critical_points.enabled) {
                    static bool  reevaluate_graph_per_frame = false;
                    static bool  enable_simplification = true;
                    static float simp_threshold = 1.0E-4f;

                    ImGui::Checkbox("Extract graph per frame", &reevaluate_graph_per_frame);
                    if (reevaluate_graph_per_frame) {
                        critical_points.raw_hash = 0;
                    }

                    // Compute a hash of the current state
                    uint64_t vol_hash = md_hash64(&vol_res, sizeof(vol_res), 0);

                    // Update volume if required
                    if (critical_points.vol_hash != vol_hash || critical_points.density_vol.tex_id == 0) {
                        critical_points.vol_hash = vol_hash;
                        init_grid(&critical_points.grid, obb.orientation, obb.min_ext, obb.max_ext, volume_resolution_samples_per_angstrom[(int)vol_res]);
                        init_volume(&critical_points.density_vol, critical_points.grid, GL_R32F);
                        if (!compute_electron_density_GPU(critical_points.density_vol.tex_id, critical_points.grid, MD_VLX_MO_TYPE_ALPHA)) {
                            MD_LOG_ERROR("Failed to compute electron density volume for critical points analysis");
                        }
                    }

                    uint64_t raw_hash = vol_hash;

                    if (critical_points.raw_hash != raw_hash) {
                        critical_points.raw_hash = raw_hash;
                        md_topo_extremum_graph_free(&critical_points.raw_graph);
                        if (!md_topo_compute_extremum_graph_GPU(&critical_points.raw_graph, critical_points.density_vol.tex_id, &critical_points.grid, 0.0f)) {
                            MD_LOG_ERROR("Failed to compute extremum graph for electron density");
                        }
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
                        if (enable_simplification) {
                            //const md_topo_extremum_graph_t& graph = critical_points.raw_graph;
                            md_allocator_i* temp_alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
                            defer{ md_arena_allocator_destroy(temp_alloc); };

                            // Extract explicit vertices and edges to work with
                            md_array(int) vertex_type       = md_array_create(int, critical_points.raw_graph.num_vertices, temp_alloc);
                            md_array(md_topo_vert_t) vertex = md_array_create(md_topo_vert_t, critical_points.raw_graph.num_vertices, temp_alloc);
                            md_array(md_topo_edge_t) edge   = md_array_create(md_topo_edge_t, critical_points.raw_graph.num_edges, temp_alloc);

                            // Construct adjacency information for each vertex
                            md_array(md_array(int)) vertex_adj = md_array_create(md_array(int), critical_points.raw_graph.num_vertices, temp_alloc);

                            md_topo_extract_vertex_types(vertex_type, md_array_size(vertex_type), &critical_points.raw_graph);
                            MEMCPY(vertex, critical_points.raw_graph.vertices, md_array_bytes(vertex));
                            MEMSET(vertex_adj, 0, md_array_bytes(vertex_adj));
                            MEMCPY(edge, critical_points.raw_graph.edges, md_array_bytes(edge));

                            for (size_t i = 0; i < md_array_size(edge); ++i) {
                                md_topo_edge_t e = edge[i];
                                md_array_push(vertex_adj[e.from], (int)e.to, temp_alloc);
                                md_array_push(vertex_adj[e.to], (int)e.from, temp_alloc);
                            }

                            // Placeholders
                            bool threshold_simplification = true;
                            bool prune_duplicate_saddles  = true;

                            if (threshold_simplification) {
                                for (int i = 0; i < md_array_size(vertex); ++i) {
                                    if (vertex[i].value < simp_threshold) {
                                        vertex_type[i] = 0;  // Mark vertex as dead
                                    }
                                }
                            }

                            if (prune_duplicate_saddles) {
                                // What we want to do here is to find all pairs of maxima that are connected by multiple saddles
                                // From this set of saddles, we want to keep only the one with the highest value, and mark the rest for removal
                                // If there is a saddle that is connected to those two maximas, but also connected to others, we want to reindex it to
                                // the 'maximum' saddle within the set

                                // Divide multi connected saddles (3+ connections) into single connections
                                for (int i = 0; i < md_array_size(vertex); ++i) {
                                    if (vertex_type[i] != MD_TOPO_SPLIT_SADDLE) continue;
                                    size_t num_adj = md_array_size(vertex_adj[i]);
                                    if (num_adj <= 2) continue;  // Only interested in saddles with 3 or more connections
                                    // For each pair of adjacent maxima, create a new saddle vertex
                                    for (size_t j = 0; j < num_adj - 1; ++j) {
                                        for (size_t k = j + 1; k < num_adj; ++k) {
                                            int max_a = vertex_adj[i][j];
                                            int max_b = vertex_adj[i][k];
                                            // Create new saddle vertex
                                            md_topo_vert_t new_saddle = vertex[i];
                                            md_array_push(vertex, new_saddle, temp_alloc);
                                            md_array_push(vertex_adj, NULL,   temp_alloc);
                                            int new_saddle_idx = (int)(md_array_size(vertex) - 1);
                                            vertex_type[new_saddle_idx] = MD_TOPO_SPLIT_SADDLE;
                                            // Update adjacency
                                            md_array_push(vertex_adj[new_saddle_idx], max_a, temp_alloc);
                                            md_array_push(vertex_adj[new_saddle_idx], max_b, temp_alloc);
                                            // Update edges
                                            md_topo_edge_t edge_a = { .from = (uint32_t)max_a, .to = (uint32_t)new_saddle_idx };
                                            md_topo_edge_t edge_b = { .from = (uint32_t)max_b, .to = (uint32_t)new_saddle_idx };
                                            md_array_push(edge, edge_a, temp_alloc);
                                            md_array_push(edge, edge_b, temp_alloc);
                                        }
                                    }
                                    // Mark original saddle as dead
                                    vertex_type[i] = 0;
                                }

                                md_array(int) saddle_list = 0;

                                // Iterate over all Maxima pairs
                                for (size_t i = 0; i < md_array_size(vertex) - 1; ++i) {
                                    if (vertex_type[i] != MD_TOPO_MAXIMUM) continue;
                                    for (size_t j = i + 1; j < md_array_size(vertex); ++j) {
                                        if (vertex_type[j] != MD_TOPO_MAXIMUM) continue;

                                        md_array_shrink(saddle_list, 0);

                                        // Find all saddles connecting the Maxima i and j
                                        for (size_t k = 0; k < md_array_size(vertex); ++k) {
                                            if (vertex_type[k] != MD_TOPO_SPLIT_SADDLE) continue;
                                            // Test all saddles to see if they connect the pair of maximas
                                            size_t num_adj = md_array_size(vertex_adj[k]);
                                            if (num_adj != 2) continue;  // Saddles must connect exactly two maxima

                                            bool connects_i = false;
                                            bool connects_j = false;
                                            for (size_t m = 0; m < num_adj; ++m) {
                                                int v = vertex_adj[k][m];
                                                if (v == i) connects_i = true;
                                                if (v == j) connects_j = true;
                                            }
                                            if (connects_i && connects_j) {
                                                md_array_push(saddle_list, (int)k, temp_alloc);
                                            }
                                        }

                                        // If we found multiple saddles, we need to prune them
                                        if (md_array_size(saddle_list) > 1) {

                                            float saddle_max = -FLT_MAX;
                                            int   saddle_max_idx = -1;

                                            // Find maximum saddle
                                            for (size_t k = 0; k < md_array_size(saddle_list); ++k) {
                                                int idx = saddle_list[k];
                                                if (vertex[idx].value > saddle_max) {
                                                    saddle_max = vertex[idx].value;
                                                    saddle_max_idx = idx;
                                                }
                                            }

                                            // Mark all other saddles as dead
                                            for (size_t k = 0; k < md_array_size(saddle_list); ++k) {
                                                int idx = saddle_list[k];
                                                if (idx == saddle_max_idx) continue;
                                                vertex_type[idx] = 0;  // Mark saddle as dead
                                            }
                                        }
                                    }
                                }
                            }

                            // Construct a new simplified graph
                            uint32_t num_vertices = 0;
                            uint32_t num_maxima = 0;
                            uint32_t num_split_saddles = 0;
                            uint32_t num_minima = 0;
                            uint32_t num_join_saddles = 0;

                            // Create a remapping from old to new vertex indices
                            md_array(int) vertex_remap = md_array_create(int, md_array_size(vertex), temp_alloc);
                            MEMSET(vertex_remap, -1, md_array_bytes(vertex_remap));

                            for (size_t i = 0; i < md_array_size(vertex_type); ++i) {
                                int new_idx = 0;
                                switch (vertex_type[i]) {
                                    case MD_TOPO_MAXIMUM:      new_idx = num_maxima++;        break;
                                    case MD_TOPO_SPLIT_SADDLE: new_idx = num_split_saddles++; break;
                                    case MD_TOPO_MINIMUM:      new_idx = num_minima++;        break;
                                    case MD_TOPO_JOIN_SADDLE:  new_idx = num_join_saddles++;  break;
                                    default: continue;
                                }
                                vertex_remap[i] = new_idx;
                                num_vertices++;
                            }

                            uint32_t offset_maxima = 0;
                            uint32_t offset_split_saddle = offset_maxima + num_maxima;
                            uint32_t offset_minima = offset_split_saddle + num_split_saddles;
                            uint32_t offset_join_saddle = offset_minima + num_minima;

                            // Adjust remapping to global indices
                            for (size_t i = 0; i < md_array_size(vertex_type); ++i) {
                                switch (vertex_type[i]) {
                                    case MD_TOPO_SPLIT_SADDLE: vertex_remap[i] += offset_split_saddle;  break;
                                    case MD_TOPO_MINIMUM:      vertex_remap[i] += offset_minima;        break;
                                    case MD_TOPO_JOIN_SADDLE:  vertex_remap[i] += offset_join_saddle;   break;
                                }
                            }

                            md_topo_extremum_graph_t& simp_graph = critical_points.simp_graph;
                            simp_graph.num_vertices         = num_vertices;
                            simp_graph.num_maxima           = num_maxima;
                            simp_graph.num_split_saddles    = num_split_saddles;
                            simp_graph.num_minima           = num_minima;
                            simp_graph.num_join_saddles     = num_join_saddles;

                            // Allocate vertex data
                            simp_graph.vertices = (md_topo_vert_t*)md_alloc(simp_graph.alloc, sizeof(md_topo_vert_t) * simp_graph.num_vertices);

                            // Fill in vertex data
                            for (size_t i = 0; i < md_array_size(vertex_type); ++i) {
                                int idx = vertex_remap[i];
                                if (idx != -1) {
                                    simp_graph.vertices[idx] = vertex[i];
                                }
                            }

                            size_t num_edges = 0;
                            for (size_t i = 0; i < md_array_size(edge); ++i) {
                                md_topo_edge_t e = edge[i];
                                if (vertex_remap[e.from] != -1 && vertex_remap[e.to] != -1) {
                                    num_edges++;
                                }
                            }

                            simp_graph.num_edges = num_edges;
                            simp_graph.edges     = (md_topo_edge_t*)md_alloc(simp_graph.alloc, sizeof(md_topo_edge_t) * simp_graph.num_edges);

                            // Fill in edge data
                            size_t edge_counter = 0;
                            for (size_t i = 0; i < md_array_size(edge); ++i) {
                                md_topo_edge_t old_edge = edge[i];
                                md_topo_edge_t new_edge;
                                new_edge.from = vertex_remap[old_edge.from];
                                new_edge.to   = vertex_remap[old_edge.to];
                                if (new_edge.from != -1 && new_edge.to != -1) {
                                    simp_graph.edges[edge_counter++] = new_edge;
                                }
                            }
                        }
                        else {
                            md_topo_extremum_graph_copy(&critical_points.simp_graph, &critical_points.raw_graph);
                        }
                    }

                    const md_topo_extremum_graph_t& graph = critical_points.simp_graph;
                    ImGui::Text("Number of Edges:           %zu", graph.num_edges);
                    ImGui::Text("Number of Critical Points: %zu", graph.num_vertices);
                    ImGui::Text("\tNumber of Maxima:        %zu", graph.num_maxima);
                    ImGui::Text("\tNumber of Split Saddles: %zu", graph.num_split_saddles);
                    ImGui::Text("\tNumber of Minima:        %zu", graph.num_minima);
                    ImGui::Text("\tNumber of Join Saddles:  %zu", graph.num_join_saddles);

                    if (ImGui::Button("Print Critical Points Info")) {
                        size_t temp_pos = md_temp_get_pos();
                        defer{ md_temp_set_pos_back(temp_pos); };
                        md_array(int) vertex_types = md_array_create(int, critical_points.simp_graph.num_vertices, md_get_temp_allocator());
                        md_topo_extract_vertex_types(vertex_types, md_array_size(vertex_types), &critical_points.simp_graph);
                        for (size_t i = 0; i < critical_points.simp_graph.num_vertices; ++i) {
                            const md_topo_vert_t& v = critical_points.simp_graph.vertices[i];
							const char* type_cstr = "Unknown";
                            switch (vertex_types[i]) {
                            case MD_TOPO_MAXIMUM:      type_cstr = "Maximum";       break;
                            case MD_TOPO_SPLIT_SADDLE: type_cstr = "Split Saddle";  break;
                            case MD_TOPO_MINIMUM:      type_cstr = "Minimum";       break;
                            case MD_TOPO_JOIN_SADDLE:  type_cstr = "Join Saddle";   break;
                            }
                            printf("Vertex %zu: Type=%s, Value=%.6f, Pos=(%.4f, %.4f, %.4f)\n", i, type_cstr, v.value, v.x, v.y, v.z);
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
            {rsp.x_unit_samples, rsp.eps,   STR_LIT("Absorption"), STR_LIT((const char*)u8"ε (L mol⁻¹ cm⁻¹)")},
            {rsp.x_unit_samples, rsp.ecd,   STR_LIT("ECD"),        STR_LIT((const char*)u8"Δε(ω) (L mol⁻¹ cm⁻¹)")},
#if 0
            {rsp.vib_x,          rsp.vib_y, STR_LIT("Vibration"),  STR_LIT("IR Intensity (km/mol)")}
#endif
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
                for (int i = 0; i < ARRAY_SIZE(properties); ++i) {
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
                export_valid = rsp.x_unit_samples && rsp.eps && rsp.ecd;

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
                            md_file_o* file = md_file_open(path, MD_FILE_WRITE | MD_FILE_BINARY);
                            if (file) {
                                const size_t written_bytes = md_file_write(file, xvg.ptr, xvg.len);
                                if (written_bytes != xvg.len) {
                                    MD_LOG_ERROR("CSV: Unexpected error, some bytes were not written");
                                }
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
    }

    void draw_rsp_linear() {
		int num_excited_states = (int)md_vlx_rsp_number_of_excited_states(vlx);
        if (ImGui::TreeNode("Absorption Analysis")) {
            if (rsp.first_plot_rot_ecd) {
                ImGui::SetNextItemOpen(true);
            }
            bool refit  = false;
            bool recalc = false;

            rsp.hovered = -1;

            const float avail_width = ImGui::GetContentRegionAvail().x;
            ImGui::PushItemWidth(MIN(avail_width, 200));
            recalc |= ImGui::SliderScalar((const char*)u8"Broadening γ HWHM (eV)", ImGuiDataType_Double, &rsp.broadening_gamma, &gamma_min, &gamma_max);
            recalc |= ImGui::Combo("Broadening mode", (int*)(&rsp.broadening_mode), broadening_mode_str, BROADENING_MODE_COUNT);
            refit  |= ImGui::Combo("X unit", (int*)(&rsp.x_unit), x_unit_full_str, X_UNIT_COUNT);
            ImGui::PopItemWidth();

            refit  |= recalc;

            const int num_peaks = (int)num_excited_states;

            const double* y_osc_peaks = md_vlx_rsp_oscillator_strengths(vlx);
            const double* y_cgs_peaks = md_vlx_rsp_rotatory_strengths(vlx);
            const double* x_abs_ev = md_vlx_rsp_absorption_ev(vlx);

            ImVec2* pixel_osc_peaks  = (ImVec2*)md_temp_push(sizeof(ImVec2) * num_peaks);
            ImVec2* pixel_cgs_peaks  = (ImVec2*)md_temp_push(sizeof(ImVec2) * num_peaks);
            ImVec2* pixel_osc_points = (ImVec2*)md_temp_push(sizeof(ImVec2) * num_peaks);
            ImVec2* pixel_cgs_points = (ImVec2*)md_temp_push(sizeof(ImVec2) * num_peaks);

            if (recalc || rsp.first_plot_rot_ecd) {
                osc_to_eps(rsp.eps, rsp.x_ev_samples, NUM_SAMPLES, y_osc_peaks, x_abs_ev, num_peaks, rsp.broadening_mode, rsp.broadening_gamma * 2);
                rot_to_eps_delta(rsp.ecd, rsp.x_ev_samples, NUM_SAMPLES, y_cgs_peaks, x_abs_ev, num_peaks, rsp.broadening_mode, rsp.broadening_gamma * 2);
                rsp.y_axis_sync_factor_eps = compute_display_multiplier_from_peaks(rsp.eps, NUM_SAMPLES, y_osc_peaks, num_peaks);
                rsp.y_axis_sync_factor_ecd = compute_display_multiplier_from_peaks(rsp.ecd, NUM_SAMPLES, y_cgs_peaks, num_peaks);
            }

            if (refit || rsp.first_plot_rot_ecd) {
                // Do unit conversions (X axis)
                convert_values(rsp.x_unit_peaks, x_abs_ev, num_peaks, rsp.x_unit);
                convert_values(rsp.x_unit_samples, rsp.x_ev_samples, NUM_SAMPLES, rsp.x_unit);
            }

            // Explicit x axis limits for linking
            static double x_min = 0.0;
            static double x_max = 10.0;

            ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
            const ImGuiTreeNodeFlags tree_flags = ImGuiTreeNodeFlags_DefaultOpen | ImGuiTreeNodeFlags_NoTreePushOnOpen;
            if (ImGui::TreeNodeEx("Absorption Plot", tree_flags)) {
                ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

                if (refit || rsp.first_plot_rot_ecd) {
                    ImPlot::SetNextAxesToFit();
                }

                if (ImPlot::BeginPlot("Absorption")) {
                    ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                    ImPlot::SetupAxis(ImAxis_X1, x_unit_full_str[rsp.x_unit]);
                    ImPlot::SetupAxis(ImAxis_Y1, "f", ImPlotAxisFlags_AuxDefault);
                    ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"ε (L mol⁻¹ cm⁻¹)");
                    ImPlot::SetupFinish();

                    ImPlot::SetAxis(ImAxis_Y2);
                    ImPlot::PlotLine("Spectrum", rsp.x_unit_samples, rsp.eps, NUM_SAMPLES);

                    ImPlot::SetAxis(ImAxis_Y1);
                    plot_peaks("Oscillator Strength", rsp.x_unit_peaks, y_osc_peaks, num_peaks, rsp.selected, rsp.hovered);

                    if (ImPlot::IsPlotHovered() && rsp.hovered != -1) {
                        double x_val = rsp.x_unit_peaks[rsp.hovered];
                        double y_val = y_osc_peaks[rsp.hovered];
                        int state = rsp.hovered + 1;
                        const char* x_label = x_unit_label_str[rsp.x_unit];
                        const char* x_unit = x_unit_short_str[rsp.x_unit];
                        if (ImGui::BeginTooltip()) {
                            ImGui::Text("State %i: %s = %.2f (%s), f = %.3f", state, x_label, x_val, x_unit, y_val);
                            ImGui::EndTooltip();
                        }
                    }

                    // adjust plot axes before endplot, to keep y-axes in sync
                    ImPlot::SyncAxesWithPadding(ImAxis_Y2, ImAxis_Y1, rsp.y_axis_sync_factor_eps);
                    ImPlot::EndPlot();
                }
            }

            ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
            if (ImGui::TreeNodeEx("ECD Plot", tree_flags)) {
                ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

                if (refit || rsp.first_plot_rot_ecd) {
                    ImPlot::SetNextAxesToFit();
                }

                if (ImPlot::BeginPlot("ECD")) {
                    ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                    ImPlot::SetupAxis(ImAxis_X1, x_unit_full_str[rsp.x_unit]);
                    ImPlot::SetupAxis(ImAxis_Y1, (const char*)u8"R (10⁻⁴⁰ cgs)", ImPlotAxisFlags_AuxDefault);
                    ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"Δε(ω) (L mol⁻¹ cm⁻¹)");
                    ImPlot::SetupFinish();

                    ImPlot::SetAxis(ImAxis_Y2);
                    ImPlot::PlotLine("Spectrum", rsp.x_unit_samples, rsp.ecd, NUM_SAMPLES);

                    ImPlot::SetAxis(ImAxis_Y1);
                    plot_peaks("Rotatory Strength", rsp.x_unit_peaks, y_cgs_peaks, num_peaks, rsp.selected, rsp.hovered);

                    if (ImPlot::IsPlotHovered() && rsp.hovered != -1) {
                        double x_val = rsp.x_unit_peaks[rsp.hovered];
                        double y_val = y_cgs_peaks[rsp.hovered];
                        int state = rsp.hovered + 1;
                        const char* x_label = x_unit_label_str[rsp.x_unit];
                        const char* x_unit  = x_unit_short_str[rsp.x_unit];
                        if (ImGui::BeginTooltip()) {
                            ImGui::Text((const char*)u8"State %i: %s = %.2f (%s), R = %.3f (10⁻⁴⁰ cgs)", state, x_label, x_val, x_unit, y_val);
                            ImGui::EndTooltip();
                        }
                    }


                    // adjust plot axes before endplot, to keep zeros aligned
                    ImPlot::SyncAxesY();
                    ImPlot::EndPlot();
                }
            }
            rsp.first_plot_rot_ecd = false;

            // Table
            static const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY |
                ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_Sortable;

            static const ImGuiTableColumnFlags columns_base_flags = ImGuiTableColumnFlags_DefaultSort;

            ImVec2 table_size = { 0, 0 };

            const double* x_values = rsp.x_unit_peaks;
            const double* f_values = y_osc_peaks;
            const double* r_values = y_cgs_peaks;

            if (ImGui::BeginTable("abs table", 4, flags, table_size, 0)) {
                ImGui::TableSetupColumn("State", columns_base_flags, 0.0f);
                ImGui::TableSetupColumn(x_unit_full_str[rsp.x_unit], columns_base_flags, 0.0f);
                ImGui::TableSetupColumn("f", columns_base_flags, 0.0f);
                ImGui::TableSetupColumn((const char*)u8"R (10⁻⁴⁰ cgs)", columns_base_flags, 0.0f);
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableHeadersRow();

                // TODO: Add sorting to the table once the vibs data structure is properly defined

                ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);
                ImGui::PushStyleColor(ImGuiCol_Header, IM_BLUE);
                bool item_hovered = false;
                for (int row_n = 0; row_n < (int)num_excited_states; row_n++) {

                    ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                    bool is_sel = row_n == rsp.selected;
                    bool is_hov = row_n == rsp.hovered;
                    bool hov_col = false;
                    ImGui::TableNextRow(ImGuiTableRowFlags_None, 0);
                    ImGui::TableNextColumn();

                    if (is_sel) {
                        ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_BLUE);
                    }
                    else {
                        ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);
                    }

                    char lable[16];
                    sprintf(lable, "%i", row_n + 1);
                    if (ImGui::Selectable(lable, is_sel || is_hov, selectable_flags)) {
                        rsp.selected = (rsp.selected == row_n) ? -1 : row_n;
                    }
                    ImGui::TableNextColumn();
                    ImGui::Text("%12.6f", x_values[row_n]);
                    ImGui::TableNextColumn();
                    ImGui::Text("%12.6f", f_values[row_n]);
                    ImGui::TableNextColumn();
                    ImGui::Text("%12.6f", r_values[row_n]);

                    ImGui::PopStyleColor(1);
                }
                if (ImGui::IsWindowHovered() && ImGui::TableGetHoveredRow() > 0) {
                    rsp.hovered = ImGui::TableGetHoveredRow() - 1;
                }
                else {
                    rsp.hovered = -1;
                }

                ImGui::PopStyleColor(2);
                ImGui::EndTable();
            }
            ImGui::TreePop();
        }
    }

    void draw_rsp_cpp() {
		static double x_min = 0.0;
		static double x_max = 10.0;

        /*
        ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
        const ImGuiTreeNodeFlags tree_flags = ImGuiTreeNodeFlags_DefaultOpen | ImGuiTreeNodeFlags_NoTreePushOnOpen;
        if (ImGui::TreeNodeEx("Absorption Plot", tree_flags)) {
            ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

            if (rsp.first_plot_rot_ecd) {
                ImPlot::SetNextAxesToFit();
            }

            if (ImPlot::BeginPlot("Absorption")) {
                ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                ImPlot::SetupAxis(ImAxis_X1, x_unit_full_str[rsp.x_unit]);
                ImPlot::SetupAxis(ImAxis_Y1, "f", ImPlotAxisFlags_AuxDefault);
                ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"ε (L mol⁻¹ cm⁻¹)");
                ImPlot::SetupFinish();

                ImPlot::SetAxis(ImAxis_Y2);
                ImPlot::PlotLine("Spectrum", rsp.x_unit_samples, rsp.eps, NUM_SAMPLES);

                ImPlot::SetAxis(ImAxis_Y1);
                plot_peaks("Oscillator Strength", rsp.x_unit_peaks, y_osc_peaks, num_peaks, rsp.selected, rsp.hovered);

                if (ImPlot::IsPlotHovered() && rsp.hovered != -1) {
                    double x_val = rsp.x_unit_peaks[rsp.hovered];
                    double y_val = y_osc_peaks[rsp.hovered];
                    int state = rsp.hovered + 1;
                    const char* x_label = x_unit_label_str[rsp.x_unit];
                    const char* x_unit  = x_unit_short_str[rsp.x_unit];
                    if (ImGui::BeginTooltip()) {
                        ImGui::Text("State %i: %s = %.2f (%s), f = %.3f", state, x_label, x_val, x_unit, y_val);
                        ImGui::EndTooltip();
                    }
                }

                // adjust plot axes before endplot, to keep y-axes in sync
                ImPlot::SyncAxesWithPadding(ImAxis_Y2, ImAxis_Y1, rsp.y_axis_sync_factor_eps);
                ImPlot::EndPlot();
            }
        }

        ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
        if (ImGui::TreeNodeEx("ECD Plot", tree_flags)) {
            ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

            if (refit || rsp.first_plot_rot_ecd) {
                ImPlot::SetNextAxesToFit();
            }

            if (ImPlot::BeginPlot("ECD")) {
                ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                ImPlot::SetupAxis(ImAxis_X1, x_unit_full_str[rsp.x_unit]);
                ImPlot::SetupAxis(ImAxis_Y1, (const char*)u8"R (10⁻⁴⁰ cgs)", ImPlotAxisFlags_AuxDefault);
                ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"Δε(ω) (L mol⁻¹ cm⁻¹)");
                ImPlot::SetupFinish();

                ImPlot::SetAxis(ImAxis_Y2);
                ImPlot::PlotLine("Spectrum", rsp.x_unit_samples, rsp.ecd, NUM_SAMPLES);

                ImPlot::SetAxis(ImAxis_Y1);
                plot_peaks("Rotatory Strength", rsp.x_unit_peaks, y_cgs_peaks, num_peaks, rsp.selected, rsp.hovered);

                if (ImPlot::IsPlotHovered() && rsp.hovered != -1) {
                    double x_val = rsp.x_unit_peaks[rsp.hovered];
                    double y_val = y_cgs_peaks[rsp.hovered];
                    int state = rsp.hovered + 1;
                    const char* x_label = x_unit_label_str[rsp.x_unit];
                    const char* x_unit  = x_unit_short_str[rsp.x_unit];
                    if (ImGui::BeginTooltip()) {
                        ImGui::Text((const char*)u8"State %i: %s = %.2f (%s), R = %.3f (10⁻⁴⁰ cgs)", state, x_label, x_val, x_unit, y_val);
                        ImGui::EndTooltip();
                    }
                }


                // adjust plot axes before endplot, to keep zeros aligned
                ImPlot::SyncAxesY();
                ImPlot::EndPlot();
            }
        }
        rsp.first_plot_rot_ecd = false;
        */
    }

    void draw_rsp_window(ApplicationState& state) {
        ScopedTemp temp_reset;

        if (!rsp.show_window) return;
        size_t num_excited_states = md_vlx_rsp_number_of_excited_states(vlx);
        size_t num_normal_modes   = md_vlx_vib_number_of_normal_modes(vlx);
		size_t num_frequencies    = md_vlx_rsp_cpp_number_of_frequencies(vlx);
        if (num_excited_states == 0 && num_normal_modes == 0 && num_frequencies == 0) return;

        ImGui::SetNextWindowSize({ 300, 350 }, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Response", &rsp.show_window, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("File")) {
                    if (ImGui::MenuItem("Export")) {
                        rsp.show_export_window = true;
                    }
                    ImGui::EndMenu();
                }
                ImGui::EndMenuBar();
            }
            
            if (num_excited_states > 0) {
                draw_rsp_linear();
			}

            if (num_frequencies > 0) {
                draw_rsp_cpp();
            }

            if (num_normal_modes > 0) {
                if (vib.first_plot) {
                    ImGui::SetNextItemOpen(true);
                }
                if (ImGui::TreeNode("Vibrational Analysis")) {
                    size_t num_atoms = md_vlx_number_of_atoms(vlx);

                    // Frequency scaling factor limits
                    static const double scale_min = 0.85;
                    static const double scale_max = 1.05;

                    bool recalc = false;

                    const float avail_width = ImGui::GetContentRegionAvail().x;
                    ImGui::PushItemWidth(MIN(avail_width, 200));

                    recalc |= ImGui::SliderFloat((const char*)u8"Broadening γ HWHM (cm⁻¹)", &vib.gamma, 1.0f, 100.0f, "%.3f", ImGuiSliderFlags_Logarithmic);
                    recalc |= ImGui::Combo("Broadening mode", (int*)(&vib.broadening_mode), broadening_mode_str, IM_ARRAYSIZE(broadening_mode_str));
                    recalc |= ImGui::SliderScalar("Scaling factor", ImGuiDataType_Double, &vib.freq_scaling_factor, &scale_min, &scale_max, "%.4f");
                    ImGui::PopItemWidth();

                    ImGui::SetItemTooltip("Frequency scaling factor");

                    const double* y_values = md_vlx_vib_ir_intensities(vlx);
                    const double* x_values_raw = md_vlx_vib_frequencies(vlx);
                    double* x_values = (double*)md_temp_push(sizeof(double) * num_normal_modes);
                    for (size_t i = 0; i < num_normal_modes; ++i) {
                        x_values[i] = x_values_raw[i] * vib.freq_scaling_factor;
                    }

                    ImVec2* pixel_peaks = (ImVec2*)md_temp_push(sizeof(ImVec2) * num_normal_modes);

                    if (vib.first_plot) {
                        // Populate x_values
                        const double x_min = x_values[0] - 100.0;
                        const double x_max = x_values[num_normal_modes - 1] + 100.0;
                        for (int i = 0; i < NUM_SAMPLES; ++i) {
                            double t = (double)i / (double)(ARRAY_SIZE(vib.x_samples) - 1);
                            double value = lerp(x_min, x_max, t);
                            vib.x_samples[i] = value;
                        }
                    }

                    if (vib.first_plot || recalc) {
                        general_broadening(vib.y_samples, vib.x_samples, NUM_SAMPLES, y_values, x_values, num_normal_modes, vib.broadening_mode, vib.gamma * 2);
                    }

                    ImGui::Checkbox("Invert X", &vib.invert_x);
                    ImGui::SameLine();
                    ImGui::Checkbox("Invert Y", &vib.invert_y);

                    ImPlotAxisFlags x_flag = vib.invert_x ? ImPlotAxisFlags_Invert : 0;
                    ImPlotAxisFlags y_flag = vib.invert_y ? ImPlotAxisFlags_Invert : 0;

                    if (ImPlot::BeginPlot("Vibrational analysis")) {
                        ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                        ImPlot::SetupAxis(ImAxis_X1, (const char*)u8"Harmonic Frequency (cm⁻¹)", x_flag);
                        ImPlot::SetupAxis(ImAxis_Y1, "IR Intensity (km/mol)", y_flag);
                        ImPlot::SetupFinish();

                        ImPlot::PlotLine("Spectrum", vib.x_samples, vib.y_samples, NUM_SAMPLES);
                        plot_peaks("IR Intensity", x_values, y_values, num_normal_modes, vib.selected, vib.hovered);

                        vib.first_plot = false;
                        ImPlot::EndPlot();
                    }

                    ImGui::PushItemWidth(MIN(avail_width, 200));

                    int idx = vib.selected + 1;
                    ImGui::SliderInt("Index", &idx, 0, (int)num_normal_modes);
                    vib.selected = CLAMP(idx - 1, -1, (int)num_normal_modes - 1);
                    ImGui::SliderFloat((const char*)"Amplitude", &vib.displacement_amp_scl, 0.25f, 2.0f);
                    ImGui::SetItemTooltip("Displacement Amplitude Scale");
                    ImGui::SliderFloat((const char*)"Frequency", &vib.displacement_freq_scl, 0.25f, 2.0f);
                    ImGui::SetItemTooltip("Displacement Frequency Scale");
                    ImGui::PopItemWidth();

                    // Table
                    static const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY |
                                                         ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_Sortable;

                    static const ImGuiTableColumnFlags columns_base_flags = ImGuiTableColumnFlags_DefaultSort;

                    ImVec2 table_size = {0, 0};

                    if (ImGui::BeginTable("vib mode table", 3, flags, table_size, 0)) {
                        ImGui::TableSetupColumn("Vibration mode", columns_base_flags, 0.0f);
                        ImGui::TableSetupColumn("Harmonic Frequency", columns_base_flags, 0.0f);
                        ImGui::TableSetupColumn("IR Intensity", columns_base_flags, 0.0f);
                        ImGui::TableSetupScrollFreeze(0, 1);
                        ImGui::TableHeadersRow();

                        // TODO: Add sorting to the table once the vibs data structure is properly defined

                        ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);
                        ImGui::PushStyleColor(ImGuiCol_Header, IM_BLUE);
                        bool item_hovered = false;
                        for (int row_n = 0; row_n < (int)num_normal_modes; row_n++) {

                            ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                            bool is_sel = row_n == vib.selected;
                            bool is_hov = row_n == vib.hovered;
                            bool hov_col = false;
                            ImGui::TableNextRow(ImGuiTableRowFlags_None, 0);
                            ImGui::TableNextColumn();

                            if (is_sel) {
                                ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_BLUE);
                            } else {
                                ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);
                            }

                            char lable[16];
                            sprintf(lable, "%i", row_n + 1);
                            if (ImGui::Selectable(lable, is_sel || is_hov, selectable_flags)) {
                                vib.selected = (vib.selected == row_n) ? -1 : row_n;
                            }
                            ImGui::TableNextColumn();
                            ImGui::Text("%12.6f", x_values[row_n]);
                            ImGui::TableNextColumn();
                            ImGui::Text("%12.6f", y_values[row_n]);

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
                        // Animation
                        vib.t += state.app.timing.delta_s * vib.displacement_freq_scl * 8.0;
                        const double scl = vib.displacement_amp_scl * 0.25 * sin(vib.t);
                        const dvec3_t* norm_modes = md_vlx_vib_normal_mode(vlx, vib.selected);

                        if (norm_modes) {
                            for (size_t i = 0; i < num_atoms; i++) {
                                state.mold.sys.atom.x[i] = (float)(atom_coord[i].x + norm_modes[i].x * scl);
                                state.mold.sys.atom.y[i] = (float)(atom_coord[i].y + norm_modes[i].y * scl);
                                state.mold.sys.atom.z[i] = (float)(atom_coord[i].z + norm_modes[i].z * scl);
                            }
                            state.mold.dirty_gpu_buffers |= MolBit_DirtyPosition;
                            vib.coord_modified = true;
                        }
                    }
                    // If all is deselected, reset coords once
                    else if (vib.coord_modified) {
                        for (size_t i = 0; i < num_atoms; i++) {
                            state.mold.sys.atom.x[i] = (float)atom_coord[i].x;
                            state.mold.sys.atom.y[i] = (float)atom_coord[i].y;
                            state.mold.sys.atom.z[i] = (float)atom_coord[i].z;
                        }
                        state.mold.dirty_gpu_buffers |= MolBit_DirtyPosition | MolBit_ClearVelocity;
                        vib.coord_modified = false;
                    }

                    ImGui::TreePop();
                }
            }
        }
        ImGui::End();

        if (rsp.show_export_window) { draw_rsp_spectra_export_window(state); }
    }

    void draw_orb_window(const ApplicationState& state) {
        if (!orb.show_window) return;
        if (num_molecular_orbitals() == 0) return;
        ImGui::SetNextWindowSize({600, 300}, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Orbital Grid", &orb.show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {

            const double* occ_alpha = md_vlx_scf_mo_occupancy(vlx, MD_VLX_MO_TYPE_ALPHA);
            const double* occ_beta  = md_vlx_scf_mo_occupancy(vlx, MD_VLX_MO_TYPE_BETA);
            const double* ene_alpha = md_vlx_scf_mo_energy(vlx, MD_VLX_MO_TYPE_ALPHA);
            const double* ene_beta  = md_vlx_scf_mo_energy(vlx, MD_VLX_MO_TYPE_BETA);

            const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();
            md_vlx_scf_type_t type = md_vlx_scf_type(vlx);

            int num_x = (type == MD_VLX_SCF_TYPE_UNRESTRICTED) ? 2 : orb.num_x;
            int num_y = orb.num_y;

            int orb_homo_idx = MAX(homo_idx[0], homo_idx[1]);
            int orb_lumo_idx = MIN(lumo_idx[0], lumo_idx[1]);
            if (type == MD_VLX_SCF_TYPE_RESTRICTED_OPENSHELL) {
                orb_lumo_idx = MAX(lumo_idx[0], lumo_idx[1]);
            }

            int num_mos = num_x * num_y;
            int beg_mo_idx = orb.mo_idx - num_mos / 2 + (num_mos % 2 == 0 ? 1 : 0);
            int window_size = num_mos;
            if (type == MD_VLX_SCF_TYPE_UNRESTRICTED) {
                beg_mo_idx = orb.mo_idx - num_y / 2 + (num_y % 2 == 0 ? 1 : 0);
                window_size = num_y;
            }

            // LEFT PANE
            {
                ImGui::BeginChild("left pane", ImVec2(300, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
                ImGui::SameLine();

                const ImVec2 outer_size = {300.f, 0.f};
                ImGui::PushItemWidth(-1);
                ImGui::BeginGroup();

                ImGui::SliderInt("##Rows", &orb.num_y, 1, 4);
                if (type == MD_VLX_SCF_TYPE_UNRESTRICTED) ImGui::PushDisabled();
                ImGui::SliderInt("##Cols", &orb.num_x, 1, 4);
                if (type == MD_VLX_SCF_TYPE_UNRESTRICTED) ImGui::PopDisabled();

                const double iso_min = 1.0e-4;
                const double iso_max = 5.0;
                double iso_val = orb.iso.values[0];
                ImGui::SliderScalar("##Iso Value", ImGuiDataType_Double, &iso_val, &iso_min, &iso_max, "%.6f", ImGuiSliderFlags_Logarithmic);
                ImGui::SetItemTooltip("Iso Value");

                orb.iso.values[0] =  (float)iso_val;
                orb.iso.values[1] = -(float)iso_val;
                orb.iso.count = 2;
                orb.iso.enabled = true;

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
                if (type == MD_VLX_SCF_TYPE_RESTRICTED_OPENSHELL && (homo_idx[0] != homo_idx[1])) {
                    btn_text = "Goto SOMO";
                }

                if (ImGui::IsWindowAppearing()) {
                    orb.scroll_to_idx = orb.mo_idx;
                }
                if (ImGui::Button(btn_text, ImVec2(-1, 0))) {
                    orb.scroll_to_idx = homo_idx[0];
                }

                int num_cols = (type == MD_VLX_SCF_TYPE_UNRESTRICTED) ? 5 : 3;

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
                    if (type == MD_VLX_SCF_TYPE_UNRESTRICTED) {
                        ImGui::TableSetupColumn("MO",                       ImGuiTableColumnFlags_DefaultSort | ImGuiTableColumnFlags_WidthFixed,         0.0f, Col_Idx);
                        ImGui::TableSetupColumn((const char*)u8"Occ. α",  ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed,  0.0f, Col_Occ_Alpha);
                        ImGui::TableSetupColumn((const char*)u8"Occ. β",  ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed,  0.0f, Col_Occ_Beta);
                        ImGui::TableSetupColumn((const char*)u8"Ene. α",  ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed,  0.0f, Col_Ene_Alpha);
                        ImGui::TableSetupColumn((const char*)u8"Ene. β",  ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed,  0.0f, Col_Ene_Beta);
                    } else {
                        ImGui::TableSetupColumn("MO",           ImGuiTableColumnFlags_DefaultSort          | ImGuiTableColumnFlags_WidthFixed,   0.0f, Col_Idx);
                        ImGui::TableSetupColumn("Occupancy",    ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed,   0.0f, Col_Occ_Alpha);
                        ImGui::TableSetupColumn("Energy",       ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed,   0.0f, Col_Ene_Alpha);
                    }
                    ImGui::TableSetupScrollFreeze(0, 1); // Make row always visible
                    ImGui::TableHeadersRow();

                    for (int n = (int)num_molecular_orbitals() - 1; n >= 0; n--) {
                        ImGui::PushID(n + 1);
                        ImGui::TableNextRow();
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
                            if (type == MD_VLX_SCF_TYPE_UNRESTRICTED) {
                                if (n == homo_idx[0] && n == orb_homo_idx) {
                                    lbl = "HOMO";
                                } else if (n == lumo_idx[0] && n == orb_lumo_idx) {
                                    lbl = "LUMO";
                                }
                            } else {
                                if (occ_beta) {
                                    occ += occ_beta[n];
                                }
                                if (n == orb_homo_idx) {
                                    lbl = (homo_idx[0] == homo_idx[1]) ? "HOMO" : "SOMO";
                                } else if (n == orb_lumo_idx) {
                                    lbl = "LUMO";
                                }
                            }
                            ImGui::Text("%.1f %s", occ, lbl);
                        } else {
                            ImGui::Text("-");
                        }
                        if (type == MD_VLX_SCF_TYPE_UNRESTRICTED) {
                            ImGui::TableNextColumn();
                            if (occ_beta) {
                                const char* lbl = "";
                                if (n == homo_idx[1] && n == orb_homo_idx) {
                                    lbl = "HOMO";
                                } else if (n == lumo_idx[1] && n == orb_lumo_idx) {
                                    lbl = "LUMO";
                                }
                                ImGui::Text("%.1f %s", occ_beta[n], lbl);
                            } else {
                                ImGui::Text("-");
                            }
                        }
                        ImGui::TableNextColumn();
                        if (ene_alpha) {
                            ImGui::Text("%.4f", ene_alpha[n]);
                        } else {
                            ImGui::Text("-");
                        }
                        if (type == MD_VLX_SCF_TYPE_UNRESTRICTED) {
                            ImGui::TableNextColumn();
                            if (ene_beta) {
                                ImGui::Text("%.4f", ene_beta[n]);
                            } else {
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
            int vol_mo_idx[16]  = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
            md_vlx_mo_type_t vol_mo_type[16] = {};
            for (int i = 0; i < num_mos; ++i) {
                if (type == MD_VLX_SCF_TYPE_UNRESTRICTED) {
                    vol_mo_idx[i] = beg_mo_idx + i / 2;
                    vol_mo_type[i] = (i & 1) ? MD_VLX_MO_TYPE_ALPHA : MD_VLX_MO_TYPE_BETA;
                } else {
                    vol_mo_idx[i] = beg_mo_idx + i;
                    vol_mo_type[i] = MD_VLX_MO_TYPE_ALPHA;
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
                        ImSwap(orb.vol_mo_idx[i],  orb.vol_mo_idx[j]);
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
                md_grid_t grid {0};
                init_grid(&grid, obb.orientation, obb.min_ext, obb.max_ext, samples_per_unit_length);

                for (int i = 0; i < num_jobs; ++i) {
                    int slot_idx = job_queue[i];
                    int mo_idx = vol_mo_idx[slot_idx];
                    md_vlx_mo_type_t mo_type = vol_mo_type[slot_idx];
                    orb.vol_mo_idx[slot_idx] = mo_idx;
                    orb.vol_mo_type[slot_idx] = mo_type;

                    if (-1 < mo_idx && mo_idx < num_molecular_orbitals()) {
                        if (task_system::task_is_running(orb.vol_task[slot_idx])) {
                            task_system::task_interrupt(orb.vol_task[slot_idx]);
                        }

                        init_volume(&orb.vol[slot_idx], grid);

                        if (use_gpu_path) {
                            compute_mo_GPU(orb.vol[slot_idx].tex_id, grid, mo_type, mo_idx, MD_GTO_EVAL_MODE_PSI);
                        } else {
                            orb.vol_task[slot_idx] = compute_mo_async(orb.vol[slot_idx].tex_id, grid, mo_type, mo_idx, MD_GTO_EVAL_MODE_PSI);
                        }
                    }
                }
            }

            // Animate camera towards targets
            const double dt = state.app.timing.delta_s;
            camera_animate(&orb.camera, orb.target.ori, orb.target.pos, orb.target.dist, dt);

            ImVec2 canvas_sz = ImGui::GetContentRegionAvail();   // Resize canvas to what's available
            canvas_sz.x = MAX(canvas_sz.x, 50.0f);
            canvas_sz.y = MAX(canvas_sz.y, 50.0f);

            // This will catch our interactions
            ImGui::InvisibleButton("canvas", canvas_sz, ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_AllowOverlap);

            // Draw border and background color
            ImGuiIO& io = ImGui::GetIO();

            ImVec2 canvas_p0 = ImGui::GetItemRectMin();
            ImVec2 canvas_p1 = ImGui::GetItemRectMax();

            ImVec2 orb_win_sz = (canvas_p1 - canvas_p0) / ImVec2((float)num_x, (float)num_y);
            orb_win_sz.x = floorf(orb_win_sz.x);
            orb_win_sz.y = floorf(orb_win_sz.y);
            canvas_p1.x = canvas_p0.x + num_x * orb_win_sz.x;
            canvas_p1.y = canvas_p0.y + num_y * orb_win_sz.y;

            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(255, 255, 255, 255));

            const int num_orbs = (int)num_molecular_orbitals();

            for (int i = 0; i < num_mos; ++i) {
                int mo_idx = beg_mo_idx + i;
                if (type == MD_VLX_SCF_TYPE_UNRESTRICTED) {
                    mo_idx = beg_mo_idx + i / 2;
                }
                int x = num_x - i % num_x - 1;
                int y = num_y - i / num_x - 1;
                ImVec2 p0 = canvas_p0 + orb_win_sz * ImVec2((float)(x+0), (float)(y+0));
                ImVec2 p1 = canvas_p0 + orb_win_sz * ImVec2((float)(x+1), (float)(y+1));
                if (-1 < mo_idx && mo_idx < num_orbs) {
                    const ImVec2 text_pos_bl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p1.y - TEXT_BASE_HEIGHT);
                    const ImVec2 text_pos_tl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p0.y + TEXT_BASE_HEIGHT * 0.25f);
                    const ImVec2 text_pos_br = ImVec2(p1.x - TEXT_BASE_HEIGHT * 0.5f, p1.y - TEXT_BASE_HEIGHT);

                    const char* lbl = "";
                    if (type == MD_VLX_SCF_TYPE_UNRESTRICTED) {
                        int j = (i & 1) ? 0 : 1;
                        if (mo_idx == orb_homo_idx && orb_homo_idx == homo_idx[j]) {
                            lbl = "(HOMO)";
                        } else if (mo_idx == orb_lumo_idx && orb_lumo_idx == lumo_idx[j]) {
                            lbl = "(LUMO)";
                        }
                    } else {
                        if (mo_idx == orb_homo_idx) {
                            lbl = (homo_idx[0] == homo_idx[1]) ? "(HOMO)" : "(SOMO)";
                        } else if (mo_idx == orb_lumo_idx) {
                            lbl = "(LUMO)";
                        }
                    }

                    char buf[32];
                    snprintf(buf, sizeof(buf), "%i %s", mo_idx + 1, lbl);
                    draw_list->AddImage((ImTextureID)(intptr_t)orb.gbuf.tex.transparency, p0, p1, { 0,1 }, { 1,0 });
                    draw_list->AddImage((ImTextureID)(intptr_t)orb.iso_tex[i], p0, p1, { 0,1 }, { 1,0 });
                    draw_list->AddText(text_pos_bl, ImColor(0,0,0), buf);

                    if (type == MD_VLX_SCF_TYPE_UNRESTRICTED) {
                        draw_list->AddText(text_pos_tl, ImColor(0, 0, 0), (i & 1) ? (const char*)u8"α" : (const char*)u8"β");
                        snprintf(buf, sizeof(buf), "%.4f", (i & 1) ? ene_alpha[mo_idx] : ene_beta[mo_idx]);
                    } else {
                        snprintf(buf, sizeof(buf), "%.4f", ene_alpha[mo_idx]);
                    }
                    float width = ImGui::CalcTextSize(buf).x;
                    draw_list->AddText(text_pos_br - ImVec2(width, 0), ImColor(0,0,0), buf);
                }
            }
            for (int x = 1; x < orb.num_x; ++x) {
                ImVec2 p0 = {canvas_p0.x + orb_win_sz.x * x, canvas_p0.y};
                ImVec2 p1 = {canvas_p0.x + orb_win_sz.x * x, canvas_p1.y};
                draw_list->AddLine(p0, p1, IM_COL32(0, 0, 0, 255));
            }
            for (int y = 1; y < orb.num_y; ++y) {
                ImVec2 p0 = {canvas_p0.x, canvas_p0.y + orb_win_sz.y * y};
                ImVec2 p1 = {canvas_p1.x, canvas_p0.y + orb_win_sz.y * y};
                draw_list->AddLine(p0, p1, IM_COL32(0, 0, 0, 255));
            }

            const bool is_hovered = ImGui::IsItemHovered();
            const bool is_active = ImGui::IsItemActive();
            const ImVec2 origin(canvas_p0.x, canvas_p0.y);  // Lock scrolled origin
            const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

            int width  = MAX(1, (int)orb_win_sz.x);
            int height = MAX(1, (int)orb_win_sz.y);

            auto& gbuf = orb.gbuf;
            if ((int)gbuf.width != width || (int)gbuf.height != height) {
                init_gbuffer(&gbuf, width, height);
                for (int i = 0; i < num_mos; ++i) {
                    gl::init_texture_2D(orb.iso_tex + i, width, height, GL_RGBA8);
                }
            }

            bool reset_hard = false;
            bool reset_view = false;
            if (is_hovered) {
                if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                    reset_view = true;
                }
            }

            if (reset_view) {
                camera_compute_optimal_view(&orb.target.pos, &orb.target.ori, &orb.target.dist, obb.orientation, obb.min_ext * BOHR_TO_ANGSTROM, obb.max_ext * BOHR_TO_ANGSTROM, orb.distance_scale);

                if (reset_hard) {
                    orb.camera.position         = orb.target.pos;
                    orb.camera.orientation      = orb.target.ori;
                    orb.camera.focus_distance   = orb.target.dist;
                }
            }

            if (is_active || is_hovered) {
                const vec2_t delta = { io.MouseDelta.x, io.MouseDelta.y };
                const vec2_t curr = {mouse_pos_in_canvas.x, mouse_pos_in_canvas.y};
                const vec2_t prev = curr - delta;

                TrackballControllerInput input = {
                    .rotate_button = is_active && ImGui::IsMouseDown(ImGuiMouseButton_Left),
                    .pan_button    = is_active && ImGui::IsMouseDown(ImGuiMouseButton_Right),
                    .dolly_button  = is_active && ImGui::IsMouseDown(ImGuiMouseButton_Middle),
                    .dolly_delta   = is_hovered ? io.MouseWheel : 0.0f,
                    .mouse_coord_prev = prev,
                    .mouse_coord_curr = curr,
                    .screen_size = {canvas_sz.x, canvas_sz.y},
                    .fov_y = orb.camera.fov_y,
                };
                camera_controller_trackball(&orb.target.pos, &orb.target.ori, &orb.target.dist, input);
            }

            if (orb.show_coordinate_system_widget) {
                float  ext = MIN(orb_win_sz.x, orb_win_sz.y) * 0.4f;
                float  pad = 20.0f;

                ImVec2 min = ImGui::GetItemRectMin() - ImGui::GetWindowPos();
                ImVec2 max = ImGui::GetItemRectMax() - ImGui::GetWindowPos();

                CoordSystemWidgetParam param = {
                    .pos = ImVec2(min.x + pad, max.y - ext - pad),
                    .size = {ext, ext},
                    .view_matrix = camera_world_to_view_matrix(orb.camera),
                    .camera_ori  = orb.target.ori,
                    .camera_pos  = orb.target.pos,
                    .camera_dist = orb.target.dist,
                };

                ImGui::DrawCoordinateSystemWidget(param);
            }

            if (gl_rep.id) {
                const float aspect_ratio = orb_win_sz.x / orb_win_sz.y;
                mat4_t view_mat = camera_world_to_view_matrix(orb.camera);
                mat4_t proj_mat = camera_perspective_projection_matrix(orb.camera, aspect_ratio);
                mat4_t inv_proj_mat = camera_inverse_perspective_projection_matrix(orb.camera, aspect_ratio);

                clear_gbuffer(&gbuf);

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
                    .shaders = state.mold.gl_shaders,
                    .draw_operations = {
                        .count = 1,
                        .ops = &draw_op
                    },
                    .view_transform = {
                        .view_matrix = (const float*)view_mat.elem,
                        .proj_matrix = (const float*)proj_mat.elem,
                    },
                };

                md_gl_draw(&draw_args);

                glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
                glClearColor(1, 1, 1, 0);
                glClear(GL_COLOR_BUFFER_BIT);

                PUSH_GPU_SECTION("Postprocessing")
                postprocessing::Descriptor postprocess_desc = {
                    .background = {
                        .color = {24.f, 24.f, 24.f},
                    },
                    .tonemapping = {
                        .enabled    = state.visuals.tonemapping.enabled,
                        .mode       = state.visuals.tonemapping.tonemapper,
                        .exposure   = state.visuals.tonemapping.exposure,
                        .gamma      = state.visuals.tonemapping.gamma,
                    },
                    .ambient_occlusion = {
                        .enabled = false
                    },
                    .depth_of_field = {
                        .enabled = false,
                    },
                    .fxaa = {
                        .enabled = true,
                    },
                    .temporal_aa = {
                        .enabled = false,
                    },
                    .sharpen = {
                        .enabled = false,
                    },
                    .input_textures = {
                        .depth          = orb.gbuf.tex.depth,
                        .color          = orb.gbuf.tex.color,
                        .normal         = orb.gbuf.tex.normal,
                        .velocity       = orb.gbuf.tex.velocity,
                    }
                };

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

                postprocessing::shade_and_postprocess(postprocess_desc, view_param);
                POP_GPU_SECTION()

                glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
                glDrawBuffer(GL_BACK);
                glDisable(GL_SCISSOR_TEST);

                if (orb.iso.enabled) {
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
                                    .volume = orb.vol[i].tex_id,
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
        }
        ImGui::End();
    }

    void draw_export_window(ApplicationState& state) {
        if (!export_state.show_window) return;
        if (ImGui::Begin("VeloxChem Export", &export_state.show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::BeginCombo("Orbital Type", electronic_structure_type_str[(int)export_state.orb_type])) {
                for (int i = 0; i < (int)ElectronicStructureType::Count; i++) {
                    bool is_selected = ((int)export_state.orb_type == i);
                    bool disabled = false;

                    if (i == (int)ElectronicStructureType::AttachmentDensity || i == (int)ElectronicStructureType::DetachmentDensity) {
                        disabled = (md_vlx_rsp_number_of_excited_states(vlx) == 0);    
                    }

                    if (disabled) ImGui::PushDisabled();
                    if (ImGui::Selectable(electronic_structure_type_str[i], is_selected)) {
                        export_state.orb_type = (ElectronicStructureType)i;
                    }
                    if (disabled) ImGui::PopDisabled();

                    if (is_selected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }

            int orb_idx = 0;
            bool show_mo_combo = false;
            bool show_lambda_combo = false;
            bool show_excited_state_combo = false;

            switch (export_state.orb_type) {
            case ElectronicStructureType::MolecularOrbital:
            case ElectronicStructureType::MolecularOrbitalDensity:
                show_mo_combo = true;
                break;
            case ElectronicStructureType::NaturalTransitionOrbitalParticle:
            case ElectronicStructureType::NaturalTransitionOrbitalHole:
            case ElectronicStructureType::NaturalTransitionOrbitalDensityParticle:
            case ElectronicStructureType::NaturalTransitionOrbitalDensityHole:
                show_excited_state_combo = true;
                show_lambda_combo = true;
                break;
            case ElectronicStructureType::AttachmentDensity:
            case ElectronicStructureType::DetachmentDensity:
                show_excited_state_combo = true;
                break;
            case ElectronicStructureType::ElectronDensity:
                break;
            default:
                break;
            }

            if (show_mo_combo) {
                md_vlx_scf_type_t type = md_vlx_scf_type(vlx);
                if (type == MD_VLX_SCF_TYPE_RESTRICTED_OPENSHELL) {
                    const char* options[2] = {"Alpha", "Beta"};
                    if (ImGui::BeginCombo("##MO_TYPE", options[export_state.mo.type])) {
                        for (int i = 0; i < ARRAY_SIZE(options); ++i) {
                            if (ImGui::Selectable(options[i])) {
                                export_state.mo.type = (md_vlx_mo_type_t)i;
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
                const double* lambdas = md_vlx_rsp_nto_lambdas(vlx, (size_t)export_state.nto.idx);
                if (lambdas) {
                    char lbl[32];
                    int i = export_state.nto.lambda_idx;
                    snprintf(lbl, sizeof(lbl), "%i (%.4f)", i + 1, lambdas[i]);
                    if (ImGui::BeginCombo("Lambda Idx", lbl)) {
                        for (i = 0; i < 4; i++) {
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
                extent = obb.max_ext - obb.min_ext;
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

                    md_file_o* file = md_file_open(path, MD_FILE_WRITE | MD_FILE_BINARY);
                    if (!file) {
                        MD_LOG_ERROR("Failed to open file for writing: '" STR_FMT "'", STR_ARG(path));
                        return;
                    }

                    const float samples_per_unit_length = volume_resolution_samples_per_angstrom[(int)export_state.resolution] * BOHR_TO_ANGSTROM;
                    md_grid_t grid = {};
                    if (export_state.use_obb) {
                        init_grid(&grid, obb.orientation, obb.min_ext, obb.max_ext, samples_per_unit_length);
                    } else {
                        init_grid(&grid, mat3_ident(), aabb.min_ext, aabb.max_ext, samples_per_unit_length);
                    }

                    defer { md_file_close(file); };

                    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(4));
                    defer { md_vm_arena_destroy(temp_arena); };

                    Volume vol = {};
                    init_volume(&vol, grid, GL_R32F);

                    task_system::ID task = task_system::INVALID_ID;

                    switch (export_state.orb_type) {
                    case ElectronicStructureType::MolecularOrbital:
                    case ElectronicStructureType::MolecularOrbitalDensity:
                    {
                        md_gto_eval_mode_t mode = (export_state.orb_type == ElectronicStructureType::MolecularOrbital) ? MD_GTO_EVAL_MODE_PSI : MD_GTO_EVAL_MODE_PSI_SQUARED;
                        if (use_gpu_path) {
                           compute_mo_GPU(vol.tex_id, grid, export_state.mo.type, export_state.mo.idx, mode);
                        } else {
                            task = compute_mo_async(vol.tex_id, grid, export_state.mo.type, export_state.mo.idx, mode);
                        }
                        break;
                    }
                    case ElectronicStructureType::NaturalTransitionOrbitalParticle:
                    case ElectronicStructureType::NaturalTransitionOrbitalHole:
                    case ElectronicStructureType::NaturalTransitionOrbitalDensityParticle:
                    case ElectronicStructureType::NaturalTransitionOrbitalDensityHole:
                    {
                        md_vlx_nto_type_t type = (export_state.orb_type == ElectronicStructureType::NaturalTransitionOrbitalParticle ||
                                                  export_state.orb_type == ElectronicStructureType::NaturalTransitionOrbitalDensityParticle)
                                                      ? MD_VLX_NTO_TYPE_PARTICLE : MD_VLX_NTO_TYPE_HOLE;
                        md_gto_eval_mode_t mode = (export_state.orb_type == ElectronicStructureType::NaturalTransitionOrbitalDensityParticle ||
                                                   export_state.orb_type == ElectronicStructureType::NaturalTransitionOrbitalDensityHole)
                                                      ? MD_GTO_EVAL_MODE_PSI_SQUARED : MD_GTO_EVAL_MODE_PSI;
                        if (use_gpu_path) {
                            compute_nto_GPU(vol.tex_id, grid, export_state.nto.idx, export_state.nto.lambda_idx, type, mode);
                        } else {
                            task = compute_nto_async(vol.tex_id, grid, export_state.nto.idx, export_state.nto.lambda_idx, type, mode);
                        }
                        break;
                    }
                    case ElectronicStructureType::AttachmentDensity:
                    case ElectronicStructureType::DetachmentDensity:
                    {
                        AttachmentDetachmentType type = (export_state.orb_type == ElectronicStructureType::AttachmentDensity)
                                                            ? AttachmentDetachmentType::Attachment
                                                            : AttachmentDetachmentType::Detachment;
                        if (use_gpu_path) {
                            compute_attachment_detachment_density_GPU(vol.tex_id, grid, export_state.nto.idx, type);
                        } else {
                            task = compute_attachment_detachment_density_async(vol.tex_id, grid, export_state.nto.idx, type);
                        }
                        break;
                    }
                    case ElectronicStructureType::ElectronDensity:
                        if (use_gpu_path) {
                            compute_electron_density_GPU(vol.tex_id, grid, MD_VLX_MO_TYPE_ALPHA);
                        } else {
                            task = compute_electron_density_async(vol.tex_id, grid, MD_VLX_MO_TYPE_ALPHA);
                        }
                        break;
                    default:
                        MD_LOG_ERROR("Unsupported export type");
                        goto done;
                    }

                    int num_samples = vol.dim[0] * vol.dim[1] * vol.dim[2];
                    size_t              natoms = md_vlx_number_of_atoms(vlx);
                    const dvec3_t* vlx_coords  = md_vlx_atom_coordinates(vlx);
                    const uint8_t* vlx_numbers = md_vlx_atomic_numbers(vlx);
                    float* data = (float*)md_vm_arena_push(temp_arena, num_samples * sizeof(float));

                    vec3_t origin = grid.origin;
                    vec3_t step_x = grid.orientation[0] * grid.spacing[0];
                    vec3_t step_y = grid.orientation[1] * grid.spacing[1];
                    vec3_t step_z = grid.orientation[2] * grid.spacing[2];

                    // Extract data from OpenGL Texture
                    glBindTexture(GL_TEXTURE_3D, vol.tex_id);
                    glGetTexImage(GL_TEXTURE_3D, 0, GL_RED, GL_FLOAT, data);
                    glBindTexture(GL_TEXTURE_3D, 0);

                    switch (file_format) {
                    case EXPORT_FILE_FORMAT_CUBE: {
                        md_file_printf(file, "Cube file generated by VIAMD\n");
                        md_file_printf(file, "%s\n", electronic_structure_type_str[(int)export_state.orb_type]);

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

                        str_t mhd_path = str_printf(temp_arena, STR_FMT ".mhd", STR_ARG(basepath));
                        str_t xyz_path = str_printf(temp_arena, STR_FMT ".xyz", STR_ARG(basepath));

                        md_file_o* mhd_file = md_file_open(mhd_path, MD_FILE_WRITE | MD_FILE_BINARY);
                        if (!mhd_file) {
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

                        md_file_close(mhd_file);

                        md_file_o* xyz_file = md_file_open(xyz_path, MD_FILE_WRITE | MD_FILE_BINARY);
                        if (!xyz_file) {
                            MD_LOG_ERROR("Failed to open .xyz file");
                            goto done;
                        }
                        defer { md_file_close(xyz_file); };
                        
                        // XYZ
                        md_file_printf(xyz_file, "%i\n", natoms);
                        md_file_printf(xyz_file, "Geometry extracted from VeloxChem dataset: '%s'\n", state.files.molecule);
                        for (int i = 0; i < natoms; ++i) {
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
        ScopedTemp temp_reset;

        int* donors         = (int*)   md_temp_push_zero(sizeof(int)    * num_groups);
        int* acceptors      = (int*)   md_temp_push_zero(sizeof(int)    * num_groups);
        double* charge_diff = (double*)md_temp_push_zero(sizeof(double) * num_groups);

        int num_donors = 0;
        int num_acceptors = 0;

        double hole_sum = 0;
        double part_sum = 0;

        for (size_t i = 0; i < num_groups; i++) {
            hole_sum += hole_charges[i];
            part_sum += part_charges[i];
        }

        for (size_t i = 0; i < num_groups; i++) {
            // percentages
            double gsCharge = hole_charges[i] / hole_sum;
            double esCharge = part_charges[i] / part_sum;

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
            total_acceptor_charge += charge_diff[acceptors[i]];
        }

        for (int i = 0; i < num_donors; i++) {
            int donor = donors[i];
            double charge_deficit = -charge_diff[donor];
            for (int j = 0; j < num_acceptors; j++) {
                int acceptor = acceptors[j];
                double contrib = charge_deficit * charge_diff[acceptor] / total_acceptor_charge;
                out_matrix[acceptor * num_groups + donor] = (float)contrib;
            }
        }
    }

    struct InteractionCanvasState {
        md_bitfield_t* highlight_mask;
        md_bitfield_t* selection_mask;

        bool out_canvas_pressed;
        bool out_open_context_menu;
    };

    struct SelectionState {
        int hovered_atom_idx;
        int hovered_bond_idx;
        md_bitfield_t* highlight_mask;
        md_bitfield_t* selection_mask;
        SingleSelectionSequence* single_selection_sequence;
        SelectionGranularity granularity;
    };

    struct ViewState {
        const Camera& camera;
        const mat4_t& MVP;
        const TrackballControllerParam& trackball_param;
        const vec3_t& picking_world_coord;
        float picking_depth;
        quat_t* target_ori;
        vec3_t* target_pos;
        float*  target_dist;
    };

    static void interaction_canvas(ImVec2 size, SelectionState& select, ViewState& view, const md_system_t& mol) {
        enum class RegionMode { Append, Remove };

        bool is_selecting = false;

        ScopedTemp temp_reset;
        md_allocator_i* temp_alloc = md_get_temp_allocator();

        bool open_context_menu = false;

        ImGuiWindow* window = ImGui::GetCurrentWindow();
        if (window) {
            bool pressed = ImGui::InvisibleButton("canvas", size, ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_AllowOverlap);
            bool hovered = ImGui::IsItemHovered();
            bool focused = ImGui::IsItemFocused();

            //ImVec2 coord = ImGui::GetMousePos() - ImGui::GetWindowPos();

            ImVec2 canvas_min = ImGui::GetItemRectMin();
            ImVec2 canvas_max = ImGui::GetItemRectMax();
            ImVec2 canvas_size = ImGui::GetItemRectSize();

            ImDrawList* dl = window->DrawList;
            ASSERT(dl);

            ImVec2 win_pos = ImGui::GetWindowPos();

            if (pressed || ImGui::IsItemActive() || ImGui::IsItemDeactivated()) {
                if (ImGui::IsKeyPressed(ImGuiMod_Shift, false)) {
                    ImGui::ResetMouseDragDelta();
                }

                if (ImGui::IsKeyDown(ImGuiMod_Shift)) {
                    RegionMode mode = RegionMode::Append;
                    if (ImGui::IsMouseDown(ImGuiMouseButton_Left) || ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                        mode = RegionMode::Append;
                    }
                    else if (ImGui::IsMouseDown(ImGuiMouseButton_Right) || ImGui::IsMouseReleased(ImGuiMouseButton_Right)) {
                        mode = RegionMode::Remove;
                    }

                    const ImVec2 ext = ImGui::GetMouseDragDelta(mode == RegionMode::Append ? ImGuiMouseButton_Left : ImGuiMouseButton_Right);
                    const ImVec2 pos = ImGui::GetMousePos() - ext;
                    const ImU32 fill_col = 0x22222222;
                    const ImU32 line_col = 0x88888888;

                    // This is relative to the window, clip selection region to canvas
                    ImVec2 sel_min = ImClamp(ImMin(pos, pos + ext), canvas_min, canvas_max);
                    ImVec2 sel_max = ImClamp(ImMax(pos, pos + ext), canvas_min, canvas_max);

                    dl->AddRectFilled(sel_min, sel_max, fill_col);
                    dl->AddRect      (sel_min, sel_max, line_col);

                    // Make it relative to the Canvas
                    sel_min = sel_min - canvas_min;
                    sel_max = sel_max - canvas_min;

                    md_bitfield_t mask = md_bitfield_create(temp_alloc);

                    if (sel_min != sel_max) {
                        md_bitfield_clear(select.highlight_mask);
                        is_selecting = true;

                        // We probably will not use a visibility mask here to filter out only visible atoms
                        // But in the case that we will transition back to that, this code snipped is left here
                        
                        //md_bitfield_iter_t it = md_bitfield_iter_create(&state.representation.visibility_mask);
                        //while (md_bitfield_iter_next(&it)) {
                        //    const uint64_t i = md_bitfield_iter_idx(&it);

                        for (size_t i = 0; i < mol.atom.count; ++i) {
                            const vec4_t p = mat4_mul_vec4(view.MVP, vec4_set(mol.atom.x[i], mol.atom.y[i], mol.atom.z[i], 1.0f));
                            const vec2_t c = {
                                ( p.x / p.w * 0.5f + 0.5f) * canvas_size.x,
                                (-p.y / p.w * 0.5f + 0.5f) * canvas_size.y,
                            };

                            // Test projected point if within selection region
                            if (sel_min.x <= c.x && c.x <= sel_max.x && sel_min.y <= c.y && c.y <= sel_max.y) {
                                md_bitfield_set_bit(&mask, i);
                            }
                        }
                        grow_mask_by_selection_granularity(&mask, select.granularity, mol);

                        if (mode == RegionMode::Append) {
                            md_bitfield_or(select.highlight_mask, select.selection_mask, &mask);
                        }
                        else if (mode == RegionMode::Remove) {
                            md_bitfield_andnot(select.highlight_mask, select.selection_mask, &mask);
                        }
                        if (pressed || ImGui::IsMouseReleased(0)) {
                            md_bitfield_copy(select.selection_mask, select.highlight_mask);
                        }
                    }
                    else if (pressed) {
                        if (select.hovered_atom_idx != -1 || select.hovered_bond_idx != -1) {
                            if (0 <= select.hovered_atom_idx && select.hovered_atom_idx < (int)mol.atom.count) {
                                if (mode == RegionMode::Append) {
                                    single_selection_sequence_push_idx(select.single_selection_sequence, select.hovered_atom_idx);
                                } else {
                                    single_selection_sequence_pop_idx (select.single_selection_sequence, select.hovered_atom_idx);
                                }
                                md_bitfield_set_bit(&mask, select.hovered_atom_idx);
                            } else if (mol.bond.pairs && 0 <= select.hovered_bond_idx && select.hovered_bond_idx < (int)mol.bond.count) {
                                md_atom_pair_t pair = mol.bond.pairs[select.hovered_bond_idx];
                                md_bitfield_set_bit(&mask, pair.idx[0]);
                                md_bitfield_set_bit(&mask, pair.idx[1]);
                            }
                            grow_mask_by_selection_granularity(&mask, select.granularity, mol);
                            if (mode == RegionMode::Append) {
                                md_bitfield_or_inplace(select.selection_mask, &mask);
                            } else {
                                md_bitfield_andnot_inplace(select.selection_mask, &mask);
                            }
                        }
                        else {
                            single_selection_sequence_clear(select.single_selection_sequence);
                            md_bitfield_clear(select.selection_mask);
                            md_bitfield_clear(select.highlight_mask);
                        }
                    }
                }
            }
            else if (ImGui::IsItemHovered() && !ImGui::IsAnyItemActive()) {
                md_bitfield_clear(select.highlight_mask);
                if (select.hovered_atom_idx != -1 || select.hovered_bond_idx != -1) {
                    if (0 <= select.hovered_atom_idx && select.hovered_atom_idx < (int)mol.atom.count) {
                        md_bitfield_set_bit(select.highlight_mask, select.hovered_atom_idx);
                    } else if (mol.bond.pairs && 0 <= select.hovered_bond_idx && select.hovered_bond_idx < (int)mol.bond.count) {
                        md_atom_pair_t pair = mol.bond.pairs[select.hovered_bond_idx];
                        md_bitfield_set_bit(select.highlight_mask, pair.idx[0]);
                        md_bitfield_set_bit(select.highlight_mask, pair.idx[1]);
                    }
                    grow_mask_by_selection_granularity(select.highlight_mask, select.granularity, mol);
                    //draw_info_window(state, select.hovered_atom_idx);
                }
            }

            if (ImGui::IsItemActive() || ImGui::IsItemHovered()) {
                if (!ImGui::IsKeyDown(ImGuiMod_Shift) && !is_selecting) {
                    const ImVec2 delta = ImGui::GetIO().MouseDelta;
                    const ImVec2 coord = ImGui::GetMousePos() - canvas_min;
                    const vec2_t mouse_delta = {delta.x, delta.y};
                    const vec2_t mouse_coord = {coord.x, coord.y};
                    const float  scroll_delta = ImGui::GetIO().MouseWheel;

                    TrackballControllerInput input;
                    input.rotate_button = ImGui::IsMouseDown(ImGuiMouseButton_Left);
                    input.pan_button    = ImGui::IsMouseDown(ImGuiMouseButton_Right);
                    input.dolly_button  = ImGui::IsMouseDown(ImGuiMouseButton_Middle);
                    input.mouse_coord_curr = mouse_coord;
                    input.mouse_coord_prev = mouse_coord - mouse_delta;
                    input.screen_size = {canvas_size.x, canvas_size.y};
                    input.dolly_delta = scroll_delta;
                    input.fov_y = view.camera.fov_y;

                    TrackballFlags flags = TrackballFlags_AnyInteractionReturnsTrue;
                    if (ImGui::IsItemActive()) {
                        flags |= TrackballFlags_EnableAllInteractions;
                    } else {
                        flags |= TrackballFlags_DollyEnabled;
                    }

                    camera_controller_trackball(view.target_pos, view.target_ori, view.target_dist, input, view.trackball_param, flags);
                }
            }
        }
    }

    void delete_group(int index) {
        vec4_t deleted_color = nto.group.color[index];
        for (size_t i = 0; i < nto.num_atoms; i++) {
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
        sprintf(nto.group.label[nto.group.count], "Group %i", (int)nto.group.count);
        nto.group.color[nto.group.count] = deleted_color;
    }

    void copy_group(int from_idx, int to_idx) {
        nto.group.color[to_idx] = nto.group.color[from_idx];
        MEMCPY(nto.group.label[to_idx], nto.group.label[from_idx], sizeof(nto.group.label[from_idx]));
        for (size_t i = 0; i < nto.num_atoms; i++) {
            if (nto.atom_group_idx[i] == to_idx) {
                nto.atom_group_idx[i] = from_idx;
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

        for (size_t i = 0; i < nto.num_atoms; i++) {
            if (nto.atom_group_idx[i] == a) {
                nto.atom_group_idx[i] = b;
            }
            else if (nto.atom_group_idx[i] == b) {
                nto.atom_group_idx[i] = a;
            }
        }
    }

    static inline void imgui_delayed_hover_tooltip(const char* text) {
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayNormal)) {
            if (ImGui::BeginTooltip()) {
                ImGui::Text(text);
                ImGui::EndTooltip();
            }
        }
    }

    void draw_nto_window(ApplicationState& state) {
        if (!nto.show_window) return;

        size_t num_excited_states = md_vlx_rsp_number_of_excited_states(vlx);

        if (num_excited_states == 0) return;

        bool open_context_menu = false;
        static bool edit_mode = false;

        ImGui::SetNextWindowSize(ImVec2(500, 300), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Transition Analysis", &nto.show_window, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
            nto.group.hovered_index = -1;
            bool viewport_hovered = false;

            if (ImGui::BeginMenuBar()) {
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
            ImGui::Text((const char*)u8"Iso Value (Ψ²)"); 
            ImGui::SliderScalar("##Iso Value", ImGuiDataType_Double, &nto.iso_val, &iso_min, &iso_max, "%.6f", ImGuiSliderFlags_Logarithmic);
            ImGui::SetItemTooltip("Iso Value");

            ImGui::Spacing();

            static const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersH | ImGuiTableFlags_SizingFixedFit;
            static const ImGuiTableColumnFlags columns_base_flags = ImGuiTableColumnFlags_NoSort;
            
            static bool hide_overlap_text = true;
            bool refresh = false;
            ImVec2 button_size(ImGui::GetFontSize() + ImGui::GetStyle().FramePadding.x, 0.f);

            ImGui::Checkbox("Edit mode", &edit_mode);
            ImGui::SameLine();
            ImGui::Checkbox("Hide overlapping text", &hide_overlap_text);

            int group_counts[MAX_NTO_GROUPS] = {0};
            for (size_t i = 0; i < nto.num_atoms; ++i) {
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
                bool item_hovered = false;
                bool show_unassigned = group_counts[0] > 0;
                int row_n = show_unassigned ? 0 : 1;
                for (; row_n < nto.group.count; row_n++) {
                    ImGui::PushID(row_n);
                    ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                    ImGui::TableNextRow(ImGuiTableRowFlags_None, 0);

                    ImGui::TableNextColumn();
                    ImGui::AlignTextToFramePadding();
                    if (edit_mode) {
                        ImGui::PushItemWidth(150.f);
                        if (row_n > 0) {
                            ImGui::InputText("##label", nto.group.label[row_n], sizeof(nto.group.label[row_n]));
                        } else {
                            ImGui::Text(nto.group.label[0]);
                        }
                        ImGui::PopItemWidth();
                    }
                    else {
                        ImGui::Selectable(nto.group.label[row_n], false, selectable_flags);
                        if (ImGui::TableGetHoveredRow() == row_n + show_unassigned) {
                            nto.group.hovered_index = (int8_t)row_n;
                            md_bitfield_clear(&state.selection.highlight_mask);
                            for (size_t j = 0; j < nto.num_atoms; j++) {
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
                            if (row_n == nto.group.count - 1) {
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
                                        for (int i = last; i >= row_n - 1; i--) {
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
                                    for (size_t i = 0; i < nto.num_atoms; i++) {
                                        if (nto.atom_group_idx[i] == row_n) {
                                            nto.atom_group_idx[i] = row_n - 1;
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
            camera_animate(&nto.camera, nto.target.ori, nto.target.pos, nto.target.dist, dt);

            ImVec2 canvas_sz = ImGui::GetContentRegionAvail();   // Resize canvas to what's available
            canvas_sz.x = MAX(canvas_sz.x, 50.0f);
            canvas_sz.y = MAX(canvas_sz.y, 50.0f);

            ImGui::Dummy(canvas_sz);

            ImVec2 cursor = ImGui::GetCursorPos();

            // Draw border and background color
            ImGuiIO& io = ImGui::GetIO();

            ImVec2 canvas_p0 = ImGui::GetItemRectMin();
            ImVec2 canvas_p1 = ImGui::GetItemRectMax();

            bool reset_view = false;

            if (rsp.selected != -1) {
                // This represents the cutoff for contributing orbitals to be part of the orbital 'grid'
                // If the occupation parameter is less than this it will not be displayed
                const double* lambda = md_vlx_rsp_nto_lambdas(vlx, rsp.selected);
                if (!lambda) {
                    MD_LOG_ERROR("No lambda information available for NTOs in veloxchem object");
                } else {
                    if (nto.sel_nto_idx != rsp.selected) {
                        nto.sel_nto_idx = rsp.selected;
                        size_t nto_idx = (size_t)rsp.selected;

                        size_t num_lambdas = 0;
                        for (size_t i = 0; i < MAX_NTO_LAMBDAS; ++i) {
                            if (lambda[i] < NTO_LAMBDA_CUTOFF_VALUE) {
                                break;
							}
							num_lambdas++;
                        }

                        init_volume(&nto.vol[NTO_Attachment], nto.grid, GL_R32F);
                        init_volume(&nto.vol[NTO_Detachment], nto.grid, GL_R32F);

                        if (use_gpu_path) {
                            compute_attachment_detachment_density_GPU(nto.vol[NTO_Attachment].tex_id, nto.grid, nto_idx, AttachmentDetachmentType::Attachment);
                            compute_attachment_detachment_density_GPU(nto.vol[NTO_Detachment].tex_id, nto.grid, nto_idx, AttachmentDetachmentType::Detachment);
                        } else {
                            nto.vol_task[NTO_Attachment] = compute_attachment_detachment_density_async(nto.vol[NTO_Attachment].tex_id, nto.grid, nto_idx, AttachmentDetachmentType::Attachment);
                            nto.vol_task[NTO_Detachment] = compute_attachment_detachment_density_async(nto.vol[NTO_Detachment].tex_id, nto.grid, nto_idx, AttachmentDetachmentType::Detachment);
                        }

                        for (size_t lambda_idx = 0; lambda_idx < num_lambdas; ++lambda_idx) {
                            size_t pi = (size_t)NTO_Part_0 + lambda_idx;
                            size_t hi = (size_t)NTO_Hole_0 + lambda_idx;

                            if (task_system::task_is_running(nto.vol_task[pi])) {
                                task_system::task_interrupt (nto.vol_task[pi]);
                            }
                            if (task_system::task_is_running(nto.vol_task[hi])) {
                                task_system::task_interrupt (nto.vol_task[hi]);
                            }

                            init_volume(&nto.vol[pi], nto.grid);
                            init_volume(&nto.vol[hi], nto.grid);

                            if (use_gpu_path) {
                                compute_nto_GPU(nto.vol[pi].tex_id, nto.grid, nto_idx, lambda_idx, MD_VLX_NTO_TYPE_PARTICLE, MD_GTO_EVAL_MODE_PSI);
                                compute_nto_GPU(nto.vol[hi].tex_id, nto.grid, nto_idx, lambda_idx, MD_VLX_NTO_TYPE_HOLE,     MD_GTO_EVAL_MODE_PSI);
                            } else {
                                nto.vol_task[pi] = compute_nto_async(nto.vol[pi].tex_id, nto.grid, nto_idx, lambda_idx, MD_VLX_NTO_TYPE_PARTICLE, MD_GTO_EVAL_MODE_PSI);
                                nto.vol_task[hi] = compute_nto_async(nto.vol[hi].tex_id, nto.grid, nto_idx, lambda_idx, MD_VLX_NTO_TYPE_HOLE,     MD_GTO_EVAL_MODE_PSI);
                            }

                        }
                        if (task_system::task_is_running(nto.seg_task[0])) {
                            task_system::task_interrupt(nto.seg_task[0]);
                        }
                        if (task_system::task_is_running(nto.seg_task[1])) {
                            task_system::task_interrupt(nto.seg_task[1]);
                        }
                    }
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
            mat4_t proj_mat     = camera_perspective_projection_matrix(nto.camera, aspect_ratio);
            mat4_t inv_view_mat = camera_view_to_world_matrix(nto.camera);
            mat4_t inv_proj_mat = camera_inverse_perspective_projection_matrix(nto.camera, aspect_ratio);

            mat4_t MVP      = proj_mat * view_mat;
            mat4_t inv_MVP  = inv_view_mat * inv_proj_mat;

            if (rsp.selected != -1) {
                SelectionState selection = {
                    .hovered_atom_idx = state.selection.atom_idx.hovered,
                    .hovered_bond_idx = state.selection.bond_idx.hovered,
                    .highlight_mask = &state.selection.highlight_mask,
                    .selection_mask = &state.selection.selection_mask,
                    .single_selection_sequence = &state.selection.single_selection_sequence,
                    .granularity = state.selection.granularity,
                };

                ViewState view = {
                    .camera = nto.camera,
                    .MVP = MVP,
                    .trackball_param = state.view.trackball_param,
                    .picking_world_coord = nto.picking.world_coord,
                    .picking_depth = nto.picking.depth,
                    .target_ori = &nto.target.ori,
                    .target_pos = &nto.target.pos,
                    .target_dist = &nto.target.dist,
                };

                ImRect hovered_canvas_rect = {};

                const int nto_target_idx[2] = {
                    NTO_Attachment,
                    NTO_Detachment,
                };

                // Draw Attachment / Detachment
                for (int i = 0; i < 2; ++i) {
                    int idx = nto_target_idx[i];

                    ImVec2 p0 = grid_p0 + win_sz * ImVec2(0.0f, (float)(i+0));
                    ImVec2 p1 = grid_p0 + win_sz * ImVec2(1.0f, (float)(i+1));
                    ImRect rect = {p0, p1};
                    ImGui::SetCursorScreenPos(p0);

                    draw_list->PushClipRect(p0, p1);
                    defer{ draw_list->PopClipRect(); };

                    draw_list->ChannelsSetCurrent(1);
                    ImGui::PushID(i);
                    interaction_canvas(p1-p0, selection, view, state.mold.sys);

                    if (ImGui::IsItemHovered()) {
                        viewport_hovered = true;
                        if (ImGui::GetIO().MouseDoubleClicked[0]) {
                            if (view.picking_depth < 1.0f) {
                                const vec3_t forward = view.camera.orientation * vec3_t{0, 0, 1};
                                nto.target.pos = view.picking_world_coord + forward * *view.target_dist;
                            } else {
                                reset_view = true;
                            }
                        }

                        if (!ImGui::IsKeyDown(ImGuiMod_Shift) && ImGui::IsMouseReleased(ImGuiMouseButton_Right) && ImGui::GetMouseDragDelta(ImGuiMouseButton_Right) == ImVec2(0,0)) {
                            open_context_menu = true;
                        }
                        hovered_canvas_rect = rect;
                    }

                    ImGui::PopID();

                    draw_list->ChannelsSetCurrent(0);

                    ImVec2 text_pos_bl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p1.y - TEXT_BASE_HEIGHT);
                    ImVec2 text_pos_tl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p0.y + TEXT_BASE_HEIGHT * 0.5f);
                    const char* lbl = ((i & 1) == 0) ? "Attachment" : "Detachment";

#if DEBUG
                    char lbls[64];
                    snprintf(lbls, sizeof(lbls), "hovered_atom_idx: %i, hovered_bond_idx: %i", state.selection.atom_idx.hovered, state.selection.bond_idx.hovered);
#endif

                    //char buf[32];
                    //snprintf(buf, sizeof(buf), (const char*)u8"λ: %.3f", nto_lambda[i / 2]);
                    draw_list->AddImage((ImTextureID)(intptr_t)nto.gbuf.tex.transparency, p0, p1, { 0,1 }, { 1,0 });
                    draw_list->AddImage((ImTextureID)(intptr_t)nto.iso_tex[idx], p0, p1, { 0,1 }, { 1,0 });
                    //draw_list->AddText(text_pos_bl, ImColor(0,0,0), buf);
                    draw_list->AddText(text_pos_tl, ImColor(0,0,0), lbl);

                    if (nto.dipole.enabled) {
                        const vec3_t mid = vec3_lerp(aabb.min_ext, aabb.max_ext, 0.5f);
                        const vec3_t ext = aabb.max_ext - aabb.min_ext;

                        const dvec3_t* magnetic_dp = md_vlx_rsp_magnetic_transition_dipole_moments(vlx);
                        const dvec3_t* electric_dp = md_vlx_rsp_electric_transition_dipole_moments(vlx);

                        vec3_t magn_vec = {(float)magnetic_dp[rsp.selected].x, (float)magnetic_dp[rsp.selected].y, (float)magnetic_dp[rsp.selected].z};
                        vec3_t elec_vec = {(float)electric_dp[rsp.selected].x, (float)electric_dp[rsp.selected].y, (float)electric_dp[rsp.selected].z};

                        magn_vec *= (float)(nto.dipole.vector_scale * BOHR_TO_ANGSTROM);
                        elec_vec *= (float)(nto.dipole.vector_scale * BOHR_TO_ANGSTROM);

                        auto proj_point = [MVP, win_sz, p0](vec3_t point) -> ImVec2 {
                            const vec4_t p = mat4_mul_vec4(MVP, vec4_from_vec3(point, 1.0));
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
                                float o_scl = 0.57735026918962576451 * d_scl;
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
                const bool is_active  = ImGui::IsItemActive();
                const ImVec2 origin(canvas_p0.x, canvas_p0.y);  // Lock scrolled origin
                const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

                int width  = MAX(1, (int)win_sz.x);
                int height = MAX(1, (int)win_sz.y);

                const int num_win = 2;

                auto& gbuf = nto.gbuf;
                if ((int)gbuf.width != width || (int)gbuf.height != height) {
                    init_gbuffer(&gbuf, width, height);
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
                    camera_compute_optimal_view(&nto.target.pos, &nto.target.ori, &nto.target.dist, obb.orientation, obb.min_ext * BOHR_TO_ANGSTROM, obb.max_ext * BOHR_TO_ANGSTROM, nto.distance_scale);
                }

                if (nto.show_coordinate_system_widget) {
                    float  ext = MIN(win_sz.x, win_sz.y) * 0.4f;
                    float  pad = 20.0f;

                    ImVec2 min = ImGui::GetItemRectMin() - ImGui::GetWindowPos();
                    ImVec2 max = ImGui::GetItemRectMax() - ImGui::GetWindowPos();

                    CoordSystemWidgetParam param = {
                        .pos = ImVec2(min.x + pad, max.y - ext - pad),
                        .size = {ext, ext},
                        .view_matrix = camera_world_to_view_matrix(nto.camera),
                        .camera_ori  = nto.target.ori,
                        .camera_pos  = nto.target.pos,
                        .camera_dist = nto.target.dist,
                    };

                    ImGui::DrawCoordinateSystemWidget(param);
                }

                clear_gbuffer(&gbuf);

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

                auto draw_rep = [](md_gl_rep_t& rep, md_gl_shaders_t& shaders, mat4_t& view_mat, mat4_t& proj_mat, uint32_t atom_mask = 0) {
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
                        .atom_mask = atom_mask,
                    };

                    md_gl_draw(&draw_args);
                };

                draw_rep(nto.gl_rep, state.mold.gl_shaders, view_mat, proj_mat);

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
                        draw_rep(nto.gl_rep, state.mold.gl_shaders_lean_and_mean, view_mat, proj_mat, AtomBit_Selected);

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
                        draw_rep(nto.gl_rep, state.mold.gl_shaders_lean_and_mean, view_mat, proj_mat, AtomBit_Highlighted);

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
                postprocessing::Descriptor postprocess_desc = {
                    .background = {
                        .color = state.visuals.background.color * state.visuals.background.intensity,
                    },
                    .tonemapping = {
                        .enabled    = state.visuals.tonemapping.enabled,
                        .mode       = state.visuals.tonemapping.tonemapper,
                        .exposure   = state.visuals.tonemapping.exposure,
                        .gamma      = state.visuals.tonemapping.gamma,
                    },
                    .ambient_occlusion = {
                        .enabled = false
                    },
                    .depth_of_field = {
                        .enabled = false,
                    },
                    .fxaa = {
                        .enabled = true,
                    },
                    .temporal_aa = {
                        .enabled = false,
                    },
                    .sharpen = {
                        .enabled = false,
                    },
                    .input_textures = {
                        .depth          = nto.gbuf.tex.depth,
                        .color          = nto.gbuf.tex.color,
                        .normal         = nto.gbuf.tex.normal,
                        .velocity       = nto.gbuf.tex.velocity,
                        .transparency   = nto.gbuf.tex.transparency,
                    }
                };

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

                postprocessing::shade_and_postprocess(postprocess_desc, view_param);
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
                        values[0] = nto.iso_val * nto.iso_val;
                    } else {
                        colors[0] = nto.col_pos;
                        colors[1] = nto.col_neg;
                        values[0] = nto.iso_val;
                        values[1] = -nto.iso_val;
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
                                .volume = nto.vol[i].tex_id,
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

                if (hovered_canvas_rect.GetArea() > 0) {
                    ImVec2 coord = ImGui::GetMousePos() - hovered_canvas_rect.Min;
                    coord.y = hovered_canvas_rect.GetSize().y - coord.y;
                    extract_picking_data(nto.picking, nto.gbuf, {coord.x, coord.y}, inv_MVP);

                    state.selection.atom_idx.hovered = INVALID_PICKING_IDX;
                    state.selection.bond_idx.hovered = INVALID_PICKING_IDX;

                    if (nto.picking.idx != INVALID_PICKING_IDX) {
                        // The index space is segmented into two parts, the first half is for atoms and the second half is for bonds
                        if (nto.picking.idx < 0x80000000) {
                            state.selection.atom_idx.hovered = nto.picking.idx;
                        } else {
                            state.selection.bond_idx.hovered = nto.picking.idx & 0x7FFFFFFF;
                        }
                    }
                }
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
            if (ImGui::BeginMenu("Assign to group")) {
                for (size_t i = 0; i < nto.group.count; i++) {
                    if (ImGui::MenuItem(nto.group.label[i])) {
                        for (size_t j = 0; j < nto.num_atoms; j++) {
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
                    for (size_t i = 0; i < nto.num_atoms; i++) {
                        if (md_bitfield_test_bit(&state.selection.selection_mask, i)) {
                            //If so, set its group to the selected group index
                            nto.atom_group_idx[i] = group_idx;
                        }
                    }
                }
            }
            ImGui::EndPopup();
        }

        // Create hash to check for changes to trigger update of colors in gl representation
        uint64_t atom_idx_hash = md_hash64(nto.atom_group_idx, sizeof(uint32_t) * nto.num_atoms, 0);
        uint64_t color_hash    = md_hash64(nto.group.color, sizeof(nto.group.color), atom_idx_hash);

        static uint64_t cur_color_hash = 0;
        if (color_hash != cur_color_hash) {
            cur_color_hash = color_hash;
            update_nto_group_colors();
        }

        // Create hash to check for changes to trigger recomputation of transition matrix
        uint64_t matrix_hash = atom_idx_hash ^ nto.sel_nto_idx ^ nto.group.count;
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

            if (nto.sel_nto_idx != -1) {
                const size_t nto_idx = (size_t)nto.sel_nto_idx;
#if 1
                {
                    ScopedTemp temp_reset;

				    size_t num_aos = md_vlx_scf_number_of_atomic_orbitals(vlx);
					size_t num_lambdas = 4;

                    // Calculate transition densities for groups
				    const double* lambdas = md_vlx_rsp_nto_lambdas(vlx, nto_idx);
				    const double* S = md_vlx_scf_overlap_matrix_data(vlx);
				    const double* a_part[4] = {0};
				    const double* a_hole[4] = {0};

                    for (size_t i = 0; i < num_lambdas; i++) {
                        double* part = (double*)md_temp_push(sizeof(double) * num_aos);
                        double* hole = (double*)md_temp_push(sizeof(double) * num_aos);
                        md_vlx_rsp_nto_extract_coefficients(part, vlx, nto_idx, i, MD_VLX_NTO_TYPE_PARTICLE);
                        md_vlx_rsp_nto_extract_coefficients(hole, vlx, nto_idx, i, MD_VLX_NTO_TYPE_HOLE);
                        a_part[i] = part;
                        a_hole[i] = hole;
                    }

                    double group_density_part[MAX_NTO_GROUPS] = {0};
                    double group_density_hole[MAX_NTO_GROUPS] = {0};

                    const int* ao_to_atom = md_vlx_ao_to_atom_idx(vlx);

					int* ao_to_group = (int*)md_temp_push(sizeof(int) * num_aos);

                    for (size_t i = 0; i < num_aos; i++) {
                        int atom_idx = ao_to_atom[i];
						ASSERT(0 <= atom_idx && (size_t)atom_idx < nto.num_atoms);
                        ao_to_group[i] = (int)nto.atom_group_idx[atom_idx];
                        if (ao_to_group[i] < 0 || ao_to_group[i] >= MAX_NTO_GROUPS) {
                            ao_to_group[i] = 0;
                        }
					}

					attribute_charge(group_density_part, ao_to_group, lambdas, a_part, S, num_lambdas, num_aos);
					attribute_charge(group_density_hole, ao_to_group, lambdas, a_hole, S, num_lambdas, num_aos);
                    for (size_t g1 = 0; g1 < nto.group.count; g1++) {
                        nto.transition_density_part[g1] = (float)group_density_part[g1];
                        nto.transition_density_hole[g1] = (float)group_density_hole[g1];
					}
					compute_transition_matrix(nto.transition_matrix, nto.group.count, nto.transition_density_hole, nto.transition_density_part);
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

                    if (compute_transition_group_values_async(&eval_attach, &seg_attach, nto.transition_density_part, nto.group.count, nto.grid, nto.atom_group_idx, nto.atom_xyzr, nto.num_atoms, nto_idx, MD_VLX_NTO_TYPE_PARTICLE, MD_GTO_EVAL_MODE_PSI_SQUARED, samples_per_unit_length) &&
                        compute_transition_group_values_async(&eval_detach, &seg_detach, nto.transition_density_hole, nto.group.count, nto.grid, nto.atom_group_idx, nto.atom_xyzr, nto.num_atoms, nto_idx, MD_VLX_NTO_TYPE_HOLE,     MD_GTO_EVAL_MODE_PSI_SQUARED, samples_per_unit_length))
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
