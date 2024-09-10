#define IMGUI_DEFINE_MATH_OPERATORS

#include <event.h>
#include <viamd.h>
#include <task_system.h>
#include <color_utils.h>

#include <md_gto.h>
#include <md_vlx.h>
#include <md_util.h>
#include <core/md_vec_math.h>
#include <core/md_log.h>
#include <core/md_arena_allocator.h>

#include <gfx/volumerender_utils.h>
#include <gfx/gl_utils.h>
#include <gfx/immediate_draw_utils.h>

#include <imgui_internal.h>
#include <imgui_widgets.h>
#include <implot_widgets.h>

#include <md_csv.h>
#include <md_xvg.h>

#define BLK_DIM 8
#define ANGSTROM_TO_BOHR 1.8897261246257702
#define BOHR_TO_ANGSTROM 0.529177210903
#define DEFAULT_SAMPLES_PER_ANGSTROM 6
#define MAX_GROUPS 16

#define IM_GREEN ImVec4{0, 1, 0, 1}
#define IM_RED ImVec4{1, 0, 0, 1}
#define IM_YELLOW ImVec4{1, 1, 0.5, 0.3}
#define IM_BLUE ImVec4{0.5, 0.5, 1, 0.3}

#define U32_MAGENTA IM_COL32(255, 0, 255, 255)

inline const ImVec4& vec_cast(const vec4_t& v) { return *(const ImVec4*)(&v); }
inline const vec4_t& vec_cast(const ImVec4& v) { return *(const vec4_t*)(&v); }
inline const ImVec2& vec_cast(const vec2_t& v) { return *(const ImVec2*)(&v); }
inline const vec2_t& vec_cast(const ImVec2& v) { return *(const vec2_t*)(&v); }

inline ImVec4& vec_cast(vec4_t& v) { return *(ImVec4*)(&v); }
inline vec4_t& vec_cast(ImVec4& v) { return *(vec4_t*)(&v); }
inline ImVec2& vec_cast(vec2_t& v) { return *(ImVec2*)(&v); }
inline vec2_t& vec_cast(ImVec2& v) { return *(vec2_t*)(&v); }

enum class VolumeRes {
    Low,
    Mid,
    High,
    Count,
};

static const float vol_res_scl[3] = {
    4.0f,
    8.0f,
    16.0f,
};

struct OBB {
    mat3_t basis;
    vec3_t min_ext;
    vec3_t max_ext;
};

// Compute an oriented bounding box (OBB) for the supplied PGTOS
OBB compute_pgto_obb(const mat3_t& PCA, const md_gto_t* pgtos, size_t num_pgtos) {
    mat4_t Ri  = mat4_from_mat3(PCA);

    // Compute min and maximum extent along the PCA axes
    vec3_t min_ext = { FLT_MAX, FLT_MAX, FLT_MAX};
    vec3_t max_ext = {-FLT_MAX,-FLT_MAX,-FLT_MAX};        

    // Transform the pgtos (x,y,z,cutoff) into the PCA frame to find the min and max extend within it
    for (size_t i = 0; i < num_pgtos; ++i) {
        vec3_t xyz = { pgtos[i].x, pgtos[i].y, pgtos[i].z };

        // The 0.9 scaling factor here is a bit arbitrary, but the cutoff-radius is computed on a value which is lower than the rendered iso-value
        // So the effective radius is a bit overestimated and thus we scale it back a bit
        float  r = pgtos[i].cutoff * 0.9f;

        vec3_t p = mat4_mul_vec3(Ri, xyz, 1.0f);
        min_ext = vec3_min(min_ext, vec3_sub_f(p, r));
        max_ext = vec3_max(max_ext, vec3_add_f(p, r));
    }

    OBB obb = {
        .basis = mat3_transpose(PCA),
        .min_ext = min_ext,
        .max_ext = max_ext,
    };

    return obb;
}

mat4_t compute_vol_mat(const OBB& obb) {
    mat4_t T = mat4_translate_vec3(obb.basis * obb.min_ext * BOHR_TO_ANGSTROM);
    mat4_t R = mat4_from_mat3(obb.basis);
    mat4_t S = mat4_scale_vec3((obb.max_ext - obb.min_ext) * BOHR_TO_ANGSTROM);
    return T * R * S;
}

// Voronoi segmentation
void grid_segment_and_attribute(float* out_group_values, size_t group_cap, const uint8_t* point_group_idx, const vec3_t* point_xyz, const float* point_r, size_t num_points, const md_grid_t* grid) {
    for (int iz = 0; iz < grid->dim[2]; ++iz) {
        for (int iy = 0; iy < grid->dim[1]; ++iy) {
            for (int ix = 0; ix < grid->dim[0]; ++ix) {
                int index = ix + iy * grid->dim[0] + iz * grid->dim[0] * grid->dim[1];
                float value = grid->data[index];

				// Skip if its does not contribute
                if (value == 0.0f) continue;

                float x = grid->origin[0] + ix * grid->step_x[0] + iy * grid->step_y[0] + iz * grid->step_z[0];
                float y = grid->origin[1] + ix * grid->step_x[1] + iy * grid->step_y[1] + iz * grid->step_z[1];
                float z = grid->origin[2] + ix * grid->step_x[2] + iy * grid->step_y[2] + iz * grid->step_z[2];

                float  min_dist = FLT_MAX;
                size_t group_idx = 0;

                // find closest point to grid point
                for (size_t i = 0; i < num_points; ++i) {
                    float px = point_xyz[i].x;
                    float py = point_xyz[i].y;
                    float pz = point_xyz[i].z;
                    float pr = point_r[i];

                    float dx = px - x;
                    float dy = py - y;
                    float dz = pz - z;

                    float dist = dx*dx + dy*dy + dz*dz - (pr*pr);
                    if (dist < min_dist) {
                        min_dist = dist;
                        group_idx = point_group_idx[i];
                    }
                }

                if (group_idx < group_cap) {
                    out_group_values[group_idx] += value;
                }
            }
        }
    }
}

struct VeloxChem : viamd::EventHandler {
    VeloxChem() { viamd::event_system_register_handler(*this); }

    md_vlx_data_t vlx {};

    // Used for clearing volumes
    uint32_t vol_fbo = 0;

    // GL representations
    //md_gl_mol_t gl_mol = {};
    md_gl_rep_t gl_rep = {};

    int homo_idx = 0;
    int lumo_idx = 0;

    // Principal Component Axes of the geometry
    mat3_t PCA = mat3_ident();
    vec3_t com = {};
    vec3_t min_aabb = {};
    vec3_t max_aabb = {};

    struct Volume {
        mat4_t tex_to_world = {};
        int dim[3] = {128, 128, 128};
        vec3_t step_size = {};
        vec3_t extent = {};
        uint32_t tex_id = 0;
    };

    struct Scf {
        bool show_window = false;
    } scf;

    struct Orb {
        bool show_window = false;
        Volume   vol[16] = {};
        int      vol_mo_idx[16] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
        uint32_t iso_tex[16] = {};
        task_system::ID vol_task[16] = {};
        int num_x  = 3;
        int num_y  = 3;
        int mo_idx = -1;
        int scroll_to_idx = -1;

        struct {
            bool enabled = true;
            size_t count = 2;
            float  values[2] = {0.05f, -0.05};
            vec4_t colors[2] = {{215.f/255.f,25.f/255.f,28.f/255.f,0.75f}, {44.f/255.f,123.f/255.f,182.f/255.f,0.75f}};

        } iso;

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
        // In practice I don't think more than 2, maybe 3 will be used in practice
        Volume   vol[8] = {};
        uint32_t iso_tex[8] = {};
        task_system::ID vol_task[8] = {};
        task_system::ID seg_task[2] = {};
        int vol_nto_idx = -1;

        size_t num_atoms = 0;
        // These are in Atomic Units (Bohr)
        vec3_t* atom_xyz = nullptr;
        float*  atom_r   = nullptr;
        uint8_t* atom_group_idx = nullptr;

        md_gl_rep_t gl_rep = {0};

        // Square transition matrix, should have dim == num_groups
        size_t transition_matrix_dim = 0;
        float* transition_matrix = nullptr;
        float* transition_density_hole = nullptr;
        float* transition_density_part = nullptr;

        struct {
            size_t count = 0;
            char   label[MAX_GROUPS][64] = {};
            vec4_t color[MAX_GROUPS] = {};
        } group;

        struct {
            bool enabled = true;
            size_t count = 2;
            float  values[2] = {0.05f, -0.05};
            vec4_t colors[2] = {{215.f/255.f,25.f/255.f,28.f/255.f,0.75f}, {44.f/255.f,123.f/255.f,182.f/255.f,0.75f}};
            vec4_t vectorColors[3] = {{0.f / 255.f, 255.f / 255.f, 255.f / 255.f, 1.f},
                                              {255.f / 255.f, 0.f / 255.f, 255.f / 255.f, 1.f},
                                              {255.f / 255.f, 255.f / 255.f, 0.f / 255.f, 1.f}};
            double vector_length = 1.0f;
            bool display_vectors = false;
            bool display_angle = false;
        } iso;

        GBuffer gbuf = {};
        Camera camera = {};
        PickingData picking = {};

        struct {
            quat_t ori = {};
            vec3_t pos = {};
            float dist = {};
        } target;

        float distance_scale = 2.0f;

        bool show_coordinate_system_widget = true;
    } nto;

    struct Rsp {
        bool show_window = false;
        bool show_export_window = false;
        int hovered = -1;
        int selected = -1;
        int focused_plot = -1;

        double* x_ev_samples;
        double* x_unit_samples;
        double* x_unit_peaks;
        //Spectra Y values, calculated from osc
        double* eps;
        //Spectra y values, calculated from rot
        double* ecd;
        double* vib_y;
        double* vib_x;
        double* vib_points;
        double* osc_points;
        double* cgs_points;
        const char* x_unit;

        bool first_plot_rot_ecd = true;
        bool first_plot_vib = true;
    } rsp;

    // Arena for persistent allocations for the veloxchem module (tied to the lifetime of the VLX object)
    md_allocator_i* arena = 0;

    size_t num_orbitals() const {
        return vlx.scf.alpha.orbitals.dim[1];
    }

    size_t num_cgtos() const {
        return vlx.scf.alpha.orbitals.dim[0];
    }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t event_idx = 0; event_idx < num_events; ++event_idx) {
            const viamd::Event& e = events[event_idx];

            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState& state = *(ApplicationState*)e.payload;
                arena = md_arena_allocator_create(state.allocator.persistent, MEGABYTES(1));
                break;
            }
            case viamd::EventType_ViamdShutdown:
                md_arena_allocator_destroy(arena);
                break;
            case viamd::EventType_ViamdFrameTick: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState& state = *(ApplicationState*)e.payload;

                draw_orb_window(state);
                draw_nto_window(state);
                draw_summary_window(state);
                draw_rsp_window(state);
                break;
            }
            case viamd::EventType_ViamdDrawMenu:
                ImGui::Checkbox("VeloxChem Summary", &scf.show_window);
                ImGui::Checkbox("VeloxChem RSP", &rsp.show_window);
                ImGui::Checkbox("VeloxChem ORB", &orb.show_window);
                ImGui::Checkbox("VeloxChem NTO", &nto.show_window);
                break;
            case viamd::EventType_ViamdRenderTransparent: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                //ApplicationState& state = *(ApplicationState*)e.payload;
                //draw_orb_volume(state);
                break;
            }
            case HASH_STR_LIT("Secret Sauce"): {
                struct Payload {
                    ApplicationState* state;
                    str_t filename;
                };
                Payload* payload = (Payload*)e.payload;
                init_from_file(payload->filename, *payload->state);
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

            case viamd::EventType_RepresentationInfoFill: {
                ASSERT(e.payload_type == viamd::EventPayloadType_RepresentationInfo);
                RepresentationInfo& info = *(RepresentationInfo*)e.payload;

                info.mo_homo_idx = homo_idx;
                info.mo_lumo_idx = lumo_idx;

                for (size_t i = 0; i < num_orbitals(); ++i) {
                    MolecularOrbital mo = {
                        .idx = (int)i,
                        .occupation = (float)vlx.scf.alpha.energies.data[i],
                        .energy = (float)vlx.scf.alpha.energies.data[i],
                    };
                    md_array_push(info.molecular_orbitals, mo, info.alloc);
                }
                
                auto push_dipole = [&info](md_vlx_dipole_moment_t vlx_dp) {
                    DipoleMoment dp = {
                        .label = str_copy(vlx_dp.ident, info.alloc),
                        .vector = vec3_set((float)vlx_dp.x, (float)vlx_dp.y, (float)vlx_dp.z),
                    };
                    md_array_push(info.dipole_moments, dp, info.alloc);
                };

                push_dipole(vlx.scf.ground_state_dipole_moment);
                for (size_t i = 0; i < vlx.rsp.num_excited_states; ++i)
                    push_dipole(vlx.rsp.electronic_transition_length[i]);
                for (size_t i = 0; i < vlx.rsp.num_excited_states; ++i)
                    push_dipole(vlx.rsp.electronic_transition_velocity[i]);
                for (size_t i = 0; i < vlx.rsp.num_excited_states; ++i)
                    push_dipole(vlx.rsp.magnetic_transition[i]);
                break;
            }
            case viamd::EventType_RepresentationComputeOrbital: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ComputeOrbital);
                ComputeOrbital& data = *(ComputeOrbital*)e.payload;

                if (!data.output_written) {
                    md_gto_eval_mode_t mode = MD_GTO_EVAL_MODE_PSI;
                    if (data.type == OrbitalType::PsiSquared) {
                       mode = MD_GTO_EVAL_MODE_PSI_SQUARED;
                    }
                    task_system::ID id = compute_mo_async(&data.tex_mat, &data.voxel_spacing, data.dst_texture, data.orbital_idx, mode, data.samples_per_angstrom);
                    data.output_written = (id != task_system::INVALID_ID);
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
        md_arena_allocator_reset(arena);
        vlx = {};
        orb = {};
        nto = {};
        rsp = {};
    }

    void init_from_file(str_t filename, ApplicationState& state) {
        str_t ext;
        if (extract_ext(&ext, filename) && str_eq_ignore_case(ext, STR_LIT("out"))) {
            MD_LOG_INFO("Attempting to load VeloxChem data from file '" STR_FMT "'", STR_ARG(filename));
            md_vlx_data_free(&vlx);
            if (md_vlx_data_parse_file(&vlx, filename, arena)) {
                MD_LOG_INFO("Successfully loaded VeloxChem data");

                if (!vol_fbo) glGenFramebuffers(1, &vol_fbo);

                // Scf
                scf.show_window = true;

                homo_idx = (int)vlx.scf.homo_idx;
                lumo_idx = (int)vlx.scf.lumo_idx;

                vec4_t min_box = vec4_set1( FLT_MAX);
                vec4_t max_box = vec4_set1(-FLT_MAX);

                nto.atom_xyz = (vec3_t*)md_arena_allocator_push(arena, sizeof(vec3_t) * vlx.geom.num_atoms);
                nto.atom_r   = (float*) md_arena_allocator_push(arena, sizeof(float)  * vlx.geom.num_atoms);
                nto.num_atoms = vlx.geom.num_atoms;

                // Compute the PCA of the provided geometry
                // This is used in determining a better fitting volume for the orbitals
                vec4_t* xyzw = (vec4_t*)md_vm_arena_push(state.allocator.frame, sizeof(vec4_t) * vlx.geom.num_atoms);
                for (size_t i = 0; i < vlx.geom.num_atoms; ++i) {
                    nto.atom_xyz[i] = vec3_set((float)vlx.geom.coord_x[i], (float)vlx.geom.coord_y[i], (float)vlx.geom.coord_z[i]) * ANGSTROM_TO_BOHR;
                    nto.atom_r[i]   = md_util_element_vdw_radius(vlx.geom.atomic_number[i]) * ANGSTROM_TO_BOHR;
                    xyzw[i] = {(float)vlx.geom.coord_x[i], (float)vlx.geom.coord_y[i], (float)vlx.geom.coord_z[i], 1.0f};
                    min_box = vec4_min(min_box, xyzw[i]);
                    max_box = vec4_max(max_box, xyzw[i]);
                }
                min_aabb = vec3_from_vec4(min_box);
                max_aabb = vec3_from_vec4(max_box);

                md_molecule_t mol = {0};
                md_vlx_molecule_init(&mol, &vlx, state.allocator.frame);
                md_util_molecule_postprocess(&mol, state.allocator.frame, MD_UTIL_POSTPROCESS_ELEMENT_BIT | MD_UTIL_POSTPROCESS_RADIUS_BIT | MD_UTIL_POSTPROCESS_BOND_BIT);
                //gl_mol = md_gl_mol_create(&mol);

                uint32_t* colors = (uint32_t*)md_vm_arena_push(state.allocator.frame, mol.atom.count * sizeof(uint32_t));
                color_atoms_cpk(colors, mol.atom.count, mol);

                gl_rep = md_gl_rep_create(state.mold.gl_mol);
                md_gl_rep_set_color(gl_rep, 0, (uint32_t)mol.atom.count, colors, 0);

                com =  md_util_com_compute_vec4(xyzw, 0, vlx.geom.num_atoms, 0);
                mat3_t C = mat3_covariance_matrix_vec4(xyzw, 0, vlx.geom.num_atoms, com);
                mat3_eigen_t eigen = mat3_eigen(C);
                PCA = mat3_extract_rotation(eigen.vectors);

                // NTO
                if (vlx.rsp.num_excited_states > 0 && vlx.rsp.nto) {

                    nto.show_window = true;
                    camera_compute_optimal_view(&nto.target.pos, &nto.target.ori, &nto.target.dist, min_aabb, max_aabb, nto.distance_scale);
					nto.atom_group_idx = (uint8_t*)md_alloc(arena, sizeof(uint8_t) * mol.atom.count);
					MEMSET(nto.atom_group_idx, 0, sizeof(uint8_t) * mol.atom.count);
					
					for (int i = 1; i < (int)ARRAY_SIZE(nto.group.color); ++i) {
                        ImVec4 color = ImPlot::GetColormapColor(i - 1, ImPlotColormap_Deep);
                        nto.group.color[i] = vec_cast(color);
                        snprintf(nto.group.label[i], sizeof(nto.group.label[i]), "Group %i", i + 1);
                    }
                    snprintf(nto.group.label[0], sizeof(nto.group.label[0]), "Unassigned");
                    nto.group.color[0] = vec4_set(0.5f, 0.5f, 0.5f, 1.0f);

                    str_t file = {};
                    extract_file(&file, filename);
                    if (str_eq_cstr(file, "tq.out")) {
                        uint8_t index_from_text[23] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
                        //nto.atom_group_idx = index_from_text;
                        nto.group.count = 2;
                        MEMCPY(nto.atom_group_idx, index_from_text, sizeof(index_from_text));
                        snprintf(nto.group.label[0], sizeof(nto.group.label[0]), "Thio");
                        snprintf(nto.group.label[1], sizeof(nto.group.label[1]), "Quin");
                    }
                    else {

                        // @TODO: Remove once proper interface is there
					    nto.group.count = 2;
					    // Assign half of the atoms to group 1
                        MEMSET(nto.atom_group_idx, 1, sizeof(uint8_t) * mol.atom.count / 2);
                    }
					nto.gl_rep = md_gl_rep_create(state.mold.gl_mol);

                    for (size_t i = 0; i < mol.atom.count; ++i) {
                        const size_t group_idx = nto.atom_group_idx[i];
                        const uint32_t color = group_idx < nto.group.count ? convert_color(nto.group.color[group_idx]) : U32_MAGENTA;
                        colors[i] = color;
                    }

                    md_gl_rep_set_color(nto.gl_rep, 0, mol.atom.count, colors, 0);
                }

                // RSP
                rsp.show_window = true;
                rsp.hovered  = -1;
                rsp.selected = -1;

                // ORB
                orb.show_window = true;
                camera_compute_optimal_view(&orb.target.pos, &orb.target.ori, &orb.target.dist, min_aabb, max_aabb, orb.distance_scale);
                orb.mo_idx = homo_idx;
                orb.scroll_to_idx = homo_idx;

            } else {
                MD_LOG_INFO("Failed to load VeloxChem data");
                reset_data();
            }
        }
    }

    task_system::ID compute_nto_async(mat4_t* out_vol_mat, vec3_t* out_voxel_spacing, uint32_t* in_out_vol_tex, size_t nto_idx, size_t lambda_idx, md_vlx_nto_type_t type, md_gto_eval_mode_t mode, float samples_per_angstrom = DEFAULT_SAMPLES_PER_ANGSTROM) {
		md_allocator_i* alloc = md_get_heap_allocator();
        size_t num_pgtos = md_vlx_nto_pgto_count(&vlx);
        md_gto_t* pgtos  = (md_gto_t*)md_alloc(alloc, sizeof(md_gto_t) * num_pgtos);

        if (!md_vlx_nto_pgto_extract(pgtos, &vlx, nto_idx, lambda_idx, type)) {
            MD_LOG_ERROR("Failed to extract NTO pgtos for nto index: %zu and lambda: %zu", nto_idx, lambda_idx);
            md_free(alloc, pgtos, sizeof(md_gto_t) * num_pgtos);
            return task_system::INVALID_ID;
        }
        md_gto_cutoff_compute(pgtos, num_pgtos, 1.0e-6);

        OBB obb = compute_pgto_obb(PCA, pgtos, num_pgtos);

        vec3_t extent = obb.max_ext - obb.min_ext;

        // Target resolution per spatial unit (We want this number of samples per Ångström in each dimension)
        // Round up to some multiple of 8 -> 8x8x8 is the chunksize we process in parallel

        // Compute required volume dimensions
        int dim[3] = {
            CLAMP(ALIGN_TO((int)(extent.x * samples_per_angstrom * BOHR_TO_ANGSTROM), 8), 8, 512),
            CLAMP(ALIGN_TO((int)(extent.y * samples_per_angstrom * BOHR_TO_ANGSTROM), 8), 8, 512),
            CLAMP(ALIGN_TO((int)(extent.z * samples_per_angstrom * BOHR_TO_ANGSTROM), 8), 8, 512),
        };

        const vec3_t stepsize = vec3_div(extent, vec3_set((float)dim[0], (float)dim[1], (float)dim[2]));
        mat4_t vol_mat = compute_vol_mat(obb);

        vec3_t step_x = obb.basis.col[0] * stepsize.x;
        vec3_t step_y = obb.basis.col[1] * stepsize.y;
        vec3_t step_z = obb.basis.col[2] * stepsize.z;
        vec3_t origin = obb.basis * (obb.min_ext + stepsize * 0.5f);

        // Init and clear volume texture
        const float zero[4] = { 0,0,0,0 };
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, vol_fbo);
        gl::init_texture_3D(in_out_vol_tex, dim[0], dim[1], dim[2], GL_R16F);
        glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, *in_out_vol_tex, 0);
        glClearBufferfv(GL_COLOR, 0, zero);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

        // WRITE OUTPUT
        *out_vol_mat = vol_mat;
        *out_voxel_spacing = stepsize * BOHR_TO_ANGSTROM;

        struct Payload {
            AsyncGridEvalArgs args;
            md_allocator_i* alloc;
            size_t mem_size;
            void* mem;
            uint32_t* tex_ptr;
        };

        size_t mem_size = sizeof(Payload) + sizeof(float) * dim[0] * dim[1] * dim[2];
        void* mem = md_alloc(alloc, mem_size);
        Payload* payload = (Payload*)mem;
        *payload = {
            .args = {
                .grid = {
                    .data = (float*)((char*)mem + sizeof(Payload)),
                    .dim = {dim[0], dim[1], dim[2]},
                    .origin = {origin.x, origin.y, origin.z},
                    .step_x = {step_x.x, step_x.y, step_x.z},
                    .step_y = {step_y.x, step_y.y, step_y.z},
                    .step_z = {step_z.x, step_z.y, step_z.z},
                },
                .pgtos = pgtos,
                .num_pgtos = num_pgtos,
                .mode = mode,
            },
            .alloc = alloc,
            .mem_size = mem_size,
            .mem = mem,
            .tex_ptr = in_out_vol_tex,
        };

        task_system::ID async_task = evaluate_pgtos_on_grid_async(&payload->args);

        // Launch task for main (render) thread to update the volume texture
        task_system::ID main_task = task_system::create_main_task(STR_LIT("##Update Volume"), [](void* user_data) {
            Payload* data = (Payload*)user_data;

            // The init here is just to ensure that the volume has not changed its dimensions during the async evaluation
            gl::init_texture_3D(data->tex_ptr, data->args.grid.dim[0], data->args.grid.dim[1], data->args.grid.dim[2], GL_R16F);
            gl::set_texture_3D_data(*data->tex_ptr, data->args.grid.data, GL_R32F);

            md_free(data->alloc, data->args.pgtos, data->args.num_pgtos * sizeof(md_gto_t));
            md_free(data->alloc, data->mem, data->mem_size);
            }, payload);

        task_system::set_task_dependency(main_task, async_task);
        task_system::enqueue_task(async_task);

        return async_task;
    }

	struct AsyncGridEvalArgs {
		md_grid_t  grid;
		md_gto_t* pgtos;
		size_t num_pgtos;
		md_gto_eval_mode_t mode;
	};

    task_system::ID compute_mo_async(mat4_t* out_vol_mat, vec3_t* out_voxel_spacing, uint32_t* in_out_vol_tex, size_t mo_idx, md_gto_eval_mode_t mode, float samples_per_angstrom = DEFAULT_SAMPLES_PER_ANGSTROM) {
		md_allocator_i* alloc = md_get_heap_allocator();
        size_t num_pgtos = md_vlx_mol_pgto_count(&vlx);
        md_gto_t* pgtos  = (md_gto_t*)md_alloc(alloc, sizeof(md_gto_t)* num_pgtos);

        if (!md_vlx_mol_pgto_extract(pgtos, &vlx, mo_idx)) {
            MD_LOG_ERROR("Failed to extract molecular pgtos for orbital index: %zu", mo_idx);
            md_free(md_get_heap_allocator(), pgtos, sizeof(md_gto_t) * num_pgtos);
            return task_system::INVALID_ID;
        }
        md_gto_cutoff_compute(pgtos, num_pgtos, 1.0e-6);

        OBB obb = compute_pgto_obb(PCA, pgtos, num_pgtos);
        vec3_t extent = obb.max_ext - obb.min_ext;

        // Target resolution per spatial unit (We want this number of samples per Ångström in each dimension)
        // Round up to some multiple of 8 -> 8x8x8 is the chunksize we process in parallel

        // Compute required volume dimensions
        int dim[3] = {
            CLAMP(ALIGN_TO((int)(extent.x * samples_per_angstrom * BOHR_TO_ANGSTROM), 8), 8, 512),
            CLAMP(ALIGN_TO((int)(extent.y * samples_per_angstrom * BOHR_TO_ANGSTROM), 8), 8, 512),
            CLAMP(ALIGN_TO((int)(extent.z * samples_per_angstrom * BOHR_TO_ANGSTROM), 8), 8, 512),
        };

        const vec3_t stepsize = vec3_div(extent, vec3_set((float)dim[0], (float)dim[1], (float)dim[2]));
        mat4_t vol_mat = compute_vol_mat(obb);

        vec3_t step_x = obb.basis.col[0] * stepsize.x;
        vec3_t step_y = obb.basis.col[1] * stepsize.y;
        vec3_t step_z = obb.basis.col[2] * stepsize.z;
        vec3_t origin = obb.basis * (obb.min_ext + stepsize * 0.5f);

        // Init and clear volume texture
        const float zero[4] = {0,0,0,0};
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, vol_fbo);
        gl::init_texture_3D(in_out_vol_tex, dim[0], dim[1], dim[2], GL_R16F);
        glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, *in_out_vol_tex, 0);
        glClearBufferfv(GL_COLOR, 0, zero);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

        // WRITE OUTPUT
        *out_vol_mat = vol_mat;
        *out_voxel_spacing = stepsize * BOHR_TO_ANGSTROM;

        struct Payload {
            AsyncGridEvalArgs args;
            uint32_t* tex_ptr;
            md_allocator_i* alloc;
            size_t alloc_size;
        };

		size_t mem_size = sizeof(Payload) + sizeof(float) * dim[0] * dim[1] * dim[2];
		void*       mem = md_alloc(alloc, mem_size);
		Payload* payload = (Payload*)mem;
		*payload = {
			.args = {
				.grid = {
					.data = (float*)((char*)mem + sizeof(Payload)),
					.dim = {dim[0], dim[1], dim[2]},
					.origin = {origin.x, origin.y, origin.z},
					.step_x = {step_x.x, step_x.y, step_x.z},
                    .step_y = {step_y.x, step_y.y, step_y.z},
                    .step_z = {step_z.x, step_z.y, step_z.z},
				},
				.pgtos = pgtos,
				.num_pgtos = num_pgtos,
				.mode = mode,
			},
			.tex_ptr = in_out_vol_tex,
			.alloc = alloc,
			.alloc_size = mem_size,
		};

		task_system::ID async_task = evaluate_pgtos_on_grid_async(&payload->args);

        // Launch task for main (render) thread to update the volume texture
        task_system::ID main_task = task_system::create_main_task(STR_LIT("##Update Volume"), [](void* user_data) {
            Payload* data = (Payload*)user_data;

            // The init here is just to ensure that the volume has not changed its dimensions during the async evaluation
            gl::init_texture_3D(data->tex_ptr, data->args.grid.dim[0], data->args.grid.dim[1], data->args.grid.dim[2], GL_R16F);
            gl::set_texture_3D_data(*data->tex_ptr, data->args.grid.data, GL_R32F);

            md_free(data->alloc, data->args.pgtos, data->args.num_pgtos * sizeof(md_gto_t));
            md_free(data->alloc, data, data->alloc_size);
            }, payload);

        task_system::set_task_dependency(main_task, async_task);
        task_system::enqueue_task(async_task);

        return async_task;
    }

    static inline ImVec4 make_highlight_color(const ImVec4& color, float factor = 1.3f) {
        // Ensure the factor is not too high, to avoid over-brightening
        factor = (factor < 1.0f) ? 1.0f : factor;

        // Increase the brightness of each component by the factor
        float r = color.x * factor;
        float g = color.y * factor;
        float b = color.z * factor;

        // Clamp the values to the range [0, 1]
        r = (r > 1.0f) ? 1.0f : r;
        g = (g > 1.0f) ? 1.0f : g;
        b = (b > 1.0f) ? 1.0f : b;

        // Return the new highlighted color with the same alpha
        return ImVec4(r, g, b, color.w);
    }

    static inline void draw_vertical_sankey_flow(ImDrawList* draw_list, ImVec2 source_pos, ImVec2 dest_pos, float thickness, ImU32 flow_color) {
        // Get the draw list from the current window

        // Define the control points for the Bezier curve
        ImVec2 p1 = ImVec2(source_pos.x + thickness / 2, source_pos.y);  // Start point (bottom-center of source node)
        ImVec2 p4 = ImVec2(dest_pos.x + thickness / 2, dest_pos.y);  // End point (top-center of destination node)

        float dist = fabs(p1.y - p4.y);
        float curve_offset = dist / 4;

        ImVec2 p2 = ImVec2(source_pos.x + thickness / 2, source_pos.y - curve_offset);  // Control point 1
        ImVec2 p3 = ImVec2(dest_pos.x + thickness / 2, dest_pos.y + curve_offset);  // Control point 2


        // Define the color and thickness for the flow
        //ImU32 flow_color = IM_COL32(100, 149, 237, 255);  // Cornflower Blue
        // ImBezierCubicCalc Use this for calculating mouse distance to curve
        // Draw the Bezier curve representing the flow
        draw_list->AddBezierCubic(p1, p2, p3, p4, flow_color, thickness, 100);
    }

    static inline void draw_aligned_text(ImDrawList* draw_list, const char* text, ImVec2 pos, ImVec2 alignment = { 0,0 }) {
        ImVec2 text_size = ImGui::CalcTextSize(text);
        ImVec2 text_pos = pos - text_size * alignment;
        draw_list->AddText(text_pos, ImGui::ColorConvertFloat4ToU32({ 0,0,0,1 }), text);
    }

    static inline void im_sankey_diagram(ImRect area, Nto* nto) {
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
        const float plot_percent = 0.8;
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
        const float bar_height = plot_area.GetHeight() * 0.05;
        int num_bars = nto->group.count;
        int num_gaps = num_bars - 1;
        float gap_size = plot_area.GetWidth() * 0.05;
        float bars_avail_width = plot_area.GetWidth() - gap_size * num_gaps;

        //Calculate bar percentages
        float start_sum = 0;
        float end_sum = 0;
        for (size_t i = 0; i < num_bars; i++) {
            start_sum += nto->transition_density_hole[i];
            end_sum += nto->transition_density_part[i];
        }
        md_array(float) start_percentages = md_array_create(float, num_bars, temp_alloc);
        md_array(float) end_percentages   = md_array_create(float, num_bars, temp_alloc);
        for (size_t i = 0; i < num_bars; i++) {
            start_percentages[i] = nto->transition_density_hole[i] / start_sum;
            end_percentages[i] = nto->transition_density_part[i] / end_sum;
        }

        //Calculate start positions
        md_array(float) start_positions = md_array_create(float, num_bars, temp_alloc);
        float cur_bottom_pos = plot_area.Min.x;
        for (int i = 0; i < num_bars; i++) {
            start_positions[i] = cur_bottom_pos;
            cur_bottom_pos += bars_avail_width * start_percentages[i] + gap_size;
        }

        //Calculate end positions
        md_array(float) end_positions = md_array_create(float, num_bars, temp_alloc);
        md_array(float) sub_end_positions = md_array_create(float, num_bars, temp_alloc);
        float cur_pos = plot_area.Min.x;
        for (int end_i = 0; end_i < num_bars; end_i++) {
            end_positions[end_i] = cur_pos;
            sub_end_positions[end_i] = cur_pos;
            cur_pos += bars_avail_width * end_percentages[end_i] + gap_size;
        }


        //Draw curves
        for (int start_i = 0; start_i < num_bars; start_i++) {
            ImVec2 start_pos = { start_positions[start_i], plot_area.Max.y - bar_height + 0.1f * bar_height };

            ImVec4 start_col = vec_cast(nto->group.color[start_i]);

            start_col.w = 0.5;
            for (int end_i = 0; end_i < num_bars; end_i++) {
                float percentage = nto->transition_matrix[end_i * num_bars + start_i];
                ImVec4 end_col = vec_cast(nto->group.color[end_i]);

                if (percentage != 0) {
                    float width = bars_avail_width * percentage;
                    ImVec2 end_pos = { sub_end_positions[end_i], plot_area.Min.y + bar_height - 0.1f * bar_height };

                    int vert_beg = draw_list->VtxBuffer.Size;
                    draw_vertical_sankey_flow(draw_list, start_pos, end_pos, width, ImGui::ColorConvertFloat4ToU32(start_col));
                    int vert_end = draw_list->VtxBuffer.Size;

                    // Apply a gradient based on y-value from start color to end color if they belong to different groups
                    if (end_i != start_i) {
                        ImVec2 grad_p0 = {start_pos.x, start_pos.y};
                        ImVec2 grad_p1 = {start_pos.x, end_pos.y};
                        ImGui::ShadeVertsLinearColorGradientKeepAlpha(draw_list, vert_beg, vert_end, grad_p0, grad_p1, ImGui::ColorConvertFloat4ToU32(start_col), ImGui::ColorConvertFloat4ToU32(end_col));
                    }
                    
                    ImVec2 midpoint = (start_pos + end_pos) * 0.5 + ImVec2{width / 2, 0};
                    char lable[16];
                    sprintf(lable, "%3.2f%%", percentage * 100);
                    if (width > ImGui::CalcTextSize(lable).x) {
                        draw_aligned_text(draw_list, lable, midpoint, { 0.5, 0.5 });
                    }

                    start_pos.x += width;
                    sub_end_positions[end_i] += width;
                }
            }
        }

        //Draw bars
        for (int i = 0; i < num_bars; i++) {
            ImVec4 bar_color = vec_cast(nto->group.color[i]);
            ImVec2 mouse_pos = ImGui::GetMousePos();

            //Calculate start
            ImVec2 start_p0 = { start_positions[i], plot_area.Max.y - bar_height};
            ImVec2 start_p1 = { start_positions[i] + bars_avail_width * start_percentages[i], plot_area.Max.y };
            ImVec2 start_midpoint = { (start_p0.x + start_p1.x) * 0.5f, start_p1.y };
            ImRect start_bar = ImRect{ start_p0, start_p1 };

            //Calculate end
            ImVec2 end_p0 = { end_positions[i], plot_area.Min.y };
            ImVec2 end_p1 = { end_positions[i] + bars_avail_width * end_percentages[i], plot_area.Min.y + bar_height };
            ImRect end_bar = ImRect{ end_p0, end_p1 };
            ImVec2 end_midpoint = { (end_p0.x + end_p1.x) * 0.5f, end_p0.y };

            if (start_bar.Contains(mouse_pos) || end_bar.Contains(mouse_pos)) {
                bar_color = make_highlight_color(bar_color);
            }

            //Draw start
            draw_list->AddRectFilled(start_p0, start_p1, ImGui::ColorConvertFloat4ToU32(bar_color));
            draw_list->AddRect(start_p0, start_p1, ImGui::ColorConvertFloat4ToU32({0,0,0,0.5}));
            char start_lable[16];
            sprintf(start_lable, "%3.2f%%", start_percentages[i] * 100);
            draw_aligned_text(draw_list, nto->group.label[i], start_midpoint, {0.5, -0.2});
            draw_aligned_text(draw_list, start_lable, start_midpoint, { 0.5, -1.2 });

            //Draw end
            draw_list->AddRectFilled(end_p0, end_p1, ImGui::ColorConvertFloat4ToU32(bar_color));
            draw_list->AddRect(end_p0, end_p1, ImGui::ColorConvertFloat4ToU32({ 0,0,0,0.5 }));
            char end_lable[16];
            sprintf(end_lable, "%3.2f%%", end_percentages[i] * 100);
            draw_aligned_text(draw_list, nto->group.label[i], end_midpoint, {0.5, 1.2});
            draw_aligned_text(draw_list, end_lable, end_midpoint, { 0.5, 2.2 });
        }
    }
    
    bool compute_nto_group_values_async(task_system::ID* out_eval_task, task_system::ID* out_segment_task, float* out_group_values, size_t num_groups, const uint8_t* point_group_idx, const vec3_t* point_xyz, const float* point_r, size_t num_points, size_t nto_idx, size_t lambda_idx, md_vlx_nto_type_t type, md_gto_eval_mode_t mode, float samples_per_angstrom = DEFAULT_SAMPLES_PER_ANGSTROM) {       
        md_allocator_i* alloc = md_get_heap_allocator();
        size_t num_pgtos = md_vlx_nto_pgto_count(&vlx);
        md_gto_t* pgtos = (md_gto_t*)md_alloc(alloc, sizeof(md_gto_t) * num_pgtos);

        if (!md_vlx_nto_pgto_extract(pgtos, &vlx, nto_idx, lambda_idx, type)) {
            MD_LOG_ERROR("Failed to extract NTO pgtos for nto index: %zu and lambda: %zu", nto_idx, lambda_idx);
            md_free(alloc, pgtos, sizeof(md_gto_t) * num_pgtos);
            return false;
        }

        md_gto_cutoff_compute(pgtos, num_pgtos, 1.0e-6);
        OBB obb = compute_pgto_obb(PCA, pgtos, num_pgtos);
        vec3_t extent = obb.max_ext - obb.min_ext;

        // Target resolution per spatial unit (We want this number of samples per Ångström in each dimension)
        // Round up to some multiple of 8 -> 8x8x8 is the chunksize we process in parallel

        // Compute required volume dimensions
        int dim[3] = {
            CLAMP(ALIGN_TO((int)(extent.x * samples_per_angstrom * BOHR_TO_ANGSTROM), 8), 8, 512),
            CLAMP(ALIGN_TO((int)(extent.y * samples_per_angstrom * BOHR_TO_ANGSTROM), 8), 8, 512),
            CLAMP(ALIGN_TO((int)(extent.z * samples_per_angstrom * BOHR_TO_ANGSTROM), 8), 8, 512),
        };

        const vec3_t stepsize = vec3_div(extent, vec3_set((float)dim[0], (float)dim[1], (float)dim[2]));
        mat4_t vol_mat = compute_vol_mat(obb);

        vec3_t step_x = obb.basis.col[0] * stepsize.x;
        vec3_t step_y = obb.basis.col[1] * stepsize.y;
        vec3_t step_z = obb.basis.col[2] * stepsize.z;
        vec3_t origin = obb.basis * (obb.min_ext + stepsize * 0.5f);

        struct Payload {
            AsyncGridEvalArgs args;

            float* dst_group_values;
			size_t num_groups;

            const uint8_t* point_group_idx;
            const vec3_t* point_xyz;
            const float* point_r;
            size_t num_points;

            md_allocator_i* alloc;
            size_t alloc_size;
        };

        size_t mem_size = sizeof(Payload) + sizeof(float) * dim[0] * dim[1] * dim[2];
        void* mem = md_alloc(alloc, mem_size);
        Payload* payload = (Payload*)mem;
        *payload = {
            .args = {
                .grid = {
                    .data = (float*)((char*)mem + sizeof(Payload)),
                    .dim = {dim[0], dim[1], dim[2]},
                    .origin = {origin.x, origin.y, origin.z},
                    .step_x = {step_x.x, step_x.y, step_x.z},
                    .step_y = {step_y.x, step_y.y, step_y.z},
                    .step_z = {step_z.x, step_z.y, step_z.z},
                },
                .pgtos = pgtos,
                .num_pgtos = num_pgtos,
                .mode = mode,
            },
			.dst_group_values = out_group_values,
            .num_groups = num_groups,
            
			.point_group_idx = point_group_idx,
			.point_xyz = point_xyz,
            .point_r = point_r,
			.num_points = num_points,
            .alloc = alloc,
            .alloc_size = mem_size,
        };

        task_system::ID eval_task = evaluate_pgtos_on_grid_async(&payload->args);

        // @TODO: This should be performed as a range task in parallel
        task_system::ID segment_task = task_system::create_pool_task(STR_LIT("##Segment Volume"), [](void* user_data) {
            Payload* data = (Payload*)user_data;

#if DEBUG
            double sum = 0.0;
            size_t len = data->args.grid.dim[0] * data->args.grid.dim[1] * data->args.grid.dim[2];
            for (size_t i = 0; i < len; ++i) {
                sum += data->args.grid.data[i];
            }
            MD_LOG_DEBUG("SUM: %g");
#endif

            MD_LOG_DEBUG("Starting segmentation of volume");
            grid_segment_and_attribute(data->dst_group_values, data->num_groups, data->point_group_idx, data->point_xyz, data->point_r, data->num_points, &data->args.grid);
            MD_LOG_DEBUG("Finished segmentation of volume");

            md_free(data->alloc, data->args.pgtos, data->args.num_pgtos * sizeof(md_gto_t));
            md_free(data->alloc, data, data->alloc_size);
        }, payload);

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
    task_system::ID evaluate_pgtos_on_grid_async(AsyncGridEvalArgs* args) {
        ASSERT(args);

        // We evaluate the in parallel over smaller NxNxN blocks
        const uint32_t num_blocks = (args->grid.dim[0] / BLK_DIM) * (args->grid.dim[1] / BLK_DIM) * (args->grid.dim[2] / BLK_DIM);
        task_system::ID async_task = task_system::create_pool_task(STR_LIT("Evaluate Orbital"), 0, num_blocks, [](uint32_t range_beg, uint32_t range_end, void* user_data, uint32_t thread_num) {
            (void)thread_num;
            AsyncGridEvalArgs* data = (AsyncGridEvalArgs*)user_data;
            MD_LOG_DEBUG("Starting async eval of orbital grid [%i][%i][%i]", data->grid.dim[0], data->grid.dim[1], data->grid.dim[2]);

            // Number of NxNxN blocks in each dimension
            int num_blk[3] = {
                data->grid.dim[0] / BLK_DIM,
                data->grid.dim[1] / BLK_DIM,
                data->grid.dim[2] / BLK_DIM,
            };

			md_grid_t* grid = &data->grid;

            size_t temp_pos = md_temp_get_pos();
            md_gto_t* sub_pgtos = (md_gto_t*)md_temp_push(sizeof(md_gto_t) * data->num_pgtos);

            for (uint32_t i = range_beg; i < range_end; ++i) {
                // Determine block index from linear input index i
                int blk_x = i % num_blk[0];
                int blk_y = (i / num_blk[0]) % num_blk[1];
                int blk_z = i / (num_blk[0] * num_blk[1]);

                int off_idx[3] = { blk_x * BLK_DIM, blk_y * BLK_DIM, blk_z * BLK_DIM };
                int len_idx[3] = { BLK_DIM, BLK_DIM, BLK_DIM };

                int beg_idx[3] = {off_idx[0], off_idx[1], off_idx[2]};
                int end_idx[3] = {off_idx[0] + len_idx[0], off_idx[1] + len_idx[1], off_idx[2] + len_idx[2]};

                // @TODO: This filtering of PGTOs can be improved by transforming the PGTO x,y,z,radius into the OBB that defines the region

                float aabb_min[3] = {
                    grid->origin[0] + beg_idx[0] * grid->step_x[0] + beg_idx[1] * grid->step_y[0] + beg_idx[2] * grid->step_z[0],
                    grid->origin[1] + beg_idx[0] * grid->step_x[1] + beg_idx[1] * grid->step_y[1] + beg_idx[2] * grid->step_z[1],
                    grid->origin[2] + beg_idx[0] * grid->step_x[2] + beg_idx[1] * grid->step_y[2] + beg_idx[2] * grid->step_z[2],
                };
                float aabb_max[3] = {
                    grid->origin[0] + end_idx[0] * grid->step_x[0] + end_idx[1] * grid->step_y[0] + end_idx[2] * grid->step_z[0],
                    grid->origin[1] + end_idx[0] * grid->step_x[1] + end_idx[1] * grid->step_y[1] + end_idx[2] * grid->step_z[1],
                    grid->origin[2] + end_idx[0] * grid->step_x[2] + end_idx[1] * grid->step_y[2] + end_idx[2] * grid->step_z[2],
                };

                size_t num_sub_pgtos = md_gto_aabb_test(sub_pgtos, aabb_min, aabb_max, data->pgtos, data->num_pgtos);
                md_gto_grid_evaluate_sub(grid, off_idx, len_idx, sub_pgtos, num_sub_pgtos, data->mode);
            }

            md_temp_set_pos_back(temp_pos);
        }, args);

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


    enum x_unit_t {
        X_UNIT_EV,
        X_UNIT_NM,
        X_UNIT_CM_INVERSE,
        X_UNIT_HARTREE,
    };

    enum broadening_mode_t {
        BROADENING_GAUSSIAN,
        BROADENING_LORENTZIAN,
    };

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

    static inline double lorentzian(double x, double x_0, double gamma) {
        double sigma = gamma / 2;
        double res = (1 / PI) * sigma / (pow((x - x_0), 2) + pow(sigma, 2));
        return res;
    }

    static inline double phys_lorentzian(double x, double x_0, double gamma, double intensity) {
        double sigma = gamma / 2;
        double res = intensity * pow(sigma, 2) / (pow((x - x_0), 2) + pow(sigma, 2));
        return res;
    }

    static inline double gaussian(double x, double x_0, double gamma) {
        double sigma = gamma / 2.3548;
        return (1 / (sigma * sqrt(2 * PI))) * exp(-(pow(x - x_0, 2) / (2 * pow(sigma, 2)))); 
    }

    //TODO: Check that phys_gaussian implementation is actually correct
    static inline double phys_gaussian(double p, double p_0, double gamma, double intensity) {
        double sigma = gamma / 2;
        double x = (p - p_0) / sigma;
        return intensity * exp(-log(2) * pow(x,2));
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

    static inline void general_broadening(double* y_out, const double* x, size_t num_samples, const double* y_peaks, const double* x_peaks, size_t num_peaks, double (*distr_func)(double x, double x_0, double gamma, double intensity), double gamma) {
        double integral = 0;
        double dist = x[1] - x[0];
        for (size_t si = 0; si < num_samples; si++) {
            double sum = 0;
            double b = 0;
            for (size_t pi = 0; pi < num_peaks; pi++) {
                b = (*distr_func)(x[si], x_peaks[pi], gamma, y_peaks[pi]);
                sum += b;
            }
            y_out[si] = sum;
            integral += y_out[si] * dist;
        }

        double i_sum = integral;
    }

    static inline void osc_to_eps(double* eps_out, const double* x, size_t num_samples, const double* osc_peaks, const double* x_peaks, size_t num_peaks, double (*distr_func)(double x, double x_0, double gamma), double gamma) {
        double c = 137.035999;
        double a_0 = 5.29177210903e-11;
        double NA = 6.02214076e23;
        double eV2au = 1 / 27.211396;

        for (size_t si = 0; si < num_samples; si++) {
            double sum = 0;
            for (size_t pi = 0; pi < num_peaks; pi++) {
                sum += (*distr_func)(x[si] * eV2au, x_peaks[pi] * eV2au, gamma * eV2au) * (osc_peaks[pi] / (x_peaks[pi] * eV2au));
            }
            double sigma = 2 * pow(PI, 2) * x[si] * eV2au * sum / c;
            double sigma_cm2 = sigma * pow(a_0, 2) * 1e4;
            eps_out[si] = sigma_cm2 * NA / (log(10) * 1e3);
        }
    }

    static inline void rot_to_eps_delta(double* eps_out, const double* x, size_t num_samples, const double* rot_peaks, const double* x_peaks, size_t num_peaks, double (*distr_func)(double x, double x_0, double gamma), double gamma) {
        double inv = 1 / (22.94);
        double eV2au = 1 / 27.211396;

        for (size_t si = 0; si < num_samples; si++) {
            double sum = 0;
            for (size_t pi = 0; pi < num_peaks; pi++) {
                sum += (*distr_func)(x[si], x_peaks[pi], gamma) * rot_peaks[pi] * x_peaks[pi];
            }
            
            eps_out[si] = sum * inv;
        }
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

        for (int i = 0; i < num_peaks; i++) {
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
                closest_idx = i;
            }
            else if (sqrt(pow(distance_x, 2) + pow(distance_y, 2)) < closest_distance) {
                closest_distance = sqrt(pow(distance_x, 2) + pow(distance_y, 2));
                closest_idx = i;
                //ImPlot::Annotation()
            }
            else if (distance_x > closest_distance){
                // We are now so far away that a closer bar will not occur, no matter the y value.
                break;
            }
        }
        return closest_distance < proxy_distance ? closest_idx : -1;
    }

    static inline void draw_bar(int id, double x, double y, double width, ImVec4 color) {
        double x1 = x - width / 2;
        double x2 = x + width / 2;
        double y1 = 0;
        ImPlot::DragRect(id, &x1, &y1, &x2, &y, color, ImPlotDragToolFlags_NoInputs);
    }

    //Calculates the maximum point and populates out_point with it
    static inline void max_points(double* out_points, const double* in_peaks, size_t num_peaks, double offset = 0.05) {
        double y_max = in_peaks[0];
        for (size_t i = 0; i < num_peaks; i++) {
            y_max = MAX(y_max, in_peaks[i]);
        }
        for (size_t i = 0; i < num_peaks; i++) {
            out_points[i] = y_max + y_max * offset;
        }
    }

    static inline bool is_all_zero(const double* array, size_t count) {
        for (size_t i = 0; i < count; i++) {
            if (array[i] != 0.0) {
                return false;
            }
        }
        return true;
    }

    /*
    static inline void osc_to_eps(double* eps_out, const double* x_peaks, const double* osc_peaks, size_t num_peaks) {
        double NA = 6.02214076e23;
        double ln10 = 2.30258509299;
        //double eps_0 = 8.854188e-12;
        double c = 137.035999; //in au
        double a_0 = 5.29177210903e-11;
        // me = 1
        double sigma = 0;
        double sigma_cm2 = 0;
        for (size_t i = 0; i < num_peaks; i++) {
            sigma = ((2 * pow(PI, 2) * x_peaks[i]) / c) * osc_peaks[i];
            sigma_cm2 = sigma * pow(a_0, 2) * 1e4;
            eps_out[i] = (sigma_cm2 * NA) / (ln10 * 1000);
        }
        //We do the broadening in the next, external step
    }
    */

    void draw_summary_window(ApplicationState& state) {
        if (!scf.show_window) { return; }
        if (vlx.scf.iter.count == 0) { return; }

        size_t temp_pos = md_temp_get_pos();
        defer {  md_temp_set_pos_back(temp_pos); };

        // We set up iterations as doubles for easier use
        double* iter = (double*)md_temp_push(sizeof(double) * vlx.scf.iter.count);
        for (int i = 0; i < vlx.scf.iter.count; ++i) {
            iter[i] = (double)vlx.scf.iter.iteration[i];
        }
        static ImPlotRect lims{ 0,1,0,1 };
        // The actual plot


        double* energy_offsets = (double*)md_temp_push(sizeof(double) * (int)vlx.scf.iter.count);
        double ref_energy = vlx.scf.iter.energy_total[vlx.scf.iter.count - 1];
        for (size_t i = 0; i < vlx.scf.iter.count; i++) {
            energy_offsets[i] = fabs(vlx.scf.iter.energy_total[i] - ref_energy);
        }
        double y1_to_y2_mult = axis_conversion_multiplier(vlx.scf.iter.gradient_norm, energy_offsets, vlx.scf.iter.count, vlx.scf.iter.count);

        ImGui::SetNextWindowSize({ 300, 350 }, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Summary", &scf.show_window)) {
            if (ImGui::TreeNode("Level of calculation")) {
                ImGui::Text("Method:");
                ImGui::Text("Basis Set: %s", (const char*)vlx.basis.ident.ptr);
                ImGui::Spacing();
                
                ImGui::TreePop();
            }
            if (ImGui::TreeNode("System Information")) {
                ImGui::Text("Num Atoms:           %6zu", vlx.geom.num_atoms);
                ImGui::Text("Num Alpha Electrons: %6zu", vlx.geom.num_alpha_electrons);
                ImGui::Text("Num Beta Electrons:  %6zu", vlx.geom.num_beta_electrons);
                ImGui::Text("Molecular Charge:    %6i", vlx.geom.molecular_charge);
                ImGui::Text("Spin Multiplicity:   %6i", vlx.geom.spin_multiplicity);
                ImGui::Spacing();
                ImGui::TreePop();
            }
            if (ImGui::TreeNode("SCF")) {
                if (ImPlot::BeginPlot("SCF")) {
                    ImPlot::SetupAxisLimits(ImAxis_X1, 1.0, (int)vlx.scf.iter.count);
                    ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);
                    ImPlot::SetupAxes("Iteration", "Gradient Norm (au)");
                    // We draw 2 y axis as "Energy total" has values in a different range then the rest of the data
                    ImPlot::SetupAxis(ImAxis_Y2, "Energy (hartree)", ImPlotAxisFlags_AuxDefault);
                    ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

                    ImPlot::PlotLine("Gradient", iter, vlx.scf.iter.gradient_norm, (int)vlx.scf.iter.count);
                    ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
                    ImPlot::PlotLine("Energy", iter, energy_offsets, (int)vlx.scf.iter.count - 1);
                    lims = ImPlot::GetPlotLimits(ImAxis_X1, ImAxis_Y1);
                    //ImPlot::PlotLine("Density Change", iter, vlx.scf.iter.density_change, (int)vlx.scf.iter.count);
                    //ImPlot::PlotLine("Energy Change", iter, vlx.scf.iter.energy_change, (int)vlx.scf.iter.count);
                    //ImPlot::PlotLine("Max Gradient", iter, vlx.scf.iter.max_gradient, (int)vlx.scf.iter.count);
                    ImPlot::EndPlot();
                }
                ImGui::Spacing();
                ImGui::Text("Total energy:              %16.10f a.u.", vlx.scf.total_energy);
                ImGui::Text("Electronic energy:         %16.10f a.u.", vlx.scf.electronic_energy);
                ImGui::Text("Nuclear repulsion energy:  %16.10f a.u.", vlx.scf.nuclear_repulsion_energy);
                ImGui::Text("Gradient norm:             %16.10f a.u.", vlx.scf.gradient_norm);
                ImGui::Spacing();
                ImGui::TreePop();
            }

            if (ImGui::TreeNode("Geometry")) {
                if (vlx.geom.num_atoms) {
                    static const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollX |
                                                         ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit;

                    static const ImGuiTableColumnFlags columns_base_flags = ImGuiTableColumnFlags_NoSort;

                    if (ImGui::BeginTable("Geometry Table", 5, flags, ImVec2(500, -1), 0)) {
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
                        for (int row_n = 0; row_n < vlx.geom.num_atoms; row_n++) {

                            ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                            bool is_sel = md_bitfield_test_bit(&state.selection.selection_mask, row_n); //If atom is selected, mark it as such
                            bool is_hov = md_bitfield_test_bit(&state.selection.highlight_mask, row_n); //If atom is hovered, mark it as such
                            bool hov_col = false;
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
                                if (state.mold.mol.atom.count > row_n) {
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
                            ImGui::Text(vlx.geom.atom_symbol[row_n].buf);
                            ImGui::TableNextColumn();
                            ImGui::Text("%12.6f", vlx.geom.coord_x[row_n]);
                            ImGui::TableNextColumn();
                            ImGui::Text("%12.6f", vlx.geom.coord_y[row_n]);
                            ImGui::TableNextColumn();
                            ImGui::Text("%12.6f", vlx.geom.coord_z[row_n]);

                            ImGui::PopStyleColor(1);
                                
                        }
                        if (!item_hovered && ImGui::IsWindowHovered()) {
                            //Makes sure that we clear the highlight if we are in this window, but don't hover an item
                            md_bitfield_clear(&state.selection.highlight_mask);
                        }

                        ImGui::PopStyleColor(2);
                        ImGui::EndTable();
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
            {rsp.x_unit_samples, rsp.eps, str_from_cstr("Absorption"), str_from_cstr((const char*)u8"ε (L mol⁻¹ cm⁻¹)")},
            {rsp.x_unit_samples, rsp.ecd, str_from_cstr("ECD"), str_from_cstr((const char*)u8"Δε(ω) (L mol⁻¹ cm⁻¹)")},
            {rsp.vib_x, rsp.vib_y, str_from_cstr("Vibration"), str_from_cstr("IR Intensity (km/mol)")}
        };

        int num_properties = ARRAY_SIZE(properties);
        
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
                for (int i = 0; i < num_properties; ++i) {
                    if (ImGui::Selectable(properties[i].lable.ptr, property_idx == i)) {
                        property_idx = i;
                    }
                }
                ImGui::EndCombo();
            }

            if (property_idx == 0 || property_idx == 1) {
                x_unit = rsp.x_unit;
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

    void draw_rsp_window(ApplicationState& state) {
        if (!rsp.show_window) return;
        if (vlx.rsp.num_excited_states == 0) return;
        // Keep track of the temp position so we can reset to it after we are done
        size_t temp_pos = md_temp_get_pos();
        defer { md_temp_set_pos_back(temp_pos); };

        static ImVec2 mouse_pos = { 0,0 };

        const char* broadening_str[] = { "Gaussian","Lorentzian" };

        const int num_samples = 1024;

        ImGui::SetNextWindowSize({ 300, 350 }, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Spectra", &rsp.show_window, ImGuiWindowFlags_MenuBar)) {
            if (ImGui::BeginMenuBar())
            {
                if (ImGui::BeginMenu("File")) {
                    char path_buf[1024] = "";
                    if (ImGui::MenuItem("Export")) {
                        rsp.show_export_window = true;
                    }
                    ImGui::EndMenu();
                }
                ImGui::EndMenuBar();
            }
            
            if (rsp.first_plot_rot_ecd) { ImGui::SetNextItemOpen(true); }
            if (ImGui::TreeNode("Absorption & ECD")) {
                bool refit1 = false;
                bool recalculate1 = false;

                static float gamma1 = 0.123;
                static x_unit_t x_unit = X_UNIT_EV;
                static broadening_mode_t broadening_mode1 = BROADENING_LORENTZIAN;
                const char* x_unit_str[] = {"Energy (eV)", "Wavelength (nm)", (const char*)u8"Wavenumber (cm⁻¹)", "Energy (hartree)"};

                recalculate1 = ImGui::SliderFloat((const char*)u8"Broadening γ HWHM (eV)", &gamma1, 0.01f, 1.0f);
                refit1 |= ImGui::Combo("Broadening mode", (int*)(&broadening_mode1), broadening_str, IM_ARRAYSIZE(broadening_str));
                refit1 |= ImGui::Combo("X unit", (int*)(&x_unit), x_unit_str, IM_ARRAYSIZE(x_unit_str));
                rsp.x_unit = x_unit_str[x_unit];

                const int num_peaks = (int)vlx.rsp.num_excited_states;
                const double* y_osc_peaks = vlx.rsp.absorption_osc_str;
                const double* y_cgs_peaks = vlx.rsp.electronic_circular_dichroism_cgs;

                if (rsp.first_plot_rot_ecd) {
                    rsp.x_ev_samples = md_array_create(double, num_samples, arena);
                    rsp.x_unit_samples = md_array_create(double, num_samples, arena);
                    rsp.x_unit_peaks = md_array_create(double, num_peaks, arena);
                    rsp.eps = md_array_create(double, num_samples, arena);
                    rsp.ecd = md_array_create(double, num_samples, arena);

                    // Populate x_values
                    const double x_min = vlx.rsp.absorption_ev[0] - 1.0;
                    const double x_max = vlx.rsp.absorption_ev[num_peaks - 1] + 1.0;
                    for (int i = 0; i < num_samples; ++i) {
                        double t = (double)i / (double)(num_samples - 1);
                        double value = lerp(x_min, x_max, t);
                        rsp.x_ev_samples[i] = value;
                    }
                }

                // double* temp_x_values  = (double*)md_temp_push(sizeof(double) * num_samples);
                // double* y_ecd_str = (double*)md_temp_push(sizeof(double) * num_samples);
                // double* y_eps_str   = (double*)md_temp_push(sizeof(double) * num_samples);

                ImVec2* pixel_osc_peaks = (ImVec2*)md_temp_push(sizeof(ImVec2) * num_peaks);
                ImVec2* pixel_cgs_peaks = (ImVec2*)md_temp_push(sizeof(ImVec2) * num_peaks); 
                ImVec2* pixel_osc_points = (ImVec2*)md_temp_push(sizeof(ImVec2) * num_peaks);
                ImVec2* pixel_cgs_points = (ImVec2*)md_temp_push(sizeof(ImVec2) * num_peaks);

                double (*distr_func)(double x, double x_o, double gamma) = 0;
                // @NOTE: Do broadening in eV
                switch (broadening_mode1) {
                    case BROADENING_GAUSSIAN:
                        distr_func = &gaussian;
                        break;
                    case BROADENING_LORENTZIAN:
                        distr_func = &lorentzian;
                        break;
                    default:
                        ASSERT(false);  // Should not happen
                        break;
                }

                if (recalculate1 || rsp.first_plot_rot_ecd) {
                    osc_to_eps(rsp.eps, rsp.x_ev_samples, num_samples, y_osc_peaks, vlx.rsp.absorption_ev, num_peaks, distr_func, gamma1 * 2);
                    rot_to_eps_delta(rsp.ecd, rsp.x_ev_samples, num_samples, y_cgs_peaks, vlx.rsp.absorption_ev, num_peaks, distr_func, gamma1 * 2);
                }

                static ImPlotRect osc_lim_constraint = {0, 0, 0, 0};
                static ImPlotRect cgs_lim_constraint = {0, 0, 0, 0};
                if (refit1 || rsp.first_plot_rot_ecd) {
                    // Do conversions
                    convert_values(rsp.x_unit_peaks, vlx.rsp.absorption_ev, num_peaks, x_unit);
                    convert_values(rsp.x_unit_samples, rsp.x_ev_samples, num_samples, x_unit);

                    osc_lim_constraint = get_plot_limits(rsp.x_unit_samples, y_osc_peaks, num_peaks, num_samples);
                    cgs_lim_constraint = get_plot_limits(rsp.x_unit_samples, y_cgs_peaks, num_peaks, num_samples);
                    if (is_all_zero(y_osc_peaks, num_peaks)) {
                        osc_lim_constraint.Y.Min = -1;
                        osc_lim_constraint.Y.Max = 1;
                    }
                    if (is_all_zero(y_cgs_peaks, num_peaks)) {
                        cgs_lim_constraint.Y.Min = -1;
                        cgs_lim_constraint.Y.Max = 1;
                    }
                }

                if (rsp.first_plot_rot_ecd) {
                    rsp.osc_points = md_array_create(double, num_peaks, arena);
                    rsp.cgs_points = md_array_create(double, num_peaks, arena);
                    max_points(rsp.osc_points, y_osc_peaks, num_peaks);
                    max_points(rsp.cgs_points, y_cgs_peaks, num_peaks);
                }

#if 1
                // Hovered display text
                /*if (rsp.hovered != -1 && rsp.focused_plot == 0) {
                    ImGui::BulletText("Hovered: %s = %f, Y = %f", x_unit_str[x_unit], (float)x_peaks[rsp.hovered], (float)y_osc_peaks[rsp.hovered]);

                }
                else if (rsp.hovered != -1 && rsp.focused_plot == 1){
                    ImGui::BulletText("Hovered: %s = %f, Y = %f", x_unit_str[x_unit], (float)x_peaks[rsp.hovered], (float)y_cgs_peaks[rsp.hovered]);
                }
                else {
                    ImGui::BulletText("Hovered:");
                }*/

                // Selected display text
                if (rsp.selected != -1) {
                    ImGui::Text((const char*)u8"Selected: State %i: Energy = %.2f eV, Wavelength = %.0f nm, f = %.3f, R = %.3f 10⁻⁴⁰ cgs",
                                rsp.selected + 1, (float)rsp.x_unit_peaks[rsp.selected], 1239.84193 / (float)rsp.x_unit_peaks[rsp.selected],
                                (float)y_osc_peaks[rsp.selected], (float)y_cgs_peaks[rsp.selected]);
                } else {
                    ImGui::Text("Selected:");
                }
#endif
                rsp.focused_plot = -1;
                if (ImPlot::BeginSubplots("##AxisLinking", 2, 1, ImVec2(-1, -1), ImPlotSubplotFlags_LinkCols)) {
                    // Absorption
                    static double osc_to_eps_mult = 1;
                    if (recalculate1 || rsp.first_plot_rot_ecd) {
                        osc_to_eps_mult = is_all_zero(y_osc_peaks, num_peaks) ? 1 : axis_conversion_multiplier(y_osc_peaks, rsp.eps, num_peaks, num_samples);
                    }

                    static ImPlotRect cur_osc_lims = {0, 1, 0, 1};
                    if (refit1 || rsp.first_plot_rot_ecd) {
                        ImPlot::SetNextAxisToFit(ImAxis_X1);
                    }
                    if (ImPlot::BeginPlot("Absorption")) {
                        ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                        ImPlot::SetupAxis(ImAxis_X1, x_unit_str[x_unit]);
                        ImPlot::SetupAxis(ImAxis_Y1, "f", ImPlotAxisFlags_AuxDefault);
                        ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"ε (L mol⁻¹ cm⁻¹)");
                        if (refit1 || rsp.first_plot_rot_ecd) {
                            ImPlot::SetupAxisLimits(ImAxis_X1, osc_lim_constraint.X.Min, osc_lim_constraint.X.Max);
                            ImPlot::SetupAxisLimits(ImAxis_Y1, osc_lim_constraint.Y.Min, osc_lim_constraint.Y.Max);
                            cur_osc_lims = osc_lim_constraint;
                        }
                        ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, osc_lim_constraint.X.Min, osc_lim_constraint.X.Max);
                        ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, osc_lim_constraint.Y.Min, osc_lim_constraint.Y.Max);
                        ImPlot::SetupAxisLimits(ImAxis_Y2, cur_osc_lims.Y.Min * osc_to_eps_mult, cur_osc_lims.Y.Max * osc_to_eps_mult,
                                                ImPlotCond_Always);
                        ImPlot::SetupFinish();

                        peaks_to_pixels(pixel_osc_peaks, rsp.x_unit_peaks, y_osc_peaks, num_peaks);
                        peaks_to_pixels(pixel_osc_points, rsp.x_unit_peaks, rsp.osc_points, num_peaks);
                        mouse_pos = ImPlot::PlotToPixels(ImPlot::GetPlotMousePos(IMPLOT_AUTO));
                        if (ImPlot::IsPlotHovered()) {
                            rsp.hovered = get_hovered_peak(mouse_pos, pixel_osc_peaks, pixel_osc_points, num_peaks);
                            rsp.focused_plot = 0;
                        }

                        // @HACK: Compute pixel width of 2 'plot' units
                        const double bar_width = ImPlot::PixelsToPlot(ImVec2(2, 0)).x - ImPlot::PixelsToPlot(ImVec2(0, 0)).x;

                        ImPlot::SetAxis(ImAxis_Y2);
                        ImPlot::PlotLine("Spectrum", rsp.x_unit_samples, rsp.eps, num_samples);
                        ImPlot::SetAxis(ImAxis_Y1);
                        ImPlot::PlotBars("Oscillator Strength", rsp.x_unit_peaks, y_osc_peaks, num_peaks, bar_width);
                        ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 3);
                        ImPlot::PlotScatter("##Peak marker", rsp.x_unit_peaks, rsp.osc_points, num_peaks);

                        // Check hovered state
                        if (rsp.hovered != -1) {
                            draw_bar(0, rsp.x_unit_peaks[rsp.hovered], y_osc_peaks[rsp.hovered], bar_width, IM_GREEN);
                            ImPlot::DragPoint(0, &rsp.x_unit_peaks[rsp.hovered], &rsp.osc_points[rsp.hovered], IM_GREEN, 4, ImPlotDragToolFlags_NoInputs);
                        }

                        // Update selected peak on click
                        if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Left) &&
                            ImPlot::IsPlotHovered()) {
                            rsp.selected = rsp.hovered == rsp.selected ? -1 : rsp.hovered;
                        }
                        // Check selected state
                        if (rsp.selected != -1) {
                            draw_bar(1, rsp.x_unit_peaks[rsp.selected], y_osc_peaks[rsp.selected], bar_width, IM_RED);
                            ImPlot::DragPoint(0, &rsp.x_unit_peaks[rsp.selected], &rsp.osc_points[rsp.selected], IM_RED, 4, ImPlotDragToolFlags_NoInputs);
                        }

                        cur_osc_lims = ImPlot::GetPlotLimits(ImAxis_X1, ImAxis_Y1);
                        ImPlot::EndPlot();
                    }

                    // Rotatory ECD
                    static double cgs_to_ecd_mult = 1;
                    if (recalculate1 || rsp.first_plot_rot_ecd) {
                        cgs_to_ecd_mult = is_all_zero(y_cgs_peaks, num_peaks) ? 1 : axis_conversion_multiplier(y_cgs_peaks, rsp.ecd, num_peaks, num_samples);
                    }
                    static ImPlotRect cur_cgs_lims = {0, 1, 0, 1};
                    if (refit1 || rsp.first_plot_rot_ecd) {
                        ImPlot::SetNextAxisToFit(ImAxis_X1);
                    }

                    if (ImPlot::BeginPlot("ECD")) {
                        ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                        ImPlot::SetupAxis(ImAxis_X1, x_unit_str[x_unit]);
                        ImPlot::SetupAxis(ImAxis_Y1, (const char*)u8"R (10⁻⁴⁰ cgs)", ImPlotAxisFlags_AuxDefault);
                        ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"Δε(ω) (L mol⁻¹ cm⁻¹)");
                        if (refit1 || rsp.first_plot_rot_ecd) {
                            ImPlot::SetupAxisLimits(ImAxis_X1, cgs_lim_constraint.X.Min, cgs_lim_constraint.X.Max);
                            ImPlot::SetupAxisLimits(ImAxis_Y1, cgs_lim_constraint.Y.Min, cgs_lim_constraint.Y.Max);
                            cur_cgs_lims = cgs_lim_constraint;
                        }
                        ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, cgs_lim_constraint.X.Min, cgs_lim_constraint.X.Max);
                        ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, cgs_lim_constraint.Y.Min, cgs_lim_constraint.Y.Max);
                        ImPlot::SetupAxisLimits(ImAxis_Y2, cur_cgs_lims.Y.Min * cgs_to_ecd_mult, cur_cgs_lims.Y.Max * cgs_to_ecd_mult,
                                                ImPlotCond_Always);
                        ImPlot::SetupFinish();

                        peaks_to_pixels(pixel_cgs_peaks, rsp.x_unit_peaks, y_cgs_peaks, num_peaks);
                        peaks_to_pixels(pixel_cgs_points, rsp.x_unit_peaks, rsp.cgs_points, num_peaks);
                        mouse_pos = ImPlot::PlotToPixels(ImPlot::GetPlotMousePos(IMPLOT_AUTO));

                        if (ImPlot::IsPlotHovered()) {
                            rsp.hovered = get_hovered_peak(mouse_pos, pixel_cgs_peaks, pixel_cgs_points, num_peaks);
                            rsp.focused_plot = 1;
                        }
                        // @HACK: Compute pixel width of 2 'plot' units

                        const double bar_width = ImPlot::PixelsToPlot(ImVec2(2, 0)).x - ImPlot::PixelsToPlot(ImVec2(0, 0)).x;

                        ImPlot::SetAxis(ImAxis_Y2);
                        ImPlot::PlotLine("Spectrum", rsp.x_unit_samples, rsp.ecd, num_samples);
                        ImPlot::SetAxis(ImAxis_Y1);
                        ImPlot::PlotBars("Rotatory Strength", rsp.x_unit_peaks, y_cgs_peaks, num_peaks, bar_width);
                        ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 3);
                        ImPlot::PlotScatter("##Peak marker", rsp.x_unit_peaks, rsp.cgs_points, num_peaks);

                        if (rsp.hovered != -1) {
                            draw_bar(2, rsp.x_unit_peaks[rsp.hovered], y_cgs_peaks[rsp.hovered], bar_width, IM_GREEN);
                            ImPlot::DragPoint(0, &rsp.x_unit_peaks[rsp.hovered], &rsp.cgs_points[rsp.hovered], IM_GREEN, 4, ImPlotDragToolFlags_NoInputs);

                        }

                        // Update selected peak on click
                        if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Left) &&
                            ImPlot::IsPlotHovered()) {
                            rsp.selected = rsp.hovered == rsp.selected ? -1 : rsp.hovered;
                        }
                        if (rsp.selected != -1) {
                            draw_bar(3, rsp.x_unit_peaks[rsp.selected], y_cgs_peaks[rsp.selected], bar_width, IM_RED);
                            ImPlot::DragPoint(0, &rsp.x_unit_peaks[rsp.selected], &rsp.cgs_points[rsp.selected], IM_RED, 4, ImPlotDragToolFlags_NoInputs);

                        }
                        cur_cgs_lims = ImPlot::GetPlotLimits(ImAxis_X1, ImAxis_Y1);
                        ImPlot::EndPlot();
                    }
                    ImPlot::EndSubplots();
                }
                rsp.first_plot_rot_ecd = false;
                ImGui::TreePop();
            }
#if 0
            if (rsp.first_plot_vib) { ImGui::SetNextItemOpen(true); }
            if (ImGui::TreeNode("Vibrational Analysis")) {
                // draw the vibrational analysis
                double har_freqs[3] = {1562.20, 3663.36, 3677.39};
                double irs[3] = {132.6605, 14.2605, 5.8974};

                double x0[3] = {0, 0, 0};
                double y0[3] = {0, 0.4272, -0.4272};
                double z0[3] = {-0.0707, 0.5612, 0.5612};

                double x1[3] = {0, 0, 0};
                double y1[3] = {0.0701, -0.5563, -0.5563};
                double z1[3] = {0, 0.4337, -0.4337};

                double x2[3] = {0, 0, 0};
                double y2[3] = {0, -0.5851, 0.5851};
                double z2[3] = {-0.0498, 0.3955, 0.3955};

                vibration_mode vib_modes[3] = {
                    {1562.20, 0, 0, 132.6605, x0, y0, z0},
                    {3663.36, 0, 0, 14.2605, x1, y1, z1},
                    {3677.39, 0, 0, 5.8974, x2, y2, z2},
                };

                // ASSERT(ARRAY_SIZE(har_freqs) == ARRAY_SIZE(irs));
                size_t num_vibs = ARRAY_SIZE(vib_modes);
                size_t num_atoms = 3;
                static int hov_vib = -1;
                static int sel_vib = -1;

                bool refit2 = false;
                bool recalculate2 = false;
                static float gamma2 = 5.0f;
                static broadening_mode_t broadening_mode2 = BROADENING_LORENTZIAN;
                recalculate2 = ImGui::SliderFloat((const char*)u8"Broadening γ HWHM (cm⁻¹)", &gamma2, 1.0f, 10.0f);
                refit2 |= ImGui::Combo("Broadening mode", (int*)(&broadening_mode2), broadening_str, IM_ARRAYSIZE(broadening_str));



                ImVec2* pixel_peaks = (ImVec2*)md_temp_push(sizeof(ImVec2) * num_vibs);
                ImVec2* pixel_points = (ImVec2*)md_temp_push(sizeof(ImVec2) * num_vibs);

                static bool coord_modified = false;
                static float amp_mult = 1;
                static float speed_mult = 1;
                static float time = 0;

                double (*distr_func)(double x, double x_o, double gamma, double intensity) = 0;
                switch (broadening_mode2) {
                    case BROADENING_GAUSSIAN:
                        distr_func = &phys_gaussian;
                        break;
                    case BROADENING_LORENTZIAN:
                        distr_func = &phys_lorentzian;
                        break;
                    default:
                        ASSERT(false);  // Should not happen
                        break;
                }

                if (rsp.first_plot_vib) {
                    rsp.vib_x = md_array_create(double, num_samples, arena);
                    rsp.vib_y = md_array_create(double, num_samples, arena);


                    // Populate x_values
                    const double x_min = har_freqs[0] - 100.0;
                    const double x_max = har_freqs[num_vibs - 1] + 100.0;
                    for (int i = 0; i < num_samples; ++i) {
                        double t = (double)i / (double)(num_samples - 1);
                        double value = lerp(x_min, x_max, t);
                        rsp.vib_x[i] = value;
                    }
                }

                if (rsp.first_plot_vib || recalculate2 || refit2) {
                    general_broadening(rsp.vib_y, rsp.vib_x, num_samples, irs, har_freqs, num_vibs, distr_func, gamma2 * 2);
                }

                if (rsp.first_plot_vib) {
                    rsp.vib_points = md_array_create(double, num_vibs, arena);
                    max_points(rsp.vib_points, irs, num_vibs);
                }
                static bool invert_x = false;
                static bool invert_y = false;
                ImGui::Checkbox("Invert X", &invert_x); ImGui::SameLine();
                ImGui::Checkbox("Invert Y", &invert_y);

                ImPlotAxisFlags x_flag = invert_x ? ImPlotAxisFlags_Invert : 0;
                ImPlotAxisFlags y_flag = invert_y ? ImPlotAxisFlags_Invert : 0;

                static ImPlotRect lim_constraint = { 0, 0, 0, 0 };
                if (refit2 || rsp.first_plot_vib) {
                    lim_constraint = get_plot_limits(rsp.vib_x, irs, num_vibs, num_samples);
                }

                if (ImPlot::BeginPlot("Vibrational analysis")) {
                    // @HACK: Compute pixel width of 2 'plot' units
                    ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                    ImPlot::SetupAxis(ImAxis_X1, (const char*)u8"Harmonic Frequency (cm⁻¹)", x_flag);
                    ImPlot::SetupAxis(ImAxis_Y1, "IR Intensity (km/mol)", y_flag);
                    if (refit2 || rsp.first_plot_vib) {
                        ImPlot::SetupAxisLimits(ImAxis_X1, lim_constraint.X.Min, lim_constraint.X.Max);
                        ImPlot::SetupAxisLimits(ImAxis_Y1, lim_constraint.Y.Min, lim_constraint.Y.Max);
                    }
                    ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, lim_constraint.X.Min, lim_constraint.X.Max);
                    ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, lim_constraint.Y.Min, lim_constraint.Y.Max);
                    ImPlot::SetupFinish();

                    ImPlot::PlotLine("Spectrum", rsp.vib_x, rsp.vib_y, num_samples);

                    const double bar_width = ImPlot::PixelsToPlot(ImVec2(2, 0)).x - ImPlot::PixelsToPlot(ImVec2(0, 0)).x;
                    ImPlot::PlotBars("IR Intensity", har_freqs, irs, (int)num_vibs, bar_width);

                    ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 3);
                    ImPlot::PlotScatter("##Peak markers", har_freqs, rsp.vib_points, (int)num_vibs);

                    peaks_to_pixels(pixel_peaks, har_freqs, irs, num_vibs);
                    peaks_to_pixels(pixel_points, har_freqs, rsp.vib_points, num_vibs);
                    mouse_pos = ImPlot::PlotToPixels(ImPlot::GetPlotMousePos(IMPLOT_AUTO));
                    if (ImPlot::IsPlotHovered()) {
                        hov_vib = get_hovered_peak(mouse_pos, pixel_peaks, pixel_points, num_vibs, invert_y);
                    }

                    // Check hovered state
                    if (hov_vib != -1) {
                        draw_bar(0, har_freqs[hov_vib], irs[hov_vib], bar_width, IM_GREEN);
                        ImPlot::DragPoint(0, &har_freqs[hov_vib], &rsp.vib_points[hov_vib], IM_GREEN, 4, ImPlotDragToolFlags_NoInputs);
                    }

                    // Update selected peak on click
                    if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Left) &&
                        ImPlot::IsPlotHovered()) {
                        sel_vib = hov_vib == sel_vib ? -1 : hov_vib;
                    }
                    // Check selected state
                    if (sel_vib != -1) {
                        draw_bar(1, har_freqs[sel_vib], irs[sel_vib], bar_width, IM_RED);
                        ImPlot::DragPoint(1, &har_freqs[sel_vib], &rsp.vib_points[sel_vib], IM_RED, 4, ImPlotDragToolFlags_NoInputs);

                        //Animation
                        time += state.app.timing.delta_s * speed_mult * 7;
                        for (size_t id = 0; id < num_atoms; id++) {
                            state.mold.mol.atom.x[id] = vlx.geom.coord_x[id] + amp_mult * 0.5 * vib_modes[sel_vib].x[id] * sin(time);
                            state.mold.mol.atom.y[id] = vlx.geom.coord_y[id] + amp_mult * 0.5 * vib_modes[sel_vib].y[id] * sin(time);
                            state.mold.mol.atom.z[id] = vlx.geom.coord_z[id] + amp_mult * 0.5 * vib_modes[sel_vib].z[id] * sin(time);
                        }
                        state.mold.dirty_buffers |= MolBit_DirtyPosition;
                        coord_modified = true;
                    }
                    // If all is deselected, reset coords once
                    else if (coord_modified) {
                        for (size_t id = 0; id < num_atoms; id++) {
                            state.mold.mol.atom.x[id] = (float)vlx.geom.coord_x[id];
                            state.mold.mol.atom.y[id] = (float)vlx.geom.coord_y[id];
                            state.mold.mol.atom.z[id] = (float)vlx.geom.coord_z[id];
                        }
                        state.mold.dirty_buffers |= MolBit_DirtyPosition | MolBit_ClearVelocity;
                        coord_modified = false;
                    }
                    rsp.first_plot_vib = false;
                    ImPlot::EndPlot();
                }

                // ImGui::Text("%i is hovered", hov_vib);
                // ImGui::Text("%f is z coord", (float)state.mold.mol.atom.z[2]);

                ImGui::SliderFloat((const char*)"Amplitude", &amp_mult, 0.2f, 2.0f);
                ImGui::SliderFloat((const char*)"Speed", &speed_mult, 0.5f, 2.0f);

                //Table
                static ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders |
                    ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_Sortable;

                static ImGuiTableColumnFlags columns_base_flags = ImGuiTableColumnFlags_DefaultSort;

                double height = ImGui::GetTextLineHeightWithSpacing();

                if (ImGui::BeginTable("table_advanced", 3, flags, ImVec2(460, height * (num_vibs + 1)), 0)) {
                    ImGui::TableSetupColumn("Vibration mode", columns_base_flags, 0.0f);
                    ImGui::TableSetupColumn("Harmonic Frequency", columns_base_flags, 0.0f);
                    ImGui::TableSetupColumn("IR Intensity", columns_base_flags, 0.0f);
                    ImGui::TableSetupScrollFreeze(0, 1);
                    ImGui::TableHeadersRow();

                    //TODO: Add sorting to the table once the vibs data structure is properly defined

                    ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);
                    ImGui::PushStyleColor(ImGuiCol_Header, IM_BLUE);
                    bool item_hovered = false;
                    for (int row_n = 0; row_n < num_vibs; row_n++) {

                        ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                        bool is_sel = row_n == sel_vib; //If atom is selected, mark it as such
                        bool is_hov = row_n == hov_vib; //If atom is hovered, mark it as such
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
                            sel_vib = sel_vib == row_n ? -1 : row_n;
                        }
                        ImGui::TableNextColumn();
                        ImGui::Text("%12.6f", har_freqs[row_n]);
                        ImGui::TableNextColumn();
                        ImGui::Text("%12.6f", irs[row_n]);

                        ImGui::PopStyleColor(1);

                    }
                    if (ImGui::IsWindowHovered() && ImGui::TableGetHoveredRow() > 0) {
                        hov_vib = ImGui::TableGetHoveredRow() - 1;
                    }
                    else {
                        hov_vib = -1;
                    }

                    ImGui::PopStyleColor(2);
                    ImGui::EndTable();
                }

                ImGui::TreePop();
            }
#endif
        }
        ImGui::End();

        if (rsp.show_export_window) { draw_rsp_spectra_export_window(state); }
    }

    void draw_orb_window(const ApplicationState& state) {
        if (!orb.show_window) return;
        if (num_orbitals() == 0) return;
        ImGui::SetNextWindowSize({600,300}, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("VeloxChem Orbital Grid", &orb.show_window)) {
#if 0
            if (vlx.geom.num_atoms) {
                if (ImGui::TreeNode("Geometry")) {
                    ImGui::Text("Num Atoms:           %6zu", vlx.geom.num_atoms);
                    ImGui::Text("Num Alpha Electrons: %6zu", vlx.geom.num_alpha_electrons);
                    ImGui::Text("Num Beta Electrons:  %6zu", vlx.geom.num_beta_electrons);
                    ImGui::Text("Molecular Charge:    %6i",  vlx.geom.molecular_charge);
                    ImGui::Text("Spin Multiplicity:   %6i",  vlx.geom.spin_multiplicity);
                    ImGui::Spacing();
                    ImGui::Text("Atom      Coord X      Coord Y      Coord Z");
                    for (size_t i = 0; i < vlx.geom.num_atoms; ++i) {
                        ImGui::Text("%4s %12.6f %12.6f %12.6f", vlx.geom.atom_symbol[i].buf, vlx.geom.coord_x[i], vlx.geom.coord_y[i], vlx.geom.coord_z[i]);
                    }
                    ImGui::TreePop();
                }
            }
#endif
            const ImVec2 outer_size = {300.f, 0.f};
            ImGui::PushItemWidth(outer_size.x);
            ImGui::BeginGroup();

            ImGui::SliderInt("##Rows", &orb.num_y, 1, 4);
            ImGui::SliderInt("##Cols", &orb.num_x, 1, 4);

            const int num_mos = orb.num_x * orb.num_y;
            const int beg_mo_idx = orb.mo_idx - num_mos / 2 + (num_mos % 2 == 0 ? 1 : 0);

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

            const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();
            enum {
                Col_Idx,
                Col_Occ,
                Col_Ene,
            };

            if (ImGui::IsWindowAppearing()) {
                orb.scroll_to_idx = orb.mo_idx;
            }
            if (ImGui::Button("Goto HOMO", ImVec2(outer_size.x,0))) {
                orb.scroll_to_idx = homo_idx;
            }

            const ImGuiTableFlags flags =
                ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable | ImGuiTableFlags_RowBg |
                ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersV | ImGuiTableFlags_NoBordersInBody | ImGuiTableFlags_ScrollY;
            if (ImGui::BeginTable("Molecular Orbitals", 3, flags, outer_size))//, ImVec2(0.0f, TEXT_BASE_HEIGHT * 15), 0.0f))
            {
                // Declare columns
                // We use the "user_id" parameter of TableSetupColumn() to specify a user id that will be stored in the sort specifications.
                // This is so our sort function can identify a column given our own identifier. We could also identify them based on their index!
                // Demonstrate using a mixture of flags among available sort-related flags:
                // - ImGuiTableColumnFlags_DefaultSort
                // - ImGuiTableColumnFlags_NoSort / ImGuiTableColumnFlags_NoSortAscending / ImGuiTableColumnFlags_NoSortDescending
                // - ImGuiTableColumnFlags_PreferSortAscending / ImGuiTableColumnFlags_PreferSortDescending
                ImGui::TableSetupColumn("Index",        ImGuiTableColumnFlags_DefaultSort          | ImGuiTableColumnFlags_WidthFixed,   0.0f, Col_Idx);
                ImGui::TableSetupColumn("Occupation",   ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed,   0.0f, Col_Occ);
                ImGui::TableSetupColumn("Energy",       ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed,   0.0f, Col_Ene);
                ImGui::TableSetupScrollFreeze(0, 1); // Make row always visible
                ImGui::TableHeadersRow();

                for (int n = (int)num_orbitals() - 1; n >= 0; n--) {
                    ImGui::PushID(n + 1);
                    ImGui::TableNextRow();
                    bool is_selected = (beg_mo_idx <= n && n < beg_mo_idx + num_mos);
                    ImGui::TableNextColumn();
                    if (orb.scroll_to_idx != -1 && n == orb.scroll_to_idx) {
                        orb.scroll_to_idx = -1;
                        ImGui::SetScrollHereY();
                    }
                    char buf[32];
                    const char* lbl = (n == homo_idx) ? " (HOMO)" : (n == lumo_idx) ? " (LUMO)" : "";
                    snprintf(buf, sizeof(buf), "%i%s", n + 1, lbl);
                    ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                    if (ImGui::Selectable(buf, is_selected, selectable_flags)) {
                        if (orb.mo_idx != n) {
                            orb.mo_idx = n;
                        }
                    }
                    ImGui::TableNextColumn();
                    ImGui::Text("%.1f", vlx.scf.alpha.occupations.data[n]);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.4f", vlx.scf.alpha.energies.data[n]);
                    ImGui::PopID();
                }

                ImGui::EndTable();
            }

            ImGui::EndGroup();
            ImGui::PopItemWidth();

            ImGui::SameLine();

            // These represent the new mo_idx we want to have in each slot
            int vol_mo_idx[16] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
            for (int i = 0; i < num_mos; ++i) {
                int mo_idx = beg_mo_idx + i;
                if (-1 < mo_idx && mo_idx < num_orbitals()) {
                    vol_mo_idx[i] = mo_idx;
                }
            }

            int job_queue[16];
            int num_jobs = 0;
            // Find and reuse volume data from existing slots (if applicable)
            // If there is no existing volume, we queue up a new job
            for (int i = 0; i < num_mos; ++i) {
                // Check if we already have that entry in the correct slot
                if (orb.vol_mo_idx[i] == vol_mo_idx[i]) continue;

                // Try to find the entry in the existing list
                bool found = false;
                for (int j = 0; j < num_mos; ++j) {
                    if (i == j) continue;
                    if (vol_mo_idx[i] == orb.vol_mo_idx[j]) {
                        // Swap to correct location
                        ImSwap(orb.vol[i], orb.vol[j]);
                        ImSwap(orb.vol_mo_idx[i], orb.vol_mo_idx[j]);
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
                const float samples_per_angstrom = 6.0f;
                for (int i = 0; i < num_jobs; ++i) {
                    int slot_idx = job_queue[i];
                    int mo_idx = vol_mo_idx[slot_idx];
                    orb.vol_mo_idx[slot_idx] = mo_idx;

                    if (-1 < mo_idx && mo_idx < num_orbitals()) {
                        if (task_system::task_is_running(orb.vol_task[slot_idx])) {
                            task_system::task_interrupt(orb.vol_task[slot_idx]);
                        }
                        orb.vol_task[slot_idx] = compute_mo_async(&orb.vol[slot_idx].tex_to_world, &orb.vol[slot_idx].step_size, &orb.vol[slot_idx].tex_id, mo_idx, MD_GTO_EVAL_MODE_PSI, samples_per_angstrom);
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

            ImVec2 orb_win_sz = (canvas_p1 - canvas_p0) / ImVec2((float)orb.num_x, (float)orb.num_y);
            orb_win_sz.x = floorf(orb_win_sz.x);
            orb_win_sz.y = floorf(orb_win_sz.y);
            canvas_p1.x = canvas_p0.x + orb.num_x * orb_win_sz.x;
            canvas_p1.y = canvas_p0.y + orb.num_y * orb_win_sz.y;

            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(255, 255, 255, 255));
            for (int i = 0; i < num_mos; ++i) {
                int mo_idx = beg_mo_idx + i;
                int x = orb.num_x - i % orb.num_x - 1;
                int y = orb.num_y - i / orb.num_x - 1;
                ImVec2 p0 = canvas_p0 + orb_win_sz * ImVec2((float)(x+0), (float)(y+0));
                ImVec2 p1 = canvas_p0 + orb_win_sz * ImVec2((float)(x+1), (float)(y+1));
                if (-1 < mo_idx && mo_idx < num_orbitals()) {
                    ImVec2 text_pos = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p1.y - TEXT_BASE_HEIGHT);
                    char buf[32];
                    const char* lbl = (mo_idx == homo_idx) ? " (HOMO)" : (mo_idx == lumo_idx) ? " (LUMO)" : "";
                    snprintf(buf, sizeof(buf), "%i%s", mo_idx + 1, lbl);
                    draw_list->AddImage((ImTextureID)(intptr_t)orb.gbuf.tex.transparency, p0, p1, { 0,1 }, { 1,0 });
                    draw_list->AddImage((ImTextureID)(intptr_t)orb.iso_tex[i], p0, p1, { 0,1 }, { 1,0 });
                    draw_list->AddText(text_pos, ImColor(0,0,0), buf);
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
                camera_compute_optimal_view(&orb.target.pos, &orb.target.ori, &orb.target.dist, min_aabb, max_aabb, orb.distance_scale);

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
                                    .model = orb.vol[i].tex_to_world,
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
                                .voxel_spacing = orb.vol[i].step_size,
                            };
                            volume::render_volume(vol_desc);
                        }
                    POP_GPU_SECTION();
                }
            }
        }
        ImGui::End();
    }

    //Calculates the transition matrix heuristic
    static inline void distribute_charges_heuristic(float* out_matrix, const size_t num_groups, const float* hole_charges, const float* particle_charges) {
        size_t temp_pos = md_temp_get_pos();
        md_allocator_i* temp_alloc = md_get_temp_allocator();
        int* donors = 0;
        int* acceptors = 0;
        float* charge_diff = 0;

        float hole_sum = 0;
        float part_sum = 0;
        for (size_t i = 0; i < num_groups; i++) {
            hole_sum += hole_charges[i];
            part_sum += particle_charges[i];
        }

        md_array(float) hole_percentages     = md_array_create(float, num_groups, temp_alloc);
        md_array(float) particle_percentages = md_array_create(float, num_groups, temp_alloc);

        for (size_t i = 0; i < num_groups; i++) {
            hole_percentages[i] = hole_charges[i] / hole_sum;
            particle_percentages[i] = particle_charges[i] / part_sum;
        }

        for (size_t i = 0; i < num_groups; i++) {
            float gsCharge = hole_percentages[i];
            float esCharge = particle_percentages[i];
            if (gsCharge > esCharge) {
                md_array_push(donors, (int)i, temp_alloc);
            }
            else {
                md_array_push(acceptors, (int)i, temp_alloc);
            }
            float diff = esCharge - gsCharge;
            out_matrix[i * num_groups + i] = MIN(gsCharge, esCharge);
            md_array_push(charge_diff, diff, temp_alloc);
        }

        int num_donors = md_array_size(donors);
        int num_acceptors = md_array_size(acceptors);

        float total_acceptor_charge = 0;
        for (size_t i = 0; i < md_array_size(acceptors); i++) {
            total_acceptor_charge += charge_diff[acceptors[i]];
        }
        for (size_t don_i = 0; don_i < md_array_size(donors); don_i++) {
            float charge_deficit = -charge_diff[donors[don_i]];
            for (size_t acc_i = 0; acc_i < md_array_size(acceptors); acc_i++) {
                float contrib = charge_deficit * charge_diff[acceptors[acc_i]] / total_acceptor_charge;
                out_matrix[acceptors[acc_i] * num_groups + donors[don_i]] = contrib;
            }
        }
        md_temp_set_pos_back(temp_pos);
    }

    //Takes the hole and particle charges of all atoms, and calculates the per group charges
    static inline void accumulate_subgroup_charges(const float* hole_charges, const float* particle_charges, size_t num_charges, size_t num_subgroups, float* ligandGSCharges, float* ligandESCharges, int* atom_subgroup_map) {
        md_allocator_i* temp_alloc = md_get_temp_allocator();
        float sumGSCharges = 0;
        float sumESCharges = 0;
        for (size_t i = 0; i < num_charges; i++) {
            sumGSCharges += hole_charges[i];
            sumESCharges += particle_charges[i];
        }

        for (size_t i = 0; i < num_charges; i++) {
            int subgroup_index = atom_subgroup_map[i];
            ligandGSCharges[subgroup_index] += hole_charges[i] / sumGSCharges;
            ligandESCharges[subgroup_index] += particle_charges[i] / sumESCharges;
        }
    }

    //Calculates the subgroup charges and the transition matrix
    static inline void compute_subgroup_charges(float* hole_charges, float* particle_charges, size_t num_charges, size_t num_subgroups, int* atom_subgroup_map) {
        md_allocator_i* temp_alloc = md_get_temp_allocator();

        //These two arrays are the group charges. They are already defined in the GroupData.
        float* ligandGSCharges = md_array_create(float, num_subgroups, temp_alloc);
        md_array(float) ligandESCharges = md_array_create(float, num_subgroups, temp_alloc);
        accumulate_subgroup_charges(hole_charges, particle_charges, num_charges, num_subgroups, ligandGSCharges, ligandESCharges, atom_subgroup_map);
        

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

    static void interaction_canvas(ImVec2 size, SelectionState& select, ViewState& view, const md_molecule_t& mol) {
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

            ImVec2 coord = ImGui::GetMousePos() - ImGui::GetWindowPos();

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
                                md_bond_pair_t pair = mol.bond.pairs[select.hovered_bond_idx];
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
                        md_bond_pair_t pair = mol.bond.pairs[select.hovered_bond_idx];
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

    static void update_picking_data(PickingData& picking, const vec2_t& coord, GBuffer& gbuffer, mat4_t inv_MVP) {
        picking = {};

#if MD_PLATFORM_OSX
        coord = coord * vec_cast(ImGui::GetIO().DisplayFramebufferScale);
#endif
        if (0.f < coord.x && coord.x < (float)gbuffer.width && 0.f < coord.y && coord.y < (float)gbuffer.height) {
            extract_picking_data(&picking.idx, &picking.depth, &gbuffer, (int)coord.x, (int)coord.y);
            const vec4_t viewport = {0, 0, (float)gbuffer.width, (float)gbuffer.height};
            picking.world_coord = mat4_unproject({coord.x, coord.y, picking.depth}, inv_MVP, viewport);
            picking.screen_coord = {coord.x, coord.y};
        }
    }

    void draw_nto_window(ApplicationState& state) {
        if (!nto.show_window) return;
        if (vlx.rsp.num_excited_states == 0) return;
        if (vlx.rsp.nto == NULL) return;

        bool open_context_menu = false;

        ImGui::SetNextWindowSize(ImVec2(500, 300), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("NTO viewer", &nto.show_window, ImGuiWindowFlags_MenuBar)) {

            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("Settings")) {
                    ImGui::Text("Orbital Colors");
                    ImGui::ColorEdit4("##Color Positive", nto.iso.colors[0].elem);
                    ImGui::SetItemTooltip("Color Positive");
                    ImGui::ColorEdit4("##Color Negative", nto.iso.colors[1].elem);
                    ImGui::SetItemTooltip("Color Negative");

                    ImGui::Text("Transition Dipole Moments");
                    ImGui::Spacing();
                    const double vector_length_min = 1.0;
                    const double vector_length_max = 10.0;
                    double vector_length_input = nto.iso.vector_length;
                    ImGui::Text("Scaling");
                    ImGui::SliderScalar("##Vector length", ImGuiDataType_Double, &vector_length_input, &vector_length_min, &vector_length_max, "%.6f",
                        ImGuiSliderFlags_Logarithmic);
                    ImGui::SetItemTooltip("Vector length");
                    nto.iso.vector_length = (float)vector_length_input;

                    bool show_vector = nto.iso.display_vectors;
                    ImGui::Checkbox("Display transition dipole moments", &show_vector);
                    nto.iso.display_vectors = show_vector;

                    bool show_angle = nto.iso.display_angle;
                    ImGui::Checkbox("Display angle", &show_angle);
                    nto.iso.display_angle = show_angle;
                    ImGui::Text("Vectors and angle Colors");
                    ImGui::ColorEdit4("##Color Electric", nto.iso.vectorColors[0].elem);
                    ImGui::SetItemTooltip("Color Electric");
                    ImGui::ColorEdit4("##Color Magnetic", nto.iso.vectorColors[1].elem);
                    ImGui::SetItemTooltip("Color Magnetic");
                    ImGui::ColorEdit4("##Color Angle", nto.iso.vectorColors[2].elem);
                    ImGui::SetItemTooltip("Color Angle");
                    ImGui::EndMenu();
                }
                ImGui::EndMenuBar();
            }

            const ImVec2 outer_size = {300.f, 0.f};
            ImGui::PushItemWidth(outer_size.x);
            ImGui::BeginGroup();

            ImGui::Text("Transition State Index");
            ImGui::Spacing();
            if (ImGui::BeginListBox("##NTO Index", outer_size)) {
                if (ImGui::IsWindowHovered()) {
                    rsp.hovered = -1;
                }
                for (int i = 0; i < (int)vlx.rsp.num_excited_states; ++i) {
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
            const double iso_min = 1.0e-4;
            const double iso_max = 5.0;
            double iso_val = nto.iso.values[0];             
            ImGui::Spacing();
            ImGui::Text("Isovalue"); 
            ImGui::SliderScalar("##Iso Value", ImGuiDataType_Double, &iso_val, &iso_min, &iso_max, "%.6f", ImGuiSliderFlags_Logarithmic);
            ImGui::SetItemTooltip("Iso Value");

            nto.iso.values[0] = (float)iso_val;
            nto.iso.values[1] = -(float)iso_val;
            nto.iso.count = 2;
            nto.iso.enabled = true;

            ImGui::Spacing();

            // @TODO: Enlist all defined groups here
            if (ImGui::BeginListBox("##Groups", outer_size)) {
                bool item_hovered = false;
                for (size_t i = 0; i < nto.group.count; i++) {
                    char color_buf[16];
                    sprintf(color_buf, "##Group-Color%i", (int)i);
                    ImGui::ColorEdit4Minimal(color_buf, nto.group.color[i].elem); 
                    ImGui::SameLine(); 
                    if (ImGui::IsKeyDown(ImGuiKey_LeftCtrl)) {
                        ImGui::Selectable(nto.group.label[i]);
                        if (ImGui::IsItemHovered()) {
                            item_hovered = true;
                            md_bitfield_clear(&state.selection.highlight_mask);
                            for (size_t j = 0; j < nto.num_atoms; j++) {
                                if (nto.atom_group_idx[j] == i) {
                                    //Add j to highlight mask
                                    md_bitfield_set_bit(&state.selection.highlight_mask, j);
                                }
                            }
                            //Select
                            if (ImGui::IsKeyDown(ImGuiKey_MouseLeft)) {
                                md_bitfield_or_inplace(&state.selection.selection_mask, &state.selection.highlight_mask);
                            }
                            //Deselect
                            else if (ImGui::IsKeyDown(ImGuiKey_MouseRight)) {
                                md_bitfield_andnot_inplace(&state.selection.selection_mask, &state.selection.highlight_mask);
                            }
                        }
                        else if (!item_hovered && ImGui::IsWindowHovered()) {
                            md_bitfield_clear(&state.selection.highlight_mask);
                        }
                    }
                    else if (i == 0) {
                        ImGui::Text(nto.group.label[i]);
                    }
                    else {
                        ImGui::InputText(color_buf, nto.group.label[i], 16);
                    }

                    ImGui::SameLine();
                    int atom_count = 0;
                    for (size_t k = 0; k < nto.num_atoms; k++) {
                        if (nto.atom_group_idx[k] == i) {
                            atom_count++;
                        }
                    }
                    ImGui::Text("Atoms %i", atom_count);
                    //if (ImGui::Selectable(nto.group.label[i])) {
                    //    for (size_t j = 0; j < nto.num_atoms; j++) {
                    //        if (nto.atom_group_idx[j] == i) {
                    //            md_bitfield_set_bit(&state.selection.selection_mask, j);
                    //            ////Selection
                    //            //if (ImGui::IsKeyDown(ImGuiKey_MouseLeft) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                    //            //    md_bitfield_set_bit(&state.selection.selection_mask, j);
                    //            //}
                    //            ////Deselect
                    //            //else if (ImGui::IsKeyDown(ImGuiKey_MouseRight) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                    //            //    md_bitfield_clear_bit(&state.selection.selection_mask, j);
                    //            //}
                    //        }
                    //    }
                    //}
                }
                ImGui::EndListBox();
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

            double nto_lambda[4] = {};
            int num_lambdas = 1;
            bool reset_view = false;

            if (rsp.selected != -1) {
                // This represents the cutoff for contributing orbitals to be part of the orbital 'grid'
                // If the occupation parameter is less than this it will not be displayed
                const double lambda_cutoff = 0.10f;
                for (size_t i = 0; i < MIN(ARRAY_SIZE(nto_lambda), vlx.rsp.nto[rsp.selected].occupations.count); ++i) {
                    nto_lambda[i] = vlx.rsp.nto[rsp.selected].occupations.data[lumo_idx + i];
                    if (nto_lambda[i] < lambda_cutoff) {
                        num_lambdas = (int)i;
                        break;
                    }
                }


				// This should be triggered by a change in the selected NTO indxex
                // And when groups are changed
                bool recompute_transition_matrix = false;

                if (nto.vol_nto_idx != rsp.selected) {
                    nto.vol_nto_idx = rsp.selected;
                    const float samples_per_angstrom = 6.0f;
                    size_t nto_idx = (size_t)rsp.selected;
                    for (int i = 0; i < num_lambdas; ++i) {
                        int pi = i * num_lambdas + 0;
                        int hi = i * num_lambdas + 1;
                        size_t lambda_idx = (size_t)i;

                        if (task_system::task_is_running(nto.vol_task[pi])) {
                            task_system::task_interrupt(nto.vol_task[pi]);
                        }
                        if (task_system::task_is_running(nto.vol_task[hi])) {
                            task_system::task_interrupt(nto.vol_task[hi]);
                        }

                        nto.vol_task[pi] = compute_nto_async(&nto.vol[pi].tex_to_world, &nto.vol[pi].step_size, &nto.vol[pi].tex_id, nto_idx, lambda_idx, MD_VLX_NTO_TYPE_PARTICLE, MD_GTO_EVAL_MODE_PSI, samples_per_angstrom);
                        nto.vol_task[hi] = compute_nto_async(&nto.vol[hi].tex_to_world, &nto.vol[hi].step_size, &nto.vol[hi].tex_id, nto_idx, lambda_idx, MD_VLX_NTO_TYPE_HOLE, MD_GTO_EVAL_MODE_PSI, samples_per_angstrom);
                    }
                    if (task_system::task_is_running(nto.seg_task[0])) {
                        task_system::task_interrupt(nto.seg_task[0]);
                    }
                    if (task_system::task_is_running(nto.seg_task[1])) {
                        task_system::task_interrupt(nto.seg_task[1]);
                    }

                    recompute_transition_matrix = true;
                }

                if (recompute_transition_matrix) {
                    const float samples_per_angstrom = 6.0f;
                    const size_t nto_idx = (size_t)rsp.selected;

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

					MEMSET(nto.transition_matrix, 0, sizeof(float) * nto.transition_matrix_dim * (nto.transition_matrix_dim + 2));

                    task_system::ID eval_part = 0;
                    task_system::ID seg_part  = 0;
                    task_system::ID eval_hole = 0;
                    task_system::ID seg_hole  = 0;

                    if (compute_nto_group_values_async(&eval_part, &seg_part, nto.transition_density_part, nto.group.count, nto.atom_group_idx, nto.atom_xyz, nto.atom_r, nto.num_atoms, nto_idx, 0, MD_VLX_NTO_TYPE_PARTICLE, MD_GTO_EVAL_MODE_PSI_SQUARED, samples_per_angstrom) &&
                        compute_nto_group_values_async(&eval_hole, &seg_hole, nto.transition_density_hole, nto.group.count, nto.atom_group_idx, nto.atom_xyz, nto.atom_r, nto.num_atoms, nto_idx, 0, MD_VLX_NTO_TYPE_HOLE,     MD_GTO_EVAL_MODE_PSI_SQUARED, samples_per_angstrom))
                    {
                        task_system::ID compute_matrix_task = task_system::create_main_task(STR_LIT("##Compute Transition Matrix"), [](void* user_data) {
                            VeloxChem::Nto* nto = (VeloxChem::Nto*)user_data;
                            distribute_charges_heuristic(nto->transition_matrix, nto->group.count, nto->transition_density_hole, nto->transition_density_part);
                        }, &nto);

                        task_system::set_task_dependency(compute_matrix_task, seg_part);
                        task_system::set_task_dependency(compute_matrix_task, seg_hole);

                        task_system::enqueue_task(eval_part);
                        task_system::enqueue_task(eval_hole);

                        nto.seg_task[0] = seg_part;
                        nto.seg_task[1] = seg_hole;
                    } else {
                        MD_LOG_DEBUG("An error occured when computing nto group values");
                    }
                }
            }

            const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();

            ImVec2 grid_p0 = canvas_p0;
            ImVec2 grid_p1 = canvas_p0 + canvas_sz * ImVec2(0.5f, 1.0f);
            ImVec2 win_sz = (grid_p1 - grid_p0) / ImVec2(1.0f, (float)(num_lambdas * 2));
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

                // Draw P / H orbitals
                for (int i = 0; i < num_lambdas * 2; ++i) {
                    ImVec2 p0 = grid_p0 + win_sz * ImVec2(0.0f, (float)(i+0));
                    ImVec2 p1 = grid_p0 + win_sz * ImVec2(1.0f, (float)(i+1));
                    ImRect rect = {p0, p1};
                    ImGui::SetCursorScreenPos(p0);

                    draw_list->ChannelsSetCurrent(1);
                    ImGui::PushID(i);
                    interaction_canvas(p1-p0, selection, view, state.mold.mol);

                    if (ImGui::IsItemHovered()) {
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
                    const char* lbl = ((i & 1) == 0) ? "Particle" : "Hole";

                    char lbls[64];
                    snprintf(lbls, sizeof(lbls), "hovered_atom_idx: %i, hovered_bond_idx: %i", state.selection.atom_idx.hovered, state.selection.bond_idx.hovered);

                    char buf[32];
                    snprintf(buf, sizeof(buf), (const char*)u8"λ: %.3f", nto_lambda[i / 2]);
                    draw_list->AddImage((ImTextureID)(intptr_t)nto.gbuf.tex.transparency, p0, p1, { 0,1 }, { 1,0 });
                    draw_list->AddImage((ImTextureID)(intptr_t)nto.iso_tex[i], p0, p1, { 0,1 }, { 1,0 });
                    draw_list->AddText(text_pos_bl, ImColor(0,0,0), buf);
                    draw_list->AddText(text_pos_tl, ImColor(0,0,0), lbls);
                    
                    const float aspect_ratio1 = win_sz.x / win_sz.y;
                    mat4_t view_mat1 = camera_world_to_view_matrix(nto.camera);
                    mat4_t proj_mat1 = camera_perspective_projection_matrix(nto.camera, aspect_ratio1);
                    int number_of_atoms = vlx.geom.num_atoms;
                    float middle_x=0;
                    float middle_y = 0;
                    float middle_z = 0;
                    for (int i = 0; i < number_of_atoms; ++i) {
                        middle_x = middle_x + (float)vlx.geom.coord_x[i];
                        middle_y = middle_y + (float)vlx.geom.coord_y[i];
                        middle_z = middle_z + (float)vlx.geom.coord_z[i];
                    }

                    middle_x=middle_x/number_of_atoms;
                    middle_y=middle_y/number_of_atoms;
                    middle_z=middle_z/number_of_atoms;
                    float farthest_x = 0;
                    float farthest_y = 0;
                    float farthest_z = 0;
                    float current_distance_x = 0;
                    float current_distance_y = 0;
                    float current_distance_z = 0;
                    for (int i = 0; i < number_of_atoms; ++i) {
                        current_distance_x = middle_x - (float)vlx.geom.coord_x[i];
                        current_distance_y = middle_y - (float)vlx.geom.coord_y[i];
                        current_distance_z = middle_z - (float)vlx.geom.coord_z[i];
                        if (vec3_length({current_distance_x, current_distance_y, current_distance_z}) > vec3_length({farthest_x, farthest_y, farthest_z})) {
                            farthest_x = current_distance_x;
                            farthest_y = current_distance_y;
                            farthest_z = current_distance_z;
                        }
                    }
                    float farthest_distance = vec3_length({farthest_x, farthest_y, farthest_z});
                    float longest_vector = 0;
                    if (vec3_length({(float)vlx.rsp.electronic_transition_length[rsp.selected].x, (float)vlx.rsp.electronic_transition_length[rsp.selected].y, (float)vlx.rsp.electronic_transition_length[rsp.selected].z}) > vec3_length({(float)vlx.rsp.magnetic_transition[rsp.selected].x, (float)vlx.rsp.magnetic_transition[rsp.selected].y, (float)vlx.rsp.magnetic_transition[rsp.selected].z}))
                    {
                        longest_vector = vec3_length({(float)vlx.rsp.electronic_transition_length[rsp.selected].x,
                                                      (float)vlx.rsp.electronic_transition_length[rsp.selected].y,
                                                      (float)vlx.rsp.electronic_transition_length[rsp.selected].z});
                    }
                    else
                    {
                        longest_vector = vec3_length({(float)vlx.rsp.magnetic_transition[rsp.selected].x,
                                                      (float)vlx.rsp.magnetic_transition[rsp.selected].y,
                                                      (float)vlx.rsp.magnetic_transition[rsp.selected].z});
                    }
                    const mat4_t mvp = proj_mat1 * view_mat1;
                    const vec4_t p = mat4_mul_vec4(mvp, {middle_x, middle_y, middle_z, 1.0});
                    const vec4_t p_electric_target =                
                        mat4_mul_vec4(mvp, {middle_x + (float)vlx.rsp.electronic_transition_length[rsp.selected].x * (float)nto.iso.vector_length*farthest_distance/longest_vector,
                                            middle_y + (float)vlx.rsp.electronic_transition_length[rsp.selected].y * (float)nto.iso.vector_length*farthest_distance/longest_vector,
                                            middle_z + (float)vlx.rsp.electronic_transition_length[rsp.selected].z * (float)nto.iso.vector_length*farthest_distance/longest_vector, 1.0f});
                    const vec2_t c = {
                        (p.x / p.w * 0.5f + 0.5f) * win_sz.x,
                        (-p.y / p.w * 0.5f + 0.5f) * win_sz.y,
                    };
                    const vec2_t c_electric_target = {
                        (p_electric_target.x / p_electric_target.w * 0.5f + 0.5f) * win_sz.x,
                        (-p_electric_target.y / p_electric_target.w * 0.5f + 0.5f) * win_sz.y,
                    };
                    float angle_electric = atan2f(c.y - c_electric_target.y, c.x - c_electric_target.x);

                    const vec2_t electric_triangle_point1 = {
                        p0.x + c_electric_target.x + cosf(angle_electric + 0.523599f) * 5,
                        p0.y + c_electric_target.y + sinf(angle_electric + 0.523599f) * 5,
                    };
                    const vec2_t electric_triangle_point2 = {
                        p0.x + c_electric_target.x + cosf(angle_electric - 0.523599f) * 5,
                        p0.y + c_electric_target.y + sinf(angle_electric - 0.523599f) * 5,
                    };
                    if (nto.iso.display_vectors) {
                            draw_list->AddLine({p0.x + c.x, p0.y + c.y}, {p0.x + c_electric_target.x, p0.y + c_electric_target.y},
                                           ImColor(nto.iso.vectorColors[0].elem[0], nto.iso.vectorColors[0].elem[1], nto.iso.vectorColors[0].elem[2],
                                                    nto.iso.vectorColors[0].elem[3]),
                                           5.0f);
                            draw_list->AddTriangle({p0.x + c_electric_target.x, p0.y + c_electric_target.y},
                                                   {electric_triangle_point1.x, electric_triangle_point1.y},
                                                   {electric_triangle_point2.x, electric_triangle_point2.y}, ImColor(nto.iso.vectorColors[0].elem[0], nto.iso.vectorColors[0].elem[1],nto.iso.vectorColors[0].elem[2],nto.iso.vectorColors[0].elem[3]), 5.0f);
                            draw_list->AddText({p0.x + c_electric_target.x, p0.y + c_electric_target.y},
                                               ImColor(nto.iso.vectorColors[0].elem[0], nto.iso.vectorColors[0].elem[1],
                                                       nto.iso.vectorColors[0].elem[2], nto.iso.vectorColors[0].elem[3]),
                                               (const char*)u8"μe");
                        }
                    const vec4_t p_magnetic_target =
                        mat4_mul_vec4(mvp, {middle_x + (float)vlx.rsp.magnetic_transition[rsp.selected].x * (float)nto.iso.vector_length*farthest_distance/longest_vector,
                                            middle_y + (float)vlx.rsp.magnetic_transition[rsp.selected].y * (float)nto.iso.vector_length*farthest_distance/longest_vector,
                                            middle_z + (float)vlx.rsp.magnetic_transition[rsp.selected].z * (float)nto.iso.vector_length*farthest_distance/longest_vector, 1.0f});
                    const vec2_t c_magnetic_target = {
                        (p_magnetic_target.x / p_magnetic_target.w * 0.5f + 0.5f) * win_sz.x,
                        (-p_magnetic_target.y / p_magnetic_target.w * 0.5f + 0.5f) * win_sz.y,
                    };
                    float angle_magnetic = atan2f(c.y - c_magnetic_target.y, c.x - c_magnetic_target.x);
                    const vec2_t magnetic_triangle_point1 = {
                        p0.x + c_magnetic_target.x + cosf(angle_magnetic + 0.523599f) * 5.0f,
                        p0.y + c_magnetic_target.y + sinf(angle_magnetic + 0.523599f) * 5.0f,
                    };
                    const vec2_t magnetic_triangle_point2 = {
                        p0.x + c_magnetic_target.x + cosf(angle_magnetic - 0.523599f) * 5.0f,
                        p0.y + c_magnetic_target.y + sinf(angle_magnetic - 0.523599f) * 5.0f,
                    };
                    vec3_t magnetic_vector_3d = {(float)vlx.rsp.magnetic_transition[rsp.selected].x, (float)vlx.rsp.magnetic_transition[rsp.selected].y, (float)vlx.rsp.magnetic_transition[rsp.selected].z};
                    vec3_t electronic_vector_3d = {(float)vlx.rsp.electronic_transition_length[rsp.selected].x, (float)vlx.rsp.electronic_transition_length[rsp.selected].y, (float)vlx.rsp.electronic_transition_length[rsp.selected].z};

                    float el_ma_dot = vec3_dot(magnetic_vector_3d, electronic_vector_3d);
                    float el_len = vec3_length(electronic_vector_3d);
                    float ma_len = vec3_length(magnetic_vector_3d);
                    float angle = acosf(el_ma_dot / (el_len * ma_len)) * (180.0 / 3.141592653589793238463);
                    char bufDPM[32];
                    const vec2_t magnetic_vector_2d = c_magnetic_target - c;
                    const vec2_t electric_vector_2d = c_electric_target - c;
                    if (nto.iso.display_vectors) {
                        draw_list->AddLine({p0.x + c.x, p0.y + c.y}, {p0.x + c_magnetic_target.x, p0.y + c_magnetic_target.y},
                                           ImColor(nto.iso.vectorColors[1].elem[0], nto.iso.vectorColors[1].elem[1], nto.iso.vectorColors[1].elem[2],
                                                   nto.iso.vectorColors[1].elem[3]),
                                           5.0f);
                        draw_list->AddTriangle({p0.x + c_magnetic_target.x, p0.y + c_magnetic_target.y},
                                               {magnetic_triangle_point1.x, magnetic_triangle_point1.y},
                                               {magnetic_triangle_point2.x, magnetic_triangle_point2.y}, ImColor(nto.iso.vectorColors[1].elem[0], nto.iso.vectorColors[1].elem[1],
                                                       nto.iso.vectorColors[1].elem[2], nto.iso.vectorColors[1].elem[3]), 5.0f);
                        draw_list->AddText({p0.x + c_magnetic_target.x, p0.y + c_magnetic_target.y},
                                           ImColor(nto.iso.vectorColors[1].elem[0], nto.iso.vectorColors[1].elem[1], nto.iso.vectorColors[1].elem[2],
                                                   nto.iso.vectorColors[1].elem[3]),
                                           (const char*)u8"μm");
                    }
                    if (nto.iso.display_angle) {
                        int num_of_seg = 10;
                        for (int segment = 0; segment < num_of_seg; segment++) {
                            vec3_t middle_point_3d0 = vec3_normalize(vec3_normalize(magnetic_vector_3d) * (num_of_seg - segment) / num_of_seg * 0.1f +
                                                                    vec3_normalize(electronic_vector_3d)*segment/num_of_seg * 0.1f) *0.03f;
                            const vec4_t p_middle_point_3d_target0 =
                                mat4_mul_vec4(mvp, {middle_x + middle_point_3d0.x * (float)nto.iso.vector_length*farthest_distance/longest_vector,
                                                    middle_y + middle_point_3d0.y * (float)nto.iso.vector_length*farthest_distance/longest_vector,
                                                    middle_z + middle_point_3d0.z * (float)nto.iso.vector_length*farthest_distance/longest_vector, 1.0f});
                            const vec2_t c_middle_point_target0 = {
                                (p_middle_point_3d_target0.x / p_middle_point_3d_target0.w * 0.5f + 0.5f) * win_sz.x,
                                (-p_middle_point_3d_target0.y / p_middle_point_3d_target0.w * 0.5f + 0.5f) * win_sz.y,
                            };
                            vec3_t middle_point_3d1 =
                                vec3_normalize(vec3_normalize(magnetic_vector_3d) * (num_of_seg - segment - 1) / num_of_seg * 0.1f +
                                                                    vec3_normalize(electronic_vector_3d)*(segment + 1 ) / num_of_seg * 0.1f) *0.03f;
                            const vec4_t p_middle_point_3d_target1 =
                                mat4_mul_vec4(mvp, {middle_x + middle_point_3d1.x * (float)nto.iso.vector_length*farthest_distance/longest_vector,
                                                    middle_y + middle_point_3d1.y * (float)nto.iso.vector_length*farthest_distance/longest_vector,
                                                    middle_z + middle_point_3d1.z * (float)nto.iso.vector_length*farthest_distance/longest_vector, 1.0f});
                            const vec2_t c_middle_point_target1 = {
                                (p_middle_point_3d_target1.x / p_middle_point_3d_target1.w * 0.5f + 0.5f) * win_sz.x,
                                (-p_middle_point_3d_target1.y / p_middle_point_3d_target1.w * 0.5f + 0.5f) * win_sz.y,
                            };
                            draw_list->AddLine({p0.x + c_middle_point_target0.x, p0.y + c_middle_point_target0.y},
                                               {p0.x + c_middle_point_target1.x, p0.y + c_middle_point_target1.y},
                                               ImColor(nto.iso.vectorColors[2].elem[0], nto.iso.vectorColors[2].elem[1],
                                                       nto.iso.vectorColors[2].elem[2], nto.iso.vectorColors[2].elem[3]),
                                               5.0f);
                        }

                        
                        snprintf(bufDPM, sizeof(bufDPM), (const char*)u8"θ=%.2f°", angle);
                        draw_list->AddText({p0.x + c.x, p0.y + c.y},
                                           ImColor(nto.iso.vectorColors[2].elem[0], nto.iso.vectorColors[2].elem[1], nto.iso.vectorColors[2].elem[2],
                                                   nto.iso.vectorColors[2].elem[3]),
                                           bufDPM);
                    }
                    

                }
                // @TODO: Draw Sankey Diagram of Transition Matrix
                {
                    ImVec2 p0 = canvas_p0 + canvas_sz * ImVec2(0.5f, 0.0f);
                    ImVec2 p1 = canvas_p1;
                    im_sankey_diagram({p0.x, p0.y, p1.x, p1.y}, &nto);
                    ImVec2 text_pos_bl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p1.y - TEXT_BASE_HEIGHT);
                    draw_list->AddText(text_pos_bl, ImColor(0, 0, 0, 255), "Transition Diagram");
                }
                // Draw grid
                {
                    ImVec2 p0 = {floorf(canvas_p0.x + canvas_sz.x * 0.5f), canvas_p0.y};
                    ImVec2 p1 = {floorf(canvas_p0.x + canvas_sz.x * 0.5f), canvas_p1.y};
                    draw_list->AddLine(p0, p1, IM_COL32(0, 0, 0, 255));
                }
                for (int i = 1; i < num_lambdas * 2; ++i) {
                    float y = floorf(canvas_p0.y + canvas_sz.y / ((float)num_lambdas * 2.0f) * i);
                    float x0 = canvas_p0.x;
                    float x1 = floorf(canvas_p0.x + canvas_sz.x * (i & 1 ? 0.5f : 1.0f));
                    draw_list->AddLine({x0, y}, {x1, y}, IM_COL32(0, 0, 0, 255));
                }

                // Draw stuff


                const bool is_hovered = ImGui::IsItemHovered();
                const bool is_active = ImGui::IsItemActive();
                const ImVec2 origin(canvas_p0.x, canvas_p0.y);  // Lock scrolled origin
                const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

                int width  = MAX(1, (int)win_sz.x);
                int height = MAX(1, (int)win_sz.y);

                int num_win = num_lambdas * 2;

                auto& gbuf = nto.gbuf;
                if ((int)gbuf.width != width || (int)gbuf.height != height) {
                    init_gbuffer(&gbuf, width, height);
                    for (int i = 0; i < num_win; ++i) {
                        gl::init_texture_2D(nto.iso_tex + i, width, height, GL_RGBA8);
                    }
                }

                if (is_hovered) {
                    if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                        reset_view = true;
                    }
                }

                if (reset_view) {
                    camera_compute_optimal_view(&nto.target.pos, &nto.target.ori, &nto.target.dist, min_aabb, max_aabb, nto.distance_scale);
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

                    if (nto.iso.enabled) {
                        PUSH_GPU_SECTION("NTO RAYCAST")
                            for (int i = 0; i < num_win; ++i) {
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
                                        .model = nto.vol[i].tex_to_world,
                                        .view  = view_mat,
                                        .proj  = proj_mat,
                                        .inv_proj = inv_proj_mat,
                                },
                                .iso = {
                                        .enabled = true,
                                        .count  = (size_t)nto.iso.count,
                                        .values = nto.iso.values,
                                        .colors = nto.iso.colors,
                                },
                                .shading = {
                                        .env_radiance = state.visuals.background.color * state.visuals.background.intensity * 0.25f,
                                        .roughness = 0.3f,
                                        .dir_radiance = {10,10,10},
                                        .ior = 1.5f,
                                },
                                .voxel_spacing = nto.vol[i].step_size,
                                };
                                volume::render_volume(vol_desc);
                            }
                        POP_GPU_SECTION();
                    }

                // Reset state
                glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
                glDrawBuffer(GL_BACK);
                glDisable(GL_SCISSOR_TEST);

                if (hovered_canvas_rect.GetArea() > 0) {
                    ImVec2 coord = ImGui::GetMousePos() - hovered_canvas_rect.Min;
                    coord.y = hovered_canvas_rect.GetSize().y - coord.y;
                    update_picking_data(nto.picking, {coord.x, coord.y}, nto.gbuf, inv_MVP);

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
        }
        ImGui::End();

        if (open_context_menu) {
            ImGui::OpenPopup("Context Menu");
        }

        if (ImGui::BeginPopup("Context Menu")) {
            if (ImGui::MenuItem("Example Button")) {

            }
            ImGui::EndPopup();
        }
    }
};
static VeloxChem instance = {};
