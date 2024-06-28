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

#define BLK_DIM 8
#define ANGSTROM_TO_BOHR 1.8897261246257702
#define BOHR_TO_ANGSTROM 0.529177210903

#define IM_GREEN ImVec4{0, 1, 0, 1}
#define IM_RED ImVec4{1, 0, 0, 1}

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

struct VeloxChem : viamd::EventHandler {
    VeloxChem() { viamd::event_system_register_handler(*this); }

    md_vlx_data_t vlx {};

    // Used for clearing volumes
    uint32_t vol_fbo = 0;

    // GL representations
    md_gl_mol_t gl_mol = {};
    md_gl_rep_t gl_rep = {};

    int homo_idx = 0;
    int lumo_idx = 0;

    // Principal Component Axes of the geometry
    mat3_t PCA = mat3_ident();
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
        int vol_nto_idx = -1;

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
        int hovered = -1;
        int selected = -1;
        int focused_plot = -1;

        double* x_ev_samples;
        double* x_unit_samples;
        double* x_unit_peaks;
        double* eps;
        double* ecd;
        double* vib_y;
        double* vib_x;
        double* vib_points;
        double* osc_points;
        double* cgs_points;
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
                    task_system::ID id = compute_mo(&data.tex_mat, &data.voxel_spacing, data.dst_texture, data.orbital_idx, mode, data.samples_per_angstrom);
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
        md_gl_mol_destroy(gl_mol);
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

                // Compute the PCA of the provided geometry
                // This is used in determining a better fitting volume for the orbitals
                vec4_t* xyzw = (vec4_t*)md_vm_arena_push(state.allocator.frame, sizeof(vec4_t) * vlx.geom.num_atoms);
                for (size_t i = 0; i < vlx.geom.num_atoms; ++i) {
                    xyzw[i] = {(float)vlx.geom.coord_x[i], (float)vlx.geom.coord_y[i], (float)vlx.geom.coord_z[i], 1.0f};
                    min_box = vec4_min(min_box, xyzw[i]);
                    max_box = vec4_max(max_box, xyzw[i]);
                }
                min_aabb = vec3_from_vec4(min_box);
                max_aabb = vec3_from_vec4(max_box);

                md_molecule_t mol = {0};
                md_vlx_molecule_init(&mol, &vlx, state.allocator.frame);
                md_util_molecule_postprocess(&mol, state.allocator.frame, MD_UTIL_POSTPROCESS_ELEMENT_BIT | MD_UTIL_POSTPROCESS_RADIUS_BIT | MD_UTIL_POSTPROCESS_BOND_BIT);
                gl_mol = md_gl_mol_create(&mol);

                uint32_t* colors = (uint32_t*)md_vm_arena_push(state.allocator.frame, mol.atom.count * sizeof(uint32_t));
                color_atoms_cpk(colors, mol.atom.count, mol);

                gl_rep = md_gl_rep_create(gl_mol);
                md_gl_rep_set_color(gl_rep, 0, (uint32_t)mol.atom.count, colors, 0);

                vec3_t com =  md_util_com_compute_vec4(xyzw, 0, vlx.geom.num_atoms, 0);
                mat3_t C = mat3_covariance_matrix_vec4(xyzw, 0, vlx.geom.num_atoms, com);
                mat3_eigen_t eigen = mat3_eigen(C);
                PCA = mat3_extract_rotation(eigen.vectors);


                // NTO
                if (vlx.rsp.num_excited_states > 0 && vlx.rsp.nto) {
                    nto.show_window = true;
                    camera_compute_optimal_view(&nto.target.pos, &nto.target.ori, &nto.target.dist, min_aabb, max_aabb, nto.distance_scale);
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

    // Compute an 'optimal' OOBB for the supplied PGTOS
    void compute_volume_basis_and_transform_pgtos(mat4_t* out_tex_mat, vec3_t* out_ext_in_angstrom, md_gto_t* pgtos, size_t num_pgtos) {
        // Compute min and maximum extent along the PCA axes
        vec4_t min_ext = { FLT_MAX, FLT_MAX, FLT_MAX, 0};
        vec4_t max_ext = {-FLT_MAX,-FLT_MAX,-FLT_MAX, 0};

        mat4_t R = mat4_from_mat3(mat3_transpose(PCA));
        mat4_t Ri = mat4_from_mat3(PCA);

        for (size_t i = 0; i < num_pgtos; ++i) {
            vec4_t p = vec4_set(pgtos[i].x, pgtos[i].y, pgtos[i].z, 1.0f);
            p = Ri * (p * BOHR_TO_ANGSTROM);
            // The 0.9 scaling factor here is a bit arbitrary, but the cutoff-radius is computed on a value which is lower than the rendered iso-value
            // So the effective radius is a bit overestimated and thus we scale it back a bit
            min_ext = vec4_min(min_ext, vec4_sub_f(p, (float)(pgtos[i].cutoff * BOHR_TO_ANGSTROM * 0.9)));
            max_ext = vec4_max(max_ext, vec4_add_f(p, (float)(pgtos[i].cutoff * BOHR_TO_ANGSTROM * 0.9)));
        }

        min_ext.w = 0.0f;
        max_ext.w = 0.0f;

        vec4_t origin = R * min_ext;
        vec4_t extent = max_ext - min_ext;

        mat4_t T = mat4_translate(origin.x, origin.y, origin.z);
        mat4_t S = mat4_scale(extent.x, extent.y, extent.z);

        // Transform pgtos into the 'model space' of the volume
        // There is a scaling factor here as well as the PGTOs are given and must be
        // Processed in a.u. (Bohr) and not Ångström, in which the volume is defined
        
        vec4_t t = vec4_mul_f(origin, ANGSTROM_TO_BOHR);
        mat4_t M = Ri * mat4_translate(-t.x, -t.y, -t.z);
        for (size_t i = 0; i < num_pgtos; ++i) {
            vec4_t c = {pgtos[i].x, pgtos[i].y, pgtos[i].z, 1.0f};
            c = M * c;
            pgtos[i].x = c.x;
            pgtos[i].y = c.y;
            pgtos[i].z = c.z;
        }

        *out_tex_mat = T * R * S;
        *out_ext_in_angstrom = vec3_from_vec4(extent);
    }

    task_system::ID compute_nto(mat4_t* out_tex_mat, vec3_t* out_voxel_spacing, uint32_t* in_out_vol_tex, size_t nto_idx, size_t lambda_idx, md_vlx_nto_type_t type, md_gto_eval_mode_t mode, float samples_per_angstrom = 8.0f) {
        size_t num_pgtos = md_vlx_nto_pgto_count(&vlx);
        md_gto_t* pgtos  = (md_gto_t*)md_alloc(md_get_heap_allocator(), sizeof(md_gto_t) * num_pgtos);

        if (!md_vlx_nto_pgto_extract(pgtos, &vlx, nto_idx, lambda_idx, type)) {
            MD_LOG_ERROR("Failed to extract NTO pgtos for nto index: %zu and lambda: %zu", nto_idx, lambda_idx);
            md_free(md_get_heap_allocator(), pgtos, sizeof(md_gto_t) * num_pgtos);
            return task_system::INVALID_ID;
        }
        md_gto_cutoff_compute(pgtos, num_pgtos, 1.0e-6);

        mat4_t tex_mat;
        vec3_t extent;

        compute_volume_basis_and_transform_pgtos(&tex_mat, &extent, pgtos, num_pgtos);

        // Target resolution per spatial unit (We want this number of samples per Ångström in each dimension)
        // Round up to some multiple of 8 -> 8x8x8 is the chunksize we process in parallel

        // Compute required volume dimensions
        int dim[3] = {
            CLAMP(ALIGN_TO((int)(extent.x * samples_per_angstrom), 8), 8, 512),
            CLAMP(ALIGN_TO((int)(extent.y * samples_per_angstrom), 8), 8, 512),
            CLAMP(ALIGN_TO((int)(extent.z * samples_per_angstrom), 8), 8, 512),
        };

        vec3_t step_size = vec3_div(extent, vec3_set((float)dim[0], (float)dim[1], (float)dim[2]));
        step_size = vec3_mul_f(step_size, (float)ANGSTROM_TO_BOHR);

        // Create texture of dim
        gl::init_texture_3D(in_out_vol_tex, dim[0], dim[1], dim[2], GL_R16F);

        // WRITE OUTPUT
        *out_tex_mat = tex_mat;
        *out_voxel_spacing = step_size;

        return async_evaluate_orbital_on_grid(in_out_vol_tex, step_size.elem, dim, pgtos, num_pgtos, mode);
    }

    task_system::ID compute_mo(mat4_t* out_tex_mat, vec3_t* out_voxel_spacing, uint32_t* in_out_vol_tex, size_t mo_idx, md_gto_eval_mode_t mode, float samples_per_angstrom = 8.0f) {
        size_t num_pgtos = md_vlx_mol_pgto_count(&vlx);
        md_gto_t* pgtos  = (md_gto_t*)md_alloc(md_get_heap_allocator(), sizeof(md_gto_t) * num_pgtos);

        if (!md_vlx_mol_pgto_extract(pgtos, &vlx, mo_idx)) {
            MD_LOG_ERROR("Failed to extract molecular pgtos for orbital index: %zu", mo_idx);
            md_free(md_get_heap_allocator(), pgtos, sizeof(md_gto_t) * num_pgtos);
            return task_system::INVALID_ID;
        }
        md_gto_cutoff_compute(pgtos, num_pgtos, 1.0e-6);

        mat4_t tex_mat;
        vec3_t extent;

        compute_volume_basis_and_transform_pgtos(&tex_mat, &extent, pgtos, num_pgtos);

        // Target resolution per spatial unit (We want this number of samples per Ångström in each dimension)
        // Round up to some multiple of 8 -> 8x8x8 is the chunksize we process in parallel

        // Compute required volume dimensions
        int dim[3] = {
            CLAMP(ALIGN_TO((int)(extent.x * samples_per_angstrom), 8), 8, 512),
            CLAMP(ALIGN_TO((int)(extent.y * samples_per_angstrom), 8), 8, 512),
            CLAMP(ALIGN_TO((int)(extent.z * samples_per_angstrom), 8), 8, 512),
        };

        vec3_t step_size = vec3_div(extent, vec3_set((float)dim[0], (float)dim[1], (float)dim[2]));
        step_size = vec3_mul_f(step_size, (float)ANGSTROM_TO_BOHR);

        // Init and clear volume texture
        const float zero[4] = {};
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, vol_fbo);
        gl::init_texture_3D(in_out_vol_tex, dim[0], dim[1], dim[2], GL_R16F);
        glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, *in_out_vol_tex, 0);
        glClearBufferfv(GL_COLOR, 0, zero);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

        // WRITE OUTPUT
        *out_tex_mat = tex_mat;
        *out_voxel_spacing = step_size;

        return async_evaluate_orbital_on_grid(in_out_vol_tex, step_size.elem, dim, pgtos, num_pgtos, mode);
    }

    // This is a bit quirky, this will take ownership of pgtos and will free them after the evaluation is complete
    task_system::ID async_evaluate_orbital_on_grid(uint32_t* tex_ptr, const float step_size[3], const int dim[3], md_gto_t* pgtos, size_t num_pgtos, md_gto_eval_mode_t mode) {
        struct Payload {
            size_t bytes;
            size_t num_pgtos;
            md_gto_t* pgtos;
            float* vol_data;
            int    vol_dim[3];
            float  step_size[3];
            uint32_t* tex_ptr;
            md_gto_eval_mode_t mode;
        };

        size_t num_vol_bytes = dim[0] * dim[1] * dim[2] * sizeof(float);
        size_t num_bytes = sizeof(Payload) + num_vol_bytes;
        void* mem = md_alloc(md_get_heap_allocator(), num_bytes);
        Payload* payload    = (Payload*)mem;
        payload->bytes      = num_bytes;
        payload->num_pgtos  = num_pgtos;
        payload->pgtos      = pgtos;
        payload->vol_data   = (float*)((char*)mem + sizeof(Payload));
        MEMSET(payload->vol_data, 0, num_vol_bytes);
        MEMCPY(payload->vol_dim, dim, sizeof(payload->vol_dim));
        MEMCPY(payload->step_size, step_size, sizeof(payload->step_size));
        payload->tex_ptr    = tex_ptr;
        payload->mode       = mode;

        // We evaluate the in parallel over smaller NxNxN blocks
        uint32_t num_blocks = (dim[0] / BLK_DIM) * (dim[1] / BLK_DIM) * (dim[2] / BLK_DIM);

        MD_LOG_DEBUG("Starting async eval of orbital grid [%i][%i][%i]", dim[0], dim[1], dim[2]);

        task_system::ID async_task = task_system::create_pool_task(STR_LIT("Evaluate Orbital"), 0, num_blocks, [](uint32_t range_beg, uint32_t range_end, void* user_data, uint32_t thread_num) {
            (void)thread_num;
            Payload* data = (Payload*)user_data;

            // Number of NxNxN blocks in each dimension
            int num_blk[3] = {
                data->vol_dim[0] / BLK_DIM,
                data->vol_dim[1] / BLK_DIM,
                data->vol_dim[2] / BLK_DIM,
            };

            const float* step = data->step_size;

            md_grid_t grid = {
                .data = data->vol_data,
                .dim  = {data->vol_dim[0], data->vol_dim[1], data->vol_dim[2]},
                // Shift origin by half a voxel to evaluate at the voxel center
                .origin = {0.5f * step[0], 0.5f * step[1], 0.5f * step[2]},
                .stepsize = {step[0], step[1], step[2]},
            };

            size_t temp_pos = md_temp_get_pos();
            md_gto_t* sub_pgtos = (md_gto_t*)md_temp_push(sizeof(md_gto_t) * data->num_pgtos);

            for (uint32_t i = range_beg; i < range_end; ++i) {
                // Determine block index from linear input index i
                int blk_x =  i % num_blk[0];
                int blk_y = (i / num_blk[0]) % num_blk[1];
                int blk_z =  i / (num_blk[0] * num_blk[1]);

                int off_idx[3] = {blk_x * BLK_DIM, blk_y * BLK_DIM, blk_z * BLK_DIM};
                int len_idx[3] = {BLK_DIM, BLK_DIM, BLK_DIM};

                float aabb_min[3] = {
                    grid.origin[0] + off_idx[0] * grid.stepsize[0],
                    grid.origin[1] + off_idx[1] * grid.stepsize[1],
                    grid.origin[2] + off_idx[2] * grid.stepsize[2],
                };
                float aabb_max[3] = {
                    grid.origin[0] + (off_idx[0] + len_idx[0]) * grid.stepsize[0],
                    grid.origin[1] + (off_idx[1] + len_idx[1]) * grid.stepsize[1],
                    grid.origin[2] + (off_idx[2] + len_idx[2]) * grid.stepsize[2],
                };

                size_t num_sub_pgtos = md_gto_aabb_test(sub_pgtos, aabb_min, aabb_max, data->pgtos, data->num_pgtos);
                md_gto_grid_evaluate_sub(&grid, off_idx, len_idx, sub_pgtos, num_sub_pgtos, data->mode);
            }

            md_temp_set_pos_back(temp_pos);
        }, payload);

        // Launch task for main (render) thread to update the volume texture
        task_system::ID main_task = task_system::create_main_task(STR_LIT("##Update Volume"), [](void* user_data) {
            Payload* data = (Payload*)user_data;
            
            // The init here is just to ensure that the volume has not changed its dimensions during the async evaluation
            gl::init_texture_3D(data->tex_ptr, data->vol_dim[0], data->vol_dim[1], data->vol_dim[2], GL_R16F);
            gl::set_texture_3D_data(*data->tex_ptr, data->vol_data, GL_R32F);

            md_free(md_get_heap_allocator(), data->pgtos, data->num_pgtos * sizeof(md_gto_t));
            md_free(md_get_heap_allocator(), data, data->bytes);
        }, payload);

        task_system::set_task_dependency(main_task, async_task);
        task_system::enqueue_task(async_task);

        return async_task;
    }

    static inline double axis_conversion_multiplier(const double* y1_array, const double* y2_array, size_t y1_array_size, size_t y2_array_size) {
        double y1_min = 0;
        double y1_max = 0;
        double y2_min = 0;
        double y2_max = 0;
        for (size_t i = 0; i < y1_array_size; i++) {
            y1_min = MIN(y1_min, y1_array[i]);
            y1_max = MAX(y1_max, y1_array[i]);
        }
        for (size_t i = 0; i < y2_array_size; i++) {
            y2_min = MIN(y2_min, y2_array[i]);
            y2_max = MAX(y2_max, y2_array[i]);
        }

        double y1_dist = fabs(y1_max - y1_min);
        double y2_dist = fabs(y2_max - y2_min);

        return y2_dist / y1_dist;
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
        ImPlotRect lim = { x_samples[0],x_samples[num_samples - 1],0,0};
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
                    static ImGuiTableFlags flags =  ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollX |
                                                    ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit;

                    static ImGuiTableColumnFlags columns_base_flags = ImGuiTableColumnFlags_NoSort;

                    if (ImGui::BeginTable("table_advanced", 5, flags, ImVec2(500, -1), 0)) {
                        ImGui::TableSetupColumn("Atom", columns_base_flags, 0.0f);
                        ImGui::TableSetupColumn("Symbol", columns_base_flags, 0.0f);
                        ImGui::TableSetupColumn("Coord X", columns_base_flags, 0.0f);
                        ImGui::TableSetupColumn("Coord Y", columns_base_flags, 0.0f);
                        ImGui::TableSetupColumn("Coord Z", columns_base_flags | ImGuiTableColumnFlags_WidthFixed, 0.0f);
                        ImGui::TableSetupScrollFreeze(0, 1);
                        ImGui::TableHeadersRow();

                        ImGui::PushStyleColor(ImGuiCol_HeaderHovered, ImVec4(1, 1, 0.5, 0.3));
                        ImGui::PushStyleColor(ImGuiCol_Header, ImVec4(0.5, 0.5, 1, 0.3));
                        bool item_hovered = false;
                        for (int row_n = 0; row_n < vlx.geom.num_atoms; row_n++) {

                            ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                            bool is_sel = md_bitfield_test_bit(&state.selection.selection_mask, row_n); //If atom is selected, mark it as such
                            bool is_hov = md_bitfield_test_bit(&state.selection.highlight_mask, row_n); //If atom is hovered, mark it as such
                            bool hov_col = false;
                            ImGui::TableNextRow(ImGuiTableRowFlags_None, 0);
                            ImGui::TableNextColumn();

                            if (is_hov) {
                                ImGui::PushStyleColor(ImGuiCol_Header, ImVec4(1, 1, 0.5, 0.3));
                            }
                            else {
                                ImGui::PushStyleColor(ImGuiCol_Header, ImVec4(0.5, 0.5, 1, 0.3));
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
        if (ImGui::Begin("Spectra", &rsp.show_window)) {
            if (ImGui::TreeNode("Absorption & ECD")) {
                bool refit1 = false;
                static bool first_plot1 = true;
                bool recalculate1 = false;

                static float gamma1 = 0.123;
                static x_unit_t x_unit = X_UNIT_EV;
                static broadening_mode_t broadening_mode1 = BROADENING_LORENTZIAN;
                const char* x_unit_str[] = {"Energy (eV)", "Wavelength (nm)", (const char*)u8"Wavenumber (cm⁻¹)", "Energy (hartree)"};

                recalculate1 = ImGui::SliderFloat((const char*)u8"Broadening γ HWHM (eV)", &gamma1, 0.01f, 1.0f);
                refit1 |= ImGui::Combo("Broadening mode", (int*)(&broadening_mode1), broadening_str, IM_ARRAYSIZE(broadening_str));
                refit1 |= ImGui::Combo("X unit", (int*)(&x_unit), x_unit_str, IM_ARRAYSIZE(x_unit_str));

                const int num_peaks = (int)vlx.rsp.num_excited_states;
                const double* y_osc_peaks = vlx.rsp.absorption_osc_str;
                const double* y_cgs_peaks = vlx.rsp.electronic_circular_dichroism_cgs;

                if (first_plot1) {
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

                if (recalculate1 || first_plot1) {
                    osc_to_eps(rsp.eps, rsp.x_ev_samples, num_samples, y_osc_peaks, vlx.rsp.absorption_ev, num_peaks, distr_func, gamma1 * 2);
                    rot_to_eps_delta(rsp.ecd, rsp.x_ev_samples, num_samples, y_cgs_peaks, vlx.rsp.absorption_ev, num_peaks, distr_func, gamma1 * 2);
                }

                static ImPlotRect osc_lim_constraint = {0, 0, 0, 0};
                static ImPlotRect cgs_lim_constraint = {0, 0, 0, 0};
                if (refit1 || first_plot1) {
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

                if (first_plot1) {
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
                    if (recalculate1 || first_plot1) {
                        osc_to_eps_mult = is_all_zero(y_osc_peaks, num_peaks) ? 1 : axis_conversion_multiplier(y_osc_peaks, rsp.eps, num_peaks, num_samples);
                    }

                    static ImPlotRect cur_osc_lims = {0, 1, 0, 1};
                    if (refit1 || first_plot1) {
                        ImPlot::SetNextAxisToFit(ImAxis_X1);
                    }
                    if (ImPlot::BeginPlot("Absorption")) {
                        ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                        ImPlot::SetupAxis(ImAxis_X1, x_unit_str[x_unit]);
                        ImPlot::SetupAxis(ImAxis_Y1, "f", ImPlotAxisFlags_AuxDefault);
                        ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"ε (L mol⁻¹ cm⁻¹)");
                        if (refit1 || first_plot1) {
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
                            rsp.selected = rsp.hovered;
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
                    if (recalculate1 || first_plot1) {
                        cgs_to_ecd_mult = is_all_zero(y_cgs_peaks, num_peaks) ? 1 : axis_conversion_multiplier(y_cgs_peaks, rsp.ecd, num_peaks, num_samples);
                    }
                    static ImPlotRect cur_cgs_lims = {0, 1, 0, 1};
                    if (refit1 || first_plot1) {
                        ImPlot::SetNextAxisToFit(ImAxis_X1);
                    }

                    if (ImPlot::BeginPlot("ECD")) {
                        ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                        ImPlot::SetupAxis(ImAxis_X1, x_unit_str[x_unit]);
                        ImPlot::SetupAxis(ImAxis_Y1, (const char*)u8"R (10⁻⁴⁰ cgs)", ImPlotAxisFlags_AuxDefault);
                        ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"Δε(ω) (L mol⁻¹ cm⁻¹)");
                        if (refit1 || first_plot1) {
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

                        if (rsp.hovered != -1 && ImPlot::IsPlotHovered()) {
                            draw_bar(2, rsp.x_unit_peaks[rsp.hovered], y_cgs_peaks[rsp.hovered], bar_width, IM_GREEN);
                            ImPlot::DragPoint(0, &rsp.x_unit_peaks[rsp.hovered], &rsp.cgs_points[rsp.hovered], IM_GREEN, 4, ImPlotDragToolFlags_NoInputs);

                        }

                        // Update selected peak on click
                        if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Left) &&
                            ImPlot::IsPlotHovered()) {
                            rsp.selected = rsp.hovered;
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
                first_plot1 = false;
                ImGui::TreePop();
            }
            if (ImGui::TreeNode("Vibrational Analysis")) {
                bool refit2 = false;
                bool recalculate2 = false;
                static bool first_plot2 = true;
                static float gamma2 = 5.0f;
                static broadening_mode_t broadening_mode2 = BROADENING_LORENTZIAN;
                recalculate2 = ImGui::SliderFloat((const char*)u8"Broadening γ HWHM (cm⁻¹)", &gamma2, 1.0f, 10.0f);
                refit2 |= ImGui::Combo("Broadening mode", (int*)(&broadening_mode2), broadening_str, IM_ARRAYSIZE(broadening_str));


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

                ImVec2* pixel_peaks = (ImVec2*)md_temp_push(sizeof(ImVec2) * num_vibs);
                ImVec2* pixel_points = (ImVec2*)md_temp_push(sizeof(ImVec2) * num_vibs);

                int hov_vib = -1;
                static int sel_vib = -1;
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

                if (first_plot2) {
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

                if (first_plot2 || recalculate2 || refit2) {
                    general_broadening(rsp.vib_y, rsp.vib_x, num_samples, irs, har_freqs, num_vibs, distr_func, gamma2 * 2);
                }

                if (first_plot2) {
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
                if (refit2 || first_plot2) {
                    lim_constraint = get_plot_limits(rsp.vib_x, irs, num_vibs, num_samples);
                }

                if (ImPlot::BeginPlot("Vibrational analysis")) {
                    // @HACK: Compute pixel width of 2 'plot' units
                    ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                    ImPlot::SetupAxis(ImAxis_X1, (const char*)u8"Harmonic Frequency (cm⁻¹)", x_flag);
                    ImPlot::SetupAxis(ImAxis_Y1, "IR Intensity (km/mol)", y_flag);
                    if (refit2 || first_plot2) {
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
                        sel_vib = hov_vib;
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
                    first_plot2 = false;
                    ImPlot::EndPlot();
                }

                // ImGui::Text("%i is hovered", hov_vib);
                // ImGui::Text("%f is z coord", (float)state.mold.mol.atom.z[2]);

                ImGui::SliderFloat((const char*)"Amplitude", &amp_mult, 0.2f, 2.0f);
                ImGui::SliderFloat((const char*)"Speed", &speed_mult, 0.5f, 2.0f);

                ImGui::TreePop();
            }
        }
        ImGui::End();
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
                        orb.vol_task[slot_idx] = compute_mo(&orb.vol[slot_idx].tex_to_world, &orb.vol[slot_idx].step_size, &orb.vol[slot_idx].tex_id, mo_idx, MD_GTO_EVAL_MODE_PSI, samples_per_angstrom);
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

    void draw_nto_window(const ApplicationState& state) {
        if (!nto.show_window) return;
        if (vlx.rsp.num_excited_states == 0) return;
        if (vlx.rsp.nto == NULL) return;

        ImGui::SetNextWindowSize(ImVec2(500, 300), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("NTO viewer", &nto.show_window, ImGuiWindowFlags_MenuBar)) {

            const ImVec2 outer_size = {300.f, 0.f};
            ImGui::PushItemWidth(outer_size.x);
            ImGui::BeginGroup();

            ImGui::Text("States");
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
            ImGui::Spacing();
            const double iso_min = 1.0e-4;
            const double iso_max = 5.0;
            double iso_val = nto.iso.values[0];
            ImGui::Text("NTO");                
            ImGui::Spacing();
            ImGui::Text("Isovalue"); 
            ImGui::SliderScalar("##Iso Value", ImGuiDataType_Double, &iso_val, &iso_min, &iso_max, "%.6f", ImGuiSliderFlags_Logarithmic);
            ImGui::SetItemTooltip("Iso Value");

            nto.iso.values[0] = (float)iso_val;
            nto.iso.values[1] = -(float)iso_val;
            nto.iso.count = 2;
            nto.iso.enabled = true;
            ImGui::Text("Orbital Colors");
            ImGui::ColorEdit4("##Color Positive", nto.iso.colors[0].elem);
            ImGui::SetItemTooltip("Color Positive");
            ImGui::ColorEdit4("##Color Negative", nto.iso.colors[1].elem);
            ImGui::SetItemTooltip("Color Negative");
            ImGui::Spacing();
            ImGui::Spacing();

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

            // @TODO: Enlist all defined groups here
            if (ImGui::BeginListBox("##Groups", outer_size)) {
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

            // This will catch our interactions
            ImGui::InvisibleButton("canvas", canvas_sz, ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_AllowOverlap);

            // Draw border and background color
            ImGuiIO& io = ImGui::GetIO();

            ImVec2 canvas_p0 = ImGui::GetItemRectMin();
            ImVec2 canvas_p1 = ImGui::GetItemRectMax();

            double nto_lambda[4] = {};

            int num_lambdas = 1;

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

                if (nto.vol_nto_idx != rsp.selected) {
                    nto.vol_nto_idx  = rsp.selected;
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

                        nto.vol_task[pi] = compute_nto(&nto.vol[pi].tex_to_world, &nto.vol[pi].step_size, &nto.vol[pi].tex_id, nto_idx, lambda_idx, MD_VLX_NTO_TYPE_PARTICLE, MD_GTO_EVAL_MODE_PSI, samples_per_angstrom);
                        nto.vol_task[hi] = compute_nto(&nto.vol[hi].tex_to_world, &nto.vol[hi].step_size, &nto.vol[hi].tex_id, nto_idx, lambda_idx, MD_VLX_NTO_TYPE_HOLE,     MD_GTO_EVAL_MODE_PSI, samples_per_angstrom);
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
            draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(255, 255, 255, 255));
            if (rsp.selected != -1) {
                // Draw P / H orbitals
                for (int i = 0; i < num_lambdas * 2; ++i) {
                    ImVec2 p0 = grid_p0 + win_sz * ImVec2(0.0f, (float)(i+0));
                    ImVec2 p1 = grid_p0 + win_sz * ImVec2(1.0f, (float)(i+1));
                    ImVec2 text_pos_bl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p1.y - TEXT_BASE_HEIGHT);
                    ImVec2 text_pos_tl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p0.y + TEXT_BASE_HEIGHT * 0.5f);
                    const char* lbl = ((i & 1) == 0) ? "Particle" : "Hole";
                    char buf[32];
                    snprintf(buf, sizeof(buf), (const char*)u8"λ: %.3f", nto_lambda[i / 2]);
                    draw_list->AddImage((ImTextureID)(intptr_t)nto.gbuf.tex.transparency, p0, p1, { 0,1 }, { 1,0 });
                    draw_list->AddImage((ImTextureID)(intptr_t)nto.iso_tex[i], p0, p1, { 0,1 }, { 1,0 });
                    draw_list->AddText(text_pos_bl, ImColor(0,0,0), buf);
                    draw_list->AddText(text_pos_tl, ImColor(0,0,0), lbl);
                    
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
                    float angle_electric = atan2(c.y - c_electric_target.y, c.x - c_electric_target.x);

                    const vec2_t electric_triangle_point1 = {
                        p0.x + c_electric_target.x + cos(angle_electric + 0.523599) * 5,
                        p0.y + c_electric_target.y + sin(angle_electric + 0.523599) * 5,
                    };
                    const vec2_t electric_triangle_point2 = {
                        p0.x + c_electric_target.x + cos(angle_electric - 0.523599) * 5,
                        p0.y + c_electric_target.y + sin(angle_electric - 0.523599) * 5,
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
                    float angle_magnetic = atan2(c.y - c_magnetic_target.y, c.x - c_magnetic_target.x);
                    const vec2_t magnetic_triangle_point1 = {
                        p0.x + c_magnetic_target.x + cos(angle_magnetic + 0.523599) * 5,
                        p0.y + c_magnetic_target.y + sin(angle_magnetic + 0.523599) * 5,
                    };
                    const vec2_t magnetic_triangle_point2 = {
                        p0.x + c_magnetic_target.x + cos(angle_magnetic - 0.523599) * 5,
                        p0.y + c_magnetic_target.y + sin(angle_magnetic - 0.523599) * 5,
                    };
                    vec3_t magnetic_vector_3d = {(float)vlx.rsp.magnetic_transition[rsp.selected].x, (float)vlx.rsp.magnetic_transition[rsp.selected].y, (float)vlx.rsp.magnetic_transition[rsp.selected].z};
                    vec3_t electronic_vector_3d = {(float)vlx.rsp.electronic_transition_length[rsp.selected].x, (float)vlx.rsp.electronic_transition_length[rsp.selected].y, (float)vlx.rsp.electronic_transition_length[rsp.selected].z};

                    float el_ma_dot = vec3_dot(magnetic_vector_3d, electronic_vector_3d);
                    float el_len = vec3_length(electronic_vector_3d);
                    float ma_len = vec3_length(magnetic_vector_3d);
                    float angle = acos(el_ma_dot / (el_len * ma_len)) * (180.0 / 3.141592653589793238463);
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
            }

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

            bool reset_view = false;
            if (is_hovered) {
                if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                    reset_view = true;
                }
            }

            if (reset_view) {
                camera_compute_optimal_view(&nto.target.pos, &nto.target.ori, &nto.target.dist, min_aabb, max_aabb, nto.distance_scale);
            }

            if (is_active || is_hovered) {
                static const TrackballControllerParam param = {
                    .min_distance = 1.0,
                    .max_distance = 1000.0,
                };

                vec2_t delta = { io.MouseDelta.x, io.MouseDelta.y };
                vec2_t curr = {mouse_pos_in_canvas.x, mouse_pos_in_canvas.y};
                vec2_t prev = curr - delta;
                float  wheel_delta = io.MouseWheel;

                TrackballControllerInput input = {
                    .rotate_button = is_active && ImGui::IsMouseDown(ImGuiMouseButton_Left),
                    .pan_button    = is_active && ImGui::IsMouseDown(ImGuiMouseButton_Right),
                    .dolly_button  = is_active && ImGui::IsMouseDown(ImGuiMouseButton_Middle),
                    .dolly_delta   = is_hovered ? wheel_delta : 0.0f,
                    .mouse_coord_prev = prev,
                    .mouse_coord_curr = curr,
                    .screen_size = {canvas_sz.x, canvas_sz.y},
                    .fov_y = nto.camera.fov_y,
                };
                camera_controller_trackball(&nto.target.pos, &nto.target.ori, &nto.target.dist, input, param);
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

            const float aspect_ratio = win_sz.x / win_sz.y;

            mat4_t view_mat = camera_world_to_view_matrix(nto.camera);
            mat4_t proj_mat = camera_perspective_projection_matrix(nto.camera, aspect_ratio);
            mat4_t inv_proj_mat = camera_inverse_perspective_projection_matrix(nto.camera, aspect_ratio);



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
                    .depth      = nto.gbuf.tex.depth,
                    .color      = nto.gbuf.tex.color,
                    .normal     = nto.gbuf.tex.normal,
                    .velocity   = nto.gbuf.tex.velocity,
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

            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
            glDrawBuffer(GL_BACK);
            glDisable(GL_SCISSOR_TEST);

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
        }
        ImGui::End();
    }
};
static VeloxChem instance = {};
