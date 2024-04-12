#define IMGUI_DEFINE_MATH_OPERATORS

#include <event.h>
#include <viamd.h>
#include <task_system.h>
#include <color_utils.h>

#include <md_gto.h>
#include <md_vlx.h>
#include <md_util.h>
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
        md_gl_representation_t gl_rep = {};
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
        md_gl_representation_t gl_rep = {};
    } nto;

    struct Rsp {
        bool show_window = false;
        int hovered = -1;
        int selected = -1;
        int focused_plot = -1;
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
                draw_scf_window();
                draw_rsp_window();
                break;
            }
            case viamd::EventType_ViamdDrawMenu:
                ImGui::Checkbox("VeloxChem Orbital", &orb.show_window);
                ImGui::Checkbox("VeloxChem NTO", &nto.show_window);
                ImGui::Checkbox("VeloxChem RSP", &rsp.show_window);
                ImGui::Checkbox("VeloxChem SCF", &scf.show_window);
                break;
            case viamd::EventType_ViamdRenderTransparent: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                //ApplicationState& state = *(ApplicationState*)e.payload;
                //draw_orb_volume(state);
                break;
            }
            case viamd::EventType_ViamdTopologyInit: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState& state = *(ApplicationState*)e.payload;
                str_t ext;
                str_t top_file = str_from_cstr(state.files.molecule);
                if (extract_ext(&ext, top_file) && str_eq_ignore_case(ext, STR_LIT("out"))) {
                    MD_LOG_INFO("Attempting to load VeloxChem data from file '" STR_FMT "'", STR_ARG(top_file));
                    md_vlx_data_free(&vlx);
                    if (md_vlx_data_parse_file(&vlx, top_file, arena)) {
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

                        vec3_t com = md_util_com_compute_vec4(xyzw, vlx.geom.num_atoms, 0);
                        mat3_t C = mat3_covariance_matrix_vec4(xyzw, 0, vlx.geom.num_atoms, com);
                        mat3_eigen_t eigen = mat3_eigen(C);
                        PCA = mat3_extract_rotation(eigen.vectors);

                        uint32_t* colors = (uint32_t*)md_vm_arena_push(state.allocator.frame, state.mold.mol.atom.count * sizeof(uint32_t));
                        color_atoms_cpk(colors, state.mold.mol.atom.count, state.mold.mol);

                        // NTO
                        if (vlx.rsp.num_excited_states > 0 && vlx.rsp.nto) {
                            nto.show_window = true;
                            md_gl_representation_init(&nto.gl_rep, &state.mold.gl_mol);
                            md_gl_representation_set_color(&nto.gl_rep, 0, (uint32_t)state.mold.mol.atom.count, colors, 0);
                            camera_compute_optimal_view(&nto.target.pos, &nto.target.ori, &nto.target.dist, min_aabb, max_aabb, nto.distance_scale);
                        }

                        // RSP
                        rsp.show_window = true;
                        rsp.hovered  = -1;
                        rsp.selected = -1;

                        // ORB
                        orb.show_window = true;
                        md_gl_representation_init(&orb.gl_rep, &state.mold.gl_mol);
                        md_gl_representation_set_color(&orb.gl_rep, 0, (uint32_t)state.mold.mol.atom.count, colors, 0);
                        camera_compute_optimal_view(&orb.target.pos, &orb.target.ori, &orb.target.dist, min_aabb, max_aabb, orb.distance_scale);
                        orb.mo_idx = homo_idx;
                        orb.scroll_to_idx = homo_idx;

                    } else {
                        MD_LOG_INFO("Failed to load VeloxChem data");
                        md_arena_allocator_reset(arena);
                        vlx = {};
                        orb = {};
                        nto = {};
                        rsp = {};
                    }
                }
                break;
            }
            case viamd::EventType_ViamdTopologyFree:
                md_gl_representation_free(&nto.gl_rep);

                md_arena_allocator_reset(arena);
                vlx = {};
                orb = {};
                nto = {};
                rsp = {};
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

        task_system::ID async_task = task_system::pool_enqueue(STR_LIT("Evaluate Orbital"), 0, num_blocks, [](uint32_t range_beg, uint32_t range_end, void* user_data, uint32_t thread_num) {
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
        task_system::ID main_task = task_system::main_enqueue(STR_LIT("##Update Volume"), [](void* user_data) {
            Payload* data = (Payload*)user_data;
            
            // The init here is just to ensure that the volume has not changed its dimensions during the async evaluation
            gl::init_texture_3D(data->tex_ptr, data->vol_dim[0], data->vol_dim[1], data->vol_dim[2], GL_R16F);
            gl::set_texture_3D_data(*data->tex_ptr, data->vol_data, GL_R32F);

            md_free(md_get_heap_allocator(), data->pgtos, data->num_pgtos * sizeof(md_gto_t));
            md_free(md_get_heap_allocator(), data, data->bytes);
        }, payload);

        task_system::set_task_dependency(main_task, async_task);

        return async_task;
    }

    static inline double axis_conversion_multiplier(const double* y1_array, const double* y2_array, size_t array_size) {
        double y1_min = y1_array[0];
        double y1_max = y1_array[0];
        double y2_min = y2_array[0];
        double y2_max = y2_array[0];
        for (size_t i = 1; i < array_size; i++) {
            y1_min = MIN(y1_min, y1_array[i]);
            y1_max = MAX(y1_max, y1_array[i]);
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

    static inline void convert_values(double* values, size_t num_values, x_unit_t unit) {
        switch (unit) {
        case X_UNIT_EV: break; // Do nothing
        case X_UNIT_NM:
            for (size_t i = 0; i < num_values; ++i) {
                values[i] = 1239.84193 / values[i];
            }
            break;
        case X_UNIT_CM_INVERSE:
            for (size_t i = 0; i < num_values; ++i) {
                values[i] = values[i] * 8065.73;
            }
            break;
        case X_UNIT_HARTREE:
            for (size_t i = 0; i < num_values; ++i) {
                values[i] = values[i] * 0.0367502;
            }
            break;
        default:
            ASSERT(false); // Should not happen
            break;
        }
    }

    static inline void broaden_gaussian(double* out_y, const double* in_x, size_t num_samples, const double* in_x_peaks, const double* in_y_peaks, size_t num_peaks, double sigma) {
        for (size_t xi = 0; xi < num_samples; xi++) {
            double tot = 0.0;
            for (size_t pi = 0; pi < num_peaks; pi++) {
                tot += in_y_peaks[pi] * exp(-(pow(in_x_peaks[pi] - in_x[xi], 2) / (2 * pow(sigma, 2))));
            }
            out_y[xi] = tot;
        }
    }

    static inline void broaden_lorentzian(double* out_y, const double* in_x, size_t num_samples, const double* in_x_peaks, const double* in_y_peaks, size_t num_peaks, double sigma) {
        for (size_t xi = 0; xi < num_samples; xi++) {
            double tot = 0.0;
            for (size_t pi = 0; pi < num_peaks; pi++) {
                tot += in_y_peaks[pi] / (1 + pow((in_x[xi] - in_x_peaks[pi]) / sigma, 2));
            }
            out_y[xi] = tot;
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

    //converts x and y peaks into pixel points in context of the current plot. Use between BeginPlot() and EndPlot()
    static inline void peaks_to_pixels(ImVec2* pixel_peaks, const double* x_peaks, const double* y_peaks, size_t num_peaks) {
        for (size_t i = 0; i < num_peaks; i++) {
            pixel_peaks[i] = ImPlot::PlotToPixels(ImPlotPoint{ x_peaks[i], y_peaks[i] });
        }
    }
    // Returns peak index closest to mouse pixel position, assumes that x-values are sorted.
    static inline int get_hovered_peak(const ImVec2 mouse_pos, const ImVec2* pixel_peaks, size_t num_peaks, double proxy_distance = 10.0) {
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
            }
            else if (y < y_min) {
                distance_y = fabs(y - y_min);
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

    static inline void cgs_to_ecd(double* ecd_out_peaks, const double* x_peaks, const double* cgs_peaks, size_t num_peaks) {
        double inv = 1 / (22.94 * PI);
        for (size_t i = 0; i < num_peaks; i++) {
            ecd_out_peaks[i] = x_peaks[i] * cgs_peaks[i] * inv;
        }
    }

    void draw_scf_window() {
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
        double y1_to_y2_mult = axis_conversion_multiplier(vlx.scf.iter.gradient_norm, energy_offsets, vlx.scf.iter.count);

        ImGui::SetNextWindowSize({ 300, 350 }, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("SCF", &scf.show_window)) {
            if (ImPlot::BeginPlot("SCF")) {
                ImPlot::SetupAxisLimits(ImAxis_X1, 1.0, (int)vlx.scf.iter.count);
                ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);
                ImPlot::SetupAxes("Iteration", "Gradient Norm (a.u.)");
                // We draw 2 y axis as "Energy total" has values in a different range then the rest of the data
                ImPlot::SetupAxis(ImAxis_Y2, "Energy (Hartree)", ImPlotAxisFlags_AuxDefault);
#if 1
                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                //ImPlot::SetupAxisScale(ImAxis_Y2, ImPlotScale_Log10);
#endif
                //ImPlot::SetupAxisLimits(ImAxis_Y2, lims.Y.Min * y1_to_y2_mult, lims.Y.Max * y1_to_y2_mult, ImPlotCond_Always);


                ImPlot::PlotLine("Gradient", iter, vlx.scf.iter.gradient_norm, (int)vlx.scf.iter.count);
                ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
                ImPlot::PlotLine("Energy", iter, energy_offsets, (int)vlx.scf.iter.count - 1);
                lims = ImPlot::GetPlotLimits(ImAxis_X1, ImAxis_Y1);
                //ImPlot::PlotLine("Density Change", iter, vlx.scf.iter.density_change, (int)vlx.scf.iter.count);
                //ImPlot::PlotLine("Energy Change", iter, vlx.scf.iter.energy_change, (int)vlx.scf.iter.count);
                //ImPlot::PlotLine("Max Gradient", iter, vlx.scf.iter.max_gradient, (int)vlx.scf.iter.count);
                ImPlot::EndPlot();
            }
        }
        ImGui::End();
    }

    void draw_rsp_window() {
        if (!rsp.show_window) return;
        if (vlx.rsp.num_excited_states == 0) return;

        // Keep track of the temp position so we can reset to it after we are done
        size_t temp_pos = md_temp_get_pos();
        defer { md_temp_set_pos_back(temp_pos); };

        static float sigma = 0.1;
        static ImVec2 mouse_pos = { 0,0 };

        const char* broadening_str[] = { "Gaussian","Lorentzian" };
        static broadening_mode_t broadening_mode = BROADENING_GAUSSIAN;

        const char* x_unit_str[] = { "eV", "nm", "cm-1", "hartree"};
        static x_unit_t x_unit = X_UNIT_EV;

        ImGui::SetNextWindowSize({ 300, 350 }, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("RSP", &rsp.show_window)) {
            bool refit = false;
            static bool first_plot = true;
            
            ImGui::SliderFloat((const char*)u8"Broadening (σ)", &sigma, 0.01f, 1.0f);
            refit |= ImGui::Combo("Broadening mode", (int*)(&broadening_mode), broadening_str, IM_ARRAYSIZE(broadening_str));
            refit |= ImGui::Combo("X unit", (int*)(&x_unit), x_unit_str, IM_ARRAYSIZE(x_unit_str));
            
            const int   num_peaks = (int)vlx.rsp.num_excited_states;
            double*       x_peaks = (double*)md_temp_push(sizeof(double) * num_peaks);
            const double* y_osc_peaks = vlx.rsp.absorption_osc_str;
            const double* y_cgs_peaks = vlx.rsp.electronic_circular_dichroism_cgs;
            for (int i = 0; i < num_peaks; ++i) {
                x_peaks[i] = vlx.rsp.absorption_ev[i];
            }

            double* y_ecd_peaks = (double*)md_temp_push(sizeof(double) * num_peaks);

            cgs_to_ecd(y_ecd_peaks, x_peaks, y_cgs_peaks, num_peaks);

            const int num_samples = 1024;
            double* x_values    = (double*)md_temp_push(sizeof(double) * num_samples);
            double* y_osc_str = (double*)md_temp_push(sizeof(double) * num_samples);
            double* y_cgs_str = (double*)md_temp_push(sizeof(double) * num_samples);
            double* y_ecd_str   = (double*)md_temp_push(sizeof(double) * num_samples);

            ImVec2*   pixel_osc_peaks = (ImVec2*)md_temp_push(sizeof(ImVec2) * num_peaks);
            ImVec2*   pixel_cgs_peaks = (ImVec2*)md_temp_push(sizeof(ImVec2) * num_peaks);

            const double x_min = vlx.rsp.absorption_ev[0] - 1.0;
            const double x_max = vlx.rsp.absorption_ev[num_peaks - 1] + 1.0;
            for (int i = 0; i < num_samples; ++i) {
                double t = (double)i / (double)(num_samples - 1);
                double value = lerp(x_min, x_max, t);
                x_values[i] = value;
            }

            // @NOTE: Do broadening in eV
            switch (broadening_mode) {
            case BROADENING_GAUSSIAN:
                broaden_gaussian(y_osc_str, x_values, num_samples, x_peaks, y_osc_peaks, num_peaks, sigma);
                broaden_gaussian(y_cgs_str, x_values, num_samples, x_peaks, y_cgs_peaks, num_peaks, sigma);
                broaden_gaussian(y_ecd_str, x_values, num_samples, x_peaks, y_ecd_peaks, num_peaks, sigma);
                break;
            case BROADENING_LORENTZIAN:
                broaden_lorentzian(y_osc_str, x_values, num_samples, x_peaks, y_osc_peaks, num_peaks, sigma);
                broaden_lorentzian(y_cgs_str, x_values, num_samples, x_peaks, y_cgs_peaks, num_peaks, sigma);
                broaden_lorentzian(y_ecd_str, x_values, num_samples, x_peaks, y_ecd_peaks, num_peaks, sigma);
                break;
            default:
                ASSERT(false); // Should not happen
                break;
            }

            // Do conversions
            convert_values(x_peaks,  num_peaks,   x_unit);
            convert_values(x_values, num_samples, x_unit);

            // Calulate constraint limits for the plot
            double y_osc_max_con = 0;
            double y_osc_min_con = 0;
            double y_cgs_max_con = 0;
            double y_cgs_min_con = 0;
            double x_max_con = x_values[num_samples - 1];
            double x_min_con = x_values[0];
            for (int i = 0; i < num_samples; i++) {
                y_osc_max_con = MAX(y_osc_max_con, y_osc_str[i]);
                y_osc_min_con = MIN(y_osc_min_con, y_osc_str[i]);
                y_cgs_max_con = MAX(y_cgs_max_con, y_cgs_str[i]);
                y_cgs_min_con = MIN(y_cgs_min_con, y_cgs_str[i]);
                x_max_con = MAX(x_max_con, x_values[i]);
                x_min_con = MIN(x_min_con, x_values[i]);
            }
            double con_lim_fac = 0.1;
            double y_osc_graph_width = y_osc_max_con - y_osc_min_con;
            double y_cgs_graph_width = y_cgs_max_con - y_cgs_min_con;
            double x_graph_width = x_max_con - x_min_con;
            y_osc_max_con += con_lim_fac * y_osc_graph_width;
            y_osc_min_con -= con_lim_fac * y_osc_graph_width;
            y_cgs_max_con += con_lim_fac * y_cgs_graph_width;
            y_cgs_min_con -= con_lim_fac * y_cgs_graph_width;
            x_max_con += con_lim_fac * x_graph_width;
            x_min_con -= con_lim_fac * x_graph_width;


#if 0
            //Hovered display text
            if (rsp.hovered != -1 && rsp.focused_plot == 0) {
                ImGui::BulletText("Hovered: %s = %f, Y = %f", x_unit_str[x_unit], (float)x_peaks[rsp.hovered], (float)y_osc_peaks[rsp.hovered]);

            }
            else if (rsp.hovered != -1 && rsp.focused_plot == 1){
                ImGui::BulletText("Hovered: %s = %f, Y = %f", x_unit_str[x_unit], (float)x_peaks[rsp.hovered], (float)y_cgs_peaks[rsp.hovered]);
            }
            else {
                ImGui::BulletText("Hovered:");
            }

            //Selected display text
            if (rsp.selected != -1) {
                ImGui::BulletText("Selected: %s = %f, Y_OSC = %f, Y_CGS = %f", x_unit_str[x_unit], (float)x_peaks[rsp.selected], (float)y_osc_peaks[rsp.selected], (float)y_osc_peaks[rsp.selected]);
            }
            else {
                ImGui::BulletText("Selected:");
            }
            ImGui::BulletText("Mouse: X = %f, Y = %f", mouse_pos.x, mouse_pos.y);
            ImGui::BulletText("Peak 0: X = %f, Y = %f", (float)pixel_osc_peaks[0].x, (float)pixel_osc_peaks[0].y);
            ImGui::BulletText("Closest Index = %i", rsp.hovered);
#endif
            rsp.focused_plot = -1;
            if (ImPlot::BeginSubplots("##AxisLinking", 2, 1, ImVec2(-1, -1), ImPlotSubplotFlags_LinkCols)) {
                if (refit || first_plot) { ImPlot::SetNextAxesToFit(); }
                // Absorption
                if (ImPlot::BeginPlot("Absorption")) {
                    ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                    ImPlot::SetupAxes(x_unit_str[x_unit], "Oscillator Strength");
                    ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, x_min_con, x_max_con);
                    ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, y_osc_min_con, y_osc_max_con);

                    peaks_to_pixels(pixel_osc_peaks, x_peaks, y_osc_peaks, num_peaks);
                    mouse_pos = ImPlot::PlotToPixels(ImPlot::GetPlotMousePos(IMPLOT_AUTO));
                    if (ImPlot::IsPlotHovered()) {
                        rsp.hovered = get_hovered_peak(mouse_pos, pixel_osc_peaks, num_peaks);
                        rsp.focused_plot = 0;
                    }

                    // @HACK: Compute pixel width of 2 'plot' units
                    const double bar_width = ImPlot::PixelsToPlot(ImVec2(2, 0)).x - ImPlot::PixelsToPlot(ImVec2(0, 0)).x;

                    ImPlot::PlotBars("Exited States", x_peaks, y_osc_peaks, num_peaks, bar_width);
                    ImPlot::PlotLine("Oscillator Strength", x_values, y_osc_str, num_samples);
                    //Check hovered state
                    if (rsp.hovered != -1) {
                        draw_bar(0, x_peaks[rsp.hovered], y_osc_peaks[rsp.hovered], bar_width, ImVec4{ 0,1,0,1 });
                    }

                    // Update selected peak on click
                    if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Left) && ImPlot::IsPlotHovered()) {
                        rsp.selected = rsp.hovered;
                    }
                    //Check selected state
                    if (rsp.selected != -1) {
                        draw_bar(1, x_peaks[rsp.selected], y_osc_peaks[rsp.selected], bar_width, ImVec4{ 1,0,0,1 });
                    }
                    
                    ImPlot::EndPlot();
                }

                if (refit || first_plot) { ImPlot::SetNextAxesToFit(); }
                // Rotatory ECD
                if (ImPlot::BeginPlot("ECD")) {
                    ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                    ImPlot::SetupAxes(x_unit_str[x_unit], "Rotatory Strength");
                    ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"ECD (L mol−1 cm−1)", ImPlotAxisFlags_AuxDefault);
                    ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, x_min_con, x_max_con);
                    ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, y_cgs_min_con, y_cgs_max_con);
                    peaks_to_pixels(pixel_cgs_peaks, x_peaks, y_cgs_peaks, num_peaks);
                    mouse_pos = ImPlot::PlotToPixels(ImPlot::GetPlotMousePos(IMPLOT_AUTO));

                    if (ImPlot::IsPlotHovered()) { 
                        rsp.hovered = get_hovered_peak(mouse_pos, pixel_cgs_peaks, num_peaks);
                        rsp.focused_plot = 1;
                    }
                    // @HACK: Compute pixel width of 2 'plot' units

                    const double bar_width = ImPlot::PixelsToPlot(ImVec2(2, 0)).x - ImPlot::PixelsToPlot(ImVec2(0, 0)).x;

                    ImPlot::PlotBars("Exited States", x_peaks, y_cgs_peaks, num_peaks, bar_width);
                    ImPlot::PlotLine("CGS", x_values, y_cgs_str, num_samples);
                    ImPlot::SetAxis(ImAxis_Y2);
                    ImPlot::PlotLine("ECD", x_values, y_ecd_str, num_samples);

                    if (rsp.hovered != -1 && ImPlot::IsPlotHovered()) {
                        draw_bar(2, x_peaks[rsp.hovered], y_cgs_peaks[rsp.hovered], bar_width, ImVec4{ 0,1,0,1 });
                    }

                    // Update selected peak on click
                    if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Left) && ImPlot::IsPlotHovered()) {
                        rsp.selected = rsp.hovered;
                    }
                    if (rsp.selected != -1) {
                        draw_bar(3, x_peaks[rsp.selected], y_cgs_peaks[rsp.selected], bar_width, ImVec4{ 1,0,0,1 });
                    }
                    ImPlot::EndPlot();
                }
                ImPlot::EndSubplots();
            }
            first_plot = false;
            /*
            constexpr str_t ABS_FILE_EXTENSION = STR_LIT("abs");
            char path_buf[2048] = "";


            if (ImGui::Button("Print absorption")) {
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, ABS_FILE_EXTENSION)) {
                    // This is where we save the absorbtion into a file
                    save_absorption({ path_buf, strnlen(path_buf, sizeof(path_buf)) }, &x_values, x_unit_str[x_unit], &y_osc_str, &y_cgs_str, 10);
                }
            }
            */
        }
        ImGui::End();
    }

    void draw_orb_window(const ApplicationState& state) {
        if (!orb.show_window) return;
        ImGui::SetNextWindowSize({600,300}, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("VeloxChem Orbital Viewer", &orb.show_window)) {
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
            draw_op.rep = &orb.gl_rep;

            md_gl_draw_args_t draw_args = {
                .shaders = &state.mold.gl_shaders,
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
                    .color = {24.f, 24.f, 24.f, 1.0f},
                },
                .bloom = {
                    .enabled = false,
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
                .temporal_reprojection = {
                    .enabled = false,
                },
                .sharpen = {
                    .enabled = true,
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
                        .voxel_spacing = orb.vol[i].step_size,
                        };
                        volume::render_volume(vol_desc);
                    }
                POP_GPU_SECTION();
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

            const double iso_min = 1.0e-4;
            const double iso_max = 5.0;
            double iso_val = nto.iso.values[0];
            ImGui::SliderScalar("##Iso Value", ImGuiDataType_Double, &iso_val, &iso_min, &iso_max, "%.6f", ImGuiSliderFlags_Logarithmic);
            ImGui::SetItemTooltip("Iso Value");

            nto.iso.values[0] =  (float)iso_val;
            nto.iso.values[1] = -(float)iso_val;
            nto.iso.count = 2;
            nto.iso.enabled = true;

            ImGui::ColorEdit4("##Color Positive", nto.iso.colors[0].elem);
            ImGui::SetItemTooltip("Color Positive");
            ImGui::ColorEdit4("##Color Negative", nto.iso.colors[1].elem);
            ImGui::SetItemTooltip("Color Negative");

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
            draw_op.rep = &nto.gl_rep;

            md_gl_draw_args_t draw_args = {
                .shaders = &state.mold.gl_shaders,
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
                    .color = {24.f, 24.f, 24.f, 1.0f},
                },
                .bloom = {
                    .enabled = false,
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
                .temporal_reprojection = {
                    .enabled = false,
                },
                .sharpen = {
                    .enabled = true,
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
                PUSH_GPU_SECTION("ORB GRID RAYCAST")
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
