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

static const char* vol_res_lbl[3] = {
    "Low",
    "Mid",
    "High",
};

static const float vol_res_scl[3] = {
    4.0f,
    8.0f,
    16.0f,
};

static int vol_res_idx = 1;
static int vol_mode = 0;

struct VeloxChem : viamd::EventHandler {
    VeloxChem() { viamd::event_system_register_handler(*this); }

    md_vlx_data_t vlx {};

    struct Volume {
        int dim[3] = {128, 128, 128};
        float* data = 0;
        size_t num_bytes = 0;
        vec3_t step_size = {};
        vec3_t extent = {};
        mat4_t model_to_world = {};
        mat4_t tex_to_world   = {};
    };

    struct Scf {
        bool show_window = false;
    } scf;

    struct Orb {
        bool show_window = false;
        bool show_volume = false;
        Volume vol = {};

        uint32_t vol_texture = 0;
        uint32_t tf_texture = 0;

        bool tf_dirty = true;

        vec3_t clip_min = {0,0,0};
        vec3_t clip_max = {1,1,1};

        struct {
            bool enabled = true;
            vec4_t color = {0,0,0,1};
        } bounding_box;

        ImPlotColormap colormap = ImPlotColormap_Plasma;
        float colormap_alpha_scale = 1.0f;

        int mo_idx = 0;
        int homo_idx = 0;
        int lumo_idx = 0;

        task_system::ID compute_volume_task = 0;

        md_array(md_gto_t) pgtos = 0;

        float density_scale = 1.0f;

        struct {
            bool enabled = true;
            size_t count = 2;
            float  values[2] = {0.05f, -0.05};
            vec4_t colors[2] = {{215.f/255.f,25.f/255.f,28.f/255.f,0.75f}, {44.f/255.f,123.f/255.f,182.f/255.f,0.75f}};
        } iso;

        struct {
            bool enabled = false;
        } dvr;
    } orb;

    struct Nto {
        bool show_window = false;
        Volume vol = {};
        uint32_t vol_tex = {};

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

        bool show_coordinate_system_widget = true;
        md_gl_representation_t gl_rep = {};
    } nto;

    struct Rsp {
        bool show_window = false;
    } rsp;

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
                ApplicationState& state = *(ApplicationState*)e.payload;
                arena = md_arena_allocator_create(state.allocator.persistent, MEGABYTES(1));
                break;
            }
            case viamd::EventType_ViamdShutdown:
                md_arena_allocator_destroy(arena);
                break;
            case viamd::EventType_ViamdFrameTick: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                draw_orb_window();
                draw_nto_window(state);
                draw_scf_window();
                draw_rsp_window();
                break;
            }
            case viamd::EventType_ViamdDrawMenu:
                ImGui::Checkbox("VeloxChem Orbital", &orb.show_window);
                ImGui::Checkbox("VeloxChem RSP", &rsp.show_window);
                ImGui::Checkbox("VeloxChem SCF", &scf.show_window);
                break;
            case viamd::EventType_ViamdPostRender: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                draw_orb_volume(state);
                break;
            }
            case viamd::EventType_ViamdTopologyInit: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                str_t ext;
                str_t top_file = str_from_cstr(state.files.molecule);
                if (extract_ext(&ext, top_file) && str_eq_ignore_case(ext, STR_LIT("out"))) {
                    MD_LOG_INFO("Attempting to load VeloxChem data from file '" STR_FMT "'", STR_ARG(top_file));
                    md_vlx_data_free(&vlx);
                    if (md_vlx_data_parse_file(&vlx, top_file, arena)) {
                        MD_LOG_INFO("Successfully loaded VeloxChem data");

                        // Scf
                        scf.show_window = true;

                        // Orb
                        orb.show_window = true;

                        for (int i = 0; i < (int)vlx.scf.alpha.occupations.count; ++i) {
                            if (vlx.scf.alpha.occupations.data[i] > 0) {
                                orb.homo_idx = i;
                            } else {
                                orb.lumo_idx = i;
                                break;
                            }
                        }

                        size_t num_pgtos = md_vlx_pgto_count(&vlx);
                        md_array_resize(orb.pgtos, num_pgtos, arena);
                        MEMSET(orb.pgtos, 0, md_array_bytes(orb.pgtos));

                        // Compute Object Oriented Bounding Box
                        // @NOTE: This is a crude approximation, by using the PCA
                        vec4_t* xyzw = (vec4_t*)md_temp_push(sizeof(vec4_t) * vlx.geom.num_atoms);
                        for (size_t i = 0; i < vlx.geom.num_atoms; ++i) {
                            xyzw[i] = {(float)vlx.geom.coord_x[i], (float)vlx.geom.coord_y[i], (float)vlx.geom.coord_z[i], 1.0f};
                        }

                        vec3_t com = md_util_com_compute_vec4(xyzw, vlx.geom.num_atoms, 0);
                        mat3_t C = mat3_covariance_matrix_vec4(xyzw, 0, vlx.geom.num_atoms, com);
                        mat3_eigen_t eigen = mat3_eigen(C);

                        mat3_t A = mat3_extract_rotation(eigen.vectors);

                        // Compute min and maximum extent along the axes
                        vec3_t min_ext = { FLT_MAX, FLT_MAX, FLT_MAX};
                        vec3_t max_ext = {-FLT_MAX,-FLT_MAX,-FLT_MAX};

                        mat4_t R = mat4_from_mat3(mat3_transpose(A));
                        mat4_t Ri = mat4_from_mat3(A);

                        for (size_t i = 0; i < vlx.geom.num_atoms; ++i) {
                            vec4_t p = mat4_mul_vec4(Ri, xyzw[i]);
                            min_ext = vec3_min(min_ext, vec3_from_vec4(p));
                            max_ext = vec3_max(max_ext, vec3_from_vec4(p));
                        }

                        const float pad = 2.5f;
                        min_ext = vec3_sub_f(min_ext, pad);
                        max_ext = vec3_add_f(max_ext, pad);

                        vec3_t origin = mat4_mul_vec3(R, min_ext, 0.0f);
                        vec3_t extent = {max_ext.x - min_ext.x, max_ext.y - min_ext.y, max_ext.z - min_ext.z};

                        mat4_t T = mat4_translate(origin.x, origin.y, origin.z);
                        mat4_t S = mat4_scale(extent.x, extent.y, extent.z);

                        orb.vol.extent = extent;
                        orb.vol.model_to_world = T * R;
                        orb.vol.tex_to_world = T * R * S;

                        update_orb_volume();
                        orb.show_volume = true;

                        // NTO
                        md_gl_representation_init(&nto.gl_rep, &state.mold.gl_mol);

                        size_t temp_pos = md_temp_get_pos();
                        uint32_t* colors = (uint32_t*)md_temp_push(state.mold.mol.atom.count * sizeof(uint32_t));
                        color_atoms_cpk(colors, state.mold.mol.atom.count, state.mold.mol);
                        md_gl_representation_set_color(&nto.gl_rep, 0, (uint32_t)state.mold.mol.atom.count, colors, 0);
                        md_temp_set_pos_back(temp_pos);
                        nto.vol.model_to_world = orb.vol.model_to_world;
                        nto.vol.tex_to_world   = orb.vol.tex_to_world;

                        vec3_t min_box = vec3_from_vec4(orb.vol.model_to_world * vec4_set(0,0,0,1));
                        vec3_t max_box = vec3_from_vec4(orb.vol.model_to_world * vec4_set(1,1,1,1));
                        camera_compute_optimal_view(&nto.target.pos, &nto.target.ori, &nto.target.dist, min_box, max_box);

                        // RSP
                        rsp.show_window = true;
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
            default:
                break;
            }
        }
    }

    void update_orb_volume() {
        if (task_system::task_is_running(orb.compute_volume_task)) {
            task_system::task_interrupt(orb.compute_volume_task);
        }
        else {
#define BLK_DIM 8
            // Target resolution per spatial unit (We want this number of samples per Ångström in each dimension)
            // Round up to some multiple of 8 -> 8x8x8 is the chunksize we process in parallel

            const float res_scale = vol_res_scl[vol_res_idx];

            int dim[3] = {
                CLAMP(ALIGN_TO((int)(orb.vol.extent.x * res_scale), 8), 8, 512),
                CLAMP(ALIGN_TO((int)(orb.vol.extent.y * res_scale), 8), 8, 512),
                CLAMP(ALIGN_TO((int)(orb.vol.extent.z * res_scale), 8), 8, 512),
            };

            size_t num_bytes = sizeof(float) * dim[0] * dim[1] * dim[2];
            if (num_bytes != orb.vol.num_bytes) {
                orb.vol.data = (float*)md_realloc(arena, orb.vol.data, orb.vol.num_bytes, num_bytes);
                ASSERT(orb.vol.data);
                orb.vol.num_bytes = num_bytes;
            }

            MEMSET(orb.vol.data, 0, orb.vol.num_bytes);
            
            MEMCPY(orb.vol.dim, dim, sizeof(orb.vol.dim));
            orb.vol.step_size = vec3_div(orb.vol.extent, vec3_set((float)orb.vol.dim[0], (float)orb.vol.dim[1], (float)orb.vol.dim[2]));

            MD_LOG_DEBUG("Created Orbital volume of dimensions [%i][%i][%i]", orb.vol.dim[0], orb.vol.dim[1], orb.vol.dim[2]);

            if (!md_vlx_extract_alpha_mo_pgtos(orb.pgtos, &vlx, orb.mo_idx)) {
                MD_LOG_ERROR("Failed to extract alpha orbital for orbital index: %i", orb.mo_idx);
                return;
            }
            md_gto_cutoff_compute(orb.pgtos, md_array_size(orb.pgtos), 1.0e-6);

            // Transform pgtos into the 'model space' of the volume
            // There is a scaling factor here as well as the PGTOs are given and must be
            // Processed in a.u. (Bohr) and not Ångström, in which the volume is defined
            mat4_t R = orb.vol.model_to_world;
            R.col[3] = vec4_set(0,0,0,1);
            R = mat4_transpose(R);

            // Ångström to Bohr
            vec4_t t = vec4_mul_f(orb.vol.model_to_world.col[3], 1.0 / 0.529177210903);
            mat4_t T = mat4_translate(-t.x, -t.y, -t.z);

            mat4_t M = R * T;
            for (size_t i = 0; i < md_array_size(orb.pgtos); ++i) {
                vec4_t c = {orb.pgtos[i].x, orb.pgtos[i].y, orb.pgtos[i].z, 1.0f};
                c = M * c;
                orb.pgtos[i].x = c.x;
                orb.pgtos[i].y = c.y;
                orb.pgtos[i].z = c.z;
            }

            // We evaluate the in parallel over smaller NxNxN blocks
            uint32_t num_blocks = (orb.vol.dim[0] / BLK_DIM) * (orb.vol.dim[1] / BLK_DIM) * (orb.vol.dim[2] / BLK_DIM);

            orb.compute_volume_task = task_system::pool_enqueue(STR_LIT("Compute Volume"), 0, num_blocks, [](uint32_t range_beg, uint32_t range_end, void* user_data) {
                Orb* orb = (Orb*)user_data;

                // Number of NxNxN blocks in each dimension
                int num_blk[3] = {
                    orb->vol.dim[0] / BLK_DIM,
                    orb->vol.dim[1] / BLK_DIM,
                    orb->vol.dim[2] / BLK_DIM,
                };

                // Conversion from Ångström to Bohr
                const float factor = 1.0 / 0.529177210903;

                // The PGTOs have been pre-transformed into the 'model space' of the volume
                vec3_t step = vec3_mul_f(orb->vol.step_size, factor);

                md_grid_t grid = {
                    .data = orb->vol.data,
                    .dim  = {orb->vol.dim[0], orb->vol.dim[1], orb->vol.dim[2]},
                    // Shift origin by half a voxel to evaluate at the voxel center
                    .origin = {0.5f * step.x, 0.5f * step.y, 0.5f * step.z},
                    .stepsize = {step.x, step.y, step.z},
                };

                for (uint32_t i = range_beg; i < range_end; ++i) {
                    // Determine block index
                    int blk_x =  i % num_blk[0];
                    int blk_y = (i / num_blk[0]) % num_blk[1];
                    int blk_z =  i / (num_blk[0] * num_blk[1]);

                    const int off_idx[3] = {blk_x * BLK_DIM, blk_y * BLK_DIM, blk_z * BLK_DIM};
                    const int len_idx[3] = {BLK_DIM, BLK_DIM, BLK_DIM};

                    md_gto_eval_mode_t eval_mode = vol_mode == 0 ? MD_GTO_EVAL_MODE_PSI : MD_GTO_EVAL_MODE_PSI_SQUARED;
                    md_gto_grid_evaluate_sub(&grid, off_idx, len_idx, orb->pgtos, md_array_size(orb->pgtos), eval_mode);
                }
            }, &this->orb);

#undef BLK_DIM

            // Launch task for main (render) thread to update the volume texture
            task_system::main_enqueue(STR_LIT("Update Volume"), [](void* user_data) {
                Orb* orb = (Orb*)user_data;
                gl::init_texture_3D(&orb->vol_texture, orb->vol.dim[0], orb->vol.dim[1], orb->vol.dim[2], GL_R32F);
                gl::set_texture_3D_data(orb->vol_texture, orb->vol.data, GL_R32F);
            }, &this->orb, orb.compute_volume_task);
        }
    }

    void draw_orb_volume(ApplicationState& state) {
        if (!orb.show_volume) return;

        // Volume model matrix expects the 

        volume::RenderDesc desc = {
            .render_target = {
                .depth = state.gbuffer.deferred.depth,
                .color = state.gbuffer.deferred.post_tonemap,
                .normal = state.gbuffer.deferred.normal,
                .width   = state.gbuffer.width,
                .height  = state.gbuffer.height,
            },
            .texture = {
                .volume = orb.vol_texture,
                .transfer_function = orb.tf_texture,
            },
            .matrix = {
                .model = orb.vol.tex_to_world,
                .view = state.view.param.matrix.current.view,
                .proj = state.view.param.matrix.current.proj,
            },
            .clip_volume = {
                .min = orb.clip_min,
                .max = orb.clip_max,
            },
            .global_scaling = {
                .density = orb.density_scale,
            },
            .iso_surface = {
                .count  = orb.iso.count,
                .values = orb.iso.values,
                .colors = orb.iso.colors,
            },
            .isosurface_enabled = orb.iso.enabled,
            .direct_volume_rendering_enabled = orb.dvr.enabled,
            .voxel_spacing = orb.vol.step_size,
        };

        volume::render_volume(desc);

        if (orb.bounding_box.enabled) {
            if (orb.vol.model_to_world != mat4_t{0}) {
                immediate::set_model_view_matrix(mat4_mul(state.view.param.matrix.current.view, orb.vol.tex_to_world));
                immediate::set_proj_matrix(state.view.param.matrix.current.proj);
                immediate::draw_box_wireframe(vec3_set1(0), vec3_set1(1), orb.bounding_box.color);
                immediate::render();
            }
        }
    }

    void draw_orb_window() {
        if (!orb.show_window) return;
        ImGui::SetNextWindowSize({300,350}, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("VeloxChem", &orb.show_window)) {
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
            if (ImGui::TreeNode("Molecular Orbitals")) {
                static const char* vol_mode_str[2] = {
                    (const char*)u8"Orbital (Ψ)",
                    (const char*)u8"Density (Ψ²)"
                };
                if (ImGui::BeginCombo("Mode", vol_mode_str[vol_mode])) {
                    for (int i = 0; i < ARRAY_SIZE(vol_mode_str); ++i) {
                        if (ImGui::Selectable(vol_mode_str[i], vol_mode == i)) {
                            if (vol_mode != i) {
                                vol_mode = i;
                                update_orb_volume();
                            }
                        }
                    }
                    ImGui::EndCombo();
                }

                if (ImGui::BeginCombo("Volume Resolution", vol_res_lbl[vol_res_idx])) {
                    for (int i = 0; i < ARRAY_SIZE(vol_res_lbl); ++i) {
                        if (ImGui::Selectable(vol_res_lbl[i], vol_res_idx == i)) {
                            if (vol_res_idx != i) {
                                vol_res_idx = i;
                                update_orb_volume();
                            }
                        }
                    }
                    ImGui::EndCombo();
                }

                if (vol_mode == 0) {
                    static double iso_val = 0.05f;
                    const  double iso_min = 1.0e-4;
                    const  double iso_max = 5.0;
                    ImGui::SliderScalar("Iso Value", ImGuiDataType_Double, &iso_val, &iso_min, &iso_max, "%.6f", ImGuiSliderFlags_Logarithmic);
                    orb.iso.values[0] =  (float)iso_val;
                    orb.iso.values[1] = -(float)iso_val;
                    orb.iso.count = 2;
                    orb.iso.enabled = true;
                    orb.dvr.enabled = false;

                    orb.density_scale = 1.0f;

                    ImGui::ColorEdit4("Color Positive", orb.iso.colors[0].elem);
                    ImGui::ColorEdit4("Color Negative", orb.iso.colors[1].elem);
                } else {
                    orb.iso.enabled = false;
                    orb.dvr.enabled = true;

                    const ImVec2 button_size = {160, 0};
                    if (ImPlot::ColormapButton(ImPlot::GetColormapName(orb.colormap), button_size, orb.colormap)) {
                        ImGui::OpenPopup("Colormap Selector");
                    }
                    if (ImGui::BeginPopup("Colormap Selector")) {
                        for (int map = 4; map < ImPlot::GetColormapCount(); ++map) {
                            if (ImPlot::ColormapButton(ImPlot::GetColormapName(map), button_size, map)) {
                                orb.colormap = map;
                                orb.tf_dirty = true;
                                ImGui::CloseCurrentPopup();
                            }
                        }
                        ImGui::EndPopup();
                    }
                    if (ImGui::SliderFloat("Alpha Scale", &orb.colormap_alpha_scale, 0.001f, 10.f, "%.3f", ImGuiSliderFlags_Logarithmic)) {
                        orb.tf_dirty = true;
                    }

                    ImGui::SliderFloat("Density Scale", &orb.density_scale, 0.001f, 10.f, "%.3f", ImGuiSliderFlags_Logarithmic);

                    // Update colormap texture
                    if (orb.tf_dirty) {
                        orb.tf_dirty = false;
                        uint32_t pixel_data[128];
                        for (size_t i = 0; i < ARRAY_SIZE(pixel_data); ++i) {
                            float t = (float)i / (float)(ARRAY_SIZE(pixel_data) - 1);
                            ImVec4 col = ImPlot::SampleColormap(t, orb.colormap);

                            // This is a small alpha ramp in the start of the TF to avoid rendering low density values.
                            col.w = MIN(160 * t*t, 0.341176f);
                            col.w = CLAMP(col.w * orb.colormap_alpha_scale, 0.0f, 1.0f);
                            pixel_data[i] = ImGui::ColorConvertFloat4ToU32(col);
                        }

                        gl::init_texture_2D(&orb.tf_texture, (int)ARRAY_SIZE(pixel_data), 1, GL_RGBA8);
                        gl::set_texture_2D_data(orb.tf_texture, pixel_data, GL_RGBA8);
                    }
                }

                const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();

                enum {
                    Col_Idx,
                    Col_Occ,
                    Col_Ene,
                    Col_Lbl,
                };

                static float table_height = 0.0f;
                if (ImGui::Button("Goto HOMO/LUMO")) {
                    // 1.5f offset = 1 for skipping the column header row and 0.5 for placing the view between homo/lumo
                    ImGui::SetNextWindowScroll(ImVec2(-1, (TEXT_BASE_HEIGHT * (orb.homo_idx + 1.5f)) - table_height * 0.5f));
                }

                const ImGuiTableFlags flags =
                    ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable | ImGuiTableFlags_RowBg |
                    ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersV | ImGuiTableFlags_NoBordersInBody | ImGuiTableFlags_ScrollY;
                if (ImGui::BeginTable("Molecular Orbitals", 4, flags))//, ImVec2(0.0f, TEXT_BASE_HEIGHT * 15), 0.0f))
                {
                    auto draw_row = [this](int n) {
                        // Display a data item
                        ImGui::PushID(n + 1);
                        ImGui::TableNextRow();
                        bool is_selected = (orb.mo_idx == n);
                        ImGui::TableNextColumn();
                        char buf[32];
                        snprintf(buf, sizeof(buf), "%i", n + 1);
                        ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                        if (ImGui::Selectable(buf, is_selected, selectable_flags)) {
                            orb.mo_idx = n;
                            update_orb_volume();
                        }
                        ImGui::TableNextColumn();
                        ImGui::Text("%.1f", vlx.scf.alpha.occupations.data[n]);
                        ImGui::TableNextColumn();
                        ImGui::Text("%.4f", vlx.scf.alpha.energies.data[n]);
                        ImGui::TableNextColumn();
                        const char* lbl = (n == orb.homo_idx) ? "HOMO" : (n == orb.lumo_idx) ? "LUMO" : "";
                        ImGui::TextUnformatted(lbl);
                        ImGui::PopID();
                        };

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
                    ImGui::TableSetupColumn("Label",        ImGuiTableColumnFlags_NoSort               | ImGuiTableColumnFlags_WidthStretch, 0.0f, Col_Lbl);
                    ImGui::TableSetupScrollFreeze(0, 1); // Make row always visible
                    ImGui::TableHeadersRow();

                    for (int n = 0; n < num_orbitals(); n++) {
                        draw_row(n);
                    }

                    ImGui::EndTable();
                    table_height = ImGui::GetItemRectSize().y;
                }
                ImGui::TreePop();
            }
        }
        ImGui::End();
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

        // The actual plot
        ImGui::SetNextWindowSize({ 300, 350 }, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("SCF", &scf.show_window)) {
            // We draw 2 plots as "Energy total" has values in a different range then the rest of the data
            if (ImPlot::BeginSubplots("##AxisLinking", 2, 1, ImVec2(-1, -1), ImPlotSubplotFlags_LinkCols)) {
                if (ImPlot::BeginPlot("SCF")) {

                    ImPlot::SetupAxisLimits(ImAxis_X1, 1.0, (int)vlx.scf.iter.count);
                    ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);
                    ImPlot::SetupAxes("Iterations", "eV");

                    ImPlot::PlotLine("Density Change", iter, vlx.scf.iter.density_change, (int)vlx.scf.iter.count);
                    ImPlot::PlotLine("Energy Change", iter, vlx.scf.iter.energy_change, (int)vlx.scf.iter.count);
                    ImPlot::PlotLine("Gradient Norm", iter, vlx.scf.iter.gradient_norm, (int)vlx.scf.iter.count);
                    ImPlot::PlotLine("Max Gradient", iter, vlx.scf.iter.max_gradient, (int)vlx.scf.iter.count);
                }
                ImPlot::EndPlot();

                if (ImPlot::BeginPlot("SCF")) {

                    ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);

                    ImPlot::PlotLine("Energy Total", iter, vlx.scf.iter.energy_total, (int)vlx.scf.iter.count);
                }
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
        }
        ImGui::End();
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

    void draw_rsp_window() {
        if (!rsp.show_window) return;
        if (vlx.rsp.num_excited_states == 0) return;

        // Keep track of the temp position so we can reset to it after we are done
        size_t temp_pos = md_temp_get_pos();
        defer { md_temp_set_pos_back(temp_pos); };

        static float sigma = 0.1;

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

            const int num_samples = 1024;
            double* x_values    = (double*)md_temp_push(sizeof(double) * num_samples);
            double* y_osc_str = (double*)md_temp_push(sizeof(double) * num_samples);
            double* y_cgs_str   = (double*)md_temp_push(sizeof(double) * num_samples);

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
                break;
            case BROADENING_LORENTZIAN:
                broaden_lorentzian(y_osc_str, x_values, num_samples, x_peaks, y_osc_peaks, num_peaks, sigma);
                broaden_lorentzian(y_cgs_str, x_values, num_samples, x_peaks, y_cgs_peaks, num_peaks, sigma);
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

            if (ImPlot::BeginSubplots("##AxisLinking", 2, 1, ImVec2(-1, -1), ImPlotSubplotFlags_LinkCols)) {
                if (refit || first_plot) { ImPlot::SetNextAxesToFit(); }
                // Absorption
                if (ImPlot::BeginPlot("Absorption")) {
                    // ImPlot::SetupAxisLimits(ImAxis_X1, 1.0, vlx.scf.iter.count);
                    ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                    ImPlot::SetupAxes(x_unit_str[x_unit], "Oscillator Strength");
                    ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, x_min_con, x_max_con);
                    ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, y_osc_min_con, y_osc_max_con);

                    // @HACK: Compute pixel width of 2 'plot' units
                    const double bar_width = ImPlot::PixelsToPlot(ImVec2(2, 0)).x - ImPlot::PixelsToPlot(ImVec2(0, 0)).x;

                    ImPlot::PlotBars("Exited States", x_peaks, y_osc_peaks, num_peaks, bar_width);
                    ImPlot::PlotLine("Oscillator Strength", x_values, y_osc_str, num_samples);

                }
                ImPlot::EndPlot();

                if (refit || first_plot) { ImPlot::SetNextAxesToFit(); }
                // Rotary ECD
                if (ImPlot::BeginPlot("ECD")) {
                    // ImPlot::SetupAxisLimits(ImAxis_X1, 1.0, vlx.scf.iter.count);
                    ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                    ImPlot::SetupAxes(x_unit_str[x_unit], "Rotary strength");
                    ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, x_min_con, x_max_con);
                    ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, y_cgs_min_con, y_cgs_max_con);

                    // @HACK: Compute pixel width of 2 'plot' units
                    const double bar_width = ImPlot::PixelsToPlot(ImVec2(2, 0)).x - ImPlot::PixelsToPlot(ImVec2(0, 0)).x;

                    ImPlot::PlotBars("Exited States", x_peaks, y_cgs_peaks, num_peaks, bar_width);
                    ImPlot::PlotLine("ECD", x_values, y_cgs_str, num_samples);

                }
                ImPlot::EndPlot();
            }
            ImPlot::EndSubplots();
            first_plot = false;
        }
        ImGui::End();
    }

    void draw_nto_window(const ApplicationState& state) {
        if (!nto.show_window) return;

        ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("NTO viewer", &nto.show_window, ImGuiWindowFlags_MenuBar)) {
            const ImVec2 button_size = {160, 0};

            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("Render")) {
                    ImGui::Checkbox("Enable Iso Surfaces", &nto.iso.enabled);
                    if (nto.iso.enabled) {
                        for (int i = 0; i < nto.iso.count; ++i) {
                            ImGui::PushID(i);
                            ImGui::SliderFloat("##Isovalue", &nto.iso.values[i], 0.0f, 10.f, "%.3f", ImGuiSliderFlags_Logarithmic);
                            ImGui::SameLine();
                            ImGui::ColorEdit4Minimal("##Color", nto.iso.colors[i].elem);
                            ImGui::PopID();
                        }
                    }
                    ImGui::EndMenu();
                }

                if (ImGui::BeginMenu("Show")) {
                    ImGui::Checkbox("Coordinate System Widget", &nto.show_coordinate_system_widget);
                    ImGui::EndMenu();
                }

                ImGui::EndMenuBar();
            }

            //update_density_volume(data);

            // Animate camera towards targets
            const float dt = 1.0f / 60.0f;
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

            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->AddImage((ImTextureID)(intptr_t)nto.gbuf.deferred.post_tonemap, canvas_p0, canvas_p1, { 0,1 }, { 1,0 });
            draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));

            const bool is_hovered = ImGui::IsItemHovered();
            const bool is_active = ImGui::IsItemActive();
            const ImVec2 origin(canvas_p0.x, canvas_p0.y);  // Lock scrolled origin
            const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

            auto& gbuf = nto.gbuf;
            if (gbuf.width != (uint32_t)canvas_sz.x || gbuf.height != (uint32_t)canvas_sz.y) {
                init_gbuffer(&gbuf, (int)canvas_sz.x, (int)canvas_sz.y);
            }

            bool reset_hard = false;
            bool reset_view = false;
            if (is_hovered) {
                if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                    reset_view = true;
                }
            }

            if (reset_view) {
                vec3_t min_box = mat4_mul_vec3(nto.vol.tex_to_world, vec3_set1(0.f), 1.f);
                vec3_t max_box = mat4_mul_vec3(nto.vol.tex_to_world, vec3_set1(1.f), 1.f);
                camera_compute_optimal_view(&nto.target.pos, &nto.target.ori, &nto.target.dist, min_box, max_box);

                if (reset_hard) {
                    nto.camera.position         = nto.target.pos;
                    nto.camera.orientation      = nto.target.ori;
                    nto.camera.focus_distance   = nto.target.dist;
                }
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
                ImVec2 win_size = ImGui::GetWindowSize();
                float  ext = MIN(win_size.x, win_size.y) * 0.2f;
                float  pad = 0.1f * ext;

                CoordSystemWidgetParam params = {
                    .pos = ImVec2(pad, win_size.y - ext - pad),
                    .size = {ext, ext},
                    .view_matrix = camera_world_to_view_matrix(nto.camera),
                    .camera_ori  = nto.target.ori,
                    .camera_pos  = nto.target.pos,
                    .camera_dist = nto.target.dist,
                };

                ImGui::DrawCoordinateSystemWidget(params);
            }

            const float aspect_ratio = canvas_sz.x / canvas_sz.y;
            mat4_t view_mat = camera_world_to_view_matrix(nto.camera);
            mat4_t proj_mat = camera_perspective_projection_matrix(nto.camera, aspect_ratio);
            mat4_t inv_proj_mat = camera_inverse_perspective_projection_matrix(nto.camera, aspect_ratio);

            clear_gbuffer(&gbuf);

            glEnable(GL_DEPTH_TEST);
            glDepthMask(GL_TRUE);

            const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
                GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_POST_TONEMAP };

            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);

            glEnable(GL_DEPTH_TEST);
            glDepthFunc(GL_LESS);

            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.deferred.fbo);
            glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
            glViewport(0, 0, gbuf.width, gbuf.height);

            md_gl_draw_op_t draw_op = {};
            draw_op.type = MD_GL_REP_BALL_AND_STICK;
            draw_op.args.ball_and_stick.ball_scale = 1.0f;
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
                    .projection_matrix = (const float*)proj_mat.elem,
                },
            };

            md_gl_draw(&draw_args);

            glDrawBuffer(GL_COLOR_ATTACHMENT_POST_TONEMAP);
            glClearColor(1, 1, 1, 0);
            glClear(GL_COLOR_BUFFER_BIT);

            if (nto.iso.enabled) {
                volume::RenderDesc vol_desc = {
                    .render_target = {
                        .depth  = nto.gbuf.deferred.depth,
                        .color  = nto.gbuf.deferred.post_tonemap,
                        .width  = nto.gbuf.width,
                        .height = nto.gbuf.height,
                    },
                    .texture = {
                        .volume = nto.vol_tex,
                    },
                    .matrix = {
                        .model = nto.vol.tex_to_world,
                        .view  = view_mat,
                        .proj  = proj_mat,
                    },
                    .iso_surface = {
                        .count  = (size_t)nto.iso.count,
                        .values = nto.iso.values,
                        .colors = nto.iso.colors,
                    },
                    .isosurface_enabled = nto.iso.enabled,
                    .voxel_spacing = nto.vol.step_size,
                };
                volume::render_volume(vol_desc);
            }

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
                .temporal_reprojection = {
                    .enabled = false,
                },
                .input_textures = {
                    .depth          = nto.gbuf.deferred.depth,
                    .color          = nto.gbuf.deferred.color,
                    .normal         = nto.gbuf.deferred.normal,
                    .velocity       = nto.gbuf.deferred.velocity,
                    .post_tonemap   = nto.gbuf.deferred.post_tonemap,
                }
            };

            ViewParam view_param = {
                .matrix = {
                    .current = {
                    .view = view_mat,
                    .proj = proj_mat,
                    .norm = view_mat,
                },
                .inverse = {
                    .proj = inv_proj_mat,
                }
                },
                .clip_planes = {
                    .near = nto.camera.near_plane,
                    .far  = nto.camera.far_plane,
                },
                .fov_y = nto.camera.fov_y,
                .resolution = {canvas_sz.x, canvas_sz.y}
            };

            postprocessing::shade_and_postprocess(postprocess_desc, view_param);
            POP_GPU_SECTION()

            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
            glDrawBuffer(GL_BACK);
        }

        ImGui::End();
    }
};
static VeloxChem instance = {};
