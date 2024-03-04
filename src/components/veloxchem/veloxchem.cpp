#define IMGUI_DEFINE_MATH_OPERATORS

#include <event.h>
#include <viamd.h>
#include <core/md_arena_allocator.h>
#include <task_system.h>

#include <core/md_log.h>
#include <md_vlx.h>

#include <gfx/volumerender_utils.h>
#include <gfx/gl_utils.h>
#include <gfx/immediate_draw_utils.h>

#include <imgui_widgets.h>
#include <implot_widgets.h>

struct VeloxChem : viamd::EventHandler {
    VeloxChem() { viamd::event_system_register_handler(*this); }

    bool show_window = false;
    md_vlx_data_t vlx {};
    int vol_dim = 128;
    mat4_t model_mat = {};

    float* vol_data = 0;
    vec3_t min_box = {};
    vec3_t max_box = {};
    vec3_t step_size = {};

    uint32_t vol_texture = 0;
    uint32_t tf_texture = 0;

    vec3_t clip_min = {0,0,0};
    vec3_t clip_max = {1,1,1};

    bool   bounding_box_enabled = true;
    vec4_t bounding_box_color = {0,0,0,1};

    int mo_idx = 0;

    task_system::ID compute_volume = 0;

    // This holds the entire mo coefficient matrix in a linear format
    // It is also transposed compared to the input such that reading a linear segment within it, you get the mo coefficients, not the ao coefficients
    md_array(double) mo_coeffs = 0;

    //double r_values[100] = {0};
    //double orb_values[5][100] = {0};
    //double r_density [3][100] = {0};

    float density_scale = 1.0f;

    struct {
        bool enabled = true;
        size_t count = 2;
        float  values[8] = {0.05f, -0.05};
        vec4_t colors[8] = {{1,0,0,1}, {0,0,1,1}};
    } iso_surface;

    md_allocator_i* arena = 0;

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
            case viamd::EventType_ViamdFrameTick:
                draw_window();
                break;
            case viamd::EventType_ViamdDrawMenu:
                ImGui::Checkbox("VeloxChem", &show_window);
                break;
            case viamd::EventType_ViamdPostRender: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                draw_volume(state);
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
                        show_window = true;

                        size_t num_rows = vlx.scf.alpha.orbitals.dim[0];
                        size_t num_cols = vlx.scf.alpha.orbitals.dim[1];

                        md_array_resize(mo_coeffs, num_rows * num_cols, arena);

                        for (size_t i = 0; i < num_rows; ++i) {
                            for (size_t j = 0; j < num_cols; ++j) {
                                mo_coeffs[i * num_rows + j] = vlx.scf.alpha.orbitals.data[j * num_cols + i];
                            }
                        }

                        vol_data = (float*)md_arena_allocator_push(arena, sizeof(float) * vol_dim * vol_dim * vol_dim);

                        min_box = vec3_set1(FLT_MAX);
                        max_box = vec3_set1(-FLT_MAX);

                        for (size_t i = 0; i < vlx.geom.num_atoms; ++i) {
                            vec3_t coord = vec3_set((float)vlx.geom.coord_x[i], (float)vlx.geom.coord_y[i], (float)vlx.geom.coord_z[i]);
                            min_box = vec3_min(min_box, coord);
                            max_box = vec3_max(max_box, coord);
                        }

                        const float pad = 2.0f;
                        min_box = vec3_sub_f(min_box, pad);
                        max_box = vec3_add_f(max_box, pad);
                        step_size = vec3_div_f(vec3_sub(max_box, min_box), (float)vol_dim);

                        model_mat = volume::compute_model_to_world_matrix(min_box, max_box);

                        update_volume();

                        /*
                        vec3_t coords[100];
                        for (size_t i = 0; i < ARRAY_SIZE(coords); ++i) {
                            double z = lerp(0.0, 2.0, i / (double)(ARRAY_SIZE(coords) - 1));
                            r_values[i] = z;
                            coords[i] = vec3_set(0,0, z * 0.529177210903);
                        }
                        for (int mo_idx = 0; mo_idx < 5; ++mo_idx) {
                            for (size_t i = 0; i < num_mo_coeffs; ++i) {
                                mo_coeffs[i] = vlx.scf.alpha.orbitals.data[i * num_cols + mo_idx];
                            }
                            md_vlx_get_mo(orb_values[mo_idx], coords, ARRAY_SIZE(coords), &vlx.geom, &vlx.basis, mo_coeffs, num_rows);
                            const double s = sign(orb_values[mo_idx][10]);
                            for (size_t i = 0; i < ARRAY_SIZE(coords); ++i) {
                                orb_values[mo_idx][i] *= s;
                            }
                        }
                        for (size_t i = 0; i < ARRAY_SIZE(r_values); ++i) {
                            const double r2 = pow(r_values[i], 2);
                            r_density[0][i] = 4 * PI * 2 * r2 *  pow(orb_values[0][i], 2);
                            r_density[1][i] = 4 * PI * 2 * r2 *  pow(orb_values[1][i], 2);
                            r_density[2][i] = 4 * PI * 2 * r2 * (pow(orb_values[2][i], 2) + pow(orb_values[3][i],2) + pow(orb_values[4][i],2));
                        }
                        */
                    } else {
                        MD_LOG_INFO("Failed to load VeloxChem data");
                        md_vlx_data_free(&vlx);
                        vlx = {};
                    }
                }
                break;
            }
            default:
                break;
            }
        }
    }

    void update_volume() {
        task_system::task_wait_for(compute_volume);

        uint32_t num_blocks = (vol_dim / 8) * (vol_dim / 8) * (vol_dim / 8);
        // We evaluate the in parallel over smaller 8x8x8 blocks

        compute_volume = task_system::pool_enqueue(STR_LIT("Compute Volume"), 0, num_blocks, [](uint32_t range_beg, uint32_t range_end, void* user_data) {
            VeloxChem* data = (VeloxChem*)user_data;

            int block_dim = data->vol_dim / 8;

            md_vlx_grid_t grid = {
                .data = data->vol_data,
                .dim  = {data->vol_dim, data->vol_dim, data->vol_dim},
                .origin = {data->min_box.x + 0.5f * data->step_size.x, data->min_box.y + 0.5f * data->step_size.y, data->min_box.z + 0.5f * data->step_size.z},
                .stepsize = {data->step_size.x, data->step_size.y, data->step_size.z},
            };

            for (uint32_t i = range_beg; i < range_end; ++i) {
                // Determine block index
                int blk_z = i & (block_dim - 1);
                int blk_y = (i / block_dim) & (block_dim - 1);
                int blk_x = i / (block_dim * block_dim);

                const int beg_idx[3] = {blk_x * 8, blk_y * 8, blk_z * 8};
                const int end_idx[3] = {blk_x * 8 + 8, blk_y * 8 + 8, blk_z * 8 + 8};

                const size_t num_mo_coeffs = data->vlx.scf.alpha.orbitals.dim[1];
                const double* mo_coeff = data->mo_coeffs + data->vlx.scf.alpha.orbitals.dim[0] * data->mo_idx;
                md_vlx_grid_evaluate_sub(&grid, beg_idx, end_idx, &data->vlx.geom, &data->vlx.basis, mo_coeff, num_mo_coeffs);
            }
        }, this);

        task_system::main_enqueue(STR_LIT("Update Volume"), [](void* user_data) {
            VeloxChem* data = (VeloxChem*)user_data;
            gl::init_texture_3D(&data->vol_texture, data->vol_dim, data->vol_dim, data->vol_dim, GL_R32F);
            gl::set_texture_3D_data(data->vol_texture, data->vol_data, GL_R32F);
        }, this, compute_volume);
    }

    void draw_volume(ApplicationState& state) {
        volume::RenderDesc desc = {
            .render_target = {
                .depth = state.gbuffer.deferred.depth,
                .color = state.gbuffer.deferred.post_tonemap,
                .normal = state.gbuffer.deferred.normal,
                .width   = state.gbuffer.width,
                .height  = state.gbuffer.height,
            },
            .texture = {
                .volume = vol_texture,
                .transfer_function = tf_texture,
            },
            .matrix = {
                .model = model_mat,
                .view = state.view.param.matrix.current.view,
                .proj = state.view.param.matrix.current.proj_jittered,
            },
            .clip_volume = {
                .min = clip_min,
                .max = clip_max,
            },
            .global_scaling = {
                .density = density_scale,
            },
            .iso_surface = {
                .count  = iso_surface.count,
                .values = iso_surface.values,
                .colors = iso_surface.colors,
            },
            .isosurface_enabled = iso_surface.enabled,
            .direct_volume_rendering_enabled = false,
            .voxel_spacing = step_size,
        };

        volume::render_volume(desc);

        if (bounding_box_enabled) {
            if (model_mat != mat4_t{0}) {
                immediate::set_model_view_matrix(mat4_mul(state.view.param.matrix.current.view, model_mat));
                immediate::set_proj_matrix(state.view.param.matrix.current.proj);
                immediate::draw_box_wireframe(vec3_set1(0), vec3_set1(1), bounding_box_color);
                immediate::render();
            }
        }
    }

    void draw_window() {
        if (!show_window) return;
        ImGui::SetNextWindowSize({300,350}, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("VeloxChem", &show_window)) {
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
            ImGui::SliderFloat("Iso Value 0", &iso_surface.values[0], -1.0f, 1.0f);
            ImGui::ColorEdit4("Color 0", iso_surface.colors[0].elem);

            ImGui::SliderFloat("Iso Value 1", &iso_surface.values[1], -1.0f, 1.0f);
            ImGui::ColorEdit4("Color 1", iso_surface.colors[1].elem);

            char buf[32];
            snprintf(buf, sizeof(buf), "%i", mo_idx);
            if (ImGui::BeginCombo("MO index", buf)) {
                for (int n = 0; n < (int)vlx.scf.alpha.orbitals.dim[0]; n++) {
                    snprintf(buf, sizeof(buf), "%i", n);
                    if (ImGui::Selectable(buf, mo_idx == n)) {
                        mo_idx = n;
                        update_volume();
                    }
                }
                ImGui::EndCombo();
            }

            /*
            if (ImPlot::BeginPlot("MO")) {
                ImPlot::SetupAxes("Distance to nucleus (Bohr)", "Orbital value (a.u.)");
                ImPlot::PlotLine("1s",  r_values, orb_values[0], ARRAY_SIZE(r_values));
                ImPlot::PlotLine("2s",  r_values, orb_values[1], ARRAY_SIZE(r_values));
                ImPlot::PlotLine("2p1", r_values, orb_values[2], ARRAY_SIZE(r_values));
                ImPlot::PlotLine("2p2", r_values, orb_values[3], ARRAY_SIZE(r_values));
                ImPlot::PlotLine("2p3", r_values, orb_values[4], ARRAY_SIZE(r_values));
                ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Density")) {
                ImPlot::SetupAxes("Distance to nucleus (Bohr)", "Radial densities (a.u.)");
                ImPlot::PlotLine("1s", r_values, r_density[0], ARRAY_SIZE(r_values));
                ImPlot::PlotLine("2s", r_values, r_density[1], ARRAY_SIZE(r_values));
                ImPlot::PlotLine("2p", r_values, r_density[2], ARRAY_SIZE(r_values));
                ImPlot::EndPlot();
            }
            */
        }
        ImGui::End();
    }
} instance;
