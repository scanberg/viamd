#define IMGUI_DEFINE_MATH_OPERATORS

#include <event.h>
#include <viamd.h>
#include <task_system.h>
#include <color_utils.h>

#include <md_gto.h>
#include <md_vlx.h>
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

    bool show_volume = false;
    bool show_window = false;
    md_vlx_data_t vlx {};

    struct Volume {
        int dim[3] = {128, 128, 128};
        float* data = 0;
        size_t num_bytes = 0;
        vec3_t min_box = {};
        vec3_t max_box = {};
        vec3_t step_size = {};
        mat4_t model_mat = {};
    };

    Volume vol = {};

    struct {
        bool show_window = false;
        Volume vol = {};
        uint32_t vol_tex = {};

        struct {
            bool enabled = true;
            static const size_t count = 2;
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

    uint32_t vol_texture = 0;
    uint32_t tf_texture = 0;

    bool tf_dirty = true;

    vec3_t clip_min = {0,0,0};
    vec3_t clip_max = {1,1,1};

    bool   bounding_box_enabled = false;
    vec4_t bounding_box_color = {0,0,0,1};

    ImPlotColormap colormap = ImPlotColormap_Plasma;
    float colormap_alpha_scale = 1.0f;

    int mo_idx = 0;
    int homo_idx = 0;
    int lumo_idx = 0;

    task_system::ID compute_volume_task = 0;

    md_gto_data_t pgto = {0};

    float density_scale = 1.0f;

    struct {
        bool enabled = true;
        size_t count = 2;
        float  values[8] = {0.05f, -0.05};
        vec4_t colors[8] = {{215.f/255.f,25.f/255.f,28.f/255.f,0.75f}, {44.f/255.f,123.f/255.f,182.f/255.f,0.75f}};
    } iso;

    struct {
        bool enabled = false;
    } dvr;

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
                draw_window();
                draw_nto_window(state);
                break;
            }
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

                        for (int i = 0; i < (int)vlx.scf.alpha.occupations.count; ++i) {
                            if (vlx.scf.alpha.occupations.data[i] > 0) {
                                homo_idx = i;
                            } else {
                                lumo_idx = i;
                                break;
                            }
                        }

                        size_t num_pgtos = md_vlx_pgto_count(&vlx);
                        md_gto_data_init(&pgto, num_pgtos, arena);

                        vec3_t min_box = vec3_set1(FLT_MAX);
                        vec3_t max_box = vec3_set1(-FLT_MAX);

                        for (size_t i = 0; i < vlx.geom.num_atoms; ++i) {
                            vec3_t coord = vec3_set((float)vlx.geom.coord_x[i], (float)vlx.geom.coord_y[i], (float)vlx.geom.coord_z[i]);
                            min_box = vec3_min(min_box, coord);
                            max_box = vec3_max(max_box, coord);
                        }

                        const float pad = 3.0f;
                        min_box = vec3_sub_f(min_box, pad);
                        max_box = vec3_add_f(max_box, pad);

                        vol.min_box = min_box;
                        vol.max_box = max_box;
                        vol.model_mat = volume::compute_model_to_world_matrix(vol.min_box, vol.max_box);

                        update_volume();
                        show_volume = true;

                        // NTO
                        md_gl_representation_init(&nto.gl_rep, &state.mold.gl_mol);

                        size_t temp_pos = md_temp_get_pos();
                        uint32_t* colors = (uint32_t*)md_temp_push(state.mold.mol.atom.count * sizeof(uint32_t));
                        color_atoms_cpk(colors, state.mold.mol.atom.count, state.mold.mol);
                        md_gl_representation_set_color(&nto.gl_rep, 0, (uint32_t)state.mold.mol.atom.count, colors, 0);
                        md_temp_set_pos_back(temp_pos);
                        nto.show_window = true;

                        nto.vol.min_box = vol.min_box;
                        nto.vol.max_box = vol.max_box;
                        nto.vol.model_mat = vol.model_mat;

                        camera_compute_optimal_view(&nto.target.pos, &nto.target.ori, &nto.target.dist, min_box, max_box);
                    } else {
                        MD_LOG_INFO("Failed to load VeloxChem data");
                        md_arena_allocator_reset(arena);
                        vol = {};
                        vlx = {};
                        pgto = {};
                        nto = {};
                        mo_idx = 0;
                        homo_idx = 0;
                        lumo_idx = 0;
                        show_volume = false;
                    }
                }
                break;
            }
            case viamd::EventType_ViamdTopologyFree:
                md_gl_representation_free(&nto.gl_rep);

                md_arena_allocator_reset(arena);
                vol = {};
                vlx = {};
                pgto = {};
                nto = {};
                mo_idx = 0;
                homo_idx = 0;
                lumo_idx = 0;
                show_volume = false;
                break;
            default:
                break;
            }
        }
    }

    void update_volume() {
        if (task_system::task_is_running(compute_volume_task)) {
            task_system::task_interrupt(compute_volume_task);
        }
        else {
#define BLK_DIM 8
            vec3_t ext_box = vec3_sub(vol.max_box, vol.min_box);

            // Target resolution per spatial unit (We want this number of samples per Ångström in each dimension)
            // Round up to some multiple of 8 -> 8x8x8 is the chunksize we process in parallel

            const float res_scale = vol_res_scl[vol_res_idx];

            int dim[3] = {
                CLAMP(ALIGN_TO((int)(ext_box.x * res_scale), 8), 8, 512),
                CLAMP(ALIGN_TO((int)(ext_box.y * res_scale), 8), 8, 512),
                CLAMP(ALIGN_TO((int)(ext_box.z * res_scale), 8), 8, 512),
            };

            size_t num_bytes = sizeof(float) * dim[0] * dim[1] * dim[2];
            if (num_bytes != vol.num_bytes) {
                vol.data = (float*)md_realloc(arena, vol.data, vol.num_bytes, num_bytes);
                ASSERT(vol.data);
                vol.num_bytes = num_bytes;
            }

            MEMSET(vol.data, 0, vol.num_bytes);
            
            MEMCPY(vol.dim, dim, sizeof(vol.dim));
            vol.step_size = vec3_div(ext_box, vec3_set((float)vol.dim[0], (float)vol.dim[1], (float)vol.dim[2]));

            MD_LOG_DEBUG("Created Orbital volume of dimensions [%i][%i][%i]", vol.dim[0], vol.dim[1], vol.dim[2]);

            if (!md_vlx_extract_alpha_mo_pgtos(&pgto, &vlx, mo_idx)) {
                MD_LOG_ERROR("Failed to extract alpha orbital for orbital index: %i", mo_idx);
                return;
            }

            uint32_t num_blocks = (vol.dim[0] / BLK_DIM) * (vol.dim[1] / BLK_DIM) * (vol.dim[2] / BLK_DIM);
            // We evaluate the in parallel over smaller NxNxN blocks

            compute_volume_task = task_system::pool_enqueue(STR_LIT("Compute Volume"), 0, num_blocks, [](uint32_t range_beg, uint32_t range_end, void* user_data) {
                VeloxChem* data = (VeloxChem*)user_data;

                int num_blk[3] = {
                    data->vol.dim[0] / BLK_DIM,
                    data->vol.dim[1] / BLK_DIM,
                    data->vol.dim[2] / BLK_DIM,
                };

                md_grid_t grid = {
                    .data = data->vol.data,
                    .dim  = {data->vol.dim[0], data->vol.dim[1], data->vol.dim[2]},
                    // Shift origin by half a voxel to evaluate at the voxel center
                    .origin = {data->vol.min_box.x + 0.5f * data->vol.step_size.x, data->vol.min_box.y + 0.5f * data->vol.step_size.y, data->vol.min_box.z + 0.5f * data->vol.step_size.z},
                    .stepsize = {data->vol.step_size.x, data->vol.step_size.y, data->vol.step_size.z},
                };

                // Conversion from Ångström to Bohr
                const float factor = 1.0 / 0.529177210903;

                grid.origin[0] *= factor;
                grid.origin[1] *= factor;
                grid.origin[2] *= factor;

                grid.stepsize[0] *= factor;
                grid.stepsize[1] *= factor;
                grid.stepsize[2] *= factor;

                for (uint32_t i = range_beg; i < range_end; ++i) {
                    // Determine block index
                    int blk_x =  i % num_blk[0];
                    int blk_y = (i / num_blk[0]) % num_blk[1];
                    int blk_z =  i / (num_blk[0] * num_blk[1]);

                    const int off_idx[3] = {blk_x * BLK_DIM, blk_y * BLK_DIM, blk_z * BLK_DIM};
                    const int len_idx[3] = {BLK_DIM, BLK_DIM, BLK_DIM};

                    md_gto_eval_mode_t eval_mode = vol_mode == 0 ? MD_GTO_EVAL_MODE_PSI : MD_GTO_EVAL_MODE_PSI_SQUARED;
                    md_gto_grid_evaluate_sub(&grid, off_idx, len_idx, &data->pgto, eval_mode);
                }
            }, this);

#undef BLK_DIM

            // Launch task for main (render) thread to update the volume texture
            task_system::main_enqueue(STR_LIT("Update Volume"), [](void* user_data) {
                VeloxChem* data = (VeloxChem*)user_data;
                gl::init_texture_3D(&data->vol_texture, data->vol.dim[0], data->vol.dim[1], data->vol.dim[2], GL_R32F);
                gl::set_texture_3D_data(data->vol_texture, data->vol.data, GL_R32F);
            }, this, compute_volume_task);
        }
    }

    void draw_volume(ApplicationState& state) {
        if (!show_volume) return;

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
                .model = vol.model_mat,
                .view = state.view.param.matrix.current.view,
                .proj = state.view.param.matrix.current.proj,
            },
            .clip_volume = {
                .min = clip_min,
                .max = clip_max,
            },
            .global_scaling = {
                .density = density_scale,
            },
            .iso_surface = {
                .count  = iso.count,
                .values = iso.values,
                .colors = iso.colors,
            },
            .isosurface_enabled = iso.enabled,
            .direct_volume_rendering_enabled = dvr.enabled,
            .voxel_spacing = vol.step_size,
        };

        volume::render_volume(desc);

        if (bounding_box_enabled) {
            if (vol.model_mat != mat4_t{0}) {
                immediate::set_model_view_matrix(mat4_mul(state.view.param.matrix.current.view, vol.model_mat));
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
                                update_volume();
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
                                update_volume();
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
                    iso.values[0] =  (float)iso_val;
                    iso.values[1] = -(float)iso_val;
                    iso.count = 2;
                    iso.enabled = true;
                    dvr.enabled = false;

                    density_scale = 1.0f;

                    ImGui::ColorEdit4("Color Positive", iso.colors[0].elem);
                    ImGui::ColorEdit4("Color Negative", iso.colors[1].elem);
                } else {
                    iso.enabled = false;
                    dvr.enabled = true;

                    const ImVec2 button_size = {160, 0};
                    if (ImPlot::ColormapButton(ImPlot::GetColormapName(colormap), button_size, colormap)) {
                        ImGui::OpenPopup("Colormap Selector");
                    }
                    if (ImGui::BeginPopup("Colormap Selector")) {
                        for (int map = 4; map < ImPlot::GetColormapCount(); ++map) {
                            if (ImPlot::ColormapButton(ImPlot::GetColormapName(map), button_size, map)) {
                                colormap = map;
                                tf_dirty = true;
                                ImGui::CloseCurrentPopup();
                            }
                        }
                        ImGui::EndPopup();
                    }
                    if (ImGui::SliderFloat("Alpha Scale", &colormap_alpha_scale, 0.001f, 10.f, "%.3f", ImGuiSliderFlags_Logarithmic)) {
                        tf_dirty = true;
                    }

                    ImGui::SliderFloat("Density Scale", &density_scale, 0.001f, 10.f, "%.3f", ImGuiSliderFlags_Logarithmic);

                    // Update colormap texture
                    if (tf_dirty) {
                        tf_dirty = false;
                        uint32_t pixel_data[128];
                        for (size_t i = 0; i < ARRAY_SIZE(pixel_data); ++i) {
                            float t = (float)i / (float)(ARRAY_SIZE(pixel_data) - 1);
                            ImVec4 col = ImPlot::SampleColormap(t, colormap);

                            // This is a small alpha ramp in the start of the TF to avoid rendering low density values.
                            col.w = MIN(160 * t*t, 0.341176f);
                            col.w = CLAMP(col.w * colormap_alpha_scale, 0.0f, 1.0f);
                            pixel_data[i] = ImGui::ColorConvertFloat4ToU32(col);
                        }

                        gl::init_texture_2D(&tf_texture, (int)ARRAY_SIZE(pixel_data), 1, GL_RGBA8);
                        gl::set_texture_2D_data(tf_texture, pixel_data, GL_RGBA8);
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
                    ImGui::SetNextWindowScroll(ImVec2(-1, (TEXT_BASE_HEIGHT * (homo_idx + 1.5f)) - table_height * 0.5f));
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
                        bool is_selected = mo_idx == n;
                        ImGui::TableNextColumn();
                        char buf[32];
                        snprintf(buf, sizeof(buf), "%i", n + 1);
                        ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                        if (ImGui::Selectable(buf, is_selected, selectable_flags)) {
                            mo_idx = n;
                            update_volume();
                        }
                        ImGui::TableNextColumn();
                        ImGui::Text("%.1f", vlx.scf.alpha.occupations.data[n]);
                        ImGui::TableNextColumn();
                        ImGui::Text("%.4f", vlx.scf.alpha.energies.data[n]);
                        ImGui::TableNextColumn();
                        const char* lbl = (n == homo_idx) ? "HOMO" : (n == lumo_idx) ? "LUMO" : "";
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
                camera_compute_optimal_view(&nto.target.pos, &nto.target.ori, &nto.target.dist, nto.vol.min_box, nto.vol.max_box);

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
                        .model = nto.vol.model_mat,
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
