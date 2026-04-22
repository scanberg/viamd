#include <core/md_arena_allocator.h>

#include <viamd_event.h>
#include <viamd.h>

#include <gfx/gl_utils.h>
#include <gfx/volumerender_utils.h>
#include <gfx/immediate_draw_utils.h>
#include <color_utils.h>

#include <imgui_internal.h>
#include <imgui_widgets.h>

#include <implot_internal.h>
#include <implot_widgets.h>

constexpr uint64_t PickingSource_DensityVolume = HASH_STR_LIT64("Density Volume");

enum LegendColorMapMode_ {
    LegendColorMapMode_Opaque,
    LegendColorMapMode_Transparent,
    LegendColorMapMode_Split,
};

struct DensityVolume : viamd::EventHandler {
    bool show_window = false;
    bool enabled = false;

    md_allocator_i* arena = nullptr;

    struct {
        bool enabled = true;
        struct {
            uint32_t id = 0;
            float alpha_scale = 1.f;
            int colormap = DEFAULT_COLORMAP;
            bool dirty = true;
            float min_val = 0.0f;
            float max_val = 1.0f;
        } tf;
    } dvr;

    struct {
        bool enabled = false;
        float values[8] = {};
        vec4_t colors[8] = {};
        size_t count = 0;
    } iso;

    struct {
        uint32_t id = 0;
        bool dirty = false;
        int  dim[3] = {0};
        float max_value = 1.f;
    } volume_texture;


    struct {
        vec3_t min = {0, 0, 0};
        vec3_t max = {1, 1, 1};
    } clip_volume;

    struct {
        bool enabled = true;
        bool checkerboard = true;
        int  colormap_mode = LegendColorMapMode_Split;
    } legend;

    vec3_t voxel_spacing = {1.0f, 1.0f, 1.0f};
    float resolution_scale = 2.0f;

    vec4_t clip_volume_color = {1,0,0,1};
    vec4_t bounding_box_color = {0,0,0,1};

    bool show_bounding_box = true;
    bool show_reference_structures = true;
    bool show_reference_ensemble = false;
    bool show_density_volume = false;
    bool show_coordinate_system_widget = true;

    bool dirty_rep = false;
    bool dirty_vol = false;

    struct {
        RepresentationType type = RepresentationType::BallAndStick;
        ColorMapping colormap = ColorMapping::Type;
        float param[4] = {1,1,1,1};
        vec4_t color = {1,1,1,1};
    } rep;

    md_array(md_gl_rep_t) gl_reps = nullptr;
    md_array(mat4_t) rep_model_mats = nullptr;
    mat4_t model_mat = {0};

    GBuffer gbuf = {0};
    PickingSurface picking_surface = {};
    Camera camera = {};
	ViewTransform target = {};

    DensityVolume() {
        viamd::event_system_register_handler(*this);
    }

    void update(ApplicationState* state) {
        md_allocator_i* temp_arena = state->allocator.frame;
        md_vm_arena_temp_t temp_scope = md_vm_arena_temp_begin(temp_arena);
        defer { md_vm_arena_temp_end(temp_scope); };

        if (dvr.tf.dirty) {
            dvr.tf.dirty = false;
            // Update colormap texture
            volume::compute_transfer_function_texture_simple(&dvr.tf.id, dvr.tf.colormap, dvr.tf.alpha_scale);
        }

        int64_t selected_property = -1;
        for (size_t i = 0; i < md_array_size(state->display_properties); ++i) {
            const DisplayProperty& dp = state->display_properties[i];
            if (dp.type == DisplayProperty::Type_Volume && dp.show_in_volume) {
                selected_property = i;
                break;
            }
        }

        const md_script_property_data_t* prop_data = 0;
        const md_script_vis_payload_o* vis_payload = 0;
        uint64_t data_fingerprint = 0;

        static int64_t s_selected_property = 0;
        if (s_selected_property != selected_property) {
            s_selected_property = selected_property;
            dirty_vol = true;
            dirty_rep = true;
        }

        if (selected_property != -1) {
            prop_data = state->display_properties[selected_property].prop_data;
            vis_payload = state->display_properties[selected_property].vis_payload;
            data_fingerprint = state->display_properties[selected_property].prop_data->fingerprint;
        }
        show_density_volume = selected_property != -1;

        static uint64_t s_script_fingerprint = 0;
        if (s_script_fingerprint != md_script_ir_fingerprint(state->script.eval_ir)) {
            s_script_fingerprint = md_script_ir_fingerprint(state->script.eval_ir);
            dirty_vol = true;
            dirty_rep = true;
        }

        static uint64_t s_data_fingerprint = 0;
        if (s_data_fingerprint != data_fingerprint) {
            s_data_fingerprint = data_fingerprint;
            dirty_vol = true;
        }

        static double s_frame = 0;
        if (s_frame != state->animation.frame) {
            s_frame = state->animation.frame;
            dirty_rep = true;
        }

        if (dirty_rep) {
            if (prop_data && vis_payload) {
                dirty_rep = false;
                size_t num_reps = 0;
                bool result = false;
                md_script_vis_t vis = {};

                if (md_script_ir_valid(state->script.eval_ir)) {
                    md_script_vis_init(&vis, state->allocator.frame);
                    md_script_vis_ctx_t ctx = {
                        .ir = state->script.eval_ir,
                        .mol = &state->mold.sys,
                        .traj = state->mold.sys.trajectory,
                    };
                    result = md_script_vis_eval_payload(&vis, vis_payload, 0, &ctx, MD_SCRIPT_VISUALIZE_SDF);
                }

                if (result) {
                    if (vis.sdf.extent) {
                        const float s = vis.sdf.extent;
                        vec3_t min_aabb = { -s, -s, -s };
                        vec3_t max_aabb = { s, s, s };
                        model_mat = volume::compute_model_to_world_matrix(min_aabb, max_aabb);
                        voxel_spacing = vec3_t{2*s / prop_data->dim[1], 2*s / prop_data->dim[2], 2*s / prop_data->dim[3]};
                    }
                    num_reps = md_array_size(vis.sdf.structures);
                }

                // We need to limit this for performance reasons
                num_reps = MIN(num_reps, 100);

                const size_t old_size = md_array_size(gl_reps);
                if (gl_reps) {
                    // Only free superflous entries
                    for (size_t i = num_reps; i < old_size; ++i) {
                        md_gl_rep_destroy(gl_reps[i]);
                    }
                }
                md_array_resize(gl_reps, num_reps, arena);
                md_array_resize(rep_model_mats, num_reps, arena);

                for (size_t i = old_size; i < num_reps; ++i) {
                    // Only init new entries
                    gl_reps[i] = md_gl_rep_create(state->mold.gl_mol);
                }

                const md_system_t& sys = state->mold.sys;
			    const size_t num_atoms = md_system_atom_count(&sys);
                const size_t num_bytes = sizeof(uint32_t) * num_atoms;
                uint32_t* colors = (uint32_t*)md_vm_arena_push(temp_arena, num_bytes);

                switch (rep.colormap) {
                case ColorMapping::Uniform:
                    color_atoms_uniform(colors, num_atoms, convert_color(rep.color));
                    break;
                case ColorMapping::Type:
                    color_atoms_type(colors, num_atoms, sys);
                    break;
                case ColorMapping::Serial:
                    color_atoms_idx(colors, num_atoms, sys);
				    break;
                case ColorMapping::CompName:
                    color_atoms_comp_name(colors, num_atoms, sys);
                    break;
                case ColorMapping::CompSeqId:
                    color_atoms_comp_seq_id(colors, num_atoms, sys);
                    break;
                case ColorMapping::InstId:
                    color_atoms_inst_id(colors, num_atoms, sys);
                    break;
                case ColorMapping::InstIndex:
                    color_atoms_inst_idx(colors, num_atoms, sys);
                    break;
                case ColorMapping::SecondaryStructure:
                    color_atoms_secondary_structure(colors, num_atoms, sys);
                    break;
                default:
                    ASSERT(false);
                    break;
                }

                for (size_t i = 0; i < num_reps; ++i) {
                    filter_colors(colors, num_atoms, &vis.sdf.structures[i]);
                    md_gl_rep_set_color(gl_reps[i], 0, (uint32_t)num_atoms, colors, 0);
                    rep_model_mats[i] = vis.sdf.matrices[i];
                }
            }
        }

        if (dirty_vol) {
            if (prop_data) {
                dirty_vol = false;
                if (!volume_texture.id) {
                    int dim[3] = { prop_data->dim[1], prop_data->dim[2], prop_data->dim[3] };
                    gl::init_texture_3D(&volume_texture.id, dim[0], dim[1], dim[2], GL_R16F);
                    MEMCPY(volume_texture.dim, dim, sizeof(dim));
                    volume_texture.max_value = prop_data->max_value;
                }
                gl::set_texture_3D_data(volume_texture.id, 0, prop_data->values, GL_R32F);
            }
        }
    }

    void draw(ApplicationState* state) {
        if (!show_window) return;

        md_allocator_i* temp_arena = state->allocator.frame;
        md_vm_arena_temp_t temp_scope = md_vm_arena_temp_begin(temp_arena);
        defer { md_vm_arena_temp_end(temp_scope); };

        ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Density Volume", &show_window, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
            const ImVec2 button_size = {160, 0};
            bool volume_changed = false;

            if (ImGui::IsWindowFocused() && ImGui::IsKeyPressed(KEY_PLAY_PAUSE, false)) {
                state->animation.mode = state->animation.mode == PlaybackMode::Playing ? PlaybackMode::Stopped : PlaybackMode::Playing;
            }

            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("Property")) {
                    int64_t selected_index = -1;
                    int64_t candidate_count = 0;
                    for (int64_t i = 0; i < (int64_t)md_array_size(state->display_properties); ++i) {
                        DisplayProperty& dp = state->display_properties[i];
                        if (dp.type != DisplayProperty::Type_Volume) continue;
                        if (!state->timeline.filter.enabled && dp.partial_evaluation) {
                            continue;
                        }
                        ImPlot::ItemIcon(dp.color); ImGui::SameLine();
                        if (ImGui::Selectable(dp.label, dp.show_in_volume)) {
                            selected_index = i;
                            volume_changed = true;
                        }
                        if (ImGui::IsItemHovered()) {
                            script_visualize_payload(state, dp.vis_payload, -1, MD_SCRIPT_VISUALIZE_DEFAULT);
                            script_set_hovered_property(state, str_from_cstr(dp.label));
                        }
                        candidate_count += 1;
                    }

                    if (candidate_count == 0) {
                        ImGui::Text("No volume properties available.");
                    }

                    // Currently we only support viewing one volume at a time.
                    // This will probably change over time but not now.
                    if (selected_index != -1) {
                        for (int64_t i = 0; i < (int64_t)md_array_size(state->display_properties); ++i) {
                            if (selected_index == i) {
                                // Toggle bool
                                state->display_properties[i].show_in_volume = !state->display_properties[i].show_in_volume;
                            } else {
                                state->display_properties[i].show_in_volume = false;
                            }
                        }
                    }
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Render")) {
                    ImGui::Checkbox("Direct Volume Rendering", &dvr.enabled);
                    if (dvr.enabled) {
                        ImGui::Indent();
                        if (ImPlot::ColormapButton(ImPlot::GetColormapName(dvr.tf.colormap), button_size, dvr.tf.colormap)) {
                            ImGui::OpenPopup("Colormap Selector");
                        }
                        if (ImGui::BeginPopup("Colormap Selector")) {
                            for (int map = 4; map < ImPlot::GetColormapCount(); ++map) {
                                if (ImPlot::ColormapButton(ImPlot::GetColormapName(map), button_size, map)) {
                                    dvr.tf.colormap = map;
                                    dvr.tf.dirty = true;
                                    ImGui::CloseCurrentPopup();
                                }
                            }
                            ImGui::EndPopup();
                        }
                        if (ImGui::SliderFloat("TF Alpha Scaling", &dvr.tf.alpha_scale, 0.001f, 10.f, "%.3f", ImGuiSliderFlags_Logarithmic)) {
                            dvr.tf.dirty = true;
                        }
                        const float tf_min = 0.0f;
                        const float tf_max = 1000.0f;
                        ImGui::SliderFloat("TF Min Value", &dvr.tf.min_val, tf_min, dvr.tf.max_val, "%.3f", ImGuiSliderFlags_Logarithmic);
                        ImGui::SliderFloat("TF Max Value", &dvr.tf.max_val, dvr.tf.min_val, tf_max, "%.3f", ImGuiSliderFlags_Logarithmic);

                        ImGui::Unindent();
                    }
                    ImGui::Checkbox("Iso Surfaces", &iso.enabled);
                    if (iso.enabled) {
                        ImGui::Indent();
                        for (int i = 0; i < (int)iso.count; ++i) {
                            ImGui::PushID(i);
                            ImGui::SliderFloat("##Isovalue", &iso.values[i], 0.0f, 10.f, "%.3f", ImGuiSliderFlags_Logarithmic);
                            if (ImGui::IsItemDeactivatedAfterEdit()) {
                                // @TODO(Robin): Sort?
                            }
                            ImGui::SameLine();
                            ImGui::ColorEdit4Minimal("##Color", iso.colors[i].elem);
                            ImGui::SameLine();
                            if (ImGui::DeleteButton(ICON_FA_XMARK)) {
                                for (int j = i; j < (int)iso.count - 1; ++j) {
                                    iso.colors[j] = iso.colors[j+1];
                                    iso.values[j] = iso.values[j+1];
                                }
                                iso.count -= 1;
                            }
                            ImGui::PopID();
                        }
                        if ((iso.count < ARRAY_SIZE(iso.values)) && ImGui::Button("Add", button_size)) {
                            size_t idx = iso.count++;
                            iso.values[idx] = 0.1f;
                            iso.colors[idx] = { 0.2f, 0.1f, 0.9f, 1.0f };
                            // @TODO(Robin): Sort?
                        }
                            ImGui::SameLine();
                        if (ImGui::Button("Clear", button_size)) {
                            iso.count = 0;
                        }
                        ImGui::Unindent();
                    }
                    ImGui::EndMenu();
                }

                if (ImGui::BeginMenu("Clip planes")) {
                    ImGui::RangeSliderFloat("x", &clip_volume.min.x, &clip_volume.max.x, 0.0f, 1.0f);
                    ImGui::RangeSliderFloat("y", &clip_volume.min.y, &clip_volume.max.y, 0.0f, 1.0f);
                    ImGui::RangeSliderFloat("z", &clip_volume.min.z, &clip_volume.max.z, 0.0f, 1.0f);
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Show")) {
                    ImGui::Checkbox("Bounding Box", &show_bounding_box);
                    if (show_bounding_box) {
                        ImGui::Indent();
                        ImGui::ColorEdit4("Color", bounding_box_color.elem);
                        ImGui::Unindent();
                    }
                    ImGui::Checkbox("Reference Structure", &show_reference_structures);
                    if (show_reference_structures) {
                        ImGui::Indent();
                        ImGui::Checkbox("Show Superimposed Structures", &show_reference_ensemble);

                        if (ImGui::BeginCombo("type", representation_type_str[(int)rep.type])) {
                            for (int i = 0; i <= (int)RepresentationType::Cartoon; ++i) {
                                if (ImGui::Selectable(representation_type_str[i], (int)rep.type == i)) {
                                    rep.type = (RepresentationType)i;
                                    dirty_rep = true;
                                }
                            }
                            ImGui::EndCombo();
                        }

                        if (ImGui::BeginCombo("color", color_mapping_str[(int)rep.colormap])) {
                            for (int i = 0; i < (int)ColorMapping::Property; ++i) {
                                if (ImGui::Selectable(color_mapping_str[i], (int)rep.type == i)) {
                                    rep.colormap = (ColorMapping)i;
                                    dirty_rep = true;
                                }
                            }
                            ImGui::EndCombo();
                        }

                        if (rep.colormap == ColorMapping::Uniform) {
                            dirty_rep |= ImGui::ColorEdit4("color", rep.color.elem, ImGuiColorEditFlags_NoInputs);
                        }

                        if (rep.type == RepresentationType::SpaceFill || rep.type == RepresentationType::Licorice) {
                            dirty_rep |= ImGui::SliderFloat("scale", &rep.param[0], 0.1f, 2.f);
                        }
                        if (rep.type == RepresentationType::Ribbons) {
                            dirty_rep |= ImGui::SliderFloat("width", &rep.param[0], 0.1f, 2.f);
                            dirty_rep |= ImGui::SliderFloat("thickness", &rep.param[1], 0.1f, 2.f);
                        }
                        if (rep.type == RepresentationType::Cartoon) {
                            dirty_rep |= ImGui::SliderFloat("coil scale",  &rep.param[0], 0.1f, 3.f);
                            dirty_rep |= ImGui::SliderFloat("sheet scale", &rep.param[1], 0.1f, 3.f);
                            dirty_rep |= ImGui::SliderFloat("helix scale", &rep.param[2], 0.1f, 3.f);
                        }
                        ImGui::Unindent();
                    }
                    ImGui::Checkbox("Legend", &legend.enabled);
                    if (legend.enabled) {
                        ImGui::Indent();
                        const char* colormap_modes[] = {"Opaque", "Transparent", "Split"};
                        if (ImGui::BeginCombo("Colormap", colormap_modes[legend.colormap_mode])) {
                            for (int i = 0; i < IM_ARRAYSIZE(colormap_modes); ++i) {
                                if (ImGui::Selectable(colormap_modes[i])) {
                                    legend.colormap_mode = i;
                                }
                            }
                            ImGui::EndCombo();
                        }
                        ImGui::Checkbox("Use Checkerboard", &legend.checkerboard);
                        if (ImGui::IsItemHovered()) {
                            ImGui::SetTooltip("Use a checkerboard background for transparent parts in the legend.");
                        }
                        ImGui::Unindent();
                    }
                    ImGui::Checkbox("Coordinate System Widget", &show_coordinate_system_widget);
                    ImGui::EndMenu();
                }

                ImGui::EndMenuBar();
            }

            update(state);

            // Animate camera towards targets
            camera_animate(&camera, target, state->app.timing.delta_s);

            const ImVec2 canvas_sz = ImMax(ImGui::GetContentRegionAvail(), ImVec2(50.0f, 50.0f));   // Resize canvas to what's available

            int width  = (int)(canvas_sz.x * ImGui::GetIO().DisplayFramebufferScale.x);
            int height = (int)(canvas_sz.y * ImGui::GetIO().DisplayFramebufferScale.y);
            if ((int)gbuf.width != width || (int)gbuf.height != height) {
                init_gbuffer(&gbuf, width, height);
            }

            // @NOTE: We cannot use the *standard* interaction surface here
            // Since we render a bunch of structures, each with its own model matrix
            // This means that the picking logic needs to be aware of all those model matrices, which is a bit of a pain to implement in a robust way.

            // This will catch our interactions
            ImGui::InvisibleButton("canvas", canvas_sz, ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_AllowOverlap);
            const ImVec2 canvas_p0 = ImGui::GetItemRectMin();
            const ImVec2 canvas_p1 = ImGui::GetItemRectMax();

            // Draw border and background color
            ImGuiIO& io = ImGui::GetIO();

            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->AddImage((ImTextureID)(intptr_t)gbuf.tex.transparency, canvas_p0, canvas_p1, { 0,1 }, { 1,0 });
            draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));

            if (dvr.enabled && legend.enabled) {
                ImVec2 canvas_ext = canvas_p1 - canvas_p0;
                ImVec2 cmap_ext = {MIN(canvas_ext.x * 0.5f, 250.0f), MIN(canvas_ext.y * 0.25f, 30.0f)};
                ImVec2 cmap_pad = {10, 10};
                ImVec2 cmap_pos = canvas_p1 - ImVec2(cmap_ext.x, cmap_ext.y) - cmap_pad;
                ImPlotColormap cmap = dvr.tf.colormap;
                ImPlotContext& gp = *ImPlot::GetCurrentContext();
                ImU32 checker_bg = IM_COL32(255, 255, 255, 255);
                ImU32 checker_fg = IM_COL32(128, 128, 128, 255);
                float checker_size = 8.0f;
                ImVec2 checker_offset = ImVec2(0,0);

                int mode = legend.colormap_mode;

                ImVec2 opaque_scl = ImVec2(1,1);
                ImVec2 transp_scl = ImVec2(0,0);

                if (mode == LegendColorMapMode_Split) {
                    opaque_scl = ImVec2(1.0f, 0.5f);
                    transp_scl = ImVec2(0.0f, 0.5f);
                }

                ImRect opaque_rect = ImRect(cmap_pos, cmap_pos + cmap_ext * opaque_scl);
                ImRect transp_rect = ImRect(cmap_pos + cmap_ext * transp_scl, cmap_pos + cmap_ext);
            
                // Opaque
                if (mode == LegendColorMapMode_Opaque || mode == LegendColorMapMode_Split) {
                    ImPlot::RenderColorBar(gp.ColormapData.GetKeys(cmap),gp.ColormapData.GetKeyCount(cmap),*draw_list,opaque_rect,false,false,!gp.ColormapData.IsQual(cmap));
                }
            
                if (mode == LegendColorMapMode_Transparent || mode == LegendColorMapMode_Split) {
                    if (legend.checkerboard) {
                        // Checkerboard
                        ImGui::DrawCheckerboard(draw_list, transp_rect.Min, transp_rect.Max, checker_bg, checker_fg, checker_size, checker_offset);
                    }
                    // Transparent
                    draw_list->AddImage((ImTextureID)(intptr_t)dvr.tf.id, transp_rect.Min, transp_rect.Max);
                }
            
                // Boarder
                draw_list->AddRect(cmap_pos, cmap_pos + cmap_ext, IM_COL32(0, 0, 0, 255), 0.0f, 0, 0.1f);
            }

            const bool is_hovered = ImGui::IsItemHovered();
            const bool is_active = ImGui::IsItemActive();
            const ImVec2 origin(canvas_p0.x, canvas_p0.y);  // Lock scrolled origin
            const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

            bool reset_hard = false;
            if (volume_changed) {
                static bool first_time = true;
                if (first_time) {
                    reset_hard = true;
                    first_time = false;
                }
            }
            bool reset_view = reset_hard;
            if (is_hovered) {
                if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                    reset_view = true;
                }
            }

            if (reset_view) {
                vec3_t min_ext = vec3_from_vec4(model_mat * vec4_set(0,0,0,1));
                vec3_t max_ext = vec3_from_vec4(model_mat * vec4_set(1,1,1,1));
                target = compute_optimal_view(min_ext, max_ext);

                if (reset_hard) {
				    camera = target;
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
                    .fov_y = camera.fov_y,
                };
                camera_controller_trackball(&target, input, param);
            }

            if (show_coordinate_system_widget) {
                ImVec2 win_size = ImGui::GetWindowSize();
                float  ext = MIN(win_size.x, win_size.y) * 0.2f;
                float  pad = 0.1f * ext;
			    ImVec2 size = { ext, ext };

			    ImGui::SetCursorScreenPos(ImVec2(pad, win_size.y - ext - pad));

                quat_t out_orientation = target.orientation;
                if (ImGui::CoordinateSystemWidget(&out_orientation, camera.orientation, size)) {
                    const vec3_t look_at = camera_get_look_at(target);
                    target.orientation = quat_normalize(out_orientation);
                    target.position = camera_position_from_look_at(look_at, target.orientation, target.distance);
                }
            }

            mat4_t view_mat = camera_world_to_view_matrix(camera);
            mat4_t proj_mat = camera_view_to_clip_matrix_persp(camera, (float)canvas_sz.x / (float)canvas_sz.y);
            mat4_t inv_proj_mat = camera_clip_to_view_matrix_persp(camera, (float)canvas_sz.x / (float)canvas_sz.y);

            PUSH_GPU_SECTION("RENDER DENSITY VOLUME");
            clear_gbuffer(&gbuf);

            const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
                GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_TRANSPARENCY };

            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);

            glEnable(GL_DEPTH_TEST);
            glDepthMask(GL_TRUE);
            glDepthFunc(GL_LESS);

            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
            glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
            glViewport(0, 0, gbuf.width, gbuf.height);
            glScissor(0, 0,  gbuf.width, gbuf.height);

            int64_t selected_property = -1;
            for (size_t i = 0; i < md_array_size(state->display_properties); ++i) {
                const DisplayProperty& dp = state->display_properties[i];
                if (dp.type == DisplayProperty::Type_Volume && dp.show_in_volume) {
                    selected_property = i;
                    break;
                }
            }

            size_t num_reps = md_array_size(gl_reps);
            if (selected_property > -1 && show_reference_structures && num_reps > 0) {
                if (!show_reference_ensemble) {
                    num_reps = 1;
                }

                md_gl_draw_op_t* draw_ops = md_array_create(md_gl_draw_op_t, num_reps, temp_arena);

                md_gl_draw_op_t op = {};
                op.type = (md_gl_rep_type_t)rep.type;
                MEMCPY(&op.args, rep.param, sizeof(op.args));
                if (op.type == MD_GL_REP_BALL_AND_STICK) {
                    op.args.ball_and_stick.color_mode = MD_GL_BOND_MODE_NEAREST;
                } else if (op.type == MD_GL_REP_LICORICE) {
                    op.args.licorice.color_mode = MD_GL_BOND_MODE_NEAREST;
                }

                for (size_t i = 0; i < num_reps; ++i) {
                    op.rep = gl_reps[i];
                    op.model_matrix = &rep_model_mats[i].elem[0][0];
                    draw_ops[i] = op;
                }

                md_gl_draw_args_t draw_args = {
                    .shaders = state->gl.shaders,
                    .draw_operations = {
                        .count = md_array_size(draw_ops),
                        .ops = draw_ops
                    },
                    .view_transform = {
                        .view_matrix = &view_mat.elem[0][0],
                        .proj_matrix = &proj_mat.elem[0][0],
                    },
                };

                md_gl_draw(&draw_args);

                if (is_hovered) {
                    mat4_t inv_MVP = mat4_mul(inv_proj_mat, mat4_transpose(view_mat));
                    const vec2_t coord = {mouse_pos_in_canvas.x, (float)gbuf.height - mouse_pos_in_canvas.y};
                    // @TODO: Reimplement picking here using the new picking api 
                    //extract_picking_data(data->picking, gbuf, coord, inv_MVP);
                    //if (data->picking.idx != INVALID_PICKING_IDX) {
                    //    draw_info_window(*data, data->picking.idx);
                    //}
                }
                glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
            }

            if (show_density_volume) {
                if (model_mat != mat4_t{ 0 }) {
                    volume::RenderDesc vol_desc = {
                        .render_target = {
                            .depth  = gbuf.tex.depth,
                            .color  = gbuf.tex.transparency,
                            .width  = gbuf.width,
                            .height = gbuf.height,
                        },
                        .texture = {
                            .density_volume = volume_texture.id,
                            .transfer_function = dvr.tf.id,
                        },
                        .matrix = {
                            .model = model_mat,
                            .view = view_mat,
                            .proj = proj_mat,
                            .inv_proj = inv_proj_mat,
                        },
                        .clip_volume = {
                            .min = clip_volume.min,
                            .max = clip_volume.max,
                        },
                        .iso = {
                            .enabled = iso.enabled,
                            .count = iso.count,
                            .values = iso.values,
                            .colors = iso.colors,
                        },
                        .dvr = {
                            .enabled = dvr.enabled,
                            .min_tf_value = dvr.tf.min_val,
                            .max_tf_value = dvr.tf.max_val,
                        },
                        .shading = {
                            .env_radiance = state->visuals.background.color * state->visuals.background.intensity,
                            .roughness = 0.3f,
                            .dir_radiance = {10,10,10},
                            .ior = 1.5f,
                            .exposure = state->visuals.tonemapping.exposure,
                            .gamma = state->visuals.tonemapping.gamma,
                        },
                        .voxel_spacing = voxel_spacing
                    
                    };
                    volume::render_volume(vol_desc);
                }
            }

            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
            glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
            glViewport(0, 0, gbuf.width, gbuf.height);
            glScissor(0, 0, gbuf.width, gbuf.height);

            if (show_bounding_box) {
                glEnable(GL_DEPTH_TEST);
                glDepthMask(GL_TRUE);
                glEnable(GL_BLEND);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

                if (model_mat != mat4_t{0}) {
                    immediate::set_model_view_matrix(mat4_mul(view_mat, model_mat));
                    immediate::set_proj_matrix(proj_mat);

                    uint32_t box_color = convert_color(bounding_box_color);
                    uint32_t clip_color = convert_color(clip_volume_color);
                    immediate::draw_box_wireframe({0,0,0}, {1,1,1}, box_color);
                    immediate::draw_box_wireframe(clip_volume.min, clip_volume.max, clip_color);

                    immediate::render();
                }
                glDisable(GL_BLEND);
            }

            PUSH_GPU_SECTION("Postprocessing")
            postprocessing::Descriptor postprocess_desc = {
                .background = {
                    .color = state->visuals.background.color * state->visuals.background.intensity,
                },
                .tonemapping = {
                    .enabled = state->visuals.tonemapping.enabled,
                    .mode = state->visuals.tonemapping.tonemapper,
                    .exposure = state->visuals.tonemapping.exposure,
                    .gamma = state->visuals.tonemapping.gamma,
                },
                .ambient_occlusion = {
                    .enabled = false,
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
                    .depth = gbuf.tex.depth,
                    .color = gbuf.tex.color,
                    .normal = gbuf.tex.normal,
                    .velocity = gbuf.tex.velocity,
                    .transparency = gbuf.tex.transparency,
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
                    .near = camera.near_plane,
                    .far = camera.far_plane,
                },
                .resolution = {canvas_sz.x, canvas_sz.y},
                .fov_y = camera.fov_y,
            };

            postprocessing::shade_and_postprocess(postprocess_desc, view_param);
            POP_GPU_SECTION()

            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
            glDrawBuffer(GL_BACK);

            POP_GPU_SECTION();
        }

        ImGui::End();
    }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                // Initialize component
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState* state = (ApplicationState*)e.payload;
                arena = md_arena_allocator_create(state->allocator.persistent, MEGABYTES(1));
                picking_surface_init(&picking_surface, PickingSource_DensityVolume);
                break;
            }
            case viamd::EventType_ViamdShutdown: {
                // Cleanup component
                md_arena_allocator_destroy(arena);
                break;
            }
            case viamd::EventType_ViamdSystemInit: {
                break;
            }
            case viamd::EventType_ViamdSystemFree: {
                md_array_shrink(gl_reps, 0);
                md_array_shrink(rep_model_mats, 0);
                model_mat = {0};
                break;
            }
            case viamd::EventType_ViamdFrameTick: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState* state = (ApplicationState*)e.payload;
                draw(state);
                break;
            }
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("Density Volume", &show_window);
                break;
            default:
                break;
            }
        }
    }
};

static DensityVolume instance;