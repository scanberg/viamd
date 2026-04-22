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

constexpr uint64_t interaction_surface_density_vol = HASH_STR_LIT64("interaction surface density volume");

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

    struct DensityVolumeRepresentation {
        mat4_t model_mat = mat4_ident();
        md_array(md_atom_idx_t) atom_indices = nullptr;
        md_gl_rep_t gl_rep = {};
        bool enabled = false;
    };

    md_array(DensityVolumeRepresentation) reps = nullptr;

    // Base model matrix to position the volume correctly in the world.
    mat4_t model_mat = mat4_ident();

    GBuffer gbuf = {};
    PickingSurface picking_surface = {};
    Camera camera = {};
	ViewTransform target = {};
    ViewTransform default_view = {};

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

        bool reset_view = false;
        static int64_t s_selected_property = 0;
        if (s_selected_property != selected_property) {
            if (s_selected_property == -1) {
                reset_view = true;
            }
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
                        vec3_t min_aabb = vec3_set1(-s);
                        vec3_t max_aabb = vec3_set1( s);
                        model_mat = volume::compute_model_to_world_matrix(min_aabb, max_aabb);
                        voxel_spacing = vec3_t{2*s / prop_data->dim[1], 2*s / prop_data->dim[2], 2*s / prop_data->dim[3]};
                        default_view = compute_optimal_view(min_aabb, max_aabb);
                        if (reset_view) {
                            target = default_view;
                            camera = default_view;
                        }
                    }
                    num_reps = md_array_size(vis.sdf.structures);
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

                // We need to limit this for performance reasons
                num_reps = MIN(num_reps, 100);

                const size_t old_size = md_array_size(reps);
                if (reps) {
                    // Only free superflous entries
                    for (size_t i = num_reps; i < old_size; ++i) {
                        md_gl_rep_destroy(reps[i].gl_rep);
                        md_array_free(reps[i].atom_indices, arena);
                    }
                }
                md_array_resize(reps, num_reps, arena);

                for (size_t i = old_size; i < num_reps; ++i) {
                    // Only init new entries
                    reps[i].gl_rep = md_gl_rep_create(state->mold.gl_mol);
                    reps[i].atom_indices = nullptr;
                }

                for (size_t i = 0; i < num_reps; ++i) {
                    filter_colors(colors, num_atoms, &vis.sdf.structures[i]);
                    md_gl_rep_set_atom_colors(reps[i].gl_rep, 0, (uint32_t)num_atoms, colors, 0);
                    reps[i].model_mat = vis.sdf.matrices[i];
                    size_t popcount = md_bitfield_popcount(&vis.sdf.structures[i]);
                    md_array_resize(reps[i].atom_indices, popcount, arena);
                    md_bitfield_iter_extract_indices(reps[i].atom_indices, popcount, md_bitfield_iter_create(&vis.sdf.structures[i]));
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

            ViewTransform reset_transform = default_view;

            const float aspect_ratio = canvas_sz.x / canvas_sz.y;
            const mat4_t world_to_view = camera_world_to_view_matrix(camera);
            const mat4_t view_to_world = camera_view_to_world_matrix(camera);

            const mat4_t view_to_clip  = camera_view_to_clip_matrix_persp(camera, aspect_ratio);
            const mat4_t clip_to_view  = camera_clip_to_view_matrix_persp(camera, aspect_ratio);

            const mat4_t world_to_clip = view_to_clip * world_to_view;
            const mat4_t clip_to_world = view_to_world * clip_to_view;

            const TrackballControllerParam trackball_param = {
                .min_distance = 1.0,
                .max_distance = 1000.0,
            };

            const InteractionSurfaceFlags surface_flags = InteractionSurfaceFlags_NoRegionSelect;
            InteractionSurfaceState surface_state = interaction_surface(interaction_surface_density_vol, vec_cast(canvas_sz), surface_flags);
            const ImVec2 canvas_p0 = ImGui::GetItemRectMin();
            const ImVec2 canvas_p1 = ImGui::GetItemRectMax();
            const ImRect canvas_rect = ImRect(canvas_p0, canvas_p1);

            if (surface_state.hovered) {
                InteractionSurfaceHitArgs args = {
                    .picking_surface = &picking_surface,
                    .picking_handler = state->picking_handler,
                    .fbo = gbuf.fbo,
                    .width = gbuf.width,
                    .height = gbuf.height,
                    .clip_to_world = clip_to_world,
                };

                PickingHit hit = {};
                interaction_surface_hit_extract(&hit, surface_state, args);

                if (hit.depth < 1.0f) {
                    reset_transform.distance = target.distance;
                    reset_transform.orientation = camera.orientation;
                    reset_transform.position = hit.world_pos + camera.orientation * vec3_set(0, 0, target.distance);                    
                }

                InteractionSurfaceEvent event = {};
                interaction_surface_event_extract(&event, surface_state, hit);

                event.clip_to_world = clip_to_world;
                event.world_to_clip = world_to_clip;

                if (event.kind == InteractionSurfaceEventKind::RegionSelect) {
                    /*
                    const md_bitfield_t* candidate_mask = &state.representation.visibility_mask;
                    if (event.selection_mode == InteractionSelectionMode::Remove) {
                        // When removing, only consider currently selected atoms as candidates for region selection
                        candidate_mask = &state->selection.selection_mask;
                    }
                    point_set_region_mask_compute(&state->selection.highlight_mask,
                        state->mold.sys.atom.x, state->mold.sys.atom.y, state->mold.sys.atom.z, state->mold.sys.atom.count,
                        candidate_mask, world_to_clip, event.region_min, event.region_max, event.surface_size);

                    grow_mask_by_selection_granularity(&state->selection.highlight_mask, state->selection.granularity, state->mold.sys);
                    if (event.region_phase == InteractionSurfaceEventPhase::Commit) {
                        // Merge highlight into selection
                        if (event.selection_mode == InteractionSelectionMode::Append) {
                            md_bitfield_or_inplace(&state->selection.selection_mask, &state->selection.highlight_mask);
                        }
                        else if (event.selection_mode == InteractionSelectionMode::Remove) {
                            md_bitfield_andnot_inplace(&state->selection.selection_mask, &state->selection.highlight_mask);
                        }
                        md_bitfield_clear(&state->selection.highlight_mask);
                    }
                    */
                } else {
                    viamd::event_system_broadcast_event(viamd::EventType_ViamdInteractionSurface, viamd::EventPayloadType_InteractionSurfaceEvent, &event);
                }
            }

            InteractionSurfaceViewTransformArgs view_args = {
                .camera = camera,
                .trackball_param = trackball_param,
                .reset_transform = reset_transform,
            };

            interaction_surface_view_transform_apply(&target, surface_state, view_args);

            // Draw border and background color
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

            if (show_coordinate_system_widget) {
                float  ext = MIN(canvas_rect.GetWidth(), canvas_rect.GetHeight()) * 0.2f;
                float  pad = 0.1f * ext;
			    ImVec2 size = { ext, ext };

			    ImGui::SetCursorScreenPos(ImVec2(canvas_p0.x + pad, canvas_p1.y - ext - pad));

                quat_t out_orientation = target.orientation;
                if (ImGui::CoordinateSystemWidget(&out_orientation, camera.orientation, size)) {
                    const vec3_t look_at = camera_get_look_at(target);
                    target.orientation = quat_normalize(out_orientation);
                    target.position = camera_position_from_look_at(look_at, target.orientation, target.distance);
                }
            }

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

            size_t num_reps = md_array_size(reps);
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
                    op.rep = reps[i].gl_rep;
                    op.model_matrix = &reps[i].model_mat.elem[0][0];
                    draw_ops[i] = op;
                }

                md_gl_draw_args_t draw_args = {
                    .shaders = state->gl.shaders,
                    .draw_operations = {
                        .count = md_array_size(draw_ops),
                        .ops = draw_ops
                    },
                    .view_transform = {
                        .view_matrix = &world_to_view.elem[0][0],
                        .proj_matrix = &view_to_clip.elem[0][0],
                    },
                    .picking_offset = {
                        .atom_base = state->picking_range_atom.beg,
                        .bond_base = state->picking_range_bond.beg,
                    }
                };

                md_gl_draw(&draw_args);

                glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
            }

            if (show_density_volume) {
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
                        .view = world_to_view,
                        .proj = view_to_clip,
                        .inv_proj = clip_to_view,
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

            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
            glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
            glViewport(0, 0, gbuf.width, gbuf.height);
            glScissor(0, 0, gbuf.width, gbuf.height);

            if (show_bounding_box && selected_property != -1) {
                glEnable(GL_DEPTH_TEST);
                glDepthMask(GL_TRUE);
                glEnable(GL_BLEND);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

                immediate::set_model_view_matrix(mat4_mul(world_to_view, model_mat));
                immediate::set_proj_matrix(view_to_clip);

                uint32_t box_color = convert_color(bounding_box_color);
                uint32_t clip_color = convert_color(clip_volume_color);
                immediate::draw_box_wireframe({0,0,0}, {1,1,1}, box_color);
                immediate::draw_box_wireframe(clip_volume.min, clip_volume.max, clip_color);

                immediate::render();

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
                        .view = world_to_view,
                        .proj = view_to_clip,
                        .norm = world_to_view,
                    },
                    .inv = {
                        .proj = clip_to_view,
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
                picking_surface_init(&picking_surface, interaction_surface_density_vol);
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
                md_array_shrink(reps, 0);
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