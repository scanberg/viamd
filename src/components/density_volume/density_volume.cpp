#include <core/md_arena_allocator.h>

#include <viamd_event.h>
#include <viamd.h>
#include <serialization_utils.h>

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

    // The volume property shown, by path; empty for none. One at a time.
    SeriesKey volume_key = {};

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

    // ## Workspace: which volume is shown, and how

    void reset_workspace_settings() {
        volume_key = {};
        dvr.enabled = true;
        dvr.tf.alpha_scale = 1.0f;
        dvr.tf.colormap = DEFAULT_COLORMAP;
        dvr.tf.min_val = 0.0f;
        dvr.tf.max_val = 1.0f;
        dvr.tf.dirty = true;
        iso.enabled = false;
        iso.count = 0;
        MEMSET(iso.values, 0, sizeof(iso.values));
        MEMSET(iso.colors, 0, sizeof(iso.colors));
        clip_volume.min = {0, 0, 0};
        clip_volume.max = {1, 1, 1};
        legend.enabled = true;
        legend.checkerboard = true;
        legend.colormap_mode = LegendColorMapMode_Split;
        resolution_scale = 2.0f;
        clip_volume_color = {1,0,0,1};
        bounding_box_color = {0,0,0,1};
        show_bounding_box = true;
        show_reference_structures = true;
        show_reference_ensemble = false;
        show_coordinate_system_widget = true;
        rep.type = RepresentationType::BallAndStick;
        rep.colormap = ColorMapping::Type;
        rep.param[0] = rep.param[1] = rep.param[2] = rep.param[3] = 1.0f;
        rep.color = {1,1,1,1};
        dirty_rep = true;
        dirty_vol = true;
    }

    void serialize(viamd::serialization_state_t& state) {
        viamd::write_section_header(state, STR_LIT("DensityVolume"));
        if (volume_key.path[0] != '\0') {
            viamd::write_str(state, STR_LIT("PropertySource"), str_from_cstr(series_source_name(volume_key.source)));
            viamd::write_str(state, STR_LIT("PropertyPath"), str_from_cstr(volume_key.path));
        }
        viamd::write_bool(state, STR_LIT("DvrEnabled"), dvr.enabled);
        viamd::write_int (state, STR_LIT("DvrColormap"), dvr.tf.colormap);
        viamd::write_flt (state, STR_LIT("DvrAlphaScale"), dvr.tf.alpha_scale);
        const float tf_range[2] = { dvr.tf.min_val, dvr.tf.max_val };
        viamd::write_flt_vec(state, STR_LIT("DvrRange"), tf_range, 2);
        viamd::write_bool(state, STR_LIT("IsoEnabled"), iso.enabled);
        viamd::write_int (state, STR_LIT("IsoCount"), (int64_t)iso.count);
        for (size_t i = 0; i < iso.count && i < ARRAY_SIZE(iso.values); ++i) {
            char key[32];
            snprintf(key, sizeof(key), "IsoValue%zu", i);
            viamd::write_flt(state, str_from_cstr(key), iso.values[i]);
            snprintf(key, sizeof(key), "IsoColor%zu", i);
            viamd::write_vec4(state, str_from_cstr(key), iso.colors[i]);
        }
        viamd::write_vec3(state, STR_LIT("ClipMin"), clip_volume.min);
        viamd::write_vec3(state, STR_LIT("ClipMax"), clip_volume.max);
        viamd::write_bool(state, STR_LIT("LegendEnabled"), legend.enabled);
        viamd::write_bool(state, STR_LIT("LegendCheckerboard"), legend.checkerboard);
        viamd::write_int (state, STR_LIT("LegendColormapMode"), legend.colormap_mode);
        viamd::write_flt (state, STR_LIT("ResolutionScale"), resolution_scale);
        viamd::write_vec4(state, STR_LIT("ClipVolumeColor"), clip_volume_color);
        viamd::write_vec4(state, STR_LIT("BoundingBoxColor"), bounding_box_color);
        viamd::write_bool(state, STR_LIT("ShowBoundingBox"), show_bounding_box);
        viamd::write_bool(state, STR_LIT("ShowReferenceStructures"), show_reference_structures);
        viamd::write_bool(state, STR_LIT("ShowReferenceEnsemble"), show_reference_ensemble);
        viamd::write_bool(state, STR_LIT("ShowCoordinateSystem"), show_coordinate_system_widget);
        viamd::write_int (state, STR_LIT("RepType"), (int)rep.type);
        viamd::write_int (state, STR_LIT("RepColorMapping"), (int)rep.colormap);
        viamd::write_flt_vec(state, STR_LIT("RepParam"), rep.param, 4);
        viamd::write_vec4(state, STR_LIT("RepColor"), rep.color);
    }

    void deserialize(viamd::deserialization_state_t& state) {
        SeriesSource source = SeriesSource_Script;
        str_t ident, arg;
        while (viamd::next_entry(ident, arg, state)) {
            if (str_eq_cstr(ident, "PropertySource")) {
                str_t name;
                viamd::extract_str(name, arg);
                series_source_from_name(&source, name);
            } else if (str_eq_cstr(ident, "PropertyPath")) {
                str_t path;
                viamd::extract_str(path, arg);
                volume_key = series_key(source, path);
            } else if (str_eq_cstr(ident, "DvrEnabled")) {
                viamd::extract_bool(dvr.enabled, arg);
            } else if (str_eq_cstr(ident, "DvrColormap")) {
                viamd::extract_int(dvr.tf.colormap, arg);
            } else if (str_eq_cstr(ident, "DvrAlphaScale")) {
                viamd::extract_flt(dvr.tf.alpha_scale, arg);
            } else if (str_eq_cstr(ident, "DvrRange")) {
                float r[2];
                if (viamd::extract_flt_vec(r, 2, arg)) { dvr.tf.min_val = r[0]; dvr.tf.max_val = r[1]; }
            } else if (str_eq_cstr(ident, "IsoEnabled")) {
                viamd::extract_bool(iso.enabled, arg);
            } else if (str_eq_cstr(ident, "IsoCount")) {
                int n;
                if (viamd::extract_int(n, arg)) iso.count = (size_t)CLAMP(n, 0, (int)ARRAY_SIZE(iso.values));
            } else if (str_eq_cstr(ident, "ClipMin")) {
                viamd::extract_vec3(clip_volume.min, arg);
            } else if (str_eq_cstr(ident, "ClipMax")) {
                viamd::extract_vec3(clip_volume.max, arg);
            } else if (str_eq_cstr(ident, "LegendEnabled")) {
                viamd::extract_bool(legend.enabled, arg);
            } else if (str_eq_cstr(ident, "LegendCheckerboard")) {
                viamd::extract_bool(legend.checkerboard, arg);
            } else if (str_eq_cstr(ident, "LegendColormapMode")) {
                viamd::extract_int(legend.colormap_mode, arg);
            } else if (str_eq_cstr(ident, "ResolutionScale")) {
                viamd::extract_flt(resolution_scale, arg);
            } else if (str_eq_cstr(ident, "ClipVolumeColor")) {
                viamd::extract_vec4(clip_volume_color, arg);
            } else if (str_eq_cstr(ident, "BoundingBoxColor")) {
                viamd::extract_vec4(bounding_box_color, arg);
            } else if (str_eq_cstr(ident, "ShowBoundingBox")) {
                viamd::extract_bool(show_bounding_box, arg);
            } else if (str_eq_cstr(ident, "ShowReferenceStructures")) {
                viamd::extract_bool(show_reference_structures, arg);
            } else if (str_eq_cstr(ident, "ShowReferenceEnsemble")) {
                viamd::extract_bool(show_reference_ensemble, arg);
            } else if (str_eq_cstr(ident, "ShowCoordinateSystem")) {
                viamd::extract_bool(show_coordinate_system_widget, arg);
            } else if (str_eq_cstr(ident, "RepType")) {
                viamd::extract_enum(rep.type, arg, (int)RepresentationType::Count);
            } else if (str_eq_cstr(ident, "RepColorMapping")) {
                viamd::extract_enum(rep.colormap, arg, (int)ColorMapping::Count);
            } else if (str_eq_cstr(ident, "RepParam")) {
                viamd::extract_flt_vec(rep.param, 4, arg);
            } else if (str_eq_cstr(ident, "RepColor")) {
                viamd::extract_vec4(rep.color, arg);
            } else {
                for (size_t i = 0; i < ARRAY_SIZE(iso.values); ++i) {
                    char key[32];
                    snprintf(key, sizeof(key), "IsoValue%zu", i);
                    if (str_eq_cstr(ident, key)) { viamd::extract_flt(iso.values[i], arg); break; }
                    snprintf(key, sizeof(key), "IsoColor%zu", i);
                    if (str_eq_cstr(ident, key)) { viamd::extract_vec4(iso.colors[i], arg); break; }
                }
            }
        }
        // Workspaces from before the two were exclusive can have both on: DVR wins, as it did by default
        if (dvr.enabled && iso.enabled) {
            iso.enabled = false;
        }
        dvr.tf.dirty = true;
        dirty_rep = true;
        dirty_vol = true;
    }

    void update(ApplicationState* state) {
        md_allocator_i* temp_arena = state->allocator.frame;
        md_temp_scope_t temp_scope = md_temp_begin_in(temp_arena);
        defer { md_temp_end(temp_scope); };

        if (dvr.tf.dirty) {
            dvr.tf.dirty = false;
            // Update colormap texture
            volume::compute_transfer_function_texture_simple(&dvr.tf.id, dvr.tf.colormap, dvr.tf.alpha_scale);
        }

        // Resolved by path every frame: the evaluation's table is rebuilt when the script recompiles
        const bool selected = volume_key.path[0] != '\0';
        SeriesVolumeView vol_view = {};
        const bool resolved = selected && series_resolve_volume(&vol_view, state, volume_key);

        const md_attribute_t* prop_attr = resolved ? vol_view.attr : NULL;
        const md_script_vis_payload_o* vis_payload = vol_view.vis_payload;
        const uint64_t data_version = vol_view.version;

        bool reset_view = false;
        static SeriesKey s_volume_key = {};
        static bool s_first = true;
        if (s_first || !series_key_equal(s_volume_key, volume_key)) {
            if (!s_first && s_volume_key.path[0] == '\0' && selected) {
                reset_view = true;
            }
            s_first = false;
            s_volume_key = volume_key;
            dirty_vol = true;
            dirty_rep = true;
        }
        show_density_volume = selected;

        static uint64_t s_script_fingerprint = 0;
        if (s_script_fingerprint != md_script_ir_fingerprint(state->script.eval_ir)) {
            s_script_fingerprint = md_script_ir_fingerprint(state->script.eval_ir);
            dirty_vol = true;
            dirty_rep = true;
        }

        static uint64_t s_data_version = 0;
        if (s_data_version != data_version) {
            s_data_version = data_version;
            dirty_vol = true;
        }

        static double s_frame = 0;
        if (s_frame != state->animation.frame) {
            s_frame = state->animation.frame;
            dirty_rep = true;
        }

        if (dirty_rep) {
            if (prop_attr && vis_payload) {
                dirty_rep = false;
                size_t num_reps = 0;
                bool result = false;
                md_script_vis_t vis = {};

                if (md_script_ir_valid(state->script.eval_ir)) {
                    md_script_vis_init(&vis, state->allocator.frame);
                    md_script_vis_ctx_t ctx = {
                        .ir = state->script.eval_ir,
                        .sys = &state->mold.sys,
                        .state = &state->mold.state
                    };
                    result = md_script_vis_eval_payload(&vis, vis_payload, 0, &ctx, MD_SCRIPT_VISUALIZE_SDF);
                }

                if (result) {
                    if (vis.sdf.extent) {
                        const float s = vis.sdf.extent;
                        vec3_t min_aabb = vec3_set1(-s);
                        vec3_t max_aabb = vec3_set1( s);
                        model_mat = volume::compute_model_to_world_matrix(min_aabb, max_aabb);
                        const uint32_t* dim = prop_attr->format.shape;
                        voxel_spacing = vec3_t{2*s / dim[0], 2*s / dim[1], 2*s / dim[2]};
                        default_view = compute_optimal_view((min_aabb + max_aabb) * 0.5f, (max_aabb - min_aabb) * 0.5f);
                        if (reset_view) {
                            target = default_view;
                            camera = default_view;
                        }
                    }
                    num_reps = md_array_size(vis.sdf.structures);
                }

                const md_system_t& sys = state->mold.sys;
			    const size_t num_atoms = md_system_atom_count(&sys);
                md_temp_scope_t temp = md_temp_begin_in(state->allocator.frame);
                defer { md_temp_end(temp); };
                const size_t num_bytes = sizeof(uint32_t) * num_atoms;
                uint32_t* colors = (uint32_t*)md_temp_alloc(temp, num_bytes);

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
                    color_atoms_secondary_structure(colors, num_atoms, sys, displayed_secondary_structure(state));
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
            if (prop_attr) {
                dirty_vol = false;
                if (!volume_texture.id) {
                    int dim[3] = { (int)prop_attr->format.shape[0], (int)prop_attr->format.shape[1], (int)prop_attr->format.shape[2] };
                    gl::init_texture_3D(&volume_texture.id, dim[0], dim[1], dim[2], GL_R16F);
                    MEMCPY(volume_texture.dim, dim, sizeof(dim));
                }
                gl::set_texture_3D_data(volume_texture.id, 0, prop_attr->data, GL_R32F);
                volume::notify_data_changed(volume_texture.id);
            }
        }
    }

    void draw(ApplicationState* state) {
        if (!show_window) return;

        md_allocator_i* temp_arena = state->allocator.frame;
        md_temp_scope_t temp_scope = md_temp_begin_in(temp_arena);
        defer { md_temp_end(temp_scope); };

        ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Density Volume", &show_window, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
            const ImVec2 button_size = {160, 0};

            if (ImGui::IsWindowFocused() && ImGui::IsKeyPressed(KEY_PLAY_PAUSE, false)) {
                state->animation.mode = state->animation.mode == PlaybackMode::Playing ? PlaybackMode::Stopped : PlaybackMode::Playing;
            }

            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("Property")) {
                    // One volume at a time: picking one replaces the other, picking it again hides it
                    int candidate_count = 0;
                    auto list = [&](SeriesSource source) {
                        series_for_each_script_property(state, source, MD_SCRIPT_PROPERTY_FLAG_VOLUME, [&](const SeriesKey& key) {
                            char label[96];
                            series_label(label, sizeof(label), state, key);
                            const bool is_selected = series_key_equal(key, volume_key);
                            ImGui::PushID(key.source);
                            ImPlot::ItemIcon(series_default_color(state, key)); ImGui::SameLine();
                            if (ImGui::Selectable(label, is_selected)) {
                                volume_key = is_selected ? SeriesKey{} : key;
                            }
                            if (ImGui::IsItemHovered()) {
                                const str_t ident = series_script_ident(key);
                                const md_script_vis_payload_o* vis = state->script.eval_ir ? md_script_ir_property_vis_payload(state->script.eval_ir, ident) : nullptr;
                                script_visualize_payload(state, vis, -1, MD_SCRIPT_VISUALIZE_DEFAULT);
                                script_set_hovered_property(state, ident);
                            }
                            ImGui::PopID();
                            candidate_count += 1;
                        });
                    };
                    list(SeriesSource_Script);
                    if (state->timeline.filter.enabled) {
                        list(SeriesSource_ScriptFiltered);
                    }

                    if (candidate_count == 0) {
                        ImGui::Text("No volume properties available.");
                    }
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Render")) {
                    // One or the other (or neither): the two are separate renderers and never mixed
                    if (ImGui::RadioButton("Off", !dvr.enabled && !iso.enabled)) {
                        dvr.enabled = false;
                        iso.enabled = false;
                    }
                    if (ImGui::RadioButton("Direct Volume Rendering", dvr.enabled)) {
                        dvr.enabled = true;
                        iso.enabled = false;
                    }
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
                    if (ImGui::RadioButton("Iso Surfaces", iso.enabled)) {
                        iso.enabled = true;
                        dvr.enabled = false;
                    }
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
                            for (int i = 0; i < (int)ColorMapping::Attribute; ++i) {
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
                gbuffer_init(&gbuf, width, height);
            }

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
            PickingHit hit = {};

            if (surface_state.hovered) {
                InteractionSurfaceHitArgs args = {
                    .picking_surface = &picking_surface,
                    .picking_handler = state->picking_handler,
                    .fbo = gbuf.fbo,
                    .width = gbuf.width,
                    .height = gbuf.height,
                    .clip_to_world = clip_to_world,
                };

                interaction_surface_hit_extract(&hit, surface_state, args);

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
            };

            InteractionSurfaceViewTransformResult view_result = interaction_surface_view_transform_apply(&target, surface_state, view_args);
            if (view_result.reset_requested) {
                ViewTransform reset_transform = default_view;
                if (hit.depth < 1.0f) {
                    reset_transform.distance = target.distance;
                    reset_transform.orientation = camera.orientation;
                    reset_transform.position = hit.world_pos + camera.orientation * vec3_set(0, 0, target.distance);
                }
                target = reset_transform;
            }

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
            gbuffer_clear(&gbuf);

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

            const bool selected_property = volume_key.path[0] != '\0';

            size_t num_reps = md_array_size(reps);
            if (selected_property && show_reference_structures && num_reps > 0) {
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

            bool iso_written = false;
            if (show_density_volume) {
                if (dvr.enabled) {
                    volume::DvrRenderDesc vol_desc = {
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
                        .tf = {
                            .min_value = dvr.tf.min_val,
                            .max_value = dvr.tf.max_val,
                        },
                    };
                    volume::render_dvr(vol_desc);
                } else if (iso.enabled) {
                    volume::IsoRenderDesc vol_desc = {
                        .render_target = {
                            .depth  = gbuf.tex.depth,
                            .color  = gbuf.tex.transparency_hdr,
                            .width  = gbuf.width,
                            .height = gbuf.height,
                            .clear_color = true,
                        },
                        .texture = {
                            .density_volume = volume_texture.id,
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
                            .count = iso.count,
                            .values = iso.values,
                            .colors = iso.colors,
                        },
                        // Lit like the compose pass lights the reference structures: env = background / 4
                        .shading = {
                            .env_radiance = state->visuals.background.color * state->visuals.background.intensity * 0.25f,
                            .roughness = 0.3f,
                            .dir_radiance = {10,10,10},
                            .ior = 1.5f,
                        },
                    };
                    iso_written = volume::render_isosurfaces(vol_desc);
                }
            }

            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
            glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
            glViewport(0, 0, gbuf.width, gbuf.height);
            glScissor(0, 0, gbuf.width, gbuf.height);

            if (show_bounding_box && selected_property) {
                glEnable(GL_DEPTH_TEST);
                glDepthMask(GL_TRUE);
                glEnable(GL_BLEND);
                // Into the transparency buffer, which holds premultiplied colour
                glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

				immediate::Queue* queue = immediate::queue_create("Density Volume Bounds");

                immediate::set_model(queue, model_mat);

                uint32_t box_color = convert_color(bounding_box_color);
                uint32_t clip_color = convert_color(clip_volume_color);
                immediate::box_wireframe(queue, {0,0,0}, {1,1,1}, box_color);
                immediate::box_wireframe(queue, clip_volume.min, clip_volume.max, clip_color);

                immediate::RenderParams params = {};
                params.view = world_to_view;
                params.proj = view_to_clip;
                immediate::render(queue, params);

				immediate::queue_destroy(queue);

                glDisable(GL_BLEND);
            }

            PUSH_GPU_SECTION("Postprocessing")
            postprocess_pipeline::Settings postprocess_settings = {};
            postprocess_pipeline::Inputs postprocess_inputs = {};

            postprocess_settings.background_color = state->visuals.background.color * state->visuals.background.intensity;
            postprocess_settings.tonemap.enabled = state->visuals.tonemapping.enabled;
            postprocess_settings.tonemap.mode = state->visuals.tonemapping.tonemapper;
            postprocess_settings.tonemap.exposure = state->visuals.tonemapping.exposure;
            postprocess_settings.tonemap.gamma = state->visuals.tonemapping.gamma;
            postprocess_settings.ssao.enabled = false;
            postprocess_settings.dof.enabled = false;
            postprocess_settings.fxaa.enabled = true;
            postprocess_settings.taa.enabled = false;
            postprocess_settings.sharpen.enabled = false;

            postprocess_inputs.depth = gbuf.tex.depth;
            postprocess_inputs.color = gbuf.tex.color;
            postprocess_inputs.normal = gbuf.tex.normal;
            postprocess_inputs.velocity = gbuf.tex.velocity;
            postprocess_inputs.transparency = gbuf.tex.transparency;
            postprocess_inputs.transparency_hdr = iso_written ? gbuf.tex.transparency_hdr : 0;

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

            postprocess_pipeline::execute(postprocess_inputs, postprocess_settings, view_param);
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
                workspace_register_window("DensityVolume", &show_window);
                break;
            }
            case viamd::EventType_ViamdSerialize:
                serialize(*(viamd::serialization_state_t*)e.payload);
                break;
            case viamd::EventType_ViamdDeserializeBegin:
                reset_workspace_settings();
                break;
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t& state = *(viamd::deserialization_state_t*)e.payload;
                if (str_eq(viamd::section_header(state), STR_LIT("DensityVolume"))) {
                    deserialize(state);
                }
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