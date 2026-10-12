#define IMGUI_DEFINE_MATH_OPERATORS

#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_vec_math.h>
#include <core/md_bitfield.h>
#include <core/md_arena_allocator.h>
#include <core/md_hash.h>
#include <core/md_os.h>

#include <md_csv.h>
#include <md_xvg.h>
#include <md_util.h>
#include <md_script.h>
#include <md_filter.h>

#include <viamd_event.h>
#include <viamd.h>
#include <event.h>

#include <serialization_utils.h>
#include <imgui_widgets.h>
#include <implot_widgets.h>
#include <implot_internal.h>
#include <components/analysis_density.h>
#include <components/analysis_density_renderer.h>
#include <cmath>
#include <memory>
#include <vector>

// The layers of the density of a series, in the order they are drawn
enum DensityLayer { DensityLayer_Map, DensityLayer_Levels, DensityLayer_Lines, DensityLayer_Count };

// What is known about a series between frames: its points in the triangle and the density of them
struct SeriesCache {
    uint32_t uid = 0;
    char path[SERIES_PATH_CAP] = "";
    uint64_t version = 0;
    uint32_t num_frames = 0;
    uint32_t count = 0;
    std::vector<vec2_t> coords;     // member major: frame f of member k at k * num_frames + f

    analysis_density::Triangle density;
    uint64_t density_hash = 0;
    analysis_density::Renderer renderer[DensityLayer_Count];
    uint64_t render_key[DensityLayer_Count] = {};
    bool rendered[DensityLayer_Count] = {};
    bool used = false;
};

struct Shapespace : viamd::EventHandler {
    // The shortcut: a selection which becomes a shapespace() property of the script
    char input[256] = "all";
    char error[256] = "";
    char property_name[128] = "shape";
    char property_error[256] = "";
    bool input_valid = false;
    uint64_t input_hash = 0;
    md_array(md_bitfield_t) bitfields = 0;
    md_bitfield_t joined_bitfield = {};

    // A property just created, added to the plot once the script has evaluated it
    char pending_name[128] = "";
    double pending_time = 0.0;

    bool show_window = false;

    // What the plot shows. Only the first subplot is used, the layout is the one of the timeline's
    PlotSubplot plot[PLOT_MAX_SUBPLOTS];

    std::vector<std::unique_ptr<SeriesCache>> caches;

    md_allocator_i* arena = 0;
    ApplicationState* app_state = 0;

    Shapespace() { viamd::event_system_register_handler(*this); }

    void clear_caches() {
        for (auto& c : caches) {
            for (auto& r : c->renderer) r.destroy();
        }
        caches.clear();
    }

    void clear_input() {
        input_valid = false;
        input_hash = 0;
        bitfields = 0;
        error[0] = '\0';
        if (arena) {
            md_arena_allocator_reset(arena);
            joined_bitfield = {};
            md_bitfield_init(&joined_bitfield, arena);
        }
    }

    void reset_settings() {
        str_copy_to_char_buf(input, sizeof(input), STR_LIT("all"));
        str_copy_to_char_buf(property_name, sizeof(property_name), STR_LIT("shape"));
        property_error[0] = '\0';
        pending_name[0] = '\0';
        plot_clear(plot, PLOT_MAX_SUBPLOTS);
        clear_caches();
        clear_input();
    }

    SeriesCache& cache_for(uint32_t uid, const char* path) {
        for (auto& c : caches) {
            if (c->uid != uid) continue;
            if (strcmp(c->path, path) != 0) {
                str_copy_to_char_buf(c->path, sizeof(c->path), str_from_cstr(path));
                c->version = 0;
                c->density_hash = 0;
            }
            return *c;
        }
        caches.push_back(std::make_unique<SeriesCache>());
        SeriesCache& c = *caches.back();
        c.uid = uid;
        str_copy_to_char_buf(c.path, sizeof(c.path), str_from_cstr(path));
        return c;
    }

    // The weights of every frame as points in the triangle
    static void update_points(SeriesCache& c, const SeriesVectorView& v) {
        if (c.version == v.version && c.num_frames == v.num_frames && c.count == v.count) return;
        c.version = v.version;
        c.num_frames = v.num_frames;
        c.count = v.count;
        c.density_hash = 0;
        c.coords.resize((size_t)v.count * v.num_frames);

        const vec2_t p[3] = {{0.0f, 0.0f}, {1.0f, 0.0f}, {0.5f, 0.86602540378f}};
        for (size_t k = 0; k < v.count; ++k) {
            for (size_t f = 0; f < v.num_frames; ++f) {
                const float* w = v.data + (f * v.count + k) * 3;
                // A frame which has not been evaluated holds zeros, not a shape: placed far outside the triangle, not as NaN
                const bool evaluated = w[0] != 0.0f || w[1] != 0.0f || w[2] != 0.0f;
                c.coords[k * v.num_frames + f] = evaluated ? p[0] * w[0] + p[1] * w[1] + p[2] * w[2] : vec2_t{-1.0e6f, -1.0e6f};
            }
        }
    }

    // The density of the members of a series over its frames, and the layers of it which the series shows
    void update_density(SeriesCache& c, const PlotSeries& s, size_t beg, size_t end) {
        struct Key {
            uint64_t version;
            uint64_t mask[4];
            uint64_t beg, end;
            float smoothing;
        } key = {};
        key.version = c.version;
        for (uint32_t k = 0; k < c.count && k < PLOT_MAX_POPULATION; ++k) {
            if (c.count == 1 || s.population_mask.test(k)) key.mask[k / 64] |= 1ULL << (k % 64);
        }
        key.beg = beg;
        key.end = end;
        key.smoothing = s.smoothing;

        const uint64_t hash = md_hash64(&key, sizeof(key), 0) | 1;
        if (hash != c.density_hash) {
            c.density_hash = hash;
            std::vector<vec2_t> points;
            size_t groups = 0;
            for (uint32_t k = 0; k < c.count; ++k) {
                if (!(c.count == 1 || (k < PLOT_MAX_POPULATION && s.population_mask.test(k)))) continue;
                points.insert(points.end(), c.coords.begin() + (size_t)k * c.num_frames, c.coords.begin() + (size_t)(k + 1) * c.num_frames);
                ++groups;
            }
            c.density.update(points.data(), c.num_frames, groups, beg, end, s.smoothing);
            for (bool& r : c.rendered) r = false;
        }

        const bool enabled[DensityLayer_Count] = { s.density_map, s.iso_levels, s.iso_lines };
        const analysis_density::Display display[DensityLayer_Count] = { analysis_density::Colormap, analysis_density::FilledLevels, analysis_density::Isolines };
        for (int layer = 0; layer < DensityLayer_Count; ++layer) {
            if (!enabled[layer]) continue;
            const struct { uint64_t density; uint32_t colormap; float color[4]; uint32_t pad; } rkey = { hash, (uint32_t)s.density_colormap, { s.color.x, s.color.y, s.color.z, s.color.w }, 0 };
            const uint64_t render_key = md_hash64(&rkey, sizeof(rkey), (uint64_t)layer) | 1;
            if (c.rendered[layer] && c.render_key[layer] == render_key) continue;
            c.rendered[layer] = c.renderer[layer].render(c.density, display[layer], s.density_colormap, s.color);
            c.render_key[layer] = render_key;
        }
    }

    // The selection of the shortcut, evaluated when it changes
    void update_input() {
        const uint64_t hash = md_hash64(input, sizeof(input), app_state->script.ir_fingerprint ^ app_state->mold.sys.topology_version) | 1;
        if (hash == input_hash) return;
        clear_input();
        input_hash = hash;
        if (!md_filter_evaluate(&bitfields, str_from_cstr(input), &app_state->mold.sys, &app_state->mold.state, app_state->script.ir, NULL, error, sizeof(error), arena)) {
            return;
        }
        if (!md_array_size(bitfields)) {
            snprintf(error, sizeof(error), "The selection contains no structures.");
            return;
        }
        for (size_t i = 0; i < md_array_size(bitfields); ++i) {
            if (!md_bitfield_popcount(&bitfields[i])) {
                snprintf(error, sizeof(error), "Selection group %zu contains no atoms.", i + 1);
                return;
            }
            md_bitfield_or_inplace(&joined_bitfield, &bitfields[i]);
        }
        input_valid = true;
    }

    void create_property() {
        property_error[0] = '\0';
        if (!md_script_identifier_name_valid(str_from_cstr(property_name))) {
            snprintf(property_error, sizeof(property_error), "Enter a valid script identifier as the property name.");
            return;
        }
        char statement[1024];
        const bool whole_system = str_eq(str_trim(str_from_cstr(input)), STR_LIT("all"));
        snprintf(statement, sizeof(statement), "%s = shapespace(%s);", property_name, whole_system ? "" : input);
        if (!script_append_statement(app_state, str_from_cstr(statement))) {
            snprintf(property_error, sizeof(property_error), "Could not create the property: check the name and selection for script errors or a duplicate identifier.");
            return;
        }
        snprintf(pending_name, sizeof(pending_name), "%s", property_name);
        pending_time = ImGui::GetTime();
    }

    // Adds a property created from the shortcut once the script has evaluated it
    void update_pending() {
        if (!pending_name[0]) return;
        char path[SERIES_PATH_CAP];
        snprintf(path, sizeof(path), "script/%s", pending_name);
        const SeriesKey key = series_key(SeriesSource_Script, str_from_cstr(path));
        SeriesVectorView v;
        if (series_resolve_vector(&v, app_state, key, 3)) {
            plot_add_vector_series(app_state, plot[0], key, false);
            pending_name[0] = '\0';
        } else if (ImGui::GetTime() - pending_time > 20.0) {
            pending_name[0] = '\0';
        }
    }

    void draw_window() {
        if (!show_window) return;

        update_input();
        update_pending();

        ImGui::SetNextWindowSize({300,350}, ImGuiCond_FirstUseEver);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(2, 2));
        defer { ImGui::PopStyleVar(1); };

        if (ImGui::Begin("Shapespace", &show_window, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
            PlotSubplot& sp = plot[0];

            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("Properties")) {
                    series_vector_properties_menu(app_state, sp, 3);
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Layers")) {
                    series_vector_layers_menu(app_state, sp, 3);
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Export", sp.count > 0)) {
                    if (ImGui::MenuItem("XVG")) export_shape_space(STR_LIT("xvg"));
                    if (ImGui::MenuItem("CSV")) export_shape_space(STR_LIT("csv"));
                    ImGui::EndMenu();
                }
                ImGui::EndMenuBar();
            }

            // The shortcut
            ImGui::SetNextItemWidth(MAX(ImGui::GetContentRegionAvail().x, 150.f));
            ImGui::InputQuery("##input", input, sizeof(input), input_valid);
            if (ImGui::IsItemHovered()) {
                if (input_valid) {
                    md_bitfield_copy(&app_state->selection.highlight_mask, &joined_bitfield);
                } else if (error[0] != '\0') {
                    ImGui::SetTooltip("%s", error);
                }
            }
            ImGui::SetNextItemWidth(MAX(60.0f, ImGui::GetContentRegionAvail().x - 145.0f));
            ImGui::InputText("##Property Name", property_name, sizeof(property_name));
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Name of the property.");
            ImGui::SameLine();
            ImGui::BeginDisabled(!input_valid);
            if (ImGui::Button("Create Property")) create_property();
            ImGui::EndDisabled();
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
                ImGui::SetTooltip("Add '%s = shapespace(...)' for this selection to the script and show it here.", property_name);
            }
            if (property_error[0]) ImGui::TextWrapped("%s", property_error);
            if (!input_valid && error[0]) ImGui::TextWrapped("%s", error);
            if (pending_name[0]) ImGui::TextDisabled("Waiting for '%s' to be evaluated...", pending_name);

            // The series, which this may add to or remove from
            series_vector_bar(app_state, sp, 3);

            // What the series are, before anything is drawn
            SeriesVectorView views[PLOT_MAX_SERIES_PER_SUBPLOT];
            bool resolved[PLOT_MAX_SERIES_PER_SUBPLOT] = {};
            char plot_id[PLOT_MAX_SERIES_PER_SUBPLOT][SERIES_PATH_CAP + 72];
            SeriesCache* series_cache[PLOT_MAX_SERIES_PER_SUBPLOT] = {};

            for (auto& c : caches) c->used = false;
            for (int j = 0; j < sp.count; ++j) {
                const PlotSeries& s = sp.series[j];
                resolved[j] = series_resolve_vector(&views[j], app_state, s.key, 3);
                series_label(views[j].label, sizeof(views[j].label), app_state, s.key);
                snprintf(plot_id[j], sizeof(plot_id[j]), "%s##%u", views[j].label, s.uid);
                if (!resolved[j]) continue;

                SeriesCache& c = cache_for(s.uid, s.key.path);
                c.used = true;
                series_cache[j] = &c;
                update_points(c, views[j]);
                if (s.iso_lines || s.iso_levels || s.density_map) {
                    update_density(c, s, 0, c.num_frames);
                } else {
                    c.density_hash = 0;
                    for (bool& r : c.rendered) r = false;
                }
            }
            for (size_t i = caches.size(); i-- > 0;) {
                if (!caches[i]->used) {
                    for (auto& r : caches[i]->renderer) r.destroy();
                    caches.erase(caches.begin() + i);
                }
            }

            ImPlot::PushStyleVar(ImPlotStyleVar_PlotPadding, ImVec2(2,2));
            defer { ImPlot::PopStyleVar(1); };

            const ImPlotFlags flags = ImPlotFlags_Equal | ImPlotFlags_NoMenus | ImPlotFlags_NoMouseText | ImPlotFlags_NoLegend;
            const ImPlotAxisFlags axis_flags = ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoTickLabels | ImPlotAxisFlags_NoTickMarks;

            const float x_reset[2] = {-0.10f, 1.10f};
            const float y_reset[2] = {-0.10f, 0.98f};
            const double zoom_constraints[2] = {0.1, 1.5};

            if (ImPlot::BeginPlot("##Shapespace Plot", ImVec2(-1,-1), flags)) {
                ImPlot::SetupAxesLimits(x_reset[0], x_reset[1], y_reset[0], y_reset[1], ImGuiCond_Appearing);
                ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, -1.0f, 2.0f);
                ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, -1.0f, 2.0f);
                ImPlot::SetupAxisZoomConstraints(ImAxis_X1, zoom_constraints[0], zoom_constraints[1]);
                ImPlot::SetupAxisZoomConstraints(ImAxis_Y1, zoom_constraints[0], zoom_constraints[1]);
                ImPlot::SetupAxes(0, 0, axis_flags, axis_flags);
                ImPlot::SetupFinish();

                const ImPlotPoint lin(0.0, 0.0);
                const ImPlotPoint pla(1.0, 0.0);
                const ImPlotPoint iso(0.5, 0.86602540378);

                const ImVec2 p0 = ImPlot::PlotToPixels(lin);
                const ImVec2 p1 = ImPlot::PlotToPixels(pla);
                const ImVec2 p2 = ImPlot::PlotToPixels(iso);

                static const char* text_sph = "Isotropic";
                static const char* text_pla = "Planar";
                static const char* text_lin = "Linear";

                const ImVec2 text_sph_offset = p2 + ImGui::CalcTextSize(text_sph) * ImVec2(-0.5f, -1.2f);
                const ImVec2 text_pla_offset = p1 + ImGui::CalcTextSize(text_pla) * ImVec2(-0.5f,  0.2f);
                const ImVec2 text_lin_offset = p0 + ImGui::CalcTextSize(text_lin) * ImVec2(-0.5f,  0.2f);

                ImPlot::PushPlotClipRect();
                ImPlot::GetPlotDrawList()->AddTriangleFilled(p0, p1, p2, IM_COL32(255,255,255,20));
                ImPlot::GetPlotDrawList()->AddTriangle(p0, p1, p2, IM_COL32(255,255,255,50));
                ImPlot::PopPlotClipRect();

                // The layers of the density of each series, below the points
                for (int j = 0; j < sp.count; ++j) {
                    const SeriesCache* c = series_cache[j];
                    if (!c) continue;
                    for (int layer = 0; layer < DensityLayer_Count; ++layer) {
                        if (c->rendered[layer]) c->renderer[layer].draw(sp.series[j].density_alpha);
                    }
                }
                ImPlot::PushPlotClipRect();

                ImPlot::GetPlotDrawList()->AddText(text_sph_offset, IM_COL32(255, 255, 255, 255), text_sph);
                ImPlot::GetPlotDrawList()->AddText(text_pla_offset, IM_COL32(255, 255, 255, 255), text_pla);
                ImPlot::GetPlotDrawList()->AddText(text_lin_offset, IM_COL32(255, 255, 255, 255), text_lin);
                ImPlot::PopPlotClipRect();

                ImPlot::PushStyleVar(ImPlotStyleVar_Marker, ImPlotMarker_Square);
                ImPlot::PushStyleColor(ImPlotCol_MarkerFill, ImVec4(0,0,0,0));
                ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, ImVec4(0,0,0,0));
                ImPlot::PlotScatter("##Hidden reset helper", x_reset, y_reset, 2);
                ImPlot::PopStyleColor(2);
                ImPlot::PopStyleVar();

                auto getter = [](int idx, void* user_data) -> ImPlotPoint {
                    const vec2_t* coords = (vec2_t*)user_data;
                    return {coords[idx].x, coords[idx].y};
                };

                int   hovered_series = -1;
                int   hovered_member = -1;
                float hovered_point_min_dist = FLT_MAX;
                uint32_t hovered_frame = 0;

                const vec2_t mouse_coord = {(float)ImPlot::GetPlotMousePos().x, (float)ImPlot::GetPlotMousePos().y};

                const ImPlotRect lim = ImPlot::GetPlotLimits();
                const float scl = (float)lim.X.Max - (float)lim.X.Min;
                const float MAX_D2 = 0.0004f * scl * scl;

                // The members of every series at the displayed frame
                for (int j = 0; j < sp.count; ++j) {
                    const PlotSeries& s = sp.series[j];
                    const char* id = plot_id[j];
                    const SeriesCache* c = series_cache[j];
                    if (!c || !s.show_points || c->num_frames == 0) {
                        ImPlot::PlotDummy(id);
                        continue;
                    }

                    const uint32_t frame = MIN((uint32_t)MAX(0.0, app_state->animation.frame + 0.5), c->num_frames - 1);
                    for (uint32_t k = 0; k < c->count; ++k) {
                        if (c->count > 1 && !s.population_mask.test(k)) continue;

                        ImVec4 col = s.color;
                        if (s.use_colormap && c->count > 1) {
                            col = ImPlot::SampleColormap((float)k / (float)(c->count - 1), s.colormap);
                            col.w *= s.colormap_alpha;
                        }
                        const vec2_t* point = c->coords.data() + (size_t)k * c->num_frames + frame;
                        ImPlot::SetNextMarkerStyle(s.marker, s.marker_size, col, 1.0f, col);
                        ImPlot::PlotScatterG(id, getter, (void*)point, 1);

                        const vec2_t delta = mouse_coord - *point;
                        const float d2 = vec2_dot(delta, delta);
                        if (d2 < MAX_D2 && d2 < hovered_point_min_dist) {
                            hovered_point_min_dist = d2;
                            hovered_series = j;
                            hovered_member = (int)k;
                            hovered_frame = frame;
                        }
                    }
                }

                if (ImPlot::IsPlotHovered()) {
                    md_bitfield_clear(&app_state->selection.highlight_mask);
                }

                if (hovered_series != -1 && resolved[hovered_series]) {
                    const SeriesVectorView& v = views[hovered_series];
                    const float* w = v.data + ((size_t)hovered_frame * v.count + hovered_member) * 3;
                    char buf[160] = "";
                    int len = 0;
                    if (v.count > 1) len += snprintf(buf, sizeof(buf), "Structure: %i, ", hovered_member + 1);
                    snprintf(buf + len, sizeof(buf) - len, "Frame: %u, lin: %.2f, plan: %.2f, iso: %.2f", hovered_frame, w[0], w[1], w[2]);
                    ImGui::SetTooltip("%s", buf);
                    script_visualize_payload(app_state, v.vis_payload, v.count > 1 ? hovered_member : -1, ~MD_SCRIPT_VISUALIZE_SDF);
                    script_set_hovered_property(app_state, str_from_cstr(v.script_ident), hovered_member);
                }

                ImPlot::EndPlot();
            }

            series_vector_drag_out(sp);
        }
        ImGui::End();
    }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& e = events[i];

            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                app_state = (ApplicationState*)e.payload;
                arena = md_arena_allocator_create(app_state->allocator.persistent, MEGABYTES(1));
                md_bitfield_init(&joined_bitfield, arena);
                workspace_register_window("ShapeSpace", &show_window);
                break;
            }
            case viamd::EventType_ViamdDeserializeBegin:
                reset_settings();
                break;
            case viamd::EventType_ViamdSystemFree:
            case viamd::EventType_ViamdTrajectoryFree:
            case viamd::EventType_ViamdSystemInit:
            case viamd::EventType_ViamdTrajectoryInit:
            case viamd::EventType_ViamdLoadData:
                clear_caches();
                clear_input();
                break;
            case viamd::EventType_ViamdShutdown:
                clear_caches();
                md_arena_allocator_destroy(arena);
                arena = 0;
                break;
            case viamd::EventType_ViamdFrameTick:
                draw_window();
                break;
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("Shapespace", &show_window);
                break;
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t& state = *(viamd::deserialization_state_t*)e.payload;
                if (str_eq(viamd::section_header(state), STR_LIT("ShapeSpace"))) {
                    str_t ident, arg;
                    while (viamd::next_entry(ident, arg, state)) {
                        if (str_eq(ident, STR_LIT("Filter"))) {
                            str_t filter;
                            viamd::extract_str(filter, arg);
                            str_copy_to_char_buf(input, sizeof(input), filter);
                        } else if (str_eq(ident, STR_LIT("PropertyName"))) {
                            viamd::extract_to_char_buf(property_name, sizeof(property_name), arg);
                        }
                    }
                } else {
                    int num_subplots = 1;
                    plot_layout_deserialize(state, STR_LIT("ShapeSpaceLayout"), STR_LIT("ShapeSpaceSeries"), app_state, plot, &num_subplots);
                }
                break;
            }
            case viamd::EventType_ViamdSerialize: {
                viamd::serialization_state_t& state = *(viamd::serialization_state_t*)e.payload;
                viamd::write_section_header(state, STR_LIT("ShapeSpace"));
                viamd::write_str(state, STR_LIT("Filter"), str_from_cstr(input));
                viamd::write_str(state, STR_LIT("PropertyName"), str_from_cstr(property_name));
                plot_layout_serialize(state, STR_LIT("ShapeSpaceLayout"), STR_LIT("ShapeSpaceSeries"), plot, 1);
                break;
            }
            default:
                break;
            }
        }
    }

    // The weights of every series of the plot, a column each
    void export_shape_space(str_t ext) {
        char path_buf[2048];
        if (!application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, ext)) return;
        str_t path = {path_buf, strnlen(path_buf, sizeof(path_buf))};
        md_array(str_t)  column_labels = 0;
        md_array(float*) column_values = 0;

        md_temp_scope_t temp_scope = md_temp_begin();
        md_allocator_i* temp_arena = md_temp_allocator(temp_scope);
        defer { md_temp_end(temp_scope); };

        // The time column below is timeline.x_values, which is held in display units, so the
        // label has to come from the same place rather than from the trajectory.
        md_unit_t time_unit = app_state->timeline.time_unit;

        str_t x_label = STR_INIT("Frame");
        if (!md_unit_is_none(time_unit)) {
            char unit_buf[64];
            md_unit_print(unit_buf, sizeof(unit_buf), time_unit);
            x_label = str_printf(temp_arena, "Time (%s)", unit_buf);
        }
        md_array_push(column_labels, x_label, temp_arena);
        md_array_push(column_values, NULL, temp_arena);

        const PlotSubplot& sp = plot[0];
        size_t num_rows = 0;
        for (int j = 0; j < sp.count; ++j) {
            SeriesVectorView v;
            if (!series_resolve_vector(&v, app_state, sp.series[j].key, 3)) continue;
            num_rows = MAX(num_rows, (size_t)v.num_frames);
            static const char* const names[3] = {"linear", "planar", "isotropic"};
            for (uint32_t k = 0; k < v.count; ++k) {
                for (int c = 0; c < 3; ++c) {
                    const str_t label = v.count > 1 ? str_printf(temp_arena, "%s [%u] %s", v.label, k + 1, names[c])
                                                    : str_printf(temp_arena, "%s %s", v.label, names[c]);
                    md_array_push(column_labels, label, temp_arena);
                    md_array_push(column_values, NULL, temp_arena);
                    float** col = &column_values[md_array_size(column_values) - 1];
                    for (uint32_t f = 0; f < v.num_frames; ++f) {
                        md_array_push(*col, v.data[((size_t)f * v.count + k) * 3 + c], temp_arena);
                    }
                }
            }
        }
        if (!num_rows) return;

        const float* time = app_state->timeline.x_values;
        for (size_t i = 0; i < num_rows; ++i) {
            md_array_push(column_values[0], time ? time[i] : (float)i, temp_arena);
        }

        const size_t num_cols = md_array_size(column_values);
        ASSERT(num_cols == md_array_size(column_labels));
        if (str_eq(ext, STR_LIT("csv"))) {
            md_csv_write_to_file(column_values, column_labels, num_cols, num_rows, path);
        } else if (str_eq(ext, STR_LIT("xvg"))) {
            str_t title = STR_INIT("Shape Space");
            str_t y_label = STR_INIT("Shape Weight Ratio");
            str_t header = md_xvg_format_header(title, x_label, y_label, md_array_size(column_labels) - 1, column_labels + 1, temp_arena);
            str_t xvg    = md_xvg_format(header, num_cols, num_rows, column_values, temp_arena);
            md_file_t file = {0};
            if (md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                md_file_write(file, xvg.ptr, xvg.len);
                md_file_close(&file);
            }
        } else {
            MD_LOG_DEBUG("Unrecognized export format");
            ASSERT(false);
        }
    }
};

static Shapespace instance = {};
