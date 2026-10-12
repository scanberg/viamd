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

// What is known about a series between frames: its points in the triangle and the density of them
struct SeriesCache {
    char path[SERIES_PATH_CAP] = "";
    uint64_t version = 0;
    uint32_t num_frames = 0;
    uint32_t count = 0;
    std::vector<vec2_t> coords;     // member major: frame f of member k at k * num_frames + f

    analysis_density::Triangle density;
    analysis_density::Renderer renderer;
    uint64_t density_hash = 0;
    bool render_dirty = true;
    bool rendered = false;
    analysis_density::Display rendered_display = analysis_density::Scatter;
    ImPlotColormap rendered_colormap = -1;
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

    analysis_density::Display display = analysis_density::Scatter;
    bool show_scatter = true;
    float smoothing = 2.0f;
    float density_opacity = 0.8f;
    ImPlotColormap density_colormap = ImPlotColormap_Viridis;

    std::vector<std::unique_ptr<SeriesCache>> caches;

    md_allocator_i* arena = 0;
    ApplicationState* app_state = 0;

    Shapespace() { viamd::event_system_register_handler(*this); }

    void clear_caches() {
        for (auto& c : caches) c->renderer.destroy();
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
        display = analysis_density::Scatter;
        show_scatter = true;
        smoothing = 2.0f;
        density_opacity = 0.8f;
        density_colormap = ImPlotColormap_Viridis;
        plot_clear(plot, PLOT_MAX_SUBPLOTS);
        clear_caches();
        clear_input();
    }

    SeriesCache& cache_for(const char* path) {
        for (auto& c : caches) {
            if (strcmp(c->path, path) == 0) return *c;
        }
        caches.push_back(std::make_unique<SeriesCache>());
        SeriesCache& c = *caches.back();
        str_copy_to_char_buf(c.path, sizeof(c.path), str_from_cstr(path));
        return c;
    }

    void frame_range(size_t num_frames, size_t& beg, size_t& end) const {
        beg = 0;
        end = num_frames;
        if (app_state->timeline.filter.enabled) {
            beg = (size_t)CLAMP(app_state->timeline.filter.beg_frame, 0.0, (double)num_frames);
            end = (size_t)CLAMP(app_state->timeline.filter.end_frame, (double)beg, (double)num_frames);
        }
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
        key.smoothing = smoothing;

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
            c.density.update(points.data(), c.num_frames, groups, beg, end, smoothing);
            c.render_dirty = true;
        }
        if (c.render_dirty || !c.rendered || c.rendered_display != display || c.rendered_colormap != density_colormap) {
            c.rendered = c.renderer.render(c.density, display, density_colormap);
            if (c.rendered) {
                c.render_dirty = false;
                c.rendered_display = display;
                c.rendered_colormap = density_colormap;
            }
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
            plot_add_series(app_state, plot[0], key);
            pending_name[0] = '\0';
        } else if (ImGui::GetTime() - pending_time > 20.0) {
            pending_name[0] = '\0';
        }
    }

    void draw_properties_menu() {
        PlotSubplot& sp = plot[0];
        int num_listed = 0;
        series_for_each_script_property(app_state, SeriesSource_Script, MD_SCRIPT_PROPERTY_FLAG_VECTOR, [&](const SeriesKey& key) {
            SeriesVectorView v;
            if (!series_resolve_vector(&v, app_state, key, 3)) return;
            plot_series_list_item(app_state, VECTOR_SERIES_DND, sp, key);
            num_listed += 1;
        });
        if (num_listed == 0) {
            ImGui::TextDisabled("No shape properties, create one below or define shapespace() in the script editor");
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
                    draw_properties_menu();
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Export", sp.count > 0)) {
                    if (ImGui::MenuItem("XVG")) export_shape_space(STR_LIT("xvg"));
                    if (ImGui::MenuItem("CSV")) export_shape_space(STR_LIT("csv"));
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Settings")) {
                    if (ImGui::BeginCombo("Display", analysis_density::display_names[display])) {
                        for (int i = 0; i < analysis_density::DisplayCount; ++i) {
                            if (ImGui::Selectable(analysis_density::display_names[i], display == i)) {
                                display = (analysis_density::Display)i;
                            }
                        }
                        ImGui::EndCombo();
                    }
                    if (display != analysis_density::Scatter) {
                        ImGui::Checkbox("Scatter overlay", &show_scatter);
                        ImGui::SliderFloat("Smoothing", &smoothing, 0.25f, 6.0f, "%.2f");
                        ImPlot::ColormapSelection("Density Color Map", &density_colormap);
                        ImGui::SliderFloat("Density Opacity", &density_opacity, 0.0f, 1.0f);
                        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Probability density relative to a uniform triangle. Independent of sample count.");
                    }
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

            // What the series are, before anything is drawn
            SeriesVectorView views[PLOT_MAX_SERIES_PER_SUBPLOT];
            bool resolved[PLOT_MAX_SERIES_PER_SUBPLOT] = {};
            char plot_id[PLOT_MAX_SERIES_PER_SUBPLOT][SERIES_PATH_CAP + 72];
            SeriesCache* series_cache[PLOT_MAX_SERIES_PER_SUBPLOT] = {};
            size_t beg[PLOT_MAX_SERIES_PER_SUBPLOT] = {};
            size_t end[PLOT_MAX_SERIES_PER_SUBPLOT] = {};

            for (auto& c : caches) c->used = false;
            for (int j = 0; j < sp.count; ++j) {
                const PlotSeries& s = sp.series[j];
                resolved[j] = series_resolve_vector(&views[j], app_state, s.key, 3);
                series_label(views[j].label, sizeof(views[j].label), app_state, s.key);
                snprintf(plot_id[j], sizeof(plot_id[j]), "%s##%s", views[j].label, s.key.path);
                if (!resolved[j]) continue;

                SeriesCache& c = cache_for(s.key.path);
                c.used = true;
                series_cache[j] = &c;
                update_points(c, views[j]);
                frame_range(c.num_frames, beg[j], end[j]);
                if (display != analysis_density::Scatter) {
                    update_density(c, s, beg[j], end[j]);
                } else {
                    c.density_hash = 0;
                }
            }
            for (size_t i = caches.size(); i-- > 0;) {
                if (!caches[i]->used) {
                    caches[i]->renderer.destroy();
                    caches.erase(caches.begin() + i);
                }
            }

            ImPlot::PushStyleVar(ImPlotStyleVar_PlotPadding, ImVec2(2,2));
            defer { ImPlot::PopStyleVar(1); };

            const ImPlotFlags flags = ImPlotFlags_Equal | ImPlotFlags_NoMenus | ImPlotFlags_NoMouseText;
            const ImPlotAxisFlags axis_flags = ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoTickLabels | ImPlotAxisFlags_NoTickMarks;

            const float x_reset[2] = {-0.10f, 1.10f};
            const float y_reset[2] = {-0.10f, 0.98f};
            const double zoom_constraints[2] = {0.1, 1.5};

            int remove_idx = -1;

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
                if (display != analysis_density::Scatter) {
                    for (int j = 0; j < sp.count; ++j) {
                        if (series_cache[j] && series_cache[j]->rendered && !series_cache[j]->render_dirty) {
                            series_cache[j]->renderer.draw(density_opacity);
                        }
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
                int   hovered_frame  = -1;
                float hovered_point_min_dist = FLT_MAX;

                const vec2_t mouse_coord = {(float)ImPlot::GetPlotMousePos().x, (float)ImPlot::GetPlotMousePos().y};

                const ImPlotRect lim = ImPlot::GetPlotLimits();
                const float scl = (float)lim.X.Max - (float)lim.X.Min;
                const float MAX_D2 = 0.0001f * scl * scl;

                for (int j = 0; j < sp.count; ++j) {
                    PlotSeries& s = sp.series[j];
                    const char* id = plot_id[j];

                    if (!resolved[j]) {
                        // Nothing behind the series right now: still a legend entry, which can be removed
                        ImPlot::PlotDummy(id);
                    } else {
                        const SeriesCache& c = *series_cache[j];
                        const bool draw_scatter = display == analysis_density::Scatter || show_scatter;
                        const size_t count = beg[j] < end[j] ? end[j] - beg[j] : 0;
                        for (uint32_t k = 0; k < c.count; ++k) {
                            if (c.count > 1 && !s.population_mask.test(k)) continue;

                            ImVec4 col = s.color;
                            if (s.use_colormap && c.count > 1) {
                                col = ImPlot::SampleColormap((float)k / (float)(c.count - 1), s.colormap);
                                col.w *= s.colormap_alpha;
                            }
                            const vec2_t* points = c.coords.data() + (size_t)k * c.num_frames + beg[j];
                            if (draw_scatter) {
                                ImPlot::SetNextMarkerStyle(s.marker, s.marker_size, col, 1.0f, col);
                            } else {
                                ImPlot::SetNextMarkerStyle(s.marker, 0.0f, ImVec4(0,0,0,0), 0.0f, ImVec4(0,0,0,0));
                            }
                            ImPlot::PlotScatterG(id, getter, (void*)points, (int)count);

                            ImPlotItem* item = ImPlot::GetItem(id);
                            if (!item || !item->Show) continue;
                            if (item->LegendHovered && hovered_series == -1) {
                                hovered_series = j;
                            }
                            for (size_t f = beg[j]; f < end[j]; ++f) {
                                const vec2_t delta = mouse_coord - c.coords[(size_t)k * c.num_frames + f];
                                const float d2 = vec2_dot(delta, delta);
                                if (d2 < MAX_D2 && d2 < hovered_point_min_dist) {
                                    hovered_point_min_dist = d2;
                                    hovered_series = j;
                                    hovered_member = (int)k;
                                    hovered_frame  = (int)f;
                                }
                            }
                        }
                    }

                    if (ImPlot::BeginLegendPopup(id)) {
                        if (ImGui::DeleteButton("Remove")) {
                            remove_idx = j;
                            ImGui::CloseCurrentPopup();
                        }
                        if (ImGui::BeginCombo("Marker", ImPlot::GetMarkerName(s.marker))) {
                            for (int m = 0; m < ImPlotMarker_COUNT; ++m) {
                                if (ImGui::Selectable(ImPlot::GetMarkerName(m), s.marker == m)) s.marker = (ImPlotMarker)m;
                            }
                            ImGui::EndCombo();
                        }
                        ImGui::SliderFloat("Marker Size", &s.marker_size, 0.1f, 10.0f, "%.2f");
                        int hovered_pop_idx = -1;
                        plot_series_style_popup(app_state, s, resolved[j] ? (int)views[j].count : 1, views[j].script_ident, views[j].vis_payload, &hovered_pop_idx);
                        ImPlot::EndLegendPopup();
                    }
                }

                if (ImPlot::IsPlotHovered()) {
                    md_bitfield_clear(&app_state->selection.highlight_mask);
                }

                // Redraw hovered index
                if (hovered_series != -1 && resolved[hovered_series]) {
                    const SeriesVectorView& v = views[hovered_series];
                    const SeriesCache& c = *series_cache[hovered_series];
                    if (hovered_frame != -1) {
                        ImPlot::PushStyleVar(ImPlotStyleVar_Marker, ImPlotMarker_Square);
                        ImPlot::PushStyleVar(ImPlotStyleVar_MarkerSize, sp.series[hovered_series].marker_size * 1.1f + 1.0f);
                        ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, ImVec4(1,1,1,1));
                        const vec2_t* point = c.coords.data() + (size_t)hovered_member * c.num_frames + hovered_frame;
                        ImPlot::PlotScatterG("##hovered idx", getter, (void*)point, 1);
                        ImPlot::PopStyleColor();
                        ImPlot::PopStyleVar(2);

                        const float* w = v.data + ((size_t)hovered_frame * v.count + hovered_member) * 3;
                        char buf[160] = "";
                        int len = 0;
                        if (v.count > 1) len += snprintf(buf, sizeof(buf), "Structure: %i, ", hovered_member + 1);
                        snprintf(buf + len, sizeof(buf) - len, "Frame: %i, lin: %.2f, plan: %.2f, iso: %.2f", hovered_frame, w[0], w[1], w[2]);
                        ImGui::SetTooltip("%s", buf);

                        if (ImGui::IsWindowFocused() && ImPlot::IsPlotHovered() && ImGui::GetIO().MouseClicked[0]) {
                            app_state->animation.frame = hovered_frame;
                        }
                    }
                    script_visualize_payload(app_state, v.vis_payload, v.count > 1 ? hovered_member : -1, ~MD_SCRIPT_VISUALIZE_SDF);
                    script_set_hovered_property(app_state, str_from_cstr(v.script_ident), hovered_member);
                }

                if (ImPlot::BeginDragDropTargetPlot()) {
                    if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(VECTOR_SERIES_DND)) {
                        if (payload->DataSize == sizeof(SeriesDragPayload)) {
                            const SeriesDragPayload* dnd = (const SeriesDragPayload*)payload->Data;
                            SeriesVectorView v;
                            if (series_resolve_vector(&v, app_state, dnd->key, 3)) {
                                plot_add_series(app_state, sp, dnd->key);
                            }
                        }
                    }
                    ImPlot::EndDragDropTarget();
                }

                ImPlot::EndPlot();
            }

            if (remove_idx != -1) {
                plot_remove_series(sp, remove_idx);
            }
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
                        } else if (str_eq(ident, STR_LIT("Display"))) {
                            viamd::extract_enum(display, arg, analysis_density::DisplayCount);
                        } else if (str_eq(ident, STR_LIT("ShowScatter"))) {
                            viamd::extract_bool(show_scatter, arg);
                        } else if (str_eq(ident, STR_LIT("Smoothing"))) {
                            viamd::extract_flt(smoothing, arg);
                            smoothing = std::isfinite(smoothing) ? CLAMP(smoothing, 0.25f, 6.0f) : 2.0f;
                        } else if (str_eq(ident, STR_LIT("DensityOpacity"))) {
                            viamd::extract_flt(density_opacity, arg);
                            density_opacity = std::isfinite(density_opacity) ? CLAMP(density_opacity, 0.0f, 1.0f) : 0.8f;
                        } else if (str_eq(ident, STR_LIT("DensityColormap"))) {
                            viamd::extract_enum(density_colormap, arg, ImPlot::GetColormapCount());
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
                viamd::write_int(state, STR_LIT("Display"), display);
                viamd::write_bool(state, STR_LIT("ShowScatter"), show_scatter);
                viamd::write_flt(state, STR_LIT("Smoothing"), smoothing);
                viamd::write_flt(state, STR_LIT("DensityOpacity"), density_opacity);
                viamd::write_int(state, STR_LIT("DensityColormap"), density_colormap);
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
