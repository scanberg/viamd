#define IMGUI_DEFINE_MATH_OPERATORS

#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_vec_math.h>
#include <core/md_bitfield.h>
#include <core/md_arena_allocator.h>
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
#include <atomic>
#include <cmath>

enum ColorMode {
    ColorMode_Uniform,
    ColorMode_Colormap,
    ColorMode_Count,
};

struct Shapespace : viamd::EventHandler {
    char input[256] = "all";
    char error[256] = "";

    bool input_valid = false;
    bool show_window = false;

    // Used to compare against the current 'app_state' if things should be reevaluated
    uint64_t eval_hash   = 0;

    // coords and weights should be of size num_frames * num_structures
    size_t num_frames = 0;
    size_t num_structures = 0;

    md_array(vec3_t) weights = 0;
    md_array(vec2_t) coords  = 0;
    md_array(vec3_t) pending_weights = 0;
    md_array(vec2_t) pending_coords = 0;
    md_array(md_bitfield_t) bitfields = 0;
    md_bitfield_t joined_bitfield = {};

    md_allocator_i* arena = 0;

    // Settings
    float marker_size = 1.5f;
    bool  use_mass = true;
    ColorMode mode = ColorMode_Colormap;
    ImPlotColormap colormap = 0;
    ImVec4 color = {0.5f,1.0f,0.5f,1.0f};
    analysis_density::Display display = analysis_density::Scatter;
    bool show_scatter = true;
    float smoothing = 2.0f;
    float density_opacity = 0.8f;
    ImPlotColormap density_colormap = ImPlotColormap_Viridis;
    // A script property (shapespace()) shown instead of the selection
    bool use_source = false;
    SeriesKey source_key = {};
    uint64_t source_version = 0;

    char property_name[128] = "shape";
    char property_error[256] = "";

    analysis_density::Triangle density;
    analysis_density::Renderer density_renderer;
    bool density_render_dirty = true;
    analysis_density::Display rendered_display = analysis_density::Scatter;
    ImPlotColormap rendered_colormap = -1;
    uint64_t data_generation = 1;
    size_t density_beg = 0;
    size_t density_end = 0;
    float density_smoothing = -1.0f;
    bool evaluation_pending = false;
    std::atomic<bool> evaluation_failed{false};

    task_system::ID evaluate_task = 0;

    ApplicationState* app_state = 0;

    Shapespace() { viamd::event_system_register_handler(*this); }

    void clear_data() {
        task_system::task_interrupt_and_wait_for(evaluate_task);
        evaluate_task = 0;
        evaluation_pending = false;
        bitfields = 0;
        weights = 0;
        coords = 0;
        pending_weights = 0;
        pending_coords = 0;
        num_frames = num_structures = 0;
        eval_hash = 0;
        input_valid = false;
        density.clear();
        density_render_dirty = true;
        ++data_generation;
        if (arena) {
            md_arena_allocator_reset(arena);
            joined_bitfield = {};
            md_bitfield_init(&joined_bitfield, arena);
        }
    }

    uint64_t evaluation_hash() const {
        uint64_t hash = md_hash64(input, sizeof(input), app_state->script.ir_fingerprint ^
            (1ULL << (uint64_t)use_mass) ^ data_generation);
        if (use_source) {
            hash = md_hash64(&source_version, sizeof(source_version), md_hash64_str(str_from_cstr(source_key.path), hash) ^ 1ULL);
        }
        return hash;
    }

    void frame_range(size_t& beg, size_t& end) const {
        beg = 0;
        end = num_frames;
        if (app_state->timeline.filter.enabled) {
            beg = (size_t)CLAMP(app_state->timeline.filter.beg_frame, 0.0, (double)num_frames);
            end = (size_t)CLAMP(app_state->timeline.filter.end_frame, (double)beg, (double)num_frames);
        }
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
        }
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
                clear_data();
                str_copy_to_char_buf(input, sizeof(input), STR_LIT("all"));
                marker_size = 1.5f;
                use_mass = true;
                mode = ColorMode_Colormap;
                colormap = 0;
                color = {0.5f,1.0f,0.5f,1.0f};
                display = analysis_density::Scatter;
                show_scatter = true;
                smoothing = 2.0f;
                density_opacity = 0.8f;
                density_colormap = ImPlotColormap_Viridis;
                str_copy_to_char_buf(property_name, sizeof(property_name), STR_LIT("shape"));
                property_error[0] = '\0';
                break;
            case viamd::EventType_ViamdSystemFree:
            case viamd::EventType_ViamdTrajectoryFree:
            case viamd::EventType_ViamdSystemInit:
            case viamd::EventType_ViamdTrajectoryInit:
            case viamd::EventType_ViamdLoadData:
                clear_data();
                break;
            case viamd::EventType_ViamdShutdown:
                clear_data();
                density_renderer.destroy();
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
                        } else if (str_eq(ident, STR_LIT("MarkerSize"))) {
                            viamd::extract_flt(marker_size, arg);
                            marker_size = std::isfinite(marker_size) ? CLAMP(marker_size, 0.01f, 10.0f) : 1.5f;
                        } else if (str_eq(ident, STR_LIT("UseMass"))) {
                            viamd::extract_bool(use_mass, arg);
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
                        } else if (str_eq(ident, STR_LIT("ColorMode"))) {
                            viamd::extract_enum(mode, arg, ColorMode_Count);
                        } else if (str_eq(ident, STR_LIT("Colormap"))) {
                            viamd::extract_enum(colormap, arg, ImPlot::GetColormapCount());
                        } else if (str_eq(ident, STR_LIT("Color"))) {
                            viamd::extract_flt_vec(&color.x, 4, arg);
                        } else if (str_eq(ident, STR_LIT("PropertyName"))) {
                            viamd::extract_to_char_buf(property_name, sizeof(property_name), arg);
                        }
                    }
                }
                break;
            }
            case viamd::EventType_ViamdSerialize: {
                viamd::serialization_state_t& state = *(viamd::serialization_state_t*)e.payload;
                viamd::write_section_header(state, STR_LIT("ShapeSpace"));
                viamd::write_str(state, STR_LIT("Filter"), str_from_cstr(input));
                viamd::write_flt(state, STR_LIT("MarkerSize"), marker_size);
                viamd::write_bool(state, STR_LIT("UseMass"), use_mass);
                viamd::write_int(state, STR_LIT("Display"), display);
                viamd::write_bool(state, STR_LIT("ShowScatter"), show_scatter);
                viamd::write_flt(state, STR_LIT("Smoothing"), smoothing);
                viamd::write_flt(state, STR_LIT("DensityOpacity"), density_opacity);
                viamd::write_int(state, STR_LIT("DensityColormap"), density_colormap);
                viamd::write_int(state, STR_LIT("ColorMode"), mode);
                viamd::write_int(state, STR_LIT("Colormap"), colormap);
                viamd::write_flt_vec(state, STR_LIT("Color"), &color.x, 4);
                viamd::write_str(state, STR_LIT("PropertyName"), str_from_cstr(property_name));
                break;
            }
            default:
                break;
            }
        }
    }

    void draw_window() {
        if (!show_window) return;

        if (evaluation_pending && !task_system::task_is_running(evaluate_task)) {
            task_system::task_wait_for(evaluate_task);
            evaluation_pending = false;
            if (eval_hash == evaluation_hash() && !evaluation_failed.load()) {
                weights = pending_weights;
                coords = pending_coords;
                density_smoothing = -1.0f;
            } else if (evaluation_failed.load()) {
                snprintf(error, sizeof(error), "Could not read all selected trajectory frames.");
                input_valid = false;
            }
        }

        if (use_source) {
            SeriesVectorView view;
            source_version = series_resolve_vector(&view, app_state, source_key, 3) ? view.version + 1 : 0;
        }

        ImGui::SetNextWindowSize({300,350}, ImGuiCond_FirstUseEver);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(2, 2));
        defer { ImGui::PopStyleVar(1); };

        if (ImGui::Begin("Shapespace", &show_window, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("Export", coords != nullptr)) {
                    if (ImGui::MenuItem("XVG")) {
                        export_shape_space(STR_LIT("xvg"));
                    }
                    if (ImGui::MenuItem("CSV")) {
                        export_shape_space(STR_LIT("csv"));
                    }
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Settings")) {
                    static constexpr float marker_min_size = 0.01f;
                    static constexpr float marker_max_size = 10.0f;
                    const char* color_names[] = {"Uniform", "Colormap"};
                    if (ImGui::BeginCombo("Scatter Color", color_names[mode])) {
                        for (int i = 0; i < (int)ColorMode_Count; ++i) {
                            if (ImGui::Selectable(color_names[i], (int)mode == i)) {
                                mode = (ColorMode)i;
                            }
                        }
                        ImGui::EndCombo();
                    }
                    if (mode == ColorMode_Colormap) {
                        ImPlot::ColormapSelection("Colormap Picker", &colormap);
                    } else if (mode == ColorMode_Uniform) {
                        ImGui::ColorPicker4("Uniform Color", &color.x);
                    }
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
                    ImGui::SliderFloat("Marker Size", &marker_size, marker_min_size, marker_max_size);
                    ImGui::BeginDisabled(use_source);
                    ImGui::Checkbox("Use Mass", &use_mass);
                    ImGui::EndDisabled();
                    ImGui::EndMenu();
                }
                SeriesKey chosen = {};
                const int pick = series_vector_menu(app_state, 3, use_source ? &source_key : nullptr, &chosen);
                if (pick == 1) { source_key = chosen; use_source = true; }
                if (pick == 2) use_source = false;
                ImGui::EndMenuBar();
            }

            SeriesKey dropped = {};
            if (series_vector_drop(app_state, 3, &dropped)) {
                source_key = dropped;
                use_source = true;
            }

            if (use_source) {
                char label[96];
                series_label(label, sizeof(label), app_state, source_key);
                ImGui::AlignTextToFramePadding();
                ImPlot::ItemIcon(series_default_color(app_state, source_key));
                ImGui::SameLine();
                ImGui::TextUnformatted(label);
                ImGui::SameLine();
                if (ImGui::SmallButton("Use selection")) use_source = false;
            } else {
                const float item_width = MAX(ImGui::GetContentRegionAvail().x, 150.f);
                ImGui::PushItemWidth(item_width);
                ImGui::InputQuery("##input", input, sizeof(input), input_valid && eval_hash == evaluation_hash());
                ImGui::PopItemWidth();
                if (ImGui::IsItemHovered()) {
                    if (input_valid && eval_hash == evaluation_hash()) {
                        md_bitfield_copy(&app_state->selection.highlight_mask, &joined_bitfield);
                    } else if (error[0] != '\0') {
                        ImGui::SetTooltip("%s", error);
                    }
                }
            }
            ImGui::SetNextItemWidth(MAX(60.0f, ImGui::GetContentRegionAvail().x - 145.0f));
            ImGui::InputText("##Property Name", property_name, sizeof(property_name));
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Name of the property.");
            ImGui::SameLine();
            ImGui::BeginDisabled(use_source || !input_valid || evaluation_pending || !coords || eval_hash != evaluation_hash());
            if (ImGui::Button("Create Property")) create_property();
            ImGui::EndDisabled();
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
                ImGui::SetTooltip("Add '%s = shapespace(...)' for this selection to the script.", property_name);
            }
            if (property_error[0]) ImGui::TextWrapped("%s", property_error);
            if (evaluation_pending) ImGui::TextUnformatted("Evaluating shape space...");
            if (!input_valid && error[0]) ImGui::TextWrapped("%s", error);

            size_t frame_beg, frame_end;
            frame_range(frame_beg, frame_end);
            if (coords && display != analysis_density::Scatter &&
                (density_smoothing != smoothing || density_beg != frame_beg || density_end != frame_end)) {
                density.update(coords, num_frames, num_structures, frame_beg, frame_end, smoothing);
                density_smoothing = smoothing;
                density_beg = frame_beg;
                density_end = frame_end;
                density_render_dirty = true;
            }
            if (coords && display != analysis_density::Scatter &&
                (density_render_dirty || rendered_display != display || rendered_colormap != density_colormap)) {
                if (density_renderer.render(density, display, density_colormap)) {
                    density_render_dirty = false;
                    rendered_display = display;
                    rendered_colormap = density_colormap;
                }
            }

            ImPlot::PushStyleVar(ImPlotStyleVar_PlotPadding, ImVec2(2,2));
            defer { ImPlot::PopStyleVar(1); };

            const ImPlotFlags flags = ImPlotFlags_Equal | ImPlotFlags_NoMenus | ImPlotFlags_NoMouseText;
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

                const ImVec2 pos_iso = ImPlot::PlotToPixels(iso);
                const ImVec2 pos_pla = ImPlot::PlotToPixels(pla);
                const ImVec2 pos_lin = ImPlot::PlotToPixels(lin);

                static const char* text_sph = "Isotropic";
                static const char* text_pla = "Planar";
                static const char* text_lin = "Linear";

                const ImVec2 text_sph_offset = pos_iso + ImGui::CalcTextSize(text_sph) * ImVec2(-0.5f, -1.2f);
                const ImVec2 text_pla_offset = pos_pla + ImGui::CalcTextSize(text_pla) * ImVec2(-0.5f,  0.2f);
                const ImVec2 text_lin_offset = pos_lin + ImGui::CalcTextSize(text_lin) * ImVec2(-0.5f,  0.2f);

                ImPlot::PushPlotClipRect();
                ImPlot::GetPlotDrawList()->AddTriangleFilled(p0, p1, p2, IM_COL32(255,255,255,20));
                ImPlot::GetPlotDrawList()->AddTriangle(p0, p1, p2, IM_COL32(255,255,255,50));
                ImPlot::PopPlotClipRect();
                if (coords && display != analysis_density::Scatter && !density_render_dirty) {
                    density_renderer.draw(density_opacity);
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

                int     hovered_structure_idx = -1;
                int     hovered_point_idx = -1;
                float   hovered_point_min_dist = FLT_MAX;

                vec2_t mouse_coord = {(float)ImPlot::GetPlotMousePos().x, (float)ImPlot::GetPlotMousePos().y};

                ImPlotRect lim = ImPlot::GetPlotLimits();
                const float scl = (float)lim.X.Max - (float)lim.X.Min;
                const float MAX_D2 = 0.0001f * scl * scl;

                ImPlot::PushStyleVar(ImPlotStyleVar_MarkerSize, marker_size);
                ImPlot::PushStyleVar(ImPlotStyleVar_Marker, ImPlotMarker_Square);
                ImPlot::PushColormap(colormap);
                for (int i = 0; coords && i < (int)num_structures; ++i) {
                    int offset = (int)num_frames * i + (int)frame_beg;
                    int count  = (int)(frame_end - frame_beg);
                    vec2_t* coordinates = coords + offset;
                    char buf[32] = "";
                    if (num_structures == 1) {
                        snprintf(buf, sizeof(buf), "##%i", i+1);
                    } else {
                        snprintf(buf, sizeof(buf), "%i", i+1);
                    }
                    const bool draw_scatter = display == analysis_density::Scatter || show_scatter;
                    if (mode == ColorMode_Uniform) ImPlot::SetNextMarkerStyle(ImPlotMarker_Square, marker_size, color, 1.0f, color);
                    if (!draw_scatter) ImPlot::SetNextMarkerStyle(ImPlotMarker_Square, 0.0f, ImVec4(0,0,0,0), 0.0f, ImVec4(0,0,0,0));
                    ImPlot::PlotScatterG(buf, getter, coordinates, count);

                    ImPlotItem* item = ImPlot::GetItem(buf);
                    if (item) {
                        if (item->LegendHovered) {
                            hovered_structure_idx = i;
                        }
                        ImPlot::PopColormap();

                        if (item->Show) {
                            for (int32_t j = offset; j < offset + count; ++j ) {
                                vec2_t delta = mouse_coord - coords[j];
                                float d2 = vec2_dot(delta, delta);
                                if (d2 < MAX_D2 && d2 < hovered_point_min_dist) {
                                    hovered_point_min_dist = d2;
                                    hovered_point_idx = j;
                                }
                            }
                        }
                    }
                }
                ImPlot::PopStyleVar(2);

                if (ImPlot::IsPlotHovered()) {
                    md_bitfield_clear(&app_state->selection.highlight_mask);
                }

                // Redraw hovered index
                ImPlot::PushStyleVar(ImPlotStyleVar_Marker, ImPlotMarker_Square);
                ImPlot::PushStyleVar(ImPlotStyleVar_MarkerSize, marker_size * 1.1f);
                ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, ImVec4(1,1,1,1));
                if (hovered_structure_idx != -1) {
                    int offset = (int)num_frames * hovered_structure_idx + (int)frame_beg;
                    int count  = (int)(frame_end - frame_beg);
                    vec2_t* coordinates = coords + offset;
                    ImPlot::PlotScatterG("##hovered structure", getter, coordinates, count);
                    if (hovered_structure_idx < (int)md_array_size(bitfields)) {
                        md_bitfield_copy(&app_state->selection.highlight_mask, &bitfields[hovered_structure_idx]);
                    }
                }
                if (hovered_point_idx != -1) {
                    vec2_t* coordinates = coords + hovered_point_idx;
                    ImPlot::PlotScatterG("##hovered idx", getter, coordinates, 1);
                    char buf[128] = "";
                    int len = 0;
                    int structure_idx = hovered_point_idx / (int)num_frames;
                    int frame_idx = hovered_point_idx % num_frames;
                    vec3_t w = weights[hovered_point_idx];
                    if (num_structures > 1) {
                        len += snprintf(buf, sizeof(buf), "Structure: %i, ", structure_idx + 1);
                    }
                    snprintf(buf + len, sizeof(buf) - len, "Frame: %i, lin: %.2f, plan: %.2f, iso: %.2f", frame_idx, w.x, w.y, w.z);
                    ImGui::SetTooltip("%s", buf);

                    if (structure_idx < (int)md_array_size(bitfields)) {
                        md_bitfield_copy(&app_state->selection.highlight_mask, &bitfields[structure_idx]);
                    }

                    if (ImGui::IsWindowFocused() && ImPlot::IsPlotHovered() && ImGui::GetIO().MouseClicked[0]) {
                        app_state->animation.frame = frame_idx;
                    }
                }
                ImPlot::PopStyleColor();
                ImPlot::PopStyleVar(2);

                ImPlot::EndPlot();
            }
        }
        ImGui::End();

        // Create a hash from everything which dictates the app_state of an evaluation to compare against
        uint64_t hash = evaluation_hash();
        if (hash != eval_hash) {
            if (task_system::task_is_running(evaluate_task)) {
                task_system::task_interrupt (evaluate_task);
            } else {
                clear_data();
                hash = evaluation_hash();
                eval_hash = hash;

                MEMSET(error, 0, sizeof(error));
                if (use_source) {
                    SeriesVectorView view;
                    if (!series_resolve_vector(&view, app_state, source_key, 3) || !view.num_frames || !view.count) {
                        snprintf(error, sizeof(error), "The property has no values.");
                        return;
                    }
                    // The values are frame major, the arrays here structure major
                    num_frames = view.num_frames;
                    num_structures = view.count;
                    md_array_resize(weights, num_frames * num_structures, arena);
                    md_array_resize(coords,  num_frames * num_structures, arena);
                    const vec2_t p[3] = {{0.0f, 0.0f}, {1.0f, 0.0f}, {0.5f, 0.86602540378f}};
                    for (size_t s = 0; s < num_structures; ++s) {
                        for (size_t f = 0; f < num_frames; ++f) {
                            const float* v = view.data + (f * num_structures + s) * 3;
                            const size_t dst = num_frames * s + f;
                            weights[dst] = vec3_t{v[0], v[1], v[2]};
                            // A frame that has not been evaluated holds zeros, not a shape
                            const bool evaluated = v[0] != 0.0f || v[1] != 0.0f || v[2] != 0.0f;
                            coords[dst] = evaluated ? p[0] * v[0] + p[1] * v[1] + p[2] * v[2] : vec2_t{NAN, NAN};
                        }
                    }
                    input_valid = true;
                    density_smoothing = -1.0f;
                    return;
                }
                if (md_filter_evaluate(&bitfields, str_from_cstr(input), &app_state->mold.sys, &app_state->mold.state, app_state->script.ir, NULL, error, sizeof(error), arena)) {
                    num_structures = md_array_size(bitfields);
                    if (!num_structures) {
                        snprintf(error, sizeof(error), "The selection contains no structures.");
                        return;
                    }
                    num_frames = run_num_frames(app_state);
                    if (!num_frames) {
                        snprintf(error, sizeof(error), "No trajectory frames are available.");
                        return;
                    }
                    for (size_t i = 0; i < num_structures; ++i) {
                        if (!md_bitfield_popcount(&bitfields[i])) {
                            snprintf(error, sizeof(error), "Selection group %zu contains no atoms.", i + 1);
                            return;
                        }
                        md_bitfield_or_inplace(&joined_bitfield, &bitfields[i]);
                    }
                    input_valid = true;
                    evaluation_pending = true;
                    evaluation_failed.store(false);
                    md_array_resize(pending_weights, num_frames * num_structures, arena);
                    md_array_resize(pending_coords,  num_frames * num_structures, arena);
                    const bool evaluate_mass = use_mass;
                    evaluate_task = task_system::create_pool_task(STR_LIT("Eval Shape Space"), (uint32_t)num_frames, [shapespace = this, evaluate_mass](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                        (void)thread_num;
                        ApplicationState* app_state = shapespace->app_state;
                        const size_t stride = ALIGN_TO(app_state->mold.sys.atom.count, 8);
                        const size_t bytes = stride * sizeof(vec3_t);
                        md_temp_scope_t temp_scope = md_temp_begin();
                        md_allocator_i* alloc = md_temp_allocator(temp_scope);
                        defer { md_temp_end(temp_scope); };

                        vec3_t* xyz = (vec3_t*)md_temp_alloc(temp_scope, bytes);
                        float* w = 0;
                        if (evaluate_mass) {
                            w = (float*)md_temp_alloc(temp_scope, stride * sizeof(float));
                            md_atom_extract_masses(w, 0, app_state->mold.sys.atom.count, &app_state->mold.sys.atom);
                        }

                        const vec2_t p[3] = {{0.0f, 0.0f}, {1.0f, 0.0f}, {0.5f, 0.86602540378f}};

                        // One context for the range, so the run's files stay open across its frames.
                        const str_t paths[] = { STR_INIT("atom/position"), STR_INIT("unitcell") };
                        md_system_extract_t* ex = md_system_extract_begin(&app_state->mold.sys, str_from_cstr(app_state->mold.run), paths, ARRAY_SIZE(paths), md_get_heap_allocator());
                        if (!ex) {
                            shapespace->evaluation_failed.store(true);
                            return;
                        }
                        defer { md_system_extract_end(ex); };

                        md_array(vec4_t) xyzw = 0;
                        for (uint32_t frame_idx = range_beg; frame_idx < range_end; ++frame_idx) {
                            md_system_state_t frame_state = {};
                            frame_state.num_atoms = app_state->mold.sys.atom.count;
                            frame_state.xyz = xyz;
                            if (!md_system_extract_frame(ex, frame_idx, &frame_state)) {
                                shapespace->evaluation_failed.store(true);
                                continue;
                            }

                            for (size_t i = 0; i < md_array_size(shapespace->bitfields); ++i) {
                                const md_bitfield_t* bf = &shapespace->bitfields[i];
                                size_t count = md_bitfield_popcount(bf);
                                md_array_resize(xyzw, count, alloc);

                                md_bitfield_iter_t iter = md_bitfield_iter_create(bf);
                                size_t dst_idx = 0;
                                while (md_bitfield_iter_next(&iter)) {
                                    const size_t src_idx = md_bitfield_iter_idx(&iter);
                                    xyzw[dst_idx++] = vec4_from_vec3(xyz[src_idx], w ? w[src_idx] : 1.0f);
                                }

                                // Covariance about the mean of the placed points, not the circular mean the
                                // placement is seeded from (their offset enters the covariance squared)
                                vec3_t com = vec3_zero();
                                md_util_deperiodize_self_vec4(xyzw, count, &frame_state.unitcell, &com);

                                const mat3_t M = mat3_covariance_matrix_vec4(xyzw, 0, count, com);
                                const vec3_t weights = md_util_shape_weights(&M);

                                dst_idx = shapespace->num_frames * i + frame_idx;
                                if (!std::isfinite(weights.x) || !std::isfinite(weights.y) || !std::isfinite(weights.z)) {
                                    shapespace->evaluation_failed.store(true);
                                }
                                shapespace->pending_weights[dst_idx] = weights;
                                shapespace->pending_coords[dst_idx] = p[0] * weights[0] + p[1] * weights[1] + p[2] * weights[2];
                            }
                        }
                    });

                    if (!evaluate_task) {
                        evaluation_pending = false;
                        input_valid = false;
                        snprintf(error, sizeof(error), "Could not schedule shape-space evaluation.");
                        return;
                    }
                    task_system::enqueue_task(evaluate_task);
                }
            }
        }
    }

    void export_shape_space(str_t ext) {
        char path_buf[2048];
        if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, ext)) {
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

            size_t num_structs = num_structures;
            size_t num_rows    = num_frames;

            for (size_t i = 0; i < num_structs; ++i) {
                md_array_push(column_labels, str_printf(temp_arena, "linear",  i + 1), temp_arena);
                md_array_push(column_values, NULL, temp_arena);
                md_array_push(column_labels, str_printf(temp_arena, "planar", i + 1), temp_arena);
                md_array_push(column_values, NULL, temp_arena);
                md_array_push(column_labels, str_printf(temp_arena, "isotropic",  i + 1), temp_arena);
                md_array_push(column_values, NULL, temp_arena);
            }

            const float* time = app_state->timeline.x_values;

            for (size_t i = 0; i < num_rows; ++i) {
                float t = time ? time[i] : (float)i;
                md_array_push(column_values[0], t, temp_arena);

                for (size_t j = 0; j < num_structs; ++j) {
                    size_t idx = num_rows * j + i;
                    vec3_t w = weights[idx];
                    md_array_push(column_values[1 + j * 3 + 0], w.x, temp_arena);
                    md_array_push(column_values[1 + j * 3 + 1], w.y, temp_arena);
                    md_array_push(column_values[1 + j * 3 + 2], w.z, temp_arena);
                }
            }

            size_t num_cols = num_structs * 3 + 1;
            ASSERT(num_cols == md_array_size(column_values));
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
    }
};

static Shapespace instance = {};
