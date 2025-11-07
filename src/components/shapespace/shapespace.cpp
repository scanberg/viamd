#define IMGUI_DEFINE_MATH_OPERATORS

#include <event.h>

#include <core/md_log.h>
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

#include <viamd.h>
#include <serialization_utils.h>
#include <imgui_widgets.h>
#include <implot_widgets.h>
#include <implot_internal.h>
#include <loader.h>

/*
    @NOTE:
    It would be nice to have a colormap applied to the points in the plot where they are graded by t=time
    However, the current implot api for scatterplot does not currently support that.
    Bummer!
*/

enum ColorMode {
    ColorMode_Uniform,
    ColorMode_Colormap,
    ColorMode_Count,
};

static const char* colormode_lbl[] = {
    "Uniform",
    "Colormap",
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
    md_array(md_bitfield_t) bitfields = 0;
    md_bitfield_t joined_bitfield = {};

    md_allocator_i* arena = 0;

    // Settings
    float marker_size = 1.5f;
    bool  use_mass = true;
    ColorMode mode = ColorMode_Colormap;
    ImPlotColormap colormap = 0;
    ImVec4 color = {0.5f,1.0f,0.5f,1.0f};

    task_system::ID evaluate_task = 0;

    ApplicationState* app_state = 0;

    Shapespace() { viamd::event_system_register_handler(*this); }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& e = events[i];

            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                app_state = (ApplicationState*)e.payload;
                arena = md_arena_allocator_create(app_state->allocator.persistent, MEGABYTES(1));
                md_bitfield_init(&joined_bitfield, arena);
                break;
            }
            case viamd::EventType_ViamdShutdown:
                task_system::task_interrupt_and_wait_for(evaluate_task);
                md_arena_allocator_destroy(arena);
                break;
            case viamd::EventType_ViamdFrameTick:
                draw_window();
                break;
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("ShapeSpace", &show_window);
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
                        } else if (str_eq(ident, STR_LIT("UseMass"))) {
                            viamd::extract_bool(use_mass, arg);
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
                break;
            }
            default:
                break;
            }
        }
    }

    void draw_window() {
        if (!show_window) return;

        ImGui::SetNextWindowSize({300,350}, ImGuiCond_FirstUseEver);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(2, 2));
        defer { ImGui::PopStyleVar(1); };

        if (ImGui::Begin("Shapespace", &show_window, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("Export")) {
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
#if 0
                    if (ImGui::BeginCombo("Color", colormode_lbl[mode])) {
                        for (int i = 0; i < (int)ColorMode_Count; ++i) {
                            if (ImGui::Selectable(colormode_lbl[i], (int)mode == i)) {
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
#endif
                    ImGui::SliderFloat("Marker Size", &marker_size, marker_min_size, marker_max_size);
                    ImGui::Checkbox("Use Mass", &use_mass);
                    ImGui::EndMenu();
                }
                ImGui::EndMenuBar();
            }

            const float item_width = MAX(ImGui::GetContentRegionAvail().x, 150.f);
            ImGui::PushItemWidth(item_width);
            ImGui::InputQuery("##input", input, sizeof(input), input_valid);
            ImGui::PopItemWidth();
            if (ImGui::IsItemHovered()) {
                if (input_valid) {
                    md_bitfield_copy(&app_state->selection.highlight_mask, &joined_bitfield);
                } else if (error[0] != '\0') {
                    ImGui::SetTooltip("%s", error);
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

                static const char* text_sph = "Spherical";
                static const char* text_pla = "Planar";
                static const char* text_lin = "Linear";

                const ImVec2 text_sph_offset = pos_iso + ImGui::CalcTextSize(text_sph) * ImVec2(-0.5f, -1.2f);
                const ImVec2 text_pla_offset = pos_pla + ImGui::CalcTextSize(text_pla) * ImVec2(-0.5f,  0.2f);
                const ImVec2 text_lin_offset = pos_lin + ImGui::CalcTextSize(text_lin) * ImVec2(-0.5f,  0.2f);

                ImPlot::PushPlotClipRect();
                ImPlot::GetPlotDrawList()->AddTriangleFilled(p0, p1, p2, IM_COL32(255,255,255,20));
                ImPlot::GetPlotDrawList()->AddTriangle(p0, p1, p2, IM_COL32(255,255,255,50));

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
                for (int i = 0; i < (int)num_structures; ++i) {
                    int offset = (int)num_frames * i;
                    int count  = (int)num_frames;
                    if (app_state->timeline.filter.enabled) {
                        offset += (int)app_state->timeline.filter.beg_frame;
                        count = (int)MAX(app_state->timeline.filter.end_frame - app_state->timeline.filter.beg_frame, 0);
                    }
                    vec2_t* coordinates = coords + offset;
                    char buf[32] = "";
                    if (num_structures == 1) {
                        snprintf(buf, sizeof(buf), "##%i", i+1);
                    } else {
                        snprintf(buf, sizeof(buf), "%i", i+1);
                    }
                    ImPlot::PlotScatterG(buf, getter, coordinates, count);

                    ImPlotItem* item = ImPlot::GetItem(buf);
                    if (item) {
                        if (item->LegendHovered) {
                            hovered_structure_idx = i;
                        }

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
                    int offset = (int)num_frames * hovered_structure_idx;
                    int count  = (int)num_frames;
                    if (app_state->timeline.filter.enabled) {
                        offset += (int)app_state->timeline.filter.beg_frame;
                        count   = MAX((int)(app_state->timeline.filter.end_frame - app_state->timeline.filter.beg_frame), 0);
                    }
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
                    len += snprintf(buf + len, sizeof(buf) - len, "Frame: %i, lin: %.2f, plan: %.2f, iso: %.2f", frame_idx, w.x, w.y, w.z);
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
        uint64_t hash = md_hash64(input, sizeof(input), app_state->script.ir_fingerprint ^ (1ULL << (uint64_t)use_mass));
        if (hash != eval_hash) {
            if (task_system::task_is_running(evaluate_task)) {
                task_system::task_interrupt (evaluate_task);
            } else {
                bitfields = 0;
                weights = 0;
                coords = 0;
                num_frames = 0;
                num_structures = 0;
                md_arena_allocator_reset(arena);
                joined_bitfield = {0};
                md_bitfield_init(&joined_bitfield, arena);

				md_trajectory_i* traj = load::traj::get_raw_trajectory(app_state->mold.traj);

                input_valid = false;
                MEMSET(error, 0, sizeof(error));
                if (md_filter_evaluate(&bitfields, str_from_cstr(input), &app_state->mold.sys, app_state->mold.sys.atom.x, app_state->mold.sys.atom.y, app_state->mold.sys.atom.z, app_state->script.ir, NULL, error, sizeof(error), arena)) {
                    eval_hash = hash;

                    input_valid = true;                    
                    num_structures = md_array_size(bitfields);
                    if (!num_structures) {
                        MD_LOG_ERROR("No structures present when attempting to populate shape space");
                        return;
                    }
                    num_frames = md_trajectory_num_frames(traj);
                    if (!num_frames) {
                        MD_LOG_ERROR("No trajectory frames present when attempting to populate shape space");
                        return;
                    }
                    for (size_t i = 0; i < num_structures; ++i) {
                        md_bitfield_or_inplace(&joined_bitfield, &bitfields[i]);
                    }

                    md_array_resize(weights, num_frames * num_structures, arena);
                    md_array_resize(coords,  num_frames * num_structures, arena);
                    MEMSET(weights, 0, md_array_bytes(weights));
                    MEMSET(coords,  0, md_array_bytes(coords));
                    evaluate_task = task_system::create_pool_task(STR_LIT("Eval Shape Space"), (uint32_t)num_frames, [shapespace = this, traj](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                        (void)thread_num;
                        ApplicationState* app_state = shapespace->app_state;
                        const size_t stride = ALIGN_TO(app_state->mold.sys.atom.count, 8);
                        const size_t bytes = stride * 3 * sizeof(float);
                        md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
                        defer { md_arena_allocator_destroy(alloc); };

                        float* coords = (float*)md_arena_allocator_push(alloc, bytes);
                        float* x = coords + stride * 0;
                        float* y = coords + stride * 1;
                        float* z = coords + stride * 2;
                        float* w = 0;
                        if (shapespace->use_mass) {
                            w = (float*)md_arena_allocator_push(alloc, stride * sizeof(float));
                            md_atom_extract_masses(w, 0, app_state->mold.sys.atom.count, &app_state->mold.sys.atom);
                        }

                        const vec2_t p[3] = {{0.0f, 0.0f}, {1.0f, 0.0f}, {0.5f, 0.86602540378f}};

                        md_array(vec4_t) xyzw = 0;
                        for (uint32_t frame_idx = range_beg; frame_idx < range_end; ++frame_idx) {
                            md_trajectory_frame_header_t header;
                            md_trajectory_load_frame(traj, frame_idx, &header, x, y, z);

                            for (size_t i = 0; i < md_array_size(shapespace->bitfields); ++i) {
                                const md_bitfield_t* bf = &shapespace->bitfields[i];
                                size_t count = md_bitfield_popcount(bf);
                                md_array_resize(xyzw, count, alloc);

                                md_bitfield_iter_t iter = md_bitfield_iter_create(bf);
                                size_t dst_idx = 0;
                                while (md_bitfield_iter_next(&iter)) {
                                    const size_t src_idx = md_bitfield_iter_idx(&iter);
                                    xyzw[dst_idx++] = vec4_set(x[src_idx], y[src_idx], z[src_idx], w ? w[src_idx] : 1.0f);
                                }

                                vec3_t com = md_util_com_compute_vec4(xyzw, 0, count, &app_state->mold.sys.unitcell);
                                md_util_deperiodize_vec4(xyzw, count, com, &app_state->mold.sys.unitcell);

                                const mat3_t M = mat3_covariance_matrix_vec4(xyzw, 0, count, com);
                                const vec3_t weights = md_util_shape_weights(&M);

                                dst_idx = shapespace->num_frames * i + frame_idx;
                                shapespace->weights[dst_idx] = weights;
                                shapespace->coords[dst_idx] = p[0] * weights[0] + p[1] * weights[1] + p[2] * weights[2];
                            }
                        }
                    });

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

            md_allocator_i* temp_arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
            defer { md_arena_allocator_destroy(temp_arena); };

            // @TODO: add unit to time (if available)
            md_unit_t time_unit = md_trajectory_time_unit(app_state->mold.traj);

            str_t x_label = STR_LIT("Frame");
            if (!md_unit_empty(time_unit)) {
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
                str_t title = STR_LIT("Shape Space");
                str_t y_label = STR_LIT("Shape Weight Ratio");
                str_t header = md_xvg_format_header(title, x_label, y_label, md_array_size(column_labels) - 1, column_labels + 1, temp_arena);
                str_t xvg    = md_xvg_format(header, num_cols, num_rows, column_values, temp_arena);
                md_file_o* file = md_file_open(path, MD_FILE_WRITE);
                md_file_write(file, xvg.ptr, xvg.len);
            } else {
                MD_LOG_DEBUG("Unrecognized export format");
                ASSERT(false);
            }
        }
    }
};

static Shapespace instance = {};
