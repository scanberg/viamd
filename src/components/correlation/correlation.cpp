#include <event.h>

#include <core/md_log.h>
#include <core/md_array.h>
#include <core/md_vec_math.h>
#include <core/md_bitfield.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

#include <md_script.h>
#include <md_util.h>

#include <viamd.h>
#include <serialization_utils.h>
#include <imgui_widgets.h>
#include <implot_widgets.h>
#include <implot_internal.h>

#include <string.h>
#include <stdlib.h>

struct PropertyInfo {
    char name[64] = "";
    int index = 0;
    size_t num_values = 0;
    size_t num_frames = 0;
    const float* data = nullptr;
    bool is_array = false;
};

struct ScatterSeries {
    md_array(float) x_data = 0;
    md_array(float) y_data = 0;
    md_array(int) frame_indices = 0;
    char name[64] = "";
    ImVec4 color = {1.0f, 1.0f, 1.0f, 1.0f};
};

struct Correlation : viamd::EventHandler {
    char error[256] = "";
    bool show_window = false;
    
    // Property selection
    int x_property_idx = -1;
    int y_property_idx = -1;
    int x_series_idx = 0;  // For array properties
    int y_series_idx = 0;  // For array properties
    
    // Available properties
    md_array(PropertyInfo) properties = 0;
    
    // Scatter plot data
    md_array(ScatterSeries) series = 0;
    
    // Interaction
    int hovered_point = -1;
    int clicked_frame = -1;
    
    md_allocator_i* arena = 0;
    ApplicationState* app_state = 0;

    Correlation() { viamd::event_system_register_handler(*this); }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& e = events[i];

            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                app_state = (ApplicationState*)e.payload;
                arena = md_arena_allocator_create(app_state->allocator.persistent, MEGABYTES(1));
                break;
            }
            case viamd::EventType_ViamdShutdown:
                md_arena_allocator_destroy(arena);
                break;
            case viamd::EventType_ViamdFrameTick:
                draw_window();
                break;
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("Correlation", &show_window);
                break;
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t& state = *(viamd::deserialization_state_t*)e.payload;
                if (str_eq(viamd::section_header(state), STR_LIT("Correlation"))) {
                    str_t ident, arg;
                    while (viamd::next_entry(ident, arg, state)) {
                        if (str_eq(ident, STR_LIT("show_window"))) {
                            viamd::extract_bool(show_window, arg);
                        }
                    }
                }
                break;
            }
            case viamd::EventType_ViamdSerialize: {
                viamd::serialization_state_t& state = *(viamd::serialization_state_t*)e.payload;
                viamd::write_section_header(state, STR_LIT("Correlation"));
                viamd::write_bool(state, STR_LIT("show_window"), show_window);
                break;
            }
            default:
                break;
            }
        }
    }

    void update_properties() {
        if (!app_state || !app_state->script.eval_ir) return;
        
        md_array_resize(properties, 0, arena);
        
        const md_script_ir_t* ir = app_state->script.eval_ir;
        const size_t num_properties = md_script_ir_property_count(ir);
        const str_t* property_names = md_script_ir_property_names(ir);
        
        for (size_t i = 0; i < num_properties; ++i) {
            PropertyInfo info = {};
            str_copy_to_char_buf(info.name, sizeof(info.name), property_names[i]);
            info.index = (int)i;
            
            // Get property data from evaluation
            if (app_state->script.full_eval) {
                const md_script_property_data_t* eval_data = md_script_eval_property_data(app_state->script.full_eval, property_names[i]);
                if (eval_data) {
                    info.num_frames = eval_data->dim[0]; // Temporal dimension
                    info.num_values = eval_data->num_values;
                    info.data = eval_data->values;
                    info.is_array = (eval_data->num_values > (size_t)eval_data->dim[0]);
                }
            }
            
            md_array_push(properties, info, arena);
        }
    }
    
    void update_scatter_data() {
        if (x_property_idx < 0 || y_property_idx < 0 || 
            x_property_idx >= (int)md_array_size(properties) || 
            y_property_idx >= (int)md_array_size(properties)) {
            return;
        }
        
        const PropertyInfo& x_prop = properties[x_property_idx];
        const PropertyInfo& y_prop = properties[y_property_idx];
        
        if (!x_prop.data || !y_prop.data || x_prop.num_frames != y_prop.num_frames) {
            return;
        }
        
        // Clear existing series
        for (size_t i = 0; i < md_array_size(series); ++i) {
            md_array_free(series[i].x_data, arena);
            md_array_free(series[i].y_data, arena);
            md_array_free(series[i].frame_indices, arena);
        }
        md_array_resize(series, 0, arena);
        
        const size_t num_frames = x_prop.num_frames;
        
        if (x_prop.is_array || y_prop.is_array) {
            // Handle array properties - create multiple series
            const size_t x_values_per_frame = x_prop.is_array ? (x_prop.num_values / x_prop.num_frames) : 1;
            const size_t y_values_per_frame = y_prop.is_array ? (y_prop.num_values / y_prop.num_frames) : 1;
            
            const size_t max_series = ImMax(x_values_per_frame, y_values_per_frame);
            
            for (size_t s = 0; s < max_series; ++s) {
                ScatterSeries scatter = {};
                snprintf(scatter.name, sizeof(scatter.name), "%s[%zu] vs %s[%zu]", 
                    x_prop.name, x_prop.is_array ? s : 0,
                    y_prop.name, y_prop.is_array ? s : 0);
                
                // Assign a color based on series index
                scatter.color = ImPlot::GetColormapColor((int)s);
                
                for (size_t f = 0; f < num_frames; ++f) {
                    float x_val, y_val;
                    
                    if (x_prop.is_array) {
                        if (s < x_values_per_frame) {
                            x_val = x_prop.data[f * x_values_per_frame + s];
                        } else {
                            continue; // Skip this series if x doesn't have enough values
                        }
                    } else {
                        x_val = x_prop.data[f];
                    }
                    
                    if (y_prop.is_array) {
                        if (s < y_values_per_frame) {
                            y_val = y_prop.data[f * y_values_per_frame + s];
                        } else {
                            continue; // Skip this series if y doesn't have enough values
                        }
                    } else {
                        y_val = y_prop.data[f];
                    }
                    
                    md_array_push(scatter.x_data, x_val, arena);
                    md_array_push(scatter.y_data, y_val, arena);
                    md_array_push(scatter.frame_indices, (int)f, arena);
                }
                
                if (md_array_size(scatter.x_data) > 0) {
                    md_array_push(series, scatter, arena);
                }
            }
        } else {
            // Simple case - single series
            ScatterSeries scatter = {};
            snprintf(scatter.name, sizeof(scatter.name), "%s vs %s", x_prop.name, y_prop.name);
            scatter.color = ImPlot::GetColormapColor(0);
            
            for (size_t f = 0; f < num_frames; ++f) {
                md_array_push(scatter.x_data, x_prop.data[f], arena);
                md_array_push(scatter.y_data, y_prop.data[f], arena);
                md_array_push(scatter.frame_indices, (int)f, arena);
            }
            
            md_array_push(series, scatter, arena);
        }
    }

    void draw_window() {
        if (!show_window) return;
        
        if (ImGui::Begin("Correlation Plot", &show_window)) {
            update_properties();
            
            // Property selection UI
            ImGui::Text("X-Axis Property:");
            ImGui::SameLine();
            if (ImGui::BeginCombo("##x_property", x_property_idx >= 0 ? properties[x_property_idx].name : "Select...")) {
                for (size_t i = 0; i < md_array_size(properties); ++i) {
                    bool is_selected = (x_property_idx == (int)i);
                    if (ImGui::Selectable(properties[i].name, is_selected)) {
                        x_property_idx = (int)i;
                    }
                    if (is_selected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }
            
            ImGui::Text("Y-Axis Property:");
            ImGui::SameLine();
            if (ImGui::BeginCombo("##y_property", y_property_idx >= 0 ? properties[y_property_idx].name : "Select...")) {
                for (size_t i = 0; i < md_array_size(properties); ++i) {
                    bool is_selected = (y_property_idx == (int)i);
                    if (ImGui::Selectable(properties[i].name, is_selected)) {
                        y_property_idx = (int)i;
                    }
                    if (is_selected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }
            
            if (x_property_idx >= 0 && y_property_idx >= 0) {
                if (ImGui::Button("Generate Plot")) {
                    update_scatter_data();
                }
                
                // Display scatter plot
                if (md_array_size(series) > 0) {
                    if (ImPlot::BeginPlot("Property Correlation", ImVec2(-1, -1))) {
                        
                        for (size_t s = 0; s < md_array_size(series); ++s) {
                            const ScatterSeries& scatter = series[s];
                            if (md_array_size(scatter.x_data) > 0) {
                                ImPlot::PushStyleColor(ImPlotCol_MarkerFill, scatter.color);
                                ImPlot::PlotScatter(scatter.name, 
                                    scatter.x_data, scatter.y_data, 
                                    (int)md_array_size(scatter.x_data));
                                ImPlot::PopStyleColor();
                                
                                // Check for hover/click on points
                                if (ImPlot::IsPlotHovered()) {
                                    // Find closest point in screen space
                                    float min_dist_sq = FLT_MAX;
                                    int closest_point = -1;
                                    
                                    for (size_t p = 0; p < md_array_size(scatter.x_data); ++p) {
                                        // Convert plot coordinates to screen space for distance calculation
                                        ImVec2 screen_pos = ImPlot::PlotToPixels(scatter.x_data[p], scatter.y_data[p]);
                                        ImPlotPoint mouse_plot = ImPlot::GetPlotMousePos();
                                        ImVec2 mouse_pixel = ImPlot::PlotToPixels(mouse_plot.x, mouse_plot.y);
                                        
                                        float dx = screen_pos.x - mouse_pixel.x;
                                        float dy = screen_pos.y - mouse_pixel.y;
                                        float dist_sq = dx * dx + dy * dy;
                                        
                                        if (dist_sq < min_dist_sq) {
                                            min_dist_sq = dist_sq;
                                            closest_point = (int)p;
                                        }
                                    }
                                    
                                    // Check if close enough to hover (within 10 pixels)
                                    if (closest_point >= 0 && min_dist_sq < 100.0f) {
                                        hovered_point = closest_point;
                                        
                                        // Show tooltip
                                        ImGui::BeginTooltip();
                                        ImGui::Text("Frame: %d", scatter.frame_indices[closest_point]);
                                        ImGui::Text("X: %.3f", scatter.x_data[closest_point]);
                                        ImGui::Text("Y: %.3f", scatter.y_data[closest_point]);
                                        ImGui::Text("Click to jump to this frame");
                                        ImGui::EndTooltip();
                                        
                                        // Handle click to jump to frame
                                        if (ImGui::IsMouseClicked(0)) {
                                            clicked_frame = scatter.frame_indices[closest_point];
                                            app_state->animation.frame = (double)clicked_frame;
                                        }
                                    }
                                }
                            }
                        }
                        
                        ImPlot::EndPlot();
                    }
                }
            }
            
            // Show error messages if any
            if (strlen(error) > 0) {
                ImGui::TextColored(ImVec4(1, 0, 0, 1), "Error: %s", error);
            }
        }
        ImGui::End();
    }
};

static Correlation instance = {};