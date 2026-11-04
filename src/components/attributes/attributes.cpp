// Attributes window: every attribute the application holds, in one table.
//
// The system's table (topology columns, the run and everything loaded along it) and the script
// evaluations' tables (script/<ident> and its siblings) are listed side by side, one row per
// attribute, with what it is, how it is stored and whether the Timelines and Distributions windows
// can show it. A row that can be shown is dragged from here into a subplot of either window, and one
// of the system that the script reads (attr()) into the Script Editor, where it goes in as the
// expression that reads it.
//
// Opened with KEY_SHOW_ATTRIBUTE_WINDOW (Ctrl+Shift+A) from anywhere, or from the Windows menu.

#include <viamd_event.h>
#include <viamd.h>
#include <event.h>
#include <plot_series.h>
#include <timeline_periods_ui.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_str.h>
#include <md_system.h>

#include <imgui.h>
#include <imgui_widgets.h>
#include <implot.h>
#include <app/IconsFontAwesome6.h>

#include <stdio.h>
#include <string.h>

namespace attributes {

static const char* type_name(md_attribute_type_t type) {
    switch (type) {
    case MD_ATTRIBUTE_TYPE_F32: return "f32";
    case MD_ATTRIBUTE_TYPE_F64: return "f64";
    case MD_ATTRIBUTE_TYPE_I8:  return "i8";
    case MD_ATTRIBUTE_TYPE_U8:  return "u8";
    case MD_ATTRIBUTE_TYPE_I16: return "i16";
    case MD_ATTRIBUTE_TYPE_U16: return "u16";
    case MD_ATTRIBUTE_TYPE_I32: return "i32";
    case MD_ATTRIBUTE_TYPE_U32: return "u32";
    case MD_ATTRIBUTE_TYPE_I64: return "i64";
    case MD_ATTRIBUTE_TYPE_U64: return "u64";
    case MD_ATTRIBUTE_TYPE_STR: return "str";
    default: return "?";
    }
}

static const char* storage_name(const md_attribute_t& attr) {
    if (md_attribute_is_alias(&attr))   return "alias";
    if (md_attribute_is_virtual(&attr)) return "virtual";
    return "resident";
}

// "System", "Script", "Script (Production)": the script over a timeline period
static void source_label(char* buf, size_t cap, const ApplicationState& state, SeriesSource source, uint32_t period) {
    const char* name = source == SeriesSource_System ? "System" : source == SeriesSource_Script ? "Script" : "?";
    if (period != 0) {
        snprintf(buf, cap, "%s (%s)", name, series_period_name(&state, period));
    } else {
        snprintf(buf, cap, "%s", name);
    }
}

// "f32 [501 x 2]", "f64 [3001] x3", "str"
static void format_string(char* buf, size_t cap, const md_attribute_format_t& fmt) {
    int len = snprintf(buf, cap, "%s", type_name(fmt.type));
    if (fmt.rank > 0 && len > 0 && (size_t)len < cap) {
        len += snprintf(buf + len, cap - len, " [");
        for (uint32_t i = 0; i < fmt.rank && (size_t)len < cap; ++i) {
            len += snprintf(buf + len, cap - len, i ? " x %u" : "%u", fmt.shape[i]);
        }
        if ((size_t)len < cap) len += snprintf(buf + len, cap - len, "]");
    }
    if (fmt.components > 1 && (size_t)len < cap) {
        snprintf(buf + len, cap - len, " x%u", fmt.components);
    }
}

// ASCII case insensitive 'needle occurs in haystack'
static bool contains_ignore_case(str_t haystack, str_t needle) {
    if (needle.len == 0) return true;
    if (needle.len > haystack.len) return false;
    for (size_t i = 0; i + needle.len <= haystack.len; ++i) {
        if (str_eq_ignore_case(str_substr(haystack, i, needle.len), needle)) return true;
    }
    return false;
}

// An attribute of a distribution or a volume of the script ("script/<ident>", "script/<ident>/weight")
static bool script_aggregate_path(const ApplicationState& state, str_t path) {
    const str_t prefix = STR_LIT("script/");
    if (!state.script.eval_ir || !str_begins_with(path, prefix)) return false;
    str_t ident = str_substr(path, prefix.len);
    size_t loc;
    if (str_find_char(&loc, ident, '/')) ident = str_substr(ident, 0, loc);
    return script::output_flags(state.script.eval_ir, ident) & (script::OUTPUT_DISTRIBUTION | script::OUTPUT_VOLUME);
}

// Every whitespace separated word of the filter occurs in the path, the label or the source
static bool matches_filter(str_t filter, str_t path, str_t label, const char* source) {
    const str_t src = str_from_cstr(source);
    size_t i = 0;
    while (i < filter.len) {
        while (i < filter.len && (filter.ptr[i] == ' ' || filter.ptr[i] == '\t')) ++i;
        size_t j = i;
        while (j < filter.len && filter.ptr[j] != ' ' && filter.ptr[j] != '\t') ++j;
        if (j > i) {
            const str_t word = {filter.ptr + i, j - i};
            if (!contains_ignore_case(path, word) && !contains_ignore_case(label, word) && !contains_ignore_case(src, word)) {
                return false;
            }
        }
        i = j;
    }
    return true;
}

struct Row {
    SeriesSource source;
    uint32_t period;        // a timeline period the script is taken over, 0 for every frame
    const md_attributes_t* table;
    const md_attribute_t* attr;  // valid for the frame: nothing publishes while the window draws
    uint32_t views;         // SeriesView_*
    const char* reason;     // why no view can show it
};

struct Attributes : viamd::EventHandler {
    bool show_window = false;
    bool plottable_only = false;
    bool focus_filter = false;
    char filter[128] = "";

    Attributes() {
        viamd::event_system_register_handler(*this);
    }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize:
                workspace_register_window("Attributes", &show_window);
                break;
            case viamd::EventType_ViamdFrameTick: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                // Global, below whatever has focus: a text field keeps its own use of the keys
                if (ImGui::Shortcut(KEY_SHOW_ATTRIBUTE_WINDOW, ImGuiInputFlags_RouteGlobal)) {
                    show_window = !show_window;
                    focus_filter = show_window;
                }
                draw(state);
                break;
            }
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("Attributes", &show_window);
                ImGui::SameLine();
                ImGui::TextDisabled("Ctrl+Shift+A");
                break;
            default:
                break;
            }
        }
    }

    void draw(ApplicationState& state) {
        if (!show_window) return;

        ImGui::SetNextWindowSize(ImVec2(900, 500), ImGuiCond_FirstUseEver);
        if (!ImGui::Begin("Attributes", &show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            ImGui::End();
            return;
        }
        if (focus_filter) {
            ImGui::SetWindowFocus();
            ImGui::SetKeyboardFocusHere();
            focus_filter = false;
        }

        ImGui::SetNextItemWidth(-ImGui::CalcTextSize("Plottable only").x - ImGui::GetFrameHeight() - ImGui::GetStyle().ItemSpacing.x * 2 - ImGui::GetStyle().ItemInnerSpacing.x);
        ImGui::InputTextWithHint("##filter", "Filter: words in the path, label or source, all must match", filter, sizeof(filter));
        ImGui::SameLine();
        ImGui::Checkbox("Plottable only", &plottable_only);

        md_temp_scope_t temp = md_temp_begin();
        md_allocator_i* alloc = md_temp_allocator(temp);

        // The tables there are right now, and the script's over each timeline period, with only what it has of its
        // own over a period: the distributions and volumes, over the period's frames. Everything else in the table is
        // the same over every frame.
        struct Source {
            SeriesSource source;
            uint32_t period;
        };
        Source sources[2 + TIMELINE_MAX_PERIODS];
        int num_sources = 0;
        sources[num_sources++] = { SeriesSource_System, 0 };
        sources[num_sources++] = { SeriesSource_Script, 0 };
        for (int i = 0; i < state.timeline.periods.count; ++i) {
            sources[num_sources++] = { SeriesSource_Script, state.timeline.periods.items[i].id };
        }
        md_array(Row) rows = 0;
        size_t total = 0;
        const str_t filter_str = str_from_cstr(filter);
        for (int si = 0; si < num_sources; ++si) {
            const SeriesSource source = sources[si].source;
            const uint32_t period = sources[si].period;
            const md_attributes_t* table = series_table(&state, source);
            if (!table) continue;
            if (period == 0) total += md_attributes_count(table);
            char src_label[96];
            source_label(src_label, sizeof(src_label), state, source, period);
            for (md_attribute_iter_t it = md_attributes_iter(table, {}); md_attributes_next(&it);) {
                const md_attribute_t& attr = *it.attr;
                if (period != 0 && !script_aggregate_path(state, attr.path)) continue;
                if (period != 0) total += 1;
                if (!matches_filter(filter_str, attr.path, attr.label, src_label)) continue;
                Row row = { source, period, table, it.attr, 0, nullptr };
                const SeriesKey key = series_key(source, attr.path, SeriesVariant_Values, period);
                row.views = series_views(&state, key, &row.reason);
                if (plottable_only && !row.views) continue;
                md_array_push(rows, row, alloc);
            }
        }

        ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
        ImGui::TextWrapped("%zu of %zu attributes. Drag a row marked " ICON_FA_CHART_LINE " or " ICON_FA_CHART_COLUMN " into a Timelines or Distributions subplot, one the script reads into the Script Editor.", md_array_size(rows), total);
        ImGui::PopStyleColor();

        const ImGuiTableFlags table_flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_Resizable |
            ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_Hideable;
        if (ImGui::BeginTable("##attributes", 8, table_flags)) {
            const float icon_w = ImGui::GetFrameHeight();
            ImGui::TableSetupScrollFreeze(0, 1);
            ImGui::TableSetupColumn(ICON_FA_CHART_LINE,   ImGuiTableColumnFlags_WidthFixed | ImGuiTableColumnFlags_NoResize, icon_w);
            ImGui::TableSetupColumn(ICON_FA_CHART_COLUMN, ImGuiTableColumnFlags_WidthFixed | ImGuiTableColumnFlags_NoResize, icon_w);
            ImGui::TableSetupColumn("Path",    ImGuiTableColumnFlags_WidthStretch, 3.0f);
            ImGui::TableSetupColumn("Label",   ImGuiTableColumnFlags_WidthStretch, 1.6f);
            ImGui::TableSetupColumn("Format",  ImGuiTableColumnFlags_WidthStretch, 1.5f);
            ImGui::TableSetupColumn("Unit",    ImGuiTableColumnFlags_WidthStretch, 0.9f);
            ImGui::TableSetupColumn("Storage", ImGuiTableColumnFlags_WidthStretch, 0.9f);
            ImGui::TableSetupColumn("Source",  ImGuiTableColumnFlags_WidthStretch, 2.6f);

            // Header, with what the two marks mean
            ImGui::TableNextRow(ImGuiTableRowFlags_Headers);
            const char* header_tips[2] = { "Can be shown in the Timelines window", "Can be shown in the Distributions window" };
            for (int c = 0; c < 8; ++c) {
                if (!ImGui::TableSetColumnIndex(c)) continue;
                ImGui::TableHeader(ImGui::TableGetColumnName(c));
                if (c < 2 && ImGui::IsItemHovered()) ImGui::SetTooltip("%s", header_tips[c]);
            }

            const ImVec4 on_color  = ImGui::GetStyleColorVec4(ImGuiCol_Text);
            const ImVec4 off_color = ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled);

            ImGuiListClipper clipper;
            clipper.Begin((int)md_array_size(rows));
            while (clipper.Step()) {
                for (int r = clipper.DisplayStart; r < clipper.DisplayEnd; ++r) {
                    const Row& row = rows[r];
                    const md_attribute_t& attr = *row.attr;
                    const SeriesKey key = series_key(row.source, attr.path, SeriesVariant_Values, row.period);
                    const bool timeline = row.views & SeriesView_Timeline;
                    const bool distribution = row.views & SeriesView_Distribution;

                    ImGui::TableNextRow();
                    ImGui::PushID(r);

                    ImGui::TableSetColumnIndex(0);
                    ImGui::TextColored(timeline ? on_color : off_color, "%s", timeline ? ICON_FA_CHART_LINE : "-");
                    ImGui::TableSetColumnIndex(1);
                    ImGui::TextColored(distribution ? on_color : off_color, "%s", distribution ? ICON_FA_CHART_COLUMN : "-");

                    // The path spans the row: it is what is hovered, dragged and right clicked
                    ImGui::TableSetColumnIndex(2);
                    char path_buf[512];
                    str_copy_to_char_buf(path_buf, sizeof(path_buf), attr.path);
                    ImGui::Selectable(path_buf, false, ImGuiSelectableFlags_SpanAllColumns);
                    const bool hovered = ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort);

                    // How the script reads it, if it does: an attribute of the system the script is compiled against
                    char expr[sizeof(path_buf) + 16];
                    const bool scriptable = row.source == SeriesSource_System && script_editor::attr_expression(expr, sizeof(expr), state.script.language, row.table, row.attr);

                    if ((row.views || scriptable) && ImGui::BeginDragDropSource()) {
                        if (row.views) {
                            char label[96];
                            series_label(label, sizeof(label), &state, key);
                            // As a timeline series when it is one, which a distribution subplot accepts
                            // as well; otherwise as what only a distribution subplot takes. The Script
                            // Editor takes either as the expression that reads it.
                            const char* dnd_type = timeline ? TIMELINE_SERIES_DND : DISTRIBUTION_SERIES_DND;
                            series_set_drag_payload(dnd_type, key, -1, label, series_default_color(&state, key));
                        } else {
                            // What only the Script Editor takes
                            ImGui::SetDragDropPayload(SCRIPT_TEXT_DND, expr, strlen(expr) + 1);
                            ImGui::TextUnformatted(expr);
                        }
                        ImGui::EndDragDropSource();
                    }

                    if (ImGui::BeginPopupContextItem("##context")) {
                        if (ImGui::MenuItem("Show in timeline", nullptr, false, timeline)) {
                            plot_add_series(&state, state.timeline.subplots[0], key);
                            state.timeline.show_window = true;
                        }
                        if (ImGui::MenuItem("Show distribution", nullptr, false, distribution)) {
                            plot_add_series(&state, state.distributions.subplots[0], key);
                            state.distributions.show_window = true;
                        }
                        ImGui::Separator();
                        if (ImGui::MenuItem("Copy path")) {
                            char buf[512];
                            snprintf(buf, sizeof(buf), STR_FMT, STR_ARG(attr.path));
                            ImGui::SetClipboardText(buf);
                        }
                        if (scriptable && ImGui::MenuItem("Copy script expression", expr)) {
                            ImGui::SetClipboardText(expr);
                        }
                        ImGui::EndPopup();
                    }


                    ImGui::TableSetColumnIndex(3);
                    if (!str_empty(attr.label)) ImGui::TextUnformatted(attr.label.ptr, attr.label.ptr + attr.label.len);

                    ImGui::TableSetColumnIndex(4);
                    {
                        char buf[96];
                        format_string(buf, sizeof(buf), attr.format);
                        ImGui::TextUnformatted(buf);
                    }

                    ImGui::TableSetColumnIndex(5);
                    {
                        char buf[32] = "";
                        md_unit_print(buf, sizeof(buf), attr.unit);
                        ImGui::TextUnformatted(buf[0] ? buf : "-");
                    }

                    ImGui::TableSetColumnIndex(6);
                    ImGui::TextUnformatted(storage_name(attr));

                    ImGui::TableSetColumnIndex(7);
                    {
                        char buf[96];
                        source_label(buf, sizeof(buf), state, row.source, row.period);
                        ImGui::TextUnformatted(buf);
                        if (row.period != 0 && ImGui::IsItemHovered()) timeline_periods_highlight(row.period);
                    }

                    if (hovered) {
                        ImGui::BeginTooltip();
                        ImGui::TextUnformatted(attr.path.ptr, attr.path.ptr + attr.path.len);
                        if (row.period != 0) {
                            ImGui::TextDisabled("Over the frames of %s", series_period_name(&state, row.period));
                        }
                        if (!str_empty(attr.description)) {
                            ImGui::PushTextWrapPos(ImGui::GetFontSize() * 30.0f);
                            ImGui::TextDisabled(STR_FMT, STR_ARG(attr.description));
                            ImGui::PopTextWrapPos();
                        }
                        if (attr.format.type == MD_ATTRIBUTE_TYPE_STR && attr.format.rank <= 1 && md_attribute_value_count(&attr.format) == 1) {
                            const str_t value = md_attribute_str(row.table, &attr, 0);
                            ImGui::Text("Value: " STR_FMT, STR_ARG(value));
                        }
                        ImGui::Separator();
                        const char* plots = timeline && distribution ? "a Timelines or Distributions subplot" :
                                            timeline ? "a Timelines subplot" : distribution ? "a Distributions subplot" : nullptr;
                        if (plots && scriptable) {
                            ImGui::Text("Drag into %s, or into the Script Editor as %s", plots, expr);
                        } else if (plots) {
                            ImGui::Text("Drag into %s", plots);
                        } else if (scriptable) {
                            ImGui::Text("Drag into the Script Editor as %s", expr);
                        }
                        if (distribution && !timeline && row.reason) {
                            ImGui::TextDisabled("Not over time: %s", row.reason);
                        } else if (!plots && row.reason) {
                            ImGui::TextDisabled("Cannot be plotted: %s", row.reason);
                        }
                        ImGui::EndTooltip();
                    }

                    ImGui::PopID();
                }
            }
            ImGui::EndTable();
        }

        md_temp_end(temp);
        ImGui::End();
    }
};

static Attributes instance;

}  // namespace attributes
