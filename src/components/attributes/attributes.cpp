// Attributes window: every attribute the application holds, in one table.
//
// The system's table (topology columns, the run and everything loaded along it) and the script
// evaluations' tables (script/<ident> and its siblings) are listed side by side, one row per
// attribute, with what it is, how it is stored and whether the Timelines and Distributions windows
// can show it. A row that can be shown is dragged from here into a subplot of either window.
//
// Opened with KEY_SHOW_ATTRIBUTE_WINDOW (Ctrl+Shift+A) from anywhere, or from the Windows menu.

#include <viamd_event.h>
#include <viamd.h>
#include <event.h>
#include <plot_series.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_str.h>
#include <md_system.h>

#include <imgui.h>
#include <imgui_widgets.h>
#include <implot.h>
#include <app/IconsFontAwesome6.h>

#include <stdio.h>

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

static const char* storage_name(md_attribute_storage_t storage) {
    switch (storage) {
    case MD_ATTRIBUTE_STORAGE_RESIDENT: return "resident";
    case MD_ATTRIBUTE_STORAGE_VIRTUAL:  return "virtual";
    case MD_ATTRIBUTE_STORAGE_ALIAS:    return "alias";
    default: return "-";
    }
}

static const char* source_name(SeriesSource source) {
    switch (source) {
    case SeriesSource_System:         return "System";
    case SeriesSource_Script:         return "Script";
    case SeriesSource_ScriptFiltered: return "Script (filtered)";
    default: return "?";
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
    const md_attributes_t* table;
    uint32_t index;         // into table->attr
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

        // The tables there are right now. The filtered evaluation only while the timeline has a filter,
        // otherwise it repeats the full one row for row.
        const SeriesSource sources[] = { SeriesSource_System, SeriesSource_Script, SeriesSource_ScriptFiltered };
        md_array(Row) rows = 0;
        size_t total = 0;
        const str_t filter_str = str_from_cstr(filter);
        for (SeriesSource source : sources) {
            if (source == SeriesSource_ScriptFiltered && !state.timeline.filter.enabled) continue;
            const md_attributes_t* table = series_table(&state, source);
            if (!table) continue;
            const size_t count = md_array_size(table->attr);
            total += count;
            for (size_t i = 0; i < count; ++i) {
                const md_attribute_t& attr = table->attr[i];
                if (!matches_filter(filter_str, attr.path, attr.label, source_name(source))) continue;
                Row row = { source, table, (uint32_t)i, 0, nullptr };
                const SeriesKey key = series_key(source, attr.path);
                row.views = series_views(&state, key, &row.reason);
                if (plottable_only && !row.views) continue;
                md_array_push(rows, row, alloc);
            }
        }

        ImGui::TextDisabled("%zu of %zu attributes. Drag a row marked " ICON_FA_CHART_LINE " or " ICON_FA_CHART_COLUMN " into a Timelines or Distributions subplot.", md_array_size(rows), total);

        const ImGuiTableFlags table_flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_Resizable |
            ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_Hideable;
        if (ImGui::BeginTable("##attributes", 8, table_flags)) {
            const float icon_w = ImGui::GetFrameHeight();
            ImGui::TableSetupScrollFreeze(0, 1);
            ImGui::TableSetupColumn(ICON_FA_CHART_LINE,   ImGuiTableColumnFlags_WidthFixed | ImGuiTableColumnFlags_NoResize, icon_w);
            ImGui::TableSetupColumn(ICON_FA_CHART_COLUMN, ImGuiTableColumnFlags_WidthFixed | ImGuiTableColumnFlags_NoResize, icon_w);
            ImGui::TableSetupColumn("Path",    ImGuiTableColumnFlags_WidthStretch, 4.0f);
            ImGui::TableSetupColumn("Label",   ImGuiTableColumnFlags_WidthStretch, 2.0f);
            ImGui::TableSetupColumn("Format",  ImGuiTableColumnFlags_WidthStretch, 1.6f);
            ImGui::TableSetupColumn("Unit",    ImGuiTableColumnFlags_WidthStretch, 1.0f);
            ImGui::TableSetupColumn("Storage", ImGuiTableColumnFlags_WidthStretch, 0.9f);
            ImGui::TableSetupColumn("Source",  ImGuiTableColumnFlags_WidthStretch, 1.0f);

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
                    const md_attribute_t& attr = row.table->attr[row.index];
                    const SeriesKey key = series_key(row.source, attr.path);
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

                    if (row.views && ImGui::BeginDragDropSource()) {
                        char label[96];
                        series_label(label, sizeof(label), &state, key);
                        // As a timeline series when it is one, which a distribution subplot accepts
                        // as well; otherwise as what only a distribution subplot takes
                        const char* dnd_type = timeline ? TIMELINE_SERIES_DND : DISTRIBUTION_SERIES_DND;
                        series_set_drag_payload(dnd_type, key, -1, label, series_default_color(&state, key));
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
                    ImGui::TextUnformatted(storage_name(attr.storage));

                    ImGui::TableSetColumnIndex(7);
                    ImGui::TextUnformatted(source_name(row.source));

                    if (hovered) {
                        ImGui::BeginTooltip();
                        ImGui::TextUnformatted(attr.path.ptr, attr.path.ptr + attr.path.len);
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
                        if (timeline && distribution) {
                            ImGui::TextUnformatted("Drag into a Timelines or Distributions subplot");
                        } else if (timeline) {
                            ImGui::TextUnformatted("Drag into a Timelines subplot");
                        } else if (distribution) {
                            ImGui::TextUnformatted("Drag into a Distributions subplot");
                            if (row.reason) ImGui::TextDisabled("Not over time: %s", row.reason);
                        } else if (row.reason) {
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
