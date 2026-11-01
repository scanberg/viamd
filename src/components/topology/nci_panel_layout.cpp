#include "nci_panel_layout.h"

#include "critical_point_distribution_core.h"

#include <color_scale.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <layer_groups.h>
#include <md_util.h>
#include <property_grid.h>
#include <representation_tree.h>

#include <imgui_internal.h>
#include <implot_widgets.h>
#include <app/IconsFontAwesome6.h>

#include <algorithm>
#include <cfloat>
#include <string>
#include <array>
#include <cstdio>
#include <cstring>
#include <memory>
#include <vector>

namespace topology {

struct NciPanelDrawCallbacks {
    void* user_data = nullptr;
    void (*draw_scatter)(void*) = nullptr;
    void (*draw_layers)(void*) = nullptr;
    void (*on_enabled)(void*) = nullptr;
};

static void draw_nci_panel_layout(bool* enabled, NciPanelDrawCallbacks callbacks) {
    static bool show_layers = true;
    static bool restore_split = false;
    static float scatter_height = 0.0f;
    if (!enabled) return;

    if (ImGui::Checkbox("NCI enabled", enabled) && callbacks.on_enabled)
        callbacks.on_enabled(callbacks.user_data);
    ImGui::SetItemTooltip("Enable NCI computation and display. Disabling retains saved objects and their settings.");
    ImGui::SameLine();
    if (ImGui::Button(show_layers ? ICON_FA_LAYER_GROUP "##nci_layers" : ICON_FA_LAYER_GROUP "##nci_layers_off")) {
        show_layers = !show_layers;
        restore_split = show_layers;
    }
    ImGui::SetItemTooltip(show_layers ? "Hide the layers" : "Show the layers");

    ImGui::BeginDisabled(!*enabled);
    const float available_height = std::max(1.0f, ImGui::GetContentRegionAvail().y);
    if (scatter_height <= 0.0f) scatter_height = available_height * 0.6f;
    const float clamped = std::clamp(scatter_height, std::min(120.0f, available_height),
                                     std::max(std::min(120.0f, available_height), available_height - 120.0f));
    // A resizable child keeps its own size: while the layers were hidden it took the full height, so the split is
    // set back explicitly, as it is when the window shrinks below it
    if (show_layers && (restore_split || clamped != scatter_height))
        ImGui::SetNextWindowSize(ImVec2(0.0f, clamped), ImGuiCond_Always);
    scatter_height = clamped;
    restore_split = false;

    ImGui::BeginChild("##NCI scatter panel", ImVec2(0, show_layers ? scatter_height : available_height),
                      ImGuiChildFlags_Borders | (show_layers ? ImGuiChildFlags_ResizeY : ImGuiChildFlags_None),
                      ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
    if (callbacks.draw_scatter) callbacks.draw_scatter(callbacks.user_data);
    ImGui::EndChild();
    if (show_layers) {
        scatter_height = ImGui::GetItemRectSize().y;
        ImGui::BeginChild("##NCI layer panel", ImVec2(0, 0), ImGuiChildFlags_Borders);
        if (callbacks.draw_layers) callbacks.draw_layers(callbacks.user_data);
        ImGui::EndChild();
    }
    ImGui::EndDisabled();
}

namespace {

bool surface_ready(const NciSurface& n) {
    return n.want_key && n.have_key == n.want_key && !n.cancelled && !n.failed;
}

// In the scatter plot, a layer is its chip colour and an interaction always the same colour
ImVec4 type_colour(const NciSurface& n) {
    return nci_type_colour(n.group);
}

ImVec4 layer_chip(int key, const std::vector<NciSurface*>& members) {
    if (key == NCI_LAYER_FULL_SYSTEM) return nci_layer_colour(0);
    if (key == NCI_LAYER_PREVIEW) return nci_layer_colour(1);
    return members.empty() ? ImVec4(0.5f, 0.5f, 0.5f, 1.0f) : members.front()->chip;
}

const char* layer_name(int key, const std::vector<NciSurface*>& members) {
    if (key == NCI_LAYER_FULL_SYSTEM) return "Full system";
    if (key == NCI_LAYER_PREVIEW) return "Preview";
    return members.empty() ? "Layer" : members.front()->name;
}

struct NciPanel {
    static constexpr uint32_t NCI_LAYER_LIST_ID = 3;
    static constexpr const char* resolution_str[] = { "Coarse", "Medium", "Fine" };
    static constexpr const char* mode_str[] = { "Uniform", "CP-adaptive" };
    static constexpr const char* mode_tooltip =
        "How the isovalue and density cutoff apply to every NCI layer, the full system included.\n"
        "Uniform: one isovalue and cutoff for every CP type; the physically consistent mode.\n"
        "CP-adaptive: around each CP, that CP type's own isovalue and cutoff. Illustrative only: it\n"
        "brings out weak interactions such as pi-stacking, but the surfaces are not mutually consistent.\n"
        "In both modes the extent around the CPs is the radius of each CP type, and one colour scale applies.";
    static constexpr size_t max_samples = 32u << 20;

    NciState& st;
    bool critical_points_enabled;

    NciPanel(NciState& state, bool cp_enabled) : st(state), critical_points_enabled(cp_enabled) {}

    // The member of the shown layer the scatter guides follow
    NciSurface* active_member() {
        const auto members = nci_layer_members(st, st.selected_layer);
        for (NciSurface* n : members) if (n->id == st.active_id) return n;
        return members.empty() ? nullptr : members.front();
    }

    // The lowest "Layer N" not taken by a saved layer; the group id stays unique and is not reused
    int free_layer_number() const {
        for (int number = 1;; ++number) {
            char name[64];
            snprintf(name, sizeof(name), "Layer %d", number);
            if (std::none_of(st.objects.begin(), st.objects.end(),
                             [&](const auto& n) { return strcmp(n->name, name) == 0; }))
                return number;
        }
    }

    std::unique_ptr<NciSurface> copy_preview(const NciSurface& source, int group_id, int number) {
        auto n = std::make_unique<NciSurface>();
        n->saved = true;
        n->id = st.next_id++;
        n->object_group = group_id;
        snprintf(n->name, sizeof(n->name), "Layer %d", number);
        n->shown = true;
        n->layer_visible = true;
        n->s_iso = source.s_iso;
        n->rho_max = source.rho_max;
        n->scale = source.scale;
        n->families = source.families;
        n->type_counts = source.type_counts;
        n->group = source.group;
        n->environments = source.environments;
        n->opacity = source.opacity;
        n->resolution = source.resolution;
        n->scope = source.scope;
        n->near_radius = source.near_radius;
        n->connected = source.connected;
        n->radius_cap = source.radius_cap;
        n->anchors = source.anchors;
        n->source_mesh = source.source_mesh;
        n->scatter = source.scatter;
        n->outline = source.outline;
        n->rho = source.rho;
        n->grid = source.grid;
        n->rho_key = source.rho_key;
        n->want_key = source.want_key;
        n->have_key = source.have_key;
        n->generation = source.generation;
        n->seconds = source.seconds;
        n->rho_seconds = source.rho_seconds;
        if (!topology::upload_nci_surface_mesh(*n)) { topology::nci_surface_release(*n); return nullptr; }
        return n;
    }

    bool preview_ready() const {
        bool any_surface = false, ready = true;
        for (const auto& n : st.previews) {
            if (!n->shown) continue;
            any_surface |= n->triangles > 0 && n->mesh.index_count > 0;
            ready &= nci_scale_member(*n) && surface_ready(*n) && !n->job && !n->editing;
        }
        return any_surface && ready;
    }

    void save_preview() {
        if (!preview_ready()) {
            MD_LOG_ERROR("NCI surface: wait for the preview to complete before saving it as a layer");
            return;
        }
        std::vector<std::unique_ptr<NciSurface>> copies;
        const int group_id = st.next_group_id++;
        const int number = free_layer_number();
        for (const auto& n : st.previews) {
            if (!n->shown) continue;
            auto saved = copy_preview(*n, group_id, number);
            if (!saved) {
                for (auto& c : copies) topology::nci_surface_release(*c);
                return;
            }
            copies.push_back(std::move(saved));
        }
        const ImVec4 chip = nci_layer_colour(1 + group_id);
        for (auto& n : copies) {
            n->chip = chip;
            st.objects.push_back(std::move(n));
        }
        st.layer_order.push_back(group_id);
        nci_select_layer(st, group_id);
        nci_update_scale(st);
    }

    void delete_layer(int key) {
        for (size_t i = 0; i < st.objects.size();) {
            if (nci_layer_key(*st.objects[i]) == key) {
                topology::nci_surface_release(*st.objects[i]);
                st.objects.erase(st.objects.begin() + i);
            } else ++i;
        }
        st.layer_order.erase(std::remove(st.layer_order.begin(), st.layer_order.end(), key), st.layer_order.end());
        st.scatter_layers.erase(key);
        if (st.selected_layer == key) nci_select_layer(st, NCI_LAYER_FULL_SYSTEM);
        nci_update_scale(st);
    }

    static void draw_environments(const std::vector<std::vector<std::string>>& environments) {
        if (environments.empty() || !ImGui::TreeNode("Chemical environments")) return;
        const auto tree = critical_point_distribution::environment_tree(environments);
        const auto draw = [&](auto&& self, uint32_t idx) -> void {
            const auto& node = tree[idx];
            ImGui::PushID((int)idx);
            if (node.children.empty()) ImGui::BulletText("%s", node.name.c_str());
            else if (ImGui::TreeNode(node.name.c_str())) {
                for (uint32_t child : node.children) self(self, child);
                ImGui::TreePop();
            }
            ImGui::PopID();
        };
        for (uint32_t child : tree[0].children) draw(draw, child);
        ImGui::TreePop();
    }

    // ---- settings: the general mode, the shared options, the colour scale and the per-type values, in one place
    void draw_settings() {
        bool editing = false;
        const bool adaptive = st.mode == NciParameterMode::CPAdaptive;
        prop::begin("NCI settings");
        prop::row("mode", mode_tooltip);
        int mode = adaptive ? 1 : 0;
        ImGui::SetNextItemWidth(prop::control_width());
        if (ImGui::Combo("##NCI mode", &mode, mode_str, (int)ARRAY_SIZE(mode_str)))
            st.mode = mode ? NciParameterMode::CPAdaptive : NciParameterMode::Uniform;
        if (!adaptive) {
            prop::row("isovalue s", "Reduced density gradient isovalue, the same for every CP type.");
            ImGui::SetNextItemWidth(prop::control_width());
            ImGui::SliderFloat("##NCI uniform s", &st.uniform.s_iso, 0.1f, 1.5f, "%.2f");
            editing |= ImGui::IsItemActive();
            prop::row("\xCF\x81 cutoff", "Density cutoff, the same for every CP type.");
            ImGui::SetNextItemWidth(prop::control_width());
            ImGui::SliderFloat("##NCI uniform rho", &st.uniform.rho_max, 0.005f, 0.2f, "< %.3f", ImGuiSliderFlags_Logarithmic);
            editing |= ImGui::IsItemActive();
        }
        prop::row("resolution", "Grid spacing for density sampling. Finer grids are smoother and take longer.");
        ImGui::SetNextItemWidth(prop::control_width());
        ImGui::Combo("##NCI resolution", &st.resolution, resolution_str, (int)ARRAY_SIZE(resolution_str));
        ImGui::SetItemTooltip("Coarse 0.3, Medium 0.2, Fine 0.13 Bohr; made coarser where that would take more than %zu M samples.",
                              max_samples >> 20);
        prop::row("connected regions", "Keep six-connected low-s regions seeded at the anchor CPs. Association depends on "
                                       "grid resolution and thresholds, not a unique physical partition.");
        ImGui::Checkbox("##NCI connected", &st.connected);
        if (!st.connected) st.radius_cap = true;
        ImGui::SameLine();
        ImGui::BeginDisabled(!st.connected);
        ImGui::Checkbox("Radius cap", &st.radius_cap);
        ImGui::EndDisabled();
        ImGui::SetItemTooltip("Limit the extent around the CPs to each type's radius. Without a cap, connected regions may merge across the system.");
        prop::row("opacity");
        ImGui::SetNextItemWidth(prop::control_width());
        ImGui::SliderFloat("##NCI opacity", &st.opacity, 0.05f, 1.0f, "%.2f");
        color_scale_draw_controls(&st.scale, st.span, md_unit_pow(md_unit_bohr_radius(), -3), "shown layer",
                                  "the 1st to 99th percentile of the shown layer", true);
        prop::end();

        ImGui::SeparatorText("Per CP type");
        const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_SizingFixedFit;
        if (ImGui::BeginTable("##NCI type settings", adaptive ? 4 : 2, flags)) {
            const float em = ImGui::GetFontSize();
            ImGui::TableSetupColumn("Type");
            ImGui::TableSetupColumn("Radius", ImGuiTableColumnFlags_WidthFixed, 6.0f * em);
            if (adaptive) {
                ImGui::TableSetupColumn("s", ImGuiTableColumnFlags_WidthFixed, 5.0f * em);
                ImGui::TableSetupColumn("\xCF\x81 cutoff", ImGuiTableColumnFlags_WidthFixed, 6.0f * em);
            }
            ImGui::TableHeadersRow();
            for (size_t f = 0; f < st.family_settings.size(); ++f) {
                if (!nci_surface::has_nci((nci_surface::Family)f)) continue;
                NciFamilySettings& fs = st.family_settings[f];
                ImGui::PushID((int)f);
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::AlignTextToFramePadding();
                ImGui::TextUnformatted(nci_surface::family_name((nci_surface::Family)f));
                if (!nci_surface::family_defaults(nci_surface::family_bit((nci_surface::Family)f)).established)
                    ImGui::SetItemTooltip("Uncalibrated type: provisional defaults");
                ImGui::TableNextColumn();
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::SliderFloat("##radius", &fs.radius, 0.5f, 12.0f, "%.1f Bohr");
                editing |= ImGui::IsItemActive();
                if (adaptive) {
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    ImGui::SliderFloat("##s", &fs.s_iso, 0.1f, 1.5f, "%.2f");
                    editing |= ImGui::IsItemActive();
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    ImGui::SliderFloat("##rho", &fs.rho_max, 0.005f, 0.2f, "< %.3f", ImGuiSliderFlags_Logarithmic);
                    editing |= ImGui::IsItemActive();
                }
                ImGui::PopID();
            }
            ImGui::EndTable();
        }
        if (ImGui::Button("Reset to defaults")) {
            for (size_t f = 0; f < st.family_settings.size(); ++f)
                st.family_settings[f] = nci_family_default_settings((nci_surface::Family)f);
            st.uniform = nci_uniform_default_settings();
            const bool legend = st.scale.show_legend;
            st.scale = nci_default_scale();
            st.scale.show_legend = legend;
            nci_update_scale(st);
        }
        st.settings_editing |= editing;
    }

    // ---- the layer list: Full system and Preview, a line, then the saved layers
    struct LayerRow {
        int key;
        const char* label;
        std::vector<NciSurface*> members;
    };

    bool draw_row(LayerRow& row, bool saved, RepresentationNode* const* nodes, size_t node_count,
                  LayerNodeRef* move_src, RepresentationDrop* move_drop) {
        bool dropped = false;
        ImGui::PushID(row.key);
        const float line_h = ImGui::GetTextLineHeight();
        const float row_h = line_h * 2.0f + 8.0f;
        const ImVec2 row_min = ImGui::GetCursorScreenPos();
        const float row_w = ImGui::GetContentRegionAvail().x;
        const ImRect rect(row_min, ImVec2(row_min.x + row_w, row_min.y + row_h));
        const bool selected = st.selected_layer == row.key;
        const bool visible = selected && !st.layer_hidden;
        ImGui::SetNextItemAllowOverlap();
        // The row: the layer in the viewport, one at a time; clicking the one shown hides it
        if (ImGui::Selectable("##row", visible, ImGuiSelectableFlags_AllowOverlap, ImVec2(row_w, row_h))) {
            if (visible) {
                st.layer_hidden = true;
                nci_apply_settings(st);
            } else nci_select_layer(st, row.key);
        }
        const bool row_hovered = ImGui::IsItemHovered();
        if (saved) {
            layer_tree_drag_source({NCI_LAYER_LIST_ID, (uint32_t)row.key, false}, row.label);
            if (layer_tree_drop_target(rect, row_min.x, {NCI_LAYER_LIST_ID, (uint32_t)row.key, false},
                                       false, nodes, node_count, move_src, move_drop))
                dropped = true;
            if (ImGui::BeginPopupContextItem("##layer menu")) {
                if (ImGui::MenuItem("Delete layer")) st.pending_delete = row.key;
                ImGui::EndPopup();
            }
        }

        ImDrawList* dl = ImGui::GetWindowDrawList();
        const float pad = 4.0f;
        const float button_w = ImGui::GetFrameHeight();
        const float text_x = row_min.x + pad + button_w + 6.0f + line_h + 8.0f;
        const ImU32 text_col = ImGui::GetColorU32(visible ? ImGuiCol_Text : ImGuiCol_TextDisabled);
        const ImU32 dim_col = ImGui::GetColorU32(ImGuiCol_TextDisabled);
        const ImVec4 clip(text_x, row_min.y, row_min.x + row_w - pad, row_min.y + row_h);

        ImGui::SetCursorScreenPos(ImVec2(row_min.x + pad, row_min.y + (row_h - button_w) * 0.5f));
        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(2, ImGui::GetStyle().FramePadding.y));
        // The eye: the layer in the scatter plot, any number of them
        const bool in_scatter = scatter_visible(row.key);
        if (ImGui::Button(in_scatter ? ICON_FA_EYE "##visible" : ICON_FA_EYE_SLASH "##visible"))
            set_scatter_visible(row.key, !in_scatter);
        ImGui::PopStyleVar();
        ImGui::PopStyleColor();
        ImGui::SetItemTooltip(in_scatter ? "Hide in the scatter plot" : "Show in the scatter plot");

        ImGui::SetCursorScreenPos(ImVec2(row_min.x + pad + button_w + 6.0f, row_min.y + (row_h - line_h) * 0.5f));
        if (saved && !row.members.empty()) {
            ImVec4 chip = row.members.front()->chip;
            if (ImGui::ColorEdit3("##chip", &chip.x, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel))
                for (NciSurface* n : row.members) n->chip = chip;
        } else {
            const ImVec4 chip = layer_chip(row.key, row.members);
            ImGui::ColorButton("##chip", chip, ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoDragDrop, ImVec2(line_h, line_h));
        }

        const std::string types = nci_layer_types(row.members);
        size_t triangles = 0;
        float progress = -1.0f;
        for (const NciSurface* n : row.members) {
            triangles += n->triangles;
            if (n->job) progress = MAX(progress, n->job->progress.load(std::memory_order_relaxed));
        }
        // As the Ramachandran layers: where the samples fall, here in the interactions' expected regions
        const RegionStats stats = region_stats(row.members);
        char expected[64] = "";
        if (stats.inside + stats.outside)
            snprintf(expected, sizeof(expected), " \xC2\xB7 %.0f %% expected", 100.0 * (double)stats.inside / (double)(stats.inside + stats.outside));
        char detail[256];
        if (progress >= 0.0f) snprintf(detail, sizeof(detail), "%s \xC2\xB7 computing %.0f %%", types.c_str(), 100.0f * progress);
        else if (visible) snprintf(detail, sizeof(detail), "%s \xC2\xB7 %zu triangles%s", types.c_str(), triangles, expected);
        else snprintf(detail, sizeof(detail), "%s%s", types.c_str(), expected);
        // CP-adaptive is not mutually consistent across parameter sets: flag the layers that mix them. Each member
        // has one set; CPs that must share it (a ring closed by a contact and that contact) are in one member.
        size_t parameter_sets = 0;
        for (const NciSurface* n : row.members) parameter_sets += n->families != 0;
        const bool illustrative = st.mode == NciParameterMode::CPAdaptive && parameter_sets > 1;
        const ImU32 warn_col = ImGui::GetColorU32(ImVec4(0.95f, 0.70f, 0.25f, visible ? 1.0f : 0.6f));
        const float icon_w = ImGui::CalcTextSize(ICON_FA_TRIANGLE_EXCLAMATION).x;
        const ImRect icon_rect(ImVec2(row_min.x + row_w - pad - icon_w, row_min.y + (row_h - line_h) * 0.5f),
                               ImVec2(row_min.x + row_w - pad, row_min.y + (row_h + line_h) * 0.5f));
        ImVec4 text_clip = clip;
        if (illustrative) {
            text_clip.z = icon_rect.Min.x - 6.0f;
            dl->AddText(icon_rect.Min, warn_col, ICON_FA_TRIANGLE_EXCLAMATION);
        }
        const float small = ImGui::GetFontSize() * 0.85f;
        const ImVec2 detail_pos(text_x, row_min.y + 4.0f + line_h + 1.0f);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize(), ImVec2(text_x, row_min.y + 4.0f), text_col, row.label, nullptr, 0.0f, &text_clip);
        dl->AddText(ImGui::GetFont(), small, detail_pos, dim_col, detail, nullptr, 0.0f, &text_clip);
        ImRect tag_rect;
        if (illustrative) {
            const char* tag = " \xC2\xB7 illustrative";
            const float detail_w = ImGui::GetFont()->CalcTextSizeA(small, FLT_MAX, 0.0f, detail).x;
            const ImVec2 tag_pos(detail_pos.x + detail_w, detail_pos.y);
            dl->AddText(ImGui::GetFont(), small, tag_pos, warn_col, tag, nullptr, 0.0f, &text_clip);
            tag_rect = ImRect(tag_pos, ImVec2(MIN(tag_pos.x + ImGui::GetFont()->CalcTextSizeA(small, FLT_MAX, 0.0f, tag).x, text_clip.z),
                                              tag_pos.y + small));
        }
        if (row_hovered && !ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
            if (illustrative && (icon_rect.Contains(ImGui::GetMousePos()) || tag_rect.Contains(ImGui::GetMousePos())))
                ImGui::SetTooltip("CP-adaptive: each CP type uses its own isovalue and \xCF\x81 cutoff.\n"
                                  "The surfaces of different types are not mutually consistent;\n"
                                  "use Uniform for a quantitative comparison.");
            else {
                ImGui::BeginTooltip();
                ImGui::TextUnformatted(row.label);
                ImGui::TextDisabled("%s", detail);
                draw_region_table(stats);
                ImGui::TextDisabled(visible ? "Click to hide it in the viewport"
                                            : "Click to show it in the viewport (one layer at a time)");
                ImGui::EndTooltip();
            }
        }
        ImGui::SetCursorScreenPos(ImVec2(row_min.x, row_min.y + row_h + 2.0f));
        ImGui::Dummy(ImVec2(0, 0));
        ImGui::PopID();
        return dropped;
    }

    // The interactions' expected regions in the scatter plot, pinned above the layers as the Ramachandran reference;
    // it has no surface. The Full system below is the CP reference.
    void draw_reference_row() {
        ImGui::PushID("##reference");
        const float line_h = ImGui::GetTextLineHeight();
        const float row_h = line_h * 2.0f + 8.0f;
        const ImVec2 row_min = ImGui::GetCursorScreenPos();
        const float row_w = ImGui::GetContentRegionAvail().x;
        const bool on = st.show_reference;
        ImGui::SetNextItemAllowOverlap();
        ImGui::Selectable("##row", false, ImGuiSelectableFlags_AllowOverlap, ImVec2(row_w, row_h));
        const bool row_hovered = ImGui::IsItemHovered();
        ImDrawList* dl = ImGui::GetWindowDrawList();
        const float pad = 4.0f;
        const float button_w = ImGui::GetFrameHeight();
        const float text_x = row_min.x + pad + button_w + 6.0f + line_h + 8.0f;
        const ImVec4 clip(text_x, row_min.y, row_min.x + row_w - pad, row_min.y + row_h);
        ImGui::SetCursorScreenPos(ImVec2(row_min.x + pad, row_min.y + (row_h - button_w) * 0.5f));
        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(2, ImGui::GetStyle().FramePadding.y));
        if (ImGui::Button(on ? ICON_FA_EYE "##visible" : ICON_FA_EYE_SLASH "##visible")) st.show_reference = !on;
        ImGui::PopStyleVar();
        ImGui::PopStyleColor();
        ImGui::SetItemTooltip(on ? "Hide the expected regions in the scatter plot" : "Show the expected regions in the scatter plot");
        const ImVec2 box(row_min.x + pad + button_w + 6.0f, row_min.y + (row_h - line_h) * 0.5f);
        dl->AddRectFilled(box, ImVec2(box.x + line_h, box.y + line_h), ImGui::GetColorU32(ImVec4(0.85f, 0.85f, 0.85f, 0.25f)));
        dl->AddRect(box, ImVec2(box.x + line_h, box.y + line_h), ImGui::GetColorU32(ImVec4(0.85f, 0.85f, 0.85f, 1.0f)));
        const ImU32 text_col = ImGui::GetColorU32(on ? ImGuiCol_Text : ImGuiCol_TextDisabled);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize(), ImVec2(text_x, row_min.y + 4.0f), text_col, "Expected regions", nullptr, 0.0f, &clip);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize() * 0.85f, ImVec2(text_x, row_min.y + 5.0f + line_h),
                    ImGui::GetColorU32(ImGuiCol_TextDisabled), "provisional, from the calibration systems", nullptr, 0.0f, &clip);
        if (row_hovered && !ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
            ImGui::BeginTooltip();
            ImGui::TextUnformatted("Expected regions");
            ImGui::TextDisabled("Where each interaction's spike is expected in the scatter plot. Provisional boxes\n"
                                "from the calibration systems, to be replaced by a reference map; densities from\n"
                                "other wave functions may shift them.");
            if (ImGui::BeginTable("##regions", 3, ImGuiTableFlags_SizingFixedFit)) {
                ImGui::TableSetupColumn("");
                ImGui::TableSetupColumn("sign(\xCE\xBB\xE2\x82\x82)\xCF\x81");
                ImGui::TableSetupColumn("s below");
                ImGui::TableHeadersRow();
                for (size_t f = 0; f < nci_surface::family_count; ++f) {
                    const nci_surface::ExpectedRegion r = nci_surface::expected_region((nci_surface::Family)f);
                    if (!r.valid) continue;
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::TextColored(nci_type_colour((nci_surface::Family)f), "%s", nci_surface::family_name((nci_surface::Family)f));
                    ImGui::TableNextColumn();
                    ImGui::Text("%+.3f to %+.3f", r.rho_min, r.rho_max);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.2f", r.s_max);
                }
                ImGui::EndTable();
            }
            ImGui::EndTooltip();
        }
        ImGui::SetCursorScreenPos(ImVec2(row_min.x, row_min.y + row_h + 2.0f));
        ImGui::Dummy(ImVec2(0, 0));
        ImGui::PopID();
    }

    // Only what needs the user: a running computation to cancel, or a failure to retry
    void draw_member_status(NciSurface& n) {
        const bool problem = !n.job && (n.cancelled || n.failed || n.rho_failed);
        if (!n.job && !problem) return;
        ImGui::PushID(n.id);
        ImGui::AlignTextToFramePadding();
        if (n.job) {
            ImGui::TextDisabled("%s: computing %.0f %%", nci_family_label(n).c_str(), 100.0f * n.job->progress.load(std::memory_order_relaxed));
            ImGui::SameLine();
            ImGui::BeginDisabled(n.job->cancel.load(std::memory_order_relaxed) != 0);
            if (ImGui::SmallButton("Cancel")) n.job->cancel.store(1, std::memory_order_relaxed);
            ImGui::EndDisabled();
        } else {
            ImGui::TextDisabled("%s: %s", nci_family_label(n).c_str(),
                                n.rho_failed ? "density evaluation failed" : n.failed ? "failed, see the log" : "cancelled");
            ImGui::SameLine();
            if (ImGui::SmallButton("Retry")) {
                n.have_key = 0;
                if (n.rho_failed) n.rho_key = 0;
                n.cancelled = n.failed = n.rho_failed = false;
            }
        }
        ImGui::PopID();
    }

    void draw_layers() {
        ImGui::AlignTextToFramePadding();
        ImGui::TextDisabled("LAYERS");
        ImGui::SameLine();
        int mode = st.mode == NciParameterMode::CPAdaptive ? 1 : 0;
        ImGui::SetNextItemWidth(9.0f * ImGui::GetFontSize());
        if (ImGui::Combo("##NCI layer mode", &mode, mode_str, (int)ARRAY_SIZE(mode_str)))
            st.mode = mode ? NciParameterMode::CPAdaptive : NciParameterMode::Uniform;
        ImGui::SetItemTooltip("%s", mode_tooltip);
        const float button = ImGui::GetFrameHeight();
        const float spacing = ImGui::GetStyle().ItemSpacing.x;
        ImGui::SameLine();
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + MAX(0.0f, ImGui::GetContentRegionAvail().x - 2.0f * button - spacing));
        ImGui::BeginDisabled(!preview_ready());
        if (ImGui::Button(ICON_FA_PLUS "##save", ImVec2(button, 0))) save_preview();
        ImGui::EndDisabled();
        ImGui::SetItemTooltip("Save the preview (the selected CPs) as a layer.\nSelect CPs and wait for the preview to complete.");
        ImGui::SameLine();
        if (ImGui::Button(ICON_FA_GEAR "##settings", ImVec2(button, 0))) ImGui::OpenPopup("##NCI settings");
        ImGui::SetItemTooltip("NCI settings: mode, isovalue, colour scale and the values of each CP type");
        if (ImGui::BeginPopup("##NCI settings")) {
            draw_settings();
            ImGui::EndPopup();
        }
        ImGui::Separator();

        draw_reference_row();
        LayerRow full = {NCI_LAYER_FULL_SYSTEM, "Full system", nci_layer_members(st, NCI_LAYER_FULL_SYSTEM)};
        draw_row(full, false, nullptr, 0, nullptr, nullptr);
        LayerRow preview = {NCI_LAYER_PREVIEW, "Preview", nci_layer_members(st, NCI_LAYER_PREVIEW)};
        if (!preview.members.empty()) draw_row(preview, false, nullptr, 0, nullptr, nullptr);

        // Saved layers, below a line, in the user's order
        std::vector<LayerRow> rows;
        for (const auto& n : st.objects) {
            const int key = nci_layer_key(*n);
            if (std::none_of(rows.begin(), rows.end(), [&](const LayerRow& r) { return r.key == key; }))
                rows.push_back({key, n->name, nci_layer_members(st, key)});
        }
        std::vector<LayerRow> ordered;
        std::vector<int> order;
        for (int key : st.layer_order) {
            auto it = std::find_if(rows.begin(), rows.end(), [&](const LayerRow& r) { return r.key == key; });
            if (it == rows.end()) continue;
            order.push_back(key);
            ordered.push_back(std::move(*it));
            rows.erase(it);
        }
        for (LayerRow& row : rows) {
            order.push_back(row.key);
            ordered.push_back(std::move(row));
        }
        st.layer_order = order;

        ImGui::SeparatorText("Saved layers");
        if (ordered.empty()) ImGui::TextDisabled("Select CPs, then " ICON_FA_PLUS " to save the preview as a layer.");
        std::vector<RepresentationNode> nodes(ordered.size());
        std::vector<RepresentationNode*> node_ptrs(ordered.size());
        for (size_t i = 0; i < ordered.size(); ++i) {
            nodes[i] = { (uint32_t)ordered[i].key, 0, (int32_t)(i * 2), false };
            node_ptrs[i] = &nodes[i];
        }
        LayerNodeRef move_src, move_target;
        RepresentationDrop move_drop = RepresentationDrop_Before;
        bool move = false;
        st.pending_delete = 0;
        for (LayerRow& row : ordered) {
            if (draw_row(row, true, node_ptrs.data(), node_ptrs.size(), &move_src, &move_drop)) {
                move = true;
                move_target = {NCI_LAYER_LIST_ID, (uint32_t)row.key, false};
            }
        }
        if (move) {
            const int src_idx = representation_tree_find(node_ptrs.data(), node_ptrs.size(), false, move_src.id);
            const int target_idx = representation_tree_find(node_ptrs.data(), node_ptrs.size(), false, move_target.id);
            if (src_idx >= 0 && target_idx >= 0 &&
                layer_tree_move(node_ptrs.data(), node_ptrs.size(), (size_t)src_idx, target_idx, move_drop)) {
                std::sort(nodes.begin(), nodes.end(), [](const RepresentationNode& a, const RepresentationNode& b) {
                    return a.order < b.order;
                });
                st.layer_order.clear();
                for (const RepresentationNode& node : nodes) st.layer_order.push_back((int)node.id);
            }
        }
        if (st.pending_delete) {
            delete_layer(st.pending_delete);
            st.pending_delete = 0;
            return;
        }

        // The selected layer: its name, CP types and the state of each type's surface
        ImGui::Separator();
        std::vector<NciSurface*> members = nci_layer_members(st, st.selected_layer);
        if (members.empty()) {
            ImGui::TextDisabled(st.selected_layer != NCI_LAYER_PREVIEW ? "Select a layer."
                                : critical_points_enabled ? "Select critical points to preview their NCI surface."
                                                          : "Tick Enabled to search for critical points.");
            return;
        }
        if (st.selected_layer > 0) {
            NciSurface& first = *members.front();
            ImGui::SetNextItemWidth(MAX(1.0f, ImGui::GetContentRegionAvail().x - ImGui::CalcTextSize("Delete").x -
                                                  2.0f * ImGui::GetStyle().FramePadding.x - spacing));
            if (ImGui::InputText("##layer name", first.name, sizeof(first.name)))
                for (size_t i = 1; i < members.size(); ++i)
                    snprintf(members[i]->name, sizeof(members[i]->name), "%s", first.name);
            ImGui::SameLine();
            if (ImGui::Button("Delete")) {
                delete_layer(st.selected_layer);
                return;
            }
        }
        ImGui::TextWrapped("Types: %s", nci_layer_types(members).c_str());
        ImGui::AlignTextToFramePadding();
        ImGui::TextDisabled("Scatter");
        ImGui::SameLine();
        draw_scatter_layer_controls(st.selected_layer);
        draw_region_table(region_stats(members));
        for (NciSurface* n : members) draw_member_status(*n);
        std::vector<std::vector<std::string>> environments;
        for (const NciSurface* n : members) environments.insert(environments.end(), n->environments.begin(), n->environments.end());
        draw_environments(environments);
    }

    // ---- what the scatter plot draws of each layer, and how its samples sit in the expected regions
    bool is_shown(int key) const { return key == st.selected_layer && !st.layer_hidden; }

    NciScatterLayer scatter_layer(int key) const {
        const auto it = st.scatter_layers.find(key);
        return it == st.scatter_layers.end() ? NciScatterLayer{} : it->second;
    }
    bool scatter_visible(int key) const { return scatter_layer(key).visible; }
    bool points_on(int key) const {
        const NciScatterLayer l = scatter_layer(key);
        return l.visible && l.points;
    }
    bool lines_on(int key) const {
        const NciScatterLayer l = scatter_layer(key);
        return l.visible && l.lines;
    }
    void set_scatter_visible(int key, bool visible) {
        NciScatterLayer l = scatter_layer(key);
        l.visible = visible;
        st.scatter_layers[key] = l;
    }

    std::vector<int> layer_keys() const {
        std::vector<int> keys = {NCI_LAYER_FULL_SYSTEM, NCI_LAYER_PREVIEW};
        keys.insert(keys.end(), st.layer_order.begin(), st.layer_order.end());
        return keys;
    }

    // The members whose samples the plot draws: the shown layer's ready ones, or another layer's last results
    std::vector<NciSurface*> plotted_members(int key) {
        std::vector<NciSurface*> out;
        for (NciSurface* n : nci_layer_members(st, key)) {
            const bool ok = is_shown(key) ? nci_scale_member(*n) && surface_ready(*n)
                                          : n->shown && n->have_key && !n->cancelled && !n->failed;
            if (ok) out.push_back(n);
        }
        return out;
    }

    struct RegionStats {
        std::array<nci_surface::RegionCounts, nci_surface::family_count> by_type = {};
        size_t inside = 0, outside = 0;
    };
    static RegionStats region_stats(const std::vector<NciSurface*>& members) {
        RegionStats stats;
        for (const NciSurface* n : members) {
            if (!n->families || !n->have_key || n->cancelled || n->failed) continue;
            const auto c = nci_surface::expected_region_counts(n->scatter, n->group);
            stats.by_type[(size_t)n->group].inside += c.inside;
            stats.by_type[(size_t)n->group].outside += c.outside;
            stats.inside += c.inside;
            stats.outside += c.outside;
        }
        return stats;
    }
    static void draw_region_table(const RegionStats& stats) {
        if (stats.inside + stats.outside == 0) return;
        if (!ImGui::BeginTable("##regions", 2, ImGuiTableFlags_SizingFixedFit)) return;
        ImGui::TableSetupColumn("");
        ImGui::TableSetupColumn("in expected region");
        ImGui::TableHeadersRow();
        for (size_t f = 0; f < stats.by_type.size(); ++f) {
            const auto& c = stats.by_type[f];
            if (c.inside + c.outside == 0) continue;
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextColored(nci_type_colour((nci_surface::Family)f), "%s", nci_surface::family_name((nci_surface::Family)f));
            ImGui::TableNextColumn();
            ImGui::Text("%.0f %%", 100.0 * (double)c.inside / (double)(c.inside + c.outside));
        }
        ImGui::EndTable();
    }

    // Points (of a layer), isolines and area in the scatter plot, of a layer or an interaction
    static bool draw_style_controls(bool* points, bool* lines, NciScatterArea* area, const char* area_tooltip) {
        bool changed = false;
        if (points) {
            changed |= ImGui::Checkbox("Points", points);
            ImGui::SameLine();
        }
        changed |= ImGui::Checkbox("Isolines", lines);
        ImGui::SameLine();
        ImGui::SetNextItemWidth(6.5f * ImGui::GetFontSize());
        if (ImGui::BeginCombo("##area", nci_scatter_area_name(*area))) {
            for (int a = 0; a < (int)NciScatterArea::Count; ++a)
                if (ImGui::Selectable(nci_scatter_area_name((NciScatterArea)a), (int)*area == a)) {
                    *area = (NciScatterArea)a;
                    changed = true;
                }
            ImGui::EndCombo();
        }
        ImGui::SetItemTooltip("%s", area_tooltip);
        return changed;
    }

    void draw_scatter_layer_controls(int key) {
        NciScatterLayer l = scatter_layer(key);
        draw_style_controls(&l.points, &l.lines, &l.area,
                            "Of the whole layer, in its colour. Filled: the inside of its outline. Heatmap: the density\n"
                            "of its samples.");
        st.scatter_layers[key] = l;
    }

    // An interaction's isolines and area, the same in every layer; or those of all the interactions at once
    void draw_type_style_controls(nci_surface::Family family, const char* name) {
        static const char* area_tooltip = "In each layer, in the interaction's colour. Filled: the inside of its outline.\n"
                                          "Heatmap: the density of its samples.";
        NciTypeStyle& t = st.type_styles[(size_t)family];
        ImGui::TextUnformatted(name);
        ImGui::PushID("type");
        draw_style_controls(nullptr, &t.lines, &t.area, area_tooltip);
        ImGui::PopID();
        ImGui::Separator();
        ImGui::TextUnformatted("All interactions");
        NciTypeStyle all = t;
        ImGui::PushID("all");
        if (draw_style_controls(nullptr, &all.lines, &all.area, area_tooltip))
            for (NciTypeStyle& other : st.type_styles) {
                other.lines = all.lines;
                other.area = all.area;
            }
        ImGui::PopID();
    }

    // Clicking the active interaction hides the guides; clicking any other makes it the active one, with its guides
    void toggle_guide(NciSurface& n) {
        NciSurface* active = st.layer_hidden ? nullptr : active_member();
        if (active == &n && st.show_guides) st.show_guides = false;
        else {
            st.active_id = n.id;
            st.show_guides = true;
        }
    }

    // Areas in the plot: the inside of the outlines of scatters filled, or their sample density
    static void begin_area(ImDrawList** dl, ImDrawListFlags* flags) {
        ImPlot::PushPlotClipRect();
        *dl = ImPlot::GetPlotDrawList();
        *flags = (*dl)->Flags;
        (*dl)->Flags &= ~ImDrawListFlags_AntiAliasedFill;   // the triangles and cells tile: no seams between them
    }
    static void end_area(ImDrawList* dl, ImDrawListFlags flags) {
        dl->Flags = flags;
        ImPlot::PopPlotClipRect();
    }
    static void fill_area(std::span<const nci_surface::Scatter* const> scatters, float rho_max, ImVec4 colour, int grow) {
        std::vector<std::array<float, 2>> triangles;
        if (!(rho_max > 0.0f) || !nci_surface::scatter_union_fill(&triangles, scatters, rho_max, grow)) return;
        ImDrawList* dl;
        ImDrawListFlags flags;
        begin_area(&dl, &flags);
        const ImU32 col = ImGui::ColorConvertFloat4ToU32(colour);
        for (size_t i = 0; i + 2 < triangles.size(); i += 3)
            dl->AddTriangleFilled(ImPlot::PlotToPixels(triangles[i][0], triangles[i][1]),
                                  ImPlot::PlotToPixels(triangles[i + 1][0], triangles[i + 1][1]),
                                  ImPlot::PlotToPixels(triangles[i + 2][0], triangles[i + 2][1]), col);
        end_area(dl, flags);
    }
    static void heatmap_area(std::span<const nci_surface::Scatter* const> scatters, float rho_max, float alpha) {
        nci_surface::Heatmap map;
        if (!(rho_max > 0.0f) || !nci_surface::scatter_heatmap(&map, scatters, rho_max) || map.bins <= 0) return;
        ImDrawList* dl;
        ImDrawListFlags flags;
        begin_area(&dl, &flags);
        const float w = 2.0f * map.rho_max / map.bins, h = nci_surface::scatter_s_max / map.bins;
        for (int y = 0; y < map.bins; ++y) for (int x = 0; x < map.bins; ++x) {
            const float v = map.values[x + (size_t)map.bins * y];
            if (v < 0.01f) continue;   // empty bins stay transparent, over the context
            ImVec4 colour = ImPlot::SampleColormap(v, ImPlotColormap_Viridis);
            colour.w = alpha;
            const ImVec2 a = ImPlot::PlotToPixels(-map.rho_max + x * w, (y + 1) * h);
            const ImVec2 b = ImPlot::PlotToPixels(-map.rho_max + (x + 1) * w, y * h);
            dl->AddRectFilled(a, b, ImGui::ColorConvertFloat4ToU32(colour));
        }
        end_area(dl, flags);
    }

    // A layer's area, of all its samples together, in its colour or on the heatmap colour map. On the histogram
    // ranges of the outlines: each surface's cutoff, the largest for a layer.
    void draw_layer_area(ImVec4 colour, const std::vector<const NciSurface*>& layer, NciScatterArea area) {
        float rho_max = 0.0f;
        std::vector<const nci_surface::Scatter*> scatters;
        for (const NciSurface* n : layer) {
            rho_max = MAX(rho_max, n->rho_max);
            scatters.push_back(&n->scatter);
        }
        if (scatters.empty()) return;
        if (area == NciScatterArea::Filled) {
            colour.w = 0.25f;
            fill_area(scatters, rho_max, colour, 1);
        } else if (area == NciScatterArea::Heatmap) heatmap_area(scatters, rho_max, 0.8f);
    }

    // ---- the scatter plot: the expected regions, the full-density context and the layers, Ramachandran-like
    void draw_scatter() {
        static bool fit_view = false;
        NciSurface& context = st.context;
        ImGui::Checkbox("Guides", &st.show_guides);
        ImGui::SetItemTooltip("The isovalue and the colour range. Drag the isovalue line to edit it: in Uniform mode\n"
                              "for every interaction, in CP-adaptive mode for the active one.");
        const bool global = st.view == nci_surface::View::Overview && !st.classified_overview;
        if (global) {
            ImGui::SameLine();
            if (ImGui::SmallButton("Global scale...")) ImGui::OpenPopup("##Global scatter scale");
        }
        if (ImGui::BeginPopup("##Global scatter scale")) {
            prop::begin("Global scatter scale");
            color_scale_draw_controls(&context.scale, context.span, md_unit_pow(md_unit_bohr_radius(), -3),
                                      "context samples", "the 1st to 99th percentile of the full density-grid scatter", true);
            prop::end();
            ImGui::EndPopup();
        }
        ImGui::SameLine();
        if (context.rho_failed) ImGui::TextDisabled("Context density failed; see the log.");
        else if (context.job || context.have_key != context.want_key) ImGui::TextDisabled("Computing context samples...");
        else if (context.cancelled) ImGui::TextDisabled("Context cancelled.");
        else if (context.failed) ImGui::TextDisabled("Context extraction failed; see the log.");
        else {
            ImGui::TextDisabled("%zu samples", context.scatter.s.size());
            ImGui::SetItemTooltip("%zu eligible grid samples; a capped volume sample, not an interaction-energy measure.",
                                  context.scatter.eligible_count);
        }
        if (context.job) {
            ImGui::SameLine();
            ImGui::BeginDisabled(context.job->cancel.load(std::memory_order_relaxed) != 0);
            if (ImGui::SmallButton("Cancel")) context.job->cancel.store(1, std::memory_order_relaxed);
            ImGui::EndDisabled();
        } else if (context.rho_failed || context.failed || context.cancelled) {
            ImGui::SameLine();
            if (ImGui::SmallButton("Retry")) {
                context.have_key = 0;
                if (context.rho_failed) context.rho_key = 0;
                context.rho_failed = context.failed = context.cancelled = false;
            }
        }
        // On the right, as in the Ramachandran plot: fitting the view and the settings
        {
            const ImGuiStyle& style = ImGui::GetStyle();
            const float bw = ImGui::CalcTextSize(ICON_FA_EXPAND).x + ImGui::CalcTextSize(ICON_FA_GEAR).x +
                             style.FramePadding.x * 4 + style.ItemSpacing.x;
            ImGui::SameLine();
            ImGui::SetCursorPosX(ImGui::GetCursorPosX() + MAX(0.0f, ImGui::GetContentRegionAvail().x - bw));
            if (ImGui::Button(ICON_FA_EXPAND "##fit")) fit_view = true;
            ImGui::SetItemTooltip("Fit the view (or double-click the plot)");
            ImGui::SameLine();
            if (ImGui::Button(ICON_FA_GEAR "##scatter settings")) ImGui::OpenPopup("##NCI scatter settings");
            ImGui::SetItemTooltip("Scatter plot settings");
            if (ImGui::BeginPopup("##NCI scatter settings")) {
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10.0f);
                ImGui::SliderFloat("Layer line width", &st.layer_outline_width, 0.5f, 10.0f, "%.1f px");
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10.0f);
                ImGui::SliderFloat("Interaction line width", &st.type_outline_width, 0.5f, 8.0f, "%.1f px");
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10.0f);
                ImGui::SliderFloat("Point size", &st.point_size, 0.5f, 6.0f, "%.1f px");
                if (ImGui::Button("Reset")) {
                    st.layer_outline_width = NciState::default_layer_outline_width;
                    st.type_outline_width = NciState::default_type_outline_width;
                    st.point_size = NciState::default_point_size;
                }
                ImGui::Separator();
                ImGui::TextDisabled("Legend: click a swatch to show or hide it, an interaction's name to make it\n"
                                    "the active one for the isovalue guide. Right-click a layer for how it is\n"
                                    "drawn.");
                ImGui::EndPopup();
            }
        }

        // The legend is a row of chips under the plot, so it never covers the data: the expected regions and the layers,
        // then the interactions
        const std::vector<NciSurface*> members = nci_layer_members(st, st.selected_layer);
        const ImVec4 layer_colour = layer_chip(st.selected_layer, members);
        const std::vector<int> keys = layer_keys();
        struct LegendEntry {
            std::string text;
            ImVec4 colour;
            enum Kind { Reference, Context, Layer, Type } kind = Layer;
            int layer = 0;
            nci_surface::Family family = nci_surface::Family::Count;
            NciSurface* surface = nullptr;   // the shown layer's surface of that interaction
            bool on = true;
            bool dim = false;
            bool new_line = false;
        };
        std::vector<LegendEntry> legend;
        legend.push_back({"Expected regions", ImVec4(0.85f, 0.85f, 0.85f, 1.0f), LegendEntry::Reference, 0,
                          nci_surface::Family::Count, nullptr, st.show_reference, !st.show_reference});
        if (surface_ready(context))
            legend.push_back({"All samples", ImVec4(0.8f, 0.8f, 0.8f, 1.0f), LegendEntry::Context, 0,
                              nci_surface::Family::Count, nullptr, st.show_context_samples, !st.show_context_samples});
        std::array<bool, nci_surface::family_count> present = {};
        std::vector<NciSurface*> shown;
        for (int key : keys) {
            const auto layer = nci_layer_members(st, key);
            if (layer.empty()) continue;
            LegendEntry e = {layer_name(key, layer), layer_chip(key, layer), LegendEntry::Layer, key};
            e.on = scatter_visible(key);
            e.dim = !e.on;
            legend.push_back(e);
            if (is_shown(key)) shown = plotted_members(key);
            if (!e.dim)
                for (const NciSurface* n : plotted_members(key)) if (n->families) present[(size_t)n->group] = true;
        }
        bool first_type = true;
        for (size_t f = 0; f < present.size(); ++f) {
            if (!present[f]) continue;
            LegendEntry e;
            e.kind = LegendEntry::Type;
            e.family = (nci_surface::Family)f;
            e.text = nci_surface::family_name(e.family);
            e.colour = nci_type_colour(e.family);
            e.on = st.type_styles[f].visible;
            e.dim = !e.on;
            e.new_line = first_type;
            first_type = false;
            for (NciSurface* n : shown) if (n->families && n->group == e.family) e.surface = n;
            legend.push_back(e);
        }
        const ImGuiStyle& style = ImGui::GetStyle();
        const float chip_h = ImGui::GetTextLineHeight();
        const float avail_w = ImGui::GetContentRegionAvail().x;
        const auto swatch_w = [&](const LegendEntry& e) { return e.kind == LegendEntry::Layer ? 2.0f * chip_h : chip_h; };
        // Chips wrapped into centred lines: each line's width and first chip
        std::vector<float> line_widths;
        std::vector<size_t> line_starts;
        for (size_t i = 0; i < legend.size(); ++i) {
            const float w = swatch_w(legend[i]) + style.ItemInnerSpacing.x + ImGui::CalcTextSize(legend[i].text.c_str()).x;
            if (line_widths.empty() || legend[i].new_line || line_widths.back() + style.ItemSpacing.x + w > avail_w) {
                line_widths.push_back(w);
                line_starts.push_back(i);
            } else line_widths.back() += style.ItemSpacing.x + w;
        }
        const float legend_h = line_widths.size() * ImGui::GetFrameHeightWithSpacing();
        if (fit_view) ImPlot::SetNextAxesToFit();
        fit_view = false;
        const bool plot_open = ImPlot::BeginPlot("##NCI scatter", ImVec2(-1, MAX(1.0f, ImGui::GetContentRegionAvail().y - legend_h)),
                                                 ImPlotFlags_NoLegend);
        if (!plot_open) return;
        ImPlot::SetupAxes("sign(\xCE\xBB\xE2\x82\x82)\xCF\x81 [bohr^-3]", "s");
        ImPlot::SetupAxisLimits(ImAxis_X1, -context.rho_max, context.rho_max, ImPlotCond_Once);
        ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, nci_surface::scatter_s_max, ImPlotCond_Once);

        struct Points {
            const nci_surface::Scatter* scatter;
            const std::vector<uint32_t>* indices;
        };
        const auto plot_bin = [&](const NciSurface& n, size_t bin, const std::vector<uint32_t>& indices, ImVec4 colour, float size) {
            if (indices.empty()) return;
            const auto getter = [](int i, void* data) {
                const auto& p = *(const Points*)data;
                const uint32_t idx = (*p.indices)[i];
                return ImPlotPoint(p.scatter->signed_rho[idx], p.scatter->s[idx]);
            };
            Points points = {&n.scatter, &indices};
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, size, colour, 0.0f, colour);
            char label[64];
            snprintf(label, sizeof(label), "##NCI %d bin %zu", n.id, bin);
            ImPlot::PlotScatterG(label, getter, &points, (int)indices.size());
        };
        // The samples of a surface, on the colour scale of its surface in the viewport; the context is neutral unless
        // the global scale applies
        const auto plot = [&](const NciSurface& n, bool neutral) {
            const auto& scatter = n.scatter;
            const ColorScale scale = n.id == -1 ? nci_effective_scale(st, st.primary) : nci_effective_scale(st, n);
            std::array<std::vector<uint32_t>, 65> bins;
            for (size_t i = 0; i < scatter.s.size(); ++i) {
                const int bin = neutral ? 0 : MIN(64, (int)(64.0f * color_scale_param(scale, scatter.signed_rho[i]) + 0.5f));
                bins[bin].push_back((uint32_t)i);
            }
            for (size_t b = 0; b < bins.size(); ++b) {
                ImVec4 colour = neutral ? ImVec4(0.8f, 0.8f, 0.8f, 0.25f)
                                        : ImPlot::SampleColormap((float)b / (float)(bins.size() - 1), scale.colormap);
                if (!neutral) colour.w = 0.7f;
                plot_bin(n, b, bins[b], colour, neutral ? 1.0f : st.point_size);
            }
        };
        if (st.show_context_samples && surface_ready(context)) plot(context, !global);

        const auto draw_outline = [&](const std::vector<std::array<float, 4>>& outline, ImVec4 color, float width) {
            ImPlot::PushPlotClipRect();
            for (const auto& edge : outline)
                ImPlot::GetPlotDrawList()->AddLine(ImPlot::PlotToPixels(edge[0], edge[1]), ImPlot::PlotToPixels(edge[2], edge[3]),
                                                  ImGui::ColorConvertFloat4ToU32(color), width);
            ImPlot::PopPlotClipRect();
        };

        // Each interaction's expected region, those of the plotted interactions (all when none is)
        if (st.show_reference) {
            const bool any = std::any_of(present.begin(), present.end(), [](bool p) { return p; });
            ImPlot::PushPlotClipRect();
            ImDrawList* dl = ImPlot::GetPlotDrawList();
            for (size_t f = 0; f < present.size(); ++f) {
                const nci_surface::ExpectedRegion r = nci_surface::expected_region((nci_surface::Family)f);
                if (!r.valid || (any && !present[f]) || !st.type_styles[f].visible) continue;
                ImVec4 colour = nci_type_colour((nci_surface::Family)f);
                const ImVec2 a = ImPlot::PlotToPixels(r.rho_min, r.s_max), b = ImPlot::PlotToPixels(r.rho_max, 0.0);
                const ImVec2 lo(MIN(a.x, b.x), MIN(a.y, b.y)), hi(MAX(a.x, b.x), MAX(a.y, b.y));
                colour.w = 0.10f;
                dl->AddRectFilled(lo, hi, ImGui::ColorConvertFloat4ToU32(colour));
                colour.w = 0.45f;
                dl->AddRect(lo, hi, ImGui::ColorConvertFloat4ToU32(colour), 0.0f, 0, 1.0f);
            }
            ImPlot::PopPlotClipRect();
        }

        NciSurface* hovered = nullptr;
        size_t hovered_sample = 0;
        float closest = 36.0f;
        const auto layer_outline = [&](int key) -> const std::vector<std::array<float, 4>>* {
            const auto it = st.layer_outlines.find(key);
            return it == st.layer_outlines.end() ? nullptr : &it->second.segments;
        };
        // A layer's area, samples, then outlines: thick for the whole layer, thin on top for each interaction
        const auto draw_layer = [&](int key) {
            if (!scatter_visible(key)) return;
            const auto layer = plotted_members(key);
            if (layer.empty()) return;
            const NciScatterLayer style_of = scatter_layer(key);
            // The interactions' isolines and areas; the points are the layer's, whatever the interactions
            const auto type_style = [&](const NciSurface& n) -> NciTypeStyle {
                if (n.families) return st.type_styles[(size_t)n.group];
                NciTypeStyle none;
                none.visible = false;
                return none;
            };
            if (style_of.area != NciScatterArea::None) {
                const std::vector<const NciSurface*> whole(layer.begin(), layer.end());
                draw_layer_area(layer_chip(key, nci_layer_members(st, key)), whole, style_of.area);
            }
            for (const NciSurface* n : layer) {
                const NciTypeStyle t = type_style(*n);
                if (t.visible && t.area != NciScatterArea::None) draw_layer_area(type_colour(*n), {n}, t.area);
            }
            for (NciSurface* n : layer) {
                if (!style_of.points) break;
                plot(*n, false);
                if (!ImPlot::IsPlotHovered()) continue;
                const ImVec2 mouse = ImGui::GetMousePos();
                for (size_t i = 0; i < n->scatter.s.size(); ++i) {
                    const ImVec2 p = ImPlot::PlotToPixels(n->scatter.signed_rho[i], n->scatter.s[i]);
                    const float d2 = (p.x - mouse.x) * (p.x - mouse.x) + (p.y - mouse.y) * (p.y - mouse.y);
                    if (d2 < closest) { closest = d2; hovered = n; hovered_sample = i; }
                }
            }
            const auto* outline = lines_on(key) ? layer_outline(key) : nullptr;
            if (outline) {
                ImVec4 colour = layer_chip(key, nci_layer_members(st, key));
                colour.w = 0.6f;
                draw_outline(*outline, colour, st.layer_outline_width);
            }
            for (const NciSurface* n : layer) {
                const NciTypeStyle t = type_style(*n);
                if (!t.visible || !t.lines) continue;
                ImVec4 colour = type_colour(*n);
                colour.w = 1.0f;
                draw_outline(n->outline, colour, st.type_outline_width);
            }
        };
        // The layer shown in the viewport on top
        for (int key : keys) if (key != st.selected_layer) draw_layer(key);
        draw_layer(st.selected_layer);

        NciSurface* active = st.layer_hidden ? nullptr : active_member();
        bool guide_hovered = false;
        if (st.show_guides && active && nci_scale_member(*active)) {
            const ColorScale scale = nci_effective_scale(st, *active);
            const bool adaptive = st.mode == NciParameterMode::CPAdaptive && active->families;
            // Uniform: one isovalue for the layer; CP-adaptive: the active interaction's
            const ImVec4 col = adaptive ? type_colour(*active) : layer_colour;
            double iso = active->s_iso;
            bool held = false, hovered_line = false;
            if (ImPlot::DragLineY(1, &iso, col, 2.0f, 0, nullptr, &hovered_line, &held)) {
                float& target = adaptive ? st.family_settings[nci_settings_index(*active)].s_iso : st.uniform.s_iso;
                target = (float)std::clamp(iso, 0.1, 1.5);
            }
            st.settings_editing |= held;
            guide_hovered = hovered_line || held;
            const double limits[] = {scale.range_beg, scale.range_end};
            ImPlot::SetNextLineStyle(col, 1.0f);
            ImPlot::PlotInfLines("##Colour limits", limits, 2);
            const double right = ImPlot::GetPlotLimits().X.Max;
            if (adaptive) ImPlot::Annotation(right, active->s_iso, col, ImVec2(-4, -4), true, "%s: s=%.2f", nci_family_label(*active).c_str(), active->s_iso);
            else ImPlot::Annotation(right, active->s_iso, col, ImVec2(-4, -4), true, "s=%.2f", active->s_iso);
            ImPlot::TagX(scale.range_beg, col, "%.4g", scale.range_beg);
            ImPlot::TagX(scale.range_end, col, "%.4g", scale.range_end);
        }
        if (hovered && !guide_hovered) {
            const int key = nci_layer_key(*hovered);
            const float rho = hovered->scatter.signed_rho[hovered_sample], s = hovered->scatter.s[hovered_sample];
            const nci_surface::ExpectedRegion region = nci_surface::expected_region(hovered->group);
            ImGui::BeginTooltip();
            ImGui::Text("%s | %s", layer_name(key, nci_layer_members(st, key)), nci_family_label(*hovered).c_str());
            ImGui::Text("sign(lambda2)rho %.5g au | s %.3f", rho, s);
            if (hovered->families && region.valid && s <= region.s_max)
                ImGui::TextDisabled(nci_surface::in_expected_region(region, rho, s) ? "In the expected region"
                                                                                     : "Outside the expected region");
            const bool guided = st.show_guides && active == hovered;
            if (is_shown(key))
                ImGui::TextDisabled(guided ? "Click to hide the guides" : "Click to make this interaction the active one for the guides");
            ImGui::EndTooltip();
            if (is_shown(key) && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) toggle_guide(*hovered);
        }
        ImPlot::EndPlot();

        const float left = ImGui::GetCursorPosX();
        for (size_t i = 0, line = 0; i < legend.size(); ++i) {
            const LegendEntry& e = legend[i];
            if (line < line_starts.size() && line_starts[line] == i) {
                ImGui::SetCursorPosX(left + MAX(0.0f, 0.5f * (avail_w - line_widths[line])));
                ++line;
            } else ImGui::SameLine();
            ImGui::PushID((int)i);
            ImVec4 colour = e.colour;
            if (!e.on) colour.w = 0.2f;
            // The swatch toggles: a square for an interaction, a box for the expected regions, a thick bar for a layer
            const ImVec2 p = ImGui::GetCursorScreenPos();
            const bool toggled = ImGui::InvisibleButton("##swatch", ImVec2(swatch_w(e), chip_h));
            ImDrawList* dl = ImGui::GetWindowDrawList();
            const ImU32 col = ImGui::ColorConvertFloat4ToU32(colour);
            if (e.kind == LegendEntry::Layer) {
                const float h = MIN(chip_h, MAX(2.0f, st.layer_outline_width));
                dl->AddRectFilled(ImVec2(p.x, p.y + 0.5f * (chip_h - h)), ImVec2(p.x + swatch_w(e), p.y + 0.5f * (chip_h + h)), col, 0.5f * h);
                ImGui::SetItemTooltip("Click: %s in the scatter plot\nRight-click: points, isolines and area",
                                      e.on ? "hide" : "show");
                if (ImGui::BeginPopupContextItem("##layer scatter")) {
                    ImGui::TextUnformatted(e.text.c_str());
                    draw_scatter_layer_controls(e.layer);
                    ImGui::EndPopup();
                }
            } else if (e.kind == LegendEntry::Reference) {
                ImVec4 fill = colour;
                fill.w *= 0.25f;
                dl->AddRectFilled(p, ImVec2(p.x + chip_h, p.y + chip_h), ImGui::ColorConvertFloat4ToU32(fill));
                dl->AddRect(p, ImVec2(p.x + chip_h, p.y + chip_h), col);
                ImGui::SetItemTooltip("%s the expected regions of the interactions\n(provisional, from the calibration systems)",
                                      e.on ? "Hide" : "Show");
            } else if (e.kind == LegendEntry::Context) {
                dl->AddCircleFilled(ImVec2(p.x + 0.5f * chip_h, p.y + 0.5f * chip_h), 0.3f * chip_h, col);
                ImGui::SetItemTooltip("%s every grid sample of the whole density, under the layers", e.on ? "Hide" : "Show");
            } else {
                dl->AddRectFilled(p, ImVec2(p.x + chip_h, p.y + chip_h), col, style.FrameRounding);
                ImGui::SetItemTooltip("Click: %s this interaction's isolines, area and expected region in every layer\n"
                                      "Right-click: its isolines and area, or those of all the interactions",
                                      e.on ? "hide" : "show");
                if (ImGui::BeginPopupContextItem("##type scatter")) {
                    draw_type_style_controls(e.family, e.text.c_str());
                    ImGui::EndPopup();
                }
            }
            if (toggled) {
                if (e.kind == LegendEntry::Type) st.type_styles[(size_t)e.family].visible = !e.on;
                else if (e.kind == LegendEntry::Reference) st.show_reference = !e.on;
                else if (e.kind == LegendEntry::Context) st.show_context_samples = !e.on;
                else {
                    set_scatter_visible(e.layer, !e.on);
                }
            }
            ImGui::SameLine(0.0f, style.ItemInnerSpacing.x);
            ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(e.dim ? ImGuiCol_TextDisabled : ImGuiCol_Text));
            if (e.surface) {
                // The name makes the interaction the active one for the isovalue guide
                const bool guided = e.surface == active && st.show_guides;
                if (ImGui::Selectable(e.text.c_str(), guided, 0, ImGui::CalcTextSize(e.text.c_str()))) toggle_guide(*e.surface);
                ImGui::SetItemTooltip(guided ? "Click to hide the isovalue and colour-range guides"
                                             : "Click to make this interaction the active one for the guides");
            } else ImGui::TextUnformatted(e.text.c_str());
            ImGui::PopStyleColor();
            ImGui::PopID();
        }
    }
};

} // namespace

void draw_nci_panel(NciState& state, bool critical_points_enabled) {
    NciPanel panel(state, critical_points_enabled);
    const topology::NciPanelDrawCallbacks callbacks = {
        &panel,
        [](void* data) { static_cast<NciPanel*>(data)->draw_scatter(); },
        [](void* data) { static_cast<NciPanel*>(data)->draw_layers(); },
        [](void* data) { static_cast<NciPanel*>(data)->st.selection_key = 0; },
    };
    draw_nci_panel_layout(&state.enabled, callbacks);
}

} // namespace topology
