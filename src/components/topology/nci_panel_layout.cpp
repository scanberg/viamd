#include "nci_panel_layout.h"

#include "critical_point_distribution_core.h"

#include <color_scale.h>
#include <core/md_common.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <density_map.h>
#include <gfx/gl_utils.h>
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

    if (!*enabled) {
        ImGui::TextDisabled("Enable NCI to show the scatter plot and layers.");
        return;
    }

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
        "How the isovalue, density cutoff and radius apply to every NCI layer, the full system included.\n"
        "Uniform: standard NCI, one isovalue, cutoff and radius around the CPs for every CP type.\n"
        "CP-adaptive: the same isovalue and cutoff, with each CP type's own radius around its CPs, so\n"
        "extended interactions such as pi-stacking are covered and compact ones stay local.\n"
        "One colour scale applies in both modes.";
    static constexpr size_t max_samples = 32u << 20;

    NciState& st;
    bool critical_points_enabled;
    std::vector<uint32_t>& critical_point_list_scope;
    std::vector<uint32_t>& critical_point_environment_selection;

    NciPanel(NciState& state, bool cp_enabled, std::vector<uint32_t>& cp_list_scope,
             std::vector<uint32_t>& cp_environment_selection)
        : st(state), critical_points_enabled(cp_enabled), critical_point_list_scope(cp_list_scope),
          critical_point_environment_selection(cp_environment_selection) {}

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

    // A layer is named by the CP families of its CPs, their branch of the CP tree: "Hydrogen bonds", "Aromatic rings",
    // "Hydrogen bonds + \xCF\x80-stacking", "Hydrogen bonds + 2 more"; numbered when another layer has the name
    std::string selection_layer_name(const std::vector<NciSurface*>& members, int number) const {
        std::vector<std::string> families;
        for (const NciSurface* member : members) {
            for (const auto& path : member->environments) {
                if (path.empty()) continue;
                const bool branched = (path[0] == "Contacts" || path[0] == "Rings") && path.size() > 1;
                const std::string& family = branched ? path[1] : path[0];
                if (std::find(families.begin(), families.end(), family) == families.end()) families.push_back(family);
            }
        }
        if (families.empty()) return "Layer " + std::to_string(number);
        std::string name = families.front();
        if (families.size() == 2) name += " + " + families[1];
        else if (families.size() > 2) name += " + " + std::to_string(families.size() - 1) + " more";
        const auto taken = [&](const std::string& candidate) {
            return std::any_of(st.objects.begin(), st.objects.end(), [&](const auto& n) { return candidate == n->name; });
        };
        if (taken(name)) {
            for (int k = 2;; ++k) {
                const std::string candidate = name + " " + std::to_string(k);
                if (!taken(candidate)) {
                    name = candidate;
                    break;
                }
            }
        }
        if (name.size() >= sizeof(NciSurface::name)) {
            size_t end = sizeof(NciSurface::name) - 1;
            while (end && (static_cast<unsigned char>(name[end]) & 0xC0) == 0x80) --end;
            name.resize(end);
        }
        return name;
    }

    std::unique_ptr<NciSurface> copy_preview(const NciSurface& source, int group_id, const std::string& name) {
        auto n = std::make_unique<NciSurface>();
        n->saved = true;
        n->id = st.next_id++;
        n->object_group = group_id;
        snprintf(n->name, sizeof(n->name), "%s", name.c_str());
        n->shown = true;
        n->layer_visible = true;
        n->s_iso = source.s_iso;
        n->scatter_s = source.scatter_s;
        n->rho_max = source.rho_max;
        n->scale = source.scale;
        n->families = source.families;
        n->type_counts = source.type_counts;
        n->critical_point_indices = source.critical_point_indices;
        n->critical_points_by_family = source.critical_points_by_family;
        n->reference_scatters = source.reference_scatters;
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
        n->reference_failed = source.reference_failed;
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
        const std::string name = selection_layer_name(nci_layer_members(st, NCI_LAYER_PREVIEW), number);
        for (const auto& n : st.previews) {
            if (!n->shown) continue;
            auto saved = copy_preview(*n, group_id, name);
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
        prop::row("isovalue s", "Reduced density gradient isovalue, the same for every CP type.");
        ImGui::SetNextItemWidth(prop::control_width());
        ImGui::SliderFloat("##NCI uniform s", &st.uniform.s_iso, 0.1f, 1.5f, "%.2f");
        editing |= ImGui::IsItemActive();
        prop::row("scatter s", "How far in s the scatter plot's connected regions reach, independent of the surface "
                               "isovalue; at the reference's 1.0 the scatter covers what the reference map does. "
                               "Without connected regions the scatter holds every sample up to s = 2 within the radius.");
        ImGui::BeginDisabled(!st.connected);
        ImGui::SetNextItemWidth(prop::control_width());
        ImGui::SliderFloat("##NCI scatter s", &st.uniform.scatter_s, 0.1f, 2.0f, "%.2f");
        editing |= ImGui::IsItemActive();
        ImGui::EndDisabled();
        prop::row("\xCF\x81 cutoff", "Density cutoff, the same for every CP type.");
        ImGui::SetNextItemWidth(prop::control_width());
        ImGui::SliderFloat("##NCI uniform rho", &st.uniform.rho_max, 0.005f, 0.2f, "< %.3f", ImGuiSliderFlags_Logarithmic);
        editing |= ImGui::IsItemActive();
        if (!adaptive) {
            prop::row("radius", "Extent around the CPs, the same for every CP type (with the radius cap).");
            ImGui::SetNextItemWidth(prop::control_width());
            ImGui::SliderFloat("##NCI uniform radius", &st.uniform.radius, 0.5f, 12.0f, "%.1f Bohr");
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
        ImGui::SetItemTooltip("Limit the extent around the CPs to the radius. Without a cap, connected regions may merge across the system.");
        prop::row("opacity");
        ImGui::SetNextItemWidth(prop::control_width());
        ImGui::SliderFloat("##NCI opacity", &st.opacity, 0.05f, 1.0f, "%.2f");
        color_scale_draw_controls(&st.scale, st.span, md_unit_pow(md_unit_bohr_radius(), -3), "shown layer",
                                  "the 1st to 99th percentile of the shown layer", true);
        prop::end();

        if (adaptive) {
            ImGui::SeparatorText("Radius per CP type");
            const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_SizingFixedFit;
            if (ImGui::BeginTable("##NCI type settings", 2, flags)) {
                const float em = ImGui::GetFontSize();
                ImGui::TableSetupColumn("Type");
                ImGui::TableSetupColumn("Radius", ImGuiTableColumnFlags_WidthFixed, 6.0f * em);
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
                        ImGui::SetItemTooltip("Uncalibrated type: provisional radius");
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    ImGui::SliderFloat("##radius", &fs.radius, 0.5f, 12.0f, "%.1f Bohr");
                    editing |= ImGui::IsItemActive();
                    ImGui::PopID();
                }
                ImGui::EndTable();
            }
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

    // ---- the rows, as the Ramachandran layers': the eye, a swatch of how it is drawn, the name and a line below it
    struct RowFrame {
        ImVec2 min;
        float w = 0.0f, h = 0.0f, line_h = 0.0f;
        ImVec2 swatch_min, swatch_max;
        ImVec2 text;
        ImVec4 clip;
        bool hovered = false, clicked = false;
    };
    // The row's selectable and its eye; the caller draws the swatch and text into the frame
    RowFrame begin_row(bool selected, bool* visible, const char* eye_tooltip) {
        RowFrame r;
        r.line_h = ImGui::GetTextLineHeight();
        r.h = r.line_h * 2.0f + 8.0f;
        r.min = ImGui::GetCursorScreenPos();
        r.w = ImGui::GetContentRegionAvail().x;
        ImGui::SetNextItemAllowOverlap();
        r.clicked = ImGui::Selectable("##row", selected, ImGuiSelectableFlags_AllowOverlap, ImVec2(r.w, r.h));
        r.hovered = ImGui::IsItemHovered(ImGuiHoveredFlags_DelayNormal | ImGuiHoveredFlags_NoSharedDelay) &&
                    !ImGui::IsMouseDragging(ImGuiMouseButton_Left);
        const float pad = 4.0f;
        ImGui::SetCursorScreenPos(ImVec2(r.min.x + pad, r.min.y + (r.h - ImGui::GetFrameHeight()) * 0.5f));
        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(2, ImGui::GetStyle().FramePadding.y));
        if (ImGui::Button(*visible ? ICON_FA_EYE "##vis" : ICON_FA_EYE_SLASH "##vis")) *visible = !*visible;
        ImGui::PopStyleVar();
        ImGui::PopStyleColor();
        ImGui::SetItemTooltip("%s", eye_tooltip);
        r.swatch_min = ImVec2(r.min.x + pad + r.line_h + 8, r.min.y + (r.h - r.line_h) * 0.5f);
        r.swatch_max = ImVec2(r.swatch_min.x + r.line_h, r.swatch_min.y + r.line_h);
        r.text = ImVec2(r.min.x + pad + r.line_h + 6 + r.line_h + 6, r.min.y + 4);
        r.clip = ImVec4(r.text.x, r.min.y, r.min.x + r.w - pad, r.min.y + r.h);
        return r;
    }
    void row_text(const RowFrame& r, const char* name, const char* detail, bool on, bool selected) {
        ImDrawList* dl = ImGui::GetWindowDrawList();
        const ImU32 text_col = ImGui::GetColorU32(on ? ImGuiCol_Text : ImGuiCol_TextDisabled);
        // On the selection's colour the disabled one is too faint to read
        const ImU32 dim_col = selected ? ImGui::GetColorU32(ImGuiCol_Text, 0.7f) : ImGui::GetColorU32(ImGuiCol_TextDisabled);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize(), r.text, text_col, name, nullptr, 0.0f, &r.clip);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize() * 0.85f, ImVec2(r.text.x, r.text.y + r.line_h + 1), dim_col, detail,
                    nullptr, 0.0f, &r.clip);
    }
    static void end_row(const RowFrame& r) {
        ImGui::SetCursorScreenPos(ImVec2(r.min.x, r.min.y + r.h + 2));
        ImGui::Dummy(ImVec2(0, 0));
    }

    // The reference's swatch: its four categories as it draws them, a quarter each
    void reference_swatch(ImDrawList* dl, ImVec2 mn, ImVec2 mx, float alpha) const {
        const NciReferenceLayer& r = st.reference;
        const ImVec2 mid((mn.x + mx.x) * 0.5f, (mn.y + mx.y) * 0.5f);
        const float s = mx.y - mn.y;
        const ImVec2 centres[4] = {ImVec2(mn.x + 0.3f * s, mn.y + 0.35f * s), ImVec2(mn.x + 0.7f * s, mn.y + 0.35f * s),
                                   ImVec2(mn.x + 0.3f * s, mn.y + 0.7f * s), ImVec2(mn.x + 0.7f * s, mn.y + 0.7f * s)};
        for (size_t c = 0; c < nci_surface::reference_category_count; ++c) {
            ImVec4 col = nci_type_colour(reference_colours[c]);
            col.w = alpha * r.opacity;
            if (r.fill) {
                ImVec4 f = col;
                f.w *= r.style == NciReferenceStyle::Combined ? 0.35f : 0.6f;
                if (r.style == NciReferenceStyle::Combined) dl->AddCircleFilled(centres[c], 0.32f * s, ImGui::ColorConvertFloat4ToU32(f));
                else {
                    const ImVec2 a(c % 2 ? mid.x : mn.x, c < 2 ? mn.y : mid.y), b(c % 2 ? mx.x : mid.x, c < 2 ? mid.y : mx.y);
                    dl->AddRectFilled(a, b, ImGui::ColorConvertFloat4ToU32(f), c == 0 ? 3.0f : 0.0f);
                }
            }
            if (r.lines || !r.fill) dl->AddCircle(centres[c], 0.2f * s, ImGui::ColorConvertFloat4ToU32(col), 0, 1.2f);
        }
    }

    bool draw_row(LayerRow& row, bool saved, RepresentationNode* const* nodes, size_t node_count,
                  LayerNodeRef* move_src, RepresentationDrop* move_drop) {
        bool dropped = false;
        ImGui::PushID(row.key);
        const bool selected = !st.reference_selected && st.selected_layer == row.key;
        const bool in_viewport = st.selected_layer == row.key && !st.layer_hidden;
        NciScatterLayer& style = scatter_style(row.key);
        const RowFrame r = begin_row(selected, &style.visible,
                                     style.visible ? "Hide it in the scatter plot" : "Show it in the scatter plot");
        // A click selects it, shown in the viewport (one layer at a time); clicking the one shown hides it there
        if (r.clicked) {
            if (selected && in_viewport) {
                st.layer_hidden = true;
                nci_apply_settings(st);
            } else select_layer(row.key);
            st.reference_selected = false;
            critical_point_environment_selection.clear();
            update_critical_point_list_scope(row.key);
        }
        if (saved) {
            const ImRect rect(r.min, ImVec2(r.min.x + r.w, r.min.y + r.h));
            layer_tree_drag_source({NCI_LAYER_LIST_ID, (uint32_t)row.key, false}, row.label);
            if (layer_tree_drop_target(rect, r.min.x, {NCI_LAYER_LIST_ID, (uint32_t)row.key, false},
                                       false, nodes, node_count, move_src, move_drop))
                dropped = true;
        }
        if (ImGui::BeginPopupContextItem("##layer menu")) {
            if (ImGui::MenuItem(in_viewport ? "Hide in the viewport" : "Show in the viewport")) {
                if (in_viewport) {
                    st.layer_hidden = true;
                    nci_apply_settings(st);
                } else select_layer(row.key);
            }
            if (saved) {
                ImGui::Separator();
                if (ImGui::MenuItem(ICON_FA_XMARK " Remove")) st.pending_delete = row.key;
            }
            ImGui::EndPopup();
        }

        ImDrawList* dl = ImGui::GetWindowDrawList();
        const DensityMapSwatch kind = style.style == NciLayerStyle::Heatmap ? DensityMapSwatch_Heatmap : DensityMapSwatch_Contours;
        ImVec4 colour = layer_colour(row.key, row.members);
        colour.w *= style.opacity;
        density_map_swatch(dl, r.swatch_min, r.swatch_max, kind, colour, style.fill, ImPlotColormap_Viridis,
                           style.visible ? 1.0f : 0.35f);

        const std::string types = nci_layer_types(row.members);
        size_t triangles = 0;
        float progress = -1.0f;
        for (const NciSurface* n : row.members) {
            triangles += n->triangles;
            if (n->job) progress = MAX(progress, n->job->progress.load(std::memory_order_relaxed));
        }
        // As the Ramachandran layers: where this layer's samples fall inside the reference regions of their categories
        const RegionStats stats = region_stats(row.members);
        char expected[64] = "";
        if (stats.inside + stats.outside)
            snprintf(expected, sizeof(expected), " \xC2\xB7 %.0f %% in reference", 100.0 * (double)stats.inside / (double)(stats.inside + stats.outside));
        char detail[256];
        if (progress >= 0.0f) snprintf(detail, sizeof(detail), "%s \xC2\xB7 computing %.0f %%", types.c_str(), 100.0f * progress);
        else if (in_viewport) snprintf(detail, sizeof(detail), "%s \xC2\xB7 %zu triangles%s", types.c_str(), triangles, expected);
        else snprintf(detail, sizeof(detail), "%s%s", types.c_str(), expected);
        row_text(r, row.label, detail, in_viewport || style.visible, selected);
        if (r.hovered) {
            ImGui::BeginTooltip();
            ImGui::TextUnformatted(row.label);
            ImGui::TextDisabled("%s", detail);
            draw_region_table(stats);
            ImGui::TextDisabled(in_viewport ? "Shown in the viewport; click again to hide it there"
                                            : "Click to show it in the viewport (one layer at a time)");
            ImGui::EndTooltip();
        }
        end_row(r);
        ImGui::PopID();
        return dropped;
    }

    // The built-in reference, a row as the Ramachandran reference's; independent of the loaded system
    void draw_reference_row() {
        ImGui::PushID("##reference");
        const bool selected = st.reference_selected;
        const RowFrame r = begin_row(selected, &st.reference.visible,
                                     st.reference.visible ? "Hide the reference in the scatter plot"
                                                          : "Show the reference in the scatter plot");
        if (r.clicked) st.reference_selected = true;
        reference_swatch(ImGui::GetWindowDrawList(), r.swatch_min, r.swatch_max, st.reference.visible ? 1.0f : 0.35f);
        const nci_surface::ReferenceMap& map = nci_surface::builtin_reference_map();
        int calculations = 0;
        for (const auto& c : map.categories) calculations += c.valid ? c.systems : 0;
        char detail[128];
        snprintf(detail, sizeof(detail), "H-bonds, rings, \xCF\x80, other \xC2\xB7 %s",
                 st.reference.style == NciReferenceStyle::Combined ? "combined" : "dominant");
        row_text(r, "Reference", detail, st.reference.visible, selected);
        if (r.hovered) {
            ImGui::BeginTooltip();
            ImGui::TextUnformatted("Reference");
            ImGui::TextDisabled("Built in, independent of the loaded system, as the Ramachandran reference: per CP category,\n"
                                "the pooled (sign(lambda2)rho, s) density of reference extractions (s <= %.1f, rho <= %.2f,\n"
                                "the category's radius) over many calculations, each weighing the same. A sample is in the\n"
                                "reference inside its category's %.0f %% area.",
                                nci_surface::reference_s_iso, nci_surface::reference_rho_max,
                                100.0f * nci_surface::reference_inside_fraction);
            if (*nci_surface::reference_data().source) ImGui::TextDisabled("%s", nci_surface::reference_data().source);
            if (ImGui::BeginTable("##regions", 2, ImGuiTableFlags_SizingFixedFit)) {
                ImGui::TableSetupColumn("");
                ImGui::TableSetupColumn("Calculations");
                ImGui::TableHeadersRow();
                for (size_t c = 0; c < map.categories.size(); ++c) {
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::TextColored(nci_type_colour(reference_colours[c]), "%s", reference_labels[c]);
                    ImGui::TableNextColumn();
                    if (map.categories[c].valid) ImGui::Text("%d", map.categories[c].systems);
                    else ImGui::TextDisabled("None");
                }
                ImGui::EndTable();
            }
            ImGui::EndTooltip();
        }
        end_row(r);
        ImGui::PopID();
    }

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

        // Below the list, as in the Ramachandran window: what is selected, its properties
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();
        if (ImGui::IsWindowFocused(ImGuiFocusedFlags_ChildWindows) && !ImGui::GetIO().WantTextInput &&
            ImGui::IsKeyPressed(ImGuiKey_Delete) && !st.reference_selected && st.selected_layer > 0) {
            delete_layer(st.selected_layer);
            return;
        }
        if (st.reference_selected) {
            draw_reference_inspector();
            return;
        }
        std::vector<NciSurface*> members = nci_layer_members(st, st.selected_layer);
        if (members.empty()) {
            ImGui::TextDisabled(st.selected_layer != NCI_LAYER_PREVIEW ? "Select a layer."
                                : critical_points_enabled ? "Select critical points to preview their NCI surface."
                                                          : "Tick Enabled to search for critical points.");
            return;
        }
        draw_layer_inspector(st.selected_layer, members);
        draw_region_table(region_stats(members));
        for (NciSurface* n : members) draw_member_status(*n);
        std::vector<std::vector<std::string>> environments;
        for (const NciSurface* n : members) environments.insert(environments.end(), n->environments.begin(), n->environments.end());
        draw_environments(environments);
    }

    static void draw_levels(float* levels, int* num_levels, const char* tooltip) {
        prop::row("levels", tooltip);
        density_map_levels_editor(levels, num_levels);
    }

    void draw_reference_inspector() {
        NciReferenceLayer& r = st.reference;
        prop::begin("nci layer");
        ImGui::PushID("##reference inspector");
        prop::row("name");
        ImGui::TextUnformatted("Reference");
        prop::row("style", "Dominant: per point the densest category, more opaque at each deeper level.\n"
                           "Combined: the categories' colours mixed where they overlap, more opaque towards each core.");
        {
            static const char* names[] = {"Dominant", "Combined"};
            static const int values[] = {(int)NciReferenceStyle::Dominant, (int)NciReferenceStyle::Combined};
            int s = (int)r.style;
            if (density_map_segmented("##style", names, values, 2, &s, nullptr, true) && s != (int)r.style) {
                r.style = (NciReferenceStyle)s;
                // Combined's graded fill reads without isolines; Dominant's levels with them
                r.lines = r.style == NciReferenceStyle::Dominant;
            }
        }
        prop::row("opacity");
        ImGui::SliderFloat("##opacity", &r.opacity, 0.0f, 1.0f, "%.2f");
        prop::row("draw");
        ImGui::Checkbox("fill", &r.fill);
        ImGui::SameLine();
        ImGui::Checkbox("lines", &r.lines);
        if (r.lines) {
            prop::row("line width");
            ImGui::SliderFloat("##line_width", &r.line_width, 0.5f, 6.0f, "%.1f px");
        }
        draw_levels(r.levels, &r.num_levels, "The areas enclose these fractions of each category's volume (highest density regions)");
        prop::hint("A sample is in the reference inside its category's %.0f %% area, whatever levels are drawn",
                   100.0f * nci_surface::reference_inside_fraction);
        ImGui::PopID();
        prop::end();
    }

    void draw_layer_inspector(int key, const std::vector<NciSurface*>& members) {
        NciScatterLayer& l = scatter_style(key);
        prop::begin("nci layer");
        ImGui::PushID(key);
        prop::row("name");
        if (key > 0) {
            NciSurface& first = *members.front();
            if (ImGui::InputText("##name", first.name, sizeof(first.name)))
                for (size_t i = 1; i < members.size(); ++i) snprintf(members[i]->name, sizeof(members[i]->name), "%s", first.name);
        } else ImGui::TextUnformatted(layer_name(key, members));
        prop::row("types");
        ImGui::TextWrapped("%s", nci_layer_types(members).c_str());
        prop::row("style", "Contours: the union of its interactions' regions at the levels.\n"
                           "Heatmap: the density of its samples.\n"
                           "The layer shown in the viewport also plots its points, on the surface colour scale.");
        {
            static const char* names[] = {"Contours", "Heatmap"};
            static const int values[] = {(int)NciLayerStyle::Contours, (int)NciLayerStyle::Heatmap};
            int s = (int)l.style;
            if (density_map_segmented("##style", names, values, 2, &s, nullptr, true)) l.style = (NciLayerStyle)s;
        }
        if (l.style == NciLayerStyle::Contours) {
            prop::row("color");
            ImVec4 colour = layer_colour(key, members);
            if (prop::color("##color", &colour.x)) {
                l.colour = colour;
                l.custom_colour = true;
                if (key > 0) for (NciSurface* n : members) n->chip = colour;
            }
        }
        prop::row("opacity");
        ImGui::SliderFloat("##opacity", &l.opacity, 0.0f, 1.0f, "%.2f");
        if (l.style == NciLayerStyle::Contours) {
            prop::row("draw");
            ImGui::Checkbox("fill", &l.fill);
            ImGui::SameLine();
            ImGui::Checkbox("lines", &l.lines);
            if (l.lines) {
                prop::row("line width");
                ImGui::SliderFloat("##line_width", &l.line_width, 0.5f, 10.0f, "%.1f px");
            }
            draw_levels(l.levels, &l.num_levels,
                        "The contours enclose these fractions of each interaction's volume, their union for the layer");
        }
        ImGui::PopID();
        prop::end();
    }

    // ---- what the scatter plot draws of each layer, and how its samples sit in the data-derived regions
    bool is_shown(int key) const { return key == st.selected_layer && !st.layer_hidden; }

    NciScatterLayer scatter_layer(int key) const {
        const auto it = st.scatter_layers.find(key);
        return it == st.scatter_layers.end() ? NciScatterLayer{} : it->second;
    }
    bool scatter_visible(int key) const { return scatter_layer(key).visible; }
    NciScatterLayer& scatter_style(int key) { return st.scatter_layers[key]; }
    ImVec4 layer_colour(int key, const std::vector<NciSurface*>& members) const {
        const NciScatterLayer l = scatter_layer(key);
        return l.custom_colour ? l.colour : layer_chip(key, members);
    }

    std::vector<int> layer_keys() const {
        std::vector<int> keys = {NCI_LAYER_FULL_SYSTEM, NCI_LAYER_PREVIEW};
        keys.insert(keys.end(), st.layer_order.begin(), st.layer_order.end());
        return keys;
    }

    void update_critical_point_list_scope(int key) {
        critical_point_list_scope.clear();
        for (const NciSurface* n : nci_layer_members(st, key))
            critical_point_list_scope.insert(critical_point_list_scope.end(), n->critical_point_indices.begin(),
                                             n->critical_point_indices.end());
        std::sort(critical_point_list_scope.begin(), critical_point_list_scope.end());
        critical_point_list_scope.erase(std::unique(critical_point_list_scope.begin(), critical_point_list_scope.end()),
                                        critical_point_list_scope.end());
    }

    void select_layer(int key) {
        nci_select_layer(st, key);
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
        std::array<nci_surface::RegionCounts, 4> by_type = {};
        size_t inside = 0, outside = 0;
    };
    static size_t reference_index(nci_surface::Family family) { return nci_surface::reference_category(family); }

    static inline const char* const reference_labels[nci_surface::reference_category_count] = {
        nci_surface::reference_category_name(0), nci_surface::reference_category_name(1),
        nci_surface::reference_category_name(2), nci_surface::reference_category_name(3)};
    static constexpr nci_surface::Family reference_colours[nci_surface::reference_category_count] = {
        nci_surface::Family::HydrogenBond, nci_surface::Family::Ring, nci_surface::Family::PiStacking,
        nci_surface::Family::OtherContact};

    // Each member's category samples (extracted with the reference settings) against their category's reference
    static RegionStats region_stats(const std::vector<NciSurface*>& members) {
        const nci_surface::ReferenceMap& map = nci_surface::builtin_reference_map();
        RegionStats stats;
        for (const NciSurface* n : members) {
            if (!n->have_key || n->cancelled || n->failed || n->reference_failed) continue;
            for (size_t family = 0; family < n->reference_scatters.size(); ++family) {
                const size_t group = reference_index((nci_surface::Family)family);
                if (group >= stats.by_type.size()) continue;
                const auto c = nci_surface::reference_counts(n->reference_scatters[family], map, group);
                stats.by_type[group].inside += c.inside;
                stats.by_type[group].outside += c.outside;
                stats.inside += c.inside;
                stats.outside += c.outside;
            }
        }
        return stats;
    }
    static void draw_region_table(const RegionStats& stats) {
        if (stats.inside + stats.outside == 0) return;
        if (!ImGui::BeginTable("##regions", 2, ImGuiTableFlags_SizingFixedFit)) return;
        ImGui::TableSetupColumn("");
        ImGui::TableSetupColumn("in reference (90 %)");
        ImGui::TableHeadersRow();
        for (size_t f = 0; f < stats.by_type.size(); ++f) {
            const auto& c = stats.by_type[f];
            if (c.inside + c.outside == 0) continue;
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextColored(nci_type_colour(reference_colours[f]), "%s", reference_labels[f]);
            ImGui::TableNextColumn();
            ImGui::Text("%.0f %%", 100.0 * (double)c.inside / (double)(c.inside + c.outside));
        }
        ImGui::EndTable();
    }

    // Isolines and area in the scatter plot, of a layer or an interaction
    static bool draw_style_controls(bool* lines, NciScatterArea* area, const char* area_tooltip) {
        bool changed = false;
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

    // An interaction's isolines and area, the same in every layer; or those of all the interactions at once
    void draw_type_style_controls(nci_surface::Family family, const char* name) {
        static const char* area_tooltip = "In each layer, in the interaction's colour. Filled: the inside of its outline.\n"
                                          "Heatmap: the density of its samples.";
        NciTypeStyle& t = st.type_styles[(size_t)family];
        ImGui::TextUnformatted(name);
        ImGui::PushID("type");
        draw_style_controls(&t.lines, &t.area, area_tooltip);
        ImGui::PopID();
        ImGui::Separator();
        ImGui::TextUnformatted("All interactions");
        NciTypeStyle all = t;
        ImGui::PushID("all");
        if (draw_style_controls(&all.lines, &all.area, area_tooltip))
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

    // The levels as the core expects them: sorted, distinct, within (0, 1). While one is dragged they may not be.
    static std::vector<float> clean_levels(const float* levels, int num_levels) {
        std::vector<float> out;
        for (int i = 0; i < CLAMP(num_levels, 0, nci_max_levels); ++i)
            if (levels[i] > 0.0f && levels[i] < 1.0f) out.push_back(levels[i]);
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    }
    // Per level, as the Ramachandran fills stacked from the outermost: n of them reach 0.8 at the core
    static float level_alpha(int levels) { return 1.0f - powf(0.2f, 1.0f / (float)MAX(1, levels)); }

    // The reference at the inspector's levels, rebuilt when they change
    const nci_surface::ReferenceMap& reference_map() {
        const std::vector<float> levels = clean_levels(st.reference.levels, st.reference.num_levels);
        uint64_t key = md_hash64(levels.data(), levels.size() * sizeof(float), 0) | 1;
        if (key != st.reference_map_key) {
            st.reference_map_key = key;
            st.reference_texture_key = 0;
            if (levels.empty() || !nci_surface::build_reference_map(&st.reference_map, nci_surface::reference_data(), levels))
                st.reference_map = {};
        }
        return st.reference_map;
    }

    // The reference fill as a texture, made per levels and style. Dominant: per pixel the densest category's colour,
    // more opaque at each level it is inside; Combined: the categories' colours mixed by their depth, more opaque with
    // depth and overlap. Empty pixels take a neighbour's colour at zero alpha so linear filtering leaves no dark fringe.
    void draw_reference_fill(const nci_surface::ReferenceMap& map) {
        constexpr int size = 384;
        constexpr float combined_alpha = 0.8f;
        constexpr size_t n_cat = nci_surface::reference_category_count;
        const uint64_t key = md_hash64(&st.reference.style, sizeof(st.reference.style), st.reference_map_key) | 1;
        if (st.reference_texture && st.reference_texture_key != key) gl::free_texture(&st.reference_texture);
        if (!st.reference_texture) {
            nci_surface::ReferenceImage image;
            nci_surface::reference_image(&image, map, size);
            const float a = level_alpha((int)map.fractions.size());
            std::vector<ImVec4> colours((size_t)size * size, ImVec4(0, 0, 0, 0));
            std::vector<bool> filled((size_t)size * size, false);
            for (size_t i = 0; i < colours.size(); ++i) {
                if (st.reference.style == NciReferenceStyle::Dominant) {
                    if (image.category[i] < 0) continue;
                    colours[i] = nci_type_colour(reference_colours[(size_t)image.category[i]]);
                    colours[i].w = 1.0f - powf(1.0f - a, (float)image.level[i]);
                } else {
                    ImVec4 sum(0, 0, 0, 0);
                    float weight = 0.0f, clear = 1.0f;
                    for (size_t c = 0; c < n_cat; ++c) {
                        const float d = image.depth[i * n_cat + c];
                        if (d <= 0.0f) continue;
                        const ImVec4 col = nci_type_colour(reference_colours[c]);
                        sum.x += d * col.x, sum.y += d * col.y, sum.z += d * col.z;
                        weight += d;
                        clear *= 1.0f - d;
                    }
                    if (weight <= 0.0f) continue;
                    colours[i] = ImVec4(sum.x / weight, sum.y / weight, sum.z / weight, combined_alpha * (1.0f - clear));
                }
                filled[i] = true;
            }
            for (int pass = 0; pass < 2; ++pass) {
                const std::vector<bool> from = filled;
                for (int j = 0; j < size; ++j) for (int i = 0; i < size; ++i) {
                    const size_t at = i + (size_t)size * j;
                    if (from[at]) continue;
                    for (int dy = -1; dy <= 1 && !filled[at]; ++dy) for (int dx = -1; dx <= 1; ++dx) {
                        const int x = i + dx, y = j + dy;
                        if (x < 0 || y < 0 || x >= size || y >= size || !from[x + (size_t)size * y]) continue;
                        colours[at] = colours[x + (size_t)size * y];
                        colours[at].w = 0.0f;
                        filled[at] = true;
                        break;
                    }
                }
            }
            std::vector<uint32_t> pixels(colours.size());
            for (size_t i = 0; i < pixels.size(); ++i) pixels[i] = ImGui::ColorConvertFloat4ToU32(colours[i]);
            if (!gl::init_texture_2D(&st.reference_texture, size, size, GL_RGBA8) ||
                !gl::set_texture_2D_data(st.reference_texture, 0, pixels.data(), GL_RGBA8)) {
                gl::free_texture(&st.reference_texture);
                return;
            }
            st.reference_texture_key = key;
        }
        // Image row 0 is s = 0, the bottom; ImPlot puts uv (0, 0) at the top left
        ImPlot::PlotImage("##NCI reference", (ImTextureID)(intptr_t)st.reference_texture, ImPlotPoint(-map.rho_extent, 0.0),
                          ImPlotPoint(map.rho_extent, nci_surface::scatter_s_max), ImVec2(0, 1), ImVec2(1, 0),
                          ImVec4(1, 1, 1, st.reference.opacity));
    }

    static void draw_triangles(const std::vector<std::array<float, 2>>& triangles, ImVec4 colour) {
        if (triangles.size() < 3) return;
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

    // An interaction's area: the inside of its outline in its colour, or the density of its samples
    void draw_type_area(ImVec4 colour, const NciSurface& n, NciScatterArea area) {
        if (!(n.rho_max > 0.0f)) return;
        const nci_surface::Scatter* scatter = &n.scatter;
        if (area == NciScatterArea::Filled) {
            colour.w = 0.25f;
            fill_area({&scatter, 1}, n.rho_max, colour, 0);
        } else if (area == NciScatterArea::Heatmap) heatmap_area({&scatter, 1}, n.rho_max, 0.8f);
    }

    // ---- the scatter plot: the expected regions, the full-density context and the layers, Ramachandran-like
    void draw_scatter() {
        static bool fit_view = false;
        NciSurface& context = st.context;
        ImGui::Checkbox("Guides", &st.show_guides);
        ImGui::SetItemTooltip("The isovalue and the colour range. Drag the isovalue line to edit it for every interaction.");
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
            const double asked = nci_resolution_spacing(context.resolution);
            const double used = MAX(context.grid.spacing[0], MAX(context.grid.spacing[1], context.grid.spacing[2]));
            const bool coarsened = used > 1.05 * asked;
            const nci_surface::ScatterDensity& density = context.scatter.density;
            ImGui::BeginGroup();
            ImGui::TextDisabled("%zu points, %.3f Bohr", context.scatter.s.size(), used);
            if (coarsened) {
                ImGui::SameLine(0.0f, 0.0f);
                ImGui::TextColored(ImVec4(0.95f, 0.65f, 0.2f, 1.0f), " (coarsened)");
            }
            ImGui::EndGroup();
            ImGui::SetItemTooltip(
                "Grid spacing %.3f Bohr (%s asks for %.3f%s).\n"
                "Points: %zu of %zu eligible grid samples, a capped subsample for display only.\n"
                "Outlines, areas, heatmaps and expected regions use the continuous density: %zu grid cells\n"
                "supersampled (%zu subsamples), %.1f Bohr^3 of eligible volume. Not an interaction-energy measure.",
                used, resolution_str[CLAMP(context.resolution, 0, (int)NciSurface::Res_Count - 1)], asked,
                coarsened ? "; the box would take too many samples" : "",
                context.scatter.s.size(), context.scatter.eligible_count, density.cells, density.subsamples,
                density.total_volume);
        }
        ImGui::SameLine();
        ImGui::BeginDisabled(st.export_scatter_in_progress);
        if (ImGui::SmallButton("Export NCI samples...")) st.export_scatter_requested = true;
        ImGui::EndDisabled();
        ImGui::SetItemTooltip("Export the shown layer's scatter and separate recomputed scatter clouds for each CP category, "
                              "with CP membership in a companion -anchors.csv file and the continuous density per\n"
                              "(sign(lambda2)rho, s) bin, in Bohr^3, in a companion -density.csv file.");
        if (st.export_scatter_in_progress) {
            ImGui::SameLine();
            ImGui::TextDisabled("Categorizing samples %.0f%%", 100.0f * st.export_scatter_progress);
        } else if (!st.export_scatter_status.empty()) {
            ImGui::SameLine();
            ImGui::TextDisabled("%s", st.export_scatter_status.c_str());
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
                ImGui::SliderFloat("Interaction line width", &st.type_outline_width, 0.5f, 8.0f, "%.1f px");
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10.0f);
                ImGui::SliderFloat("Point size", &st.point_size, 0.5f, 6.0f, "%.1f px");
                ImGui::Checkbox("Interaction isolines and areas", &st.show_type_styles);
                ImGui::SetItemTooltip("Per CP type, the isolines and areas of the layers' interactions, with a legend chip each.\n"
                                      "Off by default: the reference already shows the categories.");
                if (ImGui::Button("Reset")) {
                    st.type_outline_width = NciState::default_type_outline_width;
                    st.point_size = NciState::default_point_size;
                }
                ImGui::Separator();
                ImGui::TextDisabled("Legend: click a swatch to show or hide it, an interaction's name to make it\n"
                                    "the active one for the isovalue guide. How a layer or the reference is drawn\n"
                                    "is set in its inspector, below the layers.");
                ImGui::EndPopup();
            }
        }

        // The legend is a row of chips under the plot, so it never covers the data: the expected regions and the layers,
        // then the interactions
        const std::vector<NciSurface*> members = nci_layer_members(st, st.selected_layer);
        const ImVec4 guide_colour = layer_colour(st.selected_layer, members);
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
        legend.push_back({"Reference", ImVec4(0.85f, 0.85f, 0.85f, 1.0f), LegendEntry::Reference, 0,
                          nci_surface::Family::Count, nullptr, st.reference.visible, !st.reference.visible});
        if (surface_ready(context))
            legend.push_back({"All samples", ImVec4(0.8f, 0.8f, 0.8f, 1.0f), LegendEntry::Context, 0,
                              nci_surface::Family::Count, nullptr, st.show_context_samples, !st.show_context_samples});
        std::array<bool, nci_surface::family_count> present = {};
        std::vector<NciSurface*> shown;
        for (int key : keys) {
            const auto layer = nci_layer_members(st, key);
            if (layer.empty()) continue;
            LegendEntry e = {layer_name(key, layer), layer_colour(key, layer), LegendEntry::Layer, key};
            e.on = scatter_visible(key);
            e.dim = !e.on;
            legend.push_back(e);
            if (is_shown(key)) shown = plotted_members(key);
            if (!e.dim)
                for (const NciSurface* n : plotted_members(key)) if (n->families) present[(size_t)n->group] = true;
        }
        bool first_type = true;
        for (size_t f = 0; f < present.size(); ++f) {
            if (!present[f] || !st.show_type_styles) continue;
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
        const auto swatch_w = [&](const LegendEntry&) { return chip_h; };
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
        const nci_surface::ReferenceMap& reference = reference_map();
        float rho_extent = context.rho_max;
        if (reference.valid) rho_extent = MAX(rho_extent, reference.rho_extent);
        ImPlot::SetupAxisLimits(ImAxis_X1, -rho_extent, rho_extent, ImPlotCond_Once);
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

        // The reference is independent of the loaded system, the displayed layers and their legend visibility.
        if (st.reference.visible && reference.valid) {
            if (st.reference.fill) draw_reference_fill(reference);
            if (st.reference.lines) {
                for (size_t c = 0; c < reference.categories.size(); ++c) {
                    if (!reference.categories[c].valid) continue;
                    ImVec4 colour = nci_type_colour(reference_colours[c]);
                    colour.w = st.reference.opacity;
                    for (const auto& isoline : reference.categories[c].isolines) draw_outline(isoline, colour, st.reference.line_width);
                }
            }
        }

        NciSurface* hovered = nullptr;
        size_t hovered_sample = 0;
        float closest = 36.0f;
        const auto layer_outline = [&](int key) -> const NciState::LayerOutline* {
            const auto it = st.layer_outlines.find(key);
            return it == st.layer_outlines.end() ? nullptr : &it->second;
        };
        // A layer as its inspector sets: its contours (filled per level, stacked) or its heatmap; then, thin on top,
        // each interaction's isolines and area when shown
        const auto draw_layer = [&](int key) {
            const NciScatterLayer l = scatter_layer(key);
            if (!l.visible) return;
            const auto layer = plotted_members(key);
            const auto members_all = nci_layer_members(st, key);
            ImVec4 colour = layer_colour(key, members_all);
            const NciState::LayerOutline* outline = l.style == NciLayerStyle::Contours ? layer_outline(key) : nullptr;
            // The full system's outline also covers the saved layers, so it shows with none of its own plotted
            if (layer.empty() && key != NCI_LAYER_FULL_SYSTEM) return;
            if (outline) {
                if (l.fill) {
                    ImVec4 fill = colour;
                    fill.w = level_alpha((int)outline->fills.size()) * l.opacity;
                    for (const auto& triangles : outline->fills) draw_triangles(triangles, fill);
                }
                if (l.lines) {
                    ImVec4 line = colour;
                    line.w = l.opacity;
                    for (const auto& segments : outline->segments) draw_outline(segments, line, l.line_width);
                }
            } else if (l.style == NciLayerStyle::Heatmap && !layer.empty()) {
                float rho_max = 0.0f;
                std::vector<const nci_surface::Scatter*> scatters;
                for (const NciSurface* n : layer) {
                    rho_max = MAX(rho_max, n->rho_max);
                    scatters.push_back(&n->scatter);
                }
                heatmap_area(scatters, rho_max, 0.8f * l.opacity);
            }
            if (!st.show_type_styles) return;
            for (const NciSurface* n : layer) {
                if (!n->families) continue;
                const NciTypeStyle& t = st.type_styles[(size_t)n->group];
                if (!t.visible) continue;
                if (t.area != NciScatterArea::None) draw_type_area(type_colour(*n), *n, t.area);
                if (t.lines) {
                    ImVec4 c = type_colour(*n);
                    c.w = 1.0f;
                    draw_outline(n->outline, c, st.type_outline_width);
                }
            }
        };
        // The layer shown in the viewport on top
        for (int key : keys) if (key != st.selected_layer) draw_layer(key);
        draw_layer(st.selected_layer);
        // Point colours always come from the layer actually shown in the viewport, never from background layers.
        if (!st.layer_hidden) for (NciSurface* n : plotted_members(st.selected_layer)) {
            plot(*n, false);
            if (!ImPlot::IsPlotHovered()) continue;
            const ImVec2 mouse = ImGui::GetMousePos();
            for (size_t i = 0; i < n->scatter.s.size(); ++i) {
                const ImVec2 p = ImPlot::PlotToPixels(n->scatter.signed_rho[i], n->scatter.s[i]);
                const float d2 = (p.x - mouse.x) * (p.x - mouse.x) + (p.y - mouse.y) * (p.y - mouse.y);
                if (d2 < closest) { closest = d2; hovered = n; hovered_sample = i; }
            }
        }

        NciSurface* active = st.layer_hidden ? nullptr : active_member();
        bool guide_hovered = false;
        if (st.show_guides && active && nci_scale_member(*active)) {
            const ColorScale scale = nci_effective_scale(st, *active);
            const ImVec4 col = guide_colour;
            double iso = active->s_iso;
            bool held = false, hovered_line = false;
            if (ImPlot::DragLineY(1, &iso, col, 2.0f, 0, nullptr, &hovered_line, &held)) {
                st.uniform.s_iso = (float)std::clamp(iso, 0.1, 1.5);
            }
            st.settings_editing |= held;
            guide_hovered = hovered_line || held;
            const double limits[] = {scale.range_beg, scale.range_end};
            ImPlot::SetNextLineStyle(col, 1.0f);
            ImPlot::PlotInfLines("##Colour limits", limits, 2);
            const double right = ImPlot::GetPlotLimits().X.Max;
            ImPlot::Annotation(right, active->s_iso, col, ImVec2(-4, -4), true, "s=%.2f", active->s_iso);
            ImPlot::TagX(scale.range_beg, col, "%.4g", scale.range_beg);
            ImPlot::TagX(scale.range_end, col, "%.4g", scale.range_end);
        }
        if (hovered && !guide_hovered) {
            const int key = nci_layer_key(*hovered);
            const float rho = hovered->scatter.signed_rho[hovered_sample], s = hovered->scatter.s[hovered_sample];
            ImGui::BeginTooltip();
            ImGui::Text("%s | %s", layer_name(key, nci_layer_members(st, key)), nci_family_label(*hovered).c_str());
            ImGui::Text("sign(lambda2)rho %.5g au | s %.3f", rho, s);
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
            // The swatch toggles: a square for an interaction, a contour for expected regions, a thick bar for a layer
            const ImVec2 p = ImGui::GetCursorScreenPos();
            const bool toggled = ImGui::InvisibleButton("##swatch", ImVec2(swatch_w(e), chip_h));
            ImDrawList* dl = ImGui::GetWindowDrawList();
            const ImU32 col = ImGui::ColorConvertFloat4ToU32(colour);
            if (e.kind == LegendEntry::Layer) {
                const NciScatterLayer l = scatter_layer(e.layer);
                density_map_swatch(dl, p, ImVec2(p.x + chip_h, p.y + chip_h),
                                   l.style == NciLayerStyle::Heatmap ? DensityMapSwatch_Heatmap : DensityMapSwatch_Contours,
                                   e.colour, l.fill, ImPlotColormap_Viridis, e.on ? 1.0f : 0.35f);
                ImGui::SetItemTooltip("Click: %s the layer in the scatter plot\nPoints always follow the viewport layer.",
                                      e.on ? "hide" : "show");
            } else if (e.kind == LegendEntry::Reference) {
                reference_swatch(dl, p, ImVec2(p.x + chip_h, p.y + chip_h), e.on ? 1.0f : 0.35f);
                ImGui::SetItemTooltip("Click: %s the reference of hydrogen bonds, rings, pi stacking and other CPs",
                                      e.on ? "hide" : "show");
            } else if (e.kind == LegendEntry::Context) {
                dl->AddCircleFilled(ImVec2(p.x + 0.5f * chip_h, p.y + 0.5f * chip_h), 0.3f * chip_h, col);
                ImGui::SetItemTooltip("%s every grid sample of the whole density, under the layers", e.on ? "Hide" : "Show");
            } else {
                dl->AddRectFilled(p, ImVec2(p.x + chip_h, p.y + chip_h), col, style.FrameRounding);
                ImGui::SetItemTooltip("Click: %s this interaction's isolines and area in every layer\n"
                                      "Right-click: its isolines and area, or those of all the interactions",
                                      e.on ? "hide" : "show");
                if (ImGui::BeginPopupContextItem("##type scatter")) {
                    draw_type_style_controls(e.family, e.text.c_str());
                    ImGui::EndPopup();
                }
            }
            if (toggled) {
                if (e.kind == LegendEntry::Type) st.type_styles[(size_t)e.family].visible = !e.on;
                else if (e.kind == LegendEntry::Reference) st.reference.visible = !e.on;
                else if (e.kind == LegendEntry::Context) st.show_context_samples = !e.on;
                else scatter_style(e.layer).visible = !e.on;
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

void draw_nci_panel(NciState& state, bool critical_points_enabled, std::vector<uint32_t>& critical_point_list_scope,
                    std::vector<uint32_t>& critical_point_environment_selection) {
    NciPanel panel(state, critical_points_enabled, critical_point_list_scope, critical_point_environment_selection);
    const topology::NciPanelDrawCallbacks callbacks = {
        &panel,
        [](void* data) { static_cast<NciPanel*>(data)->draw_scatter(); },
        [](void* data) { static_cast<NciPanel*>(data)->draw_layers(); },
        [](void* data) { static_cast<NciPanel*>(data)->st.selection_key = 0; },
    };
    draw_nci_panel_layout(&state.enabled, callbacks);
}

} // namespace topology
