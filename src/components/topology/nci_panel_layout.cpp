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
#include <cmath>
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
    void (*draw_settings)(void*) = nullptr;
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
    ImGui::SameLine();
    if (ImGui::Button(ICON_FA_GEAR " NCI settings")) ImGui::OpenPopup("##NCI settings");
    if (ImGui::BeginPopup("##NCI settings")) {
        if (callbacks.draw_settings) callbacks.draw_settings(callbacks.user_data);
        ImGui::EndPopup();
    }

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
    static constexpr size_t max_samples = 32u << 20;
    static constexpr ImVec4 warning_colour = ImVec4(0.95f, 0.65f, 0.2f, 1.0f);

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

    // Prefer the deepest shared environment; mixed selections fall back to their branches.
    std::string selection_layer_name(const std::vector<NciSurface*>& members, int number) const {
        std::vector<std::string> families;
        std::vector<std::vector<std::string>> paths;
        for (const NciSurface* member : members) {
            for (const auto& path : member->environments) {
                if (path.empty()) continue;
                paths.push_back(path);
                const bool branched = (path[0] == "Contacts" || path[0] == "Rings") && path.size() > 1;
                const std::string& family = branched ? path[1] : path[0];
                if (std::find(families.begin(), families.end(), family) == families.end()) families.push_back(family);
            }
        }
        if (families.empty()) return "Layer " + std::to_string(number);
        std::string name = critical_point_distribution::selection_environment_name(paths);
        if (name.empty() || (families.size() > 1 && (name == "Contacts" || name == "Rings"))) {
            name = families.front();
            if (families.size() == 2) name += " + " + families[1];
            else if (families.size() > 2) name += " + " + std::to_string(families.size() - 1) + " more";
        }
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
        n->display_s = source.display_s;
        n->rho_max = source.rho_max;
        n->scale = source.scale;
        n->families = source.families;
        n->type_counts = source.type_counts;
        n->critical_point_indices = source.critical_point_indices;
        n->critical_points_by_family = source.critical_points_by_family;
        n->competitors = source.competitors;
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

    // Shared controls; layer-specific display stays in the layer inspector.
    void draw_settings() {
        bool editing = false;
        ImGui::SeparatorText("Surface");
        prop::begin("NCI settings");
        prop::row("surface isovalue", "The viewport surface. Scatter samples above it are drawn desaturated: in the region, "
                                      "not in the surface. Changing it does not change the scatter samples.");
        ImGui::SetNextItemWidth(prop::control_width());
        ImGui::SliderFloat("##NCI uniform s", &st.uniform.s_iso, 0.1f, st.uniform.scatter_s, "%.2f");
        editing |= ImGui::IsItemActive();
        prop::row("resolution", "Grid spacing for density sampling. Finer grids are smoother and take longer.");
        ImGui::SetNextItemWidth(prop::control_width());
        ImGui::Combo("##NCI resolution", &st.resolution, resolution_str, (int)ARRAY_SIZE(resolution_str));
        ImGui::SetItemTooltip("Coarse 0.3, Medium 0.2, Fine 0.13 Bohr; made coarser where that would take more than %zu M samples.",
                              max_samples >> 20);
        prop::row("opacity");
        ImGui::SetNextItemWidth(prop::control_width());
        ImGui::SliderFloat("##NCI opacity", &st.opacity, 0.05f, 1.0f, "%.2f");
        prop::end();
        ImGui::SeparatorText("Scatter / extraction");
        prop::begin("NCI scatter");
        prop::row("plot s range", "Display only: the scatter plot's s axis. Points above s max are shown in grey.");
        ImGui::SetNextItemWidth(prop::control_width());
        ImGui::SliderFloat("##NCI plot s", &st.uniform.plot_s, 0.2f, nci_surface::scatter_s_max, "%.2f");
        editing |= ImGui::IsItemActive();
        if (prop::section("Scatter region", false)) {
            // Tighter extraction cutoffs select a subset; the reference itself stays unchanged.
            const auto warn = [](bool beyond, const char* text) {
                if (!beyond) return;
                prop::row("");
                ImGui::PushStyleColor(ImGuiCol_Text, warning_colour);
                ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + prop::control_width());
                ImGui::TextWrapped(ICON_FA_TRIANGLE_EXCLAMATION " %s", text);
                ImGui::PopTextWrapPos();
                ImGui::PopStyleColor();
            };
            prop::row("s max", "Samples, connected regions, contours and density up to this reduced density gradient; above "
                               "it, up to the plot s range, points are still shown in grey, for display only. Also bounds the surface.");
            ImGui::SetNextItemWidth(prop::control_width());
            ImGui::SliderFloat("##NCI scatter s", &st.uniform.scatter_s, 0.2f, nci_surface::scatter_s_max, "%.2f");
            editing |= ImGui::IsItemActive();
            st.uniform.s_iso = MIN(st.uniform.s_iso, st.uniform.scatter_s);
            warn(st.uniform.scatter_s > nci_surface::extended_reference_data().s_max + 1e-4,
                 "Beyond the reference (s <= 2): no reference data above it.");
            prop::row("\xCF\x81 max", "Samples denser than this are out, both signs of lambda2. Also bounds the surface.");
            ImGui::SetNextItemWidth(prop::control_width());
            ImGui::SliderFloat("##NCI scatter rho", &st.uniform.rho_max, 0.01f, 0.15f, "%.3f au");
            editing |= ImGui::IsItemActive();
            warn(st.uniform.rho_max > nci_surface::extended_reference_data().rho_extent + 1e-6,
                 "Beyond the reference (|rho| <= 0.1): no reference data outside it.");
            prop::row("R", "Distance from the CPs the regions may reach. Also bounds the surface.");
            ImGui::SetNextItemWidth(prop::control_width());
            ImGui::SliderFloat("##NCI scatter R", &st.uniform.radius, 1.0f, 5.0f, "%.1f Bohr");
            editing |= ImGui::IsItemActive();
            warn(std::abs(st.uniform.radius - nci_surface::reference_radius) > 1e-4,
                 "The reference used R = 2.7 Bohr: points and contours are not aligned with it.");
            prop::row("");
            if (ImGui::Button("Reset extraction defaults")) {
                const NciUniformSettings defaults = nci_uniform_default_settings();
                st.uniform.scatter_s = defaults.scatter_s;
                st.uniform.rho_max = defaults.rho_max;
                st.uniform.radius = defaults.radius;
                st.uniform.s_iso = MIN(st.uniform.s_iso, st.uniform.scatter_s);
            }
            prop::hint("Scores use the selected extraction region, excluding grey display-only points.");
            prop::section_end();
        }
        prop::hint("Uniform connected regions, both signs of lambda2, nearest-CP ownership.");
        prop::end();
        ImGui::SeparatorText("Display");
        prop::begin("NCI display");
        prop::row("interaction line width");
        ImGui::SetNextItemWidth(prop::control_width());
        ImGui::SliderFloat("##NCI line width", &st.type_outline_width, 0.5f, 8.0f, "%.1f px");
        prop::row("point size");
        ImGui::SetNextItemWidth(prop::control_width());
        ImGui::SliderFloat("##NCI point size", &st.point_size, 0.5f, 6.0f, "%.1f px");
        prop::row("interaction isolines and areas");
        ImGui::Checkbox("##NCI type styles", &st.show_type_styles);
        ImGui::SetItemTooltip("Per CP type, the layers' interaction isolines and areas. Off by default: "
                              "the reference already shows the categories.");
        color_scale_draw_controls(&st.scale, st.span, md_unit_pow(md_unit_bohr_radius(), -3), "shown layer",
                                  "the 1st to 99th percentile of the shown layer", true);
        prop::end();
        ImGui::SeparatorText("Plot guides");
        ImGui::Checkbox("Surface isovalue", &st.guide_surface);
        ImGui::Checkbox("Extraction cutoff", &st.guide_cutoff);
        ImGui::Checkbox("Colour-range boundaries", &st.guide_colour_limits);
        ImGui::TextDisabled("Guides beside the plot shows/hides these without changing your choices.");

        if (ImGui::Button("Reset to defaults")) {
            st.uniform = nci_uniform_default_settings();
            st.type_outline_width = NciState::default_type_outline_width;
            st.point_size = NciState::default_point_size;
            st.show_type_styles = false;
            st.guide_surface = st.guide_cutoff = st.guide_colour_limits = true;
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
        bool hovered = false, clicked = false, context_clicked = false, remove = false;
    };
    // The row's selectable and its eye; the caller draws the swatch and text into the frame
    RowFrame begin_row(bool selected, bool* visible, const char* eye_tooltip, const LayerRow* saved = nullptr) {
        RowFrame r;
        r.line_h = ImGui::GetTextLineHeight();
        r.h = r.line_h * 2.0f + 8.0f;
        r.min = ImGui::GetCursorScreenPos();
        r.w = ImGui::GetContentRegionAvail().x;
        ImGui::SetNextItemAllowOverlap();
        r.clicked = ImGui::Selectable("##row", selected, ImGuiSelectableFlags_AllowOverlap, ImVec2(r.w, r.h));
        r.hovered = ImGui::IsItemHovered(ImGuiHoveredFlags_DelayNormal | ImGuiHoveredFlags_NoSharedDelay) &&
                    !ImGui::IsMouseDragging(ImGuiMouseButton_Left);
        r.context_clicked = ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenOverlappedByItem) &&
                            ImGui::IsMouseReleased(ImGuiMouseButton_Right);
        if (saved) layer_tree_drag_source({NCI_LAYER_LIST_ID, (uint32_t)saved->key, false}, saved->label);
        const float pad = 4.0f;
        float controls_x = r.min.x + pad;
        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(2, ImGui::GetStyle().FramePadding.y));
        if (saved) {
            ImGui::SetCursorScreenPos(ImVec2(controls_x, r.min.y + (r.h - ImGui::GetFrameHeight()) * 0.5f));
            r.remove = ImGui::Button(ICON_FA_TRASH_CAN "##delete");
            ImGui::SetItemTooltip("Delete layer");
            controls_x += ImGui::GetItemRectSize().x + 4.0f;
            if (ImGui::IsItemHovered() || ImGui::IsItemActive()) r.clicked = r.hovered = false;
        }
        ImGui::SetCursorScreenPos(ImVec2(controls_x, r.min.y + (r.h - ImGui::GetFrameHeight()) * 0.5f));
        if (ImGui::Button(*visible ? ICON_FA_EYE "##vis" : ICON_FA_EYE_SLASH "##vis")) *visible = !*visible;
        if (ImGui::IsItemHovered() || ImGui::IsItemActive()) r.clicked = r.hovered = false;
        ImGui::PopStyleVar();
        ImGui::PopStyleColor();
        ImGui::SetItemTooltip("%s", eye_tooltip);
        r.swatch_min = ImVec2(controls_x + r.line_h + 8, r.min.y + (r.h - r.line_h) * 0.5f);
        r.swatch_max = ImVec2(r.swatch_min.x + r.line_h, r.swatch_min.y + r.line_h);
        r.text = ImVec2(controls_x + r.line_h + 6 + r.line_h + 6, r.min.y + 4);
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
                                     style.visible ? "Hide it in the scatter plot" : "Show it in the scatter plot",
                                     saved ? &row : nullptr);
        if (r.remove) st.pending_delete = row.key;
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
        if (saved && !r.remove) {
            const ImRect rect(r.min, ImVec2(r.min.x + r.w, r.min.y + r.h));
            if (layer_tree_drop_target(rect, r.min.x, {NCI_LAYER_LIST_ID, (uint32_t)row.key, false},
                                       false, nodes, node_count, move_src, move_drop))
                dropped = true;
        }
        if (r.context_clicked) ImGui::OpenPopup("##layer menu");
        if (ImGui::BeginPopup("##layer menu")) {
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
        const nci_surface::ReferenceMap& map = display_reference_map();
        int calculations = 0;
        for (const auto& c : map.categories) calculations += c.valid ? c.systems : 0;
        char detail[128];
        snprintf(detail, sizeof(detail), "s <= 2, |rho| <= 0.1 \xC2\xB7 H-bonds, rings, \xCF\x80, other");
        row_text(r, "Reference", detail, st.reference.visible, selected);
        if (r.hovered) {
            ImGui::BeginTooltip();
            ImGui::TextUnformatted("Reference");
            ImGui::TextDisabled("Built in, independent of the loaded system, as the Ramachandran reference: per CP category,\n"
                                "the pooled (sign(lambda2)rho, s) density of reference extractions (s <= %.1f, rho <= %.2f,\n"
                                "%.1f Bohr from a CP, each sample to its nearest CP) over many calculations, each weighing\n"
                                "the same. Scores compare the selected extraction region against its %.0f %% area.",
                                2.0,
                                display_reference_data().rho_extent, nci_surface::reference_radius,
                                100.0f * nci_surface::reference_inside_fraction);
            if (*display_reference_data().source) ImGui::TextDisabled("%s", display_reference_data().source);
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
        const float button = ImGui::GetFrameHeight();
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + MAX(0.0f, ImGui::GetContentRegionAvail().x - button));
        ImGui::BeginDisabled(!preview_ready());
        if (ImGui::Button(ICON_FA_PLUS "##save", ImVec2(button, 0))) save_preview();
        ImGui::EndDisabled();
        ImGui::SetItemTooltip("Save the preview (the selected CPs) as a layer.\nSelect CPs and wait for the preview to complete.");
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
        prop::row("domain");
        ImGui::TextUnformatted("s <= 2, |rho| <= 0.1, R = 2.7 Bohr");
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
        prop::hint("Scores: fraction of extracted samples inside this reference's %.0f %% area",
                   100.0f * nci_surface::reference_inside_fraction);
        ImGui::PopID();
        prop::end();
        draw_reference_subfamilies();
    }

    static const nci_surface::ContactClass* guide_contact_class(const char* name) {
        for (const auto& c : nci_surface::contact_classes()) {
            if (critical_point_distribution::contact_environment_name(c.name) == name) return &c;
        }
        return nullptr;
    }

    static const char* guide_description(const nci_surface::ReferenceCategory& sub) {
        const char* branch = sub.parent[0] ? sub.parent : sub.name;
        if (const auto* c = guide_contact_class(branch)) return c->description;
        const std::string_view name = branch;
        if (name == "Aromatic rings") return "Ring critical points assigned to aromatic rings.";
        if (name == "Aliphatic rings") return "Ring critical points assigned to non-aromatic rings.";
        if (name == "Sugar rings") return "Ring critical points assigned to sugar rings.";
        if (name == "Contact-closed rings") return "Cycles closed by noncovalent contacts, such as hydrogen-bonded rings.";
        return "Cage critical points enclosed by several ring faces.";
    }

    void draw_reference_subfamilies() {
        const nci_surface::ReferenceMap& map = display_reference_map();
        if (map.subfamilies.empty()) return;
        st.reference.guides.resize(map.subfamilies.size(), 0);
        if (!ImGui::CollapsingHeader("Interaction guides")) return;
        ImGui::PushID("##subfamilies");
        if (ImGui::Button("Add interaction guide...")) {
            st.reference.guide_search[0] = '\0';
            ImGui::OpenPopup("##guide catalogue");
        }
        if (std::any_of(st.reference.guides.begin(), st.reference.guides.end(), [](uint8_t v) { return v != 0; })) {
            ImGui::SameLine();
            if (ImGui::Button("Clear")) {
                std::fill(st.reference.guides.begin(), st.reference.guides.end(), 0);
            }
        }
        ImGui::SetNextWindowSize(ImVec2(440, 460), ImGuiCond_Appearing);
        if (ImGui::BeginPopup("##guide catalogue")) {
            ImGui::InputTextWithHint("##search", "Find a name or chemistry...", st.reference.guide_search,
                                     sizeof(st.reference.guide_search));
            ImGuiTextFilter filter(st.reference.guide_search);
            bool found = false;
            ImGui::BeginChild("##catalogue", ImVec2(0, 0));
            const auto matches = [&](const auto& sub) {
                return filter.PassFilter(sub.name) || filter.PassFilter(guide_description(sub));
            };
            const auto entry = [&](size_t i, const char* label) {
                const auto& sub = map.subfamilies[i];
                ImGui::PushID((int)i);
                if (ImGui::Selectable(label, st.reference.guides[i] != 0)) {
                    st.reference.guides[i] = 2;
                    ImGui::CloseCurrentPopup();
                }
                ImGui::SetItemTooltip("%s\n%d supporting systems. Draws the 90 %% volume edge.",
                                      guide_description(sub), sub.systems);
                ImGui::PopID();
            };
            for (size_t c = 0; c < map.categories.size(); ++c) {
                bool heading = false;
                for (size_t i = 0; i < map.subfamilies.size(); ++i) {
                    const auto& sub = map.subfamilies[i];
                    if (sub.category != c || !sub.valid || sub.parent[0]) continue;
                    bool match = matches(sub), has_children = false;
                    for (const auto& child : map.subfamilies) {
                        if (child.category != c || std::string_view(child.parent) != sub.name) continue;
                        has_children = true;
                        match |= matches(child);
                    }
                    if (!match) continue;
                    if (!heading) {
                        ImGui::SeparatorText(reference_labels[c]);
                        heading = true;
                    }
                    found = true;
                    if (!has_children) {
                        entry(i, sub.name);
                        continue;
                    }
                    ImGui::PushID((int)i);
                    if (filter.IsActive()) ImGui::SetNextItemOpen(true, ImGuiCond_Always);
                    const bool open = ImGui::TreeNode(sub.name);
                    ImGui::SetItemTooltip("%s\nExpand to choose the whole class or a finer interaction.",
                                          guide_description(sub));
                    if (open) {
                        entry(i, "Whole class");
                        for (size_t j = 0; j < map.subfamilies.size(); ++j) {
                            const auto& child = map.subfamilies[j];
                            if (child.category == c && std::string_view(child.parent) == sub.name && matches(child))
                                entry(j, child.name);
                        }
                        ImGui::TreePop();
                    }
                    ImGui::PopID();
                }
            }
            if (!found) ImGui::TextDisabled("No available maps match.");
            if (ImGui::TreeNode("Other recognised interactions")) {
                ImGui::TextWrapped("Recognised by CP classification, but without enough reference data for a standalone map.");
                for (const auto& c : nci_surface::contact_classes()) {
                    bool available = false;
                    for (const auto& sub : map.subfamilies)
                        available |= sub.valid && guide_contact_class(sub.name) == &c;
                    if (available || (!filter.PassFilter(c.name) && !filter.PassFilter(c.description))) continue;
                    ImGui::TextDisabled("%s", c.name);
                    ImGui::SetItemTooltip("%s\nNo standalone map: fewer than 5 supporting systems in this reference.",
                                          c.description);
                }
                ImGui::TreePop();
            }
            ImGui::EndChild();
            ImGui::EndPopup();
        }
        for (size_t i = 0; i < map.subfamilies.size(); ++i) {
            if (!st.reference.guides[i]) continue;
            const auto& sub = map.subfamilies[i];
            ImGui::PushID((int)i);
            bool on = st.reference.guides[i] == 2;
            if (ImGui::Checkbox("##on", &on))
                st.reference.guides[i] = on ? 2 : 1;
            ImGui::SameLine();
            ImGui::TextColored(subfamily_colour(map, i), "%s", sub.name);
            ImGui::SetItemTooltip("%s\n%d supporting systems; 90 %% volume edge.\nIndependent of the category map's visibility and levels.",
                                  guide_description(sub), sub.systems);
            ImGui::SameLine();
            if (ImGui::SmallButton("x")) {
                st.reference.guides[i] = 0;
            }
            ImGui::PopID();
        }
        prop::hint("Guides show only their 90 %% isoline. Hide Reference to view guides alone.");
        ImGui::PopID();
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

    // A subfamily's colour: a shade of its category's, apart from its siblings'
    static ImVec4 subfamily_colour(const nci_surface::ReferenceMap& map, size_t subfamily) {
        const nci_surface::ReferenceCategory& sub = map.subfamilies[subfamily];
        size_t index = 0, count = 0;
        for (size_t i = 0; i < map.subfamilies.size(); ++i) {
            if (map.subfamilies[i].category != sub.category) continue;
            if (i < subfamily) ++index;
            ++count;
        }
        const ImVec4 base = nci_type_colour(reference_colours[sub.category]);
        float h, s, v;
        ImGui::ColorConvertRGBtoHSV(base.x, base.y, base.z, h, s, v);
        const float t = count > 1 ? (float)index / (float)(count - 1) : 0.5f;
        h = std::fmod(h + 0.14f * (t - 0.5f) + 1.0f, 1.0f);
        v = std::clamp(v * (1.15f - 0.45f * t), 0.0f, 1.0f);
        ImVec4 out = base;
        ImGui::ColorConvertHSVtoRGB(h, s, v, out.x, out.y, out.z);
        return out;
    }

    // Whether (rho, s) lies inside some category's outermost drawn reference level.
    bool reference_inside(float rho, float s) {
        const nci_surface::ReferenceMap& map = reference_map();
        for (size_t c = 0; c < nci_surface::reference_category_count; ++c)
            if (nci_surface::reference_level(map, c, rho, s) > 0) return true;
        return false;
    }

    // The reference at (rho, s) for a tooltip: its categories from the most likely, each with its share and the
    // smallest drawn area enclosing the point
    void draw_reference_here(float rho, float s) {
        const nci_surface::ReferenceMap& map = reference_map();
        const auto shares = nci_surface::reference_shares(map, rho, s);
        std::array<size_t, nci_surface::reference_category_count> order;
        for (size_t c = 0; c < order.size(); ++c) order[c] = c;
        std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return shares[a] > shares[b]; });
        ImGui::TextDisabled("Reference here");
        for (size_t c : order) {
            if (shares[c] < 0.005f) continue;
            ImGui::TextColored(nci_type_colour(reference_colours[c]), "%s", reference_labels[c]);
            ImGui::SameLine();
            ImGui::Text("%.0f %%", 100.0f * shares[c]);
            const int level = nci_surface::reference_level(map, c, rho, s);
            if (level > 0) {
                ImGui::SameLine();
                ImGui::TextDisabled("inside its %.0f %% area", 100.0f * map.fractions[map.fractions.size() - level]);
            }
            // Its subfamilies, by their share of the category's
            const std::vector<float> sub = nci_surface::reference_subfamily_shares(map, c, rho, s);
            std::vector<size_t> subs;
            for (size_t i = 0; i < sub.size(); ++i) if (shares[c] * sub[i] >= 0.005f) subs.push_back(i);
            std::stable_sort(subs.begin(), subs.end(), [&](size_t a, size_t b) { return sub[a] > sub[b]; });
            ImGui::Indent();
            for (size_t i : subs) {
                ImGui::TextColored(subfamily_colour(map, i), "%s", map.subfamilies[i].name);
                ImGui::SameLine();
                ImGui::TextDisabled("%.0f %%", 100.0f * shares[c] * sub[i]);
            }
            ImGui::Unindent();
        }
        if (!map.subfamilies.empty())
            ImGui::TextDisabled("Subfamilies split their category's share, each weighing the same");
    }

    // Each member's category samples from the selected extraction region against the reference.
    static RegionStats region_stats(const std::vector<NciSurface*>& members) {
        const nci_surface::ReferenceMap& map = nci_surface::extended_reference_map();
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

    const nci_surface::ReferenceData& display_reference_data() const {
        return nci_surface::extended_reference_data();
    }

    const nci_surface::ReferenceMap& display_reference_map() const {
        return nci_surface::extended_reference_map();
    }

    // The selected display reference at the inspector's levels.
    const nci_surface::ReferenceMap& reference_map() {
        const std::vector<float> levels = clean_levels(st.reference.levels, st.reference.num_levels);
        uint64_t key = md_hash64(levels.data(), levels.size() * sizeof(float), 0) | 1;
        if (key != st.reference_map_key) {
            st.reference_map_key = key;
            st.reference_texture_key = 0;
            if (levels.empty() || !nci_surface::build_reference_map(&st.reference_map, display_reference_data(), levels))
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
        static float shown_plot_s = -1.0f;
        NciSurface& context = st.context;
        ImGui::Checkbox("Guides", &st.show_guides);
        ImGui::SetItemTooltip("Show selected plot guides: surface isovalue, extraction cutoff and colour boundaries. "
                              "Choose guides in NCI settings. Drag the surface line to edit its isovalue.");
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
                ImGui::TextColored(warning_colour, " (coarsened)");
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
        // Fit remains beside the plot.
        {
            const ImGuiStyle& style = ImGui::GetStyle();
            const float bw = ImGui::CalcTextSize(ICON_FA_EXPAND).x + style.FramePadding.x * 2;
            ImGui::SameLine();
            ImGui::SetCursorPosX(ImGui::GetCursorPosX() + MAX(0.0f, ImGui::GetContentRegionAvail().x - bw));
            if (ImGui::Button(ICON_FA_EXPAND "##fit")) fit_view = true;
            ImGui::SetItemTooltip("Fit the view (or double-click the plot)");
        }

        // The legend is a row of chips under the plot, so it never covers the data: the expected regions and the layers,
        // then the interactions
        const std::vector<NciSurface*> members = nci_layer_members(st, st.selected_layer);
        const ImVec4 guide_colour = layer_colour(st.selected_layer, members);
        const std::vector<int> keys = layer_keys();
        struct LegendEntry {
            std::string text;
            ImVec4 colour;
            enum Kind { Reference, ReferenceGuide, Context, Layer, Type } kind = Layer;
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
        const auto& guide_map = display_reference_map();
        st.reference.guides.resize(guide_map.subfamilies.size(), 0);
        for (size_t i = 0; i < guide_map.subfamilies.size(); ++i) {
            if (!st.reference.guides[i]) continue;
            const bool on = st.reference.guides[i] == 2;
            legend.push_back({std::string("Ref: ") + guide_map.subfamilies[i].name, subfamily_colour(guide_map, i),
                              LegendEntry::ReferenceGuide, (int)i, nci_surface::Family::Count, nullptr, on, !on});
        }
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
        const bool reset_view = fit_view || shown_plot_s != st.uniform.plot_s;
        shown_plot_s = st.uniform.plot_s;
        if (reset_view) ImPlot::SetNextAxisToFit(ImAxis_X1);
        fit_view = false;
        const bool plot_open = ImPlot::BeginPlot("##NCI scatter", ImVec2(-1, MAX(1.0f, ImGui::GetContentRegionAvail().y - legend_h)),
                                                 ImPlotFlags_NoLegend);
        if (!plot_open) return;
        ImPlot::SetupAxes("sign(\xCE\xBB\xE2\x82\x82)\xCF\x81 [bohr^-3]", "s");
        const nci_surface::ReferenceMap& reference = reference_map();
        float rho_extent = context.rho_max;
        if (reference.valid) rho_extent = MAX(rho_extent, reference.rho_extent);
        ImPlot::SetupAxisLimits(ImAxis_X1, -rho_extent, rho_extent, ImPlotCond_Once);
        ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, st.uniform.plot_s, reset_view ? ImPlotCond_Always : ImPlotCond_Once);
        ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, 0.0, nci_surface::scatter_s_max);
        ImPlot::SetupFinish();
        if (ImPlot::IsPlotHovered() && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) fit_view = true;

        struct Points {
            const std::vector<float>* rho;
            const std::vector<float>* s;
            const std::vector<uint32_t>* indices;
        };
        const auto plot_bin = [&](const NciSurface& n, size_t bin, bool tail, const std::vector<uint32_t>& indices, ImVec4 colour,
                                  float size) {
            if (indices.empty()) return;
            const auto getter = [](int i, void* data) {
                const auto& p = *(const Points*)data;
                const uint32_t idx = (*p.indices)[i];
                return ImPlotPoint((*p.rho)[idx], (*p.s)[idx]);
            };
            Points points = tail ? Points{&n.scatter.tail.signed_rho, &n.scatter.tail.s, &indices}
                                 : Points{&n.scatter.signed_rho, &n.scatter.s, &indices};
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, size, colour, 0.0f, colour);
            char label[64];
            snprintf(label, sizeof(label), "##NCI %d bin %zu", n.id, bin);
            ImPlot::PlotScatterG(label, getter, &points, (int)indices.size());
        };
        // The samples of a surface, on the colour scale of its surface in the viewport; the context is neutral unless
        // the global scale applies. Below the surface isovalue: full colour, in the surface. Up to the scatter's s:
        // desaturated, in the analysed region but not in the surface. The display-only tail above it: plain grey.
        // Drawn tail first so the surface's own samples stay on top.
        constexpr size_t colour_bins = 65;
        const auto plot = [&](const NciSurface& n, bool neutral) {
            const auto& scatter = n.scatter;
            const ColorScale scale = n.id == -1 ? nci_effective_scale(st, st.primary) : nci_effective_scale(st, n);
            const auto colour_bin = [&](float rho) {
                return neutral ? 0 : MIN(64, (int)(64.0f * color_scale_param(scale, rho) + 0.5f));
            };
            // Tail, above the isovalue, then the surface's samples
            std::array<std::vector<uint32_t>, 3 * colour_bins> bins;
            for (size_t i = 0; i < scatter.tail.s.size(); ++i)
                bins[colour_bin(scatter.tail.signed_rho[i])].push_back((uint32_t)i);
            for (size_t i = 0; i < scatter.s.size(); ++i) {
                const bool above = scatter.s[i] > n.s_iso;
                bins[(above ? 1 : 2) * colour_bins + colour_bin(scatter.signed_rho[i])].push_back((uint32_t)i);
            }
            for (size_t b = 0; b < bins.size(); ++b) {
                const bool tail = b < colour_bins;
                const bool faded = b < 2 * colour_bins;
                const size_t c = b % colour_bins;
                ImVec4 colour = neutral ? ImVec4(0.8f, 0.8f, 0.8f, tail ? 0.08f : faded ? 0.18f : 0.35f)
                                        : ImPlot::SampleColormap((float)c / (float)(colour_bins - 1), scale.colormap);
                if (!neutral) {
                    colour.w = 0.95f;
                    if (tail) {
                        colour = ImVec4(0.5f, 0.5f, 0.5f, 0.16f);
                    } else if (faded) {
                        constexpr float keep = 0.3f;
                        colour = ImVec4(0.85f + keep * (colour.x - 0.85f), 0.85f + keep * (colour.y - 0.85f),
                                        0.85f + keep * (colour.z - 0.85f), 0.55f);
                    }
                }
                plot_bin(n, b, tail, bins[b], colour, st.point_size);
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

        // The reference domain is independent of extraction limits.
        ImPlot::PushPlotClipRect();
        {
            const float rho = reference.rho_extent;
            const float s = nci_surface::scatter_s_max;
            const ImVec2 lo = ImPlot::PlotToPixels(-rho, 0.0);
            const ImVec2 hi = ImPlot::PlotToPixels(rho, s);
            ImGui::PushClipRect(ImVec2(MIN(lo.x, hi.x), MIN(lo.y, hi.y)), ImVec2(MAX(lo.x, hi.x), MAX(lo.y, hi.y) + 1.0f), true);
        }
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
        for (size_t i = 0; i < guide_map.subfamilies.size(); ++i) {
            if (st.reference.guides[i] != 2 || !guide_map.subfamilies[i].valid) continue;
            const auto& sub = guide_map.subfamilies[i];
            if (sub.isolines.size() < 2) continue;
            draw_outline(sub.isolines[1], subfamily_colour(guide_map, i), std::max(1.5f, st.reference.line_width));
        }
        ImGui::PopClipRect();
        ImPlot::PopPlotClipRect();

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
        const bool surface_guide = st.show_guides && st.guide_surface && active && nci_scale_member(*active);
        const bool merged_guide = surface_guide && st.guide_cutoff &&
                                  std::abs(active->s_iso - st.uniform.scatter_s) < 1e-4f;
        if (st.show_guides && st.guide_cutoff && !merged_guide) {
            const ImPlotRect limits = ImPlot::GetPlotLimits();
            const double cutoff = st.uniform.scatter_s;
            const ImVec4 colour(0.75f, 0.75f, 0.75f, 0.85f);
            const ImVec2 left = ImPlot::PlotToPixels(limits.X.Min, cutoff);
            const ImVec2 right = ImPlot::PlotToPixels(limits.X.Max, cutoff);
            ImPlot::PushPlotClipRect();
            for (float x = left.x; x < right.x; x += 10.0f)
                ImPlot::GetPlotDrawList()->AddLine(ImVec2(x, left.y), ImVec2(std::min(x + 5.0f, right.x), right.y),
                                                  ImGui::ColorConvertFloat4ToU32(colour), 1.0f);
            ImPlot::PopPlotClipRect();
            ImPlot::Annotation(limits.X.Max, cutoff, colour, ImVec2(-4, -4), true, "Cutoff s = %.3g", cutoff);
        }
        if (st.show_guides && active && nci_scale_member(*active)) {
            const ColorScale scale = nci_effective_scale(st, *active);
            const ImVec4 col = guide_colour;
            double iso = active->s_iso;
            bool held = false, hovered_line = false;
            if (st.guide_surface &&
                ImPlot::DragLineY(1, &iso, col, 2.0f, ImPlotDragToolFlags_NoFit, nullptr, &hovered_line, &held)) {
                st.uniform.s_iso = (float)std::clamp(iso, 0.1, (double)st.uniform.scatter_s);
            }
            st.settings_editing |= held;
            guide_hovered = hovered_line || held;
            const double right = ImPlot::GetPlotLimits().X.Max;
            if (st.guide_surface)
                ImPlot::Annotation(right, st.uniform.s_iso, col, ImVec2(-4, -4), true,
                                   merged_guide ? "Surface / cutoff s = %.2f" : "Surface s = %.2f", st.uniform.s_iso);
            if (st.guide_colour_limits) {
                const double limits[] = {scale.range_beg, scale.range_end};
                ImPlot::SetNextLineStyle(col, 1.0f);
                ImPlot::PlotInfLines("##Colour limits", limits, 2);
                ImPlot::TagX(scale.range_beg, col, "%.4g", scale.range_beg);
                ImPlot::TagX(scale.range_end, col, "%.4g", scale.range_end);
            }
        }
        if (hovered && !guide_hovered) {
            const int key = nci_layer_key(*hovered);
            const float rho = hovered->scatter.signed_rho[hovered_sample], s = hovered->scatter.s[hovered_sample];
            ImGui::BeginTooltip();
            ImGui::Text("%s | %s", layer_name(key, nci_layer_members(st, key)), nci_family_label(*hovered).c_str());
            ImGui::Text("sign(lambda2)rho %.5g au | s %.3f", rho, s);
            if (s > hovered->s_iso)
                ImGui::TextDisabled("Above the surface isovalue (s = %.2f): in the region, not in the surface", hovered->s_iso);
            const bool guided = st.show_guides && active == hovered;
            if (is_shown(key))
                ImGui::TextDisabled(guided ? "Click to hide the guides" : "Click to make this interaction the active one for the guides");
            if (st.reference.visible && reference_inside(rho, s)) {
                ImGui::Separator();
                draw_reference_here(rho, s);
            }
            ImGui::EndTooltip();
            if (is_shown(key) && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) toggle_guide(*hovered);
        } else if (st.reference.visible && !guide_hovered && ImPlot::IsPlotHovered()) {
            const ImPlotPoint mouse = ImPlot::GetPlotMousePos();
            if (reference_inside((float)mouse.x, (float)mouse.y)) {
                ImGui::BeginTooltip();
                ImGui::TextDisabled("sign(lambda2)rho %.5g au | s %.3f", mouse.x, mouse.y);
                draw_reference_here((float)mouse.x, (float)mouse.y);
                ImGui::EndTooltip();
            }
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
            } else if (e.kind == LegendEntry::ReferenceGuide) {
                dl->AddLine(ImVec2(p.x, p.y + chip_h * 0.5f),
                            ImVec2(p.x + chip_h, p.y + chip_h * 0.5f), col, 2.0f);
                const auto& sub = guide_map.subfamilies[(size_t)e.layer];
                ImGui::SetItemTooltip("%s\n%d supporting systems; 90 %% reference isoline.\nClick to %s this guide.",
                                      guide_description(sub), sub.systems, e.on ? "hide" : "show");
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
                else if (e.kind == LegendEntry::ReferenceGuide) st.reference.guides[(size_t)e.layer] = e.on ? 1 : 2;
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
        [](void* data) { static_cast<NciPanel*>(data)->draw_settings(); },
        [](void* data) { static_cast<NciPanel*>(data)->st.selection_key = 0; },
    };
    draw_nci_panel_layout(&state.enabled, callbacks);
}

} // namespace topology
