#include <layer_groups.h>

#include <viamd.h>
#include <serialization_utils.h>
#include <property_grid.h>
#include <imgui_widgets.h>

#include <core/md_allocator.h>
#include <core/md_hash.h>

#include <imgui_internal.h>

#include <stdio.h>
#include <string.h>

// ## Groups

void layer_group_update(LayerGroup* g, const ApplicationState* app, uint64_t generation, bool follow_frame, bool want_sets) {
    ASSERT(g && app);
    md_allocator_i* heap = md_get_heap_allocator();

    uint64_t key = md_hash64_combine(generation, (g->follow_selection ? 1 : 2) | (want_sets ? 4 : 0));
    if (g->follow_selection) {
        key = md_hash64_combine(key, md_bitset_hash64(md_bitmask_view(&app->selection.selection_mask), 0));
    } else if (g->filter.text[0]) {
        key = md_hash64_combine(key, md_hash64(g->filter.text, strlen(g->filter.text), app->script.ir_fingerprint));
        key = md_hash64_combine(key, app->mold.sys.topology_version);
        if (follow_frame && g->filter.result.is_dynamic) {
            key = md_hash64_combine(key, (uint64_t)(app->animation.frame * 4096.0));
        }
    }
    if (key == g->key) return;
    g->key = key;

    for (size_t i = 0; i < md_array_size(g->sets); ++i) md_bitset_free(&g->sets[i], heap);
    md_array_shrink(g->sets, 0);
    g->valid = true;

    if (g->follow_selection) {
        md_bitset_assign(&g->atoms, md_bitmask_view(&app->selection.selection_mask), heap);
        if (want_sets && !md_bitset_empty(g->atoms)) md_array_push(g->sets, md_bitset_copy(g->atoms, heap), heap);
    } else if (g->filter.text[0]) {
        g->valid = filter_evaluate(&g->filter, &g->atoms, filter_context(app), heap);
        if (g->valid && want_sets) {
            g->valid = filter_evaluate_sets(&g->filter, &g->sets, filter_context(app), heap);
        }
    } else {
        md_bitset_free(&g->atoms, heap);
        g->atoms = {};
    }

    uint64_t h = md_hash64_combine(key, g->valid ? 0x600D : 0xBAD);
    h = md_bitset_hash64(g->atoms, h);
    for (size_t i = 0; i < md_array_size(g->sets); ++i) h = md_bitset_hash64(g->sets[i], h + i);
    g->hash = h;
}

void layer_group_free(LayerGroup* g) {
    md_allocator_i* heap = md_get_heap_allocator();
    for (size_t i = 0; i < md_array_size(g->sets); ++i) md_bitset_free(&g->sets[i], heap);
    md_array_free(g->sets, heap);
    md_bitset_free(&g->atoms, heap);
    g->sets = 0;
    g->atoms = {};
    g->key = 0;
    g->hash = 0;
}

md_filter_context_t layer_group_filter_context(const ApplicationState* app, const LayerGroup* g) {
    md_filter_context_t ctx = filter_context(app);
    if (g && layer_group_has_context(*g)) ctx.atoms = &g->atoms;
    return ctx;
}

void layer_group_label(char* buf, size_t cap, const LayerGroup& g) {
    if (g.name[0])                snprintf(buf, cap, "%s", g.name);
    else if (g.follow_selection)  snprintf(buf, cap, "Selection");
    else if (g.filter.text[0])    snprintf(buf, cap, "%s", g.filter.text);
    else                          snprintf(buf, cap, "Group");
}

bool layer_group_compose_filter(char* layer_filter, size_t cap, const LayerGroup& g) {
    if (g.follow_selection) return false;
    if (!g.filter.text[0]) return true;
    char out[sizeof(AtomFilter::text)];
    int n = 0;
    if (!layer_filter[0]) n = snprintf(out, sizeof(out), "%s", g.filter.text);
    // A context is all of the group's atoms as one: 'in' evaluates its left within each structure of its right, so the
    // right is made one (residue(2) in residue(1:5) would ask each of the five for a second residue)
    else                  n = snprintf(out, sizeof(out), "(%s) in union(%s)", layer_filter, g.filter.text);
    if (n < 0 || (size_t)n >= sizeof(out) || (size_t)n >= cap) return false;
    MEMCPY(layer_filter, out, (size_t)n + 1);
    return true;
}

// ## The tree

bool layer_tree_can_move(RepresentationNode* const nodes[], size_t count, size_t idx, int target, RepresentationDrop drop) {
    if (!representation_tree_can_move(nodes, count, idx, target, drop)) return false;
    if (!nodes[idx]->group || target < 0) return true;
    // A group stays at the top level: not into another, nor beside what is in one
    if (drop == RepresentationDrop_IntoFirst || drop == RepresentationDrop_IntoLast) return false;
    return nodes[target]->parent == 0;
}

bool layer_tree_move(RepresentationNode* nodes[], size_t count, size_t idx, int target, RepresentationDrop drop) {
    if (!layer_tree_can_move(nodes, count, idx, target, drop)) return false;
    return representation_tree_move(nodes, count, idx, target, drop);
}

int32_t layer_tree_first_order(RepresentationNode* const nodes[], size_t count, uint32_t parent) {
    if (parent == 0) return -1;     // Before everything: the orders start at 0
    const int g = representation_tree_find(nodes, count, true, parent);
    return g >= 0 ? nodes[g]->order + 1 : -1;
}

// ## UI

void layer_tree_drag_source(LayerNodeRef ref, const char* label) {
    if (ImGui::BeginDragDropSource(ImGuiDragDropFlags_None)) {
        ImGui::SetDragDropPayload(LAYER_NODE_DND, &ref, sizeof(ref));
        if (ref.group) ImGui::Text("%s %s", ICON_FA_FOLDER, label);
        else           ImGui::TextUnformatted(label);
        ImGui::EndDragDropSource();
    }
}

// Where on a row something dropped goes, by the height the mouse is at
static RepresentationDrop drop_zone(float t, bool group, bool open) {
    if (group) {
        if (t < 0.25f) return RepresentationDrop_Before;
        if (t > 0.75f) return open ? RepresentationDrop_IntoFirst : RepresentationDrop_After;
        return RepresentationDrop_IntoLast;
    }
    return t < 0.5f ? RepresentationDrop_Before : RepresentationDrop_After;
}

bool layer_tree_drop_target(const ImRect& row, float x0, LayerNodeRef target, bool open, RepresentationNode* const nodes[], size_t count,
                            LayerNodeRef* dropped, RepresentationDrop* drop_out) {
    if (!ImGui::BeginDragDropTargetCustom(row, ImGui::GetID("##drop"))) return false;
    bool delivered = false;
    const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(LAYER_NODE_DND, ImGuiDragDropFlags_AcceptBeforeDelivery | ImGuiDragDropFlags_AcceptNoDrawDefaultRect);
    if (payload && payload->DataSize == (int)sizeof(LayerNodeRef)) {
        LayerNodeRef src;
        MEMCPY(&src, payload->Data, sizeof(src));
        const int si = representation_tree_find(nodes, count, src.group, src.id);
        const int ti = representation_tree_find(nodes, count, target.group, target.id);
        const float t = (ImGui::GetMousePos().y - row.Min.y) / MAX(row.GetHeight(), 1.0f);
        const RepresentationDrop drop = drop_zone(t, target.group, open);
        if (src.window == target.window && si >= 0 && ti >= 0 && layer_tree_can_move(nodes, count, (size_t)si, ti, drop)) {
            ImDrawList* dl = ImGui::GetWindowDrawList();
            const ImU32 col = ImGui::GetColorU32(ImGuiCol_DragDropTarget);
            switch (drop) {
            case RepresentationDrop_Before:
                dl->AddLine({x0, row.Min.y}, {row.Max.x, row.Min.y}, col, 2.0f);
                break;
            case RepresentationDrop_After:
                dl->AddLine({x0, row.Max.y}, {row.Max.x, row.Max.y}, col, 2.0f);
                break;
            case RepresentationDrop_IntoFirst:
                dl->AddLine({x0 + ImGui::GetStyle().IndentSpacing, row.Max.y}, {row.Max.x, row.Max.y}, col, 2.0f);
                break;
            case RepresentationDrop_IntoLast:
                dl->AddRect(row.Min, row.Max, col, ImGui::GetStyle().FrameRounding, 0, 2.0f);
                break;
            }
            if (payload->IsDelivery()) {
                *dropped = src;
                *drop_out = drop;
                delivered = true;
            }
        }
    }
    ImGui::EndDragDropTarget();
    return delivered;
}

LayerGroupRow layer_group_row(LayerGroup& g, uint32_t window, bool selected, const char* summary, bool invalid) {
    LayerGroupRow r;
    ImGui::PushID("group");
    ImGui::PushID((int)g.node.id);

    const float h = ImGui::GetFrameHeight() + 4.0f;
    const ImVec2 row_min = ImGui::GetCursorScreenPos();
    const float row_w = ImGui::GetContentRegionAvail().x;
    r.rect = ImRect(row_min, ImVec2(row_min.x + row_w, row_min.y + h));
    r.x0 = row_min.x;

    char label[96];
    layer_group_label(label, sizeof(label), g);

    ImGui::SetNextItemAllowOverlap();
    if (ImGui::Selectable("##row", selected, ImGuiSelectableFlags_AllowOverlap | ImGuiSelectableFlags_AllowDoubleClick, ImVec2(row_w, h))) {
        r.clicked = true;
        if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) g.open = !g.open;
    }
    r.hovered = ImGui::IsItemHovered();
    layer_tree_drag_source({window, g.node.id, true}, label);
    if (ImGui::BeginPopupContextItem("##menu")) {
        r.clicked = true;
        if (ImGui::MenuItem(ICON_FA_COPY " Duplicate", nullptr, false)) r.action = LayerGroupAction_Duplicate;
        if (ImGui::MenuItem("Ungroup")) r.action = LayerGroupAction_Ungroup;
        ImGui::SetItemTooltip("Take the layers out of the group, each with a filter that keeps what it shows");
        ImGui::Separator();
        if (ImGui::MenuItem(ICON_FA_XMARK " Remove")) r.action = LayerGroupAction_Remove;
        ImGui::SetItemTooltip("Remove the group and its layers");
        ImGui::EndPopup();
    }

    const float pad = 4.0f;
    const float y_btn = row_min.y + (h - ImGui::GetFrameHeight()) * 0.5f;
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
    ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(2, ImGui::GetStyle().FramePadding.y));

    // The folder: open or not. In the accent colour it has a subset of its own, in red one that does not evaluate.
    ImGui::SetCursorScreenPos(ImVec2(row_min.x + pad, y_btn));
    const ImVec4 icon_col = invalid ? ImVec4(0.95f, 0.35f, 0.35f, 1.0f)
                          : layer_group_has_context(g) ? ImGui::GetStyleColorVec4(ImGuiCol_CheckMark)
                          : ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled);
    ImGui::PushStyleColor(ImGuiCol_Text, icon_col);
    if (ImGui::Button(g.open ? ICON_FA_FOLDER_OPEN "##open" : ICON_FA_FOLDER "##open")) g.open = !g.open;
    ImGui::PopStyleColor();
    ImGui::SetItemTooltip(g.open ? "Collapse" : "Expand");

    ImGui::SameLine(0, 2);
    if (ImGui::Button(g.visible ? ICON_FA_EYE "##vis" : ICON_FA_EYE_SLASH "##vis")) g.visible = !g.visible;
    ImGui::SetItemTooltip(g.visible ? "Hide the group" : "Show the group");
    ImGui::PopStyleVar();
    ImGui::PopStyleColor();

    // The name, and what it is of, to the right
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const float text_x = ImGui::GetItemRectMax().x + 6.0f;
    const float text_y = row_min.y + (h - ImGui::GetTextLineHeight()) * 0.5f;
    const float right = row_min.x + row_w - pad;
    const float sum_w = summary && summary[0] ? MIN(ImGui::CalcTextSize(summary).x, (right - text_x) * 0.45f) : 0.0f;
    const ImVec4 clip_name = ImVec4(text_x, row_min.y, right - sum_w - 8.0f, row_min.y + h);
    dl->AddText(ImGui::GetFont(), ImGui::GetFontSize(), ImVec2(text_x, text_y), ImGui::GetColorU32(g.visible ? ImGuiCol_Text : ImGuiCol_TextDisabled),
                label, nullptr, 0.0f, &clip_name);
    if (sum_w > 0.0f) {
        const ImVec4 clip_sum = ImVec4(right - sum_w, row_min.y, right, row_min.y + h);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize(), ImVec2(right - sum_w, text_y),
                    invalid ? IM_COL32(255, 110, 110, 255) : ImGui::GetColorU32(selected ? ImGuiCol_Text : ImGuiCol_TextDisabled, selected ? 0.7f : 1.0f),
                    summary, nullptr, 0.0f, &clip_sum);
    }

    ImGui::SetCursorScreenPos(ImVec2(row_min.x, row_min.y + h + 2));
    ImGui::Dummy(ImVec2(0, 0));
    ImGui::PopID();
    ImGui::PopID();
    return r;
}

bool layer_group_inspector(LayerGroup& g, const ApplicationState* app, const char* what) {
    bool changed = false;
    ImGui::PushID("group");
    ImGui::PushID((int)g.node.id);

    char hint[96];
    {
        char saved[sizeof(g.name)];
        MEMCPY(saved, g.name, sizeof(saved));
        g.name[0] = '\0';
        layer_group_label(hint, sizeof(hint), g);
        MEMCPY(g.name, saved, sizeof(saved));
    }
    prop::row("group");
    ImGui::InputTextWithHint("##name", hint, g.name, sizeof(g.name));

    char tip[256];
    snprintf(tip, sizeof(tip), "The %s of the layers in the group are within what this picks: their filters are evaluated within it, "
                               "as the left of 'in' in a script, and a layer without one shows all of it", what);
    prop::row("subset", tip);
    {
        const float btn_w = ImGui::GetFrameHeight() + 6;
        ImGui::SetNextItemWidth(MAX(1.0f, prop::control_width() - btn_w - 2));
        ImGui::BeginDisabled(g.follow_selection);
        const bool invalid = !g.follow_selection && !g.valid && g.filter.text[0];
        if (invalid) ImGui::PushInvalid();
        if (filter_input("##filter", &g.filter, filter_context(app), 0, nullptr, g.follow_selection ? "the selection" : "none: the group only groups")) {
            changed = true;
        }
        if (invalid) ImGui::PopInvalid();
        ImGui::EndDisabled();
        ImGui::SameLine(0, 2);
        if (g.follow_selection) ImGui::PushStyleColor(ImGuiCol_Button, ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));
        if (ImGui::Button(ICON_FA_ARROW_POINTER "##follow", ImVec2(btn_w, 0))) {
            g.follow_selection = !g.follow_selection;
            changed = true;
        }
        if (g.follow_selection) ImGui::PopStyleColor();
        ImGui::SetItemTooltip(g.follow_selection ? "Following the selection: click to use the filter" : "Follow the selection instead of the filter");
        if (!g.follow_selection && g.filter.text[0] && g.valid) {
            prop::hint("%llu atoms", (unsigned long long)md_bitset_count(g.atoms));
        }
        if (!g.follow_selection && g.filter.result.is_dynamic) {
            prop::hint("Depends on the coordinates: follows the frame only when every layer in the group is of the frame shown");
        }
    }

    ImGui::PopID();
    ImGui::PopID();
    return changed;
}

// ## Workspace

void layer_group_serialize(viamd::serialization_state_t& state, str_t section, const LayerGroup& g) {
    viamd::write_section_header(state, section);
    viamd::write_int(state, STR_LIT("Id"), (int64_t)g.node.id);
    if (g.name[0]) viamd::write_str(state, STR_LIT("Name"), str_from_cstr(g.name));
    viamd::write_bool(state, STR_LIT("Visible"), g.visible);
    viamd::write_bool(state, STR_LIT("Open"), g.open);
    if (g.filter.text[0]) viamd::write_str(state, STR_LIT("Filter"), str_from_cstr(g.filter.text));
    if (g.follow_selection) viamd::write_bool(state, STR_LIT("FollowSelection"), true);
}

void layer_group_deserialize_entry(LayerGroup* g, str_t ident, str_t arg, int* file_id) {
    if (str_eq_cstr(ident, "Id")) {
        viamd::extract_int(*file_id, arg);
    } else if (str_eq_cstr(ident, "Name")) {
        viamd::extract_to_char_buf(g->name, sizeof(g->name), arg);
    } else if (str_eq_cstr(ident, "Visible")) {
        viamd::extract_bool(g->visible, arg);
    } else if (str_eq_cstr(ident, "Open")) {
        viamd::extract_bool(g->open, arg);
    } else if (str_eq_cstr(ident, "Filter")) {
        viamd::extract_to_char_buf(g->filter.text, sizeof(g->filter.text), arg);
    } else if (str_eq_cstr(ident, "FollowSelection")) {
        viamd::extract_bool(g->follow_selection, arg);
    }
}
