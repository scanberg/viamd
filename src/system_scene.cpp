#include <system_scene.h>

#include <viamd.h>
#include <representation.h>
#include <filter_input.h>
#include <imgui_widgets.h>
#include <property_grid.h>
#include <serialization_utils.h>
#include <script_editor.h>

#include <core/md_bitset.h>
#include <core/md_str.h>

#include <imgui.h>
#include <imgui_internal.h>

#include <stdio.h>

// ---------------------------------------------------------------------------------------------------------
// The operations
// ---------------------------------------------------------------------------------------------------------

struct OperationInfo {
    SystemOperation op;
    const char* name;       // As the window lists it
    const char* key;        // As a workspace names it ([SystemScene] Operations)
    const char* icon;
    const char* tooltip;
    bool moves_atoms;       // Not applied while an assembly is shown
};

// In the order they apply in (apply_state_operations), which is the order they are listed in
static const OperationInfo operation_info[SYSTEM_OPERATION_COUNT] = {
    {SystemOperation_Center, "Center", "center", ICON_FA_CROSSHAIRS,
     "In every frame, move everything so the target is in the middle of the box (at the origin without a box),\n"
     "and keep its orientation if asked", true},
    {SystemOperation_Wrap, "Wrap into box", "wrap", ICON_FA_CUBE,
     "In every frame, move every atom that is outside the periodic box back in through the opposite side.\n"
     "Molecules across the boundary are split.", true},
    {SystemOperation_MakeWhole, "Make whole", "make_whole", ICON_FA_PUZZLE_PIECE,
     "In every frame, join the molecules that are split across the periodic box boundary,\n"
     "so that each one is drawn in one piece", true},
    {SystemOperation_TrackBonds, "Track bond changes", "track_bonds", ICON_FA_LINK,
     "In every frame, show the bonds that form and break against the bonds of the system.\n"
     "The bonds of the system stay as they are.", false},
};

static const OperationInfo& info_of(SystemOperation op) {
    for (const OperationInfo& info : operation_info) {
        if (info.op == op) return info;
    }
    ASSERT(false);
    return operation_info[0];
}

SystemOperation system_operation_at(int idx) {
    ASSERT(0 <= idx && idx < SYSTEM_OPERATION_COUNT);
    return operation_info[idx].op;
}

const char* system_operation_name(SystemOperation op) {
    return info_of(op).name;
}

// The flag that turns it on and off
static bool* enabled_flag(ApplicationState* app, SystemOperation op) {
    auto& o = app->system.operations;
    switch (op) {
    case SystemOperation_Center:     return &o.recenter;
    case SystemOperation_Wrap:       return &o.apply_pbc;
    case SystemOperation_MakeWhole:  return &o.unwrap_structures;
    case SystemOperation_TrackBonds: return &o.bond_tracking.enabled;
    }
    ASSERT(false);
    return &o.recenter;
}

bool system_operation_enabled(const ApplicationState& app, SystemOperation op) {
    const auto& o = app.system.operations;
    switch (op) {
    case SystemOperation_Center:     return o.recenter;
    case SystemOperation_Wrap:       return o.apply_pbc;
    case SystemOperation_MakeWhole:  return o.unwrap_structures;
    case SystemOperation_TrackBonds: return o.bond_tracking.enabled;
    }
    return false;
}

bool system_has_operation(const ApplicationState& app, SystemOperation op) {
    return (app.system.operations.present & op) || system_operation_enabled(app, op);
}

void system_apply_now(ApplicationState* app, bool center, bool wrap, bool make_whole) {
    ASSERT(app);
    // Its operators place the coordinates as the file gives them (system.assembly)
    if (app->system.assembly >= 0) return;
    if (apply_state_operations(app, center, wrap, make_whole)) {
        app->mold.dirty_gpu_buffers |= MolBit_DirtyPosition | MolBit_ClearVelocity;
    }
}

void system_set_operation_enabled(ApplicationState* app, SystemOperation op, bool enabled) {
    ASSERT(app);
    auto& o = app->system.operations;
    *enabled_flag(app, op) = enabled;
    if (enabled) o.present |= op;

    if (op == SystemOperation_TrackBonds) {
        o.bond_tracking.scope_dirty = true;
    } else if (enabled) {
        // The frame shown as every frame will be: all that is on, in its order, the one turned on with it
        system_apply_now(app, o.recenter, o.apply_pbc, o.unwrap_structures);
    }
}

void system_add_operation(ApplicationState* app, SystemOperation op) {
    system_set_operation_enabled(app, op, true);
}

void system_remove_operation(ApplicationState* app, SystemOperation op) {
    ASSERT(app);
    auto& o = app->system.operations;
    *enabled_flag(app, op) = false;
    o.present &= ~(uint32_t)op;

    switch (op) {
    case SystemOperation_Center:
        o.fixate_orientation = false;
        o.recenter_query.enabled = false;
        o.recenter_query.filter.text[0] = '\0';
        o.recenter_query.filter.result = {};
        md_bitset_free(&o.selection_mask, app->allocator.persistent);
        o.selection_mask = md_bitset_empty_set();
        recenter_mark_query_dirty(app);
        break;
    case SystemOperation_TrackBonds: {
        auto& bt = o.bond_tracking;
        bt.scope.text[0] = '\0';
        bt.scope.result = {};
        bt.scope_dirty = true;
        bt.tolerance = MD_BOND_DELTA_TOLERANCE;
        bt.breakable = MD_BOND_DELTA_BREAKABLE;
        break;
    }
    default:
        break;
    }
}

void system_operations_sync(ApplicationState* app) {
    ASSERT(app);
    for (const OperationInfo& info : operation_info) {
        if (system_operation_enabled(*app, info.op)) app->system.operations.present |= info.op;
    }
}

// ---------------------------------------------------------------------------------------------------------
// Copies and naming
// ---------------------------------------------------------------------------------------------------------

void system_set_assembly(ApplicationState* app, int assembly) {
    ASSERT(app);
    if (assembly == app->system.assembly) return;
    // Its copies at once (they are set every frame, update_view_param), and what is shown changes its extent: the
    // view frames it
    app->system.assembly = assembly;
    assembly_transforms(&app->view.assembly_transforms, &app->view.assembly_copy_offset, app->mold.sys, assembly, app->allocator.persistent);
    reset_view(&app->view.target, app->mold.state, &app->representation.visibility_mask, app->view.camera.fov_y, app->view.region, &app->mold.sys, assembly);
    // Room to zoom out to all of it
    app->view.trackball_param.max_distance = MAX(app->view.trackball_param.max_distance, app->view.target.distance * 10.0f);
}

void system_label(char* buf, size_t cap, const ApplicationState& app) {
    const FileEntry* file = file_list_system_file(&app.files.list);
    str_t name = {};
    if (file && extract_file(&name, str_from_cstr(file->path)) && !str_empty(name)) {
        snprintf(buf, cap, STR_FMT, STR_ARG(name));
    } else {
        snprintf(buf, cap, "System");
    }
}

// A chain of an assembly: a copy of a polymer strand, a part of the instances it copies
static size_t assembly_num_chains(const md_system_t& sys, size_t asm_idx) {
    const md_assembly_data_t& assemblies = sys.assembly;
    size_t n = 0;
    for (size_t g = assemblies.group_offset[asm_idx]; g < assemblies.group_offset[asm_idx + 1]; ++g) {
        const AssemblyGroup group = assembly_group(assemblies, g);
        size_t polymers = 0;
        for (const md_instance_idx_t* k = group.inst; k < group.inst_end; ++k) {
            const md_urange_t res = md_system_instance_residue_range(&sys, (size_t)*k);
            for (uint32_t r = res.beg; r < res.end;) {
                const md_part_idx_t part = md_system_part_find_by_residue_idx(&sys, r);
                if (part < 0) break;
                polymers += md_part_kind_is_polymer(md_system_part_kind(&sys, (size_t)part)) ? 1 : 0;
                r = md_system_part_residue_range(&sys, (size_t)part).end;
            }
        }
        n += polymers * (size_t)(group.oper_end - group.oper);
    }
    return n;
}

// ---------------------------------------------------------------------------------------------------------
// UI
// ---------------------------------------------------------------------------------------------------------

// The width of the column the icons are in, the same as the handles of the representations so that the names line up
static float icon_column_width() {
    return MAX(ImGui::CalcTextSize(ICON_FA_GRIP_VERTICAL).x,
           MAX(ImGui::CalcTextSize(ICON_FA_FOLDER).x, ImGui::CalcTextSize(ICON_FA_FOLDER_OPEN).x)) + 6.0f;
}

static float row_button_size() {
    return ImGui::GetFontSize() + 3.0f * 2;
}

// The representations were not updated while the system was hidden; brought up to what they show now
static void refresh_representations(ApplicationState* app) {
    for (size_t i = 0; i < md_array_size(app->representation.reps); ++i) {
        app->representation.reps[i].needs_update = true;
    }
    flag_dynamic_representations_as_dirty(app);
}

// The system's menu, of its header's right click and of its ⋮ button: the operations it can be given, then what is
// done once
static void draw_system_menu(ApplicationState* app) {
    const bool assembly_shown = app->system.assembly >= 0;
    const bool has_cell = md_unitcell_flags(&app->mold.state.unitcell) != 0;

    ImGui::TextDisabled("Add operation");
    bool any = false;
    for (const OperationInfo& info : operation_info) {
        if (system_has_operation(*app, info.op)) continue;
        any = true;
        char label[64];
        snprintf(label, sizeof(label), "%s %s", info.icon, info.name);
        if (ImGui::MenuItem(label)) {
            system_add_operation(app, info.op);
        }
        ImGui::SetItemTooltip("%s", info.tooltip);
    }
    if (!any) {
        ImGui::TextDisabled("  (it has all of them)");
    }

    ImGui::Separator();
    ImGui::TextDisabled("Now, to the frame shown");

    const md_bitset_t& target_mask = recenter_get_active_target_mask(app);
    const bool target = !md_bitset_empty(target_mask);
    if (ImGui::MenuItem(ICON_FA_CROSSHAIRS " Center target", nullptr, false, target && !assembly_shown)) {
        system_apply_now(app, true, false, false);
    }
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
        if (target) {
            md_bitmask_assign(&app->selection.highlight_mask, target_mask);
            ImGui::SetTooltip("Move everything so the target (highlighted) is in the middle of the box,\nor at the origin when there is no box");
        } else {
            ImGui::SetTooltip("There is nothing to center on yet: give the system the Center operation and a target");
        }
    }
    if (ImGui::MenuItem(ICON_FA_CUBE " Wrap into box", nullptr, false, has_cell && !assembly_shown)) {
        system_apply_now(app, false, true, false);
    }
    ImGui::SetItemTooltip("Move every atom that is outside the periodic box back in through the opposite side");
    if (ImGui::MenuItem(ICON_FA_PUZZLE_PIECE " Make whole", nullptr, false, !assembly_shown)) {
        system_apply_now(app, false, false, true);
    }
    ImGui::SetItemTooltip("Join molecules that are split across the periodic box boundary");
    if (ImGui::MenuItem(ICON_FA_ARROWS_ROTATE " Recompute bonds")) {
        // The whole frame nearest the animation time, or the coordinates shown when there is no run
        recompute_bonds(app, (int64_t)(app->animation.frame + 0.5));
    }
    ImGui::SetItemTooltip("Guess the covalent bonds again from the distances between atoms in this frame, for when\n"
                          "the file's bonds are wrong or a reaction should be the topology from here on.\n"
                          "Changes the bonds of the system: the bonds a file, a topology or you gave stay.");
    if (assembly_shown) {
        ImGui::TextDisabled("Centering, wrapping and making whole are not\napplied while an assembly is shown");
    }
}

static void draw_bond_tracking_details(ApplicationState* app);

// The settings of an operation, under its row
static void draw_operation_properties(ApplicationState* app, SystemOperation op) {
    auto& o = app->system.operations;
    prop::begin("system_operation");
    switch (op) {
    case SystemOperation_Center: {
        const md_bitset_t& target_mask = recenter_get_active_target_mask(app);

        prop::row("target", "What is centered: the atoms set from a selection, or those a query picks\n"
                            "(evaluated again in every frame when it depends on it)");
        bool by_query = o.recenter_query.enabled;
        if (ImGui::RadioButton("selection", !by_query)) by_query = false;
        ImGui::SameLine();
        if (ImGui::RadioButton("query", by_query)) by_query = true;
        if (by_query != o.recenter_query.enabled) {
            o.recenter_query.enabled = by_query;
            if (o.recenter) system_apply_now(app, true, o.apply_pbc, o.unwrap_structures);
        }

        if (o.recenter_query.enabled) {
            prop::row("query");
            if (filter_input("##recenter-query", &o.recenter_query.filter, app->script.language, filter_context(app))) {
                recenter_mark_query_dirty(app);
            }
        } else {
            prop::row("atoms");
            const size_t num_selected = md_bitmask_count(&app->selection.selection_mask);
            ImGui::BeginDisabled(num_selected == 0);
            if (ImGui::SmallButton("set from selection")) {
                md_bitset_assign(&o.selection_mask, md_bitmask_view(&app->selection.selection_mask), app->allocator.persistent);
                recenter_update_target_data(app);
                if (o.recenter) system_apply_now(app, true, o.apply_pbc, o.unwrap_structures);
            }
            ImGui::EndDisabled();
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
                ImGui::SetTooltip(num_selected ? "The atoms selected now, %zu of them" : "Select atoms first (%zu selected)", num_selected);
            }
            ImGui::SameLine();
            const size_t n = md_bitset_count(o.selection_mask);
            if (n) ImGui::Text("%zu atoms", n);
            else   ImGui::TextDisabled("none yet");
        }
        // The target where it is: highlighted while its rows are hovered
        if (ImGui::IsItemHovered() && !md_bitset_empty(target_mask)) {
            md_bitmask_assign(&app->selection.highlight_mask, target_mask);
        }

        prop::row("keep orientation", "Also turn everything so the target keeps the orientation it has in the first frame.\n"
                                      "Everything is wrapped into the box around the target before it is turned,\n"
                                      "and the box turns with it.");
        if (ImGui::Checkbox("##keep-orientation", &o.fixate_orientation) && o.recenter) {
            // Takes effect on the frame shown, not the next one
            system_apply_now(app, true, o.apply_pbc, o.unwrap_structures);
        }
        break;
    }
    case SystemOperation_TrackBonds: {
        auto& bt = o.bond_tracking;
        prop::row("scope", "Where bonds are tracked: the bonds with an atom among those it picks.\n"
                           "Empty for every atom; a smaller scope is faster.");
        if (filter_input("##track-bonds-scope", &bt.scope, app->script.language, filter_context(app), 0, nullptr, "all atoms")) {
            bt.scope_dirty = true;
        }
        prop::end();
        if (bt.enabled) {
            draw_bond_tracking_details(app);
        }
        return;
    }
    default:
        break;
    }
    prop::end();
}

static bool operation_has_properties(SystemOperation op) {
    return op == SystemOperation_Center || op == SystemOperation_TrackBonds;
}

// The row of an operation: its icon, its name (open for its settings, if it has any), on or off, and removing it
static void draw_operation_row(ApplicationState* app, const OperationInfo& info) {
    ImGui::PushID(info.key);

    const bool assembly_shown = app->system.assembly >= 0;
    const bool inactive = info.moves_atoms && assembly_shown;
    const bool enabled = system_operation_enabled(*app, info.op);

    const float btn = row_button_size();
    const float spacing = 2.0f;
    const ImVec2 min = ImGui::GetCursorScreenPos();
    const float avail = ImGui::GetContentRegionAvail().x;
    const float right_x = ImGui::GetCursorPosX() + avail;

    // The icon where the representations have their handles, dimmed while the operation does nothing
    const float icon_w = icon_column_width();
    ImGui::Dummy({icon_w, btn});
    {
        const ImVec2 icon_size = ImGui::CalcTextSize(info.icon);
        const ImU32 col = ImGui::GetColorU32((enabled && !inactive) ? ImGuiCol_Text : ImGuiCol_TextDisabled);
        ImGui::GetWindowDrawList()->AddText(ImFloor(min + (ImVec2(icon_w, btn) - icon_size) * 0.5f), col, info.icon);
    }
    ImGui::SameLine(0, 0);

    const float pad = 3.0f;
    bool open = false;
    ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0, pad));
    if (operation_has_properties(info.op)) {
        open = ImGui::TreeNodeEx("##node", ImGuiTreeNodeFlags_FramePadding | ImGuiTreeNodeFlags_AllowOverlap, "%s", info.name);
    } else {
        ImGui::TreeNodeEx("##node", ImGuiTreeNodeFlags_FramePadding | ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_AllowOverlap, "%s", info.name);
    }
    ImGui::PopStyleVar();
    if (ImGui::IsItemHovered()) {
        if (inactive) ImGui::SetTooltip("%s\n\nNot applied while an assembly is shown: its operators place the coordinates as the file gives them.", info.tooltip);
        else          ImGui::SetTooltip("%s", info.tooltip);
    }

    // On or off, and removing it, from the right edge
    ImGui::SameLine(right_x - (btn + spacing) * 2, spacing);
    bool on = enabled;
    ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(pad, pad));
    if (ImGui::Checkbox("##on", &on)) {
        system_set_operation_enabled(app, info.op, on);
    }
    ImGui::PopStyleVar();
    ImGui::SetItemTooltip(enabled ? "On: applied in every frame. Off keeps its settings." : "Off, its settings kept. On: applied in every frame.");
    ImGui::SameLine(0, spacing);
    bool remove = false;
    if (ImGui::DeleteButton(ICON_FA_XMARK, {btn, btn})) {
        remove = true;
    }
    ImGui::SetItemTooltip("Remove the operation and its settings");

    if (open) {
        if (inactive) ImGui::BeginDisabled();
        draw_operation_properties(app, info.op);
        if (inactive) ImGui::EndDisabled();
        ImGui::Spacing();
        ImGui::TreePop();
    }

    if (remove) {
        system_remove_operation(app, info.op);
    }
    ImGui::PopID();
}

// What the system is shown as: as it is, in its periodic images, or as one of its assemblies. One row, as a system is
// shown as one of them; none when it can only be shown as it is.
static void draw_copies_row(ApplicationState* app) {
    const md_system_t& sys = app->mold.sys;
    const bool has_cell = md_unitcell_flags(&app->mold.state.unitcell) != 0;
    const size_t num_assemblies = sys.assembly.count;
    if (!has_cell && num_assemblies == 0) return;

    ImGui::PushID("copies");
    int* n = app->system.periodic_images;
    const bool images = app->system.assembly < 0 && (n[0] | n[1] | n[2]) != 0;
    const bool shown = images || app->system.assembly >= 0;

    const float btn = row_button_size();
    const ImVec2 min = ImGui::GetCursorScreenPos();
    const float icon_w = icon_column_width();
    ImGui::Dummy({icon_w, btn});
    ImGui::GetWindowDrawList()->AddText(ImFloor(min + (ImVec2(icon_w, btn) - ImGui::CalcTextSize(ICON_FA_CUBES)) * 0.5f),
                                        ImGui::GetColorU32(shown ? ImGuiCol_Text : ImGuiCol_TextDisabled), ICON_FA_CUBES);
    ImGui::SameLine(0, 0);
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("Shown as");
    ImGui::SetItemTooltip("The system as its coordinates are, in the periodic images of its cell (the system again in the\n"
                          "neighboring cells), or as a biological assembly the file describes: the copies of its chains its\n"
                          "operators make. Each drawn by instancing; the atoms of a copy pick and select as those of the system.\n"
                          "A representation can leave itself out of the periodic images (its properties). The operators of an\n"
                          "assembly place the coordinates as the file gives them: centering, wrapping and making whole are not\n"
                          "applied while one is shown, and coordinates they moved before stay moved.");
    ImGui::SameLine();

    char preview[64] = "Coordinates as they are";
    if (app->system.assembly >= 0 && (size_t)app->system.assembly < num_assemblies) {
        const str_t id = md_assembly_id(&sys.assembly, (size_t)app->system.assembly);
        snprintf(preview, sizeof(preview), "Assembly " STR_FMT " (%zu chains)", STR_ARG(id), assembly_num_chains(sys, (size_t)app->system.assembly));
    } else if (images) {
        snprintf(preview, sizeof(preview), "Periodic images");
    }
    ImGui::SetNextItemWidth(MAX(ImGui::GetContentRegionAvail().x, 1.0f));
    if (ImGui::BeginCombo("##shown_as", preview)) {
        if (ImGui::Selectable("Coordinates as they are", !shown)) {
            system_set_assembly(app, -1);
            n[0] = n[1] = n[2] = 0;
        }
        if (has_cell && ImGui::Selectable("Periodic images", images)) {
            system_set_assembly(app, -1);
            if ((n[0] | n[1] | n[2]) == 0) n[0] = n[1] = n[2] = 1;
        }
        for (size_t k = 0; k < num_assemblies; ++k) {
            char label[64];
            const str_t id = md_assembly_id(&sys.assembly, k);
            snprintf(label, sizeof(label), "Assembly " STR_FMT " (%zu chains)##%zu", STR_ARG(id), assembly_num_chains(sys, k), k);
            if (ImGui::Selectable(label, app->system.assembly == (int)k)) {
                system_set_assembly(app, (int)k);
            }
        }
        ImGui::EndCombo();
    }

    if (images) {
        // Indented as the settings of an operation are, under its name
        ImGui::Indent();
        prop::begin("system_operation");
        prop::row("images", "This many each way along a, b and c: 1 1 1 shows the 26 cells around it");
        ImGui::SliderInt3("##periodic_images", n, 0, PERIODIC_IMAGES_MAX, "%d", ImGuiSliderFlags_AlwaysClamp);
        prop::end();
        ImGui::Unindent();
    }
    ImGui::PopID();
}

bool draw_system_header(ApplicationState* app) {
    ASSERT(app);
    ImGui::PushID("system");

    char label[128];
    system_label(label, sizeof(label), *app);

    const float btn = row_button_size();
    const float spacing = 2.0f;
    const float right_x = ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x;

    // A band across the window. What is the system's is indented under it by the caller (draw_system_operations, then
    // the representations), without the tree's own push: the two are drawn in different places.
    const ImGuiTreeNodeFlags flags = ImGuiTreeNodeFlags_Framed | ImGuiTreeNodeFlags_SpanAvailWidth | ImGuiTreeNodeFlags_AllowOverlap |
                                     ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_DefaultOpen | ImGuiTreeNodeFlags_FramePadding;
    if (!app->system.visible) ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
    const bool open = ImGui::TreeNodeEx("##system", flags, ICON_FA_ATOM "  %s", label);
    if (!app->system.visible) ImGui::PopStyleColor();
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("The system, %zu atoms: the operations on its coordinates, what it is shown as, its representations.\n"
                          "Right click or " ICON_FA_ELLIPSIS_VERTICAL ": add an operation, or apply one to the frame shown now.",
                          (size_t)app->mold.sys.atom.count);
    }
    if (ImGui::BeginPopupContextItem("##system_menu_rc")) {
        draw_system_menu(app);
        ImGui::EndPopup();
    }

    ImGui::SameLine(right_x - (btn + spacing) * 2, spacing);
    if (ImGui::Button(ICON_FA_ELLIPSIS_VERTICAL, {btn, btn})) {
        ImGui::OpenPopup("##system_menu");
    }
    ImGui::SetItemTooltip("Add an operation, or apply one to the frame shown now");
    if (ImGui::BeginPopup("##system_menu")) {
        draw_system_menu(app);
        ImGui::EndPopup();
    }
    ImGui::SameLine(0, spacing);
    if (ImGui::Button(app->system.visible ? ICON_FA_EYE : ICON_FA_EYE_SLASH, {btn, btn})) {
        app->system.visible = !app->system.visible;
        app->representation.atom_visibility_mask_dirty = true;
        if (app->system.visible) refresh_representations(app);
    }
    ImGui::SetItemTooltip("Show/Hide the system: all of its representations");

    ImGui::PopID();
    return open;
}

void draw_system_operations(ApplicationState* app) {
    ASSERT(app);
    ImGui::PushID("system");
    bool any = false;
    for (const OperationInfo& info : operation_info) {
        if (!system_has_operation(*app, info.op)) continue;
        draw_operation_row(app, info);
        any = true;
    }
    const ImVec2 before = ImGui::GetCursorScreenPos();
    draw_copies_row(app);
    any |= ImGui::GetCursorScreenPos().y != before.y;
    if (any) {
        // What is done to the system and what it is shown as, apart from how it is drawn
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();
    }
    ImGui::PopID();
}

// The state of the tracking of bond changes in the frame shown, and its advanced settings, under the operation
static void draw_bond_tracking_details(ApplicationState* app) {
    ASSERT(app);
    auto& tracking = app->system.operations.bond_tracking;

    const md_bond_delta_t* delta = bond_tracking_delta(app);
    if (delta) {
        const size_t num_broken = md_bitset_count(delta->broken);
        ImGui::Text("%zu formed, %zu broken in this frame (%.1f ms)", delta->num_formed, num_broken, tracking.time_ms);
        if (ImGui::IsItemHovered()) {
            // The atoms of the bonds that changed
            md_bitmask_clear(&app->selection.highlight_mask);
            for (size_t i = 0; i < delta->num_formed; ++i) {
                md_bitmask_set_bit(&app->selection.highlight_mask, (uint32_t)delta->formed_pair[i].idx[0]);
                md_bitmask_set_bit(&app->selection.highlight_mask, (uint32_t)delta->formed_pair[i].idx[1]);
            }
            md_bitset_iter_t it = md_bitset_iter(delta->broken);
            while (md_bitset_iter_next(&it)) {
                if (it.idx >= app->mold.sys.bond.count) break;
                const md_atom_pair_t p = app->mold.sys.bond.pairs[it.idx];
                md_bitmask_set_bit(&app->selection.highlight_mask, (uint32_t)p.idx[0]);
                md_bitmask_set_bit(&app->selection.highlight_mask, (uint32_t)p.idx[1]);
            }
            ImGui::SetTooltip("The bonds formed are drawn as sticks and the broken ones are not drawn.\n"
                              "Their atoms are highlighted.");
        }
    } else if (!tracking.scope.result.success) {
        ImGui::TextDisabled("The scope does not pick atoms: nothing is tracked");
    } else {
        ImGui::TextDisabled("Nothing to track in this system");
    }

    ImGui::SameLine();
    // VIAL has no bond changes per frame yet: its bonds() are the system's
    const bool v2 = app->script.language == script::Language::V2;
    ImGui::BeginDisabled(v2);
    if (ImGui::SmallButton("Copy to Script")) {
        // The same as a script value, to count or plot: formed(d), broken(d), bond(d)
        char args[512] = "";
        int len = 0;
        const str_t scope = str_trim(str_from_cstr(tracking.scope.text));
        if (!str_empty(scope)) len += snprintf(args + len, sizeof(args) - len, STR_FMT, STR_ARG(scope));
        if (tracking.tolerance != MD_BOND_DELTA_TOLERANCE && len < (int)sizeof(args)) {
            len += snprintf(args + len, sizeof(args) - len, "%stolerance=%.3g", len ? ", " : "", tracking.tolerance);
        }
        if (tracking.breakable != MD_BOND_DELTA_BREAKABLE && len < (int)sizeof(args)) {
            len += snprintf(args + len, sizeof(args) - len, "%sbreakable={", len ? ", " : "");
            bool first = true;
            for (int o = MD_BOND_ORIGIN_FILE; o <= MD_BOND_ORIGIN_INFERRED && len < (int)sizeof(args); ++o) {
                if (!(tracking.breakable & (1u << o))) continue;
                len += snprintf(args + len, sizeof(args) - len, "%s'%s'", first ? "" : ", ", md_bond_origin_name((md_bond_origin_t)o));
                first = false;
            }
            if (len < (int)sizeof(args)) len += snprintf(args + len, sizeof(args) - len, "}");
        }
        const str_t ident = create_unique_identifier(app->script.ir, STR_LIT("bonds_delta"), app->allocator.frame);
        char line[640];
        snprintf(line, sizeof(line), STR_FMT " = bond_delta(%s);", STR_ARG(ident), args);
        script_editor::append_line(app->editor, str_from_cstr(line));
    }
    ImGui::EndDisabled();
    if (v2) {
        ImGui::SetItemTooltip("VIAL does not track bond changes yet: bond() gives the bonds of the system");
    } else {
        ImGui::SetItemTooltip("Add the same bond changes to the script, to count them over the frames\n"
                              "(count(formed(d))) or select them (bond(d) is the bonds of each frame)");
    }

    if (ImGui::TreeNode("Advanced##track-bonds")) {
        float tolerance = tracking.tolerance;
        if (ImGui::SliderFloat("tolerance", &tolerance, 0.0f, 0.45f, "%.2f", ImGuiSliderFlags_AlwaysClamp)) {
            tracking.tolerance = tolerance;
        }
        ImGui::SetItemTooltip("A bond breaks beyond (1 + tolerance) times its cutoff and forms within (1 - tolerance) times it,\n"
                              "so that it does not flicker as it vibrates around the cutoff");
        ImGui::TextUnformatted("breakable:");
        for (int o = MD_BOND_ORIGIN_FILE; o <= MD_BOND_ORIGIN_INFERRED; ++o) {
            ImGui::SameLine();
            bool on = (tracking.breakable & (1u << o)) != 0;
            if (ImGui::Checkbox(md_bond_origin_name((md_bond_origin_t)o), &on)) {
                tracking.breakable = on ? (tracking.breakable | (1u << o)) : (tracking.breakable & ~(1u << o));
            }
        }
        ImGui::SetItemTooltip("The bonds of which origins may break. A force field topology holds throughout a\n"
                              "classical simulation, and bonds you added stay unless you say otherwise.");
        if (ImGui::SmallButton("Defaults")) {
            tracking.tolerance = MD_BOND_DELTA_TOLERANCE;
            tracking.breakable = MD_BOND_DELTA_BREAKABLE;
        }
        ImGui::TreePop();
    }
}

// ---------------------------------------------------------------------------------------------------------
// Workspace
// ---------------------------------------------------------------------------------------------------------

void system_scene_serialize(viamd::serialization_state_t& state, const ApplicationState* app, int system_id) {
    ASSERT(app);
    const auto& s = app->system;
    const auto& op = s.operations;

    viamd::write_section_header(state, STR_LIT("SystemScene"));
    viamd::write_int(state, STR_LIT("System"), system_id);
    viamd::write_bool(state, STR_LIT("Visible"), s.visible);

    // The operations it has, in their order; the settings of those alone
    char list[128] = "";
    int len = 0;
    for (const OperationInfo& info : operation_info) {
        if (!system_has_operation(*app, info.op)) continue;
        len += snprintf(list + len, sizeof(list) - len, "%s%s", len ? ", " : "", info.key);
    }
    viamd::write_str(state, STR_LIT("Operations"), str_from_cstr(list));

    if (system_has_operation(*app, SystemOperation_Center)) {
        viamd::write_bool(state, STR_LIT("Recenter"), op.recenter);
        viamd::write_bool(state, STR_LIT("FixateOrientation"), op.fixate_orientation);
        viamd::write_bool(state, STR_LIT("RecenterQueryEnabled"), op.recenter_query.enabled);
        viamd::write_str(state,  STR_LIT("RecenterQuery"), str_from_cstr(op.recenter_query.filter.text));
        if (!md_bitset_empty(op.selection_mask)) {
            viamd::write_bitset(state, STR_LIT("RecenterTarget"), &op.selection_mask);
        }
    }
    if (system_has_operation(*app, SystemOperation_Wrap)) {
        viamd::write_bool(state, STR_LIT("ApplyPbc"), op.apply_pbc);
    }
    if (system_has_operation(*app, SystemOperation_MakeWhole)) {
        viamd::write_bool(state, STR_LIT("UnwrapStructures"), op.unwrap_structures);
    }
    if (system_has_operation(*app, SystemOperation_TrackBonds)) {
        viamd::write_bool(state, STR_LIT("TrackBonds"), op.bond_tracking.enabled);
        viamd::write_str(state,  STR_LIT("TrackBondsScope"), str_from_cstr(op.bond_tracking.scope.text));
        viamd::write_flt(state,  STR_LIT("TrackBondsTolerance"), op.bond_tracking.tolerance);
        char origins[64];
        viamd::write_str(state,  STR_LIT("TrackBondsBreakable"), str_from_cstr(bond_origins_str(origins, sizeof(origins), op.bond_tracking.breakable)));
    }

    // What it is shown as
    if ((s.periodic_images[0] | s.periodic_images[1] | s.periodic_images[2]) != 0) {
        viamd::write_int_vec(state, STR_LIT("PeriodicImages"), s.periodic_images, 3);
    }
    if (s.assembly >= 0 && (size_t)s.assembly < app->mold.sys.assembly.count) {
        viamd::write_str(state, STR_LIT("Assembly"), md_assembly_id(&app->mold.sys.assembly, (size_t)s.assembly));
    }
}

bool system_scene_read_entry(ApplicationState* app, str_t ident, str_t arg) {
    ASSERT(app);
    auto& s = app->system;
    auto& op = s.operations;

    if (str_eq(ident, STR_LIT("Visible"))) {
        viamd::extract_bool(s.visible, arg);
    } else if (str_eq(ident, STR_LIT("Operations"))) {
        str_t list = {};
        if (!viamd::extract_str(list, arg)) return true;
        op.present = 0;
        while (!str_empty(list)) {
            size_t loc = list.len;
            str_find_char(&loc, list, ',');
            const str_t tok = str_trim(str_substr(list, 0, loc));
            list = loc < list.len ? str_substr(list, loc + 1) : str_t{};
            if (str_empty(tok)) continue;
            bool known = false;
            for (const OperationInfo& info : operation_info) {
                if (str_eq_cstr(tok, info.key)) {
                    op.present |= info.op;
                    known = true;
                }
            }
            if (!known) VIAMD_LOG_INFO("Workspace: no operation '" STR_FMT "' in this build of viamd, it is left out", STR_ARG(tok));
        }
    } else if (str_eq(ident, STR_LIT("Recenter"))) {
        viamd::extract_bool(op.recenter, arg);
    } else if (str_eq(ident, STR_LIT("FixateOrientation"))) {
        viamd::extract_bool(op.fixate_orientation, arg);
    } else if (str_eq(ident, STR_LIT("ApplyPbc"))) {
        viamd::extract_bool(op.apply_pbc, arg);
    } else if (str_eq(ident, STR_LIT("UnwrapStructures"))) {
        viamd::extract_bool(op.unwrap_structures, arg);
    } else if (str_eq(ident, STR_LIT("RecalcBonds")) || str_eq(ident, STR_LIT("TrackBonds"))) {
        // RecalcBonds: what a workspace from before the bond changes were tracked had, the bonds of every frame guessed anew
        viamd::extract_bool(op.bond_tracking.enabled, arg);
    } else if (str_eq(ident, STR_LIT("TrackBondsScope"))) {
        viamd::extract_to_char_buf(op.bond_tracking.scope.text, sizeof(op.bond_tracking.scope.text), arg);
        op.bond_tracking.scope_dirty = true;
    } else if (str_eq(ident, STR_LIT("TrackBondsTolerance"))) {
        float tolerance = MD_BOND_DELTA_TOLERANCE;
        if (viamd::extract_flt(tolerance, arg)) op.bond_tracking.tolerance = CLAMP(tolerance, 0.0f, 0.45f);
    } else if (str_eq(ident, STR_LIT("TrackBondsBreakable"))) {
        str_t str;
        uint32_t mask = 0;
        if (viamd::extract_str(str, arg) && bond_origins_parse(&mask, str)) op.bond_tracking.breakable = mask;
    } else if (str_eq(ident, STR_LIT("RecenterQueryEnabled"))) {
        viamd::extract_bool(op.recenter_query.enabled, arg);
    } else if (str_eq(ident, STR_LIT("RecenterQuery"))) {
        viamd::extract_to_char_buf(op.recenter_query.filter.text, sizeof(op.recenter_query.filter.text), arg);
        recenter_mark_query_dirty(app);
    } else if (str_eq(ident, STR_LIT("PeriodicImages"))) {
        // Integers; a workspace written before [SystemScene] had them as floats ([RenderSettings])
        float n[3];
        if (viamd::extract_flt_vec(n, 3, arg)) {
            for (int k = 0; k < 3; ++k) s.periodic_images[k] = CLAMP((int)n[k], 0, PERIODIC_IMAGES_MAX);
        }
    } else {
        return false;
    }
    return true;
}

void system_scene_from_legacy(ApplicationState* app) {
    ASSERT(app);
    auto& op = app->system.operations;
    system_operations_sync(app);
    if (op.recenter_query.enabled || op.recenter_query.filter.text[0] || op.fixate_orientation) op.present |= SystemOperation_Center;
    if (op.bond_tracking.scope.text[0]) op.present |= SystemOperation_TrackBonds;
}
