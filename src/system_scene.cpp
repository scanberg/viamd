#include <system_scene.h>

#include <viamd.h>
#include <representation.h>
#include <filter_input.h>
#include <imgui_widgets.h>
#include <property_grid.h>
#include <serialization_utils.h>
#include <script_editor.h>
#include <plot_series.h>
#include <display_units.h>

#include <core/md_bitset.h>
#include <core/md_str.h>
#include <core/md_path.h>

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
    MD_ASSERT(false);
    return operation_info[0];
}

SystemOperation system_operation_at(int idx) {
    MD_ASSERT(0 <= idx && idx < SYSTEM_OPERATION_COUNT);
    return operation_info[idx].op;
}

const char* system_operation_name(SystemOperation op) {
    return info_of(op).name;
}

// The flag that turns it on and off
static bool* enabled_flag(SystemContext& sys, SystemOperation op) {
    auto& o = sys.operations;
    switch (op) {
    case SystemOperation_Center:     return &o.recenter;
    case SystemOperation_Wrap:       return &o.apply_pbc;
    case SystemOperation_MakeWhole:  return &o.unwrap_structures;
    case SystemOperation_TrackBonds: return &o.bond_tracking.enabled;
    }
    MD_ASSERT(false);
    return &o.recenter;
}

bool system_operation_enabled(const SystemContext& sys, SystemOperation op) {
    const auto& o = sys.operations;
    switch (op) {
    case SystemOperation_Center:     return o.recenter;
    case SystemOperation_Wrap:       return o.apply_pbc;
    case SystemOperation_MakeWhole:  return o.unwrap_structures;
    case SystemOperation_TrackBonds: return o.bond_tracking.enabled;
    }
    return false;
}


void system_apply_now(SystemContext& sys, bool center, bool wrap, bool make_whole) {
    for (size_t k = 0; k < sys.num_instances; ++k) {
        SystemInstance& inst = *sys.instances[k];
        // Its operators place the coordinates as the file gives them (scene.assembly)
        if (inst.scene.assembly >= 0) continue;
        if (apply_state_operations(inst, center, wrap, make_whole)) {
            inst.dirty_gpu_buffers |= MolBit_DirtyPosition | MolBit_ClearVelocity;
        }
    }
}

void system_set_operation_enabled(SystemContext& sys, SystemOperation op, bool enabled) {
    auto& o = sys.operations;
    *enabled_flag(sys, op) = enabled;

    if (op == SystemOperation_TrackBonds) {
        o.bond_tracking.scope_dirty = true;
    } else if (enabled) {
        // The frame shown as every frame will be: all that is on, in its order, the one turned on with it
        system_apply_now(sys, o.recenter, o.apply_pbc, o.unwrap_structures);
    }
}

void system_reset_operations(SystemContext& sys) {
    ApplicationState* app = sys.app;
    MD_ASSERT(app);
    auto& o = sys.operations;
    o.recenter = false;
    o.fixate_orientation = false;
    o.recenter_query.enabled = false;
    o.recenter_query.filter.text[0] = '\0';
    o.recenter_query.filter.result = {};
    md_bitset_free(&o.selection_mask, app->allocator.persistent);
    o.selection_mask = md_bitset_empty_set();
    recenter_mark_query_dirty(sys);

    o.apply_pbc = false;
    o.unwrap_structures = false;

    auto& bt = o.bond_tracking;
    bt.enabled = false;
    bt.scope.text[0] = '\0';
    bt.scope.result = {};
    bt.scope_dirty = true;
    bt.tolerance = MD_BOND_DELTA_TOLERANCE;
    bt.breakable = MD_BOND_DELTA_BREAKABLE;
}

// ---------------------------------------------------------------------------------------------------------
// Copies and naming
// ---------------------------------------------------------------------------------------------------------

void instance_set_assembly(SystemInstance& inst, int assembly, bool frame) {
    SystemContext& sys = *inst.sys;
    ApplicationState* app = sys.app;
    MD_ASSERT(app);
    if (assembly == inst.scene.assembly) return;
    // Its copies at once (they are set every frame, update_view_param), and what is shown changes its extent: the
    // view frames it
    inst.scene.assembly = assembly;
    assembly_transforms(&inst.scene.assembly_transforms, &inst.scene.assembly_copy_offset, sys.data.sys, assembly, app->allocator.persistent);
    if (!frame) return;
    reset_view(&app->view.target, inst.state, &sys.representation.visibility_mask, app->view.camera.fov_y, view_aspect_ratio(*app), &sys.data.sys, assembly);
    // Room to zoom out to all of it
    app->view.trackball_param.max_distance = MD_MAX(app->view.trackball_param.max_distance, app->view.target.distance * 10.0f);
}

void system_label(char* buf, size_t cap, const SystemContext& sys) {
    const char* name = sys.name;
    snprintf(buf, cap, "%s", name[0] ? name : "System");
}

// Called after where its coordinates come from: its run by the name of its file ("run/<name>"), at the frame it keeps to,
// or the system's own coordinates
static void instance_source_label(char* buf, size_t cap, const SystemInstance& inst) {
    const SystemRun* run = instance_run(inst);
    const char* name = "Own coordinates";
    if (run) {
        name = run->path;
        if (strncmp(name, "run/", 4) == 0) name += 4;
    }
    if (run && inst.frame >= 0) snprintf(buf, cap, "%s @ frame %lld", name, (long long)inst.frame);
    else                        snprintf(buf, cap, "%s", name);
}

bool system_combo(const char* label, uint32_t* system_id, const ApplicationState& app, bool (*candidate)(const SystemContext&)) {
    auto listed = [&](const SystemContext& sc) { return system_loaded(sc) && (!candidate || candidate(sc)); };
    size_t num_listed = 0;
    for (size_t i = 0; i < app.num_systems; ++i) num_listed += listed(*app.systems[i]) ? 1 : 0;
    if (num_listed < 2) return false;
    const SystemContext& cur = system_or_first(app, *system_id);
    char preview[80];
    system_label(preview, sizeof(preview), cur);
    bool changed = false;
    if (ImGui::BeginCombo(label, preview)) {
        for (size_t i = 0; i < app.num_systems; ++i) {
            const SystemContext& sc = *app.systems[i];
            if (!listed(sc)) continue;
            char item[80];
            system_label(item, sizeof(item), sc);
            ImGui::PushID((int)sc.id);
            if (ImGui::Selectable(item, &sc == &cur) && sc.id != *system_id) {
                *system_id = sc.id;
                changed = true;
            }
            ImGui::PopID();
        }
        ImGui::EndCombo();
    }
    return changed;
}

bool instance_combo(const char* label, uint32_t* instance_id, const SystemContext& sys) {
    if (sys.num_instances < 2) return false;
    const SystemInstance& cur = instance_resolve(sys, *instance_id);
    char preview[96];
    instance_label(preview, sizeof(preview), cur);
    bool changed = false;
    if (ImGui::BeginCombo(label, preview)) {
        for (size_t k = 0; k < sys.num_instances; ++k) {
            const SystemInstance& inst = *sys.instances[k];
            char item[96];
            instance_label(item, sizeof(item), inst);
            ImGui::PushID((int)inst.id);
            const uint32_t id = k == 0 ? 0 : inst.id;
            if (ImGui::Selectable(item, &inst == &cur) && id != *instance_id) {
                *instance_id = id;
                changed = true;
            }
            ImGui::PopID();
        }
        ImGui::EndCombo();
    }
    return changed;
}

void instance_label(char* buf, size_t cap, const SystemInstance& inst) {
    if (inst.name[0]) {
        snprintf(buf, cap, "%s", inst.name);
        return;
    }
    instance_source_label(buf, cap, inst);
    // Of several called the same, the second and on by their place among them
    const SystemContext& sys = *inst.sys;
    size_t same = 0;
    for (size_t k = 0; k < sys.num_instances && sys.instances[k] != &inst; ++k) {
        char other[96];
        if (sys.instances[k]->name[0]) continue;
        instance_source_label(other, sizeof(other), *sys.instances[k]);
        if (strcmp(other, buf) == 0) same += 1;
    }
    if (same > 0) {
        const size_t len = strlen(buf);
        if (len < cap) snprintf(buf + len, cap - len, " (%zu)", same + 1);
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
    return MD_MAX(ImGui::CalcTextSize(ICON_FA_GRIP_VERTICAL).x,
           MD_MAX(ImGui::CalcTextSize(ICON_FA_FOLDER).x, ImGui::CalcTextSize(ICON_FA_FOLDER_OPEN).x)) + 6.0f;
}

static float row_button_size() {
    return ImGui::GetFontSize() + 3.0f * 2;
}

// The representations were not updated while the system was hidden; brought up to what they show now
static void refresh_representations(SystemContext& sys) {
    for (size_t i = 0; i < md_array_size(sys.representation.reps); ++i) {
        sys.representation.reps[i].needs_update = true;
    }
    flag_dynamic_representations_as_dirty(sys);
}

// The system's menu, of its header's right click and of its ⋮ button: what is done once, the details, removing it
static void draw_system_menu(SystemContext& sys, bool* remove_request, InstanceActions* inst_act) {
    ApplicationState* app = sys.app;
    // Every instance shows an assembly: there is nothing the operations move
    bool assembly_shown = true;
    for (size_t k = 0; k < sys.num_instances; ++k) assembly_shown &= sys.instances[k]->scene.assembly >= 0;
    const bool has_cell = md_unitcell_flags(&sys.data.state.unitcell) != 0;

    ImGui::TextDisabled("Now, to the frame shown");

    const md_bitset_t& target_mask = recenter_get_active_target_mask(sys);
    const bool target = !md_bitset_empty(target_mask);
    if (ImGui::MenuItem(ICON_FA_CROSSHAIRS " Center target", nullptr, false, target && !assembly_shown)) {
        system_apply_now(sys, true, false, false);
    }
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
        if (target) {
            md_bitmask_assign(&sys.selection.highlight_mask, target_mask);
            ImGui::SetTooltip("Move everything so the target (highlighted) is in the middle of the box,\nor at the origin when there is no box");
        } else {
            ImGui::SetTooltip("There is nothing to center on yet: give the Center operation a target (Operations)");
        }
    }
    if (ImGui::MenuItem(ICON_FA_CUBE " Wrap into box", nullptr, false, has_cell && !assembly_shown)) {
        system_apply_now(sys, false, true, false);
    }
    ImGui::SetItemTooltip("Move every atom that is outside the periodic box back in through the opposite side");
    if (ImGui::MenuItem(ICON_FA_PUZZLE_PIECE " Make whole", nullptr, false, !assembly_shown)) {
        system_apply_now(sys, false, false, true);
    }
    ImGui::SetItemTooltip("Join molecules that are split across the periodic box boundary");
    if (ImGui::MenuItem(ICON_FA_ARROWS_ROTATE " Recompute bonds")) {
        // The whole frame nearest the animation time, or the coordinates shown when there is no run
        recompute_bonds(sys, (int64_t)(instance_current_frame(system_instance(sys)) + 0.5));
    }
    ImGui::SetItemTooltip("Guess the covalent bonds again from the distances between atoms in this frame, for when\n"
                          "the file's bonds are wrong or a reaction should be the topology from here on.\n"
                          "Changes the bonds of the system: the bonds a file, a topology or you gave stay.");
    if (assembly_shown) {
        ImGui::TextDisabled("Centering, wrapping and making whole are not\napplied while an assembly is shown");
    }

    ImGui::Separator();
    if (ImGui::MenuItem(ICON_FA_CLONE " Add instance", nullptr, false, sys.num_instances < INSTANCE_MAX)) {
        inst_act->add_like = system_instance(sys).id;
    }
    ImGui::SetItemTooltip("Another instance of the system, as the first is: of the same coordinates, with copies of its\n"
                          "representations. Shown beside the others, of another run or frame, placed elsewhere.");
    if (ImGui::MenuItem(ICON_FA_FILE_CIRCLE_PLUS " Add File...")) {
        char path[1024] = "";
        if (application::file_dialog(path, sizeof(path), application::FileDialogFlag_Open)) {
            file_queue_push(&app->file_queue, md_str_from_cstr(path), FileFlags_None, sys.id);
        }
    }
    ImGui::SetItemTooltip("Load a file into this system: a trajectory becomes its run, energies or a topology add to it,\n"
                          "a structure replaces it");
    if (bool* show = workspace_window_flag("System")) {
        // The System window is of the first system until it has a system picker of its own
        const bool first = &sys == &first_system(*app);
        if (ImGui::MenuItem(ICON_FA_CIRCLE_INFO " Details", nullptr, false, first)) {
            *show = true;
        }
        if (first) ImGui::SetItemTooltip("The System window: its files in full, what it is made of, its trajectory and series,\nits atom and residue types");
        else       ImGui::SetItemTooltip("The System window shows the first system only, so far");
    }
    if (ImGui::MenuItem(ICON_FA_TRASH_CAN " Remove system...")) {
        *remove_request = true;
    }
    ImGui::SetItemTooltip("Unload the system, its files, operations and representations");
}

static void draw_bond_tracking_details(SystemContext& sys);

// The settings of an operation, under its row
static void draw_operation_properties(SystemContext& sys, SystemOperation op) {
    ApplicationState* app = sys.app;
    auto& o = sys.operations;
    prop::begin("system_operation");
    switch (op) {
    case SystemOperation_Center: {
        const md_bitset_t& target_mask = recenter_get_active_target_mask(sys);

        prop::row("target", "What is centered: the atoms set from a selection, or those a query picks\n"
                            "(evaluated again in every frame when it depends on it)");
        bool by_query = o.recenter_query.enabled;
        if (ImGui::RadioButton("selection", !by_query)) by_query = false;
        ImGui::SameLine();
        if (ImGui::RadioButton("query", by_query)) by_query = true;
        if (by_query != o.recenter_query.enabled) {
            o.recenter_query.enabled = by_query;
            if (o.recenter) system_apply_now(sys, true, o.apply_pbc, o.unwrap_structures);
        }

        if (o.recenter_query.enabled) {
            prop::row("query");
            if (filter_input("##recenter-query", &o.recenter_query.filter, app->script.language, filter_context(sys))) {
                recenter_mark_query_dirty(sys);
            }
        } else {
            prop::row("atoms");
            const size_t num_selected = md_bitmask_count(&sys.selection.selection_mask);
            ImGui::BeginDisabled(num_selected == 0);
            if (ImGui::SmallButton("set from selection")) {
                md_bitset_assign(&o.selection_mask, md_bitmask_view(&sys.selection.selection_mask), app->allocator.persistent);
                recenter_update_target_data(sys);
                if (o.recenter) system_apply_now(sys, true, o.apply_pbc, o.unwrap_structures);
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
            md_bitmask_assign(&sys.selection.highlight_mask, target_mask);
        }

        prop::row("keep orientation", "Also turn everything so the target keeps the orientation it has in the first frame.\n"
                                      "Everything is wrapped into the box around the target before it is turned,\n"
                                      "and the box turns with it.");
        if (ImGui::Checkbox("##keep-orientation", &o.fixate_orientation) && o.recenter) {
            // Takes effect on the frame shown, not the next one
            system_apply_now(sys, true, o.apply_pbc, o.unwrap_structures);
        }
        break;
    }
    case SystemOperation_TrackBonds: {
        auto& bt = o.bond_tracking;
        prop::row("scope", "Where bonds are tracked: the bonds with an atom among those it picks.\n"
                           "Empty for every atom; a smaller scope is faster.");
        if (filter_input("##track-bonds-scope", &bt.scope, app->script.language, filter_context(sys), 0, nullptr, "all atoms")) {
            bt.scope_dirty = true;
        }
        prop::end();
        if (bt.enabled) {
            draw_bond_tracking_details(sys);
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

// The row of an operation: its icon, its name (open for its settings, if it has any), and on or off
static void draw_operation_row(SystemContext& sys, const OperationInfo& info) {
    ImGui::PushID(info.key);

    const bool assembly_shown = system_instance(sys).scene.assembly >= 0;
    const bool inactive = info.moves_atoms && assembly_shown;
    const bool enabled = system_operation_enabled(sys, info.op);

    const float btn = row_button_size();
    const float spacing = 0.0f;
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

    // On or off, at the right edge
    ImGui::SameLine(right_x - btn, spacing);
    bool on = enabled;
    ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(pad, pad));
    if (ImGui::Checkbox("##on", &on)) {
        system_set_operation_enabled(sys, info.op, on);
    }
    ImGui::PopStyleVar();
    ImGui::SetItemTooltip(enabled ? "On: applied in every frame. Off keeps its settings." : "Off, its settings kept. On: applied in every frame.");

    if (open) {
        if (inactive) ImGui::BeginDisabled();
        draw_operation_properties(sys, info.op);
        if (inactive) ImGui::EndDisabled();
        ImGui::Spacing();
        ImGui::TreePop();
    }

    ImGui::PopID();
}

// What the instance is shown as: as it is, in its periodic images, or as one of its assemblies. One row, as an instance
// is shown as one of them; none when it can only be shown as it is.
void draw_instance_copies_row(SystemInstance& inst) {
    const SystemContext& sc = *inst.sys;
    const md_system_t& sys = sc.data.sys;
    const bool has_cell = md_unitcell_flags(&inst.state.unitcell) != 0;
    const size_t num_assemblies = sys.assembly.count;
    if (!has_cell && num_assemblies == 0) return;

    ImGui::PushID("copies");
    int* n = inst.scene.periodic_images;
    const bool images = inst.scene.assembly < 0 && (n[0] | n[1] | n[2]) != 0;
    const bool shown = images || inst.scene.assembly >= 0;

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
    if (inst.scene.assembly >= 0 && (size_t)inst.scene.assembly < num_assemblies) {
        const md_str_t id = md_assembly_id(&sys.assembly, (size_t)inst.scene.assembly);
        snprintf(preview, sizeof(preview), "Assembly " MD_STR_FMT " (%zu chains)", MD_STR_ARG(id), assembly_num_chains(sys, (size_t)inst.scene.assembly));
    } else if (images) {
        snprintf(preview, sizeof(preview), "Periodic images");
    }
    ImGui::SetNextItemWidth(MD_MAX(ImGui::GetContentRegionAvail().x, 1.0f));
    if (ImGui::BeginCombo("##shown_as", preview)) {
        if (ImGui::Selectable("Coordinates as they are", !shown)) {
            instance_set_assembly(inst, -1);
            n[0] = n[1] = n[2] = 0;
        }
        if (has_cell && ImGui::Selectable("Periodic images", images)) {
            instance_set_assembly(inst, -1);
            if ((n[0] | n[1] | n[2]) == 0) n[0] = n[1] = n[2] = 1;
        }
        for (size_t k = 0; k < num_assemblies; ++k) {
            char label[64];
            const md_str_t id = md_assembly_id(&sys.assembly, k);
            snprintf(label, sizeof(label), "Assembly " MD_STR_FMT " (%zu chains)##%zu", MD_STR_ARG(id), assembly_num_chains(sys, k), k);
            if (ImGui::Selectable(label, inst.scene.assembly == (int)k)) {
                instance_set_assembly(inst, (int)k);
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

bool draw_system_header(SystemContext& sys, bool* remove, InstanceActions* inst_act) {
    ImGui::PushID("system");

    char label[128];
    system_label(label, sizeof(label), sys);

    const float btn = row_button_size();
    const float spacing = 2.0f;
    const float right_x = ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x;

    // A band across the window. What is the system's is indented under it by the caller (draw_system_operations, then
    // the representations), without the tree's own push: the two are drawn in different places.
    const ImGuiTreeNodeFlags flags = ImGuiTreeNodeFlags_Framed | ImGuiTreeNodeFlags_SpanAvailWidth | ImGuiTreeNodeFlags_AllowOverlap |
                                     ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_DefaultOpen | ImGuiTreeNodeFlags_FramePadding;
    if (!sys.visible) ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
    const bool open = ImGui::TreeNodeEx("##system", flags, ICON_FA_ATOM "  %s", label);
    if (!sys.visible) ImGui::PopStyleColor();
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("The system, %zu atoms: the operations on its coordinates, what it is shown as, its representations.\n"
                          "Right click or " ICON_FA_ELLIPSIS_VERTICAL ": apply an operation to the frame shown now, recompute the bonds, the details.",
                          (size_t)sys.data.sys.atom.count);
    }
    bool remove_request = false;
    if (ImGui::BeginPopupContextItem("##system_menu_rc")) {
        draw_system_menu(sys, &remove_request, inst_act);
        ImGui::EndPopup();
    }

    // As the rows of the representations: show/hide, more, remove
    ImGui::SameLine(right_x - (btn + spacing) * 3, spacing);
    if (ImGui::Button(sys.visible ? ICON_FA_EYE : ICON_FA_EYE_SLASH, {btn, btn})) {
        sys.visible = !sys.visible;
        sys.representation.atom_visibility_mask_dirty = true;
        if (sys.visible) refresh_representations(sys);
    }
    ImGui::SetItemTooltip("Show/Hide the system: all of its representations");
    ImGui::SameLine(0, spacing);
    if (ImGui::Button(ICON_FA_ELLIPSIS_VERTICAL, {btn, btn})) {
        ImGui::OpenPopup("##system_menu");
    }
    ImGui::SetItemTooltip("Apply an operation to the frame shown now, recompute the bonds, the details");
    if (ImGui::BeginPopup("##system_menu")) {
        draw_system_menu(sys, &remove_request, inst_act);
        ImGui::EndPopup();
    }
    ImGui::SameLine(0, spacing);
    if (ImGui::DeleteButton(ICON_FA_XMARK, {btn, btn})) {
        remove_request = true;
    }
    ImGui::SetItemTooltip("Remove the system: asked first");

    // Asked before it is done: it cannot be taken back. Done by the caller, once the window is drawn.
    if (remove_request) {
        ImGui::OpenPopup("Remove system?##confirm");
    }
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    if (ImGui::BeginPopupModal("Remove system?##confirm", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::Text("Remove '%s'?", label);
        ImGui::TextDisabled("Its files are unloaded, and its operations and representations go with it.\n"
                            "This cannot be undone; a saved workspace still has it.");
        ImGui::Spacing();
        if (ImGui::DeleteButton("Remove")) {
            *remove = true;
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel") || ImGui::IsKeyPressed(ImGuiKey_Escape)) {
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }

    ImGui::PopID();
    return open;
}

// ---------------------------------------------------------------------------------------------------------
// Instances
// ---------------------------------------------------------------------------------------------------------

// What an instance is turned about by its placement: the centre of the system's cell, or of its own coordinates without
// one, so that turning it turns it in place
static md_float3_t instance_pivot(const SystemInstance& inst) {
    const SystemContext& sys = *inst.sys;
    if (md_unitcell_flags(&sys.data.state.unitcell) != 0) {
        md_float3x3_t A;
        md_unitcell_A_extract_float(A.elem, &sys.data.state.unitcell);
        return md_float3x3_mul_float3(A, md_float3_set(0.5f, 0.5f, 0.5f));
    }
    const md_system_state_t& st = sys.data.state;
    md_float3_t c = md_float3_zero();
    if (!st.xyz || st.num_atoms == 0) return c;
    for (size_t i = 0; i < st.num_atoms; ++i) c = md_float3_add(c, st.xyz[i]);
    return md_float3_mul1(c, 1.0f / (float)st.num_atoms);
}

// The placement as an offset and turns about x, y and z (degrees, applied in that order) about the pivot, and back
static void placement_decompose(md_float3_t* offset, md_float3_t* degrees, const md_float4x4_t& M, md_float3_t pivot) {
    // R(i, j) is M.elem[j][i]: column major
    const float r20 = M.elem[0][2], r21 = M.elem[1][2], r22 = M.elem[2][2], r10 = M.elem[0][1], r00 = M.elem[0][0];
    const float to_deg = 180.0f / 3.14159265f;
    // + 0: no -0 shown
    degrees->x = atan2f(r21, r22) * to_deg + 0.0f;
    degrees->y = asinf(MD_CLAMP(-r20, -1.0f, 1.0f)) * to_deg + 0.0f;
    degrees->z = atan2f(r10, r00) * to_deg + 0.0f;
    // M p = R (p - c) + c + offset: the translation of M is c + offset - R c
    const md_float3_t Rc = md_float4x4_mul_float3(M, pivot, 0.0f);
    const md_float3_t t  = md_float3_set(M.elem[3][0], M.elem[3][1], M.elem[3][2]);
    *offset = md_float3_add(md_float3_sub(t, pivot), Rc);
}

static md_float4x4_t placement_compose(md_float3_t offset, md_float3_t degrees, md_float3_t pivot) {
    const float to_rad = 3.14159265f / 180.0f;
    const md_quat_t q = md_quat_mul(md_quat_axis_angle(md_float3_set(0, 0, 1), degrees.z * to_rad),
                        md_quat_mul(md_quat_axis_angle(md_float3_set(0, 1, 0), degrees.y * to_rad),
                                    md_quat_axis_angle(md_float3_set(1, 0, 0), degrees.x * to_rad)));
    const md_float4x4_t R = md_float4x4_from_quat(q);
    return md_float4x4_mul(md_float4x4_translate_float3(md_float3_add(pivot, offset)), md_float4x4_mul(R, md_float4x4_translate_float3(md_float3_mul1(pivot, -1.0f))));
}

bool instance_rows_shown(const SystemContext& sys) {
    if (sys.num_instances != 1) return true;
    // Of the one instance, what the system's rows would not tell: a name, a frame it keeps to, a placement
    const SystemInstance& inst = system_instance(sys);
    return inst.name[0] || (inst.run >= 0 && inst.frame >= 0) || !md_float4x4_equal(inst.scene.placement, md_float4x4_ident()) || !inst.scene.visible;
}

void draw_run_time_offset(SystemContext& sys, const SystemRun* run) {
    const ApplicationState* app = sys.app;
    const ClockAxis a = run_clock_axis(app, sys, run);
    if (!clock_axis_timed(a)) {
        ImGui::TextDisabled("Its frames go with the first run's, frame by frame:\nits times cannot be placed beside the first run's");
        return;
    }
    // In the display unit of its times; frames numbered show as they are
    char unit[32];
    const double scl = display_units::factor_print(unit, sizeof(unit), run_time_unit(sys, run));
    char fmt[48];
    snprintf(fmt, sizeof(fmt), unit[0] ? "%%+.6g %s" : "%%+.6g", unit);
    double shown = run_time_offset(sys, run) * scl;
    const double span = a.count > 1 ? (a.times[a.count - 1] - a.times[0]) * scl : 1.0;
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
    if (ImGui::DragScalar("##time_offset", ImGuiDataType_Double, &shown, (float)MD_MAX(span / 500.0, 1e-6), nullptr, nullptr, fmt)) {
        run_set_time_offset(sys, run, shown / scl);
    }
    ImGui::SetItemTooltip("The run moved along the time shown: added to its times. Shared by every instance of the run.");
    ImGui::SameLine();
    const SystemInstance* ref = clock_reference(app);
    ImGui::BeginDisabled(!ref || instance_run(*ref) == run);
    if (ImGui::Button("Align start")) run_align_start(sys, run);
    ImGui::EndDisabled();
    ImGui::SetItemTooltip("%s", ref && instance_run(*ref) == run ? "It is the first run, which the others are aligned with"
                                                                 : "Its first frame at the first frame of the first run: the two start together");
    ImGui::SameLine();
    ImGui::BeginDisabled(run_time_offset(sys, run) == 0.0);
    if (ImGui::Button("As in its file")) run_set_time_offset(sys, run, 0.0);
    ImGui::EndDisabled();
    ImGui::SetItemTooltip("Its times as its file has them");
}

// The menu of an instance, of its row's right click and of its ⋮ button
static void draw_instance_menu(SystemInstance& inst, InstanceActions* act) {
    SystemContext& sys = *inst.sys;
    const bool first = &inst == &system_instance(sys);

    char label[96];
    instance_label(label, sizeof(label), inst);
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14);
    ImGui::InputTextWithHint("##name", label, inst.name, sizeof(inst.name));
    ImGui::SetItemTooltip("Its name; without one it is called after where its coordinates come from");

    ImGui::SeparatorText("Coordinates");
    if (ImGui::RadioButton("The system's own", inst.run < 0)) instance_show(inst, -1);
    ImGui::SetItemTooltip("Those its files gave, or the first frame of its first run");
    for (size_t r = 0; r < sys.num_runs; ++r) {
        const SystemRun& run = *sys.runs[r];
        if (run_num_frames(sys, &run) == 0) continue;   // held in its place, not opened
        const char* name = strncmp(run.path, "run/", 4) == 0 ? run.path + 4 : run.path;
        char item[160];
        snprintf(item, sizeof(item), "%s##run%zu", name, r);
        if (ImGui::RadioButton(item, inst.run == (int)r)) instance_show(inst, (int)r, inst.run == (int)r ? inst.frame : -1);
    }
    if (const SystemRun* run = instance_run(inst)) {
        const size_t num_frames = run_num_frames(sys, run);
        ImGui::Indent();
        if (ImGui::RadioButton("Follows the time shown", inst.frame < 0)) instance_show(inst, inst.run, -1);
        ImGui::SetItemTooltip(first ? "Its frame of the time shown" : "Its frame of the time shown: of the same time as the others', whatever its frames are numbered");
        ImGui::SameLine();
        if (ImGui::RadioButton("Keeps to a frame", inst.frame >= 0) && inst.frame < 0) {
            instance_show(inst, inst.run, (int64_t)MD_CLAMP(instance_current_frame(inst) + 0.5, 0.0, (double)(num_frames - 1)));
        }
        ImGui::SetItemTooltip("One frame of the run, whatever the time shown: a reference, a ghost of where it was");
        if (inst.frame >= 0 && num_frames > 0) {
            int frame = (int)inst.frame;
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14);
            if (ImGui::SliderInt("frame", &frame, 0, (int)num_frames - 1, "%d", ImGuiSliderFlags_AlwaysClamp)) {
                instance_show(inst, inst.run, frame);
            }
        }
        ImGui::TextDisabled("Its run on the time shown");
        draw_run_time_offset(sys, run);
        ImGui::Unindent();
    }

    ImGui::SeparatorText("Placement");
    {
        const md_float3_t pivot = instance_pivot(inst);
        md_float3_t offset, degrees;
        placement_decompose(&offset, &degrees, inst.scene.placement, pivot);
        bool changed = false;
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14);
        changed |= ImGui::DragFloat3("offset (" TOOLTIP_ANGSTROM ")", offset.elem, 0.1f, 0.0f, 0.0f, "%.2f");
        ImGui::SetItemTooltip("Where it is moved to, from where its coordinates put it");
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14);
        changed |= ImGui::DragFloat3("turn (deg)", degrees.elem, 0.5f, -180.0f, 180.0f, "%.1f");
        ImGui::SetItemTooltip("Turned about x, then y, then z, about the centre of its cell (of its atoms without one)");
        if (changed) inst.scene.placement = placement_compose(offset, degrees, pivot);
        if (ImGui::Button("Reset placement")) inst.scene.placement = md_float4x4_ident();
        ImGui::SetItemTooltip("Where its coordinates put it, over the other instances");
    }

    ImGui::Separator();
    if (ImGui::MenuItem(ICON_FA_COPY " Duplicate instance", nullptr, false, sys.num_instances < INSTANCE_MAX)) {
        act->add_like = inst.id;
    }
    ImGui::SetItemTooltip("Another instance like it: of the same coordinates, with copies of its representations");
    if (ImGui::MenuItem(ICON_FA_XMARK " Remove instance", nullptr, false, sys.num_instances > 1)) {
        act->remove = inst.id;
    }
    ImGui::SetItemTooltip(first ? "It and its representations. The timeline and the script follow the next one." : "It and its representations");
}

bool draw_instance_header(SystemInstance& inst, InstanceActions* act) {
    SystemContext& sys = *inst.sys;
    ImGui::PushID((int)inst.id);
    ImGui::PushID("instance");

    char label[96];
    instance_label(label, sizeof(label), inst);
    const float btn = row_button_size();
    const float spacing = 2.0f;
    const float right_x = ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x;

    const bool hidden = !inst.scene.visible || !sys.visible;
    if (hidden) ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
    const ImGuiTreeNodeFlags flags = ImGuiTreeNodeFlags_SpanAvailWidth | ImGuiTreeNodeFlags_AllowOverlap | ImGuiTreeNodeFlags_DefaultOpen |
                                     ImGuiTreeNodeFlags_FramePadding | ImGuiTreeNodeFlags_NoTreePushOnOpen;
    const bool open = ImGui::TreeNodeEx("##instance", flags, ICON_FA_CLONE "  %s", label);
    if (hidden) ImGui::PopStyleColor();
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) {
        const bool first = &inst == &system_instance(sys);
        ImGui::SetTooltip("An instance of the system: its coordinates from %s, made by the system's operations,\n"
                          "placed and drawn with representations of its own.%s\n"
                          "Right click or " ICON_FA_ELLIPSIS_VERTICAL ": where its coordinates come from, its placement, its name.",
                          inst.run < 0 ? "the system's own" : inst.frame >= 0 ? "a frame of a run" : "a run, at the time shown",
                          first ? "\nThe first: the one the timeline, the script and the tools read." : "");
    }
    if (ImGui::BeginPopupContextItem("##instance_menu_rc")) {
        draw_instance_menu(inst, act);
        ImGui::EndPopup();
    }

    ImGui::SameLine(right_x - (btn + spacing) * 3, spacing);
    if (ImGui::Button(inst.scene.visible ? ICON_FA_EYE : ICON_FA_EYE_SLASH, {btn, btn})) {
        inst.scene.visible = !inst.scene.visible;
        sys.representation.atom_visibility_mask_dirty = true;
        if (inst.scene.visible) refresh_representations(sys);
    }
    ImGui::SetItemTooltip("Show/Hide the instance: all of its representations");
    ImGui::SameLine(0, spacing);
    if (ImGui::Button(ICON_FA_ELLIPSIS_VERTICAL, {btn, btn})) {
        ImGui::OpenPopup("##instance_menu");
    }
    ImGui::SetItemTooltip("Where its coordinates come from, its placement, its name");
    if (ImGui::BeginPopup("##instance_menu")) {
        draw_instance_menu(inst, act);
        ImGui::EndPopup();
    }
    ImGui::SameLine(0, spacing);
    ImGui::BeginDisabled(sys.num_instances == 1);
    if (ImGui::DeleteButton(ICON_FA_XMARK, {btn, btn})) {
        act->remove = inst.id;
    }
    ImGui::EndDisabled();
    ImGui::SetItemTooltip(sys.num_instances == 1 ? "A system has an instance at least" : "Remove the instance and its representations");

    ImGui::PopID();
    ImGui::PopID();
    return open;
}

void apply_instance_actions(SystemContext& sys, const InstanceActions& act) {
    if (act.remove) {
        for (size_t k = 0; k < sys.num_instances; ++k) {
            if (sys.instances[k]->id == act.remove) {
                instance_remove(sys, k);
                break;
            }
        }
    }
    if (act.add_like) {
        const SystemInstance* like = instance_find(sys, act.add_like);
        if (!like) return;
        SystemInstance* inst = instance_add(sys);
        if (!inst) return;
        // Of another run, where its coordinates put it; else as like is, where it is placed
        const bool other_run = act.add_run >= 0;
        instance_show(*inst, other_run ? act.add_run : like->run, other_run ? -1 : like->frame);
        inst->scene.visible = like->scene.visible;
        if (!other_run) inst->scene.placement = like->scene.placement;
        for (int i = 0; i < 3; ++i) inst->scene.periodic_images[i] = like->scene.periodic_images[i];
        instance_set_assembly(*inst, like->scene.assembly, false);
        representations_copy_to_instance(sys, *like, *inst);
    }
}

void draw_system_operations(SystemContext& sys) {
    ImGui::PushID("system");

    // One row, folded by default, that says which are on: all of them, opted into, under it
    char on[128] = "";
    int len = 0;
    bool inactive = false;
    for (const OperationInfo& info : operation_info) {
        if (!system_operation_enabled(sys, info.op)) continue;
        len += snprintf(on + len, sizeof(on) - len, "%s%s", len ? ", " : "", info.name);
        for (size_t k = 0; k < sys.num_instances; ++k) inactive |= info.moves_atoms && sys.instances[k]->scene.assembly >= 0;
    }

    const float btn = row_button_size();
    const ImVec2 min = ImGui::GetCursorScreenPos();
    const float icon_w = icon_column_width();
    ImGui::Dummy({icon_w, btn});
    ImGui::GetWindowDrawList()->AddText(ImFloor(min + (ImVec2(icon_w, btn) - ImGui::CalcTextSize(ICON_FA_SLIDERS)) * 0.5f),
                                        ImGui::GetColorU32(len ? ImGuiCol_Text : ImGuiCol_TextDisabled), ICON_FA_SLIDERS);
    ImGui::SameLine(0, 0);
    ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0, 3.0f));
    const bool open = ImGui::TreeNodeEx("##operations", ImGuiTreeNodeFlags_FramePadding | ImGuiTreeNodeFlags_AllowOverlap, "Operations");
    ImGui::PopStyleVar();
    ImGui::SetItemTooltip("What is done to the coordinates of the system before they are drawn, in every frame, in this order.\n"
                          "The script measures the coordinates as the files give them, not as these leave them.");
    ImGui::SameLine();
    if (len) {
        if (inactive) ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.2f, 1.0f), "%s", on);
        else          ImGui::TextDisabled("%s", on);
        if (inactive) ImGui::SetItemTooltip("Centering, wrapping and making whole are not applied while an assembly is shown");
    } else {
        ImGui::TextDisabled("none");
    }

    if (open) {
        for (const OperationInfo& info : operation_info) {
            draw_operation_row(sys, info);
        }
        ImGui::TreePop();
    }

    // What the system is shown as, of its one instance, when that is not listed on a row of its own
    if (!instance_rows_shown(sys)) draw_instance_copies_row(system_instance(sys));

    // What is done to the system and what it is shown as, apart from how it is drawn
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();
    ImGui::PopID();
}

// The state of the tracking of bond changes in the frame shown, and its advanced settings, under the operation
static void draw_bond_tracking_details(SystemContext& sys) {
    ApplicationState* app = sys.app;
    MD_ASSERT(app);
    auto& tracking = sys.operations.bond_tracking;

    const md_bond_delta_t* delta = bond_tracking_delta(sys);
    if (delta) {
        const size_t num_broken = md_bitset_count(delta->broken);
        ImGui::Text("%zu formed, %zu broken in this frame (%.1f ms)", delta->num_formed, num_broken, tracking.time_ms);
        if (ImGui::IsItemHovered()) {
            // The atoms of the bonds that changed
            md_bitmask_clear(&sys.selection.highlight_mask);
            for (size_t i = 0; i < delta->num_formed; ++i) {
                md_bitmask_set_bit(&sys.selection.highlight_mask, (uint32_t)delta->formed_pair[i].idx[0]);
                md_bitmask_set_bit(&sys.selection.highlight_mask, (uint32_t)delta->formed_pair[i].idx[1]);
            }
            md_bitset_iter_t it = md_bitset_iter(delta->broken);
            while (md_bitset_iter_next(&it)) {
                if (it.idx >= sys.data.sys.bond.count) break;
                const md_atom_pair_t p = sys.data.sys.bond.pairs[it.idx];
                md_bitmask_set_bit(&sys.selection.highlight_mask, (uint32_t)p.idx[0]);
                md_bitmask_set_bit(&sys.selection.highlight_mask, (uint32_t)p.idx[1]);
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
    const bool v2 = app->script.language == script::Language::V2;
    if (ImGui::SmallButton("Copy to Script")) {
        // The same as a script value, to count or plot. v1: formed(d), broken(d), bond(d) of a bond_delta();
        // VIAL: formed() and broken() of the scope, with the tracking's tolerance and breakable origins
        char args[512] = "";
        int len = 0;
        const md_str_t scope = md_str_trim(md_str_from_cstr(tracking.scope.text));
        if (!md_str_empty(scope)) len += snprintf(args + len, sizeof(args) - len, MD_STR_FMT, MD_STR_ARG(scope));
        if (tracking.tolerance != MD_BOND_DELTA_TOLERANCE && len < (int)sizeof(args)) {
            len += snprintf(args + len, sizeof(args) - len, "%stolerance=%.3g", len ? ", " : "", tracking.tolerance);
        }
        if (tracking.breakable != MD_BOND_DELTA_BREAKABLE && len < (int)sizeof(args)) {
            len += snprintf(args + len, sizeof(args) - len, "%sbreakable=%s", len ? ", " : "", v2 ? "[" : "{");
            bool first = true;
            for (int o = MD_BOND_ORIGIN_FILE; o <= MD_BOND_ORIGIN_INFERRED && len < (int)sizeof(args); ++o) {
                if (!(tracking.breakable & (1u << o))) continue;
                len += snprintf(args + len, sizeof(args) - len, "%s'%s'", first ? "" : ", ", md_bond_origin_name((md_bond_origin_t)o));
                first = false;
            }
            if (len < (int)sizeof(args)) len += snprintf(args + len, sizeof(args) - len, "%s", v2 ? "]" : "}");
        }
        char line[640];
        if (v2) {
            const md_str_t formed = create_unique_identifier(app->script.ir, MD_STR_LIT("formed"), app->allocator.frame);
            const md_str_t broken = create_unique_identifier(app->script.ir, MD_STR_LIT("broken"), app->allocator.frame);
            snprintf(line, sizeof(line), MD_STR_FMT " = formed(%s)", MD_STR_ARG(formed), args);
            script_editor::append_line(app->editor, md_str_from_cstr(line));
            snprintf(line, sizeof(line), MD_STR_FMT " = broken(%s)", MD_STR_ARG(broken), args);
            script_editor::append_line(app->editor, md_str_from_cstr(line));
        } else {
            const md_str_t ident = create_unique_identifier(app->script.ir, MD_STR_LIT("bonds_delta"), app->allocator.frame);
            snprintf(line, sizeof(line), MD_STR_FMT " = bond_delta(%s);", MD_STR_ARG(ident), args);
            script_editor::append_line(app->editor, md_str_from_cstr(line));
        }
    }
    if (v2) {
        ImGui::SetItemTooltip("Add the same bond changes to the script, per frame: formed() and broken() with this\n"
                              "scope, tolerance and breakable origins, to count them (len(...)) or select them");
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
// Files
// ---------------------------------------------------------------------------------------------------------

// What a file of the system is: a topology, coordinates (a restart), a structure, or the data of its atoms
static const char* system_file_kind(const FileEntry& e) {
    const LoaderFlags flags = loader::type_flags(e.type);
    if (flags & LoaderFlag_Topology) return "Topology";
    if (e.type == LoaderType_AMBER_RESTART || e.type == LoaderType_NAMD_COOR) return "Coordinates";
    if (flags & LoaderFlag_System) return "Structure";
    return e.type == LoaderType_NDX ? "Groups" : "Supplement";
}

const char* system_file_role_label(const FileEntry& e, FileRole role) {
    switch (role) {
    case FileRole_System:     return system_file_kind(e);
    case FileRole_Trajectory: return "Trajectory";
    case FileRole_Supplement: return system_file_kind(e);
    case FileRole_Series:     return (e.type == LoaderType_EDR || e.type == LoaderType_AMBER_ENERGY) ? "Energies" : "Series";
    default:                  return "";
    }
}

static const char* file_role_icon(const FileEntry& e, FileRole role) {
    switch (role) {
    case FileRole_System:     return ICON_FA_ATOM;
    case FileRole_Trajectory: return ICON_FA_FILM;
    case FileRole_Supplement: return (e.type == LoaderType_ITP || e.type == LoaderType_AMBER_PRMTOP || e.type == LoaderType_PSF) ? ICON_FA_DIAGRAM_PROJECT : ICON_FA_FILE_CIRCLE_PLUS;
    case FileRole_Series:     return ICON_FA_CHART_LINE;
    default:                  return ICON_FA_FILE;
    }
}

void system_file_details(char* buf, size_t cap, const SystemContext& sys, size_t i) {
    MD_ASSERT(buf && cap > 0);
    buf[0] = '\0';
    const FileList* list = &sys.files;
    const FileEntry* e = file_list_get(list, i);
    if (!e || e->missing) return;
    switch (file_list_role(list, i)) {
    case FileRole_System: {
        // A topology alone has no coordinates, until a file of them (or a run) joins it
        snprintf(buf, cap, "%zu atoms%s%s", sys.data.sys.atom.count, md_system_is_coarse_grained(&sys.data.sys) ? ", coarse grained" : "",
                 md_system_state_has_coords(&sys.data.state) ? "" : ", no coordinates");
        break;
    }
    case FileRole_Trajectory: {
        // Of its own run: run i of the list is run i of the system
        const SystemRun* run = (e->run >= 0 && (size_t)e->run < sys.num_runs) ? sys.runs[e->run] : nullptr;
        const size_t num_frames = run_num_frames(sys, run);
        const double* times = run_frame_times(sys, run);
        if (num_frames > 0 && times && !md_unit_is_none(run_time_unit(sys, run))) {
            char unit_buf[32];
            const double scl = display_units::factor_print(unit_buf, sizeof(unit_buf), run_time_unit(sys, run));
            snprintf(buf, cap, "%zu frames, %.4g - %.4g %s", num_frames, times[0] * scl, times[num_frames - 1] * scl, unit_buf);
        } else {
            snprintf(buf, cap, "%zu frames", num_frames);
        }
        break;
    }
    case FileRole_Series: {
        // The group the file published along its run, which names its source
        const SystemRun* run = (e->run >= 0 && (size_t)e->run < sys.num_runs) ? sys.runs[e->run] : nullptr;
        md_str_t groups[64];
        const size_t num_groups = run ? MD_MIN(system_series_groups(groups, MD_COUNTOF(groups), sys, md_str_from_cstr(run->path)), MD_COUNTOF(groups)) : 0;
        const md_str_t path = md_str_from_cstr(e->path);
        for (size_t g = 0; g < num_groups; ++g) {
            if (md_str_eq(system_series_group_source(sys, groups[g]), path)) {
                snprintf(buf, cap, "%zu series", system_series_members(nullptr, 0, sys, groups[g]));
                break;
            }
        }
        break;
    }
    case FileRole_Supplement:
        snprintf(buf, cap, MD_STR_FMT, MD_STR_ARG(loader::type_name(e->type)));
        break;
    default:
        break;
    }
}

// Where a file of the system may go: into another system, or one which can stand alone (a structure, a topology, a
// trajectory) into a new one, while something of the system stays (file_list_split)
static void draw_file_move_menu(SystemContext& sc, size_t i) {
    ApplicationState* app = sc.app;
    MD_ASSERT(app);
    if (!ImGui::BeginMenu("Move to")) return;
    bool any = false;
    for (size_t k = 0; k < app->num_systems; ++k) {
        const SystemContext& other = *app->systems[k];
        if (&other == &sc || !system_loaded(other)) continue;
        char label[128];
        system_label(label, sizeof(label), other);
        ImGui::PushID((int)other.id);
        if (ImGui::MenuItem(label)) app->file_move = { sc.id, i, other.id };
        ImGui::PopID();
        any = true;
    }
    const FileEntry& e = *file_list_get(&sc.files, i);
    const PlacementKind kind = drop_file_kind(loader::type_flags(e.type));
    const bool alone = kind == PlacementKind_System || kind == PlacementKind_Topology || kind == PlacementKind_Trajectory;
    FileList rest = {}, moved = {};
    file_list_split(&rest, &moved, &sc.files, i, app->allocator.frame);
    if (alone) {
        if (any) ImGui::Separator();
        ImGui::BeginDisabled(file_list_count(&rest) == 0);
        if (ImGui::MenuItem("A new system")) app->file_move = { sc.id, i, SYSTEM_ID_NONE };
        ImGui::EndDisabled();
        ImGui::SetItemTooltip("%s", file_list_count(&rest) ? "A system of its own, and this one without it" : "It is all of this system");
    } else if (!any) {
        ImGui::TextDisabled("No other system");
    }
    size_t num_along = 0;
    for (size_t j = 0; j < file_list_count(&moved); ++j) num_along += strcmp(file_list_get(&moved, j)->path, e.path) != 0;
    if (num_along) {
        ImGui::Separator();
        ImGui::TextDisabled("With it: the %zu other files of its run", num_along);
    }
    ImGui::EndMenu();
}

void draw_system_files(SystemContext& sc, InstanceActions* inst_act) {
    const SystemContext& sys = sc;
    const FileList* list = &sys.files;
    const size_t num_files = file_list_count(list);
    if (num_files == 0) return;

    ImGui::PushID("system_files");
    const float btn = row_button_size();
    const ImVec2 min = ImGui::GetCursorScreenPos();
    const float icon_w = icon_column_width();
    const float right_x = ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x;

    ImGui::Dummy({icon_w, btn});
    ImGui::GetWindowDrawList()->AddText(ImFloor(min + (ImVec2(icon_w, btn) - ImGui::CalcTextSize(ICON_FA_FOLDER_OPEN)) * 0.5f),
                                        ImGui::GetColorU32(ImGuiCol_Text), ICON_FA_FOLDER_OPEN);
    ImGui::SameLine(0, 0);

    size_t num_missing = 0;
    for (size_t i = 0; i < num_files; ++i) num_missing += file_list_get(list, i)->missing ? 1 : 0;

    ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0, 3.0f));
    const bool open = ImGui::TreeNodeEx("##files", ImGuiTreeNodeFlags_FramePadding | ImGuiTreeNodeFlags_AllowOverlap, "Files (%zu)", num_files);
    ImGui::PopStyleVar();
    ImGui::SetItemTooltip("The files the system is built from, in the order they were loaded: its own, then those of its run");
    if (num_missing) {
        ImGui::SameLine();
        ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.2f, 1.0f), ICON_FA_TRIANGLE_EXCLAMATION " %zu not loaded", num_missing);
    }

    // The System window has them in full
    if (bool* show = workspace_window_flag("System")) {
        ImGui::SameLine(right_x - btn);
        if (ImGui::Button(ICON_FA_CIRCLE_INFO, {btn, btn})) {
            *show = true;
        }
        ImGui::SetItemTooltip("Details: the System window");
    }

    if (open) {
        for (size_t i = 0; i < num_files; ++i) {
            const FileEntry& e = *file_list_get(list, i);
            const FileRole role = file_list_role(list, i);
            ImGui::PushID((int)i);

            // The files of the run under it, one further in
            const bool of_run = e.run >= 0;
            if (of_run && role != FileRole_Trajectory) ImGui::Indent();

            const md_str_t name = md_path_filename(md_str_from_cstr(e.path));
            ImGui::AlignTextToFramePadding();
            if (e.missing) {
                ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.2f, 1.0f), ICON_FA_TRIANGLE_EXCLAMATION " " MD_STR_FMT, MD_STR_ARG(name));
            } else {
                ImGui::Text("%s " MD_STR_FMT, file_role_icon(e, role), MD_STR_ARG(name));
            }
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) {
                if (e.missing) ImGui::SetTooltip("%s, not loaded: the workspace names it, but it could not be loaded\n%s", system_file_role_label(e, role), e.path);
                else           ImGui::SetTooltip("%s\n%s", system_file_role_label(e, role), e.path);
            }
            if (!e.missing && ImGui::BeginDragDropSource(ImGuiDragDropFlags_SourceAllowNullID)) {
                const SystemFileRef ref = { sc.id, (uint32_t)i };
                ImGui::SetDragDropPayload(SYSTEM_FILE_DND, &ref, sizeof(ref));
                ImGui::Text("%s " MD_STR_FMT, file_role_icon(e, role), MD_STR_ARG(name));
                ImGui::TextDisabled("Onto a system: into it. Below the systems: into a new one.");
                ImGui::EndDragDropSource();
            }
            if (ImGui::BeginPopupContextItem("##file")) {
                if (ImGui::MenuItem("Copy path")) ImGui::SetClipboardText(e.path);
                if (!e.missing) draw_file_move_menu(sc, i);
                if (role == FileRole_Trajectory && e.run >= 0 && (size_t)e.run < sc.num_runs) {
                    ImGui::Separator();
                    SystemInstance& first = system_instance(sc);
                    if (ImGui::MenuItem("Show this run", nullptr, e.run == first.run, !e.missing && e.run != first.run)) {
                        instance_show(first, e.run);
                    }
                    ImGui::SetItemTooltip("In the first instance of the system, which the timeline and the script follow");
                    if (ImGui::MenuItem("Show in a new instance", nullptr, false, !e.missing && sc.num_instances < INSTANCE_MAX)) {
                        inst_act->add_like = first.id;
                        inst_act->add_run  = e.run;
                    }
                    ImGui::SetItemTooltip("Beside the instances there are, superposed on them, drawn as the first is");
                    if (ImGui::MenuItem("Close this run")) sc.app->run_close = { sc.id, (size_t)e.run };
                    ImGui::SetItemTooltip("Its file and the series along it are unloaded, and the instances which show it\n"
                                          "with their representations, but for the first");
                    if (!e.missing && run_num_frames(sc, sc.runs[e.run]) > 0) {
                        ImGui::SeparatorText("On the time shown");
                        draw_run_time_offset(sc, sc.runs[e.run]);
                    }
                }
                ImGui::EndPopup();
            }
            char details[128];
            system_file_details(details, sizeof(details), sc, i);
            if (details[0]) {
                ImGui::SameLine();
                ImGui::TextDisabled("%s", details);
            }
            // Of several runs, the one the system is drawn at
            bool shown = false;
            for (size_t k = 0; k < sc.num_instances; ++k) shown |= sc.instances[k]->run == e.run && sc.instances[k]->scene.visible;
            if (role == FileRole_Trajectory && sc.num_runs > 1 && shown) {
                ImGui::SameLine();
                ImGui::TextColored(ImGui::GetStyleColorVec4(ImGuiCol_CheckMark), ICON_FA_EYE " shown");
            }

            if (of_run && role != FileRole_Trajectory) ImGui::Unindent();
            ImGui::PopID();
        }
        ImGui::TreePop();
    }
    ImGui::PopID();
}

void system_remove(SystemContext& sys) {
    ApplicationState* app = sys.app;
    MD_ASSERT(app);
    // The quantum chemistry is of another, when it was of this one and another has some
    qm_system_update(app, &sys);
    // Its data, when it has any: an empty system (the first, before anything is loaded) has nothing to free
    if (system_loaded(sys) || sys.num_runs > 0) {
        interrupt_async_tasks(app);
        free_trajectory_data(app, sys);
        free_system_data(app, sys);
    }

    // What was of it: its representations, its operations, its instances and what they were shown as. The first stays,
    // as a system has one.
    remove_all_representations(sys);
    while (sys.num_instances > 1) instance_remove(sys, sys.num_instances - 1);
    sys.visible = true;
    system_reset_operations(sys);
    SystemInstance& first = system_instance(sys);
    first.name[0] = '\0';
    first.scene.visible = true;
    first.scene.placement = md_float4x4_ident();
    first.scene.periodic_images[0] = first.scene.periodic_images[1] = first.scene.periodic_images[2] = 0;
    first.scene.assembly = -1;
    md_array_shrink(first.scene.assembly_transforms, 0);
    md_array_shrink(first.scene.assembly_copy_offset, 0);
    sys.representation.atom_visibility_mask_dirty = true;
    if (system_count(*app) == 0) app->animation.time = 0.0;
}

// ---------------------------------------------------------------------------------------------------------
// Workspace
// ---------------------------------------------------------------------------------------------------------

void system_scene_serialize(viamd::serialization_state_t& state, const SystemContext& sys) {
    const auto& op = sys.operations;

    viamd::write_section_header(state, MD_STR_LIT("SystemScene"));
    viamd::write_int(state, MD_STR_LIT("System"), (int)sys.id);
    viamd::write_str(state, MD_STR_LIT("Dir"), md_str_from_cstr(sys.dir));
    viamd::write_bool(state, MD_STR_LIT("Visible"), sys.visible);

    // The operations: on or off, and their settings
    viamd::write_bool(state, MD_STR_LIT("Recenter"), op.recenter);
    viamd::write_bool(state, MD_STR_LIT("FixateOrientation"), op.fixate_orientation);
    viamd::write_bool(state, MD_STR_LIT("RecenterQueryEnabled"), op.recenter_query.enabled);
    viamd::write_str(state,  MD_STR_LIT("RecenterQuery"), md_str_from_cstr(op.recenter_query.filter.text));
    if (!md_bitset_empty(op.selection_mask)) {
        viamd::write_bitset(state, MD_STR_LIT("RecenterTarget"), &op.selection_mask);
    }
    viamd::write_bool(state, MD_STR_LIT("ApplyPbc"), op.apply_pbc);
    viamd::write_bool(state, MD_STR_LIT("UnwrapStructures"), op.unwrap_structures);
    viamd::write_bool(state, MD_STR_LIT("TrackBonds"), op.bond_tracking.enabled);
    viamd::write_str(state,  MD_STR_LIT("TrackBondsScope"), md_str_from_cstr(op.bond_tracking.scope.text));
    viamd::write_flt(state,  MD_STR_LIT("TrackBondsTolerance"), op.bond_tracking.tolerance);
    char origins[64];
    viamd::write_str(state,  MD_STR_LIT("TrackBondsBreakable"), md_str_from_cstr(bond_origins_str(origins, sizeof(origins), op.bond_tracking.breakable)));

    // Its instances, the first first
    for (size_t k = 0; k < sys.num_instances; ++k) {
        const SystemInstance& inst = *sys.instances[k];
        const SystemScene& s = inst.scene;
        viamd::write_section_header(state, MD_STR_LIT("Instance"));
        viamd::write_int(state, MD_STR_LIT("System"), (int)sys.id);
        viamd::write_int(state, MD_STR_LIT("Id"), (int)inst.id);
        if (inst.name[0]) viamd::write_str(state, MD_STR_LIT("Name"), md_str_from_cstr(inst.name));
        // The run by its index among the system's, which is that of the workspace's [Run]s of it
        viamd::write_int(state, MD_STR_LIT("Run"), inst.run);
        if (inst.run >= 0 && inst.frame >= 0) viamd::write_int(state, MD_STR_LIT("Frame"), (int)inst.frame);
        viamd::write_bool(state, MD_STR_LIT("Visible"), s.visible);
        if (!md_float4x4_equal(s.placement, md_float4x4_ident())) {
            // Column major
            viamd::write_flt_vec(state, MD_STR_LIT("Placement"), &s.placement.elem[0][0], 16);
        }
        if ((s.periodic_images[0] | s.periodic_images[1] | s.periodic_images[2]) != 0) {
            viamd::write_int_vec(state, MD_STR_LIT("PeriodicImages"), s.periodic_images, 3);
        }
        if (s.assembly >= 0 && (size_t)s.assembly < sys.data.sys.assembly.count) {
            viamd::write_str(state, MD_STR_LIT("Assembly"), md_assembly_id(&sys.data.sys.assembly, (size_t)s.assembly));
        }
    }
}

bool system_scene_read_instance_entry(SystemScene& s, md_str_t ident, md_str_t arg) {
    if (md_str_eq(ident, MD_STR_LIT("Visible"))) {
        viamd::extract_bool(s.visible, arg);
    } else if (md_str_eq(ident, MD_STR_LIT("Placement"))) {
        // Column major
        md_float4x4_t M;
        if (viamd::extract_flt_vec(&M.elem[0][0], 16, arg)) s.placement = M;
    } else if (md_str_eq(ident, MD_STR_LIT("PeriodicImages"))) {
        // Integers; a workspace written before [SystemScene] had them as floats ([RenderSettings])
        float n[3];
        if (viamd::extract_flt_vec(n, 3, arg)) {
            for (int k = 0; k < 3; ++k) s.periodic_images[k] = MD_CLAMP((int)n[k], 0, PERIODIC_IMAGES_MAX);
        }
    } else {
        return false;
    }
    return true;
}

bool system_scene_read_entry(SystemContext& sys, md_str_t ident, md_str_t arg) {
    auto& op = sys.operations;

    if (md_str_eq(ident, MD_STR_LIT("Visible"))) {
        viamd::extract_bool(sys.visible, arg);
    } else if (md_str_eq(ident, MD_STR_LIT("Dir"))) {
        md_str_t dir = {};
        if (viamd::extract_str(dir, arg)) system_dir_take(sys, dir);
    } else if (md_str_eq(ident, MD_STR_LIT("Placement")) || md_str_eq(ident, MD_STR_LIT("PeriodicImages"))) {
        // Of a workspace written before the instances ([Instance]): the system's, now its first instance's
        return system_scene_read_instance_entry(system_instance(sys).scene, ident, arg);
    } else if (md_str_eq(ident, MD_STR_LIT("Operations"))) {
        // The operations a system had, of a workspace written while they were added and removed: every system has
        // them all now, on or off by their own entries
    } else if (md_str_eq(ident, MD_STR_LIT("Recenter"))) {
        viamd::extract_bool(op.recenter, arg);
    } else if (md_str_eq(ident, MD_STR_LIT("FixateOrientation"))) {
        viamd::extract_bool(op.fixate_orientation, arg);
    } else if (md_str_eq(ident, MD_STR_LIT("ApplyPbc"))) {
        viamd::extract_bool(op.apply_pbc, arg);
    } else if (md_str_eq(ident, MD_STR_LIT("UnwrapStructures"))) {
        viamd::extract_bool(op.unwrap_structures, arg);
    } else if (md_str_eq(ident, MD_STR_LIT("RecalcBonds")) || md_str_eq(ident, MD_STR_LIT("TrackBonds"))) {
        // RecalcBonds: what a workspace from before the bond changes were tracked had, the bonds of every frame guessed anew
        viamd::extract_bool(op.bond_tracking.enabled, arg);
    } else if (md_str_eq(ident, MD_STR_LIT("TrackBondsScope"))) {
        viamd::extract_to_char_buf(op.bond_tracking.scope.text, sizeof(op.bond_tracking.scope.text), arg);
        op.bond_tracking.scope_dirty = true;
    } else if (md_str_eq(ident, MD_STR_LIT("TrackBondsTolerance"))) {
        float tolerance = MD_BOND_DELTA_TOLERANCE;
        if (viamd::extract_flt(tolerance, arg)) op.bond_tracking.tolerance = MD_CLAMP(tolerance, 0.0f, 0.45f);
    } else if (md_str_eq(ident, MD_STR_LIT("TrackBondsBreakable"))) {
        md_str_t str;
        uint32_t mask = 0;
        if (viamd::extract_str(str, arg) && bond_origins_parse(&mask, str)) op.bond_tracking.breakable = mask;
    } else if (md_str_eq(ident, MD_STR_LIT("RecenterQueryEnabled"))) {
        viamd::extract_bool(op.recenter_query.enabled, arg);
    } else if (md_str_eq(ident, MD_STR_LIT("RecenterQuery"))) {
        viamd::extract_to_char_buf(op.recenter_query.filter.text, sizeof(op.recenter_query.filter.text), arg);
        recenter_mark_query_dirty(sys);
    } else {
        return false;
    }
    return true;
}
