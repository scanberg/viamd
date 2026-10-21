#include <property_grid.h>

#include <core/md_common.h>

#include <imgui_internal.h>

#include <stdarg.h>
#include <stdio.h>

namespace prop {

// The panel being drawn into: one at a time, a grid can have many in a frame
static struct {
    bool    active = false;
    ImGuiID key = 0;            // of the grid, in the window's state storage
    float   base_x = 0.0f;      // window local x the panel starts at, which the labels are measured from
    float   right_x = 0.0f;     // window local x of the right edge of the panel
    float   label_end = 0.0f;   // window local x where the widest label of the grid ended, the frame before
    int     depth = 0;          // sections open
} g;

// The storage of a grid, by its key: [0] the frame [2] was last restarted in, [1] where the widest label of the frame
// before ended, which this frame lays out by, [2] where the widest label of this frame so far ends. Window local x,
// so that panels at different depths of indentation share the column.
static ImGuiID slot(int i) {
    return ImHashData(&i, sizeof(i), g.key);
}

static float indent_width() {
    return ImGui::GetTreeNodeToLabelSpacing();
}

static float label_gap() {
    return 2.0f * ImGui::GetStyle().ItemInnerSpacing.x;
}

void begin(const char* grid_name) {
    IM_ASSERT(!g.active && "prop::begin without prop::end");
    ImGuiStorage* st = ImGui::GetStateStorage();
    g.active = true;
    g.key    = ImHashStr(grid_name);    // by the name alone, not the ID stack: one column for the whole window
    g.depth  = 0;
    g.base_x = ImGui::GetCursorPosX();
    g.right_x = g.base_x + ImGui::GetContentRegionAvail().x;

    const int frame = ImGui::GetFrameCount();
    if (st->GetInt(slot(0), -1) != frame) {
        const float widest = st->GetFloat(slot(2), 0.0f);
        if (widest > 0.0f) st->SetFloat(slot(1), widest);
        st->SetFloat(slot(2), 0.0f);
        st->SetInt(slot(0), frame);
    }
    // Before anything has been measured: as wide as a short label
    g.label_end = st->GetFloat(slot(1), 0.0f);
    if (!(g.label_end > 0.0f)) g.label_end = g.base_x + ImGui::CalcTextSize("coloring").x;
}

void end() {
    IM_ASSERT(g.active && g.depth == 0 && "prop::end with a section open");
    g.active = false;
}

float control_x() {
    // At most a little under half the panel to the labels: a long label is cut rather than the controls squeezed
    const float max_end = g.base_x + 0.45f * MAX(g.right_x - g.base_x, 0.0f);
    return MIN(g.label_end, max_end) + label_gap();
}

float control_width() {
    return MAX(g.right_x - control_x(), 1.0f);
}

float control_width_before(const char* text) {
    const float w = (text && text[0]) ? ImGui::CalcTextSize(text).x + ImGui::GetStyle().ItemInnerSpacing.x * 2.0f : 0.0f;
    return MAX(control_width() - w, 1.0f);
}

// The header of a section: dimmed, unframed, with a rule from its name to the right edge - a divider that folds,
// beneath the rows of whatever the grid is the properties of. With 'enabled', the name is a checkbox that turns
// what the section is about on and off ('mixed': shown as neither, for what is on in part). As high as a row, so that
// the checkbox fits it.
static bool header(const char* name, bool default_open, bool* enabled, bool* toggled, const char* tooltip, bool mixed = false) {
    const ImGuiStyle& style = ImGui::GetStyle();
    ImGuiTreeNodeFlags flags = ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_SpanAvailWidth;
    if (default_open) flags |= ImGuiTreeNodeFlags_DefaultOpen;
    if (enabled)      flags |= ImGuiTreeNodeFlags_AllowOverlap;     // the checkbox over the node takes the clicks on it

    ImGui::AlignTextToFramePadding();
    ImGui::PushStyleColor(ImGuiCol_Text, style.Colors[ImGuiCol_TextDisabled]);
    const bool open = ImGui::TreeNodeEx("##section", flags, "%s", enabled ? "" : name);
    const ImVec2 mn = ImGui::GetItemRectMin();
    const ImVec2 mx = ImGui::GetItemRectMax();
    float x0 = mn.x + ImGui::GetTreeNodeToLabelSpacing() + ImGui::CalcTextSize(name).x + style.ItemInnerSpacing.x * 2.0f;
    if (enabled) {
        ImGui::SameLine();
        ImGui::SetCursorScreenPos(ImVec2(mn.x + ImGui::GetTreeNodeToLabelSpacing(), mn.y));
        if (mixed) ImGui::PushItemFlag(ImGuiItemFlags_MixedValue, true);
        const bool t = ImGui::Checkbox(name, enabled);
        if (mixed) ImGui::PopItemFlag();
        if (toggled) *toggled = t;
        if (tooltip && tooltip[0]) ImGui::SetItemTooltip("%s", tooltip);
        x0 = ImGui::GetItemRectMax().x + style.ItemInnerSpacing.x * 2.0f;
    }
    ImGui::PopStyleColor();

    const float x1 = ImGui::GetWindowPos().x - ImGui::GetScrollX() + g.right_x;
    const float y  = IM_ROUND((mn.y + mx.y) * 0.5f);
    if (x0 < x1) {
        ImGui::GetWindowDrawList()->AddLine(ImVec2(x0, y), ImVec2(x1, y), ImGui::GetColorU32(ImGuiCol_Separator));
    }
    return open;
}

static void section_open() {
    ImGui::Indent(indent_width());
    g.depth += 1;
}

bool section(const char* name, bool default_open) {
    IM_ASSERT(g.active);
    ImGui::PushID(name);
    if (!header(name, default_open, nullptr, nullptr, nullptr)) {
        ImGui::PopID();
        return false;
    }
    section_open();
    return true;
}

bool section_toggle(const char* name, bool* enabled, bool* toggled, bool default_open, const char* tooltip) {
    IM_ASSERT(g.active && enabled);
    if (toggled) *toggled = false;
    ImGui::PushID(name);
    const bool open = header(name, default_open, enabled, toggled, tooltip);
    if (!open || !*enabled) {
        ImGui::PopID();
        return false;
    }
    section_open();
    return true;
}

bool section_toggle_all(const char* name, bool* const* enabled, int count, bool* toggled, bool default_open, const char* tooltip) {
    IM_ASSERT(g.active && enabled && count > 0);
    if (toggled) *toggled = false;
    int on = 0;
    for (int i = 0; i < count; ++i) on += *enabled[i] ? 1 : 0;
    bool all = on == count;
    bool clicked = false;
    ImGui::PushID(name);
    const bool open = header(name, default_open, &all, &clicked, tooltip, on > 0 && on < count);
    if (clicked) {
        for (int i = 0; i < count; ++i) *enabled[i] = on == 0;
        if (toggled) *toggled = true;
    }
    if (!open) {
        ImGui::PopID();
        return false;
    }
    section_open();
    return true;
}

void section_end() {
    IM_ASSERT(g.active && g.depth > 0);
    ImGui::Unindent(indent_width());
    ImGui::PopID();
    g.depth -= 1;
}

void row(const char* label, const char* tooltip) {
    IM_ASSERT(g.active);
    const float x0 = ImGui::GetCursorPosX();
    const float cx = control_x();

    ImGui::AlignTextToFramePadding();
    const char* text_end = label ? ImGui::FindRenderedTextEnd(label) : nullptr;
    if (label && text_end > label) {
        // Measured for the next frame's column, where it ends in the window
        ImGuiStorage* st = ImGui::GetStateStorage();
        const float w = x0 + ImGui::CalcTextSize(label, text_end).x;
        if (w > st->GetFloat(slot(2), 0.0f)) st->SetFloat(slot(2), w);

        // Cut at the control column, should it be narrower than the label
        const ImVec2 pos = ImGui::GetCursorScreenPos();
        const float clip_x = ImGui::GetWindowPos().x - ImGui::GetScrollX() + cx - 0.5f * label_gap();
        ImGui::PushClipRect(ImVec2(pos.x, pos.y), ImVec2(clip_x, pos.y + ImGui::GetFrameHeight()), true);
        ImVec4 col = ImGui::GetStyleColorVec4(ImGuiCol_Text);
        col.w *= 0.7f;
        ImGui::TextColored(col, "%.*s", (int)(text_end - label), label);
        ImGui::PopClipRect();
        if (tooltip && tooltip[0] && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("%s", tooltip);
        }
    } else {
        ImGui::TextUnformatted("");
    }
    ImGui::SameLine(cx);
    ImGui::SetNextItemWidth(-FLT_MIN);
}

void hint(const char* fmt, ...) {
    row("");
    char buf[512];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    ImGui::PushTextWrapPos(0.0f);
    ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyle().Colors[ImGuiCol_TextDisabled]);
    ImGui::TextUnformatted(buf);
    ImGui::PopStyleColor();
    ImGui::PopTextWrapPos();
}

bool color(const char* id, float* c, bool alpha) {
    const float w = ImGui::CalcItemWidth();
    ImGui::PushID(id);
    const ImGuiColorEditFlags flags = alpha ? (ImGuiColorEditFlags_AlphaPreviewHalf | ImGuiColorEditFlags_AlphaBar) : ImGuiColorEditFlags_NoAlpha;
    const ImVec4 col(c[0], c[1], c[2], alpha ? c[3] : 1.0f);
    if (ImGui::ColorButton("##swatch", col, flags, ImVec2(w, ImGui::GetFrameHeight()))) {
        ImGui::OpenPopup("##picker");
    }
    bool changed = false;
    if (ImGui::BeginPopup("##picker")) {
        changed = alpha ? ImGui::ColorPicker4("##color", c, flags) : ImGui::ColorPicker3("##color", c, flags);
        ImGui::EndPopup();
    }
    ImGui::PopID();
    return changed;
}

}  // namespace prop
