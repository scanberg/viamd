#include <color_scale.h>

#include <display_units.h>
#include <imgui_widgets.h>

#include <core/md_common.h>

#include <imgui.h>
#include <implot.h>
#include <implot_internal.h>   // SampleColormapU32

#include <math.h>
#include <stdio.h>

uint32_t color_scale_color_u32(const ColorScale& scale, float value) {
    return ImPlot::SampleColormapU32(color_scale_param(scale, value), scale.colormap);
}

bool color_scale_draw_controls(ColorScale* scale, const ColorScaleSpan& span, md_unit_t unit, const char* span_label, const char* span_tooltip) {
    ASSERT(scale);
    bool changed = false;

    ImGui::PushID("color_scale");

    const ImVec2 button_size = {ImGui::CalcItemWidth(), 0};
    if (ImPlot::ColormapButton(ImPlot::GetColormapName(scale->colormap), button_size, scale->colormap)) {
        ImGui::OpenPopup("Colormap");
    }
    ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
    changed |= ImGui::Checkbox("legend", &scale->show_legend);
    ImGui::SetItemTooltip("Show the colour map with its range and unit in the view.\nHold Alt to move or resize it.");
    if (ImGui::BeginPopup("Colormap")) {
        // The qualitative maps are left out: their colours do not order, so they cannot carry a value
        for (int m = COLOR_SCALE_FIRST_CONTINUOUS_COLORMAP; m < ImPlot::GetColormapCount(); ++m) {
            if (ImPlot::ColormapButton(ImPlot::GetColormapName(m), button_size, m)) {
                scale->colormap = m;
                changed = true;
                ImGui::CloseCurrentPopup();
            }
        }
        ImGui::EndPopup();
    }

    // Shown in the display unit, kept in the values' own
    char unit_str[32] = "";
    const double scl = display_units::factor_print(unit_str, sizeof(unit_str), unit);
    const float fscl = (float)scl;

    bool mode_changed = false;
    mode_changed |= ImGui::Checkbox("auto range", &scale->auto_range);
    if (span_tooltip && span_tooltip[0]) {
        ImGui::SetItemTooltip("Follow the span of the values: %s", span_tooltip);
    } else {
        ImGui::SetItemTooltip("Follow the span of the values");
    }
    ImGui::SameLine();
    mode_changed |= ImGui::Checkbox("symmetric", &scale->symmetric);
    ImGui::SetItemTooltip("A range centred on zero");
    if (mode_changed) {
        if (scale->symmetric && !scale->auto_range) {
            const float r = MAX(fabsf(scale->range_beg), fabsf(scale->range_end));
            scale->range_beg = -r;
            scale->range_end =  r;
        }
        color_scale_update_range(scale, span);
        changed = true;
    }

    // The unit goes into a format string, so a '%' in it would be read as a conversion
    char unit_fmt[64] = "";
    for (size_t i = 0, j = 0; unit_str[i] && j + 2 < sizeof(unit_fmt); ++i) {
        if (unit_str[i] == '%') unit_fmt[j++] = '%';
        unit_fmt[j++] = unit_str[i];
        unit_fmt[j] = '\0';
    }
    char fmt[96];
    if (scale->symmetric) {
        snprintf(fmt, sizeof(fmt), "%%.4g %s", unit_fmt);
    } else {
        snprintf(fmt, sizeof(fmt), "%%.4g to %%.4g %s", unit_fmt);
    }
    if (scale->symmetric) {
        // From zero to half as far again as the largest magnitude, and never short of the range itself
        float r = MAX(fabsf(scale->range_beg), fabsf(scale->range_end));
        float limit = span.valid ? 1.5f * MAX(fabsf(span.min), fabsf(span.max)) : 0.0f;
        limit = MAX(limit, r);
        if (!(limit > 0.0f)) limit = 1.0f;
        float r_disp = r * fscl;
        if (ImGui::SliderFloat("range ±", &r_disp, 0.0f, limit * fscl, fmt)) {
            scale->range_beg = -r_disp / fscl;
            scale->range_end =  r_disp / fscl;
            scale->auto_range = false;
            changed = true;
        }
    } else {
        // The span padded by half its width on either side, and never short of the range itself
        float lo = scale->range_beg, hi = scale->range_end;
        if (span.valid) {
            float pad = 0.5f * (span.max - span.min);
            if (!(pad > 0.0f)) pad = MAX(0.5f * fabsf(span.max), 1.0f);
            lo = MIN(lo, span.min - pad);
            hi = MAX(hi, span.max + pad);
        }
        if (!(hi > lo)) { lo -= 1.0f; hi += 1.0f; }
        float beg = scale->range_beg * fscl;
        float end = scale->range_end * fscl;
        if (ImGui::RangeSliderFloat("range", &beg, &end, lo * fscl, hi * fscl, fmt)) {
            scale->range_beg = beg / fscl;
            scale->range_end = end / fscl;
            scale->auto_range = false;
            changed = true;
        }
    }

    if (span.valid) {
        ImGui::TextDisabled("%s: %.4g to %.4g %s", span_label ? span_label : "values", span.lo * scl, span.hi * scl, unit_str);
        if (span.lo != span.min || span.hi != span.max) {
            ImGui::SetItemTooltip("%s%sextremes %.4g to %.4g %s", span_tooltip ? span_tooltip : "", span_tooltip ? "; " : "", span.min * scl, span.max * scl, unit_str);
        } else if (span_tooltip && span_tooltip[0]) {
            ImGui::SetItemTooltip("%s", span_tooltip);
        }
    }

    ImGui::PopID();
    return changed;
}

void color_scale_draw_legend(const ColorScale& scale, md_unit_t unit, const char* label, const char* owner, int id, int slot) {
    if (scale.range_end == scale.range_beg) return;

    char unit_str[32] = "";
    const double scl = display_units::factor_print(unit_str, sizeof(unit_str), unit);
    char header[128];
    if (unit_str[0]) snprintf(header, sizeof(header), "%s (%s)", label ? label : "", unit_str);
    else             snprintf(header, sizeof(header), "%s", label ? label : "");

    const ImGuiStyle& style = ImGui::GetStyle();
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    const float width = MAX(ImGui::CalcTextSize(header).x + 2.0f * style.WindowPadding.x, 120.0f);
    const float height = 300.0f;
    // Where a legend starts out, the first time: down the right edge of the view in the order the
    // legends are drawn, a new column to the left of the last when the view is full
    const float top = 60.0f, gap = 10.0f;
    const int per_column = MAX(1, (int)((viewport->WorkSize.y - top) / (height + gap)));
    const int column = MAX(slot, 0) / per_column, row = MAX(slot, 0) % per_column;
    ImGui::SetNextWindowPos(ImVec2(viewport->WorkPos.x + viewport->WorkSize.x - (float)(column + 1) * (width + 20.0f),
                                   viewport->WorkPos.y + top + (float)row * (height + gap)), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(width, height), ImGuiCond_FirstUseEver);
    // Never narrower than its heading: the size is the user's and kept in the .ini, and what it heads
    // changes - another attribute, another unit - without the user having resized anything
    ImGui::SetNextWindowSizeConstraints(ImVec2(MAX(width, 80.0f), 120), ImVec2(MAX(width, 1000.0f), 2000));

    const bool editable = ImGui::IsKeyDown(ImGuiMod_Alt);
    ImGuiWindowFlags flags = ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoCollapse |
                             ImGuiWindowFlags_NoFocusOnAppearing;
    if (!editable) {
        flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoNavFocus | ImGuiWindowFlags_NoInputs;
    }

    // The look of the text the script draws into the view: light on a translucent dark backing,
    // which reads on any background and any colour map
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.0f, 0.0f, 0.0f, editable ? 0.65f : 0.5f));
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 1.0f, 1.0f, 1.0f));
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 5.0f);
    ImPlot::PushStyleColor(ImPlotCol_FrameBg, ImVec4(0, 0, 0, 0));

    char title[160];
    snprintf(title, sizeof(title), "Legend: %s###legend_%d", owner ? owner : "", id);
    if (ImGui::Begin(title, nullptr, flags)) {
        // Kept within the view: it grows to the right to fit a longer heading, and the view can
        // shrink under a position that was kept from a larger one
        const ImVec2 pos = ImGui::GetWindowPos(), size = ImGui::GetWindowSize();
        const ImVec2 lo = viewport->WorkPos;
        const ImVec2 hi = ImVec2(viewport->WorkPos.x + viewport->WorkSize.x - size.x, viewport->WorkPos.y + viewport->WorkSize.y - size.y);
        const ImVec2 clamped = ImVec2(MAX(lo.x, MIN(pos.x, hi.x)), MAX(lo.y, MIN(pos.y, hi.y)));
        if (clamped.x != pos.x || clamped.y != pos.y) {
            ImGui::SetWindowPos(clamped);
        }

        ImGui::TextUnformatted(header);
        const ImVec2 avail = ImGui::GetContentRegionAvail();
        char scale_id[32];
        snprintf(scale_id, sizeof(scale_id), "##scale%d", id);
        ImPlot::ColormapScale(scale_id, scale.range_beg * scl, scale.range_end * scl, ImVec2(avail.x, avail.y), "%g", ImPlotColormapScaleFlags_NoLabel, scale.colormap);
    }
    ImGui::End();

    ImPlot::PopStyleColor();
    ImGui::PopStyleVar(2);
    ImGui::PopStyleColor(2);
}
