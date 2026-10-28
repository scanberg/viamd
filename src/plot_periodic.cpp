#include <plot_periodic.h>

#include <serialization_utils.h>

#include <implot_internal.h>

#include <core/md_common.h>
#include <core/md_str.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

// Few: one per period and unit shown. Never moved, so that axes can link to them.
#define MAX_PERIODIC_VIEWS 32
static PeriodicView views[MAX_PERIODIC_VIEWS];

PeriodicView* periodic_view(const SeriesPeriodic& p, const char* unit) {
    if (!series_periodic_valid(p)) return nullptr;
    unit = unit ? unit : "";
    PeriodicView* free_view = nullptr;
    for (int i = 0; i < MAX_PERIODIC_VIEWS; ++i) {
        PeriodicView& v = views[i];
        if (!v.used) {
            if (!free_view) free_view = &v;
            continue;
        }
        if (series_periodic_equal(v.period, p) && strcmp(v.unit, unit) == 0) {
            return &v;
        }
    }
    // Full: the oldest is reused, which only ever happens to a view nothing has shown in a long while
    PeriodicView* v = free_view ? free_view : &views[0];
    *v = PeriodicView{};
    v->used = true;
    v->period = p;
    snprintf(v->unit, sizeof(v->unit), "%s", unit);
    v->min = p.seam;
    v->max = p.seam + series_periodic_length(p);
    return v;
}

void periodic_views_clear() {
    for (int i = 0; i < MAX_PERIODIC_VIEWS; ++i) {
        views[i] = PeriodicView{};
    }
}

void periodic_view_cut(PeriodicView* view, double cut) {
    if (!view) return;
    view->pending_min = cut;
    view->has_pending = true;
}

void periodic_view_cut_least_populated(PeriodicView* view, const float* density, int num_bins) {
    if (!view) return;
    const int best = periodic_least_populated(density, num_bins);
    if (best < 0) return;
    const double len = periodic_view_length(*view);
    periodic_view_cut(view, view->period.beg + (best + 0.5) * len / num_bins);
}

static int periodic_format(double value, char* buf, int size, void* user_data) {
    const PeriodicView* view = (const PeriodicView*)user_data;
    double v = periodic_wrap(value, view->period.beg, periodic_view_length(*view));
    // Within rounding of a whole period: as the beginning of it, not 1e-14 short of the end
    const double len = periodic_view_length(*view);
    if (fabs(v - view->period.end) < 1.0e-9 * len) v = view->period.beg;
    if (fabs(v) < 1.0e-9 * len) v = 0.0;
    return snprintf(buf, size, "%g", v);
}

void periodic_axis_setup(ImAxis axis, const char* label, PeriodicView* view, ImPlotAxisFlags flags) {
    ASSERT(view);
    const double len = periodic_view_length(*view);
    double span = view->max - view->min;
    if (!(span > 0.0) || span > len) span = len;

    if (view->has_pending) {
        view->min = view->pending_min;
        view->max = view->pending_min + span;
        view->has_pending = false;
    }
    // Turned around by whole periods, the cut stays in the period the values are given in, which nothing can see: the
    // values and the labels wrap alike. Not while the mouse is down, when ImPlot may be panning from where it was.
    if (!ImGui::IsMouseDown(ImGuiMouseButton_Left) && !ImGui::IsMouseDown(ImGuiMouseButton_Right) && !ImGui::IsMouseDown(ImGuiMouseButton_Middle)) {
        const double shift = view->min - periodic_wrap(view->min, view->period.beg, len);
        if (shift != 0.0) {
            view->min -= shift;
            view->max -= shift;
        }
    }

    ImPlot::SetupAxis(axis, label, flags | ImPlotAxisFlags_NoMenus);
    ImPlot::SetupAxisLinks(axis, &view->min, &view->max);
    ImPlot::SetupAxisZoomConstraints(axis, len * 1.0e-4, len);
    ImPlot::SetupAxisFormat(axis, periodic_format, view);

    // Ticks a whole fraction of the period apart, counted from where the values are given: ImPlot's own, at 50 apart,
    // would read 150, -160, -110 around the end of a dihedral's period. As far apart as the plot had room for last frame.
    const ImPlotPlot* plot = ImPlot::GetCurrentPlot();
    const bool vertical = axis >= ImAxis_Y1;
    float pixels = plot ? (vertical ? plot->PlotRect.GetHeight() : plot->PlotRect.GetWidth()) : 0.0f;
    if (pixels <= 1.0f) pixels = 300.0f;
    const double min_spacing = ImGui::GetFontSize() * (vertical ? 2.5 : 4.0);
    // The divisors of 360: a step of a whole number of degrees, and of a period of any other length the same fractions
    static const int divisions[] = {360, 180, 120, 90, 72, 60, 45, 40, 36, 30, 24, 20, 18, 12, 10, 8, 6, 4, 3, 2, 1};
    const double view_span = view->max - view->min;
    double step = len;
    for (int d : divisions) {
        const double s = len / d;
        if (s / view_span * pixels >= min_spacing) {
            step = s;
            break;
        }
    }
    const double first = view->period.beg + ceil ((view->min - view->period.beg) / step - 1.0e-9) * step;
    const double last  = view->period.beg + floor((view->max - view->period.beg) / step + 1.0e-9) * step;
    const int num_ticks = (int)floor((last - first) / step + 0.5) + 1;
    if (num_ticks >= 2) {
        ImPlot::SetupAxisTicks(axis, first, last, num_ticks);
    }
}

bool periodic_axis_popup(ImAxis axis, PeriodicView* view) {
    ASSERT(view);
    char id[48];
    snprintf(id, sizeof(id), "##periodic_axis_%d", (int)axis);
    if (ImPlot::IsAxisHovered(axis) && ImGui::IsMouseReleased(ImGuiMouseButton_Right) && ImGui::GetMouseDragDelta(ImGuiMouseButton_Right).x == 0.0f) {
        ImGui::OpenPopup(id);
    }
    if (ImPlot::IsAxisHovered(axis)) {
        ImGui::SetItemTooltip("Periodic: drag to turn it, right-click for where to cut it");
    }

    bool least = false;
    if (ImGui::BeginPopup(id)) {
        const SeriesPeriodic& p = view->period;
        const double len = periodic_view_length(*view);
        const char* u = view->unit;
        ImGui::TextDisabled("Periodic in [%g, %g[ %s", p.beg, p.end, u);
        ImGui::SeparatorText("Cut at");

        char label[96];
        auto option = [&](double cut, const char* what) {
            snprintf(label, sizeof(label), "%g %s  (%s)", periodic_wrap(cut, p.beg, len), u, what);
            const bool on = fabs(periodic_wrap(view->min - cut, -0.5 * len, len)) < 1.0e-6 * len;
            if (ImGui::MenuItem(label, nullptr, on)) periodic_view_cut(view, cut);
        };
        option(p.seam, "where values of this kind are least likely");
        if (fabs(p.seam - p.beg) > 1.0e-6 * len) {
            option(p.beg, "as the values are given");
        }
        // The other common ways of writing an angle: about 0, and from 0
        if (fabs(periodic_wrap(-0.5 * len - p.beg, -0.5 * len, len)) > 1.0e-6 * len && fabs(periodic_wrap(-0.5 * len - p.seam, -0.5 * len, len)) > 1.0e-6 * len) {
            option(-0.5 * len, "about 0");
        }
        if (fabs(periodic_wrap(p.beg, -0.5 * len, len)) > 1.0e-6 * len && fabs(periodic_wrap(p.seam, -0.5 * len, len)) > 1.0e-6 * len) {
            option(0.0, "from 0");
        }
        if (ImGui::MenuItem("Where the values shown are least populated")) {
            least = true;
        }
        ImGui::SetItemTooltip("Of the values in every plot cut here, over the frames they are shown over");

        double cut = periodic_wrap(view->min, p.beg, len);
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        if (ImGui::DragScalar("Cut", ImGuiDataType_Double, &cut, (float)(len / 360.0), nullptr, nullptr, "%.1f")) {
            periodic_view_cut(view, cut);
        }

        ImGui::Separator();
        if (ImGui::MenuItem("Show a whole period", nullptr, false, view->max - view->min < len * (1.0 - 1.0e-6))) {
            view->pending_min = view->min;
            view->has_pending = true;
            view->max = view->min + len;
        }
        ImGui::EndPopup();
    }
    return least;
}

// ## Workspace

void periodic_views_serialize(viamd::serialization_state_t& state) {
    for (int i = 0; i < MAX_PERIODIC_VIEWS; ++i) {
        const PeriodicView& v = views[i];
        if (!v.used) continue;
        viamd::write_section_header(state, STR_LIT("PeriodicView"));
        viamd::write_dbl(state, STR_LIT("Beg"),  v.period.beg);
        viamd::write_dbl(state, STR_LIT("End"),  v.period.end);
        viamd::write_dbl(state, STR_LIT("Seam"), v.period.seam);
        viamd::write_str(state, STR_LIT("Unit"), str_from_cstr(v.unit));
        viamd::write_dbl(state, STR_LIT("Min"),  v.min);
        viamd::write_dbl(state, STR_LIT("Max"),  v.max);
    }
}

bool periodic_views_deserialize(viamd::deserialization_state_t& state) {
    if (!str_eq_cstr(viamd::section_header(state), "PeriodicView")) return false;
    SeriesPeriodic p = {};
    char unit[32] = "";
    double min = 0, max = 0;
    bool has_min = false, has_max = false;
    str_t ident, arg;
    while (viamd::next_entry(ident, arg, state)) {
        if      (str_eq_cstr(ident, "Beg"))  viamd::extract_dbl(p.beg, arg);
        else if (str_eq_cstr(ident, "End"))  viamd::extract_dbl(p.end, arg);
        else if (str_eq_cstr(ident, "Seam")) viamd::extract_dbl(p.seam, arg);
        else if (str_eq_cstr(ident, "Unit")) viamd::extract_to_char_buf(unit, sizeof(unit), arg);
        else if (str_eq_cstr(ident, "Min"))  has_min = viamd::extract_dbl(min, arg);
        else if (str_eq_cstr(ident, "Max"))  has_max = viamd::extract_dbl(max, arg);
    }
    if (series_periodic_valid(p)) {
        PeriodicView* v = periodic_view(p, unit);
        if (v && has_min && has_max && max > min) {
            v->min = min;
            v->max = MIN(max, min + series_periodic_length(p));
        }
    }
    return true;
}
