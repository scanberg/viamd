#pragma once

#include <math.h>
#include <stdint.h>

#include <core/md_unit.h>

// COLOUR SCALES
//
// How scalar values become colours, wherever viamd colours by a value: the atoms of a representation
// by one of their attributes, an isosurface by a field. A scale is the colour map, the range it is
// spread over, how that range is chosen, and whether its legend is drawn over the view - the same
// controls and the same legend wherever it is used.
//
// A scale holds no values. What the values SPAN is measured by whoever has them (ColorScaleSpan),
// because what "the values" are differs from one use to the next: the atoms a representation shows,
// the samples of a field on the surface it colours.
//
// UNITS: the range is in the values' own unit, the unit the caller passes along with them. The
// controls and the legend show it in the display unit (display_units.h); nothing stored is converted.

// ImPlotColormap_* indices, so that a header need not include implot.h for a default
#define COLOR_SCALE_COLORMAP_VIRIDIS          4
#define COLOR_SCALE_COLORMAP_PLASMA           5
#define COLOR_SCALE_COLORMAP_RDBU             11    // red below zero, blue above, white between
#define COLOR_SCALE_FIRST_CONTINUOUS_COLORMAP 4     // 0-3 are qualitative (Deep, Dark, Pastel, Paired)

struct ColorScale {
    int   colormap    = COLOR_SCALE_COLORMAP_PLASMA;
    float range_beg   = 0.0f;
    float range_end   = 1.0f;
    bool  symmetric   = false;  // a range centred on zero
    bool  auto_range  = true;   // the range follows the span of the values
    bool  show_legend = false;  // the colour map with its range and unit, drawn over the view
};

// What the values span. lo / hi are what an automatic range takes: the extremes, or a robust part of
// the spread where the extremes are not representative of it. min / max are the extremes.
struct ColorScaleSpan {
    bool  valid = false;
    float min = 0.0f;
    float max = 0.0f;
    float lo  = 0.0f;
    float hi  = 0.0f;
};

static inline ColorScaleSpan color_scale_span(float min, float max) {
    ColorScaleSpan span;
    span.valid = min <= max;
    span.min = span.lo = min;
    span.max = span.hi = max;
    return span;
}

// Whether a span straddles zero: the values of a charge or a potential, which read best on a
// symmetric range and a diverging colour map
static inline bool color_scale_span_signed(const ColorScaleSpan& span) {
    return span.valid && span.min < 0.0f && span.max > 0.0f;
}

// The range an automatic scale takes from the span. Leaves a manual scale, or an empty span, alone.
static inline void color_scale_update_range(ColorScale* scale, const ColorScaleSpan& span) {
    if (!scale || !scale->auto_range || !span.valid) return;
    if (scale->symmetric) {
        const float r = fmaxf(fabsf(span.lo), fabsf(span.hi));
        scale->range_beg = -r;
        scale->range_end =  r;
    } else {
        scale->range_beg = span.lo;
        scale->range_end = span.hi;
    }
}

// Where a value falls on the scale: 0 at range_beg, 1 at range_end, clamped to [0, 1]. A range may
// run backwards, which reverses the colour map; a range of no extent puts every value at its middle.
static inline float color_scale_param(const ColorScale& scale, float value) {
    const float ext = scale.range_end - scale.range_beg;
    if (!(fabsf(ext) > 0.0f)) return 0.5f;
    const float t = (value - scale.range_beg) / ext;
    return t < 0.0f ? 0.0f : (t > 1.0f ? 1.0f : t);
}

// The colour of a value, as packed RGBA (ImU32)
uint32_t color_scale_color_u32(const ColorScale& scale, float value);

// The controls: colour map and legend, auto range and symmetric, the range, and what the values span.
// 'unit' is the values' own. 'span_label' says what the span is of ("on surface"), 'span_tooltip'
// how it was measured. Returns true when the scale changed.
bool color_scale_draw_controls(ColorScale* scale, const ColorScaleSpan& span, md_unit_t unit, const char* span_label, const char* span_tooltip);

// The legend: the colour map over its range with values, headed by 'label' and the unit, in a window
// without decoration over the view. Like the coordinate widget it takes no input unless Alt is held,
// and can then be moved and resized; where it was is kept in the .ini under 'id'. 'owner' names it in
// its title bar, which shows while Alt is held. 'slot' is its place among the legends drawn, which
// decides where it starts out.
void color_scale_draw_legend(const ColorScale& scale, md_unit_t unit, const char* label, const char* owner, int id, int slot);
