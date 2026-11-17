#pragma once

// The editor of a transfer function (transfer_function.h), an ImGui / ImPlot widget of its own, to be put wherever
// a volume's colours and opacity are set:
//
//   - a plot over the values, panned and zoomed as any ImPlot plot (horizontally; the vertical is the curve's 0 to 1)
//   - behind it, the histogram of the values in the volume, over what is in view (the caller makes it: see view)
//   - the colour map over its range, as the area under the opacity curve, and as a strip along the bottom
//   - the range as two lines, dragged; the curve's points dragged, added by ctrl + click, removed by a right click
//   - under the plot: the colour map, the opacity, the range as numbers, symmetric about 0, fitting the view
//
// It reads no volume and owns nothing: the histogram comes in, and the part of the values in view goes out.

#include <transfer_function.h>

struct TransferFunctionHistogram {
    const float* counts = nullptr;  // bins of them, evenly over [lo, hi]
    int          bins = 0;
    float        lo = 0.0f;
    float        hi = 0.0f;
};

// Kept by the caller from frame to frame, one per editor shown
struct TransferFunctionView {
    double x_min = 0.0;     // the values in view: what the histogram is wanted over next
    double x_max = 1.0;
    int    fit = 1;         // fit the view at the next frame: 1 to the range, 2 to the data (set to 0 once done)
};

// data_min, data_max: the values in the volume, what 'fit' to the data fits to (not offered when they are equal).
// hist may be NULL. width <= 0: all there is. True when tf changed.
// colors_only: the colour map over its range alone, for what takes no opacity (a plane drawn opaque): no curve, no
// opacity.
bool transfer_function_editor(const char* id, TransferFunction* tf, TransferFunctionView* view,
                              const TransferFunctionHistogram* hist, float data_min, float data_max,
                              bool colors_only = false, float width = 0.0f, float height = 150.0f);

// n RGBA float texels from t = 0 to 1: the colour map's colours, and the alpha (transfer_function_bake_alpha)
void transfer_function_bake_rgba(const TransferFunction& tf, float* rgba, int n);
