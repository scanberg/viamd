#pragma once

// A light-weight TRANSFER FUNCTION for volume rendering: a colour map gives the colour, over a range of values, and
// a curve of control points gives the opacity over the same range. What a value maps to:
//
//     t     = (value - range_beg) / (range_end - range_beg), clamped to [0, 1]
//     color = the colour map at t
//     alpha = opacity * curve(t),   curve piecewise linear through the points, constant beyond the first and last
//
// The curve is in [0, 1] and the opacity scales it, so the curve says WHERE the volume is seen and the opacity HOW
// MUCH: the opacity of a DVR is small, and a curve drawn in [0, 0.05] could not be edited. The points are in t, so
// they move with the range.
//
// This part is plain data and arithmetic (the tests link it); transfer_function_widget.h is its editor, an ImGui and
// ImPlot widget over a histogram of the values, and the baking with the colour map's colours.

#include <stdint.h>
#include <stddef.h>

#define TRANSFER_FUNCTION_MAX_POINTS 16

struct TransferFunctionPoint {
    float t = 0.0f;         // where over the range, 0 to 1
    float alpha = 0.0f;     // of the curve, 0 to 1
};

struct TransferFunction {
    int   colormap = 0;             // ImPlotColormap
    float range_beg = 0.0f;         // the values the colour map runs over, t 0 to 1 (a range may run backwards)
    float range_end = 1.0f;
    bool  symmetric = false;        // the range kept as [-r, r]: for a signed field, about its 0
    float opacity = 0.1f;           // the alpha at curve value 1: of a slab 1/150 of the volume thick, in the DVR
    int   num_points = 2;           // at least 1, in order of t
    TransferFunctionPoint points[TRANSFER_FUNCTION_MAX_POINTS] = {{0.0f, 0.0f}, {1.0f, 1.0f}};
};

// t of a value
float transfer_function_param(const TransferFunction& tf, float value);
// The value at t
float transfer_function_value(const TransferFunction& tf, float t);

// The curve at t, 0 to 1 (not scaled by the opacity)
float transfer_function_curve(const TransferFunction& tf, float t);

// The alpha, opacity * curve, at n evenly spaced t from 0 to 1, into out[i * stride] (stride in floats, 1 = packed)
void transfer_function_bake_alpha(const TransferFunction& tf, float* out, int n, int stride);

// A point at (t, alpha), in order, both clamped to [0, 1]: its index, or -1 when there is no room
int transfer_function_insert(TransferFunction* tf, float t, float alpha);
// Removes point i, unless it is the last one left: whether it did
bool transfer_function_remove(TransferFunction* tf, int i);
// Moves point i to (t, alpha), kept between its neighbours in t and in [0, 1]
void transfer_function_move(TransferFunction* tf, int i, float t, float alpha);

// The range set, as [-r, r] when it is symmetric (r the larger magnitude of the two)
void transfer_function_set_range(TransferFunction* tf, float beg, float end);

// Of everything that decides what it maps values to
uint64_t transfer_function_hash(const TransferFunction& tf);

// Curves to start from
void transfer_function_ramp(TransferFunction* tf, float from);          // 0 up to 'from', rising to 1 at t = 1
void transfer_function_ramp_magnitude(TransferFunction* tf, float from); // 0 within 'from' of the middle, 1 at both ends (a signed field)
