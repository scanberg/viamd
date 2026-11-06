#pragma once

#include <stdint.h>

// ISOLINES
// Lines on a surface coloured by a field, where the field takes given values: the contours of chosen values of
// the colour map, drawn by whichever renderer shades the surface (isomesh_utils.h, volumerender_utils.h) as part
// of its colour, so that they are lit, hidden and seen through like the rest of it.
//
// A line is found per pixel from the field's value there and its change from one pixel to the next along the
// surface (screen space derivatives on a mesh, ray differentials in the raycaster): its distance in pixels is
// |f - v| / |grad_screen f|, so a line keeps its width in pixels whatever the slope of the field or the surface.
// The lines are where the surface's colour is the colour map's at v, as they read the same values.
//
// Two kinds: the persistent ones, dark with a light rim, and one transient line (the value under the cursor of a
// legend), light with a dark rim, so either reads on any colour of a map and the two are told apart.

#define ISOLINES_MAX 8

struct IsolineDesc {
    uint32_t count = 0;                     // persistent lines
    float    values[ISOLINES_MAX] = {};     // in the field's own unit
    bool     hover = false;                 // the transient line
    float    hover_value = 0.0f;
    float    width = 1.5f;                  // of a line, in pixels of the render target
};

static inline bool isolines_active(const IsolineDesc& desc) {
    return desc.count > 0 || desc.hover;
}
