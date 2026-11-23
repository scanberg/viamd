#pragma once

// THE MOUSE IN VIAMD'S VIEWS
//
// One set of gestures for every view, 2D and 3D: a gesture means the same thing wherever it means anything, and a view
// may leave one unused. The 3D viewport is where they come from.
//
//   wheel                   zoom
//   right-drag              pan
//   left-drag               3D: turn (the viewport orbits). 2D: pan as well, as nothing turns there.
//   middle-drag             3D only: dolly, zoom (the viewport). Unused in 2D.
//   shift + left-drag       select: a region, or what is under a click. Shift + right-drag deselects; a shift + right
//                           click on nothing clears.
//   click                   pick, or go to (a frame)
//   right-click             without a drag: the context menu
//   double-click            fit, or back to all of it
//
// There is no box zoom: the wheel and double-click do that. A view with something of its own to drag with the left
// button (the timeline: ctrl-drag the current frame, shift-drag a period) keeps the left button for it while it does,
// and pans with the right only.
//
// The plots' maps are set to these once, when their contexts are made (application.cpp). A view that needs something
// else changes the map before its BeginPlot and restores it after, and says why.

#include <implot.h>
#include <implot3d.h>

// ImPlot (2D): left or right drag pans, double-click fits, the wheel zooms, a right click opens the menu. Shift is for
// selecting, which the views do themselves: ImPlot ignores the mouse while it is held.
static inline void plot_input_house(ImPlotInputMap* map) {
    ImPlot::MapInputDefault(map);
    map->Pan         = ImGuiMouseButton_Left;
    map->PanAlt      = ImGuiMouseButton_Right;
    map->Fit         = ImGuiMouseButton_Left;
    map->Select      = -1;
    map->Menu        = ImGuiMouseButton_Right;
    map->OverrideMod = ImGuiMod_Shift;
}

// ImPlot3D: left-drag turns, right-drag pans, middle-drag and the wheel zoom (its own), double-click fits, a right click
// opens the menu
static inline void plot_input_house(ImPlot3DInputMap* map) {
    *map = ImPlot3DInputMap();
    map->Rotate        = ImGuiMouseButton_Left;
    map->Pan           = ImGuiMouseButton_Right;
    map->Fit           = ImGuiMouseButton_Left;
    map->ResetRotation = -1;
    map->Menu          = ImGuiMouseButton_Right;
}
