# Density Volume overlays

Open **Windows → Density Volume**. Select properties in its **Property** menu to
add layers, or drag them from that menu or from the cube-marked rows in
**Attributes** onto the 3D canvas. Plotting drag payloads remain separate: timeline
and distribution rows still go into their respective subplots.

Up to six full/filtered scalar volume properties share one camera. Adding an
existing property selects and shows its layer rather than duplicating it. Select a
layer row to edit that property's **Render**, **Clip planes**, and **Show** settings.
The row checkbox hides it, **X** removes it, and **Up / Dn** rearrange the list.
Reordering is organizational: it does not change depth or transparency ordering.
Reference structures, their atom picking, and camera interaction remain available
for each visible, valid layer. The first added volume fits the camera; subsequent
layers do not reset it.

Each layer keeps its source/path identity, texture, transfer function, clipping,
isosurface values/colors, reference representation settings, and data version.
Missing or invalid properties remain in the list but render neither stale density
nor stale reference structures. Recompiling a script resolves the paths again;
changed grid dimensions reallocate that layer's texture.

## Rendering

A single visible rendering layer uses the existing DVR or isosurface renderer,
including the **exact isosurfaces** setting. Multiple rendering layers use one
shared ray traversal, stopped by the opaque reference-structure depth. DVR
extinction/emission is combined within each ray step; refined surface crossings
are sorted by distance and composited between the DVR segments. DVR colormap
colors are decoded to linear light and mixed with lit surfaces in the HDR
transparency target before tone mapping. Completed per-property images are **not**
blended in list order.

Overlay surfaces use the fast, sampled intersection method (half-voxel steps and
bisection refinement), not the exact single-volume cell solver. Sub-voxel features
and multiple crossings within one step can be missed. The window labels this mode
explicitly. For exact intersection of one property, hide the other rendering
layers. Coincident surfaces have no unique physical front/back ordering; separating
their isovalues avoids that ambiguity.

Workspaces write a `[DensityVolume]` header and repeated `[DensityVolumeLayer]`
sections, preserving layer order, selected layer, visibility, and per-property
settings. Legacy single-property `[DensityVolume]` sections load as one layer.

## Optional headless rendering regression

On Linux with EGL/OpenGL development headers and a surfaceless driver:

```sh
c++ -std=c++17 -Iext/mdlib/src test/gl/test_volume_overlay_gl.cpp \
    -o build/test_volume_overlay_gl $(pkg-config --cflags --libs egl gl)
LIBGL_ALWAYS_SOFTWARE=1 build/test_volume_overlay_gl
```

This standalone test compiles the actual overlay shader, then checks overlapping
DVR extinction, opaque-depth truncation, isosurface front/back ordering, and mixed
DVR/isosurface order independence. It is outside the CPU-only unit-test glob and
adds no required application dependency.
