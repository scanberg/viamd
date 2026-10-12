# Unified series model for the plotting windows and Density Volume

Branch: `uinification` (fork `mathieulinares/viamd`), target: `scanberg/viamd` `master`.
Status: draft until the mdlib changes are merged upstream (see [Merge requirements](#merge-requirements)).

## Summary

Timeline, Distributions, Ramachandran, Shape space and Density Volume now share one way of working with series:

- a **Properties** menu first,
- a list of series with a per-series style,
- drag and drop to add a series and to remove it (drag it out of the window),
- filtered series shown as separate series,
- everything saved in the workspace.

Ramachandran and Shape space can be filled by writing a script property, so showing a selection is `r = ramachandran(sel);` and nothing more.

## What changed

### Scripting (ext/mdlib)

- New procedures `ramachandran(sel)`, `shapespace(sel)`, `ramachandran_weights()` and `shapespace_weights()`, each also callable without arguments for the whole system.
- Vector properties (`float[N][2]`, `float[N][3]`) are stored for every frame and published with their number of components (`MD_SCRIPT_PROPERTY_FLAG_VECTOR`).
- Fixes found on the way: destructuring stride, transpose of a single row, visualisation of properties without arguments.
- Reference and tests: `docs/script_reference.md`, `unittest/test_script.c`.

### Ramachandran and Shape space

- Series chips under the menu bar (Properties, Layers and View or Export menus), one chip per series.
- Per-series style: points for the displayed frame or every frame, isolines, filled levels between them, density colormap, number of levels, opacity, smoothing; colour per series or per member of a population.
- The same property can be shown twice with different styles (Duplicate). Filtered evaluations appear as series of their own.
- A shortcut row creates a property from a selection and adds it to the script.
- The densities are computed in the background and debounced, and no longer restart every frame.
- Workspace sections `Ramachandran`/`ShapeSpace`, `*Layout` and `*Series`.

### Density Volume

- **Properties** first menu, layers as chips with a style popup (isosurfaces or DVR, clip planes, bounding box, reference structures), a **Layers** menu with the same settings, a **View** menu.
- Duplicates (the same volume with another style), removal by dragging a chip out, click to show or hide.
- New layers default to **one isosurface** in a colour of their own. DVR stays available, but several DVR layers are hard to read.
- The up/down reordering is gone: layers are depth sorted along the ray, so the list order never mattered.
- Old workspaces load unchanged (legacy `[DensityVolume]` and `[DensityVolumeLayer]` sections). The coordinate system widget is one window setting.

### Bug fix: lost lighting with two or more densities

The multi-layer shader built its normal matrix from the upper 3x3 of a 4x4 product, `view * transpose(inverse(model))`. The translations of the model and the view leak into that block, so the normals were wrong for any layer which is not at the origin and the isosurfaces looked flat. One density takes another path and was not affected. The matrix is now built from the two 3x3 blocks directly. `test/gl/test_volume_overlay_gl.cpp` has a headless sphere check that the brightness follows the normal (it uses an identity model, so it does not cover the translation itself).

### Consistency between windows

- One series editor (`series_vector_editor` in `src/plot_series.cpp`) is used by the chips and by the plot legends: right click a legend entry in Ramachandran or Shape space as in Timeline.
- A click on a chip toggles the visibility of the series and stays in step with its legend entry.
- Space plays and pauses in every plot window.
- The Properties menu shows the colour a series has in the window once it is there.
- Timeline and Distributions can show the same series twice: each entry of a plot has its own ImPlot id (it carries the uid of the entry), and the legend popup has Duplicate like the chips.
- Every plotting window has an **Export** entry which opens the Property Export window on its kind with the first series selected. The window gained the kinds Ramachandran and Shape Space (one column per component of each member, angles in degrees, XVG or CSV), which replace the export menu Shape space had on its own.
- Every window says the same thing where it has no series (`series_empty_hint`).

### Docs and example

- `docs/density-volume.md` describes the new layout.
- `docs/examples/showcase_gui.via` opens every window with series (it needs the aspirin-phospholipase files next to it).

## Merge requirements

`ext/mdlib` points to the fork `mathieulinares/mdlib` (branch `uinification`) because the script procedures, vector properties and fixes above are in mdlib and the author has no write access to the upstream repository. Before this can be merged into `scanberg/viamd`:

1. The mdlib commits are merged upstream, rebased on the current mdlib API (upstream master renamed parts of the script API and removed `md_bitfield` from the places this branch used it).
2. `.gitmodules` and the submodule pointer are set back to the upstream repository.

## Testing

- `cmake --build build --target viamd viamd_test` builds.
- `./build/bin/viamd_test`: all pass except `analysis_density.empty_subset_and_invalid_samples_are_empty` and `viamd_camera.default_view_of_a_planar_molecule_is_face_on`, which also fail without this branch.
- `./build/bin/md_unittest` covers the new script procedures.
- Headless GL: `test/gl/test_volume_overlay_gl.cpp` (EGL with software rendering, outside the unit-test glob).
- The GUI was tested by hand by the author with the example workspace.
