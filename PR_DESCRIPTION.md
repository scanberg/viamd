## Movie maker: keyframed camera, timeline editor and video recording

### Summary

A **movie maker** for VIAMD. Set the length of the movie, place camera keyframes and keyed look settings on a timeline, add overlays, and record to MP4 (H.264 / H.265), WebM (VP9) through ffmpeg, or to a numbered PNG sequence.

Everything is in one **Windows > Movie** editor: a timeline with draggable lanes on the left, tabbed controls (**Output**, **Timing**, **Camera**, **Looks**, **Representations**, **Overlays**) on the right, the camera keyframe table across the bottom, and shared preview playback on top. A **Scene view / Movie preview** switch gives you a free editor camera next to the movie camera, with a live picture of the movie camera while you edit. A green **Preview** button plays the movie full-viewport with every window hidden and a small control bar at the bottom.

The manual is [`docs/movie_maker.md`](docs/movie_maker.md). [`docs/examples/aspirin_phospholipase_movie.via`](docs/examples/aspirin_phospholipase_movie.via) is a finished 80 s movie using most features, and [`docs/examples/aspirin_binding_movie.via`](docs/examples/aspirin_binding_movie.via) a 60 s one with a dolly zoom, a bullet-time orbit, representation hand-overs and energy plots (their trajectory files are not in the repository).

![Current Movie editor with the aspirin binding example](docs/images/movie/overview.png)

*Native Linux capture, 12 October 2026: Scene view, camera-work lanes, live picture and the full-width camera table. The manual includes the current Output, Timing, target, look, representation and overlay inspectors, plus fullscreen Preview.*

### Features

**Recording**
- MP4 H.264, MP4 H.265 or WebM VP9 through an ffmpeg pipe (CRF setting), or PNG sequence with a button that copies the ffmpeg command.
- Window size, presets up to 8K, or custom resolution; output FPS.
- Frames are written off the render thread by `frame_sink` with a fixed buffer pool, so memory stays bounded.
- Quick tests: **Scale** (100/75/50/25 %), **Samples per frame**, **Render only a range** (PNG files keep their whole-movie numbers).
- **Pause** / **Resume**, time left, optional workspace copy (`prefix.via`) next to the movie.
- The user's view, playback and screenshot settings are restored when it ends.

**Timing**
- The movie length is the master value: changing it scales everything on the timeline.
- Movable trajectory start/end anchors (fly-over before the dynamics, hold at the end); frame pins on keys vary the speed and allow backward play; **Match Animation speed**; **Snap to frames**.

**Camera**
- Keys with per-key ease (Smooth, Ease in/out, Linear, Hold), field of view, names, spin turns, **Add Orbit**, seamless loop, **Key on Selection**, copy/paste, in a full-width, resizable keyframe table.
- Per key, **Look at**: *Fixed point* or *Track selection* (atoms followed through the trajectory), and **Focus on**: *Look-at* or *Track selection*, with transition and ease. Selections can be typed in the filter language. The blur is the keyable **Depth of field blur** look parameter.
- Three viewport paths (camera, look-at, focus) with draggable handles.
- **Keep upright** with a roll per key, so loops and spins never tilt the camera.
- The movie owns focus while its editor is open; the ordinary viewport always focuses at its look-at point (Auto), and its blur is restored when the editor closes.
- The camera path is drawn in the viewport (ticks at round times, chevrons, sight lines, cameras, spin rings, the green camera at the preview time) with handles on the eye and the look-at point of each key: click to go there, drag to edit (Ctrl moves both), Ctrl + click on the path adds a key.

**Looks and representations**
- Keyable look parameters: background, ambient occlusion, exposure, depth-of-field blur, clipping.
- Keyable representation properties: visibility (grow in / shrink away over a transition, since solid representations cannot fade), scales, tint, saturation, colors. Representations get a stable `id` that keys refer to, so reordering, duplicating and removing are safe.
- The **Representation overview** lane shows when each representation is shown, one row per system (name before the first hyphen), with block editing and swap. The Representations window is locked while recording.

**Overlays**
- Text, time stamp, scale bar, time bar, logo, images, **timeline** and **distribution** (subplots of the Timelines and Distributions windows, drawn as the movie plays, with markers), and the visualization of a script property. Overlays shown at the same time are moved apart so their plates never overlap, and marker labels that would collide go on further rows.
- Time range with fades, nine anchors, size in percent or points, color, background plate, **Duplicate**. Compact list with a Content / Timing / Appearance inspector.
- Subplots of the plot windows get saved ids and names so overlays keep finding them.

**Timeline panel**
- Ruler, trajectory, camera, lens (field of view and distance on two axes), look parameter, representation, overview and overlay lanes, toggled from one row with presets (Camera work, Scene work, All).
- Pick keys, overlay bars and overview blocks (click, Ctrl + click, box) and move, stretch, nudge, delete, copy and paste them together; an inspector sets start time and ease for the whole selection.
- Time zoom, middle-button pan, lane height / fit to window, lane selection saved in the workspace.

**Undo/redo** (Ctrl+Z, Ctrl+Y) cover everything on the movie.

**Shortcuts:** K / Shift+K add a key, Tab switches Scene view / Movie preview, Space plays, Ctrl+Z / Ctrl+Y undo and redo, from the Movie window or the viewport.

**Workspace**
- Saved under `[Movie]` (with `[MovieOverlay]` and `[MovieMarker]`), versioned with `Timeline=2`; older formats are migrated on load.
- Look parameter ids, representation property numbers, overlay kinds and lane/path option bits are stable and must never be renumbered.

### Fixes outside the movie code
- **Locale (`application.cpp`):** GTK (NativeFileDialog) switched the numeric locale, so on decimal-comma machines workspaces were written with commas. `LC_NUMERIC` is reset to `"C"`.
- **Folder picker:** `FileDialogFlag_Dir` is supported.
- **Space** plays the movie instead of the trajectory while the Movie window is open (`density_volume.cpp` honours this too).
- **`extract_flt_vec`:** accepts up to 32 values and rejects input with more values than requested, so newer keyframe formats load and older ones fall back correctly.

### New code
- `src/frame_sink.{h,cpp}`: asynchronous PNG/ffmpeg writer.
- `src/movie_keys.{h,cpp}`: parameter and representation keys, key selection and group edits, undo snapshot, time scaling, render range, snapping.
- `src/movie_overlay.{h,cpp}`: overlay maths, frame fit, migration of older overlays, default logo overlay.
- `src/plot_series.{h,cpp}`: stable subplot ids and names.
- `src/image.{h,cpp}`: in-memory image decoding for the logo (`icon/viamd.png`, baked in by CMake) and image overlays.
- `src/gfx/camera_utils.{h,cpp}`, `camera.h`: keyframe evaluation, anchors, keyed curves, `camera_aim_at`, `camera_level`, `camera_roll`.
- `Representation::id` (`viamd.h`).
- Most of the UI is in `src/main.cpp`; state and workspace I/O in `viamd.{h,cpp}`.

### Tests
New suites: `test_camera_utils`, `test_movie_keys`, `test_frame_sink`, `test_movie_overlay`, `test_image`, `test_serialization`. All 268 tests pass in a Release build on Linux (`./build/bin/viamd_test`).

### Testing done and not done
- PNG sequences and H.264 MP4s were recorded on Linux.
- The editor was checked by hand on Linux. [`MOVIE_STATUS.md`](MOVIE_STATUS.md) has developer notes and what is left out.
- Both aspirin examples were loaded with the actual trajectory and energy files. Their blur schedules and the binding generator use active aperture keys, re-authored for the current renderer. The binding script evaluated successfully; its bullet-time Preview was rendered and captured with 2 % aperture.
- Merged with `scanberg/viamd` master through `ac191f75`.
- Not tested on Windows or macOS (the ffmpeg pipe uses `_popen` / `popen`).

### Known limits
- Solid representations cannot fade; they grow in and shrink away.
- Tint and saturation keys recolor atoms every frame: slow for very large systems.
- Workspaces with separate look-at / focus tracks or follow keys are converted to per-key look-at and focus on load; the eye path of such a movie can change slightly.
- Timeline / distribution overlays draw at most six series per subplot and the first member of a population.
- Picked overlay bars and representation blocks move with picked keys but are not copied or pasted.
