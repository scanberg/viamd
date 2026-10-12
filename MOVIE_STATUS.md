# Movie feature: developer notes

Branch `video`; no pull request yet. Merged with `scanberg/viamd` master through `ac191f75`. The GUI has been checked by hand on Linux.

The user manual is [`docs/movie_maker.md`](docs/movie_maker.md); this file is for developers: where the code is, what must stay compatible, and what is left.

| File | What it is |
|---|---|
| `docs/movie_maker.md` | The user manual (also the text for the GitHub wiki, a separate repository) |
| `PR_DESCRIPTION.md` | The pull request text; remove it before merging |
| `docs/examples/fov_and_focus.via` + `.xyz` | Self-contained FOV and focus tutorial; in the older track format, so it also tests conversion on load |
| `docs/examples/aspirin_phospholipase_movie.via`, `aspirin_binding_movie.via` | Example movies (their trajectory files are not in the repository) |
| `docs/aspirin_binding_movie_article.md`, `docs/examples/tools/aspirin_movie/` | Historical case study with a current-renderer update, plus the generation and analysis tools |
| `docs/images/movie/` | Current native Linux captures using the 60 s aspirin binding example; embedded in the manual |

## Build and test

- `cmake --build build --target viamd viamd_test -j8` (run `cmake -S . -B build` after adding files), then `./build/bin/viamd_test` (268 tests; `--filter=viamd_movie_keys.*` for the 109 movie-key ones).
- The tests cover logic only: key evaluation (eye, aim and focus separately, transitions, spin, easing), conversion of older workspaces, group edits, undo, time scaling, overlay maths and the frame sink (with a fake script, so no ffmpeg is needed).

## Documentation and example verification

- Native captures taken on Linux on 12 October 2026 cover the combined editor, bottom camera table, camera-work timeline, Output, Timing, Look at / Focus on, Looks, representations, overlay Content / Appearance, distribution bins and fullscreen Preview.
- All GUI figures use `aspirin_binding_movie.via`; viewport/layout choices and the tracking popups demonstrate edits. Scientific plots in the case study are retained as camera/data diagrams, not GUI captures.
- Both aspirin workspaces load in the current application with active aperture keys (parameter 9), without obsolete blur-key warnings or unknown blur rows. The 60 s script was evaluated on the actual `.gro`, `.xtc` and `.edr`, and the bullet-time preview was rendered and inspected with 2 % aperture.
- The examples and generator no longer emit legacy parameter-5 blur keys. Apertures are re-authored in percent of view height: sharp introduction and ending, moderate dolly-zoom blur, stronger bullet-time blur. The camera/timing design is unchanged.
- The generator's blur-emission block was executed separately and matched against the shipped binding keys; both schedules, relative data paths and preservation of all non-blur workspace settings were checked. Python syntax, documentation links/anchors and all 12 movie-image references passed validation. This is not a full rerun of the scientific analysis scripts.

## Where things live

| What | Where |
|---|---|
| All movie GUI: editor window, lanes (`draw_movie_*_lane`, `draw_movie_strip`, `draw_movie_ruler`), camera path overlay (`movie_draw_camera_path`), picking and group edits (`movie_key_point`, `movie_pick_click`, `movie_group_drag`, `movie_lane_box`, `movie_selection_edit`, `movie_selection_shortcuts`), Scene view / Movie preview (`movie_scene_view`, `movie_camera_pose`, `movie_scene_fit_path`, `movie_set_scene_view`), live picture (`movie_render_pip`, `render_scene`), recording, overlay drawing | `src/main.cpp` |
| Movie state (`ApplicationState::movie`) | `src/viamd.h` |
| Workspace save and load (`[Movie]`, `[MovieOverlay]`, `[MovieMarker]`, `[Representation]`) | `src/viamd.cpp` |
| Subplot ids and names of the Timelines / Distributions windows | `src/plot_series.cpp` |
| Async frame writer | `src/frame_sink.{h,cpp}` |
| Keyed look parameters and representation keys, key selection and group edits (shift, scale, resolve, delete, copy, paste), undo snapshot, render range, snapping, time left, representation system rows | `src/movie_keys.{h,cpp}` |
| Overlay maths (fades, sizes, scale bar, time bar, elapsed curves, ticks, histograms, subplot times, frame fit, migration of older overlays) | `src/movie_overlay.{h,cpp}` |
| Logo (`icon/viamd.png`, baked into `gen/viamd_icon.inl` by CMake) and image decoding | `src/image.{h,cpp}` |
| Camera path evaluation, levelling, roll | `src/gfx/camera_utils.{h,cpp}` |

Numbers that are saved and must never be renumbered, only added to: parameter ids in `movie_param_table` (main.cpp), `RepProp` (movie_keys.h), `MovieOverlayType` (movie_overlay.h; 6 and 7 are timeline and distribution, 9 is the old one-figure version split on read), the `Lanes` bits (1, 2, 4 tracks; 8 parameter; 16 representation lane; 32 overlay; 64 overview; 256 camera; 1024 ruler; 128, 512, 2048 are the "set" flags of the last three) and the `PathOptions` bits (mask 31, default 15; bit 16 shows the path in Movie preview).

## Behaviour that must stay

- **Focus ownership.** The movie owns focus while its window is open, in Preview and while recording (`movie_focus_owned`). On entry it saves the viewport's focus depth (`editor_focus_distance`) and blur (`editor_saved_aperture`) and starts unkeyed blur at 0; on exit it restores keyed parameters, then the blur. The ordinary viewport is always Auto (focus at look-at); old global modes are converted to Auto on load and Auto is what is saved.
- **Per-key targets.** `CameraKeyframe::look_set` (0 = fixed point) and `focus_target` (`LookAt`, `Point`, `Distance`, `Selection`; Point and Distance are legacy, read-only in the GUI). A key on `LookAt` returns focus to its look-at; it does not inherit an earlier target.
- **Shortcuts outside the `WantCaptureKeyboard` guard** (so they work over the Movie window): Shift+K, Tab, Ctrl+Z / Ctrl+Y / Ctrl+Shift+Z. All are blocked while typing, recording, and (except Shift+K) in Preview.

## Known limits

- Recording has not been tried on Windows and macOS (the ffmpeg pipe uses `_popen` / `popen`).
- Only some representation properties can be keyed (not type, filter, color mapping, the other colors or electronic structure settings). Tint and saturation recolor atoms every frame: slow for very large systems.
- The global follow target of older workspaces is not part of undo and can no longer be set (a key's **Look at** / **Focus on** replace it).
- Where a frame curve turns around the trajectory slows to a stop; for a hard reversal add a key at the turn-around frame with Linear ease.
- Overlay text has one style (shadow and optional plate).
- Timeline / distribution overlays draw at most six series per subplot and the first member of a population; a timeline overlay has one axis kind for all subplots.
- Tracked path drawings use the current trajectory frame, not future-frame sampling.
- A key that looks at atoms aims at their centre; the offset that a legacy follow key kept is not preserved on conversion.
- Picked overlay bars and representation blocks move with picked keys but are not copied, pasted or (overlays) deleted by the picking shortcuts; unkeyed blocks cannot be picked.
- With **Keep upright**, a spin around a horizontal axis flips at the poles; roll eases as a number, so 170 to -170 degrees turns the long way.
- Master code, not ours: a perfectly symmetric flat molecule (exact ideal benzene) can settle about 10 degrees off face-on in the default view.

## Left for later

- Text overlays with live values (`{time}`, `{frame}`, a script property) and entrance/exit animations.
- Saved overlay templates.
- Editing representations during a recording, and a transparent background.
- Not planned: merging the Timelines and Distributions windows, 2D plots, exporting figures other than the movie.
- Upstream strays left alone: `Claude outputs/windows-arm64.yml` (tracked on master) and `TODO.md`.
