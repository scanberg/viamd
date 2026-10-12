# Making movies

VIAMD records a movie of a trajectory. You choose how long it is, place camera keyframes, key look settings (exposure, background and so on) and the properties of representations (show or hide one, change its size or color), and add overlays (text, a time stamp, a scale bar, a time bar, images, plots of trajectory properties, the visualization of a script property, the VIAMD logo). The result is an MP4 (H.264 or H.265) or WebM file encoded by ffmpeg, or a numbered PNG sequence.

![The movie tools in an earlier layout, with the camera path in the viewport](images/movie/overview.png)

*The screenshots show an earlier layout with the timeline and the controls in two windows. They are now the two panels of one window, and the settings are the same.*

Contents: [The editor](#the-editor) · [Quick start](#quick-start) · [Output](#output) · [Timing](#timing) · [Camera](#camera) · [Looks](#looks) · [Representations](#representations) · [Overlays](#overlays) · [The timeline](#the-timeline) · [Saving](#saving) · [Example](#example-workspace) · [Tips](#tips) · [Known limits](#known-limits)

## The editor

Open it with **Windows > Movie**. One window holds everything:

**Ctrl+Z** undoes movie edits with either the Movie window or viewport active;
**Ctrl+Y** or **Ctrl+Shift+Z** redoes them. While typing, these shortcuts belong
to the text field. Movie undo/redo shortcuts are disabled during recording and fullscreen Preview.

- **Left, the timeline:** lanes on a shared time axis for the trajectory frame, camera, lens (field of view and distance), looks, representations and overlays.
- **Right, the controls:** tabs **Output**, **Timing**, **Camera**, **Looks**, **Representations** and **Overlays**.
- **Bottom, across the full width:** the camera keyframe table, available from every tab. Drag its horizontal divider to resize it; the height is saved with the workspace. The upper area gives most of its width to the timeline.
- Drag the divider to resize the panels. Each scrolls on its own. The **Timeline** and **Controls** boxes hide either one (one always stays).
- **Top:** playback (**Play Preview**, **Play from start**, **Repeat**) and the preview-time scrubber.
- **Right, above the tabs:** **Undo** and **Redo** (Ctrl+Z, Ctrl+Y or Ctrl+Shift+Z) cover everything on the movie: camera keys, look and representation keys, overlays, the length and the trajectory timing. Below them the **Keyframes** section: **Add Keyframe**, **with trajectory frame**, **Key on Selection**. While recording, **Stop Recording**, **Pause** and the progress are shown there too.
- **Start Recording** is at the bottom of the **Output** tab, with a summary of frames, seconds and pixel size.

### Scene view, Movie preview and Preview

The large button at the top right (or **Tab**; hold Tab to peek at the other one) switches the viewport between two views. Tab also works when the Movie window has keyboard focus, except while typing. Neither changes the keys, the playback or the recording, and the switch is disabled while recording.

| | **Scene view** | **Movie preview** |
|---|---|---|
| Camera | Your own editor camera | The movie camera |
| Shows | The recording frame (outside dimmed), the camera path and key handles | The recording frame (outside dimmed), the overlays, and the path only if **Show path** is ticked in this mode |
| Scrubbing and playing | The green movie camera moves along the path; the trajectory plays; your view stays | The viewport is the movie camera, following the keys: a camera turning around the system looks like the system turning |

Each mode remembers its own view. The first time you enter Scene view it frames the whole camera path, seen from where you were looking; **Fit path** does that again. Both frame it, close and centred, in the frame. Both modes show the same recording frame, the same size and in the same place (in the part of the viewport the Movie window leaves free, beside, above or below it), with what is outside dimmed: what is in the frame is what a key added from the view records, so keys are aimed the same way in both modes. The view turns around the middle of the frame, and in Scene view the mouse wheel zooms towards the point under the mouse. **Add Keyframe** (K) adds a key at the preview time with what the frame shows, in both modes; **Update position** in the key table works in both too. To add a key where the path already has the camera, without changing the path's shape, right click **Add Keyframe** > **Add a key on the path**, or double click the Camera lane at that time.

Next to the Timeline and Controls boxes: **Fit path** (Scene view only), **Show path** and **Live picture** (Scene view only). **Show path** has a setting of its own in each mode: in Scene view it draws the path, in Movie preview it adds the path over the movie camera's view.

In Movie preview the frame is placed in the part of the viewport the Movie window leaves free, so the window never hides the picture, with the movie camera's view centred in the frame as in the recording (the projection is shifted with the frame; the camera and the keys do not change). Move the window and the frame moves along. When the window floats in the middle, or is moved out of the main window, the frame is centred in the whole viewport.

**Live picture** shows a small live picture of the movie camera, with overlays, at the lower right of the upper panels, above the keyframe table, so you can drag the path in the viewport and watch the shot at the same time. Right-click it for its size (**Small**, **Medium**, **Large**, **Extra large**) or to hide it. The sizes grow with the Movie window's width, limited to about 65 % of the upper area's height. It costs one extra render when the picture changes (at most about 20 times a second), so switch it off on very large systems. It has no temporal anti-aliasing or sharpening, and no property visualization.

### Preview: watch the movie

The green **Preview** button next to the switch plays the movie in the viewport as it will be recorded, with every window, the menu and the camera path hidden and black outside the frame. It always starts from the beginning of the movie. A bar at the bottom has play/pause, back to the start, the time slider (scrubbing pauses, release resumes), **Repeat** and **Exit preview**. It hides after a moment while the mouse is still and the movie plays; move the mouse to bring it back. **Space** plays and pauses, **Esc** leaves. Whenever the Movie window is up, **Space** plays and pauses the movie (as **Play Preview**), not the trajectory; in a text field it types a space. Opening the Movie window stops the trajectory if it was playing; closing it stops the movie, and **Space** plays the trajectory again. Leaving stops the playback and goes back to Scene view if that was shown. It cannot be entered while recording.

## Quick start

1. Load a trajectory and open **Windows > Movie**.
2. **Output:** click **Select Output Folder...** and pick a format.
3. **Movie length (s)** (Output or Timing): 60 s for a new movie; type another and press Enter. At first the trajectory fills the movie.
4. Set the preview time to 0, set up the view, press **K** (**Add Keyframe**).
5. Move the preview time, set up another view, press **K** again. Repeat.
6. Press **Play Preview** (or **Play from start**). Switch to **Movie preview** to see the frame and overlays, or press the green **Preview** to watch it without the windows.
7. **Start Recording** (bottom of **Output**). Esc or **Stop Recording** stops early.

When the recording ends, your own view, the animation playback and the screenshot settings are restored.

## Output

![The Movie window](images/movie/movie_window.png)

| Setting | What it does |
|---|---|
| **Select Output Folder...** | Where the movie is written. Recording is disabled until one is chosen. |
| **Filename Prefix** | `prefix.mp4`, or `prefix_00000.png`, `prefix_00001.png`, ... |
| **Format** | **MP4 H.264**, **MP4 H.265** and **WebM VP9** pipe frames straight into ffmpeg (it needs libx264, libx265 or libvpx-vp9 respectively); no image files are written. H.265 files are usually noticeably smaller and are tagged so QuickTime plays them. **PNG sequence** writes one image per frame. |
| **Quality (CRF)** | Video only. Lower is better and larger. H.264: 18 is close to lossless, 23 the default. H.265: 4 to 6 higher looks the same (default 28). VP9: 15 to 35 is usual (default 32, up to 63). |
| **ffmpeg** | Video only. The executable; a bare name is looked up on the PATH. |
| **Copy ffmpeg command** | PNG only. Copies a command that encodes the sequence into an MP4. |
| **Resolution** | **Window** (its size when recording starts), Full HD, Quad HD, 4K, 8K, or **Custom** (**Res X**, **Res Y**). |
| **Output FPS** | Frames per second. The frame count is length × FPS. |
| **Scale** | 100, 75, 50 or 25 % of the size above. Small renders much faster: use it to check a movie first. |
| **Samples per frame** | With temporal anti-aliasing: how many rendered images are averaged per frame. 0 uses the jitter sequence length. Fewer is faster, more smoother. |
| **Render only a range** | Renders the frames between two times (**Range (s)**), e.g. to redo a part. A PNG sequence keeps the frame numbers it has in the whole movie, so the files can replace the old ones. **Start at preview time** and **End at preview time** take the times from the playhead. |
| **Save a workspace copy with the movie** | Writes `prefix.via` next to the movie with the camera path, looks, overlays and settings. Opening it and recording again gives the same movie. |

While recording, **Pause** and **Resume** stop and carry on. The progress bar shows frames done and, after a few frames, the time left. The **Representations** window is locked while recording.

## Timing

![Timing and camera keyframes](images/movie/timing_and_keyframes.png)

- **Trajectory Frames:** the first and last frame used. Start after end plays backwards.
- **Movie length (s):** applied on Enter. **Everything on the timeline scales with it** (camera keys, look and representation keys, overlays and fades, the trajectory start and end), so the movie keeps its shape and only gets faster or slower. Frame numbers and spin counts are unchanged.
- **Trajectory plays (s):** the movie times at which the trajectory is at its first and last frame. Before the start it is held on the first frame (for a fly-over first), after the end on the last. These are the blue **Start** and **End** lines in the timeline, where they can also be dragged.
- A grey line shows the average speed in frames per second and relative to the Animation window. **Match Animation speed** sets the movie length so the trajectory plays at that speed.

To vary the speed, pin frames on camera keyframes (**Frame** column below): between two pins the trajectory plays at whatever speed gets from one to the other. A pin before Start or after End moves that line out to it.

## Camera

A keyframe is a camera (eye position, look-at point, distance, field of view) at a time, and it can say what it **looks at** and what is **in focus**, each without touching the rest. Between keys the camera moves smoothly; before the first and after the last it holds still.

In Movie preview, previewing and recording always follow the keys; in Scene view your own camera stays and the green movie camera moves along the path.

### Keyframes section (above the tabs)

- **Add Keyframe** (**K**): key of what the frame shows at the preview time. When there is a key at the preview time already, the new one goes **2 s later** and the preview time moves there (refused past the movie end). Right click: **Add a key on the path**. **with trajectory frame** is checked by default and also pins the frame shown now.
- **Key on Selection:** key at the preview time framing the selected atoms, from the current direction, and moves the viewport there. Molecules split by the periodic cell are framed as one.

### Camera tab

- **Add Orbit:** fields for turns, duration and axis (**Camera up**, **World X/Y/Z**) set a key at the preview time and another one *duration* later after the camera has gone round the look-at point. Refused if it would run past the movie end.
- **Seamless loop** and **Close Loop:** make the path cyclic so a repeating movie does not jump. **Close Loop** adds a final key in the first key's pose and turns the loop on.
- **Keep upright** (on for new movies), **Up**, **From view:** the camera stays level with the chosen world axis up on screen, so loops and spins cannot leave it tilted. **From view** picks the axis closest to the view's up. Keys then give only where the camera looks from and at; the tilt is the key's **Roll (deg)** (positive leans left, 0 levels), eased smoothly. A key made from the view (**K**, **Add Keyframe**, **Key on Selection**, **Update position**) takes the view's own tilt as its Roll, so it plays back as it was seen; set Roll to 0 to level it. A spin around **Camera up** turns around the up axis. Off, the keys' own tilts are used, as in old workspaces; the **Roll** column then shows each tilt with a **Level** button.
- **Paste Keyframe:** puts in the key copied from the table (see below).

### The keyframe table

![Keyframe table, key 2 waiting for an atom to be clicked](images/movie/look_at_pick.png)

The table occupies the full-width bottom section, independently of the selected controls tab. Columns resize by dragging their borders; the table scrolls sideways when wider than the window.

| Column | Meaning |
|---|---|
| **#** | The key's number. Drag it onto another row to move the key's pose there; times stay. |
| **Name** | A name of up to 23 characters ("intro", "close-up"), shown in the camera lane and saved. |
| **Time (s)** | When. The list re-sorts when you finish editing. |
| **Frame** | Tick to pin a trajectory frame to the key, then drag the number. A frame lower than the previous pin plays backwards. |
| **FOV (deg)** | Field of view. Narrow it to zoom without moving. |
| **Go to** | Moves the viewport to the key. Clicking a key's orange marker in the **Timelines** window does the same. |
| **Update** | Moves the key's eye to the viewport camera, still looking at the same point; also takes the FOV and, if pinned, the frame. |
| **Look at** | **Fixed point** or **Track selection**. Edit the point/use view aim, or track a typed selection, selected atoms or a picked atom. The eye keeps its place. See [Look at and focus on](#look-at-and-focus-on). |
| **Focus on** | **look-at** follows the movie camera's aim at this key. Choose a separate target to decouple focus, with its own transition and ease. |
| **Ease** | How the movie gets *to* this key: **Smooth** passes without stopping, **Ease in/out** slows at both ends, **Linear** is constant speed, **Hold** stays then jumps. |
| **Spin** | Extra whole turns around the look-at point on the way to this key. **...** opens the axis and **Constant speed**. |
| **Roll (deg)** | See **Keep upright**. |
| **Copy** | Remembers the key. **Paste Keyframe** puts it in at the preview time, replacing one that is there. |
| **Dup** | Copies the key to one second later. |
| **Remove** | Deletes the key. |

**Shift+K** captures a key while the Movie window or viewport has keyboard focus,
including fullscreen Preview. It is disabled during recording and while typing in a
text field. Plain **K** keeps its viewport shortcut behavior.

To edit a key: go to it, move around, **Update position**, and **Look at** / **Focus on** to choose the targets. A line below the table shows the eye, look-at point and distance at the preview time (and the focus distance with depth of field).

### Look at and focus on

Every key can say two things besides where the eye is, each on its own and independently of the other keys' times:

- **Look at:** what the camera is aimed at. A **fixed point** (the default: where the pose aims) or **atoms**, tracked through the trajectory: a selection typed in the language of the filters (`protein`, `resname("AIN")`, `resid(12)`, a named selection ...; Enter applies it, red when it is not valid or selects nothing), the current selection, or one atom picked in the viewport. An expression is evaluated again for the system and frame that are shown, so it keeps working when the system is replaced; the centre is mass-weighted and made whole across periodic boundaries. The aim moves from one key's look-at to the next with the key's easing, so a camera can look at one molecule and then turn to another while the eye travels its own way. The look-at path is the *yellow* one in the viewport.
- **Focus on:** what is sharp. Choose **Look-at** (default) or **Track selection**, independently of aim. A tracked target can be a typed selection, selected atoms or a picked atom, with **Transition** (seconds, from the key's time) and **Ease**. Choosing Track selection uses the selected atoms, or starts viewport picking if nothing is selected; expressions can be applied directly in either popup. Choosing Fixed point for look-at freezes the current tracked centre. A default focus key returns to look-at, rather than inheriting an earlier target. Explicit transitions finish no later than the next camera key. Saved legacy fixed-point/distance focus keys remain readable and retain their effect until you choose a different target. The focus path is the *magenta* one.

**How blurred the rest is** is the look parameter **Depth of field blur** (Looks tab and look-parameter lane), not part of the key: key it, ease it and smooth it like any other look parameter, with **Key the blur here** in the Focus on popup as a shortcut. Movie blur is on when that parameter is above 0, without the global depth-of-field switch; zero makes everything sharp, regardless of the focus target.

Opening Movie sets its unkeyed blur to **zero**, independently of the ordinary viewport's blur. Saved movie blur keys still apply at the preview time and during playback; closing Movie restores the ordinary viewport's blur.

Where a key is selected, **Look at** and **Focus on** act on that key only. Their targets are saved with the key and are copied, pasted, moved, stretched and deleted with it; undo covers them, and a movie's length scales their transitions.

Movie preview, the live picture, fullscreen Preview and recording always evaluate movie focus, including default look-at keys. Global controls show **Focus: driven by movie** and are read-only while the movie is open. The Scene view's editor camera keeps a separate internal depth captured on entry; movie target changes do not pull its focus. Closing the movie restores the ordinary viewport's Auto focus and saved blur settings. Atoms that are not in the current system are reported once and the key's saved point is used.

Workspaces made with the separate **Look-at** and **Focus** tracks, or with keys that followed a follow target or an atom, are converted when loaded: a camera key is added (on the path, so the movie goes the same way) where a track had a key, tracked keys look at the centre of their atoms, the focus keys become **Focus on** of the key at their time and their blur becomes keys of **Depth of field blur** with the same transitions. The eye path of such a movie can differ slightly from before.

### The camera path in the viewport

**Show path** (at the top of the Movie window) draws the path over the viewport while the editor is open (not in recorded frames). What to draw (time ticks, sight lines, cameras, spin rings) is set here in the Camera tab.

- The eye's path is blue, the look-at path yellow; each is brighter ahead of the preview time than behind.
- **Time ticks:** dots at round times (the step adapts to the movie length) with the time written at some, and chevrons for the direction. Close dots mean slow, far apart mean fast, so the easing can be read off. Yellow diamonds mark the same times on the look-at path.
- **Sight lines:** a thin line from the eye to its target at each tick. **Cameras:** a camera drawn at each key, orange at the preview time. **Spin rings:** the circle a spin goes round, with arrows and the turn count.
- **Key handles:** a round one on the eye (with the key's number and name) and a ring with a cross on what it looks at, joined by a line. White for a key that looks at a fixed point; blue for one that looks at atoms; orange at the preview time. Hovering one highlights the key in the camera lane and vice versa. **Click** goes to the key and picks it. **Drag** the eye handle to move the eye (the target stays), the ring to move the target (the eye stays), **Ctrl** to move both. The handle moves in the plane facing you; Esc cancels; a drag is one undo step.
- **Ctrl + click on the path** adds a key at that time without moving the camera.
- The **green camera** is the camera at the preview time: a green dot on the eye, a green ring (**looks at**), and a thick arrow between them. Orbit the scene to see the path from other sides.
- With keys that look at atoms, the look-at path follows them through the trajectory. It is computed a few points at a time (frames must be read), so it grows for a moment; the old path stays until the new one is ready.

### Three paths

The viewport shows the paths of the **camera** (blue), what it **looks at** (yellow) and what is **in focus** (magenta, where a key sets it), each with a handle at every key that sets it. **Drag** a handle to move that part of the key only (a tracked look-at or focus becomes a fixed point), **Ctrl** to move the eye and the aim together; release to apply, Esc cancels. **Ctrl + click** on the path adds a camera key there. The spin is made about the look-at point.

### Global depth of field (ordinary viewport)

Under **Visuals > Depth of Field**, focus is always **Auto (look-at)**: the camera's look-at/orbit point defines the sharp depth. Double-clicking an object resets that point and thus focus. There is no global Fixed mode. Old global fixed/target modes are converted to Auto on load, with a notification; movie focus keys are unaffected. Blur settings are saved with the workspace and remain separate from movie targets.
Target tracking belongs to **Movie > Camera > Focus on**, not global settings.
Legacy global target modes load as Auto with a notification. Old Looks **Focus distance**
keys are retained but inactive; use a key's **Focus on > Distance** instead.
Existing movie keys without an explicit focus now mean look-at, not inherited focus;
repeat an explicit target on later keys if it should continue through them.

**Blur Strength** uses the upstream renderer's aperture: blur of distant objects as a percentage of the view height, independent of zoom level. The movie's **Depth of field blur** keys use the same units.

### A small FOV and focus tutorial

Open [fov_and_focus.via](examples/fov_and_focus.via) and press the green **Preview** button. Its companion [fov_and_focus.xyz](examples/fov_and_focus.xyz) is a tiny synthetic trajectory, included in the repository: no aspirin files are needed. Keep the two files together.

- **0-18 s, field of view:** all atoms and the camera stay fixed; blur is off. FOV changes from **45 to 20 degrees** (objects grow), then to **70 degrees** (objects shrink), then back to 45. This is lens zoom, not camera travel. The Camera table shows the changing FOV while distance stays at 32 A.
- **18-32 s, depth of field:** FOV, camera position and aim stay fixed. A movie-local focus key turns blur on over two seconds and tracks the centre carbon. The red foreground oxygens become sharp when it reaches their depth (22.5 s), then the blue background nitrogens become sharp at their depth (31.5 s).
- **32-34 s:** focus racks from the moving carbon to a fixed point on the red oxygen over one second.
- **34-36 s:** a no-blur key fades blur to zero over one second. All atoms become sharp again.

Open **Focus on** of the keys at 18, 32 and 34 s: the key at 18 s stores the carbon (atoms), the key at 32 s the red oxygen's fixed position, and the key at 34 s goes back to what the camera looks at. The blur, 0 to 3.5 % and back, is keys of **Depth of field blur** in the Looks tab. The example is in the format of the earlier separate tracks: it is converted to keys when it is loaded.

**In short:** FOV changes how much you see and the apparent size; depth of field changes what is sharp. Narrow FOV does not by itself enable blur. Increasing **Blur Strength** strengthens depth separation, not zoom.

## Looks

![Look parameters and overlays](images/movie/looks_and_overlays.png)

Background color and intensity, ambient occlusion, exposure, depth-of-field blur, near and far clipping planes and focus distance can change during the movie.

After the upstream renderer update, occlusion radius is no longer adjustable. Old occlusion-radius and depth-of-field-blur keys are retained in workspaces but inactive, with a message when loaded: their old units cannot be converted reliably. Re-key **Depth of field blur** using the new aperture units. Other parameter IDs are unchanged.

1. Pick the parameter in the list.
2. Set it up in **Visuals** at some preview time.
3. **Key Now**.

The table lists keys with time, value, ease and **Remove**. With **Animate parameters** on, keys are followed when scrubbing, previewing and recording; afterwards the values go back to what they were.

## Representations

### Keying a property

Pick a representation and a property, set it up in the **Representations** window as it should be at the preview time, and **Key Now**.

- **Visible:** shown or hidden. Solid representations cannot fade, so they grow in or shrink away instead: at a key they start going to the new state over **Transition (s)** (default 2, set above the table; 0 is instant), their sizes scaled during that time and hidden at zero. Key it once to appear and once to go. Keys closer than the transition turn around from where it was. Electronic-structure and dipole representations cannot be scaled: they are shown fully until the transition ends and then vanish.
- **Scales** (**Radius scale**, **Ball scale**, **Bond scale**, **Width**, **Coil**, ... whatever the type has), **Tint scale** and **Saturation** move smoothly with the usual **Ease**.
- **Base color** (visible where color mapping is Uniform) and **Tint color** (visible where **Tint scale** > 0), edited in the table; in the lane their keys are squares in their color.

A key belongs to a representation by an id that survives moving, duplicating or removing representations and is saved in the workspace; removing a representation removes its keys. **Animate parameters** switches these on with the look keys. Tint and saturation recolor the atoms every frame, which is slow for very large systems; scales and visibility are cheap.

### The overview lane

The **Representations** lane (the representation overview) has one band per system (protein, ligand, water, ...). A colored block is a stretch where one representation is shown; colors identify types. A block that starts where another ends stays on the same line and the two cross-fade; only representations shown at the same time stack. **Layout... > Lines** sets how a system is drawn: **Compact** (default), one line per system with overlapping blocks sharing its height; **Stacked**, another full line only where blocks overlap; **Separate**, a line for each representation, grouped by system. Stacked and Separate make the lane grow unless **Fit to window** is on. **Block height** (×0.5 to ×2.5) makes the blocks taller or flatter, and the lane with them; it has no effect with **Fit to window**. A representation enabled without keys has a block over the whole movie.

- **Drag** a bar to move it, or an end to change when it starts or stops. Times land on frames with **Snap to frames**; a bar stays between its neighbours.
- **Double-click** empty space in a system row to add one of its existing representations (a tenth of the movie long, up to that representation's next block); it may overlap others.
- **Right-click** a block: switch that interval to another existing representation of the system, or **Remove block**. Switching keeps other systems and property keys; overlapping target blocks merge. Removing the last block leaves the representation hidden.
- Clicking a block selects its representation (and **Visible**) for the **Representation lane**.

**Systems.** Representations named alike up to the first hyphen share a row: `protein-cartoon`, `protein-cpk`, `protein-vdw` belong to `protein`. A plain `protein` joins the same row (its block is labelled with the type). Other names without a hyphen get a row of their own. Rename in the Representations window to group. Switching a block chooses an existing representation; it does not make a new selection.

**Swap.** **Swap with the next of its group** hands over at the preview time from the selected representation to the next one of its group with the transition, e.g. `protein-cartoon` shrinks away while `protein-cpk` grows in. A representation without keys is taken to have been shown (the one that goes) or hidden (the one that comes) until then.

## Overlays

Overlays are drawn in recorded frames and in **Movie preview**. In the **Overlays** tab use **Add overlay...** and pick a type. Select one in the list to edit it: **Content** (type-specific), **Timing** (range and fades), **Appearance** (position, size, colors). Clicking an overlay bar in the timeline selects it too. A new movie starts with the VIAMD logo in the top left for the whole movie; **Remove** it if unwanted (the workspace remembers).

### Common to all overlays

- On/off box, **Type**, **Duplicate** (a copy right below: handy to show something else at another time or place) and **Remove**.
- **Shown (s)**, with **Start at preview time**, **End at preview time**, **Fade in (s)** and **Fade out (s)**.
- **Position:** one of nine anchors.
- **Size:** a part of the frame height, or in points (a point is a pixel of a 1080-high frame, scaled with the frame so it looks the same at any resolution). Text: its height; logo and image: their height; timeline and distribution: the whole stack.
- **Color** and **Background** plate (none while its opacity is 0).
- Overlays shown at the same time never overlap. Those in the top and bottom rows keep their place before those in the middle row, and the centre column before the sides. One that would hit a placed one moves past it: up from the bottom row, down from the top row, the shorter way in the middle row. So a bottom-left time bar sits above a wide bottom-centre timeline, and a tall middle-right timeline moves down a little under a top-right distribution. An overlay is never pushed out of the frame; when there is no room it stays where it is. Overlays shown at different times do not push each other, and an overlay keeps the same place for its whole range.

### Kinds

- **Text:** the words you type.
- **Time stamp:** time gone since the movie started, in the timeline's unit (**N frames** if the trajectory has no times). It starts at 0 and only grows, whichever way the trajectory plays, and reads the same whenever it is shown (one appearing at 40 s reads the time since the start). A hold does not advance it. Same quantity as the time bar's **Time that has gone**.
- **Scale bar:** a bar of a known structure length, as long on the frame as it is at the distance the camera looks at, so it follows the zoom. **Length** 0 picks 1, 2 or 5 times a power of ten.
- **Logo:** the VIAMD icon on transparency. Its light grey is made black for contrast on a light picture; on a dark one give it a light **Background**. **Color** multiplies the picture (white keeps it).
- **Image:** a png or jpg (a group logo, a figure), shown like the logo with transparency kept. The path is saved relative to the workspace, so keep the file with it. **Reload** rereads it; the box turns red when it cannot be read.
- **Time bar:** how far the trajectory has gone. It fills left to right whichever way the trajectory plays, fast where it plays fast, slowly where slowed and not at all in a hold, so it shows the pace. **Width** is a part of the frame width, **Size** the text height. **Time that has gone** writes the time covered over the whole; **Speed** the multiple of the Animation window's speed. The unfilled track is neutral grey; new time bars have a translucent grey plate (**Appearance > Background**); old workspaces keep their saved plate.
- **Property:** the visualization of a script property (atoms, geometry and the labels hovering its plot shows) while the overlay is shown, with fades, in preview and in the recording. Type the identifier or pick it from the script. The script must have evaluated the property.
- **Timeline** and **Distribution:** subplots of the **Timelines** and **Distributions** windows drawn into the movie (see next section). A timeline is wide and low at the bottom center, a distribution narrow and tall at the middle right; both can be moved and resized.

### Timeline and distribution overlays

*Which subplots.* A new overlay holds the subplots of its window that have series. **Add subplot** adds one, **Up** moves it up, **Remove** takes it out; several are stacked, each with its own scale and legend. What is drawn is what the subplot shows in its window, in its colors, at most six series per subplot and the first member of a population. A subplot is found by an identity of its own, not its position, so adding or naming subplots does not change what an overlay draws. Name them in the **Subplots** menu (**Names**) of their window. A subplot that is gone or empty is reported in the overlay and not drawn.

*As the movie plays.* Plots grow with the movie (**As the movie plays**; off shows the whole plot). A line and a dot mark the current frame; **Value** writes the value there, in the plot's unit, in the legend.

- **Timeline.** The horizontal axis is by default **Elapsed time**: trajectory time covered by the movie, as in the time bar. The curve always grows to the right, also where the movie plays backward (the series is read backward; at a turn the curve goes on rightward; a hold adds nothing). **Trajectory time** is the trajectory's own time, turned around when played backward. The axis spans the whole movie, and the curve is drawn from the start up to now, so a timeline that appears late already shows what happened before, with the same labels as ever. The overlay's own range only decides when it is visible. Timelines stacked above each other share the axis, labelled under the last.
- **Distribution.** Counted over the frames played since the movie started, so its bars grow to the shape of the whole, also when it appears late. A distribution of a script that is not over frames (an RDF) is drawn as it is.
- **Number of bins** (2 to 4096) applies to all subplots of a distribution overlay without changing the Distributions window; new overlays use 128. **Use source bin counts** keeps each series' own setting (the default for older workspaces). Script distributions can only be coarsened to a divisor of the resolution they were evaluated at.

*Subplots that come in at different times.* Each subplot has **in at** and **out at** (movie seconds within the overlay's range; **stays to the end** when equal), and **In at preview time**. A subplot added while the preview is later in the overlay comes in there. It fades with the overlay's fades. A property that comes in later is drawn from the start of the movie as if it had always been there: a curve that comes in at 30 s appears with what it did up to then and keeps growing; a distribution has counted every frame so far. The plate covers only the subplots that are there, so it grows when one comes in, and the place of a subplot not there yet stays free, so the others do not move. For a different place, use another overlay.

*Titles and look.* **Titles** writes a title above each subplot: the **Title** in its row, else the subplot's name, else none. **Text (points)** and **Lines (points)** set sizes in points (0 follows the height). **Colours** gives series colors of their own instead of those in the plots, and two buttons set light text on a dark plate or dark text on a light plate. **Width** is a part of the frame width.

*Markers.* Expand **Timeline markers** under the overlay inspector. **Add marker at the preview time** or **Add one at each camera key**, then set time, label and target. **All timeline subplots** puts the line and label on every timeline; choosing a subplot puts them on that subplot only, in any timeline overlay containing it. Targets use stable subplot identities, survive renaming and save/load, and a missing one is reported, not reassigned. **Markers** in an overlay hides its markers. Marker times scale with the movie length. Each marker has a color of its own: the swatch before its time picks one, right click goes back to automatic, which gives markers different colors from a palette in time order. Its line and label in the overlays and its triangle in the ruler and overlay lane have that color. Older workspaces load with automatic colors. A label starts at its marker and runs right in the first half of the plot; in the second half it ends at its marker and runs left, so labels near the end stay readable. A label that would run into one written before (on the same subplot) goes on a lower row (up to three, as room allows); if every row is taken it goes on the other side of the marker, and only if that is taken too does it overlap, where it overlaps least. Labels never run past the plot's edges. The marker labels in the editor's ruler are placed the same way, on lines above the frame numbers.

## The timeline

![The Movie timeline](images/movie/movie_timeline.png)

The lanes share the time axis; the yellow line is the playhead. Above them, **Lanes** is a row of toggles, one per lane (lit when shown): **Ruler**, **Camera**, **Trajectory**, **Lens**, **Representations**, **Overlays**, **Rep keys** and **Look**, in the order the lanes are drawn. On the line below, **All**, **Camera work** and **Scene work** are presets and **Layout...** opens lane height, **Fit to window**, **Snap to frames**, **Show whole movie** and the representations lane's **Lines** and **Block height**. When the **Rep keys** or **Look** lanes are on, a second row picks what they show. Numbers on the lanes' axes and in the mouse readout are short: at most two decimals, without trailing zeros.

New/reset movies use **Fit to window**, the **All** lane set, and **Depth of field blur** as the Look parameter. Saved workspace layouts take precedence. Presets are:

- **Camera work:** Ruler, Camera, Trajectory, Lens.
- **Scene work:** Ruler, Trajectory, Representations, Overlays, Look; selects Depth of field blur.
- **All:** Ruler, Camera, Trajectory, Lens, Representations, Look.

Overlays and Rep keys can still be toggled individually.

### Lanes

- **Ruler:** the time axis with a tick per frame once frames are a few pixels apart, the movie's markers as yellow triangles (click to go there), and the preview time. Click or drag in it to move the preview time. Every lane has its name in its top left corner.
- **Trajectory:** the frame shown at each moment. Drag the blue **Start** and **End** lines. Orange dots are keys with a pinned frame: drag sideways for time, up or down for frame.
- **Camera lane:** the camera on one lane. Before the first key and after the last one the lane is hatched grey: the camera holds still there. The first row has keys as orange dots with number and name (yellow at the preview time). Drag sideways to change when, click to go there, right-click for **Go to** and **Remove** (and the name), double-click an empty spot to add a key there on the path, without moving the camera. Row **Look at**: a blue band over keys that look at atoms ("3 atoms", "atom 12"). Row **Focus**: a magenta band from each key that sets the focus to the next one that does, its transition darker, with what is sharp and the blur. Row **Spin**: a purple band over stretches with extra turns (`+2 x`). Row **Frame**: a green mark, with the frame number, at each key that pins a frame. A key added by double-click looks at the atoms its two neighbours look at; a spin's turns stay with the part ending at the later key.
- **Lens:** field of view (left axis, degrees, 1 to 170) and camera distance (right axis) on one plot, with a legend. Drag dots sideways for time and up or down for value; each dot belongs to its curve's axis (distance moves the camera along its line of sight).
- **Look parameter lane:** the parameter chosen at the top. Double-click adds a key, drag changes it, right-click removes.
- **Representation lane:** the property chosen with the two lists next to it. The line is the value, dots are keys; same editing. For **Visible**, double-click turns it around at that time; the first one on a representation without keys also keys its present state at 0 so it holds until then.
- **Representation overview:** see [above](#the-overview-lane).
- **Overlay lane:** overlays as bars (one row and color per kind), markers as small yellow triangles. Drag a bar to move it, an end to change when it starts or stops. Moving a whole bar moves the in and out times of its subplots; small ticks in a bar show where subplots come in. Hover for name and times.

### Picking and moving things together

Camera keys, trajectory, lens (distance and FOV), look-parameter and representation-lane keys, overlay bars and overview blocks can be picked and moved as a group.

- **Pick:** click an item (a white ring shows it). **Ctrl + click** adds or removes. **Drag on a lane's empty background** draws a box (Ctrl or Shift: adds; in the camera lane it picks by time only). A click on the background puts all down. **Pick all** (Ctrl + A) picks the camera keys and the shown look and representation keys. A camera key is one key in every lane it appears in. Clicking a key's handle in the viewport picks it too.
- **Move:** drag a picked item and all move by the same time and, in lanes with a value, by the same amount (frame, distance, FOV, parameter value, property value; colors and **Visible** move in time only). A tooltip shows how far (`5 keys: +1.20 s  +0.5`). The group stays inside the movie, keeping the spacing, and lands on frames with **Snap to frames**. Dragging an unpicked key picks it alone. One drag is one undo step. A key landing on another in the same lane replaces it on release.
- **Bars and blocks:** dragging a picked bar or block by its middle moves everything picked by the same time, and takes along what is timed inside a bar. A block stops at other blocks of its representation. Dragging an end still changes only that end. **Delete** removes picked blocks and keys (never an overlay). Unkeyed blocks cannot be picked.
- **Stretch:** the percentage box next to **Copy** stretches the time between picked items about the first one or the preview time (the list beside it). 100 % is unchanged. It stops where something would leave the movie. One undo step.
- **Inspector:** a line above the lanes says how many of each kind are picked. **first at** moves all so the first is at the typed time. The ease list sets how the movie gets to every picked key (**(mixed)** when they differ). One picked camera key can be renamed there.
- **Keyboard** (mouse over the lanes): left/right move picked keys one frame (Shift one second, Ctrl ten frames); up/down change the value of the lane last used (Shift ten steps); **Delete** removes; **Esc** puts down; **Ctrl+C** / **Ctrl+V** copy and paste with the first at the preview time. **Pick all**, **Put down**, **Delete**, **Copy** and **Paste** above the lanes do the same.

### Zoom, panning and layout

- Ctrl + wheel zooms time; **Show whole movie** resets. Because dragging on the background picks, pan with the **middle mouse button**, **Shift + wheel** or sideways scrolling.
- **Lane height** (60 to 400 px) times each lane's proportion (drag between lanes to change it). The window scrolls with the wheel when lanes do not fit; the overview, overlay and camera lanes grow to fit their rows without making others taller. **Fit to window** shares the window height among the lanes instead, and then the wheel zooms.
- The **Lanes** toggles hide or show a lane; the settings below are in **Layout...**.
- **Snap to frames** (default on): dragged keys, the playhead, the trajectory start and end, and keys added at the preview time land on a frame at the **Output FPS**.

## Saving

The movie (length, trajectory timing, camera keys with their names, look-at and focus, look and representation keys, overlays with bin counts and markers with subplot targets, path options and timeline layout) is saved in the workspace (`.via`) under `[Movie]`, with the subplot names the overlays use. Workspaces from earlier versions load and are converted. Whether Scene view or Movie preview is shown, and which editor panels are visible, are not saved.

## Example workspaces

[`examples/aspirin_phospholipase_movie.via`](examples/aspirin_phospholipase_movie.via) is an 80 s movie of aspirin in phospholipase, with the trajectory played backward. It has camera keys (the first look at the aspirin), keyed depth-of-field blur, representation keys, two markers ("Pocket shown", "Aspirin enters the pocket"), a logo, title, time stamp, scale bar, a time bar with a grey plate, four timeline subplots that come in at different times, a 128-bin distribution with a title, and three property visualization overlays. It opens at the saved preview time (about 55.45 s): set the scrubber to 0 and play.

[`examples/aspirin_binding_movie.via`](examples/aspirin_binding_movie.via) is a 60 s movie of the same data, built to show what the tools can do together. The trajectory is played in reverse, so aspirin approaches and docks. It is in five chapters:

1. **The approach** (4 to 16 s): the trajectory plays fast, the camera glides in, the waters around aspirin come in.
2. **First contact** (16 to 24 s): a **dolly zoom** at touchdown. Over 7 s the camera backs away from 14 to 80 A while the field of view narrows from 76 to 16 degrees, so aspirin keeps its size and the protein behind it seems to rise. The cartoon hands over to a spacefill surface that loses its color behind aspirin.
3. **Bullet time** (25 to 30.5 s): the trajectory is frozen at the frame with the strongest protein-aspirin Coulomb attraction (frame 288) while the camera goes once around aspirin (a spin of one turn, slightly tilted, with a tinted background and shallow depth of field). The surface changes to licorice, then the pocket residues grow in.
4. **Into the pocket** (31 to 44 s): the cartoon returns and the camera tracks aspirin down into the pocket.
5. **The calcium handshake** (47 to 57 s): in slow motion the carboxylate meets the Ca2+ ion, which swells at the strongest attraction; the waters leave and aspirin changes to spheres. The camera then pulls back for the end title.

The plots are the short-range Lennard-Jones and Coulomb terms from the `.edr` between aspirin and the protein, the Ca2+ ion and the water, and the geometry (distances to the pocket and to the ion, the number of waters around aspirin, the torsion of the carboxylate), each coming in at the time it matters, with four markers. The energies are read from the `.edr` with `attr("edr/...")` in the workspace's script.

The reasoning behind its choices (the physics of the energy terms, the placement of aspirin and of the camera, the dolly zoom, the bullet time, the representation hand-overs) is written up in [`aspirin_binding_movie_article.md`](aspirin_binding_movie_article.md).

Both refer to `aspirin-phospholipase.gro`, `.xtc` and `.edr` in their own folder, which are not in the repository. Put your copies next to the workspace and use **File > Open Workspace**.

## Tips

- Plan the length first, then place keys. If the pace is wrong, change **Movie length**: everything stays in proportion.
- For a fly-over before the dynamics, drag the blue **Start** line right.
- To slow an event, pin frames on two keys around it and move them further apart.
- If a loop or spin tilts the camera, tick **Keep upright** (use **From view** if up is not +Y).
- **Smooth** keys never stop the camera; use **Ease in/out** or **Hold** where it should rest.
- PNG sequences are safest for very long or large movies; encode them afterwards with **Copy ffmpeg command**.
- Check motion at **Scale** 25 % with few **Samples per frame**, then render at full size. To fix one part, use **Render only a range** with a PNG sequence and replace the files.
- Name representations after their selection (`protein-cartoon`, `protein-cpk`) so the overview keeps them together.

## Known limits

- Legacy global target-focus modes convert to Auto. Use a key's **Focus on** to track an object in a movie.
- A timeline or distribution overlay draws at most six series per subplot and the first member of a population, and a timeline overlay has one axis kind for all its subplots.
- Solid representations cannot fade; they grow in and shrink away.
- The Representations window is locked while recording.
- With **Keep upright**, a spin around a horizontal axis flips over at the top and bottom (the camera cannot stay level looking along the up axis). Roll eases as a number, so 170 to -170 degrees turns the long way round.
