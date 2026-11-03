# Correlation

The **Correlation** window plots one time series against another, frame by frame, in a single scatter plot. It shows how
two quantities of a trajectory move together (a distance against a hydrogen bond, an angle against the size of the
chain, ...), with the strength of the correlation and how far it can be trusted.

It works like the [Ramachandran](Ramachandran) window: a stack of **layers** you define, each drawn as points,
contours or a heat map. Unlike Ramachandran, there is one plot instead of four, no reference layer, and you choose
what goes on each axis.

![The Correlation window: layers on the left, the plot, its caption and the lag view](images/correlation/overview.png)
<!-- IMAGE PLACEHOLDER: overview.png — the whole window with the alanine demo, layer list, plot with contours, caption strip and lag view visible -->

---

## Contents

- [Opening the window](#opening-the-window)
- [Quick start](#quick-start)
- [Layers](#layers)
- [Putting series on the axes](#putting-series-on-the-axes)
- [Pairing: 1 × 1, 1 × n, n × 1, n × n](#pairing-1--1-1--n-n--1-n--n)
- [One kind of axes per plot](#one-kind-of-axes-per-plot)
- [Layer settings](#layer-settings)
- [The plot](#the-plot)
- [Periodic axes (angles)](#periodic-axes-angles)
- [The caption](#the-caption)
- [Statistics](#statistics)
- [The lag view](#the-lag-view)
- [Mouse and keyboard](#mouse-and-keyboard)
- [Settings](#settings)
- [Example: the alanine workspace](#example-the-alanine-workspace)
- [Limits](#limits)

---

## Opening the window

**Windows → Correlation**. The window and all its layers are saved with the workspace (`.via`).

---

## Quick start

1. Define a few properties over time in the script editor, for example a distance and a dihedral.
2. Open **Windows → Correlation**. A first layer is already there.
3. Drag a series from the **timeline** onto the **X** strip at the bottom of the plot, and another onto the **Y** strip
   on the left. You can also click **+** next to X or Y in the inspector and pick from the list.
4. The points appear, one per frame. The caption under the plot gives the correlation.
5. In the inspector, set **style** to *Contours* to see where the trajectory spends its time, or **frames** to a
   timeline period to look at part of it.

![Dragging a series from the timeline onto the X strip](images/correlation/drag_drop.png)
<!-- IMAGE PLACEHOLDER: drag_drop.png — a series being dragged from the timeline, the X drop strip highlighted with "Drop on X" -->

---

## Layers

The sidebar on the left (toggle it with the ![layers](images/correlation/icon_layers.png) button in the toolbar) lists
the layers. Each layer has its own series on X and Y, its own frames and its own look. They are drawn bottom up: the
first layer of the list is on top.

| On a layer card | Does |
|---|---|
| Click | Selects it and shows its settings in the inspector below the list |
| Eye | Shows or hides it |
| Right-click | **+ Series on X / Y**, **Duplicate**, **Swap X and Y**, **Remove** |
| Spinner | Its density is being computed |

**Adding layers.** The **+** button next to *LAYERS*, at the top of the list, adds a layer over the *current frame*, *all frames*, or a *period*
marked in the timeline. A new layer starts with the selected layer's series. Dropping a timeline **period** on the
plot does the same for that period, and dropping it on a layer card sets that layer's frames to the period.

![The layer list grouped by axes, with the inspector of the selected layer](images/correlation/layers.png)
<!-- IMAGE PLACEHOLDER: layers.png — the sidebar with two axes groups (e.g. "Å × Å" and "° × °"), one dimmed, and the inspector below -->

---

## Putting series on the axes

A series can be any property **over time**:

- the script's temporal properties (distances, angles, dihedrals, RMSD, ...);
- the series loaded with the trajectory (energies, temperature, ...).

There are two ways to add one to the selected layer:

- **Drag and drop** from the timeline onto the **X** or **Y** strip along the plot's axes, or onto the X / Y list in
  the inspector.
- **Pick** it: click **+** next to the X or Y list in the inspector (or right-click the layer card → **+ Series on X / Y**).
  The list has a search box and is sorted into sections (script, then each group of loaded series). Series already on
  the axis are ticked.

Click ✕ next to a series in the inspector to remove it. **Swap X and Y** exchanges the two axes.

A series is **refused**, with the reason shown, when:

- it is not over time (a distribution, for example), or is a range of values;
- its unit or periodicity differs from the series already on that axis. Put it on a layer of its own instead.

![The series picker opened from the + of the Y list](images/correlation/picker.png)
<!-- IMAGE PLACEHOLDER: picker.png — the picker popup with its search box, sections, a ticked series and a refusal in red -->

> Distributions are not supported: only series over time can be correlated.

---

## Pairing: 1 × 1, 1 × n, n × 1, n × n

A layer can hold up to 16 series on each axis. They are paired like this:

| X | Y | Pairs |
|---|---|---|
| 1 | 1 | X against Y |
| 1 | n | X against each Y (e.g. one size against several H-bonds) |
| n | 1 | each X against Y |
| n | n | X₁ against Y₁, X₂ against Y₂, ... in order (e.g. φ₄/ψ₄, φ₈/ψ₈, φ₁₂/ψ₁₂) |

Any other combination (2 against 3, ...) is refused with a message in the inspector.

**Populations.** A series with one value per member (a distance for each of several atom pairs, an RMSD per residue)
is taken member by member. Against a single value, every member is paired with it. Against another population, members
are paired by index, up to the smaller of the two.

Each pair, and each member of a single population, gets its own colour, statistics, caption entry and fit line, so a
1 × n comparison shows n correlations within **one** layer.

![A 1 × 3 layer split by pair: three densities in three colours](images/correlation/split_pairs.png)
<!-- IMAGE PLACEHOLDER: split_pairs.png — "Collapse vs H-bonds" from the demo, contours split by pair, caption with three entries -->

---

## One kind of axes per plot

The plot shows one kind of axes at a time: one unit on X and one on Y, periodic or not. Layers are therefore **grouped
by their axes** in the list (for instance *Å × Å*, *° × °*, *Å × °*), under a header naming them. The group the plot
shows is highlighted; the others are dimmed.

- **Clicking a layer** of another group turns the plot to its axes and shows that layer. The layers of the previous
  group are hidden.
- Turning on the **eye** of a layer of another group does the same.
- Clicking a group header shows its first layer.
- Within the shown group, turn on as many layers as you like to overlay them.

This keeps the axes meaningful: a distance is never drawn on the same axis as an angle.

---

## Layer settings

The inspector shows the selected layer's settings.

| Setting | Meaning |
|---|---|
| **name** | Shown in the list and the caption. Empty: named after its series |
| **X**, **Y** | The series on each axis (drop, or **+** to pick) |
| **frames** | *Current frame* (one point per pair, follows the animation), *All frames*, or a timeline *period* |
| **style** | *Points*, *Contours* of the density, or *Heatmap* of the density |
| **color** | The layer's colour |
| **opacity** | Of the whole layer |
| **draw** (contours) | *fill* and/or *lines* |
| **split** (contours) | *by pair*: one density per pair, each in its own colour. *by member*: likewise, per member of a single population (32 at most) |
| **line width** (contours) | Of the contour lines |
| **levels** (contours) | The fractions of the samples the contours enclose (highest density regions) |
| **colormap**, **log scale** (heatmap) | Colours of the heat map; log scale spans four decades below the maximum |
| **bandwidth** | Smoothing of the density, in % of the extent of the samples |
| **point size**, **color by** (points) | Colour by *layer*, by *pair* (or member), or by *time* along a colour map |
| **fit** | The least squares line of each pair (not on periodic axes) |

Under the settings, a table lists each pair with its number of samples, *r*, its 95 % interval, Spearman's ρ and the
slope.

A density needs several frames: over the current frame, or without a trajectory, the layer is drawn as points.
Densities are computed in the background; a spinner in the plot corner shows when one takes a while. Past 100 000
samples, only every n-th one is drawn as a point (noted in the corner), but all of them count in the statistics and
the densities.

![The three styles: points coloured by time, contours, heat map](images/correlation/styles.png)
<!-- IMAGE PLACEHOLDER: styles.png — the same layer in Points (color by time), Contours and Heatmap, side by side -->

---

## The plot

- The axis labels and tick labels are drawn **inside** the plot, as in Ramachandran (tick labels can be turned off in
  the settings).
- The **X** and **Y** drop strips run along the bottom and left edges.
- **Double-click** (or the ![fit](images/correlation/icon_expand.png) button) fits the view: on a linear axis, to all
  but the farthest 0.5 % of the samples either way, so a few outliers do not squash the plot.
- When both axes are of the same kind (angle against angle), they keep equal scales.

---

## Periodic axes (angles)

An axis of periodic values (dihedrals, angles in ° or rad) behaves like Ramachandran's:

- The points and densities **repeat** one period apart beyond the edges, and thin lines mark where the circle is
  **cut**.
- Panning goes round the circle; you can zoom out to at most two periods.
- **Right-click a periodic axis** to move the cut. The cut is shared with the timelines of the same values: moving it
  here moves it there, and vice versa. Panning or zooming here does not move the timelines.
- Double-click fits a whole turn starting at the cut.

![A φ/ψ layer on periodic axes, with the cut lines and repeated images](images/correlation/periodic.png)
<!-- IMAGE PLACEHOLDER: periodic.png — "Phi/psi of residue 8" from the demo, contours wrapping across the cut, the right-click cut menu open -->

---

## The caption

Under the plot, the **caption** lists each pair (or member) of the shown layers with its colour, name and correlation,
for example `r 0.82 ±0.05`. Correlations that cannot be told apart from zero are dimmed. It scrolls past three rows and
can be turned off in the settings.

**Hover an entry** to make it stand out: the other pairs fade in the plot and the lag view, and the atoms of that
pair are highlighted in the 3D viewport. Hovering also shows the full statistics in a tooltip.

![Hovering a caption entry: its pair stands out, the others fade](images/correlation/caption_hover.png)
<!-- IMAGE PLACEHOLDER: caption_hover.png — the caption strip with one entry hovered, its tooltip, and the plot dimmed except that pair -->

---

## Statistics

For each pair:

| Statistic | Meaning |
|---|---|
| **Pearson's r** | Linear correlation, −1 … 1 |
| **95 % interval of r** | By Fisher's z, over the **effective** number of samples (see below). Dimmed when it includes 0 |
| **Spearman's ρ** | Rank correlation: how monotonic the relation is, whatever its shape |
| **slope** | Of the least squares line, Y on X |

**Periodic values** use circular correlations instead of Pearson's r:

- **angle against angle**: Jammalamadaka–SenGupta's circular correlation, −1 … 1, marked **rₒ**;
- **angle against a linear value**: Mardia's circular–linear correlation, 0 … 1 (unsigned). Its tooltip also gives the
  signed correlations with the cosine and sine of the angle.

**Effective number of samples.** Successive frames are alike, and members of a population move together, so N frames
are worth fewer than N independent samples. The effective number is estimated from the lag-1 autocorrelations of X and Y
(Bretherton et al., 1999). It is what the confidence interval is based on, so long, slowly varying trajectories do not
give falsely tight intervals.

---

## The lag view

The ![lag](images/correlation/icon_lag.png) button in the toolbar opens the **lag view** under the plot. For each pair,
it plots the correlation of Y at frame *f + lag* against X at frame *f*, for lags up to a bound either way (set in the
settings, 1000 frames at most).

- A peak at a **positive** lag: Y follows X.
- A peak at a **negative** lag: Y leads X.

It needs layers over several frames. Hover a curve to read off the lag and r; hovering a caption entry highlights its
curve. For very large populations, only some members are used (noted in the corner).

![The lag view: correlation against lag for three pairs](images/correlation/lag.png)
<!-- IMAGE PLACEHOLDER: lag.png — the lag view with three curves and the hover readout -->

---

## Mouse and keyboard

| Action | Does |
|---|---|
| Drag a series from the timeline onto X / Y | Adds it to the selected layer |
| Drop a timeline period on the plot | New layer of the selected layer's series over that period |
| Drop a timeline period on a layer card | Sets that layer's frames |
| Hover a point | Tooltip with its values, member, frame and time; shows it in the viewport |
| Click a point | Goes to its frame |
| **Shift**-drag / **Shift**-click | Selects the atoms of the boxed points' series (as Ramachandran selects residues) |
| **Shift**-right-drag | Deselects them |
| Right-click on nothing | Clears the selection |
| Right-click a periodic axis | Moves its cut |
| Double-click | Fits the view |
| Hover a caption entry | Highlights its pair, its lag curve and its atoms |

---

## Settings

The ![settings](images/correlation/icon_gear.png) button at the top right:

| Setting | Meaning |
|---|---|
| **Caption** | Show the caption under the plot |
| **Tick labels** | Show the tick labels inside the plot |
| **Lag up to** | How far either way the lag view goes, in frames |

---

## Example: the alanine workspace

`datasets/1ALA-500_correlation.via` loads a 15-residue polyalanine (500 frames) that collapses and starts to fold near
frame 350, with a script and six layers ready to explore:

| Layer | Axes | What to look at |
|---|---|---|
| Collapse vs H-bonds | Å × Å | End-to-end distance against three i→i+4 H-bonds (1 × 3, split by pair): which bond forms as the chain shrinks |
| Phi/psi of residue 8 | ° × ° | A Ramachandran-like map of one residue on periodic axes |
| Phi/psi along the chain | ° × ° | φ/ψ of residues 4, 8 and 12 paired (3 × 3): does the chain fold evenly? |
| Phi 8 vs size | Å × ° | Circular–linear correlation of a dihedral with the chain's size |
| Residue RMSD vs size | Å × Å | A population (RMSD per residue) against one value |
| Residue RMSD vs size, by residue | Å × Å | The same, contours split by member: which residues move with the collapse |

Try:

1. Select each layer in turn and watch the plot switch axes.
2. Open the lag view on *Collapse vs H-bonds* to see whether the bonds form before or after the collapse.
3. Mark a period after frame 350 in the timeline and drop it on the plot to compare the folded part with the rest.
4. Shift-drag over a cluster of points to select the atoms involved.

![The alanine workspace, layer "Residue RMSD vs size, by residue"](images/correlation/alanine.png)
<!-- IMAGE PLACEHOLDER: alanine.png — the demo workspace with the by-residue layer, viewport showing a highlighted residue -->

---

## Limits

| | |
|---|---|
| Layers | 16 |
| Series per axis of a layer | 16 |
| Members a pair can be split into | 32 |
| Points drawn per layer | 100 000 (more are thinned for drawing only) |
| Lag | 1000 frames either way |
| Distributions | Not supported |
