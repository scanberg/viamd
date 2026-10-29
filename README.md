# VIAMD
Visual Interactive Analysis of Molecular Dynamics

VIAMD is an interactive analysis tool for molecular dynamics (MD) written in C/C++. VIAMD is developed at the PDC Center for High Performance Computing (KTH, Stockholm). It exposes a rudimentary scripting language for declaring operations performed across trajectory frames.
The results can then be viewed in the different windows exposed in the application. 
<p align="center">
<img src="https://github.com/scanberg/viamd/assets/38646069/5651ef62-28bc-4f41-8234-75cf9ba85612" alt="This is an overview of the viamd software" width="800"/>
</p>

## Status
[![Windows (MSVC 19)](https://github.com/scanberg/viamd/actions/workflows/windows.yml/badge.svg?branch=master)](https://github.com/scanberg/viamd/actions/workflows/windows.yml)
[![Ubuntu 22.04 (GCC 11)](https://github.com/scanberg/viamd/actions/workflows/ubuntu22.yml/badge.svg)](https://github.com/scanberg/viamd/actions/workflows/ubuntu22.yml)
[![Ubuntu 24.04 (GCC 13)](https://github.com/scanberg/viamd/actions/workflows/ubuntu24.yml/badge.svg)](https://github.com/scanberg/viamd/actions/workflows/ubuntu24.yml)
[![MacOS (Clang)](https://github.com/scanberg/viamd/actions/workflows/macos.yml/badge.svg)](https://github.com/scanberg/viamd/actions/workflows/macos.yml)

## Running VIAMD 

### Windows
For Windows, we recommend using the latest binary available on the [release page](https://github.com/scanberg/viamd/releases/).

### macOS
For macOS, we recommend using the latest binary available on the [release page](https://github.com/scanberg/viamd/releases/).

**NOTE:** To be able to run this, you will have to do some magic in a terminal to be able to run the app bundle since we do not have code signing.

```xattr -dr com.apple.quarantine viamd.app```


### Building for Ubuntu and macOS
To [build](https://github.com/scanberg/viamd/wiki/0.-Building) VIAMD on your machine, you can follow the procedure described in detail in the wiki for [Linux](https://github.com/scanberg/viamd/wiki/0.-Building#linux) and [MacOS](https://github.com/scanberg/viamd/wiki/0.-Building#mac).

On macOS, VIAMD builds without Xcode: CMake and the macOS Command Line Tools are
enough. If Apple's Metal compiler happens to be available, the Metal compute
kernels are compiled to `.metallib` during the build; otherwise VIAMD embeds the
shader source and compiles it with Metal's runtime shader compiler on first use,
which costs a fraction of a second at startup and is otherwise identical. Xcode
is only needed if you are changing the offline Metal shader compilation pipeline
itself.

## Documentation
Documentation about VIAMD is available on the github [wiki](https://github.com/scanberg/viamd/wiki). The first two chapters relate to the [visual](https://github.com/scanberg/viamd/wiki/1.-Visual) and [analysis](https://github.com/scanberg/viamd/wiki/2.-Analysis) features, respectively, where we highlight the interactive part of the software. The third chapter focuses on the VIAMD [language](https://github.com/scanberg/viamd/wiki/3.-Language) used for scripting, and the fourth chapter proposes a series of [tutorial](https://github.com/scanberg/viamd/wiki/4.-Tutorials) (under construction). 

A series of videos is available on [youtube](https://youtube.com/playlist?list=PLNx9MpJY8ffr9CeK7WefdOnuGRw_E5rSj&si=VatBHEwiL7jWyhPK).

### Critical-point density distributions
Below the search controls and the category line, the Critical Points window has
two rows of two panes: the chemical environments beside the density
distributions, and below them the list of critical points beside the details of
one. Drag the bar between the rows to share the height, and the pane borders to
share the width.

The distributions cover attractors, bond critical points, ring critical points,
and cage critical points. Each panel uses the category's marker color and an
independent logarithmic electron-density axis (au, e bohr^-3). The filled curve
is an adaptive nearest-neighbor Gaussian estimate in log10 density, normalized
to unit area per decade, not a histogram or an electron density. Categories
containing zero or negative densities use a linear axis and probability per
density unit instead. Each curve is split into stacked colored layers, one per
chemical environment (the five largest in view, the rest pooled as Other environments);
the layers sum exactly to the category's total curve, drawn as a line on top.
There is no legend: the tables' environment swatches use the layer colors.

Bars along the bottom show the actual critical points. Bars within six screen
pixels are grouped, with a count label, so degenerate values remain visible.
Bar height is directly proportional to the group count within each panel, with
the tallest bar's count displayed as a scale reference. Bars whose CPs share one
environment take its color; selected and highlighted bars use the viewport's
selection and highlight colors.

Density plots follow the Distributions window's navigation: mouse-wheel zoom,
left-drag pan, and double-click to fit the density axis. The vertical axis
auto-fits; density-axis zoom persists until refitted or new search results arrive.
Dragging to pan does not select a bar or clear the selection.

Each plot's curve is split into peaks at its valleys, the deepest first
(persistence), into at most six. The environments table lists its rows in
blocks, one per peak: a strip in the peak's colour and its number (1, 2, ...)
runs down the block, and the rows are the chemical environments within that
peak, with their count and share of it. Sorting applies within blocks. Clicking
the strip selects the whole peak; clicking a row selects that environment.

Zooming refines rather than filters. A zoomed plot draws the curve of the CPs
in view at the view's resolution (the overview's curve stays as a faint line),
and its peaks split at that curve's valleys, again into at most six in view.
The overview's peaks are never merged, so each part lies within one of them and
is numbered after it (2.1, 2.2) in shades of its colour. Peaks outside the
view keep their overview split. Each peak's CPs are named at the finest level
of detail that keeps it in at most four environments: the group ("C
aromatic"), then with what the atoms are bonded to ("C aromatic [C2H]", "H (on
N [C2H])"; a contact ring's size and composition), then with their residues
where the system has them ("in PHE", "between GLY and SER"). The plot's shading,
both tables and the viewport's environment colors and tags follow; environments
selected in the table stay selected through the change, as their parts or as
the whole they merge into. The details pane always shows the finest name. The
CSV export includes each row's peak, its density interval and the share.

Each category has an independent, optional density filter made of one or more
inclusive ranges; a point passes if it lies in any of them (their union). The
filters are the chips above the plots: + Filter (or +) adds a range over the
middle third of the plot's current view (zoom onto a peak first); drag its two
vertical handles in the plot, click the chip to type Min/Max, or press its x to
remove it. The line shows how many CPs pass. Accepted ranges are shaded; excluded
bars are dimmed, and a partially accepted group has a bright segment
proportional to its passing count. Only points whose category is shown and
whose density passes its filter appear in the viewport or take part in viewport
picking, region selection, or view fitting; paths from excluded points are
hidden too. Filtering does not recompute the search or change the topology export.

The environments table lists, one flat row each, the chemical environments of
the CPs of the shown categories that pass their filter, at the detail their
plot's view gives them: Type, Kind, Chemical environment, count, density range,
mean Laplacian and mean ellipticity.

The critical point list below it shows the CPs of the environments chosen in
the table above, else of the bars clicked, else every CP the environments table
counts (List all returns to that), with the same first columns (Type, Kind,
Chemical environment, index) and then density, Laplacian, ellipticity and
position. Both tables sort by clicking a header, Shift-click further headers to
sort by several columns; right-click a header to hide columns (Type and
Chemical environment always stay).

The Export menu at the top writes two kinds of data. Topology (JSON or VTP) is
every critical point and path, unfiltered and without chemical labels, for scripts
and ParaView. Chemical analysis (CSV) is either table as listed and sorted, every
column at full precision (hidden ones too); the environments add their CP indices,
the CP list each CP's label detail, and the CP index joins it to the topology
export.

Selection works at both levels. Click an environment to choose it, Shift-click to
add or remove one: its CPs are listed below and become the selection. Click a
bar to list and select its CPs, Shift-click to add. In the list, click a CP to
select only it, Shift-click to add or remove it, Shift-right-click to remove it;
the list itself stays. CPs picked in the viewport join the same selection, and a
plain click on empty plot space (or viewport) clears it. Hovering only
highlights, everywhere, with no tooltips.

The last CP clicked in the list (or the only selected CP) is explained in the
details pane and ringed in the viewport. Show sets what the viewport does with
the selection: Selected CPs (the default: every other CP is hidden and cannot be
picked), All CPs (the other CPs only lose saturation, as selections do elsewhere
in VIAMD), or Selected + 3 Angstrom (other CPs within 3.0 Angstrom of a selected
one keep full color, the rest are dimmed to 15% opacity). Selected CPs carry a
selection ring. With nothing selected every CP is shown. Category visibility and
density filters still apply.

Next to Show, Environment colours fills the viewport markers with their chemical
environment's colour (the one of the plots and tables) instead of their
category's, and Labels names the environment of each selected CP beside it in
the viewport (at most 200).

The details give the CP's type, density, Laplacian (and for a BCP, shared-shell
or closed-shell), Hessian eigenvalues, ellipticity, position and why it got its
label.

Every label has a kind, shown in its own sortable column of both tables
and as a prefix elsewhere ("Bond: C aromatic – C aromatic"):

- Atom: an attractor with exactly one nucleus within 0.30 Angstrom, and no other
  attractor competing for that nucleus. Labels use the existing element,
  aromaticity, double-bond, and hydrogen metadata (C aromatic, C/O carbonyl,
  O/N (H-bearing)); hydrogens are named by the atom they are bonded to (H (on C)).
- Bond (–): a BCP whose two traced attractor endpoints are atoms bonded in the topology.
- Contact (···): the same, for atoms not bonded in the topology, classified from
  elements, connectivity and geometry, and tagged intramolecular or intermolecular
  when molecules are known:
  - Hydrogen bond: H on N/O/F/S to N/O/F/S/Cl/Br/I, with D–H···A above 110 degrees.
  - X–H···π / C–H···π: such an H, or an H on C, to an atom of an aromatic ring
    (same angle).
  - C–H···X: an H on C to N/O/F/S/Cl/Br/I (same angle).
  - H···H, Halogen bond (C–X···Y above 140 degrees, X = Cl/Br/I, Y = N/O/S),
    π···π (atoms of different aromatic rings, or sp2/sp atoms bonded to one, such
    as stacked carbonyls), else Other contact.
- Bond path (~): the same, when the system has no connectivity.
- Ring: an RCP within 0.50 Angstrom of the centroid of exactly one perceived ring
  (3-7 atoms), named by aromaticity, size and composition ("Aromatic 6-ring C6").
  Otherwise it is matched to a cycle of bond paths, taken as straight lines between
  nuclei, whose polygon encloses it within 1.0 Angstrom of its plane and whose BCPs
  are all denser than it (a ring CP is the density minimum of its ring's surface).
  Cycles of 3-10 atoms within 5 Angstrom are tried first, then, if none fits, of
  up to 14 atoms within 6 Angstrom (the larger hydrogen-bonded rings of proteins).
  Matching is global: smallest and closest pairs first, each cycle claimed by one
  ring CP. Between free cycles of the same size about as close, the one whose face
  holds no atom of the others wins, else the smallest by area; if the two smallest
  differ by less than 15% the ring CP is tied. With the analytic search, a tied
  ring CP's ring surface is then traced in the background (the gradient paths
  leaving it in the plane of its two positive curvatures, bisected down to the
  narrowest sectors) and the cycle holding every nucleus it reaches wins; if none
  or several do, the ring CP stays unassigned. Bonded
  cycles are named "Bonded 8-ring"; cycles closed by contacts are named by their
  strongest contact class ("Closed by Hydrogen bond", in the order Hydrogen bond,
  Halogen bond, π···π, X–H···π, C–H···π, C–H···X, H···H, Other contact), with size,
  composition and contacts in the details. A ring CP with a
  density within 2% of a BCP on its ring is marked near-degenerate: the pair is
  close to annihilating.
- Cage: a CCP whose ring-CP paths reach at least 3 faces, all assigned rings
  ("Bonded cage, 6 faces", or "Contact cage, 6 faces" if any face is closed by a
  contact); its composition and faces are in the details.
- Non-nuclear attractor: an attractor more than 0.50 Angstrom from every nucleus,
  named by its two nearest ones.
- Unassigned: everything else, named by its topology ("Unassigned attractor",
  "Unassigned BCP", "Unassigned ring CP", "Cage CP (unassigned)"). Attractors
  0.30-0.50 Angstrom from a nucleus stay unassigned. Missing paths and ambiguous
  rings are not replaced with nearest-atom guesses.

With the analytic search, the details also give the sign of a BCP's Laplacian
(negative: shared-shell, positive: closed-shell) and its ellipticity, which the
tables have as a sortable ε column.

Contact classes name geometry and topology only: they do not establish
interaction energies or stabilization, and do not modify the molecule, other
windows, or the topology export.

## Update
If you want to stay informed about the latest updates of VIAMD, please register your email address on the [form](https://forms.gle/fAxuWob8nMLcrS5h9). 

## Citations:
* General Framework:
  * R Skånberg, I Hotz, A Ynnerman, M Linares, VIAMD: a Software for Visual Interactive Analysis of Molecular Dynamics, J. Chem. Inf. Model. 2023, 63, 23, 7382–7391 https://doi.org/10.1021/acs.jcim.3c01033
  * R Skånberg, C König, P Norman, M Linares, D Jönsson, I Hotz, A Ynnerman, VIA-MD: Visual Interactive Analysis of Molecular Dynamics, 2018, Eurographics Proceedings, p. 19–27

* Specific tool:
  * Selection tool: Robin Skånberg, Mathieu Linares, Martin Falk, Ingrid Hotz, Anders Ynnerman, MolFind-Integrated Multi-Selection Schemes for Complex Molecular Structures, 2019, The Eurographics Association, p. 17-21​
  * Shape Space and Spatial Distribution Function: Robin Skånberg, Martin Falk, Mathieu Linares, Anders Ynnerman, Ingrid Hotz, Tracking Internal Frames of Reference for Consistent Molecular Distribution Functions, 2021, IEEE Transactions on Visualization and Computer Graphics, 28 (9), 3126-3137​

## Financial Support
<p align="center">
<img src="https://raw.githubusercontent.com/scanberg/viamd/master/docs/images/serc.png"  height="75" >
 &nbsp;&nbsp;&nbsp;
<img src="https://raw.githubusercontent.com/scanberg/viamd/master/docs/images/wallenberg.jpg"  height="75" >
 &nbsp;&nbsp;&nbsp;
<img src="https://raw.githubusercontent.com/scanberg/viamd/master/docs/images/infravis.png"  height="75" >
 
</p>

* VIAMD has received financial support since 2018 from the Swedish e-Research center ([SeRC](https://e-science.se/)) and the [Wallenberg Foundation](https://www.wallenberg.org/en)
* VIAMD is currently supported as a development project in [InfraVis](https://infravis.se/)

## Acknowledgements

https://github.com/glfw/glfw

https://github.com/dougbinks/enkiTS

https://github.com/ocornut/imgui

https://github.com/epezent/implot

https://github.com/BalazsJako/ImGuiColorTextEdit

https://github.com/skaslev/gl3w

https://github.com/max0x7ba/atomic_queue

https://github.com/mlabbe/nativefiledialog

https://github.com/nothings/stb

#
<p align="center">
<img src="https://github.com/user-attachments/assets/39b69b10-88a1-43a7-9d69-68513ac4e632"  width="150" alt="This is the VIAMD logo" >
</p>
