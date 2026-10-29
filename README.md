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
(persistence), into at most six, coloured in the plot. Zooming refines them: a
zoomed plot draws the curve of the CPs in view at the view's resolution (the
overview's curve stays as a faint line), and its peaks split at that curve's
valleys, again into at most six in view, numbered after the overview's peak
they lie in (2.1, 2.2) in shades of its colour. The viewport's Colour menu can
colour the CPs by their peak.

The environments table is organised by chemistry, not by peaks: one row per
Type, Kind, Chemical environment and, for BCPs, Character, in that order. A row
never mixes characters, and the kind and environment are written once for the
character rows that follow. The Detail menu sets how finely environments are
named: the group ("C aromatic"), then what the atoms are bonded to ("C aromatic
[C2H]", "H (on N [C2H])"; a contact ring's size and composition), then their
residues where the system has them ("in PHE", "between GLY and SER"). Auto, the
default, picks for each category and kind the finest level that gives at most
six environments. The Distribution column draws each row's CPs as a small
histogram along its plot's axis, following the plot's zoom, coloured by
character (other CPs in their category's colour). Grey arrows mark CPs outside
the view, so the peaks show up within each environment without splitting the
rows. Share is the row's fraction of its category's listed CPs. Sorting applies
within categories. The details pane always shows the finest name. The CSV
export includes each row's character and share.

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

The Plot selector above the filters chooses what the plots show, their peaks
group by and new filters select by: density, the energy density H = G + V, the
ratio |V|/G, an empirical hydrogen-bond energy E(HB) (Emamian et al. 2019 by
default, Espinosa et al. 1998 selectable; only for contacts classified as
hydrogen bonds), the Mayer bond order, or an approximate interaction energy
E(int) = qA qB / R - BO / (2R) from Mulliken charges and the bond order
(kcal/mol, negative favourable; a crude stand-in for an IQA pair energy). The
energy descriptors come from the analytic search only. G is computed at each CP
from the density matrix, and V follows from the local virial theorem. The bond
orders and charges need the overlap and density matrices. Filters are kept per
variable, and a CP must pass each enabled one. CPs without a value of the
plotted variable are left out of the Distribution column, whose tooltip counts
them and gives the reason. Each BCP also gets a Character from a non-covalent
interaction (NCI) analysis of the density, the same for every interaction. The
region around the BCP, on the surface between its two atoms, where the reduced
density gradient s ≤ 0.5 and ρ ≤ 0.05 au, is split by the sign of λ2. Negative
λ2 means attractive (density accumulates), positive means repulsive (density
depleted, steric crowding). Below 0.01 au, λ2 is near 0 and its sign means
little (NCI's green band). A region mostly below that density is Weak: mostly
dispersion, which the density does not capture, so it only gets a lean from the
density-weighted repulsive share of that weak part: an attractive lean
(< 35 %), a repulsive lean (≥ 65 %), or Weak, no lean, in between. Otherwise the
density-weighted repulsive share of the denser part decides: Attractive
(< 25 %), Mixed (25–50 %) or Repulsive (≥ 50 %). A BCP with ρ > 0.05 au is
beyond NCI's range: Covalent between bonded atoms, Attractive (strong) between
atoms that are not, such as a short, strong hydrogen bond.
One with no region is Undetermined. The thresholds are conventions, and "NCI repulsive share" can be plotted and filtered. The
BCPs line above the filters hides characters, and the viewport's Colour menu
can colour BCPs by character.

The analytic mdlib evaluator also provides the local Pauli kinetic energy
density t_P = G - |grad rho|^2 / (8 rho) at arbitrary points, with G, rho and
its gradient computed in one AO pass. At a critical point t_P equals G, so
there is no duplicate CP column. This is a total-density local descriptor,
not a fragment interaction energy or a repulsion classifier. It is not
used by the NCI character rules or exposed as a surface. The analytic contact
job samples it on the accepted, projected NCI region and retains density-weighted
means of t_P and of t_P / tau_TF, where tau_TF = 3/10 (3 pi^2)^(2/3) rho^(5/3).
These are local means over a projected disk, not a 3D integral or an interaction
energy. All region samples must be valid; incomplete evaluation leaves the
means unavailable with a reason. Dense contacts outside NCI's sampling range
have no regional measurement. Nonpositive density
(including screened vacuum) gives an undefined value (NaN); small negative
roundoff values are not clamped.
Neither mean has a universal clash threshold. In particular, normalization by
tau_TF can amplify dilute regions: at fixed orbital shape, scaling the density
by c scales t_P by c but t_P / tau_TF by c^(-2/3). Hydrogen-bond and crowded-contact
values can overlap; a large normalized value alone does not establish repulsion.
The Distribution variable selector exposes **Regional Pauli (experimental)**
and **Regional Pauli / TF (experimental)**, with the usual per-variable filters
and distribution/peak colouring. Both tables offer optional columns, hidden
by default (right-click a table header to enable them), marked with an asterisk
for experimental descriptors. Details show the regional values and sample
count, or the reason they are unavailable. CSV exports include both values;
the CP export also records valid samples and status. Environment and parent
means give each CP with a valid regional measurement equal weight, rather than
pooling the samples of regions of different sizes; their tooltips and exports
give the valid CP count.
Pauli remains automatic with the analytic contact analysis. Projection batches
reuse unchanged point values and omit points already excluded by the projection;
the kinetic evaluator shares the blocked density-matrix product used for density
derivatives. The energy/contact job uses atomic cancellation through mdlib's
cancellation callback, with additional checks between batches and during regional
processing. Cancelled or stale job results are not published.
Mulliken populations also supply the electron-count trace without another
density/overlap matrix pass; population and bond-order work honour the same
job cancellation.

The NCI surface checkbox draws the NCI isosurface (Johnson et al., J. Am. Chem.
Soc. 2010) in the viewport. It is the surface s = 0.5 where ρ is below the
cutoff (0.05 au), coloured by sign(λ2)ρ: negative attractive, near-zero weak,
positive repulsive. The default Jet palette maps these to blue, green
and red, respectively. The density is sampled on a grid (Coarse 0.3, Medium 0.2 or Fine
0.13 Bohr) and differentiated by finite differences, so a finer grid gives a
smoother and more accurate surface. Scope offers Everywhere, Around shown
(the visible BCPs, respecting filters), and Around selected (the default).
Around selected uses selected CPs of any type, regardless of their visibility
or filters; selecting a ring CP therefore shows its ring's surface. With
nothing selected, no surface is drawn. The default radius is 2.7 Bohr. The radius, isovalue, cutoff, colour
range and opacity can be adjusted. Progress and Cancel stay in the Critical
Points window, without an NCI task overlay or routine completion notifications
in the viewport. Timings remain in debug logs; errors are still reported. It is drawn as a
mesh in both isosurface rendering modes. Everywhere does not need the
critical point search; the other scopes use its results.
The Colour scale section uses the same palette selector, automatic/manual
range, symmetric range and optional viewport legend as the other surface
fields. The default palette is Jet, with an automatic, zero-centred symmetric
range and a visible legend. Auto uses the 1st to 99th percentile of values on
the current scoped surface: selected regions, shown regions or Everywhere.
It updates when a new scoped surface arrives, not with camera motion; outliers
only saturate in colour. Manual keeps a fixed range for comparisons between
selections. Changing the colour scale does not recompute the surface.

The distributions pane also offers **NCI scatter**: reduced density gradient
s against sign(λ2)ρ, using the same scope, density cutoff and colour scale as
the NCI surface. It samples the density volume, not the surface vertices.
Up to 50,000 grid samples with s ≤ 2 are retained by deterministic reservoir
sampling, with the displayed and eligible counts shown above the plot.
The horizontal reference line marks the surface isovalue; changing that
isovalue leaves the scatter samples unchanged. This plot does not change the
character thresholds or apply table filters (except through Around shown).
With the CP search disabled, Everywhere still provides the plot below the
NCI controls.

The environments table is a collapsible chemistry hierarchy of the CPs of the
shown categories that pass their filters. Type is the first column and
Chemical environment the second; the tree starts at chemistry branches,
without an extra Attractor/BCP/RCP/CCP dropdown. BCPs split into Bonds, Contacts and
Bond paths with unknown connectivity. Bonds group by element pair (C-C, C-N,
etc.); Contacts by interaction family (Hydrogen bonds, pi-stacking, etc.)
and pair, with donor-H...acceptor groups for hydrogen bonds. Ring CPs group
into aromatic, aliphatic, sugar or contact-closed rings when recognized.
Detail still controls the environment names at the leaves. Every leaf repeats
its environment name, even when separate characters share a name.
Character remains a column, not a branch: neither density peaks nor the NCI
character determine the chemistry hierarchy.
Kind is omitted from this table because the branches already identify it;
it remains in the individual CP list and CSV exports.

Parent rows summarize all listed descendants, including collapsed ones:
count, share of the CP category, a character-stacked distribution, density
range and means of Laplacian, ellipticity and energy descriptors. Means use
all descendants with a value, not averages of the child means. A parent
with several characters says Multiple, not the NCI classification Mixed.
Arrows expand/collapse without selecting. Collapse all, Expand all and
Reveal selection control the tree.

The critical point list below it shows the CPs of the environments chosen in
the table above, else of the bars clicked, else every CP the environments table
counts (List all returns to that), with the same first columns (Type, Kind,
Chemical environment, index) and then density, Laplacian, ellipticity and
position. Both tables sort by clicking a header, Shift-click further headers to
sort by several columns; the environments tree sorts siblings without
separating their descendants. Right-click a header to hide columns (Type and
Chemical environment always stay).

The Export menu at the top writes two kinds of data. Topology (JSON or VTP) is
every critical point and path, unfiltered and without chemical labels, for scripts
and ParaView. Chemical analysis (CSV) is either table as listed and sorted, every
column at full precision (hidden ones too); the environments export all
filtered leaves in tree order, including collapsed ones, with a hierarchy path
and their CP indices (parent summaries are omitted to avoid double counting),
the CP list each CP's label detail, and the CP index joins it to the topology
export.

Selection works at both levels. Click an environment or a parent category to
select all its listed CPs, even in collapsed descendants. Ctrl-click or
Shift-click adds the whole group, or removes it if it is fully selected.
Hovering a parent highlights all its listed descendants. Its environments
are listed below and its CPs become the selection. Click a
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
