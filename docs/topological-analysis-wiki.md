# Topological Analysis

> Requires the development `CP` branch. This page describes the implementation at
> VIAMD revision `8a9337cf` plus the chemistry/reference corrections of
> 4 November 2026, not a released upstream version.
> Standalone wiki-page draft. Wiki assembly and publication take place after merge.

Use **Windows > Quantum Chemistry > Topological Analysis** to explore the
topology of an electronic density, chemical environments of its critical points
(CPs), and non-covalent interaction (NCI) regions. The window is available when
the loaded quantum-chemistry system supplies atoms and molecular orbitals.
For loading VeloxChem HDF5 output, see
[VIAMD for VeloxChem](https://github.com/mathieulinares/viamd/wiki/VeloxChem).

## What this window does

The window has two independently resizable columns:

- **CP:** search, four CP-category distributions, a chemical-environment tree
  when enabled, and an individual-point table.
- **NCI:** an RDG scatter plot, reference maps and guides, the full-system layer,
  a selection preview, and saved session layers.

CP search, chemical analysis, and NCI computation are separate switches.
All three are off initially. **Show CPs** only changes marker/path/label
visibility; it does not disable computation or erase a selection.
The former fixed Details pane has been removed. Hover point and environment
rows for compact information; enable optional columns or export for more detail.

## Start with a small system

1. Load a VeloxChem wavefunction-bearing HDF5 file using **File > Open File...**.
2. Open **Windows > Quantum Chemistry > Topological Analysis**.
3. Check **Enabled** in the CP column.
4. Leave the default **Analytic** search selected while its progress completes.
5. Check **Chemical environments** to enable chemical assignments and descriptors.
6. Select a chemical-environment row or individual CP.
7. Check **NCI enabled** in the NCI column.
8. Inspect the selected CPs' **Preview** layer and the scatter plot.
9. Wait for the preview computation to finish before using **+** above the layers.
10. Change one parameter at a time and inspect both the plot and real-space region.

Large systems take longer. Turning off chemical analysis keeps CP search and
NCI available, but removes chemistry-dependent information. In this mode CP
anchors may be grouped as Unassigned rather than by chemical interaction.

## Critical points: what the categories mean

A CP is a stationary point of the electron density. The three eigenvalues of
its Hessian distinguish the four nondegenerate categories:

| Category | Signature | Meaning in the density topology |
|---|---|---|
| Attractor | (3, -3) | Local density maximum; usually associated with a nucleus |
| Bond CP / BCP | (3, -1) | Saddle with two negative and one positive eigenvalues |
| Ring CP / RCP | (3, +1) | Saddle with one negative and two positive eigenvalues |
| Cage CP / CCP | (3, +3) | Local density minimum |

These are mathematical categories, not a list of chemical bonds or interaction
energies. In particular, finding a bond path does not by itself establish an
energetically stabilizing chemical bond.

The four category toggles and counts are above the density plots. They affect
which CPs are shown and listed, not which CPs were found by the search.

### CP settings

The single **CP settings** button contains:

- **Search:** Analytic or Volume method, the lower density threshold, and the
  available backend/resolution options.
- **Appearance:** marker size, paths, colours, labels, and how the viewport
  displays a selection.
- **Distributions / filters:** plotted descriptor, its options, per-variable
  filters, and BCP character filters.
- **Export:** topology JSON/VTP and environment/point CSV files.

**Enabled**, **Show CPs**, and **Chemical environments** remain outside this
popup for quick access.

The Analytic method searches the basis-function density and polishes CP
positions; the Volume method operates on sampled grid data. Precision and
available descriptors differ. The search threshold excludes low-density
regions: lowering it searches further and can take longer. Search status and
warnings matter; a finished job should not be read as an unconditional proof
that every possible CP in an arbitrary wavefunction was found.

## Chemical environments and character

Chemical assignment uses nuclei, perceived connectivity, traced CP paths and
geometry, with charge and ring information where available. It is not assigned
from the density-distribution peak containing the CP.

The environment tree distinguishes atoms, covalent bonds, nonbonded contacts,
bond paths with unknown connectivity, rings, cages and unassigned cases.
Contacts are grouped into chemical families and finer element-pair or
donor-H...acceptor environments. Ring groups distinguish aromatic, aliphatic,
sugar and contact-closed environments where recognized.

The catalogue includes hydrogen and ionic hydrogen bonds, halogen/chalcogen/
pnictogen/tetrel/triel/aerogen contacts, pi-related contacts, dihydrogen contacts,
coordination and other weak contacts. Recognition is rule-based and conservative;
an available label is not a claim that every member is stabilizing.
Ambiguous assignments retain an explicit unassigned explanation.

**Character** is a separate descriptor, not a branch defining chemical identity.
Its categories include attractive, mixed, repulsive and weak variants, plus
covalent, strong attractive and undetermined cases. These use implementation
thresholds and a separate local analysis, not a fragment energy decomposition.
The character analysis is not redefined by changing the display surface or
scatter extraction parameters.

With **Chemical environments** off, the Chemical environment and Character
columns are unavailable, chemistry-dependent filters are ignored, and CP
tooltips omit chemical character. Density filtering remains usable.

### Distributions, filtering and selection

Distributions describe the population of CP values, not the spatial electron
density. Density distributions normally use a logarithmic density coordinate.
The curves are adaptive smoothed estimates; bars represent actual CPs.
Nearby bars are grouped on screen with their count.

- Use the wheel to zoom, drag to pan, and double-click to fit.
- Choose the plotted variable under **CP settings > Distributions / filters**.
- Add filter ranges and drag their handles, or enter numerical bounds.
- Multiple ranges for a variable form a union. Enabled filters on different
  variables must all pass.
- Filters change the visible/listed subset, not the stored topological graph.
- Click an environment parent to select its listed descendants, including
  collapsed descendants. Use <kbd>Ctrl</kbd> or <kbd>Shift</kbd> to extend selection.
- Click an individual point to select it; <kbd>Shift</kbd> toggles membership.
- Hovering highlights the corresponding CPs without changing selection.
- Use **List all** and **Clear selection** to recover from a narrow list scope.

The tree summarizes descendants rather than averaging averages of child rows.
Sorting is hierarchical. Right-click a table header to expose optional columns;
**Compact columns** restores the concise column layout.

Available descriptors include density, Laplacian, ellipticity, energy-density
quantities, hydrogen-bond energy estimates, bond orders and approximate pair
energies when their input matrices are available. Missing values are not zero.
Empirical estimates have restricted domains of applicability.
Regional Pauli descriptors are explicitly experimental; they are neither
fragment interaction energies nor a validated clash score.

## NCI scatter and surfaces

The reduced density gradient is

$$s = \frac{|\nabla\rho|}{2(3\pi^2)^{1/3}\rho^{4/3}}.$$

The plot shows **s** against **sign(lambda2) rho**, with lambda2 the middle
eigenvalue of the density Hessian. Low-density, low-gradient features can reveal
non-covalent contact regions. Their real-space location is essential to
interpretation.

- Negative signed density indicates the lambda2 < 0 side, associated with
  density concentration.
- Positive signed density indicates lambda2 > 0, often associated with
  depletion/steric features.
- Near-zero signed density is the low-density region often shown in green.

Do not interpret the horizontal coordinate or a blue/red patch as an
interaction energy. Geometry, density method, spatial scope and nearby regions
all influence the plot.

### Defaults and independent controls

| Control | Default | What changing it does |
|---|---:|---|
| Scatter extraction s max | 2.0 | Changes analysed samples, contours and score |
| Extraction rho max | 0.07 au | Rejects denser samples; applies to surface and scatter |
| CP-local radius R | 2.7 Bohr | Changes the spatial scope around CP anchors |
| Surface s isovalue | 0.5 | Changes the rendered surface, not the analysed scatter |
| Plot s range | 2.0 | Changes the visible/display-context range |
| Sampling resolution | Medium | Changes the sampled grid, accuracy and cost |
| Interaction line width | 1 px | Display only |
| Scatter point size | 0.6 px | Display only |
| Layer contour width | 2 px | Display only |

Use **NCI settings > Scatter / extraction > Scatter region** for extraction
settings. Resetting extraction defaults does not redefine the reference map.
The surface cannot exceed the extraction s cutoff; lowering the cutoff below
the surface also lowers the surface isovalue.

Regions are connected, CP-anchored and radius-limited, with nearest-CP ownership.
These are VIAMD spatial filters, not requirements of conventional NCI.
There is no CP-adaptive mode or per-interaction radius.

At s max = 2 there is no above-cutoff tail within the default plot domain.
If you lower the extraction cutoff:

1. Samples below the surface isovalue use strong saturated colour.
2. Samples between the surface and extraction cutoff use pale desaturated
   colour; they are still analysed and pickable.
3. Samples above extraction cutoff, up to the plot range, use faint neutral
   grey; they are display-only and excluded from contours and scores.

Raising extraction s max includes more points in analysis and leaves fewer grey
points. Increasing plot range alone does not promote grey points into analysis.

### Plot guides

**Guides** beside the plot is a master visibility switch. Under
**NCI settings > Display > Plot guides**, choose:

- **Surface isovalue:** a solid draggable line labelled **Surface s = ...**.
- **Extraction cutoff:** a dashed draggable line labelled **Cutoff s = ...**.
  Dragging changes the analysed region; expensive recomputation waits until
  the drag is released.
- **Extraction density limits:** vertical lines at exactly +/- the configured
  rho max, with labels inside the plot. These are not automatic colour limits.

When surface and cutoff coincide, one combined line is drawn; dragging it
edits the cutoff. The surface can still be set separately in Settings.
The guide choices survive use of the master switch.
Both axes use VIAMD's inset ticks and titles. The fit button or double-click
restores the selected plot range.

The colour scale is separate: Auto follows the shown surface, and a manual
range supports visual comparisons. Changing colour mapping does not change
the extraction limits.

## Layers

- **Reference:** fixed reference data, independent of the loaded molecule.
- **Full system:** the CP-anchored regions of the system, not an unrestricted
  conventional whole-box NCI calculation.
- **Preview:** the selected CPs, updated as selection changes.
- **Saved layers:** snapshots of the preview selection, retained for the
  current loaded geometry. They are not written to a persistent layer file.

Click **+** after the preview has finished to save it. Saved layers use the
short shared chemical name where possible, such as N-H...O or Br...O.
Mixed selections use a shared branch or combined name.

Only one layer is shown in the viewport at a time. Other layers can remain as
background contours or heatmaps in the scatter.
The eye controls scatter visibility; selecting a card controls the shown
viewport layer. The inspector configures each layer's contours/heatmap, colour,
opacity, fills, line width and enclosed-volume levels.

Click the **trash can** at the left of a saved card to remove it without a
confirmation dialog. Alternatively, right-click the card and choose
**Remove**, or select it and press <kbd>Delete</kbd> while the panel has focus.
Reference, Full system and Preview are built-in and cannot be deleted.

## Reference maps, guides and percentages

The active reference spans **s = 0 to 2** and **signed rho = -0.1 to +0.1 au**.
It was regenerated from **910 calculations**: 48 VeloxTest systems,
528 S66x8 geometries, 205 usable SH250 systems, 100 IHB100 systems and
29 hand-built systems. S66x8 separations share one dimer's statistical weight.
Effective-core-potential systems unsupported by the density/topology pipeline
were excluded, not silently treated as all-electron systems.

The four pooled categories are Hydrogen bonds, Rings, pi interactions and Other.
These chemical categories can overlap substantially in scatter coordinates.
The map is not a classifier that uniquely identifies an unknown interaction.

The Reference inspector configures combined/dominant colouring, opacity,
fill/lines and highest-density contours enclosing 50%, 90% and 98% of each
category's reference volume. The same extended dataset supplies maps, chemical
guides and comparison scores; there is no original/extended selector.

Expand **Interaction guides**, then **Add interaction guide...** to browse
by category or search by chemical name. A broad class offers **Whole class**
and finer supported maps. The reference contains **18 broad maps and 89
refinements**, including donor/acceptor, intra/intermolecular, ring and
element-pair subdivisions. Guides require at least five supporting pooled
systems. Recognized classes without sufficient data are listed separately,
not presented as supported reference contours.

A selected guide draws its 90% isoline. Its visibility is independent of the
category reference fill and of the master plot-guide checkbox.
Hover map entries for descriptions and supporting-system counts.

### What “in reference (90%)” means

For each retained analysed scatter sample, VIAMD checks whether it lies inside
its CP category's fixed extended-reference 90% contour. The displayed percentage
is the fraction of these samples inside that contour. Display-only grey samples
do not enter the count; the surface isovalue does not redefine the sample set.
Changing extraction s max, rho max or radius can change the percentage.

The plotted reference density and contours are volume-based; the reported
percentage is counted from the capped retained sample cloud, not computed by
integrating the entire live continuous density grid. Do not equate it with an
exact full-grid volume fraction.

**90% describes the reference contour, not an expected score for every system
or restricted subset.** The score is not a probability that the interaction
exists, an interaction strength, or a validated structure-quality metric.
The reference was built with R = 2.7 Bohr; changing R changes comparability.
Beyond |rho| = 0.1, the reference has no data. A tighter extraction domain is
an intentional subset, not missing reference coverage.

## Export and reproducibility

Use **CP settings > Export** for:

- **Topology JSON / VTP...:** the unfiltered CP graph and paths. Positions are
  in Bohr; use VTP in ParaView.
- **Environments CSV...:** filtered leaves in hierarchy order, including
  collapsed leaves, with CP indices. Parent rows are omitted to avoid
  double-counting.
- **Critical points CSV...:** filtered/sorted CP rows with hidden descriptor
  columns at full precision and chemical assignment reasons.

**Export NCI samples...** writes the shown layer and category-specific samples,
CP-anchor membership and continuous density data in companion files.
The in-application category export is a reference-oriented recomputation;
record the actual settings written in the files rather than assuming it is
identical to a modified live scatter region.
For reference generation, the headless `viamd_nci_reference` utility supports
explicit extraction settings and per-CP density exports.

For scientific use, record the exact VIAMD/mdlib revisions, VeloxChem method,
basis and input geometry, density source, CP search threshold, grid spacing,
extraction domain, surface value and reference version. Recheck grid and cutoff
sensitivity before interpreting changes between systems.

## Theory, implementation and limitations

Conventional NCI uses the electron density and RDG field without requiring
critical-point anchors. CP-local connected extraction, chemical ownership,
the pooled category reference and its coverage percentages are VIAMD additions.
The extraction rho cutoff 0.07 and radius 2.7 Bohr are adjustable scope choices,
not universal theoretical constants. The SCF scatter s max 2 and surface
isovalue 0.5 are practical starting values, not a physical interaction boundary.

The datasets provide chemical coverage, but do not validate universal
classification or scoring accuracy. Unassigned chemistry, finite grid
resolution, basis/density-method dependence, nearby covalent density, unequal
dataset representation and unsupported ECP densities all limit interpretation.
Experimental energy/Pauli descriptors must be distinguished from the
geometrical/density information.

See the [implementation and reference-generation document](topological-analysis-implementation.md)
for algorithms, all datasets and rationales, the full map atlas, citations,
parameter choices and reproduction instructions.

Primary starting points:

- Johnson et al., *Revealing Noncovalent Interactions* (2010),
  [doi:10.1021/ja100936w](https://doi.org/10.1021/ja100936w).
- Contreras-Garcia et al., *NCIPLOT: A Program for Plotting Noncovalent
  Interaction Regions* (2011),
  [doi:10.1021/ct100641a](https://doi.org/10.1021/ct100641a).
- [Official NCIPLOT manual](https://github.com/aoterodelaroza/nciplot/blob/master/doc/nciplot-manual.pdf).

Related published pages:
[Density volume](https://github.com/mathieulinares/viamd/wiki/Density-volume),
[Selection](https://github.com/mathieulinares/viamd/wiki/Selection),
[Distributions](https://github.com/mathieulinares/viamd/wiki/Distributions),
[VIAMD for VeloxChem](https://github.com/mathieulinares/viamd/wiki/VeloxChem).
