# Interactive Quantum Chemical Analysis: VIAMD and VeloxChem 1.0 Integration

## Author List and Contributions
The author order should be:
1. Robin Skånberg
2. Gustav Eriksson
3. Tobias Pettersson
4. Mohit Sharma
5. Iulia Brumboiu
6. Josefine Hvarregaard Andersen
7. Xin Li
8. Talha Bin Masood
9. Patrick Norman
10. Mathieu Linares

Affiliations currently used in the manuscript:
* Robin Skånberg: PDC Center for High Performance Computing, KTH Royal Institute of Technology, SE-100 44 Stockholm, Sweden.
* Gustav Eriksson: Linköping University, SE-581 83 Linköping, Sweden.
* Tobias Pettersson: Linköping University, SE-581 83 Linköping, Sweden.
* Mohit Sharma: Linköping University, SE-581 83 Linköping, Sweden.
* Iulia Brumboiu: Nicolaus Copernicus University (UMK), Torun, Poland.
* Josefine Hvarregaard Andersen: UiT The Arctic University of Norway.
* Xin Li: PDC Center for High Performance Computing, KTH Royal Institute of Technology, SE-100 44 Stockholm, Sweden.
* Talha Bin Masood: Linköping University, SE-581 83 Linköping, Sweden.
* Patrick Norman: PDC Center for High Performance Computing, KTH Royal Institute of Technology, SE-100 44 Stockholm, Sweden; Division of Theoretical Chemistry and Biology, School of Engineering Sciences in Chemistry, Biotechnology and Health, KTH Royal Institute of Technology, SE-100 44 Stockholm, Sweden.
* Mathieu Linares: PDC Center for High Performance Computing, KTH Royal Institute of Technology, SE-100 44 Stockholm, Sweden.

The Author Contributions section should state:
* Robin Skånberg was the main developer of the code and contributed to the code design.
* Gustav Eriksson implemented the electronic-structure graphical interface, including interactive spectra, vibrational analysis, and transition analysis.
* Tobias Pettersson contributed to the transition-analysis functionality.
* Mohit Sharma provided the first implementation and guidelines for the critical point analysis.
* Iulia Brumboiu, Josefine Hvarregaard Andersen, and Xin Li contributed to the VeloxChem interface and provided user feedback on the graphical interface.
* Talha Bin Masood advised on electronic-structure computation, visualization, and topological analysis.
* Patrick Norman initiated the project, co-led the project, and contributed to the VeloxChem interface.
* Mathieu Linares led the project and contributed to the code design.

## I. Introduction
The introduction should motivate the article by first placing electronic-structure visualization in the broader context of quantum chemistry workflows. It should explain that orbitals, densities, spectra, excited states, and related molecular properties are often analyzed through separate programs and intermediate files.

It should then identify the limitations of this workflow: large grid files such as cube files, fragmented analysis across multiple tools, and the lack of a direct interactive connection between molecular geometry and quantum chemical data.

The introduction should present **VIAMD** as a C/C++ environment originally developed for visual interactive analysis of molecular dynamics and spatio-temporal molecular data, and then position it as a platform for interactive quantum chemical analysis.

The introduction should introduce **VeloxChem 1.0** as the target electronic-structure program for this integration and emphasize native support for VeloxChem **HDF5 (.h5)** output files in VIAMD.

Finally, it should state the scope of the article: loading VeloxChem files, inspecting calculation metadata and convergence, visualizing orbitals and densities, analyzing electronic and vibrational spectra, interpreting excited states through transition-analysis tools, and using VeloxChem input/output documentation datasets as reproducible examples.

## Abstract
The abstract should summarize the problem, contribution, implementation, and first performance results in one compact paragraph. It should state that quantum-chemical analysis is often fragmented across electronic-structure output files, cube files, molecular viewers, and plotting tools. It should then present the VIAMD--VeloxChem 1.0 integration as a direct HDF5-based route for loading and inspecting molecular structures, metadata, orbitals, densities, spectra, vibrational modes, transition-analysis data, and critical points. The abstract should mention on-demand GTO volume generation and the shared parsed VeloxChem data model, but avoid excessive implementation detail. It should close with the first CPU benchmark result at a qualitative level: millisecond-scale HDF5 parsing and interactive molecular-orbital volume generation for representative datasets.

## II. Initial Audit: Loading Data and the Summary Window
Files can be loaded via the `File > Load Data` menu or by simply **dragging and dropping** the `.h5` file directly into the spatial view. Upon loading a VeloxChem dataset, VIAMD automatically initializes a **Ball and Stick** representation for the structure and, when molecular orbital data are available, a visualization of the **HOMO orbital**.

The **Summary Window** (`Windows > VeloxChem > Summary`) should be presented as the first audit point after loading a dataset. It includes:
*   The **level of calculation** (method or functional and basis set).
*   System metadata such as number of atoms, charge, spin multiplicity, and electron count.
*   **SCF convergence** plots to verify electronic convergence.
*   **Geometry convergence**, when present in the HDF5 file, for optimization-related calculations such as geometry optimizations, scan calculations, transition-state searches, and IRC paths.
*   An **interactive geometry table** where selecting an entry highlights the corresponding atom in the 3D view.
*   **Interactive critical point analysis** computed by VIAMD from the electron density reconstructed from the molecular orbitals; the critical point table is also linked to the 3D molecular view.

The figure for this section should use `figure/summary.png` and emphasize metadata, convergence diagnostics, the interactive geometry table, and the interactive critical point table.

Installation details should not be treated as a main article section. The availability section should briefly mention that precompiled binaries are provided for Windows and macOS on both Intel and Apple Silicon architectures, while Linux build instructions are available on the VIAMD wiki.

## III. Visualizing the Electronic Landscape
This section should present the electronic-structure visualization workflow as two complementary tools: electronic representations in the main spatial view and the dedicated **Orbital Grid** window.

In the main view, users can add electronic volumetric representations and choose the volume source from the data available in the VeloxChem file and derived by VIAMD. Relevant sources include:
*   Molecular orbitals.
*   Orbital densities.
*   Total electron density.
*   Differences between compatible densities, such as alpha--beta density differences.
*   Natural Transition Orbitals (NTOs).
*   Attachment and detachment densities for TDDFT calculations.
*   Attachment--detachment density differences for TDDFT states, including charge-transfer states.
*   Differences between selected NTO-related quantities.
*   Density-related properties from complex polarization propagator (CPP) calculations.

The text should explain that these representations can be inspected together with the molecular geometry and adjusted through visual parameters such as isovalue, opacity, colors, and phase convention.

This point should connect to the spectroscopy and transition-analysis sections: a TDDFT charge-transfer state can later be illustrated with an attachment--detachment density-difference representation.

The **Orbital Grid** window (`Windows > VeloxChem > Orbital Grid`) should be described as a compact comparison view for multiple orbitals or density-derived quantities. Users can compare orbital shapes, phases, energies, occupations, and localization patterns side by side. The text should also mention support for ROHF and UHF calculations, with UHF orbitals split into alpha and beta channels for separate inspection.

The representation figure for this section should be expanded to four panels. This is a stronger visual choice than only showing a generic representation view because it demonstrates the range of representation types supported by the VeloxChem integration. The four panels should be:
*   `figure/representation-homo.png`: the HOMO representation generated upon loading a VeloxChem dataset.
*   `figure/representation-density-dipole-pna.png`: a green total-electron-density isosurface shown together with a molecular dipole or PNA representation.
*   `figure/representation-density-resp.png`: an electron-density representation shown with RESP-based atom coloring.
*   `figure/representation-alpha-beta.png`: an alpha--beta density-difference representation for a UHF calculation.

The manuscript text should emphasize the conceptual progression across the four panels: automatic orbital inspection, density plus directional descriptor, density plus atom-centered electrostatic coloring, and spin-density difference. This keeps the figure focused on representation diversity rather than on file generation.

The Orbital Grid should still be described in the text and should remain as a large standalone figure using `figure/orbital-grid.png`, because it represents a different inspection mode: compact comparison of multiple orbitals or density-derived quantities rather than a single spatial representation. It should not be compressed into the four-panel representation figure, since the grid contains many small orbital panels and interface controls that need space to remain legible.

## IV. Spectroscopy: Electronic and Vibrational Responses
This section is explicitly split into electronic and vibrational spectroscopy. It focuses on representative workflows rather than covering every VeloxChem example, while still noting that other spectroscopies are available.

The opening paragraph should explain that spectroscopy continues the electronic-structure visualization workflow by linking spectral data to molecular and electronic representations. The article should focus on **TDDFT electronic spectroscopy** and **IR vibrational spectroscopy**, while mentioning that other workflows are available, such as optical activity, CPP-based response calculations, Raman-type vibrational spectroscopy, and X-ray spectroscopies.

### Electronic spectroscopy from TDDFT
The **Response Window** (`Windows > VeloxChem > Response`) should be introduced with TDDFT spectra as the main electronic-spectroscopy example. Users can apply **Lorentzian or Gaussian broadening**, adjust the broadening width, and export spectral data to **.xvg or .csv** formats.

The TDDFT figure should use `figure/tddft.png` as a single composite figure. It should show the interactive TDDFT spectrum on the left and, for the fourth excited state of the TQ system, the difference between attachment and detachment densities on the right. The purpose is to showcase a charge-transfer state and prepare the transition-analysis discussion that follows.

### Vibrational spectroscopy from IR calculations
VIAMD should present IR as the representative vibrational-spectroscopy example, also accessed through the **Response Window**. The section should describe plotting vibrational frequencies/intensities and **animating specific vibrational modes** directly within the spatial view. The text should emphasize the link between a spectral peak and the corresponding nuclear displacement pattern.

The IR workflow will remain text-only for now. Do not include `figure/ir-spectrum.png` unless a clear static representation of vibrational motion is produced later, such as a ghosted multi-frame overlay, displacement arrows, or paired positive/negative normal-mode geometries.

## V. Advanced Transition Analysis: Charge Transfer and Transition Moments
This section should describe the transition-analysis workflow as following the **VALET** approach for visual analysis of electronic densities and transitions in molecules, citing `MasoodCGF2021`.

The analysis should be described as based on **attachment and detachment densities**, not primarily on NTOs. NTOs can remain part of the broader visualization workflow, but the transition-analysis section should focus on:
*   A selected excited state from the **Response Window**.
*   **Detachment density**, showing where electron density is removed during excitation.
*   **Attachment density**, showing where electron density accumulates during excitation.
*   Attachment--detachment difference representations for highlighting charge redistribution.
*   The fourth excited state of the TQ system as a charge-transfer example connected to the TDDFT spectroscopy section.

The settings should be described as allowing users to change attachment/detachment colors and control transition-analysis display options. The section should also mention display of:
*   The electric transition dipole moment vector.
*   The magnetic transition dipole moment vector.
*   The angle between the electric and magnetic transition dipole moments.

For fragment analysis, users can define molecular subgroups/fragments and analyze how attachment and detachment densities are distributed over them. This supports quantification of **local excitation versus charge-transfer** character.

The figure for this section should use `figure/transition-analysis.png` and show the transition-analysis interface with attachment/detachment density visualization, transition moment vectors, and charge-transfer analysis controls.

## VI. Implementation
This section should be placed after the transition-analysis section so that the article first presents the user-facing workflow and then explains how it is implemented.

The implementation section should be grounded in the VIAMD and mdlib code. It should explain that VeloxChem HDF5 files are registered by the VIAMD loader when VeloxChem support is enabled, and that the generic loader initializes the molecular system while the VeloxChem component parses the full dataset into an `md_vlx_t` object.

The text should describe the mdlib parser as reading the core molecular data and the `scf`, `rsp`, and `opt` HDF5 groups when they are present. These groups provide SCF orbitals and densities, response data including excited-state and vibrational spectra, and optimization coordinates and energies. Do not describe `vib` as a separate HDF5 group; vibrational information belongs to the response workflow. Atomic and density properties stored in the file are also made available to VIAMD.

The section should explain the GTO reconstruction path: mdlib resolves the basis-set identifier, reads the basis file from packaged VIAMD resources or from the dataset directory, normalizes the basis, builds an AO remapping table, permutes SCF coefficient and matrix data into canonical shell order, and extracts an `md_gto_basis_t`. Molecular-orbital coefficients and SCF density matrices are then used directly by VIAMD for volume evaluation. Restricted calculations share alpha/beta storage, while unrestricted calculations keep separate alpha and beta data.

The text should describe on-demand volume generation. VIAMD hashes the electronic representation state—source, spin, orbital or excited-state index, transition-density component, NTO component, resolution, and frame-dependent state—and reevaluates only when the hash changes. Volumes are sampled on a molecular grid and stored in 3D OpenGL textures. Orbitals and NTOs are evaluated from AO coefficient vectors; electron densities, transition densities, and density properties are evaluated from AO-space matrices. The GPU path packs atoms, coefficients, and matrices into mdlib GPU buffers; otherwise the OpenGL GTO grid-evaluation path is used.

For transition analysis, the section should state that attachment/detachment densities are computed from VeloxChem TDDFT response solution vectors in mdlib. The implementation builds `T = Z - Y`, forms attachment and detachment matrices in virtual and occupied MO spaces, transforms them to AO matrices, symmetrizes them, and returns attachment, detachment, or difference matrices to VIAMD.

The section should close by connecting implementation to the user-facing workflow: Summary queries metadata, SCF history, optimization data, and geometry from `md_vlx_t`; critical-point analysis is computed in VIAMD from an evaluated SCF density volume using mdlib topology routines; Response uses excitation, oscillator/rotatory strength, transition-dipole, CPP, and vibrational arrays; Orbital Grid uses the same coefficient accessors and GTO volume evaluation as the main representation system.

## VII. Performance
This section should follow the implementation section and support the main claim that direct HDF5-based analysis enables interactive workflows.

The section should explain the performance argument without inventing benchmark numbers. It should first emphasize that VIAMD avoids the repeated generation, storage, and loading of intermediate cube files. In a conventional workflow, each orbital, density, excited-state descriptor, or density difference often requires a separate grid file at a fixed resolution. In the VIAMD--VeloxChem workflow, the HDF5 file remains the central data source and volumes are reconstructed on demand.

The section should distinguish between cheap interactive rendering changes and expensive volume recomputation. Changing colors, opacity, signs, or visibility should be described as a rendering update, while changing orbital index, spin channel, excited state, transition-density component, or grid resolution should be described as triggering a new GTO-based volume evaluation. The text should mention that VIAMD hashes the representation state to avoid unnecessary reevaluation.

The performance factors to discuss are:
*   Number of atoms.
*   AO basis size.
*   Grid resolution.
*   Representation type.
*   Hardware and availability of the GPU path.

The section should explain that molecular orbitals and NTOs evaluate coefficient vectors over the grid, whereas electron densities, density properties, attachment densities, and detachment densities evaluate AO-space matrices. Critical-point analysis should be described as more expensive because it combines density evaluation with topological extraction.

The text should state that when GPU support is available, VIAMD packs atom positions, coefficients, and density matrices into mdlib GPU buffers for GTO volume evaluation; otherwise the OpenGL GTO grid-evaluation path is used. Avoid numerical speedup claims until actual benchmarks are collected.

The first-draft benchmark protocol now uses a dedicated `viamd_perf_benchmark` target in `article/performance/viamd_perf_benchmark.c`. The target links against mdlib and measures CPU-side operations on existing VeloxChem test datasets. The current configuration is a Release build with VeloxChem/HDF5 support enabled, AVX/AVX2/FMA enabled, and the mdlib GPU path disabled. The first benchmark uses a $64^3$ grid and nine repeats, reporting median wall time.

Measured first-draft median timings:
*   H2O (`h2o.h5`, 3 atoms, 24 AOs, 1 state): HDF5 parse 1.3 ms, basis extraction 1.2 ms, molecular-orbital volume 3.9 ms, transition-density extraction 1.6 ms.
*   Molecule (`mol.h5`, 26 atoms, 229 AOs, 1 state): HDF5 parse 3.2 ms, basis extraction 3.0 ms, molecular-orbital volume 31.7 ms, transition-density extraction 30.3 ms.
*   TQ SCF (`tq.scf.results.h5`, 23 atoms, 254 AOs, no response states): HDF5 parse 3.7 ms, basis extraction 3.4 ms, molecular-orbital volume 35.9 ms.
*   Acrolein RSP (`acro-rsp.h5`, 8 atoms, 76 AOs, 5 states): HDF5 parse 2.4 ms, basis extraction 2.2 ms, molecular-orbital volume 20.9 ms, transition-density extraction 2.9 ms.

The section includes one generated figure:
*   `figure/performance-timings.png`: bar chart of the first-draft CPU benchmark timings for HDF5 parsing, basis extraction, molecular-orbital volume generation, and transition-density extraction.

## VIII. Workspace Management
All analysis settings—including camera position, representations, and custom subgroup definitions—can be saved in a **.via workspace file**. This allows for the immediate restoration of a complex quantum chemical analysis session.

This section should explain that the HDF5 file remains the source of the quantum-chemical data while the `.via` workspace stores how the data are inspected in VIAMD. The text should mention restoration of camera views, molecular and electronic representations, isovalues/colors, selected spectra or states, and custom molecular subgroups/fragments. The section should emphasize reproducibility, collaboration, and continuing complex analysis sessions without manually rebuilding the same visualization state.

## IX. Conclusions
The conclusion should synthesize the paper without introducing new results. It should restate the main contribution: direct VIAMD support for VeloxChem 1.0 HDF5 output. It should summarize the user-facing capabilities—loading, Summary, orbitals/densities, Orbital Grid, Response spectra, vibrational modes, transition analysis, critical points, and workspace restoration. It should then connect these features to the implementation: shared parsed VeloxChem data, GTO-basis reconstruction, AO remapping, and on-demand volume generation. It should close with future work: larger benchmark systems, GPU-enabled measurements, and broader documentation of additional VeloxChem response/property workflows.

## X. Availability, Documentation, and Datasets
This section should state that VIAMD is available under the MIT license at `https://github.com/scanberg/viamd`, with precompiled binaries for Windows and macOS and Linux build instructions in the documentation. It should mention active development and point users to the VIAMD wiki for documentation and tutorials. The TREXIO sentence should be integrated into the availability paragraph as an interoperability note, not placed between dataset entries.

The Datasets subsection should state that the article uses the VeloxChem input/output file collection from `https://veloxchem.org/docs/input-files/`. These examples provide Python inputs, text inputs, log outputs, and HDF5 outputs. The HDF5 files are regenerated regularly from the VeloxChem master branch and are used by the VIAMD developers as a practical feature-test collection for the VeloxChem interface.

The manuscript should include a compact feature matrix that follows the order of the VeloxChem input/output page. The table should use `Example` and `HDF5 file` as the identifying columns and should be placed on a landscape page for readability, with tilted/slim headers for feature columns. The text before the table should explicitly define the left-side Common column as the shared VIAMD starting point for VeloxChem HDF5 files: direct file loading, molecular-structure inspection, atom selection and highlighting through interactive tables, Summary-window metadata and convergence review, standard molecular representations, molecular-orbital/density/density-difference representations where electronic-structure data are present, Orbital Grid access where orbital data are present, molecular vectors or axes when present, density topology and critical-point inspection where density data allow it, and workspace saving/restoration. Additional columns should mark only the main dataset-specific VIAMD workflows with bullets rather than prose: Optimisation, Response, transition analysis, and Atomic props. The Optimisation column should be used for explicit optimization/path workflows such as S0/S1 optimization, frozen-coordinate examples, coordinate scans, TS searches, and IRC paths; it should not be marked for the reference-state SCF, ROHF, and UHF optimization rows unless those rows are intentionally used as optimization/path examples. The Atomic props. column should be restricted to atom-centered property datasets such as ESP, RESP, Boltzmann-weighted RESP, and LoProp. Do not create separate columns for molecular vectors or density topology because these are common/general VIAMD capabilities. Do not create a separate vibrational-window column; vibrational spectra and normal-mode interactions belong under the same Response workflow. The current table is a draft and should be verified against the monthly regenerated HDF5 files before final submission. The final public tutorial page should associate each file with the exact VIAMD features that can be inspected from it.

The table should include all major page categories in order: reference states, Hamiltonian, environment, potential energy surfaces, UV/Vis absorption/emission, first-order properties, polarizability, optical activity/dichroism, vibrational spectroscopies, weak interactions, localized properties, multi-photon interactions, X-ray spectroscopies, and general response functions. Rows with no HDF5 output listed on the VeloxChem page can be retained as log-file references and marked as not current VIAMD HDF5 visualization targets.

The paragraph after the dataset table should not feel like a detached figure-summary note. It should close the Datasets subsection by connecting the VeloxChem documentation collection and the matrix to reproducibility, regression testing, and future tutorial material.

## XI. Supporting Information
The Supporting Information should be aligned with the current article rather than the older generic VIAMD manuscript. It should include or point to example workspaces, benchmark source and CSV data, additional screenshots if needed, dataset provenance, and any practical limitations. The existing generic bullets for representation types, script types, script glossary, limitations, and example workspace can remain only if they are still relevant to this VeloxChem-focused article.

## XII. Final polish checklist
Before submission or sharing a full draft, check the following:
*   Add or guard the graphical table-of-contents image `figure/toc.png`.
*   Keep the IR workflow text-only unless a clear static representation of vibrational motion is produced later.
*   Replace future-tense captions such as "will combine" with final descriptive captions.
*   Replace "first-draft" wording in final-facing performance text if the benchmark is kept as an article result.
*   Confirm that "Transition density" in the timing figure is described as matrix extraction unless full transition-density volume generation is benchmarked.
*   Verify all bibliography keys and figure files during PDF compilation.
*   Finalize keywords, abbreviations, acknowledgements, dataset descriptions, and supporting-information text.

## XIII. Literature and citation plan
The manuscript now has a first literature pass. Keep the citations grouped at conceptual entry points rather than citing every implementation sentence. The Introduction should cite both classical molecular viewers and state-of-the-art molecular-visualization or wavefunction-analysis tools, including VMD, PyMOL, Molden, Chimera, JSmol, Mol*, Avogadro, Multiwfn, and the molecular-visualization survey. The Summary-window discussion of density critical points should cite QTAIM/electron-density topology via Bader. The Electronic Landscape section should cite NTOs and attachment/detachment density analysis when introducing excited-state volumetric descriptors, and it should cite ESP/RESP/LoProp references where atom-centered electrostatic coloring is discussed. The Spectroscopy section should cite TDDFT and response theory through Runge--Gross, Casida, and VeloxChem. The Transition Analysis section should cite VALET plus foundational excited-state analysis references for attachment/detachment densities, NTOs, charge-transfer descriptors, and TheoDORE-style excited-state analysis. The Implementation section should cite HDF5 when motivating the direct hierarchical data container. Further additions, if needed, should focus on CPP/X-ray spectroscopy only if those workflows become more than dataset-table entries.

## XIV. References
1. R. Skånberg, I. Hotz, A. Ynnerman, M. Linares, *VIAMD: a Software for Visual Interactive Analysis of Molecular Dynamics*, J. Chem. Inf. Model. 2023, 63, 23, 7382–7391.
2. R. Skånberg, C. König, P. Norman, M. Linares, D. Jönsson, I. Hotz, A. Ynnerman, *VIA-MD: Visual Interactive Analysis of Molecular Dynamics*, 2018, Eurographics Proceedings, p. 19–27.
3. R. Skånberg, M. Linares, M. Falk, I. Hotz, A. Ynnerman, *MolFind-Integrated Multi-Selection Schemes for Complex Molecular Structures*, 2019, The Eurographics Association, p. 17-21.
4. R. Skånberg, M. Falk, M. Linares, A. Ynnerman, I. Hotz, *Tracking Internal Frames of Reference for Consistent Molecular Distribution Functions*, 2021, IEEE Transactions on Visualization and Computer Graphics, 28 (9), 3126-3137.
