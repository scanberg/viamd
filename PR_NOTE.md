# PR Note: QM Source Integration, TREXIO Support, and AO Ordering

## Summary

This PR extends VIAMD's quantum/electronic-structure file support and makes the user-facing QM workflow more consistent across source formats.

The main changes are:

- add TREXIO file support with summary and orbital-grid visualization
- improve Molden AO-order handling and simplify its UI
- treat VeloxChem, Molden, and TREXIO as QM/electronic-structure sources internally
- keep the GUI source-specific, so users still see `VeloxChem`, `Molden`, or `TREXIO`
- add validation examples and loader tests for Molden
- document the AO-ordering policy for future maintenance

## Molden

Molden support was cleaned up around AO ordering and validation.

### User-facing changes

- Molden still appears as `Molden` in the GUI.
- The Molden summary window is now source-focused and less verbose.
- AO ordering is automatic for normal use.
- Manual AO-order override remains available under `Source Details / Advanced`.

### AO-order handling

Molden files vary by producer, especially for p shells and SP-derived p functions. The default is now automatic resolution:

- GANSU-generated files resolve to `Standard P Order`.
- files with only SP-derived p shells resolve to `SP-Derived P Order`.
- otherwise Molden defaults to `Standard P Order`.

Manual modes are kept only as an advanced/debug override:

- `Auto`
- `Source Order`
- `Standard P Order`
- `SP-Derived P Order`

The Molden component resets its settings between files so an override from one file does not leak into the next file.

### Validation

The Molden loader test was updated to current mdlib APIs and extended with new examples:

- `h2o_sto3g.molden`
- `h2_sto3g.molden`
- `Ammonia_NH3.molden`
- `Anthracene_C14H10.molden`

All intended Molden examples parse and load successfully.

The NTO Molden example is intentionally excluded from this PR.

## TREXIO

This PR adds optional TREXIO support behind `VIAMD_ENABLE_TREXIO`.

### Build integration

- CMake detects TREXIO through `pkg-config` or a direct header/library search.
- If TREXIO is unavailable, TREXIO support is disabled cleanly.
- The TREXIO library directories are passed to the linker to handle installations where `pkg-config` exposes dependent libraries such as HDF5 separately.

### Loader integration

TREXIO files are treated as QM system sources. Detection supports:

- `.trexio`
- `.h5`
- `.hdf5`

For HDF5-like extensions, the loader probes the file before selecting the TREXIO loader so it does not steal non-TREXIO HDF5 files such as VeloxChem files.

### GUI and representations

TREXIO exposes:

- `TREXIO Summary`
- `TREXIO Orbital Grid`
- molecular-orbital electronic-structure representations

The GUI shows the source as `TREXIO`, not generic `QM`.

### TREXIO example

A proper ORCA-derived TREXIO water example is included under `datasets/trexio_examples/`.

This file is independent of the VeloxChem reference path and is intended as the TREXIO validation example.

## AO-ordering policy

The central rule is:

> Source-specific parsers/adapters must map AO/MO coefficients into VIAMD/mdlib-compatible order before evaluation.

The rendering and representation path should not need to know about producer-specific AO quirks.

### VIAMD/mdlib internal order

The internal GTO evaluation path expects shell-contiguous ordering sorted by:

```text
angular momentum -> atom -> radial shell/function -> shell component
```

### VeloxChem

VeloxChem raw data has its own ordering. mdlib's VeloxChem reader remaps the raw VeloxChem AO/MO data into VIAMD/mdlib order before orbital evaluation.

### Molden

Molden ordering is producer-dependent. The Molden component resolves this automatically by default and exposes manual override only under advanced source details.

### TREXIO

TREXIO spherical AO order follows the official convention:

```text
m = 0,+1,-1,+2,-2,...
```

The TREXIO component maps this official spherical order into VIAMD/mdlib shell-contiguous order. The default TREXIO mode must remain `TREXIO spherical`.

`Source / mdlib order` and `Source p: px,py,pz` are debug overrides for unusual or synthetic files only.

## QM layer / UI organization

The implementation now treats VeloxChem, Molden, and TREXIO as sources of QM/electronic-structure data. A small shared UI helper centralizes common source-window labels and source display.

Important UI decision:

- internally, these formats share the QM/electronic-structure concept
- externally, the GUI continues to show the concrete source name

So users see:

- `VeloxChem Summary`, `VeloxChem Orbital Grid`
- `Molden Summary`, `Molden Orbital Grid`
- `TREXIO Summary`, `TREXIO Orbital Grid`

Source-specific details and AO-order overrides belong under `Source Details / Advanced`.

## Validation performed

- full VIAMD build succeeds with TREXIO enabled
- Molden parser validation passed for the intended examples
- Molden loader validation passed for H2O, H2, ammonia, and anthracene examples
- TREXIO water example was generated from ORCA JSON using official TREXIO tooling and validated through the TREXIO pathe concrete source name


## Out of scope

- Molden NTO support
- one-off helper tools for VLX-to-TREXIO conversion or TREXIO probing
- broad unification of all source-specific backends into a single canonical data structure

The current PR keeps source-specific parsing/remapping localized and unifies the GUI/representation behavior where it is safe to do so.
