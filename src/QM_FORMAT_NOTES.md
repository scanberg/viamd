# QM Source Integration Notes

VIAMD treats VeloxChem (`.h5`), Molden (`.molden`/`.mold`), and TREXIO (`.trexio`/`.h5`/`.hdf5`) files as sources of quantum-mechanical/electronic-structure data.

The GUI still shows the concrete source name (`VeloxChem`, `Molden`, or `TREXIO`) so users know where the data came from, but the runtime capabilities are shared conceptually:

- summary information
- orbital metadata
- orbital grid preview
- electronic-structure representations

## Runtime model

Each source-specific component is responsible for parsing and normalizing its own format. Shared rendering and representation code should receive data in VIAMD/mdlib-compatible order.

The intended boundary is:

```text
source file -> source adapter/parser -> canonical mdlib-compatible basis + MO coefficients -> shared rendering/evaluation path
```

Format-specific ordering rules should not leak into representation evaluation.

## AO ordering policy

### VIAMD/mdlib internal order

The internal GTO basis/evaluation path expects shell-contiguous ordering sorted by:

```text
angular momentum -> atom -> radial shell/function -> shell component
```

All source formats should map their AO/MO coefficients into this order before evaluation.

### VeloxChem

VeloxChem raw data uses a different AO layout. The mdlib VeloxChem reader remaps this into VIAMD/mdlib order before the component evaluates molecular orbitals.

### Molden

Molden has no fully reliable universal AO ordering convention across producers. The Molden component therefore defaults to automatic ordering resolution and only exposes manual overrides under source-specific advanced details.

Current automatic behavior:

- GANSU-generated files resolve to `Standard P Order`.
- files with only SP-derived p shells resolve to `SP-Derived P Order`.
- otherwise default to `Standard P Order`.

Manual overrides are for debugging unusual producers and should normally remain on `Auto`.

### TREXIO

TREXIO spherical AO ordering follows the official convention:

```text
m = 0,+1,-1,+2,-2,...
```

The TREXIO component maps this into VIAMD/mdlib shell-contiguous order. The default must remain `TREXIO spherical`; source/mdlib order is only a debug override for synthetic data already written in mdlib order.

## GUI policy

The visible GUI labels should show the concrete source name, not the generic internal concept:

- `VeloxChem Summary`, `VeloxChem Orbital Grid`
- `Molden Summary`, `Molden Orbital Grid`
- `TREXIO Summary`, `TREXIO Orbital Grid`

Source-specific details and manual AO-order overrides belong under `Source Details / Advanced`.

## PR validation checklist

Before merging changes in this area:

- build `viamd`
- load at least one Molden file with orbitals
- load the ORCA-derived TREXIO water example
- confirm orbital-grid generation works for both Molden and TREXIO
- confirm VeloxChem still exposes its existing response/NTO/export actions
- keep NTO-specific Molden behavior out of scope unless explicitly re-enabled
