# TREXIO Integration Guide for VIAMD

This document describes the implementation of TREXIO file support in VIAMD.

## Overview

TREXIO (https://github.com/TREX-CoE/trexio) is an open-source file format and library for quantum chemistry data, developed by the TREX Centre of Excellence. It provides a standardized way to store and exchange data such as molecular geometries, basis sets, molecular orbitals, and other quantum chemical properties.

## Current Implementation Status

### ✅ Completed (Phase 1 - Infrastructure)

1. **CMake Integration**
   - Added `VIAMD_ENABLE_TREXIO` CMake option (default: OFF)
   - Created `cmake/FindTREXIO.cmake` module to locate TREXIO library
   - Updated main `CMakeLists.txt` to conditionally compile TREXIO support
   - Added TREXIO component to build system

2. **TREXIO Component Skeleton**
   - Created `src/components/trexio/trexio.cpp`
   - Implemented basic file reading using TREXIO C API
   - Added component menu and info dialogs
   - Reads basic metadata: atom count, electron count, basis set info, MO count

3. **Loader System Registration**
   - Updated `src/loader.cpp` to recognize `.trexio` file extension
   - Added TREXIO to loader enumerations and tables
   - Prepared infrastructure for full loader integration

4. **Documentation**
   - Updated `README.md` with build instructions
   - Documented TREXIO library prerequisites
   - Provided installation examples for various platforms

### 🚧 TODO (Phase 2 - Full Integration)

To complete TREXIO support and enable full molecule/wavefunction visualization, the following work is needed:

#### 1. Implement mdlib TREXIO Loader

Create `ext/mdlib/src/md_trexio.c` and `ext/mdlib/src/md_trexio.h` following the pattern of `md_vlx.c`:

**Header Structure (`md_trexio.h`):**
```c
#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <md_types.h>
#include <core/md_str.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_molecule_t;
struct md_molecule_loader_i;

typedef struct md_trexio_t md_trexio_t;

// Core API
struct md_trexio_t* md_trexio_create(struct md_allocator_i* alloc);
void md_trexio_destroy(struct md_trexio_t* trexio);
bool md_trexio_parse_file(struct md_trexio_t* trexio, str_t filename);

// Molecule data accessors
size_t md_trexio_number_of_atoms(const struct md_trexio_t* trexio);
size_t md_trexio_number_of_electrons(const struct md_trexio_t* trexio);
const dvec3_t* md_trexio_atom_coordinates(const struct md_trexio_t* trexio);
const uint8_t* md_trexio_atomic_numbers(const struct md_trexio_t* trexio);

// Basis set accessors
size_t md_trexio_basis_shell_num(const struct md_trexio_t* trexio);
// ... add more as needed

// Molecule loader interface
md_molecule_loader_i* md_trexio_molecule_api(void);

#ifdef __cplusplus
}
#endif
```

**Implementation (`md_trexio.c`):**
```c
#include "md_trexio.h"
#include <trexio.h>
#include <md_molecule.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>

// Internal structure
struct md_trexio_t {
    md_allocator_i* alloc;
    
    // Atom data
    size_t num_atoms;
    dvec3_t* atom_coords;
    uint8_t* atomic_numbers;
    
    // Electron data
    size_t num_electrons;
    size_t num_alpha;
    size_t num_beta;
    
    // Basis set data
    size_t num_shells;
    // ... add basis set arrays
    
    // MO data
    size_t num_mo;
    // ... add MO coefficient arrays
};

// Implementation of functions...
static bool trexio_mol_init_from_file(md_molecule_t* mol, str_t filename, 
                                      const void* arg, md_allocator_i* alloc) {
    (void)arg;
    md_trexio_t* trexio = md_trexio_create(md_get_heap_allocator());
    
    bool success = false;
    if (md_trexio_parse_file(trexio, filename)) {
        success = md_trexio_molecule_init(mol, trexio, alloc);
    }
    
    md_trexio_destroy(trexio);
    return success;
}

static md_molecule_loader_i trexio_loader = {
    NULL,  // init_from_str
    trexio_mol_init_from_file
};

md_molecule_loader_i* md_trexio_molecule_api(void) {
    return &trexio_loader;
}
```

#### 2. Update src/loader.cpp

Replace the NULL placeholder with actual loader:
```cpp
#if MD_TREXIO
#include <md_trexio.h>
#endif

// In mol_loader_api array:
#if MD_TREXIO
    md_trexio_molecule_api(),
#endif
```

#### 3. Enhance TREXIO Component

Expand `src/components/trexio/trexio.cpp` to include:
- Orbital visualization (following VeloxChem pattern)
- Basis set display
- Electronic structure information
- Interactive property plots
- Volume rendering for orbitals

#### 4. Testing

Create test cases with sample TREXIO files:
- Benzene molecule (from TREXIO test suite)
- Water molecule (minimal example)
- Test both HDF5 and text backends
- Verify coordinate conversion (Bohr → Angstrom)
- Validate basis set parsing
- Check MO coefficient loading

## TREXIO Data Mapping

| TREXIO Group | TREXIO Field | VIAMD Structure | Notes |
|--------------|--------------|-----------------|-------|
| nucleus | num | mol->atom.count | Number of atoms |
| nucleus | charge | mol->atom.element | Atomic numbers (Z) |
| nucleus | coord | mol->atom.x/y/z | Coordinates (Bohr → Å) |
| nucleus | label | mol->atom.type | Element symbols |
| electron | num | N/A | Total electrons |
| electron | up_num | N/A | Alpha electrons |
| electron | dn_num | N/A | Beta electrons |
| basis | type | N/A | Basis set type (GTO) |
| basis | shell_num | gto_data | Number of shells |
| basis | prim_num | gto_data | Number of primitives |
| basis | nucleus_index | gto_data | Shell centers |
| basis | shell_ang_mom | gto_data | Angular momentum (s,p,d,...) |
| basis | exponent | gto_data | Gaussian exponents |
| basis | coefficient | gto_data | Contraction coefficients |
| mo | num | N/A | Number of MOs |
| mo | coefficient | MO arrays | MO coefficients |
| mo | energy | MO arrays | Orbital energies |
| mo | occupation | MO arrays | Occupations |

## Building with TREXIO Support

### Prerequisites

**Ubuntu/Debian:**
```bash
sudo apt-get install libtrexio-dev
```

**From source:**
```bash
git clone https://github.com/TREX-CoE/trexio.git
cd trexio
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local
make -j4
make install
```

### Build VIAMD

```bash
cd viamd
git submodule update --init --recursive
cmake -DVIAMD_ENABLE_TREXIO=ON .
make -j4
```

If TREXIO is in a custom location:
```bash
cmake -DVIAMD_ENABLE_TREXIO=ON -DTREXIO_ROOT=/path/to/trexio .
```

## Example Usage

Once fully implemented, users will be able to:

1. Load TREXIO files via File → Open or drag-and-drop
2. View molecular structure with proper atom types and coordinates
3. Visualize molecular orbitals (HOMO, LUMO, etc.)
4. Display basis set information
5. Plot orbital energies and occupation numbers
6. Export visualizations and data

## References

- TREXIO Documentation: https://trex-coe.github.io/trexio/
- TREXIO GitHub: https://github.com/TREX-CoE/trexio
- TREXIO Paper: https://doi.org/10.1063/5.0148161
- VeloxChem Integration (similar pattern): `src/components/veloxchem/veloxchem.cpp`

## Contributing

To contribute to TREXIO support:

1. Implement mdlib loader (Phase 2, Step 1)
2. Test with real TREXIO files
3. Add visualization features
4. Document any TREXIO-specific considerations
5. Submit pull request with tests

## Notes

- TREXIO coordinates are typically in Bohr; convert to Angstrom for VIAMD
- TREXIO supports both HDF5 and text backends; handle both
- Some TREXIO data is optional; handle missing fields gracefully
- Follow VeloxChem integration pattern for consistency
- Maintain backward compatibility (TREXIO support optional)
