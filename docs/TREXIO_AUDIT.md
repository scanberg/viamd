# TREXIO Support Audit - Phase 0

**Date:** November 17, 2025  
**Repository:** scanberg/viamd  
**Audited by:** GitHub Copilot  

## Executive Summary

This document provides an audit of TREXIO support in the viamd project, specifically examining whether the mdlib dependency (ext/mdlib) contains TREXIO C-bindings or if external TREXIO library linking is required.

**Key Finding:** mdlib does **NOT** currently contain TREXIO support. Integration will require linking against an external TREXIO library.

---

## Audit Methodology

1. Initialized all git submodules including ext/mdlib
2. Searched for md_trexio.* files in ext/mdlib/src
3. Searched for TREXIO references in source code and CMake configuration
4. Examined available test data files
5. Reviewed current file format support in mdlib

---

## Detailed Findings

### 1. TREXIO C-Bindings in ext/mdlib

**Status:** ❌ NOT PRESENT

- No `md_trexio.c` or `md_trexio.h` files exist in `ext/mdlib/src/`
- No TREXIO-related symbols found in the codebase
- No references to TREXIO in any source files

**Currently Supported File Formats in mdlib:**
- PDB (Protein Data Bank) - `md_pdb.c/h`
- XYZ - `md_xyz.c/h`
- GRO (GROMACS) - `md_gro.c/h`
- TRR (GROMACS Trajectory) - `md_trr.c/h`
- XTC (GROMACS Compressed Trajectory) - `md_xtc.c/h`
- EDR (GROMACS Energy) - `md_edr.c/h`
- LAMMPS - `md_lammps.c/h`
- mmCIF - `md_mmcif.c/h`
- CUBE - `md_cube.c/h`
- CSV - `md_csv.c/h`
- XVG - `md_xvg.c/h`
- VLX (VeloxChem) - `md_vlx.c/h` (optional, controlled by MD_ENABLE_VLX)
- GTO - `md_gto.c/h`
- SMILES - `md_smiles.c/h`

### 2. CMake Build Configuration

**Status:** ❌ NO TREXIO SUPPORT

- No `FindTREXIO.cmake` module exists
- No TREXIO-related options in CMakeLists.txt
- No TREXIO package dependencies configured

**Current External Dependencies in mdlib:**
- HDF5 (optional, used for VeloxChem support)
- NetCDF (optional)
- OpenGL (required for rendering)

**Relevant CMake configuration from ext/mdlib/CMakeLists.txt:**
```cmake
option(MD_ENABLE_VLX "Enable Veloxchem Support" OFF)

if (${MD_ENABLE_VLX})
    set(MD_ENABLE_HDF5 TRUE)
endif()

if (MD_ENABLE_HDF5)
    find_package(HDF5 REQUIRED)
    set(MD_LIBS ${MD_LIBS} HDF5::HDF5)
endif()
```

This pattern shows how external libraries are integrated into mdlib.

### 3. TREXIO Test/Example Files

**Status:** ❌ NONE AVAILABLE

- No `.trexio` files found in `datasets/` directory
- No TREXIO files in `ext/mdlib/test_data/` directory
- Test data includes: PDB, XYZ, GRO, TRR, XTC, EDR, LAMMPS, CUBE, XVG, and VeloxChem HDF5 files

**Available Test Data Formats:**
```
ext/mdlib/test_data/
├── *.pdb (PDB files)
├── *.xyz (XYZ files)
├── *.gro (GROMACS structure)
├── *.trr (GROMACS trajectory)
├── *.xtc (GROMACS compressed trajectory)
├── *.edr (GROMACS energy)
├── *.lammpstrj (LAMMPS trajectory)
├── *.cube (Gaussian CUBE)
├── *.xvg (XVG plot data)
├── *.arc (Tinker archive)
└── vlx/*.h5 (VeloxChem HDF5 files)
```

---

## TREXIO Background

TREXIO (https://github.com/TREX-CoE/trexio) is a file format and library designed for quantum chemistry calculations. It provides:

- Efficient storage of quantum chemistry data (basis sets, molecular orbitals, integrals, etc.)
- Support for both HDF5 and text-based backends
- Standardized format for quantum chemistry codes in the TREX (Targeting Real Chemical Accuracy at the Exascale) project

---

## Integration Plan Recommendations

### Option 1: Link Against External TREXIO Library (RECOMMENDED)

**Approach:**
1. Add TREXIO as an optional dependency similar to VeloxChem/HDF5
2. Create `FindTREXIO.cmake` or use pkg-config to locate system TREXIO installation
3. Add CMake option `MD_ENABLE_TREXIO` (similar to `MD_ENABLE_VLX`)
4. Implement `md_trexio.c/h` in mdlib to wrap TREXIO C API
5. Link against system-installed TREXIO library

**Advantages:**
- Leverages official TREXIO library with full feature support
- Reduced maintenance burden (bug fixes/features come from upstream)
- Ensures compatibility with TREXIO file format specification
- Can utilize TREXIO's optimized I/O implementations

**Disadvantages:**
- Adds external dependency that users must install
- May require different installation procedures per platform
- Version compatibility considerations

**Implementation Steps:**
```cmake
# In ext/mdlib/CMakeLists.txt
option(MD_ENABLE_TREXIO "Enable TREXIO Support" OFF)

if (MD_ENABLE_TREXIO)
    find_package(TREXIO REQUIRED)
    set(SRC_FILES ${SRC_FILES} src/md_trexio.c src/md_trexio.h)
    set(MD_LIBS ${MD_LIBS} TREXIO::TREXIO)
    set(MD_DEFINES ${MD_DEFINES} MD_TREXIO)
endif()
```

### Option 2: Vendor TREXIO Library

**Approach:**
1. Include TREXIO source as a git submodule in `ext/mdlib/ext/trexio`
2. Build TREXIO as part of mdlib build process
3. Implement md_trexio wrapper as above

**Advantages:**
- No external dependency for users
- Guaranteed version compatibility
- Simplified build process

**Disadvantages:**
- Larger repository size
- Must maintain TREXIO submodule
- May complicate packaging for distributions

### Recommended File Structure After Integration

```
ext/mdlib/
├── src/
│   ├── md_trexio.c         # New: TREXIO format support
│   ├── md_trexio.h         # New: TREXIO API declarations
│   └── ... (existing files)
├── cmake/
│   ├── FindTREXIO.cmake    # New: Locate TREXIO library
│   └── ... (existing files)
├── test_data/
│   ├── example.trexio/     # New: Sample TREXIO file(s)
│   └── ... (existing files)
└── CMakeLists.txt          # Modified: Add TREXIO option
```

### API Design Considerations

The `md_trexio.h` API should follow the pattern of other mdlib format handlers:

```c
// Example API structure (to be refined during implementation)

// Open/close functions
bool md_trexio_open(const char* filename, md_trexio_context_t* ctx);
void md_trexio_close(md_trexio_context_t* ctx);

// Read molecular structure
bool md_trexio_read_structure(md_trexio_context_t* ctx, md_molecule_t* mol);

// Read basis set information (if applicable)
bool md_trexio_read_basis(md_trexio_context_t* ctx, md_basis_t* basis);

// Read molecular orbitals (if applicable)
bool md_trexio_read_orbitals(md_trexio_context_t* ctx, md_orbitals_t* orbitals);
```

---

## Required Resources for Integration

### Dependencies
- **TREXIO Library**: v2.0+ (recommended latest stable)
  - Available via: GitHub, conda-forge, package managers
  - Requires: HDF5 (for HDF5 backend support)

### Test Files
- Need to acquire or generate minimal TREXIO test files
- Suggested sources:
  - TREXIO repository examples
  - Quantum chemistry codes that export TREXIO (e.g., Quantum Package, GAMESS)
  - Manually create minimal test files using TREXIO API

### Documentation
- TREXIO specification: https://trex-coe.github.io/trexio/
- API reference: https://trex-coe.github.io/trexio/trexio.html
- File format documentation

---

## Next Steps (Phase 1 - Implementation)

1. **Decision Point**: Choose integration approach (Option 1 vs Option 2)
2. **Acquire Test Files**: Obtain minimal TREXIO files for testing
3. **Create FindTREXIO.cmake**: CMake module to locate TREXIO library
4. **Add CMake Options**: Integrate TREXIO option into build system
5. **Implement md_trexio.c/h**: Create TREXIO format support module
6. **Add Tests**: Unit tests for TREXIO file reading
7. **Update Documentation**: User guide for TREXIO support
8. **CI/CD Updates**: Ensure TREXIO dependency available in build pipelines

---

## Conclusion

The current state of viamd/mdlib does **NOT** include TREXIO support. To add TREXIO functionality, the project will need to:

1. Link against the external TREXIO library (recommended approach)
2. Create new `md_trexio.c/h` modules following existing mdlib patterns
3. Add CMake configuration for optional TREXIO support
4. Acquire or generate test files for validation

The integration should follow the established pattern used for VeloxChem support, treating TREXIO as an optional feature that can be enabled via CMake configuration.

---

## References

- TREXIO GitHub: https://github.com/TREX-CoE/trexio
- TREXIO Documentation: https://trex-coe.github.io/trexio/
- mdlib Repository: https://github.com/scanberg/mdlib
- viamd Repository: https://github.com/scanberg/viamd

---

**Audit Status:** ✅ Complete  
**Phase:** 0 - Preparation/Quick Audit  
**Next Phase:** Implementation Planning
