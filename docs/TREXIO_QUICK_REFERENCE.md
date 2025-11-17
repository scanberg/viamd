# TREXIO Support - Quick Reference

## Current Status: ❌ NOT IMPLEMENTED

### Key Findings

1. **mdlib TREXIO Bindings:** Not present
   - No `md_trexio.c` or `md_trexio.h` files
   - No TREXIO symbols in the codebase

2. **Build System:** No TREXIO support
   - No FindTREXIO.cmake
   - No TREXIO-related CMake options

3. **Test Files:** None available
   - No .trexio files in test_data directories

## How mdlib and viamd Locate Libraries

### Current Pattern for Optional Dependencies

mdlib uses a CMake-based pattern for optional file format support:

```cmake
# Example: VeloxChem support (from ext/mdlib/CMakeLists.txt)
option(MD_ENABLE_VLX "Enable Veloxchem Support" OFF)

if (${MD_ENABLE_VLX})
    set(MD_ENABLE_HDF5 TRUE)
endif()

if (MD_ENABLE_HDF5)
    find_package(HDF5 REQUIRED)
    set(MD_LIBS ${MD_LIBS} HDF5::HDF5)
endif()
```

### Integration Approach for TREXIO (Recommended)

**External Library Linking** (following VeloxChem pattern):

1. **System-wide installation:** TREXIO library installed via package manager or from source
2. **CMake detection:** FindTREXIO.cmake module locates the library
3. **Optional feature:** Controlled by `MD_ENABLE_TREXIO` CMake option
4. **Link-time:** mdlib links against TREXIO::TREXIO target

```cmake
# Proposed addition to ext/mdlib/CMakeLists.txt
option(MD_ENABLE_TREXIO "Enable TREXIO Support" OFF)

if (MD_ENABLE_TREXIO)
    find_package(TREXIO REQUIRED)
    set(SRC_FILES ${SRC_FILES} src/md_trexio.c src/md_trexio.h)
    set(MD_LIBS ${MD_LIBS} TREXIO::TREXIO)
    set(MD_DEFINES ${MD_DEFINES} MD_TREXIO)
endif()
```

## Currently Supported File Formats

mdlib currently supports these formats (built-in):
- PDB, XYZ, GRO, TRR, XTC, EDR, LAMMPS, mmCIF, CUBE, CSV, XVG, GTO, SMILES
- VLX (optional, requires HDF5)

## Test Files Available

Located in `ext/mdlib/test_data/`:
- Various molecular structure files (PDB, XYZ, GRO, etc.)
- Trajectory files (TRR, XTC, LAMMPS)
- Energy files (EDR, XVG)
- VeloxChem files (HDF5)
- **No TREXIO files currently**

## What's Needed for TREXIO Integration

1. **External TREXIO library** (https://github.com/TREX-CoE/trexio)
   - Can be installed via conda, package managers, or from source
   - Requires HDF5 for HDF5 backend support

2. **New mdlib files:**
   - `ext/mdlib/src/md_trexio.c` - Implementation
   - `ext/mdlib/src/md_trexio.h` - API declarations
   - `ext/mdlib/cmake/FindTREXIO.cmake` - CMake module

3. **Test data:**
   - Minimal TREXIO example files for validation

4. **CMake changes:**
   - Add `MD_ENABLE_TREXIO` option
   - Add find_package(TREXIO) logic
   - Update source file lists and link libraries

## See Also

- Full audit report: [docs/TREXIO_AUDIT.md](TREXIO_AUDIT.md)
- TREXIO project: https://github.com/TREX-CoE/trexio
- TREXIO documentation: https://trex-coe.github.io/trexio/
