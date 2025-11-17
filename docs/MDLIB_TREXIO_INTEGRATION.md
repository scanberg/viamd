# TREXIO Integration - mdlib Changes

## Overview

This directory contains patches and instructions for integrating TREXIO support into the mdlib library.

## Files

- `mdlib_trexio.patch` - Git patch file containing all mdlib changes for TREXIO support
- `TREXIO_INTEGRATION_DESIGN.md` - Technical design document (in /tmp for review)
- `TREXIO_SUPPORT.md` - User documentation

## Applying the mdlib Patch

### Option 1: Apply Patch to mdlib Submodule (Recommended for Development)

```bash
cd ext/mdlib
git checkout -b add-trexio-support
git apply ../../docs/mdlib_trexio.patch
git add .
git commit -m "Add TREXIO file format support"
```

### Option 2: Apply Patch Directly (For Testing)

```bash
cd ext/mdlib
patch -p1 < ../../docs/mdlib_trexio.patch
```

### Option 3: Manual File Addition

Copy the following files to `ext/mdlib/src/`:
- `md_trexio.h` - TREXIO header file
- `md_trexio.c` - TREXIO implementation

Then apply the CMakeLists.txt changes manually:

1. Add option: `option(MD_ENABLE_TREXIO "Enable TREXIO Support" OFF)`
2. Add TREXIO library detection and linking
3. Add md_trexio.c and md_trexio.h to source files when enabled

## mdlib Changes Summary

### New Files

**src/md_trexio.h**
- Public API for TREXIO file loading
- System loader interface
- Accessor functions for quantum chemistry data
- Follows same pattern as md_vlx.h

**src/md_trexio.c**
- Implementation of TREXIO file parsing
- Reads nucleus, basis, MO, and electron data
- Converts TREXIO data to md_system_t
- Handles both HDF5 and text backends
- ~600 lines of C code

### Modified Files

**CMakeLists.txt**
- Added `MD_ENABLE_TREXIO` option
- Added TREXIO library detection
- Added md_trexio.c/h to source files when enabled
- Added MD_TREXIO preprocessor definition

## Integration Status

### ✅ Completed
- Core TREXIO file parsing
- Nucleus data extraction
- Basis set data reading
- MO data reading
- System loader interface
- Conditional compilation support

### 🚧 In Progress
- GTO extraction implementation
- Basis set to md_gto_data_t conversion
- MO to GTO conversion

### 📋 Future Work
- Excited state data
- Response properties
- Density matrix extraction
- Unit tests

## Dependencies

The mdlib TREXIO support requires:

1. **TREXIO library** (>= 2.0.0)
   - Install via package manager or from source
   - See: https://github.com/TREX-CoE/trexio

2. **HDF5 library** (optional, for HDF5 backend)
   - Recommended for performance
   - Required for `.h5` file support

## Contributing to mdlib

If you maintain or have access to the mdlib repository, consider:

1. **Submit Pull Request**:
   - Fork https://github.com/scanberg/mdlib
   - Create feature branch
   - Apply this patch
   - Submit PR with detailed description

2. **Coordinate with Maintainer**:
   - Discuss integration approach
   - Ensure compatibility with mdlib design
   - Follow mdlib coding standards

3. **Alternative: Local Fork**:
   - Maintain a local fork with TREXIO support
   - Update .gitmodules to point to your fork
   - Keep synchronized with upstream

## Testing

After applying the patch, test with:

```bash
cd ext/mdlib
mkdir build && cd build
cmake -DMD_ENABLE_TREXIO=ON ..
make
# Run unit tests if available
ctest
```

## Patch Contents

The patch includes:
- Complete md_trexio.h header (122 lines)
- Complete md_trexio.c implementation (600+ lines)
- CMakeLists.txt modifications (20 lines added)

Total additions: ~750 lines of code

## Notes

1. **Submodule Consideration**: mdlib is a git submodule. Changes here should ideally be contributed upstream to the mdlib repository.

2. **Detached HEAD**: The submodule is currently in detached HEAD state at commit `06da5a3`. You may need to checkout a branch before making changes.

3. **Build System**: The TREXIO support is optional and disabled by default. It requires explicit `-DMD_ENABLE_TREXIO=ON` flag.

4. **Compatibility**: The implementation follows the same patterns as existing loaders (md_vlx, md_pdb, etc.) for consistency.

5. **Error Handling**: Comprehensive error handling with MD_LOG_ERROR and MD_LOG_WARN for debugging.

## License

All mdlib changes follow the mdlib license (same as VIAMD project).

## Support

For questions about the mdlib integration:
- VIAMD Issues: https://github.com/scanberg/viamd/issues
- mdlib Repository: https://github.com/scanberg/mdlib

For questions about TREXIO:
- TREXIO Issues: https://github.com/TREX-CoE/trexio/issues
