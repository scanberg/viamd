# TREXIO Integration - Known Issues and Required Fixes

## Current Status

TREXIO support in VIAMD has been implemented but **is not yet fully functional** due to API compatibility issues between the TREXIO loader code and the current mdlib API.

## Known Issues

### 1. Compilation Errors Due to mdlib API Changes

**Status:** ❌ **BLOCKING** - Prevents compilation with TREXIO enabled

**Description:**
The TREXIO loader implementation (`md_trexio.c`) uses an outdated mdlib API. Specifically, the `md_atom_type_data_t` structure no longer has an `element` field, which is referenced in the TREXIO loader code.

**Error Messages:**
```
error: 'md_atom_type_data_t' has no member named 'element'
    if (sys->atom.type.element && sys->atom.type.element[j] == element) {
                      ^
```

**Location:** `ext/mdlib/src/md_trexio.c`, function `md_trexio_system_init()`

**Impact:** Cannot build VIAMD with `-DVIAMD_ENABLE_TREXIO=ON`

**Required Fix:**
- Update md_trexio.c to use current mdlib API
- Investigate how element information is now stored in `md_atom_type_data_t`
- Update atom type assignment logic accordingly
- Regenerate the mdlib_trexio patch after fixing

### 2. Type Warnings - int32_t vs int64_t

**Status:** ⚠️  **Non-blocking** - Compiles with warnings but may cause issues

**Description:**
TREXIO API uses `int32_t` for counts (nucleus_num, mo_num, electron counts) but the md_trexio structure uses `int64_t`.

**Warning Messages:**
```
warning: passing argument 2 of 'trexio_read_nucleus_num' from incompatible pointer type
warning: passing argument 2 of 'trexio_read_mo_num' from incompatible pointer type
warning: passing argument 2 of 'trexio_read_electron_up_num' from incompatible pointer type
```

**Impact:** 
- On 64-bit systems, usually works but technically undefined behavior
- Could cause issues with very large values or on different architectures

**Required Fix:**
- Use temporary `int32_t` variables in TREXIO read calls
- Cast to `int64_t` for storage in md_trexio_t structure
- Or change md_trexio_t structure to use `int32_t` if sizes permit

### 3. Corrupted Patch Files

**Status:** ✅ **FIXED** - Workaround implemented

**Description:**
The main patch file (`mdlib_trexio.patch`) appears to be corrupted at line 919, causing `git apply` to fail.

**Solution Implemented:**
- Patch application script now uses `mdlib_trexio_original.patch` instead
- Added `FindTREXIO.cmake` to mdlib cmake directory
- Created `mdlib_trexio_updated.patch` with FindTREXIO module included

**Remaining Work:**
- Test the updated patch thoroughly
- Update mdlib_trexio.c to fix API compatibility
- Regenerate final clean patch file

## What Works

- ✅ Documentation is accurate and complete
- ✅ TREXIO library detection via CMake FindTREXIO module
- ✅ Patch application script is idempotent and robust
- ✅ Build instructions are correct (when TREXIO support is disabled)
- ✅ Test data files exist and are valid
- ✅ UI component structure is in place (`src/components/trexio/trexio.cpp`)

## What Doesn't Work

- ❌ Building with `-DVIAMD_ENABLE_TREXIO=ON` fails due to compilation errors
- ❌ TREXIO file loading cannot be tested (loader doesn't compile)
- ❌ Integration tests cannot run
- ❌ Molecular orbital visualization is incomplete

## Required Steps to Fix

### Priority 1: Fix Compilation (Required for basic functionality)

1. **Update md_trexio.c to use current mdlib API**
   - File: `ext/mdlib/src/md_trexio.c`
   - Function: `md_trexio_system_init()`
   - Lines: ~470-495 (element assignment logic)
   
2. **Investigate current mdlib atom type API**
   ```bash
   cd ext/mdlib/src
   grep -r "md_atom_type_data_t" *.h
   # Find how element information is currently stored
   ```

3. **Fix type mismatches**
   - Use `int32_t` temporaries for TREXIO API calls
   - Document why int64_t is used in md_trexio_t if there's a reason

4. **Test compilation**
   ```bash
   ./scripts/apply_mdlib_trexio_patch.sh
   mkdir build && cd build
   cmake -DVIAMD_ENABLE_TREXIO=ON ..
   make
   ```

5. **Regenerate patch**
   ```bash
   cd ext/mdlib
   git add src/md_trexio.* cmake/FindTREXIO.cmake
   git diff --cached HEAD > ../../docs/mdlib_trexio_final.patch
   # Also capture CMakeLists.txt changes
   ```

### Priority 2: Verify Functionality (After compilation works)

1. **Test with sample files**
   ```bash
   ./viamd ../test_data/h2_molecule.trexio
   ./viamd ../test_data/h2o_molecule.trexio
   ```

2. **Run unit tests** (if they exist)
   ```bash
   cd build
   ctest -V -R trexio
   ```

3. **Test with HDF5 format files**
   - Requires HDF5-enabled TREXIO installation
   - Test with PySCF-generated files

### Priority 3: Complete Integration (After basic loading works)

1. **Complete molecular orbital visualization**
   - Update `src/components/trexio/trexio.cpp`
   - Implement grid-based orbital rendering
   - Add UI controls for orbital selection

2. **Add more data group support**
   - Basis set information display
   - Electron density visualization
   - Additional quantum chemistry properties

3. **Performance optimization**
   - Asynchronous file loading
   - Caching strategies
   - Large file handling

## Testing Checklist

Once compilation is fixed, verify:

- [ ] VIAMD builds successfully with `-DVIAMD_ENABLE_TREXIO=ON`
- [ ] Can load H2 test file without crashes
- [ ] Can load H2O test file without crashes
- [ ] Atomic coordinates are displayed correctly
- [ ] Atom counts match expected values
- [ ] Electron counts are correct
- [ ] Can open both .trexio (text) and .h5 (HDF5) formats
- [ ] Error handling works for invalid files
- [ ] TREXIO summary window displays correct information
- [ ] Build works on Linux (tested: Ubuntu 22.04, Ubuntu 24.04)
- [ ] Build works on macOS
- [ ] Build works on Windows (if tested)

## Documentation Status

All user-facing documentation has been updated and is accurate:

- ✅ README.md - Build instructions, prerequisites
- ✅ docs/TREXIO_SUPPORT.md - Complete user guide
- ✅ docs/PHASE5_COMPONENT.md - Component architecture
- ✅ scripts/apply_mdlib_trexio_patch.sh - Working patch script
- ✅ test_data/README.md - Test file documentation

Documentation correctly states that:
- TREXIO must be installed separately (not auto-downloaded)
- HDF5 is optional but recommended
- Current limitations and known issues

## Recommendations for Next Developer

1. **Start with API investigation**: Don't modify md_trexio.c blindly. First understand how the current mdlib API handles element information.

2. **Check mdlib changelog/commits**: Look for commits that changed `md_atom_type_data_t` structure to understand the reason for the change.

3. **Consider using md_vlx.c as reference**: VeloxChem loader (md_vlx.c) should use the current API and can serve as a template.

4. **Test incrementally**: Fix one compilation error at a time, test after each fix.

5. **Update patch carefully**: The patch must apply cleanly to mdlib at commit `06da5a3`. Test patch application on fresh mdlib checkout.

6. **Keep documentation updated**: Update this file and other docs as you fix issues.

## References

- mdlib repository: https://github.com/scanberg/mdlib
- Current mdlib commit: `06da5a3` (vlx fixes)
- TREXIO library: https://github.com/TREX-CoE/trexio
- TREXIO version tested: 2.6.0
- Related documentation:
  - TREXIO_IMPLEMENTATION_STATUS.md
  - PHASE4_SUMMARY.md
  - test_data/TESTING_GUIDE.md
  - test_data/VALIDATION_RECIPE.md

## Contact

For questions about TREXIO integration:
- VIAMD Issues: https://github.com/scanberg/viamd/issues
- TREXIO Issues: https://github.com/TREX-CoE/trexio/issues
