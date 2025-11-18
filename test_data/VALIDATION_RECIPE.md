# TREXIO Test Validation Recipe

## Quick Validation

This recipe provides step-by-step instructions to validate VIAMD's TREXIO loader implementation.

## Prerequisites

- CMake >= 3.20
- C/C++ compiler (gcc, clang, or MSVC)
- Optional: HDF5 library for HDF5 backend support
- Optional: Python with `trexio` and `pyscf` packages for generating test files

## Step 1: Build VIAMD with TREXIO Support

```bash
# Clone and initialize submodules
cd /path/to/viamd
git submodule update --init --recursive

# Create build directory
mkdir build && cd build

# Configure with TREXIO enabled
cmake -DVIAMD_ENABLE_TREXIO=ON ..

# Build
make -j$(nproc)
```

## Step 2: Run Unit Tests

```bash
# From build directory
./bin/md_unittest --filter=trexio.*
```

**Expected Results:**
- ✅ `trexio.create_destroy` - PASS
- ✅ `trexio.disabled_compilation` - PASS
- ⚠️ `trexio.parse_h2_text` - Currently blocked (see Known Issues)
- ⚠️ `trexio.parse_h2o_text` - Currently blocked (see Known Issues)
- ⚠️ `trexio.system_init_h2` - Currently blocked (see Known Issues)
- ⚠️ `trexio.system_loader` - Currently blocked (see Known Issues)

## Step 3: Validate Test Data

### Check Test Files Exist

```bash
cd ../test_data
ls -l h2_molecule.trexio/ h2o_molecule.trexio/ ch4_molecule.trexio/
```

**Expected Files:**
Each `.trexio` directory should contain:
- `nucleus.txt` - Atomic coordinates and charges
- `electron.txt` - Electron configuration
- `metadata.txt` - Package information

### Validate File Contents

#### H2 Molecule (h2_molecule.trexio)
```bash
# Check atom count
grep -A 1 "^num$" h2_molecule.trexio/nucleus.txt
# Expected: 2

# Check charges
grep -A 2 "^charge$" h2_molecule.trexio/nucleus.txt
# Expected: 1.0, 1.0 (two hydrogen atoms)

# Check labels
grep -A 2 "^label$" h2_molecule.trexio/nucleus.txt
# Expected: H, H

# Check electrons
cat h2_molecule.trexio/electron.txt
# Expected: up_num=1, dn_num=1
```

#### H2O Molecule (h2o_molecule.trexio)
```bash
# Check atom count
grep -A 1 "^num$" h2o_molecule.trexio/nucleus.txt
# Expected: 3

# Check charges
grep -A 3 "^charge$" h2o_molecule.trexio/nucleus.txt
# Expected: 8.0 (O), 1.0 (H), 1.0 (H)

# Check electrons
cat h2o_molecule.trexio/electron.txt
# Expected: up_num=5, dn_num=5 (total 10 electrons)
```

#### CH4 Molecule (ch4_molecule.trexio)
```bash
# Check atom count
grep -A 1 "^num$" ch4_molecule.trexio/nucleus.txt
# Expected: 5 (1 C + 4 H)

# Check charges
grep -A 5 "^charge$" ch4_molecule.trexio/nucleus.txt
# Expected: 6.0 (C), 1.0, 1.0, 1.0, 1.0 (4 H)
```

## Step 4: Low-Level Parsing Function Tests

These tests are included in the unit test suite (`test_trexio.c`):

### Test: TREXIO Object Creation/Destruction
- **Function:** `md_trexio_create()`, `md_trexio_destroy()`
- **Validation:** Object is created and destroyed without memory leaks
- **Status:** ✅ PASSING

### Test: Atom Count Parsing
- **Function:** `md_trexio_number_of_atoms()`
- **Validation:** Returns correct count: 2 for H2, 3 for H2O, 5 for CH4
- **Status:** ⚠️ BLOCKED (file format issue)

### Test: Coordinate Parsing
- **Function:** `md_trexio_atom_coordinates()`
- **Validation:** 
  - Returns coordinates in Angstrom (converted from Bohr)
  - H2: Second atom at ~0.741 Angstrom from first
- **Status:** ⚠️ BLOCKED (file format issue)

### Test: Charge Parsing
- **Function:** `md_trexio_atomic_charges()`
- **Validation:** Returns correct charges for each atom type
- **Status:** ⚠️ BLOCKED (file format issue)

### Test: AO/MO Dimensions
- **Function:** `md_trexio_number_of_aos()`, `md_trexio_mo_num()`
- **Validation:** Returns 0 for geometry-only files (not yet implemented for full quantum files)
- **Status:** Not yet tested (requires files with basis set data)

### Test: System Initialization
- **Function:** `md_trexio_system_init()`
- **Validation:** Populates `md_system_t` structure correctly
- **Status:** ⚠️ BLOCKED (file format issue)

## Step 5: Manual Verification (Optional)

### Load Test Files in VIAMD GUI

```bash
cd build
./viamd ../test_data/h2_molecule.trexio
```

**Expected Behavior:**
- VIAMD loads without errors
- Two hydrogen atoms visible
- Bond length approximately 0.74 Angstrom
- Atoms colored according to element (hydrogen typically white/gray)

**Known Limitation:** Currently blocked due to TREXIO text format compatibility issue.

## Known Issues

### TREXIO Text Format Compatibility

**Issue:** TREXIO 2.6.0 library returns `TREXIO_READONLY` error when opening text format files created manually.

**Error Message:**
```
TREXIO: Failed to open file: Read-only file
```

**Root Cause:** The text format files may be missing required metadata or version information that TREXIO 2.6.0 expects.

**Workarounds:**

1. **Generate files with TREXIO library:**
   ```bash
   pip install trexio pyscf
   cd test_data
   python3 create_pyscf_trexio.py
   ```

2. **Use HDF5 backend:**
   ```bash
   cmake -DVIAMD_ENABLE_TREXIO=ON -DENABLE_HDF5=ON ..
   make
   # Generate .h5 files instead of .trexio text format
   ```

3. **Use files from quantum chemistry codes:**
   - Export from Quantum Package: `qp_run export_trexio molecule.trexio`
   - Use PySCF with TREXIO bindings
   - Use other supported codes (CP2K, FHI-aims, etc.)

## Test Data Sources

### Current Test Files

All test files are located in `test_data/` directory:

| File | Description | Atoms | Electrons | Source |
|------|-------------|-------|-----------|--------|
| `h2_molecule.trexio` | Hydrogen molecule | 2 (H, H) | 2 | Created with `create_test_trexio.py` |
| `h2o_molecule.trexio` | Water molecule | 3 (O, H, H) | 10 | Created with `create_test_trexio.py` |
| `ch4_molecule.trexio` | Methane molecule | 5 (C, 4×H) | 10 | Created with `create_test_trexio.py` |

### Reference Data

Coordinates are specified in Bohr (atomic units) in TREXIO files and should be converted to Angstrom for display:

- Conversion factor: 1 Bohr = 0.529177210903 Angstrom
- H2 bond length: 1.4 Bohr ≈ 0.741 Angstrom
- H2O bond length: ~1.8 Bohr ≈ 0.95 Angstrom
- H2O bond angle: ~104.5 degrees

## Acceptance Criteria

- [x] **Sample TREXIO files collected:** 3 files (H2, H2O, CH4) ✅
- [x] **Test recipe documented:** This file + TESTING_GUIDE.md ✅
- [x] **Unit tests created:** 6 tests in `ext/mdlib/unittest/test_trexio.c` ✅
- [x] **Low-level parsing functions tested:** Tests implemented ✅
- [ ] **All tests passing:** 2/6 passing, 4/6 blocked by file format issue ⚠️
- [x] **Documentation complete:** README.md, TESTING_GUIDE.md, this recipe ✅

## Next Steps

1. **Resolve TREXIO text format compatibility:**
   - Investigate TREXIO 2.6.0 text format requirements
   - Update test file generation script
   - Or switch to HDF5 format for testing

2. **Generate quantum chemistry test files:**
   - Create files with basis set data (geometry + AO)
   - Create files with MO data (geometry + basis + MOs)
   - Test with `create_pyscf_trexio.py`

3. **Manual GUI testing:**
   - Load each test file in VIAMD
   - Verify geometry rendering
   - Check atom labels and bonds
   - Screenshot for documentation

4. **Performance testing:**
   - Measure file loading time
   - Compare text vs HDF5 backend
   - Test with larger molecules (50+ atoms)

## Support and References

- **TREXIO Documentation:** https://trex-coe.github.io/trexio/
- **TREXIO GitHub:** https://github.com/TREX-CoE/trexio
- **Issue Tracking:** See repository issues tagged with `TREXIO` and `testing`
- **Test Scripts:**
  - `test_data/create_test_trexio.py` - Simple text format generator
  - `test_data/create_pyscf_trexio.py` - Full quantum chemistry data generator
  - `test_data/validate_trexio.sh` - Validation script
  - `test_data/build_and_test.sh` - Automated build and test
