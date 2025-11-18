# TREXIO Test Data

This directory contains test TREXIO files for validating VIAMD's TREXIO support.

## Quick Links

- 📘 **[Complete Testing Guide](TESTING_GUIDE.md)** - Comprehensive testing documentation
- ✅ **[Validation Recipe](VALIDATION_RECIPE.md)** - Step-by-step validation instructions
- 🔨 **[Build Script](build_and_test.sh)** - Automated build and test
- ✅ **[Validation Script](validate_trexio.sh)** - Quick setup verification
- 🧪 **[PySCF Generator](create_pyscf_trexio.py)** - Create files with quantum chemistry data

## Test Status

**Current Status:** Phase 4 Complete - Test infrastructure ready

- ✅ Sample TREXIO files: 3 files (H2, H2O, CH4)
- ✅ Unit tests created: 6 tests in `ext/mdlib/unittest/test_trexio.c`
- ✅ Documentation complete: TESTING_GUIDE.md, VALIDATION_RECIPE.md
- ⚠️ Test execution: 2/6 passing, 4/6 blocked by file format compatibility

**Known Issue:** TREXIO 2.6.0 cannot open manually-created text format files.
See [VALIDATION_RECIPE.md](VALIDATION_RECIPE.md#known-issues) for details and workarounds.

## Test Files

### Text Format (.trexio directories)

- **h2_molecule.trexio** - Simple H2 molecule (2 atoms)
- **h2o_molecule.trexio** - Water molecule (3 atoms)
- **ch4_molecule.trexio** - Methane molecule (5 atoms)

### File Format

TREXIO text format uses a directory structure with .txt files for each data group.
This format is human-readable and doesn't require HDF5 library.

## Usage

To test VIAMD with these files:

```bash
# Build VIAMD with TREXIO support
cd /home/runner/work/viamd/viamd
cd ext/mdlib && git apply ../../docs/mdlib_trexio.patch
cd ../..
mkdir build && cd build
cmake -DVIAMD_ENABLE_TREXIO=ON ..
make

# Load a test file
./viamd ../test_data/h2_molecule.trexio
```

## Quick Validation

Run the validation script to verify everything is set up correctly:

```bash
cd test_data
./validate_trexio.sh
```

This will check:
- mdlib TREXIO files are in place
- Test data files exist and are valid
- Unit tests are present

## Data Validation

Each test file contains:
- Nucleus data (coordinates in Bohr, charges, labels)
- Electron configuration (up/down spin numbers)
- Metadata (package version, description)

Coordinates should be converted from Bohr to Angstrom when loaded in VIAMD.

## Test File Generation

To regenerate test files:

```bash
cd test_data
python3 create_test_trexio.py
```

## Unit Tests

Unit tests are located in `ext/mdlib/unittest/test_trexio.c` and include:
- `create_destroy` - Test TREXIO object creation/destruction
- `parse_h2_text` - Test parsing H2 molecule
- `parse_h2o_text` - Test parsing water molecule
- `system_init_h2` - Test md_system_t initialization
- `system_loader` - Test system loader interface

Run tests with:
```bash
cd build
ctest -V -R trexio
```
