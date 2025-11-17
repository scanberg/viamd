# TREXIO Testing Guide

## Phase 4: Testing and Validation

This guide covers comprehensive testing of VIAMD's TREXIO support.

## Prerequisites

### Required
- C/C++ compiler (gcc, clang, or MSVC)
- CMake (>= 3.20)
- Make or Ninja

### Optional but Recommended
- TREXIO library (>= 2.0.0) - for HDF5 backend
- HDF5 library (>= 1.8) - for HDF5 backend
- PySCF - for generating test files
- Python 3 with TREXIO bindings

## Quick Start

### 1. Validate Setup

```bash
cd test_data
./validate_trexio.sh
```

This checks that all necessary files are in place.

### 2. Build VIAMD with TREXIO Support

**Option A: Automated Build Script**
```bash
cd test_data
./build_and_test.sh
```

**Option B: Manual Build**
```bash
# Apply mdlib patch (if not already done)
cd ext/mdlib
git apply ../../docs/mdlib_trexio.patch

# Configure and build
cd ../..
mkdir build && cd build
cmake -DVIAMD_ENABLE_TREXIO=ON ..
make -j$(nproc)

# Run tests
ctest -V -R trexio
```

### 3. Test with Sample Files

```bash
# From build directory
./viamd ../test_data/h2_molecule.trexio
./viamd ../test_data/h2o_molecule.trexio
./viamd ../test_data/ch4_molecule.trexio
```

## Test Data

### Included Test Files (Text Format)

The repository includes three pre-generated TREXIO text format files:

1. **h2_molecule.trexio** - Hydrogen molecule (H2)
   - 2 atoms
   - 2 electrons (1 up, 1 down)
   - Simple diatomic for basic testing

2. **h2o_molecule.trexio** - Water molecule (H2O)
   - 3 atoms (O, H, H)
   - 10 electrons (5 up, 5 down)
   - Tests bent molecular geometry

3. **ch4_molecule.trexio** - Methane molecule (CH4)
   - 5 atoms (C, 4×H)
   - 10 electrons (5 up, 5 down)
   - Tests tetrahedral geometry

These files use TREXIO text format and don't require HDF5.

### Generating Files with PySCF

For more realistic test cases with full quantum chemistry data:

```bash
# Install dependencies
pip install pyscf trexio

# Generate test files
cd test_data
python3 create_pyscf_trexio.py
```

This creates:
- `h2_pyscf.h5` - H2 with SCF calculation data
- `h2o_pyscf.h5` - H2O with SCF calculation data

These files include:
- Molecular orbital coefficients
- MO energies
- Basis set information
- SCF energies

## Unit Tests

### Running Unit Tests

```bash
cd build
ctest -V -R trexio
```

### Test Coverage

The unit tests (`ext/mdlib/unittest/test_trexio.c`) include:

1. **create_destroy** - TREXIO object lifecycle
   - Tests object creation and destruction
   - Verifies memory management

2. **parse_h2_text** - H2 molecule parsing
   - Validates atom count (2)
   - Checks atomic charges (both 1.0)
   - Verifies coordinates
   - Checks Bohr to Angstrom conversion
   - Validates electron configuration

3. **parse_h2o_text** - Water molecule parsing
   - Validates 3 atoms
   - Checks charges (O=8, H=1, H=1)
   - Verifies 10 electrons

4. **system_init_h2** - System initialization
   - Tests md_system_t creation from TREXIO
   - Validates coordinate conversion to float
   - Checks atom count consistency

5. **system_loader** - Loader interface
   - Verifies loader function pointer
   - Tests init function availability

6. **Conditional compilation test**
   - Ensures tests pass when TREXIO is disabled

## Testing Different Backends

### Text Backend (No HDF5)

```bash
cmake -DVIAMD_ENABLE_TREXIO=ON -DENABLE_HDF5=OFF ..
make
./viamd ../test_data/h2_molecule.trexio
```

**Advantages:**
- No HDF5 dependency
- Human-readable files
- Easy debugging

**Limitations:**
- Slower I/O
- Larger file sizes
- Text precision limits

### HDF5 Backend

```bash
cmake -DVIAMD_ENABLE_TREXIO=ON -DENABLE_HDF5=ON ..
make
./viamd ../test_data/h2_pyscf.h5
```

**Advantages:**
- Fast I/O
- Compact storage
- Full precision

**Limitations:**
- Requires HDF5 library
- Binary format (not human-readable)

## Testing with Real Quantum Chemistry Files

### From PySCF

```python
from pyscf import gto, scf
import trexio

mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='6-31g')
mf = scf.RHF(mol).run()

with trexio.File('molecule.h5', 'w', trexio.TREXIO_HDF5) as f:
    trexio.write_nucleus_num(f, mol.natm)
    trexio.write_nucleus_charge(f, mol.atom_charges())
    trexio.write_nucleus_coord(f, mol.atom_coords())
    # ... write more data
```

### From Quantum Package

Quantum Package can export TREXIO files directly:
```bash
qp_run export_trexio molecule.trexio
```

### From Other Codes

Many quantum chemistry codes support TREXIO:
- CP2K
- FHI-aims
- Dirac
- TurboRVB

## Validation Checklist

- [ ] Build succeeds with TREXIO enabled
- [ ] Build succeeds with TREXIO disabled (stubs work)
- [ ] Unit tests pass
- [ ] Can load .trexio text format files
- [ ] Can load .h5 HDF5 format files (if HDF5 enabled)
- [ ] Coordinates display correctly (Bohr → Angstrom)
- [ ] Atomic charges are correct
- [ ] Atom labels display properly
- [ ] Electron configuration is accurate
- [ ] Multiple molecules load without errors
- [ ] No memory leaks (run with valgrind if available)

## Troubleshooting

### Build Fails with "TREXIO library not found"

**Solution:**
```bash
# Option 1: Install via conda
conda install -c conda-forge trexio

# Option 2: Build from source
git clone https://github.com/TREX-CoE/trexio.git
cd trexio
./configure --prefix=$HOME/.local
make && make install

# Option 3: Specify TREXIO location
cmake -DVIAMD_ENABLE_TREXIO=ON -DTREXIO_ROOT=/path/to/trexio ..
```

### Unit Tests Fail

**Check:**
1. Test data files exist in `test_data/`
2. Paths in test code are correct (relative paths from build dir)
3. TREXIO library is properly installed
4. Run with verbose output: `ctest -V -R trexio`

### Cannot Load .h5 Files

**Solution:**
Ensure HDF5 support is enabled:
```bash
cmake -DVIAMD_ENABLE_TREXIO=ON -DENABLE_HDF5=ON ..
```

### Coordinates Appear Incorrect

**Check:**
- TREXIO stores coordinates in Bohr
- md_trexio.c should convert to Angstrom
- Verify conversion factor: 1 Bohr = 0.529177 Angstrom

## Performance Testing

### File Loading Speed

```bash
time ./viamd ../test_data/h2_molecule.trexio
```

Compare text vs HDF5 backends.

### Memory Usage

```bash
# Linux
/usr/bin/time -v ./viamd molecule.trexio

# macOS
/usr/bin/time -l ./viamd molecule.trexio
```

## Continuous Integration

For CI/CD pipelines:

```yaml
# .github/workflows/test_trexio.yml
name: TREXIO Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
        with:
          submodules: recursive
      
      - name: Install TREXIO
        run: |
          conda install -c conda-forge trexio
      
      - name: Build and Test
        run: |
          cd test_data
          ./build_and_test.sh
```

## Advanced Testing

### Fuzzing Test Files

Create malformed TREXIO files to test error handling:

```python
# Create invalid TREXIO file
import os
os.makedirs('invalid.trexio', exist_ok=True)
with open('invalid.trexio/nucleus.txt', 'w') as f:
    f.write('num\n')
    f.write('invalid_number\n')  # Should be integer
```

Test that VIAMD handles gracefully:
```bash
./viamd invalid.trexio  # Should show error, not crash
```

### Large Molecule Tests

Test scalability with larger systems:
- 50+ atoms
- 1000+ basis functions
- Multiple MB file size

### Concurrent Loading

Test thread safety (if applicable):
```bash
for i in {1..10}; do
    ./viamd molecule.trexio &
done
wait
```

## Reporting Issues

When reporting TREXIO-related issues, include:

1. VIAMD version and commit
2. TREXIO library version
3. Build configuration (CMake output)
4. Sample TREXIO file (if possible)
5. Error messages and logs
6. Operating system and compiler

## References

- [TREXIO Documentation](https://trex-coe.github.io/trexio/)
- [TREXIO GitHub](https://github.com/TREX-CoE/trexio)
- [VIAMD Documentation](docs/TREXIO_SUPPORT.md)
- [PySCF Documentation](https://pyscf.org/)

## Summary

Phase 4 testing ensures VIAMD's TREXIO support is:
- ✅ Correctly implemented
- ✅ Well tested
- ✅ Compatible with quantum chemistry codes
- ✅ Robust and reliable
- ✅ Ready for production use
