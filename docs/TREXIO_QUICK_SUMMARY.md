# TREXIO Support - Quick Summary

**Status**: Ready for Integration  
**Date**: 2025-11-17

## TL;DR

- ❌ TREXIO symbols are **NOT** in mdlib upstream
- ✅ Patch adds TREXIO support to mdlib (md_trexio.c/h)
- ✅ TREXIO library auto-downloaded via CMake FetchContent
- ✅ 3 test files ready (H2, H2O, CH4)
- ✅ Complete documentation exists

## Quick Integration

```bash
# 1. Apply patch
./scripts/apply_mdlib_trexio_patch.sh

# 2. Build with TREXIO
mkdir build && cd build
cmake -DVIAMD_ENABLE_TREXIO=ON ..
make

# 3. Test
./viamd ../test_data/h2_molecule.trexio
```

## How TREXIO Symbols are Provided

### NOT Internal to mdlib
- Upstream mdlib at `ext/mdlib` does NOT contain md_trexio.*
- No TREXIO code in scanberg/mdlib repository

### Integration Method: Patch + External Library

**Step 1: Patch adds integration code**
- Patch file: `docs/mdlib_trexio.patch`
- Adds: `src/md_trexio.c` (672 lines)
- Adds: `src/md_trexio.h` (117 lines)
- Modifies: `CMakeLists.txt` (TREXIO build config)

**Step 2: CMake downloads TREXIO library**
- Uses `FetchContent` to download TREXIO v2.6.0 from GitHub releases
- URL: https://github.com/TREX-CoE/trexio/releases/download/v2.6.0/trexio-2.6.0.tar.gz
- Builds TREXIO as static library
- Links into mdlib
- No system TREXIO installation required

**Step 3: Symbol resolution**
```
viamd → links mdlib → includes md_trexio.c → calls TREXIO C library
```

## Test Files Available

All in `test_data/` directory, TREXIO text format:

| File | Atoms | Description |
|------|-------|-------------|
| h2_molecule.trexio | 2 (H, H) | Simple diatomic |
| h2o_molecule.trexio | 3 (O, H, H) | Water molecule |
| ch4_molecule.trexio | 5 (C, 4×H) | Methane |

**Data included:**
- Atomic coordinates (Bohr → Ångström conversion)
- Nuclear charges
- Element labels
- Electron counts (up/down spin)

## Documentation

| File | Purpose |
|------|---------|
| `docs/TREXIO_AUDIT_AND_INTEGRATION_PLAN.md` | **Complete audit** (450+ lines) |
| `docs/TREXIO_SUPPORT.md` | User guide |
| `docs/MDLIB_TREXIO_INTEGRATION.md` | Integration details |
| `test_data/TESTING_GUIDE.md` | Testing instructions |
| `test_data/README.md` | Quick start |

## Build Options

```cmake
# Enable TREXIO support
cmake -DVIAMD_ENABLE_TREXIO=ON ..

# With HDF5 backend (recommended for .h5 files)
cmake -DVIAMD_ENABLE_TREXIO=ON -DENABLE_HDF5=ON ..

# Combined with VeloxChem
cmake -DVIAMD_ENABLE_TREXIO=ON -DVIAMD_ENABLE_VELOXCHEM=ON ..
```

## File Format Support

| Extension | Backend | HDF5 Required | Status |
|-----------|---------|---------------|--------|
| `.trexio` | TEXT | No | ✅ Ready |
| `.h5` | HDF5 | Yes | ✅ Ready (with HDF5) |

## Implementation Status

### ✅ Complete
- Core file parsing
- Nucleus data extraction
- System initialization
- Loader registration
- Build integration
- Test infrastructure

### ⏸️ Incomplete
- GTO extraction (placeholder)
- MO visualization
- UI component
- Unit tests in mdlib

## Next Steps

1. **For basic molecule loading**: Just apply patch and build
2. **For orbital visualization**: Implement GTO extraction functions
3. **For full features**: Complete UI component in `src/components/trexio/`

## Dependencies

**Required:**
- CMake ≥ 3.20
- C compiler
- Internet (for TREXIO download)

**Optional:**
- HDF5 library (for .h5 file support)
- Python + trexio (for generating test files)

## Support

- Full audit: `docs/TREXIO_AUDIT_AND_INTEGRATION_PLAN.md`
- Issues: https://github.com/scanberg/viamd/issues
- TREXIO docs: https://trex-coe.github.io/trexio/

---

**Ready to integrate**: ✅ All prerequisites met, documentation complete, test files available.
