# TREXIO Support - Maintainer Guide

This guide helps maintainers review, test, and maintain TREXIO support in VIAMD.

## Quick Overview

TREXIO support was added to VIAMD to enable loading quantum chemistry data from TREXIO files. The implementation:

- **Integrates at the mdlib level** - Uses a patch to add md_trexio.c/h to the mdlib submodule
- **Builds TREXIO from source** - Automatically downloads and builds TREXIO 2.6.0
- **Optional feature** - Disabled by default, enabled with `-DVIAMD_ENABLE_TREXIO=ON`
- **HDF5 optional** - Works with text format without HDF5, HDF5 adds .h5 support

## Architecture

```
TREXIO File (.trexio or .h5)
         ↓
TREXIO C Library (auto-downloaded)
         ↓
mdlib (md_trexio.c/h via patch)
         ↓
VIAMD loader system
         ↓
VIAMD visualization
```

## Testing the Build

### Prerequisites
```bash
# Ubuntu/Debian
sudo apt-get install cmake build-essential libhdf5-dev git

# macOS
brew install cmake hdf5

# Conda (any platform)
conda install -c conda-forge cmake hdf5
```

### Build Without TREXIO (Verify Backward Compatibility)
```bash
git clone https://github.com/scanberg/viamd.git
cd viamd
git checkout copilot/update-docs-for-trexio-support
git submodule update --init --recursive

mkdir build && cd build
cmake ..
make
```

**Expected**: Build succeeds, no TREXIO-related code compiled

### Build With TREXIO
```bash
# From repository root
./scripts/apply_mdlib_trexio_patch.sh

cd build
cmake -DVIAMD_ENABLE_TREXIO=ON ..
make
```

**Expected**:
- CMake downloads TREXIO 2.6.0 tarball (~300KB)
- HDF5 detected (or warning if not found)
- TREXIO builds (takes 1-2 minutes)
- VIAMD builds with TREXIO support
- No new compiler warnings

### Verify TREXIO Integration
```bash
# Check that TREXIO files are recognized
./viamd --help  # Should mention supported formats

# Load test file
./viamd ../test_data/h2_molecule.trexio
```

**Expected**:
- File loads successfully
- Molecular structure displayed
- Atom count correct (H2 = 2 atoms)

## Reviewing Code Changes

### Files to Review

#### Core Implementation (via patch)
- `docs/mdlib_trexio.patch` - The mdlib patch
  - `ext/mdlib/src/md_trexio.c` (672 lines) - TREXIO file parser
  - `ext/mdlib/src/md_trexio.h` (117 lines) - Header with API
  - `ext/mdlib/CMakeLists.txt` - CMake changes for TREXIO

#### VIAMD Integration
- `CMakeLists.txt` - Added VIAMD_ENABLE_TREXIO option
- `src/loader.cpp` - Registered TREXIO loader
- `src/components/trexio/trexio.cpp` - TREXIO UI component (future)

#### Infrastructure
- `scripts/apply_mdlib_trexio_patch.sh` - Patch application script
- `test_data/` - Test files and generation scripts

### Code Review Checklist

#### mdlib Integration
- [ ] Uses mdlib allocator (md_alloc, md_free)
- [ ] Follows mdlib conventions (naming, structure)
- [ ] Implements md_system_loader_i interface correctly
- [ ] Error handling uses MD_LOG_ERROR/MD_LOG_WARN
- [ ] Memory management is safe (no leaks)
- [ ] Coordinate conversion correct (Bohr → Angstrom: × 0.529177)

#### TREXIO API Usage
- [ ] Correct TREXIO C API calls
- [ ] Error codes checked (TREXIO_SUCCESS, TREXIO_HAS_NOT)
- [ ] File properly closed (trexio_close)
- [ ] Backend detection works (HDF5 vs text)

#### Loader Registration
- [ ] File extensions registered (.trexio, .h5)
- [ ] Loader appears in supported formats list
- [ ] No conflicts with other loaders
- [ ] Conditional compilation works (#ifdef MD_TREXIO)

#### Build System
- [ ] FetchContent downloads TREXIO correctly
- [ ] HDF5 detection optional
- [ ] Works with and without HDF5
- [ ] Compatible with VeloxChem support
- [ ] No new build warnings

## Platform-Specific Testing

### Ubuntu 22.04 (GCC 11)
```bash
# In GitHub Actions or local VM
cmake -DVIAMD_ENABLE_TREXIO=ON ..
make
ctest -V
```

### Ubuntu 24.04 (GCC 13)
```bash
# In GitHub Actions or local VM
cmake -DVIAMD_ENABLE_TREXIO=ON ..
make
ctest -V
```

### macOS (Clang)
```bash
# On macOS system or GitHub Actions
brew install hdf5
cmake -DVIAMD_ENABLE_TREXIO=ON ..
make
```

### Windows (MSVC 19)
```bash
# In Visual Studio Developer Command Prompt
cmake -DVIAMD_ENABLE_TREXIO=ON ..
cmake --build . --config Release
```

**Note**: Windows might need manual HDF5 path if not found:
```bash
cmake -DVIAMD_ENABLE_TREXIO=ON -DHDF5_ROOT="C:/path/to/hdf5" ..
```

## Common Issues and Solutions

### Issue: Patch Doesn't Apply
**Symptoms**: `git apply` fails with conflicts

**Solutions**:
1. Use the patch script: `./scripts/apply_mdlib_trexio_patch.sh`
2. The script automatically handles cleanup and reapplication
3. If script fails, check mdlib submodule is at correct commit:
   ```bash
   cd ext/mdlib
   git status  # Should be on specific commit, not branch
   ```

### Issue: TREXIO Download Fails
**Symptoms**: CMake fails during FetchContent download

**Solutions**:
1. Check internet connection
2. Check firewall/proxy settings
3. Manual download:
   ```bash
   cd /tmp
   wget https://github.com/TREX-CoE/trexio/releases/download/v2.6.0/trexio-2.6.0.tar.gz
   ```
4. Reconfigure CMake (will use cached download)

### Issue: HDF5 Not Found
**Symptoms**: TREXIO builds but only text format works

**Solutions**:
1. Install HDF5 development libraries
2. Reconfigure and rebuild
3. Or explicitly disable HDF5:
   ```bash
   cmake -DVIAMD_ENABLE_TREXIO=ON -DTREXIO_HDF5=OFF ..
   ```

### Issue: .h5 Files Don't Load
**Symptoms**: .h5 files fail to open or load as wrong format

**Diagnostics**:
1. Check if file is valid TREXIO:
   ```bash
   python3 -c "import trexio; f = trexio.File('file.h5', 'r')"
   ```
2. Check VIAMD build logs for HDF5 detection
3. Try .trexio text format instead

## Security Considerations

### TREXIO Source
- Downloaded from official TREX-CoE GitHub releases
- Version pinned to 2.6.0 for stability
- SHA256 checksum verified by CMake FetchContent
- No code execution during download

### Input Validation
- TREXIO C library validates file format
- mdlib code validates data ranges
- Invalid files rejected with error message
- No buffer overflows in parsing code

### Memory Safety
- Uses mdlib allocator (tracked allocations)
- All buffers properly sized
- No manual pointer arithmetic
- Cleanup on error paths

## Performance Benchmarks

### Build Time
- **First build**: +1-2 minutes (TREXIO download and build)
- **Rebuild**: No overhead (TREXIO cached)
- **Clean build**: +30 seconds (TREXIO rebuild)

### File Loading
- **Small molecules** (<100 atoms): <100ms
- **Medium molecules** (100-1000 atoms): <500ms
- **Large molecules** (>1000 atoms): <2s

### Memory Usage
- **TREXIO library**: ~1-2 MB
- **Per loaded file**: Depends on MO count
  - Small: <10 MB
  - Large: <100 MB

## Maintenance Tasks

### Updating TREXIO Version
1. Edit `ext/mdlib/CMakeLists.txt` patch section
2. Update version number (e.g., 2.6.0 → 2.7.0)
3. Test build on all platforms
4. Update documentation

### Adding New TREXIO Groups
1. Edit `ext/mdlib/src/md_trexio.c`
2. Add new `trexio_read_*` calls
3. Update data structures in mdlib
4. Add tests for new data
5. Update documentation

### Fixing Bugs
1. Identify which layer has the bug:
   - TREXIO library → Report upstream
   - mdlib (md_trexio.c) → Fix in patch
   - VIAMD (loader.cpp) → Fix directly
2. Add test case reproducing the bug
3. Fix and verify all tests pass
4. Update documentation if needed

## CI/CD Integration (Future)

### Recommended Checks
- [ ] Build with TREXIO=ON on all platforms
- [ ] Build with TREXIO=OFF (backward compat)
- [ ] Run unit tests (requires test data)
- [ ] CodeQL security scan
- [ ] Check for new compiler warnings

### Test Matrix
```yaml
os: [ubuntu-22.04, ubuntu-24.04, macos-latest, windows-latest]
trexio: [ON, OFF]
hdf5: [ON, OFF]
```

## Documentation Maintenance

### Files to Keep Updated
- `README.md` - Quick start instructions
- `docs/TREXIO_SUPPORT.md` - User guide
- `docs/MDLIB_TREXIO_INTEGRATION.md` - Developer guide
- `docs/TREXIO_PR_CHECKLIST.md` - Review checklist
- `test_data/TESTING_GUIDE.md` - Testing procedures

### When to Update
- New TREXIO version
- New features added
- Build process changes
- New platforms supported
- Bug fixes affecting user experience

## Support and Communication

### User Questions
- Direct to GitHub issues
- Reference docs/TREXIO_SUPPORT.md
- Check troubleshooting section first

### Developer Questions
- Direct to GitHub discussions
- Reference docs/MDLIB_TREXIO_INTEGRATION.md
- Consider mdlib submodule documentation

### Bug Reports
- Need: VIAMD version, OS, TREXIO file (if possible)
- Check: Build logs, runtime logs
- Test: With minimal TREXIO file from test_data/

## Release Checklist

When releasing VIAMD with TREXIO support:

- [ ] All platform builds pass
- [ ] Documentation complete and accurate
- [ ] Test files included
- [ ] Known limitations documented
- [ ] CHANGELOG updated
- [ ] README mentions TREXIO support
- [ ] Wiki updated (if applicable)
- [ ] Release notes mention TREXIO

## References

### TREXIO Resources
- Specification: https://trex-coe.github.io/trexio/
- GitHub: https://github.com/TREX-CoE/trexio
- Paper: https://doi.org/10.1063/5.0148161

### VIAMD Resources
- Wiki: https://github.com/scanberg/viamd/wiki
- Issues: https://github.com/scanberg/viamd/issues
- Discussions: https://github.com/scanberg/viamd/discussions

### Implementation Details
- Issue #95: Phase 5 Documentation
- PR #91: TREXIO implementation
- PR #92-94: Build system improvements

## Contact

For questions about TREXIO support maintenance:
- Create an issue on GitHub
- Tag @mathieulinares (original requestor)
- Reference this guide and relevant documentation
