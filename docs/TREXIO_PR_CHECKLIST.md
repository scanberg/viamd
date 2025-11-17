# TREXIO Support PR Review Checklist

This document provides a comprehensive checklist for reviewing TREXIO support in VIAMD, addressing the requirements from Phase 5 (Issue #95).

## Documentation Review

### User Documentation
- [x] **README.md updated** with TREXIO build instructions
  - Prerequisites clearly stated (HDF5 optional)
  - Build steps documented
  - CMake options explained (`-DVIAMD_ENABLE_TREXIO=ON`)
- [x] **docs/TREXIO_SUPPORT.md** provides comprehensive user guide
  - What is TREXIO and why use it
  - Supported features and limitations
  - Build instructions with examples
  - Usage examples and workflows
  - Troubleshooting section
  - File format details
- [x] **Loader limitations documented**
  - Read-only access (no TREXIO writing)
  - Limited to nucleus, basis, MO, electron groups
  - No trajectory support
  - Basic visualization only (advanced features in progress)
- [x] **Error handling documented**
  - TREXIO download failures
  - HDF5 backend unavailable
  - File format detection issues

### Developer Documentation
- [x] **docs/MDLIB_TREXIO_INTEGRATION.md** explains mdlib integration
- [x] **docs/IMPLEMENTATION_SUMMARY.md** provides technical overview
- [x] **docs/PHASE5_COMPONENT.md** documents component architecture
- [x] **docs/mdlib_trexio.patch** available for applying mdlib changes
- [x] **scripts/README.md** documents helper scripts
- [x] **Code comments** in key areas:
  - src/components/trexio/trexio.cpp
  - src/loader.cpp (TREXIO registration)
  - mdlib patch includes comments

### Testing Documentation
- [x] **test_data/TESTING_GUIDE.md** explains testing process
- [x] **test_data/README.md** documents test files
- [x] **docs/PHASE4_TESTING_REPORT.md** summarizes test results

## Build System Review

### CMake Configuration
- [x] **VIAMD_ENABLE_TREXIO** option added to CMakeLists.txt
- [x] **Default OFF** - feature is optional
- [x] **Automatic TREXIO download** - uses FetchContent
- [x] **HDF5 detection** - optional dependency handled gracefully
- [x] **Conditional compilation** - uses `#ifdef MD_TREXIO`
- [x] **Compatible with VeloxChem** - both can be enabled simultaneously

### Patch Application
- [x] **scripts/apply_mdlib_trexio_patch.sh** - idempotent script
- [x] **Manual patch instructions** documented as fallback
- [x] **No trailing whitespace** in patch file
- [x] **Patch applies cleanly** to mdlib submodule

## Platform Testing

### Linux
- [ ] **Ubuntu 22.04** - Build succeeds with TREXIO=ON
  - [ ] Test with GCC 11
  - [ ] HDF5 detection works
  - [ ] TREXIO download succeeds
  - [ ] Sample files load correctly
- [ ] **Ubuntu 24.04** - Build succeeds with TREXIO=ON
  - [ ] Test with GCC 13
  - [ ] No build warnings introduced

### macOS
- [ ] **macOS (Clang)** - Build succeeds with TREXIO=ON
  - [ ] HDF5 from Homebrew detected
  - [ ] No compilation errors
  - [ ] Sample files load correctly

### Windows
- [ ] **Windows (MSVC 19)** - Build succeeds with TREXIO=ON
  - [ ] HDF5 detection works (if available)
  - [ ] Text-only backend works without HDF5
  - [ ] Sample files load correctly

## Functional Testing

### File Loading
- [ ] **.trexio (text format)** files load correctly
  - [ ] Atomic coordinates displayed
  - [ ] Element types correct
  - [ ] Molecule structure visualized
- [ ] **.h5 (HDF5 format)** files load correctly
  - [ ] Requires HDF5 backend
  - [ ] Falls back to text format if HDF5 unavailable
- [ ] **Error handling** works as expected
  - [ ] Invalid files rejected gracefully
  - [ ] Missing data handled appropriately
  - [ ] User-friendly error messages

### Test Files
- [ ] **h2_molecule.trexio** - Simple H2 loads
- [ ] **h2o_molecule.trexio** - Water molecule loads
- [ ] **ch4_molecule.trexio** - Methane loads
- [ ] **PySCF-generated files** - Real quantum chemistry data loads

### Integration
- [ ] **Loader registration** - TREXIO appears in supported formats
- [ ] **File extension detection** - .trexio and .h5 recognized
- [ ] **No conflicts** with other loaders (VeloxChem, PDB, etc.)
- [ ] **Component UI** - TREXIO component accessible if built

## Code Quality

### Implementation
- [x] **mdlib integration** - md_trexio.c/h follows mdlib patterns
- [x] **Memory management** - Uses mdlib allocator
- [x] **Error handling** - Comprehensive logging
- [x] **Coordinate conversion** - Bohr to Angstrom with correct constants
- [x] **API usage** - Correct TREXIO C API calls

### Code Review
- [x] **Follows existing patterns** - Similar to VeloxChem integration
- [x] **Minimal changes** - Only necessary modifications
- [x] **Backward compatible** - No breaking changes
- [x] **Well-commented** - Key sections explained
- [x] **No memory leaks** - Proper cleanup

### Security
- [ ] **CodeQL scan** - No new vulnerabilities
- [ ] **Input validation** - File data validated
- [ ] **Buffer safety** - No buffer overflows
- [ ] **Dependency security** - TREXIO from official source

## CI/CD Integration

### GitHub Actions
- [ ] **Builds pass** on all platforms (Linux, macOS, Windows)
- [ ] **No new warnings** introduced
- [ ] **Tests pass** (when CI tests added)
- [ ] **Code coverage** acceptable

### Performance
- [ ] **Build time** - TREXIO download and build time acceptable
- [ ] **File loading** - Performance adequate for typical files
- [ ] **Memory usage** - No excessive memory consumption

## User Experience

### Expected Behavior
1. **Enable TREXIO support**
   - Run patch script
   - Configure with `-DVIAMD_ENABLE_TREXIO=ON`
   - Build completes successfully
   - TREXIO download automatic

2. **Load TREXIO file**
   - Open File → Open or drag-and-drop
   - Select .trexio or .h5 file
   - Molecule structure displays
   - Atom count and types correct

3. **Error scenarios**
   - Invalid file: Clear error message
   - Missing HDF5: Graceful degradation to text format
   - Corrupted data: Error logged, user notified

### Log Output
- [ ] **Informative messages** during build
- [ ] **Helpful warnings** if HDF5 not found
- [ ] **Clear errors** for file loading failures
- [ ] **No spam** - Only relevant messages

## Acceptance Criteria Checklist

✅ = Complete | ⚠️ = Needs Testing | ❌ = Not Done

### From Issue #95:

1. **Update README/developer docs**
   - ✅ README.md updated with TREXIO section
   - ✅ How to enable TREXIO documented
   - ✅ Prerequisites listed (HDF5 optional)
   - ✅ Build instructions provided
   - ✅ CMake options explained

2. **Document loader limitations**
   - ✅ Read-only operation
   - ✅ Limited group support
   - ✅ No trajectory support
   - ✅ Current visualization limitations
   - ✅ Expected user experience described

3. **Document expected errors**
   - ✅ Download failures
   - ✅ HDF5 backend issues
   - ✅ File format detection
   - ✅ Invalid file handling
   - ✅ Troubleshooting guide provided

4. **Review PR: builds on Linux/Mac**
   - ⚠️ Ubuntu 22.04 build verification needed
   - ⚠️ Ubuntu 24.04 build verification needed
   - ⚠️ macOS build verification needed
   - ⚠️ Windows build verification needed

5. **Review PR: tests pass**
   - ⚠️ Unit tests need to be run (requires TREXIO C library)
   - ✅ Test files created and documented
   - ✅ Test generation scripts provided
   - ✅ Testing guide available

6. **Review PR: log output is reasonable**
   - ⚠️ Build log review needed
   - ⚠️ Runtime log review needed
   - ✅ Error messages documented

7. **Review PR: code well-commented**
   - ✅ mdlib patch includes comments
   - ✅ Loader registration commented
   - ✅ Component code documented
   - ✅ API boundaries explained
   - ✅ Error boundaries handled

## Remaining Work

### Critical (Required for merge)
- [ ] Verify builds on Ubuntu 22.04/24.04
- [ ] Verify builds on macOS
- [ ] Verify builds on Windows (if applicable)
- [ ] Run and document test results
- [ ] CodeQL security scan

### Nice to Have (Future improvements)
- [ ] Add CI/CD automation
- [ ] Wiki page with detailed examples
- [ ] Video tutorial
- [ ] More test files
- [ ] Performance benchmarks

## Reviewer Notes

### Things to Check
1. Does the build work without TREXIO support? (backward compatibility)
2. Does the build work with TREXIO support enabled?
3. Does HDF5 detection work correctly?
4. Are error messages clear and helpful?
5. Is the patch application idempotent?
6. Do test files load correctly?
7. Are there any memory leaks?
8. Is the code maintainable?

### Potential Issues
- TREXIO download might fail on restricted networks
- HDF5 detection might need manual configuration on some systems
- .h5 extension conflict with VeloxChem needs clear documentation
- First-time build might be slow due to TREXIO compilation

## Sign-off

- [ ] Documentation is clear and complete
- [ ] Build instructions are accurate
- [ ] Limitations are well-documented
- [ ] Error handling is comprehensive
- [ ] Code quality is acceptable
- [ ] Tests demonstrate functionality
- [ ] Ready for maintainer review

## References

- Issue #95: https://github.com/scanberg/viamd/issues/95
- PR #91: https://github.com/scanberg/viamd/pull/91
- TREXIO Docs: https://trex-coe.github.io/trexio/
- VIAMD Wiki: https://github.com/scanberg/viamd/wiki
