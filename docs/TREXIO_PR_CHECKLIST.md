# TREXIO Support - PR Review Checklist

## Documentation Review

### User Documentation

- [x] **README.md**
  - [x] Accurate build instructions (no false claims about automatic download)
  - [x] Clear prerequisites section (TREXIO + HDF5)
  - [x] Multiple installation methods documented
  - [x] Status warning about current non-functional state
  - [x] Link to troubleshooting/known issues

- [x] **docs/TREXIO_SUPPORT.md**
  - [x] Complete user guide with all sections
  - [x] Prerequisites clearly stated
  - [x] Step-by-step build instructions
  - [x] Usage examples and workflows
  - [x] Comprehensive troubleshooting section
  - [x] Limitations and known issues documented
  - [x] Status warning at top of document

- [x] **docs/TREXIO_KNOWN_ISSUES.md** (NEW)
  - [x] Lists all known compilation/runtime issues
  - [x] Explains API compatibility problems
  - [x] Provides fix recommendations
  - [x] Testing checklist for future work
  - [x] References to related documentation

### Developer Documentation

- [x] **docs/TREXIO_SUPPORT.md - Contributing Section**
  - [x] Code structure explained
  - [x] How to extend TREXIO support
  - [x] Patch maintenance instructions
  - [x] Build system integration details
  - [x] Testing procedures

- [x] **Existing Implementation Docs**
  - [x] TREXIO_IMPLEMENTATION_STATUS.md (from Phase 1-4)
  - [x] PHASE4_SUMMARY.md (testing phase)
  - [x] PHASE5_COMPONENT.md (UI component)
  - [x] test_data/TESTING_GUIDE.md
  - [x] test_data/VALIDATION_RECIPE.md

### Build System Documentation

- [x] **scripts/apply_mdlib_trexio_patch.sh**
  - [x] Uses correct patch file (mdlib_trexio_original.patch)
  - [x] Idempotent (can run multiple times)
  - [x] Good error messages
  - [x] Troubleshooting hints

- [x] **CMake Integration**
  - [x] FindTREXIO.cmake module created
  - [x] Proper find_package usage
  - [x] Clear error messages when TREXIO not found
  - [x] PKG_CONFIG_PATH usage documented

## Code Review

### Build System

- [x] **CMakeLists.txt (root)**
  - [x] VIAMD_ENABLE_TREXIO option exists
  - [x] Properly propagates to mdlib

- [ ] **ext/mdlib/CMakeLists.txt** (via patch)
  - [ ] MD_ENABLE_TREXIO option added
  - [ ] find_package(TREXIO) called correctly
  - [ ] Proper error handling
  - [ ] Links TREXIO library
  - [ ] Adds MD_TREXIO definition
  - **Status:** ⚠️ Not tested - compilation fails before this

- [x] **ext/mdlib/cmake/FindTREXIO.cmake** (NEW)
  - [x] Tries pkg-config first
  - [x] Falls back to find_path/find_library
  - [x] Sets all required variables
  - [x] Creates imported target
  - **Status:** ✅ Created but not fully tested

### Source Code

- [ ] **ext/mdlib/src/md_trexio.c** (via patch)
  - [ ] Implements md_system_loader_i interface
  - [ ] Reads nucleus data correctly
  - [ ] Handles errors gracefully
  - [ ] Memory management is correct
  - **Status:** ❌ Does not compile - API compatibility issues

- [ ] **ext/mdlib/src/md_trexio.h** (via patch)
  - [ ] Public API is clean
  - [ ] Well documented
  - [ ] Follows mdlib patterns
  - **Status:** ⚠️ Header OK, implementation has issues

- [x] **src/components/trexio/trexio.cpp**
  - [x] Basic structure in place
  - [x] UI windows defined
  - [ ] Fully functional (blocked by loader issues)
  - **Status:** ⚠️ Stub implementation, needs working loader

- [x] **src/loader.cpp**
  - [x] TREXIO loader registered
  - [x] File extension detection
  - [x] Conditional compilation with #if MD_TREXIO
  - **Status:** ✅ Integration code looks correct

### Test Data

- [x] **test_data/**
  - [x] Sample .trexio files exist (H2, H2O, CH4)
  - [x] File creation scripts present
  - [x] Validation scripts available
  - [x] Documentation complete
  - **Status:** ✅ Test infrastructure ready

## Testing

### Build Tests

- [x] **Without TREXIO**
  - [x] VIAMD builds successfully
  - [x] No TREXIO-related errors
  - [x] Binary runs correctly
  - **Tested:** ✅ Ubuntu 22.04 (in CI)

- [ ] **With TREXIO** 
  - [ ] CMake configures successfully
  - [ ] All files compile without errors
  - [ ] Linking succeeds
  - [ ] Binary runs
  - **Status:** ❌ FAILS - compilation errors in md_trexio.c

### Functionality Tests

- [ ] **File Loading**
  - [ ] Can load .trexio text format
  - [ ] Can load .h5 HDF5 format
  - [ ] Correct atom counts
  - [ ] Correct coordinates
  - [ ] Error handling works
  - **Status:** ⏸️ Cannot test - build fails

- [ ] **UI Components**
  - [ ] TREXIO summary window displays
  - [ ] Shows correct file information
  - [ ] No crashes
  - **Status:** ⏸️ Cannot test - build fails

### Platform Tests

- [ ] **Linux**
  - [ ] Ubuntu 22.04 - Build + Run
  - [ ] Ubuntu 24.04 - Build + Run
  - **Status:** ⏸️ Build fails on all Linux

- [ ] **macOS**
  - [ ] Build successful
  - [ ] Run successful
  - **Status:** ⏸️ Not tested

- [ ] **Windows**
  - [ ] Build successful (MSVC)
  - [ ] Run successful
  - **Status:** ⏸️ Not tested

## CI/CD

- [ ] **GitHub Actions Workflows**
  - [ ] Test TREXIO build on Ubuntu 22.04
  - [ ] Test TREXIO build on Ubuntu 24.04
  - [ ] Test TREXIO build on macOS
  - [ ] Test TREXIO build on Windows
  - **Status:** ⚠️ Should be added when build works

## Documentation Completeness

### For Users

- [x] How to install TREXIO prerequisites
- [x] How to build VIAMD with TREXIO support
- [x] What file formats are supported
- [x] How to use TREXIO features
- [x] What to do when things go wrong
- [x] What limitations exist
- [x] Current status and known issues

### For Developers

- [x] How TREXIO integration is structured
- [x] How to extend TREXIO support
- [x] How to maintain the mdlib patch
- [x] How to test changes
- [x] What known issues exist
- [x] What needs to be fixed
- [x] How to regenerate patches

### For Reviewers

- [x] What has been implemented
- [x] What works and what doesn't
- [x] What needs to be fixed before merge
- [x] Testing procedures
- [x] Review checklist (this document)

## Acceptance Criteria (from Issue)

Based on the original Phase 5 issue requirements:

### Documentation

- [x] **README/developer docs updated**
  - [x] How to enable TREXIO
  - [x] Prerequisites clearly stated
  - [x] Build instructions accurate
  - [x] Developer guide for maintaining

- [x] **Loader limitations documented**
  - [x] What data is read vs ignored
  - [x] What file formats work
  - [x] Performance considerations
  - [x] Known issues and workarounds

- [x] **Expected user experience/errors documented**
  - [x] What users should see when loading files
  - [x] Error messages explained
  - [x] Troubleshooting guide
  - [x] When to expect issues

### PR Readiness

- [ ] **Builds on Linux** ❌ FAILS
- [ ] **Builds on macOS** ⏸️ NOT TESTED
- [ ] **Builds on Windows** ⏸️ NOT TESTED  
- [ ] **Tests pass** ⏸️ CANNOT RUN
- [ ] **Log output reasonable** ⏸️ CANNOT TEST
- [x] **New code well-commented** ✅ (existing code)

### Review Feedback

- [ ] **PR submitted and reviewed**
  - Status: Documentation is review-ready, code is not
- [ ] **Feedback addressed**
  - Status: Waiting for initial review

## Blocking Issues

### Must Fix Before Merge

1. **md_trexio.c compilation errors**
   - Priority: CRITICAL
   - Impact: Cannot build with TREXIO support
   - See: docs/TREXIO_KNOWN_ISSUES.md

2. **API compatibility with current mdlib**
   - Priority: CRITICAL
   - Impact: Code uses outdated mdlib structures
   - Requires: Investigation of current mdlib API

### Should Fix Before Merge

1. **int32_t vs int64_t warnings**
   - Priority: HIGH
   - Impact: May cause issues on some platforms
   - Simple to fix with temporary variables

2. **Complete testing**
   - Priority: HIGH
   - Impact: Unknown if loader actually works
   - Blocked by compilation issues

### Nice to Have

1. **Molecular orbital visualization**
   - Priority: MEDIUM
   - Status: UI stub exists
   - Can be done in future PR

2. **Extended data group support**
   - Priority: LOW
   - Can be added incrementally

## Recommendations

### For Current PR

**Option 1: Documentation-Only PR (RECOMMENDED)**
- Merge current PR with documentation updates only
- Create follow-up issue for fixing compilation
- Benefits:
  - Documentation improvements are valuable independently
  - Sets accurate expectations for users
  - Provides clear roadmap for future work

**Option 2: Fix and Complete**
- Fix md_trexio.c API compatibility issues
- Complete testing
- Merge fully functional implementation
- Benefits:
  - Feature is complete and usable
- Drawbacks:
  - Requires more time and mdlib API investigation
  - May need multiple iterations

### For Follow-Up Work

1. Create issue: "Fix TREXIO loader API compatibility"
   - Use TREXIO_KNOWN_ISSUES.md as basis
   - Assign to developer familiar with mdlib

2. Create issue: "Complete TREXIO molecular orbital visualization"
   - After loader works
   - Implement UI features

3. Add TREXIO builds to CI when working
   - Start with one platform
   - Expand to all platforms

## Sign-Off Checklist

Before requesting final review:

- [x] All documentation reviewed and accurate
- [x] Build instructions tested and work (for disabled TREXIO)
- [x] Known issues documented clearly
- [x] Status warnings prominent in all docs
- [ ] Code compiles (❌ BLOCKED)
- [ ] Tests pass (⏸️ BLOCKED)
- [x] PR description complete
- [x] Reviewer guidance provided

## Notes for Reviewer

This PR primarily contains **documentation improvements** for TREXIO support. The actual TREXIO loader implementation has API compatibility issues that prevent compilation.

**What to focus on:**
1. Documentation accuracy and completeness
2. User-facing instructions and troubleshooting
3. Developer guidance for future fixes
4. Whether the status warnings are prominent enough

**What NOT to focus on:**
1. The non-working loader code (already documented in known issues)
2. Missing features (documented in limitations)
3. Lack of CI integration (premature until code works)

**Recommendation:** Approve this PR for the documentation improvements, then address code issues in a follow-up PR/issue.
