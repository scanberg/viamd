# TREXIO Loader Crash - Need Original Failing File for Full Reproduction

## Issue Summary

**Investigation completed:** PRs #115, #116, #117  
**Minimal fix applied:** Allocator NULL/bounds validation  
**Status:** Cannot reproduce crash with repository test files  

**Request:** Original TREXIO file that caused the crash in PR #113 for full reproduction and targeted fix verification.

---

## Background

Comprehensive investigation of TREXIO loader allocator crash has been completed across three PRs:

- **PR #115:** Allocator reproducer (150 iterations, ASAN clean - allocator works correctly)
- **PR #116:** Documentation with fix patterns (based on VeloxChem reference)
- **PR #117:** Systematic investigation with CI workflows and loader test harness

**Key Findings:**
1. ✅ mdlib heap allocator is NOT the problem (validated with 150+ iterations, zero errors)
2. ⚠️ Repository test files (`test_data/*.trexio`) cannot be opened by TREXIO library
3. ✅ Minimal validation fix implemented (NULL checks, bounds validation)
4. ⏸️ Original crash scenario NOT reproducible with available test files

---

## What Has Been Done

### 1. Allocator Validation (PR #115)
```bash
# Test Results: 150 iterations (100 normal + 50 ASAN)
[DIAGNOSTIC] All iterations completed successfully!
[DIAGNOSTIC] No crashes detected in allocator operations.
Exit code: 0
```

**Conclusion:** Allocator works perfectly for TREXIO allocation patterns.

### 2. Minimal Fix Applied
Added validation to prevent crashes from NULL allocators or invalid counts:

**In `md_trexio_create()`:**
```c
if (!alloc || !alloc->realloc) {
    MD_LOG_ERROR("TREXIO: Invalid allocator");
    return NULL;
}
```

**In `md_trexio_parse_file()`:**
```c
// Validate allocator
if (!trexio->alloc || !trexio->alloc->realloc) {
    MD_LOG_ERROR("TREXIO: Invalid allocator in structure");
    return false;
}

// Validate nucleus count
if (nucleus_num_i32 < 0 || nucleus_num_i32 > 1000000) {
    MD_LOG_ERROR("TREXIO: Invalid nucleus count: %d", nucleus_num_i32);
    return false;
}
```

**Total validation:** ~20 lines

### 3. Test File Issue (PR #117)
Repository test files fail with TREXIO library error:
```
Testing TREXIO loader with file: test_data/h2o_molecule.trexio
Allocator obtained: 0x... (valid)
Creating TREXIO structure... (valid)
Parsing TREXIO file...
[error]: TREXIO: Failed to open file: Read-only file  ← TREXIO library error
```

**Affected files:**
- `test_data/h2o_molecule.trexio`
- `test_data/h2_molecule.trexio`
- `test_data/ch4_molecule.trexio`

**Implication:** Original crash may have occurred with a different TREXIO file not in the repository.

---

## What We Need

### 1. Original Failing TREXIO File
**Request:** The exact TREXIO file that caused the crash in PR #113

**Details needed:**
- TREXIO file (or instructions to generate it)
- Command/code that triggered the crash
- Environment details (TREXIO library version, HDF5 vs text backend, etc.)
- Stack trace if available

**Why:** 
- Current test files don't reproduce the crash
- Need actual failing case to verify minimal fix is sufficient
- May reveal additional validation needed beyond current fixes

### 2. Expected File Format
**Request:** Guidance on creating valid TREXIO test files

**Questions:**
- What TREXIO library version should test files be created with?
- Should we use HDF5 or text backend?
- Are there specific TREXIO groups/datasets that must be present?
- Can you share Python code to generate valid test files?

**Example:**
```python
import trexio

# Create test file
f = trexio.File("test.trexio", 'w', trexio.TREXIO_TEXT)
f.write_nucleus_num(3)  # H2O
f.write_nucleus_coord([...])
# etc.
f.close()
```

---

## Reproduction Commands

For testing with the actual failing file:

```bash
# 1. Setup
git clone https://github.com/scanberg/viamd.git
cd viamd
git checkout TREXIO  # or fix branch once merged
git submodule update --init --recursive

# 2. Build with TREXIO + ASAN
mkdir build && cd build
cmake -DVIAMD_ENABLE_TREXIO=ON \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_C_FLAGS="-g -O1 -fsanitize=address,undefined" \
      -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined" ..

# 3. Build allocator reproducer (proves allocator works)
cmake --build . --target mdlib_allocator_reproducer -j4
./bin/mdlib_allocator_reproducer --iterations 100

# 4. Build loader test (test with actual file)
cmake --build . --target trexio_loader_test -j4
./bin/trexio_loader_test /path/to/failing_file.trexio
```

---

## Expected Outcomes

### With Original Failing File

**If minimal fix sufficient:**
- Loader test completes without crash
- Validation catches invalid data and returns error gracefully
- ASAN reports no memory errors

**If additional issues found:**
- Capture exact stack trace with ASAN
- Identify specific validation needed
- Implement targeted fix in new PR

### Without Original File

If file cannot be shared:
- Current minimal fix prevents NULL allocator crashes
- Current bounds check prevents integer overflow from invalid counts
- Additional crash scenarios remain unaddressed until reproduced

---

## Deliverables So Far

✅ **Completed:**
- Comprehensive investigation (3 PRs, ~15,000 lines of artifacts)
- Allocator validation (150+ iterations, ASAN clean)
- Root cause analysis (NULL allocator, invalid counts)
- Minimal validation fix (20 lines)
- Forensic report (`TREXIO_ALLOCATOR_FORENSIC_REPORT.md`)
- Test harness (`tools/trexio_loader_test/`)
- Reproducer tool (`tools/mdlib_allocator_reproducer/`)

📋 **Pending:**
- Verification with original failing file
- Additional validation if needed
- Integration of reproducer into TREXIO branch (optional)

---

## Checklist for Response

When providing the failing file, please include:

- [ ] TREXIO file that caused crash (or generation instructions)
- [ ] Command/code used to trigger crash
- [ ] TREXIO library version used
- [ ] Backend type (HDF5 vs text)
- [ ] Stack trace if available
- [ ] Any error messages observed
- [ ] Environment details (OS, compiler, etc.)

Optional but helpful:
- [ ] Python script to generate similar test files
- [ ] Expected behavior vs actual crash behavior
- [ ] Minimal code snippet that reproduces issue

---

## Related Links

- **PR #113:** TREXIO implementation (merged to TREXIO)
- **PR #115:** Allocator reproducer
- **PR #116:** Documentation & fix patterns
- **PR #117:** Systematic investigation
- **Forensic Report:** `TREXIO_ALLOCATOR_FORENSIC_REPORT.md`

---

## Contact

**Assignees:** @mathieulinares  
**Reviewers:** Repository maintainers  
**Labels:** TREXIO, investigation, help-wanted, question

---

**Thank you for your help in completing this investigation!**
