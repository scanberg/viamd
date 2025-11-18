# TREXIO Allocator Fix - Implementation Checklist

Use this checklist to apply the allocator initialization fix to the TREXIO implementation.

## Pre-Implementation

- [ ] Locate TREXIO implementation files (md_trexio.c, md_trexio.h)
- [ ] Locate TREXIO test files (test_trexio.c)
- [ ] Read TREXIO_ALLOCATOR_FIX.md for context
- [ ] Read TREXIO_IMPLEMENTATION_REFERENCE.md for patterns
- [ ] Backup current implementation

## Code Changes

### md_trexio.c - Parse Function

- [ ] Add NULL check for `trexio` parameter at function start
- [ ] Add NULL check for `alloc` parameter
- [ ] Add NULL check for `alloc->realloc` function pointer
- [ ] Add logging for each validation failure
- [ ] Verify early return on validation failure

Example location: Start of `md_trexio_parse_file()`:
```c
if (!trexio) {
    MD_LOG_ERROR("NULL trexio structure");
    return false;
}

if (!alloc) {
    MD_LOG_ERROR("NULL allocator - cannot parse TREXIO file");
    return false;
}

if (!alloc->realloc) {
    MD_LOG_ERROR("Allocator has NULL realloc function");
    return false;
}
```

### md_trexio.c - Create Function

- [ ] Add NULL check for `alloc` parameter
- [ ] Add NULL check for `alloc->realloc`
- [ ] Check allocation result before using pointer
- [ ] Add error logging

Example location: Start of `md_trexio_create()`:
```c
if (!alloc) {
    MD_LOG_ERROR("NULL allocator passed to md_trexio_create");
    return NULL;
}

if (!alloc->realloc) {
    MD_LOG_ERROR("Allocator missing realloc function");
    return NULL;
}

md_trexio_t* trexio = (md_trexio_t*)md_alloc(alloc, sizeof(md_trexio_t));
if (!trexio) {
    MD_LOG_ERROR("Failed to allocate TREXIO structure");
    return NULL;
}
```

### md_trexio.c - Destroy Function

- [ ] Add NULL checks for both parameters
- [ ] Early return if either is NULL
- [ ] Check array pointers before freeing

Example:
```c
void md_trexio_destroy(md_trexio_t* trexio, md_allocator_i* alloc) {
    if (!trexio || !alloc) return;
    
    if (trexio->atom_coordinates) {
        md_free(alloc, trexio->atom_coordinates, /* size */);
    }
    // ... other arrays ...
    
    md_free(alloc, trexio, sizeof(md_trexio_t));
}
```

### md_trexio.c - System Loader Interface

- [ ] Add validation in `trexio_sys_init_from_file()`
- [ ] Add validation in `trexio_sys_init_from_str()` (if implemented)
- [ ] Add validation in `trexio_sys_free()`

Example location: Start of `trexio_sys_init_from_file()`:
```c
static bool trexio_sys_init_from_file(
    struct md_system_o* inst,
    str_t filename,
    struct md_allocator_i* alloc,
    uint32_t flags
) {
    if (!alloc || !alloc->realloc) {
        MD_LOG_ERROR("Invalid allocator in TREXIO system loader");
        return false;
    }
    
    // ... implementation ...
}
```

### md_trexio.c - All Other Functions

- [ ] Review every function that accepts `md_allocator_i*`
- [ ] Add validation to each
- [ ] Ensure consistent error messages
- [ ] Log function name in error messages

## Test Changes

### test_trexio.c - Setup

- [ ] Include correct headers (`utest.h`, `md_trexio.h`, `core/md_allocator.h`)
- [ ] Remove any global allocator variables
- [ ] Each test gets its own allocator

### test_trexio.c - Allocator Validation Test

- [ ] Add test: `UTEST(trexio, allocator_validation)`
- [ ] Verify `md_get_heap_allocator()` returns non-NULL
- [ ] Verify `alloc->realloc` is non-NULL
- [ ] Test basic allocation works

```c
UTEST(trexio, allocator_validation) {
    md_allocator_i* alloc = md_get_heap_allocator();
    ASSERT_NE(alloc, NULL);
    ASSERT_NE(alloc->realloc, NULL);
    
    void* test = md_alloc(alloc, 64);
    ASSERT_NE(test, NULL);
    md_free(alloc, test, 64);
}
```

### test_trexio.c - Create/Destroy Test

- [ ] Get allocator: `md_allocator_i* alloc = md_get_heap_allocator();`
- [ ] Assert allocator is not NULL
- [ ] Pass allocator to create function
- [ ] Assert result is not NULL
- [ ] Pass allocator to destroy function

```c
UTEST(trexio, create_and_destroy) {
    md_allocator_i* alloc = md_get_heap_allocator();
    ASSERT_NE(alloc, NULL);
    
    md_trexio_t* trexio = md_trexio_create(alloc);
    ASSERT_NE(trexio, NULL);
    
    md_trexio_destroy(trexio, alloc);
}
```

### test_trexio.c - Parse Test

- [ ] Get allocator first
- [ ] Create TREXIO structure with allocator
- [ ] Pass allocator to parse function
- [ ] Verify parse succeeded
- [ ] Check loaded data
- [ ] Clean up with allocator

```c
UTEST(trexio, parse_water) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    md_trexio_t* trexio = md_trexio_create(alloc);
    ASSERT_NE(trexio, NULL);
    
    bool result = md_trexio_parse_file(
        trexio,
        STR_LIT(MD_UNITTEST_DATA_DIR "/trexio/water.trexio"),
        alloc
    );
    ASSERT_TRUE(result);
    
    EXPECT_GT(md_trexio_num_atoms(trexio), 0);
    
    md_trexio_destroy(trexio, alloc);
}
```

### test_trexio.c - NULL Allocator Test

- [ ] Add test: `UTEST(trexio, null_allocator_handling)`
- [ ] Verify NULL allocator is rejected gracefully
- [ ] Verify no crash occurs

```c
UTEST(trexio, null_allocator_handling) {
    md_trexio_t* trexio = md_trexio_create(NULL);
    EXPECT_EQ(NULL, trexio);  // Should return NULL, not crash
}
```

## CMake Changes

### ext/mdlib/unittest/CMakeLists.txt

- [ ] Find section with TREXIO test configuration
- [ ] Verify `target_link_libraries` includes `mdlib`
- [ ] Verify TREXIO includes are correct
- [ ] Verify test is only built when `MD_ENABLE_TREXIO` is ON

```cmake
if (MD_ENABLE_TREXIO)
    add_executable(test_trexio test_trexio.c)
    target_link_libraries(test_trexio PRIVATE mdlib)  # Critical!
    target_include_directories(test_trexio PRIVATE ${TREXIO_INCLUDE_DIRS})
    target_link_libraries(test_trexio PRIVATE ${TREXIO_LIBRARIES})
    add_test(NAME test_trexio COMMAND test_trexio)
endif()
```

## Build and Test

### Build

- [ ] Clean build directory: `rm -rf build && mkdir build && cd build`
- [ ] Configure with TREXIO: `cmake -DMD_ENABLE_TREXIO=ON ..`
- [ ] Build: `make -j4`
- [ ] Verify no compilation errors
- [ ] Verify no warnings about allocator usage

### Test Execution

- [ ] Run test binary: `./unittest/test_trexio`
- [ ] Verify all tests pass
- [ ] Check for memory leaks (if tools available): `valgrind ./unittest/test_trexio`
- [ ] Verify error messages appear correctly for NULL allocator test

### Integration Test

- [ ] Try loading actual TREXIO file with system loader
- [ ] Verify geometry data loads correctly
- [ ] Verify atom count matches file
- [ ] Verify atomic numbers are correct
- [ ] Verify coordinates are in Angstrom (converted from Bohr)

## Verification

### Code Review

- [ ] All functions that use allocator have validation
- [ ] Error messages are clear and helpful
- [ ] Early returns prevent use of invalid allocator
- [ ] Logging includes function names
- [ ] Code follows VeloxChem pattern

### Testing

- [ ] Allocator validation test passes
- [ ] Create/destroy test passes
- [ ] Parse test passes
- [ ] NULL allocator test passes (returns NULL, doesn't crash)
- [ ] No memory leaks detected
- [ ] Integration with viamd works

### Documentation

- [ ] Update TREXIO_IMPLEMENTATION_STATUS.md with fix status
- [ ] Note any deviations from reference implementation
- [ ] Document any additional changes needed

## Completion

- [ ] All code changes committed
- [ ] All tests passing
- [ ] Documentation updated
- [ ] PR ready for review

## Notes

Use this space to track issues or decisions:

```
Example:
- Added extra validation in trexio_read_basis_data() 
- Test file location: ext/mdlib/unittest/data/trexio/water.trexio
- Modified CMakeLists.txt to copy test data
```

## Success Criteria

✅ The fix is successful when:
1. All unit tests pass
2. NULL allocator returns error (not crash)
3. Parse function loads geometry correctly
4. No memory leaks
5. Error messages are clear
6. Code follows VeloxChem pattern
7. Documentation is updated
