# TREXIO Allocator Initialization Fix

## Problem Summary

The TREXIO loader implementation encounters an allocator crash during `md_alloc()` in the parse function, despite:
- Direct TREXIO API tests working perfectly
- The implementation being complete and following the VeloxChem pattern

This indicates an **allocator initialization issue in the test context** rather than a TREXIO implementation problem.

## Root Cause Analysis

The mdlib heap allocator (`md_get_heap_allocator()`) is a static allocator defined in `ext/mdlib/src/core/md_allocator.c`:

```c
static struct md_allocator_i _heap_allocator = {
    NULL,
    realloc_internal,
};

md_allocator_i* md_get_heap_allocator(void) {
    return &_heap_allocator;
}
```

This allocator **should always work** and requires no initialization. However, crashes can occur if:

1. **The allocator pointer is NULL** when passed to `md_alloc()`
2. **The allocator symbols are not linked** into the test binary
3. **The function pointer is invalid** due to incorrect struct initialization

## Common Causes and Solutions

### 1. NULL Allocator Pointer

**Problem**: The allocator parameter passed to a function is NULL.

**Fix**: Always check allocator validity before use:

```c
bool md_trexio_parse_file(md_trexio_t* trexio, str_t filename, md_allocator_i* alloc) {
    // Add NULL check
    if (!alloc) {
        MD_LOG_ERROR("NULL allocator passed to TREXIO parse function");
        return false;
    }
    
    // Verify allocator has required function
    if (!alloc->realloc) {
        MD_LOG_ERROR("Allocator missing realloc function");
        return false;
    }
    
    // Now safe to use
    void* data = md_alloc(alloc, size);
    // ...
}
```

### 2. Allocator Not Initialized in Tests

**Problem**: Test code doesn't properly get the allocator.

**VeloxChem Pattern** (from `test_vlx.c`):
```c
UTEST(vlx, vlx_parse) {
    md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());  // ✓ Correct
    bool result = md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.out"));
    // ...
    md_vlx_destroy(vlx);
}
```

**Correct TREXIO Test Pattern**:
```c
UTEST(trexio, trexio_parse) {
    md_allocator_i* alloc = md_get_heap_allocator();  // Get allocator first
    ASSERT_NE(alloc, NULL);  // Verify it's not NULL
    
    md_trexio_t* trexio = md_trexio_create(alloc);
    ASSERT_NE(trexio, NULL);
    
    bool result = md_trexio_parse_file(trexio, STR_LIT(MD_UNITTEST_DATA_DIR "/trexio/h2o.trexio"), alloc);
    ASSERT_TRUE(result);
    
    md_trexio_destroy(trexio, alloc);
}
```

### 3. System Loader Interface

The `md_system_loader_i` interface expects specific function signatures. Ensure they match:

```c
// Correct signature
static bool trexio_sys_init_from_file(
    struct md_system_o* inst,
    str_t filename,
    struct md_allocator_i* alloc,
    uint32_t flags
) {
    // MUST check allocator before use
    if (!alloc || !alloc->realloc) {
        MD_LOG_ERROR("Invalid allocator in TREXIO loader");
        return false;
    }
    
    // Cast inst to TREXIO type
    md_trexio_t* trexio = (md_trexio_t*)inst;
    
    // Use allocator safely
    void* data = md_alloc(alloc, size);
    // ...
}
```

### 4. CMake Linking Issues

**Problem**: Test binary doesn't link allocator symbols.

**Fix**: Ensure test CMakeLists.txt includes mdlib:

```cmake
if (MD_ENABLE_TREXIO)
    add_executable(test_trexio test_trexio.c)
    target_link_libraries(test_trexio PRIVATE mdlib)  # Must link mdlib!
    target_include_directories(test_trexio PRIVATE ${TREXIO_INCLUDE_DIRS})
    target_link_libraries(test_trexio PRIVATE ${TREXIO_LIBRARIES})
endif()
```

## Implementation Checklist

- [ ] Add NULL checks for allocator parameter in all TREXIO functions
- [ ] Verify allocator->realloc is not NULL before calling md_alloc
- [ ] Follow VeloxChem test pattern: get allocator, verify not NULL, pass to functions
- [ ] Ensure test binary links with mdlib
- [ ] Add logging for allocator validation failures
- [ ] Test with both heap and arena allocators to verify generality

## Testing Pattern

```c
// test_trexio.c
#include "utest.h"
#include <md_trexio.h>
#include <core/md_allocator.h>

UTEST(trexio, allocator_initialization) {
    // Verify heap allocator works
    md_allocator_i* alloc = md_get_heap_allocator();
    ASSERT_NE(alloc, NULL);
    ASSERT_NE(alloc->realloc, NULL);
    
    // Test basic allocation
    void* test_mem = md_alloc(alloc, 64);
    ASSERT_NE(test_mem, NULL);
    md_free(alloc, test_mem, 64);
}

UTEST(trexio, parse_with_heap_allocator) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    md_trexio_t* trexio = md_trexio_create(alloc);
    ASSERT_NE(trexio, NULL);
    
    // Parse test file
    bool result = md_trexio_parse_file(
        trexio, 
        STR_LIT(MD_UNITTEST_DATA_DIR "/trexio/water.trexio"),
        alloc
    );
    ASSERT_TRUE(result);
    
    // Verify data loaded
    EXPECT_GT(md_trexio_num_atoms(trexio), 0);
    
    md_trexio_destroy(trexio, alloc);
}

UTEST(trexio, parse_with_arena_allocator) {
    md_allocator_i* arena = md_arena_allocator_create(
        md_get_heap_allocator(), 
        MEGABYTES(1)
    );
    
    md_trexio_t* trexio = md_trexio_create(arena);
    ASSERT_NE(trexio, NULL);
    
    bool result = md_trexio_parse_file(
        trexio, 
        STR_LIT(MD_UNITTEST_DATA_DIR "/trexio/water.trexio"),
        arena
    );
    ASSERT_TRUE(result);
    
    md_arena_allocator_destroy(arena);
}
```

## Quick Fix Template

If the TREXIO loader already exists, add this at the start of parse functions:

```c
bool md_trexio_parse_file(md_trexio_t* trexio, str_t filename, md_allocator_i* alloc) {
    // === ADD THESE CHECKS ===
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
    // === END CHECKS ===
    
    // Existing parse code...
    void* nucleus_charge = md_alloc(alloc, sizeof(double) * num_atoms);
    // ...
}
```

## References

- `ext/mdlib/src/core/md_allocator.c` - Allocator implementation
- `ext/mdlib/unittest/test_vlx.c` - Reference test pattern
- `ext/mdlib/unittest/test_allocator.c` - Allocator tests
- `ext/mdlib/src/md_vlx.c` - VeloxChem loader (reference implementation)
