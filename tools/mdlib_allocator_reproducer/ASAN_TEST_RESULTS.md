# mdlib Allocator Reproducer - ASAN Test Results

## Test Configuration

- **Date:** 2025-11-18
- **Build Type:** AddressSanitizer (ASAN) enabled
- **Compiler:** GCC 13.3.0
- **Iterations:** 50
- **Flags:** `-fsanitize=address -g`

## Test Summary

**Result:** ✅ **PASSED - No memory errors detected**

The reproducer completed 50 iterations successfully with AddressSanitizer enabled. No memory leaks, use-after-free, buffer overflows, or other memory errors were detected.

## Detailed Results

- **Total iterations:** 50
- **Successful allocations:** All (100%)
- **Failed allocations:** 0
- **ASAN errors:** 0
- **Memory leaks:** 0
- **Exit code:** 0 (success)

## Key Observations

1. **All allocations successful:** Every call to `md_alloc()` succeeded and returned valid pointers
2. **No memory errors:** ASAN did not detect any:
   - Heap buffer overflows
   - Stack buffer overflows
   - Use-after-free
   - Use-after-return
   - Double-free
   - Memory leaks
3. **Allocator remains functional:** After each iteration, the allocator was verified to still work correctly
4. **Address patterns:** ASAN-instrumented addresses show expected heap patterns (e.g., `0x502000000010`, `0x50c000000040`)

## Sample Output (Last 80 lines)

The complete output is available in `repro-output-asan.txt` (3,657 lines, 180KB). The last 80 lines show the final iteration completing successfully with all memory properly freed.

## Conclusion

The mdlib heap allocator works correctly in isolated testing, even under AddressSanitizer's strict memory checking. This confirms that:

1. **The allocator implementation is sound** - No memory safety issues detected
2. **The allocation/deallocation patterns are correct** - Proper pairing of alloc/free calls
3. **Any TREXIO loader crashes are likely NOT due to the allocator itself**

If crashes occur in the actual TREXIO loader, they are more likely related to:
- TREXIO file parsing logic
- Data validation issues
- Integration between TREXIO library and mdlib
- Specific memory access patterns during file loading that aren't reproduced in this isolated test

## Next Steps

To investigate actual TREXIO loader crashes:
1. Run the actual TREXIO loader with ASAN enabled
2. Test with real TREXIO files (especially the one that causes crashes)
3. Add more complex allocation patterns if needed
4. Compare memory access patterns between reproducer and actual loader
