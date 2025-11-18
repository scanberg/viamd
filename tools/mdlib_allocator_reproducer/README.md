# mdlib Allocator Reproducer

A minimal, self-contained diagnostic tool to investigate mdlib allocator crashes observed during TREXIO parsing.

## Purpose

This tool reproduces the exact `md_get_heap_allocator()` / `md_alloc()` call pattern used by the TREXIO loader to help diagnose allocator-related crashes. It emits detailed diagnostic logging including:

- Allocator pointer and internal state
- Allocation sizes and addresses
- Return values from allocator operations
- Verification that the allocator remains functional after each operation

## Building

The reproducer is built as part of the viamd project:

```bash
# From the viamd root directory
mkdir -p build && cd build
cmake ..
make mdlib_allocator_reproducer
```

The tool will be located at `build/tools/mdlib_allocator_reproducer/mdlib_allocator_reproducer`.

## Running

### Using the Helper Script (Recommended)

The easiest way to run the reproducer and capture output is using the provided script:

```bash
# From the viamd root directory
./scripts/run_trexio_allocator_repro.sh
```

This script will:
- Build the reproducer (or use an existing build with `SKIP_BUILD=1`)
- Run 100 iterations by default
- Capture all output to `tools/mdlib_allocator_reproducer/repro-output.txt`
- Display a summary of results

**Environment variables:**
- `ITERATIONS=N` - Number of iterations to run (default: 100)
- `ENABLE_ASAN=1` - Enable AddressSanitizer for memory error detection
- `SKIP_BUILD=1` - Skip the build step and use existing build
- `BUILD_DIR=/path/to/build` - Specify build directory (default: ./build)

**Examples:**

```bash
# Run with 1000 iterations
ITERATIONS=1000 ./scripts/run_trexio_allocator_repro.sh

# Run with AddressSanitizer
ENABLE_ASAN=1 ./scripts/run_trexio_allocator_repro.sh

# Use existing build and run 500 iterations
SKIP_BUILD=1 ITERATIONS=500 ./scripts/run_trexio_allocator_repro.sh
```

### Basic Usage

Run a single iteration of the allocation pattern:

```bash
./build/tools/mdlib_allocator_reproducer/mdlib_allocator_reproducer
```

### Multiple Iterations

Test with multiple iterations to catch intermittent issues:

```bash
./build/tools/mdlib_allocator_reproducer/mdlib_allocator_reproducer --iterations 100
```

### Using CMake Target

You can also use the provided CMake target to run the reproducer:

```bash
cd build
make run_allocator_reproducer
```

This runs the reproducer with 100 iterations by default.

## Command Line Options

- `--iterations N` - Run N iterations of the allocation pattern (default: 1)
- `--help` - Show help message

## What It Tests

The reproducer mimics the TREXIO loader's allocation pattern:

1. **Get heap allocator** - Calls `md_get_heap_allocator()` (same as TREXIO loader)
2. **Create main structure** - Allocates a structure similar to `md_trexio_t`
3. **Simulate TREXIO allocations** - Performs allocations similar to parsing TREXIO data:
   - Nucleus coordinates array
   - Nucleus charges array
   - Nucleus label strings (array of strings)
   - Basis set shell data
   - Basis set primitive data (exponents)
   - Molecular orbital energies
4. **Free all data** - Frees everything in reverse order
5. **Verify allocator** - Tests that the allocator remains functional after all operations

## Diagnostic Output

The tool produces detailed diagnostic output for each step:

```
[DIAGNOSTIC] mdlib Allocator Reproducer Started
[DIAGNOSTIC] Iterations: 1
[DIAGNOSTIC] ========================================
[DIAGNOSTIC] 
--- Iteration 1/1 ---
[DIAGNOSTIC] Calling md_get_heap_allocator()...
[DIAGNOSTIC]   [After md_get_heap_allocator()]
[DIAGNOSTIC]     Allocator address: 0x...
[DIAGNOSTIC]     Allocator inst:    0x...
[DIAGNOSTIC]     Allocator realloc: 0x...
[DIAGNOSTIC]     Testing allocator with small allocation...
[DIAGNOSTIC]     Test allocation successful at: 0x...
[DIAGNOSTIC]     Test free completed
...
```

## Integration with CI

This tool can be integrated into CI pipelines to catch allocator regressions:

```yaml
- name: Test mdlib allocator
  run: |
    ./scripts/run_trexio_allocator_repro.sh
    
- name: Upload reproducer output
  uses: actions/upload-artifact@v3
  with:
    name: allocator-repro-output
    path: tools/mdlib_allocator_reproducer/repro-output.txt
```

Or run directly:

```yaml
- name: Test mdlib allocator
  run: |
    cd build
    ./tools/mdlib_allocator_reproducer/mdlib_allocator_reproducer --iterations 1000
```

## Troubleshooting

### Build Errors

If you get errors about missing mdlib headers:
- Ensure you've initialized submodules: `git submodule update --init --recursive`
- Ensure you're building from the viamd root directory

### Runtime Crashes

If the reproducer crashes:
1. Check the diagnostic output to see which allocation failed
2. Note the allocator addresses and state before the crash
3. Try running with fewer iterations to isolate the issue
4. Check if the crash is deterministic or intermittent

## Platform Support

This tool is platform-neutral and should work on:
- Linux (GCC, Clang)
- macOS (Clang)
- Windows (MSVC, MinGW)

The tool uses only standard C (C11) and mdlib APIs, with no platform-specific dependencies.

## Related Files

- Main implementation: `main.c`
- Build configuration: `CMakeLists.txt`
- This documentation: `README.md`

## Future Enhancements

Potential improvements for future iterations:

- Add TREXIO file reading to test with real data
- Add memory leak detection integration (Valgrind, ASan)
- Add stress testing with random allocation sizes
- Add multi-threaded testing for thread-local allocators
- Add comparison between heap and temp allocators
