#!/bin/bash

# Script to build and run the mdlib allocator reproducer with TREXIO support enabled
# Captures output for diagnostic purposes

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIAMD_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-${VIAMD_ROOT}/build}"
OUTPUT_FILE="${VIAMD_ROOT}/tools/mdlib_allocator_reproducer/repro-output.txt"
ITERATIONS=${ITERATIONS:-100}
ENABLE_ASAN=${ENABLE_ASAN:-0}
SKIP_BUILD=${SKIP_BUILD:-0}

echo "============================================"
echo "TREXIO Allocator Reproducer Build & Run"
echo "============================================"
echo "Build directory: ${BUILD_DIR}"
echo "Output file: ${OUTPUT_FILE}"
echo "Iterations: ${ITERATIONS}"
echo "ASAN enabled: ${ENABLE_ASAN}"
echo ""

echo "============================================"
echo "TREXIO Allocator Reproducer Build & Run"
echo "============================================"
echo "Build directory: ${BUILD_DIR}"
echo "Output file: ${OUTPUT_FILE}"
echo "Iterations: ${ITERATIONS}"
echo "ASAN enabled: ${ENABLE_ASAN}"
echo "Skip build: ${SKIP_BUILD}"
echo ""

# Clean or create build directory
if [ "${SKIP_BUILD}" = "0" ]; then
    if [ -d "${BUILD_DIR}" ]; then
        echo "Using existing build directory..."
    else
        echo "Creating build directory..."
        mkdir -p "${BUILD_DIR}"
    fi
    
    # Configure CMake with TREXIO enabled
    echo "Configuring CMake with TREXIO support..."
    cd "${BUILD_DIR}"
    
    CMAKE_FLAGS="-DVIAMD_ENABLE_TREXIO=ON"
    
    # Add ASAN flags if requested
    if [ "${ENABLE_ASAN}" = "1" ]; then
        echo "Enabling AddressSanitizer..."
        CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_C_FLAGS=-fsanitize=address -DCMAKE_CXX_FLAGS=-fsanitize=address -DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address"
    fi
    
    cmake ${CMAKE_FLAGS} "${VIAMD_ROOT}" || {
        echo "ERROR: CMake configuration failed"
        echo "You may need to install dependencies (e.g., OpenGL, X11)"
        echo "Or set SKIP_BUILD=1 to use an existing build"
        exit 1
    }
    
    # Build only the reproducer target
    echo ""
    echo "Building mdlib_allocator_reproducer target..."
    cmake --build . --target mdlib_allocator_reproducer -j$(nproc 2>/dev/null || echo 4) || {
        echo "ERROR: Build failed"
        exit 1
    }
else
    echo "Skipping build (SKIP_BUILD=1)"
    cd "${BUILD_DIR}"
fi

# Verify the binary exists
REPRO_BIN="${BUILD_DIR}/bin/mdlib_allocator_reproducer"
if [ ! -f "${REPRO_BIN}" ]; then
    echo "ERROR: Reproducer binary not found at ${REPRO_BIN}"
    exit 1
fi

# Run the reproducer and capture output
echo ""
echo "Running reproducer with ${ITERATIONS} iterations..."
echo "Output will be saved to: ${OUTPUT_FILE}"
echo ""

# Create output directory if needed
mkdir -p "$(dirname "${OUTPUT_FILE}")"

# Run and capture both stdout and stderr
{
    echo "============================================"
    echo "TREXIO Allocator Reproducer Output"
    echo "============================================"
    echo "Date: $(date)"
    echo "Host: $(hostname)"
    echo "OS: $(uname -a)"
    echo "Compiler: $(cmake --version | head -1)"
    echo "Iterations: ${ITERATIONS}"
    echo "ASAN: ${ENABLE_ASAN}"
    echo "============================================"
    echo ""
    
    "${REPRO_BIN}" --iterations "${ITERATIONS}"
    REPRO_EXIT=$?
    
    echo ""
    echo "============================================"
    echo "Reproducer Exit Code: ${REPRO_EXIT}"
    echo "============================================"
    
    exit ${REPRO_EXIT}
} > "${OUTPUT_FILE}" 2>&1

REPRO_EXIT=$?

# Display results
echo ""
echo "============================================"
echo "Results"
echo "============================================"
echo "Exit code: ${REPRO_EXIT}"
echo "Output saved to: ${OUTPUT_FILE}"
echo ""

if [ ${REPRO_EXIT} -eq 0 ]; then
    echo "✓ Reproducer completed successfully"
    echo ""
    echo "Last 30 lines of output:"
    tail -30 "${OUTPUT_FILE}" 2>/dev/null || echo "(Could not display output)"
else
    echo "✗ Reproducer failed with exit code ${REPRO_EXIT}"
    echo ""
    echo "Last 50 lines of output:"
    tail -50 "${OUTPUT_FILE}" 2>/dev/null || echo "(Could not display output)"
fi

echo ""
echo "Full output available at: ${OUTPUT_FILE}"
echo "============================================"

exit ${REPRO_EXIT}
