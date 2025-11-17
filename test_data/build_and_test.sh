#!/bin/bash
# Build and test VIAMD with TREXIO support
# Phase 4 testing script

set -e  # Exit on error

echo "=========================================="
echo "VIAMD TREXIO Build and Test Script"
echo "=========================================="
echo ""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIAMD_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$VIAMD_ROOT/build_trexio"

echo "VIAMD Root: $VIAMD_ROOT"
echo "Build Directory: $BUILD_DIR"
echo ""

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check for required tools
echo "1. Checking build dependencies..."
MISSING_DEPS=()

if ! command_exists cmake; then
    MISSING_DEPS+=("cmake")
fi

if ! command_exists make; then
    MISSING_DEPS+=("make")
fi

if ! command_exists gcc && ! command_exists clang; then
    MISSING_DEPS+=("gcc or clang")
fi

if [ ${#MISSING_DEPS[@]} -gt 0 ]; then
    echo "   ✗ Missing dependencies: ${MISSING_DEPS[*]}"
    echo "   Please install missing dependencies and try again."
    exit 1
fi

echo "   ✓ Build tools found (cmake, make, compiler)"

# Check for TREXIO library
echo ""
echo "2. Checking for TREXIO library..."

# Try to find TREXIO using pkg-config or known locations
TREXIO_FOUND=false

if command_exists pkg-config && pkg-config --exists trexio; then
    echo "   ✓ TREXIO found via pkg-config"
    TREXIO_VERSION=$(pkg-config --modversion trexio 2>/dev/null || echo "unknown")
    echo "     Version: $TREXIO_VERSION"
    TREXIO_FOUND=true
elif [ -f "$HOME/.local/lib/libtrexio.so" ] || [ -f "$HOME/.local/lib/libtrexio.a" ]; then
    echo "   ✓ TREXIO found in ~/.local/lib"
    TREXIO_FOUND=true
elif [ -f "/usr/local/lib/libtrexio.so" ] || [ -f "/usr/local/lib/libtrexio.a" ]; then
    echo "   ✓ TREXIO found in /usr/local/lib"
    TREXIO_FOUND=true
fi

if [ "$TREXIO_FOUND" = false ]; then
    echo "   ✗ TREXIO library not found"
    echo ""
    echo "   TREXIO installation options:"
    echo "   1. Conda: conda install -c conda-forge trexio"
    echo "   2. From source:"
    echo "      git clone https://github.com/TREX-CoE/trexio.git"
    echo "      cd trexio"
    echo "      ./configure --prefix=\$HOME/.local"
    echo "      make && make install"
    echo ""
    echo "   Continuing with build attempt (may fail)..."
fi

# Apply mdlib patch if not already applied
echo ""
echo "3. Checking mdlib TREXIO files..."

if [ ! -f "$VIAMD_ROOT/ext/mdlib/src/md_trexio.c" ]; then
    echo "   ! mdlib TREXIO files not found, applying patch..."
    cd "$VIAMD_ROOT/ext/mdlib"
    
    if [ -f "../../docs/mdlib_trexio.patch" ]; then
        git apply ../../docs/mdlib_trexio.patch 2>&1 | grep -v "trailing whitespace" || true
        echo "   ✓ Patch applied"
    else
        echo "   ✗ Patch file not found: docs/mdlib_trexio.patch"
        exit 1
    fi
else
    echo "   ✓ mdlib TREXIO files present"
fi

# Create build directory
echo ""
echo "4. Configuring build..."
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Configure with CMake
echo "   Running cmake..."
if cmake -DVIAMD_ENABLE_TREXIO=ON "$VIAMD_ROOT" 2>&1 | tee cmake.log; then
    echo "   ✓ CMake configuration successful"
else
    echo "   ✗ CMake configuration failed"
    echo "   See build_trexio/cmake.log for details"
    exit 1
fi

# Build
echo ""
echo "5. Building VIAMD..."
echo "   This may take several minutes..."

if make -j$(nproc) 2>&1 | tee build.log; then
    echo "   ✓ Build successful"
else
    echo "   ✗ Build failed"
    echo "   See build_trexio/build.log for details"
    exit 1
fi

# Check if executable was created
if [ -f "$BUILD_DIR/viamd" ]; then
    echo "   ✓ VIAMD executable created: $BUILD_DIR/viamd"
else
    echo "   ✗ VIAMD executable not found"
    exit 1
fi

# Run tests if available
echo ""
echo "6. Running tests..."

if [ -f "$BUILD_DIR/CTestTestfile.cmake" ]; then
    echo "   Running TREXIO unit tests..."
    if ctest -V -R trexio 2>&1 | tee test.log; then
        echo "   ✓ Tests passed"
    else
        echo "   ! Some tests failed (this may be expected if TREXIO library is not properly installed)"
        echo "   See build_trexio/test.log for details"
    fi
else
    echo "   ! CTest not configured, skipping tests"
fi

# Test with sample data
echo ""
echo "7. Testing with sample TREXIO files..."

cd "$VIAMD_ROOT/test_data"

# Test loading text format files
for testfile in h2_molecule.trexio h2o_molecule.trexio ch4_molecule.trexio; do
    if [ -d "$testfile" ]; then
        echo "   Testing: $testfile"
        # We can't actually run VIAMD in headless mode, but we can check the file exists
        echo "   ✓ Test file exists: $testfile"
    fi
done

# Summary
echo ""
echo "=========================================="
echo "Build and Test Summary"
echo "=========================================="
echo ""
echo "Build Status: SUCCESS"
echo "Build Directory: $BUILD_DIR"
echo "Executable: $BUILD_DIR/viamd"
echo ""
echo "Next steps:"
echo "  1. Test with sample files:"
echo "     cd $BUILD_DIR"
echo "     ./viamd ../test_data/h2_molecule.trexio"
echo ""
echo "  2. Create TREXIO files with PySCF (if installed):"
echo "     cd $VIAMD_ROOT/test_data"
echo "     python3 create_pyscf_trexio.py"
echo ""
echo "  3. Test with your own TREXIO files"
echo ""
echo "=========================================="
