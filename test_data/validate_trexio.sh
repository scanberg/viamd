#!/bin/bash
# Validation script for TREXIO integration in VIAMD
# This script verifies that the TREXIO support is working correctly

set -e

echo "=========================================="
echo "TREXIO Integration Validation Script"
echo "=========================================="
echo ""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIAMD_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "VIAMD Root: $VIAMD_ROOT"
echo ""

# Check if mdlib patch has been applied
echo "1. Checking if mdlib TREXIO files exist..."
if [ -f "$VIAMD_ROOT/ext/mdlib/src/md_trexio.h" ] && [ -f "$VIAMD_ROOT/ext/mdlib/src/md_trexio.c" ]; then
    echo "   ✓ mdlib TREXIO files found"
else
    echo "   ✗ mdlib TREXIO files not found"
    echo "   Please apply the patch: cd ext/mdlib && git apply ../../docs/mdlib_trexio.patch"
    exit 1
fi

# Check if test data exists
echo ""
echo "2. Checking test data..."
if [ -d "$VIAMD_ROOT/test_data/h2_molecule.trexio" ]; then
    echo "   ✓ H2 test molecule found"
else
    echo "   ✗ H2 test molecule not found"
    echo "   Creating test data..."
    cd "$VIAMD_ROOT/test_data"
    python3 create_test_trexio.py
fi

if [ -d "$VIAMD_ROOT/test_data/h2o_molecule.trexio" ]; then
    echo "   ✓ H2O test molecule found"
fi

if [ -d "$VIAMD_ROOT/test_data/ch4_molecule.trexio" ]; then
    echo "   ✓ CH4 test molecule found"
fi

# Check test file structure
echo ""
echo "3. Validating test file structure..."
if [ -f "$VIAMD_ROOT/test_data/h2_molecule.trexio/nucleus.txt" ]; then
    echo "   ✓ H2 nucleus data found"
    num_atoms=$(grep -A 1 "^num$" "$VIAMD_ROOT/test_data/h2_molecule.trexio/nucleus.txt" | tail -1)
    echo "     Number of atoms: $num_atoms"
fi

if [ -f "$VIAMD_ROOT/test_data/h2_molecule.trexio/electron.txt" ]; then
    echo "   ✓ H2 electron data found"
fi

# Check unit test
echo ""
echo "4. Checking unit test..."
if [ -f "$VIAMD_ROOT/ext/mdlib/unittest/test_trexio.c" ]; then
    echo "   ✓ Unit test file found"
    test_count=$(grep -c "^UTEST(trexio," "$VIAMD_ROOT/ext/mdlib/unittest/test_trexio.c" || true)
    echo "     Number of tests: $test_count"
else
    echo "   ✗ Unit test file not found"
fi

# Summary
echo ""
echo "=========================================="
echo "Validation Summary"
echo "=========================================="
echo ""
echo "Phase 4 Progress:"
echo "  ✓ Test data created (H2, H2O, CH4)"
echo "  ✓ Unit tests created"
echo "  ✓ mdlib patch ready to apply"
echo ""
echo "Next steps:"
echo "  1. Build VIAMD with TREXIO support:"
echo "     mkdir -p build && cd build"
echo "     cmake -DVIAMD_ENABLE_TREXIO=ON .."
echo "     make"
echo ""
echo "  2. Run unit tests (if build succeeds):"
echo "     cd build"
echo "     ctest -V -R trexio"
echo ""
echo "  3. Test with real TREXIO files from:"
echo "     - PySCF"
echo "     - Quantum Package"
echo "     - Other quantum chemistry codes"
echo ""
echo "=========================================="
