#!/bin/bash

# Script to apply mdlib TREXIO patch in an idempotent way
# This script can be run multiple times without issues

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIAMD_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
MDLIB_DIR="${VIAMD_ROOT}/ext/mdlib"
PATCH_FILE="${VIAMD_ROOT}/docs/mdlib_trexio_original.patch"

echo "Applying mdlib TREXIO patch..."

# Check if mdlib directory exists
if [ ! -d "${MDLIB_DIR}" ]; then
    echo "Error: mdlib directory not found at ${MDLIB_DIR}"
    echo "Please run 'git submodule update --init --recursive' first"
    exit 1
fi

# Check if patch file exists
if [ ! -f "${PATCH_FILE}" ]; then
    echo "Error: Patch file not found at ${PATCH_FILE}"
    exit 1
fi

cd "${MDLIB_DIR}"

# Check if patch is already applied by looking for the TREXIO option in CMakeLists.txt
if grep -q "MD_ENABLE_TREXIO" CMakeLists.txt 2>/dev/null && \
   [ -f "src/md_trexio.c" ] && [ -f "src/md_trexio.h" ]; then
    echo "Patch appears to be already applied. Skipping."
    exit 0
fi

# Reset any partial application of the patch
if [ -f "src/md_trexio.c" ] || [ -f "src/md_trexio.h" ]; then
    echo "Cleaning up partial patch application..."
    rm -f src/md_trexio.c src/md_trexio.h
fi

# Revert CMakeLists.txt if it has uncommitted changes
if git diff --quiet CMakeLists.txt; then
    : # No changes, nothing to do
else
    echo "Reverting uncommitted changes to CMakeLists.txt..."
    git checkout -- CMakeLists.txt
fi

# Apply the patch
echo "Applying patch from: ${PATCH_FILE}"
if git apply "${PATCH_FILE}"; then
    echo "Patch applied successfully!"
    echo ""
    echo "TREXIO support files added:"
    echo "  - CMakeLists.txt (modified)"
    echo "  - src/md_trexio.c (new)"
    echo "  - src/md_trexio.h (new)"
else
    echo "Error: Failed to apply patch"
    echo ""
    echo "Troubleshooting:"
    echo "1. Ensure mdlib is at commit 06da5a3:"
    echo "   cd ${MDLIB_DIR} && git log --oneline -1"
    echo ""
    echo "2. Try resetting mdlib and applying again:"
    echo "   cd ${MDLIB_DIR} && git reset --hard HEAD && git clean -fd"
    echo "   cd ${VIAMD_ROOT} && ./scripts/apply_mdlib_trexio_patch.sh"
    exit 1
fi
