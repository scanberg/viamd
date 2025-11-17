# VIAMD Scripts

This directory contains helper scripts for building and maintaining VIAMD.

## apply_mdlib_trexio_patch.sh

Applies the mdlib TREXIO patch in an idempotent way. This script can be run multiple times safely.

### Usage

```bash
./scripts/apply_mdlib_trexio_patch.sh
```

### Features

- **Idempotent**: Can be run multiple times without issues
- **Auto-cleanup**: Automatically cleans up partial patch applications
- **Error recovery**: If you manually removed files or made changes, just run the script again

### Scenarios Handled

1. **Fresh application**: Applies the patch for the first time
2. **Already applied**: Detects that the patch is already applied and skips
3. **Partial application**: Cleans up and reapplies if some files are missing
4. **Manual changes**: Reverts changes and reapplies the patch cleanly

### Troubleshooting

If you encounter any issues with the patch:

1. Simply run the script again: `./scripts/apply_mdlib_trexio_patch.sh`
2. The script will automatically detect the state and fix any issues

For manual patch application, see the alternative instructions in `docs/TREXIO_SUPPORT.md`.
