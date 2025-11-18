# TREXIO GUI Headless Test

This directory contains automated tests for the TREXIO GUI panel.

## Overview

The test suite validates:
- TREXIO file loading
- Metadata extraction
- Grid computation correctness
- File format compliance

## Running Tests Locally

### Prerequisites

- Python 3.6+
- viamd built with `-DVIAMD_ENABLE_TREXIO=ON`
- Test TREXIO files in `test_data/`

### Basic Test Run

```bash
cd tools/trexio/trexio_gui_test
python3 run_tests.py
```

### With Xvfb (Headless)

For headless testing (e.g., in CI):

```bash
# Install Xvfb
sudo apt-get install xvfb

# Run tests
xvfb-run -a python3 run_tests.py
```

## Test Structure

### Files

- `run_tests.py` - Main test runner
- `expected_output.json` - Reference values for validation
- `output/` - Test results and artifacts (created during test run)

### Test Cases

1. **File Existence Tests**
   - Verify test data files exist
   - Check TREXIO directory structure

2. **Metadata Tests**
   - Validate atom counts
   - Check basis set information
   - Verify MO data presence

3. **Grid Format Tests**
   - Check magic number (0x4754524F)
   - Validate version
   - Verify header structure

## Output

Test results are saved to `output/test_results.json`:

```json
{
  "tests": [
    {
      "test": "File exists: h2_molecule.trexio",
      "passed": true,
      "details": "File found: ..."
    }
  ],
  "passed": 8,
  "failed": 0,
  "timestamp": "2024-11-18 12:34:56"
}
```

## CI Integration

The tests are automatically run in CI via `.github/workflows/trexio-gui-test.yml`.

Artifacts uploaded:
- `test_results.json` - Test outcomes
- `trexio-debug-*.log` - Debug logs (if verbose logging enabled)
- `orbital_*.grid` - Generated grid files

## Extending Tests

To add new test cases:

1. Edit `run_tests.py` and add test methods
2. Update `expected_output.json` with reference values
3. Run locally to verify
4. Commit changes

Example test method:

```python
def test_my_feature(self):
    test_name = "My Feature Test"
    self.log(f"Running: {test_name}")
    
    # Your test logic here
    passed = True  # or False
    
    result = {
        'test': test_name,
        'passed': passed,
        'details': 'Test details'
    }
    
    self.results['tests'].append(result)
    if passed:
        self.results['passed'] += 1
    else:
        self.results['failed'] += 1
    
    return passed
```

## Troubleshooting

### Tests fail with "TREXIO library not found"

- Ensure TREXIO is installed: `pkg-config --modversion trexio`
- Check LD_LIBRARY_PATH includes TREXIO library path
- Run `sudo ldconfig` after installing TREXIO

### "No display" errors

- Use `xvfb-run -a` to create virtual display
- Or set `DISPLAY=:99` with Xvfb running

### Grid files not generated

- Check viamd built successfully with TREXIO support
- Verify test data files are valid
- Review debug logs in `tools/trexio/debug-*.log`

## Performance Benchmarks

Expected computation times (32³ grid):

| Molecule | Atoms | MOs | Time (approx) |
|----------|-------|-----|---------------|
| H2       | 2     | 2   | < 5s          |
| H2O      | 3     | 24  | < 30s         |
| CH4      | 5     | 40+ | < 60s         |

Times measured on CI runners (2 CPU cores, 7GB RAM).
