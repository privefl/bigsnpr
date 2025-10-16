# Test Import Error Fix

## Issue

When trying to run the test file directly, you encountered:
```
ImportError: attempted relative import with no known parent package
```

## Root Cause

The test file (`ldpred2_python/tests/test_ldpred2.py`) was using relative imports:
```python
from ldpred2_inf import ldpred2_inf
from ldpred2_grid import ldpred2_grid
```

When you run a Python file directly as a script (e.g., `python3 ldpred2_python/tests/test_ldpred2.py`), Python doesn't recognize it as part of a package, so relative imports fail.

## Solution

### 1. Fixed the Test File

Changed from relative imports to absolute imports:

**Before:**
```python
from ldpred2_inf import ldpred2_inf
from ldpred2_grid import ldpred2_grid
from utils import standardize_betas
```

**After:**
```python
# Add parent directory to path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

# Use absolute imports
from ldpred2_python.ldpred2_inf import ldpred2_inf
from ldpred2_python.ldpred2_grid import ldpred2_grid
from ldpred2_python.utils import standardize_betas
```

### 2. Created Test Runners

To make running tests easier, three test runners were created:

#### Option 1: Main Test Runner (Recommended)
```bash
python3 run_tests.py
```

**Location**: `run_tests.py` in project root

**Features**:
- Automatically sets up Python path
- Checks for dependencies
- Discovers and runs all tests
- Provides clear output
- Returns proper exit codes

#### Option 2: Simple Test Runner
```bash
python3 ldpred2_python/tests/run_tests_simple.py
```

**Location**: `ldpred2_python/tests/run_tests_simple.py`

**Features**:
- Can be run from any directory
- Explicitly imports all test classes
- Good for debugging specific tests

#### Option 3: Unittest Discovery
```bash
python3 -m unittest discover ldpred2_python/tests
```

**Built-in Python method**:
- Uses Python's unittest discovery
- No custom runner needed
- Standard approach

## How to Run Tests

### From Project Root

```bash
cd /path/to/bigsnpr

# Recommended
python3 run_tests.py

# Or use unittest
python3 -m unittest discover ldpred2_python/tests

# Or use the simple runner
python3 ldpred2_python/tests/run_tests_simple.py
```

### From Tests Directory

```bash
cd /path/to/bigsnpr/ldpred2_python/tests

# Use the simple runner
python3 run_tests_simple.py

# Or use unittest with parent directory
cd ../..
python3 -m unittest discover ldpred2_python/tests
```

### With pytest (if installed)

```bash
pip install pytest
pytest ldpred2_python/tests/
```

## What NOT to Do

❌ **Don't run the test file directly:**
```bash
# This will cause ImportError
python3 ldpred2_python/tests/test_ldpred2.py
```

❌ **Don't use relative imports when running as script:**
```python
# This fails when run directly
from ..ldpred2_inf import ldpred2_inf
```

## Understanding Python Imports

### Why This Happens

Python has two ways to run files:

1. **As a script**: `python3 file.py`
   - `__name__` is `"__main__"`
   - File is not part of a package
   - Relative imports don't work

2. **As a module**: `python3 -m package.module`
   - `__name__` is `"package.module"`
   - File is part of a package
   - Relative imports work

### The Fix

The test file now:
1. Adds the project root to `sys.path`
2. Uses absolute imports (`from ldpred2_python.module import ...`)
3. Works when run directly or as part of a package

## Verification

To verify the fix works:

```bash
# 1. Check dependencies are installed
python3 -c "import numpy, scipy, pandas; print('Dependencies OK')"

# 2. Run the test runner
python3 run_tests.py

# 3. Should see output like:
#    ✓ Dependencies found
#    test_standardize_betas (test_ldpred2.TestUtils) ... ok
#    test_ld_scores (test_ldpred2.TestUtils) ... ok
#    ...
#    ✓ All tests passed!
```

## Test Structure

The test file includes these test classes:

1. **TestUtils** - Tests for utility functions
   - `test_standardize_betas`
   - `test_ld_scores`
   - `test_check_dimensions`

2. **TestLDpred2Inf** - Tests for infinitesimal model
   - `test_basic_functionality`
   - `test_reproducibility`

3. **TestLDpred2Grid** - Tests for grid search
   - `test_gibbs_one_basic`
   - `test_grid_basic`
   - `test_reproducibility`

4. **TestLDpred2Auto** - Tests for automatic inference
   - `test_auto_basic`
   - `test_multiple_chains`

5. **TestDataGeneration** - Tests for synthetic data
   - `test_ld_matrix_generation`
   - `test_gwas_sumstats_generation`
   - `test_complete_dataset`

## Additional Notes

### Running Specific Tests

```bash
# Run a specific test class
python3 -m unittest ldpred2_python.tests.test_ldpred2.TestUtils

# Run a specific test method
python3 -m unittest ldpred2_python.tests.test_ldpred2.TestUtils.test_standardize_betas
```

### Verbose Output

```bash
# More detailed output
python3 -m unittest discover -v ldpred2_python/tests
```

### Quiet Output

```bash
# Less output
python3 -m unittest discover -q ldpred2_python/tests
```

## Summary

The `ImportError: attempted relative import with no known parent package` error has been fixed by:

1. ✅ Converting relative imports to absolute imports in test file
2. ✅ Adding proper path setup in test file
3. ✅ Creating dedicated test runners (`run_tests.py`, `run_tests_simple.py`)
4. ✅ Documenting correct ways to run tests
5. ✅ Updating all documentation with proper test commands

**Always use one of the test runners or `python3 -m unittest discover` instead of running the test file directly.**

