# LDpred2 Python - Troubleshooting Guide

## Installation Issues

### Error: `ModuleNotFoundError: No module named 'numpy'`

**Problem**: Required dependencies are not installed.

**Solution**:
```bash
# Option 1: Install all dependencies at once
pip install -r requirements.txt

# Option 2: Install individually
pip install numpy scipy pandas scikit-learn matplotlib

# Option 3: Use the installation script
bash install_dependencies.sh
```

**For Python 3.12+**, ensure you're using compatible versions:
```bash
pip install numpy>=1.26.0 scipy>=1.11.0 pandas>=2.0.0 scikit-learn>=1.3.0
```

### Error: `NameError: name 'Union' is not defined`

**Problem**: This was a type hint compatibility issue with Python 3.12+.

**Solution**: This has been fixed in the latest version. The code now properly imports `Union` and `Optional` from the `typing` module where needed.

If you still see this error:
1. Make sure you have the latest version of the files
2. Check that the imports at the top of the file include:
   ```python
   from typing import Union, Optional, List, Dict
   ```

### Error: `pip: command not found`

**Problem**: pip is not installed or not in PATH.

**Solution**:
```bash
# Try using python3 -m pip instead
python3 -m pip install numpy scipy pandas scikit-learn

# Or install pip first
python3 -m ensurepip --upgrade
```

### Error: Permission denied when installing packages

**Problem**: Trying to install to system Python without permissions.

**Solution**:
```bash
# Option 1: Install for user only
pip install --user numpy scipy pandas scikit-learn

# Option 2: Use a virtual environment (recommended)
python3 -m venv ldpred2_env
source ldpred2_env/bin/activate  # On Windows: ldpred2_env\Scripts\activate
pip install -r requirements.txt
```

## Runtime Errors

### Error: `ValueError: df_beta must have column 'beta'`

**Problem**: Input DataFrame is missing required columns.

**Solution**: Ensure your summary statistics DataFrame has these columns:
```python
df_beta = pd.DataFrame({
    'beta': [...],      # Effect sizes
    'beta_se': [...],   # Standard errors
    'n_eff': [...]      # Effective sample size
})
```

### Error: `ValueError: All beta_se values must be positive`

**Problem**: Some standard errors are zero or negative.

**Solution**: Filter out invalid entries:
```python
df_beta = df_beta[df_beta['beta_se'] > 0]
```

### Error: `ValueError: Correlation matrix shape doesn't match`

**Problem**: LD matrix dimensions don't match number of variants in summary stats.

**Solution**: Ensure both have the same number of variants:
```python
m = len(df_beta)
assert corr.shape == (m, m), f"Expected {m}x{m}, got {corr.shape}"
```

### Warning: `Gibbs sampler diverged - returning NaN`

**Problem**: The Gibbs sampler detected divergence (instability).

**Solutions**:
1. Use shrinkage on the correlation matrix:
   ```python
   results = ldpred2_auto(corr, df_beta, shrink_corr=0.95)
   ```

2. Prevent sign jumps:
   ```python
   results = ldpred2_auto(corr, df_beta, allow_jump_sign=False)
   ```

3. Increase burn-in period:
   ```python
   results = ldpred2_auto(corr, df_beta, burn_in=1000)
   ```

4. Check for LD mismatch between reference panel and GWAS

### Error: `MemoryError`

**Problem**: Not enough RAM for the LD matrix.

**Solutions**:
1. Ensure LD matrix is in sparse format:
   ```python
   from scipy import sparse
   corr = sparse.csr_matrix(corr_dense)
   ```

2. Reduce the number of variants

3. Use a machine with more RAM

## Python Version Issues

### Python 3.12+ Compatibility

The code has been updated for Python 3.12+ compatibility. Key changes:
- Proper imports of `Union`, `Optional` from `typing`
- Compatible with stricter type checking

If you encounter type hint errors:
```bash
# Verify Python version
python3 --version

# Should be 3.7 or later
```

### Python 3.6 or earlier

**Problem**: Code requires Python 3.7+

**Solution**: Upgrade Python:
```bash
# Check current version
python3 --version

# Install Python 3.7+ from python.org or use conda
conda install python=3.10
```

## Import Errors

### Error: `ModuleNotFoundError: No module named 'ldpred2_python'`

**Problem**: Package not in Python path.

**Solutions**:
1. Run from the correct directory:
   ```bash
   cd /path/to/bigsnpr
   python3 examples/example_basic.py
   ```

2. Install the package:
   ```bash
   pip install -e .
   ```

3. Add to Python path:
   ```python
   import sys
   sys.path.insert(0, '/path/to/bigsnpr')
   from ldpred2_python import ldpred2_inf
   ```

### Error: `ImportError: attempted relative import with no known parent package`

**Problem**: Test file is being run directly with relative imports.

**Solution**: Use one of the provided test runners:

```bash
# Option 1: Use the main test runner (from project root)
python3 run_tests.py

# Option 2: Use the simple test runner
python3 ldpred2_python/tests/run_tests_simple.py

# Option 3: Use unittest discovery
python3 -m unittest discover ldpred2_python/tests

# Option 4: Use pytest (if installed)
pytest ldpred2_python/tests/
```

**Why this happens**: Python's import system requires proper package structure. Running `python3 ldpred2_python/tests/test_ldpred2.py` directly doesn't work because it tries to use relative imports without the parent package being available.

## Performance Issues

### Code is very slow

**Possible causes and solutions**:

1. **Large number of iterations**:
   - Reduce `num_iter` for testing
   - Use fewer parameter combinations in grid search

2. **Dense LD matrix**:
   - Convert to sparse format
   - Increase sparsity threshold when computing LD

3. **Large number of variants**:
   - This is expected - LDpred2 scales with m²
   - Consider using a subset for testing

## Data Issues

### Poor predictive performance

**Possible causes**:

1. **LD mismatch**: Reference panel ancestry differs from GWAS
   - Use `shrink_corr=0.95` or lower
   - Ensure reference panel matches GWAS ancestry

2. **Insufficient GWAS power**: Sample size too small
   - LDpred2 requires well-powered GWAS
   - Check N*h²/p is sufficiently large

3. **Wrong hyperparameters**:
   - Try LDpred2-auto instead of fixed parameters
   - Test wider range of p and h2 values

### LDpred2-auto estimates seem wrong

**Diagnostics**:

1. Check convergence:
   ```python
   import matplotlib.pyplot as plt
   plt.plot(results[0]['path_h2_est'])
   plt.title('h2 convergence')
   plt.show()
   ```

2. Run multiple chains:
   ```python
   results = ldpred2_auto(corr, df_beta, vec_p_init=[0.001, 0.01, 0.1, 0.5])
   ```

3. Check for divergence:
   ```python
   for i, res in enumerate(results):
       if np.any(np.isnan(res['beta_est'])):
           print(f"Chain {i+1} diverged")
   ```

## Getting Help

If you're still having issues:

1. **Check the documentation**:
   - README_MAIN.md
   - QUICKSTART.md
   - IMPLEMENTATION_NOTES.md

2. **Run the verification script**:
   ```bash
   python3 verify_installation.py
   ```

3. **Check the original LDpred2 resources**:
   - Tutorial: https://privefl.github.io/bigsnpr/articles/LDpred2.html
   - Paper: https://doi.org/10.1093/bioinformatics/btaa1029

4. **Common debugging steps**:
   ```python
   # Check data shapes
   print(f"LD matrix: {corr.shape}")
   print(f"Summary stats: {len(df_beta)}")
   
   # Check data types
   print(f"LD matrix type: {type(corr)}")
   print(f"LD matrix format: {corr.format if hasattr(corr, 'format') else 'N/A'}")
   
   # Check for NaN/Inf
   print(f"NaN in beta: {df_beta['beta'].isna().sum()}")
   print(f"Inf in beta: {np.isinf(df_beta['beta']).sum()}")
   ```

## Quick Fixes Checklist

- [ ] Python 3.7+ installed
- [ ] All dependencies installed (`pip install -r requirements.txt`)
- [ ] Running from correct directory
- [ ] Input data has required columns (beta, beta_se, n_eff)
- [ ] LD matrix is sparse format
- [ ] LD matrix dimensions match summary stats
- [ ] No NaN or Inf values in input data
- [ ] Sufficient memory available

