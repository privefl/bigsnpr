# LDpred2 Python - Quick Reference Card

## Installation

```bash
# Install dependencies
pip install numpy scipy pandas scikit-learn matplotlib

# Or use requirements file
pip install -r requirements.txt

# Or use the installation script
bash install_dependencies.sh

# Optional: Install package
pip install -e .
```

## Verification

```bash
# Verify installation
python3 verify_installation.py

# Check Python version
python3 --version  # Should be 3.7+

# Check dependencies
python3 -c "import numpy, scipy, pandas; print('OK')"
```

## Running Tests

```bash
# Recommended: Use test runner
python3 run_tests.py

# Alternative: Use unittest
python3 -m unittest discover ldpred2_python/tests

# Alternative: Simple runner
python3 ldpred2_python/tests/run_tests_simple.py

# With pytest (if installed)
pytest ldpred2_python/tests/
```

**❌ Don't do this:**
```bash
# This causes ImportError
python3 ldpred2_python/tests/test_ldpred2.py
```

## Running Examples

```bash
# Basic example (all three methods)
python3 examples/example_basic.py

# Complete workflow with visualization
python3 examples/example_workflow.py
```

## Basic Usage

### LDpred2-inf (Infinitesimal)

```python
from ldpred2_python import ldpred2_inf
import pandas as pd
from scipy import sparse

corr = sparse.load_npz('ld_matrix.npz')
df_beta = pd.read_csv('sumstats.csv')  # needs: beta, beta_se, n_eff

beta = ldpred2_inf(corr, df_beta, h2=0.5)
```

### LDpred2-grid (Parameter Search)

```python
from ldpred2_python import ldpred2_grid

params = [
    {'p': 0.01, 'h2': 0.5, 'sparse': False},
    {'p': 0.1, 'h2': 0.5, 'sparse': True}
]

beta_grid = ldpred2_grid(corr, df_beta, params, 
                         burn_in=50, num_iter=100)
```

### LDpred2-auto (Automatic)

```python
from ldpred2_python import ldpred2_auto

results = ldpred2_auto(
    corr, df_beta,
    h2_init=0.5,
    vec_p_init=[0.01, 0.1, 0.5],
    burn_in=500,
    num_iter=200,
    sparse=True
)

# Access results
beta = results[0]['beta_est']
h2 = results[0]['h2_est']
p = results[0]['p_est']
```

## Data Requirements

### Summary Statistics DataFrame

Must have these columns:
```python
df_beta = pd.DataFrame({
    'beta': [...],      # Effect sizes
    'beta_se': [...],   # Standard errors (must be > 0)
    'n_eff': [...]      # Effective sample size
})
```

### LD Correlation Matrix

```python
from scipy import sparse

# Should be sparse format
corr = sparse.csr_matrix(corr_dense)

# Or load from file
corr = sparse.load_npz('ld_matrix.npz')

# Must be square and match df_beta length
assert corr.shape == (len(df_beta), len(df_beta))
```

## Common Issues

### Import Errors

```bash
# ModuleNotFoundError: No module named 'numpy'
pip install numpy scipy pandas scikit-learn

# ImportError: attempted relative import
# Use test runners, don't run test file directly
python3 run_tests.py
```

### Type Hint Errors (Python 3.12)

```bash
# NameError: name 'Union' is not defined
# This has been fixed - update your files
# See PYTHON312_FIX.md for details
```

### Divergence Warnings

```python
# Use shrinkage
results = ldpred2_auto(corr, df_beta, shrink_corr=0.95)

# Prevent sign jumps
results = ldpred2_auto(corr, df_beta, allow_jump_sign=False)

# Increase burn-in
results = ldpred2_auto(corr, df_beta, burn_in=1000)
```

## File Structure

```
bigsnpr/
├── ldpred2_python/          # Main package
│   ├── ldpred2_inf.py       # Infinitesimal model
│   ├── ldpred2_grid.py      # Grid search
│   ├── ldpred2_auto.py      # Automatic inference
│   ├── mle_alpha.py         # MLE optimization
│   ├── utils.py             # Utilities
│   ├── generate_test_data.py # Test data
│   └── tests/
│       ├── test_ldpred2.py  # Unit tests
│       └── run_tests_simple.py
├── examples/
│   ├── example_basic.py     # Basic usage
│   └── example_workflow.py  # Complete workflow
├── run_tests.py             # Test runner
├── verify_installation.py   # Installation check
├── requirements.txt         # Dependencies
└── Documentation/
    ├── README_MAIN.md       # Main README
    ├── QUICKSTART.md        # Quick start
    ├── TROUBLESHOOTING.md   # Troubleshooting
    ├── PYTHON312_FIX.md     # Python 3.12 fix
    └── TEST_IMPORT_FIX.md   # Test import fix
```

## Key Parameters

### LDpred2-inf
- `h2`: SNP heritability (0 to 1)

### LDpred2-grid
- `p`: Proportion of causal variants (0 to 1)
- `h2`: SNP heritability (0 to 1)
- `sparse`: Use sparse prior (True/False)
- `burn_in`: Burn-in iterations (default: 50)
- `num_iter`: Sampling iterations (default: 100)

### LDpred2-auto
- `h2_init`: Initial h2 value
- `vec_p_init`: Initial p values (list for multiple chains)
- `burn_in`: Burn-in iterations (default: 500)
- `num_iter`: Sampling iterations (default: 200)
- `sparse`: Use sparse prior (default: False)
- `shrink_corr`: LD shrinkage (default: 1.0, try 0.95 if diverging)
- `allow_jump_sign`: Allow sign changes (default: True)

## Performance Tips

1. **Use sparse matrices**: `sparse.csr_matrix(corr)`
2. **Start with small iterations**: Test with `num_iter=10`
3. **Use sparse prior**: `sparse=True` for faster computation
4. **Reduce variants**: Test on subset first
5. **Multiple chains**: Use `vec_p_init=[0.01, 0.1, 0.5]`

## Validation Against R

To confirm identical performance with R reference:

```bash
# 1. Generate R reference (requires R + bigsnpr)
Rscript validation/generate_r_reference.R

# 2. Run Python implementation
python3 validation/run_python_implementation.py

# 3. Compare results
python3 validation/compare_with_r.py \
  --r-dir validation/r_results \
  --py-dir validation/py_results
```

See `VALIDATION_SUMMARY.md` for details.

## Getting Help

1. **Documentation**:
   - `README_MAIN.md` - Overview
   - `QUICKSTART.md` - Usage guide
   - `TROUBLESHOOTING.md` - Common issues
   - `IMPLEMENTATION_NOTES.md` - Technical details
   - `VALIDATION_GUIDE.md` - Validation against R

2. **Verification**:
   ```bash
   python3 verify_installation.py
   ```

3. **Original LDpred2**:
   - Tutorial: https://privefl.github.io/bigsnpr/articles/LDpred2.html
   - Paper: https://doi.org/10.1093/bioinformatics/btaa1029

## Quick Debugging

```python
# Check data
print(f"LD shape: {corr.shape}")
print(f"Summary stats: {len(df_beta)}")
print(f"Columns: {df_beta.columns.tolist()}")

# Check for issues
print(f"NaN in beta: {df_beta['beta'].isna().sum()}")
print(f"Zero beta_se: {(df_beta['beta_se'] <= 0).sum()}")
print(f"LD format: {corr.format if hasattr(corr, 'format') else 'dense'}")

# Check results
print(f"NaN in output: {np.isnan(beta).sum()}")
print(f"Inf in output: {np.isinf(beta).sum()}")
```

## Version Info

- **Python**: 3.7+ (tested up to 3.12)
- **NumPy**: >= 1.20.0
- **SciPy**: >= 1.7.0
- **Pandas**: >= 1.3.0
- **scikit-learn**: >= 0.24.0

## License & Citation

Please cite the original LDpred2 paper:

Privé, F., Arbel, J., & Vilhjálmsson, B. J. (2020). LDpred2: better, faster, stronger. *Bioinformatics*, 36(22-23), 5424-5431.

