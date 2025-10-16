# LDpred2 Python - Quick Start Guide

## Installation

### Prerequisites

Ensure you have Python 3.7 or later installed. Then install the required dependencies:

```bash
pip install -r requirements.txt
```

Or install individual packages:

```bash
pip install numpy scipy pandas scikit-learn matplotlib
```

**Note for Python 3.12+**: The code has been updated to be compatible with Python 3.12's stricter type hint requirements.

### Install the package (optional)

```bash
pip install -e .
```

### Verify Installation

```bash
python3 verify_installation.py
```

### Run Tests

To run the full test suite:

```bash
# Option 1: Use the test runner (recommended)
python3 run_tests.py

# Option 2: Use unittest discovery
python3 -m unittest discover ldpred2_python/tests
```

**Important**: Do not run the test file directly (`python3 ldpred2_python/tests/test_ldpred2.py`), as this will cause `ImportError: attempted relative import with no known parent package`. Always use one of the test runners above.

## Basic Usage

### 1. LDpred2-inf (Infinitesimal Model)

The simplest model that assumes all variants are causal:

```python
import pandas as pd
from scipy import sparse
from ldpred2_python import ldpred2_inf

# Load your data
corr = sparse.load_npz('ld_matrix.npz')  # LD correlation matrix
df_beta = pd.read_csv('sumstats.csv')    # GWAS summary statistics

# df_beta must have columns: beta, beta_se, n_eff
# Run LDpred2-inf
beta_inf = ldpred2_inf(corr, df_beta, h2=0.5)

# Compute polygenic scores
pgs = genotypes @ beta_inf
```

### 2. LDpred2-grid (Grid Search)

Test multiple hyperparameter combinations:

```python
from ldpred2_python import ldpred2_grid
from ldpred2_python.utils import create_grid_params

# Define parameter grid
p_seq = [0.001, 0.01, 0.1, 1.0]      # Proportion causal
h2_seq = [0.3, 0.5, 0.7]              # Heritability
sparse_seq = [False, True]            # Sparse thresholding

params = create_grid_params(p_seq, h2_seq, sparse_seq)

# Run LDpred2-grid
beta_grid = ldpred2_grid(
    corr, df_beta, params,
    burn_in=50,
    num_iter=100,
    seed=42
)

# beta_grid is a matrix: (n_variants x n_param_combinations)
# Evaluate each column on a validation set to choose best parameters
```

### 3. LDpred2-auto (Automatic Inference)

Automatically infer hyperparameters:

```python
from ldpred2_python import ldpred2_auto

# Run multiple chains for quality control
results = ldpred2_auto(
    corr, df_beta,
    h2_init=0.5,                    # Initial heritability
    vec_p_init=[0.01, 0.1, 0.5],   # Multiple initial p values
    burn_in=500,
    num_iter=200,
    sparse=True,                    # Also compute sparse solution
    shrink_corr=0.95,              # Recommended for robustness
    seed=42
)

# results is a list (one per p_init value)
for res in results:
    print(f"Estimated h2: {res['h2_est']:.3f}")
    print(f"Estimated p: {res['p_est']:.4f}")
    print(f"Estimated alpha: {res['alpha_est']:.3f}")
    
    # Use the estimated betas
    beta_auto = res['beta_est']
    pgs = genotypes @ beta_auto
```

## Running Examples

### Basic Example

```bash
cd examples
python example_basic.py
```

This demonstrates all three LDpred2 variants on synthetic data.

### Complete Workflow

```bash
cd examples
python example_workflow.py
```

This shows a realistic workflow including:
- Data loading and QC
- Running all LDpred2 methods
- Comparing performance
- Visualizing convergence

## Running Tests

```bash
cd ldpred2_python/tests
python test_ldpred2.py
```

Or using pytest:

```bash
pytest ldpred2_python/tests/
```

## Data Format Requirements

### Summary Statistics (df_beta)

Must be a pandas DataFrame with these columns:

- `beta`: Effect size estimates from GWAS
- `beta_se`: Standard errors of effect sizes
- `n_eff`: Effective sample size (or total N for continuous traits)

Example:
```python
df_beta = pd.DataFrame({
    'beta': [0.01, -0.02, 0.005, ...],
    'beta_se': [0.001, 0.001, 0.001, ...],
    'n_eff': [10000, 10000, 10000, ...]
})
```

### LD Correlation Matrix (corr)

Must be a scipy sparse matrix (CSR or CSC format) containing:
- Pairwise correlations between variants
- Typically computed from a reference panel
- Should match the variants in df_beta

Example:
```python
from scipy import sparse

# Load from file
corr = sparse.load_npz('ld_matrix.npz')

# Or create from dense matrix
corr_dense = np.corrcoef(genotypes.T)
corr = sparse.csr_matrix(corr_dense)
```

## Generating Test Data

```python
from ldpred2_python.generate_test_data import generate_complete_dataset

data = generate_complete_dataset(
    m=1000,           # Number of variants
    n_gwas=10000,     # GWAS sample size
    n_test=1000,      # Test sample size
    h2=0.5,           # True heritability
    p_causal=0.05,    # Proportion of causal variants
    seed=42
)

corr = data['corr']
df_beta = data['df_beta']
G_test = data['G_test']
y_test = data['y_test']
```

## Tips and Best Practices

### 1. Choosing the Right Method

- **LDpred2-inf**: Fast, works well when many variants are causal
- **LDpred2-grid**: Requires validation set, good when you can tune parameters
- **LDpred2-auto**: No validation set needed, automatically infers parameters

### 2. LDpred2-auto Recommendations

- Run multiple chains with different `p_init` values for QC
- Use `shrink_corr=0.95` for robustness (especially with LD mismatch)
- Check convergence by plotting `path_h2_est` and `path_p_est`
- Discard chains that diverge (contain NaN values)

### 3. Parameter Bounds

You can constrain parameter estimates:

```python
results = ldpred2_auto(
    corr, df_beta,
    h2_init=0.5,
    p_bounds=(1e-5, 0.5),        # Constrain p to [1e-5, 0.5]
    alpha_bounds=(-1.5, 0.5),    # Constrain alpha
    ...
)
```

### 4. Computational Considerations

- Larger `burn_in` and `num_iter` improve accuracy but take longer
- For LDpred2-grid, test fewer parameter combinations first
- Use `seed` parameter for reproducibility

## Troubleshooting

### Divergence Issues

If LDpred2-auto diverges (returns NaN):
- Try `shrink_corr=0.95` or lower
- Use `allow_jump_sign=False`
- Increase `burn_in`
- Check for LD mismatch between reference and GWAS

### Memory Issues

For very large datasets:
- Ensure LD matrix is sparse
- Process in chunks if needed
- Use sparse output: `sparse=True` in parameters

### Poor Performance

- Check that variants are properly matched between datasets
- Verify LD matrix quality
- Try different hyperparameter ranges
- Ensure sufficient GWAS sample size

## References

- Privé, F., Arbel, J., & Vilhjálmsson, B. J. (2020). LDpred2: better, faster, stronger. *Bioinformatics*, 36(22-23), 5424-5431.
- Original R implementation: https://github.com/privefl/bigsnpr
- Tutorial: https://privefl.github.io/bigsnpr/articles/LDpred2.html

