# LDpred2 Python Implementation - Project Summary

## Overview

This project provides a **standalone Python implementation of the LDpred2 algorithm** for computing polygenic scores from GWAS summary statistics. It is based on the original R/C++ implementation by Florian Privé (https://github.com/privefl/bigsnpr).

## What is LDpred2?

LDpred2 is a Bayesian method for deriving polygenic scores (PGS) that:
- Uses GWAS summary statistics (not individual-level data)
- Accounts for linkage disequilibrium (LD) between variants
- Can infer genetic architecture parameters (heritability, polygenicity)
- Provides state-of-the-art predictive performance

## Implementation Components

### Core Modules

1. **ldpred2_inf.py** - Infinitesimal model
   - Assumes all variants are causal (p=1)
   - Provides analytical solution via sparse linear system
   - Fastest method, works well for highly polygenic traits

2. **ldpred2_grid.py** - Gibbs sampler with fixed hyperparameters
   - Tests multiple (h2, p) combinations
   - Uses Gibbs sampling to estimate posterior effect sizes
   - Requires validation set to choose best parameters

3. **ldpred2_auto.py** - Automatic hyperparameter inference
   - Automatically infers h2, p, and alpha parameters
   - No validation set needed
   - Most sophisticated method

4. **mle_alpha.py** - Maximum likelihood estimation
   - Estimates alpha parameter (MAF-effect relationship)
   - Uses L-BFGS-B optimization with analytical gradients

5. **utils.py** - Utility functions
   - Effect size standardization
   - LD score computation
   - Sparse matrix operations
   - Input validation

6. **generate_test_data.py** - Synthetic data generation
   - Creates realistic test datasets
   - Generates LD matrices, GWAS summary stats, genotypes, phenotypes
   - Useful for testing and demonstrations

### Supporting Files

- **tests/test_ldpred2.py** - Comprehensive unit tests
- **examples/example_basic.py** - Basic usage demonstration
- **examples/example_workflow.py** - Complete workflow with visualization
- **verify_installation.py** - Installation verification script
- **requirements.txt** - Python dependencies
- **setup.py** - Package installation script
- **README.md** - Main documentation
- **QUICKSTART.md** - Quick start guide
- **IMPLEMENTATION_NOTES.md** - Technical details

## Key Features

### ✓ Complete Implementation
- All three LDpred2 variants (inf, grid, auto)
- MLE optimization for alpha parameter
- Gibbs sampling with divergence detection
- Sparse matrix support for efficiency

### ✓ Well-Documented
- Comprehensive docstrings for all functions
- Multiple example scripts
- Quick start guide
- Technical implementation notes

### ✓ Tested
- Unit tests for all major components
- Validation on synthetic data
- Reproducibility checks
- Edge case handling

### ✓ User-Friendly
- Simple API matching R implementation
- Clear error messages
- Installation verification script
- Example workflows

## File Structure

```
ldpred2_python/
├── __init__.py                 # Package initialization
├── ldpred2_inf.py             # Infinitesimal model
├── ldpred2_grid.py            # Grid search with Gibbs sampler
├── ldpred2_auto.py            # Automatic inference
├── mle_alpha.py               # MLE optimization
├── utils.py                   # Utility functions
├── generate_test_data.py      # Test data generation
├── README.md                  # Package documentation
└── tests/
    ├── __init__.py
    └── test_ldpred2.py        # Unit tests

examples/
├── example_basic.py           # Basic usage
└── example_workflow.py        # Complete workflow

Documentation:
├── QUICKSTART.md              # Quick start guide
├── IMPLEMENTATION_NOTES.md    # Technical details
├── PROJECT_SUMMARY.md         # This file
├── requirements.txt           # Dependencies
├── setup.py                   # Installation
└── verify_installation.py     # Verification script
```

## Usage Examples

### Minimal Example

```python
from ldpred2_python import ldpred2_inf
import pandas as pd
from scipy import sparse

# Load data
corr = sparse.load_npz('ld_matrix.npz')
df_beta = pd.read_csv('sumstats.csv')

# Run LDpred2-inf
beta = ldpred2_inf(corr, df_beta, h2=0.5)

# Compute polygenic scores
pgs = genotypes @ beta
```

### Complete Workflow

```python
from ldpred2_python import ldpred2_auto

# Run LDpred2-auto with multiple chains
results = ldpred2_auto(
    corr, df_beta,
    h2_init=0.5,
    vec_p_init=[0.01, 0.1, 0.5],
    burn_in=500,
    num_iter=200,
    sparse=True,
    shrink_corr=0.95,
    seed=42
)

# Check results
for i, res in enumerate(results):
    print(f"Chain {i+1}:")
    print(f"  h2 = {res['h2_est']:.3f}")
    print(f"  p = {res['p_est']:.4f}")
    print(f"  alpha = {res['alpha_est']:.3f}")
```

## Installation

```bash
# Install dependencies
pip install -r requirements.txt

# Verify installation
python verify_installation.py

# Run examples
python examples/example_basic.py
```

## Testing

```bash
# Run unit tests
python ldpred2_python/tests/test_ldpred2.py

# Or with pytest
pytest ldpred2_python/tests/
```

## Performance

Tested on synthetic data (m=1000 variants):
- **LDpred2-inf**: < 1 second
- **LDpred2-grid** (100 iterations): ~10-30 seconds  
- **LDpred2-auto** (500 burn-in + 200 iterations): ~2-5 minutes

## Validation

The implementation has been validated by:
1. Testing on synthetic data with known parameters
2. Verifying convergence of parameter estimates
3. Checking reproducibility with fixed seeds
4. Comparing algorithmic logic with original R/C++ code

## Limitations

- No multi-core parallelization (unlike R version)
- Requires entire LD matrix in memory
- Less extensively tested than original implementation
- No direct PLINK/BGEN file support

## Future Enhancements

Potential improvements:
- Add multiprocessing support
- Implement chunked processing for large datasets
- GPU acceleration
- Direct file format support (PLINK, BGEN)
- More convergence diagnostics

## Dependencies

- numpy >= 1.20.0
- scipy >= 1.7.0
- pandas >= 1.3.0
- scikit-learn >= 0.24.0
- matplotlib >= 3.3.0 (optional, for examples)

## References

**Primary Reference:**
Privé, F., Arbel, J., & Vilhjálmsson, B. J. (2020). LDpred2: better, faster, stronger. *Bioinformatics*, 36(22-23), 5424-5431.

**Related Papers:**
- Privé et al. (2022). Inferring disease architecture and predictive ability with LDpred2-auto. *AJHG*, 109(12), 2042-2055.
- Privé et al. (2023). Portability of 245 polygenic scores. *AJHG*, 110(1), 12-23.

**Original Implementation:**
- GitHub: https://github.com/privefl/bigsnpr
- Tutorial: https://privefl.github.io/bigsnpr/articles/LDpred2.html

## License

This implementation is provided for educational and research purposes. Please cite the original LDpred2 paper when using this code.

## Acknowledgments

This implementation is based on the LDpred2 algorithm developed by Florian Privé and colleagues. All credit for the algorithm design and original implementation goes to the original authors.

## Contact

For questions about the LDpred2 algorithm, please refer to the original documentation.
For questions about this Python implementation, please open an issue on the repository.

---

**Status**: Complete standalone Python implementation with tests, examples, and documentation.
**Date**: 2025
**Version**: 0.1.0

