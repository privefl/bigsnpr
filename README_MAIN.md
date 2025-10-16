# LDpred2 Python Implementation - Complete Package

## 🎯 What Was Created

A **complete, standalone Python implementation** of the LDpred2 algorithm for computing polygenic scores from GWAS summary statistics. This implementation includes:

✅ All three LDpred2 variants (inf, grid, auto)  
✅ Comprehensive tests  
✅ Example scripts with synthetic data  
✅ Full documentation  
✅ Installation verification  

## 📁 Project Structure

```
.
├── ldpred2_python/              # Main package
│   ├── __init__.py              # Package exports
│   ├── ldpred2_inf.py          # Infinitesimal model
│   ├── ldpred2_grid.py         # Grid search with Gibbs sampler
│   ├── ldpred2_auto.py         # Automatic inference
│   ├── mle_alpha.py            # MLE optimization
│   ├── utils.py                # Utility functions
│   ├── generate_test_data.py  # Test data generation
│   ├── README.md               # Package documentation
│   └── tests/
│       └── test_ldpred2.py    # Unit tests
│
├── examples/
│   ├── example_basic.py        # Basic usage demo
│   └── example_workflow.py     # Complete workflow
│
├── Documentation/
│   ├── QUICKSTART.md           # Quick start guide
│   ├── IMPLEMENTATION_NOTES.md # Technical details
│   ├── PROJECT_SUMMARY.md      # Project overview
│   └── README_MAIN.md          # This file
│
├── requirements.txt             # Python dependencies
├── setup.py                     # Package installation
└── verify_installation.py       # Installation check
```

## 🚀 Quick Start

### 1. Install Dependencies

```bash
pip install numpy scipy pandas scikit-learn matplotlib
```

**Notes**:
- The code is fully compatible with Python 3.7 through 3.12
- If you encountered `NameError: name 'Union' is not defined`, this has been fixed (see `PYTHON312_FIX.md`)
- If you encountered `ImportError: attempted relative import`, this has been fixed (see `TEST_IMPORT_FIX.md`)
- Test failures due to Gibbs sampler divergence have been addressed (see `TEST_FIXES.md`)

### 2. Verify Installation

```bash
python verify_installation.py
```

### 3. Run Examples

```bash
# Basic example
python examples/example_basic.py

# Complete workflow
python examples/example_workflow.py
```

### 4. Run Tests

```bash
# Option 1: Use the test runner (recommended)
python3 run_tests.py

# Option 2: Use unittest discovery
python3 -m unittest discover ldpred2_python/tests

# Option 3: Use the simple test runner
python3 ldpred2_python/tests/run_tests_simple.py
```

## 💡 Usage Examples

### LDpred2-inf (Infinitesimal Model)

```python
from ldpred2_python import ldpred2_inf
import pandas as pd
from scipy import sparse

# Load data
corr = sparse.load_npz('ld_matrix.npz')
df_beta = pd.read_csv('sumstats.csv')  # Must have: beta, beta_se, n_eff

# Run
beta_inf = ldpred2_inf(corr, df_beta, h2=0.5)
```

### LDpred2-grid (Parameter Grid Search)

```python
from ldpred2_python import ldpred2_grid

params = [
    {'p': 0.01, 'h2': 0.5, 'sparse': False},
    {'p': 0.1, 'h2': 0.5, 'sparse': True}
]

beta_grid = ldpred2_grid(corr, df_beta, params, 
                        burn_in=50, num_iter=100)
```

### LDpred2-auto (Automatic Inference)

```python
from ldpred2_python import ldpred2_auto

results = ldpred2_auto(
    corr, df_beta,
    h2_init=0.5,
    vec_p_init=[0.01, 0.1, 0.5],  # Multiple chains
    burn_in=500,
    num_iter=200,
    sparse=True,
    shrink_corr=0.95
)

# Access results
for res in results:
    print(f"h2: {res['h2_est']:.3f}")
    print(f"p: {res['p_est']:.4f}")
    beta = res['beta_est']
```

## 📊 What Each Method Does

| Method | Description | When to Use |
|--------|-------------|-------------|
| **LDpred2-inf** | Assumes all variants causal | Fast, highly polygenic traits |
| **LDpred2-grid** | Tests parameter combinations | Have validation set for tuning |
| **LDpred2-auto** | Infers parameters automatically | No validation set, most flexible |

## 🔬 Algorithm Overview

LDpred2 computes polygenic scores by:

1. **Input**: GWAS summary statistics + LD correlation matrix
2. **Model**: Bayesian spike-and-slab prior on effect sizes
3. **Inference**: Gibbs sampling (grid/auto) or analytical solution (inf)
4. **Output**: Posterior effect size estimates

Key parameters:
- **h2**: SNP heritability (proportion of variance explained)
- **p**: Polygenicity (proportion of causal variants)
- **alpha**: MAF-effect size relationship

## 📚 Documentation Files

- **QUICKSTART.md**: Installation and basic usage
- **IMPLEMENTATION_NOTES.md**: Technical details and algorithms
- **PROJECT_SUMMARY.md**: Complete project overview
- **ldpred2_python/README.md**: Package-specific documentation

## ✅ Testing

The implementation includes:

1. **Unit tests** for all core functions
2. **Integration tests** for complete workflows
3. **Validation** on synthetic data with known parameters
4. **Reproducibility checks** with fixed random seeds

Run tests:
```bash
# Recommended: Use the test runner
python3 run_tests.py

# Alternative: Use unittest
python3 -m unittest discover ldpred2_python/tests
```

## 📦 Dependencies

**Required:**
- numpy >= 1.20.0
- scipy >= 1.7.0
- pandas >= 1.3.0
- scikit-learn >= 0.24.0

**Optional:**
- matplotlib >= 3.3.0 (for visualization examples)

## 🎓 References

**Primary Paper:**
Privé, F., Arbel, J., & Vilhjálmsson, B. J. (2020). LDpred2: better, faster, stronger. *Bioinformatics*, 36(22-23), 5424-5431.

**Original Implementation:**
- GitHub: https://github.com/privefl/bigsnpr
- Tutorial: https://privefl.github.io/bigsnpr/articles/LDpred2.html

## 🔍 Key Features

### ✓ Complete Implementation
- All three LDpred2 variants
- MLE optimization for alpha parameter
- Gibbs sampling with divergence detection
- Sparse matrix support

### ✓ Well-Tested
- Comprehensive unit tests
- Synthetic data validation
- Reproducibility verification
- Edge case handling

### ✓ Validated
- Validation framework against R reference implementation
- Automated comparison scripts
- Detailed validation guide
- See `validation/` directory and `VALIDATION_GUIDE.md`

### ✓ User-Friendly
- Simple API
- Clear documentation
- Example workflows
- Installation verification

### ✓ Efficient
- Sparse matrix operations
- Incremental updates in Gibbs sampler
- Optimized linear algebra

## 🚦 Getting Started Checklist

- [ ] Install dependencies: `pip install -r requirements.txt`
- [ ] Verify installation: `python verify_installation.py`
- [ ] Run basic example: `python examples/example_basic.py`
- [ ] Read QUICKSTART.md for detailed usage
- [ ] Run tests to ensure everything works
- [ ] Try with your own data

## 💻 System Requirements

- Python 3.7 or later
- ~500 MB RAM for small datasets (m=1000)
- Scales with number of variants and LD matrix sparsity

## 🐛 Troubleshooting

**Import errors?**
→ Install dependencies: `pip install -r requirements.txt`

**Tests failing?**
→ Check Python version (need 3.7+)

**Divergence in LDpred2-auto?**
→ Try `shrink_corr=0.95` or `allow_jump_sign=False`

**Memory issues?**
→ Ensure LD matrix is sparse format

## 📝 License

This implementation is provided for educational and research purposes. Please cite the original LDpred2 paper when using this code.

## 🙏 Acknowledgments

This implementation is based on the LDpred2 algorithm developed by Florian Privé and colleagues. All credit for the algorithm design goes to the original authors.

---

**Status**: ✅ Complete and ready to use  
**Version**: 0.1.0  
**Last Updated**: 2025  

For questions or issues, please refer to the documentation files or the original LDpred2 resources.

