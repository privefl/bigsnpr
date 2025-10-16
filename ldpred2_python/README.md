# LDpred2 Python Implementation

A standalone Python implementation of the LDpred2 algorithm for computing polygenic scores from GWAS summary statistics.

## Overview

LDpred2 is a Bayesian method for deriving polygenic scores based on GWAS summary statistics and a reference panel for linkage disequilibrium (LD). This implementation provides three variants:

1. **LDpred2-inf**: Infinitesimal model assuming all variants are causal (analytical solution)
2. **LDpred2-grid**: Gibbs sampler with fixed hyperparameters (h2, p) 
3. **LDpred2-auto**: Automatic inference of hyperparameters (h2, p, alpha)

## Installation

```bash
pip install numpy scipy pandas scikit-learn
```

## Quick Start

```python
import numpy as np
import pandas as pd
from ldpred2_python import ldpred2_inf, ldpred2_grid, ldpred2_auto
from ldpred2_python.utils import load_ld_matrix, standardize_betas

# Load your data
corr = load_ld_matrix('ld_matrix.npz')  # Sparse LD correlation matrix
df_beta = pd.read_csv('sumstats.csv')   # GWAS summary statistics

# LDpred2-inf (infinitesimal model)
beta_inf = ldpred2_inf(corr, df_beta, h2=0.5)

# LDpred2-grid (with hyperparameter grid)
params = {'p': [0.001, 0.01, 0.1, 1.0], 'h2': [0.3, 0.5, 0.7], 'sparse': [False, True]}
beta_grid = ldpred2_grid(corr, df_beta, params, burn_in=50, num_iter=100)

# LDpred2-auto (automatic hyperparameter inference)
results = ldpred2_auto(corr, df_beta, h2_init=0.5, burn_in=500, num_iter=200)
beta_auto = results['beta_est']
```

## Input Data Format

### Summary Statistics (df_beta)
Required columns:
- `beta`: Effect size estimates
- `beta_se`: Standard errors of effect sizes
- `n_eff`: Effective sample size (or total sample size for continuous traits)

### LD Correlation Matrix (corr)
- Sparse correlation matrix (scipy.sparse format)
- Typically computed from a reference panel of the same ancestry as the GWAS
- Should match the variants in the summary statistics

## Algorithm Details

### LDpred2-inf
Assumes an infinitesimal model where all variants contribute to the trait. Solves:
```
beta_inf = (R + (m/(h2*N)) * I)^-1 * beta_hat
```
where R is the LD correlation matrix, m is the number of variants, h2 is heritability, and N is sample size.

### LDpred2-grid
Uses Gibbs sampling to estimate posterior effect sizes under a spike-and-slab prior:
- p: proportion of causal variants
- h2: SNP heritability
- Iterates through variants, sampling from posterior distributions

### LDpred2-auto
Extends LDpred2-grid by automatically inferring hyperparameters:
- p: polygenicity (proportion of causal variants)
- h2: SNP heritability
- alpha: parameter controlling relationship between MAF and effect sizes

## References

- Privé, F., Arbel, J., & Vilhjálmsson, B. J. (2020). LDpred2: better, faster, stronger. Bioinformatics, 36(22-23), 5424-5431.
- Original R implementation: https://github.com/privefl/bigsnpr

## License

This is an educational implementation based on the published LDpred2 algorithm.

