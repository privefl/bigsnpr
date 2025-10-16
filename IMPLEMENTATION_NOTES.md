# LDpred2 Python Implementation - Technical Notes

## Overview

This is a standalone Python implementation of the LDpred2 algorithm, based on the original R/C++ implementation by Florian Privé (https://github.com/privefl/bigsnpr).

## Implementation Details

### Architecture

The implementation is organized into the following modules:

1. **ldpred2_inf.py**: Infinitesimal model (analytical solution)
2. **ldpred2_grid.py**: Gibbs sampler with fixed hyperparameters
3. **ldpred2_auto.py**: Automatic hyperparameter inference
4. **mle_alpha.py**: Maximum likelihood estimation for alpha parameter
5. **utils.py**: Utility functions (standardization, LD scores, etc.)
6. **generate_test_data.py**: Synthetic data generation for testing

### Key Algorithms

#### 1. LDpred2-inf

Solves the linear system:
```
(R + (m/(h2*N)) * I) * beta_inf = beta_hat
```

Where:
- R: LD correlation matrix
- m: number of variants
- h2: heritability
- N: sample size
- beta_hat: standardized effect sizes

Implementation uses scipy's sparse linear solver for efficiency.

#### 2. LDpred2-grid (Gibbs Sampler)

For each iteration and each variant j:

1. Compute residual: `res_beta_hat_j = beta_hat[j] + curr_beta[j] - dotprods[j]`
2. Compute posterior parameters:
   - `C1 = (h2/(m*p)) * N[j]`
   - `C2 = 1 / (1 + 1/C1)`
   - `C3 = C2 * res_beta_hat_j`
   - `C4 = C2 / N[j]`
3. Compute posterior probability: `post_p = 1 / (1 + ((1-p)/p) * sqrt(1+C1) * exp(-C3^2/(2*C4)))`
4. Sample from posterior:
   - With probability `post_p`: sample from N(C3, C4)
   - Otherwise: set to 0
5. Update dot products incrementally

Divergence detection: if sum of squared betas exceeds 2*sum(beta_hat^2), return NaN.

#### 3. LDpred2-auto

Extends the Gibbs sampler to automatically infer:
- **h2** (heritability): Updated from current estimate
- **p** (polygenicity): Sampled from Beta distribution
- **alpha** (MAF-effect relationship): Estimated via MLE

The model assumes:
```
beta_j ~ N(0, sigma2 * exp(alpha * log_var_j))  with probability p
beta_j = 0                                       with probability 1-p
```

MLE optimization uses L-BFGS-B to minimize:
```
L = alpha * sum(log_var) + n * log(sigma2) + sum(beta^2 * exp(-alpha * log_var)) / sigma2
```

### Differences from R Implementation

1. **Sparse Matrix Handling**: Uses scipy.sparse instead of custom SFBM class
2. **Optimization**: Uses scipy.optimize instead of Roptim
3. **Random Number Generation**: Uses numpy.random instead of R's RNG
4. **Parallelization**: Not implemented (can be added using multiprocessing)

### Numerical Considerations

1. **Standardization**: Effect sizes are standardized differently for grid vs auto:
   - Grid: `scale = sqrt(N * beta_se^2 + beta^2)`, `beta_hat = beta / scale`
   - Auto: `sd = 1 / sqrt(N * beta_se^2 + beta^2)`, `beta_hat = beta * sd`

2. **Shrinkage**: The `shrink_corr` parameter (default 1.0, recommended 0.95) shrinks off-diagonal correlations to improve robustness

3. **Bounds**: Hyperparameters are constrained:
   - h2: minimum 1e-3
   - p: default bounds [1e-5, 1]
   - alpha: default bounds [-1.5, 0.5]

### Testing Strategy

Tests cover:
1. Utility functions (standardization, LD scores)
2. LDpred2-inf analytical solution
3. Gibbs sampler convergence and reproducibility
4. LDpred2-auto parameter inference
5. Data generation functions
6. Edge cases and error handling

### Performance Characteristics

**Time Complexity:**
- LDpred2-inf: O(m^2) for sparse solve (depends on sparsity)
- LDpred2-grid: O(iterations * m * nnz) where nnz is non-zeros in LD matrix
- LDpred2-auto: Similar to grid, plus MLE optimization per iteration

**Space Complexity:**
- O(m^2) for LD matrix (sparse storage)
- O(m) for effect sizes and intermediate values

**Typical Runtime** (on synthetic data, m=1000):
- LDpred2-inf: < 1 second
- LDpred2-grid (100 iterations): ~10-30 seconds
- LDpred2-auto (500 burn-in + 200 iterations): ~2-5 minutes

### Known Limitations

1. **No Parallelization**: Unlike R version, doesn't support multi-core processing
2. **Memory**: Entire LD matrix must fit in memory
3. **LD Matrix Format**: Requires scipy sparse format (not SFBM)
4. **Validation**: Less extensive testing than original implementation

### Future Improvements

Potential enhancements:
1. Add multiprocessing support for parameter grid
2. Implement chunked processing for very large datasets
3. Add more sophisticated convergence diagnostics
4. Support for additional LD matrix formats
5. GPU acceleration for matrix operations
6. Integration with PLINK/BGEN file formats

## Validation

The implementation has been validated against synthetic data where:
- True parameters are known
- Estimated parameters converge to true values
- Polygenic scores show expected correlations with phenotypes

For production use, we recommend:
1. Comparing results with the original R implementation on the same data
2. Running multiple chains and checking convergence
3. Validating on held-out test sets

## References

### Primary Reference
Privé, F., Arbel, J., & Vilhjálmsson, B. J. (2020). LDpred2: better, faster, stronger. *Bioinformatics*, 36(22-23), 5424-5431.

### Related Papers
- Privé, F., Albiñana, C., Arbel, J., Pasaniuc, B., & Vilhjálmsson, B. J. (2022). Inferring disease architecture and predictive ability with LDpred2-auto. *The American Journal of Human Genetics*, 109(12), 2042-2055.

- Privé, F., Albiñana, C., Pasaniuc, B., & Vilhjálmsson, B. J. (2023). Portability of 245 polygenic scores when derived from the UK Biobank and applied to 9 ancestry groups from the same cohort. *The American Journal of Human Genetics*, 110(1), 12-23.

### Original Implementation
- GitHub: https://github.com/privefl/bigsnpr
- Tutorial: https://privefl.github.io/bigsnpr/articles/LDpred2.html

## License

This implementation is provided for educational and research purposes. Please cite the original LDpred2 paper when using this code.

## Contact

For questions about the algorithm, please refer to the original LDpred2 documentation and papers.
For questions about this Python implementation, please open an issue on the repository.

