# Validating Python Implementation Against R Reference

This guide explains how to confirm that the Python implementation performs identically (or near-identically) to the original R implementation in the bigsnpr package.

## Overview

Perfect numerical identity is difficult to achieve due to:
1. **Random number generators**: R and Python use different RNG algorithms
2. **Floating-point precision**: Minor differences in computation order
3. **Linear algebra libraries**: Different BLAS/LAPACK implementations
4. **Optimization algorithms**: L-BFGS-B implementations may differ slightly

However, we can validate that the implementations are **functionally equivalent** by comparing:
- Algorithm logic and flow
- Parameter estimates (within tolerance)
- Prediction accuracy on the same data
- Convergence behavior

## Validation Approaches

### Approach 1: Cross-Language Comparison (Recommended)

Compare outputs on the same input data with fixed random seeds.

#### Step 1: Prepare Shared Test Data

Create test data in R and export for Python:

```r
# In R
library(bigsnpr)
set.seed(42)

# Generate test data
m <- 1000
n <- 5000

# Create genotype matrix
G <- matrix(rbinom(n * m, 2, 0.3), n, m)

# Simulate phenotype
h2 <- 0.5
p_causal <- 0.1
n_causal <- round(m * p_causal)
causal_idx <- sample(m, n_causal)
true_beta <- numeric(m)
true_beta[causal_idx] <- rnorm(n_causal, sd = sqrt(h2 / n_causal))

y <- G %*% true_beta + rnorm(n, sd = sqrt(1 - h2))

# Compute GWAS summary statistics
gwas <- big_univLinReg(as_FBM(G), y)
df_beta <- data.frame(
  beta = gwas$estim,
  beta_se = gwas$std.err,
  n_eff = n
)

# Compute LD matrix
corr <- cor(G)
# Make sparse
corr[abs(corr) < 0.05] <- 0
corr_sparse <- Matrix::Matrix(corr, sparse = TRUE)

# Export data
write.csv(df_beta, "test_data/df_beta.csv", row.names = FALSE)
saveRDS(corr_sparse, "test_data/corr_matrix.rds")
saveRDS(true_beta, "test_data/true_beta.rds")
saveRDS(G, "test_data/genotypes.rds")
saveRDS(y, "test_data/phenotypes.rds")
```

#### Step 2: Run R Implementation

```r
# Run LDpred2-inf
beta_inf_r <- snp_ldpred2_inf(corr_sparse, df_beta, h2 = 0.5)

# Run LDpred2-grid
params_grid <- expand.grid(
  p = c(0.01, 0.1),
  h2 = c(0.3, 0.5),
  sparse = c(FALSE, TRUE)
)
beta_grid_r <- snp_ldpred2_grid(corr_sparse, df_beta, params_grid,
                                burn_in = 50, num_iter = 100)

# Run LDpred2-auto
auto_r <- snp_ldpred2_auto(corr_sparse, df_beta,
                           h2_init = 0.5,
                           vec_p_init = c(0.01, 0.1),
                           burn_in = 500,
                           num_iter = 200,
                           allow_jump_sign = FALSE,
                           shrink_corr = 0.95)

# Export results
saveRDS(beta_inf_r, "test_data/beta_inf_r.rds")
saveRDS(beta_grid_r, "test_data/beta_grid_r.rds")
saveRDS(auto_r, "test_data/auto_r.rds")
```

#### Step 3: Run Python Implementation

```python
import numpy as np
import pandas as pd
from scipy import sparse, io
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

# Load data exported from R
df_beta = pd.read_csv("test_data/df_beta.csv")

# Load sparse matrix from R
ro.r('corr <- readRDS("test_data/corr_matrix.rds")')
# Convert to Python sparse matrix (requires rpy2)
# Or manually convert via CSV

# Run Python implementations
from ldpred2_python import ldpred2_inf, ldpred2_grid, ldpred2_auto

# LDpred2-inf
beta_inf_py = ldpred2_inf(corr, df_beta, h2=0.5)

# LDpred2-grid
params = [
    {'p': 0.01, 'h2': 0.3, 'sparse': False},
    {'p': 0.01, 'h2': 0.3, 'sparse': True},
    {'p': 0.01, 'h2': 0.5, 'sparse': False},
    {'p': 0.01, 'h2': 0.5, 'sparse': True},
    {'p': 0.1, 'h2': 0.3, 'sparse': False},
    {'p': 0.1, 'h2': 0.3, 'sparse': True},
    {'p': 0.1, 'h2': 0.5, 'sparse': False},
    {'p': 0.1, 'h2': 0.5, 'sparse': True},
]
beta_grid_py = ldpred2_grid(corr, df_beta, params,
                            burn_in=50, num_iter=100, seed=42)

# LDpred2-auto
auto_py = ldpred2_auto(corr, df_beta,
                       h2_init=0.5,
                       vec_p_init=[0.01, 0.1],
                       burn_in=500,
                       num_iter=200,
                       allow_jump_sign=False,
                       shrink_corr=0.95,
                       seed=42)

# Save results
np.save("test_data/beta_inf_py.npy", beta_inf_py)
np.save("test_data/beta_grid_py.npy", beta_grid_py)
import pickle
with open("test_data/auto_py.pkl", "wb") as f:
    pickle.dump(auto_py, f)
```

#### Step 4: Compare Results

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Load R results (converted to numpy)
beta_inf_r = np.load("test_data/beta_inf_r.npy")
beta_inf_py = np.load("test_data/beta_inf_py.npy")

# Compare LDpred2-inf (should be nearly identical)
correlation = pearsonr(beta_inf_r, beta_inf_py)[0]
print(f"LDpred2-inf correlation: {correlation:.6f}")
print(f"Max absolute difference: {np.max(np.abs(beta_inf_r - beta_inf_py)):.6e}")
print(f"Mean absolute difference: {np.mean(np.abs(beta_inf_r - beta_inf_py)):.6e}")

# Plot comparison
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.scatter(beta_inf_r, beta_inf_py, alpha=0.5)
plt.plot([beta_inf_r.min(), beta_inf_r.max()], 
         [beta_inf_r.min(), beta_inf_r.max()], 'r--')
plt.xlabel("R implementation")
plt.ylabel("Python implementation")
plt.title(f"LDpred2-inf (r={correlation:.4f})")

plt.subplot(1, 2, 2)
plt.hist(beta_inf_r - beta_inf_py, bins=50)
plt.xlabel("Difference (R - Python)")
plt.ylabel("Count")
plt.title("Distribution of differences")
plt.tight_layout()
plt.savefig("validation_ldpred2_inf.png")
```

### Approach 2: Prediction Accuracy Comparison

Compare predictive performance on held-out data:

```python
# Load test genotypes and phenotypes
G_test = np.load("test_data/genotypes_test.npy")
y_test = np.load("test_data/phenotypes_test.npy")

# Compute PGS with both implementations
pgs_r = G_test @ beta_inf_r
pgs_py = G_test @ beta_inf_py

# Compare prediction accuracy
from sklearn.metrics import r2_score

r2_r = r2_score(y_test, pgs_r)
r2_py = r2_score(y_test, pgs_py)

print(f"R² (R implementation): {r2_r:.4f}")
print(f"R² (Python implementation): {r2_py:.4f}")
print(f"Difference: {abs(r2_r - r2_py):.6f}")
```

### Approach 3: Algorithm Component Testing

Test individual components in isolation:

#### Test 1: Standardization

```python
# Python
from ldpred2_python.utils import standardize_betas
beta_hat_py, scale_py = standardize_betas(df_beta)

# R
# beta_hat_r <- with(df_beta, beta / sqrt(n_eff * beta_se^2 + beta^2))
# scale_r <- sqrt(n_eff * beta_se^2 + beta^2)

# Compare
print(f"Standardization correlation: {pearsonr(beta_hat_r, beta_hat_py)[0]:.8f}")
```

#### Test 2: LD Scores

```python
from ldpred2_python.utils import compute_ld_scores

ld_scores_py = compute_ld_scores(corr)

# R: ld_scores_r <- Matrix::colSums(corr^2)

print(f"LD scores correlation: {pearsonr(ld_scores_r, ld_scores_py)[0]:.8f}")
```

#### Test 3: MLE Optimization

```python
from ldpred2_python.mle_alpha import mle_alpha

# Use same test data
result_py = mle_alpha(par_init, ind_causal, beta_hat, ld_scores, n_vec)

# Compare with R optim() results
print(f"Alpha estimate - R: {alpha_r:.4f}, Python: {result_py[0]:.4f}")
print(f"Sigma2 estimate - R: {sigma2_r:.4f}, Python: {result_py[1]:.4f}")
```

### Approach 4: Convergence Behavior

Compare convergence paths for LDpred2-auto:

```python
import matplotlib.pyplot as plt

# Extract convergence paths
h2_path_r = auto_r[[0]]['path_h2_est']
h2_path_py = auto_py[0]['path_h2_est']

plt.figure(figsize=(12, 4))

plt.subplot(1, 3, 1)
plt.plot(h2_path_r, label='R', alpha=0.7)
plt.plot(h2_path_py, label='Python', alpha=0.7)
plt.xlabel('Iteration')
plt.ylabel('h² estimate')
plt.legend()
plt.title('h² convergence')

plt.subplot(1, 3, 2)
plt.plot(auto_r[[0]]['path_p_est'], label='R', alpha=0.7)
plt.plot(auto_py[0]['path_p_est'], label='Python', alpha=0.7)
plt.xlabel('Iteration')
plt.ylabel('p estimate')
plt.legend()
plt.title('p convergence')

plt.subplot(1, 3, 3)
plt.plot(auto_r[[0]]['path_alpha_est'], label='R', alpha=0.7)
plt.plot(auto_py[0]['path_alpha_est'], label='Python', alpha=0.7)
plt.xlabel('Iteration')
plt.ylabel('α estimate')
plt.legend()
plt.title('α convergence')

plt.tight_layout()
plt.savefig("validation_convergence.png")
```

## Expected Results

### LDpred2-inf
- **Correlation**: > 0.9999 (should be nearly identical)
- **Max difference**: < 1e-6
- **Why**: Deterministic algorithm, no randomness

### LDpred2-grid
- **Correlation**: > 0.95 (may vary due to RNG)
- **Prediction R² difference**: < 0.01
- **Why**: Stochastic Gibbs sampling, different RNGs

### LDpred2-auto
- **Final h² difference**: < 0.05
- **Final p difference**: < 0.02
- **Correlation of beta estimates**: > 0.90
- **Why**: Stochastic + optimization, may converge to slightly different local optima

## Validation Checklist

- [ ] LDpred2-inf produces nearly identical results (correlation > 0.9999)
- [ ] LDpred2-grid produces similar prediction accuracy (R² difference < 0.01)
- [ ] LDpred2-auto converges to similar parameter estimates
- [ ] Standardization functions produce identical results
- [ ] LD score computation matches exactly
- [ ] MLE optimization produces similar estimates
- [ ] Both implementations handle edge cases similarly
- [ ] Both implementations detect divergence similarly
- [ ] Prediction accuracy on held-out data is comparable

## Automated Validation Script

I'll create a script that automates this comparison:

```bash
# Run validation
python3 validation/run_validation.py --r-results path/to/r_results --output validation_report.html
```

This will:
1. Load R and Python results
2. Compute all comparison metrics
3. Generate plots
4. Create an HTML report with pass/fail for each test

## Known Differences

### Acceptable Differences

1. **Random number generation**: Different RNG algorithms mean stochastic methods won't match exactly
2. **Floating-point order**: Minor differences in computation order can cause small numerical differences
3. **Optimization convergence**: L-BFGS-B may take different paths to similar solutions

### Unacceptable Differences

1. **LDpred2-inf differs significantly**: Indicates bug in linear solver or standardization
2. **Prediction accuracy differs substantially**: Indicates algorithmic error
3. **Divergence behavior differs**: May indicate numerical stability issues
4. **Parameter estimates consistently biased**: Indicates systematic error

## Troubleshooting

### If LDpred2-inf doesn't match:
- Check standardization formulas
- Verify sparse matrix format
- Check linear solver (try dense matrix)
- Verify diagonal addition formula

### If LDpred2-grid/auto don't converge similarly:
- Check Gibbs sampling loop logic
- Verify posterior probability calculation
- Check dot product update logic
- Verify MLE optimization

### If predictions differ:
- Ensure same genotype matrix
- Check effect size scaling
- Verify matrix multiplication order

## Conclusion

Perfect numerical identity is not expected due to RNG and implementation differences, but **functional equivalence** should be demonstrated through:

1. ✅ LDpred2-inf: Near-perfect correlation (> 0.9999)
2. ✅ LDpred2-grid/auto: Similar prediction accuracy
3. ✅ Parameter estimates: Within reasonable tolerance
4. ✅ Convergence behavior: Similar patterns

If these criteria are met, the Python implementation can be considered validated against the R reference.

