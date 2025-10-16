#!/usr/bin/env Rscript
# Generate reference results using R bigsnpr implementation
#
# This script creates test data and runs the R implementation of LDpred2
# to generate reference results for validation.
#
# Usage:
#   Rscript validation/generate_r_reference.R

library(bigsnpr)
library(Matrix)

cat("="*70, "\n")
cat("Generating R Reference Results for LDpred2\n")
cat("="*70, "\n\n")

# Set seed for reproducibility
set.seed(42)

# Create output directory
dir.create("validation/r_results", showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Generate Test Data
# ============================================================================

cat("Generating test data...\n")

m <- 1000  # Number of variants
n <- 5000  # Sample size

# Generate genotype matrix (0, 1, 2)
cat("  - Genotype matrix (", n, "x", m, ")...\n", sep="")
G <- matrix(rbinom(n * m, 2, 0.3), n, m)

# Simulate phenotype
cat("  - Simulating phenotype...\n")
h2_true <- 0.5
p_causal_true <- 0.1
n_causal <- round(m * p_causal_true)

# Select causal variants
causal_idx <- sample(m, n_causal)
true_beta <- numeric(m)
true_beta[causal_idx] <- rnorm(n_causal, sd = sqrt(h2_true / n_causal))

# Generate phenotype
genetic_value <- G %*% true_beta
environmental_noise <- rnorm(n, sd = sqrt(1 - h2_true))
y <- genetic_value + environmental_noise

cat("  - True h²:", h2_true, "\n")
cat("  - True p:", p_causal_true, "\n")
cat("  - Number of causal variants:", n_causal, "\n")

# Compute GWAS summary statistics
cat("  - Computing GWAS summary statistics...\n")

# Simple linear regression for each variant
compute_sumstats <- function(G, y) {
  m <- ncol(G)
  n <- nrow(G)
  
  beta <- numeric(m)
  beta_se <- numeric(m)
  
  for (j in 1:m) {
    fit <- lm(y ~ G[, j])
    beta[j] <- coef(fit)[2]
    beta_se[j] <- summary(fit)$coefficients[2, 2]
  }
  
  data.frame(
    beta = beta,
    beta_se = beta_se,
    n_eff = n
  )
}

df_beta <- compute_sumstats(G, y)

# Compute LD correlation matrix
cat("  - Computing LD matrix...\n")
corr <- cor(G)

# Make sparse (threshold at 0.05)
cat("  - Making LD matrix sparse...\n")
corr[abs(corr) < 0.05] <- 0
corr_sparse <- Matrix(corr, sparse = TRUE)

cat("  - LD matrix sparsity:", 
    round(100 * (1 - nnzero(corr_sparse) / (m * m)), 2), "%\n")

# Save test data
cat("\nSaving test data...\n")
write.csv(df_beta, "validation/r_results/df_beta.csv", row.names = FALSE)
saveRDS(corr_sparse, "validation/r_results/corr_matrix.rds")
saveRDS(true_beta, "validation/r_results/true_beta.rds")
saveRDS(G, "validation/r_results/genotypes.rds")
saveRDS(y, "validation/r_results/phenotypes.rds")

# Also save as numpy-compatible format
library(reticulate)
np <- import("numpy", convert = FALSE)

# Save dense correlation matrix for Python
corr_dense <- as.matrix(corr_sparse)
np$save("validation/r_results/corr_matrix.npy", corr_dense)

cat("  ✓ Test data saved\n")

# ============================================================================
# Run LDpred2-inf
# ============================================================================

cat("\n")
cat("="*70, "\n")
cat("Running LDpred2-inf\n")
cat("="*70, "\n")

beta_inf <- snp_ldpred2_inf(corr_sparse, df_beta, h2 = 0.5)

cat("  - Number of estimates:", length(beta_inf), "\n")
cat("  - Range: [", min(beta_inf), ",", max(beta_inf), "]\n", sep="")

# Save results
saveRDS(beta_inf, "validation/r_results/beta_inf_r.rds")
np$save("validation/r_results/beta_inf_r.npy", beta_inf)

cat("  ✓ LDpred2-inf results saved\n")

# ============================================================================
# Run LDpred2-grid
# ============================================================================

cat("\n")
cat("="*70, "\n")
cat("Running LDpred2-grid\n")
cat("="*70, "\n")

# Define parameter grid
params_grid <- expand.grid(
  p = c(0.01, 0.1),
  h2 = c(0.3, 0.5),
  sparse = c(FALSE, TRUE)
)

cat("  - Testing", nrow(params_grid), "parameter combinations\n")
cat("  - Burn-in: 50, Iterations: 100\n")

beta_grid <- snp_ldpred2_grid(
  corr_sparse, 
  df_beta, 
  params_grid,
  burn_in = 50,
  num_iter = 100,
  ncores = 1
)

cat("  - Result dimensions:", nrow(beta_grid), "x", ncol(beta_grid), "\n")

# Save results
saveRDS(beta_grid, "validation/r_results/beta_grid_r.rds")
np$save("validation/r_results/beta_grid_r.npy", as.matrix(beta_grid))
write.csv(params_grid, "validation/r_results/params_grid.csv", row.names = FALSE)

cat("  ✓ LDpred2-grid results saved\n")

# ============================================================================
# Run LDpred2-auto
# ============================================================================

cat("\n")
cat("="*70, "\n")
cat("Running LDpred2-auto\n")
cat("="*70, "\n")

cat("  - Initial h²: 0.5\n")
cat("  - Initial p: [0.01, 0.1]\n")
cat("  - Burn-in: 500, Iterations: 200\n")
cat("  - Shrinkage: 0.95\n")

auto_result <- snp_ldpred2_auto(
  corr_sparse,
  df_beta,
  h2_init = 0.5,
  vec_p_init = c(0.01, 0.1),
  burn_in = 500,
  num_iter = 200,
  allow_jump_sign = FALSE,
  shrink_corr = 0.95,
  ncores = 1
)

cat("  - Number of chains:", length(auto_result), "\n")

# Print results for each chain
for (i in seq_along(auto_result)) {
  cat("\n  Chain", i, ":\n")
  cat("    h² estimate:", auto_result[[i]]$h2_est, "\n")
  cat("    p estimate:", auto_result[[i]]$p_est, "\n")
  cat("    α estimate:", auto_result[[i]]$alpha_est, "\n")
  
  # Check for divergence
  if (any(is.na(auto_result[[i]]$beta_est))) {
    cat("    ⚠ WARNING: Chain diverged (contains NaN)\n")
  }
}

# Save results
saveRDS(auto_result, "validation/r_results/auto_r.rds")

# Save as pickle for Python
library(reticulate)
pickle <- import("pickle")
py_auto <- r_to_py(auto_result)
pickle$dump(py_auto, file("validation/r_results/auto_r.pkl", "wb"))

cat("\n  ✓ LDpred2-auto results saved\n")

# ============================================================================
# Summary
# ============================================================================

cat("\n")
cat("="*70, "\n")
cat("Reference Data Generation Complete\n")
cat("="*70, "\n\n")

cat("Files saved to: validation/r_results/\n\n")

cat("Test data:\n")
cat("  - df_beta.csv (summary statistics)\n")
cat("  - corr_matrix.rds / .npy (LD matrix)\n")
cat("  - true_beta.rds (true effect sizes)\n")
cat("  - genotypes.rds (genotype matrix)\n")
cat("  - phenotypes.rds (phenotypes)\n\n")

cat("R results:\n")
cat("  - beta_inf_r.rds / .npy (LDpred2-inf)\n")
cat("  - beta_grid_r.rds / .npy (LDpred2-grid)\n")
cat("  - auto_r.rds / .pkl (LDpred2-auto)\n\n")

cat("Next steps:\n")
cat("  1. Run Python implementation: python3 validation/run_python_implementation.py\n")
cat("  2. Compare results: python3 validation/compare_with_r.py \\\n")
cat("       --r-dir validation/r_results --py-dir validation/py_results\n\n")

cat("="*70, "\n")

