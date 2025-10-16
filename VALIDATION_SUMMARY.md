# Validation Against R Reference - Quick Summary

## How to Confirm Identical Performance

To validate that the Python implementation performs identically to the R reference implementation, follow this 3-step process:

### Step 1: Generate R Reference Results

```bash
Rscript validation/generate_r_reference.R
```

**What it does:**
- Creates synthetic test data (1000 variants, 5000 samples)
- Runs R's bigsnpr implementation (LDpred2-inf, grid, auto)
- Saves results to `validation/r_results/`

**Requirements:**
- R with bigsnpr package installed
- ~5-10 minutes runtime

### Step 2: Run Python Implementation

```bash
python3 validation/run_python_implementation.py
```

**What it does:**
- Loads the same test data
- Runs Python implementation
- Saves results to `validation/py_results/`

**Requirements:**
- Python with numpy, scipy, pandas
- ~5-10 minutes runtime

### Step 3: Compare Results

```bash
python3 validation/compare_with_r.py \
  --r-dir validation/r_results \
  --py-dir validation/py_results \
  --output-dir validation/comparison
```

**What it does:**
- Compares R and Python outputs
- Generates comparison plots
- Creates validation report

**Output:**
- `validation_report.md` - Summary with pass/fail
- Comparison plots for each method
- Detailed metrics in JSON format

## Expected Results

### ✅ LDpred2-inf (Deterministic)

```
Correlation: > 0.9999
Max difference: < 1e-5
Status: PASSED
```

**Why near-perfect:** No randomness, only numerical precision differences.

### ✅ LDpred2-grid (Stochastic)

```
Correlation: > 0.95
Prediction R² difference: < 0.01
Status: PASSED
```

**Why not perfect:** Different RNG algorithms in R and Python.

### ✅ LDpred2-auto (Stochastic + Optimization)

```
h² difference: < 0.1
p difference: < 0.05
Beta correlation: > 0.90
Status: PASSED
```

**Why not perfect:** Stochastic sampling + may converge to different local optima.

## What "Identical Performance" Means

### Functional Equivalence ✅

The implementations are **functionally equivalent** if:

1. **LDpred2-inf produces nearly identical results** (correlation > 0.9999)
   - This is deterministic, so should match almost exactly
   
2. **LDpred2-grid/auto produce similar prediction accuracy**
   - Same R² on held-out data (within 0.01)
   - Similar parameter estimates (within tolerance)
   
3. **Both handle edge cases similarly**
   - Divergence detection
   - Parameter bounds
   - Error handling

### NOT Perfect Numerical Identity ❌

We do NOT expect:
- Exact bit-for-bit identical outputs for stochastic methods
- Identical random number sequences (R and Python use different RNGs)
- Identical optimization paths (L-BFGS-B may differ slightly)

## Quick Validation Checklist

Run through this checklist to confirm validation:

- [ ] **Step 1 complete:** R reference results generated
- [ ] **Step 2 complete:** Python results generated  
- [ ] **Step 3 complete:** Comparison completed
- [ ] **LDpred2-inf:** Correlation > 0.9999 ✅
- [ ] **LDpred2-grid:** Correlation > 0.95 ✅
- [ ] **LDpred2-auto:** Parameter estimates within tolerance ✅
- [ ] **Plots generated:** Visual inspection shows agreement ✅
- [ ] **Report created:** validation_report.md shows PASSED ✅

## Alternative: Manual Spot Check

If you don't want to run the full validation, you can do a quick manual check:

```python
# Load R and Python results
import numpy as np

beta_r = np.load("validation/r_results/beta_inf_r.npy")
beta_py = np.load("validation/py_results/beta_inf_py.npy")

# Quick comparison
from scipy.stats import pearsonr
print(f"Correlation: {pearsonr(beta_r, beta_py)[0]:.8f}")
print(f"Max diff: {np.max(np.abs(beta_r - beta_py)):.2e}")

# Should see:
# Correlation: 0.99999xxx
# Max diff: x.xxe-06
```

## Validation on Real Data

For production use, validate on real GWAS data:

1. **Get public GWAS summary statistics**
   - E.g., from GWAS Catalog, UK Biobank, etc.

2. **Compute LD from reference panel**
   - Use 1000 Genomes or similar
   - Match ancestry to GWAS

3. **Run both implementations**
   - Same input data
   - Same parameters

4. **Compare prediction accuracy**
   - On held-out test samples
   - R² should be within 0.01

## Troubleshooting

### Large differences in LDpred2-inf

**Problem:** Correlation < 0.999

**Likely cause:** Bug in standardization or linear solver

**Fix:** Review `ldpred2_inf.py` and `utils.py`

### Large differences in LDpred2-grid/auto

**Problem:** Correlation < 0.90

**Likely cause:** Bug in Gibbs sampler or MLE

**Fix:** Review `ldpred2_grid.py`, `ldpred2_auto.py`, `mle_alpha.py`

### Divergence in Python but not R

**Problem:** Python returns NaN, R doesn't

**Likely cause:** Numerical stability issue

**Fix:** Add shrinkage (`shrink_corr=0.95`), check LD matrix

## Documentation

For detailed information, see:

- **`VALIDATION_GUIDE.md`** - Comprehensive validation guide
- **`validation/README.md`** - Validation scripts documentation
- **`IMPLEMENTATION_NOTES.md`** - Technical implementation details

## Conclusion

The Python implementation is validated against the R reference if:

1. ✅ LDpred2-inf: Near-perfect correlation (> 0.9999)
2. ✅ LDpred2-grid/auto: Similar prediction accuracy
3. ✅ Parameter estimates: Within reasonable tolerance
4. ✅ Edge cases: Handled similarly

**This demonstrates functional equivalence**, which is the appropriate standard given different RNG and optimization implementations.

Perfect numerical identity is neither expected nor required for validation.

