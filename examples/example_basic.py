"""
Basic example of using LDpred2 Python implementation.

This script demonstrates:
1. Generating synthetic test data
2. Running LDpred2-inf (infinitesimal model)
3. Running LDpred2-grid with parameter grid
4. Running LDpred2-auto for automatic inference
5. Evaluating polygenic scores
"""

import numpy as np
import pandas as pd
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from ldpred2_python import ldpred2_inf, ldpred2_grid, ldpred2_auto
from ldpred2_python.generate_test_data import generate_complete_dataset
from ldpred2_python.utils import create_grid_params, seq_log


def compute_pgs(G, beta):
    """Compute polygenic scores."""
    return G @ beta


def evaluate_pgs(pgs, y):
    """Evaluate polygenic score performance."""
    # Correlation
    corr = np.corrcoef(pgs, y)[0, 1]
    
    # R-squared
    r2 = corr ** 2
    
    return {'correlation': corr, 'r2': r2}


def main():
    print("=" * 70)
    print("LDpred2 Python Implementation - Basic Example")
    print("=" * 70)
    
    # Generate synthetic data
    print("\n1. Generating synthetic data...")
    print("-" * 70)
    
    data = generate_complete_dataset(
        m=500,           # 500 variants
        n_gwas=10000,    # GWAS sample size
        n_test=1000,     # Test sample size
        h2=0.5,          # True heritability
        p_causal=0.05,   # 5% causal variants
        seed=42
    )
    
    corr = data['corr']
    df_beta = data['df_beta']
    G_test = data['G_test']
    y_test = data['y_test']
    true_beta = data['true_beta']
    
    print(f"Number of variants: {len(df_beta)}")
    print(f"True heritability: {data['true_h2']}")
    print(f"True proportion causal: {data['true_p']}")
    print(f"Test sample size: {len(y_test)}")
    
    # LDpred2-inf
    print("\n2. Running LDpred2-inf (infinitesimal model)...")
    print("-" * 70)
    
    beta_inf = ldpred2_inf(corr, df_beta, h2=0.5)
    pgs_inf = compute_pgs(G_test, beta_inf)
    perf_inf = evaluate_pgs(pgs_inf, y_test)
    
    print(f"LDpred2-inf correlation: {perf_inf['correlation']:.4f}")
    print(f"LDpred2-inf R²: {perf_inf['r2']:.4f}")
    
    # LDpred2-grid
    print("\n3. Running LDpred2-grid with parameter grid...")
    print("-" * 70)
    
    # Create parameter grid
    p_seq = [0.001, 0.01, 0.05, 0.1, 0.5, 1.0]
    h2_seq = [0.3, 0.5, 0.7]
    sparse_seq = [False, True]
    
    params = create_grid_params(p_seq, h2_seq, sparse_seq)
    print(f"Testing {len(params)} parameter combinations...")
    
    beta_grid = ldpred2_grid(
        corr, df_beta, params,
        burn_in=50,
        num_iter=100,
        seed=42
    )
    
    # Evaluate all parameter combinations
    best_r2 = 0
    best_idx = 0
    
    for i in range(beta_grid.shape[1]):
        pgs = compute_pgs(G_test, beta_grid[:, i])
        perf = evaluate_pgs(pgs, y_test)
        
        if perf['r2'] > best_r2:
            best_r2 = perf['r2']
            best_idx = i
    
    best_params = params[best_idx]
    print(f"\nBest parameters: p={best_params['p']}, h2={best_params['h2']}, "
          f"sparse={best_params['sparse']}")
    print(f"Best LDpred2-grid R²: {best_r2:.4f}")
    
    # LDpred2-auto
    print("\n4. Running LDpred2-auto (automatic inference)...")
    print("-" * 70)
    
    # Run multiple chains for quality control
    results = ldpred2_auto(
        corr, df_beta,
        h2_init=0.5,
        vec_p_init=[0.01, 0.1, 0.5],  # Multiple initial values
        burn_in=200,
        num_iter=100,
        sparse=True,
        shrink_corr=0.95,  # Recommended for robustness
        seed=42
    )
    
    print(f"Ran {len(results)} chains")
    
    # Evaluate each chain
    for i, res in enumerate(results):
        if np.any(np.isnan(res['beta_est'])):
            print(f"Chain {i+1}: DIVERGED")
            continue
        
        pgs = compute_pgs(G_test, res['beta_est'])
        perf = evaluate_pgs(pgs, y_test)
        
        print(f"\nChain {i+1} (p_init={res['p_init']}):")
        print(f"  Estimated h2: {res['h2_est']:.4f} (true: {data['true_h2']:.4f})")
        print(f"  Estimated p: {res['p_est']:.4f} (true: {data['true_p']:.4f})")
        print(f"  Estimated alpha: {res['alpha_est']:.4f}")
        print(f"  Correlation: {perf['correlation']:.4f}")
        print(f"  R²: {perf['r2']:.4f}")
        
        # Also evaluate sparse solution
        if 'beta_est_sparse' in res:
            pgs_sparse = compute_pgs(G_test, res['beta_est_sparse'])
            perf_sparse = evaluate_pgs(pgs_sparse, y_test)
            n_nonzero = np.sum(res['beta_est_sparse'] != 0)
            print(f"  Sparse solution: {n_nonzero} non-zero effects, R²={perf_sparse['r2']:.4f}")
    
    # Summary
    print("\n5. Summary")
    print("=" * 70)
    print(f"True heritability: {data['true_h2']:.4f}")
    print(f"True proportion causal: {data['true_p']:.4f}")
    print(f"\nLDpred2-inf R²: {perf_inf['r2']:.4f}")
    print(f"LDpred2-grid best R²: {best_r2:.4f}")
    
    # Find best auto result
    best_auto_r2 = 0
    for res in results:
        if not np.any(np.isnan(res['beta_est'])):
            pgs = compute_pgs(G_test, res['beta_est'])
            perf = evaluate_pgs(pgs, y_test)
            if perf['r2'] > best_auto_r2:
                best_auto_r2 = perf['r2']
    
    print(f"LDpred2-auto best R²: {best_auto_r2:.4f}")
    
    print("\n" + "=" * 70)
    print("Example completed successfully!")
    print("=" * 70)


if __name__ == '__main__':
    main()

