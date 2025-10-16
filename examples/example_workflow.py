"""
Complete workflow example for LDpred2.

This demonstrates a realistic workflow:
1. Load GWAS summary statistics
2. Load or compute LD matrix
3. Match variants between datasets
4. Run LDpred2 methods
5. Evaluate and compare results
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from ldpred2_python import ldpred2_inf, ldpred2_grid, ldpred2_auto
from ldpred2_python.generate_test_data import generate_complete_dataset
from ldpred2_python.utils import compute_ld_scores


def plot_convergence(results, save_path=None):
    """Plot convergence of LDpred2-auto."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    for i, res in enumerate(results):
        if np.any(np.isnan(res['beta_est'])):
            continue
        
        # Plot h2 path
        axes[0].plot(res['path_h2_est'], label=f"Chain {i+1}", alpha=0.7)
        
        # Plot p path
        axes[1].plot(res['path_p_est'], label=f"Chain {i+1}", alpha=0.7)
        
        # Plot alpha path
        if not np.all(np.isnan(res['path_alpha_est'])):
            axes[2].plot(res['path_alpha_est'], label=f"Chain {i+1}", alpha=0.7)
    
    axes[0].set_xlabel('Iteration')
    axes[0].set_ylabel('h² estimate')
    axes[0].set_title('Heritability convergence')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].set_xlabel('Iteration')
    axes[1].set_ylabel('p estimate')
    axes[1].set_title('Polygenicity convergence')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    axes[2].set_xlabel('Iteration')
    axes[2].set_ylabel('α estimate')
    axes[2].set_title('Alpha convergence')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Convergence plot saved to {save_path}")
    else:
        plt.show()
    
    plt.close()


def plot_grid_performance(params, performances, save_path=None):
    """Plot LDpred2-grid performance across parameters."""
    # Extract p and h2 values
    p_vals = [p['p'] for p in params]
    h2_vals = [p['h2'] for p in params]
    
    # Create scatter plot
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(p_vals, h2_vals, c=performances, 
                         s=100, cmap='viridis', alpha=0.7)
    plt.colorbar(scatter, label='R²')
    plt.xlabel('p (proportion causal)')
    plt.ylabel('h² (heritability)')
    plt.title('LDpred2-grid performance across parameters')
    plt.xscale('log')
    plt.grid(True, alpha=0.3)
    
    # Mark best point
    best_idx = np.argmax(performances)
    plt.scatter([p_vals[best_idx]], [h2_vals[best_idx]], 
               marker='*', s=500, c='red', edgecolors='black', 
               label='Best', zorder=10)
    plt.legend()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Grid performance plot saved to {save_path}")
    else:
        plt.show()
    
    plt.close()


def main():
    print("=" * 70)
    print("LDpred2 Complete Workflow Example")
    print("=" * 70)
    
    # Step 1: Generate/load data
    print("\nStep 1: Loading data...")
    print("-" * 70)
    
    data = generate_complete_dataset(
        m=1000,
        n_gwas=15000,
        n_test=2000,
        h2=0.4,
        p_causal=0.02,
        seed=123
    )
    
    corr = data['corr']
    df_beta = data['df_beta']
    G_test = data['G_test']
    y_test = data['y_test']
    
    print(f"Loaded {len(df_beta)} variants")
    print(f"GWAS sample size: {df_beta['n_eff'].iloc[0]:.0f}")
    print(f"Test sample size: {len(y_test)}")
    
    # Step 2: Compute LD scores (for QC)
    print("\nStep 2: Computing LD scores...")
    print("-" * 70)
    
    ld_scores = compute_ld_scores(corr)
    print(f"Mean LD score: {np.mean(ld_scores):.2f}")
    print(f"Median LD score: {np.median(ld_scores):.2f}")
    
    # Step 3: Run LDpred2-inf
    print("\nStep 3: Running LDpred2-inf...")
    print("-" * 70)
    
    beta_inf = ldpred2_inf(corr, df_beta, h2=0.4)
    pgs_inf = G_test @ beta_inf
    r2_inf = np.corrcoef(pgs_inf, y_test)[0, 1] ** 2
    print(f"LDpred2-inf R²: {r2_inf:.4f}")
    
    # Step 4: Run LDpred2-grid
    print("\nStep 4: Running LDpred2-grid...")
    print("-" * 70)
    
    # Define grid
    from ldpred2_python.utils import seq_log
    p_seq = seq_log(1e-4, 1, 10)
    h2_seq = [0.2, 0.4, 0.6]
    
    params = []
    for p in p_seq:
        for h2 in h2_seq:
            params.append({'p': p, 'h2': h2, 'sparse': False})
    
    print(f"Testing {len(params)} parameter combinations...")
    
    beta_grid = ldpred2_grid(corr, df_beta, params, 
                            burn_in=50, num_iter=100, seed=123)
    
    # Evaluate all
    performances = []
    for i in range(beta_grid.shape[1]):
        pgs = G_test @ beta_grid[:, i]
        r2 = np.corrcoef(pgs, y_test)[0, 1] ** 2
        performances.append(r2)
    
    best_idx = np.argmax(performances)
    best_r2_grid = performances[best_idx]
    best_params = params[best_idx]
    
    print(f"Best parameters: p={best_params['p']:.4f}, h2={best_params['h2']:.2f}")
    print(f"Best R²: {best_r2_grid:.4f}")
    
    # Plot grid performance
    plot_grid_performance(params, performances, 
                         save_path='ldpred2_grid_performance.png')
    
    # Step 5: Run LDpred2-auto
    print("\nStep 5: Running LDpred2-auto...")
    print("-" * 70)
    
    results = ldpred2_auto(
        corr, df_beta,
        h2_init=0.5,
        vec_p_init=[0.001, 0.01, 0.1],
        burn_in=300,
        num_iter=200,
        sparse=True,
        shrink_corr=0.95,
        seed=123
    )
    
    # Evaluate chains
    print("\nChain results:")
    best_r2_auto = 0
    best_chain = None
    
    for i, res in enumerate(results):
        if np.any(np.isnan(res['beta_est'])):
            print(f"  Chain {i+1}: DIVERGED")
            continue
        
        pgs = G_test @ res['beta_est']
        r2 = np.corrcoef(pgs, y_test)[0, 1] ** 2
        
        print(f"  Chain {i+1}: h2={res['h2_est']:.3f}, p={res['p_est']:.4f}, "
              f"alpha={res['alpha_est']:.3f}, R²={r2:.4f}")
        
        if r2 > best_r2_auto:
            best_r2_auto = r2
            best_chain = res
    
    # Plot convergence
    plot_convergence(results, save_path='ldpred2_auto_convergence.png')
    
    # Step 6: Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nTrue parameters:")
    print(f"  h² = {data['true_h2']:.3f}")
    print(f"  p = {data['true_p']:.4f}")
    
    print(f"\nPredictive performance (R²):")
    print(f"  LDpred2-inf:  {r2_inf:.4f}")
    print(f"  LDpred2-grid: {best_r2_grid:.4f}")
    print(f"  LDpred2-auto: {best_r2_auto:.4f}")
    
    if best_chain is not None:
        print(f"\nBest LDpred2-auto estimates:")
        print(f"  h² = {best_chain['h2_est']:.3f}")
        print(f"  p = {best_chain['p_est']:.4f}")
        print(f"  α = {best_chain['alpha_est']:.3f}")
    
    print("\n" + "=" * 70)
    print("Workflow completed successfully!")
    print("=" * 70)


if __name__ == '__main__':
    main()

