#!/usr/bin/env python3
"""
Compare Python implementation with R reference implementation.

This script loads results from both implementations and compares them.

Usage:
    python3 validation/compare_with_r.py --r-dir path/to/r_results --py-dir path/to/py_results
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from pathlib import Path
import json


def compare_vectors(vec_r, vec_py, name, tolerance=1e-6):
    """Compare two vectors and return comparison metrics."""
    # Remove NaN values for comparison
    mask = ~(np.isnan(vec_r) | np.isnan(vec_py))
    vec_r_clean = vec_r[mask]
    vec_py_clean = vec_py[mask]
    
    if len(vec_r_clean) == 0:
        return {
            'name': name,
            'correlation': np.nan,
            'max_abs_diff': np.nan,
            'mean_abs_diff': np.nan,
            'rmse': np.nan,
            'passed': False,
            'message': 'All values are NaN'
        }
    
    correlation, _ = pearsonr(vec_r_clean, vec_py_clean)
    max_abs_diff = np.max(np.abs(vec_r_clean - vec_py_clean))
    mean_abs_diff = np.mean(np.abs(vec_r_clean - vec_py_clean))
    rmse = np.sqrt(np.mean((vec_r_clean - vec_py_clean) ** 2))
    
    # Determine if passed
    passed = correlation > 0.95 and max_abs_diff < tolerance
    
    return {
        'name': name,
        'correlation': correlation,
        'max_abs_diff': max_abs_diff,
        'mean_abs_diff': mean_abs_diff,
        'rmse': rmse,
        'passed': passed,
        'n_compared': len(vec_r_clean),
        'n_total': len(vec_r)
    }


def plot_comparison(vec_r, vec_py, name, output_path):
    """Create comparison plots."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    # Remove NaN for plotting
    mask = ~(np.isnan(vec_r) | np.isnan(vec_py))
    vec_r_clean = vec_r[mask]
    vec_py_clean = vec_py[mask]
    
    if len(vec_r_clean) == 0:
        plt.close()
        return
    
    # Scatter plot
    axes[0].scatter(vec_r_clean, vec_py_clean, alpha=0.5, s=10)
    lim_min = min(vec_r_clean.min(), vec_py_clean.min())
    lim_max = max(vec_r_clean.max(), vec_py_clean.max())
    axes[0].plot([lim_min, lim_max], [lim_min, lim_max], 'r--', lw=2)
    axes[0].set_xlabel('R implementation')
    axes[0].set_ylabel('Python implementation')
    corr = pearsonr(vec_r_clean, vec_py_clean)[0]
    axes[0].set_title(f'{name}\n(r = {corr:.6f})')
    axes[0].grid(True, alpha=0.3)
    
    # Difference histogram
    diff = vec_r_clean - vec_py_clean
    axes[1].hist(diff, bins=50, edgecolor='black', alpha=0.7)
    axes[1].axvline(0, color='r', linestyle='--', lw=2)
    axes[1].set_xlabel('Difference (R - Python)')
    axes[1].set_ylabel('Count')
    axes[1].set_title(f'Mean diff: {np.mean(diff):.2e}\nStd: {np.std(diff):.2e}')
    axes[1].grid(True, alpha=0.3)
    
    # Bland-Altman plot
    mean_val = (vec_r_clean + vec_py_clean) / 2
    axes[2].scatter(mean_val, diff, alpha=0.5, s=10)
    axes[2].axhline(0, color='r', linestyle='--', lw=2)
    axes[2].axhline(np.mean(diff), color='g', linestyle='--', lw=1, label='Mean diff')
    axes[2].axhline(np.mean(diff) + 1.96*np.std(diff), color='orange', 
                    linestyle='--', lw=1, label='±1.96 SD')
    axes[2].axhline(np.mean(diff) - 1.96*np.std(diff), color='orange', 
                    linestyle='--', lw=1)
    axes[2].set_xlabel('Mean of R and Python')
    axes[2].set_ylabel('Difference (R - Python)')
    axes[2].set_title('Bland-Altman Plot')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def compare_ldpred2_inf(r_dir, py_dir, output_dir):
    """Compare LDpred2-inf results."""
    print("\n" + "="*70)
    print("Comparing LDpred2-inf")
    print("="*70)
    
    # Load results
    beta_r = np.load(r_dir / "beta_inf_r.npy")
    beta_py = np.load(py_dir / "beta_inf_py.npy")
    
    # Compare
    result = compare_vectors(beta_r, beta_py, "LDpred2-inf", tolerance=1e-4)
    
    # Print results
    print(f"Correlation: {result['correlation']:.8f}")
    print(f"Max absolute difference: {result['max_abs_diff']:.2e}")
    print(f"Mean absolute difference: {result['mean_abs_diff']:.2e}")
    print(f"RMSE: {result['rmse']:.2e}")
    print(f"Status: {'✅ PASSED' if result['passed'] else '❌ FAILED'}")
    
    # Plot
    plot_comparison(beta_r, beta_py, "LDpred2-inf", 
                   output_dir / "comparison_ldpred2_inf.png")
    
    return result


def compare_ldpred2_grid(r_dir, py_dir, output_dir):
    """Compare LDpred2-grid results."""
    print("\n" + "="*70)
    print("Comparing LDpred2-grid")
    print("="*70)
    
    # Load results (assuming multiple parameter sets)
    beta_grid_r = np.load(r_dir / "beta_grid_r.npy")
    beta_grid_py = np.load(py_dir / "beta_grid_py.npy")
    
    results = []
    
    # Compare each parameter set
    n_params = beta_grid_r.shape[1] if beta_grid_r.ndim > 1 else 1
    
    if beta_grid_r.ndim == 1:
        beta_grid_r = beta_grid_r.reshape(-1, 1)
        beta_grid_py = beta_grid_py.reshape(-1, 1)
    
    for i in range(n_params):
        result = compare_vectors(beta_grid_r[:, i], beta_grid_py[:, i], 
                                f"LDpred2-grid (param set {i+1})", 
                                tolerance=0.01)
        results.append(result)
        
        print(f"\nParameter set {i+1}:")
        print(f"  Correlation: {result['correlation']:.6f}")
        print(f"  Max abs diff: {result['max_abs_diff']:.2e}")
        print(f"  Status: {'✅ PASSED' if result['passed'] else '❌ FAILED'}")
    
    # Plot first parameter set
    if n_params > 0:
        plot_comparison(beta_grid_r[:, 0], beta_grid_py[:, 0], 
                       "LDpred2-grid (param set 1)",
                       output_dir / "comparison_ldpred2_grid.png")
    
    return results


def compare_ldpred2_auto(r_dir, py_dir, output_dir):
    """Compare LDpred2-auto results."""
    print("\n" + "="*70)
    print("Comparing LDpred2-auto")
    print("="*70)
    
    # Load results
    import pickle
    
    with open(r_dir / "auto_r.pkl", "rb") as f:
        auto_r = pickle.load(f)
    
    with open(py_dir / "auto_py.pkl", "rb") as f:
        auto_py = pickle.load(f)
    
    results = []
    
    # Compare each chain
    n_chains = len(auto_py)
    
    for i in range(n_chains):
        print(f"\nChain {i+1}:")
        
        # Compare beta estimates
        beta_r = auto_r[i]['beta_est']
        beta_py = auto_py[i]['beta_est']
        
        result_beta = compare_vectors(beta_r, beta_py, 
                                     f"LDpred2-auto chain {i+1} (beta)",
                                     tolerance=0.05)
        
        print(f"  Beta correlation: {result_beta['correlation']:.6f}")
        
        # Compare parameter estimates
        h2_r = auto_r[i]['h2_est']
        h2_py = auto_py[i]['h2_est']
        p_r = auto_r[i]['p_est']
        p_py = auto_py[i]['p_est']
        alpha_r = auto_r[i].get('alpha_est', np.nan)
        alpha_py = auto_py[i].get('alpha_est', np.nan)
        
        print(f"  h² - R: {h2_r:.4f}, Python: {h2_py:.4f}, Diff: {abs(h2_r - h2_py):.4f}")
        print(f"  p  - R: {p_r:.4f}, Python: {p_py:.4f}, Diff: {abs(p_r - p_py):.4f}")
        if not np.isnan(alpha_r) and not np.isnan(alpha_py):
            print(f"  α  - R: {alpha_r:.4f}, Python: {alpha_py:.4f}, Diff: {abs(alpha_r - alpha_py):.4f}")
        
        # Check if parameters are close
        params_close = (abs(h2_r - h2_py) < 0.1 and 
                       abs(p_r - p_py) < 0.05)
        
        result_beta['params_close'] = params_close
        result_beta['h2_diff'] = abs(h2_r - h2_py)
        result_beta['p_diff'] = abs(p_r - p_py)
        
        results.append(result_beta)
        
        # Plot first chain
        if i == 0:
            plot_comparison(beta_r, beta_py, 
                           f"LDpred2-auto chain {i+1}",
                           output_dir / "comparison_ldpred2_auto.png")
    
    return results


def generate_report(all_results, output_path):
    """Generate a summary report."""
    with open(output_path, 'w') as f:
        f.write("# LDpred2 Python vs R Validation Report\n\n")
        
        # Summary
        all_passed = all(r.get('passed', False) for results in all_results.values() 
                        for r in (results if isinstance(results, list) else [results]))
        
        f.write("## Summary\n\n")
        if all_passed:
            f.write("✅ **All tests PASSED**\n\n")
        else:
            f.write("❌ **Some tests FAILED**\n\n")
        
        # Detailed results
        for method, results in all_results.items():
            f.write(f"## {method}\n\n")
            
            if not isinstance(results, list):
                results = [results]
            
            for result in results:
                f.write(f"### {result['name']}\n\n")
                f.write(f"- **Correlation**: {result['correlation']:.8f}\n")
                f.write(f"- **Max absolute difference**: {result['max_abs_diff']:.2e}\n")
                f.write(f"- **Mean absolute difference**: {result['mean_abs_diff']:.2e}\n")
                f.write(f"- **RMSE**: {result['rmse']:.2e}\n")
                f.write(f"- **Status**: {'✅ PASSED' if result['passed'] else '❌ FAILED'}\n")
                
                if 'h2_diff' in result:
                    f.write(f"- **h² difference**: {result['h2_diff']:.4f}\n")
                    f.write(f"- **p difference**: {result['p_diff']:.4f}\n")
                
                f.write("\n")
        
        # Interpretation
        f.write("## Interpretation\n\n")
        f.write("- **LDpred2-inf**: Should have correlation > 0.9999 (deterministic)\n")
        f.write("- **LDpred2-grid**: Should have correlation > 0.95 (stochastic)\n")
        f.write("- **LDpred2-auto**: Parameter estimates should be within tolerance\n")


def main():
    parser = argparse.ArgumentParser(description="Compare Python and R implementations")
    parser.add_argument("--r-dir", type=str, required=True,
                       help="Directory containing R results")
    parser.add_argument("--py-dir", type=str, required=True,
                       help="Directory containing Python results")
    parser.add_argument("--output-dir", type=str, default="validation_output",
                       help="Output directory for plots and reports")
    
    args = parser.parse_args()
    
    r_dir = Path(args.r_dir)
    py_dir = Path(args.py_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    print("="*70)
    print("LDpred2 Python vs R Validation")
    print("="*70)
    
    all_results = {}
    
    # Compare each method
    try:
        result_inf = compare_ldpred2_inf(r_dir, py_dir, output_dir)
        all_results['LDpred2-inf'] = result_inf
    except Exception as e:
        print(f"Error comparing LDpred2-inf: {e}")
    
    try:
        results_grid = compare_ldpred2_grid(r_dir, py_dir, output_dir)
        all_results['LDpred2-grid'] = results_grid
    except Exception as e:
        print(f"Error comparing LDpred2-grid: {e}")
    
    try:
        results_auto = compare_ldpred2_auto(r_dir, py_dir, output_dir)
        all_results['LDpred2-auto'] = results_auto
    except Exception as e:
        print(f"Error comparing LDpred2-auto: {e}")
    
    # Generate report
    generate_report(all_results, output_dir / "validation_report.md")
    
    # Save results as JSON
    with open(output_dir / "validation_results.json", 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print("\n" + "="*70)
    print("Validation complete!")
    print(f"Results saved to: {output_dir}")
    print("="*70)


if __name__ == "__main__":
    main()

