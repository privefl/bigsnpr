#!/usr/bin/env python3
"""
Verify that LDpred2 Python implementation is correctly installed and working.

This script performs basic checks without requiring the full test suite.
"""

import sys
import os

def check_imports():
    """Check that all required packages can be imported."""
    print("Checking imports...")
    
    required = ['numpy', 'scipy', 'pandas']
    optional = ['matplotlib', 'sklearn']
    
    missing = []
    for pkg in required:
        try:
            __import__(pkg)
            print(f"  ✓ {pkg}")
        except ImportError:
            print(f"  ✗ {pkg} - MISSING")
            missing.append(pkg)
    
    for pkg in optional:
        try:
            __import__(pkg)
            print(f"  ✓ {pkg} (optional)")
        except ImportError:
            print(f"  ○ {pkg} (optional) - not installed")
    
    if missing:
        print(f"\nERROR: Missing required packages: {', '.join(missing)}")
        print("Install with: pip install -r requirements.txt")
        return False
    
    return True


def check_ldpred2_modules():
    """Check that LDpred2 modules can be imported."""
    print("\nChecking LDpred2 modules...")
    
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    
    modules = [
        'ldpred2_python.ldpred2_inf',
        'ldpred2_python.ldpred2_grid',
        'ldpred2_python.ldpred2_auto',
        'ldpred2_python.utils',
        'ldpred2_python.mle_alpha',
        'ldpred2_python.generate_test_data'
    ]
    
    for mod in modules:
        try:
            __import__(mod)
            print(f"  ✓ {mod}")
        except ImportError as e:
            print(f"  ✗ {mod} - ERROR: {e}")
            return False
    
    return True


def run_basic_test():
    """Run a basic functionality test."""
    print("\nRunning basic functionality test...")
    
    try:
        import numpy as np
        import pandas as pd
        from scipy import sparse
        
        from ldpred2_python import ldpred2_inf, ldpred2_grid, ldpred2_auto
        from ldpred2_python.generate_test_data import generate_complete_dataset
        
        # Generate small test dataset
        print("  Generating test data (m=50)...")
        data = generate_complete_dataset(
            m=50, n_gwas=1000, n_test=100,
            h2=0.3, p_causal=0.1, seed=42
        )
        
        corr = data['corr']
        df_beta = data['df_beta']
        
        # Test LDpred2-inf
        print("  Testing LDpred2-inf...")
        beta_inf = ldpred2_inf(corr, df_beta, h2=0.3)
        assert len(beta_inf) == 50
        assert np.all(np.isfinite(beta_inf))
        print("    ✓ LDpred2-inf works")
        
        # Test LDpred2-grid
        print("  Testing LDpred2-grid...")
        params = {'p': 0.1, 'h2': 0.3, 'sparse': False}
        beta_grid = ldpred2_grid(corr, df_beta, params, 
                                burn_in=5, num_iter=10, seed=42)
        assert beta_grid.shape == (50, 1)
        print("    ✓ LDpred2-grid works")
        
        # Test LDpred2-auto
        print("  Testing LDpred2-auto...")
        results = ldpred2_auto(corr, df_beta, h2_init=0.3,
                              burn_in=10, num_iter=10, seed=42)
        assert len(results) == 1
        assert 'beta_est' in results[0]
        assert len(results[0]['beta_est']) == 50
        print("    ✓ LDpred2-auto works")
        
        print("\n✓ All basic tests passed!")
        return True
        
    except Exception as e:
        print(f"\n✗ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Main verification routine."""
    print("=" * 70)
    print("LDpred2 Python Implementation - Installation Verification")
    print("=" * 70)
    
    # Check imports
    if not check_imports():
        print("\n" + "=" * 70)
        print("FAILED: Missing required packages")
        print("=" * 70)
        return 1
    
    # Check modules
    if not check_ldpred2_modules():
        print("\n" + "=" * 70)
        print("FAILED: Could not import LDpred2 modules")
        print("=" * 70)
        return 1
    
    # Run basic test
    if not run_basic_test():
        print("\n" + "=" * 70)
        print("FAILED: Basic functionality test failed")
        print("=" * 70)
        return 1
    
    # Success
    print("\n" + "=" * 70)
    print("SUCCESS: Installation verified!")
    print("=" * 70)
    print("\nNext steps:")
    print("  1. Run examples: python examples/example_basic.py")
    print("  2. Run full tests: python ldpred2_python/tests/test_ldpred2.py")
    print("  3. Read QUICKSTART.md for usage instructions")
    print("=" * 70)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

