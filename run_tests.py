#!/usr/bin/env python3
"""
Test runner for LDpred2 Python implementation.

This script properly sets up the Python path and runs all tests.
"""

import sys
import os
import unittest

# Add the project root to Python path
project_root = os.path.dirname(os.path.abspath(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

def main():
    """Run all tests."""
    print("=" * 70)
    print("Running LDpred2 Python Tests")
    print("=" * 70)
    print()
    
    # Check if dependencies are installed
    try:
        import numpy
        import scipy
        import pandas
        print("✓ Dependencies found")
    except ImportError as e:
        print(f"✗ Missing dependency: {e}")
        print("\nPlease install dependencies first:")
        print("  pip install -r requirements.txt")
        return 1
    
    print()
    
    # Discover and run tests
    loader = unittest.TestLoader()
    start_dir = os.path.join(project_root, 'ldpred2_python', 'tests')
    suite = loader.discover(start_dir, pattern='test_*.py')
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print()
    print("=" * 70)
    if result.wasSuccessful():
        print("✓ All tests passed!")
        print("=" * 70)
        return 0
    else:
        print("✗ Some tests failed")
        print(f"  Failures: {len(result.failures)}")
        print(f"  Errors: {len(result.errors)}")
        print("=" * 70)
        return 1

if __name__ == '__main__':
    sys.exit(main())

