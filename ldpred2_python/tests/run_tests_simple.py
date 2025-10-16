#!/usr/bin/env python3
"""
Simple test runner that can be run from the tests directory.

Usage:
    cd ldpred2_python/tests
    python3 run_tests_simple.py
    
Or from project root:
    python3 ldpred2_python/tests/run_tests_simple.py
"""

import sys
import os

# Add project root to path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, '../..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# Now we can import and run the tests
if __name__ == '__main__':
    import unittest
    
    print("=" * 70)
    print("LDpred2 Python - Simple Test Runner")
    print("=" * 70)
    print()
    
    # Check dependencies
    try:
        import numpy
        import scipy
        import pandas
        print("✓ Dependencies installed")
    except ImportError as e:
        print(f"✗ Missing dependency: {e}")
        print("\nInstall with: pip install numpy scipy pandas scikit-learn")
        sys.exit(1)
    
    print()
    print("Running tests...")
    print("-" * 70)
    
    # Import the test module
    from test_ldpred2 import (TestUtils, TestLDpred2Inf, TestLDpred2Grid,
                              TestLDpred2Auto, TestDataGeneration)
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestUtils))
    suite.addTests(loader.loadTestsFromTestCase(TestLDpred2Inf))
    suite.addTests(loader.loadTestsFromTestCase(TestLDpred2Grid))
    suite.addTests(loader.loadTestsFromTestCase(TestLDpred2Auto))
    suite.addTests(loader.loadTestsFromTestCase(TestDataGeneration))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Summary
    print()
    print("=" * 70)
    if result.wasSuccessful():
        print("SUCCESS: All tests passed!")
        print(f"  Ran {result.testsRun} tests")
    else:
        print("FAILURE: Some tests failed")
        print(f"  Tests run: {result.testsRun}")
        print(f"  Failures: {len(result.failures)}")
        print(f"  Errors: {len(result.errors)}")
    print("=" * 70)
    
    sys.exit(0 if result.wasSuccessful() else 1)

