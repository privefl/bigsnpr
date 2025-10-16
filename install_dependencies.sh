#!/bin/bash
# Install dependencies for LDpred2 Python implementation

echo "Installing LDpred2 Python dependencies..."
echo "=========================================="

# Check Python version
python3 --version

echo ""
echo "Installing required packages..."
python3 -m pip install numpy scipy pandas scikit-learn

echo ""
echo "Installing optional packages..."
python3 -m pip install matplotlib

echo ""
echo "=========================================="
echo "Installation complete!"
echo ""
echo "To verify installation, run:"
echo "  python3 verify_installation.py"

