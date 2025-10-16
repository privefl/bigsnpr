"""
Setup script for LDpred2 Python implementation.
"""

from setuptools import setup, find_packages

with open("ldpred2_python/README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="ldpred2-python",
    version="0.1.0",
    author="LDpred2 Python Contributors",
    description="A standalone Python implementation of the LDpred2 algorithm",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "pandas>=1.3.0",
        "scikit-learn>=0.24.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "matplotlib>=3.3.0",
        ],
    },
)

