"""
Generate synthetic test data for LDpred2.

This module creates realistic synthetic data for testing and demonstration.
"""

import numpy as np
import pandas as pd
from scipy import sparse
from typing import Tuple, Dict


def generate_ld_matrix(m: int, block_size: int = 50,
                      within_block_corr: float = 0.5,
                      seed: int = None) -> sparse.spmatrix:
    """
    Generate a sparse LD correlation matrix with block structure.

    This creates a positive definite correlation matrix by using
    exponentially decaying correlations within blocks.

    Parameters
    ----------
    m : int
        Number of variants
    block_size : int
        Size of LD blocks
    within_block_corr : float
        Average correlation within blocks (should be < 0.5 for stability)
    seed : int, optional
        Random seed

    Returns
    -------
    corr : sparse matrix
        Sparse LD correlation matrix (positive definite)
    """
    if seed is not None:
        np.random.seed(seed)

    # Ensure within_block_corr is not too high to maintain positive definiteness
    within_block_corr = min(within_block_corr, 0.4)

    # Create block diagonal structure
    n_blocks = (m + block_size - 1) // block_size

    # Build correlation matrix
    row_ind = []
    col_ind = []
    data = []

    for block_idx in range(n_blocks):
        start = block_idx * block_size
        end = min(start + block_size, m)
        block_m = end - start

        # Add off-diagonal correlations that decay with distance
        for i in range(block_m):
            for j in range(i + 1, block_m):
                dist = j - i
                # Correlation decays exponentially with distance
                # Use smaller decay to ensure positive definiteness
                corr_val = within_block_corr * np.exp(-dist / 5.0)

                # Only add if correlation is significant
                if abs(corr_val) > 0.01:
                    # Add to sparse matrix
                    row_ind.append(start + i)
                    col_ind.append(start + j)
                    data.append(corr_val)

                    # Symmetric
                    row_ind.append(start + j)
                    col_ind.append(start + i)
                    data.append(corr_val)

        # Add diagonal
        for i in range(block_m):
            row_ind.append(start + i)
            col_ind.append(start + i)
            data.append(1.0)

    # Create sparse matrix
    corr = sparse.csr_matrix((data, (row_ind, col_ind)), shape=(m, m))

    return corr


def generate_gwas_sumstats(m: int, 
                          h2: float = 0.5,
                          p_causal: float = 0.01,
                          n_samples: int = 10000,
                          seed: int = None) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Generate synthetic GWAS summary statistics.
    
    Parameters
    ----------
    m : int
        Number of variants
    h2 : float
        True heritability
    p_causal : float
        Proportion of causal variants
    n_samples : int
        Sample size
    seed : int, optional
        Random seed
        
    Returns
    -------
    df_beta : pd.DataFrame
        Summary statistics with columns: beta, beta_se, n_eff
    true_beta : np.ndarray
        True effect sizes (for validation)
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Determine causal variants
    n_causal = int(m * p_causal)
    causal_idx = np.random.choice(m, size=n_causal, replace=False)
    
    # Generate true effect sizes
    true_beta = np.zeros(m)
    
    # Effect sizes for causal variants
    # Variance per causal variant
    var_per_causal = h2 / n_causal
    true_beta[causal_idx] = np.random.randn(n_causal) * np.sqrt(var_per_causal)
    
    # Generate observed effect sizes with noise
    # Standard error depends on sample size
    beta_se = np.ones(m) / np.sqrt(n_samples)
    
    # Add sampling noise
    observed_beta = true_beta + np.random.randn(m) * beta_se
    
    # Create dataframe
    df_beta = pd.DataFrame({
        'beta': observed_beta,
        'beta_se': beta_se,
        'n_eff': np.ones(m) * n_samples
    })
    
    return df_beta, true_beta


def generate_genotype_matrix(n_samples: int, m: int, 
                            maf_range: Tuple[float, float] = (0.05, 0.5),
                            seed: int = None) -> np.ndarray:
    """
    Generate synthetic genotype matrix.
    
    Parameters
    ----------
    n_samples : int
        Number of individuals
    m : int
        Number of variants
    maf_range : tuple
        Range of minor allele frequencies
    seed : int, optional
        Random seed
        
    Returns
    -------
    G : np.ndarray
        Genotype matrix (n_samples x m)
        Values are 0, 1, or 2 (number of minor alleles)
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Generate random MAFs
    mafs = np.random.uniform(maf_range[0], maf_range[1], size=m)
    
    # Generate genotypes
    G = np.zeros((n_samples, m), dtype=int)
    
    for j in range(m):
        # Sample genotypes based on Hardy-Weinberg equilibrium
        p = mafs[j]
        probs = [(1-p)**2, 2*p*(1-p), p**2]  # P(0), P(1), P(2)
        G[:, j] = np.random.choice([0, 1, 2], size=n_samples, p=probs)
    
    return G


def generate_phenotypes(G: np.ndarray, true_beta: np.ndarray,
                       h2: float = 0.5, seed: int = None) -> np.ndarray:
    """
    Generate phenotypes from genotypes and true effects.
    
    Parameters
    ----------
    G : np.ndarray
        Genotype matrix (n_samples x m)
    true_beta : np.ndarray
        True effect sizes (length m)
    h2 : float
        Heritability
    seed : int, optional
        Random seed
        
    Returns
    -------
    y : np.ndarray
        Phenotypes (length n_samples)
    """
    if seed is not None:
        np.random.seed(seed)
    
    n_samples = G.shape[0]
    
    # Genetic component
    y_genetic = G @ true_beta
    
    # Standardize genetic component
    y_genetic = (y_genetic - np.mean(y_genetic)) / np.std(y_genetic)
    
    # Environmental component
    var_env = (1 - h2) / h2
    y_env = np.random.randn(n_samples) * np.sqrt(var_env)
    
    # Total phenotype
    y = y_genetic + y_env
    
    return y


def generate_complete_dataset(m: int = 1000,
                             n_gwas: int = 10000,
                             n_test: int = 1000,
                             h2: float = 0.5,
                             p_causal: float = 0.01,
                             seed: int = None) -> Dict:
    """
    Generate a complete synthetic dataset for testing LDpred2.
    
    Parameters
    ----------
    m : int
        Number of variants
    n_gwas : int
        GWAS sample size
    n_test : int
        Test set sample size
    h2 : float
        Heritability
    p_causal : float
        Proportion of causal variants
    seed : int, optional
        Random seed
        
    Returns
    -------
    data : dict
        Dictionary containing:
        - corr: LD correlation matrix
        - df_beta: GWAS summary statistics
        - G_test: test genotype matrix
        - y_test: test phenotypes
        - true_beta: true effect sizes
        - true_h2: true heritability
        - true_p: true proportion causal
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Generate LD matrix
    print(f"Generating LD matrix for {m} variants...")
    corr = generate_ld_matrix(m, seed=seed)
    
    # Generate GWAS summary statistics
    print(f"Generating GWAS summary statistics (N={n_gwas})...")
    df_beta, true_beta = generate_gwas_sumstats(
        m, h2=h2, p_causal=p_causal, n_samples=n_gwas, seed=seed
    )
    
    # Generate test genotypes
    print(f"Generating test genotypes (N={n_test})...")
    G_test = generate_genotype_matrix(n_test, m, seed=seed)
    
    # Generate test phenotypes
    print("Generating test phenotypes...")
    y_test = generate_phenotypes(G_test, true_beta, h2=h2, seed=seed)
    
    data = {
        'corr': corr,
        'df_beta': df_beta,
        'G_test': G_test,
        'y_test': y_test,
        'true_beta': true_beta,
        'true_h2': h2,
        'true_p': p_causal
    }
    
    print("Dataset generation complete!")
    return data

