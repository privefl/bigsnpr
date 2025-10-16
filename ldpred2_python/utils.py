"""
Utility functions for LDpred2 implementation.
"""

import numpy as np
import pandas as pd
from scipy import sparse
from typing import Tuple


def validate_df_beta(df_beta: pd.DataFrame) -> None:
    """
    Validate that df_beta has required columns.
    
    Parameters
    ----------
    df_beta : pd.DataFrame
        Summary statistics dataframe
        
    Raises
    ------
    ValueError
        If required columns are missing or invalid
    """
    required_cols = ['beta', 'beta_se', 'n_eff']
    for col in required_cols:
        if col not in df_beta.columns:
            raise ValueError(f"df_beta must have column '{col}'")
    
    if not np.all(df_beta['beta_se'] > 0):
        raise ValueError("All beta_se values must be positive")
    
    if not np.all(df_beta['n_eff'] > 0):
        raise ValueError("All n_eff values must be positive")


def standardize_betas(df_beta: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    """
    Standardize effect sizes for LDpred2.
    
    Parameters
    ----------
    df_beta : pd.DataFrame
        Summary statistics with columns: beta, beta_se, n_eff
        
    Returns
    -------
    beta_hat : np.ndarray
        Standardized effect sizes
    scale : np.ndarray
        Scaling factors (to convert back to original scale)
    """
    validate_df_beta(df_beta)
    
    N = df_beta['n_eff'].values
    beta = df_beta['beta'].values
    beta_se = df_beta['beta_se'].values
    
    # Scale factor: sqrt(N * beta_se^2 + beta^2)
    scale = np.sqrt(N * beta_se**2 + beta**2)
    beta_hat = beta / scale
    
    return beta_hat, scale


def standardize_betas_auto(df_beta: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    """
    Standardize effect sizes for LDpred2-auto (different scaling).
    
    Parameters
    ----------
    df_beta : pd.DataFrame
        Summary statistics with columns: beta, beta_se, n_eff
        
    Returns
    -------
    beta_hat : np.ndarray
        Standardized effect sizes
    sd : np.ndarray
        Standard deviation (inverse of scale)
    """
    validate_df_beta(df_beta)
    
    N = df_beta['n_eff'].values
    beta = df_beta['beta'].values
    beta_se = df_beta['beta_se'].values
    
    # For auto: sd = 1 / sqrt(N * beta_se^2 + beta^2)
    sd = 1.0 / np.sqrt(N * beta_se**2 + beta**2)
    beta_hat = beta * sd
    
    return beta_hat, sd


def compute_ld_scores(corr: sparse.spmatrix) -> np.ndarray:
    """
    Compute LD scores from correlation matrix.
    
    LD score for variant j = sum_k(r_jk^2)
    
    Parameters
    ----------
    corr : sparse matrix
        Sparse LD correlation matrix
        
    Returns
    -------
    ld_scores : np.ndarray
        LD scores for each variant
    """
    # Convert to CSC format for efficient column operations
    corr_csc = sparse.csc_matrix(corr)
    
    # Compute sum of squared correlations for each column
    ld_scores = np.array(corr_csc.power(2).sum(axis=0)).flatten()
    
    return ld_scores


def load_ld_matrix(filepath: str) -> sparse.spmatrix:
    """
    Load sparse LD correlation matrix from file.
    
    Parameters
    ----------
    filepath : str
        Path to .npz file containing sparse matrix
        
    Returns
    -------
    corr : sparse matrix
        Sparse LD correlation matrix
    """
    loaded = sparse.load_npz(filepath)
    return loaded


def save_ld_matrix(corr: sparse.spmatrix, filepath: str) -> None:
    """
    Save sparse LD correlation matrix to file.
    
    Parameters
    ----------
    corr : sparse matrix
        Sparse LD correlation matrix
    filepath : str
        Path to save .npz file
    """
    sparse.save_npz(filepath, corr)


def check_dimensions(corr: sparse.spmatrix, df_beta: pd.DataFrame) -> None:
    """
    Check that dimensions match between correlation matrix and summary stats.
    
    Parameters
    ----------
    corr : sparse matrix
        LD correlation matrix
    df_beta : pd.DataFrame
        Summary statistics
        
    Raises
    ------
    ValueError
        If dimensions don't match
    """
    m = len(df_beta)
    if corr.shape[0] != m or corr.shape[1] != m:
        raise ValueError(
            f"Correlation matrix shape {corr.shape} doesn't match "
            f"number of variants {m}"
        )


def compute_mean_ld(corr: sparse.spmatrix) -> float:
    """
    Compute mean LD score.
    
    Parameters
    ----------
    corr : sparse matrix
        LD correlation matrix
        
    Returns
    -------
    mean_ld : float
        Mean LD score
    """
    ld_scores = compute_ld_scores(corr)
    return np.mean(ld_scores)


def sparse_solve_symmetric(A: sparse.spmatrix, b: np.ndarray, 
                           add_to_diag: float = 0.0) -> np.ndarray:
    """
    Solve symmetric sparse linear system (A + add_to_diag*I) x = b.
    
    Parameters
    ----------
    A : sparse matrix
        Symmetric sparse matrix
    b : np.ndarray
        Right-hand side vector
    add_to_diag : float
        Value to add to diagonal
        
    Returns
    -------
    x : np.ndarray
        Solution vector
    """
    from scipy.sparse.linalg import spsolve
    
    if add_to_diag != 0.0:
        # Add to diagonal
        n = A.shape[0]
        A_modified = A + sparse.eye(n) * add_to_diag
    else:
        A_modified = A
    
    # Solve the system
    x = spsolve(A_modified, b)
    
    return x


def create_grid_params(p_seq, h2_seq, sparse_seq):
    """
    Create grid of parameters for LDpred2-grid.
    
    Parameters
    ----------
    p_seq : list
        Sequence of p values (proportion causal)
    h2_seq : list
        Sequence of h2 values (heritability)
    sparse_seq : list
        Sequence of sparse flags (bool)
        
    Returns
    -------
    params : list of dict
        List of parameter combinations
    """
    import itertools
    
    params = []
    for p, h2, sparse in itertools.product(p_seq, h2_seq, sparse_seq):
        params.append({'p': p, 'h2': h2, 'sparse': sparse})
    
    return params


def seq_log(min_val: float, max_val: float, length: int) -> np.ndarray:
    """
    Create log-spaced sequence.
    
    Parameters
    ----------
    min_val : float
        Minimum value
    max_val : float
        Maximum value
    length : int
        Number of points
        
    Returns
    -------
    seq : np.ndarray
        Log-spaced sequence
    """
    return np.exp(np.linspace(np.log(min_val), np.log(max_val), length))

