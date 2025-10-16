"""
LDpred2-inf: Infinitesimal model implementation.

This module implements the infinitesimal model which assumes all variants
are causal (p = 1). It provides an analytical solution.
"""

import numpy as np
import pandas as pd
from scipy import sparse

from .utils import standardize_betas, check_dimensions, sparse_solve_symmetric


def ldpred2_inf(corr: sparse.spmatrix, 
                df_beta: pd.DataFrame, 
                h2: float) -> np.ndarray:
    """
    LDpred2 infinitesimal model.
    
    Assumes all variants are causal (p = 1) and provides an analytical solution
    by solving the linear system:
        beta_inf = (R + (m/(h2*N)) * I)^-1 * beta_hat
    
    where:
        R is the LD correlation matrix
        m is the number of variants
        h2 is the SNP heritability
        N is the (effective) sample size
        beta_hat are the standardized effect sizes
    
    Parameters
    ----------
    corr : sparse matrix
        Sparse LD correlation matrix (m x m)
    df_beta : pd.DataFrame
        Summary statistics with columns:
        - beta: effect size estimates
        - beta_se: standard errors
        - n_eff: effective sample size
    h2 : float
        Heritability estimate (must be positive)
        
    Returns
    -------
    beta_inf : np.ndarray
        Vector of effect sizes under infinitesimal model
        
    Raises
    ------
    ValueError
        If h2 is not positive or dimensions don't match
        
    Examples
    --------
    >>> import numpy as np
    >>> import pandas as pd
    >>> from scipy import sparse
    >>> from ldpred2_python import ldpred2_inf
    >>> 
    >>> # Create example data
    >>> m = 100
    >>> corr = sparse.eye(m, format='csr')  # Identity matrix (no LD)
    >>> df_beta = pd.DataFrame({
    ...     'beta': np.random.randn(m) * 0.01,
    ...     'beta_se': np.ones(m) * 0.001,
    ...     'n_eff': np.ones(m) * 10000
    ... })
    >>> 
    >>> # Run LDpred2-inf
    >>> beta_inf = ldpred2_inf(corr, df_beta, h2=0.5)
    >>> print(beta_inf.shape)
    (100,)
    
    References
    ----------
    Privé, F., Arbel, J., & Vilhjálmsson, B. J. (2020). 
    LDpred2: better, faster, stronger. 
    Bioinformatics, 36(22-23), 5424-5431.
    """
    # Validate inputs
    if h2 <= 0:
        raise ValueError("h2 must be positive")
    
    check_dimensions(corr, df_beta)
    
    # Standardize effect sizes
    beta_hat, scale = standardize_betas(df_beta)
    
    # Get sample sizes
    N = df_beta['n_eff'].values
    
    # Number of variants
    m = len(df_beta)
    
    # Compute diagonal addition: m / (h2 * N)
    # This is added to the diagonal of the correlation matrix
    add_to_diag = m / (h2 * N)
    
    # Solve the linear system for each variant
    # Since N can vary by variant, we need to handle this carefully
    # The system is: (R + diag(m/(h2*N))) * beta_inf = beta_hat
    
    if np.allclose(N, N[0]):
        # All sample sizes are the same - simpler case
        diag_add = m / (h2 * N[0])
        beta_inf_scaled = sparse_solve_symmetric(corr, beta_hat, add_to_diag=diag_add)
    else:
        # Sample sizes vary - need to add different values to diagonal
        diag_values = m / (h2 * N)
        corr_modified = corr + sparse.diags(diag_values, format='csr')
        from scipy.sparse.linalg import spsolve
        beta_inf_scaled = spsolve(corr_modified, beta_hat)
    
    # Convert back to original scale
    beta_inf = beta_inf_scaled * scale
    
    return beta_inf

