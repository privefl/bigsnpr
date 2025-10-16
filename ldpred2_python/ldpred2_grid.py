"""
LDpred2-grid: Gibbs sampler with fixed hyperparameters.

This module implements the Gibbs sampling algorithm for LDpred2 with
fixed hyperparameters (h2 and p).
"""

import numpy as np
import pandas as pd
from scipy import sparse
from typing import Union, List, Dict, Optional
import warnings

from .utils import standardize_betas, check_dimensions


def ldpred2_gibbs_one(corr: sparse.spmatrix,
                      beta_hat: np.ndarray,
                      n_vec: np.ndarray,
                      h2: float,
                      p: float,
                      sparse_output: bool = False,
                      burn_in: int = 50,
                      num_iter: int = 100,
                      seed: Optional[int] = None) -> np.ndarray:
    """
    Single Gibbs sampler run for LDpred2.
    
    Implements the core Gibbs sampling algorithm that iterates through variants,
    computing posterior probabilities and sampling effect sizes from their
    posterior distributions.
    
    Parameters
    ----------
    corr : sparse matrix
        Sparse LD correlation matrix (m x m)
    beta_hat : np.ndarray
        Standardized effect sizes (length m)
    n_vec : np.ndarray
        Sample sizes for each variant (length m)
    h2 : float
        Heritability
    p : float
        Proportion of causal variants
    sparse_output : bool
        If True, set effects to 0 if posterior probability < p
    burn_in : int
        Number of burn-in iterations
    num_iter : int
        Number of iterations after burn-in
    seed : int, optional
        Random seed for reproducibility
        
    Returns
    -------
    avg_beta : np.ndarray
        Average effect sizes across iterations (after burn-in)
        Returns NaN if divergence detected
    """
    if seed is not None:
        np.random.seed(seed)
    
    m = len(beta_hat)
    
    # Initialize current beta estimates
    curr_beta = np.zeros(m)
    avg_beta = np.zeros(m)
    
    # Dot products: corr @ curr_beta
    # We'll update this incrementally
    dotprods = np.zeros(m)
    
    # Precompute constants
    h2_per_var = h2 / (m * p)
    inv_odd_p = (1 - p) / p
    
    # Divergence detection threshold
    gap0 = 2 * np.sum(beta_hat ** 2)
    
    # Convert to CSR for efficient row access
    corr_csr = sparse.csr_matrix(corr)
    
    # Total iterations
    total_iter = burn_in + num_iter
    
    for k in range(-burn_in, num_iter):
        gap = 0.0
        
        for j in range(m):
            # Compute residual beta_hat for variant j
            # res_beta_hat_j = beta_hat[j] - (dotprods[j] - curr_beta[j])
            res_beta_hat_j = beta_hat[j] + curr_beta[j] - dotprods[j]
            
            # Compute posterior parameters
            C1 = h2_per_var * n_vec[j]
            C2 = 1.0 / (1.0 + 1.0 / C1)
            C3 = C2 * res_beta_hat_j
            C4 = C2 / n_vec[j]
            
            # Posterior probability of being causal
            post_p_j = 1.0 / (1.0 + inv_odd_p * np.sqrt(1.0 + C1) * 
                             np.exp(-C3 * C3 / C4 / 2.0))
            
            # Store old value for update
            old_beta = curr_beta[j]
            diff = -old_beta
            
            if sparse_output and (post_p_j < p):
                # Force to zero if posterior prob is low
                curr_beta[j] = 0.0
            else:
                # Sample from posterior
                if np.random.rand() < post_p_j:
                    # Sample from normal distribution
                    curr_beta[j] = np.random.normal(C3, np.sqrt(C4))
                    diff += curr_beta[j]
                    gap += curr_beta[j] ** 2
                else:
                    curr_beta[j] = 0.0
                
                # Accumulate average (only after burn-in)
                if k >= 0:
                    avg_beta[j] += C3 * post_p_j
            
            # Update dot products if beta changed
            if diff != 0:
                # dotprods += diff * corr[:, j]
                # For sparse matrix, extract column j
                col_j = corr_csr[:, j].toarray().flatten()
                dotprods += diff * col_j
        
        # Check for divergence
        if gap > gap0:
            warnings.warn("Gibbs sampler diverged - returning NaN")
            return np.full(m, np.nan)
    
    # Return average over iterations
    return avg_beta / num_iter


def ldpred2_grid(corr: sparse.spmatrix,
                df_beta: pd.DataFrame,
                grid_param: Union[Dict, List[Dict]],
                burn_in: int = 50,
                num_iter: int = 100,
                seed: Optional[int] = None) -> np.ndarray:
    """
    LDpred2-grid: Run Gibbs sampler over a grid of hyperparameters.
    
    Parameters
    ----------
    corr : sparse matrix
        Sparse LD correlation matrix
    df_beta : pd.DataFrame
        Summary statistics with columns: beta, beta_se, n_eff
    grid_param : dict or list of dict
        Grid of parameters. Each dict should have:
        - p: proportion of causal variants
        - h2: heritability
        - sparse: whether to use sparse thresholding
    burn_in : int
        Number of burn-in iterations
    num_iter : int
        Number of iterations after burn-in
    seed : int, optional
        Random seed for reproducibility
        
    Returns
    -------
    beta_grid : np.ndarray
        Matrix of effect sizes (m x n_params)
        Each column corresponds to one parameter combination
        
    Examples
    --------
    >>> import numpy as np
    >>> import pandas as pd
    >>> from scipy import sparse
    >>> from ldpred2_python import ldpred2_grid
    >>> 
    >>> # Create example data
    >>> m = 100
    >>> corr = sparse.eye(m, format='csr')
    >>> df_beta = pd.DataFrame({
    ...     'beta': np.random.randn(m) * 0.01,
    ...     'beta_se': np.ones(m) * 0.001,
    ...     'n_eff': np.ones(m) * 10000
    ... })
    >>> 
    >>> # Define parameter grid
    >>> params = [
    ...     {'p': 0.01, 'h2': 0.5, 'sparse': False},
    ...     {'p': 0.1, 'h2': 0.5, 'sparse': True}
    ... ]
    >>> 
    >>> # Run LDpred2-grid
    >>> beta_grid = ldpred2_grid(corr, df_beta, params, burn_in=10, num_iter=20)
    >>> print(beta_grid.shape)
    (100, 2)
    """
    check_dimensions(corr, df_beta)
    
    # Standardize effect sizes
    beta_hat, scale = standardize_betas(df_beta)
    N = df_beta['n_eff'].values
    
    # Convert single dict to list
    if isinstance(grid_param, dict):
        grid_param = [grid_param]
    
    m = len(df_beta)
    n_params = len(grid_param)
    
    # Initialize output matrix
    beta_grid = np.zeros((m, n_params))
    
    # Run Gibbs sampler for each parameter combination
    for i, params in enumerate(grid_param):
        p = params['p']
        h2 = params['h2']
        sparse_output = params.get('sparse', False)
        
        # Set seed for this run (for reproducibility)
        run_seed = None if seed is None else seed + i
        
        beta_scaled = ldpred2_gibbs_one(
            corr=corr,
            beta_hat=beta_hat,
            n_vec=N,
            h2=h2,
            p=p,
            sparse_output=sparse_output,
            burn_in=burn_in,
            num_iter=num_iter,
            seed=run_seed
        )
        
        # Convert back to original scale
        beta_grid[:, i] = beta_scaled * scale
    
    return beta_grid

