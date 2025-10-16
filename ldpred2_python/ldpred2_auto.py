"""
LDpred2-auto: Automatic hyperparameter inference.

This module implements LDpred2-auto which automatically infers the
hyperparameters h2, p, and alpha during Gibbs sampling.
"""

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import beta as beta_dist
from typing import Dict, List, Optional, Union
import warnings

from .utils import standardize_betas_auto, check_dimensions, compute_mean_ld
from .mle_alpha import mle_alpha
from .ldpred2_grid import ldpred2_gibbs_one


MIN_H2 = 1e-3


def ldpred2_gibbs_auto(corr: sparse.spmatrix,
                       beta_hat: np.ndarray,
                       n_vec: np.ndarray,
                       log_var: np.ndarray,
                       h2_init: float,
                       p_init: float,
                       burn_in: int = 500,
                       num_iter: int = 200,
                       shrink_corr: float = 1.0,
                       use_mle: bool = True,
                       p_bounds: tuple = (1e-5, 1.0),
                       alpha_bounds: tuple = (-1.5, 0.5),
                       mean_ld: Optional[float] = None,
                       allow_jump_sign: bool = True,
                       seed: Optional[int] = None) -> Dict:
    """
    LDpred2-auto Gibbs sampler with automatic hyperparameter inference.
    
    Parameters
    ----------
    corr : sparse matrix
        Sparse LD correlation matrix
    beta_hat : np.ndarray
        Standardized effect sizes
    n_vec : np.ndarray
        Sample sizes
    log_var : np.ndarray
        Log variance (related to allele frequencies)
    h2_init : float
        Initial heritability
    p_init : float
        Initial proportion of causal variants
    burn_in : int
        Number of burn-in iterations
    num_iter : int
        Number of iterations after burn-in
    shrink_corr : float
        Shrinkage factor for correlation matrix (0.95 recommended for robustness)
    use_mle : bool
        Whether to use MLE for alpha estimation
    p_bounds : tuple
        Bounds for p parameter
    alpha_bounds : tuple
        Bounds for alpha parameter
    mean_ld : float, optional
        Mean LD score (computed if not provided)
    allow_jump_sign : bool
        Whether to allow effect signs to change between iterations
    seed : int, optional
        Random seed
        
    Returns
    -------
    results : dict
        Dictionary containing:
        - beta_est: estimated effect sizes
        - postp_est: posterior probabilities
        - corr_est: estimated correlations
        - path_p_est: path of p estimates
        - path_h2_est: path of h2 estimates
        - path_alpha_est: path of alpha estimates
    """
    if seed is not None:
        np.random.seed(seed)
    
    m = len(beta_hat)
    
    # Compute mean LD if not provided
    if mean_ld is None:
        mean_ld = compute_mean_ld(corr)
    
    # Initialize
    curr_beta = np.zeros(m)
    avg_beta = np.zeros(m)
    avg_postp = np.zeros(m)
    avg_beta_hat = np.zeros(m)
    
    dotprods = np.zeros(m)
    
    # Initialize hyperparameters
    h2 = max(h2_init, MIN_H2)
    p = np.clip(p_init, p_bounds[0], p_bounds[1])
    
    # MLE parameters: [alpha_plus_one, sigma2]
    par_mle = np.array([0.0, h2 / (m * p)])
    
    # Paths for tracking convergence
    num_iter_tot = burn_in + num_iter
    p_est = np.full(num_iter_tot, np.nan)
    h2_est = np.full(num_iter_tot, np.nan)
    alpha_est = np.full(num_iter_tot, np.nan)
    
    # Divergence detection
    gap0 = 2 * np.sum(beta_hat ** 2)
    
    # Convert to CSR for efficient operations
    corr_csr = sparse.csr_matrix(corr)
    
    cur_h2_est = 0.0
    
    for k in range(num_iter_tot):
        inv_odd_p = (1 - p) / p
        alpha_plus_one = par_mle[0]
        sigma2 = par_mle[1]
        gap = 0.0
        
        ind_causal = []
        
        for j in range(m):
            # Compute residual
            dotprod = dotprods[j]
            res_beta_hat_j = beta_hat[j] - shrink_corr * (dotprod - curr_beta[j])
            
            # Scale by frequency if using MLE
            if use_mle:
                scale_freq = np.exp(alpha_plus_one * log_var[j])
            else:
                scale_freq = 1.0
            
            # Posterior parameters
            C1 = scale_freq * sigma2 * n_vec[j]
            C2 = 1.0 / (1.0 + 1.0 / C1)
            C3 = C2 * res_beta_hat_j
            C4 = C2 / n_vec[j]
            
            # Posterior probability
            postp = 1.0 / (1.0 + inv_odd_p * np.sqrt(1.0 + C1) * 
                          np.exp(-C3 * C3 / C4 / 2.0))
            
            prev_beta = curr_beta[j]
            dotprod_shrunk = shrink_corr * dotprod + (1 - shrink_corr) * prev_beta
            
            # Accumulate averages (after burn-in)
            if k >= burn_in:
                avg_postp[j] += postp
                avg_beta[j] += C3 * postp
                avg_beta_hat[j] += dotprod_shrunk
            
            # Sample new beta
            diff = -prev_beta
            if np.random.rand() < postp:
                samp_beta = np.random.normal(C3, np.sqrt(C4))
                
                # Check sign jump if not allowed
                if (not allow_jump_sign) and (samp_beta * prev_beta < 0):
                    curr_beta[j] = 0.0
                else:
                    curr_beta[j] = samp_beta
                    diff += samp_beta
                    ind_causal.append(j)
                    gap += samp_beta ** 2
            else:
                curr_beta[j] = 0.0
            
            # Update dot products
            if diff != 0:
                cur_h2_est += diff * (2 * dotprod_shrunk + diff)
                col_j = corr_csr[:, j].toarray().flatten()
                dotprods += diff * col_j
        
        # Check divergence
        if gap > gap0:
            warnings.warn("LDpred2-auto diverged")
            avg_beta[:] = np.nan
            avg_postp[:] = np.nan
            avg_beta_hat[:] = np.nan
            break
        
        # Update hyperparameters
        nb_causal = len(ind_causal)
        
        # Update p using beta distribution
        p = beta_dist.rvs(1 + nb_causal / mean_ld, 1 + (m - nb_causal) / mean_ld)
        p = np.clip(p, p_bounds[0], p_bounds[1])
        
        # Update h2
        h2 = max(cur_h2_est, MIN_H2)
        
        # Update alpha and sigma2 using MLE
        if use_mle and len(ind_causal) > 0:
            try:
                par_mle = mle_alpha(
                    par_init=par_mle,
                    ind_causal=ind_causal,
                    log_var=log_var,
                    curr_beta=curr_beta,
                    alpha_bounds=alpha_bounds,
                    bootstrap=True,
                    verbose=False
                )
            except:
                # If MLE fails, use simple update
                par_mle[1] = h2 / (m * p)
        else:
            par_mle[1] = h2 / (m * p)
        
        # Store paths
        p_est[k] = p
        h2_est[k] = h2
        if use_mle:
            alpha_est[k] = par_mle[0] - 1.0
    
    # Compute final estimates
    results = {
        'beta_est': avg_beta / num_iter,
        'postp_est': avg_postp / num_iter,
        'corr_est': avg_beta_hat / num_iter,
        'path_p_est': p_est,
        'path_h2_est': h2_est,
        'path_alpha_est': alpha_est
    }
    
    return results


def ldpred2_auto(corr: sparse.spmatrix,
                df_beta: pd.DataFrame,
                h2_init: float,
                vec_p_init: Union[float, List[float]] = 0.1,
                burn_in: int = 500,
                num_iter: int = 200,
                sparse: bool = False,
                shrink_corr: float = 1.0,
                use_mle: bool = True,
                p_bounds: tuple = (1e-5, 1.0),
                alpha_bounds: tuple = (-1.5, 0.5),
                allow_jump_sign: bool = True,
                seed: Optional[int] = None) -> List[Dict]:
    """
    LDpred2-auto: Automatic hyperparameter inference.
    
    Parameters
    ----------
    corr : sparse matrix
        Sparse LD correlation matrix
    df_beta : pd.DataFrame
        Summary statistics with columns: beta, beta_se, n_eff
    h2_init : float
        Initial heritability estimate
    vec_p_init : float or list of float
        Initial value(s) for p (proportion causal)
        Multiple values run multiple chains for QC
    burn_in : int
        Number of burn-in iterations
    num_iter : int
        Number of iterations after burn-in
    sparse : bool
        Whether to also compute sparse solution
    shrink_corr : float
        Shrinkage for correlation matrix (0.95 recommended)
    use_mle : bool
        Whether to use MLE for alpha
    p_bounds : tuple
        Bounds for p
    alpha_bounds : tuple
        Bounds for alpha
    allow_jump_sign : bool
        Whether to allow sign changes
    seed : int, optional
        Random seed
        
    Returns
    -------
    results_list : list of dict
        List of results (one per p_init value)
        Each dict contains beta_est, postp_est, h2_est, p_est, alpha_est, etc.
    """
    check_dimensions(corr, df_beta)
    
    # Standardize (different from grid!)
    beta_hat, sd = standardize_betas_auto(df_beta)
    N = df_beta['n_eff'].values
    
    # Compute log_var (2 * log(sd))
    log_var = 2 * np.log(sd)
    
    # Compute mean LD
    mean_ld = compute_mean_ld(corr)
    
    # Convert single p_init to list
    if not isinstance(vec_p_init, list):
        vec_p_init = [vec_p_init]
    
    results_list = []
    
    for i, p_init in enumerate(vec_p_init):
        run_seed = None if seed is None else seed + i
        
        results = ldpred2_gibbs_auto(
            corr=corr,
            beta_hat=beta_hat,
            n_vec=N,
            log_var=log_var,
            h2_init=h2_init,
            p_init=p_init,
            burn_in=burn_in,
            num_iter=num_iter,
            shrink_corr=shrink_corr,
            use_mle=use_mle,
            p_bounds=p_bounds,
            alpha_bounds=(alpha_bounds[0] + 1, alpha_bounds[1] + 1),
            mean_ld=mean_ld,
            allow_jump_sign=allow_jump_sign,
            seed=run_seed
        )
        
        # Convert back to original scale
        results['beta_est'] = results['beta_est'] / sd
        
        # Compute summary statistics
        results['h2_est'] = np.nanmean(results['path_h2_est'][burn_in:])
        results['p_est'] = np.nanmean(results['path_p_est'][burn_in:])
        results['alpha_est'] = np.nanmean(results['path_alpha_est'][burn_in:])
        results['h2_init'] = h2_init
        results['p_init'] = p_init
        
        # Optionally compute sparse solution
        if sparse and not np.isnan(results['h2_est']):
            beta_sparse = ldpred2_gibbs_one(
                corr=corr,
                beta_hat=beta_hat,
                n_vec=N,
                h2=results['h2_est'],
                p=results['p_est'],
                sparse_output=True,
                burn_in=50,
                num_iter=100,
                seed=run_seed
            )
            results['beta_est_sparse'] = beta_sparse / sd
        
        results_list.append(results)
    
    return results_list

