"""
Maximum Likelihood Estimation for alpha parameter.

This module implements MLE optimization for the alpha parameter that controls
the relationship between allele frequencies and effect sizes.
"""

import numpy as np
from scipy.optimize import minimize
from typing import Tuple, List


class MLEAlpha:
    """
    Maximum Likelihood Estimator for alpha parameter.
    
    The model assumes:
        beta_j^2 ~ sigma2 * exp(alpha * log_var_j)
    
    where log_var_j is related to the allele frequency.
    
    The negative log-likelihood is:
        L = alpha * sum(log_var) + n * log(sigma2) + sum(beta^2 * exp(-alpha * log_var)) / sigma2
    """
    
    def __init__(self, ind_causal: List[int], log_var: np.ndarray, 
                 curr_beta: np.ndarray, bootstrap: bool = False):
        """
        Initialize MLE optimizer.
        
        Parameters
        ----------
        ind_causal : list of int
            Indices of causal variants
        log_var : np.ndarray
            Log variance values for all variants
        curr_beta : np.ndarray
            Current beta estimates for all variants
        bootstrap : bool
            Whether to use bootstrap sampling
        """
        self.nb = len(ind_causal)
        
        # Extract values for causal variants only
        if bootstrap:
            # Bootstrap: sample with replacement
            indices = np.random.choice(self.nb, size=self.nb, replace=True)
            sampled_ind = [ind_causal[i] for i in indices]
        else:
            sampled_ind = ind_causal
        
        self.a = np.array([log_var[j] for j in sampled_ind])
        self.b = np.array([curr_beta[j] ** 2 for j in sampled_ind])
        self.sum_a = np.sum(self.a)
    
    def objective(self, par: np.ndarray) -> float:
        """
        Compute negative log-likelihood.
        
        Parameters
        ----------
        par : np.ndarray
            Parameters [alpha_plus_one, sigma2]
            
        Returns
        -------
        nll : float
            Negative log-likelihood
        """
        alpha_plus_one = par[0]
        sigma2 = par[1]
        
        # Compute sum_c = sum(b * exp(-alpha * a))
        sum_c = np.sum(self.b * np.exp(-alpha_plus_one * self.a))
        
        # Negative log-likelihood
        nll = alpha_plus_one * self.sum_a + self.nb * np.log(sigma2) + sum_c / sigma2
        
        return nll
    
    def gradient(self, par: np.ndarray) -> np.ndarray:
        """
        Compute gradient of negative log-likelihood.
        
        Parameters
        ----------
        par : np.ndarray
            Parameters [alpha_plus_one, sigma2]
            
        Returns
        -------
        grad : np.ndarray
            Gradient vector
        """
        alpha_plus_one = par[0]
        sigma2 = par[1]
        
        # Compute intermediate values
        c_k = self.b * np.exp(-alpha_plus_one * self.a)
        sum_c = np.sum(c_k)
        sum_ac = np.sum(self.a * c_k)
        
        # Gradient
        grad = np.array([
            self.sum_a - sum_ac / sigma2,
            (self.nb - sum_c / sigma2) / sigma2
        ])
        
        return grad


def mle_alpha(par_init: np.ndarray,
              ind_causal: List[int],
              log_var: np.ndarray,
              curr_beta: np.ndarray,
              alpha_bounds: Tuple[float, float],
              bootstrap: bool = False,
              verbose: bool = False) -> np.ndarray:
    """
    Perform MLE optimization for alpha parameter.
    
    Parameters
    ----------
    par_init : np.ndarray
        Initial parameters [alpha_plus_one, sigma2]
    ind_causal : list of int
        Indices of causal variants
    log_var : np.ndarray
        Log variance values
    curr_beta : np.ndarray
        Current beta estimates
    alpha_bounds : tuple
        Bounds for alpha parameter (min, max)
    bootstrap : bool
        Whether to use bootstrap sampling
    verbose : bool
        Whether to print optimization progress
        
    Returns
    -------
    par : np.ndarray
        Optimized parameters [alpha_plus_one, sigma2]
    """
    # Create MLE object
    mle = MLEAlpha(ind_causal, log_var, curr_beta, bootstrap)
    
    # Set bounds
    # alpha_plus_one bounds: alpha_bounds
    # sigma2 bounds: [par_init[1] / 2, par_init[1] * 2]
    bounds = [
        (alpha_bounds[0], alpha_bounds[1]),
        (par_init[1] / 2.0, par_init[1] * 2.0)
    ]
    
    # Optimize using L-BFGS-B
    result = minimize(
        fun=mle.objective,
        x0=par_init,
        method='L-BFGS-B',
        jac=mle.gradient,
        bounds=bounds,
        options={'disp': verbose}
    )
    
    if not result.success and verbose:
        print(f"Warning: Optimization did not converge: {result.message}")
    
    return result.x

