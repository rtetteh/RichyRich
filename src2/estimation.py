import numpy as np
import math
from scipy.optimize import minimize
from scipy.special import erf
from typing import Tuple, Optional, List

class ParameterEstimator:
    
    @staticmethod
    def estimate_alpha_improved(above_tau_data: np.ndarray, tau: float, method: str = 'combined') -> float:
        """
        Improved alpha estimation for Pareto distribution.
        """
        if len(above_tau_data) == 0:
            return np.nan
        
        # Remove extreme outliers
        Q99 = np.percentile(above_tau_data, 99)
        filtered_data = above_tau_data[above_tau_data <= Q99]
        
        if len(filtered_data) < 5:
            filtered_data = above_tau_data
            
        if method == 'mle_robust':
            return ParameterEstimator._alpha_mle_robust(filtered_data, tau)
        elif method == 'hill':
            return ParameterEstimator._alpha_hill_estimator(filtered_data, tau)
        elif method == 'method_moments':
            return ParameterEstimator._alpha_method_of_moments(filtered_data, tau)
        elif method == 'combined':
            return ParameterEstimator._alpha_combined(filtered_data, tau)
        else:
            return len(filtered_data) / np.sum(np.log(filtered_data / tau))

    @staticmethod
    def _alpha_mle_robust(data: np.ndarray, tau: float) -> float:
        log_ratios = np.log(data / tau)
        log_ratios = log_ratios[log_ratios < 20]
        
        if len(log_ratios) == 0:
            return np.nan
            
        median_log = np.median(log_ratios)
        mad = np.median(np.abs(log_ratios - median_log))
        threshold = 3 * mad
        robust_log_ratios = log_ratios[np.abs(log_ratios - median_log) <= threshold]
        
        if len(robust_log_ratios) < 3:
            robust_log_ratios = log_ratios
            
        return len(robust_log_ratios) / np.sum(robust_log_ratios)

    @staticmethod
    def _alpha_hill_estimator(data: np.ndarray, tau: float) -> float:
        log_ratios = np.log(data / tau)
        n = len(log_ratios)
        alpha_hill = 1.0 / np.mean(log_ratios)
        bias_correction = 1 - 1/(2*n) + 1/(12*n**2)
        return alpha_hill * bias_correction

    @staticmethod
    def _alpha_method_of_moments(data: np.ndarray, tau: float) -> float:
        ratios = data / tau
        sample_mean = np.mean(ratios)
        
        if sample_mean <= 1:
            return len(data) / np.sum(np.log(ratios))
            
        alpha_mom = sample_mean / (sample_mean - 1)
        n = len(data)
        bias_factor = (n - 1) / (n - 2) if n > 2 else 1
        return alpha_mom * bias_factor

    @staticmethod
    def _alpha_combined(data: np.ndarray, tau: float) -> float:
        estimates = []
        
        alpha_mle = ParameterEstimator._alpha_mle_robust(data, tau)
        if not np.isnan(alpha_mle) and alpha_mle > 0:
            estimates.append(alpha_mle)
            
        alpha_hill = ParameterEstimator._alpha_hill_estimator(data, tau)
        if not np.isnan(alpha_hill) and alpha_hill > 0:
            estimates.append(alpha_hill)
            
        alpha_mom = ParameterEstimator._alpha_method_of_moments(data, tau)
        if not np.isnan(alpha_mom) and 1 < alpha_mom < 20:
            estimates.append(alpha_mom)
            
        if not estimates:
            return np.nan
            
        if len(estimates) >= 2:
            weights = [0.5, 0.3, 0.2][:len(estimates)]
            return np.average(estimates, weights=weights)
        else:
            return estimates[0]

    @staticmethod
    def estimate_parameters(tau: float, data: np.ndarray, alpha_method: str = 'combined') -> Tuple[float, float, float]:
        """
        Estimate mu, sigma^2, and alpha.
        
        Returns:
            Tuple[float, float, float]: (mu, sigma^2, alpha)
        """
        below = data[data < tau]
        above = data[data >= tau]
        
        if len(below) == 0 or len(above) == 0:
            return np.nan, np.nan, np.nan
            
        log_below = np.log(below)

        def neg_loglik(theta):
            mu, sigma2 = theta
            if sigma2 <= 0:
                return np.inf
            
            sigma = math.sqrt(sigma2)
            ll = -0.5 * np.sum((log_below - mu) ** 2 / sigma2) - len(log_below) * np.log(sigma) - 0.5 * len(log_below) * math.log(2 * math.pi)
            
            z = (math.log(tau) - mu) / sigma
            trunc_prob = 0.5 * (1 + erf(z / math.sqrt(2)))
            if trunc_prob <= 0:
                return np.inf
            
            ll = ll - len(log_below) * math.log(trunc_prob)
            return -ll

        mu0 = np.mean(log_below)
        sigma20 = np.var(log_below, ddof=1)
        bounds = [(-np.inf, np.inf), (1e-12, None)]
        
        try:
            res = minimize(neg_loglik, x0=[mu0, sigma20], method='L-BFGS-B', bounds=bounds)
            if res.success:
                mu_hat, sigma2_hat = res.x
            else:
                mu_hat, sigma2_hat = mu0, sigma20
        except Exception:
            mu_hat, sigma2_hat = mu0, sigma20

        alpha_hat = ParameterEstimator.estimate_alpha_improved(above, tau, method=alpha_method)
        
        return mu_hat, sigma2_hat, alpha_hat
