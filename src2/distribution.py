import numpy as np
from dataclasses import dataclass
from typing import Optional
from scipy import stats

@dataclass
class HybridParams:
    mu: float
    sigma: float  # This is the standard deviation (sqrt of variance)
    alpha: float
    tau: float
    tau_perc: float = 0.82

@dataclass
class DeviatedParams:
    mu: float
    sigma: float  # This is the standard deviation (sqrt of variance)
    alpha: float
    tau: float
    tau_perc: float = 0.82
    delta: float = 0.1  # fraction of sample from log-Cauchy

class HybridLogNormalPareto:
    """
    Hybrid distribution: LogNormal body + Pareto tail.
    """
    
    @staticmethod
    def generate_sample(n: int, params: HybridParams, seed: Optional[int] = None) -> np.ndarray:
        """
        Generate hybrid sample with EXACT percentile placement.
        
        Args:
            n (int): Sample size.
            params (HybridParams): Distribution parameters.
            seed (int, optional): Random seed.
            
        Returns:
            np.ndarray: Generated data.
        """
        if seed is not None:
            rng = np.random.RandomState(seed)
        else:
            rng = np.random

        n1 = int(np.floor(params.tau_perc * n))
        n2 = n - n1

        # Step 1: Generate sufficient lognormal data < tau (strictly less than)
        dat1 = []
        max_attempts = 20
        attempts = 0
        
        while len(dat1) < n1 and attempts < max_attempts:
            draws = max(10 * n1, 10000)
            # Note: params.sigma is expected to be std dev here
            v1 = np.exp(rng.normal(loc=params.mu, scale=params.sigma, size=draws))
            keep = v1[v1 < params.tau]
            
            if keep.size > 0:
                to_take = min(keep.size, n1 - len(dat1))
                dat1.extend(keep[:to_take].tolist())
            
            attempts += 1
        
        # Ensure we have exactly n1 values
        if len(dat1) < n1:
            remaining = n1 - len(dat1)
            fill_values = params.tau - rng.uniform(1e-4, 1e-2, remaining)
            dat1.extend(fill_values.tolist())
        
        dat1 = np.array(dat1[:n1])

        # Step 2: Generate Pareto tail >= tau
        # X = tau / U^(1/alpha)
        U = rng.uniform(size=n2)
        v2 = params.tau / (U ** (1.0 / params.alpha))
        v2 = np.maximum(v2, params.tau)

        # Step 3: Combine data
        data = np.concatenate([dat1, v2])
        
        # Step 4: FORCE exact percentile by strategic placement
        sorted_data = np.sort(data)
        
        percentile_position = params.tau_perc * (n - 1)
        lower_index = int(np.floor(percentile_position))
        upper_index = int(np.ceil(percentile_position))
        
        if lower_index == upper_index:
            sorted_data[lower_index] = params.tau
        else:
            sorted_data[lower_index] = params.tau
            sorted_data[upper_index] = params.tau
        
        # Step 5: Shuffle
        rng.shuffle(sorted_data)
        
        return sorted_data


class DeviatedHybridLogNormalPareto:
    """
    Deviated hybrid distribution: LogNormal body + (1-delta)*Pareto + delta*log-Cauchy tail.
    """
    
    @staticmethod
    def generate_log_cauchy_sample(n: int, m: float, s: float, tau: float, seed: Optional[int] = None) -> np.ndarray:
        """
        Generate log-Cauchy sample > tau.
        
        Parameters m and s estimated from log of Pareto tail.
        """
        if seed is not None:
            rng = np.random.RandomState(seed)
        else:
            rng = np.random
        
        # Generate Cauchy: m + s * tan(pi * (U - 0.5))
        U = rng.rand(n)
        cauchy_vals = m + s * np.tan(np.pi * (U - 0.5))
        Y = np.exp(cauchy_vals)
        # Ensure > tau
        Y = np.maximum(Y, tau + 1e-4)
        return Y
    
    @staticmethod
    def generate_sample(n: int, params: DeviatedParams, seed: Optional[int] = None) -> np.ndarray:
        """
        Generate deviated hybrid sample.
        
        Args:
            n (int): Sample size.
            params (DeviatedParams): Distribution parameters.
            seed (int, optional): Random seed.
            
        Returns:
            np.ndarray: Generated data.
        """
        if seed is not None:
            rng = np.random.RandomState(seed)
        else:
            rng = np.random

        n1 = int(np.floor(params.tau_perc * n))
        n2 = n - n1

        # Step 1: Generate sufficient lognormal data < tau (strictly less than)
        dat1 = []
        max_attempts = 20
        attempts = 0
        
        while len(dat1) < n1 and attempts < max_attempts:
            draws = max(10 * n1, 10000)
            v1 = np.exp(rng.normal(loc=params.mu, scale=params.sigma, size=draws))
            keep = v1[v1 < params.tau]
            
            if keep.size > 0:
                to_take = min(keep.size, n1 - len(dat1))
                dat1.extend(keep[:to_take].tolist())
            
            attempts += 1
        
        # Ensure we have exactly n1 values
        if len(dat1) < n1:
            remaining = n1 - len(dat1)
            fill_values = params.tau - rng.uniform(1e-4, 1e-2, remaining)
            dat1.extend(fill_values.tolist())
        
        dat1 = np.array(dat1[:n1])

        # Step 2: Generate full Pareto tail
        k = 1.0 / params.alpha
        sigma_pareto = params.tau / params.alpha
        theta = params.tau
        v2 = stats.genpareto.rvs(c=k, scale=sigma_pareto, loc=theta, size=n2, random_state=rng)

        # Step 3: Generate deviations
        d = int(np.floor(params.delta * n))
        if d > 0:
            # Estimate log-Cauchy parameters from v2
            log_v2 = np.log(v2)
            m = np.median(log_v2)
            s = np.median(np.abs(log_v2 - m))
            
            # Generate log-Cauchy deviations
            dev_v2 = DeviatedHybridLogNormalPareto.generate_log_cauchy_sample(d, m, s, params.tau, seed)
            
            # Replace the last d observations in v2
            v2 = np.concatenate([v2[:-d], dev_v2])

        # Step 3: Combine data
        data = np.concatenate([dat1, v2])
        
        # Step 4: FORCE exact percentile by strategic placement
        sorted_data = np.sort(data)
        
        percentile_position = params.tau_perc * (n - 1)
        lower_index = int(np.floor(percentile_position))
        upper_index = int(np.ceil(percentile_position))
        
        if lower_index == upper_index:
            sorted_data[lower_index] = params.tau
        else:
            sorted_data[lower_index] = params.tau
            sorted_data[upper_index] = params.tau
        
        # Step 5: Shuffle
        rng.shuffle(sorted_data)
        
        return sorted_data
