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
    delta: float = 0.05  # Contamination fraction (δ) - fraction of sample from misspecified Pareto
    taumis: float = 0.05  # Tau deviation parameter (τ_d = τ × (1 + taumis))

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
    Deviated hybrid distribution: LogNormal body + Pareto tail with tau deviation contamination.
    Uses delta as contamination fraction (δ) where δ*n observations are randomly replaced
    with Pareto draws using tau_d = tau * (1 + taumis) instead of the true tau.
    """
    
    @staticmethod
    def generate_sample(n: int, params: DeviatedParams, seed: Optional[int] = None) -> np.ndarray:
        """
        Generate deviated hybrid sample with tau deviation contamination.
        
        A fraction δ (params.delta) of the sample is randomly selected and replaced
        with Pareto draws using tau_d = tau * (1 + taumis) instead of the true tau.
        
        Args:
            n (int): Sample size.
            params (DeviatedParams): Distribution parameters.
            seed (int, optional): Random seed for reproducibility.
            
        Returns:
            np.ndarray: Generated data with contamination.
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

        # Step 2: Generate Pareto tail >= tau
        v2 = stats.genpareto.rvs(
            c=1.0 / params.alpha,
            loc=params.tau,
            scale=params.tau / params.alpha,
            size=n2,
            random_state=rng
        )

        # Step 3: Apply tau deviation contamination to Pareto tail subsample
        # Generate deviations from Pareto subsample - replace last d observations
        d = int(np.floor(params.delta * n))  # number of misspecified observations
        if d > 0:
            tau_d = params.tau * (1.0 + params.taumis)  # misspecified tau parameter
            
            # Generate misspecified Pareto draws (size 2*d to ensure enough > tau)
            devv2 = stats.genpareto.rvs(
                c=1.0 / params.alpha,
                loc=tau_d,
                scale=tau_d / params.alpha,
                size=2 * d,
                random_state=rng
            )
            
            # Create v2dev: first (n2-d) from v2, last d from devv2 (only those > tau)
            v2dev = np.zeros(n2)
            v2dev[:n2 - d] = v2[:n2 - d]
            
            # Fill last d positions with valid devv2 values (> tau)
            c1 = n2 - d  # Start index for replacement
            i = 0  # Index in devv2
            while c1 < n2 and i < len(devv2):
                if devv2[i] > params.tau:
                    v2dev[c1] = devv2[i]
                    c1 += 1
                i += 1
            
            # If we don't have enough valid values, fill remaining with the last valid value
            while c1 < n2:
                v2dev[c1] = devv2[i-1] if i > 0 else tau_d
                c1 += 1
        else:
            v2dev = v2

        # Step 4: Combine data
        data = np.concatenate([dat1, v2dev])
        
        # Step 5: FORCE exact percentile by strategic placement
        sorted_data = np.sort(data)
        
        percentile_position = params.tau_perc * (n - 1)
        lower_index = int(np.floor(percentile_position))
        upper_index = int(np.ceil(percentile_position))
        
        if lower_index == upper_index:
            sorted_data[lower_index] = params.tau
        else:
            sorted_data[lower_index] = params.tau
            sorted_data[upper_index] = params.tau
        
        # Step 6: Shuffle
        rng.shuffle(sorted_data)
        
        return sorted_data
