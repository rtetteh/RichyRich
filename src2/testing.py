import numpy as np
from typing import List, Tuple, Optional
from .utils import nearest_positive_definite
from .statistics import AuxiliaryStatistics
from .distribution import HybridLogNormalPareto, HybridParams
from .estimation import ParameterEstimator

class HypothesisTest:
    
    @staticmethod
    def calculate_test_statistic_1(observed_psi: np.ndarray, bar_psi: np.ndarray) -> float:
        """Test 1 (unweighted): Euclidean distance."""
        return np.linalg.norm(np.array(observed_psi) - np.array(bar_psi)) ** 2

    @staticmethod
    def calculate_test_statistic_1opt(observed_psi: np.ndarray, bar_psi: np.ndarray, optimal_weight: np.ndarray) -> float:
        """Test 1opt (optimally weighted)."""
        observed_psi = np.array(observed_psi).reshape(-1)
        bar_psi = np.array(bar_psi).reshape(-1)
        diff = observed_psi - bar_psi
        weighted_diff = np.dot(optimal_weight, diff)
        return np.linalg.norm(weighted_diff) ** 2

    @staticmethod
    def compute_kernel(matrix: np.ndarray, sample_size: int) -> np.ndarray:
        """Compute kernel matrix for auxiliary statistics."""
        if isinstance(matrix, list):
            matrix = np.array(matrix)
        
        if matrix.ndim == 1:
            matrix = matrix.reshape(-1, 1)
        elif matrix.ndim == 0:
            matrix = np.array([[matrix]])
        
        moment_order = matrix.shape[1]
        kernel = np.zeros((moment_order, moment_order))
        means = np.mean(matrix, axis=0)
        
        # Vectorized computation
        diffs = matrix - means
        kernel = (diffs.T @ diffs) / sample_size
        return kernel

    @staticmethod
    def select_optimal_beta(
        N: int, 
        c_values: np.ndarray, 
        moment_order: int, 
        S: int, 
        R: int, 
        true_data: np.ndarray, 
        best_tau: float,
        tau_grid: np.ndarray,
        seed: int
    ) -> float:
        """
        Select optimal regularization parameter beta using cross-validation.
        
        This method splits the data into training and testing sets to evaluate
        which beta value minimizes the test statistic on unseen data.
        
        Args:
            N (int): Sample size.
            c_values (np.ndarray): Array of candidate c values (beta = c / N^(1/3)).
            moment_order (int): The moment order being analyzed.
            S (int): Number of simulated auxiliary statistics.
            R (int): Number of bootstrap replications.
            true_data (np.ndarray): The original dataset.
            best_tau (float): The optimal threshold determined for the full dataset.
            tau_grid (np.ndarray): Grid of percentiles for tau evaluation.
            seed (int): Random seed for reproducibility.
            
        Returns:
            float: The optimal beta value.
        """
        
        # Save random state
        saved_random_state = np.random.get_state()
        
        try:
            np.random.seed(seed + 9999)
            
            # Estimate initial parameters
            mu1, sigma2_1, alpha_hat = ParameterEstimator.estimate_parameters(best_tau, true_data)
            sigma1 = np.sqrt(sigma2_1)
            
            params = HybridParams(mu=mu1, sigma=sigma1, alpha=alpha_hat, tau=best_tau)
            
            # Split data based on best_tau
            n_below = int(np.sum(true_data < best_tau))
            n_above = int(np.sum(true_data >= best_tau))
            
            data_below_tau = HybridLogNormalPareto.generate_sample(n_below, params)
            data_above_tau = HybridLogNormalPareto.generate_sample(n_above, params)
            
            # Create training/testing splits
            train1_size = int(2 * len(data_below_tau) / 3)
            train2_size = int(2 * len(data_above_tau) / 3)
            
            train1, test1 = data_below_tau[:train1_size], data_below_tau[train1_size:]
            train2, test2 = data_above_tau[:train2_size], data_above_tau[train2_size:]
            
            train_data = np.concatenate([train1, train2])
            test_data = np.concatenate([test1, test2])
            train_size = len(train_data)
            test_size = len(test_data)
            
            # Evaluate tau on training data
            tau_values_local = np.percentile(train_data, tau_grid)
            best_tau2 = HypothesisTest.evaluate_tau(train_data, tau_values_local, moment_order, S=200)
            
            if best_tau2 is None:
                # Fallback if no suitable tau found on training set
                best_tau2 = best_tau

            best_beta = None
            min_norm_test = float('inf')
            
            # Pre-compute testing stats
            psi_test = AuxiliaryStatistics.compute_vector(test_data, moment_order)
            
            # OPTIMIZATION: Pre-compute simulated data ONCE for all c_values
            # The parameters (mu_train, etc.) depend on best_tau2, which is fixed for this loop.
            mu_train, sigma2_train, alpha_train = ParameterEstimator.estimate_parameters(best_tau2, train_data)
            sigma_train = np.sqrt(sigma2_train)
            params_train = HybridParams(mu=mu_train, sigma=sigma_train, alpha=alpha_train, tau=best_tau2)
            
            # Generate S simulated datasets for training (ONCE)
            psi_s_train = np.array([
                AuxiliaryStatistics.compute_vector(
                    HybridLogNormalPareto.generate_sample(train_size, params_train), 
                    moment_order
                ) for _ in range(S)
            ])
            psi_bar_s_train = np.mean(psi_s_train, axis=0)
            
            # Generate R bootstrap samples for kernel (ONCE)
            psi_r_train = np.array([
                AuxiliaryStatistics.compute_vector(
                    HybridLogNormalPareto.generate_sample(train_size, params), # Uses original params
                    moment_order
                ) for _ in range(R)
            ])
            
            # Compute kernel (ONCE)
            kernel_train = HypothesisTest.compute_kernel(psi_r_train - psi_bar_s_train, train_size)
            
            # Pre-compute K_Train (ONCE)
            integral_operator_train = (1 + 1/S) * kernel_train
            K_Train = (1 + 1/S) * integral_operator_train
            squared_integral_operator_train = np.dot(K_Train, K_Train)
            
            # Generate S simulated datasets for testing (ONCE)
            psi_s_test = np.array([
                AuxiliaryStatistics.compute_vector(
                    HybridLogNormalPareto.generate_sample(test_size, params_train), 
                    moment_order
                ) for _ in range(S)
            ])
            psi_bar_s_test = np.mean(psi_s_test, axis=0)
            z_test = psi_test - psi_bar_s_test
            
            # Generate R bootstrap samples for testing kernel (ONCE)
            psi_r_test = np.array([
                AuxiliaryStatistics.compute_vector(
                    HybridLogNormalPareto.generate_sample(test_size, params), 
                    moment_order
                ) for _ in range(R)
            ])
            
            # Compute testing kernel (ONCE)
            kernel_test = HypothesisTest.compute_kernel(psi_r_test - psi_bar_s_test, test_size)
            integral_operator_test = (1 + 1/S) * kernel_test
            K_T = (1 + 1/S) * integral_operator_test
            squared_integral_operator_test = np.dot(K_T, K_T)

            # Loop over c_values (Now very fast)
            for c in c_values:
                beta = c / (N ** (1/3))
                
                # Training Phase: Compute optimal weights
                pinv_matrix_train = np.linalg.pinv(squared_integral_operator_train + beta * np.eye(moment_order))
                optimal_weight_train = np.dot(pinv_matrix_train, K_Train)
                optimal_weight_train = optimal_weight_train / np.abs(optimal_weight_train).max()
                
                eigenvalues_train = np.linalg.eigvals(optimal_weight_train)
                if not np.all(eigenvalues_train > 0):
                    optimal_weight_train = nearest_positive_definite(optimal_weight_train)
                
                # Testing Phase: Compute optimal weights for test set
                pinv_matrix_test = np.linalg.pinv(squared_integral_operator_test + beta * np.eye(moment_order))
                optimal_weight_test = np.dot(pinv_matrix_test, K_T)
                optimal_weight_test = optimal_weight_test / np.abs(optimal_weight_test).max()
                
                eigenvalues_test = np.linalg.eigvals(optimal_weight_test)
                if not np.all(eigenvalues_test > 0):
                    optimal_weight_test = nearest_positive_definite(optimal_weight_test)
                
                norm_test = np.linalg.norm(optimal_weight_test @ z_test)
                
                if norm_test < min_norm_test:
                    min_norm_test = norm_test
                    best_beta = beta
            
            return best_beta if best_beta is not None else c_values[0] / (N ** (1/3))

        finally:
            np.random.set_state(saved_random_state)

    @staticmethod
    def evaluate_tau(data: np.ndarray, tau_values: np.ndarray, current_order: int, S: int = 200) -> Optional[float]:
        """Evaluate optimal tau."""
        best_tau = None
        min_Q = np.inf
        
        for tau in tau_values:
            try:
                mu, sigma2, alpha_hat = ParameterEstimator.estimate_parameters(tau, data)
                
                if np.isnan(mu) or np.isnan(sigma2) or np.isnan(alpha_hat):
                    continue
                
                params = HybridParams(mu=mu, sigma=np.sqrt(sigma2), alpha=alpha_hat, tau=tau)
                
                psi_obs = AuxiliaryStatistics.compute(data, current_order)
                
                psi_sim = []
                for _ in range(S):
                    sim_data = HybridLogNormalPareto.generate_sample(len(data), params)
                    psi_sim.append(AuxiliaryStatistics.compute(sim_data, current_order))
                
                psi_mean = np.mean(psi_sim)
                Q = (psi_obs - psi_mean) ** 2
                
                if Q < min_Q:
                    min_Q = Q
                    best_tau = tau
            except Exception:
                continue
                
        return best_tau
