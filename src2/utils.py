import numpy as np
from scipy.linalg import eigh

def nearest_positive_definite(A, epsilon=1e-6):
    """
    Ensure matrix A is strictly positive definite.
    
    Args:
        A (np.ndarray): Input matrix.
        epsilon (float): Minimum eigenvalue threshold.
        
    Returns:
        np.ndarray: Positive definite matrix.
    """
    # Make symmetric
    B = (A + A.T) / 2
    # Eigenvalue decomposition
    eigvals, eigvecs = eigh(B)
    # Shift eigenvalues to ensure positive definiteness
    eigvals[eigvals < epsilon] = epsilon
    # Reconstruct the matrix
    A_pd = eigvecs @ np.diag(eigvals) @ eigvecs.T
    return A_pd
