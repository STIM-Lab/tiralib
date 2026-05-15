import numpy as np


def diff_norm(A, B, alpha=1.0):
    """
    Calculate the Frobenius norm of the difference matrix between two tensor fields, after scaling B by alpha.
    Args:
        A:     First tensor field (e.g., ground truth), shape (H, W, 2, 2)
        B:     Second tensor field (e.g., test case), shape (H, W, 2, 2)
        alpha: Scaling factor applied to B before comparison
    Returns:
        Frobenius norm of the difference between A and alpha * B
    """
    Diff = A - alpha * B
    
    axes = (A.ndim - 2, A.ndim - 1)
    Norms = np.linalg.norm(Diff, axis=axes)
    return Norms
    

def mean_diff_norm(A, B, alpha=1.0):
    """
    Calculate the mean Frobenius norm of the difference matrix between two tensor fields, after scaling B by alpha.
    Args:
        A:     First tensor field (e.g., ground truth), shape (H, W, 2, 2)
        B:     Second tensor field (e.g., test case), shape (H, W, 2, 2)
        alpha: Scaling factor applied to B before comparison
    Returns:
        Mean Frobenius norm of the difference between A and alpha * B
    """
    Norms = diff_norm(A, B, alpha)
    
    return np.mean(Norms)


## Scaled Invarient Mean Frobenius Norm - SIMF
def simf(A, B):
    alpha = _scale_factor(A, B)
    return mean_diff_norm(A, B, alpha)


## Scaled Invarient Frobenius Norm - SIF
def sif(A, B):
    alpha = _scale_factor(A, B)
    return diff_norm(A, B, alpha)


## Calculate the scale factor that minimizes the mean Frobenius norm of: MEAN(||A - alpha * B||)
def _scale_factor(A, B):
    num = np.sum(A * B)          # sum of element-wise products
    den = np.sum(B * B)          # sum of squared elements of B (sum of squared Frobenius norms)

    # avoid division by zero
    if den < 1e-8:
        return 0.0
    
    alpha = num / den
    return alpha
