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


## Mask-restricted SIMF: same scale-invariant Frobenius comparison but the mean
# is taken only over pixels where the mask is True. Use this when the phantom
# is mostly empty background — global SIMF dilutes the structure-recovery
# signal with low-error background pixels and drives optimizers toward "don't
# filter" solutions.
def masked_simf(A, B, mask=None):
    if mask is None:
        mask = _support_mask(A)
    if not np.any(mask):
        return 0.0
    alpha = _scale_factor(A, B)
    norms = diff_norm(A, B, alpha)
    return float(np.mean(norms[mask]))


## Scaled Invarient Frobenius Norm - SIF
def sif(A, B):
    alpha = _scale_factor(A, B)
    return diff_norm(A, B, alpha)


def _dominant_eigvec(T):
    """Return the eigenvector of the largest-magnitude eigenvalue at each pixel.
    Shape: (..., 2). Sign-disambiguated so v[..., 0] >= 0 (sticks are sign-free).
    """
    evals, evecs = np.linalg.eigh(T)
    mags = np.abs(evals)
    idx = np.argmax(mags, axis=-1)                # (...,)
    take = idx[..., None, None]                   # (...,1,1)
    # eigh returns eigenvectors as columns; pick the column at `take`.
    v = np.take_along_axis(evecs, take.repeat(2, axis=-2), axis=-1)[..., 0]  # (...,2)
    # Disambiguate sign so undirected-line comparisons are well defined.
    flip = v[..., 0] < 0
    v[flip] = -v[flip]
    return v


def _support_mask(A, eps=1e-8):
    return np.linalg.norm(A, axis=(-2, -1)) > eps


## Mean angular error (radians) between dominant eigenvectors of A and B,
# evaluated only at pixels where A has nonzero structure. Sticks are undirected
# so angles are folded into [0, pi/2].
def angular_error(A, B, mask=None):
    if mask is None:
        mask = _support_mask(A)
    if not np.any(mask):
        return 0.0
    vA = _dominant_eigvec(A)
    vB = _dominant_eigvec(B)
    dot = np.abs(np.sum(vA * vB, axis=-1))        # |cos angle| — undirected
    dot = np.clip(dot, 0.0, 1.0)
    ang = np.arccos(dot)                          # [0, pi/2]
    return float(np.mean(ang[mask]))


def _eccentricity_field(T):
    evals = np.linalg.eigh(T)[0]
    a = np.abs(evals[..., 1])
    b = np.abs(evals[..., 0])
    safe = np.where(a > 1e-12, a, 1.0)
    ratio = np.clip(b / safe, 0.0, 1.0)
    return np.sqrt(1.0 - ratio * ratio)           # 1 = stick, 0 = isotropic


## Mean absolute eccentricity difference on GT support — captures whether the
# method preserves stick-vs-plate character (e.g. at junctions, where blur
# destroys anisotropy by averaging incompatible orientations).
def anisotropy_error(A, B, mask=None):
    if mask is None:
        mask = _support_mask(A)
    if not np.any(mask):
        return 0.0
    eA = _eccentricity_field(A)
    eB = _eccentricity_field(B)
    return float(np.mean(np.abs(eA[mask] - eB[mask])))


## Mean alignment-weighted magnitude of B at pixels where the GT line had a
# carved-out gap. Higher = more successful gap bridging.
#   score = mean( ||B||_F * |<vB, gt_dir>| ) over gap pixels
# where gt_dir is the unit stick direction the line should have at the gap.
def gap_bridging_score(B, gap_mask, gt_direction):
    if not np.any(gap_mask):
        return 0.0
    gt = np.asarray(gt_direction, dtype=float)
    gt = gt / (np.linalg.norm(gt) + 1e-12)
    vB = _dominant_eigvec(B)
    mag = np.linalg.norm(B, axis=(-2, -1))
    align = np.abs(vB @ gt)                       # |cos| with GT direction
    return float(np.mean((mag * align)[gap_mask]))


## Calculate the scale factor that minimizes the mean Frobenius norm of: MEAN(||A - alpha * B||)
def _scale_factor(A, B):
    num = np.sum(A * B)          # sum of element-wise products
    den = np.sum(B * B)          # sum of squared elements of B (sum of squared Frobenius norms)

    # avoid division by zero
    if den < 1e-8:
        return 0.0
    
    alpha = num / den
    return alpha
