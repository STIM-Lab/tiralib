"""Reconstruction: convert a voted 2D tensor field into a binary mask.

Generic (works on any (H, W, 2, 2) tensor field), so it lives in tiralib rather
than in a vessel-specific repo. The first method is native tensor-voting
saliency-thresholding; level-set and SymPR variants are stubbed for later.

See the tk-vesselvote design spec
(docs/superpowers/specs/2026-06-08-2d-reconstruction-design.md).
"""

import numpy as np


def stick_saliency(T):
    """Per-pixel stick saliency s = |lambda1| - |lambda0| of a tensor field.

    lambda1, lambda0 are the larger/smaller-magnitude eigenvalues. High saliency
    marks confident, oriented (stick-like) structure; near-zero marks isotropic
    (ball/plate) or empty regions. Sign is irrelevant (sticks are undirected).

    Args:
        T: tensor field, shape (H, W, 2, 2).

    Returns:
        s: float64 saliency map, shape (H, W), all >= 0.
    """
    evals = np.linalg.eigh(T)[0]                 # ascending by value
    mags = np.abs(evals)
    mags.sort(axis=-1)                           # ascending by magnitude
    return mags[..., 1] - mags[..., 0]


def _otsu_threshold(values):
    """Otsu's threshold on a 1-D array of nonnegative values.

    Standard between-class-variance maximization over a 256-bin histogram.
    Returns the threshold in the original value units.
    """
    vmax = float(values.max())
    if vmax <= 0:
        return 0.0
    hist, edges = np.histogram(values, bins=256, range=(0.0, vmax))
    hist = hist.astype(np.float64)
    total = hist.sum()
    if total == 0:
        return 0.0
    # Bin centers.
    centers = (edges[:-1] + edges[1:]) / 2.0
    w0 = np.cumsum(hist)                         # weight of background class
    w1 = total - w0                              # weight of foreground class
    # Cumulative means.
    cum = np.cumsum(hist * centers)
    mu_total = cum[-1]
    # Avoid division by zero where a class is empty.
    with np.errstate(divide="ignore", invalid="ignore"):
        mu0 = np.where(w0 > 0, cum / w0, 0.0)
        mu1 = np.where(w1 > 0, (mu_total - cum) / w1, 0.0)
    between = w0 * w1 * (mu0 - mu1) ** 2
    return float(centers[np.argmax(between)])


def saliency_threshold(T, fill=False, threshold="otsu", close_radius=0):
    """Reconstruct a binary mask from a voted tensor field by thresholding stick
    saliency.

    Args:
        T:            tensor field, shape (H, W, 2, 2).
        fill:         if False (default), the mask is the thresholded saliency
                      directly -- use when the field fires on the region you want
                      (interior/medial fields). If True, treat thresholded pixels
                      as boundary curves: optionally close gaps, then flood-fill
                      the enclosed regions -- use when the field fires on vessel
                      walls/edges.
        threshold:    "otsu" (default; unsupervised, fits the training-free story)
                      or an explicit float applied to the saliency map.
        close_radius: morphological-closing radius (pixels) applied before the
                      flood fill when fill=True; 0 disables it. Ignored when
                      fill=False.

    Returns:
        mask: bool array, shape (H, W).
    """
    s = stick_saliency(T)

    if threshold == "otsu":
        thr = _otsu_threshold(s.ravel())
    else:
        thr = float(threshold)

    binary = s > thr

    if not fill:
        return binary

    return _fill_enclosed(binary, close_radius)


def _fill_enclosed(boundary, close_radius):
    """Close small gaps in a boundary mask, then fill enclosed background.

    A pixel is "enclosed" if it is background that cannot reach the image border
    through background. Implemented by flood-filling the background from the
    border and inverting: anything the flood did not reach is interior.
    """
    from scipy.ndimage import binary_closing, binary_fill_holes

    b = boundary.copy()
    if close_radius and close_radius > 0:
        r = int(close_radius)
        size = 2 * r + 1
        struct = np.ones((size, size), dtype=bool)
        b = binary_closing(b, structure=struct)

    # binary_fill_holes fills regions of background fully enclosed by True.
    return binary_fill_holes(b)


# --- Future methods (stubbed; see design spec) ------------------------------

def levelset_reconstruct(T, **kwargs):
    """DEFERRED (Option B): tensor-driven level-set / active-contour evolution.

    The voted field supplies the edge-stopping / advection terms; the zero level
    set converges to the vessel boundary, giving a closed contour = mask.
    """
    raise NotImplementedError(
        "levelset_reconstruct is a future method (Option B); see the 2026-06-08 "
        "reconstruction design spec."
    )


def sympr_reconstruct(T, **kwargs):
    """DEFERRED (Option C): 2D Symmetrized Poisson Reconstruction.

    Minimize the quartic energy ||grad chi (.) grad chi - T||^2 + alpha sum chi^2
    by coarse-to-fine coordinate descent, with the voted field injected directly
    as the target tensor field T (skipping the paper's lossy outer-product splat).
    Reference C++ (mkazhdan/SymmetricPoissonReconstruction) is 3D-hardcoded, so
    this is a from-scratch 2D Python implementation.
    """
    raise NotImplementedError(
        "sympr_reconstruct is a future method (Option C); see the 2026-06-08 "
        "reconstruction design spec."
    )
