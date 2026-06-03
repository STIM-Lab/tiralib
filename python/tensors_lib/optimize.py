import numpy as np
import math
import scipy.ndimage
from scipy.optimize import minimize_scalar, differential_evolution, dual_annealing
from . import metrics
import tensors as ts


def _w_from_sigma(sigma1, sigma2):
    sigmax = max(float(sigma1), float(sigma2))
    pad = int(3 * sigmax)
    w = 2 * pad + 1
    if w < 3: w = 3
    if w % 2 == 0: w += 1
    return w

def tk_vote2_fast(field, sigma1=3, sigma2=0, power=1):
    """
    Optimized TK tensor voting (N=0 analytical plate, stick+plate).

    Equivalent to ts.tk_vote2(field, sigma1, sigma2, power, plate=True, N=0)
    but with per-offset geometry precomputed once as w×w grid arrays, removing
    all arctan2/sqrt/exp calls from the per-pixel inner work.

    Structure mirrors the reference tk_vote2 exactly:
      - Precompute Dx, Dy, c1, c2, R components, PF as w×w arrays (done once).
      - Per-pixel loop computes only polynomial operations (qTd, power, scale),
        then adds into a padded slice with direct slice += (no np.add.at).

    Parameters
    ----------
    field  : ndarray shape (H, W, 2, 2) — input symmetric tensor field
    sigma1 : float — lateral decay sigma
    sigma2 : float — axial decay sigma
    power  : int   — refinement exponent

    Returns
    -------
    ndarray shape (H, W, 2, 2) — voted tensor field
    """
    H, W = field.shape[:2]

    # ── 1. Batch eigendecomposition ──────────────────────────────────────────
    evals, evecs = np.linalg.eigh(field)      # evals: (H,W,2), evecs: (H,W,2,2)
    mags  = np.abs(evals)
    idx   = np.argsort(mags, axis=-1)         # sort by magnitude ascending

    sorted_evals = np.take_along_axis(evals, idx, axis=-1)
    sorted_mags  = np.take_along_axis(mags,  idx, axis=-1)

    # sorted eigenvectors: (H,W,2,2) — last axis is eigenvector index
    E = np.zeros_like(evecs)
    E[..., 0, :] = np.take_along_axis(evecs[..., 0, :], idx, axis=-1)
    E[..., 1, :] = np.take_along_axis(evecs[..., 1, :], idx, axis=-1)

    # ── 2. Window geometry — precomputed w×w grids ───────────────────────────
    sigmax = max(sigma1, sigma2)
    pad    = int(3 * sigmax)
    w      = 2 * pad + 1

    # Grid of offsets: x[k] = k - pad for k in 0..w-1
    # meshgrid(x, x) gives X0 varying along cols (fast/j), X1 along rows (slow/i)
    x = np.linspace(-pad, pad, w)
    X0, X1 = np.meshgrid(x, x)   # X0 = col-offset (dj), X1 = row-offset (di)

    L_sq = X0 ** 2 + X1 ** 2
    L    = np.sqrt(L_sq)

    # Decay arrays (shape w×w). Convention: exp(-ℓ²/σ²) = 1 at ℓ=0 regardless of σ.
    d1 = np.exp(-L_sq / sigma1 ** 2) if sigma1 > 0 else np.zeros_like(L)
    d2 = np.exp(-L_sq / sigma2 ** 2) if sigma2 > 0 else np.zeros_like(L)
    d1[pad, pad] = 1.0
    d2[pad, pad] = 1.0

    # Normalization for stick vote
    num     = math.pi * math.factorial(2 * power)
    den     = 2 ** (2 * power) * (math.factorial(power) ** 2)
    eta_val = den / (num * (sigma1 ** 2 + sigma2 ** 2))

    # Normalized direction (Dx=col-axis=X0/L, Dy=row-axis=X1/L)
    Dx = np.divide(X0, L, out=np.zeros_like(X0), where=L != 0)
    Dy = np.divide(X1, L, out=np.zeros_like(X1), where=L != 0)

    # Reflection matrix components — all w×w
    R11 = 1.0 - 2.0 * Dx ** 2
    R12 = -2.0 * Dx * Dy
    R22 = 1.0 - 2.0 * Dy ** 2

    # Analytical plate field (N=0) — w×w (precomputed once, identical to reference)
    c_norm   = 1.0 / (sigma1 ** 2 + sigma2 ** 2)
    ALPHA    = np.arctan2(X1, X0)   # arctan2(row_offset, col_offset)
    TWO_A    = 2.0 * ALPHA
    COS_2A   = np.cos(TWO_A)
    SIN_2A   = np.sin(TWO_A)
    M00 = 0.25 * (COS_2A + 2.0)
    M01 = 0.25 * SIN_2A
    M11 = 0.25 * (2.0 - COS_2A)
    PF = np.zeros((w, w, 2, 2))
    PF[:, :, 0, 0] = c_norm * (d1 * (1.0 - M00) + d2 * M00)
    PF[:, :, 0, 1] = c_norm * (d1 * (0.0 - M01) + d2 * M01)
    PF[:, :, 1, 0] = PF[:, :, 0, 1]
    PF[:, :, 1, 1] = c_norm * (d1 * (1.0 - M11) + d2 * M11)

    # ── 3. Allocate padded output ──────────────────────────────────────────────
    VF = np.zeros((H + 2 * pad, W + 2 * pad, 2, 2), dtype=field.dtype)

    # ── 4. Pixel loop — only polynomial inner work, direct slice += ───────────
    for i in range(H):
        for j in range(W):
            l0_mag = sorted_mags[i, j, 0]
            l1_mag = sorted_mags[i, j, 1]
            l0     = sorted_evals[i, j, 0]
            l1     = sorted_evals[i, j, 1]

            if l0_mag == 0.0 and l1_mag == 0.0:
                continue  # zero tensor — nothing to vote

            # Largest eigenvector (stick direction)
            qx = E[i, j, 0, 1]
            qy = E[i, j, 1, 1]

            vf_slice = VF[i:i + w, j:j + w]

            # ── Stick vote ─────────────────────────────────────────────────────
            scale_stick = (l1_mag - l0_mag) * np.sign(l1)
            if scale_stick != 0.0:
                # q · d — only polynomial operations over w×w, no trig
                qTd        = qx * Dx + qy * Dy
                qTd2       = qTd * qTd
                sin2_theta = 1.0 - qTd2
                cos2_theta = qTd2

                if power == 1:
                    DECAY = d1 * sin2_theta + d2 * cos2_theta
                else:
                    DECAY = d1 * np.power(sin2_theta, power) + d2 * np.power(cos2_theta, power)

                # R * q — pure linear (precomputed R11, R12, R22)
                Rq_x = R11 * qx + R12 * qy
                Rq_y = R12 * qx + R22 * qy

                factor = scale_stick * eta_val * DECAY   # w×w scalar field

                vf_slice[:, :, 0, 0] += factor * (Rq_x * Rq_x)
                vf_slice[:, :, 0, 1] += factor * (Rq_x * Rq_y)
                vf_slice[:, :, 1, 0] += factor * (Rq_x * Rq_y)
                vf_slice[:, :, 1, 1] += factor * (Rq_y * Rq_y)

            # ── Plate vote (analytical N=0) ────────────────────────────────────
            scale_plate = l0
            if scale_plate != 0.0:
                vf_slice[:, :, 0, 0] += scale_plate * PF[:, :, 0, 0]
                vf_slice[:, :, 0, 1] += scale_plate * PF[:, :, 0, 1]
                vf_slice[:, :, 1, 0] += scale_plate * PF[:, :, 0, 1]
                vf_slice[:, :, 1, 1] += scale_plate * PF[:, :, 1, 1]

    return VF[pad:pad + H, pad:pad + W, :, :]


def _grid_then_refine(objective, bounds, n_grid=30):
    """
    Coarse uniform grid scan followed by a bounded Brent refinement.

    The SIMF landscape is not guaranteed to be unimodal: the optimal scale
    factor alpha = <A,B>/<B,B> varies with the parameter, so the effective
    error surface can have shoulders or multiple local dips. A pure
    minimize_scalar (Brent's method) assumes unimodality and will stop at
    the first local minimum it encounters.

    Strategy:
      1. Evaluate the objective at n_grid uniformly spaced points.
      2. Identify the best point. Bracket it with its two neighbours.
      3. Run bounded Brent on that narrow bracket for fast, accurate refinement.
    """
    lo, hi = bounds
    grid = np.linspace(lo, hi, n_grid)
    vals = np.array([objective(x) for x in grid])
    best_idx = int(np.argmin(vals))

    # Build a tight bracket around the best grid point (clamp to [lo, hi]).
    bracket_lo = grid[max(best_idx - 1, 0)]
    bracket_hi = grid[min(best_idx + 1, n_grid - 1)]

    # If the best point is at a boundary the bracket collapses — fall back to
    # the full interval so Brent has room to move.
    if bracket_lo >= bracket_hi:
        bracket_lo, bracket_hi = lo, hi

    res = minimize_scalar(objective, bounds=(bracket_lo, bracket_hi), method='bounded')
    return res


# Minimizes SIMF(G, blur(I, σ)) over scalar σ.
#
# Previously used minimize_scalar alone (Brent's bounded method), which assumes
# a unimodal landscape and can stop at a local dip.  Replaced with a coarse
# grid scan to locate the global basin, followed by a tight Brent refinement.
def optimize_blur(G, I, bounds=(0.1, 50.0)):
    def objective(sg):
        return metrics.simf(G, scipy.ndimage.gaussian_filter(I, sigma=(sg, sg, 0, 0), mode="constant"))

    # res = minimize_scalar(objective, bounds=bounds, method='bounded')
    res = _grid_then_refine(objective, bounds, n_grid=40)
    B_opt = scipy.ndimage.gaussian_filter(I, sigma=(res.x, res.x, 0, 0), mode="constant")
    return res.x, B_opt, res.fun


# Minimizes SIMF(G, structure(gaussian_filter(image, σ))) over scalar σ.
#
# Image-domain pre-blur: unlike optimize_blur (which smooths a tensor FIELD),
# this smooths the IMAGE before the structure tensor is built, denoising the
# gradients that form it. structure_fn maps an image -> (H,W,2,2) tensor field
# (e.g. tensors.structure). sigma=0 means no blur (raw image).
def optimize_preblur(G, image, structure_fn, bounds=(0.0, 8.0)):
    def vote(sg):
        img = scipy.ndimage.gaussian_filter(image, sg) if sg > 0 else image
        return structure_fn(img)

    def objective(sg):
        B = vote(sg)
        if B.shape != G.shape:
            return np.inf
        return metrics.simf(G, B)

    res = _grid_then_refine(objective, bounds, n_grid=40)
    return res.x, vote(res.x), res.fun


# Minimizes SIMF(G, ATV(I, σ)) over scalar σ.
#
# Same fix as optimize_blur: replaced bare minimize_scalar with grid scan +
# Brent refinement to handle a potentially non-unimodal SIMF landscape.
def optimize_atv(G, I, bounds=(0.5, 20.0), cuda=False):
    # cuda=True votes on the GPU (tensors_lib.cuda_wrappers.atv_vote2_cuda),
    # which is ~1000x faster than the CPU ts.atv_vote2 the grid scan calls
    # ~30 times — the difference between a sub-second and a ~20min optimize.
    if cuda:
        from tensors_lib.cuda_wrappers import atv_vote2_cuda
        vote = lambda sg: atv_vote2_cuda(I, sigma=sg)
    else:
        vote = lambda sg: ts.atv_vote2(I, sigma=sg)

    def objective(sigma):
        B = vote(sigma)
        if B.shape != G.shape:
            return np.inf
        return metrics.simf(G, B)

    # res = minimize_scalar(objective, bounds=bounds, method='bounded')
    res = _grid_then_refine(objective, bounds, n_grid=30)
    return res.x, vote(res.x), res.fun


# Minimizes SIMF(G, TK(I, σ1, σ2, power, plate=True)) over (σ1, σ2, power)
# using differential evolution.
def optimize_tk(G, I, bounds=None, seed=42, cuda=False):
    if bounds is None:
        bounds = [(1.0, 30.0), (0.0, 20.0), (1, 30)]

    if cuda:
        from tensors_lib.cuda_wrappers import (
            _field_to_flat4, _build_nb_numpy, tk_vote2_cuda_prebuilt_flat,
        )
        # Pre-convert I once; (sigma1, sigma2, power) vary per evaluation so the
        # neighbor table must be rebuilt each time — but using vectorized numpy
        # (_build_nb_numpy) instead of the C++ scalar loop inside the extension.
        I_flat    = _field_to_flat4(I)
        out_dtype = I.dtype

        def vote_fn_cuda(s1, s2, pw):
            nb = _build_nb_numpy(s1, s2, pw, 20, _w_from_sigma(s1, s2))
            return tk_vote2_cuda_prebuilt_flat(
                I_flat, nb, s1, s2, pw, plate=True, N=20, out_dtype=out_dtype)

        def objective(params):
            if not np.all(np.isfinite(params)):
                return np.inf
            s1, s2, pw = float(params[0]), float(params[1]), max(1, int(round(params[2])))
            B = vote_fn_cuda(s1, s2, pw)
            if B.shape != G.shape:
                return np.inf
            val = metrics.simf(G, B)
            return val if np.isfinite(val) else np.inf

        res = differential_evolution(objective, bounds=bounds, seed=seed,
                                     maxiter=500, popsize=15, tol=1e-5, polish=True)
        s1, s2, pw = res.x[0], res.x[1], max(1, int(round(res.x[2])))
        return s1, s2, pw, vote_fn_cuda(s1, s2, pw), res.fun
    else:
        def objective(params):
            if not np.all(np.isfinite(params)):
                return np.inf
            s1, s2, pw = float(params[0]), float(params[1]), max(1, int(round(params[2])))
            B = ts.tk_vote2(I, sigma1=s1, sigma2=s2, power=pw, plate=True, N=20)
            if B.shape != G.shape:
                return np.inf
            val = metrics.simf(G, B)
            return val if np.isfinite(val) else np.inf

        res = differential_evolution(objective, bounds=bounds, seed=seed,
                                     maxiter=500, popsize=30, tol=1e-5, polish=True)
        s1, s2, pw = res.x[0], res.x[1], max(1, int(round(res.x[2])))
        return s1, s2, pw, ts.tk_vote2(I, sigma1=s1, sigma2=s2, power=pw, plate=True, N=20), res.fun


# Minimizes the same TK objective using dual annealing — a hybrid of fast simulated annealing
# and classical annealing.
def optimize_tk_annealing(G, I, bounds=None, seed=42):
    if bounds is None:
        bounds = [(1.0, 30.0), (0.0, 20.0), (1, 30), (0, 1)]

    def objective(params):
        s1, s2, pw = float(params[0]), float(params[1]), max(1, int(round(params[2])))
        plate = bool(round(params[3]))
        return metrics.simf(G, ts.tk_vote2(I, sigma1=s1, sigma2=s2, power=pw, plate=plate, N=20))

    res = dual_annealing(objective, bounds=bounds, seed=seed, maxiter=150)
    s1, s2, pw, plate = res.x[0], res.x[1], max(1, int(round(res.x[2]))), bool(round(res.x[3]))
    return s1, s2, pw, plate, ts.tk_vote2(I, sigma1=s1, sigma2=s2, power=pw, plate=plate, N=20), res.fun
