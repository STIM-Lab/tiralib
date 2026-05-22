"""
cuda_wrappers.py

Optional CUDA-accelerated tensor voting wrappers.

Drop-in replacements for ts.tk_vote2() and ts.atv_vote2() that call the
pybind11 CUDA extensions when available, and fall back to the pure-Python
implementations automatically when the extensions have not been built.

Usage
-----
    from tensors_lib.cuda_wrappers import tk_vote2_cuda, atv_vote2_cuda, CUDA_TK_AVAILABLE, CUDA_ATV_AVAILABLE

Building the extensions
-----------------------
    cd python/cuda_ext && mkdir build && cd build
    cmake .. -DCMAKE_PREFIX_PATH=/path/to/venv
    make -j$(nproc)
    # The resulting .so files land in the build/ dir; copy them to python/cuda_ext/
    # or set PYTHONPATH to include that directory.

Field layout convention
-----------------------
  - tensors.py and optimize.py use (H, W, 2, 2) float64 arrays.
  - The CUDA extension expects (H, W, 4) float32 C-contiguous arrays,
    where the four values per pixel are [T00, T01, T10, T11].
  - _field_to_flat4 and _flat4_to_field handle the conversion.
"""

import sys
import os
import numpy as np
import tensors as ts

# ---------------------------------------------------------------------------
# Locate and import CUDA extensions
# ---------------------------------------------------------------------------
_CUDA_EXT_DIR = os.path.join(os.path.dirname(__file__), '..', 'cuda_ext')
_abs_ext_dir = os.path.abspath(_CUDA_EXT_DIR)
if _abs_ext_dir not in sys.path:
    sys.path.insert(0, _abs_ext_dir)

try:
    import tensorvote_tk as _cuda_tk_mod
    CUDA_TK_AVAILABLE = True
except ImportError:
    _cuda_tk_mod = None
    CUDA_TK_AVAILABLE = False

try:
    import tensorvote_atv as _cuda_atv_mod
    CUDA_ATV_AVAILABLE = True
except ImportError:
    _cuda_atv_mod = None
    CUDA_ATV_AVAILABLE = False


# ---------------------------------------------------------------------------
# Neighbor2D dtype — must match tira::Neighbor2D exactly (48 bytes):
#   int32 du, dv  |  float32 dx, dy, c1, c2, Ra, Rb, Rc, Va, Vb, Vc
# ---------------------------------------------------------------------------

_NB_DTYPE = np.dtype([
    ('du', np.int32),   ('dv', np.int32),
    ('dx', np.float32), ('dy', np.float32),
    ('c1', np.float32), ('c2', np.float32),
    ('Ra', np.float32), ('Rb', np.float32), ('Rc', np.float32),
    ('Va', np.float32), ('Vb', np.float32), ('Vc', np.float32),
])  # itemsize == 48


_TV_EPSILON = 1e-12


def _build_nb_numpy(sigma1, sigma2, power, samples, w):
    """
    Build a tira::Neighbor2D table using vectorized numpy.

    Matches cpu::build_neighbors2d exactly but replaces the O(w^2 * N) scalar
    C++ loop with N vectorized numpy passes over all w^2 offsets — much faster
    for large sigma or power.

    Parameters
    ----------
    sigma1  : float — lateral sigma
    sigma2  : float — axial sigma (0 = no axial channel)
    power   : int   — refinement exponent
    samples : int   — plate numerical-integration samples (matches 'N' in Python)
    w       : int   — window side length (odd)

    Returns
    -------
    ndarray, shape (w*w,), dtype _NB_DTYPE — contiguous structured array
    """
    hw = w // 2

    # All (dv, du) pairs — dv = row offset, du = col offset
    # Outer loop dv, inner loop du (matches cpu::build_neighbors2d order)
    dv_g, du_g = np.meshgrid(np.arange(-hw, hw + 1, dtype=np.int32),
                              np.arange(-hw, hw + 1, dtype=np.int32),
                              indexing='ij')
    dv = dv_g.ravel()
    du = du_g.ravel()

    du_f = du.astype(np.float32)
    dv_f = dv.astype(np.float32)
    l2   = du_f ** 2 + dv_f ** 2
    L    = np.sqrt(l2)
    nz   = L > 0.0

    dx = np.where(nz, du_f / np.where(nz, L, 1.0), 0.0).astype(np.float32)
    dy = np.where(nz, dv_f / np.where(nz, L, 1.0), 0.0).astype(np.float32)

    if sigma1 > 0:
        c1 = np.exp(-(l2 / sigma1 ** 2)).astype(np.float32)
    else:
        c1 = np.where(l2 > _TV_EPSILON, 0.0, 1.0).astype(np.float32)

    if sigma2 > 0:
        c2 = np.exp(-(l2 / sigma2 ** 2)).astype(np.float32)
    else:
        c2 = np.where(l2 > _TV_EPSILON, 0.0, 1.0).astype(np.float32)

    Ra = (1.0 - 2.0 * dx * dx).astype(np.float32)
    Rb = (-2.0 * dx * dy).astype(np.float32)
    Rc = (1.0 - 2.0 * dy * dy).astype(np.float32)

    # Plate kernel: vectorized numerical integration matching tk_plate_numerical.
    # Result = (2/N) * sum_{n=0}^{N-1} tk_stickvote_2d(du, dv, theta_n)
    n_nb   = len(du)
    Va_acc = np.zeros(n_nb, dtype=np.float64)
    Vb_acc = np.zeros(n_nb, dtype=np.float64)
    Vc_acc = np.zeros(n_nb, dtype=np.float64)

    Ra_f64 = Ra.astype(np.float64)
    Rb_f64 = Rb.astype(np.float64)
    Rc_f64 = Rc.astype(np.float64)
    c1_f64 = c1.astype(np.float64)
    c2_f64 = c2.astype(np.float64)
    dx_f64 = dx.astype(np.float64)
    dy_f64 = dy.astype(np.float64)

    dtheta = np.pi / samples
    for n in range(samples):
        theta = n * dtheta
        qx    = np.cos(theta)
        qy    = np.sin(theta)

        qTd  = qx * dx_f64 + qy * dy_f64
        qTd2 = qTd * qTd
        t1   = 1.0 - qTd2
        t2   = qTd2

        if power == 1:
            eta = c1_f64 * t1 + c2_f64 * t2
        else:
            eta = c1_f64 * (t1 ** power) + c2_f64 * (t2 ** power)

        Rq_x = Ra_f64 * qx + Rb_f64 * qy
        Rq_y = Rb_f64 * qx + Rc_f64 * qy

        Va_acc += eta * Rq_x * Rq_x
        Vb_acc += eta * Rq_x * Rq_y
        Vc_acc += eta * Rq_y * Rq_y

    Va = (Va_acc * 2.0 / samples).astype(np.float32)
    Vb = (Vb_acc * 2.0 / samples).astype(np.float32)
    Vc = (Vc_acc * 2.0 / samples).astype(np.float32)

    nb = np.empty(n_nb, dtype=_NB_DTYPE)
    nb['du'] = du;  nb['dv'] = dv
    nb['dx'] = dx;  nb['dy'] = dy
    nb['c1'] = c1;  nb['c2'] = c2
    nb['Ra'] = Ra;  nb['Rb'] = Rb;  nb['Rc'] = Rc
    nb['Va'] = Va;  nb['Vb'] = Vb;  nb['Vc'] = Vc
    return np.ascontiguousarray(nb)


# ---------------------------------------------------------------------------
# Internal layout converters
# ---------------------------------------------------------------------------

def _field_to_flat4(field):
    """
    Convert (H, W, 2, 2) tensor field (any float dtype) to (H, W, 4) float32.
    Layout: [T00, T01, T10, T11] per pixel, C-contiguous.
    """
    H, W = field.shape[:2]
    flat = np.empty((H, W, 4), dtype=np.float32)
    flat[:, :, 0] = field[:, :, 0, 0]
    flat[:, :, 1] = field[:, :, 0, 1]
    flat[:, :, 2] = field[:, :, 1, 0]
    flat[:, :, 3] = field[:, :, 1, 1]
    return np.ascontiguousarray(flat)


def _flat4_to_field(flat, out_dtype=np.float64):
    """
    Convert (H, W, 4) float32 flat layout back to (H, W, 2, 2).
    """
    H, W = flat.shape[:2]
    field = np.empty((H, W, 2, 2), dtype=out_dtype)
    field[:, :, 0, 0] = flat[:, :, 0]
    field[:, :, 0, 1] = flat[:, :, 1]
    field[:, :, 1, 0] = flat[:, :, 2]
    field[:, :, 1, 1] = flat[:, :, 3]
    return field


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def tk_vote2_cuda_flat(flat_in, H, W, sigma1, sigma2=0.0, power=1, plate=True, N=20, out_dtype=np.float64):
    """
    Fast-path CUDA TK tensor voting for repeated calls with the same field.

    Accepts a pre-converted (H, W, 4) float32 C-contiguous array (as produced by
    _field_to_flat4), bypassing the dtype conversion and array copy on every call.
    Use this inside optimization loops where the input field is fixed but (sigma1,
    sigma2, power) vary across evaluations.

    Falls back to ts.tk_vote2() when the CUDA extension is not available.

    Parameters
    ----------
    flat_in  : ndarray, shape (H, W, 4), float32, C-contiguous — pre-converted field
    H, W     : int — field dimensions (must match flat_in.shape[:2])
    sigma1   : float — orthogonal (lateral) sigma
    sigma2   : float — axial (parallel) sigma; 0 = no parallel channel
    power    : int   — refinement exponent (>= 1)
    plate    : bool  — include plate votes
    N        : int   — plate numerical integration samples (0 -> uses 10 internally)
    out_dtype: dtype — dtype for the returned (H, W, 2, 2) array

    Returns
    -------
    ndarray, shape (H, W, 2, 2), dtype out_dtype
    """
    if not CUDA_TK_AVAILABLE:
        # Reconstruct a float64 field for the CPU fallback
        field = _flat4_to_field(flat_in, out_dtype=np.float64)
        return ts.tk_vote2(field, sigma1=sigma1, sigma2=sigma2, power=power, plate=plate, N=N)

    samples = max(1, int(N)) if N > 0 else 10
    if not plate:
        samples = 1

    flat_out = _cuda_tk_mod.tensorvote_tk(
        flat_in,
        float(sigma1),
        float(sigma2),
        int(power),
        0,          # w=0: auto-compute from sigma
        True,       # stick always enabled
        bool(plate),
        samples,
    )
    return _flat4_to_field(flat_out, out_dtype=out_dtype)


def tk_vote2_cuda_prebuilt_flat(flat_in, nb, sigma1, sigma2=0.0, power=1,
                                plate=True, N=20, out_dtype=np.float64):
    """
    TK voting using a neighbor table built externally by _build_nb_numpy().

    Passes the pre-built table directly to the CUDA extension, bypassing the
    internal cpu::build_neighbors2d call that dominates at large sigma/power.

    Parameters
    ----------
    flat_in  : ndarray, shape (H, W, 4), float32, C-contiguous
    H, W     : int — field dimensions
    nb       : ndarray, shape (w*w,), dtype _NB_DTYPE — from _build_nb_numpy()
    sigma1   : float
    sigma2   : float
    power    : int
    plate    : bool
    N        : int — plate samples (must match what was used to build nb)
    out_dtype: dtype

    Returns
    -------
    ndarray, shape (H, W, 2, 2), dtype out_dtype
    """
    if not CUDA_TK_AVAILABLE:
        field = _flat4_to_field(flat_in, out_dtype=np.float64)
        return ts.tk_vote2(field, sigma1=sigma1, sigma2=sigma2, power=power, plate=plate, N=N)

    samples = max(1, int(N)) if N > 0 else 10
    if not plate:
        samples = 1

    sigmax = max(float(sigma1), float(sigma2))
    pad    = int(3 * sigmax)
    w      = 2 * pad + 1
    if w < 3:  w = 3
    if w % 2 == 0: w += 1

    nb_bytes = np.frombuffer(nb.tobytes(), dtype=np.uint8)

    flat_out = _cuda_tk_mod.tensorvote_tk_prebuilt(
        flat_in,
        nb_bytes,
        int(len(nb)),
        int(w),
        float(sigma1),
        float(sigma2),
        int(power),
        True,
        bool(plate),
        samples,
    )
    return _flat4_to_field(flat_out, out_dtype=out_dtype)


def tk_vote2_cuda(field, sigma1, sigma2=0.0, power=1, plate=True, N=20):
    """
    CUDA-accelerated TK tensor voting.

    Drop-in replacement for ts.tk_vote2(field, sigma1, sigma2, power, plate, N).
    Falls back to ts.tk_vote2() when the CUDA extension is not available.

    Parameters
    ----------
    field  : ndarray, shape (H, W, 2, 2) — input tensor field
    sigma1 : float — orthogonal (lateral) sigma
    sigma2 : float — axial (parallel) sigma; 0 = no parallel channel
    power  : int   — refinement exponent (>= 1)
    plate  : bool  — include plate votes
    N      : int   — plate numerical integration samples.
                     N=0 in the Python ts.tk_vote2 triggers an analytical closed-form;
                     the CUDA extension does not implement that analytical shortcut — it
                     uses numerical integration with samples>=1.  For N=0 inputs we pass
                     samples=10 to the CUDA extension as a reasonable default.
                     Use N>=1 for faithful reproduction of the numerical path.

    Returns
    -------
    ndarray, shape (H, W, 2, 2)
    """
    if not CUDA_TK_AVAILABLE:
        return ts.tk_vote2(field, sigma1=sigma1, sigma2=sigma2, power=power, plate=plate, N=N)

    out_dtype = field.dtype
    flat_in = _field_to_flat4(field)

    # Map Python N to CUDA samples:
    #   N=0 analytical -> use 10 samples (CUDA has no closed-form shortcut)
    #   N>0 numerical  -> use N samples
    if plate:
        samples = max(1, int(N)) if N > 0 else 10
    else:
        samples = 1  # value is ignored when plate=False, but must be >= 1

    flat_out = _cuda_tk_mod.tensorvote_tk(
        flat_in,
        float(sigma1),
        float(sigma2),
        int(power),
        0,          # w=0: auto-compute as 2*floor(3*max(sigma1,sigma2))+1
        True,       # stick always enabled
        bool(plate),
        samples,
    )
    return _flat4_to_field(flat_out, out_dtype=out_dtype)


def atv_vote2_cuda(field, sigma):
    """
    CUDA-accelerated ATV tensor voting.

    Drop-in replacement for ts.atv_vote2(field, sigma).
    Falls back to ts.atv_vote2() when the CUDA extension is not available.

    Parameters
    ----------
    field : ndarray, shape (H, W, 2, 2)
    sigma : float — isotropic decay sigma

    Returns
    -------
    ndarray, shape (H, W, 2, 2)
    """
    if not CUDA_ATV_AVAILABLE:
        return ts.atv_vote2(field, sigma=sigma)

    out_dtype = field.dtype
    flat_in = _field_to_flat4(field)
    flat_out = _cuda_atv_mod.tensorvote_atv(flat_in, float(sigma))
    return _flat4_to_field(flat_out, out_dtype=out_dtype)
