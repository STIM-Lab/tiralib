import numpy

def stick_field(direction, size):
    """Create a field of stick tensors in the specified orientation.

    Parameters
    ----------
    direction : array_like
        Orientation of the stick tensor, specified as a 1D vector.
    size : array_like
        Shape of the tensor field to generate, as a 1D vector.

    Returns
    -------
    Field : ndarray
        Tensor field of shape ``(*size, dim, dim)`` filled with the outer
        product of ``direction`` with itself.
    """
    direction = numpy.asarray(direction)
    if direction.ndim != 1:
        raise ValueError("direction must be a 1D vector")
        
    size = numpy.asarray(size)
    if size.ndim != 1:
        raise ValueError("size must be a 1D vector")
    
    if len(size) != len(direction):
        raise ValueError("direction and size must be the same dimension")
        
    dimension = len(direction)
    
    field_dimension = numpy.append(size, (dimension, dimension))    
    Field = numpy.zeros(field_dimension)
    
    T = numpy.outer(direction, direction)
    
    Field[..., :, :] = T
    
    return Field


def tensor_grid(N, cell_size=1, noise=0.0):
    """Create a grid of alternating horizontal/vertical stick tensors.

    Even cells get ``T = [[1, 0], [0, 0]]`` (horizontal sticks); odd cells get
    ``T = [[0, 0], [0, 1]]`` (vertical sticks).

    Parameters
    ----------
    N : int
        Grid side length; output shape is ``(N, N, 2, 2)``.
    cell_size : int, optional
        Size of the alternating blocks. Smaller values produce higher frequency
        patterns. Default is 1.
    noise : float, optional
        Standard deviation of symmetric Gaussian noise added to each tensor.
        Default is 0.0 (no noise).

    Returns
    -------
    T : ndarray
        Tensor field of shape ``(N, N, 2, 2)`` with dtype ``float32``.
    """
    T = numpy.zeros((N, N, 2, 2), dtype=numpy.float32)
    
    i, j = numpy.indices((N, N))                            # create coordinate grids
    
    mask = ((i // cell_size) + (j // cell_size)) % 2 == 0   # true: even blocks, false: odd blocks
    
    # horizontal sticks in True mask regions
    T[mask, 0, 0] = 1.0  # Ixx
    T[mask, 1, 1] = 0.0  # Iyy
    T[mask, 0, 1] = 0.0  # Ixy
    T[mask, 1, 0] = 0.0  
    
    # vertical sticks in False mask regions
    T[~mask, 0, 0] = 0.0  # Ixx
    T[~mask, 1, 1] = 1.0  # Iyy
    T[~mask, 0, 1] = 0.0  # Ixy
    T[~mask, 1, 0] = 0.0
    
    # add optional symmetric Gaussian noise
    if noise != 0:
        noise_tensor = numpy.random.normal(0, noise, T.shape)
        noise_tensor[:, :, 0, 1] = noise_tensor[:, :, 1, 0]
        T += noise_tensor

    return T


def _stick_outer(theta):
    """Compute the outer product of a stick direction vector."""
    v = numpy.array([numpy.cos(theta), numpy.sin(theta)])
    return numpy.outer(v, v)


def crossing(shape, th1=0.0, th2=numpy.pi / 2, width=1):
    """True 2D crossing: two bars with orientations that sum at the junction.

    At the intersection, ``T = T(th1) + T(th2)`` — a legitimate plate-like tensor
    that TK is designed to recover via plate voting. (Blur cannot reconstruct this
    from either single-orientation neighborhood.)

    Parameters
    ----------
    shape : tuple of int
        Output shape as ``(rows, cols)``.
    th1 : float, optional
        Stick angle in radians for the horizontal bar. Default is 0.0.
    th2 : float, optional
        Stick angle in radians for the vertical bar. Default is ``pi/2``.
    width : int, optional
        Bar thickness in pixels. Default is 1.

    Returns
    -------
    T : ndarray
        Tensor field of shape ``(rows, cols, 2, 2)`` with dtype ``float32``.
    """
    rows, cols = shape
    T = numpy.zeros((rows, cols, 2, 2), dtype=numpy.float32)
    T1 = _stick_outer(th1).astype(numpy.float32)
    T2 = _stick_outer(th2).astype(numpy.float32)

    half = width // 2
    rc, cc = rows // 2, cols // 2
    r_lo, r_hi = max(rc - half, 0), min(rc + half + (width % 2), rows)
    c_lo, c_hi = max(cc - half, 0), min(cc + half + (width % 2), cols)

    # Horizontal bar gets th1, vertical bar gets th2; intersection gets the sum.
    bar_h_mask = numpy.zeros((rows, cols), dtype=bool)
    bar_v_mask = numpy.zeros((rows, cols), dtype=bool)
    bar_h_mask[r_lo:r_hi, :] = True
    bar_v_mask[:, c_lo:c_hi] = True

    T[bar_h_mask] = T1
    T[bar_v_mask] = T2
    junction = bar_h_mask & bar_v_mask
    T[junction] = T1 + T2

    return T


def gapped_line(shape, theta=numpy.pi / 2, gap_frac=0.25, width=1):
    """Stripe of co-oriented sticks with a centered gap of zero tensors.

    Parameters
    ----------
    shape : tuple of int
        Output shape as ``(rows, cols)``.
    theta : float, optional
        Stored stick direction in radians. Default is ``pi/2``.
    gap_frac : float, optional
        Fraction of the row axis carved out as a gap in the middle.
        Default is 0.25.
    width : int, optional
        Line thickness in pixels. Default is 1.

    Returns
    -------
    T : ndarray
        Tensor field of shape ``(rows, cols, 2, 2)`` with dtype ``float32``.
    gap_mask : ndarray of bool
        ``True`` at pixels that lie on the ground-truth line but were zeroed
        out, used to compute the gap-bridging score.
    """
    rows, cols = shape
    T = numpy.zeros((rows, cols, 2, 2), dtype=numpy.float32)
    Ts = _stick_outer(theta).astype(numpy.float32)

    # The line runs through the center along the axis perpendicular to "across-width".
    # For simplicity we draw a vertical bar (constant column band) and let theta set
    # the stored stick direction; this keeps the demo geometry trivial while still
    # testing direction recovery.
    half = width // 2
    cc = cols // 2
    c_lo, c_hi = max(cc - half, 0), min(cc + half + (width % 2), cols)

    line_mask = numpy.zeros((rows, cols), dtype=bool)
    line_mask[:, c_lo:c_hi] = True
    T[line_mask] = Ts

    # Carve a gap out of the middle along the row axis.
    gap_len = int(gap_frac * rows)
    g_lo = (rows - gap_len) // 2
    g_hi = g_lo + gap_len
    gap_mask = numpy.zeros((rows, cols), dtype=bool)
    gap_mask[g_lo:g_hi, c_lo:c_hi] = True
    T[gap_mask] = 0.0

    return T, gap_mask


def arc(shape, radius=None, width=1, arc_span=(0.0, 2 * numpy.pi)):
    """Circular arc of stick tensors, each storing the local tangent direction.

    Tests recovery of smoothly varying orientation — blur averages across the
    normal direction and dulls anisotropy, TK reinforces along the tangent.

    Parameters
    ----------
    shape : tuple of int
        Output shape as ``(rows, cols)``.
    radius : float, optional
        Arc radius in pixels. Defaults to ``0.4 * min(rows, cols)``.
    width : int, optional
        Arc thickness in pixels. Default is 1.
    arc_span : tuple of float, optional
        ``(a0, a1)`` angular span in radians. If ``a0 < a1`` the arc covers
        ``[a0, a1]``; otherwise it wraps through ``±pi``. Default is the full
        circle ``(0, 2*pi)``.

    Returns
    -------
    T : ndarray
        Tensor field of shape ``(rows, cols, 2, 2)`` with dtype ``float32``.
    """
    rows, cols = shape
    if radius is None:
        radius = 0.4 * min(rows, cols)
    T = numpy.zeros((rows, cols, 2, 2), dtype=numpy.float32)

    rc, cc = rows / 2.0, cols / 2.0
    yy, xx = numpy.indices((rows, cols))
    dy = yy - rc
    dx = xx - cc
    r = numpy.sqrt(dx * dx + dy * dy)
    phi = numpy.arctan2(dy, dx)

    half = max(width / 2.0, 0.5)
    on_arc = numpy.abs(r - radius) <= half

    a0, a1 = arc_span
    if a0 < a1:
        in_span = (phi >= a0) & (phi <= a1)
    else:
        in_span = (phi >= a0) | (phi <= a1)
    mask = on_arc & in_span

    # Tangent direction at each pixel: perpendicular to radial = (-sin phi, cos phi)
    tx = -numpy.sin(phi)
    ty = numpy.cos(phi)
    T[mask, 0, 0] = (tx * tx)[mask]
    T[mask, 0, 1] = (tx * ty)[mask]
    T[mask, 1, 0] = (tx * ty)[mask]
    T[mask, 1, 1] = (ty * ty)[mask]

    return T


def parallel_lines(shape, theta=numpy.pi / 2, separation=6, width=1):
    """Two parallel stripes of co-oriented sticks.

    Tests whether the method merges nearby structures. Blur with sigma
    comparable to ``separation`` will fuse them; TK with σ₂ along the stick
    axis should keep them distinct.

    Parameters
    ----------
    shape : tuple of int
        Output shape as ``(rows, cols)``.
    theta : float, optional
        Stored stick direction in radians. Default is ``pi/2``.
    separation : int, optional
        Distance in pixels between the two stripe centers. Default is 6.
    width : int, optional
        Stripe thickness in pixels. Default is 1.

    Returns
    -------
    T : ndarray
        Tensor field of shape ``(rows, cols, 2, 2)`` with dtype ``float32``.
    """
    rows, cols = shape
    T = numpy.zeros((rows, cols, 2, 2), dtype=numpy.float32)
    Ts = _stick_outer(theta).astype(numpy.float32)

    half = width // 2
    cc = cols // 2
    offset = separation // 2

    for sign in (-1, +1):
        c_center = cc + sign * offset
        c_lo = max(c_center - half, 0)
        c_hi = min(c_center + half + (width % 2), cols)
        T[:, c_lo:c_hi] = Ts

    return T
