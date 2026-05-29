import numpy


def dropout(Field, fraction, seed=None):
    """Randomly delete tensors from the field given the fraction to remove.

    Parameters
    ----------
    Field : ndarray
        Original tensor field. A copy is returned; the input is not modified.
    fraction : float
        Fraction of tensors in ``[0, 1]`` that will be randomly zeroed out.
    seed : int, optional
        Seed for the random number generator. Passing the same seed (with the
        same field shape and fraction) reproduces the same dropout pattern.

    Returns
    -------
    Field : ndarray
        Copy of the input with a random subset of tensors set to zero.
    """
    Field = numpy.array(Field, copy=True)
    rng = numpy.random.default_rng(seed)
    mask = rng.random(Field.shape[:-2]) < fraction
    Field[mask] = 0

    return Field


def add_isotropic_background(field, density=0.02, magnitude=0.5, seed=None):
    """Sprinkle small isotropic tensors at random empty pixels.

    False-positive background that TK can suppress (low anisotropy -> low
    stick-vote contribution) but Gaussian blur and ATV cannot distinguish from
    real structure.

    Parameters
    ----------
    field : ndarray
        Input tensor field of shape ``(..., 2, 2)``. A copy is returned; the
        input is not modified.
    density : float, optional
        Probability that any given empty pixel is filled with an isotropic
        tensor. Default is 0.02.
    magnitude : float, optional
        Magnitude of the inserted isotropic tensor ``magnitude * I``.
        Default is 0.5.
    seed : int or None, optional
        Seed for the random number generator. Default is ``None``.

    Returns
    -------
    field : ndarray
        Copy of the input with isotropic tensors inserted at a random subset
        of previously-empty pixels.
    """
    field = numpy.array(field, copy=True)
    rng = numpy.random.default_rng(seed)

    norms = numpy.linalg.norm(field, axis=(-2, -1))
    empty = norms < 1e-8

    mask = rng.random(field.shape[:-2]) < density
    mask &= empty

    I2 = numpy.eye(2, dtype=field.dtype) * magnitude
    field[mask] = I2

    return field


def add_dual_noise(T, s_theta, s_eig1=None, s_eig2=None, seed=None):
    """"
    Add Gaussian noise to the orientation and eigenvalues of a tensor field T, ensuring the result remains PSD.
    Reconstruct the PSD tensor from the noisy eigenvectors and eigenvalues.
    """
    if s_eig1 == None:
        s_eig1 = s_theta
    if s_eig2 == None:
        s_eig2 = s_eig1

    rng = numpy.random.default_rng(seed)
    H, W = T.shape[:2]
    evals, evecs = numpy.linalg.eigh(T)  # evals: (H,W,2), evecs columns are eigenvectors

    # add independent Gaussian noise to each eigenvalue, clip to keep PSD
    l0 = numpy.maximum(0.0, evals[:, :, 0] + rng.normal(0.0, s_eig1, (H, W)))
    l1 = numpy.maximum(0.0, evals[:, :, 1] + rng.normal(0.0, s_eig2, (H, W)))

    # rotate eigenvectors by a noisy angle
    theta = numpy.arctan2(evecs[:, :, 1, 0], evecs[:, :, 0, 0])
    theta += rng.normal(0.0, s_theta, (H, W))

    cos_t = numpy.cos(theta)
    sin_t = numpy.sin(theta)

    # reconstruct symmetric PSD tensor: T = l0*v0*v0^T + l1*v1*v1^T -> where v0=[cos,sin], v1=[-sin,cos]
    T_noisy = numpy.zeros_like(T)
    T_noisy[:, :, 0, 0] = l0 * cos_t**2 + l1 * sin_t**2
    T_noisy[:, :, 1, 1] = l0 * sin_t**2 + l1 * cos_t**2
    T_noisy[:, :, 0, 1] = (l0 - l1) * cos_t * sin_t
    T_noisy[:, :, 1, 0] = T_noisy[:, :, 0, 1]

    return T_noisy


def add_noise(T, sigma, max_retries=5, seed=None):
    """Add Gaussian noise to a tensor field, retrying until the result is PSD.

    For each pixel, independent noise is sampled and added to all three unique
    components (T00, T01, T11). If the noisy tensor at any pixel is not PSD
    (i.e. has a negative eigenvalue), that pixel is resampled up to max_retries
    times before falling back to the original noiseless tensor at that pixel.
    """
    rng = numpy.random.default_rng(seed)
    H, W = T.shape[:2]
    T_noisy = T.copy()

    for i in range(H):
        for j in range(W):
            t = T[i, j]
            for _ in range(max_retries):
                noise = rng.normal(0.0, sigma, (2, 2))
                noise = (noise + noise.T) / 2  # keep symmetric
                candidate = t + noise
                if numpy.all(numpy.linalg.eigvalsh(candidate) >= 0):
                    T_noisy[i, j] = candidate
                    break
            # If no valid sample found, keep the original tensor (already PSD)

    return T_noisy
