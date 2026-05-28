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
