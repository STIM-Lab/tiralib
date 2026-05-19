import numpy

## Create a field of stick tensors in the specified orientation.
#
# @direction is the orientation of the stick tensor specified as a vector
# @size is the size of the tensor field to generate
def stick_field(direction, size):
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

## Randomly delete tensors from the field given the fraction to remove
#
# @Field is the original tensor field (a copy will be returned)
# @fraction is the fraction of tensors [0, 1] that will be randomly deleted
def dropout(Field, fraction):
    
    Field = numpy.array(Field, copy=True)
    mask = numpy.random.rand(*Field.shape[:-2]) < fraction
    Field[mask] = 0
    
    return Field


## create a grid of alternating horizontal/vertical stick tensors
# even cells -> T = [[1, 0], [0, 0]]
# odd cells  -> T = [[0, 0], [0, 1]]
# cell_size controls the size of the alternating blocks -> smaller, higher frequency
def tensor_grid(N, cell_size=1, noise=0.0):
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
    