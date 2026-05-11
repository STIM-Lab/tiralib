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
    