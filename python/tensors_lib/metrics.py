import numpy


## Calculate the Frobenius norm of the difference matrix between two tensor fields
def diff_norm(A, B):
    Diff = A - B
    
    axes = (A.ndim - 2, A.ndim - 1)
    Norms = numpy.linalg.norm(Diff, axis=axes)
    return Norms
    

## Calculate the mean Frobenius norm of the difference matrix
#
# @A is the first matrix field to compare (usually the ground truth)
# @B is the second matrix field to compare (usually the test case)
def mean_diff_norm(A, B):
    
    Norms = diff_norm(A, B)
    
    return numpy.mean(Norms)

## Calculate the scale factor that minimizes the mean Frobenius norm of: MEAN(||A - alpha * B||)
def min_mean_diff_norm(A, B):
    
    return alpha
