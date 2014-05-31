# util.py

import numpy
from numpy import dot, eye

EPS = 1.e-4 # tolerance for floating-point comparison

def is_unitary(A):
    """Indicates whether or not the matrix A is unitary."""
    return numpy.max(abs(dot(A, A.conj().T) - eye(A.shape[0]))) < EPS


