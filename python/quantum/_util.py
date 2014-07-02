# util.py

from math import floor, log, sqrt
import numpy
from numpy import dot, eye

EPS = 1.e-4 # tolerance for floating-point comparison

def is_unitary(A):
    """Indicates whether or not the matrix A is unitary."""
    return numpy.max(abs(dot(A, A.conj().T) - eye(A.shape[0]))) < EPS


def gcd(a, b):
    """Returns the greatest common divisor of a and b, using the Euclidean
    algorithm."""
    if a <= 0 or b <= 0:
        raise ValueError('Arguments must be positive integers')

    while b != 0:
        tmp = b
        b = a % b
        a = tmp
    return a


def is_prime(n):
    """Tests whether n is prime.

    This uses a naive algorithm that would be useless for numbers of substantial
    size, but is fine for the purposes of QAS, which are understanding and
    simulating the use of quantum computing rather than actually solving
    difficult problem instances.
    """
    if  n <= 1:
        raise ValueError('n must be greater than 1')

    # trial  division
    for k in range(2, int(floor(sqrt(n))) + 1):
        if n % k == 0:
            return False
    return True



def int_root(n):
    """Checks whether n is the k-th power of an integer for some k > 1. If so,
    it returns the root r such that r^k = n. Otherwise, returns -1."""
    if n <= 0:
        raise ValueError('n must be positive')
    k = 2
    while k <= log(n, 2):
        root = n**(1. / k)
        iroot = round(root)
        if iroot**k == n:
            return iroot
        k += 1
    return -1


