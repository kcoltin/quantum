# algorithms.py
#
# Contains various algorithms for quantum computing.

from cmath import exp
from math import pi, sqrt
from numpy import outer

from quantum_gate import QuantumGate

def qft(q):
    """Performs the quantum Fourier transform on a system of qubits."""
    # Make quantum gate F, the DFT matrix of dimension N := 2^n, with a unitary
    # normalization constant of 1/sqrt(N)
    N = 2**q.n()
    indices = range(N) # row and column indices
    F = QuantumGate(exp(-2.j*pi/N)**outer(indices, indices) / sqrt(N))
    # apply F to q
    F * q



