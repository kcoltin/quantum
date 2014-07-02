# gate_factory.py

"""Contains factory functions for creating quantum gates."""

from cmath import exp
from math import sqrt
from numpy import array, eye, ones, zeros

from quantum.quantum_gate import QuantumGate, tensor_power

def x_gate():
    """Factory method for Pauli-X gate (i.e. "not gate")"""
    return QuantumGate(array([[0., 1.],
                               [1., 0.]]))


def y_gate():
    """Factory method for Pauli-Y gate"""
    return QuantumGate(array([[0., -1.j],
                               [1.j, 0.]]))


def z_gate():
    """Factory method for Pauli-Z gate"""
    return QuantumGate(array([[1., 0.],
                               [0., -1.]]))


def hadamard_gate(n=1):
    """Factory method for Hadamard gate.

    Args:
        n: Number of qubits the gate operates on. Default is 1.
    """
    return tensor_power(QuantumGate(array([[1./sqrt(2), 1./sqrt(2)],
                                            [1./sqrt(2), -1./sqrt(2)]])), n)


def phase_shift_gate(theta):
    """Factory method for phase-shift gate.

    Args:
        theta: Phase
    """
    return QuantumGate(array([[1., 0.],
                              [0., exp(theta * 1.j)]]))


def swap_gate():
    """Factory method for swap gate"""
    return QuantumGate(array([[1., 0., 0., 0.],
                               [0., 0., 1., 0.],
                               [0., 1., 0., 0.],
                               [0., 0., 0., 1.]]))


def cnot_gate():
    """Factory method for controlled-not gate"""
    return QuantumGate(array([[1., 0., 0., 0.],
                               [0., 1., 0., 0.],
                               [0., 0., 0., 1.],
                               [0., 0., 1., 0.]]))


def controlled_gate(U):
    """Factory method for controlled-U gate.

    Args:
        U: A quantum gate operating on n qubits.
    Returns:
        The controlled-U gate operates on n+1 qubits and is of the form
        |I 0|
        |0 U|
        where U is any gate operating on n-bits and I is the 2^n x 2^n identity
        matrix.
    """
    G = eye(2**(U.n() + 1))
    G[2**U.n():,2**U.n():] = U.matrix()
    return QuantumGate(G)


def toffoli_gate():
    """Factory method for Toffoli gate a.k.a. controlled-controlled-not gate."""
    return controlled_gate(cnot_gate())


def fredkin_gate():
    """Factory method for Fredkin gate, a.k.a. controlled-swap gate."""
    return controlled_gate(swap_gate())


def function_gate(f, m, k):
    """Factory method for a gate that implements a classical function.

    Args:
        f: A function that takes a sequence of bits of length m and returns
        another sequence of bits of length k. The bit sequences are represented
        by the integers they represent in binary: e.g. "00101" would be the
        integer 5.
        m: Length of bits for the input to f.
        k: Length of bits for the output to f.
    Returns:
        The resulting gate, Uf, takes as its input a qubit system of length m+k.
        Call the first m qubits in the system "x" and the last k qubits "y":
        then the action of Uf on a pure state |x, y> is to map it to
        |x, f(x) (+) y> where (+) represents bitwise addition mod 2 (equivalent
        to the XOR operation).  It is easily shown that the map
        |x, y> --> |x, f(x) (+) y> is a bijection, which implies that Uf is a
        permutation matrix (and is therefore unitary).
    """
    # the gate operates on an (m+k)-qubit system
    G = QuantumGate(zeros((2**(m+k), 2**(m+k))))
    # filter that zeros out all but the last k bits
    FILTER = (1 << k) - 1

    for j in range(G.matrix().shape[1]):
        y = j & FILTER # y is the last k bits of the jth state
        x = j >> k # x is the state ignoring the last k bits
        # The gate maps |x,y> to |x,y + f(x)>, where + is bitwise addition mod 2
        i = (x << k) + (y ^ f(x))
        G[i,j] = 1.

    return G


def grover_diffusion_operator(n):
    """Factory method for Grover diffusion operator.

    Every non-diagonal entry is 1/2^(n-1); diagonal entries are that value minus
    one.
    Note: at least one source I've see defines the operator as the negative of
    this matrix; either definition works for Grover's algorithm.

    Args:
        n: Number of qubits for the gate to operate it.
    """
    return QuantumGate(1. / 2.**(n - 1) * ones((2**n, 2**n)) - eye(2**n))



