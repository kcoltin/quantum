# quantum_gate.py

"""Contains the QuantumGate class and several closely related functions, for
representing and working with quantum logic gates."""

from math import log
import numpy
from numpy import array, copy, eye, linalg

from quantum.qubit_system import QubitSystem
from quantum._util import EPS, is_unitary

class QuantumGate(object):
    """Represents a quantum logic gate.

    The class is a thin wrapper around a Numpy array that represents the action
    of the gate.

    Attributes:
        __matrix: The unitary matrix representing the action of the gate on a
        system of qubits. A quantum gate that operates on n qubits has a matrix
        of dimension 2^n x 2^n.
    """

    def __init__(self, matrix):
        """Constructor. Creates a quantum gate with the given matrix.

        Args:
            matrix: A 2^n x 2^n array-like object, where n is the number of
            quantum bits the gate should operate on.

        Raises:
            ValueError if "matrix" has invalid dimensions for an quantum gate
            matrix.
        """
        # test that matrix is square and with dimensions that are a power of 2
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] \
           or abs(log(matrix.shape[0], 2) % 1) > EPS:
            raise ValueError('invalid matrix dimensions')
        self.__matrix = array(matrix)


    def __repr__(self):
        return 'QuantumGate(' + repr(self.__matrix) + ')'


    def __str__(self):
        return 'QuantumGate(' + str(self.__matrix) + ')'


    def copy(self):
        """Deep copy"""
        return QuantumGate(self.__matrix)


    def n(self):
        """Number of qubits that the quantum gate acts on."""
        return int(round(log(self.__matrix.shape[0], 2)))


    # Functions for accessing and setting the matrix
    def matrix(self):
        """Returns the matrix that defines the action of the quantum gate."""
        return self.__matrix

    def __getitem__(self, key):
        return self.__matrix[key]

    def __setitem__(self, key, value):
        self.__matrix[key] = value


    def act(self, qubits, index=1):
        """Performs the action of a quantum gate on a qubit system.

        Operates in-place on the system "qubits", so the original system is
        changed by interaction with the gate. This avoids violations of the no-
        cloning theorem.

        Args:
            qubits: A system of qubits for the gate to act on.
            index: Starting index of the first qubit in the system for the gate
            to act one. E.g. if the gate acts on two qubits and index = 3, then
            the gate would act on the 3rd and 4th qubits in the system (where
            qubits are indexed starting at one).

        Raises:
            ValueError if there is a dimension mismatch between the gate and the
            qubits to be operated on.
            RuntimeError if the matrix is not unitary.
        """
        if index <= 0 or index + self.n() - 1 > qubits.n():
            raise ValueError('Dimension mismatch with gate and qubit system')
        if not is_unitary(self.__matrix):
            raise RuntimeError('Non-unitary matrix')

        # construct a matrix to operate only on the desired qubits
        G = copy(self.__matrix)
        if index > 1:
            G = numpy.kron(eye(2**(index - 1)), G)
        if index + self.n() - 1 < qubits.n():
            G = numpy.kron(G, eye(2**(qubits.n() - (index + self.n() -1))))

        qubits._QubitSystem__coeffs = numpy.dot(G, qubits._QubitSystem__coeffs)


    def H(self):
        """Returns a new quantum gate whose matrix is the conjugate transpose of
        the current gate's matrix."""
        return QuantumGate(self.__matrix.conj().T)

    def abs(self):
        """Returns a new quantum gate whose matrix is the (entrywise) absolute
        value of the current gate's matrix."""
        return QuantumGate(abs(self.__matrix))

    def __pos__(self):
        """Returns a copy of the object unchanged."""
        return copy(self)

    def __neg__(self):
        """Negates the gate's matrix entrywise."""
        return QuantumGate(-self.__matrix)

    def __add__(self, other):
        """Returns a quantum gate created by adding the matrices of two gates,
        via ordinary matrix addition."""
        return QuantumGate(self.__matrix + other.__matrix)

    def __sub__(self, other):
        """Returns a quantum gate created by subtracting the matrices of two
        gates, via ordinary matrix subtraction."""
        return QuantumGate(self.__matrix - other.__matrix)

    # Entrywise matrix multiplication - also shorthand for self.act(other) where
    # other is a qubit system.
    def __mul__(self, other):
        """Returns a quantum gate created by multiplying the matrices of two
        gates, via entrywise multiplication."""
        if isinstance(other, QubitSystem):
            return self.act(other)
        else:
            return QuantumGate(self.__matrix - other.__matrix)

    def __pow__(self, e):
        """Returns a quantum gate created by raising the matrix of the gate to a
        power, using entrywise exponentiation."""
        return QuantumGate(self.__matrix**e)

    # Comparison operators
    def __lt__(self, other):
        return self.__matrix < other.__matrix
    def __le__(self, other):
        return self.__matrix <= other.__matrix
    def __eq__(self, other):
        return self.__matrix == other.__matrix
    def __ne__(self, other):
        return self.__matrix != other.__matrix
    def __gt__(self, other):
        return self.__matrix > other.__matrix
    def __ge__(self, other):
        return self.__matrix >= other.__matrix



def dot(A, B):
    """Returns a quantum gate created by multiplying the matrices of two gates,
    via matrix multiplication."""
    return QuantumGate(numpy.dot(A.matrix(), B.matrix()))

def matrix_power(G, e):
    """Returns a quantum gate created by raising the matrix of the gate to an
    integer power, using matrix exponentiation."""
    return QuantumGate(linalg.matrix_power(G.matrix(), e))

def kron(A, B):
    """Returns a quantum gate created by taking the Kronecker product (tensor
    product) of two other gates."""
    return QuantumGate(numpy.kron(A.matrix(), B.matrix()))

def tensor_power(G, e):
    """Returns the gate that results from taking the tensor product of the given
    gate with itself e times, analogous to raising the gate to an exponent e."""
    if e < 1 or e % 1 != 0:
        raise ValueError('Exponent must be a nonnegative integer')

    M = copy(G.matrix())
    for _ in range(1, e):
        M = numpy.kron(M, G.matrix())
    return QuantumGate(M)





