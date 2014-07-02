# algorithms.py

"""Contains various quantum computing algorithms."""

from cmath import exp
from math import asin, ceil, log, pi, sqrt
from numpy import eye, outer

from quantum.gate_factory import function_gate, grover_diffusion_operator,\
                                 hadamard_gate
from quantum.quantum_gate import kron, QuantumGate
from quantum.qubit_system import QubitSystem
from quantum._util import gcd, int_root, is_prime

def qft(q):
    """Performs the quantum Fourier transform on a system of qubits."""
    # Make quantum gate F, the DFT matrix of dimension N := 2^n, with a unitary
    # normalization constant of 1/sqrt(N)
    N = 2**q.n()
    indices = range(N) # row and column indices
    F = QuantumGate(exp(-2.j*pi/N)**outer(indices, indices) / sqrt(N))
    # apply F to q
    F * q


def grover_search(match_text, lst):
    """Grover's quantum algorithm for searching.

    Args:
      match_text: Text to find in the list lst.
      lst: List of strings to search to find a string matching match_text.

    Returns:
      The index i of the item such that lst[i] is the same string as
        match_text. The lines must match exactly; it is not enough for the text
        to be contained in the line. If two or more lines match, it will only
        return one of the line numbers. Returns -1 if no matching line is found,
        i.e. the algorithm fails to find a solution.
    """
    if len(lst) <= 0:
        raise ValueError('List must be of positive length')

    n = len(lst)
    N = int(ceil(log(n, 2))) # number of qubits needed
    Hn = hadamard_gate(N)
    Ui = _search_oracle(match_text, lst)
    Ud = grover_diffusion_operator(N)

    MAX_ITER = 50
    count = 0
    index = n
    # Repeat until a solution is found or the iteration limit is reached
    while count < MAX_ITER and (index >= n or lst[index] != match_text):
        q = QubitSystem(N) # system of log2(n) bits in state |0>
        # apply Hadamard gate to create uniform superposition of basis states
        Hn * q

        for _ in range(_r(2**N)):
            Ui * q # apply operator that flips the sign of the matching index
            Ud * q # apply Grover's diffusion operator
        index = q.measure()
        count += 1
    return index if index < n and lst[index] == match_text else -1


def grover_invert(f, y, n):
    """Grover's algorithm for inverting a general function f that maps a
    sequence of n bits (represented as an int whose binary representation is the
    bit sequence) to another sequence of bits.

    Args:
      f: Function to invert
      y: Value of the function at which to evaluate the inverse.

    Returns:
      The input x such that f(x) = y. If more than one input suffices, it
        returns one at random. If no input suffices, returns -1.
    """
    if n <= 0:
        raise ValueError('n must be positive')

    Hn = hadamard_gate(n)
    Ui = _function_oracle(f, y, n)
    Ud = grover_diffusion_operator(n)

    MAX_ITER = 50
    count = 0
    x = None
    # Repeat until a solution is found or the iteration limit is reached
    while count < MAX_ITER and (x is None or f(x) != y):
        q = QubitSystem(n) # system of n bits in state |0>
        # apply Hadamard gate to create uniform superposition of basis states
        Hn * q

        for _ in range(_r(2**n)):
            Ui * q # apply operator that flips the sign of the matching index
            Ud * q # apply Grover's diffusion operator
        x = q.measure()
        count += 1
    return x if f(x) == y else -1


def shor_factor(n):
    """Shor's factorization algorithm.

    Args:
      n: Integer >=2 to factor.

    Returns:
      If n is composite, a non-trivial factor of n. If n is prime, returns 1.

    Raises:
      ValueError if n is <= 1.
    """
    if n <= 1:
        raise ValueError('n must be at least 2')
    if is_prime(n):
        return 1
    if n % 2 == 0:
        return 2 # even numbers > 2 are trivial

    # Need to check that n is not a power of an integer for algorithm to work
    root = int_root(n)
    if root != -1:
        return root

    # choose m s.t. n^2 <= 2^m < 2*n^2
    # log2(n^2) <= m <= log2(2 * n^2) = 1 + log2(n^2)
    m = ceil(log(n**2, 2))

    ny = ceil(log(n - 1, 2)) # number of qubits in output of function f
    I = QuantumGate(eye(2**ny))
    H = kron(hadamard_gate(m), I)

    MAX_ITER = 10
    niter = 0
    while True:
        a = n - 1 # arbitrary integer coprime to n

        # Initialize a system of qubits long enough to represent the integers 0
        # to 2^m - 1 alongside an integer up to n - 1, then apply Hadamard gate
        # to the first m qubits to create a uniform superposition.
        q = QubitSystem(m + ny)
        H * q

        # Apply the function f(x) = a^x (mod n) to the system
        f = lambda x: pow(a, x, n)
        Uf = function_gate(f, m, ny)
        Uf * q

        # Find the period of f via quantum Fourier transform
        qft(q)
        r = 1. / q.measure() # period = 1 / frequency

        niter += 1
        if niter >= MAX_ITER or (r % 2 == 0 and (a**(r / 2)) % n != -1):
            break
    return gcd(a**(r / 2) + 1, n)



# Creates the "quantum oracle/black box" gate used by the Grover search
# algorithm.

# Args:
#   match_text: Text to find in the list lst.
#   lst: List of strings to be searched to find a string matching match_text.

# Returns:
#   A gate that maps a state |k> to  -|k> if the kth item of "list" matches
#     the text match_text, and maps it to itself otherwise.
def _search_oracle(match_text, lst):
    n = len(lst)
    N = int(ceil(log(n, 2))) # number of qubits needed
    gate = QuantumGate(eye(2**N)) # identity gate
    for i in range(n):
        if lst[i] == match_text:
            gate[i][i] = -1.
    return gate


# Creates the "quantum oracle/black box" gate used by the Grover function
# inversion algorithm.

# Args:
#   f: Function to invert.
#   y: Value of the function at which to evaluate the inverse.

# Returns:
#   A gate that maps a state |x> to  -|x> if f(x) = y, and maps it to itself
#     otherwise.
def _function_oracle(f, y, n):
    gate = QuantumGate(eye(2**n)) # identity gate
    for i in range(2**n):
        if f(i) == y:
            gate[i][i] = -1.
    return gate



# Function returning the optimal number of iterations for Grover's
# algorithm. It is important to stop at exactly this many iterations or else
# future iterations may actually lower the probability of measuring the
# correct answer. See http://www.quantiki.org/wiki/Grover's_search_algorithm.
def _r(n):
    theta = asin(1. / sqrt(n))
    return int(round((pi / theta - 2.) / 4.))










