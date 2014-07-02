# qubit_system.py

"""Contains the class QubitSystem for representing a system of one or more
qubits."""

from math import log, sqrt
from numpy import arange, array, cumsum, kron, nonzero, random, zeros

class QubitSystem(object):
    """Represents a system of one or more qubits.

    Attributes:
        __coeffs: The coefficients on each of the basis states. __coeffs[i] is
        the coefficient on the basis state that is the binary representation of
        the integer i. For example, a 2-qubit system
        a|00> + b|01> + c|10> + d|11> would have coeffs = [a, b, c, d].
    """

    def __init__(self, n=1, state=0):
        """Constructor for a system of zero or more qubits initialized to a pure
        state.

        Args:
            n: Number of qubits in the system. Must be a nonnegative integer.
            Default 1.
            state: Pure state that the system is initialized to. The system is
            initialized to the pure state that is the binary representation of
            the integer "state".

        Raises:
            ValueError if illogical values are passed for n or state.
        """
        if n < 0 or n % 1 != 0:
            raise ValueError('n must be a nonnegative integer')
        if state < 0 or state >= 2**n or state % 1 != 0:
            raise ValueError('state must be an integer between 0 and 2^n - 1')

        self.__coeffs = zeros(2**n)
        self.__coeffs[state] = 1.


    def __repr__(self):
        return '<qubit_system.QubitSystem object at ' + hex(id(self)) + '>'


    def __str__(self):
        return 'QubitSystem of ' + str(self.n()) + ' quantum bits'


    def copy(self):
        """You can't make a copy of a QubitSystem: this enforces the no-cloning
        theorem."""
        raise TypeError('Cannot copy qubits: violation of no-cloning theorem')


    def n(self):
        """Returns the number of qubits in the system."""
        return int(round(log(len(self.__coeffs), 2)))


    def measure(self, bit_index=None):
        """Measures the state of the qubit system.

        Args:
            bit_index: If supplied, measures the state of the qubit given by
            bit_index. (One-based indexing is used for qubits, so 1 is the first
            qubit and so on.) If omitted, measures the state of all qubits in
            the system.

        Returns:
            If bit_index is supplied, returns either 0 or 1 depending on the
            measured value of the bit_index-th qubit. If bit_index is absent,
            returns an integer whose binary representation is the values of all
            the qubits. E.g., if the system is measured to be in the state 0110,
            returns 5.

        Raises:
            ValueError if bit_index is out of range.
        """
        if bit_index is not None and (bit_index < 1 or bit_index > self.n()):
            raise ValueError('bit_index out of range')

        # Randomly observe a state, based on the probability amplitudes
        rand = random.random()
        cumprobs = cumsum(abs(self.__coeffs)**2)
        state = int(min(nonzero(cumprobs >= rand)[0])) # Note 1

        if bit_index is None:
            # collapse to the observed state
            self.__coeffs = zeros(len(self.__coeffs))
            self.__coeffs[state] = 1.
            return state
        else:
            # value of bit_index-th qubit
            bit = (state >> (self.n() - bit_index)) & 1
            # collapse only the bit_index-th qubit to the observed value
            self.__collapse_qubit(bit_index, bit)
            return bit


    def smeasure(self):
        """Measures the state of the qubit system as a binary string.

        Returns:
            A string of 0's and 1's corresponding to the observed values of each
            qubit. E.g., if the system is observed in the state 011, returns
            "011".
        """
        state = bin(self.measure())[2:] # state as binary int sans leading zeros
        return state.zfill(self.n())


    def __collapse_qubit(self, bit_index, observed_value):
        """Collapses the system of qubits to those states where the bit given by
        bit_index has the given observed value.

        E.g., if q is a two-qubit system, then calling q.__collapse_qubit(2,1)
        would collapse it to a superposition of the states |01> and |11>.

        Args:
            bit_index: Index of qubit to collapse into a particular state.
            observed_value: The value of the bit_index-th bit that was observed,
            either 0 or 1.
        """
        # Interval over which the values 0 and 1 alternate for the bit given by
        # bit_index. For example, if bit_index = 2 and n = 3, then interval
        # would equal 2, because the bit is zero for 2 states, then one for 2
        # states, then zero for 2, and then one for 2:
        # |000> |001> |010> |011> |100> |101> |110> |111>
        interval = 2**(self.n() - bit_index)

        # Set the probabilities of all states where the bit given by bit_index
        # is NOT equal to observed_value to zero.
        # Iterate over the first state in each "block" of consecutive states
        # where the bit given by bit_index is NOT equal to observed_value. Note
        # that the operation "observed_value ^ 1" maps 0 to 1 and vice versa.
        for i in range((observed_value ^ 1) * interval, self.n(), 2 * interval):
            self.__coeffs[i+arange(interval)] = 0.

        # Then, normalize the remaining states such that their squared
        # coefficients sum to one.
        self.__coeffs /= sqrt(sum(abs(self.__coeffs)**2))


    def __mul__(self, other):
        """Combines two systems of qubits.

        On return, the coefficients on the current system of qubits are replaced
        by the tensor product of the two sets of coefficients. The system
        "other" is "erased" in that it is resized to a system of size zero. This
        represents the fact that the two systems are combined into a single
        system of separable states, and avoids violating the no-cloning theorem.

        Args:
            other: Another system of qubits to join with the current system.
        """
        self.__coeffs = kron(self.__coeffs, other.__coeffs)
        other.__coeffs = array([])



# NOTES
# 1. The int() cast is to change the output of nonzero() from int64, which can
#    cause weird errors - noticed when calling pow(a, b, c) where b is an int64
#    and a, c are ints.




