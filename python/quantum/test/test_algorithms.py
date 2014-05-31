# test_algorithms.py

import unittest
from algorithms import qft
from gate_factory import hadamard_gate, phase_shift_gate
from math import pi, sqrt
import numpy
from numpy.fft import fft
from quantum_gate import tensor_power
from qubit_system import QubitSystem
import util

class Test_algorithms(unittest.TestCase):
    # Test quantum Fourier transform
    def test_qft(self):
        N = 4
        q = QubitSystem(N)
        # put q into a non-uniform mixed state
        hadamard_gate(N) * q
        tensor_power(phase_shift_gate(pi), N) * q

        # classical FFT of coefficients, normalized to be unitary like the QFT
        yhat = fft(q._QubitSystem__coeffs) / sqrt(2**N)
        qft(q)
        self.assertTrue(numpy.max(numpy.abs(q._QubitSystem__coeffs - yhat)) \
                        < util.EPS)


    # Tests Grover's search algorithm
#    def test_grover_search(self):
#        # List of strings to search
#        LIST = ['John', 'Paul', 'George', 'Ringo', 'Paul']
#
#        # Test that it can correctly find each string in the list
#        for item in LIST:
#            self.assertEqual(LIST[grover_search(item, LIST)], item)
#        # Test that it returns -1 when the item is missing
#        self.assertEqual(grover_search('Yoko', LIST), -1)
#
#
#    # Test Grover's function inversion algorithm
#    def test_grover_inverstion(self):
#        # function for testing: inverse of discrete log problem
#        f = lambda x: pow(13, x, 8)
#        n = 3 # number of bits to allow as input to f
#
#        # Test that it can solve the problem for each possible y value
#        yrange = map(f, range(2**n))
#        for y in yrange:
#            self.assertEqual(f(grover_invert(f, y, n)), y)
#        # Test that it returns -1 when the value is not in the range
#        self.assertEqual(grover_invert(f, 3, n), -1)
#
#
#    # Test Shor's algorithm
#    def test_shor(self):
#        composites = [4, 6, 9, 10, 12]
#        factors = map(shor_factor, composites)
#        # test that the factors are indeed valid factors
#        self.assertEqual(map(mod, composites, factors), [0] * len(composites))
#        # test that the factors are non-trivial
#        self.assertTrue(all(array(factors) != ones(len(composites))))
#        self.assertTrue(all(array(factors) != array(composites)))
#
#        # Test that it returns 1 for prime numbers
#        primes = [2, 3, 7, 29]
#        self.assertEqual(map(shor_factor, primes), [1] * len(primes))



if __name__ == '__main__':
    unittest.main()





