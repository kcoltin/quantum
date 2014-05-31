# test_qubit_system.py

import sys
sys.path.append('../..')

import unittest

import quantum
from quantum.gate_factory import hadamard_gate
from quantum.qubit_system import QubitSystem

class Test_qubit_system(unittest.TestCase):
    # Test measurement of a qubit system
    def test_measure(self):
        # First, test measurement of a system in a pure state
        STATE = '01010'
        N = len(STATE)
        q = QubitSystem(N, int('0b' + STATE, 2))
        self.assertEqual(q.measure(), int('0b' + STATE, 2))
        self.assertEqual(q.smeasure(), STATE)
        for i in range(N):
            self.assertEqual(q.measure(i + 1), int(STATE[i]))

        # Test probabilistic measurement: repeatedly measure superpositions of
        # basis states and test that the outcome varies randomly.
        N = 2
        NTRIALS = 100
        # number of 1's in the first and second bit, respectively
        nsuccess1 = 0
        nsuccess2 = 0
        for i in range(NTRIALS):
            q = QubitSystem(N) # set the qubits to |00>
            hadamard_gate(N) * q # equal superposition of four basis states
            nsuccess1 += q.measure(1)
            nsuccess2 += q.measure(2)
        # Test that repeated measurement gives the same result
        state = q.measure()
        self.assertEqual(q.measure(), state)
        self.assertEqual(q.measure(), state)
        # Test that the number of 1's that appeared in the first and second
        # qubits is approximately half the total (within approx. 99% confidence
        # interval)
        CUTOFF = NTRIALS * 0.15
        self.assertTrue(abs(nsuccess1 - NTRIALS / 2.) < CUTOFF)
        self.assertTrue(abs(nsuccess2 - NTRIALS / 2.) < CUTOFF)


    # Test computation of the tensor product of two qubit systems.
    def test_tensor_product(self):
        # Test that |1> (x) |0> = |10>
        q0 = QubitSystem(1, 0)
        q1 = QubitSystem(1, 1)
        q1 * q0
        self.assertEqual(q1.measure(), 0b10)

        # Test that |0100> (x) |101> = |0100101>
        q4 = QubitSystem(4, 0b0100)
        q5 = QubitSystem(3, 0b101)
        q4 * q5
        self.assertEqual(q4.measure(), 0b0100101)




if __name__ == '__main__':
    unittest.main()



