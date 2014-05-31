# test_gate_factory.py

import sys
sys.path.append('../..')

import unittest
from numpy import array

import quantum
from quantum.gate_factory import controlled_gate, fredkin_gate, function_gate,\
                         hadamard_gate, toffoli_gate, x_gate
from quantum.quantum_gate import QuantumGate
from quantum.qubit_system import QubitSystem
from quantum._util import EPS

class Test_gate_factory(unittest.TestCase):
    # Tests the factory method for producing a Pauli X-gate (aka "not" gate).
    # This # is sufficient for testing all the "basic" factory methods since
    # their # internal workings are identical except for the values of the
    # particular matrix.
    def test_x_gate(self):
        X = x_gate()
        self.assertTrue((X.matrix() == array([[0., 1.], [1., 0.]])).all())

        # Test that applying the X-gate to a qubit in state |0> converts it to
        # state |1>, and vice versa
        q = QubitSystem()
        X * q
        self.assertEqual(q.measure(), 1)
        X * q
        self.assertEqual(q.measure(), 0)


    # Test Hadamard gate
    def test_hadamard_gate(self):
        N = 3
        # Matrix with 1's for all positive entries of three-qubit Hadamard gate
        is_pos = array([[1, 1, 1, 1, 1, 1, 1, 1],
                      [1, 0, 1, 0, 1, 0, 1, 0],
                      [1, 1, 0, 0, 1, 1, 0, 0],
                      [1, 0, 0, 1, 1, 0, 0, 1],
                      [1, 1, 1, 1, 0, 0, 0, 0],
                      [1, 0, 1, 0, 0, 1, 0, 1],
                      [1, 1, 0, 0, 0, 0, 1, 1],
                      [1, 0, 0, 1, 0, 1, 1, 0]])
        H = hadamard_gate(N)
        ENTRY_VAL = 2.**(-N / 2.) # abs val of each entry
        H_exp = ENTRY_VAL * (2 * is_pos - 1) # expected H
        self.assertTrue((abs(H.matrix() - H_exp) < EPS).all())


    # Tests functions for producing controlled gates
    def test_controlled_gate(self):
        # gate that permutes the coefficients on a two-qubit system
        P = QuantumGate(array([[0., 0., 0., 1.],
                                [1., 0., 0., 0.],
                                [0., 1., 0., 0.],
                                [0., 0., 1., 0.]]))
        q1 = QubitSystem(3, 0b011)
        q2 = QubitSystem(3, 0b101)
        q3 = QubitSystem(3, 0b101)

        # controlled gate permutes the second two bits iff the first bit is one
        CG = controlled_gate(P)
        CG * q1
        self.assertEqual(q1.measure(), 0b011)
        CG * q2
        self.assertEqual(q2.measure(), 0b110)

        # test Toffoli gate, a.k.a. controlled-controlled-not gate
        T = toffoli_gate()
        T * q1
        self.assertEqual(q1.measure(), 0b011)
        T * q2
        self.assertEqual(q2.measure(), 0b111)

        # Test Fredkin gate, a.k.a. controlled swap gate
        F = fredkin_gate()
        F * q1
        self.assertEqual(q1.measure(), 0b011)
        F * q3
        self.assertEqual(q3.measure(), 0b110)


    # Tests the factory method for producing a gate that implements a classical
    # function.
    def test_function_gate(self):
        f = lambda x: x >> 2
        m, k = 4, 2
        Uf = function_gate(f, m, k)

        # Initial state - M-bit binary string for 6 plus random K-bit register
        q = QubitSystem(m + k, 0b011010)
        # Expected output state: first M bits unchanged, followed by initial
        # string added to f(x) mod 2 (where x = first M bits)
        Uf * q
        self.assertEqual(q.measure(), 0b011011)



if __name__ == '__main__':
    unittest.main()



