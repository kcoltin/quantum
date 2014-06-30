# test_quantum_gate.py

import sys
sys.path.append('../..')

import unittest
from numpy import array, eye

import quantum
from quantum.gate_factory import cnot_gate, swap_gate, x_gate, y_gate, z_gate
from quantum import quantum_gate
from quantum.quantum_gate import dot, QuantumGate
from quantum.qubit_system import QubitSystem

class Test_quantum_gate(unittest.TestCase):
    # Tests the action of multiplying a quantum gate by a qubit system.
    def test_act(self):
        q1 = QubitSystem(2)
        q2 = QubitSystem(2, 0b10)
        SWAP = swap_gate()
        CNOT = cnot_gate()

        # SWAP * |00> = |00>
        SWAP * q1
        self.assertEqual(q1.measure(), 0b00)

        # SWAP * |10> = |01>, and vice versa
        SWAP * q2
        self.assertEqual(q2.measure(), 0b01)
        SWAP * q2
        self.assertEqual(q2.measure(), 0b10)

        # CNOT * |00> = |00>
        CNOT * q1
        self.assertEqual(q1.measure(), 0b00)

        # CNOT * |10> = |11> and vice versa
        CNOT * q2
        self.assertEqual(q2.measure(), 0b11)
        CNOT * q2
        self.assertEqual(q2.measure(), 0b10)

        # Test applying a gate to a subset of the qubits in the system
        q = QubitSystem(4) # state 0000
        X = x_gate()

        X.act(q, 3)
        self.assertEqual(q.measure(), 0b0010)
        X.act(q, 4)
        self.assertEqual(q.measure(), 0b0011)
        SWAP.act(q, 2)
        self.assertEqual(q.measure(), 0b0101)
        SWAP * q
        self.assertEqual(q.measure(), 0b1001)


    # Tests basic linear algebra operations on gates. Note: most linear algebra
    # operations are implemented with simple one-liners so it's not really
    # necessary to test them all; this just tests a few representative ones.
    def test_linalg(self):
        A = QuantumGate(array([[0., 1.],
                               [0., 0.]]))
        B = QuantumGate(array([[0., 0.],
                               [1., 0.]]))
        self.assertTrue((A.H().matrix() == B.matrix()).all())
        X = A + B
        self.assertTrue((X.matrix() == x_gate().matrix()).all())
        self.assertTrue((dot(X, X).matrix() == eye(2)).all())


    # Tests the tensor product operation
    def test_tensor_product(self):
        Y = y_gate()
        Z = z_gate()
        # tensor product of the Y and Z gates
        YZ = QuantumGate(array([[0., 0., -1.j, 0.],
                                 [0., 0., 0., 1.j],
                                 [1.j, 0., 0., 0.],
                                 [0., -1.j, 0., 0.]]))
        self.assertTrue((quantum_gate.kron(Y, Z) == YZ).all())


if __name__ == '__main__':
    unittest.main()

