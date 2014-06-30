// Tests for the functionality of the QuantumGate class.

#include "gate_factory.h"
#include "quantum_gate.h" 
#include "qubit_system.h" 
#include "util.h"
#include <gtest/gtest.h>
using std::complex;
using std::string;
using namespace quantum_algorithm_simulator; 

// Tolerance for floating-point comparison
static const double EPS = 1.e-6; 

// Tests the action of multiplying a quantum gate by a qubit system. 
TEST (QuantumGateTest, TestAct) {
	int N = 2; 
	QubitSystem q1(N);
	QubitSystem q2(N, 2);
	const QuantumGate SWAP = qgates::swap_gate(); 
	const QuantumGate CNOT = qgates::cnot_gate(); 

	// SWAP * |00> = |00>
	SWAP * q1; 
	EXPECT_EQ("00", q1.smeasure()); 

	// SWAP * |10> = |01>, and vice versa
	SWAP * q2; 
	EXPECT_EQ("01", q2.smeasure());
	SWAP * q2; 
	EXPECT_EQ("10", q2.smeasure());
	
	// CNOT * |00> = |00> 
	CNOT * q1; 
	EXPECT_EQ("00", q1.smeasure()); 
	
	// CNOT * |10> = |11> and vice versa
	CNOT * q2; 
	EXPECT_EQ("11", q2.smeasure());
	CNOT * q2; 
	EXPECT_EQ("10", q2.smeasure()); 

	// Test applying a gate to a subset of the qubits in the system
	N = 4; 
	QubitSystem q(N); // state 0000
	const QuantumGate X = qgates::x_gate(); 

	X.act(&q, 3); 
	EXPECT_EQ("0010", q.smeasure()); 
	X.act(&q, 4); 
	EXPECT_EQ("0011", q.smeasure()); 
	SWAP.act(&q, 2);  
	EXPECT_EQ("0101", q.smeasure()); 
	SWAP * q; 
	EXPECT_EQ("1001", q.smeasure()); 
}


// Tests the operation of taking the conjugate transpose of a quantum gate to
// obtain the inverse gate. 
TEST (QuantumGateTest, TestTranspose) { 
	const QuantumGate gate= qgates::controlled_gate(qgates::phase_shift_gate(PI));
	const QuantumGate gateH = gate.H(); 

	// Test that gateH is the conjugate transpose
	for (int i = 0; i < pow(2, gate.N()); i++) { 
		for (int j = 0; j < pow(2, gate.N()); j++) { 
			EXPECT_EQ(conj(gate(i,j)), gateH(j,i)); 
		}
	}

	// Test that gateH is the inverse of the original gate
	const int state_index = rdunif(0, pow(2, gate.N()) - 1); 
	QubitSystem q(gate.N(), state_index); 
	const string state_orig = q.smeasure(); 

	gate * q; 
	gateH * q; 
	EXPECT_EQ(state_orig, q.smeasure()); 

	// Test applying gates in the reverse order
	gateH * q; 
	gate * q; 
	EXPECT_EQ(state_orig, q.smeasure()); 
}

// Tests the action of adding two quantum gates. Tests both the versions 
// G += H and G + H.
TEST (QuantumGateTest, TestSum) {
	// Test that |0><0| (x) I + |1><1| (x) X equals the controlled-not gate
	const QubitSystem q0(1); 
	const QubitSystem q1(1, 1); 
	const QuantumGate I(1); 
	const QuantumGate X = qgates::x_gate(); 
	const QuantumGate CN = qgates::cnot_gate(); 

	static const complex<double> mat00[] = {1., 0., 0., 0.}; 
	QuantumGate gate00(1, mat00);
	static const complex<double> mat11[] = {0., 0., 0., 1.}; 
	QuantumGate gate11(1, mat11);
	
	QuantumGate g1 = kron(gate00, I); 
	QuantumGate g2 = kron(gate11, X); 
	g1 += g2; 

	for (int i = 0; i < 4; i++) { 
		for (int j = 0; j < 4; j++) {
			EXPECT_LT(abs(g1(i,j) - CN(i,j)), EPS); 
		}
	}

	static const double THETA = PI / 3.;
	static const complex<double> thetai(0., THETA);
	// Test that phase shift gate plus Z gate equals this matrix
	static const complex<double> MATRIX[][2] = {{2., 0.},
	                                            {0., exp(thetai) - 1.}};  
	const QuantumGate R = qgates::phase_shift_gate(THETA); 
	const QuantumGate Z = qgates::z_gate(); 
	const QuantumGate g_sum = R + Z; 
	
	for (int i = 0; i < 2; i++) { 
		for (int j = 0; j < 2; j++) {
			EXPECT_LT(abs(g_sum(i,j) - MATRIX[i][j]), EPS); 
		}
	}
}


// Tests the operation of taking the matrix product of two quantum gates. 
TEST (QuantumGateTest, TestMatrixProduct) {
	// The product of the Y and Z gates is the gate with this matrix
	static const complex<double> YZ_MATRIX[][2] = {{ 0.,       {0., 1.}},
	                                               {{0., 1.},  0.}};
	QuantumGate Y = qgates::y_gate(); 
	const QuantumGate Z = qgates::z_gate(); 
	const QuantumGate YZ = Y * Z; 
	Y *= Z; // Test *= form as well
	static const int N = 1; 

	// Test correct dimension
	EXPECT_EQ(N, YZ.N()); 
	EXPECT_EQ(N, Y.N()); 

	for (int i = 0; i < pow(2, N); i++) {
		for (int j = 0; j < pow(2, N); j++) {
			EXPECT_LT(abs(YZ(i,j) - YZ_MATRIX[i][j]), EPS); 
			EXPECT_LT(abs(Y(i,j) - YZ_MATRIX[i][j]), EPS); 
		}
	}

	// Test power operation, i.e. multiplying a gate with itself multiple times
	QuantumGate X = qgates::x_gate(); 
	QuantumGate X4 = X^4; 
	X ^= 4; // test ^= form too
	double expected; 

	// X^4 is the identity
	for (int i = 0; i < pow(2, X4.N()); i++) {
		for (int j = 0; j < pow(2, X4.N()); j++) {
			expected = i == j; 
			EXPECT_LT(abs(X(i,j) - expected), EPS); 
			EXPECT_LT(abs(X4(i,j) - expected), EPS); 
		}
	}
}



// Test the tensor product of two quantum gates.
TEST (QuantumGateTest, TestTensorProduct) {
	// The product of the Y and Z gates is the gate with this matrix
	static const complex<double> YZ_MATRIX[][4] = {{0., 0., {0., -1.}, 0.},
	                                               {0., 0., 0., {0., 1.}},
	                                               {{0., 1.}, 0., 0., 0.},
	                                               {0., {0., -1.}, 0., 0.}};
	QuantumGate Y = qgates::y_gate(); 
	const QuantumGate Z = qgates::z_gate(); 
	const QuantumGate YZ = kron(Y, Z); 
	static const int N = 2; 

	// Test correct dimension
	EXPECT_EQ(N, YZ.N()); 

	for (int i = 0; i < pow(2, N); i++) {
		for (int j = 0; j < pow(2, N); j++) {
			EXPECT_LT(abs(YZ(i,j) - YZ_MATRIX[i][j]), EPS); 
		}
	}

	// Test power operation, i.e. tensoring a gate with itself multiple times
	QuantumGate X = qgates::x_gate(); 
	QuantumGate X3 = tensor_pow(X, 3); 
	double expected; 

	// X^3 is all zeros except for ones on the reverse diagonal
	for (int i = 0; i < pow(2, X3.N()); i++) {
		for (int j = 0; j < pow(2, X3.N()); j++) {
			expected = j == pow(2, X3.N()) - i - 1;
			EXPECT_LT(abs(X3(i,j) - expected), EPS); 
		}
	}
}




