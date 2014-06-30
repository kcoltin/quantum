// Tests for the functionality of the QuantumGate class.

#include "gate_factory.h" 
#include "quantum_gate.h" 
#include "qubit_system.h" 
#include <gtest/gtest.h>
using std::complex;
using std::string;
using namespace quantum_algorithm_simulator; 

// Tolerance for floating-point comparison
static const double EPS = 1.e-6; 

// Tests the factory method for producing a Pauli X-gate (aka "not" gate). This
// is sufficient for testing all the "basic" factory methods since their 
// internal workings are identical except for the values of the particular 
// matrix. 
TEST (GateFactoryTest, TestXGate) {
	static const int N = 1; 
	QubitSystem q(N); 
	QuantumGate X = qgates::x_gate(); 

	// Test correct dimension
	EXPECT_EQ(N, X.N()); 

	// Test that applying the X-gate to a qubit in state |0> converts it to 
	// state |1>, and vice versa
	X * q;  
	int *state = q.ameasure(); 
	EXPECT_EQ(1, state[0]); // the state should be |1>
	delete [] state; 

	X * q; 
	state = q.ameasure(); 
	EXPECT_EQ(0, state[0]); // the state should be |0>
	delete [] state; 
}


// Test Hadamard gate
TEST (GateFactoryTest, TestHadamardGate) {
	static const int N = 3; 
	// Matrix with 1's for all positive entries of three-qubit Hadamard gate
	static const int is_pos[] = {1, 1, 1, 1, 1, 1, 1, 1, 
	                             1, 0, 1, 0, 1, 0, 1, 0,
	                             1, 1, 0, 0, 1, 1, 0, 0,
	                             1, 0, 0, 1, 1, 0, 0, 1, 
	                             1, 1, 1, 1, 0, 0, 0, 0, 
	                             1, 0, 1, 0, 0, 1, 0, 1, 
	                             1, 1, 0, 0, 0, 0, 1, 1, 
	                             1, 0, 0, 1, 0, 1, 1, 0}; 
	const QuantumGate H = qgates::hadamard_gate(N); 
	const double ENTRY_VAL = pow(2., - N / 2.); // abs val of each entry 
	int sign; 
	for (int i = 0; i < pow(2, N); i++) { 
		for (int j = 0; j < pow(2, N); j++) { 
			sign = is_pos[(int) (i * pow(2, N) + j)] == 1 ? 1. : -1.;   
			EXPECT_LT(abs(H(i,j) - sign * ENTRY_VAL), EPS); 
		}
	}
}


// Tests the factory method for producing a controlled gate. 
TEST (GateFactoryTest, TestControlledGate) {
	// Gate that permutes the coefficients on a 2-qubit system
	static const int N = 2; 
	static const complex<double> matrix[] = {0., 0., 0., 1.,
	                                         1., 0., 0., 0.,
	                                         0., 1., 0., 0.,
	                                         0., 0., 1., 0.};
	const QuantumGate P(N, matrix);
	const QuantumGate CG = qgates::controlled_gate(P); 

	// Systems of three qubits, in states |011> and |101>
	QubitSystem q1(N + 1, 3); 
	QubitSystem q2(N + 1, 5); 
	QubitSystem q3(N + 1, 5); 

	// The controlled gate permutes the second two bits iff the first bit is one.

	// Test that applying the gate when the first bit is zero has no effect
	CG * q1;  
	EXPECT_EQ("011", q1.smeasure());

	// Test that the gate permutes the coefficients when the first bit is one, 
	// sending the coefficient of |101> to |110>
	CG * q2;  
	EXPECT_EQ("110", q2.smeasure()); 

	// Test Toffoli gate, a.k.a. controlled-controlled-not gate
	QuantumGate T = qgates::toffoli_gate(); 
	// Test that the gate has no effect on |011> and sends |110> to |111> 
	T * q1; 
	EXPECT_EQ("011", q1.smeasure()); 
	T * q2; 
	EXPECT_EQ("111", q2.smeasure()); 

	// Test Fredkin gate, a.k.a. controlled swap gate
	QuantumGate F = qgates::fredkin_gate(); 
	// Test that the gate has no effect on |011> and sends |101> to |110> 
	F * q1; 
	EXPECT_EQ("011", q1.smeasure()); 
	F * q3; 
	EXPECT_EQ("110", q3.smeasure()); 
}


// Tests the factory method for producing a gate that implements a classical 
// function. 
// Function f shifts the bit sequence x to the right by two, equivalent to 
// dividing by 4 and rounding down to the nearest integer.  
static int f (int x) { int y = x >> 2; return y; } 
TEST (GateFactoryTest, TestFunctionGate) {
	static const int M = 4;
	static const int K = 2;
	static const int N = M + K; 
	QuantumGate Uf = qgates::function_gate(f, M, K); 

	// Initial state - M-bit binary string for 6 plus random K-bit register 
	static const string state0 = "011010"; 
	QubitSystem q(N, bin_to_int(state0)); 

	// Expected output state - first M bits unchanged, followed by initial string
	// added to f(x) mod 2 (where x = first M bits)
	Uf * q; 
	EXPECT_EQ("011011", q.smeasure()); 
}






