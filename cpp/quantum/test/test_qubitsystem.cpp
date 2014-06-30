// Tests for the functionality of the QubitSystem class.

#include "gate_factory.h" 
#include "quantum_gate.h" 
#include "qubit_system.h" 
#include "util.h" 
#include <gtest/gtest.h>
using std::complex;
using std::string;
using namespace quantum_algorithm_simulator; 

// Tests measurement of a qubit system 
TEST (QubitSystemTest, TestMeasure) {
	// First, test measurement of a system in a pure state
	QubitSystem q(5, 10); // state |01010>
	static const string EXPECTED_STATE = "01010"; 
	EXPECT_EQ(EXPECTED_STATE, q.smeasure()); 
	int state = q.measure(); 
	int *astate = q.ameasure(); 
	for (int i = 0; i < q.N(); i++) { 
		EXPECT_EQ(EXPECTED_STATE.at(i) - '0', get_bit(state, q.N() - i - 1)); 
		EXPECT_EQ(EXPECTED_STATE.at(i) - '0', q.measure(i+1)); 
		EXPECT_EQ(EXPECTED_STATE.at(i) - '0', astate[i]); 
	}
	delete [] astate; 

	// Test probabilistic measurement: repeatedly measure superpositions of basis
	// states and test that the outcome varies randomly
	static const int N = 2; 
	QuantumGate H = qgates::hadamard_gate(N);
	static const int NTRIALS = 100; 
	int nsuccess1 = 0, nsuccess2 = 0; // number of 1's in the first and second bit
	for (int i = 0; i < NTRIALS; i++) { 
		q.init(N); // set the qubits to |00>
		H * q; // equal superposition of four basis states
		nsuccess1 += q.measure(1); 
		nsuccess2 += q.measure(2); 
	}
	// Test that repeated measurement gives the same result
	state = q.measure(); 
	astate = q.ameasure(); 
	EXPECT_EQ(get_bit(state, 1), q.measure(1)); 
	EXPECT_EQ(astate[0], q.measure(1)); 
	EXPECT_EQ(get_bit(state, 0), q.measure(2)); 
	EXPECT_EQ(astate[1], q.measure(2)); 

	// Test that the number of 1's that appeared in the first and second qubits 
	// is approximately half the total (within approx. 99% confidence interval) 
	static const double CUTOFF = NTRIALS * 0.15;  
	EXPECT_LE(abs(nsuccess1 - NTRIALS / 2.), CUTOFF); 
	EXPECT_LE(abs(nsuccess2 - NTRIALS / 2.), CUTOFF); 

	delete [] astate; 
}



// Tests computation of the tensor product of two qubit systems. 
TEST (QubitSystemTest, TestTensorProduct) {
	// Test that |1> (x) |0> = |10>
	static const int N = 1; 
	QubitSystem q1(N, 1); 
	QubitSystem q2(N); 

	q1 * q2; 
	EXPECT_EQ("10", q1.smeasure()); 

	// Test that |0100> (x) |101> = |0100101>
	static const int N1 = 4; 
	static const int N2 = 3; 
	q1.init(N1, 4); 
	q2.init(N2, 5); 

	q1 * q2; 
	EXPECT_EQ("0100101", q1.smeasure()); 
}



