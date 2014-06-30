// gate_factory.cpp 
// 
// Functions in namespace qgates: factory functions for creating quantum gates.

#include "gate_factory.h"
#include "quantum_gate.h"
#include "util.h"
#include <functional>
#include <stdexcept>
using std::complex; 

namespace quantum_algorithm_simulator {
namespace qgates {

// Factory method for Pauli-X gate (i.e. "not gate")
QuantumGate x_gate () {
	static const complex<double> vals[] = {0., 1., 
	                                       1., 0.}; 
	QuantumGate gate(1, vals); 
	return gate; 
}

// Factory method for Pauli-Y gate 
QuantumGate y_gate () {
	static const complex<double> vals[] = {0.,      {0., -1.}, 
	                                       {0, 1.}, 0.}; 
	QuantumGate gate(1, vals); 
	return gate; 
}

// Factory method for Pauli-Z gate 
QuantumGate z_gate () {
	static const complex<double> vals[] = {1., 0., 
	                                       0., -1.}; 
	QuantumGate gate(1, vals); 
	return gate; 
}

// Factory method for Hadamard ("square root of not") gate
QuantumGate hadamard_gate () {
	return hadamard_gate(1); 
}


// Hadamard gate repeated n times, to operate on n qubits. 
QuantumGate hadamard_gate (int n) {
	static const complex<double> vals[] = {1./sqrt(2), 1./sqrt(2), 
	                                       1./sqrt(2), -1./sqrt(2)}; 
	static const QuantumGate H(1, vals); 
	return tensor_pow(H, n); 
}


// Factory method for phase-shift gate
QuantumGate phase_shift_gate (double theta) {
	const complex<double> THETAI(0., theta); // Note 1 
	complex<double> vals[] = {1., 0., 
	                          0., exp(THETAI)}; 
	QuantumGate gate(1, vals); 
	return gate; 
}


// Factory method for swap gate
QuantumGate swap_gate () {
	static const complex<double> vals[] = {1., 0., 0., 0., 
	                                       0., 0., 1., 0., 
	                                       0., 1., 0., 0., 
	                                       0., 0., 0., 1.}; 
	QuantumGate gate(2, vals); 
	return gate; 
}

// Factory method for controlled not gate
QuantumGate cnot_gate () {
	static const complex<double> vals[] = {1., 0., 0., 0., 
	                                       0., 1., 0., 0., 
	                                       0., 0., 0., 1., 
	                                       0., 0., 1., 0.}; 
	QuantumGate gate(2, vals); 
	return gate; 
}

// Factory method for controlled-U gate, where U is an arbitrary gate operating
// on n qubits. 
// The controlled-U gate operates on n+1 qubits and is of the form
// |I 0|
// |0 U|
// where U is any gate operating on n-bits and I is the 2^n x 2^n identity 
// matrix. 
QuantumGate controlled_gate (const QuantumGate &U) {
	// Make n+1 bit identity gate, i.e. 2^(n+1) x 2^(n+1) identity matrix
	QuantumGate cgate(U.N() + 1); 
	int m = (int) pow(2, U.N()); // m := 2^n

	// Copy the matrix of U into the lower righthand quadrant of the new gate
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			cgate.set(m+i, m+j, U(i,j)); 
		}
	}

	return cgate; 
}


// Toffoli gate, a.k.a. controlled-controlled-not gate. 
QuantumGate toffoli_gate () {
	return controlled_gate(cnot_gate()); 
}


// Fredkin gate, a.k.a. controlled swap gate. 
QuantumGate fredkin_gate () {
	return controlled_gate(swap_gate()); 
}


// Factory method for producing a gate that implements a classical function f. 
// f should be a function that takes a sequence of bits of length m and returns 
// another sequence of bits of length k. The bit sequences are represented by 
// the integers they represent in binary: e.g. "00101" would be the integer 5.  
// 
// The resulting gate, Uf, takes as its input a qubit system of length m + k. 
// Call the first m qubits in the system "x" and the last k qubits "y": then 
// the action of Uf on a pure state |x, y> is to map it to |x, f(x) (+) y> where
// (+) represents bitwise addition mod 2 (equivalent to the XOR operation).  
// It is easily shown that the map |x, y> --> |x, f(x) (+) y> is a bijection, 
// which implies that Uf is a permutation matrix (and is therefore unitary).  
QuantumGate function_gate (int (*f) (int), int m, int k) {
	auto g = [f] (int x) { return f(x); }; 
	return function_gate(g, m, k); 
} 

// Version of function_gate taking a C++11 functional rather than a function 
// pointer. This version implements the actual functionality.  
QuantumGate function_gate (std::function<int (int)> f, int m, int k) {
	unsigned int n = m + k; // gate operates on an (m + k)-qubit system
	complex<double> **matrix = czeros(pow(2, n), pow(2, n)); 
	unsigned int x, y, i; 
	// Filter that zeroes out all but the last k bits 
	const unsigned int FILTER = (1 << k) - 1; 

	// Iterate over columns of matrix
	for (unsigned int j = 0; j < pow(2, n); j++) {
		y = j & FILTER; // y is the last k bits of the jth state
		x = j >> k; // x is state ignoring the last k bits
		
		// Gate maps |x, y> to |x, y (+) f(x)>, where (+) is bitwise addition mod 2
		i = (x << k) + (y ^ (unsigned int) f((int) x)); 
		matrix[i][j] = 1.; 
	}

	QuantumGate gate(n, matrix); 
	return gate; 
}


// Grover diffusion operator, operating on n qubits. 
// Every non-diagonal entry is 1/2^(n-1); diagonal entries are that value minus
// one. 
// Note: at least one source I've see defines the operator as the negative of
// this matrix; either definition works for Grover's algorithm.
QuantumGate grover_diffusion_operator (int n) {
	const double VAL = 1. / pow(2, n - 1); 
	QuantumGate gate(n); 

	for (int i = 0; i < pow(2, n); i++) { 
		for (int j = 0; j < pow(2, n); j++) { 
			gate.set(i, j, VAL - (i == j)); 
		}
	}

	return gate; 
} 


}
} 


/* NOTES
1. e^theta*i must be computed outside of the function exp(), to avoid dumb 
   errors with the complex exponential function. 


*/



