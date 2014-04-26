// algorithms.cpp
// Implements quantum algorithms.

#include "algorithms.h"
#include "gate_factory.h"
#include "quantum.h"
#include "util.h"
#include <cmath>
#include <complex>
#include <stdexcept>
using std::string;
using std::complex; 

namespace quantum_algorithm_simulator {

static QuantumGate search_oracle (const string &match_text, const string *list, 
                                  int n);
static QuantumGate function_oracle (int (*f) (int), int y, int n);  
static int r (int n); 


// Performs the quantum Fourier transform on a system of qubits.
void qft (QubitSystem *q) {
	// Make quantum gate F, the discrete Fourier transform matrix of dimension 
	// 2^n, with a unitary normalization constant of 1/sqrt(N). 
	const int N = pow(2, q->N());
	const complex<double> TWOPIIN(0., 2. * PI / N); // Note 1
	const complex<double> omega = exp(TWOPIIN); 
	complex<double> **matrix = new_cmat(N, N); 

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			matrix[i][j] = pow(omega, i * j) / sqrt(N);
		}
	}
	QuantumGate F(N, matrix, true); 

	// Apply F to q
	F * (*q);
}


// Grover's quantum algorithm for searching. 
// Searches an array of strings of length n to find a string matching the given 
// string. 
// Returns the index i of the item such that list[i] is the same string as 
// match_text. 
// The lines must match exactly; it is not enough for the text to be contained 
// in the line. If two or more lines match, it will only return one of the line 
// numbers. Returns -1 if no matching line is found, i.e. the algorithm fails to
// find a solution. 
int grover_search (const string &match_text, const string *list, int n) {
	if (n <= 0) throw std::invalid_argument("n must be positive");

	const int N = ceil(log2(n)); // number of qubits needed
	const QuantumGate Hn = qgates::hadamard_gate(N);
	const QuantumGate Ui = search_oracle(match_text, list, n);
	const QuantumGate Ud = qgates::grover_diffusion_operator(N); 

	int index; 
	static const int MAX_ITER = 20; 
	int count = 0; 
	// Repeat until a solution is found or the iteration limit is reached
	do {
		QubitSystem q(N); // system of log2(n) bits in state |0>
		// apply Hadamard gate to create uniform superposition of all basis states
		Hn * q; 

		for (int i = 0; i < r(pow(2., N)); i++) { 
			Ui * q; // apply operator that flips the sign of the matching index
			Ud * q; // apply Grover's diffusion operator 
		}

		index = q.measure(); 
		count++;
	} while (count < MAX_ITER && (index >= n || list[index] != match_text)); 

	return index < n && list[index] == match_text ? index : -1; 
}


// Grover's algorithm for inverting a general function f that maps a sequence of
// n bits (represented as an int whose binary representation is the bit
// sequence) to another sequence of bits. Given f and a value y, returns the 
// input x such that f(x) = y. If more than one input suffices, it returns one
// at random. If no input suffices, returns -1. 
int grover_invert (int (*f) (int), int y, int n) { 
	if (n <= 0) throw std::invalid_argument("n must be positive");

	const QuantumGate Hn = qgates::hadamard_gate(n);
	const QuantumGate Ui = function_oracle(f, y, n);
	const QuantumGate Ud = qgates::grover_diffusion_operator(n); 

	int x; 
	static const int MAX_ITER = 20; 
	int count = 0; 
	// Repeat until a solution is found or the iteration limit is reached
	do {
		QubitSystem q(n); // system of n bits in state |0>
		// apply Hadamard gate to create uniform superposition of all basis states
		Hn * q; 

		for (int i = 0; i < r(pow(2., n)); i++) { 
			Ui * q; // apply operator that flips the sign of the matching input
			Ud * q; // apply Grover's diffusion operator 
		}

		x = q.measure(); 
		count++;
	} while (count < MAX_ITER && f(x) != y); 

	return f(x) == y ? x : -1; 
}

// Creates the "quantum oracle/black box" gate used by the Grover search 
// algorithm. The gate maps a state |k> to  -|k> if the kth item of "list" 
// matches the text match_text, and maps it to itself otherwise.  
static QuantumGate search_oracle (const string &match_text, const string *list, 
                                  int n) {
	QuantumGate gate(ceil(log2(n))); // identity gate

	for (int i = 0; i < n; i++) { 
		if (list[i] == match_text) { 
			gate.set(i, i, -1.); 
		}
	}

	return gate; 
}


// Creates the "quantum oracle/black box" gate used by the Grover function 
// inversion algorithm. The gate maps a state |x> to  -|x> if f(x) = y, and 
// maps it to itself otherwise.  
static QuantumGate function_oracle (int (*f) (int), int y, int n) { 
	QuantumGate gate(n); // identity gate

	for (int i = 0; i < n; i++) { 
		if (f(i) == y) {
			gate.set(i, i, -1.); 
		}
	}

	return gate; 
}


// Function returning the optimal number of iterations for the algorithm. It is
// important to stop at exactly this many iterations or else future iterations
// may actually lower the probability of measuring the correct answer.  
// See http://www.quantiki.org/wiki/Grover's_search_algorithm.
static int r (int n) { 
	double theta = asin(1. / sqrt(n)); 
	return (int) round((PI / theta - 2.) / 4.); 
}


} // end namespace 


// NOTES
// 1. This constant must be computed outside of the function exp(), to avoid 
// dumb errors with the complex exponential function. 

