// algorithms.cpp
// Implements quantum algorithms.

#include "algorithms.h"
#include "gate_factory.h"
#include "quantum_gate.h"
#include "qubit_system.h"
#include "util.h"
#include <cmath>
#include <complex>
#include <functional>
#include <stdexcept>
using std::string;
using std::complex; 
using std::runtime_error; 

namespace quantum_algorithm_simulator {

static QuantumGate search_oracle (const string &match_text, const string *list, 
                                  int n);
static QuantumGate function_oracle (int (*f) (int), int y, int n);  
static int r (int n); 
static int modexp (int b, int e, int n); 


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
	F.act(q);
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
	static const int MAX_ITER = 50; 
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
	static const int MAX_ITER = 50; 
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


// Uses Shor's algorithm to return a factor of the integer n. If n is composite,
// returns a non-trivial factor (i.e. a factor other than n or 1). Otherwise,
// if n is prime, returns 1.  
int shor_factor (int n) {
	if (n <= 1) throw runtime_error("n must be at least 2"); 
	if (is_prime(n)) return 1; 
	if (n % 2 == 0) return 2; // even numbers > 2 are trivial
	// Need to check that n is not the power of an integer for algorithm to work
	int root = int_root(n); 
	if (root != -1) return root; 

	// choose m s.t. n^2 <= 2^m < 2*n^2
	// log2(n^2) <= m <= log2(2 * n^2) = 1 + log2(n^2) 
	int m = ceil(log2(pow(n, 2)));

	const int ny = ceil(log2(n - 1)); // number of qubits in output of function f
	const QuantumGate I(ny); 
	const QuantumGate H = kron(qgates::hadamard_gate(m), I); 

	QubitSystem q; 
	int a, r; 
	const int MAX_ITER = 10; 
	int niter = 0; 

	do {
		a = n - 1; // arbitrary integer coprime to n

		// Initialize a system of qubits long enough to represent the integers 0 to
		// 2^m - 1 alongside an integer up to n - 1, then apply Hadamard gate to the
		// first m qubits to create a uniform superposition.
		q.init(m + ny); 
		H * q; 

		// Apply the function f(x) = a^x (mod n) to the system
		auto f = [a, n] (int x) { return modexp(a, x, n); }; 
		const QuantumGate Uf = qgates::function_gate(f, m, ny); 
		Uf * q; 

		// Find the period of f via quantum Fourier transform
		qft(&q); 
		r = 1. / q.measure(); // period = 1 / frequency

		niter++; 
	} while (niter < MAX_ITER && (r % 2 != 0 || ((int) pow(a, r / 2)) % n == -1));

	return gcd(pow(a, r / 2) + 1, n); 
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

	for (int i = 0; i < pow(2, n); i++) { 
		if (f(i) == y) {
			gate.set(i, i, -1.); 
		}
	}

	return gate; 
}


// Function returning the optimal number of iterations for Grover's algorithm.
// It is important to stop at exactly this many iterations or else future
// iterations may actually lower the probability of measuring the correct
// answer. See http://www.quantiki.org/wiki/Grover's_search_algorithm.
static int r (int n) { 
	double theta = asin(1. / sqrt(n)); 
	return (int) round((PI / theta - 2.) / 4.); 
}

// Performs modular exponentiation. Given integers b, e, and n, returns 
// y = b^e (mod n).
static int modexp (int b, int e, int n) {
	return (int) pow(b, e) % n; 
}


} // end namespace 


// NOTES
// 1. This constant must be computed outside of the function exp(), to avoid 
// dumb errors with the complex exponential function. 

