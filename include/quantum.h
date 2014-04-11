#ifndef QUANTUM_H
#define QUANTUM_H

#include <complex>

class QuantumGate; 

// Represents a system of one or more qubits. 
//
// The number of qubits is given by N, and the coefficients on each of the
// basis states are given by coeffs. For example, a 2-qubit system would be of
// the form 
// a|00> + b|01> + c|10> + d|11>.
// In this case N = 2 and coeffs would be the array {a, b, c, d}. 
//
// Indices of bits are always one-based, but indices of states are always zero-
// based. 
//
// The indices are listed in coeffs as a vector, corresponding to the entries
// of an N-dimensional array stored in row-major order. That is,
// the coefficient on the basis state |b_1,b_2,b_3,...,b_N> (where each b is 
// 0 or 1) is given by 
// coeffs[N*b_1 + (N-1)*b_2 + (N-2)*b_3 + ... + 2*b_N-1 + b_N].  
// Equivalently, the state whose coefficient is given by coeffs[i] is the state
// whose representation in 0's and 1's is the binary representation of the 
// integer i: e.g., coeffs[3] will always represent the state |00...011> (where
// the number of leading zeros depends on N). 
class QubitSystem {
	private:
	int n; // Number of qubits in the system
	std::complex<double> *coeffs; // Coefficients on each of the basis states 

	QubitSystem (const QubitSystem &orig);
	QubitSystem & operator= (const QubitSystem &orig); 
	int get_observed_state (); 
	void collapse (int state); 
	void collapse (int bit_index, int observed_value); 

	public: 
	QubitSystem (); 
	QubitSystem (int n); 
	QubitSystem (int n, int state); 
	QubitSystem (int n, const std::string state); 
	void init (); 
	void init (int n); 
	void init (int n, int state); 
	void init (int n, const std::string state); 
	~QubitSystem (); 

	int N () const; 
	int measure (); 
	int measure (int bit_index); 
	std::string smeasure (); 
	int * ameasure (); 

	friend void operator* (const QuantumGate &gate, QubitSystem &qubits);  
	friend void operator* (QubitSystem &q1, QubitSystem &q2);
}; 



// Represents a quantum logic gate.
class QuantumGate {
	private:
	// Number of qubits that the gate operates on. A gate operating on n qubits is
	// represented as a 2^n x 2^n matrix. 
	int n; 
	// unitary matrix representing the action of the gate
	std::complex<double> **matrix; 

	public:
	QuantumGate (int n);
	QuantumGate (int n, std::complex<double> **matrix); 
	QuantumGate (int n, std::complex<double> **matrix, bool byref); 
	QuantumGate (int n, const std::complex<double> *vector); 
	QuantumGate (const QuantumGate &orig);
	QuantumGate & operator= (const QuantumGate &orig); 
	~QuantumGate (); 
	
	int N () const; 
	void set (int i, int j, std::complex<double> val); 
	std::complex<double> operator() (int i, int j) const; 
	QuantumGate H () const; 
	QuantumGate & operator+= (const QuantumGate &other); 
	QuantumGate & operator*= (const QuantumGate &other); 
	QuantumGate & operator^= (int e); 
	QuantumGate & operator%= (const QuantumGate &other); 

	friend QuantumGate operator+ (const QuantumGate &g1, const QuantumGate &g2); 
	friend void operator* (const QuantumGate &gate, QubitSystem &qubits);  
	friend QuantumGate operator* (const QuantumGate &g1, const QuantumGate &g2); 
	friend QuantumGate operator^ (const QuantumGate &g1, int e); 
	friend QuantumGate operator% (const QuantumGate &g1, const QuantumGate &g2); 
	friend QuantumGate tensor_pow (const QuantumGate &g, int e); 
};


#endif 



