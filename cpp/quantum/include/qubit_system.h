#ifndef QUBIT_SYSTEM_H
#define QUBIT_SYSTEM_H

#include "quantum_gate.h"
#include <complex>

namespace quantum_algorithm_simulator {

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
	QubitSystem (int n=1, int state=0); 
	void init (int n=1, int state=0); 
	~QubitSystem (); 

	int N () const; 
	int measure (); 
	int measure (int bit_index); 
	std::string smeasure (); 
	int * ameasure (); 

	friend void QuantumGate::act (QubitSystem *q) const; 
	friend void QuantumGate::act (QubitSystem *q, int index) const; 
	friend void operator* (const QuantumGate &gate, QubitSystem &qubits);  
	friend void operator* (QubitSystem &q1, QubitSystem &q2);
}; 


int bin_to_int (const std::string bin); 


}

#endif 



