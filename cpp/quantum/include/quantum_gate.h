#ifndef QUANTUM_GATE_H
#define QUANTUM_GATE_H

#include <complex>

namespace quantum_algorithm_simulator {

class QubitSystem;

// Represents a quantum logic gate.
class QuantumGate {
	private:
	// Number of qubits that the gate operates on. A gate operating on n qubits is
	// represented as a 2^n x 2^n matrix. 
	int n; 
	// unitary matrix representing the action of the gate
	std::complex<double> **matrix; 

	public:
	QuantumGate (int n, std::complex<double> **matrix=NULL, bool byref=true); 
	QuantumGate (int n, const std::complex<double> *vector); 
	QuantumGate (const QuantumGate &orig);
	QuantumGate & operator= (const QuantumGate &orig); 
	~QuantumGate (); 
	
	int N () const; 
	void set (int i, int j, std::complex<double> val); 
	std::complex<double> operator() (int i, int j) const; 
	QuantumGate H () const; 
	void act (QubitSystem *q) const; 
	void act (QubitSystem *q, int index) const; 
	QuantumGate & operator+= (const QuantumGate &other); 
	QuantumGate & operator-= (const QuantumGate &other); 
	QuantumGate & operator*= (const QuantumGate &other); 
	QuantumGate & operator^= (int e); 

	friend void operator* (const QuantumGate &gate, QubitSystem &qubits);  
	friend QuantumGate operator+ (const QuantumGate &g1, const QuantumGate &g2); 
	friend QuantumGate operator- (const QuantumGate &g1, const QuantumGate &g2); 
	friend QuantumGate operator* (const QuantumGate &g1, const QuantumGate &g2); 
	friend QuantumGate operator^ (const QuantumGate &g1, int e); 
	friend QuantumGate kron (const QuantumGate &g1, const QuantumGate &g2); 
	friend QuantumGate tensor_pow (const QuantumGate &g, int e); 
};


}

#endif 



