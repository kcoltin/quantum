// quantum_gate.cpp

#include "quantum_gate.h" 
#include "qubit_system.h" 
#include "util.h"
#include <cmath>
#include <stdexcept>
using std::complex; 
using std::string;
using std::length_error;
using std::runtime_error;
using std::invalid_argument;

namespace quantum_algorithm_simulator {

// Constructs a new quantum gate operating on an n-qubit system, where the 
// unitary matrix describing the action of the gate is given. If byref is
// false, it assigns the given matrix by reference rather than copying its 
// entries by value. 
// By default, the matrix is the identity matrix and byref is true.
QuantumGate::QuantumGate (int n, complex<double> **matrix, bool byref) {
	if (n <= 0) throw invalid_argument("n must be positive"); 
	this->n = n; 
	if (matrix == NULL) {
		matrix = eyec(pow(2, n), pow(2, n));
		byref = true; // pointless to make a copy of a locally-allocated matrix
	}
	this->matrix = byref ? matrix : copym(matrix, pow(2, n), pow(2, n)); 
} 


// Constructs a new quantum gate operating on an n-qubit system, where the 
// unitary matrix describing the action of the gate is given by a vector of
// values in row-major order.  
QuantumGate::QuantumGate (int n, const complex<double> *vector) {
	if (n <= 0) throw invalid_argument("n must be positive"); 
	this->n = n; 
	this->matrix = vec_to_mat(vector, pow(2, n), pow(2, n));
}

// Copy constructor
QuantumGate::QuantumGate (const QuantumGate &orig) {
	this->n = orig.n; 
	this->matrix = copym(orig.matrix, pow(2, this->n)); 
}

// Copy assignment operator
QuantumGate & QuantumGate::operator= (const QuantumGate &orig) {
	// Do nothing if this == orig
	if (this == &orig) return *this;

	// If the n's are unequal, delete and reallocate matrix
	if (this->n != orig.n) {
		this->n = orig.n; 
		deletem(this->matrix, pow(2, this->n));  
		this->matrix = new_cmat(pow(2, this->n), pow(2, this->n)); 
	}

	// Copy coefficients in place
	copym(orig.matrix, this->matrix, (int) pow(2, this->n)); 

	return *this; 
}


// Destructor
QuantumGate::~QuantumGate () {
	deletem(this->matrix, pow(2, this->n)); 
}


// Returns the number of qubits the gate operates on
int QuantumGate::N () const { 
	return this->n; 
}

// Sets the i, j entry of the matrix to the specified value. 
void QuantumGate::set (int i, int j, complex<double> val) { 
	if (i < 0 || pow(2, this->n) <= i || j < 0 || pow(2, this->n) <= j) { 
		throw invalid_argument("index out of range"); 
	}

	this->matrix[i][j] = val; 
}



// Accesses the i,j entry of the matrix. Can only be used as an accessor - 
// cannot be used to set the value.
complex<double> QuantumGate::operator() (int i, int j) const {
	return this->matrix[i][j]; 
}


// Returns the conjugate transpose of the quantum gate, equivalent to the 
// inverse of the gate (since gate matrices are unitary). 
QuantumGate QuantumGate::H () const {
	QuantumGate gate(this->n, conj_transpose(this->matrix, pow(2, this->n)));
	return gate; 
}


// Overloaded/default version of act, applies to first n qubits where n is the
// number of qubits that q acts on. 
void QuantumGate::act (QubitSystem *q) const {
	this->act(q, 1); 
}


// Performs the action of a quantum gate on a qubit system.  
// The coefficients of the product are the result of the left matrix 
// multiplication of the gate by the system q, starting with the qubit in 
// position INDEX (where 1, not 0, is the first bit). E.g., if the gate G
// operates on 2 qubits, then G.act(q, 2) would operate on the 2nd and 3rd
// qubits of the system, e.g. on the 1 bits of a system in state |01100>. 
// 
// It operates IN PLACE on the argument q, so the qubit system is changed
// by operation of the gate. This setup, as opposed to returning a new system as
// the product, prevents violations of the no-cloning theorem. 
void QuantumGate::act (QubitSystem *q, int index) const {
	// Throw error if the dimensions don't fit
	if (index <= 0 || index + this->n - 1 > q->N()) { 
		throw invalid_argument("Dimension mismatch with gate and qubit system"); 
	}
	// Throw error if matrix is not unitary
	if (!is_unitary(this->matrix, pow(2, this->n))) {
		throw runtime_error("Cannot perform non-unitary operation on a qubit");
	}

	complex<double> **mat, **tmp, **I; 

	if (index == 1 && this->n == q->N()) { 
		prod_ip(this->matrix, q->coeffs, pow(2, this->n)); 
	}
	else if (index > 1 && index + this->n - 1 == q->N()) { 
		I = eyec(pow(2, index - 1), pow(2, index - 1)); 
		mat = kronecker(I, this->matrix, pow(2, index - 1), pow(2, index - 1), 
		                pow(2, this->n), pow(2, this->n)); 
		prod_ip(mat, q->coeffs, pow(2, q->N())); 

		deletem(mat, pow(2, q->N())); 
		deletem(I, pow(2, index - 1)); 
	}
	else if (index == 1 && this->n < q->N()) { 
		I = eyec(pow(2, q->N() - this->n), pow(2, q->N() - this->n)); 
		mat = kronecker(this->matrix, I, pow(2, this->n), pow(2, this->n), 
		                pow(2, q->N() - this->n), pow(2, q->N() - this->n)); 
		prod_ip(mat, q->coeffs, pow(2, q->N())); 

		deletem(mat, pow(2, q->N())); 
		deletem(I, pow(2, q->N() - this->n)); 
	}
	else { // if (index > 1 && index + this->n - 1 < q->N())  
		I = eyec(pow(2, index - 1), pow(2, index - 1)); 
		tmp = kronecker(I, this->matrix, pow(2, index - 1), pow(2, index - 1), 
		                pow(2, this->n), pow(2, this->n)); 

		deletem(I, pow(2, index - 1)); 

		I = eyec(pow(2, q->N() - (index + this->n - 1)), 
		         pow(2, q->N() - (index + this->n - 1)));
		mat = kronecker(tmp, I, pow(2, this->n + index - 1),
		                pow(2, this->n + index - 1),
		                pow(2, q->N() - (index + this->n - 1)), 
		                pow(2, q->N() - (index + this->n - 1))); 
		prod_ip(mat, q->coeffs, pow(2, q->N())); 

		deletem(mat, pow(2, q->N())); 
		deletem(tmp, pow(2, this->n + index - 1)); 
		deletem(I, pow(2, q->N() - (index + this->n - 1))); 
	}
}


// Sum of two quantum gates - just the matrix sum of the two matrices. 
QuantumGate & QuantumGate::operator+= (const QuantumGate &other) {
	// Throw error if the lengths differ
	if (this->n != other.n) {
		throw length_error("Cannot add quantum gates of unequal dimension");
	}

	sum_ip(this->matrix, other.matrix, pow(2, this->n)); 
	return *this; 
}


// Difference of two quantum gates - just the difference of the two matrices. 
QuantumGate & QuantumGate::operator-= (const QuantumGate &other) {
	// Throw error if the lengths differ
	if (this->n != other.n) {
		throw length_error("Cannot subtract quantum gates of unequal dimension");
	}

	diff_ip(this->matrix, other.matrix, pow(2, this->n)); 
	return *this; 
}

// Multiplies two gates via standard matrix multiplication. The action of the
// returned gate is equivalent to the action of the righthand gate followed by
// that of the lefthand gate.
QuantumGate & QuantumGate::operator*= (const QuantumGate &other) {
	// Throw error if the lengths differ
	if (this->n != other.n) { 
		throw length_error("Cannot multiply quantum gates of unequal dimension");
	}

	prod_ip(this->matrix, other.matrix, pow(2, this->n));
	return *this; 
}


// Returns the gate that results from taking the matrix product of the given
// gate with itself e times, that is, raising the matrix to the exponent 
QuantumGate & QuantumGate::operator^= (int e) {
	if (e < 0) throw invalid_argument("exponent must be nonnegative"); 
	
	mpow_ip(this->matrix, e, pow(2, this->n)); 
	return *this; 
}

// Sum of two quantum gates - just the matrix sum of the two matrices. 
QuantumGate operator+ (const QuantumGate &g1, const QuantumGate &g2) {
	// Throw error if the lengths differ
	if (g1.n != g2.n) {
		throw length_error("Cannot add quantum gates of unequal dimension");
	}

	QuantumGate g_out(g1.n, sum(g1.matrix, g2.matrix, pow(2, g1.n))); 
	return g_out;
}

// Difference of two quantum gates - just the difference of the two matrices. 
QuantumGate operator- (const QuantumGate &g1, const QuantumGate &g2) {
	// Throw error if the lengths differ
	if (g1.n != g2.n) {
		throw length_error("Cannot subtract quantum gates of unequal dimension");
	}

	QuantumGate g_out(g1.n, diff(g1.matrix, g2.matrix, pow(2, g1.n))); 
	return g_out;
}

// Multiplies two gates via standard matrix multiplication. The action of the
// returned gate is equivalent to the action of the righthand gate followed by
// that of the lefthand gate.
QuantumGate operator* (const QuantumGate &g1, const QuantumGate &g2) {
	// Throw error if the lengths differ
	if (g1.n != g2.n) { 
		throw length_error("Cannot multiply quantum gates of unequal dimension");
	}

	QuantumGate product(g1.n, prod(g1.matrix, g2.matrix, pow(2, g1.n)));
	return product;
}

// Returns the gate that results from taking the matrix product of the given
// gate with itself e times, that is, raising the matrix to the exponent 
QuantumGate operator^ (const QuantumGate &gate, int e) {
	if (e < 0) throw invalid_argument("exponent must be nonnegative"); 
	
	QuantumGate g_out(gate.n, mpow(gate.matrix, e, pow(2, gate.n))); 
	return g_out; 
}

	
// Returns the tensor product, i.e. the Kronecker product, of two quantum gates.
QuantumGate kron (const QuantumGate &g1, const QuantumGate &g2) {
	QuantumGate product(g1.n + g2.n, kronecker(g1.matrix, g2.matrix, pow(2, g1.n),
	                                   pow(2, g1.n), pow(2, g2.n), pow(2, g2.n)));
	return product; 
}


// Returns the gate that results from taking the tensor product of the given
// gate with itself e times, analogous to raising g to an exponent e.  
QuantumGate tensor_pow (const QuantumGate &g, int e) {
	if (e <= 0) throw invalid_argument("exponent must be positive"); 

	const int pow2n = pow(2, g.n); 
	int nprod = pow2n; // dimension of matrix "product"
	complex<double> **product = copym(g.matrix, nprod); 
	complex<double> **tmp;

	for (int i = 2; i <= e; i++) { 
		tmp = copym(product, nprod); 
		deletem(product, nprod); 
		product = kronecker(tmp, g.matrix, nprod, nprod, pow2n, pow2n); 
		deletem(tmp, nprod); 
		nprod *= pow2n; 
	}

	QuantumGate g_out(g.n * e, product); 
	return g_out; 
} 

} 
