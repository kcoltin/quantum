// quantum.cpp
// 
// Contains functions and data structures for simulating a quantum computer.

#include "quantum.h" 
#include "util.h"
#include <cmath>
#include <stdexcept>
using std::complex; 
using std::string;
using std::length_error;
using std::runtime_error;
using std::invalid_argument;

namespace quantum_algorithm_simulator {

// The constructors just call the corresponding initializer
QubitSystem::QubitSystem () {
	this->coeffs = NULL; 
	this->init(); 
}

QubitSystem::QubitSystem (int n) {
	this->coeffs = NULL; 
	this->init(n); 
}

QubitSystem::QubitSystem (int n, int state) {
	this->coeffs = NULL; 
	this->init(n, state); 
}

QubitSystem::QubitSystem (int n, const string state) { 
	this->coeffs = NULL; 
	this->init(n, state); 
}

// Default initializer. Creates a single qubit in state |0>.
void QubitSystem::init () {
	this->n = 1; 
	delete [] this->coeffs; 
	this->coeffs = new complex<double>[2]; 
	this->collapse(0); // set system to pure state |0>
}

// Initializer for an n-qubit system. The coefficients are initialized to 
// [1, 0, 0, ..., 0]. 
void QubitSystem::init (int n) {
	this->n = n; 
	delete [] this->coeffs; 
	this->coeffs = new complex<double>[(int) pow(2, n)]; 
	this->collapse(0); // set system to pure state |0000....0>
}

// Initializer for an n-qubit system. The coefficients are initialized such that
// the system is in the pure state given by state. 
void QubitSystem::init (int n, int state) {
	this->n = n; 
	delete [] this->coeffs; 
	this->coeffs = new complex<double>[(int) pow(2, n)]; 
	this->collapse(state); // set system to pure state 
}

// Initializer for an n-qubit system. The coefficients are initialized such that
// the system is in the pure state represented by the binary string STATE. 
void QubitSystem::init (int n, const string state) { 
	this->n = n; 
	delete [] this->coeffs; 
	this->coeffs = new complex<double>[(int) pow(2, n)]; 
	this->collapse(bin_to_int(state)); 
}

// Destructor
QubitSystem::~QubitSystem () {
	delete [] this->coeffs;
}


// Returns the number of qubits in the system
int QubitSystem::N () const { 
	return this->n; 
}


// Measures the system, returning the integer whose binary representation 
// is the sequence of measured qubits. E.g., if |0110> is measured, returns 6.
int QubitSystem::measure () {
	int state = this->get_observed_state(); 
	this->collapse(state); 
	return state;
}


// Returns the observed value resulting from measuring the ith qubit in the
// system. For example, in the system 0.5|00> + 0.5|01>, measure(2) would return
// 0 or 1, each with 50% probability. 
// Measuring the system necessarily collapses it to a linear combination of pure
// states that are compatible with the observation. 
int QubitSystem::measure (int bit_index) {
	int state = this->get_observed_state();
	int bit = get_bit(state, this->n - bit_index); 

	this->collapse(bit_index, bit); 
	return bit;
}


// Like measure(), but returns the state as a string. E.g. 1 (from measure()) or
// {0, 1} (from ameasure()) would become "01". 
string QubitSystem::smeasure () {
	int state = this->measure(); 
	return int_to_bin(state, this->n); 
}


// Returns a vector of observed values resulting when all qubits in a system are
// measured. For example, in the system 0.5|00> + 0.5|01>, measure() would 
// return {0, 0} or {0, 1}, each with 50% probability. 
// Measuring the system necessarily collapses it to the observed pure state.
int * QubitSystem::ameasure () {
	int state = this->measure(); 

	int *astate = new int[this->n];
	for (int i = 0; i < this->n; i++) { 
		astate[i] = get_bit(state, this->n - i -1);  
	}

	return astate;
}

// Alternate notation for gate.act(qubits)
void operator* (const QuantumGate &gate, QubitSystem &qubits){
	gate.act(&qubits); 
} 


// Computes the tensor product of two systems of qubits. 
// On return, the system q1 is overwritten with the tensor product of q1 and q2.
// The system q2 is "erased" in that it is resized to a system of size zero. 
// This represents the fact that the two systems are combined into a single 
// system of separable states. 
void operator* (QubitSystem &q1, QubitSystem &q2) {
	// Number of qubits in the product system
	int N = q1.n + q2.n; 
	// coefficients of the product system
	complex<double> *coeffs = new complex<double>[(int) pow(2, N)]; 
	int index; 

	// The coefficients of the tensor product are essentially the Kronecker 
	// product of the two vectors q1.coeffs and q2.coeffs thought of as n x 1 
	// matrices. 
	for (int i = 0; i < pow(2, q1.n); i++) {
		for (int j = 0; j < pow(2, q2.n); j++) {
			index = (int) i * pow(2, q2.n) + j;
			coeffs[index] = q1.coeffs[i] * q2.coeffs[j]; 
		}
	}

	// Move coeffs to be the coefficients of q1 
	q1.n = N; 
	delete [] q1.coeffs; 
	q1.coeffs = coeffs; 

	// Erase q2
	q2.n = 0; 
	delete [] q2.coeffs; 
	q2.coeffs = NULL; // avoids memory errors when destructor is called
}


// Copy constructor. This is private so as to not allow programs to violate
// the no-cloning theorem. 
QubitSystem::QubitSystem (const QubitSystem &orig) {
	this->n = orig.n; 
	this->coeffs = copyv(orig.coeffs, pow(2, this->n)); 
}

// Copy assignment operator. This is private so as to not allow programs to 
// violate the no-cloning theorem.
QubitSystem & QubitSystem::operator= (const QubitSystem &orig) {
	// Do nothing if this == orig
	if (this == &orig) return *this;

	// If the n's are unequal, delete and reallocate coeffs
	if (this->n != orig.n) {
		this->n = orig.n; 
		delete [] this->coeffs; 
		this->coeffs = new complex<double>[(int) pow(2, this->n)]; 
	}

	// Copy coefficients in place
	copyv(orig.coeffs, this->coeffs, (int) pow(2, this->n)); 
	
	return *this; 
}

// Collapses the system of qubits to a single pure state. It sets the 
// coefficients on the state given by state to one, and all other 
// coefficients to zero.
void QubitSystem::collapse (int state) {
	for (int i = 0; i < pow(2, this->n); i++) {
		this->coeffs[i] = 0.;
	}
	this->coeffs[state] = 1.;
}


// Collapses a system of qubits to those states where the bit given by bit_index
// has the given observed value. E.g., if q is a 2-qubit system, then calling
// q.collapse(2, 1) would collapse it down to a combination of the basis states
// |01> and |11>. 
void QubitSystem::collapse (int bit_index, int observed_value) {
	// Interval over which the values 0 and 1 alternate for the bit given by 
	// bit_index. For example, if bit_index = 2 and n = 3, then interval would
	// equal 2, because the bit is zero for 2 states, then one for 2 states, then
	// zero for 2, and then one for 2: 
	// |000> |001> |010> |011> |100> |101> |110> |111>
	int interval = pow(2, this->n - bit_index); 

	// Set the probabilities of all states where the bit given by bit_index is NOT
	// equal to observed_value to zero.
	// Iterate over the first state in each "block" of consecutive states where
	// the bit given by bit_index is NOT equal to observed_value.  
	// Note: the operator "observed_value ^ 1" maps 0 to 1 and vice versa. 
	for (int i = (observed_value ^ 1) * interval; i < this->n; i += interval * 2){
		// Then, interate over each state within this block
		for (int j = 0; j < interval; j++) {
			this->coeffs[i+j] = 0.;
		}
	}

	// Then, normalize the remaining states such that their squared coefficients
	// sum to one. 
	normalize (this->coeffs, this->n); 
}


// Randomly observes a particular state of the system, based on the probability
// magnitudes. Returns the index (between 0 and 2^n-1) of the state that is
// observed. 
int QubitSystem::get_observed_state () {
	double randval = runif(); 
	double total_prob = pow(abs(this->coeffs[0]), 2); 
	int index = 0; 
	while (total_prob < randval && index < pow(2, this->n) - 1) {
		index++; 
		total_prob += pow(abs(this->coeffs[index]), 2);
	}

	return index; 
}




// Default constructor for a quantum gate operating on n qubits. Initialized to
// the identity matrix. 
QuantumGate::QuantumGate (int n) {
	this->n = n; 
	this->matrix = eyec(pow(2, n), pow(2, n)); 
}


// Constructs a new quantum gate operating on an n-qubit system, where the 
// unitary matrix describing the action of the gate is given. Assigns the matrix
// by reference. 
QuantumGate::QuantumGate (int n, complex<double> **matrix) {
	this->n = n; 
	this->matrix = matrix; 
}


// Constructs a new quantum gate operating on an n-qubit system. If byref is
// false, it assigns the given matrix by reference rather than copying its 
// entries by value. 
QuantumGate::QuantumGate (int n, complex<double> **matrix, bool byref) {
	this->n = n; 
	this->matrix = byref ? matrix : copym(matrix, pow(2, n), pow(2, n)); 
} 


// Constructs a new quantum gate operating on an n-qubit system, where the 
// unitary matrix describing the action of the gate is given by a vector of
// values in row-major order.  
QuantumGate::QuantumGate (int n, const complex<double> *vector) {
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
// the product, a) is more realistic, and b) prevents cloning via multiplying a
// qubit by (for example) the identity gate. 
void QuantumGate::act (QubitSystem *q, int index) const {
	// TODO: for now, assumes index == 1 and q->n == this->n; change that

	// Throw error if the lengths differ
	if (q->N() != this->n) {
		throw length_error("Cannot operate on qubit system of different dimension");
	}
	// Throw error if matrix is not unitary
	if (!is_unitary(this->matrix, pow(2, this->n))) {
		throw runtime_error("Cannot perform non-unitary operation on a qubit");
	}

	prod_ip(this->matrix, q->coeffs, pow(2, this->n)); 
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



// Returns the tensor product, i.e. the Kronecker product, of two quantum gates.
QuantumGate & QuantumGate::operator%= (const QuantumGate &other) {
	complex<double> **matrix = kron(this->matrix, other.matrix, pow(2, this->n),
	                           pow(2, this->n), pow(2, other.n), pow(2, other.n));
	deletem(this->matrix, pow(2, this->n)); 
	this->n += other.n; 
	this->matrix = matrix; 
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
QuantumGate operator% (const QuantumGate &g1, const QuantumGate &g2) {
	QuantumGate product(g1.n + g2.n, kron(g1.matrix, g2.matrix, pow(2, g1.n), 
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
		product = kron(tmp, g.matrix, nprod, nprod, pow2n, pow2n); 
		deletem(tmp, nprod); 
		nprod *= pow2n; 
	}

	QuantumGate g_out(g.n * e, product); 
	return g_out; 
} 

} 
