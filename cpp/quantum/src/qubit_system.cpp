// quantum.cpp
// 
// Contains functions and data structures for simulating a quantum computer.

#include "qubit_system.h" 
#include "quantum_gate.h" 
#include "util.h"
#include <cmath>
#include <stdexcept>
using std::complex; 
using std::string;
using std::length_error;
using std::runtime_error;
using std::invalid_argument;

namespace quantum_algorithm_simulator {

// Constructor just calls the initializer
QubitSystem::QubitSystem (int n, int state) {
	this->coeffs = NULL; 
	this->init(n, state); 
}


// Initializer for an n-qubit system. The coefficients are initialized such that
// the system is in the pure state given by state. 
//
// Defaults to a single-qubit system in state |0>.
//
// Examples of initializing a system in state |010>:
//   q.init(3, 2)
//   or
//   q.init(3, bin_to_int("010"))
void QubitSystem::init (int n, int state) {
	this->n = n; 
	delete [] this->coeffs; 
	this->coeffs = new complex<double>[(int) pow(2, n)]; 

	if (state >= pow(2, n)) throw invalid_argument("invalid state");
	this->collapse(state); // set system to pure state 
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


// Combines two systems of qubits. 
// On return, the coefficients of the system q1 are overwritten with the tensor 
// product of the coefficients of  q1 and q2.
// The system q2 is "erased" in that it is resized to a system of size zero. 
// This represents the fact that the two systems are combined into a single 
// system of separable states, and avoids violating the no-cloning theorem.  
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
	// Note that the operation "observed_value ^ 1" maps 0 to 1 and vice versa. 
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


// Converts a string of 0s and 1s to the integer it represents in binary form. 
// E.g. if str is the string "1010", then bin_to_int(str) would return 10. 
// This is here, rather than in util.cpp where it might seem to belong better,
// so that it can be publically visible and used to do things like
// q = QubitSystem(3, bin_to_int("011")).
int bin_to_int (string bin) {
	int n = 0; 
	int len = bin.length(); 
	int bit; 

	for (int i = 0; i < len; i++) { 
		bit = bin[len-i-1];  
		if (bit == '1') { 
			n += 1 << i; // equals 2^i
		} 
		else if (bit != '0') {
			throw invalid_argument("Argument must be string of 0s and 1s only");
		}
	}

	return n; 
} 


} 
