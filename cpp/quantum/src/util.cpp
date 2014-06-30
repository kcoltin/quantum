#include "util.h"
#include <cmath>
#include <stdexcept>
#include <gsl/gsl_rng.h>
#include <sys/time.h>
using std::complex; 
using std::string;
using std::invalid_argument; 

const double PI = 4. * atan(1.); 

// Static function declarations
static gsl_rng * init_runif ();


// Allocates an m x n complex matrix. 
complex<double> ** new_cmat (int m, int n) {
	complex<double> **A = new complex<double>*[m]; 
	for (int i = 0; i < m; i++) {
		A[i] = new complex<double>[n]; 
	}
	return A; 
} 


// Makes an m x n complex matrix of zeros. 
complex<double> ** czeros (int m, int n) {
	complex<double> **A = new complex<double>*[m]; 

	for (int i = 0; i < m; i++) {
		A[i] = new complex<double>[n]; 

		for (int j = 0; j < n; j++) {
			A[i][j] = 0.;
		}
	}

	return A; 
} 


// Makes an m x n complex identity matrix. 
complex<double> ** eyec (int m, int n) {
	complex<double> **I = new complex<double>*[m]; 

	for (int i = 0; i < m; i++) {
		I[i] = new complex<double>[n]; 

		for (int j = 0; j < n; j++) {
			I[i][j] = i == j; 
		}
	}

	return I; 
} 


// Makes a copy of an n-dimensional complex vector
complex<double> * copyv (const complex<double> *orig, int n) {
	complex<double> *vnew = new complex<double>[n]; 
	for (int i = 0; i < n; i++) {
		vnew[i] = orig[i];
	}
	return vnew; 
}


// Makes a copy of an n-dimensional complex vector, copying the contents of 
// "orig" into "vnew".
void copyv (const complex<double> *orig, complex<double> *vnew, int n) {
	for (int i = 0; i < n; i++) {
		vnew[i] = orig[i];
	}
}


// Makes a copy of an m x n complex matrix
complex<double> ** copym (complex<double> **orig, int m, int n) {
	complex<double> **mnew = new_cmat(m, n); 
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			mnew[i][j] = orig[i][j]; 
		}
	}
	return mnew;
}

// Overloaded version for a square matrix. 
complex<double> ** copym (complex<double> **orig, int n) {
	return copym(orig, n, n);
}

// Makes a copy of an n-dimensional complex vector, copying the contents of 
// "orig" into "mnew".
void copym (complex<double> **orig, complex<double> **mnew, int m, int n) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			mnew[i][j] = orig[i][j]; 
		}
	}
}

// Overloaded version for a square matrix. 
void copym (complex<double> **orig, complex<double> **mnew, int n) {
	copym(orig, mnew, n, n);
}





// Copies the values of a complex vector with entries in row-major order into a 
// matrix.
complex<double> ** vec_to_mat (const complex<double> *x, int m, int n) {
	complex<double> **A = new_cmat(m, n); 
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = x[n*i+j];
		}
	}
	return A; 
}



// Deletes memory for a complex matrix with leading dimension m. 
void deletem (complex<double> **A, int m) {
	for (int i = 0; i < m; i++) {
		delete [] A[i]; 
	}
	delete [] A;
}


// Returns the conjugate transpose of a matrix.
complex<double> ** conj_transpose (complex<double> **A, int m, int n) {
	complex<double> **B = new_cmat(n, m); 
	for (int i = 0; i < m; i++) { 
		for (int j = 0; j < n; j++) { 
			B[j][i] = conj(A[i][j]); 
		}
	}
	return B; 
}

// Overloaded version for a square n x n matrix. 
complex<double> ** conj_transpose (complex<double> **A, int n) {
	return conj_transpose(A, n, n); 
}

// Returns the sum of two n x n complex matrices. 
complex<double> ** sum (complex<double> **A, complex<double> **B, int n) {
	complex<double> **C = new_cmat(n, n); 
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) { 
			C[i][j] = A[i][j] + B[i][j]; 
		}
	}

	return C; 
}
 
// Adds two n x n complex matrices, operating on the first matrix in place. That
// is, on return the matrix A will be the sum of the original A and B.
void sum_ip (complex<double> **A, complex<double> **B, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) { 
			A[i][j] += B[i][j]; 
		}
	}
} 

// Returns the difference of two n x n complex matrices. 
complex<double> ** diff (complex<double> **A, complex<double> **B, int n) {
	complex<double> **C = new_cmat(n, n); 
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) { 
			C[i][j] = A[i][j] - B[i][j]; 
		}
	}

	return C; 
}

// Subtracts two n x n complex matrices, operating on the first matrix in place.
// That is, on return the matrix A will be the difference of the original A
// minus B.
void diff_ip (complex<double> **A, complex<double> **B, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) { 
			A[i][j] -= B[i][j]; 
		}
	}
} 



// Inner product of two complex vectors; equivalent to the matrix product 
// x'y, where x' is the conjugate transpose. 
complex<double> inner_prod (const complex<double> *x, const complex<double> *y, 
                            int n) {
	complex<double> ip = 0.; 
	for (int i = 0; i < n; i++) { 
		ip += conj(x[i]) * y[i]; 
	}
	return ip; 
}

// Outer product of two complex vectors; equivalent to the matrix product 
// xy', where y' is the conjugate transpose. 
complex<double> ** outer_prod (const complex<double> *x, 
                               const complex<double> *y, int n) {
	complex<double> **op = new_cmat(n, n); 

	for (int i = 0; i < n; i++) { 
		for (int j = 0; j < n; j++) { 
			op[i][j] = x[i] * conj(y[j]); 
		}
	}

	return op; 
}


// Multiplies an m x n complex matrix by an n x 1 complex vector. 
complex<double> * prod (complex<double> **A, const complex<double> *x, int m, 
                        int n) {
	complex<double> *y = new complex<double>[m]; 
	for (int i = 0; i < m; i++) {
		y[i] = 0.; 
		for (int j = 0; j < n; j++) { 
			y[i] += A[i][j] * x[j]; 
		}
	}
	return y; 
}

// Multiplies an m x n complex matrix by an n x p complex matrix. 
complex<double> ** prod (complex<double> **A, complex<double> **B, int m, int n,
                         int p) { 
	complex<double> **C = new_cmat(m, p); 
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < p; j++) { 
			C[i][j] = 0.; 
			for (int k = 0; k < n; k++) { 
				C[i][j] += A[i][k] * B[k][j]; 
			}
		}
	}
	return C; 
}


// Overloaded version of prod() for a square matrix. 
complex<double> * prod (complex<double> **A, const complex<double> *x, int n) {
	return prod(A, x, n, n);
}


// Overloaded version of prod when both matrices are square n x n
complex<double> ** prod (complex<double> **A, complex<double> **B, int n) { 
	return prod(A, B, n, n, n); 
} 

// Multiplies an n x n complex matrix by an n x 1 complex vector, operating on
// the vector in place. That is, on return the vector x will be the product of
// A times the original x.  
void prod_ip (complex<double> **A, complex<double> *x, int n) {
	complex<double> *tmp = copyv(x, n); 
	for (int i = 0; i < n; i++) {
		x[i] = 0.; 
		for (int j = 0; j < n; j++) { 
			x[i] += A[i][j] * tmp[j]; 
		}
	}
	delete [] tmp; 
}

// Multiplies an m x n complex matrix A by an n x n complex matrix B, operating
// on the matrix A in place.  
void prod_ip (complex<double> **A, complex<double> **B, int m, int n) { 
	complex<double> **tmp = copym(A, m, n); 
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) { 
			A[i][j] = 0.; 
			for (int k = 0; k < n; k++) { 
				A[i][j] += tmp[i][k] * B[k][j]; 
			}
		}
	}
	deletem(tmp, m); 
}


// Overloaded version of prod_ip for when A is square n x n 
void prod_ip (complex<double> **A, complex<double> **B, int n) {
	prod_ip(A, B, n, n); 
}



// Returns the Kronecker product of two complex matrices, an m x n matrix A and
// a p x q matrix B. 
complex<double> ** kronecker (complex<double> **A, complex<double> **B, int m,
                              int n, int p, int q) {
	complex<double> **product = new_cmat(m * p, n * q); 
	int i, j; 

	for (int i1 = 0; i1 < m; i1++) {
		for (int j1 = 0; j1 < n; j1++) {
			for (int i2 = 0; i2 < p; i2++) {
				for (int j2 = 0; j2 < q; j2++) {
					i = i1 * p + i2; 
					j = j1 * q + j2; 
					product[i][j] = A[i1][j1] * B[i2][j2];
				}
			}
		}
	}
	
	return product; 
}


// Raises an n x n matrix to the power e, where e is a nonnegative integer. 
// Uses the naive method of repeated multiplication. 
complex<double> ** mpow (complex<double> **A, int e, int n) { 
	if (e < 0) throw invalid_argument("Exponent must be nonnegative"); 

	complex<double> **X = eyec(n, n); 
	
	for (int i = 1; i <= e; i++) { 
		prod_ip(X, A, n); 
	}

	return X; 
}


// Raises an n x n matrix to the power e, where e is a nonnegative integer. 
// Uses the naive method of repeated multiplication. 
// Edits the matrix A in place. 
void mpow_ip (complex<double> **A, int e, int n) { 
	if (e < 0) throw invalid_argument("Exponent must be nonnegative"); 

	complex<double> **tmp = copym(A, n); 
	for (int i = 0; i < n; i++) { 
		for (int j = 0; j < n; j++) { 
			A[i][j] = i == j; 
		}
	}
	
	for (int i = 1; i <= e; i++) { 
		prod_ip(A, tmp, n); 
	}

	deletem(tmp, n); 
}



// Tests whether an n x n complex matrix is unitary.
bool is_unitary (complex<double> **U, int n) {
	static const double EPS = 1.e-4; // tolerance for floating-point comparison
	complex<double> actual; 
	double expected; 

	// U is unitary iff UU* = I, where U* is the conjugate transpose of U. 
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			actual = inner_prod(U[j], U[i], n); // (i,j) entry of UU*
			expected = i == j ? 1. : 0.; // corresponding entry of identity matrix
			if (abs(actual - expected) > EPS) {
				return false;
			}
		}
	}

	return true; 
}


// Returns the greatest common divisor of a and b using the Euclidean algorithm.
int gcd (int a, int b) {
	if (a <= 0 || b <= 0) 
		throw invalid_argument ("arguments must be positive integers");

	int tmp;

	while (b != 0) {
		tmp = b;
		b = a % b;
		a = tmp;
	}

	return a;
}


// Tests whether n is prime.
// This uses a naive algorithm that would be useless for numbers of substantial
// size, but is fine for the purposes of QAS, which are understanding and 
// simulating the use of quantum computing rather than actually solving 
// difficult problem instances. 
bool is_prime (int n) {
	if (n <= 1) throw invalid_argument("n must be greater than 1");

	// trial division
	for (int k = 2; k <= floor(sqrt(n)); k++) { 
		if (n % k == 0) { 
			return false; 
		}
	}
	return true; 
}


// Checks whether n is the k-th power of an integer for some k > 1. If so, it 
// returns the root r such that r^k = n. Otherwise, it returns -1.
int int_root (int n) {
	if (n <= 0) throw invalid_argument("n must be positive");

	int k = 2;
	double root;
	int iroot;

	while (k <= log2(n)) {
		root = pow(n, 1. / k);
		iroot = round(root);
		if (round(pow(iroot, k)) == n) {
			return iroot;
		}
		k++;
	}

	return -1;
}


// Uniform random number generator. Returns a random number uniformly 
// distributed on [0, 1). 
double runif () {
	static bool is_seeded = false;
	static gsl_rng *rng;

	if (!is_seeded) {
		rng = init_runif ();
		is_seeded = true;
	}

	return gsl_rng_uniform(rng);
}

// Discrete uniform random number generator. Returns an integer uniformly 
// distributed between a and b, inclusive. 
int rdunif (int a, int b) {
	if (a > b) throw invalid_argument("a must be less than or equal to b");
	return (int) floor(a + runif() * (b - a + 1));
}


// Seeds the runif() function with a pseudo-random starting point, and 
// initializes a random number generator. 
// This function is automatically called by runif the first time that runif is
// called in a program, so it never needs to be called by the user. 
static gsl_rng * init_runif () {
	const gsl_rng_type *type;
	gsl_rng *rng;
	struct timeval time;

	gsl_rng_env_setup();
	type = gsl_rng_default;
	rng = gsl_rng_alloc(type);

	// Set seed: uses the "timeval" structure from sys/time.h to seed time based 
	// on microseconds rather than just seconds, as the traditional way of seeding
	// srand() does. 
	gettimeofday(&time, NULL);
	gsl_rng_set(rng, (unsigned long int) time.tv_usec * time.tv_sec); 

	// Note: rng is alloc-ed but never freed in either this function or in runif -
	// it stays in memory throughout the duration of the program. This isn't 
	// really a memory leak problem because there is only a single instance of it 
	// (since this function is only called once per program, through runif()). 

	return rng;
}


// Converts an integer to a string that represents that integer in binary form.
// len is the length of the string to return. 
// Example: int_to_bin(5, 4) = "0101". 
string int_to_bin (int n, int len) {
	if (n >= pow(2, len)) throw invalid_argument("len is too small");

	string str; 
	for (int i = 0; i < len; i++) { 
		str += get_bit(n, len - i - 1) + '0'; 
	}

	return str; 
}


// Returns the bit, zero or one, in the index-th position from the right, 
// starting at zero. That is, it returns the bit in the position with place 
// value 2^index in binary. For example:
// get_bit(1, 0) = 1
// get_bit(1, 1) = 0
// get_bit(2, 0) = 0
// get_bit(2, 1) = 1
// get_bit(5, 2) = 1 
// get_bit(7, 2) = 1 
// get_bit(7, 3) = 0 
int get_bit (int n, int index) { 
	static const int MASK = 1; 
	return (n >> index) & MASK; 
}


// Normalizes the coefficients in coeffs, a vector of length 2^n, such that 
// the squares of the absolute values of each entry sum to one. 
void normalize (complex<double> *coeffs, int n) {
	double sum = 0.; 

	for (int i = 0; i < pow(2, n); i++) {
		sum += pow(abs(coeffs[i]), 2); 
	}
	
	double factor = sqrt(sum); // normalization factor
	for (int i = 0; i < pow(2, n); i++) {
		coeffs[i] /= factor;
	}
}








