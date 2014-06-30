#include "algorithms.h" 
#include <gtest/gtest.h>
using std::string;
using namespace quantum_algorithm_simulator; 

// Tests Grover's search algorithm
TEST (AlgorithmsTest, TestGroverSearch) {
	// List of strings to search 
	static const string LIST[] = {"John", "Paul", "George", "Ringo", "Paul"}; 
	static const int LIST_LENGTH = 5; 
	int index; 
	
	// Test that it can correctly find each string in the list
	for (int i = 0; i < LIST_LENGTH; i++) { 
		index = grover_search(LIST[i], LIST, LIST_LENGTH); 
		EXPECT_GE(index, 0); 
		EXPECT_LT(index, LIST_LENGTH); 
		if (0 <= index && index < LIST_LENGTH) { 
			EXPECT_EQ(LIST[i], LIST[index]); 
		}
	}

	// Test that it returns -1 when the item is missing
	index = grover_search("Yoko", LIST, LIST_LENGTH); 
	EXPECT_EQ(-1, index); 
}



// Function to use to test Grover's algorithm: inverse of discrete log problem
static int f (int x) { return ((int) pow(13, x)) % 8; }

// Tests Grover's function inversion algorithm
TEST (AlgorithmsTest, TestGroverInversion) {
	static const int N = 3; // number of bits to allow as input to f
	int x; 
	
	// Test that it can correctly solve the problem for each possible y value
	for (int i = 0; i < pow(2, N); i++) { 
		x = grover_invert(f, f(i), N); 
		EXPECT_GE(x, 0); 
		if (x >= 0) {
			EXPECT_EQ(f(i), f(x)); 
		}
	}

	// Test that it returns -1 when the item is missing
	x = grover_invert(f, 3, N); 
	EXPECT_EQ(-1, x); 
}


// Test Shor's algorithm
TEST (AlgorithmsTest, TestShor) {
	int composites[] = {4, 6, 9, 10, 12}; 
	const int LENGTH = 5; 

	for (int i = 0; i < LENGTH; i++) {
		int n = composites[i]; 
		int p = shor_factor(n); 
		// test that p is a factor
		EXPECT_EQ(0, n % p); 
		// test that p is not trivial
		EXPECT_NE(1, p); 
		EXPECT_NE(n, p); 
	}

	// Test that it returns 1 for prime numbers
	int primes[] = {2, 3, 7, 29}; 
	const int PRIMES_LENGTH = 4; 
	for (int i = 0; i < PRIMES_LENGTH; i++) { 
		int p = primes[i]; 
		EXPECT_EQ(1, shor_factor(p)); 
	}
}






