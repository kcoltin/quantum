// Tests for the functionality of the QubitSystem class.

#include "algorithms.h" 
#include "quantum.h" 
#include <gtest/gtest.h>
using std::string;

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

