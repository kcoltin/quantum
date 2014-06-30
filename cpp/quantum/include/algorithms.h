#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include "qubit_system.h"
#include <string>

namespace quantum_algorithm_simulator {

void qft (QubitSystem *q); 
int grover_search (const std::string &match_text, const std::string *list, 
                   int n); 
int grover_invert (int (*f) (int), int y, int n); 
int shor_factor (int n); 

}

#endif 


