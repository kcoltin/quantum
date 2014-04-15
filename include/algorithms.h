#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <string>

namespace quantum_algorithm_simulator {

int grover_search (const std::string &match_text, const std::string *list, 
                   int n); 
int grover_invert (int (*f) (int), int y, int n); 
int shor_factor (int n); 

}

#endif 


