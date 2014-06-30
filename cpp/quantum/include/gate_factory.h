// gate_factory.h 

#ifndef GATE_FACTORY_H
#define GATE_FACTORY_H

#include "quantum_gate.h"
#include <functional>

namespace quantum_algorithm_simulator {

// Namespace containing factory methods for specific types of gates
namespace qgates {

QuantumGate x_gate (); 
QuantumGate y_gate (); 
QuantumGate z_gate (); 
QuantumGate hadamard_gate (); 
QuantumGate hadamard_gate (int n); 
QuantumGate phase_shift_gate (double theta); 
QuantumGate swap_gate ();
QuantumGate cnot_gate (); 
QuantumGate controlled_gate (const QuantumGate &U); 
QuantumGate toffoli_gate (); 
QuantumGate fredkin_gate (); 
QuantumGate function_gate (int (*f) (int), int m, int k); 
QuantumGate function_gate (std::function<int (int)> f, int m, int k);
QuantumGate grover_diffusion_operator (int n); 

}

}

#endif 
