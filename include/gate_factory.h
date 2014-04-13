// gate_factory.h 

#ifndef GATE_FACTORY_H
#define GATE_FACTORY_H

#include "quantum.h"

namespace quantum_algorithm_simulator {

// Namespace containing factory methods for specific types of gates
namespace qgates {

QuantumGate x_gate (); 
QuantumGate y_gate (); 
QuantumGate z_gate (); 
QuantumGate hadamard_gate (); 
QuantumGate phase_shift_gate (double theta); 
QuantumGate swap_gate ();
QuantumGate cnot_gate (); 
QuantumGate controlled_gate (const QuantumGate &U); 
QuantumGate toffoli_gate (); 
QuantumGate fredkin_gate (); 
QuantumGate function_gate (unsigned int (*f) (unsigned int), int m, int k); 
QuantumGate grover_diffusion_operator (int n); 

}

}

#endif 
