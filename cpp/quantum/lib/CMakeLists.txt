set(quantum_SRCS 
    ${quantum_SOURCE_DIR}/src/algorithms.cpp
    ${quantum_SOURCE_DIR}/src/gate_factory.cpp
    ${quantum_SOURCE_DIR}/src/quantum_gate.cpp
    ${quantum_SOURCE_DIR}/src/qubit_system.cpp
    ${quantum_SOURCE_DIR}/src/util.cpp
   )
set(LIBRARY_OUTPUT_PATH ${quantum_SOURCE_DIR}/lib)

include_directories(${quantum_SOURCE_DIR}/include)
add_library(quantum SHARED ${quantum_SRCS})


