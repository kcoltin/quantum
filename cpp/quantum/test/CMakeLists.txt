set(test_SRCS 
    ${quantum_SOURCE_DIR}/test/run_tests.cpp
    ${quantum_SOURCE_DIR}/test/test_gate_factory.cpp
    ${quantum_SOURCE_DIR}/test/test_qubitsystem.cpp
    ${quantum_SOURCE_DIR}/test/test_quantumgate.cpp
    ${quantum_SOURCE_DIR}/test/test_algorithms.cpp
   )
set(EXECUTABLE_OUTPUT_PATH ${quantum_SOURCE_DIR}/test)

include_directories(${quantum_SOURCE_DIR}/src)
include_directories(${quantum_SOURCE_DIR}/include)

add_executable(run_tests ${test_SRCS})

target_link_libraries(run_tests quantum gsl blas gtest pthread) 

