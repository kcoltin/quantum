project(quantum)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE debug)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_subdirectory(lib)
add_subdirectory(test)


