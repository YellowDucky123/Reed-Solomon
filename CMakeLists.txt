cmake_minimum_required(VERSION 3.10)
project(LibffExample CXX)

set(CMAKE_CXX_STANDARD 14)

# Define your executable
add_executable(my_project main.cpp)

# Link libff and its core dependencies directly
target_link_libraries(my_project ff gmp gmpxx ssl crypto)
