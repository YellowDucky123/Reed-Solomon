cmake_minimum_required(VERSION 3.14)
project(ReedSolomon CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# End-to-end demo executable (header-only library, so just the one source).
add_executable(rs_demo demo.cpp)

# Unit tests via GoogleTest.
enable_testing()
find_package(GTest REQUIRED)

add_executable(rs_tests
	tests/field_test.cpp
	tests/poly_test.cpp
	tests/reed_solomon_test.cpp
)
target_link_libraries(rs_tests GTest::gtest GTest::gtest_main pthread)

add_test(NAME rs_tests COMMAND rs_tests)
