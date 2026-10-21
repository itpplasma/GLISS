cmake_minimum_required(VERSION 3.24)
project(free_boundary_proof LANGUAGES Fortran)
find_package(fortsym CONFIG REQUIRED)
add_executable(prove_free_boundary prove_free_boundary.f90)
target_link_libraries(prove_free_boundary PRIVATE fortsym::fortsym)
enable_testing()
add_test(NAME strict_free_boundary_identities
    COMMAND prove_free_boundary ${CMAKE_CURRENT_BINARY_DIR}/identities.json)
set_tests_properties(strict_free_boundary_identities PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
