cmake_minimum_required(VERSION 3.24)
project(weighted_derivative_proof LANGUAGES Fortran)
find_package(fortsym CONFIG REQUIRED)
add_executable(weighted_derivative_proof weighted_derivative_proof.f90)
target_link_libraries(weighted_derivative_proof PRIVATE fortsym::fortsym)
target_compile_options(weighted_derivative_proof PRIVATE
    -Wall -Wextra -Werror=array-temporaries -fcheck=all)
enable_testing()
add_test(NAME strict_weighted_derivative_identities
    COMMAND weighted_derivative_proof ${CMAKE_CURRENT_BINARY_DIR}/strict-identities.json)
set_tests_properties(strict_weighted_derivative_identities PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
