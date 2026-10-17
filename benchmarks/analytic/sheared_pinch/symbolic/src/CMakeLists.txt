cmake_minimum_required(VERSION 3.24)
project(sheared_pinch_proof LANGUAGES Fortran)
find_package(fortsym CONFIG REQUIRED)
add_executable(sheared_pinch_proof sheared_pinch_proof.f90)
target_link_libraries(sheared_pinch_proof PRIVATE fortsym::fortsym)
target_compile_options(sheared_pinch_proof PRIVATE
    -Wall -Wextra -Werror=array-temporaries -fcheck=all)
enable_testing()
add_test(NAME strict_sheared_pinch_identities
    COMMAND sheared_pinch_proof ${CMAKE_CURRENT_BINARY_DIR}/strict-identities.json)
set_tests_properties(strict_sheared_pinch_identities PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
