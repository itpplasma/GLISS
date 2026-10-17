cmake_minimum_required(VERSION 3.24)
project(resonant_pinch_proof LANGUAGES Fortran)
find_package(fortsym CONFIG REQUIRED)
add_executable(resonant_pinch_proof resonant_pinch_proof.f90)
target_link_libraries(resonant_pinch_proof PRIVATE fortsym::fortsym)
target_compile_options(resonant_pinch_proof PRIVATE
    -Wall -Wextra -Werror=array-temporaries -fcheck=all)
enable_testing()
add_test(NAME strict_resonant_pinch_identities
    COMMAND resonant_pinch_proof ${CMAKE_CURRENT_BINARY_DIR}/strict-identities.json)
set_tests_properties(strict_resonant_pinch_identities PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
