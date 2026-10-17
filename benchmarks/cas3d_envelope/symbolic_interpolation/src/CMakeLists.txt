cmake_minimum_required(VERSION 3.24)
project(axis_interpolation_proof LANGUAGES Fortran)
find_package(fortsym CONFIG REQUIRED)
add_executable(axis_interpolation_proof axis_interpolation_proof.f90)
target_link_libraries(axis_interpolation_proof PRIVATE fortsym::fortsym)
target_compile_options(axis_interpolation_proof PRIVATE
    -Wall -Wextra -Werror=array-temporaries -fcheck=all)
enable_testing()
add_test(NAME strict_axis_interpolation_identities
    COMMAND axis_interpolation_proof ${CMAKE_CURRENT_BINARY_DIR}/strict-identities.json)
