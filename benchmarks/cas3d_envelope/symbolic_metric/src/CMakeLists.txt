cmake_minimum_required(VERSION 3.24)
project(metric_averaging_proof LANGUAGES Fortran)
find_package(fortsym CONFIG REQUIRED)
add_executable(metric_averaging_proof metric_averaging_proof.f90)
target_link_libraries(metric_averaging_proof PRIVATE fortsym::fortsym)
target_compile_options(metric_averaging_proof PRIVATE
    -Wall -Wextra -Werror=array-temporaries -fcheck=all)
enable_testing()
add_test(NAME strict_metric_averaging_identity
    COMMAND metric_averaging_proof ${CMAKE_CURRENT_BINARY_DIR}/strict-identities.json)
