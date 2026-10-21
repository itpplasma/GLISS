cmake_minimum_required(VERSION 3.22)
project(ideal_mhd_batches LANGUAGES Fortran)
find_package(OpenMP REQUIRED COMPONENTS Fortran)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})
set(shared_root ${CMAKE_CURRENT_SOURCE_DIR}/../..)
add_library(mhd_batch STATIC ${shared_root}/src/generated/ideal_mhd_batches.f90
    ${shared_root}/src/generated/ideal_mhd_energy_generated.f90
    ${shared_root}/src/physical_constants.f90
    ${shared_root}/src/two_component_kernel.f90
    ${shared_root}/src/three_component_kernel.f90
    ${shared_root}/src/perpendicular_kinetic_kernel.f90)
target_link_libraries(mhd_batch PUBLIC OpenMP::OpenMP_Fortran)
target_include_directories(mhd_batch PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
add_executable(benchmark_mhd_batches ${shared_root}/test/test_ideal_mhd_batches.f90)
target_link_libraries(benchmark_mhd_batches PRIVATE mhd_batch)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(mhd_batch PRIVATE -O3 -fno-fast-math -ffp-contract=off
        -finline-limit=100000)
endif()
enable_testing()
add_test(NAME mhd_batch_oracles COMMAND benchmark_mhd_batches)
option(MHD_BATCH_CUDA "Build the CUDA batch oracle and timing harness" OFF)
if(MHD_BATCH_CUDA)
    if(NOT CMAKE_CUDA_ARCHITECTURES)
        message(FATAL_ERROR "Set CMAKE_CUDA_ARCHITECTURES explicitly")
    endif()
    enable_language(CUDA)
    add_executable(benchmark_mhd_batches_cuda gpu.cu reference.f90)
    target_link_libraries(benchmark_mhd_batches_cuda PRIVATE mhd_batch)
    target_include_directories(benchmark_mhd_batches_cuda PRIVATE ${shared_root}/src/generated)
    set_target_properties(benchmark_mhd_batches_cuda PROPERTIES CUDA_STANDARD 17)
    target_compile_options(benchmark_mhd_batches_cuda PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:-O3>)
    add_test(NAME mhd_batch_cuda_oracles COMMAND benchmark_mhd_batches_cuda 257 3)
endif()
