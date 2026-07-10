cmake_minimum_required(VERSION 3.22)

project(gliss VERSION 0.1.0 LANGUAGES Fortran)

include(CTest)
find_package(LAPACK REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF_C REQUIRED IMPORTED_TARGET netcdf)

add_library(gliss
    src/gvec_cas3d_adapter.f90
    src/netcdf_c_bindings.f90
    src/netcdf_c_api.f90
    src/gvec_cas3d_netcdf.f90
    src/gvec_cas3d_reader.f90
    src/gvec_cas3d_reconstruction.f90
    src/gvec_cas3d_integrals.f90
    src/gvec_cas3d_types.f90
    src/local_mode_model.f90
    src/radial_bspline_basis.f90
    src/symmetric_eigensolver.f90)
target_link_libraries(gliss PUBLIC
    LAPACK::LAPACK
    PkgConfig::NETCDF_C)
set_target_properties(gliss PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_include_directories(gliss PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>")

add_executable(gliss_demo app/gliss_demo.f90)
target_link_libraries(gliss_demo PRIVATE gliss)
set_target_properties(gliss_demo PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gvec_cas3d_inspect app/gvec_cas3d_inspect.f90)
target_link_libraries(gvec_cas3d_inspect PRIVATE gliss)
set_target_properties(gvec_cas3d_inspect PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gvec_cas3d_compare_pointwise
    app/gvec_cas3d_compare_pointwise.f90)
target_link_libraries(gvec_cas3d_compare_pointwise PRIVATE gliss)
set_target_properties(gvec_cas3d_compare_pointwise PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gvec_cas3d_integrate app/gvec_cas3d_integrate.f90)
target_link_libraries(gvec_cas3d_integrate PRIVATE gliss)
set_target_properties(gvec_cas3d_integrate PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

if(BUILD_TESTING)
    add_executable(test_local_mode test/test_local_mode.f90)
    target_link_libraries(test_local_mode PRIVATE gliss)
    add_test(NAME local_mode COMMAND test_local_mode)

    add_executable(test_eigen_sensitivity test/test_eigen_sensitivity.f90)
    target_link_libraries(test_eigen_sensitivity PRIVATE gliss)
    add_test(NAME eigen_sensitivity COMMAND test_eigen_sensitivity)

    add_executable(test_gvec_cas3d_adapter
        test/test_gvec_cas3d_adapter.f90)
    target_link_libraries(test_gvec_cas3d_adapter PRIVATE gliss)
    add_test(NAME gvec_cas3d_adapter COMMAND test_gvec_cas3d_adapter)

    add_executable(test_gvec_cas3d_reader
        test/test_gvec_cas3d_reader.f90)
    target_link_libraries(test_gvec_cas3d_reader PRIVATE gliss)
    add_test(NAME gvec_cas3d_reader COMMAND test_gvec_cas3d_reader)

    add_executable(test_radial_bspline_basis
        test/test_radial_bspline_basis.f90)
    target_link_libraries(test_radial_bspline_basis PRIVATE gliss)
    add_test(NAME radial_bspline_basis COMMAND test_radial_bspline_basis)

    add_executable(test_gvec_cas3d_reconstruction
        test/test_gvec_cas3d_reconstruction.f90)
    target_link_libraries(test_gvec_cas3d_reconstruction PRIVATE gliss)
    add_test(NAME gvec_cas3d_reconstruction
        COMMAND test_gvec_cas3d_reconstruction)

    add_executable(test_gvec_cas3d_integrals test/test_gvec_cas3d_integrals.f90)
    target_link_libraries(test_gvec_cas3d_integrals PRIVATE gliss)
    add_test(NAME gvec_cas3d_integrals COMMAND test_gvec_cas3d_integrals)
endif()

option(GLISS_ENABLE_ENZYME "Build Enzyme derivative checks" OFF)
if(GLISS_ENABLE_ENZYME)
    include(cmake/EnzymeFortran.cmake)
    add_enzyme_fortran_test(enzyme_drive_gradient
        SOURCES
            src/local_mode_model.f90
            test/enzyme_drive_gradient.f90)
endif()
