cmake_minimum_required(VERSION 3.22)

project(gliss VERSION 0.1.0 LANGUAGES Fortran)

include(CTest)
find_package(LAPACK REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF_C REQUIRED IMPORTED_TARGET netcdf)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-Wall -Wextra -Wno-compare-reals)
endif()

add_library(gliss
    src/gliss_capi.f90
    src/gvec_cas3d_adapter.f90
    src/netcdf_c_bindings.f90
    src/netcdf_c_api.f90
    src/gvec_cas3d_netcdf.f90
    src/gvec_cas3d_reader.f90
    src/gvec_cas3d_reconstruction.f90
    src/gvec_cas3d_integrals.f90
    src/gvec_cas3d_types.f90
    src/block_tridiagonal.f90
    src/compressible_stiffness_assembly.f90
    src/compressible_stiffness_family_assembly.f90
    src/generalized_block_solver.f90
    src/eigenvalue_tracking.f90
    src/dynamic_family_layout.f90
    src/family_point_assembly.f90
    src/local_mode_model.f90
    src/mass_density_policy.f90
    src/phase_assembly_policy.f90
    src/physical_mass_assembly.f90
    src/physical_mass_family_assembly.f90
    src/physical_mass_kernel.f90
    src/nonuniform_derivative.f90
    src/fourier_phase_kind.f90
    src/mode_topology.f90
    src/trial_space_topology.f90
    src/phase_factor_topology.f90
    src/physical_constants.f90
    src/radial_space_policy.f90
    src/symmetric_pivot_inertia.f90
    src/terpsichore_topology.f90
    src/three_component_kernel.f90
    src/two_component_kernel.f90
    src/variable_block_tridiagonal.f90
    src/variable_generalized_solver.f90
    src/newcomb_limit.f90
    src/family_assembly.f90
    src/field_profile_identities.f90
    src/mercier_diagnostic.f90
    src/radial_bspline_basis.f90
    src/symmetric_eigensolver.f90)
target_link_libraries(gliss PUBLIC
    LAPACK::LAPACK
    PkgConfig::NETCDF_C)
set_target_properties(gliss PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_include_directories(gliss PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>")

add_library(gliss_c SHARED src/gliss_capi.f90)

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

add_executable(gliss_newcomb app/gliss_newcomb.f90)
target_link_libraries(gliss_newcomb PRIVATE gliss)
set_target_properties(gliss_newcomb PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_mercier app/gliss_mercier.f90)
target_link_libraries(gliss_mercier PRIVATE gliss)
set_target_properties(gliss_mercier PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_field_profiles app/gliss_field_profiles.f90)
target_link_libraries(gliss_field_profiles PRIVATE gliss)
set_target_properties(gliss_field_profiles PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_family app/gliss_family.f90)
target_link_libraries(gliss_family PRIVATE gliss)
set_target_properties(gliss_family PROPERTIES
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

    add_executable(test_radial_space_policy
        test/test_radial_space_policy.f90)
    target_link_libraries(test_radial_space_policy PRIVATE gliss)
    add_test(NAME radial_space_policy COMMAND test_radial_space_policy)

    add_executable(test_gvec_cas3d_reconstruction
        test/test_gvec_cas3d_reconstruction.f90)
    target_link_libraries(test_gvec_cas3d_reconstruction PRIVATE gliss)
    add_test(NAME gvec_cas3d_reconstruction
        COMMAND test_gvec_cas3d_reconstruction)

    add_library(gliss_test_support test/cylinder_fixture.f90)
    target_link_libraries(gliss_test_support PUBLIC gliss)

    add_executable(test_export_assembly test/test_export_assembly.f90)
    target_link_libraries(test_export_assembly PRIVATE gliss_test_support)
    add_test(NAME export_assembly COMMAND test_export_assembly)

    add_executable(test_family_assembly test/test_family_assembly.f90)
    target_link_libraries(test_family_assembly PRIVATE gliss)
    add_test(NAME family_assembly COMMAND test_family_assembly)

    add_executable(test_generalized_block_solver
        test/test_generalized_block_solver.f90)
    target_link_libraries(test_generalized_block_solver PRIVATE gliss)
    add_test(NAME generalized_block_solver COMMAND test_generalized_block_solver)

    add_executable(test_newcomb_limit test/test_newcomb_limit.f90)
    target_link_libraries(test_newcomb_limit PRIVATE gliss)
    add_test(NAME newcomb_limit COMMAND test_newcomb_limit)

    add_executable(test_two_component_kernel
        test/test_two_component_kernel.f90)
    target_link_libraries(test_two_component_kernel PRIVATE gliss)
    add_test(NAME two_component_kernel COMMAND test_two_component_kernel)

    add_executable(test_three_component_kernel
        test/test_three_component_kernel.f90)
    target_link_libraries(test_three_component_kernel PRIVATE gliss)
    add_test(NAME three_component_kernel COMMAND test_three_component_kernel)

    add_executable(test_physical_mass_kernel
        test/test_physical_mass_kernel.f90)
    target_link_libraries(test_physical_mass_kernel PRIVATE gliss)
    add_test(NAME physical_mass_kernel COMMAND test_physical_mass_kernel)

    add_executable(test_mass_density_policy
        test/test_mass_density_policy.f90)
    target_link_libraries(test_mass_density_policy PRIVATE gliss)
    add_test(NAME mass_density_policy COMMAND test_mass_density_policy)

    add_executable(test_physical_mass_assembly
        test/test_physical_mass_assembly.f90)
    target_link_libraries(test_physical_mass_assembly PRIVATE gliss)
    add_test(NAME physical_mass_assembly COMMAND test_physical_mass_assembly)

    add_executable(test_compressible_stiffness_assembly
        test/test_compressible_stiffness_assembly.f90)
    target_link_libraries(test_compressible_stiffness_assembly PRIVATE gliss)
    add_test(NAME compressible_stiffness_assembly
        COMMAND test_compressible_stiffness_assembly)

    add_executable(test_compressible_stiffness_family_assembly
        test/test_compressible_stiffness_family_assembly.f90)
    target_link_libraries(test_compressible_stiffness_family_assembly PRIVATE
        gliss)
    add_test(NAME compressible_stiffness_family_assembly
        COMMAND test_compressible_stiffness_family_assembly)

    add_executable(test_dynamic_family_layout
        test/test_dynamic_family_layout.f90)
    target_link_libraries(test_dynamic_family_layout PRIVATE gliss)
    add_test(NAME dynamic_family_layout COMMAND test_dynamic_family_layout)

    add_executable(test_variable_block_tridiagonal
        test/test_variable_block_tridiagonal.f90)
    target_link_libraries(test_variable_block_tridiagonal PRIVATE gliss)
    add_test(NAME variable_block_tridiagonal
        COMMAND test_variable_block_tridiagonal)

    add_executable(test_variable_block_factorization
        test/test_variable_block_factorization.f90)
    target_link_libraries(test_variable_block_factorization PRIVATE gliss)
    add_test(NAME variable_block_factorization
        COMMAND test_variable_block_factorization)

    add_executable(test_variable_generalized_solver
        test/test_variable_generalized_solver.f90)
    target_link_libraries(test_variable_generalized_solver PRIVATE gliss)
    add_test(NAME variable_generalized_solver
        COMMAND test_variable_generalized_solver)

    add_executable(test_physical_mass_family_assembly
        test/test_physical_mass_family_assembly.f90)
    target_link_libraries(test_physical_mass_family_assembly PRIVATE gliss)
    add_test(NAME physical_mass_family_assembly
        COMMAND test_physical_mass_family_assembly)

    add_executable(test_zero_family_homogeneous_spectrum
        test/test_zero_family_homogeneous_spectrum.f90)
    target_link_libraries(test_zero_family_homogeneous_spectrum PRIVATE gliss)
    add_test(NAME zero_family_homogeneous_spectrum
        COMMAND test_zero_family_homogeneous_spectrum)

    add_executable(test_mode_topology test/test_mode_topology.f90)
    target_link_libraries(test_mode_topology PRIVATE gliss)
    add_test(NAME mode_topology COMMAND test_mode_topology)

    add_executable(test_trial_space_topology
        test/test_trial_space_topology.f90)
    target_link_libraries(test_trial_space_topology PRIVATE gliss)
    add_test(NAME trial_space_topology COMMAND test_trial_space_topology)

    add_executable(test_phase_factor_topology
        test/test_phase_factor_topology.f90)
    target_link_libraries(test_phase_factor_topology PRIVATE gliss)
    add_test(NAME phase_factor_topology COMMAND test_phase_factor_topology)

    add_executable(test_terpsichore_topology
        test/test_terpsichore_topology.f90)
    target_link_libraries(test_terpsichore_topology PRIVATE gliss)
    add_test(NAME terpsichore_topology COMMAND test_terpsichore_topology)

    add_executable(test_nonuniform_derivative
        test/test_nonuniform_derivative.f90)
    target_link_libraries(test_nonuniform_derivative PRIVATE gliss)
    add_test(NAME nonuniform_derivative COMMAND test_nonuniform_derivative)

    add_executable(test_mercier_diagnostic test/test_mercier_diagnostic.f90)
    target_link_libraries(test_mercier_diagnostic PRIVATE gliss_test_support)
    add_test(NAME mercier_diagnostic COMMAND test_mercier_diagnostic)

    add_executable(test_field_profile_identities
        test/test_field_profile_identities.f90)
    target_link_libraries(test_field_profile_identities PRIVATE
        gliss_test_support)
    add_test(NAME field_profile_identities
        COMMAND test_field_profile_identities)

    add_executable(test_gvec_cas3d_integrals test/test_gvec_cas3d_integrals.f90)
    target_link_libraries(test_gvec_cas3d_integrals PRIVATE gliss)
    add_test(NAME gvec_cas3d_integrals COMMAND test_gvec_cas3d_integrals)

    add_executable(test_suydam_marginal test/test_suydam_marginal.f90)
    target_link_libraries(test_suydam_marginal PRIVATE gliss_test_support)
    add_test(NAME suydam_marginal_m4 COMMAND test_suydam_marginal 4)
    add_test(NAME suydam_marginal_m6 COMMAND test_suydam_marginal 6)

    add_executable(test_eigenvalue_tracking
        test/test_eigenvalue_tracking.f90)
    target_link_libraries(test_eigenvalue_tracking PRIVATE gliss)
    add_test(NAME eigenvalue_tracking COMMAND test_eigenvalue_tracking)
endif()

option(GLISS_ENABLE_ENZYME "Build Enzyme derivative checks" OFF)
if(GLISS_ENABLE_ENZYME)
    include(cmake/EnzymeFortran.cmake)
    add_enzyme_fortran_test(enzyme_drive_gradient
        SOURCES
            src/physical_constants.f90
            src/local_mode_model.f90
            test/enzyme_drive_gradient.f90)
    add_enzyme_fortran_test(enzyme_phase_assembly_gradient
        SOURCES
            src/fourier_phase_kind.f90
            src/phase_factor_topology.f90
            src/radial_space_policy.f90
            src/two_component_kernel.f90
            src/family_point_assembly.f90
            test/enzyme_phase_fixture.f90
            test/enzyme_phase_assembly_gradient.f90)
    add_enzyme_fortran_test(enzyme_physical_mass_gradient
        SOURCES
            src/physical_mass_kernel.f90
            test/enzyme_physical_mass_gradient.f90)
    add_enzyme_fortran_test(enzyme_physical_mass_assembly_gradient
        SOURCES
            src/fourier_phase_kind.f90
            src/mass_density_policy.f90
            src/phase_assembly_policy.f90
            src/phase_factor_topology.f90
            src/physical_mass_kernel.f90
            src/radial_space_policy.f90
            src/physical_mass_assembly.f90
            test/enzyme_physical_mass_assembly_fixture.f90
            test/enzyme_physical_mass_assembly_gradient.f90)
    add_enzyme_fortran_test(enzyme_physical_mass_family_gradient
        SOURCES
            src/fourier_phase_kind.f90
            src/trial_space_topology.f90
            src/dynamic_family_layout.f90
            src/mass_density_policy.f90
            src/phase_assembly_policy.f90
            src/phase_factor_topology.f90
            src/physical_mass_kernel.f90
            src/radial_space_policy.f90
            src/physical_mass_assembly.f90
            src/physical_mass_family_assembly.f90
            test/enzyme_physical_mass_family_fixture.f90
            test/enzyme_physical_mass_family_gradient.f90)
    add_enzyme_fortran_test(enzyme_three_component_gradient
        SOURCES
            src/physical_constants.f90
            src/two_component_kernel.f90
            src/three_component_kernel.f90
            test/enzyme_three_component_gradient.f90)
    add_enzyme_fortran_test(enzyme_compressible_stiffness_gradient
        SOURCES
            src/fourier_phase_kind.f90
            src/phase_assembly_policy.f90
            src/phase_factor_topology.f90
            src/physical_constants.f90
            src/radial_space_policy.f90
            src/two_component_kernel.f90
            src/three_component_kernel.f90
            src/compressible_stiffness_assembly.f90
            test/enzyme_compressible_stiffness_fixture.f90
            test/enzyme_compressible_stiffness_gradient.f90)
    add_enzyme_fortran_test(enzyme_compressible_stiffness_family_gradient
        SOURCES
            src/fourier_phase_kind.f90
            src/trial_space_topology.f90
            src/dynamic_family_layout.f90
            src/phase_assembly_policy.f90
            src/phase_factor_topology.f90
            src/physical_constants.f90
            src/radial_space_policy.f90
            src/two_component_kernel.f90
            src/three_component_kernel.f90
            src/compressible_stiffness_assembly.f90
            src/compressible_stiffness_family_assembly.f90
            test/enzyme_compressible_stiffness_family_fixture.f90
            test/enzyme_compressible_stiffness_family_gradient.f90)
endif()
