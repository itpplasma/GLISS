cmake_minimum_required(VERSION 3.22)

project(gliss VERSION 0.0.1 LANGUAGES C Fortran)

include(CTest)
set(BLA_VENDOR "Generic" CACHE STRING "BLAS/LAPACK implementation")
find_package(LAPACK REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF_C REQUIRED IMPORTED_TARGET netcdf)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-Wall>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wextra>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-compare-reals>)
    if(SKBUILD)
        add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-ffile-prefix-map=${CMAKE_CURRENT_SOURCE_DIR}=.>)
    endif()
endif()

add_library(gliss
    src/gvec_cas3d_adapter.f90
    src/netcdf_c_bindings.f90
    src/netcdf_c_api.f90
    src/gvec_cas3d_netcdf.f90
    src/gvec_cas3d_reader.f90
    src/gvec_cas3d_writer.f90
    src/gvec_cas3d_reconstruction.f90
    src/gvec_cas3d_integrals.f90
    src/gvec_cas3d_types.f90
    src/block_tridiagonal.f90
    src/compressible_geometry.f90
    src/compressible_stiffness_validation.f90
    src/compatible_compressible_stiffness_assembly.f90
    src/compressible_stiffness_assembly.f90
    src/compressible_stiffness_family_assembly.f90
    src/generalized_block_solver.f90
    src/eigenvalue_tracking.f90
    src/dynamic_family_layout.f90
    src/family_point_assembly.f90
    src/helical_cylinder_limit.f90
    src/local_mode_model.f90
    src/mass_density_policy.f90
    src/phase_assembly_policy.f90
    src/compatible_physical_mass_assembly.f90
    src/physical_mass_assembly.f90
    src/physical_mass_family_assembly.f90
    src/physical_mass_kernel.f90
    src/nonuniform_derivative.f90
    src/fourier_phase_kind.f90
    src/mode_topology.f90
    src/trial_space_topology.f90
    src/phase_factor_topology.f90
    src/physical_constants.f90
    src/radial_cubic_spline.f90
    src/axis_regular_harmonic_spline.f90
    src/cartesian_harmonic_spline.f90
    src/compatible_family_point_assembly.f90
    src/compatible_radial_quadrature.f90
    src/compatible_operator_trace_types.f90
    src/compatible_problem_assembly_support.f90
    src/field_periodic_cartesian.f90
    src/cartesian_primitive_geometry.f90
    src/primitive_geometry_grid.f90
    src/primitive_equilibrium_spline.f90
    src/radial_feec_complex.f90
    src/radial_space_policy.f90
    src/perpendicular_kinetic_kernel.f90
    src/symmetric_pivot_inertia.f90
    src/terpsichore_model_policy.f90
    src/terpsichore_normalization.f90
    src/terpsichore_pair_average.f90
    src/terpsichore_matrix_fixture.f90
    src/terpsichore_noninteracting_coefficients.f90
    src/terpsichore_noninteracting_stiffness.f90
    src/terpsichore_pseudoplasma_fixture.f90
    src/terpsichore_pseudoplasma_stiffness.f90
    src/terpsichore_pseudoplasma_coupling.f90
    src/terpsichore_eigen_diagnostics.f90
    src/terpsichore_reduced_layout.f90
    src/terpsichore_reduced_mass.f90
    src/terpsichore_reduced_mass_adapter.f90
    src/terpsichore_reduced_mass_family_assembly.f90
    src/terpsichore_solution_fixture.f90
    src/terpsichore_topology.f90
    src/terpsichore_fixed_boundary_spectrum.f90
    src/three_component_kernel.f90
    src/two_component_kernel.f90
    src/boundary/starwall_cylinder_limit.f90
    src/boundary/starwall_mesh_geometry.f90
    src/boundary/starwall_ideal_vacuum.f90
    src/boundary/starwall_fourier_coupling.f90
    src/vacuum_schur.f90
    src/variable_block_tridiagonal.f90
    src/variable_generalized_equilibration.f90
    src/stable_reduction.f90
    src/fixed_boundary_solver_controls.f90
    src/variable_generalized_solver.f90
    src/variable_spectrum_analysis.f90
    src/dense_spectrum_support.f90
    src/dense_generalized_inverse_iteration.f90
    src/fixed_boundary_energy.f90
    src/fixed_boundary_eigen_bracket.f90
    src/fixed_boundary_input_validation.f90
    src/fixed_boundary_spectrum.f90
    src/newcomb_limit.f90
    src/family_assembly.f90
    src/two_component_artificial_problem.f90
    src/two_component_spectrum.f90
    src/axisymmetric_spectrum.f90
    src/field_profile_identities.f90
    src/export_surface_geometry.f90
    src/primitive_kernel_geometry.f90
    src/compatible_two_component_problem.f90
    src/compatible_three_component_problem.f90
    src/mercier_diagnostic.f90
    src/radial_bspline_basis.f90
    src/symmetric_eigensolver.f90)
target_link_libraries(gliss PUBLIC
    LAPACK::LAPACK
    PkgConfig::NETCDF_C)
set_target_properties(gliss PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
    POSITION_INDEPENDENT_CODE ON)
target_include_directories(gliss PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>")

add_library(gliss_c SHARED
    src/gliss_c_abi_support.f90
    src/gliss_c_contexts.f90
    src/gliss_equilibrium_io_capi.f90
    src/gliss_capi.f90
    src/gliss_spectrum_capi.f90
    src/gliss_solver_controls_capi.f90
    src/gliss_energy_capi.f90
    src/gliss_terpsichore_capi.f90
    src/gliss_axisymmetric_capi.f90
    src/gliss_marginality_capi.f90
    src/gliss_full_spectrum_capi.f90)
add_dependencies(gliss_c gliss)
target_link_libraries(gliss_c PRIVATE gliss)
install(TARGETS gliss_c
    LIBRARY DESTINATION gliss
    RUNTIME DESTINATION gliss)
if(SKBUILD)
    install(FILES include/gliss.h DESTINATION gliss/include)
else()
    install(FILES include/gliss.h DESTINATION include)
endif()

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

add_executable(gliss_local_mode_diagnostic
    app/gliss_local_mode_diagnostic.f90)
target_link_libraries(gliss_local_mode_diagnostic PRIVATE gliss)
set_target_properties(gliss_local_mode_diagnostic PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_vacuum_schur_diagnostic
    app/gliss_vacuum_schur_diagnostic.f90)
target_link_libraries(gliss_vacuum_schur_diagnostic PRIVATE gliss)
set_target_properties(gliss_vacuum_schur_diagnostic PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_starwall_diagnostic
    app/gliss_starwall_diagnostic.f90)
target_link_libraries(gliss_starwall_diagnostic PRIVATE gliss)
set_target_properties(gliss_starwall_diagnostic PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_mercier app/gliss_mercier.f90)
target_link_libraries(gliss_mercier PRIVATE gliss)
set_target_properties(gliss_mercier PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_field_profiles app/gliss_field_profiles.f90)
target_link_libraries(gliss_field_profiles PRIVATE gliss)
set_target_properties(gliss_field_profiles PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_primitive_geometry_trace
    app/gliss_primitive_geometry_trace.f90)
target_link_libraries(gliss_primitive_geometry_trace PRIVATE gliss)
set_target_properties(gliss_primitive_geometry_trace PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_family app/gliss_family.f90)
target_link_libraries(gliss_family PRIVATE gliss)
set_target_properties(gliss_family PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_spectrum app/gliss_spectrum.f90)
target_link_libraries(gliss_spectrum PRIVATE gliss)
set_target_properties(gliss_spectrum PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_compatible_spectrum app/gliss_compatible_spectrum.f90)
target_link_libraries(gliss_compatible_spectrum PRIVATE gliss)
set_target_properties(gliss_compatible_spectrum PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_compatible_marginality
    app/gliss_compatible_marginality.f90)
target_link_libraries(gliss_compatible_marginality PRIVATE gliss)
set_target_properties(gliss_compatible_marginality PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_compatible_operator_trace
    app/compatible_operator_geometry_trace.f90
    app/gliss_compatible_operator_trace.f90)
target_link_libraries(gliss_compatible_operator_trace PRIVATE gliss)
set_target_properties(gliss_compatible_operator_trace PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_axisymmetric app/gliss_axisymmetric.f90)
target_link_libraries(gliss_axisymmetric PRIVATE gliss)
set_target_properties(gliss_axisymmetric PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_axisymmetric_fkg app/gliss_axisymmetric_fkg.f90)
target_link_libraries(gliss_axisymmetric_fkg PRIVATE gliss)
set_target_properties(gliss_axisymmetric_fkg PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_terpsichore_eigen
    app/gliss_terpsichore_eigen.f90)
target_link_libraries(gliss_terpsichore_eigen PRIVATE gliss)
set_target_properties(gliss_terpsichore_eigen PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gvec_cas3d_integrate app/gvec_cas3d_integrate.f90)
target_link_libraries(gvec_cas3d_integrate PRIVATE gliss)
set_target_properties(gvec_cas3d_integrate PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_terpsichore_reduced_mass
    app/gliss_terpsichore_reduced_mass.f90)
target_link_libraries(gliss_terpsichore_reduced_mass PRIVATE gliss)
set_target_properties(gliss_terpsichore_reduced_mass PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_terpsichore_potential_inspect
    app/gliss_terpsichore_potential_inspect.f90)
target_link_libraries(gliss_terpsichore_potential_inspect PRIVATE gliss)
set_target_properties(gliss_terpsichore_potential_inspect PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_terpsichore_noninteracting_stiffness
    app/gliss_terpsichore_noninteracting_stiffness.f90)
target_link_libraries(gliss_terpsichore_noninteracting_stiffness PRIVATE gliss)
set_target_properties(gliss_terpsichore_noninteracting_stiffness PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

add_executable(gliss_terpsichore_pseudoplasma_stiffness
    app/gliss_terpsichore_pseudoplasma_stiffness.f90)
target_link_libraries(gliss_terpsichore_pseudoplasma_stiffness PRIVATE gliss)
set_target_properties(gliss_terpsichore_pseudoplasma_stiffness PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fo/bin")

if(BUILD_TESTING)
    function(add_spectrum_cli_rejection name expected)
        add_test(NAME ${name} COMMAND ${CMAKE_COMMAND}
            -DPROGRAM=$<TARGET_FILE:gliss_spectrum>
            "-DEXPECTED=gliss_spectrum: ${expected}"
            "-DARGUMENTS=${ARGN}"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)
    endfunction()

    add_spectrum_cli_rejection(spectrum_cli_bad_gamma
        "GAMMA must be a number" missing.nc bad 1.0 1.0 1,1)
    add_spectrum_cli_rejection(spectrum_cli_negative_gamma
        "GAMMA must be nonnegative" missing.nc -1.0 1.0 1.0 1,1)
    add_spectrum_cli_rejection(spectrum_cli_nonfinite_floor
        "FLOOR must be finite" missing.nc 1.6666666666666667 1.0 NaN 1,1)
    add_spectrum_cli_rejection(spectrum_cli_zero_density
        "DENSITY must be positive" missing.nc 1.6666666666666667 0.0 1.0 1,1)
    add_spectrum_cli_rejection(spectrum_cli_negative_m
        "poloidal mode m must be nonnegative"
        missing.nc 1.6666666666666667 1.0 1.0 -1,1)
    add_spectrum_cli_rejection(spectrum_cli_duplicate_mode
        "duplicate mode" missing.nc 1.6666666666666667 1.0 1.0 1,1 1,1)
    add_spectrum_cli_rejection(spectrum_cli_incomplete_family
        "--family requires exactly INDEX MMAX NMAX"
        missing.nc 1.6666666666666667 1.0 1.0 --family 0 8)
    add_spectrum_cli_rejection(spectrum_cli_bad_quadrature
        "quadrature RULE must be midpoint"
        missing.nc 1.6666666666666667 1.0 1.0 --quadrature bad 1,1)
    add_spectrum_cli_rejection(spectrum_cli_rejects_gauss2
        "quadrature RULE must be midpoint"
        missing.nc 1.6666666666666667 1.0 1.0 --quadrature gauss2 1,1)
    add_spectrum_cli_rejection(spectrum_cli_incomplete_quadrature
        "--quadrature requires RULE and a mode selector"
        missing.nc 1.6666666666666667 1.0 1.0 --quadrature)

    function(add_axisymmetric_cli_rejection name expected)
        add_test(NAME ${name} COMMAND ${CMAKE_COMMAND}
            -DPROGRAM=$<TARGET_FILE:gliss_axisymmetric>
            "-DEXPECTED=gliss_axisymmetric: ${expected}"
            "-DARGUMENTS=${ARGN}"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)
    endfunction()

    add_axisymmetric_cli_rejection(axisymmetric_cli_bad_n
        "N must be an integer" missing.nc bad 4)
    add_axisymmetric_cli_rejection(axisymmetric_cli_zero_n
        "N must be positive" missing.nc 0 4)
    add_axisymmetric_cli_rejection(axisymmetric_cli_zero_mmax
        "MMAX must be positive" missing.nc 1 0)
    add_axisymmetric_cli_rejection(axisymmetric_cli_bad_option
        "the optional argument must be --quadrature"
        missing.nc 1 4 --bad gauss2)
    add_axisymmetric_cli_rejection(axisymmetric_cli_bad_quadrature
        "quadrature RULE must be midpoint"
        missing.nc 1 4 --quadrature bad)
    add_axisymmetric_cli_rejection(axisymmetric_cli_rejects_gauss2
        "quadrature RULE must be midpoint"
        missing.nc 1 4 --quadrature gauss2)

    function(add_axisymmetric_fkg_cli_rejection name expected)
        add_test(NAME ${name} COMMAND ${CMAKE_COMMAND}
            -DPROGRAM=$<TARGET_FILE:gliss_axisymmetric_fkg>
            "-DEXPECTED=gliss_axisymmetric_fkg: ${expected}"
            "-DARGUMENTS=${ARGN}"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)
    endfunction()

    add_axisymmetric_fkg_cli_rejection(axisymmetric_fkg_cli_bad_n
        "N must be an integer" missing.nc bad 4)
    add_axisymmetric_fkg_cli_rejection(axisymmetric_fkg_cli_zero_mmax
        "MMAX must be positive" missing.nc 1 0)
    add_axisymmetric_cli_rejection(axisymmetric_cli_bad_final_option
        "the final option must be --count-only"
        missing.nc 1 4 --quadrature midpoint --bad)

    add_executable(test_local_mode test/test_local_mode.f90)
    target_link_libraries(test_local_mode PRIVATE gliss)
    add_test(NAME local_mode COMMAND test_local_mode)
    add_test(NAME local_mode_diagnostic COMMAND gliss_local_mode_diagnostic)

    add_executable(test_starwall_ideal_vacuum
        test/boundary/test_starwall_ideal_vacuum.f90)
    target_link_libraries(test_starwall_ideal_vacuum PRIVATE gliss)
    add_test(NAME starwall_ideal_vacuum COMMAND test_starwall_ideal_vacuum)
    add_executable(test_starwall_fourier_coupling
        test/boundary/test_starwall_fourier_coupling.f90)
    target_link_libraries(test_starwall_fourier_coupling PRIVATE gliss)
    add_test(NAME starwall_fourier_coupling
        COMMAND test_starwall_fourier_coupling)
    add_test(NAME starwall_diagnostic_cli_extra COMMAND ${CMAKE_COMMAND}
        -DPROGRAM=$<TARGET_FILE:gliss_starwall_diagnostic>
        "-DEXPECTED=gliss_starwall_diagnostic: no arguments are accepted"
        "-DARGUMENTS=extra"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)
    add_test(NAME local_mode_diagnostic_cli_extra COMMAND ${CMAKE_COMMAND}
        -DPROGRAM=$<TARGET_FILE:gliss_local_mode_diagnostic>
        "-DEXPECTED=gliss_local_mode_diagnostic: no arguments are accepted"
        "-DARGUMENTS=extra"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)

    add_executable(test_vacuum_schur test/test_vacuum_schur.f90)
    target_link_libraries(test_vacuum_schur PRIVATE gliss)
    add_test(NAME vacuum_schur COMMAND test_vacuum_schur)
    add_test(NAME vacuum_schur_diagnostic COMMAND gliss_vacuum_schur_diagnostic)
    add_test(NAME vacuum_schur_diagnostic_cli_extra COMMAND ${CMAKE_COMMAND}
        -DPROGRAM=$<TARGET_FILE:gliss_vacuum_schur_diagnostic>
        "-DEXPECTED=gliss_vacuum_schur_diagnostic: no arguments are accepted"
        "-DARGUMENTS=extra"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)

    add_executable(test_helical_cylinder_limit
        test/test_helical_cylinder_limit.f90)
    target_link_libraries(test_helical_cylinder_limit PRIVATE gliss)
    add_test(NAME helical_cylinder_limit COMMAND test_helical_cylinder_limit)

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

    add_executable(test_gvec_cas3d_writer
        test/test_gvec_cas3d_writer.f90)
    target_link_libraries(test_gvec_cas3d_writer PRIVATE gliss_test_support)
    add_test(NAME gvec_cas3d_writer COMMAND test_gvec_cas3d_writer)

    add_executable(test_radial_bspline_basis
        test/test_radial_bspline_basis.f90)
    target_link_libraries(test_radial_bspline_basis PRIVATE gliss)
    add_test(NAME radial_bspline_basis COMMAND test_radial_bspline_basis)

    add_executable(test_radial_feec_complex
        test/test_radial_feec_complex.f90)
    target_link_libraries(test_radial_feec_complex PRIVATE gliss)
    add_test(NAME radial_feec_complex COMMAND test_radial_feec_complex)

    add_executable(test_compatible_radial_quadrature
        test/test_compatible_radial_quadrature.f90)
    target_link_libraries(test_compatible_radial_quadrature PRIVATE gliss)
    add_test(NAME compatible_radial_quadrature
        COMMAND test_compatible_radial_quadrature)

    add_executable(test_compatible_problem_assembly_support
        test/test_compatible_problem_assembly_support.f90)
    target_link_libraries(test_compatible_problem_assembly_support PRIVATE gliss)
    add_test(NAME compatible_problem_assembly_support
        COMMAND test_compatible_problem_assembly_support)

    function(add_compatible_cli_rejection name target expected)
        add_test(NAME ${name} COMMAND ${CMAKE_COMMAND}
            -DPROGRAM=$<TARGET_FILE:${target}>
            "-DEXPECTED=${expected}"
            "-DARGUMENTS=${ARGN}"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)
    endfunction()

    add_compatible_cli_rejection(compatible_spectrum_cli_missing
        gliss_compatible_spectrum
        "gliss_compatible_spectrum: missing required arguments")
    add_compatible_cli_rejection(compatible_spectrum_cli_nonfinite_gamma
        gliss_compatible_spectrum
        "gliss_compatible_spectrum: GAMMA must be finite"
        missing.nc 2 8 8 nan 1 1 1 1,0)
    add_compatible_cli_rejection(compatible_spectrum_cli_duplicate_mode
        gliss_compatible_spectrum
        "gliss_compatible_spectrum: duplicate mode"
        missing.nc 2 8 8 1 1 1 1 1,0 1,0)
    add_compatible_cli_rejection(compatible_marginality_cli_missing
        gliss_compatible_marginality
        "gliss_compatible_marginality: missing required arguments")
    add_compatible_cli_rejection(compatible_marginality_cli_nonfinite_floor
        gliss_compatible_marginality
        "gliss_compatible_marginality: FLOOR must be finite"
        missing.nc 2 8 8 nan 1 1,0)
    add_compatible_cli_rejection(compatible_marginality_cli_missing_density
        gliss_compatible_marginality
        "gliss_compatible_marginality: --physical-density requires a value"
        missing.nc 2 8 8 1e-8 1 --physical-density=)
    add_compatible_cli_rejection(compatible_marginality_cli_nonfinite_density
        gliss_compatible_marginality
        "gliss_compatible_marginality: physical density must be finite"
        missing.nc 2 8 8 1e-8 1 --physical-density=nan 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_duplicate_density
        gliss_compatible_marginality
        "gliss_compatible_marginality: duplicate --physical-density"
        missing.nc 2 8 8 1e-8 1 --physical-density=1
        --physical-density=1 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_nonfinite_m0_power
        gliss_compatible_marginality
        "gliss_compatible_marginality: m=0 stored power must be finite"
        missing.nc 2 8 8 1e-8 1 --m0-stored-power=nan 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_duplicate_option
        gliss_compatible_marginality
        "gliss_compatible_marginality: duplicate --m0-stored-power"
        missing.nc 2 8 8 1e-8 1 --m0-stored-power=0.5
        --m0-stored-power=0.5 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_bad_power_count
        gliss_compatible_marginality
        "gliss_compatible_marginality: stored powers requires one value per mode"
        missing.nc 2 8 8 1e-8 1 --stored-powers=0,0 0,1 1,1 2,1)
    add_compatible_cli_rejection(compatible_marginality_cli_nonfinite_power
        gliss_compatible_marginality
        "gliss_compatible_marginality: stored powers must be finite"
        missing.nc 2 8 8 1e-8 1 --stored-powers=0,nan,0 0,1 1,1 2,1)
    add_compatible_cli_rejection(compatible_marginality_cli_power_conflict
        gliss_compatible_marginality
        "gliss_compatible_marginality: --m0-stored-power conflicts"
        missing.nc 2 8 8 1e-8 1 --m0-stored-power=0
        --stored-powers=0,0,0 0,1 1,1 2,1)
    add_compatible_cli_rejection(compatible_marginality_cli_bad_eta_power_count
        gliss_compatible_marginality
        "gliss_compatible_marginality: eta stored powers requires one value per mode"
        missing.nc 2 8 8 1e-8 1 --eta-stored-powers=0,0 0,1 1,1 2,1)
    add_compatible_cli_rejection(compatible_marginality_cli_nonfinite_eta_power
        gliss_compatible_marginality
        "gliss_compatible_marginality: eta stored powers must be finite"
        missing.nc 2 8 8 1e-8 1 --eta-stored-powers=0,nan,0 0,1 1,1 2,1)
    add_compatible_cli_rejection(compatible_marginality_cli_duplicate_eta_power
        gliss_compatible_marginality
        "gliss_compatible_marginality: duplicate --eta-stored-powers"
        missing.nc 2 8 8 1e-8 1 --eta-stored-powers=0,0,0
        --eta-stored-powers=0,0,0 0,1 1,1 2,1)
    add_compatible_cli_rejection(compatible_marginality_cli_bad_bracket_order
        gliss_compatible_marginality
        "gliss_compatible_marginality: negative bracket must obey"
        missing.nc 2 8 8 1e-8 1 --negative-bracket=1e-6,1e-8,1e-3 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_bad_bracket_shape
        gliss_compatible_marginality
        "gliss_compatible_marginality: negative bracket requires"
        missing.nc 2 8 8 1e-8 1 --negative-bracket=1e-8,1e-6 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_bad_eigen_bracket_order
        gliss_compatible_marginality
        "gliss_compatible_marginality: eigenvalue bracket must obey"
        missing.nc 2 8 8 1e-8 1 --eigenvalue-bracket=1,1,1e-3 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_bad_eigen_bracket_shape
        gliss_compatible_marginality
        "gliss_compatible_marginality: eigenvalue bracket requires"
        missing.nc 2 8 8 1e-8 1 --eigenvalue-bracket=0,1 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_nonfinite_eigen_bracket
        gliss_compatible_marginality
        "gliss_compatible_marginality: eigenvalue bracket lower must be finite"
        missing.nc 2 8 8 1e-8 1 --eigenvalue-bracket=nan,1,1e-3 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_duplicate_eigen_bracket
        gliss_compatible_marginality
        "gliss_compatible_marginality: duplicate --eigenvalue-bracket"
        missing.nc 2 8 8 1e-8 1 --eigenvalue-bracket=0,1,1e-3
        --eigenvalue-bracket=0,1,1e-3 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_conflicting_brackets
        gliss_compatible_marginality
        "gliss_compatible_marginality: conflicting eigenvalue bracket options"
        missing.nc 2 8 8 1e-8 1 --negative-bracket=1e-8,1e-6,1e-3
        --eigenvalue-bracket=0,1,1e-3 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_missing_count_shifts
        gliss_compatible_marginality
        "gliss_compatible_marginality: --count-shifts requires"
        missing.nc 2 8 8 1e-8 1 --count-shifts=)
    add_compatible_cli_rejection(compatible_marginality_cli_duplicate_count_shifts
        gliss_compatible_marginality
        "gliss_compatible_marginality: duplicate --count-shifts"
        missing.nc 2 8 8 1e-8 1 --count-shifts=0,1
        --count-shifts=0,1 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_empty_count_shift
        gliss_compatible_marginality
        "gliss_compatible_marginality: count shifts require nonempty"
        missing.nc 2 8 8 1e-8 1 --count-shifts=0,,1 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_nonfinite_count_shift
        gliss_compatible_marginality
        "gliss_compatible_marginality: count shift must be finite"
        missing.nc 2 8 8 1e-8 1 --count-shifts=0,nan 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_negative_count_shift
        gliss_compatible_marginality
        "gliss_compatible_marginality: count shifts must be nonnegative"
        missing.nc 2 8 8 1e-8 1 --count-shifts=-1,0 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_unsorted_count_shifts
        gliss_compatible_marginality
        "gliss_compatible_marginality: count shifts must be strictly increasing"
        missing.nc 2 8 8 1e-8 1 --count-shifts=0,1,1 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_count_bracket_conflict
        gliss_compatible_marginality
        "gliss_compatible_marginality: conflicting eigenvalue count options"
        missing.nc 2 8 8 1e-8 1 --eigenvalue-bracket=0,1,1e-3
        --count-shifts=0,1 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_count_inertia_conflict
        gliss_compatible_marginality
        "gliss_compatible_marginality: --count-shifts conflicts with --inertia-only"
        missing.nc 2 8 8 1e-8 1 --count-shifts=0,1 --inertia-only 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_missing_eigenprofile_index
        gliss_compatible_marginality
        "gliss_compatible_marginality: --eigenprofile-index requires"
        missing.nc 2 8 8 1e-8 1 --eigenprofile-index= 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_bad_eigenprofile_index
        gliss_compatible_marginality
        "gliss_compatible_marginality: eigenprofile index must be positive"
        missing.nc 2 8 8 1e-8 1 --eigenprofile-index=0 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_noninteger_eigenprofile_index
        gliss_compatible_marginality
        "gliss_compatible_marginality: eigenprofile index must be an integer"
        missing.nc 2 8 8 1e-8 1 --eigenprofile-index=one 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_duplicate_eigenprofile
        gliss_compatible_marginality
        "gliss_compatible_marginality: duplicate --eigenprofile-index"
        missing.nc 2 8 8 1e-8 1 --eigenprofile-index=1
        --eigenprofile-index=2 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_profile_points_without_profile
        gliss_compatible_marginality
        "gliss_compatible_marginality: --profile-points requires --eigenprofile-index"
        missing.nc 2 8 8 1e-8 1 --profile-points=8 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_bad_profile_points
        gliss_compatible_marginality
        "gliss_compatible_marginality: profile points must lie in [8,10000]"
        missing.nc 2 8 8 1e-8 1 --eigenprofile-index=1
        --profile-points=7 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_missing_profile_points
        gliss_compatible_marginality
        "gliss_compatible_marginality: --profile-points requires an integer"
        missing.nc 2 8 8 1e-8 1 --eigenprofile-index=1
        --profile-points= 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_noninteger_profile_points
        gliss_compatible_marginality
        "gliss_compatible_marginality: profile points must be an integer"
        missing.nc 2 8 8 1e-8 1 --eigenprofile-index=1
        --profile-points=many 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_duplicate_profile_points
        gliss_compatible_marginality
        "gliss_compatible_marginality: duplicate --profile-points"
        missing.nc 2 8 8 1e-8 1 --eigenprofile-index=1
        --profile-points=8 --profile-points=9 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_profile_count_conflict
        gliss_compatible_marginality
        "gliss_compatible_marginality: --eigenprofile-index conflicts with --count-shifts"
        missing.nc 2 8 8 1e-8 1 --eigenprofile-index=1
        --count-shifts=0,1 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_profile_missing_bracket
        gliss_compatible_marginality
        "gliss_compatible_marginality: --eigenprofile-index requires --eigenvalue-bracket"
        missing.nc 2 8 8 1e-8 1 --eigenprofile-index=1 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_profile_negative_conflict
        gliss_compatible_marginality
        "gliss_compatible_marginality: --eigenprofile-index conflicts with --negative-bracket"
        missing.nc 2 8 8 1e-8 1 --eigenprofile-index=1
        --negative-bracket=1e-3,1e-2,1e-3 0,1)
    add_compatible_cli_rejection(compatible_marginality_cli_profile_inertia_conflict
        gliss_compatible_marginality
        "gliss_compatible_marginality: --eigenprofile-index conflicts with --inertia-only"
        missing.nc 2 8 8 1e-8 1 --eigenprofile-index=1 --inertia-only 0,1)
    add_compatible_cli_rejection(compatible_trace_cli_missing
        gliss_compatible_operator_trace
        "gliss_compatible_operator_trace: missing required arguments")
    add_compatible_cli_rejection(compatible_trace_cli_bad_theta
        gliss_compatible_operator_trace
        "gliss_compatible_operator_trace: NTHETA and NZETA must be at least 8"
        missing.nc 2 7 8 1 1,0)
    add_compatible_cli_rejection(compatible_trace_cli_bad_mode
        gliss_compatible_operator_trace
        "gliss_compatible_operator_trace: modes must be given once as m,n"
        missing.nc 2 8 8 1 broken)
    add_compatible_cli_rejection(compatible_trace_cli_missing_first_cells
        gliss_compatible_operator_trace
        "gliss_compatible_operator_trace: --first-cells requires a value"
        missing.nc 2 8 8 1 --first-cells=)
    add_compatible_cli_rejection(compatible_trace_cli_bad_first_cells
        gliss_compatible_operator_trace
        "gliss_compatible_operator_trace: --first-cells must be between"
        missing.nc 2 8 8 1 --first-cells=0 0,1)
    add_compatible_cli_rejection(compatible_trace_cli_unknown_option
        gliss_compatible_operator_trace
        "gliss_compatible_operator_trace: unknown option"
        missing.nc 2 8 8 1 --cells=1 0,1)

    add_executable(test_radial_cubic_spline
        test/test_radial_cubic_spline.f90)
    target_link_libraries(test_radial_cubic_spline PRIVATE gliss)
    add_test(NAME radial_cubic_spline COMMAND test_radial_cubic_spline)

    add_executable(test_axis_regular_harmonic_spline
        test/test_axis_regular_harmonic_spline.f90)
    target_link_libraries(test_axis_regular_harmonic_spline PRIVATE gliss)
    add_test(NAME axis_regular_harmonic_spline
        COMMAND test_axis_regular_harmonic_spline)

    add_executable(test_cartesian_primitive_geometry
        test/test_cartesian_primitive_geometry.f90)
    target_link_libraries(test_cartesian_primitive_geometry PRIVATE gliss)
    add_test(NAME cartesian_primitive_geometry
        COMMAND test_cartesian_primitive_geometry)

    add_executable(test_cartesian_harmonic_spline
        test/test_cartesian_harmonic_spline.f90)
    target_link_libraries(test_cartesian_harmonic_spline PRIVATE gliss)
    add_test(NAME cartesian_harmonic_spline
        COMMAND test_cartesian_harmonic_spline)

    add_executable(test_field_periodic_cartesian
        test/test_field_periodic_cartesian.f90)
    target_link_libraries(test_field_periodic_cartesian PRIVATE gliss)
    add_test(NAME field_periodic_cartesian
        COMMAND test_field_periodic_cartesian)

    add_test(NAME primitive_geometry_trace_usage
        COMMAND ${CMAKE_COMMAND}
            -DCOMMAND_PATH=$<TARGET_FILE:gliss_primitive_geometry_trace>
            "-DEXPECTED=usage: gliss_primitive_geometry_trace"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_failure.cmake)

    add_test(NAME primitive_geometry_trace_invalid_theta
        COMMAND ${CMAKE_COMMAND}
            -DCOMMAND_PATH=$<TARGET_FILE:gliss_primitive_geometry_trace>
            "-DEXPECTED=NTHETA must be an integer of at least 8"
            -DARG1=missing.nc -DARG2=bad -DARG3=8
            -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_failure.cmake)

    add_test(NAME primitive_geometry_trace_invalid_zeta
        COMMAND ${CMAKE_COMMAND}
            -DCOMMAND_PATH=$<TARGET_FILE:gliss_primitive_geometry_trace>
            "-DEXPECTED=NZETA must be an integer of at least 8"
            -DARG1=missing.nc -DARG2=8 -DARG3=7
            -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_failure.cmake)

    add_test(NAME primitive_geometry_trace_missing_export
        COMMAND ${CMAKE_COMMAND}
            -DCOMMAND_PATH=$<TARGET_FILE:gliss_primitive_geometry_trace>
            "-DEXPECTED=export could not be read"
            -DARG1=missing.nc -DARG2=8 -DARG3=8
            -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_failure.cmake)

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

    add_executable(test_gliss_capi test/test_gliss_capi.f90)
    target_link_libraries(test_gliss_capi PRIVATE gliss_c gliss_test_support)
    add_test(NAME gliss_capi COMMAND test_gliss_capi)

    add_executable(test_gliss_spectrum_capi
        test/test_gliss_spectrum_capi.f90)
    target_link_libraries(test_gliss_spectrum_capi PRIVATE
        gliss_c gliss_test_support)
    add_test(NAME gliss_spectrum_capi COMMAND test_gliss_spectrum_capi)

    add_executable(test_gliss_marginality_capi
        test/test_gliss_marginality_capi.f90)
    target_link_libraries(test_gliss_marginality_capi PRIVATE
        gliss_c gliss_test_support)
    add_test(NAME gliss_marginality_capi COMMAND test_gliss_marginality_capi)

    add_executable(test_gliss_header cmake/test_gliss_header.c)
    target_include_directories(test_gliss_header PRIVATE include)
    target_link_libraries(test_gliss_header PRIVATE gliss_c)
    add_test(NAME gliss_header COMMAND test_gliss_header)

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

    add_executable(test_dense_generalized_inverse_iteration
        test/test_dense_generalized_inverse_iteration.f90)
    target_link_libraries(test_dense_generalized_inverse_iteration PRIVATE gliss)
    add_test(NAME dense_generalized_inverse_iteration
        COMMAND test_dense_generalized_inverse_iteration)

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

    add_executable(test_perpendicular_kinetic_kernel
        test/test_perpendicular_kinetic_kernel.f90)
    target_link_libraries(test_perpendicular_kinetic_kernel PRIVATE gliss)
    add_test(NAME perpendicular_kinetic_kernel
        COMMAND test_perpendicular_kinetic_kernel)

    add_executable(test_mass_density_policy
        test/test_mass_density_policy.f90)
    target_link_libraries(test_mass_density_policy PRIVATE gliss)
    add_test(NAME mass_density_policy COMMAND test_mass_density_policy)

    add_executable(test_physical_mass_assembly
        test/test_physical_mass_assembly.f90)
    target_link_libraries(test_physical_mass_assembly PRIVATE gliss)
    add_test(NAME physical_mass_assembly COMMAND test_physical_mass_assembly)

    add_executable(test_compatible_vector_evaluation
        test/test_compatible_vector_evaluation.f90)
    target_link_libraries(test_compatible_vector_evaluation PRIVATE gliss)
    add_test(NAME compatible_vector_evaluation
        COMMAND test_compatible_vector_evaluation)

    add_executable(test_compressible_stiffness_assembly
        test/test_compressible_stiffness_assembly.f90)
    target_link_libraries(test_compressible_stiffness_assembly PRIVATE gliss)
    add_test(NAME compressible_stiffness_assembly
        COMMAND test_compressible_stiffness_assembly)

    add_executable(test_compatible_three_component_assembly
        test/test_compatible_three_component_assembly.f90)
    target_link_libraries(test_compatible_three_component_assembly PRIVATE gliss)
    add_test(NAME compatible_three_component_assembly
        COMMAND test_compatible_three_component_assembly)

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

    add_executable(test_variable_spectrum_analysis
        test/test_variable_spectrum_analysis.f90)
    target_link_libraries(test_variable_spectrum_analysis PRIVATE gliss)
    add_test(NAME variable_spectrum_analysis
        COMMAND test_variable_spectrum_analysis)

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

    add_executable(test_cylinder_physical_spectrum
        test/test_cylinder_physical_spectrum.f90)
    target_link_libraries(test_cylinder_physical_spectrum PRIVATE
        gliss_test_support)
    add_test(NAME cylinder_physical_spectrum
        COMMAND test_cylinder_physical_spectrum)

    add_executable(test_fixed_boundary_spectrum
        test/test_fixed_boundary_spectrum.f90)
    target_link_libraries(test_fixed_boundary_spectrum PRIVATE
        gliss_test_support)
    add_test(NAME fixed_boundary_spectrum
        COMMAND test_fixed_boundary_spectrum)

    add_executable(test_two_component_spectrum
        test/test_two_component_spectrum.f90)
    target_link_libraries(test_two_component_spectrum PRIVATE
        gliss_test_support)
    add_test(NAME two_component_spectrum
        COMMAND test_two_component_spectrum)

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

    add_executable(test_terpsichore_model_policy
        test/test_terpsichore_model_policy.f90)
    target_link_libraries(test_terpsichore_model_policy PRIVATE gliss)
    add_test(NAME terpsichore_model_policy COMMAND test_terpsichore_model_policy)

    add_executable(test_terpsichore_reduced_mass
        test/test_terpsichore_reduced_mass.f90)
    target_link_libraries(test_terpsichore_reduced_mass PRIVATE gliss)
    add_test(NAME terpsichore_reduced_mass
        COMMAND test_terpsichore_reduced_mass)

    add_executable(test_terpsichore_reduced_mass_adapter
        test/test_terpsichore_reduced_mass_adapter.f90)
    target_link_libraries(test_terpsichore_reduced_mass_adapter PRIVATE gliss)
    add_test(NAME terpsichore_reduced_mass_adapter
        COMMAND test_terpsichore_reduced_mass_adapter)

    add_executable(test_terpsichore_noninteracting_stiffness
        test/test_terpsichore_noninteracting_stiffness.f90)
    target_link_libraries(test_terpsichore_noninteracting_stiffness PRIVATE gliss)
    add_test(NAME terpsichore_noninteracting_stiffness
        COMMAND test_terpsichore_noninteracting_stiffness)

    add_executable(test_terpsichore_pseudoplasma_stiffness
        test/test_terpsichore_pseudoplasma_stiffness.f90)
    target_link_libraries(test_terpsichore_pseudoplasma_stiffness PRIVATE gliss)
    add_test(NAME terpsichore_pseudoplasma_stiffness
        COMMAND test_terpsichore_pseudoplasma_stiffness)
    add_test(NAME terpsichore_pseudoplasma_stiffness_cli_missing
        COMMAND ${CMAKE_COMMAND}
        -DPROGRAM=$<TARGET_FILE:gliss_terpsichore_pseudoplasma_stiffness>
        "-DEXPECTED=usage: gliss_terpsichore_pseudoplasma_stiffness FORT.24"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)
    add_test(NAME terpsichore_pseudoplasma_stiffness_cli_extra
        COMMAND ${CMAKE_COMMAND}
        -DPROGRAM=$<TARGET_FILE:gliss_terpsichore_pseudoplasma_stiffness>
        "-DEXPECTED=usage: gliss_terpsichore_pseudoplasma_stiffness FORT.24"
        "-DARGUMENTS=one;two"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)

    add_executable(test_terpsichore_solution_fixture
        test/test_terpsichore_solution_fixture.f90)
    target_link_libraries(test_terpsichore_solution_fixture PRIVATE gliss)
    add_test(NAME terpsichore_solution_fixture
        COMMAND test_terpsichore_solution_fixture)

    add_executable(test_terpsichore_eigen_diagnostics
        test/test_terpsichore_eigen_diagnostics.f90)
    target_link_libraries(test_terpsichore_eigen_diagnostics PRIVATE gliss)
    add_test(NAME terpsichore_eigen_diagnostics
        COMMAND test_terpsichore_eigen_diagnostics)
    add_test(NAME terpsichore_eigen_cli_missing
        COMMAND ${CMAKE_COMMAND}
        -DPROGRAM=$<TARGET_FILE:gliss_terpsichore_eigen>
        "-DEXPECTED=usage: gliss_terpsichore_eigen FORT.23 IVAC [FORT.24]"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)
    add_test(NAME terpsichore_eigen_cli_bad_ivac
        COMMAND ${CMAKE_COMMAND}
        -DPROGRAM=$<TARGET_FILE:gliss_terpsichore_eigen>
        "-DEXPECTED=IVAC must be a nonnegative decimal integer"
        "-DARGUMENTS=unused;1x"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)
    add_test(NAME terpsichore_eigen_cli_negative_ivac
        COMMAND ${CMAKE_COMMAND}
        -DPROGRAM=$<TARGET_FILE:gliss_terpsichore_eigen>
        "-DEXPECTED=IVAC must be a nonnegative decimal integer"
        "-DARGUMENTS=unused;-1"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)
    add_test(NAME terpsichore_eigen_cli_missing_vacuum
        COMMAND ${CMAKE_COMMAND}
        -DPROGRAM=$<TARGET_FILE:gliss_terpsichore_eigen>
        "-DEXPECTED=IVAC>0 requires FORT.24"
        "-DARGUMENTS=unused;1"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)
    add_test(NAME terpsichore_eigen_cli_extra_fixed
        COMMAND ${CMAKE_COMMAND}
        -DPROGRAM=$<TARGET_FILE:gliss_terpsichore_eigen>
        "-DEXPECTED=IVAC=0 does not accept FORT.24"
        "-DARGUMENTS=unused;0;unused"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)

    add_executable(test_terpsichore_normalization
        test/test_terpsichore_normalization.f90)
    target_link_libraries(test_terpsichore_normalization PRIVATE gliss)
    add_test(NAME terpsichore_normalization
        COMMAND test_terpsichore_normalization)

    add_executable(test_terpsichore_reduced_mass_family_assembly
        test/test_terpsichore_reduced_mass_family_assembly.f90)
    target_link_libraries(test_terpsichore_reduced_mass_family_assembly PRIVATE gliss)
    add_test(NAME terpsichore_reduced_mass_family_assembly
        COMMAND test_terpsichore_reduced_mass_family_assembly)

    add_executable(test_nonuniform_derivative
        test/test_nonuniform_derivative.f90)
    target_link_libraries(test_nonuniform_derivative PRIVATE gliss)
    add_test(NAME nonuniform_derivative COMMAND test_nonuniform_derivative)

    add_executable(test_mercier_diagnostic test/test_mercier_diagnostic.f90)
    target_link_libraries(test_mercier_diagnostic PRIVATE gliss_test_support)
    add_test(NAME mercier_diagnostic COMMAND test_mercier_diagnostic)

    add_executable(test_mercier_gradient test/test_mercier_gradient.f90)
    target_link_libraries(test_mercier_gradient PRIVATE gliss)
    add_test(NAME mercier_gradient COMMAND test_mercier_gradient)

    add_executable(test_mercier_full_gradient
        test/test_mercier_full_gradient.f90)
    target_link_libraries(test_mercier_full_gradient PRIVATE gliss)
    add_test(NAME mercier_full_gradient COMMAND test_mercier_full_gradient)

    add_executable(test_mercier_fluxslope_gradient
        test/test_mercier_fluxslope_gradient.f90)
    target_link_libraries(test_mercier_fluxslope_gradient PRIVATE gliss)
    add_test(NAME mercier_fluxslope_gradient
        COMMAND test_mercier_fluxslope_gradient)

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
    add_test(NAME suydam_marginal_cli_bad_mode COMMAND ${CMAKE_COMMAND}
        -DPROGRAM=$<TARGET_FILE:test_suydam_marginal>
        "-DEXPECTED=test_suydam_marginal: MODE must be 4 or 6"
        "-DARGUMENTS=5"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)
    add_test(NAME suydam_marginal_cli_extra COMMAND ${CMAKE_COMMAND}
        -DPROGRAM=$<TARGET_FILE:test_suydam_marginal>
        "-DEXPECTED=test_suydam_marginal: usage:"
        "-DARGUMENTS=4;out.csv;extra"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_cli_rejection.cmake)

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
    add_enzyme_fortran_test(enzyme_helical_cylinder_gradient
        SOURCES
            src/helical_cylinder_limit.f90
            test/enzyme_helical_cylinder_gradient.f90)
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
    add_enzyme_fortran_test(enzyme_perpendicular_kinetic_gradient
        SOURCES
            src/perpendicular_kinetic_kernel.f90
            test/enzyme_perpendicular_kinetic_gradient.f90)
    add_enzyme_fortran_test(enzyme_terpsichore_reduced_mass_gradient
        SOURCES
            src/terpsichore_reduced_mass.f90
            test/enzyme_terpsichore_reduced_mass_gradient.f90)
    add_enzyme_fortran_test(enzyme_terpsichore_normalization_gradient
        SOURCES
            src/terpsichore_normalization.f90
            src/terpsichore_reduced_mass.f90
            test/enzyme_terpsichore_normalization_fixture.f90
            test/enzyme_terpsichore_normalization_gradient.f90)
    add_enzyme_fortran_test(enzyme_terpsichore_reduced_mass_family_gradient
        SOURCES
            src/fourier_phase_kind.f90
            src/trial_space_topology.f90
            src/dynamic_family_layout.f90
            src/terpsichore_reduced_layout.f90
            src/terpsichore_reduced_mass.f90
            src/terpsichore_reduced_mass_family_assembly.f90
            test/enzyme_terpsichore_reduced_mass_family_fixture.f90
            test/enzyme_terpsichore_reduced_mass_family_gradient.f90)
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
