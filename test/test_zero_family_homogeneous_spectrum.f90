program test_zero_family_homogeneous_spectrum
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compressible_stiffness_family_assembly, only: &
        assemble_compressible_family_stiffness
    use dynamic_family_layout, only: dynamic_family_layout_t
    use mass_density_policy, only: mass_density_profile_t
    use phase_assembly_policy, only: phase_assembly_transformed
    use physical_constants, only: vacuum_permeability
    use physical_mass_family_assembly, only: assemble_physical_family_mass
    use radial_space_policy, only: radial_space_config_t
    use symmetric_eigensolver, only: solve_symmetric_generalized
    use variable_block_tridiagonal, only: pack_variable_blocks, &
        variable_block_tridiagonal_t
    use variable_generalized_solver, only: variable_generalized_inertia, &
        variable_generalized_ok
    implicit none

    real(dp), parameter :: density_kg_m3 = 2.0_dp
    real(dp), parameter :: gamma_pressure_pa = 3.0_dp
    real(dp), parameter :: flux_t_slope = 1.0_dp
    real(dp), parameter :: current_i = -1.0_dp
    real(dp), parameter :: signed_sqrtg = -1.0_dp
    real(dp), parameter :: bmag = 1.0_dp
    integer, parameter :: meshes(4) = [8, 16, 32, 64]
    real(dp) :: errors(size(meshes))
    integer :: mesh

    call check_fixture_conventions()
    do mesh = 1, size(meshes)
        call check_homogeneous_mesh(meshes(mesh), errors(mesh))
    end do
    call require(all(errors(2:) < errors(:size(errors) - 1)), &
        "fixed-endpoint manufactured branch does not converge monotonically")
    call require(all(errors(:size(errors) - 1) / errors(2:) > 3.5_dp), &
        "fixed-endpoint manufactured branch is not second-order convergent")

    write (*, "(a)") "PASS"

contains

    subroutine check_homogeneous_mesh(intervals, continuum_error)
        integer, intent(in) :: intervals
        real(dp), intent(out) :: continuum_error
        real(dp), allocatable :: fields(:, :, :, :), drive(:, :, :)
        real(dp), allocatable :: jacobian_s(:, :, :), jacobian_t(:, :, :)
        real(dp), allocatable :: jacobian_z(:, :, :), gamma_p(:, :, :)
        real(dp), allocatable :: stiffness(:, :), mass(:, :)
        real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
        type(dynamic_family_layout_t) :: layout, mass_layout
        type(mass_density_profile_t) :: density
        type(radial_space_config_t) :: radial_space
        integer :: info

        call build_homogeneous_fields(intervals, fields, drive, jacobian_s, &
            jacobian_t, jacobian_z, gamma_p)
        call assemble_compressible_family_stiffness(fields, drive, &
            jacobian_s, jacobian_t, jacobian_z, gamma_p, [0, 0], [0, 0], &
            [1, 2], [0.0_dp, 0.0_dp], 1, radial_space, &
            1.0_dp / real(intervals, dp), phase_assembly_transformed, &
            stiffness, layout, info)
        call require(info == 0, "homogeneous N=0 stiffness assembly failed")
        call require(layout%normal_unknowns == intervals - 1 .and. &
            layout%eta_unknowns == intervals .and. &
            layout%mu_unknowns == intervals, &
            "homogeneous N=0 component dimensions are wrong")
        density%s = [0.0_dp, 1.0_dp]
        density%kilograms_per_cubic_metre = &
            [density_kg_m3, density_kg_m3]
        call assemble_physical_family_mass(fields, density, [0, 0], [0, 0], &
            [1, 2], [0.0_dp, 0.0_dp], 1, radial_space, &
            1.0_dp / real(intervals, dp), phase_assembly_transformed, mass, &
            mass_layout, info)
        call require(info == 0, "homogeneous N=0 mass assembly failed")
        call require(mass_layout%total_unknowns == layout%total_unknowns, &
            "homogeneous N=0 K/M layouts differ")
        call check_spectrum(stiffness, mass, intervals, eigenvalues, &
            eigenvectors, info)
        call certify_zero_cluster(stiffness, mass, intervals, eigenvalues)
        continuum_error = abs(eigenvalues(2 * intervals + 1) &
            - fixed_endpoint_fast_eigenvalue())
    end subroutine check_homogeneous_mesh

    subroutine check_spectrum(stiffness, mass, intervals, eigenvalues, &
            eigenvectors, info)
        real(dp), intent(in) :: stiffness(:, :), mass(:, :)
        integer, intent(in) :: intervals
        real(dp), allocatable, intent(out) :: eigenvalues(:)
        real(dp), allocatable, intent(out) :: eigenvectors(:, :)
        integer, intent(out) :: info
        real(dp) :: expected, scale
        integer :: mode, offset

        call solve_symmetric_generalized(stiffness, mass, eigenvalues, &
            eigenvectors, info)
        call require(info == 0, "homogeneous N=0 dense solve failed")
        offset = 2 * intervals
        call require(size(eigenvalues) == 3 * intervals - 1, &
            "homogeneous N=0 dimension is wrong")
        scale = discrete_fast_eigenvalue(intervals, 1)
        call require(maxval(abs(eigenvalues(:offset))) < 1.0e-12_dp * scale, &
            "homogeneous N=0 structural nullity is wrong")
        do mode = 1, intervals - 1
            expected = discrete_fast_eigenvalue(intervals, mode)
            call require(abs(eigenvalues(offset + mode) - expected) &
                < 1.0e-11_dp * expected, &
                "homogeneous N=0 discrete fast eigenvalue is wrong")
        end do
    end subroutine check_spectrum

    subroutine certify_zero_cluster(stiffness, mass, intervals, eigenvalues)
        real(dp), intent(in) :: stiffness(:, :), mass(:, :)
        integer, intent(in) :: intervals
        real(dp), intent(in) :: eigenvalues(:)
        type(variable_block_tridiagonal_t) :: k_blocks, m_blocks
        integer :: widths(1), cluster, count, info
        real(dp) :: scale, floor, first_positive, window

        cluster = 2 * intervals
        widths(1) = 3 * intervals - 1
        call pack_variable_blocks(stiffness, widths, k_blocks, info)
        call require(info == 0, "zero-family stiffness packing failed")
        call pack_variable_blocks(mass, widths, m_blocks, info)
        call require(info == 0, "zero-family mass packing failed")

        scale = discrete_fast_eigenvalue(intervals, 1)
        floor = 1.0e-7_dp * scale
        first_positive = eigenvalues(cluster + 1)
        window = 0.25_dp * (discrete_fast_eigenvalue(intervals, 2) &
            - first_positive)

        call variable_generalized_inertia(k_blocks, m_blocks, -floor, &
            count, info)
        call require(info == variable_generalized_ok .and. count == 0, &
            "N=0 spectrum has a mode below the arithmetic floor")
        call variable_generalized_inertia(k_blocks, m_blocks, floor, &
            count, info)
        call require(info == variable_generalized_ok .and. count == cluster, &
            "N=0 zero cluster is not counted at the floor")
        call variable_generalized_inertia(k_blocks, m_blocks, &
            first_positive - window, count, info)
        call require(info == variable_generalized_ok .and. count == cluster, &
            "lowest positive mode is not bracketed above the window base")
        call variable_generalized_inertia(k_blocks, m_blocks, &
            first_positive + window, count, info)
        call require(info == variable_generalized_ok &
            .and. count == cluster + 1, &
            "window does not pin the lowest positive mode")
    end subroutine certify_zero_cluster

    subroutine build_homogeneous_fields(intervals, fields, drive, &
            jacobian_s, jacobian_t, jacobian_z, gamma_p)
        integer, intent(in) :: intervals
        real(dp), allocatable, intent(out) :: fields(:, :, :, :)
        real(dp), allocatable, intent(out) :: drive(:, :, :)
        real(dp), allocatable, intent(out) :: jacobian_s(:, :, :)
        real(dp), allocatable, intent(out) :: jacobian_t(:, :, :)
        real(dp), allocatable, intent(out) :: jacobian_z(:, :, :)
        real(dp), allocatable, intent(out) :: gamma_p(:, :, :)

        allocate (fields(1, 1, 13, intervals), source=0.0_dp)
        allocate (drive(1, 1, intervals), source=0.0_dp)
        allocate (jacobian_s(1, 1, intervals), source=0.0_dp)
        allocate (jacobian_t(1, 1, intervals), source=0.0_dp)
        allocate (jacobian_z(1, 1, intervals), source=0.0_dp)
        allocate (gamma_p(1, 1, intervals), source=gamma_pressure_pa)
        fields(:, :, 1, :) = flux_t_slope
        fields(:, :, 5, :) = current_i
        fields(:, :, 7, :) = signed_sqrtg
        fields(:, :, 8, :) = bmag
        fields(:, :, 9, :) = 1.0_dp
    end subroutine build_homogeneous_fields

    subroutine check_fixture_conventions()
        real(dp) :: flux_p_slope, current_j

        flux_p_slope = 0.0_dp
        current_j = 0.0_dp
        call require(flux_t_slope * current_i + flux_p_slope * current_j &
            == bmag**2 * signed_sqrtg, &
            "homogeneous fixture violates the signed Boozer identity")
    end subroutine check_fixture_conventions

    pure function discrete_fast_eigenvalue(intervals, mode) result(value)
        integer, intent(in) :: intervals, mode
        real(dp) :: value, angle, inverse_step

        angle = acos(-1.0_dp) * real(mode, dp) / real(intervals, dp)
        inverse_step = real(intervals, dp)
        value = fast_speed_squared() * 4.0_dp * inverse_step**2 &
            * tan(0.5_dp * angle)**2
    end function discrete_fast_eigenvalue

    pure function fixed_endpoint_fast_eigenvalue() result(value)
        real(dp) :: value

        value = fast_speed_squared() * acos(-1.0_dp)**2
    end function fixed_endpoint_fast_eigenvalue

    pure function fast_speed_squared() result(value)
        real(dp) :: value

        value = (1.0_dp / vacuum_permeability + gamma_pressure_pa) &
            / density_kg_m3
    end function fast_speed_squared

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_zero_family_homogeneous_spectrum
