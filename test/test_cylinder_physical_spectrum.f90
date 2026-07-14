program test_cylinder_physical_spectrum
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use compressible_geometry, only: build_compressible_geometry, &
        compressible_geometry_invalid_input, compressible_geometry_ok
    use compressible_stiffness_family_assembly, only: &
        assemble_compressible_family_stiffness
    use cylinder_fixture, only: create_cylinder_fixture
    use dynamic_family_layout, only: build_dynamic_block_permutation, &
        dynamic_family_layout_t, dynamic_layout_ok
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mass_density_policy, only: mass_density_profile_t
    use mercier_diagnostic, only: build_kernel_geometry, mercier_ok
    use phase_assembly_policy, only: phase_assembly_transformed
    use physical_mass_family_assembly, only: assemble_physical_family_mass
    use radial_space_policy, only: radial_space_config_t
    use symmetric_eigensolver, only: solve_symmetric_generalized
    use variable_block_tridiagonal, only: pack_permuted_variable_blocks, &
        variable_block_ok, variable_block_tridiagonal_t
    use variable_generalized_solver, only: &
        iterate_variable_generalized_eigenvalue, variable_generalized_ok, &
        variable_generalized_inertia
    implicit none

    real(dp), parameter :: density_kg_m3 = 2.0_dp
    real(dp), parameter :: adiabatic_index = 5.0_dp / 3.0_dp
    integer, parameter :: n_angle = 16
    integer, parameter :: meshes(3) = [32, 64, 128]
    ! frozen references from derivations/cylinder_compressional_spectrum.wl
    ! (theta pinch mode (3,1); Suydam member (4,4) at pressure scale 1;
    ! non-resonant member (2,1) at pressure scale 1/4; rho = 2 kg/m^3,
    ! gamma = 5/3)
    real(dp), parameter :: theta_slow_point = 9.257320410606333_dp
    real(dp), parameter :: theta_alfven_point = 44209.70641441537_dp
    real(dp), parameter :: theta_fast_lowest = 2.814092222219040e7_dp
    real(dp), parameter :: theta_fast_second = 1.023131542357869e8_dp
    real(dp), parameter :: screw_stable_slow_edge = 17.44662897569245_dp
    real(dp), parameter :: screw_unstable_omega2 = -180.4232327098574_dp

    real(dp) :: slow_errors(size(meshes)), fast_errors(size(meshes))
    real(dp) :: second_errors(size(meshes)), growth_errors(size(meshes))
    real(dp) :: slow_values(size(meshes)), fast_values(size(meshes))
    real(dp) :: second_values(size(meshes))
    real(dp) :: growth_values(size(meshes)), growth_extrapolated
    character(len=1024) :: output_path
    integer :: mesh
    logical :: output_requested

    call validate_arguments()
    call check_bridge_rejections()
    do mesh = 1, size(meshes)
        call check_theta_pinch(meshes(mesh), slow_errors(mesh), &
            fast_errors(mesh), second_errors(mesh), slow_values(mesh), &
            fast_values(mesh), second_values(mesh))
        call check_screw_stable(meshes(mesh))
        call check_screw_unstable(meshes(mesh), growth_values(mesh), &
            mesh == size(meshes))
    end do
    call check_unsafe_quadrature_rejection()
    growth_errors = abs(growth_values - screw_unstable_omega2) &
        / abs(screw_unstable_omega2)
    write (error_unit, "(a, 3es13.5)") "slow-point deviations ", slow_errors
    write (error_unit, "(a, 3es13.5)") "fast errors           ", fast_errors
    write (error_unit, "(a, 3es13.5)") "growth errors         ", growth_errors
    ! The exact slow band is 3.3e-7 wide relative above the slow point.
    ! The midpoint radial quadrature lets marginally resolved slow
    ! modes dip BELOW the essential infimum, and how far the worst
    ! mode dips is arithmetic sensitive across the near-degenerate
    ! cluster (measured 1.7e-4/2.2e-3/2.9e-2 at 32/64/128 locally,
    ! 2.3e-4/2.4e-3/7.0e-2 on CI); the envelope freezes the defect
    ! with cross-platform margin.  Exact quadrature (E5) and exact
    ! axis spaces (E2) are the queued fixes; the fast branch and the
    ! growth rate below carry the convergence claims.
    call require(all(slow_errors < [1.0e-3_dp, 1.0e-2_dp, 1.5e-1_dp]), &
        "theta-pinch slow edge leaves the frozen defect envelope")
    call require_convergent(fast_errors, 3.5_dp, &
        "theta-pinch fast branch is not second-order convergent")
    call require_convergent(second_errors, 3.5_dp, &
        "theta-pinch second fast mode is not second-order convergent")
    call require_convergent(growth_errors, 3.5_dp, &
        "screw-pinch growth rate is not second-order convergent")
    call require(fast_errors(size(meshes)) < 2.0e-4_dp, &
        "theta-pinch fast eigenvalue misses the Bessel reference")
    ! the resonant-layer constant is large (measured 61/16/3.9 percent
    ! at 32/64/128) but exactly second order (ratios 3.91 and 3.98);
    ! Richardson extrapolation of the two finest meshes lands on the
    ! shooting reference.
    call require(growth_errors(size(meshes)) < 6.0e-2_dp, &
        "screw-pinch growth rate misses the shooting reference")
    growth_extrapolated = (4.0_dp * growth_values(size(meshes)) &
        - growth_values(size(meshes) - 1)) / 3.0_dp
    write (error_unit, "(a, es13.5)") "extrapolated growth error ", &
        abs(growth_extrapolated - screw_unstable_omega2) &
        / abs(screw_unstable_omega2)
    call require(abs(growth_extrapolated - screw_unstable_omega2) &
        < 1.0e-2_dp * abs(screw_unstable_omega2), &
        "extrapolated growth rate misses the shooting reference")
    call write_results(slow_values, fast_values, second_values, growth_values)

    write (*, "(a)") "PASS"

contains

    subroutine assemble_pair(filename, surfaces, pressure_scale, &
            poloidal_scale, mode_m, mode_n, stiffness, mass, layout, &
            quadrature_points, reject_quadrature)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: surfaces, mode_m, mode_n
        real(dp), intent(in) :: pressure_scale, poloidal_scale
        real(dp), allocatable, intent(out) :: stiffness(:, :), mass(:, :)
        type(dynamic_family_layout_t), intent(out) :: layout
        integer, intent(in), optional :: quadrature_points
        logical, intent(in), optional :: reject_quadrature
        type(gvec_cas3d_equilibrium_t) :: equilibrium
        type(dynamic_family_layout_t) :: mass_layout
        type(mass_density_profile_t) :: density
        type(radial_space_config_t) :: radial_space
        real(dp), allocatable :: fields(:, :, :, :), drive(:, :, :)
        real(dp), allocatable :: jacobian_s(:, :, :), jacobian_t(:, :, :)
        real(dp), allocatable :: jacobian_z(:, :, :), gamma_p(:, :, :)
        real(dp) :: step
        integer :: info

        if (present(quadrature_points)) then
            radial_space%quadrature_points = quadrature_points
        end if
        call create_cylinder_fixture(filename, surfaces=surfaces, &
            pressure_scale=pressure_scale, poloidal_scale=poloidal_scale)
        call read_gvec_cas3d_file(filename, equilibrium, info)
        call require(info == reader_ok, "cylinder fixture was rejected")
        call delete_file(filename)
        call build_kernel_geometry(equilibrium, n_angle, n_angle, fields, &
            drive, info)
        call require(info == mercier_ok, "kernel geometry build failed")
        call build_compressible_geometry(equilibrium, n_angle, n_angle, &
            adiabatic_index, jacobian_s, jacobian_t, jacobian_z, gamma_p, &
            info)
        call require(info == compressible_geometry_ok, &
            "compressible geometry build failed")
        step = 1.0_dp / real(surfaces, dp)
        call assemble_compressible_family_stiffness(fields, drive, &
            jacobian_s, jacobian_t, jacobian_z, gamma_p, [mode_m], &
            [mode_n], [1], [0.0_dp], 1, radial_space, step, &
            phase_assembly_transformed, stiffness, layout, info)
        if (present(reject_quadrature)) then
            if (reject_quadrature) then
                call require(info /= 0, &
                    "unsafe interpolated quadrature was accepted")
                return
            end if
        end if
        call require(info == 0, "compressible stiffness assembly failed")
        density%s = [0.0_dp, 1.0_dp]
        density%kilograms_per_cubic_metre = [density_kg_m3, density_kg_m3]
        call assemble_physical_family_mass(fields, density, [mode_m], &
            [mode_n], [1], [0.0_dp], 1, radial_space, step, &
            phase_assembly_transformed, mass, mass_layout, info)
        call require(info == 0, "physical mass assembly failed")
        call require(mass_layout%total_unknowns == layout%total_unknowns, &
            "stiffness and mass layouts differ")
    end subroutine assemble_pair

    subroutine check_theta_pinch(surfaces, slow_error, fast_error, &
            second_error, slow_value, fast_value, second_value)
        integer, intent(in) :: surfaces
        real(dp), intent(out) :: slow_error, fast_error, second_error
        real(dp), intent(out) :: slow_value, fast_value, second_value
        type(dynamic_family_layout_t) :: layout
        real(dp), allocatable :: stiffness(:, :), mass(:, :)
        real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
        real(dp) :: gap_bottom, gap_top, fast_first, fast_second
        integer :: info, band_count, i

        ! mode (3,1): xi ~ r^2 = a^2 s at the axis is exactly
        ! representable under the CAS3D xi(0) = 0 rule; m = 1, 2 hit
        ! the documented axis-rule limitation (exact spaces are E2).
        call assemble_pair("theta_pinch_spectrum.nc", surfaces, 1.0_dp, &
            0.0_dp, 3, 1, stiffness, mass, layout)
        call solve_symmetric_generalized(stiffness, mass, eigenvalues, &
            eigenvectors, info)
        call require(info == 0, "theta-pinch dense solve failed")
        call require(all(eigenvalues > 0.0_dp), &
            "theta-pinch spectrum is not positive")
        slow_value = eigenvalues(1)
        slow_error = abs(slow_value - theta_slow_point) &
            / theta_slow_point

        ! the reference branch separations are factors of about 4.8e3
        ! (slow to Alfven) and 6.4e2 (Alfven to fast), so decade-wide
        ! windows detect pollution without pinning the a-priori
        ! unknown discrete cluster widths.
        gap_bottom = 2.0_dp * theta_slow_point
        gap_top = 0.5_dp * theta_alfven_point
        call require(count(eigenvalues > gap_bottom &
            .and. eigenvalues < gap_top) == 0, &
            "spectral pollution in the slow-Alfven gap")
        band_count = count(abs(eigenvalues - theta_alfven_point) &
            < 0.2_dp * theta_alfven_point)
        call require(band_count >= 1, &
            "Alfven cluster is missing from the discrete spectrum")

        fast_first = 0.0_dp
        fast_second = 0.0_dp
        do i = 1, size(eigenvalues) - 1
            if (eigenvalues(i) > 0.5_dp * theta_fast_lowest) then
                fast_first = eigenvalues(i)
                fast_second = eigenvalues(i + 1)
                exit
            end if
        end do
        call require(fast_first > 0.0_dp, "fast branch is missing")
        call require(count(eigenvalues > 2.0_dp * theta_alfven_point &
            .and. eigenvalues < 0.5_dp * theta_fast_lowest) == 0, &
            "spectral pollution in the Alfven-fast gap")
        fast_error = abs(fast_first - theta_fast_lowest) &
            / theta_fast_lowest
        second_error = abs(fast_second - theta_fast_second) &
            / theta_fast_second
        fast_value = fast_first
        second_value = fast_second
    end subroutine check_theta_pinch

    subroutine validate_arguments()
        integer :: io_status

        if (command_argument_count() > 1) then
            call fail_cli("usage: test_cylinder_physical_spectrum " // &
                "[OUTPUT_CSV]")
        end if
        output_requested = command_argument_count() == 1
        if (.not. output_requested) return
        call get_command_argument(1, output_path, status=io_status)
        if (io_status /= 0) call fail_cli("output path is unavailable")
        if (len_trim(output_path) == 0) call fail_cli("output path is empty")
    end subroutine validate_arguments

    subroutine write_results(slow, fast, second, growth)
        real(dp), intent(in) :: slow(:), fast(:), second(:), growth(:)
        integer :: i, io_status, unit

        if (.not. output_requested) return
        open (newunit=unit, file=trim(output_path), status="replace", &
            action="write", iostat=io_status)
        if (io_status /= 0) call fail_cli("output CSV cannot be opened")
        write (unit, "(a)") "n_radial,theta_slow_s_minus_2," // &
            "theta_fast_first_s_minus_2,theta_fast_second_s_minus_2," // &
            "screw_unstable_omega2_s_minus_2"
        do i = 1, size(meshes)
            write (unit, "(i0,4(',',es24.16))", iostat=io_status) &
                meshes(i), slow(i), fast(i), second(i), growth(i)
            if (io_status /= 0) call fail_cli("output CSV write failed")
        end do
        close (unit, iostat=io_status)
        if (io_status /= 0) call fail_cli("output CSV close failed")
    end subroutine write_results

    subroutine fail_cli(message)
        character(len=*), intent(in) :: message

        write (error_unit, "(a)") &
            "test_cylinder_physical_spectrum: " // message
        flush (error_unit)
        stop 2
    end subroutine fail_cli

    subroutine check_screw_stable(surfaces)
        integer, intent(in) :: surfaces
        type(dynamic_family_layout_t) :: layout
        real(dp), allocatable :: stiffness(:, :), mass(:, :)
        real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
        integer :: info

        ! resonant (4,4) member: iota crosses 1, so the slow continuum
        ! touches zero and the sign of near-zero discrete modes is
        ! arithmetic noise; stability is asserted above a self-scaled
        ! roundoff floor.
        call assemble_pair("screw_stable_spectrum.nc", surfaces, 0.25_dp, &
            1.0_dp, 4, 4, stiffness, mass, layout)
        call solve_symmetric_generalized(stiffness, mass, eigenvalues, &
            eigenvectors, info)
        call require(info == 0, "stable screw dense solve failed")
        call require(all(eigenvalues > -1.0e-10_dp &
            * maxval(abs(eigenvalues))), &
            "Suydam-stable member is unstable beyond the roundoff floor")

        ! non-resonant (2,1) member: F stays positive, so no discrete
        ! mode may fall below the slow continuum edge (variational
        ! bound; 0.98 covers quadrature nonvariationality).
        call assemble_pair("screw_nonres_spectrum.nc", surfaces, 0.25_dp, &
            1.0_dp, 2, 1, stiffness, mass, layout)
        call solve_symmetric_generalized(stiffness, mass, eigenvalues, &
            eigenvectors, info)
        call require(info == 0, "non-resonant screw dense solve failed")
        call require(all(eigenvalues > 0.0_dp), &
            "non-resonant stable member has a negative eigenvalue")
        call require(count(eigenvalues &
            < 0.98_dp * screw_stable_slow_edge) == 0, &
            "stable member has a mode below the slow continuum edge")
    end subroutine check_screw_stable

    subroutine check_screw_unstable(surfaces, growth_value, certify)
        integer, intent(in) :: surfaces
        real(dp), intent(out) :: growth_value
        logical, intent(in) :: certify
        type(dynamic_family_layout_t) :: layout
        real(dp), allocatable :: stiffness(:, :), mass(:, :)
        real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
        integer :: info

        call assemble_pair("screw_unstable_spectrum.nc", surfaces, 1.0_dp, &
            1.0_dp, 4, 4, stiffness, mass, layout)
        call solve_symmetric_generalized(stiffness, mass, eigenvalues, &
            eigenvectors, info)
        call require(info == 0, "unstable screw dense solve failed")
        call require(eigenvalues(1) < 0.0_dp, &
            "Suydam-unstable member is not unstable")
        growth_value = eigenvalues(1)
        if (certify) call certify_unstable_pair(stiffness, mass, layout, &
            eigenvalues)
    end subroutine check_screw_unstable

    subroutine certify_unstable_pair(stiffness, mass, layout, eigenvalues)
        real(dp), intent(in) :: stiffness(:, :), mass(:, :)
        type(dynamic_family_layout_t), intent(in) :: layout
        real(dp), intent(in) :: eigenvalues(:)
        type(variable_block_tridiagonal_t) :: k_blocks, m_blocks
        real(dp), allocatable :: vector(:)
        integer, allocatable :: permutation(:), widths(:)
        real(dp) :: eigenvalue, residual, resolution
        integer :: info, negatives

        call build_dynamic_block_permutation(layout, widths, permutation, &
            info)
        call require(info == dynamic_layout_ok, &
            "block permutation build failed")
        call pack_permuted_variable_blocks(stiffness, permutation, widths, &
            k_blocks, info)
        call require(info == variable_block_ok, "stiffness packing failed")
        call pack_permuted_variable_blocks(mass, permutation, widths, &
            m_blocks, info)
        call require(info == variable_block_ok, "mass packing failed")
        call variable_generalized_inertia(k_blocks, m_blocks, 0.0_dp, &
            negatives, info)
        call require(info == variable_generalized_ok, &
            "variable-block inertia failed")
        call require(negatives == count(eigenvalues < 0.0_dp), &
            "variable-block inertia disagrees with the dense oracle")
        call iterate_variable_generalized_eigenvalue(k_blocks, m_blocks, &
            1.05_dp * eigenvalues(1), eigenvalue, vector, residual, &
            resolution, info)
        call require(info == variable_generalized_ok, &
            "variable-block inverse iteration failed")
        call require(abs(eigenvalue - eigenvalues(1)) &
            <= residual + resolution + 1.0e-9_dp * abs(eigenvalues(1)), &
            "certified eigenvalue disagrees with the dense oracle")
    end subroutine certify_unstable_pair

    subroutine check_unsafe_quadrature_rejection()
        type(dynamic_family_layout_t) :: layout
        real(dp), allocatable :: stiffness(:, :), mass(:, :)

        call assemble_pair("rejected_quadrature.nc", meshes(1), 1.0_dp, &
            0.0_dp, 3, 1, stiffness, mass, layout, quadrature_points=2, &
            reject_quadrature=.true.)
    end subroutine check_unsafe_quadrature_rejection

    subroutine check_bridge_rejections()
        type(gvec_cas3d_equilibrium_t) :: equilibrium
        real(dp), allocatable :: jacobian_s(:, :, :), jacobian_t(:, :, :)
        real(dp), allocatable :: jacobian_z(:, :, :), gamma_p(:, :, :)
        integer :: info

        call create_cylinder_fixture("bridge_reject.nc", surfaces=16)
        call read_gvec_cas3d_file("bridge_reject.nc", equilibrium, info)
        call require(info == reader_ok, "rejection fixture was rejected")
        call delete_file("bridge_reject.nc")
        call build_compressible_geometry(equilibrium, 4, n_angle, &
            adiabatic_index, jacobian_s, jacobian_t, jacobian_z, gamma_p, &
            info)
        call require(info == compressible_geometry_invalid_input, &
            "coarse angular grid was not rejected")
        call build_compressible_geometry(equilibrium, n_angle, n_angle, &
            -1.0_dp, jacobian_s, jacobian_t, jacobian_z, gamma_p, info)
        call require(info == compressible_geometry_invalid_input, &
            "negative adiabatic index was not rejected")
        ! edge spline undershoot within the relative noise floor is
        ! accepted and clamped out of gamma*p; undershoot beyond the
        ! floor is rejected.
        equilibrium%pressure(size(equilibrium%pressure)) = &
            -1.0e-4_dp * maxval(equilibrium%pressure)
        call build_compressible_geometry(equilibrium, n_angle, n_angle, &
            adiabatic_index, jacobian_s, jacobian_t, jacobian_z, gamma_p, &
            info)
        call require(info == compressible_geometry_ok, &
            "sub-floor pressure undershoot was rejected")
        call require(all(gamma_p(:, :, size(equilibrium%pressure)) &
            == 0.0_dp), "pressure undershoot was not clamped")
        equilibrium%pressure(1) = -0.01_dp * maxval(equilibrium%pressure)
        call build_compressible_geometry(equilibrium, n_angle, n_angle, &
            adiabatic_index, jacobian_s, jacobian_t, jacobian_z, gamma_p, &
            info)
        call require(info == compressible_geometry_invalid_input, &
            "negative pressure beyond the floor was not rejected")
        ! vacuum equilibria carry exactly zero pressure
        equilibrium%pressure = 0.0_dp
        call build_compressible_geometry(equilibrium, n_angle, n_angle, &
            adiabatic_index, jacobian_s, jacobian_t, jacobian_z, gamma_p, &
            info)
        call require(info == compressible_geometry_ok, &
            "vacuum zero pressure was rejected")
        call require(all(gamma_p == 0.0_dp), &
            "vacuum gamma*p is not identically zero")
    end subroutine check_bridge_rejections

    subroutine require_convergent(errors, ratio, message)
        real(dp), intent(in) :: errors(:), ratio
        character(len=*), intent(in) :: message

        call require(all(errors(2:) < errors(:size(errors) - 1)), message)
        call require(all(errors(:size(errors) - 1) / errors(2:) > ratio), &
            message)
    end subroutine require_convergent

    subroutine delete_file(filename)
        character(len=*), intent(in) :: filename
        integer :: unit

        open (newunit=unit, file=filename, status="old")
        close (unit, status="delete")
    end subroutine delete_file

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_cylinder_physical_spectrum
