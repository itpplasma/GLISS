! Conforming axis space of the compatible FEEC problems (#35).
!
! A Cartesian-regular |m|=1 displacement has xi^s ~ a s^(1/2) and
! eta ~ b s^(-1/2) with b fixed by a. The discrete space ties the first
! eta coefficient to the first normal coefficient with the factor kappa of
! compatible_axis_regularity. The oracle is the kernel itself: on the
! analytic screw-pinch fixture the energy density of the tied pair stays
! bounded as s -> 0 in the incompressible and the compressible kernel and
! both parities, while the untied normal function alone grows like 1/s (a
! logarithmically divergent energy). On the toroidal Solov'ev export
! (shifted, elongated axis) the same cancellation holds to the regularity
! of the exported Boozer chart: the tied density at s = 1e-8 is below
! 1e-8 of the untied one. The quadrature and block-compaction helpers are
! checked against exact values.
program test_axis_regularity
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use compatible_axis_regularity, only: axis_flux_slope, &
        axis_regularity_ok, axis_tie_factor
    use compatible_block_storage, only: compatible_block_ok, &
        initialize_compatible_block_pencil
    use compatible_compressible_stiffness_assembly, only: &
        assemble_compatible_compressible_stiffness_surface, &
        compatible_stiffness_term_count
    use compatible_family_point_assembly, only: &
        assemble_compatible_transformed_surface
    use cylinder_fixture, only: create_cylinder_fixture
    use compatible_radial_quadrature, only: build_axis_quadrature, &
        build_gauss_legendre, compatible_quadrature_ok
    use export_surface_geometry, only: build_angular_grids
    use gvec_cas3d_reader, only: read_gvec_cas3d_file, reader_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use phase_assembly_policy, only: phase_assembly_transformed
    use primitive_equilibrium_spline, only: fit_primitive_equilibrium, &
        primitive_equilibrium_ok, primitive_equilibrium_spline_t
    use primitive_kernel_geometry, only: evaluate_primitive_kernel_surface, &
        primitive_kernel_ok
    use radial_feec_complex, only: build_radial_feec_complex, &
        radial_feec_complex_t, radial_feec_ok
    use variable_block_tridiagonal, only: variable_block_tridiagonal_t
    implicit none

    character(len=1024) :: directory
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(primitive_equilibrium_spline_t) :: spline
    type(radial_feec_complex_t) :: complex
    real(dp) :: flux_slope
    integer :: info, parity

    call get_command_argument(1, directory)
    if (len_trim(directory) == 0) error stop 'supply the fixture directory'
    ! Degree one on unit cells: the first retained H1 function is N = s near
    ! the axis, so the tie factor applies to xi^s = s^(-1/2) N.
    call build_radial_feec_complex([0.0_dp, 1.0_dp, 2.0_dp], 1, .true., &
        .true., complex, info)
    call require(info == radial_feec_ok, 'FEEC complex failed')
    call create_cylinder_fixture('axis_regularity_cylinder.nc')
    call load('axis_regularity_cylinder.nc')
    do parity = 1, 2
        call check_bounded_density(parity, &
            axis_tie_factor(complex, parity, flux_slope), .true.)
    end do
    call load(trim(directory) // '/solovev_q1.035.nc')
    do parity = 1, 2
        call check_bounded_density(parity, &
            axis_tie_factor(complex, parity, flux_slope), .false.)
    end do
    call check_quadrature()
    call check_block_compaction()
    write (*, '(a)') 'axis regularity checks passed'

contains

    subroutine load(path)
        character(len=*), intent(in) :: path

        call read_gvec_cas3d_file(path, equilibrium, info)
        call require(info == reader_ok, 'fixture read failed')
        call fit_primitive_equilibrium(equilibrium, spline, info)
        call require(info == primitive_equilibrium_ok, 'spline fit failed')
        call axis_flux_slope(spline, flux_slope, info)
        call require(info == axis_regularity_ok, 'axis flux slope failed')
    end subroutine load

    subroutine check_bounded_density(parity, kappa, analytic)
        integer, intent(in) :: parity
        real(dp), intent(in) :: kappa
        logical, intent(in) :: analytic
        real(dp) :: tied(2, 3), untied(2, 3)
        integer :: level

        do level = 1, 3
            call densities(parity, kappa, 10.0_dp**(-2 * level - 2), &
                tied(:, level), untied(:, level))
        end do
        ! s drops by 100 per level: the untied normal function alone grows
        ! by the same factor, the tied densities converge.
        call require(all(untied(:, 3) >= 50.0_dp * untied(:, 2)), &
            'untied |m|=1 compression does not diverge like 1/s')
        if (analytic) then
            call require(all(abs(tied(:, 3) - tied(:, 2)) &
                <= 1.0e-2_dp * abs(tied(:, 2))), &
                'tied |m|=1 energy density is not bounded at the axis')
        else
            call require(abs(tied(1, 3) - tied(1, 2)) &
                <= 1.0e-2_dp * abs(tied(1, 2)), &
                'tied |m|=1 incompressible density is not bounded')
            call require(all(abs(tied(:, 3)) <= 1.0e-8_dp * untied(:, 3)), &
                'tied |m|=1 compression does not cancel at the axis')
        end if
    end subroutine check_bounded_density

    subroutine densities(parity, kappa, s, tied, untied)
        integer, intent(in) :: parity
        real(dp), intent(in) :: kappa, s
        real(dp), intent(out) :: tied(2), untied(2)
        real(dp), allocatable :: theta(:), zeta(:), fields(:, :, :)
        real(dp), allocatable :: drive(:, :), jacobian_s(:, :)
        real(dp), allocatable :: jacobian_t(:, :), jacobian_z(:, :)
        real(dp), allocatable :: gamma_p(:, :)
        real(dp) :: incompressible(2, 2), terms2(2, 2, 4)
        real(dp) :: compressible(3, 3), terms3(3, 3, &
            compatible_stiffness_term_count)
        real(dp) :: h1(1, 1), dh1(1, 1), l2(1, 1), eta(1, 1), zero(1, 1)
        real(dp) :: pressure, pair(2), triple(3), image2(2), image3(3)
        integer :: status, trial_m(1), trial_n(1), parities(1)

        call build_angular_grids(64, 8, theta, zeta)
        call evaluate_primitive_kernel_surface(spline, s, theta, zeta, &
            fields, drive, status, jacobian_s, jacobian_t, jacobian_z, &
            pressure)
        call require(status == primitive_kernel_ok, 'kernel surface failed')
        h1 = sqrt(s)
        dh1 = 0.5_dp / sqrt(s)
        l2 = 1.0_dp / sqrt(s)
        eta = l2
        zero = 0.0_dp
        trial_m = 1
        trial_n = 1
        parities = parity
        incompressible = 0.0_dp
        terms2 = 0.0_dp
        call assemble_compatible_transformed_surface(fields, drive, trial_m, &
            trial_n, parities, equilibrium%field_periods, h1, dh1, l2, &
            incompressible, status, terms2)
        call require(status == 0, 'incompressible kernel failed')
        pair(1) = 1.0_dp
        pair(2) = kappa
        image2 = matmul(incompressible, pair)
        tied(1) = dot_product(pair, image2)
        untied(1) = incompressible(1, 1)
        allocate (gamma_p, mold=drive)
        gamma_p = 5.0_dp / 3.0_dp * pressure
        compressible = 0.0_dp
        terms3 = 0.0_dp
        call assemble_compatible_compressible_stiffness_surface(fields, &
            drive, jacobian_s, jacobian_t, jacobian_z, gamma_p, trial_m, &
            trial_n, parities, equilibrium%field_periods, h1, dh1, eta, zero, &
            1.0_dp, phase_assembly_transformed, compressible, status, terms3)
        call require(status == 0, 'compressible kernel failed')
        triple(1) = 1.0_dp
        triple(2) = kappa
        triple(3) = 0.0_dp
        image3 = matmul(compressible, triple)
        tied(2) = dot_product(triple, image3)
        untied(2) = compressible(1, 1)
    end subroutine densities

    subroutine check_quadrature()
        real(dp), allocatable :: nodes(:), weights(:)
        real(dp) :: exact, total
        integer :: power, points, status

        do points = 1, 24
            call build_gauss_legendre(points, nodes, weights, status)
            call require(status == compatible_quadrature_ok, &
                'Gauss-Legendre rule failed')
            do power = 0, 2 * points - 1
                exact = merge(2.0_dp / real(power + 1, dp), 0.0_dp, &
                    modulo(power, 2) == 0)
                total = sum(weights * nodes**power)
                call require(abs(total - exact) <= 64.0_dp &
                    * epsilon(1.0_dp), 'Gauss-Legendre rule is not exact')
            end do
        end do
        ! The axis rule in sqrt(s) integrates s^(k/2) on [0, h] exactly for
        ! 2 points -> k + 1 <= 2 * 2 * points - 1.
        points = 7
        call build_axis_quadrature(0.25_dp, points, nodes, weights, status)
        call require(status == compatible_quadrature_ok, 'axis rule failed')
        do power = 0, 2 * points - 2
            exact = 0.25_dp**(0.5_dp * power + 1.0_dp) &
                / (0.5_dp * power + 1.0_dp)
            total = sum(weights * nodes**(0.5_dp * power))
            call require(abs(total - exact) <= 64.0_dp * epsilon(1.0_dp) &
                * exact, 'axis rule is not exact for half-integer powers')
        end do
        call require(all(nodes > 0.0_dp .and. nodes < 0.25_dp), &
            'axis rule nodes leave the axis element')
    end subroutine check_quadrature

    subroutine check_block_compaction()
        type(variable_block_tridiagonal_t) :: stiffness, mass
        integer, allocatable :: block_index(:), local_index(:)
        logical :: eliminated(22)
        integer :: status

        ! Degree 2 on three cells: 5 H1 and 6 L2 functions, two normal and
        ! two eta trials -> blocks [H1 1:2, L2 1:2], [H1 3:4, L2 3:4],
        ! [H1 5, L2 5:6]; widths 8, 8, 6 over 22 unknowns before removing
        ! unknowns 11 and 12 (eta function 1 of both trials, block 1
        ! locals 5 and 6).
        eliminated = .false.
        eliminated([11, 12]) = .true.
        call initialize_compatible_block_pencil(5, 6, 2, 2, 2, stiffness, &
            mass, block_index, local_index, status)
        call require(status == compatible_block_ok, 'block pencil failed')
        call require(all(stiffness%widths == [8, 8, 6]), &
            'full block widths differ')
        call initialize_compatible_block_pencil(5, 6, 2, 2, 2, stiffness, &
            mass, block_index, local_index, status, eliminated(1:15))
        call require(status /= compatible_block_ok, &
            'mismatched elimination mask was accepted')
        call initialize_compatible_block_pencil(5, 6, 2, 2, 2, stiffness, &
            mass, block_index, local_index, status, eliminated)
        call require(status == compatible_block_ok, &
            'compacted block pencil failed')
        call require(all(stiffness%widths == [6, 8, 6]), &
            'compacted block widths differ')
        call require(size(block_index) == 20, 'compacted size differs')
        ! Normal unknowns 1:4 keep block 1 locals 1:4; the first kept eta
        ! unknowns (full 13, 14 = eta function 2) move to locals 5, 6.
        call require(all(block_index(1:4) == 1) &
            .and. all(local_index(1:4) == [1, 2, 3, 4]), &
            'normal locals changed')
        call require(all(block_index(11:12) == 1) &
            .and. all(local_index(11:12) == [5, 6]), &
            'compacted eta locals differ')
        call require(all(block_index(13:14) == 2) &
            .and. all(local_index(13:14) == [5, 6]), &
            'second block eta locals differ')
    end subroutine check_block_compaction

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (*, '(a)') 'FAIL: ' // message
            error stop 1
        end if
    end subroutine require

end program test_axis_regularity
