program test_terpsichore_noninteracting_stiffness
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use dynamic_family_layout, only: dynamic_family_layout_t
    use terpsichore_matrix_fixture, only: terpsichore_matrix_fixture_t
    use terpsichore_noninteracting_coefficients, only: &
        build_terpsichore_noninteracting_coefficients, &
        build_terpsichore_noninteracting_coefficients_direct, &
        terpsichore_coefficients_ok
    use terpsichore_noninteracting_stiffness, only: &
        assemble_terpsichore_noninteracting_fixed_boundary_stiffness, &
        terpsichore_noninteracting_ok
    implicit none

    type(terpsichore_matrix_fixture_t) :: fixture
    type(dynamic_family_layout_t) :: layout
    real(dp), allocatable :: expected(:, :), stiffness(:, :)
    real(dp), allocatable :: fast(:, :, :, :), direct(:, :, :, :)
    integer :: info

    call build_constant_fixture(fixture)
    call assemble_terpsichore_noninteracting_fixed_boundary_stiffness( &
        fixture, stiffness, layout, info)
    call require(info == terpsichore_noninteracting_ok, &
        "constant non-interacting fixture was rejected")
    call build_expected(expected)
    call require(layout%total_unknowns == 5, &
        "non-interacting fixed-boundary layout is wrong")
    call require(maxval(abs(stiffness - expected)) < 2.0e-12_dp, &
        "non-interacting analytical matrix is wrong")
    call build_terpsichore_noninteracting_coefficients(fixture, fast, &
        info)
    call require(info == terpsichore_coefficients_ok, &
        "transform coefficient path failed")
    call build_terpsichore_noninteracting_coefficients_direct(fixture, &
        direct, info)
    call require(info == terpsichore_coefficients_ok, &
        "direct coefficient oracle failed")
    call require(maxval(abs(fast - direct)) < 1.0e-13_dp &
        * max(1.0_dp, maxval(abs(direct))), &
        "transform coefficients disagree with the point-pair oracle")
    fixture%legacy_modelk = 2
    call assemble_terpsichore_noninteracting_fixed_boundary_stiffness( &
        fixture, stiffness, layout, info)
    call require(info == terpsichore_noninteracting_ok, &
        "physical-norm non-interacting potential was rejected")
    call require(maxval(abs(stiffness - expected)) < 2.0e-12_dp, &
        "kinetic norm changed the non-interacting potential")
    fixture%legacy_modelk = 1
    call assemble_terpsichore_noninteracting_fixed_boundary_stiffness( &
        fixture, stiffness, layout, info)
    call require(info /= terpsichore_noninteracting_ok, &
        "Kruskal-Oberman potential entered the non-interacting kernel")
    fixture%legacy_modelk = 0
    fixture%sigma_b(1, 1) = tiny(1.0_dp)
    call assemble_terpsichore_noninteracting_fixed_boundary_stiffness( &
        fixture, stiffness, layout, info)
    call require(info /= terpsichore_noninteracting_ok, &
        "numerically singular sigma_b entered the coefficient kernel")
    fixture%sigma_b(1, 1) = -1.0_dp
    fixture%current_j(1) = tiny(1.0_dp)
    call assemble_terpsichore_noninteracting_fixed_boundary_stiffness( &
        fixture, stiffness, layout, info)
    call require(info /= terpsichore_noninteracting_ok, &
        "numerically singular current flux entered the stiffness kernel")
    write (*, "(a)") "PASS"

contains

    subroutine build_constant_fixture(value)
        type(terpsichore_matrix_fixture_t), intent(out) :: value

        value%intervals = 3
        value%poloidal_points = 8
        value%toroidal_points = 1
        value%stability_periods = 1
        value%field_periods = 1
        value%modes = 1
        value%legacy_modelk = 0
        value%parity = 0.0_dp
        value%current_factor = 1.0_dp
        allocate (value%s(0:3))
        value%s = [0.0_dp, 0.2_dp, 0.6_dp, 1.0_dp]
        value%flux_p_slope = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
        value%flux_t_slope = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
        allocate (value%flux_p_curve(4), source=0.0_dp)
        allocate (value%flux_t_curve(4), source=0.0_dp)
        allocate (value%current_i(3), source=0.0_dp)
        allocate (value%current_j(3), source=1.0_dp)
        allocate (value%pressure_slope(3), source=0.0_dp)
        value%radial_power = [0.0_dp]
        value%mode_m = [1]
        value%mode_n = [0]
        allocate (value%signed_bjac(8, 0:3), source=-1.0_dp)
        value%signed_bjac(:, 0) = 0.0_dp
        allocate (value%signed_bjac_radial(8, 0:3), source=0.0_dp)
        allocate (value%sigma_b_s(8, 0:3), source=0.0_dp)
        allocate (value%metric_ss_over_jacobian(8, 0:3), source=-1.0_dp)
        allocate (value%metric_st_over_jacobian(8, 0:3), source=0.0_dp)
        allocate (value%metric_tt_over_jacobian(8, 0:3), source=-2.0_dp)
        allocate (value%sigma_b(8, 3), source=-1.0_dp)
        allocate (value%parallel_current(8, 3), source=0.0_dp)
    end subroutine build_constant_fixture

    subroutine build_expected(matrix)
        real(dp), allocatable, intent(out) :: matrix(:, :)
        real(dp) :: diagonal(3), off_diagonal, radial_step(3)
        integer :: interval

        radial_step = [0.2_dp, 0.4_dp, 0.4_dp]
        do interval = 1, 3
            diagonal(interval) = 3.0_dp / radial_step(interval) &
                + 0.75_dp * radial_step(interval)
        end do
        off_diagonal = -3.0_dp / radial_step(2) &
            + 0.75_dp * radial_step(2)
        allocate (matrix(5, 5), source=0.0_dp)
        matrix(1, 1) = diagonal(1) + diagonal(2)
        matrix(2, 2) = diagonal(2) + diagonal(3)
        matrix(1, 2) = off_diagonal
        matrix(2, 1) = off_diagonal
        do interval = 1, 3
            matrix(2 + interval, 2 + interval) = 3.0_dp * radial_step(interval)
        end do
    end subroutine build_expected

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") message
        error stop 1
    end subroutine require

end program test_terpsichore_noninteracting_stiffness
