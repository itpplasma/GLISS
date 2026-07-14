program test_radial_cubic_spline
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use radial_cubic_spline, only: build_radial_cubic_spline_grid, &
        evaluate_radial_cubic_spline, fit_radial_cubic_spline, &
        radial_cubic_spline_grid_t, radial_cubic_spline_invalid, &
        radial_cubic_spline_ok, radial_cubic_spline_t
    implicit none

    real(dp), parameter :: nodes(6) = [0.04_dp, 0.17_dp, 0.36_dp, &
        0.58_dp, 0.79_dp, 0.96_dp]
    real(dp), parameter :: points(9) = [0.0_dp, 0.04_dp, 0.11_dp, &
        0.36_dp, 0.5_dp, 0.79_dp, 0.9_dp, 0.96_dp, 1.0_dp]
    type(radial_cubic_spline_grid_t) :: grid
    type(radial_cubic_spline_t) :: spline
    real(dp) :: derivative, value
    integer :: i, info

    call build_radial_cubic_spline_grid(nodes, 0.0_dp, 1.0_dp, grid, info)
    call require(info == radial_cubic_spline_ok, "valid grid was rejected")
    call fit_radial_cubic_spline(grid, polynomial(nodes), spline, info)
    call require(info == radial_cubic_spline_ok, "cubic fit failed")
    do i = 1, size(points)
        call evaluate_radial_cubic_spline(grid, spline, points(i), value, &
            derivative, info)
        call require(info == radial_cubic_spline_ok, "evaluation failed")
        call require(abs(value - polynomial(points(i))) < 2.0e-13_dp, &
            "cubic value is not exact")
        call require(abs(derivative - polynomial_slope(points(i))) &
            < 2.0e-12_dp, "cubic derivative is not exact")
    end do
    call check_continuity(grid)
    call check_mathematica_fixture()
    call check_invalid_inputs(grid, spline)
    write (*, "(a)") "PASS"

contains

    elemental function polynomial(x) result(value)
        real(dp), intent(in) :: x
        real(dp) :: value

        value = 1.0_dp - 2.0_dp * x + 3.0_dp * x**2 - 0.5_dp * x**3
    end function polynomial

    elemental function polynomial_slope(x) result(value)
        real(dp), intent(in) :: x
        real(dp) :: value

        value = -2.0_dp + 6.0_dp * x - 1.5_dp * x**2
    end function polynomial_slope

    subroutine check_continuity(valid_grid)
        type(radial_cubic_spline_grid_t), intent(in) :: valid_grid
        type(radial_cubic_spline_t) :: smooth
        real(dp) :: left_slope, right_slope, second, value
        real(dp) :: first_third, last_third, scale
        integer :: interval, status

        call fit_radial_cubic_spline(valid_grid, sin(3.0_dp * nodes), smooth, &
            status)
        call require(status == radial_cubic_spline_ok, "smooth fit failed")
        scale = max(1.0_dp, maxval(abs(smooth%second_derivatives)))
        first_third = (smooth%second_derivatives(2) &
            - smooth%second_derivatives(1)) / valid_grid%intervals(1)
        last_third = (smooth%second_derivatives(3) &
            - smooth%second_derivatives(2)) / valid_grid%intervals(2)
        call require(abs(first_third - last_third) < 1.0e-12_dp * scale, &
            "left not-a-knot condition failed")
        first_third = (smooth%second_derivatives(size(nodes) - 1) &
            - smooth%second_derivatives(size(nodes) - 2)) &
            / valid_grid%intervals(size(nodes) - 2)
        last_third = (smooth%second_derivatives(size(nodes)) &
            - smooth%second_derivatives(size(nodes) - 1)) &
            / valid_grid%intervals(size(nodes) - 1)
        call require(abs(first_third - last_third) < 1.0e-12_dp * scale, &
            "right not-a-knot condition failed")
        do interval = 1, size(nodes) - 2
            left_slope = (smooth%values(interval + 1) &
                - smooth%values(interval)) / valid_grid%intervals(interval) &
                + valid_grid%intervals(interval) &
                * (smooth%second_derivatives(interval) &
                + 2.0_dp * smooth%second_derivatives(interval + 1)) / 6.0_dp
            right_slope = (smooth%values(interval + 2) &
                - smooth%values(interval + 1)) &
                / valid_grid%intervals(interval + 1) &
                - valid_grid%intervals(interval + 1) &
                * (2.0_dp * smooth%second_derivatives(interval + 1) &
                + smooth%second_derivatives(interval + 2)) / 6.0_dp
            call require(abs(left_slope - right_slope) < 1.0e-12_dp * scale, &
                "first derivative is discontinuous")
            call evaluate_radial_cubic_spline(valid_grid, smooth, &
                nodes(interval + 1), value, left_slope, status, second)
            call require(status == radial_cubic_spline_ok .and. &
                abs(second - smooth%second_derivatives(interval + 1)) &
                < 1.0e-13_dp * scale, "second derivative is discontinuous")
        end do
    end subroutine check_continuity

    subroutine check_mathematica_fixture()
        real(dp), parameter :: fixture_nodes(6) = [0.04_dp, 0.17_dp, &
            0.36_dp, 0.58_dp, 0.79_dp, 0.96_dp]
        real(dp), parameter :: fixture_values(6) = [1.0_dp, -2.0_dp, &
            3.0_dp, 0.5_dp, -1.0_dp, 2.0_dp]
        real(dp), parameter :: query(3) = [0.0_dp, 0.5_dp, 1.0_dp]
        real(dp), parameter :: expected_value(3) = [ &
            4.6555397308640102_dp, 2.1105063103881973_dp, &
            3.3526179251273755_dp]
        real(dp), parameter :: expected_slope(3) = [ &
            -112.13643631143828_dp, -19.278249093642739_dp, &
            37.130255510230877_dp]
        type(radial_cubic_spline_grid_t) :: fixture_grid
        type(radial_cubic_spline_t) :: fixture_spline
        real(dp) :: actual_slope, actual_value
        integer :: i, status

        call build_radial_cubic_spline_grid(fixture_nodes, 0.0_dp, 1.0_dp, &
            fixture_grid, status)
        call require(status == radial_cubic_spline_ok, &
            "Mathematica fixture grid failed")
        call fit_radial_cubic_spline(fixture_grid, fixture_values, &
            fixture_spline, status)
        call require(status == radial_cubic_spline_ok, &
            "Mathematica fixture fit failed")
        do i = 1, size(query)
            call evaluate_radial_cubic_spline(fixture_grid, fixture_spline, &
                query(i), actual_value, actual_slope, status)
            call require(status == radial_cubic_spline_ok, &
                "Mathematica fixture evaluation failed")
            call require(abs(actual_value - expected_value(i)) &
                < 2.0e-13_dp * max(1.0_dp, abs(expected_value(i))), &
                "Mathematica fixture value differs")
            call require(abs(actual_slope - expected_slope(i)) &
                < 2.0e-13_dp * max(1.0_dp, abs(expected_slope(i))), &
                "Mathematica fixture slope differs")
        end do
    end subroutine check_mathematica_fixture

    subroutine check_invalid_inputs(valid_grid, valid_spline)
        type(radial_cubic_spline_grid_t), intent(in) :: valid_grid
        type(radial_cubic_spline_t), intent(in) :: valid_spline
        type(radial_cubic_spline_grid_t) :: local_grid
        type(radial_cubic_spline_t) :: local_spline
        real(dp) :: local_nodes(6), local_values(6), result, slope
        integer :: status

        local_nodes = nodes
        local_nodes(3) = local_nodes(2)
        call build_radial_cubic_spline_grid(local_nodes, 0.0_dp, 1.0_dp, &
            local_grid, status)
        call require(status == radial_cubic_spline_invalid, &
            "repeated radial node was accepted")
        local_nodes = nodes
        local_nodes(2) = ieee_value(0.0_dp, ieee_quiet_nan)
        call build_radial_cubic_spline_grid(local_nodes, 0.0_dp, 1.0_dp, &
            local_grid, status)
        call require(status == radial_cubic_spline_invalid, &
            "nonfinite radial node was accepted")
        call build_radial_cubic_spline_grid(nodes, 0.1_dp, 1.0_dp, &
            local_grid, status)
        call require(status == radial_cubic_spline_invalid, &
            "sample outside reconstruction domain was accepted")
        call fit_radial_cubic_spline(valid_grid, polynomial(nodes(:5)), &
            local_spline, status)
        call require(status == radial_cubic_spline_invalid, &
            "wrong value count was accepted")
        local_values = polynomial(nodes)
        local_values(4) = ieee_value(0.0_dp, ieee_quiet_nan)
        call fit_radial_cubic_spline(valid_grid, local_values, local_spline, &
            status)
        call require(status == radial_cubic_spline_invalid, &
            "nonfinite value was accepted")
        call evaluate_radial_cubic_spline(valid_grid, valid_spline, -0.01_dp, &
            result, slope, status)
        call require(status == radial_cubic_spline_invalid, &
            "coordinate outside domain was accepted")
    end subroutine check_invalid_inputs

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_radial_cubic_spline
