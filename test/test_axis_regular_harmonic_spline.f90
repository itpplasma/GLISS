program test_axis_regular_harmonic_spline
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use axis_regular_harmonic_spline, only: &
        axis_regular_harmonic_field_t, axis_regular_harmonic_invalid, &
        axis_regular_harmonic_ok, evaluate_axis_regular_harmonics, &
        fit_axis_regular_harmonics
    use radial_cubic_spline, only: build_radial_cubic_spline_grid, &
        radial_cubic_spline_grid_t, radial_cubic_spline_ok
    implicit none

    real(dp), parameter :: nodes(6) = [0.04_dp, 0.17_dp, 0.36_dp, &
        0.58_dp, 0.79_dp, 0.96_dp]
    type(radial_cubic_spline_grid_t) :: grid
    integer :: info

    call build_radial_cubic_spline_grid(nodes, 0.0_dp, 1.0_dp, grid, info)
    call require(info == radial_cubic_spline_ok, "grid construction failed")
    call check_manufactured_jet(grid)
    call check_axis_limits(grid)
    call check_mathematica_fixture(grid)
    call check_invalid_inputs(grid)
    write (*, "(a)") "PASS"

contains

    subroutine check_manufactured_jet(valid_grid)
        type(radial_cubic_spline_grid_t), intent(in) :: valid_grid
        integer, parameter :: modes(6) = [0, 1, -2, 3, 4, 5]
        real(dp), parameter :: queries(3) = [0.02_dp, 0.43_dp, 1.0_dp]
        type(axis_regular_harmonic_field_t) :: field
        real(dp) :: samples(size(nodes), size(modes))
        real(dp) :: values(size(modes)), slopes(size(modes))
        real(dp) :: seconds(size(modes))
        real(dp) :: expected, exponent
        integer :: column, query, status

        call build_samples(modes, samples)
        call fit_axis_regular_harmonics(valid_grid, modes, samples, field, &
            status)
        call require(status == axis_regular_harmonic_ok, &
            "manufactured fit failed")
        do query = 1, size(queries)
            call evaluate_axis_regular_harmonics(valid_grid, field, &
                queries(query), values, slopes, seconds, status)
            call require(status == axis_regular_harmonic_ok, &
                "manufactured evaluation failed")
            do column = 1, size(modes)
                exponent = 0.5_dp * real(abs(modes(column)), dp)
                expected = queries(query)**exponent &
                    * quotient(queries(query), column)
                call require(close(values(column), expected), &
                    "manufactured value differs")
                expected = queries(query)**exponent &
                    * (quotient_s(queries(query), column) + exponent &
                    * quotient(queries(query), column) / queries(query))
                call require(close(slopes(column), expected), &
                    "manufactured slope differs")
                expected = queries(query)**exponent &
                    * (quotient_ss(queries(query), column) &
                    + 2.0_dp * exponent &
                    * quotient_s(queries(query), column) / queries(query) &
                    + exponent * (exponent - 1.0_dp) &
                    * quotient(queries(query), column) / queries(query)**2)
                call require(close(seconds(column), expected), &
                    "manufactured second derivative differs")
            end do
        end do
    end subroutine check_manufactured_jet

    subroutine check_axis_limits(valid_grid)
        type(radial_cubic_spline_grid_t), intent(in) :: valid_grid
        integer, parameter :: modes(4) = [0, -2, 4, 5]
        type(axis_regular_harmonic_field_t) :: field
        real(dp) :: samples(size(nodes), size(modes))
        real(dp) :: values(size(modes)), slopes(size(modes))
        real(dp) :: seconds(size(modes))
        integer :: status

        call build_samples(modes, samples)
        call fit_axis_regular_harmonics(valid_grid, modes, samples, field, &
            status)
        call evaluate_axis_regular_harmonics(valid_grid, field, 0.0_dp, &
            values, slopes, seconds, status)
        call require(status == axis_regular_harmonic_ok, &
            "regular axis evaluation failed")
        call require(close(values(1), quotient(0.0_dp, 1)), &
            "m=0 axis value differs")
        call require(close(slopes(1), quotient_s(0.0_dp, 1)), &
            "m=0 axis slope differs")
        call require(close(seconds(1), quotient_ss(0.0_dp, 1)), &
            "m=0 axis second derivative differs")
        call require(close(slopes(2), quotient(0.0_dp, 2)), &
            "absolute m=2 axis slope differs")
        call require(close(seconds(2), 2.0_dp * quotient_s(0.0_dp, 2)), &
            "absolute m=2 axis second derivative differs")
        call require(close(seconds(3), 2.0_dp * quotient(0.0_dp, 3)), &
            "absolute m=4 axis second derivative differs")
        call require(all(values(2:) == 0.0_dp) &
            .and. all(slopes(3:) == 0.0_dp) &
            .and. seconds(4) == 0.0_dp, "zero axis limits differ")
    end subroutine check_axis_limits

    subroutine check_mathematica_fixture(valid_grid)
        type(radial_cubic_spline_grid_t), intent(in) :: valid_grid
        type(axis_regular_harmonic_field_t) :: field
        real(dp) :: samples(size(nodes), 1), values(1), slopes(1), seconds(1)
        integer :: status

        samples(:, 1) = sqrt(nodes) &
            * (1.0_dp + 2.0_dp * nodes - 3.0_dp * nodes**2 + nodes**3)
        call fit_axis_regular_harmonics(valid_grid, [1], samples, field, &
            status)
        call evaluate_axis_regular_harmonics(valid_grid, field, 0.25_dp, &
            values, slopes, seconds, status)
        call require(status == axis_regular_harmonic_ok, &
            "Mathematica fixture evaluation failed")
        call require(close(values(1), 85.0_dp / 128.0_dp), &
            "Mathematica fixture value differs")
        call require(close(slopes(1), 107.0_dp / 64.0_dp), &
            "Mathematica fixture slope differs")
        call require(close(seconds(1), -113.0_dp / 32.0_dp), &
            "Mathematica fixture second derivative differs")
    end subroutine check_mathematica_fixture

    subroutine check_invalid_inputs(valid_grid)
        type(radial_cubic_spline_grid_t), intent(in) :: valid_grid
        type(radial_cubic_spline_grid_t) :: shifted_grid
        type(axis_regular_harmonic_field_t) :: field
        real(dp) :: samples(size(nodes), 2), values(2), slopes(2), seconds(2)
        integer :: status

        call build_samples([0, 1], samples)
        call fit_axis_regular_harmonics(valid_grid, [0], samples, field, &
            status)
        call require(status == axis_regular_harmonic_invalid, &
            "mismatched modes were accepted")
        call fit_axis_regular_harmonics(valid_grid, [-huge(0) - 1, 1], &
            samples, field, status)
        call require(status == axis_regular_harmonic_invalid, &
            "unrepresentable absolute mode was accepted")
        samples(2, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
        call fit_axis_regular_harmonics(valid_grid, [0, 1], samples, field, &
            status)
        call require(status == axis_regular_harmonic_invalid, &
            "nonfinite samples were accepted")
        call build_samples([0, 1], samples)
        call build_radial_cubic_spline_grid(nodes, 0.01_dp, 1.0_dp, &
            shifted_grid, status)
        call fit_axis_regular_harmonics(shifted_grid, [0, 1], samples, field, &
            status)
        call require(status == axis_regular_harmonic_invalid, &
            "non-axis reconstruction domain was accepted")
        call fit_axis_regular_harmonics(valid_grid, [0, 1], samples, field, &
            status)
        call evaluate_axis_regular_harmonics(valid_grid, field, 0.0_dp, &
            values, slopes, seconds, status)
        call require(status == axis_regular_harmonic_invalid, &
            "singular m=1 axis jet was accepted")
        call require(all(values == 0.0_dp) .and. all(slopes == 0.0_dp) &
            .and. all(seconds == 0.0_dp), "failed evaluation changed outputs")
        call evaluate_axis_regular_harmonics(valid_grid, field, -0.01_dp, &
            values, slopes, seconds, status)
        call require(status == axis_regular_harmonic_invalid, &
            "coordinate outside domain was accepted")
        call evaluate_axis_regular_harmonics(valid_grid, field, 0.5_dp, &
            values(:1), slopes, seconds, status)
        call require(status == axis_regular_harmonic_invalid, &
            "wrong output extent was accepted")
    end subroutine check_invalid_inputs

    subroutine build_samples(modes, samples)
        integer, intent(in) :: modes(:)
        real(dp), intent(out) :: samples(:, :)
        integer :: column

        do column = 1, size(modes)
            samples(:, column) = nodes**(0.5_dp * real(abs(modes(column)), dp)) &
                * quotient(nodes, column)
        end do
    end subroutine build_samples

    elemental function quotient(s, column) result(value)
        real(dp), intent(in) :: s
        integer, intent(in) :: column
        real(dp) :: value

        value = 1.0_dp + 0.2_dp * column &
            + (0.5_dp * column - 1.0_dp) * s &
            + (0.25_dp - 0.1_dp * column) * s**2 &
            + 0.05_dp * column * s**3
    end function quotient

    elemental function quotient_s(s, column) result(value)
        real(dp), intent(in) :: s
        integer, intent(in) :: column
        real(dp) :: value

        value = 0.5_dp * column - 1.0_dp &
            + 2.0_dp * (0.25_dp - 0.1_dp * column) * s &
            + 0.15_dp * column * s**2
    end function quotient_s

    elemental function quotient_ss(s, column) result(value)
        real(dp), intent(in) :: s
        integer, intent(in) :: column
        real(dp) :: value

        value = 2.0_dp * (0.25_dp - 0.1_dp * column) &
            + 0.3_dp * column * s
    end function quotient_ss

    function close(actual, expected) result(matches)
        real(dp), intent(in) :: actual, expected
        logical :: matches

        matches = abs(actual - expected) &
            <= 3.0e-11_dp * max(1.0_dp, abs(expected))
    end function close

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        write (error_unit, "(a)") "FAIL: " // message
        error stop 1
    end subroutine require

end program test_axis_regular_harmonic_spline
