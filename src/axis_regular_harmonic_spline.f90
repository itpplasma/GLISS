module axis_regular_harmonic_spline
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use radial_cubic_spline, only: evaluate_radial_cubic_spline_field, &
        fit_radial_cubic_spline_field, radial_cubic_spline_allocation_error, &
        radial_cubic_spline_field_t, radial_cubic_spline_grid_t, &
        radial_cubic_spline_ok
    implicit none
    private

    integer, parameter, public :: axis_regular_harmonic_ok = 0
    integer, parameter, public :: axis_regular_harmonic_invalid = -1
    integer, parameter, public :: axis_regular_harmonic_allocation_error = -2

    type, public :: axis_regular_harmonic_field_t
        integer, allocatable :: poloidal_modes(:)
        type(radial_cubic_spline_field_t) :: quotient
    end type axis_regular_harmonic_field_t

    public :: evaluate_axis_regular_harmonics
    public :: fit_axis_regular_harmonics

contains

    subroutine fit_axis_regular_harmonics(grid, poloidal_modes, values, &
            field, info)
        type(radial_cubic_spline_grid_t), intent(in) :: grid
        integer, intent(in) :: poloidal_modes(:)
        real(dp), intent(in) :: values(:, :)
        type(axis_regular_harmonic_field_t), intent(out) :: field
        integer, intent(out) :: info
        real(dp), allocatable :: quotient_values(:, :)
        real(dp) :: exponent
        integer :: allocation_status, column, spline_info

        info = axis_regular_harmonic_invalid
        if (.not. grid_is_axis_half_grid(grid)) return
        if (size(values, 1) /= size(grid%nodes)) return
        if (size(values, 2) < 1) return
        if (size(poloidal_modes) /= size(values, 2)) return
        if (.not. all(ieee_is_finite(values))) return
        if (.not. modes_are_valid(poloidal_modes)) return
        allocate (field%poloidal_modes(size(poloidal_modes)), &
            quotient_values(size(values, 1), size(values, 2)), &
            stat=allocation_status)
        if (allocation_status /= 0) then
            info = axis_regular_harmonic_allocation_error
            return
        end if
        field%poloidal_modes = poloidal_modes
        do column = 1, size(values, 2)
            exponent = 0.5_dp * real(abs(poloidal_modes(column)), dp)
            quotient_values(:, column) = values(:, column) &
                / grid%nodes**exponent
        end do
        if (.not. all(ieee_is_finite(quotient_values))) return
        call fit_radial_cubic_spline_field(grid, quotient_values, &
            field%quotient, spline_info)
        if (spline_info == radial_cubic_spline_allocation_error) then
            info = axis_regular_harmonic_allocation_error
            return
        end if
        if (spline_info /= radial_cubic_spline_ok) return
        info = axis_regular_harmonic_ok
    end subroutine fit_axis_regular_harmonics

    subroutine evaluate_axis_regular_harmonics(grid, field, coordinate, &
            values, derivatives, second_derivatives, info)
        type(radial_cubic_spline_grid_t), intent(in) :: grid
        type(axis_regular_harmonic_field_t), intent(in) :: field
        real(dp), intent(in) :: coordinate
        real(dp), intent(out) :: values(:), derivatives(:)
        real(dp), intent(out) :: second_derivatives(:)
        integer, intent(out) :: info
        real(dp) :: quotient(size(values)), quotient_derivative(size(values))
        real(dp) :: quotient_second(size(values))
        integer :: spline_info

        values = 0.0_dp
        derivatives = 0.0_dp
        second_derivatives = 0.0_dp
        info = axis_regular_harmonic_invalid
        if (.not. grid_is_axis_half_grid(grid)) return
        if (.not. field_is_valid(field)) return
        if (size(values) /= size(field%poloidal_modes)) return
        if (size(derivatives) /= size(field%poloidal_modes)) return
        if (size(second_derivatives) /= size(field%poloidal_modes)) return
        if (.not. ieee_is_finite(coordinate)) return
        if (coordinate < grid%domain_min .or. &
            coordinate > grid%domain_max) return
        if (coordinate == 0.0_dp) then
            if (axis_jet_is_singular(field%poloidal_modes)) return
        end if
        call evaluate_radial_cubic_spline_field(grid, field%quotient, &
            coordinate, quotient, quotient_derivative, quotient_second, &
            spline_info)
        if (spline_info /= radial_cubic_spline_ok) return
        if (coordinate == 0.0_dp) then
            call evaluate_axis_limits(field%poloidal_modes, quotient, &
                quotient_derivative, quotient_second, values, derivatives, &
                second_derivatives)
        else
            call apply_axis_factors(field%poloidal_modes, coordinate, &
                quotient, quotient_derivative, quotient_second, values, &
                derivatives, second_derivatives)
        end if
        if (.not. all(ieee_is_finite(values)) .or. &
            .not. all(ieee_is_finite(derivatives)) .or. &
            .not. all(ieee_is_finite(second_derivatives))) then
            values = 0.0_dp
            derivatives = 0.0_dp
            second_derivatives = 0.0_dp
            return
        end if
        info = axis_regular_harmonic_ok
    end subroutine evaluate_axis_regular_harmonics

    subroutine apply_axis_factors(poloidal_modes, coordinate, quotient, &
            quotient_derivative, quotient_second, values, derivatives, &
            second_derivatives)
        integer, intent(in) :: poloidal_modes(:)
        real(dp), intent(in) :: coordinate
        real(dp), intent(in) :: quotient(:), quotient_derivative(:)
        real(dp), intent(in) :: quotient_second(:)
        real(dp), intent(out) :: values(:), derivatives(:)
        real(dp), intent(out) :: second_derivatives(:)
        real(dp) :: exponent, factor
        integer :: column

        do column = 1, size(poloidal_modes)
            exponent = 0.5_dp * real(abs(poloidal_modes(column)), dp)
            factor = coordinate**exponent
            values(column) = factor * quotient(column)
            derivatives(column) = factor * (quotient_derivative(column) &
                + exponent * quotient(column) / coordinate)
            second_derivatives(column) = factor * (quotient_second(column) &
                + 2.0_dp * exponent * quotient_derivative(column) / coordinate &
                + exponent * (exponent - 1.0_dp) * quotient(column) &
                / coordinate**2)
        end do
    end subroutine apply_axis_factors

    subroutine evaluate_axis_limits(poloidal_modes, quotient, &
            quotient_derivative, quotient_second, values, derivatives, &
            second_derivatives)
        integer, intent(in) :: poloidal_modes(:)
        real(dp), intent(in) :: quotient(:), quotient_derivative(:)
        real(dp), intent(in) :: quotient_second(:)
        real(dp), intent(out) :: values(:), derivatives(:)
        real(dp), intent(out) :: second_derivatives(:)
        integer :: column, mode

        values = 0.0_dp
        derivatives = 0.0_dp
        second_derivatives = 0.0_dp
        do column = 1, size(poloidal_modes)
            mode = abs(poloidal_modes(column))
            select case (mode)
            case (0)
                values(column) = quotient(column)
                derivatives(column) = quotient_derivative(column)
                second_derivatives(column) = quotient_second(column)
            case (2)
                derivatives(column) = quotient(column)
                second_derivatives(column) = 2.0_dp &
                    * quotient_derivative(column)
            case (4)
                second_derivatives(column) = 2.0_dp * quotient(column)
            end select
        end do
    end subroutine evaluate_axis_limits

    function grid_is_axis_half_grid(grid) result(valid)
        type(radial_cubic_spline_grid_t), intent(in) :: grid
        logical :: valid

        valid = allocated(grid%nodes)
        if (.not. valid) return
        valid = grid%domain_min == 0.0_dp .and. size(grid%nodes) >= 4
        if (.not. valid) return
        valid = all(grid%nodes > 0.0_dp)
    end function grid_is_axis_half_grid

    function modes_are_valid(poloidal_modes) result(valid)
        integer, intent(in) :: poloidal_modes(:)
        logical :: valid

        valid = size(poloidal_modes) > 0
        if (.not. valid) return
        valid = all(poloidal_modes >= -huge(poloidal_modes))
    end function modes_are_valid

    function axis_jet_is_singular(poloidal_modes) result(singular)
        integer, intent(in) :: poloidal_modes(:)
        logical :: singular

        singular = any(abs(poloidal_modes) == 1) &
            .or. any(abs(poloidal_modes) == 3)
    end function axis_jet_is_singular

    function field_is_valid(field) result(valid)
        type(axis_regular_harmonic_field_t), intent(in) :: field
        logical :: valid

        valid = allocated(field%poloidal_modes)
        if (.not. valid) return
        valid = modes_are_valid(field%poloidal_modes)
        if (.not. valid) return
        valid = allocated(field%quotient%values)
        if (.not. valid) return
        valid = size(field%quotient%values, 2) == size(field%poloidal_modes)
    end function field_is_valid

end module axis_regular_harmonic_spline
