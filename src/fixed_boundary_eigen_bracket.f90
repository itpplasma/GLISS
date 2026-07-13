module fixed_boundary_eigen_bracket
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use variable_block_tridiagonal, only: variable_block_tridiagonal_t
    use variable_generalized_solver, only: variable_generalized_inertia, &
        variable_generalized_ok
    implicit none
    private

    integer, parameter, public :: fixed_boundary_bracket_ok = 0
    integer, parameter, public :: fixed_boundary_bracket_error = -1

    public :: bracket_lowest_negative

contains

    subroutine bracket_lowest_negative(stiffness, mass, zero_floor, shift, &
            interval, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: zero_floor
        real(dp), intent(out) :: shift, interval
        integer, intent(out) :: info
        real(dp) :: lower, upper, middle
        integer :: count, iteration

        info = fixed_boundary_bracket_error
        lower = -2.0_dp * zero_floor
        do iteration = 1, 200
            call variable_generalized_inertia(stiffness, mass, lower, count, &
                info)
            if (info /= variable_generalized_ok) then
                lower = lower * (1.0_dp + 1.0e-8_dp)
                cycle
            end if
            if (count == 0) exit
            lower = 2.0_dp * lower
        end do
        if (iteration > 200) then
            info = fixed_boundary_bracket_error
            return
        end if
        upper = -zero_floor
        do iteration = 1, 200
            middle = 0.5_dp * (lower + upper)
            if (upper - lower <= 1.0e-9_dp * abs(middle) &
                + 1.0e-3_dp * zero_floor) exit
            call variable_generalized_inertia(stiffness, mass, middle, count, &
                info)
            if (info /= variable_generalized_ok) then
                middle = middle * (1.0_dp + 1.0e-8_dp)
                cycle
            end if
            if (count == 0) then
                lower = middle
            else
                upper = middle
            end if
        end do
        if (iteration > 200) then
            info = fixed_boundary_bracket_error
            return
        end if
        shift = lower
        interval = upper - lower
        info = fixed_boundary_bracket_ok
    end subroutine bracket_lowest_negative

end module fixed_boundary_eigen_bracket
