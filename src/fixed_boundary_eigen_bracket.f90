module fixed_boundary_eigen_bracket
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fixed_boundary_solver_controls, only: fixed_boundary_solver_controls_t
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
            interval, info, controls)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: zero_floor
        real(dp), intent(out) :: shift, interval
        integer, intent(out) :: info
        type(fixed_boundary_solver_controls_t), intent(in), optional :: controls
        type(fixed_boundary_solver_controls_t) :: stopping
        real(dp) :: lower, upper, middle
        integer :: count, iteration

        info = fixed_boundary_bracket_error
        stopping = fixed_boundary_solver_controls_t()
        if (present(controls)) stopping = controls
        lower = -2.0_dp * zero_floor
        do iteration = 1, stopping%bracket_iteration_limit
            call variable_generalized_inertia(stiffness, mass, lower, count, &
                info)
            if (info /= variable_generalized_ok) then
                lower = lower * (1.0_dp + 1.0e-8_dp)
                cycle
            end if
            if (count == 0) exit
            lower = 2.0_dp * lower
        end do
        if (iteration > stopping%bracket_iteration_limit) then
            info = fixed_boundary_bracket_error
            return
        end if
        upper = -zero_floor
        do iteration = 1, stopping%bracket_iteration_limit
            middle = 0.5_dp * (lower + upper)
            if (upper - lower <= stopping%negative_bracket_relative &
                * abs(middle) + stopping%negative_bracket_floor &
                * zero_floor) exit
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
        if (iteration > stopping%bracket_iteration_limit) then
            info = fixed_boundary_bracket_error
            return
        end if
        shift = lower
        interval = upper - lower
        info = fixed_boundary_bracket_ok
    end subroutine bracket_lowest_negative

end module fixed_boundary_eigen_bracket
