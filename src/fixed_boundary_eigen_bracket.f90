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
    integer, parameter, public :: fixed_boundary_bracket_expansion_error = -12
    integer, parameter, public :: fixed_boundary_bracket_probe_error = -13
    integer, parameter, public :: fixed_boundary_bracket_refinement_error = -14

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
                lower = 2.0_dp * lower
                cycle
            end if
            if (count == 0) exit
            lower = 2.0_dp * lower
        end do
        if (iteration > stopping%bracket_iteration_limit) then
            info = fixed_boundary_bracket_expansion_error
            return
        end if
        upper = -zero_floor
        do iteration = 1, stopping%bracket_iteration_limit
            middle = 0.5_dp * (lower + upper)
            if (upper - lower <= stopping%negative_bracket_relative &
                * abs(middle) + stopping%negative_bracket_floor &
                * zero_floor) exit
            call bounded_inertia_probe(stiffness, mass, lower, upper, &
                middle, count, info)
            if (info /= fixed_boundary_bracket_ok) return
            if (count == 0) then
                lower = middle
            else
                upper = middle
            end if
        end do
        if (iteration > stopping%bracket_iteration_limit) then
            info = fixed_boundary_bracket_refinement_error
            return
        end if
        shift = lower
        interval = upper - lower
        info = fixed_boundary_bracket_ok
    end subroutine bracket_lowest_negative

    subroutine bounded_inertia_probe(stiffness, mass, lower, upper, probe, &
            count, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: lower, upper
        real(dp), intent(inout) :: probe
        integer, intent(out) :: count, info
        real(dp) :: candidate, delta, origin, scale
        integer :: attempt

        origin = probe
        scale = max(1.0_dp, abs(origin))
        delta = min(16.0_dp * epsilon(1.0_dp) * scale, &
            0.25_dp * (upper - lower))
        do attempt = 0, 15
            if (attempt == 0) then
                candidate = origin
            else if (modulo(attempt, 2) == 1) then
                candidate = max(lower, origin - delta)
            else
                candidate = min(upper, origin + delta)
                delta = min(16.0_dp * delta, 0.25_dp * (upper - lower))
            end if
            call variable_generalized_inertia(stiffness, mass, candidate, &
                count, info)
            if (info == variable_generalized_ok) then
                probe = candidate
                info = fixed_boundary_bracket_ok
                return
            end if
        end do
        info = fixed_boundary_bracket_probe_error
    end subroutine bounded_inertia_probe

end module fixed_boundary_eigen_bracket
