module fixed_boundary_eigen_bracket
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fixed_boundary_solver_controls, only: fixed_boundary_solver_controls_t
    use variable_block_tridiagonal, only: variable_block_tridiagonal_t, &
        variable_pencil_scale
    use variable_generalized_solver, only: validate_variable_pencil, &
        variable_generalized_inertia, variable_generalized_ok
    implicit none
    private

    integer, parameter, public :: fixed_boundary_bracket_ok = 0
    integer, parameter, public :: fixed_boundary_bracket_error = -1
    integer, parameter, public :: fixed_boundary_bracket_expansion_error = -12
    integer, parameter, public :: fixed_boundary_bracket_probe_error = -13
    integer, parameter, public :: fixed_boundary_bracket_refinement_error = -14

    public :: bracket_lowest_negative, bounded_inertia_probe
    public :: bracket_lowest_positive, certify_lowest_bracket

    real(dp), parameter :: expansion_factor = 16.0_dp

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
        call validate_variable_pencil(stiffness, mass, info)
        if (info /= variable_generalized_ok) then
            info = fixed_boundary_bracket_error
            return
        end if
        ! Expand geometrically by a factor of 16 per probe; the last probe
        ! with negative eigenvalues below it is the upper end of the bracket.
        upper = -zero_floor
        lower = -2.0_dp * zero_floor
        do iteration = 1, stopping%bracket_iteration_limit
            call variable_generalized_inertia(stiffness, mass, lower, count, &
                info, validated=.true.)
            if (info == variable_generalized_ok) then
                if (count == 0) exit
                upper = lower
            end if
            lower = expansion_factor * lower
        end do
        if (iteration > stopping%bracket_iteration_limit) then
            info = fixed_boundary_bracket_expansion_error
            return
        end if
        do iteration = 1, stopping%bracket_iteration_limit
            middle = 0.5_dp * (lower + upper)
            if (upper - lower <= stopping%negative_bracket_relative &
                * abs(middle) + stopping%negative_bracket_floor &
                * zero_floor) exit
            call bounded_inertia_probe(stiffness, mass, lower, upper, &
                middle, count, info, validated=.true.)
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

    ! Bisect the first-positive inertia bracket with the same stopping rule
    ! as the negative branch, so the reported interval is a refined
    ! inertia enclosure rather than the coarse spectrum-summary bracket.
    subroutine bracket_lowest_positive(stiffness, mass, lower, upper, shift, &
            interval, info, controls)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: lower, upper
        real(dp), intent(out) :: shift, interval
        integer, intent(out) :: info
        type(fixed_boundary_solver_controls_t), intent(in), optional :: controls
        type(fixed_boundary_solver_controls_t) :: stopping
        real(dp) :: below, above, middle
        integer :: base_count, count, iteration

        stopping = fixed_boundary_solver_controls_t()
        if (present(controls)) stopping = controls
        shift = lower
        interval = upper - lower
        below = lower
        above = upper
        call validate_variable_pencil(stiffness, mass, info)
        if (info /= variable_generalized_ok) then
            info = fixed_boundary_bracket_error
            return
        end if
        call bounded_inertia_probe(stiffness, mass, below, above, below, &
            base_count, info, validated=.true.)
        if (info /= fixed_boundary_bracket_ok) return
        do iteration = 1, stopping%bracket_iteration_limit
            middle = 0.5_dp * (below + above)
            if (above - below <= stopping%negative_bracket_relative &
                * abs(middle)) exit
            call bounded_inertia_probe(stiffness, mass, below, above, &
                middle, count, info, validated=.true.)
            if (info /= fixed_boundary_bracket_ok) return
            if (count == base_count) then
                below = middle
            else
                above = middle
            end if
        end do
        if (iteration > stopping%bracket_iteration_limit) then
            info = fixed_boundary_bracket_refinement_error
            return
        end if
        shift = below
        interval = above - below
        info = fixed_boundary_bracket_ok
    end subroutine bracket_lowest_positive

    subroutine bounded_inertia_probe(stiffness, mass, lower, upper, probe, &
            count, info, validated)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: lower, upper
        real(dp), intent(inout) :: probe
        integer, intent(out) :: count, info
        logical, intent(in), optional :: validated
        real(dp) :: candidate, delta, origin, scale
        integer :: attempt

        origin = probe
        scale = max(abs(origin), variable_pencil_scale(stiffness, mass))
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
                count, info, validated)
            if (info == variable_generalized_ok) then
                probe = candidate
                info = fixed_boundary_bracket_ok
                return
            end if
        end do
        info = fixed_boundary_bracket_probe_error
    end subroutine bounded_inertia_probe

    ! Inertia bracket [lower, upper] = [eigenvalue - margin, eigenvalue +
    ! margin] of the eigenvalue above base_count others: base_count
    ! eigenvalues lie below lower and at least one more below upper. Two
    ! factorizations certify an inverse-iteration eigenvalue whose residual
    ! bound is below margin, instead of bisecting to that width.
    subroutine certify_lowest_bracket(stiffness, mass, base_count, &
            eigenvalue, margin, lower, upper, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        integer, intent(in) :: base_count
        real(dp), intent(in) :: eigenvalue, margin
        real(dp), intent(out) :: lower, upper
        integer, intent(out) :: info
        integer :: count

        info = fixed_boundary_bracket_error
        lower = eigenvalue - margin
        upper = eigenvalue + margin
        if (.not. (margin > 0.0_dp) .or. .not. (upper > lower)) return
        call variable_generalized_inertia(stiffness, mass, lower, count, &
            info, validated=.true.)
        if (info /= variable_generalized_ok .or. count /= base_count) then
            info = fixed_boundary_bracket_error
            return
        end if
        call variable_generalized_inertia(stiffness, mass, upper, count, &
            info, validated=.true.)
        if (info /= variable_generalized_ok .or. count <= base_count) then
            info = fixed_boundary_bracket_error
            return
        end if
        info = fixed_boundary_bracket_ok
    end subroutine certify_lowest_bracket

end module fixed_boundary_eigen_bracket
