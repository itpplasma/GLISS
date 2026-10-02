module variable_spectrum_analysis
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use variable_block_tridiagonal, only: variable_block_tridiagonal_t, &
        variable_pencil_scale
    use variable_generalized_solver, only: validate_variable_pencil, &
        variable_generalized_inertia, variable_generalized_ok
    implicit none
    private

    integer, parameter, public :: variable_spectrum_ok = 0
    integer, parameter, public :: variable_spectrum_invalid = -1
    integer, parameter, public :: variable_spectrum_no_convergence = -2

    type, public :: variable_spectrum_summary_t
        real(dp) :: zero_floor = 0.0_dp
        ! Counts use lambda < -zero_floor and the closed floor band
        ! -zero_floor <= lambda <= zero_floor.
        integer :: negative_count = 0
        integer :: zero_count = 0
        logical :: has_positive = .false.
        ! The first positive cluster lies in [lower, upper).
        real(dp) :: first_positive_lower = 0.0_dp
        real(dp) :: first_positive_upper = 0.0_dp
        integer :: first_positive_cluster_count = 0
    end type variable_spectrum_summary_t

    public :: analyze_variable_spectrum

contains

    ! lowest_only: the caller needs only the lowest eigenvalue. The first
    ! positive cluster is then sought only when no eigenvalue is negative,
    ! and only isolated (base count at the lower end, more at the upper)
    ! rather than refined to the floor width.
    subroutine analyze_variable_spectrum(stiffness, mass, zero_floor, &
            summary, info, lowest_only)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: zero_floor
        type(variable_spectrum_summary_t), intent(out) :: summary
        integer, intent(out) :: info
        logical, intent(in), optional :: lowest_only
        real(dp) :: lower_shift, upper_shift
        integer :: count_above, count_below, dimension
        logical :: isolate

        summary = variable_spectrum_summary_t()
        info = variable_spectrum_invalid
        if (.not. ieee_is_finite(zero_floor) .or. zero_floor <= 0.0_dp) return
        if (zero_floor > 0.125_dp * huge(zero_floor)) return
        ! Every probe below factors a shift of this validated pencil.
        call validate_variable_pencil(stiffness, mass, info)
        if (info /= variable_generalized_ok) then
            info = variable_spectrum_invalid
            return
        end if
        lower_shift = -zero_floor
        upper_shift = zero_floor
        call directed_inertia(stiffness, mass, lower_shift, -1.0_dp, &
            count_below, info)
        if (info /= variable_spectrum_ok) return
        call directed_inertia(stiffness, mass, upper_shift, 1.0_dp, &
            count_above, info)
        if (info /= variable_spectrum_ok .or. count_above < count_below) then
            info = variable_spectrum_invalid
            return
        end if
        dimension = sum(stiffness%widths)
        summary%zero_floor = zero_floor
        summary%negative_count = count_below
        summary%zero_count = count_above - count_below
        if (count_above == dimension) then
            info = variable_spectrum_ok
            return
        end if
        isolate = .false.
        if (present(lowest_only)) isolate = lowest_only
        if (isolate .and. count_below > 0) then
            info = variable_spectrum_ok
            return
        end if
        call bracket_first_positive(stiffness, mass, zero_floor, count_above, &
            summary, info, isolate)
    end subroutine analyze_variable_spectrum

    subroutine bracket_first_positive(stiffness, mass, zero_floor, base_count, &
            summary, info, isolate)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: zero_floor
        integer, intent(in) :: base_count
        type(variable_spectrum_summary_t), intent(inout) :: summary
        integer, intent(out) :: info
        logical, intent(in) :: isolate
        real(dp) :: lower, upper
        integer :: upper_count

        lower = zero_floor
        upper = 2.0_dp * zero_floor
        call expand_positive_bracket(stiffness, mass, base_count, lower, &
            upper, upper_count, info, merge(16.0_dp, 2.0_dp, isolate))
        if (info /= variable_spectrum_ok) return
        if (.not. isolate) then
            call refine_positive_bracket(stiffness, mass, zero_floor, &
                base_count, lower, upper, upper_count, info)
            if (info /= variable_spectrum_ok) return
        end if
        summary%has_positive = .true.
        summary%first_positive_lower = lower
        summary%first_positive_upper = upper
        summary%first_positive_cluster_count = upper_count - base_count
    end subroutine bracket_first_positive

    ! Grow upper by factor until an eigenvalue lies below it; lower is the
    ! last probe that had none.
    subroutine expand_positive_bracket(stiffness, mass, base_count, lower, &
            upper, upper_count, info, factor)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        integer, intent(in) :: base_count
        real(dp), intent(inout) :: lower, upper
        integer, intent(out) :: upper_count, info
        real(dp), intent(in) :: factor
        integer :: iteration

        do iteration = 1, 4096
            call upward_inertia(stiffness, mass, upper, upper_count, info)
            if (info /= variable_spectrum_ok) return
            if (upper_count > base_count) return
            if (upper > huge(upper) / (2.0_dp * factor)) then
                info = variable_spectrum_no_convergence
                return
            end if
            lower = upper
            upper = factor * upper
        end do
        info = variable_spectrum_no_convergence
    end subroutine expand_positive_bracket

    subroutine refine_positive_bracket(stiffness, mass, tolerance, base_count, &
            lower, upper, upper_count, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: base_count
        real(dp), intent(inout) :: lower, upper
        integer, intent(inout) :: upper_count
        integer, intent(out) :: info
        real(dp) :: midpoint
        integer :: count, iteration

        do iteration = 1, 4096
            if (upper - lower <= tolerance) then
                info = variable_spectrum_ok
                return
            end if
            midpoint = lower + 0.5_dp * (upper - lower)
            if (midpoint <= lower .or. midpoint >= upper) then
                info = variable_spectrum_no_convergence
                return
            end if
            call variable_generalized_inertia(stiffness, mass, midpoint, &
                count, info, validated=.true.)
            if (info /= variable_generalized_ok) then
                call resolve_singular_probe(stiffness, mass, midpoint, &
                    base_count, lower, upper, upper_count, info)
                if (info /= variable_spectrum_ok) return
                cycle
            end if
            if (count > base_count) then
                upper = midpoint
                upper_count = count
            else
                lower = midpoint
            end if
        end do
        info = variable_spectrum_no_convergence
    end subroutine refine_positive_bracket

    subroutine upward_inertia(stiffness, mass, shift, count, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(inout) :: shift
        integer, intent(out) :: count, info

        call directed_inertia(stiffness, mass, shift, 1.0_dp, count, info)
    end subroutine upward_inertia

    subroutine directed_inertia(stiffness, mass, shift, direction, count, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(inout) :: shift
        real(dp), intent(in) :: direction
        integer, intent(out) :: count, info
        real(dp) :: origin, delta, left, right
        integer :: attempt, left_count, left_info, right_count, right_info

        origin = shift
        call variable_generalized_inertia(stiffness, mass, shift, count, info, &
            validated=.true.)
        if (info == variable_generalized_ok) then
            info = variable_spectrum_ok
            return
        end if
        shift = nearest(shift, direction)
        call variable_generalized_inertia(stiffness, mass, shift, count, info, &
            validated=.true.)
        if (info == variable_generalized_ok) then
            info = variable_spectrum_ok
            return
        end if
        ! An unresolved radial pivot can affect more than one ULP. Probe
        ! both sides with increasing separation; equal reliable counts
        ! certify the original shift without silently moving the floor
        ! across a physical eigenvalue. Different counts remain unresolved.
        shift = origin
        delta = min(16.0_dp * epsilon(1.0_dp) &
            * max(abs(origin), variable_pencil_scale(stiffness, mass)), &
            0.25_dp * abs(origin))
        do attempt = 1, 12
            left = origin - delta
            right = origin + delta
            call variable_generalized_inertia(stiffness, mass, left, &
                left_count, left_info, validated=.true.)
            call variable_generalized_inertia(stiffness, mass, right, &
                right_count, right_info, validated=.true.)
            if (left_info == variable_generalized_ok &
                .and. right_info == variable_generalized_ok) then
                if (left_count == right_count) then
                    count = left_count
                    info = variable_spectrum_ok
                    return
                end if
                exit
            end if
            if (delta >= 0.25_dp * abs(origin)) exit
            delta = min(16.0_dp * delta, 0.25_dp * abs(origin))
        end do
        count = -1
        info = variable_spectrum_invalid
    end subroutine directed_inertia

    subroutine resolve_singular_probe(stiffness, mass, shift, base_count, &
            lower, upper, upper_count, info)
        type(variable_block_tridiagonal_t), intent(in) :: stiffness, mass
        real(dp), intent(in) :: shift
        integer, intent(in) :: base_count
        real(dp), intent(inout) :: lower, upper
        integer, intent(inout) :: upper_count
        integer, intent(out) :: info
        real(dp) :: left, right
        integer :: left_count, right_count

        left = nearest(shift, -1.0_dp)
        right = nearest(shift, 1.0_dp)
        call variable_generalized_inertia(stiffness, mass, left, left_count, &
            info, validated=.true.)
        if (info /= variable_generalized_ok) then
            info = variable_spectrum_invalid
            return
        end if
        call variable_generalized_inertia(stiffness, mass, right, right_count, &
            info, validated=.true.)
        if (info /= variable_generalized_ok .or. left_count < base_count &
            .or. right_count < left_count) then
            info = variable_spectrum_invalid
            return
        end if
        if (left_count > base_count) then
            upper = left
            upper_count = left_count
        else if (right_count > base_count) then
            lower = left
            upper = right
            upper_count = right_count
        else
            lower = right
        end if
        info = variable_spectrum_ok
    end subroutine resolve_singular_probe

end module variable_spectrum_analysis
