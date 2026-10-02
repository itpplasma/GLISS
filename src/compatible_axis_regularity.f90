module compatible_axis_regularity
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fourier_phase_kind, only: phase_sine
    use primitive_equilibrium_spline, only: primitive_equilibrium_spline_t
    use radial_cubic_spline, only: evaluate_radial_cubic_spline_field, &
        radial_cubic_spline_ok
    use radial_feec_complex, only: radial_feec_complex_t
    implicit none
    private

    ! A displacement that is smooth in Cartesian coordinates near the
    ! magnetic axis has, in every regular flux chart
    ! x = x_axis(zeta) + sqrt(s) (c(zeta) cos(theta) + d(zeta) sin(theta))
    ! + O(s), a |m|=1 part
    !     xi^s = 2 sqrt(s) (u cos(theta) + v sin(theta)),
    !     xi^theta = (v cos(theta) - u sin(theta)) / sqrt(s)
    ! (theta in radians), independent of the axis shape c, d. The leading
    ! coefficients of xi^s ~ a s^(1/2) and eta ~ Phi' xi^theta ~ b s^(-1/2)
    ! are therefore tied: the s^(-1/2) part of the compression
    ! d(xi^s)/ds + d(eta)/dtheta / Phi' vanishes. The FEEC space stores
    ! a in the first retained H1 function and b in the first L2 function as
    ! independent coefficients; their untied combination has a
    ! logarithmically divergent compression energy and lies outside the
    ! energy space, so the discrete space was not conforming at the axis.
    ! The conforming subspace eliminates the first eta coefficient of each
    ! |m|=1 trial: eta_1 = kappa xi_1 with
    !     kappa = -sigma N'(0) Phi'(0) / (4 pi),
    ! sigma = +1 when the normal component is the cosine family and -1 when
    ! it is the sine family, N'(0) the slope of the first retained H1
    ! function and Phi'(0) the kernel toroidal-flux slope sqrt(g) B^zeta at
    ! s = 0 (theta has period one). For |m| >= 2 the untied leading
    ! coefficients have finite energy; they need no constraint.

    integer, parameter, public :: axis_regularity_ok = 0
    integer, parameter, public :: axis_regularity_invalid = -1

    type, public :: axis_tie_t
        integer :: full_unknowns = 0
        integer :: reduced_unknowns = 0
        integer :: eliminated_count = 0
        integer, allocatable :: target(:)
        real(dp), allocatable :: factor(:)
        logical, allocatable :: eliminated(:)
    end type axis_tie_t

    real(dp), parameter :: pi = acos(-1.0_dp)

    public :: axis_flux_slope
    public :: axis_tie_factor
    public :: build_axis_tie
    public :: build_trial_axis_tie
    public :: requires_axis_tie
    public :: tie_local_map

contains

    pure logical function requires_axis_tie(mode_m, stored_power) &
            result(required)
        integer, intent(in) :: mode_m
        real(dp), intent(in) :: stored_power

        ! Only the regular-axis representation s^(-1/2) of |m|=1 carries
        ! the tied leading coefficients.
        required = abs(mode_m) == 1 .and. stored_power == 0.5_dp
    end function requires_axis_tie

    subroutine axis_flux_slope(spline, slope, info)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        real(dp), intent(out) :: slope
        integer, intent(out) :: info
        real(dp) :: values(3), slopes(3), seconds(3)
        integer :: local_info

        ! The kernel field sqrt(g) B^zeta equals -dPhi/ds of the profile
        ! spline (cartesian_primitive_geometry), evaluated at the axis.
        slope = 0.0_dp
        info = axis_regularity_invalid
        call evaluate_radial_cubic_spline_field(spline%radial_grid, &
            spline%profiles, 0.0_dp, values, slopes, seconds, local_info)
        if (local_info /= radial_cubic_spline_ok) return
        if (.not. ieee_is_finite(slopes(1)) .or. slopes(1) == 0.0_dp) return
        slope = -slopes(1)
        info = axis_regularity_ok
    end subroutine axis_flux_slope

    pure function axis_tie_factor(complex, parity, flux_slope) result(kappa)
        type(radial_feec_complex_t), intent(in) :: complex
        integer, intent(in) :: parity
        real(dp), intent(in) :: flux_slope
        real(dp) :: kappa, sigma, slope

        ! The first retained H1 function of the open knot vector is the
        ! only one with a nonzero slope at s = 0, and the first L2
        ! function is the only one with a nonzero value (equal to one).
        slope = real(complex%h1_degree, dp) &
            / (complex%h1_knots(complex%h1_degree + 2) - complex%h1_knots(2))
        sigma = 1.0_dp
        if (parity == phase_sine) sigma = -1.0_dp
        kappa = -sigma * slope * flux_slope / (4.0_dp * pi)
    end function axis_tie_factor

    subroutine build_axis_tie(full_unknowns, normal_index, eta_index, &
            kappa, tie, info)
        integer, intent(in) :: full_unknowns
        integer, intent(in) :: normal_index(:), eta_index(:)
        real(dp), intent(in) :: kappa(:)
        type(axis_tie_t), intent(out) :: tie
        integer, intent(out) :: info
        logical, allocatable :: eliminated(:)
        integer :: full, pair, reduced

        info = axis_regularity_invalid
        if (full_unknowns < 1) return
        if (size(eta_index) /= size(normal_index)) return
        if (size(kappa) /= size(normal_index)) return
        if (any(normal_index < 1) .or. any(normal_index > full_unknowns)) &
            return
        if (any(eta_index < 1) .or. any(eta_index > full_unknowns)) return
        if (.not. all(ieee_is_finite(kappa))) return
        allocate (eliminated(full_unknowns), source=.false.)
        do pair = 1, size(eta_index)
            if (eliminated(eta_index(pair))) return
            eliminated(eta_index(pair)) = .true.
        end do
        if (any(eliminated(normal_index))) return
        call move_alloc(eliminated, tie%eliminated)
        allocate (tie%target(full_unknowns), source=0)
        allocate (tie%factor(full_unknowns), source=1.0_dp)
        reduced = 0
        do full = 1, full_unknowns
            if (tie%eliminated(full)) cycle
            reduced = reduced + 1
            tie%target(full) = reduced
        end do
        do pair = 1, size(eta_index)
            tie%target(eta_index(pair)) = tie%target(normal_index(pair))
            tie%factor(eta_index(pair)) = kappa(pair)
        end do
        tie%full_unknowns = full_unknowns
        tie%reduced_unknowns = reduced
        tie%eliminated_count = size(eta_index)
        info = axis_regularity_ok
    end subroutine build_axis_tie

    subroutine build_trial_axis_tie(spline, complex, mode_m, parity, &
            stored_power, normal_rank, eta_rank, full_unknowns, conforming, &
            tie, info)
        type(primitive_equilibrium_spline_t), intent(in) :: spline
        type(radial_feec_complex_t), intent(in) :: complex
        integer, intent(in) :: mode_m(:), parity(:)
        real(dp), intent(in) :: stored_power(:)
        integer, intent(in) :: normal_rank(:), eta_rank(:), full_unknowns
        logical, intent(in) :: conforming
        type(axis_tie_t), intent(out) :: tie
        integer, intent(out) :: info
        integer, allocatable :: normal_index(:), eta_index(:)
        real(dp), allocatable :: kappa(:)
        real(dp) :: flux_slope
        integer :: pairs, trial

        ! Global layout of the compatible problems: the normal unknown of
        ! H1 function b and trial t is (b - 1) * N + rank_n(t), and the first
        ! eta unknown of trial t follows all normal unknowns at
        ! h1_dofs * N + rank_eta(t), with N active normal trials.
        info = axis_regularity_invalid
        pairs = 0
        if (conforming) then
            do trial = 1, size(mode_m)
                if (requires_axis_tie(mode_m(trial), stored_power(trial)) &
                    .and. normal_rank(trial) > 0 .and. eta_rank(trial) > 0) &
                    pairs = pairs + 1
            end do
        end if
        allocate (normal_index(pairs), eta_index(pairs), kappa(pairs))
        if (pairs > 0) then
            if (.not. complex%left_trace .or. complex%h1_dofs < 1) return
            if (complex%h1_basis_index(1) /= 2) return
            call axis_flux_slope(spline, flux_slope, info)
            if (info /= axis_regularity_ok) return
            pairs = 0
            do trial = 1, size(mode_m)
                if (.not. requires_axis_tie(mode_m(trial), &
                    stored_power(trial))) cycle
                if (normal_rank(trial) < 1 .or. eta_rank(trial) < 1) cycle
                pairs = pairs + 1
                normal_index(pairs) = normal_rank(trial)
                eta_index(pairs) = complex%h1_dofs * count(normal_rank > 0) &
                    + eta_rank(trial)
                kappa(pairs) = axis_tie_factor(complex, parity(trial), &
                    flux_slope)
            end do
        end if
        call build_axis_tie(full_unknowns, normal_index, eta_index, kappa, &
            tie, info)
    end subroutine build_trial_axis_tie

    pure subroutine tie_local_map(tie, map, reduced_map, scale)
        type(axis_tie_t), intent(in) :: tie
        integer, intent(in) :: map(:)
        integer, intent(out) :: reduced_map(:)
        real(dp), intent(out) :: scale(:)
        integer :: column

        ! A local column with full index g contributes factor(g) times its
        ! basis function to reduced unknown target(g). The reduced local
        ! matrix is diag(scale) A diag(scale), scattered through the
        ! reduced map; repeated targets add, which forms T^T A T.
        do column = 1, size(map)
            reduced_map(column) = 0
            scale(column) = 0.0_dp
            if (map(column) == 0) cycle
            reduced_map(column) = tie%target(map(column))
            scale(column) = tie%factor(map(column))
        end do
    end subroutine tie_local_map

end module compatible_axis_regularity
