program test_beta_resonance
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    use export_surface_geometry, only: surface_data_t, &
        solve_beta_derivatives_modes, mu0, mercier_ok
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: mercier_surface_terms, mercier_d_terms_t, &
        mercier_gradient_t
    implicit none
    integer, parameter :: n = 16
    type(surface_data_t) :: surface
    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(mercier_d_terms_t) :: terms
    type(mercier_gradient_t) :: gradient
    real(dp) :: pressure_derivative, toroidal_derivative, poloidal_derivative
    real(dp) :: theta(n), zeta(n), rhs(n, n), delta, flux, plus(n, n), minus(n, n), h
    real(dp) :: difference(n, n)
    real(dp), allocatable :: beta(:, :), bt(:, :), bz(:, :)
    integer :: i, j, info, k

    allocate (surface%jacobian(n, n), surface%b_theta(n, n), surface%b_zeta(n, n))
    do i = 1, n
        theta(i) = real(i - 1, dp) / real(n, dp)
        zeta(i) = theta(i)
    end do
    do j = 1, n
        do i = 1, n
            rhs(i, j) = 0.25_dp * cos(2.0_dp * acos(-1.0_dp) * (theta(i) - zeta(j)))
        end do
    end do
    call require_finite(rhs)
    surface%jacobian = 1.0_dp + rhs
    surface%b_theta = 1.0_dp / surface%jacobian
    surface%b_zeta = surface%b_theta
    equilibrium%field_periods = 1
    equilibrium%poloidal_modes = [0, 1]
    equilibrium%toroidal_modes = [0, 1, -1]
    allocate (surface%mod_b(n, n), surface%area_element(n, n))
    surface%mod_b = 1.0_dp
    surface%area_element = 1.0_dp
    call check_omitted_status(1.0_dp / mu0, -1.0_dp)
    ! The primal zero forcing is solvable, but its pressure tangent is not.
    call check_omitted_status(0.0_dp, 0.0_dp)
    call solve(1.0_dp, 1.0_dp)
    if (info == mercier_ok) error stop 'incompatible rational resonance accepted'
    call solve(1.0_dp + epsilon(1.0_dp), 1.0_dp)
    if (info == mercier_ok) error stop 'unresolved denominator accepted'
    ! Analytic beta = A sin(2 pi(theta-zeta))/(2 pi (p-t)).
    ! Both sides of the former clipping threshold must retain that exact map.
    do k = 1, 4
        delta = (-1.0_dp)**k * 10.0_dp**(-8 - k)
        call solve(1.0_dp + delta, 1.0_dp)
        if (info /= mercier_ok) error stop 'resolved near resonance rejected'
        difference = ((1.0_dp + delta) - 1.0_dp) * bt - rhs
        if (finite_max_abs(difference) > 1.0e-12_dp) &
            error stop 'near resonance analytic derivative'
        difference = bt + bz
        if (finite_max_abs(difference) > 1.0e-12_dp * finite_max_abs(bt)) &
            error stop 'opposite angular derivatives'
    end do
    ! Centered directional derivative of the exact inverse denominator on
    ! either side of resonance; the O((h/delta)**2) error is below 2e-6.
    h = 1.0e-7_dp
    do k = -1, 1, 2
        delta = real(k, dp) * 1.0e-4_dp
        call solve(1.0_dp + delta + h, 1.0_dp)
        if (info /= mercier_ok) error stop 'positive tangent sample rejected'
        plus = bt
        call require_finite(plus)
        call solve(1.0_dp + delta - h, 1.0_dp)
        if (info /= mercier_ok) error stop 'negative tangent sample rejected'
        minus = bt
        call require_finite(minus)
        difference = (plus - minus) * delta**2 / (2.0_dp * h) + rhs
        if (finite_max_abs(difference) &
            > 2.0e-6_dp * finite_max_abs(rhs)) error stop 'flux directional derivative'
    end do
    ! Scale invariance: flux rescaling inversely rescales beta derivatives.
    do k = -1, 1
        flux = 10.0_dp**(150 * k)
        call solve(2.0_dp * flux, flux)
        if (info /= mercier_ok) error stop 'finite reference scale rejected'
        difference = flux * bt - rhs
        if (finite_max_abs(difference) > 1.0e-12_dp) error stop 'flux scaling'
    end do
    call solve(0.0_dp, 0.0_dp)
    if (info == mercier_ok) error stop 'zero flux accepted'
    call solve(ieee_value(0.0_dp, ieee_quiet_nan), 1.0_dp)
    if (info == mercier_ok) error stop 'nonfinite flux accepted'
    surface%jacobian = 1.0_dp
    surface%b_theta = 1.0_dp
    surface%b_zeta = 1.0_dp
    call solve(1.0_dp, 1.0_dp)
    if (info /= mercier_ok) error stop 'compatible gauge rejected'
    if (finite_max_abs(beta) /= 0.0_dp) error stop 'nonzero gauge'
    ! The contract is explicitly mean projected: a constant force-balance
    ! defect has zero projected derivative but a nonzero unprojected residual.
    surface%jacobian = 2.0_dp
    surface%b_theta = 0.5_dp
    surface%b_zeta = 0.5_dp
    call solve(1.0_dp, 1.0_dp)
    if (info /= mercier_ok) error stop 'mean projection rejected'
    difference = bt + bz
    if (finite_max_abs(difference) > 1.0e-12_dp) error stop 'mean projection derivative'
    difference = bt + bz - 1.0_dp
    if (finite_max_abs(difference) < 0.99_dp) error stop 'mean residual hidden'
contains
    subroutine solve(p, t)
        real(dp), intent(in) :: p, t
        call solve_beta_derivatives_modes([0, 1], [0, 1, -1], surface, theta, zeta, &
            0.0_dp, -1.0_dp, 1.0_dp / mu0, p, t, beta, bt, bz, info=info)
        if (info == mercier_ok) then
            call require_finite(beta)
            call require_finite(bt)
            call require_finite(bz)
        end if
    end subroutine solve
    subroutine check_omitted_status(pressure, covariant_slope)
        real(dp), intent(in) :: pressure, covariant_slope
        terms = mercier_d_terms_t(7.0_dp, 7.0_dp, 7.0_dp, 7.0_dp, 7.0_dp)
        gradient = mercier_gradient_t(7.0_dp, 7.0_dp)
        pressure_derivative = 7.0_dp
        toroidal_derivative = 7.0_dp
        poloidal_derivative = 7.0_dp
        call mercier_surface_terms(equilibrium, surface, theta, zeta, &
            1.0_dp, 1.0_dp, 0.0_dp, covariant_slope, 1.0_dp, 0.0_dp, &
            1.0_dp, 0.0_dp, pressure, 1.0_dp, terms, gradient, &
            pressure_derivative, beta, toroidal_derivative, poloidal_derivative)
        if (ieee_is_finite(terms%shear)) error stop 'failed shear retained'
        if (ieee_is_finite(terms%current)) error stop 'failed current retained'
        if (ieee_is_finite(terms%well)) error stop 'failed well retained'
        if (ieee_is_finite(terms%geodesic)) error stop 'failed geodesic retained'
        if (ieee_is_finite(terms%mercier)) error stop 'failed Mercier retained'
        if (ieee_is_finite(gradient%d_iota_slope)) error stop 'failed iota gradient retained'
        if (ieee_is_finite(gradient%d_pressure_slope)) &
            error stop 'failed pressure gradient retained'
        if (ieee_is_finite(pressure_derivative)) error stop 'failed pressure derivative retained'
        if (ieee_is_finite(toroidal_derivative)) error stop 'failed toroidal derivative retained'
        if (ieee_is_finite(poloidal_derivative)) error stop 'failed poloidal derivative retained'
        if (any(ieee_is_finite(beta))) error stop 'failed beta retained'
    end subroutine check_omitted_status

    subroutine require_finite(values)
        real(dp), intent(in) :: values(:, :)
        if (.not. all(ieee_is_finite(values))) error stop 'nonfinite resonance evidence'
    end subroutine require_finite

    function finite_max_abs(values) result(value)
        real(dp), intent(in) :: values(:, :)
        real(dp) :: value
        call require_finite(values)
        value = maxval(abs(values))
    end function finite_max_abs
end program test_beta_resonance
