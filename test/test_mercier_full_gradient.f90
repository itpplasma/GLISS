program test_mercier_full_gradient
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use gvec_cas3d_types, only: gvec_cas3d_equilibrium_t
    use mercier_diagnostic, only: mercier_surface_terms, surface_data_t, &
        mercier_d_terms_t, mercier_gradient_t
    implicit none

    real(dp), parameter :: covariant_theta = 0.7_dp
    real(dp), parameter :: covariant_zeta = -1.3_dp
    real(dp), parameter :: covariant_theta_slope = 0.5_dp
    real(dp), parameter :: covariant_zeta_slope = -0.4_dp
    real(dp), parameter :: flux_slope = 1.2_dp
    real(dp), parameter :: flux_curvature = 0.35_dp
    real(dp), parameter :: volume_slope = 2.0_dp
    real(dp), parameter :: volume_curvature = 0.6_dp
    real(dp), parameter :: iota_slope = 0.7_dp
    real(dp), parameter :: pressure_slope = -0.9_dp
    real(dp), parameter :: tol = 1.0e-7_dp

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(surface_data_t) :: surface
    real(dp), allocatable :: theta(:), zeta(:)
    real(dp), allocatable :: beta_values(:, :)
    type(mercier_d_terms_t) :: d_terms
    type(mercier_gradient_t) :: grad
    real(dp) :: analytic, explicit_only, implicit_part
    real(dp) :: fd, rel_err, implicit_fraction, h

    call build_fixture(equilibrium, surface, theta, zeta)

    call mercier_surface_terms(equilibrium, surface, theta, zeta, &
        covariant_theta, covariant_zeta, covariant_theta_slope, &
        covariant_zeta_slope, flux_slope, flux_curvature, volume_slope, &
        volume_curvature, pressure_slope, iota_slope, d_terms, grad, &
        analytic, beta_values)
    explicit_only = grad%d_pressure_slope
    implicit_part = analytic - explicit_only

    ! d_mercier is a quadratic polynomial in pressure_slope (beta is linear
    ! in it; the well and geodesic terms are at most quadratic), so central
    ! differencing carries no truncation error. Its only error is rounding,
    ! which shrinks as the step grows because the mu0-small sensitivity sits
    ! on an O(1) d_mercier; a unit step keeps cancellation below tolerance.
    h = max(1.0_dp, abs(pressure_slope))
    fd = (d_mercier_of(pressure_slope + h) &
        - d_mercier_of(pressure_slope - h)) / (2.0_dp * h)
    rel_err = abs(fd - analytic) / max(abs(analytic), tiny(1.0_dp))
    implicit_fraction = abs(fd - explicit_only) / max(abs(analytic), &
        tiny(1.0_dp))

    write (*, "(a, es24.16)") "d_mercier               = ", d_terms%mercier
    write (*, "(a, es24.16)") "analytic full gradient  = ", analytic
    write (*, "(a, es24.16)") "explicit-only gradient  = ", explicit_only
    write (*, "(a, es24.16)") "beta-implicit part      = ", implicit_part
    write (*, "(a, es24.16)") "finite difference       = ", fd
    write (*, "(a, es24.16)") "relative error (full)   = ", rel_err
    write (*, "(a, es24.16)") "fd vs explicit-only      = ", implicit_fraction

    call require(rel_err <= tol, &
        "full pressure_slope gradient disagrees with finite difference")
    call require(implicit_fraction >= 100.0_dp * tol, &
        "finite difference must resolve the beta-implicit contribution")
    write (*, "(a)") "PASS"

contains

    function d_mercier_of(p) result(value)
        real(dp), intent(in) :: p
        real(dp) :: value
        type(mercier_d_terms_t) :: local_terms
        type(mercier_gradient_t) :: local_grad
        real(dp) :: local_full
        real(dp), allocatable :: local_beta(:, :)

        call mercier_surface_terms(equilibrium, surface, theta, zeta, &
            covariant_theta, covariant_zeta, covariant_theta_slope, &
            covariant_zeta_slope, flux_slope, flux_curvature, &
            volume_slope, volume_curvature, p, iota_slope, &
            local_terms, local_grad, local_full, local_beta)
        value = local_terms%mercier
    end function d_mercier_of

    subroutine build_fixture(equilibrium, surface, theta, zeta)
        type(gvec_cas3d_equilibrium_t), intent(out) :: equilibrium
        type(surface_data_t), intent(out) :: surface
        real(dp), allocatable, intent(out) :: theta(:), zeta(:)
        integer, parameter :: n_theta = 16, n_zeta = 16
        real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)
        integer :: j, k
        real(dp) :: tp, zp

        equilibrium%field_periods = 1
        equilibrium%poloidal_modes = [0, 1, 2]
        equilibrium%toroidal_modes = [0, 1, 2]
        allocate (theta(n_theta), zeta(n_zeta))
        do j = 1, n_theta
            theta(j) = real(j - 1, dp) / real(n_theta, dp)
        end do
        do k = 1, n_zeta
            zeta(k) = real(k - 1, dp) / real(n_zeta, dp)
        end do
        allocate (surface%jacobian(n_theta, n_zeta))
        allocate (surface%mod_b(n_theta, n_zeta))
        allocate (surface%b_theta(n_theta, n_zeta))
        allocate (surface%b_zeta(n_theta, n_zeta))
        allocate (surface%area_element(n_theta, n_zeta))
        do k = 1, n_zeta
            do j = 1, n_theta
                tp = two_pi * theta(j)
                zp = two_pi * zeta(k)
                surface%jacobian(j, k) = 2.0_dp + 0.3_dp * cos(tp) &
                    + 0.2_dp * cos(zp) + 0.1_dp * cos(tp - zp)
                surface%mod_b(j, k) = 1.5_dp + 0.2_dp * cos(tp) &
                    + 0.15_dp * cos(zp - tp)
                surface%b_theta(j, k) = 0.4_dp + 0.1_dp * sin(tp) &
                    + 0.05_dp * cos(zp)
                surface%b_zeta(j, k) = 0.9_dp + 0.15_dp * cos(zp) &
                    + 0.05_dp * sin(tp - zp)
                surface%area_element(j, k) = 1.0_dp &
                    + 0.25_dp * cos(tp) * cos(zp)
            end do
        end do
    end subroutine build_fixture

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, "(a)") message
            error stop 1
        end if
    end subroutine require

end program test_mercier_full_gradient
