program test_mercier_fluxslope_gradient
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
    real(dp), parameter :: tol = 1.0e-6_dp

    type(gvec_cas3d_equilibrium_t) :: equilibrium
    type(surface_data_t) :: surface
    real(dp), allocatable :: theta(:), zeta(:)
    real(dp), allocatable :: beta_values(:, :)
    type(mercier_d_terms_t) :: d_terms
    type(mercier_gradient_t) :: grad
    real(dp) :: full, poloidal_base
    real(dp) :: analytic_toroidal, analytic_poloidal
    real(dp) :: fd_toroidal, fd_poloidal, h, rel_err

    call build_fixture(equilibrium, surface, theta, zeta)
    poloidal_base = sum(surface%jacobian * surface%b_theta) &
        / real(size(surface%jacobian), dp)

    call mercier_surface_terms(equilibrium, surface, theta, zeta, &
        covariant_theta, covariant_zeta, covariant_theta_slope, &
        covariant_zeta_slope, flux_slope, flux_curvature, volume_slope, &
        volume_curvature, pressure_slope, iota_slope, d_terms, grad, &
        full, beta_values, analytic_toroidal, analytic_poloidal, &
        poloidal_flux_slope_override=poloidal_base)

    ! Toroidal flux slope: re-run grad_psi, the beta solve, and the
    ! integrals at the perturbed toroidal slope while holding the poloidal
    ! slope fixed.
    h = 1.0e-5_dp * abs(flux_slope)
    fd_toroidal = (d_mercier_toroidal(flux_slope + h) &
        - d_mercier_toroidal(flux_slope - h)) / (2.0_dp * h)
    rel_err = abs(fd_toroidal - analytic_toroidal) &
        / max(abs(analytic_toroidal), tiny(1.0_dp))
    write (*, "(a, es24.16)") "d_mercier                 = ", d_terms%mercier
    write (*, "(a, es24.16)") "analytic d/d toroidal     = ", analytic_toroidal
    write (*, "(a, es24.16)") "finite diff toroidal      = ", fd_toroidal
    write (*, "(a, es24.16)") "relative error toroidal   = ", rel_err
    call require(rel_err <= tol, &
        "toroidal_flux_slope gradient disagrees with finite difference")

    ! Poloidal flux slope enters only through the beta denominator; the
    ! finite difference re-solves beta at the perturbed poloidal slope.
    h = 1.0e-5_dp * abs(poloidal_base)
    fd_poloidal = (d_mercier_poloidal(poloidal_base + h) &
        - d_mercier_poloidal(poloidal_base - h)) / (2.0_dp * h)
    rel_err = abs(fd_poloidal - analytic_poloidal) &
        / max(abs(analytic_poloidal), tiny(1.0_dp))
    write (*, "(a, es24.16)") "analytic d/d poloidal     = ", analytic_poloidal
    write (*, "(a, es24.16)") "finite diff poloidal      = ", fd_poloidal
    write (*, "(a, es24.16)") "relative error poloidal   = ", rel_err
    call require(rel_err <= tol, &
        "poloidal_flux_slope gradient disagrees with finite difference")
    call require(abs(analytic_poloidal) > 1.0e-9_dp, &
        "poloidal gradient must resolve a nonzero beta-implicit response")

    write (*, "(a)") "PASS"

contains

    function d_mercier_toroidal(phi_t) result(value)
        real(dp), intent(in) :: phi_t
        real(dp) :: value
        type(mercier_d_terms_t) :: local_terms
        type(mercier_gradient_t) :: local_grad
        real(dp) :: local_full, local_toroidal, local_poloidal
        real(dp), allocatable :: local_beta(:, :)

        call mercier_surface_terms(equilibrium, surface, theta, zeta, &
            covariant_theta, covariant_zeta, covariant_theta_slope, &
            covariant_zeta_slope, phi_t, flux_curvature, volume_slope, &
            volume_curvature, pressure_slope, iota_slope, local_terms, &
            local_grad, local_full, local_beta, local_toroidal, &
            local_poloidal, poloidal_flux_slope_override=poloidal_base)
        value = local_terms%mercier
    end function d_mercier_toroidal

    function d_mercier_poloidal(phi_p) result(value)
        real(dp), intent(in) :: phi_p
        real(dp) :: value
        type(mercier_d_terms_t) :: local_terms
        type(mercier_gradient_t) :: local_grad
        real(dp) :: local_full, local_toroidal, local_poloidal
        real(dp), allocatable :: local_beta(:, :)

        call mercier_surface_terms(equilibrium, surface, theta, zeta, &
            covariant_theta, covariant_zeta, covariant_theta_slope, &
            covariant_zeta_slope, flux_slope, flux_curvature, volume_slope, &
            volume_curvature, pressure_slope, iota_slope, local_terms, &
            local_grad, local_full, local_beta, local_toroidal, &
            local_poloidal, poloidal_flux_slope_override=phi_p)
        value = local_terms%mercier
    end function d_mercier_poloidal

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

end program test_mercier_fluxslope_gradient
